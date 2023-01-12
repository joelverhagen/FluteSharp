﻿using System.Drawing;
using FluteSharp;
using static FluteSharp.LookUpTable.Constants;
using static FluteSharp.LookUpTable;

internal class Program
{
    private const int MAXD = 150;

    private static void Main(string[] args)
    {
        Initialize();

        string? line;
        var points = new List<Point>();
        while ((line = Console.In.ReadLine()) != null)
        {
            var pieces = line.Split(' ', StringSplitOptions.RemoveEmptyEntries);
            points.Add(new Point(int.Parse(pieces[0]), int.Parse(pieces[1])));
        }

        var flutetree = Flute(points, 3);
        Console.WriteLine($"FLUTE wirelength = {flutetree.length}");
    }

    public static Tree Flute(IReadOnlyList<Point> points, int acc)
    {
        if (points.Count == 1)
        {
            return new Tree
            {
                deg = 1,
                length = 0,
                branch = new[]
                {
                    new Branch
                    {
                        x = points[0].X,
                        y = points[0].Y,
                        n = 0,
                    },
                }
            };
        }
        else if (points.Count == 2)
        {
            return new Tree
            {
                deg = 2,
                length = GetManhattanDistance(points[0], points[1]),
                branch = new[] 
                {
                    new Branch
                    {
                        x = points[0].X,
                        y = points[0].Y,
                        n = 1,
                    },
                    new Branch
                    {
                        x = points[1].X,
                        y = points[1].Y,
                        n = 1,
                    },
                }
            };
        }

        var pointsByX = points.OrderBy(p => p.X).ToList();
        var pointsAndIndexByY = pointsByX.Select((p, i) => (Point: p, Index: i)).OrderBy(p => p.Point.Y).ToList();

        var xs = new int[points.Count];
        var ys = new int[points.Count];
        var s = new int[points.Count];

        for (var i = 0; i < points.Count; i++)
        {
            xs[i] = pointsByX[i].X;
            ys[i] = pointsAndIndexByY[i].Point.Y;
            s[i] = pointsAndIndexByY[i].Index;
        }

        return flutes(points.Count, xs, ys, s, acc);
    }

    private static Tree flutes(int d, ReadOnlySpan<int> xs, ReadOnlySpan<int> ys, ReadOnlySpan<int> s, int acc)
    {
        if (d <= D)
        {
            return flutes_LD(d, xs, ys, s, acc);
        }
        else if (d <= D1(acc))
        {
            return flutes_MD(d, xs, ys, s, acc);
        }
        else
        {
            return flutes_HD(d, xs, ys, s, acc);
        }
    }

    private static Tree flutes_LMD(int d, ReadOnlySpan<int> xs, ReadOnlySpan<int> ys, ReadOnlySpan<int> s, int acc)
    {
        if (d <= D)
        {
            return flutes_LD(d, xs, ys, s, acc);
        }
        else
        {
            return flutes_MD(d, xs, ys, s, acc);
        }
    }
    

    private static Tree flutes_LD(int d, ReadOnlySpan<int> xs, ReadOnlySpan<int> ys, ReadOnlySpan<int> s, int acc)
    {
        int k, pi, i, j;
        Csoln[] rlistarr;
        int rlisti = 0;
        Csoln bestrlist;
        int[] dd = new int[2*D-2];  // 0..D-2 for v, D-1..2*D-3 for h
        int minl, sum;
        int[] l = new int[MPOWV + 1];
        int hflip;
        Tree t = new Tree();

        t.deg = d;
        t.branch = new Branch[2 * d - 2];
        if (d == 2) {
            minl = xs[1]-xs[0]+ys[1]-ys[0];
            t.branch[0].x = xs[s[0]];
            t.branch[0].y = ys[0];
            t.branch[0].n = 1;
            t.branch[1].x = xs[s[1]];
            t.branch[1].y = ys[1];
            t.branch[1].n = 1;
        }
        else if (d == 3) {
            minl = xs[2]-xs[0]+ys[2]-ys[0];
            t.branch[0].x = xs[s[0]];
            t.branch[0].y = ys[0];
            t.branch[0].n = 3;
            t.branch[1].x = xs[s[1]];
            t.branch[1].y = ys[1];
            t.branch[1].n = 3;
            t.branch[2].x = xs[s[2]];
            t.branch[2].y = ys[2];
            t.branch[2].n = 3;
            t.branch[3].x = xs[1];
            t.branch[3].y = ys[1];
            t.branch[3].n = 3;
        }
        else {
            k = 0;
            if (s[0] < s[2]) k++;
            if (s[1] < s[2]) k++;
        
            for (i=3; i<=d-1; i++) {  // p0=0 always, skip i=1 for symmetry
                pi = s[i];
                for (j=d-1; j>i; j--)
                    if (s[j] < s[i])
                        pi--;
                k = pi + (i+1)*k;
            }
        
            if (k < numgrp[d]) { // no horizontal flip
                hflip = 0;
                for (i=1; i<=d-3; i++) {
                    dd[i]=ys[i+1]-ys[i];
                    dd[d-1+i]=xs[i+1]-xs[i];
                }
            }
            else {
                hflip = 1;
                k=2*numgrp[d]-1-k;
                for (i=1; i<=d-3; i++) {
                    dd[i]=ys[i+1]-ys[i];
                    dd[d-1+i]=xs[d-1-i]-xs[d-2-i];
                }
            }
        
            minl = l[0] = xs[d-1]-xs[0]+ys[d-1]-ys[0];
            rlistarr = LUT[d,k];
            for (i=0; rlistarr[0].seg[i]>0; i++)
                minl += dd[rlistarr[0].seg[i]];
            bestrlist = rlistarr[0];
            l[1] = minl;
            j = 2;
            while (j <= numsoln[d,k]) {
                rlisti++;
                sum = l[rlistarr[rlisti].parent];
                for (i=0; rlistarr[rlisti].seg[i]>0; i++)
                    sum += dd[rlistarr[rlisti].seg[i]];
                for (i=10; rlistarr[rlisti].seg[i]>0; i--)
                    sum -= dd[rlistarr[rlisti].seg[i]];
                if (sum < minl) {
                    minl = sum;
                    bestrlist = rlistarr[rlisti];
                }
                l[j++] = sum;
            }
        
            t.branch[0].x = xs[s[0]];
            t.branch[0].y = ys[0];
            t.branch[1].x = xs[s[1]];
            t.branch[1].y = ys[1];
            for (i=2; i<d-2; i++) {
                t.branch[i].x = xs[s[i]];
                t.branch[i].y = ys[i];
                t.branch[i].n = bestrlist.neighbor[i];
            }
            t.branch[d-2].x = xs[s[d-2]];
            t.branch[d-2].y = ys[d-2];
            t.branch[d-1].x = xs[s[d-1]];
            t.branch[d-1].y = ys[d-1];
            if (hflip > 0) {
                if (s[1] < s[0]) {
                    t.branch[0].n = bestrlist.neighbor[1];
                    t.branch[1].n = bestrlist.neighbor[0];
                }
                else {
                    t.branch[0].n = bestrlist.neighbor[0];
                    t.branch[1].n = bestrlist.neighbor[1];
                }
                if (s[d-1] < s[d-2]) {
                    t.branch[d-2].n = bestrlist.neighbor[d-1];
                    t.branch[d-1].n = bestrlist.neighbor[d-2];
                }
                else {
                    t.branch[d-2].n = bestrlist.neighbor[d-2];
                    t.branch[d-1].n = bestrlist.neighbor[d-1];
                }
                for (i=d; i<2*d-2; i++) {
                    t.branch[i].x = xs[d-1-bestrlist.rowcol[i-d]%16];
                    t.branch[i].y = ys[bestrlist.rowcol[i-d]/16];
                    t.branch[i].n = bestrlist.neighbor[i];
                    }
            }
            else {  // !hflip
                if (s[0] < s[1]) {
                    t.branch[0].n = bestrlist.neighbor[1];
                    t.branch[1].n = bestrlist.neighbor[0];
                }
                else {
                    t.branch[0].n = bestrlist.neighbor[0];
                    t.branch[1].n = bestrlist.neighbor[1];
                }
                if (s[d-2] < s[d-1]) {
                    t.branch[d-2].n = bestrlist.neighbor[d-1];
                    t.branch[d-1].n = bestrlist.neighbor[d-2];
                }
                else {
                    t.branch[d-2].n = bestrlist.neighbor[d-2];
                    t.branch[d-1].n = bestrlist.neighbor[d-1];
                }
                for (i=d; i<2*d-2; i++) {
                    t.branch[i].x = xs[bestrlist.rowcol[i-d]%16];
                    t.branch[i].y = ys[bestrlist.rowcol[i-d]/16];
                    t.branch[i].n = bestrlist.neighbor[i];
                }
            }
        }
        t.length = minl;
    
        return t;
    }

    private static Tree flutes_MD(int d, ReadOnlySpan<int> xs, ReadOnlySpan<int> ys, ReadOnlySpan<int> s, int acc)
    {
        int[] x1 = new int[MAXD], x2 = new int[MAXD], y1 = new int[MAXD], y2 = new int[MAXD];
        int[] si = new int[MAXD], s1 = new int[MAXD], s2 = new int[MAXD];
        float[] score = new float[2 * MAXD], penalty = new float[MAXD];
        float pnlty, dx, dy;
        int ll, minl, coord1, coord2;
        int i, r, p, maxbp, bestbp = 0, bp, nbp, ub, lb, n1, n2, nn1 = 0, nn2 = 0, newacc;
        Tree t, t1, t2, bestt1, bestt2;
        int ms, mins, maxs, minsi, maxsi;
        int[] distx = new int[MAXD], disty = new int[MAXD];
        int xydiff;

        if (s[0] < s[d-1]) {
            ms = Math.Max(s[0], s[1]);
            for (i=2; i<=ms; i++)
                ms = Math.Max(ms, s[i]);
            if (ms <= d-3) {
                for (i=0; i<=ms; i++) {
                    x1[i] = xs[i];
                    y1[i] = ys[i];
                    s1[i] = s[i];
                }
                x1[ms+1] = xs[ms]; 
                y1[ms+1] = ys[ms]; 
                s1[ms+1] = ms+1;
            
                s2[0] = 0;
                for (i=1; i<=d-1-ms; i++)
                    s2[i] = s[i+ms]-ms;

                t1 = flutes_LMD(ms+2, x1, y1, s1, acc);
                t2 = flutes_LMD(d-ms, xs.Slice(ms), ys.Slice(ms), s2, acc);
                t = dmergetree(t1, t2);
            
                return t;
            }
        }
        else {  // (s[0] > s[d-1])
            ms = Math.Min(s[0], s[1]);
            for (i=2; i<=d-1-ms; i++)
                ms = Math.Min(ms, s[i]);
            if (ms >= 2) {
                x1[0] = xs[ms];
                y1[0] = ys[0];
                s1[0] = s[0]-ms+1;
                for (i=1; i<=d-1-ms; i++) {
                    x1[i] = xs[i+ms-1];
                    y1[i] = ys[i];
                    s1[i] = s[i]-ms+1;
                }
                x1[d-ms] = xs[d-1];
                y1[d-ms] = ys[d-1-ms];
                s1[d-ms] = 0;
            
                s2[0] = ms;
                for (i=1; i<=ms; i++)
                    s2[i] = s[i+d-1-ms];
            
                t1 = flutes_LMD(d+1-ms, x1, y1, s1, acc);
                t2 = flutes_LMD(ms+1, xs, ys.Slice(d-1-ms), s2, acc);
                t = dmergetree(t1, t2);
            
                return t;
            }
        }

    // Find inverse si[] of s[]
        for (r=0; r<d; r++)
            si[s[r]] = r;
    
    // Determine breaking directions and positions dp[]
        lb=(d-2*acc+2)/4;
        if (lb < 2) lb = 2;
        ub=d-1-lb;

        // Compute scores    
        const float AA = 0.6f;  // 2.0*BB
        const float BB = 0.3f;
        float CC = 7.4f/((d+10.0f)*(d-3.0f));
        float DD = 4.8f/(d-1);

        // Compute penalty[]    
        dx = CC*(xs[d-2]-xs[1]);
        dy = CC*(ys[d-2]-ys[1]);
        for (r = d/2, pnlty = 0; r>=2; r--, pnlty += dx)
        {
            penalty[r] = pnlty;
            penalty[d - 1 - r] = pnlty;
        }

        penalty[1] = pnlty;
        penalty[d-2] = pnlty;
        penalty[0] = pnlty;
        penalty[d-1] = pnlty; 
        
        for (r = d/2-1, pnlty = dy; r>=2; r--, pnlty += dy)
        {
            penalty[s[r]] += pnlty;
            penalty[s[d - 1 - r]] += pnlty;
        }
            
        penalty[s[1]] += pnlty;
        penalty[s[d-2]] += pnlty;
        penalty[s[0]] += pnlty;
        penalty[s[d-1]] += pnlty;
    //#define CC 0.16
    //#define v(r) ((r==0||r==1||r==d-2||r==d-1) ? d-3 : abs(d-1-r-r))
    //    for (r=0; r<d; r++)
    //        penalty[r] = v(r)*dx + v(si[r])*dy;

        // Compute distx[], disty[]
        xydiff = (xs[d-1] - xs[0]) - (ys[d-1] - ys[0]);
        if (s[0] < s[1])
        {
            mins = s[0];
            maxs = s[1];
        }
        else
        {
            mins = s[1];
            maxs = s[0];
        }
        if (si[0] < si[1])
        {
            minsi = si[0];
            maxsi = si[1];
        }
        else
        {
            minsi = si[1];
            maxsi = si[0];
        }
        for (r=2; r<=ub; r++) {
            if (s[r] < mins)
                mins = s[r];
            else if (s[r] > maxs)
                maxs = s[r];
            distx[r] = xs[maxs] - xs[mins];
            if (si[r] < minsi)
                minsi = si[r];
            else if (si[r] > maxsi)
                maxsi = si[r];
            disty[r] = ys[maxsi] - ys[minsi] + xydiff;
        }

        if (s[d - 2] < s[d - 1])
        {
            mins = s[d - 2];
            maxs = s[d - 1];
        }

        else
        {
            mins = s[d - 1];
            maxs = s[d - 2];
        }

        if (si[d-2] < si[d-1])
        {
            minsi = si[d - 2];
            maxsi = si[d - 1];
        }
        else
        {
            minsi = si[d - 1];
            maxsi = si[d - 2];
        }

        for (r=d-3; r>=lb; r--) {
            if (s[r] < mins)
                mins = s[r];
            else if (s[r] > maxs)
                maxs = s[r];
            distx[r] += xs[maxs] - xs[mins];
            if (si[r] < minsi)
                minsi = si[r];
            else if (si[r] > maxsi)
                maxsi = si[r];
            disty[r] += ys[maxsi] - ys[minsi];
        }

        nbp=0;
        for (r=lb; r<=ub; r++) {
            if (si[r]<=1)
                score[nbp] = (xs[r+1] - xs[r-1]) - penalty[r]
                    - AA*(ys[2]-ys[1]) - DD*disty[r];
            else if (si[r]>=d-2)
                score[nbp] = (xs[r+1] - xs[r-1]) - penalty[r]
                    - AA*(ys[d-2]-ys[d-3]) - DD*disty[r];
            else score[nbp] = (xs[r+1] - xs[r-1]) - penalty[r]
                     - BB*(ys[si[r]+1]-ys[si[r]-1]) - DD*disty[r];
            nbp++;
        
            if (s[r]<=1)
                score[nbp] = (ys[r+1] - ys[r-1]) - penalty[s[r]]
                    - AA*(xs[2]-xs[1]) - DD*distx[r];
            else if (s[r]>=d-2)
                score[nbp] = (ys[r+1] - ys[r-1]) - penalty[s[r]]
                    - AA*(xs[d-2]-xs[d-3]) - DD*distx[r];
            else score[nbp] = (ys[r+1] - ys[r-1]) - penalty[s[r]]
                     - BB*(xs[s[r]+1]-xs[s[r]-1]) - DD*distx[r];
            nbp++;
        }

        if (acc <= 3)
            newacc = 1;
        else {
            newacc = acc/2;
            if (acc >= nbp) acc = nbp-1;
        }
    
        minl = int.MaxValue;
        bestt1 = new Tree();
        bestt2 = new Tree();
        bestt1.branch = bestt2.branch = Array.Empty<Branch>();
        for (i=0; i<acc; i++) {
            maxbp = 0;
            for (bp=1; bp<nbp; bp++)
                if (score[maxbp] < score[bp]) maxbp = bp;
            score[maxbp] = -9e9f;

            p = BreakPt(maxbp, lb);
    // Breaking in p
            if (BreakInX(maxbp)) {  // break in x
                n1 = n2 = 0;
                for (r=0; r<d; r++) {
                    if (s[r] < p) {
                        s1[n1] = s[r];
                        y1[n1] = ys[r];
                        n1++;
                    }
                    else if (s[r] > p) {
                        s2[n2] = s[r]-p;
                        y2[n2] = ys[r];
                        n2++;
                    }
                    else { // if (s[r] == p)  i.e.,  r = si[p]
                        s1[n1] = p;  s2[n2] = 0;
                        y1[n1] = y2[n2] = ys[r];
                        nn1 = n1;  nn2 = n2;
                        n1++;  n2++;
                    }
                }

                t1 = flutes_LMD(p+1, xs, y1, s1, newacc);
                t2 = flutes_LMD(d-p, xs.Slice(p), y2, s2, newacc);
                ll = t1.length + t2.length;
                coord1 = t1.branch[t1.branch[nn1].n].y;
                coord2 = t2.branch[t2.branch[nn2].n].y;
                if (t2.branch[nn2].y > Math.Max(coord1, coord2))
                    ll -= t2.branch[nn2].y - Math.Max(coord1, coord2);
                else if (t2.branch[nn2].y < Math.Min(coord1, coord2))
                    ll -= Math.Min(coord1, coord2) - t2.branch[nn2].y;
            }
            else {  // if (!BreakInX(maxbp))
                n1 = n2 = 0;
                for (r=0; r<d; r++) {
                    if (si[r] < p) {
                        s1[si[r]] = n1;
                        x1[n1] = xs[r];
                        n1++;
                    }
                    else if (si[r] > p) {
                        s2[si[r]-p] = n2;
                        x2[n2] = xs[r];
                        n2++;
                    }
                    else { // if (si[r] == p)  i.e.,  r = s[p]
                        s1[p] = n1;  s2[0] = n2;
                        x1[n1] = x2[n2] = xs[r];
                        n1++;  n2++;
                    }
                }

                t1 = flutes_LMD(p+1, x1, ys, s1, newacc);
                t2 = flutes_LMD(d-p, x2, ys.Slice(p), s2, newacc);
                ll = t1.length + t2.length;
                coord1 = t1.branch[t1.branch[p].n].x;
                coord2 = t2.branch[t2.branch[0].n].x;
                if (t2.branch[0].x > Math.Max(coord1, coord2))
                    ll -= t2.branch[0].x - Math.Max(coord1, coord2);
                else if (t2.branch[0].x < Math.Min(coord1, coord2))
                    ll -= Math.Min(coord1, coord2) - t2.branch[0].x;
            }
            if (minl > ll) {
                minl = ll;
                bestt1 = t1;
                bestt2 = t2;
                bestbp = maxbp;
            }
        }

        if (BreakInX(bestbp)) {
            t = hmergetree(bestt1, bestt2, s);
        } else {
            t = vmergetree(bestt1, bestt2);
        }
    

        return t;
    }

    private static int BreakPt(int bp, int lb)
    {
        return ((bp) / 2 + lb);
    }

    private static bool BreakInX(int bp)
    {
        return ((bp) % 2 == 0);
    }



    private static Tree hmergetree(Tree t1, Tree t2, ReadOnlySpan<int> s)
    {
        int i, prev, curr, next, extra, offset1, offset2;
        int p, ii = 0, n1, n2, nn1 = 0, nn2 = 0;
        int coord1, coord2;
        Tree t = new Tree();

        t.deg = t1.deg + t2.deg - 1;
        t.length = t1.length + t2.length;
        t.branch = new Branch[2 * t.deg - 2];
        offset1 = t2.deg - 1;
        offset2 = 2 * t1.deg - 3;

        p = t1.deg - 1;
        n1 = n2 = 0;
        for (i = 0; i < t.deg; i++)
        {
            if (s[i] < p)
            {
                t.branch[i].x = t1.branch[n1].x;
                t.branch[i].y = t1.branch[n1].y;
                t.branch[i].n = t1.branch[n1].n + offset1;
                n1++;
            }
            else if (s[i] > p)
            {
                t.branch[i].x = t2.branch[n2].x;
                t.branch[i].y = t2.branch[n2].y;
                t.branch[i].n = t2.branch[n2].n + offset2;
                n2++;
            }
            else
            {
                t.branch[i].x = t2.branch[n2].x;
                t.branch[i].y = t2.branch[n2].y;
                t.branch[i].n = t2.branch[n2].n + offset2;
                nn1 = n1; nn2 = n2; ii = i;
                n1++; n2++;
            }
        }
        for (i = t.deg; i <= t.deg + t1.deg - 3; i++)
        {
            t.branch[i].x = t1.branch[i - offset1].x;
            t.branch[i].y = t1.branch[i - offset1].y;
            t.branch[i].n = t1.branch[i - offset1].n + offset1;
        }
        for (i = t.deg + t1.deg - 2; i <= 2 * t.deg - 4; i++)
        {
            t.branch[i].x = t2.branch[i - offset2].x;
            t.branch[i].y = t2.branch[i - offset2].y;
            t.branch[i].n = t2.branch[i - offset2].n + offset2;
        }
        extra = 2 * t.deg - 3;
        coord1 = t1.branch[t1.branch[nn1].n].y;
        coord2 = t2.branch[t2.branch[nn2].n].y;
        if (t2.branch[nn2].y > Math.Max(coord1, coord2))
        {
            t.branch[extra].y = Math.Max(coord1, coord2);
            t.length -= t2.branch[nn2].y - t.branch[extra].y;
        }
        else if (t2.branch[nn2].y < Math.Min(coord1, coord2))
        {
            t.branch[extra].y = Math.Min(coord1, coord2);
            t.length -= t.branch[extra].y - t2.branch[nn2].y;
        }
        else t.branch[extra].y = t2.branch[nn2].y;
        t.branch[extra].x = t2.branch[nn2].x;
        t.branch[extra].n = t.branch[ii].n;
        t.branch[ii].n = extra;

        prev = extra;
        curr = t1.branch[nn1].n + offset1;
        next = t.branch[curr].n;
        while (curr != next)
        {
            t.branch[curr].n = prev;
            prev = curr;
            curr = next;
            next = t.branch[curr].n;
        }
        t.branch[curr].n = prev;

        return t;
    }

    private static Tree vmergetree(Tree t1, Tree t2)
    {
        int i, prev, curr, next, extra, offset1, offset2;
        int coord1, coord2;
        Tree t = new Tree();

        t.deg = t1.deg + t2.deg - 1;
        t.length = t1.length + t2.length;
        t.branch = new Branch[2 * t.deg - 2];
        offset1 = t2.deg - 1;
        offset2 = 2 * t1.deg - 3;

        for (i = 0; i <= t1.deg - 2; i++)
        {
            t.branch[i].x = t1.branch[i].x;
            t.branch[i].y = t1.branch[i].y;
            t.branch[i].n = t1.branch[i].n + offset1;
        }
        for (i = t1.deg - 1; i <= t.deg - 1; i++)
        {
            t.branch[i].x = t2.branch[i - t1.deg + 1].x;
            t.branch[i].y = t2.branch[i - t1.deg + 1].y;
            t.branch[i].n = t2.branch[i - t1.deg + 1].n + offset2;
        }
        for (i = t.deg; i <= t.deg + t1.deg - 3; i++)
        {
            t.branch[i].x = t1.branch[i - offset1].x;
            t.branch[i].y = t1.branch[i - offset1].y;
            t.branch[i].n = t1.branch[i - offset1].n + offset1;
        }
        for (i = t.deg + t1.deg - 2; i <= 2 * t.deg - 4; i++)
        {
            t.branch[i].x = t2.branch[i - offset2].x;
            t.branch[i].y = t2.branch[i - offset2].y;
            t.branch[i].n = t2.branch[i - offset2].n + offset2;
        }
        extra = 2 * t.deg - 3;
        coord1 = t1.branch[t1.branch[t1.deg - 1].n].x;
        coord2 = t2.branch[t2.branch[0].n].x;
        if (t2.branch[0].x > Math.Max(coord1, coord2))
        {
            t.branch[extra].x = Math.Max(coord1, coord2);
            t.length -= t2.branch[0].x - t.branch[extra].x;
        }
        else if (t2.branch[0].x < Math.Min(coord1, coord2))
        {
            t.branch[extra].x = Math.Min(coord1, coord2);
            t.length -= t.branch[extra].x - t2.branch[0].x;
        }
        else t.branch[extra].x = t2.branch[0].x;
        t.branch[extra].y = t2.branch[0].y;
        t.branch[extra].n = t.branch[t1.deg - 1].n;
        t.branch[t1.deg - 1].n = extra;

        prev = extra;
        curr = t1.branch[t1.deg - 1].n + offset1;
        next = t.branch[curr].n;
        while (curr != next)
        {
            t.branch[curr].n = prev;
            prev = curr;
            curr = next;
            next = t.branch[curr].n;
        }
        t.branch[curr].n = prev;

        return t;
    }


    private static Tree dmergetree(Tree t1, Tree t2)
    {
        int i, d, prev, curr, next, offset1, offset2;
        Tree t = new Tree();

        t.deg = d = t1.deg + t2.deg - 2;
        t.length = t1.length + t2.length;
        t.branch = new Branch[2 * d - 2];
        offset1 = t2.deg - 2;
        offset2 = 2 * t1.deg - 4;

        for (i = 0; i <= t1.deg - 2; i++)
        {
            t.branch[i].x = t1.branch[i].x;
            t.branch[i].y = t1.branch[i].y;
            t.branch[i].n = t1.branch[i].n + offset1;
        }
        for (i = t1.deg - 1; i <= d - 1; i++)
        {
            t.branch[i].x = t2.branch[i - t1.deg + 2].x;
            t.branch[i].y = t2.branch[i - t1.deg + 2].y;
            t.branch[i].n = t2.branch[i - t1.deg + 2].n + offset2;
        }
        for (i = d; i <= d + t1.deg - 3; i++)
        {
            t.branch[i].x = t1.branch[i - offset1].x;
            t.branch[i].y = t1.branch[i - offset1].y;
            t.branch[i].n = t1.branch[i - offset1].n + offset1;
        }
        for (i = d + t1.deg - 2; i <= 2 * d - 3; i++)
        {
            t.branch[i].x = t2.branch[i - offset2].x;
            t.branch[i].y = t2.branch[i - offset2].y;
            t.branch[i].n = t2.branch[i - offset2].n + offset2;
        }

        prev = t2.branch[0].n + offset2;
        curr = t1.branch[t1.deg - 1].n + offset1;
        next = t.branch[curr].n;
        while (curr != next)
        {
            t.branch[curr].n = prev;
            prev = curr;
            curr = next;
            next = t.branch[curr].n;
        }
        t.branch[curr].n = prev;

        return t;
    }

    private static Tree flutes_HD(int d, ReadOnlySpan<int> xs, ReadOnlySpan<int> ys, ReadOnlySpan<int> s, int acc)
    {
        throw new NotImplementedException();
    }

    private static double TAU(int A)
    {
        return 8 + 1.3 * A;
    }

    private static int D1(int A)
    {
        return 25 + 120 / (A * A);
    }

    private static int D2(int A)
    {
        return A <= 6 ? 500 : 75 + 5 * A;
    }

    private static Tree Flutes_HD(int[] xs, int[] ys, int[] s, int acc)
    {
        throw new NotImplementedException();
    }

    private static int GetManhattanDistance(Point a, Point b)
    {
        return Math.Abs(a.X - b.X) + Math.Abs(a.Y - b.Y);
    }
}