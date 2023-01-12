using System.Drawing;
using FluteSharp;
using static FluteSharp.LookUpTable.Constants;
using static FluteSharp.LookUpTable;
using System.Xml.Linq;

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

        var hd = new HighDegreeContext();
        return flutes(hd, points.Count, xs, ys, s, acc);
    }


    private class HighDegreeContext
    {
        public int FIRST_ROUND { get; set; } = 2;
        public int EARLY_QUIT_CRITERIA { get; set; } = 1;
        public int D3 { get; set; } = int.MaxValue;
    }

    private static Tree flutes(HighDegreeContext hd, int d, ReadOnlySpan<int> xs, ReadOnlySpan<int> ys, ReadOnlySpan<int> s, int acc)
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
            return flutes_HD(hd,d, xs, ys, s, acc);
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

    private static Tree flutes_HD(HighDegreeContext hd, int d, ReadOnlySpan<int> xs, ReadOnlySpan<int> ys, ReadOnlySpan<int> s, int acc)
    {
        int i, A, orig_D3;
        Tree t;
        //DTYPE *dist[MAXD], *dist_base;
        int[][] dist;
        int[] dist_base;
        int threshold, threshold_x, threshold_y;
        int best_round, min_node1, min_node2;
        int[] nb;
        int prev_len;

        //Chris
        if (d <= D2(acc))
        {
            if (acc <= 6)
            {
                hd.FIRST_ROUND = 0;
                A = acc;
            }
            else
            {
                hd.FIRST_ROUND = acc - 6;
                A = 6 + ((acc - 5) / 4) * 2;  // Even A is better
            }
            hd.EARLY_QUIT_CRITERIA = (int)(0.75 * hd.FIRST_ROUND + 0.5);

            dist_base = new int[d];
            dist = new int[d][];
            nb = new int[d] (int**)malloc(d * sizeof(int*));
            for (i = 0; i < d; i++)
            {
                dist[i] = &(dist_base[i * d]);
                nb[i] = (int*)malloc(DEFAULT_QSIZE * sizeof(int));
                nb[i][0] = DEFAULT_QSIZE;
                nb[i][1] = 2; // queue head
            }

            t = flute_mr(d, xs, ys, s, A, FIRST_ROUND,
                 dist, &threshold_x, &threshold_y, &threshold,
                 &best_round, &min_node1,
                 &min_node2, nb);
        }
        else
        {
            A = acc;
            orig_D3 = hd.D3;
            if (orig_D3 >= int.MinValue && d > 1000)
            {
                hd.D3 = (d <= 10000) ? 1000 : 10000;
            }
            t = flute_am(hd, d, xs, ys, s, A,
                 out threshold_x, out threshold_y, out threshold);
            if (orig_D3 >= int.MinValue)
            {
                hd.D3 = orig_D3;
            }
        }

        return t;
    }


    static Tree flute_am(HighDegreeContext hd, int d, ReadOnlySpan<int> xs, ReadOnlySpan<int> ys, ReadOnlySpan<int> s, int acc,
              out int threshold_x, out int threshold_y, out int threshold)
    {
        int i, j, k, m, n, itr, node1, node2;
        int smallest_gap, gap;
        Tree t, t0;
        Tree[] subtree;
        int prev_effort;
        /*
        int num_subtree, subroot[MAXPART], suproot[MAXPART], isSuperRoot[MAXD];
        int tree_id[MAXD], tid, tree_size[MAXD], edges[2*MAXD];
        int idx[MAXPART], offset[MAXPART], *order[MAXT],
              order_base[MAXD+10]; //order_base[MAXT*MAXD];
        DTYPE x[MAXD+MAXPART], y[MAXD+MAXPART];
        int new_s[MAXD+MAXPART], si[MAXD], xmap[MAXD+MAXPART];
        */
        int[] x, y;
        int num_subtree;
        int[] subroot = new int[3], suproot = new int[3], isSuperRoot;
        int[] tree_id, tree_size, edges;
        int tid;
        int[] idx = new int[3], offset = new int[3], order_base;
        Memory<int>[] order = new Memory<int>[3];
        int[] new_s, si, xmap;

        int maxd = d + 1;
        x = new int[maxd + 3];
        y = new int[maxd + 3];
        isSuperRoot = new int[maxd];
        tree_id = new int[maxd];
        tree_size = new int[maxd];
        edges = new int[maxd * 2];
        order_base = new int[maxd + 10];
        new_s = new int[maxd + 3];
        si = new int[maxd];
        xmap = new int[maxd + 3];

        /*
        for (i=0; i<MAXT; i++) {
          order[i] = &(order_base[i*MAXD]);
        }
        */

        for (i = 0; i < d; i++)
        {
            x[i] = xs[s[i]];
            y[i] = ys[i];
            isSuperRoot[i] = 0;
        }

        build_rmst((long)d, x, y, edges, tree_id);

        for (i = 0; i < d; i++)
        {
            tree_size[i] = 1;  // the node itself
        }

        suproot[0] = subroot[0] = edges[0];
        num_subtree = 1;

        /*
        for (i=2*d-3; i>=0; ) {
          node2 = edges[i--];
          node1 = edges[i--];
          j = tree_size[node1]+tree_size[node2];
          //if (j >= d/2) {
          if (j >= d/2 && num_subtree<2) {
            isSuperRoot[node1] = 1;
            suproot[num_subtree] = node1;
            subroot[num_subtree++] = node2;
          } else {
            tree_size[node1] = j;
          }
        }
        */

        for (i = 2 * d - 3; i >= 0;)
        {
            node2 = edges[i--];
            node1 = edges[i--];
            tree_size[node1] += tree_size[node2];
        }

        j = 0;
        smallest_gap = Math.Abs(tree_size[j] - (d / 2));
        for (i = 1; i < d; i++)
        {
            gap = Math.Abs(tree_size[i] - (d / 2));
            if (gap < smallest_gap)
            {
                j = i;
                smallest_gap = gap;
            }
        }

        for (i = 2 * d - 3; i >= 0;)
        {
            node2 = edges[i--];
            node1 = edges[i--];
            if (node2 == j)
            {
                isSuperRoot[node1] = 1;
                suproot[num_subtree] = node1;
                subroot[num_subtree++] = node2;
                tree_size[subroot[0]] -= tree_size[j];
                break;
            }
        }

        //assert(num_subtree<=MAXT);

        for (i = 1; i < num_subtree; i++)
        {
            tree_id[subroot[i]] = i + 1;
            tree_size[subroot[i]] += 1;  // to account for the link to parent tree
        }

        for (i = 0; i < 2 * d - 2;)
        {
            node1 = edges[i++];
            node2 = edges[i++];
            if (tree_id[node2] == 1)
            { // non-root node
                tree_id[node2] = tree_id[node1];
            }
        }

        // Find inverse si[] of s[]
        for (i = 0; i < d; i++)
            si[s[i]] = i;

        offset[1] = 0;
        for (i = 1; i < num_subtree; i++)
        {
            offset[i + 1] = offset[i] + tree_size[subroot[i - 1]];
        }
        //assert(offset[num_subtree]==d+num_subtree-1-tree_size[subroot[num_subtree-1]]);

        for (i = 0; i < num_subtree; i++)
        {
            order[i] = order_base.AsMemory().Slice(offset[i + 1]);
        }

        for (i = 0; i <= num_subtree; i++)
            idx[i] = 0;

        for (i = 0; i < d; i++)
        {
            tid = tree_id[si[i]];
            j = idx[tid]++;
            x[offset[tid] + j] = xs[i];
            xmap[i] = j;

            if (isSuperRoot[si[i]] > 0)
            {
                for (k = 1; k < num_subtree; k++)
                {
                    if (suproot[k] == si[i])
                    {
                        tid = k + 1;
                        j = idx[tid]++;
                        x[offset[tid] + j] = xs[i];
                        xmap[d - 1 + tid] = j;
                    }
                }
            }
        }

        for (i = 0; i <= num_subtree; i++)
            idx[i] = 0;

        for (i = 0; i < d; i++)
        {
            tid = tree_id[i];
            j = idx[tid]++;
            y[offset[tid] + j] = ys[i];
            new_s[offset[tid] + j] = xmap[s[i]];
            order[tid - 1].Span[j] = i;

            if (isSuperRoot[i] > 0)
            {
                for (k = 1; k < num_subtree; k++)
                {
                    if (suproot[k] == i)
                    {
                        tid = k + 1;
                        j = idx[tid]++;
                        y[offset[tid] + j] = ys[i];
                        new_s[offset[tid] + j] = xmap[d - 1 + tid];
                        order[tid - 1].Span[j] = i;
                    }
                }
            }
        }

        subtree = new Tree[num_subtree];
        for (i = 1; i <= num_subtree; i++)
        {
            if (tree_size[subroot[i - 1]] <= 1)
            {
                subtree[i - 1].deg = 0;
                continue;
            }

            t = flutes(hd, tree_size[subroot[i - 1]], x.AsSpan().Slice(offset[i]), y.AsSpan().Slice(offset[i]),
                    new_s.AsSpan().Slice(offset[i]), acc);
            subtree[i - 1] = t;
        }

        for (i = 1; i < num_subtree; i++)
        {
            //assert(tree_id[suproot[i]] != tree_id[subroot[i]]);

            t = xmergetree(subtree[tree_id[suproot[i]] - 1],
                   subtree[tree_id[subroot[i]] - 1],
                   order[tree_id[suproot[i]] - 1],
                   order[tree_id[subroot[i]] - 1],
                   xs[s[suproot[i]]], ys[suproot[i]]);

            subtree[tree_id[subroot[i]] - 1].deg = 0;
            subtree[tree_id[suproot[i]] - 1] = t;
        }

        t0 = subtree[0];

        t = t0;

        return t;
    }

    static Tree xmergetree(Tree t1, Tree t2, Memory<int> order1, Memory<int> order2,
        int cx, int cy)
    {
        int i, num, cnt, order_by_x = 1;
        Tree t;
        TreeNode* tn1, *tn2, *n1, *p1, **nodes;
        dl_t list_of_nodes = dl_alloc();
        DTYPE threshold_x, threshold_y;
        DTYPE min_x, max_x, max_len, len, gain;

        if (t1.deg <= 0)
        {
            for (i = 0; i < t2.deg; i++)
            {
                order1.Span[i] = order2.Span[i];
            }
            return t2;
        }
        else if (t2.deg <= 0)
        {
            return t1;
        }

        redirect(t1, cx, cy);
        redirect(t2, cx, cy);

        curr_mark = 0;
        tn1 = createRootedTree(t1, order1, 1, list_of_nodes);
        tn2 = createRootedTree(t2, order2, 2, list_of_nodes);

        num = dl_length(list_of_nodes);
        nodes = (TreeNode**)malloc(sizeof(TreeNode*) * num);
        i = 0;
        dl_forall(TreeNode *, list_of_nodes, n1) {
            nodes[i++] = n1;
        }
        dl_endfor;
        dl_clear(list_of_nodes);

        qsort(nodes, num, sizeof(TreeNode*), cmpNodeByYX);

        max_len = 0;
        min_x = max_x = nodes[0]->x;
        for (i = 0; i < num; i++)
        {
            n1 = nodes[i];
            p1 = n1->parent;
            if (p1)
            {
                len = ADIFF(n1->x, p1->x) + ADIFF(n1->y, p1->y);
                if (len > max_len)
                {
                    max_len = len;
                }
            }
            if (n1->x < min_x)
            {
                min_x = n1->x;
            }
            else if (n1->x > max_x)
            {
                max_x = n1->x;
            }
        }

        threshold_x = (max_x - min_x) / 4;
        threshold_y = (nodes[num - 1]->y - nodes[0]->y) / 4;

        threshold_x = min(threshold_x, max_len);
        threshold_y = min(threshold_y, max_len);

        for (cnt = (t1.deg + t2.deg) / 2; cnt > 0; cnt--)
        {
            gain = (order_by_x) ?
              exchange_branches_order_x(num, nodes, threshold_x, threshold_y, max_len) :
              exchange_branches_order_y(num, nodes, threshold_x, threshold_y, max_len);

            //assert(gain>=0);

            if (gain <= 0 && !order_by_x)
            {
                break;
            }
            if (cnt > 1)
            {
                collect_nodes(tn1, list_of_nodes);
                num = dl_length(list_of_nodes);
                if (num <= 1)
                {
                    break;
                }

                collect_nodes(tn2, list_of_nodes);
                if (dl_length(list_of_nodes) - num <= 1)
                {
                    break;
                }

                free(nodes);
                num = dl_length(list_of_nodes);
                nodes = (TreeNode**)malloc(sizeof(TreeNode*) * num);
                i = 0;
                dl_forall(TreeNode *, list_of_nodes, n1) {
                    nodes[i++] = n1;
                }
                dl_endfor;
                dl_clear(list_of_nodes);

                if (order_by_x)
                {
                    order_by_x = 0;
                    qsort(nodes, num, sizeof(TreeNode*), cmpNodeByXY);
                }
                else
                {
                    order_by_x = 1;
                    qsort(nodes, num, sizeof(TreeNode*), cmpNodeByYX);
                }
            }
        }

        dl_free(list_of_nodes);
        free(nodes);

        t = mergeRootedTrees(tn1, tn2, order1);

        free(t1.branch);
        free(t2.branch);

        return t;
    }

    private static TreeNode* createRootedTree(Tree t, int* order, int id, dl_t list_of_nodes)
    {
        int i, dd, n;
        TreeNode* root = 0, **nodes, *p;

        dd = t.deg * 2 - 2;
        nodes = (TreeNode**)malloc(sizeof(TreeNode*) * dd);
        for (i = 0; i < dd; i++)
        {
            nodes[i] = (TreeNode*)malloc(sizeof(TreeNode));
            nodes[i]->mark = curr_mark;
            nodes[i]->children = dl_alloc();
        }

        curr_mark++;
        for (i = 0; i < dd; i++)
        {
            nodes[i]->mark = curr_mark;
            n = t.branch[i].n;
            if (i == n)
            {
                if (i < t.deg)
                {
                    //assert(root==0);
                    nodes[i]->parent = 0;
                    root = nodes[i];
                }
                else
                {  /* must be redundant */
                    dl_free(nodes[i]->children);
                    free(nodes[i]);
                    nodes[i] = 0;
                    continue;
                }
            }
            else
            {
                p = nodes[n];
                nodes[i]->parent = p;
                dl_append(TreeNode *, p->children, nodes[i]);
            }
            nodes[i]->order = (i < t.deg) ? order[i] : -1;
            nodes[i]->id = id;
            nodes[i]->x = t.branch[i].x;
            nodes[i]->y = t.branch[i].y;

            /* len will be computed in update_subtree 
            nodes[i]->blen =
              ADIFF(t.branch[i].x, t.branch[n].x)+ADIFF(t.branch[i].y, t.branch[n].y);

            nodes[i]->e = nodes[i];
            nodes[i]->len =
              ADIFF(t.branch[i].x, t.branch[n].x)+ADIFF(t.branch[i].y, t.branch[n].y);
            */

            dl_append(TreeNode *, list_of_nodes, nodes[i]);
        }

        //assert(root);

        update_subtree(root, 0);

        for (i = 0; i < dd; i++)
        {
            if (nodes[i] && nodes[i]->mark != curr_mark)
            {
                dl_free(nodes[i]->children);
                free(nodes[i]);
            }
        }

        free(nodes);
        return root;
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

    private static int GetManhattanDistance(Point a, Point b)
    {
        return Math.Abs(a.X - b.X) + Math.Abs(a.Y - b.Y);
    }
}
