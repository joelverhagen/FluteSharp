using System.Drawing;

namespace Knapcode.FluteSharp;

public class FLUTE
{
    private readonly LookUpTable _lut;
    private readonly int _maxD;

    public FLUTE(LookUpTable lut) : this(lut, maxD: 150)
    {
    }

    public FLUTE(LookUpTable lut, int maxD)
    {
        _lut = lut;
        _maxD = maxD;
    }

    public Tree Execute(IReadOnlyList<Point> points)
    {
        return Execute(points, accuracy: 3);
    }

    public Tree Execute(IReadOnlyList<Point> points, int accuracy)
    {
        if (points.Count == 1)
        {
            return new Tree
            {
                Deg = 1,
                Length = 0,
                Branch = new[]
                {
                    new Branch
                    {
                        X = points[0].X,
                        Y = points[0].Y,
                        N = 0,
                    },
                }
            };
        }
        else if (points.Count == 2)
        {
            return new Tree
            {
                Deg = 2,
                Length = GetManhattanDistance(points[0], points[1]),
                Branch = new[] 
                {
                    new Branch
                    {
                        X = points[0].X,
                        Y = points[0].Y,
                        N = 1,
                    },
                    new Branch
                    {
                        X = points[1].X,
                        Y = points[1].Y,
                        N = 1,
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

        return flutes(points.Count, xs, ys, s, accuracy);
    }

    private static int GetManhattanDistance(Point a, Point b)
    {
        return Math.Abs(a.X - b.X) + Math.Abs(a.Y - b.Y);
    }

    private Tree flutes(int d, ReadOnlySpan<int> xs, ReadOnlySpan<int> ys, ReadOnlySpan<int> s, int acc)
    {
        if (d <= _lut.D)
        {
            return flutes_LD(d, xs, ys, s, acc);
        }
        else
        {
            return flutes_MD(d, xs, ys, s, acc);
        }
    }

    private Tree flutes_LD(int d, ReadOnlySpan<int> xs, ReadOnlySpan<int> ys, ReadOnlySpan<int> s, int acc)
    {
        int k, pi, i, j;
        Csoln[] rlistarr;
        int rlisti = 0;
        Csoln bestrlist;
        int[] dd = new int[2*_lut.D-2];  // 0..D-2 for v, D-1..2*D-3 for h
        int minl, sum;
        int[] l = new int[_lut.MPOWV + 1];
        int hflip;
        Tree t = new Tree();

        t.Deg = d;
        t.Branch = new Branch[2 * d - 2];
        if (d == 2) {
            minl = xs[1]-xs[0]+ys[1]-ys[0];
            t.Branch[0].X = xs[s[0]];
            t.Branch[0].Y = ys[0];
            t.Branch[0].N = 1;
            t.Branch[1].X = xs[s[1]];
            t.Branch[1].Y = ys[1];
            t.Branch[1].N = 1;
        }
        else if (d == 3) {
            minl = xs[2]-xs[0]+ys[2]-ys[0];
            t.Branch[0].X = xs[s[0]];
            t.Branch[0].Y = ys[0];
            t.Branch[0].N = 3;
            t.Branch[1].X = xs[s[1]];
            t.Branch[1].Y = ys[1];
            t.Branch[1].N = 3;
            t.Branch[2].X = xs[s[2]];
            t.Branch[2].Y = ys[2];
            t.Branch[2].N = 3;
            t.Branch[3].X = xs[1];
            t.Branch[3].Y = ys[1];
            t.Branch[3].N = 3;
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
        
            if (k < LookUpTable.numgrp[d]) { // no horizontal flip
                hflip = 0;
                for (i=1; i<=d-3; i++) {
                    dd[i]=ys[i+1]-ys[i];
                    dd[d-1+i]=xs[i+1]-xs[i];
                }
            }
            else {
                hflip = 1;
                k=2*LookUpTable.numgrp[d]-1-k;
                for (i=1; i<=d-3; i++) {
                    dd[i]=ys[i+1]-ys[i];
                    dd[d-1+i]=xs[d-1-i]-xs[d-2-i];
                }
            }
        
            minl = l[0] = xs[d-1]-xs[0]+ys[d-1]-ys[0];
            rlistarr = _lut.LUT[d,k];
            for (i=0; rlistarr[0].seg[i]>0; i++)
                minl += dd[rlistarr[0].seg[i]];
            bestrlist = rlistarr[0];
            l[1] = minl;
            j = 2;
            while (j <= _lut.numsoln[d,k]) {
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
        
            t.Branch[0].X = xs[s[0]];
            t.Branch[0].Y = ys[0];
            t.Branch[1].X = xs[s[1]];
            t.Branch[1].Y = ys[1];
            for (i=2; i<d-2; i++) {
                t.Branch[i].X = xs[s[i]];
                t.Branch[i].Y = ys[i];
                t.Branch[i].N = bestrlist.neighbor[i];
            }
            t.Branch[d-2].X = xs[s[d-2]];
            t.Branch[d-2].Y = ys[d-2];
            t.Branch[d-1].X = xs[s[d-1]];
            t.Branch[d-1].Y = ys[d-1];
            if (hflip > 0) {
                if (s[1] < s[0]) {
                    t.Branch[0].N = bestrlist.neighbor[1];
                    t.Branch[1].N = bestrlist.neighbor[0];
                }
                else {
                    t.Branch[0].N = bestrlist.neighbor[0];
                    t.Branch[1].N = bestrlist.neighbor[1];
                }
                if (s[d-1] < s[d-2]) {
                    t.Branch[d-2].N = bestrlist.neighbor[d-1];
                    t.Branch[d-1].N = bestrlist.neighbor[d-2];
                }
                else {
                    t.Branch[d-2].N = bestrlist.neighbor[d-2];
                    t.Branch[d-1].N = bestrlist.neighbor[d-1];
                }
                for (i=d; i<2*d-2; i++) {
                    t.Branch[i].X = xs[d-1-bestrlist.rowcol[i-d]%16];
                    t.Branch[i].Y = ys[bestrlist.rowcol[i-d]/16];
                    t.Branch[i].N = bestrlist.neighbor[i];
                    }
            }
            else {  // !hflip
                if (s[0] < s[1]) {
                    t.Branch[0].N = bestrlist.neighbor[1];
                    t.Branch[1].N = bestrlist.neighbor[0];
                }
                else {
                    t.Branch[0].N = bestrlist.neighbor[0];
                    t.Branch[1].N = bestrlist.neighbor[1];
                }
                if (s[d-2] < s[d-1]) {
                    t.Branch[d-2].N = bestrlist.neighbor[d-1];
                    t.Branch[d-1].N = bestrlist.neighbor[d-2];
                }
                else {
                    t.Branch[d-2].N = bestrlist.neighbor[d-2];
                    t.Branch[d-1].N = bestrlist.neighbor[d-1];
                }
                for (i=d; i<2*d-2; i++) {
                    t.Branch[i].X = xs[bestrlist.rowcol[i-d]%16];
                    t.Branch[i].Y = ys[bestrlist.rowcol[i-d]/16];
                    t.Branch[i].N = bestrlist.neighbor[i];
                }
            }
        }
        t.Length = minl;
    
        return t;
    }

    private Tree flutes_MD(int d, ReadOnlySpan<int> xs, ReadOnlySpan<int> ys, ReadOnlySpan<int> s, int acc)
    {
        int[] x1 = new int[_maxD], x2 = new int[_maxD], y1 = new int[_maxD], y2 = new int[_maxD];
        int[] si = new int[_maxD], s1 = new int[_maxD], s2 = new int[_maxD];
        float[] score = new float[2 * _maxD], penalty = new float[_maxD];
        float pnlty, dx, dy;
        int ll, minl, coord1, coord2;
        int i, r, p, maxbp, bestbp = 0, bp, nbp, ub, lb, n1, n2, nn1 = 0, nn2 = 0, newacc;
        Tree t, t1, t2, bestt1, bestt2;
        int ms, mins, maxs, minsi, maxsi;
        int[] distx = new int[_maxD], disty = new int[_maxD];
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

                t1 = flutes(ms+2, x1, y1, s1, acc);
                t2 = flutes(d-ms, xs.Slice(ms), ys.Slice(ms), s2, acc);
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
            
                t1 = flutes(d+1-ms, x1, y1, s1, acc);
                t2 = flutes(ms+1, xs, ys.Slice(d-1-ms), s2, acc);
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
        bestt1.Branch = bestt2.Branch = Array.Empty<Branch>();
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

                t1 = flutes(p+1, xs, y1, s1, newacc);
                t2 = flutes(d-p, xs.Slice(p), y2, s2, newacc);
                ll = t1.Length + t2.Length;
                coord1 = t1.Branch[t1.Branch[nn1].N].Y;
                coord2 = t2.Branch[t2.Branch[nn2].N].Y;
                if (t2.Branch[nn2].Y > Math.Max(coord1, coord2))
                    ll -= t2.Branch[nn2].Y - Math.Max(coord1, coord2);
                else if (t2.Branch[nn2].Y < Math.Min(coord1, coord2))
                    ll -= Math.Min(coord1, coord2) - t2.Branch[nn2].Y;
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

                t1 = flutes(p+1, x1, ys, s1, newacc);
                t2 = flutes(d-p, x2, ys.Slice(p), s2, newacc);
                ll = t1.Length + t2.Length;
                coord1 = t1.Branch[t1.Branch[p].N].X;
                coord2 = t2.Branch[t2.Branch[0].N].X;
                if (t2.Branch[0].X > Math.Max(coord1, coord2))
                    ll -= t2.Branch[0].X - Math.Max(coord1, coord2);
                else if (t2.Branch[0].X < Math.Min(coord1, coord2))
                    ll -= Math.Min(coord1, coord2) - t2.Branch[0].X;
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

    private static Tree dmergetree(Tree t1, Tree t2)
    {
        int i, d, prev, curr, next, offset1, offset2;
        Tree t = new Tree();

        t.Deg = d = t1.Deg + t2.Deg - 2;
        t.Length = t1.Length + t2.Length;
        t.Branch = new Branch[2 * d - 2];
        offset1 = t2.Deg - 2;
        offset2 = 2 * t1.Deg - 4;

        for (i = 0; i <= t1.Deg - 2; i++)
        {
            t.Branch[i].X = t1.Branch[i].X;
            t.Branch[i].Y = t1.Branch[i].Y;
            t.Branch[i].N = t1.Branch[i].N + offset1;
        }
        for (i = t1.Deg - 1; i <= d - 1; i++)
        {
            t.Branch[i].X = t2.Branch[i - t1.Deg + 2].X;
            t.Branch[i].Y = t2.Branch[i - t1.Deg + 2].Y;
            t.Branch[i].N = t2.Branch[i - t1.Deg + 2].N + offset2;
        }
        for (i = d; i <= d + t1.Deg - 3; i++)
        {
            t.Branch[i].X = t1.Branch[i - offset1].X;
            t.Branch[i].Y = t1.Branch[i - offset1].Y;
            t.Branch[i].N = t1.Branch[i - offset1].N + offset1;
        }
        for (i = d + t1.Deg - 2; i <= 2 * d - 3; i++)
        {
            t.Branch[i].X = t2.Branch[i - offset2].X;
            t.Branch[i].Y = t2.Branch[i - offset2].Y;
            t.Branch[i].N = t2.Branch[i - offset2].N + offset2;
        }

        prev = t2.Branch[0].N + offset2;
        curr = t1.Branch[t1.Deg - 1].N + offset1;
        next = t.Branch[curr].N;
        while (curr != next)
        {
            t.Branch[curr].N = prev;
            prev = curr;
            curr = next;
            next = t.Branch[curr].N;
        }
        t.Branch[curr].N = prev;

        return t;
    }

    private static Tree hmergetree(Tree t1, Tree t2, ReadOnlySpan<int> s)
    {
        int i, prev, curr, next, extra, offset1, offset2;
        int p, ii = 0, n1, n2, nn1 = 0, nn2 = 0;
        int coord1, coord2;
        Tree t = new Tree();

        t.Deg = t1.Deg + t2.Deg - 1;
        t.Length = t1.Length + t2.Length;
        t.Branch = new Branch[2 * t.Deg - 2];
        offset1 = t2.Deg - 1;
        offset2 = 2 * t1.Deg - 3;

        p = t1.Deg - 1;
        n1 = n2 = 0;
        for (i = 0; i < t.Deg; i++)
        {
            if (s[i] < p)
            {
                t.Branch[i].X = t1.Branch[n1].X;
                t.Branch[i].Y = t1.Branch[n1].Y;
                t.Branch[i].N = t1.Branch[n1].N + offset1;
                n1++;
            }
            else if (s[i] > p)
            {
                t.Branch[i].X = t2.Branch[n2].X;
                t.Branch[i].Y = t2.Branch[n2].Y;
                t.Branch[i].N = t2.Branch[n2].N + offset2;
                n2++;
            }
            else
            {
                t.Branch[i].X = t2.Branch[n2].X;
                t.Branch[i].Y = t2.Branch[n2].Y;
                t.Branch[i].N = t2.Branch[n2].N + offset2;
                nn1 = n1; nn2 = n2; ii = i;
                n1++; n2++;
            }
        }
        for (i = t.Deg; i <= t.Deg + t1.Deg - 3; i++)
        {
            t.Branch[i].X = t1.Branch[i - offset1].X;
            t.Branch[i].Y = t1.Branch[i - offset1].Y;
            t.Branch[i].N = t1.Branch[i - offset1].N + offset1;
        }
        for (i = t.Deg + t1.Deg - 2; i <= 2 * t.Deg - 4; i++)
        {
            t.Branch[i].X = t2.Branch[i - offset2].X;
            t.Branch[i].Y = t2.Branch[i - offset2].Y;
            t.Branch[i].N = t2.Branch[i - offset2].N + offset2;
        }
        extra = 2 * t.Deg - 3;
        coord1 = t1.Branch[t1.Branch[nn1].N].Y;
        coord2 = t2.Branch[t2.Branch[nn2].N].Y;
        if (t2.Branch[nn2].Y > Math.Max(coord1, coord2))
        {
            t.Branch[extra].Y = Math.Max(coord1, coord2);
            t.Length -= t2.Branch[nn2].Y - t.Branch[extra].Y;
        }
        else if (t2.Branch[nn2].Y < Math.Min(coord1, coord2))
        {
            t.Branch[extra].Y = Math.Min(coord1, coord2);
            t.Length -= t.Branch[extra].Y - t2.Branch[nn2].Y;
        }
        else t.Branch[extra].Y = t2.Branch[nn2].Y;
        t.Branch[extra].X = t2.Branch[nn2].X;
        t.Branch[extra].N = t.Branch[ii].N;
        t.Branch[ii].N = extra;

        prev = extra;
        curr = t1.Branch[nn1].N + offset1;
        next = t.Branch[curr].N;
        while (curr != next)
        {
            t.Branch[curr].N = prev;
            prev = curr;
            curr = next;
            next = t.Branch[curr].N;
        }
        t.Branch[curr].N = prev;

        return t;
    }

    private static Tree vmergetree(Tree t1, Tree t2)
    {
        int i, prev, curr, next, extra, offset1, offset2;
        int coord1, coord2;
        Tree t = new Tree();

        t.Deg = t1.Deg + t2.Deg - 1;
        t.Length = t1.Length + t2.Length;
        t.Branch = new Branch[2 * t.Deg - 2];
        offset1 = t2.Deg - 1;
        offset2 = 2 * t1.Deg - 3;

        for (i = 0; i <= t1.Deg - 2; i++)
        {
            t.Branch[i].X = t1.Branch[i].X;
            t.Branch[i].Y = t1.Branch[i].Y;
            t.Branch[i].N = t1.Branch[i].N + offset1;
        }
        for (i = t1.Deg - 1; i <= t.Deg - 1; i++)
        {
            t.Branch[i].X = t2.Branch[i - t1.Deg + 1].X;
            t.Branch[i].Y = t2.Branch[i - t1.Deg + 1].Y;
            t.Branch[i].N = t2.Branch[i - t1.Deg + 1].N + offset2;
        }
        for (i = t.Deg; i <= t.Deg + t1.Deg - 3; i++)
        {
            t.Branch[i].X = t1.Branch[i - offset1].X;
            t.Branch[i].Y = t1.Branch[i - offset1].Y;
            t.Branch[i].N = t1.Branch[i - offset1].N + offset1;
        }
        for (i = t.Deg + t1.Deg - 2; i <= 2 * t.Deg - 4; i++)
        {
            t.Branch[i].X = t2.Branch[i - offset2].X;
            t.Branch[i].Y = t2.Branch[i - offset2].Y;
            t.Branch[i].N = t2.Branch[i - offset2].N + offset2;
        }
        extra = 2 * t.Deg - 3;
        coord1 = t1.Branch[t1.Branch[t1.Deg - 1].N].X;
        coord2 = t2.Branch[t2.Branch[0].N].X;
        if (t2.Branch[0].X > Math.Max(coord1, coord2))
        {
            t.Branch[extra].X = Math.Max(coord1, coord2);
            t.Length -= t2.Branch[0].X - t.Branch[extra].X;
        }
        else if (t2.Branch[0].X < Math.Min(coord1, coord2))
        {
            t.Branch[extra].X = Math.Min(coord1, coord2);
            t.Length -= t.Branch[extra].X - t2.Branch[0].X;
        }
        else t.Branch[extra].X = t2.Branch[0].X;
        t.Branch[extra].Y = t2.Branch[0].Y;
        t.Branch[extra].N = t.Branch[t1.Deg - 1].N;
        t.Branch[t1.Deg - 1].N = extra;

        prev = extra;
        curr = t1.Branch[t1.Deg - 1].N + offset1;
        next = t.Branch[curr].N;
        while (curr != next)
        {
            t.Branch[curr].N = prev;
            prev = curr;
            curr = next;
            next = t.Branch[curr].N;
        }
        t.Branch[curr].N = prev;

        return t;
    }


}
