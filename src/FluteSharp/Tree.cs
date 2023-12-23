using System;
using System.Collections.Generic;

namespace Knapcode.FluteSharp;

public class Tree
{
    public int Deg { get; set; }
    public int Length { get; set; }
    public Branch[] Branch { get; set; } = Array.Empty<Branch>();

    public Dictionary<Point, HashSet<Point>> GetNeighbors()
    {
        var output = new Dictionary<Point, HashSet<Point>>();

        HashSet<Point> GetOrAddPoint(Point point)
        {
            if (!output!.TryGetValue(point, out var neighbors))
            {
                neighbors = new HashSet<Point>();
                output.Add(point, neighbors);
            }

            return neighbors;
        }


        for (int i = 0; i < Branch.Length; i++)
        {
            var branch = Branch[i];
            if (branch.N == i)
            {
                continue;
            }

            var nextBranch = Branch[branch.N];

            var point = new Point(branch.X, branch.Y);
            var nextPoint = new Point(nextBranch.X, nextBranch.Y);

            if (!point.Equals(nextPoint))
            {
                var neighbors = GetOrAddPoint(point);
                var nextNeighbors = GetOrAddPoint(nextPoint);

                neighbors.Add(nextPoint);
                nextNeighbors.Add(point);
            }
        }

        return output;
    }
}
