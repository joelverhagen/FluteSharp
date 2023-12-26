using System;
using System.Collections.Generic;

namespace Knapcode.FluteSharp;

public class Tree
{
    public int Deg { get; set; }
    public int Length { get; set; }
    public Branch[] Branch { get; set; } = Array.Empty<Branch>();

    public Dictionary<Point, List<Point>> GetNeighbors()
    {
        var pointToNeighbors = new Dictionary<Point, List<Point>>();
        var pointToAdded = new Dictionary<Point, Dictionary<Point, bool>>();

        (List<Point>, Dictionary<Point, bool>) GetOrAddPoint(Point point)
        {
            Dictionary<Point, bool> added;
            if (!pointToNeighbors!.TryGetValue(point, out var neighbors))
            {
                neighbors = new List<Point>();
                added = new Dictionary<Point, bool>();
                pointToNeighbors.Add(point, neighbors);
                pointToAdded.Add(point, added);
            }
            else
            {
                added = pointToAdded[point];
            }

            return (neighbors, added);
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
                (var neighbors, var added) = GetOrAddPoint(point);
                (var nextNeighbors, var nextAdded) = GetOrAddPoint(nextPoint);

#if NET6_0_OR_GREATER
                if (added.TryAdd(nextPoint, true))
                {
                    neighbors.Add(nextPoint);
                }

                if (nextAdded.TryAdd(nextPoint, true))
                {
                    nextNeighbors.Add(nextPoint);
                }
#else
                if (!added.ContainsKey(nextPoint))
                {
                    added.Add(nextPoint, true);
                    neighbors.Add(nextPoint);
                }

                if (!nextAdded.ContainsKey(nextPoint))
                {
                    nextAdded.Add(nextPoint, true);
                    nextNeighbors.Add(nextPoint);
                }
#endif
            }
        }

        return pointToNeighbors;
    }
}
