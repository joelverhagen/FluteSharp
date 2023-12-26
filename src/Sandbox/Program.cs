using System;
using System.Collections.Generic;

namespace Knapcode.FluteSharp;

public class Program
{
    public static void Main()
    {
        var lut = LookUpTable.Degree6;

        var points = new List<Point>
        {
            new Point(94, 18),
            new Point(6, 19),
            new Point(65, 27),
            new Point(98, 27),
            new Point(72, 29),
            new Point(38, 45),
            new Point(50, 67),
            new Point(75, 69),
            new Point(95, 75),
            new Point(21, 96),
        };

        var flute = new FLUTE(lut);
        var flutetree = flute.Execute(points);

        // print all of the paths
        var paths = new Dictionary<string, bool>();
        for (int i = 0; i < flutetree.Branch.Length; i++)
        {
            var current = flutetree.Branch[i];

            while (true)
            {
                var next = flutetree.Branch[current.N];

                var path = $"{current.X} {current.Y} {next.X} {next.Y}";

                if (paths.TryAdd(path, true))
                {
                    Console.WriteLine(path);
                }

                if (current.N == next.N)
                {
                    break;
                }

                current = next;
            }
        }
    }
}