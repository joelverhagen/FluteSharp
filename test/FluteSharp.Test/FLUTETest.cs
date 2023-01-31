using System.Collections.Immutable;
using System.Drawing;
using System.Text;

namespace Knapcode.FluteSharp;

public class FLUTETest
{
    private static readonly bool GenerateTestData = false;
    private const string InputFileName = "input.txt";

    private static readonly IDictionary<int, FLUTE> FLUTEs = Enumerable.Range(4, 6).ToDictionary(d => d, GetFLUTE);

    private static FLUTE GetFLUTE(int d)
    {
        using var powvStream = File.OpenRead(Path.Combine("data", $"POWV{d}.dat"));
        using var postStream = File.OpenRead(Path.Combine("data", $"POST{d}.dat"));

        var lut = new LookUpTable(d, powvStream, postStream);
        var flute = new FLUTE(lut);

        return flute;
    }

    [Theory]
    [MemberData(nameof(AllTestCaseNamesAndLookupTableSizes))]
    public void DataDrivenTest(string testDataCase, int d)
    {
        // Arrange
        var testDataDir = GetTestDataDir();
        var inputPath = Path.Combine(testDataDir, testDataCase, InputFileName);
        var input = File.ReadAllText(inputPath);
        var points = ReadGrid(input);
        if (GenerateTestData)
        {
            var normalizedInput = PrintGrid(points);
            if (normalizedInput != input)
            {
                File.WriteAllText(inputPath, normalizedInput);
            }
        }

        // Act
        var tree = FLUTEs[d].Execute(points);

        // Assert
        var neighbors = tree.GetNeighbors();
        var actualSolution = PrintGrid(points, neighbors);
        var solutionPath = Path.Combine(testDataDir, testDataCase, $"solution-flute-d-{d}.txt");
        if (GenerateTestData)
        {
            File.WriteAllText(solutionPath, actualSolution);
        }
        else
        {
            var expectedSolution = File.ReadAllText(solutionPath);
            Assert.Equal(expectedSolution, actualSolution);
        }
    }

    public static IEnumerable<object[]> AllTestCaseNamesAndLookupTableSizes
    {
        get
        {
            foreach (var testCase in TestCaseNames)
            {
                foreach (var d in FLUTEs.Keys)
                {
                    yield return new object[] { testCase, d };
                }
            }
        }
    }

    public static IEnumerable<string> TestCaseNames
    {
        get
        {
            var testDataDir = GetTestDataDir();
            var inputPaths = Directory.EnumerateFiles(testDataDir, InputFileName, SearchOption.AllDirectories);
            foreach (var inputPath in inputPaths)
            {
                var testDir = Path.GetDirectoryName(inputPath)!;
                var testName = testDir.Substring(testDataDir.Length + 1).Replace('\\', '/');
                yield return testName;
            }
        }
    }

    private static string GetTestDataDir()
    {
        var current = Directory.GetCurrentDirectory();
        while (current != null && !Directory.EnumerateFiles(current, "FluteSharp.sln").Any())
        {
            current = Path.GetDirectoryName(current);
        }

        if (current == null)
        {
            throw new InvalidOperationException($"Could not find the repository root by probing upwards for FluteSharp.sln. Current directory: {Directory.GetCurrentDirectory()}");
        }

        return Path.GetFullPath(Path.Combine(current, "test", "FluteSharp.Test", "test-data"));
    }

    public static List<Point> ReadGrid(string grid)
    {
        using var reader = new StringReader(grid);
        string? line;
        var y = 0;
        var points = new List<Point>();
        while ((line = reader.ReadLine()) != null)
        {
            if (string.IsNullOrWhiteSpace(line))
            {
                continue;
            }

            for (var x = 0; x < line.Length; x++)
            {
                if (line[x] == 'p')
                {
                    points.Add(new Point(x, y));
                }
            }

            y++;
        }

        return points;
    }

    private static string PrintGrid(IReadOnlyList<Point> points)
    {
        return PrintGrid(points, points.ToDictionary(p => p, p => new HashSet<Point>()));
    }

    private static string PrintGrid(IReadOnlyList<Point> points, Dictionary<Point, HashSet<Point>> neighbors)
    {
        var maxX = points.Concat(neighbors.Keys).Max(p => p.X);
        var maxY = points.Concat(neighbors.Keys).Max(p => p.Y);

        var pointToSymbol = new Dictionary<Point, char>();
        foreach (var point in points)
        {
            pointToSymbol.Add(point, 'p'); // input point
        }

        foreach (var point in neighbors.Keys)
        {
            if (!pointToSymbol.ContainsKey(point))
            {
                pointToSymbol.Add(point, 's'); // Steiner point
            }
        }

        var grid = new char[maxX + 1, maxY + 1];

        for (var y = 0; y <= maxY; y++)
        {
            for (var x = 0; x <= maxX; x++)
            {
                grid[x, y] = pointToSymbol.TryGetValue(new Point(x, y), out var symbol) ? symbol : '.';
            }
        }

        foreach ((var from, var tos) in neighbors)
        {
            foreach (var to in tos)
            {
                (var a, var b) = from.X > to.X || from.X == to.X && from.Y > to.Y ? (to, from) : (from, to);

                var tx = a.X > b.X ? -1 : 1;
                var dx = Math.Abs(a.X - b.X);
                for (var x = 0; x < dx; x++)
                {
                    var point = new Point(a.X + (x * tx), a.Y);
                    if (!pointToSymbol.ContainsKey(point))
                    {
                        grid[point.X, point.Y] = '+';
                    }
                }

                var ty = a.Y > b.Y ? -1 : 1;
                var dy = Math.Abs(a.Y - b.Y);
                for (var y = 0; y < dy; y++)
                {
                    var point = new Point(b.X, a.Y + (y * ty));
                    if (!pointToSymbol.ContainsKey(point))
                    {
                        grid[point.X, point.Y] = '+';
                    }
                }
            }
        }

        var sb = new StringBuilder();


        for (var y = 0; y < grid.GetLength(1); y++)
        {
            for (var x = 0; x < grid.GetLength(0); x++)
            {
                sb.Append(grid[x, y]);
            }

            sb.AppendLine();
        }

        return sb.ToString();
    }
}