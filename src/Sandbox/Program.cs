using Knapcode.FluteSharp;

var d = 9;

using var powvStream = File.OpenRead(Path.Combine("data", $"POWV{d}.dat"));
using var postStream = File.OpenRead(Path.Combine("data", $"POST{d}.dat"));

var lut = new LookUpTable(d, powvStream, postStream);

string? line;
var points = new List<Point>();
while ((line = Console.In.ReadLine()) != null)
{
    var pieces = line.Split(' ', StringSplitOptions.RemoveEmptyEntries);
    points.Add(new Point(int.Parse(pieces[0]), int.Parse(pieces[1])));
}

var flute = new FLUTE(lut);
var flutetree = flute.Execute(points);

// print all of the paths
var paths = new HashSet<string>();
for (int i = 0; i < flutetree.Branch.Length; i++)
{
    var current = flutetree.Branch[i];

    while (true)
    {
        var next = flutetree.Branch[current.N];

        var path = $"{current.X} {current.Y} {next.X} {next.Y}";

        if (paths.Add(path))
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
