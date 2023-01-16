using System.Drawing;
using Knapcode.FluteSharp;

var lut = new LookUpTable();

string? line;
var points = new List<Point>();
while ((line = Console.In.ReadLine()) != null)
{
    var pieces = line.Split(' ', StringSplitOptions.RemoveEmptyEntries);
    points.Add(new Point(int.Parse(pieces[0]), int.Parse(pieces[1])));
}

var flute = new FLUTE(lut);
var flutetree = flute.Execute(points);
Console.WriteLine($"FLUTE wirelength = {flutetree.Length}");
