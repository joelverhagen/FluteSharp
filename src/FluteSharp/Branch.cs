using System.Diagnostics;

namespace Knapcode.FluteSharp;

[DebuggerDisplay("X = {X}, Y = {Y}, N = {N}")]
public struct Branch
{
    public int X { get; set; }
    public int Y { get; set; }
    public int N { get; set; }
}
