namespace Knapcode.FluteSharp;

public class Tree
{
    public int Deg { get; set; }
    public int Length { get; set; }
    public Branch[] Branch { get; set; } = Array.Empty<Branch>();
}
