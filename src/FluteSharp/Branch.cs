using System;

namespace Knapcode.FluteSharp;

public struct Branch : IEquatable<Branch>
{
    public int X;
    public int Y;
    public int N;

    public override bool Equals(object? obj)
    {
        return obj is Branch branch && Equals(branch);
    }

    public bool Equals(Branch other)
    {
        return X == other.X &&
               Y == other.Y &&
               N == other.N;
    }

    public override int GetHashCode()
    {
        return HashCode.Combine(X, Y, N);
    }

    public static bool operator ==(Branch left, Branch right)
    {
        return left.Equals(right);
    }

    public static bool operator !=(Branch left, Branch right)
    {
        return !(left == right);
    }
}
