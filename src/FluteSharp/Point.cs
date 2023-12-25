using System;

namespace Knapcode.FluteSharp;

public struct Point : IEquatable<Point>
{
    public Point(int x, int y)
    {
        X = x;
        Y = y;
    }

    public int X;
    public int Y;

    public override bool Equals(object? obj)
    {
        return obj is Point point && Equals(point);
    }

    public bool Equals(Point other)
    {
        return X == other.X &&
               Y == other.Y;
    }

    public override int GetHashCode()
    {
        return HashCode.Combine(X, Y);
    }

    public static bool operator ==(Point left, Point right)
    {
        return left.Equals(right);
    }

    public static bool operator !=(Point left, Point right)
    {
        return !(left == right);
    }
}