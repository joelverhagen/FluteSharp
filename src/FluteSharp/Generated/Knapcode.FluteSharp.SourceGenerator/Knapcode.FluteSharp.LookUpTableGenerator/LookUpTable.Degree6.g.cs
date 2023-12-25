﻿// <auto-generated/>

using System;

namespace Knapcode.FluteSharp;

public partial class LookUpTable
{
    private static readonly byte[] POWV6 = new byte[]
    {
        100, 61, 52, 10, 49, 10, 48, 36, 10, 48, 48, 10, 50, 10, 48, 52, 36, 10, 48, 49, 36, 10, 48, 50, 10, 48, 50, 10, 48, 50,
        10, 100, 61, 53, 10, 49, 10, 48, 36, 10, 48, 48, 10, 50, 10, 48, 53, 36, 10, 48, 50, 36, 10, 51, 10, 48, 53, 54, 36, 10,
        48, 50, 54, 36, 10, 48, 49, 50, 36, 10, 48, 51, 10, 48, 50, 10, 48, 50, 10, 48, 50, 10, 48, 51, 10, 48, 51, 10, 48, 51,
        10, 48, 51, 10, 51, 10, 48, 53, 54, 36, 10, 48, 49, 54, 36, 10, 48, 49, 50, 36, 10, 48, 49, 50, 10, 48, 49, 50, 10, 48,
        51, 10, 48, 51, 10, 48, 49, 50, 10, 50, 10, 48, 53, 36, 10, 48, 49, 36, 10, 48, 49, 56, 10, 50, 10, 48, 54, 36, 10, 48,
        49, 36, 10, 48, 50, 48, 10, 51, 10, 48, 53, 54, 36, 10, 48, 49, 53, 36, 10, 48, 49, 50, 36, 10, 50, 10, 48, 53, 54, 36,
        10, 48, 49, 50, 36, 10, 48, 50, 51, 10, 48, 50, 50, 10, 48, 50, 50, 10, 48, 50, 50, 10, 48, 50, 51, 10, 48, 50, 51, 10,
        100, 61, 54, 10, 49, 10, 48, 36, 10, 48, 48, 10, 50, 10, 48, 54, 36, 10, 48, 51, 36, 10, 51, 10, 48, 54, 55, 36, 10, 48,
        55, 51, 36, 10, 48, 50, 51, 36, 10, 52, 10, 48, 54, 55, 56, 36, 10, 49, 51, 36, 54, 10, 50, 50, 36, 55, 10, 51, 49, 36,
        56, 10, 48, 52, 10, 48, 50, 10, 48, 50, 10, 48, 50, 10, 48, 51, 10, 48, 52, 10, 48, 52, 10, 48, 51, 10, 48, 51, 10, 51,
        10, 48, 54, 55, 36, 10, 48, 50, 55, 36, 10, 48, 50, 51, 36, 10, 48, 49, 52, 10, 53, 10, 48, 54, 55, 56, 36, 10, 49, 50,
        36, 54, 10, 50, 51, 36, 55, 10, 51, 49, 54, 36, 56, 10, 52, 50, 36, 54, 10, 48, 49, 54, 10, 48, 52, 10, 48, 52, 10, 48,
        49, 54, 10, 53, 10, 48, 54, 55, 56, 36, 10, 49, 50, 36, 54, 10, 50, 49, 36, 55, 10, 51, 54, 51, 36, 56, 10, 52, 50, 36,
        54, 10, 48, 50, 49, 10, 48, 50, 49, 10, 48, 52, 10, 48, 52, 10, 48, 49, 54, 10, 48, 50, 49, 10, 51, 10, 48, 54, 55, 36,
        10, 48, 50, 55, 36, 10, 48, 49, 50, 36, 10, 48, 50, 56, 10, 50, 10, 48, 55, 36, 10, 48, 50, 36, 10, 48, 51, 48, 10, 51,
        10, 48, 54, 55, 36, 10, 48, 54, 50, 36, 10, 48, 50, 51, 36, 10, 50, 10, 48, 54, 55, 36, 10, 48, 50, 51, 36, 10, 52, 10,
        48, 54, 55, 56, 36, 10, 48, 50, 51, 56, 36, 10, 50, 49, 54, 36, 56, 10, 51, 50, 36, 54, 10, 48, 51, 52, 10, 48, 51, 50,
        10, 48, 51, 50, 10, 48, 51, 50, 10, 48, 51, 51, 10, 48, 51, 52, 10, 48, 51, 52, 10, 48, 51, 51, 10, 48, 51, 51, 10, 48,
        49, 52, 10, 48, 49, 52, 10, 48, 49, 54, 10, 48, 49, 54, 10, 48, 51, 52, 10, 48, 51, 52, 10, 48, 49, 54, 10, 48, 50, 49,
        10, 48, 50, 49, 10, 48, 50, 49, 10, 48, 51, 52, 10, 48, 51, 52, 10, 48, 49, 54, 10, 48, 50, 49, 10, 48, 50, 56, 10, 48,
        50, 56, 10, 51, 10, 48, 55, 56, 36, 10, 48, 50, 56, 36, 10, 48, 49, 50, 36, 10, 48, 54, 48, 10, 53, 10, 48, 54, 55, 56,
        36, 10, 49, 50, 36, 55, 10, 50, 49, 36, 56, 10, 50, 51, 36, 54, 10, 51, 50, 51, 36, 54, 10, 56, 10, 48, 54, 55, 55, 56,
        36, 10, 49, 49, 36, 55, 10, 49, 51, 36, 55, 10, 50, 50, 36, 56, 10, 51, 50, 36, 54, 10, 52, 51, 36, 55, 10, 54, 56, 36,
        54, 10, 54, 50, 36, 54, 10, 55, 10, 48, 54, 55, 55, 56, 36, 10, 49, 49, 36, 55, 10, 49, 51, 36, 55, 10, 51, 50, 36, 54,
        10, 52, 49, 36, 55, 10, 53, 54, 36, 56, 10, 54, 50, 36, 54, 10, 48, 54, 52, 10, 48, 54, 50, 10, 48, 54, 50, 10, 48, 54,
        50, 10, 48, 54, 51, 10, 48, 54, 52, 10, 48, 54, 52, 10, 48, 54, 51, 10, 48, 54, 51, 10, 53, 10, 48, 49, 55, 56, 36, 10,
        49, 50, 36, 56, 10, 50, 51, 36, 55, 10, 49, 54, 55, 36, 49, 10, 52, 51, 36, 55, 10, 48, 55, 52, 10, 52, 10, 48, 49, 55,
        56, 36, 10, 48, 49, 50, 51, 36, 10, 49, 54, 55, 36, 49, 10, 51, 51, 36, 55, 10, 48, 55, 54, 10, 48, 54, 52, 10, 48, 54,
        52, 10, 48, 55, 54, 10, 52, 10, 48, 54, 55, 56, 36, 10, 49, 49, 36, 54, 10, 50, 50, 36, 55, 10, 51, 51, 36, 56, 10, 48,
        56, 49, 10, 48, 56, 49, 10, 48, 54, 52, 10, 48, 54, 52, 10, 48, 55, 54, 10, 48, 56, 49, 10, 51, 10, 48, 54, 55, 36, 10,
        48, 49, 55, 36, 10, 48, 49, 50, 36, 10, 48, 56, 56, 10, 48, 54, 48, 10, 48, 54, 48, 10, 48, 54, 50, 10, 48, 54, 51, 10,
        53, 10, 48, 54, 55, 51, 36, 10, 49, 50, 36, 54, 10, 50, 49, 36, 55, 10, 49, 55, 56, 36, 51, 10, 52, 49, 36, 55, 10, 48,
        57, 52, 10, 48, 54, 50, 10, 48, 54, 50, 10, 48, 54, 50, 10, 48, 54, 51, 10, 48, 57, 52, 10, 48, 57, 52, 10, 48, 54, 51,
        10, 48, 54, 51, 10, 48, 55, 52, 10, 48, 55, 52, 10, 53, 10, 48, 54, 55, 51, 36, 10, 49, 49, 36, 54, 10, 50, 56, 36, 51,
        10, 50, 50, 36, 55, 10, 49, 55, 56, 36, 51, 10, 48, 49, 48, 54, 10, 48, 57, 52, 10, 48, 57, 52, 10, 48, 49, 48, 54, 10,
        52, 10, 48, 54, 56, 36, 10, 48, 54, 51, 36, 10, 48, 49, 56, 36, 10, 48, 49, 51, 36, 10, 48, 49, 49, 49, 10, 48, 49, 49,
        49, 10, 48, 57, 52, 10, 48, 57, 52, 10, 48, 49, 48, 54, 10, 48, 49, 49, 49, 10, 50, 10, 48, 54, 36, 10, 48, 49, 36, 10,
        48, 49, 49, 56, 10, 50, 10, 48, 56, 36, 10, 48, 49, 36, 10, 48, 49, 50, 48, 10, 52, 10, 48, 54, 56, 36, 10, 48, 49, 54,
        36, 10, 48, 51, 56, 36, 10, 48, 49, 51, 36, 10, 53, 10, 48, 49, 54, 55, 36, 10, 49, 51, 36, 54, 10, 50, 56, 36, 49, 10,
        50, 50, 36, 55, 10, 49, 55, 56, 36, 49, 10, 52, 10, 48, 55, 51, 56, 36, 10, 48, 49, 50, 51, 36, 10, 49, 54, 55, 36, 51,
        10, 51, 49, 36, 55, 10, 48, 49, 50, 52, 10, 48, 49, 50, 50, 10, 48, 49, 50, 50, 10, 48, 49, 50, 50, 10, 48, 49, 50, 51,
        10, 48, 49, 50, 52, 10, 48, 49, 50, 52, 10, 48, 49, 50, 51, 10, 48, 49, 50, 51, 10, 53, 10, 48, 49, 54, 55, 36, 10, 49,
        50, 36, 54, 10, 50, 51, 36, 55, 10, 49, 55, 56, 36, 49, 10, 52, 51, 36, 55, 10, 48, 49, 51, 52, 10, 55, 10, 48, 54, 55,
        55, 56, 36, 10, 49, 49, 36, 55, 10, 49, 51, 36, 55, 10, 50, 50, 36, 54, 10, 52, 51, 36, 55, 10, 53, 54, 36, 56, 10, 54,
        50, 36, 54, 10, 48, 49, 51, 54, 10, 48, 49, 50, 52, 10, 48, 49, 50, 52, 10, 48, 49, 51, 54, 10, 52, 10, 48, 54, 55, 56,
        36, 10, 48, 49, 50, 56, 36, 10, 50, 54, 51, 36, 56, 10, 51, 50, 36, 54, 10, 48, 49, 52, 49, 10, 48, 49, 52, 49, 10, 48,
        49, 50, 52, 10, 48, 49, 50, 52, 10, 48, 49, 51, 54, 10, 48, 49, 52, 49, 10, 50, 10, 48, 54, 55, 36, 10, 48, 49, 50, 36,
        10, 48, 49, 52, 56, 10, 51, 10, 48, 55, 56, 36, 10, 48, 49, 55, 36, 10, 48, 49, 50, 36, 10, 48, 49, 53, 48, 10, 52, 10,
        48, 54, 55, 56, 36, 10, 49, 49, 36, 56, 10, 50, 50, 36, 55, 10, 51, 51, 36, 54, 10, 52, 10, 48, 49, 54, 55, 36, 10, 48,
        49, 50, 51, 36, 10, 49, 55, 56, 36, 49, 10, 51, 51, 36, 55, 10, 54, 10, 48, 54, 55, 55, 56, 36, 10, 49, 49, 36, 55, 10,
        49, 51, 36, 55, 10, 48, 49, 54, 50, 51, 36, 10, 52, 56, 36, 54, 10, 52, 50, 36, 54, 10, 48, 49, 53, 52, 10, 48, 49, 53,
        50, 10, 48, 49, 53, 50, 10, 48, 49, 53, 50, 10, 48, 49, 53, 51, 10, 48, 49, 53, 52, 10, 48, 49, 53, 52, 10, 48, 49, 53,
        51, 10, 48, 49, 53, 51, 10, 48, 49, 51, 52, 10, 48, 49, 51, 52, 10, 48, 49, 51, 54, 10, 48, 49, 51, 54, 10, 48, 49, 53,
        52, 10, 48, 49, 53, 52, 10, 48, 49, 51, 54, 10, 48, 49, 52, 49, 10, 48, 49, 52, 49, 10, 48, 49, 52, 49, 10, 48, 49, 53,
        52, 10, 48, 49, 53, 52, 10, 48, 49, 51, 54, 10, 48, 49, 52, 49, 10, 48, 49, 52, 56, 10, 48, 49, 52, 56, 10
    };

    private static readonly byte[] POST6 = new byte[]
    {
        100, 61, 52, 10, 73, 88, 68, 85, 85, 10, 72, 88, 68, 85, 85, 10, 72, 73, 84, 84, 85, 10, 100, 61, 53, 10, 74, 89, 104, 85,
        103, 118, 102, 10, 74, 88, 104, 85, 103, 118, 119, 10, 74, 88, 89, 85, 103, 103, 119, 10, 73, 89, 105, 85, 103, 118, 119, 10, 73, 88,
        89, 85, 103, 103, 119, 10, 88, 89, 90, 118, 87, 86, 119, 10, 73, 89, 105, 85, 103, 118, 119, 10, 72, 73, 105, 101, 87, 118, 119, 10,
        88, 89, 90, 117, 87, 102, 102, 10, 72, 88, 106, 85, 103, 118, 119, 10, 72, 73, 106, 101, 87, 118, 119, 10, 73, 89, 104, 85, 103, 118,
        119, 10, 73, 74, 104, 101, 103, 118, 101, 10, 73, 89, 105, 85, 103, 118, 119, 10, 73, 74, 105, 101, 103, 119, 87, 10, 88, 89, 90, 117,
        118, 86, 119, 10, 73, 89, 105, 85, 103, 118, 119, 10, 88, 89, 90, 117, 103, 86, 119, 10, 100, 61, 54, 10, 75, 90, 105, 120, 102, 120,
        153, 120, 136, 10, 75, 90, 104, 120, 102, 120, 153, 120, 153, 10, 75, 90, 104, 105, 102, 120, 152, 121, 153, 10, 75, 89, 105, 121, 102, 120,
        153, 120, 153, 10, 75, 89, 104, 105, 102, 120, 152, 121, 153, 10, 75, 104, 105, 106, 102, 135, 151, 152, 153, 10, 74, 89, 105, 121, 102, 120,
        153, 120, 136, 10, 74, 89, 104, 105, 102, 120, 152, 121, 153, 10, 74, 104, 105, 106, 102, 135, 151, 152, 153, 10, 89, 90, 91, 104, 135, 105,
        137, 119, 118, 10, 75, 89, 105, 121, 102, 120, 153, 120, 153, 10, 75, 88, 89, 121, 102, 119, 153, 136, 153, 10, 75, 104, 105, 106, 102, 119,
        152, 152, 136, 10, 74, 89, 105, 121, 102, 120, 153, 120, 136, 10, 74, 88, 89, 121, 102, 119, 153, 136, 153, 10, 74, 104, 105, 106, 102, 119,
        152, 152, 153, 10, 89, 90, 91, 105, 135, 105, 137, 119, 118, 10, 104, 105, 106, 107, 152, 102, 151, 120, 136, 10, 73, 89, 105, 122, 102, 120,
        153, 120, 153, 10, 73, 88, 89, 122, 102, 119, 153, 136, 136, 10, 88, 89, 90, 122, 135, 102, 153, 120, 136, 10, 89, 105, 106, 107, 150, 103,
        152, 119, 120, 10, 104, 105, 106, 107, 151, 102, 152, 119, 120, 10, 73, 89, 105, 123, 102, 120, 153, 120, 136, 10, 73, 88, 89, 123, 102, 119,
        153, 136, 136, 10, 88, 89, 90, 123, 135, 102, 153, 120, 136, 10, 75, 89, 105, 120, 102, 120, 153, 120, 153, 10, 75, 89, 90, 120, 102, 120,
        153, 137, 121, 10, 75, 89, 105, 121, 102, 120, 153, 120, 153, 10, 75, 89, 90, 121, 102, 120, 153, 137, 121, 10, 75, 104, 105, 106, 102, 121,
        135, 152, 153, 10, 75, 89, 105, 121, 102, 120, 153, 120, 153, 10, 75, 104, 105, 106, 102, 120, 151, 152, 153, 10, 74, 89, 105, 121, 102, 120,
        153, 120, 136, 10, 74, 104, 105, 106, 102, 120, 151, 152, 153, 10, 89, 90, 91, 105, 135, 105, 137, 119, 118, 10, 104, 105, 106, 107, 152, 103,
        150, 120, 153, 10, 74, 90, 106, 120, 102, 120, 153, 120, 136, 10, 74, 89, 90, 120, 102, 120, 153, 136, 135, 10, 89, 90, 91, 120, 135, 104,
        153, 120, 134, 10, 74, 90, 106, 121, 102, 120, 153, 120, 136, 10, 74, 89, 90, 121, 102, 120, 153, 136, 135, 10, 89, 90, 91, 121, 135, 104,
        153, 120, 134, 10, 74, 104, 105, 106, 102, 121, 135, 152, 153, 10, 104, 105, 106, 107, 152, 105, 118, 120, 153, 10, 74, 90, 106, 122, 102, 120,
        153, 120, 153, 10, 89, 90, 106, 122, 118, 104, 153, 120, 136, 10, 73, 89, 105, 106, 102, 121, 152, 120, 153, 10, 89, 90, 91, 122, 134, 104,
        153, 119, 119, 10, 73, 88, 89, 106, 102, 121, 151, 136, 136, 10, 89, 105, 106, 107, 150, 105, 135, 120, 136, 10, 88, 89, 90, 106, 135, 105,
        150, 120, 153, 10, 104, 105, 106, 107, 151, 105, 134, 120, 153, 10, 74, 90, 106, 122, 102, 120, 153, 120, 153, 10, 89, 90, 106, 122, 118, 104,
        153, 120, 136, 10, 73, 89, 105, 106, 102, 121, 152, 120, 153, 10, 73, 88, 89, 106, 102, 121, 151, 136, 136, 10, 88, 89, 90, 106, 135, 105,
        150, 120, 153, 10, 89, 105, 106, 107, 150, 104, 151, 120, 136, 10, 104, 105, 106, 107, 151, 104, 150, 120, 153, 10, 88, 90, 106, 122, 118, 104,
        153, 120, 136, 10, 88, 90, 91, 122, 134, 104, 153, 119, 119, 10, 88, 105, 106, 107, 150, 105, 135, 120, 136, 10, 74, 90, 106, 122, 102, 120,
        153, 120, 153, 10, 73, 89, 105, 106, 102, 121, 152, 120, 153, 10, 88, 90, 106, 122, 118, 104, 153, 120, 136, 10, 88, 105, 106, 107, 150, 104,
        151, 120, 136, 10, 74, 90, 106, 122, 102, 120, 153, 120, 153, 10, 73, 89, 105, 106, 102, 121, 152, 120, 153, 10, 73, 89, 105, 122, 102, 120,
        153, 120, 153, 10, 72, 73, 105, 122, 118, 104, 153, 120, 136, 10, 88, 89, 90, 122, 134, 103, 153, 120, 136, 10, 88, 105, 106, 107, 150, 103,
        152, 119, 120, 10, 73, 89, 105, 123, 102, 120, 153, 120, 136, 10, 72, 73, 105, 123, 118, 104, 153, 120, 136, 10, 88, 89, 90, 123, 134, 103,
        153, 120, 136, 10, 73, 89, 105, 107, 102, 121, 152, 120, 153, 10, 73, 88, 89, 107, 102, 121, 151, 136, 136, 10, 88, 89, 90, 107, 135, 105,
        150, 120, 153, 10, 74, 90, 106, 122, 102, 120, 153, 119, 120, 10, 89, 90, 106, 122, 118, 104, 153, 120, 136, 10, 73, 89, 105, 107, 102, 121,
        152, 120, 153, 10, 88, 90, 106, 107, 118, 105, 152, 119, 120, 10, 88, 90, 106, 122, 118, 104, 153, 120, 136, 10, 88, 105, 106, 107, 134, 105,
        151, 120, 136, 10, 74, 90, 106, 122, 102, 120, 153, 119, 120, 10, 73, 89, 106, 122, 102, 120, 153, 119, 120, 10, 73, 89, 106, 107, 102, 121,
        152, 119, 120, 10, 72, 73, 106, 122, 118, 104, 153, 120, 153, 10, 88, 89, 106, 107, 118, 105, 152, 119, 120, 10, 72, 88, 106, 123, 102, 120,
        153, 120, 136, 10, 72, 73, 106, 123, 118, 104, 153, 120, 136, 10, 74, 90, 105, 120, 102, 120, 153, 120, 153, 10, 74, 75, 105, 120, 118, 120,
        153, 134, 136, 10, 74, 90, 105, 121, 102, 120, 153, 120, 136, 10, 74, 75, 105, 121, 118, 120, 153, 134, 136, 10, 74, 90, 104, 105, 102, 120,
        152, 119, 151, 10, 90, 91, 104, 105, 118, 120, 152, 150, 153, 10, 89, 91, 105, 121, 118, 120, 153, 134, 136, 10, 89, 91, 104, 105, 118, 120,
        152, 150, 153, 10, 74, 90, 104, 106, 102, 120, 152, 121, 153, 10, 91, 104, 105, 106, 104, 103, 151, 152, 153, 10, 74, 90, 106, 122, 102, 120,
        153, 120, 153, 10, 74, 90, 104, 106, 102, 120, 152, 121, 153, 10, 89, 90, 91, 104, 134, 121, 137, 119, 118, 10, 74, 90, 106, 122, 102, 120,
        153, 120, 153, 10, 89, 90, 105, 121, 118, 120, 153, 134, 136, 10, 89, 91, 105, 121, 118, 120, 153, 134, 136, 10, 88, 89, 91, 121, 134, 134,
        153, 121, 121, 10, 91, 104, 105, 106, 103, 103, 152, 152, 136, 10, 74, 90, 106, 122, 102, 120, 153, 120, 153, 10, 74, 90, 105, 106, 102, 120,
        152, 121, 153, 10, 74, 90, 106, 122, 102, 120, 153, 120, 153, 10, 89, 90, 105, 121, 118, 120, 153, 134, 136, 10, 74, 90, 105, 106, 102, 120,
        152, 121, 153, 10, 88, 89, 90, 121, 134, 134, 153, 120, 135, 10, 90, 104, 105, 106, 103, 103, 152, 152, 153, 10, 89, 90, 91, 105, 134, 121,
        137, 119, 118, 10, 104, 105, 106, 107, 150, 134, 151, 120, 136, 10, 73, 89, 105, 122, 102, 120, 153, 120, 153, 10, 88, 89, 90, 122, 134, 118,
        153, 120, 136, 10, 89, 105, 106, 107, 150, 103, 152, 119, 120, 10, 104, 105, 106, 107, 150, 118, 152, 119, 120, 10, 73, 89, 105, 123, 102, 120,
        153, 120, 136, 10, 88, 89, 90, 123, 134, 118, 153, 120, 136, 10, 74, 90, 106, 120, 102, 120, 153, 120, 136, 10, 74, 75, 106, 120, 118, 120,
        153, 134, 136, 10, 89, 90, 91, 120, 134, 135, 153, 120, 134, 10, 74, 90, 106, 121, 102, 120, 153, 120, 136, 10, 74, 75, 106, 121, 118, 120,
        153, 134, 136, 10, 89, 90, 91, 121, 134, 135, 153, 120, 134, 10, 91, 104, 105, 106, 103, 105, 135, 152, 153, 10, 89, 91, 105, 121, 118, 120,
        153, 134, 136, 10, 91, 104, 105, 106, 103, 104, 151, 152, 153, 10, 74, 90, 106, 122, 102, 120, 153, 120, 153, 10, 74, 90, 105, 106, 102, 120,
        152, 121, 153, 10, 74, 90, 106, 122, 102, 120, 153, 120, 153, 10, 89, 90, 105, 121, 118, 120, 153, 134, 136, 10, 74, 90, 105, 106, 102, 120,
        152, 121, 153, 10, 89, 90, 91, 105, 134, 121, 137, 119, 118, 10, 90, 104, 105, 106, 103, 104, 151, 152, 153, 10, 104, 105, 106, 107, 150, 135,
        150, 120, 153, 10
    };

    public static LookUpTable Degree6 => new LookUpTable(6, POWV6, POST6);
}
