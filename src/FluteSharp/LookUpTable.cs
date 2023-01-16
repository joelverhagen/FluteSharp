/* Source: https://home.engineering.iastate.edu/~cnchu/Flute.java
 * The Java-based implementation of the LUT reader is by Stefan Mücke. Only the LUT reader is copied from his Java implementation and ported to C#.
 */

/**
 * READ THIS LICENSE AGREEMENT CAREFULLY BEFORE USING THIS PRODUCT. BY USING
 * THIS PRODUCT YOU INDICATE YOUR ACCEPTANCE OF THE TERMS OF THE FOLLOWING
 * AGREEMENT. THESE TERMS APPLY TO YOU AND ANY SUBSEQUENT LICENSEE OF THIS
 * PRODUCT.
 * 
 * License Agreement for FLUTE
 * 
 * Copyright (c) 2004 by Dr. Chris C. N. Chu
 * All rights reserved
 * 
 * Copyright (c) 2012 by Stefan Mücke -- Java port
 * All rights reserved
 * 
 * ATTRIBUTION ASSURANCE LICENSE (adapted from the original BSD license)
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the conditions below are
 * met. These conditions require a modest attribution to Dr. Chris C. N. Chu
 * (the "Author").
 * 
 * 1. Redistributions of the source code, with or without modification (the
 *    "Code"), must be accompanied by any documentation and, each time
 *    the resulting executable program or a program dependent thereon is
 *    launched, a prominent display (e.g., splash screen or banner text) of
 *    the Author's attribution information, which includes:
 *    (a) Dr. Chris C. N. Chu ("AUTHOR"),
 *    (b) Iowa State University ("PROFESSIONAL IDENTIFICATION"), and
 *    (c) http://home.eng.iastate.edu/~cnchu/ ("URL").
 * 
 * 2. Users who intend to use the Code for commercial purposes will notify
 *    Author prior to such commercial use.
 * 
 * 3. Neither the name nor any trademark of the Author may be used to
 *    endorse or promote products derived from this software without
 *    specific prior written permission.
 * 
 * 4. Users are entirely responsible, to the exclusion of the Author and any
 *    other persons, for compliance with (1) regulations set by owners or
 *    administrators of employed equipment, (2) licensing terms of any other
 *    software, and (3) local, national, and international regulations
 *    regarding use, including those regarding import, export, and use of
 *    encryption software.
 * 
 * THIS FREE SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS" AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR OR ANY CONTRIBUTOR BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, EFFECTS OF UNAUTHORIZED OR MALICIOUS
 * NETWORK ACCESS; PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

namespace Knapcode.FluteSharp;

public class LookUpTable
{
    /// <summary>
    /// The default file name for the "potentially optimal wirelength vector" data.
    /// </summary>
    public const string POWVFile = "POWV9.dat";

    /// <summary>
    /// The default file name for the "potentially optimal Steiner tree" data.
    /// </summary>
    public const string POSTFile = "POST9.dat";

    /// <summary>
    /// Max. # of groups, 9! = 362880
    /// </summary>
    private const int MGROUP = 362880 / 4;

    /// <summary>
    /// Maximum degree supported by the look-up table.
    /// </summary>
    internal readonly int D = 9;

    /// <summary>
    /// Max. # of POWVs per group
    /// </summary>
    internal readonly int MPOWV = 79;

    internal static readonly IReadOnlyList<int> numgrp = new[] { 0, 0, 0, 0, 6, 30, 180, 1260, 10080, 90720 };

    internal readonly Csoln[,][] LUT;

    internal readonly int[,] numsoln;

    public LookUpTable()
    {
        LUT = new Csoln[D + 1, MGROUP][]; // storing 4 .. D
        numsoln = new int[D + 1, MGROUP];

        Initialize();
    }
    private void Initialize()
    {
        using var powvStream = File.OpenRead(POWVFile);
        using var postStream = File.OpenRead(POSTFile);

        readLUT(powvStream, postStream);
    }

    private void readLUT(Stream powvStream, Stream postStream)
    {
        char[] charnum = new char[256];
        char[] lineBuf = new char[32];
        int linep;
        char c;
        int[] number = new int[1];

        for (int i = 0; i <= 255; i++)
        {
            if ('0' <= i && i <= '9')
                charnum[i] = (char)(i - '0');
            else if (i >= 'A')
                charnum[i] = (char)(i - 'A' + 10);
            else
                // if (i=='$' || i=='\n' || ... )
                charnum[i] = (char)0;
        }

        for (int d = 4; d <= D; d++)
        {
            // d=%d\n
            linep = readLine(powvStream, lineBuf);
            linep = scanString(lineBuf, linep, "d=");
            linep = scanNumber(lineBuf, linep, number);
            d = number[0];
            scanEOL(lineBuf, linep);

            // d=%d\n
            linep = readLine(postStream, lineBuf);
            linep = scanString(lineBuf, linep, "d=");
            linep = scanNumber(lineBuf, linep, number);
            d = number[0];
            scanEOL(lineBuf, linep);

            for (int k = 0; k < numgrp[d]; k++)
            {
                int ns = charnum[powvStream.ReadByte() & 0xff];

                if (ns == 0)
                { // same as some previous group
                  // %d\n
                    linep = readLine(powvStream, lineBuf);
                    linep = scanNumber(lineBuf, linep, number);
                    int kk = number[0];
                    scanEOL(lineBuf, linep);
                    numsoln[d, k] = numsoln[d, kk];
                    LUT[d, k] = LUT[d, kk];
                }
                else
                {
                    powvStream.ReadByte(); // '\n'
                    numsoln[d, k] = ns;
                    Csoln[] p = new Csoln[ns];
                    for (int i = 0; i < ns; i++)
                    {
                        p[i] = new Csoln(D);
                    }
                    int poffset = 0; // C# workaround for C-style pointer arithmetic on 'p'
                    LUT[d, k] = p;
                    for (int i = 1; i <= ns; i++)
                    {
                        linep = readLine(powvStream, lineBuf);
                        p[poffset].parent = (byte)charnum[lineBuf[linep++]];
                        int j = 0;
                        while ((p[poffset].seg[j++] = (byte)charnum[lineBuf[linep++]]) != 0)
                        {
                        }
                        j = 10;
                        while ((p[poffset].seg[j--] = (byte)charnum[lineBuf[linep++]]) != 0)
                        {
                        }
                        int nn = 2 * d - 2;
                        readChars(postStream, lineBuf, d - 2);
                        linep = 0;
                        for (j = d; j < nn; j++)
                        {
                            c = charnum[lineBuf[linep++]];
                            p[poffset].rowcol[j - d] = (byte)c;
                        }
                        readChars(postStream, lineBuf, nn / 2 + 1); // last char \n
                        linep = 0;
                        for (j = 0; j < nn;)
                        {
                            c = lineBuf[linep++];
                            p[poffset].neighbor[j++] = (byte)(c / 16);
                            p[poffset].neighbor[j++] = (byte)(c % 16);
                        }
                        poffset++;
                    }
                }
            }
        }
    }

    private int readLine(Stream ins, char[] buf)
    {
        int c;
        int i = 0;
        while ((c = ins.ReadByte()) != -1)
        {
            if (c == '\n')
            {
                buf[i] = '\n';
                break;
            }
            buf[i++] = (char)(c & 0xff);
        }
        return 0;
    }

    private void readChars(Stream ins, char[] buf, int count)
    {
        for (int i = 0; i < count; i++)
        {
            buf[i] = (char)(ins.ReadByte() & 0xff);
        }
    }

    private int scanNumber(char[] buf, int offset, int[] result)
    {
        int c = buf[offset++];
        if (!isDigit(c))
            throw new InvalidDataException("Reading error. Expected digit but got '" + c + "'.");
        int number = c - '0';
        while (isDigit(c = buf[offset++]))
        {
            number = number * 10 + (c - '0');
        }
        result[0] = number;
        return offset - 1;
    }

    private int scanString(char[] buf, int offset, string str)
    {
        int len = str.Length;
        for (int i = 0; i < len; i++)
        {
            int c = buf[offset + i];
            if (c != str[i])
                throw new InvalidDataException("Reading error. Expected '"
                        + str[i]
                        + "' but got '"
                        + c
                        + "'.");
        }
        return offset + len;
    }

    private void scanEOL(char[] buf, int offset)
    {
        if (buf[offset] != '\n')
            throw new InvalidDataException("Reading error. Expected EOL but got '" + buf[offset] + "'.");
    }

    private bool isDigit(int c)
    {
        return '0' <= c && c <= '9';
    }
}
