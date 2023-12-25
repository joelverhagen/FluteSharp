/* Source: https://home.engineering.iastate.edu/~cnchu/Flute.java
 * The Java-based implementation of the LUT reader is by Stefan Mücke. Only the LUT reader is copied from his Java implementation and ported to C#.
 */

using System;
using System.Collections.Generic;
using System.IO;

namespace Knapcode.FluteSharp;

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

public partial class LookUpTable
{
    /// <summary>
    /// Max. # of groups
    /// </summary>
    private readonly int MGROUP;

    /// <summary>
    /// Maximum degree supported by the look-up table.
    /// </summary>
    internal readonly int D;

    /// <summary>
    /// Max. # of POWVs per group
    /// </summary>
    internal readonly int MPOWV;

    internal static readonly IReadOnlyList<int> numgrp = new[] { 0, 0, 0, 0, 6, 30, 180, 1260, 10080, 90720 };

    internal readonly Csoln[,][] LUT;

    internal readonly int[,] numsoln;

#if USE_STREAMS
    public LookUpTable(int d, Stream powvStream, Stream postStream)
        : this(d, powvStream.ReadByte, postStream.ReadByte)
    {
    }
#endif

    public LookUpTable(int d, string powvPath, string postPath)
        : this(d, File.ReadAllBytes(powvPath), File.ReadAllBytes(postPath))
    {
    }

    public LookUpTable(int d, byte[] powv, byte[] post)
        : this(d, ByteReader.NewReadByte(powv), ByteReader.NewReadByte(post))
    {
    }

    public LookUpTable(int d, Func<int> readBytePowv, Func<int> readBytePost)
    {
        if (d < 4 || d > 9)
        {
            throw new ArgumentOutOfRangeException(nameof(d), d, "The degree for the look-up table must be between 4 and 9, inclusive.");
        }
        
        if (d <= 7)
        {
            MGROUP = 5040 / 4; // 7! = 5040
            MPOWV = 15;
        }
        else if (d == 8)
        {
            MGROUP = 40320 / 4; // 8! = 40320
            MPOWV = 33;
        }
        else if (d == 9)
        {
            MGROUP = 362880 / 4; // 9! = 362880
            MPOWV = 79;
        }
        else
        {
            throw new NotImplementedException();
        }

        D = d;
        LUT = new Csoln[D + 1, MGROUP][]; // storing 4 .. D
        numsoln = new int[D + 1, MGROUP];

        readLUT(readBytePowv, readBytePost);
    }

    private void readLUT(Func<int> readBytePowv, Func<int> readBytePost)
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
            linep = readLine(readBytePowv, lineBuf);
            linep = scanString(lineBuf, linep, "d=");
            linep = scanNumber(lineBuf, linep, number);
            d = number[0];
            scanEOL(lineBuf, linep);

            // d=%d\n
            linep = readLine(readBytePost, lineBuf);
            linep = scanString(lineBuf, linep, "d=");
            linep = scanNumber(lineBuf, linep, number);
            d = number[0];
            scanEOL(lineBuf, linep);

            for (int k = 0; k < numgrp[d]; k++)
            {
                int ns = charnum[readBytePowv() & 0xff];

                if (ns == 0)
                { // same as some previous group
                  // %d\n
                    linep = readLine(readBytePowv, lineBuf);
                    linep = scanNumber(lineBuf, linep, number);
                    int kk = number[0];
                    scanEOL(lineBuf, linep);
                    numsoln[d, k] = numsoln[d, kk];
                    LUT[d, k] = LUT[d, kk];
                }
                else
                {
                    readBytePowv(); // '\n'
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
                        linep = readLine(readBytePowv, lineBuf);
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
                        readChars(readBytePost, lineBuf, d - 2);
                        linep = 0;
                        for (j = d; j < nn; j++)
                        {
                            c = charnum[lineBuf[linep++]];
                            p[poffset].rowcol[j - d] = (byte)c;
                        }
                        readChars(readBytePost, lineBuf, nn / 2 + 1); // last char \n
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

    private int readLine(Func<int> readByte, char[] buf)
    {
        int c;
        int i = 0;
        while ((c = readByte()) != -1)
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

    private void readChars(Func<int> readByte, char[] buf, int count)
    {
        for (int i = 0; i < count; i++)
        {
            buf[i] = (char)(readByte() & 0xff);
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

    private class ByteReader
    {
        private readonly byte[] _data;
        private bool _eof;
        private int _index;

        public static Func<int> NewReadByte(byte[] data)
        {
            var iterator = new ByteReader(data);
            return iterator.ReadByte;
        }

        public ByteReader(byte[] data)
        {
            _data = data;
            _eof = data.Length == 0;
        }

        public int ReadByte()
        {
            if (_eof)
            {
                return -1;
            }

            var b = _data[_index];

            _index++;
            if (_index >= _data.Length)
            {
                _eof = true;
            }

            return b;
        }
    }
}
