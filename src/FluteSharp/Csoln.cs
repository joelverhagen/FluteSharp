﻿/* Source: https://home.engineering.iastate.edu/~cnchu/Flute.java
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

internal class Csoln
{
    public Csoln(int d)
    {
        seg = new byte[11];
        rowcol = new byte[d - 2];
        neighbor = new byte[2 * d - 2];
    }

    public byte parent;
    public byte[] seg;  // Add: 0..i, Sub: j..10; seg[i+1]=seg[j-1]=0
    public byte[] rowcol;  // row = rowcol[]/16, col = rowcol[]%16, 
    public byte[] neighbor;
}
