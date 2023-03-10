# FluteSharp

[![Latest version](https://img.shields.io/nuget/v/Knapcode.FluteSharp)](https://www.nuget.org/packages/Knapcode.FluteSharp) ![Build](https://github.com/joelverhagen/FluteSharp/workflows/Build/badge.svg)

A .NET port of [the FLUTE algorithm](https://home.engineering.iastate.edu/~cnchu/flute.html). It is a based on the [FLUTE 3.1 C implementation](https://home.engineering.iastate.edu/~cnchu/flute-3.1/flute-3.1.tgz). Note that I have not ported the high degree optimizations that Chris added in FLUTE 3.0 and 3.1 (e.g. the `flutes_HD` method). For my own purposes I think I can get by with the simpler implementation.

**Dr. Chris C. N. Chu** is the author of the algorithm and the associated C implementation. He is a professor at **Iowa State University**. His home page URL is **http://home.engineering.iastate.edu/~cnchu/**.

The Java-based implementation of the LUT reader is by Stefan Mücke. Only the LUT reader is copied from [his Java implementation](https://home.engineering.iastate.edu/~cnchu/Flute.java) and ported to C#.

## Install

```plaintext
dotnet add package Knapcode.FluteSharp
```

## Example

See the sample in [`src/Sandbox/Program.cs`](https://github.com/joelverhagen/FluteSharp/blob/main/src/Sandbox/Program.cs)

## FLUTE license

This license is copied from the [FLUTE website](https://home.engineering.iastate.edu/~cnchu/flute.html#License).

```plaintext
READ THIS LICENSE AGREEMENT CAREFULLY BEFORE USING THIS PRODUCT. BY USING THIS PRODUCT YOU INDICATE YOUR ACCEPTANCE OF THE TERMS OF THE FOLLOWING AGREEMENT. THESE TERMS APPLY TO YOU AND ANY SUBSEQUENT LICENSEE OF THIS PRODUCT.

License Agreement for FLUTE

Copyright (c) 2004 by Dr. Chris C. N. Chu
All rights reserved

ATTRIBUTION ASSURANCE LICENSE (adapted from the original BSD license) Redistribution and use in source and binary forms, with or without modification, are permitted provided that the conditions below are met. These conditions require a modest attribution to Dr. Chris C. N. Chu (the "Author").

Redistributions of the source code, with or without modification (the "Code"), must be accompanied by any documentation and, each time the resulting executable program or a program dependent thereon is launched, a prominent display (e.g., splash screen or banner text) of the Author's attribution information, which includes:
(a) Dr. Chris C. N. Chu ("AUTHOR"),
(b) Iowa State University ("PROFESSIONAL IDENTIFICATION"), and
(c) http://home.engineering.iastate.edu/~cnchu/ ("URL").
Users who intend to use the Code for commercial purposes will notify Author prior to such commercial use.
Neither the name nor any trademark of the Author may be used to endorse or promote products derived from this software without specific prior written permission.
Users are entirely responsible, to the exclusion of the Author and any other persons, for compliance with (1) regulations set by owners or administrators of employed equipment, (2) licensing terms of any other software, and (3) local, national, and international regulations regarding use, including those regarding import, export, and use of encryption software.
THIS FREE SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR ANY CONTRIBUTOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, EFFECTS OF UNAUTHORIZED OR MALICIOUS NETWORK ACCESS; PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```