#pragma once
/*
Copyright (c) 2014, Ben Buhrow
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

*/

#include <stdint.h>

// a Fermat PRP test on 1 64-bit input
int fermat_prp_64x1(uint64_t n);
// a Fermat PRP test on 1 128-bit input
int fermat_prp_128x1(uint64_t* n);
// a Fermat PRP test on 8x 52-bit inputs
uint8_t fermat_prp_52x8(uint64_t* n);
// a Fermat PRP test on 8x 104-bit inputs
uint8_t fermat_prp_104x8(uint64_t* n);
// a Miller-Rabin SPRP test on 8x 52-bit inputs using base 2
uint8_t MR_2sprp_52x8(uint64_t* n);
// a Miller-Rabin SPRP test on 8x 104-bit inputs using base 2
uint8_t MR_2sprp_104x8(uint64_t* n);
// a Miller-Rabin SPRP test on 8x 104-bit inputs using an
// independent arbitrary 52-bit base on each
uint8_t MR_sprp_104x8(uint64_t* n, uint64_t* bases);
// a Miller-Rabin SPRP test on 1 104-bit input using 8x
// different bases
uint8_t MR_sprp_104x8base(uint64_t* n, uint64_t* one, uint64_t* bases);

// test routine for the above
int test_tinyprp();












