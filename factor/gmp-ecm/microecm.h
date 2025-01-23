/*
Copyright (c) 2014, Ben Buhrow and (c) 2022, Jeff Hurchalla.
All rights reserved.
This software is provided "as is," without warranty of any kind,
express or implied.  In no event shall the author or contributors
be held liable for any damages arising in any way from the use of
this software.
The contents of this file are DUAL-LICENSED.  You may modify and/or
redistribute this software according to the terms of one of the
following two licenses (at your option):


LICENSE 1 ("FreeBSD")
---------------------
Copyright (c) 2014, Ben Buhrow and (c) 2022, Jeff Hurchalla.
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


LICENSE 2 (MPL 2.0)
-------------------
Copyright (c) 2014, Ben Buhrow and (c) 2022, Jeff Hurchalla.
All rights reserved.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
*/


#ifndef MICROECM_GETFACTOR_UECM_H_INCLUDED
#define MICROECM_GETFACTOR_UECM_H_INCLUDED

#include <stdint.h>


#ifdef __cplusplus   // C compilers skip this ifdef section
extern "C" {
#endif


// getfactor_uecm() returns 1 if it is unable to find a factor of q64, and
// otherwise it returns a single factor of q64.
//
// if the input is known to have no small factors, set is_arbitrary=0, 
// otherwise, set is_arbitrary=1 and a few curves targetting small factors
// will be run prior to the standard sequence of curves for the input size.
//  
// Prior to your first call of getfactor_uecm(), set *pran = 0  (or set it to
// some other arbitrary value); after that, don't change *pran.
// FYI: *pran is used within microecm.c by a random number generator, and it
// holds the current value of a pseudo random sequence.  Your first assigment
// to *pran seeds the sequence, and after seeding it you don't want to change
// *pran, since that would restart the sequence.

uint64_t getfactor_uecm(uint64_t q64, int is_arbitrary, uint64_t *pran);
void getfactor_uecm_x8_list(uint64_t* q64, uint64_t* f64, uint32_t num_in, uint64_t* pran);
uint64_t getfactor_upm1(uint64_t q64, uint32_t b1);

int prp_uecm(uint64_t n);

#ifdef __cplusplus
}
#endif


#endif
