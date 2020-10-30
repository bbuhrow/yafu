/*
Copyright (c) 2019, Ben Buhrow
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

/* Elliptic Curve Method: toplevel and stage 1 routines.

Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011,
2012, 2016 Paul Zimmermann, Alexander Kruppa, Cyril Bouvier, David Cleaver.

This file is part of the ECM Library.

The ECM Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The ECM Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the ECM Library; see the file COPYING.LIB.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#include "avx_ecm.h"
#include "soe.h"
#include "arith.h"

//#define D 1155
//#define U 8
//#define R 483
//#define L 2 * U

//#define DEBUG 1

// local functions
void vec_add(vec_monty_t *mdata, ecm_work *work, ecm_pt *Pin, ecm_pt *Pout);
void vec_duplicate(vec_monty_t *mdata, ecm_work *work, vec_bignum_t *insum, vec_bignum_t *indiff, ecm_pt *P);
void next_pt_vec(vec_monty_t *mdata, ecm_work *work, ecm_pt *P, uint64_t c);
void euclid(vec_monty_t *mdata, ecm_work *work, ecm_pt *P, uint64_t c);
void vec_prac(vec_monty_t *mdata, ecm_work *work, ecm_pt *P, uint64_t c);
void vec_build_one_curve(thread_data_t *tdata, mpz_t X, mpz_t Z, mpz_t A, uint64_t sigma);
void build_one_curve_param1(thread_data_t *tdata, mpz_t X, mpz_t Z, mpz_t X2, mpz_t Z2, 
    mpz_t A, uint64_t sigma);
void vec_ecm_stage1(vec_monty_t *mdata, ecm_work *work, ecm_pt *P, base_t b1, base_t *primes, int verbose);
void vec_ecm_stage2_init(ecm_pt *P, vec_monty_t *mdata, ecm_work *work, base_t *primes, int verbose);
void vec_ecm_stage2_pair(ecm_pt *P, vec_monty_t *mdata, ecm_work *work, base_t *primes, int verbose);
int vec_check_factor(mpz_t Z, mpz_t n, mpz_t f);

// GMP-ECM batch stage1
unsigned int compute_s(mpz_t s, uint64_t *primes, uint64_t B1);
int vec_ecm_stage1_batch(mpz_t f, ecm_work *work, ecm_pt *P, vec_bignum_t * A, 
    vec_bignum_t * n, uint64_t B1, mpz_t s, vec_monty_t *mdata);

// support parallel MMIX by requiring the lcg's state to be passed in
uint64_t par_lcg_rand(uint64_t* lcg_state)
{
    // Knuth's MMIX LCG
    *lcg_state = 6364136223846793005ULL * *lcg_state + 1442695040888963407ULL;
    return *lcg_state;
}

#define INV_2_POW_32 2.3283064365386962890625e-10
// Knuth's 64 bit MMIX LCG, using a global 64 bit state variable.
uint32_t spRandp(uint64_t * lcg_state, uint32_t lower, uint32_t upper)
{
    // advance the state of the LCG and return the appropriate result
    *lcg_state = par_lcg_rand(lcg_state);
    return lower + (uint32_t)(
        (double)(upper - lower) * (double)(*lcg_state >> 32) * INV_2_POW_32);
}

// a map of the 479 integers relatively prime to 2310. plus 0, 1, and 2310.
// The map maps into the Pb array of stored elliptic delta points: [d]Q.
// 0 maps to 0 which is used as scratch space.
int rprime_map[2311] = {0, 1,
2, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 3, 0, 0, 0, 4, 0, 5, 0,
0, 0, 6, 0, 0, 0, 0, 0, 7, 0,
8, 0, 0, 0, 0, 0, 9, 0, 0, 0,
10, 0, 11, 0, 0, 0, 12, 0, 0, 0,
0, 0, 13, 0, 0, 0, 0, 0, 14, 0,
15, 0, 0, 0, 0, 0, 16, 0, 0, 0,
17, 0, 18, 0, 0, 0, 0, 0, 19, 0,
0, 0, 20, 0, 0, 0, 0, 0, 21, 0,
0, 0, 0, 0, 0, 0, 22, 0, 0, 0,
23, 0, 24, 0, 0, 0, 25, 0, 26, 0,
0, 0, 27, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 28, 0, 0, 0,
29, 0, 0, 0, 0, 0, 30, 0, 31, 0,
0, 0, 0, 0, 0, 0, 0, 0, 32, 0,
33, 0, 0, 0, 0, 0, 34, 0, 0, 0,
0, 0, 35, 0, 0, 0, 36, 0, 37, 0,
0, 0, 38, 0, 0, 0, 0, 0, 39, 0,
40, 0, 0, 0, 0, 0, 0, 0, 0, 0,
41, 0, 42, 0, 0, 0, 43, 0, 44, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
45, 0, 0, 0, 0, 0, 0, 0, 0, 0,
46, 0, 47, 0, 0, 0, 48, 0, 49, 0,
0, 0, 50, 0, 0, 0, 0, 0, 51, 0,
52, 0, 0, 0, 0, 0, 53, 0, 0, 0,
54, 0, 0, 0, 0, 0, 55, 0, 0, 0,
0, 0, 56, 0, 0, 0, 0, 0, 57, 0,
58, 0, 0, 0, 0, 0, 59, 0, 0, 0,
60, 0, 61, 0, 0, 0, 0, 0, 62, 0,
0, 0, 63, 0, 0, 0, 0, 0, 64, 0,
0, 0, 0, 0, 0, 0, 65, 0, 0, 0,
66, 0, 67, 0, 0, 0, 68, 0, 0, 0,
0, 0, 69, 0, 0, 0, 0, 0, 0, 0,
70, 0, 0, 0, 0, 0, 71, 0, 0, 0,
0, 0, 0, 0, 0, 0, 72, 0, 73, 0,
0, 0, 74, 0, 0, 0, 0, 0, 75, 0,
76, 0, 0, 0, 0, 0, 77, 0, 0, 0,
0, 0, 78, 0, 0, 0, 79, 0, 80, 0,
0, 0, 81, 0, 0, 0, 0, 0, 82, 0,
83, 0, 0, 0, 0, 0, 84, 0, 0, 0,
85, 0, 86, 0, 0, 0, 0, 0, 87, 0,
0, 0, 0, 0, 0, 0, 0, 0, 88, 0,
89, 0, 0, 0, 0, 0, 0, 0, 0, 0,
90, 0, 91, 0, 0, 0, 92, 0, 93, 0,
0, 0, 94, 0, 0, 0, 0, 0, 95, 0,
0, 0, 0, 0, 0, 0, 96, 0, 0, 0,
97, 0, 98, 0, 0, 0, 99, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 100, 0,
101, 0, 0, 0, 0, 0, 102, 0, 0, 0,
103, 0, 104, 0, 0, 0, 0, 0, 105, 0,
0, 0, 106, 0, 0, 0, 0, 0, 107, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
108, 0, 109, 0, 0, 0, 110, 0, 111, 0,
0, 0, 112, 0, 0, 0, 0, 0, 0, 0,
113, 0, 0, 0, 0, 0, 114, 0, 0, 0,
115, 0, 0, 0, 0, 0, 116, 0, 117, 0,
0, 0, 118, 0, 0, 0, 0, 0, 119, 0,
120, 0, 0, 0, 0, 0, 121, 0, 0, 0,
0, 0, 0, 0, 0, 0, 122, 0, 123, 0,
0, 0, 124, 0, 0, 0, 0, 0, 125, 0,
126, 0, 0, 0, 0, 0, 127, 0, 0, 0,
128, 0, 129, 0, 0, 0, 130, 0, 131, 0,
0, 0, 0, 0, 0, 0, 0, 0, 132, 0,
133, 0, 0, 0, 0, 0, 0, 0, 0, 0,
134, 0, 135, 0, 0, 0, 136, 0, 0, 0,
0, 0, 137, 0, 0, 0, 0, 0, 138, 0,
139, 0, 0, 0, 0, 0, 140, 0, 0, 0,
0, 0, 141, 0, 0, 0, 142, 0, 0, 0,
0, 0, 143, 0, 0, 0, 0, 0, 144, 0,
145, 0, 0, 0, 0, 0, 146, 0, 0, 0,
147, 0, 148, 0, 0, 0, 0, 0, 149, 0,
0, 0, 150, 0, 0, 0, 0, 0, 151, 0,
0, 0, 0, 0, 0, 0, 152, 0, 0, 0,
153, 0, 154, 0, 0, 0, 0, 0, 155, 0,
0, 0, 156, 0, 0, 0, 0, 0, 0, 0,
157, 0, 0, 0, 0, 0, 158, 0, 0, 0,
159, 0, 0, 0, 0, 0, 160, 0, 161, 0,
0, 0, 162, 0, 0, 0, 0, 0, 163, 0,
0, 0, 0, 0, 0, 0, 164, 0, 0, 0,
0, 0, 165, 0, 0, 0, 166, 0, 167, 0,
0, 0, 0, 0, 0, 0, 0, 0, 168, 0,
169, 0, 0, 0, 0, 0, 170, 0, 0, 0,
171, 0, 172, 0, 0, 0, 173, 0, 174, 0,
0, 0, 0, 0, 0, 0, 0, 0, 175, 0,
176, 0, 0, 0, 0, 0, 0, 0, 0, 0,
177, 0, 178, 0, 0, 0, 179, 0, 180, 0,
0, 0, 181, 0, 0, 0, 0, 0, 0, 0,
182, 0, 0, 0, 0, 0, 183, 0, 0, 0,
184, 0, 185, 0, 0, 0, 186, 0, 0, 0,
0, 0, 187, 0, 0, 0, 0, 0, 188, 0,
189, 0, 0, 0, 0, 0, 190, 0, 0, 0,
191, 0, 0, 0, 0, 0, 0, 0, 192, 0,
0, 0, 193, 0, 0, 0, 0, 0, 194, 0,
0, 0, 0, 0, 0, 0, 195, 0, 0, 0,
196, 0, 197, 0, 0, 0, 198, 0, 199, 0,
0, 0, 200, 0, 0, 0, 0, 0, 0, 0,
201, 0, 0, 0, 0, 0, 202, 0, 0, 0,
203, 0, 0, 0, 0, 0, 204, 0, 0, 0,
0, 0, 205, 0, 0, 0, 0, 0, 206, 0,
207, 0, 0, 0, 0, 0, 208, 0, 0, 0,
0, 0, 209, 0, 0, 0, 210, 0, 211, 0,
0, 0, 212, 0, 0, 0, 0, 0, 213, 0,
214, 0, 0, 0, 0, 0, 215, 0, 0, 0,
216, 0, 217, 0, 0, 0, 218, 0, 219, 0,
0, 0, 0, 0, 0, 0, 0, 0, 220, 0,
221, 0, 0, 0, 0, 0, 0, 0, 0, 0,
222, 0, 223, 0, 0, 0, 0, 0, 224, 0,
0, 0, 225, 0, 0, 0, 0, 0, 226, 0,
227, 0, 0, 0, 0, 0, 228, 0, 0, 0,
229, 0, 230, 0, 0, 0, 231, 0, 0, 0,
0, 0, 232, 0, 0, 0, 0, 0, 233, 0,
0, 0, 0, 0, 0, 0, 234, 0, 0, 0,
235, 0, 236, 0, 0, 0, 0, 0, 237, 0,
0, 0, 0, 0, 0, 0, 0, 0, 238, 0,
0, 0, 0, 0, 0, 0, 239, 0, 0, 0,
240, 0, 241, 0, 0, 0, 242, 0, 243, 0,
0, 0, 244, 0, 0, 0, 0, 0, 0, 0,
245, 0, 0, 0, 0, 0, 0, 0, 0, 0,
246, 0, 0, 0, 0, 0, 247, 0, 248, 0,
0, 0, 249, 0, 0, 0, 0, 0, 0, 0,
250, 0, 0, 0, 0, 0, 251, 0, 0, 0,
0, 0, 252, 0, 0, 0, 253, 0, 254, 0,
0, 0, 255, 0, 0, 0, 0, 0, 256, 0,
257, 0, 0, 0, 0, 0, 258, 0, 0, 0,
259, 0, 0, 0, 0, 0, 260, 0, 261, 0,
0, 0, 0, 0, 0, 0, 0, 0, 262, 0,
263, 0, 0, 0, 0, 0, 0, 0, 0, 0,
264, 0, 265, 0, 0, 0, 266, 0, 267, 0,
0, 0, 268, 0, 0, 0, 0, 0, 269, 0,
270, 0, 0, 0, 0, 0, 271, 0, 0, 0,
272, 0, 273, 0, 0, 0, 274, 0, 0, 0,
0, 0, 275, 0, 0, 0, 0, 0, 276, 0,
277, 0, 0, 0, 0, 0, 278, 0, 0, 0,
0, 0, 279, 0, 0, 0, 0, 0, 280, 0,
0, 0, 281, 0, 0, 0, 0, 0, 282, 0,
0, 0, 0, 0, 0, 0, 283, 0, 0, 0,
284, 0, 285, 0, 0, 0, 286, 0, 287, 0,
0, 0, 288, 0, 0, 0, 0, 0, 0, 0,
289, 0, 0, 0, 0, 0, 290, 0, 0, 0,
291, 0, 0, 0, 0, 0, 0, 0, 292, 0,
0, 0, 293, 0, 0, 0, 0, 0, 294, 0,
295, 0, 0, 0, 0, 0, 296, 0, 0, 0,
0, 0, 297, 0, 0, 0, 298, 0, 299, 0,
0, 0, 300, 0, 0, 0, 0, 0, 301, 0,
0, 0, 0, 0, 0, 0, 302, 0, 0, 0,
303, 0, 304, 0, 0, 0, 305, 0, 306, 0,
0, 0, 0, 0, 0, 0, 0, 0, 307, 0,
308, 0, 0, 0, 0, 0, 0, 0, 0, 0,
309, 0, 310, 0, 0, 0, 311, 0, 312, 0,
0, 0, 313, 0, 0, 0, 0, 0, 314, 0,
315, 0, 0, 0, 0, 0, 0, 0, 0, 0,
316, 0, 317, 0, 0, 0, 318, 0, 0, 0,
0, 0, 319, 0, 0, 0, 0, 0, 0, 0,
320, 0, 0, 0, 0, 0, 321, 0, 0, 0,
322, 0, 323, 0, 0, 0, 0, 0, 324, 0,
0, 0, 325, 0, 0, 0, 0, 0, 326, 0,
0, 0, 0, 0, 0, 0, 327, 0, 0, 0,
328, 0, 0, 0, 0, 0, 329, 0, 330, 0,
0, 0, 331, 0, 0, 0, 0, 0, 0, 0,
332, 0, 0, 0, 0, 0, 333, 0, 0, 0,
334, 0, 0, 0, 0, 0, 335, 0, 336, 0,
0, 0, 337, 0, 0, 0, 0, 0, 338, 0,
339, 0, 0, 0, 0, 0, 340, 0, 0, 0,
0, 0, 341, 0, 0, 0, 342, 0, 0, 0,
0, 0, 343, 0, 0, 0, 0, 0, 344, 0,
345, 0, 0, 0, 0, 0, 346, 0, 0, 0,
0, 0, 347, 0, 0, 0, 348, 0, 349, 0,
0, 0, 0, 0, 0, 0, 0, 0, 350, 0,
351, 0, 0, 0, 0, 0, 0, 0, 0, 0,
352, 0, 353, 0, 0, 0, 354, 0, 355, 0,
0, 0, 356, 0, 0, 0, 0, 0, 357, 0,
358, 0, 0, 0, 0, 0, 359, 0, 0, 0,
360, 0, 361, 0, 0, 0, 0, 0, 0, 0,
0, 0, 362, 0, 0, 0, 0, 0, 363, 0,
364, 0, 0, 0, 0, 0, 365, 0, 0, 0,
366, 0, 367, 0, 0, 0, 0, 0, 368, 0,
0, 0, 369, 0, 0, 0, 0, 0, 370, 0,
0, 0, 0, 0, 0, 0, 371, 0, 0, 0,
372, 0, 373, 0, 0, 0, 374, 0, 375, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
376, 0, 0, 0, 0, 0, 377, 0, 0, 0,
378, 0, 0, 0, 0, 0, 379, 0, 380, 0,
0, 0, 381, 0, 0, 0, 0, 0, 382, 0,
383, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 384, 0, 0, 0, 385, 0, 386, 0,
0, 0, 387, 0, 0, 0, 0, 0, 0, 0,
388, 0, 0, 0, 0, 0, 389, 0, 0, 0,
390, 0, 391, 0, 0, 0, 392, 0, 393, 0,
0, 0, 0, 0, 0, 0, 0, 0, 394, 0,
395, 0, 0, 0, 0, 0, 0, 0, 0, 0,
396, 0, 0, 0, 0, 0, 397, 0, 398, 0,
0, 0, 399, 0, 0, 0, 0, 0, 400, 0,
401, 0, 0, 0, 0, 0, 402, 0, 0, 0,
403, 0, 404, 0, 0, 0, 405, 0, 0, 0,
0, 0, 406, 0, 0, 0, 0, 0, 407, 0,
408, 0, 0, 0, 0, 0, 409, 0, 0, 0,
410, 0, 411, 0, 0, 0, 0, 0, 0, 0,
0, 0, 412, 0, 0, 0, 0, 0, 413, 0,
0, 0, 0, 0, 0, 0, 414, 0, 0, 0,
0, 0, 415, 0, 0, 0, 416, 0, 417, 0,
0, 0, 418, 0, 0, 0, 0, 0, 0, 0,
419, 0, 0, 0, 0, 0, 420, 0, 0, 0,
421, 0, 0, 0, 0, 0, 422, 0, 423, 0,
0, 0, 424, 0, 0, 0, 0, 0, 425, 0,
426, 0, 0, 0, 0, 0, 427, 0, 0, 0,
0, 0, 428, 0, 0, 0, 0, 0, 429, 0,
0, 0, 430, 0, 0, 0, 0, 0, 431, 0,
432, 0, 0, 0, 0, 0, 433, 0, 0, 0,
434, 0, 435, 0, 0, 0, 436, 0, 437, 0,
0, 0, 0, 0, 0, 0, 0, 0, 438, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
439, 0, 440, 0, 0, 0, 441, 0, 442, 0,
0, 0, 0, 0, 0, 0, 0, 0, 443, 0,
444, 0, 0, 0, 0, 0, 445, 0, 0, 0,
446, 0, 447, 0, 0, 0, 448, 0, 0, 0,
0, 0, 449, 0, 0, 0, 0, 0, 450, 0,
451, 0, 0, 0, 0, 0, 0, 0, 0, 0,
452, 0, 453, 0, 0, 0, 0, 0, 454, 0,
0, 0, 455, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 456, 0, 0, 0,
457, 0, 458, 0, 0, 0, 459, 0, 460, 0,
0, 0, 461, 0, 0, 0, 0, 0, 0, 0,
462, 0, 0, 0, 0, 0, 463, 0, 0, 0,
464, 0, 0, 0, 0, 0, 465, 0, 466, 0,
0, 0, 467, 0, 0, 0, 0, 0, 468, 0,
469, 0, 0, 0, 0, 0, 470, 0, 0, 0,
0, 0, 471, 0, 0, 0, 472, 0, 473, 0,
0, 0, 474, 0, 0, 0, 0, 0, 475, 0,
476, 0, 0, 0, 0, 0, 477, 0, 0, 0,
478, 0, 479, 0, 0, 0, 480, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 481, 482};

int rprime_map_4620[4621] = { 0,
1, 2, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 3, 0, 0, 0, 4, 0, 5, 0,
0, 0, 6, 0, 0, 0, 0, 0, 7, 0,
8, 0, 0, 0, 0, 0, 9, 0, 0, 0,
10, 0, 11, 0, 0, 0, 12, 0, 0, 0,
0, 0, 13, 0, 0, 0, 0, 0, 14, 0,
15, 0, 0, 0, 0, 0, 16, 0, 0, 0,
17, 0, 18, 0, 0, 0, 0, 0, 19, 0,
0, 0, 20, 0, 0, 0, 0, 0, 21, 0,
0, 0, 0, 0, 0, 0, 22, 0, 0, 0,
23, 0, 24, 0, 0, 0, 25, 0, 26, 0,
0, 0, 27, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 28, 0, 0, 0,
29, 0, 0, 0, 0, 0, 30, 0, 31, 0,
0, 0, 0, 0, 0, 0, 0, 0, 32, 0,
33, 0, 0, 0, 0, 0, 34, 0, 0, 0,
0, 0, 35, 0, 0, 0, 36, 0, 37, 0,
0, 0, 38, 0, 0, 0, 0, 0, 39, 0,
40, 0, 0, 0, 0, 0, 0, 0, 0, 0,
41, 0, 42, 0, 0, 0, 43, 0, 44, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
45, 0, 0, 0, 0, 0, 0, 0, 0, 0,
46, 0, 47, 0, 0, 0, 48, 0, 49, 0,
0, 0, 50, 0, 0, 0, 0, 0, 51, 0,
52, 0, 0, 0, 0, 0, 53, 0, 0, 0,
54, 0, 0, 0, 0, 0, 55, 0, 0, 0,
0, 0, 56, 0, 0, 0, 0, 0, 57, 0,
58, 0, 0, 0, 0, 0, 59, 0, 0, 0,
60, 0, 61, 0, 0, 0, 0, 0, 62, 0,
0, 0, 63, 0, 0, 0, 0, 0, 64, 0,
0, 0, 0, 0, 0, 0, 65, 0, 0, 0,
66, 0, 67, 0, 0, 0, 68, 0, 0, 0,
0, 0, 69, 0, 0, 0, 0, 0, 0, 0,
70, 0, 0, 0, 0, 0, 71, 0, 0, 0,
0, 0, 0, 0, 0, 0, 72, 0, 73, 0,
0, 0, 74, 0, 0, 0, 0, 0, 75, 0,
76, 0, 0, 0, 0, 0, 77, 0, 0, 0,
0, 0, 78, 0, 0, 0, 79, 0, 80, 0,
0, 0, 81, 0, 0, 0, 0, 0, 82, 0,
83, 0, 0, 0, 0, 0, 84, 0, 0, 0,
85, 0, 86, 0, 0, 0, 0, 0, 87, 0,
0, 0, 0, 0, 0, 0, 0, 0, 88, 0,
89, 0, 0, 0, 0, 0, 0, 0, 0, 0,
90, 0, 91, 0, 0, 0, 92, 0, 93, 0,
0, 0, 94, 0, 0, 0, 0, 0, 95, 0,
0, 0, 0, 0, 0, 0, 96, 0, 0, 0,
97, 0, 98, 0, 0, 0, 99, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 100, 0,
101, 0, 0, 0, 0, 0, 102, 0, 0, 0,
103, 0, 104, 0, 0, 0, 0, 0, 105, 0,
0, 0, 106, 0, 0, 0, 0, 0, 107, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
108, 0, 109, 0, 0, 0, 110, 0, 111, 0,
0, 0, 112, 0, 0, 0, 0, 0, 0, 0,
113, 0, 0, 0, 0, 0, 114, 0, 0, 0,
115, 0, 0, 0, 0, 0, 116, 0, 117, 0,
0, 0, 118, 0, 0, 0, 0, 0, 119, 0,
120, 0, 0, 0, 0, 0, 121, 0, 0, 0,
0, 0, 0, 0, 0, 0, 122, 0, 123, 0,
0, 0, 124, 0, 0, 0, 0, 0, 125, 0,
126, 0, 0, 0, 0, 0, 127, 0, 0, 0,
128, 0, 129, 0, 0, 0, 130, 0, 131, 0,
0, 0, 0, 0, 0, 0, 0, 0, 132, 0,
133, 0, 0, 0, 0, 0, 0, 0, 0, 0,
134, 0, 135, 0, 0, 0, 136, 0, 0, 0,
0, 0, 137, 0, 0, 0, 0, 0, 138, 0,
139, 0, 0, 0, 0, 0, 140, 0, 0, 0,
0, 0, 141, 0, 0, 0, 142, 0, 0, 0,
0, 0, 143, 0, 0, 0, 0, 0, 144, 0,
145, 0, 0, 0, 0, 0, 146, 0, 0, 0,
147, 0, 148, 0, 0, 0, 0, 0, 149, 0,
0, 0, 150, 0, 0, 0, 0, 0, 151, 0,
0, 0, 0, 0, 0, 0, 152, 0, 0, 0,
153, 0, 154, 0, 0, 0, 0, 0, 155, 0,
0, 0, 156, 0, 0, 0, 0, 0, 0, 0,
157, 0, 0, 0, 0, 0, 158, 0, 0, 0,
159, 0, 0, 0, 0, 0, 160, 0, 161, 0,
0, 0, 162, 0, 0, 0, 0, 0, 163, 0,
0, 0, 0, 0, 0, 0, 164, 0, 0, 0,
0, 0, 165, 0, 0, 0, 166, 0, 167, 0,
0, 0, 0, 0, 0, 0, 0, 0, 168, 0,
169, 0, 0, 0, 0, 0, 170, 0, 0, 0,
171, 0, 172, 0, 0, 0, 173, 0, 174, 0,
0, 0, 0, 0, 0, 0, 0, 0, 175, 0,
176, 0, 0, 0, 0, 0, 0, 0, 0, 0,
177, 0, 178, 0, 0, 0, 179, 0, 180, 0,
0, 0, 181, 0, 0, 0, 0, 0, 0, 0,
182, 0, 0, 0, 0, 0, 183, 0, 0, 0,
184, 0, 185, 0, 0, 0, 186, 0, 0, 0,
0, 0, 187, 0, 0, 0, 0, 0, 188, 0,
189, 0, 0, 0, 0, 0, 190, 0, 0, 0,
191, 0, 0, 0, 0, 0, 0, 0, 192, 0,
0, 0, 193, 0, 0, 0, 0, 0, 194, 0,
0, 0, 0, 0, 0, 0, 195, 0, 0, 0,
196, 0, 197, 0, 0, 0, 198, 0, 199, 0,
0, 0, 200, 0, 0, 0, 0, 0, 0, 0,
201, 0, 0, 0, 0, 0, 202, 0, 0, 0,
203, 0, 0, 0, 0, 0, 204, 0, 0, 0,
0, 0, 205, 0, 0, 0, 0, 0, 206, 0,
207, 0, 0, 0, 0, 0, 208, 0, 0, 0,
0, 0, 209, 0, 0, 0, 210, 0, 211, 0,
0, 0, 212, 0, 0, 0, 0, 0, 213, 0,
214, 0, 0, 0, 0, 0, 215, 0, 0, 0,
216, 0, 217, 0, 0, 0, 218, 0, 219, 0,
0, 0, 0, 0, 0, 0, 0, 0, 220, 0,
221, 0, 0, 0, 0, 0, 0, 0, 0, 0,
222, 0, 223, 0, 0, 0, 0, 0, 224, 0,
0, 0, 225, 0, 0, 0, 0, 0, 226, 0,
227, 0, 0, 0, 0, 0, 228, 0, 0, 0,
229, 0, 230, 0, 0, 0, 231, 0, 0, 0,
0, 0, 232, 0, 0, 0, 0, 0, 233, 0,
0, 0, 0, 0, 0, 0, 234, 0, 0, 0,
235, 0, 236, 0, 0, 0, 0, 0, 237, 0,
0, 0, 0, 0, 0, 0, 0, 0, 238, 0,
0, 0, 0, 0, 0, 0, 239, 0, 0, 0,
240, 0, 241, 0, 0, 0, 242, 0, 243, 0,
0, 0, 244, 0, 0, 0, 0, 0, 0, 0,
245, 0, 0, 0, 0, 0, 0, 0, 0, 0,
246, 0, 0, 0, 0, 0, 247, 0, 248, 0,
0, 0, 249, 0, 0, 0, 0, 0, 0, 0,
250, 0, 0, 0, 0, 0, 251, 0, 0, 0,
0, 0, 252, 0, 0, 0, 253, 0, 254, 0,
0, 0, 255, 0, 0, 0, 0, 0, 256, 0,
257, 0, 0, 0, 0, 0, 258, 0, 0, 0,
259, 0, 0, 0, 0, 0, 260, 0, 261, 0,
0, 0, 0, 0, 0, 0, 0, 0, 262, 0,
263, 0, 0, 0, 0, 0, 0, 0, 0, 0,
264, 0, 265, 0, 0, 0, 266, 0, 267, 0,
0, 0, 268, 0, 0, 0, 0, 0, 269, 0,
270, 0, 0, 0, 0, 0, 271, 0, 0, 0,
272, 0, 273, 0, 0, 0, 274, 0, 0, 0,
0, 0, 275, 0, 0, 0, 0, 0, 276, 0,
277, 0, 0, 0, 0, 0, 278, 0, 0, 0,
0, 0, 279, 0, 0, 0, 0, 0, 280, 0,
0, 0, 281, 0, 0, 0, 0, 0, 282, 0,
0, 0, 0, 0, 0, 0, 283, 0, 0, 0,
284, 0, 285, 0, 0, 0, 286, 0, 287, 0,
0, 0, 288, 0, 0, 0, 0, 0, 0, 0,
289, 0, 0, 0, 0, 0, 290, 0, 0, 0,
291, 0, 0, 0, 0, 0, 0, 0, 292, 0,
0, 0, 293, 0, 0, 0, 0, 0, 294, 0,
295, 0, 0, 0, 0, 0, 296, 0, 0, 0,
0, 0, 297, 0, 0, 0, 298, 0, 299, 0,
0, 0, 300, 0, 0, 0, 0, 0, 301, 0,
0, 0, 0, 0, 0, 0, 302, 0, 0, 0,
303, 0, 304, 0, 0, 0, 305, 0, 306, 0,
0, 0, 0, 0, 0, 0, 0, 0, 307, 0,
308, 0, 0, 0, 0, 0, 0, 0, 0, 0,
309, 0, 310, 0, 0, 0, 311, 0, 312, 0,
0, 0, 313, 0, 0, 0, 0, 0, 314, 0,
315, 0, 0, 0, 0, 0, 0, 0, 0, 0,
316, 0, 317, 0, 0, 0, 318, 0, 0, 0,
0, 0, 319, 0, 0, 0, 0, 0, 0, 0,
320, 0, 0, 0, 0, 0, 321, 0, 0, 0,
322, 0, 323, 0, 0, 0, 0, 0, 324, 0,
0, 0, 325, 0, 0, 0, 0, 0, 326, 0,
0, 0, 0, 0, 0, 0, 327, 0, 0, 0,
328, 0, 0, 0, 0, 0, 329, 0, 330, 0,
0, 0, 331, 0, 0, 0, 0, 0, 0, 0,
332, 0, 0, 0, 0, 0, 333, 0, 0, 0,
334, 0, 0, 0, 0, 0, 335, 0, 336, 0,
0, 0, 337, 0, 0, 0, 0, 0, 338, 0,
339, 0, 0, 0, 0, 0, 340, 0, 0, 0,
0, 0, 341, 0, 0, 0, 342, 0, 0, 0,
0, 0, 343, 0, 0, 0, 0, 0, 344, 0,
345, 0, 0, 0, 0, 0, 346, 0, 0, 0,
0, 0, 347, 0, 0, 0, 348, 0, 349, 0,
0, 0, 0, 0, 0, 0, 0, 0, 350, 0,
351, 0, 0, 0, 0, 0, 0, 0, 0, 0,
352, 0, 353, 0, 0, 0, 354, 0, 355, 0,
0, 0, 356, 0, 0, 0, 0, 0, 357, 0,
358, 0, 0, 0, 0, 0, 359, 0, 0, 0,
360, 0, 361, 0, 0, 0, 0, 0, 0, 0,
0, 0, 362, 0, 0, 0, 0, 0, 363, 0,
364, 0, 0, 0, 0, 0, 365, 0, 0, 0,
366, 0, 367, 0, 0, 0, 0, 0, 368, 0,
0, 0, 369, 0, 0, 0, 0, 0, 370, 0,
0, 0, 0, 0, 0, 0, 371, 0, 0, 0,
372, 0, 373, 0, 0, 0, 374, 0, 375, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
376, 0, 0, 0, 0, 0, 377, 0, 0, 0,
378, 0, 0, 0, 0, 0, 379, 0, 380, 0,
0, 0, 381, 0, 0, 0, 0, 0, 382, 0,
383, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 384, 0, 0, 0, 385, 0, 386, 0,
0, 0, 387, 0, 0, 0, 0, 0, 0, 0,
388, 0, 0, 0, 0, 0, 389, 0, 0, 0,
390, 0, 391, 0, 0, 0, 392, 0, 393, 0,
0, 0, 0, 0, 0, 0, 0, 0, 394, 0,
395, 0, 0, 0, 0, 0, 0, 0, 0, 0,
396, 0, 0, 0, 0, 0, 397, 0, 398, 0,
0, 0, 399, 0, 0, 0, 0, 0, 400, 0,
401, 0, 0, 0, 0, 0, 402, 0, 0, 0,
403, 0, 404, 0, 0, 0, 405, 0, 0, 0,
0, 0, 406, 0, 0, 0, 0, 0, 407, 0,
408, 0, 0, 0, 0, 0, 409, 0, 0, 0,
410, 0, 411, 0, 0, 0, 0, 0, 0, 0,
0, 0, 412, 0, 0, 0, 0, 0, 413, 0,
0, 0, 0, 0, 0, 0, 414, 0, 0, 0,
0, 0, 415, 0, 0, 0, 416, 0, 417, 0,
0, 0, 418, 0, 0, 0, 0, 0, 0, 0,
419, 0, 0, 0, 0, 0, 420, 0, 0, 0,
421, 0, 0, 0, 0, 0, 422, 0, 423, 0,
0, 0, 424, 0, 0, 0, 0, 0, 425, 0,
426, 0, 0, 0, 0, 0, 427, 0, 0, 0,
0, 0, 428, 0, 0, 0, 0, 0, 429, 0,
0, 0, 430, 0, 0, 0, 0, 0, 431, 0,
432, 0, 0, 0, 0, 0, 433, 0, 0, 0,
434, 0, 435, 0, 0, 0, 436, 0, 437, 0,
0, 0, 0, 0, 0, 0, 0, 0, 438, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
439, 0, 440, 0, 0, 0, 441, 0, 442, 0,
0, 0, 0, 0, 0, 0, 0, 0, 443, 0,
444, 0, 0, 0, 0, 0, 445, 0, 0, 0,
446, 0, 447, 0, 0, 0, 448, 0, 0, 0,
0, 0, 449, 0, 0, 0, 0, 0, 450, 0,
451, 0, 0, 0, 0, 0, 0, 0, 0, 0,
452, 0, 453, 0, 0, 0, 0, 0, 454, 0,
0, 0, 455, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 456, 0, 0, 0,
457, 0, 458, 0, 0, 0, 459, 0, 460, 0,
0, 0, 461, 0, 0, 0, 0, 0, 0, 0,
462, 0, 0, 0, 0, 0, 463, 0, 0, 0,
464, 0, 0, 0, 0, 0, 465, 0, 466, 0,
0, 0, 467, 0, 0, 0, 0, 0, 468, 0,
469, 0, 0, 0, 0, 0, 470, 0, 0, 0,
0, 0, 471, 0, 0, 0, 472, 0, 473, 0,
0, 0, 474, 0, 0, 0, 0, 0, 475, 0,
476, 0, 0, 0, 0, 0, 477, 0, 0, 0,
478, 0, 479, 0, 0, 0, 480, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 481, 0,
482, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 483, 0, 0, 0, 484, 0, 485, 0,
0, 0, 486, 0, 0, 0, 0, 0, 487, 0,
488, 0, 0, 0, 0, 0, 489, 0, 0, 0,
490, 0, 491, 0, 0, 0, 492, 0, 0, 0,
0, 0, 493, 0, 0, 0, 0, 0, 494, 0,
495, 0, 0, 0, 0, 0, 496, 0, 0, 0,
497, 0, 498, 0, 0, 0, 0, 0, 499, 0,
0, 0, 500, 0, 0, 0, 0, 0, 501, 0,
0, 0, 0, 0, 0, 0, 502, 0, 0, 0,
503, 0, 504, 0, 0, 0, 505, 0, 506, 0,
0, 0, 507, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 508, 0, 0, 0,
509, 0, 0, 0, 0, 0, 510, 0, 511, 0,
0, 0, 0, 0, 0, 0, 0, 0, 512, 0,
513, 0, 0, 0, 0, 0, 514, 0, 0, 0,
0, 0, 515, 0, 0, 0, 516, 0, 517, 0,
0, 0, 518, 0, 0, 0, 0, 0, 519, 0,
520, 0, 0, 0, 0, 0, 0, 0, 0, 0,
521, 0, 522, 0, 0, 0, 523, 0, 524, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
525, 0, 0, 0, 0, 0, 0, 0, 0, 0,
526, 0, 527, 0, 0, 0, 528, 0, 529, 0,
0, 0, 530, 0, 0, 0, 0, 0, 531, 0,
532, 0, 0, 0, 0, 0, 533, 0, 0, 0,
534, 0, 0, 0, 0, 0, 535, 0, 0, 0,
0, 0, 536, 0, 0, 0, 0, 0, 537, 0,
538, 0, 0, 0, 0, 0, 539, 0, 0, 0,
540, 0, 541, 0, 0, 0, 0, 0, 542, 0,
0, 0, 543, 0, 0, 0, 0, 0, 544, 0,
0, 0, 0, 0, 0, 0, 545, 0, 0, 0,
546, 0, 547, 0, 0, 0, 548, 0, 0, 0,
0, 0, 549, 0, 0, 0, 0, 0, 0, 0,
550, 0, 0, 0, 0, 0, 551, 0, 0, 0,
0, 0, 0, 0, 0, 0, 552, 0, 553, 0,
0, 0, 554, 0, 0, 0, 0, 0, 555, 0,
556, 0, 0, 0, 0, 0, 557, 0, 0, 0,
0, 0, 558, 0, 0, 0, 559, 0, 560, 0,
0, 0, 561, 0, 0, 0, 0, 0, 562, 0,
563, 0, 0, 0, 0, 0, 564, 0, 0, 0,
565, 0, 566, 0, 0, 0, 0, 0, 567, 0,
0, 0, 0, 0, 0, 0, 0, 0, 568, 0,
569, 0, 0, 0, 0, 0, 0, 0, 0, 0,
570, 0, 571, 0, 0, 0, 572, 0, 573, 0,
0, 0, 574, 0, 0, 0, 0, 0, 575, 0,
0, 0, 0, 0, 0, 0, 576, 0, 0, 0,
577, 0, 578, 0, 0, 0, 579, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 580, 0,
581, 0, 0, 0, 0, 0, 582, 0, 0, 0,
583, 0, 584, 0, 0, 0, 0, 0, 585, 0,
0, 0, 586, 0, 0, 0, 0, 0, 587, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
588, 0, 589, 0, 0, 0, 590, 0, 591, 0,
0, 0, 592, 0, 0, 0, 0, 0, 0, 0,
593, 0, 0, 0, 0, 0, 594, 0, 0, 0,
595, 0, 0, 0, 0, 0, 596, 0, 597, 0,
0, 0, 598, 0, 0, 0, 0, 0, 599, 0,
600, 0, 0, 0, 0, 0, 601, 0, 0, 0,
0, 0, 0, 0, 0, 0, 602, 0, 603, 0,
0, 0, 604, 0, 0, 0, 0, 0, 605, 0,
606, 0, 0, 0, 0, 0, 607, 0, 0, 0,
608, 0, 609, 0, 0, 0, 610, 0, 611, 0,
0, 0, 0, 0, 0, 0, 0, 0, 612, 0,
613, 0, 0, 0, 0, 0, 0, 0, 0, 0,
614, 0, 615, 0, 0, 0, 616, 0, 0, 0,
0, 0, 617, 0, 0, 0, 0, 0, 618, 0,
619, 0, 0, 0, 0, 0, 620, 0, 0, 0,
0, 0, 621, 0, 0, 0, 622, 0, 0, 0,
0, 0, 623, 0, 0, 0, 0, 0, 624, 0,
625, 0, 0, 0, 0, 0, 626, 0, 0, 0,
627, 0, 628, 0, 0, 0, 0, 0, 629, 0,
0, 0, 630, 0, 0, 0, 0, 0, 631, 0,
0, 0, 0, 0, 0, 0, 632, 0, 0, 0,
633, 0, 634, 0, 0, 0, 0, 0, 635, 0,
0, 0, 636, 0, 0, 0, 0, 0, 0, 0,
637, 0, 0, 0, 0, 0, 638, 0, 0, 0,
639, 0, 0, 0, 0, 0, 640, 0, 641, 0,
0, 0, 642, 0, 0, 0, 0, 0, 643, 0,
0, 0, 0, 0, 0, 0, 644, 0, 0, 0,
0, 0, 645, 0, 0, 0, 646, 0, 647, 0,
0, 0, 0, 0, 0, 0, 0, 0, 648, 0,
649, 0, 0, 0, 0, 0, 650, 0, 0, 0,
651, 0, 652, 0, 0, 0, 653, 0, 654, 0,
0, 0, 0, 0, 0, 0, 0, 0, 655, 0,
656, 0, 0, 0, 0, 0, 0, 0, 0, 0,
657, 0, 658, 0, 0, 0, 659, 0, 660, 0,
0, 0, 661, 0, 0, 0, 0, 0, 0, 0,
662, 0, 0, 0, 0, 0, 663, 0, 0, 0,
664, 0, 665, 0, 0, 0, 666, 0, 0, 0,
0, 0, 667, 0, 0, 0, 0, 0, 668, 0,
669, 0, 0, 0, 0, 0, 670, 0, 0, 0,
671, 0, 0, 0, 0, 0, 0, 0, 672, 0,
0, 0, 673, 0, 0, 0, 0, 0, 674, 0,
0, 0, 0, 0, 0, 0, 675, 0, 0, 0,
676, 0, 677, 0, 0, 0, 678, 0, 679, 0,
0, 0, 680, 0, 0, 0, 0, 0, 0, 0,
681, 0, 0, 0, 0, 0, 682, 0, 0, 0,
683, 0, 0, 0, 0, 0, 684, 0, 0, 0,
0, 0, 685, 0, 0, 0, 0, 0, 686, 0,
687, 0, 0, 0, 0, 0, 688, 0, 0, 0,
0, 0, 689, 0, 0, 0, 690, 0, 691, 0,
0, 0, 692, 0, 0, 0, 0, 0, 693, 0,
694, 0, 0, 0, 0, 0, 695, 0, 0, 0,
696, 0, 697, 0, 0, 0, 698, 0, 699, 0,
0, 0, 0, 0, 0, 0, 0, 0, 700, 0,
701, 0, 0, 0, 0, 0, 0, 0, 0, 0,
702, 0, 703, 0, 0, 0, 0, 0, 704, 0,
0, 0, 705, 0, 0, 0, 0, 0, 706, 0,
707, 0, 0, 0, 0, 0, 708, 0, 0, 0,
709, 0, 710, 0, 0, 0, 711, 0, 0, 0,
0, 0, 712, 0, 0, 0, 0, 0, 713, 0,
0, 0, 0, 0, 0, 0, 714, 0, 0, 0,
715, 0, 716, 0, 0, 0, 0, 0, 717, 0,
0, 0, 0, 0, 0, 0, 0, 0, 718, 0,
0, 0, 0, 0, 0, 0, 719, 0, 0, 0,
720, 0, 721, 0, 0, 0, 722, 0, 723, 0,
0, 0, 724, 0, 0, 0, 0, 0, 0, 0,
725, 0, 0, 0, 0, 0, 0, 0, 0, 0,
726, 0, 0, 0, 0, 0, 727, 0, 728, 0,
0, 0, 729, 0, 0, 0, 0, 0, 0, 0,
730, 0, 0, 0, 0, 0, 731, 0, 0, 0,
0, 0, 732, 0, 0, 0, 733, 0, 734, 0,
0, 0, 735, 0, 0, 0, 0, 0, 736, 0,
737, 0, 0, 0, 0, 0, 738, 0, 0, 0,
739, 0, 0, 0, 0, 0, 740, 0, 741, 0,
0, 0, 0, 0, 0, 0, 0, 0, 742, 0,
743, 0, 0, 0, 0, 0, 0, 0, 0, 0,
744, 0, 745, 0, 0, 0, 746, 0, 747, 0,
0, 0, 748, 0, 0, 0, 0, 0, 749, 0,
750, 0, 0, 0, 0, 0, 751, 0, 0, 0,
752, 0, 753, 0, 0, 0, 754, 0, 0, 0,
0, 0, 755, 0, 0, 0, 0, 0, 756, 0,
757, 0, 0, 0, 0, 0, 758, 0, 0, 0,
0, 0, 759, 0, 0, 0, 0, 0, 760, 0,
0, 0, 761, 0, 0, 0, 0, 0, 762, 0,
0, 0, 0, 0, 0, 0, 763, 0, 0, 0,
764, 0, 765, 0, 0, 0, 766, 0, 767, 0,
0, 0, 768, 0, 0, 0, 0, 0, 0, 0,
769, 0, 0, 0, 0, 0, 770, 0, 0, 0,
771, 0, 0, 0, 0, 0, 0, 0, 772, 0,
0, 0, 773, 0, 0, 0, 0, 0, 774, 0,
775, 0, 0, 0, 0, 0, 776, 0, 0, 0,
0, 0, 777, 0, 0, 0, 778, 0, 779, 0,
0, 0, 780, 0, 0, 0, 0, 0, 781, 0,
0, 0, 0, 0, 0, 0, 782, 0, 0, 0,
783, 0, 784, 0, 0, 0, 785, 0, 786, 0,
0, 0, 0, 0, 0, 0, 0, 0, 787, 0,
788, 0, 0, 0, 0, 0, 0, 0, 0, 0,
789, 0, 790, 0, 0, 0, 791, 0, 792, 0,
0, 0, 793, 0, 0, 0, 0, 0, 794, 0,
795, 0, 0, 0, 0, 0, 0, 0, 0, 0,
796, 0, 797, 0, 0, 0, 798, 0, 0, 0,
0, 0, 799, 0, 0, 0, 0, 0, 0, 0,
800, 0, 0, 0, 0, 0, 801, 0, 0, 0,
802, 0, 803, 0, 0, 0, 0, 0, 804, 0,
0, 0, 805, 0, 0, 0, 0, 0, 806, 0,
0, 0, 0, 0, 0, 0, 807, 0, 0, 0,
808, 0, 0, 0, 0, 0, 809, 0, 810, 0,
0, 0, 811, 0, 0, 0, 0, 0, 0, 0,
812, 0, 0, 0, 0, 0, 813, 0, 0, 0,
814, 0, 0, 0, 0, 0, 815, 0, 816, 0,
0, 0, 817, 0, 0, 0, 0, 0, 818, 0,
819, 0, 0, 0, 0, 0, 820, 0, 0, 0,
0, 0, 821, 0, 0, 0, 822, 0, 0, 0,
0, 0, 823, 0, 0, 0, 0, 0, 824, 0,
825, 0, 0, 0, 0, 0, 826, 0, 0, 0,
0, 0, 827, 0, 0, 0, 828, 0, 829, 0,
0, 0, 0, 0, 0, 0, 0, 0, 830, 0,
831, 0, 0, 0, 0, 0, 0, 0, 0, 0,
832, 0, 833, 0, 0, 0, 834, 0, 835, 0,
0, 0, 836, 0, 0, 0, 0, 0, 837, 0,
838, 0, 0, 0, 0, 0, 839, 0, 0, 0,
840, 0, 841, 0, 0, 0, 0, 0, 0, 0,
0, 0, 842, 0, 0, 0, 0, 0, 843, 0,
844, 0, 0, 0, 0, 0, 845, 0, 0, 0,
846, 0, 847, 0, 0, 0, 0, 0, 848, 0,
0, 0, 849, 0, 0, 0, 0, 0, 850, 0,
0, 0, 0, 0, 0, 0, 851, 0, 0, 0,
852, 0, 853, 0, 0, 0, 854, 0, 855, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
856, 0, 0, 0, 0, 0, 857, 0, 0, 0,
858, 0, 0, 0, 0, 0, 859, 0, 860, 0,
0, 0, 861, 0, 0, 0, 0, 0, 862, 0,
863, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 864, 0, 0, 0, 865, 0, 866, 0,
0, 0, 867, 0, 0, 0, 0, 0, 0, 0,
868, 0, 0, 0, 0, 0, 869, 0, 0, 0,
870, 0, 871, 0, 0, 0, 872, 0, 873, 0,
0, 0, 0, 0, 0, 0, 0, 0, 874, 0,
875, 0, 0, 0, 0, 0, 0, 0, 0, 0,
876, 0, 0, 0, 0, 0, 877, 0, 878, 0,
0, 0, 879, 0, 0, 0, 0, 0, 880, 0,
881, 0, 0, 0, 0, 0, 882, 0, 0, 0,
883, 0, 884, 0, 0, 0, 885, 0, 0, 0,
0, 0, 886, 0, 0, 0, 0, 0, 887, 0,
888, 0, 0, 0, 0, 0, 889, 0, 0, 0,
890, 0, 891, 0, 0, 0, 0, 0, 0, 0,
0, 0, 892, 0, 0, 0, 0, 0, 893, 0,
0, 0, 0, 0, 0, 0, 894, 0, 0, 0,
0, 0, 895, 0, 0, 0, 896, 0, 897, 0,
0, 0, 898, 0, 0, 0, 0, 0, 0, 0,
899, 0, 0, 0, 0, 0, 900, 0, 0, 0,
901, 0, 0, 0, 0, 0, 902, 0, 903, 0,
0, 0, 904, 0, 0, 0, 0, 0, 905, 0,
906, 0, 0, 0, 0, 0, 907, 0, 0, 0,
0, 0, 908, 0, 0, 0, 0, 0, 909, 0,
0, 0, 910, 0, 0, 0, 0, 0, 911, 0,
912, 0, 0, 0, 0, 0, 913, 0, 0, 0,
914, 0, 915, 0, 0, 0, 916, 0, 917, 0,
0, 0, 0, 0, 0, 0, 0, 0, 918, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
919, 0, 920, 0, 0, 0, 921, 0, 922, 0,
0, 0, 0, 0, 0, 0, 0, 0, 923, 0,
924, 0, 0, 0, 0, 0, 925, 0, 0, 0,
926, 0, 927, 0, 0, 0, 928, 0, 0, 0,
0, 0, 929, 0, 0, 0, 0, 0, 930, 0,
931, 0, 0, 0, 0, 0, 0, 0, 0, 0,
932, 0, 933, 0, 0, 0, 0, 0, 934, 0,
0, 0, 935, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 936, 0, 0, 0,
937, 0, 938, 0, 0, 0, 939, 0, 940, 0,
0, 0, 941, 0, 0, 0, 0, 0, 0, 0,
942, 0, 0, 0, 0, 0, 943, 0, 0, 0,
944, 0, 0, 0, 0, 0, 945, 0, 946, 0,
0, 0, 947, 0, 0, 0, 0, 0, 948, 0,
949, 0, 0, 0, 0, 0, 950, 0, 0, 0,
0, 0, 951, 0, 0, 0, 952, 0, 953, 0,
0, 0, 954, 0, 0, 0, 0, 0, 955, 0,
956, 0, 0, 0, 0, 0, 957, 0, 0, 0,
958, 0, 959, 0, 0, 0, 960, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 961, 962};

int mulcnt[256];
int sqrcnt[256];

#include "threadpool.h"

void vec_ecm_sync(void *vptr);
void vec_ecm_dispatch(void *vptr);
void vec_ecm_build_curve_work_fcn(void *vptr);
void vec_ecm_stage1_work_fcn(void *vptr);
void vec_ecm_stage2_work_fcn(void *vptr);


void
mpres_get_z(mpz_t R, const mpz_t S, mpz_t modulus, vec_monty_t *mdata)
{
    // not fast, but gets the input S out of Mrep.
    mpz_set(mdata->gmp_t1, S);
    //ecm_redc_basecase(R, mdata->gmp_t1, modulus);
    mpz_mul(mdata->gmp_t2, mdata->nhat, S);
    mpz_tdiv_r_2exp(mdata->gmp_t2, mdata->gmp_t2, MAXBITS);
    mpz_mul(mdata->gmp_t2, modulus, mdata->gmp_t2);
    mpz_add(mdata->gmp_t2, S, mdata->gmp_t2);
    mpz_tdiv_q_2exp(mdata->gmp_t2, mdata->gmp_t2, MAXBITS);
    mpz_mod(R, mdata->gmp_t2, modulus);
}

/* R <- S / 2^n mod modulus. Does not need to be fast. */
void mpres_div_2exp(mpz_t R, mpz_t S, const unsigned int n, mpz_t modulus)
{
    int i;

    if (n == 0)
    {
        mpz_set(R, S);
        return;
    }

    if (mpz_odd_p(S))
    {
        mpz_add(R, S, modulus);
        mpz_tdiv_q_2exp(R, R, 1);
    }
    else
        mpz_tdiv_q_2exp(R, S, 1);

    for (i = n; i > 1; i--)
    {
        if (mpz_odd_p(R))
        {
            mpz_add(R, R, modulus);
        }
        mpz_tdiv_q_2exp(R, R, 1);
    }

}

/* R <- S * 2^k mod modulus */
void mpres_mul_2exp(mpz_t R, mpz_t S, const unsigned long k, mpz_t modulus)
{
    mpz_mul_2exp(R, S, k);
    /* This is the same for all methods: just reduce with original modulus */
    mpz_mod(R, R, modulus);
}

void mpres_add_ui(mpz_t R, mpz_t S, const unsigned long n, mpz_t modulus)
{
    mpz_set_ui(R, n);
    mpz_mul_2exp(R, R, MAXBITS);
    mpz_add(R, R, S);
    mpz_mod(R, R, modulus);
}


void vec_ecm_sync(void *vptr)
{
    tpool_t *tpdata = (tpool_t *)vptr;
    thread_data_t *udata = (thread_data_t *)tpdata->user_data;
    uint32_t tid = tpdata->tindex;

    if (udata[tid].ecm_phase == 0)
    {
        udata[tid].phase_done = 1;
    }
    else if (udata->ecm_phase == 1)
    {
        udata[tid].phase_done = 1;
    }
    else if (udata->ecm_phase == 2)
    {
        udata[tid].phase_done = 1;
    }
    else if (udata->ecm_phase == 3)
    {
        udata[tid].phase_done = 1;
    }
    
    return;
}

void vec_ecm_dispatch(void *vptr)
{
    tpool_t *tpdata = (tpool_t *)vptr;
    thread_data_t *udata = (thread_data_t *)tpdata->user_data;
    uint32_t tid = tpdata->tindex;
    
    if (udata[tid].ecm_phase == 0)
    {        
        if (udata[tid].phase_done == 0)
        {
            tpdata->work_fcn_id = 0;
        }
        else
        {
            tpdata->work_fcn_id = tpdata->num_work_fcn;
        }
    }
    else if (udata[tid].ecm_phase == 1)
    {
        if (udata[tid].phase_done == 0)
        {
            tpdata->work_fcn_id = 1;
        }
        else
        {
            tpdata->work_fcn_id = tpdata->num_work_fcn;
        }
    }
    else if (udata[tid].ecm_phase == 2)
    {
        if (udata[tid].phase_done == 0)
        {
            tpdata->work_fcn_id = 2;
        }
        else
        {
            tpdata->work_fcn_id = tpdata->num_work_fcn;
        }
    }
    else if (udata[tid].ecm_phase == 3)
    {
        if (udata[tid].phase_done == 0)
        {
            tpdata->work_fcn_id = 3;
        }
        else
        {
            tpdata->work_fcn_id = tpdata->num_work_fcn;
        }
    }

    return;
}

void vec_ecm_stage1_work_fcn(void *vptr)
{
    tpool_t *tpdata = (tpool_t *)vptr;
    thread_data_t *udata = (thread_data_t *)tpdata->user_data;
    uint32_t tid = tpdata->tindex;

    vec_ecm_stage1(udata[tid].mdata, udata[tid].work, udata[tid].P, 
        udata[tid].b1, NULL, udata[tid].verbose && (tid == 0));

    return;
}

void vec_ecm_stage2_init_work_fcn(void *vptr)
{
    tpool_t *tpdata = (tpool_t *)vptr;
    thread_data_t *udata = (thread_data_t *)tpdata->user_data;
    uint32_t tid = tpdata->tindex;

    vec_ecm_stage2_init(udata[tid].P, udata[tid].mdata, 
        udata[tid].work, NULL, udata[tid].verbose && (tid == 0));

    return;
}

void vec_ecm_stage2_work_fcn(void *vptr)
{
    tpool_t *tpdata = (tpool_t *)vptr;
    thread_data_t *udata = (thread_data_t *)tpdata->user_data;
    uint32_t tid = tpdata->tindex;

    vec_ecm_stage2_pair(udata[tid].P, udata[tid].mdata,
        udata[tid].work, NULL, udata[tid].verbose && (tid == 0));

    return;
}

//#define PARAM1

void vec_ecm_build_curve_work_fcn(void *vptr)
{
    tpool_t *tpdata = (tpool_t *)vptr;
    thread_data_t *tdata = (thread_data_t *)tpdata->user_data;
    uint32_t tid = tpdata->tindex;
    mpz_t X, Z, A; // , t, n, one, r;
#ifdef PARAM1
    mpz_t X2, Z2;
#endif
    int i;

    mpz_init(X);
    mpz_init(Z);
    mpz_init(A);
#ifdef PARAM1
    mpz_init(X2);
    mpz_init(Z2);
#endif

	for (i = 0; i < VECLEN; i++)
    {
        int j;
#ifdef PARAM1
        build_one_curve_param1(&tdata[tid], X, Z, X2, Z2, A, 0); // 65953669); // 36875098);
        insert_mpz_to_vec(tdata[tid].P->X, X, i);
        insert_mpz_to_vec(tdata[tid].P->Z, Z, i);
        insert_mpz_to_vec(tdata[tid].work->pt2.X, X2, i);
        insert_mpz_to_vec(tdata[tid].work->pt2.Z, Z2, i);
        insert_mpz_to_vec(tdata[tid].work->s, A, i);
#else
        // 8689346476060549ULL
        vec_build_one_curve(&tdata[tid], X, Z, A, 0); // 6710676431370252287 + i);

        insert_mpz_to_vec(tdata[tid].P->X, X, i);
        insert_mpz_to_vec(tdata[tid].P->Z, Z, i);
        insert_mpz_to_vec(tdata[tid].work->s, A, i);
#endif
        tdata[tid].sigma[i] = tdata[tid].work->sigma;
    }

    //print_vechexbignum(tdata[tid].P->X, "POINT x: ");
    //print_vechexbignum(tdata[tid].P->Z, "POINT z: ");
    //print_vechexbignum(tdata[tid].work->s, "POINT a: ");

    tdata[tid].work->diff1->size = tdata[tid].work->n->size;
    tdata[tid].work->diff2->size = tdata[tid].work->n->size;
    tdata[tid].work->sum1->size  = tdata[tid].work->n->size;
    tdata[tid].work->sum2->size  = tdata[tid].work->n->size;
    tdata[tid].work->pt1.X->size = tdata[tid].work->n->size;
    tdata[tid].work->pt1.Z->size = tdata[tid].work->n->size;
    tdata[tid].work->pt2.X->size = tdata[tid].work->n->size;
    tdata[tid].work->pt2.Z->size = tdata[tid].work->n->size;
    tdata[tid].work->tt1->size   = tdata[tid].work->n->size;
    tdata[tid].work->tt2->size   = tdata[tid].work->n->size;
    tdata[tid].work->tt3->size   = tdata[tid].work->n->size;
    tdata[tid].work->tt4->size   = tdata[tid].work->n->size;

    mpz_clear(X);
    mpz_clear(Z);
    mpz_clear(A);

#ifdef PARAM1
    mpz_clear(X2);
    mpz_clear(Z2);
#endif

    return;
}

void vec_ecm_work_init(ecm_work *work)
{
	int i, j, m, sz = 4 * VECLEN * NWORDS;
	uint32_t U = work->U;
	uint32_t L = work->L;
	uint32_t D = work->D;
	uint32_t R = work->R;

	work->stg1Add = 0;
	work->stg1Doub = 0;

	work->diff1 = vecInit();
	work->diff2 = vecInit();
	work->sum1 = vecInit();
	work->sum2 = vecInit();
	work->tt1 = vecInit();
	work->tt2 = vecInit();
	work->tt3 = vecInit();
	work->tt4 = vecInit();
	work->tt5 = vecInit();
	work->s = vecInit();
	work->n = vecInit();

	work->Paprod = (vec_bignum_t **)malloc((2 * L) * sizeof(vec_bignum_t *));
	work->Pa = (ecm_pt *)malloc((2 * L) * sizeof(ecm_pt));
	for (j = 0; j < (2 * L); j++)
	{
		work->Paprod[j] = vecInit();
		vec_ecm_pt_init(&work->Pa[j]);
	}

	work->Pad = (ecm_pt *)malloc(sizeof(ecm_pt));
    vec_ecm_pt_init(work->Pad);

	// build an array to hold values of f(b).
	// will need to be of size U*R, allowed values will be
	// up to a multiple of U of the residues mod D
	work->Pb = (ecm_pt *)malloc(U * (R + 1) * sizeof(ecm_pt));
	work->Pbprod = (vec_bignum_t **)malloc(U * (R + 1) * sizeof(vec_bignum_t *));
	for (j = 0; j < U * (R + 1); j++)
	{
        vec_ecm_pt_init(&work->Pb[j]);
		work->Pbprod[j] = vecInit();
	}

	work->marks = (uint8_t *)malloc(2 * U * D * sizeof(uint8_t));
	work->nmarks = (uint8_t *)malloc(U * D * sizeof(uint8_t));
	work->map = (uint32_t *)calloc(U * (D + 1) + 3, sizeof(uint32_t));

	work->map[0] = 0;
	work->map[1] = 1;
	work->map[2] = 2;
	m = 3;
	for (i = 0; i < U; i++)
	{
		if (i == 0)
			j = 3;
		else
			j = 1;

		for (; j < D; j++)
		{
			if (spGCD(j, D) == 1)
			{
				work->map[i*D + j] = m++;
			}
			else
			{
				work->map[i*D + j] = 0;
			}
		}

		if (i == 0)
			work->map[i*D + j] = m++;

	}

	work->Qmap = (uint32_t *)malloc(D * sizeof(uint32_t));
	work->Qrmap = (uint32_t *)malloc(D * sizeof(uint32_t));

	for (j = 0, i = 0; i < D; i++)
	{
		if (spGCD(i, D) == 1)
		{
			work->Qmap[i] = j;
			work->Qrmap[j++] = i;
		}
		else
		{
			work->Qmap[i] = (uint32_t)-1;
		}
	}

	for (i = j; i < D; i++)
	{
		work->Qrmap[i] = (uint32_t)-1;
	}

	work->Q = (Queue_t **)malloc(2 * j * sizeof(Queue_t *));
	for (i = 0; i < 2 * j; i++)
	{
		work->Q[i] = newQueue(D);
	}


	/*
	j = 0;
	printf("residue map of %u\n", D);
	for (i = 0; i <= U * D; i++, j++)
	{
		printf("%u, ", work->map[i]);
		if (j % 10 == 0)
			printf("\n");
	}
	printf("\n");
	*/

	vec_ecm_pt_init(&work->pt1);
	vec_ecm_pt_init(&work->pt2);
	vec_ecm_pt_init(&work->pt3);
	vec_ecm_pt_init(&work->pt4);
	vec_ecm_pt_init(&work->pt5);

	work->stg2acc = vecInit();

	return;
}

void vec_ecm_work_free(ecm_work *work)
{
	int i;
	uint32_t U = work->U;
	uint32_t L = work->L;
	uint32_t D = work->D;
	uint32_t R = work->R;

	vecFree(work->diff1);
	vecFree(work->diff2);
	vecFree(work->sum1);
	vecFree(work->sum2);
	vecFree(work->tt1);
	vecFree(work->tt2);
	vecFree(work->tt3);
	vecFree(work->tt4);
	vecFree(work->tt5);
	vecFree(work->s);
	vecFree(work->n);
	vec_ecm_pt_free(&work->pt1);
	vec_ecm_pt_free(&work->pt2);
	vec_ecm_pt_free(&work->pt3);
	vec_ecm_pt_free(&work->pt4);
	vec_ecm_pt_free(&work->pt5);

	for (i = 0; i < (2 * L); i++)
	{
		vec_ecm_pt_free(&work->Pa[i]);
		vecFree(work->Paprod[i]);
	}

	vec_ecm_pt_free(work->Pad);
	free(work->Pa);
	free(work->Pad);
	free(work->marks);
	free(work->nmarks);
	free(work->map);
	
	for (i = 0; i < U * (R + 1); i++)
	{
		vec_ecm_pt_free(&work->Pb[i]);
		vecFree(work->Pbprod[i]);
	}

	for (i = 0; i < 2 * (R - 3); i++)
	{
		clearQueue(work->Q[i]);
		free(work->Q[i]);
	}
	free(work->Q);
	free(work->Qmap);
	free(work->Qrmap);

	free(work->Pbprod);
	free(work->Paprod);
	free(work->Pb);
	vecFree(work->stg2acc);

	return;
}

void vec_ecm_pt_init(ecm_pt *pt)
{
	pt->X = vecInit();
	pt->Z = vecInit();
}

void vec_ecm_pt_free(ecm_pt *pt)
{
	vecFree(pt->X);
	vecFree(pt->Z);
}

void vec_add(vec_monty_t *mdata, ecm_work *work, ecm_pt *Pin, ecm_pt *Pout)
{
    // compute:
    //x+ = z- * [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
    //z+ = x- * [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
    // where:
    //x- = original x
    //z- = original z
    // given the sums and differences of the original points (stored in work structure).

    vecmulmod_ptr(work->diff1, work->sum2, work->tt1, work->n, work->tt4, mdata);	//U
    vecmulmod_ptr(work->sum1, work->diff2, work->tt2, work->n, work->tt4, mdata);	//V

    vecaddsubmod_ptr(work->tt1, work->tt2, work->tt3, work->tt4, work->n);
    vecsqrmod_ptr(work->tt3, work->tt1, work->n, work->tt5, mdata);					//(U + V)^2
    vecsqrmod_ptr(work->tt4, work->tt2, work->n, work->tt5, mdata);					//(U - V)^2

    // choosing the initial point Pz0 = 1 means that z_p-q = 1 and this mul isn't necessary...
    // but that involves a different way to initialize curves, so for now
    // we can't assume Z=1
	if (Pin->X == Pout->X)
	{
		base_t *swap;
		vecmulmod_ptr(work->tt1, Pin->Z, Pout->Z, work->n, work->tt4, mdata);		//Z * (U + V)^2
		vecmulmod_ptr(work->tt2, Pin->X, Pout->X, work->n, work->tt4, mdata);		//x * (U - V)^2
		swap = Pout->Z->data;
		Pout->Z->data = Pout->X->data;
		Pout->X->data = swap;
	}
	else
	{
		vecmulmod_ptr(work->tt1, Pin->Z, Pout->X, work->n, work->tt4, mdata);		//Z * (U + V)^2
		vecmulmod_ptr(work->tt2, Pin->X, Pout->Z, work->n, work->tt4, mdata);		//x * (U - V)^2
	}
	work->stg1Add++;
    return;
}

void vec_duplicate(vec_monty_t *mdata, ecm_work *work, vec_bignum_t *insum, vec_bignum_t *indiff, ecm_pt *P)
{
    vecsqrmod_ptr(indiff, work->tt1, work->n, work->tt4, mdata);			    // V=(x1 - z1)^2
    vecsqrmod_ptr(insum, work->tt2, work->n, work->tt4, mdata);			        // U=(x1 + z1)^2
    vecmulmod_ptr(work->tt1, work->tt2, P->X, work->n, work->tt4, mdata);       // x=U*V

    vecsubmod_ptr(work->tt2, work->tt1, work->tt3, work->n);	                // w = U-V
    vecmulmod_ptr(work->tt3, work->s, work->tt2, work->n, work->tt4, mdata);    // t = (A+2)/4 * w
    vecaddmod_ptr(work->tt2, work->tt1, work->tt2, work->n);                    // t = t + V
    vecmulmod_ptr(work->tt2, work->tt3, P->Z, work->n, work->tt4, mdata);       // Z = t*w
	work->stg1Doub++;
    return;
}

double vec_getEcost(uint64_t d, uint64_t e)
{
	int doub = 0, add = 0;

	while (d > 0)
	{
		if ((e / 2) < d)
		{
			d = e - d;
		}
		else if ((d < (e / 4)) && ((e & 1) == 0))
		{
			e = e / 2;
			doub++;
			add++;
		}
		else
		{
			e = e - d;
			add++;
		}

	}
	return (doub + add) * 2 * 0.75 + add * 4 + doub * 3;
}

int * getEseq(uint64_t d, uint64_t e)
{
	uint64_t target = e;
	int doub = 0, add = 0;
	int seqLen = 64;
	int *seq = (int *)malloc(seqLen * sizeof(int));
	int it = 1;
	seq[0] = 0;

	d = e - d;

	while (d > 0)
	{
		if ((e / 2) < d)
		{
#ifdef DEBUG
			printf("[3]: d = %u\n", e - d);
#endif
			seq[it] = 3;
			d = e - d;
		}
		else if ((d < (e / 4)) && ((e & 1) == 0))
		{
#ifdef DEBUG
			printf("[1]: d = %u, e = %u\n", d, e);
			printf("[1]: d = %u, e = %u\n", d, e - d);
#endif
			e = e / 2;
			doub++;
			add++;
			seq[it] = 1;
		}
		else
		{
#ifdef DEBUG
			printf("[2]: d = %u, e = %u\n", d, e);
#endif
			e = e - d;
			add++;
			seq[it] = 2;
		}

		it++;
		if (it > seqLen)
		{
			seqLen *= 2;
			seq = (int *)realloc(seq, seqLen * sizeof(int));
		}
	}

	// reverse the sequence so it's constructive starting from 0,1
	seqLen = it;
	for (it = 0; it < seqLen / 2; it++)
	{
		uint32_t tmp = seq[it];
		seq[it] = seq[seqLen - it - 1];
		seq[seqLen - it - 1] = tmp;
	}

	// for fun.
#ifdef DEBUG
	if (1)
#else
	if (0)
#endif
	{
		uint32_t d = 0;
		uint32_t x = 1;
		it = 0;

		printf("target is %u, sequence length is %u\n", target, seqLen);
		while (seq[it] != 0)
		{
			if (seq[it] == 1)
			{
				printf("[1]: double to get %u and add %d to get %u\n", x * 2, x - d, x + x - d);
				x *= 2;
			}
			else if (seq[it] == 2)
			{
				x += d;
				printf("[2]: add %d to get %u\n", d, x);
			}
			else if (seq[it] == 3)
			{
				d = x - d;
				printf("[3]: d is now %u\n", d);
			}
			it++;
		}
	}

	return seq;
}

#define ADD 6.0
#define DUP 5.0

static double
lucas_cost(uint64_t n, double v)
{
	uint64_t d, e, r;
	double c; /* cost */

	d = n;
	r = (uint64_t)((double)d * v + 0.5);
	if (r >= n)
		return (ADD * (double)n);
	d = n - r;
	e = 2 * r - n;
	c = DUP + ADD; /* initial duplicate and final addition */
	while (d != e)
	{
		if (d < e)
		{
			r = d;
			d = e;
			e = r;
		}
		if (d - e <= e / 4 && ((d + e) % 3) == 0)
		{ /* condition 1 */
			d = (2 * d - e) / 3;
			e = (e - d) / 2;
			c += 3.0 * ADD; /* 3 additions */
		}
		else if (d - e <= e / 4 && (d - e) % 6 == 0)
		{ /* condition 2 */
			d = (d - e) / 2;
			c += ADD + DUP; /* one addition, one duplicate */
		}
		else if ((d + 3) / 4 <= e)
		{ /* condition 3 */
			d -= e;
			c += ADD; /* one addition */
		}
		else if ((d + e) % 2 == 0)
		{ /* condition 4 */
			d = (d - e) / 2;
			c += ADD + DUP; /* one addition, one duplicate */
		}
		/* now d+e is odd */
		else if (d % 2 == 0)
		{ /* condition 5 */
			d /= 2;
			c += ADD + DUP; /* one addition, one duplicate */
		}
		/* now d is odd and e is even */
		else if (d % 3 == 0)
		{ /* condition 6 */
			d = d / 3 - e;
			c += 3.0 * ADD + DUP; /* three additions, one duplicate */
		}
		else if ((d + e) % 3 == 0)
		{ /* condition 7 */
			d = (d - 2 * e) / 3;
			c += 3.0 * ADD + DUP; /* three additions, one duplicate */
		}
		else if ((d - e) % 3 == 0)
		{ /* condition 8 */
			d = (d - e) / 3;
			c += 3.0 * ADD + DUP; /* three additions, one duplicate */
		}
		else /* necessarily e is even: catches all cases */
		{ /* condition 9 */
			e /= 2;
			c += ADD + DUP; /* one addition, one duplicate */
		}
	}

	return c;
}

void vec_prac(vec_monty_t *mdata, ecm_work *work, ecm_pt *P, uint64_t c)
{
	uint64_t d, e, r;
	double cmin, cost;
	int i;
	vec_bignum_t *s1, *s2, *d1, *d2;
	base_t *sw_x, *sw_z;

#define NV 10  
	/* 1/val[0] = the golden ratio (1+sqrt(5))/2, and 1/val[i] for i>0
	   is the real number whose continued fraction expansion is all 1s
	   except for a 2 in i+1-st place */
	static double val[NV] =
	{ 0.61803398874989485, 0.72360679774997897, 0.58017872829546410,
	  0.63283980608870629, 0.61242994950949500, 0.62018198080741576,
	  0.61721461653440386, 0.61834711965622806, 0.61791440652881789,
	  0.61807966846989581 };

	/* chooses the best value of v */
	for (d = 0, cmin = ADD * (double)c; d < NV; d++)
	{
		cost = lucas_cost(c, val[d]);
		if (cost < cmin)
		{
			cmin = cost;
			i = d;
		}
	}

	d = c;
	r = (uint64_t)((double)d * val[i] + 0.5);

	s1 = work->sum1;
	s2 = work->sum2;
	d1 = work->diff1;
	d2 = work->diff2;

	/* first iteration always begins by Condition 3, then a swap */
	d = c - r;
	e = 2 * r - c;

	// mpres_set(xB, xA, n);
	// mpres_set(zB, zA, n); /* B=A */
	// mpres_set(xC, xA, n);
	// mpres_set(zC, zA, n); /* C=A */
	// duplicate(xA, zA, xA, zA, n, b, u, v, w); /* A = 2*A */

	// the first one is always a doubling
	// point1 is [1]P
	vecCopy(P->X, work->pt1.X);
	vecCopy(P->Z, work->pt1.Z);
	vecCopy(P->X, work->pt2.X);
	vecCopy(P->Z, work->pt2.Z);
	vecCopy(P->X, work->pt3.X);
	vecCopy(P->Z, work->pt3.Z);
	vecsubmod_ptr(work->pt1.X, work->pt1.Z, d1, work->n);
	vecaddmod_ptr(work->pt1.X, work->pt1.Z, s1, work->n);

	// point2 is [2]P
	vec_duplicate(mdata, work, s1, d1, &work->pt1);

	while (d != e)
	{
		if (d < e)
		{
			r = d;
			d = e;
			e = r;
			//mpres_swap(xA, xB, n);
			//mpres_swap(zA, zB, n);
			sw_x = work->pt1.X->data;
			sw_z = work->pt1.Z->data;
			work->pt1.X->data = work->pt2.X->data;
			work->pt1.Z->data = work->pt2.Z->data;
			work->pt2.X->data = sw_x;
			work->pt2.Z->data = sw_z;
		}
		/* do the first line of Table 4 whose condition qualifies */
		if (d - e <= e / 4 && ((d + e) % 3) == 0)
		{ /* condition 1 */
			d = (2 * d - e) / 3;
			e = (e - d) / 2;
            
            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s1, d1, work->n);
            vecaddsubmod_ptr(work->pt2.X, work->pt2.Z, s2, d2, work->n);

			vec_add(mdata, work, &work->pt3, &work->pt4); // T = A + B (C)

            vecaddsubmod_ptr(work->pt4.X, work->pt4.Z, s1, d1, work->n);
            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s2, d2, work->n);

			vec_add(mdata, work, &work->pt2, &work->pt5); // T2 = T + A (B)

            vecaddsubmod_ptr(work->pt2.X, work->pt2.Z, s1, d1, work->n);
            vecaddsubmod_ptr(work->pt4.X, work->pt4.Z, s2, d2, work->n);

			vec_add(mdata, work, &work->pt1, &work->pt2); // B = B + T (A)

			//add3(xT, zT, xA, zA, xB, zB, xC, zC, n, u, v, w); /* T = f(A,B,C) */
			//add3(xT2, zT2, xT, zT, xA, zA, xB, zB, n, u, v, w); /* T2 = f(T,A,B) */
			//add3(xB, zB, xB, zB, xT, zT, xA, zA, n, u, v, w); /* B = f(B,T,A) */
			//mpres_swap(xA, xT2, n);
			//mpres_swap(zA, zT2, n); /* swap A and T2 */

			sw_x = work->pt1.X->data;
			sw_z = work->pt1.Z->data;
			work->pt1.X->data = work->pt5.X->data;
			work->pt1.Z->data = work->pt5.Z->data;
			work->pt5.X->data = sw_x;
			work->pt5.Z->data = sw_z;

		}
		else if (d - e <= e / 4 && (d - e) % 6 == 0)
		{ /* condition 2 */
			d = (d - e) / 2;
			
            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s1, d1, work->n);
            vecaddsubmod_ptr(work->pt2.X, work->pt2.Z, s2, d2, work->n);

			vec_add(mdata, work, &work->pt3, &work->pt2);		// B = A + B (C)
			vec_duplicate(mdata, work, s1, d1, &work->pt1);		// A = 2A

			//add3(xB, zB, xA, zA, xB, zB, xC, zC, n, u, v, w); /* B = f(A,B,C) */
			//duplicate(xA, zA, xA, zA, n, b, u, v, w); /* A = 2*A */

		}
		else if ((d + 3) / 4 <= e)
		{ /* condition 3 */
			d -= e;
			
            vecaddsubmod_ptr(work->pt2.X, work->pt2.Z, s1, d1, work->n);
            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s2, d2, work->n);

			vec_add(mdata, work, &work->pt3, &work->pt4);		// T = B + A (C)
			//add3(xT, zT, xB, zB, xA, zA, xC, zC, n, u, v, w); /* T = f(B,A,C) */
			
			/* circular permutation (B,T,C) */
			//tmp = xB;
			//xB = xT;
			//xT = xC;
			//xC = tmp;
			//tmp = zB;
			//zB = zT;
			//zT = zC;
			//zC = tmp;

			sw_x = work->pt2.X->data;
			sw_z = work->pt2.Z->data;
			work->pt2.X->data = work->pt4.X->data;
			work->pt2.Z->data = work->pt4.Z->data;
			work->pt4.X->data = work->pt3.X->data;
			work->pt4.Z->data = work->pt3.Z->data;
			work->pt3.X->data = sw_x;
			work->pt3.Z->data = sw_z;

		}
		else if ((d + e) % 2 == 0)
		{ /* condition 4 */
			d = (d - e) / 2;

            vecaddsubmod_ptr(work->pt2.X, work->pt2.Z, s1, d1, work->n);
            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s2, d2, work->n);

			vec_add(mdata, work, &work->pt3, &work->pt2);		// B = B + A (C)
			vec_duplicate(mdata, work, s2, d2, &work->pt1);		// A = 2A
			
			//add3(xB, zB, xB, zB, xA, zA, xC, zC, n, u, v, w); /* B = f(B,A,C) */
			//duplicate(xA, zA, xA, zA, n, b, u, v, w); /* A = 2*A */
		}
		/* now d+e is odd */
		else if (d % 2 == 0)
		{ /* condition 5 */
			d /= 2;
			
            vecaddsubmod_ptr(work->pt3.X, work->pt3.Z, s1, d1, work->n);
            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s2, d2, work->n);

			vec_add(mdata, work, &work->pt2, &work->pt3);		// C = C + A (B)
			vec_duplicate(mdata, work, s2, d2, &work->pt1);		// A = 2A

			//add3(xC, zC, xC, zC, xA, zA, xB, zB, n, u, v, w); /* C = f(C,A,B) */
			//duplicate(xA, zA, xA, zA, n, b, u, v, w); /* A = 2*A */
		}
		/* now d is odd, e is even */
		else if (d % 3 == 0)
		{ /* condition 6 */
			d = d / 3 - e;

            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s1, d1, work->n);

			vec_duplicate(mdata, work, s1, d1, &work->pt4);		// T = 2A

            vecaddsubmod_ptr(work->pt2.X, work->pt2.Z, s2, d2, work->n);

			vec_add(mdata, work, &work->pt3, &work->pt5);		// T2 = A + B (C)
			
            vecaddsubmod_ptr(work->pt4.X, work->pt4.Z, s1, d1, work->n);
            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s2, d2, work->n);

			vec_add(mdata, work, &work->pt1, &work->pt1);		// A = T + A (A)

            vecaddsubmod_ptr(work->pt5.X, work->pt5.Z, s2, d2, work->n);

			vec_add(mdata, work, &work->pt3, &work->pt4);		// T = T + T2 (C)

			//duplicate(xT, zT, xA, zA, n, b, u, v, w); /* T = 2*A */
			//add3(xT2, zT2, xA, zA, xB, zB, xC, zC, n, u, v, w); /* T2 = f(A,B,C) */
			//add3(xA, zA, xT, zT, xA, zA, xA, zA, n, u, v, w); /* A = f(T,A,A) */
			//add3(xT, zT, xT, zT, xT2, zT2, xC, zC, n, u, v, w); /* T = f(T,T2,C) */

			/* circular permutation (C,B,T) */
			//tmp = xC;
			//xC = xB;
			//xB = xT;
			//xT = tmp;
			//tmp = zC;
			//zC = zB;
			//zB = zT;
			//zT = tmp;

			sw_x = work->pt3.X->data;
			sw_z = work->pt3.Z->data;
			work->pt3.X->data = work->pt2.X->data;
			work->pt3.Z->data = work->pt2.Z->data;
			work->pt2.X->data = work->pt4.X->data;
			work->pt2.Z->data = work->pt4.Z->data;
			work->pt4.X->data = sw_x;
			work->pt4.Z->data = sw_z;

		}
		else if ((d + e) % 3 == 0)
		{ /* condition 7 */
			d = (d - 2 * e) / 3;

            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s1, d1, work->n);
            vecaddsubmod_ptr(work->pt2.X, work->pt2.Z, s2, d2, work->n);

			vec_add(mdata, work, &work->pt3, &work->pt4);		// T = A + B (C)

            vecaddsubmod_ptr(work->pt4.X, work->pt4.Z, s1, d1, work->n);
            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s2, d2, work->n);

			vec_add(mdata, work, &work->pt2, &work->pt2);		// B = T + A (B)

			vec_duplicate(mdata, work, s2, d2, &work->pt4);		// T = 2A

            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s1, d1, work->n);
            vecaddsubmod_ptr(work->pt4.X, work->pt4.Z, s2, d2, work->n);

			vec_add(mdata, work, &work->pt1, &work->pt1);		// A = A + T (A) = 3A

			//add3(xT, zT, xA, zA, xB, zB, xC, zC, n, u, v, w); /* T = f(A,B,C) */
			//add3(xB, zB, xT, zT, xA, zA, xB, zB, n, u, v, w); /* B = f(T,A,B) */
			//duplicate(xT, zT, xA, zA, n, b, u, v, w);
			//add3(xA, zA, xA, zA, xT, zT, xA, zA, n, u, v, w); /* A = 3*A */
		}
		else if ((d - e) % 3 == 0)
		{ /* condition 8 */
			d = (d - e) / 3;
            // in a run with B1=100M, this condition only executed a few times.

            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s1, d1, work->n);
            vecaddsubmod_ptr(work->pt2.X, work->pt2.Z, s2, d2, work->n);

			vec_add(mdata, work, &work->pt3, &work->pt4);		// T = A + B (C)

            vecaddsubmod_ptr(work->pt3.X, work->pt3.Z, s1, d1, work->n);
            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s2, d2, work->n);

			vec_add(mdata, work, &work->pt2, &work->pt3);		// C = C + A (B)

			//add3(xT, zT, xA, zA, xB, zB, xC, zC, n, u, v, w); /* T = f(A,B,C) */
			//add3(xC, zC, xC, zC, xA, zA, xB, zB, n, u, v, w); /* C = f(A,C,B) */
			//mpres_swap(xB, xT, n);
			//mpres_swap(zB, zT, n); /* swap B and T */
			sw_x = work->pt2.X->data;
			sw_z = work->pt2.Z->data;
			work->pt2.X->data = work->pt4.X->data;
			work->pt2.Z->data = work->pt4.Z->data;
			work->pt4.X->data = sw_x;
			work->pt4.Z->data = sw_z;

            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s2, d2, work->n);

			vec_duplicate(mdata, work, s2, d2, &work->pt4);		// T = 2A

            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s1, d1, work->n);
            vecaddsubmod_ptr(work->pt4.X, work->pt4.Z, s2, d2, work->n);

			vec_add(mdata, work, &work->pt1, &work->pt1);		// A = A + T (A) = 3A

			//duplicate(xT, zT, xA, zA, n, b, u, v, w);
			//add3(xA, zA, xA, zA, xT, zT, xA, zA, n, u, v, w); /* A = 3*A */
		}
		else /* necessarily e is even here */
		{ /* condition 9 */
			e /= 2;
            // in a run with B1=100M, this condition never executed.

            vecaddsubmod_ptr(work->pt3.X, work->pt3.Z, s1, d1, work->n);
            vecaddsubmod_ptr(work->pt2.X, work->pt2.Z, s2, d2, work->n);

			vec_add(mdata, work, &work->pt1, &work->pt3);		// C = C + B (A)
			vec_duplicate(mdata, work, s2, d2, &work->pt2);		// B = 2B

			//add3(xC, zC, xC, zC, xB, zB, xA, zA, n, u, v, w); /* C = f(C,B,A) */
			//duplicate(xB, zB, xB, zB, n, b, u, v, w); /* B = 2*B */
		}
	}

	vecsubmod_ptr(work->pt1.X, work->pt1.Z, d1, work->n);
	vecaddmod_ptr(work->pt1.X, work->pt1.Z, s1, work->n);
	vecsubmod_ptr(work->pt2.X, work->pt2.Z, d2, work->n);
	vecaddmod_ptr(work->pt2.X, work->pt2.Z, s2, work->n);

	vec_add(mdata, work, &work->pt3, P);		// A = A + B (C)

	//add3(xA, zA, xA, zA, xB, zB, xC, zC, n, u, v, w);

	if (d != 1)
	{
		printf("problem: d != 1\n");
	}

	return;

}

void euclid(vec_monty_t *mdata, ecm_work *work, ecm_pt *P, uint64_t c)
{
	uint64_t startd;
	uint64_t numd = 10;
	uint64_t bestd = 0;
	int *seq;
	double best;
	double cost;
	int i;
	uint64_t d, x;
	vec_bignum_t *x1, *z1, *x2, *z2, *x3, *z3, *x4, *z4, *s1, *s2, *d1, *d2;
	base_t *sw_x, *sw_z;

	if (c == 1)
	{
		return;
	}
	
	// thank you gmp-ecm
	static double val[10] =
	{ 0.61803398874989485, 0.72360679774997897, 0.58017872829546410,
	  0.63283980608870629, 0.61242994950949500, 0.62018198080741576,
	  0.61721461653440386, 0.61834711965622806, 0.61791440652881789,
	  0.61807966846989581 };

	i = 0;
	best = 999999.0;

	//while (numd > 0)
	while (i < numd)
	{
		//uint64_t d = startd - i;
		uint64_t d = (uint64_t)((double)c * val[i]);
		uint64_t e = c;

		if (spGCD(d, e) != 1)
		{
			i++;
			continue;
		}

		cost = vec_getEcost(d, e);
		if (cost < best)
		{
			best = cost;
			bestd = d;
		}
		i++;
		numd--;
	}

	if ((c & 1) == 0)
	{
		printf("input should not be even\n");
	}

	seq = getEseq(bestd, c);

	// now follow the sequence
	x1 = work->pt1.X;
	z1 = work->pt1.Z;
	x2 = work->pt2.X;
	z2 = work->pt2.Z;
	x3 = work->pt3.X;
	z3 = work->pt3.Z;
	x4 = work->pt4.X;
	z4 = work->pt4.Z;
	s1 = work->sum1;
	s2 = work->sum2;
	d1 = work->diff1;
	d2 = work->diff2;

	// the first one is always a doubling
	// point1 is [1]P
	vecCopy(P->X, x1);
	vecCopy(P->Z, z1);
	vecCopy(P->X, x3);
	vecCopy(P->Z, z3);
	vecsubmod_ptr(P->X, P->Z, d1, work->n);
	vecaddmod_ptr(P->X, P->Z, s1, work->n);

	d = 1;
	x = 2;

	// point2 is [2]P
	vec_duplicate(mdata, work, s1, d1, &work->pt2);

#ifdef DEBUG
	printf("target is %lu\n", c);
	printf("pt1 holds [1]P\n");
	print_vechex(work->pt1.X->data, 0, NWORDS, "");
	print_vechex(work->pt1.Z->data, 0, NWORDS, "");
	printf("pt2 holds [2]P\n");
	print_vechex(work->pt2.X->data, 0, NWORDS, "");
	print_vechex(work->pt2.Z->data, 0, NWORDS, "");
#endif

	if (c == 2)
	{
		free(seq);
		vecCopy(x2, P->X);
		vecCopy(z2, P->Z);
		return;
	}
	
	startd = 0;
	i = 2;
	while (seq[i])
	{
		// both add and dup require the sum and difference 
		// of the X and Z coords of the two points we are tracking.
		vecaddsubmod_ptr(x2, z2, s2, d2, work->n);
		vecaddsubmod_ptr(x1, z1, s1, d1, work->n);

		if (seq[i] == 1)
		{
			// add point1 to point2, store in point1
			// double point2
			vec_add(mdata, work, &work->pt3, &work->pt1);
			vec_duplicate(mdata, work, s2, d2, &work->pt2);

			if (seq[i + 1] == 3)
			{
				d = x + d;
				x *= 2;
				i++;
			}
			else if (seq[i + 1] == 2)
			{
				d = x - d;
				x *= 2;

				// swap the "initial" point held in pt3 with the 
				// output point held in pt1
				sw_x = x3->data;
				sw_z = z3->data;
				x3->data = x1->data;
				z3->data = z1->data;
				x1->data = sw_x;
				z1->data = sw_z;
			}
			else
			{
				d = x + d;
				x *= 2;
			}

#ifdef DEBUG
			printf("[1]: pt1 holds [%d]P = ", d);
			print_vechex(work->pt1.X->data, 0, NWORDS, "");
			print_vechex(work->pt1.Z->data, 0, NWORDS, "");
			printf("[1]: pt2 holds [%d]P = ", x);
			print_vechex(work->pt2.X->data, 0, NWORDS, "");
			print_vechex(work->pt2.Z->data, 0, NWORDS, "");
#endif

		}
		else if(seq[i] == 2)
		{
			// add point1 (implicit in add/sub inputs) to point2, store in point2.
			// what is now in point 2 will either be copied to 
			// the new initial point (pt3) or to point 1.
			vecCopy(x2, x4);
			vecCopy(z2, z4);
			
#ifdef DEBUG
			printf("[2]: tmp holds [%d]P = ", x);
			print_vechex(work->pt4.X->data, 0, NWORDS, "");
			print_vechex(work->pt4.Z->data, 0, NWORDS, "");
#endif

			vec_add(mdata, work, &work->pt3, &work->pt2);
			x = x + d;

#ifdef DEBUG
			printf("[2]: pt2 holds [%d]P = ", x);
			print_vechex(work->pt2.X->data, 0, NWORDS, "");
			print_vechex(work->pt2.Z->data, 0, NWORDS, "");
#endif

			if (seq[i + 1] == 2)
			{
				vecCopy(x4, x3);
				vecCopy(z4, z3);

#ifdef DEBUG
				printf("no swap\n");
				printf("[2]: pt3 holds [%d]P = ", x - d);
				print_vechex(work->pt3.X->data, 0, NWORDS, "");
				print_vechex(work->pt3.Z->data, 0, NWORDS, "");
				printf("[2]: pt1 holds [%d]P = ", d);
				print_vechex(work->pt1.X->data, 0, NWORDS, "");
				print_vechex(work->pt1.Z->data, 0, NWORDS, "");
#endif
				
				
			}
			else
			{
				// change point1
				vecCopy(x1, x3);
				vecCopy(z1, z3);
				vecCopy(x4, x1);
				vecCopy(z4, z1);

#ifdef DEBUG
				printf("next is a swap\n");
				printf("[2]: pt3 holds [%d]P = ", d);
				print_vechex(work->pt3.X->data, 0, NWORDS, "");
				print_vechex(work->pt3.Z->data, 0, NWORDS, "");
				printf("[2]: pt1 holds [%d]P = ", x - d);
				print_vechex(work->pt1.X->data, 0, NWORDS, "");
				print_vechex(work->pt1.Z->data, 0, NWORDS, "");
#endif

				d = x - d;

				if (seq[i + 1] == 3)
					i++;
			}


		}

		i++;

#ifdef DEBUG
		printf("==================\n");
#endif
	}

	if (x != c)
	{
		printf("expected %lu, euclid returned %lu\n", c, x);
		exit(1);
	}

	// copy out answer
	vecCopy(x2, P->X);
	vecCopy(z2, P->Z);

#ifdef DEBUG
	//exit(1);
#endif

	free(seq);
	return;
}

void next_pt_vec(vec_monty_t *mdata, ecm_work *work, ecm_pt *P, uint64_t c)
{
	uint64_t mask;
	vec_bignum_t *x1, *z1, *x2, *z2, *s1, *s2, *d1, *d2;
	uint64_t e, d;

	x1 = work->pt1.X;
	z1 = work->pt1.Z;
	x2 = work->pt2.X;
	z2 = work->pt2.Z;
	s1 = work->sum1;
	s2 = work->sum2;
	d1 = work->diff1;
	d2 = work->diff2;

	if (c == 1)
	{
		return;
	}

	vecCopy(P->X, x1);
	vecCopy(P->Z, z1);
	vecsubmod_ptr(P->X, P->Z, d1, work->n);
	vecaddmod_ptr(P->X, P->Z, s1, work->n);
	vec_duplicate(mdata, work, s1, d1, &work->pt2);

	//mulcnt[tid] += 3;
	//sqrcnt[tid] += 2;

	if (c == 2)
	{
		vecCopy(x2, P->X);
		vecCopy(z2, P->Z);
		return;
	}

	// find the first '1' bit then skip it
#ifdef __INTEL_COMPILER
	mask = 1ULL << (64 - _lzcnt_u64((uint64_t)c) - 2);
#elif defined(__GNUC__)
	mask = 1ULL << (64 - __builtin_clzll((uint64_t)c) - 2);
#elif defined _MSC_VER
    mask = 1ULL << (64 - _lzcnt_u64((uint64_t)c) - 2);
#endif

	//goal is to compute x_c, z_c using montgomery's addition
	//and multiplication formula's for x and z
	//the procedure will be similar to p+1 lucas chain formula's
	//but this time we are simultaneously incrementing both x and z
	//rather than just Vn.  In each bit of the binary expansion of
	//c, we need just one addition and one duplication for both x and z.
	//compute loop
	//for each bit of M to the right of the most significant bit
	d = 1;
	e = 2;
	while (mask > 0)
	{
		vecaddsubmod_ptr(x2, z2, s2, d2, work->n);
		vecaddsubmod_ptr(x1, z1, s1, d1, work->n);

		//if the bit is 1
		if (c & mask)
		{
			//add x1,z1, duplicate x2,z2
			vec_add(mdata, work, P, &work->pt1);
			vec_duplicate(mdata, work, s2, d2, &work->pt2);
			d = d + e;
			e *= 2;
		}
		else
		{
			//add x2,z2, duplicate x1,z1
			vec_add(mdata, work, P, &work->pt2);
			vec_duplicate(mdata, work, s1, d1, &work->pt1);
			e = e + d;
			d *= 2;
		}

		//mulcnt[tid] += 7;
		//sqrcnt[tid] += 4;

		mask >>= 1;
	}
	if (d != c)
	{
		printf("expected %lu, ladder returned %lu\n", c, d);
		exit(1);
	}
	vecCopy(x1, P->X);
	vecCopy(z1, P->Z);

	return;
}

void array_mul(uint64_t *primes, int num_p, mpz_t piprimes)
{
	mpz_t *p;
	int i;
	int alloc;
	uint64_t stg1 = (uint64_t)STAGE1_MAX;

	if (num_p & 1)
	{
		p = (mpz_t *)malloc((num_p + 1) * sizeof(mpz_t));
		alloc = num_p + 1;
	}
	else
	{
		p = (mpz_t *)malloc((num_p + 0) * sizeof(mpz_t));
		alloc = num_p;
	}

	for (i = 0; i < num_p; i++)
	{
		uint64_t c = 1;
		uint64_t q;

		q = primes[i];
		do {
			c *= q;
		} while ((c * q) < stg1);

		mpz_init(p[i]);
		mpz_set_ui(p[i], q);
	}

	if (num_p & 1)
	{
		mpz_init(p[i]);
		mpz_set_ui(p[i], 1);
		num_p++;
	}

	while (num_p != 1)
	{
		for (i = 0; i < num_p / 2; i++)
		{
			mpz_mul(p[i], p[i], p[i + num_p / 2]);
		}

		num_p /= 2;

		if ((num_p > 1) && (num_p & 1))
		{
			mpz_set_ui(p[num_p], 1);
			num_p++;
		}
	}

	mpz_set(piprimes, p[0]);

	printf("piprimes has %lu bits\n\n", mpz_sizeinbase(piprimes, 2));

	for (i = 0; i < alloc; i++)
	{
		mpz_clear(p[i]);
	}
	free(p);

	return;
}

void vececm(thread_data_t *tdata)
{
	//attempt to factor n with the elliptic curve method
	//following brent and montgomery's papers, and CP's book
    tpool_t *tpool_data;
    uint32_t threads = tdata[0].total_threads;
    base_t retval;
	base_t i, j;
	int curve;
	int tid;
	FILE *save;
	char fname[80];
	char *wstr;
	int found = 0;
    int result;
    uint64_t num_found;
	vec_bignum_t *one = vecInit();
    mpz_t gmpt, gmpn;
    int verbose = tdata[0].verbose;
    mpz_t tmp_factor;
    uint32_t curves_run = 0;

	// timing variables
	struct timeval stopt;	// stop time of this job
	struct timeval startt;	// start time of this job
	struct timeval fullstopt;	// stop time of this job
	struct timeval fullstartt;	// start time of this job
	double t_time;

	gettimeofday(&fullstartt, NULL);
    mpz_init(gmpt);
    mpz_init(gmpn);
    mpz_init(tmp_factor);

    // in this function, gmpn is only used to for screen output and to 
    // check factors.  So if this is a Mersenne input, use the original
    // input number.  (all math is still done relative to the base Mersenne.)
    if (tdata[0].mdata->isMersenne)
    {
        extract_bignum_from_vec_to_mpz(gmpn, tdata[0].mdata->vnhat, 0, NWORDS);
    }
    else
    {
        extract_bignum_from_vec_to_mpz(gmpn, tdata[0].mdata->n, 0, NWORDS);
    }

	for (i = 0; i < VECLEN; i++)
	{
		one->data[i] = 1;
	}
    one->size = 1;
    
    tpool_data = tpool_setup(tdata[0].total_threads, NULL, NULL, &vec_ecm_sync,
        &vec_ecm_dispatch, tdata);

    for (i = 0; i < threads; i++)
        tpool_data[i].debug = 0;

    tpool_add_work_fcn(tpool_data, &vec_ecm_build_curve_work_fcn);
    tpool_add_work_fcn(tpool_data, &vec_ecm_stage1_work_fcn);
    tpool_add_work_fcn(tpool_data, &vec_ecm_stage2_init_work_fcn);
    tpool_add_work_fcn(tpool_data, &vec_ecm_stage2_work_fcn);

	for (curve = 0; curve < tdata[0].curves; curve += VECLEN)
	{
        uint64_t p;

        if (verbose > 1)
            printf("\n");

        if (verbose >= 0)
        {
            printf("ecm: %d/%d curves on C%d @ B1=%lu, B2=100*B1\r",
                (curve + VECLEN) * threads, tdata[0].curves * threads,
                (int)gmp_base10(gmpn), STAGE1_MAX);
        }
        else if (verbose > 1)
        {
            printf("ecm: %d/%d curves on C%d @ B1=%lu, B2=100*B1\n",
                (curve + VECLEN) * threads, tdata[0].curves * threads,
                (int)gmp_base10(gmpn), STAGE1_MAX);
        }

		gettimeofday(&startt, NULL);
        
        // parallel curve building        
        for (i = 0; i < threads; i++)
        {
			tdata[i].work->stg1Add = 0;
			tdata[i].work->stg1Doub = 0;
			tdata[i].work->last_pid = 0;
            tdata[i].phase_done = 0;
            tdata[i].ecm_phase = 0;
        }
        tpool_go(tpool_data);

		gettimeofday(&stopt, NULL);
		t_time = yafu_difftime (&startt, &stopt);
		
        if (verbose > 1)
		    printf("Building curves took %1.4f seconds.\n",t_time);

		// parallel stage 1
		gettimeofday(&startt, NULL);

        for (p = 0; p < STAGE1_MAX; p += PRIME_RANGE)
        {
			// first condition will fetch more primes if B1 is large
			// and thus stage 1 takes several iterations of PRIME_RANGE.
			// second condition resets primes array, if it is modified
			// by stage 2, before starting a new batch of curves.
			if ((tdata[0].work->last_pid == NUM_P) || ((p == 0) && (P_MIN != 2)))
			{
				if (PRIMES != NULL) { free(PRIMES); PRIMES = NULL; };
				PRIMES = GetPRIMESRange(spSOEprimes, szSOEp, NULL,
					p, MIN(STAGE2_MAX + 1000, p + (uint64_t)PRIME_RANGE), &num_found);
				NUM_P = num_found;
				P_MIN = PRIMES[0];
				P_MAX = PRIMES[NUM_P - 1];

                if (verbose > 1)
				    printf("found %lu primes in range [%lu : %lu]\n", NUM_P, P_MIN, P_MAX);
			}

            for (i = 0; i < threads; i++)
            {
                tdata[i].phase_done = 0;
                tdata[i].ecm_phase = 1;
            }

            if (verbose > 1)
            {
                printf("commencing Stage 1 @ prime %lu\n", P_MIN);
            }

            // start the stage 1 thread pool
            tpool_go(tpool_data);

            if (p < STAGE1_MAX)
            {
                // record results for this batch of primes
                if (tdata[0].save_b1) // && (save_intermediate)
                {
                    sprintf(fname, "save_b1_intermediate.txt");
                    save = fopen(fname, "a");
                }

                gettimeofday(&stopt, NULL);
                t_time = yafu_difftime(&startt, &stopt);
                if (verbose > 1)
                    printf("Stage 1 current elapsed time is %1.4f seconds\n", t_time);

                for (j = 0; j < threads; j++)
                {

                    // GMP-ECM wants X/Z.
                    // or, equivalently, X and Z listed separately.
                    vecmulmod_ptr(tdata[j].P->X, one, tdata[j].work->tt4, tdata[j].work->n,
                        tdata[j].work->tt2, tdata[j].mdata);

                    vecmulmod_ptr(tdata[j].P->Z, one, tdata[j].work->tt3, tdata[j].work->n,
                        tdata[j].work->tt2, tdata[j].mdata);

                    //print_vechexbignum(tdata[j].P->Z, "Z point\n");

                    for (i = 0; i < VECLEN; i++)
                    {
                        extract_bignum_from_vec_to_mpz(gmpt, tdata[j].P->Z, i, NWORDS);

                        if (mpz_cmp_ui(gmpt, 0) == 0)
                        {
                            printf("something failed: tid = %d, vec = %d has zero result\n", j, i);
                        }

                        result = vec_check_factor(gmpt, gmpn, tmp_factor);
                        if (result == 1)
                        {

                            //FILE *out = fopen("ecm_results.txt", "a");

                            if (verbose > 1)
                            {
                                gmp_printf("\nfound factor %Zd in stage 1 in thread %d, vec position %d, with sigma = ",
                                    tmp_factor, j, i);
                                printf("%"PRIu64"\n", tdata[j].sigma[i]);
                            }

                            //if (out != NULL)
                            //{
                            //	gmp_fprintf(out, "\nfound factor %Zd in stage 1 at curve %d, "
                            //		"in thread %d, vec position %d, with sigma = ",
                            //        tmp_factor, threads * curve + j * VECLEN + i, j, i);
                            //    fprintf(out, "%"PRIu64"\n", tdata[j].sigma[i]);
                            //	fclose(out);
                            //}
                            //
                            //fflush(stdout);
                            found = 1;

                            if (tdata[j].numfactors == 0)
                            {
                                tdata[j].factors = (avx_ecm_factor_t*)malloc(1 * sizeof(avx_ecm_factor_t));
                            }
                            else
                            {
                                tdata[j].factors = (avx_ecm_factor_t*)realloc(tdata[j].factors,
                                    (tdata[j].numfactors + 1) * sizeof(avx_ecm_factor_t));
                            }

                            tdata[j].factors[tdata[j].numfactors].thread_id = j;
                            tdata[j].factors[tdata[j].numfactors].stg_id = 1;
                            tdata[j].factors[tdata[j].numfactors].sigma = tdata[j].sigma[i];
                            tdata[j].factors[tdata[j].numfactors].vec_id = i;
                            tdata[j].factors[tdata[j].numfactors].curve_id =
                                (curve * threads) + (j * VECLEN) + i;
                            mpz_init(tdata[j].factors[tdata[j].numfactors].factor);
                            mpz_set(tdata[j].factors[tdata[j].numfactors++].factor, tmp_factor);

                            //gmp_printf("thread %d now has %d factors in stg1: %Zd\n",
                            //    j, tdata[j].numfactors, tmp_factor);
                        }

                        if (tdata[0].save_b1)
                        {
                            fprintf(save, "METHOD=ECM; SIGMA=%"PRIu64"; B1=%"PRIu64"; ",
                                tdata[j].sigma[i], PRIMES[tdata[j].work->last_pid - 1]);
                            gmp_fprintf(save, "N=0x%Zx; ", gmpn);

                            extract_bignum_from_vec_to_mpz(gmpt, tdata[j].work->tt4, i, NWORDS);
                            gmp_fprintf(save, "X=0x%Zx; ", gmpt);

                            extract_bignum_from_vec_to_mpz(gmpt, tdata[j].work->tt3, i, NWORDS);
                            gmp_fprintf(save, "Z=0x%Zx; PROGRAM=AVX-ECM;\n", gmpt);
                        }
                    }

                }

                if (tdata[0].save_b1)
                {
                    fclose(save);
                }
            }
        }

        // record results for this batch of primes
        if (tdata[0].save_b1) // && (save_intermediate)
        {
            sprintf(fname, "save_b1.txt");
            save = fopen(fname, "a");
        }

        for (j = 0; j < threads; j++)
        {

            // GMP-ECM wants X/Z.
            // or, equivalently, X and Z listed separately.
            vecmulmod_ptr(tdata[j].P->X, one, tdata[j].work->tt4, tdata[j].work->n,
                tdata[j].work->tt2, tdata[j].mdata);

            vecmulmod_ptr(tdata[j].P->Z, one, tdata[j].work->tt3, tdata[j].work->n,
                tdata[j].work->tt2, tdata[j].mdata);

            //print_vechexbignum(tdata[j].P->Z, "Z point\n");

            for (i = 0; i < VECLEN; i++)
            {
                extract_bignum_from_vec_to_mpz(gmpt, tdata[j].P->Z, i, NWORDS);

                if (mpz_cmp_ui(gmpt, 0) == 0)
                {
                    printf("something failed: tid = %d, vec = %d has zero result\n", j, i);
                }

                result = vec_check_factor(gmpt, gmpn, tmp_factor);
                if (result == 1)
                {

                    //FILE *out = fopen("ecm_results.txt", "a");

                    if (verbose > 1)
                    {
                        gmp_printf("\nfound factor %Zd in stage 1 in thread %d, vec position %d, with sigma = ",
                            tmp_factor, j, i);
                        printf("%"PRIu64"\n", tdata[j].sigma[i]);
                    }

                    //if (out != NULL)
                    //{
                    //	gmp_fprintf(out, "\nfound factor %Zd in stage 1 at curve %d, "
                    //		"in thread %d, vec position %d, with sigma = ",
                    //        tmp_factor, threads * curve + j * VECLEN + i, j, i);
                    //    fprintf(out, "%"PRIu64"\n", tdata[j].sigma[i]);
                    //	fclose(out);
                    //}
                    //
                    //fflush(stdout);
                    found = 1;

                    if (tdata[j].numfactors == 0)
                    {
                        tdata[j].factors = (avx_ecm_factor_t*)malloc(1 * sizeof(avx_ecm_factor_t));
                    }
                    else
                    {
                        tdata[j].factors = (avx_ecm_factor_t*)realloc(tdata[j].factors,
                            (tdata[j].numfactors + 1) * sizeof(avx_ecm_factor_t));
                    }

                    tdata[j].factors[tdata[j].numfactors].thread_id = j;
                    tdata[j].factors[tdata[j].numfactors].stg_id = 1;
                    tdata[j].factors[tdata[j].numfactors].sigma = tdata[j].sigma[i];
                    tdata[j].factors[tdata[j].numfactors].vec_id = i;
                    tdata[j].factors[tdata[j].numfactors].curve_id =
                        (curve * threads) + (j * VECLEN) + i;
                    mpz_init(tdata[j].factors[tdata[j].numfactors].factor);
                    mpz_set(tdata[j].factors[tdata[j].numfactors++].factor, tmp_factor);

                    //gmp_printf("thread %d now has %d factors in stg1: %Zd\n",
                    //    j, tdata[j].numfactors, tmp_factor);
                }

                if (tdata[0].save_b1)
                {
                    fprintf(save, "METHOD=ECM; SIGMA=%"PRIu64"; B1=%"PRIu64"; ",
                        tdata[j].sigma[i], PRIMES[tdata[j].work->last_pid - 1]);
                    gmp_fprintf(save, "N=0x%Zx; ", gmpn);

                    extract_bignum_from_vec_to_mpz(gmpt, tdata[j].work->tt4, i, NWORDS);
                    gmp_fprintf(save, "X=0x%Zx; ", gmpt);

                    extract_bignum_from_vec_to_mpz(gmpt, tdata[j].work->tt3, i, NWORDS);
                    gmp_fprintf(save, "Z=0x%Zx; PROGRAM=AVX-ECM;\n", gmpt);
                }
            }

        }

        if (tdata[0].save_b1)
        {
            fclose(save);
        }

        gettimeofday(&stopt, NULL);
        t_time = yafu_difftime(&startt, &stopt);
        if (verbose > 1)
            printf("Stage 1 took %1.4f seconds\n", t_time);
	
        curves_run += VECLEN * threads;

        // always stop when a factor is found
        //if (found)
        //	break;


        if (DO_STAGE2)
        {
            uint64_t last_p = PRIMES[tdata[0].work->last_pid];

            // parallel stage 2
            gettimeofday(&startt, NULL);

			// stage 2 parallel init
            for (i = 0; i < threads; i++)
            {
                tdata[i].phase_done = 0;
                tdata[i].ecm_phase = 2;
            }
            tpool_go(tpool_data);

            for (; last_p < STAGE2_MAX; )
            {
				if (last_p == P_MAX)
				{
					if (PRIMES != NULL) { free(PRIMES); PRIMES = NULL; };
					PRIMES = GetPRIMESRange(spSOEprimes, szSOEp, NULL,
						last_p, MIN(last_p + (uint64_t)PRIME_RANGE, STAGE2_MAX + 1000), &num_found);
					NUM_P = num_found;
					P_MIN = PRIMES[0];
					P_MAX = PRIMES[NUM_P - 1];

					for (i = 0; i < threads; i++)
					{
						tdata[i].work->last_pid = 1;
					}

                    if (verbose > 1)
					    printf("found %lu primes in range [%lu : %lu]\n", NUM_P, P_MIN, P_MAX);
				}

                for (i = 0; i < threads; i++)
                {
                    tdata[i].phase_done = 0;
                    tdata[i].ecm_phase = 3;
                }
                tpool_go(tpool_data);

                //printf("stage 2 finished at P=%lu (maxP = %lu) (id %u of %lu)\n",
                //    PRIMES[tdata[0].work->last_pid], P_MAX, tdata[0].work->last_pid, NUM_P);

				if (tdata[0].work->last_pid == NUM_P)
				{
					// we ended at the last prime we cached.  Check if we
					// need to do more.  
					if (STAGE2_MAX == PRIME_RANGE)
						break;
					else
						last_p = P_MAX;
				}
				else
				{
					last_p = PRIMES[tdata[0].work->last_pid];
				}
            }

            gettimeofday(&stopt, NULL);
            t_time = yafu_difftime(&startt, &stopt);
            if (verbose > 1)
            {
                printf("Stage 2 took %1.4f seconds\n", t_time);
                printf("performed %d pair-multiplies for %lu primes in stage 2\n",
                    tdata[0].work->paired, tdata[0].work->numprimes);
                printf("performed %u point-additions and %u point-doubles in stage 2\n",
                    tdata[0].work->ptadds + tdata[0].work->stg1Add, tdata[0].work->stg1Doub);
            }

            for (j = 0; j < threads; j++)
            {
				for (i = 0; i < VECLEN; i++)
                {
                    extract_bignum_from_vec_to_mpz(gmpt, tdata[j].work->stg2acc, i, NWORDS);
                    result = vec_check_factor(gmpt, gmpn, tmp_factor);

                    if (mpz_cmp_ui(gmpt, 0) == 0)
                    {
                        printf("something failed: tid = %d, vec = %d has zero result\n", j, i);
                    }

                    if (result == 1)
                    {
						//FILE *out = fopen("ecm_results.txt", "a");

                        if (verbose > 1)
                        {
                            gmp_printf("\nfound factor %Zd in stage 2 in thread %d, vec position %d, with sigma = ",
                                tmp_factor, j, i);
                            printf("%"PRIu64"\n", tdata[j].sigma[i]);
                        }

						//if (out != NULL)
						//{
						//	gmp_fprintf(out, "\nfound factor %Zd in stage 1 at curve %d, "
						//		"in thread %d, vec position %d, with sigma = ",
                        //        tmp_factor, threads * curve + j * VECLEN + i, j, i);
                        //    fprintf(out, "%"PRIu64"\n", tdata[j].sigma[i]);
						//	fclose(out);
						//}
                        //
                        //fflush(stdout);
                        found = 1;

                        if (tdata[j].numfactors == 0)
                        {
                            tdata[j].factors = (avx_ecm_factor_t*)malloc(1 * sizeof(avx_ecm_factor_t));
                        }
                        else
                        {
                            tdata[j].factors = (avx_ecm_factor_t*)realloc(tdata[j].factors,
                                (tdata[j].numfactors + 1) * sizeof(avx_ecm_factor_t));
                        }

                        tdata[j].factors[tdata[j].numfactors].thread_id = j;
                        tdata[j].factors[tdata[j].numfactors].stg_id = 2;
                        tdata[j].factors[tdata[j].numfactors].sigma = tdata[j].sigma[i];
                        tdata[j].factors[tdata[j].numfactors].vec_id = i;
                        tdata[j].factors[tdata[j].numfactors].curve_id =
                            (curve * threads) + (j * VECLEN) + i;
                        mpz_init(tdata[j].factors[tdata[j].numfactors].factor);
                        mpz_set(tdata[j].factors[tdata[j].numfactors++].factor, tmp_factor);

                        //gmp_printf("thread %d now has %d factors in stg2: %Zd\n", 
                        //    j, tdata[j].numfactors, tmp_factor);
                    }
                }
            }
        }

        if (verbose >= 0)
        {
            printf("ecm: %d/%d curves on C%d @ B1=%lu, B2=100*B1\r",
                (curve + VECLEN) * threads, tdata[0].curves * threads,
                (int)gmp_base10(gmpn), STAGE1_MAX);
        }

		if (found)
			break;

	}

    if (verbose > 1)
    {
        printf("ecm: %d/%d curves on C%d @ B1=%lu, B2=100*B1\n",
            curves_run, tdata[0].curves * threads,
            (int)gmp_base10(gmpn), STAGE1_MAX);
    }
    else
    {
        printf("\n");
    }

    tdata[0].curves = curves_run;

	gettimeofday(&fullstopt, NULL);
	t_time = yafu_difftime(&fullstartt, &fullstopt);
    if (verbose > 1)
	    printf("Process took %1.4f seconds.\n", t_time);

	vecClear(one);
    mpz_clear(gmpn);
    mpz_clear(gmpt);
    mpz_clear(tmp_factor);
	return;
}

//#define PRINT_DEBUG

void vec_build_one_curve(thread_data_t *tdata, mpz_t X, mpz_t Z, mpz_t A, uint64_t sigma)
{
    vec_monty_t *mdata = tdata->mdata;
    ecm_work *work = tdata->work;
    uint32_t tid = tdata->tid;
    ecm_pt *P = &work->pt1;

    mpz_t n, u, v, t1, t2, t3, t4;
    mpz_init(n);
    mpz_init(u);
    mpz_init(v);
    mpz_init(t1);
    mpz_init(t2);
    mpz_init(t3);
    mpz_init(t4);

    extract_bignum_from_vec_to_mpz(n, mdata->n, 0, NWORDS);

    if (sigma == 0)
    {
        do
        {
            work->sigma = spRandp(&tdata->lcg_state, 6, VEC_MAXDIGIT);  //lcg_rand(&tdata->lcg_state);
        } while (work->sigma < 6);
    }
    else
    {
        work->sigma = sigma;
    }

    //sigma = work->sigma = 1632562926;
    //sigma = 269820583;
    //sigma = 50873471;
    //sigma = 444711979;		// both good
    //work->sigma = sigma;
    //printf("thread %d running curve on sigma = %"PRIu64"\n", tid, work->sigma);

#ifdef PRINT_DEBUG
    printf("sigma = %lu\n", work->sigma);
#endif

    // v = 4*sigma
    mpz_set_ui(v, work->sigma);
    mpz_mul_2exp(v, v, 2);

#ifdef PRINT_DEBUG
    gmp_printf("v = %Zx\n", v);
#endif

    // u = sigma^2 - 5
    mpz_set_ui(u, work->sigma);
    mpz_mul(u, u, u);
    mpz_sub_ui(u, u, 5);

#ifdef PRINT_DEBUG
    gmp_printf("u = sigma^2 - 5 = %Zx\n", u);
#endif

    // x = u^3
    mpz_mul(X, u, u);
    mpz_mul(X, X, u);
    mpz_tdiv_r(X, X, n);

#ifdef PRINT_DEBUG
    gmp_printf("u^3 = %Zx\n", X);
#endif

    // z = v^3
    mpz_mul(Z, v, v);
    mpz_mul(Z, Z, v);
    mpz_tdiv_r(Z, Z, n);

#ifdef PRINT_DEBUG
    gmp_printf("v^3 = %Zx\n", Z);
#endif


    // compute parameter A
    // (v-u)
    if (mpz_cmp(u, v) > 0)
    {
        mpz_sub(t1, v, u);
        mpz_add(t1, t1, n);
    }
    else
    {
        mpz_sub(t1, v, u);
    }

    // (v-u)^3
    mpz_mul(t2, t1, t1);
    mpz_tdiv_r(t2, t2, n);
    mpz_mul(t4, t2, t1);
    mpz_tdiv_r(t4, t4, n);

    // 3u + v
    mpz_mul_ui(t1, u, 3);
    mpz_add(t3, t1, v);
    mpz_tdiv_r(t3, t3, n);

    // a = (v-u)^3 * (3u + v)
    mpz_mul(t1, t3, t4);
    mpz_tdiv_r(t1, t1, n);

#ifdef PRINT_DEBUG
    gmp_printf("(v-u)^3 = %Zx\n", t4);
    gmp_printf("(3u + v) = %Zx\n", t3);
    gmp_printf("a = %Zx\n", t1);
#endif

    if (0)
    {
        // This is how gmp-ecm does things since sometimes they want
        // the Montgomery parameter A.  We always use Suyama's parameterization
        // so we just go ahead and build (A+2)/4.
        // 4*x*v
        mpz_mul_ui(t2, X, 4);
        mpz_mul(t4, t2, v);
        mpz_tdiv_r(t4, t4, n);

        // 4*x*v * z
        mpz_mul(t3, t4, Z);
        mpz_tdiv_r(t3, t3, n);

        /* u = 1 / (v^3 * 4*u^3*v) */
        mpz_invert(t2, t3, n);

        /* v = z^(-1) (mod n)  = 1 / v^3   */
        mpz_mul(v, t2, t4);   
        mpz_tdiv_r(v, v, n);

        /* x = x * z^(-1)      = u^3 / v^3 */
        mpz_mul(X, X, v);
        mpz_tdiv_r(X, X, n);

        /* v = b^(-1) (mod n)  = 1 / 4*u^3*v */
        mpz_mul(v, t2, Z);
        mpz_tdiv_r(v, v, n);

        /* t = ((v-u)^3 * (3*u+v)) / 4*u^3*v */
        mpz_mul(t1, t1, v);       
        mpz_tdiv_r(t1, t1, n);

        /* A = ((v-u)^3 * (3*u+v)) / 4*u^3*v - 2*/
        mpz_sub_ui(A, t1, 2);
        if (mpz_sgn(A) < 0)
        {
            mpz_add(A, A, n);
        }

        mpz_mul_2exp(X, X, DIGITBITS * NWORDS);
        mpz_tdiv_r(X, X, n);
        mpz_mul_2exp(Z, Z, DIGITBITS * NWORDS);
        mpz_tdiv_r(Z, Z, n);
        mpz_mul_2exp(A, A, DIGITBITS * NWORDS);
        mpz_tdiv_r(A, A, n);

        //gmp_printf("X/Z = %Zx\n", X);
        //gmp_printf("Z = %Zx\n", Z);
        //gmp_printf("A = %Zx\n", A);
        //
        //printf("(A+2)*B/4\n");

        mpz_set_ui(t1, 2);
        mpz_mul_2exp(t1, t1, DIGITBITS * NWORDS);
        mpz_add(A, A, t1);
        mpz_tdiv_r(A, A, n);

        if (mpz_odd_p(A))
        {
            mpz_add(A, A, n);
        }
        mpz_tdiv_q_2exp(A, A, 1);

        if (mpz_odd_p(A))
        {
            mpz_add(A, A, n);
        }
        mpz_tdiv_q_2exp(A, A, 1);

        //gmp_printf("A = %Zx\n", A);
        //exit(1);
    }
    else
    {
        // We always use Suyama's parameterization
        // so we just go ahead and build (A+2)/4.

        // 16*u^3*v
        mpz_mul_ui(t2, X, 16);
        mpz_mul(t4, t2, v);
        mpz_tdiv_r(t4, t4, n);

#ifdef PRINT_DEBUG
        gmp_printf("16*u^3*v = %Zx\n", t4);
#endif

        // accomplish the division by multiplying by the modular inverse
        // of the denom
        mpz_invert(t2, t4, n);

#ifdef PRINT_DEBUG
        gmp_printf("inv = %Zx\n", t2);
#endif

        // t1 = b = (v - u)^3 * (3*u + v) / 16u^3v)
        mpz_mul(A, t1, t2);
        mpz_tdiv_r(A, A, n);

#ifdef PRINT_DEBUG
        gmp_printf("b = %Zx\n", A);
#endif

        mpz_invert(t1, Z, n);
        mpz_mul(X, X, t1);
        mpz_set_ui(Z, 1);

        if (mdata->isMersenne)
        {
            // into Monty rep
            mpz_mul_2exp(X, X, DIGITBITS * NWORDS);
            mpz_tdiv_r(X, X, n);
            mpz_mul_2exp(Z, Z, DIGITBITS * NWORDS);
            mpz_tdiv_r(Z, Z, n);
            mpz_mul_2exp(A, A, DIGITBITS * NWORDS);
            mpz_tdiv_r(A, A, n);
        }
        else
        {
            mpz_tdiv_r(X, X, n);
            mpz_tdiv_r(Z, Z, n);
            mpz_tdiv_r(A, A, n);
        }
        

        //gmp_printf("X = %Zx\n", X);
        //gmp_printf("Z = %Zx\n", Z);
        //gmp_printf("A = %Zx\n", A);
        //exit(1);
    }
	
    mpz_clear(n);
    mpz_clear(u);
    mpz_clear(v);
    mpz_clear(t1);
    mpz_clear(t2);
    mpz_clear(t3);
    mpz_clear(t4);

	return;
}

//#define PRINT_DEBUG

void build_one_curve_param1(thread_data_t *tdata, mpz_t X, mpz_t Z, 
    mpz_t X2, mpz_t Z2, mpz_t A, uint64_t sigma)
{
    vec_monty_t *mdata = tdata->mdata;
    ecm_work *work = tdata->work;
    uint32_t tid = tdata->tid;
    ecm_pt *P = &work->pt1;
    int i;
    mpz_t n, v, d, dorig;
    mpz_init(n);
    mpz_init(v);
    mpz_init(d);
    mpz_init(dorig);

    extract_bignum_from_vec_to_mpz(n, mdata->n, 0, NWORDS);

    if (sigma == 0)
    {
        do
        {
            work->sigma = spRandp(&tdata->lcg_state, 6, VEC_MAXDIGIT) & 0X3FFFFFF;  //lcg_rand(&tdata->lcg_state) & 0x3ffffff;
        } while (work->sigma < 6);
    }
    else
    {
        work->sigma = sigma;
    }
    
#ifdef PRINT_DEBUG
    printf("sigma = %lu\n", work->sigma);
#endif

    // v = sigma^2
    mpz_set_ui(v, work->sigma);
    mpz_mul(v, v, v);

    /* A=4*d-2 with d = sigma^2/2^GMP_NUMB_BITS*/
      /* Compute d = sigma^2/2^GMP_NUMB_BITS */
    for (i = 0; i < 52; i++)
    {
        if (mpz_tstbit(v, 0) == 1)
            mpz_add(v, v, n);
        mpz_div_2exp(v, v, 1);
    }

    mpz_mod(v, v, n);

#ifdef PRINT_DEBUG
    gmp_printf("n = %Zd\n", n);
#endif
#ifdef PRINT_DEBUG
    gmp_printf("d = %Zd\n", v);
#endif
    mpz_set(dorig, v);

    /* TODO add d!=-1/8*/
    if (mpz_sgn(v) == 0 || mpz_cmp_ui(v, 1) == 0)
    {
        printf("parameter cannot be 0 or 1\n");
        exit(1);
    }

    mpz_mul_2exp(v, v, 2);           /* 4d */
    mpz_sub_ui(v, v, 2);             /* 4d-2 */

    // copy and convert to Montgomery representation
#ifdef PRINT_DEBUG
    gmp_printf("4d-2 = %Zd\n", v);
#endif
    //mpres_set_z(A, v, n);
    mpz_mul_2exp(A, v, MAXBITS);
    mpz_mod(A, A, n);
#ifdef PRINT_DEBUG
    gmp_printf("A = %Zd\n", A);
#endif

    /* Compute d=(A+2)/4 from A and d'=B*d thus d' = 2^(GMP_NUMB_BITS-2)*(A+2) */
    mpres_get_z(d, A, n, mdata);
#ifdef PRINT_DEBUG
    gmp_printf("d (mpres_get_z) = %Zd\n", d);
#endif
    mpz_add_ui(d, d, 2);
#ifdef PRINT_DEBUG
    gmp_printf("d (mpz_add_ui) = %Zd\n", d);
#endif
    mpz_mul_2exp(d, d, 52 - 2);
#ifdef PRINT_DEBUG
    gmp_printf("d (mpz_mul_2exp) = %Zd\n", d);
#endif
    mpz_mod(d, d, n);
#ifdef PRINT_DEBUG
    gmp_printf("d (mpz_mod) = %Zd\n", d);
#endif

    mpz_set(A, d);

#ifdef PRINT_DEBUG
    gmp_printf("d' = (A+2)*2^(GMP_NUMB_BITS-2) = %Zd\n", d);
#endif
#ifdef PRINT_DEBUG
    gmp_printf("dorig = %Zd\n", dorig);
#endif

    //mpres_set_ui(X, 2, n);
    mpz_set_ui(X, 2);
    mpz_mul_2exp(X, X, MAXBITS);
    mpz_mod(X, X, n);

#ifdef PRINT_DEBUG
    gmp_printf("X = %Zd\n", X);
#endif

    //mpres_set_ui(Z, 1, n);
    mpz_set_ui(Z, 1);
    mpz_mul_2exp(Z, Z, MAXBITS);
    mpz_tdiv_r(Z, Z, n);

#ifdef PRINT_DEBUG
    gmp_printf("Z = %Zd\n", Z);
#endif

    /* Compute 2P : no need to duplicate P, the coordinates are simple. */
    mpz_set_ui(X2, 9);
    mpz_mul_2exp(X2, X2, MAXBITS);
    mpz_tdiv_r(X2, X2, n);

#ifdef PRINT_DEBUG
    gmp_printf("X2 = %Zd\n", X2);
#endif

    mpz_set(Z2, d);
    mpz_mul_2exp(Z2, Z2, MAXBITS);
    mpz_tdiv_r(Z2, Z2, n);
    mpres_div_2exp(Z2, Z2, GMP_NUMB_BITS, n);

    mpres_mul_2exp(Z2, Z2, 6, n);
    mpres_add_ui(Z2, Z2, 8, n); /* P2 <- 2P = (9 : : 64d+8) */

#ifdef PRINT_DEBUG
    gmp_printf("Z2 = %Zd\n", Z2);
#endif

    mpz_clear(n);
    mpz_clear(v);
    mpz_clear(d);

    return;
}

//#define TESTMUL
void vec_ecm_stage1(vec_monty_t *mdata, ecm_work *work, ecm_pt *P, base_t b1, base_t *primes, int verbose)
{
	int i;
	uint64_t q;
	uint64_t stg1 = (uint64_t)STAGE1_MAX;

#ifdef PARAM1
    {
        mpz_t s, f;
        // timing variables
        struct timeval stopt;	// stop time of this job
        struct timeval startt;	// start time of this job
        double t_time;

        gettimeofday(&startt, NULL);
        mpz_init(s);
        mpz_init(f);
        work->last_pid = compute_s(s, PRIMES, stg1);
        gettimeofday(&stopt, NULL);
        t_time = yafu_difftime(&startt, &stopt);
        printf("Built product of %lu bits in %1.0f ms\n", mpz_sizeinbase(s, 2), t_time * 1000);
        ecm_stage1_batch(f, work, P, work->s, work->n, stg1, s, mdata);
        mpz_clear(s);
        mpz_clear(f);
        return;
    }
#endif

	// handle the only even case 
	q = 2;
	while (q < STAGE1_MAX)
	{
		vecsubmod_ptr(P->X, P->Z, work->diff1, work->n);
		vecaddmod_ptr(P->X, P->Z, work->sum1, work->n);
		vec_duplicate(mdata, work, work->sum1, work->diff1, P);
		q *= 2;
	}

	for (i = 1; (i < NUM_P) && ((uint32_t)PRIMES[i] < STAGE1_MAX); i++)
	{
		uint64_t c = 1;
	
		q = PRIMES[i];
		do {
			vec_prac(mdata, work, P, q);
			c *= q;
		} while ((c * q) < stg1);
	
#ifdef SKYLAKEX
		if ((verbose > 1) && ((i & 8191) == 0))
#else
        if ((verbose > 1) && ((i & 511) == 0))
#endif
		{
			printf("accumulating prime %lu\r", q);
			fflush(stdout);
		}
	}

	work->last_pid = i;

    if (verbose > 1)
	{
		printf("\nStage 1 completed at prime %lu with %u point-adds and %u point-doubles\n", 
			PRIMES[i-1], work->stg1Add, work->stg1Doub);
		fflush(stdout);
	}
	return;
}

void vec_ecm_stage2_init(ecm_pt *P, vec_monty_t *mdata, ecm_work *work, base_t *primes, int verbose)
{
	// run Montgomery's PAIR algorithm.  
	uint32_t D = work->D;
	uint32_t R = work->R;
	uint32_t w = D;
	uint32_t U = work->U;
	uint32_t L = work->L;
	uint32_t umax = U * w;
	int i, j, k, pid;
	
	uint32_t ainc = 2 * D;
	uint32_t ascale = 2 * D;
	uint32_t amin = work->amin = (STAGE1_MAX + w) / ascale;
	uint32_t s;
	uint32_t a;
	
	uint32_t numR;
	uint32_t u, ap;
	int q, mq;
	int debug = 0;

	uint32_t *rprime_map_U = work->map;
	ecm_pt *Pa = work->Pa;
	vec_bignum_t **Paprod = work->Paprod;
	vec_bignum_t **Pbprod = work->Pbprod;
	ecm_pt *Pb = work->Pb;
	ecm_pt *Pd;
	vec_bignum_t *acc = work->stg2acc;


    if (verbose > 1)
		printf("\n");

	work->paired = 0;
	work->numprimes = 0;
	work->ptadds = 0;
	work->stg1Add = 0;
	work->stg1Doub = 0;

	//stage 2 init
	//Q = P = result of stage 1
	//compute [d]Q for 0 < d <= D
	Pd = &Pb[rprime_map_U[w]];

	// [1]Q
	vecCopy(P->Z, Pb[1].Z);
	vecCopy(P->X, Pb[1].X);
	vecmulmod_ptr(Pb[1].X, Pb[1].Z, Pbprod[1], work->n, work->tt4, mdata);

	// [2]Q
	vecCopy(P->Z, Pb[2].Z);
	vecCopy(P->X, Pb[2].X);
    vecaddsubmod_ptr(P->X, P->Z, work->sum1, work->diff1, work->n);
	//vecsubmod_ptr(P->X, P->Z, work->diff1, work->n);
	//vecaddmod_ptr(P->X, P->Z, work->sum1, work->n);
	vec_duplicate(mdata, work, work->sum1, work->diff1, &Pb[2]);
	vecmulmod_ptr(Pb[2].X, Pb[2].Z, Pbprod[2], work->n, work->tt4, mdata);

	// [3]Q, [4]Q, ... [D]Q
	// because of the way we pick 'a' and 'D', we'll only need the 479 points that 
	// are relatively prime to D.  We need to keep a few points during the running 
	// computation i.e., ... j, j-1, j-2. The lookup table rprime_map maps the 
	// integer j into the relatively prime indices that we store. We also store
	// points 1, 2, and D, and storage point 0 is used as scratch for a total of
	// 483 stored points.

	vecCopy(Pb[1].X, work->pt2.X);
	vecCopy(Pb[1].Z, work->pt2.Z);
	vecCopy(Pb[2].X, work->pt1.X);
	vecCopy(Pb[2].Z, work->pt1.Z);

	for (j = 3; j <= U * D; j++)
	{
		ecm_pt *P1 = &work->pt1;			// Sd - 1
		ecm_pt *P2 = &Pb[1];				// S1
		ecm_pt *P3 = &work->pt2;			// Sd - 2
		ecm_pt *Pout = &Pb[rprime_map_U[j]];	// Sd

		// vecAdd:
		//x+ = z- * [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
		//z+ = x- * [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
		//x- = original x
		//z- = original z

		// compute Sd from Sd-1 + S1, requiring Sd-1 - S1 = Sd-2
        vecaddsubmod_ptr(P1->X, P1->Z, work->sum1, work->diff1, work->n);
        //vecsubmod_ptr(P1->X, P1->Z, work->diff1, work->n);
		//vecaddmod_ptr(P1->X, P1->Z, work->sum1, work->n);
        vecaddsubmod_ptr(P2->X, P2->Z, work->sum2, work->diff2, work->n);
		//vecsubmod_ptr(P2->X, P2->Z, work->diff2, work->n);
		//vecaddmod_ptr(P2->X, P2->Z, work->sum2, work->n);

		vecmulmod_ptr(work->diff1, work->sum2, work->tt1, work->n, work->tt4, mdata);	//U
		vecmulmod_ptr(work->sum1, work->diff2, work->tt2, work->n, work->tt4, mdata);	//V

        vecaddsubmod_ptr(work->tt1, work->tt2, Pout->X, Pout->Z, work->n);
		//vecaddmod_ptr(work->tt1, work->tt2, Pout->X, work->n);		//U + V
		//vecsubmod_ptr(work->tt1, work->tt2, Pout->Z, work->n);		//U - V
		vecsqrmod_ptr(Pout->X, work->tt1, work->n, work->tt4, mdata);					//(U + V)^2
		vecsqrmod_ptr(Pout->Z, work->tt2, work->n, work->tt4, mdata);					//(U - V)^2

		// if gcd(j,D) != 1, Pout maps to scratch space (Pb[0])
		vecmulmod_ptr(work->tt1, P3->Z, Pout->X, work->n, work->tt4, mdata);			//Z * (U + V)^2
		vecmulmod_ptr(work->tt2, P3->X, Pout->Z, work->n, work->tt4, mdata);			//x * (U - V)^2

		//store Pb[j].X * Pb[j].Z as well
		vecmulmod_ptr(Pout->X, Pout->Z, Pbprod[rprime_map_U[j]],
			work->n, work->tt4, mdata);

		work->ptadds++;

		// advance
		vecCopy(P1->X, P3->X);
		vecCopy(P1->Z, P3->Z);
		vecCopy(Pout->X, P1->X);
		vecCopy(Pout->Z, P1->Z);
	}

	//printf("B table generated to umax = %d\n", U * D);

	// Pd = [2w]Q
	vecCopy(P->Z, Pd->Z);
	vecCopy(P->X, Pd->X);
	next_pt_vec(mdata, work, Pd, ainc);
	//vec_prac(mdata, work, Pd, ainc);

	//first a value: first multiple of D greater than B1
	work->A = amin * ascale;

	//initialize info needed for giant step
	vecCopy(P->Z, Pa[0].Z);
	vecCopy(P->X, Pa[0].X);
	next_pt_vec(mdata, work, &Pa[0], work->A);
	//vec_prac(mdata, work, &Pa[0], work->A);

	//and Paprod
	vecmulmod_ptr(Pa[0].X, Pa[0].Z, Paprod[0], work->n, work->tt4, mdata);
	if (verbose & (debug == 2))
		printf("Pa[0] = [%lu]Q\n", work->A);

	vecCopy(P->Z, work->Pad->Z);
	vecCopy(P->X, work->Pad->X);
	next_pt_vec(mdata, work, work->Pad, work->A - ainc);
	//vec_prac(mdata, work, work->Pad, work->A - ainc);

	if (verbose & (debug == 2))
		printf("Pad = [%lu]Q\n", work->A - ainc);

	vecaddmod_ptr(Pa[0].X, Pa[0].Z, work->sum1, work->n);
	vecaddmod_ptr(Pd->X, Pd->Z, work->sum2, work->n);
	vecsubmod_ptr(Pa[0].X, Pa[0].Z, work->diff1, work->n);
	vecsubmod_ptr(Pd->X, Pd->Z, work->diff2, work->n);
	vec_add(mdata, work, work->Pad, &Pa[1]);
	vecmulmod_ptr(Pa[1].X, Pa[1].Z, Paprod[1], work->n, work->tt4, mdata);

	work->A += ainc;
	if (verbose & (debug == 2))
		printf("Pa[1] = [%lu]Q\n", work->A + ainc);

	for (i = 2; i < 2 * L; i++)
	{
		//giant step - use the addition formula for ECM
		//Pa + Pd
		//x+ = z- * [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
		//z+ = x- * [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
		//x- = [a-d]x
		//z- = [a-d]z
		//vecaddmod_ptr(Pa[i - 1].X, Pa[i - 1].Z, work->sum1, work->n);
		//vecaddmod_ptr(Pd->X, Pd->Z, work->sum2, work->n);
		//vecsubmod_ptr(Pa[i - 1].X, Pa[i - 1].Z, work->diff1, work->n);
		//vecsubmod_ptr(Pd->X, Pd->Z, work->diff2, work->n);
        vecaddsubmod_ptr(Pa[i - 1].X, Pa[i - 1].Z, work->sum1, work->diff1, work->n);
        vecaddsubmod_ptr(Pd->X, Pd->Z, work->sum2, work->diff2, work->n);
		vec_add(mdata, work, &Pa[i - 2], &Pa[i]);

		work->A += ainc;
		work->ptadds++;
		if (verbose & (debug == 2))
			printf("Pa[%d] = [%lu]Q\n", i, work->A + i * ainc);

		//and Paprod
		vecmulmod_ptr(Pa[i].X, Pa[i].Z, Paprod[i], work->n, work->tt4, mdata);
	}

	if (verbose & (debug == 2))
		printf("A table generated to 2 * L = %d\n", 2 * L);

	// initialize accumulator
    vecCopy(mdata->one, acc);

	return;
}

// batch inversions:
// From: H. Shacham and D. Boneh, "Improving SSL Handshake Performance via Batching"
// A general batched - inversion algorithm proceeds, in three phases, as follows.
// first, set A1 = x1 and Ai = xi * A(i-1) so that Ai = prod(j=1,i,xj).
// then, invert An and store in Bn and set Bi = x(i+1) * B(i+1) for i < n.
// Now we have Bi = prod(j=1,i,xj^-1).
// finally, set C1 = B1 and Ci = A(i-1) * B(i) for i > 1.
// Then Ci = xi^-1 for i > 1.
// each phase takes n-1 multiplications so we have 3n-3 total multiplications
// and one inversion mod N.

// Now, how do we use this to speed up stage 2?
// In projective coordinates (Px:Pz) = (Qx:Qz) does not
// imply that Px=Qx.  We need to cancel the Z-coordinates.
// to cancel a z-coord we can multiply a point by z^-1, then
// we have (Px/Pz:1).
// the cross product we compute and accumulate is:
// (Ax - Bx) * (Az + Bz) + BxBz - AxAz
// AxBz - BxAz
// In the baby-step giant-step algorithm for stage 2 we can normalize
// all pre-computed Pb's by the batch-inversion process.
// Likewise if we precompute all Pa's we can normalize in the same way.
// (If there are too many Pa's then maybe in batches at a time?)
// Then we simply work with the normalized values and directly 
// accumulate the cross products (Xr/Zr*1 - Xd/Zd*1)?
// does the final accumulation then need another inversion?



// _TEST
#define CROSS_PRODUCT_TEST \
vecsubmod_ptr(Pa[pa].X, Pb[rprime_map_U[pb]].X, work->tt1, work->n);          \
vecmulmod_ptr(acc, work->tt1, acc, work->n, work->tt4, mdata);        

// pre-computing the sum/diff multiplies is not efficient - many
// of the sum/diff product combinations will never be used
#define CROSS_PRODUCT \
vecsubmod_ptr(Pa[pa].X, Pb[rprime_map_U[pb]].X, work->tt1, work->n);          \
vecaddmod_ptr(Pa[pa].Z, Pb[rprime_map_U[pb]].Z, work->tt2, work->n);          \
vecmulmod_ptr(work->tt1, work->tt2, work->tt3, work->n, work->tt4, mdata);    \
vecaddmod_ptr(work->tt3, Pbprod[rprime_map_U[pb]], work->tt1, work->n);       \
vecsubmod_ptr(work->tt1, Paprod[pa], work->tt2, work->n);                     \
vecmulmod_ptr(acc, work->tt2, acc, work->n, work->tt4, mdata);                


void vec_ecm_stage2_pair(ecm_pt *P, vec_monty_t *mdata, ecm_work *work, base_t *primes, int verbose)
{
	// run Montgomery's PAIR algorithm.  
	uint32_t D = work->D;
	uint32_t R = work->R;
	uint32_t w = D;
	uint32_t U = work->U;
	uint32_t L = work->L;
	uint32_t umax = U * w;
	int i, j, k, pid;
	Queue_t **Q = work->Q;
	uint32_t ainc = 2 * D;
	uint32_t ascale = 2 * D;
	uint32_t amin = work->amin;
	uint64_t s;
	uint32_t a;
	uint32_t *map = work->Qmap;
	uint32_t *rmap = work->Qrmap;
	uint32_t numR = R - 3;
	uint32_t u, ap;
	int q, mq;
	int debug = 0;

	uint32_t *rprime_map_U = work->map;
	ecm_pt *Pa = work->Pa;
	vec_bignum_t **Paprod = work->Paprod;
	vec_bignum_t **Pbprod = work->Pbprod;
	ecm_pt *Pb = work->Pb;
	ecm_pt *Pd = &Pb[rprime_map_U[w]];
	vec_bignum_t *acc = work->stg2acc;
    uint32_t aranges = 0;

	if (verbose > 1)
		printf("\n");

	pid = work->last_pid;
    if (verbose > 1)
	{
		printf("commencing stage 2 at p=%lu, A=%u\n"
			"w = %u, R = %u, L = %u, umax = %u, amin = %u\n",
			PRIMES[pid], amin * ascale, w, numR, L, umax, amin);
	}

	while ((pid < NUM_P) && (PRIMES[pid] < STAGE2_MAX))
	{
		work->numprimes++;

		s = PRIMES[pid++];
		a = (s + w) / ascale;

		if ((verbose > 1) && ((pid & 32767) == 0))
		{
			printf("accumulating prime %lu\r", PRIMES[pid]);
			fflush(stdout);
		}

		// new range of a
		while (a >= (amin + L))
		{
			uint32_t oldmin = amin;

			amin = amin + L - U;
			
			if (verbose & (debug == 2))
				printf("dumping tables from %u to %u\n", oldmin, amin-1); fflush(stdout);

			for (i = 0; i < numR; i++)
			{
				while (Q[i]->len > 0)
				{
					if (peekqueue(Q[i]) < amin)
					{
						//printf("(%u,%u)\n", (2 * dequeue(Q[i])), rmap[numR - i - 1]);
						int pa = dequeue(Q[i]) - oldmin;
						int pb = rmap[numR - i - 1];

						if (verbose & debug) {
							if (((((oldmin + pa) * ascale) + pb) == 3266117) ||
								((((oldmin + pa) * ascale) - pb) == 3266117))
							{
								printf("5: accumulating (%u,%u) = 3266117\n", pa, pb); fflush(stdout);
								printf("a  = %u\n", a);
								printf("ap  = %u\n", ap);
								printf("q  = %d\n", q);
								printf("mq = %d\n", mq);
								printf("u = %d\n", u);
								printf("amin = %d\n", amin);
							}
						}

						// accumulate the cross product  (zimmerman syntax).
						// page 342 in C&P
                        CROSS_PRODUCT;

						work->paired++;
					}
					else
					{
						break;
					}
				}
			}

			for (i = numR; i < 2 * numR; i++)
			{
				while (Q[i]->len > 0)
				{
					if (peekqueue(Q[i]) < amin)
					{
						//printf("(%u,%u)\n", (2 * dequeue(Q[i])), rmap[i - numR]);
						int pa = dequeue(Q[i]) - oldmin;
						int pb = rmap[i - numR];

						if (verbose & debug) {
							if (((((oldmin + pa) * ascale) + pb) == 3266117) ||
								((((oldmin + pa) * ascale) - pb) == 3266117))
							{
								printf("6: accumulating (%u,%u) = 3266117\n", pa, pb); fflush(stdout);
								printf("a  = %u\n", a);
								printf("ap  = %u\n", ap);
								printf("q  = %d\n", q);
								printf("mq = %d\n", mq);
								printf("u = %d\n", u);
								printf("amin = %d\n", amin);
							}
						}

						//printf("accumulate (%u,%u)\n", pa, pb); fflush(stdout);
						// accumulate the cross product  (zimmerman syntax).
						// page 342 in C&P
                        CROSS_PRODUCT;
						work->paired++;
					}
					else
					{
						break;
					}
				}
			}

			// shift out uneeded A's
			j = 0;
			for (i = (amin - oldmin); i < (2 * L); i++, j++)
			{
				vecCopy(Pa[i].X, Pa[j].X);
				vecCopy(Pa[i].Z, Pa[j].Z);
				vecCopy(Paprod[i], Paprod[j]);

				if (verbose & (debug == 2))
					printf("Pa[%d] = [%d]Q\n", j, oldmin * ascale + i * ainc);
			}

			// make new A's using the last two points
			//printf("making new A's from %d to %d\n", (2 * L) - (amin - oldmin), (2 * L));
			for (i = (2 * L) - (amin - oldmin); i < (2 * L); i++)
			{
				//giant step - use the addition formula for ECM
				//Pa + Pd
				//x+ = z- * [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
				//z+ = x- * [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
				//x- = [a-d]x
				//z- = [a-d]z
				//vecaddmod_ptr(Pa[i - 1].X, Pa[i - 1].Z, work->sum1, work->n);
				//vecaddmod_ptr(Pd->X, Pd->Z, work->sum2, work->n);
				//vecsubmod_ptr(Pa[i - 1].X, Pa[i - 1].Z, work->diff1, work->n);
				//vecsubmod_ptr(Pd->X, Pd->Z, work->diff2, work->n);
                vecaddsubmod_ptr(Pa[i - 1].X, Pa[i - 1].Z, work->sum1, work->diff1, work->n);
                vecaddsubmod_ptr(Pd->X, Pd->Z, work->sum2, work->diff2, work->n);
				vec_add(mdata, work, &Pa[i - 2], &Pa[i]);

				//and Paprod
				vecmulmod_ptr(Pa[i].X, Pa[i].Z, Paprod[i], work->n, work->tt4, mdata);

				work->A += ainc;
				
				if (verbose & (debug == 2))
					printf("Pa[%d] = [%lu]Q\n", i, work->A);

				work->ptadds++;
			}

            aranges++;

			if (verbose & (debug == 2))
				printf("amin is now %u\n", amin);  fflush(stdout);
		}

		q = s - a * ascale;
		mq = q * -1;

		do
		{
			if (mq < 0)
			{
				if (Q[numR - map[abs(mq)] - 1]->len > 0)
				{
					ap = dequeue(Q[numR - map[abs(mq)] - 1]);
					//printf("dequeued %u from Q[%u](%d)\n", ap, R - map[abs(mq)] - 1, mq);

					if (ap == 0)
					{
						printf("dequeued %u from Q[%u](%d)\n", ap, numR - map[abs(mq)] - 1, mq);
						printf("a = %u\n", a);
						printf("ap = %u\n", ap);
						printf("s = %lu\n", s);
						fflush(stdout);
					}

					u = w * (a - ap) + q;

					if (u > umax)
					{
						//printf("(%u,%u)\n", 2 * ap, abs(q));
						int pa = ap - amin;
						int pb = abs(q);

						if (verbose & debug) {
							if (((((amin + pa) * ascale) + pb) == 3266117) ||
								((((amin + pa) * ascale) - pb) == 3266117))
							{
								printf("1: accumulating (%u,%u) = 3266117\n", pa, pb); fflush(stdout);
								printf("a  = %u\n", a);
								printf("ap  = %u\n", ap);
								printf("q  = %d\n", q);
								printf("mq = %d\n", mq);
								printf("amin = %d\n", amin);
							}
						}

						//printf("accumulate (%u,%u)\n", pa, pb); fflush(stdout);
						// accumulate the cross product  (zimmerman syntax).
						// page 342 in C&P

						if ((pb < 0) || (pb >= (U * D)))
						{
							printf("invalid pb = %d\n", pb);
						}

						if (rprime_map_U[pb] == 0)
						{
							printf("invalid distance %d\n", pb);
							printf("accumulate (%u,%u)\n", pa, pb); fflush(stdout);
							exit(-1);
						}

						if ((pa < 0) || (pa >= 2 * L))
						{
							printf("invalid Pa[%d]\n", pa);
							printf("accumulate (%u,%u)\n", pa, pb); fflush(stdout);
							exit(-1);
						}

                        CROSS_PRODUCT;

					}
					else
					{
						//printf("(%u,%u)\n", (a + ap), u);
						int pa = (a + ap) - 2 * amin;
						int pb = u;
						//printf("accumulate (%u,%u)\n", pa, pb); fflush(stdout);
						// accumulate the cross product  (zimmerman syntax).
						// page 342 in C&P
						if (verbose & debug) {
							if (((((2 * amin + pa) * w) + pb) == 3266117) ||
								((((2 * amin + pa) * w) - pb) == 3266117))
							{
								printf("2: accumulating (%u,%u) = 3266117\n", pa, pb); fflush(stdout);
								printf("a  = %u\n", a);
								printf("ap  = %u\n", ap);
								printf("q  = %d\n", q);
								printf("mq = %d\n", mq);
								printf("u = %d\n", u);
								printf("map[u] = %u\n", rprime_map_U[pb]);
								printf("amin = %d\n", amin);
							}
						}

						if (pa & 1)
						{
							pa = (pa - 1) / 2;
							pb -= D;
							if (pb < 0)
							{
								pa++;
								pb += 2 * D;
							}
						}
						else
						{
							pa >>= 1;
						}

						if ((pb < 0) || (pb >= (U * D)))
						{
							printf("invalid pb = %d\n", pb);
						}

						if (rprime_map_U[pb] == 0)
						{
							printf("invalid distance %d\n", pb);
							printf("accumulate (%u,%u)\n", pa, pb); fflush(stdout);
							exit(-1);
						}

						if ((pa < 0) || (pa >= 2 * L))
						{
							printf("invalid Pa[%d]\n", pa);
							printf("accumulate (%u,%u)\n", pa, pb); fflush(stdout);
							exit(-1);
						}

                        CROSS_PRODUCT;
					}
					work->paired++;
				}
				else
				{
					if (q < 0)
					{
						//printf("queueing %u in Q[%u](%d)\n", a, R - map[abs(q)] - 1, q);
						enqueue(Q[numR - map[abs(q)] - 1], a);
					}
					else
					{
						//printf("queueing %u in Q[%u](%d)\n", a, R + map[abs(q)], q);
						enqueue(Q[numR + map[abs(q)]], a);
					}
					u = 0;
				}
			}
			else if (mq > 0)
			{
				if (Q[numR + map[abs(mq)]]->len > 0)
				{
					ap = dequeue(Q[numR + map[abs(mq)]]);
					//printf("dequeued %u from Q[%u](%d)\n", ap, R + map[abs(mq)],mq);

					if (ap == 0)
					{
						printf("dequeued %u from Q[%u](%d)\n", ap, numR + map[abs(mq)], mq);
						printf("a = %u\n", a);
						printf("ap = %u\n", ap);
						printf("s = %lu\n", s);
						fflush(stdout);
					}

					u = w * (a - ap) + q;

					if (u > umax)
					{
						//printf("(%u,%u)\n", 2 * ap, abs(q));
						int pa = ap - amin;
						int pb = abs(q);
						//printf("accumulate (%u,%u)\n", pa, pb); fflush(stdout);
						// accumulate the cross product  (zimmerman syntax).
						// page 342 in C&P

						if (verbose & debug) {
							if (((((amin + pa) * ascale) + pb) == 3266117) ||
								((((amin + pa) * ascale) - pb) == 3266117))
							{
								printf("3: accumulating (%u,%u) = 3266117\n", pa, pb); fflush(stdout);
								printf("a  = %u\n", a);
								printf("ap  = %u\n", ap);
								printf("q  = %d\n", q);
								printf("mq = %d\n", mq);
								printf("amin = %d\n", amin);
							}
						}

						if ((pb < 0) || (pb >= (U * D)))
						{
							printf("invalid pb = %d\n", pb);
						}

						if (rprime_map_U[pb] == 0)
						{
							printf("invalid distance %d\n", pb);
							printf("accumulate (%u,%u)\n", pa, pb); fflush(stdout);
							exit(-1);
						}

						if ((pa < 0) || (pa >= 2 * L))
						{
							printf("invalid Pa[%d]\n", pa);
							printf("accumulate (%u,%u)\n", pa, pb); fflush(stdout);
							exit(-1);
						}

                        CROSS_PRODUCT;
					}
					else
					{
						//printf("(%u,%u)\n", (a + ap), u);
						int pa = (a + ap) - 2 * amin;
						int pb = u;
						//printf("accumulate (%u,%u)\n", pa, pb); fflush(stdout);
						// accumulate the cross product  (zimmerman syntax).
						// page 342 in C&P
						if (verbose & debug) {
							if (((((2 * amin + pa) * w) + pb) == 3266117) ||
								((((2 * amin + pa) * w) - pb) == 3266117))
							{
								printf("4: accumulating (%u,%u) = 3266117\n", pa, pb); fflush(stdout);
								printf("a  = %u\n", a);
								printf("ap  = %u\n", ap);
								printf("q  = %d\n", q);
								printf("mq = %d\n", mq);
								printf("u = %d\n", u);
								printf("amin = %d\n", amin);
							}
						}

						if (pa & 1)
						{
							pa = (pa - 1) / 2;
							pb -= D;
							if (pb < 0)
							{
								pa++;
								pb += 2 * D;
							}
						}
						else
						{
							pa >>= 1;
						}

						if ((pb < 0) || (pb >= (U * D)))
						{
							printf("invalid pb = %d\n", pb);
						}

						if (rprime_map_U[pb] == 0)
						{
							printf("invalid distance %d\n", pb);
							printf("accumulate (%u,%u)\n", pa, pb); fflush(stdout);
							exit(-1);
						}

						if ((pa < 0) || (pa >= 2 * L))
						{
							printf("invalid Pa[%d]\n", pa);
							printf("accumulate (%u,%u)\n", pa, pb); fflush(stdout);
							exit(-1);
						}

                        CROSS_PRODUCT;
					}
					work->paired++;
				}
				else
				{
					if (q < 0)
					{
						//printf("queueing %u in Q[%u](%d)\n", a, R - map[abs(q)] - 1, q);
						enqueue(Q[numR - map[abs(q)] - 1], a);
					}
					else
					{
						//printf("queueing %u in Q[%u](%d)\n", a, R + map[abs(q)], q);
						enqueue(Q[numR + map[abs(q)]], a);
					}
					u = 0;
				}
			}

		} while (u > umax);

	}

	//printf("main loop done, flushing Queues\n"); fflush(stdout);
	j = 0;

	for (i = 0; i < numR; i++)
	{
		while (Q[i]->len > 0)
		{
			//printf("(%u,%u)\n", (2 * dequeue(Q[i])), rmap[numR - i - 1]);
			int pa = dequeue(Q[i]) - amin;
			int pb = rmap[numR - i - 1];
			//printf("accumulate (%u,%u)\n", pa, pb); fflush(stdout);
			// accumulate the cross product  (zimmerman syntax).
			// page 342 in C&P
            CROSS_PRODUCT;
			work->paired++;
		}
	}
	for (i = numR; i < 2 * numR; i++)
	{
		while (Q[i]->len > 0)
		{
			//printf("(%u,%u)\n", (2 * dequeue(Q[i])), rmap[i - numR]);
			int pa = dequeue(Q[i]) - amin;
			int pb = rmap[i - numR];
			//printf("accumulate (%u,%u)\n", pa, pb); fflush(stdout);
			// accumulate the cross product  (zimmerman syntax).
			// page 342 in C&P
            CROSS_PRODUCT;
			work->paired++;
		}
	}

    //printf("computed %u new a ranges (%d each)\n", aranges, L - U);

	work->amin = amin;
	work->last_pid = pid;

	return;
}


#ifdef PARAM1


/* ECM stage 1 in batch mode, for initial point (x:z) with small coordinates,
   such that x and z fit into a mp_limb_t.
   For example we can start with (x=2:y=1) with the curve by^2 = x^3 + ax^2 + x
   with a = 4d-2 and b=16d+2, then we have to multiply by d=(a+2)/4 in the
   duplicates.
   With the change of variable x=b*X, y=b*Y, this curve becomes:
   Y^2 = X^3 + a/b*X^2 + 1/b^2*X.
*/

#define MAX_HEIGHT 32

#if defined(_WIN32)
/* Due to a limitation in GMP on 64-bit Windows, should also
    affect 32-bit Windows, sufficient memory cannot be allocated
    for the batch product s when using primes larger than the following */
#define MAX_B1_BATCH 3124253146UL
#else
/* nth_prime(2^(MAX_HEIGHT-1)) */
#define MAX_B1_BATCH 50685770167ULL
#endif

unsigned int compute_s(mpz_t s, uint64_t *primes, uint64_t B1)
{
    mpz_t acc[MAX_HEIGHT]; /* To accumulate products of prime powers */
    mpz_t ppz;
    unsigned int i, j, it;
    uint64_t pi = 2, pp, maxpp, qi;

    for (j = 0; j < MAX_HEIGHT; j++)
        mpz_init(acc[j]); /* sets acc[j] to 0 */
    mpz_init(ppz);

    i = 0;
    while (pi <= B1)
    {
        pp = qi = pi;
        maxpp = B1 / qi;

        while (pp <= maxpp)
            pp *= qi;

        mpz_set_ui(ppz, pp);

        if ((i & 1) == 0)
            mpz_set(acc[0], ppz);
        else
            mpz_mul(acc[0], acc[0], ppz);

        j = 0;
        /* We have accumulated i+1 products so far. If bits 0..j of i are all
           set, then i+1 is a multiple of 2^(j+1). */
        while ((i & (1 << j)) != 0)
        {
            /* we use acc[MAX_HEIGHT-1] as 0-sentinel below, thus we need
               j+1 < MAX_HEIGHT-1 */
            if ((i & (1 << (j + 1))) == 0) /* i+1 is not multiple of 2^(j+2),
                                              thus add[j+1] is "empty" */
                mpz_swap(acc[j + 1], acc[j]); /* avoid a copy with mpz_set */
            else
                mpz_mul(acc[j + 1], acc[j + 1], acc[j]); /* accumulate in acc[j+1] */
            mpz_set_ui(acc[j], 1);
            j++;
        }

        i++;
        pi = primes[i];
    }
    it = i;

    for (mpz_set(s, acc[0]), j = 1; mpz_cmp_ui(acc[j], 0) != 0; j++)
        mpz_mul(s, s, acc[j]);

    for (i = 0; i < MAX_HEIGHT; i++)
        mpz_clear(acc[i]);
    mpz_clear(ppz);

    return it;
}



/* R <- S*m/B mod modulus where m fits in a mp_limb_t.
   Here S (w in dup_add_batch1) is the result of a subtraction,
   thus with the notations from http://www.loria.fr/~zimmerma/papers/norm.pdf
   we have S < 2 \alpha N.
   Then R < (2 \alpha N \beta + \beta N) = (2 \alpha + 1) N.
   This result R is used in an addition with u being the result of a squaring
   thus u < \alpha N, which gives a result < (3 \alpha + 1) N.
   Finally this result is used in a multiplication with another operand less
   than 2 \alpha N, thus we want:
   ((2 \alpha) (3 \alpha + 1) N^2 + \beta N)/\beta \leq \alpha N, i.e.,
   2 \alpha (3 \alpha + 1) \varepsilon + 1 \leq \alpha
   This implies \varepsilon \leq 7/2 - sqrt(3)/2 ~ 0.0359, in which case
   we can take \alpha = 2/3*sqrt(3)+1 ~ 2.1547.
   In that case no adjustment is needed in mpresn_mul_1.
   However we prefer to keep the adjustment here, to allow a larger set of
   inputs (\varepsilon \leq 1/16 = 0.0625 instead of 0.0359).
*/
void
vecmulmod_1(vec_bignum_t * S, uint64_t *m, vec_bignum_t * R, vec_bignum_t * modulus, vec_bignum_t * t, vec_monty_t *mdata)
{
    //mp_ptr t1 = PTR(modulus->temp1);
    //mp_ptr t2 = PTR(modulus->temp2);
    //mp_size_t n = ABSIZ(modulus->orig_modulus);
    //mp_limb_t q;
    //
    //{
    //    t1[n] = mpn_mul_1(t1, PTR(S), n, m);
    //    q = t1[0] * modulus->Nprim[0];
    //    t2[n] = mpn_mul_1(t2, PTR(modulus->orig_modulus), n, q);
    //
    //    q = mpn_add_n(PTR(R), t1 + 1, t2 + 1, n);
    //    q += mpn_add_1(PTR(R), PTR(R), n, t1[0] != 0);
    //
    //    while (q != 0)
    //        q -= mpn_sub_n(PTR(R), PTR(R), PTR(modulus->orig_modulus), n);
    //}
    vecmulmod52_1(S, m, R, modulus, t, mdata);
}


/* (x1:z1) <- 2(x1:z1)
   (x2:z2) <- (x1:z1) + (x2:z2)
   assume (x2:z2) - (x1:z1) = (2:1)
   Uses 4 modular multiplies and 4 modular squarings.
   Inputs are x1, z1, x2, z2, d, n.
   Use two auxiliary variables: t, w (it seems using one only is not possible
   if all mpresn_mul and mpresn_sqr calls don't overlap input and output).

   In the batch 1 mode, we pass d_prime such that the actual d is d_prime/beta.
   Since beta is a square, if d_prime is a square (on 64-bit machines),
   so is d.
   In mpresn_mul_1, we multiply by d_prime = beta*d and divide by beta.
*/
static void
dup_add_batch1(vec_bignum_t *x1, vec_bignum_t *z1, vec_bignum_t *x2, vec_bignum_t *z2,
    vec_bignum_t *t, vec_bignum_t *w, vec_bignum_t *s, uint64_t *d_prime, vec_bignum_t *n, vec_monty_t *mdata)
{
    //bignum *t2;
    //t2 = vecInit();
    //memcpy(t2->data, d_prime, VECLEN * sizeof(base_t));

    /* active: x1 z1 x2 z2 */
    //mpresn_addsub(w, z1, x1, z1, n); /* w = x1+z1, z1 = x1-z1 */
    vecaddsubmod_ptr(x1, z1, w, z1, n);

    /* active: w z1 x2 z2 */
    //mpresn_addsub(x1, x2, x2, z2, n); /* x1 = x2+z2, x2 = x2-z2 */
    vecaddsubmod_ptr(x2, z2, x1, x2, n);
    /* active: w z1 x1 x2 */

    //mpresn_mul(z2, w, x2, n); /* w = (x1+z1)(x2-z2) */
    vecmulmod_ptr(w, x2, z2, n, s, mdata);
    /* active: w z1 x1 z2 */
    //mpresn_mul(x2, z1, x1, n); /* x2 = (x1-z1)(x2+z2) */
    vecmulmod_ptr(z1, x1, x2, n, s, mdata);
    /* active: w z1 x2 z2 */
    //mpresn_sqr(t, z1, n);    /* t = (x1-z1)^2 */
    vecsqrmod_ptr(z1, t, n, s, mdata);
    /* active: w t x2 z2 */
    //mpresn_sqr(z1, w, n);    /* z1 = (x1+z1)^2 */
    vecsqrmod_ptr(w, z1, n, s, mdata);
    /* active: z1 t x2 z2 */

    //mpresn_mul(x1, z1, t, n); /* xdup = (x1+z1)^2 * (x1-z1)^2 */
    vecmulmod_ptr(z1, t, x1, n, s, mdata);
    /* active: x1 z1 t x2 z2 */

    //mpresn_sub(w, z1, t, n);   /* w = (x1+z1)^2 - (x1-z1)^2 */
    vecsubmod_ptr(z1, t, w, n);
    /* active: x1 w t x2 z2 */

    //mpresn_mul_1(z1, w, d_prime, n); /* z1 = d * ((x1+z1)^2 - (x1-z1)^2) */
    vecmulmod_1(w, d_prime, z1, n, s, mdata);
    //vecmulmod_ptr(w, t2, z1, n, s, mdata);
    /* active: x1 z1 w t x2 z2 */

    //mpresn_add(t, t, z1, n);  /* t = (x1-z1)^2 - d* ((x1+z1)^2 - (x1-z1)^2) */
    vecaddmod_ptr(z1, t, t, n);
    //vecsubmod_ptr(z1, t, t, n);

    /* active: x1 w t x2 z2 */
    //mpresn_mul(z1, w, t, n); /* zdup = w * [(x1-z1)^2 - d* ((x1+z1)^2 - (x1-z1)^2)] */
    vecmulmod_ptr(w, t, z1, n, s, mdata);
    /* active: x1 z1 x2 z2 */

    //mpresn_addsub(w, z2, x2, z2, n);
    vecaddsubmod_ptr(x2, z2, w, z2, n);
    /* active: x1 z1 w z2 */

    //mpresn_sqr(x2, w, n);
    vecsqrmod_ptr(w, x2, n, s, mdata);
    /* active: x1 z1 x2 z2 */
    //mpresn_sqr(w, z2, n);
    vecsqrmod_ptr(z2, w, n, s, mdata);
    /* active: x1 z1 x2 w */
    //mpresn_add(z2, w, w, n);
    vecaddmod_ptr(w, w, z2, n);

    //vecFree(t2);
}



/* Input: x is initial point
          A is curve parameter in Montgomery's form:
          g*y^2*z = x^3 + a*x^2*z + x*z^2
          n is the number to factor
      B1 is the stage 1 bound
   Output: If a factor is found, it is returned in x.
           Otherwise, x contains the x-coordinate of the point computed
           in stage 1 (with z coordinate normalized to 1).
       B1done is set to B1 if stage 1 completed normally,
       or to the largest prime processed if interrupted, but never
       to a smaller value than B1done was upon function entry.
   Return value: ECM_FACTOR_FOUND_STEP1 if a factor, otherwise
           ECM_NO_FACTOR_FOUND
*/
/*
For now we don't take into account go stop_asap and chkfilename
*/
int vec_ecm_stage1_batch(mpz_t f, ecm_work *work, ecm_pt *P, vec_bignum_t * A, 
    vec_bignum_t * n, uint64_t B1, mpz_t s, vec_monty_t *mdata)
{
    uint64_t *d_1;
    mpz_t d_2;

    vec_bignum_t *x1, *z1, *x2, *z2;
    uint64_t i;
    vec_bignum_t *t, *u;
    int ret = 0;

    x1 = vecInit();
    z1 = vecInit();
    x2 = vecInit();
    z2 = vecInit();
    t  = vecInit();
    u  = vecInit();

    d_1 = (uint64_t *)xmalloc_align(VECLEN * sizeof(uint64_t));

    vecCopy(P->X, x1);
    vecCopy(P->Z, z1);
    vecCopy(work->pt2.X, x2);
    vecCopy(work->pt2.Z, z2);
    memcpy(d_1, work->s->data, VECLEN * sizeof(base_t));

    //printf("d_1 = ");
    //for (i = 0; i < VECLEN; i++)
    //    printf("%lu, ", d_1[i]);
    //printf("\n");

    /* initialize P */
    //mpres_set(x1, x, n);
    //mpres_set_ui(z1, 1, n); /* P1 <- 1P */

    /* Compute d=(A+2)/4 from A and d'=B*d thus d' = 2^(GMP_NUMB_BITS-2)*(A+2) */
    //mpres_get_z(u, A, n);
    //mpz_add_ui(u, u, 2);
    //mpz_mul_2exp(u, u, GMP_NUMB_BITS - 2);
    //mpres_set_z_for_gcd(u, u, n); /* reduces u mod n */
    //if (mpz_size(u) > 1)
    //{
    //    mpres_get_z(u, A, n);
    //    gmp_fprintf(stderr,
    //        "Error, 2^%d*(A+2) should fit in a mp_limb_t, A=%Zd\n",
    //        GMP_NUMB_BITS - 2, u);
    //    return -1;
    //}
    //d_1 = mpz_getlimbn(u, 0);

    ///* Compute 2P : no need to duplicate P, the coordinates are simple. */
    //mpres_set_ui(x2, 9, n);
    ///* here d = d_1 / GMP_NUMB_BITS */
    ///* warning: mpres_set_ui takes an unsigned long which has only 32 bits
    //    on Windows, while d_1 might have 64 bits */
    //if (mpz_size(u) != 1 || mpz_getlimbn(u, 0) != d_1_1)
    //{
    //    printf("problem with mpres_set_ui\n");
    //    exit(-2);
    //}
    //mpres_set_z(z2, u, n);
    //mpres_div_2exp(z2, z2, GMP_NUMB_BITS, n);
    //
    //mpres_mul_2exp(z2, z2, 6, n);
    //mpres_add_ui(z2, z2, 8, n); /* P2 <- 2P = (9 : : 64d+8) */

    /* invariant: if j represents the upper bits of s,
       then P1 = j*P and P2=(j+1)*P */

    //mpresn_pad(x1, n);
    //mpresn_pad(z1, n);
    //mpresn_pad(x2, n);
    //mpresn_pad(z2, n);

    /* now perform the double-and-add ladder */
    for (i = mpz_sizeinbase(s, 2) - 1; i-- > 0;)
    {
        if (mpz_tstbit(s, i) == 0) /* (j,j+1) -> (2j,2j+1) */
            /* P2 <- P1+P2    P1 <- 2*P1 */
            dup_add_batch1(x1, z1, x2, z2, t, u, mdata->mtmp1, d_1, n, mdata);
        else /* (j,j+1) -> (2j+1,2j+2) */
            /* P1 <- P1+P2     P2 <- 2*P2 */
            dup_add_batch1(x2, z2, x1, z1, t, u, mdata->mtmp1, d_1, n, mdata);
    }

    //mpresn_unpad(x1);
    //mpresn_unpad(z1);

    // this sets x = x1/z1, but our vececm expects
    // them to be separate.
    //if (!mpres_invert(u, z1, n)) /* Factor found? */
    //{
    //    mpres_gcd(f, z1, n);
    //    ret = 1;
    //}
    //mpres_mul(x, x1, u, n);
    //vecmulmod_ptr(x1, u, x, mdata->mtmp1, n, mdata);

    vecCopy(x1, P->X);
    vecCopy(z1, P->Z);

    vecFree(x1);
    vecFree(z1);
    vecFree(x2);
    vecFree(z2);
    vecFree(t);
    vecFree(u);
    align_free(d_1);

    return ret;
}

#endif

int vec_check_factor(mpz_t Z, mpz_t n, mpz_t f)
{
    //gmp_printf("checking point Z = %Zx against input N = %Zx\n", Z, n);
    mpz_gcd(f, Z, n);

	if (mpz_cmp_ui(f, 1) > 0)
	{
		if (mpz_cmp(f, n) == 0)
		{
            mpz_set_ui(f, 0);
			return 0;
		}
		return 1;
	}
	return 0;
}


