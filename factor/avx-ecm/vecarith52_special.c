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


Copyright (c) 2018 by The Mayo Clinic, though its Special Purpose
 Processor Development Group (SPPDG). All Rights Reserved Worldwide.
 Licensed under the Apache License, Version 2.0 (the "License"); you may
 not use this file except in compliance with the License. You may obtain
 a copy of the License at http://www.apache.org/licenses/LICENSE-2.0.
Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied,
           including conditions of title, non-infringement, merchantability,
          or fitness for a particular purpose
 See the License for the specific language governing permissions and
 limitations under the License.
This file is a snapshot of a work in progress, originated by Mayo
 Clinic SPPDG.
*/

#include "avx_ecm.h"

#ifdef _MSC_VER
//#define USE_AVX512F
#endif

#ifdef USE_AVX512F
#include <immintrin.h>

//#define PRINT_DEBUG 2
//#define NPRINT_DEBUG 0

// test cases:
// is stage 2 broken for mersennes? - yes looks that way

// ./yafu "ecm(2^823-1,8)" -v -v -B1ecm 3000000 -sigma 4248094344 finds factor 1460915248436556406607
// ./yafu "ecm(2^823-1,8)" -v -v -B1ecm 3000000 -sigma 2102985689 finds both factors 1460915248436556406607 and 1534086200463688788034864584049
// ./yafu "ecm((2^1063+1)/3,8)" -v -v -B1ecm 3000000 -sigma 3299028348 finds both factors 114584129081 and 26210488518118323164267329859
// ./yafu "ecm((2^827+1)/3,8)" -v -v -B1ecm 28000000 -B2ecm 760000000 -sigma 4029008539 fails to find factor!
// ./yafu "ecm((2^773-1)/6864241/9461521,8)" -v -v -B1ecm 2500000 -B2ecm 1800000000 -sigma 8170945836124664 fails to find factor!


void vecmulmod52_mersenne_416(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    int i, j;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t NBLOCKS = mdata->NBLOCKS;
    // needed in loops
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19
    __m512i c00, c01, c02, c03, c04, c05, c06, c07,
        c08, c09, c10, c11, c12, c13, c14, c15;

    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry;

    // deal with the sign
    c->size = NWORDS;
    c->signmask = a->signmask ^ b->signmask;

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

#ifdef DEBUG_MERSENNE
    printf("commencing vecmulmod52_mersenne\n");
    print_vechexbignum(a, "input a: ");
    print_vechexbignum(b, "input b: ");
    print_vechexbignum(s, "input s: ");
#endif

    // first half mul
    //for (i = 0; i < NBLOCKS; i++)
    i = 0;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (j = i; j > 0; j--)
        {
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum
        a0 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 0) * VECLEN);
        a1 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + 3 * VECLEN);
        b1 = _mm512_load_epi64(b->data + 2 * VECLEN);
        b2 = _mm512_load_epi64(b->data + 1 * VECLEN);
        b3 = _mm512_load_epi64(b->data + 0 * VECLEN);

        {
            prod5_ld = _mm512_cvtepu64_pd(a0);
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(a2);
            prod3_ld = _mm512_cvtepu64_pd(b2);
            prod4_ld = _mm512_cvtepu64_pd(b3);

            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod1_hd));  // a1 * b3 -> to te2/3
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));  // a2 * b3 -> to te4/5 
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a0 * b2 -> to te2/3
            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod4_hd));  // a0 * b3 -> to te0/1

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod1_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod4_ld));
        }

        {
            prod5_ld = _mm512_cvtepu64_pd(a0);
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(b0);
            prod3_ld = _mm512_cvtepu64_pd(b1);
            prod4_ld = _mm512_cvtepu64_pd(b2);

            prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_hd = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod2_hd));
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod3_hd));
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod4_hd));

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_ld = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod2_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod3_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod4_ld));
        }

        {
            prod1_ld = _mm512_cvtepu64_pd(a2);
            prod2_ld = _mm512_cvtepu64_pd(a3);
            prod3_ld = _mm512_cvtepu64_pd(b2);
            prod4_ld = _mm512_cvtepu64_pd(b3);

            prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod2_hd));

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

            prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod2_ld));
        }

        // subtract out all of the bias at once.  these are
        // counts of how many times we have mul_accumulated into each column
        // during the block a*b and final triangular a*b.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 2,
            i * 4 + 3,
            i * 4 + 4);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 2,
            i * 4 + 3,
            i * 4 + 4);


        // now, do a carry-propagating column summation and store the results.
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        // i = 3, a3..0 = c[0..3]
        // i = 4, a3..0 = c[4..7]
        // i = 5, a3..0 = c[8..11]
        c00 = a3;
        c01 = a2;
        c02 = a1;
        c03 = a0;
    }

    i = 1;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (j = i; j > 0; j--)
        {
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum
        a0 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 0) * VECLEN);
        a1 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + 3 * VECLEN);
        b1 = _mm512_load_epi64(b->data + 2 * VECLEN);
        b2 = _mm512_load_epi64(b->data + 1 * VECLEN);
        b3 = _mm512_load_epi64(b->data + 0 * VECLEN);

        {
            prod5_ld = _mm512_cvtepu64_pd(a0);
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(a2);
            prod3_ld = _mm512_cvtepu64_pd(b2);
            prod4_ld = _mm512_cvtepu64_pd(b3);

            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod1_hd));  // a1 * b3 -> to te2/3
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));  // a2 * b3 -> to te4/5 
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a0 * b2 -> to te2/3
            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod4_hd));  // a0 * b3 -> to te0/1

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod1_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod4_ld));
        }

        {
            prod5_ld = _mm512_cvtepu64_pd(a0);
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(b0);
            prod3_ld = _mm512_cvtepu64_pd(b1);
            prod4_ld = _mm512_cvtepu64_pd(b2);

            prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_hd = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod2_hd));
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod3_hd));
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod4_hd));

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_ld = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod2_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod3_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod4_ld));
        }

        {
            prod1_ld = _mm512_cvtepu64_pd(a2);
            prod2_ld = _mm512_cvtepu64_pd(a3);
            prod3_ld = _mm512_cvtepu64_pd(b2);
            prod4_ld = _mm512_cvtepu64_pd(b3);

            prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod2_hd));

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

            prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod2_ld));
        }

        // subtract out all of the bias at once.  these are
        // counts of how many times we have mul_accumulated into each column
        // during the block a*b and final triangular a*b.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 2,
            i * 4 + 3,
            i * 4 + 4);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 2,
            i * 4 + 3,
            i * 4 + 4);


        // now, do a carry-propagating column summation and store the results.
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        // i = 3, a3..0 = c[0..3]
        // i = 4, a3..0 = c[4..7]
        // i = 5, a3..0 = c[8..11]
        c04 = a3;
        c05 = a2;
        c06 = a1;
        c07 = a0;
    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(s, "after lo half:");
#endif

    // second half mul
    //for (i = NBLOCKS; i < 2 * NBLOCKS; i++)
    i = 2;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (j = i - NBLOCKS + 1; j < NBLOCKS; j++)
        {
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }


        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN);

        {
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(a2);
            prod3_ld = _mm512_cvtepu64_pd(b0);
            prod4_ld = _mm512_cvtepu64_pd(b1);

            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));  // a1*b0 -> t0/1
            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod4_hd));  // a2*b1 -> t0/1
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a2*b0 -> t2/3

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod4_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
        }

        {
            prod1_ld = _mm512_cvtepu64_pd(a3);
            prod2_ld = _mm512_cvtepu64_pd(b2);
            prod3_ld = _mm512_cvtepu64_pd(b1);
            prod4_ld = _mm512_cvtepu64_pd(b0);

            prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));  // a3*b2 -> t0/1
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a3*b1 -> t2/3
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod4_hd));  // a3*b0 -> t4/5

            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod4_ld));
        }


        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            (2 * NBLOCKS - i - 1) * 4 + 3,
            (2 * NBLOCKS - i - 1) * 4 + 2,
            (2 * NBLOCKS - i - 1) * 4 + 1,
            (2 * NBLOCKS - i - 1) * 4 + 0);
        SUB_BIAS_LO(
            (2 * NBLOCKS - i - 1) * 4 + 3,
            (2 * NBLOCKS - i - 1) * 4 + 2,
            (2 * NBLOCKS - i - 1) * 4 + 1,
            (2 * NBLOCKS - i - 1) * 4 + 0);

        // final column accumulation
        {
            ACCUM_SHIFT(a3, te0, te1);
            ACCUM_SHIFT(a2, te2, te3);
            ACCUM_SHIFT(a1, te4, te5);
            ACCUM_SHIFT(a0, te6, te7);

            a3 = _mm512_and_epi64(vlmask, a3);
            a2 = _mm512_and_epi64(vlmask, a2);
            a1 = _mm512_and_epi64(vlmask, a1);
            a0 = _mm512_and_epi64(vlmask, a0);

            // i = 3, a3..0 = c[0..3]
            // i = 4, a3..0 = c[4..7]
            // i = 5, a3..0 = c[8..11]
            //_mm512_store_epi64(c->data + (i * BLOCKWORDS + 0) * VECLEN, a3);
            //_mm512_store_epi64(c->data + (i * BLOCKWORDS + 1) * VECLEN, a2);
            //_mm512_store_epi64(c->data + (i * BLOCKWORDS + 2) * VECLEN, a1);
            //_mm512_store_epi64(c->data + (i * BLOCKWORDS + 3) * VECLEN, a0);

            c08 = a3;
            c09 = a2;
            c10 = a1;
            c11 = a0;

        }

    }

    i = 3;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (j = i - NBLOCKS + 1; j < NBLOCKS; j++)
        {
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }


        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN);

        {
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(a2);
            prod3_ld = _mm512_cvtepu64_pd(b0);
            prod4_ld = _mm512_cvtepu64_pd(b1);

            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));  // a1*b0 -> t0/1
            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod4_hd));  // a2*b1 -> t0/1
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a2*b0 -> t2/3

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod4_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
        }

        {
            prod1_ld = _mm512_cvtepu64_pd(a3);
            prod2_ld = _mm512_cvtepu64_pd(b2);
            prod3_ld = _mm512_cvtepu64_pd(b1);
            prod4_ld = _mm512_cvtepu64_pd(b0);

            prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));  // a3*b2 -> t0/1
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a3*b1 -> t2/3
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod4_hd));  // a3*b0 -> t4/5

            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod4_ld));
        }


        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            (2 * NBLOCKS - i - 1) * 4 + 3,
            (2 * NBLOCKS - i - 1) * 4 + 2,
            (2 * NBLOCKS - i - 1) * 4 + 1,
            (2 * NBLOCKS - i - 1) * 4 + 0);
        SUB_BIAS_LO(
            (2 * NBLOCKS - i - 1) * 4 + 3,
            (2 * NBLOCKS - i - 1) * 4 + 2,
            (2 * NBLOCKS - i - 1) * 4 + 1,
            (2 * NBLOCKS - i - 1) * 4 + 0);

        // final column accumulation
        {
            ACCUM_SHIFT(a3, te0, te1);
            ACCUM_SHIFT(a2, te2, te3);
            ACCUM_SHIFT(a1, te4, te5);
            ACCUM_SHIFT(a0, te6, te7);

            a3 = _mm512_and_epi64(vlmask, a3);
            a2 = _mm512_and_epi64(vlmask, a2);
            a1 = _mm512_and_epi64(vlmask, a1);
            a0 = _mm512_and_epi64(vlmask, a0);

            // i = 3, a3..0 = c[0..3]
            // i = 4, a3..0 = c[4..7]
            // i = 5, a3..0 = c[8..11]
            //_mm512_store_epi64(c->data + (i * BLOCKWORDS + 0) * VECLEN, a3);
            //_mm512_store_epi64(c->data + (i * BLOCKWORDS + 1) * VECLEN, a2);
            //_mm512_store_epi64(c->data + (i * BLOCKWORDS + 2) * VECLEN, a1);
            //_mm512_store_epi64(c->data + (i * BLOCKWORDS + 3) * VECLEN, a0);

            c12 = a3;
            c13 = a2;
            c14 = a1;
            c15 = a0;

        }

    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(s, "after hi half:");
#endif

    // reduce by adding hi to lo.  first right shift hi into output.
    __m512i vbshift, vbpshift;
    int bshift = mdata->nbits % 52;
    int wshift = mdata->nbits / 52;

    if (mdata->use_vnbits)
    {
        // for different lanes computing with different N.
        // right now that is only for avxppm1
        int bshift_a[VECLEN];
        for (i = 0; i < VECLEN; i++)
        {
            bshift_a[i] = mdata->vnbits[i] % 52;
        }

        vbshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[7])) - 1ULL);

        vbpshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])),
            (1ULL << (uint64_t)(bshift_a[1])),
            (1ULL << (uint64_t)(bshift_a[2])),
            (1ULL << (uint64_t)(bshift_a[3])),
            (1ULL << (uint64_t)(bshift_a[4])),
            (1ULL << (uint64_t)(bshift_a[5])),
            (1ULL << (uint64_t)(bshift_a[6])),
            (1ULL << (uint64_t)(bshift_a[7])));

        vec_bignum52_mask_rshift_vn(s, c, mdata->vnbits, 0xff);
    }
    else
    {
        _mm512_store_epi64(s->data + 0 * VECLEN, c00);
        _mm512_store_epi64(s->data + 1 * VECLEN, c01);
        _mm512_store_epi64(s->data + 2 * VECLEN, c02);
        _mm512_store_epi64(s->data + 3 * VECLEN, c03);
        _mm512_store_epi64(s->data + 4 * VECLEN, c04);
        _mm512_store_epi64(s->data + 5 * VECLEN, c05);
        _mm512_store_epi64(s->data + 6 * VECLEN, c06);
        _mm512_store_epi64(s->data + 7 * VECLEN, c07);
        _mm512_store_epi64(s->data + 8 * VECLEN, c08);
        _mm512_store_epi64(s->data + 9 * VECLEN, c09);
        _mm512_store_epi64(s->data + 10 * VECLEN, c10);
        _mm512_store_epi64(s->data + 11 * VECLEN, c11);
        _mm512_store_epi64(s->data + 12 * VECLEN, c12);
        _mm512_store_epi64(s->data + 13 * VECLEN, c13);
        _mm512_store_epi64(s->data + 14 * VECLEN, c14);
        _mm512_store_epi64(s->data + 15 * VECLEN, c15);

        vbshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL);
        vbpshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)));
        vec_bignum52_mask_rshift_n(s, c, mdata->nbits, 0xff);


        //int wshift = mdata->nbits / 52;
        //int bshift = mdata->nbits % 52;
        //__m512i lowmask = _mm512_set1_epi64((1ULL << (uint64_t)bshift) - 1ULL);
        //
        //// split off the high half into the top words.  
        //c15 = _mm512_or_epi64(_mm512_slli_epi64(c15, (DIGITBITS - bshift)), _mm512_srli_epi64(c14, bshift));
        //c14 = _mm512_or_epi64(_mm512_slli_epi64(c14, (DIGITBITS - bshift)), _mm512_srli_epi64(c13, bshift));
        //c13 = _mm512_or_epi64(_mm512_slli_epi64(c13, (DIGITBITS - bshift)), _mm512_srli_epi64(c12, bshift));
        //c12 = _mm512_or_epi64(_mm512_slli_epi64(c12, (DIGITBITS - bshift)), _mm512_srli_epi64(c11, bshift));
        //c11 = _mm512_or_epi64(_mm512_slli_epi64(c11, (DIGITBITS - bshift)), _mm512_srli_epi64(c10, bshift));
        //c10 = _mm512_or_epi64(_mm512_slli_epi64(c10, (DIGITBITS - bshift)), _mm512_srli_epi64(c09, bshift));
        //c09 = _mm512_or_epi64(_mm512_slli_epi64(c09, (DIGITBITS - bshift)), _mm512_srli_epi64(c08, bshift));
        //c08 = _mm512_or_epi64(_mm512_slli_epi64(c08, (DIGITBITS - bshift)), _mm512_srli_epi64(c07, bshift));
    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(c, "hi part:");
    print_vechexbignum(s, "lo part:");
#endif

    if (mdata->isMersenne > 1)
    {
        // this number has been identified as a "pseudo"-Mersenne: 2^n-c, with
        // 'c' small.  Fast reduction is still possible.

        // multiply c * hi
        b0 = _mm512_set1_epi64(mdata->isMersenne);
        vecmul52_1(c, b0, c, n, s, mdata);

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after mul_1 with c:");
#endif

        // clear any hi-bits in the lo result
        a1 = _mm512_load_epi64(s->data + wshift * VECLEN);
        _mm512_store_epi64(s->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(s, "cleared low part:");
#endif

        // add c * hi + lo
        scarry = 0;
        for (i = 0; i <= NWORDS; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            b0 = _mm512_load_epi64(s->data + i * VECLEN);
            a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "lo + hi * c:");
#endif

        // now we need to do it again, but the hi part is just a single word so
        // we only have a single-precision vecmul
        //vec_bignum52_mask_rshift_n(c, s, mdata->nbits, 0xff);
        //a0 = _mm512_load_epi64(s->data);

        // It's possible these hi bits are located in two words
        a0 = _mm512_load_epi64(c->data + wshift * VECLEN);
        a0 = _mm512_srli_epi64(a0, bshift);
        a1 = _mm512_load_epi64(c->data + (wshift + 1) * VECLEN);
        a0 = _mm512_or_epi64(a0, _mm512_slli_epi64(a1, 52 - bshift));

        //print_regvechex(a0, 0, "hi part:");

        b0 = _mm512_set1_epi64(mdata->isMersenne);

        //print_regvechex(b0, 0, "c:");

        VEC_MUL_LOHI_PD(a0, b0, acc_e0, acc_e1);

        // clear any hi-bits now that we have the multiplier
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_set1_epi64(0));
        }

        // now add that into the lo part again.  Here we add in 
        // the two words of the hi-mul result.
        b0 = _mm512_load_epi64(c->data + 0 * VECLEN);
        b1 = _mm512_load_epi64(c->data + 1 * VECLEN);
        b2 = zero;
        b3 = zero;
        acc_e0 = _mm512_add_epi64(acc_e0, b0);
        VEC_CARRYPROP_LOHI(acc_e0, b2);         // lo word and carry for 1st add.
        acc_e1 = _mm512_add_epi64(acc_e1, b1);
        VEC_CARRYPROP_LOHI(acc_e1, b3);         // lo word and carry for 2nd add.
        acc_e1 = _mm512_add_epi64(acc_e1, b2);  // add 1st carry to 2nd lo word.
        VEC_CARRYPROP_LOHI(acc_e1, b3);         // that could have a carry, add it to b3.

        // save first two result words 
        _mm512_store_epi64(c->data + 0 * VECLEN, acc_e0);
        _mm512_store_epi64(c->data + 1 * VECLEN, acc_e1);

        // propagate the carry through the rest of the lo result.
        b2 = _mm512_load_epi64(c->data + 2 * VECLEN);
        scarry = 0;
        a0 = _mm512_adc_epi52(b2, scarry, b3, &scarry);
        _mm512_store_epi64(c->data + 2 * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i = 3; (i < wshift) && (scarry > 0); i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "lo + hi * c:");
#endif

        // it is unlikely, but possible that we still have a final carry to propagate
        if (scarry > 0)
        {
            printf("unlikely add\n");
            b0 = _mm512_set1_epi64(mdata->isMersenne);
            a0 = _mm512_load_epi64(c->data + 0 * VECLEN);
            a0 = _mm512_adc_epi52(a0, 0, b0, &scarry);
            _mm512_store_epi64(c->data + 0 * VECLEN, a0);

            i = 1;
            while (scarry > 0)
            {
                a1 = _mm512_load_epi64(c->data + i * VECLEN);
                a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
                _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
                i++;
            }
        }

        // finally clear any hi-bits in the result
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif

    }
    else if (mdata->isMersenne > 0)
    {
        // now add the low part into the high.
        scarry = 0;
        for (i = 0; i < wshift; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            b0 = _mm512_load_epi64(s->data + i * VECLEN);
            a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

        //c00 = _mm512_adc_epi52(c00, scarry, c08, &scarry);
        //c01 = _mm512_adc_epi52(c01, scarry, c09, &scarry);
        //c02 = _mm512_adc_epi52(c02, scarry, c10, &scarry);
        //c03 = _mm512_adc_epi52(c03, scarry, c11, &scarry);
        //c04 = _mm512_adc_epi52(c04, scarry, c12, &scarry);
        //c05 = _mm512_adc_epi52(c05, scarry, c13, &scarry);
        //c06 = _mm512_adc_epi52(c06, scarry, c14, &scarry);
        //c07 = _mm512_adc_epi52(c07, scarry, c15, &scarry);

        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(s->data + i * VECLEN);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add it back in.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);

        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        //scarry = _mm512_test_epi64_mask(c07, vbpshift);
        //c00 = _mm512_addcarry_epi52(c00, scarry, &scarry); if (scarry == 0) goto carryadddone;
        //c01 = _mm512_addcarry_epi52(c01, scarry, &scarry); if (scarry == 0) goto carryadddone;
        //c02 = _mm512_addcarry_epi52(c02, scarry, &scarry); if (scarry == 0) goto carryadddone;
        //c03 = _mm512_addcarry_epi52(c03, scarry, &scarry); if (scarry == 0) goto carryadddone;
        //c04 = _mm512_addcarry_epi52(c04, scarry, &scarry); if (scarry == 0) goto carryadddone;
        //c05 = _mm512_addcarry_epi52(c05, scarry, &scarry); if (scarry == 0) goto carryadddone;
        //c06 = _mm512_addcarry_epi52(c06, scarry, &scarry); if (scarry == 0) goto carryadddone;
        //c07 = _mm512_addcarry_epi52(c07, scarry, &scarry); if (scarry == 0) goto carryadddone;

    carryadddone:


        // clear the potential hi-bit
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

        //c07 = _mm512_and_epi64(vbshift, c07);

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_set1_epi64(0));
        }

        //_mm512_store_epi64(c->data + 0 * VECLEN, c00);
        //_mm512_store_epi64(c->data + 1 * VECLEN, c01);
        //_mm512_store_epi64(c->data + 2 * VECLEN, c02);
        //_mm512_store_epi64(c->data + 3 * VECLEN, c03);
        //_mm512_store_epi64(c->data + 4 * VECLEN, c04);
        //_mm512_store_epi64(c->data + 5 * VECLEN, c05);
        //_mm512_store_epi64(c->data + 6 * VECLEN, c06);
        //_mm512_store_epi64(c->data + 7 * VECLEN, c07);

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif
    }
    else
    {
        // now subtract the hi part from the lo
        scarry = 0;
        for (i = 0; i < wshift; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);   // hi
            b0 = _mm512_load_epi64(s->data + i * VECLEN);   // lo
            a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(s->data + i * VECLEN);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add 1.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // and add one to the hi-bit.  This should resolve the borrow bit.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        //a1 = _mm512_add_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift))), a1);
        //_mm512_store_epi64(c->data + wshift * VECLEN,
        //    _mm512_and_epi64(_mm512_set1_epi64(0xfffffffffffffULL), a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif

        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));
    }

    c->size = NWORDS;
    return;
}

void vecsqrmod52_mersenne_416(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    // 8x sqr:
    // input 8 bignums in the even lanes of a.
    // output 8 squaremod bignums in the even lanes of c.
    int i, j, k;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t NBLOCKS = mdata->NBLOCKS;
    vec_bignum_t* b = a;

    // needed in loops
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19
    __m512i c00, c01, c02, c03, c04, c05, c06, c07,
        c08, c09, c10, c11, c12, c13, c14, c15;

#ifndef IFMA
    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
#endif

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry;

    // deal with the sign
    c->size = NWORDS;
    c->signmask = 0;

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

#ifdef DEBUG_MERSENNE
    printf("commencing vecsqrmod52_mersenne\n");
    print_vechexbignum(a, "input a: ");
    print_vechexbignum(b, "input b: ");
    print_vechexbignum(s, "input s: ");
#endif

    // first half sqr
    i = 0;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (k = 0, j = i; j > (i + 1) / 2; j--, k++)
        {
            // when i = 0, j = 0, j > 0 --> 0 iterations
            // when i = 1, j = 1, j > 1 --> 0 iterations
            // when i = 2, j = 2, j > 1 --> 1 iteration @ a[3..0], b[5..11]
            // when i = 3, j = 3, j > 2 --> 1 iteration @ a[3..0], b[9..15]
            // when i = 4, j = 4, j > 2 --> 2 iteration @ a[3..0], b[13..19] and a[7..4], b[9..15]
            // when i = 5, j = 5, j > 3 --> 2 iteration @ a[3..0], b[17..23] and a[7..4], b[13..19]
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }


        {
            // i even
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);

#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a0, a1, te2, te3);
            VEC_MUL_ACCUM_LOHI_PD(a0, a2, te4, te5);
            VEC_MUL_ACCUM_LOHI_PD(a0, a3, te6, te7);
            VEC_MUL_ACCUM_LOHI_PD(a1, a2, te6, te7);
#else
            //prod1_e = _mm512_mul_epu32(a0, a1);    // te2
            //prod2_e = _mm512_mul_epu32(a0, a2);    // te4
            //prod3_e = _mm512_mul_epu32(a0, a3);    // te6
            {
                prod1_ld = _mm512_cvtepu64_pd(a0);
                prod2_ld = _mm512_cvtepu64_pd(a1);
                prod3_ld = _mm512_cvtepu64_pd(a2);
                prod4_ld = _mm512_cvtepu64_pd(a3);

                prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a1 -> to te2/3
                prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a2 -> to te4/5
                prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a3 -> to te6/7

                te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod2_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod3_hd));
                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod4_hd));

                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
                prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

                prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod2_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod3_ld));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod4_ld));
            }

            //prod1_e = _mm512_mul_epu32(a1, a2);    // te6
            {
                prod1_ld = _mm512_cvtepu64_pd(a1);
                prod2_ld = _mm512_cvtepu64_pd(a2);

                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * a2 -> to te6/7

                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            }

            // all terms so far need to be doubled.  
            // but to do that we first need to remove all bias
            // that has been accumulated so far.
            SUB_BIAS_HI(
                k * 4 + 0,
                k * 4 + 1,
                k * 4 + 1,
                k * 4 + 2);
            SUB_BIAS_LO(
                k * 4 + 0,
                k * 4 + 1,
                k * 4 + 1,
                k * 4 + 2);
#endif

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a0, a0, te0, te1);
            VEC_MUL_ACCUM_LOHI_PD(a1, a1, te4, te5);
#else
            // finally, accumulate the two non-doubled terms.
            //prod1_e = _mm512_mul_epu32(a0, a0);
            //prod2_e = _mm512_mul_epu32(a1, a1);
            {
                prod1_ld = _mm512_cvtepu64_pd(a0);
                prod2_ld = _mm512_cvtepu64_pd(a1);

                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a0 -> to te0/1
                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * a1 -> to te4/5

                te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
            }
#endif

        }

#ifndef IFMA
        // need to remove bias from the two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            1,
            0,
            1,
            0);
        SUB_BIAS_LO(
            1,
            0,
            1,
            0);
#endif

        // final column accumulation
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        // i = 3, a3..0 = c[0..3]
        // i = 4, a3..0 = c[4..7]
        // i = 5, a3..0 = c[8..11]
        c00 = a3;
        c01 = a2;
        c02 = a1;
        c03 = a0;
    }

    i = 1;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (k = 0, j = i; j > (i + 1) / 2; j--, k++)
        {
            // when i = 0, j = 0, j > 0 --> 0 iterations
            // when i = 1, j = 1, j > 1 --> 0 iterations
            // when i = 2, j = 2, j > 1 --> 1 iteration @ a[3..0], b[5..11]
            // when i = 3, j = 3, j > 2 --> 1 iteration @ a[3..0], b[9..15]
            // when i = 4, j = 4, j > 2 --> 2 iteration @ a[3..0], b[13..19] and a[7..4], b[9..15]
            // when i = 5, j = 5, j > 3 --> 2 iteration @ a[3..0], b[17..23] and a[7..4], b[13..19]
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }


        {
            // i odd
            // for 384-bit inputs when i == 1, j = 0, a = {3,2,1,0} and b = {2,3,4,5,6,7}
            // for 512-bit inputs when i == 3, j = 1, a = {7,6,5,4} and b = {6,7,8,9,a,b}
            // for 512-bit inputs when i == 1, j = 0, a = {3,2,1,0} and b = {2,3,4,5,6,7}

            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);


#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a2, b2, te0, te1);
            VEC_MUL_ACCUM_LOHI_PD(a1, b2, te2, te3);
            VEC_MUL_ACCUM_LOHI_PD(a1, b3, te4, te5);
            VEC_MUL_ACCUM_LOHI_PD(a0, b3, te6, te7);
#else
            //prod1_e = _mm512_mul_epu32(a2, b2);   // te0
            //prod2_e = _mm512_mul_epu32(a1, b2);   // te2
            //prod3_e = _mm512_mul_epu32(a1, b3);   // te4
            //prod4_e = _mm512_mul_epu32(a0, b3);   // te6
            //ACCUM_4X_DOUBLED_PROD;
            {
                prod5_ld = _mm512_cvtepu64_pd(a0);
                prod1_ld = _mm512_cvtepu64_pd(a1);
                prod2_ld = _mm512_cvtepu64_pd(a2);
                prod3_ld = _mm512_cvtepu64_pd(b2);
                prod4_ld = _mm512_cvtepu64_pd(b3);

                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b2 -> to te0/1 
                prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b2 -> to te2/3
                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b3 -> to te4/5
                prod4_hd = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b3 -> to te6/7

                te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod1_hd));
                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod4_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
                prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod4_ld = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod1_ld));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod4_ld));
            }
#endif


#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a3, b3, te0, te1);
            VEC_MUL_ACCUM_LOHI_PD(a2, b3, te2, te3);
            VEC_MUL_ACCUM_LOHI_PD(a2, b4, te4, te5);
            VEC_MUL_ACCUM_LOHI_PD(a1, b4, te6, te7);
#else
            //prod1_e = _mm512_mul_epu32(a3, b3);   // te0
            //prod2_e = _mm512_mul_epu32(a2, b3);   // te2
            //prod3_e = _mm512_mul_epu32(a2, b4);   // te4
            //prod4_e = _mm512_mul_epu32(a1, b4);   // te6
            //ACCUM_4X_DOUBLED_PROD;
            {
                prod5_ld = _mm512_cvtepu64_pd(a1);
                prod1_ld = _mm512_cvtepu64_pd(a2);
                prod2_ld = _mm512_cvtepu64_pd(a3);
                prod3_ld = _mm512_cvtepu64_pd(b3);
                prod4_ld = _mm512_cvtepu64_pd(b4);

                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * b3 -> to te0/1
                prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b3 -> to te2/3
                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b4 -> to te4/5
                prod4_hd = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b4 -> to te6/7

                te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod1_hd));
                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod4_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
                prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod4_ld = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod1_ld));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod4_ld));
            }
#endif


#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a3, b4, te2, te3);
            VEC_MUL_ACCUM_LOHI_PD(a3, b5, te4, te5);
            VEC_MUL_ACCUM_LOHI_PD(a2, b5, te6, te7);
            VEC_MUL_ACCUM_LOHI_PD(a3, b6, te6, te7);
#else
            //prod1_e = _mm512_mul_epu32(a3, b4);  // te2
            //prod2_e = _mm512_mul_epu32(a3, b5);  // te4
            //prod3_e = _mm512_mul_epu32(a2, b5);  // te6
            {
                prod1_ld = _mm512_cvtepu64_pd(a2);
                prod2_ld = _mm512_cvtepu64_pd(a3);
                prod3_ld = _mm512_cvtepu64_pd(b4);
                prod4_ld = _mm512_cvtepu64_pd(b5);

                prod3_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * b4 -> to te2/3
                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * b5 -> to te4/5
                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b5 -> to te6/7

                te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));
                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);

                prod3_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            }

            //prod1_e = _mm512_mul_epu32(a3, b6);  // te6
            {
                prod1_ld = _mm512_cvtepu64_pd(a3);
                prod2_ld = _mm512_cvtepu64_pd(b6);

                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * b6 -> to te6/7

                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            }

            // all terms so far need to be doubled.  
            // but to do that we first need to remove all bias
            // that has been accumulated so far.
            SUB_BIAS_HI(
                k * 4 + 2,
                k * 4 + 3,
                k * 4 + 3,
                k * 4 + 4);
            SUB_BIAS_LO(
                k * 4 + 2,
                k * 4 + 3,
                k * 4 + 3,
                k * 4 + 4);
#endif

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a1, a1, te0, te1);
            VEC_MUL_ACCUM_LOHI_PD(a0, a0, te4, te5);
#else
            // finally, accumulate the two non-doubled terms.
            //prod1_e = _mm512_mul_epu32(a1, a1);    // te0
            //prod2_e = _mm512_mul_epu32(a0, a0);    // te4
            {
                prod1_ld = _mm512_cvtepu64_pd(a0);
                prod2_ld = _mm512_cvtepu64_pd(a1);

                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * a1 -> to te0/1
                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a0 -> to te4/5

                te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod1_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod1_ld));
            }
#endif

        }


#ifndef IFMA
        // need to remove bias from the two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            1,
            0,
            1,
            0);
        SUB_BIAS_LO(
            1,
            0,
            1,
            0);
#endif

        // final column accumulation
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        // i = 3, a3..0 = c[0..3]
        // i = 4, a3..0 = c[4..7]
        // i = 5, a3..0 = c[8..11]
        c04 = a3;
        c05 = a2;
        c06 = a1;
        c07 = a0;
    }
#ifdef DEBUG_MERSENNE
    print_vechexbignum(s, "after lo half:");
#endif

    // second half sqr
    i = 0;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (k = 0, j = 0; j < (NBLOCKS - i - 1) / 2; j++, k++)
        {
            // Compute a solid block (all matching terms are in the lower
            // half triangle of the expansion).
            a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
            a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
            a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);
            a3 = _mm512_load_epi64(a->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);

            b0 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // The final block shape depends on the parity of i and NBLOCKS
        // Block shape 1 (small upper triangle) if 'i' is odd and NBLOCKS is even,
        // or if 'i' is even and NBLOCKS is odd.
        // Block shape 2 if 'i' is odd and NBLOCKS is odd or if 'i' is even
        // and NBLOCKS is even.
        // NBLOCKS is even
        {
            // i even, block shape 1.
            {

                // always a continuation of the full-block loop, so use the same 
                // loading pattern.  Only now we don't need as many b-terms.
                a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);		// {f, b}
                a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);		// {e, a}
                a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);		// {d, 9}
                a3 = _mm512_load_epi64(a->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);		// {c, 8}

                b0 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 1) * VECLEN); // {9, 5}
                b1 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 2) * VECLEN);	// {a, 6}
                b2 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 3) * VECLEN);	// {b, 7}
                b3 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 4) * VECLEN);	// {c, 8}
                b4 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 5) * VECLEN); // {d, 9}

                //prod1_e = _mm512_mul_epu32(a0, b0);    // te0
                //prod2_e = _mm512_mul_epu32(a0, b1);    // te2
                //prod3_e = _mm512_mul_epu32(a0, b2);    // te4
                //prod4_e = _mm512_mul_epu32(a0, b3);    // te6
                //ACCUM_4X_DOUBLED_PROD;
                VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);

                //prod1_e = _mm512_mul_epu32(a1, b1);    // te0
                //prod2_e = _mm512_mul_epu32(a1, b2);    // te2
                //prod3_e = _mm512_mul_epu32(a1, b3);    // te4
                //prod4_e = _mm512_mul_epu32(a1, b4);    // te6
                //ACCUM_4X_DOUBLED_PROD;
                VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a2, b2, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a2, b3, te2, te3);
#else
                //prod1_e = _mm512_mul_epu32(a2, b2);    // te0
                //prod2_e = _mm512_mul_epu32(a2, b3);    // te2
                {
                    prod1_ld = _mm512_cvtepu64_pd(a2);
                    prod2_ld = _mm512_cvtepu64_pd(b2);
                    prod3_ld = _mm512_cvtepu64_pd(b3);

                    prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b2 -> to te0/1
                    prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b3 -> to te2/3

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));

                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);

                    prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                }

                // all terms so far need to be doubled.  
                // but to do that we first need to remove all bias
                // that has been accumulated so far.
                SUB_BIAS_HI(
                    j * 4 + 3,
                    j * 4 + 3,
                    j * 4 + 2,
                    j * 4 + 2);
                SUB_BIAS_LO(
                    j * 4 + 3,
                    j * 4 + 3,
                    j * 4 + 2,
                    j * 4 + 2);
#endif

                // now double
                te1 = _mm512_slli_epi64(te1, 1);
                te3 = _mm512_slli_epi64(te3, 1);
                te5 = _mm512_slli_epi64(te5, 1);
                te7 = _mm512_slli_epi64(te7, 1);
                te0 = _mm512_slli_epi64(te0, 1);
                te2 = _mm512_slli_epi64(te2, 1);
                te4 = _mm512_slli_epi64(te4, 1);
                te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a3, a3, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a2, a2, te4, te5);
#else
                // finally, accumulate the two non-doubled terms.
                //prod1_e = _mm512_mul_epu32(a3, a3);   // te0
                //prod2_e = _mm512_mul_epu32(a2, a2);   // te4
                {
                    prod1_ld = _mm512_cvtepu64_pd(a3);
                    prod2_ld = _mm512_cvtepu64_pd(a2);

                    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * a3 -> to te0/1
                    prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * a2 -> to te4/5

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));
                    te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));

                    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
                    te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
                }
#endif
            }
        }

#ifndef IFMA
        // need to remove bias from the two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            1,
            0,
            1,
            0);
        SUB_BIAS_LO(
            1,
            0,
            1,
            0);
#endif

        // final column accumulation
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        // i = 3, a3..0 = c[0..3]
        // i = 4, a3..0 = c[4..7]
        // i = 5, a3..0 = c[8..11]
        c08 = a3;
        c09 = a2;
        c10 = a1;
        c11 = a0;


    }

    i = 1;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (k = 0, j = 0; j < (NBLOCKS - i - 1) / 2; j++, k++)
        {
            // Compute a solid block (all matching terms are in the lower
            // half triangle of the expansion).
            a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
            a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
            a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);
            a3 = _mm512_load_epi64(a->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);

            b0 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // The final block shape depends on the parity of i and NBLOCKS
        // Block shape 1 (small upper triangle) if 'i' is odd and NBLOCKS is even,
        // or if 'i' is even and NBLOCKS is odd.
        // Block shape 2 if 'i' is odd and NBLOCKS is odd or if 'i' is even
        // and NBLOCKS is even.
        // NBLOCKS is even
        {
            // i odd, block shape 1.
            {
                // always a continuation of the full-block loop, so use the same 
                // loading pattern.  Only now we don't need as many b-terms.
                a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);  // {f, b}
                a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);  // {e, a}
                a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);  // {d, 9}

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a0, a2, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a0, a1, te2, te3);
#else
                //k == 0;
                //prod1_e = _mm512_mul_epu32(a0, a2);   // te0
                //prod2_e = _mm512_mul_epu32(a0, a1);   // te2
                {
                    prod1_ld = _mm512_cvtepu64_pd(a0);
                    prod2_ld = _mm512_cvtepu64_pd(a1);
                    prod3_ld = _mm512_cvtepu64_pd(a2);

                    prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a2 -> to te0/1
                    prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a1 -> to te2/3

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod3_hd));
                    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod2_hd));

                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);

                    prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod3_ld));
                    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod2_ld));
                }

                // all terms so far need to be doubled.  
                // but to do that we first need to remove all bias
                // that has been accumulated so far.
                SUB_BIAS_HI(
                    j * 4 + 1,
                    j * 4 + 1,
                    j * 4 + 0,
                    j * 4 + 0);
                SUB_BIAS_LO(
                    j * 4 + 1,
                    j * 4 + 1,
                    j * 4 + 0,
                    j * 4 + 0);
#endif

                // now double
                te1 = _mm512_slli_epi64(te1, 1);
                te3 = _mm512_slli_epi64(te3, 1);
                te5 = _mm512_slli_epi64(te5, 1);
                te7 = _mm512_slli_epi64(te7, 1);
                te0 = _mm512_slli_epi64(te0, 1);
                te2 = _mm512_slli_epi64(te2, 1);
                te4 = _mm512_slli_epi64(te4, 1);
                te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a1, a1, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a0, a0, te4, te5);
#else
                // finally, accumulate the two non-doubled terms.
                //prod1_e = _mm512_mul_epu32(a1, a1);       // te0
                //prod2_e = _mm512_mul_epu32(a0, a0);       // te4
                {
                    prod1_ld = _mm512_cvtepu64_pd(a1);
                    prod2_ld = _mm512_cvtepu64_pd(a0);

                    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b1 -> to te0/1
                    prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b0 -> to te4/5

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));
                    te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));

                    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
                    te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
                }
#endif
            }

        }

#ifndef IFMA
        // need to remove bias from the two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            1,
            0,
            1,
            0);
        SUB_BIAS_LO(
            1,
            0,
            1,
            0);
#endif

        // final column accumulation
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        // i = 3, a3..0 = c[0..3]
        // i = 4, a3..0 = c[4..7]
        // i = 5, a3..0 = c[8..11]
        c12 = a3;
        c13 = a2;
        c14 = a1;
        c15 = a0;


    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(s, "after hi half:");
#endif

    // reduce by adding hi to lo.  first right shift hi into output.
    __m512i vbshift, vbpshift;
    int bshift = mdata->nbits % 52;
    int wshift = mdata->nbits / 52;

    if (mdata->use_vnbits)
    {
        int bshift_a[VECLEN];
        for (i = 0; i < VECLEN; i++)
        {
            bshift_a[i] = mdata->vnbits[i] % 52;
        }

        //vbshift = _mm512_set_epi64(
        //    (1ULL << (uint64_t)(bshift_a[7])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[0])) - 1ULL);
        //
        //vbpshift = _mm512_set_epi64(
        //    (1ULL << (uint64_t)(bshift_a[7])),
        //    (1ULL << (uint64_t)(bshift_a[6])),
        //    (1ULL << (uint64_t)(bshift_a[5])),
        //    (1ULL << (uint64_t)(bshift_a[4])),
        //    (1ULL << (uint64_t)(bshift_a[3])),
        //    (1ULL << (uint64_t)(bshift_a[2])),
        //    (1ULL << (uint64_t)(bshift_a[1])),
        //    (1ULL << (uint64_t)(bshift_a[0])));

        vbshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[7])) - 1ULL);

        vbpshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])),
            (1ULL << (uint64_t)(bshift_a[1])),
            (1ULL << (uint64_t)(bshift_a[2])),
            (1ULL << (uint64_t)(bshift_a[3])),
            (1ULL << (uint64_t)(bshift_a[4])),
            (1ULL << (uint64_t)(bshift_a[5])),
            (1ULL << (uint64_t)(bshift_a[6])),
            (1ULL << (uint64_t)(bshift_a[7])));

        vec_bignum52_mask_rshift_vn(s, c, mdata->vnbits, 0xff);
    }
    else
    {
        _mm512_store_epi64(s->data + 0 * VECLEN, c00);
        _mm512_store_epi64(s->data + 1 * VECLEN, c01);
        _mm512_store_epi64(s->data + 2 * VECLEN, c02);
        _mm512_store_epi64(s->data + 3 * VECLEN, c03);
        _mm512_store_epi64(s->data + 4 * VECLEN, c04);
        _mm512_store_epi64(s->data + 5 * VECLEN, c05);
        _mm512_store_epi64(s->data + 6 * VECLEN, c06);
        _mm512_store_epi64(s->data + 7 * VECLEN, c07);
        _mm512_store_epi64(s->data + 8 * VECLEN, c08);
        _mm512_store_epi64(s->data + 9 * VECLEN, c09);
        _mm512_store_epi64(s->data + 10 * VECLEN, c10);
        _mm512_store_epi64(s->data + 11 * VECLEN, c11);
        _mm512_store_epi64(s->data + 12 * VECLEN, c12);
        _mm512_store_epi64(s->data + 13 * VECLEN, c13);
        _mm512_store_epi64(s->data + 14 * VECLEN, c14);
        _mm512_store_epi64(s->data + 15 * VECLEN, c15);

        vbshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL);
        vbpshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)));
        vec_bignum52_mask_rshift_n(s, c, mdata->nbits, 0xff);
    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(c, "hi part:");
    print_vechexbignum(s, "lo part:");
#endif

    if (mdata->isMersenne > 1)
    {
        // this number has been identified as a "pseudo"-Mersenne: 2^n-c, with
        // 'c' small.  Fast reduction is still possible.

        // multiply c * hi
        b0 = _mm512_set1_epi64(mdata->isMersenne);
        vecmul52_1(c, b0, c, n, s, mdata);

        // clear any hi-bits in the lo result
        a1 = _mm512_load_epi64(s->data + wshift * VECLEN);
        _mm512_store_epi64(s->data + wshift * VECLEN,
            _mm512_and_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL), a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(s, "cleared low part:");
#endif

        // add c * hi + lo
        scarry = 0;
        for (i = 0; i <= NWORDS; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            b0 = _mm512_load_epi64(s->data + i * VECLEN);
            a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "lo + hi * c:");
#endif

        // now we need to do it again, but the hi part is just a single word so
        // we only have a single-precision vecmul.  It's possible these hi bits
        // are located in two words
        //vec_bignum52_mask_rshift_n(c, s, mdata->nbits, 0xff);
        //a0 = _mm512_load_epi64(s->data);

        a0 = _mm512_load_epi64(c->data + wshift * VECLEN);
        a0 = _mm512_srli_epi64(a0, bshift);
        a1 = _mm512_load_epi64(c->data + (wshift + 1) * VECLEN);
        a0 = _mm512_or_epi64(a0, _mm512_slli_epi64(a1, 52 - bshift));

        b0 = _mm512_set1_epi64(mdata->isMersenne);
        VEC_MUL_LOHI_PD(a0, b0, acc_e0, acc_e1);

        // clear any hi-bits now that we have the multiplier
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL), a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_set1_epi64(0));
        }

        // now add that into the lo part again.  Here we add in 
        // the two words of the hi-mul result.
        b0 = _mm512_load_epi64(c->data + 0 * VECLEN);
        b1 = _mm512_load_epi64(c->data + 1 * VECLEN);
        b2 = zero;
        b3 = zero;
        acc_e0 = _mm512_add_epi64(acc_e0, b0);
        VEC_CARRYPROP_LOHI(acc_e0, b2);         // lo word and carry for 1st add.
        acc_e1 = _mm512_add_epi64(acc_e1, b1);
        VEC_CARRYPROP_LOHI(acc_e1, b3);         // lo word and carry for 2nd add.
        acc_e1 = _mm512_add_epi64(acc_e1, b2);  // add 1st carry to 2nd lo word.
        VEC_CARRYPROP_LOHI(acc_e1, b3);         // that could have a carry, add it to b3.

        // save first two result words 
        _mm512_store_epi64(c->data + 0 * VECLEN, acc_e0);
        _mm512_store_epi64(c->data + 1 * VECLEN, acc_e1);

        // propagate the carry through the rest of the lo result.
        b2 = _mm512_load_epi64(c->data + 2 * VECLEN);
        scarry = 0;
        a0 = _mm512_adc_epi52(b2, scarry, b3, &scarry);
        _mm512_store_epi64(c->data + 2 * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i = 3; (i < wshift) && (scarry > 0); i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "lo + hi * c:");
#endif

        if (scarry > 0)
        {
            // it is unlikely, but possible that we still have a final carry to propagate
            b0 = _mm512_set1_epi64(mdata->isMersenne);
            a0 = _mm512_load_epi64(c->data + 0 * VECLEN);
            a0 = _mm512_adc_epi52(a0, 0, b0, &scarry);
            _mm512_store_epi64(c->data + 0 * VECLEN, a0);

            i = 1;
            while (scarry > 0)
            {
                a1 = _mm512_load_epi64(c->data + i * VECLEN);
                a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
                _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
                i++;
            }
        }

        // finally clear any hi-bits in the result
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL), a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
#endif

    }
    else if (mdata->isMersenne > 0)
    {
        // now add the low part into the high.
        scarry = 0;
        for (i = 0; i < wshift; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            b0 = _mm512_load_epi64(s->data + i * VECLEN);
            a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(s->data + i * VECLEN);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add it back in.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // clear the potential hi-bit
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif
    }
    else
    {
        // now subtract the hi part from the lo
        scarry = 0;
        for (i = 0; i < wshift; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);   // hi
            b0 = _mm512_load_epi64(s->data + i * VECLEN);   // lo
            a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(s->data + i * VECLEN);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));

        // zero any remaining lo words.
        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after sub:");
#endif

        // if there was a carry, add 1 to the lo bit.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // and add one to the hi-bit.  This should resolve the borrow bit.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        //a1 = _mm512_add_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift))), a1);
        //_mm512_store_epi64(c->data + wshift * VECLEN,
        //    _mm512_and_epi64(_mm512_set1_epi64(0xfffffffffffffULL), a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
#endif

        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after bit clear:");
        exit(1);
#endif

    }

    c->size = NWORDS;
    return;
}

void vecmulmod52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    int i, j;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t NBLOCKS = mdata->NBLOCKS;

    //if (NWORDS == 8)
    //{
    //    vecmulmod52_mersenne_416(a, b, c, n, s, mdata);
    //    return;
    //}


    // needed in loops
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19

#ifndef IFMA
    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
#endif

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry;

    // deal with the sign
    c->size = NWORDS;
    c->signmask = a->signmask ^ b->signmask;

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

#ifdef DEBUG_MERSENNE
    printf("commencing vecmulmod52_mersenne\n");
    print_vechexbignum(a, "input a: ");
    print_vechexbignum(b, "input b: ");
    print_vechexbignum(s, "input s: ");
#endif

    // first half mul
    for (i = 0; i < NBLOCKS; i++)
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (j = i; j > 0; j--)
        {
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum
        a0 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 0) * VECLEN);
        a1 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + 3 * VECLEN);
        b1 = _mm512_load_epi64(b->data + 2 * VECLEN);
        b2 = _mm512_load_epi64(b->data + 1 * VECLEN);
        b3 = _mm512_load_epi64(b->data + 0 * VECLEN);

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a0, b3, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a1, b3, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a0, b2, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a2, b3, te4, te5);
#else
        // ======
        //VEC_MUL_ACCUM_LOHI_PD(a0, b3, te0, te1);
        //VEC_MUL_ACCUM_LOHI_PD(a1, b3, te2, te3);
        //VEC_MUL_ACCUM_LOHI_PD(a0, b2, te2, te3);
        //VEC_MUL_ACCUM_LOHI_PD(a2, b3, te4, te5);
        {
            prod5_ld = _mm512_cvtepu64_pd(a0);
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(a2);
            prod3_ld = _mm512_cvtepu64_pd(b2);
            prod4_ld = _mm512_cvtepu64_pd(b3);

            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod1_hd));  // a1 * b3 -> to te2/3
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));  // a2 * b3 -> to te4/5 
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a0 * b2 -> to te2/3
            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod4_hd));  // a0 * b3 -> to te0/1

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod1_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod4_ld));
        }
#endif

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a1, b2, te4, te5);
        VEC_MUL_ACCUM_LOHI_PD(a0, b1, te4, te5);
        VEC_MUL_ACCUM_LOHI_PD(a1, b1, te6, te7);
        VEC_MUL_ACCUM_LOHI_PD(a0, b0, te6, te7);
#else
        // ======
        //VEC_MUL_ACCUM_LOHI_PD(a1, b2, te4, te5);
        //VEC_MUL_ACCUM_LOHI_PD(a0, b1, te4, te5);
        //VEC_MUL_ACCUM_LOHI_PD(a1, b1, te6, te7);
        //VEC_MUL_ACCUM_LOHI_PD(a0, b0, te6, te7);
        {
            prod5_ld = _mm512_cvtepu64_pd(a0);
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(b0);
            prod3_ld = _mm512_cvtepu64_pd(b1);
            prod4_ld = _mm512_cvtepu64_pd(b2);

            prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_hd = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod2_hd));
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod3_hd));
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod4_hd));

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_ld = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod2_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod3_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod4_ld));
        }
#endif

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a3, b3, te6, te7);
        VEC_MUL_ACCUM_LOHI_PD(a2, b2, te6, te7);
#else
        //VEC_MUL_ACCUM_LOHI_PD(a3, b3, te6, te7);
        //VEC_MUL_ACCUM_LOHI_PD(a2, b2, te6, te7);
        {
            prod1_ld = _mm512_cvtepu64_pd(a2);
            prod2_ld = _mm512_cvtepu64_pd(a3);
            prod3_ld = _mm512_cvtepu64_pd(b2);
            prod4_ld = _mm512_cvtepu64_pd(b3);

            prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod2_hd));

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

            prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod2_ld));
        }

        // subtract out all of the bias at once.  these are
        // counts of how many times we have mul_accumulated into each column
        // during the block a*b and final triangular a*b.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 2,
            i * 4 + 3,
            i * 4 + 4);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 2,
            i * 4 + 3,
            i * 4 + 4);
#endif



        // now, do a carry-propagating column summation and store the results.
        {
            j = 0;
            // accumulate this column-sum:
            // carry propagate low to high.
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 1;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te2);
            acc_e1 = _mm512_add_epi64(acc_e1, te3);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 2;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te4);
            acc_e1 = _mm512_add_epi64(acc_e1, te5);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 3;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te6);
            acc_e1 = _mm512_add_epi64(acc_e1, te7);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;
        }
    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(s, "after lo half:");
#endif

    // second half mul
    for (i = NBLOCKS; i < 2 * NBLOCKS; i++)
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (j = i - NBLOCKS + 1; j < NBLOCKS; j++)
        {
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }


        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN);

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a1, b0, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b1, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b0, te2, te3);
#else
        // ======
        //VEC_MUL_ACCUM_LOHI_PD(a1, b0, te0, te1);
        //VEC_MUL_ACCUM_LOHI_PD(a2, b1, te0, te1);
        //VEC_MUL_ACCUM_LOHI_PD(a2, b0, te2, te3);
        {
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(a2);
            prod3_ld = _mm512_cvtepu64_pd(b0);
            prod4_ld = _mm512_cvtepu64_pd(b1);

            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));  // a1*b0 -> t0/1
            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod4_hd));  // a2*b1 -> t0/1
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a2*b0 -> t2/3

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod4_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
        }
#endif

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a3, b2, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a3, b1, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a3, b0, te4, te5);
#else
        //VEC_MUL_ACCUM_LOHI_PD(a3, b2, te0, te1);
        //VEC_MUL_ACCUM_LOHI_PD(a3, b1, te2, te3);
        //VEC_MUL_ACCUM_LOHI_PD(a3, b0, te4, te5);
        {
            prod1_ld = _mm512_cvtepu64_pd(a3);
            prod2_ld = _mm512_cvtepu64_pd(b2);
            prod3_ld = _mm512_cvtepu64_pd(b1);
            prod4_ld = _mm512_cvtepu64_pd(b0);

            prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));  // a3*b2 -> t0/1
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a3*b1 -> t2/3
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod4_hd));  // a3*b0 -> t4/5

            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod4_ld));
        }


        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            (2 * NBLOCKS - i - 1) * 4 + 3,
            (2 * NBLOCKS - i - 1) * 4 + 2,
            (2 * NBLOCKS - i - 1) * 4 + 1,
            (2 * NBLOCKS - i - 1) * 4 + 0);
        SUB_BIAS_LO(
            (2 * NBLOCKS - i - 1) * 4 + 3,
            (2 * NBLOCKS - i - 1) * 4 + 2,
            (2 * NBLOCKS - i - 1) * 4 + 1,
            (2 * NBLOCKS - i - 1) * 4 + 0);
#endif

        // final column accumulation
        {
            j = 0;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a3 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN, 
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 1;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te2);
            acc_e1 = _mm512_add_epi64(acc_e1, te3);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a2 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 2;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te4);
            acc_e1 = _mm512_add_epi64(acc_e1, te5);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a1 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 3;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te6);
            acc_e1 = _mm512_add_epi64(acc_e1, te7);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a0 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            a3 = _mm512_and_epi64(vlmask, a3);
            a2 = _mm512_and_epi64(vlmask, a2);
            a1 = _mm512_and_epi64(vlmask, a1);
            a0 = _mm512_and_epi64(vlmask, a0);

            _mm512_store_epi64(s->data + (i * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + 3) * VECLEN, a0);

        }
    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(s, "after hi half:");
#endif

    // reduce by adding hi to lo.  first right shift hi into output.
    __m512i vbshift, vbpshift;
    int bshift = mdata->nbits % 52;
    int wshift = mdata->nbits / 52;
#ifdef DEBUG_MERSENNE
    mpz_t ga, gb;
    mpz_init(ga);
    mpz_init(gb);
    extract_bignum_from_vec_to_mpz(ga, a, 0, NWORDS);
    extract_bignum_from_vec_to_mpz(gb, b, 0, NWORDS);
    mpz_mul(ga, ga, gb);
#endif

    if (mdata->use_vnbits)
    {
        int bshift_a[VECLEN];
        for (i = 0; i < VECLEN; i++)
        {
            bshift_a[i] = mdata->vnbits[i] % 52;
        }

        //vbshift = _mm512_set_epi64(
        //    (1ULL << (uint64_t)(bshift_a[7])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[0])) - 1ULL);
        //
        //vbpshift = _mm512_set_epi64(
        //    (1ULL << (uint64_t)(bshift_a[7])),
        //    (1ULL << (uint64_t)(bshift_a[6])),
        //    (1ULL << (uint64_t)(bshift_a[5])),
        //    (1ULL << (uint64_t)(bshift_a[4])),
        //    (1ULL << (uint64_t)(bshift_a[3])),
        //    (1ULL << (uint64_t)(bshift_a[2])),
        //    (1ULL << (uint64_t)(bshift_a[1])),
        //    (1ULL << (uint64_t)(bshift_a[0])));

        vbshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[7])) - 1ULL);

        vbpshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])),
            (1ULL << (uint64_t)(bshift_a[1])),
            (1ULL << (uint64_t)(bshift_a[2])),
            (1ULL << (uint64_t)(bshift_a[3])),
            (1ULL << (uint64_t)(bshift_a[4])),
            (1ULL << (uint64_t)(bshift_a[5])),
            (1ULL << (uint64_t)(bshift_a[6])),
            (1ULL << (uint64_t)(bshift_a[7])));

        vec_bignum52_mask_rshift_vn(s, c, mdata->vnbits, 0xff);
    }
    else
    {
        vbshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL);
        vbpshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)));
        vec_bignum52_mask_rshift_n(s, c, mdata->nbits, 0xff);
    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(c, "hi part:");
    print_vechexbignum(s, "lo part:");
#endif

    if (mdata->isMersenne > 1)
    {
        // this number has been identified as a "pseudo"-Mersenne: 2^n-c, with
        // 'c' small.  Fast reduction is still possible.

        // multiply c * hi
        b0 = _mm512_set1_epi64(mdata->isMersenne);
        vecmul52_1(c, b0, c, n, s, mdata);

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after mul_1 with c:");
#endif

        // clear any hi-bits in the lo result
        a1 = _mm512_load_epi64(s->data + wshift * VECLEN);
        _mm512_store_epi64(s->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(s, "cleared low part:");
#endif

        // add c * hi + lo
        scarry = 0;
        for (i = 0; i <= NWORDS; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            b0 = _mm512_load_epi64(s->data + i * VECLEN);
            
            //a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            __m512i t = _mm512_add_epi64(a1, b0);
            t = _mm512_add_epi64(t, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpgt_epu64_mask(t, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "lo + hi * c:");
#endif

        // now we need to do it again, but the hi part is just a single word so
        // we only have a single-precision vecmul
        //vec_bignum52_mask_rshift_n(c, s, mdata->nbits, 0xff);
        //a0 = _mm512_load_epi64(s->data);

        // It's possible these hi bits are located in two words
        a0 = _mm512_load_epi64(c->data + wshift * VECLEN);
        a0 = _mm512_srli_epi64(a0, bshift);
        a1 = _mm512_load_epi64(c->data + (wshift + 1) * VECLEN);
        a0 = _mm512_or_epi64(a0, _mm512_slli_epi64(a1, 52 - bshift));

        //print_regvechex(a0, 0, "hi part:");

        b0 = _mm512_set1_epi64(mdata->isMersenne);

        //print_regvechex(b0, 0, "c:");

        VEC_MUL_LOHI_PD(a0, b0, acc_e0, acc_e1);

        // clear any hi-bits now that we have the multiplier
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_set1_epi64(0));
        }

        // now add that into the lo part again.  Here we add in 
        // the two words of the hi-mul result.
        b0 = _mm512_load_epi64(c->data + 0 * VECLEN);
        b1 = _mm512_load_epi64(c->data + 1 * VECLEN);
        b2 = zero;
        b3 = zero;
        acc_e0 = _mm512_add_epi64(acc_e0, b0);
        VEC_CARRYPROP_LOHI(acc_e0, b2);         // lo word and carry for 1st add.
        acc_e1 = _mm512_add_epi64(acc_e1, b1);
        VEC_CARRYPROP_LOHI(acc_e1, b3);         // lo word and carry for 2nd add.
        acc_e1 = _mm512_add_epi64(acc_e1, b2);  // add 1st carry to 2nd lo word.
        VEC_CARRYPROP_LOHI(acc_e1, b3);         // that could have a carry, add it to b3.

        // save first two result words 
        _mm512_store_epi64(c->data + 0 * VECLEN, acc_e0);
        _mm512_store_epi64(c->data + 1 * VECLEN, acc_e1);

        // propagate the carry through the rest of the lo result.
        b2 = _mm512_load_epi64(c->data + 2 * VECLEN);
        scarry = 0;
        a0 = _mm512_adc_epi52(b2, scarry, b3, &scarry);
        _mm512_store_epi64(c->data + 2 * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i = 3; (i < wshift) && (scarry > 0); i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            
            //a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            __m512i t = _mm512_add_epi64(a1, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpeq_epu64_mask(a1, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "lo + hi * c:");
#endif

        // it is unlikely, but possible that we still have a final carry to propagate
        if (scarry > 0)
        {
            printf("unlikely add\n");
            b0 = _mm512_set1_epi64(mdata->isMersenne);
            a0 = _mm512_load_epi64(c->data + 0 * VECLEN);
            a0 = _mm512_adc_epi52(a0, 0, b0, &scarry);
            _mm512_store_epi64(c->data + 0 * VECLEN, a0);

            i = 1;
            while (scarry > 0)
            {
                a1 = _mm512_load_epi64(c->data + i * VECLEN);
                
                //a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
                __m512i t = _mm512_add_epi64(a1, _mm512_maskz_set1_epi64(scarry, 1));
                scarry = _mm512_cmpeq_epu64_mask(a1, vlmask);
                a0 = _mm512_and_epi64(t, vlmask);

                _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
                i++;
            }
        }

        // finally clear any hi-bits in the result
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif

    }
    else if (mdata->isMersenne > 0)
    {
        // now add the low part into the high.
        scarry = 0;
        for (i = 0; i < wshift; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            b0 = _mm512_load_epi64(s->data + i * VECLEN);
            
            //a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            __m512i t = _mm512_add_epi64(a1, b0);
            t = _mm512_add_epi64(t, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpgt_epu64_mask(t, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(s->data + i * VECLEN);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add it back in.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            
            //a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            __m512i t = _mm512_add_epi64(a1, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpeq_epu64_mask(a1, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // clear the potential hi-bit
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");


        mpz_t gn, gc;
        mpz_init(gn);
        mpz_init(gc);
        mpz_set_ui(gn, 1);
        mpz_mul_2exp(gn, gn, mdata->nbits);
        mpz_sub_ui(gn, gn, 1);
        mpz_mod(gc, ga, gn);

        extract_bignum_from_vec_to_mpz(gb, c, 0, NWORDS);

        if (mpz_cmp(gc, gb) != 0)
        {
            gmp_printf("a = %Zd\n", ga);
            gmp_printf("n = %Zd\n", gn);

            gmp_printf("b = %Zd\n", gb);
            gmp_printf("c = %Zd\n", gc);
            exit(1);
        }


        mpz_clear(gc);
        mpz_clear(gn);


        exit(1);
#endif
        

    }
    else
    {
        // now subtract the hi part from the lo
        scarry = 0;
        __mmask8 carryout = 0;
        for (i = 0; i < wshift; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);   // hi
            b0 = _mm512_load_epi64(s->data + i * VECLEN);   // lo
            
            //a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
            __m512i t = _mm512_sub_epi64(b0, a1);
            carryout = _mm512_cmpgt_epu64_mask(a1, b0);
            __m512i t2 = _mm512_sub_epi64(t, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_kor(carryout, _mm512_cmpgt_epu64_mask(t2, t));
            a0 = _mm512_and_epi64(t2, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(s->data + i * VECLEN);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add 1.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            
            //a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            __m512i t = _mm512_add_epi64(a1, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpeq_epu64_mask(a1, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // and add one to the hi-bit.  This should resolve the borrow bit.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        //a1 = _mm512_add_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift))), a1);
        //_mm512_store_epi64(c->data + wshift * VECLEN,
        //    _mm512_and_epi64(_mm512_set1_epi64(0xfffffffffffffULL), a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif

        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));
    }

#ifdef DEBUG_MERSENNE
    mpz_clear(ga);
    mpz_clear(gb);
#endif

    c->size = NWORDS;
    return;
}

void vecsqrmod52_mersenne(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    vecmulmod52_mersenne(a, a, c, n, s, mdata);
    return;

    // 8x sqr:
    // input 8 bignums in the even lanes of a.
    // output 8 squaremod bignums in the even lanes of c.
    int i, j, k;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t NBLOCKS = mdata->NBLOCKS;
    vec_bignum_t* b = a;

    //if (NWORDS == 8)
    //{
    //    vecsqrmod52_mersenne_416(a, c, n, s, mdata);
    //    return;
    //}

    // needed in loops
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19

#ifndef IFMA
    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
#endif

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry;

    // deal with the sign
    c->size = NWORDS;
    c->signmask = 0;

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

#ifdef DEBUG_MERSENNE
    printf("commencing vecsqrmod52_mersenne\n");
    print_vechexbignum(a, "input a: ");
    print_vechexbignum(b, "input b: ");
    print_vechexbignum(s, "input s: ");
#endif

    // first half sqr
    for (i = 0; i < NBLOCKS; i++)
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (k = 0, j = i; j > (i + 1) / 2; j--, k++)
        {
            // when i = 0, j = 0, j > 0 --> 0 iterations
            // when i = 1, j = 1, j > 1 --> 0 iterations
            // when i = 2, j = 2, j > 1 --> 1 iteration @ a[3..0], b[5..11]
            // when i = 3, j = 3, j > 2 --> 1 iteration @ a[3..0], b[9..15]
            // when i = 4, j = 4, j > 2 --> 2 iteration @ a[3..0], b[13..19] and a[7..4], b[9..15]
            // when i = 5, j = 5, j > 3 --> 2 iteration @ a[3..0], b[17..23] and a[7..4], b[13..19]
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        if (i & 1)
        {
            // i odd
            // for 384-bit inputs when i == 1, j = 0, a = {3,2,1,0} and b = {2,3,4,5,6,7}
            // for 512-bit inputs when i == 3, j = 1, a = {7,6,5,4} and b = {6,7,8,9,a,b}
            // for 512-bit inputs when i == 1, j = 0, a = {3,2,1,0} and b = {2,3,4,5,6,7}

            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);


#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a2, b2, te0, te1);
            VEC_MUL_ACCUM_LOHI_PD(a1, b2, te2, te3);
            VEC_MUL_ACCUM_LOHI_PD(a1, b3, te4, te5);
            VEC_MUL_ACCUM_LOHI_PD(a0, b3, te6, te7);
#else
            //prod1_e = _mm512_mul_epu32(a2, b2);   // te0
            //prod2_e = _mm512_mul_epu32(a1, b2);   // te2
            //prod3_e = _mm512_mul_epu32(a1, b3);   // te4
            //prod4_e = _mm512_mul_epu32(a0, b3);   // te6
            //ACCUM_4X_DOUBLED_PROD;
            {
                prod5_ld = _mm512_cvtepu64_pd(a0);
                prod1_ld = _mm512_cvtepu64_pd(a1);
                prod2_ld = _mm512_cvtepu64_pd(a2);
                prod3_ld = _mm512_cvtepu64_pd(b2);
                prod4_ld = _mm512_cvtepu64_pd(b3);

                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b2 -> to te0/1 
                prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b2 -> to te2/3
                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b3 -> to te4/5
                prod4_hd = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b3 -> to te6/7

                te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod1_hd));
                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod4_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
                prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod4_ld = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod1_ld));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod4_ld));
            }
#endif


#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a3, b3, te0, te1);
            VEC_MUL_ACCUM_LOHI_PD(a2, b3, te2, te3);
            VEC_MUL_ACCUM_LOHI_PD(a2, b4, te4, te5);
            VEC_MUL_ACCUM_LOHI_PD(a1, b4, te6, te7);
#else
            //prod1_e = _mm512_mul_epu32(a3, b3);   // te0
            //prod2_e = _mm512_mul_epu32(a2, b3);   // te2
            //prod3_e = _mm512_mul_epu32(a2, b4);   // te4
            //prod4_e = _mm512_mul_epu32(a1, b4);   // te6
            //ACCUM_4X_DOUBLED_PROD;
            {
                prod5_ld = _mm512_cvtepu64_pd(a1);
                prod1_ld = _mm512_cvtepu64_pd(a2);
                prod2_ld = _mm512_cvtepu64_pd(a3);
                prod3_ld = _mm512_cvtepu64_pd(b3);
                prod4_ld = _mm512_cvtepu64_pd(b4);

                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * b3 -> to te0/1
                prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b3 -> to te2/3
                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b4 -> to te4/5
                prod4_hd = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b4 -> to te6/7

                te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod1_hd));
                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod4_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
                prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod4_ld = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod1_ld));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod4_ld));
            }
#endif


#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a3, b4, te2, te3);
            VEC_MUL_ACCUM_LOHI_PD(a3, b5, te4, te5);
            VEC_MUL_ACCUM_LOHI_PD(a2, b5, te6, te7);
            VEC_MUL_ACCUM_LOHI_PD(a3, b6, te6, te7);
#else
            //prod1_e = _mm512_mul_epu32(a3, b4);  // te2
            //prod2_e = _mm512_mul_epu32(a3, b5);  // te4
            //prod3_e = _mm512_mul_epu32(a2, b5);  // te6
            {
                prod1_ld = _mm512_cvtepu64_pd(a2);
                prod2_ld = _mm512_cvtepu64_pd(a3);
                prod3_ld = _mm512_cvtepu64_pd(b4);
                prod4_ld = _mm512_cvtepu64_pd(b5);

                prod3_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * b4 -> to te2/3
                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * b5 -> to te4/5
                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b5 -> to te6/7

                te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));
                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);

                prod3_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            }

            //prod1_e = _mm512_mul_epu32(a3, b6);  // te6
            {
                prod1_ld = _mm512_cvtepu64_pd(a3);
                prod2_ld = _mm512_cvtepu64_pd(b6);

                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * b6 -> to te6/7

                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            }

            // all terms so far need to be doubled.  
            // but to do that we first need to remove all bias
            // that has been accumulated so far.
            SUB_BIAS_HI(
                k * 4 + 2,
                k * 4 + 3,
                k * 4 + 3,
                k * 4 + 4);
            SUB_BIAS_LO(
                k * 4 + 2,
                k * 4 + 3,
                k * 4 + 3,
                k * 4 + 4);
#endif

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a1, a1, te0, te1);
            VEC_MUL_ACCUM_LOHI_PD(a0, a0, te4, te5);
#else
            // finally, accumulate the two non-doubled terms.
            //prod1_e = _mm512_mul_epu32(a1, a1);    // te0
            //prod2_e = _mm512_mul_epu32(a0, a0);    // te4
            {
                prod1_ld = _mm512_cvtepu64_pd(a0);
                prod2_ld = _mm512_cvtepu64_pd(a1);

                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * a1 -> to te0/1
                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a0 -> to te4/5

                te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod1_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod1_ld));
            }
#endif

        }
        else
        {
            // i even
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);

#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a0, a1, te2, te3);
            VEC_MUL_ACCUM_LOHI_PD(a0, a2, te4, te5);
            VEC_MUL_ACCUM_LOHI_PD(a0, a3, te6, te7);
            VEC_MUL_ACCUM_LOHI_PD(a1, a2, te6, te7);
#else
            //prod1_e = _mm512_mul_epu32(a0, a1);    // te2
            //prod2_e = _mm512_mul_epu32(a0, a2);    // te4
            //prod3_e = _mm512_mul_epu32(a0, a3);    // te6
            {
                prod1_ld = _mm512_cvtepu64_pd(a0);
                prod2_ld = _mm512_cvtepu64_pd(a1);
                prod3_ld = _mm512_cvtepu64_pd(a2);
                prod4_ld = _mm512_cvtepu64_pd(a3);

                prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a1 -> to te2/3
                prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a2 -> to te4/5
                prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a3 -> to te6/7

                te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod2_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod3_hd));
                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod4_hd));

                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
                prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

                prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod2_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod3_ld));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod4_ld));
            }

            //prod1_e = _mm512_mul_epu32(a1, a2);    // te6
            {
                prod1_ld = _mm512_cvtepu64_pd(a1);
                prod2_ld = _mm512_cvtepu64_pd(a2);

                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * a2 -> to te6/7

                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            }

            // all terms so far need to be doubled.  
            // but to do that we first need to remove all bias
            // that has been accumulated so far.
            SUB_BIAS_HI(
                k * 4 + 0,
                k * 4 + 1,
                k * 4 + 1,
                k * 4 + 2);
            SUB_BIAS_LO(
                k * 4 + 0,
                k * 4 + 1,
                k * 4 + 1,
                k * 4 + 2);
#endif

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a0, a0, te0, te1);
            VEC_MUL_ACCUM_LOHI_PD(a1, a1, te4, te5);
#else
            // finally, accumulate the two non-doubled terms.
            //prod1_e = _mm512_mul_epu32(a0, a0);
            //prod2_e = _mm512_mul_epu32(a1, a1);
            {
                prod1_ld = _mm512_cvtepu64_pd(a0);
                prod2_ld = _mm512_cvtepu64_pd(a1);

                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a0 -> to te0/1
                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * a1 -> to te4/5

                te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
            }
#endif

        }

#ifndef IFMA
        // need to remove bias from the two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            1,
            0,
            1,
            0);
        SUB_BIAS_LO(
            1,
            0,
            1,
            0);
#endif

        // final column accumulation
        {
            j = 0;
            // accumulate this column-sum:
            // carry propagate low to high.
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 1;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te2);
            acc_e1 = _mm512_add_epi64(acc_e1, te3);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 2;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te4);
            acc_e1 = _mm512_add_epi64(acc_e1, te5);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 3;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te6);
            acc_e1 = _mm512_add_epi64(acc_e1, te7);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;
        }
    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(s, "after lo half:");
#endif

    // second half sqr
    for (i = 0; i < NBLOCKS; i++)
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (k = 0, j = 0; j < (NBLOCKS - i - 1) / 2; j++, k++)
        {
            // Compute a solid block (all matching terms are in the lower
            // half triangle of the expansion).
            a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
            a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
            a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);
            a3 = _mm512_load_epi64(a->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);

            b0 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // The final block shape depends on the parity of i and NBLOCKS
        // Block shape 1 (small upper triangle) if 'i' is odd and NBLOCKS is even,
        // or if 'i' is even and NBLOCKS is odd.
        // Block shape 2 if 'i' is odd and NBLOCKS is odd or if 'i' is even
        // and NBLOCKS is even.
        if (NBLOCKS & 1)		// NBLOCKS is odd
        {
            // i odd, block shape 2.
            if (i & 1)
            {
                // always a continuation of the full-block loop, so use the same 
                // loading pattern.  Only now we don't need as many b-terms.
                a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
                a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
                a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);
                a3 = _mm512_load_epi64(a->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);

                b0 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 1) * VECLEN);
                b1 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 2) * VECLEN);
                b2 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 3) * VECLEN);
                b3 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 4) * VECLEN);
                b4 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 5) * VECLEN);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a0, b0, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a1, b2, te2, te3);
                VEC_MUL_ACCUM_LOHI_PD(a1, b1, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a0, b1, te2, te3);
#else
                //prod1_e = _mm512_mul_epu32(a0, b0);        // te0
                //prod1_e = _mm512_mul_epu32(a1, b1);        // te0
                //prod1_e = _mm512_mul_epu32(a0, b1);        // te2
                //prod1_e = _mm512_mul_epu32(a1, b2);        // te2
                {
                    prod5_ld = _mm512_cvtepu64_pd(a0);
                    prod1_ld = _mm512_cvtepu64_pd(a1);
                    prod2_ld = _mm512_cvtepu64_pd(b0);
                    prod3_ld = _mm512_cvtepu64_pd(b1);
                    prod4_ld = _mm512_cvtepu64_pd(b2);

                    prod2_hd = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b0 -> to te0/1
                    prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b2 -> to te2/3
                    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b1 -> to te0/1
                    prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b1 -> to te2/3

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod4_hd));
                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));
                    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));

                    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
                    prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

                    prod2_ld = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod4_ld));
                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
                    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                }
#endif

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a2, b2, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a2, b3, te2, te3);
#else
                //prod1_e = _mm512_mul_epu32(a2, b2);        // te0
                //prod1_e = _mm512_mul_epu32(a2, b3);        // te2
                {
                    prod1_ld = _mm512_cvtepu64_pd(a2);
                    prod2_ld = _mm512_cvtepu64_pd(b2);
                    prod3_ld = _mm512_cvtepu64_pd(b3);

                    prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b2 -> to te0/1
                    prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b3 -> to te2/3

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));

                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);

                    prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                }
#endif

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a0, b2, te4, te5);
                VEC_MUL_ACCUM_LOHI_PD(a1, b4, te6, te7);
                VEC_MUL_ACCUM_LOHI_PD(a1, b3, te4, te5);
                VEC_MUL_ACCUM_LOHI_PD(a0, b3, te6, te7);
#else
                // prod1_e = _mm512_mul_epu32(a0, b2);       // te4
                // prod1_e = _mm512_mul_epu32(a1, b3);       // te4
                // prod1_e = _mm512_mul_epu32(a0, b3);       // te6
                // prod1_e = _mm512_mul_epu32(a1, b4);       // te6
                {
                    prod5_ld = _mm512_cvtepu64_pd(a0);
                    prod1_ld = _mm512_cvtepu64_pd(a1);
                    prod2_ld = _mm512_cvtepu64_pd(b2);
                    prod3_ld = _mm512_cvtepu64_pd(b3);
                    prod4_ld = _mm512_cvtepu64_pd(b4);

                    prod2_hd = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b2 -> to te4/5
                    prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b4 -> to te6/7
                    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b3 -> to te4/5
                    prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b3 -> to te6/7

                    te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));
                    te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod4_hd));
                    te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod1_hd));
                    te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod3_hd));

                    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
                    prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

                    prod2_ld = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
                    te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod4_ld));
                    te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod1_ld));
                    te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod3_ld));
                }

                // all terms so far need to be doubled.  
                // but to do that we first need to remove all bias
                // that has been accumulated so far.
                SUB_BIAS_HI(
                    k * 4 + 3,
                    k * 4 + 3,
                    k * 4 + 2,
                    k * 4 + 2);
                SUB_BIAS_LO(
                    k * 4 + 3,
                    k * 4 + 3,
                    k * 4 + 2,
                    k * 4 + 2);
#endif

                // now double
                te1 = _mm512_slli_epi64(te1, 1);
                te3 = _mm512_slli_epi64(te3, 1);
                te5 = _mm512_slli_epi64(te5, 1);
                te7 = _mm512_slli_epi64(te7, 1);
                te0 = _mm512_slli_epi64(te0, 1);
                te2 = _mm512_slli_epi64(te2, 1);
                te4 = _mm512_slli_epi64(te4, 1);
                te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a3, a3, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a2, a2, te4, te5);
#else
                // finally the two non-doubled terms.
                //prod1_e = _mm512_mul_epu32(a3, a3);    // te0
                //prod1_e = _mm512_mul_epu32(a2, a2);    // te4
                {
                    prod1_ld = _mm512_cvtepu64_pd(a3);
                    prod2_ld = _mm512_cvtepu64_pd(a2);

                    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * a3 -> to te0/1
                    prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * a2 -> to te4/5

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));
                    te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));

                    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
                    te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
                }
#endif

            }
            else
            {
                // i even, block shape 1.
                // always a continuation of the full-block loop, so use the same 
                // loading pattern.  Only now we don't need as many b-terms.
                a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
                a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
                a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a0, a2, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a0, a1, te2, te3);
#else
                //k == 0;
                //prod1_e = _mm512_mul_epu32(a0, a2);      // te0
                //prod1_e = _mm512_mul_epu32(a0, a1);      // te2
                {
                    prod1_ld = _mm512_cvtepu64_pd(a0);
                    prod2_ld = _mm512_cvtepu64_pd(a1);
                    prod3_ld = _mm512_cvtepu64_pd(a2);

                    prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a2 -> to te0/1
                    prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a1 -> to te2/3

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod3_hd));
                    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod2_hd));

                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);

                    prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod3_ld));
                    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod2_ld));
                }

                // all terms so far need to be doubled.  
                // but to do that we first need to remove all bias
                // that has been accumulated so far.
                SUB_BIAS_HI(
                    k * 4 + 1,
                    k * 4 + 1,
                    k * 4 + 0,
                    k * 4 + 0);
                SUB_BIAS_LO(
                    k * 4 + 1,
                    k * 4 + 1,
                    k * 4 + 0,
                    k * 4 + 0);
#endif

                // now double
                te1 = _mm512_slli_epi64(te1, 1);
                te3 = _mm512_slli_epi64(te3, 1);
                te5 = _mm512_slli_epi64(te5, 1);
                te7 = _mm512_slli_epi64(te7, 1);
                te0 = _mm512_slli_epi64(te0, 1);
                te2 = _mm512_slli_epi64(te2, 1);
                te4 = _mm512_slli_epi64(te4, 1);
                te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a1, a1, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a0, a0, te4, te5);
#else
                // finally the two non-doubled terms.
                //prod1_e = _mm512_mul_epu32(a1, a1);    // te0
                //prod1_e = _mm512_mul_epu32(a0, a0);    // te4
                {
                    prod1_ld = _mm512_cvtepu64_pd(a1);
                    prod2_ld = _mm512_cvtepu64_pd(a0);

                    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b1 -> to te0/1
                    prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b0 -> to te4/5

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));
                    te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));

                    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
                    te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
                }
#endif

            }

        }
        else				// NBLOCKS is even
        {
            // i odd, block shape 1.
            if (i & 1)
            {
                // always a continuation of the full-block loop, so use the same 
                // loading pattern.  Only now we don't need as many b-terms.
                a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);  // {f, b}
                a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);  // {e, a}
                a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);  // {d, 9}

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a0, a2, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a0, a1, te2, te3);
#else
                //k == 0;
                //prod1_e = _mm512_mul_epu32(a0, a2);   // te0
                //prod2_e = _mm512_mul_epu32(a0, a1);   // te2
                {
                    prod1_ld = _mm512_cvtepu64_pd(a0);
                    prod2_ld = _mm512_cvtepu64_pd(a1);
                    prod3_ld = _mm512_cvtepu64_pd(a2);

                    prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a2 -> to te0/1
                    prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a1 -> to te2/3

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod3_hd));
                    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod2_hd));

                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);

                    prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod3_ld));
                    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod2_ld));
                }

                // all terms so far need to be doubled.  
                // but to do that we first need to remove all bias
                // that has been accumulated so far.
                SUB_BIAS_HI(
                    j * 4 + 1,
                    j * 4 + 1,
                    j * 4 + 0,
                    j * 4 + 0);
                SUB_BIAS_LO(
                    j * 4 + 1,
                    j * 4 + 1,
                    j * 4 + 0,
                    j * 4 + 0);
#endif

                // now double
                te1 = _mm512_slli_epi64(te1, 1);
                te3 = _mm512_slli_epi64(te3, 1);
                te5 = _mm512_slli_epi64(te5, 1);
                te7 = _mm512_slli_epi64(te7, 1);
                te0 = _mm512_slli_epi64(te0, 1);
                te2 = _mm512_slli_epi64(te2, 1);
                te4 = _mm512_slli_epi64(te4, 1);
                te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a1, a1, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a0, a0, te4, te5);
#else
                // finally, accumulate the two non-doubled terms.
                //prod1_e = _mm512_mul_epu32(a1, a1);       // te0
                //prod2_e = _mm512_mul_epu32(a0, a0);       // te4
                {
                    prod1_ld = _mm512_cvtepu64_pd(a1);
                    prod2_ld = _mm512_cvtepu64_pd(a0);

                    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b1 -> to te0/1
                    prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b0 -> to te4/5

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));
                    te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));

                    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
                    te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
                }
#endif
            }
            else
            {
                // i even, block shape 1.
                // always a continuation of the full-block loop, so use the same 
                // loading pattern.  Only now we don't need as many b-terms.
                a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);		// {f, b}
                a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);		// {e, a}
                a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);		// {d, 9}
                a3 = _mm512_load_epi64(a->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);		// {c, 8}

                b0 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 1) * VECLEN); // {9, 5}
                b1 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 2) * VECLEN);	// {a, 6}
                b2 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 3) * VECLEN);	// {b, 7}
                b3 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 4) * VECLEN);	// {c, 8}
                b4 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 5) * VECLEN); // {d, 9}

                //prod1_e = _mm512_mul_epu32(a0, b0);    // te0
                //prod2_e = _mm512_mul_epu32(a0, b1);    // te2
                //prod3_e = _mm512_mul_epu32(a0, b2);    // te4
                //prod4_e = _mm512_mul_epu32(a0, b3);    // te6
                //ACCUM_4X_DOUBLED_PROD;
                VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);

                //prod1_e = _mm512_mul_epu32(a1, b1);    // te0
                //prod2_e = _mm512_mul_epu32(a1, b2);    // te2
                //prod3_e = _mm512_mul_epu32(a1, b3);    // te4
                //prod4_e = _mm512_mul_epu32(a1, b4);    // te6
                //ACCUM_4X_DOUBLED_PROD;
                VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a2, b2, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a2, b3, te2, te3);
#else
                //prod1_e = _mm512_mul_epu32(a2, b2);    // te0
                //prod2_e = _mm512_mul_epu32(a2, b3);    // te2
                {
                    prod1_ld = _mm512_cvtepu64_pd(a2);
                    prod2_ld = _mm512_cvtepu64_pd(b2);
                    prod3_ld = _mm512_cvtepu64_pd(b3);

                    prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b2 -> to te0/1
                    prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b3 -> to te2/3

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));

                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);

                    prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                }

                // all terms so far need to be doubled.  
                // but to do that we first need to remove all bias
                // that has been accumulated so far.
                SUB_BIAS_HI(
                    j * 4 + 3,
                    j * 4 + 3,
                    j * 4 + 2,
                    j * 4 + 2);
                SUB_BIAS_LO(
                    j * 4 + 3,
                    j * 4 + 3,
                    j * 4 + 2,
                    j * 4 + 2);
#endif

                // now double
                te1 = _mm512_slli_epi64(te1, 1);
                te3 = _mm512_slli_epi64(te3, 1);
                te5 = _mm512_slli_epi64(te5, 1);
                te7 = _mm512_slli_epi64(te7, 1);
                te0 = _mm512_slli_epi64(te0, 1);
                te2 = _mm512_slli_epi64(te2, 1);
                te4 = _mm512_slli_epi64(te4, 1);
                te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a3, a3, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a2, a2, te4, te5);
#else
                // finally, accumulate the two non-doubled terms.
                //prod1_e = _mm512_mul_epu32(a3, a3);   // te0
                //prod2_e = _mm512_mul_epu32(a2, a2);   // te4
                {
                    prod1_ld = _mm512_cvtepu64_pd(a3);
                    prod2_ld = _mm512_cvtepu64_pd(a2);

                    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * a3 -> to te0/1
                    prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * a2 -> to te4/5

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));
                    te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));

                    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
                    te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
                }
#endif
            }
        }

#ifndef IFMA
        // need to remove bias from the two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            1,
            0,
            1,
            0);
        SUB_BIAS_LO(
            1,
            0,
            1,
            0);
#endif

        // final column accumulation
        {
            j = 0;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a3 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN, 
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 1;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te2);
            acc_e1 = _mm512_add_epi64(acc_e1, te3);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a2 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 2;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te4);
            acc_e1 = _mm512_add_epi64(acc_e1, te5);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a1 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 3;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te6);
            acc_e1 = _mm512_add_epi64(acc_e1, te7);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a0 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            a3 = _mm512_and_epi64(vlmask, a3);
            a2 = _mm512_and_epi64(vlmask, a2);
            a1 = _mm512_and_epi64(vlmask, a1);
            a0 = _mm512_and_epi64(vlmask, a0);

            _mm512_store_epi64(s->data + (i * BLOCKWORDS + NWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + NWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + NWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + NWORDS + 3) * VECLEN, a0);

        }


    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(s, "after hi half:");
#endif

    // reduce by adding hi to lo.  first right shift hi into output.
    __m512i vbshift, vbpshift;
    int bshift = mdata->nbits % 52;
    int wshift = mdata->nbits / 52;

    if (mdata->use_vnbits)
    {
        int bshift_a[VECLEN];
        for (i = 0; i < VECLEN; i++)
        {
            bshift_a[i] = mdata->vnbits[i] % 52;
        }

        //vbshift = _mm512_set_epi64(
        //    (1ULL << (uint64_t)(bshift_a[7])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[0])) - 1ULL);
        //
        //vbpshift = _mm512_set_epi64(
        //    (1ULL << (uint64_t)(bshift_a[7])),
        //    (1ULL << (uint64_t)(bshift_a[6])),
        //    (1ULL << (uint64_t)(bshift_a[5])),
        //    (1ULL << (uint64_t)(bshift_a[4])),
        //    (1ULL << (uint64_t)(bshift_a[3])),
        //    (1ULL << (uint64_t)(bshift_a[2])),
        //    (1ULL << (uint64_t)(bshift_a[1])),
        //    (1ULL << (uint64_t)(bshift_a[0])));

        vbshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[7])) - 1ULL);

        vbpshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])),
            (1ULL << (uint64_t)(bshift_a[1])),
            (1ULL << (uint64_t)(bshift_a[2])),
            (1ULL << (uint64_t)(bshift_a[3])),
            (1ULL << (uint64_t)(bshift_a[4])),
            (1ULL << (uint64_t)(bshift_a[5])),
            (1ULL << (uint64_t)(bshift_a[6])),
            (1ULL << (uint64_t)(bshift_a[7])));

        vec_bignum52_mask_rshift_vn(s, c, mdata->vnbits, 0xff);
    }
    else
    {
        vbshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL);
        vbpshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)));
        vec_bignum52_mask_rshift_n(s, c, mdata->nbits, 0xff);
    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(c, "hi part:");
    print_vechexbignum(s, "lo part:");
#endif

    if (mdata->isMersenne > 1)
    {
        // this number has been identified as a "pseudo"-Mersenne: 2^n-c, with
        // 'c' small.  Fast reduction is still possible.

        // multiply c * hi
        b0 = _mm512_set1_epi64(mdata->isMersenne);
        vecmul52_1(c, b0, c, n, s, mdata);

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after mul_1 with c:");
#endif

        // clear any hi-bits in the lo result
        a1 = _mm512_load_epi64(s->data + wshift * VECLEN);
        _mm512_store_epi64(s->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(s, "cleared low part:");
#endif

        // add c * hi + lo
        scarry = 0;
        for (i = 0; i <= NWORDS; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            b0 = _mm512_load_epi64(s->data + i * VECLEN);

            //a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            __m512i t = _mm512_add_epi64(a1, b0);
            t = _mm512_add_epi64(t, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpgt_epu64_mask(t, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "lo + hi * c:");
#endif

        // now we need to do it again, but the hi part is just a single word so
        // we only have a single-precision vecmul
        //vec_bignum52_mask_rshift_n(c, s, mdata->nbits, 0xff);
        //a0 = _mm512_load_epi64(s->data);

        // It's possible these hi bits are located in two words
        a0 = _mm512_load_epi64(c->data + wshift * VECLEN);
        a0 = _mm512_srli_epi64(a0, bshift);
        a1 = _mm512_load_epi64(c->data + (wshift + 1) * VECLEN);
        a0 = _mm512_or_epi64(a0, _mm512_slli_epi64(a1, 52 - bshift));

        //print_regvechex(a0, 0, "hi part:");

        b0 = _mm512_set1_epi64(mdata->isMersenne);

        //print_regvechex(b0, 0, "c:");

        VEC_MUL_LOHI_PD(a0, b0, acc_e0, acc_e1);

        // clear any hi-bits now that we have the multiplier
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_set1_epi64(0));
        }

        // now add that into the lo part again.  Here we add in 
        // the two words of the hi-mul result.
        b0 = _mm512_load_epi64(c->data + 0 * VECLEN);
        b1 = _mm512_load_epi64(c->data + 1 * VECLEN);
        b2 = zero;
        b3 = zero;
        acc_e0 = _mm512_add_epi64(acc_e0, b0);
        VEC_CARRYPROP_LOHI(acc_e0, b2);         // lo word and carry for 1st add.
        acc_e1 = _mm512_add_epi64(acc_e1, b1);
        VEC_CARRYPROP_LOHI(acc_e1, b3);         // lo word and carry for 2nd add.
        acc_e1 = _mm512_add_epi64(acc_e1, b2);  // add 1st carry to 2nd lo word.
        VEC_CARRYPROP_LOHI(acc_e1, b3);         // that could have a carry, add it to b3.

        // save first two result words 
        _mm512_store_epi64(c->data + 0 * VECLEN, acc_e0);
        _mm512_store_epi64(c->data + 1 * VECLEN, acc_e1);

        // propagate the carry through the rest of the lo result.
        b2 = _mm512_load_epi64(c->data + 2 * VECLEN);
        scarry = 0;
        a0 = _mm512_adc_epi52(b2, scarry, b3, &scarry);
        _mm512_store_epi64(c->data + 2 * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i = 3; (i < wshift) && (scarry > 0); i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);

            //a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            __m512i t = _mm512_add_epi64(a1, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpeq_epu64_mask(a1, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "lo + hi * c:");
#endif

        // it is unlikely, but possible that we still have a final carry to propagate
        if (scarry > 0)
        {
            printf("unlikely add\n");
            b0 = _mm512_set1_epi64(mdata->isMersenne);
            a0 = _mm512_load_epi64(c->data + 0 * VECLEN);
            a0 = _mm512_adc_epi52(a0, 0, b0, &scarry);
            _mm512_store_epi64(c->data + 0 * VECLEN, a0);

            i = 1;
            while (scarry > 0)
            {
                a1 = _mm512_load_epi64(c->data + i * VECLEN);

                //a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
                __m512i t = _mm512_add_epi64(a1, _mm512_maskz_set1_epi64(scarry, 1));
                scarry = _mm512_cmpeq_epu64_mask(a1, vlmask);
                a0 = _mm512_and_epi64(t, vlmask);

                _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
                i++;
            }
        }

        // finally clear any hi-bits in the result
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif

    }
    else if (mdata->isMersenne > 0)
    {
        // now add the low part into the high.
        scarry = 0;
        for (i = 0; i < wshift; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            b0 = _mm512_load_epi64(s->data + i * VECLEN);

            //a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            __m512i t = _mm512_add_epi64(a1, b0);
            t = _mm512_add_epi64(t, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpgt_epu64_mask(t, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(s->data + i * VECLEN);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add it back in.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);

            //a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            __m512i t = _mm512_add_epi64(a1, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpeq_epu64_mask(a1, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // clear the potential hi-bit
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif
    }
    else
    {
        // now subtract the hi part from the lo
        scarry = 0;
        __mmask8 carryout = 0;
        for (i = 0; i < wshift; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);   // hi
            b0 = _mm512_load_epi64(s->data + i * VECLEN);   // lo

            //a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
            __m512i t = _mm512_sub_epi64(b0, a1);
            carryout = _mm512_cmpgt_epu64_mask(a1, b0);
            __m512i t2 = _mm512_sub_epi64(t, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_kor(carryout, _mm512_cmpgt_epu64_mask(t2, t));
            a0 = _mm512_and_epi64(t2, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(s->data + i * VECLEN);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add 1.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);

            //a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            __m512i t = _mm512_add_epi64(a1, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpeq_epu64_mask(a1, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // and add one to the hi-bit.  This should resolve the borrow bit.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        //a1 = _mm512_add_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift))), a1);
        //_mm512_store_epi64(c->data + wshift * VECLEN,
        //    _mm512_and_epi64(_mm512_set1_epi64(0xfffffffffffffULL), a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif

        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));
    }

    c->size = NWORDS;
    return;
}

void vecaddmod52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata)
{
    // assumptions:
    // a, b, c are of length VECLEN * NWORDS
    // a, b, c, and n are aligned
    // a and b are both positive
    // n is the montgomery base
    int i;
    uint32_t NWORDS = a->WORDS_ALLOC;
    __mmask8 carry = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    __m512i lomask = _mm512_set1_epi64(0xfffffffffffffULL);
    __m512i zero = _mm512_set1_epi64(0);

    __m512i vbshift, vbpshift;
    int bshift = mdata->nbits % 52;
    int wshift = mdata->nbits / 52;

    if (mdata->use_vnbits)
    {
        int bshift_a[VECLEN];
        for (i = 0; i < VECLEN; i++)
        {
            bshift_a[i] = mdata->vnbits[i] % 52;
        }

        //vbshift = _mm512_set_epi64(
        //    (1ULL << (uint64_t)(bshift_a[7])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[0])) - 1ULL);
        //
        //vbpshift = _mm512_set_epi64(
        //    (1ULL << (uint64_t)(bshift_a[7])),
        //    (1ULL << (uint64_t)(bshift_a[6])),
        //    (1ULL << (uint64_t)(bshift_a[5])),
        //    (1ULL << (uint64_t)(bshift_a[4])),
        //    (1ULL << (uint64_t)(bshift_a[3])),
        //    (1ULL << (uint64_t)(bshift_a[2])),
        //    (1ULL << (uint64_t)(bshift_a[1])),
        //    (1ULL << (uint64_t)(bshift_a[0])));

        vbshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[7])) - 1ULL);

        vbpshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])),
            (1ULL << (uint64_t)(bshift_a[1])),
            (1ULL << (uint64_t)(bshift_a[2])),
            (1ULL << (uint64_t)(bshift_a[3])),
            (1ULL << (uint64_t)(bshift_a[4])),
            (1ULL << (uint64_t)(bshift_a[5])),
            (1ULL << (uint64_t)(bshift_a[6])),
            (1ULL << (uint64_t)(bshift_a[7])));
    }
    else
    {
        vbshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL);
        vbpshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)));
    }

    // add
    for (i = 0; i < NWORDS; i++)
    {
        avec = _mm512_load_epi64(a->data + i * VECLEN);
        bvec = _mm512_load_epi64(b->data + i * VECLEN);

        //cvec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        __m512i t = _mm512_add_epi64(avec, bvec);
        t = _mm512_add_epi64(t, _mm512_maskz_set1_epi64(carry, 1));
        carry = _mm512_cmpgt_epu64_mask(t, lomask);
        cvec = _mm512_and_epi64(t, lomask);

        _mm512_store_epi64(c->data + i * VECLEN, cvec);
    }

    // check for a carry.
    avec = _mm512_load_epi64(c->data + wshift * VECLEN);
    carry = _mm512_test_epi64_mask(avec, vbpshift);

    if (mdata->isMersenne > 0)
    {
        // the modulo is just the low part plus 1 (the carry, if present).
        cvec = _mm512_load_epi64(c->data + 0 * VECLEN);
        //bvec = _mm512_addcarry_epi52(cvec, carry, &carry);

        bvec = _mm512_set1_epi64(mdata->isMersenne);
        bvec = _mm512_mask_adc_epi52(cvec, carry, 0, bvec, &carry);
        _mm512_store_epi64(c->data + 0 * VECLEN, bvec);

        for (i = 1; (i < NWORDS) && (carry > 0); i++)
        {
            cvec = _mm512_load_epi64(c->data + i * VECLEN);

            //bvec = _mm512_addcarry_epi52(cvec, carry, &carry);
            __m512i t = _mm512_add_epi64(cvec, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpeq_epu64_mask(cvec, lomask);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(c->data + i * VECLEN, bvec);
        }
    }
    else
    {
        // the modulo is just the low part minus 1 (the carry, if present).
        cvec = _mm512_load_epi64(c->data + 0 * VECLEN);
        bvec = _mm512_set1_epi64(-mdata->isMersenne);
        bvec = _mm512_mask_sbb_epi52(cvec, carry, 0, bvec, &carry);
        _mm512_store_epi64(c->data + 0 * VECLEN, bvec);

        for (i = 1; (i < NWORDS) && (carry > 0); i++)
        {
            cvec = _mm512_load_epi64(c->data + i * VECLEN);
            
            //bvec = _mm512_subborrow_epi52(cvec, carry, &carry);
            __m512i t = _mm512_sub_epi64(cvec, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpeq_epu64_mask(cvec, zero);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(c->data + i * VECLEN, bvec);
        }
    }

    // clear the potential hi-bit
    avec = _mm512_load_epi64(c->data + wshift * VECLEN);
    _mm512_store_epi64(c->data + wshift * VECLEN,
        _mm512_and_epi64(vbshift, avec));

    return;
}

void vecsubmod52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata)
{
    // assumptions:
    // a, b, c are of length VECLEN * NWORDS
    // s1 is of length VECLEN
    // a, b, c, n, and s1 are aligned
    // a and b are both positive
    // a >= b
    // n is the montgomery base
    int i;
    uint32_t NWORDS = a->WORDS_ALLOC;
    __mmask8 carry = 0;
    __mmask8 carryout = 0;
    __mmask8 mask = 0;
    __mmask8 mask2 = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    __m512i nvec;
    __m512i lomask = _mm512_set1_epi64(0xfffffffffffffULL);
    __m512i zero = _mm512_set1_epi64(0);

    __m512i vbshift;
    int bshift = mdata->nbits % 52;
    int wshift = mdata->nbits / 52;

    if (mdata->use_vnbits)
    {
        int bshift_a[VECLEN];
        for (i = 0; i < VECLEN; i++)
        {
            bshift_a[i] = mdata->vnbits[i] % 52;
        }

        //vbshift = _mm512_set_epi64(
        //    (1ULL << (uint64_t)(bshift_a[7])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[0])) - 1ULL);
        //

        vbshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[7])) - 1ULL);
    }
    else
    {
        vbshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL);
    }

    // subtract
    carry = 0;
    //for (i = 0; i < NWORDS; i++)
    for (i = 0; i <= wshift; i++)
    {
        avec = _mm512_load_epi64(a->data + i * VECLEN);
        bvec = _mm512_load_epi64(b->data + i * VECLEN);
        
        //cvec = _mm512_sbb_epi52(avec, carry, bvec, &carry);
        __m512i t = _mm512_sub_epi64(avec, bvec);
        carryout = _mm512_cmpgt_epu64_mask(bvec, avec);
        __m512i t2 = _mm512_sub_epi64(t, _mm512_maskz_set1_epi64(carry, 1));
        carry = _mm512_kor(carryout, _mm512_cmpgt_epu64_mask(t2, t));
        cvec = _mm512_and_epi64(t2, lomask);

        _mm512_store_epi64(c->data + i * VECLEN, cvec);
    }

    if (mdata->isMersenne > 0)
    {
        // if we had a final carry, then b was bigger than a so we need to re-add n.
        nvec = _mm512_set1_epi64(0x10000000000000ULL);
        nvec = _mm512_sub_epi64(nvec, _mm512_set1_epi64(mdata->isMersenne));
        cvec = _mm512_load_epi64(c->data + 0 * VECLEN);
        bvec = _mm512_mask_adc_epi52(cvec, carry, 0, nvec, &carry);
        _mm512_store_epi64(c->data + 0 * VECLEN, bvec);
        nvec = lomask;
        for (i = 1; i <= wshift; i++)
        {
            cvec = _mm512_load_epi64(c->data + i * VECLEN);
            
            //bvec = _mm512_mask_adc_epi52(cvec, mask, carry, nvec, &carry);
            __m512i t = _mm512_add_epi64(cvec, nvec);
            t = _mm512_mask_add_epi64(cvec, mask, t, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpgt_epu64_mask(t, lomask);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(c->data + i * VECLEN, bvec);
        }
    }
    else
    {
        // if we had a final carry, then b was bigger than a so we need to re-add n.
        mask = carry;
        carry = 0;
        nvec = _mm512_set1_epi64(1);
        cvec = _mm512_load_epi64(c->data + 0 * VECLEN);
        bvec = _mm512_mask_adc_epi52(cvec, mask, carry, nvec, &carry);
        _mm512_store_epi64(c->data + 0 * VECLEN, bvec);
        for (i = 1; (i <= wshift) && (carry > 0); i++)
        {
            cvec = _mm512_load_epi64(c->data + i * VECLEN);
            
            //bvec = _mm512_addcarry_epi52(cvec, carry, &carry);
            __m512i t = _mm512_add_epi64(cvec, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpeq_epu64_mask(cvec, lomask);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(c->data + i * VECLEN, bvec);
        }
    }

    //for (; i < NWORDS; i++)
    //{
    //    _mm512_store_epi64(c->data + i * VECLEN, _mm512_set1_epi64(0));
    //}

    // clear the potential hi-bit
    avec = _mm512_load_epi64(c->data + wshift * VECLEN);
    _mm512_store_epi64(c->data + wshift * VECLEN,
        _mm512_and_epi64(vbshift, avec));

    return;
}

void vec_simul_addsub52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* sum, vec_bignum_t* diff,
    vec_monty_t* mdata)
{
    // assumptions:
    // a, b, c are of length VECLEN * NWORDS
    // a, b, c, and n are aligned
    // a and b are both positive
    // n is the montgomery base
    // produce sum = a + b and diff = a - b at the same time which
    // saves 3N loads (only have to load a,b, and n once)
    int i;
    uint32_t NWORDS = a->WORDS_ALLOC;
    __mmask8 carry = 0;
    __mmask8 carryout = 0;
    __mmask8 borrow = 0;
    __mmask8 bmask = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    __m512i nvec;
    __m512i lomask = _mm512_set1_epi64(0xfffffffffffffULL);
    __m512i zero = _mm512_set1_epi64(0);

    __m512i vbshift, vbpshift;
    int bshift = mdata->nbits % 52;
    int wshift = mdata->nbits / 52;

    if (mdata->use_vnbits)
    {
        int bshift_a[VECLEN];
        for (i = 0; i < VECLEN; i++)
        {
            bshift_a[i] = mdata->vnbits[i] % 52;
        }

        //vbshift = _mm512_set_epi64(
        //    (1ULL << (uint64_t)(bshift_a[7])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[0])) - 1ULL);
        //
        //vbpshift = _mm512_set_epi64(
        //    (1ULL << (uint64_t)(bshift_a[7])),
        //    (1ULL << (uint64_t)(bshift_a[6])),
        //    (1ULL << (uint64_t)(bshift_a[5])),
        //    (1ULL << (uint64_t)(bshift_a[4])),
        //    (1ULL << (uint64_t)(bshift_a[3])),
        //    (1ULL << (uint64_t)(bshift_a[2])),
        //    (1ULL << (uint64_t)(bshift_a[1])),
        //    (1ULL << (uint64_t)(bshift_a[0])));

        vbshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[7])) - 1ULL);

        vbpshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])),
            (1ULL << (uint64_t)(bshift_a[1])),
            (1ULL << (uint64_t)(bshift_a[2])),
            (1ULL << (uint64_t)(bshift_a[3])),
            (1ULL << (uint64_t)(bshift_a[4])),
            (1ULL << (uint64_t)(bshift_a[5])),
            (1ULL << (uint64_t)(bshift_a[6])),
            (1ULL << (uint64_t)(bshift_a[7])));
    }
    else
    {
        vbshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL);
        vbpshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)));
    }


    //for (i = 0; i < NWORDS; i++)
    for (i = 0; i <= wshift; i++)
    {
        // add/sub
        avec = _mm512_load_epi64(a->data + i * VECLEN);
        bvec = _mm512_load_epi64(b->data + i * VECLEN);
        
        //cvec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        //cvec = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
        _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, avec, bvec);
        
        _mm512_store_epi64(sum->data + i * VECLEN, avec);
        _mm512_store_epi64(diff->data + i * VECLEN, bvec);
    }

    bmask = borrow;     // result too small, need to add n

    // check for a carry.
    avec = _mm512_load_epi64(sum->data + wshift * VECLEN);
    carry = _mm512_test_epi64_mask(avec, vbpshift);

    // deal with the sum:
    if (mdata->isMersenne > 0)
    {
        // sum produced a carry, subtract n, which is all FF's, so
        // equivalent to adding 1.
        // the modulo is just the low part plus 1 (the carry, if present).
        cvec = _mm512_load_epi64(sum->data + 0 * VECLEN);
        //bvec = _mm512_addcarry_epi52(cvec, carry, &carry);
        bvec = _mm512_set1_epi64(mdata->isMersenne);
        bvec = _mm512_mask_adc_epi52(cvec, carry, 0, bvec, &carry);
        _mm512_store_epi64(sum->data + 0 * VECLEN, bvec);

        for (i = 1; (i < NWORDS) && (carry > 0); i++)
        {
            cvec = _mm512_load_epi64(sum->data + i * VECLEN);
            
            //bvec = _mm512_addcarry_epi52(cvec, carry, &carry);
            __m512i t = _mm512_add_epi64(cvec, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpeq_epu64_mask(cvec, lomask);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(sum->data + i * VECLEN, bvec);
        }
    }
    else
    {
        // the modulo is just the low part minus 1 (the carry, if present).
        cvec = _mm512_load_epi64(sum->data + 0 * VECLEN);
        //bvec = _mm512_subborrow_epi52(cvec, carry, &carry);
        //mm512_mask_sbb_epi52(__m512i a, __mmask8 m, __mmask8 c, __m512i b, __mmask8* cout)
        bvec = _mm512_set1_epi64(-mdata->isMersenne);
        bvec = _mm512_mask_sbb_epi52(cvec, carry, 0, bvec, &carry);
        _mm512_store_epi64(sum->data + 0 * VECLEN, bvec);

        for (i = 1; (i < NWORDS) && (carry > 0); i++)
        {
            cvec = _mm512_load_epi64(sum->data + i * VECLEN);
            
            //bvec = _mm512_subborrow_epi52(cvec, carry, &carry);
            __m512i t = _mm512_sub_epi64(cvec, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpeq_epu64_mask(cvec, zero);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(sum->data + i * VECLEN, bvec);
        }
    }

    // clear the potential hi-bit
    avec = _mm512_load_epi64(sum->data + wshift * VECLEN);
    _mm512_store_epi64(sum->data + wshift * VECLEN,
        _mm512_and_epi64(vbshift, avec));

    // deal with the diff:
    if (mdata->isMersenne > 0)
    {
        // if we had a final carry, then b was bigger than a so we need to re-add n.
        nvec = _mm512_set1_epi64(0x10000000000000ULL);
        nvec = _mm512_sub_epi64(nvec, _mm512_set1_epi64(mdata->isMersenne));
        cvec = _mm512_load_epi64(diff->data + 0 * VECLEN);
        bvec = _mm512_mask_adc_epi52(cvec, bmask, 0, nvec, &carry);
        _mm512_store_epi64(diff->data + 0 * VECLEN, bvec);
        nvec = _mm512_set1_epi64(0xfffffffffffffULL);
        for (i = 1; i <= wshift; i++)
        {
            cvec = _mm512_load_epi64(diff->data + i * VECLEN);

            //bvec = _mm512_mask_adc_epi52(cvec, bmask, carry, nvec, &carry);
            __m512i t = _mm512_add_epi64(cvec, nvec);
            t = _mm512_mask_add_epi64(cvec, bmask, t, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpgt_epu64_mask(t, lomask);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(diff->data + i * VECLEN, bvec);
        }
    }
    else
    {
        // diff produced a carry, add n, which is all 1000...0001, so
        // equivalent to adding 1.
        // if we had a final carry, then b was bigger than a so we need to re-add n.
        carry = 0;
        nvec = _mm512_set1_epi64(1);
        cvec = _mm512_load_epi64(diff->data + 0 * VECLEN);
        bvec = _mm512_mask_adc_epi52(cvec, bmask, carry, nvec, &carry);
        _mm512_store_epi64(diff->data + 0 * VECLEN, bvec);
        for (i = 1; (i <= wshift) && (carry > 0); i++)
        {
            cvec = _mm512_load_epi64(diff->data + i * VECLEN);

            //bvec = _mm512_addcarry_epi52(cvec, carry, &carry);
            __m512i t = _mm512_add_epi64(cvec, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpeq_epu64_mask(cvec, lomask);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(diff->data + i * VECLEN, bvec);
        }
    }

    //for (; i < NWORDS; i++)
    //{
    //    _mm512_store_epi64(diff->data + i * VECLEN, _mm512_set1_epi64(0));
    //}

    // clear the potential hi-bit
    avec = _mm512_load_epi64(diff->data + wshift * VECLEN);
    _mm512_store_epi64(diff->data + wshift * VECLEN,
        _mm512_and_epi64(vbshift, avec));


    return;
}

#define DEBUG_MERSENNE
void vecmod_mersenne(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* s, vec_monty_t* mdata)
{
    int i, j;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t NBLOCKS = mdata->NBLOCKS;
    // needed in loops
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19

#ifndef IFMA
    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
#endif

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry;

    // deal with the sign
    c->size = NWORDS;
    c->signmask = a->signmask;

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;


    // reduce by adding hi to lo.  first right shift hi into output.
    // for Mersenne's, nbits should be set to the exponent of the number.
    __m512i vbshift, vbpshift;
    int bshift = mdata->nbits % 52;
    int wshift = mdata->nbits / 52;

    mpz_t ga;
    mpz_init(ga);
    extract_bignum_from_vec_to_mpz(ga, a, 0, 2 * NWORDS);

    if (mdata->use_vnbits)
    {
        int bshift_a[VECLEN];
        for (i = 0; i < VECLEN; i++)
        {
            bshift_a[i] = mdata->vnbits[i] % 52;
        }

        vbshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[7])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[0])) - 1ULL);

        vbpshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[7])),
            (1ULL << (uint64_t)(bshift_a[6])),
            (1ULL << (uint64_t)(bshift_a[5])),
            (1ULL << (uint64_t)(bshift_a[4])),
            (1ULL << (uint64_t)(bshift_a[3])),
            (1ULL << (uint64_t)(bshift_a[2])),
            (1ULL << (uint64_t)(bshift_a[1])),
            (1ULL << (uint64_t)(bshift_a[0])));

        vec_bignum52_mask_rshift_vn(s, c, mdata->vnbits, 0xff);
    }
    else
    {
        vbshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL);
        vbpshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)));
        vec_bignum52_mask_rshift_n(a, c, mdata->nbits, 0xff);
    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(c, "hi part:");
    print_vechexbignum(a, "lo part:");
#endif

    if (mdata->isMersenne > 1)
    {
        printf("vecmod_mersenne does not yet handle pseudo-mersenne inputs\n");
        exit(1);

        /*
        // this number has been identified as a "pseudo"-Mersenne: 2^n-c, with
        // 'c' small.  Fast reduction is still possible.

        // multiply c * hi
        b0 = _mm512_set1_epi64(mdata->isMersenne);
        vecmul52_1(c, b0, c, n, s, mdata);

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after mul_1 with c:");
#endif

        // clear any hi-bits in the lo result
        a1 = _mm512_load_epi64(s->data + wshift * VECLEN);
        _mm512_store_epi64(s->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(s, "cleared low part:");
#endif

        // add c * hi + lo
        scarry = 0;
        for (i = 0; i <= NWORDS; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            b0 = _mm512_load_epi64(s->data + i * VECLEN);
            a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "lo + hi * c:");
#endif

        // now we need to do it again, but the hi part is just a single word so
        // we only have a single-precision vecmul
        //vec_bignum52_mask_rshift_n(c, s, mdata->nbits, 0xff);
        //a0 = _mm512_load_epi64(s->data);

        // It's possible these hi bits are located in two words
        a0 = _mm512_load_epi64(c->data + wshift * VECLEN);
        a0 = _mm512_srli_epi64(a0, bshift);
        a1 = _mm512_load_epi64(c->data + (wshift + 1) * VECLEN);
        a0 = _mm512_or_epi64(a0, _mm512_slli_epi64(a1, 52 - bshift));

        //print_regvechex(a0, 0, "hi part:");

        b0 = _mm512_set1_epi64(mdata->isMersenne);

        //print_regvechex(b0, 0, "c:");

        VEC_MUL_LOHI_PD(a0, b0, acc_e0, acc_e1);

        // clear any hi-bits now that we have the multiplier
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_set1_epi64(0));
        }

        // now add that into the lo part again.  Here we add in
        // the two words of the hi-mul result.
        b0 = _mm512_load_epi64(c->data + 0 * VECLEN);
        b1 = _mm512_load_epi64(c->data + 1 * VECLEN);
        b2 = zero;
        b3 = zero;
        acc_e0 = _mm512_add_epi64(acc_e0, b0);
        VEC_CARRYPROP_LOHI(acc_e0, b2);         // lo word and carry for 1st add.
        acc_e1 = _mm512_add_epi64(acc_e1, b1);
        VEC_CARRYPROP_LOHI(acc_e1, b3);         // lo word and carry for 2nd add.
        acc_e1 = _mm512_add_epi64(acc_e1, b2);  // add 1st carry to 2nd lo word.
        VEC_CARRYPROP_LOHI(acc_e1, b3);         // that could have a carry, add it to b3.

        // save first two result words
        _mm512_store_epi64(c->data + 0 * VECLEN, acc_e0);
        _mm512_store_epi64(c->data + 1 * VECLEN, acc_e1);

        // propagate the carry through the rest of the lo result.
        b2 = _mm512_load_epi64(c->data + 2 * VECLEN);
        scarry = 0;
        a0 = _mm512_adc_epi52(b2, scarry, b3, &scarry);
        _mm512_store_epi64(c->data + 2 * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i = 3; (i < wshift) && (scarry > 0); i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "lo + hi * c:");
#endif

        // it is unlikely, but possible that we still have a final carry to propagate
        if (scarry > 0)
        {
            printf("unlikely add\n");
            b0 = _mm512_set1_epi64(mdata->isMersenne);
            a0 = _mm512_load_epi64(c->data + 0 * VECLEN);
            a0 = _mm512_adc_epi52(a0, 0, b0, &scarry);
            _mm512_store_epi64(c->data + 0 * VECLEN, a0);

            i = 1;
            while (scarry > 0)
            {
                a1 = _mm512_load_epi64(c->data + i * VECLEN);
                a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
                _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
                i++;
            }
        }

        // finally clear any hi-bits in the result
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif

        */
    }
    else if (mdata->isMersenne > 0)
    {
        // now add the low part into the high.
        scarry = 0;
        for (i = 0; i < wshift; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            b0 = _mm512_load_epi64(a->data + i * VECLEN);
            a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(a->data + i * VECLEN);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add it back in.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // clear the potential hi-bit
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
#endif

        mpz_t gn, gb, gc;
        mpz_init(gn);
        mpz_init(gb);
        mpz_init(gc);
        mpz_set_ui(gn, 1);
        mpz_mul_2exp(gn, gn, mdata->nbits);
        mpz_mod(gc, ga, gn);

        extract_bignum_from_vec_to_mpz(gb, c, 0, NWORDS);

        if (mpz_cmp(gc, gb) != 0)
        {
            gmp_printf("a = %Zd\n", ga);
            gmp_printf("n = %Zd\n", gn);

            gmp_printf("b = %Zd\n", gb);
            gmp_printf("c = %Zd\n", gc);
            exit(1);
        }

        
        mpz_clear(gb);
        mpz_clear(gc);
        mpz_clear(gn);

    }
    else
    {
        // now subtract the hi part from the lo
        scarry = 0;
        for (i = 0; i < wshift; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);   // hi
            b0 = _mm512_load_epi64(a->data + i * VECLEN);   // lo
            a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(a->data + i * VECLEN);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add 1.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // and add one to the hi-bit.  This should resolve the borrow bit.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        //a1 = _mm512_add_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift))), a1);
        //_mm512_store_epi64(c->data + wshift * VECLEN,
        //    _mm512_and_epi64(_mm512_set1_epi64(0xfffffffffffffULL), a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif

        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));
    }

    mpz_clear(ga);
    c->size = NWORDS;

}


#define and64 _mm512_and_epi64
#define storeu64 _mm512_store_epi64
#define add64 _mm512_add_epi64
#define set64 _mm512_set1_epi64
#define srli64 _mm512_srli_epi64
#define loadu64 _mm512_load_epi64
#define DIGIT_SIZE 52
#define DIGIT_MASK 0x000fffffffffffffULL
#define MB_WIDTH 8
#define SIMD_BYTES 64

__inline __m512i fma52lo(__m512i a, __m512i b, __m512i c)
{
    //return _mm512_madd52lo_epu64(a, b, c);
    __m512i t = _mm512_and_si512(_mm512_mullo_epi64(b, c), _mm512_set1_epi64(0x000fffffffffffffull));
    return _mm512_add_epi64(a, t);
}
__inline __m512i fma52hi(__m512i a, __m512i b, __m512i c, __m512d dbias, __m512i vbias1)
{
    //return _mm512_madd52hi_epu64(a, b, c);
    __m512d prod1_ld = _mm512_cvtepu64_pd(b);
    __m512d prod2_ld = _mm512_cvtepu64_pd(c);
    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
    return _mm512_add_epi64(a, _mm512_sub_epi64(castepu(prod1_ld), vbias1)); \
}

__mmask8 base_abssub_52(uint64_t* a, uint64_t* b, uint64_t* c, __mmask8 sm, int wordsa, int wordsb)
{
    // if sign-mask 'sm' is set, b > a, so we sub b-a instead of a-b.
    int i;
    __m512i avec, bvec, cvec, tmp;
    __mmask8 carry = 0;

    for (i = 0; i < MIN(wordsa, wordsb); i++)
    {
        avec = _mm512_load_epi64(a + i * VECLEN);
        bvec = _mm512_load_epi64(b + i * VECLEN);

        tmp = avec;
        tmp = _mm512_mask_mov_epi64(tmp, ~sm, bvec);

        avec = _mm512_mask_mov_epi64(avec, sm, bvec);
        bvec = _mm512_mask_mov_epi64(bvec, sm, tmp);

        cvec = _mm512_sbb_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c + i * VECLEN, cvec);
    }

    if (wordsa > wordsb)
    {
        avec = _mm512_load_epi64(a + i * VECLEN);
        bvec = _mm512_setzero_si512();

        tmp = avec;
        tmp = _mm512_mask_mov_epi64(tmp, ~sm, bvec);

        avec = _mm512_mask_mov_epi64(avec, sm, bvec);
        bvec = _mm512_mask_mov_epi64(bvec, sm, tmp);

        cvec = _mm512_sbb_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c + i * VECLEN, cvec);
    }
    else if (wordsb > wordsa)
    {
        avec = _mm512_setzero_si512();
        bvec = _mm512_load_epi64(b + i * VECLEN);

        tmp = avec;
        tmp = _mm512_mask_mov_epi64(tmp, ~sm, bvec);

        avec = _mm512_mask_mov_epi64(avec, sm, bvec);
        bvec = _mm512_mask_mov_epi64(bvec, sm, tmp);

        cvec = _mm512_sbb_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c + i * VECLEN, cvec);
    }

    return carry;
}

__mmask8 base_add_52(uint64_t* a, uint64_t* b, uint64_t* c, int wordsa, int wordsb)
{
    int i;
    __m512i avec, bvec, cvec;
    __mmask8 carry = 0;

    for (i = 0; i < MIN(wordsa, wordsb); i++)
    {
        avec = _mm512_load_epi64(a + i * VECLEN);
        bvec = _mm512_load_epi64(b + i * VECLEN);
        cvec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c + i * VECLEN, cvec);
    }

    if (wordsa > wordsb)
    {
        avec = _mm512_load_epi64(a + i * VECLEN);
        bvec = _mm512_setzero_si512();
        cvec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c + i * VECLEN, cvec);
    }
    else if (wordsb > wordsa)
    {
        avec = _mm512_setzero_si512();
        bvec = _mm512_load_epi64(b + i * VECLEN);
        cvec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c + i * VECLEN, cvec);
    }

    return carry;
}

uint32_t vec_gte2_52(uint64_t* u, uint64_t* v, int sz)
{
    // decide if each of the vec_bignums in vec 'u' is >=
    // the corresponding vec_bignum in vec 'v'.
    // return a mask of results.
    int i;
    __mmask8 mdecided = 0x00;
    __mmask8 mgte = 0;

    for (i = sz - 1; i >= 0; --i)
    {
        __m512i a = _mm512_load_epi64(u + i * VECLEN);
        __m512i b = _mm512_load_epi64(v + i * VECLEN);

        mgte |= _mm512_mask_cmp_epu64_mask(~mdecided, a, b, _MM_CMPINT_GT);
        mdecided |= mgte | _mm512_mask_cmp_epu64_mask(~mdecided, a, b, _MM_CMPINT_LT);

        if (mdecided == 0xff)
            break;
    }

    //equal if still undecided
    mgte |= ~mdecided;

    return (uint32_t)mgte;
}

void kcombine2(uint64_t* c, uint64_t* z1, uint64_t* s1, __mmask8 sm, int words, int hiwords)
{
    int halfwords = words / 2;
    int i;
    // combine all of the last addition steps and utilize the
    // reduced radix instead of carryflag emulation.
    // z1 += clo
    // z1 += chi
    // c += z1 * base^halfwords
    // becomes
    // c = c + clo * base^halfwords + chi * base^halfwords + neg(z1) * base^halfwords
    __m512i zero = _mm512_set1_epi64(0);
    __m512i vone = _mm512_set1_epi64(1);
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i v1, v2, v3, vc = zero;
    __mmask8 carry = 0;

    for (i = 0; i < 2 * hiwords; i++)
    {
        v1 = _mm512_load_epi64(c + i * VECLEN);
        v2 = _mm512_load_epi64(c + (words + i) * VECLEN);
        _mm512_store_epi64(s1 + i * VECLEN, _mm512_add_epi64(v1, v2));
    }

    for (; i < words; i++)
    {
        v1 = _mm512_load_epi64(c + i * VECLEN);
        _mm512_store_epi64(s1 + i * VECLEN, v1);
    }

    i = halfwords;
    v1 = _mm512_load_epi64(c + i * VECLEN);
    v2 = _mm512_load_epi64(s1 + (i - halfwords) * VECLEN);
    v3 = _mm512_load_epi64(z1 + (i - halfwords) * VECLEN);

    v1 = _mm512_add_epi64(v1, v2);
    v3 = _mm512_mask_xor_epi64(v3, sm, v3, vlmask);
    v1 = _mm512_mask_add_epi64(v1, sm, v1, vone);

    v3 = _mm512_add_epi64(v3, vc);
    v1 = _mm512_add_epi64(v1, v3);

    // carryprop
    vc = _mm512_srli_epi64(v1, 52);
    v1 = _mm512_and_epi64(vlmask, v1);

    _mm512_store_epi64(c + i * VECLEN, v1);

    for (i = halfwords + 1; i < halfwords + words; i++)
    {
        v1 = _mm512_load_epi64(c + i * VECLEN);
        v2 = _mm512_load_epi64(s1 + (i - halfwords) * VECLEN);
        v3 = _mm512_load_epi64(z1 + (i - halfwords) * VECLEN);

        v1 = _mm512_add_epi64(v1, v2);
        v3 = _mm512_mask_xor_epi64(v3, sm, v3, vlmask);

        v3 = _mm512_add_epi64(v3, vc);
        v1 = _mm512_add_epi64(v1, v3);

        // carryprop
        vc = _mm512_srli_epi64(v1, 52);
        v1 = _mm512_and_epi64(vlmask, v1);

        _mm512_store_epi64(c + i * VECLEN, v1);
    }

    // add in the final hi-word carry then propagate any
    // additional carries through the final words of the output.
    // z1 was negative, so we have to subtract a sign bit from
    // this word.
    v1 = _mm512_load_epi64(c + i * VECLEN);
    v1 = _mm512_addsetc_epi52(v1, vc, &carry);
    v1 = _mm512_mask_sub_epi64(v1, sm, v1, vone);
    _mm512_store_epi64(c + i * VECLEN, v1);
    i++;

    while (carry)
    {
        v1 = _mm512_load_epi64(c + i * VECLEN);
        v1 = _mm512_addcarry_epi52(v1, carry, &carry);
        _mm512_store_epi64(c + i * VECLEN, v1);
        i++;
    }

    return;
}

void kcombine2_lo(uint64_t* c, uint64_t* z1, uint64_t* s1, __mmask8 sm, int words, int hiwords)
{
    int halfwords = words / 2;
    int i;
    // combine all of the last addition steps and utilize the
    // reduced radix instead of carryflag emulation.
    // z1 += clo
    // z1 += chi
    // c += z1 * base^halfwords
    // becomes
    // c = c + clo * base^halfwords + chi * base^halfwords + neg(z1) * base^halfwords
    __m512i zero = _mm512_set1_epi64(0);
    __m512i vone = _mm512_set1_epi64(1);
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i v1, v2, v3, vc = zero;
    __mmask8 carry = 0;

    for (i = 0; i < halfwords; i++)
    {
        v1 = _mm512_load_epi64(c + i * VECLEN);
        v2 = _mm512_load_epi64(c + (words + i) * VECLEN);
        _mm512_store_epi64(s1 + i * VECLEN, _mm512_add_epi64(v1, v2));
    }

    i = halfwords;
    v1 = _mm512_load_epi64(c + i * VECLEN);
    v2 = _mm512_load_epi64(s1 + (i - halfwords) * VECLEN);
    v3 = _mm512_load_epi64(z1 + (i - halfwords) * VECLEN);

    v1 = _mm512_add_epi64(v1, v2);
    v3 = _mm512_mask_xor_epi64(v3, sm, v3, vlmask);
    v1 = _mm512_mask_add_epi64(v1, sm, v1, vone);

    v3 = _mm512_add_epi64(v3, vc);
    v1 = _mm512_add_epi64(v1, v3);

    // carryprop
    vc = _mm512_srli_epi64(v1, 52);
    v1 = _mm512_and_epi64(vlmask, v1);

    _mm512_store_epi64(c + i * VECLEN, v1);

    for (i = halfwords + 1; i < words; i++)
    {
        v1 = _mm512_load_epi64(c + i * VECLEN);
        v2 = _mm512_load_epi64(s1 + (i - halfwords) * VECLEN);
        v3 = _mm512_load_epi64(z1 + (i - halfwords) * VECLEN);

        v1 = _mm512_add_epi64(v1, v2);
        v3 = _mm512_mask_xor_epi64(v3, sm, v3, vlmask);

        v3 = _mm512_add_epi64(v3, vc);
        v1 = _mm512_add_epi64(v1, v3);

        // carryprop
        vc = _mm512_srli_epi64(v1, 52);
        v1 = _mm512_and_epi64(vlmask, v1);

        _mm512_store_epi64(c + i * VECLEN, v1);
    }

    return;
}

#ifdef IFMA
#define VEC_MUL4_ACCUMX(x, b0, b1, b2, b3, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4) \
    lo1 = _mm512_madd52lo_epu64(lo1, x, b0); \
    lo2 = _mm512_madd52lo_epu64(lo2, x, b1); \
    lo3 = _mm512_madd52lo_epu64(lo3, x, b2); \
    lo4 = _mm512_madd52lo_epu64(lo4, x, b3); \
    hi1 = _mm512_madd52hi_epu64(hi1, x, b0); \
    hi2 = _mm512_madd52hi_epu64(hi2, x, b1); \
    hi3 = _mm512_madd52hi_epu64(hi3, x, b2); \
    hi4 = _mm512_madd52hi_epu64(hi4, x, b3);

#define SUB_BIAS_HI4(bias1, bias2, bias3, bias4) {}
#define SUB_BIAS_LO4(bias1, bias2, bias3, bias4) {}
#define SUB_BIAS_HI3(bias1, bias2, bias3) {}
#define SUB_BIAS_LO3(bias1, bias2, bias3) {}

#else
#define VEC_MUL4_ACCUMX(x, b0, b1, b2, b3, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4) \
    prod5_ld = _mm512_cvtepu64_pd(x);                                                                           \
    prod1_ld = _mm512_cvtepu64_pd(b0);                                                                          \
    prod2_ld = _mm512_cvtepu64_pd(b1);                                                                          \
    prod3_ld = _mm512_cvtepu64_pd(b2);                                                                          \
    prod4_ld = _mm512_cvtepu64_pd(b3);                                                                          \
    prod1_hd = _mm512_fmadd_round_pd(prod5_ld, prod1_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));      \
    prod2_hd = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));      \
    prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));      \
    prod4_hd = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));      \
    hi1 = _mm512_add_epi64(hi1, _mm512_castpd_si512(prod1_hd));                                                 \
    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);                                            \
    hi2 = _mm512_add_epi64(hi2, _mm512_castpd_si512(prod2_hd));                                                 \
    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);                                            \
    hi3 = _mm512_add_epi64(hi3, _mm512_castpd_si512(prod3_hd));                                                 \
    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);                                            \
    hi4 = _mm512_add_epi64(hi4, _mm512_castpd_si512(prod4_hd));                                                 \
    prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);                                            \
    prod1_ld = _mm512_fmadd_round_pd(prod5_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
    prod2_ld = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
    prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
    prod4_ld = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
    lo1 = _mm512_add_epi64(lo1, _mm512_castpd_si512(prod1_ld));                                                 \
    lo2 = _mm512_add_epi64(lo2, _mm512_castpd_si512(prod2_ld));                                                 \
    lo3 = _mm512_add_epi64(lo3, _mm512_castpd_si512(prod3_ld));                                                 \
    lo4 = _mm512_add_epi64(lo4, _mm512_castpd_si512(prod4_ld));


#define SUB_BIAS_HI4(bias1, bias2, bias3, bias4) \
    A0 = _mm512_set1_epi64((bias1) * 0x467);         \
    A1 = _mm512_set1_epi64((bias2) * 0x467);         \
    A2 = _mm512_set1_epi64((bias3) * 0x467);         \
    A3 = _mm512_set1_epi64((bias4) * 0x467);         \
    A0 = _mm512_slli_epi64(A0, 52);            \
    A1 = _mm512_slli_epi64(A1, 52);            \
    A2 = _mm512_slli_epi64(A2, 52);            \
    A3 = _mm512_slli_epi64(A3, 52);            \
    hi0 = _mm512_sub_epi64(hi0, A0); \
    hi1 = _mm512_sub_epi64(hi1, A1); \
    hi2 = _mm512_sub_epi64(hi2, A2); \
    hi3 = _mm512_sub_epi64(hi3, A3);


#define SUB_BIAS_LO4(bias1, bias2, bias3, bias4) \
    A0 = _mm512_set1_epi64((bias1) * 0x433);         \
    A1 = _mm512_set1_epi64((bias2) * 0x433);         \
    A2 = _mm512_set1_epi64((bias3) * 0x433);         \
    A3 = _mm512_set1_epi64((bias4) * 0x433);         \
    A0 = _mm512_slli_epi64(A0, 52);            \
    A1 = _mm512_slli_epi64(A1, 52);            \
    A2 = _mm512_slli_epi64(A2, 52);            \
    A3 = _mm512_slli_epi64(A3, 52);            \
    lo0 = _mm512_sub_epi64(lo0, A0); \
    lo1 = _mm512_sub_epi64(lo1, A1); \
    lo2 = _mm512_sub_epi64(lo2, A2); \
    lo3 = _mm512_sub_epi64(lo3, A3);

#define SUB_BIAS_HI3(bias1, bias2, bias3) \
    A0 = _mm512_set1_epi64((bias1) * 0x467);         \
    A1 = _mm512_set1_epi64((bias2) * 0x467);         \
    A2 = _mm512_set1_epi64((bias3) * 0x467);         \
    A0 = _mm512_slli_epi64(A0, 52);            \
    A1 = _mm512_slli_epi64(A1, 52);            \
    A2 = _mm512_slli_epi64(A2, 52);            \
    hi4 = _mm512_sub_epi64(hi4, A0); \
    hi5 = _mm512_sub_epi64(hi5, A1); \
    hi6 = _mm512_sub_epi64(hi6, A2);


#define SUB_BIAS_LO3(bias1, bias2, bias3) \
    A0 = _mm512_set1_epi64((bias1) * 0x433);         \
    A1 = _mm512_set1_epi64((bias2) * 0x433);         \
    A2 = _mm512_set1_epi64((bias3) * 0x433);         \
    A0 = _mm512_slli_epi64(A0, 52);            \
    A1 = _mm512_slli_epi64(A1, 52);            \
    A2 = _mm512_slli_epi64(A2, 52);            \
    lo4 = _mm512_sub_epi64(lo4, A0); \
    lo5 = _mm512_sub_epi64(lo5, A1); \
    lo6 = _mm512_sub_epi64(lo6, A2);
#endif


//#define DEBUG_VECMUL

void vecmul52_cios(uint64_t* a, uint64_t* b, uint64_t* c, int n)
{
    int i, j, k;

    __m512i lo0, lo1, lo2, lo3, lo4, lo5, lo6;
    __m512i hi0, hi1, hi2, hi3, hi4, hi5, hi6;
    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);
    __m512i zero = set64(0), B, A, C, plo, carry;
    __m512i MASK = set64(DIGIT_MASK);

    int full_blocks = n / BLOCKWORDS;

#ifdef DEBUG_VECMUL
    printf("commencing vecmul52_cios with WORDS=%d and NBLOCKS=%d\n", n, full_blocks);
    print_vechex(a, 0, n, "input a: ");
    print_vechex(b, 0, n, "input b: ");
    print_vechex(c, 0, 2 * n, "output c: ");
#endif

    for (i = 0; i < 2 * n; i++) {
        _mm512_store_epi64(c + i * VECLEN, _mm512_setzero_si512());
    }
    
    for (j = 0; j < full_blocks; j++) {
        __m512i B0 = loadu64(b + (j * BLOCKWORDS + 0) * VECLEN);
        __m512i B1 = loadu64(b + (j * BLOCKWORDS + 1) * VECLEN);
        __m512i B2 = loadu64(b + (j * BLOCKWORDS + 2) * VECLEN);
        __m512i B3 = loadu64(b + (j * BLOCKWORDS + 3) * VECLEN);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            __m512i A0 = loadu64(a + (i * BLOCKWORDS + 0) * VECLEN);
            __m512i A1 = loadu64(a + (i * BLOCKWORDS + 1) * VECLEN);
            __m512i A2 = loadu64(a + (i * BLOCKWORDS + 2) * VECLEN);
            __m512i A3 = loadu64(a + (i * BLOCKWORDS + 3) * VECLEN);

            VEC_MUL4_ACCUMX(A0, B0, B1, B2, B3, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, B1, B2, B3, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, B1, B2, B3, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, B1, B2, B3, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

            A0 = loadu64(c + ((i + j) * BLOCKWORDS + 0) * VECLEN);
            A1 = loadu64(c + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            A2 = loadu64(c + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            A3 = loadu64(c + ((i + j) * BLOCKWORDS + 3) * VECLEN);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS + 0) * VECLEN, 
                _mm512_and_si512(lo0, MASK));
            
            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS + 1) * VECLEN,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS + 2) * VECLEN,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS + 3) * VECLEN,
                _mm512_and_si512(lo3, MASK));

            hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

            // subtract out all of the bias at once.
            SUB_BIAS_HI3(3, 2, 1);
            SUB_BIAS_LO3(3, 2, 1);

            A0 = loadu64(c + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            A1 = loadu64(c + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            A2 = loadu64(c + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            A3 = loadu64(c + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS + 4) * VECLEN,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS + 5) * VECLEN,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS + 6) * VECLEN,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS + 7) * VECLEN,
                _mm512_and_si512(hi6, MASK));

            if (((i + j) * BLOCKWORDS + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i + j) * BLOCKWORDS + 8) * VECLEN, 
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i + j) * BLOCKWORDS + 8) * VECLEN), 
                        _mm512_srli_epi64(hi6, 52)));
            }

#ifdef DEBUG_VECMUL
            print_vechex(c, 0, 2 * n, "end of block iteration:");
#endif
        }
    }

    if (n - (full_blocks * BLOCKWORDS) == 3)
    {
        j = full_blocks * BLOCKWORDS;
        __m512i B0 = loadu64(b + (j + 0) * VECLEN);
        __m512i B1 = loadu64(b + (j + 1) * VECLEN);
        __m512i B2 = loadu64(b + (j + 2) * VECLEN);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            __m512i A0 = loadu64(a + (i * BLOCKWORDS + 0) * VECLEN);
            __m512i A1 = loadu64(a + (i * BLOCKWORDS + 1) * VECLEN);
            __m512i A2 = loadu64(a + (i * BLOCKWORDS + 2) * VECLEN);
            __m512i A3 = loadu64(a + (i * BLOCKWORDS + 3) * VECLEN);

            VEC_MUL4_ACCUMX(A0, B0, B1, B2, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, B1, B2, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, B1, B2, zero, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, B1, B2, zero, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

            A0 = loadu64(c + ((i) * BLOCKWORDS + j + 0) * VECLEN);
            A1 = loadu64(c + ((i) * BLOCKWORDS + j + 1) * VECLEN);
            A2 = loadu64(c + ((i) * BLOCKWORDS + j + 2) * VECLEN);
            A3 = loadu64(c + ((i) * BLOCKWORDS + j + 3) * VECLEN);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 0) * VECLEN,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 1) * VECLEN,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 2) * VECLEN,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 3) * VECLEN,
                _mm512_and_si512(lo3, MASK));

            hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

            // subtract out all of the bias at once.
            SUB_BIAS_HI3(3, 2, 1);
            SUB_BIAS_LO3(3, 2, 1);

            A0 = loadu64(c + ((i) * BLOCKWORDS + j + 4) * VECLEN);
            A1 = loadu64(c + ((i) * BLOCKWORDS + j + 5) * VECLEN);
            A2 = loadu64(c + ((i) * BLOCKWORDS + j + 6) * VECLEN);
            A3 = loadu64(c + ((i) * BLOCKWORDS + j + 7) * VECLEN);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 4) * VECLEN,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 5) * VECLEN,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 6) * VECLEN,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 7) * VECLEN,
                _mm512_and_si512(hi6, MASK));

            if (((i)*BLOCKWORDS + j + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 8) * VECLEN,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i)*BLOCKWORDS + j + 8) * VECLEN),
                        _mm512_srli_epi64(hi6, 52)));
            }
        }

        j = full_blocks * BLOCKWORDS;
        B0 = loadu64(a + (j + 0) * VECLEN);
        B1 = loadu64(a + (j + 1) * VECLEN);
        B2 = loadu64(a + (j + 2) * VECLEN);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            __m512i A0 = loadu64(b + (i * BLOCKWORDS + 0) * VECLEN);
            __m512i A1 = loadu64(b + (i * BLOCKWORDS + 1) * VECLEN);
            __m512i A2 = loadu64(b + (i * BLOCKWORDS + 2) * VECLEN);
            __m512i A3 = loadu64(b + (i * BLOCKWORDS + 3) * VECLEN);

            VEC_MUL4_ACCUMX(A0, B0, B1, B2, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, B1, B2, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, B1, B2, zero, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, B1, B2, zero, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

            A0 = loadu64(c + ((i)*BLOCKWORDS + j + 0) * VECLEN);
            A1 = loadu64(c + ((i)*BLOCKWORDS + j + 1) * VECLEN);
            A2 = loadu64(c + ((i)*BLOCKWORDS + j + 2) * VECLEN);
            A3 = loadu64(c + ((i)*BLOCKWORDS + j + 3) * VECLEN);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 0) * VECLEN,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 1) * VECLEN,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 2) * VECLEN,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 3) * VECLEN,
                _mm512_and_si512(lo3, MASK));

            hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

            // subtract out all of the bias at once.
            SUB_BIAS_HI3(3, 2, 1);
            SUB_BIAS_LO3(3, 2, 1);

            A0 = loadu64(c + ((i)*BLOCKWORDS + j + 4) * VECLEN);
            A1 = loadu64(c + ((i)*BLOCKWORDS + j + 5) * VECLEN);
            A2 = loadu64(c + ((i)*BLOCKWORDS + j + 6) * VECLEN);
            A3 = loadu64(c + ((i)*BLOCKWORDS + j + 7) * VECLEN);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 4) * VECLEN,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 5) * VECLEN,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 6) * VECLEN,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 7) * VECLEN,
                _mm512_and_si512(hi6, MASK));

            if (((i)*BLOCKWORDS + j + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 8) * VECLEN,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i)*BLOCKWORDS + j + 8) * VECLEN),
                        _mm512_srli_epi64(hi6, 52)));
            }
        }

#ifdef DEBUG_VECMUL
        print_vechex(c, 0, 2 * n, "end of last partial block:");
#endif

        j = full_blocks * BLOCKWORDS;
        i = full_blocks * BLOCKWORDS;
        lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
        hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

        B0 = loadu64(b + (j + 0) * VECLEN);
        B1 = loadu64(b + (j + 1) * VECLEN);
        B2 = loadu64(b + (j + 2) * VECLEN);
        __m512i A0 = loadu64(a + (i + 0) * VECLEN);
        __m512i A1 = loadu64(a + (i + 1) * VECLEN);
        __m512i A2 = loadu64(a + (i + 2) * VECLEN);
        __m512i A3 = zero;

        VEC_MUL4_ACCUMX(A0, B0, B1, B2, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
        VEC_MUL4_ACCUMX(A1, B0, B1, B2, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
        VEC_MUL4_ACCUMX(A2, B0, B1, B2, zero, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);

        // subtract out all of the bias at once.
        SUB_BIAS_HI4(1, 2, 3, 3);
        SUB_BIAS_LO4(1, 2, 3, 3);

        A0 = loadu64(c + ((i + j) + 0) * VECLEN);
        A1 = loadu64(c + ((i + j) + 1) * VECLEN);
        A2 = loadu64(c + ((i + j) + 2) * VECLEN);
        A3 = loadu64(c + ((i + j) + 3) * VECLEN);

        lo0 = _mm512_add_epi64(lo0, A0);
        lo1 = _mm512_add_epi64(lo1, A1);
        lo2 = _mm512_add_epi64(lo2, A2);
        lo3 = _mm512_add_epi64(lo3, A3);

        lo1 = _mm512_add_epi64(lo1, hi0);
        lo2 = _mm512_add_epi64(lo2, hi1);
        lo3 = _mm512_add_epi64(lo3, hi2);

        _mm512_store_epi64(c + ((i + j) + 0) * VECLEN,
            _mm512_and_si512(lo0, MASK));

        hi0 = _mm512_srli_epi64(lo0, 52);
        lo1 = _mm512_add_epi64(lo1, hi0);

        _mm512_store_epi64(c + ((i + j) + 1) * VECLEN,
            _mm512_and_si512(lo1, MASK));

        hi1 = _mm512_srli_epi64(lo1, 52);
        lo2 = _mm512_add_epi64(lo2, hi1);

        _mm512_store_epi64(c + ((i + j) + 2) * VECLEN,
            _mm512_and_si512(lo2, MASK));

        hi2 = _mm512_srli_epi64(lo2, 52);
        lo3 = _mm512_add_epi64(lo3, hi2);

        _mm512_store_epi64(c + ((i + j) + 3) * VECLEN,
            _mm512_and_si512(lo3, MASK));

        hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

        // subtract out all of the bias at once.
        SUB_BIAS_HI3(2, 1, 0);
        SUB_BIAS_LO3(2, 1, 0);

        A0 = loadu64(c + ((i + j) + 4) * VECLEN);
        A1 = loadu64(c + ((i + j) + 5) * VECLEN);

        lo4 = _mm512_add_epi64(lo4, A0);
        lo5 = _mm512_add_epi64(lo5, A1);

        lo4 = _mm512_add_epi64(lo4, hi3);
        lo5 = _mm512_add_epi64(lo5, hi4);

        _mm512_store_epi64(c + ((i + j) + 4) * VECLEN,
            _mm512_and_si512(lo4, MASK));

        hi4 = _mm512_srli_epi64(lo4, 52);
        lo5 = _mm512_add_epi64(lo5, hi4);

        _mm512_store_epi64(c + ((i + j) + 5) * VECLEN,
            _mm512_and_si512(lo5, MASK));
    }

    if (n - (full_blocks * BLOCKWORDS) == 2)
    {
        j = full_blocks * BLOCKWORDS;
        __m512i B0 = loadu64(b + (j + 0) * VECLEN);
        __m512i B1 = loadu64(b + (j + 1) * VECLEN);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            __m512i A0 = loadu64(a + (i * BLOCKWORDS + 0) * VECLEN);
            __m512i A1 = loadu64(a + (i * BLOCKWORDS + 1) * VECLEN);
            __m512i A2 = loadu64(a + (i * BLOCKWORDS + 2) * VECLEN);
            __m512i A3 = loadu64(a + (i * BLOCKWORDS + 3) * VECLEN);

            VEC_MUL4_ACCUMX(A0, B0, B1, zero, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, B1, zero, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, B1, zero, zero, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, B1, zero, zero, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

#ifdef DEBUG_VECMUL
            print_regvechex(lo0, 0, "lo0: ");
            print_regvechex(lo1, 0, "lo1: ");
            print_regvechex(lo2, 0, "lo2: ");
            print_regvechex(lo3, 0, "lo3: ");

            print_regvechex(hi0, 0, "hi0: ");
            print_regvechex(hi1, 0, "hi1: ");
            print_regvechex(hi2, 0, "hi2: ");
            print_regvechex(hi3, 0, "hi3: ");
#endif

            A0 = loadu64(c + ((i) * BLOCKWORDS + j + 0) * VECLEN);
            A1 = loadu64(c + ((i) * BLOCKWORDS + j + 1) * VECLEN);
            A2 = loadu64(c + ((i) * BLOCKWORDS + j + 2) * VECLEN);
            A3 = loadu64(c + ((i) * BLOCKWORDS + j + 3) * VECLEN);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 0) * VECLEN,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 1) * VECLEN,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 2) * VECLEN,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 3) * VECLEN,
                _mm512_and_si512(lo3, MASK));

            hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

            // subtract out all of the bias at once.
            SUB_BIAS_HI3(3, 2, 1);
            SUB_BIAS_LO3(3, 2, 1);

#ifdef DEBUG_VECMUL
            print_regvechex(lo4, 0, "lo4: ");
            print_regvechex(lo5, 0, "lo5: ");
            print_regvechex(lo6, 0, "lo6: ");
            print_regvechex(hi4, 0, "hi4: ");
            print_regvechex(hi5, 0, "hi5: ");
            print_regvechex(hi6, 0, "hi6: ");
#endif

            A0 = loadu64(c + ((i) * BLOCKWORDS + j + 4) * VECLEN);
            A1 = loadu64(c + ((i) * BLOCKWORDS + j + 5) * VECLEN);
            A2 = loadu64(c + ((i) * BLOCKWORDS + j + 6) * VECLEN);
            A3 = loadu64(c + ((i) * BLOCKWORDS + j + 7) * VECLEN);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 4) * VECLEN,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 5) * VECLEN,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 6) * VECLEN,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 7) * VECLEN,
                _mm512_and_si512(hi6, MASK));

            if (((i)*BLOCKWORDS + j + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 8) * VECLEN,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i)*BLOCKWORDS + j + 8) * VECLEN),
                        _mm512_srli_epi64(hi6, 52)));
            }
        }

        B0 = loadu64(a + (j + 0) * VECLEN);
        B1 = loadu64(a + (j + 1) * VECLEN);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            __m512i A0 = loadu64(b + (i * BLOCKWORDS + 0) * VECLEN);
            __m512i A1 = loadu64(b + (i * BLOCKWORDS + 1) * VECLEN);
            __m512i A2 = loadu64(b + (i * BLOCKWORDS + 2) * VECLEN);
            __m512i A3 = loadu64(b + (i * BLOCKWORDS + 3) * VECLEN);

            VEC_MUL4_ACCUMX(A0, B0, B1, zero, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, B1, zero, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, B1, zero, zero, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, B1, zero, zero, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

            A0 = loadu64(c + ((i)*BLOCKWORDS + j + 0) * VECLEN);
            A1 = loadu64(c + ((i)*BLOCKWORDS + j + 1) * VECLEN);
            A2 = loadu64(c + ((i)*BLOCKWORDS + j + 2) * VECLEN);
            A3 = loadu64(c + ((i)*BLOCKWORDS + j + 3) * VECLEN);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 0) * VECLEN,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 1) * VECLEN,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 2) * VECLEN,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 3) * VECLEN,
                _mm512_and_si512(lo3, MASK));

            hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

            // subtract out all of the bias at once.
            SUB_BIAS_HI3(3, 2, 1);
            SUB_BIAS_LO3(3, 2, 1);

            A0 = loadu64(c + ((i)*BLOCKWORDS + j + 4) * VECLEN);
            A1 = loadu64(c + ((i)*BLOCKWORDS + j + 5) * VECLEN);
            A2 = loadu64(c + ((i)*BLOCKWORDS + j + 6) * VECLEN);
            A3 = loadu64(c + ((i)*BLOCKWORDS + j + 7) * VECLEN);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 4) * VECLEN,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 5) * VECLEN,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 6) * VECLEN,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 7) * VECLEN,
                _mm512_and_si512(hi6, MASK));

            if (((i)*BLOCKWORDS + j + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 8) * VECLEN,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i)*BLOCKWORDS + j + 8) * VECLEN),
                        _mm512_srli_epi64(hi6, 52)));
            }
        }

        j = full_blocks * BLOCKWORDS;
        i = full_blocks * BLOCKWORDS;
        lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
        hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

        B0 = loadu64(b + (j + 0) * VECLEN);
        B1 = loadu64(b + (j + 1) * VECLEN);
        __m512i A0 = loadu64(a + (i + 0) * VECLEN);
        __m512i A1 = loadu64(a + (i + 1) * VECLEN);
        __m512i A2 = zero;
        __m512i A3 = zero;

        VEC_MUL4_ACCUMX(A0, B0, B1, zero, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
        VEC_MUL4_ACCUMX(A1, B0, B1, zero, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);

        // subtract out all of the bias at once.
        SUB_BIAS_HI4(1, 2, 2, 2);
        SUB_BIAS_LO4(1, 2, 2, 2);

        A0 = loadu64(c + ((i + j) + 0) * VECLEN);
        A1 = loadu64(c + ((i + j) + 1) * VECLEN);
        A2 = loadu64(c + ((i + j) + 2) * VECLEN);
        A3 = loadu64(c + ((i + j) + 3) * VECLEN);

        lo0 = _mm512_add_epi64(lo0, A0);
        lo1 = _mm512_add_epi64(lo1, A1);
        lo2 = _mm512_add_epi64(lo2, A2);
        lo3 = _mm512_add_epi64(lo3, A3);

        lo1 = _mm512_add_epi64(lo1, hi0);
        lo2 = _mm512_add_epi64(lo2, hi1);
        lo3 = _mm512_add_epi64(lo3, hi2);

        _mm512_store_epi64(c + ((i + j) + 0) * VECLEN,
            _mm512_and_si512(lo0, MASK));

        hi0 = _mm512_srli_epi64(lo0, 52);
        lo1 = _mm512_add_epi64(lo1, hi0);

        _mm512_store_epi64(c + ((i + j) + 1) * VECLEN,
            _mm512_and_si512(lo1, MASK));

        hi1 = _mm512_srli_epi64(lo1, 52);
        lo2 = _mm512_add_epi64(lo2, hi1);

        _mm512_store_epi64(c + ((i + j) + 2) * VECLEN,
            _mm512_and_si512(lo2, MASK));

        hi2 = _mm512_srli_epi64(lo2, 52);
        lo3 = _mm512_add_epi64(lo3, hi2);

        _mm512_store_epi64(c + ((i + j) + 3) * VECLEN,
            _mm512_and_si512(lo3, MASK));
    }

    if (n - (full_blocks * BLOCKWORDS) == 1)
    {
        j = full_blocks * BLOCKWORDS;
        __m512i B0 = loadu64(b + (j + 0) * VECLEN);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            __m512i A0 = loadu64(a + (i * BLOCKWORDS + 0) * VECLEN);
            __m512i A1 = loadu64(a + (i * BLOCKWORDS + 1) * VECLEN);
            __m512i A2 = loadu64(a + (i * BLOCKWORDS + 2) * VECLEN);
            __m512i A3 = loadu64(a + (i * BLOCKWORDS + 3) * VECLEN);

            VEC_MUL4_ACCUMX(A0, B0, zero, zero, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, zero, zero, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, zero, zero, zero, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, zero, zero, zero, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

            A0 = loadu64(c + ((i)*BLOCKWORDS + j + 0) * VECLEN);
            A1 = loadu64(c + ((i)*BLOCKWORDS + j + 1) * VECLEN);
            A2 = loadu64(c + ((i)*BLOCKWORDS + j + 2) * VECLEN);
            A3 = loadu64(c + ((i)*BLOCKWORDS + j + 3) * VECLEN);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 0) * VECLEN,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 1) * VECLEN,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 2) * VECLEN,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 3) * VECLEN,
                _mm512_and_si512(lo3, MASK));

            hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

            // subtract out all of the bias at once.
            SUB_BIAS_HI3(3, 2, 1);
            SUB_BIAS_LO3(3, 2, 1);

#ifdef DEBUG_VECMUL
            print_regvechex(lo4, 0, "lo4: ");
            print_regvechex(lo5, 0, "lo5: ");
            print_regvechex(lo6, 0, "lo6: ");
            print_regvechex(hi4, 0, "hi4: ");
            print_regvechex(hi5, 0, "hi5: ");
            print_regvechex(hi6, 0, "hi6: ");
#endif

            A0 = loadu64(c + ((i)*BLOCKWORDS + j + 4) * VECLEN);
            A1 = loadu64(c + ((i)*BLOCKWORDS + j + 5) * VECLEN);
            A2 = loadu64(c + ((i)*BLOCKWORDS + j + 6) * VECLEN);
            A3 = loadu64(c + ((i)*BLOCKWORDS + j + 7) * VECLEN);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 4) * VECLEN,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 5) * VECLEN,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 6) * VECLEN,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 7) * VECLEN,
                _mm512_and_si512(hi6, MASK));

            if (((i)*BLOCKWORDS + j + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 8) * VECLEN,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i)*BLOCKWORDS + j + 8) * VECLEN),
                        _mm512_srli_epi64(hi6, 52)));
            }
        }

        B0 = loadu64(a + (j + 0) * VECLEN);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            __m512i A0 = loadu64(b + (i * BLOCKWORDS + 0) * VECLEN);
            __m512i A1 = loadu64(b + (i * BLOCKWORDS + 1) * VECLEN);
            __m512i A2 = loadu64(b + (i * BLOCKWORDS + 2) * VECLEN);
            __m512i A3 = loadu64(b + (i * BLOCKWORDS + 3) * VECLEN);

            VEC_MUL4_ACCUMX(A0, B0, zero, zero, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, zero, zero, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, zero, zero, zero, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, zero, zero, zero, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

            A0 = loadu64(c + ((i)*BLOCKWORDS + j + 0) * VECLEN);
            A1 = loadu64(c + ((i)*BLOCKWORDS + j + 1) * VECLEN);
            A2 = loadu64(c + ((i)*BLOCKWORDS + j + 2) * VECLEN);
            A3 = loadu64(c + ((i)*BLOCKWORDS + j + 3) * VECLEN);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 0) * VECLEN,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 1) * VECLEN,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 2) * VECLEN,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 3) * VECLEN,
                _mm512_and_si512(lo3, MASK));

            hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

            // subtract out all of the bias at once.
            SUB_BIAS_HI3(3, 2, 1);
            SUB_BIAS_LO3(3, 2, 1);

            A0 = loadu64(c + ((i)*BLOCKWORDS + j + 4) * VECLEN);
            A1 = loadu64(c + ((i)*BLOCKWORDS + j + 5) * VECLEN);
            A2 = loadu64(c + ((i)*BLOCKWORDS + j + 6) * VECLEN);
            A3 = loadu64(c + ((i)*BLOCKWORDS + j + 7) * VECLEN);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 4) * VECLEN,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 5) * VECLEN,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 6) * VECLEN,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 7) * VECLEN,
                _mm512_and_si512(hi6, MASK));

            if (((i)*BLOCKWORDS + j + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 8) * VECLEN,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i)*BLOCKWORDS + j + 8) * VECLEN),
                        _mm512_srli_epi64(hi6, 52)));
            }
        }

        j = full_blocks * BLOCKWORDS;
        i = full_blocks * BLOCKWORDS;
        lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
        hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

        B0 = loadu64(b + (j + 0) * VECLEN);
        __m512i A0 = loadu64(a + (i + 0) * VECLEN);
        __m512i A1 = zero;
        __m512i A2 = zero;
        __m512i A3 = zero;

        VEC_MUL4_ACCUMX(A0, B0, zero, zero, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);

        // subtract out all of the bias at once.
        SUB_BIAS_HI4(1, 1, 1, 1);
        SUB_BIAS_LO4(1, 1, 1, 1);

        A0 = loadu64(c + ((i + j) + 0) * VECLEN);
        A1 = loadu64(c + ((i + j) + 1) * VECLEN);

        lo0 = _mm512_add_epi64(lo0, A0);
        lo1 = _mm512_add_epi64(lo1, A1);

        lo1 = _mm512_add_epi64(lo1, hi0);

        _mm512_store_epi64(c + ((i + j) + 0) * VECLEN,
            _mm512_and_si512(lo0, MASK));

        hi0 = _mm512_srli_epi64(lo0, 52);
        lo1 = _mm512_add_epi64(lo1, hi0);

        _mm512_store_epi64(c + ((i + j) + 1) * VECLEN,
            _mm512_and_si512(lo1, MASK));
    }

    if (0) {
        for (j = full_blocks * BLOCKWORDS; j < n; j++) {
            B = loadu64(b + j * VECLEN);

            carry = zero;
            for (i = 0; i < n; i++) {
                A = loadu64(a + i * VECLEN);
                C = loadu64(c + (i + j) * VECLEN);

                plo = fma52lo(C, A, B);
                plo = _mm512_add_epi64(plo, carry);
                carry = fma52hi(_mm512_srli_epi64(plo, 52), A, B, dbias, vbias1);

                _mm512_store_epi64(c + (i + j) * VECLEN, _mm512_and_si512(plo, MASK));
            }
            _mm512_store_epi64(c + (i + j) * VECLEN, _mm512_add_epi64(
                loadu64(c + (i + j) * VECLEN), _mm512_and_si512(carry, MASK)));
        }
    }

#ifdef DEBUG_VECMUL
    print_vechex(c, 0, 2 * n, "end of multiplication:");
#endif

#ifdef DEBUG_VECMUL
    exit(1);
#endif
    return;
}

void vecmul52_cios_lo(uint64_t* a, uint64_t* b, uint64_t* c, int n)
{
    int i, j, k;

    __m512i lo0, lo1, lo2, lo3, lo4, lo5, lo6;
    __m512i hi0, hi1, hi2, hi3, hi4, hi5, hi6;
    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);
    __m512i zero = set64(0), B, A, C, plo, carry;
    __m512i MASK = set64(DIGIT_MASK);

    int full_blocks = n / BLOCKWORDS;

#ifdef DEBUG_VECMUL
    printf("commencing vecmul52_cios with WORDS=%d and NBLOCKS=%d\n", n, full_blocks);
    print_vechex(a, 0, n, "input a: ");
    print_vechex(b, 0, n, "input b: ");
    print_vechex(c, 0, 2 * n, "output c: ");
#endif

    for (i = 0; i < 2 * n; i++) {
        _mm512_store_epi64(c + i * VECLEN, _mm512_setzero_si512());
    }

    for (j = 0; j < full_blocks; j++) {
        __m512i B0 = loadu64(b + (j * BLOCKWORDS + 0) * VECLEN);
        __m512i B1 = loadu64(b + (j * BLOCKWORDS + 1) * VECLEN);
        __m512i B2 = loadu64(b + (j * BLOCKWORDS + 2) * VECLEN);
        __m512i B3 = loadu64(b + (j * BLOCKWORDS + 3) * VECLEN);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            if ((i + j) * BLOCKWORDS > n)
                continue;

            __m512i A0 = loadu64(a + (i * BLOCKWORDS + 0) * VECLEN);
            __m512i A1 = loadu64(a + (i * BLOCKWORDS + 1) * VECLEN);
            __m512i A2 = loadu64(a + (i * BLOCKWORDS + 2) * VECLEN);
            __m512i A3 = loadu64(a + (i * BLOCKWORDS + 3) * VECLEN);

            VEC_MUL4_ACCUMX(A0, B0, B1, B2, B3, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, B1, B2, B3, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, B1, B2, B3, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, B1, B2, B3, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

            A0 = loadu64(c + ((i + j) * BLOCKWORDS + 0) * VECLEN);
            A1 = loadu64(c + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            A2 = loadu64(c + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            A3 = loadu64(c + ((i + j) * BLOCKWORDS + 3) * VECLEN);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS + 0) * VECLEN,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS + 1) * VECLEN,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS + 2) * VECLEN,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS + 3) * VECLEN,
                _mm512_and_si512(lo3, MASK));

            hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

            // subtract out all of the bias at once.
            SUB_BIAS_HI3(3, 2, 1);
            SUB_BIAS_LO3(3, 2, 1);

            A0 = loadu64(c + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            A1 = loadu64(c + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            A2 = loadu64(c + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            A3 = loadu64(c + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS + 4) * VECLEN,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS + 5) * VECLEN,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS + 6) * VECLEN,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS + 7) * VECLEN,
                _mm512_and_si512(hi6, MASK));

            if (((i + j) * BLOCKWORDS + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i + j) * BLOCKWORDS + 8) * VECLEN,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i + j) * BLOCKWORDS + 8) * VECLEN),
                        _mm512_srli_epi64(hi6, 52)));
            }

#ifdef DEBUG_VECMUL
            print_vechex(c, 0, 2 * n, "end of block iteration:");
#endif
        }
    }

    if (n - (full_blocks * BLOCKWORDS) == 3)
    {
        j = full_blocks * BLOCKWORDS;
        __m512i B0 = loadu64(b + (j + 0) * VECLEN);
        __m512i B1 = loadu64(b + (j + 1) * VECLEN);
        __m512i B2 = loadu64(b + (j + 2) * VECLEN);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            if ((i * BLOCKWORDS + j) > n)
                continue;

            __m512i A0 = loadu64(a + (i * BLOCKWORDS + 0) * VECLEN);
            __m512i A1 = loadu64(a + (i * BLOCKWORDS + 1) * VECLEN);
            __m512i A2 = loadu64(a + (i * BLOCKWORDS + 2) * VECLEN);
            __m512i A3 = loadu64(a + (i * BLOCKWORDS + 3) * VECLEN);

            VEC_MUL4_ACCUMX(A0, B0, B1, B2, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, B1, B2, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, B1, B2, zero, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, B1, B2, zero, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

            A0 = loadu64(c + ((i)*BLOCKWORDS + j + 0) * VECLEN);
            A1 = loadu64(c + ((i)*BLOCKWORDS + j + 1) * VECLEN);
            A2 = loadu64(c + ((i)*BLOCKWORDS + j + 2) * VECLEN);
            A3 = loadu64(c + ((i)*BLOCKWORDS + j + 3) * VECLEN);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 0) * VECLEN,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 1) * VECLEN,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 2) * VECLEN,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 3) * VECLEN,
                _mm512_and_si512(lo3, MASK));

            hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

            // subtract out all of the bias at once.
            SUB_BIAS_HI3(3, 2, 1);
            SUB_BIAS_LO3(3, 2, 1);

            A0 = loadu64(c + ((i)*BLOCKWORDS + j + 4) * VECLEN);
            A1 = loadu64(c + ((i)*BLOCKWORDS + j + 5) * VECLEN);
            A2 = loadu64(c + ((i)*BLOCKWORDS + j + 6) * VECLEN);
            A3 = loadu64(c + ((i)*BLOCKWORDS + j + 7) * VECLEN);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 4) * VECLEN,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 5) * VECLEN,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 6) * VECLEN,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 7) * VECLEN,
                _mm512_and_si512(hi6, MASK));

            if (((i)*BLOCKWORDS + j + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 8) * VECLEN,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i)*BLOCKWORDS + j + 8) * VECLEN),
                        _mm512_srli_epi64(hi6, 52)));
            }
        }

        j = full_blocks * BLOCKWORDS;
        B0 = loadu64(a + (j + 0) * VECLEN);
        B1 = loadu64(a + (j + 1) * VECLEN);
        B2 = loadu64(a + (j + 2) * VECLEN);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            if ((i * BLOCKWORDS + j) > n)
                continue;

            __m512i A0 = loadu64(b + (i * BLOCKWORDS + 0) * VECLEN);
            __m512i A1 = loadu64(b + (i * BLOCKWORDS + 1) * VECLEN);
            __m512i A2 = loadu64(b + (i * BLOCKWORDS + 2) * VECLEN);
            __m512i A3 = loadu64(b + (i * BLOCKWORDS + 3) * VECLEN);

            VEC_MUL4_ACCUMX(A0, B0, B1, B2, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, B1, B2, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, B1, B2, zero, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, B1, B2, zero, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

            A0 = loadu64(c + ((i)*BLOCKWORDS + j + 0) * VECLEN);
            A1 = loadu64(c + ((i)*BLOCKWORDS + j + 1) * VECLEN);
            A2 = loadu64(c + ((i)*BLOCKWORDS + j + 2) * VECLEN);
            A3 = loadu64(c + ((i)*BLOCKWORDS + j + 3) * VECLEN);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 0) * VECLEN,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 1) * VECLEN,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 2) * VECLEN,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 3) * VECLEN,
                _mm512_and_si512(lo3, MASK));

            hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

            // subtract out all of the bias at once.
            SUB_BIAS_HI3(3, 2, 1);
            SUB_BIAS_LO3(3, 2, 1);

            A0 = loadu64(c + ((i)*BLOCKWORDS + j + 4) * VECLEN);
            A1 = loadu64(c + ((i)*BLOCKWORDS + j + 5) * VECLEN);
            A2 = loadu64(c + ((i)*BLOCKWORDS + j + 6) * VECLEN);
            A3 = loadu64(c + ((i)*BLOCKWORDS + j + 7) * VECLEN);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 4) * VECLEN,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 5) * VECLEN,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 6) * VECLEN,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 7) * VECLEN,
                _mm512_and_si512(hi6, MASK));

            if (((i)*BLOCKWORDS + j + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 8) * VECLEN,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i)*BLOCKWORDS + j + 8) * VECLEN),
                        _mm512_srli_epi64(hi6, 52)));
            }
        }
    }

    if (n - (full_blocks * BLOCKWORDS) == 2)
    {
        j = full_blocks * BLOCKWORDS;
        __m512i B0 = loadu64(b + (j + 0) * VECLEN);
        __m512i B1 = loadu64(b + (j + 1) * VECLEN);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            if ((i * BLOCKWORDS + j) > n)
                continue;

            __m512i A0 = loadu64(a + (i * BLOCKWORDS + 0) * VECLEN);
            __m512i A1 = loadu64(a + (i * BLOCKWORDS + 1) * VECLEN);
            __m512i A2 = loadu64(a + (i * BLOCKWORDS + 2) * VECLEN);
            __m512i A3 = loadu64(a + (i * BLOCKWORDS + 3) * VECLEN);

            VEC_MUL4_ACCUMX(A0, B0, B1, zero, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, B1, zero, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, B1, zero, zero, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, B1, zero, zero, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

#ifdef DEBUG_VECMUL
            print_regvechex(lo0, 0, "lo0: ");
            print_regvechex(lo1, 0, "lo1: ");
            print_regvechex(lo2, 0, "lo2: ");
            print_regvechex(lo3, 0, "lo3: ");

            print_regvechex(hi0, 0, "hi0: ");
            print_regvechex(hi1, 0, "hi1: ");
            print_regvechex(hi2, 0, "hi2: ");
            print_regvechex(hi3, 0, "hi3: ");
#endif

            A0 = loadu64(c + ((i)*BLOCKWORDS + j + 0) * VECLEN);
            A1 = loadu64(c + ((i)*BLOCKWORDS + j + 1) * VECLEN);
            A2 = loadu64(c + ((i)*BLOCKWORDS + j + 2) * VECLEN);
            A3 = loadu64(c + ((i)*BLOCKWORDS + j + 3) * VECLEN);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 0) * VECLEN,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 1) * VECLEN,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 2) * VECLEN,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 3) * VECLEN,
                _mm512_and_si512(lo3, MASK));

            hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

            // subtract out all of the bias at once.
            SUB_BIAS_HI3(3, 2, 1);
            SUB_BIAS_LO3(3, 2, 1);

#ifdef DEBUG_VECMUL
            print_regvechex(lo4, 0, "lo4: ");
            print_regvechex(lo5, 0, "lo5: ");
            print_regvechex(lo6, 0, "lo6: ");
            print_regvechex(hi4, 0, "hi4: ");
            print_regvechex(hi5, 0, "hi5: ");
            print_regvechex(hi6, 0, "hi6: ");
#endif

            A0 = loadu64(c + ((i)*BLOCKWORDS + j + 4) * VECLEN);
            A1 = loadu64(c + ((i)*BLOCKWORDS + j + 5) * VECLEN);
            A2 = loadu64(c + ((i)*BLOCKWORDS + j + 6) * VECLEN);
            A3 = loadu64(c + ((i)*BLOCKWORDS + j + 7) * VECLEN);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 4) * VECLEN,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 5) * VECLEN,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 6) * VECLEN,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 7) * VECLEN,
                _mm512_and_si512(hi6, MASK));

            if (((i)*BLOCKWORDS + j + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 8) * VECLEN,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i)*BLOCKWORDS + j + 8) * VECLEN),
                        _mm512_srli_epi64(hi6, 52)));
            }
        }

        B0 = loadu64(a + (j + 0) * VECLEN);
        B1 = loadu64(a + (j + 1) * VECLEN);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            if ((i * BLOCKWORDS + j) > n)
                continue;

            __m512i A0 = loadu64(b + (i * BLOCKWORDS + 0) * VECLEN);
            __m512i A1 = loadu64(b + (i * BLOCKWORDS + 1) * VECLEN);
            __m512i A2 = loadu64(b + (i * BLOCKWORDS + 2) * VECLEN);
            __m512i A3 = loadu64(b + (i * BLOCKWORDS + 3) * VECLEN);

            VEC_MUL4_ACCUMX(A0, B0, B1, zero, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, B1, zero, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, B1, zero, zero, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, B1, zero, zero, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

            A0 = loadu64(c + ((i)*BLOCKWORDS + j + 0) * VECLEN);
            A1 = loadu64(c + ((i)*BLOCKWORDS + j + 1) * VECLEN);
            A2 = loadu64(c + ((i)*BLOCKWORDS + j + 2) * VECLEN);
            A3 = loadu64(c + ((i)*BLOCKWORDS + j + 3) * VECLEN);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 0) * VECLEN,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 1) * VECLEN,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 2) * VECLEN,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 3) * VECLEN,
                _mm512_and_si512(lo3, MASK));

            hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

            // subtract out all of the bias at once.
            SUB_BIAS_HI3(3, 2, 1);
            SUB_BIAS_LO3(3, 2, 1);

            A0 = loadu64(c + ((i)*BLOCKWORDS + j + 4) * VECLEN);
            A1 = loadu64(c + ((i)*BLOCKWORDS + j + 5) * VECLEN);
            A2 = loadu64(c + ((i)*BLOCKWORDS + j + 6) * VECLEN);
            A3 = loadu64(c + ((i)*BLOCKWORDS + j + 7) * VECLEN);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 4) * VECLEN,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 5) * VECLEN,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 6) * VECLEN,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 7) * VECLEN,
                _mm512_and_si512(hi6, MASK));

            if (((i)*BLOCKWORDS + j + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 8) * VECLEN,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i)*BLOCKWORDS + j + 8) * VECLEN),
                        _mm512_srli_epi64(hi6, 52)));
            }
        }
    }

    if (n - (full_blocks * BLOCKWORDS) == 1)
    {
        j = full_blocks * BLOCKWORDS;
        __m512i B0 = loadu64(b + (j + 0) * VECLEN);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            if ((i * BLOCKWORDS + j) > n)
                continue;

            __m512i A0 = loadu64(a + (i * BLOCKWORDS + 0) * VECLEN);
            __m512i A1 = loadu64(a + (i * BLOCKWORDS + 1) * VECLEN);
            __m512i A2 = loadu64(a + (i * BLOCKWORDS + 2) * VECLEN);
            __m512i A3 = loadu64(a + (i * BLOCKWORDS + 3) * VECLEN);

            VEC_MUL4_ACCUMX(A0, B0, zero, zero, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, zero, zero, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, zero, zero, zero, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, zero, zero, zero, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

            A0 = loadu64(c + ((i)*BLOCKWORDS + j + 0) * VECLEN);
            A1 = loadu64(c + ((i)*BLOCKWORDS + j + 1) * VECLEN);
            A2 = loadu64(c + ((i)*BLOCKWORDS + j + 2) * VECLEN);
            A3 = loadu64(c + ((i)*BLOCKWORDS + j + 3) * VECLEN);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 0) * VECLEN,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 1) * VECLEN,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 2) * VECLEN,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 3) * VECLEN,
                _mm512_and_si512(lo3, MASK));

            hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

            // subtract out all of the bias at once.
            SUB_BIAS_HI3(3, 2, 1);
            SUB_BIAS_LO3(3, 2, 1);

#ifdef DEBUG_VECMUL
            print_regvechex(lo4, 0, "lo4: ");
            print_regvechex(lo5, 0, "lo5: ");
            print_regvechex(lo6, 0, "lo6: ");
            print_regvechex(hi4, 0, "hi4: ");
            print_regvechex(hi5, 0, "hi5: ");
            print_regvechex(hi6, 0, "hi6: ");
#endif

            A0 = loadu64(c + ((i)*BLOCKWORDS + j + 4) * VECLEN);
            A1 = loadu64(c + ((i)*BLOCKWORDS + j + 5) * VECLEN);
            A2 = loadu64(c + ((i)*BLOCKWORDS + j + 6) * VECLEN);
            A3 = loadu64(c + ((i)*BLOCKWORDS + j + 7) * VECLEN);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 4) * VECLEN,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 5) * VECLEN,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 6) * VECLEN,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 7) * VECLEN,
                _mm512_and_si512(hi6, MASK));

            if (((i)*BLOCKWORDS + j + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 8) * VECLEN,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i)*BLOCKWORDS + j + 8) * VECLEN),
                        _mm512_srli_epi64(hi6, 52)));
            }
        }

        B0 = loadu64(a + (j + 0) * VECLEN);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            if ((i * BLOCKWORDS + j) > n)
                continue;

            __m512i A0 = loadu64(b + (i * BLOCKWORDS + 0) * VECLEN);
            __m512i A1 = loadu64(b + (i * BLOCKWORDS + 1) * VECLEN);
            __m512i A2 = loadu64(b + (i * BLOCKWORDS + 2) * VECLEN);
            __m512i A3 = loadu64(b + (i * BLOCKWORDS + 3) * VECLEN);

            VEC_MUL4_ACCUMX(A0, B0, zero, zero, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, zero, zero, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, zero, zero, zero, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, zero, zero, zero, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

            A0 = loadu64(c + ((i)*BLOCKWORDS + j + 0) * VECLEN);
            A1 = loadu64(c + ((i)*BLOCKWORDS + j + 1) * VECLEN);
            A2 = loadu64(c + ((i)*BLOCKWORDS + j + 2) * VECLEN);
            A3 = loadu64(c + ((i)*BLOCKWORDS + j + 3) * VECLEN);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 0) * VECLEN,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 1) * VECLEN,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 2) * VECLEN,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 3) * VECLEN,
                _mm512_and_si512(lo3, MASK));

            hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

            // subtract out all of the bias at once.
            SUB_BIAS_HI3(3, 2, 1);
            SUB_BIAS_LO3(3, 2, 1);

            A0 = loadu64(c + ((i)*BLOCKWORDS + j + 4) * VECLEN);
            A1 = loadu64(c + ((i)*BLOCKWORDS + j + 5) * VECLEN);
            A2 = loadu64(c + ((i)*BLOCKWORDS + j + 6) * VECLEN);
            A3 = loadu64(c + ((i)*BLOCKWORDS + j + 7) * VECLEN);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 4) * VECLEN,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 5) * VECLEN,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 6) * VECLEN,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 7) * VECLEN,
                _mm512_and_si512(hi6, MASK));

            if (((i)*BLOCKWORDS + j + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i)*BLOCKWORDS + j + 8) * VECLEN,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i)*BLOCKWORDS + j + 8) * VECLEN),
                        _mm512_srli_epi64(hi6, 52)));
            }
        }

    }


#ifdef DEBUG_VECMUL
    print_vechex(c, 0, 2 * n, "end of multiplication:");
#endif

#ifdef DEBUG_VECMUL
    exit(1);
#endif
    return;
}

#define _mm512_madd52lo_epu64_(r, a, b, c, o) \
    { \
        r=fma52lo(a, b, _mm512_loadu_si512((__m512i*)(((char*)c)+o))); \
    }
#define _mm512_madd52hi_epu64_(r, a, b, c, d, v, o) \
    { \
        r=fma52hi(a, b, _mm512_loadu_si512((__m512i*)(((char*)c)+o)), d, v); \
    }
#define fma52lo_mem(r, a, b, c, o) _mm512_madd52lo_epu64_(r, a, b, c, o)
#define fma52hi_mem(r, a, b, c, d, v, o) _mm512_madd52hi_epu64_(r, a, b, c, d, v, o)

#define fma52lo_mem_x4(a1, a2, a3, a4, b, c, o, m) { \
    __m512i t1 = _mm512_loadu_si512((__m512i*)(((char*)c) + o + 0 * SIMD_BYTES)); \
    __m512i t2 = _mm512_loadu_si512((__m512i*)(((char*)c) + o + 1 * SIMD_BYTES)); \
    __m512i t3 = _mm512_loadu_si512((__m512i*)(((char*)c) + o + 2 * SIMD_BYTES)); \
    __m512i t4 = _mm512_loadu_si512((__m512i*)(((char*)c) + o + 3 * SIMD_BYTES)); \
    t1 = _mm512_mullo_epi64(b, t1); \
    t2 = _mm512_mullo_epi64(b, t2); \
    t3 = _mm512_mullo_epi64(b, t3); \
    t4 = _mm512_mullo_epi64(b, t4); \
    t1 = _mm512_and_si512(t1, m); \
    t2 = _mm512_and_si512(t2, m); \
    t3 = _mm512_and_si512(t3, m); \
    t4 = _mm512_and_si512(t4, m); \
    a1 = _mm512_add_epi64(a1, t1); \
    a2 = _mm512_add_epi64(a2, t2); \
    a3 = _mm512_add_epi64(a3, t3); \
    a4 = _mm512_add_epi64(a4, t4); \
}

#define fma52hi_mem_x4(a1, a2, a3, a4, b, c, d, v, o) { \
    __m512d bd = _mm512_cvtepu64_pd(b); \
    __m512d cd1 = _mm512_cvtepu64_pd(_mm512_loadu_si512((__m512i*)(((char*)c) + o + 0 * SIMD_BYTES))); \
    __m512d cd2 = _mm512_cvtepu64_pd(_mm512_loadu_si512((__m512i*)(((char*)c) + o + 1 * SIMD_BYTES))); \
    __m512d cd3 = _mm512_cvtepu64_pd(_mm512_loadu_si512((__m512i*)(((char*)c) + o + 2 * SIMD_BYTES))); \
    __m512d cd4 = _mm512_cvtepu64_pd(_mm512_loadu_si512((__m512i*)(((char*)c) + o + 3 * SIMD_BYTES))); \
    cd1 = _mm512_fmadd_round_pd(bd, cd1, d, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
    cd2 = _mm512_fmadd_round_pd(bd, cd2, d, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
    cd3 = _mm512_fmadd_round_pd(bd, cd3, d, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
    cd4 = _mm512_fmadd_round_pd(bd, cd4, d, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
    a1 = _mm512_add_epi64(a1, _mm512_sub_epi64(castepu(cd1), v)); \
    a2 = _mm512_add_epi64(a2, _mm512_sub_epi64(castepu(cd2), v)); \
    a3 = _mm512_add_epi64(a3, _mm512_sub_epi64(castepu(cd3), v)); \
    a4 = _mm512_add_epi64(a4, _mm512_sub_epi64(castepu(cd4), v)); \
}

#define fma52hi_mem_x4_shift(a1, a2, a3, a4, a5, b, c, d, v, o) { \
    __m512d bd = _mm512_cvtepu64_pd(b); \
    __m512d cd1 = _mm512_cvtepu64_pd(_mm512_loadu_si512((__m512i*)(((char*)c) + o + 0 * SIMD_BYTES))); \
    __m512d cd2 = _mm512_cvtepu64_pd(_mm512_loadu_si512((__m512i*)(((char*)c) + o + 1 * SIMD_BYTES))); \
    __m512d cd3 = _mm512_cvtepu64_pd(_mm512_loadu_si512((__m512i*)(((char*)c) + o + 2 * SIMD_BYTES))); \
    __m512d cd4 = _mm512_cvtepu64_pd(_mm512_loadu_si512((__m512i*)(((char*)c) + o + 3 * SIMD_BYTES))); \
    cd1 = _mm512_fmadd_round_pd(bd, cd1, d, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
    cd2 = _mm512_fmadd_round_pd(bd, cd2, d, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
    cd3 = _mm512_fmadd_round_pd(bd, cd3, d, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
    cd4 = _mm512_fmadd_round_pd(bd, cd4, d, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
    a1 = _mm512_add_epi64(a2, _mm512_sub_epi64(castepu(cd1), v)); \
    a2 = _mm512_add_epi64(a3, _mm512_sub_epi64(castepu(cd2), v)); \
    a3 = _mm512_add_epi64(a4, _mm512_sub_epi64(castepu(cd3), v)); \
    a4 = _mm512_add_epi64(a5, _mm512_sub_epi64(castepu(cd4), v)); \
}



void vecmulmod52_fixed624_cios(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    int i, j, k;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t NBLOCKS = mdata->NBLOCKS;

//#ifndef IFMA
    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(set64(0x4670000000000000ULL));
    __m512i vbias1 = set64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = set64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = set64(0x4330000000000000ULL);  // 31
//#endif

    // needed after loops
    __m512i zero = set64(0);
    __m512i MASK = set64(DIGIT_MASK);
    __m512i K = loadu64(mdata->vrho);
    __mmask8 scarry2;
    __mmask8 scarry;

    // deal with the sign
    c->size = NWORDS;
    c->signmask = a->signmask ^ b->signmask;

    __m512i res00, res01, res02, res03, res04, res05, res06, res07, res08, res09,
        res10, res11;

    res00 = res01 = res02 = res03 = res04 = res05 = res06 = res07 = res08 =
        res09 = res10 = res11 = zero;

    uint64_t* inpB_mb = b->data;
    uint64_t* inpA_mb = a->data;
    uint64_t* inpM_mb = n->data;
    uint64_t* out_mb = c->data;

    int itr;
    for (itr = 0; itr < 12; itr++) {
        __m512i Yi;
        __m512i Bi = loadu64(inpB_mb);
        inpB_mb += MB_WIDTH;

        fma52lo_mem_x4(res00, res01, res02, res03, Bi, inpA_mb, SIMD_BYTES * 0, MASK);
        fma52lo_mem_x4(res04, res05, res06, res07, Bi, inpA_mb, SIMD_BYTES * 4, MASK);
        fma52lo_mem_x4(res08, res09, res10, res11, Bi, inpA_mb, SIMD_BYTES * 8, MASK);

        Yi = fma52lo(zero, res00, K);

        fma52lo_mem_x4(res00, res01, res02, res03, Yi, inpM_mb, SIMD_BYTES * 0, MASK);
        fma52lo_mem_x4(res04, res05, res06, res07, Yi, inpM_mb, SIMD_BYTES * 4, MASK);
        fma52lo_mem_x4(res08, res09, res10, res11, Yi, inpM_mb, SIMD_BYTES * 8, MASK);

        res00 = srli64(res00, DIGIT_SIZE);
        res01 = add64(res01, res00);

        fma52hi_mem_x4_shift(res00, res01, res02, res03, res04, Bi, inpA_mb, dbias, vbias1, SIMD_BYTES * 0);
        fma52hi_mem_x4_shift(res04, res05, res06, res07, res08, Bi, inpA_mb, dbias, vbias1, SIMD_BYTES * 4);
        fma52hi_mem_x4_shift(res08, res09, res10, res11, zero, Bi, inpA_mb, dbias, vbias1, SIMD_BYTES * 8);

        fma52hi_mem_x4(res00, res01, res02, res03, Yi, inpM_mb, dbias, vbias1, SIMD_BYTES * 0);
        fma52hi_mem_x4(res04, res05, res06, res07, Yi, inpM_mb, dbias, vbias1, SIMD_BYTES * 4);
        fma52hi_mem_x4(res08, res09, res10, res11, Yi, inpM_mb, dbias, vbias1, SIMD_BYTES * 8);
    }


    // Normalization
    __m512i T = zero;

    {
        T = srli64(res00, DIGIT_SIZE);
        res00 = and64(res00, MASK);
        storeu64(out_mb + MB_WIDTH * 0, res00);
        res01 = add64(res01, T);
        T = srli64(res01, DIGIT_SIZE);
        res01 = and64(res01, MASK);
        storeu64(out_mb + MB_WIDTH * 1, res01);
        res02 = add64(res02, T);
        T = srli64(res02, DIGIT_SIZE);
        res02 = and64(res02, MASK);
        storeu64(out_mb + MB_WIDTH * 2, res02);
        res03 = add64(res03, T);
        T = srli64(res03, DIGIT_SIZE);
        res03 = and64(res03, MASK);
        storeu64(out_mb + MB_WIDTH * 3, res03);
        res04 = add64(res04, T);
        T = srli64(res04, DIGIT_SIZE);
        res04 = and64(res04, MASK);
        storeu64(out_mb + MB_WIDTH * 4, res04);
        res05 = add64(res05, T);
        T = srli64(res05, DIGIT_SIZE);
        res05 = and64(res05, MASK);
        storeu64(out_mb + MB_WIDTH * 5, res05);
        res06 = add64(res06, T);
        T = srli64(res06, DIGIT_SIZE);
        res06 = and64(res06, MASK);
        storeu64(out_mb + MB_WIDTH * 6, res06);
        res07 = add64(res07, T);
        T = srli64(res07, DIGIT_SIZE);
        res07 = and64(res07, MASK);
        storeu64(out_mb + MB_WIDTH * 7, res07);
        res08 = add64(res08, T);
        T = srli64(res08, DIGIT_SIZE);
        res08 = and64(res08, MASK);
        storeu64(out_mb + MB_WIDTH * 8, res08);
        res09 = add64(res09, T);
        T = srli64(res09, DIGIT_SIZE);
        res09 = and64(res09, MASK);
        storeu64(out_mb + MB_WIDTH * 9, res09);
        res10 = add64(res10, T);
        T = srli64(res10, DIGIT_SIZE);
        res10 = and64(res10, MASK);
        storeu64(out_mb + MB_WIDTH * 10, res10);
        res11 = add64(res11, T);
        T = srli64(res11, DIGIT_SIZE);
        res11 = and64(res11, MASK);
        storeu64(out_mb + MB_WIDTH * 11, res11);
    }

    //printf("lane 0 cios: \n");
    //for (i = 0; i < NWORDS; i++)
    //{
    //    printf("%016lx", out_mb[MB_WIDTH * i + 0]);
    //}
    //printf("\n");
    //
    //print_regvechex64(T, 0, "carry reg");

    // AMM doesn't require further reduction mod N
    if (0)
    {
        scarry2 = _mm512_cmp_epu64_mask(T, zero, _MM_CMPINT_EQ);

        if (scarry2 != 0xff)
        {
            printf("carry word: %08x\n", scarry2);
        }

        // subtract n from tmp
        scarry = 0;
        for (i = 0; i < NWORDS; i++)
        {
            __m512i a1 = _mm512_load_epi64(s->data + i * VECLEN);
            __m512i b0 = _mm512_load_epi64(n->data + i * VECLEN);
            __m512i a0 = _mm512_sbb_epi52(a1, scarry, b0, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(MASK, a0));
        }

        // negate any final borrows if there was also a final carry.
        scarry &= scarry2;

        // if there was a final borrow, we didn't need to do the subtraction after all.
        // replace with original results based on final borrow mask.
        for (i = NWORDS - 1; i >= 0; i--)
        {
            __m512i b0 = _mm512_load_epi64(s->data + i * VECLEN);
            _mm512_mask_store_epi64(c->data + i * VECLEN, scarry, b0);
        }
    }

    //printf("lane 0 cios: \n");
    //for (i = 0; i < NWORDS; i++)
    //{
    //    printf("%016lx", c->data[MB_WIDTH * i + 0]);
    //}
    //printf("\n");

    c->size = 12;
    return;
}

void vecmul52_lo(uint64_t* a, uint64_t* b, uint64_t* c, int words)
{
    int i, j;
    uint32_t NBLOCKS;
    // needed in loops
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19

#ifndef IFMA
    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
#endif

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry;

    __m512i arestore[3];
    __m512i brestore[3];

    NBLOCKS = words / BLOCKWORDS;
    if ((words % BLOCKWORDS) > 0)
    {
        NBLOCKS++;

#ifdef DEBUG_VECMUL
        printf("commencing vecmul52 with WORDS=%d and NBLOCKS=%d\n", words, NBLOCKS);
        print_vechex(a, 0, NBLOCKS * BLOCKWORDS, "input a: ");
        print_vechex(b, 0, NBLOCKS * BLOCKWORDS, "input b: ");
        print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS, "output c: ");
#endif

        for (i = words, j = 0; i < NBLOCKS * BLOCKWORDS; i++, j++)
        {
            arestore[j] = _mm512_load_epi64(a + i * VECLEN);
            brestore[j] = _mm512_load_epi64(b + i * VECLEN);
            _mm512_store_epi64(a + i * VECLEN, zero);
            _mm512_store_epi64(b + i * VECLEN, zero);
        }
    }

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

#ifdef DEBUG_VECMUL
    printf("after padding top words\n");
    print_vechex(a, 0, NBLOCKS * BLOCKWORDS, "input a: ");
    print_vechex(b, 0, NBLOCKS * BLOCKWORDS, "input b: ");
    print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS, "output c: ");
#endif

    // first half mul
    for (i = 0; i < NBLOCKS; i++)
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (j = i; j > 0; j--)
        {
            a0 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum
        a0 = _mm512_load_epi64(a + (i * BLOCKWORDS + 0) * VECLEN);
        a1 = _mm512_load_epi64(a + (i * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a + (i * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a + (i * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b + 3 * VECLEN);
        b1 = _mm512_load_epi64(b + 2 * VECLEN);
        b2 = _mm512_load_epi64(b + 1 * VECLEN);
        b3 = _mm512_load_epi64(b + 0 * VECLEN);

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a0, b3, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a1, b3, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a0, b2, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a2, b3, te4, te5);
#else
        // ======
        //VEC_MUL_ACCUM_LOHI_PD(a0, b3, te0, te1);
        //VEC_MUL_ACCUM_LOHI_PD(a1, b3, te2, te3);
        //VEC_MUL_ACCUM_LOHI_PD(a0, b2, te2, te3);
        //VEC_MUL_ACCUM_LOHI_PD(a2, b3, te4, te5);
        {
            prod5_ld = _mm512_cvtepu64_pd(a0);
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(a2);
            prod3_ld = _mm512_cvtepu64_pd(b2);
            prod4_ld = _mm512_cvtepu64_pd(b3);

            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod1_hd));  // a1 * b3 -> to te2/3
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));  // a2 * b3 -> to te4/5 
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a0 * b2 -> to te2/3
            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod4_hd));  // a0 * b3 -> to te0/1

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod1_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod4_ld));
        }
#endif

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a1, b2, te4, te5);
        VEC_MUL_ACCUM_LOHI_PD(a0, b1, te4, te5);
        VEC_MUL_ACCUM_LOHI_PD(a1, b1, te6, te7);
        VEC_MUL_ACCUM_LOHI_PD(a0, b0, te6, te7);
#else
        // ======
        //VEC_MUL_ACCUM_LOHI_PD(a1, b2, te4, te5);
        //VEC_MUL_ACCUM_LOHI_PD(a0, b1, te4, te5);
        //VEC_MUL_ACCUM_LOHI_PD(a1, b1, te6, te7);
        //VEC_MUL_ACCUM_LOHI_PD(a0, b0, te6, te7);
        {
            prod5_ld = _mm512_cvtepu64_pd(a0);
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(b0);
            prod3_ld = _mm512_cvtepu64_pd(b1);
            prod4_ld = _mm512_cvtepu64_pd(b2);

            prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_hd = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod2_hd));
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod3_hd));
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod4_hd));

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_ld = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod2_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod3_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod4_ld));
        }
#endif

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a3, b3, te6, te7);
        VEC_MUL_ACCUM_LOHI_PD(a2, b2, te6, te7);
#else
        //VEC_MUL_ACCUM_LOHI_PD(a3, b3, te6, te7);
        //VEC_MUL_ACCUM_LOHI_PD(a2, b2, te6, te7);
        {
            prod1_ld = _mm512_cvtepu64_pd(a2);
            prod2_ld = _mm512_cvtepu64_pd(a3);
            prod3_ld = _mm512_cvtepu64_pd(b2);
            prod4_ld = _mm512_cvtepu64_pd(b3);

            prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod2_hd));

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

            prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod2_ld));
        }

        // subtract out all of the bias at once.  these are
        // counts of how many times we have mul_accumulated into each column
        // during the block a*b and final triangular a*b.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 2,
            i * 4 + 3,
            i * 4 + 4);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 2,
            i * 4 + 3,
            i * 4 + 4);
#endif



        // now, do a carry-propagating column summation and store the results.
        {
            j = 0;
            // accumulate this column-sum:
            // carry propagate low to high.
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(c + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 1;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te2);
            acc_e1 = _mm512_add_epi64(acc_e1, te3);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(c + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 2;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te4);
            acc_e1 = _mm512_add_epi64(acc_e1, te5);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(c + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 3;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te6);
            acc_e1 = _mm512_add_epi64(acc_e1, te7);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(c + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;
        }
    }

#ifdef DEBUG_VECMUL
    print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS, "after lo half:");
#endif

    for (i = words, j = 0; i < NBLOCKS * BLOCKWORDS; i++, j++)
    {
        _mm512_store_epi64(a + i * VECLEN, arestore[j]);
        _mm512_store_epi64(b + i * VECLEN, brestore[j]);
    }


#ifdef DEBUG_VECMUL
    printf("after restoring inputs and outputs\n");
    print_vechex(a, 0, NBLOCKS * BLOCKWORDS, "input a: ");
    print_vechex(b, 0, NBLOCKS * BLOCKWORDS, "input b: ");
    print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS, "output c: ");
    exit(1);
#endif

    return;
}

void vecmul52_n(uint64_t* a, uint64_t* b, uint64_t* c, int words)
{
    int i, j;
    uint32_t NBLOCKS;
    // needed in loops
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19

#ifndef IFMA
    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
#endif

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry;

    NBLOCKS = words / BLOCKWORDS;
    if ((words % BLOCKWORDS) > 0)
    {
        NBLOCKS++;

#ifdef DEBUG_VECMUL
        printf("commencing vecmul52 with WORDS=%d and NBLOCKS=%d\n", words, NBLOCKS);
        print_vechex(a, 0, NBLOCKS * BLOCKWORDS, "input a: ");
        print_vechex(b, 0, NBLOCKS * BLOCKWORDS, "input b: ");
        print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS, "output c: ");
#endif
    }

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

    // first half mul
    for (i = 0; i < NBLOCKS; i++)
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (j = i; j > 0; j--)
        {
            a0 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum
        a0 = _mm512_load_epi64(a + (i * BLOCKWORDS + 0) * VECLEN);
        a1 = _mm512_load_epi64(a + (i * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a + (i * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a + (i * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b + 3 * VECLEN);
        b1 = _mm512_load_epi64(b + 2 * VECLEN);
        b2 = _mm512_load_epi64(b + 1 * VECLEN);
        b3 = _mm512_load_epi64(b + 0 * VECLEN);

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a0, b3, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a1, b3, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a0, b2, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a2, b3, te4, te5);
#else
        // ======
        //VEC_MUL_ACCUM_LOHI_PD(a0, b3, te0, te1);
        //VEC_MUL_ACCUM_LOHI_PD(a1, b3, te2, te3);
        //VEC_MUL_ACCUM_LOHI_PD(a0, b2, te2, te3);
        //VEC_MUL_ACCUM_LOHI_PD(a2, b3, te4, te5);
        {
            prod5_ld = _mm512_cvtepu64_pd(a0);
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(a2);
            prod3_ld = _mm512_cvtepu64_pd(b2);
            prod4_ld = _mm512_cvtepu64_pd(b3);

            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod1_hd));  // a1 * b3 -> to te2/3
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));  // a2 * b3 -> to te4/5 
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a0 * b2 -> to te2/3
            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod4_hd));  // a0 * b3 -> to te0/1

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod1_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod4_ld));
        }
#endif

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a1, b2, te4, te5);
        VEC_MUL_ACCUM_LOHI_PD(a0, b1, te4, te5);
        VEC_MUL_ACCUM_LOHI_PD(a1, b1, te6, te7);
        VEC_MUL_ACCUM_LOHI_PD(a0, b0, te6, te7);
#else
        // ======
        //VEC_MUL_ACCUM_LOHI_PD(a1, b2, te4, te5);
        //VEC_MUL_ACCUM_LOHI_PD(a0, b1, te4, te5);
        //VEC_MUL_ACCUM_LOHI_PD(a1, b1, te6, te7);
        //VEC_MUL_ACCUM_LOHI_PD(a0, b0, te6, te7);
        {
            prod5_ld = _mm512_cvtepu64_pd(a0);
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(b0);
            prod3_ld = _mm512_cvtepu64_pd(b1);
            prod4_ld = _mm512_cvtepu64_pd(b2);

            prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_hd = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod2_hd));
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod3_hd));
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod4_hd));

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_ld = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod2_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod3_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod4_ld));
        }
#endif

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a3, b3, te6, te7);
        VEC_MUL_ACCUM_LOHI_PD(a2, b2, te6, te7);
#else
        //VEC_MUL_ACCUM_LOHI_PD(a3, b3, te6, te7);
        //VEC_MUL_ACCUM_LOHI_PD(a2, b2, te6, te7);
        {
            prod1_ld = _mm512_cvtepu64_pd(a2);
            prod2_ld = _mm512_cvtepu64_pd(a3);
            prod3_ld = _mm512_cvtepu64_pd(b2);
            prod4_ld = _mm512_cvtepu64_pd(b3);

            prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod2_hd));

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

            prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod2_ld));
        }

        // subtract out all of the bias at once.  these are
        // counts of how many times we have mul_accumulated into each column
        // during the block a*b and final triangular a*b.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 2,
            i * 4 + 3,
            i * 4 + 4);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 2,
            i * 4 + 3,
            i * 4 + 4);
#endif



        // now, do a carry-propagating column summation and store the results.
        {
            j = 0;
            // accumulate this column-sum:
            // carry propagate low to high.
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(c + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 1;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te2);
            acc_e1 = _mm512_add_epi64(acc_e1, te3);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(c + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 2;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te4);
            acc_e1 = _mm512_add_epi64(acc_e1, te5);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(c + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 3;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te6);
            acc_e1 = _mm512_add_epi64(acc_e1, te7);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(c + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;
        }
    }

#ifdef DEBUG_VECMUL
    print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS, "after lo half:");
#endif

    // second half mul
    for (i = NBLOCKS; i < 2 * NBLOCKS; i++)
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (j = i - NBLOCKS + 1; j < NBLOCKS; j++)
        {
            a0 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }


        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b + ((NBLOCKS * BLOCKWORDS) - 1) * VECLEN);
        b1 = _mm512_load_epi64(b + ((NBLOCKS * BLOCKWORDS) - 2) * VECLEN);
        b2 = _mm512_load_epi64(b + ((NBLOCKS * BLOCKWORDS) - 3) * VECLEN);

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a1, b0, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b1, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b0, te2, te3);
#else
        // ======
        //VEC_MUL_ACCUM_LOHI_PD(a1, b0, te0, te1);
        //VEC_MUL_ACCUM_LOHI_PD(a2, b1, te0, te1);
        //VEC_MUL_ACCUM_LOHI_PD(a2, b0, te2, te3);
        {
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(a2);
            prod3_ld = _mm512_cvtepu64_pd(b0);
            prod4_ld = _mm512_cvtepu64_pd(b1);

            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));  // a1*b0 -> t0/1
            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod4_hd));  // a2*b1 -> t0/1
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a2*b0 -> t2/3

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod4_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
        }
#endif

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a3, b2, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a3, b1, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a3, b0, te4, te5);
#else
        //VEC_MUL_ACCUM_LOHI_PD(a3, b2, te0, te1);
        //VEC_MUL_ACCUM_LOHI_PD(a3, b1, te2, te3);
        //VEC_MUL_ACCUM_LOHI_PD(a3, b0, te4, te5);
        {
            prod1_ld = _mm512_cvtepu64_pd(a3);
            prod2_ld = _mm512_cvtepu64_pd(b2);
            prod3_ld = _mm512_cvtepu64_pd(b1);
            prod4_ld = _mm512_cvtepu64_pd(b0);

            prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));  // a3*b2 -> t0/1
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a3*b1 -> t2/3
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod4_hd));  // a3*b0 -> t4/5

            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod4_ld));
        }


        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            (2 * NBLOCKS - i - 1) * 4 + 3,
            (2 * NBLOCKS - i - 1) * 4 + 2,
            (2 * NBLOCKS - i - 1) * 4 + 1,
            (2 * NBLOCKS - i - 1) * 4 + 0);
        SUB_BIAS_LO(
            (2 * NBLOCKS - i - 1) * 4 + 3,
            (2 * NBLOCKS - i - 1) * 4 + 2,
            (2 * NBLOCKS - i - 1) * 4 + 1,
            (2 * NBLOCKS - i - 1) * 4 + 0);
#endif

        // final column accumulation
        {
            j = 0;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a3 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN, 
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 1;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te2);
            acc_e1 = _mm512_add_epi64(acc_e1, te3);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a2 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 2;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te4);
            acc_e1 = _mm512_add_epi64(acc_e1, te5);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a1 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 3;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te6);
            acc_e1 = _mm512_add_epi64(acc_e1, te7);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a0 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            a3 = _mm512_and_epi64(vlmask, a3);
            a2 = _mm512_and_epi64(vlmask, a2);
            a1 = _mm512_and_epi64(vlmask, a1);
            a0 = _mm512_and_epi64(vlmask, a0);

            _mm512_store_epi64(c + (i * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(c + (i * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(c + (i * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(c + (i * BLOCKWORDS + 3) * VECLEN, a0);

        }
    }

#ifdef DEBUG_VECMUL
    print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS, "after hi half:");
#endif

    return;
}

void vecmul52_np(uint64_t* a, uint64_t* b, uint64_t* c, int words)
{
    int i, j;
    uint32_t NBLOCKS;
    // needed in loops
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19

#ifndef IFMA
    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
#endif

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry;

    __m512i ahi[BLOCKWORDS];
    __m512i bhi[BLOCKWORDS];

    NBLOCKS = words / BLOCKWORDS;
    if ((words % BLOCKWORDS) == 0)
    {
        printf("error only call vecmul52_np when 4 does not divide n\n");
    }

    for (i = NBLOCKS * BLOCKWORDS, j = 0; i < words; i++, j++)
    {
        ahi[j] = _mm512_load_epi64(a + i * VECLEN);
        bhi[j] = _mm512_load_epi64(b + i * VECLEN);
    }
    for (; j < BLOCKWORDS; j++)
    {
        ahi[j] = zero;
        bhi[j] = zero;
    }
    NBLOCKS++;


    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

    // first half mul
    for (i = 0; i < NBLOCKS; i++)
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (j = i; j > 0; j--)
        {

            a0 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            
            if (j == (NBLOCKS-1))
            {
                b3 = bhi[0];
                b4 = bhi[1];
                b5 = bhi[2];
                b6 = bhi[3];
            }
            else
            {
                b3 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
                b4 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
                b5 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
                b6 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 7) * VECLEN);
            }

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum
        if (i == (NBLOCKS-1))
        {
            a0 = ahi[0];
            a1 = ahi[1];
            a2 = ahi[2];
            a3 = ahi[3];
        }
        else
        {
            a0 = _mm512_load_epi64(a + (i * BLOCKWORDS + 0) * VECLEN);
            a1 = _mm512_load_epi64(a + (i * BLOCKWORDS + 1) * VECLEN);
            a2 = _mm512_load_epi64(a + (i * BLOCKWORDS + 2) * VECLEN);
            a3 = _mm512_load_epi64(a + (i * BLOCKWORDS + 3) * VECLEN);
        }

        b0 = _mm512_load_epi64(b + 3 * VECLEN);
        b1 = _mm512_load_epi64(b + 2 * VECLEN);
        b2 = _mm512_load_epi64(b + 1 * VECLEN);
        b3 = _mm512_load_epi64(b + 0 * VECLEN);

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a0, b3, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a1, b3, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a0, b2, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a2, b3, te4, te5);
#else
        // ======
        //VEC_MUL_ACCUM_LOHI_PD(a0, b3, te0, te1);
        //VEC_MUL_ACCUM_LOHI_PD(a1, b3, te2, te3);
        //VEC_MUL_ACCUM_LOHI_PD(a0, b2, te2, te3);
        //VEC_MUL_ACCUM_LOHI_PD(a2, b3, te4, te5);
        {
            prod5_ld = _mm512_cvtepu64_pd(a0);
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(a2);
            prod3_ld = _mm512_cvtepu64_pd(b2);
            prod4_ld = _mm512_cvtepu64_pd(b3);

            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod1_hd));  // a1 * b3 -> to te2/3
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));  // a2 * b3 -> to te4/5 
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a0 * b2 -> to te2/3
            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod4_hd));  // a0 * b3 -> to te0/1

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod1_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod4_ld));
        }
#endif

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a1, b2, te4, te5);
        VEC_MUL_ACCUM_LOHI_PD(a0, b1, te4, te5);
        VEC_MUL_ACCUM_LOHI_PD(a1, b1, te6, te7);
        VEC_MUL_ACCUM_LOHI_PD(a0, b0, te6, te7);
#else
        // ======
        //VEC_MUL_ACCUM_LOHI_PD(a1, b2, te4, te5);
        //VEC_MUL_ACCUM_LOHI_PD(a0, b1, te4, te5);
        //VEC_MUL_ACCUM_LOHI_PD(a1, b1, te6, te7);
        //VEC_MUL_ACCUM_LOHI_PD(a0, b0, te6, te7);
        {
            prod5_ld = _mm512_cvtepu64_pd(a0);
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(b0);
            prod3_ld = _mm512_cvtepu64_pd(b1);
            prod4_ld = _mm512_cvtepu64_pd(b2);

            prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_hd = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod2_hd));
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod3_hd));
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod4_hd));

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_ld = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod2_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod3_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod4_ld));
        }
#endif

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a3, b3, te6, te7);
        VEC_MUL_ACCUM_LOHI_PD(a2, b2, te6, te7);
#else
        //VEC_MUL_ACCUM_LOHI_PD(a3, b3, te6, te7);
        //VEC_MUL_ACCUM_LOHI_PD(a2, b2, te6, te7);
        {
            prod1_ld = _mm512_cvtepu64_pd(a2);
            prod2_ld = _mm512_cvtepu64_pd(a3);
            prod3_ld = _mm512_cvtepu64_pd(b2);
            prod4_ld = _mm512_cvtepu64_pd(b3);

            prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod2_hd));

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

            prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod2_ld));
        }

        // subtract out all of the bias at once.  these are
        // counts of how many times we have mul_accumulated into each column
        // during the block a*b and final triangular a*b.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 2,
            i * 4 + 3,
            i * 4 + 4);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 2,
            i * 4 + 3,
            i * 4 + 4);
#endif

        // now, do a carry-propagating column summation and store the results.
        {
            j = 0;
            // accumulate this column-sum:
            // carry propagate low to high.
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(c + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 1;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te2);
            acc_e1 = _mm512_add_epi64(acc_e1, te3);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(c + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 2;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te4);
            acc_e1 = _mm512_add_epi64(acc_e1, te5);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(c + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 3;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te6);
            acc_e1 = _mm512_add_epi64(acc_e1, te7);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(c + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;
        }
    }

#ifdef DEBUG_VECMUL
    print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS, "after lo half:");
#endif

    // second half mul
    for (i = NBLOCKS; i < 2 * NBLOCKS; i++)
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (j = i - NBLOCKS + 1; j < NBLOCKS; j++)
        {
            if (j == (i - NBLOCKS + 1)) // the first j-iteration
            {
                a0 = ahi[3];
                a1 = ahi[2];
                a2 = ahi[1];
                a3 = ahi[0];
            }
            else
            {
                a0 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 3) * VECLEN);
                a1 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 2) * VECLEN);
                a2 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 1) * VECLEN);
                a3 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 0) * VECLEN);
            }

            b0 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            if (j == (NBLOCKS - 1))
            {
                b3 = bhi[0];
                b4 = bhi[1];
                b5 = bhi[2];
                b6 = bhi[3];
            }
            else
            {
                b3 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
                b4 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
                b5 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
                b6 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 7) * VECLEN);
            }

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum (a * b)
        if (i == (2 * NBLOCKS - 1))
        {
            a1 = ahi[1];
            a2 = ahi[2];
            a3 = ahi[3];
        }
        else
        {
            a1 = _mm512_load_epi64(a + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
            a2 = _mm512_load_epi64(a + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
            a3 = _mm512_load_epi64(a + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);
        }

        b0 = bhi[3];
        b1 = bhi[2];
        b2 = bhi[1];

        //b0 = _mm512_load_epi64(b + ((NBLOCKS * BLOCKWORDS) - 1) * VECLEN);
        //b1 = _mm512_load_epi64(b + ((NBLOCKS * BLOCKWORDS) - 2) * VECLEN);
        //b2 = _mm512_load_epi64(b + ((NBLOCKS * BLOCKWORDS) - 3) * VECLEN);

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a1, b0, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b1, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b0, te2, te3);
#else
        // ======
        //VEC_MUL_ACCUM_LOHI_PD(a1, b0, te0, te1);
        //VEC_MUL_ACCUM_LOHI_PD(a2, b1, te0, te1);
        //VEC_MUL_ACCUM_LOHI_PD(a2, b0, te2, te3);
        {
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(a2);
            prod3_ld = _mm512_cvtepu64_pd(b0);
            prod4_ld = _mm512_cvtepu64_pd(b1);

            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));  // a1*b0 -> t0/1
            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod4_hd));  // a2*b1 -> t0/1
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a2*b0 -> t2/3

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod4_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
        }
#endif

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a3, b2, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a3, b1, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a3, b0, te4, te5);
#else
        //VEC_MUL_ACCUM_LOHI_PD(a3, b2, te0, te1);
        //VEC_MUL_ACCUM_LOHI_PD(a3, b1, te2, te3);
        //VEC_MUL_ACCUM_LOHI_PD(a3, b0, te4, te5);
        {
            prod1_ld = _mm512_cvtepu64_pd(a3);
            prod2_ld = _mm512_cvtepu64_pd(b2);
            prod3_ld = _mm512_cvtepu64_pd(b1);
            prod4_ld = _mm512_cvtepu64_pd(b0);

            prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));  // a3*b2 -> t0/1
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a3*b1 -> t2/3
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod4_hd));  // a3*b0 -> t4/5

            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod4_ld));
        }


        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            (2 * NBLOCKS - i - 1) * 4 + 3,
            (2 * NBLOCKS - i - 1) * 4 + 2,
            (2 * NBLOCKS - i - 1) * 4 + 1,
            (2 * NBLOCKS - i - 1) * 4 + 0);
        SUB_BIAS_LO(
            (2 * NBLOCKS - i - 1) * 4 + 3,
            (2 * NBLOCKS - i - 1) * 4 + 2,
            (2 * NBLOCKS - i - 1) * 4 + 1,
            (2 * NBLOCKS - i - 1) * 4 + 0);
#endif

        // final column accumulation
        {
            j = 0;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a3 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN, 
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 1;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te2);
            acc_e1 = _mm512_add_epi64(acc_e1, te3);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a2 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 2;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te4);
            acc_e1 = _mm512_add_epi64(acc_e1, te5);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a1 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 3;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te6);
            acc_e1 = _mm512_add_epi64(acc_e1, te7);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a0 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            a3 = _mm512_and_epi64(vlmask, a3);
            a2 = _mm512_and_epi64(vlmask, a2);
            a1 = _mm512_and_epi64(vlmask, a1);
            a0 = _mm512_and_epi64(vlmask, a0);

            if ((i * BLOCKWORDS + 0) < words)
                _mm512_store_epi64(c + (i * BLOCKWORDS + 0) * VECLEN, a3);
            if ((i * BLOCKWORDS + 1) < words)
                _mm512_store_epi64(c + (i * BLOCKWORDS + 1) * VECLEN, a2);
            if ((i * BLOCKWORDS + 2) < words)
                _mm512_store_epi64(c + (i * BLOCKWORDS + 2) * VECLEN, a1);
            if ((i * BLOCKWORDS + 3) < words)
                _mm512_store_epi64(c + (i * BLOCKWORDS + 3) * VECLEN, a0);

        }
    }

#ifdef DEBUG_VECMUL
    print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS, "after hi half:");
#endif


    return;
}

void vecsqr52_n(uint64_t* a, uint64_t* c, int words)
{
    int i, j, k;
    uint64_t* b = a;
    uint32_t NBLOCKS;
    __m512i i0, i1;
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19

#ifndef IFMA
    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
#endif

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry;

    __m512i arestore[3];
    __m512i crestore[6];

    NBLOCKS = words / BLOCKWORDS;
    if ((words % BLOCKWORDS) > 0)
    {
        NBLOCKS++;

#ifdef DEBUG_VECMUL
        printf("commencing vecsqr52 with WORDS=%d and NBLOCKS=%d\n", words, NBLOCKS);
        print_vechex(a, 0, NBLOCKS * BLOCKWORDS, "input a: ");
        print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS, "output c: ");
#endif

        for (i = words, j = 0; i < NBLOCKS * BLOCKWORDS; i++, j++)
        {
            arestore[j] = _mm512_load_epi64(a + i * VECLEN);
            _mm512_store_epi64(a + i * VECLEN, zero);
        }

        for (i = 2 * words, j = 0; i < 2 * NBLOCKS * BLOCKWORDS; i++, j++)
        {
            crestore[j] = _mm512_load_epi64(c + i * VECLEN);
            _mm512_store_epi64(c + i * VECLEN, zero);
        }
    }

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

#ifdef DEBUG_VECMUL
    printf("after padding top words\n");
    print_vechex(a, 0, NBLOCKS * BLOCKWORDS, "input a: ");
    print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS, "output c: ");
#endif

    // first half sqr
    for (i = 0; i < NBLOCKS; i++)
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (k = 0, j = i; j > (i + 1) / 2; j--, k++)
        {
            // when i = 0, j = 0, j > 0 --> 0 iterations
            // when i = 1, j = 1, j > 1 --> 0 iterations
            // when i = 2, j = 2, j > 1 --> 1 iteration @ a[3..0], b[5..11]
            // when i = 3, j = 3, j > 2 --> 1 iteration @ a[3..0], b[9..15]
            // when i = 4, j = 4, j > 2 --> 2 iteration @ a[3..0], b[13..19] and a[7..4], b[9..15]
            // when i = 5, j = 5, j > 3 --> 2 iteration @ a[3..0], b[17..23] and a[7..4], b[13..19]
            a0 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        if (i & 1)
        {
            // i odd
            // for 384-bit inputs when i == 1, j = 0, a = {3,2,1,0} and b = {2,3,4,5,6,7}
            // for 512-bit inputs when i == 3, j = 1, a = {7,6,5,4} and b = {6,7,8,9,a,b}
            // for 512-bit inputs when i == 1, j = 0, a = {3,2,1,0} and b = {2,3,4,5,6,7}

            a0 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b1 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a2, b2, te0, te1);
            VEC_MUL_ACCUM_LOHI_PD(a1, b2, te2, te3);
            VEC_MUL_ACCUM_LOHI_PD(a1, b3, te4, te5);
            VEC_MUL_ACCUM_LOHI_PD(a0, b3, te6, te7);
#else
            //prod1_e = _mm512_mul_epu32(a2, b2);   // te0
            //prod2_e = _mm512_mul_epu32(a1, b2);   // te2
            //prod3_e = _mm512_mul_epu32(a1, b3);   // te4
            //prod4_e = _mm512_mul_epu32(a0, b3);   // te6
            //ACCUM_4X_DOUBLED_PROD;
            {
                prod5_ld = _mm512_cvtepu64_pd(a0);
                prod1_ld = _mm512_cvtepu64_pd(a1);
                prod2_ld = _mm512_cvtepu64_pd(a2);
                prod3_ld = _mm512_cvtepu64_pd(b2);
                prod4_ld = _mm512_cvtepu64_pd(b3);

                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b2 -> to te0/1 
                prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b2 -> to te2/3
                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b3 -> to te4/5
                prod4_hd = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b3 -> to te6/7

                te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod1_hd));
                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod4_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
                prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod4_ld = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod1_ld));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod4_ld));
            }
#endif


#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a3, b3, te0, te1);
            VEC_MUL_ACCUM_LOHI_PD(a2, b3, te2, te3);
            VEC_MUL_ACCUM_LOHI_PD(a2, b4, te4, te5);
            VEC_MUL_ACCUM_LOHI_PD(a1, b4, te6, te7);
#else
            //prod1_e = _mm512_mul_epu32(a3, b3);   // te0
            //prod2_e = _mm512_mul_epu32(a2, b3);   // te2
            //prod3_e = _mm512_mul_epu32(a2, b4);   // te4
            //prod4_e = _mm512_mul_epu32(a1, b4);   // te6
            //ACCUM_4X_DOUBLED_PROD;
            {
                prod5_ld = _mm512_cvtepu64_pd(a1);
                prod1_ld = _mm512_cvtepu64_pd(a2);
                prod2_ld = _mm512_cvtepu64_pd(a3);
                prod3_ld = _mm512_cvtepu64_pd(b3);
                prod4_ld = _mm512_cvtepu64_pd(b4);

                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * b3 -> to te0/1
                prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b3 -> to te2/3
                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b4 -> to te4/5
                prod4_hd = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b4 -> to te6/7

                te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod1_hd));
                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod4_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
                prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod4_ld = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod1_ld));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod4_ld));
            }

#endif

#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a3, b4, te2, te3);
            VEC_MUL_ACCUM_LOHI_PD(a3, b5, te4, te5);
            VEC_MUL_ACCUM_LOHI_PD(a2, b5, te6, te7);
            VEC_MUL_ACCUM_LOHI_PD(a3, b6, te6, te7);
#else
            //prod1_e = _mm512_mul_epu32(a3, b4);  // te2
            //prod2_e = _mm512_mul_epu32(a3, b5);  // te4
            //prod3_e = _mm512_mul_epu32(a2, b5);  // te6
            {
                prod1_ld = _mm512_cvtepu64_pd(a2);
                prod2_ld = _mm512_cvtepu64_pd(a3);
                prod3_ld = _mm512_cvtepu64_pd(b4);
                prod4_ld = _mm512_cvtepu64_pd(b5);
                prod5_ld = _mm512_cvtepu64_pd(b6);

                prod3_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * b4 -> to te2/3
                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * b5 -> to te4/5
                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b5 -> to te6/7
                prod4_hd = _mm512_fmadd_round_pd(prod2_ld, prod5_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * b6 -> to te6/7

                te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));
                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod4_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
                prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

                prod5_ld = _mm512_fmadd_round_pd(prod2_ld, prod5_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod3_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod5_ld));
            }

            //prod1_e = _mm512_mul_epu32(a3, b6);  // te6
            //{
            //    prod1_ld = _mm512_cvtepu64_pd(a3);
            //    prod2_ld = _mm512_cvtepu64_pd(b6);
            //
            //    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
            //        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * b6 -> to te6/7
            //
            //    te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
            //    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            //    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            //    te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            //}

            // all terms so far need to be doubled.  
            // but to do that we first need to remove all bias
            // that has been accumulated so far.
            SUB_BIAS_HI(
                k * 4 + 2,
                k * 4 + 3,
                k * 4 + 3,
                k * 4 + 4);
            SUB_BIAS_LO(
                k * 4 + 2,
                k * 4 + 3,
                k * 4 + 3,
                k * 4 + 4);
#endif

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a1, a1, te0, te1);
            VEC_MUL_ACCUM_LOHI_PD(a0, a0, te4, te5);
#else
            // finally, accumulate the two non-doubled terms.
            //prod1_e = _mm512_mul_epu32(a1, a1);    // te0
            //prod2_e = _mm512_mul_epu32(a0, a0);    // te4
            {
                prod1_ld = _mm512_cvtepu64_pd(a0);
                prod2_ld = _mm512_cvtepu64_pd(a1);

                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * a1 -> to te0/1
                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a0 -> to te4/5

                te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod1_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod1_ld));
            }
#endif
        }
        else
        {
            // i even
            a0 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 0) * VECLEN);
            a1 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a2 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a3 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS + 3) * VECLEN);

#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a0, a1, te2, te3);
            VEC_MUL_ACCUM_LOHI_PD(a0, a2, te4, te5);
            VEC_MUL_ACCUM_LOHI_PD(a0, a3, te6, te7);
            VEC_MUL_ACCUM_LOHI_PD(a1, a2, te6, te7);
#else
            //prod1_e = _mm512_mul_epu32(a0, a1);    // te2
            //prod2_e = _mm512_mul_epu32(a0, a2);    // te4
            //prod3_e = _mm512_mul_epu32(a0, a3);    // te6
            {
                prod1_ld = _mm512_cvtepu64_pd(a0);
                prod2_ld = _mm512_cvtepu64_pd(a1);
                prod3_ld = _mm512_cvtepu64_pd(a2);
                prod4_ld = _mm512_cvtepu64_pd(a3);

                prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a1 -> to te2/3
                prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a2 -> to te4/5
                prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a3 -> to te6/7
                prod1_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * a2 -> to te6/7

                te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod2_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod3_hd));
                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod4_hd));
                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));

                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
                prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);
                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);

                prod5_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod2_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod3_ld));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod4_ld));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod5_ld));
            }

            ////prod1_e = _mm512_mul_epu32(a1, a2);    // te6
            //{
            //    prod1_ld = _mm512_cvtepu64_pd(a1);
            //    prod2_ld = _mm512_cvtepu64_pd(a2);
            //
            //    
            //
            //    te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
            //    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            //    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            //    te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            //}

            // all terms so far need to be doubled.  
            // but to do that we first need to remove all bias
            // that has been accumulated so far.
            SUB_BIAS_HI(
                k * 4 + 0,
                k * 4 + 1,
                k * 4 + 1,
                k * 4 + 2);
            SUB_BIAS_LO(
                k * 4 + 0,
                k * 4 + 1,
                k * 4 + 1,
                k * 4 + 2);

#endif

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a0, a0, te0, te1);
            VEC_MUL_ACCUM_LOHI_PD(a1, a1, te4, te5);
#else
            // finally, accumulate the two non-doubled terms.
            //prod1_e = _mm512_mul_epu32(a0, a0);
            //prod2_e = _mm512_mul_epu32(a1, a1);
            {
                prod1_ld = _mm512_cvtepu64_pd(a0);
                prod2_ld = _mm512_cvtepu64_pd(a1);

                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a0 -> to te0/1
                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * a1 -> to te4/5

                te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
            }
#endif
        }

#ifndef IFMA
        // need to remove bias from two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            1,
            0,
            1,
            0);
        SUB_BIAS_LO(
            1,
            0,
            1,
            0);
#endif


        // now, do a carry-propagating column summation and store the results.
        {
            j = 0;
            // accumulate this column-sum:
            // carry propagate low to high.
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(c + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 1;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te2);
            acc_e1 = _mm512_add_epi64(acc_e1, te3);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(c + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 2;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te4);
            acc_e1 = _mm512_add_epi64(acc_e1, te5);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(c + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 3;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te6);
            acc_e1 = _mm512_add_epi64(acc_e1, te7);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(c + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;
        }
    }

#ifdef DEBUG_VECMUL
    print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS, "after lo half:");
#endif

    // second half sqr
    for (i = 0; i < NBLOCKS; i++)
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (k = 0, j = 0; j < (NBLOCKS - i - 1) / 2; j++, k++)
        {
            // Compute a solid block (all matching terms are in the lower
            // half triangle of the expansion).
            a0 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS) - 1 - j * BLOCKWORDS) * VECLEN);
            a1 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS) - 2 - j * BLOCKWORDS) * VECLEN);
            a2 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS) - 3 - j * BLOCKWORDS) * VECLEN);
            a3 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS) - 4 - j * BLOCKWORDS) * VECLEN);

            b0 = _mm512_load_epi64(b + ((j + i) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b + ((j + i) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b + ((j + i) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b + ((j + i) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b + ((j + i) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b + ((j + i) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b + ((j + i) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }


        if (NBLOCKS & 1)		// NBLOCKS is odd
        {
            // i odd, block shape 2.
            if (i & 1)
            {
                // always a continuation of the full-block loop, so use the same 
                // loading pattern.  Only now we don't need as many b-terms.
                a0 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS) - 1 - j * BLOCKWORDS) * VECLEN);
                a1 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS) - 2 - j * BLOCKWORDS) * VECLEN);
                a2 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS) - 3 - j * BLOCKWORDS) * VECLEN);
                a3 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS) - 4 - j * BLOCKWORDS) * VECLEN);

                b0 = _mm512_load_epi64(b + (j * BLOCKWORDS + i * BLOCKWORDS + 1) * VECLEN);
                b1 = _mm512_load_epi64(b + (j * BLOCKWORDS + i * BLOCKWORDS + 2) * VECLEN);
                b2 = _mm512_load_epi64(b + (j * BLOCKWORDS + i * BLOCKWORDS + 3) * VECLEN);
                b3 = _mm512_load_epi64(b + (j * BLOCKWORDS + i * BLOCKWORDS + 4) * VECLEN);
                b4 = _mm512_load_epi64(b + (j * BLOCKWORDS + i * BLOCKWORDS + 5) * VECLEN);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a0, b0, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a1, b2, te2, te3);
                VEC_MUL_ACCUM_LOHI_PD(a1, b1, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a0, b1, te2, te3);
#else
                //prod1_e = _mm512_mul_epu32(a0, b0);        // te0
                //prod1_e = _mm512_mul_epu32(a1, b1);        // te0
                //prod1_e = _mm512_mul_epu32(a0, b1);        // te2
                //prod1_e = _mm512_mul_epu32(a1, b2);        // te2
                {
                    prod5_ld = _mm512_cvtepu64_pd(a0);
                    prod1_ld = _mm512_cvtepu64_pd(a1);
                    prod2_ld = _mm512_cvtepu64_pd(b0);
                    prod3_ld = _mm512_cvtepu64_pd(b1);
                    prod4_ld = _mm512_cvtepu64_pd(b2);

                    prod2_hd = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b0 -> to te0/1
                    prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b2 -> to te2/3
                    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b1 -> to te0/1
                    prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b1 -> to te2/3

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod4_hd));
                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));
                    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));

                    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
                    prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

                    prod2_ld = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod4_ld));
                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
                    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                }

#endif

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a2, b2, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a2, b3, te2, te3);
#else
                //prod1_e = _mm512_mul_epu32(a2, b2);        // te0
                //prod1_e = _mm512_mul_epu32(a2, b3);        // te2
                {
                    prod1_ld = _mm512_cvtepu64_pd(a2);
                    prod2_ld = _mm512_cvtepu64_pd(b2);
                    prod3_ld = _mm512_cvtepu64_pd(b3);

                    prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b2 -> to te0/1
                    prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b3 -> to te2/3

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));

                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);

                    prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                }
#endif

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a0, b2, te4, te5);
                VEC_MUL_ACCUM_LOHI_PD(a1, b4, te6, te7);
                VEC_MUL_ACCUM_LOHI_PD(a1, b3, te4, te5);
                VEC_MUL_ACCUM_LOHI_PD(a0, b3, te6, te7);
#else
                // prod1_e = _mm512_mul_epu32(a0, b2);       // te4
                // prod1_e = _mm512_mul_epu32(a1, b3);       // te4
                // prod1_e = _mm512_mul_epu32(a0, b3);       // te6
                // prod1_e = _mm512_mul_epu32(a1, b4);       // te6
                {
                    prod5_ld = _mm512_cvtepu64_pd(a0);
                    prod1_ld = _mm512_cvtepu64_pd(a1);
                    prod2_ld = _mm512_cvtepu64_pd(b2);
                    prod3_ld = _mm512_cvtepu64_pd(b3);
                    prod4_ld = _mm512_cvtepu64_pd(b4);

                    prod2_hd = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b2 -> to te4/5
                    prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b4 -> to te6/7
                    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b3 -> to te4/5
                    prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b3 -> to te6/7

                    te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));
                    te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod4_hd));
                    te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod1_hd));
                    te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod3_hd));

                    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
                    prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

                    prod2_ld = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
                    te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod4_ld));
                    te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod1_ld));
                    te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod3_ld));
                }


                // all terms so far need to be doubled.  
                // but to do that we first need to remove all bias
                // that has been accumulated so far.
                SUB_BIAS_HI(
                    k * 4 + 3,
                    k * 4 + 3,
                    k * 4 + 2,
                    k * 4 + 2);
                SUB_BIAS_LO(
                    k * 4 + 3,
                    k * 4 + 3,
                    k * 4 + 2,
                    k * 4 + 2);

#endif

                // now double
                te1 = _mm512_slli_epi64(te1, 1);
                te3 = _mm512_slli_epi64(te3, 1);
                te5 = _mm512_slli_epi64(te5, 1);
                te7 = _mm512_slli_epi64(te7, 1);
                te0 = _mm512_slli_epi64(te0, 1);
                te2 = _mm512_slli_epi64(te2, 1);
                te4 = _mm512_slli_epi64(te4, 1);
                te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a3, a3, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a2, a2, te4, te5);
#else
                // finally the two non-doubled terms.
                //prod1_e = _mm512_mul_epu32(a3, a3);    // te0
                //prod1_e = _mm512_mul_epu32(a2, a2);    // te4
                {
                    prod1_ld = _mm512_cvtepu64_pd(a3);
                    prod2_ld = _mm512_cvtepu64_pd(a2);

                    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * a3 -> to te0/1
                    prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * a2 -> to te4/5

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));
                    te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));

                    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
                    te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
                }
#endif
            }
            else
            {
                // i even, block shape 1.
                // always a continuation of the full-block loop, so use the same 
                // loading pattern.  Only now we don't need as many b-terms.
                a0 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS) - 1 - j * BLOCKWORDS) * VECLEN);
                a1 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS) - 2 - j * BLOCKWORDS) * VECLEN);
                a2 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS) - 3 - j * BLOCKWORDS) * VECLEN);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a0, a2, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a0, a1, te2, te3);
#else
                //k == 0;
                //prod1_e = _mm512_mul_epu32(a0, a2);      // te0
                //prod1_e = _mm512_mul_epu32(a0, a1);      // te2
                {
                    prod1_ld = _mm512_cvtepu64_pd(a0);
                    prod2_ld = _mm512_cvtepu64_pd(a1);
                    prod3_ld = _mm512_cvtepu64_pd(a2);

                    prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a2 -> to te0/1
                    prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a1 -> to te2/3

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod3_hd));
                    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod2_hd));

                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);

                    prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod3_ld));
                    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod2_ld));
                }

                // all terms so far need to be doubled.  
                // but to do that we first need to remove all bias
                // that has been accumulated so far.
                SUB_BIAS_HI(
                    k * 4 + 1,
                    k * 4 + 1,
                    k * 4 + 0,
                    k * 4 + 0);
                SUB_BIAS_LO(
                    k * 4 + 1,
                    k * 4 + 1,
                    k * 4 + 0,
                    k * 4 + 0);

#endif

                // now double
                te1 = _mm512_slli_epi64(te1, 1);
                te3 = _mm512_slli_epi64(te3, 1);
                te5 = _mm512_slli_epi64(te5, 1);
                te7 = _mm512_slli_epi64(te7, 1);
                te0 = _mm512_slli_epi64(te0, 1);
                te2 = _mm512_slli_epi64(te2, 1);
                te4 = _mm512_slli_epi64(te4, 1);
                te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a1, a1, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a0, a0, te4, te5);
#else
                // finally the two non-doubled terms.
                //prod1_e = _mm512_mul_epu32(a1, a1);    // te0
                //prod1_e = _mm512_mul_epu32(a0, a0);    // te4
                {
                    prod1_ld = _mm512_cvtepu64_pd(a1);
                    prod2_ld = _mm512_cvtepu64_pd(a0);

                    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b1 -> to te0/1
                    prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b0 -> to te4/5

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));
                    te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));

                    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
                    te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
                }



#endif
            }

        }
        else				// NBLOCKS is even
        {
            // i odd, block shape 1.
            if (i & 1)
            {
                // always a continuation of the full-block loop, so use the same 
                // loading pattern.  Only now we don't need as many b-terms.
                a0 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS) - 1 - j * BLOCKWORDS) * VECLEN);  // {f, b}
                a1 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS) - 2 - j * BLOCKWORDS) * VECLEN);  // {e, a}
                a2 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS) - 3 - j * BLOCKWORDS) * VECLEN);  // {d, 9}

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a0, a2, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a0, a1, te2, te3);
#else
                //k == 0;
                //prod1_e = _mm512_mul_epu32(a0, a2);   // te0
                //prod2_e = _mm512_mul_epu32(a0, a1);   // te2
                {
                    prod1_ld = _mm512_cvtepu64_pd(a0);
                    prod2_ld = _mm512_cvtepu64_pd(a1);
                    prod3_ld = _mm512_cvtepu64_pd(a2);

                    prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a2 -> to te0/1
                    prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a1 -> to te2/3

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod3_hd));
                    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod2_hd));

                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);

                    prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod3_ld));
                    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod2_ld));
                }

                // all terms so far need to be doubled.  
                // but to do that we first need to remove all bias
                // that has been accumulated so far.
                SUB_BIAS_HI(
                    j * 4 + 1,
                    j * 4 + 1,
                    j * 4 + 0,
                    j * 4 + 0);
                SUB_BIAS_LO(
                    j * 4 + 1,
                    j * 4 + 1,
                    j * 4 + 0,
                    j * 4 + 0);

#endif

                // now double
                te1 = _mm512_slli_epi64(te1, 1);
                te3 = _mm512_slli_epi64(te3, 1);
                te5 = _mm512_slli_epi64(te5, 1);
                te7 = _mm512_slli_epi64(te7, 1);
                te0 = _mm512_slli_epi64(te0, 1);
                te2 = _mm512_slli_epi64(te2, 1);
                te4 = _mm512_slli_epi64(te4, 1);
                te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a1, a1, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a0, a0, te4, te5);
#else
                // finally, accumulate the two non-doubled terms.
                //prod1_e = _mm512_mul_epu32(a1, a1);       // te0
                //prod2_e = _mm512_mul_epu32(a0, a0);       // te4
                {
                    prod1_ld = _mm512_cvtepu64_pd(a1);
                    prod2_ld = _mm512_cvtepu64_pd(a0);

                    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b1 -> to te0/1
                    prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b0 -> to te4/5

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));
                    te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));

                    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
                    te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
                }
#endif
            }
            else
            {
                // i even, block shape 1.
                // always a continuation of the full-block loop, so use the same 
                // loading pattern.  Only now we don't need as many b-terms.
                a0 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS) - 1 - j * BLOCKWORDS) * VECLEN);		// {f, b}
                a1 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS) - 2 - j * BLOCKWORDS) * VECLEN);		// {e, a}
                a2 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS) - 3 - j * BLOCKWORDS) * VECLEN);		// {d, 9}
                a3 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS) - 4 - j * BLOCKWORDS) * VECLEN);		// {c, 8}

                b0 = _mm512_load_epi64(b + (j * BLOCKWORDS + i * BLOCKWORDS + 1) * VECLEN); // {9, 5}
                b1 = _mm512_load_epi64(b + (j * BLOCKWORDS + i * BLOCKWORDS + 2) * VECLEN);	// {a, 6}
                b2 = _mm512_load_epi64(b + (j * BLOCKWORDS + i * BLOCKWORDS + 3) * VECLEN);	// {b, 7}
                b3 = _mm512_load_epi64(b + (j * BLOCKWORDS + i * BLOCKWORDS + 4) * VECLEN);	// {c, 8}
                b4 = _mm512_load_epi64(b + (j * BLOCKWORDS + i * BLOCKWORDS + 5) * VECLEN); // {d, 9}

                //prod1_e = _mm512_mul_epu32(a0, b0);    // te0
                //prod2_e = _mm512_mul_epu32(a0, b1);    // te2
                //prod3_e = _mm512_mul_epu32(a0, b2);    // te4
                //prod4_e = _mm512_mul_epu32(a0, b3);    // te6
                //ACCUM_4X_DOUBLED_PROD;
                VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);

                //prod1_e = _mm512_mul_epu32(a1, b1);    // te0
                //prod2_e = _mm512_mul_epu32(a1, b2);    // te2
                //prod3_e = _mm512_mul_epu32(a1, b3);    // te4
                //prod4_e = _mm512_mul_epu32(a1, b4);    // te6
                //ACCUM_4X_DOUBLED_PROD;
                VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a2, b2, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a2, b3, te2, te3);
#else
                //prod1_e = _mm512_mul_epu32(a2, b2);    // te0
                //prod2_e = _mm512_mul_epu32(a2, b3);    // te2
                {
                    prod1_ld = _mm512_cvtepu64_pd(a2);
                    prod2_ld = _mm512_cvtepu64_pd(b2);
                    prod3_ld = _mm512_cvtepu64_pd(b3);

                    prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b2 -> to te0/1
                    prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b3 -> to te2/3

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));

                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);

                    prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                }

                // all terms so far need to be doubled.  
                // but to do that we first need to remove all bias
                // that has been accumulated so far.
                SUB_BIAS_HI(
                    j * 4 + 3,
                    j * 4 + 3,
                    j * 4 + 2,
                    j * 4 + 2);
                SUB_BIAS_LO(
                    j * 4 + 3,
                    j * 4 + 3,
                    j * 4 + 2,
                    j * 4 + 2);

#endif

                // now double
                te1 = _mm512_slli_epi64(te1, 1);
                te3 = _mm512_slli_epi64(te3, 1);
                te5 = _mm512_slli_epi64(te5, 1);
                te7 = _mm512_slli_epi64(te7, 1);
                te0 = _mm512_slli_epi64(te0, 1);
                te2 = _mm512_slli_epi64(te2, 1);
                te4 = _mm512_slli_epi64(te4, 1);
                te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a3, a3, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a2, a2, te4, te5);
#else
                // finally, accumulate the two non-doubled terms.
                //prod1_e = _mm512_mul_epu32(a3, a3);   // te0
                //prod2_e = _mm512_mul_epu32(a2, a2);   // te4
                {
                    prod1_ld = _mm512_cvtepu64_pd(a3);
                    prod2_ld = _mm512_cvtepu64_pd(a2);

                    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * a3 -> to te0/1
                    prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * a2 -> to te4/5

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));
                    te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));

                    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
                    te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
                }
#endif
            }

        }

#ifndef IFMA
        // need to remove bias from the non-doubled terms.
        SUB_BIAS_HI(
            1,
            0,
            1,
            0);
        SUB_BIAS_LO(
            1,
            0,
            1,
            0);
#endif

        // final column accumulation
        {
            j = 0;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a3 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN, 
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 1;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te2);
            acc_e1 = _mm512_add_epi64(acc_e1, te3);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a2 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 2;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te4);
            acc_e1 = _mm512_add_epi64(acc_e1, te5);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a1 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 3;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te6);
            acc_e1 = _mm512_add_epi64(acc_e1, te7);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a0 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            a3 = _mm512_and_epi64(vlmask, a3);
            a2 = _mm512_and_epi64(vlmask, a2);
            a1 = _mm512_and_epi64(vlmask, a1);
            a0 = _mm512_and_epi64(vlmask, a0);

            _mm512_store_epi64(c + ((NBLOCKS + i) * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(c + ((NBLOCKS + i) * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(c + ((NBLOCKS + i) * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(c + ((NBLOCKS + i) * BLOCKWORDS + 3) * VECLEN, a0);

        }
    }

#ifdef DEBUG_VECMUL
    print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS, "after hi half:");
#endif

    for (i = words, j = 0; i < NBLOCKS * BLOCKWORDS; i++, j++)
    {
        _mm512_store_epi64(a + i * VECLEN, arestore[j]);
    }

    for (i = 2 * words, j = 0; i < 2 * NBLOCKS * BLOCKWORDS; i++, j++)
    {
        _mm512_store_epi64(c + i * VECLEN, crestore[j]);
    }

#ifdef DEBUG_VECMUL
    printf("after restoring inputs and outputs\n");
    print_vechex(a, 0, NBLOCKS* BLOCKWORDS, "input a: ");
    print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS, "output c: ");
    //exit(1);
#endif

    return;
}

void kmuln(uint64_t* a, uint64_t* b, uint64_t* c, uint64_t* scratch, int words);
void kmuln_lo(uint64_t* a, uint64_t* b, uint64_t* c, uint64_t* scratch, int words);
void ksqrn(uint64_t* a, uint64_t* c, uint64_t* scratch, int words);

#define veckmul(a, b, c, s, words)				\
  do {									\
    if (words <= 20) {						\
        if ((words & 0x3) == 0) { vecmul52_n(a, b, c, words); } \
        else { vecmul52_cios(a, b, c, words); } }     \
    else								\
        kmuln(a, b, c, s, words);				\
  } while (0);

#define veckmul_lo(a, b, c, s, words)				\
  do {									\
    if (words <= 20) {						\
        if ((words & 0x3) == 0) { vecmul52_lo(a, b, c, words); } \
        else { vecmul52_cios_lo(a, b, c, words); } }     \
    else								\
        kmuln_lo(a, b, c, s, words);				\
  } while (0);

#define vecksqr(a, c, s, words)				\
  do {									\
    if (words <= 20) {						\
        if ((words & 0x3) == 0) { vecsqr52_n(a, c, words); } \
        else { vecmul52_cios(a, a, c, words); } }     \
    else								\
        ksqrn(a, c, s, words);				\
  } while (0);

#define vectmul(a, b, c, s, words)				\
  do {									\
    if (words <= 20) {						\
        if ((words & 0x3) == 0) { vecmul52_n(a, b, c, words); } \
        else { vecmul52_cios(a, b, c, words); } }     \
    else								\
        tmuln(a, b, c, s, words);				\
  } while (0);

void tmuln(uint64_t* a, uint64_t* b, uint64_t* c, uint64_t* scratch, int words)
{
    int w0, w1, w2;

    w0 = w1 = words / 3;
    w2 = words - w1 - w0;
    int wmax = MAX(w0, w2);
    int scratchwords = wmax + 2;

    //printf("Splitting phase\n");
    uint64_t* a0 = a;
    uint64_t* a1 = a + w0 * VECLEN;
    uint64_t* a2 = a + (w0 + w1) * VECLEN;
    uint64_t* b0 = b;
    uint64_t* b1 = b + w1 * VECLEN;
    uint64_t* b2 = b + (w0 + w1) * VECLEN;
    uint64_t* p1 = scratch + 0 * scratchwords;
    uint64_t* pm1 = scratch + 1 * scratchwords;
    uint64_t* pm2 = scratch + 2 * scratchwords;
    uint64_t* q1 = scratch + 3 * scratchwords;
    uint64_t* qm1 = scratch + 4 * scratchwords;
    uint64_t* qm2 = scratch + 5 * scratchwords;
    __mmask8 m1;
    __mmask8 m2;
    __mmask8 m3;
    __mmask8 m4;

    //printf("Evaluation phase\n");
    //mpz_set(p0, a0);
    //mpz_set(pinf, a2);
    //mpz_add(pm1, a0, a2);
    //mpz_add(p1, pm1, a1);
    //mpz_sub(pm1, pm1, a1);
    //mpz_add(pm2, pm1, a2);
    //mpz_mul_2exp(pm2, pm2, 1);
    //mpz_sub(pm2, pm2, a0);
    base_add_52(a0, a2, pm1, w0, w2);
    base_add_52(p1, pm1, a1, wmax, w1);
    m1 = vec_gte2_52(a1, pm1, wmax);
    base_abssub_52(pm1, a1, pm1, m1, wmax, w1);
    base_add_52(pm1, a2, pm2, wmax, w2);
    base_add_52(pm2, pm2, pm2, wmax, wmax);
    m2 = vec_gte2_52(a0, pm2, wmax);
    base_abssub_52(pm2, a0, pm2, m2, wmax, w0);

    //mpz_set(q0, b0);
    //mpz_set(qinf, b2);
    //mpz_add(qm1, b0, b2);
    //mpz_add(q1, qm1, b1);
    //mpz_sub(qm1, qm1, b1);
    //mpz_add(qm2, qm1, b2);
    //mpz_mul_2exp(qm2, qm2, 1);
    //mpz_sub(qm2, qm2, b0);
    base_add_52(b0, b2, qm1, w0, w2);
    base_add_52(q1, qm1, b1, wmax, w1);
    m3 = vec_gte2_52(b1, qm1, wmax);
    base_abssub_52(qm1, b1, qm1, m3, wmax, w1);
    base_add_52(qm1, b2, qm2, wmax, w2);
    base_add_52(qm2, qm2, qm2, wmax, wmax);
    m4 = vec_gte2_52(b0, qm2, wmax);
    base_abssub_52(qm2, b0, qm2, m4, wmax, w0);

    //printf("Pointwise multiplication phase\n");
    //mpz_mul(r0, p0, q0);
    //mpz_mul(r1, p1, q1);
    //mpz_mul(rm1, pm1, qm1);
    //mpz_mul(rm2, pm2, qm2);
    //mpz_mul(rinf, pinf, qinf);

    //printf("Interpolation phase\n");
    //mpz_set(i0, r0);
    //mpz_set(i4, rinf);
    //mpz_sub(i3, rm2, r1);
    //mpz_tdiv_q_ui(i3, i3, 3);
    //mpz_sub(i1, r1, rm1);
    //mpz_tdiv_q_2exp(i1, i1, 1);
    //mpz_sub(i2, rm1, r0);
    //mpz_sub(i3, i2, i3);
    //mpz_tdiv_q_2exp(i3, i3, 1);
    //mpz_mul_2exp(rinf, rinf, 1);
    //mpz_add(i3, i3, rinf);
    //mpz_add(i2, i2, i1);
    //mpz_sub(i2, i2, i4);
    //mpz_sub(i1, i1, i3);

    //printf("Recomposition phase\n");
    //mpz_set(out, i0);
    //mpz_mul(i1, i1, B);
    //mpz_mul(i2, i2, B);
    //mpz_mul(i2, i2, B);
    //mpz_mul(i3, i3, B);
    //mpz_mul(i3, i3, B);
    //mpz_mul(i3, i3, B);
    //mpz_mul(i4, i4, B);
    //mpz_mul(i4, i4, B);
    //mpz_mul(i4, i4, B);
    //mpz_mul(i4, i4, B);
    //mpz_add(out, out, i1);
    //mpz_add(out, out, i2);
    //mpz_add(out, out, i3);
    //mpz_add(out, out, i4);



    return;
}


void kmuln(uint64_t* a, uint64_t* b, uint64_t* c, uint64_t* scratch, int words)
{
    int lowords, hiwords, i;

    if (words & 1)
    {
        hiwords = words / 2;
        lowords = hiwords + 1;
    }
    else
    {
        hiwords = lowords = words / 2;
    }

    uint64_t* a0 = a;
    uint64_t* a1 = a + lowords * VECLEN;
    uint64_t* b0 = b;
    uint64_t* b1 = b + lowords * VECLEN;
    uint64_t* z1;
    uint64_t* s1;
    uint64_t* s2;
    __mmask8 m1;
    __mmask8 m2;

    s1 = scratch;
    s2 = scratch + (lowords)*VECLEN;
    z1 = scratch + (lowords * 2) * VECLEN;

    veckmul(a0, b0, c, scratch, lowords);
    veckmul(a1, b1, c + 2 * lowords * VECLEN, scratch, hiwords);
    // result is negative if the longer subtractand does not have a leading digit
    // AND if the rest of the comparison is true.
    m1 = hiwords < lowords ? (_mm512_cmpeq_epu64_mask(_mm512_load_epi64(a0 + lowords * VECLEN),
        _mm512_setzero_si512()) & vec_gte2_52(a1, a0, hiwords)) : vec_gte2_52(a1, a0, hiwords);
    // result is negative if the longer subtractand has a leading digit
    // OR if the rest of the comparison is true.
    m2 = hiwords < lowords ? (_mm512_cmpgt_epu64_mask(_mm512_load_epi64(b0 + lowords * VECLEN),
        _mm512_setzero_si512()) | vec_gte2_52(b0, b1, hiwords)) : vec_gte2_52(b0, b1, hiwords);
    base_abssub_52(a0, a1, s1, m1, lowords, hiwords);
    base_abssub_52(b1, b0, s2, m2, hiwords, lowords);
    veckmul(s1, s2, z1, scratch + (4 * lowords) * VECLEN, lowords);

    kcombine2(c, z1, s1, m1 ^ m2, 2 * lowords, hiwords);
    return;
}

void kmuln_lo(uint64_t* a, uint64_t* b, uint64_t* c, uint64_t* scratch, int words)
{
    int lowords, hiwords, i;

    if (words & 1)
    {
        hiwords = words / 2;
        lowords = hiwords + 1;
    }
    else
    {
        hiwords = lowords = words / 2;
    }

    uint64_t* a0 = a;
    uint64_t* a1 = a + lowords * VECLEN;
    uint64_t* b0 = b;
    uint64_t* b1 = b + lowords * VECLEN;
    uint64_t* z1;
    uint64_t* s1;
    uint64_t* s2;
    __mmask8 m1;
    __mmask8 m2;

    s1 = scratch;
    s2 = scratch + (lowords)*VECLEN;
    z1 = scratch + (lowords * 2) * VECLEN;

    veckmul(a0, b0, c, scratch, lowords);
    veckmul_lo(a1, b1, c + 2 * lowords * VECLEN, scratch, hiwords);
    // result is negative if the longer subtractand does not have a leading digit
    // AND if the rest of the comparison is true.
    m1 = hiwords < lowords ? (_mm512_cmpeq_epu64_mask(_mm512_load_epi64(a0 + lowords * VECLEN),
        _mm512_setzero_si512()) & vec_gte2_52(a1, a0, hiwords)) : vec_gte2_52(a1, a0, hiwords);
    // result is negative if the longer subtractand has a leading digit
    // OR if the rest of the comparison is true.
    m2 = hiwords < lowords ? (_mm512_cmpgt_epu64_mask(_mm512_load_epi64(b0 + lowords * VECLEN),
        _mm512_setzero_si512()) | vec_gte2_52(b0, b1, hiwords)) : vec_gte2_52(b0, b1, hiwords);
    base_abssub_52(a0, a1, s1, m1, lowords, hiwords);
    base_abssub_52(b1, b0, s2, m2, hiwords, lowords);
    veckmul_lo(s1, s2, z1, scratch + (4 * lowords) * VECLEN, lowords);

    kcombine2_lo(c, z1, s1, m1 ^ m2, 2 * lowords, hiwords);
    return;
}

void ksqrn(uint64_t* a, uint64_t* c, uint64_t* scratch, int words)
{
    int lowords, hiwords;

    if (words & 1)
    {
        hiwords = words / 2;
        lowords = hiwords + 1;
    }
    else
    {
        hiwords = lowords = words / 2;
    }
    uint64_t* a0 = a;
    uint64_t* a1 = a + lowords * VECLEN;
    uint64_t* z1;
    uint64_t* s1;
    uint64_t* s2;
    __mmask8 m1;
    __mmask8 m2;
    int i;

    if (a == c)
    {
        printf("error output pointer cannot be the same as input pointer\n");
        exit(1);
    }

    s1 = scratch;
    s2 = scratch + (lowords)*VECLEN;
    z1 = scratch + (lowords * 2) * VECLEN;

    //printf("ksqr lo len %d words: ", lowords);
    //for (i = lowords - 1; i >= 0; i--)
    //{
    //    printf("%013llx", a0[i * VECLEN]);
    //}
    //printf("\n");

    vecksqr(a0, c, scratch, lowords);

    //printf("result: ");
    //for (i = lowords * 2 - 1; i >= 0; i--)
    //{
    //    printf("%013llx", c[i * VECLEN]);
    //}
    //printf("\n");
    //
    //printf("ksqr hi len %d words: ", hiwords);
    //for (i = hiwords - 1; i >= 0; i--)
    //{
    //    printf("%013llx", a1[i * VECLEN]);
    //}
    //printf("\n");

    vecksqr(a1, c + 2 * lowords * VECLEN, scratch, hiwords);

    //printf("result: ");
    //for (i = lowords * 2 - 1; i >= 0; i--)
    //{
    //    printf("%013llx", c[(2 * lowords + i) * VECLEN]);
    //}
    //printf("\n");
    //
    //printf("a0: ");
    //for (i = lowords - 1; i >= 0; i--)
    //{
    //    printf("%013llx", a0[i * VECLEN]);
    //}
    //printf("\n");
    //
    //printf("a1: ");
    //for (i = hiwords - 1; i >= 0; i--)
    //{
    //    printf("%013llx", a1[i * VECLEN]);
    //}
    //printf("\n");

    // result is negative if the longer subtractand does not have a leading digit
    // AND if the rest of the comparison is true.
    m1 = hiwords < lowords ? (_mm512_cmpeq_epu64_mask(_mm512_load_epi64(a0 + lowords * VECLEN),
        _mm512_setzero_si512()) & vec_gte2_52(a1, a0, hiwords)) : vec_gte2_52(a1, a0, hiwords);
    base_abssub_52(a0, a1, s1, m1, lowords, hiwords);

    //printf("ksqr absdiff len %d words: ", lowords);
    //for (i = lowords - 1; i >= 0; i--)
    //{
    //    printf("%013llx", s1[i * VECLEN]);
    //}
    //printf("\n");

    vecksqr(s1, z1, scratch + (4 * lowords) * VECLEN, lowords);

    //printf("combine at offset %d words\n", lowords);
    //
    //printf("z0=0x");
    //for (i = lowords * 2 - 1; i >= 0; i--)
    //{
    //    printf("%013llx", c[i * VECLEN]);
    //}
    //printf(";\n");
    //
    //printf("z1=0x");
    //for (i = lowords * 2 - 1; i >= 0; i--)
    //{
    //    printf("%013llx", z1[i * VECLEN]);
    //}
    //printf(";\n");
    //
    //printf("z2=0x");
    //for (i = lowords * 2 - 1; i >= 0; i--)
    //{
    //    printf("%013llx", c[(2 * lowords + i) * VECLEN]);
    //}
    //printf(";\n");
    //
    //printf("s1=0x");
    //for (i = lowords * 2 - 1; i >= 0; i--)
    //{
    //    printf("%013llx", s1[i * VECLEN]);
    //}
    //printf(";\n");

    kcombine2(c, z1, s1, 0xff, 2 * lowords, hiwords);

    //printf("result: ");
    //for (i = lowords * 4 - 1; i >= 0; i--)
    //{
    //    printf("%013llx", c[i * VECLEN]);
    //}
    //printf("\n");

    return;
}

void vecksqr_mersenne(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    vecksqr(a->data, s->data, mdata->mtmp4->data, mdata->NWORDS);
    vecmod_mersenne(s, c, mdata->mtmp2, mdata);
    return;
}

void veckmul_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    veckmul(a->data, b->data, s->data, mdata->mtmp4->data, mdata->NWORDS);
    vecmod_mersenne(s, c, mdata->mtmp4, mdata);
    return;
}

uint32_t vec_bignum52_special_add(uint64_t* a, uint64_t* b, uint64_t* c, int words)
{
    // assumptions:
    // a, b, c are of length VECLEN * NWORDS
    // s1 is of length VECLEN
    // a, b, c, and s1 are aligned
    // a and b are both positive
    // a >= b
    int i;
    __mmask8 carry = 0xff;  // because the lower half add in karatsuba always carries.
    __m512i avec;
    __m512i bvec;
    __m512i cvec;

    // subtract the selected elements ('1' in the mask)
    for (i = 0; i < words; i++)
    {
        avec = _mm512_load_epi64(a + i * VECLEN);
        bvec = _mm512_load_epi64(b + i * VECLEN);
        cvec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c + i * VECLEN,
            _mm512_and_epi64(_mm512_set1_epi64(VEC_MAXDIGIT), cvec));
    }

    if (carry)
    {
        // subtract any final borrows that exceed the size of b.
        _mm512_mask_store_epi64(c + i * VECLEN, carry, _mm512_set1_epi64(1));
    }

    return carry;
}

void vecksqr_redc(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    vecksqr(a->data, mdata->mtmp1->data, mdata->mtmp4->data, mdata->NWORDS);

    veckmul_lo(mdata->mtmp1->data, mdata->vnhat->data, mdata->mtmp2->data, mdata->mtmp4->data, mdata->NWORDS);
    veckmul(mdata->mtmp2->data, n->data, mdata->mtmp3->data, mdata->mtmp4->data, mdata->NWORDS);

    mdata->mtmp1->size = mdata->NWORDS * 2;
    mdata->mtmp3->size = mdata->NWORDS * 2;
    uint32_t m = vec_bignum52_special_add(mdata->mtmp1->data + mdata->NWORDS * VECLEN,
        mdata->mtmp3->data + mdata->NWORDS * VECLEN, c->data, mdata->NWORDS);

#ifndef USE_AMM
    //m |= vec_gte52(c, n);
    c->WORDS_ALLOC = c->size = mdata->NWORDS;
    //vec_bignum52_mask_sub(c, n, c, m);
#endif
    return;
}

void veckmul_redc(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    veckmul(a->data, b->data, mdata->mtmp1->data, mdata->mtmp4->data, mdata->NWORDS);

    veckmul_lo(mdata->mtmp1->data, mdata->vnhat->data, mdata->mtmp2->data, mdata->mtmp4->data, mdata->NWORDS);
    veckmul(mdata->mtmp2->data, n->data, mdata->mtmp3->data, mdata->mtmp4->data, mdata->NWORDS);
    mdata->mtmp1->size = mdata->NWORDS * 2;
    mdata->mtmp3->size = mdata->NWORDS * 2;
    uint32_t m = vec_bignum52_special_add(mdata->mtmp1->data + mdata->NWORDS * VECLEN,
        mdata->mtmp3->data + mdata->NWORDS * VECLEN, c->data, mdata->NWORDS);
    
#ifndef USE_AMM
    //m |= vec_gte52(c, n);
    c->WORDS_ALLOC = c->size = mdata->NWORDS;
    //vec_bignum52_mask_sub(c, n, c, m);
#endif
    return;
}

void vectmul_redc(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    vectmul(a->data, b->data, mdata->mtmp1->data, mdata->mtmp4->data, mdata->NWORDS);

    vectmul(mdata->mtmp1->data, mdata->vnhat->data, mdata->mtmp2->data, mdata->mtmp4->data, mdata->NWORDS);
    vectmul(mdata->mtmp2->data, n->data, mdata->mtmp3->data, mdata->mtmp4->data, mdata->NWORDS);
    mdata->mtmp1->size = mdata->NWORDS * 2;
    mdata->mtmp3->size = mdata->NWORDS * 2;
    uint32_t m = vec_bignum52_special_add(mdata->mtmp1->data + mdata->NWORDS * VECLEN,
        mdata->mtmp3->data + mdata->NWORDS * VECLEN, c->data, mdata->NWORDS);

#ifndef USE_AMM
    //m |= vec_gte52(c, n);
    c->WORDS_ALLOC = c->size = mdata->NWORDS;
    //vec_bignum52_mask_sub(c, n, c, m);
#endif
    return;
}

void vectsqr_redc(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    vectmul_redc(a, a, c, n, s, mdata);
    return;
}
#endif