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
#include "monty.h"

#ifdef USE_AVX512F
#include <immintrin.h>

//#define PRINT_DEBUG 2
//#define NPRINT_DEBUG 0

void vecmulmod52_mersenne_416(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    int i, j;
    int NWORDS = (int)mdata->NWORDS;
    int NBLOCKS = (int)mdata->NBLOCKS;
    // needed in loops
    __m512i a0, a1;
    UNUSED_VAR __m512i a2, a3;                                   // 4
    __m512i b0;
    UNUSED_VAR __m512i b1, b2, b3, b4, b5, b6;                  // 11
    UNUSED_VAR __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19
    __m512i c00, c01, c02, c03, c04, c05, c06, c07,
        c08, c09, c10, c11, c12, c13, c14, c15;

    UNUSED_VAR __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    UNUSED_VAR __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    UNUSED_VAR __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    UNUSED_VAR __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    UNUSED_VAR __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    UNUSED_VAR __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry;

    // deal with the sign
    c->size = (int)(uint32_t)NWORDS;
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
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 3) * VECLEN_I);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 2) * VECLEN_I);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 1) * VECLEN_I);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 0) * VECLEN_I);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 1) * VECLEN_I);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 2) * VECLEN_I);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 3) * VECLEN_I);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 4) * VECLEN_I);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 5) * VECLEN_I);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 6) * VECLEN_I);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 7) * VECLEN_I);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum
        a0 = _mm512_load_epi64(a->data + (i * BLOCKWORDS_I + 0) * VECLEN_I);
        a1 = _mm512_load_epi64(a->data + (i * BLOCKWORDS_I + 1) * VECLEN_I);
        a2 = _mm512_load_epi64(a->data + (i * BLOCKWORDS_I + 2) * VECLEN_I);
        a3 = _mm512_load_epi64(a->data + (i * BLOCKWORDS_I + 3) * VECLEN_I);

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
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 3) * VECLEN_I);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 2) * VECLEN_I);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 1) * VECLEN_I);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 0) * VECLEN_I);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 1) * VECLEN_I);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 2) * VECLEN_I);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 3) * VECLEN_I);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 4) * VECLEN_I);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 5) * VECLEN_I);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 6) * VECLEN_I);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 7) * VECLEN_I);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum
        a0 = _mm512_load_epi64(a->data + (i * BLOCKWORDS_I + 0) * VECLEN_I);
        a1 = _mm512_load_epi64(a->data + (i * BLOCKWORDS_I + 1) * VECLEN_I);
        a2 = _mm512_load_epi64(a->data + (i * BLOCKWORDS_I + 2) * VECLEN_I);
        a3 = _mm512_load_epi64(a->data + (i * BLOCKWORDS_I + 3) * VECLEN_I);

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
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 3) * VECLEN_I);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 2) * VECLEN_I);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 1) * VECLEN_I);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 0) * VECLEN_I);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 1) * VECLEN_I);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 2) * VECLEN_I);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 3) * VECLEN_I);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 4) * VECLEN_I);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 5) * VECLEN_I);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 6) * VECLEN_I);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 7) * VECLEN_I);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }


        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS_I + 1) * VECLEN_I);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS_I + 2) * VECLEN_I);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS_I + 3) * VECLEN_I);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN_I);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN_I);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN_I);

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
            //_mm512_store_epi64(c->data + (i * BLOCKWORDS_I + 0) * VECLEN_I, a3);
            //_mm512_store_epi64(c->data + (i * BLOCKWORDS_I + 1) * VECLEN_I, a2);
            //_mm512_store_epi64(c->data + (i * BLOCKWORDS_I + 2) * VECLEN_I, a1);
            //_mm512_store_epi64(c->data + (i * BLOCKWORDS_I + 3) * VECLEN_I, a0);

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
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 3) * VECLEN_I);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 2) * VECLEN_I);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 1) * VECLEN_I);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 0) * VECLEN_I);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 1) * VECLEN_I);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 2) * VECLEN_I);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 3) * VECLEN_I);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 4) * VECLEN_I);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 5) * VECLEN_I);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 6) * VECLEN_I);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 7) * VECLEN_I);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }


        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS_I + 1) * VECLEN_I);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS_I + 2) * VECLEN_I);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS_I + 3) * VECLEN_I);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN_I);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN_I);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN_I);

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
            //_mm512_store_epi64(c->data + (i * BLOCKWORDS_I + 0) * VECLEN_I, a3);
            //_mm512_store_epi64(c->data + (i * BLOCKWORDS_I + 1) * VECLEN_I, a2);
            //_mm512_store_epi64(c->data + (i * BLOCKWORDS_I + 2) * VECLEN_I, a1);
            //_mm512_store_epi64(c->data + (i * BLOCKWORDS_I + 3) * VECLEN_I, a0);

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
            (long long)((1ULL << (uint64_t)(bshift_a[0])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[1])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[2])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[3])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[4])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[5])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[6])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[7])) - 1ULL));

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

        vbshift = _mm512_set1_epi64((long long)((1ULL << (uint64_t)(bshift)) - 1ULL));
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
        a1 = _mm512_load_epi64(s->data + wshift * VECLEN_I);
        _mm512_store_epi64(s->data + wshift * VECLEN_I,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN_I, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(s, "cleared low part:");
#endif

        // add c * hi + lo
        scarry = 0;
        for (i = 0; i <= NWORDS; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
            b0 = _mm512_load_epi64(s->data + i * VECLEN_I);
            a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "lo + hi * c:");
#endif

        // now we need to do it again, but the hi part is just a single word so
        // we only have a single-precision vecmul
        //vec_bignum52_mask_rshift_n(c, s, mdata->nbits, 0xff);
        //a0 = _mm512_load_epi64(s->data);

        // It's possible these hi bits are located in two words
        a0 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        a0 = _mm512_srli_epi64(a0, (unsigned int)bshift);
        a1 = _mm512_load_epi64(c->data + (wshift + 1) * VECLEN_I);
        a0 = _mm512_or_epi64(a0, _mm512_slli_epi64(a1, (unsigned int)(52 - bshift)));

        //print_regvechex(a0, 0, "hi part:");

        b0 = _mm512_set1_epi64(mdata->isMersenne);

        //print_regvechex(b0, 0, "c:");

        VEC_MUL_LOHI_PD(a0, b0, acc_e0, acc_e1);

        // clear any hi-bits now that we have the multiplier
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        _mm512_store_epi64(c->data + wshift * VECLEN_I,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_set1_epi64(0));
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
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
            a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
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
                a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
                a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
                _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
                i++;
            }
        }

        // finally clear any hi-bits in the result
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        _mm512_store_epi64(c->data + wshift * VECLEN_I,
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
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
            b0 = _mm512_load_epi64(s->data + i * VECLEN_I);
            a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
        }

        //c00 = _mm512_adc_epi52(c00, scarry, c08, &scarry);
        //c01 = _mm512_adc_epi52(c01, scarry, c09, &scarry);
        //c02 = _mm512_adc_epi52(c02, scarry, c10, &scarry);
        //c03 = _mm512_adc_epi52(c03, scarry, c11, &scarry);
        //c04 = _mm512_adc_epi52(c04, scarry, c12, &scarry);
        //c05 = _mm512_adc_epi52(c05, scarry, c13, &scarry);
        //c06 = _mm512_adc_epi52(c06, scarry, c14, &scarry);
        //c07 = _mm512_adc_epi52(c07, scarry, c15, &scarry);

        a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
        b0 = _mm512_load_epi64(s->data + i * VECLEN_I);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN_I, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add it back in.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);

        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
            a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
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



        // clear the potential hi-bit
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        _mm512_store_epi64(c->data + wshift * VECLEN_I,
            _mm512_and_epi64(vbshift, a1));

        //c07 = _mm512_and_epi64(vbshift, c07);

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_set1_epi64(0));
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
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);   // hi
            b0 = _mm512_load_epi64(s->data + i * VECLEN_I);   // lo
            a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
        b0 = _mm512_load_epi64(s->data + i * VECLEN_I);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN_I, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add 1.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
            a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // and add one to the hi-bit.  This should resolve the borrow bit.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        //a1 = _mm512_add_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift))), a1);
        //_mm512_store_epi64(c->data + wshift * VECLEN_I,
        //    _mm512_and_epi64(_mm512_set1_epi64(0xfffffffffffffULL), a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif

        _mm512_store_epi64(c->data + wshift * VECLEN_I,
            _mm512_and_epi64(vbshift, a1));
    }

    c->size = (int)(uint32_t)NWORDS;
    return;
}

void vecsqrmod52_mersenne_416(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    // 8x sqr:
    // input 8 bignums in the even lanes of a.
    // output 8 squaremod bignums in the even lanes of c.
    int i, j;
    UNUSED_VAR int k;
    int NWORDS = (int)mdata->NWORDS;
    int NBLOCKS = (int)mdata->NBLOCKS;
    vec_bignum_t* b = a;

    // needed in loops
    __m512i a0, a1;
    UNUSED_VAR __m512i a2, a3;                                   // 4
    __m512i b0;
    UNUSED_VAR __m512i b1, b2, b3, b4, b5, b6;                  // 11
    UNUSED_VAR __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19
    __m512i c00, c01, c02, c03, c04, c05, c06, c07,
        c08, c09, c10, c11, c12, c13, c14, c15;

#ifndef IFMA
    UNUSED_VAR __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    UNUSED_VAR __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    UNUSED_VAR __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    UNUSED_VAR __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    UNUSED_VAR __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    UNUSED_VAR __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
#endif

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry;

    // deal with the sign
    c->size = (int)(uint32_t)NWORDS;
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
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 3) * VECLEN_I);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 2) * VECLEN_I);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 1) * VECLEN_I);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 0) * VECLEN_I);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 1) * VECLEN_I);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 2) * VECLEN_I);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 3) * VECLEN_I);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 4) * VECLEN_I);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 5) * VECLEN_I);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 6) * VECLEN_I);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 7) * VECLEN_I);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }


        {
            // i even
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 0) * VECLEN_I);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 1) * VECLEN_I);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 2) * VECLEN_I);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 3) * VECLEN_I);

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
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 3) * VECLEN_I);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 2) * VECLEN_I);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 1) * VECLEN_I);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 0) * VECLEN_I);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 1) * VECLEN_I);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 2) * VECLEN_I);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 3) * VECLEN_I);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 4) * VECLEN_I);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 5) * VECLEN_I);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 6) * VECLEN_I);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 7) * VECLEN_I);

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

            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 3) * VECLEN_I);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 2) * VECLEN_I);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 1) * VECLEN_I);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 0) * VECLEN_I);

            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 2) * VECLEN_I);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 3) * VECLEN_I);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 4) * VECLEN_I);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 5) * VECLEN_I);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 6) * VECLEN_I);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 7) * VECLEN_I);


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
            a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS_I) * VECLEN_I);
            a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS_I) * VECLEN_I);
            a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS_I) * VECLEN_I);
            a3 = _mm512_load_epi64(a->data + (NWORDS - 4 - j * BLOCKWORDS_I) * VECLEN_I);

            b0 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS_I + 1) * VECLEN_I);
            b1 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS_I + 2) * VECLEN_I);
            b2 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS_I + 3) * VECLEN_I);
            b3 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS_I + 4) * VECLEN_I);
            b4 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS_I + 5) * VECLEN_I);
            b5 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS_I + 6) * VECLEN_I);
            b6 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS_I + 7) * VECLEN_I);

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
                a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS_I) * VECLEN_I);		// {f, b}
                a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS_I) * VECLEN_I);		// {e, a}
                a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS_I) * VECLEN_I);		// {d, 9}
                a3 = _mm512_load_epi64(a->data + (NWORDS - 4 - j * BLOCKWORDS_I) * VECLEN_I);		// {c, 8}

                b0 = _mm512_load_epi64(b->data + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 1) * VECLEN_I); // {9, 5}
                b1 = _mm512_load_epi64(b->data + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 2) * VECLEN_I);	// {a, 6}
                b2 = _mm512_load_epi64(b->data + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 3) * VECLEN_I);	// {b, 7}
                b3 = _mm512_load_epi64(b->data + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 4) * VECLEN_I);	// {c, 8}
                b4 = _mm512_load_epi64(b->data + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 5) * VECLEN_I); // {d, 9}

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
            a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS_I) * VECLEN_I);
            a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS_I) * VECLEN_I);
            a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS_I) * VECLEN_I);
            a3 = _mm512_load_epi64(a->data + (NWORDS - 4 - j * BLOCKWORDS_I) * VECLEN_I);

            b0 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS_I + 1) * VECLEN_I);
            b1 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS_I + 2) * VECLEN_I);
            b2 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS_I + 3) * VECLEN_I);
            b3 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS_I + 4) * VECLEN_I);
            b4 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS_I + 5) * VECLEN_I);
            b5 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS_I + 6) * VECLEN_I);
            b6 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS_I + 7) * VECLEN_I);

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
                a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS_I) * VECLEN_I);  // {f, b}
                a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS_I) * VECLEN_I);  // {e, a}
                a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS_I) * VECLEN_I);  // {d, 9}

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
            (long long)((1ULL << (uint64_t)(bshift_a[0])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[1])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[2])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[3])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[4])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[5])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[6])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[7])) - 1ULL));

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

        vbshift = _mm512_set1_epi64((long long)((1ULL << (uint64_t)(bshift)) - 1ULL));
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
        a1 = _mm512_load_epi64(s->data + wshift * VECLEN_I);
        _mm512_store_epi64(s->data + wshift * VECLEN_I,
            _mm512_and_epi64(_mm512_set1_epi64((long long)((1ULL << (uint64_t)(bshift)) - 1ULL)), a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN_I, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(s, "cleared low part:");
#endif

        // add c * hi + lo
        scarry = 0;
        for (i = 0; i <= NWORDS; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
            b0 = _mm512_load_epi64(s->data + i * VECLEN_I);
            a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "lo + hi * c:");
#endif

        // now we need to do it again, but the hi part is just a single word so
        // we only have a single-precision vecmul.  It's possible these hi bits
        // are located in two words
        //vec_bignum52_mask_rshift_n(c, s, mdata->nbits, 0xff);
        //a0 = _mm512_load_epi64(s->data);

        a0 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        a0 = _mm512_srli_epi64(a0, (unsigned int)bshift);
        a1 = _mm512_load_epi64(c->data + (wshift + 1) * VECLEN_I);
        a0 = _mm512_or_epi64(a0, _mm512_slli_epi64(a1, (unsigned int)(52 - bshift)));

        b0 = _mm512_set1_epi64(mdata->isMersenne);
        VEC_MUL_LOHI_PD(a0, b0, acc_e0, acc_e1);

        // clear any hi-bits now that we have the multiplier
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        _mm512_store_epi64(c->data + wshift * VECLEN_I,
            _mm512_and_epi64(_mm512_set1_epi64((long long)((1ULL << (uint64_t)(bshift)) - 1ULL)), a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_set1_epi64(0));
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
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
            a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
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
                a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
                a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
                _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
                i++;
            }
        }

        // finally clear any hi-bits in the result
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        _mm512_store_epi64(c->data + wshift * VECLEN_I,
            _mm512_and_epi64(_mm512_set1_epi64((long long)((1ULL << (uint64_t)(bshift)) - 1ULL)), a1));

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
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
            b0 = _mm512_load_epi64(s->data + i * VECLEN_I);
            a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
        b0 = _mm512_load_epi64(s->data + i * VECLEN_I);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN_I, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add it back in.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
            a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // clear the potential hi-bit
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        _mm512_store_epi64(c->data + wshift * VECLEN_I,
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
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);   // hi
            b0 = _mm512_load_epi64(s->data + i * VECLEN_I);   // lo
            a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
        b0 = _mm512_load_epi64(s->data + i * VECLEN_I);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));

        // zero any remaining lo words.
        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN_I, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after sub:");
#endif

        // if there was a carry, add 1 to the lo bit.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
            a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // and add one to the hi-bit.  This should resolve the borrow bit.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        //a1 = _mm512_add_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift))), a1);
        //_mm512_store_epi64(c->data + wshift * VECLEN_I,
        //    _mm512_and_epi64(_mm512_set1_epi64(0xfffffffffffffULL), a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
#endif

        _mm512_store_epi64(c->data + wshift * VECLEN_I,
            _mm512_and_epi64(vbshift, a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after bit clear:");
        exit(1);
#endif

    }

    c->size = (int)(uint32_t)NWORDS;
    return;
}

void vecmulmod52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    int i, j;
    int NWORDS = (int)mdata->NWORDS;
    int NBLOCKS = (int)mdata->NBLOCKS;

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
    UNUSED_VAR __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    UNUSED_VAR __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    UNUSED_VAR __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    UNUSED_VAR __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    UNUSED_VAR __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    UNUSED_VAR __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
#endif

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry;

    // deal with the sign
    c->size = (int)(uint32_t)NWORDS;
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
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 3) * VECLEN_I);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 2) * VECLEN_I);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 1) * VECLEN_I);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 0) * VECLEN_I);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 1) * VECLEN_I);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 2) * VECLEN_I);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 3) * VECLEN_I);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 4) * VECLEN_I);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 5) * VECLEN_I);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 6) * VECLEN_I);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 7) * VECLEN_I);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum
        a0 = _mm512_load_epi64(a->data + (i * BLOCKWORDS_I + 0) * VECLEN_I);
        a1 = _mm512_load_epi64(a->data + (i * BLOCKWORDS_I + 1) * VECLEN_I);
        a2 = _mm512_load_epi64(a->data + (i * BLOCKWORDS_I + 2) * VECLEN_I);
        a3 = _mm512_load_epi64(a->data + (i * BLOCKWORDS_I + 3) * VECLEN_I);

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
            _mm512_store_epi64(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I, acc_e0);

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
            _mm512_store_epi64(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I, acc_e0);

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
            _mm512_store_epi64(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I, acc_e0);

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
            _mm512_store_epi64(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I, acc_e0);

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
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 3) * VECLEN_I);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 2) * VECLEN_I);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 1) * VECLEN_I);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 0) * VECLEN_I);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 1) * VECLEN_I);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 2) * VECLEN_I);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 3) * VECLEN_I);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 4) * VECLEN_I);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 5) * VECLEN_I);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 6) * VECLEN_I);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 7) * VECLEN_I);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }


        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS_I + 1) * VECLEN_I);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS_I + 2) * VECLEN_I);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS_I + 3) * VECLEN_I);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN_I);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN_I);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN_I);

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
            //_mm512_store_pd(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I, 
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
            //_mm512_store_pd(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I,
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
            //_mm512_store_pd(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I,
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
            //_mm512_store_pd(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            a3 = _mm512_and_epi64(vlmask, a3);
            a2 = _mm512_and_epi64(vlmask, a2);
            a1 = _mm512_and_epi64(vlmask, a1);
            a0 = _mm512_and_epi64(vlmask, a0);

            _mm512_store_epi64(s->data + (i * BLOCKWORDS_I + 0) * VECLEN_I, a3);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS_I + 1) * VECLEN_I, a2);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS_I + 2) * VECLEN_I, a1);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS_I + 3) * VECLEN_I, a0);

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
            (long long)((1ULL << (uint64_t)(bshift_a[0])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[1])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[2])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[3])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[4])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[5])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[6])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[7])) - 1ULL));

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
        vbshift = _mm512_set1_epi64((long long)((1ULL << (uint64_t)(bshift)) - 1ULL));
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
        a1 = _mm512_load_epi64(s->data + wshift * VECLEN_I);
        _mm512_store_epi64(s->data + wshift * VECLEN_I,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN_I, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(s, "cleared low part:");
#endif

        // add c * hi + lo
        scarry = 0;
        for (i = 0; i <= NWORDS; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
            b0 = _mm512_load_epi64(s->data + i * VECLEN_I);
            
            //a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            __m512i t = _mm512_add_epi64(a1, b0);
            t = _mm512_add_epi64(t, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpgt_epu64_mask(t, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "lo + hi * c:");
#endif

        // now we need to do it again, but the hi part is just a single word so
        // we only have a single-precision vecmul
        //vec_bignum52_mask_rshift_n(c, s, mdata->nbits, 0xff);
        //a0 = _mm512_load_epi64(s->data);

        // It's possible these hi bits are located in two words
        a0 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        a0 = _mm512_srli_epi64(a0, (unsigned int)bshift);
        a1 = _mm512_load_epi64(c->data + (wshift + 1) * VECLEN_I);
        a0 = _mm512_or_epi64(a0, _mm512_slli_epi64(a1, (unsigned int)(52 - bshift)));

        //print_regvechex(a0, 0, "hi part:");

        b0 = _mm512_set1_epi64(mdata->isMersenne);

        //print_regvechex(b0, 0, "c:");

        VEC_MUL_LOHI_PD(a0, b0, acc_e0, acc_e1);

        // clear any hi-bits now that we have the multiplier
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        _mm512_store_epi64(c->data + wshift * VECLEN_I,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_set1_epi64(0));
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
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
            
            //a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            __m512i t = _mm512_add_epi64(a1, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpeq_epu64_mask(a1, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
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
                a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
                
                //a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
                __m512i t = _mm512_add_epi64(a1, _mm512_maskz_set1_epi64(scarry, 1));
                scarry = _mm512_cmpeq_epu64_mask(a1, vlmask);
                a0 = _mm512_and_epi64(t, vlmask);

                _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
                i++;
            }
        }

        // finally clear any hi-bits in the result
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        _mm512_store_epi64(c->data + wshift * VECLEN_I,
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
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
            b0 = _mm512_load_epi64(s->data + i * VECLEN_I);
            
            //a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            __m512i t = _mm512_add_epi64(a1, b0);
            t = _mm512_add_epi64(t, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpgt_epu64_mask(t, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
        b0 = _mm512_load_epi64(s->data + i * VECLEN_I);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN_I, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add it back in.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
            
            //a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            __m512i t = _mm512_add_epi64(a1, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpeq_epu64_mask(a1, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // clear the potential hi-bit
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        _mm512_store_epi64(c->data + wshift * VECLEN_I,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_set1_epi64(0));
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
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);   // hi
            b0 = _mm512_load_epi64(s->data + i * VECLEN_I);   // lo
            
            //a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
            __m512i t = _mm512_sub_epi64(b0, a1);
            carryout = _mm512_cmpgt_epu64_mask(a1, b0);
            __m512i t2 = _mm512_sub_epi64(t, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = (__mmask8)_mm512_kor(carryout, _mm512_cmpgt_epu64_mask(t2, t));
            a0 = _mm512_and_epi64(t2, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
        b0 = _mm512_load_epi64(s->data + i * VECLEN_I);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN_I, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add 1.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
            
            //a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            __m512i t = _mm512_add_epi64(a1, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpeq_epu64_mask(a1, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // and add one to the hi-bit.  This should resolve the borrow bit.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        //a1 = _mm512_add_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift))), a1);
        //_mm512_store_epi64(c->data + wshift * VECLEN_I,
        //    _mm512_and_epi64(_mm512_set1_epi64(0xfffffffffffffULL), a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif

        _mm512_store_epi64(c->data + wshift * VECLEN_I,
            _mm512_and_epi64(vbshift, a1));
    }

#ifdef DEBUG_MERSENNE
    mpz_clear(ga);
    mpz_clear(gb);
#endif

    c->size = (int)(uint32_t)NWORDS;
    return;
}

void vecsqrmod52_mersenne(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    vecmulmod52_mersenne(a, a, c, n, s, mdata);
    return;

#ifdef USE_NATIVE_SQRMOD_MERSENNE
    // 8x sqr:
    // input 8 bignums in the even lanes of a.
    // output 8 squaremod bignums in the even lanes of c.
    int i, j;
    int k;
    int NWORDS = (int)mdata->NWORDS;
    int NBLOCKS = (int)mdata->NBLOCKS;
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
    UNUSED_VAR __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    UNUSED_VAR __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    UNUSED_VAR __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    UNUSED_VAR __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    UNUSED_VAR __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    UNUSED_VAR __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
#endif

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry;

    // deal with the sign
    c->size = (uint32_t)NWORDS;
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
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 3) * VECLEN_I);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 2) * VECLEN_I);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 1) * VECLEN_I);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 0) * VECLEN_I);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 1) * VECLEN_I);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 2) * VECLEN_I);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 3) * VECLEN_I);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 4) * VECLEN_I);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 5) * VECLEN_I);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 6) * VECLEN_I);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 7) * VECLEN_I);

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

            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 3) * VECLEN_I);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 2) * VECLEN_I);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 1) * VECLEN_I);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 0) * VECLEN_I);

            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 2) * VECLEN_I);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 3) * VECLEN_I);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 4) * VECLEN_I);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 5) * VECLEN_I);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 6) * VECLEN_I);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS_I + 7) * VECLEN_I);


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
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 0) * VECLEN_I);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 1) * VECLEN_I);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 2) * VECLEN_I);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS_I + 3) * VECLEN_I);

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
            _mm512_store_epi64(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I, acc_e0);

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
            _mm512_store_epi64(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I, acc_e0);

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
            _mm512_store_epi64(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I, acc_e0);

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
            _mm512_store_epi64(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I, acc_e0);

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
            a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS_I) * VECLEN_I);
            a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS_I) * VECLEN_I);
            a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS_I) * VECLEN_I);
            a3 = _mm512_load_epi64(a->data + (NWORDS - 4 - j * BLOCKWORDS_I) * VECLEN_I);

            b0 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS_I + 1) * VECLEN_I);
            b1 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS_I + 2) * VECLEN_I);
            b2 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS_I + 3) * VECLEN_I);
            b3 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS_I + 4) * VECLEN_I);
            b4 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS_I + 5) * VECLEN_I);
            b5 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS_I + 6) * VECLEN_I);
            b6 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS_I + 7) * VECLEN_I);

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
                a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS_I) * VECLEN_I);
                a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS_I) * VECLEN_I);
                a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS_I) * VECLEN_I);
                a3 = _mm512_load_epi64(a->data + (NWORDS - 4 - j * BLOCKWORDS_I) * VECLEN_I);

                b0 = _mm512_load_epi64(b->data + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 1) * VECLEN_I);
                b1 = _mm512_load_epi64(b->data + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 2) * VECLEN_I);
                b2 = _mm512_load_epi64(b->data + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 3) * VECLEN_I);
                b3 = _mm512_load_epi64(b->data + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 4) * VECLEN_I);
                b4 = _mm512_load_epi64(b->data + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 5) * VECLEN_I);

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
                a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS_I) * VECLEN_I);
                a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS_I) * VECLEN_I);
                a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS_I) * VECLEN_I);

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
                a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS_I) * VECLEN_I);  // {f, b}
                a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS_I) * VECLEN_I);  // {e, a}
                a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS_I) * VECLEN_I);  // {d, 9}

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
                a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS_I) * VECLEN_I);		// {f, b}
                a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS_I) * VECLEN_I);		// {e, a}
                a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS_I) * VECLEN_I);		// {d, 9}
                a3 = _mm512_load_epi64(a->data + (NWORDS - 4 - j * BLOCKWORDS_I) * VECLEN_I);		// {c, 8}

                b0 = _mm512_load_epi64(b->data + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 1) * VECLEN_I); // {9, 5}
                b1 = _mm512_load_epi64(b->data + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 2) * VECLEN_I);	// {a, 6}
                b2 = _mm512_load_epi64(b->data + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 3) * VECLEN_I);	// {b, 7}
                b3 = _mm512_load_epi64(b->data + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 4) * VECLEN_I);	// {c, 8}
                b4 = _mm512_load_epi64(b->data + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 5) * VECLEN_I); // {d, 9}

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
            //_mm512_store_pd(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I, 
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
            //_mm512_store_pd(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I,
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
            //_mm512_store_pd(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I,
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
            //_mm512_store_pd(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            a3 = _mm512_and_epi64(vlmask, a3);
            a2 = _mm512_and_epi64(vlmask, a2);
            a1 = _mm512_and_epi64(vlmask, a1);
            a0 = _mm512_and_epi64(vlmask, a0);

            _mm512_store_epi64(s->data + (i * BLOCKWORDS_I + NWORDS + 0) * VECLEN_I, a3);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS_I + NWORDS + 1) * VECLEN_I, a2);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS_I + NWORDS + 2) * VECLEN_I, a1);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS_I + NWORDS + 3) * VECLEN_I, a0);

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
        a1 = _mm512_load_epi64(s->data + wshift * VECLEN_I);
        _mm512_store_epi64(s->data + wshift * VECLEN_I,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN_I, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(s, "cleared low part:");
#endif

        // add c * hi + lo
        scarry = 0;
        for (i = 0; i <= NWORDS; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
            b0 = _mm512_load_epi64(s->data + i * VECLEN_I);

            //a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            __m512i t = _mm512_add_epi64(a1, b0);
            t = _mm512_add_epi64(t, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpgt_epu64_mask(t, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "lo + hi * c:");
#endif

        // now we need to do it again, but the hi part is just a single word so
        // we only have a single-precision vecmul
        //vec_bignum52_mask_rshift_n(c, s, mdata->nbits, 0xff);
        //a0 = _mm512_load_epi64(s->data);

        // It's possible these hi bits are located in two words
        a0 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        a0 = _mm512_srli_epi64(a0, bshift);
        a1 = _mm512_load_epi64(c->data + (wshift + 1) * VECLEN_I);
        a0 = _mm512_or_epi64(a0, _mm512_slli_epi64(a1, 52 - bshift));

        //print_regvechex(a0, 0, "hi part:");

        b0 = _mm512_set1_epi64(mdata->isMersenne);

        //print_regvechex(b0, 0, "c:");

        VEC_MUL_LOHI_PD(a0, b0, acc_e0, acc_e1);

        // clear any hi-bits now that we have the multiplier
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        _mm512_store_epi64(c->data + wshift * VECLEN_I,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_set1_epi64(0));
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
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);

            //a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            __m512i t = _mm512_add_epi64(a1, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpeq_epu64_mask(a1, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
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
                a1 = _mm512_load_epi64(c->data + i * VECLEN_I);

                //a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
                __m512i t = _mm512_add_epi64(a1, _mm512_maskz_set1_epi64(scarry, 1));
                scarry = _mm512_cmpeq_epu64_mask(a1, vlmask);
                a0 = _mm512_and_epi64(t, vlmask);

                _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
                i++;
            }
        }

        // finally clear any hi-bits in the result
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        _mm512_store_epi64(c->data + wshift * VECLEN_I,
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
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
            b0 = _mm512_load_epi64(s->data + i * VECLEN_I);

            //a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            __m512i t = _mm512_add_epi64(a1, b0);
            t = _mm512_add_epi64(t, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpgt_epu64_mask(t, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
        b0 = _mm512_load_epi64(s->data + i * VECLEN_I);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN_I, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add it back in.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);

            //a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            __m512i t = _mm512_add_epi64(a1, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpeq_epu64_mask(a1, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // clear the potential hi-bit
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        _mm512_store_epi64(c->data + wshift * VECLEN_I,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_set1_epi64(0));
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
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);   // hi
            b0 = _mm512_load_epi64(s->data + i * VECLEN_I);   // lo

            //a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
            __m512i t = _mm512_sub_epi64(b0, a1);
            carryout = _mm512_cmpgt_epu64_mask(a1, b0);
            __m512i t2 = _mm512_sub_epi64(t, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = (__mmask8)_mm512_kor(carryout, _mm512_cmpgt_epu64_mask(t2, t));
            a0 = _mm512_and_epi64(t2, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
        b0 = _mm512_load_epi64(s->data + i * VECLEN_I);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN_I, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add 1.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);

            //a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            __m512i t = _mm512_add_epi64(a1, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpeq_epu64_mask(a1, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // and add one to the hi-bit.  This should resolve the borrow bit.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        //a1 = _mm512_add_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift))), a1);
        //_mm512_store_epi64(c->data + wshift * VECLEN_I,
        //    _mm512_and_epi64(_mm512_set1_epi64(0xfffffffffffffULL), a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif

        _mm512_store_epi64(c->data + wshift * VECLEN_I,
            _mm512_and_epi64(vbshift, a1));
    }

    c->size = (uint32_t)NWORDS;
    return;
#endif /* USE_NATIVE_SQRMOD_MERSENNE */
}

void vecaddmod52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata)
{
    // assumptions:
    // a, b, c are of length VECLEN_I * NWORDS
    // a, b, c, and n are aligned
    // a and b are both positive
    // n is the montgomery base
    int i;
    int NWORDS = (int)a->WORDS_ALLOC;
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
            (long long)((1ULL << (uint64_t)(bshift_a[0])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[1])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[2])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[3])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[4])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[5])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[6])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[7])) - 1ULL));

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
        vbshift = _mm512_set1_epi64((long long)((1ULL << (uint64_t)(bshift)) - 1ULL));
        vbpshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)));
    }

    // add
    carry = 0;
    for (i = 0; i < wshift; i++)
    {
        avec = _mm512_load_epi64(b->data + i * VECLEN_I);
        bvec = _mm512_load_epi64(a->data + i * VECLEN_I);
        avec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(lomask, avec));
    }

    avec = _mm512_load_epi64(b->data + i * VECLEN_I);
    bvec = _mm512_load_epi64(a->data + i * VECLEN_I);
    bvec = _mm512_and_epi64(vbshift, bvec);
    avec = _mm512_adc_epi52(avec, carry, bvec, &carry);
    _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(lomask, avec));

    for (i++; i < NWORDS; i++)
    {
        _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_set1_epi64(0));
    }

    // check for a carry.
    avec = _mm512_load_epi64(c->data + wshift * VECLEN_I);
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
            cvec = _mm512_load_epi64(c->data + i * VECLEN_I);

            //bvec = _mm512_addcarry_epi52(cvec, carry, &carry);
            __m512i t = _mm512_add_epi64(cvec, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpeq_epu64_mask(cvec, lomask);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(c->data + i * VECLEN_I, bvec);
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
            cvec = _mm512_load_epi64(c->data + i * VECLEN_I);
            
            //bvec = _mm512_subborrow_epi52(cvec, carry, &carry);
            __m512i t = _mm512_sub_epi64(cvec, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpeq_epu64_mask(cvec, zero);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(c->data + i * VECLEN_I, bvec);
        }
    }

    // clear the potential hi-bit
    avec = _mm512_load_epi64(c->data + wshift * VECLEN_I);
    _mm512_store_epi64(c->data + wshift * VECLEN_I,
        _mm512_and_epi64(vbshift, avec));

    return;
}

void vecsubmod52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata)
{
    // assumptions:
    // a, b, c are of length VECLEN_I * NWORDS
    // s1 is of length VECLEN_I
    // a, b, c, n, and s1 are aligned
    // a and b are both positive
    // a >= b
    // n is the montgomery base
    int i;
    __mmask8 carry = 0;
    __mmask8 carryout = 0;
    __mmask8 mask = 0;
    UNUSED_VAR __mmask8 mask2 = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    __m512i nvec;
    __m512i lomask = _mm512_set1_epi64(0xfffffffffffffULL);
    UNUSED_VAR __m512i zero = _mm512_set1_epi64(0);

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
            (long long)((1ULL << (uint64_t)(bshift_a[0])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[1])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[2])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[3])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[4])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[5])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[6])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[7])) - 1ULL));
    }
    else
    {
        vbshift = _mm512_set1_epi64((long long)((1ULL << (uint64_t)(bshift)) - 1ULL));
    }

    // subtract
    carry = 0;
    //for (i = 0; i < NWORDS; i++)
    for (i = 0; i <= wshift; i++)
    {
        avec = _mm512_load_epi64(a->data + i * VECLEN_I);
        bvec = _mm512_load_epi64(b->data + i * VECLEN_I);

        //cvec = _mm512_sbb_epi52(avec, carry, bvec, &carry);
        __m512i t = _mm512_sub_epi64(avec, bvec);
        carryout = _mm512_cmpgt_epu64_mask(bvec, avec);
        __m512i t2 = _mm512_sub_epi64(t, _mm512_maskz_set1_epi64(carry, 1));
        carry = (__mmask8)_mm512_kor(carryout, _mm512_cmpgt_epu64_mask(t2, t));
        cvec = _mm512_and_epi64(t2, lomask);

        _mm512_store_epi64(c->data + i * VECLEN_I, cvec);
    }

    if (mdata->isMersenne > 0)
    {
        // if we had a final carry, then b was bigger than a so we need to re-add n.
        // save the borrow now: every limb needs to add n's corresponding word,
        // and the limb-0 adc below overwrites carry with its own overflow.
        mask = carry;
        nvec = _mm512_set1_epi64(0x10000000000000ULL);
        nvec = _mm512_sub_epi64(nvec, _mm512_set1_epi64(mdata->isMersenne));
        cvec = _mm512_load_epi64(c->data + 0 * VECLEN);
        bvec = _mm512_mask_adc_epi52(cvec, carry, 0, nvec, &carry);
        _mm512_store_epi64(c->data + 0 * VECLEN, bvec);
        nvec = lomask;
        for (i = 1; i <= wshift; i++)
        {
            cvec = _mm512_load_epi64(c->data + i * VECLEN_I);

            //bvec = _mm512_mask_adc_epi52(cvec, mask, carry, nvec, &carry);
            __m512i t = _mm512_add_epi64(cvec, nvec);
            t = _mm512_mask_add_epi64(cvec, mask, t, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpgt_epu64_mask(t, lomask);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(c->data + i * VECLEN_I, bvec);
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
            cvec = _mm512_load_epi64(c->data + i * VECLEN_I);

            //bvec = _mm512_addcarry_epi52(cvec, carry, &carry);
            __m512i t = _mm512_add_epi64(cvec, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpeq_epu64_mask(cvec, lomask);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(c->data + i * VECLEN_I, bvec);
        }
    }

    //for (; i < NWORDS; i++)
    //{
    //    _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_set1_epi64(0));
    //}

    // clear the potential hi-bit
    avec = _mm512_load_epi64(c->data + wshift * VECLEN_I);
    _mm512_store_epi64(c->data + wshift * VECLEN_I,
        _mm512_and_epi64(vbshift, avec));

    return;
}

void vecsubmod52_mersenne2(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata)
{
    // assumptions:
    // a, b, c are of length VECLEN_I * NWORDS
    // s1 is of length VECLEN_I
    // a, b, c, n, and s1 are aligned
    // a and b are both positive
    // a >= b
    // n is the montgomery base
    int i;
    int NWORDS = (int)a->WORDS_ALLOC;
    __mmask8 carry = 0;
    __mmask8 mask = 0;
    UNUSED_VAR __mmask8 mask2 = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    __m512i nvec;
    __m512i lomask = _mm512_set1_epi64(0xfffffffffffffULL);
    UNUSED_VAR __m512i zero = _mm512_set1_epi64(0);

    __m512i vbshift;
    __m512i vbpshift;
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
            (long long)((1ULL << (uint64_t)(bshift_a[0])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[1])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[2])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[3])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[4])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[5])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[6])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[7])) - 1ULL));
    }
    else
    {
        vbshift = _mm512_set1_epi64((long long)((1ULL << (uint64_t)(bshift)) - 1ULL));
        vbpshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)));
    }

    // subtract
    carry = 0;
    for (i = 0; i < wshift; i++)
    {
        avec = _mm512_load_epi64(b->data + i * VECLEN_I);
        bvec = _mm512_load_epi64(a->data + i * VECLEN_I);
        avec = _mm512_sbb_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(lomask, avec));
    }

    avec = _mm512_load_epi64(b->data + i * VECLEN_I);
    bvec = _mm512_load_epi64(a->data + i * VECLEN_I);
    bvec = _mm512_and_epi64(vbshift, bvec);
    avec = _mm512_sbb_epi52(avec, carry, bvec, &carry);
    _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(lomask, avec));

    for (i++; i < NWORDS; i++)
    {
        _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_set1_epi64(0));
    }

    // check for a carry.
    avec = _mm512_load_epi64(c->data + wshift * VECLEN_I);
    carry = _mm512_test_epi64_mask(avec, vbpshift);

    if (mdata->isMersenne > 0)
    {
        // if we had a final carry, then b was bigger than a so we need to re-add n.
        // save the borrow now: every limb needs to add n's corresponding word,
        // and the limb-0 adc below overwrites carry with its own overflow.
        mask = carry;
        nvec = _mm512_set1_epi64(0x10000000000000ULL);
        nvec = _mm512_sub_epi64(nvec, _mm512_set1_epi64(mdata->isMersenne));
        cvec = _mm512_load_epi64(c->data + 0 * VECLEN);
        bvec = _mm512_mask_adc_epi52(cvec, carry, 0, nvec, &carry);
        _mm512_store_epi64(c->data + 0 * VECLEN, bvec);
        nvec = lomask;
        for (i = 1; i <= wshift; i++)
        {
            cvec = _mm512_load_epi64(c->data + i * VECLEN_I);
            
            //bvec = _mm512_mask_adc_epi52(cvec, mask, carry, nvec, &carry);
            __m512i t = _mm512_add_epi64(cvec, nvec);
            t = _mm512_mask_add_epi64(cvec, mask, t, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpgt_epu64_mask(t, lomask);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(c->data + i * VECLEN_I, bvec);
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
            cvec = _mm512_load_epi64(c->data + i * VECLEN_I);
            
            //bvec = _mm512_addcarry_epi52(cvec, carry, &carry);
            __m512i t = _mm512_add_epi64(cvec, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpeq_epu64_mask(cvec, lomask);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(c->data + i * VECLEN_I, bvec);
        }
    }

    //for (; i < NWORDS; i++)
    //{
    //    _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_set1_epi64(0));
    //}

    // clear the potential hi-bit
    avec = _mm512_load_epi64(c->data + wshift * VECLEN_I);
    _mm512_store_epi64(c->data + wshift * VECLEN_I,
        _mm512_and_epi64(vbshift, avec));

    return;
}

void vec_simul_addsub52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* sum, vec_bignum_t* diff,
    vec_monty_t* mdata)
{
    //vecaddmod52_mersenne(a, b, sum, mdata);
    //vecsubmod52_mersenne(a, b, diff, mdata);
    //return;

    // assumptions:
    // a, b, c are of length VECLEN_I * NWORDS
    // a, b, c, and n are aligned
    // a and b are both positive
    // n is the montgomery base
    // produce sum = a + b and diff = a - b at the same time which
    // saves 3N loads (only have to load a,b, and n once)
    int i;
    int NWORDS = (int)a->WORDS_ALLOC;
    __mmask8 carry = 0;
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
            (long long)((1ULL << (uint64_t)(bshift_a[0])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[1])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[2])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[3])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[4])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[5])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[6])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[7])) - 1ULL));

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
        vbshift = _mm512_set1_epi64((long long)((1ULL << (uint64_t)(bshift)) - 1ULL));
        vbpshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)));
    }


    //for (i = 0; i < NWORDS; i++)
    for (i = 0; i <= wshift; i++)
    {
        // add/sub
        avec = _mm512_load_epi64(a->data + i * VECLEN_I);
        bvec = _mm512_load_epi64(b->data + i * VECLEN_I);
        
        //cvec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        //cvec = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
        _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, avec, bvec);
        
        _mm512_store_epi64(sum->data + i * VECLEN_I, avec);
        _mm512_store_epi64(diff->data + i * VECLEN_I, bvec);
    }

    bmask = borrow;     // result too small, need to add n

    // check for a carry.
    avec = _mm512_load_epi64(sum->data + wshift * VECLEN_I);
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
            cvec = _mm512_load_epi64(sum->data + i * VECLEN_I);
            
            //bvec = _mm512_addcarry_epi52(cvec, carry, &carry);
            __m512i t = _mm512_add_epi64(cvec, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpeq_epu64_mask(cvec, lomask);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(sum->data + i * VECLEN_I, bvec);
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
            cvec = _mm512_load_epi64(sum->data + i * VECLEN_I);
            
            //bvec = _mm512_subborrow_epi52(cvec, carry, &carry);
            __m512i t = _mm512_sub_epi64(cvec, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpeq_epu64_mask(cvec, zero);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(sum->data + i * VECLEN_I, bvec);
        }
    }

    // clear the potential hi-bit
    avec = _mm512_load_epi64(sum->data + wshift * VECLEN_I);
    _mm512_store_epi64(sum->data + wshift * VECLEN_I,
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
            cvec = _mm512_load_epi64(diff->data + i * VECLEN_I);

            //bvec = _mm512_mask_adc_epi52(cvec, bmask, carry, nvec, &carry);
            __m512i t = _mm512_add_epi64(cvec, nvec);
            t = _mm512_mask_add_epi64(cvec, bmask, t, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpgt_epu64_mask(t, lomask);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(diff->data + i * VECLEN_I, bvec);
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
            cvec = _mm512_load_epi64(diff->data + i * VECLEN_I);

            //bvec = _mm512_addcarry_epi52(cvec, carry, &carry);
            __m512i t = _mm512_add_epi64(cvec, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpeq_epu64_mask(cvec, lomask);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(diff->data + i * VECLEN_I, bvec);
        }
    }

    //for (; i < NWORDS; i++)
    //{
    //    _mm512_store_epi64(diff->data + i * VECLEN_I, _mm512_set1_epi64(0));
    //}

    // clear the potential hi-bit
    avec = _mm512_load_epi64(diff->data + wshift * VECLEN_I);
    _mm512_store_epi64(diff->data + wshift * VECLEN_I,
        _mm512_and_epi64(vbshift, avec));


    return;
}

//#define DEBUG_MERSENNE
void vecmod_mersenne(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* s, vec_monty_t* mdata)
{
    int i;
    int NWORDS = (int)mdata->NWORDS;
    // needed in loops
    __m512i a0, a1;                                     // 4
    UNUSED_VAR __m512i b0, b1, b2, b3;                         // 11

#ifndef IFMA
    UNUSED_VAR __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    UNUSED_VAR __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    UNUSED_VAR __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    UNUSED_VAR __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    UNUSED_VAR __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    UNUSED_VAR __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
#endif

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    UNUSED_VAR __m512i acc_e0, acc_e1, acc_e2;
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry;

    // deal with the sign
    c->size = (int)(uint32_t)NWORDS;
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

#ifdef DEBUG_MERSENNE
    mpz_t ga;
    mpz_init(ga);
    extract_bignum_from_vec_to_mpz(ga, a, 0, 2 * NWORDS);
#endif

    if (mdata->use_vnbits)
    {
        int bshift_a[VECLEN];
        for (i = 0; i < VECLEN; i++)
        {
            bshift_a[i] = mdata->vnbits[i] % 52;
        }

        vbshift = _mm512_set_epi64(
            (long long)((1ULL << (uint64_t)(bshift_a[7])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[6])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[5])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[4])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[3])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[2])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[1])) - 1ULL),
            (long long)((1ULL << (uint64_t)(bshift_a[0])) - 1ULL));

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
        vbshift = _mm512_set1_epi64((long long)((1ULL << (uint64_t)(bshift)) - 1ULL));
        vbpshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)));
        vec_bignum52_mask_rshift_n(a, c, mdata->nbits, 0xff);
    }

#ifdef DEBUG_MERSENNE
    //print_vechexbignum(c, "hi part:");
    //print_vechexbignum(a, "lo part:");
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
        a1 = _mm512_load_epi64(s->data + wshift * VECLEN_I);
        _mm512_store_epi64(s->data + wshift * VECLEN_I,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN_I, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(s, "cleared low part:");
#endif

        // add c * hi + lo
        scarry = 0;
        for (i = 0; i <= NWORDS; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
            b0 = _mm512_load_epi64(s->data + i * VECLEN_I);
            a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "lo + hi * c:");
#endif

        // now we need to do it again, but the hi part is just a single word so
        // we only have a single-precision vecmul
        //vec_bignum52_mask_rshift_n(c, s, mdata->nbits, 0xff);
        //a0 = _mm512_load_epi64(s->data);

        // It's possible these hi bits are located in two words
        a0 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        a0 = _mm512_srli_epi64(a0, bshift);
        a1 = _mm512_load_epi64(c->data + (wshift + 1) * VECLEN_I);
        a0 = _mm512_or_epi64(a0, _mm512_slli_epi64(a1, 52 - bshift));

        //print_regvechex(a0, 0, "hi part:");

        b0 = _mm512_set1_epi64(mdata->isMersenne);

        //print_regvechex(b0, 0, "c:");

        VEC_MUL_LOHI_PD(a0, b0, acc_e0, acc_e1);

        // clear any hi-bits now that we have the multiplier
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        _mm512_store_epi64(c->data + wshift * VECLEN_I,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_set1_epi64(0));
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
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
            a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
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
                a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
                a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
                _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
                i++;
            }
        }

        // finally clear any hi-bits in the result
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        _mm512_store_epi64(c->data + wshift * VECLEN_I,
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
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
            b0 = _mm512_load_epi64(a->data + i * VECLEN_I);
            a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
        b0 = _mm512_load_epi64(a->data + i * VECLEN_I);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN_I, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        //print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add it back in.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
            a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // clear the potential hi-bit
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        _mm512_store_epi64(c->data + wshift * VECLEN_I,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        //print_vechexbignum(c, "after carry add:");
#endif

#ifdef DEBUG_MERSENNE
        mpz_t gn, gb, gc;
        mpz_init(gn);
        mpz_init(gb);
        mpz_init(gc);
        mpz_set_ui(gn, 1);
        mpz_mul_2exp(gn, gn, mdata->nbits);
        mpz_sub_ui(gn, gn, 1);
        mpz_mod(gc, ga, gn);
        
        extract_bignum_from_vec_to_mpz(gb, c, 0, NWORDS);
        
        if (mpz_cmp(gc, gb) != 0)
        {
            //gmp_printf("a = %Zd\n", ga);
            //gmp_printf("n = %Zd\n", gn);
            //
            //gmp_printf("b = %Zd\n", gb);
            //gmp_printf("c = %Zd\n", gc);
            exit(1);
        }
        
        
        mpz_clear(gb);
        mpz_clear(gc);
        mpz_clear(gn);
#endif

    }
    else
    {
        // now subtract the hi part from the lo
        scarry = 0;
        for (i = 0; i < wshift; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);   // hi
            b0 = _mm512_load_epi64(a->data + i * VECLEN_I);   // lo
            a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
        b0 = _mm512_load_epi64(a->data + i * VECLEN_I);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN_I, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add 1.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
            a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // and add one to the hi-bit.  This should resolve the borrow bit.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN_I);
        //a1 = _mm512_add_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift))), a1);
        //_mm512_store_epi64(c->data + wshift * VECLEN_I,
        //    _mm512_and_epi64(_mm512_set1_epi64(0xfffffffffffffULL), a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif

        _mm512_store_epi64(c->data + wshift * VECLEN_I,
            _mm512_and_epi64(vbshift, a1));
    }

#ifdef DEBUG_MERSENNE
    mpz_clear(ga);
#endif
    c->size = (int)(uint32_t)NWORDS;

}

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
        avec = _mm512_load_epi64(a + i * VECLEN_I);
        bvec = _mm512_load_epi64(b + i * VECLEN_I);

        tmp = avec;
        tmp = _mm512_mask_mov_epi64(tmp, ~sm, bvec);

        avec = _mm512_mask_mov_epi64(avec, sm, bvec);
        bvec = _mm512_mask_mov_epi64(bvec, sm, tmp);

        cvec = _mm512_sbb_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c + i * VECLEN_I, cvec);
    }

    if (wordsa > wordsb)
    {
        avec = _mm512_load_epi64(a + i * VECLEN_I);
        bvec = _mm512_setzero_si512();

        tmp = avec;
        tmp = _mm512_mask_mov_epi64(tmp, ~sm, bvec);

        avec = _mm512_mask_mov_epi64(avec, sm, bvec);
        bvec = _mm512_mask_mov_epi64(bvec, sm, tmp);

        cvec = _mm512_sbb_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c + i * VECLEN_I, cvec);
    }
    else if (wordsb > wordsa)
    {
        avec = _mm512_setzero_si512();
        bvec = _mm512_load_epi64(b + i * VECLEN_I);

        tmp = avec;
        tmp = _mm512_mask_mov_epi64(tmp, ~sm, bvec);

        avec = _mm512_mask_mov_epi64(avec, sm, bvec);
        bvec = _mm512_mask_mov_epi64(bvec, sm, tmp);

        cvec = _mm512_sbb_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c + i * VECLEN_I, cvec);
    }

    return carry;
}

__mmask8 base_absaddsub_52(uint64_t* a, uint64_t* b, uint64_t* c, 
    __mmask8 addmask, __mmask8 abssubmask, int wordsa, int wordsb)
{
    // if addmask is set then both a and b are negative so we add them 
    // and write any carries.
    // if abssubmask is set then b > a, so we subtract b - a instead of
    // the normal a - b (if neither mask is set).
    // calling code should make sure the two masks are mutually exclusive.
    int i;
    __m512i avec, bvec, cvec, tmp;
    __mmask8 scarry = 0, acarry = 0;

    for (i = 0; i < MIN(wordsa, wordsb); i++)
    {
        avec = _mm512_load_epi64(a + i * VECLEN_I);
        bvec = _mm512_load_epi64(b + i * VECLEN_I);

        // swap words in the lanes we are abs-subbing.
        tmp = avec;
        tmp = _mm512_mask_mov_epi64(tmp, ~abssubmask, bvec);
        avec = _mm512_mask_mov_epi64(avec, abssubmask, bvec);
        bvec = _mm512_mask_mov_epi64(bvec, abssubmask, tmp);

        // sub everything we are not adding.  some of these are swapped
        // abssubs and some are normal subs.
        cvec = _mm512_mask_sbb_epi52(avec, ~addmask, scarry, bvec, &scarry);

        // add the indicated lanes then merge them into the output vector.
        avec = _mm512_mask_adc_epi52(avec, addmask, acarry, bvec, &acarry);
        cvec = _mm512_mask_mov_epi64(cvec, addmask, avec);

        _mm512_store_epi64(c + i * VECLEN_I, cvec);
    }

    if (wordsa > wordsb)
    {
        avec = _mm512_load_epi64(a + i * VECLEN_I);
        bvec = _mm512_setzero_si512();

        // swap words in the lanes we are abs-subbing.
        tmp = avec;
        tmp = _mm512_mask_mov_epi64(tmp, ~abssubmask, bvec);
        avec = _mm512_mask_mov_epi64(avec, abssubmask, bvec);
        bvec = _mm512_mask_mov_epi64(bvec, abssubmask, tmp);

        // sub everything we are not adding.  some of these are swapped
        // abssubs and some are normal subs.
        cvec = _mm512_mask_sbb_epi52(avec, ~addmask, scarry, bvec, &scarry);

        // add the indicated lanes then merge them into the output vector.
        avec = _mm512_mask_adc_epi52(avec, addmask, acarry, bvec, &acarry);
        cvec = _mm512_mask_mov_epi64(cvec, addmask, avec);

        _mm512_store_epi64(c + i * VECLEN_I, cvec);
        i++;
    }
    else if (wordsb > wordsa)
    {
        avec = _mm512_setzero_si512();
        bvec = _mm512_load_epi64(b + i * VECLEN_I);

        // swap words in the lanes we are abs-subbing.
        tmp = avec;
        tmp = _mm512_mask_mov_epi64(tmp, ~abssubmask, bvec);
        avec = _mm512_mask_mov_epi64(avec, abssubmask, bvec);
        bvec = _mm512_mask_mov_epi64(bvec, abssubmask, tmp);

        // sub everything we are not adding.  some of these are swapped
        // abssubs and some are normal subs.
        cvec = _mm512_mask_sbb_epi52(avec, ~addmask, scarry, bvec, &scarry);

        // add the indicated lanes then merge them into the output vector.
        avec = _mm512_mask_adc_epi52(avec, addmask, acarry, bvec, &acarry);
        cvec = _mm512_mask_mov_epi64(cvec, addmask, avec);

        _mm512_store_epi64(c + i * VECLEN_I, cvec);
        i++;
    }

    // write any final addition carries.
    _mm512_store_epi64(c + i * VECLEN_I, _mm512_mask_mov_epi64(
        _mm512_setzero_si512(), acarry, _mm512_set1_epi64(1)));

    return acarry;
}

__mmask8 base_add_52(uint64_t* a, uint64_t* b, uint64_t* c, int wordsa, int wordsb)
{
    int i;
    __m512i avec, bvec, cvec;
    __mmask8 carry = 0;

    for (i = 0; i < MIN(wordsa, wordsb); i++)
    {
        avec = _mm512_load_epi64(a + i * VECLEN_I);
        bvec = _mm512_load_epi64(b + i * VECLEN_I);
        cvec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c + i * VECLEN_I, cvec);
    }

    if (wordsa > wordsb)
    {
        avec = _mm512_load_epi64(a + i * VECLEN_I);
        bvec = _mm512_setzero_si512();
        cvec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c + i * VECLEN_I, cvec);
    }
    else if (wordsb > wordsa)
    {
        avec = _mm512_setzero_si512();
        bvec = _mm512_load_epi64(b + i * VECLEN_I);
        cvec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c + i * VECLEN_I, cvec);
    }

    return carry;
}

__mmask8 base_add_52_wc(uint64_t* a, uint64_t* b, uint64_t* c, int wordsa, int wordsb)
{
    // add and write the carry.
    int i;
    __m512i avec, bvec, cvec;
    __mmask8 carry = 0;

    for (i = 0; i < MIN(wordsa, wordsb); i++)
    {
        avec = _mm512_load_epi64(a + i * VECLEN_I);
        bvec = _mm512_load_epi64(b + i * VECLEN_I);
        cvec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c + i * VECLEN_I, cvec);
    }

    if (wordsa > wordsb)
    {
        avec = _mm512_load_epi64(a + i * VECLEN_I);
        bvec = _mm512_setzero_si512();
        cvec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c + i * VECLEN_I, cvec);
        i++;
    }
    else if (wordsb > wordsa)
    {
        avec = _mm512_setzero_si512();
        bvec = _mm512_load_epi64(b + i * VECLEN_I);
        cvec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c + i * VECLEN_I, cvec);
        i++;
    }

    _mm512_store_epi64(c + i * VECLEN_I, _mm512_mask_mov_epi64(
        _mm512_setzero_si512(), carry, _mm512_set1_epi64(1)));

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
        __m512i a = _mm512_load_epi64(u + i * VECLEN_I);
        __m512i b = _mm512_load_epi64(v + i * VECLEN_I);

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
        v1 = _mm512_load_epi64(c + i * VECLEN_I);
        v2 = _mm512_load_epi64(c + (words + i) * VECLEN_I);
        _mm512_store_epi64(s1 + i * VECLEN_I, _mm512_add_epi64(v1, v2));
    }

    for (; i < words; i++)
    {
        v1 = _mm512_load_epi64(c + i * VECLEN_I);
        _mm512_store_epi64(s1 + i * VECLEN_I, v1);
    }

    i = halfwords;
    v1 = _mm512_load_epi64(c + i * VECLEN_I);
    v2 = _mm512_load_epi64(s1 + (i - halfwords) * VECLEN_I);
    v3 = _mm512_load_epi64(z1 + (i - halfwords) * VECLEN_I);

    v1 = _mm512_add_epi64(v1, v2);
    v3 = _mm512_mask_xor_epi64(v3, sm, v3, vlmask);
    v1 = _mm512_mask_add_epi64(v1, sm, v1, vone);

    v3 = _mm512_add_epi64(v3, vc);
    v1 = _mm512_add_epi64(v1, v3);

    // carryprop
    vc = _mm512_srli_epi64(v1, 52);
    v1 = _mm512_and_epi64(vlmask, v1);

    _mm512_store_epi64(c + i * VECLEN_I, v1);

    for (i = halfwords + 1; i < halfwords + words; i++)
    {
        v1 = _mm512_load_epi64(c + i * VECLEN_I);
        v2 = _mm512_load_epi64(s1 + (i - halfwords) * VECLEN_I);
        v3 = _mm512_load_epi64(z1 + (i - halfwords) * VECLEN_I);

        v1 = _mm512_add_epi64(v1, v2);
        v3 = _mm512_mask_xor_epi64(v3, sm, v3, vlmask);

        v3 = _mm512_add_epi64(v3, vc);
        v1 = _mm512_add_epi64(v1, v3);

        // carryprop
        vc = _mm512_srli_epi64(v1, 52);
        v1 = _mm512_and_epi64(vlmask, v1);

        _mm512_store_epi64(c + i * VECLEN_I, v1);
    }

    // add in the final hi-word carry then propagate any
    // additional carries through the final words of the output.
    // z1 was negative, so we have to subtract a sign bit from
    // this word.
    v1 = _mm512_load_epi64(c + i * VECLEN_I);
    v1 = _mm512_addsetc_epi52(v1, vc, &carry);
    v1 = _mm512_mask_sub_epi64(v1, sm, v1, vone);
    _mm512_store_epi64(c + i * VECLEN_I, v1);
    i++;

    while (carry)
    {
        v1 = _mm512_load_epi64(c + i * VECLEN_I);
        v1 = _mm512_addcarry_epi52(v1, carry, &carry);
        _mm512_store_epi64(c + i * VECLEN_I, v1);
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

    for (i = 0; i < halfwords; i++)
    {
        v1 = _mm512_load_epi64(c + i * VECLEN_I);
        v2 = _mm512_load_epi64(c + (words + i) * VECLEN_I);
        _mm512_store_epi64(s1 + i * VECLEN_I, _mm512_add_epi64(v1, v2));
    }

    i = halfwords;
    v1 = _mm512_load_epi64(c + i * VECLEN_I);
    v2 = _mm512_load_epi64(s1 + (i - halfwords) * VECLEN_I);
    v3 = _mm512_load_epi64(z1 + (i - halfwords) * VECLEN_I);

    v1 = _mm512_add_epi64(v1, v2);
    v3 = _mm512_mask_xor_epi64(v3, sm, v3, vlmask);
    v1 = _mm512_mask_add_epi64(v1, sm, v1, vone);

    v3 = _mm512_add_epi64(v3, vc);
    v1 = _mm512_add_epi64(v1, v3);

    // carryprop
    vc = _mm512_srli_epi64(v1, 52);
    v1 = _mm512_and_epi64(vlmask, v1);

    _mm512_store_epi64(c + i * VECLEN_I, v1);

    for (i = halfwords + 1; i < words; i++)
    {
        v1 = _mm512_load_epi64(c + i * VECLEN_I);
        v2 = _mm512_load_epi64(s1 + (i - halfwords) * VECLEN_I);
        v3 = _mm512_load_epi64(z1 + (i - halfwords) * VECLEN_I);

        v1 = _mm512_add_epi64(v1, v2);
        v3 = _mm512_mask_xor_epi64(v3, sm, v3, vlmask);

        v3 = _mm512_add_epi64(v3, vc);
        v1 = _mm512_add_epi64(v1, v3);

        // carryprop
        vc = _mm512_srli_epi64(v1, 52);
        v1 = _mm512_and_epi64(vlmask, v1);

        _mm512_store_epi64(c + i * VECLEN_I, v1);
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
    int i, j;

    __m512i lo0, lo1, lo2, lo3, lo4, lo5, lo6;
    __m512i hi0, hi1, hi2, hi3, hi4, hi5, hi6;
    UNUSED_VAR __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;
    UNUSED_VAR __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);
    UNUSED_VAR __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);
    UNUSED_VAR __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);
    __m512i zero = set64(0), B, A, C, plo, carry;
    __m512i MASK = set64(DIGIT_MASK);

    int full_blocks = n / BLOCKWORDS_I;

#ifdef DEBUG_VECMUL
    printf("commencing vecmul52_cios with WORDS=%d and NBLOCKS=%d\n", n, full_blocks);
    print_vechex(a, 0, n, "input a: ");
    print_vechex(b, 0, n, "input b: ");
    print_vechex(c, 0, 2 * n, "output c: ");
#endif

    for (i = 0; i < 2 * n; i++) {
        _mm512_store_epi64(c + i * VECLEN_I, _mm512_setzero_si512());
    }
    
    for (j = 0; j < full_blocks; j++) {
        __m512i B0 = loadu64(b + (j * BLOCKWORDS_I + 0) * VECLEN_I);
        __m512i B1 = loadu64(b + (j * BLOCKWORDS_I + 1) * VECLEN_I);
        __m512i B2 = loadu64(b + (j * BLOCKWORDS_I + 2) * VECLEN_I);
        __m512i B3 = loadu64(b + (j * BLOCKWORDS_I + 3) * VECLEN_I);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            __m512i A0 = loadu64(a + (i * BLOCKWORDS_I + 0) * VECLEN_I);
            __m512i A1 = loadu64(a + (i * BLOCKWORDS_I + 1) * VECLEN_I);
            __m512i A2 = loadu64(a + (i * BLOCKWORDS_I + 2) * VECLEN_I);
            __m512i A3 = loadu64(a + (i * BLOCKWORDS_I + 3) * VECLEN_I);

            VEC_MUL4_ACCUMX(A0, B0, B1, B2, B3, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, B1, B2, B3, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, B1, B2, B3, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, B1, B2, B3, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

            A0 = loadu64(c + ((i + j) * BLOCKWORDS_I + 0) * VECLEN_I);
            A1 = loadu64(c + ((i + j) * BLOCKWORDS_I + 1) * VECLEN_I);
            A2 = loadu64(c + ((i + j) * BLOCKWORDS_I + 2) * VECLEN_I);
            A3 = loadu64(c + ((i + j) * BLOCKWORDS_I + 3) * VECLEN_I);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS_I + 0) * VECLEN_I, 
                _mm512_and_si512(lo0, MASK));
            
            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS_I + 1) * VECLEN_I,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS_I + 2) * VECLEN_I,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS_I + 3) * VECLEN_I,
                _mm512_and_si512(lo3, MASK));

            hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

            // subtract out all of the bias at once.
            SUB_BIAS_HI3(3, 2, 1);
            SUB_BIAS_LO3(3, 2, 1);

            A0 = loadu64(c + ((i + j) * BLOCKWORDS_I + 4) * VECLEN_I);
            A1 = loadu64(c + ((i + j) * BLOCKWORDS_I + 5) * VECLEN_I);
            A2 = loadu64(c + ((i + j) * BLOCKWORDS_I + 6) * VECLEN_I);
            A3 = loadu64(c + ((i + j) * BLOCKWORDS_I + 7) * VECLEN_I);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS_I + 4) * VECLEN_I,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS_I + 5) * VECLEN_I,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS_I + 6) * VECLEN_I,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS_I + 7) * VECLEN_I,
                _mm512_and_si512(hi6, MASK));

            if (((i + j) * BLOCKWORDS_I + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i + j) * BLOCKWORDS_I + 8) * VECLEN_I, 
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i + j) * BLOCKWORDS_I + 8) * VECLEN_I), 
                        _mm512_srli_epi64(hi6, 52)));
            }

#ifdef DEBUG_VECMUL
            print_vechex(c, 0, 2 * n, "end of block iteration:");
#endif
        }
    }

    if (n - (full_blocks * BLOCKWORDS_I) == 3)
    {
        j = full_blocks * BLOCKWORDS_I;
        __m512i B0 = loadu64(b + (j + 0) * VECLEN_I);
        __m512i B1 = loadu64(b + (j + 1) * VECLEN_I);
        __m512i B2 = loadu64(b + (j + 2) * VECLEN_I);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            __m512i A0 = loadu64(a + (i * BLOCKWORDS_I + 0) * VECLEN_I);
            __m512i A1 = loadu64(a + (i * BLOCKWORDS_I + 1) * VECLEN_I);
            __m512i A2 = loadu64(a + (i * BLOCKWORDS_I + 2) * VECLEN_I);
            __m512i A3 = loadu64(a + (i * BLOCKWORDS_I + 3) * VECLEN_I);

            VEC_MUL4_ACCUMX(A0, B0, B1, B2, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, B1, B2, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, B1, B2, zero, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, B1, B2, zero, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

            A0 = loadu64(c + ((i) * BLOCKWORDS_I + j + 0) * VECLEN_I);
            A1 = loadu64(c + ((i) * BLOCKWORDS_I + j + 1) * VECLEN_I);
            A2 = loadu64(c + ((i) * BLOCKWORDS_I + j + 2) * VECLEN_I);
            A3 = loadu64(c + ((i) * BLOCKWORDS_I + j + 3) * VECLEN_I);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 0) * VECLEN_I,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 1) * VECLEN_I,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 2) * VECLEN_I,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 3) * VECLEN_I,
                _mm512_and_si512(lo3, MASK));

            hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

            // subtract out all of the bias at once.
            SUB_BIAS_HI3(3, 2, 1);
            SUB_BIAS_LO3(3, 2, 1);

            A0 = loadu64(c + ((i) * BLOCKWORDS_I + j + 4) * VECLEN_I);
            A1 = loadu64(c + ((i) * BLOCKWORDS_I + j + 5) * VECLEN_I);
            A2 = loadu64(c + ((i) * BLOCKWORDS_I + j + 6) * VECLEN_I);
            A3 = loadu64(c + ((i) * BLOCKWORDS_I + j + 7) * VECLEN_I);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 4) * VECLEN_I,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 5) * VECLEN_I,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 6) * VECLEN_I,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 7) * VECLEN_I,
                _mm512_and_si512(hi6, MASK));

            if (((i)*BLOCKWORDS_I + j + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 8) * VECLEN_I,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i)*BLOCKWORDS_I + j + 8) * VECLEN_I),
                        _mm512_srli_epi64(hi6, 52)));
            }
        }

        j = full_blocks * BLOCKWORDS_I;
        B0 = loadu64(a + (j + 0) * VECLEN_I);
        B1 = loadu64(a + (j + 1) * VECLEN_I);
        B2 = loadu64(a + (j + 2) * VECLEN_I);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            __m512i A0 = loadu64(b + (i * BLOCKWORDS_I + 0) * VECLEN_I);
            __m512i A1 = loadu64(b + (i * BLOCKWORDS_I + 1) * VECLEN_I);
            __m512i A2 = loadu64(b + (i * BLOCKWORDS_I + 2) * VECLEN_I);
            __m512i A3 = loadu64(b + (i * BLOCKWORDS_I + 3) * VECLEN_I);

            VEC_MUL4_ACCUMX(A0, B0, B1, B2, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, B1, B2, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, B1, B2, zero, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, B1, B2, zero, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

            A0 = loadu64(c + ((i)*BLOCKWORDS_I + j + 0) * VECLEN_I);
            A1 = loadu64(c + ((i)*BLOCKWORDS_I + j + 1) * VECLEN_I);
            A2 = loadu64(c + ((i)*BLOCKWORDS_I + j + 2) * VECLEN_I);
            A3 = loadu64(c + ((i)*BLOCKWORDS_I + j + 3) * VECLEN_I);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 0) * VECLEN_I,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 1) * VECLEN_I,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 2) * VECLEN_I,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 3) * VECLEN_I,
                _mm512_and_si512(lo3, MASK));

            hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

            // subtract out all of the bias at once.
            SUB_BIAS_HI3(3, 2, 1);
            SUB_BIAS_LO3(3, 2, 1);

            A0 = loadu64(c + ((i)*BLOCKWORDS_I + j + 4) * VECLEN_I);
            A1 = loadu64(c + ((i)*BLOCKWORDS_I + j + 5) * VECLEN_I);
            A2 = loadu64(c + ((i)*BLOCKWORDS_I + j + 6) * VECLEN_I);
            A3 = loadu64(c + ((i)*BLOCKWORDS_I + j + 7) * VECLEN_I);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 4) * VECLEN_I,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 5) * VECLEN_I,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 6) * VECLEN_I,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 7) * VECLEN_I,
                _mm512_and_si512(hi6, MASK));

            if (((i)*BLOCKWORDS_I + j + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 8) * VECLEN_I,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i)*BLOCKWORDS_I + j + 8) * VECLEN_I),
                        _mm512_srli_epi64(hi6, 52)));
            }
        }

#ifdef DEBUG_VECMUL
        print_vechex(c, 0, 2 * n, "end of last partial block:");
#endif

        j = full_blocks * BLOCKWORDS_I;
        i = full_blocks * BLOCKWORDS_I;
        lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
        hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

        B0 = loadu64(b + (j + 0) * VECLEN_I);
        B1 = loadu64(b + (j + 1) * VECLEN_I);
        B2 = loadu64(b + (j + 2) * VECLEN_I);
        __m512i A0 = loadu64(a + (i + 0) * VECLEN_I);
        __m512i A1 = loadu64(a + (i + 1) * VECLEN_I);
        __m512i A2 = loadu64(a + (i + 2) * VECLEN_I);
        __m512i A3 = zero;

        VEC_MUL4_ACCUMX(A0, B0, B1, B2, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
        VEC_MUL4_ACCUMX(A1, B0, B1, B2, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
        VEC_MUL4_ACCUMX(A2, B0, B1, B2, zero, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);

        // subtract out all of the bias at once.
        SUB_BIAS_HI4(1, 2, 3, 3);
        SUB_BIAS_LO4(1, 2, 3, 3);

        A0 = loadu64(c + ((i + j) + 0) * VECLEN_I);
        A1 = loadu64(c + ((i + j) + 1) * VECLEN_I);
        A2 = loadu64(c + ((i + j) + 2) * VECLEN_I);
        A3 = loadu64(c + ((i + j) + 3) * VECLEN_I);

        lo0 = _mm512_add_epi64(lo0, A0);
        lo1 = _mm512_add_epi64(lo1, A1);
        lo2 = _mm512_add_epi64(lo2, A2);
        lo3 = _mm512_add_epi64(lo3, A3);

        lo1 = _mm512_add_epi64(lo1, hi0);
        lo2 = _mm512_add_epi64(lo2, hi1);
        lo3 = _mm512_add_epi64(lo3, hi2);

        _mm512_store_epi64(c + ((i + j) + 0) * VECLEN_I,
            _mm512_and_si512(lo0, MASK));

        hi0 = _mm512_srli_epi64(lo0, 52);
        lo1 = _mm512_add_epi64(lo1, hi0);

        _mm512_store_epi64(c + ((i + j) + 1) * VECLEN_I,
            _mm512_and_si512(lo1, MASK));

        hi1 = _mm512_srli_epi64(lo1, 52);
        lo2 = _mm512_add_epi64(lo2, hi1);

        _mm512_store_epi64(c + ((i + j) + 2) * VECLEN_I,
            _mm512_and_si512(lo2, MASK));

        hi2 = _mm512_srli_epi64(lo2, 52);
        lo3 = _mm512_add_epi64(lo3, hi2);

        _mm512_store_epi64(c + ((i + j) + 3) * VECLEN_I,
            _mm512_and_si512(lo3, MASK));

        hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

        // subtract out all of the bias at once.
        SUB_BIAS_HI3(2, 1, 0);
        SUB_BIAS_LO3(2, 1, 0);

        A0 = loadu64(c + ((i + j) + 4) * VECLEN_I);
        A1 = loadu64(c + ((i + j) + 5) * VECLEN_I);

        lo4 = _mm512_add_epi64(lo4, A0);
        lo5 = _mm512_add_epi64(lo5, A1);

        lo4 = _mm512_add_epi64(lo4, hi3);
        lo5 = _mm512_add_epi64(lo5, hi4);

        _mm512_store_epi64(c + ((i + j) + 4) * VECLEN_I,
            _mm512_and_si512(lo4, MASK));

        hi4 = _mm512_srli_epi64(lo4, 52);
        lo5 = _mm512_add_epi64(lo5, hi4);

        _mm512_store_epi64(c + ((i + j) + 5) * VECLEN_I,
            _mm512_and_si512(lo5, MASK));
    }

    if (n - (full_blocks * BLOCKWORDS_I) == 2)
    {
        j = full_blocks * BLOCKWORDS_I;
        __m512i B0 = loadu64(b + (j + 0) * VECLEN_I);
        __m512i B1 = loadu64(b + (j + 1) * VECLEN_I);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            __m512i A0 = loadu64(a + (i * BLOCKWORDS_I + 0) * VECLEN_I);
            __m512i A1 = loadu64(a + (i * BLOCKWORDS_I + 1) * VECLEN_I);
            __m512i A2 = loadu64(a + (i * BLOCKWORDS_I + 2) * VECLEN_I);
            __m512i A3 = loadu64(a + (i * BLOCKWORDS_I + 3) * VECLEN_I);

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

            A0 = loadu64(c + ((i) * BLOCKWORDS_I + j + 0) * VECLEN_I);
            A1 = loadu64(c + ((i) * BLOCKWORDS_I + j + 1) * VECLEN_I);
            A2 = loadu64(c + ((i) * BLOCKWORDS_I + j + 2) * VECLEN_I);
            A3 = loadu64(c + ((i) * BLOCKWORDS_I + j + 3) * VECLEN_I);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 0) * VECLEN_I,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 1) * VECLEN_I,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 2) * VECLEN_I,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 3) * VECLEN_I,
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

            A0 = loadu64(c + ((i) * BLOCKWORDS_I + j + 4) * VECLEN_I);
            A1 = loadu64(c + ((i) * BLOCKWORDS_I + j + 5) * VECLEN_I);
            A2 = loadu64(c + ((i) * BLOCKWORDS_I + j + 6) * VECLEN_I);
            A3 = loadu64(c + ((i) * BLOCKWORDS_I + j + 7) * VECLEN_I);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 4) * VECLEN_I,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 5) * VECLEN_I,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 6) * VECLEN_I,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 7) * VECLEN_I,
                _mm512_and_si512(hi6, MASK));

            if (((i)*BLOCKWORDS_I + j + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 8) * VECLEN_I,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i)*BLOCKWORDS_I + j + 8) * VECLEN_I),
                        _mm512_srli_epi64(hi6, 52)));
            }
        }

        B0 = loadu64(a + (j + 0) * VECLEN_I);
        B1 = loadu64(a + (j + 1) * VECLEN_I);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            __m512i A0 = loadu64(b + (i * BLOCKWORDS_I + 0) * VECLEN_I);
            __m512i A1 = loadu64(b + (i * BLOCKWORDS_I + 1) * VECLEN_I);
            __m512i A2 = loadu64(b + (i * BLOCKWORDS_I + 2) * VECLEN_I);
            __m512i A3 = loadu64(b + (i * BLOCKWORDS_I + 3) * VECLEN_I);

            VEC_MUL4_ACCUMX(A0, B0, B1, zero, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, B1, zero, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, B1, zero, zero, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, B1, zero, zero, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

            A0 = loadu64(c + ((i)*BLOCKWORDS_I + j + 0) * VECLEN_I);
            A1 = loadu64(c + ((i)*BLOCKWORDS_I + j + 1) * VECLEN_I);
            A2 = loadu64(c + ((i)*BLOCKWORDS_I + j + 2) * VECLEN_I);
            A3 = loadu64(c + ((i)*BLOCKWORDS_I + j + 3) * VECLEN_I);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 0) * VECLEN_I,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 1) * VECLEN_I,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 2) * VECLEN_I,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 3) * VECLEN_I,
                _mm512_and_si512(lo3, MASK));

            hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

            // subtract out all of the bias at once.
            SUB_BIAS_HI3(3, 2, 1);
            SUB_BIAS_LO3(3, 2, 1);

            A0 = loadu64(c + ((i)*BLOCKWORDS_I + j + 4) * VECLEN_I);
            A1 = loadu64(c + ((i)*BLOCKWORDS_I + j + 5) * VECLEN_I);
            A2 = loadu64(c + ((i)*BLOCKWORDS_I + j + 6) * VECLEN_I);
            A3 = loadu64(c + ((i)*BLOCKWORDS_I + j + 7) * VECLEN_I);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 4) * VECLEN_I,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 5) * VECLEN_I,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 6) * VECLEN_I,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 7) * VECLEN_I,
                _mm512_and_si512(hi6, MASK));

            if (((i)*BLOCKWORDS_I + j + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 8) * VECLEN_I,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i)*BLOCKWORDS_I + j + 8) * VECLEN_I),
                        _mm512_srli_epi64(hi6, 52)));
            }
        }

        j = full_blocks * BLOCKWORDS_I;
        i = full_blocks * BLOCKWORDS_I;
        lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
        hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

        B0 = loadu64(b + (j + 0) * VECLEN_I);
        B1 = loadu64(b + (j + 1) * VECLEN_I);
        __m512i A0 = loadu64(a + (i + 0) * VECLEN_I);
        __m512i A1 = loadu64(a + (i + 1) * VECLEN_I);
        __m512i A2 = zero;
        __m512i A3 = zero;

        VEC_MUL4_ACCUMX(A0, B0, B1, zero, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
        VEC_MUL4_ACCUMX(A1, B0, B1, zero, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);

        // subtract out all of the bias at once.
        SUB_BIAS_HI4(1, 2, 2, 2);
        SUB_BIAS_LO4(1, 2, 2, 2);

        A0 = loadu64(c + ((i + j) + 0) * VECLEN_I);
        A1 = loadu64(c + ((i + j) + 1) * VECLEN_I);
        A2 = loadu64(c + ((i + j) + 2) * VECLEN_I);
        A3 = loadu64(c + ((i + j) + 3) * VECLEN_I);

        lo0 = _mm512_add_epi64(lo0, A0);
        lo1 = _mm512_add_epi64(lo1, A1);
        lo2 = _mm512_add_epi64(lo2, A2);
        lo3 = _mm512_add_epi64(lo3, A3);

        lo1 = _mm512_add_epi64(lo1, hi0);
        lo2 = _mm512_add_epi64(lo2, hi1);
        lo3 = _mm512_add_epi64(lo3, hi2);

        _mm512_store_epi64(c + ((i + j) + 0) * VECLEN_I,
            _mm512_and_si512(lo0, MASK));

        hi0 = _mm512_srli_epi64(lo0, 52);
        lo1 = _mm512_add_epi64(lo1, hi0);

        _mm512_store_epi64(c + ((i + j) + 1) * VECLEN_I,
            _mm512_and_si512(lo1, MASK));

        hi1 = _mm512_srli_epi64(lo1, 52);
        lo2 = _mm512_add_epi64(lo2, hi1);

        _mm512_store_epi64(c + ((i + j) + 2) * VECLEN_I,
            _mm512_and_si512(lo2, MASK));

        hi2 = _mm512_srli_epi64(lo2, 52);
        lo3 = _mm512_add_epi64(lo3, hi2);

        _mm512_store_epi64(c + ((i + j) + 3) * VECLEN_I,
            _mm512_and_si512(lo3, MASK));
    }

    if (n - (full_blocks * BLOCKWORDS_I) == 1)
    {
        j = full_blocks * BLOCKWORDS_I;
        __m512i B0 = loadu64(b + (j + 0) * VECLEN_I);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            __m512i A0 = loadu64(a + (i * BLOCKWORDS_I + 0) * VECLEN_I);
            __m512i A1 = loadu64(a + (i * BLOCKWORDS_I + 1) * VECLEN_I);
            __m512i A2 = loadu64(a + (i * BLOCKWORDS_I + 2) * VECLEN_I);
            __m512i A3 = loadu64(a + (i * BLOCKWORDS_I + 3) * VECLEN_I);

            VEC_MUL4_ACCUMX(A0, B0, zero, zero, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, zero, zero, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, zero, zero, zero, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, zero, zero, zero, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

            A0 = loadu64(c + ((i)*BLOCKWORDS_I + j + 0) * VECLEN_I);
            A1 = loadu64(c + ((i)*BLOCKWORDS_I + j + 1) * VECLEN_I);
            A2 = loadu64(c + ((i)*BLOCKWORDS_I + j + 2) * VECLEN_I);
            A3 = loadu64(c + ((i)*BLOCKWORDS_I + j + 3) * VECLEN_I);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 0) * VECLEN_I,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 1) * VECLEN_I,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 2) * VECLEN_I,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 3) * VECLEN_I,
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

            A0 = loadu64(c + ((i)*BLOCKWORDS_I + j + 4) * VECLEN_I);
            A1 = loadu64(c + ((i)*BLOCKWORDS_I + j + 5) * VECLEN_I);
            A2 = loadu64(c + ((i)*BLOCKWORDS_I + j + 6) * VECLEN_I);
            A3 = loadu64(c + ((i)*BLOCKWORDS_I + j + 7) * VECLEN_I);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 4) * VECLEN_I,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 5) * VECLEN_I,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 6) * VECLEN_I,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 7) * VECLEN_I,
                _mm512_and_si512(hi6, MASK));

            if (((i)*BLOCKWORDS_I + j + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 8) * VECLEN_I,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i)*BLOCKWORDS_I + j + 8) * VECLEN_I),
                        _mm512_srli_epi64(hi6, 52)));
            }
        }

        B0 = loadu64(a + (j + 0) * VECLEN_I);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            __m512i A0 = loadu64(b + (i * BLOCKWORDS_I + 0) * VECLEN_I);
            __m512i A1 = loadu64(b + (i * BLOCKWORDS_I + 1) * VECLEN_I);
            __m512i A2 = loadu64(b + (i * BLOCKWORDS_I + 2) * VECLEN_I);
            __m512i A3 = loadu64(b + (i * BLOCKWORDS_I + 3) * VECLEN_I);

            VEC_MUL4_ACCUMX(A0, B0, zero, zero, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, zero, zero, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, zero, zero, zero, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, zero, zero, zero, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

            A0 = loadu64(c + ((i)*BLOCKWORDS_I + j + 0) * VECLEN_I);
            A1 = loadu64(c + ((i)*BLOCKWORDS_I + j + 1) * VECLEN_I);
            A2 = loadu64(c + ((i)*BLOCKWORDS_I + j + 2) * VECLEN_I);
            A3 = loadu64(c + ((i)*BLOCKWORDS_I + j + 3) * VECLEN_I);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 0) * VECLEN_I,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 1) * VECLEN_I,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 2) * VECLEN_I,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 3) * VECLEN_I,
                _mm512_and_si512(lo3, MASK));

            hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

            // subtract out all of the bias at once.
            SUB_BIAS_HI3(3, 2, 1);
            SUB_BIAS_LO3(3, 2, 1);

            A0 = loadu64(c + ((i)*BLOCKWORDS_I + j + 4) * VECLEN_I);
            A1 = loadu64(c + ((i)*BLOCKWORDS_I + j + 5) * VECLEN_I);
            A2 = loadu64(c + ((i)*BLOCKWORDS_I + j + 6) * VECLEN_I);
            A3 = loadu64(c + ((i)*BLOCKWORDS_I + j + 7) * VECLEN_I);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 4) * VECLEN_I,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 5) * VECLEN_I,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 6) * VECLEN_I,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 7) * VECLEN_I,
                _mm512_and_si512(hi6, MASK));

            if (((i)*BLOCKWORDS_I + j + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 8) * VECLEN_I,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i)*BLOCKWORDS_I + j + 8) * VECLEN_I),
                        _mm512_srli_epi64(hi6, 52)));
            }
        }

        j = full_blocks * BLOCKWORDS_I;
        i = full_blocks * BLOCKWORDS_I;
        lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
        hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

        B0 = loadu64(b + (j + 0) * VECLEN_I);
        __m512i A0 = loadu64(a + (i + 0) * VECLEN_I);
        __m512i A1 = zero;
        UNUSED_VAR __m512i A2 = zero;
        UNUSED_VAR __m512i A3 = zero;

        VEC_MUL4_ACCUMX(A0, B0, zero, zero, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);

        // subtract out all of the bias at once.
        SUB_BIAS_HI4(1, 1, 1, 1);
        SUB_BIAS_LO4(1, 1, 1, 1);

        A0 = loadu64(c + ((i + j) + 0) * VECLEN_I);
        A1 = loadu64(c + ((i + j) + 1) * VECLEN_I);

        lo0 = _mm512_add_epi64(lo0, A0);
        lo1 = _mm512_add_epi64(lo1, A1);

        lo1 = _mm512_add_epi64(lo1, hi0);

        _mm512_store_epi64(c + ((i + j) + 0) * VECLEN_I,
            _mm512_and_si512(lo0, MASK));

        hi0 = _mm512_srli_epi64(lo0, 52);
        lo1 = _mm512_add_epi64(lo1, hi0);

        _mm512_store_epi64(c + ((i + j) + 1) * VECLEN_I,
            _mm512_and_si512(lo1, MASK));
    }

    if (0) {
        for (j = full_blocks * BLOCKWORDS_I; j < n; j++) {
            B = loadu64(b + j * VECLEN_I);

            carry = zero;
            for (i = 0; i < n; i++) {
                A = loadu64(a + i * VECLEN_I);
                C = loadu64(c + (i + j) * VECLEN_I);

                plo = fma52lo(C, A, B);
                plo = _mm512_add_epi64(plo, carry);
                carry = fma52hi(_mm512_srli_epi64(plo, 52), A, B, dbias, vbias1);

                _mm512_store_epi64(c + (i + j) * VECLEN_I, _mm512_and_si512(plo, MASK));
            }
            _mm512_store_epi64(c + (i + j) * VECLEN_I, _mm512_add_epi64(
                loadu64(c + (i + j) * VECLEN_I), _mm512_and_si512(carry, MASK)));
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
    int i, j;

    __m512i lo0, lo1, lo2, lo3, lo4, lo5, lo6;
    __m512i hi0, hi1, hi2, hi3, hi4, hi5, hi6;
    UNUSED_VAR __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;
    UNUSED_VAR __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;
    UNUSED_VAR __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    UNUSED_VAR __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);
    UNUSED_VAR __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);
    UNUSED_VAR __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);
    __m512i zero = set64(0);
    __m512i MASK = set64(DIGIT_MASK);

    int full_blocks = n / BLOCKWORDS_I;

#ifdef DEBUG_VECMUL
    printf("commencing vecmul52_cios with WORDS=%d and NBLOCKS=%d\n", n, full_blocks);
    print_vechex(a, 0, n, "input a: ");
    print_vechex(b, 0, n, "input b: ");
    print_vechex(c, 0, 2 * n, "output c: ");
#endif

    for (i = 0; i < 2 * n; i++) {
        _mm512_store_epi64(c + i * VECLEN_I, _mm512_setzero_si512());
    }

    for (j = 0; j < full_blocks; j++) {
        __m512i B0 = loadu64(b + (j * BLOCKWORDS_I + 0) * VECLEN_I);
        __m512i B1 = loadu64(b + (j * BLOCKWORDS_I + 1) * VECLEN_I);
        __m512i B2 = loadu64(b + (j * BLOCKWORDS_I + 2) * VECLEN_I);
        __m512i B3 = loadu64(b + (j * BLOCKWORDS_I + 3) * VECLEN_I);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            if ((i + j) * BLOCKWORDS_I > n)
                continue;

            __m512i A0 = loadu64(a + (i * BLOCKWORDS_I + 0) * VECLEN_I);
            __m512i A1 = loadu64(a + (i * BLOCKWORDS_I + 1) * VECLEN_I);
            __m512i A2 = loadu64(a + (i * BLOCKWORDS_I + 2) * VECLEN_I);
            __m512i A3 = loadu64(a + (i * BLOCKWORDS_I + 3) * VECLEN_I);

            VEC_MUL4_ACCUMX(A0, B0, B1, B2, B3, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, B1, B2, B3, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, B1, B2, B3, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, B1, B2, B3, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

            A0 = loadu64(c + ((i + j) * BLOCKWORDS_I + 0) * VECLEN_I);
            A1 = loadu64(c + ((i + j) * BLOCKWORDS_I + 1) * VECLEN_I);
            A2 = loadu64(c + ((i + j) * BLOCKWORDS_I + 2) * VECLEN_I);
            A3 = loadu64(c + ((i + j) * BLOCKWORDS_I + 3) * VECLEN_I);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS_I + 0) * VECLEN_I,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS_I + 1) * VECLEN_I,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS_I + 2) * VECLEN_I,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS_I + 3) * VECLEN_I,
                _mm512_and_si512(lo3, MASK));

            hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

            // subtract out all of the bias at once.
            SUB_BIAS_HI3(3, 2, 1);
            SUB_BIAS_LO3(3, 2, 1);

            A0 = loadu64(c + ((i + j) * BLOCKWORDS_I + 4) * VECLEN_I);
            A1 = loadu64(c + ((i + j) * BLOCKWORDS_I + 5) * VECLEN_I);
            A2 = loadu64(c + ((i + j) * BLOCKWORDS_I + 6) * VECLEN_I);
            A3 = loadu64(c + ((i + j) * BLOCKWORDS_I + 7) * VECLEN_I);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS_I + 4) * VECLEN_I,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS_I + 5) * VECLEN_I,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS_I + 6) * VECLEN_I,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i + j) * BLOCKWORDS_I + 7) * VECLEN_I,
                _mm512_and_si512(hi6, MASK));

            if (((i + j) * BLOCKWORDS_I + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i + j) * BLOCKWORDS_I + 8) * VECLEN_I,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i + j) * BLOCKWORDS_I + 8) * VECLEN_I),
                        _mm512_srli_epi64(hi6, 52)));
            }

#ifdef DEBUG_VECMUL
            print_vechex(c, 0, 2 * n, "end of block iteration:");
#endif
        }
    }

    if (n - (full_blocks * BLOCKWORDS_I) == 3)
    {
        j = full_blocks * BLOCKWORDS_I;
        __m512i B0 = loadu64(b + (j + 0) * VECLEN_I);
        __m512i B1 = loadu64(b + (j + 1) * VECLEN_I);
        __m512i B2 = loadu64(b + (j + 2) * VECLEN_I);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            if ((i * BLOCKWORDS_I + j) > n)
                continue;

            __m512i A0 = loadu64(a + (i * BLOCKWORDS_I + 0) * VECLEN_I);
            __m512i A1 = loadu64(a + (i * BLOCKWORDS_I + 1) * VECLEN_I);
            __m512i A2 = loadu64(a + (i * BLOCKWORDS_I + 2) * VECLEN_I);
            __m512i A3 = loadu64(a + (i * BLOCKWORDS_I + 3) * VECLEN_I);

            VEC_MUL4_ACCUMX(A0, B0, B1, B2, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, B1, B2, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, B1, B2, zero, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, B1, B2, zero, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

            A0 = loadu64(c + ((i)*BLOCKWORDS_I + j + 0) * VECLEN_I);
            A1 = loadu64(c + ((i)*BLOCKWORDS_I + j + 1) * VECLEN_I);
            A2 = loadu64(c + ((i)*BLOCKWORDS_I + j + 2) * VECLEN_I);
            A3 = loadu64(c + ((i)*BLOCKWORDS_I + j + 3) * VECLEN_I);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 0) * VECLEN_I,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 1) * VECLEN_I,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 2) * VECLEN_I,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 3) * VECLEN_I,
                _mm512_and_si512(lo3, MASK));

            hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

            // subtract out all of the bias at once.
            SUB_BIAS_HI3(3, 2, 1);
            SUB_BIAS_LO3(3, 2, 1);

            A0 = loadu64(c + ((i)*BLOCKWORDS_I + j + 4) * VECLEN_I);
            A1 = loadu64(c + ((i)*BLOCKWORDS_I + j + 5) * VECLEN_I);
            A2 = loadu64(c + ((i)*BLOCKWORDS_I + j + 6) * VECLEN_I);
            A3 = loadu64(c + ((i)*BLOCKWORDS_I + j + 7) * VECLEN_I);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 4) * VECLEN_I,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 5) * VECLEN_I,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 6) * VECLEN_I,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 7) * VECLEN_I,
                _mm512_and_si512(hi6, MASK));

            if (((i)*BLOCKWORDS_I + j + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 8) * VECLEN_I,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i)*BLOCKWORDS_I + j + 8) * VECLEN_I),
                        _mm512_srli_epi64(hi6, 52)));
            }
        }

        j = full_blocks * BLOCKWORDS_I;
        B0 = loadu64(a + (j + 0) * VECLEN_I);
        B1 = loadu64(a + (j + 1) * VECLEN_I);
        B2 = loadu64(a + (j + 2) * VECLEN_I);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            if ((i * BLOCKWORDS_I + j) > n)
                continue;

            __m512i A0 = loadu64(b + (i * BLOCKWORDS_I + 0) * VECLEN_I);
            __m512i A1 = loadu64(b + (i * BLOCKWORDS_I + 1) * VECLEN_I);
            __m512i A2 = loadu64(b + (i * BLOCKWORDS_I + 2) * VECLEN_I);
            __m512i A3 = loadu64(b + (i * BLOCKWORDS_I + 3) * VECLEN_I);

            VEC_MUL4_ACCUMX(A0, B0, B1, B2, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, B1, B2, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, B1, B2, zero, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, B1, B2, zero, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

            A0 = loadu64(c + ((i)*BLOCKWORDS_I + j + 0) * VECLEN_I);
            A1 = loadu64(c + ((i)*BLOCKWORDS_I + j + 1) * VECLEN_I);
            A2 = loadu64(c + ((i)*BLOCKWORDS_I + j + 2) * VECLEN_I);
            A3 = loadu64(c + ((i)*BLOCKWORDS_I + j + 3) * VECLEN_I);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 0) * VECLEN_I,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 1) * VECLEN_I,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 2) * VECLEN_I,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 3) * VECLEN_I,
                _mm512_and_si512(lo3, MASK));

            hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

            // subtract out all of the bias at once.
            SUB_BIAS_HI3(3, 2, 1);
            SUB_BIAS_LO3(3, 2, 1);

            A0 = loadu64(c + ((i)*BLOCKWORDS_I + j + 4) * VECLEN_I);
            A1 = loadu64(c + ((i)*BLOCKWORDS_I + j + 5) * VECLEN_I);
            A2 = loadu64(c + ((i)*BLOCKWORDS_I + j + 6) * VECLEN_I);
            A3 = loadu64(c + ((i)*BLOCKWORDS_I + j + 7) * VECLEN_I);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 4) * VECLEN_I,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 5) * VECLEN_I,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 6) * VECLEN_I,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 7) * VECLEN_I,
                _mm512_and_si512(hi6, MASK));

            if (((i)*BLOCKWORDS_I + j + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 8) * VECLEN_I,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i)*BLOCKWORDS_I + j + 8) * VECLEN_I),
                        _mm512_srli_epi64(hi6, 52)));
            }
        }
    }

    if (n - (full_blocks * BLOCKWORDS_I) == 2)
    {
        j = full_blocks * BLOCKWORDS_I;
        __m512i B0 = loadu64(b + (j + 0) * VECLEN_I);
        __m512i B1 = loadu64(b + (j + 1) * VECLEN_I);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            if ((i * BLOCKWORDS_I + j) > n)
                continue;

            __m512i A0 = loadu64(a + (i * BLOCKWORDS_I + 0) * VECLEN_I);
            __m512i A1 = loadu64(a + (i * BLOCKWORDS_I + 1) * VECLEN_I);
            __m512i A2 = loadu64(a + (i * BLOCKWORDS_I + 2) * VECLEN_I);
            __m512i A3 = loadu64(a + (i * BLOCKWORDS_I + 3) * VECLEN_I);

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

            A0 = loadu64(c + ((i)*BLOCKWORDS_I + j + 0) * VECLEN_I);
            A1 = loadu64(c + ((i)*BLOCKWORDS_I + j + 1) * VECLEN_I);
            A2 = loadu64(c + ((i)*BLOCKWORDS_I + j + 2) * VECLEN_I);
            A3 = loadu64(c + ((i)*BLOCKWORDS_I + j + 3) * VECLEN_I);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 0) * VECLEN_I,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 1) * VECLEN_I,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 2) * VECLEN_I,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 3) * VECLEN_I,
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

            A0 = loadu64(c + ((i)*BLOCKWORDS_I + j + 4) * VECLEN_I);
            A1 = loadu64(c + ((i)*BLOCKWORDS_I + j + 5) * VECLEN_I);
            A2 = loadu64(c + ((i)*BLOCKWORDS_I + j + 6) * VECLEN_I);
            A3 = loadu64(c + ((i)*BLOCKWORDS_I + j + 7) * VECLEN_I);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 4) * VECLEN_I,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 5) * VECLEN_I,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 6) * VECLEN_I,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 7) * VECLEN_I,
                _mm512_and_si512(hi6, MASK));

            if (((i)*BLOCKWORDS_I + j + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 8) * VECLEN_I,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i)*BLOCKWORDS_I + j + 8) * VECLEN_I),
                        _mm512_srli_epi64(hi6, 52)));
            }
        }

        B0 = loadu64(a + (j + 0) * VECLEN_I);
        B1 = loadu64(a + (j + 1) * VECLEN_I);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            if ((i * BLOCKWORDS_I + j) > n)
                continue;

            __m512i A0 = loadu64(b + (i * BLOCKWORDS_I + 0) * VECLEN_I);
            __m512i A1 = loadu64(b + (i * BLOCKWORDS_I + 1) * VECLEN_I);
            __m512i A2 = loadu64(b + (i * BLOCKWORDS_I + 2) * VECLEN_I);
            __m512i A3 = loadu64(b + (i * BLOCKWORDS_I + 3) * VECLEN_I);

            VEC_MUL4_ACCUMX(A0, B0, B1, zero, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, B1, zero, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, B1, zero, zero, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, B1, zero, zero, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

            A0 = loadu64(c + ((i)*BLOCKWORDS_I + j + 0) * VECLEN_I);
            A1 = loadu64(c + ((i)*BLOCKWORDS_I + j + 1) * VECLEN_I);
            A2 = loadu64(c + ((i)*BLOCKWORDS_I + j + 2) * VECLEN_I);
            A3 = loadu64(c + ((i)*BLOCKWORDS_I + j + 3) * VECLEN_I);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 0) * VECLEN_I,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 1) * VECLEN_I,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 2) * VECLEN_I,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 3) * VECLEN_I,
                _mm512_and_si512(lo3, MASK));

            hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

            // subtract out all of the bias at once.
            SUB_BIAS_HI3(3, 2, 1);
            SUB_BIAS_LO3(3, 2, 1);

            A0 = loadu64(c + ((i)*BLOCKWORDS_I + j + 4) * VECLEN_I);
            A1 = loadu64(c + ((i)*BLOCKWORDS_I + j + 5) * VECLEN_I);
            A2 = loadu64(c + ((i)*BLOCKWORDS_I + j + 6) * VECLEN_I);
            A3 = loadu64(c + ((i)*BLOCKWORDS_I + j + 7) * VECLEN_I);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 4) * VECLEN_I,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 5) * VECLEN_I,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 6) * VECLEN_I,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 7) * VECLEN_I,
                _mm512_and_si512(hi6, MASK));

            if (((i)*BLOCKWORDS_I + j + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 8) * VECLEN_I,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i)*BLOCKWORDS_I + j + 8) * VECLEN_I),
                        _mm512_srli_epi64(hi6, 52)));
            }
        }
    }

    if (n - (full_blocks * BLOCKWORDS_I) == 1)
    {
        j = full_blocks * BLOCKWORDS_I;
        __m512i B0 = loadu64(b + (j + 0) * VECLEN_I);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            if ((i * BLOCKWORDS_I + j) > n)
                continue;

            __m512i A0 = loadu64(a + (i * BLOCKWORDS_I + 0) * VECLEN_I);
            __m512i A1 = loadu64(a + (i * BLOCKWORDS_I + 1) * VECLEN_I);
            __m512i A2 = loadu64(a + (i * BLOCKWORDS_I + 2) * VECLEN_I);
            __m512i A3 = loadu64(a + (i * BLOCKWORDS_I + 3) * VECLEN_I);

            VEC_MUL4_ACCUMX(A0, B0, zero, zero, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, zero, zero, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, zero, zero, zero, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, zero, zero, zero, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

            A0 = loadu64(c + ((i)*BLOCKWORDS_I + j + 0) * VECLEN_I);
            A1 = loadu64(c + ((i)*BLOCKWORDS_I + j + 1) * VECLEN_I);
            A2 = loadu64(c + ((i)*BLOCKWORDS_I + j + 2) * VECLEN_I);
            A3 = loadu64(c + ((i)*BLOCKWORDS_I + j + 3) * VECLEN_I);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 0) * VECLEN_I,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 1) * VECLEN_I,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 2) * VECLEN_I,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 3) * VECLEN_I,
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

            A0 = loadu64(c + ((i)*BLOCKWORDS_I + j + 4) * VECLEN_I);
            A1 = loadu64(c + ((i)*BLOCKWORDS_I + j + 5) * VECLEN_I);
            A2 = loadu64(c + ((i)*BLOCKWORDS_I + j + 6) * VECLEN_I);
            A3 = loadu64(c + ((i)*BLOCKWORDS_I + j + 7) * VECLEN_I);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 4) * VECLEN_I,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 5) * VECLEN_I,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 6) * VECLEN_I,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 7) * VECLEN_I,
                _mm512_and_si512(hi6, MASK));

            if (((i)*BLOCKWORDS_I + j + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 8) * VECLEN_I,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i)*BLOCKWORDS_I + j + 8) * VECLEN_I),
                        _mm512_srli_epi64(hi6, 52)));
            }
        }

        B0 = loadu64(a + (j + 0) * VECLEN_I);

        for (i = 0; i < full_blocks; i++) {
            lo0 = lo1 = lo2 = lo3 = lo4 = lo5 = lo6 = zero;
            hi0 = hi1 = hi2 = hi3 = hi4 = hi5 = hi6 = zero;

            if ((i * BLOCKWORDS_I + j) > n)
                continue;

            __m512i A0 = loadu64(b + (i * BLOCKWORDS_I + 0) * VECLEN_I);
            __m512i A1 = loadu64(b + (i * BLOCKWORDS_I + 1) * VECLEN_I);
            __m512i A2 = loadu64(b + (i * BLOCKWORDS_I + 2) * VECLEN_I);
            __m512i A3 = loadu64(b + (i * BLOCKWORDS_I + 3) * VECLEN_I);

            VEC_MUL4_ACCUMX(A0, B0, zero, zero, zero, lo0, lo1, lo2, lo3, hi0, hi1, hi2, hi3);
            VEC_MUL4_ACCUMX(A1, B0, zero, zero, zero, lo1, lo2, lo3, lo4, hi1, hi2, hi3, hi4);
            VEC_MUL4_ACCUMX(A2, B0, zero, zero, zero, lo2, lo3, lo4, lo5, hi2, hi3, hi4, hi5);
            VEC_MUL4_ACCUMX(A3, B0, zero, zero, zero, lo3, lo4, lo5, lo6, hi3, hi4, hi5, hi6);

            // subtract out all of the bias at once.
            SUB_BIAS_HI4(1, 2, 3, 4);
            SUB_BIAS_LO4(1, 2, 3, 4);

            A0 = loadu64(c + ((i)*BLOCKWORDS_I + j + 0) * VECLEN_I);
            A1 = loadu64(c + ((i)*BLOCKWORDS_I + j + 1) * VECLEN_I);
            A2 = loadu64(c + ((i)*BLOCKWORDS_I + j + 2) * VECLEN_I);
            A3 = loadu64(c + ((i)*BLOCKWORDS_I + j + 3) * VECLEN_I);

            lo0 = _mm512_add_epi64(lo0, A0);
            lo1 = _mm512_add_epi64(lo1, A1);
            lo2 = _mm512_add_epi64(lo2, A2);
            lo3 = _mm512_add_epi64(lo3, A3);

            lo1 = _mm512_add_epi64(lo1, hi0);
            lo2 = _mm512_add_epi64(lo2, hi1);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 0) * VECLEN_I,
                _mm512_and_si512(lo0, MASK));

            hi0 = _mm512_srli_epi64(lo0, 52);
            lo1 = _mm512_add_epi64(lo1, hi0);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 1) * VECLEN_I,
                _mm512_and_si512(lo1, MASK));

            hi1 = _mm512_srli_epi64(lo1, 52);
            lo2 = _mm512_add_epi64(lo2, hi1);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 2) * VECLEN_I,
                _mm512_and_si512(lo2, MASK));

            hi2 = _mm512_srli_epi64(lo2, 52);
            lo3 = _mm512_add_epi64(lo3, hi2);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 3) * VECLEN_I,
                _mm512_and_si512(lo3, MASK));

            hi3 = _mm512_add_epi64(hi3, _mm512_srli_epi64(lo3, 52));

            // subtract out all of the bias at once.
            SUB_BIAS_HI3(3, 2, 1);
            SUB_BIAS_LO3(3, 2, 1);

            A0 = loadu64(c + ((i)*BLOCKWORDS_I + j + 4) * VECLEN_I);
            A1 = loadu64(c + ((i)*BLOCKWORDS_I + j + 5) * VECLEN_I);
            A2 = loadu64(c + ((i)*BLOCKWORDS_I + j + 6) * VECLEN_I);
            A3 = loadu64(c + ((i)*BLOCKWORDS_I + j + 7) * VECLEN_I);

            lo4 = _mm512_add_epi64(lo4, A0);
            lo5 = _mm512_add_epi64(lo5, A1);
            lo6 = _mm512_add_epi64(lo6, A2);

            lo4 = _mm512_add_epi64(lo4, hi3);
            lo5 = _mm512_add_epi64(lo5, hi4);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 4) * VECLEN_I,
                _mm512_and_si512(lo4, MASK));

            hi4 = _mm512_srli_epi64(lo4, 52);
            lo5 = _mm512_add_epi64(lo5, hi4);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 5) * VECLEN_I,
                _mm512_and_si512(lo5, MASK));

            hi5 = _mm512_srli_epi64(lo5, 52);
            lo6 = _mm512_add_epi64(lo6, hi5);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 6) * VECLEN_I,
                _mm512_and_si512(lo6, MASK));

            hi6 = _mm512_add_epi64(hi6, _mm512_srli_epi64(lo6, 52));
            hi6 = _mm512_add_epi64(A3, hi6);

            _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 7) * VECLEN_I,
                _mm512_and_si512(hi6, MASK));

            if (((i)*BLOCKWORDS_I + j + 8) < 2 * n)
            {
                _mm512_store_epi64(c + ((i)*BLOCKWORDS_I + j + 8) * VECLEN_I,
                    _mm512_add_epi64(
                        _mm512_load_epi64(c + ((i)*BLOCKWORDS_I + j + 8) * VECLEN_I),
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
    int i;
    int NWORDS = (int)mdata->NWORDS;

//#ifndef IFMA
    __m512d dbias = _mm512_castsi512_pd(set64(0x4670000000000000ULL));
    __m512i vbias1 = set64(0x4670000000000000ULL);  // 31
    UNUSED_VAR __m512i vbias2 = set64(0x4670000000000001ULL);  // 31
    UNUSED_VAR __m512i vbias3 = set64(0x4330000000000000ULL);  // 31
//#endif

    // needed after loops
    UNUSED_VAR __m512i zero = set64(0);
    __m512i MASK = set64(DIGIT_MASK);
    __m512i K = loadu64(mdata->vrho);
    __mmask8 scarry2;
    __mmask8 scarry;

    // deal with the sign
    c->size = (int)(uint32_t)NWORDS;
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
            __m512i a1 = _mm512_load_epi64(s->data + i * VECLEN_I);
            __m512i b0 = _mm512_load_epi64(n->data + i * VECLEN_I);
            __m512i a0 = _mm512_sbb_epi52(a1, scarry, b0, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(MASK, a0));
        }

        // negate any final borrows if there was also a final carry.
        scarry &= scarry2;

        // if there was a final borrow, we didn't need to do the subtraction after all.
        // replace with original results based on final borrow mask.
        for (i = NWORDS - 1; i >= 0; i--)
        {
            __m512i b0 = _mm512_load_epi64(s->data + i * VECLEN_I);
            _mm512_mask_store_epi64(c->data + i * VECLEN_I, scarry, b0);
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
    UNUSED_VAR int NBLOCKS;
    // needed in loops
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19

#ifndef IFMA
    UNUSED_VAR __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    UNUSED_VAR __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    UNUSED_VAR __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    UNUSED_VAR __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    UNUSED_VAR __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    UNUSED_VAR __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
#endif

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i zero = _mm512_set1_epi64(0);
    UNUSED_VAR __mmask8 scarry;

    __m512i arestore[3];
    __m512i brestore[3];

    NBLOCKS = words / BLOCKWORDS_I;
    if ((words % BLOCKWORDS_I) > 0)
    {
        NBLOCKS++;

#ifdef DEBUG_VECMUL
        printf("commencing vecmul52 with WORDS=%d and NBLOCKS=%d\n", words, NBLOCKS);
        print_vechex(a, 0, NBLOCKS * BLOCKWORDS_I, "input a: ");
        print_vechex(b, 0, NBLOCKS * BLOCKWORDS_I, "input b: ");
        print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS_I, "output c: ");
#endif

        for (i = words, j = 0; i < NBLOCKS * BLOCKWORDS_I; i++, j++)
        {
            arestore[j] = _mm512_load_epi64(a + i * VECLEN_I);
            brestore[j] = _mm512_load_epi64(b + i * VECLEN_I);
            _mm512_store_epi64(a + i * VECLEN_I, zero);
            _mm512_store_epi64(b + i * VECLEN_I, zero);
        }
    }

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

#ifdef DEBUG_VECMUL
    printf("after padding top words\n");
    print_vechex(a, 0, NBLOCKS * BLOCKWORDS_I, "input a: ");
    print_vechex(b, 0, NBLOCKS * BLOCKWORDS_I, "input b: ");
    print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS_I, "output c: ");
#endif

    // first half mul
    for (i = 0; i < NBLOCKS; i++)
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (j = i; j > 0; j--)
        {
            a0 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 3) * VECLEN_I);
            a1 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 2) * VECLEN_I);
            a2 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 1) * VECLEN_I);
            a3 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 0) * VECLEN_I);

            b0 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 1) * VECLEN_I);
            b1 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 2) * VECLEN_I);
            b2 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 3) * VECLEN_I);
            b3 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 4) * VECLEN_I);
            b4 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 5) * VECLEN_I);
            b5 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 6) * VECLEN_I);
            b6 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 7) * VECLEN_I);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum
        a0 = _mm512_load_epi64(a + (i * BLOCKWORDS_I + 0) * VECLEN_I);
        a1 = _mm512_load_epi64(a + (i * BLOCKWORDS_I + 1) * VECLEN_I);
        a2 = _mm512_load_epi64(a + (i * BLOCKWORDS_I + 2) * VECLEN_I);
        a3 = _mm512_load_epi64(a + (i * BLOCKWORDS_I + 3) * VECLEN_I);

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
            _mm512_store_epi64(c + (i * BLOCKWORDS_I + j) * VECLEN_I, acc_e0);

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
            _mm512_store_epi64(c + (i * BLOCKWORDS_I + j) * VECLEN_I, acc_e0);

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
            _mm512_store_epi64(c + (i * BLOCKWORDS_I + j) * VECLEN_I, acc_e0);

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
            _mm512_store_epi64(c + (i * BLOCKWORDS_I + j) * VECLEN_I, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;
        }
    }

#ifdef DEBUG_VECMUL
    print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS_I, "after lo half:");
#endif

    for (i = words, j = 0; i < NBLOCKS * BLOCKWORDS_I; i++, j++)
    {
        _mm512_store_epi64(a + i * VECLEN_I, arestore[j]);
        _mm512_store_epi64(b + i * VECLEN_I, brestore[j]);
    }


#ifdef DEBUG_VECMUL
    printf("after restoring inputs and outputs\n");
    print_vechex(a, 0, NBLOCKS * BLOCKWORDS_I, "input a: ");
    print_vechex(b, 0, NBLOCKS * BLOCKWORDS_I, "input b: ");
    print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS_I, "output c: ");
    exit(1);
#endif

    return;
}

void vecmul52_n(uint64_t* a, uint64_t* b, uint64_t* c, int words)
{
    int i, j;
    UNUSED_VAR int NBLOCKS;
    // needed in loops
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19

#ifndef IFMA
    UNUSED_VAR __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    UNUSED_VAR __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    UNUSED_VAR __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    UNUSED_VAR __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    UNUSED_VAR __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    UNUSED_VAR __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
#endif

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i zero = _mm512_set1_epi64(0);

    NBLOCKS = words / BLOCKWORDS_I;
    if ((words % BLOCKWORDS_I) > 0)
    {
        NBLOCKS++;

#ifdef DEBUG_VECMUL
        printf("commencing vecmul52 with WORDS=%d and NBLOCKS=%d\n", words, NBLOCKS);
        print_vechex(a, 0, NBLOCKS * BLOCKWORDS_I, "input a: ");
        print_vechex(b, 0, NBLOCKS * BLOCKWORDS_I, "input b: ");
        print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS_I, "output c: ");
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
            a0 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 3) * VECLEN_I);
            a1 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 2) * VECLEN_I);
            a2 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 1) * VECLEN_I);
            a3 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 0) * VECLEN_I);

            b0 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 1) * VECLEN_I);
            b1 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 2) * VECLEN_I);
            b2 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 3) * VECLEN_I);
            b3 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 4) * VECLEN_I);
            b4 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 5) * VECLEN_I);
            b5 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 6) * VECLEN_I);
            b6 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 7) * VECLEN_I);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum
        a0 = _mm512_load_epi64(a + (i * BLOCKWORDS_I + 0) * VECLEN_I);
        a1 = _mm512_load_epi64(a + (i * BLOCKWORDS_I + 1) * VECLEN_I);
        a2 = _mm512_load_epi64(a + (i * BLOCKWORDS_I + 2) * VECLEN_I);
        a3 = _mm512_load_epi64(a + (i * BLOCKWORDS_I + 3) * VECLEN_I);

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
            _mm512_store_epi64(c + (i * BLOCKWORDS_I + j) * VECLEN_I, acc_e0);

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
            _mm512_store_epi64(c + (i * BLOCKWORDS_I + j) * VECLEN_I, acc_e0);

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
            _mm512_store_epi64(c + (i * BLOCKWORDS_I + j) * VECLEN_I, acc_e0);

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
            _mm512_store_epi64(c + (i * BLOCKWORDS_I + j) * VECLEN_I, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;
        }
    }

#ifdef DEBUG_VECMUL
    print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS_I, "after lo half:");
#endif

    // second half mul
    for (i = NBLOCKS; i < 2 * NBLOCKS; i++)
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (j = i - NBLOCKS + 1; j < NBLOCKS; j++)
        {
            a0 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 3) * VECLEN_I);
            a1 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 2) * VECLEN_I);
            a2 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 1) * VECLEN_I);
            a3 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 0) * VECLEN_I);

            b0 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 1) * VECLEN_I);
            b1 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 2) * VECLEN_I);
            b2 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 3) * VECLEN_I);
            b3 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 4) * VECLEN_I);
            b4 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 5) * VECLEN_I);
            b5 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 6) * VECLEN_I);
            b6 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 7) * VECLEN_I);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }


        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a + ((i - NBLOCKS) * BLOCKWORDS_I + 1) * VECLEN_I);
        a2 = _mm512_load_epi64(a + ((i - NBLOCKS) * BLOCKWORDS_I + 2) * VECLEN_I);
        a3 = _mm512_load_epi64(a + ((i - NBLOCKS) * BLOCKWORDS_I + 3) * VECLEN_I);

        b0 = _mm512_load_epi64(b + ((NBLOCKS * BLOCKWORDS_I) - 1) * VECLEN_I);
        b1 = _mm512_load_epi64(b + ((NBLOCKS * BLOCKWORDS_I) - 2) * VECLEN_I);
        b2 = _mm512_load_epi64(b + ((NBLOCKS * BLOCKWORDS_I) - 3) * VECLEN_I);

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
            //_mm512_store_pd(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I, 
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
            //_mm512_store_pd(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I,
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
            //_mm512_store_pd(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I,
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
            //_mm512_store_pd(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            a3 = _mm512_and_epi64(vlmask, a3);
            a2 = _mm512_and_epi64(vlmask, a2);
            a1 = _mm512_and_epi64(vlmask, a1);
            a0 = _mm512_and_epi64(vlmask, a0);

            _mm512_store_epi64(c + (i * BLOCKWORDS_I + 0) * VECLEN_I, a3);
            _mm512_store_epi64(c + (i * BLOCKWORDS_I + 1) * VECLEN_I, a2);
            _mm512_store_epi64(c + (i * BLOCKWORDS_I + 2) * VECLEN_I, a1);
            _mm512_store_epi64(c + (i * BLOCKWORDS_I + 3) * VECLEN_I, a0);

        }
    }

#ifdef DEBUG_VECMUL
    print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS_I, "after hi half:");
#endif

    return;
}

void vecmul52_np(uint64_t* a, uint64_t* b, uint64_t* c, int words)
{
    int i, j;
    UNUSED_VAR int NBLOCKS;
    // needed in loops
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19

#ifndef IFMA
    UNUSED_VAR __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    UNUSED_VAR __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    UNUSED_VAR __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    UNUSED_VAR __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    UNUSED_VAR __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    UNUSED_VAR __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
#endif

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i zero = _mm512_set1_epi64(0);
    UNUSED_VAR __mmask8 scarry;

    __m512i ahi[BLOCKWORDS_I];
    __m512i bhi[BLOCKWORDS_I];

    NBLOCKS = words / BLOCKWORDS_I;
    if ((words % BLOCKWORDS_I) == 0)
    {
        printf("error only call vecmul52_np when 4 does not divide n\n");
    }

    for (i = NBLOCKS * BLOCKWORDS_I, j = 0; i < words; i++, j++)
    {
        ahi[j] = _mm512_load_epi64(a + i * VECLEN_I);
        bhi[j] = _mm512_load_epi64(b + i * VECLEN_I);
    }
    for (; j < BLOCKWORDS_I; j++)
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

            a0 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 3) * VECLEN_I);
            a1 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 2) * VECLEN_I);
            a2 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 1) * VECLEN_I);
            a3 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 0) * VECLEN_I);

            b0 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 1) * VECLEN_I);
            b1 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 2) * VECLEN_I);
            b2 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 3) * VECLEN_I);
            
            if (j == (NBLOCKS-1))
            {
                b3 = bhi[0];
                b4 = bhi[1];
                b5 = bhi[2];
                b6 = bhi[3];
            }
            else
            {
                b3 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 4) * VECLEN_I);
                b4 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 5) * VECLEN_I);
                b5 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 6) * VECLEN_I);
                b6 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 7) * VECLEN_I);
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
            a0 = _mm512_load_epi64(a + (i * BLOCKWORDS_I + 0) * VECLEN_I);
            a1 = _mm512_load_epi64(a + (i * BLOCKWORDS_I + 1) * VECLEN_I);
            a2 = _mm512_load_epi64(a + (i * BLOCKWORDS_I + 2) * VECLEN_I);
            a3 = _mm512_load_epi64(a + (i * BLOCKWORDS_I + 3) * VECLEN_I);
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
            _mm512_store_epi64(c + (i * BLOCKWORDS_I + j) * VECLEN_I, acc_e0);

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
            _mm512_store_epi64(c + (i * BLOCKWORDS_I + j) * VECLEN_I, acc_e0);

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
            _mm512_store_epi64(c + (i * BLOCKWORDS_I + j) * VECLEN_I, acc_e0);

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
            _mm512_store_epi64(c + (i * BLOCKWORDS_I + j) * VECLEN_I, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;
        }
    }

#ifdef DEBUG_VECMUL
    print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS_I, "after lo half:");
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
                a0 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 3) * VECLEN_I);
                a1 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 2) * VECLEN_I);
                a2 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 1) * VECLEN_I);
                a3 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 0) * VECLEN_I);
            }

            b0 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 1) * VECLEN_I);
            b1 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 2) * VECLEN_I);
            b2 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 3) * VECLEN_I);
            if (j == (NBLOCKS - 1))
            {
                b3 = bhi[0];
                b4 = bhi[1];
                b5 = bhi[2];
                b6 = bhi[3];
            }
            else
            {
                b3 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 4) * VECLEN_I);
                b4 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 5) * VECLEN_I);
                b5 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 6) * VECLEN_I);
                b6 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 7) * VECLEN_I);
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
            a1 = _mm512_load_epi64(a + ((i - NBLOCKS) * BLOCKWORDS_I + 1) * VECLEN_I);
            a2 = _mm512_load_epi64(a + ((i - NBLOCKS) * BLOCKWORDS_I + 2) * VECLEN_I);
            a3 = _mm512_load_epi64(a + ((i - NBLOCKS) * BLOCKWORDS_I + 3) * VECLEN_I);
        }

        b0 = bhi[3];
        b1 = bhi[2];
        b2 = bhi[1];

        //b0 = _mm512_load_epi64(b + ((NBLOCKS * BLOCKWORDS_I) - 1) * VECLEN_I);
        //b1 = _mm512_load_epi64(b + ((NBLOCKS * BLOCKWORDS_I) - 2) * VECLEN_I);
        //b2 = _mm512_load_epi64(b + ((NBLOCKS * BLOCKWORDS_I) - 3) * VECLEN_I);

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
            //_mm512_store_pd(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I, 
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
            //_mm512_store_pd(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I,
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
            //_mm512_store_pd(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I,
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
            //_mm512_store_pd(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            a3 = _mm512_and_epi64(vlmask, a3);
            a2 = _mm512_and_epi64(vlmask, a2);
            a1 = _mm512_and_epi64(vlmask, a1);
            a0 = _mm512_and_epi64(vlmask, a0);

            if ((i * BLOCKWORDS_I + 0) < words)
                _mm512_store_epi64(c + (i * BLOCKWORDS_I + 0) * VECLEN_I, a3);
            if ((i * BLOCKWORDS_I + 1) < words)
                _mm512_store_epi64(c + (i * BLOCKWORDS_I + 1) * VECLEN_I, a2);
            if ((i * BLOCKWORDS_I + 2) < words)
                _mm512_store_epi64(c + (i * BLOCKWORDS_I + 2) * VECLEN_I, a1);
            if ((i * BLOCKWORDS_I + 3) < words)
                _mm512_store_epi64(c + (i * BLOCKWORDS_I + 3) * VECLEN_I, a0);

        }
    }

#ifdef DEBUG_VECMUL
    print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS_I, "after hi half:");
#endif


    return;
}

void vecsqr52_n(uint64_t* a, uint64_t* c, int words)
{
    int i, j, k;
    uint64_t* b = a;
    UNUSED_VAR int NBLOCKS;
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19

#ifndef IFMA
    UNUSED_VAR __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    UNUSED_VAR __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    UNUSED_VAR __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    UNUSED_VAR __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    UNUSED_VAR __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    UNUSED_VAR __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
#endif

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i zero = _mm512_set1_epi64(0);
    UNUSED_VAR __mmask8 scarry;

    __m512i arestore[3];
    __m512i crestore[6];

    NBLOCKS = words / BLOCKWORDS_I;
    if ((words % BLOCKWORDS_I) > 0)
    {
        NBLOCKS++;

#ifdef DEBUG_VECMUL
        printf("commencing vecsqr52 with WORDS=%d and NBLOCKS=%d\n", words, NBLOCKS);
        print_vechex(a, 0, NBLOCKS * BLOCKWORDS_I, "input a: ");
        print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS_I, "output c: ");
#endif

        for (i = words, j = 0; i < NBLOCKS * BLOCKWORDS_I; i++, j++)
        {
            arestore[j] = _mm512_load_epi64(a + i * VECLEN_I);
            _mm512_store_epi64(a + i * VECLEN_I, zero);
        }

        for (i = 2 * words, j = 0; i < 2 * NBLOCKS * BLOCKWORDS_I; i++, j++)
        {
            crestore[j] = _mm512_load_epi64(c + i * VECLEN_I);
            _mm512_store_epi64(c + i * VECLEN_I, zero);
        }
    }

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

#ifdef DEBUG_VECMUL
    printf("after padding top words\n");
    print_vechex(a, 0, NBLOCKS * BLOCKWORDS_I, "input a: ");
    print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS_I, "output c: ");
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
            a0 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 3) * VECLEN_I);
            a1 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 2) * VECLEN_I);
            a2 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 1) * VECLEN_I);
            a3 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 0) * VECLEN_I);

            b0 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 1) * VECLEN_I);
            b1 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 2) * VECLEN_I);
            b2 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 3) * VECLEN_I);
            b3 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 4) * VECLEN_I);
            b4 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 5) * VECLEN_I);
            b5 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 6) * VECLEN_I);
            b6 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 7) * VECLEN_I);

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

            a0 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 3) * VECLEN_I);
            a1 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 2) * VECLEN_I);
            a2 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 1) * VECLEN_I);
            a3 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 0) * VECLEN_I);

            b1 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 2) * VECLEN_I);
            b2 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 3) * VECLEN_I);
            b3 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 4) * VECLEN_I);
            b4 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 5) * VECLEN_I);
            b5 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 6) * VECLEN_I);
            b6 = _mm512_load_epi64(b + ((j - 1) * BLOCKWORDS_I + 7) * VECLEN_I);

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
            a0 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 0) * VECLEN_I);
            a1 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 1) * VECLEN_I);
            a2 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 2) * VECLEN_I);
            a3 = _mm512_load_epi64(a + ((i - j) * BLOCKWORDS_I + 3) * VECLEN_I);

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
            _mm512_store_epi64(c + (i * BLOCKWORDS_I + j) * VECLEN_I, acc_e0);

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
            _mm512_store_epi64(c + (i * BLOCKWORDS_I + j) * VECLEN_I, acc_e0);

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
            _mm512_store_epi64(c + (i * BLOCKWORDS_I + j) * VECLEN_I, acc_e0);

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
            _mm512_store_epi64(c + (i * BLOCKWORDS_I + j) * VECLEN_I, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;
        }
    }

#ifdef DEBUG_VECMUL
    print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS_I, "after lo half:");
#endif

    // second half sqr
    for (i = 0; i < NBLOCKS; i++)
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (k = 0, j = 0; j < (NBLOCKS - i - 1) / 2; j++, k++)
        {
            // Compute a solid block (all matching terms are in the lower
            // half triangle of the expansion).
            a0 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS_I) - 1 - j * BLOCKWORDS_I) * VECLEN_I);
            a1 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS_I) - 2 - j * BLOCKWORDS_I) * VECLEN_I);
            a2 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS_I) - 3 - j * BLOCKWORDS_I) * VECLEN_I);
            a3 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS_I) - 4 - j * BLOCKWORDS_I) * VECLEN_I);

            b0 = _mm512_load_epi64(b + ((j + i) * BLOCKWORDS_I + 1) * VECLEN_I);
            b1 = _mm512_load_epi64(b + ((j + i) * BLOCKWORDS_I + 2) * VECLEN_I);
            b2 = _mm512_load_epi64(b + ((j + i) * BLOCKWORDS_I + 3) * VECLEN_I);
            b3 = _mm512_load_epi64(b + ((j + i) * BLOCKWORDS_I + 4) * VECLEN_I);
            b4 = _mm512_load_epi64(b + ((j + i) * BLOCKWORDS_I + 5) * VECLEN_I);
            b5 = _mm512_load_epi64(b + ((j + i) * BLOCKWORDS_I + 6) * VECLEN_I);
            b6 = _mm512_load_epi64(b + ((j + i) * BLOCKWORDS_I + 7) * VECLEN_I);

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
                a0 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS_I) - 1 - j * BLOCKWORDS_I) * VECLEN_I);
                a1 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS_I) - 2 - j * BLOCKWORDS_I) * VECLEN_I);
                a2 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS_I) - 3 - j * BLOCKWORDS_I) * VECLEN_I);
                a3 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS_I) - 4 - j * BLOCKWORDS_I) * VECLEN_I);

                b0 = _mm512_load_epi64(b + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 1) * VECLEN_I);
                b1 = _mm512_load_epi64(b + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 2) * VECLEN_I);
                b2 = _mm512_load_epi64(b + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 3) * VECLEN_I);
                b3 = _mm512_load_epi64(b + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 4) * VECLEN_I);
                b4 = _mm512_load_epi64(b + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 5) * VECLEN_I);

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
                a0 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS_I) - 1 - j * BLOCKWORDS_I) * VECLEN_I);
                a1 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS_I) - 2 - j * BLOCKWORDS_I) * VECLEN_I);
                a2 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS_I) - 3 - j * BLOCKWORDS_I) * VECLEN_I);

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
                a0 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS_I) - 1 - j * BLOCKWORDS_I) * VECLEN_I);  // {f, b}
                a1 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS_I) - 2 - j * BLOCKWORDS_I) * VECLEN_I);  // {e, a}
                a2 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS_I) - 3 - j * BLOCKWORDS_I) * VECLEN_I);  // {d, 9}

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
                a0 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS_I) - 1 - j * BLOCKWORDS_I) * VECLEN_I);		// {f, b}
                a1 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS_I) - 2 - j * BLOCKWORDS_I) * VECLEN_I);		// {e, a}
                a2 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS_I) - 3 - j * BLOCKWORDS_I) * VECLEN_I);		// {d, 9}
                a3 = _mm512_load_epi64(a + ((NBLOCKS * BLOCKWORDS_I) - 4 - j * BLOCKWORDS_I) * VECLEN_I);		// {c, 8}

                b0 = _mm512_load_epi64(b + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 1) * VECLEN_I); // {9, 5}
                b1 = _mm512_load_epi64(b + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 2) * VECLEN_I);	// {a, 6}
                b2 = _mm512_load_epi64(b + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 3) * VECLEN_I);	// {b, 7}
                b3 = _mm512_load_epi64(b + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 4) * VECLEN_I);	// {c, 8}
                b4 = _mm512_load_epi64(b + (j * BLOCKWORDS_I + i * BLOCKWORDS_I + 5) * VECLEN_I); // {d, 9}

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
            //_mm512_store_pd(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I, 
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
            //_mm512_store_pd(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I,
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
            //_mm512_store_pd(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I,
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
            //_mm512_store_pd(s->data + (i * BLOCKWORDS_I + j) * VECLEN_I,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            a3 = _mm512_and_epi64(vlmask, a3);
            a2 = _mm512_and_epi64(vlmask, a2);
            a1 = _mm512_and_epi64(vlmask, a1);
            a0 = _mm512_and_epi64(vlmask, a0);

            _mm512_store_epi64(c + ((NBLOCKS + i) * BLOCKWORDS_I + 0) * VECLEN_I, a3);
            _mm512_store_epi64(c + ((NBLOCKS + i) * BLOCKWORDS_I + 1) * VECLEN_I, a2);
            _mm512_store_epi64(c + ((NBLOCKS + i) * BLOCKWORDS_I + 2) * VECLEN_I, a1);
            _mm512_store_epi64(c + ((NBLOCKS + i) * BLOCKWORDS_I + 3) * VECLEN_I, a0);

        }
    }

#ifdef DEBUG_VECMUL
    print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS_I, "after hi half:");
#endif

    for (i = words, j = 0; i < NBLOCKS * BLOCKWORDS_I; i++, j++)
    {
        _mm512_store_epi64(a + i * VECLEN_I, arestore[j]);
    }

    for (i = 2 * words, j = 0; i < 2 * NBLOCKS * BLOCKWORDS_I; i++, j++)
    {
        _mm512_store_epi64(c + i * VECLEN_I, crestore[j]);
    }

#ifdef DEBUG_VECMUL
    printf("after restoring inputs and outputs\n");
    print_vechex(a, 0, NBLOCKS* BLOCKWORDS_I, "input a: ");
    print_vechex(c, 0, 2 * NBLOCKS * BLOCKWORDS_I, "output c: ");
    //exit(1);
#endif

    return;
}

void kmuln(uint64_t* a, uint64_t* b, uint64_t* c, uint64_t* scratch, int words);
void kmuln_lo(uint64_t* a, uint64_t* b, uint64_t* c, uint64_t* scratch, int words);
void ksqrn(uint64_t* a, uint64_t* c, uint64_t* scratch, int words);
/* Base-case threshold of the veckmul dispatch: widths <= this go straight to
   the block kernels (vecmul52_n at %4==0, else vecmul52_cios), wider widths
   recurse into kmuln.  Measured on the dev Xeon (plans/toom3_optimize.md
   4c): base beats 1-level karatsuba at every width 21..28, so 28 beats the
   old 20 -- but the balance is machine-dependent.  The standalone harness
   (plans/bench/thresh/, mode "thresh") measures the optimum per machine;
   override here with -DVECKMUL_BASE_THRESH=n when the sweep says so.
   NOTE: vecksqr and vectmul keep their own (20) thresholds -- ksqrn-vs-base
   is mixed in 21..32 and toom-3 is hopeless below ~30 words. */
#ifndef VECKMUL_BASE_THRESH
#define VECKMUL_BASE_THRESH 28
#endif
#define veckmul(a, b, c, s, words)				\
  do {									\
    if (words <= VECKMUL_BASE_THRESH) {			\
        if ((words & 0x3) == 0) { vecmul52_n(a, b, c, (int)words); } \
        else { vecmul52_cios(a, b, c, (int)words); } }     \
    else								\
        kmuln(a, b, c, s, (int)words);				\
  } while (0);

#define veckmul_lo(a, b, c, s, words)				\
  do {									\
    if (words <= VECKMUL_BASE_THRESH) {			\
        if ((words & 0x3) == 0) { vecmul52_lo(a, b, c, (int)words); } \
        else { vecmul52_cios_lo(a, b, c, (int)words); } }     \
    else								\
        kmuln_lo(a, b, c, s, (int)words);				\
  } while (0);

#define vecksqr(a, c, s, words)				\
  do {									\
    if (words <= 20) {						\
        if ((words & 0x3) == 0) { vecsqr52_n(a, c, (int)words); } \
        else { vecmul52_cios(a, a, c, (int)words); } }     \
    else								\
        ksqrn(a, c, s, (int)words);				\
  } while (0);

#define vectmul(a, b, c, s, words)				\
  do {									\
    if (words <= 20) {						\
        if ((words & 0x3) == 0) { vecmul52_n(a, b, c, (int)words); } \
        else { vecmul52_cios(a, b, c, (int)words); } }     \
    else								\
        tmuln(a, b, c, s, (int)words);				\
  } while (0);

/* exact in-place division by 3 of each vec_bignum in u (L limbs each).
   LSB-first Hensel division: q_i = (r_i * inv3) mod B, then subtract
   3*q_i*B^i -- by the divisibility invariant this touches only limb i
   (which becomes zero, no borrow) and limb i+1 (hi = 3*q_i >> 52, no
   borrow because the remaining value is a non-negative multiple of
   B^(i+1)).  Requires 3 | u.  Sign of u must be tracked by the caller. */
#define TOOM_INV3 3002399751580331ULL   /* (2^53+1)/3: 3*INV3 = 2^53+1 = 1 mod 2^52 */

void vec_hdiv3_52(uint64_t* u, int L)
{
    int i;
    __m512i inv3 = _mm512_set1_epi64(TOOM_INV3);
    __m512i lomask = _mm512_set1_epi64(0xfffffffffffffULL);
    UNUSED_VAR __m512i r, q, prod, lo, hi, n;

    for (i = 0; i < L; i++)
    {
        r = _mm512_load_epi64(u + i * VECLEN_I);
        q = _mm512_and_epi64(lomask, _mm512_mullo_epi64(r, inv3));
        prod = _mm512_add_epi64(q, _mm512_add_epi64(q, q));      /* 3*q < 2^54 */
        lo = _mm512_and_epi64(lomask, prod);
        hi = _mm512_srli_epi64(prod, 52);                        /* 0, 1 or 2 */
#ifdef DEBUG_MERSENNE
        if (_mm512_kor(0, _mm512_mask_cmpneq_epu64_mask(0xff, r, lo)) != 0)
        {
            printf("tmuln: hdiv3 limb %d not divisible by 3\n", i);
            exit(1);
        }
#endif
        _mm512_store_epi64(u + i * VECLEN_I, q);   /* write quotient limb; residual limb i is exactly 0 */
        if (i + 1 < L)
        {
            n = _mm512_load_epi64(u + (i + 1) * VECLEN_I);
#ifdef DEBUG_MERSENNE
            if (_mm512_kor(0, _mm512_mask_cmpgt_epu64_mask(0xff, hi, n)) != 0)
            {
                printf("tmuln: hdiv3 unexpected borrow at limb %d\n", i);
                exit(1);
            }
#endif
            _mm512_store_epi64(u + (i + 1) * VECLEN_I, _mm512_sub_epi64(n, hi));
        }
    }
    return;
}

/* exact in-place 1-bit right shift of each vec_bignum in u (L limbs).
   Requires u even.  Sign of u must be tracked by the caller. */
void vec_rshift1_52(uint64_t* u, int L)
{
    int i;
    __m512i lomask = _mm512_set1_epi64(0xfffffffffffffULL);
    __m512i one = _mm512_set1_epi64(1);
    __m512i lo, hi, r;

    for (i = 0; i < L; i++)
    {
        lo = _mm512_load_epi64(u + i * VECLEN_I);
        hi = (i + 1 < L) ? _mm512_load_epi64(u + (i + 1) * VECLEN_I)
                         : _mm512_setzero_si512();
        r = _mm512_or_epi64(_mm512_srli_epi64(lo, 1),
                            _mm512_slli_epi64(_mm512_and_epi64(hi, one), 51));
#ifdef DEBUG_MERSENNE
        /* exact /2 requires the whole value to be even, i.e. limb 0 even */
        if (i == 0 && _mm512_kor(0, _mm512_mask_test_epi64_mask(0xff, lo, one)) != 0)
        {
            printf("tmuln: rshift1 of odd value\n");
            exit(1);
        }
#endif
        _mm512_store_epi64(u + i * VECLEN_I, _mm512_and_epi64(lomask, r));
    }
    return;
}

/* signed add of two vec_bignum arrays with per-lane sign masks:
   c = (a, aneg) + (b, bneg), where a set neg bit means that lane's value
   is negative.  c may alias a or b.  Returns the per-lane sign mask of
   the result (set = negative; zero-magnitude lanes may keep a stale sign
   bit -- canonicalize before testing).  W = limb width; c needs W+1
   limbs of space (a carry limb is written at index W). */
__mmask8 base_sadd_52(uint64_t* a, __mmask8 aneg, uint64_t* b, __mmask8 bneg,
    uint64_t* c, int W)
{
    __mmask8 opp = aneg ^ bneg;
    __mmask8 addmask = ~opp;
    __mmask8 abssubmask = (__mmask8)(opp & vec_gte2_52(b, a, W));
    UNUSED_VAR __mmask8 acarry = base_absaddsub_52(a, b, c, addmask, abssubmask, W, W);
#ifdef DEBUG_MERSENNE
    if (acarry)
    {
        printf("base_sadd_52: unexpected carry out of %d limbs\n", W);
        exit(1);
    }
#endif
    return (aneg & addmask) | (opp & ((abssubmask & bneg) | (~abssubmask & aneg)));
}

/* fused add + subtract (GMP's mpn_add_n_sub_n pattern, used by its toom
   evaluators): one bottom-up pass writes
   sum = a + b  (carry limb written at sum[W], like base_add_52_wc) and
   diff = a - b  (raw difference; negative in the lanes of the returned
   borrow mask, i.e. where b > a).  a and b need W limbs; sum needs W+1,
   diff W.  Follow with base_neg_52(diff, borrow, W) to get the
   magnitude |a - b| in the borrow lanes (magnitude + sign-mask pair, as
   base_abssub_52 would produce).  One full-width sweep fewer than the
   add + gte + abssub triple. */
__mmask8 base_addsub_52(uint64_t* a, uint64_t* b, uint64_t* sum, uint64_t* diff, int W)
{
    int i;
    __m512i avec, bvec, svec, dvec;
    __mmask8 scarry = 0, sborrow = 0;

    for (i = 0; i < W; i++)
    {
        avec = _mm512_load_epi64(a + i * VECLEN_I);
        bvec = _mm512_load_epi64(b + i * VECLEN_I);

        svec = _mm512_adc_epi52(avec, scarry, bvec, &scarry);
        dvec = _mm512_sbb_epi52(avec, sborrow, bvec, &sborrow);

        _mm512_store_epi64(sum + i * VECLEN_I, svec);
        _mm512_store_epi64(diff + i * VECLEN_I, dvec);
    }
    _mm512_store_epi64(sum + i * VECLEN_I, _mm512_mask_mov_epi64(
        _mm512_setzero_si512(), scarry, _mm512_set1_epi64(1)));

    return sborrow;     /* set where b > a: the difference is negative */
}

/* conditional negate: in the lanes of sm, dst = -(dst) mod B^W (two's
   complement over the W limbs); the other lanes pass through unchanged.
   Turns the raw difference of base_addsub_52 into the magnitude of the
   per-lane signed value. */
void base_neg_52(uint64_t* dst, __mmask8 sm, int W)
{
    int i;
    __m512i lomask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i d, t, carry;
    __m512i czero = _mm512_maskz_set1_epi64(sm, 1);     /* +1 in sm lanes */

    carry = czero;
    for (i = 0; i < W; i++)
    {
        d = _mm512_load_epi64(dst + i * VECLEN_I);
        /* ~limb + carry: ~d is lomask ^ d for limbs < B */
        t = _mm512_add_epi64(_mm512_xor_epi64(d, lomask), carry);
        _mm512_store_epi64(dst + i * VECLEN_I, _mm512_mask_mov_epi64(d, sm,
            _mm512_and_epi64(lomask, t)));
        carry = _mm512_srli_epi64(t, 52);
    }

    return;
}

/* signed sublsh1 (GMP's as2 = 2*(a2 + as1) - a0 fusion in toom33_mul.c,
   threaded through the per-lane sign masks tmuln's magnitude+sign
   convention needs): dst = 2*(dst, dsign) - (src, all lanes negative).
   The signed sum is accumulated in reduced radix over W+1 limbs: each
   negated term enters complemented (XOR each limb with B-1) with its +1
   at limb 0, the B^(W+1) multiples land in the carry dropped past limb
   W, and a lane's result is negative exactly when its dropped carry
   equals the number of complemented terms minus one (= dsign).  The
   caller fixes the negative lanes up with base_neg_52(dst, mask, W+1).
   |result| < B^W is the caller's guarantee; dst needs W+1 limbs. */
__mmask8 base_sublsh1_52(uint64_t* dst, __mmask8 dsign, uint64_t* src, int W)
{
    int i;
    __m512i lomask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i zero = _mm512_setzero_si512();
    __m512i one = _mm512_set1_epi64(1);
    __m512i dsignv = _mm512_maskz_set1_epi64(dsign, 1);  /* value 1 in dsign lanes */
    __m512i d, dig, t, v, svec, c2 = zero, cp = zero;
    __mmask8 m;

    for (i = 0; i <= W; i++)
    {
        /* double the dst digit (limb W holds the incoming carry limb) */
        d = _mm512_load_epi64(dst + i * VECLEN_I);
        d = _mm512_add_epi64(d, d);
        d = _mm512_add_epi64(d, c2);
        dig = _mm512_and_epi64(lomask, d);
        c2 = _mm512_srli_epi64(d, 52);

        /* +2*dst, complemented in the dsign lanes */
        t = _mm512_mask_mov_epi64(dig, dsign, _mm512_xor_epi64(dig, lomask));

        /* - src (complemented; zero-extended past its W limbs) */
        svec = (i < W) ? _mm512_load_epi64(src + i * VECLEN_I) : zero;
        t = _mm512_add_epi64(t, _mm512_xor_epi64(svec, lomask));

        /* the two complement +1's at limb 0 (src always, dst in dsign) */
        v = _mm512_add_epi64(t, cp);
        if (i == 0)
        {
            v = _mm512_add_epi64(v, one);
            v = _mm512_add_epi64(v, dsignv);
        }

        _mm512_store_epi64(dst + i * VECLEN_I, _mm512_and_epi64(lomask, v));
        cp = _mm512_srli_epi64(v, 52);
    }

    /* negative exactly where the dropped carry equals dsign */
    m = _mm512_cmpeq_epu64_mask(cp, dsignv);
    return m;
}

/* signed addlsh1 (GMP's mpn_sublsh1_n_ip1 pattern in
   toom_interpolate_5pts.c, threaded through the per-lane sign masks):
   dst = (dst, dsign) + 2*(src, all lanes non-negative).  src is NOT
   modified.  Same reduced-radix scheme as base_sublsh1_52: the signed
   sum runs over W+1 limbs, complemented terms enter with their +1 at
   limb 0, and the B^(W+1) multiples land in the carry dropped past limb
   W; a dsign lane's result is negative exactly when its dropped carry
   is zero.  The caller fixes the negative lanes up with
   base_neg_52(dst, mask, W+1).  |result| < B^W is the caller's
   guarantee; dst needs W+1 limbs. */
__mmask8 base_addlsh1s_52(uint64_t* dst, __mmask8 dsign, uint64_t* src, int W)
{
    int i;
    __m512i lomask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i zero = _mm512_setzero_si512();
    UNUSED_VAR __m512i one = _mm512_set1_epi64(1);
    __m512i d, dig, t, v, svec, c2 = zero, cp = zero;
    __mmask8 m;

    for (i = 0; i <= W; i++)
    {
        /* double the src digit on the fly */
        svec = _mm512_load_epi64(src + i * VECLEN_I);
        d = _mm512_add_epi64(svec, svec);
        d = _mm512_add_epi64(d, c2);
        dig = _mm512_and_epi64(lomask, d);
        c2 = _mm512_srli_epi64(d, 52);

        /* + dst, complemented in the dsign lanes */
        d = _mm512_load_epi64(dst + i * VECLEN_I);
        t = _mm512_mask_mov_epi64(d, dsign, _mm512_xor_epi64(d, lomask));

        v = _mm512_add_epi64(t, dig);
        v = _mm512_add_epi64(v, cp);
        if (i == 0)
            v = _mm512_add_epi64(v, _mm512_maskz_set1_epi64(dsign, 1));

        _mm512_store_epi64(dst + i * VECLEN_I, _mm512_and_epi64(lomask, v));
        cp = _mm512_srli_epi64(v, 52);
    }

    /* a dsign lane is negative exactly when its dropped carry is zero */
    m = dsign & _mm512_cmpeq_epu64_mask(cp, zero);
    return m;
}

/* signed add of two vec_bignum arrays with per-lane sign masks, followed
   by an exact division by two written directly to c (shift-on-write):
   c = ((a, aneg) + (b, bneg)) / 2, where the sum must be even.  Fuses
   base_sadd_52 + vec_rshift1_52 and saves the separate full-width rshift
   sweep: limb i of the result only needs bit 0 of the sum's limb i+1,
   which is computed as the pipeline's next step.  c needs W limbs of
   space; the sum must be < 2*B^W (result fits W limbs).  Sign of the
   result returned as in base_sadd_52. */
__mmask8 base_sadd_rshift1_52(uint64_t* a, __mmask8 aneg, uint64_t* b, __mmask8 bneg,
    uint64_t* c, int W)
{
    int i;
    __m512i lomask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i one = _mm512_set1_epi64(1);
    __m512i avec, bvec, tmp, cur, nxt;
    __mmask8 opp = aneg ^ bneg;
    __mmask8 addmask = ~opp;
    __mmask8 abssubmask = (__mmask8)(opp & vec_gte2_52(b, a, W));
    __mmask8 scarry = 0, acarry = 0;

    /* pipeline: cur = signed sum of limb i, nxt = signed sum of limb i+1
       (or the sum's carry limb at the top); limb i of the halved result
       is (cur >> 1) | (nxt bit 0) << 51. */
    avec = _mm512_load_epi64(a);
    bvec = _mm512_load_epi64(b);

    /* swap words in the lanes we are abs-subbing */
    tmp = avec;
    avec = _mm512_mask_mov_epi64(avec, abssubmask, bvec);
    bvec = _mm512_mask_mov_epi64(bvec, abssubmask, tmp);

    cur = _mm512_mask_sbb_epi52(avec, ~addmask, scarry, bvec, &scarry);
    cur = _mm512_mask_mov_epi64(cur, addmask,
        _mm512_mask_adc_epi52(avec, addmask, acarry, bvec, &acarry));

    for (i = 0; i < W; i++)
    {
        if (i + 1 < W)
        {
            avec = _mm512_load_epi64(a + (i + 1) * VECLEN_I);
            bvec = _mm512_load_epi64(b + (i + 1) * VECLEN_I);

            tmp = avec;
            avec = _mm512_mask_mov_epi64(avec, abssubmask, bvec);
            bvec = _mm512_mask_mov_epi64(bvec, abssubmask, tmp);

            nxt = _mm512_mask_sbb_epi52(avec, ~addmask, scarry, bvec, &scarry);
            nxt = _mm512_mask_mov_epi64(nxt, addmask,
                _mm512_mask_adc_epi52(avec, addmask, acarry, bvec, &acarry));
        }
        else
        {
            /* the sum's carry limb: 1 in the lanes that added with carry */
            nxt = _mm512_maskz_set1_epi64(acarry, 1);
        }

        _mm512_store_epi64(c + i * VECLEN_I, _mm512_and_epi64(lomask,
            _mm512_or_epi64(_mm512_srli_epi64(cur, 1),
                _mm512_slli_epi64(_mm512_and_epi64(nxt, one), 51))));
        cur = nxt;
    }

    return (aneg & addmask) | (opp & ((abssubmask & bneg) | (~abssubmask & aneg)));
}

/* signed 3-input add with per-lane sign masks:
   c = (a, aneg) + (b, bneg) + (d, dneg).  Sums the positive and negative
   groups separately in one bottom-up pass (into caller-provided
   temporaries p and n, each with W+1 limbs of space), compares the group
   magnitudes (top-down, carry limbs included), then subtracts the
   smaller group from the larger -- one full-width sweep fewer than two
   base_sadd_52 calls.  c may alias a, b, or d.  Returns the sign mask of
   the result (set = negative; zero-magnitude lanes may keep a stale sign
   bit -- canonicalize before testing). */
__mmask8 base_sadd3_52(uint64_t* a, __mmask8 aneg, uint64_t* b, __mmask8 bneg,
    uint64_t* d, __mmask8 dneg, uint64_t* c, int W, uint64_t* p, uint64_t* n)
{
    int i;
    __m512i lomask = _mm512_set1_epi64(0x000fffffffffffffULL);
    UNUSED_VAR __m512i one = _mm512_set1_epi64(1);
    __m512i av, pacc, nacc, pc, nc, m1, m2, m3, pcarry, ncarry;
    __m512i ones = _mm512_set1_epi64(~0LL);
    __mmask8 ngt, plt, m;

    /* per-lane membership masks (all-ones in the lanes where the term
       belongs to the positive group; membership is per-lane and constant
       across limbs, so non-member lanes never accumulate and their carry
       is always zero) */
    m1 = _mm512_maskz_mov_epi64(~aneg, ones);
    m2 = _mm512_maskz_mov_epi64(~bneg, ones);
    m3 = _mm512_maskz_mov_epi64(~dneg, ones);

    pcarry = _mm512_setzero_si512();
    ncarry = _mm512_setzero_si512();
    for (i = 0; i < W; i++)
    {
        av = _mm512_load_epi64(a + i * VECLEN_I);
        pacc = _mm512_and_epi64(av, m1);
        nacc = _mm512_andnot_epi64(m1, av);

        av = _mm512_load_epi64(b + i * VECLEN_I);
        pacc = _mm512_add_epi64(pacc, _mm512_and_epi64(av, m2));
        nacc = _mm512_add_epi64(nacc, _mm512_andnot_epi64(m2, av));

        av = _mm512_load_epi64(d + i * VECLEN_I);
        pacc = _mm512_add_epi64(pacc, _mm512_and_epi64(av, m3));
        nacc = _mm512_add_epi64(nacc, _mm512_andnot_epi64(m3, av));

        /* fold in the previous limb's carry once (0 in non-member lanes) */
        pacc = _mm512_add_epi64(pacc, pcarry);
        nacc = _mm512_add_epi64(nacc, ncarry);

        _mm512_store_epi64(p + i * VECLEN_I, _mm512_and_epi64(lomask, pacc));
        pcarry = _mm512_srli_epi64(pacc, 52);

        _mm512_store_epi64(n + i * VECLEN_I, _mm512_and_epi64(lomask, nacc));
        ncarry = _mm512_srli_epi64(nacc, 52);
    }

    _mm512_store_epi64(p + i * VECLEN_I, pcarry);
    _mm512_store_epi64(n + i * VECLEN_I, ncarry);

    /* compare the group sums (including their carry limbs):
       m = (n >= p) per lane; where set the result is negative */
    pc = _mm512_load_epi64(p + i * VECLEN_I);
    nc = _mm512_load_epi64(n + i * VECLEN_I);
    ngt = _mm512_cmpgt_epu64_mask(nc, pc);
    plt = _mm512_cmpgt_epu64_mask(pc, nc);
    m = (__mmask8)(ngt | (~plt & vec_gte2_52(n, p, W)));

    base_absaddsub_52(p, n, c, 0, m, W, W);     /* |p - n|, direction m */
    return m;
}

/* optional per-phase cycle instrumentation for tmuln (zero cost when
   TMULN_PROFILE is not defined).  Accumulates rdtsc deltas for the eval,
   pointwise, interpolation and recomposition phases across calls;
   tmuln_profile_reset() zeroes the accumulators and tmuln_profile_report()
   prints the split. */
#ifdef TMULN_PROFILE
#include <x86intrin.h>
static unsigned long long tmuln_prof_calls = 0;
static unsigned long long tmuln_prof_setup = 0;
static unsigned long long tmuln_prof_eval = 0;
static unsigned long long tmuln_prof_pw = 0;
static unsigned long long tmuln_prof_interp = 0;
static unsigned long long tmuln_prof_recomp = 0;
static unsigned long long tmuln_prof_total = 0;

void tmuln_profile_reset(void)
{
    tmuln_prof_calls = 0;
    tmuln_prof_setup = 0;
    tmuln_prof_eval = 0;
    tmuln_prof_pw = 0;
    tmuln_prof_interp = 0;
    tmuln_prof_recomp = 0;
    tmuln_prof_total = 0;
    return;
}

void tmuln_profile_report(void)
{
    if (tmuln_prof_calls == 0)
    {
        printf("tmuln_profile: no calls\n");
        return;
    }
    printf("tmuln_profile: calls=%llu total=%.1f cyc/op\n"
        "  setup=%.1f (%.1f%%) eval=%.1f (%.1f%%) pointwise=%.1f (%.1f%%) "
        "interp=%.1f (%.1f%%) recomp=%.1f (%.1f%%)\n",
        tmuln_prof_calls, (double)tmuln_prof_total / tmuln_prof_calls,
        (double)tmuln_prof_setup / tmuln_prof_calls,
        100.0 * tmuln_prof_setup / tmuln_prof_total,
        (double)tmuln_prof_eval / tmuln_prof_calls,
        100.0 * tmuln_prof_eval / tmuln_prof_total,
        (double)tmuln_prof_pw / tmuln_prof_calls,
        100.0 * tmuln_prof_pw / tmuln_prof_total,
        (double)tmuln_prof_interp / tmuln_prof_calls,
        100.0 * tmuln_prof_interp / tmuln_prof_total,
        (double)tmuln_prof_recomp / tmuln_prof_calls,
        100.0 * tmuln_prof_recomp / tmuln_prof_total);
    return;
}
#endif

void tmuln(uint64_t* a, uint64_t* b, uint64_t* c, uint64_t* scratch, int words)
{
    /*
     * Toom-3 multiplication of words-limb vec_bignums, VECLEN at a time.
     * Split: a = a0 + a1*X + a2*X^2 (w0 = w1 = words/3 limbs each, X = B^w0),
     * evaluate at 0, 1, -1, -2, inf, multiply pointwise (Karatsuba), then
     * interpolate and recompose.  All signed arithmetic uses magnitudes plus
     * one sign bit per lane, exactly like kcombine2.  The algorithm was
     * validated with a scalar per-lane transcription against schoolbook
     * multiplication (1M+ cases, plans/toom3_finish.md).
     *
     * Scratch layout (offsets in vector-words; one 'word' = VECLEN_I limbs):
     *   [0, K)                     recursive Karatsuba scratch (K = 4*W+8)
     *   a0p..b2p (6 * W)           padded copies of the six split pieces
     *   pm1, p1, pm2, qm1, q1, qm2 evaluation points (W+1 each)
     *   r0, r1, rm1, rm2, rinf     pointwise products (P+1 each, P = 2*W)
     *   i1, i2, i3                 interpolation (P+1 each)
     * The region above K is zeroed once, before any recursion; the
     * recursive calls only touch [0, K).
     */
    int i;
    int w0 = words / 3, w1 = words / 3, w2 = words - 2 * w0;
    int wmax = MAX(w0, w2) + 1;
    /* eval width: the evaluation values are bounded by 9*B^w2 < B^wmax
       (|p(-2)| < 8*B^w2 + B^w0), so wmax limbs + the primitives' carry
       limb suffice; the old +2 guard pads the sweeps but is never needed
       for the wc carries.  P = 2*W covers the interpolation values, which
       are bounded by 3*B^(P/2-1) < B^P. */
    int W = wmax;                 /* eval padded width */
    int P = 2 * W;                /* pointwise/interpolation padded width */
    int A = 2 * words + 2;        /* recomposition accumulator width */
    int K = 4 * W + 8;            /* recursive karatsuba scratch (vector-words) */
    int o_a0p, o_a1p, o_a2p, o_b0p, o_b1p, o_b2p;
    int o_pm1, o_p1, o_pm2, o_qm1, o_q1, o_qm2;
    int o_r0, o_r1, o_rm1, o_rm2, o_rinf;
    int o_i1, o_i2, o_i3;
    int T;
    uint64_t *a0p, *a1p, *a2p, *b0p, *b1p, *b2p;
    uint64_t *pm1, *p1, *pm2, *qm1, *q1, *qm2;
    uint64_t *r0, *r1, *rm1, *rm2, *rinf;
    uint64_t *i1, *i2, *i3;
    __mmask8 m1, m2, m3, m4, tsign, qsign, rm1sign, rm2sign;
    __mmask8 i1sign, i2sign, i3sign;
    __m512i zero = _mm512_setzero_si512();
#ifdef TMULN_PROFILE
    unsigned long long tp0, tp1, tp2, tp3, tp4, tp5;
    tp0 = __rdtsc();
#endif

    o_a0p = K;
    o_a1p = o_a0p + W;
    o_a2p = o_a1p + W;
    o_b0p = o_a2p + W;
    o_b1p = o_b0p + W;
    o_b2p = o_b1p + W;
    o_pm1 = o_b2p + W;
    o_p1 = o_pm1 + W + 1;
    o_pm2 = o_p1 + W + 1;
    o_qm1 = o_pm2 + W + 1;
    o_q1 = o_qm1 + W + 1;
    o_qm2 = o_q1 + W + 1;
    o_r0 = o_qm2 + W + 1;
    o_r1 = o_r0 + P + 1;
    o_rm1 = o_r1 + P + 1;
    o_rm2 = o_rm1 + P + 1;
    o_rinf = o_rm2 + P + 1;
    o_i1 = o_rinf + P + 1;
    o_i2 = o_i1 + P + 1;
    o_i3 = o_i2 + P + 1;
    T = o_i3 + P + 1;

    /* mtmp4 is sized in vec_monty_alloc from the same worst-case formula:
       toomT = 34*(words/3+5) + 23 + 2*words + 2 vector-words, guaranteed
       via vecInit((toomT+1)/2) -> 2*ceil(toomT/2)+4 >= toomT+4. */
    {
        int toomWg = words / 3 + 5;
        int toomTg = 34 * toomWg + 23 + 2 * words + 2;
        if (T > toomTg + 4)
        {
            printf("tmuln: scratch overflow: need %d have %d\n", T, toomTg + 4);
            exit(1);
        }
    }

    a0p = scratch + o_a0p * VECLEN_I;
    a1p = scratch + o_a1p * VECLEN_I;
    a2p = scratch + o_a2p * VECLEN_I;
    b0p = scratch + o_b0p * VECLEN_I;
    b1p = scratch + o_b1p * VECLEN_I;
    b2p = scratch + o_b2p * VECLEN_I;
    pm1 = scratch + o_pm1 * VECLEN_I;
    p1 = scratch + o_p1 * VECLEN_I;
    pm2 = scratch + o_pm2 * VECLEN_I;
    qm1 = scratch + o_qm1 * VECLEN_I;
    q1 = scratch + o_q1 * VECLEN_I;
    qm2 = scratch + o_qm2 * VECLEN_I;
    r0 = scratch + o_r0 * VECLEN_I;
    r1 = scratch + o_r1 * VECLEN_I;
    rm1 = scratch + o_rm1 * VECLEN_I;
    rm2 = scratch + o_rm2 * VECLEN_I;
    rinf = scratch + o_rinf * VECLEN_I;
    i1 = scratch + o_i1 * VECLEN_I;
    i2 = scratch + o_i2 * VECLEN_I;
    i3 = scratch + o_i3 * VECLEN_I;

    /* zero everything above the recursive scratch; safe here because the
        pointwise recursive calls below only touch [0, K). */
    memset(scratch + K * VECLEN_I, 0, (unsigned long)((T - K) * VECLEN_I) * sizeof(uint64_t));

    /* padded copies of the six split pieces */
    memcpy(a0p, a, (unsigned long)(w0 * VECLEN_I) * sizeof(uint64_t));
    memcpy(a1p, a + w0 * VECLEN_I, (unsigned long)(w1 * VECLEN_I) * sizeof(uint64_t));
    memcpy(a2p, a + 2 * w0 * VECLEN_I, (unsigned long)(w2 * VECLEN_I) * sizeof(uint64_t));
    memcpy(b0p, b, (unsigned long)(w0 * VECLEN_I) * sizeof(uint64_t));
    memcpy(b1p, b + w0 * VECLEN_I, (unsigned long)(w1 * VECLEN_I) * sizeof(uint64_t));
    memcpy(b2p, b + 2 * w0 * VECLEN_I, (unsigned long)(w2 * VECLEN_I) * sizeof(uint64_t));

#ifdef TMULN_PROFILE
    tp1 = __rdtsc();    /* end of setup (memset + split copies) */
#endif

    /* ---- evaluation phase ----
       p(0) = a0, p(inf) = a2, p(1) = a0+a1+a2,
       p(-1) = a0+a2-a1, p(-2) = 2*(p(-1)+a2)-a0  (signs kept per lane).
       The +1/-1 pair is one fused add+sub pass (GMP's mpn_add_n_sub_n
       pattern) plus a conditional negate for the borrow lanes. */
    base_add_52_wc(a0p, a2p, pm1, W, W);            /* pm1 = a0 + a2 */
    m1 = base_addsub_52(pm1, a1p, p1, pm1, W);      /* p1 = gp+a1; pm1 = gp-a1 */
    if (m1)
        base_neg_52(pm1, m1, W);                    /* pm1 = |p(-1)|, sign m1 */

    tsign = base_sadd_52(pm1, m1, a2p, 0, pm2, W);  /* pm2 = p(-1) + a2 (signed) */
    m2 = base_sublsh1_52(pm2, tsign, a0p, W);       /* p(-2) = 2*pm2 - a0 */
    if (m2)
        base_neg_52(pm2, m2, W + 1);                /* magnitude in the negative lanes */

    base_add_52_wc(b0p, b2p, qm1, W, W);            /* qm1 = b0 + b2 */
    base_add_52_wc(qm1, b1p, q1, W, W);             /* q1 = q(1) */
    m3 = (__mmask8)vec_gte2_52(b1p, qm1, W);                  /* b1 >= qm1: q(-1) <= 0 */
    base_abssub_52(qm1, b1p, qm1, m3, W, W);        /* qm1 = |q(-1)|, sign m3 */

    qsign = base_sadd_52(qm1, m3, b2p, 0, qm2, W);
    m4 = base_sublsh1_52(qm2, qsign, b0p, W);       /* q(-2) = 2*qm2 - b0 */
    if (m4)
        base_neg_52(qm2, m4, W + 1);

#ifdef TMULN_PROFILE
    tp2 = __rdtsc();    /* end of evaluation */
#endif

    /* ---- pointwise multiplication (Karatsuba above the base case) ---- */
    veckmul(a0p, b0p, r0, scratch, W);
    veckmul(p1, q1, r1, scratch, W);
    veckmul(pm1, qm1, rm1, scratch, W);
    veckmul(pm2, qm2, rm2, scratch, W);
    veckmul(a2p, b2p, rinf, scratch, W);
    rm1sign = m1 ^ m3;      /* product of signed values */
    rm2sign = m2 ^ m4;
#ifdef TMULN_PROFILE
    tp3 = __rdtsc();    /* end of pointwise */
#endif

    /* ---- interpolation (mirrors the proven scalar sequence) ----
       i3 = (rm2 - r1)/3;  i1 = (r1 - rm1)/2;  i2 = rm1 - r0;
       i3 = (i2 - i3)/2 + 2*rinf;  i2 = i2 + i1 - rinf;  i1 = i1 - i3.
       Fused primitives: the /2 shifts are folded into the preceding
       signed adds (shift-on-write), the doubling of rinf skips its
       (always equal) magnitude comparison, and i2 += i1 - rinf is a
       single 3-input signed pass. */
    i3sign = base_sadd_52(rm2, rm2sign, r1, 0xff, i3, P);   /* i3 = rm2 - r1 */
    vec_hdiv3_52(i3, P);                                    /* /3, exact */

    i1sign = base_sadd_rshift1_52(r1, 0, rm1, ~rm1sign, i1, P);  /* i1 = (r1 - rm1)/2 */

    i2sign = base_sadd_52(rm1, rm1sign, r0, 0xff, i2, P);   /* i2 = rm1 - r0 */

    i3sign = base_sadd_rshift1_52(i2, i2sign, i3, ~i3sign, i3, P);  /* i3 = (i2 - i3)/2 */

    i3sign = base_addlsh1s_52(i3, i3sign, rinf, P);         /* i3 += 2*rinf (rinf unmodified) */
    if (i3sign)
        base_neg_52(i3, i3sign, P + 1);                     /* magnitude in the negative lanes */

    i2sign = base_sadd3_52(i2, i2sign, i1, i1sign, rinf, 0xff, i2, P, r1, rm2);  /* i2 += i1 - rinf (temps r1/rm2 are dead here; r0 still needed by the recomposition) */
    i1sign = base_sadd_52(i1, i1sign, i3, ~i3sign, i1, P);  /* i1 -= i3 */
#ifdef TMULN_PROFILE
    tp4 = __rdtsc();    /* end of interpolation */
#endif

    /* ---- recomposition ----
       c = i0 + i1*B^w0 + i2*B^2w0 + i3*B^3w0 + i4*B^4w0, with i4 = the
       original rinf (2*w2 significant limbs: 4*w0 + 2*w2 = 2*words exactly).
       Instead of four (memset + offset copy + full-width signed add)
       cycles, accumulate everything in one reduced-radix pass, kcombine2
       style: each signed term is added as-is or complemented per lane
       (XOR each limb with B-1, plus a +1 at the term's bottom limb for
       the negated lanes).  A negated P-limb term at offset o represents
       -term*B^o mod B^A with complement digits (B-1) in [o+P, A), so the
       loop below adds B-1 in the negated lanes past each term's window.
       The sum of all digit contributions is then exactly the true value
       plus (number of negated lanes)*B^A, and since the product is
       non-negative and < B^(2*words) < B^A, limbs [0, A) of the exact
       digit-wise sum are the product (the B^A multiples land in the
       dropped carry out of limb A-1). */
    memcpy(c, r0, (unsigned long)(P * VECLEN_I) * sizeof(uint64_t));         /* acc = i0 */
    memset(c + P * VECLEN_I, 0, (unsigned long)((A - P) * VECLEN_I) * sizeof(uint64_t));

    {
        __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
        __m512i vone = _mm512_set1_epi64(1);
        __m512i vB1 = _mm512_set1_epi64(0x000fffffffffffffULL);  /* B-1 */
        __m512i v, t, vc = zero;
        uint64_t *p1_in = i1, *p2 = i2, *p3 = i3, *p4 = rinf;

        for (i = 0; i < A; i++)
        {
            v = _mm512_add_epi64(_mm512_load_epi64(c + i * VECLEN_I), vc);

            /* i1*B^w0, negated in lanes of i1sign */
            if (i >= w0)
            {
                if (i < w0 + P)
                {
                    t = _mm512_load_epi64(p1_in);
                    p1_in += VECLEN_I;
                    t = _mm512_mask_xor_epi64(t, i1sign, t, vlmask);
                    v = _mm512_add_epi64(v, t);
                    if (i == w0)
                        v = _mm512_mask_add_epi64(v, i1sign, v, vone);
                }
                else
                {
                    v = _mm512_mask_add_epi64(v, i1sign, v, vB1);
                }
            }
            /* i2*B^2w0, negated in lanes of i2sign */
            if (i >= 2 * w0)
            {
                if (i < 2 * w0 + P)
                {
                    t = _mm512_load_epi64(p2);
                    p2 += VECLEN_I;
                    t = _mm512_mask_xor_epi64(t, i2sign, t, vlmask);
                    v = _mm512_add_epi64(v, t);
                    if (i == 2 * w0)
                        v = _mm512_mask_add_epi64(v, i2sign, v, vone);
                }
                else
                {
                    v = _mm512_mask_add_epi64(v, i2sign, v, vB1);
                }
            }
            /* i3*B^3w0, negated in lanes of i3sign */
            if (i >= 3 * w0)
            {
                if (i < 3 * w0 + P)
                {
                    t = _mm512_load_epi64(p3);
                    p3 += VECLEN_I;
                    t = _mm512_mask_xor_epi64(t, i3sign, t, vlmask);
                    v = _mm512_add_epi64(v, t);
                    if (i == 3 * w0)
                        v = _mm512_mask_add_epi64(v, i3sign, v, vone);
                }
                else
                {
                    v = _mm512_mask_add_epi64(v, i3sign, v, vB1);
                }
            }
            /* rinf*B^4w0, always non-negative */
            if (i >= 4 * w0)
            {
                v = _mm512_add_epi64(v, _mm512_load_epi64(p4));
                p4 += VECLEN_I;
            }

            /* carry-propagate and store this limb of the exact sum */
            _mm512_store_epi64(c + i * VECLEN_I, _mm512_and_epi64(vlmask, v));
            vc = _mm512_srli_epi64(v, 52);
        }
    }

    /* the product is non-negative and < B^(2*words): limbs [2*words, A)
       of the exact sum must be zero */
#ifdef DEBUG_MERSENNE
    {
        __mmask8 hi_limb = 0;
        for (i = 2 * words; i < A; i++)
            hi_limb |= _mm512_cmpneq_epu64_mask(_mm512_load_epi64(c + i * VECLEN_I), zero);
        if (hi_limb)
        {
            printf("tmuln: product exceeds 2*words limbs\n");
            exit(1);
        }
    }
#endif
#ifdef TMULN_PROFILE
    tp5 = __rdtsc();
    tmuln_prof_calls++;
    tmuln_prof_setup += tp1 - tp0;
    tmuln_prof_eval += tp2 - tp1;
    tmuln_prof_pw += tp3 - tp2;
    tmuln_prof_interp += tp4 - tp3;
    tmuln_prof_recomp += tp5 - tp4;
    tmuln_prof_total += tp5 - tp0;
#endif

    return;
}


void kmuln(uint64_t* a, uint64_t* b, uint64_t* c, uint64_t* scratch, int words)
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
    uint64_t* a1 = a + lowords * VECLEN_I;
    uint64_t* b0 = b;
    uint64_t* b1 = b + lowords * VECLEN_I;
    uint64_t* z1;
    uint64_t* s1;
    uint64_t* s2;
    __mmask8 m1;
    __mmask8 m2;

    s1 = scratch;
    s2 = scratch + (lowords)*VECLEN_I;
    z1 = scratch + (lowords * 2) * VECLEN_I;

    veckmul(a0, b0, c, scratch, lowords);
    veckmul(a1, b1, c + 2 * lowords * VECLEN_I, scratch, hiwords);
    // result is negative if the longer subtractand does not have a leading digit
    // AND if the rest of the comparison is true.
    m1 = (__mmask8)(hiwords < lowords ? (_mm512_cmpeq_epu64_mask(_mm512_load_epi64(a0 + hiwords * VECLEN_I),
        _mm512_setzero_si512()) & vec_gte2_52(a1, a0, hiwords)) : vec_gte2_52(a1, a0, hiwords));
    // result is negative if the longer subtractand has a leading digit
    // OR if the rest of the comparison is true.
    m2 = (__mmask8)(hiwords < lowords ? (_mm512_cmpgt_epu64_mask(_mm512_load_epi64(b0 + hiwords * VECLEN_I),
        _mm512_setzero_si512()) | vec_gte2_52(b0, b1, hiwords)) : vec_gte2_52(b0, b1, hiwords));
    base_abssub_52(a0, a1, s1, m1, lowords, hiwords);
    base_abssub_52(b1, b0, s2, m2, hiwords, lowords);
    veckmul(s1, s2, z1, scratch + (4 * lowords) * VECLEN_I, lowords);

    kcombine2(c, z1, s1, m1 ^ m2, 2 * lowords, hiwords);
    return;
}

void kmuln_lo(uint64_t* a, uint64_t* b, uint64_t* c, uint64_t* scratch, int words)
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
    uint64_t* a1 = a + lowords * VECLEN_I;
    uint64_t* b0 = b;
    uint64_t* b1 = b + lowords * VECLEN_I;
    uint64_t* z1;
    uint64_t* s1;
    uint64_t* s2;
    __mmask8 m1;
    __mmask8 m2;

    s1 = scratch;
    s2 = scratch + (lowords)*VECLEN_I;
    z1 = scratch + (lowords * 2) * VECLEN_I;

    veckmul(a0, b0, c, scratch, lowords);
    veckmul_lo(a1, b1, c + 2 * lowords * VECLEN_I, scratch, hiwords);
    // result is negative if the longer subtractand does not have a leading digit
    // AND if the rest of the comparison is true.
    m1 = (__mmask8)(hiwords < lowords ? (_mm512_cmpeq_epu64_mask(_mm512_load_epi64(a0 + hiwords * VECLEN_I),
        _mm512_setzero_si512()) & vec_gte2_52(a1, a0, hiwords)) : vec_gte2_52(a1, a0, hiwords));
    // result is negative if the longer subtractand has a leading digit
    // OR if the rest of the comparison is true.
    m2 = (__mmask8)(hiwords < lowords ? (_mm512_cmpgt_epu64_mask(_mm512_load_epi64(b0 + hiwords * VECLEN_I),
        _mm512_setzero_si512()) | vec_gte2_52(b0, b1, hiwords)) : vec_gte2_52(b0, b1, hiwords));
    base_abssub_52(a0, a1, s1, m1, lowords, hiwords);
    base_abssub_52(b1, b0, s2, m2, hiwords, lowords);
    veckmul_lo(s1, s2, z1, scratch + (4 * lowords) * VECLEN_I, lowords);

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
    uint64_t* a1 = a + lowords * VECLEN_I;
    uint64_t* z1;
    uint64_t* s1;
    UNUSED_VAR uint64_t* s2;
    __mmask8 m1;
    UNUSED_VAR int i;

    if (a == c)
    {
        printf("error output pointer cannot be the same as input pointer\n");
        exit(1);
    }

    s1 = scratch;
    s2 = scratch + (lowords)*VECLEN_I;
    z1 = scratch + (lowords * 2) * VECLEN_I;

    //printf("ksqr lo len %d words: ", lowords);
    //for (i = lowords - 1; i >= 0; i--)
    //{
    //    printf("%013llx", a0[i * VECLEN_I]);
    //}
    //printf("\n");

    vecksqr(a0, c, scratch, lowords);

    //printf("result: ");
    //for (i = lowords * 2 - 1; i >= 0; i--)
    //{
    //    printf("%013llx", c[i * VECLEN_I]);
    //}
    //printf("\n");
    //
    //printf("ksqr hi len %d words: ", hiwords);
    //for (i = hiwords - 1; i >= 0; i--)
    //{
    //    printf("%013llx", a1[i * VECLEN_I]);
    //}
    //printf("\n");

    vecksqr(a1, c + 2 * lowords * VECLEN_I, scratch, hiwords);

    //printf("result: ");
    //for (i = lowords * 2 - 1; i >= 0; i--)
    //{
    //    printf("%013llx", c[(2 * lowords + i) * VECLEN_I]);
    //}
    //printf("\n");
    //
    //printf("a0: ");
    //for (i = lowords - 1; i >= 0; i--)
    //{
    //    printf("%013llx", a0[i * VECLEN_I]);
    //}
    //printf("\n");
    //
    //printf("a1: ");
    //for (i = hiwords - 1; i >= 0; i--)
    //{
    //    printf("%013llx", a1[i * VECLEN_I]);
    //}
    //printf("\n");

    // result is negative if the longer subtractand does not have a leading digit
    // AND if the rest of the comparison is true.
    m1 = (__mmask8)(hiwords < lowords ? (_mm512_cmpeq_epu64_mask(_mm512_load_epi64(a0 + hiwords * VECLEN_I),
        _mm512_setzero_si512()) & vec_gte2_52(a1, a0, hiwords)) : vec_gte2_52(a1, a0, hiwords));
    base_abssub_52(a0, a1, s1, m1, lowords, hiwords);

    //printf("ksqr absdiff len %d words: ", lowords);
    //for (i = lowords - 1; i >= 0; i--)
    //{
    //    printf("%013llx", s1[i * VECLEN_I]);
    //}
    //printf("\n");

    vecksqr(s1, z1, scratch + (4 * lowords) * VECLEN_I, lowords);

    //printf("combine at offset %d words\n", lowords);
    //
    //printf("z0=0x");
    //for (i = lowords * 2 - 1; i >= 0; i--)
    //{
    //    printf("%013llx", c[i * VECLEN_I]);
    //}
    //printf(";\n");
    //
    //printf("z1=0x");
    //for (i = lowords * 2 - 1; i >= 0; i--)
    //{
    //    printf("%013llx", z1[i * VECLEN_I]);
    //}
    //printf(";\n");
    //
    //printf("z2=0x");
    //for (i = lowords * 2 - 1; i >= 0; i--)
    //{
    //    printf("%013llx", c[(2 * lowords + i) * VECLEN_I]);
    //}
    //printf(";\n");
    //
    //printf("s1=0x");
    //for (i = lowords * 2 - 1; i >= 0; i--)
    //{
    //    printf("%013llx", s1[i * VECLEN_I]);
    //}
    //printf(";\n");

    kcombine2(c, z1, s1, 0xff, 2 * lowords, hiwords);

    //printf("result: ");
    //for (i = lowords * 4 - 1; i >= 0; i--)
    //{
    //    printf("%013llx", c[i * VECLEN_I]);
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
    // a, b, c are of length VECLEN_I * NWORDS
    // s1 is of length VECLEN_I
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
        avec = _mm512_load_epi64(a + i * VECLEN_I);
        bvec = _mm512_load_epi64(b + i * VECLEN_I);
        cvec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c + i * VECLEN_I,
            _mm512_and_epi64(_mm512_set1_epi64(VEC_MAXDIGIT), cvec));
    }

    if (carry)
    {
        // subtract any final borrows that exceed the size of b.
        _mm512_mask_store_epi64(c + i * VECLEN_I, carry, _mm512_set1_epi64(1));
    }

    return carry;
}

void vecksqr_redc(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    vecksqr(a->data, mdata->mtmp1->data, mdata->mtmp4->data, mdata->NWORDS);

    veckmul_lo(mdata->mtmp1->data, mdata->vnhat->data, mdata->mtmp2->data, mdata->mtmp4->data, mdata->NWORDS);
    veckmul(mdata->mtmp2->data, n->data, mdata->mtmp3->data, mdata->mtmp4->data, mdata->NWORDS);

    mdata->mtmp1->size = (int)(mdata->NWORDS * 2);
    mdata->mtmp3->size = (int)(mdata->NWORDS * 2);
    uint32_t m = vec_bignum52_special_add(mdata->mtmp1->data + mdata->NWORDS * VECLEN_I,
        mdata->mtmp3->data + mdata->NWORDS * VECLEN_I, c->data, (int)mdata->NWORDS);

#ifndef USE_AMM
    c->WORDS_ALLOC = (uint32_t)(c->size = (int)mdata->NWORDS);
    
    // Need to handle an overflow (result > R = 2^MAXBITS)
    // subtract n from result
    __mmask8 scarry = 0;
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    int i;
    for (i = 0; i < mdata->NWORDS; i++)
    {
        __m512i a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
        __m512i b0 = _mm512_load_epi64(n->data + i * VECLEN_I);
        __m512i a0 = _mm512_mask_sbb_epi52(a1, (__mmask8)m, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
    }
#endif
    return;
}

void veckmul_redc(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    veckmul(a->data, b->data, mdata->mtmp1->data, mdata->mtmp4->data, mdata->NWORDS);

    veckmul_lo(mdata->mtmp1->data, mdata->vnhat->data, mdata->mtmp2->data, mdata->mtmp4->data, mdata->NWORDS);
    veckmul(mdata->mtmp2->data, n->data, mdata->mtmp3->data, mdata->mtmp4->data, mdata->NWORDS);
    mdata->mtmp1->size = (int)(mdata->NWORDS * 2);
    mdata->mtmp3->size = (int)(mdata->NWORDS * 2);
    uint32_t m = vec_bignum52_special_add(mdata->mtmp1->data + mdata->NWORDS * VECLEN_I,
        mdata->mtmp3->data + mdata->NWORDS * VECLEN_I, c->data, (int)mdata->NWORDS);
    
#ifndef USE_AMM
    
    c->WORDS_ALLOC = (uint32_t)(c->size = (int)mdata->NWORDS);
    
    // Need to handle an overflow (result > R = 2^MAXBITS)
    // subtract n from result
    __mmask8 scarry = 0;
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    int i;
    for (i = 0; i < mdata->NWORDS; i++)
    {
        __m512i a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
        __m512i b0 = _mm512_load_epi64(n->data + i * VECLEN_I);
        __m512i a0 = _mm512_mask_sbb_epi52(a1, (__mmask8)m, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
    }

#endif
    return;
}

void vectmul_redc(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    vectmul(a->data, b->data, mdata->mtmp1->data, mdata->mtmp4->data, mdata->NWORDS);

    vectmul(mdata->mtmp1->data, mdata->vnhat->data, mdata->mtmp2->data, mdata->mtmp4->data, mdata->NWORDS);
    vectmul(mdata->mtmp2->data, n->data, mdata->mtmp3->data, mdata->mtmp4->data, mdata->NWORDS);
    mdata->mtmp1->size = (int)(mdata->NWORDS * 2);
    mdata->mtmp3->size = (int)(mdata->NWORDS * 2);
    uint32_t m = vec_bignum52_special_add(mdata->mtmp1->data + mdata->NWORDS * VECLEN_I,
        mdata->mtmp3->data + mdata->NWORDS * VECLEN_I, c->data, (int)mdata->NWORDS);

#ifndef USE_AMM
    c->WORDS_ALLOC = (uint32_t)(c->size = (int)mdata->NWORDS);

    // Need to handle an overflow (result > R = 2^MAXBITS)
    // subtract n from result
    __mmask8 scarry = 0;
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    int i;
    for (i = 0; i < mdata->NWORDS; i++)
    {
        __m512i a1 = _mm512_load_epi64(c->data + i * VECLEN_I);
        __m512i b0 = _mm512_load_epi64(n->data + i * VECLEN_I);
        __m512i a0 = _mm512_mask_sbb_epi52(a1, (__mmask8)m, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN_I, _mm512_and_epi64(vlmask, a0));
    }
#endif
    return;
}

void vectsqr_redc(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    vectmul_redc(a, a, c, n, s, mdata);
    return;
}


#else

void veckmul_redc(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    return;
}

void vecksqr_redc(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    return;
}

#endif
