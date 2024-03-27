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
#include <math.h>

#ifdef _MSC_VER
//#define USE_AVX512F
#endif

#ifdef USE_AVX512F
#include <immintrin.h>

#define USE_AMM 1


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

__m512i __inline _mm512_mask_sbb_src_epi52(__m512i src, __m512i a, __mmask8 m, __mmask8 c, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_mask_sub_epi64(src, m, a, b);
    *cout = _mm512_mask_cmpgt_epu64_mask(m, b, a);
    __m512i t2 = _mm512_mask_sub_epi64(src, m, t, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_kor(*cout, _mm512_mask_cmpgt_epu64_mask(m, t2, t));
    t2 = _mm512_and_epi64(t2, _mm512_set1_epi64(0xfffffffffffffULL));
    return t2;
}

void vecmul52(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata)
{
    int i, j;
    uint32_t NWORDS;
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

    NWORDS = MAX(a->size, b->size);
    NBLOCKS = NWORDS / BLOCKWORDS + ((NWORDS % BLOCKWORDS) > 0);
    NWORDS = NBLOCKS * BLOCKWORDS;

    //printf("base case multiply with blocksize %d\n", NBLOCKS);

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

#ifdef DEBUG_VECMUL
    printf("commencing vecmul52 with NWORDS=%d and NBLOCKS=%d\n", NWORDS, NBLOCKS);
    print_vechexbignum(a, "input a: ");
    print_vechexbignum(b, "input b: ");
    print_vechexbignum(c, "input c: ");
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
            _mm512_store_epi64(c->data + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

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
            _mm512_store_epi64(c->data + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

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
            _mm512_store_epi64(c->data + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

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
            _mm512_store_epi64(c->data + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;
        }
    }

#ifdef DEBUG_VECMUL
    print_vechexbignum(c, "after lo half:");
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

            _mm512_store_epi64(c->data + (i * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(c->data + (i * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(c->data + (i * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(c->data + (i * BLOCKWORDS + 3) * VECLEN, a0);

        }
    }

#ifdef DEBUG_VECMUL
    print_vechexbignum(c, "after hi half:");
#endif

    c->size = 2 * NWORDS;
    return;
}

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

    c->size = NWORDS;

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

void kcombine(uint64_t* c, uint64_t* z1, uint64_t* s1, __mmask8 sm, int words)
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

    // from halfwords to words, we have to copy 'c'
    // to a scratch buffer because we need to add those words
    // back into c later on.
    __m512i zero = _mm512_set1_epi64(0);
    __m512i vone = _mm512_set1_epi64(1);
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i v1, v2, v3, v4, vc = zero;
    __mmask8 carry = 0;

    i = halfwords;
    v1 = _mm512_load_epi64(c + i * VECLEN);
    v2 = _mm512_load_epi64(c + (i - halfwords) * VECLEN);   // lo part
    v3 = _mm512_load_epi64(c + (i + halfwords) * VECLEN);   // hi part
    v4 = _mm512_load_epi64(z1 + (i - halfwords) * VECLEN);  // diff part
    _mm512_store_epi64(s1 + (i - halfwords) * VECLEN, v1);

    v2 = _mm512_add_epi64(v2, v3);
    v4 = _mm512_mask_xor_epi64(v4, sm, v4, vlmask);
    v1 = _mm512_mask_add_epi64(v1, sm, v1, vone);
    v2 = _mm512_add_epi64(v2, v4);
    v1 = _mm512_add_epi64(v1, v2);

    // carryprop
    vc = _mm512_srli_epi64(v1, 52);
    v1 = _mm512_and_epi64(vlmask, v1);

    _mm512_store_epi64(c + i * VECLEN, v1);

    for (i = halfwords + 1; i < words; i++)
    {
        v1 = _mm512_load_epi64(c + i * VECLEN);
        v2 = _mm512_load_epi64(c + (i - halfwords) * VECLEN);
        v3 = _mm512_load_epi64(c + (i + halfwords) * VECLEN);
        v4 = _mm512_load_epi64(z1 + (i - halfwords) * VECLEN);
        // store the output word that we will be overwritting this iteration
        _mm512_store_epi64(s1 + (i - halfwords) * VECLEN, v1);  

        v2 = _mm512_add_epi64(v2, v3);
        v4 = _mm512_mask_xor_epi64(v4, sm, v4, vlmask);

        v4 = _mm512_add_epi64(v4, vc);
        v1 = _mm512_add_epi64(v1, v2);

        v1 = _mm512_add_epi64(v1, v4);

        // carryprop
        vc = _mm512_srli_epi64(v1, 52);
        v1 = _mm512_and_epi64(vlmask, v1);

        _mm512_store_epi64(c + i * VECLEN, v1);
    }

    // stop storing inputs to scratch and instead use
    // the previously stored scratch words as inputs.
    for (i = words; i < words + halfwords; i++)
    {
        v1 = _mm512_load_epi64(c + i * VECLEN);
        v2 = _mm512_load_epi64(s1 + (i - words) * VECLEN);
        v3 = _mm512_load_epi64(c + (i + halfwords) * VECLEN);
        v4 = _mm512_load_epi64(z1 + (i - halfwords) * VECLEN);

        v2 = _mm512_add_epi64(v2, v3);
        v4 = _mm512_mask_xor_epi64(v4, sm, v4, vlmask);

        v4 = _mm512_add_epi64(v4, vc);
        v1 = _mm512_add_epi64(v1, v2);

        v1 = _mm512_add_epi64(v1, v4);

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

    for ( ; i < words; i++)
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


void vecmul52_cios(uint64_t* a, uint64_t* b, uint64_t* c, int n)
{
    int i, j, k;

    __m512d dbias = _mm512_castsi512_pd(set64(0x4670000000000000ULL));
    __m512i vbias1 = set64(0x4670000000000000ULL);
    __m512i zero = set64(0), B, A, C, plo, carry;
    __m512i MASK = set64(DIGIT_MASK);

    for (i = 0; i < 2 * n; i++) {
        _mm512_store_epi64(c + i * VECLEN, _mm512_setzero_si512());
    }

    for (j = 0; j < n; j++) {
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
        _mm512_store_epi64(c + (i + j) * VECLEN, _mm512_and_si512(carry, MASK));
    }

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

    __m512i arestore[3];
    __m512i brestore[3];

    NBLOCKS = words / BLOCKWORDS;
    if ((words % BLOCKWORDS) > 0)
    {
        NBLOCKS++;
        for (i = words, j = 0; i < NBLOCKS * BLOCKWORDS; i++, j++)
        {
            arestore[j] = _mm512_load_epi64(a + i * VECLEN);
            brestore[j] = _mm512_load_epi64(b + i * VECLEN);
            _mm512_store_epi64(a + i * VECLEN, zero);
            _mm512_store_epi64(b + i * VECLEN, zero);
        }
    }

#ifdef DEBUG_MUL
    printf("vecmul52 with blocksize %d, nblocks %d, nwords %d\n",
        BLOCKWORDS, NBLOCKS, BASE_WORDS);
#endif

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

#ifdef DEBUG_VECMUL
    printf("commencing vecmul52 with NWORDS=%d and NBLOCKS=%d\n", NWORDS, NBLOCKS);
    print_vechexbignum(a, "input a: ");
    print_vechexbignum(b, "input b: ");
    print_vechexbignum(c, "input c: ");
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
    print_vechexbignum(c, "after lo half:");
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
    print_vechexbignum(c, "after hi half:");
#endif

    for (i = words, j = 0; i < NBLOCKS * BLOCKWORDS; i++, j++)
    {
        _mm512_store_epi64(a + i * VECLEN, arestore[j]);
        _mm512_store_epi64(b + i * VECLEN, brestore[j]);
    }

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

    NBLOCKS = words / BLOCKWORDS;
    if ((words % BLOCKWORDS) > 0)
    {
        NBLOCKS++;
        for (i = words, j = 0; i < NBLOCKS * BLOCKWORDS; i++, j++)
        {
            arestore[j] = _mm512_load_epi64(a + i * VECLEN);
            _mm512_store_epi64(a + i * VECLEN, zero);
        }
    }

#ifdef DEBUG_MUL
    printf("vecmul52 with blocksize %d, nblocks %d, nwords %d\n",
        BLOCKWORDS, NBLOCKS, BASE_WORDS);
#endif

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

#ifdef DEBUG_VECMUL
    printf("commencing vecmul52 with NWORDS=%d and NBLOCKS=%d\n", NWORDS, NBLOCKS);
    print_vechexbignum(a, "input a: ");
    print_vechexbignum(b, "input b: ");
    print_vechexbignum(c, "input c: ");
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
    print_vechexbignum(c, "after lo half:");
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
    print_vechexbignum(c, "after hi half:");
#endif

    for (i = words, j = 0; i < NBLOCKS * BLOCKWORDS; i++, j++)
    {
        _mm512_store_epi64(a + i * VECLEN, arestore[j]);
    }

    return;
}

void vec_bignum52_mask_sub_n(uint64_t* a, uint64_t* b, uint64_t* c, int words, uint32_t wmask)
{
    // assumptions:
    // a, b, c are of length VECLEN * NWORDS
    // s1 is of length VECLEN
    // a, b, c, and s1 are aligned
    // a and b are both positive
    // a >= b
    int i;
    __mmask8 carry = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;

    if (wmask == 0)
        return;

    // subtract the selected elements ('1' in the mask)
    carry = 0;
    for (i = 0; i < words; i++)
    {
        avec = _mm512_load_epi64(a + i * VECLEN);
        bvec = _mm512_load_epi64(b + i * VECLEN);
        cvec = _mm512_sbb_epi52(avec, carry, bvec, &carry);
        _mm512_mask_store_epi64(c + i * VECLEN, (__mmask8)wmask,
            _mm512_and_epi64(_mm512_set1_epi64(VEC_MAXDIGIT), cvec));
    }

    if (carry)
    {
        // subtract any final borrows that exceed the size of b.
        _mm512_mask_store_epi64(c + i * VECLEN, (__mmask8)wmask & carry, _mm512_set1_epi64(0));
    }

    return;
}

void kmuln(uint64_t* a, uint64_t* b, uint64_t* c, uint64_t* scratch, int words);
void ksqrn(uint64_t* a, uint64_t* c, uint64_t* scratch, int words);

#define veckmul(a, b, c, s, words)				\
  do {									\
    if (words <= 20)						\
      vecmul52_cios(a, b, c, words);      \
    else								\
        kmuln(a, b, c, s, words);				\
  } while (0);

#define vecksqr(a, c, s, words)				\
  do {									\
    if (words <= 20)						\
      vecmul52_cios(a, a, c, words); \
    else								\
        ksqrn(a, c, s, words);				\
  } while (0);



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
    z1 = scratch + (lowords*2)*VECLEN;

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

    //if (hiwords < lowords)
    //{
    //    for (i = hiwords * 2; i < lowords * 2; i++)
    //    {
    //        _mm512_store_si512(c + (2 * lowords + i) * VECLEN, _mm512_setzero_si512());
    //    }
    //}

    kcombine2(c, z1, s1, m1 ^ m2, 2 * lowords, hiwords);
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
    z1 = scratch + (lowords*2)*VECLEN;

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

    //if (hiwords < lowords)
    //{
    //    for (i = hiwords * 2; i < lowords * 2; i++)
    //    {
    //        _mm512_store_si512(c + (2 * lowords + i) * VECLEN, _mm512_setzero_si512());
    //    }
    //}

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
    int i;
    //printf("ksqr len %d words: ", mdata->NWORDS);
    //for (i = mdata->NWORDS - 1; i >= 0; i--)
    //{
    //    printf("%013llx", a->data[i * VECLEN]);
    //}
    //printf("\n");

    vecksqr(a->data, mdata->mtmp1->data, mdata->mtmp4->data, mdata->NWORDS);

    //printf("result: ");
    //for (i = 2 * mdata->NWORDS - 1; i >= 0; i--)
    //{
    //    printf("%013llx", mdata->mtmp1->data[i * VECLEN]);
    //}
    //printf("\n");
    //
    //exit(1);
    veckmul(mdata->mtmp1->data, mdata->vnhat->data, mdata->mtmp2->data, mdata->mtmp4->data, mdata->NWORDS);
    veckmul(mdata->mtmp2->data, n->data, mdata->mtmp3->data, mdata->mtmp4->data, mdata->NWORDS);

    mdata->mtmp1->size = mdata->NWORDS * 2;
    mdata->mtmp3->size = mdata->NWORDS * 2;
    uint32_t m = vec_bignum52_special_add(mdata->mtmp1->data + mdata->NWORDS * VECLEN,
        mdata->mtmp3->data + mdata->NWORDS * VECLEN, c->data, mdata->NWORDS);
    m |= vec_gte52(c, n);
    c->WORDS_ALLOC = c->size = mdata->NWORDS;
    vec_bignum52_mask_sub(c, n, c, m);

    return;
}

void veckmul_redc(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    veckmul(a->data, b->data, mdata->mtmp1->data, mdata->mtmp4->data, mdata->NWORDS);

    veckmul(mdata->mtmp1->data, mdata->vnhat->data, mdata->mtmp2->data, mdata->mtmp4->data, mdata->NWORDS);
    veckmul(mdata->mtmp2->data, n->data, mdata->mtmp3->data, mdata->mtmp4->data, mdata->NWORDS);
    mdata->mtmp1->size = mdata->NWORDS * 2;
    mdata->mtmp3->size = mdata->NWORDS * 2;
    uint32_t m = vec_bignum52_special_add(mdata->mtmp1->data + mdata->NWORDS * VECLEN,
        mdata->mtmp3->data + mdata->NWORDS * VECLEN, c->data, mdata->NWORDS);
    m |= vec_gte52(c, n);
    c->WORDS_ALLOC = c->size = mdata->NWORDS;
    vec_bignum52_mask_sub(c, n, c, m);

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

void vecmulmod52_fixed1040_bfips(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    int i, j, k;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t NBLOCKS = mdata->NBLOCKS;

    // needed in loops
    __m512i i0, i1;
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19
    __m512i c00, c01, c02, c03, c04, c05, c06, c07,
        c08, c09, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19;

    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i nhatvec_e = _mm512_load_epi64(mdata->vrho);
    __m512i hiword = _mm512_set1_epi64(0x000000000000001);
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry_e = 0;
    __mmask8 scarry2;
    __mmask8 scarry;

    // deal with the sign
    c->size = NWORDS;
    c->signmask = a->signmask ^ b->signmask;

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

    // first half mul
    i = 0;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        // finish each triangluar shaped column sum
        a0 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 0) * VECLEN);
        a1 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + 3 * VECLEN);
        b1 = _mm512_load_epi64(b->data + 2 * VECLEN);
        b2 = _mm512_load_epi64(b->data + 1 * VECLEN);
        b3 = _mm512_load_epi64(b->data + 0 * VECLEN);

        VEC_MUL_MUL4_A(a0, a1, a2, b2, b3);
        VEC_MUL_MUL4_B(a0, a1, b0, b1, b2);
        VEC_MUL_MUL2_A(a2, a3, b2, b3);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);
        SUB_BIAS_LO(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);

        // now, column by column, add in the s*n contribution and reduce to 
        // a single 64+x bit accumulator while storing the intermediate product
        // 's' as we go.
        {
            CARRYPROP_ACCUM_0(c00);
            CARRYPROP_ACCUM_1(c01, c00);
            CARRYPROP_ACCUM_2(c02, c00, c01);
            CARRYPROP_ACCUM_3(c03, c00, c01, c02);
        }
    }

    i = 1;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        j = 1;
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

        VEC_MUL_MUL4_A(a0, a1, a2, b2, b3);
        VEC_MUL_MUL4_B(a0, a1, b0, b1, b2);
        VEC_MUL_MUL2_A(a2, a3, b2, b3);

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        j = 0;
        {
            // accumulate s * n
            b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c03, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c02, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c01, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c00, b3, b4, b5, b6);
        }

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);
        SUB_BIAS_LO(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);

        // now, column by column, add in the s*n contribution and reduce to 
        // a single 64+x bit accumulator while storing the intermediate product
        // 's' as we go.
        {
            CARRYPROP_ACCUM_0(c04);
            CARRYPROP_ACCUM_1(c05, c04);
            CARRYPROP_ACCUM_2(c06, c04, c05);
            CARRYPROP_ACCUM_3(c07, c04, c05, c06);
        }
    }

    i = 2;
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

        VEC_MUL_MUL4_A(a0, a1, a2, b2, b3);
        VEC_MUL_MUL4_B(a0, a1, b0, b1, b2);
        VEC_MUL_MUL2_A(a2, a3, b2, b3);

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        //for (j = 0; j < i; j++)
        j = 0;
        b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c03, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c02, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c01, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c00, b3, b4, b5, b6);

        j = 1;
        b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c07, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c06, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c05, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c04, b3, b4, b5, b6);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);
        SUB_BIAS_LO(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);

        // now, column by column, add in the s*n contribution and reduce to 
        // a single 64+x bit accumulator while storing the intermediate product
        // 's' as we go.
        {
            CARRYPROP_ACCUM_0(c08);
            CARRYPROP_ACCUM_1(c09, c08);
            CARRYPROP_ACCUM_2(c10, c08, c09);
            CARRYPROP_ACCUM_3(c11, c08, c09, c10);
        }
    }

    i = 3;
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

        VEC_MUL_MUL4_A(a0, a1, a2, b2, b3);
        VEC_MUL_MUL4_B(a0, a1, b0, b1, b2);
        VEC_MUL_MUL2_A(a2, a3, b2, b3);

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        //for (j = 0; j < i; j++)
        j = 0;
        b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c03, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c02, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c01, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c00, b3, b4, b5, b6);

        j = 1;
        b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c07, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c06, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c05, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c04, b3, b4, b5, b6);

        j = 2;
        b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c11, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c10, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c09, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c08, b3, b4, b5, b6);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);
        SUB_BIAS_LO(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);

        // now, column by column, add in the s*n contribution and reduce to 
        // a single 64+x bit accumulator while storing the intermediate product
        // 's' as we go.
        {
            CARRYPROP_ACCUM_0(c12);
            CARRYPROP_ACCUM_1(c13, c12);
            CARRYPROP_ACCUM_2(c14, c12, c13);
            CARRYPROP_ACCUM_3(c15, c12, c13, c14);
        }
    }

    i = 4;
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

        VEC_MUL_MUL4_A(a0, a1, a2, b2, b3);
        VEC_MUL_MUL4_B(a0, a1, b0, b1, b2);
        VEC_MUL_MUL2_A(a2, a3, b2, b3);

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        //for (j = 0; j < i; j++)
        j = 0;
        b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c03, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c02, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c01, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c00, b3, b4, b5, b6);

        j = 1;
        b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c07, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c06, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c05, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c04, b3, b4, b5, b6);

        j = 2;
        b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c11, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c10, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c09, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c08, b3, b4, b5, b6);

        j = 3;
        b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c15, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c14, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c13, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c12, b3, b4, b5, b6);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);
        SUB_BIAS_LO(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);

        // now, column by column, add in the s*n contribution and reduce to 
        // a single 64+x bit accumulator while storing the intermediate product
        // 's' as we go.
        {
            CARRYPROP_ACCUM_0(c16);
            CARRYPROP_ACCUM_1(c17, c16);
            CARRYPROP_ACCUM_2(c18, c16, c17);
            CARRYPROP_ACCUM_3(c19, c16, c17, c18);
        }
    }

    // second half mul
    //for (i = NBLOCKS; i < 2 * NBLOCKS; i++)
    i = 5;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        //for (j = i - NBLOCKS + 1; j < NBLOCKS; j++)
        j = 1;
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

        // accumulate s * n
        b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c19, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c18, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c17, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c16, b3, b4, b5, b6);

        j = 2;
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

        // accumulate s * n
        b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c15, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c14, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c13, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c12, b3, b4, b5, b6);

        j = 3;
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

        // accumulate s * n
        b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c11, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c10, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c09, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c08, b3, b4, b5, b6);

        j = 4;
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

        // accumulate s * n
        b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c07, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c06, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c05, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c04, b3, b4, b5, b6);


        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN);

        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        a1 = c01;
        a2 = c02;
        a3 = c03;

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);
        SUB_BIAS_LO(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);

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
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN, a0);

        }
    }

    i = 6;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        j = 2;
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

        // accumulate s * n
        b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c19, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c18, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c17, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c16, b3, b4, b5, b6);

        j = 3;
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

        // accumulate s * n
        b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c15, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c14, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c13, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c12, b3, b4, b5, b6);

        j = 4;
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

        // accumulate s * n
        b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c11, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c10, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c09, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c08, b3, b4, b5, b6);



        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN);

        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // i = 3, a1..3 = c[1..3]
        // i = 4, a1..3 = c[5..7]
        // i = 5, a1..3 = c[9..11]
        a1 = c05;
        a2 = c06;
        a3 = c07;

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);
        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);
        SUB_BIAS_LO(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);

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
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN, a0);

        }
    }

    i = 7;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        //for (j = i - NBLOCKS + 1; j < NBLOCKS; j++)
        j = 3;
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

        // accumulate s * n
        b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c19, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c18, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c17, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c16, b3, b4, b5, b6);

        j = 4;
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

        // accumulate s * n
        b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c15, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c14, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c13, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c12, b3, b4, b5, b6);

        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN);

        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // i = 3, a1..3 = c[1..3]
        // i = 4, a1..3 = c[5..7]
        // i = 5, a1..3 = c[9..11]
        a1 = c09;
        a2 = c10;
        a3 = c11;

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);
        SUB_BIAS_LO(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);

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
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN, a0);

        }
    }

    i = 8;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        j = 4;
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

        // accumulate s * n
        b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c19, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c18, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c17, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c16, b3, b4, b5, b6);

        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN);

        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // i = 3, a1..3 = c[1..3]
        // i = 4, a1..3 = c[5..7]
        // i = 5, a1..3 = c[9..11]
        a1 = c13;
        a2 = c14;
        a3 = c15;

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);
        SUB_BIAS_LO(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);

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
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN, a0);

        }
    }

    i = 9;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN);

        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // i = 3, a1..3 = c[1..3]
        // i = 4, a1..3 = c[5..7]
        // i = 5, a1..3 = c[9..11]
        a1 = c17;
        a2 = c18;
        a3 = c19;

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);
        SUB_BIAS_LO(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);

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
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN, a0);

        }
    }
    _mm512_store_epi64(c->data + NWORDS * VECLEN, zero);

    c->size = NWORDS;
    return;
}

void vecmulmod52_fixed832_bfips(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    int i, j, k;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t NBLOCKS = mdata->NBLOCKS;

    // needed in loops
    __m512i i0, i1;
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
    __m512i nhatvec_e = _mm512_load_epi64(mdata->vrho);
    __m512i hiword = _mm512_set1_epi64(0x000000000000001);
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry_e = 0;
    __mmask8 scarry2;
    __mmask8 scarry;

    // deal with the sign
    c->size = NWORDS;
    c->signmask = a->signmask ^ b->signmask;

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

    // first half mul
    i = 0;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        // finish each triangluar shaped column sum
        a0 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 0) * VECLEN);
        a1 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + 3 * VECLEN);
        b1 = _mm512_load_epi64(b->data + 2 * VECLEN);
        b2 = _mm512_load_epi64(b->data + 1 * VECLEN);
        b3 = _mm512_load_epi64(b->data + 0 * VECLEN);

        VEC_MUL_MUL4_A(a0, a1, a2, b2, b3);
        VEC_MUL_MUL4_B(a0, a1, b0, b1, b2);
        VEC_MUL_MUL2_A(a2, a3, b2, b3);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);
        SUB_BIAS_LO(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);

        // now, column by column, add in the s*n contribution and reduce to 
        // a single 64+x bit accumulator while storing the intermediate product
        // 's' as we go.
        {
            CARRYPROP_ACCUM_0(c00);
            CARRYPROP_ACCUM_1(c01, c00);
            CARRYPROP_ACCUM_2(c02, c00, c01);
            CARRYPROP_ACCUM_3(c03, c00, c01, c02);
        }
    }

    i = 1;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        j = 1;
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

        VEC_MUL_MUL4_A(a0, a1, a2, b2, b3);
        VEC_MUL_MUL4_B(a0, a1, b0, b1, b2);
        VEC_MUL_MUL2_A(a2, a3, b2, b3);

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        j = 0;
        {
            // accumulate s * n
            b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c03, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c02, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c01, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c00, b3, b4, b5, b6);
        }

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);
        SUB_BIAS_LO(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);

        // now, column by column, add in the s*n contribution and reduce to 
        // a single 64+x bit accumulator while storing the intermediate product
        // 's' as we go.
        {
            CARRYPROP_ACCUM_0(c04);
            CARRYPROP_ACCUM_1(c05, c04);
            CARRYPROP_ACCUM_2(c06, c04, c05);
            CARRYPROP_ACCUM_3(c07, c04, c05, c06);
        }
    }

    i = 2;
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

        VEC_MUL_MUL4_A(a0, a1, a2, b2, b3);
        VEC_MUL_MUL4_B(a0, a1, b0, b1, b2);
        VEC_MUL_MUL2_A(a2, a3, b2, b3);

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        //for (j = 0; j < i; j++)
        j = 0;
        b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c03, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c02, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c01, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c00, b3, b4, b5, b6);

        j = 1;
        b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c07, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c06, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c05, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c04, b3, b4, b5, b6);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);
        SUB_BIAS_LO(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);

        // now, column by column, add in the s*n contribution and reduce to 
        // a single 64+x bit accumulator while storing the intermediate product
        // 's' as we go.
        {
            CARRYPROP_ACCUM_0(c08);
            CARRYPROP_ACCUM_1(c09, c08);
            CARRYPROP_ACCUM_2(c10, c08, c09);
            CARRYPROP_ACCUM_3(c11, c08, c09, c10);
        }
    }

    i = 3;
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

        VEC_MUL_MUL4_A(a0, a1, a2, b2, b3);
        VEC_MUL_MUL4_B(a0, a1, b0, b1, b2);
        VEC_MUL_MUL2_A(a2, a3, b2, b3);

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        //for (j = 0; j < i; j++)
        j = 0;
        b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c03, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c02, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c01, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c00, b3, b4, b5, b6);

        j = 1;
        b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c07, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c06, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c05, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c04, b3, b4, b5, b6);

        j = 2;
        b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c11, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c10, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c09, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c08, b3, b4, b5, b6);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);
        SUB_BIAS_LO(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);

        // now, column by column, add in the s*n contribution and reduce to 
        // a single 64+x bit accumulator while storing the intermediate product
        // 's' as we go.
        {
            CARRYPROP_ACCUM_0(c12);
            CARRYPROP_ACCUM_1(c13, c12);
            CARRYPROP_ACCUM_2(c14, c12, c13);
            CARRYPROP_ACCUM_3(c15, c12, c13, c14);
        }
    }

    // second half mul
    //for (i = NBLOCKS; i < 2 * NBLOCKS; i++)
    i = 4;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        //for (j = i - NBLOCKS + 1; j < NBLOCKS; j++)
        j = 1;
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

        // accumulate s * n
        b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c15, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c14, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c13, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c12, b3, b4, b5, b6);

        j = 2;
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

        // accumulate s * n
        b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c11, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c10, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c09, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c08, b3, b4, b5, b6);

        j = 3;
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

        // accumulate s * n
        b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c07, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c06, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c05, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c04, b3, b4, b5, b6);


        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN);

        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        a1 = c01;
        a2 = c02;
        a3 = c03;

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);
        SUB_BIAS_LO(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);

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
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN, a0);

        }
    }

    i = 5;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        j = 2;
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

        // accumulate s * n
        b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c15, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c14, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c13, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c12, b3, b4, b5, b6);

        j = 3;
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

        // accumulate s * n
        b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c11, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c10, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c09, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c08, b3, b4, b5, b6);



        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN);

        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // i = 3, a1..3 = c[1..3]
        // i = 4, a1..3 = c[5..7]
        // i = 5, a1..3 = c[9..11]
        a1 = c05;
        a2 = c06;
        a3 = c07;

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);
        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);
        SUB_BIAS_LO(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);

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
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN, a0);

        }
    }

    i = 6;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        //for (j = i - NBLOCKS + 1; j < NBLOCKS; j++)
        j = 3;
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

        // accumulate s * n
        b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c15, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c14, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c13, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c12, b3, b4, b5, b6);

        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN);

        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // i = 3, a1..3 = c[1..3]
        // i = 4, a1..3 = c[5..7]
        // i = 5, a1..3 = c[9..11]
        a1 = c09;
        a2 = c10;
        a3 = c11;

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);
        SUB_BIAS_LO(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);

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
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN, a0);

        }
    }

    i = 7;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN);

        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // i = 3, a1..3 = c[1..3]
        // i = 4, a1..3 = c[5..7]
        // i = 5, a1..3 = c[9..11]
        a1 = c13;
        a2 = c14;
        a3 = c15;

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);
        SUB_BIAS_LO(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);

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
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN, a0);

        }
    }

    _mm512_store_epi64(c->data + NWORDS * VECLEN, zero);

    c->size = NWORDS;
    return;
}

void vecmulmod52_fixed624_bfips(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    int i, j, k;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t NBLOCKS = mdata->NBLOCKS;

    // needed in loops
    __m512i i0, i1;
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19
    __m512i c00, c01, c02, c03, c04, c05, c06, c07, c08, c09, c10, c11;

    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i nhatvec_e = _mm512_load_epi64(mdata->vrho);
    __m512i hiword = _mm512_set1_epi64(0x000000000000001);
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry_e = 0;
    __mmask8 scarry2;
    __mmask8 scarry;

    // deal with the sign
    c->size = NWORDS;
    c->signmask = a->signmask ^ b->signmask;

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;
    
    // first half mul
    i = 0;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        // finish each triangluar shaped column sum
        a0 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 0) * VECLEN);
        a1 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + 3 * VECLEN);
        b1 = _mm512_load_epi64(b->data + 2 * VECLEN);
        b2 = _mm512_load_epi64(b->data + 1 * VECLEN);
        b3 = _mm512_load_epi64(b->data + 0 * VECLEN);

        VEC_MUL_MUL4_A(a0, a1, a2, b2, b3);
        VEC_MUL_MUL4_B(a0, a1, b0, b1, b2);
        VEC_MUL_MUL2_A(a2, a3, b2, b3);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);
        SUB_BIAS_LO(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);

        // now, column by column, add in the s*n contribution and reduce to 
        // a single 64+x bit accumulator while storing the intermediate product
        // 's' as we go.
        {
            CARRYPROP_ACCUM_0(c00);
            CARRYPROP_ACCUM_1(c01, c00);
            CARRYPROP_ACCUM_2(c02, c00, c01);
            CARRYPROP_ACCUM_3(c03, c00, c01, c02);
        }
    }

    i = 1;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        j = 1;
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

        VEC_MUL_MUL4_A(a0, a1, a2, b2, b3);
        VEC_MUL_MUL4_B(a0, a1, b0, b1, b2);
        VEC_MUL_MUL2_A(a2, a3, b2, b3);

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        j = 0;
        {
            // accumulate s * n
            b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c03, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c02, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c01, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c00, b3, b4, b5, b6);
        }

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);
        SUB_BIAS_LO(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);

        // now, column by column, add in the s*n contribution and reduce to 
        // a single 64+x bit accumulator while storing the intermediate product
        // 's' as we go.
        {
            CARRYPROP_ACCUM_0(c04);
            CARRYPROP_ACCUM_1(c05, c04);
            CARRYPROP_ACCUM_2(c06, c04, c05);
            CARRYPROP_ACCUM_3(c07, c04, c05, c06);
        }
    }

    i = 2;
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

        VEC_MUL_MUL4_A(a0, a1, a2, b2, b3);
        VEC_MUL_MUL4_B(a0, a1, b0, b1, b2);
        VEC_MUL_MUL2_A(a2, a3, b2, b3);

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        j = 0;
        b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c03, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c02, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c01, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c00, b3, b4, b5, b6);

        j = 1;
        b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c07, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c06, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c05, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c04, b3, b4, b5, b6);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);
        SUB_BIAS_LO(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);

        // now, column by column, add in the s*n contribution and reduce to 
        // a single 64+x bit accumulator while storing the intermediate product
        // 's' as we go.
        {
            CARRYPROP_ACCUM_0(c08);
            CARRYPROP_ACCUM_1(c09, c08);
            CARRYPROP_ACCUM_2(c10, c08, c09);
            CARRYPROP_ACCUM_3(c11, c08, c09, c10);
        }
    }

    // second half mul
    i = 3;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        j = 1;
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

        // accumulate s * n
        b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c11, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c10, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c09, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c08, b3, b4, b5, b6);

        j = 2;
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

        // accumulate s * n
        b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c07, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c06, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c05, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c04, b3, b4, b5, b6);


        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN);

        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        a1 = c01;
        a2 = c02;
        a3 = c03;

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);
        SUB_BIAS_LO(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);

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
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN, a0);

        }
    }

    i = 4;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        j = 2;
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

        // accumulate s * n
        b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c11, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c10, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c09, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c08, b3, b4, b5, b6);



        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN);

        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // i = 3, a1..3 = c[1..3]
        // i = 4, a1..3 = c[5..7]
        // i = 5, a1..3 = c[9..11]
        a1 = c05;
        a2 = c06;
        a3 = c07;

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);
        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);
        SUB_BIAS_LO(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);

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
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN, a0);

        }
    }

    i = 5;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN);

        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // i = 3, a1..3 = c[1..3]
        // i = 4, a1..3 = c[5..7]
        // i = 5, a1..3 = c[9..11]
        a1 = c09;
        a2 = c10;
        a3 = c11;

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);
        SUB_BIAS_LO(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);

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
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN, a0);

        }
    }

    _mm512_store_epi64(c->data + NWORDS * VECLEN, zero);

    c->size = NWORDS;
    return;
}

void vecmulmod52_fixed416_bfips(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    int i, j, k;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t NBLOCKS = mdata->NBLOCKS;

    // needed in loops
    __m512i i0, i1;
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19
    __m512i c00, c01, c02, c03, c04, c05, c06, c07;

    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i nhatvec_e = _mm512_load_epi64(mdata->vrho);
    __m512i hiword = _mm512_set1_epi64(0x000000000000001);
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry_e = 0;
    __mmask8 scarry2;
    __mmask8 scarry;

    // deal with the sign
    c->size = NWORDS;
    c->signmask = a->signmask ^ b->signmask;

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

    // first half mul
    i = 0;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        //for (j = i; j > 0; j--)

        // finish each triangluar shaped column sum
        a0 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 0) * VECLEN);
        a1 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + 3 * VECLEN);
        b1 = _mm512_load_epi64(b->data + 2 * VECLEN);
        b2 = _mm512_load_epi64(b->data + 1 * VECLEN);
        b3 = _mm512_load_epi64(b->data + 0 * VECLEN);

        VEC_MUL_MUL4_A(a0, a1, a2, b2, b3);
        VEC_MUL_MUL4_B(a0, a1, b0, b1, b2);
        VEC_MUL_MUL2_A(a2, a3, b2, b3);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);
        SUB_BIAS_LO(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);

        // now, column by column, add in the s*n contribution and reduce to 
        // a single 64+x bit accumulator while storing the intermediate product
        // 's' as we go.
        {
            CARRYPROP_ACCUM_0(c00);
            CARRYPROP_ACCUM_1(c01, c00);
            CARRYPROP_ACCUM_2(c02, c00, c01);
            CARRYPROP_ACCUM_3(c03, c00, c01, c02);
        }
    }

    i = 1;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        //for (j = i; j > 0; j--)
        j = 1;
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

        VEC_MUL_MUL4_A(a0, a1, a2, b2, b3);
        VEC_MUL_MUL4_B(a0, a1, b0, b1, b2);
        VEC_MUL_MUL2_A(a2, a3, b2, b3);

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        //for (j = 0; j < i; j++)
        j = 0;
        {
            // accumulate s * n
            b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c03, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c02, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c01, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c00, b3, b4, b5, b6);
        }

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);
        SUB_BIAS_LO(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);

        // now, column by column, add in the s*n contribution and reduce to 
        // a single 64+x bit accumulator while storing the intermediate product
        // 's' as we go.
        {
            CARRYPROP_ACCUM_0(c04);
            CARRYPROP_ACCUM_1(c05, c04);
            CARRYPROP_ACCUM_2(c06, c04, c05);
            CARRYPROP_ACCUM_3(c07, c04, c05, c06);
        }
    }

    // second half mul
    i = 2;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        //for (j = i - NBLOCKS + 1; j < NBLOCKS; j++)
        j = 1;
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

        // accumulate s * n
        b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c07, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c06, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c05, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c04, b3, b4, b5, b6);

        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN);

        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        a1 = c01;
        a2 = c02;
        a3 = c03;

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);
        SUB_BIAS_LO(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);

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
            //_mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 0) * VECLEN, a3);
            //_mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN, a2);
            //_mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN, a1);
            //_mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN, a0);

            c00 = a3;
            c01 = a2;
            c02 = a1;
            c03 = a0;

        }
    }

    i = 3;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        //for (j = i - NBLOCKS + 1; j < NBLOCKS; j++)

        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN);

        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // i = 3, a1..3 = c[1..3]
        // i = 4, a1..3 = c[5..7]
        // i = 5, a1..3 = c[9..11]
        a1 = c05;
        a2 = c06;
        a3 = c07;

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);
        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);
        SUB_BIAS_LO(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);

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
            c04 = a3;
            c05 = a2;
            c06 = a1;
            c07 = a0;


            //_mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 0) * VECLEN, a3);
            //_mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN, a2);
            //_mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN, a1);
            //_mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN, a0);

        }
    }

#ifdef USE_AMM
    _mm512_store_epi64(c->data + 0 * VECLEN, c00);
    _mm512_store_epi64(c->data + 1 * VECLEN, c01);
    _mm512_store_epi64(c->data + 2 * VECLEN, c02);
    _mm512_store_epi64(c->data + 3 * VECLEN, c03);
    _mm512_store_epi64(c->data + 4 * VECLEN, c04);
    _mm512_store_epi64(c->data + 5 * VECLEN, c05);
    _mm512_store_epi64(c->data + 6 * VECLEN, c06);
    _mm512_store_epi64(c->data + 7 * VECLEN, c07);
    _mm512_store_epi64(c->data + NWORDS * VECLEN, zero);
#else

    __m512i cvec;
    __m512i nvec;
    __m512i bvec;

    // compare
    scarry = 0;	    // sub mask
    scarry2 = 0;	// keep looking mask

    cvec = c07;
    nvec = _mm512_load_epi64(mdata->n->data + 7 * VECLEN);
    // compare those that have not already been decided using the mask
    scarry |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_GT);
    scarry2 |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_LT);

    // decided all of them, stop comparing.
    if (scarry2 == 0xff) goto sub;

    cvec = c06;
    nvec = _mm512_load_epi64(mdata->n->data + 6 * VECLEN);
    // compare those that have not already been decided using the mask
    scarry |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_GT);
    scarry2 |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_LT);

    // decided all of them, stop comparing.
    if (scarry2 == 0xff) goto sub;

    cvec = c05;
    nvec = _mm512_load_epi64(mdata->n->data + 5 * VECLEN);
    // compare those that have not already been decided using the mask
    scarry |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_GT);
    scarry2 |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_LT);

    // decided all of them, stop comparing.
    if (scarry2 == 0xff) goto sub;

    cvec = c04;
    nvec = _mm512_load_epi64(mdata->n->data + 4 * VECLEN);
    // compare those that have not already been decided using the mask
    scarry |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_GT);
    scarry2 |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_LT);

    // decided all of them, stop comparing.
    if (scarry2 == 0xff) goto sub;

    cvec = c03;
    nvec = _mm512_load_epi64(mdata->n->data + 3 * VECLEN);
    // compare those that have not already been decided using the mask
    scarry |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_GT);
    scarry2 |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_LT);

    // decided all of them, stop comparing.
    if (scarry2 == 0xff) goto sub;

    cvec = c02;
    nvec = _mm512_load_epi64(mdata->n->data + 2 * VECLEN);
    // compare those that have not already been decided using the mask
    scarry |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_GT);
    scarry2 |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_LT);

    // decided all of them, stop comparing.
    if (scarry2 == 0xff) goto sub;

    cvec = c01;
    nvec = _mm512_load_epi64(mdata->n->data + 1 * VECLEN);
    // compare those that have not already been decided using the mask
    scarry |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_GT);
    scarry2 |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_LT);

    // decided all of them, stop comparing.
    if (scarry2 == 0xff) goto sub;

    cvec = c00;
    nvec = _mm512_load_epi64(mdata->n->data + 0 * VECLEN);
    // compare those that have not already been decided using the mask
    scarry |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_GT);
    scarry2 |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_LT);

    // decided all of them, stop comparing.
    if (scarry2 == 0xff) goto sub;

    // check for equal as well by flipping mask bits that have still
    // not been decided (i.e., are equal)
    scarry |= (~scarry2);

sub:

    if (scarry == 0) goto done;
    
    // subtract n from c when c is not less than n, as indicated by a 1 bit in mask
    scarry2 = 0;

    nvec = _mm512_load_epi64(mdata->n->data + 0 * VECLEN);
    bvec = _mm512_mask_sbb_src_epi52(zero, c00, scarry, scarry2, nvec, &scarry2);
    _mm512_store_epi64(c->data + 0 * VECLEN, bvec);

    nvec = _mm512_load_epi64(mdata->n->data + 1 * VECLEN);
    bvec = _mm512_mask_sbb_src_epi52(zero, c01, scarry, scarry2, nvec, &scarry2);
    _mm512_store_epi64(c->data + 1 * VECLEN, bvec);

    nvec = _mm512_load_epi64(mdata->n->data + 2 * VECLEN);
    bvec = _mm512_mask_sbb_src_epi52(zero, c02, scarry, scarry2, nvec, &scarry2);
    _mm512_store_epi64(c->data + 2 * VECLEN, bvec);

    nvec = _mm512_load_epi64(mdata->n->data + 3 * VECLEN);
    bvec = _mm512_mask_sbb_src_epi52(zero, c03, scarry, scarry2, nvec, &scarry2);
    _mm512_store_epi64(c->data + 3 * VECLEN, bvec);

    nvec = _mm512_load_epi64(mdata->n->data + 4 * VECLEN);
    bvec = _mm512_mask_sbb_src_epi52(zero, c04, scarry, scarry2, nvec, &scarry2);
    _mm512_store_epi64(c->data + 4 * VECLEN, bvec);

    nvec = _mm512_load_epi64(mdata->n->data + 5 * VECLEN);
    bvec = _mm512_mask_sbb_src_epi52(zero, c05, scarry, scarry2, nvec, &scarry2);
    _mm512_store_epi64(c->data + 5 * VECLEN, bvec);

    nvec = _mm512_load_epi64(mdata->n->data + 6 * VECLEN);
    bvec = _mm512_mask_sbb_src_epi52(zero, c06, scarry, scarry2, nvec, &scarry2);
    _mm512_store_epi64(c->data + 6 * VECLEN, bvec);

    nvec = _mm512_load_epi64(mdata->n->data + 7 * VECLEN);
    bvec = _mm512_mask_sbb_src_epi52(zero, c07, scarry, scarry2, nvec, &scarry2);
    _mm512_store_epi64(c->data + 7 * VECLEN, bvec);

done:

#endif
    

    c->size = NWORDS;
    return;
}

#ifdef IFMA
void vecmulmod52_fixed624_bfips_ifma(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    int i, j, k;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t NBLOCKS = mdata->NBLOCKS;

    // needed in loops
    __m512i i0, i1;
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19
    __m512i c00, c01, c02, c03, c04, c05, c06, c07, c08, c09, c10, c11;

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i nhatvec_e = _mm512_load_epi64(mdata->vrho);
    __m512i hiword = _mm512_set1_epi64(0x000000000000001);
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry_e = 0;
    __mmask8 scarry2;
    __mmask8 scarry;

    // deal with the sign
    c->size = NWORDS;
    c->signmask = a->signmask ^ b->signmask;

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

    // first half mul
    i = 0;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        // finish each triangluar shaped column sum
        a0 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 0) * VECLEN);
        a1 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + 3 * VECLEN);
        b1 = _mm512_load_epi64(b->data + 2 * VECLEN);
        b2 = _mm512_load_epi64(b->data + 1 * VECLEN);
        b3 = _mm512_load_epi64(b->data + 0 * VECLEN);

        VEC_MUL_ACCUM_LOHI_PD(a0, b3, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a1, b3, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a0, b2, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a2, b3, te4, te5);

        VEC_MUL_ACCUM_LOHI_PD(a1, b2, te4, te5);
        VEC_MUL_ACCUM_LOHI_PD(a0, b1, te4, te5);
        VEC_MUL_ACCUM_LOHI_PD(a1, b1, te6, te7);
        VEC_MUL_ACCUM_LOHI_PD(a0, b0, te6, te7);

        VEC_MUL_ACCUM_LOHI_PD(a3, b3, te6, te7);
        VEC_MUL_ACCUM_LOHI_PD(a2, b2, te6, te7);

        // now, column by column, add in the s*n contribution and reduce to 
        // a single 64+x bit accumulator while storing the intermediate product
        // 's' as we go.
        {
            j = 0;
            // accumulate this column-sum
            // now carry propagate low to high.
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);

            //_mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            c00 = a0;
            b0 = _mm512_load_epi64(n->data + 0 * VECLEN);

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            b1 = _mm512_load_epi64(n->data + 1 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c00, b1, acc_e0, acc_e1);

            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            //_mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            c01 = a0;

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            b1 = _mm512_load_epi64(n->data + 2 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c00, b1, acc_e0, acc_e1);
            b1 = _mm512_load_epi64(n->data + 1 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c01, b1, acc_e0, acc_e1);

            // these are bottleneck high-latency sequentially dependent instructions...
            // not sure what can be done about it.
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            //_mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            c02 = a0;

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            b1 = _mm512_load_epi64(n->data + 3 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c00, b1, acc_e0, acc_e1);
            b1 = _mm512_load_epi64(n->data + 2 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c01, b1, acc_e0, acc_e1);
            b1 = _mm512_load_epi64(n->data + 1 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c02, b1, acc_e0, acc_e1);

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            //_mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            c03 = a0;

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;
        }
    }

    i = 1;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        j = 1;
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

        VEC_MUL_ACCUM_LOHI_PD(a0, b3, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a1, b3, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a0, b2, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a2, b3, te4, te5);

        VEC_MUL_ACCUM_LOHI_PD(a1, b2, te4, te5);
        VEC_MUL_ACCUM_LOHI_PD(a0, b1, te4, te5);
        VEC_MUL_ACCUM_LOHI_PD(a1, b1, te6, te7);
        VEC_MUL_ACCUM_LOHI_PD(a0, b0, te6, te7);

        VEC_MUL_ACCUM_LOHI_PD(a3, b3, te6, te7);
        VEC_MUL_ACCUM_LOHI_PD(a2, b2, te6, te7);


        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        j = 0;
        {
            // accumulate s * n
            b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c03, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c02, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c01, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c00, b3, b4, b5, b6);
        }


        // now, column by column, add in the s*n contribution and reduce to 
        // a single 64+x bit accumulator while storing the intermediate product
        // 's' as we go.
        {
            j = 0;
            // accumulate this column-sum
            // now carry propagate low to high.
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);

            //_mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            c04 = a0;
            b0 = _mm512_load_epi64(n->data + 0 * VECLEN);

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            b1 = _mm512_load_epi64(n->data + 1 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c04, b1, acc_e0, acc_e1);

            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            //_mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            c05 = a0;

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            b1 = _mm512_load_epi64(n->data + 2 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c04, b1, acc_e0, acc_e1);
            b1 = _mm512_load_epi64(n->data + 1 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c05, b1, acc_e0, acc_e1);

            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            //_mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            c06 = a0;

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            b1 = _mm512_load_epi64(n->data + 3 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c04, b1, acc_e0, acc_e1);
            b1 = _mm512_load_epi64(n->data + 2 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c05, b1, acc_e0, acc_e1);
            b1 = _mm512_load_epi64(n->data + 1 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c06, b1, acc_e0, acc_e1);

            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            //_mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            c07 = a0;

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;
        }
    }

    i = 2;
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

        VEC_MUL_ACCUM_LOHI_PD(a0, b3, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a1, b3, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a0, b2, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a2, b3, te4, te5);


        VEC_MUL_ACCUM_LOHI_PD(a1, b2, te4, te5);
        VEC_MUL_ACCUM_LOHI_PD(a0, b1, te4, te5);
        VEC_MUL_ACCUM_LOHI_PD(a1, b1, te6, te7);
        VEC_MUL_ACCUM_LOHI_PD(a0, b0, te6, te7);


        VEC_MUL_ACCUM_LOHI_PD(a3, b3, te6, te7);
        VEC_MUL_ACCUM_LOHI_PD(a2, b2, te6, te7);


        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        j = 0;
        b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c03, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c02, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c01, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c00, b3, b4, b5, b6);

        j = 1;
        b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c07, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c06, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c05, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c04, b3, b4, b5, b6);


        // now, column by column, add in the s*n contribution and reduce to 
        // a single 64+x bit accumulator while storing the intermediate product
        // 's' as we go.
        {
            j = 0;
            // accumulate this column-sum
            // now carry propagate low to high.
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);

            //_mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            c08 = a0;
            b0 = _mm512_load_epi64(n->data + 0 * VECLEN);

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            b1 = _mm512_load_epi64(n->data + 1 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c08, b1, acc_e0, acc_e1);

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            //_mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            c09 = a0;

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            b1 = _mm512_load_epi64(n->data + 2 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c08, b1, acc_e0, acc_e1);
            b1 = _mm512_load_epi64(n->data + 1 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c09, b1, acc_e0, acc_e1);

            // these are bottleneck high-latency sequentially dependent instructions...
            // not sure what can be done about it.
            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            //_mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            c10 = a0;

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            b1 = _mm512_load_epi64(n->data + 3 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c08, b1, acc_e0, acc_e1);
            b1 = _mm512_load_epi64(n->data + 2 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c09, b1, acc_e0, acc_e1);
            b1 = _mm512_load_epi64(n->data + 1 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c10, b1, acc_e0, acc_e1);

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            //_mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            c11 = a0;

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;
        }
    }

    // second half mul
    i = 3;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        j = 1;
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

        // accumulate s * n
        b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c11, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c10, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c09, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c08, b3, b4, b5, b6);

        j = 2;
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

        // accumulate s * n
        b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c07, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c06, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c05, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c04, b3, b4, b5, b6);


        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN);

        VEC_MUL_ACCUM_LOHI_PD(a1, b0, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b1, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b0, te2, te3);


        VEC_MUL_ACCUM_LOHI_PD(a3, b2, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a3, b1, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a3, b0, te4, te5);


        a1 = c01;
        a2 = c02;
        a3 = c03;

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        // a1*b0 -> t0/1
        // a2*b1 -> t0/1
        // a2*b0 -> t2/3
        VEC_MUL_ACCUM_LOHI_PD(a1, b0, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b1, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b0, te2, te3);


        VEC_MUL_ACCUM_LOHI_PD(a3, b2, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a3, b1, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a3, b0, te4, te5);


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

            // i = 3, a3..0 = c[0..3]
            // i = 4, a3..0 = c[4..7]
            // i = 5, a3..0 = c[8..11]
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN, a0);

        }
    }

    i = 4;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        j = 2;
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

        // accumulate s * n
        b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
        b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
        b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
        b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
        b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

        VEC_MUL4_ACCUM(c11, b0, b1, b2, b3);
        VEC_MUL4_ACCUM(c10, b1, b2, b3, b4);
        VEC_MUL4_ACCUM(c09, b2, b3, b4, b5);
        VEC_MUL4_ACCUM(c08, b3, b4, b5, b6);



        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN);

        VEC_MUL_ACCUM_LOHI_PD(a1, b0, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b1, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b0, te2, te3);



        VEC_MUL_ACCUM_LOHI_PD(a3, b2, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a3, b1, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a3, b0, te4, te5);

        // i = 3, a1..3 = c[1..3]
        // i = 4, a1..3 = c[5..7]
        // i = 5, a1..3 = c[9..11]
        a1 = c05;
        a2 = c06;
        a3 = c07;

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        // a1*b0 -> t0/1
        // a2*b1 -> t0/1
        // a2*b0 -> t2/3
        VEC_MUL_ACCUM_LOHI_PD(a1, b0, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b1, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b0, te2, te3);


        VEC_MUL_ACCUM_LOHI_PD(a3, b2, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a3, b1, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a3, b0, te4, te5);


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

            // i = 3, a3..0 = c[0..3]
            // i = 4, a3..0 = c[4..7]
            // i = 5, a3..0 = c[8..11]
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN, a0);

        }
    }

    i = 5;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN);

        VEC_MUL_ACCUM_LOHI_PD(a1, b0, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b1, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b0, te2, te3);


        VEC_MUL_ACCUM_LOHI_PD(a3, b2, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a3, b1, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a3, b0, te4, te5);


        // i = 3, a1..3 = c[1..3]
        // i = 4, a1..3 = c[5..7]
        // i = 5, a1..3 = c[9..11]
        a1 = c09;
        a2 = c10;
        a3 = c11;

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        // a1*b0 -> t0/1
        // a2*b1 -> t0/1
        // a2*b0 -> t2/3
        VEC_MUL_ACCUM_LOHI_PD(a1, b0, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b1, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b0, te2, te3);


        VEC_MUL_ACCUM_LOHI_PD(a3, b2, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a3, b1, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a3, b0, te4, te5);

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

            // i = 3, a3..0 = c[0..3]
            // i = 4, a3..0 = c[4..7]
            // i = 5, a3..0 = c[8..11]
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(c->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN, a0);

        }
    }

    _mm512_store_epi64(c->data + NWORDS * VECLEN, zero);

    c->size = NWORDS;
    return;
}
void vecsqrmod52_fixed624_bfips_ifma(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_monty_t* mdata)
{
    // 8x sqr:
    // input 8 bignums in the even lanes of a.
    // output 8 squaremod bignums in the even lanes of c.
    int i, j, k;
    vec_bignum_t* b = a;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t NBLOCKS = mdata->NBLOCKS;

    // needed in loops
    __m512i i0, i1;
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19
    __m512i c00, c01, c02, c03, c04, c05, c06, c07, c08, c09, c10, c11;

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
    __m512i nhatvec_e = _mm512_load_epi64(mdata->vrho);
    __m512i hiword = _mm512_set1_epi64(0x000000000000001);
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry_e = 0;
    __mmask8 scarry2;
    __mmask8 scarry;

    // deal with the sign
    c->size = NWORDS;
    c->signmask = 0;

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

    // first half sqr
    i = 0;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        //for (k = 0, j = i; j > (i + 1) / 2; j--, k++)
        k = 0;
        j = 0;
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
            VEC_SQR_MUL4_C(a0, a1, a2, a3);


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
            VEC_SQR_MUL2_C(a0, a1);

#endif
        }

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations

#ifndef IFMA
        // need to remove bias from the s*n loop and the 
        // two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);
#endif

        // final monty column accumulation
        {
            j = 0;
            // accumulate this column-sum
            // now carry propagate low to high.
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);

            //_mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            c00 = a0;
            b0 = _mm512_load_epi64(n->data + 0 * VECLEN);

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            b1 = _mm512_load_epi64(n->data + 1 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c00, b1, acc_e0, acc_e1);

            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            //_mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            c01 = a0;

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            b1 = _mm512_load_epi64(n->data + 2 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c00, b1, acc_e0, acc_e1);
            b1 = _mm512_load_epi64(n->data + 1 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c01, b1, acc_e0, acc_e1);

            // these are bottleneck high-latency sequentially dependent instructions...
            // not sure what can be done about it.
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            //_mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            c02 = a0;

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            b1 = _mm512_load_epi64(n->data + 3 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c00, b1, acc_e0, acc_e1);
            b1 = _mm512_load_epi64(n->data + 2 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c01, b1, acc_e0, acc_e1);
            b1 = _mm512_load_epi64(n->data + 1 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c02, b1, acc_e0, acc_e1);

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            //_mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            c03 = a0;

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;
        }
    }

    i = 1;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        //for (k = 0, j = i; j > (i + 1) / 2; j--, k++)
        k = 0;
        j = 1;
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

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        j = 0;
        {
            // accumulate s * n
            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c03, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c02, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c01, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c00, b3, b4, b5, b6);
        }

#ifndef IFMA
        // need to remove bias from the s*n loop and the 
        // two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);
#endif

        // final monty column accumulation
        {
            j = 0;
            // accumulate this column-sum
            // now carry propagate low to high.
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);

            //_mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            c04 = a0;
            b0 = _mm512_load_epi64(n->data + 0 * VECLEN);

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            b1 = _mm512_load_epi64(n->data + 1 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c04, b1, acc_e0, acc_e1);

            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            //_mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            c05 = a0;

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            b1 = _mm512_load_epi64(n->data + 2 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c04, b1, acc_e0, acc_e1);
            b1 = _mm512_load_epi64(n->data + 1 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c05, b1, acc_e0, acc_e1);

            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            //_mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            c06 = a0;

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            b1 = _mm512_load_epi64(n->data + 3 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c04, b1, acc_e0, acc_e1);
            b1 = _mm512_load_epi64(n->data + 2 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c05, b1, acc_e0, acc_e1);
            b1 = _mm512_load_epi64(n->data + 1 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c06, b1, acc_e0, acc_e1);

            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            //_mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            c07 = a0;

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;
        }
    }

    i = 2;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        //for (k = 0, j = i; j > (i + 1) / 2; j--, k++)
        j = 2;
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

        k = 1;
        j = 1;
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

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        j = 0;
        {
            // accumulate s * n
            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c03, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c02, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c01, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c00, b3, b4, b5, b6);
        }

        j = 1;
        {
            // accumulate s * n
            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c07, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c06, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c05, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c04, b3, b4, b5, b6);
        }

#ifndef IFMA
        // need to remove bias from the s*n loop and the 
        // two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);
#endif

        // final monty column accumulation
        {
            j = 0;
            // accumulate this column-sum
            // now carry propagate low to high.
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);

            //_mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            c08 = a0;
            b0 = _mm512_load_epi64(n->data + 0 * VECLEN);

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            b1 = _mm512_load_epi64(n->data + 1 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c08, b1, acc_e0, acc_e1);

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            //_mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            c09 = a0;

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            b1 = _mm512_load_epi64(n->data + 2 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c08, b1, acc_e0, acc_e1);
            b1 = _mm512_load_epi64(n->data + 1 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c09, b1, acc_e0, acc_e1);

            // these are bottleneck high-latency sequentially dependent instructions...
            // not sure what can be done about it.
            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            //_mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            c10 = a0;

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            b1 = _mm512_load_epi64(n->data + 3 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c08, b1, acc_e0, acc_e1);
            b1 = _mm512_load_epi64(n->data + 2 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c09, b1, acc_e0, acc_e1);
            b1 = _mm512_load_epi64(n->data + 1 * VECLEN);
            VEC_MUL_ACCUM_LOHI_PD(c10, b1, acc_e0, acc_e1);

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            //_mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            c11 = a0;

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;
        }
    }

    //printf("fixed624: \n");
    //print_regvechex(c11, 0, "0:");
    //print_regvechex(c10, 0, "0:");
    //print_regvechex(c09, 0, "0:");
    //print_regvechex(c08, 0, "0:");
    //print_regvechex(c07, 0, "0:");
    //print_regvechex(c06, 0, "0:");
    //print_regvechex(c05, 0, "0:");
    //print_regvechex(c04, 0, "0:");
    //print_regvechex(c03, 0, "0:");
    //print_regvechex(c02, 0, "0:");
    //print_regvechex(c01, 0, "0:");
    //print_regvechex(c00, 0, "0:");
    //printf("\n");

    // second half sqr
    i = 0;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        //for (k = 0, j = 0; j < (NBLOCKS - i - 1) / 2; j++, k++)
        j = 0;
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
        // NBLOCKS is odd
        j = 1;
        k = 1;
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


        // the s*n terms.  No more doubling past here.
        //for (j = 0; j < NBLOCKS - 1 - i; j++)
        j = 0;
        {
            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c11, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c10, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c09, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c08, b3, b4, b5, b6);
        }

        j = 1;
        {
            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c07, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c06, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c05, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c04, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum (s * n)
        a1 = c01;
        a2 = c02;
        a3 = c03;

        j = 2;
        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a1, b0, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b1, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b0, te2, te3);
#else
        // finish each triangluar shaped column sum (s * n)
        // a1*b0 -> t0/1
        // a2*b1 -> t0/1
        // a2*b0 -> t2/3
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

        //VEC_MUL_ACCUM_LOHI(a3, b2, te0, te1);
        //VEC_MUL_ACCUM_LOHI(a3, b1, te2, te3);
        //VEC_MUL_ACCUM_LOHI(a3, b0, te4, te5);
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
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);
        SUB_BIAS_LO(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);

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

            _mm512_store_epi64(c->data + (i * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(c->data + (i * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(c->data + (i * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(c->data + (i * BLOCKWORDS + 3) * VECLEN, a0);
        }


    }

    i = 1;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        // The final block shape depends on the parity of i and NBLOCKS
        // Block shape 1 (small upper triangle) if 'i' is odd and NBLOCKS is even,
        // or if 'i' is even and NBLOCKS is odd.
        // Block shape 2 if 'i' is odd and NBLOCKS is odd or if 'i' is even
        // and NBLOCKS is even.
        // NBLOCKS is odd

        //for (k = 0, j = 0; j < (NBLOCKS - i - 1) / 2; j++, k++)

        // i odd, block shape 2.
        k = 0;
        j = 0;
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


        // the s*n terms.  No more doubling past here.
        j = 0;
        {
            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c11, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c10, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c09, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c08, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum (s * n)
        a1 = c05;
        a2 = c06;
        a3 = c07;

        j = 1;
        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a1, b0, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b1, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b0, te2, te3);
#else
        // finish each triangluar shaped column sum (s * n)
        // a1*b0 -> t0/1
        // a2*b1 -> t0/1
        // a2*b0 -> t2/3
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

        //VEC_MUL_ACCUM_LOHI(a3, b2, te0, te1);
        //VEC_MUL_ACCUM_LOHI(a3, b1, te2, te3);
        //VEC_MUL_ACCUM_LOHI(a3, b0, te4, te5);
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
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);
        SUB_BIAS_LO(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);

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

            _mm512_store_epi64(c->data + (i * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(c->data + (i * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(c->data + (i * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(c->data + (i * BLOCKWORDS + 3) * VECLEN, a0);
        }


    }

    i = 2;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        // The final block shape depends on the parity of i and NBLOCKS
        // Block shape 1 (small upper triangle) if 'i' is odd and NBLOCKS is even,
        // or if 'i' is even and NBLOCKS is odd.
        // Block shape 2 if 'i' is odd and NBLOCKS is odd or if 'i' is even
        // and NBLOCKS is even.
        // NBLOCKS is odd

        //for (k = 0, j = 0; j < (NBLOCKS - i - 1) / 2; j++, k++)
        k = 0;
        j = 0;
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


        // finish each triangluar shaped column sum (s * n)
        a1 = c09;
        a2 = c10;
        a3 = c11;

        j = 0;
        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a1, b0, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b1, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b0, te2, te3);
#else
        // finish each triangluar shaped column sum (s * n)
        // a1*b0 -> t0/1
        // a2*b1 -> t0/1
        // a2*b0 -> t2/3
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

        //VEC_MUL_ACCUM_LOHI(a3, b2, te0, te1);
        //VEC_MUL_ACCUM_LOHI(a3, b1, te2, te3);
        //VEC_MUL_ACCUM_LOHI(a3, b0, te4, te5);
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
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);
        SUB_BIAS_LO(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);

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

            _mm512_store_epi64(c->data + (i * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(c->data + (i * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(c->data + (i * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(c->data + (i * BLOCKWORDS + 3) * VECLEN, a0);
        }


    }

    //printf("fixed624: \n");
    //for (i = NWORDS - 1; i >= 0; i--)
    //{
    //    printf("%016lx", c->data[i * VECLEN + 0]);
    //}
    //printf("\n");

    c->size = NWORDS;
    return;
}

#endif

void vecmulmod52_avxecm(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    int i, j, k;
    int NWORDS = mdata->NWORDS;
    int NBLOCKS = mdata->NBLOCKS;

    // needed in loops
    __m512i i0, i1;
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19
    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i nhatvec_e = _mm512_load_epi64(mdata->vrho);
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry2;
    __mmask8 scarry;

    // deal with the sign
    c->size = NWORDS;
    c->signmask = a->signmask ^ b->signmask;

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

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        for (j = 0; j < i; j++)
        {
            // accumulate s * n
            a0 = _mm512_load_epi64(s->data + ((j + 1) * BLOCKWORDS - 1) * VECLEN);
            a1 = _mm512_load_epi64(s->data + ((j + 1) * BLOCKWORDS - 2) * VECLEN);
            a2 = _mm512_load_epi64(s->data + ((j + 1) * BLOCKWORDS - 3) * VECLEN);
            a3 = _mm512_load_epi64(s->data + ((j + 1) * BLOCKWORDS - 4) * VECLEN);

            b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);
        SUB_BIAS_LO(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);

        // now, column by column, add in the s*n contribution and reduce to 
        // a single 64+x bit accumulator while storing the intermediate product
        // 's' as we go.
        {
            j = 0;
            // accumulate this column-sum
            // now carry propagate low to high.
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);

            _mm512_store_epi64(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);
            b0 = _mm512_load_epi64(n->data + 0 * VECLEN);

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            for (k = 0; k < j; k++)
            {
                a0 = _mm512_load_epi64(s->data + (i * BLOCKWORDS + k) * VECLEN);
                b1 = _mm512_load_epi64(n->data + (j - k) * VECLEN);

                VEC_MUL_ACCUM_LOHI_PD(a0, b1, acc_e0, acc_e1);
            }

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            for (k = 0; k < j; k++)
            {
                a0 = _mm512_load_epi64(s->data + (i * BLOCKWORDS + k) * VECLEN);
                b1 = _mm512_load_epi64(n->data + (j - k) * VECLEN);

                VEC_MUL_ACCUM_LOHI_PD(a0, b1, acc_e0, acc_e1);
            }

            // these are bottleneck high-latency sequentially dependent instructions...
            // not sure what can be done about it.
            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            for (k = 0; k < j; k++)
            {
                a0 = _mm512_load_epi64(s->data + (i * BLOCKWORDS + k) * VECLEN);
                b1 = _mm512_load_epi64(n->data + (j - k) * VECLEN);

                VEC_MUL_ACCUM_LOHI_PD(a0, b1, acc_e0, acc_e1);
            }

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;
        }
    }

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

            // accumulate s * n
            a0 = _mm512_load_epi64(s->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(s->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(s->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(s->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

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

        // finish each triangluar shaped column sum (s * n)
        // a1*b0 -> t0/1
        // a2*b1 -> t0/1
        // a2*b0 -> t2/3
        {
            a1 = _mm512_load_epi64(s->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
            a2 = _mm512_load_epi64(s->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
            a3 = _mm512_load_epi64(s->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

            b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

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
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);
        SUB_BIAS_LO(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);

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

            _mm512_store_epi64(s->data + ((i - NBLOCKS) * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(s->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(s->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(s->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN, a0);

        }
    }

    a0 = acc_e0;
    scarry2 = _mm512_cmp_epu64_mask(a0, zero, _MM_CMPINT_EQ);

    // subtract n from tmp
    scarry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        a1 = _mm512_load_epi64(s->data + i * VECLEN);
        b0 = _mm512_load_epi64(n->data + i * VECLEN);
        a0 = _mm512_sbb_epi52(a1, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
    }

    // negate any final borrows if there was also a final carry.
    scarry &= scarry2;

    // if there was a final borrow, we didn't need to do the subtraction after all.
    // replace with original results based on final borrow mask.
    for (i = NWORDS - 1; i >= 0; i--)
    {
        b0 = _mm512_load_epi64(s->data + i * VECLEN);
        _mm512_mask_store_epi64(c->data + i * VECLEN, scarry, b0);
    }

    c->size = NWORDS;
    return;
}

void vecmulmod52(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    int i, j, k;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t NBLOCKS = mdata->NBLOCKS;

    // needed in loops
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
    __m512i nhatvec_e = _mm512_load_epi64(mdata->vrho);
    //__m512i hiword = _mm512_set1_epi64(0x000000000000001);
    __m512i zero = _mm512_set1_epi64(0);
    //__mmask8 scarry_e = 0;
    __mmask8 scarry2;
    __mmask8 scarry;

#ifdef USE_AMM
    uint64_t* outdata = c->data;
#else
    uint64_t* outdata = s->data;
#endif

    // deal with the sign
    c->size = NWORDS;
    c->signmask = a->signmask ^ b->signmask;

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
        //VEC_MUL_ACCUM_LOHI(a0, b3, te0, te1);
        //VEC_MUL_ACCUM_LOHI(a1, b3, te2, te3);
        //VEC_MUL_ACCUM_LOHI(a0, b2, te2, te3);
        //VEC_MUL_ACCUM_LOHI(a2, b3, te4, te5);
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
        //VEC_MUL_ACCUM_LOHI(a1, b2, te4, te5);
        //VEC_MUL_ACCUM_LOHI(a0, b1, te4, te5);
        //VEC_MUL_ACCUM_LOHI(a1, b1, te6, te7);
        //VEC_MUL_ACCUM_LOHI(a0, b0, te6, te7);
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
        //VEC_MUL_ACCUM_LOHI(a3, b3, te6, te7);
        //VEC_MUL_ACCUM_LOHI(a2, b2, te6, te7);
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
#endif

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        for (j = 0; j < i; j++)
        {
            // accumulate s * n
            a0 = _mm512_load_epi64(outdata + ((j + 1) * BLOCKWORDS - 1) * VECLEN);
            a1 = _mm512_load_epi64(outdata + ((j + 1) * BLOCKWORDS - 2) * VECLEN);
            a2 = _mm512_load_epi64(outdata + ((j + 1) * BLOCKWORDS - 3) * VECLEN);
            a3 = _mm512_load_epi64(outdata + ((j + 1) * BLOCKWORDS - 4) * VECLEN);

            b0 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

#ifndef IFMA

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);
        SUB_BIAS_LO(
            i * 8 + 1,
            i * 8 + 2,
            i * 8 + 3,
            i * 8 + 4);
#endif

        // now, column by column, add in the s*n contribution and reduce to 
        // a single 64+x bit accumulator while storing the intermediate product
        // 's' as we go.
        {
            j = 0;
            // accumulate this column-sum
            // now carry propagate low to high.
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);

            _mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            b0 = _mm512_load_epi64(n->data + 0 * VECLEN);

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            for (k = 0; k < j; k++)
            {
                a0 = _mm512_load_epi64(outdata + (i * BLOCKWORDS + k) * VECLEN);
                b1 = _mm512_load_epi64(n->data + (j - k) * VECLEN);

                VEC_MUL_ACCUM_LOHI_PD(a0, b1, acc_e0, acc_e1);
            }

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            _mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            for (k = 0; k < j; k++)
            {
                a0 = _mm512_load_epi64(outdata + (i * BLOCKWORDS + k) * VECLEN);
                b1 = _mm512_load_epi64(n->data + (j - k) * VECLEN);

                VEC_MUL_ACCUM_LOHI_PD(a0, b1, acc_e0, acc_e1);
            }

            // these are bottleneck high-latency sequentially dependent instructions...
            // not sure what can be done about it.
            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            _mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            for (k = 0; k < j; k++)
            {
                a0 = _mm512_load_epi64(outdata + (i * BLOCKWORDS + k) * VECLEN);
                b1 = _mm512_load_epi64(n->data + (j - k) * VECLEN);

                VEC_MUL_ACCUM_LOHI_PD(a0, b1, acc_e0, acc_e1);
            }

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            _mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;
        }
    }

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

            // accumulate s * n
            // i = 3, j = 1, a0..3 = c[11..8]
            // i = 3, j = 2, a0..3 = c[7..4]
            // i = 4, j = 2, a0..3 = c[11..8]
            a0 = _mm512_load_epi64(outdata + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(outdata + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(outdata + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(outdata + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

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
        //VEC_MUL_ACCUM_LOHI(a1, b0, te0, te1);
        //VEC_MUL_ACCUM_LOHI(a2, b1, te0, te1);
        //VEC_MUL_ACCUM_LOHI(a2, b0, te2, te3);
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
        //VEC_MUL_ACCUM_LOHI(a3, b2, te0, te1);
        //VEC_MUL_ACCUM_LOHI(a3, b1, te2, te3);
        //VEC_MUL_ACCUM_LOHI(a3, b0, te4, te5);
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
#endif

        // i = 3, a1..3 = c[1..3]
        // i = 4, a1..3 = c[5..7]
        // i = 5, a1..3 = c[9..11]
        a1 = _mm512_load_epi64(outdata + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(outdata + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(outdata + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        // a1*b0 -> t0/1
        // a2*b1 -> t0/1
        // a2*b0 -> t2/3
#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a1, b0, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b1, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b0, te2, te3);
#else
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
        //VEC_MUL_ACCUM_LOHI(a3, b2, te0, te1);
        //VEC_MUL_ACCUM_LOHI(a3, b1, te2, te3);
        //VEC_MUL_ACCUM_LOHI(a3, b0, te4, te5);
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
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);
        SUB_BIAS_LO(
            (2 * NBLOCKS - i - 1) * 8 + 6,
            (2 * NBLOCKS - i - 1) * 8 + 4,
            (2 * NBLOCKS - i - 1) * 8 + 2,
            (2 * NBLOCKS - i - 1) * 8 + 0);

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

            // i = 3, a3..0 = c[0..3]
            // i = 4, a3..0 = c[4..7]
            // i = 5, a3..0 = c[8..11]
            _mm512_store_epi64(outdata + ((i - NBLOCKS) * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(outdata + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(outdata + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(outdata + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN, a0);
        }
    }


#ifdef USE_AMM
    //for (i = NWORDS - 1; i >= 0; i--)
    //{
    //    b0 = _mm512_load_epi64(s->data + i * VECLEN);
    //    _mm512_store_epi64(c->data + i * VECLEN, b0);
    //}
#else
    a0 = acc_e0;
    scarry2 = _mm512_cmp_epu64_mask(a0, zero, _MM_CMPINT_EQ);

    // subtract n from tmp
    scarry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        a1 = _mm512_load_epi64(s->data + i * VECLEN);
        b0 = _mm512_load_epi64(n->data + i * VECLEN);
        a0 = _mm512_sbb_epi52(a1, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
    }

    // negate any final borrows if there was also a final carry.
    scarry &= scarry2;

    // if there was a final borrow, we didn't need to do the subtraction after all.
    // replace with original results based on final borrow mask.
    for (i = NWORDS - 1; i >= 0; i--)
    {
        b0 = _mm512_load_epi64(s->data + i * VECLEN);
        _mm512_mask_store_epi64(c->data + i * VECLEN, scarry, b0);
    }
#endif

    c->size = NWORDS;
    return;
}

void vecredc52_base(vec_bignum_t *a, vec_bignum_t *c, vec_bignum_t *n, vec_bignum_t *s, vec_monty_t*mdata)
{
    // taking a number 'a' out of Montgomery representation: REDC on an
    // input of size NWORDS.
    int i, j, k;
    uint32_t NWORDS = n->WORDS_ALLOC;
    // needed in loops
    __m512i i0, i1;
    __m512i a0, a1, carry, q;                                     // 4
    __m512i b0;
    __m512d prod1_hd, prod2_hd;                 // 23
    __m512d prod1_ld, prod2_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1;
    __m512i nhatvec_e = _mm512_load_epi64(mdata->vrho);
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry2;
    __mmask8 scarry;

    //{
    //    q = t1[0] * nhatvec_e;
    //    t2[n] = mpn_mul_1(t2, n, NWORDS, q);
    //
    //    q = mpn_add_n(c, t1 + 1, t2 + 1, NWORDS);
    //    q += mpn_add_1(c, c, NWORDS, t1[0] != 0);
    //
    //    while (q != 0)
    //        q -= mpn_sub_n(c, c, n, NWORDS);
    //}

    a0 = _mm512_load_epi64(a->data);
    _mm512_mullo_epi52(q, nhatvec_e, a0);

    carry = zero;
    a0 = _mm512_load_epi64(n->data);
    b0 = _mm512_load_epi64(a->data);
    VEC_MUL_LOHI_PD(a0, q, acc_e0, acc_e1);
    acc_e0 = _mm512_add_epi64(acc_e0, carry);
    acc_e0 = _mm512_add_epi64(acc_e0, b0);
    acc_e0 = _mm512_srli_epi64(acc_e0, 52);
    acc_e1 = _mm512_add_epi64(acc_e1, acc_e0);
    carry = acc_e1;
    for (i = 1; i < NWORDS; i++)
    {
        a0 = _mm512_load_epi64(n->data + i * VECLEN);
        b0 = _mm512_load_epi64(a->data + i * VECLEN);
        VEC_MUL_LOHI_PD(a0, q, acc_e0, acc_e1);
        acc_e0 = _mm512_add_epi64(acc_e0, carry);
        acc_e0 = _mm512_add_epi64(acc_e0, b0);
        carry = acc_e1;
        _mm512_store_epi64(s->data + (i - 1) * VECLEN, acc_e0);
    }
    _mm512_store_epi64(s->data + (i - 1) * VECLEN, carry);
    acc_e0 = _mm512_srli_epi64(carry, 52);

    a0 = acc_e0;
    scarry2 = _mm512_cmp_epu64_mask(a0, zero, _MM_CMPINT_EQ);

    // subtract n from tmp (need to more than once?)
    scarry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        a1 = _mm512_load_epi64(s->data + i * VECLEN);
        b0 = _mm512_load_epi64(n->data + i * VECLEN);
        a0 = _mm512_sbb_epi52(a1, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
    }

    // negate any final borrows if there was also a final carry.
    scarry &= scarry2;

    // if there was a final borrow, we didn't need to do the subtraction after all.
    // replace with original results based on final borrow mask.
    for (i = NWORDS - 1; i >= 0; i--)
    {
        b0 = _mm512_load_epi64(s->data + i * VECLEN);
        _mm512_mask_store_epi64(c->data + i * VECLEN, scarry, b0);
    }

    c->size = NWORDS;
    return;
}

void vecmulmod52_1(vec_bignum_t *a, base_t *b, vec_bignum_t *c, vec_bignum_t *n, vec_bignum_t *s, vec_monty_t*mdata)
{
    int i, j, k;
    uint32_t NWORDS = n->WORDS_ALLOC;
    // needed in loops
    __m512i i0, i1;
    __m512i a0, a1, carry, q;                                     // 4
    __m512i b0;
    __m512d prod1_hd, prod2_hd;                 // 23
    __m512d prod1_ld, prod2_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1;
    __m512i nhatvec_e = _mm512_load_epi64(mdata->vrho);
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry2;
    __mmask8 scarry;

    // deal with the sign
    c->size = NWORDS;
    c->signmask = a->signmask;

    //{
    //    t1[n] = mpn_mul_1(t1, a, NWORDS, b);
    //    q = t1[0] * nhatvec_e;
    //    t2[n] = mpn_mul_1(t2, n, NWORDS, q);
    //
    //    q = mpn_add_n(c, t1 + 1, t2 + 1, NWORDS);
    //    q += mpn_add_1(c, c, NWORDS, t1[0] != 0);
    //
    //    while (q != 0)
    //        q -= mpn_sub_n(c, c, n, NWORDS);
    //}

    carry = zero;
    b0 = _mm512_load_epi64(b);
    a0 = _mm512_load_epi64(a->data);
    VEC_MUL_LOHI_PD(a0, b0, acc_e0, acc_e1);
    acc_e0 = _mm512_add_epi64(acc_e0, carry);
    carry = acc_e1;
    _mm512_mullo_epi52(q, nhatvec_e, a0);
    _mm512_store_epi64(s->data, acc_e0);
    
    for (i = 1; i < NWORDS; i++)
    {
        a0 = _mm512_load_epi64(a->data + i * VECLEN);
        VEC_MUL_LOHI_PD(a0, b0, acc_e0, acc_e1);
        acc_e0 = _mm512_add_epi64(acc_e0, carry);
        carry = acc_e1;
        _mm512_store_epi64(s->data + i * VECLEN, acc_e0);
    }
    _mm512_store_epi64(s->data + i * VECLEN, carry);

    carry = zero;
    a0 = _mm512_load_epi64(n->data);
    b0 = _mm512_load_epi64(s->data);
    VEC_MUL_LOHI_PD(a0, q, acc_e0, acc_e1);
    acc_e0 = _mm512_add_epi64(acc_e0, carry);
    acc_e0 = _mm512_add_epi64(acc_e0, b0);
    acc_e0 = _mm512_srli_epi64(acc_e0, 52);
    acc_e1 = _mm512_add_epi64(acc_e1, acc_e0);
    carry = acc_e1;
    for (i = 1; i < NWORDS; i++)
    {
        a0 = _mm512_load_epi64(n->data + i * VECLEN);
        b0 = _mm512_load_epi64(s->data + i * VECLEN);
        VEC_MUL_LOHI_PD(a0, q, acc_e0, acc_e1);
        acc_e0 = _mm512_add_epi64(acc_e0, carry);
        acc_e0 = _mm512_add_epi64(acc_e0, b0);
        carry = acc_e1;
        _mm512_store_epi64(s->data + (i - 1) * VECLEN, acc_e0);
    }
    _mm512_store_epi64(s->data + (i - 1) * VECLEN, carry);
    acc_e0 = _mm512_srli_epi64(carry, 52);

    a0 = acc_e0;

    if (_mm512_cmp_epu64_mask(a0, _mm512_set1_epi64(1), _MM_CMPINT_GT))
    {
        printf("carry greater than 1!\n");
    }

    scarry2 = _mm512_cmp_epu64_mask(a0, zero, _MM_CMPINT_EQ);

    // subtract n from tmp (need to more than once?)
    scarry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        a1 = _mm512_load_epi64(s->data + i * VECLEN);
        b0 = _mm512_load_epi64(n->data + i * VECLEN);
        a0 = _mm512_sbb_epi52(a1, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
    }

    // negate any final borrows if there was also a final carry.
    scarry &= scarry2;

    // if there was a final borrow, we didn't need to do the subtraction after all.
    // replace with original results based on final borrow mask.
    for (i = NWORDS - 1; i >= 0; i--)
    {
        b0 = _mm512_load_epi64(s->data + i * VECLEN);
        _mm512_mask_store_epi64(c->data + i * VECLEN, scarry, b0);
    }

    c->size = NWORDS;
    return;
}

void vecmul52_1(vec_bignum_t* a, __m512i b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    // unsigned multiply a * b, where b is a single base_t digit.
    // resulting in a NWORDS+1 size result.
    int i;
    uint32_t NWORDS = n->WORDS_ALLOC;
    // needed in loops
    __m512i a0, carry;                                     // 4
    __m512d prod1_hd;                 // 23
    __m512d prod1_ld, prod2_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1;
    __m512i zero = _mm512_set1_epi64(0);

    carry = zero;
    a0 = _mm512_load_epi64(a->data);
    VEC_MUL_LOHI_PD(a0, b, acc_e0, acc_e1);
    carry = acc_e1;
    _mm512_store_epi64(c->data, acc_e0);

    for (i = 1; i < NWORDS; i++)
    {
        a0 = _mm512_load_epi64(a->data + i * VECLEN);
        VEC_MUL_LOHI_PD(a0, b, acc_e0, acc_e1);
        acc_e0 = _mm512_add_epi64(acc_e0, carry);
        VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
        carry = acc_e1;
        _mm512_store_epi64(c->data + i * VECLEN, acc_e0);
    }
    _mm512_store_epi64(c->data + i * VECLEN, carry);

    c->size = NWORDS + 1;
    return;
}

void vecsqrmod52_fixed1040_bfips(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    // 8x sqr:
    // input 8 bignums in the even lanes of a.
    // output 8 squaremod bignums in the even lanes of c.
    int i, j, k;
    vec_bignum_t* b = a;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t NBLOCKS = mdata->NBLOCKS;

    // needed in loops
    __m512i i0, i1;
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19
    __m512i c00, c01, c02, c03, c04, c05, c06, c07, c08, c09, c10, c11,
        c12, c13, c14, c15, c16, c17, c18, c19;


    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31


    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i nhatvec_e = _mm512_load_epi64(mdata->vrho);
    __m512i hiword = _mm512_set1_epi64(0x000000000000001);
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry_e = 0;
    __mmask8 scarry2;
    __mmask8 scarry;

    // deal with the sign
    c->size = NWORDS;
    c->signmask = 0;

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

    // first half sqr
    i = 0;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        //for (k = 0, j = i; j > (i + 1) / 2; j--, k++)
        k = 0;
        j = 0;

        // i even
        {          
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);

            VEC_SQR_MUL4_C(a0, a1, a2, a3);


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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally, accumulate the two non-doubled terms.
            VEC_SQR_MUL2_C(a0, a1);
        }

        // need to remove bias from the s*n loop and the 
        // two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);

        // final monty column accumulation
        CARRYPROP_ACCUM_0(c00);
        CARRYPROP_ACCUM_1(c01, c00);
        CARRYPROP_ACCUM_2(c02, c00, c01);
        CARRYPROP_ACCUM_3(c03, c00, c01, c02);
    }

    i = 1;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        //for (k = 0, j = i; j > (i + 1) / 2; j--, k++)
        k = 0;
        j = 1;

        // i odd
        {
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

            VEC_SQR_MUL4_A(a0, a1, a2, b2, b3);
            VEC_SQR_MUL4_A(a1, a2, a3, b3, b4);
            VEC_SQR_MUL4_B(a2, a3, b4, b5, b6);

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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally, accumulate the two non-doubled terms.
            VEC_SQR_MUL2(a0, a1);

        }

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        j = 0;
        {
            // accumulate s * n
            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c03, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c02, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c01, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c00, b3, b4, b5, b6);
        }

        // need to remove bias from the s*n loop and the 
        // two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);

        // final monty column accumulation
        CARRYPROP_ACCUM_0(c04);
        CARRYPROP_ACCUM_1(c05, c04);
        CARRYPROP_ACCUM_2(c06, c04, c05);
        CARRYPROP_ACCUM_3(c07, c04, c05, c06);
    }

    i = 2;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (k = 0, j = i; j > (i + 1) / 2; j--, k++)
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

        // i even
        {
            
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);

            VEC_SQR_MUL4_C(a0, a1, a2, a3);

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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally, accumulate the two non-doubled terms.
            VEC_SQR_MUL2_C(a0, a1);
        }

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        j = 0;
        {
            // accumulate s * n
            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c03, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c02, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c01, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c00, b3, b4, b5, b6);
        }

        j = 1;
        {
            // accumulate s * n
            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c07, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c06, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c05, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c04, b3, b4, b5, b6);
        }

        // need to remove bias from the s*n loop and the 
        // two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);

        // final monty column accumulation
        CARRYPROP_ACCUM_0(c08);
        CARRYPROP_ACCUM_1(c09, c08);
        CARRYPROP_ACCUM_2(c10, c08, c09);
        CARRYPROP_ACCUM_3(c11, c08, c09, c10);
    }

    i = 3;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (k = 0, j = i; j > (i + 1) / 2; j--, k++)
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

        // i odd
        {
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

            VEC_SQR_MUL4_A(a0, a1, a2, b2, b3);
            VEC_SQR_MUL4_A(a1, a2, a3, b3, b4);
            VEC_SQR_MUL4_B(a2, a3, b4, b5, b6);

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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally, accumulate the two non-doubled terms.
            VEC_SQR_MUL2(a0, a1);

        }

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        j = 0;
        {
            // accumulate s * n
            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c03, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c02, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c01, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c00, b3, b4, b5, b6);
        }

        j = 1;
        {
            // accumulate s * n
            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c07, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c06, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c05, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c04, b3, b4, b5, b6);
        }

        j = 2;
        {
            // accumulate s * n
            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c11, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c10, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c09, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c08, b3, b4, b5, b6);
        }

        // need to remove bias from the s*n loop and the 
        // two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);

        // final monty column accumulation
        CARRYPROP_ACCUM_0(c12);
        CARRYPROP_ACCUM_1(c13, c12);
        CARRYPROP_ACCUM_2(c14, c12, c13);
        CARRYPROP_ACCUM_3(c15, c12, c13, c14);
    }

    i = 4;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (k = 0, j = i; j > (i + 1) / 2; j--, k++)
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

        // i even
        {

            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);

            VEC_SQR_MUL4_C(a0, a1, a2, a3);

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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally, accumulate the two non-doubled terms.
            VEC_SQR_MUL2_C(a0, a1);
        }

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        j = 0;
        {
            // accumulate s * n
            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c03, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c02, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c01, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c00, b3, b4, b5, b6);
        }

        j = 1;
        {
            // accumulate s * n
            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c07, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c06, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c05, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c04, b3, b4, b5, b6);
        }

        j = 2;
        {
            // accumulate s * n
            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c11, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c10, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c09, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c08, b3, b4, b5, b6);
        }

        j = 3;
        {
            // accumulate s * n
            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c15, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c14, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c13, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c12, b3, b4, b5, b6);
        }

        // need to remove bias from the s*n loop and the 
        // two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);

        // final monty column accumulation
        CARRYPROP_ACCUM_0(c16);
        CARRYPROP_ACCUM_1(c17, c16);
        CARRYPROP_ACCUM_2(c18, c16, c17);
        CARRYPROP_ACCUM_3(c19, c16, c17, c18);
    }

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

        // NBLOCKS is odd
        // i even, block shape 1.
        {
            
            // always a continuation of the full-block loop, so use the same 
            // loading pattern.  Only now we don't need as many b-terms.
            a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
            a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
            a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);

            VEC_SQR_MUL2_D(a0, a1, a2);

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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally the two non-doubled terms.
            VEC_SQR_MUL2_C(a1, a0);
        }


        // the s*n terms.  No more doubling past here.
        //for (j = 0; j < NBLOCKS - 1 - i; j++)
        j = 0;
        {
            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c19, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c18, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c17, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c16, b3, b4, b5, b6);
        }

        j = 1;
        {
            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c15, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c14, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c13, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c12, b3, b4, b5, b6);
        }

        j = 2;
        {
            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c11, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c10, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c09, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c08, b3, b4, b5, b6);
        }

        j = 3;
        {
            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c07, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c06, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c05, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c04, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum (s * n)
        a1 = c01;
        a2 = c02;
        a3 = c03;

        j = 4;
        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);
        SUB_BIAS_LO(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);

        // final column accumulation
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 0) * VECLEN, a3);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 1) * VECLEN, a2);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 2) * VECLEN, a1);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 3) * VECLEN, a0);

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

        // NBLOCKS is odd
        // i odd, block shape 2.
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

            VEC_SQR_MUL4_D(a0, a1, b0, b1, b2, te0, te1, te2, te3);
            VEC_SQR_MUL2_B(a2, b2, b3);
            VEC_SQR_MUL4_D(a0, a1, b2, b3, b4, te4, te5, te6, te7);

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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally the two non-doubled terms.
            VEC_SQR_MUL2_C(a3, a2);

        }


        // the s*n terms.  No more doubling past here.
        //for (j = 0; j < NBLOCKS - 1 - i; j++)
        j = 0;
        {
            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c19, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c18, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c17, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c16, b3, b4, b5, b6);
        }

        j = 1;
        {
            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c15, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c14, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c13, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c12, b3, b4, b5, b6);
        }

        j = 2;
        {
            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c11, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c10, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c09, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c08, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum (s * n)
        a1 = c05;
        a2 = c06;
        a3 = c07;

        j = 3;
        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);
        SUB_BIAS_LO(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);

        // final column accumulation
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 0) * VECLEN, a3);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 1) * VECLEN, a2);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 2) * VECLEN, a1);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 3) * VECLEN, a0);

    }

    i = 2;
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

        // NBLOCKS is odd
        // i even, block shape 1.
        {
            
            // always a continuation of the full-block loop, so use the same 
            // loading pattern.  Only now we don't need as many b-terms.
            a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
            a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
            a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);

            VEC_SQR_MUL2_D(a0, a1, a2);

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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally the two non-doubled terms.
            VEC_SQR_MUL2_C(a1, a0);
        }

        //for (j = 0; j < NBLOCKS - 1 - i; j++)
        j = 0;
        {
            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c19, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c18, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c17, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c16, b3, b4, b5, b6);
        }

        j = 1;
        {
            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c15, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c14, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c13, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c12, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum (s * n)
        a1 = c09;
        a2 = c10;
        a3 = c11;

        j = 2;
        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);
        SUB_BIAS_LO(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);

        // final column accumulation
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 0) * VECLEN, a3);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 1) * VECLEN, a2);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 2) * VECLEN, a1);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 3) * VECLEN, a0);

    }

    i = 3;
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

        // NBLOCKS is odd
        // i odd, block shape 2.
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

            VEC_SQR_MUL4_D(a0, a1, b0, b1, b2, te0, te1, te2, te3);
            VEC_SQR_MUL2_B(a2, b2, b3);
            VEC_SQR_MUL4_D(a0, a1, b2, b3, b4, te4, te5, te6, te7);

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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally the two non-doubled terms.
            VEC_SQR_MUL2_C(a3, a2);

        }


        // the s*n terms.  No more doubling past here.
        //for (j = 0; j < NBLOCKS - 1 - i; j++)
        j = 0;
        {
            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c19, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c18, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c17, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c16, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum (s * n)
        a1 = c13;
        a2 = c14;
        a3 = c15;

        j = 1;
        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);
        SUB_BIAS_LO(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);

        // final column accumulation
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 0) * VECLEN, a3);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 1) * VECLEN, a2);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 2) * VECLEN, a1);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 3) * VECLEN, a0);

    }

    i = 4;
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

        // NBLOCKS is odd
        // i even, block shape 1.
        {

            // always a continuation of the full-block loop, so use the same 
            // loading pattern.  Only now we don't need as many b-terms.
            a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
            a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
            a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);

            VEC_SQR_MUL2_D(a0, a1, a2);

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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally the two non-doubled terms.
            VEC_SQR_MUL2_C(a1, a0);
        }

        //for (j = 0; j < NBLOCKS - 1 - i; j++)

        // finish each triangluar shaped column sum (s * n)
        a1 = c17;
        a2 = c18;
        a3 = c19;

        j = 0;
        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);
        SUB_BIAS_LO(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);

        // final column accumulation
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 0) * VECLEN, a3);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 1) * VECLEN, a2);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 2) * VECLEN, a1);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 3) * VECLEN, a0);

    }
    c->size = NWORDS;
    return;
}

void vecsqrmod52_fixed832_bfips(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    // 8x sqr:
    // input 8 bignums in the even lanes of a.
    // output 8 squaremod bignums in the even lanes of c.
    int i, j, k;
    vec_bignum_t* b = a;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t NBLOCKS = mdata->NBLOCKS;

    // needed in loops
    __m512i i0, i1;
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
    __m512i nhatvec_e = _mm512_load_epi64(mdata->vrho);
    __m512i hiword = _mm512_set1_epi64(0x000000000000001);
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry_e = 0;
    __mmask8 scarry2;
    __mmask8 scarry;

    // deal with the sign
    c->size = NWORDS;
    c->signmask = 0;

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

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

        // i even
        {
            
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);

            VEC_SQR_MUL4_C(a0, a1, a2, a3);

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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally, accumulate the two non-doubled terms.
            VEC_SQR_MUL2_C(a0, a1);
        }

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        //for (j = 0; j < i; j++)

        // need to remove bias from the s*n loop and the 
        // two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);

        // final monty column accumulation
        CARRYPROP_ACCUM_0(c00);
        CARRYPROP_ACCUM_1(c01, c00);
        CARRYPROP_ACCUM_2(c02, c00, c01);
        CARRYPROP_ACCUM_3(c03, c00, c01, c02);
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

        // i odd
        {
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

            VEC_SQR_MUL4_A(a0, a1, a2, b2, b3);
            VEC_SQR_MUL4_A(a1, a2, a3, b3, b4);
            VEC_SQR_MUL4_B(a2, a3, b4, b5, b6);

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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally, accumulate the two non-doubled terms.
            VEC_SQR_MUL2(a0, a1);
        }


        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        //for (j = 0; j < i; j++)
        j = 0;
        {
            // accumulate s * n
            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c03, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c02, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c01, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c00, b3, b4, b5, b6);
        }
        j = 1;

        // need to remove bias from the s*n loop and the 
        // two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);

        // final monty column accumulation
        CARRYPROP_ACCUM_0(c04);
        CARRYPROP_ACCUM_1(c05, c04);
        CARRYPROP_ACCUM_2(c06, c04, c05);
        CARRYPROP_ACCUM_3(c07, c04, c05, c06);
    }

    i = 2;
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

        // i even
        {

            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);

            VEC_SQR_MUL4_C(a0, a1, a2, a3);

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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally, accumulate the two non-doubled terms.
            VEC_SQR_MUL2_C(a0, a1);
        }

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        //for (j = 0; j < i; j++)
        j = 0;
        {
            // accumulate s * n
            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c03, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c02, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c01, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c00, b3, b4, b5, b6);
        }

        j = 1;
        {
            // accumulate s * n
            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c07, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c06, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c05, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c04, b3, b4, b5, b6);
        }
        j = 2;

        // need to remove bias from the s*n loop and the 
        // two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);

        // final monty column accumulation
        CARRYPROP_ACCUM_0(c08);
        CARRYPROP_ACCUM_1(c09, c08);
        CARRYPROP_ACCUM_2(c10, c08, c09);
        CARRYPROP_ACCUM_3(c11, c08, c09, c10);
    }

    i = 3;
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

        // i odd
        {
            
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

            VEC_SQR_MUL4_A(a0, a1, a2, b2, b3);
            VEC_SQR_MUL4_A(a1, a2, a3, b3, b4);
            VEC_SQR_MUL4_B(a2, a3, b4, b5, b6);

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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally, accumulate the two non-doubled terms.
            VEC_SQR_MUL2(a0, a1);
        }


        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        //for (j = 0; j < i; j++)
        j = 0;
        {
            // accumulate s * n
            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c03, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c02, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c01, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c00, b3, b4, b5, b6);
        }

        j = 1;
        {
            // accumulate s * n
            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c07, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c06, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c05, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c04, b3, b4, b5, b6);
        }
        j = 2;
        {
            // accumulate s * n
            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c11, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c10, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c09, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c08, b3, b4, b5, b6);
        }
        j = 3;
        // need to remove bias from the s*n loop and the 
        // two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);

        // final monty column accumulation
        CARRYPROP_ACCUM_0(c12);
        CARRYPROP_ACCUM_1(c13, c12);
        CARRYPROP_ACCUM_2(c14, c12, c13);
        CARRYPROP_ACCUM_3(c15, c12, c13, c14);
    }

    //printf("fixed416: \n");
    //for (i = NWORDS - 1; i >= 0; i--)
    //{
    //    printf("%016lx", s->data[i * VECLEN + 0]);
    //}
    //printf("\n");

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


            VEC_SQR_MUL2_B(a2, b2, b3);

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


            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally, accumulate the two non-doubled terms.
            VEC_SQR_MUL2_C(a3, a2);

        }

        // the s*n terms.  No more doubling past here.
        //for (j = 0; j < NBLOCKS - 1 - i; j++)
        j = 0;
        {
            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c15, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c14, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c13, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c12, b3, b4, b5, b6);
        }
        j = 1;
        {
            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c11, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c10, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c09, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c08, b3, b4, b5, b6);
        }
        j = 2;
        {
            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c07, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c06, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c05, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c04, b3, b4, b5, b6);
        }
        j = 3;

        // finish each triangluar shaped column sum (s * n)
        a1 = c01;
        a2 = c02;
        a3 = c03;

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        VEC_SQR_MUL3_A(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);
        SUB_BIAS_LO(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);

        // final column accumulation
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 0) * VECLEN, a3);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 1) * VECLEN, a2);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 2) * VECLEN, a1);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 3) * VECLEN, a0);
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
        // i odd, block shape 1.
        {
            
            // always a continuation of the full-block loop, so use the same 
            // loading pattern.  Only now we don't need as many b-terms.
            a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);  // {f, b}
            a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);  // {e, a}
            a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);  // {d, 9}

            VEC_SQR_MUL2_D(a0, a1, a2);

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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally, accumulate the two non-doubled terms.
            VEC_SQR_MUL2_C(a1, a0);
        }

        // the s*n terms.  No more doubling past here.
        //for (j = 0; j < NBLOCKS - 1 - i; j++)
        j = 0;
        {
            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c15, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c14, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c13, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c12, b3, b4, b5, b6);
        }
        j = 1;
        {
            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c11, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c10, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c09, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c08, b3, b4, b5, b6);
        }
        j = 2;

        // finish each triangluar shaped column sum (s * n)
        a1 = c05;
        a2 = c06;
        a3 = c07;

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        VEC_SQR_MUL3_A(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);
        SUB_BIAS_LO(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);

        // final column accumulation
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 0) * VECLEN, a3);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 1) * VECLEN, a2);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 2) * VECLEN, a1);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 3) * VECLEN, a0);
    }

    i = 2;
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


            VEC_SQR_MUL2_B(a2, b2, b3);

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


            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally, accumulate the two non-doubled terms.
            VEC_SQR_MUL2_C(a3, a2);

        }

        // the s*n terms.  No more doubling past here.
        //for (j = 0; j < NBLOCKS - 1 - i; j++)
        j = 0;
        {
            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c15, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c14, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c13, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c12, b3, b4, b5, b6);
        }
        j = 1;

        // finish each triangluar shaped column sum (s * n)
        a1 = c09;
        a2 = c10;
        a3 = c11;

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        VEC_SQR_MUL3_A(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);
        SUB_BIAS_LO(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);

        // final column accumulation
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 0) * VECLEN, a3);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 1) * VECLEN, a2);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 2) * VECLEN, a1);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 3) * VECLEN, a0);
    }

    i = 3;
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
        // i odd, block shape 1.
        {
            
            // always a continuation of the full-block loop, so use the same 
            // loading pattern.  Only now we don't need as many b-terms.
            a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);  // {f, b}
            a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);  // {e, a}
            a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);  // {d, 9}

            VEC_SQR_MUL2_D(a0, a1, a2);

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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally, accumulate the two non-doubled terms.
            VEC_SQR_MUL2_C(a1, a0);
        }

        // the s*n terms.  No more doubling past here.
        j = 0;

        // finish each triangluar shaped column sum (s * n)
        a1 = c13;
        a2 = c14;
        a3 = c15;

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        VEC_SQR_MUL3_A(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);
        SUB_BIAS_LO(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);

        // final column accumulation
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 0) * VECLEN, a3);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 1) * VECLEN, a2);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 2) * VECLEN, a1);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 3) * VECLEN, a0);
    }

    //printf("fixed416: \n");
    //for (i = NWORDS - 1; i >= 0; i--)
    //{
    //    printf("%016lx", c->data[i * VECLEN + 0]);
    //}
    //printf("\n");

    c->size = NWORDS;
    return;
}

void vecsqrmod52_fixed624_bfips(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    // 8x sqr:
    // input 8 bignums in the even lanes of a.
    // output 8 squaremod bignums in the even lanes of c.
    int i, j, k;
    vec_bignum_t* b = a;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t NBLOCKS = mdata->NBLOCKS;

    // needed in loops
    __m512i i0, i1;
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19
    __m512i c00, c01, c02, c03, c04, c05, c06, c07, c08, c09, c10, c11;


    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31


    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i nhatvec_e = _mm512_load_epi64(mdata->vrho);
    __m512i hiword = _mm512_set1_epi64(0x000000000000001);
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry_e = 0;
    __mmask8 scarry2;
    __mmask8 scarry;

    // deal with the sign
    c->size = NWORDS;
    c->signmask = 0;

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

    // first half sqr
    i = 0;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        //for (k = 0, j = i; j > (i + 1) / 2; j--, k++)
        k = 0;
        j = 0;
        {
            // i even
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);

            VEC_SQR_MUL4_C(a0, a1, a2, a3);
            

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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally, accumulate the two non-doubled terms.
            VEC_SQR_MUL2_C(a0, a1);
        }

        // need to remove bias from the s*n loop and the 
        // two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);

        // final monty column accumulation
        CARRYPROP_ACCUM_0(c00);
        CARRYPROP_ACCUM_1(c01, c00);
        CARRYPROP_ACCUM_2(c02, c00, c01);
        CARRYPROP_ACCUM_3(c03, c00, c01, c02);
    }

    i = 1;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        //for (k = 0, j = i; j > (i + 1) / 2; j--, k++)
        k = 0;
        j = 1;
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

            VEC_SQR_MUL4_A(a0, a1, a2, b2, b3);
            VEC_SQR_MUL4_A(a1, a2, a3, b3, b4);
            VEC_SQR_MUL4_B(a2, a3, b4, b5, b6);
            
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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally, accumulate the two non-doubled terms.
            VEC_SQR_MUL2(a0, a1);
            
        }
        
        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        j = 0;
        {
            // accumulate s * n
            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c03, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c02, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c01, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c00, b3, b4, b5, b6);
        }

#ifndef IFMA
        // need to remove bias from the s*n loop and the 
        // two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);
#endif

        // final monty column accumulation
        CARRYPROP_ACCUM_0(c04);
        CARRYPROP_ACCUM_1(c05, c04);
        CARRYPROP_ACCUM_2(c06, c04, c05);
        CARRYPROP_ACCUM_3(c07, c04, c05, c06);
    }

    i = 2;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        //for (k = 0, j = i; j > (i + 1) / 2; j--, k++)
        j = 2;
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

        k = 1;
        j = 1;
        {
            // i even
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);

            VEC_SQR_MUL4_C(a0, a1, a2, a3);

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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally, accumulate the two non-doubled terms.
            VEC_SQR_MUL2_C(a0, a1);
        }

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        j = 0;
        {
            // accumulate s * n
            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c03, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c02, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c01, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c00, b3, b4, b5, b6);
        }

        j = 1;
        {
            // accumulate s * n
            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c07, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c06, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c05, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c04, b3, b4, b5, b6);
        }

        // need to remove bias from the s*n loop and the 
        // two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);

        // final monty column accumulation
        CARRYPROP_ACCUM_0(c08);
        CARRYPROP_ACCUM_1(c09, c08);
        CARRYPROP_ACCUM_2(c10, c08, c09);
        CARRYPROP_ACCUM_3(c11, c08, c09, c10);
    }

    // second half sqr
    i = 0;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        //for (k = 0, j = 0; j < (NBLOCKS - i - 1) / 2; j++, k++)
        j = 0;
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

        // NBLOCKS is odd
        j = 1;
        k = 1;
        {
            // i even, block shape 1.
            // always a continuation of the full-block loop, so use the same 
            // loading pattern.  Only now we don't need as many b-terms.
            a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
            a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
            a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);

            VEC_SQR_MUL2_D(a0, a1, a2);
            
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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally the two non-doubled terms.
            VEC_SQR_MUL2_C(a1, a0);
        }


        // the s*n terms.  No more doubling past here.
        //for (j = 0; j < NBLOCKS - 1 - i; j++)
        j = 0;
        {
            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c11, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c10, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c09, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c08, b3, b4, b5, b6);
        }

        j = 1;
        {
            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c07, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c06, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c05, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c04, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum (s * n)
        a1 = c01;
        a2 = c02;
        a3 = c03;

        j = 2;
        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);
        
        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);
        SUB_BIAS_LO(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);

        // final column accumulation
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 0) * VECLEN, a3);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 1) * VECLEN, a2);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 2) * VECLEN, a1);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 3) * VECLEN, a0);

    }

    i = 1;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        // NBLOCKS is odd
        //for (k = 0, j = 0; j < (NBLOCKS - i - 1) / 2; j++, k++)
        // i odd, block shape 2.
        k = 0;
        j = 0;
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

            VEC_SQR_MUL4_D(a0, a1, b0, b1, b2, te0, te1, te2, te3);
            VEC_SQR_MUL2_B(a2, b2, b3);
            VEC_SQR_MUL4_D(a0, a1, b2, b3, b4, te4, te5, te6, te7);

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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally the two non-doubled terms.
            VEC_SQR_MUL2_C(a3, a2);
            
        }
        

        // the s*n terms.  No more doubling past here.
        j = 0;
        {
            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c11, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c10, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c09, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c08, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum (s * n)
        a1 = c05;
        a2 = c06;
        a3 = c07;

        j = 1;
        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);
        SUB_BIAS_LO(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);

        // final column accumulation
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 0) * VECLEN, a3);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 1) * VECLEN, a2);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 2) * VECLEN, a1);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 3) * VECLEN, a0);

    }

    i = 2;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        // The final block shape depends on the parity of i and NBLOCKS
        // Block shape 1 (small upper triangle) if 'i' is odd and NBLOCKS is even,
        // or if 'i' is even and NBLOCKS is odd.
        // Block shape 2 if 'i' is odd and NBLOCKS is odd or if 'i' is even
        // and NBLOCKS is even.
        // NBLOCKS is odd

        //for (k = 0, j = 0; j < (NBLOCKS - i - 1) / 2; j++, k++)
        k = 0;
        j = 0;
        {
            // i even, block shape 1.
            // always a continuation of the full-block loop, so use the same 
            // loading pattern.  Only now we don't need as many b-terms.
            a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
            a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
            a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);

            VEC_SQR_MUL2_D(a0, a1, a2);

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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally the two non-doubled terms.
            VEC_SQR_MUL2_C(a1, a0);
        }


        // finish each triangluar shaped column sum (s * n)
        a1 = c09;
        a2 = c10;
        a3 = c11;

        j = 0;
        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        VEC_SQR_MUL3_C(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);
        SUB_BIAS_LO(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);

        // final column accumulation
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 0) * VECLEN, a3);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 1) * VECLEN, a2);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 2) * VECLEN, a1);
        _mm512_store_epi64(c->data + (i * BLOCKWORDS + 3) * VECLEN, a0);

    }

    c->size = NWORDS;
    return;
}

void vecsqrmod52_fixed416_bfips(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    // 8x sqr:
    // input 8 bignums in the even lanes of a.
    // output 8 squaremod bignums in the even lanes of c.
    int i, j, k;
    vec_bignum_t* b = a;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t NBLOCKS = mdata->NBLOCKS;

    // needed in loops
    __m512i i0, i1;
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19
    __m512i c00, c01, c02, c03, c04, c05, c06, c07;

    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i nhatvec_e = _mm512_load_epi64(mdata->vrho);
    __m512i hiword = _mm512_set1_epi64(0x000000000000001);
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry_e = 0;
    __mmask8 scarry2;
    __mmask8 scarry;

    // deal with the sign
    c->size = NWORDS;
    c->signmask = 0;

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

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

            VEC_SQR_MUL4_C(a0, a1, a2, a3);

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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally, accumulate the two non-doubled terms.
            VEC_SQR_MUL2_C(a0, a1);
        }

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        //for (j = 0; j < i; j++)

        // need to remove bias from the s*n loop and the 
        // two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);

        // final monty column accumulation
        CARRYPROP_ACCUM_0(c00);
        CARRYPROP_ACCUM_1(c01, c00);
        CARRYPROP_ACCUM_2(c02, c00, c01);
        CARRYPROP_ACCUM_3(c03, c00, c01, c02);
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

            VEC_SQR_MUL4_A(a0, a1, a2, b2, b3);
            VEC_SQR_MUL4_A(a1, a2, a3, b3, b4);
            VEC_SQR_MUL4_B(a2, a3, b4, b5, b6);

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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally, accumulate the two non-doubled terms.
            VEC_SQR_MUL2(a0, a1);
        }
        

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        //for (j = 0; j < i; j++)
        j = 0;
        {
            // accumulate s * n
            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c03, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c02, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c01, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c00, b3, b4, b5, b6);
        }
        j = 1;

        // need to remove bias from the s*n loop and the 
        // two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);

        // final monty column accumulation
        CARRYPROP_ACCUM_0(c04);
        CARRYPROP_ACCUM_1(c05, c04);
        CARRYPROP_ACCUM_2(c06, c04, c05);
        CARRYPROP_ACCUM_3(c07, c04, c05, c06);
    }

    //printf("fixed416: \n");
    //for (i = NWORDS - 1; i >= 0; i--)
    //{
    //    printf("%016lx", s->data[i * VECLEN + 0]);
    //}
    //printf("\n");

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


            VEC_SQR_MUL2_B(a2, b2, b3);

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


            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally, accumulate the two non-doubled terms.
            VEC_SQR_MUL2_C(a3, a2);

        }

        // the s*n terms.  No more doubling past here.
        //for (j = 0; j < NBLOCKS - 1 - i; j++)
        j = 0;
        {
            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(c07, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(c06, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(c05, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(c04, b3, b4, b5, b6);
        }
        j = 1;

        // finish each triangluar shaped column sum (s * n)
        a1 = c01;
        a2 = c02;
        a3 = c03;

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        VEC_SQR_MUL3_A(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);

        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);
        SUB_BIAS_LO(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);

        // final column accumulation
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        //_mm512_store_epi64(c->data + (i * BLOCKWORDS + 0) * VECLEN, a3);
        //_mm512_store_epi64(c->data + (i * BLOCKWORDS + 1) * VECLEN, a2);
        //_mm512_store_epi64(c->data + (i * BLOCKWORDS + 2) * VECLEN, a1);
        //_mm512_store_epi64(c->data + (i * BLOCKWORDS + 3) * VECLEN, a0);
        c00 = a3;
        c01 = a2;
        c02 = a1;
        c03 = a0;
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
            // always a continuation of the full-block loop, so use the same 
            // loading pattern.  Only now we don't need as many b-terms.
            a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);  // {f, b}
            a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);  // {e, a}
            a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);  // {d, 9}

            VEC_SQR_MUL2_D(a0, a1, a2);

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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

            // finally, accumulate the two non-doubled terms.
            VEC_SQR_MUL2_C(a1, a0);
        }

        // the s*n terms.  No more doubling past here.
        j = 0;

        // finish each triangluar shaped column sum (s * n)
        a1 = c05;
        a2 = c06;
        a3 = c07;

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        VEC_SQR_MUL3_A(a1, a2, b0, b1);
        VEC_SQR_MUL3_B(a3, b2, b1, b0);
        
        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);
        SUB_BIAS_LO(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);

        // final column accumulation
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        //_mm512_store_epi64(c->data + (i * BLOCKWORDS + 0) * VECLEN, a3);
        //_mm512_store_epi64(c->data + (i * BLOCKWORDS + 1) * VECLEN, a2);
        //_mm512_store_epi64(c->data + (i * BLOCKWORDS + 2) * VECLEN, a1);
        //_mm512_store_epi64(c->data + (i * BLOCKWORDS + 3) * VECLEN, a0);
        c04 = a3;
        c05 = a2;
        c06 = a1;
        c07 = a0;
    }

    //printf("fixed416: \n");
    //for (i = NWORDS - 1; i >= 0; i--)
    //{
    //    printf("%016lx", c->data[i * VECLEN + 0]);
    //}
    //printf("\n");

#ifdef USE_AMM
    _mm512_store_epi64(c->data + 0 * VECLEN, c00);
    _mm512_store_epi64(c->data + 1 * VECLEN, c01);
    _mm512_store_epi64(c->data + 2 * VECLEN, c02);
    _mm512_store_epi64(c->data + 3 * VECLEN, c03);
    _mm512_store_epi64(c->data + 4 * VECLEN, c04);
    _mm512_store_epi64(c->data + 5 * VECLEN, c05);
    _mm512_store_epi64(c->data + 6 * VECLEN, c06);
    _mm512_store_epi64(c->data + 7 * VECLEN, c07);
    _mm512_store_epi64(c->data + NWORDS * VECLEN, zero);
#else

    __m512i cvec;
    __m512i nvec;
    __m512i bvec;

    // compare
    scarry = 0;	    // sub mask
    scarry2 = 0;	// keep looking mask

    cvec = c07;
    nvec = _mm512_load_epi64(mdata->n->data + 7 * VECLEN);
    // compare those that have not already been decided using the mask
    scarry |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_GT);
    scarry2 |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_LT);

    // decided all of them, stop comparing.
    if (scarry2 == 0xff) goto sub;

    cvec = c06;
    nvec = _mm512_load_epi64(mdata->n->data + 6 * VECLEN);
    // compare those that have not already been decided using the mask
    scarry |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_GT);
    scarry2 |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_LT);

    // decided all of them, stop comparing.
    if (scarry2 == 0xff) goto sub;

    cvec = c05;
    nvec = _mm512_load_epi64(mdata->n->data + 5 * VECLEN);
    // compare those that have not already been decided using the mask
    scarry |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_GT);
    scarry2 |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_LT);

    // decided all of them, stop comparing.
    if (scarry2 == 0xff) goto sub;

    cvec = c04;
    nvec = _mm512_load_epi64(mdata->n->data + 4 * VECLEN);
    // compare those that have not already been decided using the mask
    scarry |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_GT);
    scarry2 |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_LT);

    // decided all of them, stop comparing.
    if (scarry2 == 0xff) goto sub;

    cvec = c03;
    nvec = _mm512_load_epi64(mdata->n->data + 3 * VECLEN);
    // compare those that have not already been decided using the mask
    scarry |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_GT);
    scarry2 |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_LT);

    // decided all of them, stop comparing.
    if (scarry2 == 0xff) goto sub;

    cvec = c02;
    nvec = _mm512_load_epi64(mdata->n->data + 2 * VECLEN);
    // compare those that have not already been decided using the mask
    scarry |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_GT);
    scarry2 |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_LT);

    // decided all of them, stop comparing.
    if (scarry2 == 0xff) goto sub;

    cvec = c01;
    nvec = _mm512_load_epi64(mdata->n->data + 1 * VECLEN);
    // compare those that have not already been decided using the mask
    scarry |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_GT);
    scarry2 |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_LT);

    // decided all of them, stop comparing.
    if (scarry2 == 0xff) goto sub;

    cvec = c00;
    nvec = _mm512_load_epi64(mdata->n->data + 0 * VECLEN);
    // compare those that have not already been decided using the mask
    scarry |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_GT);
    scarry2 |= _mm512_mask_cmp_epu64_mask(~scarry2, cvec, nvec, _MM_CMPINT_LT);

    // decided all of them, stop comparing.
    if (scarry2 == 0xff) goto sub;

    // check for equal as well by flipping mask bits that have still
    // not been decided (i.e., are equal)
    scarry |= (~scarry2);

sub:

    if (scarry == 0) goto done;

    // subtract n from c when c is not less than n, as indicated by a 1 bit in mask
    scarry2 = 0;

    nvec = _mm512_load_epi64(mdata->n->data + 0 * VECLEN);
    bvec = _mm512_mask_sbb_src_epi52(zero, c00, scarry, scarry2, nvec, &scarry2);
    _mm512_store_epi64(c->data + 0 * VECLEN, bvec);

    nvec = _mm512_load_epi64(mdata->n->data + 1 * VECLEN);
    bvec = _mm512_mask_sbb_src_epi52(zero, c01, scarry, scarry2, nvec, &scarry2);
    _mm512_store_epi64(c->data + 1 * VECLEN, bvec);

    nvec = _mm512_load_epi64(mdata->n->data + 2 * VECLEN);
    bvec = _mm512_mask_sbb_src_epi52(zero, c02, scarry, scarry2, nvec, &scarry2);
    _mm512_store_epi64(c->data + 2 * VECLEN, bvec);

    nvec = _mm512_load_epi64(mdata->n->data + 3 * VECLEN);
    bvec = _mm512_mask_sbb_src_epi52(zero, c03, scarry, scarry2, nvec, &scarry2);
    _mm512_store_epi64(c->data + 3 * VECLEN, bvec);

    nvec = _mm512_load_epi64(mdata->n->data + 4 * VECLEN);
    bvec = _mm512_mask_sbb_src_epi52(zero, c04, scarry, scarry2, nvec, &scarry2);
    _mm512_store_epi64(c->data + 4 * VECLEN, bvec);

    nvec = _mm512_load_epi64(mdata->n->data + 5 * VECLEN);
    bvec = _mm512_mask_sbb_src_epi52(zero, c05, scarry, scarry2, nvec, &scarry2);
    _mm512_store_epi64(c->data + 5 * VECLEN, bvec);

    nvec = _mm512_load_epi64(mdata->n->data + 6 * VECLEN);
    bvec = _mm512_mask_sbb_src_epi52(zero, c06, scarry, scarry2, nvec, &scarry2);
    _mm512_store_epi64(c->data + 6 * VECLEN, bvec);

    nvec = _mm512_load_epi64(mdata->n->data + 7 * VECLEN);
    bvec = _mm512_mask_sbb_src_epi52(zero, c07, scarry, scarry2, nvec, &scarry2);
    _mm512_store_epi64(c->data + 7 * VECLEN, bvec);

done:

#endif

    c->size = NWORDS;
    return;
}

void vecsqrmod52_mul(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    vecmulmod52(a, a, c, n, s, mdata);
    return;
}

void vecsqrmod52(vec_bignum_t *a, vec_bignum_t *c, vec_bignum_t *n, vec_bignum_t *s, vec_monty_t*mdata)
{
    // 8x sqr:
    // input 8 bignums in the even lanes of a.
    // output 8 squaremod bignums in the even lanes of c.
    int i, j, k;
    vec_bignum_t *b = a;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t NBLOCKS = mdata->NBLOCKS;

    // needed in loops
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
    __m512i nhatvec_e = _mm512_load_epi64(mdata->vrho);
    __m512i hiword = _mm512_set1_epi64(0x000000000000001);
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry_e = 0;
    __mmask8 scarry2;
    __mmask8 scarry;

#ifdef USE_AMM
    uint64_t* outdata = c->data;
#else
    uint64_t* outdata = s->data;
#endif

    // deal with the sign
    c->size = NWORDS;
    c->signmask = 0;

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

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

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        for (j = 0; j < i; j++)
        {
            // accumulate s * n
            a0 = _mm512_load_epi32(outdata + ((j + 1) * BLOCKWORDS - 1) * VECLEN);
            a1 = _mm512_load_epi32(outdata + ((j + 1) * BLOCKWORDS - 2) * VECLEN);
            a2 = _mm512_load_epi32(outdata + ((j + 1) * BLOCKWORDS - 3) * VECLEN);
            a3 = _mm512_load_epi32(outdata + ((j + 1) * BLOCKWORDS - 4) * VECLEN);

            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

#ifndef IFMA
        // need to remove bias from the s*n loop and the 
        // two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);
#endif

        // final monty column accumulation
        {
            // now, column by column, add in the s*n contribution and reduce to 
            // a single 64+x bit accumulator while storing the intermediate product
            // 's' as we go.
            j = 0;
            // accumulate this column-sum
            // now carry propagate low to high.
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            _mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);
            b0 = _mm512_load_epi64(n->data + 0 * VECLEN);

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            for (k = 0; k < j; k++)
            {
                a0 = _mm512_load_epi64(outdata + (i * BLOCKWORDS + k) * VECLEN);
                b1 = _mm512_load_epi64(n->data + (j - k) * VECLEN);

                VEC_MUL_ACCUM_LOHI_PD(a0, b1, acc_e0, acc_e1);
            }

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            _mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            for (k = 0; k < j; k++)
            {
                a0 = _mm512_load_epi64(outdata + (i * BLOCKWORDS + k) * VECLEN);
                b1 = _mm512_load_epi64(n->data + (j - k) * VECLEN);

                VEC_MUL_ACCUM_LOHI_PD(a0, b1, acc_e0, acc_e1);
            }

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            _mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            for (k = 0; k < j; k++)
            {
                a0 = _mm512_load_epi64(outdata + (i * BLOCKWORDS + k) * VECLEN);
                b1 = _mm512_load_epi64(n->data + (j - k) * VECLEN);

                VEC_MUL_ACCUM_LOHI_PD(a0, b1, acc_e0, acc_e1);
            }

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            _mm512_store_epi64(outdata + (i * BLOCKWORDS + j) * VECLEN, a0);

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;
        }
    }

    //printf("normal624: \n");
    //for (i = NWORDS - 1; i >= 0; i--)
    //{
    //    printf("%016lx", s->data[i * VECLEN + 0]);
    //}
    //printf("\n");

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

                b0 = _mm512_load_epi64(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 1) * VECLEN);
                b1 = _mm512_load_epi64(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 2) * VECLEN);
                b2 = _mm512_load_epi64(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 3) * VECLEN);
                b3 = _mm512_load_epi64(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 4) * VECLEN);
                b4 = _mm512_load_epi64(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 5) * VECLEN);

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

                b0 = _mm512_load_epi64(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 1) * VECLEN); // {9, 5}
                b1 = _mm512_load_epi64(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 2) * VECLEN);	// {a, 6}
                b2 = _mm512_load_epi64(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 3) * VECLEN);	// {b, 7}
                b3 = _mm512_load_epi64(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 4) * VECLEN);	// {c, 8}
                b4 = _mm512_load_epi64(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 5) * VECLEN); // {d, 9}

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

        // the s*n terms.  No more doubling past here.
        for (j = 0; j < NBLOCKS - 1 - i; j++)
        {
            a0 = _mm512_load_epi64(outdata + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
            a1 = _mm512_load_epi64(outdata + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
            a2 = _mm512_load_epi64(outdata + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);
            a3 = _mm512_load_epi64(outdata + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);

            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum (s * n)
        a1 = _mm512_load_epi64(outdata + (i * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(outdata + (i * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(outdata + (i * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a1, b0, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b1, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b0, te2, te3);
#else
        // finish each triangluar shaped column sum (s * n)
        // a1*b0 -> t0/1
        // a2*b1 -> t0/1
        // a2*b0 -> t2/3
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

        //VEC_MUL_ACCUM_LOHI(a3, b2, te0, te1);
        //VEC_MUL_ACCUM_LOHI(a3, b1, te2, te3);
        //VEC_MUL_ACCUM_LOHI(a3, b0, te4, te5);
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
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);
        SUB_BIAS_LO(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);

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

            _mm512_store_epi64(outdata + (i * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(outdata + (i * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(outdata + (i * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(outdata + (i * BLOCKWORDS + 3) * VECLEN, a0);
        }


    }

#ifdef USE_AMM
    //for (i = NWORDS - 1; i >= 0; i--)
    //{
    //    b0 = _mm512_load_epi64(s->data + i * VECLEN);
    //    _mm512_store_epi64(c->data + i * VECLEN, b0);
    //}
#else
    a0 = acc_e0;
    scarry2 = _mm512_cmp_epu64_mask(a0, zero, _MM_CMPINT_EQ);

    // subtract n from tmp
    scarry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        a1 = _mm512_load_epi64(s->data + i * VECLEN);
        b0 = _mm512_load_epi64(n->data + i * VECLEN);
        a0 = _mm512_sbb_epi52(a1, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
    }

    // negate any final borrows if there was also a final carry.
    scarry &= scarry2;

    // if there was a final borrow, we didn't need to do the subtraction after all.
    // replace with original results based on final borrow mask.
    for (i = NWORDS - 1; i >= 0; i--)
    {
        b0 = _mm512_load_epi64(s->data + i * VECLEN);
        _mm512_mask_store_epi64(c->data + i * VECLEN, scarry, b0);
    }
#endif

    //printf("normal624: \n");
    //for (i = NWORDS - 1; i >= 0; i--)
    //{
    //    printf("%016lx", c->data[i * VECLEN + 0]);
    //}
    //printf("\n");
    //exit(0);

    c->size = NWORDS;
    return;
}

void vecsqrmod52_avxecm(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    // 8x sqr:
    // input 8 bignums in the even lanes of a.
    // output 8 squaremod bignums in the even lanes of c.
    int i, j, k;
    vec_bignum_t* b = a;
    int NWORDS = mdata->NWORDS;
    int NBLOCKS = mdata->NBLOCKS;

    // needed in loops
    __m512i i0, i1;
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19
    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i nhatvec_e = _mm512_load_epi64(mdata->vrho);
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry2;
    __mmask8 scarry;

    // deal with the sign
    c->size = NWORDS;
    c->signmask = 0;

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

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

        }
        else
        {
            // i even
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);

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

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

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
        }

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        for (j = 0; j < i; j++)
        {
            // accumulate s * n
            a0 = _mm512_load_epi32(s->data + ((j + 1) * BLOCKWORDS - 1) * VECLEN);
            a1 = _mm512_load_epi32(s->data + ((j + 1) * BLOCKWORDS - 2) * VECLEN);
            a2 = _mm512_load_epi32(s->data + ((j + 1) * BLOCKWORDS - 3) * VECLEN);
            a3 = _mm512_load_epi32(s->data + ((j + 1) * BLOCKWORDS - 4) * VECLEN);

            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // need to remove bias from the s*n loop and the 
        // two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 0,
            i * 4 + 1,
            i * 4 + 0);

        // final monty column accumulation
        {
            // now, column by column, add in the s*n contribution and reduce to 
            // a single 64+x bit accumulator while storing the intermediate product
            // 's' as we go.
            j = 0;
            // accumulate this column-sum
            // now carry propagate low to high.
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);
            b0 = _mm512_load_epi64(n->data + 0 * VECLEN);

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            for (k = 0; k < j; k++)
            {
                a0 = _mm512_load_epi64(s->data + (i * BLOCKWORDS + k) * VECLEN);
                b1 = _mm512_load_epi64(n->data + (j - k) * VECLEN);

                VEC_MUL_ACCUM_LOHI_PD(a0, b1, acc_e0, acc_e1);
            }

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            for (k = 0; k < j; k++)
            {
                a0 = _mm512_load_epi64(s->data + (i * BLOCKWORDS + k) * VECLEN);
                b1 = _mm512_load_epi64(n->data + (j - k) * VECLEN);

                VEC_MUL_ACCUM_LOHI_PD(a0, b1, acc_e0, acc_e1);
            }

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

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

            for (k = 0; k < j; k++)
            {
                a0 = _mm512_load_epi64(s->data + (i * BLOCKWORDS + k) * VECLEN);
                b1 = _mm512_load_epi64(n->data + (j - k) * VECLEN);

                VEC_MUL_ACCUM_LOHI_PD(a0, b1, acc_e0, acc_e1);
            }

            //a0 = _mm512_and_epi64(vlmask, _mm512_mullo_epi64(nhatvec_e, acc_e0));
            _mm512_mullo_epi52(a0, nhatvec_e, acc_e0);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

            // add in the final product
            VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;
        }
    }

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

                // now double
                te1 = _mm512_slli_epi64(te1, 1);
                te3 = _mm512_slli_epi64(te3, 1);
                te5 = _mm512_slli_epi64(te5, 1);
                te7 = _mm512_slli_epi64(te7, 1);
                te0 = _mm512_slli_epi64(te0, 1);
                te2 = _mm512_slli_epi64(te2, 1);
                te4 = _mm512_slli_epi64(te4, 1);
                te6 = _mm512_slli_epi64(te6, 1);

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
            }
            else
            {
                // i even, block shape 1.
                // always a continuation of the full-block loop, so use the same 
                // loading pattern.  Only now we don't need as many b-terms.
                a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
                a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
                a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);

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

                // now double
                te1 = _mm512_slli_epi64(te1, 1);
                te3 = _mm512_slli_epi64(te3, 1);
                te5 = _mm512_slli_epi64(te5, 1);
                te7 = _mm512_slli_epi64(te7, 1);
                te0 = _mm512_slli_epi64(te0, 1);
                te2 = _mm512_slli_epi64(te2, 1);
                te4 = _mm512_slli_epi64(te4, 1);
                te6 = _mm512_slli_epi64(te6, 1);

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

                // now double
                te1 = _mm512_slli_epi64(te1, 1);
                te3 = _mm512_slli_epi64(te3, 1);
                te5 = _mm512_slli_epi64(te5, 1);
                te7 = _mm512_slli_epi64(te7, 1);
                te0 = _mm512_slli_epi64(te0, 1);
                te2 = _mm512_slli_epi64(te2, 1);
                te4 = _mm512_slli_epi64(te4, 1);
                te6 = _mm512_slli_epi64(te6, 1);

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

                // now double
                te1 = _mm512_slli_epi64(te1, 1);
                te3 = _mm512_slli_epi64(te3, 1);
                te5 = _mm512_slli_epi64(te5, 1);
                te7 = _mm512_slli_epi64(te7, 1);
                te0 = _mm512_slli_epi64(te0, 1);
                te2 = _mm512_slli_epi64(te2, 1);
                te4 = _mm512_slli_epi64(te4, 1);
                te6 = _mm512_slli_epi64(te6, 1);

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
            }
        }

        // the s*n terms.  No more doubling past here.
        for (j = 0; j < NBLOCKS - 1 - i; j++)
        {
            a0 = _mm512_load_epi64(s->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
            a1 = _mm512_load_epi64(s->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
            a2 = _mm512_load_epi64(s->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);
            a3 = _mm512_load_epi64(s->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);

            b0 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum (s * n)
        a1 = _mm512_load_epi64(s->data + (i * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(s->data + (i * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(s->data + (i * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(n->data + (NWORDS - 3) * VECLEN);

        // finish each triangluar shaped column sum (s * n)
        // a1*b0 -> t0/1
        // a2*b1 -> t0/1
        // a2*b0 -> t2/3
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
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);
        SUB_BIAS_LO(
            j * 4 + 4,
            j * 4 + 2,
            j * 4 + 2,
            j * 4 + 0);

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

    a0 = acc_e0;
    scarry2 = _mm512_cmp_epu64_mask(a0, zero, _MM_CMPINT_EQ);

    // subtract n from tmp
    scarry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        a1 = _mm512_load_epi64(s->data + i * VECLEN);
        b0 = _mm512_load_epi64(n->data + i * VECLEN);
        a0 = _mm512_sbb_epi52(a1, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
    }

    // negate any final borrows if there was also a final carry.
    scarry &= scarry2;

    // if there was a final borrow, we didn't need to do the subtraction after all.
    // replace with original results based on final borrow mask.
    for (i = NWORDS - 1; i >= 0; i--)
    {
        b0 = _mm512_load_epi64(s->data + i * VECLEN);
        _mm512_mask_store_epi64(c->data + i * VECLEN, scarry, b0);
    }

    c->size = NWORDS;
    return;
}

void vecaddmod52_avxecm(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata)
{
    // assumptions:
    // a, b, c are of length VECLEN * NWORDS
    // a, b, c, and n are aligned
    // a and b are both positive
    // n is the montgomery base
    int i;
    int NWORDS = mdata->NWORDS;
    __mmask8 carry = 0;
    __mmask8 mask = 0;
    __mmask8 mask2 = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    __m512i nvec;

    // add
    for (i = 0; i < NWORDS; i++)
    {
        avec = _mm512_load_epi64(a->data + i * VECLEN);
        bvec = _mm512_load_epi64(b->data + i * VECLEN);
        cvec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c->data + i * VECLEN, cvec);
    }

    mask = carry;	// sub mask
    mask2 = mask;	// keep looking mask

    // compare.  the initial mask is equal to the addition carry
    // because if the most significant word has a carry, then the
    // result is bigger than n.
    for (i = NWORDS - 1; i >= 0; i--)
    {
        cvec = _mm512_load_epi64(c->data + i * VECLEN);
        nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
        // compare those that have not already been decided using the mask
        mask |= _mm512_mask_cmp_epu64_mask(~mask2, cvec, nvec, _MM_CMPINT_GT);
        mask2 |= _mm512_mask_cmp_epu64_mask(~mask2, cvec, nvec, _MM_CMPINT_LT);

        // decided all of them, stop comparing.
        if (mask2 == 0xff)
        {
            break;
        }
    }

    // check for equal as well by flipping mask bits that have still
    // not been decided (i.e., are equal)
    mask |= (~mask2);

    // subtract n from c when c is not less than n, as indicated by a 1 bit in mask
    carry = 0;
    for (i = 0; (i < NWORDS) && (mask > 0); i++)
    {
        cvec = _mm512_load_epi64(c->data + i * VECLEN);
        nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
        bvec = _mm512_mask_sbb_epi52(cvec, mask, carry, nvec, &carry);
        _mm512_store_epi64(c->data + i * VECLEN, bvec);
    }
    return;
}

void vecaddmod52(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, vec_monty_t* mdata)
{
    // assumptions:
    // a, b, c are of length VECLEN * NWORDS
    // a, b, c, and n are aligned
    // a and b are both positive
    // n is the montgomery base
    int i;
    uint32_t NWORDS = mdata->NWORDS; // a->WORDS_ALLOC;
    __mmask8 carry = 0;
    __mmask8 carryout = 0;
    __mmask8 mask = 0;
    __mmask8 mask2 = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    __m512i nvec;
    __m512i lomask = _mm512_set1_epi64(0xfffffffffffffULL);

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

    mask = carry;	// sub mask
    mask2 = mask;	// keep looking mask

#ifndef USE_AMM
    // compare.  the initial mask is equal to the addition carry
    // because if the most significant word has a carry, then the
    // result is bigger than n.
    for (i = NWORDS - 1; i >= 0; i--)
    {
        cvec = _mm512_load_epi64(c->data + i * VECLEN);
        nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
        // compare those that have not already been decided using the mask
        mask |= _mm512_mask_cmp_epu64_mask(~mask2, cvec, nvec, _MM_CMPINT_GT);
        mask2 |= _mm512_mask_cmp_epu64_mask(~mask2, cvec, nvec, _MM_CMPINT_LT);

        // decided all of them, stop comparing.
        if (mask2 == 0xff)
        {
            break;
        }
    }

    // check for equal as well by flipping mask bits that have still
    // not been decided (i.e., are equal)
    mask |= (~mask2);
#endif

    // subtract n from c when c is not less than n, as indicated by a 1 bit in mask
    carry = 0;
    for (i = 0; (i < NWORDS) && (mask > 0); i++)
    {
        cvec = _mm512_load_epi64(c->data + i * VECLEN);
        nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
        
        //bvec = _mm512_mask_sbb_epi52(cvec, mask, carry, nvec, &carry);
        __m512i t = _mm512_mask_sub_epi64(cvec, mask, cvec, nvec);
        carryout = _mm512_mask_cmpgt_epu64_mask(mask, nvec, cvec);
        __m512i t2 = _mm512_mask_sub_epi64(cvec, mask, t, _mm512_maskz_set1_epi64(carry, 1));
        carry = _mm512_kor(carryout, _mm512_cmpgt_epu64_mask(t2, t));
        bvec = _mm512_and_epi64(t2, lomask);

        _mm512_store_epi64(c->data + i * VECLEN, bvec);
    }
    return;
}

void vecsubmod52_avxecm(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata)
{
    // assumptions:
    // a, b, c are of length VECLEN * NWORDS
    // s1 is of length VECLEN
    // a, b, c, n, and s1 are aligned
    // a and b are both positive
    // a >= b
    // n is the montgomery base
    int i;
    int NWORDS = mdata->NWORDS;

    __mmask8 carry = 0;
    __mmask8 mask = 0;
    __m512i nvec;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;

    // subtract
    carry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        avec = _mm512_load_epi64(a->data + i * VECLEN);
        bvec = _mm512_load_epi64(b->data + i * VECLEN);
        cvec = _mm512_sbb_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c->data + i * VECLEN, cvec);
    }

    // if we had a final carry, then b was bigger than a so we need to re-add n.
    mask = carry;
    carry = 0;
    for (i = 0; (i < NWORDS) && (mask > 0); i++)
    {
        avec = _mm512_load_epi64(c->data + i * VECLEN);
        nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
        cvec = _mm512_mask_adc_epi52(avec, mask, carry, nvec, &carry);
        _mm512_store_epi64(c->data + i * VECLEN, cvec);
    }
    return;
}

void vecsubmod52(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, vec_monty_t* mdata)
{
    // assumptions:
    // a, b, c are of length VECLEN * NWORDS
    // s1 is of length VECLEN
    // a, b, c, n, and s1 are aligned
    // a and b are both positive
    // a >= b
    // n is the montgomery base
    int i;
    uint32_t NWORDS = mdata->NWORDS; // a->WORDS_ALLOC;
    __mmask8 carry = 0;
    __mmask8 carryout = 0;
    __mmask8 mask = 0;
    __mmask8 mask2 = 0;
    __m512i nvec;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    __m512i lomask = _mm512_set1_epi64(0xfffffffffffffULL);

    // subtract
    carry = 0;
    for (i = 0; i < NWORDS; i++)
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

    // if we had a final carry, then b was bigger than a so we need to re-add n.
    mask = carry;
    carry = 0;
    for (i = 0; (i < NWORDS) && (mask > 0); i++)
    {
        avec = _mm512_load_epi64(c->data + i * VECLEN);
        nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
        //cvec = _mm512_mask_adc_epi52(avec, mask, carry, nvec, &carry);

        __m512i t = _mm512_mask_add_epi64(avec, mask, avec, nvec);
        t = _mm512_mask_add_epi64(avec, mask, t, _mm512_maskz_set1_epi64(carry, 1));
        carry = _mm512_mask_cmpgt_epu64_mask(mask, t, lomask);
        cvec = _mm512_and_epi64(t, lomask);

        _mm512_store_epi64(c->data + i * VECLEN, cvec);
    }
    return;
}

void vecsignedaddmod52(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, vec_bignum_t *n)
{
    // assumptions:
    // a, b, c are of length VECLEN * NWORDS
    // a, b, c, and n are aligned
    // a and b are both positive
    // n is the montgomery base
    int i;
    uint32_t NWORDS = a->WORDS_ALLOC;
    __mmask8 borrow = 0;
    __mmask8 carry = 0;
    __mmask8 mask = 0;
    __mmask8 mask2 = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    __m512i nvec;

    // this is nominally an add, but if a or b are negative then
    // it's a subtract.  Since the vector can contain both cases
    // at the same time we have to perform both operations with 
    // the appropriate masks.
    // if both positive or negative then we add.
    mask = a->signmask ^ b->signmask;

    // add
    for (i = 0; i < NWORDS; i++)
    {
        avec = _mm512_load_epi64(a->data + i * VECLEN);
        bvec = _mm512_load_epi64(b->data + i * VECLEN);
        // this won't work, because the sbb will replace the lanes
        // where mask=1 with avec.
        cvec = _mm512_mask_adc_epi52(avec, mask, carry, bvec, &carry);
        cvec = _mm512_mask_sbb_epi52(avec, ~mask, borrow, bvec, &borrow);
        _mm512_store_epi64(c->data + i * VECLEN, cvec);
    }

    mask = carry;	// sub mask
    mask2 = mask;	// keep looking mask

    // compare.  the initial mask is equal to the addition carry
    // because if the most significant word has a carry, then the
    // result is bigger than n.
    for (i = NWORDS - 1; i >= 0; i--)
    {
        cvec = _mm512_load_epi64(c->data + i * VECLEN);
        nvec = _mm512_load_epi64(n->data + i * VECLEN);
        // compare those that have not already been decided using the mask
        mask |= _mm512_mask_cmp_epu64_mask(~mask2, cvec, nvec, _MM_CMPINT_GT);
        mask2 |= _mm512_mask_cmp_epu64_mask(~mask2, cvec, nvec, _MM_CMPINT_LT);

        // decided all of them, stop comparing.
        if (mask2 == 0xff)
        {
            break;
        }
    }

    // check for equal as well by flipping mask bits that have still
    // not been decided (i.e., are equal)
    mask |= (~mask2);

    // subtract n from c when c is not less than n, as indicated by a 1 bit in mask
    carry = 0;
    for (i = 0; (i < NWORDS) && (mask > 0); i++)
    {
        cvec = _mm512_load_epi64(c->data + i * VECLEN);
        nvec = _mm512_load_epi64(n->data + i * VECLEN);
        bvec = _mm512_mask_sbb_epi52(cvec, mask, carry, nvec, &carry);
        _mm512_store_epi64(c->data + i * VECLEN, bvec);
    }
    return;
}

void vec_simul_addsub52_fixed1040(vec_bignum_t* a, vec_bignum_t* b, 
    vec_bignum_t* sum, vec_bignum_t* diff,
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
    uint32_t NWORDS = mdata->NWORDS; // a->WORDS_ALLOC;
    __mmask8 carry = 0;
    __mmask8 borrow = 0;
    __mmask8 cmask = 0;
    __mmask8 bmask = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    __m512i nvec;
    __m512i lomask = _mm512_set1_epi64(0xfffffffffffff);
    __m512i
        s00, s01, s02, s03, s04, s05, s06, s07,
        s08, s09, s10, s11, s12, s13, s14, s15,
        s16, s17, s18, s19;
    __m512i
        d00, d01, d02, d03, d04, d05, d06, d07,
        d08, d09, d10, d11, d12, d13, d14, d15,
        d16, d17, d18, d19;

    // add
    avec = _mm512_load_epi64(a->data + 0 * VECLEN);
    bvec = _mm512_load_epi64(b->data + 0 * VECLEN);
    //s00 = _mm512_adc_epi52(avec, carry, bvec, &carry);
    //d00 = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
    _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, s00, d00)
    avec = _mm512_load_epi64(a->data + 1 * VECLEN);
    bvec = _mm512_load_epi64(b->data + 1 * VECLEN);
    //s01 = _mm512_adc_epi52(avec, carry, bvec, &carry);
    //d01 = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
    _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, s01, d01)
    avec = _mm512_load_epi64(a->data + 2 * VECLEN);
    bvec = _mm512_load_epi64(b->data + 2 * VECLEN);
    //s02 = _mm512_adc_epi52(avec, carry, bvec, &carry);
    //d02 = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
    _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, s02, d02)
    avec = _mm512_load_epi64(a->data + 3 * VECLEN);
    bvec = _mm512_load_epi64(b->data + 3 * VECLEN);
    //s03 = _mm512_adc_epi52(avec, carry, bvec, &carry);
    //d03 = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
    _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, s03, d03)

    if (NWORDS == 4) goto addsubnormalize;

    avec = _mm512_load_epi64(a->data + 4 * VECLEN);
    bvec = _mm512_load_epi64(b->data + 4 * VECLEN);
    //s04 = _mm512_adc_epi52(avec, carry, bvec, &carry);
    //d04 = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
    _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, s04, d04)
    avec = _mm512_load_epi64(a->data + 5 * VECLEN);
    bvec = _mm512_load_epi64(b->data + 5 * VECLEN);
    //s05 = _mm512_adc_epi52(avec, carry, bvec, &carry);
    //d05 = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
    _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, s05, d05)
    avec = _mm512_load_epi64(a->data + 6 * VECLEN);
    bvec = _mm512_load_epi64(b->data + 6 * VECLEN);
    //s06 = _mm512_adc_epi52(avec, carry, bvec, &carry);
    //d06 = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
    _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, s06, d06)
    avec = _mm512_load_epi64(a->data + 7 * VECLEN);
    bvec = _mm512_load_epi64(b->data + 7 * VECLEN);
    //s07 = _mm512_adc_epi52(avec, carry, bvec, &carry);
    //d07 = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
    _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, s07, d07)

    if (NWORDS == 8) goto addsubnormalize;

    avec = _mm512_load_epi64(a->data + 8 * VECLEN);
    bvec = _mm512_load_epi64(b->data + 8 * VECLEN);
    //s08 = _mm512_adc_epi52(avec, carry, bvec, &carry);
    //d08 = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
    _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, s08, d08)
    avec = _mm512_load_epi64(a->data + 9 * VECLEN);
    bvec = _mm512_load_epi64(b->data + 9 * VECLEN);
    //s09 = _mm512_adc_epi52(avec, carry, bvec, &carry);
    //d09 = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
    _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, s09, d09)
    avec = _mm512_load_epi64(a->data + 10 * VECLEN);
    bvec = _mm512_load_epi64(b->data + 10 * VECLEN);
    //s10 = _mm512_adc_epi52(avec, carry, bvec, &carry);
    //d10 = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
    _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, s10, d10)
    avec = _mm512_load_epi64(a->data + 11 * VECLEN);
    bvec = _mm512_load_epi64(b->data + 11 * VECLEN);
    //s11 = _mm512_adc_epi52(avec, carry, bvec, &carry);
    //d11 = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
    _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, s11, d11)

    if (NWORDS == 12) goto addsubnormalize;

    avec = _mm512_load_epi64(a->data + 12 * VECLEN);
    bvec = _mm512_load_epi64(b->data + 12 * VECLEN);
    //s12 = _mm512_adc_epi52(avec, carry, bvec, &carry);
    //d12 = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
    _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, s12, d12)
    avec = _mm512_load_epi64(a->data + 13 * VECLEN);
    bvec = _mm512_load_epi64(b->data + 13 * VECLEN);
    //s13 = _mm512_adc_epi52(avec, carry, bvec, &carry);
    //d13 = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
    _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, s13, d13)
    avec = _mm512_load_epi64(a->data + 14 * VECLEN);
    bvec = _mm512_load_epi64(b->data + 14 * VECLEN);
    //s14 = _mm512_adc_epi52(avec, carry, bvec, &carry);
    //d14 = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
    _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, s14, d14)
    avec = _mm512_load_epi64(a->data + 15 * VECLEN);
    bvec = _mm512_load_epi64(b->data + 15 * VECLEN);
    //s15 = _mm512_adc_epi52(avec, carry, bvec, &carry);
    //d15 = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
    _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, s15, d15)
    
    if (NWORDS == 16) goto addsubnormalize;

    avec = _mm512_load_epi64(a->data + 16 * VECLEN);
    bvec = _mm512_load_epi64(b->data + 16 * VECLEN);
    //s16 = _mm512_adc_epi52(avec, carry, bvec, &carry);
    //d16 = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
    _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, s16, d16)
    avec = _mm512_load_epi64(a->data + 17 * VECLEN);
    bvec = _mm512_load_epi64(b->data + 17 * VECLEN);
    //s17 = _mm512_adc_epi52(avec, carry, bvec, &carry);
    //d17 = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
    _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, s17, d17)
    avec = _mm512_load_epi64(a->data + 18 * VECLEN);
    bvec = _mm512_load_epi64(b->data + 18 * VECLEN);
    //s18 = _mm512_adc_epi52(avec, carry, bvec, &carry);
    //d18 = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
    _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, s18, d18)
    avec = _mm512_load_epi64(a->data + 19 * VECLEN);
    bvec = _mm512_load_epi64(b->data + 19 * VECLEN);
    //s19 = _mm512_adc_epi52(avec, carry, bvec, &carry);
    //d19 = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
    _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, s19, d19)

addsubnormalize:

    cmask = carry;	    // result too big, need to subtract n
    bmask = borrow;     // result too small, need to add n

    carry = 0;
    borrow = 0;

    i = 0;
    nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
    //avec = _mm512_mask_sbb_epi52(s00, cmask, borrow, nvec, &borrow);
    //bvec = _mm512_mask_adc_epi52(d00, bmask, carry, nvec, &carry);
    _mm512_mask_adcsbb_epi52(d00, s00, bmask, cmask, carry, borrow, nvec, lomask, bvec, avec);
    _mm512_store_epi64(sum->data + i * VECLEN, avec);
    _mm512_store_epi64(diff->data + i * VECLEN, bvec);

    i = 1;
    nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
    //avec = _mm512_mask_sbb_epi52(s01, cmask, borrow, nvec, &borrow);
    //bvec = _mm512_mask_adc_epi52(d01, bmask, carry, nvec, &carry);
    _mm512_mask_adcsbb_epi52(d01, s01, bmask, cmask, carry, borrow, nvec, lomask, bvec, avec);
    _mm512_store_epi64(sum->data + i * VECLEN, avec);
    _mm512_store_epi64(diff->data + i * VECLEN, bvec);

    i = 2;
    nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
    //avec = _mm512_mask_sbb_epi52(s02, cmask, borrow, nvec, &borrow);
    //bvec = _mm512_mask_adc_epi52(d02, bmask, carry, nvec, &carry);
    _mm512_mask_adcsbb_epi52(d02, s02, bmask, cmask, carry, borrow, nvec, lomask, bvec, avec);
    _mm512_store_epi64(sum->data + i * VECLEN, avec);
    _mm512_store_epi64(diff->data + i * VECLEN, bvec);

    i = 3;
    nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
    //avec = _mm512_mask_sbb_epi52(s03, cmask, borrow, nvec, &borrow);
    //bvec = _mm512_mask_adc_epi52(d03, bmask, carry, nvec, &carry);
    _mm512_mask_adcsbb_epi52(d03, s03, bmask, cmask, carry, borrow, nvec, lomask, bvec, avec);
    _mm512_store_epi64(sum->data + i * VECLEN, avec);
    _mm512_store_epi64(diff->data + i * VECLEN, bvec);

    if (NWORDS == 4) goto addsubdone;

    i = 4;
    nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
    //avec = _mm512_mask_sbb_epi52(s04, cmask, borrow, nvec, &borrow);
    //bvec = _mm512_mask_adc_epi52(d04, bmask, carry, nvec, &carry);
    _mm512_mask_adcsbb_epi52(d04, s04, bmask, cmask, carry, borrow, nvec, lomask, bvec, avec);
    _mm512_store_epi64(sum->data + i * VECLEN, avec);
    _mm512_store_epi64(diff->data + i * VECLEN, bvec);

    i = 5;
    nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
    //avec = _mm512_mask_sbb_epi52(s05, cmask, borrow, nvec, &borrow);
    //bvec = _mm512_mask_adc_epi52(d05, bmask, carry, nvec, &carry);
    _mm512_mask_adcsbb_epi52(d05, s05, bmask, cmask, carry, borrow, nvec, lomask, bvec, avec);
    _mm512_store_epi64(sum->data + i * VECLEN, avec);
    _mm512_store_epi64(diff->data + i * VECLEN, bvec);

    i = 6;
    nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
    //avec = _mm512_mask_sbb_epi52(s06, cmask, borrow, nvec, &borrow);
    //bvec = _mm512_mask_adc_epi52(d06, bmask, carry, nvec, &carry);
    _mm512_mask_adcsbb_epi52(d06, s06, bmask, cmask, carry, borrow, nvec, lomask, bvec, avec);
    _mm512_store_epi64(sum->data + i * VECLEN, avec);
    _mm512_store_epi64(diff->data + i * VECLEN, bvec);

    i = 7;
    nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
    //avec = _mm512_mask_sbb_epi52(s07, cmask, borrow, nvec, &borrow);
    //bvec = _mm512_mask_adc_epi52(d07, bmask, carry, nvec, &carry);
    _mm512_mask_adcsbb_epi52(d07, s07, bmask, cmask, carry, borrow, nvec, lomask, bvec, avec);
    _mm512_store_epi64(sum->data + i * VECLEN, avec);
    _mm512_store_epi64(diff->data + i * VECLEN, bvec);

    if (NWORDS == 8) goto addsubdone;

    i = 8;
    nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
    //avec = _mm512_mask_sbb_epi52(s08, cmask, borrow, nvec, &borrow);
    //bvec = _mm512_mask_adc_epi52(d08, bmask, carry, nvec, &carry);
    _mm512_mask_adcsbb_epi52(d08, s08, bmask, cmask, carry, borrow, nvec, lomask, bvec, avec);
    _mm512_store_epi64(sum->data + i * VECLEN, avec);
    _mm512_store_epi64(diff->data + i * VECLEN, bvec);

    i = 9;
    nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
    //avec = _mm512_mask_sbb_epi52(s09, cmask, borrow, nvec, &borrow);
    //bvec = _mm512_mask_adc_epi52(d09, bmask, carry, nvec, &carry);
    _mm512_mask_adcsbb_epi52(d09, s09, bmask, cmask, carry, borrow, nvec, lomask, bvec, avec);
    _mm512_store_epi64(sum->data + i * VECLEN, avec);
    _mm512_store_epi64(diff->data + i * VECLEN, bvec);

    i = 10;
    nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
    //avec = _mm512_mask_sbb_epi52(s10, cmask, borrow, nvec, &borrow);
    //bvec = _mm512_mask_adc_epi52(d10, bmask, carry, nvec, &carry);
    _mm512_mask_adcsbb_epi52(d10, s10, bmask, cmask, carry, borrow, nvec, lomask, bvec, avec);
    _mm512_store_epi64(sum->data + i * VECLEN, avec);
    _mm512_store_epi64(diff->data + i * VECLEN, bvec);

    i = 11;
    nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
    //avec = _mm512_mask_sbb_epi52(s11, cmask, borrow, nvec, &borrow);
    //bvec = _mm512_mask_adc_epi52(d11, bmask, carry, nvec, &carry);
    _mm512_mask_adcsbb_epi52(d11, s11, bmask, cmask, carry, borrow, nvec, lomask, bvec, avec);
    _mm512_store_epi64(sum->data + i * VECLEN, avec);
    _mm512_store_epi64(diff->data + i * VECLEN, bvec);

    if (NWORDS == 12) goto addsubdone;

    i = 12;
    nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
    //avec = _mm512_mask_sbb_epi52(s12, cmask, borrow, nvec, &borrow);
    //bvec = _mm512_mask_adc_epi52(d12, bmask, carry, nvec, &carry);
    _mm512_mask_adcsbb_epi52(d12, s12, bmask, cmask, carry, borrow, nvec, lomask, bvec, avec);
    _mm512_store_epi64(sum->data + i * VECLEN, avec);
    _mm512_store_epi64(diff->data + i * VECLEN, bvec);

    i = 13;
    nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
    //avec = _mm512_mask_sbb_epi52(s13, cmask, borrow, nvec, &borrow);
    //bvec = _mm512_mask_adc_epi52(d13, bmask, carry, nvec, &carry);
    _mm512_mask_adcsbb_epi52(d13, s13, bmask, cmask, carry, borrow, nvec, lomask, bvec, avec);
    _mm512_store_epi64(sum->data + i * VECLEN, avec);
    _mm512_store_epi64(diff->data + i * VECLEN, bvec);

    i = 14;
    nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
    //avec = _mm512_mask_sbb_epi52(s14, cmask, borrow, nvec, &borrow);
    //bvec = _mm512_mask_adc_epi52(d14, bmask, carry, nvec, &carry);
    _mm512_mask_adcsbb_epi52(d14, s14, bmask, cmask, carry, borrow, nvec, lomask, bvec, avec);
    _mm512_store_epi64(sum->data + i * VECLEN, avec);
    _mm512_store_epi64(diff->data + i * VECLEN, bvec);

    i = 15;
    nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
    //avec = _mm512_mask_sbb_epi52(s15, cmask, borrow, nvec, &borrow);
    //bvec = _mm512_mask_adc_epi52(d15, bmask, carry, nvec, &carry);
    _mm512_mask_adcsbb_epi52(d15, s15, bmask, cmask, carry, borrow, nvec, lomask, bvec, avec);
    _mm512_store_epi64(sum->data + i * VECLEN, avec);
    _mm512_store_epi64(diff->data + i * VECLEN, bvec);

    if (NWORDS == 16) goto addsubdone;

    i = 16;
    nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
    //avec = _mm512_mask_sbb_epi52(s16, cmask, borrow, nvec, &borrow);
    //bvec = _mm512_mask_adc_epi52(d16, bmask, carry, nvec, &carry);
    _mm512_mask_adcsbb_epi52(d16, s16, bmask, cmask, carry, borrow, nvec, lomask, bvec, avec);
    _mm512_store_epi64(sum->data + i * VECLEN, avec);
    _mm512_store_epi64(diff->data + i * VECLEN, bvec);

    i = 17;
    nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
    //avec = _mm512_mask_sbb_epi52(s17, cmask, borrow, nvec, &borrow);
    //bvec = _mm512_mask_adc_epi52(d17, bmask, carry, nvec, &carry);
    _mm512_mask_adcsbb_epi52(d17, s17, bmask, cmask, carry, borrow, nvec, lomask, bvec, avec);
    _mm512_store_epi64(sum->data + i * VECLEN, avec);
    _mm512_store_epi64(diff->data + i * VECLEN, bvec);

    i = 18;
    nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
    //avec = _mm512_mask_sbb_epi52(s18, cmask, borrow, nvec, &borrow);
    //bvec = _mm512_mask_adc_epi52(d18, bmask, carry, nvec, &carry);
    _mm512_mask_adcsbb_epi52(d18, s18, bmask, cmask, carry, borrow, nvec, lomask, bvec, avec);
    _mm512_store_epi64(sum->data + i * VECLEN, avec);
    _mm512_store_epi64(diff->data + i * VECLEN, bvec);

    i = 19;
    nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
    //avec = _mm512_mask_sbb_epi52(s19, cmask, borrow, nvec, &borrow);
    //bvec = _mm512_mask_adc_epi52(d19, bmask, carry, nvec, &carry);
    _mm512_mask_adcsbb_epi52(d19, s19, bmask, cmask, carry, borrow, nvec, lomask, bvec, avec);
    _mm512_store_epi64(sum->data + i * VECLEN, avec);
    _mm512_store_epi64(diff->data + i * VECLEN, bvec);

addsubdone:

    return;
}

void vec_simul_addsub52_avxecm(vec_bignum_t* a, vec_bignum_t* b, 
    vec_bignum_t* sum, vec_bignum_t* diff, vec_monty_t* mdata)
{
    // assumptions:
    // a, b, c are of length VECLEN * NWORDS
    // a, b, c, and n are aligned
    // a and b are both positive
    // n is the montgomery base
    // produce sum = a + b and diff = a - b at the same time which
    // saves 3N loads (only have to load a,b, and n once)
    int i;
    int NWORDS = 16;

    __mmask8 carry = 0;
    __mmask8 borrow = 0;
    __mmask8 cmask = 0;
    __mmask8 cmask2 = 0;
    __mmask8 bmask = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    __m512i nvec;


    for (i = 0; i < NWORDS; i++)
    {
        // add
        avec = _mm512_load_epi64(a->data + i * VECLEN);
        bvec = _mm512_load_epi64(b->data + i * VECLEN);
        cvec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(sum->data + i * VECLEN, cvec);

        // sub
        cvec = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
        _mm512_store_epi64(diff->data + i * VECLEN, cvec);

        //t1 = _mm512_add_epi64(avec, bvec);
        //t3 = _mm512_sub_epi64(avec, bvec);
        //t4 = _mm512_sub_epi64(t3, _mm512_maskz_set1_epi64(borrow, 1));
        //t1 = _mm512_add_epi64(t1, _mm512_maskz_set1_epi64(carry, 1));
        //borrow = _mm512_cmpgt_epu64_mask(bvec, avec);
        //carry = _mm512_cmpgt_epu64_mask(t1, vmask);
        //borrow = _mm512_kor(borrow, _mm512_cmpgt_epu64_mask(t4, t3));
        //t1 = _mm512_and_epi64(t1, vmask);
        //t4 = _mm512_and_epi64(t4, vmask);
        //
        //_mm512_store_epi64(sum->data + i * VECLEN, t1);
        //_mm512_store_epi64(diff->data + i * VECLEN, t4);
    }

    cmask = carry;	    // result too big, need to subtract n
    bmask = borrow;     // result too small, need to add n
    cmask2 = cmask;	    // keep looking mask for add

    // compare.  the initial mask is equal to the addition carry
    // because if the most significant word has a carry, then the
    // result is bigger than n.
    for (i = NWORDS - 1; i >= 0; i--)
    {
        cvec = _mm512_load_epi64(sum->data + i * VECLEN);
        nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
        // compare those that have not already been decided using the mask
        cmask |= _mm512_mask_cmp_epu64_mask(~cmask2, cvec, nvec, _MM_CMPINT_GT);
        cmask2 |= _mm512_mask_cmp_epu64_mask(~cmask2, cvec, nvec, _MM_CMPINT_LT);

        // decided all of them, stop comparing.
        if (cmask2 == 0xff)
        {
            break;
        }
    }

    // check for equal as well by flipping mask bits that have still
    // not been decided (i.e., are equal)
    cmask |= (~cmask2);

    carry = 0;
    borrow = 0;
    for (i = 0; i < NWORDS; i++)
    {
        // conditional sub
        cvec = _mm512_load_epi64(sum->data + i * VECLEN);
        nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
        bvec = _mm512_mask_sbb_epi52(cvec, cmask, borrow, nvec, &borrow);
        _mm512_store_epi64(sum->data + i * VECLEN, bvec);

        // conditional add
        cvec = _mm512_load_epi64(diff->data + i * VECLEN);
        bvec = _mm512_mask_adc_epi52(cvec, bmask, carry, nvec, &carry);
        _mm512_store_epi64(diff->data + i * VECLEN, bvec);
    }

    return;
}

void vec_simul_addsub52(vec_bignum_t *a, vec_bignum_t *b, 
    vec_bignum_t *sum, vec_bignum_t *diff, 
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
    uint32_t NWORDS = mdata->NWORDS; // a->WORDS_ALLOC;
    __mmask8 carry = 0;
    __mmask8 borrow = 0;
    __mmask8 cmask = 0;
    __mmask8 cmask2 = 0;
    __mmask8 bmask = 0;
    __mmask8 bmask2 = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    __m512i nvec;
    __m512i lomask = _mm512_set1_epi64(0xfffffffffffff);

    for (i = 0; i < NWORDS; i++)
    {
        // add/sub
        avec = _mm512_load_epi64(a->data + i * VECLEN);
        bvec = _mm512_load_epi64(b->data + i * VECLEN);

        _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, avec, bvec);

        _mm512_store_epi64(sum->data + i * VECLEN, avec);
        _mm512_store_epi64(diff->data + i * VECLEN, bvec);
    }

    cmask = carry;	    // result too big, need to subtract n
    bmask = borrow;     // result too small, need to add n
    cmask2 = cmask;	    // keep looking mask for add

#ifndef USE_AMM
    // compare.  the initial mask is equal to the addition carry
    // because if the most significant word has a carry, then the
    // result is bigger than n.
    for (i = NWORDS - 1; i >= 0; i--)
    {
        cvec = _mm512_load_epi64(sum->data + i * VECLEN);
        nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);
        // compare those that have not already been decided using the mask
        cmask |= _mm512_mask_cmp_epu64_mask(~cmask2, cvec, nvec, _MM_CMPINT_GT);
        cmask2 |= _mm512_mask_cmp_epu64_mask(~cmask2, cvec, nvec, _MM_CMPINT_LT);

        // decided all of them, stop comparing.
        if (cmask2 == 0xff)
        {
            break;
        }
    }

    // check for equal as well by flipping mask bits that have still
    // not been decided (i.e., are equal)
    cmask |= (~cmask2);
#endif

    carry = 0;
    borrow = 0;

    for (i = 0; i < NWORDS; i++)
    {
        // conditional sub/add
        avec = _mm512_load_epi64(sum->data + i * VECLEN);
        bvec = _mm512_load_epi64(diff->data + i * VECLEN);
        nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN);

        _mm512_mask_adcsbb_epi52(bvec, avec, bmask, cmask, carry, borrow, nvec, lomask, bvec, avec);

        _mm512_store_epi64(sum->data + i * VECLEN, avec);
        _mm512_store_epi64(diff->data + i * VECLEN, bvec);
    }

    return;
}



#else

void vecsqrmod52(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    return;
}

void vecmulmod52(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    return;
}

void vecaddmod52(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata)
{
    return;
}

void vecsubmod52(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata)
{
    return;
}

void vec_simul_addsub52(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* sum, vec_bignum_t* diff, vec_monty_t* mdata)
{
    return;
}

void vecsqrmod52_mersenne(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    return;
}

void vecmulmod52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    return;
}

void vecaddmod52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata)
{
    return;
}

void vecsubmod52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata)
{
    return;
}

void vec_simul_addsub52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* sum, vec_bignum_t* diff, vec_monty_t* mdata)
{
    return;
}

void vecmulmod52_fixed1040_bfips(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    return;
}

void vecmulmod52_fixed832_bfips(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    return;
}

void vecmulmod52_fixed624_bfips(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    return;
}

void vecmulmod52_fixed416_bfips(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    return;
}

void vecsqrmod52_fixed1040_bfips(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    return;
}

void vecsqrmod52_fixed832_bfips(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    return;
}

void vecsqrmod52_fixed624_bfips(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    return;
}

void vecsqrmod52_fixed416_bfips(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    return;
}

void vec_simul_addsub52_fixed1040(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* sum, vec_bignum_t* diff,
    vec_monty_t* mdata)
{
    return;
}

#endif
