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


__m512i __inline _mm512_mask_sbb_src_epi52(__m512i src, __m512i a, __mmask8 m, __mmask8 c, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_mask_sub_epi64(src, m, a, b);
    *cout = _mm512_mask_cmpgt_epu64_mask(m, b, a);
    __m512i t2 = _mm512_mask_sub_epi64(src, m, t, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_kor(*cout, _mm512_mask_cmpgt_epu64_mask(m, t2, t));
    t2 = _mm512_and_epi64(t2, _mm512_set1_epi64(0xfffffffffffffULL));
    return t2;
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
    
#ifdef USE_AMM

    // Need to check for an overflow (result > R = 2^MAXBITS)
    a0 = acc_e0;
    scarry2 = _mm512_cmp_epu64_mask(a0, zero, _MM_CMPINT_GT);

    if (scarry2)
    {
        printf("mul > R, reducing\n");
        // subtract n from result
        scarry = 0;
        for (i = 0; i < NWORDS; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            b0 = _mm512_load_epi64(n->data + i * VECLEN);
            a0 = _mm512_mask_sbb_epi52(a1, scarry2, scarry, b0, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }
    }
    //_mm512_store_epi64(c->data + NWORDS * VECLEN, zero);

#else
    // Need to check for an overflow (result > R = 2^MAXBITS) and
    // need to check if result > modulus
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

    //_mm512_store_epi64(c->data + NWORDS * VECLEN, zero);


#ifdef USE_AMM

    // Need to check for an overflow (result > R = 2^MAXBITS)
    a0 = acc_e0;
    scarry2 = _mm512_cmp_epu64_mask(a0, zero, _MM_CMPINT_GT);

    // subtract n from result
    scarry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(n->data + i * VECLEN);
        a0 = _mm512_mask_sbb_epi52(a1, scarry2, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
    }

#else
    // Need to check for an overflow (result > R = 2^MAXBITS) and
    // need to check if result > modulus
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

    //_mm512_store_epi64(c->data + NWORDS * VECLEN, zero);


#ifdef USE_AMM

    // Need to check for an overflow (result > R = 2^MAXBITS)
    a0 = acc_e0;
    scarry2 = _mm512_cmp_epu64_mask(a0, zero, _MM_CMPINT_GT);

    // subtract n from result
    scarry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(n->data + i * VECLEN);
        a0 = _mm512_mask_sbb_epi52(a1, scarry2, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
    }

#else
    // Need to check for an overflow (result > R = 2^MAXBITS) and
    // need to check if result > modulus
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

    // Need to check for an overflow (result > R = 2^MAXBITS)
    a0 = acc_e0;
    scarry2 = _mm512_cmp_epu64_mask(a0, zero, _MM_CMPINT_GT);

    // subtract n from result
    scarry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(n->data + i * VECLEN);
        a0 = _mm512_mask_sbb_epi52(a1, scarry2, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
    }

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

//#define DO_INTERMEDIATE_CARRYPROP
void vecmulmod52(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    int i, j, k;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t NBLOCKS = mdata->NBLOCKS;

    // needed in loops
    __m512i i0, i1;
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7, te8;             // 19

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

    if ((a == c) || (b == c))
        outdata = s->data;

#else
    uint64_t* outdata = s->data;
#endif

    //mpz_t ga, gb;
    //mpz_init(ga);
    //mpz_init(gb);
    //extract_bignum_from_vec_to_mpz(ga, a, 4, NWORDS);
    //extract_bignum_from_vec_to_mpz(gb, b, 4, NWORDS);

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
        int biascount = 0;
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = te8 = zero;

        for (j = i, k = 0; j > 0; j--, k++)
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

            biascount += 4;

#ifdef DO_INTERMEDIATE_CARRYPROP
            // after 12 accumulations we are in danger of overflowing the 
            // redundant bits.  Although in reality not every addition will 
            // produce a carry bit.  So every 3-4 blocks we should carry propagate.
            if (k > 2)
            {
                SUB_BIAS_HI(biascount, biascount, biascount, biascount);
                SUB_BIAS_LO(biascount, biascount, biascount, biascount);

                te2 = _mm512_add_epi64(te2, _mm512_srli_epi64(te0, 52));
                te3 = _mm512_add_epi64(te3, _mm512_srli_epi64(te1, 52));
                te4 = _mm512_add_epi64(te4, _mm512_srli_epi64(te2, 52));
                te5 = _mm512_add_epi64(te5, _mm512_srli_epi64(te3, 52));
                te6 = _mm512_add_epi64(te6, _mm512_srli_epi64(te4, 52));
                te7 = _mm512_add_epi64(te7, _mm512_srli_epi64(te5, 52));
                te7 = _mm512_add_epi64(te7, _mm512_srli_epi64(te6, 52));
                te8 = _mm512_add_epi64(te8, _mm512_srli_epi64(te7, 52));

                te0 = _mm512_and_epi64(te0, vlmask);
                te1 = _mm512_and_epi64(te1, vlmask);
                te2 = _mm512_and_epi64(te2, vlmask);
                te3 = _mm512_and_epi64(te3, vlmask);
                te4 = _mm512_and_epi64(te4, vlmask);
                te5 = _mm512_and_epi64(te5, vlmask);
                te6 = _mm512_and_epi64(te6, vlmask);
                te7 = _mm512_and_epi64(te7, vlmask);

                k = 0;
                biascount = 0;
            }
#endif
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

#ifdef DO_INTERMEDIATE_CARRYPROP

#ifndef IFMA

        // subtract out all of the bias.
        SUB_BIAS_HI(
            biascount + 1,
            biascount + 2,
            biascount + 3,
            biascount + 4);
        SUB_BIAS_LO(
            biascount + 1,
            biascount + 2,
            biascount + 3,
            biascount + 4);

        biascount = 0;
#endif

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

            biascount += 4;
#ifdef DO_INTERMEDIATE_CARRYPROP   

            // after 12 accumulations we are in danger of overflowing the 
            // redundant bits.  Although in reality not every addition will 
            // produce a carry bit.  So every 3-4 blocks we should carry propagate.
            if (k > 2)
            {
                SUB_BIAS_HI(biascount, biascount, biascount, biascount);
                SUB_BIAS_LO(biascount, biascount, biascount, biascount);

                te2 = _mm512_add_epi64(te2, _mm512_srli_epi64(te0, 52));
                te3 = _mm512_add_epi64(te3, _mm512_srli_epi64(te1, 52));
                te4 = _mm512_add_epi64(te4, _mm512_srli_epi64(te2, 52));
                te5 = _mm512_add_epi64(te5, _mm512_srli_epi64(te3, 52));
                te6 = _mm512_add_epi64(te6, _mm512_srli_epi64(te4, 52));
                te7 = _mm512_add_epi64(te7, _mm512_srli_epi64(te5, 52));
                te7 = _mm512_add_epi64(te7, _mm512_srli_epi64(te6, 52));
                te8 = _mm512_add_epi64(te8, _mm512_srli_epi64(te7, 52));

                te0 = _mm512_and_epi64(te0, vlmask);
                te1 = _mm512_and_epi64(te1, vlmask);
                te2 = _mm512_and_epi64(te2, vlmask);
                te3 = _mm512_and_epi64(te3, vlmask);
                te4 = _mm512_and_epi64(te4, vlmask);
                te5 = _mm512_and_epi64(te5, vlmask);
                te6 = _mm512_and_epi64(te6, vlmask);
                te7 = _mm512_and_epi64(te7, vlmask);

                k = 0;
                biascount = 0;
            }
#endif
        }

#ifdef DO_INTERMEDIATE_CARRYPROP
#ifndef IFMA

        // subtract out all of the bias at once.
        SUB_BIAS_HI(biascount, biascount, biascount, biascount);
        SUB_BIAS_LO(biascount, biascount, biascount, biascount);

#endif
#else
#ifndef IFMA

        // subtract out all of the bias at once.
        SUB_BIAS_HI(biascount + 1, biascount + 2, biascount + 3, biascount + 4);
        SUB_BIAS_LO(biascount + 1, biascount + 2, biascount + 3, biascount + 4);

#endif
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
            acc_e2 = _mm512_add_epi64(te8, _mm512_srli_epi64(acc_e1, 52));
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
        int biascount = 0;
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = te8 = zero;

        for (j = i - NBLOCKS + 1, k = 0; j < NBLOCKS; j++, k += 2)
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

            biascount += 8;
#ifdef DO_INTERMEDIATE_CARRYPROP

            // after 12 accumulations we are in danger of overflowing the 
            // redundant bits.  Although in reality not every addition will 
            // produce a carry bit.  So every 3-4 blocks we should carry propagate.
            //if (k > 3)
            {
                SUB_BIAS_HI(biascount, biascount, biascount, biascount);
                SUB_BIAS_LO(biascount, biascount, biascount, biascount);

                te2 = _mm512_add_epi64(te2, _mm512_srli_epi64(te0, 52));
                te3 = _mm512_add_epi64(te3, _mm512_srli_epi64(te1, 52));
                te4 = _mm512_add_epi64(te4, _mm512_srli_epi64(te2, 52));
                te5 = _mm512_add_epi64(te5, _mm512_srli_epi64(te3, 52));
                te6 = _mm512_add_epi64(te6, _mm512_srli_epi64(te4, 52));
                te7 = _mm512_add_epi64(te7, _mm512_srli_epi64(te5, 52));
                te7 = _mm512_add_epi64(te7, _mm512_srli_epi64(te6, 52));
                te8 = _mm512_add_epi64(te8, _mm512_srli_epi64(te7, 52));

                te0 = _mm512_and_epi64(te0, vlmask);
                te1 = _mm512_and_epi64(te1, vlmask);
                te2 = _mm512_and_epi64(te2, vlmask);
                te3 = _mm512_and_epi64(te3, vlmask);
                te4 = _mm512_and_epi64(te4, vlmask);
                te5 = _mm512_and_epi64(te5, vlmask);
                te6 = _mm512_and_epi64(te6, vlmask);
                te7 = _mm512_and_epi64(te7, vlmask);

                k = 0;
                biascount = 0;
            }
#endif
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
            biascount + 6,
            biascount + 4,
            biascount + 2,
            biascount + 0);
        SUB_BIAS_LO(
            biascount + 6,
            biascount + 4,
            biascount + 2,
            biascount + 0);

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
            acc_e2 = _mm512_add_epi64(te8, _mm512_srli_epi64(acc_e1, 52));
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

    if ((a == c) || (b == c))
        vecCopy(s, c);

    // Need to check for an overflow (result > R = 2^MAXBITS)
    a0 = acc_e0;
    scarry2 = _mm512_cmp_epu64_mask(a0, zero, _MM_CMPINT_GT);

    // subtract n from result
    scarry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(n->data + i * VECLEN);
        a0 = _mm512_mask_sbb_epi52(a1, scarry2, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
    }

#else

    if ((a == c) || (b == c))
        vecCopy(s, c);

    // Need to check for an overflow (result > R = 2^MAXBITS) and
    // need to check if result > modulus
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

    //mpz_t gt, gm, gn, gtmp1;
    //mpz_init(gt);
    //mpz_init(gm);
    //mpz_init(gn);
    //mpz_init(gtmp1);
    //
    //extract_bignum_from_vec_to_mpz(gn, n, 4, NWORDS);
    //mpz_mul(gt, ga, gb);
    //mpz_mod_2exp(gtmp1, gt, NWORDS * 52);
    //mpz_mul(gm, gtmp1, mdata->nhat);
    //mpz_mod_2exp(gm, gm, NWORDS * 52);
    //mpz_mul(gtmp1, gm, gn);
    //mpz_add(gt, gt, gtmp1);
    //mpz_tdiv_q_2exp(gt, gt, NWORDS * 52);
    //
    //extract_bignum_from_vec_to_mpz(gtmp1, c, 4, NWORDS);
    //if (mpz_cmp(gtmp1, gt) != 0)
    //{
    //    printf("modmul check failed, rho = %016llx\n", mdata->vrho[4]);
    //    gmp_printf("a = %Zd\n", ga);
    //    gmp_printf("b = %Zd\n", gb);
    //    gmp_printf("n = %Zd\n", gn);
    //    gmp_printf("c = %Zd\n", gtmp1);
    //    gmp_printf("g = %Zd\n", gt);
    //    exit(1);
    //}
    //else
    //{
    //    printf("modmul check success\n");
    //}
    //mpz_clear(ga);
    //mpz_clear(gb);
    //mpz_clear(gt);
    //mpz_clear(gn);
    //mpz_clear(gm);
    //mpz_clear(gtmp1);

    c->size = NWORDS;
    return;
}

void vecsqrmod52_fixed1040_bfips(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    vecmulmod52_fixed1040_bfips(a, a, c, n, s, mdata);
    return;

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

#ifdef USE_AMM

    // Need to check for an overflow (result > R = 2^MAXBITS)
    a0 = acc_e0;
    scarry2 = _mm512_cmp_epu64_mask(a0, zero, _MM_CMPINT_GT);

    // subtract n from result
    scarry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(n->data + i * VECLEN);
        a0 = _mm512_mask_sbb_epi52(a1, scarry2, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
    }

    //_mm512_store_epi64(c->data + NWORDS * VECLEN, zero);

#else
    // Need to check for an overflow (result > R = 2^MAXBITS) and
    // need to check if result > modulus
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

#ifdef USE_AMM

// Need to check for an overflow (result > R = 2^MAXBITS)
    a0 = acc_e0;
    scarry2 = _mm512_cmp_epu64_mask(a0, zero, _MM_CMPINT_GT);

    // subtract n from result
    scarry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(n->data + i * VECLEN);
        a0 = _mm512_mask_sbb_epi52(a1, scarry2, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
    }

#else
    // Need to check for an overflow (result > R = 2^MAXBITS) and
    // need to check if result > modulus
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

#ifdef USE_AMM

    // Need to check for an overflow (result > R = 2^MAXBITS)
    a0 = acc_e0;
    scarry2 = _mm512_cmp_epu64_mask(a0, zero, _MM_CMPINT_GT);

    // subtract n from result
    scarry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(n->data + i * VECLEN);
        a0 = _mm512_mask_sbb_epi52(a1, scarry2, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
    }

#else
    // Need to check for an overflow (result > R = 2^MAXBITS) and
    // need to check if result > modulus
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

    // Need to check for an overflow (result > R = 2^MAXBITS)
    a0 = acc_e0;
    scarry2 = _mm512_cmp_epu64_mask(a0, zero, _MM_CMPINT_GT);

    // subtract n from result
    scarry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(n->data + i * VECLEN);
        a0 = _mm512_mask_sbb_epi52(a1, scarry2, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
    }

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

    if (a == c)
        outdata = s->data;

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

    if (a == c)
        vecCopy(s, c);

    // Need to check for an overflow (result > R = 2^MAXBITS)
    a0 = acc_e0;
    scarry2 = _mm512_cmp_epu64_mask(a0, zero, _MM_CMPINT_GT);

    // subtract n from result
    scarry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(n->data + i * VECLEN);
        a0 = _mm512_mask_sbb_epi52(a1, scarry2, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
    }


#else
    // Need to check for an overflow (result > R = 2^MAXBITS) and
    // need to check if result > modulus
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

    //printf("vecsub:\n");
    //for (i = NWORDS - 1; i >=0; i--)
    //    printf("%lx", c->data[i * VECLEN]);
    //printf("\n");

    // if we did not have a final carry, then we still are not in
    // positive territory.  add n again.
    mask = (~carry) & mask;
    //if (mask & 1) printf("vecsub still negative, adding n again\n");
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
    vecaddmod52(a, b, sum, mdata);
    vecsubmod52(a, b, diff, mdata);
    return;

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

    //printf("vecaddsub diff:\n");
    //for (i = 7; i >= 0; i--)
    //    printf("%lx", diff->data[i * VECLEN]);
    //printf("\n");

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

void vecmodexp52(vec_bignum_t* d, vec_bignum_t* b, vec_bignum_t* e, vec_bignum_t* m,
    vec_bignum_t* s, vec_bignum_t* one, vec_monty_t* mdata)
{

    return;
}


#endif
