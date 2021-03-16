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
#define USE_AVX512F
#endif

#ifdef USE_AVX512F
#include <immintrin.h>

//#define PRINT_DEBUG 2
//#define NPRINT_DEBUG 0

#ifndef SKYLAKEX
#define _mm512_mullo_epi64 _mm512_mullox_epi64
// also need an answer for _mm512_cvtepu64_pd, probably among other things,
// because KNL does not support AVX512DQ
#define _mm512_cvtepu64_pd(x) _mm512_sub_pd(_mm512_castsi512_pd(_mm512_add_epi64((x), _mm512_set1_epi64(0x4330000000000000ULL))), _mm512_set1_pd(4503599627370496.))
#define _mm512_cvtpd_epu64(x) _mm512_sub_epi64(_mm512_castpd_si512(_mm512_add_pd((x), _mm512_set1_pd(4503599627370496.))), _mm512_set1_epi64(0x4330000000000000ULL))
#endif

void print_vechex(base_t *a, int v, int n, const char *pre);
void print_regvechex(__m512i a, int v, const char *pre);
void print_regvechex64(__m512i a, int v, const char *pre);
void print_regvechexrange(__m512i a, int v1, int v2, const char *pre);
void vec_bignum52_mask_rshift_n(vec_bignum_t* u, vec_bignum_t* v, int n, uint32_t wmask);
void vecmul52_1(vec_bignum_t* a, __m512i b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata);

// ---------------------------------------------------------------------
// emulated instructions
// ---------------------------------------------------------------------
__m512i __inline _mm512_addsetc_epi52(__m512i a, __m512i b, __mmask8 *cout)
{
    __m512i t = _mm512_add_epi64(a, b);
    *cout = _mm512_cmpgt_epu64_mask(t, _mm512_set1_epi64(0xfffffffffffffULL));
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}

__m512i __inline _mm512_adc_epi52(__m512i a, __mmask8 c, __m512i b, __mmask8 *cout)
{
    __m512i t = _mm512_add_epi64(a, b);
    t = _mm512_add_epi64(t, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_cmpgt_epu64_mask(t, _mm512_set1_epi64(0xfffffffffffffULL));
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}

__m512i __inline _mm512_mask_adc_epi52(__m512i a, __mmask8 m, __mmask8 c, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_add_epi64(a, b);
    t = _mm512_mask_add_epi64(a, m, t, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_cmpgt_epu64_mask(t, _mm512_set1_epi64(0xfffffffffffffULL));
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}

__m512i __inline _mm512_addcarry_epi52(__m512i a, __mmask8 c, __mmask8* cout)
{
    __m512i t = _mm512_add_epi64(a, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_cmpeq_epu64_mask(a, _mm512_set1_epi64(0xfffffffffffffULL));
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}

__m512i __inline _mm512_subborrow_epi52(__m512i a, __mmask8 c, __mmask8* cout)
{
    __m512i t = _mm512_sub_epi64(a, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_cmpeq_epu64_mask(a, _mm512_set1_epi64(0));
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}

__m512i __inline _mm512_sbb_epi52(__m512i a, __mmask8 c, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_sub_epi64(a, b);
    *cout = _mm512_cmpgt_epu64_mask(b, a);
    __m512i t2 = _mm512_sub_epi64(t, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_kor(*cout, _mm512_cmpgt_epu64_mask(t2, t));
    t2 = _mm512_and_epi64(t2, _mm512_set1_epi64(0xfffffffffffffULL));
    return t2;
}

__m512i __inline _mm512_mask_sbb_epi52(__m512i a, __mmask8 m, __mmask8 c, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_mask_sub_epi64(a, m, a, b);
    *cout = _mm512_mask_cmpgt_epu64_mask(m, b, a);
    __m512i t2 = _mm512_mask_sub_epi64(a, m, t, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_kor(*cout, _mm512_mask_cmpgt_epu64_mask(m, t2, t));
    t2 = _mm512_and_epi64(t2, _mm512_set1_epi64(0xfffffffffffffULL));
    return t2;
}


//#define PRINT_DEBUG 0
//#define SPRINT_DEBUG 0


#define castpd _mm512_castsi512_pd
#define castepu _mm512_castpd_si512

#ifdef IFMA

#define _mm512_mullo_epi52(c, a, b) \
    c = _mm512_madd52lo_epu64(_mm512_set1_epi64(0), a, b);

#define VEC_MUL_ACCUM_LOHI_PD(a, b, lo, hi) \
    lo = _mm512_madd52lo_epu64(lo, a, b); \
    hi = _mm512_madd52hi_epu64(hi, a, b);

#define VEC_MUL_LOHI_PD(a, b, lo, hi) \
    lo = _mm512_madd52lo_epu64(_mm512_set1_epi64(0), a, b); \
    hi = _mm512_madd52hi_epu64(_mm512_set1_epi64(0), a, b);

#define VEC_CARRYPROP_LOHI(lo, hi) \
	a0 = _mm512_srli_epi64(lo, 52);	\
	hi = _mm512_add_epi64(hi, a0);		\
	lo = _mm512_and_epi64(vlmask, lo);

#define VEC_MUL4_ACCUM(x, b0, b1, b2, b3) \
    te0 = _mm512_madd52lo_epu64(te0, x, b0); \
    te2 = _mm512_madd52lo_epu64(te2, x, b1); \
    te4 = _mm512_madd52lo_epu64(te4, x, b2); \
    te6 = _mm512_madd52lo_epu64(te6, x, b3); \
    te1 = _mm512_madd52hi_epu64(te1, x, b0); \
    te3 = _mm512_madd52hi_epu64(te3, x, b1); \
    te5 = _mm512_madd52hi_epu64(te5, x, b2); \
    te7 = _mm512_madd52hi_epu64(te7, x, b3);

#define SUB_BIAS_HI(bias1, bias2, bias3, bias4) {}
#define SUB_BIAS_LO(bias1, bias2, bias3, bias4) {}

#else

#define _mm512_mullo_epi52(c, a, b) \
    i0 = _mm512_srli_epi64(a, 32); \
    i1 = _mm512_srli_epi64(b, 32); \
    i0 = _mm512_mul_epu32(b, i0); \
    i1 = _mm512_mul_epu32(a, i1); \
    c = _mm512_mul_epu32(a, b); \
    i0 = _mm512_slli_epi64(i0, 32); \
    i1 = _mm512_slli_epi64(i1, 32); \
    c = _mm512_add_epi64(c, i0); \
    c = _mm512_add_epi64(c, i1); \
    c = _mm512_and_epi64(c, vlmask); 

//#define _mm512_mullo_epi52(c, a, b) \
//    prod1_hd = _mm512_cvtepu64_pd(a); \
//    prod2_hd = _mm512_cvtepu64_pd(b); \
//    prod3_hd = _mm512_fmadd_round_pd(prod1_hd, prod2_hd, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
//    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd); \
//    prod3_hd = _mm512_fmadd_round_pd(prod1_hd, prod2_hd, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
//    c = _mm512_sub_epi64(_mm512_castpd_si512(prod3_hd), _mm512_set1_epi64(0x4330000000000000ULL));

#define VEC_MUL_ACCUM_LOHI_PD(a, b, lo, hi) \
	prod1_ld = _mm512_cvtepu64_pd(a);		\
	prod2_ld = _mm512_cvtepu64_pd(b);		\
    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
    hi = _mm512_add_epi64(hi, _mm512_sub_epi64(castepu(prod1_hd), vbias1)); \
    prod1_hd = _mm512_sub_pd(castpd(vbias2), prod1_hd); \
	prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
	lo = _mm512_add_epi64(lo, _mm512_sub_epi64(castepu(prod1_ld), vbias3));

#define VEC_MUL_LOHI_PD(a, b, lo, hi) \
	prod1_ld = _mm512_cvtepu64_pd(a);		\
	prod2_ld = _mm512_cvtepu64_pd(b);		\
    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
    hi = _mm512_sub_epi64(castepu(prod1_hd), vbias1); \
    prod1_hd = _mm512_sub_pd(castpd(vbias2), prod1_hd); \
	prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
	lo = _mm512_sub_epi64(castepu(prod1_ld), vbias3);

#define VEC_CARRYPROP_LOHI(lo, hi) \
	a0 = _mm512_srli_epi64(lo, 52);	\
	hi = _mm512_add_epi64(hi, a0);		\
	lo = _mm512_and_epi64(vlmask, lo);

#define VEC_MUL4_ACCUM(x, b0, b1, b2, b3) \
    prod5_ld = _mm512_cvtepu64_pd(x);                                                                           \
    prod1_ld = _mm512_cvtepu64_pd(b0);                                                                          \
    prod2_ld = _mm512_cvtepu64_pd(b1);                                                                          \
    prod3_ld = _mm512_cvtepu64_pd(b2);                                                                          \
    prod4_ld = _mm512_cvtepu64_pd(b3);                                                                          \
    prod1_hd = _mm512_fmadd_round_pd(prod5_ld, prod1_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));      \
    prod2_hd = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));      \
    prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));      \
    prod4_hd = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));      \
    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));                                                 \
    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);                                            \
    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod2_hd));                                                 \
    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);                                            \
    te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod3_hd));                                                 \
    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);                                            \
    te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod4_hd));                                                 \
    prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);                                            \
    prod1_ld = _mm512_fmadd_round_pd(prod5_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
    prod2_ld = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
    prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
    prod4_ld = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));                                                 \
    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod2_ld));                                                 \
    te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod3_ld));                                                 \
    te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod4_ld));

// avoid the slow _mm512_mullo_epi64 by using _mm512_mul_epu32 and a 32-bit shift.
// works because the bias's low 32-bits are all zero.
#define SUB_BIAS_HI(bias1, bias2, bias3, bias4) \
    b0 = _mm512_set1_epi64((bias1) * 0x467);         \
    b1 = _mm512_set1_epi64((bias2) * 0x467);         \
    b2 = _mm512_set1_epi64((bias3) * 0x467);         \
    b3 = _mm512_set1_epi64((bias4) * 0x467);         \
    b0 = _mm512_slli_epi64(b0, 52);            \
    b1 = _mm512_slli_epi64(b1, 52);            \
    b2 = _mm512_slli_epi64(b2, 52);            \
    b3 = _mm512_slli_epi64(b3, 52);            \
    te1 = _mm512_sub_epi64(te1, b0); \
    te3 = _mm512_sub_epi64(te3, b1); \
    te5 = _mm512_sub_epi64(te5, b2); \
    te7 = _mm512_sub_epi64(te7, b3);

#define SUB_BIAS_LO(bias1, bias2, bias3, bias4) \
    b0 = _mm512_set1_epi64((bias1) * 0x433);         \
    b1 = _mm512_set1_epi64((bias2) * 0x433);         \
    b2 = _mm512_set1_epi64((bias3) * 0x433);         \
    b3 = _mm512_set1_epi64((bias4) * 0x433);         \
    b0 = _mm512_slli_epi64(b0, 52);            \
    b1 = _mm512_slli_epi64(b1, 52);            \
    b2 = _mm512_slli_epi64(b2, 52);            \
    b3 = _mm512_slli_epi64(b3, 52);            \
    te0 = _mm512_sub_epi64(te0, b0); \
    te2 = _mm512_sub_epi64(te2, b1); \
    te4 = _mm512_sub_epi64(te4, b2); \
    te6 = _mm512_sub_epi64(te6, b3);

#endif

//#define DEBUGLANE 0
//#define DEBUG_MERSENNE

void vecmulmod52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
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
    vec_bignum52_mask_rshift_n(s, c, mdata->nbits, 0xff);

    int bshift = mdata->nbits % 52;
    int wshift = mdata->nbits / 52;

#ifdef DEBUG_MERSENNE
    print_vechexbignum52(c, "hi part:");
    print_vechexbignum52(s, "lo part:");
#endif

    if (mdata->isMersenne > 1)
    {
        // this number has been identified as a "pseudo"-Mersenne: 2^n-c, with
        // 'c' small.  Fast reduction is still possible.

        // multiply c * hi
        b0 = _mm512_set1_epi64(mdata->isMersenne);
        vecmul52_1(c, b0, c, n, s, mdata);

#ifdef DEBUG_MERSENNE
        print_vechexbignum52(c, "after mul_1 with c:");
#endif

        // clear any hi-bits in the lo result
        a1 = _mm512_load_epi64(s->data + wshift * VECLEN);
        _mm512_store_epi64(s->data + wshift * VECLEN,
            _mm512_and_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL), a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum52(s, "cleared low part:");
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
        print_vechexbignum52(c, "lo + hi * c:");
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
        print_vechexbignum52(c, "lo + hi * c:");
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
            _mm512_and_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL), a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum52(c, "after carry add:");
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

        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(s->data + i * VECLEN);
        b0 = _mm512_and_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL), b0);
        a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum52(c, "after add:");
#endif

        // if there was a carry, add it back in.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        scarry = _mm512_test_epi64_mask(a1, _mm512_set1_epi64((1ULL << (uint64_t)bshift)));
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
            _mm512_and_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL), a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum52(c, "after carry add:");
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
        b0 = _mm512_and_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL), b0);
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
        scarry = _mm512_test_epi64_mask(a1, _mm512_set1_epi64((1ULL << (uint64_t)bshift)));
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
            _mm512_and_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL), a1));
    }

    c->size = NWORDS;
    return;
}

void vecsqrmod52_mersenne(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
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
    vec_bignum52_mask_rshift_n(s, c, mdata->nbits, 0xff);

    int bshift = mdata->nbits % 52;
    int wshift = mdata->nbits / 52;

#ifdef DEBUG_MERSENNE
    print_vechexbignum52(c, "hi part:");
    print_vechexbignum52(s, "lo part:");
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
        print_vechexbignum52(s, "cleared low part:");
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
        print_vechexbignum52(c, "lo + hi * c:");
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
        print_vechexbignum52(c, "lo + hi * c:");
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
        print_vechexbignum52(c, "after carry add:");
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
        b0 = _mm512_and_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL), b0);
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
        scarry = _mm512_test_epi64_mask(a1, _mm512_set1_epi64((1ULL << (uint64_t)bshift)));
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
            _mm512_and_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL), a1));

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
        b0 = _mm512_and_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL), b0);
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
        scarry = _mm512_test_epi64_mask(a1, _mm512_set1_epi64((1ULL << (uint64_t)bshift)));
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
            _mm512_and_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL), a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after bit clear:");
        exit(1);
#endif

    }

    c->size = NWORDS;
    return;
}

void vecmulmod52(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata)
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

        a1 = _mm512_load_epi64(s->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(s->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(s->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

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

void vecsqrmod52_test(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
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

void vecaddmod52(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, vec_monty_t* mdata)
{
    // assumptions:
    // a, b, c are of length VECLEN * NWORDS
    // a, b, c, and n are aligned
    // a and b are both positive
    // n is the montgomery base
    int i;
    uint32_t NWORDS = a->WORDS_ALLOC;
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
    int bshift = mdata->nbits % 52;
    int wshift = mdata->nbits / 52;

    // add
    for (i = 0; i < NWORDS; i++)
    {
        avec = _mm512_load_epi64(a->data + i * VECLEN);
        bvec = _mm512_load_epi64(b->data + i * VECLEN);
        cvec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c->data + i * VECLEN, cvec);
    }

    // check for a carry.
    avec = _mm512_load_epi64(c->data + wshift * VECLEN);
    carry = _mm512_test_epi64_mask(avec, _mm512_set1_epi64((1ULL << (uint64_t)bshift)));

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
            bvec = _mm512_addcarry_epi52(cvec, carry, &carry);
            _mm512_store_epi64(c->data + i * VECLEN, bvec);
        }
    }
    else
    {
        // the modulo is just the low part minus 1 (the carry, if present).
        cvec = _mm512_load_epi64(c->data + 0 * VECLEN);
        //bvec = _mm512_subborrow_epi52(cvec, carry, &carry);
        //mm512_mask_sbb_epi52(__m512i a, __mmask8 m, __mmask8 c, __m512i b, __mmask8* cout)
        bvec = _mm512_set1_epi64(-mdata->isMersenne);
        bvec = _mm512_mask_sbb_epi52(cvec, carry, 0, bvec, &carry);
        _mm512_store_epi64(c->data + 0 * VECLEN, bvec);

        for (i = 1; (i < NWORDS) && (carry > 0); i++)
        {
            cvec = _mm512_load_epi64(c->data + i * VECLEN);
            bvec = _mm512_subborrow_epi52(cvec, carry, &carry);
            _mm512_store_epi64(c->data + i * VECLEN, bvec);
        }
    }

    // clear the potential hi-bit
    avec = _mm512_load_epi64(c->data + wshift * VECLEN);
    _mm512_store_epi64(c->data + wshift * VECLEN,
        _mm512_and_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL), avec));

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
    uint32_t NWORDS = a->WORDS_ALLOC;
    __mmask8 carry = 0;
    __mmask8 mask = 0;
    __mmask8 mask2 = 0;
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
    __mmask8 mask = 0;
    __mmask8 mask2 = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    __m512i nvec;
    int bshift = mdata->nbits % 52;
    int wshift = mdata->nbits / 52;

    // subtract
    carry = 0;
    //for (i = 0; i < NWORDS; i++)
    for (i = 0; i <= wshift; i++)
    {
        avec = _mm512_load_epi64(a->data + i * VECLEN);
        bvec = _mm512_load_epi64(b->data + i * VECLEN);
        cvec = _mm512_sbb_epi52(avec, carry, bvec, &carry);
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
        nvec = _mm512_set1_epi64(0xfffffffffffffULL);
        for (i = 1; i <= wshift; i++)
        {
            cvec = _mm512_load_epi64(c->data + i * VECLEN);
            bvec = _mm512_mask_adc_epi52(cvec, mask, carry, nvec, &carry);
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
            bvec = _mm512_addcarry_epi52(cvec, carry, &carry);
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
        _mm512_and_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL), avec));

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

void vec_simul_addsub52(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *sum, vec_bignum_t *diff, 
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
    __mmask8 borrow = 0;
    __mmask8 cmask = 0;
    __mmask8 cmask2 = 0;
    __mmask8 bmask = 0;
    __mmask8 bmask2 = 0;
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
    __mmask8 borrow = 0;
    __mmask8 bmask = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    __m512i nvec;
    int bshift = mdata->nbits % 52;
    int wshift = mdata->nbits / 52;

    //for (i = 0; i < NWORDS; i++)
    for (i = 0; i <= wshift; i++)
    {
        // add
        avec = _mm512_load_epi64(a->data + i * VECLEN);
        bvec = _mm512_load_epi64(b->data + i * VECLEN);
        cvec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(sum->data + i * VECLEN, cvec);

        // sub
        cvec = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
        _mm512_store_epi64(diff->data + i * VECLEN, cvec);
    }

    bmask = borrow;     // result too small, need to add n

    // check for a carry.
    avec = _mm512_load_epi64(sum->data + wshift * VECLEN);
    carry = _mm512_test_epi64_mask(avec, _mm512_set1_epi64((1ULL << (uint64_t)bshift)));

    if (mdata->isMersenne > 0)
    {
        // the modulo is just the low part plus 1 (the carry, if present).
        cvec = _mm512_load_epi64(sum->data + 0 * VECLEN);
        //bvec = _mm512_addcarry_epi52(cvec, carry, &carry);
        bvec = _mm512_set1_epi64(mdata->isMersenne);
        bvec = _mm512_mask_adc_epi52(cvec, carry, 0, bvec, &carry);
        _mm512_store_epi64(sum->data + 0 * VECLEN, bvec);

        for (i = 1; (i < NWORDS) && (carry > 0); i++)
        {
            cvec = _mm512_load_epi64(sum->data + i * VECLEN);
            bvec = _mm512_addcarry_epi52(cvec, carry, &carry);
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
            bvec = _mm512_subborrow_epi52(cvec, carry, &carry);
            _mm512_store_epi64(sum->data + i * VECLEN, bvec);
        }
    }

    // clear the potential hi-bit
    avec = _mm512_load_epi64(sum->data + wshift * VECLEN);
    _mm512_store_epi64(sum->data + wshift * VECLEN,
        _mm512_and_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL), avec));

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
            bvec = _mm512_mask_adc_epi52(cvec, bmask, carry, nvec, &carry);
            _mm512_store_epi64(diff->data + i * VECLEN, bvec);
        }
    }
    else
    {
        // if we had a final carry, then b was bigger than a so we need to re-add n.
        carry = 0;
        nvec = _mm512_set1_epi64(1);
        cvec = _mm512_load_epi64(diff->data + 0 * VECLEN);
        bvec = _mm512_mask_adc_epi52(cvec, bmask, carry, nvec, &carry);
        _mm512_store_epi64(diff->data + 0 * VECLEN, bvec);
        for (i = 1; (i <= wshift) && (carry > 0); i++)
        {
            cvec = _mm512_load_epi64(diff->data + i * VECLEN);
            bvec = _mm512_addcarry_epi52(cvec, carry, &carry);
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
        _mm512_and_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL), avec));

    return;
}

void print_vechex(base_t *a, int v, int n, const char *pre)
{
    // print n hexdigits of the v'th position of vec bignum a
    int i;
    if (pre != NULL)
        printf("%s: ", pre);
    for (i = n - 1; i >= 0; i--)
    {
        printf("%013"PRIx64"", a[v + i * VECLEN]);
    }
    printf("\n");
    return;
}

void print_regvechex(__m512i a, int v, const char *pre)
{
    // print n hexdigits of the v'th position of vec bignum a
    base_t aa[VECLEN];

    _mm512_store_epi64(aa, a);

    if (pre != NULL)
        printf("%s: ", pre);

    printf("%08"PRIx64"", aa[v]);
    if (pre != NULL)
        printf("\n");
    return;
}

void print_regvechexrange(__m512i a, int v1, int v2, const char *pre)
{
    // print n hexdigits of the v'th position of vec bignum a
    base_t aa[VECLEN];
    int i;

    _mm512_store_epi64(aa, a);

    if (pre != NULL)
        printf("%s: ", pre);

    for (i = v2; i >= v1; i--)
        printf("%08"PRIx64"", aa[i]);

    if (pre != NULL)
        printf("\n");
    return;
}

void print_hexbignum(vec_bignum_t *a, const char *pre)
{
    // print n hexdigits of the bignum a
    int i;
    uint32_t NWORDS = a->WORDS_ALLOC;
    if (pre != NULL)
        printf("%s", pre);

    for (i = NWORDS - 1; i >= 0; i--)
    {
        printf("%013"PRIx64"", a->data[i]);
    }
    printf("\n");
    return;
}

void print_vechexbignum(vec_bignum_t *a, const char *pre)
{
    // print n hexdigits of the bignum a
    int i, j;
    uint32_t NWORDS = a->WORDS_ALLOC;
    if (pre != NULL)
        printf("%s", pre);

    for (j = 0; j < VECLEN; j++)
    {
#ifdef DEBUGLANE
        if (j != DEBUGLANE)
            continue;
#endif
        for (i = 2 * NWORDS - 1; i >= 0; i--)
        {
            printf("%016llx", a->data[i * VECLEN + j]);
        }
        printf("\n");
    }
    return;
}

void print_regvechex64(__m512i a, int v, const char *pre)
{
    // print n hexdigits of the v'th position of vec bignum a
    uint64_t aa[VECLEN];

    _mm512_store_epi64(aa, a);

    if (pre != NULL)
        printf("%s: ", pre);

    printf("%016"PRIx64"", aa[v]);
    if (pre != NULL)
        printf("\n");
    return;
}

uint32_t vec_gte52(vec_bignum_t * u, vec_bignum_t * v)
{
    // decide if each of the bignums in vec 'u' is >=
    // the corresponding bignum in vec 'v'.
    // return a mask of results.
    int i;
    uint32_t NWORDS = u->WORDS_ALLOC;
    __mmask8 mdecided = 0;
    __mmask8 mgte = 0;

    for (i = NWORDS - 1; i >= 0; --i)
    {
        __m512i a = _mm512_load_epi64(u->data + i * VECLEN);
        __m512i b = _mm512_load_epi64(v->data + i * VECLEN);

        mgte |= _mm512_mask_cmp_epu64_mask(~mdecided, a, b, _MM_CMPINT_GT);
        mdecided = mdecided | _mm512_mask_cmp_epu64_mask(~mdecided, a, b, _MM_CMPINT_LT);

        if (mdecided == 0xff)
            break;
    }

    //equal if still undecided
    mgte |= ~mdecided;

    return (uint32_t)mgte;
}

uint32_t vec_eq52(base_t * u, base_t * v, int sz)
{
    // decide if each of the bignums in vec 'u' is >=
    // the corresponding bignum in vec 'v'.
    // return a mask of results.
    int i;
    __mmask8 meq = 0xff;

    for (i = sz - 1; i >= 0; --i)
    {
        __m512i a = _mm512_load_epi64(u + i * VECLEN);
        __m512i b = _mm512_load_epi64(v + i * VECLEN);

        meq = _mm512_mask_cmp_epu64_mask(meq, a, b, _MM_CMPINT_EQ);

        if (meq == 0)
            break;
    }

    return (uint32_t)meq;
}

uint32_t vec_bignum52_mask_lshift_1(vec_bignum_t * u, uint32_t wmask)
{
    // return the left shift of bignum u by 1
    int i;
    uint32_t NWORDS = u->WORDS_ALLOC;
    __m512i nextcarry;
    __m512i carry = _mm512_set1_epi64(0);
    __m512i highmask = _mm512_set1_epi64(VEC_MAXDIGIT);
    __m512i word;

    for (i = 0; i < NWORDS; i++)
    {
        word = _mm512_load_epi64(u->data + i * VECLEN);
        nextcarry = _mm512_srli_epi64(word, (DIGITBITS - 1));
        _mm512_mask_store_epi64(u->data + i * VECLEN, (__mmask8)wmask,
            _mm512_and_epi64(highmask, _mm512_or_epi64(_mm512_slli_epi64(word, 1), carry)));
        carry = nextcarry;
    }

    _mm512_mask_store_epi64(u->data + i * VECLEN, (__mmask8)wmask,
        _mm512_and_epi64(highmask, _mm512_or_epi64(_mm512_slli_epi64(word, 1), carry)));

    // return an overflow mask
    return wmask & _mm512_cmp_epi64_mask(carry, _mm512_set1_epi64(0), _MM_CMPINT_GT);
}

uint32_t vec_bignum52_mask_lshift_n(vec_bignum_t * u, int n, uint32_t wmask)
{
    // return the left shift of bignum u by n bits
    // n is assumed less than DIGITBITS
    int i;
    uint32_t NWORDS = u->WORDS_ALLOC;
    __m512i nextcarry;
    __m512i carry = _mm512_set1_epi64(0);
    __m512i highmask = _mm512_set1_epi64(VEC_MAXDIGIT);
    __m512i word;

    if (n > DIGITBITS)
    {
        printf("error, vec_bignum52_mask_lshift_n expects n < %d\n", DIGITBITS);
        exit(0);
    }

    for (i = 0; i < NWORDS; i++)
    {
        word = _mm512_load_epi64(u->data + i * VECLEN);
        nextcarry = _mm512_srli_epi64(word, (DIGITBITS - n));
        _mm512_mask_store_epi64(u->data + i * VECLEN, (__mmask8)wmask,
            _mm512_and_epi64(highmask, _mm512_or_epi64(_mm512_slli_epi64(word, n), carry)));
        carry = nextcarry;
    }

    _mm512_mask_store_epi64(u->data + i * VECLEN, (__mmask8)wmask,
        _mm512_and_epi64(highmask, _mm512_or_epi64(_mm512_slli_epi64(word, n), carry)));

    // return an overflow mask
    return wmask & _mm512_cmp_epi64_mask(carry, _mm512_set1_epi64(0), _MM_CMPINT_GT);
}

void vec_bignum52_mask_rshift_n(vec_bignum_t* u, vec_bignum_t* v, int n, uint32_t wmask)
{
    // return the right shift of bignum u by n bits
    int i;
    uint32_t NWORDS = u->WORDS_ALLOC;
    __m512i nextcarry;
    __m512i carry = _mm512_set1_epi64(0);
    __m512i lowmask;
    __m512i word;
    int wshift = n / 52;
    int bshift = n % 52;

    lowmask = _mm512_set1_epi64((1ULL << (uint64_t)bshift) - 1ULL);

    for (i = 2 * NWORDS - 1; (i - wshift) >= 0; i--)
    {
        word = _mm512_load_epi64(u->data + i * VECLEN);
        nextcarry = _mm512_slli_epi64(_mm512_and_epi64(word, lowmask), (DIGITBITS - bshift));
        _mm512_mask_store_epi64(v->data + (i - wshift) * VECLEN, (__mmask8)wmask,
            _mm512_or_epi64(_mm512_srli_epi64(word, bshift), carry));
        carry = nextcarry;
    }

    return;
}

void vec_bignum52_mask_rshift_1(vec_bignum_t * u, uint32_t wmask)
{
    // return the right shift of bignum u by 1
    int i;
    uint32_t NWORDS = u->WORDS_ALLOC;
    __m512i nextcarry;
    __m512i carry = _mm512_set1_epi64(0);
    __m512i lowmask = _mm512_set1_epi64(1);
    __m512i word;

    //carry = 0;
    //for (i = sb - 1; i >= 0; --i)
    //{
    //    nextcarry = (b->data[i] & mask) << y;
    //    a->data[i] = b->data[i] >> x | carry;
    //    carry = nextcarry;
    //}

    for (i = NWORDS - 1; i >= 0; i--)
    {
        word = _mm512_load_epi64(u->data + i * VECLEN);
        nextcarry = _mm512_slli_epi64(_mm512_and_epi64(word, lowmask), (DIGITBITS - 1));
        _mm512_mask_store_epi64(u->data + i * VECLEN, (__mmask8)wmask,
            _mm512_or_epi64(_mm512_srli_epi64(word, 1), carry));
        carry = nextcarry;
    }

    return;
}

void vec_bignum52_mask_sub(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, uint32_t wmask)
{
    // assumptions:
    // a, b, c are of length VECLEN * NWORDS
    // s1 is of length VECLEN
    // a, b, c, and s1 are aligned
    // a and b are both positive
    // a >= b
    int i;
    uint32_t NWORDS = a->WORDS_ALLOC;
    __mmask8 carry = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;

    if (wmask == 0)
        return;

    // subtract the selected elements ('1' in the mask)
    carry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        avec = _mm512_load_epi64(a->data + i * VECLEN);
        bvec = _mm512_load_epi64(b->data + i * VECLEN);
        cvec = _mm512_sbb_epi52(avec, carry, bvec, &carry);
        _mm512_mask_store_epi64(c->data + i * VECLEN, (__mmask8)wmask,
            _mm512_and_epi64(_mm512_set1_epi64(VEC_MAXDIGIT), cvec));
    }

    if (carry)
    {
        // subtract any final borrows that exceed the size of b.
        _mm512_mask_store_epi64(c->data + i * VECLEN, (__mmask8)wmask & carry, _mm512_set1_epi64(0));
    }

    return;
}

void vec_bignum52_add_1(vec_bignum_t *a, base_t *b, vec_bignum_t *c)
{
    // assumptions:
    // a, c are of length VECLEN * NWORDS
    // b is vector of 1-word base_t's
    // s1 is of length VECLEN
    // a, b, c, and s1 are aligned
    // a and b are both positive
    // a >= b
    int i;
    uint32_t NWORDS = a->WORDS_ALLOC;
    __mmask8 carry = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;

    bvec = _mm512_load_epi64(b);
    carry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        avec = _mm512_load_epi64(a->data + i * VECLEN);
        cvec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c->data + i * VECLEN, 
            _mm512_and_epi64(_mm512_set1_epi64(VEC_MAXDIGIT), cvec));
    }

    _mm512_store_epi64(c->data + i * VECLEN, _mm512_maskz_set1_epi64(carry, 1));
    
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



#endif
