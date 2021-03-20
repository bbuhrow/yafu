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

//#ifdef _MSC_VER
//#define USE_AVX512F
//#endif

#ifdef USE_AVX512F
#include <immintrin.h>

void vec_bignum_mask_rshift_n(vec_bignum_t* u, vec_bignum_t* v, int n, uint32_t wmask);

// ---------------------------------------------------------------------
// emulated instructions
// ---------------------------------------------------------------------
__m512i __inline _mm512_mulhi_epu32(__m512i a, __m512i b)
{
    __m512i t1 = _mm512_shuffle_epi32(a, 0xB1);
    __m512i t2 = _mm512_shuffle_epi32(b, 0xB1);
    __m512i evens = _mm512_mul_epu32(a, b);
    __m512i odds = _mm512_mul_epu32(t1, t2);
    //return _mm512_mask_mov_epi32(_mm512_shuffle_epi32(evens, 0xB1), 0xaaaa, odds);
    return _mm512_mask_mov_epi32(odds, 0x5555, _mm512_shuffle_epi32(evens, 0xB1));
}

__m512i __inline _mm512_mask_adc_epi32(__m512i a, __mmask16 m, __mmask16 c, __m512i b, __mmask16 *cout)
{
    __m512i t = _mm512_add_epi32(a, b);
    *cout = _mm512_cmplt_epu32_mask(t, a);
    __m512i t2 = _mm512_mask_add_epi32(a, m, t, _mm512_maskz_set1_epi32(c, 1));
    *cout = _mm512_kor(*cout, _mm512_mask_cmplt_epu32_mask(m, t2, t));
    return t2;
}

__m512i __inline _mm512_adc_epi32_test1(__m512i a, __mmask16 c, __m512i b, __mmask16 *cout)
{
    __m512i t = _mm512_add_epi32(a, b);
    *cout = _mm512_cmplt_epu32_mask(t, a);
    __m512i t2 = _mm512_add_epi32(t, _mm512_maskz_set1_epi32(c, 1));
    *cout = _mm512_kor(*cout, _mm512_cmplt_epu32_mask(t2, t));
    return t2;
}

__m512i __inline _mm512_adc_epi32_test2(__m512i a, __mmask16 c, __m512i b, __mmask16 *cout)
{
    // looks like a slightly improved data dependency chain... 
    // but it tested slower for 1024-b inputs...
    __m512i t = _mm512_add_epi32(a, b);
    __mmask16 gt0 = _mm512_kor(_mm512_test_epi32_mask(b, b), c);

    t = _mm512_add_epi32(t, _mm512_maskz_set1_epi32(c, 1));
    *cout = _mm512_kand(_mm512_cmple_epu32_mask(t, a), gt0);
    return t;
}

__m512i __inline _mm512_adc_epi32(__m512i a, __mmask16 c, __m512i b, __mmask16 *cout)
{
    __m512i t = _mm512_add_epi32(a, b);
    t = _mm512_add_epi32(t, _mm512_maskz_set1_epi32(c, 1));
    *cout = _mm512_cmplt_epu32_mask(t, a) | (_mm512_cmpeq_epu32_mask(t, a) & c);
    return t;
}

__m512i __inline _mm512_addcarry_epi32(__m512i a, __mmask16 c, __mmask16 *cout)
{
    __m512i t = _mm512_add_epi32(a, _mm512_maskz_set1_epi32(c, 1));
    *cout = _mm512_cmplt_epu32_mask(t, a);
    return t;
}

__m512i __inline _mm512_subborrow_epi32(__m512i a, __mmask16 c, __mmask16 *cout)
{
    __m512i t = _mm512_sub_epi32(a, _mm512_maskz_set1_epi32(c, 1));
    *cout = _mm512_cmpeq_epu32_mask(a, _mm512_setzero_epi32());
    return t;
}

__m512i __inline _mm512_mask_sbb_epi32(__m512i a, __mmask16 m, __mmask16 c, __m512i b, __mmask16 *cout)
{
    __m512i t = _mm512_sub_epi32(a, b);
    *cout = _mm512_cmpgt_epu32_mask(t, a);
    __m512i t2 = _mm512_mask_sub_epi32(a, m, t, _mm512_maskz_set1_epi32(c, 1));
    *cout = _mm512_kor(*cout, _mm512_cmpgt_epu32_mask(t2, t));
    return t2;
}

__m512i __inline _mm512_sbb_epi32(__m512i a, __mmask16 c, __m512i b, __mmask16 *cout)
{
    __m512i t = _mm512_sub_epi32(a, b);
    *cout = _mm512_cmpgt_epu32_mask(t, a);
    __m512i t2 = _mm512_sub_epi32(t, _mm512_maskz_set1_epi32(c, 1));
    *cout = _mm512_kor(*cout, _mm512_cmpgt_epu32_mask(t2, t));
    return t2;
}

__m512i __inline _mm512_sbb_epi64(__m512i a, __mmask8 c, __m512i b, __mmask8 *cout)
{
    __m512i t = _mm512_sub_epi64(a, b);
    *cout = _mm512_cmpgt_epu64_mask(t, a);
    __m512i t2 = _mm512_sub_epi64(t, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_kor(*cout, _mm512_cmpgt_epu64_mask(t2, t));
    return t2;
}

__m512i __inline _mm512_addsetc_epi32(__m512i a, __m512i b, __mmask16 *cout)
{
    __m512i t = _mm512_add_epi32(a, b);
    *cout = _mm512_cmplt_epu32_mask(t, a);
    return t;
}

__m512i __inline _mm512_subsetc_epi32(__m512i a, __m512i b, __mmask16 *cout)
{
    __m512i t = _mm512_sub_epi32(a, b);
    *cout = _mm512_cmpgt_epu32_mask(b, a);
    return t;
}

__inline void _mm512_epi32_to_eo64(__m512i a, __m512i *e64, __m512i *o64)
{
    *e64 = _mm512_maskz_mov_epi32(0x5555, a);
    *o64 = _mm512_maskz_mov_epi32(0x5555, _mm512_shuffle_epi32(a, 0xB1));
    return;
}

__inline __m512i _mm512_eo64lo_to_epi32(__m512i e64, __m512i o64)
{
    return _mm512_mask_blend_epi32(0xAAAA, e64, _mm512_shuffle_epi32(o64, 0xB1));
}

__inline __m512i _mm512_eo64hi_to_epi32(__m512i e64, __m512i o64)
{
    return _mm512_mask_blend_epi32(0xAAAA, _mm512_shuffle_epi32(e64, 0xB1), o64);
}

__inline void _mm512_mul_eo64_epi32(__m512i a, __m512i b, __m512i *e64, __m512i *o64)
{
    // multiply the 16-element 32-bit vectors a and b to produce two 8-element
    // 64-bit vector products e64 and o64, where e64 is the even elements
    // of a*b and o64 is the odd elements of a*b
    //__m512i t1 = _mm512_shuffle_epi32(a, 0xB1);
    //__m512i t2 = _mm512_shuffle_epi32(b, 0xB1);

    //_mm512_shuffle_epi32(a, 0xB1);
    //_mm512_shuffle_epi32(b, 0xB1);
    *e64 = _mm512_mul_epu32(a, b);
    *o64 = _mm512_mul_epu32(_mm512_shuffle_epi32(a, 0xB1), _mm512_shuffle_epi32(b, 0xB1));

    return;
}

#define _mm512_iseven_epi32(x) \
    _mm512_cmp_epi32_mask(_mm512_setzero_epi32(), _mm512_and_epi32((x), _mm512_set1_epi32(1)), _MM_CMPINT_EQ)

#define _mm512_isodd_epi32(x) \
    _mm512_cmp_epi32_mask(_mm512_set1_epi32(1), _mm512_and_epi32((x), _mm512_set1_epi32(1)), _MM_CMPINT_EQ)


#define ACCUM_EO_PROD2(sum_e, sum_o, carry_e, carry_o, in_e, in_o) \
	sum_e = _mm512_add_epi64(sum_e, in_e);	\
	sum_o = _mm512_add_epi64(sum_o, in_o);	\
	scarry_e1 = _mm512_cmplt_epu64_mask(sum_e, in_e);	\
	scarry_o1 = _mm512_cmplt_epu64_mask(sum_o, in_o);	\
	carry_e = _mm512_mask_add_epi64(carry_e, scarry_e1, hiword, carry_e);	\
	carry_o = _mm512_mask_add_epi64(carry_o, scarry_o1, hiword, carry_o);

#define ACCUM_EO_PROD(sum_e, sum_o, carry_e, carry_o) \
	sum_e = _mm512_add_epi64(sum_e, prod1_e);	\
	sum_o = _mm512_add_epi64(sum_o, prod1_o);	\
	scarry_e1 = _mm512_cmplt_epu64_mask(sum_e, prod1_e);	\
	scarry_o1 = _mm512_cmplt_epu64_mask(sum_o, prod1_o);	\
	carry_e = _mm512_mask_add_epi64(carry_e, scarry_e1, hiword, carry_e);	\
	carry_o = _mm512_mask_add_epi64(carry_o, scarry_o1, hiword, carry_o);
#define ACCUM_DOUBLED_EO_PROD(sum_e, sum_o, carry_e, carry_o) \
	sum_e = _mm512_add_epi64(sum_e, prod1_e);	\
	sum_o = _mm512_add_epi64(sum_o, prod1_o);	\
	scarry_e1 = _mm512_cmplt_epu64_mask(sum_e, prod1_e);	\
	scarry_o1 = _mm512_cmplt_epu64_mask(sum_o, prod1_o);	\
	carry_e = _mm512_mask_add_epi64(carry_e, scarry_e1, hiword2, carry_e);	\
	carry_o = _mm512_mask_add_epi64(carry_o, scarry_o1, hiword2, carry_o);



void vecmulmod(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata)
{
    int i, j, k;
    uint32_t NWORDS = n->WORDS_ALLOC;
    uint32_t NBLOCKS = n->WORDS_ALLOC / BLOCKWORDS;
    __m512i a0, a1, a2, a3;
    __m512i b0, b1, b2, b3, b4, b5, b6;
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;
    __m512i to0, to1, to2, to3, to4, to5, to6, to7;

    __m512i acc_e0;
    __m512i acc_o0;

    __m512i acc_e1;
    __m512i acc_o1;
    // 31

    __m512i nhatvec_e = _mm512_load_epi32(mdata->vrho);
    __m512i nhatvec_o = _mm512_shuffle_epi32(nhatvec_e, 0xB1);;

    __m512i prod1_e;
    __m512i prod1_o;

    __m512i hiword = _mm512_set1_epi64(0x000000100000000);
    __m512i zero = _mm512_set1_epi64(0);
    // 37

    __mmask8 scarry_e1 = 0;
    __mmask8 scarry_o1 = 0;
    __mmask8 scarry_e = 0;
    __mmask8 scarry_o = 0;
    __mmask16 scarry2;
    __mmask16 scarry;

    // zero the accumulator
    acc_e0 = acc_o0 = acc_e1 = acc_o1 = zero;


    // first half mul
    for (i = 0; i < NBLOCKS; i++)
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;
        to0 = to1 = to2 = to3 = to4 = to5 = to6 = to7 = zero;

        for (j = i; j > 0; j--)
        {
            a0 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);


            //k == 0;
            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a2, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a3, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            //k == 1;
            _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a3, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            //k == 2;
            _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a2, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a3, b5, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            //k == 3;
            _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a2, b5, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a3, b6, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

        }

        // finish each triangular shaped column sum
        a0 = _mm512_load_epi32(a->data + (i * BLOCKWORDS + 0) * VECLEN);
        a1 = _mm512_load_epi32(a->data + (i * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi32(a->data + (i * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi32(a->data + (i * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi32(b->data + 0 * VECLEN);
        b1 = _mm512_load_epi32(b->data + 1 * VECLEN);
        b2 = _mm512_load_epi32(b->data + 2 * VECLEN);
        b3 = _mm512_load_epi32(b->data + 3 * VECLEN);

        // ======
        _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        // ======
        _mm512_mul_eo64_epi32(a1, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te2, to2, te3, to3);

        _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te2, to2, te3, to3);

        // ======
        _mm512_mul_eo64_epi32(a2, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te4, to4, te5, to5);

        _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te4, to4, te5, to5);

        _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te4, to4, te5, to5);

        // ======
        _mm512_mul_eo64_epi32(a3, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te6, to6, te7, to7);

        _mm512_mul_eo64_epi32(a2, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te6, to6, te7, to7);

        _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te6, to6, te7, to7);

        _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te6, to6, te7, to7);

        // for those 's' we have already accumulated, compute the
       // block s*n accumulations
        for (j = i; j > 0; j--)
        {
            a0 = _mm512_load_epi32(s->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi32(s->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi32(s->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi32(s->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            //k == 0;
            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a2, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a3, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            //k == 1;
            _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a3, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            //k == 2;
            _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a2, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a3, b5, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            //k == 3;
            _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a2, b5, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a3, b6, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

        }

        // now, column by column, add in the s*n contribution and reduce to 
       // a single 64+x bit accumulator while storing the intermediate product
       // 's' as we go.
        j = 0;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te0, to0);
        acc_e1 = _mm512_add_epi64(acc_e1, te1);
        acc_o1 = _mm512_add_epi64(acc_o1, to1);

        prod1_e = _mm512_mul_epu32(nhatvec_e, acc_e0);
        prod1_o = _mm512_mul_epu32(nhatvec_o, acc_o0);
        a0 = _mm512_eo64lo_to_epi32(prod1_e, prod1_o);

        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        b0 = _mm512_load_epi32(n->data + 0 * VECLEN);
        _mm512_mul_eo64_epi32(b0, a0, &prod1_e, &prod1_o);

        // add in the final product
        ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 1;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te2, to2);
        acc_e1 = _mm512_add_epi64(acc_e1, te3);
        acc_o1 = _mm512_add_epi64(acc_o1, to3);

        for (k = 0; k < j; k++)
        {
            a0 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + k) * VECLEN);
            b0 = _mm512_load_epi32(n->data + (j - k) * VECLEN);

            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);
        }

        prod1_e = _mm512_mul_epu32(nhatvec_e, acc_e0);
        prod1_o = _mm512_mul_epu32(nhatvec_o, acc_o0);
        a0 = _mm512_eo64lo_to_epi32(prod1_e, prod1_o);

        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        b0 = _mm512_load_epi32(n->data + 0 * VECLEN);
        _mm512_mul_eo64_epi32(b0, a0, &prod1_e, &prod1_o);

        // add in the final product
        ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 2;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te4, to4);
        acc_e1 = _mm512_add_epi64(acc_e1, te5);
        acc_o1 = _mm512_add_epi64(acc_o1, to5);

        for (k = 0; k < j; k++)
        {
            a0 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + k) * VECLEN);
            b0 = _mm512_load_epi32(n->data + (j - k) * VECLEN);

            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);
        }

        prod1_e = _mm512_mul_epu32(nhatvec_e, acc_e0);
        prod1_o = _mm512_mul_epu32(nhatvec_o, acc_o0);
        a0 = _mm512_eo64lo_to_epi32(prod1_e, prod1_o);

        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        b0 = _mm512_load_epi32(n->data + 0 * VECLEN);
        _mm512_mul_eo64_epi32(b0, a0, &prod1_e, &prod1_o);

        // add in the final product
        ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 3;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te6, to6);
        acc_e1 = _mm512_add_epi64(acc_e1, te7);
        acc_o1 = _mm512_add_epi64(acc_o1, to7);

        for (k = 0; k < j; k++)
        {
            a0 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + k) * VECLEN);
            b0 = _mm512_load_epi32(n->data + (j - k) * VECLEN);

            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);
        }

        prod1_e = _mm512_mul_epu32(nhatvec_e, acc_e0);
        prod1_o = _mm512_mul_epu32(nhatvec_o, acc_o0);
        a0 = _mm512_eo64lo_to_epi32(prod1_e, prod1_o);

        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        b0 = _mm512_load_epi32(n->data + 0 * VECLEN);
        _mm512_mul_eo64_epi32(b0, a0, &prod1_e, &prod1_o);

        // add in the final product
        ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

    }

    // second half mul
    for (i = NBLOCKS; i < 2 * NBLOCKS; i++)
    {

        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;
        to0 = to1 = to2 = to3 = to4 = to5 = to6 = to7 = zero;

        for (j = i - NBLOCKS + 1; j < NBLOCKS; j++)
        {
            a0 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            // accumulate a * b
            //k == 0;
            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a2, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a3, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            //k == 1;
            _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a3, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            //k == 2;
            _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a2, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a3, b5, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            //k == 3;
            _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a2, b5, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a3, b6, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);


            // accumulate s * n
            a0 = _mm512_load_epi32(s->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi32(s->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi32(s->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi32(s->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            //k == 0;
            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a2, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a3, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            //k == 1;
            _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a3, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            //k == 2;
            _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a2, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a3, b5, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            //k == 3;
            _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a2, b5, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a3, b6, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

        }

        // finish each triangular shaped column sum (a * b)
        a1 = _mm512_load_epi32(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi32(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi32(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi32(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi32(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi32(b->data + (NWORDS - 3) * VECLEN);

        // ======
        _mm512_mul_eo64_epi32(a1, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        _mm512_mul_eo64_epi32(a2, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        _mm512_mul_eo64_epi32(a3, b2, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        // ======
        _mm512_mul_eo64_epi32(a2, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te2, to2, te3, to3);

        _mm512_mul_eo64_epi32(a3, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te2, to2, te3, to3);

        // ======
        _mm512_mul_eo64_epi32(a3, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te4, to4, te5, to5);

        // finish each triangular shaped column sum (s * n)
        a1 = _mm512_load_epi32(s->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi32(s->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi32(s->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi32(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi32(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi32(n->data + (NWORDS - 3) * VECLEN);

        // ======
        _mm512_mul_eo64_epi32(a1, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        _mm512_mul_eo64_epi32(a2, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        _mm512_mul_eo64_epi32(a3, b2, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        // ======
        _mm512_mul_eo64_epi32(a2, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te2, to2, te3, to3);

        _mm512_mul_eo64_epi32(a3, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te2, to2, te3, to3);

        // ======
        _mm512_mul_eo64_epi32(a3, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te4, to4, te5, to5);

        j = 0;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te0, to0);
        acc_e1 = _mm512_add_epi64(acc_e1, te1);
        acc_o1 = _mm512_add_epi64(acc_o1, to1);

        // store the low-word final result
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + ((i - NBLOCKS) * BLOCKWORDS + j) * VECLEN, a0);

        // and shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 1;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te2, to2);
        acc_e1 = _mm512_add_epi64(acc_e1, te3);
        acc_o1 = _mm512_add_epi64(acc_o1, to3);

        // store the low-word final result
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + ((i - NBLOCKS) * BLOCKWORDS + j) * VECLEN, a0);

        // and shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 2;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te4, to4);
        acc_e1 = _mm512_add_epi64(acc_e1, te5);
        acc_o1 = _mm512_add_epi64(acc_o1, to5);

        // store the low-word final result
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + ((i - NBLOCKS) * BLOCKWORDS + j) * VECLEN, a0);

        // and shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 3;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te6, to6);
        acc_e1 = _mm512_add_epi64(acc_e1, te7);
        acc_o1 = _mm512_add_epi64(acc_o1, to7);

        // store the low-word final result
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + ((i - NBLOCKS) * BLOCKWORDS + j) * VECLEN, a0);

        // and shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;
    }

    a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
    scarry2 = _mm512_cmp_epu32_mask(a0, zero, _MM_CMPINT_EQ);

    // subtract n from tmp
    scarry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        a1 = _mm512_load_epi32(s->data + i * VECLEN);
        b0 = _mm512_load_epi32(n->data + i * VECLEN);
        a0 = _mm512_sbb_epi32(a1, scarry, b0, &scarry);
        _mm512_store_epi32(c->data + i * VECLEN, a0);
    }

    // negate any final borrows if there was also a final carry.
    scarry &= scarry2;

    // if there was a final borrow, we didn't need to do the subtraction after all.
    // replace with original results based on final borrow mask.
    for (i = NWORDS - 1; i >= 0; i--)
    {
        b0 = _mm512_load_epi32(s->data + i * VECLEN);
        _mm512_mask_store_epi32(c->data + i * VECLEN, scarry, b0);
    }

    c->size = NWORDS;
    return;
}

void vecsqrmod(vec_bignum_t *a, vec_bignum_t *c, vec_bignum_t *n, vec_bignum_t *s, vec_monty_t*mdata)
{
    int i, j, k;
    vec_bignum_t *b = a;
    uint32_t NWORDS = n->WORDS_ALLOC;
    uint32_t NBLOCKS = n->WORDS_ALLOC / BLOCKWORDS;
    __m512i a0, a1, a2, a3;
    __m512i b0, b1, b2, b3, b4, b5, b6;
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;
    __m512i to0, to1, to2, to3, to4, to5, to6, to7;

    __m512i acc_e0;
    __m512i acc_o0;

    __m512i acc_e1;
    __m512i acc_o1;
    // 31

    __m512i nhatvec_e = _mm512_load_epi32(mdata->vrho); // _mm512_set1_epi32(nhat);
    __m512i nhatvec_o = _mm512_shuffle_epi32(nhatvec_e, 0xB1);;

    __m512i prod1_e;
    __m512i prod1_o;

    __m512i hiword = _mm512_set1_epi64(0x000000100000000);
    __m512i hiword2 = _mm512_set1_epi64(0x000000200000000);
    __m512i zero = _mm512_set1_epi64(0);
    // 37

    __mmask8 scarry_e1 = 0;
    __mmask8 scarry_o1 = 0;
    __mmask8 scarry_e = 0;
    __mmask8 scarry_o = 0;
    __mmask16 scarry2;
    __mmask16 scarry;

    // zero the accumulator
    acc_e0 = acc_o0 = acc_e1 = acc_o1 = zero;


    // first half sqr
    for (i = 0; i < NBLOCKS; i++)
    {

        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;
        to0 = to1 = to2 = to3 = to4 = to5 = to6 = to7 = zero;

        if (i & 1)
        {
            __mmask8 scarry_e1 = 0;
            __mmask8 scarry_o1 = 0;

            // i odd
            for (j = 0; j < (i - 1) / 2; j++)
            {
                // for 384-bit inputs NBLOCKS=3 and this loop doesn't run at all.
                // i=0 no, even
                // i=1 no
                // i=2 no, even


                // for 1024-bit inputs NBLOCKS=8 and this loop runs 6 times over all i.
                // i=0 no, even
                // i=1 no
                // i=2 no, even
                // i=3 once
                // i=4 no, even
                // i=5 twice
                // i=6 no, even
                // i=7 thrice
                // with the doubling trick we trade 96 instructions for each j-loop iteration
                // and the 72 instructions for each odd i, after the j-loop,
                // for the 32 instructions for each odd i.  That saves 6*96+4*72-4*32=736 instructions.

                //hips_block_mul_type3(a->data + ((j + 1) * BLOCKWORDS) * VECLEN,
                //    a->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS) * VECLEN, t_e, t_o);

                a0 = _mm512_load_epi32(a->data + ((j + 1) * BLOCKWORDS - 1) * VECLEN);
                a1 = _mm512_load_epi32(a->data + ((j + 1) * BLOCKWORDS - 2) * VECLEN);
                a2 = _mm512_load_epi32(a->data + ((j + 1) * BLOCKWORDS - 3) * VECLEN);
                a3 = _mm512_load_epi32(a->data + ((j + 1) * BLOCKWORDS - 4) * VECLEN);

                b0 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS + 1) * VECLEN);
                b1 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS + 2) * VECLEN);
                b2 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS + 3) * VECLEN);
                b3 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS + 4) * VECLEN);
                b4 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS + 5) * VECLEN);
                b5 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS + 6) * VECLEN);
                b6 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS + 7) * VECLEN);

                // save independent sum/carry words for each product-column in the block.
                // uses 11 register inputs, 16 register outputs, and 3 aux vectors.
                // since all terms in this loop are doubled, we do the doubling
                // after the loop with a left shift.

                //k == 0;
                _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);          // a-1, b+1
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);          // a-2, b+2
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a2, b2, &prod1_e, &prod1_o);          // a-3, b+3
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a3, b3, &prod1_e, &prod1_o);          // a-4, b+4
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                //k == 1;
                _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);          // a-1, b+2
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);          // a-2, b+3
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);          // a-3, b+4
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                _mm512_mul_eo64_epi32(a3, b4, &prod1_e, &prod1_o);          // a-4, b+5
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                //k == 2;
                _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);          // a-1, b+3
                ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

                _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);          // a-2, b+4
                ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

                _mm512_mul_eo64_epi32(a2, b4, &prod1_e, &prod1_o);          // a-3, b+5
                ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

                _mm512_mul_eo64_epi32(a3, b5, &prod1_e, &prod1_o);          // a-4, b+6
                ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

                //k == 3;
                _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);          // a-1, b+4
                ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

                _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);          // a-2, b+5
                ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

                _mm512_mul_eo64_epi32(a2, b5, &prod1_e, &prod1_o);          // a-3, b+6
                ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

                _mm512_mul_eo64_epi32(a3, b6, &prod1_e, &prod1_o);          // a-4, b+7
                ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            }

            // for 384-bit inputs when i == 1, j = 0, a = {3,2,1,0} and b = {2,3,4,5,6,7}
                        // for 512-bit inputs when i == 3, j = 1, a = {7,6,5,4} and b = {6,7,8,9,a,b}
                        // for 512-bit inputs when i == 1, j = 0, a = {3,2,1,0} and b = {2,3,4,5,6,7}

            a0 = _mm512_load_epi32(a->data + (i * BLOCKWORDS - j * BLOCKWORDS - 1) * VECLEN);
            a1 = _mm512_load_epi32(a->data + (i * BLOCKWORDS - j * BLOCKWORDS - 2) * VECLEN);
            a2 = _mm512_load_epi32(a->data + (i * BLOCKWORDS - j * BLOCKWORDS - 3) * VECLEN);
            a3 = _mm512_load_epi32(a->data + (i * BLOCKWORDS - j * BLOCKWORDS - 4) * VECLEN);

            b1 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + 7) * VECLEN);

            // save independent sum/carry words for each product-column in the block.
            // uses 11 register inputs, 16 register outputs, and 3 aux vectors.
            //k == 0;
            _mm512_mul_eo64_epi32(a2, b2, &prod1_e, &prod1_o);
            //te0 = _mm512_addsetc_epi64(te0, prod1_e, &scarry_e);
            //to0 = _mm512_addsetc_epi64(to0, prod1_o, &scarry_o);
            //te1 = _mm512_mask_add_epi64(te1, scarry_e, hiword2, te1);
            //to1 = _mm512_mask_add_epi64(to1, scarry_o, hiword2, to1);
            ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a3, b3, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

            //k == 1;
            _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);
            //te2 = _mm512_addsetc_epi64(te2, prod1_e, &scarry_e);
            //to2 = _mm512_addsetc_epi64(to2, prod1_o, &scarry_o);
            //te3 = _mm512_mask_add_epi64(te3, scarry_e, hiword2, te3);
            //to3 = _mm512_mask_add_epi64(to3, scarry_o, hiword2, to3);

            _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a3, b4, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            //k == 2;
            _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);
            //te4 = _mm512_addsetc_epi64(te4, prod1_e, &scarry_e);
            //to4 = _mm512_addsetc_epi64(to4, prod1_o, &scarry_o);
            //te5 = _mm512_mask_add_epi64(te5, scarry_e, hiword2, te5);
            //to5 = _mm512_mask_add_epi64(to5, scarry_o, hiword2, to5);
            ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a2, b4, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a3, b5, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

            //k == 3;
            _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);
            //te6 = _mm512_addsetc_epi64(te6, prod1_e, &scarry_e);
            //to6 = _mm512_addsetc_epi64(to6, prod1_o, &scarry_o);
            //te7 = _mm512_mask_add_epi64(te7, scarry_e, hiword2, te7);
            //to7 = _mm512_mask_add_epi64(to7, scarry_o, hiword2, to7);

            _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a2, b5, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a3, b6, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            // all terms so far need to be doubled.  Do that all at once with these
            // left shifts.  We need the high-bit of the low word to end up in the
            // 32nd bit position (because the high words are offset by that much
            // for easier combining later on).  _mm512_maskz_srli_epi32 allows us
            // to do the right shift simultaneously with clearing out the low bits.
            te1 = _mm512_or_epi64(te1, _mm512_maskz_srli_epi32(0xaaaa, te0, 31));
            to1 = _mm512_or_epi64(to1, _mm512_maskz_srli_epi32(0xaaaa, to0, 31));
            te0 = _mm512_slli_epi64(te0, 1);
            to0 = _mm512_slli_epi64(to0, 1);

            te3 = _mm512_or_epi64(te3, _mm512_maskz_srli_epi32(0xaaaa, te2, 31));
            to3 = _mm512_or_epi64(to3, _mm512_maskz_srli_epi32(0xaaaa, to2, 31));
            te2 = _mm512_slli_epi64(te2, 1);
            to2 = _mm512_slli_epi64(to2, 1);

            te5 = _mm512_or_epi64(te5, _mm512_maskz_srli_epi32(0xaaaa, te4, 31));
            to5 = _mm512_or_epi64(to5, _mm512_maskz_srli_epi32(0xaaaa, to4, 31));
            te4 = _mm512_slli_epi64(te4, 1);
            to4 = _mm512_slli_epi64(to4, 1);

            te7 = _mm512_or_epi64(te7, _mm512_maskz_srli_epi32(0xaaaa, te6, 31));
            to7 = _mm512_or_epi64(to7, _mm512_maskz_srli_epi32(0xaaaa, to6, 31));
            te6 = _mm512_slli_epi64(te6, 1);
            to6 = _mm512_slli_epi64(to6, 1);

            // finally the two non-doubled terms.
            _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);
            //te0 = _mm512_addsetc_epi64(te0, prod1_e, &scarry_e);
            //to0 = _mm512_addsetc_epi64(to0, prod1_o, &scarry_o);
            //te1 = _mm512_mask_add_epi64(te1, scarry_e, hiword, te1);
            //to1 = _mm512_mask_add_epi64(to1, scarry_o, hiword, to1);

            _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);
            //te4 = _mm512_addsetc_epi64(te4, prod1_e, &scarry_e);
            //to4 = _mm512_addsetc_epi64(to4, prod1_o, &scarry_o);
            //te5 = _mm512_mask_add_epi64(te5, scarry_e, hiword, te5);
            //to5 = _mm512_mask_add_epi64(to5, scarry_o, hiword, to5);

        }
        else
        {
            __mmask8 scarry_e1 = 0;
            __mmask8 scarry_o1 = 0;

            // i even
            for (j = 0; j < i / 2; j++)
            {
                // for 384-bit inputs NBLOCKS=3 and this loop runs once. 
                // i=0: 0 times
                // i=1: no, odd
                // i=2: 1 time

                // for 1024-bit inputs NBLOCKS=8 and this loop runs 6 times over all i.
                // i=0 0 times
                // i=1 no, odd
                // i=2 1 time
                // i=3 no, odd
                // i=4 2 times
                // i=5 no, odd
                // i=6 3 times
                // i=7 no, odd
                // with the doubling trick we trade 96 instructions for each j-loop iteration
                // and the 24 instructions for each even i, after the j-loop,
                // for the 40 instructions once per even i.  That saves 6*96+4*24-4*40=512 instructions.


                //hips_block_mul_type3(a->data + ((j + 1) * BLOCKWORDS) * VECLEN,
                //    a->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS) * VECLEN, t_e, t_o);

                // when i = 2, j = 0, a = {3,2,1,0}, b = {5,6,7,8,9,a,b}
                a0 = _mm512_load_epi32(a->data + ((j + 1) * BLOCKWORDS - 1) * VECLEN);
                a1 = _mm512_load_epi32(a->data + ((j + 1) * BLOCKWORDS - 2) * VECLEN);
                a2 = _mm512_load_epi32(a->data + ((j + 1) * BLOCKWORDS - 3) * VECLEN);
                a3 = _mm512_load_epi32(a->data + ((j + 1) * BLOCKWORDS - 4) * VECLEN);

                b0 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j
                    * BLOCKWORDS + 1) * VECLEN);
                b1 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j
                    * BLOCKWORDS + 2) * VECLEN);
                b2 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j
                    * BLOCKWORDS + 3) * VECLEN);
                b3 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j
                    * BLOCKWORDS + 4) * VECLEN);
                b4 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j
                    * BLOCKWORDS + 5) * VECLEN);
                b5 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j
                    * BLOCKWORDS + 6) * VECLEN);
                b6 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j
                    * BLOCKWORDS + 7) * VECLEN);

                // save independent sum/carry words for each product-column in the block.
                // uses 11 register inputs, 16 register outputs, and 3 aux vectors.
                //k == 0;
                _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);          // a-1, b+1
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);          // a-2, b+2
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a2, b2, &prod1_e, &prod1_o);          // a-3, b+3
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a3, b3, &prod1_e, &prod1_o);          // a-4, b+4
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                //k == 1;
                _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);          // a-1, b+2
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);          // a-2, b+3
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);          // a-3, b+4
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                _mm512_mul_eo64_epi32(a3, b4, &prod1_e, &prod1_o);          // a-4, b+5
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                //k == 2;
                _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);          // a-1, b+3
                ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

                _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);          // a-2, b+4
                ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

                _mm512_mul_eo64_epi32(a2, b4, &prod1_e, &prod1_o);          // a-3, b+5
                ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

                _mm512_mul_eo64_epi32(a3, b5, &prod1_e, &prod1_o);          // a-4, b+6
                ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

                //k == 3;
                _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);          // a-1, b+4
                ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

                _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);          // a-2, b+5
                ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

                _mm512_mul_eo64_epi32(a2, b5, &prod1_e, &prod1_o);          // a-3, b+6
                ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

                _mm512_mul_eo64_epi32(a3, b6, &prod1_e, &prod1_o);          // a-4, b+7
                ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);
            }

            a0 = _mm512_load_epi32(a->data + (i / 2 * BLOCKWORDS + 0) * VECLEN);
            a1 = _mm512_load_epi32(a->data + (i / 2 * BLOCKWORDS + 1) * VECLEN);
            a2 = _mm512_load_epi32(a->data + (i / 2 * BLOCKWORDS + 2) * VECLEN);
            a3 = _mm512_load_epi32(a->data + (i / 2 * BLOCKWORDS + 3) * VECLEN);

            // save independent sum/carry words for each product-column in the block.
            // uses 11 register inputs, 16 register outputs, and 3 aux vectors.

            //k == 1;
            _mm512_mul_eo64_epi32(a0, a1, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            //k == 2;
            _mm512_mul_eo64_epi32(a0, a2, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

            //k == 3;
            _mm512_mul_eo64_epi32(a0, a3, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a1, a2, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            // all terms so far need to be doubled.  Do that all at once with these
            // left shifts.
            te1 = _mm512_or_epi64(te1, _mm512_maskz_srli_epi32(0xaaaa, te0, 31));
            to1 = _mm512_or_epi64(to1, _mm512_maskz_srli_epi32(0xaaaa, to0, 31));
            te0 = _mm512_slli_epi64(te0, 1);
            to0 = _mm512_slli_epi64(to0, 1);

            te3 = _mm512_or_epi64(te3, _mm512_maskz_srli_epi32(0xaaaa, te2, 31));
            to3 = _mm512_or_epi64(to3, _mm512_maskz_srli_epi32(0xaaaa, to2, 31));
            te2 = _mm512_slli_epi64(te2, 1);
            to2 = _mm512_slli_epi64(to2, 1);

            te5 = _mm512_or_epi64(te5, _mm512_maskz_srli_epi32(0xaaaa, te4, 31));
            to5 = _mm512_or_epi64(to5, _mm512_maskz_srli_epi32(0xaaaa, to4, 31));
            te4 = _mm512_slli_epi64(te4, 1);
            to4 = _mm512_slli_epi64(to4, 1);

            te7 = _mm512_or_epi64(te7, _mm512_maskz_srli_epi32(0xaaaa, te6, 31));
            to7 = _mm512_or_epi64(to7, _mm512_maskz_srli_epi32(0xaaaa, to6, 31));
            te6 = _mm512_slli_epi64(te6, 1);
            to6 = _mm512_slli_epi64(to6, 1);

            // finally the two non-doubled terms.
            _mm512_mul_eo64_epi32(a0, a0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a1, a1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

        }

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        for (j = 0; j < i; j++)
        {
            __mmask8 scarry_e1 = 0;
            __mmask8 scarry_o1 = 0;

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

            //k == 0;
            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a2, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a3, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            //k == 1;
            _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a3, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            //k == 2;
            _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a2, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a3, b5, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            //k == 3;
            _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a2, b5, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a3, b6, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

        }		// now, column by column, add in the s*n contribution and reduce to 
        // a single 64+x bit accumulator while storing the intermediate product
        // 's' as we go.
        j = 0;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te0, to0);
        acc_e1 = _mm512_add_epi64(acc_e1, te1);
        acc_o1 = _mm512_add_epi64(acc_o1, to1);

        prod1_e = _mm512_mul_epu32(nhatvec_e, acc_e0);
        prod1_o = _mm512_mul_epu32(nhatvec_o, acc_o0);
        a0 = _mm512_eo64lo_to_epi32(prod1_e, prod1_o);

        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        b0 = _mm512_load_epi32(n->data + 0 * VECLEN);
        _mm512_mul_eo64_epi32(b0, a0, &prod1_e, &prod1_o);

        // add in the final product
        ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 1;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te2, to2);
        acc_e1 = _mm512_add_epi64(acc_e1, te3);
        acc_o1 = _mm512_add_epi64(acc_o1, to3);

        for (k = 0; k < j; k++)
        {
            a0 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + k) * VECLEN);
            b0 = _mm512_load_epi32(n->data + (j - k) * VECLEN);

            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);
        }

        prod1_e = _mm512_mul_epu32(nhatvec_e, acc_e0);
        prod1_o = _mm512_mul_epu32(nhatvec_o, acc_o0);
        a0 = _mm512_eo64lo_to_epi32(prod1_e, prod1_o);

        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        b0 = _mm512_load_epi32(n->data + 0 * VECLEN);
        _mm512_mul_eo64_epi32(b0, a0, &prod1_e, &prod1_o);

        // add in the final product
        ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 2;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te4, to4);
        acc_e1 = _mm512_add_epi64(acc_e1, te5);
        acc_o1 = _mm512_add_epi64(acc_o1, to5);

        for (k = 0; k < j; k++)
        {
            a0 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + k) * VECLEN);
            b0 = _mm512_load_epi32(n->data + (j - k) * VECLEN);

            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);
        }

        prod1_e = _mm512_mul_epu32(nhatvec_e, acc_e0);
        prod1_o = _mm512_mul_epu32(nhatvec_o, acc_o0);
        a0 = _mm512_eo64lo_to_epi32(prod1_e, prod1_o);

        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        b0 = _mm512_load_epi32(n->data + 0 * VECLEN);
        _mm512_mul_eo64_epi32(b0, a0, &prod1_e, &prod1_o);

        // add in the final product
        ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 3;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te6, to6);
        acc_e1 = _mm512_add_epi64(acc_e1, te7);
        acc_o1 = _mm512_add_epi64(acc_o1, to7);

        for (k = 0; k < j; k++)
        {
            a0 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + k) * VECLEN);
            b0 = _mm512_load_epi32(n->data + (j - k) * VECLEN);

            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);
        }

        prod1_e = _mm512_mul_epu32(nhatvec_e, acc_e0);
        prod1_o = _mm512_mul_epu32(nhatvec_o, acc_o0);
        a0 = _mm512_eo64lo_to_epi32(prod1_e, prod1_o);

        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        b0 = _mm512_load_epi32(n->data + 0 * VECLEN);
        _mm512_mul_eo64_epi32(b0, a0, &prod1_e, &prod1_o);

        // add in the final product
        ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

    }

    // second half sqr
    for (i = 0; i < NBLOCKS; i++)
    {


        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;
        to0 = to1 = to2 = to3 = to4 = to5 = to6 = to7 = zero;


        for (j = 0; j < (NBLOCKS - i - 1) / 2; j++)
        {
            __mmask8 scarry_e1 = 0;
            __mmask8 scarry_o1 = 0;

            // Compute a solid block (all matching terms are in the lower
            // half triangle of the expansion).

            //hips_block_mul_type3(a->data + (i * BLOCKWORDS + (j + 2) * BLOCKWORDS) * VECLEN,
            //    a->data + (NWORDS - 2 * BLOCKWORDS - j * BLOCKWORDS) * VECLEN, t_e, t_o);
            a0 = _mm512_load_epi32(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
            a1 = _mm512_load_epi32(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
            a2 = _mm512_load_epi32(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);
            a3 = _mm512_load_epi32(a->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);

            b0 = _mm512_load_epi32(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 7) * VECLEN);

            // save independent sum/carry words for each product-column in the block.
            // uses 11 register inputs, 16 register outputs, and 3 aux vectors.
            //k == 0;
            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);          // a-1, b+1
            ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);          // a-2, b+2
            ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a2, b2, &prod1_e, &prod1_o);          // a-3, b+3
            ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a3, b3, &prod1_e, &prod1_o);          // a-4, b+4
            ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

            //k == 1;
            _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);          // a-1, b+2
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);          // a-2, b+3
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);          // a-3, b+4
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a3, b4, &prod1_e, &prod1_o);          // a-4, b+5
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            //k == 2;
            _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);          // a-1, b+3
            ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);          // a-2, b+4
            ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a2, b4, &prod1_e, &prod1_o);          // a-3, b+5
            ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a3, b5, &prod1_e, &prod1_o);          // a-4, b+6
            ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

            //k == 3;
            _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);          // a-1, b+4
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);          // a-2, b+5
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a2, b5, &prod1_e, &prod1_o);          // a-3, b+6
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a3, b6, &prod1_e, &prod1_o);          // a-4, b+7
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);
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
                a0 = _mm512_load_epi32(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
                a1 = _mm512_load_epi32(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
                a2 = _mm512_load_epi32(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);
                a3 = _mm512_load_epi32(a->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);

                b0 = _mm512_load_epi32(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 1) * VECLEN);
                b1 = _mm512_load_epi32(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 2) * VECLEN);
                b2 = _mm512_load_epi32(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 3) * VECLEN);
                b3 = _mm512_load_epi32(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 4) * VECLEN);
                b4 = _mm512_load_epi32(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 5) * VECLEN);

                // save independent sum/carry words for each product-column in the block.
                // uses 11 register inputs, 16 register outputs, and 3 aux vectors.
                //k == 0;
                _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a2, b2, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                //k == 1;
                _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                //k == 2;
                _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

                _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

                //k == 3;
                _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

                _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

                // all terms so far need to be doubled.  Do that all at once with these
                // left shifts.
                te1 = _mm512_or_epi64(te1, _mm512_maskz_srli_epi32(0xaaaa, te0, 31));
                to1 = _mm512_or_epi64(to1, _mm512_maskz_srli_epi32(0xaaaa, to0, 31));
                te0 = _mm512_slli_epi64(te0, 1);
                to0 = _mm512_slli_epi64(to0, 1);

                te3 = _mm512_or_epi64(te3, _mm512_maskz_srli_epi32(0xaaaa, te2, 31));
                to3 = _mm512_or_epi64(to3, _mm512_maskz_srli_epi32(0xaaaa, to2, 31));
                te2 = _mm512_slli_epi64(te2, 1);
                to2 = _mm512_slli_epi64(to2, 1);

                te5 = _mm512_or_epi64(te5, _mm512_maskz_srli_epi32(0xaaaa, te4, 31));
                to5 = _mm512_or_epi64(to5, _mm512_maskz_srli_epi32(0xaaaa, to4, 31));
                te4 = _mm512_slli_epi64(te4, 1);
                to4 = _mm512_slli_epi64(to4, 1);

                te7 = _mm512_or_epi64(te7, _mm512_maskz_srli_epi32(0xaaaa, te6, 31));
                to7 = _mm512_or_epi64(to7, _mm512_maskz_srli_epi32(0xaaaa, to6, 31));
                te6 = _mm512_slli_epi64(te6, 1);
                to6 = _mm512_slli_epi64(to6, 1);

                // finally the two non-doubled terms.
                _mm512_mul_eo64_epi32(a3, a3, &prod1_e, &prod1_o);
                ACCUM_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a2, a2, &prod1_e, &prod1_o);
                ACCUM_EO_PROD(te4, to4, te5, to5);
            }
            else
            {
                // i even, block shape 1.
                // always a continuation of the full-block loop, so use the same 
                // loading pattern.  Only now we don't need as many b-terms.
                a0 = _mm512_load_epi32(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
                a1 = _mm512_load_epi32(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
                a2 = _mm512_load_epi32(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);

                //k == 0;
                _mm512_mul_eo64_epi32(a0, a2, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                //k == 1;
                _mm512_mul_eo64_epi32(a0, a1, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                te1 = _mm512_or_epi64(te1, _mm512_maskz_srli_epi32(0xaaaa, te0, 31));
                to1 = _mm512_or_epi64(to1, _mm512_maskz_srli_epi32(0xaaaa, to0, 31));
                te0 = _mm512_slli_epi64(te0, 1);
                to0 = _mm512_slli_epi64(to0, 1);

                te3 = _mm512_or_epi64(te3, _mm512_maskz_srli_epi32(0xaaaa, te2, 31));
                to3 = _mm512_or_epi64(to3, _mm512_maskz_srli_epi32(0xaaaa, to2, 31));
                te2 = _mm512_slli_epi64(te2, 1);
                to2 = _mm512_slli_epi64(to2, 1);

                // technically only have to do these two if j > 0 (so that
                // they are non-zero from full-block loop iterations).
                te5 = _mm512_or_epi64(te5, _mm512_maskz_srli_epi32(0xaaaa, te4, 31));
                to5 = _mm512_or_epi64(to5, _mm512_maskz_srli_epi32(0xaaaa, to4, 31));
                te4 = _mm512_slli_epi64(te4, 1);
                to4 = _mm512_slli_epi64(to4, 1);

                te7 = _mm512_or_epi64(te7, _mm512_maskz_srli_epi32(0xaaaa, te6, 31));
                to7 = _mm512_or_epi64(to7, _mm512_maskz_srli_epi32(0xaaaa, to6, 31));
                te6 = _mm512_slli_epi64(te6, 1);
                to6 = _mm512_slli_epi64(to6, 1);

                // finally the two non-doubled terms.
                _mm512_mul_eo64_epi32(a1, a1, &prod1_e, &prod1_o);
                ACCUM_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a0, a0, &prod1_e, &prod1_o);
                ACCUM_EO_PROD(te4, to4, te5, to5);
            }
        }

        else
        {
            // NBLOCKS is even
            // i odd, block shape 1.
            if (i & 1)
            {
                // always a continuation of the full-block loop, so use the same 
                // loading pattern.  Only now we don't need as many b-terms.
                a0 = _mm512_load_epi32(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);  // {f, b}
                a1 = _mm512_load_epi32(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);  // {e, a}
                a2 = _mm512_load_epi32(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);  // {d, 9}

                //k == 0;
                _mm512_mul_eo64_epi32(a0, a2, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                //k == 1;
                _mm512_mul_eo64_epi32(a0, a1, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                te1 = _mm512_or_epi64(te1, _mm512_maskz_srli_epi32(0xaaaa, te0, 31));
                to1 = _mm512_or_epi64(to1, _mm512_maskz_srli_epi32(0xaaaa, to0, 31));
                te0 = _mm512_slli_epi64(te0, 1);
                to0 = _mm512_slli_epi64(to0, 1);

                te3 = _mm512_or_epi64(te3, _mm512_maskz_srli_epi32(0xaaaa, te2, 31));
                to3 = _mm512_or_epi64(to3, _mm512_maskz_srli_epi32(0xaaaa, to2, 31));
                te2 = _mm512_slli_epi64(te2, 1);
                to2 = _mm512_slli_epi64(to2, 1);

                // technically only have to do these two if j > 0 (so that
                // they are non-zero from full-block loop iterations).
                te5 = _mm512_or_epi64(te5, _mm512_maskz_srli_epi32(0xaaaa, te4, 31));
                to5 = _mm512_or_epi64(to5, _mm512_maskz_srli_epi32(0xaaaa, to4, 31));
                te4 = _mm512_slli_epi64(te4, 1);
                to4 = _mm512_slli_epi64(to4, 1);

                te7 = _mm512_or_epi64(te7, _mm512_maskz_srli_epi32(0xaaaa, te6, 31));
                to7 = _mm512_or_epi64(to7, _mm512_maskz_srli_epi32(0xaaaa, to6, 31));
                te6 = _mm512_slli_epi64(te6, 1);
                to6 = _mm512_slli_epi64(to6, 1);

                // finally the two non-doubled terms.
                _mm512_mul_eo64_epi32(a1, a1, &prod1_e, &prod1_o);
                ACCUM_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a0, a0, &prod1_e, &prod1_o);
                ACCUM_EO_PROD(te4, to4, te5, to5);
            }
            else
            {
                // i even, block shape 1.
                // always a continuation of the full-block loop, so use the same 
                // loading pattern.  Only now we don't need as many b-terms.
                a0 = _mm512_load_epi32(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);		// {f, b}
                a1 = _mm512_load_epi32(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);		// {e, a}
                a2 = _mm512_load_epi32(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);		// {d, 9}
                a3 = _mm512_load_epi32(a->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);		// {c, 8}

                b0 = _mm512_load_epi32(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 1) * VECLEN); // {9, 5}
                b1 = _mm512_load_epi32(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 2) * VECLEN);	// {a, 6}
                b2 = _mm512_load_epi32(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 3) * VECLEN);	// {b, 7}
                b3 = _mm512_load_epi32(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 4) * VECLEN);	// {c, 8}
                b4 = _mm512_load_epi32(b->data + (j*BLOCKWORDS + i * BLOCKWORDS + 5) * VECLEN); // {d, 9}

                // save independent sum/carry words for each product-column in the block.
                // uses 11 register inputs, 16 register outputs, and 3 aux vectors.
                //k == 0;
                _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a2, b2, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                //k == 1;
                _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                //k == 2;
                _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

                _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

                //k == 3;
                _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

                _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

                // all terms so far need to be doubled.  Do that all at once with these
                // left shifts.
                te1 = _mm512_or_epi64(te1, _mm512_maskz_srli_epi32(0xaaaa, te0, 31));
                to1 = _mm512_or_epi64(to1, _mm512_maskz_srli_epi32(0xaaaa, to0, 31));
                te0 = _mm512_slli_epi64(te0, 1);
                to0 = _mm512_slli_epi64(to0, 1);

                te3 = _mm512_or_epi64(te3, _mm512_maskz_srli_epi32(0xaaaa, te2, 31));
                to3 = _mm512_or_epi64(to3, _mm512_maskz_srli_epi32(0xaaaa, to2, 31));
                te2 = _mm512_slli_epi64(te2, 1);
                to2 = _mm512_slli_epi64(to2, 1);

                te5 = _mm512_or_epi64(te5, _mm512_maskz_srli_epi32(0xaaaa, te4, 31));
                to5 = _mm512_or_epi64(to5, _mm512_maskz_srli_epi32(0xaaaa, to4, 31));
                te4 = _mm512_slli_epi64(te4, 1);
                to4 = _mm512_slli_epi64(to4, 1);

                te7 = _mm512_or_epi64(te7, _mm512_maskz_srli_epi32(0xaaaa, te6, 31));
                to7 = _mm512_or_epi64(to7, _mm512_maskz_srli_epi32(0xaaaa, to6, 31));
                te6 = _mm512_slli_epi64(te6, 1);
                to6 = _mm512_slli_epi64(to6, 1);

                // finally the two non-doubled terms.
                _mm512_mul_eo64_epi32(a3, a3, &prod1_e, &prod1_o);
                ACCUM_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a2, a2, &prod1_e, &prod1_o);
                ACCUM_EO_PROD(te4, to4, te5, to5);
            }
        }

        // the s*n term.  No more doubling past here.
        for (j = 0; j < NBLOCKS - 1 - i; j++)
        {
            __mmask8 scarry_e1 = 0;
            __mmask8 scarry_o1 = 0;

            a0 = _mm512_load_epi32(s->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
            a1 = _mm512_load_epi32(s->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
            a2 = _mm512_load_epi32(s->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);
            a3 = _mm512_load_epi32(s->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);

            b0 = _mm512_load_epi32(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            //k == 0;
            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a2, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a3, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            //k == 1;
            _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a3, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            //k == 2;
            _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a2, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a3, b5, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            //k == 3;
            _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a2, b5, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a3, b6, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);
        }

        // finish each triangluar shaped column sum (s * n)
        a1 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi32(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi32(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi32(n->data + (NWORDS - 3) * VECLEN);

        // ======
        _mm512_mul_eo64_epi32(a1, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        _mm512_mul_eo64_epi32(a2, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        _mm512_mul_eo64_epi32(a3, b2, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        // ======
        _mm512_mul_eo64_epi32(a2, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te2, to2, te3, to3);

        _mm512_mul_eo64_epi32(a3, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te2, to2, te3, to3);

        // ======
        _mm512_mul_eo64_epi32(a3, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te4, to4, te5, to5);

        j = 0;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te0, to0);
        acc_e1 = _mm512_add_epi64(acc_e1, te1);
        acc_o1 = _mm512_add_epi64(acc_o1, to1);

        // store the low-word final result
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        // and shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 1;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te2, to2);
        acc_e1 = _mm512_add_epi64(acc_e1, te3);
        acc_o1 = _mm512_add_epi64(acc_o1, to3);

        // store the low-word final result
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        // and shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 2;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te4, to4);
        acc_e1 = _mm512_add_epi64(acc_e1, te5);
        acc_o1 = _mm512_add_epi64(acc_o1, to5);

        // store the low-word final result
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        // and shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 3;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te6, to6);
        acc_e1 = _mm512_add_epi64(acc_e1, te7);
        acc_o1 = _mm512_add_epi64(acc_o1, to7);

        // store the low-word final result
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        // and shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;
    }

    a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
    scarry2 = _mm512_cmp_epu32_mask(a0, zero, _MM_CMPINT_EQ);

    // subtract n from tmp
    scarry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        a1 = _mm512_load_epi32(s->data + i * VECLEN);
        b0 = _mm512_load_epi32(n->data + i * VECLEN);
        a0 = _mm512_sbb_epi32(a1, scarry, b0, &scarry);
        _mm512_store_epi32(c->data + i * VECLEN, a0);
    }

    // negate any final borrows if there was also a final carry.
    scarry &= scarry2;

    // if there was a final borrow, we didn't need to do the subtraction after all.
    // replace with original results based on final borrow mask.
    for (i = NWORDS - 1; i >= 0; i--)
    {
        b0 = _mm512_load_epi32(s->data + i * VECLEN);
        _mm512_mask_store_epi32(c->data + i * VECLEN, scarry, b0);
    }

    c->size = NWORDS;
    return;
}

void vecmulmod_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    int i, j, k;
    uint32_t NWORDS = n->WORDS_ALLOC;
    uint32_t NBLOCKS = n->WORDS_ALLOC / BLOCKWORDS;
    __m512i a0, a1, a2, a3;
    __m512i b0, b1, b2, b3, b4, b5, b6;
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;
    __m512i to0, to1, to2, to3, to4, to5, to6, to7;

    __m512i acc_e0;
    __m512i acc_o0;

    __m512i acc_e1;
    __m512i acc_o1;
    // 31

    __m512i prod1_e;
    __m512i prod1_o;

    __m512i hiword = _mm512_set1_epi64(0x000000100000000);
    __m512i zero = _mm512_set1_epi64(0);
    // 37

    __mmask8 scarry_e1 = 0;
    __mmask8 scarry_o1 = 0;
    __mmask8 scarry_e = 0;
    __mmask8 scarry_o = 0;
    __mmask16 scarry2;
    __mmask16 scarry;

    // zero the accumulator
    acc_e0 = acc_o0 = acc_e1 = acc_o1 = zero;


    // first half mul
    for (i = 0; i < NBLOCKS; i++)
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;
        to0 = to1 = to2 = to3 = to4 = to5 = to6 = to7 = zero;

        for (j = i; j > 0; j--)
        {
            a0 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);


            //k == 0;
            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a2, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a3, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            //k == 1;
            _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a3, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            //k == 2;
            _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a2, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a3, b5, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            //k == 3;
            _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a2, b5, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a3, b6, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

        }

        // finish each triangular shaped column sum
        a0 = _mm512_load_epi32(a->data + (i * BLOCKWORDS + 0) * VECLEN);
        a1 = _mm512_load_epi32(a->data + (i * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi32(a->data + (i * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi32(a->data + (i * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi32(b->data + 0 * VECLEN);
        b1 = _mm512_load_epi32(b->data + 1 * VECLEN);
        b2 = _mm512_load_epi32(b->data + 2 * VECLEN);
        b3 = _mm512_load_epi32(b->data + 3 * VECLEN);

        // ======
        _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        // ======
        _mm512_mul_eo64_epi32(a1, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te2, to2, te3, to3);

        _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te2, to2, te3, to3);

        // ======
        _mm512_mul_eo64_epi32(a2, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te4, to4, te5, to5);

        _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te4, to4, te5, to5);

        _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te4, to4, te5, to5);

        // ======
        _mm512_mul_eo64_epi32(a3, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te6, to6, te7, to7);

        _mm512_mul_eo64_epi32(a2, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te6, to6, te7, to7);

        _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te6, to6, te7, to7);

        _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te6, to6, te7, to7);

        // now, do a carry-propagating column summation and store the results.
        j = 0;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te0, to0);
        acc_e1 = _mm512_add_epi64(acc_e1, te1);
        acc_o1 = _mm512_add_epi64(acc_o1, to1);

        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 1;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te2, to2);
        acc_e1 = _mm512_add_epi64(acc_e1, te3);
        acc_o1 = _mm512_add_epi64(acc_o1, to3);

        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 2;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te4, to4);
        acc_e1 = _mm512_add_epi64(acc_e1, te5);
        acc_o1 = _mm512_add_epi64(acc_o1, to5);

        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 3;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te6, to6);
        acc_e1 = _mm512_add_epi64(acc_e1, te7);
        acc_o1 = _mm512_add_epi64(acc_o1, to7);

        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

    }

    // second half mul
    for (i = NBLOCKS; i < 2 * NBLOCKS; i++)
    {

        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;
        to0 = to1 = to2 = to3 = to4 = to5 = to6 = to7 = zero;

        for (j = i - NBLOCKS + 1; j < NBLOCKS; j++)
        {
            a0 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            // accumulate a * b
            //k == 0;
            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a2, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a3, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            //k == 1;
            _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a3, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            //k == 2;
            _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a2, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a3, b5, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            //k == 3;
            _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a2, b5, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a3, b6, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);
        }

        // finish each triangular shaped column sum (a * b)
        a1 = _mm512_load_epi32(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi32(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi32(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi32(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi32(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi32(b->data + (NWORDS - 3) * VECLEN);

        // ======
        _mm512_mul_eo64_epi32(a1, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        _mm512_mul_eo64_epi32(a2, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        _mm512_mul_eo64_epi32(a3, b2, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        // ======
        _mm512_mul_eo64_epi32(a2, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te2, to2, te3, to3);

        _mm512_mul_eo64_epi32(a3, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te2, to2, te3, to3);

        // ======
        _mm512_mul_eo64_epi32(a3, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te4, to4, te5, to5);


        j = 0;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te0, to0);
        acc_e1 = _mm512_add_epi64(acc_e1, te1);
        acc_o1 = _mm512_add_epi64(acc_o1, to1);

        // store the low-word final result
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + ((i) * BLOCKWORDS + j) * VECLEN, a0);

        // and shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 1;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te2, to2);
        acc_e1 = _mm512_add_epi64(acc_e1, te3);
        acc_o1 = _mm512_add_epi64(acc_o1, to3);

        // store the low-word final result
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + ((i) * BLOCKWORDS + j) * VECLEN, a0);

        // and shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 2;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te4, to4);
        acc_e1 = _mm512_add_epi64(acc_e1, te5);
        acc_o1 = _mm512_add_epi64(acc_o1, to5);

        // store the low-word final result
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + ((i) * BLOCKWORDS + j) * VECLEN, a0);

        // and shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 3;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te6, to6);
        acc_e1 = _mm512_add_epi64(acc_e1, te7);
        acc_o1 = _mm512_add_epi64(acc_o1, to7);

        // store the low-word final result
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + ((i) * BLOCKWORDS + j) * VECLEN, a0);

        // and shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;
    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(s, "after hi half:");
#endif

    // reduce by adding hi to lo.  first right shift hi into output.
    vec_bignum_mask_rshift_n(s, c, mdata->nbits, 0xffff);

    int bshift = mdata->nbits % 32;
    int wshift = mdata->nbits / 32;

#ifdef DEBUG_MERSENNE
    print_vechexbignum(c, "hi part:");
    print_vechexbignum(s, "lo part:");
#endif

    // now add the low part into the high.
    scarry = 0;
    for (i = 0; i < wshift; i++)
    {
        a1 = _mm512_load_epi32(c->data + i * VECLEN);
        b0 = _mm512_load_epi32(s->data + i * VECLEN);
        a0 = _mm512_adc_epi32(a1, scarry, b0, &scarry);
        _mm512_store_epi32(c->data + i * VECLEN, a0);
    }

    a1 = _mm512_load_epi32(c->data + i * VECLEN);
    b0 = _mm512_load_epi32(s->data + i * VECLEN);
    b0 = _mm512_and_epi32(_mm512_set1_epi32((1U << (uint32_t)(bshift)) - 1U), b0);
    a0 = _mm512_adc_epi32(a1, scarry, b0, &scarry);
    _mm512_store_epi32(c->data + i * VECLEN, a0);

    for (i++; i < NWORDS; i++)
    {
        _mm512_store_epi32(s->data + i * VECLEN, _mm512_set1_epi32(0));
    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(c, "after add:");
#endif

    // if there was a carry, add it back in.
    a1 = _mm512_load_epi32(c->data + wshift * VECLEN);
    scarry = _mm512_test_epi32_mask(a1, _mm512_set1_epi32((1 << (uint32_t)bshift)));
    i = 0;
    while (scarry > 0)
    {
        a1 = _mm512_load_epi32(c->data + i * VECLEN);
        a0 = _mm512_addcarry_epi32(a1, scarry, &scarry);
        _mm512_store_epi32(c->data + i * VECLEN, a0);
        i++;
    }

    // clear the potential hi-bit
    a1 = _mm512_load_epi32(c->data + wshift * VECLEN);
    _mm512_store_epi32(c->data + wshift * VECLEN,
        _mm512_and_epi32(_mm512_set1_epi32((1 << (uint32_t)(bshift)) - 1), a1));

#ifdef DEBUG_MERSENNE
    print_vechexbignum(c, "after carry add:");
    exit(1);
#endif

    c->size = NWORDS;
    return;
}

void vecsqrmod_mersenne(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    vecmulmod_mersenne(a, a, c, n, s, mdata);
    return;
}

#if 0
void vecsqrmod_mersenne_full(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    int i, j, k;
    vec_bignum_t* b = a;

    __m512i a0, a1, a2, a3;
    __m512i b0, b1, b2, b3, b4, b5, b6;
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;
    __m512i to0, to1, to2, to3, to4, to5, to6, to7;

    __m512i acc_e0;
    __m512i acc_o0;

    __m512i acc_e1;
    __m512i acc_o1;
    // 31

    __m512i prod1_e;
    __m512i prod1_o;

    __m512i hiword = _mm512_set1_epi64(0x000000100000000);
    __m512i hiword2 = _mm512_set1_epi64(0x000000200000000);
    __m512i zero = _mm512_set1_epi64(0);
    // 37

    __mmask8 scarry_e1 = 0;
    __mmask8 scarry_o1 = 0;
    __mmask8 scarry_e = 0;
    __mmask8 scarry_o = 0;
    __mmask16 scarry2;
    __mmask16 scarry;

    // zero the accumulator
    acc_e0 = acc_o0 = acc_e1 = acc_o1 = zero;


    // first half sqr
    for (i = 0; i < NBLOCKS; i++)
    {

        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;
        to0 = to1 = to2 = to3 = to4 = to5 = to6 = to7 = zero;

        if (i & 1)
        {
            __mmask8 scarry_e1 = 0;
            __mmask8 scarry_o1 = 0;

            // i odd
            for (j = 0; j < (i - 1) / 2; j++)
            {
                // for 384-bit inputs NBLOCKS=3 and this loop doesn't run at all.
                // i=0 no, even
                // i=1 no
                // i=2 no, even


                // for 1024-bit inputs NBLOCKS=8 and this loop runs 6 times over all i.
                // i=0 no, even
                // i=1 no
                // i=2 no, even
                // i=3 once
                // i=4 no, even
                // i=5 twice
                // i=6 no, even
                // i=7 thrice
                // with the doubling trick we trade 96 instructions for each j-loop iteration
                // and the 72 instructions for each odd i, after the j-loop,
                // for the 32 instructions for each odd i.  That saves 6*96+4*72-4*32=736 instructions.

                //hips_block_mul_type3(a->data + ((j + 1) * BLOCKWORDS) * VECLEN,
                //    a->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS) * VECLEN, t_e, t_o);

                a0 = _mm512_load_epi32(a->data + ((j + 1) * BLOCKWORDS - 1) * VECLEN);
                a1 = _mm512_load_epi32(a->data + ((j + 1) * BLOCKWORDS - 2) * VECLEN);
                a2 = _mm512_load_epi32(a->data + ((j + 1) * BLOCKWORDS - 3) * VECLEN);
                a3 = _mm512_load_epi32(a->data + ((j + 1) * BLOCKWORDS - 4) * VECLEN);

                b0 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS + 1) * VECLEN);
                b1 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS + 2) * VECLEN);
                b2 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS + 3) * VECLEN);
                b3 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS + 4) * VECLEN);
                b4 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS + 5) * VECLEN);
                b5 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS + 6) * VECLEN);
                b6 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS + 7) * VECLEN);

                // save independent sum/carry words for each product-column in the block.
                // uses 11 register inputs, 16 register outputs, and 3 aux vectors.
                // since all terms in this loop are doubled, we do the doubling
                // after the loop with a left shift.

                //k == 0;
                _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);          // a-1, b+1
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);          // a-2, b+2
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a2, b2, &prod1_e, &prod1_o);          // a-3, b+3
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a3, b3, &prod1_e, &prod1_o);          // a-4, b+4
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                //k == 1;
                _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);          // a-1, b+2
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);          // a-2, b+3
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);          // a-3, b+4
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                _mm512_mul_eo64_epi32(a3, b4, &prod1_e, &prod1_o);          // a-4, b+5
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                //k == 2;
                _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);          // a-1, b+3
                ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

                _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);          // a-2, b+4
                ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

                _mm512_mul_eo64_epi32(a2, b4, &prod1_e, &prod1_o);          // a-3, b+5
                ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

                _mm512_mul_eo64_epi32(a3, b5, &prod1_e, &prod1_o);          // a-4, b+6
                ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

                //k == 3;
                _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);          // a-1, b+4
                ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

                _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);          // a-2, b+5
                ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

                _mm512_mul_eo64_epi32(a2, b5, &prod1_e, &prod1_o);          // a-3, b+6
                ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

                _mm512_mul_eo64_epi32(a3, b6, &prod1_e, &prod1_o);          // a-4, b+7
                ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            }

            // for 384-bit inputs when i == 1, j = 0, a = {3,2,1,0} and b = {2,3,4,5,6,7}
                        // for 512-bit inputs when i == 3, j = 1, a = {7,6,5,4} and b = {6,7,8,9,a,b}
                        // for 512-bit inputs when i == 1, j = 0, a = {3,2,1,0} and b = {2,3,4,5,6,7}

            a0 = _mm512_load_epi32(a->data + (i * BLOCKWORDS - j * BLOCKWORDS - 1) * VECLEN);
            a1 = _mm512_load_epi32(a->data + (i * BLOCKWORDS - j * BLOCKWORDS - 2) * VECLEN);
            a2 = _mm512_load_epi32(a->data + (i * BLOCKWORDS - j * BLOCKWORDS - 3) * VECLEN);
            a3 = _mm512_load_epi32(a->data + (i * BLOCKWORDS - j * BLOCKWORDS - 4) * VECLEN);

            b1 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + 7) * VECLEN);

            // save independent sum/carry words for each product-column in the block.
            // uses 11 register inputs, 16 register outputs, and 3 aux vectors.
            //k == 0;
            _mm512_mul_eo64_epi32(a2, b2, &prod1_e, &prod1_o);
            //te0 = _mm512_addsetc_epi64(te0, prod1_e, &scarry_e);
            //to0 = _mm512_addsetc_epi64(to0, prod1_o, &scarry_o);
            //te1 = _mm512_mask_add_epi64(te1, scarry_e, hiword2, te1);
            //to1 = _mm512_mask_add_epi64(to1, scarry_o, hiword2, to1);
            ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a3, b3, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

            //k == 1;
            _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);
            //te2 = _mm512_addsetc_epi64(te2, prod1_e, &scarry_e);
            //to2 = _mm512_addsetc_epi64(to2, prod1_o, &scarry_o);
            //te3 = _mm512_mask_add_epi64(te3, scarry_e, hiword2, te3);
            //to3 = _mm512_mask_add_epi64(to3, scarry_o, hiword2, to3);

            _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a3, b4, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            //k == 2;
            _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);
            //te4 = _mm512_addsetc_epi64(te4, prod1_e, &scarry_e);
            //to4 = _mm512_addsetc_epi64(to4, prod1_o, &scarry_o);
            //te5 = _mm512_mask_add_epi64(te5, scarry_e, hiword2, te5);
            //to5 = _mm512_mask_add_epi64(to5, scarry_o, hiword2, to5);
            ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a2, b4, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a3, b5, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

            //k == 3;
            _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);
            //te6 = _mm512_addsetc_epi64(te6, prod1_e, &scarry_e);
            //to6 = _mm512_addsetc_epi64(to6, prod1_o, &scarry_o);
            //te7 = _mm512_mask_add_epi64(te7, scarry_e, hiword2, te7);
            //to7 = _mm512_mask_add_epi64(to7, scarry_o, hiword2, to7);

            _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a2, b5, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a3, b6, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            // all terms so far need to be doubled.  Do that all at once with these
            // left shifts.  We need the high-bit of the low word to end up in the
            // 32nd bit position (because the high words are offset by that much
            // for easier combining later on).  _mm512_maskz_srli_epi32 allows us
            // to do the right shift simultaneously with clearing out the low bits.
            te1 = _mm512_or_epi64(te1, _mm512_maskz_srli_epi32(0xaaaa, te0, 31));
            to1 = _mm512_or_epi64(to1, _mm512_maskz_srli_epi32(0xaaaa, to0, 31));
            te0 = _mm512_slli_epi64(te0, 1);
            to0 = _mm512_slli_epi64(to0, 1);

            te3 = _mm512_or_epi64(te3, _mm512_maskz_srli_epi32(0xaaaa, te2, 31));
            to3 = _mm512_or_epi64(to3, _mm512_maskz_srli_epi32(0xaaaa, to2, 31));
            te2 = _mm512_slli_epi64(te2, 1);
            to2 = _mm512_slli_epi64(to2, 1);

            te5 = _mm512_or_epi64(te5, _mm512_maskz_srli_epi32(0xaaaa, te4, 31));
            to5 = _mm512_or_epi64(to5, _mm512_maskz_srli_epi32(0xaaaa, to4, 31));
            te4 = _mm512_slli_epi64(te4, 1);
            to4 = _mm512_slli_epi64(to4, 1);

            te7 = _mm512_or_epi64(te7, _mm512_maskz_srli_epi32(0xaaaa, te6, 31));
            to7 = _mm512_or_epi64(to7, _mm512_maskz_srli_epi32(0xaaaa, to6, 31));
            te6 = _mm512_slli_epi64(te6, 1);
            to6 = _mm512_slli_epi64(to6, 1);

            // finally the two non-doubled terms.
            _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);
            //te0 = _mm512_addsetc_epi64(te0, prod1_e, &scarry_e);
            //to0 = _mm512_addsetc_epi64(to0, prod1_o, &scarry_o);
            //te1 = _mm512_mask_add_epi64(te1, scarry_e, hiword, te1);
            //to1 = _mm512_mask_add_epi64(to1, scarry_o, hiword, to1);

            _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);
            //te4 = _mm512_addsetc_epi64(te4, prod1_e, &scarry_e);
            //to4 = _mm512_addsetc_epi64(to4, prod1_o, &scarry_o);
            //te5 = _mm512_mask_add_epi64(te5, scarry_e, hiword, te5);
            //to5 = _mm512_mask_add_epi64(to5, scarry_o, hiword, to5);

        }
        else
        {
            __mmask8 scarry_e1 = 0;
            __mmask8 scarry_o1 = 0;

            // i even
            for (j = 0; j < i / 2; j++)
            {
                // for 384-bit inputs NBLOCKS=3 and this loop runs once. 
                // i=0: 0 times
                // i=1: no, odd
                // i=2: 1 time

                // for 1024-bit inputs NBLOCKS=8 and this loop runs 6 times over all i.
                // i=0 0 times
                // i=1 no, odd
                // i=2 1 time
                // i=3 no, odd
                // i=4 2 times
                // i=5 no, odd
                // i=6 3 times
                // i=7 no, odd
                // with the doubling trick we trade 96 instructions for each j-loop iteration
                // and the 24 instructions for each even i, after the j-loop,
                // for the 40 instructions once per even i.  That saves 6*96+4*24-4*40=512 instructions.


                //hips_block_mul_type3(a->data + ((j + 1) * BLOCKWORDS) * VECLEN,
                //    a->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS) * VECLEN, t_e, t_o);

                // when i = 2, j = 0, a = {3,2,1,0}, b = {5,6,7,8,9,a,b}
                a0 = _mm512_load_epi32(a->data + ((j + 1) * BLOCKWORDS - 1) * VECLEN);
                a1 = _mm512_load_epi32(a->data + ((j + 1) * BLOCKWORDS - 2) * VECLEN);
                a2 = _mm512_load_epi32(a->data + ((j + 1) * BLOCKWORDS - 3) * VECLEN);
                a3 = _mm512_load_epi32(a->data + ((j + 1) * BLOCKWORDS - 4) * VECLEN);

                b0 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j
                    * BLOCKWORDS + 1) * VECLEN);
                b1 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j
                    * BLOCKWORDS + 2) * VECLEN);
                b2 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j
                    * BLOCKWORDS + 3) * VECLEN);
                b3 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j
                    * BLOCKWORDS + 4) * VECLEN);
                b4 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j
                    * BLOCKWORDS + 5) * VECLEN);
                b5 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j
                    * BLOCKWORDS + 6) * VECLEN);
                b6 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j
                    * BLOCKWORDS + 7) * VECLEN);

                // save independent sum/carry words for each product-column in the block.
                // uses 11 register inputs, 16 register outputs, and 3 aux vectors.
                //k == 0;
                _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);          // a-1, b+1
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);          // a-2, b+2
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a2, b2, &prod1_e, &prod1_o);          // a-3, b+3
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a3, b3, &prod1_e, &prod1_o);          // a-4, b+4
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                //k == 1;
                _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);          // a-1, b+2
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);          // a-2, b+3
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);          // a-3, b+4
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                _mm512_mul_eo64_epi32(a3, b4, &prod1_e, &prod1_o);          // a-4, b+5
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                //k == 2;
                _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);          // a-1, b+3
                ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

                _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);          // a-2, b+4
                ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

                _mm512_mul_eo64_epi32(a2, b4, &prod1_e, &prod1_o);          // a-3, b+5
                ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

                _mm512_mul_eo64_epi32(a3, b5, &prod1_e, &prod1_o);          // a-4, b+6
                ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

                //k == 3;
                _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);          // a-1, b+4
                ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

                _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);          // a-2, b+5
                ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

                _mm512_mul_eo64_epi32(a2, b5, &prod1_e, &prod1_o);          // a-3, b+6
                ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

                _mm512_mul_eo64_epi32(a3, b6, &prod1_e, &prod1_o);          // a-4, b+7
                ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);
            }

            a0 = _mm512_load_epi32(a->data + (i / 2 * BLOCKWORDS + 0) * VECLEN);
            a1 = _mm512_load_epi32(a->data + (i / 2 * BLOCKWORDS + 1) * VECLEN);
            a2 = _mm512_load_epi32(a->data + (i / 2 * BLOCKWORDS + 2) * VECLEN);
            a3 = _mm512_load_epi32(a->data + (i / 2 * BLOCKWORDS + 3) * VECLEN);

            // save independent sum/carry words for each product-column in the block.
            // uses 11 register inputs, 16 register outputs, and 3 aux vectors.

            //k == 1;
            _mm512_mul_eo64_epi32(a0, a1, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            //k == 2;
            _mm512_mul_eo64_epi32(a0, a2, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

            //k == 3;
            _mm512_mul_eo64_epi32(a0, a3, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a1, a2, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            // all terms so far need to be doubled.  Do that all at once with these
            // left shifts.
            te1 = _mm512_or_epi64(te1, _mm512_maskz_srli_epi32(0xaaaa, te0, 31));
            to1 = _mm512_or_epi64(to1, _mm512_maskz_srli_epi32(0xaaaa, to0, 31));
            te0 = _mm512_slli_epi64(te0, 1);
            to0 = _mm512_slli_epi64(to0, 1);

            te3 = _mm512_or_epi64(te3, _mm512_maskz_srli_epi32(0xaaaa, te2, 31));
            to3 = _mm512_or_epi64(to3, _mm512_maskz_srli_epi32(0xaaaa, to2, 31));
            te2 = _mm512_slli_epi64(te2, 1);
            to2 = _mm512_slli_epi64(to2, 1);

            te5 = _mm512_or_epi64(te5, _mm512_maskz_srli_epi32(0xaaaa, te4, 31));
            to5 = _mm512_or_epi64(to5, _mm512_maskz_srli_epi32(0xaaaa, to4, 31));
            te4 = _mm512_slli_epi64(te4, 1);
            to4 = _mm512_slli_epi64(to4, 1);

            te7 = _mm512_or_epi64(te7, _mm512_maskz_srli_epi32(0xaaaa, te6, 31));
            to7 = _mm512_or_epi64(to7, _mm512_maskz_srli_epi32(0xaaaa, to6, 31));
            te6 = _mm512_slli_epi64(te6, 1);
            to6 = _mm512_slli_epi64(to6, 1);

            // finally the two non-doubled terms.
            _mm512_mul_eo64_epi32(a0, a0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a1, a1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

        }

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        for (j = 0; j < i; j++)
        {
            __mmask8 scarry_e1 = 0;
            __mmask8 scarry_o1 = 0;

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

            //k == 0;
            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a2, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a3, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            //k == 1;
            _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a3, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            //k == 2;
            _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a2, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a3, b5, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            //k == 3;
            _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a2, b5, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a3, b6, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

        }		// now, column by column, add in the s*n contribution and reduce to 
        // a single 64+x bit accumulator while storing the intermediate product
        // 's' as we go.
        j = 0;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te0, to0);
        acc_e1 = _mm512_add_epi64(acc_e1, te1);
        acc_o1 = _mm512_add_epi64(acc_o1, to1);

        prod1_e = _mm512_mul_epu32(nhatvec_e, acc_e0);
        prod1_o = _mm512_mul_epu32(nhatvec_o, acc_o0);
        a0 = _mm512_eo64lo_to_epi32(prod1_e, prod1_o);

        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        b0 = _mm512_load_epi32(n->data + 0 * VECLEN);
        _mm512_mul_eo64_epi32(b0, a0, &prod1_e, &prod1_o);

        // add in the final product
        ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 1;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te2, to2);
        acc_e1 = _mm512_add_epi64(acc_e1, te3);
        acc_o1 = _mm512_add_epi64(acc_o1, to3);

        for (k = 0; k < j; k++)
        {
            a0 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + k) * VECLEN);
            b0 = _mm512_load_epi32(n->data + (j - k) * VECLEN);

            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);
        }

        prod1_e = _mm512_mul_epu32(nhatvec_e, acc_e0);
        prod1_o = _mm512_mul_epu32(nhatvec_o, acc_o0);
        a0 = _mm512_eo64lo_to_epi32(prod1_e, prod1_o);

        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        b0 = _mm512_load_epi32(n->data + 0 * VECLEN);
        _mm512_mul_eo64_epi32(b0, a0, &prod1_e, &prod1_o);

        // add in the final product
        ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 2;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te4, to4);
        acc_e1 = _mm512_add_epi64(acc_e1, te5);
        acc_o1 = _mm512_add_epi64(acc_o1, to5);

        for (k = 0; k < j; k++)
        {
            a0 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + k) * VECLEN);
            b0 = _mm512_load_epi32(n->data + (j - k) * VECLEN);

            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);
        }

        prod1_e = _mm512_mul_epu32(nhatvec_e, acc_e0);
        prod1_o = _mm512_mul_epu32(nhatvec_o, acc_o0);
        a0 = _mm512_eo64lo_to_epi32(prod1_e, prod1_o);

        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        b0 = _mm512_load_epi32(n->data + 0 * VECLEN);
        _mm512_mul_eo64_epi32(b0, a0, &prod1_e, &prod1_o);

        // add in the final product
        ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 3;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te6, to6);
        acc_e1 = _mm512_add_epi64(acc_e1, te7);
        acc_o1 = _mm512_add_epi64(acc_o1, to7);

        for (k = 0; k < j; k++)
        {
            a0 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + k) * VECLEN);
            b0 = _mm512_load_epi32(n->data + (j - k) * VECLEN);

            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);
        }

        prod1_e = _mm512_mul_epu32(nhatvec_e, acc_e0);
        prod1_o = _mm512_mul_epu32(nhatvec_o, acc_o0);
        a0 = _mm512_eo64lo_to_epi32(prod1_e, prod1_o);

        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        b0 = _mm512_load_epi32(n->data + 0 * VECLEN);
        _mm512_mul_eo64_epi32(b0, a0, &prod1_e, &prod1_o);

        // add in the final product
        ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

    }

    // second half sqr
    for (i = 0; i < NBLOCKS; i++)
    {


        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;
        to0 = to1 = to2 = to3 = to4 = to5 = to6 = to7 = zero;


        for (j = 0; j < (NBLOCKS - i - 1) / 2; j++)
        {
            __mmask8 scarry_e1 = 0;
            __mmask8 scarry_o1 = 0;

            // Compute a solid block (all matching terms are in the lower
            // half triangle of the expansion).

            //hips_block_mul_type3(a->data + (i * BLOCKWORDS + (j + 2) * BLOCKWORDS) * VECLEN,
            //    a->data + (NWORDS - 2 * BLOCKWORDS - j * BLOCKWORDS) * VECLEN, t_e, t_o);
            a0 = _mm512_load_epi32(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
            a1 = _mm512_load_epi32(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
            a2 = _mm512_load_epi32(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);
            a3 = _mm512_load_epi32(a->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);

            b0 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 7) * VECLEN);

            // save independent sum/carry words for each product-column in the block.
            // uses 11 register inputs, 16 register outputs, and 3 aux vectors.
            //k == 0;
            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);          // a-1, b+1
            ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);          // a-2, b+2
            ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a2, b2, &prod1_e, &prod1_o);          // a-3, b+3
            ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a3, b3, &prod1_e, &prod1_o);          // a-4, b+4
            ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

            //k == 1;
            _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);          // a-1, b+2
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);          // a-2, b+3
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);          // a-3, b+4
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a3, b4, &prod1_e, &prod1_o);          // a-4, b+5
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            //k == 2;
            _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);          // a-1, b+3
            ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);          // a-2, b+4
            ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a2, b4, &prod1_e, &prod1_o);          // a-3, b+5
            ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a3, b5, &prod1_e, &prod1_o);          // a-4, b+6
            ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

            //k == 3;
            _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);          // a-1, b+4
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);          // a-2, b+5
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a2, b5, &prod1_e, &prod1_o);          // a-3, b+6
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a3, b6, &prod1_e, &prod1_o);          // a-4, b+7
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);
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
                a0 = _mm512_load_epi32(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
                a1 = _mm512_load_epi32(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
                a2 = _mm512_load_epi32(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);
                a3 = _mm512_load_epi32(a->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);

                b0 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 1) * VECLEN);
                b1 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 2) * VECLEN);
                b2 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 3) * VECLEN);
                b3 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 4) * VECLEN);
                b4 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 5) * VECLEN);

                // save independent sum/carry words for each product-column in the block.
                // uses 11 register inputs, 16 register outputs, and 3 aux vectors.
                //k == 0;
                _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a2, b2, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                //k == 1;
                _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                //k == 2;
                _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

                _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

                //k == 3;
                _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

                _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

                // all terms so far need to be doubled.  Do that all at once with these
                // left shifts.
                te1 = _mm512_or_epi64(te1, _mm512_maskz_srli_epi32(0xaaaa, te0, 31));
                to1 = _mm512_or_epi64(to1, _mm512_maskz_srli_epi32(0xaaaa, to0, 31));
                te0 = _mm512_slli_epi64(te0, 1);
                to0 = _mm512_slli_epi64(to0, 1);

                te3 = _mm512_or_epi64(te3, _mm512_maskz_srli_epi32(0xaaaa, te2, 31));
                to3 = _mm512_or_epi64(to3, _mm512_maskz_srli_epi32(0xaaaa, to2, 31));
                te2 = _mm512_slli_epi64(te2, 1);
                to2 = _mm512_slli_epi64(to2, 1);

                te5 = _mm512_or_epi64(te5, _mm512_maskz_srli_epi32(0xaaaa, te4, 31));
                to5 = _mm512_or_epi64(to5, _mm512_maskz_srli_epi32(0xaaaa, to4, 31));
                te4 = _mm512_slli_epi64(te4, 1);
                to4 = _mm512_slli_epi64(to4, 1);

                te7 = _mm512_or_epi64(te7, _mm512_maskz_srli_epi32(0xaaaa, te6, 31));
                to7 = _mm512_or_epi64(to7, _mm512_maskz_srli_epi32(0xaaaa, to6, 31));
                te6 = _mm512_slli_epi64(te6, 1);
                to6 = _mm512_slli_epi64(to6, 1);

                // finally the two non-doubled terms.
                _mm512_mul_eo64_epi32(a3, a3, &prod1_e, &prod1_o);
                ACCUM_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a2, a2, &prod1_e, &prod1_o);
                ACCUM_EO_PROD(te4, to4, te5, to5);
            }
            else
            {
                // i even, block shape 1.
                // always a continuation of the full-block loop, so use the same 
                // loading pattern.  Only now we don't need as many b-terms.
                a0 = _mm512_load_epi32(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
                a1 = _mm512_load_epi32(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
                a2 = _mm512_load_epi32(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);

                //k == 0;
                _mm512_mul_eo64_epi32(a0, a2, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                //k == 1;
                _mm512_mul_eo64_epi32(a0, a1, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                te1 = _mm512_or_epi64(te1, _mm512_maskz_srli_epi32(0xaaaa, te0, 31));
                to1 = _mm512_or_epi64(to1, _mm512_maskz_srli_epi32(0xaaaa, to0, 31));
                te0 = _mm512_slli_epi64(te0, 1);
                to0 = _mm512_slli_epi64(to0, 1);

                te3 = _mm512_or_epi64(te3, _mm512_maskz_srli_epi32(0xaaaa, te2, 31));
                to3 = _mm512_or_epi64(to3, _mm512_maskz_srli_epi32(0xaaaa, to2, 31));
                te2 = _mm512_slli_epi64(te2, 1);
                to2 = _mm512_slli_epi64(to2, 1);

                // technically only have to do these two if j > 0 (so that
                // they are non-zero from full-block loop iterations).
                te5 = _mm512_or_epi64(te5, _mm512_maskz_srli_epi32(0xaaaa, te4, 31));
                to5 = _mm512_or_epi64(to5, _mm512_maskz_srli_epi32(0xaaaa, to4, 31));
                te4 = _mm512_slli_epi64(te4, 1);
                to4 = _mm512_slli_epi64(to4, 1);

                te7 = _mm512_or_epi64(te7, _mm512_maskz_srli_epi32(0xaaaa, te6, 31));
                to7 = _mm512_or_epi64(to7, _mm512_maskz_srli_epi32(0xaaaa, to6, 31));
                te6 = _mm512_slli_epi64(te6, 1);
                to6 = _mm512_slli_epi64(to6, 1);

                // finally the two non-doubled terms.
                _mm512_mul_eo64_epi32(a1, a1, &prod1_e, &prod1_o);
                ACCUM_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a0, a0, &prod1_e, &prod1_o);
                ACCUM_EO_PROD(te4, to4, te5, to5);
            }
        }

        else
        {
            // NBLOCKS is even
            // i odd, block shape 1.
            if (i & 1)
            {
                // always a continuation of the full-block loop, so use the same 
                // loading pattern.  Only now we don't need as many b-terms.
                a0 = _mm512_load_epi32(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);  // {f, b}
                a1 = _mm512_load_epi32(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);  // {e, a}
                a2 = _mm512_load_epi32(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);  // {d, 9}

                //k == 0;
                _mm512_mul_eo64_epi32(a0, a2, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                //k == 1;
                _mm512_mul_eo64_epi32(a0, a1, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                te1 = _mm512_or_epi64(te1, _mm512_maskz_srli_epi32(0xaaaa, te0, 31));
                to1 = _mm512_or_epi64(to1, _mm512_maskz_srli_epi32(0xaaaa, to0, 31));
                te0 = _mm512_slli_epi64(te0, 1);
                to0 = _mm512_slli_epi64(to0, 1);

                te3 = _mm512_or_epi64(te3, _mm512_maskz_srli_epi32(0xaaaa, te2, 31));
                to3 = _mm512_or_epi64(to3, _mm512_maskz_srli_epi32(0xaaaa, to2, 31));
                te2 = _mm512_slli_epi64(te2, 1);
                to2 = _mm512_slli_epi64(to2, 1);

                // technically only have to do these two if j > 0 (so that
                // they are non-zero from full-block loop iterations).
                te5 = _mm512_or_epi64(te5, _mm512_maskz_srli_epi32(0xaaaa, te4, 31));
                to5 = _mm512_or_epi64(to5, _mm512_maskz_srli_epi32(0xaaaa, to4, 31));
                te4 = _mm512_slli_epi64(te4, 1);
                to4 = _mm512_slli_epi64(to4, 1);

                te7 = _mm512_or_epi64(te7, _mm512_maskz_srli_epi32(0xaaaa, te6, 31));
                to7 = _mm512_or_epi64(to7, _mm512_maskz_srli_epi32(0xaaaa, to6, 31));
                te6 = _mm512_slli_epi64(te6, 1);
                to6 = _mm512_slli_epi64(to6, 1);

                // finally the two non-doubled terms.
                _mm512_mul_eo64_epi32(a1, a1, &prod1_e, &prod1_o);
                ACCUM_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a0, a0, &prod1_e, &prod1_o);
                ACCUM_EO_PROD(te4, to4, te5, to5);
            }
            else
            {
                // i even, block shape 1.
                // always a continuation of the full-block loop, so use the same 
                // loading pattern.  Only now we don't need as many b-terms.
                a0 = _mm512_load_epi32(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);		// {f, b}
                a1 = _mm512_load_epi32(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);		// {e, a}
                a2 = _mm512_load_epi32(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);		// {d, 9}
                a3 = _mm512_load_epi32(a->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);		// {c, 8}

                b0 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 1) * VECLEN); // {9, 5}
                b1 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 2) * VECLEN);	// {a, 6}
                b2 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 3) * VECLEN);	// {b, 7}
                b3 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 4) * VECLEN);	// {c, 8}
                b4 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 5) * VECLEN); // {d, 9}

                // save independent sum/carry words for each product-column in the block.
                // uses 11 register inputs, 16 register outputs, and 3 aux vectors.
                //k == 0;
                _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a2, b2, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

                //k == 1;
                _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

                //k == 2;
                _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

                _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

                //k == 3;
                _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

                _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);
                ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

                // all terms so far need to be doubled.  Do that all at once with these
                // left shifts.
                te1 = _mm512_or_epi64(te1, _mm512_maskz_srli_epi32(0xaaaa, te0, 31));
                to1 = _mm512_or_epi64(to1, _mm512_maskz_srli_epi32(0xaaaa, to0, 31));
                te0 = _mm512_slli_epi64(te0, 1);
                to0 = _mm512_slli_epi64(to0, 1);

                te3 = _mm512_or_epi64(te3, _mm512_maskz_srli_epi32(0xaaaa, te2, 31));
                to3 = _mm512_or_epi64(to3, _mm512_maskz_srli_epi32(0xaaaa, to2, 31));
                te2 = _mm512_slli_epi64(te2, 1);
                to2 = _mm512_slli_epi64(to2, 1);

                te5 = _mm512_or_epi64(te5, _mm512_maskz_srli_epi32(0xaaaa, te4, 31));
                to5 = _mm512_or_epi64(to5, _mm512_maskz_srli_epi32(0xaaaa, to4, 31));
                te4 = _mm512_slli_epi64(te4, 1);
                to4 = _mm512_slli_epi64(to4, 1);

                te7 = _mm512_or_epi64(te7, _mm512_maskz_srli_epi32(0xaaaa, te6, 31));
                to7 = _mm512_or_epi64(to7, _mm512_maskz_srli_epi32(0xaaaa, to6, 31));
                te6 = _mm512_slli_epi64(te6, 1);
                to6 = _mm512_slli_epi64(to6, 1);

                // finally the two non-doubled terms.
                _mm512_mul_eo64_epi32(a3, a3, &prod1_e, &prod1_o);
                ACCUM_EO_PROD(te0, to0, te1, to1);

                _mm512_mul_eo64_epi32(a2, a2, &prod1_e, &prod1_o);
                ACCUM_EO_PROD(te4, to4, te5, to5);
            }
        }

        // the s*n term.  No more doubling past here.
        for (j = 0; j < NBLOCKS - 1 - i; j++)
        {
            __mmask8 scarry_e1 = 0;
            __mmask8 scarry_o1 = 0;

            a0 = _mm512_load_epi32(s->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
            a1 = _mm512_load_epi32(s->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
            a2 = _mm512_load_epi32(s->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);
            a3 = _mm512_load_epi32(s->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);

            b0 = _mm512_load_epi32(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi32(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi32(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);

            //k == 0;
            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a2, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a3, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            //k == 1;
            _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a3, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te2, to2, te3, to3);

            //k == 2;
            _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a2, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a3, b5, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);

            //k == 3;
            _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a2, b5, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a3, b6, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te6, to6, te7, to7);
        }

        // finish each triangluar shaped column sum (s * n)
        a1 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi32(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi32(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi32(n->data + (NWORDS - 3) * VECLEN);

        // ======
        _mm512_mul_eo64_epi32(a1, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        _mm512_mul_eo64_epi32(a2, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        _mm512_mul_eo64_epi32(a3, b2, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        // ======
        _mm512_mul_eo64_epi32(a2, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te2, to2, te3, to3);

        _mm512_mul_eo64_epi32(a3, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te2, to2, te3, to3);

        // ======
        _mm512_mul_eo64_epi32(a3, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te4, to4, te5, to5);

        j = 0;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te0, to0);
        acc_e1 = _mm512_add_epi64(acc_e1, te1);
        acc_o1 = _mm512_add_epi64(acc_o1, to1);

        // store the low-word final result
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        // and shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 1;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te2, to2);
        acc_e1 = _mm512_add_epi64(acc_e1, te3);
        acc_o1 = _mm512_add_epi64(acc_o1, to3);

        // store the low-word final result
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        // and shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 2;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te4, to4);
        acc_e1 = _mm512_add_epi64(acc_e1, te5);
        acc_o1 = _mm512_add_epi64(acc_o1, to5);

        // store the low-word final result
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        // and shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 3;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te6, to6);
        acc_e1 = _mm512_add_epi64(acc_e1, te7);
        acc_o1 = _mm512_add_epi64(acc_o1, to7);

        // store the low-word final result
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        // and shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;
    }

    a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
    scarry2 = _mm512_cmp_epu32_mask(a0, zero, _MM_CMPINT_EQ);

    // subtract n from tmp
    scarry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        a1 = _mm512_load_epi32(s->data + i * VECLEN);
        b0 = _mm512_load_epi32(n->data + i * VECLEN);
        a0 = _mm512_sbb_epi32(a1, scarry, b0, &scarry);
        _mm512_store_epi32(c->data + i * VECLEN, a0);
    }

    // negate any final borrows if there was also a final carry.
    scarry &= scarry2;

    // if there was a final borrow, we didn't need to do the subtraction after all.
    // replace with original results based on final borrow mask.
    for (i = NWORDS - 1; i >= 0; i--)
    {
        b0 = _mm512_load_epi32(s->data + i * VECLEN);
        _mm512_mask_store_epi32(c->data + i * VECLEN, scarry, b0);
    }

    c->size = NWORDS;
    return;
}
#endif

void vec_simul_addsub(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *sum, vec_bignum_t *diff, vec_monty_t* mdata)
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
    __mmask16 carry = 0;
    __mmask16 borrow = 0;
    __mmask16 cmask = 0;
    __mmask16 cmask2 = 0;
    __mmask16 bmask = 0;
    __mmask16 bmask2 = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    __m512i nvec;


    for (i = 0; i < NWORDS; i++)
    {
        // add
        avec = _mm512_load_epi32(a->data + i * VECLEN);
        bvec = _mm512_load_epi32(b->data + i * VECLEN);
        cvec = _mm512_adc_epi32(avec, carry, bvec, &carry);
        _mm512_store_epi32(sum->data + i * VECLEN, cvec);

        // sub
        cvec = _mm512_sbb_epi32(avec, borrow, bvec, &borrow);
        _mm512_store_epi32(diff->data + i * VECLEN, cvec);
    }

    cmask = carry;	    // result too big, need to subtract n
    bmask = borrow;     // result too small, need to add n
    cmask2 = cmask;	    // keep looking mask for add

    // compare.  the initial mask is equal to the addition carry
    // because if the most significant word has a carry, then the
    // result is bigger than n.
    for (i = NWORDS - 1; i >= 0; i--)
    {
        cvec = _mm512_load_epi32(sum->data + i * VECLEN);
        nvec = _mm512_load_epi32(mdata->n->data + i * VECLEN);
        // compare those that have not already been decided using the mask
        cmask |= _mm512_mask_cmp_epu32_mask(~cmask2, cvec, nvec, _MM_CMPINT_GT);
        cmask2 |= _mm512_mask_cmp_epu32_mask(~cmask2, cvec, nvec, _MM_CMPINT_LT);

        // decided all of them, stop comparing.
        if (cmask2 == 0xffff)
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
        cvec = _mm512_load_epi32(sum->data + i * VECLEN);
        nvec = _mm512_load_epi32(mdata->n->data + i * VECLEN);
        bvec = _mm512_mask_sbb_epi32(cvec, cmask, borrow, nvec, &borrow);
        _mm512_store_epi32(sum->data + i * VECLEN, bvec);

        // conditional add
        cvec = _mm512_load_epi32(diff->data + i * VECLEN);
        bvec = _mm512_mask_adc_epi32(cvec, bmask, carry, nvec, &carry);
        _mm512_store_epi32(diff->data + i * VECLEN, bvec);
    }

    return;
}

void vecaddmod(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, vec_monty_t* mdata)
{
    // assumptions:
    // a, b, c are of length VECLEN * NWORDS
    // a, b, c, and n are aligned
    // a and b are both positive
    // n is the montgomery base
    int i;
    uint32_t NWORDS = a->WORDS_ALLOC;
    __mmask16 carry = 0;
    __mmask16 mask = 0;
    __mmask16 mask2 = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    __m512i nvec;

    // add
    for (i = 0; i < NWORDS; i++)
    {
        avec = _mm512_load_epi32(a->data + i * VECLEN);
        bvec = _mm512_load_epi32(b->data + i * VECLEN);
        cvec = _mm512_adc_epi32(avec, carry, bvec, &carry);
        _mm512_store_epi32(c->data + i * VECLEN, cvec);
    }

    mask = carry;	// sub mask
    mask2 = mask;	// keep looking mask

    // compare.  the initial mask is equal to the addition carry
    // because if the most significant word has a carry, then the
    // result is bigger than n.
    for (i = NWORDS - 1; i >= 0; i--)
    {
        cvec = _mm512_load_epi32(c->data + i * VECLEN);
        nvec = _mm512_load_epi32(mdata->n->data + i * VECLEN);
        // compare those that have not already been decided using the mask
        mask |= _mm512_mask_cmp_epu32_mask(~mask2, cvec, nvec, _MM_CMPINT_GT);
        mask2 |= _mm512_mask_cmp_epu32_mask(~mask2, cvec, nvec, _MM_CMPINT_LT);

        // decided all of them, stop comparing.
        if (mask2 == 0xffff)
        {
            break;
        }
    }

    // check for equal as well by flipping mask bits that have still
    // not been decided (i.e., are equal)
    mask |= (~mask2);

    // subtract n from c when c is not less than n, as indicated by a 1 bit in mask
    carry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        cvec = _mm512_load_epi32(c->data + i * VECLEN);
        nvec = _mm512_load_epi32(mdata->n->data + i * VECLEN);
        bvec = _mm512_mask_sbb_epi32(cvec, mask, carry, nvec, &carry);
        _mm512_store_epi32(c->data + i * VECLEN, bvec);
    }

    return;
}

void vecsubmod(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, vec_monty_t* mdata)
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
    __mmask16 carry = 0;
    __mmask16 mask = 0;
    __mmask16 mask2 = 0;
    __m512i nvec;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;

    // subtract
    carry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        avec = _mm512_load_epi32(a->data + i * VECLEN);
        bvec = _mm512_load_epi32(b->data + i * VECLEN);
        cvec = _mm512_sbb_epi32(avec, carry, bvec, &carry);
        _mm512_store_epi32(c->data + i * VECLEN, cvec);
    }

    // if we had a final carry, then b was bigger than a so we need to re-add n.
    mask = carry;
    carry = 0;
    for (i = 0; (i < NWORDS) && (mask > 0); i++)
    {
        avec = _mm512_load_epi32(c->data + i * VECLEN);
        nvec = _mm512_load_epi32(mdata->n->data + i * VECLEN);
        cvec = _mm512_mask_adc_epi32(avec, mask, carry, nvec, &carry);
        _mm512_store_epi32(c->data + i * VECLEN, cvec);
    }

    return;
}

void vecaddmod_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata)
{
    // assumptions:
    // a, b, c are of length VECLEN * NWORDS
    // a, b, c, and n are aligned
    // a and b are both positive
    // n is the montgomery base
    int i;
    uint32_t NWORDS = a->WORDS_ALLOC;
    __mmask16 carry = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    int bshift = mdata->nbits % 32;
    int wshift = mdata->nbits / 32;

    // add
    for (i = 0; i < NWORDS; i++)
    {
        avec = _mm512_load_epi32(a->data + i * VECLEN);
        bvec = _mm512_load_epi32(b->data + i * VECLEN);
        cvec = _mm512_adc_epi32(avec, carry, bvec, &carry);
        _mm512_store_epi32(c->data + i * VECLEN, cvec);
    }

    // check for a carry.
    avec = _mm512_load_epi32(c->data + wshift * VECLEN);
    carry = _mm512_test_epi32_mask(avec, _mm512_set1_epi32((1 << bshift)));

    // the modulo is just the low part plus 1 (the carry, if present).
    cvec = _mm512_load_epi32(c->data + 0 * VECLEN);
    bvec = _mm512_addcarry_epi32(cvec, carry, &carry);
    _mm512_store_epi32(c->data + 0 * VECLEN, bvec);

    for (i = 1; (i < NWORDS) && (carry > 0); i++)
    {
        cvec = _mm512_load_epi32(c->data + i * VECLEN);
        bvec = _mm512_addcarry_epi32(cvec, carry, &carry);
        _mm512_store_epi32(c->data + i * VECLEN, bvec);
    }

    // clear the potential hi-bit
    avec = _mm512_load_epi32(c->data + wshift * VECLEN);
    _mm512_store_epi32(c->data + wshift * VECLEN,
        _mm512_and_epi32(_mm512_set1_epi64((1 << (bshift)) - 1), avec));

    return;
}

void vecsubmod_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata)
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
    __mmask16 carry = 0;
    __mmask16 mask = 0;
    __m512i nvec;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    int bshift = mdata->nbits % 32;
    int wshift = mdata->nbits / 32;

    // subtract
    carry = 0;
    for (i = 0; i <= wshift; i++)
    {
        avec = _mm512_load_epi32(a->data + i * VECLEN);
        bvec = _mm512_load_epi32(b->data + i * VECLEN);
        cvec = _mm512_sbb_epi32(avec, carry, bvec, &carry);
        _mm512_store_epi32(c->data + i * VECLEN, cvec);
    }

    // if we had a final carry, then b was bigger than a so we need to re-add n.
    mask = carry;
    carry = 0;
    nvec = _mm512_set1_epi32(0xffffffff);
    for (i = 0; i <= wshift; i++)
    {
        cvec = _mm512_load_epi32(c->data + i * VECLEN);
        bvec = _mm512_mask_adc_epi32(cvec, mask, carry, nvec, &carry);
        _mm512_store_epi32(c->data + i * VECLEN, bvec);
    }

    // clear the potential hi-bit
    avec = _mm512_load_epi32(c->data + wshift * VECLEN);
    _mm512_store_epi32(c->data + wshift * VECLEN,
        _mm512_and_epi32(_mm512_set1_epi64((1 << (bshift)) - 1), avec));

    return;
}

void vec_simul_addsub_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* sum, vec_bignum_t* diff,
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
    __mmask16 carry = 0;
    __mmask16 borrow = 0;
    __mmask16 bmask = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    __m512i nvec;
    int bshift = mdata->nbits % 32;
    int wshift = mdata->nbits / 32;

    for (i = 0; i <= wshift; i++)
    {
        // add
        avec = _mm512_load_epi32(a->data + i * VECLEN);
        bvec = _mm512_load_epi32(b->data + i * VECLEN);
        cvec = _mm512_adc_epi32(avec, carry, bvec, &carry);
        _mm512_store_epi32(sum->data + i * VECLEN, cvec);

        // sub
        cvec = _mm512_sbb_epi32(avec, borrow, bvec, &borrow);
        _mm512_store_epi32(diff->data + i * VECLEN, cvec);
    }

    bmask = borrow;     // result too small, need to add n

    // check for a carry.
    avec = _mm512_load_epi32(sum->data + wshift * VECLEN);
    carry = _mm512_test_epi32_mask(avec, _mm512_set1_epi32((1 << bshift)));

    // the modulo is just the low part plus 1 (the carry, if present).
    cvec = _mm512_load_epi32(sum->data + 0 * VECLEN);
    bvec = _mm512_addcarry_epi32(cvec, carry, &carry);
    _mm512_store_epi32(sum->data + 0 * VECLEN, bvec);

    for (i = 1; (i < NWORDS) && (carry > 0); i++)
    {
        cvec = _mm512_load_epi32(sum->data + i * VECLEN);
        bvec = _mm512_addcarry_epi32(cvec, carry, &carry);
        _mm512_store_epi32(sum->data + i * VECLEN, bvec);
    }

    // clear the potential hi-bit
    avec = _mm512_load_epi32(sum->data + wshift * VECLEN);
    _mm512_store_epi32(sum->data + wshift * VECLEN,
        _mm512_and_epi32(_mm512_set1_epi32((1 << (bshift)) - 1), avec));

    carry = 0;
    nvec = _mm512_set1_epi32(0xffffffffULL);
    for (i = 0; i <= wshift; i++)
    {
        // conditional add
        cvec = _mm512_load_epi32(diff->data + i * VECLEN);
        bvec = _mm512_mask_adc_epi32(cvec, bmask, carry, nvec, &carry);
        _mm512_store_epi32(diff->data + i * VECLEN, bvec);
    }

    // clear the potential hi-bit
    avec = _mm512_load_epi32(diff->data + wshift * VECLEN);
    _mm512_store_epi32(diff->data + wshift * VECLEN,
        _mm512_and_epi32(_mm512_set1_epi32((1 << (bshift)) - 1), avec));

    return;
}

uint32_t vec_gte(vec_bignum_t * u, vec_bignum_t * v)
{
    // decide if each of the bignums in vec 'u' is >=
    // the corresponding bignum in vec 'v'.
    // return a mask of results.
    int i;
    uint32_t NWORDS = u->WORDS_ALLOC;
    __mmask16 mdecided = 0;
    __mmask16 mgte = 0;

    for (i = NWORDS - 1; i >= 0; --i)
    {
        __m512i a = _mm512_load_epi32(u->data + i * VECLEN);
        __m512i b = _mm512_load_epi32(v->data + i * VECLEN);

        mgte |= _mm512_mask_cmp_epu32_mask(~mdecided, a, b, _MM_CMPINT_GT);
        mdecided = mdecided | _mm512_mask_cmp_epu32_mask(~mdecided, a, b, _MM_CMPINT_LT);

        if (mdecided == 0xffff)
            break;
    }

    //equal if still undecided
    mgte |= ~mdecided;

    return (uint32_t)mgte;
}

uint32_t vec_eq(base_t * u, base_t * v, int sz)
{
    // decide if each of the bignums in vec 'u' is >=
    // the corresponding bignum in vec 'v'.
    // return a mask of results.
    int i;
    __mmask16 meq = 0xffff;

    for (i = sz - 1; i >= 0; --i)
    {
        __m512i a = _mm512_load_epi32(u + i * VECLEN);
        __m512i b = _mm512_load_epi32(v + i * VECLEN);

        meq = _mm512_mask_cmp_epu32_mask(meq, a, b, _MM_CMPINT_EQ);

        if (meq == 0)
            break;
    }

    return (uint32_t)meq;
}

uint32_t vec_bignum_mask_lshift_1(vec_bignum_t * u, uint32_t wmask)
{
    // return the left shift of bignum u by 1
    int i;
    uint32_t NWORDS = u->WORDS_ALLOC;
    __m512i nextcarry;
    __m512i carry = _mm512_setzero_epi32();
    __m512i highmask = _mm512_set1_epi32(0x80000000);
    __m512i word;

    for (i = 0; i < NWORDS; i++)
    {
        word = _mm512_load_epi32(u->data + i * VECLEN);
        // _mm512_and_epi32(word, highmask)  // not necessary to mask as the shift zero extends.
        nextcarry = _mm512_srli_epi32(word, 31);
        _mm512_mask_store_epi32(u->data + i * VECLEN, (__mmask16)wmask,
            _mm512_or_epi32(_mm512_slli_epi32(word, 1), carry));
        carry = nextcarry;
    }

    _mm512_mask_store_epi32(u->data + i * VECLEN, (__mmask16)wmask,
        _mm512_or_epi32(_mm512_slli_epi32(word, 1), carry));

    // return an overflow mask
    return wmask & _mm512_cmp_epi32_mask(carry, _mm512_setzero_epi32(), _MM_CMPINT_GT);
}

void vec_bignum_mask_rshift_n(vec_bignum_t* u, vec_bignum_t* v, int n, uint32_t wmask)
{
    // return the right shift of bignum u by n bits
    int i;
    uint32_t NWORDS = u->WORDS_ALLOC;
    __m512i nextcarry;
    __m512i carry = _mm512_set1_epi32(0);
    __m512i lowmask;
    __m512i word;
    int wshift = n / 32;
    int bshift = n % 32;

    lowmask = _mm512_set1_epi32((1ULL << (uint32_t)bshift) - 1ULL);

    for (i = 2 * NWORDS - 1; (i - wshift) >= 0; i--)
    {
        word = _mm512_load_epi32(u->data + i * VECLEN);
        nextcarry = _mm512_slli_epi32(_mm512_and_epi32(word, lowmask), (DIGITBITS - bshift));
        _mm512_mask_store_epi64(v->data + (i - wshift) * VECLEN, (__mmask16)wmask,
            _mm512_or_epi32(_mm512_srli_epi32(word, bshift), carry));
        carry = nextcarry;
    }

    return;
}

void vec_bignum_mask_rshift_1(vec_bignum_t * u, uint32_t wmask)
{
    // return the right shift of bignum u by 1
    int i;
    uint32_t NWORDS = u->WORDS_ALLOC;
    __m512i nextcarry;
    __m512i carry = _mm512_setzero_epi32();
    __m512i lowmask = _mm512_set1_epi32(0x00000001);
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
        word = _mm512_load_epi32(u->data + i * VECLEN);
        nextcarry = _mm512_slli_epi32(_mm512_and_epi32(word, lowmask), 31);
        _mm512_mask_store_epi32(u->data + i * VECLEN, (__mmask16)wmask,
            _mm512_or_epi32(_mm512_srli_epi32(word, 1), carry));
        carry = nextcarry;
    }

    return;
}

void vec_bignum_mask_sub(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, uint32_t wmask)
{
    // assumptions:
    // a, b, c are of length VECLEN * NWORDS
    // s1 is of length VECLEN
    // a, b, c, and s1 are aligned
    // a and b are both positive
    // a >= b
    int i;
    uint32_t NWORDS = a->WORDS_ALLOC;
    __mmask16 carry = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;

    if (wmask == 0)
        return;

    // subtract the selected elements ('1' in the mask)
    carry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        avec = _mm512_load_epi32(a->data + i * VECLEN);
        bvec = _mm512_load_epi32(b->data + i * VECLEN);
        cvec = _mm512_sbb_epi32(avec, carry, bvec, &carry);
        _mm512_mask_store_epi32(c->data + i * VECLEN, (__mmask16)wmask, cvec);
    }

    if (carry)
    {
        // subtract any final borrows that exceed the size of b.
        _mm512_mask_store_epi32(c->data + i * VECLEN, (__mmask16)wmask & carry, _mm512_setzero_epi32());
    }

    return;
}

#else

void vecsqrmod(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    return;
}

void vecmulmod(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    return;
}

void vecaddmod(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata)
{
    return;
}

void vecsubmod(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata)
{
    return;
}

void vec_simul_addsub(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* sum, vec_bignum_t* diff, vec_monty_t* mdata)
{
    return;
}

#endif

