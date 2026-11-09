/*
Copyright (c) 2026, Ben Buhrow
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
*/

/*
 * AVX2 versions of AVX-ECM's radix-2^52 add / sub routines (builds without
 * AVX-512).  Same 8-lane vec_bignum_t layout, processed as two 4-lane
 * halves, and the same results as the AVX-512 routines in vecarith52.c
 * (USE_AMM semantics):
 *   vecaddmod52_avx2:  c = a + b, minus N in the lanes whose sum carried
 *                      out of the top word
 *   vecsubmod52_avx2:  c = a - b, plus N in the lanes that borrowed, plus N
 *                      again where that add did not carry out
 *   vec_simul_addsub52_avx2: sum = a + b (minus N on carry) and
 *                      diff = a - b (plus N once on borrow)
 * Digits are 52-bit in 64-bit lanes, so carries are plain shifts and
 * borrows are sign tests.  The output may alias an input.
 */

#include "avx_ecm.h"

#if defined(USE_AVX2) && !defined(USE_AVX512F)
#include <immintrin.h>

#define LDH(p, i, h) _mm256_load_si256((const __m256i *)(const void *)((p) + (size_t)(i) * 8 + (h)))
#define STH(p, i, h, v) _mm256_store_si256((__m256i *)(void *)((p) + (size_t)(i) * 8 + (h)), (v))

static inline int any_set(__m256i m)
{
    return !_mm256_testz_si256(m, m);
}

// c (lanes in m) += N, carry chain; returns the lanes that carried out
static inline __m256i add_n_masked(uint64_t *c, const uint64_t *n, __m256i m,
    int NWORDS, int h)
{
    const __m256i lomask = _mm256_set1_epi64x(0xfffffffffffffLL);
    __m256i carry = _mm256_setzero_si256();
    int i;
    for (i = 0; i < NWORDS; i++)
    {
        __m256i t = _mm256_add_epi64(_mm256_add_epi64(LDH(c, i, h),
            _mm256_and_si256(LDH(n, i, h), m)), carry);
        carry = _mm256_srli_epi64(t, 52);
        STH(c, i, h, _mm256_and_si256(t, lomask));
    }
    return _mm256_cmpgt_epi64(carry, _mm256_setzero_si256());
}

// c (lanes in m) -= N, borrow chain
static inline void sub_n_masked(uint64_t *c, const uint64_t *n, __m256i m,
    int NWORDS, int h)
{
    const __m256i lomask = _mm256_set1_epi64x(0xfffffffffffffLL);
    const __m256i zero = _mm256_setzero_si256();
    __m256i borrow = zero;      // 0 or -1
    int i;
    for (i = 0; i < NWORDS; i++)
    {
        __m256i t = _mm256_add_epi64(_mm256_sub_epi64(LDH(c, i, h),
            _mm256_and_si256(LDH(n, i, h), m)), borrow);
        borrow = _mm256_cmpgt_epi64(zero, t);
        STH(c, i, h, _mm256_and_si256(t, lomask));
    }
}

void vecaddmod52_avx2(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, vec_monty_t *mdata)
{
    const __m256i lomask = _mm256_set1_epi64x(0xfffffffffffffLL);
    const int NWORDS = (int)mdata->NWORDS;
    int h, i;
    for (h = 0; h < 8; h += 4)
    {
        __m256i carry = _mm256_setzero_si256(), m;
        for (i = 0; i < NWORDS; i++)
        {
            __m256i t = _mm256_add_epi64(_mm256_add_epi64(LDH(a->data, i, h),
                LDH(b->data, i, h)), carry);
            carry = _mm256_srli_epi64(t, 52);
            STH(c->data, i, h, _mm256_and_si256(t, lomask));
        }
        m = _mm256_cmpgt_epi64(carry, _mm256_setzero_si256());
        if (any_set(m))
            sub_n_masked(c->data, mdata->n->data, m, NWORDS, h);
    }
}

void vecsubmod52_avx2(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, vec_monty_t *mdata)
{
    const __m256i lomask = _mm256_set1_epi64x(0xfffffffffffffLL);
    const __m256i zero = _mm256_setzero_si256();
    const int NWORDS = (int)mdata->NWORDS;
    int h, i;
    for (h = 0; h < 8; h += 4)
    {
        __m256i borrow = zero, m;
        for (i = 0; i < NWORDS; i++)
        {
            __m256i t = _mm256_add_epi64(_mm256_sub_epi64(LDH(a->data, i, h),
                LDH(b->data, i, h)), borrow);
            borrow = _mm256_cmpgt_epi64(zero, t);
            STH(c->data, i, h, _mm256_and_si256(t, lomask));
        }
        m = borrow;                     // lanes that went negative
        if (any_set(m))
        {
            __m256i co = add_n_masked(c->data, mdata->n->data, m, NWORDS, h);
            m = _mm256_andnot_si256(co, m);     // still negative: add N again
            if (any_set(m))
            {
                co = add_n_masked(c->data, mdata->n->data, m, NWORDS, h);
                // and a third time (b up to a + 3N; the Edwards doubling
                // can need it, see vecsubmod52 in vecarith52.c)
                m = _mm256_andnot_si256(co, m);
                if (any_set(m))
                    add_n_masked(c->data, mdata->n->data, m, NWORDS, h);
            }
        }
    }
}

void vec_simul_addsub52_avx2(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *sum,
    vec_bignum_t *diff, vec_monty_t *mdata)
{
    const __m256i lomask = _mm256_set1_epi64x(0xfffffffffffffLL);
    const __m256i zero = _mm256_setzero_si256();
    const int NWORDS = (int)mdata->NWORDS;
    int h, i;
    for (h = 0; h < 8; h += 4)
    {
        __m256i carry = zero, borrow = zero, cm, bm;
        for (i = 0; i < NWORDS; i++)
        {
            __m256i x = LDH(a->data, i, h), y = LDH(b->data, i, h);
            __m256i ts = _mm256_add_epi64(_mm256_add_epi64(x, y), carry);
            __m256i td = _mm256_add_epi64(_mm256_sub_epi64(x, y), borrow);
            carry = _mm256_srli_epi64(ts, 52);
            borrow = _mm256_cmpgt_epi64(zero, td);
            STH(sum->data, i, h, _mm256_and_si256(ts, lomask));
            STH(diff->data, i, h, _mm256_and_si256(td, lomask));
        }
        cm = _mm256_cmpgt_epi64(carry, zero);   // sum too big: - N
        bm = borrow;                            // diff negative: + N (once)
        if (any_set(cm))
            sub_n_masked(sum->data, mdata->n->data, cm, NWORDS, h);
        if (any_set(bm))
            add_n_masked(diff->data, mdata->n->data, bm, NWORDS, h);
    }
}

#endif
