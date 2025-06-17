/*
MIT License

Copyright (c) 2021 Ben Buhrow

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "soe.h"
#include "soe_impl.h"
#include "ytools.h"
#include <stdint.h>
#if defined(_MSC_VER) && defined(__clang__)
#include <x86intrin.h>
#else
#include <immintrin.h>
#endif
#include "threadpool.h"

uint64_t count_8_bytes(soe_staticdata_t* sdata,
    uint64_t pcount, uint64_t byte_offset);
uint64_t count_8_bytes_bmi2(soe_staticdata_t* sdata,
    uint64_t pcount, uint64_t byte_offset);

typedef __uint128_t uint128_t;

typedef struct
{
    uint64_t data[2][8];
} vec_u104_t;

#define and64 _mm512_and_epi64
#define store64 _mm512_store_epi64
#define storeu64 _mm512_storeu_epi64
#define mstoreu64 _mm512_mask_storeu_epi64
#define storeu512 _mm512_storeu_si512
#define add64 _mm512_add_epi64
#define sub64 _mm512_sub_epi64
#define set64 _mm512_set1_epi64
#define srli64 _mm512_srli_epi64
#define load64 _mm512_load_epi64
#define loadu64 _mm512_loadu_epi64
#define loadu512 _mm512_loadu_si512
#define castpd _mm512_castsi512_pd
#define castepu _mm512_castpd_si512

#define POS_VARIANT
#define LO(x) ((uint64_t)((x) & 0xffffffffffffffffull))
#define HI(x) ((uint64_t)((x) >> 64))


// static __m512i lo52mask;


#ifdef IFMA

#define FORCE_INLINE __inline

FORCE_INLINE __m512i mul52hi(__m512i b, __m512i c)
{
    return _mm512_madd52hi_epu64(_mm512_set1_epi64(0), c, b);
}

FORCE_INLINE __m512i mul52lo(__m512i b, __m512i c)
{
    return _mm512_madd52lo_epu64(_mm512_set1_epi64(0), c, b);
}

FORCE_INLINE void mul52lohi(__m512i b, __m512i c, __m512i* l, __m512i* h)
{
    *l = _mm512_madd52lo_epu64(_mm512_set1_epi64(0), c, b);
    *h = _mm512_madd52hi_epu64(_mm512_set1_epi64(0), c, b);
    return;
}

#else

static __m512d dbias;
static __m512i vbias1;
static __m512i vbias2;
static __m512i vbias3;

__inline __m512i mul52lo(__m512i b, __m512i c)
{
    return _mm512_and_si512(_mm512_mullo_epi64(b, c), _mm512_set1_epi64(0x000fffffffffffffull));
}
__inline __m512i mul52hi(__m512i b, __m512i c)
{
    __m512d prod1_ld = _mm512_cvtepu64_pd(b);
    __m512d prod2_ld = _mm512_cvtepu64_pd(c);
    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
    return _mm512_sub_epi64(castepu(prod1_ld), vbias1);
}
__inline void mul52lohi(__m512i b, __m512i c, __m512i* l, __m512i* h)
{
    __m512d prod1_ld = _mm512_cvtepu64_pd(b);
    __m512d prod2_ld = _mm512_cvtepu64_pd(c);
    __m512d prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
    *h = _mm512_sub_epi64(castepu(prod1_hd), vbias1);
    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
    *l = _mm512_castpd_si512(prod1_ld);
    *l = _mm512_and_si512(*l, lo52mask);
    *h = _mm512_and_si512(*h, lo52mask);
    return;
}

#endif



#ifdef IFMA
#define _mm512_mullo_epi52(c, a, b) \
    c = _mm512_madd52lo_epu64(_mm512_set1_epi64(0), a, b);

#define VEC_MUL_ACCUM_LOHI_PD(a, b, lo, hi) \
    lo = _mm512_madd52lo_epu64(lo, a, b); \
    hi = _mm512_madd52hi_epu64(hi, a, b);
#else

#define VEC_MUL_ACCUM_LOHI_PD(a, b, lo, hi) \
	prod1_ld = _mm512_cvtepu64_pd(a);		\
	prod2_ld = _mm512_cvtepu64_pd(b);		\
    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
    hi = _mm512_add_epi64(hi, _mm512_sub_epi64(castepu(prod1_hd), vbias1)); \
    prod1_hd = _mm512_sub_pd(castpd(vbias2), prod1_hd); \
	prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
	lo = _mm512_add_epi64(lo, _mm512_sub_epi64(castepu(prod1_ld), vbias3));

#define VEC_MUL2_ACCUM_LOHI_PD(c, a, b, lo1, hi1, lo2, hi2) \
	prod1_ld = _mm512_cvtepu64_pd(a);		\
	prod2_ld = _mm512_cvtepu64_pd(b);		\
	prod3_ld = _mm512_cvtepu64_pd(c);		\
    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
	prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
    hi1 = _mm512_add_epi64(hi1, _mm512_sub_epi64(castepu(prod1_hd), vbias1)); \
	hi2 = _mm512_add_epi64(hi2, _mm512_sub_epi64(castepu(prod2_hd), vbias1)); \
    prod1_hd = _mm512_sub_pd(castpd(vbias2), prod1_hd); \
	prod2_hd = _mm512_sub_pd(castpd(vbias2), prod2_hd); \
	prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
	prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
	lo1 = _mm512_add_epi64(lo1, _mm512_sub_epi64(castepu(prod1_ld), vbias3)); \
	lo2 = _mm512_add_epi64(lo2, _mm512_sub_epi64(castepu(prod2_ld), vbias3));

#define _mm512_mullo_epi52(c, a, b) \
    c = _mm512_and_si512(_mm512_mullo_epi64(a, b), _mm512_set1_epi64(0x000fffffffffffffull));
#endif




#define carryprop(lo, hi, mask) \
	{ __m512i carry = _mm512_srli_epi64(lo, 52);	\
	hi = _mm512_add_epi64(hi, carry);		\
	lo = _mm512_and_epi64(mask, lo); }

__m512i __inline my_mm512_addsetc_epi52(__m512i a, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_add_epi64(a, b);
    *cout = _mm512_cmpgt_epu64_mask(t, _mm512_set1_epi64(0xfffffffffffffULL));
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}
__m512i __inline my_mm512_mask_addsetc_epi52(__m512i c, __mmask8 mask, __m512i a, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_add_epi64(a, b);
    *cout = _mm512_mask_cmpgt_epu64_mask(mask, t, _mm512_set1_epi64(0xfffffffffffffULL));
    t = _mm512_mask_and_epi64(c, mask, t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}
__m512i __inline my_mm512_subsetc_epi52(__m512i a, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_sub_epi64(a, b);
    *cout = _mm512_cmpgt_epu64_mask(b, a);
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}
__m512i __inline my_mm512_mask_subsetc_epi52(__m512i c, __mmask8 mask, __m512i a, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_sub_epi64(a, b);
    *cout = _mm512_mask_cmpgt_epu64_mask(mask, b, a);
    t = _mm512_mask_and_epi64(c, mask, t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}
__m512i __inline my_mm512_adc_epi52(__m512i a, __mmask8 c, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_add_epi64(a, b);
    t = _mm512_add_epi64(t, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_cmpgt_epu64_mask(t, _mm512_set1_epi64(0xfffffffffffffULL));
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}
__m512i __inline my_mm512_mask_adc_epi52(__m512i a, __mmask8 m, __mmask8 c, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_add_epi64(a, b);
    t = _mm512_mask_add_epi64(a, m, t, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_cmpgt_epu64_mask(t, _mm512_set1_epi64(0xfffffffffffffULL));
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}
__m512i __inline my_mm512_sbb_epi52(__m512i a, __mmask8 c, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_sub_epi64(a, b);
    *cout = _mm512_cmpgt_epu64_mask(b, a);
    __m512i t2 = _mm512_sub_epi64(t, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_kor(*cout, _mm512_cmpgt_epu64_mask(t2, t));
    t2 = _mm512_and_epi64(t2, _mm512_set1_epi64(0xfffffffffffffULL));
    return t2;
}
__m512i __inline my_mm512_mask_sbb_epi52(__m512i a, __mmask8 m, __mmask8 c, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_mask_sub_epi64(a, m, a, b);
    *cout = _mm512_mask_cmpgt_epu64_mask(m, b, a);
    __m512i t2 = _mm512_mask_sub_epi64(a, m, t, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_kor(*cout, _mm512_mask_cmpgt_epu64_mask(m, t2, t));
    t2 = _mm512_and_epi64(t2, _mm512_set1_epi64(0xfffffffffffffULL));
    return t2;
}

void printvec(char* msg, __m512i v)
{
    ALIGNED_MEM uint64_t m[8];
    store64(m, v);
    int i;
    printf("%s: ", msg);
    for (i = 0; i < 8; i++)
        printf("%016lx ", m[i]);
    printf("\n");
    return;
}


void mulredc52_mask_add_vec(__m512i* c0, __mmask8 addmsk, __m512i a0, __m512i b0, __m512i n0, __m512i vrho)
{
    // CIOS modular multiplication with normal (negative) single-word nhat
    __m512i m;
    __m512i t0, t1, C1;

#ifndef IFMA
    __m512d prod1_hd, prod2_hd;
    __m512d prod1_ld, prod2_ld;
    __m512i i0, i1;
#endif

    __m512i zero = _mm512_set1_epi64(0);
    __m512i one = _mm512_set1_epi64(1);
    __mmask8 scarry2;
    __mmask8 scarry;

    t0 = t1 = C1 = zero;

    VEC_MUL_ACCUM_LOHI_PD(a0, b0, t0, t1);

    // m0
    m = mul52lo(t0, vrho);

    VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, t1);

    // propagate t0's carry before we throw it away.
    // it is almost always exactly 1, but not always, for example when m = 0;
    t0 = _mm512_add_epi64(t1, _mm512_srli_epi64(t0, 52));

    __mmask8 bmsk = _mm512_cmpge_epu64_mask(t0, n0);
    t0 = _mm512_mask_sub_epi64(t0, bmsk, t0, n0);

    // conditional addmod (double result)
    t0 = _mm512_mask_slli_epi64(t0, addmsk, t0, 1);
    bmsk = _mm512_mask_cmpgt_epu64_mask(addmsk, t0, n0);
    *c0 = _mm512_mask_sub_epi64(t0, bmsk & addmsk, t0, n0);
    // _mm512_and_epi64(t0, lo52mask);

    return;
}
void mask_mulredc104_vec(__m512i* c1, __m512i* c0, __mmask8 mulmsk,
    __m512i a1, __m512i a0, __m512i b1, __m512i b0, __m512i n1, __m512i n0, __m512i vrho)
{
    // CIOS modular multiplication with normal (negative) single-word nhat
    __m512i m;
    __m512i t0, t1, t2, t3, C1, C2;

#ifndef IFMA
    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
    int biascount = 0;
    __m512i i0, i1;
#endif

    __m512i zero = _mm512_set1_epi64(0);
    __m512i one = _mm512_set1_epi64(1);
    __m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
    __mmask8 scarry2;
    __mmask8 scarry;

    t0 = t1 = t2 = t3 = C1 = C2 = zero;

    VEC_MUL_ACCUM_LOHI_PD(a0, b0, t0, t1);
    VEC_MUL_ACCUM_LOHI_PD(a1, b0, t1, t2);
    //VEC_MUL2_ACCUM_LOHI_PD(b0, a0, a1, t0, t1, C1, t2);
    //t1 = _mm512_add_epi64(t1, C1);

    // m0
    m = mul52lo(t0, vrho);

    VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
    VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, C2);
    //VEC_MUL2_ACCUM_LOHI_PD(m, n0, n1, t0, C1, t1, C2);

    t1 = _mm512_add_epi64(t1, C1);
    t2 = _mm512_add_epi64(t2, C2);
    // we throw t0 away after this so first propagate its carry.
    t0 = _mm512_add_epi64(t1, _mm512_srli_epi64(t0, 52));
    t1 = t2;
    t2 = C1 = zero;

    VEC_MUL_ACCUM_LOHI_PD(a0, b1, t0, C1);
    VEC_MUL_ACCUM_LOHI_PD(a1, b1, t1, t2);
    //VEC_MUL2_ACCUM_LOHI_PD(b1, a0, a1, t0, C1, t1, t2);

    t1 = _mm512_add_epi64(t1, C1);
    C1 = C2 = zero;

    // m1
    m = mul52lo(t0, vrho);

    VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
    VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, t2);
    //VEC_MUL2_ACCUM_LOHI_PD(m, n0, n1, t0, C1, t1, t2);

    t1 = _mm512_add_epi64(t1, C1);

    // final carryprop
    carryprop(t0, t1, lo52mask);
    carryprop(t1, t2, lo52mask);
    carryprop(t2, C2, lo52mask);

    scarry = _mm512_cmp_epu64_mask(C2, zero, _MM_CMPINT_GT);

    if (scarry > 0) {
        // conditionally subtract when needed (AMM - only on overflow)
        C1 = _mm512_mask_set1_epi64(zero, _mm512_cmpgt_epi64_mask(n0, t1), 1);
        t1 = _mm512_mask_sub_epi64(t1, scarry, t1, n0);
        t2 = _mm512_mask_sub_epi64(t2, scarry, t2, n1);
        t2 = _mm512_mask_sub_epi64(t2, scarry, t2, C1);
    }

    // on Zen4-epyc it is slower to do this:
    // *c0 = _mm512_mask_and_epi64(a0, mulmsk, lo52mask, t1);
    // *c1 = _mm512_mask_and_epi64(a1, mulmsk, lo52mask, t2);

    // than this:
    *c0 = _mm512_and_epi64(lo52mask, t1);
    *c1 = _mm512_and_epi64(lo52mask, t2);
    *c0 = _mm512_mask_mov_epi64(*c0, ~mulmsk, a0);
    *c1 = _mm512_mask_mov_epi64(*c1, ~mulmsk, a1);

    return;
}
void sqrredc104_vec(__m512i* c1, __m512i* c0,
    __m512i a1, __m512i a0, __m512i n1, __m512i n0, __m512i vrho)
{
    // CIOS modular multiplication with normal (negative) single-word nhat
    __m512i m;
    __m512i t0, t1, t2, t3, C1, C2, sqr_lo, sqr_hi;

#ifndef IFMA
    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
    int biascount = 0;
    __m512i i0, i1;
#endif

    __m512i zero = _mm512_set1_epi64(0);
    __m512i one = _mm512_set1_epi64(1);
    __m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
    __mmask8 scarry2;
    __mmask8 scarry;

    t0 = t1 = t2 = t3 = C1 = C2 = sqr_lo = sqr_hi = zero;

    VEC_MUL_ACCUM_LOHI_PD(a1, a0, sqr_lo, sqr_hi);
    t1 = sqr_lo;
    t2 = sqr_hi;
    VEC_MUL_ACCUM_LOHI_PD(a0, a0, t0, t1);

    // m0
    m = mul52lo(t0, vrho);

    VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
    VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, C2);

    t1 = _mm512_add_epi64(t1, C1);
    t2 = _mm512_add_epi64(t2, C2);
    // we throw t0 away after this so first propagate its carry.
    // t0 = _mm512_add_epi64(t1, one);
    // t0 = _mm512_mask_add_epi64(t1, _mm512_cmpgt_epu64_mask(m, zero), t1, one);
    t0 = _mm512_add_epi64(t1, _mm512_srli_epi64(t0, 52));
    t1 = t2;
    t2 = C1 = C2 = zero;

    VEC_MUL_ACCUM_LOHI_PD(a1, a1, t1, t2);

    t0 = _mm512_add_epi64(t0, sqr_lo);
    t1 = _mm512_add_epi64(t1, sqr_hi);

    // m1
    m = mul52lo(t0, vrho);

    VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
    VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, t2);

    t1 = _mm512_add_epi64(t1, C1);

    // final carryprop
    carryprop(t0, t1, lo52mask);
    carryprop(t1, t2, lo52mask);
    carryprop(t2, C2, lo52mask);

    scarry = _mm512_cmp_epu64_mask(C2, zero, _MM_CMPINT_GT);

    if (scarry > 0) {
        // conditionally subtract when needed (AMM - only on overflow)
        __mmask8 bmsk;
        bmsk = _mm512_cmpgt_epi64_mask(n0, t1);
        t1 = _mm512_mask_sub_epi64(t1, scarry, t1, n0);
        t2 = _mm512_mask_sub_epi64(t2, scarry, t2, n1);
        t2 = _mm512_mask_sub_epi64(t2, scarry, t2, _mm512_mask_set1_epi64(zero, bmsk, 1));
    }
    *c0 = _mm512_and_epi64(lo52mask, t1);
    *c1 = _mm512_and_epi64(lo52mask, t2);

    return;
}
void mask_sqrredc104_vec(__m512i* c1, __m512i* c0, __mmask8 mulmsk,
    __m512i a1, __m512i a0, __m512i n1, __m512i n0, __m512i vrho)
{
    // CIOS modular multiplication with normal (negative) single-word nhat
    __m512i m;
    __m512i t0, t1, t2, C3, C1, C2, sqr_lo, sqr_hi;

#ifndef IFMA
    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
    int biascount = 0;
    __m512i i0, i1;
#endif

    __m512i zero = _mm512_set1_epi64(0);
    __m512i one = _mm512_set1_epi64(1);
    __m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
    __mmask8 scarry2;
    __mmask8 scarry;

    t0 = t1 = t2 = C1 = C2 = C3 = sqr_lo = sqr_hi = zero;

    VEC_MUL_ACCUM_LOHI_PD(a1, a0, sqr_lo, sqr_hi);
    VEC_MUL_ACCUM_LOHI_PD(a0, a0, t0, t1);

    // m0
    t1 += sqr_lo;
    m = mul52lo(t0, vrho);

    VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
    VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, C2);

    t1 = _mm512_add_epi64(t1, C1);
    t2 = _mm512_add_epi64(sqr_hi, C2);
    // we throw t0 away after this so first propagate its carry.
    // t0 = _mm512_add_epi64(t1, one);
    // t0 = _mm512_mask_add_epi64(t1, _mm512_cmpgt_epu64_mask(m, zero), t1, one);
    t0 = _mm512_add_epi64(t1, _mm512_srli_epi64(t0, 52));
    t1 = t2;
    t2 = C1 = C2 = zero;

    VEC_MUL_ACCUM_LOHI_PD(a1, a1, t1, t2);

    t0 = _mm512_add_epi64(t0, sqr_lo);
    t1 = _mm512_add_epi64(t1, sqr_hi);

    // m1
    m = mul52lo(t0, vrho);

    VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
    VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, t2);

    t1 = _mm512_add_epi64(t1, C1);

    // final carryprop
    carryprop(t0, t1, lo52mask);
    carryprop(t1, t2, lo52mask);

    carryprop(t2, C2, lo52mask);
    scarry = _mm512_cmp_epu64_mask(C2, zero, _MM_CMPINT_GT);

    //scarry = _mm512_cmpge_epu64_mask(t2, n1);

    if (scarry > 0) {
        // conditionally subtract when needed (AMM - only on overflow)
        C1 = _mm512_mask_set1_epi64(zero, _mm512_cmpgt_epi64_mask(n0, t1), 1);
        t1 = _mm512_mask_sub_epi64(t1, scarry, t1, n0);
        t2 = _mm512_mask_sub_epi64(t2, scarry, t2, n1);
        t2 = _mm512_mask_sub_epi64(t2, scarry, t2, C1);
    }

    // on Zen4-epyc it is slower to do this:
    // *c0 = _mm512_mask_and_epi64(a0, mulmsk, lo52mask, t1);
    // *c1 = _mm512_mask_and_epi64(a1, mulmsk, lo52mask, t2);

    // than this:
    *c0 = _mm512_and_epi64(lo52mask, t1);
    *c1 = _mm512_and_epi64(lo52mask, t2);
    *c0 = _mm512_mask_mov_epi64(*c0, ~mulmsk, a0);
    *c1 = _mm512_mask_mov_epi64(*c1, ~mulmsk, a1);

    return;
}
void mask_sqrredc104_vec_pos(__m512i* c1, __m512i* c0, __mmask8 mulmsk,
    __m512i a1, __m512i a0, __m512i n1, __m512i n0, __m512i vrho)
{
    // CIOS modular multiplication with positive variant n'
    __m512i m;
    __m512i t0, t1, t2, C3, C1, C2, sqr_lo, sqr_hi;

#ifndef IFMA
    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
    int biascount = 0;
    __m512i i0, i1;
#endif

    __m512i zero = _mm512_set1_epi64(0);
    __m512i one = _mm512_set1_epi64(1);
    __m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
    __mmask8 scarry2;
    __mmask8 scarry;

    t0 = t1 = t2 = C1 = C2 = C3 = sqr_lo = sqr_hi = zero;

    VEC_MUL_ACCUM_LOHI_PD(a1, a0, sqr_lo, sqr_hi);
    VEC_MUL_ACCUM_LOHI_PD(a0, a0, t0, t1);

    // m0
    t1 += sqr_lo;

    // note, we leave rho = 0 - rho so that we get -m,
    // and thus the muladd becomes mulsub, since there is
    // no fmsub52 in avx512-ifma.
    m = mul52lo(t0, vrho);

    VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
    VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, C2);

    t1 = _mm512_sub_epi64(t1, C1);
    t2 = _mm512_add_epi64(sqr_hi, C2);
    t0 = t1;
    t1 = t2;
    t2 = C1 = C2 = zero;

    VEC_MUL_ACCUM_LOHI_PD(a1, a1, t1, t2);

    t0 = _mm512_add_epi64(t0, sqr_lo);
    t1 = _mm512_add_epi64(t1, sqr_hi);

    // m1
    m = mul52lo(t0, vrho);

    VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
    VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, t2);

    t1 = _mm512_sub_epi64(t1, C1);

    // final carryprop
    carryprop(t0, t1, lo52mask);
    carryprop(t1, t2, lo52mask);
    carryprop(t2, C2, lo52mask);
    scarry = _mm512_cmp_epu64_mask(C2, zero, _MM_CMPINT_GT);

    if (scarry > 0) {
        // conditionally add
        t1 = _mm512_mask_add_epi64(t1, scarry, t1, n0);
        t2 = _mm512_mask_add_epi64(t2, scarry, t2, n1);
        t2 = _mm512_mask_add_epi64(t2, scarry, t2, _mm512_srli_epi64(t1, 52));
    }

    // on Zen4-epyc it is slower to do this:
    // *c0 = _mm512_mask_and_epi64(a0, mulmsk, lo52mask, t1);
    // *c1 = _mm512_mask_and_epi64(a1, mulmsk, lo52mask, t2);

    // than this:
    *c0 = _mm512_and_epi64(lo52mask, t1);
    *c1 = _mm512_and_epi64(lo52mask, t2);
    *c0 = _mm512_mask_mov_epi64(*c0, ~mulmsk, a0);
    *c1 = _mm512_mask_mov_epi64(*c1, ~mulmsk, a1);

    return;
}
void mask_sqrredc104_exact_vec(__m512i* c1, __m512i* c0, __mmask8 mulmsk,
    __m512i a1, __m512i a0, __m512i n1, __m512i n0, __m512i vrho)
{
    // CIOS modular multiplication with normal (negative) single-word nhat
    __m512i m;
    __m512i t0, t1, t2, C3, C1, C2, sqr_lo, sqr_hi;

#ifndef IFMA
    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
    int biascount = 0;
    __m512i i0, i1;
#endif

    __m512i zero = _mm512_set1_epi64(0);
    __m512i one = _mm512_set1_epi64(1);
    __m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
    __mmask8 scarry2;
    __mmask8 scarry;

    t0 = t1 = t2 = C1 = C2 = C3 = sqr_lo = sqr_hi = zero;

    VEC_MUL_ACCUM_LOHI_PD(a1, a0, sqr_lo, sqr_hi);
    t1 = sqr_lo;
    t2 = sqr_hi;
    VEC_MUL_ACCUM_LOHI_PD(a0, a0, t0, t1);

    // m0
    m = mul52lo(t0, vrho);

    VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
    VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, C2);

    t1 = _mm512_add_epi64(t1, C1);
    t2 = _mm512_add_epi64(t2, C2);
    // we throw t0 away after this so first propagate its carry.
    // t0 = _mm512_add_epi64(t1, one);
    // t0 = _mm512_mask_add_epi64(t1, _mm512_cmpgt_epu64_mask(m, zero), t1, one);
    t0 = _mm512_add_epi64(t1, _mm512_srli_epi64(t0, 52));
    t1 = t2;
    t2 = C1 = C2 = zero;

    VEC_MUL_ACCUM_LOHI_PD(a1, a1, t1, t2);

    t0 = _mm512_add_epi64(t0, sqr_lo);
    t1 = _mm512_add_epi64(t1, sqr_hi);

    // m1
    m = mul52lo(t0, vrho);

    VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
    VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, t2);

    t1 = _mm512_add_epi64(t1, C1);

    // final carryprop
    carryprop(t0, t1, lo52mask);
    carryprop(t1, t2, lo52mask);

    //carryprop(t2, C2, lo52mask);
    //scarry = _mm512_cmp_epu64_mask(C2, zero, _MM_CMPINT_GT);

    scarry = _mm512_cmpge_epu64_mask(t2, n1);
    //scarry |= (_mm512_cmpeq_epu64_mask(t2, n1) & _mm512_cmpgt_epu64_mask(t1, n0));

    if (scarry > 0) {
        // conditionally subtract when result >= n
        C1 = _mm512_mask_set1_epi64(zero, _mm512_cmpgt_epi64_mask(n0, t1), 1);
        t1 = _mm512_mask_sub_epi64(t1, scarry, t1, n0);
        t2 = _mm512_mask_sub_epi64(t2, scarry, t2, n1);
        t2 = _mm512_mask_sub_epi64(t2, scarry, t2, C1);
    }

    // on Zen4-epyc it is slower to do this:
    // *c0 = _mm512_mask_and_epi64(a0, mulmsk, lo52mask, t1);
    // *c1 = _mm512_mask_and_epi64(a1, mulmsk, lo52mask, t2);

    // than this:
    *c0 = _mm512_and_epi64(lo52mask, t1);
    *c1 = _mm512_and_epi64(lo52mask, t2);
    *c0 = _mm512_mask_mov_epi64(*c0, ~mulmsk, a0);
    *c1 = _mm512_mask_mov_epi64(*c1, ~mulmsk, a1);

    return;
}
void addmod104_x8(__m512i* c1, __m512i* c0, __m512i a1, __m512i a0,
    __m512i b1, __m512i b0, __m512i n1, __m512i n0)
{
    __m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);

    // add
    __mmask8 bmsk;
    //a0 = my_mm512_addsetc_epi52(a0, b0, &bmsk);
    //a1 = my_mm512_adc_epi52(a1, bmsk, b1, &bmsk);
    a0 = _mm512_add_epi64(a0, b0);
    a1 = _mm512_add_epi64(a1, b1);
    a1 = _mm512_add_epi64(a1, _mm512_srli_epi64(a0, 52));
    a0 = _mm512_and_epi64(a0, lo52mask);

    // compare
    //__mmask8 msk = bmsk | _mm512_cmpgt_epu64_mask(a1, n1);
    __mmask8 msk = _mm512_cmpgt_epu64_mask(a1, n1);
    msk |= (_mm512_cmpeq_epu64_mask(a1, n1) & _mm512_cmpge_epu64_mask(a0, n0));

    // conditionally subtract N
    *c0 = my_mm512_mask_subsetc_epi52(a0, msk, a0, n0, &bmsk);
    *c1 = my_mm512_mask_sbb_epi52(a1, msk, bmsk, n1, &bmsk);
    // *c0 = _mm512_mask_sub_epi64(a0, msk, a0, n0);
    // *c1 = _mm512_mask_sub_epi64(a1, msk, a1, n1);
    // *c1 = _mm512_mask_sub_epi64(*c1, msk, *c1, _mm512_srli_epi64(*c0, 63));
    return;
}
void mask_addmod104_x8(__m512i* c1, __m512i* c0, __mmask8 addmsk,
    __m512i a1, __m512i a0, __m512i b1, __m512i b0, __m512i n1, __m512i n0)
{
    // add
    __mmask8 bmsk;
    a0 = my_mm512_mask_addsetc_epi52(a0, addmsk, a0, b0, &bmsk);
    a1 = my_mm512_mask_adc_epi52(a1, addmsk, bmsk, b1, &bmsk);

    // compare
    __mmask8 msk = bmsk | _mm512_cmpgt_epu64_mask(a1, n1);
    msk |= (_mm512_cmpeq_epu64_mask(a1, n1) & _mm512_cmpge_epu64_mask(a0, n0));

    // conditionally subtract N
    *c0 = my_mm512_mask_subsetc_epi52(a0, addmsk & msk, a0, n0, &bmsk);
    *c1 = my_mm512_mask_sbb_epi52(a1, addmsk & msk, bmsk, n1, &bmsk);
    return;
}
void mask_dblmod104_x8(__m512i* c1, __m512i* c0, __mmask8 addmsk,
    __m512i a1, __m512i a0, __m512i n1, __m512i n0)
{
    __m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);

    // add
    __mmask8 bmsk;
    //a0 = my_mm512_mask_addsetc_epi52(a0, addmsk, a0, b0, &bmsk);
    //a1 = my_mm512_mask_adc_epi52(a1, addmsk, bmsk, b1, &bmsk);

    a0 = _mm512_mask_slli_epi64(a0, addmsk, a0, 1);
    a1 = _mm512_mask_slli_epi64(a1, addmsk, a1, 1);
    // when doubling, it is safe to check both carries before adding
    // in the previous carry, because the shift makes room for
    // the previous carry.  So either the upper word shift generates
    // a carry or doesn't, the addition won't cause one.
    a1 = _mm512_add_epi64(a1, _mm512_srli_epi64(a0, 52));
    a0 = _mm512_and_epi64(lo52mask, a0);

    // compare
    __mmask8 msk = _mm512_cmpgt_epu64_mask(a1, n1);
    msk |= (_mm512_cmpeq_epu64_mask(a1, n1) & _mm512_cmpge_epu64_mask(a0, n0));

    // conditionally subtract N
    *c0 = my_mm512_mask_subsetc_epi52(a0, addmsk & msk, a0, n0, &bmsk);
    *c1 = my_mm512_mask_sbb_epi52(a1, addmsk & msk, bmsk, n1, &bmsk);
    return;
}
void mask_redsub104_x8(__m512i* c1, __m512i* c0, __mmask8 addmsk,
    __m512i a1, __m512i a0, __m512i n1, __m512i n0)
{
    __mmask8 bmsk;

    // compare
    __mmask8 msk = _mm512_cmpgt_epu64_mask(a1, n1);
    msk |= (_mm512_cmpeq_epu64_mask(a1, n1) & _mm512_cmpge_epu64_mask(a0, n0));

    // conditionally subtract N
    *c0 = my_mm512_mask_subsetc_epi52(a0, addmsk & msk, a0, n0, &bmsk);
    *c1 = my_mm512_mask_sbb_epi52(a1, addmsk & msk, bmsk, n1, &bmsk);
    return;
}
void redsub104_x8(__m512i* c1, __m512i* c0,
    __m512i a1, __m512i a0, __m512i n1, __m512i n0)
{
    __mmask8 bmsk;

    // compare
    __mmask8 msk = _mm512_cmpgt_epu64_mask(a1, n1);
    msk |= (_mm512_cmpeq_epu64_mask(a1, n1) & _mm512_cmpge_epu64_mask(a0, n0));

    // conditionally subtract N
    *c0 = my_mm512_mask_subsetc_epi52(a0, msk, a0, n0, &bmsk);
    *c1 = my_mm512_mask_sbb_epi52(a1, msk, bmsk, n1, &bmsk);

    return;
}
void submod104_x8(__m512i* c1, __m512i* c0, __m512i a1, __m512i a0,
    __m512i b1, __m512i b0, __m512i n1, __m512i n0)
{
    // compare
    __mmask8 msk = _mm512_cmplt_epu64_mask(a1, b1);
    msk |= _mm512_cmpeq_epu64_mask(a1, b1) & _mm512_cmplt_epu64_mask(a0, b0);

    // subtract
    __mmask8 bmsk;
    a0 = my_mm512_subsetc_epi52(a0, b0, &bmsk);
    a1 = my_mm512_sbb_epi52(a1, bmsk, b1, &bmsk);

    // conditionally add N
    *c0 = my_mm512_mask_addsetc_epi52(a0, msk, a0, n0, &bmsk);
    *c1 = my_mm512_mask_adc_epi52(a1, msk, bmsk, n1, &bmsk);
    return;
}
static __inline __m512i multiplicative_inverse104_x8(uint64_t* a)
{
    //    assert(a%2 == 1);  // the inverse (mod 2<<64) only exists for odd values
    __m512i x0, x1, x2, x3, x4, x5, y, n, i0, i1;
    __m512i three = _mm512_set1_epi64(3), two = _mm512_set1_epi64(2), one = _mm512_set1_epi64(1);
    __m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);

    n = _mm512_load_epi64(a);
    _mm512_mullo_epi52(x0, n, three);
    x0 = _mm512_xor_epi64(x0, two);
    _mm512_mullo_epi52(y, n, x0);
    y = _mm512_sub_epi64(one, y);
    y = _mm512_and_epi64(lo52mask, y);

    x1 = _mm512_add_epi64(y, one);
    _mm512_mullo_epi52(x1, x0, x1);
    _mm512_mullo_epi52(y, y, y);

    x2 = _mm512_add_epi64(y, one);
    _mm512_mullo_epi52(x2, x1, x2);
    _mm512_mullo_epi52(y, y, y);

    x3 = _mm512_add_epi64(y, one);
    _mm512_mullo_epi52(x3, x2, x3);
    _mm512_mullo_epi52(y, y, y);

    x4 = _mm512_add_epi64(y, one);
    _mm512_mullo_epi52(x4, x3, x4);

    return x4;
}

#if defined(INTEL_COMPILER) || defined(INTEL_LLVM_COMPILER)
#define ROUNDING_MODE (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)
#else
#define ROUNDING_MODE _MM_FROUND_CUR_DIRECTION
#endif

#ifdef _MSC_VER
#  define FORCE_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER) || defined (__INTEL_LLVM_COMPILER)
#  define FORCE_INLINE inline __attribute__((always_inline))
#else
#  define FORCE_INLINE __inline
#endif


// a Miller-Rabin SPRP test on 8x 104-bit inputs using base 2
uint8_t MR_2sprp_104x8(uint64_t* n)
{
    // assumes has no small factors.  assumes n >= 54 bits.
    // assumes n is a list of 8 104-bit integers (16 52-bit words)
    // in the format: 8 lo-words, 8 hi-words.
    // do a Miller-Rabin sprp test using base 2.
    __m512i vrho = multiplicative_inverse104_x8(n);
    __m512i mone[2];
    vec_u104_t unity;
    __m512i nv[2];
    __m512i dv[2];
    __m512i rv[2];
    __m512i bv[2];
    __m512i n1v[2];
    __m512i tv[2];
    __m512i m;
    __m512i zerov = _mm512_setzero_si512();
    __m512i onev = _mm512_set1_epi64(1);
    __m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
    uint64_t tmp = 0;

    vrho = _mm512_and_epi64(_mm512_sub_epi64(zerov, vrho), lo52mask);

    nv[0] = _mm512_load_epi64(&n[0]);
    nv[1] = _mm512_load_epi64(&n[8]);

    // the 128-bit division we do one at a time
    int i;
    for (i = 0; i < 8; i++)
    {
        uint128_t mod = ((uint128_t)n[i + 8] << 52) + n[i];
        uint128_t one = (uint128_t)1 << 104;
        one %= mod;

        unity.data[0][i] = (uint64_t)one & 0xfffffffffffffULL;
        unity.data[1][i] = (uint64_t)(one >> 52) & 0xfffffffffffffULL;
    }

    mone[0] = _mm512_load_epi64(unity.data[0]);
    mone[1] = _mm512_load_epi64(unity.data[1]);

    // compute d and tzcnt
    submod104_x8(&n1v[1], &n1v[0], nv[1], nv[0], zerov, onev, nv[1], nv[0]);

    __mmask8 done = 0;
    dv[1] = n1v[1];
    dv[0] = n1v[0];
    __m512i tzcntv = zerov;
    while (done != 0xff)
    {
        __m512i c = _mm512_mask_slli_epi64(dv[1], ~done, dv[1], 51);
        dv[0] = _mm512_mask_srli_epi64(dv[0], ~done, dv[0], 1);
        dv[0] = _mm512_mask_or_epi64(dv[0], ~done, c, dv[0]);
        dv[1] = _mm512_mask_srli_epi64(dv[1], ~done, dv[1], 1);
        tzcntv = _mm512_mask_add_epi64(tzcntv, ~done, tzcntv, onev);
        done = done | _mm512_cmpeq_epi64_mask(_mm512_and_epi64(dv[0], onev), onev);
    }
    dv[0] = _mm512_and_epi64(dv[0], lo52mask);

    // penultimate-hi-bit mask based on d
    m = _mm512_sub_epi64(_mm512_set1_epi64(62), _mm512_lzcnt_epi64(dv[1]));
    m = _mm512_sllv_epi64(_mm512_set1_epi64(1), m);
    m = _mm512_mask_set1_epi64(m, _mm512_cmple_epi64_mask(dv[1], onev), 0);

    // we know the first bit is set and the first squaring is of unity,
    // so we can do the first iteration manually (and hence the penultimate mask bit)
    addmod104_x8(&rv[1], &rv[0], mone[1], mone[0], mone[1], mone[0], nv[1], nv[0]);

    __mmask8 protect = _mm512_cmpgt_epi64_mask(_mm512_srli_epi64(n1v[1], 49), zerov);

    // compute b^d
    if (protect)
    {
        done = _mm512_cmpeq_epu64_mask(m, zerov);
        while (done != 0xff)
        {
            __mmask8 bitcmp = _mm512_test_epi64_mask(m, dv[1]);

            mask_sqrredc104_exact_vec(&rv[1], &rv[0], ~done, rv[1], rv[0], nv[1], nv[0], vrho);
            mask_dblmod104_x8(&rv[1], &rv[0], (~done) & bitcmp, rv[1], rv[0], nv[1], nv[0]);

            m = _mm512_srli_epi64(m, 1);
            done = _mm512_cmpeq_epu64_mask(m, zerov);
        }
    }
    else
    {
        done = _mm512_cmpeq_epu64_mask(m, zerov);
        while (done != 0xff)
        {
            __mmask8 bitcmp = _mm512_test_epi64_mask(m, dv[1]);

            mask_sqrredc104_vec(&rv[1], &rv[0], ~done, rv[1], rv[0], nv[1], nv[0], vrho);
            mask_dblmod104_x8(&rv[1], &rv[0], (~done) & bitcmp, rv[1], rv[0], nv[1], nv[0]);

            m = _mm512_srli_epi64(m, 1);
            done = _mm512_cmpeq_epu64_mask(m, zerov);
        }
    }

    m = _mm512_sub_epi64(_mm512_set1_epi64(62), _mm512_lzcnt_epi64(dv[0]));
    m = _mm512_sllv_epi64(_mm512_set1_epi64(1), m);
    m = _mm512_mask_set1_epi64(m, _mm512_cmpge_epi64_mask(dv[1], onev), 1ULL << 51);

    if (protect)
    {
        done = _mm512_cmpeq_epu64_mask(m, zerov);
        while (done != 0xff)
        {
            __mmask8 bitcmp = _mm512_test_epi64_mask(m, dv[0]);

            mask_sqrredc104_exact_vec(&rv[1], &rv[0], ~done, rv[1], rv[0], nv[1], nv[0], vrho);
            mask_dblmod104_x8(&rv[1], &rv[0], (~done) & bitcmp, rv[1], rv[0], nv[1], nv[0]);

            m = _mm512_srli_epi64(m, 1);
            done = _mm512_cmpeq_epu64_mask(m, zerov);
        }
    }
    else
    {
        done = _mm512_cmpeq_epu64_mask(m, zerov);
        while (done != 0xff)
        {
            __mmask8 bitcmp = _mm512_test_epi64_mask(m, dv[0]);

            mask_sqrredc104_vec(&rv[1], &rv[0], ~done, rv[1], rv[0], nv[1], nv[0], vrho);
            mask_dblmod104_x8(&rv[1], &rv[0], (~done) & bitcmp, rv[1], rv[0], nv[1], nv[0]);

            m = _mm512_srli_epi64(m, 1);
            done = _mm512_cmpeq_epu64_mask(m, zerov);
        }
    }

    // AMM possibly needs a final correction by n
    addmod104_x8(&rv[1], &rv[0], zerov, zerov, rv[1], rv[0], nv[1], nv[0]);

    // check current result == 1
    __mmask8 is1prp = _mm512_cmpeq_epu64_mask(rv[1], mone[1]) &
        _mm512_cmpeq_epu64_mask(rv[0], mone[0]);

    // now compute b^(2^s*d) and check for congruence to -1 as we go.
    // check while tzcnt is > 1 for all inputs or all are already not prp.
    done = is1prp;
    __mmask8 ism1prp = 0;

    submod104_x8(&n1v[1], &n1v[0], zerov, zerov, mone[1], mone[0], nv[1], nv[0]);

    while (done != 0xff)
    {
        tzcntv = _mm512_mask_sub_epi64(tzcntv, ~done, tzcntv, onev);

        // prp by -1 check
        ism1prp = (_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[1], n1v[1]) &
            _mm512_mask_cmpeq_epu64_mask(~is1prp, rv[0], n1v[0]));

        is1prp |= ism1prp;	// stop checking it if we've found a prp criteria.
        done = (is1prp | ism1prp);

        sqrredc104_vec(&rv[1], &rv[0], rv[1], rv[0], nv[1], nv[0], vrho);
        addmod104_x8(&rv[1], &rv[0], zerov, zerov, rv[1], rv[0], nv[1], nv[0]);

        // definitely not prp by 1 check, stop checking
        done |= (_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[1], zerov) &
            _mm512_mask_cmpeq_epu64_mask(~is1prp, rv[0], onev));

        done |= _mm512_mask_cmple_epu64_mask(~done, tzcntv, onev);
    }

    addmod104_x8(&rv[1], &rv[0], zerov, zerov, rv[1], rv[0], nv[1], nv[0]);

    // check current result == m-1
    ism1prp |= (_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[1], n1v[1]) &
        _mm512_mask_cmpeq_epu64_mask(~is1prp, rv[0], n1v[0]));

    return (is1prp | ism1prp);
}

uint8_t fermat_prp_104x8(uint64_t* n)
{
    // assumes has no small factors.  assumes n >= 54 bits.
    // assumes n is a list of 8 104-bit integers (16 52-bit words)
    // in the format: 8 lo-words, 8 hi-words.
    // do a base-2 fermat prp test on each using LR binexp.
    __m512i vrho = multiplicative_inverse104_x8(n);
    vec_u104_t unity;
    __m512i nvec[2];
    __m512i evec[2];
    __m512i m;
    __m512i zero = _mm512_setzero_si512();
    __m512i one = _mm512_set1_epi64(1);
    __m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
    uint64_t tmp = 0;

    vrho = _mm512_and_epi64(_mm512_sub_epi64(zero, vrho), lo52mask);

    nvec[0] = loadu64(&n[0]);
    nvec[1] = loadu64(&n[8]);
    submod104_x8(&evec[1], &evec[0],
        nvec[1], nvec[0], zero, one, nvec[1], nvec[0]);

    // the 128-bit division we do the slow way
    int i;
    for (i = 0; i < 8; i++)
    {
        uint128_t mod = ((uint128_t)n[i + 8] << 52) + n[i];
        uint128_t t = (uint128_t)1 << 104;
        t %= mod;

        unity.data[0][i] = (uint64_t)t & 0xfffffffffffffULL;
        unity.data[1][i] = (uint64_t)(t >> 52);
    }

    // penultimate-hi-bit mask
    m = _mm512_sub_epi64(_mm512_set1_epi64(62), _mm512_lzcnt_epi64(evec[1]));
    m = _mm512_sllv_epi64(_mm512_set1_epi64(1), m);
    m = _mm512_mask_set1_epi64(m, _mm512_cmple_epu64_mask(evec[1], one), 0);

    __mmask8 protect = _mm512_cmpgt_epi64_mask(_mm512_srli_epi64(evec[1], 49), zero);

    // we know the first bit is set and the first squaring is of unity,
    // so we can do the first iteration manually with no squaring.
    // Note: the first 5 iterations can be done much more cheaply in
    // single precision and then converted into montgomery representation,
    // but that would require a 208-bit division; not worth it.
    __m512i r0 = loadu64(unity.data[0]);
    __m512i r1 = loadu64(unity.data[1]);

    addmod104_x8(&r1, &r0, r1, r0, r1, r0, nvec[1], nvec[0]);

    __mmask8 done;
    if (protect)
    {
        done = _mm512_cmpeq_epu64_mask(m, zero);
        while (done != 0xff)
        {
            __mmask8 bitcmp = _mm512_test_epi64_mask(m, evec[1]);

            mask_sqrredc104_exact_vec(&r1, &r0, ~done, r1, r0, nvec[1], nvec[0], vrho);
            mask_dblmod104_x8(&r1, &r0, (~done) & bitcmp, r1, r0, nvec[1], nvec[0]);

            m = _mm512_srli_epi64(m, 1);
            done = _mm512_cmpeq_epu64_mask(m, zero);
        }
    }
    else
    {

        __mmask8 done = _mm512_cmpeq_epu64_mask(m, zero);
        while (done != 0xff)
        {
            __mmask8 bitcmp = _mm512_test_epi64_mask(m, evec[1]);

            mask_sqrredc104_vec(&r1, &r0, ~done, r1, r0, nvec[1], nvec[0], vrho);
            mask_dblmod104_x8(&r1, &r0, (~done) & bitcmp, r1, r0, nvec[1], nvec[0]);

            m = _mm512_srli_epi64(m, 1);
            done = _mm512_cmpeq_epu64_mask(m, zero);
        }
    }

    m = _mm512_sub_epi64(_mm512_set1_epi64(62), _mm512_lzcnt_epi64(evec[0]));
    m = _mm512_sllv_epi64(_mm512_set1_epi64(1), m);
    m = _mm512_mask_set1_epi64(m, _mm512_cmpge_epu64_mask(evec[1], one), 1ULL << 51);

    if (protect)
    {
        done = _mm512_cmpeq_epu64_mask(m, zero);
        while (done != 0xff)
        {
            __mmask8 bitcmp = _mm512_test_epi64_mask(m, evec[0]);

            mask_sqrredc104_exact_vec(&r1, &r0, ~done, r1, r0, nvec[1], nvec[0], vrho);
            mask_dblmod104_x8(&r1, &r0, (~done) & bitcmp, r1, r0, nvec[1], nvec[0]);

            m = _mm512_srli_epi64(m, 1);
            done = _mm512_cmpeq_epu64_mask(m, zero);
        }
    }
    else
    {

        __mmask8 done = _mm512_cmpeq_epu64_mask(m, zero);
        while (done != 0xff)
        {
            __mmask8 bitcmp = _mm512_test_epi64_mask(m, evec[0]);

            mask_sqrredc104_vec(&r1, &r0, ~done, r1, r0, nvec[1], nvec[0], vrho);
            mask_dblmod104_x8(&r1, &r0, (~done) & bitcmp, r1, r0, nvec[1], nvec[0]);

            m = _mm512_srli_epi64(m, 1);
            done = _mm512_cmpeq_epu64_mask(m, zero);
        }
    }

    // AMM possibly needs a final correction by n
    addmod104_x8(&r1, &r0, r1, r0, zero, zero, nvec[1], nvec[0]);

    uint8_t isprp =
        _mm512_cmpeq_epu64_mask(loadu64(unity.data[0]), r0) &
        _mm512_cmpeq_epu64_mask(loadu64(unity.data[1]), r1);

    return isprp;
}

void count_twins_dispatch(void* vptr)
{
    tpool_t* tdata = (tpool_t*)vptr;
    soe_userdata_t* t = (soe_userdata_t*)tdata->user_data;
    soe_staticdata_t* sdata = t->sdata;

    // launch one range of computation for each thread.  don't really
    // need a threadpool for this, but the infrastructure is there...
    if (sdata->sync_count < sdata->THREADS)
    {
        tdata->work_fcn_id = 0;
        sdata->sync_count++;
    }
    else
    {
        tdata->work_fcn_id = tdata->num_work_fcn;
    }

    return;
}

void count_twins_work_fcn(void* vptr)
{
    tpool_t* tdata = (tpool_t*)vptr;
    soe_userdata_t* udata = (soe_userdata_t*)tdata->user_data;
    soe_staticdata_t* sdata = udata->sdata;
    thread_soedata_t* t = &udata->ddata[tdata->tindex];
    int i;

    uint64_t count = 0;
    uint64_t numchunks = (t->stopid - t->startid);


    // counting twins and other prime constellations requires
    // all of the sieve lines in memory for reording.
#if defined(USE_BMI2) || defined(USE_AVX512F)
    if ((sdata->has_bmi2) && (sdata->numclasses <= 270) && (!(sdata->numclasses == 96))) 
        // && (!(sdata->analysis > 1)))
    {
        for (i = 0; i < numchunks; i++)
        {
            count = count_8_bytes_bmi2(sdata, count, (uint64_t)(t->startid + i) * 8);

            // if searching for prime constellations, then we need to look at the last 
            // flags of this block of 8 bytes and the first flags of the next one.
            // we do this by loading the trailing bits of this block into a carry
            // register.  The rest is handled by the block analysis function.
            if ((sdata->analysis == 2) && (sdata->is_main_sieve == 1))
            {
                uint8_t* lastline = sdata->lines[sdata->numclasses - 1];
                uint8_t* firstline = sdata->lines[0];

                if ((i + 1) < numchunks)
                {
                    uint8_t lastflag = lastline[(uint64_t)(t->startid + i) * 8 + 7] & 0x80;
                    uint8_t firstflag = firstline[(uint64_t)(t->startid + i) * 8 + 8] & 0x1;

                    if (lastflag && firstflag)
                    {
                        count++;
                    }
                }
                else if ((i + 1) < sdata->numlinebytes)
                {
                    // this thread is done but if there is more data after
                    // this thread's chunk then check between thread boundaries.
                    uint8_t lastflag = lastline[(uint64_t)(t->startid + i) * 8 + 7] & 0x80;
                    uint8_t firstflag = firstline[(uint64_t)(t->startid + i) * 8 + 8] & 0x1;

                    if (lastflag && firstflag)
                    {
                        count++;
                    }
                }
            }
        }
    }
    else
    {
        // if we don't have BMI2, or an unsupported-by-bmi2-numclasses, 
        // then we will end up here.
        for (i = 0; i < numchunks; i++)
        {
            count = count_8_bytes(sdata, count, (uint64_t)(t->startid + i) * 8);

            // if searching for prime constellations, then we need to look at the last 
            // flags of this block of 8 bytes and the first flags of the next one.
            // we do this by loading the trailing bits of this block into a carry
            // register.  The rest is handled by the block analysis function.
            if ((sdata->analysis == 2) && (sdata->is_main_sieve == 1))
            {
                uint8_t* lastline = sdata->lines[sdata->numclasses - 1];
                uint8_t* firstline = sdata->lines[0];

                if ((i + 1) < numchunks)
                {
                    uint8_t lastflag = lastline[(uint64_t)(t->startid + i) * 8 + 7] & 0x80;
                    uint8_t firstflag = firstline[(uint64_t)(t->startid + i) * 8 + 8] & 0x1;

                    if (lastflag && firstflag)
                    {
                        count++;
                    }
                }
                else if ((i + 1) < sdata->numlinebytes)
                {
                    // this thread is done but if there is more data after
                    // this thread's chunk then check between thread boundaries.
                    uint8_t lastflag = lastline[(uint64_t)(t->startid + i) * 8 + 7] & 0x80;
                    uint8_t firstflag = firstline[(uint64_t)(t->startid + i) * 8 + 8] & 0x1;

                    if (lastflag && firstflag)
                    {
                        count++;
                    }
                }
            }
        }
    }
#else
    for (i = 0; i < numchunks; i++)
    {
        count = count_8_bytes(sdata, count, (uint64_t)(t->startid + i) * 8);

        // if searching for prime constellations, then we need to look at the last 
        // flags of this block of 8 bytes and the first flags of the next one.
        // we do this by loading the trailing bits of this block into a carry
        // register.  The rest is handled by the block analysis function.
        if ((sdata->analysis == 2) && (sdata->is_main_sieve == 1))
        {
            uint8_t* lastline = sdata->lines[sdata->numclasses - 1];
            uint8_t* firstline = sdata->lines[0];

            if ((i + 1) < numchunks)
            {
                uint8_t lastflag = lastline[(uint64_t)(t->startid + i) * 8 + 7] & 0x80;
                uint8_t firstflag = firstline[(uint64_t)(t->startid + i) * 8 + 8] & 0x1;

                if (lastflag && firstflag)
                {
                    count++;
                }
            }
            else if ((t->startid + i + 1) < sdata->numlinebytes)
            {
                // this thread is done but if there is more data after
                // this thread's chunk then check between thread boundaries.
                uint8_t lastflag = lastline[(uint64_t)(t->startid + i) * 8 + 7] & 0x80;
                uint8_t firstflag = firstline[(uint64_t)(t->startid + i) * 8 + 8] & 0x1;

                if (lastflag && firstflag)
                {
                    count++;
                }
            }
        }
    }
#endif

    t->linecount = count;

    return;
}

uint64_t count_line(soe_staticdata_t *sdata, uint32_t current_line)
{
	// extract stuff from the thread data structure
	uint8_t *line = sdata->lines[current_line];
	uint64_t numlinebytes = sdata->numlinebytes;
	uint64_t lowlimit = sdata->lowlimit;
	uint64_t prodN = sdata->prodN;
	uint8_t *flagblock = line;
	uint64_t i, it = 0;
    uint64_t stopcount;
	int ix;
	int done, kx;
	uint64_t prime;
    uint8_t* masks = sdata->masks;
    uint8_t* nmasks = sdata->nmasks;

    if (sdata->sieve_range)
    {
        // if we are counting primes in a range with large offset 
        // then we need to compute the actual values of candidate primes
        // and run PRP tests on them.
        mpz_t tmpz;
        mpz_init(tmpz);
        uint64_t* flagblock64 = (uint64_t*)line;
        uint64_t next_report = 1000;
        int num_cand = 0;
        uint64_t numchunks = sdata->numlinebytes / 8;

#ifdef USE_AVX512F
        
        uint64_t candidates[64];   // space for worst-case 64 prime flags to test;
        int j;
        i = 0;

        do
        {
            while ((num_cand < 8) && (i < numchunks))
            {
                uint64_t x = flagblock64[i];

                while (x > 0)
                {
                    uint32_t idx = _trail_zcnt64(x);
                    uint64_t value = (i * 64 + idx) * prodN + sdata->rclass[current_line];
                    if ((value >= sdata->orig_llimit) && (value <= sdata->orig_hlimit))
                    {
                        candidates[num_cand++] = value;
                    }
                    x = _reset_lsb64(x);
                }

                if (num_cand >= 64)
                {
                    printf("WARNING: too many candidates\n\n");
                }

                // next chunk
                i++;
            }

            ALIGNED_MEM uint64_t n8[16];
            uint8_t loc_msk;
            uint8_t gmp_msk = 0;

            for (j = 0, loc_msk = 0; j < MIN(8, num_cand); j++)
            {
                mpz_add_ui(tmpz, sdata->offset, candidates[j]);

                loc_msk |= (1ull << j);

                n8[j] = mpz_get_ui(tmpz) & 0xfffffffffffffull;
                mpz_tdiv_q_2exp(tmpz, tmpz, 52);
                n8[j + 8] = mpz_get_ui(tmpz) & 0xfffffffffffffull;
            }

            if (num_cand > 0)
            {
                // fill any remaining spots with repeats of the last one:
                // not sure what MR_2sprp will do with garbage in a slot.
                for (; j < 8; j++)
                {
                    n8[j] = n8[0];
                    n8[j + 8] = n8[8];
                }

                uint8_t prpmask = loc_msk;
                switch (sdata->witnesses)
                {
                case 1: prpmask &= fermat_prp_104x8(n8); break;
                case 2: prpmask &= MR_2sprp_104x8(n8); break;
                default: prpmask &= MR_2sprp_104x8(n8); break;
                }

                for (j = 0; j < num_cand; j++)
                {
                    if (prpmask & (1ull << j))
                    {
                        it++;
                    }

                    candidates[j] = candidates[j + MIN(8, num_cand)];
                }
                num_cand -= MIN(8, num_cand);

                if (num_cand < 0)
                {
                    print("WARNING: num_cand < 0\n");
                }
            }

            if ((i > next_report) && (sdata->VFLAG > 0))
            {
                printf("PRP progress: %d%%\r",
                    (int)((double)(i) / (double)(numchunks) * 100.0));
                fflush(stdout);
                next_report = i + 1000;
            }

        } while (i < numchunks);


#else

        for (i = 0; i < numchunks; i++)
        {
            uint64_t x = flagblock64[i];
            int j;

            for (j = 0; j < 64; j++)
            {
                if (x & (1ull << j))
                {
                    uint64_t value = (i * 64 + j) * prodN + sdata->rclass[current_line];
                    mpz_add_ui(tmpz, sdata->offset, value);

                    if ((value >= sdata->orig_llimit) && (value <= sdata->orig_hlimit))
                    {
                        if (mpz_probab_prime_p(tmpz, 1))
                        {
                            it++;
                        }
                    }
                }
            }

            if ((i > next_report) && (sdata->VFLAG > 0))
            {
                printf("PRP progress: %d%%\r",
                    (int)((double)(i) / (double)(numchunks) * 100.0));
                fflush(stdout);
                next_report = i + 1000;
            }
        }


#endif

        mpz_clear(tmpz);
        
    }
    else
    {

#ifdef USE_AVX2

        __m256i v5, v3, v0f, v3f;
        uint32_t* tmp;

        v5 = _mm256_set1_epi32(0x55555555);
        v3 = _mm256_set1_epi32(0x33333333);
        v0f = _mm256_set1_epi32(0x0F0F0F0F);
        v3f = _mm256_set1_epi32(0x0000003F);
        tmp = (uint32_t*)xmalloc_align(8 * sizeof(uint32_t));

        uint64_t numchunks = (sdata->orig_hlimit - lowlimit) / (512 * prodN) + 1;

        stopcount = numchunks * 2; // i / 32;
        for (i = 0; i < stopcount; i += 2)
        {
            __m256i t1, t2, t3, t4;
            __m256i x = _mm256_load_si256((__m256i*)(&flagblock[32 * i]));
            __m256i y = _mm256_load_si256((__m256i*)(&flagblock[32 * i + 32]));
            t1 = _mm256_srli_epi64(x, 1);
            t3 = _mm256_srli_epi64(y, 1);
            t1 = _mm256_and_si256(t1, v5);
            t3 = _mm256_and_si256(t3, v5);
            x = _mm256_sub_epi64(x, t1);
            y = _mm256_sub_epi64(y, t3);
            t1 = _mm256_and_si256(x, v3);
            t3 = _mm256_and_si256(y, v3);
            t2 = _mm256_srli_epi64(x, 2);
            t4 = _mm256_srli_epi64(y, 2);
            t2 = _mm256_and_si256(t2, v3);
            t4 = _mm256_and_si256(t4, v3);
            x = _mm256_add_epi64(t2, t1);
            y = _mm256_add_epi64(t4, t3);
            t1 = _mm256_srli_epi64(x, 4);
            t3 = _mm256_srli_epi64(y, 4);
            x = _mm256_add_epi64(x, t1);
            y = _mm256_add_epi64(y, t3);
            x = _mm256_and_si256(x, v0f);
            y = _mm256_and_si256(y, v0f);
            t1 = _mm256_srli_epi64(x, 8);
            t3 = _mm256_srli_epi64(y, 8);
            x = _mm256_add_epi64(x, t1);
            y = _mm256_add_epi64(y, t3);
            t1 = _mm256_srli_epi64(x, 16);
            t3 = _mm256_srli_epi64(y, 16);
            x = _mm256_add_epi64(x, t1);
            y = _mm256_add_epi64(y, t3);
            t1 = _mm256_srli_epi64(x, 32);
            t3 = _mm256_srli_epi64(y, 32);
            x = _mm256_add_epi64(x, t1);
            y = _mm256_add_epi64(y, t3);
            x = _mm256_and_si256(x, v3f);
            y = _mm256_and_si256(y, v3f);
            _mm256_store_si256((__m256i*)tmp, x);
            it += tmp[0] + tmp[2] + tmp[4] + tmp[6];
            _mm256_store_si256((__m256i*)tmp, y);
            it += tmp[0] + tmp[2] + tmp[4] + tmp[6];

        }

        align_free(tmp);

#else

        // process 64 bits at a time by using Warren's algorithm
        uint64_t numchunks = (sdata->orig_hlimit - lowlimit) / (64 * prodN) + 1;
        uint64_t* flagblock64 = (uint64_t*)line;

        for (i = 0; i < numchunks; i++)
        {
            /* Convert to 64-bit unsigned integer */
            uint64_t x = flagblock64[i];

            /*  Employ bit population counter algorithm from Henry S. Warren's
                *  "Hacker's Delight" book, chapter 5.   Added one more shift-n-add
                *  to accomdate 64 bit values.
                */

            x = x - ((x >> 1) & 0x5555555555555555ULL);
            x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
            x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
            x = x + (x >> 8);
            x = x + (x >> 16);
            x = x + (x >> 32);

            it += (x & 0x000000000000003FULL);

        }

#endif

    }

	return it;
}

uint64_t count_twins(soe_staticdata_t* sdata, thread_soedata_t* thread_data)
{
    // compute twins using all of the sieved lines we have stored
    uint32_t pcount = 0;
    uint64_t i;
    int j;
    uint32_t range, lastid;

    //timing
    double t;
    struct timeval tstart, tstop;

    // threading structures
    tpool_t* tpool_data;
    soe_userdata_t udata;

    if (sdata->VFLAG > 1)
    {
        gettimeofday(&tstart, NULL);
    }

    // total number of 64-bit chunks needed
    uint32_t numchunks = (sdata->orig_hlimit - sdata->lowlimit) / (64 * sdata->prodN) + 1;
    range = numchunks / sdata->THREADS;
    lastid = 0;

    if (sdata->VFLAG > 1)
    {
        printf("counting twins from %lu to %lu\n", sdata->orig_llimit, sdata->orig_hlimit);
    }

    // divvy up the line bytes.  Unlike when counting primes,
    // here threading is by groups of bytes, over all classes.
    // this is necessary to order the primes.
    for (i = 0; i < sdata->THREADS; i++)
    {
        thread_soedata_t* t = thread_data + i;

        t->sdata = *sdata;
        t->startid = lastid;
        t->stopid = t->startid + range;

        if (i == (sdata->THREADS - 1))
        {
            t->stopid = numchunks;
        }
        lastid = t->stopid;

        if (sdata->VFLAG > 2)
        {
            printf("thread %d counting twins from byte offset %u to %u\n",
                (int)i, t->startid * 8, t->stopid * 8);
        }
    }

    udata.sdata = sdata;
    udata.ddata = thread_data;
    tpool_data = tpool_setup(sdata->THREADS, NULL, NULL, NULL,
        &count_twins_dispatch, &udata);

    if (sdata->THREADS == 1)
    {
        thread_data->linecount = 0;
        count_twins_work_fcn(tpool_data);
    }
    else
    {
        sdata->sync_count = 0;
        tpool_add_work_fcn(tpool_data, &count_twins_work_fcn);
        tpool_go(tpool_data);
    }
    free(tpool_data);

    // now combine all of the results
    if (sdata->THREADS > 1)
    {
        pcount = 0;
        for (j = 0; j < sdata->THREADS; j++)
        {
            thread_soedata_t* t = thread_data + j;

            if (t->linecount == 0)
            {
                continue;
            }

            if (sdata->VFLAG > 2)
            {
                printf("adding %" PRIu64 " twins found in thread %d\n", t->linecount, j);
            }

            pcount += t->linecount;
        }
    }
    else
    {
        pcount = thread_data[0].linecount;
    }

    if (sdata->VFLAG > 1)
    {
        gettimeofday(&tstop, NULL);

        t = ytools_difftime(&tstart, &tstop);

        if (sdata->VFLAG > 2)
        {
            printf("time to count twins = %1.4f\n", t);
        }
    }

    return pcount;
}

uint64_t count_twins_nonthreaded(soe_staticdata_t* sdata)
{
    int i;
    uint64_t lowlimit = sdata->lowlimit;
    uint64_t count = 0;
    uint64_t prodN = sdata->prodN;
    uint64_t numchunks = (sdata->orig_hlimit - lowlimit) / (64 * prodN) + 1;

	// counting twins and other prime constellations requires
	// all of the sieve lines in memory for reording.
#if defined(USE_BMI2) || defined(USE_AVX512F)
    if ((sdata->has_bmi2) && (sdata->numclasses <= 48))
    {
        for (i = 0; i < numchunks; i++)
        {
            count = count_8_bytes_bmi2(sdata, count, (uint64_t)i * 8);

            // if searching for prime constellations, then we need to look at the last 
            // flags of this block of 8 bytes and the first flags of the next one.
            // we do this by loading the trailing bits of this block into a carry
            // register.  The rest is handled by the block analysis function.
            // For the first block in a multi-threaded run, we should also
            // be receiving carry data from the previous thread in order to
            // maintain continuity across threads.

            if ((sdata->analysis == 2) && (sdata->is_main_sieve == 1))
            {
                uint8_t* lastline = sdata->lines[sdata->numclasses - 1];
                uint8_t* firstline = sdata->lines[0];

                if ((i + 1) < numchunks)
                {
                    uint8_t lastflag = lastline[i * 8 + 7] & 0x80;
                    uint8_t firstflag = firstline[i * 8 + 8] & 0x1;

                    if (lastflag && firstflag)
                    {
                        count++;
                    }
                }
                else if ((i + 1) < sdata->numlinebytes)
                {
                    // this thread is done but if there is more data after
                    // this thread's chunk then check between thread boundaries.
                    uint8_t lastflag = lastline[i * 8 + 7] & 0x80;
                    uint8_t firstflag = firstline[i * 8 + 8] & 0x1;

                    if (lastflag && firstflag)
                    {
                        count++;
                    }
                }
            }
        }
    }
    else
    {
        // if we don't have BMI2, or numclasses > 48, then we should end up here.
        for (i = 0; i < numchunks; i++)
        {
            count = count_8_bytes(sdata, count, (uint64_t)i * 8);

            // if searching for prime constellations, then we need to look at the last 
            // flags of this block of 8 bytes and the first flags of the next one.
            // we do this by loading the trailing bits of this block into a carry
            // register.  The rest is handled by the block analysis function.
            // For the first block in a multi-threaded run, we should also
            // be receiving carry data from the previous thread in order to
            // maintain continuity across threads.

            if ((sdata->analysis == 2) && (sdata->is_main_sieve == 1))
            {
                uint8_t* lastline = sdata->lines[sdata->numclasses - 1];
                uint8_t* firstline = sdata->lines[0];

                if ((i + 1) < numchunks)
                {
                    uint8_t lastflag = lastline[i * 8 + 7] & 0x80;
                    uint8_t firstflag = firstline[i * 8 + 8] & 0x1;

                    if (lastflag && firstflag)
                    {
                        count++;
                    }
                }
                else if ((i + 1) < sdata->numlinebytes)
                {
                    // this thread is done but if there is more data after
                    // this thread's chunk then check between thread boundaries.
                    uint8_t lastflag = lastline[i * 8 + 7] & 0x80;
                    uint8_t firstflag = firstline[i * 8 + 8] & 0x1;

                    if (lastflag && firstflag)
                    {
                        count++;
                    }
                }
            }
        }
    }
#else
    for (i = 0; i < numchunks; i++)
    {
        count = count_8_bytes(sdata, count, (uint64_t)i * 8);

        // if searching for prime constellations, then we need to look at the last 
        // flags of this block of 8 bytes and the first flags of the next one.
        // we do this by loading the trailing bits of this block into a carry
        // register.  The rest is handled by the block analysis function.
        // For the first block in a multi-threaded run, we should also
        // be receiving carry data from the previous thread in order to
        // maintain continuity across threads.

        if ((sdata->analysis == 2) && (sdata->is_main_sieve == 1))
        {
            uint8_t* lastline = sdata->lines[sdata->numclasses - 1];
            uint8_t* firstline = sdata->lines[0];

            if ((i + 1) < numchunks)
            {
                uint8_t lastflag = lastline[i * 8 + 7] & 0x80;
                uint8_t firstflag = firstline[i * 8 + 8] & 0x1;

                if (lastflag && firstflag)
                {
                    count++;
                }
            }
            else if ((i + 1) < sdata->numlinebytes)
            {
                // this thread is done but if there is more data after
                // this thread's chunk then check between thread boundaries.
                uint8_t lastflag = lastline[i * 8 + 7] & 0x80;
                uint8_t firstflag = firstline[i * 8 + 8] & 0x1;

                if (lastflag && firstflag)
                {
                    count++;
                }
            }
        }
    }
#endif

	return count;
}

void count_line_special(thread_soedata_t *thread_data)
{
	//extract stuff from the thread data structure
	soe_staticdata_t *sdata = &thread_data->sdata;
	uint32_t current_line = thread_data->current_line;
	uint8_t *line = sdata->lines[current_line];	
	uint64_t numlinebytes = sdata->numlinebytes;
	uint64_t lowlimit = sdata->lowlimit;
	uint64_t prodN = sdata->prodN;
	uint64_t *flagblock64 = (uint64_t *)line;
	uint64_t i, k, it, lower, upper;
	int ix;
	int64_t start, stop;
    uint8_t* masks = sdata->masks;
    uint8_t* nmasks = sdata->nmasks;

	//zero out any bits below the requested range
	for (i=lowlimit + sdata->rclass[current_line], ix=0; i < sdata->orig_llimit; i += prodN, ix++)
		line[ix >> 3] &= masks[ix & 7];
	
	//and any high bits above the requested range
	for (i=sdata->highlimit + sdata->rclass[current_line] - prodN, ix=0; i > sdata->orig_hlimit; i -= prodN, ix++)
		line[numlinebytes - 1 - (ix >> 3)] &= masks[7 - (ix & 7)];

	//count each block of 1e9
	lower = sdata->orig_llimit;
	upper = lower;
	k = 0;
	thread_data->linecount = 0;
	while (upper != sdata->orig_hlimit)
	{
		//set the bounds for the next batch
		upper = upper + 1000000000; 
		if (upper > sdata->orig_hlimit)
			upper = sdata->orig_hlimit;

		//find the starting byte number.  first find the number of bits between the current lower
		//limit and the start of the line.
		start = (int64_t)((lower - lowlimit) / prodN);

		//we'll be counting in 64 bit chunks, so compute how many 64 bit chunks this is
		start /= 64;
		
		//start a little before the range, to account for rounding errors
		start -= 2;

		if (start < 0) start = 0;

		//find the stopping byte number: first find the number of bits between the current upper
		//limit and the start of the line.
		stop = (int64_t)((upper - lowlimit) / prodN);

		//we'll be counting in 64 bit chunks, so compute how many 64 bit chunks this is
		stop /= 64;
		
		//stop a little after the range, to account for rounding errors
		stop += 2;

		if (stop > (numlinebytes >> 3)) stop = (numlinebytes >> 3);

		//count these bytes
		it = 0;
		for (ix = start; ix < stop; ix++)
		{
			/* Convert to 64-bit unsigned integer */    
			uint64_t x = flagblock64[ix];
		    
			/*  Employ bit population counter algorithm from Henry S. Warren's
			 *  "Hacker's Delight" book, chapter 5.   Added one more shift-n-add
			 *  to accomdate 64 bit values.
			 */
			x = x - ((x >> 1) & 0x5555555555555555ULL);
			x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
			x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
			x = x + (x >> 8);
			x = x + (x >> 16);
			x = x + (x >> 32);

			it += (x & 0x000000000000003FULL);
		}

		//then correct the counts
		//zero out any bits below the requested range
		for (i= (start * 64) * prodN + sdata->rclass[current_line] + lowlimit, ix=0; i < lower; i += prodN, ix++)
		{
			if (line[(ix >> 3) + (start << 3)] & nmasks[ix & 7])
				it--;
		}
		
		//and any high bits above the requested range
		for (i=(stop * 64) * prodN + sdata->rclass[current_line] - prodN + lowlimit, ix=0; i > upper; i -= prodN, ix++)
		{
			if (line[(stop << 3) - 1 - (ix >> 3)] & nmasks[7 - (ix & 7)])
				it--;

		}
		
		//add the count to the special array
		thread_data->ddata.special_count[k] = it;
		thread_data->linecount += it;
		k++;
		lower = upper;
	}

}

uint64_t count_8_bytes(soe_staticdata_t* sdata,
    uint64_t pcount, uint64_t byte_offset)
{
    uint32_t current_line;
    // re-ordering queues supporting arbitrary residue classes.
    uint64_t** pqueues; // [64][48]
    uint32_t pcounts[64];
    int i, j;
    uint32_t nc = sdata->numclasses;
    uint64_t lowlimit = sdata->lowlimit;
    uint64_t prodN = sdata->prodN;
    uint8_t** lines = sdata->lines;
    uint64_t olow = sdata->orig_llimit;
    uint64_t ohigh = sdata->orig_hlimit;
    int GLOBAL_OFFSET = sdata->GLOBAL_OFFSET;

    if ((byte_offset & 32767) == 0)
    {
        if (sdata->VFLAG > 1)
        {
            printf("computing: %d%%\r", (int)
                ((double)byte_offset / (double)(sdata->numlinebytes) * 100.0));
            fflush(stdout);
        }
    }

    pqueues = (uint64_t**)xmalloc(64 * sizeof(uint64_t*));
    for (i = 0; i < 64; i++)
    {
        pqueues[i] = (uint64_t*)xmalloc(sdata->numclasses * sizeof(uint64_t));
    }

    // Compute the primes using ctz on the 64-bit words but push the results
    // into 64 different queues depending on the bit position.  Then
    // we pull from the queues in order while storing into the primes array.
    // This time the bottleneck is mostly in the queue-based sorting
    // and associated memory operations, so we don't bother with
    // switching between branch-free inner loops or not.
    memset(pcounts, 0, 64 * sizeof(uint32_t));

    lowlimit += byte_offset * 8 * prodN;
    for (current_line = 0; current_line < nc; current_line++)
    {
        uint64_t* line64 = (uint64_t*)lines[current_line];
        uint64_t flags64 = line64[byte_offset / 8];

        while (flags64 > 0)
        {
            uint64_t pos = _trail_zcnt64(flags64);
            uint64_t prime = lowlimit + pos * prodN + sdata->rclass[current_line];

            if ((prime >= olow) && (prime <= ohigh))
            {
                pqueues[pos][pcounts[pos]] = prime;
                pcounts[pos]++;
            }
            flags64 ^= (1ULL << pos);
        }
    }

    for (i = 0; i < 64; i++)
    {

        // search for twins and only load the leading element/prime.
        // if depth-based sieving then these are candidate twins.
        if (pcounts[i] > 0)
        {

            for (j = 0; j < pcounts[i] - 1; j++)
            {
                if ((pqueues[i][j + 1] - pqueues[i][j]) == 2)
                {
                    pcount++;
                }
            }
            if (i < 63)
            {
                if (pcounts[i + 1] > 0)
                {
                    if ((pqueues[i + 1][0] - pqueues[i][j]) == 2)
                    {
                        pcount++;
                    }
                }
            }
        }

    }

    for (i = 0; i < 64; i++)
    {
        free(pqueues[i]);
    }
    free(pqueues);

    return pcount;
}



#if defined(USE_BMI2) || defined(USE_AVX512F)

__inline uint64_t interleave_pdep2x32(uint32_t x1, uint32_t x2)
{
    return _pdep_u64(x1, 0x5555555555555555)
        | _pdep_u64(x2, 0xaaaaaaaaaaaaaaaa);
}

__inline uint64_t interleave_pdep_8x8(uint8_t x1,
    uint8_t x2,
    uint8_t x3,
    uint8_t x4,
    uint8_t x5,
    uint8_t x6,
    uint8_t x7,
    uint8_t x8)
{
    return _pdep_u64(x1, 0x0101010101010101ull) |
        _pdep_u64(x2, 0x0202020202020202ull) |
        _pdep_u64(x3, 0x0404040404040404ull) |
        _pdep_u64(x4, 0x0808080808080808ull) |
        _pdep_u64(x5, 0x1010101010101010ull) |
        _pdep_u64(x6, 0x2020202020202020ull) |
        _pdep_u64(x7, 0x4040404040404040ull) |
        _pdep_u64(x8, 0x8080808080808080ull);
}

#define BIT0 0x1
#define BIT1 0x2
#define BIT2 0x4
#define BIT3 0x8
#define BIT4 0x10
#define BIT5 0x20
#define BIT6 0x40
#define BIT7 0x80

#define BYTE_TO_BINARY_PATTERN "%c%c%c%c%c%c%c%c"
#define BYTE_TO_BINARY(byte)  \
  ((byte) & 0x01 ? '1' : '0'), \
  ((byte) & 0x02 ? '1' : '0'), \
  ((byte) & 0x04 ? '1' : '0'), \
  ((byte) & 0x08 ? '1' : '0'), \
  ((byte) & 0x10 ? '1' : '0'), \
  ((byte) & 0x20 ? '1' : '0'), \
  ((byte) & 0x40 ? '1' : '0'), \
  ((byte) & 0x80 ? '1' : '0') 

uint64_t count_8_bytes_bmi2(soe_staticdata_t* sdata,
    uint64_t pcount, uint64_t byte_offset)
{
    uint32_t nc = sdata->numclasses;
    uint64_t lowlimit = sdata->lowlimit;
    uint8_t** lines = sdata->lines;

    if ((byte_offset & 32767) == 0)
    {
        if (sdata->VFLAG > 1)
        {
            printf("computing: %d%%\r", (int)
                ((double)byte_offset / (double)(sdata->numlinebytes) * 100.0));
            fflush(stdout);
        }
    }

    // AVX2 version, new instructions help quite a bit:
    // use _pdep_u64 to align/interleave bits from multiple bytes, 
    // _blsr_u64 to clear the last set bit, and depending on the 
    // number of residue classes, AVX2 vector load/store operations.

    // here is the 2 line version
    if (nc == 2)
    {
        int i;
        uint32_t last_bit = 0;
        uint32_t* lines32a = (uint32_t*)lines[0];
        uint32_t* lines32b = (uint32_t*)lines[1];
        
        // align the current bytes in next 2 residue classes
        for (i = 0; i < 2; i++)
        {
            uint64_t aligned_flags;

            aligned_flags = interleave_pdep2x32(
                lines32a[byte_offset / 4 + i],
                lines32b[byte_offset / 4 + i]);


            // alternate bits encode potential primes
            // in residue classes 1 and 5.  So twins can 
            // only exist with flags in class 5 followed by 1.
            uint64_t twins = aligned_flags & (aligned_flags >> 1);

            twins &= 0xaaaaaaaaaaaaaaaaULL;

            pcount += _mm_popcnt_u64(twins);
            pcount += (last_bit & aligned_flags);

            last_bit = (aligned_flags >> 63);
        }
    }
    else if (nc == 8)
    {
        int i;
        uint32_t last_bit = 0;

        // align the current bytes in next 8 residue classes
        for (i = 0; i < 8; i++)
        {
            uint64_t aligned_flags;

            aligned_flags = interleave_pdep_8x8(lines[0][byte_offset + i],
                lines[1][byte_offset + i],
                lines[2][byte_offset + i],
                lines[3][byte_offset + i],
                lines[4][byte_offset + i],
                lines[5][byte_offset + i],
                lines[6][byte_offset + i],
                lines[7][byte_offset + i]);

            uint64_t twins = aligned_flags & (aligned_flags >> 1);

            twins &= 0x9494949494949494ULL;

            pcount += _mm_popcnt_u64(twins);
            pcount += (last_bit & aligned_flags);
            
            last_bit = (aligned_flags >> 63);
        }
    }
    else if (nc == 30)
    {
        int i;
        uint32_t last_bit = 0;

        for (i = 0; i < 8; i++)
        {
            uint64_t aligned_flags1;
            uint64_t aligned_flags2;
            uint64_t aligned_flags3;
            uint64_t aligned_flags4;

            // first we partially order the lines such
            // that each 64-bit flags contains 8 ordered
            // bytes for a set of 8 classes.
            aligned_flags1 = interleave_pdep_8x8(
                lines[0][byte_offset + i],
                lines[1][byte_offset + i],
                lines[2][byte_offset + i],
                lines[3][byte_offset + i],
                lines[4][byte_offset + i],
                lines[5][byte_offset + i],
                lines[6][byte_offset + i],
                lines[7][byte_offset + i]);

            aligned_flags2 = interleave_pdep_8x8(
                lines[8 + 0][byte_offset + i],
                lines[8 + 1][byte_offset + i],
                lines[8 + 2][byte_offset + i],
                lines[8 + 3][byte_offset + i],
                lines[8 + 4][byte_offset + i],
                lines[8 + 5][byte_offset + i],
                lines[8 + 6][byte_offset + i],
                lines[8 + 7][byte_offset + i]);

            aligned_flags3 = interleave_pdep_8x8(
                lines[16 + 0][byte_offset + i],
                lines[16 + 1][byte_offset + i],
                lines[16 + 2][byte_offset + i],
                lines[16 + 3][byte_offset + i],
                lines[16 + 4][byte_offset + i],
                lines[16 + 5][byte_offset + i],
                lines[16 + 6][byte_offset + i],
                lines[16 + 7][byte_offset + i]);

            aligned_flags4 = interleave_pdep_8x8(
                lines[24 + 0][byte_offset + i],
                lines[24 + 1][byte_offset + i],
                lines[24 + 2][byte_offset + i],
                lines[24 + 3][byte_offset + i],
                lines[24 + 4][byte_offset + i],
                lines[24 + 5][byte_offset + i],
                0 , 
                0);


            // aligned_flags1 contains: b0(c0-7), b1(c0-7), b2(c0-7), ... b7(c0-7)
            // aligned_flags2 contains: b0(c8-15), b1(c8-15), b2(c8-15), ... b7(c8-15)
            // aligned_flags3 contains: b0(c16-23), b1(c16-23), b2(c16-23), ... b7(c16-23)
            // aligned_flags4 contains: b0(c24-29), b1(c24-29), b2(c24-29), ... b7(c24-29)

            // now, shuffle the bytes within the partially ordered chunks.
            uint8_t unordered_bytes[48];
            uint8_t ordered_bytes[64];
            uint64_t* unordered64 = (uint64_t*)unordered_bytes;
            uint64_t* ordered64 = (uint64_t*)ordered_bytes;
            uint32_t* ordered32 = (uint32_t*)ordered_bytes;
            unordered64[0] = aligned_flags1;
            unordered64[1] = aligned_flags2;
            unordered64[2] = aligned_flags3;
            unordered64[3] = aligned_flags4;

            ordered_bytes[0] = unordered_bytes[0];
            ordered_bytes[1] = unordered_bytes[8];
            ordered_bytes[2] = unordered_bytes[16];
            ordered_bytes[3] = unordered_bytes[24];
            ordered_bytes[4] = unordered_bytes[1];
            ordered_bytes[5] = unordered_bytes[9];
            ordered_bytes[6] = unordered_bytes[17];
            ordered_bytes[7] = unordered_bytes[25];

            ordered_bytes[8] = unordered_bytes[2];
            ordered_bytes[9] = unordered_bytes[10];
            ordered_bytes[10] = unordered_bytes[18];
            ordered_bytes[11] = unordered_bytes[26];
            ordered_bytes[12] = unordered_bytes[3];
            ordered_bytes[13] = unordered_bytes[11];
            ordered_bytes[14] = unordered_bytes[19];
            ordered_bytes[15] = unordered_bytes[27];

            ordered_bytes[16] = unordered_bytes[4];
            ordered_bytes[17] = unordered_bytes[12];
            ordered_bytes[18] = unordered_bytes[20];
            ordered_bytes[19] = unordered_bytes[28];
            ordered_bytes[20] = unordered_bytes[5];
            ordered_bytes[21] = unordered_bytes[13];
            ordered_bytes[22] = unordered_bytes[21];
            ordered_bytes[23] = unordered_bytes[29];

            ordered_bytes[24] = unordered_bytes[6];
            ordered_bytes[25] = unordered_bytes[14];
            ordered_bytes[26] = unordered_bytes[22];
            ordered_bytes[27] = unordered_bytes[30];
            ordered_bytes[28] = unordered_bytes[7];
            ordered_bytes[29] = unordered_bytes[15];
            ordered_bytes[30] = unordered_bytes[23];
            ordered_bytes[31] = unordered_bytes[31];

            // rotate 2 bits for continuity between
            // the two 30-bit sections.
            ordered32[0] |= ((ordered32[1] & 0x3) << 30);
            ordered32[2] |= ((ordered32[3] & 0x3) << 30);
            ordered32[4] |= ((ordered32[5] & 0x3) << 30);
            ordered32[6] |= ((ordered32[7] & 0x3) << 30);
            ordered32[1] >>= 2;
            ordered32[3] >>= 2;
            ordered32[5] >>= 2;
            ordered32[7] >>= 2;

            aligned_flags1 = ordered64[0];
            aligned_flags2 = ordered64[1];
            aligned_flags3 = ordered64[2];
            aligned_flags4 = ordered64[3];

            // twin residues mod 210:
            // 1 11 13 17 19 29 31 41 43 59 61 71 73 101 103 107 109 137 139 149 151 167 169 179 181 191 193 197 199 209
            // binary mask: 0101 0101 0101 0101 0101 0101 0101 01|01 0101 0101 0101 0101 0101 0101 0100
            // hex mask: AAAAAAAAAAAAAA20

            uint64_t twins = aligned_flags1 & (aligned_flags1 >> 1);
            twins &= 0x02AAAAAAAAAAAAAAull;
            pcount += _mm_popcnt_u64(twins);
            pcount += (last_bit & aligned_flags1);  // from previous chunk
            last_bit = (aligned_flags1 >> 59);      // generate carry

            twins = aligned_flags2 & (aligned_flags2 >> 1);
            twins &= 0x02AAAAAAAAAAAAAAull;
            pcount += _mm_popcnt_u64(twins);        
            pcount += (last_bit & aligned_flags2);  // from previous chunk
            last_bit = (aligned_flags2 >> 59);      // generate carry

            twins = aligned_flags3 & (aligned_flags3 >> 1);
            twins &= 0x02AAAAAAAAAAAAAAull;
            pcount += _mm_popcnt_u64(twins);
            pcount += (last_bit & aligned_flags3);  // from previous chunk
            last_bit = (aligned_flags3 >> 59);      // generate carry

            twins = aligned_flags4 & (aligned_flags4 >> 1);
            twins &= 0x02AAAAAAAAAAAAAAull;
            pcount += _mm_popcnt_u64(twins);
            pcount += (last_bit & aligned_flags4);  // from previous chunk
            last_bit = (aligned_flags4 >> 59);      // generate carry


        }
    }
    else if (nc == 48)
    {
        int i;
        uint32_t last_bit = 0;

        for (i = 0; i < 8; i++)
        {
            uint64_t aligned_flags1;
            uint64_t aligned_flags2;
            uint64_t aligned_flags3;
            uint64_t aligned_flags4;
            uint64_t aligned_flags5;
            uint64_t aligned_flags6;

            // first we partially order the lines such
            // that each 64-bit flags contains 8 ordered
            // bytes for a set of 8 classes.
            aligned_flags1 = interleave_pdep_8x8(
                lines[0][byte_offset + i],
                lines[1][byte_offset + i],
                lines[2][byte_offset + i],
                lines[3][byte_offset + i],
                lines[4][byte_offset + i],
                lines[5][byte_offset + i],
                lines[6][byte_offset + i],
                lines[7][byte_offset + i]);

            aligned_flags2 = interleave_pdep_8x8(
                lines[8 + 0][byte_offset + i],
                lines[8 + 1][byte_offset + i],
                lines[8 + 2][byte_offset + i],
                lines[8 + 3][byte_offset + i],
                lines[8 + 4][byte_offset + i],
                lines[8 + 5][byte_offset + i],
                lines[8 + 6][byte_offset + i],
                lines[8 + 7][byte_offset + i]);

            aligned_flags3 = interleave_pdep_8x8(
                lines[16 + 0][byte_offset + i],
                lines[16 + 1][byte_offset + i],
                lines[16 + 2][byte_offset + i],
                lines[16 + 3][byte_offset + i],
                lines[16 + 4][byte_offset + i],
                lines[16 + 5][byte_offset + i],
                lines[16 + 6][byte_offset + i],
                lines[16 + 7][byte_offset + i]);

            aligned_flags4 = interleave_pdep_8x8(
                lines[24 + 0][byte_offset + i],
                lines[24 + 1][byte_offset + i],
                lines[24 + 2][byte_offset + i],
                lines[24 + 3][byte_offset + i],
                lines[24 + 4][byte_offset + i],
                lines[24 + 5][byte_offset + i],
                lines[24 + 6][byte_offset + i],
                lines[24 + 7][byte_offset + i]);

            aligned_flags5 = interleave_pdep_8x8(
                lines[32 + 0][byte_offset + i],
                lines[32 + 1][byte_offset + i],
                lines[32 + 2][byte_offset + i],
                lines[32 + 3][byte_offset + i],
                lines[32 + 4][byte_offset + i],
                lines[32 + 5][byte_offset + i],
                lines[32 + 6][byte_offset + i],
                lines[32 + 7][byte_offset + i]);

            aligned_flags6 = interleave_pdep_8x8(
                lines[40 + 0][byte_offset + i],
                lines[40 + 1][byte_offset + i],
                lines[40 + 2][byte_offset + i],
                lines[40 + 3][byte_offset + i],
                lines[40 + 4][byte_offset + i],
                lines[40 + 5][byte_offset + i],
                lines[40 + 6][byte_offset + i],
                lines[40 + 7][byte_offset + i]);

            // aligned_flags1 contains: b0(c0-7), b1(c0-7), b2(c0-7), ... b7(c0-7)
            // aligned_flags2 contains: b0(c8-15), b1(c8-15), b2(c8-15), ... b7(c8-15)
            // aligned_flags3 contains: b0(c16-23), b1(c16-23), b2(c16-23), ... b7(c16-23)
            // aligned_flags4 contains: b0(c24-31), b1(c24-31), b2(c24-31), ... b7(c24-31)
            // aligned_flags5 contains: b0(c32-39), b1(c32-39), b2(c32-39), ... b7(c32-39)
            // aligned_flags6 contains: b0(c40-47), b1(c40-47), b2(c40-47), ... b7(c40-47)
             
            // now, shuffle the bytes within the partially ordered chunks.
            uint8_t unordered_bytes[48];
            uint8_t ordered_bytes[64];
            uint64_t* unordered64 = (uint64_t*)unordered_bytes;
            uint64_t* ordered64 = (uint64_t*)ordered_bytes;
            unordered64[0] = aligned_flags1;
            unordered64[1] = aligned_flags2;
            unordered64[2] = aligned_flags3;
            unordered64[3] = aligned_flags4;
            unordered64[4] = aligned_flags5;
            unordered64[5] = aligned_flags6;

#if 1
            // minimum masks steps:
            ordered_bytes[0] = unordered_bytes[0];
            ordered_bytes[1] = unordered_bytes[8];
            ordered_bytes[2] = unordered_bytes[16];
            ordered_bytes[3] = unordered_bytes[24];
            ordered_bytes[4] = unordered_bytes[32];
            ordered_bytes[5] = unordered_bytes[40];
            ordered_bytes[6] = unordered_bytes[1];
            ordered_bytes[7] = unordered_bytes[9];

            ordered_bytes[8] = unordered_bytes[17];
            ordered_bytes[9] = unordered_bytes[25];
            ordered_bytes[10] = unordered_bytes[33];
            ordered_bytes[11] = unordered_bytes[41];
            ordered_bytes[12] = unordered_bytes[02];
            ordered_bytes[13] = unordered_bytes[10];
            ordered_bytes[14] = unordered_bytes[18];
            ordered_bytes[15] = unordered_bytes[26];

            ordered_bytes[16] = unordered_bytes[34];
            ordered_bytes[17] = unordered_bytes[42];
            ordered_bytes[18] = unordered_bytes[3];
            ordered_bytes[19] = unordered_bytes[11];
            ordered_bytes[20] = unordered_bytes[19];
            ordered_bytes[21] = unordered_bytes[27];
            ordered_bytes[22] = unordered_bytes[35];
            ordered_bytes[23] = unordered_bytes[43];

            ordered_bytes[24] = unordered_bytes[4];
            ordered_bytes[25] = unordered_bytes[12];
            ordered_bytes[26] = unordered_bytes[20];
            ordered_bytes[27] = unordered_bytes[28];
            ordered_bytes[28] = unordered_bytes[36];
            ordered_bytes[29] = unordered_bytes[44];
            ordered_bytes[30] = unordered_bytes[5];
            ordered_bytes[31] = unordered_bytes[13];

            ordered_bytes[32] = unordered_bytes[21];
            ordered_bytes[33] = unordered_bytes[29];
            ordered_bytes[34] = unordered_bytes[37];
            ordered_bytes[35] = unordered_bytes[45];
            ordered_bytes[36] = unordered_bytes[6];
            ordered_bytes[37] = unordered_bytes[14];
            ordered_bytes[38] = unordered_bytes[22];
            ordered_bytes[39] = unordered_bytes[30];

            ordered_bytes[40] = unordered_bytes[38];
            ordered_bytes[41] = unordered_bytes[46];
            ordered_bytes[42] = unordered_bytes[7];
            ordered_bytes[43] = unordered_bytes[15];
            ordered_bytes[44] = unordered_bytes[23];
            ordered_bytes[45] = unordered_bytes[31];
            ordered_bytes[46] = unordered_bytes[39];
            ordered_bytes[47] = unordered_bytes[47];

            aligned_flags1 = ordered64[0];
            aligned_flags2 = ordered64[1];
            aligned_flags3 = ordered64[2];
            aligned_flags4 = ordered64[3];
            aligned_flags5 = ordered64[4];
            aligned_flags6 = ordered64[5];



            // twin-mask for 48 classes.
            // residue classes mod 210:
            // 1                     0
            // 11 13 *               1 0
            // 17 19 *               1 0
            // 23                    0
            // 29 31 *               1 0
            // 37                    0
            // 41 43 *               1 0
            // 47 53                 0 0
            // 59 61 *               1 0
            // 67                    0          -- 224A
            // 71 73 *               1 0
            // 79 83 89 97           0 0 0 0
            // 101 103 *             1 0
            // 107 109 *             1 0
            // 113 121 127 131       0 0 0 0
            // 137 139 *             1 0        -- 4141
            // 143                   0
            // 149 151 *             1 0
            // 157 163               0 0 
            // 167 169 *             1 0 
            // 173                   0
            // 179 181 *             1 0
            // 187                   0
            // 191 193 *             1 0
            // 197 199 *             1 0
            // 209 *                 1          -- A922

            // 48-bit mask = 0xA9224141224A;
            // repeat that 6 times and split into 64-bit chunks
            // to go with the ordered flags:

            uint64_t twins = aligned_flags1 & (aligned_flags1 >> 1);
            twins &= 0x224AA9224141224Aull;
            pcount += _mm_popcnt_u64(twins);
            pcount += (last_bit & aligned_flags1);  // from previous chunk
                                                    // no carry generated

            twins = aligned_flags2 & (aligned_flags2 >> 1);
            twins &= 0x4141224AA9224141ull;
            pcount += _mm_popcnt_u64(twins);        // no carry generated

            twins = aligned_flags3 & (aligned_flags3 >> 1);
            twins &= 0xA9224141224AA922ull;
            pcount += _mm_popcnt_u64(twins);
            last_bit = (aligned_flags3 >> 63);      // generate carry

            twins = aligned_flags4 & (aligned_flags4 >> 1);
            twins &= 0x224AA9224141224Aull;
            pcount += _mm_popcnt_u64(twins);
            pcount += (last_bit & aligned_flags4);  // from previous chunk
                                                    // no carry generated

            twins = aligned_flags5 & (aligned_flags5 >> 1);
            twins &= 0x4141224AA9224141ull;
            pcount += _mm_popcnt_u64(twins);        // no carry generated

            twins = aligned_flags6 & (aligned_flags6 >> 1);
            twins &= 0xA9224141224AA922ull;
            pcount += _mm_popcnt_u64(twins);
            last_bit = (aligned_flags6 >> 63);      // generate carry

#else

            uint64_t aligned_flags7;
            uint64_t aligned_flags8;

            ordered_bytes[0] = unordered_bytes[0];
            ordered_bytes[1] = unordered_bytes[8];
            ordered_bytes[2] = unordered_bytes[16];
            ordered_bytes[3] = unordered_bytes[24];
            ordered_bytes[4] = unordered_bytes[32];
            ordered_bytes[5] = unordered_bytes[40];
            ordered_bytes[6] = 0;
            ordered_bytes[7] = 0;

            ordered_bytes[8] = unordered_bytes[1];
            ordered_bytes[9] = unordered_bytes[9];
            ordered_bytes[10] = unordered_bytes[17];
            ordered_bytes[11] = unordered_bytes[25];
            ordered_bytes[12] = unordered_bytes[33];
            ordered_bytes[13] = unordered_bytes[41];
            ordered_bytes[14] = 0;
            ordered_bytes[15] = 0;

            ordered_bytes[16] = unordered_bytes[2];
            ordered_bytes[17] = unordered_bytes[10];
            ordered_bytes[18] = unordered_bytes[18];
            ordered_bytes[19] = unordered_bytes[26];
            ordered_bytes[20] = unordered_bytes[34];
            ordered_bytes[21] = unordered_bytes[42];
            ordered_bytes[22] = 0;
            ordered_bytes[23] = 0;

            ordered_bytes[24] = unordered_bytes[3];
            ordered_bytes[25] = unordered_bytes[11];
            ordered_bytes[26] = unordered_bytes[19];
            ordered_bytes[27] = unordered_bytes[27];
            ordered_bytes[28] = unordered_bytes[35];
            ordered_bytes[29] = unordered_bytes[43];
            ordered_bytes[30] = 0;
            ordered_bytes[31] = 0;

            ordered_bytes[32] = unordered_bytes[4];
            ordered_bytes[33] = unordered_bytes[12];
            ordered_bytes[34] = unordered_bytes[20];
            ordered_bytes[35] = unordered_bytes[28];
            ordered_bytes[36] = unordered_bytes[36];
            ordered_bytes[37] = unordered_bytes[44];
            ordered_bytes[38] = 0;
            ordered_bytes[39] = 0;

            ordered_bytes[40] = unordered_bytes[5];
            ordered_bytes[41] = unordered_bytes[13];
            ordered_bytes[42] = unordered_bytes[21];
            ordered_bytes[43] = unordered_bytes[29];
            ordered_bytes[44] = unordered_bytes[37];
            ordered_bytes[45] = unordered_bytes[45];
            ordered_bytes[46] = 0;
            ordered_bytes[47] = 0;

            ordered_bytes[48] = unordered_bytes[6];
            ordered_bytes[49] = unordered_bytes[14];
            ordered_bytes[50] = unordered_bytes[22];
            ordered_bytes[51] = unordered_bytes[30];
            ordered_bytes[52] = unordered_bytes[38];
            ordered_bytes[53] = unordered_bytes[46];
            ordered_bytes[54] = 0;
            ordered_bytes[55] = 0;

            ordered_bytes[56] = unordered_bytes[7];
            ordered_bytes[57] = unordered_bytes[15];
            ordered_bytes[58] = unordered_bytes[23];
            ordered_bytes[59] = unordered_bytes[31];
            ordered_bytes[60] = unordered_bytes[39];
            ordered_bytes[61] = unordered_bytes[47];
            ordered_bytes[62] = 0;
            ordered_bytes[63] = 0;

            aligned_flags1 = ordered64[0];
            aligned_flags2 = ordered64[1];
            aligned_flags3 = ordered64[2];
            aligned_flags4 = ordered64[3];
            aligned_flags5 = ordered64[4];
            aligned_flags6 = ordered64[5];
            aligned_flags7 = ordered64[6];
            aligned_flags8 = ordered64[7];

            // twin-mask for 48 classes.
            // residue classes mod 210:
            // 1                     0
            // 11 13 *               1 0
            // 17 19 *               1 0
            // 23                    0
            // 29 31 *               1 0
            // 37                    0
            // 41 43 *               1 0
            // 47 53                 0 0
            // 59 61 *               1 0
            // 67                    0
            // 71 73 *               1 0
            // 79 83 89 97           0 0 0 0
            // 101 103 *             1 0
            // 107 109 *             1 0
            // 113 121 127 131       0 0 0 0
            // 137 139 *             1 0
            // 143                   0
            // 149 151 *             1 0
            // 157 163               0 0 
            // 167 169 *             1 0 
            // 173                   0
            // 179 181 *             1 0
            // 187                   0
            // 191 193 *             1 0
            // 197 199 *             1 0
            // 209 *                 1

            // can skip classes: 5 8 11 12 15 18 19 20 21 26 27 28 29 32 35 36 39 42 

            // 48-bit mask = 0xA9224141224A;

            uint64_t twins = aligned_flags1 & (aligned_flags1 >> 1);
            twins &= 0xA9224141224Aull;
            pcount += _mm_popcnt_u64(twins);
            pcount += (last_bit & aligned_flags1);
            last_bit = (aligned_flags1 >> 47);

            twins = aligned_flags2 & (aligned_flags2 >> 1);
            twins &= 0xA9224141224Aull;
            pcount += _mm_popcnt_u64(twins);
            pcount += (last_bit & aligned_flags2);
            last_bit = (aligned_flags2 >> 47);

            twins = aligned_flags3 & (aligned_flags3 >> 1);
            twins &= 0xA9224141224Aull;
            pcount += _mm_popcnt_u64(twins);
            pcount += (last_bit & aligned_flags3);
            last_bit = (aligned_flags3 >> 47);

            twins = aligned_flags4 & (aligned_flags4 >> 1);
            twins &= 0xA9224141224Aull;
            pcount += _mm_popcnt_u64(twins);
            pcount += (last_bit & aligned_flags4);
            last_bit = (aligned_flags4 >> 47);

            twins = aligned_flags5 & (aligned_flags5 >> 1);
            twins &= 0xA9224141224Aull;
            pcount += _mm_popcnt_u64(twins);
            pcount += (last_bit & aligned_flags5);
            last_bit = (aligned_flags5 >> 47);

            twins = aligned_flags6 & (aligned_flags6 >> 1);
            twins &= 0xA9224141224Aull;
            pcount += _mm_popcnt_u64(twins);
            pcount += (last_bit & aligned_flags6);
            last_bit = (aligned_flags6 >> 47);

            twins = aligned_flags7 & (aligned_flags7 >> 1);
            twins &= 0xA9224141224Aull;
            pcount += _mm_popcnt_u64(twins);
            pcount += (last_bit & aligned_flags7);
            last_bit = (aligned_flags7 >> 47);

            twins = aligned_flags8 & (aligned_flags8 >> 1);
            twins &= 0xA9224141224Aull;
            pcount += _mm_popcnt_u64(twins);
            pcount += (last_bit & aligned_flags8);
            last_bit = (aligned_flags8 >> 47);
#endif
            
        }
    }
    else if (nc == 270)
    {
        /*
        
        twin residues mod 2310 (separated into groups of 64):
        1 17 19 29 31 41 43 59 61 71 73 101 103 107 109 137 139 149 151 167 169 179 181 191 193 
        197 199 221 223 227 229 239 241 269 271 281 283 311 313 347 349 359 361 377 379 389 391 
        401 403 419 421 431 433 437 439 461 463 479 481 491 493 521 523 527 

        529 557 559 569 571 587 589 599 601 611 613 617 619 629 631 641 643 659 661 689 691 701 
        703 731 733 767 769 797 799 809 811 821 823 827 829 839 841 851 853 857 859 881 883 899 
        901 941 943 947 949 989 991 1007 1009 1019 1021 1031 1033 1037 1039 1049 1051 1061 1063 1079 

        1081 1091 1093 1121 1123 1151 1153 1157 1159 1187 1189 1217 1219 1229 1231 1247 1249 1259 
        1261 1271 1273 1277 1279 1289 1291 1301 1303 1319 1321 1361 1363 1367 1369 1409 1411 1427 
        1429 1451 1453 1457 1459 1469 1471 1481 1483 1487 1489 1499 1501 1511 1513 1541 1543 1577 
        1579 1607 1609 1619 1621 1649 1651 1667 1669 1679 

        1681 1691 1693 1697 1699 1709 1711 1721 1723 1739 1741 1751 1753 1781 1783 1787 1789 1817 
        1819 1829 1831 1847 1849 1871 1873 1877 1879 1889 1891 1907 1909 1919 1921 1931 1933 1949 
        1951 1961 1963 1997 1999 2027 2029 2039 2041 2069 2071 2081 2083 2087 2089 2111 2113 2117 
        2119 2129 2131 2141 2143 2159 2161 2171 2173 2201 

        2203 2207 2209 2237 2239 2249 2251 2267 2269 2279 2281 2291 2293 2309

        */
        
        int i;
        uint32_t last_bit = 0;

        // sort this column of 270x64 bits
        for (i = 0; i < 8; i++)
        {
            // take one byte-column at a time (270x8 bits)
            uint64_t aligned_flags[34];
            int j;

            // first we partially order the lines such
            // that each 64-bit flags contains 8 ordered
            // bytes for a set of 8 classes.
            for (j = 0; j < 33; j++)
            {
                aligned_flags[j] = interleave_pdep_8x8(
                    lines[(j * 8) + 0][byte_offset + i],
                    lines[(j * 8) + 1][byte_offset + i],
                    lines[(j * 8) + 2][byte_offset + i],
                    lines[(j * 8) + 3][byte_offset + i],
                    lines[(j * 8) + 4][byte_offset + i],
                    lines[(j * 8) + 5][byte_offset + i],
                    lines[(j * 8) + 6][byte_offset + i],
                    lines[(j * 8) + 7][byte_offset + i]);
            }

            aligned_flags[j] = interleave_pdep_8x8(
                lines[(j * 8) + 0][byte_offset + i],
                lines[(j * 8) + 1][byte_offset + i],
                lines[(j * 8) + 2][byte_offset + i],
                lines[(j * 8) + 3][byte_offset + i],
                lines[(j * 8) + 4][byte_offset + i],
                lines[(j * 8) + 5][byte_offset + i],
                0,
                0);

            // now, shuffle the bytes within the partially ordered chunks.
            // every 270 bits there will be a discontinuity between columns
            // of ordered bits.  Rather than do the somewhat expensive
            // rotation like in the 30-class case we just group each 270-bit
            // column into its own set of 5 64-bit variables and then
            // carry propagate between the sets for each of the 8 columns 
            // in this byte_offset.  5 * 8 = 40 bytes per set. 
            uint8_t *unordered_bytes = (uint8_t*)aligned_flags;
            uint8_t ordered_bytes[320];
            uint64_t* ordered64 = (uint64_t*)ordered_bytes;
           
            for (j = 0; j < 8; j++)
            {
                ordered_bytes[(j * 40) + 0 ] = unordered_bytes[j + (0  * 8)];
                ordered_bytes[(j * 40) + 1 ] = unordered_bytes[j + (1  * 8)];
                ordered_bytes[(j * 40) + 2 ] = unordered_bytes[j + (2  * 8)];
                ordered_bytes[(j * 40) + 3 ] = unordered_bytes[j + (3  * 8)];
                ordered_bytes[(j * 40) + 4 ] = unordered_bytes[j + (4  * 8)];
                ordered_bytes[(j * 40) + 5 ] = unordered_bytes[j + (5  * 8)];
                ordered_bytes[(j * 40) + 6 ] = unordered_bytes[j + (6  * 8)];
                ordered_bytes[(j * 40) + 7 ] = unordered_bytes[j + (7  * 8)];
                ordered_bytes[(j * 40) + 8 ] = unordered_bytes[j + (8  * 8)];
                ordered_bytes[(j * 40) + 9 ] = unordered_bytes[j + (9  * 8)];
                ordered_bytes[(j * 40) + 10] = unordered_bytes[j + (10 * 8)];
                ordered_bytes[(j * 40) + 11] = unordered_bytes[j + (11 * 8)];
                ordered_bytes[(j * 40) + 12] = unordered_bytes[j + (12 * 8)];
                ordered_bytes[(j * 40) + 13] = unordered_bytes[j + (13 * 8)];
                ordered_bytes[(j * 40) + 14] = unordered_bytes[j + (14 * 8)];
                ordered_bytes[(j * 40) + 15] = unordered_bytes[j + (15 * 8)];
                ordered_bytes[(j * 40) + 16] = unordered_bytes[j + (16 * 8)];
                ordered_bytes[(j * 40) + 17] = unordered_bytes[j + (17 * 8)];
                ordered_bytes[(j * 40) + 18] = unordered_bytes[j + (18 * 8)];
                ordered_bytes[(j * 40) + 19] = unordered_bytes[j + (19 * 8)];
                ordered_bytes[(j * 40) + 20] = unordered_bytes[j + (20 * 8)];
                ordered_bytes[(j * 40) + 21] = unordered_bytes[j + (21 * 8)];
                ordered_bytes[(j * 40) + 22] = unordered_bytes[j + (22 * 8)];
                ordered_bytes[(j * 40) + 23] = unordered_bytes[j + (23 * 8)];
                ordered_bytes[(j * 40) + 24] = unordered_bytes[j + (24 * 8)];
                ordered_bytes[(j * 40) + 25] = unordered_bytes[j + (25 * 8)];
                ordered_bytes[(j * 40) + 26] = unordered_bytes[j + (26 * 8)];
                ordered_bytes[(j * 40) + 27] = unordered_bytes[j + (27 * 8)];
                ordered_bytes[(j * 40) + 28] = unordered_bytes[j + (28 * 8)];
                ordered_bytes[(j * 40) + 29] = unordered_bytes[j + (29 * 8)];
                ordered_bytes[(j * 40) + 30] = unordered_bytes[j + (30 * 8)];
                ordered_bytes[(j * 40) + 31] = unordered_bytes[j + (31 * 8)];
                ordered_bytes[(j * 40) + 32] = unordered_bytes[j + (32 * 8)];
                ordered_bytes[(j * 40) + 33] = unordered_bytes[j + (33 * 8)];
            }

            for (j = 0; j < 8; j++)
            {
                // count the twins
                uint64_t twins = ordered64[(j * 5) + 0] & (ordered64[(j * 5) + 0] >> 1);
                twins &= 0x2AAAAAAAAAAAAAAAull;
                pcount += _mm_popcnt_u64(twins);
                pcount += (last_bit & ordered64[(j * 5) + 0]);  // from previous chunk
                last_bit = (ordered64[(j * 5) + 0] >> 63);      // generate carry

                twins = ordered64[(j * 5) + 1] & (ordered64[(j * 5) + 1] >> 1);
                twins &= 0x2AAAAAAAAAAAAAAAull;
                pcount += _mm_popcnt_u64(twins);
                pcount += (last_bit & ordered64[(j * 5) + 1]);  // from previous chunk
                last_bit = (ordered64[(j * 5) + 1] >> 63);      // generate carry

                twins = ordered64[(j * 5) + 2] & (ordered64[(j * 5) + 2] >> 1);
                twins &= 0x2AAAAAAAAAAAAAAAull;
                pcount += _mm_popcnt_u64(twins);
                pcount += (last_bit & ordered64[(j * 5) + 2]);  // from previous chunk
                last_bit = (ordered64[(j * 5) + 2] >> 63);      // generate carry

                twins = ordered64[(j * 5) + 3] & (ordered64[(j * 5) + 3] >> 1);
                twins &= 0x2AAAAAAAAAAAAAAAull;
                pcount += _mm_popcnt_u64(twins);
                pcount += (last_bit & ordered64[(j * 5) + 3]);  // from previous chunk
                last_bit = (ordered64[(j * 5) + 3] >> 63);      // generate carry

                // 14 bits left in this last word and, similar to the previous words,
                // we mask off the last bit because it gets carried to the next chunk.
                twins = ordered64[(j * 5) + 4] & (ordered64[(j * 5) + 4] >> 1);
                twins &= 0x0AAAull;
                pcount += _mm_popcnt_u64(twins);
                pcount += (last_bit & ordered64[(j * 5) + 4]);          // from previous chunk
                last_bit = ((ordered64[(j * 5) + 4] & 0x2000) >> 13);   // generate carry
            }
        }

    }
    else
    {
        uint64_t** pqueues; // [64][48]
        uint32_t pcounts[64];
        int i, j;
        uint64_t prodN = sdata->prodN;
        uint32_t current_line;

        // ordering the bits becomes inefficient with 48 lines because
        // they would need to be dispersed over too great a distance.
        // instead we compute the primes as before but push the results
        // into 64 different queues depending on the bit position.  Then
        // we pull from the queues in order while storing into the primes array.
        // This time the bottleneck is mostly in the queue-based sorting
        // and associated memory operations, so we don't bother with
        // switching between branch-free inner loops or not.
        pqueues = (uint64_t**)xmalloc(64 * sizeof(uint64_t*));
        for (i = 0; i < 64; i++)
        {
            pqueues[i] = (uint64_t*)xmalloc(sdata->numclasses * sizeof(uint64_t));
        }

        // Compute the primes using ctz on the 64-bit words but push the results
        // into 64 different queues depending on the bit position.  Then
        // we pull from the queues in order while storing into the primes array.
        // This time the bottleneck is mostly in the queue-based sorting
        // and associated memory operations, so we don't bother with
        // switching between branch-free inner loops or not.
        memset(pcounts, 0, 64 * sizeof(uint32_t));

        lowlimit += byte_offset * 8 * prodN;
        for (current_line = 0; current_line < nc; current_line++)
        {
            uint64_t* line64 = (uint64_t*)lines[current_line];
            uint64_t flags64 = line64[byte_offset / 8];

            while (flags64 > 0)
            {
                uint64_t pos = _trail_zcnt64(flags64);
                uint64_t prime = lowlimit + pos * prodN + sdata->rclass[current_line];

                pqueues[pos][pcounts[pos]] = prime;
                pcounts[pos]++;
                flags64 ^= (1ULL << pos);
            }
        }

        for (i = 0; i < 64; i++)
        {
            // search for twins and only load the leading element/prime.
            // if depth-based sieving then these are candidate twins.
            if (pcounts[i] > 0)
            {
                for (j = 0; j < pcounts[i] - 1; j++)
                {
                    if ((pqueues[i][j + 1] - pqueues[i][j]) == 2)
                    {
                        pcount++;
                    }
                }
                if (i < 63)
                {
                    if (pcounts[i + 1] > 0)
                    {
                        if ((pqueues[i + 1][0] - pqueues[i][j]) == 2)
                        {
                            pcount++;
                        }
                    }
                }
            }

        }

        for (i = 0; i < 64; i++)
        {
            free(pqueues[i]);
        }
        free(pqueues);


    }

    return pcount;
}

#endif
