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


#include "microecm.h"
#include "monty.h"

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#if defined(_MSC_VER)
#  ifndef _WIN64
#    error "64 bit compilation mode is required for MSVC"
#  endif
#  include <immintrin.h>
#  include <intrin.h>
#endif


#ifdef USE_AVX512F
#  include <immintrin.h>
#endif


// Using the inline asm in this file can increase performance by ~20-25%
// (surprisingly).  Hence these macros are defined by default.
#if defined(__x86_64__) || defined(_M_X64)
#  define MICRO_ECM_ALT_MULREDC_USE_INLINE_ASM_X86
#endif

// We rarely will ever want to use debugging printfs
//#define MICRO_ECM_VERBOSE_PRINTF




#ifdef _MSC_VER
#  define MICRO_ECM_FORCE_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
#  define MICRO_ECM_FORCE_INLINE inline __attribute__((always_inline))
#else
#  define MICRO_ECM_FORCE_INLINE __inline
#endif


#ifdef MICRO_ECM_VERBOSE_PRINTF
#  include <stdio.h>
// these globals will cause data races, if we run multithreaded
    uint32_t stg1Doub;
    uint32_t stg1Add;
    uint32_t ptadds;
#endif


typedef struct
{
    uint64_t X;
    uint64_t Z;
} uecm_pt;

#ifdef USE_AVX512F
typedef struct
{
    __m512i X;
    __m512i Z;
} uecm_pt_x8;
#endif


static const uint32_t map[60] = {
    0, 1, 2, 0, 0, 0, 0, 3, 0, 0,
    0, 4, 0, 5, 0, 0, 0, 6, 0, 7,
    0, 0, 0, 8, 0, 0, 0, 0, 0, 9,
    0, 10, 0, 0, 0, 0, 0, 11, 0, 0,
    0, 12, 0, 13, 0, 0, 0, 14, 0, 15,
    0, 0, 0, 16, 0, 0, 0, 0, 0, 17 };

static const double INV_2_POW_32 = 1.0 / (double)((uint64_t)(1) << 32);

#define NUMP 801
static const int primes[NUMP] = {
2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
79, 83, 89, 97, 101, 103, 107, 109, 113, 127,
131, 137, 139, 149, 151, 157, 163, 167, 173, 179,
181, 191, 193, 197, 199, 211, 223, 227, 229, 233,
239, 241, 251, 257, 263, 269, 271, 277, 281, 283,
293, 307, 311, 313, 317, 331, 337, 347, 349, 353,
359, 367, 373, 379, 383, 389, 397, 401, 409, 419,
421, 431, 433, 439, 443, 449, 457, 461, 463, 467,
479, 487, 491, 499, 503, 509, 521, 523, 541, 547,
557, 563, 569, 571, 577, 587, 593, 599, 601, 607,
613, 617, 619, 631, 641, 643, 647, 653, 659, 661,
673, 677, 683, 691, 701, 709, 719, 727, 733, 739,
743, 751, 757, 761, 769, 773, 787, 797, 809, 811,
821, 823, 827, 829, 839, 853, 857, 859, 863, 877,
881, 883, 887, 907, 911, 919, 929, 937, 941, 947,
953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019,
1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087,
1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153,
1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229,
1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297,
1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381,
1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453,
1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523,
1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597,
1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663,
1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741,
1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823,
1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901,
1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993,
1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063,
2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, 2131,
2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221,
2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, 2293,
2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 2371,
2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437,
2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539,
2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621,
2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 2689,
2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741, 2749,
2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 2833,
2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 2909,
2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 3001,
3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 3083,
3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, 3187,
3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 3259,
3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343,
3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, 3433,
3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517,
3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571, 3581,
3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659,
3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, 3733,
3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, 3823,
3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 3911,
3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, 4001,
4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, 4073,
4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, 4153,
4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, 4241,
4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, 4327,
4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409, 4421,
4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, 4507,
4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, 4591,
4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, 4663,
4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, 4759,
4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 4861,
4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, 4943,
4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, 5009,
5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, 5099,
5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, 5189,
5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279, 5281,
5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, 5393,
5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, 5449,
5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, 5527,
5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, 5641,
5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, 5701,
5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, 5801,
5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, 5861,
5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, 5953,
5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, 6067,
6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133, 6143,
};

static uint32_t uecm_lcg_rand_32B(uint32_t lower, uint32_t upper, uint64_t *ploc_lcg)
{
    *ploc_lcg = 6364136223846793005ULL * (*ploc_lcg) + 1442695040888963407ULL;
    return lower + (uint32_t)(
        (double)(upper - lower) * (double)((*ploc_lcg) >> 32) * INV_2_POW_32);
}



/* --- The following two functions are written by Jeff Hurchalla, Copyright 2022 --- */

// for this algorithm, see https://jeffhurchalla.com/2022/04/28/montgomery-redc-using-the-positive-inverse-mod-r/
MICRO_ECM_FORCE_INLINE static uint64_t uecm_mulredc_alt(uint64_t x, uint64_t y, uint64_t N, uint64_t invN)
{
#if defined(_MSC_VER)
    uint64_t T_hi;
    uint64_t T_lo = _umul128(x, y, &T_hi);
    uint64_t m = T_lo * invN;
    uint64_t mN_hi = __umulh(m, N);
#else
    __uint128_t prod = (__uint128_t)x * y;
    uint64_t T_hi = (uint64_t)(prod >> 64);
    uint64_t T_lo = (uint64_t)(prod);
    uint64_t m = T_lo * invN;
    __uint128_t mN = (__uint128_t)m * N;
    uint64_t mN_hi = (uint64_t)(mN >> 64);
#endif
    uint64_t tmp = T_hi + N;
#if defined(MICRO_ECM_ALT_MULREDC_USE_INLINE_ASM_X86) && !defined(_MSC_VER)
    __asm__ (
        "subq %[mN_hi], %[tmp] \n\t"    /* tmp = T_hi + N - mN_hi */
        "subq %[mN_hi], %[T_hi] \n\t"   /* T_hi = T_hi - mN_hi */
        "cmovaeq %[T_hi], %[tmp] \n\t"  /* tmp = (T_hi >= mN_hi) ? T_hi : tmp */
        : [tmp]"+&r"(tmp), [T_hi]"+&r"(T_hi)
        : [mN_hi]"r"(mN_hi)
        : "cc");
    uint64_t result = tmp;
#else
    tmp = tmp - mN_hi;
    uint64_t result = T_hi - mN_hi;
    result = (T_hi < mN_hi) ? tmp : result;
#endif
   return result;
}

// for this algorithm, see https://jeffhurchalla.com/2022/04/25/a-faster-multiplicative-inverse-mod-a-power-of-2/
static uint64_t uecm_multiplicative_inverse(uint64_t a)
{
//    assert(a%2 == 1);  // the inverse (mod 2<<64) only exists for odd values
    uint64_t x0 = (3*a)^2;
    uint64_t y = 1 - a*x0;
    uint64_t x1 = x0*(1 + y);
    y *= y;
    uint64_t x2 = x1*(1 + y);
    y *= y;
    uint64_t x3 = x2*(1 + y);
    y *= y;
    uint64_t x4 = x3*(1 + y);
    return x4;
}

/* --- end Hurchalla functions --- */


// full strength mul/sqr redc
MICRO_ECM_FORCE_INLINE static uint64_t uecm_mulredc(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
{
    return uecm_mulredc_alt(x, y, n, 0 - nhat);
}
MICRO_ECM_FORCE_INLINE static uint64_t uecm_sqrredc(uint64_t x, uint64_t n, uint64_t nhat)
{
    return uecm_mulredc_alt(x, x, n, 0 - nhat);
}

MICRO_ECM_FORCE_INLINE static void uecm_uadd(uint64_t rho, uint64_t n, const uecm_pt P1, const uecm_pt P2,
    const uecm_pt Pin, uecm_pt *Pout)
{
    // compute:
    //x+ = z- * [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
    //z+ = x- * [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
    // where:
    //x- = original x
    //z- = original z
    // given the sums and differences of the original points
    uint64_t diff1 = submod(P1.X, P1.Z, n);
    uint64_t sum1 = addmod(P1.X, P1.Z, n);
    uint64_t diff2 = submod(P2.X, P2.Z, n);
    uint64_t sum2 = addmod(P2.X, P2.Z, n);

    uint64_t tt1 = uecm_mulredc(diff1, sum2, n, rho); //U
    uint64_t tt2 = uecm_mulredc(sum1, diff2, n, rho); //V

    uint64_t tt3 = addmod(tt1, tt2, n);
    uint64_t tt4 = submod(tt1, tt2, n);
    tt1 = uecm_sqrredc(tt3, n, rho);   //(U + V)^2
    tt2 = uecm_sqrredc(tt4, n, rho);   //(U - V)^2

    uint64_t tmpx = uecm_mulredc(tt1, Pin.Z, n, rho);     //Z * (U + V)^2
    uint64_t tmpz = uecm_mulredc(tt2, Pin.X, n, rho);     //x * (U - V)^2
    Pout->X = tmpx;
    Pout->Z = tmpz;

#ifdef MICRO_ECM_VERBOSE_PRINTF
    stg1Add++;
#endif
    return;
}

MICRO_ECM_FORCE_INLINE static void uecm_udup(uint64_t s, uint64_t rho, uint64_t n,
    uint64_t insum, uint64_t indiff, uecm_pt *P)
{
    uint64_t tt1 = uecm_sqrredc(indiff, n, rho);          // U=(x1 - z1)^2
    uint64_t tt2 = uecm_sqrredc(insum, n, rho);           // V=(x1 + z1)^2
    P->X = uecm_mulredc(tt1, tt2, n, rho);         // x=U*V

    uint64_t tt3 = submod(tt2, tt1, n);          // w = V-U
    tt2 = uecm_mulredc(tt3, s, n, rho);      // w = (A+2)/4 * w
    tt2 = addmod(tt2, tt1, n);          // w = w + U
    P->Z = uecm_mulredc(tt2, tt3, n, rho);         // Z = w*(V-U)
#ifdef MICRO_ECM_VERBOSE_PRINTF
    stg1Doub++;
#endif
    return;
}


#ifdef USE_AVX512F
#define uecm_and64 _mm512_and_epi64
#define uecm_storeu64 _mm512_store_epi64
#define uecm_add64 _mm512_add_epi64
#define uecm_sub64 _mm512_sub_epi64
#define uecm_set64 _mm512_set1_epi64
#define uecm_srli64 _mm512_srli_epi64
#define uecm_loadu64 _mm512_load_epi64
#define uecm_castpd _mm512_castsi512_pd
#define uecm_castepu _mm512_castpd_si512
#define uecm_DIGIT_SIZE 52
#define uecm_DIGIT_MASK 0x000fffffffffffffULL
#define uecm_MB_WIDTH 8
#define uecm_SIMD_BYTES 64

static __m512d dbias;
static __m512i vbias1;
static __m512i vbias2;
static __m512i lo52mask;

MICRO_ECM_FORCE_INLINE static __m512i uecm_mul52lo(__m512i b, __m512i c)
{
    return _mm512_and_si512(_mm512_mullo_epi64(b, c), _mm512_set1_epi64(0x000fffffffffffffull));
}
MICRO_ECM_FORCE_INLINE static __m512i uecm_mul52hi(__m512i b, __m512i c)
{
    __m512d prod1_ld = _mm512_cvtepu64_pd(b);
    __m512d prod2_ld = _mm512_cvtepu64_pd(c);
    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
    return _mm512_sub_epi64(uecm_castepu(prod1_ld), vbias1);
}
MICRO_ECM_FORCE_INLINE static void uecm_mul52lohi(__m512i b, __m512i c, __m512i* l, __m512i* h)
{
    __m512d prod1_ld = _mm512_cvtepu64_pd(b);
    __m512d prod2_ld = _mm512_cvtepu64_pd(c);
    __m512d prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
    *h = _mm512_sub_epi64(uecm_castepu(prod1_hd), vbias1);
    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
    *l = _mm512_castpd_si512(prod1_ld);
    return;
}
MICRO_ECM_FORCE_INLINE static __m512i uecm_mulredc_x8(__m512i x, __m512i y, __m512i N, __m512i invN)
{
    // invN is the positive variant = 0 - nhat (the standard negative inverse)
    __m512i T_hi; // = uecm_mul52hi(x, y, dbias, vbias1);
    __m512i T_lo; // = uecm_mul52lo(x, y);
    uecm_mul52lohi(x, y, &T_lo, &T_hi);
    __m512i m = _mm512_and_si512(_mm512_mullo_epi64(T_lo, invN), lo52mask);
    __m512i mN_hi = uecm_mul52hi(m, N);
    __m512i tmp = uecm_add64(T_hi, N);
    tmp = uecm_sub64(tmp, mN_hi);
    __m512i result = uecm_sub64(T_hi, mN_hi);
    __mmask8 msk = _mm512_cmplt_epu64_mask(T_hi, mN_hi);
    return _mm512_mask_mov_epi64(result, msk, tmp);
}
MICRO_ECM_FORCE_INLINE static __m512i uecm_sqrredc_x8(__m512i x, __m512i N, __m512i invN)
{
    // invN is the positive variant = 0 - nhat (the standard negative inverse)
    __m512i T_hi; // = uecm_mul52hi(x, y, dbias, vbias1);
    __m512i T_lo; // = uecm_mul52lo(x, y);
    uecm_mul52lohi(x, x, &T_lo, &T_hi);
    __m512i m = _mm512_and_si512(_mm512_mullo_epi64(T_lo, invN), lo52mask);
    __m512i mN_hi = uecm_mul52hi(m, N);
    __m512i tmp = uecm_add64(T_hi, N);
    tmp = uecm_sub64(tmp, mN_hi);
    __m512i result = uecm_sub64(T_hi, mN_hi);
    __mmask8 msk = _mm512_cmplt_epu64_mask(T_hi, mN_hi);
    return _mm512_mask_mov_epi64(result, msk, tmp);
}
MICRO_ECM_FORCE_INLINE static __m512i uecm_addmod_x8(__m512i x, __m512i y, __m512i N)
{
    __m512i t = uecm_add64(x, y);
    __mmask8 m = _mm512_cmpge_epu64_mask(t, N);
    return _mm512_mask_sub_epi64(t, m, t, N);
}
MICRO_ECM_FORCE_INLINE static __m512i uecm_submod_x8(__m512i x, __m512i y, __m512i N)
{
    __m512i t = uecm_sub64(x, y);
    __mmask8 m = _mm512_cmpgt_epu64_mask(y, x);
    return _mm512_mask_add_epi64(t, m, t, N);
}

static uint64_t uecm_multiplicative_inverse52(uint64_t a)
{
    //    assert(a%2 == 1);  // the inverse (mod 2<<64) only exists for odd values
    uint64_t x0 = (3 * a) ^ 2;
    uint64_t y = 1 - a * x0;
    uint64_t x1 = x0 * (1 + y);
    y *= y;
    uint64_t x2 = x1 * (1 + y);
    y *= y;
    uint64_t x3 = x2 * (1 + y);
    y *= y;
    uint64_t x4 = x3 * (1 + y);
    return x4 & 0x000fffffffffffffull;
}

MICRO_ECM_FORCE_INLINE uint64_t submod52(uint64_t a, uint64_t b, uint64_t n)
{
    uint64_t r0;
    if (_subborrow_u64(0, a, b, &r0))
        r0 += n;
    return r0;
}
MICRO_ECM_FORCE_INLINE uint64_t addmod52(uint64_t x, uint64_t y, uint64_t n)
{
    // FYI: The clause above often compiles with a branch in MSVC.
    // The statement below often compiles without a branch (uses cmov) in MSVC.
    return (x >= n - y) ? x - (n - y) : x + y;
}
MICRO_ECM_FORCE_INLINE static uint64_t uecm_mulredc52(uint64_t x, uint64_t y, uint64_t N, uint64_t invN)
{
#if defined(_MSC_VER)
    uint64_t T_hi;
    uint64_t T_lo = _umul128(x, y, &T_hi);
    uint64_t m = T_lo * invN;
    uint64_t mN_hi = __umulh(m, N);
#else
    __uint128_t prod = (__uint128_t)x * y;
    uint64_t T_hi = (uint64_t)(prod >> 52);
    uint64_t T_lo = (uint64_t)(prod & 0x000fffffffffffffull);
    uint64_t m = (T_lo * invN) & 0x000fffffffffffffull;
    __uint128_t mN = (__uint128_t)m * N;
    uint64_t mN_hi = (uint64_t)(mN >> 52);
#endif
    uint64_t tmp = T_hi + N;
#if defined(MICRO_ECM_ALT_MULREDC_USE_INLINE_ASM_X86) && !defined(_MSC_VER)
    __asm__(
        "subq %[mN_hi], %[tmp] \n\t"    /* tmp = T_hi + N - mN_hi */
        "subq %[mN_hi], %[T_hi] \n\t"   /* T_hi = T_hi - mN_hi */
        "cmovaeq %[T_hi], %[tmp] \n\t"  /* tmp = (T_hi >= mN_hi) ? T_hi : tmp */
        : [tmp] "+&r"(tmp), [T_hi]"+&r"(T_hi)
        : [mN_hi] "r"(mN_hi)
        : "cc");
    uint64_t result = tmp;
#else
    tmp = tmp - mN_hi;
    uint64_t result = T_hi - mN_hi;
    result = (T_hi < mN_hi) ? tmp : result;
#endif
    return result & 0x000fffffffffffffull;
}
MICRO_ECM_FORCE_INLINE static uint64_t uecm_sqrredc52(uint64_t x, uint64_t N, uint64_t invN)
{
    return uecm_mulredc52(x, x, N, invN);
}

MICRO_ECM_FORCE_INLINE static void uecm_uadd_x8(__m512i rho, __m512i n, 
    const uecm_pt_x8 P1, const uecm_pt_x8 P2,
    const uecm_pt_x8 Pin, uecm_pt_x8* Pout)
{
    // compute:
    //x+ = z- * [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
    //z+ = x- * [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
    // where:
    //x- = original x
    //z- = original z
    // given the sums and differences of the original points
    __m512i diff1 = uecm_submod_x8(P1.X, P1.Z, n);
    __m512i sum1 = uecm_addmod_x8(P1.X, P1.Z, n);
    __m512i diff2 = uecm_submod_x8(P2.X, P2.Z, n);
    __m512i sum2 = uecm_addmod_x8(P2.X, P2.Z, n);

    __m512i tt1 = uecm_mulredc_x8(diff1, sum2, n, rho); //U
    __m512i tt2 = uecm_mulredc_x8(sum1, diff2, n, rho); //V

    __m512i tt3 = uecm_addmod_x8(tt1, tt2, n);
    __m512i tt4 = uecm_submod_x8(tt1, tt2, n);
    tt1 = uecm_sqrredc_x8(tt3, n, rho);   //(U + V)^2
    tt2 = uecm_sqrredc_x8(tt4, n, rho);   //(U - V)^2

    __m512i tmpx = uecm_mulredc_x8(tt1, Pin.Z, n, rho);     //Z * (U + V)^2
    __m512i tmpz = uecm_mulredc_x8(tt2, Pin.X, n, rho);     //x * (U - V)^2
    Pout->X = tmpx;
    Pout->Z = tmpz;

    return;
}

MICRO_ECM_FORCE_INLINE static void uecm_udup_x8(__m512i s, __m512i rho, __m512i n,
    __m512i insum, __m512i indiff, uecm_pt_x8* P)
{
    __m512i tt1 = uecm_sqrredc_x8(indiff, n, rho);          // U=(x1 - z1)^2
    __m512i tt2 = uecm_sqrredc_x8(insum, n, rho);           // V=(x1 + z1)^2
    P->X = uecm_mulredc_x8(tt1, tt2, n, rho);         // x=U*V

    __m512i tt3 = uecm_submod_x8(tt2, tt1, n);          // w = V-U
    tt2 = uecm_mulredc_x8(tt3, s, n, rho);      // w = (A+2)/4 * w
    tt2 = uecm_addmod_x8(tt2, tt1, n);          // w = w + U
    P->Z = uecm_mulredc_x8(tt2, tt3, n, rho);         // Z = w*(V-U)
    return;
}

static void uecm_uprac70_x8(__m512i rho, __m512i n, uecm_pt_x8* P, __m512i s)
{
    __m512i s1, s2, d1, d2;
    __m512i swp;
    int i;
    static const uint8_t steps[116] = {
        0,6,0,6,0,6,0,4,6,0,4,6,0,4,4,6,
        0,4,4,6,0,5,4,6,0,3,3,4,6,0,3,5,
        4,6,0,3,4,3,4,6,0,5,5,4,6,0,5,3,
        3,4,6,0,3,3,4,3,4,6,0,5,3,3,3,3,
        3,3,3,3,4,3,3,4,6,0,5,4,3,3,4,6,
        0,3,4,3,5,4,6,0,5,3,3,3,4,6,0,5,
        4,3,5,4,6,0,5,5,3,3,4,6,0,4,3,3,
        3,5,4,6 };

    uecm_pt_x8 pt1, pt2, pt3;
    for (i = 0; i < 116; i++)
    {
        if (steps[i] == 0)
        {
            pt1.X = pt2.X = pt3.X = P->X;
            pt1.Z = pt2.Z = pt3.Z = P->Z;

            d1 = uecm_submod_x8(pt1.X, pt1.Z, n);
            s1 = uecm_addmod_x8(pt1.X, pt1.Z, n);
            uecm_udup_x8(s, rho, n, s1, d1, &pt1);
        }
        else if (steps[i] == 3)
        {
            // integrate step 4 followed by swap(1,2)
            uecm_pt_x8 pt4;
            uecm_uadd_x8(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

            swp = pt1.X;
            pt1.X = pt4.X;
            pt4.X = pt3.X;
            pt3.X = pt2.X;
            pt2.X = swp;
            swp = pt1.Z;
            pt1.Z = pt4.Z;
            pt4.Z = pt3.Z;
            pt3.Z = pt2.Z;
            pt2.Z = swp;
        }
        else if (steps[i] == 4)
        {
            uecm_pt_x8 pt4;
            uecm_uadd_x8(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

            swp = pt2.X;
            pt2.X = pt4.X;
            pt4.X = pt3.X;
            pt3.X = swp;
            swp = pt2.Z;
            pt2.Z = pt4.Z;
            pt4.Z = pt3.Z;
            pt3.Z = swp;
        }
        else if (steps[i] == 5)
        {
            d2 = uecm_submod_x8(pt1.X, pt1.Z, n);
            s2 = uecm_addmod_x8(pt1.X, pt1.Z, n);

            uecm_uadd_x8(rho, n, pt2, pt1, pt3, &pt2);        // B = B + A (C)
            uecm_udup_x8(s, rho, n, s2, d2, &pt1);        // A = 2A
        }
        else if (steps[i] == 6)
        {
            uecm_uadd_x8(rho, n, pt1, pt2, pt3, P);     // A = A + B (C)
        }
    }
    return;
}

static void uecm_uprac85_x8(__m512i rho, __m512i n, uecm_pt_x8* P, __m512i s)
{
    __m512i s1, s2, d1, d2;
    __m512i swp;
    int i;
    static const uint8_t steps[146] = {
        0,6,0,6,0,6,0,6,0,4,
        6,0,4,6,0,4,4,6,0,4,
        4,6,0,5,4,6,0,3,3,4,
        6,0,3,5,4,6,0,3,4,3,
        4,6,0,5,5,4,6,0,5,3,
        3,4,6,0,3,3,4,3,4,6,
        0,4,3,4,3,5,3,3,3,3,
        3,3,3,3,4,6,0,3,3,3,
        3,3,3,3,3,3,4,3,4,3,
        4,6,0,3,4,3,5,4,6,0,
        5,3,3,3,4,6,0,5,4,3,
        5,4,6,0,4,3,3,3,5,4,
        6,0,4,3,5,3,3,4,6,0,
        3,3,3,3,5,4,6,0,3,3,
        3,4,3,3,4,6 };

    uecm_pt_x8 pt1, pt2, pt3;
    for (i = 0; i < 146; i++)
    {
        if (steps[i] == 0)
        {
            pt1.X = pt2.X = pt3.X = P->X;
            pt1.Z = pt2.Z = pt3.Z = P->Z;

            d1 = uecm_submod_x8(pt1.X, pt1.Z, n);
            s1 = uecm_addmod_x8(pt1.X, pt1.Z, n);
            uecm_udup_x8(s, rho, n, s1, d1, &pt1);
        }
        else if (steps[i] == 3)
        {
            // integrate step 4 followed by swap(1,2)
            uecm_pt_x8 pt4;
            uecm_uadd_x8(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

            swp = pt1.X;
            pt1.X = pt4.X;
            pt4.X = pt3.X;
            pt3.X = pt2.X;
            pt2.X = swp;
            swp = pt1.Z;
            pt1.Z = pt4.Z;
            pt4.Z = pt3.Z;
            pt3.Z = pt2.Z;
            pt2.Z = swp;
        }
        else if (steps[i] == 4)
        {
            uecm_pt_x8 pt4;
            uecm_uadd_x8(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

            swp = pt2.X;
            pt2.X = pt4.X;
            pt4.X = pt3.X;
            pt3.X = swp;
            swp = pt2.Z;
            pt2.Z = pt4.Z;
            pt4.Z = pt3.Z;
            pt3.Z = swp;
        }
        else if (steps[i] == 5)
        {
            d2 = uecm_submod_x8(pt1.X, pt1.Z, n);
            s2 = uecm_addmod_x8(pt1.X, pt1.Z, n);

            uecm_uadd_x8(rho, n, pt2, pt1, pt3, &pt2);        // B = B + A (C)
            uecm_udup_x8(s, rho, n, s2, d2, &pt1);        // A = 2A
        }
        else if (steps[i] == 6)
        {
            uecm_uadd_x8(rho, n, pt1, pt2, pt3, P);     // A = A + B (C)
        }
    }
    return;
}

static void uecm_uprac_x8(__m512i rho, __m512i n, uecm_pt_x8* P, uint64_t c, double v, __m512i s)
{
    uint64_t d, e, r;
    int i;
    __m512i s1, s2, d1, d2;
    __m512i swp;

    // we require c != 0
    int shift = _trail_zcnt64(c);
    c = c >> shift;

    d = c;
    r = (uint64_t)((double)d * v + 0.5);

    d = c - r;
    e = 2 * r - c;

    uecm_pt_x8 pt1, pt2, pt3;

    // the first one is always a doubling
    // point1 is [1]P
    pt1.X = pt2.X = pt3.X = P->X;
    pt1.Z = pt2.Z = pt3.Z = P->Z;

    d1 = uecm_submod_x8(pt1.X, pt1.Z, n);
    s1 = uecm_addmod_x8(pt1.X, pt1.Z, n);

    // point2 is [2]P
    uecm_udup_x8(s, rho, n, s1, d1, &pt1);

    while (d != e)
    {
        if (d < e)
        {
            r = d;
            d = e;
            e = r;
            swp = pt1.X;
            pt1.X = pt2.X;
            pt2.X = swp;
            swp = pt1.Z;
            pt1.Z = pt2.Z;
            pt2.Z = swp;
        }
        if (d - e <= e / 4 && ((d + e) % 3) == 0)
        {
            d = (2 * d - e) / 3;
            e = (e - d) / 2;

            uecm_pt_x8 pt4;
            uecm_uadd_x8(rho, n, pt1, pt2, pt3, &pt4); // T = A + B (C)
            uecm_pt_x8 pt5;
            uecm_uadd_x8(rho, n, pt4, pt1, pt2, &pt5); // T2 = T + A (B)
            uecm_uadd_x8(rho, n, pt2, pt4, pt1, &pt2); // B = B + T (A)

            swp = pt1.X;
            pt1.X = pt5.X;
            pt5.X = swp;
            swp = pt1.Z;
            pt1.Z = pt5.Z;
            pt5.Z = swp;
        }
        else if (d - e <= e / 4 && (d - e) % 6 == 0)
        {
            d = (d - e) / 2;

            d1 = uecm_submod_x8(pt1.X, pt1.Z, n);
            s1 = uecm_addmod_x8(pt1.X, pt1.Z, n);

            uecm_uadd_x8(rho, n, pt1, pt2, pt3, &pt2);        // B = A + B (C)
            uecm_udup_x8(s, rho, n, s1, d1, &pt1);        // A = 2A
        }
        else if ((d + 3) / 4 <= e)
        {
            d -= e;

            uecm_pt_x8 pt4;
            uecm_uadd_x8(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

            swp = pt2.X;
            pt2.X = pt4.X;
            pt4.X = pt3.X;
            pt3.X = swp;
            swp = pt2.Z;
            pt2.Z = pt4.Z;
            pt4.Z = pt3.Z;
            pt3.Z = swp;
        }
        else if ((d + e) % 2 == 0)
        {
            d = (d - e) / 2;

            d2 = uecm_submod_x8(pt1.X, pt1.Z, n);
            s2 = uecm_addmod_x8(pt1.X, pt1.Z, n);

            uecm_uadd_x8(rho, n, pt2, pt1, pt3, &pt2);        // B = B + A (C)
            uecm_udup_x8(s, rho, n, s2, d2, &pt1);        // A = 2A
        }
        else if (d % 2 == 0)
        {
            d /= 2;

            d2 = uecm_submod_x8(pt1.X, pt1.Z, n);
            s2 = uecm_addmod_x8(pt1.X, pt1.Z, n);

            uecm_uadd_x8(rho, n, pt3, pt1, pt2, &pt3);        // C = C + A (B)
            uecm_udup_x8(s, rho, n, s2, d2, &pt1);        // A = 2A
        }
        else if (d % 3 == 0)
        {
            d = d / 3 - e;

            d1 = uecm_submod_x8(pt1.X, pt1.Z, n);
            s1 = uecm_addmod_x8(pt1.X, pt1.Z, n);

            uecm_pt_x8 pt4;
            uecm_udup_x8(s, rho, n, s1, d1, &pt4);        // T = 2A
            uecm_pt_x8 pt5;
            uecm_uadd_x8(rho, n, pt1, pt2, pt3, &pt5);        // T2 = A + B (C)
            uecm_uadd_x8(rho, n, pt4, pt1, pt1, &pt1);        // A = T + A (A)
            uecm_uadd_x8(rho, n, pt4, pt5, pt3, &pt4);        // T = T + T2 (C)

            swp = pt3.X;
            pt3.X = pt2.X;
            pt2.X = pt4.X;
            pt4.X = swp;
            swp = pt3.Z;
            pt3.Z = pt2.Z;
            pt2.Z = pt4.Z;
            pt4.Z = swp;

        }
        else if ((d + e) % 3 == 0)
        {
            d = (d - 2 * e) / 3;

            uecm_pt_x8 pt4;
            uecm_uadd_x8(rho, n, pt1, pt2, pt3, &pt4);        // T = A + B (C)


            d2 = uecm_submod_x8(pt1.X, pt1.Z, n);
            s2 = uecm_addmod_x8(pt1.X, pt1.Z, n);
            uecm_uadd_x8(rho, n, pt4, pt1, pt2, &pt2);        // B = T + A (B)
            uecm_udup_x8(s, rho, n, s2, d2, &pt4);        // T = 2A
            uecm_uadd_x8(rho, n, pt1, pt4, pt1, &pt1);        // A = A + T (A) = 3A
        }
        else if ((d - e) % 3 == 0)
        {
            d = (d - e) / 3;

            uecm_pt_x8 pt4;
            uecm_uadd_x8(rho, n, pt1, pt2, pt3, &pt4);        // T = A + B (C)

            d2 = uecm_submod_x8(pt1.X, pt1.Z, n);
            s2 = uecm_addmod_x8(pt1.X, pt1.Z, n);
            uecm_uadd_x8(rho, n, pt3, pt1, pt2, &pt3);        // C = C + A (B)

            swp = pt2.X;
            pt2.X = pt4.X;
            pt4.X = swp;
            swp = pt2.Z;
            pt2.Z = pt4.Z;
            pt4.Z = swp;

            uecm_udup_x8(s, rho, n, s2, d2, &pt4);        // T = 2A
            uecm_uadd_x8(rho, n, pt1, pt4, pt1, &pt1);        // A = A + T (A) = 3A
        }
        else
        {
            e /= 2;

            d2 = uecm_submod_x8(pt2.X, pt2.Z, n);
            s2 = uecm_addmod_x8(pt2.X, pt2.Z, n);

            uecm_uadd_x8(rho, n, pt3, pt2, pt1, &pt3);        // C = C + B (A)
            uecm_udup_x8(s, rho, n, s2, d2, &pt2);        // B = 2B
        }
    }
    uecm_uadd_x8(rho, n, pt1, pt2, pt3, P);     // A = A + B (C)

    for (i = 0; i < shift; i++)
    {
        d1 = uecm_submod_x8(P->X, P->Z, n);
        s1 = uecm_addmod_x8(P->X, P->Z, n);
        uecm_udup_x8(s, rho, n, s1, d1, P);     // P = 2P
    }
    return;
}

#endif

static void uecm_uprac70(uint64_t rho, uint64_t n, uecm_pt *P, uint64_t s)
{
    uint64_t s1, s2, d1, d2;
    uint64_t swp;
    int i;
    static const uint8_t steps[116] = {
        0,6,0,6,0,6,0,4,6,0,4,6,0,4,4,6,
        0,4,4,6,0,5,4,6,0,3,3,4,6,0,3,5,
        4,6,0,3,4,3,4,6,0,5,5,4,6,0,5,3,
        3,4,6,0,3,3,4,3,4,6,0,5,3,3,3,3,
        3,3,3,3,4,3,3,4,6,0,5,4,3,3,4,6,
        0,3,4,3,5,4,6,0,5,3,3,3,4,6,0,5,
        4,3,5,4,6,0,5,5,3,3,4,6,0,4,3,3,
        3,5,4,6 };

    uecm_pt pt1, pt2, pt3;
    for (i = 0; i < 116; i++)
    {
        if (steps[i] == 0)
        {
            pt1.X = pt2.X = pt3.X = P->X;
            pt1.Z = pt2.Z = pt3.Z = P->Z;

            d1 = submod(pt1.X, pt1.Z, n);
            s1 = addmod(pt1.X, pt1.Z, n);
            uecm_udup(s, rho, n, s1, d1, &pt1);
        }
        else if (steps[i] == 3)
        {
            // integrate step 4 followed by swap(1,2)
            uecm_pt pt4;
            uecm_uadd(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

            swp = pt1.X;
            pt1.X = pt4.X;
            pt4.X = pt3.X;
            pt3.X = pt2.X;
            pt2.X = swp;
            swp = pt1.Z;
            pt1.Z = pt4.Z;
            pt4.Z = pt3.Z;
            pt3.Z = pt2.Z;
            pt2.Z = swp;
        }
        else if (steps[i] == 4)
        {
            uecm_pt pt4;
            uecm_uadd(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

            swp = pt2.X;
            pt2.X = pt4.X;
            pt4.X = pt3.X;
            pt3.X = swp;
            swp = pt2.Z;
            pt2.Z = pt4.Z;
            pt4.Z = pt3.Z;
            pt3.Z = swp;
        }
        else if (steps[i] == 5)
        {
            d2 = submod(pt1.X, pt1.Z, n);
            s2 = addmod(pt1.X, pt1.Z, n);

            uecm_uadd(rho, n, pt2, pt1, pt3, &pt2);        // B = B + A (C)
            uecm_udup(s, rho, n, s2, d2, &pt1);        // A = 2A
        }
        else if (steps[i] == 6)
        {
            uecm_uadd(rho, n, pt1, pt2, pt3, P);     // A = A + B (C)
        }
    }
    return;
}

static void uecm_uprac85(uint64_t rho, uint64_t n, uecm_pt *P, uint64_t s)
{
    uint64_t s1, s2, d1, d2;
    uint64_t swp;
    int i;
    static const uint8_t steps[146] = {
        0,6,0,6,0,6,0,6,0,4,
        6,0,4,6,0,4,4,6,0,4,
        4,6,0,5,4,6,0,3,3,4,
        6,0,3,5,4,6,0,3,4,3,
        4,6,0,5,5,4,6,0,5,3,
        3,4,6,0,3,3,4,3,4,6,
        0,4,3,4,3,5,3,3,3,3,
        3,3,3,3,4,6,0,3,3,3,
        3,3,3,3,3,3,4,3,4,3,
        4,6,0,3,4,3,5,4,6,0,
        5,3,3,3,4,6,0,5,4,3,
        5,4,6,0,4,3,3,3,5,4,
        6,0,4,3,5,3,3,4,6,0,
        3,3,3,3,5,4,6,0,3,3,
        3,4,3,3,4,6 };

    uecm_pt pt1, pt2, pt3;
    for (i = 0; i < 146; i++)
    {
        if (steps[i] == 0)
        {
            pt1.X = pt2.X = pt3.X = P->X;
            pt1.Z = pt2.Z = pt3.Z = P->Z;

            d1 = submod(pt1.X, pt1.Z, n);
            s1 = addmod(pt1.X, pt1.Z, n);
            uecm_udup(s, rho, n, s1, d1, &pt1);
        }
        else if (steps[i] == 3)
        {
            // integrate step 4 followed by swap(1,2)
            uecm_pt pt4;
            uecm_uadd(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

            swp = pt1.X;
            pt1.X = pt4.X;
            pt4.X = pt3.X;
            pt3.X = pt2.X;
            pt2.X = swp;
            swp = pt1.Z;
            pt1.Z = pt4.Z;
            pt4.Z = pt3.Z;
            pt3.Z = pt2.Z;
            pt2.Z = swp;
        }
        else if (steps[i] == 4)
        {
            uecm_pt pt4;
            uecm_uadd(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

            swp = pt2.X;
            pt2.X = pt4.X;
            pt4.X = pt3.X;
            pt3.X = swp;
            swp = pt2.Z;
            pt2.Z = pt4.Z;
            pt4.Z = pt3.Z;
            pt3.Z = swp;
        }
        else if (steps[i] == 5)
        {
            d2 = submod(pt1.X, pt1.Z, n);
            s2 = addmod(pt1.X, pt1.Z, n);

            uecm_uadd(rho, n, pt2, pt1, pt3, &pt2);        // B = B + A (C)
            uecm_udup(s, rho, n, s2, d2, &pt1);        // A = 2A
        }
        else if (steps[i] == 6)
        {
            uecm_uadd(rho, n, pt1, pt2, pt3, P);     // A = A + B (C)
        }
    }
    return;
}

static void uecm_uprac(uint64_t rho, uint64_t n, uecm_pt *P, uint64_t c, double v, uint64_t s)
{
    uint64_t d, e, r;
    int i;
    uint64_t s1, s2, d1, d2;
    uint64_t swp;

    // we require c != 0
    int shift = _trail_zcnt64(c);
    c = c >> shift;

    d = c;
    r = (uint64_t)((double)d * v + 0.5);

    d = c - r;
    e = 2 * r - c;

    uecm_pt pt1, pt2, pt3;

    // the first one is always a doubling
    // point1 is [1]P
    pt1.X = pt2.X = pt3.X = P->X;
    pt1.Z = pt2.Z = pt3.Z = P->Z;

    d1 = submod(pt1.X, pt1.Z, n);
    s1 = addmod(pt1.X, pt1.Z, n);

    // point2 is [2]P
    uecm_udup(s, rho, n, s1, d1, &pt1);

    while (d != e)
    {
        if (d < e)
        {
            r = d;
            d = e;
            e = r;
            swp = pt1.X;
            pt1.X = pt2.X;
            pt2.X = swp;
            swp = pt1.Z;
            pt1.Z = pt2.Z;
            pt2.Z = swp;
        }
        if (d - e <= e / 4 && ((d + e) % 3) == 0)
        {
            d = (2 * d - e) / 3;
            e = (e - d) / 2;

            uecm_pt pt4;
            uecm_uadd(rho, n, pt1, pt2, pt3, &pt4); // T = A + B (C)
            uecm_pt pt5;
            uecm_uadd(rho, n, pt4, pt1, pt2, &pt5); // T2 = T + A (B)
            uecm_uadd(rho, n, pt2, pt4, pt1, &pt2); // B = B + T (A)

            swp = pt1.X;
            pt1.X = pt5.X;
            pt5.X = swp;
            swp = pt1.Z;
            pt1.Z = pt5.Z;
            pt5.Z = swp;
        }
        else if (d - e <= e / 4 && (d - e) % 6 == 0)
        {
            d = (d - e) / 2;

            d1 = submod(pt1.X, pt1.Z, n);
            s1 = addmod(pt1.X, pt1.Z, n);

            uecm_uadd(rho, n, pt1, pt2, pt3, &pt2);        // B = A + B (C)
            uecm_udup(s, rho, n, s1, d1, &pt1);        // A = 2A
        }
        else if ((d + 3) / 4 <= e)
        {
            d -= e;

            uecm_pt pt4;
            uecm_uadd(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

            swp = pt2.X;
            pt2.X = pt4.X;
            pt4.X = pt3.X;
            pt3.X = swp;
            swp = pt2.Z;
            pt2.Z = pt4.Z;
            pt4.Z = pt3.Z;
            pt3.Z = swp;
        }
        else if ((d + e) % 2 == 0)
        {
            d = (d - e) / 2;

            d2 = submod(pt1.X, pt1.Z, n);
            s2 = addmod(pt1.X, pt1.Z, n);

            uecm_uadd(rho, n, pt2, pt1, pt3, &pt2);        // B = B + A (C)
            uecm_udup(s, rho, n, s2, d2, &pt1);        // A = 2A
        }
        else if (d % 2 == 0)
        {
            d /= 2;

            d2 = submod(pt1.X, pt1.Z, n);
            s2 = addmod(pt1.X, pt1.Z, n);

            uecm_uadd(rho, n, pt3, pt1, pt2, &pt3);        // C = C + A (B)
            uecm_udup(s, rho, n, s2, d2, &pt1);        // A = 2A
        }
        else if (d % 3 == 0)
        {
            d = d / 3 - e;

            d1 = submod(pt1.X, pt1.Z, n);
            s1 = addmod(pt1.X, pt1.Z, n);

            uecm_pt pt4;
            uecm_udup(s, rho, n, s1, d1, &pt4);        // T = 2A
            uecm_pt pt5;
            uecm_uadd(rho, n, pt1, pt2, pt3, &pt5);        // T2 = A + B (C)
            uecm_uadd(rho, n, pt4, pt1, pt1, &pt1);        // A = T + A (A)
            uecm_uadd(rho, n, pt4, pt5, pt3, &pt4);        // T = T + T2 (C)

            swp = pt3.X;
            pt3.X = pt2.X;
            pt2.X = pt4.X;
            pt4.X = swp;
            swp = pt3.Z;
            pt3.Z = pt2.Z;
            pt2.Z = pt4.Z;
            pt4.Z = swp;

        }
        else if ((d + e) % 3 == 0)
        {
            d = (d - 2 * e) / 3;

            uecm_pt pt4;
            uecm_uadd(rho, n, pt1, pt2, pt3, &pt4);        // T = A + B (C)


            d2 = submod(pt1.X, pt1.Z, n);
            s2 = addmod(pt1.X, pt1.Z, n);
            uecm_uadd(rho, n, pt4, pt1, pt2, &pt2);        // B = T + A (B)
            uecm_udup(s, rho, n, s2, d2, &pt4);        // T = 2A
            uecm_uadd(rho, n, pt1, pt4, pt1, &pt1);        // A = A + T (A) = 3A
        }
        else if ((d - e) % 3 == 0)
        {
            d = (d - e) / 3;

            uecm_pt pt4;
            uecm_uadd(rho, n, pt1, pt2, pt3, &pt4);        // T = A + B (C)

            d2 = submod(pt1.X, pt1.Z, n);
            s2 = addmod(pt1.X, pt1.Z, n);
            uecm_uadd(rho, n, pt3, pt1, pt2, &pt3);        // C = C + A (B)

            swp = pt2.X;
            pt2.X = pt4.X;
            pt4.X = swp;
            swp = pt2.Z;
            pt2.Z = pt4.Z;
            pt4.Z = swp;

            uecm_udup(s, rho, n, s2, d2, &pt4);        // T = 2A
            uecm_uadd(rho, n, pt1, pt4, pt1, &pt1);        // A = A + T (A) = 3A
        }
        else
        {
            e /= 2;

            d2 = submod(pt2.X, pt2.Z, n);
            s2 = addmod(pt2.X, pt2.Z, n);

            uecm_uadd(rho, n, pt3, pt2, pt1, &pt3);        // C = C + B (A)
            uecm_udup(s, rho, n, s2, d2, &pt2);        // B = 2B
        }
    }
    uecm_uadd(rho, n, pt1, pt2, pt3, P);     // A = A + B (C)

    for (i = 0; i < shift; i++)
    {
        d1 = submod(P->X, P->Z, n);
        s1 = addmod(P->X, P->Z, n);
        uecm_udup(s, rho, n, s1, d1, P);     // P = 2P
    }
    return;
}


// jeff: "likely_gcd" is probably always the correct gcd, but I didn't add this
// parameter by using any proof; it's conceivable it might be wrong sometimes.
// See comments within the function.
inline static uint64_t uecm_modinv_64(uint64_t a, uint64_t p, uint64_t* plikely_gcd)
{
    /* thanks to the folks at www.mersenneforum.org */

    uint64_t ps1, ps2, parity, dividend, divisor, rem, q, t;

    q = 1;
    rem = a;
    dividend = p;
    divisor = a;
    ps1 = 1;
    ps2 = 0;
    parity = 0;

    while (divisor > 1) {
        rem = dividend - divisor;
        t = rem - divisor;
        if (rem >= divisor) {
            q += ps1; rem = t; t -= divisor;
            if (rem >= divisor) {
                q += ps1; rem = t; t -= divisor;
                if (rem >= divisor) {
                    q += ps1; rem = t; t -= divisor;
                    if (rem >= divisor) {
                        q += ps1; rem = t; t -= divisor;
                        if (rem >= divisor) {
                            q += ps1; rem = t; t -= divisor;
                            if (rem >= divisor) {
                                q += ps1; rem = t; t -= divisor;
                                if (rem >= divisor) {
                                    q += ps1; rem = t; t -= divisor;
                                    if (rem >= divisor) {
                                        q += ps1; rem = t;
                                        if (rem >= divisor) {
                                            q = dividend / divisor;
                                            rem = dividend % divisor;
                                            q *= ps1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        q += ps2;
        parity = ~parity;
        dividend = divisor;
        divisor = rem;
        ps2 = ps1;
        ps1 = q;
    }

    // jeff: added "likely_gcd".  this function seems to be a variant on the
    // extended euclidean algorithm, and thus the gcd likely equals the dividend
    // as calculated below.  However I'm doing this by analogy and educated
    // guess, not by proof.  It appears to work in all tests, so I suspect it is
    // correct, but it's possible this could be wrong in some cases.
    if (divisor == 1)
        dividend = divisor;
    *plikely_gcd = dividend;


    if (parity == 0)
        return ps1;
    else
        return p - ps1;
}

static uint64_t uecm_build(uecm_pt *P, uint64_t rho, uint64_t n, uint64_t *ploc_lcg, uint64_t* ps, uint64_t five, uint64_t Rsqr)
{
    uint64_t t1, t2, t3, t4;
    uint64_t u, v;

    uint32_t sigma = uecm_lcg_rand_32B(7, (uint32_t)-1, ploc_lcg);

    u = uecm_mulredc((uint64_t)sigma, Rsqr, n, rho);  // to_monty(sigma)

    //printf("sigma = %" PRIu64 ", u = %" PRIu64 ", n = %" PRIu64 "\n", sigma, u, n);

    v = addmod(u, u, n);
    v = addmod(v, v, n);            // 4*sigma

    //printf("v = %" PRIu64 "\n", v);

    u = uecm_sqrredc(u, n, rho);
    t1 = five;

    //printf("monty(5) = %" PRIu64 "\n", t1);

    u = submod(u, t1, n);           // sigma^2 - 5

    //printf("u = %" PRIu64 "\n", u);

    t1 = uecm_mulredc(u, u, n, rho);
    uint64_t tmpx = uecm_mulredc(t1, u, n, rho);  // u^3

    uint64_t v2 = addmod(v, v, n);             // 2*v
    uint64_t v4 = addmod(v2, v2, n);           // 4*v
    uint64_t v8 = addmod(v4, v4, n);           // 8*v
    uint64_t v16 = addmod(v8, v8, n);          // 16*v
    uint64_t t5 = uecm_mulredc(v16, tmpx, n, rho);    // 16*u^3*v

    t1 = uecm_mulredc(v, v, n, rho);
    uint64_t tmpz = uecm_mulredc(t1, v, n, rho);  // v^3

    //compute parameter A
    t1 = submod(v, u, n);           // (v - u)
    t2 = uecm_sqrredc(t1, n, rho);
    t4 = uecm_mulredc(t2, t1, n, rho);   // (v - u)^3

    t1 = addmod(u, u, n);           // 2u
    t2 = addmod(u, v, n);           // u + v
    t3 = addmod(t1, t2, n);         // 3u + v

    t1 = uecm_mulredc(t3, t4, n, rho);   // a = (v-u)^3 * (3u + v)

    // u holds the denom (jeff note: isn't it t5 that has the denom?)
    // t1 holds the numer
    // accomplish the division by multiplying by the modular inverse
    t2 = 1;
    t5 = uecm_mulredc(t5, t2, n, rho);   // take t5 out of monty rep

    uint64_t likely_gcd;
    t3 = uecm_modinv_64(t5, n, &likely_gcd);

    t3 = uecm_mulredc(t3, Rsqr, n, rho); // to_monty(t3)
    *ps = uecm_mulredc(t3, t1, n, rho);

    P->X = tmpx;
    P->Z = tmpz;

    return likely_gcd;
}

#ifdef USE_AVX512F


static __m512i uecm_build_x8(uecm_pt_x8* P, __m512i rho, __m512i n, uint64_t* ploc_lcg,
    uint64_t* likely_gcd, __m512i five, __m512i Rsqr)
{
    __m512i t1, t2, t3, t4;
    __m512i u, v, s;
    int i;
    uint64_t s8[8];
    uint32_t sigma = uecm_lcg_rand_32B(7, (uint32_t)-1, ploc_lcg);

    for (i = 0; i < 8; i++)
    {
        s8[i] = sigma;
        sigma = uecm_lcg_rand_32B(7, (uint32_t)-1, ploc_lcg);
    }
    s = _mm512_loadu_epi64(s8);

    u = uecm_mulredc_x8(s, Rsqr, n, rho);  // to_monty(sigma)

    //printf("sigma = %" PRIu64 ", u = %" PRIu64 ", n = %" PRIu64 "\n", sigma, u, n);

    v = uecm_addmod_x8(u, u, n);
    v = uecm_addmod_x8(v, v, n);            // 4*sigma

    //printf("v = %" PRIu64 "\n", v);

    u = uecm_sqrredc_x8(u, n, rho);
    t1 = five;

    //printf("monty(5) = %" PRIu64 "\n", t1);

    u = uecm_submod_x8(u, t1, n);           // sigma^2 - 5

    //printf("u = %" PRIu64 "\n", u);

    t1 = uecm_mulredc_x8(u, u, n, rho);
    __m512i tmpx = uecm_mulredc_x8(t1, u, n, rho);  // u^3

    __m512i v2 = uecm_addmod_x8(v, v, n);             // 2*v
    __m512i v4 = uecm_addmod_x8(v2, v2, n);           // 4*v
    __m512i v8 = uecm_addmod_x8(v4, v4, n);           // 8*v
    __m512i v16 = uecm_addmod_x8(v8, v8, n);          // 16*v
    __m512i t5 = uecm_mulredc_x8(v16, tmpx, n, rho);    // 16*u^3*v

    t1 = uecm_mulredc_x8(v, v, n, rho);
    __m512i tmpz = uecm_mulredc_x8(t1, v, n, rho);  // v^3

    //compute parameter A
    t1 = uecm_submod_x8(v, u, n);           // (v - u)
    t2 = uecm_sqrredc_x8(t1, n, rho);
    t4 = uecm_mulredc_x8(t2, t1, n, rho);   // (v - u)^3

    t1 = uecm_addmod_x8(u, u, n);           // 2u
    t2 = uecm_addmod_x8(u, v, n);           // u + v
    t3 = uecm_addmod_x8(t1, t2, n);         // 3u + v

    t1 = uecm_mulredc_x8(t3, t4, n, rho);   // a = (v-u)^3 * (3u + v)

    // u holds the denom (jeff note: isn't it t5 that has the denom?)
    // t1 holds the numer
    // accomplish the division by multiplying by the modular inverse
    t2 = uecm_set64(1);
    t5 = uecm_mulredc_x8(t5, t2, n, rho);   // take t5 out of monty rep

    uint64_t t3_64[8];
    uint64_t t5_64[8];
    uint64_t n_64[8];
    _mm512_storeu_epi64(t5_64, t5);
    _mm512_storeu_epi64(n_64, n);
    for (i = 0; i < 8; i++)
    {
        t3_64[i] = uecm_modinv_64(t5_64[i], n_64[i], &likely_gcd[i]);
    }

    t3 = _mm512_loadu_epi64(t3_64);

    t3 = uecm_mulredc_x8(t3, Rsqr, n, rho); // to_monty(t3)
    __m512i ps = uecm_mulredc_x8(t3, t1, n, rho);

    P->X = tmpx;
    P->Z = tmpz;

    return ps;
}

static uint64_t uecm_build52(uecm_pt* P, uint64_t rho, uint64_t n, uint64_t* ploc_lcg, uint64_t* ps, uint64_t five, uint64_t Rsqr)
{
    uint64_t t1, t2, t3, t4;
    uint64_t u, v;

    uint32_t sigma = uecm_lcg_rand_32B(7, (uint32_t)-1, ploc_lcg);

    u = uecm_mulredc52((uint64_t)sigma, Rsqr, n, rho);  // to_monty(sigma)

    //printf("sigma = %" PRIu64 ", u = %" PRIu64 ", n = %" PRIu64 "\n", sigma, u, n);

    v = addmod52(u, u, n);
    v = addmod52(v, v, n);            // 4*sigma

    //printf("v = %" PRIu64 "\n", v);

    u = uecm_sqrredc52(u, n, rho);
    t1 = five;

    //printf("monty(5) = %" PRIu64 "\n", t1);

    u = submod52(u, t1, n);           // sigma^2 - 5

    //printf("u = %" PRIu64 "\n", u);

    t1 = uecm_mulredc52(u, u, n, rho);
    uint64_t tmpx = uecm_mulredc52(t1, u, n, rho);  // u^3

    uint64_t v2 = addmod52(v, v, n);             // 2*v
    uint64_t v4 = addmod52(v2, v2, n);           // 4*v
    uint64_t v8 = addmod52(v4, v4, n);           // 8*v
    uint64_t v16 = addmod52(v8, v8, n);          // 16*v
    uint64_t t5 = uecm_mulredc52(v16, tmpx, n, rho);    // 16*u^3*v

    t1 = uecm_mulredc52(v, v, n, rho);
    uint64_t tmpz = uecm_mulredc52(t1, v, n, rho);  // v^3

    //compute parameter A
    t1 = submod52(v, u, n);           // (v - u)
    t2 = uecm_sqrredc52(t1, n, rho);
    t4 = uecm_mulredc52(t2, t1, n, rho);   // (v - u)^3

    t1 = addmod52(u, u, n);           // 2u
    t2 = addmod52(u, v, n);           // u + v
    t3 = addmod52(t1, t2, n);         // 3u + v

    t1 = uecm_mulredc52(t3, t4, n, rho);   // a = (v-u)^3 * (3u + v)

    // u holds the denom (jeff note: isn't it t5 that has the denom?)
    // t1 holds the numer
    // accomplish the division by multiplying by the modular inverse
    t2 = 1;
    t5 = uecm_mulredc52(t5, t2, n, rho);   // take t5 out of monty rep

    uint64_t likely_gcd;
    t3 = uecm_modinv_64(t5, n, &likely_gcd);

    t3 = uecm_mulredc52(t3, Rsqr, n, rho); // to_monty(t3)
    *ps = uecm_mulredc52(t3, t1, n, rho);

    P->X = tmpx;
    P->Z = tmpz;

    return likely_gcd;
}

#endif

static int uecm_check_factor(uint64_t Z, uint64_t n, uint64_t* f)
{
    int status = 0;
    *f = bingcd64(n, Z);

    if (*f > 1)
    {
        if (*f == n)
        {
            *f = 0;
            status = 0;
        }
        else
        {
            status = 1;
        }
    }
    return status;
}

static void uecm_stage1(uint64_t rho, uint64_t n, uecm_pt *P, uint64_t stg1, uint64_t s)
{
    uint64_t q;

    // handle the only even case
    q = 2;
    while (q < stg1 * 4)  // jeff: multiplying by 4 improves perf ~1%
    {
        uint64_t diff1 = submod(P->X, P->Z, n);
        uint64_t sum1 = addmod(P->X, P->Z, n);
        uecm_udup(s, rho, n, sum1, diff1, P);
        q *= 2;
    }

    if (stg1 == 27)
    {
        uecm_uprac(rho, n, P, 3, 0.61803398874989485, s);
        uecm_uprac(rho, n, P, 3, 0.61803398874989485, s);
        uecm_uprac(rho, n, P, 3, 0.61803398874989485, s);
        uecm_uprac(rho, n, P, 5, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 5, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 7, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 11, 0.580178728295464130, s);
        uecm_uprac(rho, n, P, 13, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 17, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 19, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 23, 0.522786351415446049, s);
    }
    else if (stg1 == 47)
    {
        // jeff: improved perf slightly by using one more uprac for 3,
        // and removing uprac for 47.
        uecm_uprac(rho, n, P, 3, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 3, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 3, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 3, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 5, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 5, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 7, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 11, 0.580178728295464130, s);
        uecm_uprac(rho, n, P, 13, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 17, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 19, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 23, 0.522786351415446049, s);
        uecm_uprac(rho, n, P, 29, 0.548409048446403258, s);
        uecm_uprac(rho, n, P, 31, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 37, 0.580178728295464130, s);
        uecm_uprac(rho, n, P, 41, 0.548409048446403258, s);
        uecm_uprac(rho, n, P, 43, 0.618033988749894903, s);
//        uecm_uprac(rho, n, P, 47, 0.548409048446403258, s);
    }
    else if (stg1 == 59)
    {   // jeff: probably stg1 of 59 would benefit from similar changes
        // as stg1 of 47 above, but I didn't bother. Stg1 of 59 seems to
        // always perform worse than stg1 of 47, so there doesn't seem
        // to be any reason to ever use stg1 of 59.
        uecm_uprac(rho, n, P, 3, 0.61803398874989485, s);
        uecm_uprac(rho, n, P, 3, 0.61803398874989485, s);
        uecm_uprac(rho, n, P, 3, 0.61803398874989485, s);
        uecm_uprac(rho, n, P, 5, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 5, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 7, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 7, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 11, 0.580178728295464130, s);
        uecm_uprac(rho, n, P, 13, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 17, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 19, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 23, 0.522786351415446049, s);
        uecm_uprac(rho, n, P, 29, 0.548409048446403258, s);
        uecm_uprac(rho, n, P, 31, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 1961, 0.552936068843375, s);   // 37 * 53
        uecm_uprac(rho, n, P, 41, 0.548409048446403258, s);
        uecm_uprac(rho, n, P, 43, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 47, 0.548409048446403258, s);
        uecm_uprac(rho, n, P, 59, 0.548409048446403258, s);
    }
    else if (stg1 == 70)
    {
        // call prac with best ratio found in deep search.
        // some composites are cheaper than their
        // constituent primes.
        uecm_uprac70(rho, n, P, s);
    }
    else // if (stg1 >= 85)
    {
        uecm_uprac85(rho, n, P, s);

        if (stg1 == 85)
        {
            uecm_uprac(rho, n, P, 61, 0.522786351415446049, s);
        }
        else
        {
            uecm_uprac(rho, n, P, 5, 0.618033988749894903, s);
            uecm_uprac(rho, n, P, 11, 0.580178728295464130, s);
//            uecm_uprac(rho, n, P, 61, 0.522786351415446049, s);
            uecm_uprac(rho, n, P, 89, 0.618033988749894903, s);
            uecm_uprac(rho, n, P, 97, 0.723606797749978936, s);
            uecm_uprac(rho, n, P, 101, 0.556250337855490828, s);
            uecm_uprac(rho, n, P, 107, 0.580178728295464130, s);
            uecm_uprac(rho, n, P, 109, 0.548409048446403258, s);
            uecm_uprac(rho, n, P, 113, 0.618033988749894903, s);

            if (stg1 == 125)
            {
                // jeff: moved 61 to here
                uecm_uprac(rho, n, P, 61, 0.522786351415446049, s);
                uecm_uprac(rho, n, P, 103, 0.632839806088706269, s);
            }
            else
            {
                uecm_uprac(rho, n, P, 7747, 0.552188778811121, s); // 61 x 127
                uecm_uprac(rho, n, P, 131, 0.618033988749894903, s);
                uecm_uprac(rho, n, P, 14111, 0.632839806088706, s);  // 103 x 137
                uecm_uprac(rho, n, P, 20989, 0.620181980807415, s);  // 139 x 151
                uecm_uprac(rho, n, P, 157, 0.640157392785047019, s);
                uecm_uprac(rho, n, P, 163, 0.551390822543526449, s);

                if (stg1 == 165)
                {
                    uecm_uprac(rho, n, P, 149, 0.580178728295464130, s);
                }
                else
                {
                    uecm_uprac(rho, n, P, 13, 0.618033988749894903, s);
                    uecm_uprac(rho, n, P, 167, 0.580178728295464130, s);
                    uecm_uprac(rho, n, P, 173, 0.612429949509495031, s);
                    uecm_uprac(rho, n, P, 179, 0.618033988749894903, s);
                    uecm_uprac(rho, n, P, 181, 0.551390822543526449, s);
                    uecm_uprac(rho, n, P, 191, 0.618033988749894903, s);
                    uecm_uprac(rho, n, P, 193, 0.618033988749894903, s);
                    uecm_uprac(rho, n, P, 29353, 0.580178728295464, s);  // 149 x 197
                    uecm_uprac(rho, n, P, 199, 0.551390822543526449, s);
                }
            }
        }
    }
    return;
}

#ifdef USE_AVX512F

static void uecm_stage1_x8(__m512i rho, __m512i n, uecm_pt_x8* P, uint64_t stg1, __m512i s)
{
    uint64_t q;

    // handle the only even case
    q = 2;
    while (q < stg1 * 4)  // jeff: multiplying by 4 improves perf ~1%
    {
        __m512i diff1 = uecm_submod_x8(P->X, P->Z, n);
        __m512i sum1 = uecm_addmod_x8(P->X, P->Z, n);
        uecm_udup_x8(s, rho, n, sum1, diff1, P);
        q *= 2;
    }

    if (stg1 == 27)
    {
        uecm_uprac_x8(rho, n, P, 3, 0.61803398874989485, s);
        uecm_uprac_x8(rho, n, P, 3, 0.61803398874989485, s);
        uecm_uprac_x8(rho, n, P, 3, 0.61803398874989485, s);
        uecm_uprac_x8(rho, n, P, 5, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 5, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 7, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 11, 0.580178728295464130, s);
        uecm_uprac_x8(rho, n, P, 13, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 17, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 19, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 23, 0.522786351415446049, s);
    }
    else if (stg1 == 47)
    {
        // jeff: improved perf slightly by using one more uprac for 3,
        // and removing uprac for 47.
        uecm_uprac_x8(rho, n, P, 3, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 3, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 3, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 3, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 5, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 5, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 7, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 11, 0.580178728295464130, s);
        uecm_uprac_x8(rho, n, P, 13, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 17, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 19, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 23, 0.522786351415446049, s);
        uecm_uprac_x8(rho, n, P, 29, 0.548409048446403258, s);
        uecm_uprac_x8(rho, n, P, 31, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 37, 0.580178728295464130, s);
        uecm_uprac_x8(rho, n, P, 41, 0.548409048446403258, s);
        uecm_uprac_x8(rho, n, P, 43, 0.618033988749894903, s);
        //        uecm_uprac(rho, n, P, 47, 0.548409048446403258, s);
    }
    else if (stg1 == 59)
    {   // jeff: probably stg1 of 59 would benefit from similar changes
        // as stg1 of 47 above, but I didn't bother. Stg1 of 59 seems to
        // always perform worse than stg1 of 47, so there doesn't seem
        // to be any reason to ever use stg1 of 59.
        uecm_uprac_x8(rho, n, P, 3, 0.61803398874989485, s);
        uecm_uprac_x8(rho, n, P, 3, 0.61803398874989485, s);
        uecm_uprac_x8(rho, n, P, 3, 0.61803398874989485, s);
        uecm_uprac_x8(rho, n, P, 5, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 5, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 7, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 7, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 11, 0.580178728295464130, s);
        uecm_uprac_x8(rho, n, P, 13, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 17, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 19, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 23, 0.522786351415446049, s);
        uecm_uprac_x8(rho, n, P, 29, 0.548409048446403258, s);
        uecm_uprac_x8(rho, n, P, 31, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 1961, 0.552936068843375, s);   // 37 * 53
        uecm_uprac_x8(rho, n, P, 41, 0.548409048446403258, s);
        uecm_uprac_x8(rho, n, P, 43, 0.618033988749894903, s);
        uecm_uprac_x8(rho, n, P, 47, 0.548409048446403258, s);
        uecm_uprac_x8(rho, n, P, 59, 0.548409048446403258, s);
    }
    else if (stg1 == 70)
    {
        // call prac with best ratio found in deep search.
        // some composites are cheaper than their
        // constituent primes.
        uecm_uprac70_x8(rho, n, P, s);
    }
    else // if (stg1 >= 85)
    {
        uecm_uprac85_x8(rho, n, P, s);

        if (stg1 == 85)
        {
            uecm_uprac_x8(rho, n, P, 61, 0.522786351415446049, s);
        }
        else
        {
            uecm_uprac_x8(rho, n, P, 5, 0.618033988749894903, s);
            uecm_uprac_x8(rho, n, P, 11, 0.580178728295464130, s);
            //            uecm_uprac(rho, n, P, 61, 0.522786351415446049, s);
            uecm_uprac_x8(rho, n, P, 89, 0.618033988749894903, s);
            uecm_uprac_x8(rho, n, P, 97, 0.723606797749978936, s);
            uecm_uprac_x8(rho, n, P, 101, 0.556250337855490828, s);
            uecm_uprac_x8(rho, n, P, 107, 0.580178728295464130, s);
            uecm_uprac_x8(rho, n, P, 109, 0.548409048446403258, s);
            uecm_uprac_x8(rho, n, P, 113, 0.618033988749894903, s);

            if (stg1 == 125)
            {
                // jeff: moved 61 to here
                uecm_uprac_x8(rho, n, P, 61, 0.522786351415446049, s);
                uecm_uprac_x8(rho, n, P, 103, 0.632839806088706269, s);
            }
            else
            {
                uecm_uprac_x8(rho, n, P, 7747, 0.552188778811121, s); // 61 x 127
                uecm_uprac_x8(rho, n, P, 131, 0.618033988749894903, s);
                uecm_uprac_x8(rho, n, P, 14111, 0.632839806088706, s);  // 103 x 137
                uecm_uprac_x8(rho, n, P, 20989, 0.620181980807415, s);  // 139 x 151
                uecm_uprac_x8(rho, n, P, 157, 0.640157392785047019, s);
                uecm_uprac_x8(rho, n, P, 163, 0.551390822543526449, s);

                if (stg1 == 165)
                {
                    uecm_uprac_x8(rho, n, P, 149, 0.580178728295464130, s);
                }
                else
                {
                    uecm_uprac_x8(rho, n, P, 13, 0.618033988749894903, s);
                    uecm_uprac_x8(rho, n, P, 167, 0.580178728295464130, s);
                    uecm_uprac_x8(rho, n, P, 173, 0.612429949509495031, s);
                    uecm_uprac_x8(rho, n, P, 179, 0.618033988749894903, s);
                    uecm_uprac_x8(rho, n, P, 181, 0.551390822543526449, s);
                    uecm_uprac_x8(rho, n, P, 191, 0.618033988749894903, s);
                    uecm_uprac_x8(rho, n, P, 193, 0.618033988749894903, s);
                    uecm_uprac_x8(rho, n, P, 29353, 0.580178728295464, s);  // 149 x 197
                    uecm_uprac_x8(rho, n, P, 199, 0.551390822543526449, s);
                }
            }
        }
    }
    return;
}

#endif

#ifdef MICRO_ECM_VERBOSE_PRINTF
static void uecm_stage2_pair(int b1, int b2, 
    int Astart, int D, int AD, int *p)
{
    // a very simple pairing procedure.  
    int pid = 0;
    int A;
    int i;
    int* pairs;
    int steps[1000];
    int s = 0;

    pairs = (int*)malloc(D * sizeof(int));

    printf("commencing pair for A0=%d, D=%d, b1=%d, b2=%d\n",
        Astart, D, b1, b2);

    while (p[pid] < b1) {
        pid++;
    }

    printf("now at p=%d\n", p[pid]);

    for (i = 0; i < D; i++)
    {
        pairs[i] = 0;
    }

    A = Astart;

    while (p[pid] < b2)
    {
        int b = A - p[pid];

        if ((p[pid] - A) > (D / 2))
        {
            A += D;
            printf("new A = %d\n", A);
            for (i = 0; i < D / 2; i++)
            {
                pairs[i] = 0;
            }
            steps[s++] = 0;
            continue;
        }

        if (b > 0)
        {
            if (b > D / 2)
            {
                printf("unable to pair %d\n", primes[pid]);
                pid++;
                continue;
            }

            printf("new pair %d +/- %d (%d:%d)\n", A, b, A - b, A + b);
            pairs[b] = 1;
            steps[s++] = b;
        }
        else
        {
            if (pairs[abs(b)] == 0)
            {
                printf("leftover pair %d +/- %d (%d:%d)\n", A, b, A - b, A + b);
                steps[s++] = abs(b);
            }
            
        }

        pid++;
    }

    free(pairs);

    printf("steps[%d] = {", s);
    for (i = 0; i < s; i++)
    {
        if (i % 16 == 0) printf("\n");
        printf("%d,", steps[i]);
    }
    printf("\b};\n");
    return;
}
#endif



// pre-paired sequences for various B1 and B2 = 25*B1
static const int numb1_27 = 76;
static const uint8_t b1_27[76] = {
1,7,13,17,23,29,31,41,43,47,49,0,59,47,43,41,
37,31,29,19,13,7,1,11,23,0,59,53,43,41,37,31,
23,17,11,7,1,19,29,49,0,53,49,47,43,31,23,19,
11,7,1,13,37,59,0,59,53,43,37,31,29,23,17,13,
11,1,47,0,59,49,41,31,23,17,11,7 };

static const int numb1_47 = 121;
static const uint8_t b1_47[121] = {
1,7,13,17,23,29,31,41,43,47,49,0,59,47,43,41,
37,31,29,19,13,7,1,11,23,0,59,53,43,41,37,31,
23,17,11,7,1,19,29,49,0,53,49,47,43,31,23,19,
11,7,1,13,37,59,0,59,53,43,37,31,29,23,17,13,
11,1,47,0,59,49,41,31,23,17,11,7,1,19,37,47,
0,59,49,47,43,41,31,17,13,11,7,37,0,53,49,43,
37,23,19,13,7,1,29,31,41,59,0,59,49,47,41,23,
19,17,13,7,1,43,53,0,59 };

static const int numb1_70 = 175;
static const uint8_t b1_70[175] = {
31,41,43,47,49,0,59,47,43,41,37,31,29,19,13,7,
1,11,23,0,59,53,43,41,37,31,23,17,11,7,1,19,
29,49,0,53,49,47,43,31,23,19,11,7,1,13,37,59,
0,59,53,43,37,31,29,23,17,13,11,1,47,0,59,49,
41,31,23,17,11,7,1,19,37,47,0,59,49,47,43,41,
31,17,13,11,7,37,0,53,49,43,37,23,19,13,7,1,
29,31,41,59,0,59,49,47,41,23,19,17,13,7,1,43,
53,0,59,49,43,37,29,17,13,7,1,19,47,53,0,59,
53,49,47,43,31,29,23,11,17,0,47,43,41,37,31,23,
19,17,11,1,13,29,53,0,59,47,41,37,31,23,19,11,
7,17,29,0,53,47,43,41,17,13,11,1,23,31,37 };

static const int numb1_85 = 225;
static const uint8_t b1_85[225] = {
    1,53,49,47,43,41,37,23,19,13,11,1,7,17,29,31,0,59,47,43,41,37,31,29,19,13,7,1,11,23,0,59,53,43,41,37,
    31,23,17,11,7,1,19,29,49,0,53,49,47,43,31,23,19,11,7,1,13,37,59,0,59,53,43,37,31,29,23,17,13,11,1,47,
    0,59,49,41,31,23,17,11,7,1,19,37,47,0,59,49,47,43,41,31,17,13,11,7,37,0,53,49,43,37,23,19,13,7,1,29,
    31,41,59,0,59,49,47,41,23,19,17,13,7,1,43,53,0,59,49,43,37,29,17,13,7,1,19,47,53,0,59,53,49,47,43,31,
    29,23,11,17,0,47,43,41,37,31,23,19,17,11,1,13,29,53,0,59,47,41,37,31,23,19,11,7,17,29,0,53,47,43,41,
    17,13,11,1,23,31,37,49,0,53,47,43,41,29,19,7,1,17,31,37,49,59,0,49,43,37,19,17,1,23,29,47,53,0,59,53,
    43,41,31,17,7,1,11,13,19,29 };

static const int numb1_125 = 319;
static const uint8_t b1_125[319] = {
    23,19,13,11,1,7,17,29,31,0,59,47,43,41,37,31,29,19,13,7,1,11,23,0,59,53,43,41,37,31,23,17,11,7,1,19,
    29,49,0,53,49,47,43,31,23,19,11,7,1,13,37,59,0,59,53,43,37,31,29,23,17,13,11,1,47,0,59,49,41,31,23,
    17,11,7,1,19,37,47,0,59,49,47,43,41,31,17,13,11,7,37,0,53,49,43,37,23,19,13,7,1,29,31,41,59,0,59,49,
    47,41,23,19,17,13,7,1,43,53,0,59,49,43,37,29,17,13,7,1,19,47,53,0,59,53,49,47,43,31,29,23,11,17,0,47,
    43,41,37,31,23,19,17,11,1,13,29,53,0,59,47,41,37,31,23,19,11,7,17,29,0,53,47,43,41,17,13,11,1,23,31,
    37,49,0,53,47,43,41,29,19,7,1,17,31,37,49,59,0,49,43,37,19,17,1,23,29,47,53,0,59,53,43,41,31,17,7,1,
    11,13,19,29,0,59,53,49,47,37,29,11,13,17,23,31,0,59,43,41,37,29,23,17,13,1,31,47,0,59,53,49,47,41,37,
    31,19,13,7,11,17,29,43,0,47,29,19,11,7,1,41,43,59,0,53,49,37,23,13,11,7,1,17,19,29,41,43,59,0,59,49,
    41,37,23,13,1,7,11,29,43,47,53,0,59,53,49,31,23,13,7,1,17,29,43,47,0,59,31,29,19,11,7,37,49,53 };

static const int numb1_165 = 425;
static const uint8_t b1_165[425] = {
    13,7,1,11,19,47,59,0,59,49,43,37,31,29,23,19,17,7,11,13,47,53,0,53,47,41,37,31,23,19,11,1,13,29,43,
    59,0,53,49,41,37,31,19,17,1,7,23,29,47,59,0,59,53,47,43,41,29,19,17,13,7,1,23,31,49,0,53,47,41,37,29,
    23,19,11,7,17,31,43,49,59,0,47,43,41,37,23,19,17,13,7,11,29,53,0,53,49,43,37,29,23,11,7,1,13,19,31,41,
    0,53,49,47,43,37,31,23,17,11,13,41,0,59,47,43,37,31,29,23,11,1,17,19,41,0,59,53,19,13,7,1,29,43,47,49,
    0,53,49,47,41,29,19,17,13,11,7,1,23,31,43,59,0,53,49,41,37,23,19,13,11,7,1,17,43,47,0,47,43,41,31,19,
    17,7,1,13,37,49,0,59,49,37,29,13,1,7,11,17,19,41,47,53,0,49,47,31,29,7,1,13,17,19,23,37,59,0,47,37,31,
    19,17,13,11,1,29,41,43,53,0,59,41,17,13,7,1,19,23,31,47,49,53,0,59,53,47,43,31,29,7,1,11,17,37,41,49,
    0,49,43,37,23,19,13,1,7,17,0,59,49,41,37,31,29,23,1,11,13,53,0,53,43,41,37,29,23,17,13,11,7,1,19,31,49,
    0,53,43,31,29,23,19,17,1,13,37,41,59,0,53,43,37,31,23,13,1,17,29,59,0,59,49,41,37,23,19,11,1,7,29,0,59,
    43,17,13,11,1,7,23,29,37,41,49,0,49,47,43,41,29,1,7,13,19,23,31,59,0,59,49,47,31,29,13,7,37,41,43,0,49,
    41,29,23,13,11,7,1,17,19,31,43,53,0,53,47,43,37,29,23,17,1,11,13,31,41,49,59,0,53,47,41,19,13,11,1,17,
    23,43,0,53,49,47,37,23,19,11,7,17,29,31,43,0,53,31,19,17,13,7,1,29,37,59 };

static const int numb1_205 = 511;
static const uint8_t b1_205[511] = {
    1,23,41,0,59,53,49,47,37,23,19,17,13,1,7,29,43,0,53,49,41,31,29,19,17,11,7,1,13,37,59,0,49,47,29,23,
    13,7,1,17,31,37,43,0,59,49,47,43,37,31,29,17,13,7,1,11,19,53,0,59,53,49,41,37,23,13,1,11,17,19,29,43,
    47,0,53,49,47,43,23,19,11,1,7,17,37,41,0,59,53,41,37,31,29,19,17,11,1,13,43,47,0,53,47,41,19,17,7,1,
    11,23,31,43,59,0,59,53,41,31,13,11,7,1,17,29,37,0,49,43,37,29,11,1,13,17,19,23,41,0,59,49,47,43,41,37,
    31,19,7,1,13,23,29,53,0,53,49,43,41,37,31,29,23,13,7,17,19,47,59,0,49,47,37,29,23,17,11,7,13,19,31,41,
    53,0,59,43,29,23,19,17,13,11,1,41,0,59,37,31,23,17,13,11,7,1,19,29,43,53,0,49,47,43,41,31,19,17,1,7,11,
    13,23,0,47,43,37,29,13,11,7,1,17,19,23,31,59,0,59,37,31,29,23,19,13,1,7,11,41,47,53,0,53,49,43,31,23,
    17,13,41,59,0,59,53,31,19,17,1,7,11,23,37,47,49,0,59,53,47,43,41,37,31,23,19,17,11,1,0,59,53,49,47,31,
    17,13,7,1,11,29,37,0,53,43,31,17,13,7,1,29,41,49,0,53,49,41,29,23,11,7,1,19,31,47,0,47,43,41,29,23,19,
    7,1,11,49,0,59,31,29,23,17,11,7,1,13,41,43,0,59,43,37,17,1,7,11,13,19,41,49,0,59,53,43,41,37,31,29,23,
    13,11,1,47,0,59,53,47,31,19,17,13,1,7,11,29,37,43,49,0,49,43,41,31,17,13,7,11,23,37,53,0,53,49,41,23,
    19,13,11,7,1,17,37,59,0,49,47,43,37,31,29,23,1,7,41,0,59,43,41,37,31,17,13,11,7,47,49,0,59,49,47,37,31,
    29,19,17,7,1,0,53,47,37,19,13,1,11,31,41,0,49,47,37,23,17,13,11,7,19,31,53,0,59,53,47,29,13,11,7,1,23,
    41,0,49,47,41,37,19,11,13,17,23,29,31,43,0,59,29,19,13,1,41,43,47,53,0,59,53,43,41,37,23,17,11,7,1,13,
    29,49 };


static uint64_t uecm_stage2(uecm_pt *P, uint64_t rho, uint64_t n, uint32_t stg1_max, uint64_t s, uint64_t unityval)
{
    int b;
    int i, j, k;
    uecm_pt Pa1;
    uecm_pt *Pa = &Pa1;
    uecm_pt Pb[18];
    uecm_pt Pd1;
    uecm_pt *Pd = &Pd1;
    const uint8_t *barray = 0;
    int numb;

#ifdef MICRO_ECM_VERBOSE_PRINTF
    ptadds = 0;
    stg1Doub = 0;
    stg1Add = 0;
#endif

    // this function has been written for MICRO_ECM_PARAM_D of 60, so you
    // probably don't want to change it.
    const int MICRO_ECM_PARAM_D = 60;

    //stage 2 init
    //Q = P = result of stage 1
    //compute [d]Q for 0 < d <= MICRO_ECM_PARAM_D

    uint64_t Pbprod[18];

    // [1]Q
    Pb[1] = *P;
    Pbprod[1] = uecm_mulredc(Pb[1].X, Pb[1].Z, n, rho);

    // [2]Q
    uint64_t diff1 = submod(P->X, P->Z, n);
    uint64_t sum1 = addmod(P->X, P->Z, n);
    uecm_udup(s, rho, n, sum1, diff1, &Pb[2]);
//    Pbprod[2] = uecm_mulredc(Pb[2].X, Pb[2].Z, n, rho);    // never used

    /*
    Let D = MICRO_ECM_PARAM_D.

    D is small in tinyecm, so it is straightforward to just enumerate the needed
    points.  We can do it efficiently with two progressions mod 6.
    Pb[0] = scratch
    Pb[1] = [1]Q;
    Pb[2] = [2]Q;
    Pb[3] = [7]Q;   prog2
    Pb[4] = [11]Q;  prog1
    Pb[5] = [13]Q;  prog2
    Pb[6] = [17]Q;  prog1
    Pb[7] = [19]Q;  prog2
    Pb[8] = [23]Q;  prog1
    Pb[9] = [29]Q;  prog1
    Pb[10] = [30 == D]Q;   // shouldn't this be [31]Q?
    Pb[11] = [31]Q; prog2   // shouldn't this be [37]Q?
    Pb[12] = [37]Q; prog2   // shouldn't this be [41]Q?
    Pb[13] = [41]Q; prog1   // [43]Q?
    Pb[14] = [43]Q; prog2   // [47]Q?
    Pb[15] = [47]Q; prog1   // [49]Q?
    Pb[16] = [49]Q; prog2   // [53]Q?
    Pb[17] = [53]Q; prog1   // [59]Q?
    // Pb[18] = [59]Q; prog1   // [60]Q?   note: we can remove this line I believe.  Pb[18] is never set, and never used.  Therefore I changed the definition of Pb above to have only 18 elements.

    two progressions with total of 17 adds to get 15 values of Pb.
    6 + 5(1) -> 11 + 6(5) -> 17 + 6(11) -> 23 + 6(17) -> 29 + 6(23) -> 35 + 6(29) -> 41 + 6(35) -> 47 + 6(41) -> 53 + 6(47) -> 59
    6 + 1(5) -> 7  + 6(1) -> 13 + 6(7)  -> 19 + 6(13) -> 25 + 6(19) -> 31 + 6(25) -> 37 + 6(31) -> 43 + 6(37) -> 49

    we also need [2D]Q = [60]Q
    to get [60]Q we just need one more add:
    compute [60]Q from [31]Q + [29]Q([2]Q), all of which we
    have after the above progressions are computed.

    we also need [A]Q = [((B1 + D) / (2D) * 2 D]Q
    which is equal to the following for various common B1:
    B1      [x]Q
    65      [120]Q
    85      [120]Q
    125     [180]Q      note: according to the A[Q] formula above, wouldn't this be [120]Q?  ( I suspect maybe the formula is supposed to be [((B1 + D) / D) * D]Q )
    165     [180]Q      note: same as above.
    205     [240]Q

    and we need [x-D]Q as well, for the above [x]Q.
    So far we are getting [x]Q and [x-2D]Q each from prac(x,Q).
    There is a better way using progressions of [2D]Q
    [120]Q = 2*[60]Q
    [180]Q = [120]Q + [60]Q([60]Q)
    [240]Q = 2*[120]Q
    [300]Q = [240]Q + [60]Q([180]Q)
    ...
    etc.

    */

    uecm_pt pt5, pt6;

    // Calculate all Pb: the following is specialized for MICRO_ECM_PARAM_D=60
    // [2]Q + [1]Q([1]Q) = [3]Q
    uecm_uadd(rho, n, Pb[1], Pb[2], Pb[1], &Pb[3]);        // <-- temporary

    // 2*[3]Q = [6]Q
    diff1 = submod(Pb[3].X, Pb[3].Z, n);
    sum1 = addmod(Pb[3].X, Pb[3].Z, n);
    uecm_udup(s, rho, n, sum1, diff1, &pt6);   // pt6 = [6]Q

    // [3]Q + [2]Q([1]Q) = [5]Q
    uecm_uadd(rho, n, Pb[3], Pb[2], Pb[1], &pt5);    // <-- pt5 = [5]Q
    Pb[3] = pt5;

    // [6]Q + [5]Q([1]Q) = [11]Q
    uecm_uadd(rho, n, pt6, pt5, Pb[1], &Pb[4]);    // <-- [11]Q

    i = 3;
    k = 4;
    j = 5;
    while ((j + 12) < MICRO_ECM_PARAM_D)
    {
        // [j+6]Q + [6]Q([j]Q) = [j+12]Q
        uecm_uadd(rho, n, pt6, Pb[k], Pb[i], &Pb[map[j + 12]]);
        i = k;
        k = map[j + 12];
        j += 6;
    }

    // [6]Q + [1]Q([5]Q) = [7]Q
    uecm_uadd(rho, n, pt6, Pb[1], pt5, &Pb[3]);    // <-- [7]Q
    i = 1;
    k = 3;
    j = 1;
    while ((j + 12) < MICRO_ECM_PARAM_D)
    {
        // [j+6]Q + [6]Q([j]Q) = [j+12]Q
        uecm_uadd(rho, n, pt6, Pb[k], Pb[i], &Pb[map[j + 12]]);
        i = k;
        k = map[j + 12];
        j += 6;
    }

    // Pd = [2w]Q
    // [31]Q + [29]Q([2]Q) = [60]Q
    uecm_uadd(rho, n, Pb[9], Pb[10], Pb[2], Pd);   // <-- [60]Q

#ifdef MICRO_ECM_VERBOSE_PRINTF
    ptadds++;
#endif

    // temporary - make [4]Q
    uecm_pt pt4;
    diff1 = submod(Pb[2].X, Pb[2].Z, n);
    sum1 = addmod(Pb[2].X, Pb[2].Z, n);
    uecm_udup(s, rho, n, sum1, diff1, &pt4);   // pt4 = [4]Q


    // make all of the Pbprod's
    for (i = 3; i < 18; i++)
    {
        Pbprod[i] = uecm_mulredc(Pb[i].X, Pb[i].Z, n, rho);
    }


    //initialize info needed for giant step
    uecm_pt Pad;

    // Pd = [w]Q
    // [17]Q + [13]Q([4]Q) = [30]Q
    uecm_uadd(rho, n, Pb[map[17]], Pb[map[13]], pt4, &Pad);    // <-- [30]Q

    // [60]Q + [30]Q([30]Q) = [90]Q
    uecm_uadd(rho, n, *Pd, Pad, Pad, Pa);

    uecm_pt pt90 = *Pa;   // set pt90 = [90]Q
    uecm_pt pt60 = *Pd;   // set pt60 = [60]Q

    // [90]Q + [30]Q([60]Q) = [120]Q
    uecm_uadd(rho, n, *Pa, Pad, *Pd, Pd);

    // [120]Q + [30]Q([90]Q) = [150]Q
    uecm_uadd(rho, n, *Pd, Pad, *Pa, Pa);


    //initialize accumulator
    uint64_t acc = unityval;

    // adjustment of Pa and Pad for particular B1.
    // Currently we have Pa=150, Pd=120, Pad=30

    if (stg1_max < 70)
    {
        // first process the appropriate b's with A=90
        static const int steps27[16] = {59,53,49,47,43,37,31,29,23,19,17,11,7,1,13,41 };
        static const int steps47[15] = {43,37,31,29,23,19,17,11,7,1,13,41,47,49,59 };
        const int* steps;
        int numsteps;
        if (stg1_max == 27)
        {
            steps = steps27;
            numsteps = 16;
        }
        else // if (stg1_max == 47)
        {
            steps = steps47;
            numsteps = 15;
        }

        uint64_t pt90prod = uecm_mulredc(pt90.X, pt90.Z, n, rho);

        for (i = 0; i < numsteps; i++)
        {
            b = steps[i];
            // accumulate the cross product  (zimmerman syntax).
            // page 342 in C&P
            uint64_t tt1 = submod(pt90.X, Pb[map[b]].X, n);
            uint64_t tt2 = addmod(pt90.Z, Pb[map[b]].Z, n);
            uint64_t tt3 = uecm_mulredc(tt1, tt2, n, rho);
            tt1 = addmod(tt3, Pbprod[map[b]], n);
            tt2 = submod(tt1, pt90prod, n);

            uint64_t tmp = uecm_mulredc(acc, tt2, n, rho);
            if (tmp == 0)
                break;
            acc = tmp;
        }
    }
    else if (stg1_max == 70)
    {
        // first process these b's with A=120
        static const int steps[15] = { 49,47,41,37,31,23,19,17,13,11,7,29,43,53,59 };
        // we currently have Pd=120

        uint64_t pdprod = uecm_mulredc(Pd->X, Pd->Z, n, rho);

        for (i = 0; i < 15; i++)
        {
            b = steps[i];
            // accumulate the cross product  (zimmerman syntax).
            // page 342 in C&P
            uint64_t tt1 = submod(Pd->X, Pb[map[b]].X, n);
            uint64_t tt2 = addmod(Pd->Z, Pb[map[b]].Z, n);
            uint64_t tt3 = uecm_mulredc(tt1, tt2, n, rho);
            tt1 = addmod(tt3, Pbprod[map[b]], n);
            tt2 = submod(tt1, pdprod, n);

            uint64_t tmp = uecm_mulredc(acc, tt2, n, rho);
            if (tmp == 0)
                break;
            acc = tmp;
        }
    }
    else if (stg1_max == 165)
    {
        // Currently we have Pa=150, Pd=120, Pad=30,  and pt60=60, pt90=90
        // Need Pa = 180, Pd = 120, Pad = 60
// either of these should be fine
#if 0
        // [150]Q + [30]Q([120]Q) = [180]Q
        uecm_uadd(rho, n, *Pa, Pad, *Pd, Pa);
#else
        diff1 = submod(pt90.X, pt90.Z, n);
        sum1 = addmod(pt90.X, pt90.Z, n);
        uecm_udup(s, rho, n, sum1, diff1, Pa);
#endif
        Pad = pt60;
        // have pa = 180, pd = 120, pad = 60
    }
    else if (stg1_max == 205)
    {
        // Currently we have Pa=150, Pd=120, Pad=30,  and pt60=60, pt90=90
        // need Pa = 210, Pd = 120, Pad = 90

        // [120]Q + [90]Q([30]Q) = [210]Q
        uecm_uadd(rho, n, *Pd, pt90, Pad, Pa);

        Pad = pt90;
    }

    //initialize Paprod
    uint64_t Paprod = uecm_mulredc(Pa->X, Pa->Z, n, rho);

    if (stg1_max == 27)
    {
        barray = b1_27;
        numb = numb1_27;
    }
    else if (stg1_max == 47)
    {
        barray = b1_47;
        numb = numb1_47;
    }
    else if (stg1_max <= 70)
    {
        barray = b1_70;
        numb = numb1_70;
    }
    else if (stg1_max == 85)
    {
        barray = b1_85;
        numb = numb1_85;
    }
    else if (stg1_max == 125)
    {
        barray = b1_125;
        numb = numb1_125;
    }
    else if (stg1_max == 165)
    {
        barray = b1_165;
        numb = numb1_165;
    }
    else if (stg1_max == 205)
    {
        barray = b1_205;
        numb = numb1_205;
    }

    for (i = 0; i < numb; i++)
    {
        if (barray[i] == 0)
        {
            //giant step - use the addition formula for ECM
            uecm_pt point = *Pa;

            //Pa + Pd
            uecm_uadd(rho, n, *Pa, *Pd, Pad, Pa);

            //Pad holds the previous Pa
            Pad = point;

            //and Paprod
            Paprod = uecm_mulredc(Pa->X, Pa->Z, n, rho);

            i++;
        }

        //we accumulate XrZd - XdZr = (Xr - Xd) * (Zr + Zd) + XdZd - XrZr
        //in CP notation, Pa -> (Xr,Zr), Pb -> (Xd,Zd)

        b = barray[i];
        // accumulate the cross product  (zimmerman syntax).
        // page 342 in C&P
        uint64_t tt1 = submod(Pa->X, Pb[map[b]].X, n);
        uint64_t tt2 = addmod(Pa->Z, Pb[map[b]].Z, n);
        uint64_t tt3 = uecm_mulredc(tt1, tt2, n, rho);
        tt1 = addmod(tt3, Pbprod[map[b]], n);
        tt2 = submod(tt1, Paprod, n);

        uint64_t tmp = uecm_mulredc(acc, tt2, n, rho);
        if (tmp == 0)
            break;
        acc = tmp;
    }

    return acc;
}

#ifdef USE_AVX512F

static __m512i uecm_stage2_x8(uecm_pt_x8* P, __m512i rho, __m512i n,
    uint32_t stg1_max, __m512i s, __m512i unityval)
{
    int b;
    int i, j, k;
    uecm_pt_x8 Pa1;
    uecm_pt_x8* Pa = &Pa1;
    uecm_pt_x8 Pb[18];
    uecm_pt_x8 Pd1;
    uecm_pt_x8* Pd = &Pd1;
    const uint8_t* barray = 0;
    int numb;

#ifdef MICRO_ECM_VERBOSE_PRINTF
    ptadds = 0;
    stg1Doub = 0;
    stg1Add = 0;
#endif

    // this function has been written for MICRO_ECM_PARAM_D of 60, so you
    // probably don't want to change it.
    const int MICRO_ECM_PARAM_D = 60;

    //stage 2 init
    //Q = P = result of stage 1
    //compute [d]Q for 0 < d <= MICRO_ECM_PARAM_D

    __m512i Pbprod[18];

    // [1]Q
    Pb[1] = *P;
    Pbprod[1] = uecm_mulredc_x8(Pb[1].X, Pb[1].Z, n, rho);

    // [2]Q
    __m512i diff1 = uecm_submod_x8(P->X, P->Z, n);
    __m512i sum1 = uecm_addmod_x8(P->X, P->Z, n);
    uecm_udup_x8(s, rho, n, sum1, diff1, &Pb[2]);
    //    Pbprod[2] = uecm_mulredc(Pb[2].X, Pb[2].Z, n, rho);    // never used

        /*
        Let D = MICRO_ECM_PARAM_D.

        D is small in tinyecm, so it is straightforward to just enumerate the needed
        points.  We can do it efficiently with two progressions mod 6.
        Pb[0] = scratch
        Pb[1] = [1]Q;
        Pb[2] = [2]Q;
        Pb[3] = [7]Q;   prog2
        Pb[4] = [11]Q;  prog1
        Pb[5] = [13]Q;  prog2
        Pb[6] = [17]Q;  prog1
        Pb[7] = [19]Q;  prog2
        Pb[8] = [23]Q;  prog1
        Pb[9] = [29]Q;  prog1
        Pb[10] = [30 == D]Q;   // shouldn't this be [31]Q?
        Pb[11] = [31]Q; prog2   // shouldn't this be [37]Q?
        Pb[12] = [37]Q; prog2   // shouldn't this be [41]Q?
        Pb[13] = [41]Q; prog1   // [43]Q?
        Pb[14] = [43]Q; prog2   // [47]Q?
        Pb[15] = [47]Q; prog1   // [49]Q?
        Pb[16] = [49]Q; prog2   // [53]Q?
        Pb[17] = [53]Q; prog1   // [59]Q?
        // Pb[18] = [59]Q; prog1   // [60]Q?   note: we can remove this line I believe.  Pb[18] is never set, and never used.  Therefore I changed the definition of Pb above to have only 18 elements.

        two progressions with total of 17 adds to get 15 values of Pb.
        6 + 5(1) -> 11 + 6(5) -> 17 + 6(11) -> 23 + 6(17) -> 29 + 6(23) -> 35 + 6(29) -> 41 + 6(35) -> 47 + 6(41) -> 53 + 6(47) -> 59
        6 + 1(5) -> 7  + 6(1) -> 13 + 6(7)  -> 19 + 6(13) -> 25 + 6(19) -> 31 + 6(25) -> 37 + 6(31) -> 43 + 6(37) -> 49

        we also need [2D]Q = [60]Q
        to get [60]Q we just need one more add:
        compute [60]Q from [31]Q + [29]Q([2]Q), all of which we
        have after the above progressions are computed.

        we also need [A]Q = [((B1 + D) / (2D) * 2 D]Q
        which is equal to the following for various common B1:
        B1      [x]Q
        65      [120]Q
        85      [120]Q
        125     [180]Q      note: according to the A[Q] formula above, wouldn't this be [120]Q?  ( I suspect maybe the formula is supposed to be [((B1 + D) / D) * D]Q )
        165     [180]Q      note: same as above.
        205     [240]Q

        and we need [x-D]Q as well, for the above [x]Q.
        So far we are getting [x]Q and [x-2D]Q each from prac(x,Q).
        There is a better way using progressions of [2D]Q
        [120]Q = 2*[60]Q
        [180]Q = [120]Q + [60]Q([60]Q)
        [240]Q = 2*[120]Q
        [300]Q = [240]Q + [60]Q([180]Q)
        ...
        etc.

        */

    uecm_pt_x8 pt5, pt6;

    // Calculate all Pb: the following is specialized for MICRO_ECM_PARAM_D=60
    // [2]Q + [1]Q([1]Q) = [3]Q
    uecm_uadd_x8(rho, n, Pb[1], Pb[2], Pb[1], &Pb[3]);        // <-- temporary

    // 2*[3]Q = [6]Q
    diff1 = uecm_submod_x8(Pb[3].X, Pb[3].Z, n);
    sum1 = uecm_addmod_x8(Pb[3].X, Pb[3].Z, n);
    uecm_udup_x8(s, rho, n, sum1, diff1, &pt6);   // pt6 = [6]Q

    // [3]Q + [2]Q([1]Q) = [5]Q
    uecm_uadd_x8(rho, n, Pb[3], Pb[2], Pb[1], &pt5);    // <-- pt5 = [5]Q
    Pb[3] = pt5;

    // [6]Q + [5]Q([1]Q) = [11]Q
    uecm_uadd_x8(rho, n, pt6, pt5, Pb[1], &Pb[4]);    // <-- [11]Q

    i = 3;
    k = 4;
    j = 5;
    while ((j + 12) < MICRO_ECM_PARAM_D)
    {
        // [j+6]Q + [6]Q([j]Q) = [j+12]Q
        uecm_uadd_x8(rho, n, pt6, Pb[k], Pb[i], &Pb[map[j + 12]]);
        i = k;
        k = map[j + 12];
        j += 6;
    }

    // [6]Q + [1]Q([5]Q) = [7]Q
    uecm_uadd_x8(rho, n, pt6, Pb[1], pt5, &Pb[3]);    // <-- [7]Q
    i = 1;
    k = 3;
    j = 1;
    while ((j + 12) < MICRO_ECM_PARAM_D)
    {
        // [j+6]Q + [6]Q([j]Q) = [j+12]Q
        uecm_uadd_x8(rho, n, pt6, Pb[k], Pb[i], &Pb[map[j + 12]]);
        i = k;
        k = map[j + 12];
        j += 6;
    }

    // Pd = [2w]Q
    // [31]Q + [29]Q([2]Q) = [60]Q
    uecm_uadd_x8(rho, n, Pb[9], Pb[10], Pb[2], Pd);   // <-- [60]Q

#ifdef MICRO_ECM_VERBOSE_PRINTF
    ptadds++;
#endif

    // temporary - make [4]Q
    uecm_pt_x8 pt4;
    diff1 = uecm_submod_x8(Pb[2].X, Pb[2].Z, n);
    sum1 = uecm_addmod_x8(Pb[2].X, Pb[2].Z, n);
    uecm_udup_x8(s, rho, n, sum1, diff1, &pt4);   // pt4 = [4]Q


    // make all of the Pbprod's
    for (i = 3; i < 18; i++)
    {
        Pbprod[i] = uecm_mulredc_x8(Pb[i].X, Pb[i].Z, n, rho);
    }


    //initialize info needed for giant step
    uecm_pt_x8 Pad;

    // Pd = [w]Q
    // [17]Q + [13]Q([4]Q) = [30]Q
    uecm_uadd_x8(rho, n, Pb[map[17]], Pb[map[13]], pt4, &Pad);    // <-- [30]Q

    // [60]Q + [30]Q([30]Q) = [90]Q
    uecm_uadd_x8(rho, n, *Pd, Pad, Pad, Pa);

    uecm_pt_x8 pt90 = *Pa;   // set pt90 = [90]Q
    uecm_pt_x8 pt60 = *Pd;   // set pt60 = [60]Q

    // [90]Q + [30]Q([60]Q) = [120]Q
    uecm_uadd_x8(rho, n, *Pa, Pad, *Pd, Pd);

    // [120]Q + [30]Q([90]Q) = [150]Q
    uecm_uadd_x8(rho, n, *Pd, Pad, *Pa, Pa);


    //initialize accumulator
    __m512i acc = unityval;

    // adjustment of Pa and Pad for particular B1.
    // Currently we have Pa=150, Pd=120, Pad=30

    if (stg1_max < 70)
    {
        // first process the appropriate b's with A=90
        static const int steps27[16] = { 59,53,49,47,43,37,31,29,23,19,17,11,7,1,13,41 };
        static const int steps47[15] = { 43,37,31,29,23,19,17,11,7,1,13,41,47,49,59 };
        const int* steps;
        int numsteps;
        if (stg1_max == 27)
        {
            steps = steps27;
            numsteps = 16;
        }
        else // if (stg1_max == 47)
        {
            steps = steps47;
            numsteps = 15;
        }

        __m512i pt90prod = uecm_mulredc_x8(pt90.X, pt90.Z, n, rho);

        for (i = 0; i < numsteps; i++)
        {
            b = steps[i];
            // accumulate the cross product  (zimmerman syntax).
            // page 342 in C&P
            __m512i tt1 = uecm_submod_x8(pt90.X, Pb[map[b]].X, n);
            __m512i tt2 = uecm_addmod_x8(pt90.Z, Pb[map[b]].Z, n);
            __m512i tt3 = uecm_mulredc_x8(tt1, tt2, n, rho);
            tt1 = uecm_addmod_x8(tt3, Pbprod[map[b]], n);
            tt2 = uecm_submod_x8(tt1, pt90prod, n);

            __m512i tmp = uecm_mulredc_x8(acc, tt2, n, rho);
            __mmask8 m8 = _mm512_cmpgt_epi64_mask(tmp, _mm512_setzero_si512());
            acc = _mm512_mask_mov_epi64(acc, m8, tmp);
        }
    }
    else if (stg1_max == 70)
    {
        // first process these b's with A=120
        static const int steps[15] = { 49,47,41,37,31,23,19,17,13,11,7,29,43,53,59 };
        // we currently have Pd=120

        __m512i pdprod = uecm_mulredc_x8(Pd->X, Pd->Z, n, rho);

        for (i = 0; i < 15; i++)
        {
            b = steps[i];
            // accumulate the cross product  (zimmerman syntax).
            // page 342 in C&P
            __m512i tt1 = uecm_submod_x8(Pd->X, Pb[map[b]].X, n);
            __m512i tt2 = uecm_addmod_x8(Pd->Z, Pb[map[b]].Z, n);
            __m512i tt3 = uecm_mulredc_x8(tt1, tt2, n, rho);
            tt1 = uecm_addmod_x8(tt3, Pbprod[map[b]], n);
            tt2 = uecm_submod_x8(tt1, pdprod, n);

            __m512i tmp = uecm_mulredc_x8(acc, tt2, n, rho);
            __mmask8 m8 = _mm512_cmpgt_epi64_mask(tmp, _mm512_setzero_si512());
            acc = _mm512_mask_mov_epi64(acc, m8, tmp);
        }
    }
    else if (stg1_max == 165)
    {
        // Currently we have Pa=150, Pd=120, Pad=30,  and pt60=60, pt90=90
        // Need Pa = 180, Pd = 120, Pad = 60
// either of these should be fine
#if 0
        // [150]Q + [30]Q([120]Q) = [180]Q
        uecm_uadd_x8(rho, n, *Pa, Pad, *Pd, Pa);
#else
        diff1 = uecm_submod_x8(pt90.X, pt90.Z, n);
        sum1 = uecm_addmod_x8(pt90.X, pt90.Z, n);
        uecm_udup_x8(s, rho, n, sum1, diff1, Pa);
#endif
        Pad = pt60;
        // have pa = 180, pd = 120, pad = 60
    }
    else if (stg1_max == 205)
    {
        // Currently we have Pa=150, Pd=120, Pad=30,  and pt60=60, pt90=90
        // need Pa = 210, Pd = 120, Pad = 90

        // [120]Q + [90]Q([30]Q) = [210]Q
        uecm_uadd_x8(rho, n, *Pd, pt90, Pad, Pa);

        Pad = pt90;
    }

    //initialize Paprod
    __m512i Paprod = uecm_mulredc_x8(Pa->X, Pa->Z, n, rho);

    if (stg1_max == 27)
    {
        barray = b1_27;
        numb = numb1_27;
    }
    else if (stg1_max == 47)
    {
        barray = b1_47;
        numb = numb1_47;
    }
    else if (stg1_max <= 70)
    {
        barray = b1_70;
        numb = numb1_70;
    }
    else if (stg1_max == 85)
    {
        barray = b1_85;
        numb = numb1_85;
    }
    else if (stg1_max == 125)
    {
        barray = b1_125;
        numb = numb1_125;
    }
    else if (stg1_max == 165)
    {
        barray = b1_165;
        numb = numb1_165;
    }
    else if (stg1_max == 205)
    {
        barray = b1_205;
        numb = numb1_205;
    }

    for (i = 0; i < numb; i++)
    {
        if (barray[i] == 0)
        {
            //giant step - use the addition formula for ECM
            uecm_pt_x8 point = *Pa;

            //Pa + Pd
            uecm_uadd_x8(rho, n, *Pa, *Pd, Pad, Pa);

            //Pad holds the previous Pa
            Pad = point;

            //and Paprod
            Paprod = uecm_mulredc_x8(Pa->X, Pa->Z, n, rho);

            i++;
        }

        //we accumulate XrZd - XdZr = (Xr - Xd) * (Zr + Zd) + XdZd - XrZr
        //in CP notation, Pa -> (Xr,Zr), Pb -> (Xd,Zd)

        b = barray[i];
        // accumulate the cross product  (zimmerman syntax).
        // page 342 in C&P
        __m512i tt1 = uecm_submod_x8(Pa->X, Pb[map[b]].X, n);
        __m512i tt2 = uecm_addmod_x8(Pa->Z, Pb[map[b]].Z, n);
        __m512i tt3 = uecm_mulredc_x8(tt1, tt2, n, rho);
        tt1 = uecm_addmod_x8(tt3, Pbprod[map[b]], n);
        tt2 = uecm_submod_x8(tt1, Paprod, n);

        __m512i tmp = uecm_mulredc_x8(acc, tt2, n, rho);
        __mmask8 m8 = _mm512_cmpgt_epi64_mask(tmp, _mm512_setzero_si512());
        acc = _mm512_mask_mov_epi64(acc, m8, tmp);
    }

    return acc;
}


#endif

static int microecm(uint64_t n, uint64_t *f, uint32_t B1, uint32_t B2,
                     uint32_t curves, uint64_t *ploc_lcg)
{
    //attempt to factor n with the elliptic curve method
    //following brent and montgomery's papers, and CP's book
    int curve;
    int found = 0;
    int result;
    uecm_pt P;
    uint64_t tmp1;

    uint64_t rho = (uint64_t)0 - uecm_multiplicative_inverse(n);

    uint32_t stg1_max = B1;
//    uint32_t stg2_max = B2;

//    uint64_t unityval = u64div(1, n);
    // Let R = 2^64.  We can see R%n ≡ (R-n)%n  (mod n)
    uint64_t unityval = ((uint64_t)0 - n) % n;   // unityval ≡ R  (mod n)

    uint64_t two = addmod(unityval, unityval, n);
    uint64_t four = addmod(two, two, n);
    uint64_t five = addmod(unityval, four, n);
    uint64_t eight = addmod(four, four, n);
    uint64_t sixteen = addmod(eight, eight, n);
    uint64_t two_8 = uecm_sqrredc(sixteen, n, rho);   // R*2^8         (mod n)
    uint64_t two_16 = uecm_sqrredc(two_8, n, rho);    // R*2^16        (mod n)
    uint64_t two_32 = uecm_sqrredc(two_16, n, rho);   // R*2^32        (mod n)
    uint64_t Rsqr = uecm_sqrredc(two_32, n, rho);     // R*2^64 ≡ R*R  (mod n)

    uint64_t s;

    int use_prebuilt_curves = 1;
    uint32_t smallsigma[8] = { 11, 61, 56, 81, 83, 7, 30, 51 };
    uint64_t Z[8] = { 85184, 14526784, 11239424, 34012224, 36594368, 21952, 1728000, 8489664 };
    uint64_t X[8] = { 1560896ULL, 51312965696ULL, 30693697091ULL, 281784327616ULL,
        326229015104ULL, 85184ULL, 716917375ULL, 17495004736ULL };
//    uint64_t u_v[8] = { 373248ULL, 41854210048ULL, 24566036643ULL, 242037319168ULL,
//        281268868608ULL, 4096ULL, 465484375ULL, 13686220288ULL };
//    uint64_t u3_v[8] = { 392, 11392, 9617, 19992, 20984, 160, 2805, 7992 };
    // the negt1 array = u_v[i] * u3_v[i].
    uint64_t negt1[8] = { 146313216ULL, 476803160866816ULL, 236251574395731ULL, 4838810084806656ULL,
                          5902145938870272ULL, 655360ULL, 1305683671875ULL, 109380272541696ULL };
    uint64_t d[8] = { 1098870784ULL, 200325818077184ULL, 110006210374144ULL, 1460769954361344ULL,
        1732928528232448ULL, 38162432ULL, 1376481360000ULL, 57103695458304ULL };

    uint64_t likely_gcd = 1;

    *f = 1;
    for (curve = 0; (uint32_t)curve < curves; curve++)
    {
#ifdef MICRO_ECM_VERBOSE_PRINTF
        stg1Add = 0;
        stg1Doub = 0;
            printf("commencing curve %d of %u\n", curve, curves);
#endif


        if ((curve < 8) && (use_prebuilt_curves))
        {
            //ubuild(&P, rho, &work, goodsigma[curve]); // sigma);
            //sigma = smallsigma[curve];
            // lookup point
            P.X = X[curve];
            P.Z = Z[curve];
            // some computation left to do for S parameter for this 'n'

            uint64_t num, uvc, u3vc;

            uint64_t dem = d[curve];
            // jeff note:  uecm_modinv_64(dem, n) appears to require dem < n,
            // so I now take the remainder to achieve dem < n.

            // This is a faster way to compute dem = dem % n, even if the CPU
            // has extremely fast division (as present in many new CPUs).
            dem = uecm_mulredc(dem, unityval, n, rho);

            num = uecm_mulredc(negt1[curve], Rsqr, n, rho);     // to Monty rep.
               // The mulredc postcondition guarantees  num < n.
            num = n - num;

            dem = uecm_modinv_64(dem, n, &likely_gcd);

            dem = uecm_mulredc(dem, Rsqr, n, rho);              // to Monty rep.
            s = uecm_mulredc(num, dem, n, rho);

            P.X = uecm_mulredc(P.X, Rsqr, n, rho);              // to Monty rep.
            P.Z = uecm_mulredc(P.Z, Rsqr, n, rho);              // to Monty rep.
        }
        else
        {
            likely_gcd = uecm_build(&P, rho, n, ploc_lcg, &s, five, Rsqr);
        }

        if (likely_gcd > 1)
        {
            // If the gcd gave us a factor, we're done.  If not, since gcd != 1
            // the inverse calculated in uecm_build would have bogus, and so this
            // curve is probably set up for failure (hence we continue).
            if (likely_gcd == n || n % likely_gcd != 0)
                continue;
            *f = likely_gcd;
            break;
        }
        
        
#ifdef MICRO_ECM_VERBOSE_PRINTF
        {
            printf("curve parameters:\n");
            printf("\tn = %" PRIu64 "\n", n);
            printf("\trho = %" PRIu64 "\n", rho);
            printf("\tx = %" PRIx64 "\n", P.X);
            printf("\tz = %" PRIx64 "\n", P.Z);
            printf("\tb = %" PRIx64 "\n", s);
        }
#endif
        uecm_stage1(rho, n, &P, (uint64_t)stg1_max, s);
        result = uecm_check_factor(P.Z, n, &tmp1);

#ifdef MICRO_ECM_VERBOSE_PRINTF
        {
            printf("after stage1: P = %" PRIx64 ", %" PRIx64 "\n", P.X, P.Z);
        }
#endif
        if (result == 1)
        {
#ifdef MICRO_ECM_VERBOSE_PRINTF
                printf("\nfound factor %" PRIx64 " in stage 1\n", tmp1);
#endif
            *f = tmp1;
            break;
        }

        if (B2 > B1)
        {
            uint64_t stg2acc = uecm_stage2(&P, rho, n, stg1_max, s, unityval);

#ifdef MICRO_ECM_VERBOSE_PRINTF
            {
                printf("after stage2: A = %" PRIx64 "\n", stg2acc);
            }
                uint32_t paired = 0;
                uint64_t numprimes = 0;
                printf("performed %d pair-multiplies for %" PRIu64 " primes in stage 2\n",
                    paired, numprimes);
                printf("performed %u point-additions and %u point-doubles in stage 2\n",
                    ptadds + stg1Add, stg1Doub);
#endif
            result = uecm_check_factor(stg2acc, n, &tmp1);

            if (result == 1)
            {
#ifdef MICRO_ECM_VERBOSE_PRINTF
                    printf("\nfound factor %" PRIx64 " in stage 2\n", tmp1);
#endif
                *f = tmp1;
                break;
            }
        }
    }

    return curve;
}

#ifdef USE_AVX512F

static void microecm_x8_list(uint64_t* n64, uint64_t* f, uint32_t B1, uint32_t B2,
    uint32_t curves, uint32_t num_in, uint64_t* ploc_lcg)
{
    //attempt to factor n with the elliptic curve method
    //following brent and montgomery's papers, and CP's book
    int result;
    uint8_t msk = 0;
    //uecm_pt P;
    uint64_t uarray[8], rarray[8], tarray[8], oarray[8];
    uint64_t narray[8], carray[8], fivea[8], Rsqra[8];
    //uint64_t xarray[8], zarray[8], sarray[8];
    //uint64_t likely_gcd;
    uint32_t stg1_max = B1;
    __m512i v_u, v_s, v_n, v_f, v_r, v_c, v_o;
    uecm_pt_x8 v_P;
    int i;
    uint32_t j = 0;

    // default (failed) factors
    for (i = 0; i < num_in; i++)
    {
        f[i] = 1;
    }

    v_u = v_n = v_r = v_s = v_c = v_o = v_P.X = v_P.Z = _mm512_setzero_si512();;

    while (1)
    {
        uint8_t lmsk = 0;

        // this keeps the vectors full
        for (i = 0; i < 8; i++)
        {
            if (((msk & (1 << i)) == 0) && (j < num_in))
            {
                // compute things we need for this new n.
                uint64_t unityval = (1ULL << 52) % n64[j];   // unityval ≡ R  (mod n)
                uint64_t rho = uecm_multiplicative_inverse52(n64[j]);
                uint64_t two = addmod52(unityval, unityval, n64[j]);
                uint64_t four = addmod52(two, two, n64[j]);
                uint64_t five = addmod52(unityval, four, n64[j]);
                uint64_t eight = addmod52(four, four, n64[j]);
                uint64_t sixteen = addmod52(eight, eight, n64[j]);

                uint64_t two_5 = addmod52(sixteen, sixteen, n64[j]);
                uint64_t two_8 = uecm_sqrredc52(sixteen, n64[j], rho);        // R*2^8         (mod n)
                uint64_t two_13 = uecm_mulredc52(two_8, two_5, n64[j], rho);  // R*2^13         (mod n)
                uint64_t two_26 = uecm_sqrredc52(two_13, n64[j], rho);        // R*2^26        (mod n)
                uint64_t Rsqr = uecm_sqrredc52(two_26, n64[j], rho);          // R*2^52 ≡ R*R  (mod n)

                //uint64_t two_8 = uecm_sqrredc52(sixteen, n64[j], rho);        // R*2^8         (mod n)
                //uint64_t two_16 = uecm_sqrredc52(two_8, n64[j], rho);         // R*2^16        (mod n)
                //uint64_t two_32 = uecm_sqrredc52(two_16, n64[j], rho);        // R*2^32        (mod n)
                //uint64_t two_48 = uecm_mulredc52(two_16, two_32, n64[j], rho);// R*2^48        (mod n)
                //uint64_t Rsqr = uecm_mulredc52(two_48, sixteen, n64[j], rho); // R*2^52 ≡ R*R  (mod n)
                uint64_t s;

                // mark this index for loading into vectors
                lmsk |= (1 << i);

                // load the info we need into our vectors.
                uarray[i] = unityval;
                rarray[i] = rho;
                narray[i] = n64[j];
                carray[i] = 0;
                oarray[i] = j;
                fivea[i] = five;
                Rsqra[i] = Rsqr;

                j++;
            }

            //// build a curve for this index
            //uint64_t s;
            //likely_gcd = uecm_build52(&P, rarray[i], narray[i], ploc_lcg,
            //    &s, fivea[i], Rsqra[i]);
            //
            //sarray[i] = s;
            //xarray[i] = P.X;
            //zarray[i] = P.Z;
            //
            //while (likely_gcd > 1)
            //{
            //    // If the gcd gave us a factor, we're done.  If not, since gcd != 1
            //    // the inverse calculated in uecm_build would have bogus, and so this
            //    // curve is probably set up for failure (hence we continue).
            //    if (likely_gcd == narray[i] || narray[i] % likely_gcd != 0)
            //    {
            //        likely_gcd = uecm_build52(&P, rarray[i], narray[i], ploc_lcg,
            //            &s, fivea[i], Rsqra[i]);
            //
            //        sarray[i] = s;
            //        xarray[i] = P.X;
            //        zarray[i] = P.Z;
            //    }
            //    else
            //    {
            //        f[oarray[i]] = likely_gcd;
            //        //msk &= (~(1 << i));
            //        i--;
            //        break;
            //    }
            //}
        }

        msk |= lmsk;
        // these are conditionally loaded depending on if
        // the input is newly loaded.
        v_u = _mm512_mask_loadu_epi64(v_u, lmsk, uarray);
        v_n = _mm512_mask_loadu_epi64(v_n, lmsk, narray);
        v_r = _mm512_mask_loadu_epi64(v_r, lmsk, rarray);
        v_c = _mm512_mask_loadu_epi64(v_c, lmsk, carray);
        v_o = _mm512_mask_loadu_epi64(v_o, lmsk, oarray);

        // build the curves
        v_s = uecm_build_x8(&v_P, v_r, v_n, ploc_lcg, tarray,
            _mm512_loadu_epi64(fivea), _mm512_loadu_epi64(Rsqra));
        
        // process any factors or errors resulting from modinv
        __m512i lgcd = _mm512_loadu_epi64(tarray);
        __mmask8 m1 = _mm512_cmplt_epi64_mask(lgcd, v_n);
        __mmask8 m2 = _mm512_cmpgt_epi64_mask(lgcd, uecm_set64(1));
        __m512i rem = _mm512_rem_epi64(v_n, lgcd);
        __mmask8 m3 = _mm512_cmpeq_epi64_mask(rem, _mm512_setzero_si512());
        _mm512_mask_store_epi64(tarray, m2 & msk & m1 & m3, lgcd);
        msk &= (~(m2 & m1 & m3));
        
        for (i = 0; i < 8; i++)
        {
            if (((m1 & m2 & m3) & (1 << i)) > 0)
                f[oarray[i]] = tarray[i];
        }

        // these are always loaded for this new curve
        //v_s = _mm512_loadu_epi64(sarray);
        //v_P.X = _mm512_loadu_epi64(xarray);
        //v_P.Z = _mm512_loadu_epi64(zarray);

        uecm_stage1_x8(v_r, v_n, &v_P, (uint64_t)stg1_max, v_s);

        _mm512_storeu_epi64(tarray, v_P.Z);
        for (i = 0; i < 8; i++)
        {
            if (msk & (1 << i))
            {
                uint64_t t;
                result = uecm_check_factor(tarray[i], narray[i], &t);

                if (result == 1)
                {
                    f[oarray[i]] = t;
                    msk &= (~(1 << i));
                }
            }
        }

        if (B2 > B1)
        {
            __m512i stg2acc = uecm_stage2_x8(&v_P, v_r, v_n, stg1_max, v_s, v_u);

            _mm512_storeu_epi64(tarray, stg2acc);

            for (i = 0; i < 8; i++)
            {
                if (msk & (1 << i))
                {
                    uint64_t t;
                    result = uecm_check_factor(tarray[i], narray[i], &t);

                    if (result == 1)
                    {
                        f[oarray[i]] = t;
                        msk &= (~(1 << i));
                    }
                }
            }
        }

        // if we've exhaused the curves count for anything in
        // the vector, flag it for replacement.
        v_c = _mm512_add_epi64(v_c, uecm_set64(1));
        msk &= _mm512_cmplt_epi64_mask(v_c, uecm_set64(curves));

        if ((j == num_in) && (msk == 0))
            break;
    }

    return;
}

#endif

static uint64_t uecm_dispatch(uint64_t n, int targetBits, int arbitrary, uint64_t *ploc_lcg)
{
    int B1, curves;
    uint64_t f64 = 1;
#ifdef MICRO_ECM_VERBOSE_PRINTF
    //uecm_stage2_pair(70, 180, 120, 120, 30, primes);
    //uecm_stage2_pair(180, 25 * 70, 150, 120, 30, primes);

    //uecm_stage2_pair(47, 150, 90, 120, 30, primes);
    //uecm_stage2_pair(150, 25 * 47, 150, 120, 30, primes);
    
    //uecm_stage2_pair(30, 150, 90, 120, 30, primes);
    //uecm_stage2_pair(150, 25 * 30, 150, 120, 30, primes);
    //exit(0);
#endif
    if (arbitrary)
    {
        // try fast attempts to find possible small factors.
        {
            B1 = 47;
            curves = 1;
            microecm(n, &f64, B1, 25 * B1, curves, ploc_lcg);
            if (f64 > 1)
                return f64;
        }
        {
            B1 = 70;
            curves = 1;
            microecm(n, &f64, B1, 25 * B1, curves, ploc_lcg);
            if (f64 > 1)
                return f64;
        }
        if (targetBits > 58)
        {
            B1 = 125;
            curves = 1;
            microecm(n, &f64, B1, 25 * B1, curves, ploc_lcg);
            if (f64 > 1)
                return f64;
        }
    }

    if (targetBits <= 40)
    {
        B1 = 27;
        curves = 32;
        microecm(n, &f64, B1, 25 * B1, curves, ploc_lcg);
    }
    else if (targetBits <= 44)
    {
        B1 = 47;
        curves = 32;
        microecm(n, &f64, B1, 25 * B1, curves, ploc_lcg);
    }
    else if (targetBits <= 48)
    {
        // multi-thread issue here...
        //f64 = LehmanFactor(n, 0, 0, 0);
        B1 = 70;
        curves = 32;
        microecm(n, &f64, B1, 25 * B1, curves, ploc_lcg);
    }
    else if (targetBits <= 52)
    {
        B1 = 85;
        curves = 32;
        microecm(n, &f64, B1, 25 * B1, curves, ploc_lcg);
    }
    else if (targetBits <= 58)
    {
        B1 = 125;
        curves = 32;
        microecm(n, &f64, B1, 25 * B1, curves, ploc_lcg);
    }
    else if (targetBits <= 62)
    {
        B1 = 165;
        curves = 42;
        microecm(n, &f64, B1, 25 * B1, curves, ploc_lcg);
    }
    else if (targetBits <= 64)
    {
        B1 = 205;
        curves = 42;
        microecm(n, &f64, B1, 25 * B1, curves, ploc_lcg);
    }

    return f64;
}

#ifdef USE_AVX512F

static void uecm_dispatch_x8_list(uint64_t* n, uint64_t* f,
    int targetBits, uint32_t num_in, uint64_t* ploc_lcg)
{
    int B1, curves;

    if (targetBits <= 40)
    {
        B1 = 27;
        curves = 32;
        microecm_x8_list(n, f, B1, 25 * B1, curves, num_in, ploc_lcg);
    }
    else if (targetBits <= 44)
    {
        B1 = 47;
        curves = 32;
        microecm_x8_list(n, f, B1, 25 * B1, curves, num_in, ploc_lcg);
    }
    else if (targetBits <= 48)
    {
        B1 = 70;
        curves = 32;
        microecm_x8_list(n, f, B1, 25 * B1, curves, num_in, ploc_lcg);
    }
    else if (targetBits <= 52)
    {
        B1 = 85;
        curves = 32;
        microecm_x8_list(n, f, B1, 25 * B1, curves, num_in, ploc_lcg);
    }
    else if (targetBits <= 58)
    {
        B1 = 125;
        curves = 32;
        microecm_x8_list(n, f, B1, 25 * B1, curves, num_in, ploc_lcg);
    }
    else if (targetBits <= 62)
    {
        B1 = 165;
        curves = 42;
        microecm_x8_list(n, f, B1, 25 * B1, curves, num_in, ploc_lcg);
    }
    else if (targetBits <= 64)
    {
        B1 = 205;
        curves = 42;
        microecm_x8_list(n, f, B1, 25 * B1, curves, num_in, ploc_lcg);
    }

    return;
}

#endif

static int uecm_get_bits(uint64_t n)
{
    int i = 0;
    while (n != 0)
    {
        n >>= 1;
        i++;
    }
    return i;
}


// getfactor_uecm() returns 1 if unable to find a factor of q64,
// Otherwise it returns a factor of q64.
// 
// if the input is known to have no small factors, set is_arbitrary=0, 
// otherwise, set is_arbitrary=1 and a few curves targetting small factors
// will be run prior to the standard sequence of curves for the input size.
//  
// Prior to your first call of getfactor_uecm(), set *pran = 0  (or set it to
// some other arbitrary value); after that, don't change *pran.
// FYI: *pran is used within this file by a random number generator, and it
// holds the current value of a pseudo random sequence.  Your first assigment
// to *pran seeds the sequence, and after seeding it you don't want to
// change *pran, since that would restart the sequence.
uint64_t getfactor_uecm(uint64_t q64, int is_arbitrary, uint64_t *pran)
{
    if (q64 % 2 == 0)
        return 2;
    int bits = uecm_get_bits(q64);
    return uecm_dispatch(q64, bits, is_arbitrary, pran);
}

#ifdef USE_AVX512F

void getfactor_uecm_x8_list(uint64_t* q64, uint64_t* f64, uint32_t num_in, uint64_t* pran)
{
    dbias = _mm512_castsi512_pd(uecm_set64(0x4670000000000000ULL));
    vbias1 = uecm_set64(0x4670000000000000ULL);
    vbias2 = uecm_set64(0x4670000000000001ULL);
    lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);

    int bits = uecm_get_bits(q64[0]);       // assume all the same size
    uecm_dispatch_x8_list(q64, f64, bits, num_in, pran);
    return;
}

#else

void getfactor_uecm_x8_list(uint64_t* q64, uint64_t* f64, uint32_t num_in, uint64_t* pran)
{
    int bits = uecm_get_bits(q64[0]);       // assume all the same size
    int i;

    for (i = 0; i < num_in; i++)
    {
        f64[i] = uecm_dispatch(q64[i], bits, 0, pran);
    }
    return;
}

#endif
