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
#include "arith.h"
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

#include <stdio.h>


#ifdef USE_AVX512F
#  include <immintrin.h>
#endif


// Using the inline asm in this file can increase performance by ~20-25%
// (surprisingly).  Hence these macros are defined by default.
#if defined(__x86_64__) || defined(_M_X64)
#  define MICRO_PM1_ALT_MULREDC_USE_INLINE_ASM_X86
#endif

// We rarely will ever want to use debugging printfs
//#define MICRO_PM1_VERBOSE_PRINTF




#ifdef _MSC_VER
#  define MICRO_PM1_FORCE_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER) || defined (__INTEL_LLVM_COMPILER)
#  define MICRO_PM1_FORCE_INLINE inline __attribute__((always_inline))
#else
#  define MICRO_PM1_FORCE_INLINE __inline
#endif


#ifdef MICRO_PM1_VERBOSE_PRINTF
#  include <stdio.h>
// these globals will cause data races, if we run multithreaded
    uint32_t stg1Doub;
    uint32_t stg1Add;
    uint32_t ptadds;
#endif

static const uint32_t map[60] = {
    0, 1, 2, 0, 0, 0, 0, 3, 0, 0,
    0, 4, 0, 5, 0, 0, 0, 6, 0, 7,
    0, 0, 0, 8, 0, 0, 0, 0, 0, 9,
    0, 10, 0, 0, 0, 0, 0, 11, 0, 0,
    0, 12, 0, 13, 0, 0, 0, 14, 0, 15,
    0, 0, 0, 16, 0, 0, 0, 0, 0, 17 };

static const double INV_2_POW_32 = 1.0 / (double)((uint64_t)(1) << 32);

#define NUMP 1438
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
6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, 6229,
6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, 6311,
6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 6373,
6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, 6481,
6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, 6577,
6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, 6679,
6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, 6763,
6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, 6841,
6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 6947,
6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997, 7001,
7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, 7109,
7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, 7211,
7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, 7307,
7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, 7417,
7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, 7507,
7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, 7573,
7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, 7649,
7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, 7727,
7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, 7841,
7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919, 7927,
7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017, 8039,
8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111, 8117,
8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219, 8221,
8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291, 8293,
8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387, 8389,
8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501, 8513,
8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597, 8599,
8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677, 8681,
8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, 8747,
8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831, 8837,
8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929, 8933,
8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011, 9013,
9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109, 9127,
9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199, 9203,
9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283, 9293,
9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377, 9391,
9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439, 9461,
9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533, 9539,
9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631, 9643,
9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733, 9739,
9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811, 9817,
9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887, 9901,
9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973, 10007, 10009,
10037, 10039, 10061, 10067, 10069, 10079, 10091, 10093, 10099, 10103, 10111,
10133, 10139, 10141, 10151, 10159, 10163, 10169, 10177, 10181, 10193, 10211,
10223, 10243, 10247, 10253, 10259, 10267, 10271, 10273, 10289, 10301, 10303,
10313, 10321, 10331, 10333, 10337, 10343, 10357, 10369, 10391, 10399, 10427,
10429, 10433, 10453, 10457, 10459, 10463, 10477, 10487, 10499, 10501, 10513,
10529, 10531, 10559, 10567, 10589, 10597, 10601, 10607, 10613, 10627, 10631,
10639, 10651, 10657, 10663, 10667, 10687, 10691, 10709, 10711, 10723, 10729,
10733, 10739, 10753, 10771, 10781, 10789, 10799, 10831, 10837, 10847, 10853,
10859, 10861, 10867, 10883, 10889, 10891, 10903, 10909, 10937, 10939, 10949,
10957, 10973, 10979, 10987, 10993, 11003, 11027, 11047, 11057, 11059, 11069,
11071, 11083, 11087, 11093, 11113, 11117, 11119, 11131, 11149, 11159, 11161,
11171, 11173, 11177, 11197, 11213, 11239, 11243, 11251, 11257, 11261, 11273,
11279, 11287, 11299, 11311, 11317, 11321, 11329, 11351, 11353, 11369, 11383,
11393, 11399, 11411, 11423, 11437, 11443, 11447, 11467, 11471, 11483, 11489,
11491, 11497, 11503, 11519, 11527, 11549, 11551, 11579, 11587, 11593, 11597,
11617, 11621, 11633, 11657, 11677, 11681, 11689, 11699, 11701, 11717, 11719,
11731, 11743, 11777, 11779, 11783, 11789, 11801, 11807, 11813, 11821, 11827,
11831, 11833, 11839, 11863, 11867, 11887, 11897, 11903, 11909, 11923, 11927,
11933, 11939, 11941, 11953, 11959, 11969, 11971, 11981, 11987, };


static uint32_t upm1_lcg_rand_32B(uint32_t lower, uint32_t upper, uint64_t *ploc_lcg)
{
    *ploc_lcg = 6364136223846793005ULL * (*ploc_lcg) + 1442695040888963407ULL;
    return lower + (uint32_t)(
        (double)(upper - lower) * (double)((*ploc_lcg) >> 32) * INV_2_POW_32);
}

__inline uint64_t upm1_submod(uint64_t a, uint64_t b, uint64_t n);
__inline uint64_t upm1_addmod(uint64_t a, uint64_t b, uint64_t n);

#if defined(MICRO_PM1_ALT_MULREDC_USE_INLINE_ASM_X86) && !defined(_MSC_VER)

MICRO_PM1_FORCE_INLINE uint64_t upm1_submod(uint64_t a, uint64_t b, uint64_t n)
{
    __asm__(
        "xorq %%r8, %%r8 \n\t"
        "subq %1, %0 \n\t"
        "cmovc %2, %%r8 \n\t"
        "addq %%r8, %0 \n\t"
        : "+r"(a)
        : "r"(b), "r"(n)
        : "r8", "cc");

    return a;
}

MICRO_PM1_FORCE_INLINE uint64_t upm1_addmod(uint64_t x, uint64_t y, uint64_t n)
{
    uint64_t t = x - n;
    x += y;
    __asm__("add %2, %1\n\t"
        "cmovc %1, %0\n\t"
        :"+r" (x), "+&r" (t)
        : "r" (y)
        : "cc"
    );
    return x;
}

#else

MICRO_PM1_FORCE_INLINE uint64_t upm1_submod(uint64_t a, uint64_t b, uint64_t n)
{
    uint64_t r0;
    if (_subborrow_u64(0, a, b, &r0))
        r0 += n;
    return r0;
}

MICRO_PM1_FORCE_INLINE uint64_t upm1_addmod(uint64_t x, uint64_t y, uint64_t n)
{
#if 0
    uint64_t r;
    uint64_t tmp = x - n;
    uint8_t c = _addcarry_u64(0, tmp, y, &r);
    return (c) ? r : x + y;
#else
    // FYI: The clause above often compiles with a branch in MSVC.
    // The statement below often compiles without a branch (uses cmov) in MSVC.
    return (x >= n - y) ? x - (n - y) : x + y;
#endif
}

#endif

/* --- The following two functions are written by Jeff Hurchalla, Copyright 2022 --- */

// for this algorithm, see https://jeffhurchalla.com/2022/04/28/montgomery-redc-using-the-positive-inverse-mod-r/
MICRO_PM1_FORCE_INLINE static uint64_t upm1_mulredc_alt(uint64_t x, uint64_t y, uint64_t N, uint64_t invN)
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
#if defined(MICRO_PM1_ALT_MULREDC_USE_INLINE_ASM_X86) && !defined(_MSC_VER)
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
static uint64_t upm1_multiplicative_inverse(uint64_t a)
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

static uint64_t upm1_bingcd64(uint64_t u, uint64_t v)
{
#if 1
    if (u == 0) {
        return v;
    }
    if (v != 0) {
        int j = _trail_zcnt64(v);
        v = (uint64_t)(v >> j);
        while (1) {
            uint64_t tmp = u;
            uint64_t sub1 = (uint64_t)(v - tmp);
            uint64_t sub2 = (uint64_t)(tmp - v);
            if (tmp == v)
                break;
            u = (tmp >= v) ? v : tmp;
            v = (tmp >= v) ? sub2 : sub1;
            // For the line below, the standard way to write this algorithm
            // would have been to use _trail_zcnt64(v)  (instead of
            // _trail_zcnt64(sub1)).  However, as pointed out by
            // https://gmplib.org/manual/Binary-GCD, "in twos complement the
            // number of low zero bits on u-v is the same as v-u, so counting or
            // testing can begin on u-v without waiting for abs(u-v) to be
            // determined."  Hence we are able to use sub1 for the argument.
            // By removing the dependency on abs(u-v), the CPU can execute
            // _trail_zcnt64() at the same time as abs(u-v).
            j = _trail_zcnt64(sub1);
            v = (uint64_t)(v >> j);
        }
    }
    return u;
#else
    // For reference, or if in the future we need to allow an even u,
    // this version allows u to be even or odd.
    if (u == 0) {
        return v;
    }
    if (v != 0) {
        int i = _trail_zcnt64(u);
        int j = _trail_zcnt64(v);
        u = (uint64_t)(u >> i);
        v = (uint64_t)(v >> j);
        int k = (i < j) ? i : j;
        while (1) {
            uint64_t tmp = u;
            uint64_t sub1 = (uint64_t)(v - tmp);
            uint64_t sub2 = (uint64_t)(tmp - v);
            if (tmp == v)
                break;
            u = (tmp >= v) ? v : tmp;
            v = (tmp >= v) ? sub2 : sub1;
            // For the line below, the standard way to write this algorithm
            // would have been to use _trail_zcnt64(v)  (instead of
            // _trail_zcnt64(sub1)).  However, as pointed out by
            // https://gmplib.org/manual/Binary-GCD, "in twos complement the
            // number of low zero bits on u-v is the same as v-u, so counting or
            // testing can begin on u-v without waiting for abs(u-v) to be
            // determined."  Hence we are able to use sub1 for the argument.
            // By removing the dependency on abs(u-v), the CPU can execute
            // _trail_zcnt64() at the same time as abs(u-v).
            j = _trail_zcnt64(sub1);
            v = (uint64_t)(v >> j);
        }
        u = (uint64_t)(u << k);
    }
    return u;
#endif
}

/* --- end Hurchalla functions --- */


// full strength mul/sqr redc
MICRO_PM1_FORCE_INLINE static uint64_t upm1_mulredc(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
{
    return upm1_mulredc_alt(x, y, n, 0 - nhat);
}
MICRO_PM1_FORCE_INLINE static uint64_t upm1_sqrredc(uint64_t x, uint64_t n, uint64_t nhat)
{
    return upm1_mulredc_alt(x, x, n, 0 - nhat);
}

static int upm1_check_factor(uint64_t Z, uint64_t n, uint64_t* f)
{
    int status = 0;
    *f = upm1_bingcd64(n, Z);

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

static int ewin100[34] = { // 170 muls
        12, 12, 14, 3, 12, 7, 13, 6, 12, 0, 15, 4, 1, 8, 7, 3, 0, 14, 13, 6, 
        5, 6, 15, 13, 0, 11, 0, 14, 13, 3, 8, 8, 12, 0 };
static int ewin333[119] = { // 595 muls
    2, 5, 1, 6, 12, 5, 1, 5, 0, 15, 3, 13, 2, 4, 2, 0, 4, 13, 11, 9, 4, 5, 
    4, 13, 7, 15, 0, 11, 10, 7, 5, 4, 7, 0, 14, 11, 12, 10, 12, 4, 11, 2, 
    5, 2, 10, 10, 7, 3, 14, 11, 8, 0, 15, 2, 2, 3, 10, 11, 6, 9, 2, 8, 15, 
    4, 12, 13, 14, 13, 0, 7, 3, 3, 12, 3, 9, 8, 4, 6, 15, 0, 3, 9, 11, 14, 
    5, 7, 4, 4, 11, 14, 8, 6, 11, 14, 0, 12, 13, 4, 12, 12, 11, 6, 13, 0, 
    10, 8, 13, 6, 13, 15, 2, 5, 14, 13, 11, 11, 15, 0, 0 };
static int ewin666[243] = { // 1215 muls
    6, 15, 12, 1, 15, 8, 2, 2, 13, 7, 9, 2, 13, 0, 6, 10, 14, 1, 4, 14, 0, 
    10, 9, 6, 0, 9, 11, 0, 11, 2, 5, 5, 12, 14, 9, 11, 7, 11, 7, 9, 2, 1, 1, 
    2, 1, 9, 14, 3, 3, 1, 0, 2, 12, 15, 0, 10, 11, 14, 10, 9, 11, 0, 13, 6, 
    5, 10, 7, 14, 7, 4, 5, 3, 11, 15, 0, 9, 1, 2, 12, 13, 0, 9, 6, 10, 14, 0, 
    11, 8, 7, 4, 15, 10, 10, 6, 8, 11, 9, 8, 1, 2, 7, 4, 0, 9, 6, 7, 11, 6, 
    14, 0, 1, 12, 8, 11, 13, 5, 4, 0, 13, 14, 2, 0, 8, 12, 2, 8, 3, 10, 14, 11, 
    0, 15, 7, 9, 1, 7, 10, 3, 8, 11, 9, 8, 10, 9, 6, 0, 1, 13, 2, 8, 14, 0, 8, 
    13, 0, 0, 6, 13, 5, 15, 4, 6, 3, 0, 8, 15, 3, 14, 12, 11, 5, 6, 14, 1, 4, 
    6, 9, 14, 1, 0, 3, 15, 6, 8, 0, 0, 8, 4, 10, 0, 1, 9, 9, 15, 14, 14, 8, 
    12, 14, 14, 13, 8, 6, 7, 6, 11, 9, 6, 7, 7, 6, 8, 15, 1, 15, 2, 13, 9, 5, 
    7, 10, 10, 11, 11, 2, 15, 4, 3, 10, 6, 9, 2, 14, 8, 14, 5, 6, 2, 15, 12, 
    10, 0, 0 };
static int ewin1000[360] = { // 1800 muls
    3, 11, 15, 13, 0, 1, 12, 9, 9, 13, 3, 2, 15, 9, 11, 15, 2, 13, 7, 4, 7, 8, 
    10, 1, 3, 7, 10, 11, 3, 9, 12, 4, 14, 5, 7, 0, 14, 9, 6, 0, 0, 12, 1, 10, 
    3, 8, 6, 1, 4, 3, 12, 15, 5, 11, 14, 14, 5, 1, 4, 2, 9, 2, 0, 7, 7, 14, 2, 
    12, 10, 3, 4, 0, 10, 9, 2, 2, 3, 2, 4, 2, 3, 12, 5, 11, 7, 11, 0, 14, 2, 
    4, 2, 7, 11, 3, 5, 1, 6, 10, 0, 2, 1, 11, 0, 10, 11, 4, 10, 0, 1, 7, 14, 1, 
    14, 7, 2, 5, 0, 6, 12, 0, 7, 15, 1, 12, 15, 15, 6, 1, 5, 13, 0, 12, 4, 8, 
    12, 8, 13, 14, 5, 12, 2, 5, 15, 3, 8, 6, 12, 11, 14, 1, 14, 2, 8, 1, 3, 13, 
    2, 1, 14, 10, 15, 14, 9, 15, 5, 1, 2, 14, 11, 11, 9, 8, 7, 10, 0, 11, 6, 
    11, 3, 3, 11, 12, 15, 11, 7, 12, 0, 14, 11, 4, 7, 3, 10, 15, 13, 4, 6, 12, 
    14, 11, 14, 13, 8, 8, 5, 9, 2, 15, 4, 1, 11, 7, 6, 11, 8, 6, 2, 1, 12, 10, 
    12, 11, 3, 15, 1, 9, 1, 14, 5, 15, 4, 12, 3, 8, 10, 9, 7, 4, 15, 11, 9, 1, 
    1, 13, 7, 9, 1, 13, 6, 9, 0, 11, 14, 15, 12, 9, 3, 15, 8, 14, 15, 10, 7, 
    11, 15, 11, 15, 15, 14, 5, 15, 0, 1, 10, 2, 13, 15, 5, 3, 10, 14, 13, 11, 
    7, 12, 4, 14, 1, 5, 14, 8, 10, 14, 3, 13, 12, 13, 4, 2, 9, 13, 2, 2, 14, 
    11, 3, 7, 5, 9, 1, 7, 10, 5, 15, 14, 2, 4, 2, 4, 4, 4, 0, 5, 11, 0, 10, 5, 
    1, 4, 13, 0, 13, 10, 5, 9, 6, 2, 7, 0, 10, 10, 9, 8, 14, 12, 8, 6, 0, 8, 5, 
    14, 8, 9, 4, 12, 1, 7, 10, 0, 0 };

static uint64_t upm1_stage1(uint64_t rho, uint64_t n, uint64_t one, uint64_t stg1)
{
    int i;
    uint64_t P = one;
    uint64_t g[16];

    g[0] = one;
    for (i = 1; i < 16; i++)
        g[i] = upm1_addmod(g[i-1], g[i-1], n);

    switch (stg1)
    {
    case 100:
        for (i = 0; i < 34; i++) {
            P = upm1_sqrredc(P, n, rho);
            P = upm1_sqrredc(P, n, rho);
            P = upm1_sqrredc(P, n, rho);
            P = upm1_sqrredc(P, n, rho);
            if (ewin100[i] > 0) P = upm1_mulredc(P, g[ewin100[i]], n, rho);
        }

        break;
    case 333:
        for (i = 0; i < 119; i++) {
            P = upm1_sqrredc(P, n, rho);
            P = upm1_sqrredc(P, n, rho);
            P = upm1_sqrredc(P, n, rho);
            P = upm1_sqrredc(P, n, rho);
            if (ewin333[i] > 0) P = upm1_mulredc(P, g[ewin333[i]], n, rho);
        }

        break;
    case 666:
        for (i = 0; i < 243; i++) {
            P = upm1_sqrredc(P, n, rho);
            P = upm1_sqrredc(P, n, rho);
            P = upm1_sqrredc(P, n, rho);
            P = upm1_sqrredc(P, n, rho);
            if (ewin666[i] > 0) P = upm1_mulredc(P, g[ewin666[i]], n, rho);
        }

        break;
    case 1000:
        for (i = 0; i < 360; i++) {
            P = upm1_sqrredc(P, n, rho);
            P = upm1_sqrredc(P, n, rho);
            P = upm1_sqrredc(P, n, rho);
            P = upm1_sqrredc(P, n, rho);
            if (ewin1000[i] > 0) P = upm1_mulredc(P, g[ewin1000[i]], n, rho);
        }

        break;
    }
    return P;
}

void upm1_generate_window_plan(uint64_t B1)
{
    int i, j;
    mpz_t e;
    uint64_t q;
    mpz_init(e);

    uint64_t p[168] = { 2, 3, 5, 7,
        11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
        101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179,
        181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269,
        271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367,
        373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461,
        463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571,
        577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661,
        673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773,
        787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883,
        887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997 };
    int g[1024];

    mpz_set_ui(e, 1);
    for (i = 0; i < 168; i++)
    {
        if (p[i] > B1)
            break;

        q = p[i];
        while ((q * p[i]) < B1)
        {
            q *= p[i];
        }
        mpz_mul_ui(e, e, q);
        gmp_printf("e*=%lu=%Zx\n", q, e);
    }

    // now we have the exponent, generate a LR-kary windowing plan.
    // assume window size 4.
    j = 0;
    while (mpz_cmp_ui(e, 0) > 0)
    {
        g[j++] =  mpz_get_ui(e) & 0xf;
        mpz_tdiv_q_2exp(e, e, 4);
    }
    printf("%d: \n", j);

    while (j > 0)
        printf("%d,", g[--j]);
    printf("\n");

    mpz_clear(e);
    return;
}

#ifdef MICRO_PM1_VERBOSE_PRINTF
static void upm1_stage2_pair(int b1, int b2, 
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

static uint64_t upm1_expRL(uint64_t P, uint64_t n, uint64_t rho, uint64_t one, uint32_t m)
{
    uint64_t s = P;
    P = one;

    while (m > 0)
    {
        if (m & 1)
            P = upm1_mulredc(P, s, n, rho);
        s = upm1_sqrredc(s, n, rho);
        m >>= 1;
    }
    return P;
}

static const uint32_t upm1_map[60] = {
    0, 1, 2, 0, 0, 0, 0, 3, 0, 0,
    0, 4, 0, 5, 0, 0, 0, 6, 0, 7,
    0, 0, 0, 8, 0, 0, 0, 0, 0, 9,
    0, 10, 0, 0, 0, 0, 0, 11, 0, 0,
    0, 12, 0, 13, 0, 0, 0, 14, 0, 15,
    0, 0, 0, 16, 0, 0, 0, 0, 0, 17 };

/* baby-steps, giant-steps for b1 = 100, w = 60
q0 = 120 */
static int num_b1_100_steps = 382;
static uint8_t b1_100_steps[382] = {
19, 17, 13, 11, 7, 0,
53, 49, 43, 41, 31, 29, 23, 17, 13, 7, 1, 0,
59, 49, 47, 43, 41, 29, 17, 13, 11, 7, 1, 0,
59, 49, 43, 37, 31, 29, 23, 19, 17, 7, 0,
53, 49, 47, 43, 29, 23, 13, 11, 7, 1, 0,
53, 47, 41, 37, 31, 23, 19, 11, 1, 0,
59, 49, 47, 41, 37, 31, 23, 19, 17, 13, 1, 0,
53, 49, 41, 37, 31, 19, 17, 0,
59, 53, 43, 37, 31, 29, 23, 13, 7, 1, 0,
59, 53, 47, 43, 41, 29, 19, 17, 13, 7, 1, 0,
59, 47, 43, 37, 29, 19, 11, 1, 0,
53, 47, 41, 37, 29, 23, 19, 11, 7, 0,
53, 43, 31, 29, 19, 17, 13, 11, 1, 0,
47, 43, 41, 37, 23, 19, 17, 13, 0,
53, 49, 41, 31, 23, 19, 13, 7, 0,
53, 49, 43, 37, 29, 23, 11, 7, 1, 0,
59, 49, 47, 41, 31, 29, 19, 17, 11, 0,
53, 49, 47, 43, 37, 31, 23, 17, 11, 0,
49, 47, 37, 29, 19, 13, 7, 0,
59, 47, 43, 37, 31, 29, 23, 11, 1, 0,
43, 41, 37, 31, 29, 23, 19, 17, 13, 1, 0,
59, 53, 19, 13, 7, 0,
59, 41, 31, 17, 13, 11, 7, 1, 0,
53, 49, 47, 41, 29, 19, 17, 13, 11, 7, 1, 0,
49, 37, 29, 17, 11, 7, 1, 0,
53, 49, 41, 37, 23, 19, 13, 11, 7, 1, 0,
59, 53, 43, 23, 17, 13, 11, 0,
47, 43, 41, 31, 19, 17, 7, 0,
59, 53, 47, 41, 23, 17, 13, 11, 0,
59, 49, 37, 29, 13, 0,
59, 53, 49, 47, 43, 41, 31, 19, 13, 7, 0,
49, 47, 31, 29, 7, 1, 0,
53, 47, 43, 41, 37, 29, 23, 13, 11, 1, 0,
47, 37, 31, 19, 17, 13, 11, 1, 0,
49, 47, 31, 29, 23, 19, 17, 7, 0,
59, 41, 17, 13, 7, 0,
59, 43, 41, 37, 29, 13, 11, 7, 0,
59, 53, 47, 43, 31, 29, 7, 1, 0,
59, 53, 49, 43, 29, 23, 19, 17, 11, 7, 1, 0,
49, 43, 37, 23, 19, 13, 1, 0,
53, 47, 43};


 /* baby-steps, giant-steps for b1 = 333, w = 60
 q0 = 360
 */
static int num_b1_333_steps = 1110;
static uint8_t b1_333_steps[1110] = {
23, 13, 11, 7, 1, 0,
53, 47, 41, 37, 31, 23, 19, 11, 1, 0,
59, 49, 47, 41, 37, 31, 23, 19, 17, 13, 1, 0,
53, 49, 41, 37, 31, 19, 17, 0,
59, 53, 43, 37, 31, 29, 23, 13, 7, 1, 0,
59, 53, 47, 43, 41, 29, 19, 17, 13, 7, 1, 0,
59, 47, 43, 37, 29, 19, 11, 1, 0,
53, 47, 41, 37, 29, 23, 19, 11, 7, 0,
53, 43, 31, 29, 19, 17, 13, 11, 1, 0,
47, 43, 41, 37, 23, 19, 17, 13, 0,
53, 49, 41, 31, 23, 19, 13, 7, 0,
53, 49, 43, 37, 29, 23, 11, 7, 1, 0,
59, 49, 47, 41, 31, 29, 19, 17, 11, 0,
53, 49, 47, 43, 37, 31, 23, 17, 11, 0,
49, 47, 37, 29, 19, 13, 7, 0,
59, 47, 43, 37, 31, 29, 23, 11, 1, 0,
43, 41, 37, 31, 29, 23, 19, 17, 13, 1, 0,
59, 53, 19, 13, 7, 0,
59, 41, 31, 17, 13, 11, 7, 1, 0,
53, 49, 47, 41, 29, 19, 17, 13, 11, 7, 1, 0,
49, 37, 29, 17, 11, 7, 1, 0,
53, 49, 41, 37, 23, 19, 13, 11, 7, 1, 0,
59, 53, 43, 23, 17, 13, 11, 0,
47, 43, 41, 31, 19, 17, 7, 0,
59, 53, 47, 41, 23, 17, 13, 11, 0,
59, 49, 37, 29, 13, 0,
59, 53, 49, 47, 43, 41, 31, 19, 13, 7, 0,
49, 47, 31, 29, 7, 1, 0,
53, 47, 43, 41, 37, 29, 23, 13, 11, 1, 0,
47, 37, 31, 19, 17, 13, 11, 1, 0,
49, 47, 31, 29, 23, 19, 17, 7, 0,
59, 41, 17, 13, 7, 0,
59, 43, 41, 37, 29, 13, 11, 7, 0,
59, 53, 47, 43, 31, 29, 7, 1, 0,
59, 53, 49, 43, 29, 23, 19, 17, 11, 7, 1, 0,
49, 43, 37, 23, 19, 13, 1, 0,
53, 47, 43, 17, 0,
59, 49, 41, 37, 31, 29, 23, 1, 0,
49, 47, 31, 23, 19, 7, 0,
53, 43, 41, 37, 29, 23, 17, 13, 11, 7, 1, 0,
53, 49, 47, 41, 31, 29, 19, 11, 7, 0,
53, 43, 31, 29, 23, 19, 17, 1, 0,
47, 43, 37, 29, 23, 19, 1, 0,
53, 43, 37, 31, 23, 13, 1, 0,
47, 43, 37, 31, 29, 1, 0,
59, 49, 41, 37, 23, 19, 11, 0,
59, 53, 41, 37, 31, 11, 1, 0,
59, 43, 17, 13, 11, 0,
59, 53, 49, 37, 31, 23, 19, 11, 0,
49, 47, 43, 41, 29, 1, 0,
59, 53, 47, 41, 37, 31, 29, 17, 13, 1, 0,
59, 49, 47, 31, 29, 13, 7, 0,
47, 31, 23, 19, 17, 13, 11, 0,
49, 41, 29, 23, 13, 11, 7, 1, 0,
59, 53, 43, 41, 29, 19, 17, 7, 0,
53, 47, 43, 37, 29, 23, 17, 1, 0,
49, 47, 43, 29, 23, 19, 11, 1, 0,
53, 47, 41, 19, 13, 11, 1, 0,
47, 43, 37, 19, 17, 7, 0,
53, 49, 47, 37, 23, 19, 11, 0,
53, 49, 43, 41, 37, 31, 29, 17, 13, 0,
53, 31, 19, 17, 13, 7, 1, 0,
59, 53, 31, 29, 23, 7, 1, 0,
49, 47, 41, 29, 13, 11, 7, 1, 0,
47, 43, 41, 23, 0,
59, 49, 43, 41, 31, 29, 19, 17, 7, 1, 0,
59, 49, 47, 37, 31, 23, 0,
53, 43, 41, 31, 23, 17, 7, 0,
49, 43, 31, 19, 17, 0,
59, 53, 49, 43, 37, 19, 17, 7, 0,
53, 47, 43, 41, 37, 13, 11, 0,
59, 53, 37, 29, 23, 17, 0,
59, 43, 41, 37, 31, 29, 23, 17, 7, 1, 0,
49, 37, 19, 17, 11, 7, 0,
49, 41, 17, 13, 11, 7, 1, 0,
59, 47, 43, 29, 0,
59, 49, 43, 31, 17, 11, 1, 0,
49, 47, 43, 37, 29, 23, 13, 11, 7, 0,
53, 47, 41, 37, 31, 29, 19, 17, 1, 0,
49, 41, 23, 19, 13, 1, 0,
59, 53, 47, 41, 13, 7, 0,
53, 49, 41, 31, 23, 11, 0,
53, 49, 47, 43, 19, 7, 1, 0,
59, 43, 37, 31, 17, 7, 0,
53, 49, 19, 13, 7, 1, 0,
53, 47, 43, 41, 29, 23, 19, 17, 11, 0,
49, 43, 41, 37, 19, 17, 13, 1, 0,
59, 53, 49, 23, 17, 11, 7, 0,
59, 49, 17, 1, 0,
59, 53, 49, 47, 43, 41, 31, 17, 11, 7, 0,
59, 49, 43, 23, 19, 17, 11, 0,
41, 37, 29, 19, 13, 7, 0,
59, 53, 41, 37, 31, 29, 23, 19, 13, 11, 1, 0,
59, 43, 37, 17, 13, 1, 0,
47, 19, 13, 0,
53, 49, 31, 23, 17, 13, 7, 0,
53, 47, 41, 31, 29, 19, 7, 0,
59, 49, 47, 37, 29, 17, 7, 0,
43, 41, 37, 29, 23, 19, 11, 0,
53, 43, 37, 31, 29, 23, 13, 1, 0,
59, 49, 43, 37, 31, 23, 17, 7, 1, 0,
59, 53, 47, 41, 31, 23, 0,
59, 53, 31, 29, 11, 7, 0,
59, 49, 19, 11, 0,
53, 49, 47, 37, 31, 29, 23, 19, 1, 0,
53, 41, 23, 7, 1, 0,
59, 47, 41, 31, 29, 19, 17, 11, 1, 0,
47, 43, 19, 17, 1, 0,
59, 49, 47, 37, 17, 13, 11, 7, 0,
59, 43, 37, 31, 29, 17, 1, 0,
53, 49, 43, 13, 11, 1, 0,
59, 53, 49, 43, 37, 29, 23, 19, 7, 1, 0,
53, 41, 37, 23, 11, 1, 0,
37, 31, 19, 13, 11, 0,
49, 41, 23, 13, 7, 0,
53, 49, 47, 41, 31, 23, 17, 13, 7, 0,
37, 23, 13, 11, 0,
59, 49, 47, 31, 29, 11, 0,
47, 29, 23, 7, 0,
49, 43, 41, 23, 19, 13, 11, 1, 0,
53, 43, 37, 31, 23, 19, 13, 11, 1, 0,
59, 47, 43, 37, 31, 29, 17, 13, 0,
59, 41, 37, 31, 11, 7, 0,
59, 53, 49, 41, 37, 23, 17, 13, 0,
59, 47, 43, 41, 11, 7, 0,
43, 37, 31, 19, 7, 0,
53, 47, 43, 41, 37, 19, 13, 1, 0,
53, 47, 43, 31, 29, 17, 0,
47, 31, 29, 23, 1, 0,
47, 41, 31, 19, 13, 11, 7, 0,
59, 49, 43, 37, 13, 0,
59, 53, 49, 41, 29, 11, 1, 0,
59, 49, 47, 43, 37, 17, 11, 7, 0,
53, 49, 47, 43, 29, 23};


static uint64_t upm1_stage2(uint64_t P, uint64_t rho, uint64_t n, uint32_t b1, uint64_t unityval)
{
    int w = 60;
    uint64_t d[32], six, five, pw, pgiant;
    int i, j;
    uint64_t acc;
    uint32_t b2 = 25 * b1;

    d[1] = P;
    d[2] = upm1_sqrredc(P, n, rho);
    six = upm1_mulredc(P, d[2], n, rho);       
    six = upm1_sqrredc(six, n, rho);          // P^6

    // 1, 7, 13, 19, 25, 31, 37, 43
    // unnecessary powers will be mapped to scratch d[0].
    j = 1;
    while ((j + 6) < 48)
    {
        d[upm1_map[j+6]] = upm1_mulredc(d[upm1_map[j]], six, n, rho);
        j += 6;
    }

    // 11, 17, 23, 29, 35, 41, 47, 53, 59
    // unnecessary powers will be mapped to scratch d[0].
    five = upm1_sqrredc(d[2], n, rho);
    five = upm1_mulredc(five, d[1], n, rho);
    d[upm1_map[11]] = upm1_mulredc(five, six, n, rho);
    j = 11;
    while ((j + 6) < 60)
    {
        d[upm1_map[j + 6]] = upm1_mulredc(d[upm1_map[j]], six, n, rho);
        j += 6;
    }

#if 0
    // baby-steps, giant-steps
    printf("/* baby-steps, giant-steps for b1 = %u, w = 60\n", b1);
    i = 1;
    uint32_t p = primes[i];
    while (p < b1) p = primes[i++];
    uint32_t q = w;
    while (q < b1) q += w;
    printf("q0 = %u\n", q);
    printf("b1_100_steps[] = {\n");
    j = 0;
    while (p < b2)
    {
        if (p < q)
        {
            printf("%u, ", q - p);
            j++;
        }
        else
        {
            q += w;
            printf("0,\n");
            printf("%u, ", q - p);
            j += 2;
        }
        p = primes[i++];
    }
    printf("};\n %d entries */ \n", j);
    
    b1 = 333;
    b2 = b1 * 25;
    printf("/* baby-steps, giant-steps for b1 = %u, w = 60\n", b1);
    i = 1;
    p = primes[i];
    while (p < b1) p = primes[i++];
    q = w;
    while (q < b1) q += w;
    printf("q0 = %u\n", q);
    printf("b1_333_steps[] = {\n");
    j = 0;
    while (p < b2)
    {
        if (p < q)
        {
            printf("%u, ", q - p);
            j++;
        }
        else
        {
            q += w;
            printf("0,\n");
            printf("%u, ", q - p);
            j += 2;
        }
        p = primes[i++];
    }
    printf("};\n %d entries */ \n", j);

    exit(1);
#endif

    pw = upm1_mulredc(d[upm1_map[59]], P, n, rho);      // assumes w=60
    pgiant = pw;
    i = w;
    while (i < b1)
    {
        pgiant = upm1_mulredc(pgiant, pw, n, rho);
        i += w;
    }

    acc = unityval;

    switch (b1)
    {
    case 100:
        
        for (i = 0; i < num_b1_100_steps; i++)
        {
            if (b1_100_steps[i] == 0)
            {
                pgiant = upm1_mulredc(pgiant, pw, n, rho);
                i++;
            }
            acc = upm1_mulredc(acc, upm1_submod(pgiant, d[upm1_map[b1_100_steps[i]]], n), n, rho);
        }
        break;
    case 333:
        
        for (i = 0; i < num_b1_333_steps; i++)
        {
            if (b1_333_steps[i] == 0)
            {
                pgiant = upm1_mulredc(pgiant, pw, n, rho);
                i++;
            }
            acc = upm1_mulredc(acc, upm1_submod(pgiant, d[upm1_map[b1_333_steps[i]]], n), n, rho);
        }
        break;
    }

    return acc;
}

static uint64_t micropm1(uint64_t n, uint32_t B1, uint32_t B2)
{
    //attempt to factor n with the elliptic curve method
    //following brent and montgomery's papers, and CP's book
    int result;
    uint64_t stg1_res, f = 0, q;
    uint64_t tmp1;
    uint64_t rho = (uint64_t)0 - upm1_multiplicative_inverse(n);

    // Let R = 2^64.  We can see R%n ≡ (R-n)%n  (mod n)
    uint64_t unityval = ((uint64_t)0 - n) % n;   // unityval == R  (mod n)

    stg1_res = upm1_stage1(rho, n, unityval, B1);
    q = upm1_mulredc(1, stg1_res, n, rho);
    result = upm1_check_factor(q-1, n, &tmp1);

    if (result == 1)
    {
        f = tmp1;
    }
    else if (B2 > B1)
    {
        uint64_t stg2acc = upm1_stage2(stg1_res, rho, n, B1, unityval);
        q = upm1_mulredc(1, stg2acc, n, rho);

        result = upm1_check_factor(q, n, &tmp1);

        if (result == 1)
        {
            f = tmp1;
        }
    }

    return f;
}

static uint64_t upm1_dispatch(uint64_t n, int targetBits, uint32_t b1)
{
    int B1 = 333;

    //upm1_generate_window_plan(100);
    //upm1_generate_window_plan(333);
    //upm1_generate_window_plan(666);
    //upm1_generate_window_plan(1000);
    //exit(0);

    if (b1 > 0)
    {
        return micropm1(n, b1, 25 * b1);
    }
    else
    {
        if (targetBits < 50) B1 = 100;
        else B1 = 333;

        return micropm1(n, B1, 25 * B1);
    }
}

static int upm1_get_bits(uint64_t n)
{
    int i = 0;
    while (n != 0)
    {
        n >>= 1;
        i++;
    }
    return i;
}


// getfactor_upm1() returns 1 if unable to find a factor of q64,
// Otherwise it returns a factor of q64.
// 
// if the input is known to have no small factors, set is_arbitrary=0, 
// otherwise, set is_arbitrary=1 and a few curves targetting small factors
// will be run prior to the standard sequence of curves for the input size.
//  
// Prior to your first call of getfactor_upm1(), set *pran = 0  (or set it to
// some other arbitrary value); after that, don't change *pran.
// FYI: *pran is used within this file by a random number generator, and it
// holds the current value of a pseudo random sequence.  Your first assigment
// to *pran seeds the sequence, and after seeding it you don't want to
// change *pran, since that would restart the sequence.
uint64_t getfactor_upm1(uint64_t q64, uint32_t b1)
{
    if (q64 % 2 == 0)
        return 2;

    if (b1 > 0)
    {
        return upm1_dispatch(q64, 0, b1);
    }
    else
    {
        int bits = upm1_get_bits(q64);
        return upm1_dispatch(q64, bits, 0);
    }
}


