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
*/

#include "gmp.h"
#include <stdint.h>
#include <stdio.h>
#include "microecm.h"
#include "if.h"
#include "gmp-aux.h"

#define D 120

#if defined( _MSC_VER) || defined(AVX512_ECM)
#define USE_AVX512F
#endif

//#define DEBUG 1
#define base_t uint64_t
#define base_digits 2

#ifdef USE_AVX512F
#include <immintrin.h>
#endif


typedef struct
{
	uint64_t base[2];
} u128_t;

typedef struct
{
	uint64_t data[2][8];
} vec_u104_t;

typedef struct
{
	base_t X[base_digits];
	base_t Z[base_digits];
} tinyecm_pt;

typedef struct
{
	vec_u104_t X;
	vec_u104_t Z;
} tecm_pt_x8;

typedef struct
{
	base_t sum1[base_digits];
	base_t diff1[base_digits];
	base_t sum2[base_digits];
	base_t diff2[base_digits];
	base_t tt1[base_digits];
	base_t tt2[base_digits];
	base_t tt3[base_digits];
	base_t tt4[base_digits];
	base_t tt5[base_digits];
	base_t s[base_digits];
	base_t n[base_digits];
	tinyecm_pt pt1;
	tinyecm_pt pt2;
	tinyecm_pt pt3;
	tinyecm_pt pt4;
	tinyecm_pt pt5;
	uint32_t sigma;

	tinyecm_pt Pa;
	tinyecm_pt Pd;
	tinyecm_pt Pad;
	tinyecm_pt Pb[20];
	base_t Paprod[base_digits];
	base_t Pbprod[20][base_digits];
	
	base_t stg2acc[base_digits];
	uint32_t stg1Add;
	uint32_t stg1Doub;
	uint32_t paired;
	uint32_t ptadds;
	uint64_t numprimes;
	uint64_t A;
	uint32_t last_pid;
	uint32_t amin;

	uint32_t stg1_max;
	uint32_t stg2_max;

} tinyecm_work;

typedef struct
{
	vec_u104_t sum1;
	vec_u104_t diff1;
	vec_u104_t sum2;
	vec_u104_t diff2;
	vec_u104_t tt1;
	vec_u104_t tt2;
	vec_u104_t tt3;
	vec_u104_t tt4;
	vec_u104_t tt5;
	vec_u104_t s;
	vec_u104_t n;
	tecm_pt_x8 pt1;
	tecm_pt_x8 pt2;
	tecm_pt_x8 pt3;
	tecm_pt_x8 pt4;
	tecm_pt_x8 pt5;
	uint32_t sigma[8];

	tecm_pt_x8 Pa;
	tecm_pt_x8 Pd;
	tecm_pt_x8 Pad;
	tecm_pt_x8 Pb[20];
	vec_u104_t Paprod;
	vec_u104_t Pbprod[20];

	vec_u104_t stg2acc;
	uint32_t stg1Add;
	uint32_t stg1Doub;
	uint32_t paired;
	uint32_t ptadds;
	uint64_t numprimes;
	uint64_t A;
	uint32_t last_pid;
	uint32_t amin;

	uint32_t stg1_max;
	uint32_t stg2_max;

} tinyecm_work_x8;

static const uint32_t map[60] = {
	0, 1, 2, 0, 0, 0, 0, 3, 0, 0,
	0, 4, 0, 5, 0, 0, 0, 6, 0, 7,
	0, 0, 0, 8, 0, 0, 0, 0, 0, 9,
	0, 10, 0, 0, 0, 0, 0, 11, 0, 0,
	0, 12, 0, 13, 0, 0, 0, 14, 0, 15,
	0, 0, 0, 16, 0, 0, 0, 0, 0, 17 };

#define NUMP 801
static const int tecm_primes[NUMP] = {
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

/********************* 128-bit Montgomery arith **********************/
typedef struct
{
    uint64_t r[2];
    uint64_t n[2];
    uint64_t np[2];
    uint64_t nhat[2];
    uint64_t rhat[2];
    uint64_t rmask[2];
    uint64_t one[2];
    uint64_t mtmp1[2];
    uint64_t mtmp2[2];
    uint64_t mtmp3[2];
    uint64_t mtmp4[2];
    uint64_t rho;
} monty128_t;

#define GCC_ASM64X
#define USE_MULX

void to_monty128(monty128_t* mdata, uint64_t* x);
void monty128_init(monty128_t* mdata, uint64_t* n);
void mulmod128(uint64_t* u, uint64_t* v, uint64_t* w, monty128_t* mdata);
void sqrmod128(uint64_t* u, uint64_t* w, monty128_t* mdata);
void addmod128(uint64_t* u, uint64_t* v, uint64_t* w, uint64_t* n);
void submod128(uint64_t* u, uint64_t* v, uint64_t* w, uint64_t* n);
__inline uint64_t _umul128(uint64_t x, uint64_t y, uint64_t* hi);

#ifdef GCC_ASM64X
__inline uint64_t _umul128(uint64_t x, uint64_t y, uint64_t* hi)
{
    __asm__(
        "mulq %3	\n\t"
        : "=&a"(x), "=&d"(y)
        : "0"(x), "1"(y)
        : "cc"
    );

    *hi = y;
    return x;
}
#endif

/********************* 128-bit Montgomery arith **********************/
void to_monty128(monty128_t* mdata, uint64_t* x)
{
    //given a number x in normal (hexadecimal) representation, 
    //find its montgomery representation

    //this uses some precomputed monty constants
    //xhat = (x * r) mod n
    // = x * R^2 / R mod n
    // = REDC(x * R^2)
    mpz_t m;
    mpz_t n;

    mpz_init(m);
    mpz_init(n);

    mpz_set_ull(m, x[1]);
    mpz_mul_2exp(m, m, 64);
    mpz_add_ull(m, m, x[0]);

    mpz_set_ull(n, mdata->n[1]);
    mpz_mul_2exp(n, n, 64);
    mpz_add_ull(n, n, mdata->n[0]);

    // implied R = 2^128
    mpz_mul_2exp(m, m, 128);
    mpz_mod(m, m, n);

    x[0] = mpz_get_ull(m);
    mpz_tdiv_q_2exp(m, m, 64);
    x[1] = mpz_get_ull(m);

    mpz_clear(m);
    mpz_clear(n);

    return;
}

void monty128_init(monty128_t* mdata, uint64_t* n)
{
    //for a input modulus n, initialize constants for 
    //montogomery representation
    //this assumes that n is relatively prime to 2, i.e. is odd.	
    uint64_t b = n[0];
    uint64_t x;

    mdata->n[0] = n[0];
    mdata->n[1] = n[1];

    // invert (odd) n mod 2^64
    x = (((b + 2) & 4) << 1) + b; // here x*a==1 mod 2**4
    x *= 2 - b * x;               // here x*a==1 mod 2**8
    x *= 2 - b * x;               // here x*a==1 mod 2**16
    x *= 2 - b * x;               // here x*a==1 mod 2**32         
    x *= 2 - b * x;               // here x*a==1 mod 2**64

    mdata->rho = (uint64_t)((uint64_t)0 - ((uint64_t)x));

    mdata->one[0] = 1;
    mdata->one[1] = 0;
    to_monty128(mdata, mdata->one);

    return;
}

void ciosFullMul128x(uint64_t* u, uint64_t* v, uint64_t rho, uint64_t* n, uint64_t* w)
{
#if defined( USE_MULX ) && defined(GCC_ASM64X)
    // requires mulx in BMI2 (via the AVX2 macro) and GCC_ASM64 syntax

    __asm__(
        "movq %0, %%r10	\n\t"			/* u ptr in r10 */
        "movq %2, %%r11	\n\t"			/* w ptr in r11 */
        "movq 0(%1), %%r9	\n\t"		/* v[0] ptr in r9 */

        /* begin s += u * v */
        "movq 0(%%r10), %%rdx	\n\t"   /* ready to multiply by u[0]  */
        "mulx %%r9, %%r12, %%r14 \n\t"  /* r14 = HI(u[0] * v)         */
        "addq %%r12, 0(%%r11) \n\t"     /* w[0] = w[0] + LO(u[0] * v) */

        "movq 8(%%r10), %%rdx	\n\t"   /* ready to multiply by u[1]  */
        "mulx %%r9, %%r12, %%r13 \n\t"  /* r13 = HI(u[1] * v)         */
        "adcq %%r14, %%r12 \n\t"        /* r12 = HI(u[0] * v) + LO(u[1] * v) + prevcarry */
        "adcq $0, %%r13 \n\t"           /* r13 = HI(u[1] * v) + prevcarry                */
        "addq %%r12, 8(%%r11) \n\t"		/* w[1] = w[1] + HI(u[0] * v) + LO(u[1] * v)*/
        "adcq %%r13, 16(%%r11) \n\t"	/* w[2] = w[2] + HI(u[1] * v) + prevcarry */

        "movq 0(%%r11), %%rdx	\n\t"   /* ready to multiply by w[0]  */
        "mulx %4, %%r9, %%r14	\n\t"   /* m = rho * w[0]         */
        "movq %3, %%r10	\n\t"			/* n ptr in r10 */

        /* begin s = (s + n * m) >> 64 */
        "movq 0(%%r10), %%rdx	\n\t"   /* ready to multiply by n[0]  */
        "mulx %%r9, %%r12, %%r14 \n\t"  /* r14 = HI(n[0] * m)         */
        "addq 0(%%r11), %%r12  \n\t"    /* r12 = w[0] (could be rdx) + LO(n[0] * m) */

        /* r12 should be 0 here */

        "movq 8(%%r10), %%rdx	\n\t"   /* ready to multiply by n[1]  */
        "mulx %%r9, %%r12, %%r13 \n\t"  /* r13 = HI(n[1] * m)         */
        "adcq %%r14, %%r12 \n\t"        /* r12 = HI(n[0] * m) + LO(n[1] * m) + prevcarry */
        "adcq $0, %%r13 \n\t"           /* r13 = HI(n[1] * m) + prevcarry                */
        "xorq %%r14, %%r14 \n\t"
        "addq 8(%%r11), %%r12  \n\t"
        "movq %%r12, 0(%%r11)	\n\t"	/* w[0] = w[1] + HI(n[0] * m) + LO(n[1] * m) + prevcarry */

        "adcq 16(%%r11), %%r13 \n\t"
        "movq %%r13, 8(%%r11)	\n\t"	/* w[1] = w[2] + HI(n[1] * m) + prevcarry */
        "adcq $0, %%r14 \n\t"
        "movq %%r14, 16(%%r11)	\n\t"   /* w[2] = carry out */

        /* round 2 */

        "movq %0, %%r10	\n\t"			/* u ptr in r10 */
        "movq 8(%1), %%r9	\n\t"		/* v[1] ptr in r9 */

        /* begin s += u * v */
        "movq 0(%%r10), %%rdx	\n\t"   /* ready to multiply by u[0]  */
        "mulx %%r9, %%r12, %%r14 \n\t"  /* r14 = HI(u[0] * v)         */
        "addq %%r12, 0(%%r11) \n\t"     /* w[0] = w[0] + LO(u[0] * v) */

        "movq 8(%%r10), %%rdx	\n\t"   /* ready to multiply by u[1]  */
        "mulx %%r9, %%r12, %%r13 \n\t"  /* r13 = HI(u[1] * v)         */
        "adcq %%r14, %%r12 \n\t"        /* r12 = HI(u[0] * v) + LO(u[1] * v) + prevcarry */
        "adcq $0, %%r13 \n\t"           /* r13 = HI(u[1] * v) + prevcarry                */
        "addq %%r12, 8(%%r11) \n\t"		/* w[1] = w[1] + HI(u[0] * v) + LO(u[1] * v)*/
        "adcq %%r13, 16(%%r11) \n\t"	/* w[2] = w[2] + HI(u[1] * v) + prevcarry */

        "movq 0(%%r11), %%rdx	\n\t"   /* ready to multiply by w[0]  */
        "mulx %4, %4, %%r14	\n\t"       /* m = rho * w[0]         */
        "movq %3, %%r10	\n\t"			/* n ptr in r10 */

        /* begin s = (s + n * m) >> 64 */
        "movq 0(%%r10), %%rdx	\n\t"   /* ready to multiply by n[0]  */
        "mulx %4, %%r12, %%r14	\n\t"   /* r14 = HI(n[0] * m)         */
        "addq 0(%%r11), %%r12  \n\t"    /* r12 = w[0] (could be rdx) + LO(n[0] * m) */

        /* r12 should be 0 here */

        "movq 8(%%r10), %%rdx	\n\t"   /* ready to multiply by n[1]  */
        "mulx %4, %%r12, %%r13	\n\t"   /* r13 = HI(n[1] * m)         */
        "adcq %%r14, %%r12 \n\t"        /* r12 = HI(n[0] * m) + LO(n[1] * m) + prevcarry */
        "adcq $0, %%r13 \n\t"           /* r13 = HI(n[1] * m) + prevcarry                */
        "xorq %%r14, %%r14 \n\t"
        "addq 8(%%r11), %%r12  \n\t"
        "movq %%r12, 0(%%r11)	\n\t"	/* w[0] = w[1] + HI(n[0] * m) + LO(n[1] * m) + prevcarry */

        "adcq 16(%%r11), %%r13 \n\t"
        "movq %%r13, 8(%%r11)	\n\t"	/* w[1] = w[2] + HI(n[1] * m) + prevcarry */
        "adcq $0, %%r14 \n\t"
        "movq %%r14, 16(%%r11)	\n\t"   /* w[2] = carry out */



        :
    : "r"(u), "r"(v), "r"(w), "r"(n), "r"(rho)
        : "r9", "r10", "rdx", "r11", "r12", "r13", "r14", "cc", "memory");

#else


#endif
    return;
}

// already defined within mingw64/msys2
#if 0 //defined( GCC_ASM64X ) && !defined(__MINGW32__) && !defined(__INTEL_COMPILER) && !defined(__INTEL_LLVM_COMPILER)

__inline uint8_t _addcarry_u64(uint64_t x, uint8_t w, uint64_t y, uint64_t* sum)
{
    uint64_t s, c;
    s = y;
    c = 0;

    __asm__("movq %2, %%rax		\n\t"
        "addq %3, %%rax		\n\t"
        "adcq $0, %5		\n\t"
        "addq %%rax, %4		\n\t"
        "adcq $0, %5		\n\t"
        : "=r"(s), "=r"(c)
        : "r"(x), "r"((uint64_t)w), "0"(s), "1"(c)
        : "rax", "memory", "cc");

    *sum = s;
    return c;
}

#endif

void mulmod128(uint64_t* u, uint64_t* v, uint64_t* w, monty128_t* mdata)
{
    // integrate multiply and reduction steps, alternating
    // between iterations of the outer loops.
    uint64_t s[3];

    s[0] = 0;
    s[1] = 0;
    s[2] = 0;

#if defined( USE_MULX ) && defined(GCC_ASM64X)
    // requires mulx in BMI2 (via the AVX2 macro) and GCC_ASM64 syntax
    ciosFullMul128x(u, v, mdata->rho, mdata->n, s);

    if ((s[2]) || (s[1] > mdata->n[1]) || ((s[1] == mdata->n[1]) && (s[0] > mdata->n[0])))
    {
        __asm__(
            "movq %4, %%r11 \n\t"
            "movq %0, 0(%%r11) \n\t"
            "movq %1, 8(%%r11) \n\t"
            "subq %2, 0(%%r11) \n\t"
            "sbbq %3, 8(%%r11) \n\t"
            :
        : "r"(s[0]), "r"(s[1]), "r"(mdata->n[0]), "r"(mdata->n[1]), "r"(w)
            : "r11", "cc", "memory");
    }
    else
    {
        w[0] = s[0];
        w[1] = s[1];
    }
#else
    // TODO: implement portable u128 x u128 modular multiplication
    uint64_t t[3], U, c2, c3;
    uint8_t c1, c4;

    // z = 0
    // for (i = 0; i < t; i++)
    // {
    //     u = (z0 + xi * y0) * -m’ mod b
    //     z = (z + xi * y + u * m) / b
    // }
    // if (x >= m) { x -= m; }

    //printf("nhat = %" PRIx64 "\n", mdata->rho);
    //printf("a = %" PRIx64 ",%" PRIx64 "\n", u[1], u[0]);
    //printf("b = %" PRIx64 ",%" PRIx64 "\n", v[1], v[0]);
    //printf("n = %" PRIx64 ",%" PRIx64 "\n", mdata->n[1], mdata->n[0]);

    //i = 0;
    //u = (z0 + x0 * y0) * nhat mod b;
    t[0] = _umul128(u[0], v[0], &t[1]);
    U = t[0] * mdata->rho;

    //printf("t[0] = %" PRIx64 ", t[1] = %" PRIx64 ", u = %" PRIx64 "\n", t[0], t[1], U);

    //j = 0;
    //z = (z + x0 * y0 + u * m0) / b;
    c1 = _addcarry_u64(0, _umul128(U, mdata->n[0], &c2), t[0], &t[0]);      // c1,c2 apply to t1
    t[2] = _addcarry_u64(c1, t[1], c2, &t[1]);                              // c1 applies to t2

    //j = 1;
    //z = (z + x0 * y1 + u * m1) / b;
    c1 = _addcarry_u64(0, _umul128(u[0], v[1], &c2), t[1], &t[1]);          // c1,c2 apply to t1
    c4 = _addcarry_u64(0, _umul128(U, mdata->n[1], &c3), t[1], &t[1]);      // c1,c2 apply to t1
    c1 = _addcarry_u64(c1, t[2], c2, &t[2]);                                // c1 applies to t2
    c4 = _addcarry_u64(c4, t[2], c3, &t[2]);                                // c4 applies to t2
    t[0] = t[1];                                                            // divide by b
    t[1] = t[2];
    t[2] = c1 + c4;

    //printf("t = %" PRIx64 ",%" PRIx64 ",%" PRIx64 "\n", t[2], t[1], t[0]);

    //i = 0;
    //u = (z + x1 * y0) * nhat mod b;
    c1 = _addcarry_u64(0, _umul128(u[1], v[0], &c2), t[0], &t[0]);          // c1,c2 apply to t1
    c1 = _addcarry_u64(c1, t[1], c2, &t[1]);                                // c1 applies to t2
    t[2] += c1;
    U = t[0] * mdata->rho;

    //printf("t[0] = %" PRIx64 ", u = %" PRIx64 "\n", t[0], U);

    //j = 0;
    //z = (z + x1 * y0 + u * m0) / b;
    c1 = _addcarry_u64(0, _umul128(U, mdata->n[0], &c2), t[0], &t[0]);      // c1,c2 apply to t1
    c1 = _addcarry_u64(c1, t[1], c2, &t[1]);                                // c1 applies to t2
    t[2] += c1;

    //j = 1;
    //z = (z + x1 * y1 + u * m1) / b;
    c1 = _addcarry_u64(0, _umul128(u[1], v[1], &c2), t[1], &t[1]);          // c1,c2 apply to t1
    c4 = _addcarry_u64(0, _umul128(U, mdata->n[1], &c3), t[1], &t[1]);      // c1,c2 apply to t1
    c1 = _addcarry_u64(c1, t[2], c2, &t[2]);                                // c1 applies to t2
    c4 = _addcarry_u64(c4, t[2], c3, &t[2]);                                // c4 applies to t2

    //printf("t = %" PRIx64 ",%" PRIx64 ",%" PRIx64 "\n", t[2], t[1], t[0]);

    w[0] = t[1];
    w[1] = t[2];

    //exit(1);

    return;



#endif

    return;
}

void sqrmod128(uint64_t* u, uint64_t* w, monty128_t* mdata)
{
    // integrate multiply and reduction steps, alternating
    // between iterations of the outer loops.
    uint64_t s[3];

    s[0] = 0;
    s[1] = 0;
    s[2] = 0;

#if defined( USE_MULX ) && defined(GCC_ASM64X)
    // requires mulx in BMI2 (via the AVX2 macro) and GCC_ASM64 syntax
    ciosFullMul128x(u, u, mdata->rho, mdata->n, s);

    if ((s[2]) || (s[1] > mdata->n[1]) || ((s[1] == mdata->n[1]) && (s[0] > mdata->n[0])))
    {
        __asm__(
            "movq %4, %%r11 \n\t"
            "movq %0, 0(%%r11) \n\t"
            "movq %1, 8(%%r11) \n\t"
            "subq %2, 0(%%r11) \n\t"
            "sbbq %3, 8(%%r11) \n\t"
            :
        : "r"(s[0]), "r"(s[1]), "r"(mdata->n[0]), "r"(mdata->n[1]), "r"(w)
            : "r11", "cc", "memory");
    }
    else
    {
        w[0] = s[0];
        w[1] = s[1];
    }

#else

    mulmod128(u, u, w, mdata);

#endif
    return;
}

void addmod128(uint64_t* a, uint64_t* b, uint64_t* w, uint64_t* n)
{
#if defined(GCC_ASM64X)
    // requires GCC_ASM64 syntax
    w[1] = a[1];
    w[0] = a[0];
    __asm__(
        "movq %0, %%r8 \n\t"
        "movq %1, %%r9 \n\t"
        "subq %4, %%r8 \n\t"		/* t = x - n */
        "sbbq %5, %%r9 \n\t"
        "addq %2, %0 \n\t"			/* x += y */
        "adcq %3, %1 \n\t"
        "addq %2, %%r8 \n\t"		/* t += y */
        "adcq %3, %%r9 \n\t"
        "cmovc %%r8, %0 \n\t"
        "cmovc %%r9, %1 \n\t"
        : "+r"(w[0]), "+r"(w[1])
        : "r"(b[0]), "r"(b[1]), "r"(n[0]), "r"(n[1])
        : "r8", "r9", "cc", "memory");

#else

    uint8_t c;
    uint64_t t[2];
    c = _addcarry_u64(0, a[0], b[0], &t[0]);
    c = _addcarry_u64(c, a[1], b[1], &t[1]);
    if (c || (t[1] > n[1]) || ((t[1] == n[1]) && (t[0] > n[0])))
    {
        c = _subborrow_u64(0, t[0], n[0], &w[0]);
        c = _subborrow_u64(c, t[1], n[1], &w[1]);
    }
    else
    {
        w[0] = t[0];
        w[1] = t[1];
    }


#endif
    return;
}

void submod128(uint64_t* a, uint64_t* b, uint64_t* w, uint64_t* n)
{
#if defined(GCC_ASM64X)
    // requires GCC_ASM64 syntax
    __asm__(
        "movq %6, %%r11 \n\t"
        "xorq %%r8, %%r8 \n\t"
        "xorq %%r9, %%r9 \n\t"
        "movq %0, 0(%%r11) \n\t"
        "movq %1, 8(%%r11) \n\t"
        "subq %2, 0(%%r11) \n\t"
        "sbbq %3, 8(%%r11) \n\t"
        "cmovc %4, %%r8 \n\t"
        "cmovc %5, %%r9 \n\t"
        "addq %%r8, 0(%%r11) \n\t"
        "adcq %%r9, 8(%%r11) \n\t"
        "1: \n\t"
        :
    : "r"(a[0]), "r"(a[1]), "r"(b[0]), "r"(b[1]), "r"(n[0]), "r"(n[1]), "r"(w)
        : "r8", "r9", "r11", "cc", "memory");

#else

    uint8_t c;
    uint64_t t[2];
    c = _subborrow_u64(0, a[0], b[0], &t[0]);
    c = _subborrow_u64(c, a[1], b[1], &t[1]);
    if (c)
    {
        c = _addcarry_u64(0, t[0], n[0], &w[0]);
        c = _addcarry_u64(c, t[1], n[1], &w[1]);
    }
    else
    {
        w[0] = t[0];
        w[1] = t[1];
    }

#endif

    return;
}

// local functions
void add(monty128_t *mdata, tinyecm_work *work, tinyecm_pt *P1, tinyecm_pt *P2, 
	tinyecm_pt *Pin, tinyecm_pt *Pout);
void duplicate(monty128_t *mdata, tinyecm_work *work, uint64_t * insum, uint64_t * indiff, tinyecm_pt *P);
void prac(monty128_t *mdata, tinyecm_work *work, tinyecm_pt *P, uint64_t c, double v);
int check_factor(uint64_t * Z, uint64_t * n, uint64_t * f);
void build_one_curve(tinyecm_pt *P, monty128_t *mdata, 
	tinyecm_work *work, uint32_t sigma, uint64_t * lcg_state, int verbose);

void ecm_stage1(monty128_t *mdata, tinyecm_work *work, tinyecm_pt *P);
void ecm_stage2(tinyecm_pt *P, monty128_t *mdata, tinyecm_work *work);

double bench_curves;
double bench_stg1;
double bench_stg2;

void u128_to_mpz(uint64_t *in, mpz_t out)
{
	mpz_set_ull(out, in[1]);
	mpz_mul_2exp(out, out, 64);
	mpz_add_ull(out, out, in[0]);
	return;
}

void mpz_to_u128(mpz_t in, uint64_t *out)
{
    out[0] = mpz_get_ull(in);
    mpz_tdiv_q_2exp(in, in, 64);
    out[1] = mpz_get_ull(in);
	
    // restore input
    mpz_mul_2exp(in, in, 64);
	mpz_add_ull(in, in, out[0]);

	return;
}

__inline void copy128(uint64_t *src, uint64_t *dest)
{
	dest[0] = src[0];
	dest[1] = src[1];
	return;
}

__inline void swap128(uint64_t *a, uint64_t *b)
{
	uint64_t tmp[2];
	tmp[0] = a[0];
	tmp[1] = a[1];
	a[0] = b[0];
	a[1] = b[1];
	b[0] = tmp[0];
	b[1] = tmp[1];
	return;
}

__inline void rot128(uint64_t *a, uint64_t *b, uint64_t *c)
{
	uint64_t tmp[2];
	tmp[0] = a[0];
	tmp[1] = a[1];
	a[0] = b[0];
	a[1] = b[1];
	b[0] = c[0];
	b[1] = c[1];
	c[0] = tmp[0];
	c[1] = tmp[1];
	return;
}

void tinyecm_work_init(tinyecm_work *work)
{
	work->stg1Add = 0;
	work->stg1Doub = 0;

	return;
}

void tinyecm_work_free(tinyecm_work *work)
{
	return;
}

void add(monty128_t *mdata, tinyecm_work *work, tinyecm_pt *P1, tinyecm_pt *P2, 
	tinyecm_pt *Pin, tinyecm_pt *Pout)
{
	// compute:
	//x+ = z- * [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
	//z+ = x- * [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
	// where:
	//x- = original x
	//z- = original z
	addmod128(P1->X, P1->Z, work->sum1, work->n);
	submod128(P1->X, P1->Z, work->diff1, work->n);
	addmod128(P2->X, P2->Z, work->sum2, work->n);
	submod128(P2->X, P2->Z, work->diff2, work->n);

	mulmod128(work->diff1, work->sum2, work->tt1, mdata);	//U
	mulmod128(work->sum1, work->diff2, work->tt2, mdata);	//V

	addmod128(work->tt1, work->tt2, work->tt3, work->n);
	submod128(work->tt1, work->tt2, work->tt4, work->n);
	sqrmod128(work->tt3, work->tt1, mdata);					//(U + V)^2
	sqrmod128(work->tt4, work->tt2, mdata);					//(U - V)^2

	// choosing the initial point Pz0 = 1 means that z_p-q = 1 and this mul isn't necessary...
	// but that involves a different way to initialize curves, so for now
	// we can't assume Z=1
	if (Pin->X == Pout->X)
	{
		mulmod128(work->tt1, Pin->Z, Pout->Z, mdata);		//Z * (U + V)^2
		mulmod128(work->tt2, Pin->X, Pout->X, mdata);		//x * (U - V)^2
		swap128(Pout->Z, Pout->X);
	}
	else
	{
		mulmod128(work->tt1, Pin->Z, Pout->X, mdata);		//Z * (U + V)^2
		mulmod128(work->tt2, Pin->X, Pout->Z, mdata);		//x * (U - V)^2
	}
	work->stg1Add++;
	return;
}

void duplicate(monty128_t *mdata, tinyecm_work *work, 
	uint64_t * insum, uint64_t * indiff, tinyecm_pt *P)
{
	sqrmod128(indiff, work->tt1, mdata);						// U=(x1 - z1)^2
	sqrmod128(insum,  work->tt2, mdata);						// V=(x1 + z1)^2
	mulmod128(work->tt1, work->tt2, P->X, mdata);				// x=U*V

	submod128(work->tt2, work->tt1, work->tt3, mdata->n);	    // w = V-U
	mulmod128(work->tt3, work->s, work->tt2, mdata);			// w = (A+2)/4 * w
	addmod128(work->tt2, work->tt1, work->tt2, mdata->n);       // w = w + U
	mulmod128(work->tt2, work->tt3, P->Z, mdata);				// Z = w*(V-U)
	work->stg1Doub++;
	return;
}

#define ADD 6.0
#define DUP 5.0
#define NV 10  

double getEcost(uint64_t d, uint64_t e)
{
	int doub = 0, add = 0;

	while (d > 0)
	{
		if ((e / 2) < d)
		{
			d = e - d;
		}
		else if ((d < (e / 4)) && ((e & 1) == 0))
		{
			e = e / 2;
			doub++;
			add++;
		}
		else
		{
			e = e - d;
			add++;
		}

	}
	return (doub + add) * 2 + add * 4 + doub * 3;
}

static double lucas_cost(uint64_t n, double v)
{
	uint64_t d, e, r;
	double c; /* cost */

	d = n;
	r = (uint64_t)((double)d * v + 0.5);
	if (r >= n)
		return (ADD * (double)n);
	d = n - r;
	e = 2 * r - n;
	c = DUP + ADD; /* initial duplicate and final addition */
	while (d != e)
	{
		if (d < e)
		{
			r = d;
			d = e;
			e = r;
		}
		if (d - e <= e / 4 && ((d + e) % 3) == 0)
		{ /* condition 1 */
			d = (2 * d - e) / 3;
			e = (e - d) / 2;
			c += 3.0 * ADD; /* 3 additions */
		}
		else if (d - e <= e / 4 && (d - e) % 6 == 0)
		{ /* condition 2 */
			d = (d - e) / 2;
			c += ADD + DUP; /* one addition, one duplicate */
		}
		else if ((d + 3) / 4 <= e)
		{ /* condition 3 */
			d -= e;
			c += ADD; /* one addition */
		}
		else if ((d + e) % 2 == 0)
		{ /* condition 4 */
			d = (d - e) / 2;
			c += ADD + DUP; /* one addition, one duplicate */
		}
		/* now d+e is odd */
		else if (d % 2 == 0)
		{ /* condition 5 */
			d /= 2;
			c += ADD + DUP; /* one addition, one duplicate */
		}
		/* now d is odd and e is even */
		else if (d % 3 == 0)
		{ /* condition 6 */
			d = d / 3 - e;
			c += 3.0 * ADD + DUP; /* three additions, one duplicate */
		}
		else if ((d + e) % 3 == 0)
		{ /* condition 7 */
			d = (d - 2 * e) / 3;
			c += 3.0 * ADD + DUP; /* three additions, one duplicate */
		}
		else if ((d - e) % 3 == 0)
		{ /* condition 8 */
			d = (d - e) / 3;
			c += 3.0 * ADD + DUP; /* three additions, one duplicate */
		}
		else /* necessarily e is even: catches all cases */
		{ /* condition 9 */
			e /= 2;
			c += ADD + DUP; /* one addition, one duplicate */
		}
	}

	if (d != 1)
	{
		c = 9999999.;
	}

	return c;
}

void prac70(monty128_t *mdata, tinyecm_work *work, tinyecm_pt *P)
{
	uint64_t *s1, *s2, *d1, *d2;
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

	s1 = work->sum1;
	s2 = work->sum2;
	d1 = work->diff1;
	d2 = work->diff2;

	for (i = 0; i < 116; i++)
	{
		if (steps[i] == 0)
		{
			work->pt1.X[0] = work->pt2.X[0] = work->pt3.X[0] = P->X[0];
			work->pt1.X[1] = work->pt2.X[1] = work->pt3.X[1] = P->X[1];
			work->pt1.Z[0] = work->pt2.Z[0] = work->pt3.Z[0] = P->Z[0];
			work->pt1.Z[1] = work->pt2.Z[1] = work->pt3.Z[1] = P->Z[1];

			submod128(work->pt1.X, work->pt1.Z, d1, work->n);
			addmod128(work->pt1.X, work->pt1.Z, s1, work->n);

			// point2 is [2]P
			duplicate(mdata, work, s1, d1, &work->pt1);
		}
		else if (steps[i] == 3)
		{
			// integrate step 4 followed by swap(1,2)
			add(mdata, work, &work->pt2, &work->pt1, &work->pt3, &work->pt4);		// T = B + A (C)
			rot128(work->pt2.X, work->pt4.X, work->pt3.X);
			rot128(work->pt2.Z, work->pt4.Z, work->pt3.Z);
			swap128(work->pt1.X, work->pt2.X);
			swap128(work->pt1.Z, work->pt2.Z);
		}
		else if (steps[i] == 4)
		{
			add(mdata, work, &work->pt2, &work->pt1, &work->pt3, &work->pt4);		// T = B + A (C)
			rot128(work->pt2.X, work->pt4.X, work->pt3.X);
			rot128(work->pt2.Z, work->pt4.Z, work->pt3.Z);
		}
		else if (steps[i] == 5)
		{
			add(mdata, work, &work->pt2, &work->pt1, &work->pt3, &work->pt2);		// B = B + A (C)

			submod128(work->pt1.X, work->pt1.Z, d2, work->n);
			addmod128(work->pt1.X, work->pt1.Z, s2, work->n);
			duplicate(mdata, work, s2, d2, &work->pt1);		// A = 2A
		}
		else if (steps[i] == 6)
		{
			add(mdata, work, &work->pt1, &work->pt2, &work->pt3, P);		// A = A + B (C)
		}

	}

	return;

}

void prac85(monty128_t *mdata, tinyecm_work *work, tinyecm_pt *P)
{
	uint64_t *s1, *s2, *d1, *d2;
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

	s1 = work->sum1;
	s2 = work->sum2;
	d1 = work->diff1;
	d2 = work->diff2;

	for (i = 0; i < 146; i++)
	{
		if (steps[i] == 0)
		{
			work->pt1.X[0] = work->pt2.X[0] = work->pt3.X[0] = P->X[0];
			work->pt1.X[1] = work->pt2.X[1] = work->pt3.X[1] = P->X[1];
			work->pt1.Z[0] = work->pt2.Z[0] = work->pt3.Z[0] = P->Z[0];
			work->pt1.Z[1] = work->pt2.Z[1] = work->pt3.Z[1] = P->Z[1];

			submod128(work->pt1.X, work->pt1.Z, d1, work->n);
			addmod128(work->pt1.X, work->pt1.Z, s1, work->n);

			// point2 is [2]P
			duplicate(mdata, work, s1, d1, &work->pt1);
		}
		else if (steps[i] == 3)
		{
			// integrate step 4 followed by swap(1,2)
			add(mdata, work, &work->pt2, &work->pt1, &work->pt3, &work->pt4);		// T = B + A (C)
			rot128(work->pt2.X, work->pt4.X, work->pt3.X);
			rot128(work->pt2.Z, work->pt4.Z, work->pt3.Z);
			swap128(work->pt1.X, work->pt2.X);
			swap128(work->pt1.Z, work->pt2.Z);
		}
		else if (steps[i] == 4)
		{
			add(mdata, work, &work->pt2, &work->pt1, &work->pt3, &work->pt4);		// T = B + A (C)
			rot128(work->pt2.X, work->pt4.X, work->pt3.X);
			rot128(work->pt2.Z, work->pt4.Z, work->pt3.Z);
		}
		else if (steps[i] == 5)
		{
			add(mdata, work, &work->pt2, &work->pt1, &work->pt3, &work->pt2);		// B = B + A (C)

			submod128(work->pt1.X, work->pt1.Z, d2, work->n);
			addmod128(work->pt1.X, work->pt1.Z, s2, work->n);
			duplicate(mdata, work, s2, d2, &work->pt1);		// A = 2A
		}
		else if (steps[i] == 6)
		{
			add(mdata, work, &work->pt1, &work->pt2, &work->pt3, P);		// A = A + B (C)
		}
	}

	return;

}

void prac_good(monty128_t *mdata, tinyecm_work *work, tinyecm_pt *P, uint64_t c, double v_in)
{
	uint64_t d, e, r;
	double cmin, cost, v;
	int i;
	uint64_t *s1, *s2, *d1, *d2;
	uint32_t *sw_x, *sw_z;
	int shift = 0;

	while ((c & 1) == 0)
	{
		shift++;
		c >>= 1;
	}

	/* 1/val[0] = the golden ratio (1+sqrt(5))/2, and 1/val[i] for i>0
	   is the real number whose continued fraction expansion is all 1s
	   except for a 2 in i+1-st place */
	if (v_in < 0.0001)
	{
		static double val[NV] =
		{ 0.61803398874989485, 0.72360679774997897, 0.58017872829546410,
		  0.63283980608870629, 0.61242994950949500, 0.62018198080741576,
		  0.61721461653440386, 0.61834711965622806, 0.61791440652881789,
		  0.61807966846989581 };

		/* chooses the best value of v */
		for (d = 0, cmin = ADD * (double)c; d < NV; d++)
		{
			cost = lucas_cost(c, val[d]);
			if (cost < cmin)
			{
				cmin = cost;
				i = d;
			}
		}
		v = val[i];
	}
	else
	{
		v = v_in;
	}

	d = c;
	r = (uint64_t)((double)d * v + 0.5);

	s1 = work->sum1;
	s2 = work->sum2;
	d1 = work->diff1;
	d2 = work->diff2;

	/* first iteration always begins by Condition 3, then a swap */
	d = c - r;
	if ((c >= (2 * r)) || (r >= c))
	{
		printf("problem\n");
		exit(1);
	}
	e = 2 * r - c;

	// mpres_set(xB, xA, n);
	// mpres_set(zB, zA, n); /* B=A */
	// mpres_set(xC, xA, n);
	// mpres_set(zC, zA, n); /* C=A */
	// duplicate(xA, zA, xA, zA, n, b, u, v, w); /* A = 2*A */

	// the first one is always a doubling
	// point1 is [1]P
	work->pt1.X[0] = work->pt2.X[0] = work->pt3.X[0] = P->X[0];
	work->pt1.X[1] = work->pt2.X[1] = work->pt3.X[1] = P->X[1];
	work->pt1.Z[0] = work->pt2.Z[0] = work->pt3.Z[0] = P->Z[0];
	work->pt1.Z[1] = work->pt2.Z[1] = work->pt3.Z[1] = P->Z[1];

	submod128(work->pt1.X, work->pt1.Z, d1, work->n);
	addmod128(work->pt1.X, work->pt1.Z, s1, work->n);

	// point2 is [2]P
	duplicate(mdata, work, s1, d1, &work->pt1);

	while (d != e)
	{
		if (d < e)
		{
			r = d;
			d = e;
			e = r;
			//mpres_swap(xA, xB, n);
			//mpres_swap(zA, zB, n);
			swap128(work->pt1.X, work->pt2.X);
			swap128(work->pt1.Z, work->pt2.Z);
		}
		/* do the first line of Table 4 whose condition qualifies */
		if (d - e <= e / 4 && ((d + e) % 3) == 0)
		{ /* condition 1 */
			d = (2 * d - e) / 3;
			e = (e - d) / 2;

			add(mdata, work, &work->pt1, &work->pt2, &work->pt3, &work->pt4); // T = A + B (C)
			add(mdata, work, &work->pt4, &work->pt1, &work->pt2, &work->pt5); // T2 = T + A (B)
			add(mdata, work, &work->pt2, &work->pt4, &work->pt1, &work->pt2); // B = B + T (A)

			//add3(xT, zT, xA, zA, xB, zB, xC, zC, n, u, v, w); /* T = f(A,B,C) */
			//add3(xT2, zT2, xT, zT, xA, zA, xB, zB, n, u, v, w); /* T2 = f(T,A,B) */
			//add3(xB, zB, xB, zB, xT, zT, xA, zA, n, u, v, w); /* B = f(B,T,A) */
			//mpres_swap(xA, xT2, n);
			//mpres_swap(zA, zT2, n); /* swap A and T2 */
			swap128(work->pt1.X, work->pt5.X);
			swap128(work->pt1.Z, work->pt5.Z);
		}
		else if (d - e <= e / 4 && (d - e) % 6 == 0)
		{ /* condition 2 */
			d = (d - e) / 2;

			add(mdata, work, &work->pt1, &work->pt2, &work->pt3, &work->pt2);		// B = A + B (C)

			submod128(work->pt1.X, work->pt1.Z, d1, work->n);
			addmod128(work->pt1.X, work->pt1.Z, s1, work->n);
			duplicate(mdata, work, s1, d1, &work->pt1);		// A = 2A

			//add3(xB, zB, xA, zA, xB, zB, xC, zC, n, u, v, w); /* B = f(A,B,C) */
			//duplicate(xA, zA, xA, zA, n, b, u, v, w); /* A = 2*A */

		}
		else if ((d + 3) / 4 <= e)
		{ /* condition 3 */
			d -= e;

			add(mdata, work, &work->pt2, &work->pt1, &work->pt3, &work->pt4);		// T = B + A (C)
			//add3(xT, zT, xB, zB, xA, zA, xC, zC, n, u, v, w); /* T = f(B,A,C) */

			/* circular permutation (B,T,C) */
			//tmp = xB;
			//xB = xT;
			//xT = xC;
			//xC = tmp;
			//tmp = zB;
			//zB = zT;
			//zT = zC;
			//zC = tmp;
			rot128(work->pt2.X, work->pt4.X, work->pt3.X);
			rot128(work->pt2.Z, work->pt4.Z, work->pt3.Z);
		}
		else if ((d + e) % 2 == 0)
		{ /* condition 4 */
			d = (d - e) / 2;

			add(mdata, work, &work->pt2, &work->pt1, &work->pt3, &work->pt2);		// B = B + A (C)

			submod128(work->pt1.X, work->pt1.Z, d2, work->n);
			addmod128(work->pt1.X, work->pt1.Z, s2, work->n);
			duplicate(mdata, work, s2, d2, &work->pt1);		// A = 2A

			//add3(xB, zB, xB, zB, xA, zA, xC, zC, n, u, v, w); /* B = f(B,A,C) */
			//duplicate(xA, zA, xA, zA, n, b, u, v, w); /* A = 2*A */
		}
		/* now d+e is odd */
		else if (d % 2 == 0)
		{ /* condition 5 */
			d /= 2;

			add(mdata, work, &work->pt3, &work->pt1, &work->pt2, &work->pt3);		// C = C + A (B)

			submod128(work->pt1.X, work->pt1.Z, d2, work->n);
			addmod128(work->pt1.X, work->pt1.Z, s2, work->n);
			duplicate(mdata, work, s2, d2, &work->pt1);		// A = 2A

			//add3(xC, zC, xC, zC, xA, zA, xB, zB, n, u, v, w); /* C = f(C,A,B) */
			//duplicate(xA, zA, xA, zA, n, b, u, v, w); /* A = 2*A */
		}
		/* now d is odd, e is even */
		else if (d % 3 == 0)
		{ /* condition 6 */
			d = d / 3 - e;

			submod128(work->pt1.X, work->pt1.Z, d1, work->n);
			addmod128(work->pt1.X, work->pt1.Z, s1, work->n);
			duplicate(mdata, work, s1, d1, &work->pt4);		// T = 2A

			add(mdata, work, &work->pt1, &work->pt2, &work->pt3, &work->pt5);		// T2 = A + B (C)
			add(mdata, work, &work->pt4, &work->pt1, &work->pt1, &work->pt1);		// A = T + A (A)
			add(mdata, work, &work->pt4, &work->pt5, &work->pt3, &work->pt4);		// T = T + T2 (C)

			//duplicate(xT, zT, xA, zA, n, b, u, v, w); /* T = 2*A */
			//add3(xT2, zT2, xA, zA, xB, zB, xC, zC, n, u, v, w); /* T2 = f(A,B,C) */
			//add3(xA, zA, xT, zT, xA, zA, xA, zA, n, u, v, w); /* A = f(T,A,A) */
			//add3(xT, zT, xT, zT, xT2, zT2, xC, zC, n, u, v, w); /* T = f(T,T2,C) */

			/* circular permutation (C,B,T) */
			//tmp = xC;
			//xC = xB;
			//xB = xT;
			//xT = tmp;
			//tmp = zC;
			//zC = zB;
			//zB = zT;
			//zT = tmp;
			rot128(work->pt3.X, work->pt2.X, work->pt4.X);
			rot128(work->pt3.Z, work->pt2.Z, work->pt4.Z);
		}
		else if ((d + e) % 3 == 0)
		{ /* condition 7 */
			d = (d - 2 * e) / 3;

			add(mdata, work, &work->pt1, &work->pt2, &work->pt3, &work->pt4);		// T = A + B (C)
			add(mdata, work, &work->pt4, &work->pt1, &work->pt2, &work->pt2);		// B = T + A (B)

			submod128(work->pt1.X, work->pt1.Z, d2, work->n);
			addmod128(work->pt1.X, work->pt1.Z, s2, work->n);
			duplicate(mdata, work, s2, d2, &work->pt4);		// T = 2A
			add(mdata, work, &work->pt1, &work->pt4, &work->pt1, &work->pt1);		// A = A + T (A) = 3A

			//add3(xT, zT, xA, zA, xB, zB, xC, zC, n, u, v, w); /* T = f(A,B,C) */
			//add3(xB, zB, xT, zT, xA, zA, xB, zB, n, u, v, w); /* B = f(T,A,B) */
			//duplicate(xT, zT, xA, zA, n, b, u, v, w);
			//add3(xA, zA, xA, zA, xT, zT, xA, zA, n, u, v, w); /* A = 3*A */
		}
		else if ((d - e) % 3 == 0)
		{ /* condition 8 */
			d = (d - e) / 3;

			add(mdata, work, &work->pt1, &work->pt2, &work->pt3, &work->pt4);		// T = A + B (C)
			add(mdata, work, &work->pt3, &work->pt1, &work->pt2, &work->pt3);		// C = C + A (B)

			//add3(xT, zT, xA, zA, xB, zB, xC, zC, n, u, v, w); /* T = f(A,B,C) */
			//add3(xC, zC, xC, zC, xA, zA, xB, zB, n, u, v, w); /* C = f(A,C,B) */
			//mpres_swap(xB, xT, n);
			//mpres_swap(zB, zT, n); /* swap B and T */
			swap128(work->pt2.X, work->pt4.X);
			swap128(work->pt2.Z, work->pt4.Z);

			submod128(work->pt1.X, work->pt1.Z, d1, work->n);
			addmod128(work->pt1.X, work->pt1.Z, s1, work->n);
			duplicate(mdata, work, s1, d1, &work->pt4);		// T = 2A
			add(mdata, work, &work->pt1, &work->pt4, &work->pt1, &work->pt1);		// A = A + T (A) = 3A

			//duplicate(xT, zT, xA, zA, n, b, u, v, w);
			//add3(xA, zA, xA, zA, xT, zT, xA, zA, n, u, v, w); /* A = 3*A */
		}
		else /* necessarily e is even here */
		{ /* condition 9 */
			e /= 2;

			add(mdata, work, &work->pt3, &work->pt2, &work->pt1, &work->pt3);		// C = C + B (A)

			submod128(work->pt2.X, work->pt2.Z, d2, work->n);
			addmod128(work->pt2.X, work->pt2.Z, s2, work->n);
			duplicate(mdata, work, s2, d2, &work->pt2);		// B = 2B

			//add3(xC, zC, xC, zC, xB, zB, xA, zA, n, u, v, w); /* C = f(C,B,A) */
			//duplicate(xB, zB, xB, zB, n, b, u, v, w); /* B = 2*B */
		}
	}

	add(mdata, work, &work->pt1, &work->pt2, &work->pt3, P);		// A = A + B (C)
	//add3(xA, zA, xA, zA, xB, zB, xC, zC, n, u, v, w);

	for (i = 0; i < shift; i++)
	{
		submod128(P->X, P->Z, d1, work->n);
		addmod128(P->X, P->Z, s1, work->n);
		duplicate(mdata, work, s1, d1, P);		// P = 2P
	}

	if (d != 1)
	{
		printf("problem: d != 1\n");
	}

	return;

}

void prac(monty128_t *mdata, tinyecm_work *work, tinyecm_pt *P, uint64_t c, double v)
{
	uint64_t d, e, r;
	int i;
	uint64_t *s1, *s2, *d1, *d2;
	uint64_t swp;

	d = c;
	r = (uint64_t)((double)d * v + 0.5);

	s1 = work->sum1;
	s2 = work->sum2;
	d1 = work->diff1;
	d2 = work->diff2;

	d = c - r;
	e = 2 * r - c;

	// the first one is always a doubling
	// point1 is [1]P
	work->pt1.X[0] = work->pt2.X[0] = work->pt3.X[0] = P->X[0];
	work->pt1.X[1] = work->pt2.X[1] = work->pt3.X[1] = P->X[1];
	work->pt1.Z[0] = work->pt2.Z[0] = work->pt3.Z[0] = P->Z[0];
	work->pt1.Z[1] = work->pt2.Z[1] = work->pt3.Z[1] = P->Z[1];

	submod128(work->pt1.X, work->pt1.Z, d1, work->n);
	addmod128(work->pt1.X, work->pt1.Z, s1, work->n);

	// point2 is [2]P
	duplicate(mdata, work, s1, d1, &work->pt1);

	while (d != e)
	{
		if (d < e)
		{
			r = d;
			d = e;
			e = r;
			swap128(work->pt1.X, work->pt2.X);
			swap128(work->pt1.Z, work->pt2.Z);
		}

		if ((d + 3) / 4 <= e)
		{
			d -= e;

			add(mdata, work, &work->pt2, &work->pt1, &work->pt3, &work->pt4);		// T = B + A (C)
			rot128(work->pt2.X, work->pt4.X, work->pt3.X);
			rot128(work->pt2.Z, work->pt4.Z, work->pt3.Z);
		}
		else if ((d + e) % 2 == 0)
		{
			d = (d - e) / 2;

			add(mdata, work, &work->pt2, &work->pt1, &work->pt3, &work->pt2);		// B = B + A (C)

			submod128(work->pt1.X, work->pt1.Z, d2, work->n);
			addmod128(work->pt1.X, work->pt1.Z, s2, work->n);
			duplicate(mdata, work, s2, d2, &work->pt1);		// A = 2A
		}
		else
		{
			// empirically, tiny B1 values only need the above prac cases.
			// just in case, fall back on this.
			printf("unhandled case in prac\n");
			exit(1);
		}
	}

	add(mdata, work, &work->pt1, &work->pt2, &work->pt3, P);		// A = A + B (C)

	return;

}

static const double INV_2_POW_32 = 1.0 / (double)((uint64_t)(1) << 32);

static uint32_t tecm_lcg_rand_32B(uint32_t lower, uint32_t upper, uint64_t* ploc_lcg)
{
    *ploc_lcg = 6364136223846793005ULL * (*ploc_lcg) + 1442695040888963407ULL;
    return lower + (uint32_t)(
        (double)(upper - lower) * (double)((*ploc_lcg) >> 32) * INV_2_POW_32);
}

void build_one_curve(tinyecm_pt *P, monty128_t *mdata, 
	tinyecm_work *work, uint32_t sigma, uint64_t * lcg_state, int verbose)
{
	base_t t1[2], t2[2], t3[2], t4[2], t5[2], s[3];
	base_t u[2], v[2], n[2];
	mpz_t gmpt, gmpn;
	mpz_init(gmpt);
	mpz_init(gmpn);

	n[0] = mdata->n[0];
	n[1] = mdata->n[1];

	if (verbose)
		printf("n = %016lx%016lx\n", n[1], n[0]);

	if (sigma == 0)
	{
		work->sigma = tecm_lcg_rand_32B(7, (uint32_t)-1, lcg_state);
	}
	else
	{
		work->sigma = sigma;
	}
	sigma = work->sigma;

	u[0] = sigma;
	u[1] = 0;
	to_monty128(mdata, u);
	
	if (verbose)
		printf("monty(sigma) = %016lx%016lx\n", u[1], u[0]);

	t1[0] = 4;
	t1[1] = 0;
	to_monty128(mdata, t1);

	if (verbose)
		printf("monty(4) = %016lx%016lx\n", t1[1], t1[0]);

	mulmod128(u, t1, v, mdata);		// v = 4*sigma

	if (verbose)
		printf("v = 4*sigma = %016lx%016lx\n", v[1], v[0]);

	sqrmod128(u, u, mdata);
	t1[0] = 5;
	t1[1] = 0;
	to_monty128(mdata, t1);
	submod128(u, t1, u, mdata->n);		// u = sigma^2 - 5

	if (verbose)
		printf("u = sigma^2 - 5 = %016lx%016lx\n", u[1], u[0]);

	sqrmod128(u, t1, mdata);
	mulmod128(t1, u, P->X, mdata);	// x = u^3

	if (verbose)
		printf("x = u^3 = %016lx%016lx\n", P->X[1], P->X[0]);

	sqrmod128(v, t1, mdata);
	mulmod128(t1, v, P->Z, mdata);	// z = v^3

	if (verbose)
		printf("z = v^3 = %016lx%016lx\n", P->Z[1], P->Z[0]);

	//compute parameter A
	submod128(v, u, t1, mdata->n);		// (v - u)
	sqrmod128(t1, t2, mdata);
	mulmod128(t2, t1, t4, mdata);	// (v - u)^3

	if (verbose)
		printf("(v - u)^3 = %016lx%016lx\n", t4[1], t4[0]);

	t1[0] = 3;
	t1[1] = 0;
	to_monty128(mdata, t1);
	mulmod128(t1, u, t2, mdata);		// 3u
	addmod128(t2, v, t3, mdata->n);		// 3u + v

	if (verbose)
		printf("3u + v = %016lx%016lx\n", t3[1], t3[0]);

	mulmod128(t3, t4, t1, mdata);	// a = (v-u)^3 * (3u + v)

	if (verbose)
		printf("a = (v-u)^3 * (3u + v) = %016lx%016lx\n", t1[1], t1[0]);

	t2[0] = 16;
	t2[1] = 0;
	to_monty128(mdata, t2);
	mulmod128(P->X, t2, t3, mdata);	// 16*u^3
	mulmod128(t3, v, t4, mdata);		// 16*u^3*v

	if (verbose)
		printf("16*u^3*v = %016lx%016lx\n", t4[1], t4[0]);

	// u holds the denom, t1 holds the numer
	// accomplish the division by multiplying by the modular inverse
	t2[0] = 1;
	t2[1] = 0;
	mulmod128(t4, t2, t4, mdata);	// take t4 out of monty rep
	mulmod128(t1, t2, t1, mdata);	// take t1 out of monty rep

	u128_to_mpz(t4, gmpt);
	u128_to_mpz(n, gmpn);
	mpz_invert(gmpt, gmpt, gmpn);		// gmpt = t4^-1 mod n
	
	mpz_to_u128(gmpt, t3);
	
	if (verbose)
		printf("1/16*u^3*v = %016lx%016lx\n", t3[1], t3[0]);

	u128_to_mpz(t1, gmpn);
	mpz_mul(gmpt, gmpt, gmpn);			// gmpt = t1 * t4 
	u128_to_mpz(n, gmpn);
	mpz_tdiv_r(gmpt, gmpt, gmpn);		// gmpt = t1 * t4 % n
	mpz_to_u128(gmpt, work->s);
	to_monty128(mdata, work->s);

	u128_to_mpz(t4, gmpn);
	u128_to_mpz(t3, gmpt);
	mpz_mul(gmpt, gmpt, gmpn);
	u128_to_mpz(n, gmpn);
	mpz_tdiv_r(gmpt, gmpt, gmpn);

	//if (mpz_cmp_ui(gmpt, 1) != 0)
	//{
	//	gmp_printf("inversion produced result %Zd from %Zd\n", gmpt, gmpn);
	//	//exit(1);
	//}

	
	// t1 = b = (v - u)^3 * (3*u + v) / 16u^3v
	//mulmod128(t3, t1, work->s, s, mdata);

	mpz_clear(gmpt);
	mpz_clear(gmpn);
	return;
}

void tinyecm(mpz_t n, mpz_t f, uint32_t B1, uint32_t B2, uint32_t curves,
    uint64_t* lcg_state, int verbose)
{
	//attempt to factor n with the elliptic curve method
	//following brent and montgomery's papers, and CP's book
	base_t retval;
	base_t i, j;
	int curve;
	int tid;
	char *wstr;
	int found = 0;
	int result;
	uint64_t num_found;
	tinyecm_work work;
	tinyecm_pt P;
	monty128_t mdata;
	uint64_t n128[2];
	uint32_t sigma;

	mpz_to_u128(n, n128);
	tinyecm_work_init(&work);
	monty128_init(&mdata, n128);
	copy128(n128, work.n);
	work.stg1_max = B1;
	work.stg2_max = B2;

	mpz_set_ull(f, 1);
	for (curve = 0; curve < curves; curve++)
	{
		work.stg1Add = 0;
		work.stg1Doub = 0;
		work.last_pid = 0;

		if (verbose)
			printf("commencing curve %d of %u\n", curve, curves);

		sigma = 0;
		build_one_curve(&P, &mdata, &work, sigma, lcg_state, verbose);

		if (verbose)
		{
			printf("curve parameters:\n\tsigma = %u\n", work.sigma);
			printf("\tn = %016lx%016lx\n", work.n[1], work.n[0]);
			printf("\tx = %016lx%016lx\n", P.X[1], P.X[0]);
			printf("\tz = %016lx%016lx\n", P.Z[1], P.Z[0]);
			printf("\tb = %016lx%016lx\n", work.s[1], work.s[0]);
		}

		ecm_stage1(&mdata, &work, &P);

		result = check_factor(P.Z, mdata.n, mdata.mtmp1);
			
		if (result == 1)
		{
			if (verbose)
				printf("\nfound factor %016lx%016lx in stage 1 with sigma = %u\n",
					mdata.mtmp1[1], mdata.mtmp1[0], sigma);

			u128_to_mpz(mdata.mtmp1, f);
			break;
		}

		if (B2 > B1)
		{
			ecm_stage2(&P, &mdata, &work);
			result = check_factor(work.stg2acc, mdata.n, mdata.mtmp1);

			if (result == 1)
			{
				if (verbose)
					printf("\nfound factor %016lx%016lx in stage 2 with sigma = %u\n",
						mdata.mtmp1[1], mdata.mtmp1[0], sigma);

				u128_to_mpz(mdata.mtmp1, f);
				break;
			}
		}
	}

	return;
}


void ecm_stage1(monty128_t *mdata, tinyecm_work *work, tinyecm_pt *P)
{
	int i;
	uint64_t q;
	uint64_t stg1 = (uint64_t)work->stg1_max;

	// handle the only even case 
	q = 2;
	while (q < stg1 * 4)  // jeff: multiplying by 4 improves perf ~1%
	{
		submod128(P->X, P->Z, work->diff1, work->n);
		addmod128(P->X, P->Z, work->sum1, work->n);
		duplicate(mdata, work, work->sum1, work->diff1, P);
		q *= 2;
	}

	if (stg1 == 27)
	{
		prac(mdata, work, P, 3, 0.61803398874989485);
		prac(mdata, work, P, 3, 0.61803398874989485);
		prac(mdata, work, P, 3, 0.61803398874989485);
		prac(mdata, work, P, 5, 0.618033988749894903);
		prac(mdata, work, P, 5, 0.618033988749894903);
		prac(mdata, work, P, 7, 0.618033988749894903);
		prac(mdata, work, P, 11, 0.580178728295464130);
		prac(mdata, work, P, 13, 0.618033988749894903);
		prac(mdata, work, P, 17, 0.618033988749894903);
		prac(mdata, work, P, 19, 0.618033988749894903);
		prac(mdata, work, P, 23, 0.522786351415446049);
	}
	else if (stg1 == 47)
	{
		// jeff: improved perf slightly by using one more uprac for 3,
		// and removing uprac for 47.
		prac(mdata, work, P, 3, 0.618033988749894903);
		prac(mdata, work, P, 3, 0.618033988749894903);
		prac(mdata, work, P, 3, 0.618033988749894903);
		prac(mdata, work, P, 3, 0.618033988749894903);
		prac(mdata, work, P, 5, 0.618033988749894903);
		prac(mdata, work, P, 5, 0.618033988749894903);
		prac(mdata, work, P, 7, 0.618033988749894903);
		prac(mdata, work, P, 11, 0.580178728295464130);
		prac(mdata, work, P, 13, 0.618033988749894903);
		prac(mdata, work, P, 17, 0.618033988749894903);
		prac(mdata, work, P, 19, 0.618033988749894903);
		prac(mdata, work, P, 23, 0.522786351415446049);
		prac(mdata, work, P, 29, 0.548409048446403258);
		prac(mdata, work, P, 31, 0.618033988749894903);
		prac(mdata, work, P, 37, 0.580178728295464130);
		prac(mdata, work, P, 41, 0.548409048446403258);
		prac(mdata, work, P, 43, 0.618033988749894903);
		//        tecm_uprac(mdata, work, P, 47, 0.548409048446403258);
	}
	else if (stg1 == 59)
	{   // jeff: probably stg1 of 59 would benefit from similar changes
		// as stg1 of 47 above, but I didn't bother. Stg1 of 59 seems to
		// always perform worse than stg1 of 47, so there doesn't seem
		// to be any reason to ever use stg1 of 59.
		prac(mdata, work, P, 3, 0.61803398874989485);
		prac(mdata, work, P, 3, 0.61803398874989485);
		prac(mdata, work, P, 3, 0.61803398874989485);
		prac(mdata, work, P, 5, 0.618033988749894903);
		prac(mdata, work, P, 5, 0.618033988749894903);
		prac(mdata, work, P, 7, 0.618033988749894903);
		prac(mdata, work, P, 7, 0.618033988749894903);
		prac(mdata, work, P, 11, 0.580178728295464130);
		prac(mdata, work, P, 13, 0.618033988749894903);
		prac(mdata, work, P, 17, 0.618033988749894903);
		prac(mdata, work, P, 19, 0.618033988749894903);
		prac(mdata, work, P, 23, 0.522786351415446049);
		prac(mdata, work, P, 29, 0.548409048446403258);
		prac(mdata, work, P, 31, 0.618033988749894903);
		prac(mdata, work, P, 1961, 0.552936068843375);   // 37 * 53
		prac(mdata, work, P, 41, 0.548409048446403258);
		prac(mdata, work, P, 43, 0.618033988749894903);
		prac(mdata, work, P, 47, 0.548409048446403258);
		prac(mdata, work, P, 59, 0.548409048446403258);
	}
	else if (stg1 == 70)
	{
		prac70(mdata, work, P);
		i = 19;
	}
	else // if (stg1 >= 85)
	{
		prac85(mdata, work, P);

		if (stg1 == 85)
		{
			prac(mdata, work, P, 61, 0.522786351415446049);
		}
		else
		{
			prac(mdata, work, P, 5, 0.618033988749894903);
			prac(mdata, work, P, 11, 0.580178728295464130);
			//            uecm_uprac(mdata, work, P, 61, 0.522786351415446049);
			prac(mdata, work, P, 89, 0.618033988749894903);
			prac(mdata, work, P, 97, 0.723606797749978936);
			prac(mdata, work, P, 101, 0.556250337855490828);
			prac(mdata, work, P, 107, 0.580178728295464130);
			prac(mdata, work, P, 109, 0.548409048446403258);
			prac(mdata, work, P, 113, 0.618033988749894903);

			if (stg1 == 125)
			{
				// jeff: moved 61 to here
				prac(mdata, work, P, 61, 0.522786351415446049);
				prac(mdata, work, P, 103, 0.632839806088706269);
			}
			else
			{
				prac(mdata, work, P, 7747, 0.552188778811121); // 61 x 127
				prac(mdata, work, P, 131, 0.618033988749894903);
				prac(mdata, work, P, 14111, 0.632839806088706);  // 103 x 137
				prac(mdata, work, P, 20989, 0.620181980807415);  // 139 x 151
				prac(mdata, work, P, 157, 0.640157392785047019);
				prac(mdata, work, P, 163, 0.551390822543526449);

				if (stg1 == 165)
				{
					prac(mdata, work, P, 149, 0.580178728295464130);
				}
				else
				{
					prac(mdata, work, P, 13, 0.618033988749894903);
					prac(mdata, work, P, 167, 0.580178728295464130);
					prac(mdata, work, P, 173, 0.612429949509495031);
					prac(mdata, work, P, 179, 0.618033988749894903);
					prac(mdata, work, P, 181, 0.551390822543526449);
					prac(mdata, work, P, 191, 0.618033988749894903);
					prac(mdata, work, P, 193, 0.618033988749894903);
					prac(mdata, work, P, 29353, 0.580178728295464);  // 149 x 197
					prac(mdata, work, P, 199, 0.551390822543526449);
				}
			}
		}
	}
	return;
}

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


void ecm_stage2(tinyecm_pt* P, monty128_t* mdata, tinyecm_work* work )
{
	int b;
	int i, j, k;
	tinyecm_pt Pa1;
	tinyecm_pt* Pa = &Pa1;
	tinyecm_pt Pb[18];
	tinyecm_pt Pd1;
	tinyecm_pt* Pd = &Pd1;
	const uint8_t* barray = 0;
	int numb;
	uint32_t stg1_max = work->stg1_max;

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

	uint64_t Pbprod[18][2];

	// [1]Q
	//Pb[1] = *P;
	copy128(P->X, Pb[1].X);
	copy128(P->Z, Pb[1].Z);
	mulmod128(Pb[1].X, Pb[1].Z, Pbprod[1], mdata);

	// [2]Q
	submod128(P->X, P->Z, work->diff1, work->n);
	addmod128(P->X, P->Z, work->sum1, work->n);
	duplicate(mdata, work, work->sum1, work->diff1, &Pb[2]);

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

	tinyecm_pt pt5, pt6;

	// Calculate all Pb: the following is specialized for MICRO_ECM_PARAM_D=60
	// [2]Q + [1]Q([1]Q) = [3]Q
	add(mdata, work, &Pb[1], &Pb[2], &Pb[1], &Pb[3]);        // <-- temporary

	// 2*[3]Q = [6]Q
	submod128(Pb[3].X, Pb[3].Z, work->diff1, work->n);
	addmod128(Pb[3].X, Pb[3].Z, work->sum1, work->n);
	duplicate(mdata, work, work->sum1, work->diff1, &pt6);   // pt6 = [6]Q

	// [3]Q + [2]Q([1]Q) = [5]Q
	add(mdata, work, &Pb[3], &Pb[2], &Pb[1], &pt5);    // <-- pt5 = [5]Q
	//Pb[3] = pt5;
	copy128(pt5.X, Pb[3].X);
	copy128(pt5.Z, Pb[3].Z);

	// [6]Q + [5]Q([1]Q) = [11]Q
	add(mdata, work, &pt6, &pt5, &Pb[1], &Pb[4]);    // <-- [11]Q

	i = 3;
	k = 4;
	j = 5;
	while ((j + 12) < MICRO_ECM_PARAM_D)
	{
		// [j+6]Q + [6]Q([j]Q) = [j+12]Q
		add(mdata, work, &pt6, &Pb[k], &Pb[i], &Pb[map[j + 12]]);
		i = k;
		k = map[j + 12];
		j += 6;
	}

	// [6]Q + [1]Q([5]Q) = [7]Q
	add(mdata, work, &pt6, &Pb[1], &pt5, &Pb[3]);    // <-- [7]Q
	i = 1;
	k = 3;
	j = 1;
	while ((j + 12) < MICRO_ECM_PARAM_D)
	{
		// [j+6]Q + [6]Q([j]Q) = [j+12]Q
		add(mdata, work, &pt6, &Pb[k], &Pb[i], &Pb[map[j + 12]]);
		i = k;
		k = map[j + 12];
		j += 6;
	}

	// Pd = [2w]Q
	// [31]Q + [29]Q([2]Q) = [60]Q
	add(mdata, work, &Pb[9], &Pb[10], &Pb[2], Pd);   // <-- [60]Q

#ifdef MICRO_ECM_VERBOSE_PRINTF
	ptadds++;
#endif

	// temporary - make [4]Q
	tinyecm_pt pt4;
	submod128(Pb[2].X, Pb[2].Z, work->diff1, work->n);
	addmod128(Pb[2].X, Pb[2].Z, work->sum1, work->n);
	duplicate(mdata, work, work->sum1, work->diff1, &pt4);   // pt4 = [4]Q


	// make all of the Pbprod's
	for (i = 3; i < 18; i++)
	{
		 mulmod128(Pb[i].X, Pb[i].Z, Pbprod[i], mdata);
	}

	//initialize info needed for giant step
	tinyecm_pt Pad;

	// Pd = [w]Q
	// [17]Q + [13]Q([4]Q) = [30]Q
	add(mdata, work, &Pb[map[17]], &Pb[map[13]], &pt4, &Pad);    // <-- [30]Q

	// [60]Q + [30]Q([30]Q) = [90]Q
	add(mdata, work, Pd, &Pad, &Pad, Pa);

	tinyecm_pt pt90;   // set pt90 = [90]Q
	tinyecm_pt pt60;   // set pt60 = [60]Q

	copy128(Pa->X, pt90.X);
	copy128(Pd->X, pt60.X);
	copy128(Pa->Z, pt90.Z);
	copy128(Pd->Z, pt60.Z);

	// [90]Q + [30]Q([60]Q) = [120]Q
	add(mdata, work, Pa, &Pad, Pd, Pd);

	// [120]Q + [30]Q([90]Q) = [150]Q
	add(mdata, work, Pd, &Pad, Pa, Pa);

	//initialize accumulator
	uint64_t acc[2];
	copy128(mdata->one, acc);

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

		uint64_t pt90prod[2];
		mulmod128(pt90.X, pt90.Z, pt90prod, mdata);

		for (i = 0; i < numsteps; i++)
		{
			b = steps[i];
			// accumulate the cross product  (zimmerman syntax).
			// page 342 in C&P
			uint64_t tt1[2];
			uint64_t tt2[2];
			submod128(pt90.X, Pb[map[b]].X, tt1, work->n);
			addmod128(pt90.Z, Pb[map[b]].Z, tt2, work->n);
			uint64_t tt3[2];
			mulmod128(tt1, tt2, tt3, mdata);
			addmod128(tt3, Pbprod[map[b]], tt1, work->n);
			submod128(tt1, pt90prod, tt2, work->n);

			uint64_t tmp[2];
			mulmod128(acc, tt2, tmp, mdata);
			if ((tmp[0] == 0) && (tmp[1] == 0))
				break;
			copy128(tmp, acc);
		}
	}
	else if (stg1_max == 70)
	{
		// first process these b's with A=120
		static const int steps[15] = { 49,47,41,37,31,23,19,17,13,11,7,29,43,53,59 };
		// we currently have Pd=120

		uint64_t pdprod[2];
		mulmod128(Pd->X, Pd->Z, pdprod, mdata);

		for (i = 0; i < 15; i++)
		{
			b = steps[i];
			// accumulate the cross product  (zimmerman syntax).
			// page 342 in C&P
			uint64_t tt1[2];
			uint64_t tt2[2];
			submod128(Pd->X, Pb[map[b]].X, tt1, work->n);
			addmod128(Pd->Z, Pb[map[b]].Z, tt2, work->n);
			uint64_t tt3[2];
			mulmod128(tt1, tt2, tt3, mdata);
			addmod128(tt3, Pbprod[map[b]], tt1, work->n);
			submod128(tt1, pdprod, tt2, work->n);

			uint64_t tmp[2];
			mulmod128(acc, tt2, tmp, mdata);
			if ((tmp[0] == 0) && (tmp[1] == 0))
				break;
			copy128(tmp, acc);
		}
	}
	else if (stg1_max == 165)
	{
		// Currently we have Pa=150, Pd=120, Pad=30,  and pt60=60, pt90=90
		// Need Pa = 180, Pd = 120, Pad = 60
		// either of these should be fine
		submod128(pt90.X, pt90.Z, work->diff1, work->n);
		addmod128(pt90.X, pt90.Z, work->sum1, work->n);
		duplicate(mdata, work, work->sum1, work->diff1, Pa);

		//Pad = pt60;
		copy128(pt60.X, Pad.X);
		copy128(pt60.Z, Pad.Z);
		// have pa = 180, pd = 120, pad = 60
	}
	else if (stg1_max == 205)
	{
		// Currently we have Pa=150, Pd=120, Pad=30,  and pt60=60, pt90=90
		// need Pa = 210, Pd = 120, Pad = 90

		// [120]Q + [90]Q([30]Q) = [210]Q
		add(mdata, work, Pd, &pt90, &Pad, Pa);

		//Pad = pt90;
		copy128(pt90.X, Pad.X);
		copy128(pt90.Z, Pad.Z);
	}

	//initialize Paprod
	uint64_t Paprod[2];
	mulmod128(Pa->X, Pa->Z, Paprod, mdata);

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
			tinyecm_pt point;
			copy128(Pa->X, point.X);
			copy128(Pa->Z, point.Z);

			//Pa + Pd
			add(mdata, work, Pa, Pd, &Pad, Pa);

			//Pad holds the previous Pa
			copy128(point.X, Pad.X);
			copy128(point.Z, Pad.Z);

			//and Paprod
			mulmod128(Pa->X, Pa->Z, Paprod, mdata);

			i++;
		}

		//we accumulate XrZd - XdZr = (Xr - Xd) * (Zr + Zd) + XdZd - XrZr
		//in CP notation, Pa -> (Xr,Zr), Pb -> (Xd,Zd)

		b = barray[i];
		// accumulate the cross product  (zimmerman syntax).
		// page 342 in C&P
		uint64_t tt1[2];
		uint64_t tt2[2];
		submod128(Pa->X, Pb[map[b]].X, tt1, work->n);
		addmod128(Pa->Z, Pb[map[b]].Z, tt2, work->n);
		uint64_t tt3[2];
		mulmod128(tt1, tt2, tt3, mdata);
		addmod128(tt3, Pbprod[map[b]], tt1, work->n);
		submod128(tt1, Paprod, tt2, work->n);

		uint64_t tmp[2];
		mulmod128(acc, tt2, tmp, mdata);
		if ((tmp[0] == 0) && (tmp[1] == 0))
			break;
		copy128(tmp, acc);
	}

	copy128(acc, work->stg2acc);
	return;
}

int check_factor(uint64_t * Z, uint64_t * n, uint64_t * f)
{
	int status;
	mpz_t gmp_z;
	mpz_t gmp_n;
	mpz_t gmp_f;

	mpz_init(gmp_z);
	mpz_init(gmp_n);
	mpz_init(gmp_f);

	u128_to_mpz(Z, gmp_z);
	u128_to_mpz(n, gmp_n);

	mpz_gcd(gmp_f, gmp_z, gmp_n);

	status = 0;
	if (mpz_cmp_ui(gmp_f, 1) > 0)
	{
		if (mpz_cmp(gmp_f, gmp_n) == 0)
		{
			// check for a square?
			uint32_t t2 = mpz_get_ull(gmp_f) & 31;
			if (t2 == 0 || t2 == 1 || t2 == 4 ||
				t2 == 9 || t2 == 16 || t2 == 17 || t2 == 25)
			{
				mpz_sqrt(gmp_z, gmp_f);
				mpz_mul(gmp_n, gmp_z, gmp_z);
				if (mpz_cmp(gmp_n, gmp_f) == 0)
				{
					mpz_set(gmp_f, gmp_z);
					status = 1;
				}
				else
				{
					mpz_set_ull(gmp_f, 0);
					status = 0;
				}
			}
			else
			{
				mpz_set_ull(gmp_f, 0);
				status = 0;
			}
		}
		else
		{
			status = 1;
		}
	}
	
	mpz_to_u128(gmp_f, f);

	mpz_clear(gmp_f);
	mpz_clear(gmp_n);
	mpz_clear(gmp_z);

	return status;
}

static int tpm1_ewin100[34] = { // 170 muls
		12, 12, 14, 3, 12, 7, 13, 6, 12, 0, 15, 4, 1, 8, 7, 3, 0, 14, 13, 6,
		5, 6, 15, 13, 0, 11, 0, 14, 13, 3, 8, 8, 12, 0 };
static int tpm1_ewin333[119] = { // 595 muls
	2, 5, 1, 6, 12, 5, 1, 5, 0, 15, 3, 13, 2, 4, 2, 0, 4, 13, 11, 9, 4, 5,
	4, 13, 7, 15, 0, 11, 10, 7, 5, 4, 7, 0, 14, 11, 12, 10, 12, 4, 11, 2,
	5, 2, 10, 10, 7, 3, 14, 11, 8, 0, 15, 2, 2, 3, 10, 11, 6, 9, 2, 8, 15,
	4, 12, 13, 14, 13, 0, 7, 3, 3, 12, 3, 9, 8, 4, 6, 15, 0, 3, 9, 11, 14,
	5, 7, 4, 4, 11, 14, 8, 6, 11, 14, 0, 12, 13, 4, 12, 12, 11, 6, 13, 0,
	10, 8, 13, 6, 13, 15, 2, 5, 14, 13, 11, 11, 15, 0, 0 };

static void tpm1_stage1(monty128_t* mdata, uint64_t* P, uint64_t stg1)
{
	int i;
	uint64_t g[16][2];

	g[0][0] = P[0] = mdata->one[0];
	g[0][1] = P[1] = mdata->one[1];

	for (i = 1; i < 16; i++)
		addmod128(g[i - 1], g[i - 1], g[i], mdata->n);

	switch (stg1)
	{
	case 100:
		for (i = 0; i < 34; i++) {
			sqrmod128(P, P, mdata);
			sqrmod128(P, P, mdata);
			sqrmod128(P, P, mdata);
			sqrmod128(P, P, mdata);
			if (tpm1_ewin100[i] > 0) mulmod128(P, g[tpm1_ewin100[i]], P, mdata);
		}
		break;
	case 333:
		for (i = 0; i < 119; i++) {
			sqrmod128(P, P, mdata);
			sqrmod128(P, P, mdata);
			sqrmod128(P, P, mdata);
			sqrmod128(P, P, mdata);
			if (tpm1_ewin333[i] > 0) mulmod128(P, g[tpm1_ewin333[i]], P, mdata);
		}
		break;
	}
	return;
}

static void tpm1_expRL(monty128_t* mdata, uint64_t* P, uint64_t* in, uint32_t m)
{
	uint64_t s[2];
	P[0] = mdata->one[0];
	P[1] = mdata->one[1];
	s[0] = in[0];
	s[1] = in[1];

	while (m > 0)
	{
		if (m & 1)
			mulmod128(P, s, P, mdata);
		sqrmod128(s, s, mdata);
		m >>= 1;
	}
	return;
}

static const uint32_t tpm1_map[60] = {
	0, 1, 2, 0, 0, 0, 0, 3, 0, 0,
	0, 4, 0, 5, 0, 0, 0, 6, 0, 7,
	0, 0, 0, 8, 0, 0, 0, 0, 0, 9,
	0, 10, 0, 0, 0, 0, 0, 11, 0, 0,
	0, 12, 0, 13, 0, 0, 0, 14, 0, 15,
	0, 0, 0, 16, 0, 0, 0, 0, 0, 17 };

/* baby-steps, giant-steps pairing for b1 = 100, w = 60
q0 = 120 */
static int tpm1_num_b1_100_pairs = 279;
static uint8_t tpm1_b1_100_pairs[279] = {
19, 17, 13, 11, 7, 29, 31, 37, 43, 47, 53, 59, 0,
59, 49, 47, 43, 41, 29, 17, 13, 11, 7, 1, 23, 31, 37, 53, 0,
53, 49, 47, 43, 29, 23, 13, 11, 7, 1, 19, 37, 41, 59, 0,
59, 49, 47, 41, 37, 31, 23, 19, 17, 13, 1, 7, 11, 29, 43, 0,
59, 53, 43, 37, 31, 29, 23, 13, 7, 1, 17, 19, 41, 47, 0,
59, 47, 43, 37, 29, 19, 11, 1, 7, 13, 23, 31, 41, 49, 53, 0,
53, 43, 31, 29, 19, 17, 13, 11, 1, 23, 37, 41, 47, 0,
53, 49, 41, 31, 23, 19, 13, 7, 11, 17, 37, 59, 0,
59, 49, 47, 41, 31, 29, 19, 17, 11, 7, 13, 23, 37, 43, 0,
49, 47, 37, 29, 19, 13, 7, 1, 17, 23, 31, 59, 0,
43, 41, 37, 31, 29, 23, 19, 17, 13, 1, 7, 47, 53, 0,
59, 41, 31, 17, 13, 11, 7, 1, 19, 43, 47, 49, 53, 0,
49, 37, 29, 17, 11, 7, 1, 19, 23, 41, 47, 53, 59, 0,
59, 53, 43, 23, 17, 13, 11, 19, 29, 41, 0,
59, 53, 47, 41, 23, 17, 13, 11, 1, 31, 0,
59, 53, 49, 47, 43, 41, 31, 19, 13, 7, 11, 29, 0,
53, 47, 43, 41, 37, 29, 23, 13, 11, 1, 49, 59, 0,
49, 47, 31, 29, 23, 19, 17, 7, 1, 43, 53, 0,
59, 43, 41, 37, 29, 13, 11, 7, 1, 17, 31, 53, 0,
59, 53, 49, 43, 29, 23, 19, 17, 11, 7, 1, 37, 41, 47, 0,
53, 47, 43 };

/* baby-steps, giant-steps pairing for b1 = 333, w = 60
q0 = 360
*/
static int tpm1_num_b1_333_pairs = 837;
static uint8_t tpm1_b1_333_pairs[837] = {
23, 13, 11, 7, 1, 19, 29, 37, 41, 49, 59, 0,
59, 49, 47, 41, 37, 31, 23, 19, 17, 13, 1, 7, 11, 29, 43, 0,
59, 53, 43, 37, 31, 29, 23, 13, 7, 1, 17, 19, 41, 47, 0,
59, 47, 43, 37, 29, 19, 11, 1, 7, 13, 23, 31, 41, 49, 53, 0,
53, 43, 31, 29, 19, 17, 13, 11, 1, 23, 37, 41, 47, 0,
53, 49, 41, 31, 23, 19, 13, 7, 11, 17, 37, 59, 0,
59, 49, 47, 41, 31, 29, 19, 17, 11, 7, 13, 23, 37, 43, 0,
49, 47, 37, 29, 19, 13, 7, 1, 17, 23, 31, 59, 0,
43, 41, 37, 31, 29, 23, 19, 17, 13, 1, 7, 47, 53, 0,
59, 41, 31, 17, 13, 11, 7, 1, 19, 43, 47, 49, 53, 0,
49, 37, 29, 17, 11, 7, 1, 19, 23, 41, 47, 53, 59, 0,
59, 53, 43, 23, 17, 13, 11, 19, 29, 41, 0,
59, 53, 47, 41, 23, 17, 13, 11, 1, 31, 0,
59, 53, 49, 47, 43, 41, 31, 19, 13, 7, 11, 29, 0,
53, 47, 43, 41, 37, 29, 23, 13, 11, 1, 49, 59, 0,
49, 47, 31, 29, 23, 19, 17, 7, 1, 43, 53, 0,
59, 43, 41, 37, 29, 13, 11, 7, 1, 17, 31, 53, 0,
59, 53, 49, 43, 29, 23, 19, 17, 11, 7, 1, 37, 41, 47, 0,
53, 47, 43, 17, 1, 11, 19, 23, 29, 31, 37, 59, 0,
49, 47, 31, 23, 19, 7, 17, 37, 43, 53, 59, 0,
53, 49, 47, 41, 31, 29, 19, 11, 7, 17, 37, 43, 59, 0,
47, 43, 37, 29, 23, 19, 1, 7, 17, 59, 0,
47, 43, 37, 31, 29, 1, 11, 19, 23, 41, 49, 0,
59, 53, 41, 37, 31, 11, 1, 17, 43, 47, 49, 0,
59, 53, 49, 37, 31, 23, 19, 11, 13, 17, 0,
59, 53, 47, 41, 37, 31, 29, 17, 13, 1, 11, 0,
47, 31, 23, 19, 17, 13, 11, 37, 49, 53, 59, 0,
59, 53, 43, 41, 29, 19, 17, 7, 13, 23, 31, 37, 0,
49, 47, 43, 29, 23, 19, 11, 1, 7, 13, 41, 59, 0,
47, 43, 37, 19, 17, 7, 11, 13, 23, 41, 49, 0,
53, 49, 43, 41, 37, 31, 29, 17, 13, 7, 47, 59, 0,
59, 53, 31, 29, 23, 7, 1, 11, 13, 19, 47, 49, 0,
47, 43, 41, 23, 1, 11, 17, 19, 29, 31, 53, 59, 0,
59, 49, 47, 37, 31, 23, 7, 17, 19, 29, 43, 53, 0,
49, 43, 31, 19, 17, 1, 7, 11, 23, 41, 53, 0,
53, 47, 43, 41, 37, 13, 11, 1, 7, 23, 31, 0,
59, 43, 41, 37, 31, 29, 23, 17, 7, 1, 11, 49, 53, 0,
49, 41, 17, 13, 11, 7, 1, 31, 0,
59, 49, 43, 31, 17, 11, 1, 13, 23, 37, 47, 53, 0,
53, 47, 41, 37, 31, 29, 19, 17, 1, 11, 59, 0,
59, 53, 47, 41, 13, 7, 11, 19, 29, 37, 49, 0,
53, 49, 47, 43, 19, 7, 1, 17, 23, 29, 0,
53, 49, 19, 13, 7, 1, 17, 31, 37, 41, 43, 0,
49, 43, 41, 37, 19, 17, 13, 1, 7, 11, 53, 0,
59, 49, 17, 1, 7, 11, 13, 19, 29, 43, 53, 0,
59, 49, 43, 23, 19, 17, 11, 31, 41, 47, 53, 0,
59, 53, 41, 37, 31, 29, 23, 19, 13, 11, 1, 17, 43, 47, 0,
47, 19, 13, 7, 11, 29, 37, 43, 53, 0,
53, 47, 41, 31, 29, 19, 7, 1, 11, 13, 23, 43, 0,
43, 41, 37, 29, 23, 19, 11, 7, 17, 31, 47, 59, 0,
59, 49, 43, 37, 31, 23, 17, 7, 1, 13, 19, 29, 0,
59, 53, 31, 29, 11, 7, 1, 41, 49, 0,
53, 49, 47, 37, 31, 29, 23, 19, 1, 7, 59, 0,
59, 47, 41, 31, 29, 19, 17, 11, 1, 13, 43, 0,
59, 49, 47, 37, 17, 13, 11, 7, 1, 23, 29, 31, 43, 0,
53, 49, 43, 13, 11, 1, 7, 17, 23, 31, 37, 41, 59, 0,
53, 41, 37, 23, 11, 1, 29, 47, 49, 0,
49, 41, 23, 13, 7, 11, 19, 29, 37, 43, 47, 53, 0,
37, 23, 13, 11, 1, 29, 31, 49, 0,
47, 29, 23, 7, 11, 17, 19, 37, 41, 49, 59, 0,
53, 43, 37, 31, 23, 19, 13, 11, 1, 17, 29, 47, 0,
59, 41, 37, 31, 11, 7, 1, 19, 23, 43, 47, 0,
59, 47, 43, 41, 11, 7, 17, 23, 29, 53, 0,
53, 47, 43, 41, 37, 19, 13, 1, 7, 17, 29, 31, 0,
47, 31, 29, 23, 1, 13, 19, 41, 49, 53, 0,
59, 49, 43, 37, 13, 1, 7, 11, 19, 31, 0,
59, 49, 47, 43, 37, 17, 11, 7, 13, 31 };

static void tpm1_stage2_pair(monty128_t* mdata, uint64_t* P, uint64_t* acc, uint32_t b1)
{
	int w = 60;
	uint64_t d[32][2], six[2], x12[2], xmid[2], x24[2], x25[2];
	uint64_t x36[2], x60[2], x72[2], five[2], pw[2], pgiant[2];
	int i, j;
	uint32_t b2 = 25 * b1;

	// we accumulate f(vw) - f(u) where f(n) = n^2, so that
	// we can pair together primes vw+/-u
	// see: P. L. Montgomery, "Speeding the Pollard and Elliptic Curve Methods
	// of Factorization," Mathematics of Computation, Vol 48, No. 177, 1987
	// 
	// b^(n+h)^2 = b^(n^2 + 2h(n) + h^2) = (b^(n^2))*(b^(2hn))*(b^(h^2))

	// u=1, u^2=1, b^u^2 = P
	copy128(P, d[1]);
	sqrmod128(P, d[2], mdata);
	mulmod128(P, d[2], six, mdata);
	sqrmod128(six, six, mdata);                // P^6
	sqrmod128(six, x12, mdata);                // P^12
	sqrmod128(x12, x24, mdata);                // P^24
	mulmod128(x24, x12, x36, mdata);           // P^36
	mulmod128(x24, x36, x60, mdata);           // P^60
	sqrmod128(x36, x72, mdata);                // P^72
	sqrmod128(d[2], five, mdata);
	mulmod128(five, d[1], five, mdata);        // P^5
	mulmod128(x24, d[1], x25, mdata);          // P^25


	// P^7^2 = P^(1+6)^2 = P^(1^2) * P^(12*1) * P^36
	// P^13^2 = P^(7+6)^2 = P^(7^2) * P^(12*7) * P^36
	// P^19^2 = P^(13+6)^2 = P^(13^2) * P^(12*13) * P^36
	// ...

	// 1, 7, 13, 19, 25, 31, 37, 43, 49
	// unnecessary powers will be mapped to scratch d[0].
	j = 1;
	copy128(x12, xmid);
	while ((j + 6) < 50)
	{
		mulmod128(d[tpm1_map[j]], x36, d[tpm1_map[j + 6]], mdata);
		mulmod128(d[tpm1_map[j + 6]], xmid, d[tpm1_map[j + 6]], mdata);
		mulmod128(xmid, x72, xmid, mdata);
		j += 6;
	}

	// P^11^2 = P^(5+6)^2 = P^(5^2) * P^(12*5) * P^36
	// P^17^2 = P^(11+6)^2 = P^(11^2) * P^(12*11) * P^36
	// P^23^2 = P^(17+6)^2 = P^(17^2) * P^(12*17) * P^36
	// ...

	// 11, 17, 23, 29, 35, 41, 47, 53, 59
	// unnecessary powers will be mapped to scratch d[0].
	mulmod128(x25, x36, d[tpm1_map[11]], mdata);
	mulmod128(d[tpm1_map[11]], x60, d[tpm1_map[11]], mdata);
	mulmod128(x60, x72, xmid, mdata);
	j = 11;
	while ((j + 6) < 60)
	{
		mulmod128(d[tpm1_map[j]], x36, d[tpm1_map[j + 6]], mdata);
		mulmod128(d[tpm1_map[j + 6]], xmid, d[tpm1_map[j + 6]], mdata);
		mulmod128(xmid, x72, xmid, mdata);
		j += 6;
	}

	// P^(2w)^2, assumes w=60
	// P^(120^2) = P^(14440)
	tpm1_expRL(mdata, pw, P, 120 * 120);
	uint64_t x14400[2];
	copy128(pw, x14400);
	uint64_t x240x120[2];
	sqrmod128(pw, x240x120, mdata);
	copy128(x240x120, xmid);
	// P^(240^2) = P^(120+120)^2 = P^(120^2) * P^(240*120) * P^(14400)
	// P^(360^2) = P^(240+120)^2 = P^(240^2) * P^(240*240) * P^(14400)
	// P^(480^2) = P^(360+120)^2 = P^(360^2) * P^(240*360) * P^(14400)
	// ...

	copy128(pw, pgiant);
	i = 2 * w;
	while (i < b1)
	{
		mulmod128(pgiant, x14400, pgiant, mdata);
		mulmod128(pgiant, xmid, pgiant, mdata);
		mulmod128(xmid, x240x120, xmid, mdata);
		i += 2 * w;
	}

	copy128(mdata->one, acc);
	switch (b1)
	{
	case 100:

		for (i = 0; i < tpm1_num_b1_100_pairs; i++)
		{
			if (tpm1_b1_100_pairs[i] == 0)
			{
				mulmod128(pgiant, x14400, pgiant, mdata);
				mulmod128(pgiant, xmid, pgiant, mdata);
				mulmod128(xmid, x240x120, xmid, mdata);
				i++;
			}

			// if we happen to pick up all factors of the input with this
			// prime pair, then the modular reduction will become 0.
			// testing for this is simple and allows us to avoid gcd == n.
			uint64_t tmp[2], sub[2];
			submod128(pgiant, d[tpm1_map[tpm1_b1_100_pairs[i]]], sub, mdata->n);
			mulmod128(acc, sub, tmp, mdata);
			if ((tmp[0] == 0) && (tmp[1] == 0)) break;
			else copy128(tmp, acc);
		}
		break;
	case 333:

		for (i = 0; i < tpm1_num_b1_333_pairs; i++)
		{
			if (tpm1_b1_333_pairs[i] == 0)
			{
				mulmod128(pgiant, x14400, pgiant, mdata);
				mulmod128(pgiant, xmid, pgiant, mdata);
				mulmod128(xmid, x240x120, xmid, mdata);
				i++;
			}

			uint64_t tmp[2], sub[2];
			submod128(pgiant, d[tpm1_map[tpm1_b1_333_pairs[i]]], sub, mdata->n);
			mulmod128(acc, sub, tmp, mdata);
			if ((tmp[0] == 0) && (tmp[1] == 0)) break;
			else copy128(tmp, acc);
		}
		break;

	}

	return;
}

static void tinypm1(mpz_t n, mpz_t f, uint32_t B1, uint32_t B2)
{
	//attempt to factor n with the elliptic curve method
	//following brent and montgomery's papers, and CP's book
	int result;
	uint64_t stg1_res[2], q[2], n128[2];
	uint64_t tmp1[2];
	monty128_t mdata;

	mpz_to_u128(n, n128);
	monty128_init(&mdata, n128);

	tmp1[0] = 1;
	tmp1[1] = 0;
	tpm1_stage1(&mdata, stg1_res, B1);
	mulmod128(tmp1, stg1_res, q, &mdata);
	submod128(q, tmp1, q, n128);
	result = check_factor(q, n128, tmp1);

	if (result == 1)
	{
		u128_to_mpz(tmp1, f);
	}
	else if (B2 > B1)
	{
		uint64_t stg2acc[2];
		tpm1_stage2_pair(&mdata, stg1_res, stg2acc, B1);

		tmp1[0] = 1;
		tmp1[1] = 0;
		mulmod128(tmp1, stg2acc, q, &mdata);
		result = check_factor(q, n128, tmp1);

		if (result == 1)
		{
			u128_to_mpz(tmp1, f);
		}
	}

	return;
}

#ifdef USE_AVX512F

/********************* 104-bit Vector Montgomery arith **********************/
#define VECLEN 8

#define tecm_and64 _mm512_and_epi64
#define tecm_storeu64 _mm512_storeu_epi64
#define tecm_add64 _mm512_add_epi64
#define tecm_sub64 _mm512_sub_epi64
#define tecm_set64 _mm512_set1_epi64
#define tecm_srli64 _mm512_srli_epi64
#define tecm_loadu64 _mm512_loadu_epi64
#define tecm_castpd _mm512_castsi512_pd
#define tecm_castepu _mm512_castpd_si512
#define tecm_DIGIT_SIZE 52
#define tecm_DIGIT_MASK 0x000fffffffffffffULL
#define tecm_MB_WIDTH 8
#define tecm_SIMD_BYTES 64

static __m512d dbias;
static __m512i vbias1;
static __m512i vbias2;
static __m512i lo52mask;

typedef struct
{
	vec_u104_t r;
	vec_u104_t n;
	vec_u104_t np;
	vec_u104_t nhat;
	vec_u104_t rhat;
	vec_u104_t rmask;
	vec_u104_t one;
	vec_u104_t mtmp1;
	vec_u104_t mtmp2;
	vec_u104_t mtmp3;
	vec_u104_t mtmp4;
	vec_u104_t rho;
} monty104_x8_t;

#define GCC_VERSION (__GNUC__ * 10000 \
                     + __GNUC_MINOR__ * 100 \
                     + __GNUC_PATCHLEVEL__)

/* Test for GCC > 3.2.0 */
#if GCC_VERSION < 95000
#define _mm512_loadu_epi64 _mm512_load_epi64
#define _mm512_storeu_epi64 _mm512_store_epi64
#endif

#ifdef IFMA

MICRO_ECM_FORCE_INLINE static __m512i tecm_mul52hi(__m512i b, __m512i c)
{
	return _mm512_madd52hi_epu64(_mm512_set1_epi64(0), c, b);
}

MICRO_ECM_FORCE_INLINE static __m512i tecm_mul52lo(__m512i b, __m512i c)
{
	return _mm512_madd52lo_epu64(_mm512_set1_epi64(0), c, b);
}

MICRO_ECM_FORCE_INLINE static void tecm_mul52lohi(__m512i b, __m512i c, __m512i* l, __m512i* h)
{
	*l = _mm512_madd52lo_epu64(_mm512_set1_epi64(0), c, b);
	*h = _mm512_madd52hi_epu64(_mm512_set1_epi64(0), c, b);
	return;
}

#else



__inline static __m512i tecm_mul52lo(__m512i b, __m512i c)
{
	return _mm512_and_si512(_mm512_mullo_epi64(b, c), _mm512_set1_epi64(0x000fffffffffffffull));
}
__inline static __m512i tecm_mul52hi(__m512i b, __m512i c)
{
	__m512d prod1_ld = _mm512_cvtepu64_pd(b);
	__m512d prod2_ld = _mm512_cvtepu64_pd(c);
	prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	return _mm512_sub_epi64(tecm_castepu(prod1_ld), vbias1);
}
__inline static void tecm_mul52lohi(__m512i b, __m512i c, __m512i* l, __m512i* h)
{
	__m512d prod1_ld = _mm512_cvtepu64_pd(b);
	__m512d prod2_ld = _mm512_cvtepu64_pd(c);
	__m512d prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	*h = _mm512_sub_epi64(tecm_castepu(prod1_hd), vbias1);
	prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
	prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	*l = _mm512_castpd_si512(prod1_ld);
	*l = _mm512_and_si512(*l, lo52mask);
	*h = _mm512_and_si512(*h, lo52mask);
	return;
}

#endif

#define tecm_carryprop(lo, hi, mask) \
	{ __m512i a0 = _mm512_srli_epi64(lo, 52);	\
	hi = _mm512_add_epi64(hi, a0);		\
	lo = _mm512_and_epi64(mask, lo); }

__m512i tecm_mm512_addsetc_epi52(__m512i a, __m512i b, __mmask8* cout)
{
	__m512i t = _mm512_add_epi64(a, b);
	*cout = _mm512_cmpgt_epu64_mask(t, _mm512_set1_epi64(0xfffffffffffffULL));
	t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
	return t;
}
__m512i tecm_mm512_mask_addsetc_epi52(__m512i c, __mmask8 mask, __m512i a, __m512i b, __mmask8* cout)
{
	__m512i t = _mm512_add_epi64(a, b);
	*cout = _mm512_mask_cmpgt_epu64_mask(mask, t, _mm512_set1_epi64(0xfffffffffffffULL));
	t = _mm512_mask_and_epi64(c, mask, t, _mm512_set1_epi64(0xfffffffffffffULL));
	return t;
}
__m512i tecm_mm512_subsetc_epi52(__m512i a, __m512i b, __mmask8* cout)
{
	__m512i t = _mm512_sub_epi64(a, b);
	*cout = _mm512_cmpgt_epu64_mask(b, a);
	t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
	return t;
}
__m512i tecm_mm512_mask_subsetc_epi52(__m512i c, __mmask8 mask, __m512i a, __m512i b, __mmask8* cout)
{
	__m512i t = _mm512_sub_epi64(a, b);
	*cout = _mm512_mask_cmpgt_epu64_mask(mask, b, a);
	t = _mm512_mask_and_epi64(c, mask, t, _mm512_set1_epi64(0xfffffffffffffULL));
	return t;
}
__m512i tecm_mm512_adc_epi52(__m512i a, __mmask8 c, __m512i b, __mmask8* cout)
{
	__m512i t = _mm512_add_epi64(a, b);
	t = _mm512_add_epi64(t, _mm512_maskz_set1_epi64(c, 1));
	*cout = _mm512_cmpgt_epu64_mask(t, _mm512_set1_epi64(0xfffffffffffffULL));
	t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
	return t;
}
__m512i tecm_mm512_mask_adc_epi52(__m512i a, __mmask8 m, __mmask8 c, __m512i b, __mmask8* cout)
{
	__m512i t = _mm512_add_epi64(a, b);
	t = _mm512_mask_add_epi64(a, m, t, _mm512_maskz_set1_epi64(c, 1));
	*cout = _mm512_cmpgt_epu64_mask(t, _mm512_set1_epi64(0xfffffffffffffULL));
	t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
	return t;
}
__m512i tecm_mm512_addcarry_epi52(__m512i a, __mmask8 c, __mmask8* cout)
{
	__m512i t = _mm512_add_epi64(a, _mm512_maskz_set1_epi64(c, 1));
	*cout = _mm512_cmpeq_epu64_mask(a, _mm512_set1_epi64(0xfffffffffffffULL));
	t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
	return t;
}
__m512i tecm_mm512_subborrow_epi52(__m512i a, __mmask8 c, __mmask8* cout)
{
	__m512i t = _mm512_sub_epi64(a, _mm512_maskz_set1_epi64(c, 1));
	*cout = _mm512_cmpeq_epu64_mask(a, _mm512_set1_epi64(0));
	t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
	return t;
}
__m512i tecm_mm512_sbb_epi52(__m512i a, __mmask8 c, __m512i b, __mmask8* cout)
{
	__m512i t = _mm512_sub_epi64(a, b);
	*cout = _mm512_cmpgt_epu64_mask(b, a);
	__m512i t2 = _mm512_sub_epi64(t, _mm512_maskz_set1_epi64(c, 1));
	*cout = _mm512_kor(*cout, _mm512_cmpgt_epu64_mask(t2, t));
	t2 = _mm512_and_epi64(t2, _mm512_set1_epi64(0xfffffffffffffULL));
	return t2;
}
__m512i tecm_mm512_mask_sbb_epi52(__m512i a, __mmask8 m, __mmask8 c, __m512i b, __mmask8* cout)
{
	__m512i t = _mm512_mask_sub_epi64(a, m, a, b);
	*cout = _mm512_mask_cmpgt_epu64_mask(m, b, a);
	__m512i t2 = _mm512_mask_sub_epi64(a, m, t, _mm512_maskz_set1_epi64(c, 1));
	*cout = _mm512_kor(*cout, _mm512_mask_cmpgt_epu64_mask(m, t2, t));
	t2 = _mm512_and_epi64(t2, _mm512_set1_epi64(0xfffffffffffffULL));
	return t2;
}

__inline static void tecm_mulredc104_x8(vec_u104_t* p, vec_u104_t* x, vec_u104_t* y, vec_u104_t* N, vec_u104_t* invN)
{
	// invN is the positive variant = 0 - nhat (the standard negative inverse)
	__m512i T0;
	__m512i T1;
	__m512i T2;
	__m512i T3;
	__m512i t0, t1, t2, t3;
	__m512i mN0;
	__m512i mN1;
	__m512i mN2;
	__m512i mN3;
	__m512i m0;
	__m512i m1;
	__m512i lomask52 = _mm512_set1_epi64(0x000fffffffffffffULL);

	__m512i x0 = _mm512_loadu_si512(x->data[0]);
	__m512i x1 = _mm512_loadu_si512(x->data[1]);
	__m512i y0 = _mm512_loadu_si512(y->data[0]);
	__m512i y1 = _mm512_loadu_si512(y->data[1]);
	__m512i invN0 = _mm512_loadu_si512(invN->data[0]);
	__m512i invN1 = _mm512_loadu_si512(invN->data[1]);
	__m512i N0 = _mm512_loadu_si512(N->data[0]);
	__m512i N1 = _mm512_loadu_si512(N->data[1]);
	__m512i N2;
	__m512i zero = _mm512_setzero_si512();

	// T = x*y
	//tecm_mul52lohi(x0, y0, &T0, &T1);
	//tecm_mul52lohi(x0, y1, &t0, &t1);
	//tecm_mul52lohi(x1, y0, &t2, &t3);
	//tecm_mul52lohi(x1, y1, &T2, &T3);

	__m512d x0_ld = _mm512_cvtepu64_pd(x0);
	__m512d y0_ld = _mm512_cvtepu64_pd(y0);
	__m512d x1_ld = _mm512_cvtepu64_pd(x1);
	__m512d y1_ld = _mm512_cvtepu64_pd(y1);

	__m512d prod1_hd = _mm512_fmadd_round_pd(x0_ld, y0_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	__m512d prod2_hd = _mm512_fmadd_round_pd(x0_ld, y1_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	__m512d prod3_hd = _mm512_fmadd_round_pd(x1_ld, y0_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	__m512d prod4_hd = _mm512_fmadd_round_pd(x1_ld, y1_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

	T1 = _mm512_sub_epi64(tecm_castepu(prod1_hd), vbias1);
	t1 = _mm512_sub_epi64(tecm_castepu(prod2_hd), vbias1);
	t3 = _mm512_sub_epi64(tecm_castepu(prod3_hd), vbias1);
	T3 = _mm512_sub_epi64(tecm_castepu(prod4_hd), vbias1);

	prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
	prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
	prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
	prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

	T0 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x0_ld, y0_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));
	t0 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x0_ld, y1_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));
	t2 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x1_ld, y0_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));
	T2 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x1_ld, y1_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));

	T0 = _mm512_and_si512(T0, lo52mask);
	T1 = _mm512_and_si512(T1, lo52mask);
	t0 = _mm512_and_si512(t0, lo52mask);
	t1 = _mm512_and_si512(t1, lo52mask);
	T2 = _mm512_and_si512(T2, lo52mask);
	T3 = _mm512_and_si512(T3, lo52mask);
	t2 = _mm512_and_si512(t2, lo52mask);
	t3 = _mm512_and_si512(t3, lo52mask);

	T1 = _mm512_add_epi64(T1, t0);
	T2 = _mm512_add_epi64(T2, t1);
	T1 = _mm512_add_epi64(T1, t2);
	T2 = _mm512_add_epi64(T2, t3);
	tecm_carryprop(T1, T2, lomask52);
	tecm_carryprop(T2, T3, lomask52);

	// get the low 104 bits of m = Tlo * invN
	//tecm_mul52lohi(T0, invN0, &m0, &m1);

	x0_ld = _mm512_cvtepu64_pd(T0);
	y0_ld = _mm512_cvtepu64_pd(invN0);
	prod1_hd = _mm512_fmadd_round_pd(x0_ld, y0_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	m1 = _mm512_sub_epi64(tecm_castepu(prod1_hd), vbias1);
	prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
	m0 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x0_ld, y0_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));
	m0 = _mm512_and_si512(m0, lo52mask);
	m1 = _mm512_and_si512(m1, lo52mask);

	m1 = _mm512_add_epi64(m1, tecm_mul52lo(T0, invN1));
	m1 = _mm512_add_epi64(m1, tecm_mul52lo(T1, invN0));
	m1 = _mm512_and_si512(m1, lomask52);

	// get the high 104 bits of m * N
	//tecm_mul52lohi(m0, N0, &mN0, &mN1);
	//tecm_mul52lohi(m0, N1, &t0, &t1);
	//tecm_mul52lohi(m1, N0, &t2, &t3);
	//tecm_mul52lohi(m1, N1, &mN2, &mN3);

	x0_ld = _mm512_cvtepu64_pd(m0);
	y0_ld = _mm512_cvtepu64_pd(N0);
	x1_ld = _mm512_cvtepu64_pd(m1);
	y1_ld = _mm512_cvtepu64_pd(N1);

	prod1_hd = _mm512_fmadd_round_pd(x0_ld, y0_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	prod2_hd = _mm512_fmadd_round_pd(x0_ld, y1_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	prod3_hd = _mm512_fmadd_round_pd(x1_ld, y0_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	prod4_hd = _mm512_fmadd_round_pd(x1_ld, y1_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

	mN1 = _mm512_sub_epi64(tecm_castepu(prod1_hd), vbias1);
	t1 = _mm512_sub_epi64(tecm_castepu(prod2_hd), vbias1);
	t3 = _mm512_sub_epi64(tecm_castepu(prod3_hd), vbias1);
	mN3 = _mm512_sub_epi64(tecm_castepu(prod4_hd), vbias1);

	prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
	prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
	prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
	prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

	mN0 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x0_ld, y0_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));
	t0 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x0_ld, y1_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));
	t2 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x1_ld, y0_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));
	mN2 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x1_ld, y1_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));

	mN0 = _mm512_and_si512(mN0, lo52mask);
	mN1 = _mm512_and_si512(mN1, lo52mask);
	t0 = _mm512_and_si512(t0, lo52mask);
	t1 = _mm512_and_si512(t1, lo52mask);
	mN2 = _mm512_and_si512(mN2, lo52mask);
	mN3 = _mm512_and_si512(mN3, lo52mask);
	t2 = _mm512_and_si512(t2, lo52mask);
	t3 = _mm512_and_si512(t3, lo52mask);

	mN1 = _mm512_add_epi64(mN1, t0);
	mN2 = _mm512_add_epi64(mN2, t1);
	mN1 = _mm512_add_epi64(mN1, t2);
	mN2 = _mm512_add_epi64(mN2, t3);
	tecm_carryprop(mN1, mN2, lomask52);
	tecm_carryprop(mN2, mN3, lomask52);

	// result is (T - mN) >> 104
	// unless T_hi < mN_hi
	// then return (T_hi + N - mN_hi).

	// compare
	__mmask8 msk = _mm512_cmplt_epu64_mask(T3, mN3);
	msk |= _mm512_cmpeq_epu64_mask(T3, mN3) & _mm512_cmplt_epu64_mask(T2, mN2);

	// subtract
	__mmask8 bmsk;
	T2 = tecm_mm512_subsetc_epi52(T2, mN2, &bmsk);
	T3 = tecm_mm512_sbb_epi52(T3, bmsk, mN3, &bmsk);

	// conditionally add N
	T2 = tecm_mm512_mask_addsetc_epi52(T2, msk, T2, N0, &bmsk);
	T3 = tecm_mm512_mask_adc_epi52(T3, msk, bmsk, N1, &bmsk);

	_mm512_storeu_si512(p->data[0], _mm512_and_epi64(T2, lo52mask));
	_mm512_storeu_si512(p->data[1], _mm512_and_epi64(T3, lo52mask));
	return;
}
__inline static void tecm_sqrredc104_x8(vec_u104_t* p, vec_u104_t* x, vec_u104_t* N, vec_u104_t* invN)
{
	tecm_mulredc104_x8(p, x, x, N, invN);
}
__inline static void tecm_addmod104_x8(vec_u104_t* z, vec_u104_t* x, vec_u104_t* y, vec_u104_t* N)
{
	__m512i x0 = _mm512_loadu_si512(x->data[0]);
	__m512i x1 = _mm512_loadu_si512(x->data[1]);
	__m512i y0 = _mm512_loadu_si512(y->data[0]);
	__m512i y1 = _mm512_loadu_si512(y->data[1]);
	__m512i N0 = _mm512_loadu_si512(N->data[0]);
	__m512i N1 = _mm512_loadu_si512(N->data[1]);

	// add
	__mmask8 bmsk;
	x0 = tecm_mm512_addsetc_epi52(x0, y0, &bmsk);
	x1 = tecm_mm512_adc_epi52(x1, bmsk, y1, &bmsk);

	// compare
	__mmask8 msk = bmsk | _mm512_cmpgt_epu64_mask(x1, N1);
	msk |= (_mm512_cmpeq_epu64_mask(x1, N1) & _mm512_cmpge_epu64_mask(x0, N0));

	// conditionally subtract N
	x0 = tecm_mm512_mask_subsetc_epi52(x0, msk, x0, N0, &bmsk);
	x1 = tecm_mm512_mask_sbb_epi52(x1, msk, bmsk, N1, &bmsk);

	_mm512_storeu_si512(z->data[0], _mm512_and_epi64(x0, lo52mask));
	_mm512_storeu_si512(z->data[1], _mm512_and_epi64(x1, lo52mask));
	return;
}
__inline static void tecm_submod104_x8(vec_u104_t* z, vec_u104_t* x, vec_u104_t* y, vec_u104_t* N)
{
	__m512i x0 = _mm512_loadu_si512(x->data[0]);
	__m512i x1 = _mm512_loadu_si512(x->data[1]);
	__m512i y0 = _mm512_loadu_si512(y->data[0]);
	__m512i y1 = _mm512_loadu_si512(y->data[1]);
	__m512i N0 = _mm512_loadu_si512(N->data[0]);
	__m512i N1 = _mm512_loadu_si512(N->data[1]);

	// compare
	__mmask8 msk = _mm512_cmplt_epu64_mask(x1, y1);
	msk |= _mm512_cmpeq_epu64_mask(x1, y1) & _mm512_cmplt_epu64_mask(x0, y0);

	// subtract
	__mmask8 bmsk;
	x0 = tecm_mm512_subsetc_epi52(x0, y0, &bmsk);
	x1 = tecm_mm512_sbb_epi52(x1, bmsk, y1, &bmsk);

	// conditionally add N
	x0 = tecm_mm512_mask_addsetc_epi52(x0, msk, x0, N0, &bmsk);
	x1 = tecm_mm512_mask_adc_epi52(x1, msk, bmsk, N1, &bmsk);

	_mm512_storeu_si512(z->data[0], _mm512_and_epi64(x0, lo52mask));
	_mm512_storeu_si512(z->data[1], _mm512_and_epi64(x1, lo52mask));
	return;
}
__inline static void tecm_submulredc104_x8(vec_u104_t* p, vec_u104_t* w, vec_u104_t* x,
	vec_u104_t* y, vec_u104_t* N, vec_u104_t* invN)
{
	// invN is the positive variant = 0 - nhat (the standard negative inverse)
	__m512i T0;
	__m512i T1;
	__m512i T2;
	__m512i T3;
	__m512i t0, t1, t2, t3;
	__m512i mN0;
	__m512i mN1;
	__m512i mN2;
	__m512i mN3;
	__m512i m0;
	__m512i m1;
	__m512i lomask52 = _mm512_set1_epi64(0x000fffffffffffffULL);

	__m512i x0 = _mm512_loadu_si512(x->data[0]);
	__m512i x1 = _mm512_loadu_si512(x->data[1]);
	__m512i w0 = _mm512_loadu_si512(w->data[0]);
	__m512i w1 = _mm512_loadu_si512(w->data[1]);
	__m512i y0 = _mm512_loadu_si512(y->data[0]);
	__m512i y1 = _mm512_loadu_si512(y->data[1]);
	__m512i invN0 = _mm512_loadu_si512(invN->data[0]);
	__m512i invN1 = _mm512_loadu_si512(invN->data[1]);
	__m512i N0 = _mm512_loadu_si512(N->data[0]);
	__m512i N1 = _mm512_loadu_si512(N->data[1]);
	__m512i N2;
	__m512i zero = _mm512_setzero_si512();

	// compare
	__mmask8 msk = _mm512_cmplt_epu64_mask(w1, x1);
	msk |= _mm512_cmpeq_epu64_mask(w1, x1) & _mm512_cmplt_epu64_mask(w0, x0);

	// subtract
	__mmask8 bmsk;
	x0 = tecm_mm512_subsetc_epi52(w0, x0, &bmsk);
	x1 = tecm_mm512_sbb_epi52(w1, bmsk, x1, &bmsk);

	// conditionally add N
	x0 = tecm_mm512_mask_addsetc_epi52(x0, msk, x0, N0, &bmsk);
	x1 = tecm_mm512_mask_adc_epi52(x1, msk, bmsk, N1, &bmsk);
	x0 = _mm512_and_epi64(x0, lo52mask);
	x1 = _mm512_and_epi64(x1, lo52mask);

	// T = x*y
	//tecm_mul52lohi(x0, y0, &T0, &T1);
	//tecm_mul52lohi(x0, y1, &t0, &t1);
	//tecm_mul52lohi(x1, y0, &t2, &t3);
	//tecm_mul52lohi(x1, y1, &T2, &T3);

	__m512d x0_ld = _mm512_cvtepu64_pd(x0);
	__m512d y0_ld = _mm512_cvtepu64_pd(y0);
	__m512d x1_ld = _mm512_cvtepu64_pd(x1);
	__m512d y1_ld = _mm512_cvtepu64_pd(y1);

	__m512d prod1_hd = _mm512_fmadd_round_pd(x0_ld, y0_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	__m512d prod2_hd = _mm512_fmadd_round_pd(x0_ld, y1_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	__m512d prod3_hd = _mm512_fmadd_round_pd(x1_ld, y0_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	__m512d prod4_hd = _mm512_fmadd_round_pd(x1_ld, y1_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

	T1 = _mm512_sub_epi64(tecm_castepu(prod1_hd), vbias1);
	t1 = _mm512_sub_epi64(tecm_castepu(prod2_hd), vbias1);
	t3 = _mm512_sub_epi64(tecm_castepu(prod3_hd), vbias1);
	T3 = _mm512_sub_epi64(tecm_castepu(prod4_hd), vbias1);

	prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
	prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
	prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
	prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

	T0 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x0_ld, y0_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));
	t0 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x0_ld, y1_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));
	t2 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x1_ld, y0_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));
	T2 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x1_ld, y1_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));

	T0 = _mm512_and_si512(T0, lo52mask);
	T1 = _mm512_and_si512(T1, lo52mask);
	t0 = _mm512_and_si512(t0, lo52mask);
	t1 = _mm512_and_si512(t1, lo52mask);
	T2 = _mm512_and_si512(T2, lo52mask);
	T3 = _mm512_and_si512(T3, lo52mask);
	t2 = _mm512_and_si512(t2, lo52mask);
	t3 = _mm512_and_si512(t3, lo52mask);

	T1 = _mm512_add_epi64(T1, t0);
	T2 = _mm512_add_epi64(T2, t1);
	T1 = _mm512_add_epi64(T1, t2);
	T2 = _mm512_add_epi64(T2, t3);
	tecm_carryprop(T1, T2, lomask52);
	tecm_carryprop(T2, T3, lomask52);

	// get the low 104 bits of m = Tlo * invN
	//tecm_mul52lohi(T0, invN0, &m0, &m1);

	x0_ld = _mm512_cvtepu64_pd(T0);
	y0_ld = _mm512_cvtepu64_pd(invN0);
	prod1_hd = _mm512_fmadd_round_pd(x0_ld, y0_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	m1 = _mm512_sub_epi64(tecm_castepu(prod1_hd), vbias1);
	prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
	m0 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x0_ld, y0_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));
	m0 = _mm512_and_si512(m0, lo52mask);
	m1 = _mm512_and_si512(m1, lo52mask);

	m1 = _mm512_add_epi64(m1, tecm_mul52lo(T0, invN1));
	m1 = _mm512_add_epi64(m1, tecm_mul52lo(T1, invN0));
	m1 = _mm512_and_si512(m1, lomask52);

	// get the high 104 bits of m * N
	//tecm_mul52lohi(m0, N0, &mN0, &mN1);
	//tecm_mul52lohi(m0, N1, &t0, &t1);
	//tecm_mul52lohi(m1, N0, &t2, &t3);
	//tecm_mul52lohi(m1, N1, &mN2, &mN3);

	x0_ld = _mm512_cvtepu64_pd(m0);
	y0_ld = _mm512_cvtepu64_pd(N0);
	x1_ld = _mm512_cvtepu64_pd(m1);
	y1_ld = _mm512_cvtepu64_pd(N1);

	prod1_hd = _mm512_fmadd_round_pd(x0_ld, y0_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	prod2_hd = _mm512_fmadd_round_pd(x0_ld, y1_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	prod3_hd = _mm512_fmadd_round_pd(x1_ld, y0_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	prod4_hd = _mm512_fmadd_round_pd(x1_ld, y1_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

	mN1 = _mm512_sub_epi64(tecm_castepu(prod1_hd), vbias1);
	t1 = _mm512_sub_epi64(tecm_castepu(prod2_hd), vbias1);
	t3 = _mm512_sub_epi64(tecm_castepu(prod3_hd), vbias1);
	mN3 = _mm512_sub_epi64(tecm_castepu(prod4_hd), vbias1);

	prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
	prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
	prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
	prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

	mN0 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x0_ld, y0_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));
	t0 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x0_ld, y1_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));
	t2 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x1_ld, y0_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));
	mN2 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x1_ld, y1_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));

	mN0 = _mm512_and_si512(mN0, lo52mask);
	mN1 = _mm512_and_si512(mN1, lo52mask);
	t0 = _mm512_and_si512(t0, lo52mask);
	t1 = _mm512_and_si512(t1, lo52mask);
	mN2 = _mm512_and_si512(mN2, lo52mask);
	mN3 = _mm512_and_si512(mN3, lo52mask);
	t2 = _mm512_and_si512(t2, lo52mask);
	t3 = _mm512_and_si512(t3, lo52mask);

	mN1 = _mm512_add_epi64(mN1, t0);
	mN2 = _mm512_add_epi64(mN2, t1);
	mN1 = _mm512_add_epi64(mN1, t2);
	mN2 = _mm512_add_epi64(mN2, t3);
	tecm_carryprop(mN1, mN2, lomask52);
	tecm_carryprop(mN2, mN3, lomask52);

	// result is (T - mN) >> 104
	// unless T_hi < mN_hi
	// then return (T_hi + N - mN_hi).

	// compare
	msk = _mm512_cmplt_epu64_mask(T3, mN3);
	msk |= _mm512_cmpeq_epu64_mask(T3, mN3) & _mm512_cmplt_epu64_mask(T2, mN2);

	// subtract
	bmsk;
	T2 = tecm_mm512_subsetc_epi52(T2, mN2, &bmsk);
	T3 = tecm_mm512_sbb_epi52(T3, bmsk, mN3, &bmsk);

	// conditionally add N
	T2 = tecm_mm512_mask_addsetc_epi52(T2, msk, T2, N0, &bmsk);
	T3 = tecm_mm512_mask_adc_epi52(T3, msk, bmsk, N1, &bmsk);

	_mm512_storeu_si512(p->data[0], _mm512_and_epi64(T2, lo52mask));
	_mm512_storeu_si512(p->data[1], _mm512_and_epi64(T3, lo52mask));
	return;
}
__inline static void tecm_addmulredc104_x8(vec_u104_t* p, vec_u104_t* w, vec_u104_t* x,
	vec_u104_t* y, vec_u104_t* N, vec_u104_t* invN)
{
	// invN is the positive variant = 0 - nhat (the standard negative inverse)
	__m512i T0;
	__m512i T1;
	__m512i T2;
	__m512i T3;
	__m512i t0, t1, t2, t3;
	__m512i mN0;
	__m512i mN1;
	__m512i mN2;
	__m512i mN3;
	__m512i m0;
	__m512i m1;
	__m512i lomask52 = _mm512_set1_epi64(0x000fffffffffffffULL);

	__m512i x0 = _mm512_loadu_si512(x->data[0]);
	__m512i x1 = _mm512_loadu_si512(x->data[1]);
	__m512i w0 = _mm512_loadu_si512(w->data[0]);
	__m512i w1 = _mm512_loadu_si512(w->data[1]);
	__m512i y0 = _mm512_loadu_si512(y->data[0]);
	__m512i y1 = _mm512_loadu_si512(y->data[1]);
	__m512i invN0 = _mm512_loadu_si512(invN->data[0]);
	__m512i invN1 = _mm512_loadu_si512(invN->data[1]);
	__m512i N0 = _mm512_loadu_si512(N->data[0]);
	__m512i N1 = _mm512_loadu_si512(N->data[1]);
	__m512i N2;
	__m512i zero = _mm512_setzero_si512();

	// add
	__mmask8 bmsk;
	x0 = tecm_mm512_addsetc_epi52(w0, x0, &bmsk);
	x1 = tecm_mm512_adc_epi52(w1, bmsk, x1, &bmsk);

	// compare
	__mmask8 msk = bmsk | _mm512_cmpgt_epu64_mask(x1, N1);
	msk |= (_mm512_cmpeq_epu64_mask(x1, N1) & _mm512_cmpge_epu64_mask(x0, N0));

	// conditionally subtract N
	x0 = tecm_mm512_mask_subsetc_epi52(x0, msk, x0, N0, &bmsk);
	x1 = tecm_mm512_mask_sbb_epi52(x1, msk, bmsk, N1, &bmsk);
	x0 = _mm512_and_epi64(x0, lo52mask);
	x1 = _mm512_and_epi64(x1, lo52mask);

	// T = x*y
	//tecm_mul52lohi(x0, y0, &T0, &T1);
	//tecm_mul52lohi(x0, y1, &t0, &t1);
	//tecm_mul52lohi(x1, y0, &t2, &t3);
	//tecm_mul52lohi(x1, y1, &T2, &T3);

	__m512d x0_ld = _mm512_cvtepu64_pd(x0);
	__m512d y0_ld = _mm512_cvtepu64_pd(y0);
	__m512d x1_ld = _mm512_cvtepu64_pd(x1);
	__m512d y1_ld = _mm512_cvtepu64_pd(y1);

	__m512d prod1_hd = _mm512_fmadd_round_pd(x0_ld, y0_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	__m512d prod2_hd = _mm512_fmadd_round_pd(x0_ld, y1_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	__m512d prod3_hd = _mm512_fmadd_round_pd(x1_ld, y0_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	__m512d prod4_hd = _mm512_fmadd_round_pd(x1_ld, y1_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

	T1 = _mm512_sub_epi64(tecm_castepu(prod1_hd), vbias1);
	t1 = _mm512_sub_epi64(tecm_castepu(prod2_hd), vbias1);
	t3 = _mm512_sub_epi64(tecm_castepu(prod3_hd), vbias1);
	T3 = _mm512_sub_epi64(tecm_castepu(prod4_hd), vbias1);

	prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
	prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
	prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
	prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

	T0 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x0_ld, y0_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));
	t0 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x0_ld, y1_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));
	t2 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x1_ld, y0_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));
	T2 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x1_ld, y1_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));

	T0 = _mm512_and_si512(T0, lo52mask);
	T1 = _mm512_and_si512(T1, lo52mask);
	t0 = _mm512_and_si512(t0, lo52mask);
	t1 = _mm512_and_si512(t1, lo52mask);
	T2 = _mm512_and_si512(T2, lo52mask);
	T3 = _mm512_and_si512(T3, lo52mask);
	t2 = _mm512_and_si512(t2, lo52mask);
	t3 = _mm512_and_si512(t3, lo52mask);

	T1 = _mm512_add_epi64(T1, t0);
	T2 = _mm512_add_epi64(T2, t1);
	T1 = _mm512_add_epi64(T1, t2);
	T2 = _mm512_add_epi64(T2, t3);
	tecm_carryprop(T1, T2, lomask52);
	tecm_carryprop(T2, T3, lomask52);

	// get the low 104 bits of m = Tlo * invN
	//tecm_mul52lohi(T0, invN0, &m0, &m1);

	x0_ld = _mm512_cvtepu64_pd(T0);
	y0_ld = _mm512_cvtepu64_pd(invN0);
	prod1_hd = _mm512_fmadd_round_pd(x0_ld, y0_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	m1 = _mm512_sub_epi64(tecm_castepu(prod1_hd), vbias1);
	prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
	m0 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x0_ld, y0_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));
	m0 = _mm512_and_si512(m0, lo52mask);
	m1 = _mm512_and_si512(m1, lo52mask);

	m1 = _mm512_add_epi64(m1, tecm_mul52lo(T0, invN1));
	m1 = _mm512_add_epi64(m1, tecm_mul52lo(T1, invN0));
	m1 = _mm512_and_si512(m1, lomask52);

	// get the high 104 bits of m * N
	//tecm_mul52lohi(m0, N0, &mN0, &mN1);
	//tecm_mul52lohi(m0, N1, &t0, &t1);
	//tecm_mul52lohi(m1, N0, &t2, &t3);
	//tecm_mul52lohi(m1, N1, &mN2, &mN3);

	x0_ld = _mm512_cvtepu64_pd(m0);
	y0_ld = _mm512_cvtepu64_pd(N0);
	x1_ld = _mm512_cvtepu64_pd(m1);
	y1_ld = _mm512_cvtepu64_pd(N1);

	prod1_hd = _mm512_fmadd_round_pd(x0_ld, y0_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	prod2_hd = _mm512_fmadd_round_pd(x0_ld, y1_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	prod3_hd = _mm512_fmadd_round_pd(x1_ld, y0_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	prod4_hd = _mm512_fmadd_round_pd(x1_ld, y1_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

	mN1 = _mm512_sub_epi64(tecm_castepu(prod1_hd), vbias1);
	t1 = _mm512_sub_epi64(tecm_castepu(prod2_hd), vbias1);
	t3 = _mm512_sub_epi64(tecm_castepu(prod3_hd), vbias1);
	mN3 = _mm512_sub_epi64(tecm_castepu(prod4_hd), vbias1);

	prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
	prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
	prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
	prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

	mN0 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x0_ld, y0_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));
	t0 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x0_ld, y1_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));
	t2 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x1_ld, y0_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));
	mN2 = _mm512_castpd_si512(
		_mm512_fmadd_round_pd(x1_ld, y1_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));

	mN0 = _mm512_and_si512(mN0, lo52mask);
	mN1 = _mm512_and_si512(mN1, lo52mask);
	t0 = _mm512_and_si512(t0, lo52mask);
	t1 = _mm512_and_si512(t1, lo52mask);
	mN2 = _mm512_and_si512(mN2, lo52mask);
	mN3 = _mm512_and_si512(mN3, lo52mask);
	t2 = _mm512_and_si512(t2, lo52mask);
	t3 = _mm512_and_si512(t3, lo52mask);

	mN1 = _mm512_add_epi64(mN1, t0);
	mN2 = _mm512_add_epi64(mN2, t1);
	mN1 = _mm512_add_epi64(mN1, t2);
	mN2 = _mm512_add_epi64(mN2, t3);
	tecm_carryprop(mN1, mN2, lomask52);
	tecm_carryprop(mN2, mN3, lomask52);

	// result is (T - mN) >> 104
	// unless T_hi < mN_hi
	// then return (T_hi + N - mN_hi).

	// compare
	msk = _mm512_cmplt_epu64_mask(T3, mN3);
	msk |= _mm512_cmpeq_epu64_mask(T3, mN3) & _mm512_cmplt_epu64_mask(T2, mN2);

	// subtract
	bmsk;
	T2 = tecm_mm512_subsetc_epi52(T2, mN2, &bmsk);
	T3 = tecm_mm512_sbb_epi52(T3, bmsk, mN3, &bmsk);

	// conditionally add N
	T2 = tecm_mm512_mask_addsetc_epi52(T2, msk, T2, N0, &bmsk);
	T3 = tecm_mm512_mask_adc_epi52(T3, msk, bmsk, N1, &bmsk);

	_mm512_storeu_si512(p->data[0], _mm512_and_epi64(T2, lo52mask));
	_mm512_storeu_si512(p->data[1], _mm512_and_epi64(T3, lo52mask));
	return;
}

static __inline void copyvec104(vec_u104_t* dest, vec_u104_t* src)
{
	_mm512_storeu_si512(dest->data[0], _mm512_loadu_si512(src->data[0]));
	_mm512_storeu_si512(dest->data[1], _mm512_loadu_si512(src->data[1]));
	return;
}

static __inline void maskcopyvec104(__mmask8 mask, vec_u104_t* dest, vec_u104_t* src)
{
	_mm512_mask_storeu_epi64(dest->data[0], mask, _mm512_loadu_si512(src->data[0]));
	_mm512_mask_storeu_epi64(dest->data[1], mask, _mm512_loadu_si512(src->data[1]));
	return;
}

static __inline __mmask8 iszerovec(vec_u104_t* src)
{
	return _mm512_cmpeq_epu64_mask(
		_mm512_and_epi64(_mm512_loadu_si512(src->data[0]), lo52mask), _mm512_setzero_si512()) &
		_mm512_cmpeq_epu64_mask(
			_mm512_and_epi64(_mm512_loadu_si512(src->data[1]), lo52mask), _mm512_setzero_si512());
}

static void tecm_multiplicative_inverse104(uint64_t* inv, uint64_t a)
{
	//    assert(a%2 == 1);  // the inverse (mod 2<<64) only exists for odd values
	uint64_t x0 = (3 * a) ^ 2;		// 2^4
	uint64_t y = 1 - a * x0;
	uint64_t x1 = x0 * (1 + y);		// 2^8
	y *= y;
	uint64_t x2 = x1 * (1 + y);		// 2^16
	y *= y;
	uint64_t x3 = x2 * (1 + y);		// 2^32
	y *= y;
	uint64_t x4 = x3 * (1 + y);		// 2^64

	// y *= y
	uint64_t yhi, ylo;
	ylo = _umul128(y, y, &yhi);

	// uint64_t x5 = x4 * (1 + y);	// 2^128
	uint64_t x5hi, x5lo;
	if (ylo == 0xffffffffffffffffULL)
	{
		x5lo = 0;
		x5hi = _umul128(yhi + 1, x4, &ylo);
	}
	else
	{
		x5lo = _umul128(1 + ylo, x4, &x5hi);
		x5hi += _umul128(yhi, x4, &ylo);
	}

	inv[0] = x5lo;
	inv[1] = x5hi & 0x000000ffffffffffULL;
	return;
}

__inline static void tecm_uadd_x8(vec_u104_t* rho, vec_u104_t* n,
	tecm_pt_x8 P1, tecm_pt_x8 P2,
	tecm_pt_x8 Pin, tecm_pt_x8* Pout)
{
	// compute:
	//x+ = z- * [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
	//z+ = x- * [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
	// where:
	//x- = original x
	//z- = original z
	// given the sums and differences of the original points
	vec_u104_t diff1, sum1, diff2, sum2, tmp;
	tecm_submod104_x8(&diff1, &P1.X, &P1.Z, n);
	tecm_addmod104_x8(&sum1, &P1.X, &P1.Z, n);
	tecm_submod104_x8(&diff2, &P2.X, &P2.Z, n);
	tecm_addmod104_x8(&sum2, &P2.X, &P2.Z, n);

	vec_u104_t tt1, tt2, tt3, tt4;
	tecm_mulredc104_x8(&tt1, &diff1, &sum2, n, rho); //U
	tecm_mulredc104_x8(&tt2, &sum1, &diff2, n, rho); //V

	tecm_addmod104_x8(&tt3, &tt1, &tt2, n);
	tecm_submod104_x8(&tt4, &tt1, &tt2, n);
	tecm_sqrredc104_x8(&tt1, &tt3, n, rho);   //(U + V)^2
	tecm_sqrredc104_x8(&tt2, &tt4, n, rho);   //(U - V)^2

	tecm_mulredc104_x8(&tmp, &tt1, &Pin.Z, n, rho);     //Z * (U + V)^2
	tecm_mulredc104_x8(&Pout->Z, &tt2, &Pin.X, n, rho);     //x * (U - V)^2
	copyvec104(&Pout->X, &tmp);

	return;
}

__inline static void tecm_udup_x8(vec_u104_t* s, vec_u104_t* rho, vec_u104_t* n,
	vec_u104_t* insum, vec_u104_t* indiff, tecm_pt_x8* P)
{
	vec_u104_t tt1, tt2, tt3;
	tecm_sqrredc104_x8(&tt1, indiff, n, rho);          // U=(x1 - z1)^2
	tecm_sqrredc104_x8(&tt2, insum, n, rho);           // V=(x1 + z1)^2
	tecm_mulredc104_x8(&P->X, &tt1, &tt2, n, rho);         // x=U*V

	tecm_submod104_x8(&tt3, &tt2, &tt1, n);          // w = V-U
	tecm_mulredc104_x8(&tt2, &tt3, s, n, rho);      // w = (A+2)/4 * w
	tecm_addmod104_x8(&tt2, &tt2, &tt1, n);          // w = w + U
	tecm_mulredc104_x8(&P->Z, &tt2, &tt3, n, rho);         // Z = w*(V-U)
	//tecm_addmulredc104_x8(&P->Z, &tt2, &tt1, &tt3, n, rho);         // Z = w*(V-U)
	return;
}

static __inline void copypt104(tecm_pt_x8* dest, tecm_pt_x8* src)
{
	_mm512_storeu_si512(dest->X.data[0], _mm512_loadu_si512(src->X.data[0]));
	_mm512_storeu_si512(dest->X.data[1], _mm512_loadu_si512(src->X.data[1]));
	_mm512_storeu_si512(dest->Z.data[0], _mm512_loadu_si512(src->Z.data[0]));
	_mm512_storeu_si512(dest->Z.data[1], _mm512_loadu_si512(src->Z.data[1]));
	return;
}

static void tecm_prac70_x8(vec_u104_t* rho, vec_u104_t* n, tecm_pt_x8* P, vec_u104_t* s)
{
	vec_u104_t s1, s2, d1, d2;
	tecm_pt_x8 swp;
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

	tecm_pt_x8 pt1, pt2, pt3;
	for (i = 0; i < 116; i++)
	{
		if (steps[i] == 0)
		{
			copypt104(&pt1, P);
			copypt104(&pt2, P);
			copypt104(&pt3, P);

			tecm_submod104_x8(&d1, &pt1.X, &pt1.Z, n);
			tecm_addmod104_x8(&s1, &pt1.X, &pt1.Z, n);
			tecm_udup_x8(s, rho, n, &s1, &d1, &pt1);
		}
		else if (steps[i] == 3)
		{
			// integrate step 4 followed by swap(1,2)
			tecm_pt_x8 pt4;
			tecm_uadd_x8(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

			copypt104(&swp, &pt1);
			copypt104(&pt1, &pt4);
			copypt104(&pt4, &pt3);
			copypt104(&pt3, &pt2);
			copypt104(&pt2, &swp);
		}
		else if (steps[i] == 4)
		{
			tecm_pt_x8 pt4;
			tecm_uadd_x8(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

			copypt104(&swp, &pt2);
			copypt104(&pt2, &pt4);
			copypt104(&pt4, &pt3);
			copypt104(&pt3, &swp);
		}
		else if (steps[i] == 5)
		{
			tecm_submod104_x8(&d2, &pt1.X, &pt1.Z, n);
			tecm_addmod104_x8(&s2, &pt1.X, &pt1.Z, n);

			tecm_uadd_x8(rho, n, pt2, pt1, pt3, &pt2);        // B = B + A (C)
			tecm_udup_x8(s, rho, n, &s2, &d2, &pt1);        // A = 2A
		}
		else if (steps[i] == 6)
		{
			tecm_uadd_x8(rho, n, pt1, pt2, pt3, P);     // A = A + B (C)
		}
	}
	return;
}

static void tecm_prac85_x8(vec_u104_t* rho, vec_u104_t* n, tecm_pt_x8* P, vec_u104_t* s)
{
	vec_u104_t s1, s2, d1, d2;
	tecm_pt_x8 swp;
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

	tecm_pt_x8 pt1, pt2, pt3;
	for (i = 0; i < 146; i++)
	{
		if (steps[i] == 0)
		{
			copypt104(&pt1, P);
			copypt104(&pt2, P);
			copypt104(&pt3, P);

			tecm_submod104_x8(&d1, &pt1.X, &pt1.Z, n);
			tecm_addmod104_x8(&s1, &pt1.X, &pt1.Z, n);
			tecm_udup_x8(s, rho, n, &s1, &d1, &pt1);
		}
		else if (steps[i] == 3)
		{
			// integrate step 4 followed by swap(1,2)
			tecm_pt_x8 pt4;
			tecm_uadd_x8(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

			copypt104(&swp, &pt1);
			copypt104(&pt1, &pt4);
			copypt104(&pt4, &pt3);
			copypt104(&pt3, &pt2);
			copypt104(&pt2, &swp);
		}
		else if (steps[i] == 4)
		{
			tecm_pt_x8 pt4;
			tecm_uadd_x8(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

			copypt104(&swp, &pt2);
			copypt104(&pt2, &pt4);
			copypt104(&pt4, &pt3);
			copypt104(&pt3, &swp);
		}
		else if (steps[i] == 5)
		{
			tecm_submod104_x8(&d2, &pt1.X, &pt1.Z, n);
			tecm_addmod104_x8(&s2, &pt1.X, &pt1.Z, n);

			tecm_uadd_x8(rho, n, pt2, pt1, pt3, &pt2);        // B = B + A (C)
			tecm_udup_x8(s, rho, n, &s2, &d2, &pt1);        // A = 2A
		}
		else if (steps[i] == 6)
		{
			tecm_uadd_x8(rho, n, pt1, pt2, pt3, P);     // A = A + B (C)
		}
	}
	return;
}

static void tecm_prac_x8(vec_u104_t* rho, vec_u104_t* n, tecm_pt_x8* P, uint64_t c, double v, vec_u104_t* s)
{
	uint64_t d, e, r;
	int i;
	vec_u104_t s1, s2, d1, d2;
	tecm_pt_x8 swp;

	// we require c != 0
	int shift = 0;
	
	//= _trail_zcnt64(c);
	//c = c >> shift;
	while ((c & 1) == 0) {
		c >>= 1;
		shift++;
	}

	d = c;
	r = (uint64_t)((double)d * v + 0.5);

	d = c - r;
	e = 2 * r - c;

	tecm_pt_x8 pt1, pt2, pt3;

	// the first one is always a doubling
	// point1 is [1]P
	copypt104(&pt1, P);
	copypt104(&pt2, P);
	copypt104(&pt3, P);

	tecm_submod104_x8(&d1, &pt1.X, &pt1.Z, n);
	tecm_addmod104_x8(&s1, &pt1.X, &pt1.Z, n);

	// point2 is [2]P
	tecm_udup_x8(s, rho, n, &s1, &d1, &pt1);

	while (d != e)
	{
		if (d < e)
		{
			r = d;
			d = e;
			e = r;
			copypt104(&swp, &pt1);
			copypt104(&pt1, &pt2);
			copypt104(&pt2, &swp);
		}
		if (d - e <= e / 4 && ((d + e) % 3) == 0)
		{
			d = (2 * d - e) / 3;
			e = (e - d) / 2;

			tecm_pt_x8 pt4;
			tecm_uadd_x8(rho, n, pt1, pt2, pt3, &pt4); // T = A + B (C)
			tecm_pt_x8 pt5;
			tecm_uadd_x8(rho, n, pt4, pt1, pt2, &pt5); // T2 = T + A (B)
			tecm_uadd_x8(rho, n, pt2, pt4, pt1, &pt2); // B = B + T (A)

			copypt104(&swp, &pt1);
			copypt104(&pt1, &pt5);
			copypt104(&pt5, &swp);
		}
		else if (d - e <= e / 4 && (d - e) % 6 == 0)
		{
			d = (d - e) / 2;

			tecm_submod104_x8(&d1, &pt1.X, &pt1.Z, n);
			tecm_addmod104_x8(&s1, &pt1.X, &pt1.Z, n);

			tecm_uadd_x8(rho, n, pt1, pt2, pt3, &pt2);        // B = A + B (C)
			tecm_udup_x8(s, rho, n, &s1, &d1, &pt1);        // A = 2A
		}
		else if ((d + 3) / 4 <= e)
		{
			d -= e;

			tecm_pt_x8 pt4;
			tecm_uadd_x8(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

			copypt104(&swp, &pt2);
			copypt104(&pt2, &pt4);
			copypt104(&pt4, &pt3);
			copypt104(&pt3, &swp);
		}
		else if ((d + e) % 2 == 0)
		{
			d = (d - e) / 2;

			tecm_submod104_x8(&d2, &pt1.X, &pt1.Z, n);
			tecm_addmod104_x8(&s2, &pt1.X, &pt1.Z, n);

			tecm_uadd_x8(rho, n, pt2, pt1, pt3, &pt2);        // B = B + A (C)
			tecm_udup_x8(s, rho, n, &s2, &d2, &pt1);        // A = 2A
		}
		else if (d % 2 == 0)
		{
			d /= 2;

			tecm_submod104_x8(&d2, &pt1.X, &pt1.Z, n);
			tecm_addmod104_x8(&s2, &pt1.X, &pt1.Z, n);

			tecm_uadd_x8(rho, n, pt3, pt1, pt2, &pt3);        // C = C + A (B)
			tecm_udup_x8(s, rho, n, &s2, &d2, &pt1);        // A = 2A
		}
		else if (d % 3 == 0)
		{
			d = d / 3 - e;

			tecm_submod104_x8(&d1, &pt1.X, &pt1.Z, n);
			tecm_addmod104_x8(&s1, &pt1.X, &pt1.Z, n);

			tecm_pt_x8 pt4;
			tecm_udup_x8(s, rho, n, &s1, &d1, &pt4);        // T = 2A
			tecm_pt_x8 pt5;
			tecm_uadd_x8(rho, n, pt1, pt2, pt3, &pt5);        // T2 = A + B (C)
			tecm_uadd_x8(rho, n, pt4, pt1, pt1, &pt1);        // A = T + A (A)
			tecm_uadd_x8(rho, n, pt4, pt5, pt3, &pt4);        // T = T + T2 (C)

			copypt104(&swp, &pt3);
			copypt104(&pt3, &pt2);
			copypt104(&pt2, &pt4);
			copypt104(&pt4, &swp);
		}
		else if ((d + e) % 3 == 0)
		{
			d = (d - 2 * e) / 3;

			tecm_pt_x8 pt4;
			tecm_uadd_x8(rho, n, pt1, pt2, pt3, &pt4);        // T = A + B (C)


			tecm_submod104_x8(&d2, &pt1.X, &pt1.Z, n);
			tecm_addmod104_x8(&s2, &pt1.X, &pt1.Z, n);
			tecm_uadd_x8(rho, n, pt4, pt1, pt2, &pt2);        // B = T + A (B)
			tecm_udup_x8(s, rho, n, &s2, &d2, &pt4);        // T = 2A
			tecm_uadd_x8(rho, n, pt1, pt4, pt1, &pt1);        // A = A + T (A) = 3A
		}
		else if ((d - e) % 3 == 0)
		{
			d = (d - e) / 3;

			tecm_pt_x8 pt4;
			tecm_uadd_x8(rho, n, pt1, pt2, pt3, &pt4);        // T = A + B (C)

			tecm_submod104_x8(&d2, &pt1.X, &pt1.Z, n);
			tecm_addmod104_x8(&s2, &pt1.X, &pt1.Z, n);
			tecm_uadd_x8(rho, n, pt3, pt1, pt2, &pt3);        // C = C + A (B)

			copypt104(&swp, &pt2);
			copypt104(&pt2, &pt4);
			copypt104(&pt4, &swp);

			tecm_udup_x8(s, rho, n, &s2, &d2, &pt4);        // T = 2A
			tecm_uadd_x8(rho, n, pt1, pt4, pt1, &pt1);        // A = A + T (A) = 3A
		}
		else
		{
			e /= 2;

			tecm_submod104_x8(&d2, &pt2.X, &pt2.Z, n);
			tecm_addmod104_x8(&s2, &pt2.X, &pt2.Z, n);

			tecm_uadd_x8(rho, n, pt3, pt2, pt1, &pt3);        // C = C + B (A)
			tecm_udup_x8(s, rho, n, &s2, &d2, &pt2);        // B = 2B
		}
	}
	tecm_uadd_x8(rho, n, pt1, pt2, pt3, P);     // A = A + B (C)

	for (i = 0; i < shift; i++)
	{
		tecm_submod104_x8(&d1, &P->X, &P->Z, n);
		tecm_addmod104_x8(&s1, &P->X, &P->Z, n);

		tecm_udup_x8(s, rho, n, &s1, &d1, P);     // P = 2P
	}
	return;
}

void ecm_stage1_x8(monty104_x8_t* mdata, tinyecm_work_x8* work, tecm_pt_x8* P)
{
	int i;
	uint64_t q;
	uint64_t stg1 = (uint64_t)work->stg1_max;
	vec_u104_t* rho = &mdata->rho;
	vec_u104_t* n = &mdata->n;
	vec_u104_t* s = &work->s;

	// handle the only even case 
	q = 2;
	while (q < stg1 * 4)  // jeff: multiplying by 4 improves perf ~1%
	{
		tecm_submod104_x8(&work->diff1, &P->X, &P->Z, n);
		tecm_addmod104_x8(&work->sum1, &P->X, &P->Z, n);
		tecm_udup_x8(s, rho, n, &work->sum1, &work->diff1, P);
		q *= 2;
	}

	if (stg1 == 27)
	{
		tecm_prac_x8(rho, n, P, 3, 0.61803398874989485, s);
		tecm_prac_x8(rho, n, P, 3, 0.61803398874989485, s);
		tecm_prac_x8(rho, n, P, 3, 0.61803398874989485, s);
		tecm_prac_x8(rho, n, P, 5, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 5, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 7, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 11, 0.580178728295464130, s);
		tecm_prac_x8(rho, n, P, 13, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 17, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 19, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 23, 0.522786351415446049, s);
	}
	else if (stg1 == 47)
	{
		// jeff: improved perf slightly by using one more uprac for 3,
		// and removing uprac for 47.
		tecm_prac_x8(rho, n, P, 3, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 3, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 3, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 3, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 5, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 5, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 7, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 11, 0.580178728295464130, s);
		tecm_prac_x8(rho, n, P, 13, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 17, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 19, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 23, 0.522786351415446049, s);
		tecm_prac_x8(rho, n, P, 29, 0.548409048446403258, s);
		tecm_prac_x8(rho, n, P, 31, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 37, 0.580178728295464130, s);
		tecm_prac_x8(rho, n, P, 41, 0.548409048446403258, s);
		tecm_prac_x8(rho, n, P, 43, 0.618033988749894903, s);
		//        tecm_uprac(mdata, work, P, 47, 0.548409048446403258);
	}
	else if (stg1 == 59)
	{   // jeff: probably stg1 of 59 would benefit from similar changes
		// as stg1 of 47 above, but I didn't bother. Stg1 of 59 seems to
		// always perform worse than stg1 of 47, so there doesn't seem
		// to be any reason to ever use stg1 of 59.
		tecm_prac_x8(rho, n, P, 3, 0.61803398874989485, s);
		tecm_prac_x8(rho, n, P, 3, 0.61803398874989485, s);
		tecm_prac_x8(rho, n, P, 3, 0.61803398874989485, s);
		tecm_prac_x8(rho, n, P, 5, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 5, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 7, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 7, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 11, 0.580178728295464130, s);
		tecm_prac_x8(rho, n, P, 13, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 17, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 19, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 23, 0.522786351415446049, s);
		tecm_prac_x8(rho, n, P, 29, 0.548409048446403258, s);
		tecm_prac_x8(rho, n, P, 31, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 1961, 0.552936068843375, s);   // 37 * 53
		tecm_prac_x8(rho, n, P, 41, 0.548409048446403258, s);
		tecm_prac_x8(rho, n, P, 43, 0.618033988749894903, s);
		tecm_prac_x8(rho, n, P, 47, 0.548409048446403258, s);
		tecm_prac_x8(rho, n, P, 59, 0.548409048446403258, s);
	}
	else if (stg1 == 70)
	{
		tecm_prac70_x8(rho, n, P, s);
		i = 19;
	}
	else // if (stg1 >= 85)
	{
		tecm_prac85_x8(rho, n, P, s);

		if (stg1 == 85)
		{
			tecm_prac_x8(rho, n, P, 61, 0.522786351415446049, s);
		}
		else
		{
			tecm_prac_x8(rho, n, P, 5, 0.618033988749894903, s);
			tecm_prac_x8(rho, n, P, 11, 0.580178728295464130, s);
			tecm_prac_x8(rho, n, P, 89, 0.618033988749894903, s);
			tecm_prac_x8(rho, n, P, 97, 0.723606797749978936, s);
			tecm_prac_x8(rho, n, P, 101, 0.556250337855490828, s);
			tecm_prac_x8(rho, n, P, 107, 0.580178728295464130, s);
			tecm_prac_x8(rho, n, P, 109, 0.548409048446403258, s);
			tecm_prac_x8(rho, n, P, 113, 0.618033988749894903, s);

			if (stg1 == 125)
			{
				// jeff: moved 61 to here
				tecm_prac_x8(rho, n, P, 61, 0.522786351415446049, s);
				tecm_prac_x8(rho, n, P, 103, 0.632839806088706269, s);
			}
			else
			{
				tecm_prac_x8(rho, n, P, 7747, 0.552188778811121, s); // 61 x 127
				tecm_prac_x8(rho, n, P, 131, 0.618033988749894903, s);
				tecm_prac_x8(rho, n, P, 14111, 0.632839806088706, s);  // 103 x 137
				tecm_prac_x8(rho, n, P, 20989, 0.620181980807415, s);  // 139 x 151
				tecm_prac_x8(rho, n, P, 157, 0.640157392785047019, s);
				tecm_prac_x8(rho, n, P, 163, 0.551390822543526449, s);

				if (stg1 == 165)
				{
					tecm_prac_x8(rho, n, P, 149, 0.580178728295464130, s);
				}
				else
				{
					tecm_prac_x8(rho, n, P, 13, 0.618033988749894903, s);
					tecm_prac_x8(rho, n, P, 167, 0.580178728295464130, s);
					tecm_prac_x8(rho, n, P, 173, 0.612429949509495031, s);
					tecm_prac_x8(rho, n, P, 179, 0.618033988749894903, s);
					tecm_prac_x8(rho, n, P, 181, 0.551390822543526449, s);
					tecm_prac_x8(rho, n, P, 191, 0.618033988749894903, s);
					tecm_prac_x8(rho, n, P, 193, 0.618033988749894903, s);
					tecm_prac_x8(rho, n, P, 29353, 0.580178728295464, s);  // 149 x 197
					tecm_prac_x8(rho, n, P, 199, 0.551390822543526449, s);
				}
			}
		}
	}
	return;
}

void ecm_stage2_x8(tecm_pt_x8* P, monty104_x8_t* mdata, tinyecm_work_x8* work)
{
	int b;
	int i, j, k;
	tecm_pt_x8 Pa1;
	tecm_pt_x8* Pa = &Pa1;
	tecm_pt_x8 Pb[18];
	tecm_pt_x8 Pd1;
	tecm_pt_x8* Pd = &Pd1;
	const uint8_t* barray = 0;
	int numb;
	uint32_t stg1_max = work->stg1_max;
	vec_u104_t* N = &mdata->n;
	vec_u104_t* rho = &mdata->rho;
	vec_u104_t* sum = &work->sum1;
	vec_u104_t* diff = &work->diff1;

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

	vec_u104_t Pbprod[18];

	// [1]Q
	//Pb[1] = *P;
	copypt104(&Pb[1], P);
	tecm_mulredc104_x8(&Pbprod[1], &Pb[1].X, &Pb[1].Z, N, rho);

	// [2]Q
	tecm_submod104_x8(diff, &P->X, &P->Z, N);
	tecm_addmod104_x8(sum, &P->X, &P->Z, N);
	tecm_udup_x8(&work->s, rho, N, sum, diff, &Pb[2]);

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

	tecm_pt_x8 pt5, pt6;

	// Calculate all Pb: the following is specialized for MICRO_ECM_PARAM_D=60
	// [2]Q + [1]Q([1]Q) = [3]Q
	tecm_uadd_x8(rho, N, Pb[1], Pb[2], Pb[1], &Pb[3]);        // <-- temporary

	// 2*[3]Q = [6]Q
	tecm_submod104_x8(diff, &Pb[3].X, &Pb[3].Z, N);
	tecm_addmod104_x8(sum, &Pb[3].X, &Pb[3].Z, N);
	tecm_udup_x8(&work->s, rho, N, sum, diff, &pt6);   // pt6 = [6]Q

	// [3]Q + [2]Q([1]Q) = [5]Q
	tecm_uadd_x8(rho, N, Pb[3], Pb[2], Pb[1], &pt5);    // <-- pt5 = [5]Q
	//Pb[3] = pt5;
	copypt104(&Pb[3], &pt5);

	// [6]Q + [5]Q([1]Q) = [11]Q
	tecm_uadd_x8(rho, N, pt6, pt5, Pb[1], &Pb[4]);    // <-- [11]Q

	i = 3;
	k = 4;
	j = 5;
	while ((j + 12) < MICRO_ECM_PARAM_D)
	{
		// [j+6]Q + [6]Q([j]Q) = [j+12]Q
		tecm_uadd_x8(rho, N, pt6, Pb[k], Pb[i], &Pb[map[j + 12]]);
		i = k;
		k = map[j + 12];
		j += 6;
	}

	// [6]Q + [1]Q([5]Q) = [7]Q
	tecm_uadd_x8(rho, N, pt6, Pb[1], pt5, &Pb[3]);    // <-- [7]Q
	i = 1;
	k = 3;
	j = 1;
	while ((j + 12) < MICRO_ECM_PARAM_D)
	{
		// [j+6]Q + [6]Q([j]Q) = [j+12]Q
		tecm_uadd_x8(rho, N, pt6, Pb[k], Pb[i], &Pb[map[j + 12]]);
		i = k;
		k = map[j + 12];
		j += 6;
	}

	// Pd = [2w]Q
	// [31]Q + [29]Q([2]Q) = [60]Q
	tecm_uadd_x8(rho, N, Pb[9], Pb[10], Pb[2], Pd);   // <-- [60]Q

#ifdef MICRO_ECM_VERBOSE_PRINTF
	ptadds++;
#endif

	// temporary - make [4]Q
	tecm_pt_x8 pt4;
	tecm_submod104_x8(diff, &Pb[2].X, &Pb[2].Z, N);
	tecm_addmod104_x8(sum, &Pb[2].X, &Pb[2].Z, N);
	tecm_udup_x8(&work->s, rho, N, sum, diff, &pt4);   // pt4 = [4]Q


	// make all of the Pbprod's
	for (i = 3; i < 18; i++)
	{
		tecm_mulredc104_x8(&Pbprod[i], &Pb[i].X, &Pb[i].Z, N, rho);
	}

	//initialize info needed for giant step
	tecm_pt_x8 Pad;

	// Pd = [w]Q
	// [17]Q + [13]Q([4]Q) = [30]Q
	tecm_uadd_x8(rho, N, Pb[map[17]], Pb[map[13]], pt4, &Pad);    // <-- [30]Q

	// [60]Q + [30]Q([30]Q) = [90]Q
	tecm_uadd_x8(rho, N, *Pd, Pad, Pad, Pa);

	tecm_pt_x8 pt90;   // set pt90 = [90]Q
	tecm_pt_x8 pt60;   // set pt60 = [60]Q

	copypt104(&pt90, Pa);
	copypt104(&pt60, Pd);

	// [90]Q + [30]Q([60]Q) = [120]Q
	tecm_uadd_x8(rho, N, *Pa, Pad, *Pd, Pd);

	// [120]Q + [30]Q([90]Q) = [150]Q
	tecm_uadd_x8(rho, N, *Pd, Pad, *Pa, Pa);

	//initialize accumulator
	vec_u104_t acc;
	copyvec104(&acc, &mdata->one);

	// adjustment of Pa and Pad for particular B1.
	// Currently we have Pa=150, Pd=120, Pad=30
	__mmask8 msk = 0xff;

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

		vec_u104_t pt90prod;
		tecm_mulredc104_x8(&pt90prod, &pt90.X, &pt90.Z, N, rho);

		for (i = 0; i < numsteps; i++)
		{
			b = steps[i];
			// accumulate the cross product  (zimmerman syntax).
			// page 342 in C&P
			vec_u104_t tt1;
			vec_u104_t tt2;
			tecm_submod104_x8(&tt1, &pt90.X, &Pb[map[b]].X, N);
			tecm_addmod104_x8(&tt2, &pt90.Z, &Pb[map[b]].Z, N);
			vec_u104_t tt3;
			tecm_mulredc104_x8(&tt3, &tt1, &tt2, N, rho);
			tecm_addmod104_x8(&tt1, &tt3, &Pbprod[map[b]], N);
			tecm_submod104_x8(&tt2, &tt1, &pt90prod, N);

			vec_u104_t tmp;
			tecm_mulredc104_x8(&tmp, &acc, &tt2, N, rho);
			__mmask8 iszero = iszerovec(&tmp);
			msk &= (~iszero);
			maskcopyvec104(msk, &acc, &tmp);
		}
	}
	else if (stg1_max == 70)
	{
		// first process these b's with A=120
		static const int steps[15] = { 49,47,41,37,31,23,19,17,13,11,7,29,43,53,59 };
		// we currently have Pd=120

		vec_u104_t pdprod;
		tecm_mulredc104_x8(&pdprod, &Pd->X, &Pd->Z, N, rho);

		for (i = 0; i < 15; i++)
		{
			b = steps[i];
			// accumulate the cross product  (zimmerman syntax).
			// page 342 in C&P
			vec_u104_t tt1;
			vec_u104_t tt2;
			tecm_submod104_x8(&tt1, &Pd->X, &Pb[map[b]].X, N);
			tecm_addmod104_x8(&tt2, &Pd->Z, &Pb[map[b]].Z, N);
			vec_u104_t tt3;
			tecm_mulredc104_x8(&tt3, &tt1, &tt2, N, rho);
			tecm_addmod104_x8(&tt1, &tt3, &Pbprod[map[b]], N);
			tecm_submod104_x8(&tt2, &tt1, &pdprod, N);

			vec_u104_t tmp;
			tecm_mulredc104_x8(&tmp, &acc, &tt2, N, rho);
			__mmask8 iszero = iszerovec(&tmp);
			msk &= (~iszero);
			maskcopyvec104(msk, &acc, &tmp);
		}
	}
	else if (stg1_max == 165)
	{
		// Currently we have Pa=150, Pd=120, Pad=30,  and pt60=60, pt90=90
		// Need Pa = 180, Pd = 120, Pad = 60
		// either of these should be fine
		tecm_submod104_x8(diff, &pt90.X, &pt90.Z, N);
		tecm_addmod104_x8(sum, &pt90.X, &pt90.Z, N);
		tecm_udup_x8(&work->s, rho, N, sum, diff, Pa);

		//Pad = pt60;
		copypt104(&Pad, &pt60);
		// have pa = 180, pd = 120, pad = 60
	}
	else if (stg1_max == 205)
	{
		// Currently we have Pa=150, Pd=120, Pad=30,  and pt60=60, pt90=90
		// need Pa = 210, Pd = 120, Pad = 90

		// [120]Q + [90]Q([30]Q) = [210]Q
		tecm_uadd_x8(rho, N, *Pd, pt90, Pad, Pa);

		//Pad = pt90;
		copypt104(&Pad, &pt90);
	}

	//initialize Paprod
	vec_u104_t Paprod;
	tecm_mulredc104_x8(&Paprod, &Pa->X, &Pa->Z, N, rho);

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
			tecm_pt_x8 point;
			copypt104(&point, Pa);

			//Pa + Pd
			tecm_uadd_x8(rho, N, *Pa, *Pd, Pad, Pa);

			//Pad holds the previous Pa
			copypt104(&Pad, &point);

			//and Paprod
			tecm_mulredc104_x8(&Paprod, &Pa->X, &Pa->Z, N, rho);

			i++;
		}

		//we accumulate XrZd - XdZr = (Xr - Xd) * (Zr + Zd) + XdZd - XrZr
		//in CP notation, Pa -> (Xr,Zr), Pb -> (Xd,Zd)

		b = barray[i];
		// accumulate the cross product  (zimmerman syntax).
		// page 342 in C&P
		vec_u104_t tt1;
		vec_u104_t tt2;
		tecm_submod104_x8(&tt1, &Pa->X, &Pb[map[b]].X, N);
		tecm_addmod104_x8(&tt2, &Pa->Z, &Pb[map[b]].Z, N);
		vec_u104_t tt3;
		tecm_mulredc104_x8(&tt3, &tt1, &tt2, N, rho);
		tecm_addmod104_x8(&tt1, &tt3, &Pbprod[map[b]], N);
		tecm_submod104_x8(&tt2, &tt1, &Paprod, N);

		vec_u104_t tmp;
		tecm_mulredc104_x8(&tmp, &acc, &tt2, N, rho);
		__mmask8 iszero = iszerovec(&tmp);
		msk &= (~iszero);
		maskcopyvec104(msk, &acc, &tmp);
	}

	copyvec104(&work->stg2acc, &acc);
	return;
}

void to_monty104_x8(vec_u104_t* x, vec_u104_t* N)
{
	mpz_t m;
	mpz_t n;

	mpz_init(m);
	mpz_init(n);

	int i;
	for (i = 0; i < 8; i++)
	{
		mpz_set_ull(n, N->data[1][i]);
		mpz_mul_2exp(n, n, 52);
		mpz_add_ull(n, n, N->data[0][i]);

		mpz_set_ull(m, x->data[1][i]);
		mpz_mul_2exp(m, m, 52);
		mpz_add_ull(m, m, x->data[0][i]);

		mpz_mul_2exp(m, m, 104);
		mpz_mod(m, m, n);
		x->data[0][i] = mpz_get_ull(m) & 0x000fffffffffffff;
		mpz_tdiv_q_2exp(m, m, 52);
		x->data[1][i] = mpz_get_ull(m) & 0x000fffffffffffff;
	}

	mpz_clear(m);
	mpz_clear(n);

	return;
}

void u104_to_mpz(vec_u104_t* x, int lane, mpz_t gx)
{
	mpz_set_ull(gx, x->data[1][lane]);
	mpz_mul_2exp(gx, gx, 52);
	mpz_add_ull(gx, gx, x->data[0][lane]);
	return;
}

void mpz_to_u104(mpz_t gx, vec_u104_t* x, int lane)
{
	x->data[0][lane] = mpz_get_ull(gx) & 0x000fffffffffffff;
	mpz_tdiv_q_2exp(gx, gx, 52);
	x->data[1][lane] = mpz_get_ull(gx) & 0x000fffffffffffff;
	// restore input
	mpz_mul_2exp(gx, gx, 52);
	mpz_add_ull(gx, gx, x->data[0][lane]);
	return;
}

void build_curves_104_x8(tecm_pt_x8* P, monty104_x8_t* mdata,
	tinyecm_work_x8* work, uint64_t* lcg_state, int verbose)
{
	vec_u104_t t1, t2, t3, t4, t5, s, * rho = &mdata->rho;
	vec_u104_t u, v, * n = &mdata->n, sigma;
	mpz_t gmpt, gmpn, gmpg, gmpt2;
	mpz_init(gmpt);
	mpz_init(gmpn);
	mpz_init(gmpg);
	mpz_init(gmpt2);
	int i;

	for (i = 0; i < 8; i++)
	{
		work->sigma[i] = sigma.data[0][i] = tecm_lcg_rand_32B(7, (uint32_t)-1, lcg_state);
		sigma.data[1][i] = 0;
		t1.data[0][i] = 4;
		t1.data[1][i] = 0;
		t2.data[0][i] = 5;
		t2.data[1][i] = 0;
		t3.data[0][i] = 3;
		t3.data[1][i] = 0;
		t5.data[0][i] = 16;
		t5.data[1][i] = 0;
	}

	to_monty104_x8(&sigma, n);
	to_monty104_x8(&t1, n);
	to_monty104_x8(&t2, n);
	to_monty104_x8(&t3, n);
	to_monty104_x8(&t5, n);

	verbose = 0;

	tecm_mulredc104_x8(&v, &sigma, &t1, n, rho);	// v = 4*sigma
	if (verbose) printf("v[0] = %016lx%016lx\n", v.data[1][0], v.data[0][0]);
	tecm_sqrredc104_x8(&u, &sigma, n, rho);
	tecm_submod104_x8(&u, &u, &t2, n);				// u = sigma^2 - 5

	tecm_sqrredc104_x8(&t1, &u, n, rho);
	tecm_mulredc104_x8(&P->X, &t1, &u, n, rho);		// x = u^3
	if (verbose) printf("x[0] = %016lx%016lx\n", P->X.data[1][0], P->X.data[0][0]);

	tecm_sqrredc104_x8(&t1, &v, n, rho);
	tecm_mulredc104_x8(&P->Z, &t1, &v, n, rho);		// z = v^3
	if (verbose) printf("z[0] = %016lx%016lx\n", P->Z.data[1][0], P->Z.data[0][0]);

	//compute parameter A
	tecm_submod104_x8(&t1, &v, &u, n);				// (v - u)
	if (verbose) printf("v-u = %016lx%016lx\n", t1.data[1][0], t1.data[0][0]);
	tecm_sqrredc104_x8(&t2, &t1, n, rho);
	tecm_mulredc104_x8(&t4, &t1, &t2, n, rho);		// (v - u)^3
	if (verbose) printf("(v-u)^3 = %016lx%016lx\n", t4.data[1][0], t4.data[0][0]);
	tecm_mulredc104_x8(&t2, &u, &t3, n, rho);		// 3u
	if (verbose) printf("3u = %016lx%016lx\n", t2.data[1][0], t2.data[0][0]);
	if (verbose) printf("v = %016lx%016lx\n", v.data[1][0], v.data[0][0]);
	tecm_addmod104_x8(&t3, &t2, &v, n);				// 3u + v
	if (verbose) printf("3u+v = %016lx%016lx\n", t3.data[1][0], t3.data[0][0]);
	tecm_mulredc104_x8(&t1, &t4, &t3, n, rho);		// a = (v-u)^3 * (3u + v)

	tecm_mulredc104_x8(&t3, &P->X, &t5, n, rho);	// 16*u^3
	tecm_mulredc104_x8(&t4, &v, &t3, n, rho);		// 16*u^3*v

	if (verbose) printf("a = %016lx%016lx\n", t1.data[1][0], t1.data[0][0]);
	if (verbose) printf("d = %016lx%016lx\n", t4.data[1][0], t4.data[0][0]);

	// u holds the denom, t1 holds the numer
	// accomplish the division by multiplying by the modular inverse
	_mm512_storeu_si512(t2.data[0], _mm512_set1_epi64(1));
	_mm512_storeu_si512(t2.data[1], _mm512_set1_epi64(0));
	tecm_mulredc104_x8(&t4, &t2, &t4, n, rho);	// take t4 out of monty rep
	tecm_mulredc104_x8(&t1, &t2, &t1, n, rho);	// take t1 out of monty rep

	for (i = 0; i < 8; i++)
	{
		u104_to_mpz(&t4, i, gmpt);
		u104_to_mpz(n, i, gmpn);

		// a = a4
		// b = n
		// then s == a^-1 mod b
		// and t == b^-1 mod a
		mpz_gcdext(gmpg, gmpt, NULL, gmpt, gmpn);

		if ((mpz_cmp_ui(gmpg, 1) > 0) && (mpz_cmp(gmpg, gmpn) < 0))
		{
			gmp_printf("gcdext found factor %Zd of %Zd while building curves\n",
				gmpg, gmpn);
		}
		
		//mpz_invert(gmpt, gmpt, gmpn);
		//int invcode = mpz_invert(gmpt, gmpt, gmpn);
		//if (invcode == 0)
		//{
		//	printf("inverse does not exist!\n");		// gmpt = t4^-1 mod n
		//}

		u104_to_mpz(&t1, i, gmpn);
		mpz_mul(gmpt, gmpt, gmpn);			// gmpt = t1 * t4 
		u104_to_mpz(n, i, gmpn);
		mpz_tdiv_r(gmpt, gmpt, gmpn);		// gmpt = t1 * t4 % n
		mpz_to_u104(gmpt, &work->s, i);
	}

	to_monty104_x8(&work->s, n);

	mpz_clear(gmpt);
	mpz_clear(gmpn);
	mpz_clear(gmpg);
	mpz_clear(gmpt2);
	return;
}

void tinyecm_x8_list(uint64_t* n, uint64_t* f, uint32_t B1, uint32_t B2, uint32_t curves,
	uint32_t num_in, uint64_t* lcg_state, int verbose)
{
	// attempt to factor n with the elliptic curve method
	// following brent and montgomery's papers, and CP's book
	base_t i, j = 0;
	int curve;
	int tid;
	char* wstr;
	int found = 0;
	int result;
	uint64_t num_found;
	tinyecm_work_x8 work;
	tecm_pt_x8 P;
	monty104_x8_t mdata;
	uint32_t sigma;
	mpz_t gN, gR, gT;

	mpz_init(gN);
	mpz_init(gR);
	mpz_init(gT);

	mpz_set_ull(gR, 1);
	mpz_mul_2exp(gR, gR, 104);

	mpz_set_ull(gN, n[1]);
	mpz_mul_2exp(gN, gN, 52);
	mpz_add_ull(gN, gN, n[0]);

	mpz_mod(gT, gR, gN);

	mpz_invert(gN, gN, gR);

	for (i = 0; i < 8; i++)
	{
		mpz_to_u104(gT, &mdata.one, i);
		mpz_to_u104(gN, &mdata.rho, i);

		work.n.data[0][i] = mdata.n.data[0][i] = n[0];
		work.n.data[1][i] = mdata.n.data[1][i] = n[1];
	}

	f[0] = 1;
	f[1] = 0;
	work.stg1_max = B1;
	work.stg2_max = B2;

	//printf("attempting tecm on %016lx%016lx\n", n[1], n[0]);

	for (curve = 0; curve < curves; curve += 8)
	{
		int done = 0;

		build_curves_104_x8(&P, &mdata, &work, lcg_state, 0);
		ecm_stage1_x8(&mdata, &work, &P);

		for (i = 0; i < 8; i++)
		{
			{
				uint64_t pz[2];
				uint64_t ni[2];
				uint64_t fi[2];

				pz[0] = P.Z.data[1][i] << 52;
				pz[0] |= P.Z.data[0][i];
				pz[1] = P.Z.data[1][i] >> 12;

				ni[0] = mdata.n.data[1][i] << 52;
				ni[0] |= mdata.n.data[0][i];
				ni[1] = mdata.n.data[1][i] >> 12;

				result = check_factor(pz, ni, fi);

				if (result == 1)
				{
					f[0] = fi[0];
					f[1] = fi[1];
					done = 1;
					break;

					//if (done)
					//{
					//	printf("other (ignored) factor %016lx%016lx found in lane %d in stage 1\n",
					//		fi[1], fi[0], i);
					//	printf("original factor is %016lx%016lx\n",
					//		f[1], f[0]);
					//}
					//else
					//{
					//	f[0] = fi[0];
					//	f[1] = fi[1];
					//	if (f[1])
					//		printf("large factor %016lx%016lx found in lane %d in stage 1\n",
					//			f[1], f[0], i);
					//	done = 1;
					//}
				}
			}
		}
		if (done) break;

		if (B2 > B1)
		{
			ecm_stage2_x8(&P, &mdata, &work);

			for (i = 0; i < 8; i++)
			{
				{
					uint64_t acc[2];
					uint64_t ni[2];
					uint64_t fi[2];

					acc[0] = work.stg2acc.data[1][i] << 52;
					acc[0] |= work.stg2acc.data[0][i];
					acc[1] = work.stg2acc.data[1][i] >> 12;

					ni[0] = mdata.n.data[1][i] << 52;
					ni[0] |= mdata.n.data[0][i];
					ni[1] = mdata.n.data[1][i] >> 12;

					result = check_factor(acc, ni, fi);

					if (result == 1)
					{
						f[0] = fi[0];
						f[1] = fi[1];
						done = 1;
						break;
					}
				}
			}
		}
		if (done) break;
	}
done:


	mpz_clear(gN);
	mpz_clear(gR);
	mpz_clear(gT);
	return;
}


#endif


int tecm_dispatch(mpz_t n, mpz_t f, int targetBits, uint64_t* ploc_lcg)
{
	int B1, curves;

	if (targetBits == 0)
	{
		// try fast attempts to find possible small factors.
		{
			B1 = 47;
			curves = 1;
			tinyecm(n, f, B1, B1 * 25, curves, ploc_lcg, 0);
			if (mpz_get_ull(f) > 1)
				return 1;
		}
		{
			B1 = 70;
			curves = 1;
			tinyecm(n, f, B1, 25 * B1, curves, ploc_lcg, 0);
			if (mpz_get_ull(f) > 1)
				return 1;
		}
		{
			B1 = 125;
			curves = 1;
			tinyecm(n, f, B1, 25 * B1, curves, ploc_lcg, 0);
			if (mpz_get_ull(f) > 1)
				return 1;
		}
	}

	if (targetBits <= 20)
	{
		B1 = 27;
		curves = 32;
		tinyecm(n, f, B1, 25 * B1, curves, ploc_lcg, 0);
	}
	else if (targetBits <= 22)
	{
		B1 = 47;
		curves = 32;
		tinyecm(n, f, B1, 25 * B1, curves, ploc_lcg, 0);
	}
	else if (targetBits <= 24)
	{
		B1 = 70;
		curves = 32;
		tinyecm(n, f, B1, 25 * B1, curves, ploc_lcg, 0);
	}
	else if (targetBits <= 26)
	{
		B1 = 85;
		curves = 32;
		tinyecm(n, f, B1, 25 * B1, curves, ploc_lcg, 0);
	}
	else if (targetBits <= 29)
	{
		B1 = 125;
		curves = 32;
		tinyecm(n, f, B1, 25 * B1, curves, ploc_lcg, 0);
	}
	else if (targetBits <= 31)
	{
		B1 = 165;
		curves = 42;
		tinyecm(n, f, B1, 25 * B1, curves, ploc_lcg, 0);
	}
	else //if (targetBits <= 32)
	{
		B1 = 205;
		curves = 42;
		tinyecm(n, f, B1, 25 * B1, curves, ploc_lcg, 0);
	}

	if (mpz_get_ull(f) > 1)
		return 1;
	else
		return 0;
}


static int tecm_get_bits(mpz_t n)
{
	return mpz_sizeinbase(n, 2);
}

// getfactor_tecm() returns 0 if unable to find a factor of n,
// Otherwise it returns 1 and a factor of n in argument f.
// 
// if the input is known to have no small factors, set is_arbitrary=0, 
// otherwise, set is_arbitrary=1 and a few curves targetting small factors
// will be run prior to the standard sequence of curves for the input size.
//  
// Prior to your first call of getfactor_tecm(), set *pran = 0  (or set it to
// some other arbitrary value); after that, don't change *pran.
// FYI: *pran is used within this file by a random number generator, and it
// holds the current value of a pseudo random sequence.  Your first assigment
// to *pran seeds the sequence, and after seeding it you don't want to
// change *pran, since that would restart the sequence.

int getfactor_tpm1(mpz_t n, mpz_t f, uint32_t b1)
{
	if (mpz_even_p(n))
	{
		mpz_set_ull(f, 2);
		return 1;
	}

	tinypm1(n, f, b1, 25 * b1);
	if (mpz_get_ull(f) > 1)
		return 1;
	else
		return 0;
}

#ifdef USE_AVX512F

void tecm_testmath()
{
	int num = 1000000;
	int bits = 32;
	vec_u104_t x, y, n, r, z;
	int i, j;
	mpz_t gx, gy, gm, gi, gz, gn, gt, gr, glomask;
	gmp_randstate_t gmp_randstate;
	gmp_randinit_default(gmp_randstate);

	mpz_init(gx);
	mpz_init(gy);
	mpz_init(gn);
	mpz_init(gt);
	mpz_init(gr);
	mpz_init(gz);
	mpz_init(gm);
	mpz_init(gi);
	mpz_init(glomask);

	mpz_set_ull(gr, 1);
	mpz_mul_2exp(gr, gr, 104);
	mpz_sub_ui(glomask, gr, 1);

	for (bits = 30; bits <= 100; bits += 5)
	{
		printf("bits = %u\n", bits);
		printf("\tmodmul\n");
		for (i = 0; i < num; i++)
		{
			for (j = 0; j < 8; j++)
			{
				mpz_urandomb(gx, gmp_randstate, bits);
				mpz_setbit(gx, bits - 1);
				mpz_urandomb(gy, gmp_randstate, bits);
				mpz_setbit(gy, bits - 1);
				mpz_urandomb(gn, gmp_randstate, bits);
				mpz_setbit(gn, bits - 1);
				if (mpz_even_p(gn))
					mpz_sub_ui(gn, gn, 1);
				mpz_invert(gt, gn, gr);

				mpz_to_u104(gx, &x, j);
				mpz_to_u104(gy, &y, j);
				mpz_to_u104(gn, &n, j);
				mpz_to_u104(gt, &r, j);
			}

			to_monty104_x8(&x, &n);
			to_monty104_x8(&y, &n);
			tecm_mulredc104_x8(&z, &x, &y, &n, &r);

			for (j = 0; j < 8; j++)
			{
				u104_to_mpz(&x, j, gx);
				u104_to_mpz(&y, j, gy);
				u104_to_mpz(&n, j, gn);
				u104_to_mpz(&z, j, gz);

				mpz_invert(gi, gn, gr);

				// get the low 104 bits of m = Tlo * invN
				// get the high 104 bits of m * N
				// result is (T - mN) >> 104
				// unless T_hi < mN_hi
				// then return (T_hi + N - mN_hi).
				mpz_mul(gt, gx, gy);
				mpz_and(gm, gt, glomask);
				mpz_mul(gm, gm, gi);
				mpz_and(gm, gm, glomask);
				mpz_mul(gm, gm, gn);
				mpz_tdiv_q_2exp(gt, gt, 104);
				mpz_tdiv_q_2exp(gm, gm, 104);
				if (mpz_cmp(gm, gt) > 0)
				{
					mpz_add(gt, gt, gn);
					mpz_sub(gt, gt, gm);
				}
				else
				{
					mpz_sub(gt, gt, gm);
				}

				if (mpz_cmp(gt, gz) != 0)
				{
					gmp_printf("mul failed: x=%Zd, y=%Zd, n=%Zd, r=%Zd, z=%Zd != %Zd\n",
						gx, gy, gn, gr, gt, gz);
					exit(1);
				}
			}
		}
		printf("\tmodadd\n");
		for (i = 0; i < num; i++)
		{
			for (j = 0; j < 8; j++)
			{
				mpz_urandomb(gx, gmp_randstate, bits);
				mpz_setbit(gx, bits - 1);
				mpz_urandomb(gy, gmp_randstate, bits);
				mpz_setbit(gy, bits - 1);
				mpz_urandomb(gn, gmp_randstate, bits);
				mpz_setbit(gn, bits - 1);
				if (mpz_even_p(gn))
					mpz_sub_ui(gn, gn, 1);

				mpz_to_u104(gx, &x, j);
				mpz_to_u104(gy, &y, j);
				mpz_to_u104(gn, &n, j);
			}

			to_monty104_x8(&x, &n);
			to_monty104_x8(&y, &n);
			tecm_addmod104_x8(&z, &x, &y, &n);

			for (j = 0; j < 8; j++)
			{
				u104_to_mpz(&x, j, gx);
				u104_to_mpz(&y, j, gy);
				u104_to_mpz(&n, j, gn);
				u104_to_mpz(&z, j, gz);

				mpz_add(gt, gx, gy);

				if (mpz_cmp(gt, gn) >= 0)
				{
					mpz_sub(gt, gt, gn);
				}

				if (mpz_cmp(gt, gz) != 0)
				{
					gmp_printf("add failed: x=%Zd, y=%Zd, n=%Zd, r=%Zd, z=%Zd != %Zd\n",
						gx, gy, gn, gr, gt, gz);
					exit(1);
				}
			}
		}
		printf("\tmodsub\n");
		for (i = 0; i < num; i++)
		{
			for (j = 0; j < 8; j++)
			{
				mpz_urandomb(gx, gmp_randstate, bits);
				mpz_setbit(gx, bits - 1);
				mpz_urandomb(gy, gmp_randstate, bits);
				mpz_setbit(gy, bits - 1);
				mpz_urandomb(gn, gmp_randstate, bits);
				mpz_setbit(gn, bits - 1);
				if (mpz_even_p(gn))
					mpz_sub_ui(gn, gn, 1);

				mpz_to_u104(gx, &x, j);
				mpz_to_u104(gy, &y, j);
				mpz_to_u104(gn, &n, j);
			}

			to_monty104_x8(&x, &n);
			to_monty104_x8(&y, &n);
			tecm_submod104_x8(&z, &x, &y, &n);

			for (j = 0; j < 8; j++)
			{
				u104_to_mpz(&x, j, gx);
				u104_to_mpz(&y, j, gy);
				u104_to_mpz(&n, j, gn);
				u104_to_mpz(&z, j, gz);

				if (mpz_cmp(gy, gx) > 0)
				{
					mpz_add(gt, gx, gn);
					mpz_sub(gt, gt, gy);
				}
				else
				{
					mpz_sub(gt, gx, gy);
				}

				if (mpz_cmp(gt, gz) != 0)
				{
					gmp_printf("sub failed: x=%Zd, y=%Zd, n=%Zd, r=%Zd, z=%Zd != %Zd\n",
						gx, gy, gn, gr, gt, gz);
					exit(1);
				}
			}
		}
	}




}

int tecm_dispatch_x8_list(mpz_t gn, mpz_t gf, int targetBits, uint64_t* ploc_lcg)
{
	int B1, curves;

	uint64_t n[2];
	uint64_t f[2];
	int i;

	n[0] = gn->_mp_d[0] & 0x000fffffffffffffull; //mpz_get_ull(gn) & 0x000fffffffffffffull;
	mpz_tdiv_q_2exp(gf, gn, 52);
	n[1] = gf->_mp_d[0] & 0x000fffffffffffffull; //mpz_get_ull(gf) & 0x000fffffffffffffull;

	//printf("in tecm_dispatch, N = %013llx%013llx\n", n[1], n[0]);

	mpz_set_ull(gf, 1);

	if (targetBits <= 20)
	{
		B1 = 27;
		curves = 32;
		tinyecm_x8_list(n, f, B1, 25 * B1, curves, curves, ploc_lcg, 0);
	}
	else if (targetBits <= 22)
	{
		B1 = 47;
		curves = 32;
		tinyecm_x8_list(n, f, B1, 25 * B1, curves, curves, ploc_lcg, 0);
	}
	else if (targetBits <= 24)
	{
		B1 = 70;
		curves = 32;
		tinyecm_x8_list(n, f, B1, 25 * B1, curves, curves, ploc_lcg, 0);
	}
	else if (targetBits <= 26)
	{
		B1 = 85;
		curves = 32;
		tinyecm_x8_list(n, f, B1, 25 * B1, curves, curves, ploc_lcg, 0);
	}
	else if (targetBits <= 29)
	{
		B1 = 125;
		curves = 32;
		tinyecm_x8_list(n, f, B1, 25 * B1, curves, curves, ploc_lcg, 0);
	}
	else if (targetBits <= 31)
	{
		B1 = 165;
		curves = 48;
		tinyecm_x8_list(n, f, B1, 25 * B1, curves, curves, ploc_lcg, 0);
	}
	else //if (targetBits <= 32)
	{
		B1 = 205;
		curves = 48;
		tinyecm_x8_list(n, f, B1, 25 * B1, curves, curves, ploc_lcg, 0);
	}

	mpz_set_ull(gf, f[1]);
	mpz_mul_2exp(gf, gf, 64);
	mpz_add_ull(gf, gf, f[0]);

	if (mpz_get_ull(gf) > 1)
		return 1;
	else
		return 0;
}

int getfactor_tecm(mpz_t n, mpz_t f, int target_bits, uint64_t* pran)
{
	dbias = _mm512_castsi512_pd(tecm_set64(0x4670000000000000ULL));
	vbias1 = tecm_set64(0x4670000000000000ULL);
	vbias2 = tecm_set64(0x4670000000000001ULL);
	lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);

	//tecm_testmath();
	//return;

	if (mpz_sizeinbase(n, 2) > 104)
		printf("warning: n is too large (%d bits)\n", mpz_sizeinbase(n, 2));

	return tecm_dispatch_x8_list(n, f, target_bits, pran);
	//return tecm_dispatch(n, f, target_bits, pran);
}



#else

int getfactor_tecm(mpz_t n, mpz_t f, int target_bits, uint64_t* pran)
{
	if (mpz_even_p(n))
	{
		mpz_set_ull(f, 2);
		return 1;
	}

	return tecm_dispatch(n, f, target_bits, pran);
}

#endif
