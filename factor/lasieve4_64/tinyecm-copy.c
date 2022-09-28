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

#define D 120

//#define DEBUG 1
#define base_t uint64_t
#define base_digits 2

typedef struct
{
	uint64_t base[2];
} u128_t;

typedef struct
{
	base_t X[base_digits];
	base_t Z[base_digits];
} tinyecm_pt;

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

    mpz_set_ui(m, x[1]);
    mpz_mul_2exp(m, m, 64);
    mpz_add_ui(m, m, x[0]);

    mpz_set_ui(n, mdata->n[1]);
    mpz_mul_2exp(n, n, 64);
    mpz_add_ui(n, n, mdata->n[0]);

    // implied R = 2^128
    mpz_mul_2exp(m, m, 128);
    mpz_mod(m, m, n);

    x[0] = mpz_get_ui(m);
    mpz_tdiv_q_2exp(m, m, 64);
    x[1] = mpz_get_ui(m);

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
#if defined( GCC_ASM64X ) && !defined(__MINGW32__)

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
	mpz_set_ui(out, in[1]);
	mpz_mul_2exp(out, out, 64);
	mpz_add_ui(out, out, in[0]);
	return;
}

void mpz_to_u128(mpz_t in, uint64_t *out)
{
    out[0] = mpz_get_ui(in);
    mpz_tdiv_q_2exp(in, in, 64);
    out[1] = mpz_get_ui(in);

    // restore input
    mpz_mul_2exp(in, in, 64);
    mpz_add_ui(in, in, out[0]);

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

	mpz_set_ui(f, 1);
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
		//        uecm_uprac(mdata, work, P, 47, 0.548409048446403258);
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
			mpz_set_ui(gmp_f, 0);
			status = 0;
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

void tinyecm_test(int sizeb, int num, int type)
{
	u128_t *inputs;
	int i;

	if (sizeb > 128)
	{
		printf("size must be <= 128 bits\n");
		return;
	}

	if (num > (1 << 30))
	{
		printf("please be reasonable\n");
		return;
	}

	inputs = (u128_t *)malloc(num * sizeof(u128_t));

	// build the requested input type, quantity and size
	for (i = 0; i < num; i++)
	{


	}

	// test them, building up timings in globals


	// report results

	return;
}


int tecm_dispatch(mpz_t n, mpz_t f, int targetBits, uint64_t* ploc_lcg)
{
	int B1, curves;

	//uecm_stage2_pair(70, 180, 120, 120, 30, primes);
	//uecm_stage2_pair(180, 25 * 70, 150, 120, 30, primes);

	//uecm_stage2_pair(47, 150, 90, 120, 30, primes);
	//uecm_stage2_pair(150, 25 * 47, 150, 120, 30, primes);


	//uecm_stage2_pair(30, 150, 90, 120, 30, primes);
	//uecm_stage2_pair(150, 25 * 30, 150, 120, 30, primes);
	//exit(0);

	if (targetBits == 0)
	{
		// try fast attempts to find possible small factors.
		{
			B1 = 47;
			curves = 1;
			tinyecm(n, f, B1, B1 * 25, curves, ploc_lcg, 0);
			if (mpz_get_ui(f) > 1)
				return 1;
		}
		{
			B1 = 70;
			curves = 1;
			tinyecm(n, f, B1, 25 * B1, curves, ploc_lcg, 0);
			if (mpz_get_ui(f) > 1)
				return 1;
		}
		{
			B1 = 125;
			curves = 1;
			tinyecm(n, f, B1, 25 * B1, curves, ploc_lcg, 0);
			if (mpz_get_ui(f) > 1)
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

	if (mpz_get_ui(f) > 1)
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
int getfactor_tecm(mpz_t n, mpz_t f, int target_bits, uint64_t* pran)
{
	if (mpz_even_p(n))
	{
		mpz_set_ui(f, 2);
		return 1;
	}

	return tecm_dispatch(n, f, target_bits, pran);
}

