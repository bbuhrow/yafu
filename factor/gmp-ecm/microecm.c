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

#include "monty.h"
#include "stdint.h"
#include "ytools.h"

#define D 60

//#define DEBUG 1
typedef struct
{
	uint64_t X;
	uint64_t Z;
} uecm_pt;

typedef struct
{
	uint64_t sum1;
	uint64_t diff1;
	uint64_t sum2;
	uint64_t diff2;
	uint64_t tt1;
	uint64_t tt2;
	uint64_t tt3;
	uint64_t tt4;
	uint64_t tt5;
	uint64_t s;
	uint64_t n;
	uecm_pt pt1;
	uecm_pt pt2;
	uecm_pt pt3;
	uecm_pt pt4;
	uecm_pt pt5;
	uint32_t sigma;

	uecm_pt Pa;
	uecm_pt Pd;
	uecm_pt Pad;
	uecm_pt Pb[20];
	uint64_t Paprod;
	uint64_t Pbprod[20];

	uint64_t stg2acc;
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

} uecm_work;

static const uint32 map[61] = {
	0, 1, 2, 0, 0, 0, 0, 3, 0, 0,
	0, 4, 0, 5, 0, 0, 0, 6, 0, 7,
	0, 0, 0, 8, 0, 0, 0, 0, 0, 9,
	0, 10, 0, 0, 0, 0, 0, 11, 0, 0,
	0, 12, 0, 13, 0, 0, 0, 14, 0, 15,
	0, 0, 0, 16, 0, 0, 0, 0, 0, 17,
	18 };

static const uint32_t primes[801] = {
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

static uint8 marks[D];
static uint8 nmarks[D];

#ifdef _MSC_VERasdf

#include <immintrin.h>
#include <intrin.h>

__inline uint64 u64div(uint64 c, uint64 n)
{
    uint64 r;
    mpz_t a;
    mpz_init(a);
    mpz_set_ui(a, c);
    mpz_mul_2exp(a, a, 64);
    r = mpz_tdiv_ui(a, n);
    mpz_clear(a);

    // first available in Visual Studio 2019
    //_udiv128(c, 0, n, &r);

    return r;
}

__inline uint64 mulredc(uint64 x, uint64 y, uint64 n, uint64 nhat)
{
    uint64 th, tl, u, ah, al;
    tl = _umul128(x, y, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    tl = _addcarry_u64(0, al, tl, &al);
    th = _addcarry_u64(tl, th, ah, &x);
    x = ((x >= n) || th) ? x - n : x;
    return x;
}

__inline uint64 sqrredc(uint64 x, uint64 n, uint64 nhat)
{
    uint64 th, tl, u, ah, al;
    tl = _umul128(x, x, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    tl = _addcarry_u64(0, al, tl, &al);
    th = _addcarry_u64(tl, th, ah, &x);
    x = ((x >= n) || th) ? x - n : x;
    return x;
}


__inline uint64 submod(uint64 a, uint64 b, uint64 n)
{
    uint64 r0;
    if (_subborrow_u64(0, a, b, &r0))
        r0 += n;
    return r0;
}

__inline uint64 addmod(uint64 x, uint64 y, uint64 n)
{
    uint64 r;
    uint8 c = _addcarry_u64(0, x, y, &r);

    if (c || (r >= n))
        r -= n;
    return r;
}


#endif

// local functions
void uadd(uint64_t rho, uecm_work *work, uecm_pt *P1, uecm_pt *p2,
	uecm_pt *Pin, uecm_pt *Pout);
void udup(uint64_t rho, uecm_work *work, uint64 insum, uint64 indiff, uecm_pt *P);
void uprac(uint64_t rho, uecm_work *work, uecm_pt *P, uint64_t c, double v);
int ucheck_factor(uint64 Z, uint64 n, uint64 * f);
void ubuild(uecm_pt *P, uint64_t rho, uecm_work *work, uint32_t sigma, int verbose);
void uecm_stage1(uint64_t rho, uecm_work *work, uecm_pt *P);
void uecm_stage2(uecm_pt *P, uint64_t rho, uecm_work *work);
uint32 lcg_rand_32B(uint32 lower, uint32 upper);

#define INV_2_POW_32 0.00000000023283064365386962890625

uint32 lcg_rand_32B(uint32 lower, uint32 upper)
{
    LCG_STATE = 6364136223846793005ULL * LCG_STATE + 1442695040888963407ULL;
    return lower + (fp_digit)(
        (double)(upper - lower) * (double)(LCG_STATE >> 32) * INV_2_POW_32);
}

void uadd(uint64_t rho, uecm_work *work, uecm_pt *P1, uecm_pt *P2, 
	uecm_pt *Pin, uecm_pt *Pout)
{
	// compute:
	//x+ = z- * [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
	//z+ = x- * [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
	// where:
	//x- = original x
	//z- = original z
	// given the sums and differences of the original points (stored in work structure).
	work->diff1 = submod(P1->X, P1->Z, work->n);
	work->sum1 = addmod(P1->X, P1->Z, work->n);
	work->diff2 = submod(P2->X, P2->Z, work->n);
	work->sum2 = addmod(P2->X, P2->Z, work->n);

	work->tt1 = mulredc(work->diff1, work->sum2, work->n, rho);	//U
	work->tt2 = mulredc(work->sum1, work->diff2, work->n, rho);	//V

	work->tt3 = addmod(work->tt1, work->tt2, work->n);
	work->tt4 = submod(work->tt1, work->tt2, work->n);
	work->tt1 = sqrredc(work->tt3, work->n, rho);	//(U + V)^2
	work->tt2 = sqrredc(work->tt4, work->n, rho);	//(U - V)^2

	if (Pin == Pout)
	{
		uint64_t tmp;
		Pout->Z = mulredc(work->tt1, Pin->Z, work->n, rho);		//Z * (U + V)^2
		Pout->X = mulredc(work->tt2, Pin->X, work->n, rho);		//x * (U - V)^2
		tmp = Pout->Z;
		Pout->Z = Pout->X;
		Pout->X = tmp;
	}
	else
	{
		Pout->X = mulredc(work->tt1, Pin->Z, work->n, rho);		//Z * (U + V)^2
		Pout->Z = mulredc(work->tt2, Pin->X, work->n, rho);		//x * (U - V)^2
	}
	work->stg1Add++;
	return;
}

void udup(uint64_t rho, uecm_work *work, 
	uint64 insum, uint64 indiff, uecm_pt *P)
{
	work->tt1 = sqrredc(indiff, work->n, rho);			// U=(x1 - z1)^2
	work->tt2 = sqrredc(insum, work->n, rho);			// V=(x1 + z1)^2
	P->X = mulredc(work->tt1, work->tt2, work->n, rho);			// x=U*V

	work->tt3 = submod(work->tt2, work->tt1, work->n);			// w = V-U
	work->tt2 = mulredc(work->tt3, work->s, work->n, rho);		// w = (A+2)/4 * w
	work->tt2 = addmod(work->tt2, work->tt1, work->n);			// w = w + U
	P->Z = mulredc(work->tt2, work->tt3, work->n, rho);			// Z = w*(V-U)
	work->stg1Doub++;
	return;
}

void uprac70(uint64_t rho, uecm_work *work, uecm_pt *P)
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

    for (i = 0; i < 116; i++)
    {
        if (steps[i] == 0)
        {
            work->pt1.X = work->pt2.X = work->pt3.X = P->X;
            work->pt1.Z = work->pt2.Z = work->pt3.Z = P->Z;

            d1 = submod(work->pt1.X, work->pt1.Z, work->n);
            s1 = addmod(work->pt1.X, work->pt1.Z, work->n);
            udup(rho, work, s1, d1, &work->pt1);
        }
        else if (steps[i] == 3)
        {
            // integrate step 4 followed by swap(1,2)
            uadd(rho, work, &work->pt2, &work->pt1, &work->pt3, &work->pt4);		// T = B + A (C)

            swp = work->pt1.X;
            work->pt1.X = work->pt4.X;
            work->pt4.X = work->pt3.X;
            work->pt3.X = work->pt2.X;
            work->pt2.X = swp;
            swp = work->pt1.Z;
            work->pt1.Z = work->pt4.Z;
            work->pt4.Z = work->pt3.Z;
            work->pt3.Z = work->pt2.Z;
            work->pt2.Z = swp;
        }
        else if (steps[i] == 4)
        {
            uadd(rho, work, &work->pt2, &work->pt1, &work->pt3, &work->pt4);		// T = B + A (C)

            swp = work->pt2.X;
            work->pt2.X = work->pt4.X;
            work->pt4.X = work->pt3.X;
            work->pt3.X = swp;
            swp = work->pt2.Z;
            work->pt2.Z = work->pt4.Z;
            work->pt4.Z = work->pt3.Z;
            work->pt3.Z = swp;
        }
        else if (steps[i] == 5)
        {
            d2 = submod(work->pt1.X, work->pt1.Z, work->n);
            s2 = addmod(work->pt1.X, work->pt1.Z, work->n);

            uadd(rho, work, &work->pt2, &work->pt1, &work->pt3, &work->pt2);		// B = B + A (C)
            udup(rho, work, s2, d2, &work->pt1);		// A = 2A
        }
        else if (steps[i] == 6)
        {
            uadd(rho, work, &work->pt1, &work->pt2, &work->pt3, P);		// A = A + B (C)
        }

    }

    return;

}

void uprac85(uint64_t rho, uecm_work *work, uecm_pt *P)
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

    for (i = 0; i < 146; i++)
    {
        if (steps[i] == 0)
        {
            work->pt1.X = work->pt2.X = work->pt3.X = P->X;
            work->pt1.Z = work->pt2.Z = work->pt3.Z = P->Z;

            d1 = submod(work->pt1.X, work->pt1.Z, work->n);
            s1 = addmod(work->pt1.X, work->pt1.Z, work->n);
            udup(rho, work, s1, d1, &work->pt1);
        }
        else if (steps[i] == 3)
        {
            // integrate step 4 followed by swap(1,2)
            uadd(rho, work, &work->pt2, &work->pt1, &work->pt3, &work->pt4);		// T = B + A (C)

            swp = work->pt1.X;
            work->pt1.X = work->pt4.X;
            work->pt4.X = work->pt3.X;
            work->pt3.X = work->pt2.X;
            work->pt2.X = swp;
            swp = work->pt1.Z;
            work->pt1.Z = work->pt4.Z;
            work->pt4.Z = work->pt3.Z;
            work->pt3.Z = work->pt2.Z;
            work->pt2.Z = swp;
        }
        else if (steps[i] == 4)
        {
            uadd(rho, work, &work->pt2, &work->pt1, &work->pt3, &work->pt4);		// T = B + A (C)

            swp = work->pt2.X;
            work->pt2.X = work->pt4.X;
            work->pt4.X = work->pt3.X;
            work->pt3.X = swp;
            swp = work->pt2.Z;
            work->pt2.Z = work->pt4.Z;
            work->pt4.Z = work->pt3.Z;
            work->pt3.Z = swp;
        }
        else if (steps[i] == 5)
        {
            d2 = submod(work->pt1.X, work->pt1.Z, work->n);
            s2 = addmod(work->pt1.X, work->pt1.Z, work->n);

            uadd(rho, work, &work->pt2, &work->pt1, &work->pt3, &work->pt2);		// B = B + A (C)
            udup(rho, work, s2, d2, &work->pt1);		// A = 2A
        }
        else if (steps[i] == 6)
        {
            uadd(rho, work, &work->pt1, &work->pt2, &work->pt3, P);		// A = A + B (C)
        }

    }

    return;

}

void uprac(uint64_t rho, uecm_work *work, uecm_pt *P, uint64_t c, double v)
{
	uint64_t d, e, r;
	int i;
	uint64 s1, s2, d1, d2;
	int shift = 0;
	uint64_t swp;

	while ((c & 1) == 0)
	{
		shift++;
		c >>= 1;
	}
	
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
	work->pt1.X = work->pt2.X = work->pt3.X = P->X;
	work->pt1.Z = work->pt2.Z = work->pt3.Z = P->Z;

	d1 = submod(work->pt1.X, work->pt1.Z, work->n);
	s1 = addmod(work->pt1.X, work->pt1.Z, work->n);

	// point2 is [2]P
	udup(rho, work, s1, d1, &work->pt1);

	while (d != e)
	{
		if (d < e)
		{
			r = d;
			d = e;
			e = r;
			swp = work->pt1.X;
			work->pt1.X = work->pt2.X;
			work->pt2.X = swp;
			swp = work->pt1.Z;
			work->pt1.Z = work->pt2.Z;
			work->pt2.Z = swp;
		}
		if (d - e <= e / 4 && ((d + e) % 3) == 0)
		{
			d = (2 * d - e) / 3;
			e = (e - d) / 2;

			uadd(rho, work, &work->pt1, &work->pt2, &work->pt3, &work->pt4); // T = A + B (C)
			uadd(rho, work, &work->pt4, &work->pt1, &work->pt2, &work->pt5); // T2 = T + A (B)
			uadd(rho, work, &work->pt2, &work->pt4, &work->pt1, &work->pt2); // B = B + T (A)

			swp = work->pt1.X;
			work->pt1.X = work->pt5.X;
			work->pt5.X = swp;
			swp = work->pt1.Z;
			work->pt1.Z = work->pt5.Z;
			work->pt5.Z = swp;
		}
		else if (d - e <= e / 4 && (d - e) % 6 == 0)
		{
			d = (d - e) / 2;

			d1 = submod(work->pt1.X, work->pt1.Z, work->n);
			s1 = addmod(work->pt1.X, work->pt1.Z, work->n);

			uadd(rho, work, &work->pt1, &work->pt2, &work->pt3, &work->pt2);		// B = A + B (C)
			udup(rho, work, s1, d1, &work->pt1);		// A = 2A
		}
		else if ((d + 3) / 4 <= e)
		{
			d -= e;

			uadd(rho, work, &work->pt2, &work->pt1, &work->pt3, &work->pt4);		// T = B + A (C)

			swp = work->pt2.X;
			work->pt2.X = work->pt4.X;
			work->pt4.X = work->pt3.X;
			work->pt3.X = swp;
			swp = work->pt2.Z;
			work->pt2.Z = work->pt4.Z;
			work->pt4.Z = work->pt3.Z;
			work->pt3.Z = swp;
		}
		else if ((d + e) % 2 == 0)
		{
			d = (d - e) / 2;

			d2 = submod(work->pt1.X, work->pt1.Z, work->n);
			s2 = addmod(work->pt1.X, work->pt1.Z, work->n);

			uadd(rho, work, &work->pt2, &work->pt1, &work->pt3, &work->pt2);		// B = B + A (C)
			udup(rho, work, s2, d2, &work->pt1);		// A = 2A
		}
		else if (d % 2 == 0)
		{
			d /= 2;

			d2 = submod(work->pt1.X, work->pt1.Z, work->n);
			s2 = addmod(work->pt1.X, work->pt1.Z, work->n);

			uadd(rho, work, &work->pt3, &work->pt1, &work->pt2, &work->pt3);		// C = C + A (B)
			udup(rho, work, s2, d2, &work->pt1);		// A = 2A
		}
		else if (d % 3 == 0)
		{
			d = d / 3 - e;

			d1 = submod(work->pt1.X, work->pt1.Z, work->n);
			s1 = addmod(work->pt1.X, work->pt1.Z, work->n);

			udup(rho, work, s1, d1, &work->pt4);		// T = 2A
			uadd(rho, work, &work->pt1, &work->pt2, &work->pt3, &work->pt5);		// T2 = A + B (C)
			uadd(rho, work, &work->pt4, &work->pt1, &work->pt1, &work->pt1);		// A = T + A (A)
			uadd(rho, work, &work->pt4, &work->pt5, &work->pt3, &work->pt4);		// T = T + T2 (C)

			swp = work->pt3.X;
			work->pt3.X = work->pt2.X;
			work->pt2.X = work->pt4.X;
			work->pt4.X = swp;
			swp = work->pt3.Z;
			work->pt3.Z = work->pt2.Z;
			work->pt2.Z = work->pt4.Z;
			work->pt4.Z = swp;

		}
		else if ((d + e) % 3 == 0)
		{
			d = (d - 2 * e) / 3;

			uadd(rho, work, &work->pt1, &work->pt2, &work->pt3, &work->pt4);		// T = A + B (C)


			d2 = submod(work->pt1.X, work->pt1.Z, work->n);
			s2 = addmod(work->pt1.X, work->pt1.Z, work->n);
			uadd(rho, work, &work->pt4, &work->pt1, &work->pt2, &work->pt2);		// B = T + A (B)
			udup(rho, work, s2, d2, &work->pt4);		// T = 2A
			uadd(rho, work, &work->pt1, &work->pt4, &work->pt1, &work->pt1);		// A = A + T (A) = 3A
		}
		else if ((d - e) % 3 == 0)
		{
			d = (d - e) / 3;

			uadd(rho, work, &work->pt1, &work->pt2, &work->pt3, &work->pt4);		// T = A + B (C)

			d2 = submod(work->pt1.X, work->pt1.Z, work->n);
			s2 = addmod(work->pt1.X, work->pt1.Z, work->n);
			uadd(rho, work, &work->pt3, &work->pt1, &work->pt2, &work->pt3);		// C = C + A (B)

			swp = work->pt2.X;
			work->pt2.X = work->pt4.X;
			work->pt4.X = swp;
			swp = work->pt2.Z;
			work->pt2.Z = work->pt4.Z;
			work->pt4.Z = swp;

			udup(rho, work, s2, d2, &work->pt4);		// T = 2A
			uadd(rho, work, &work->pt1, &work->pt4, &work->pt1, &work->pt1);		// A = A + T (A) = 3A
		}
		else
		{
			e /= 2;

			d2 = submod(work->pt2.X, work->pt2.Z, work->n);
			s2 = addmod(work->pt2.X, work->pt2.Z, work->n);

			uadd(rho, work, &work->pt3, &work->pt2, &work->pt1, &work->pt3);		// C = C + B (A)
			udup(rho, work, s2, d2, &work->pt2);		// B = 2B
		}
	}

	uadd(rho, work, &work->pt1, &work->pt2, &work->pt3, P);		// A = A + B (C)

	for (i = 0; i < shift; i++)
	{
		d1 = submod(P->X, P->Z, work->n);
		s1 = addmod(P->X, P->Z, work->n);
		udup(rho, work, s1, d1, P);		// P = 2P
	}

	return;

}

uint64_t modinv_64(uint64_t a, uint64_t p) {

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

	if (parity == 0)
		return ps1;
	else
		return p - ps1;
}

void ubuild(uecm_pt *P, uint64_t rho, 
	uecm_work *work, uint32_t sigma, int verbose)
{
	uint64_t t1, t2, t3, t4;
	uint64_t u, v, n;

	n = work->n;

	if (sigma == 0)
	{
		work->sigma = lcg_rand_32B(7, (uint32)-1);
	}
	else
	{
		work->sigma = sigma;
	}
    sigma = work->sigma; // = 1073169389;

    u = (uint64)sigma % n;
	u = u64div(u, n);               // to_monty(sigma)

    //printf("sigma = %" PRIu64 ", u = %" PRIu64 ", n = %" PRIu64 "\n", sigma, u, n);
	
    v = addmod(u, u, n);
    v = addmod(v, v, n);            // v = 4*sigma

    //printf("v = %" PRIu64 "\n", v);

	u = mulredc(u, u, n, rho);
	t1 = 5;
	t1 = u64div(t1, n);             // to_monty(5)

    //printf("monty(5) = %" PRIu64 "\n", t1);

	u = submod(u, t1, n);			// u = sigma^2 - 5

    //printf("u = %" PRIu64 "\n", u);

	t1 = mulredc(u, u, n, rho);
	P->X = mulredc(t1, u, n, rho);	// x = u^3

	t1 = mulredc(v, v, n, rho);
	P->Z = mulredc(t1, v, n, rho);	// z = v^3

	//compute parameter A
	t1 = submod(v, u, n);			// (v - u)
	t2 = mulredc(t1, t1, n, rho);
	t4 = mulredc(t2, t1, n, rho);	// (v - u)^3

	t1 = 3;
	t1 = u64div(t1, n);             // to_monty(3)
	t2 = mulredc(t1, u, n, rho);	// 3u
	t3 = addmod(t2, v, n);			// 3u + v

	t1 = mulredc(t3, t4, n, rho);	// a = (v-u)^3 * (3u + v)

	t2 = 16;
	t2 = u64div(t2, n);             // to_monty(16)
	t3 = mulredc(P->X, t2, n, rho);	// 16*u^3
	t4 = mulredc(t3, v, n, rho);	// 16*u^3*v

	// u holds the denom, t1 holds the numer
	// accomplish the division by multiplying by the modular inverse
	t2 = 1;
	t4 = mulredc(t4, t2, n, rho);	// take t4 out of monty rep
	t1 = mulredc(t1, t2, n, rho);	// take t1 out of monty rep

	t3 = modinv_64(t4, n);
#ifdef _MSC_VER
    {
        mpz_t gmpt;
        mpz_init(gmpt);
        mpz_set_ui(gmpt, t3);
        mpz_mul_ui(gmpt, gmpt, t1);
        work->s = mpz_tdiv_ui(gmpt, n);
        mpz_clear(gmpt);
    }
#else
	spMulMod(t3, t1, n, &work->s);
#endif
	work->s = u64div(work->s, n);   // to_monty(s)

	return;
}

void microecm(uint64_t n, uint64_t *f, uint32 B1, uint32 B2, 
	uint32 curves, int verbose)
{
	//attempt to factor n with the elliptic curve method
	//following brent and montgomery's papers, and CP's book
	uint64_t retval;
	uint64_t i, j;
	int curve;
	int tid;
	char *wstr;
	int found = 0;
	int result;
	uint64_t num_found;
	uecm_work work;
	uecm_pt P;
	uint32 sigma;
	uint64_t rho, x, tmp1;

	x = (((n + 2) & 4) << 1) + n; // here x*a==1 mod 2**4
	x *= 2 - n * x;               // here x*a==1 mod 2**8
	x *= 2 - n * x;               // here x*a==1 mod 2**16
	x *= 2 - n * x;               // here x*a==1 mod 2**32         
	x *= 2 - n * x;               // here x*a==1 mod 2**64
	rho = (uint64)0 - x;
	work.n = n;
	work.stg1_max = B1;
	work.stg2_max = B2;

	*f = 1;
	for (curve = 0; curve < curves; curve++)
	{
		uint64_t p;

		work.stg1Add = 0;
		work.stg1Doub = 0;
		work.last_pid = 0;

		if (verbose)
			printf("commencing curve %d of %u\n", curve, curves);

		sigma = 0;
		ubuild(&P, rho, &work, sigma, verbose);

		if (verbose)
		{
			printf("curve parameters:\n\tsigma = %u\n", work.sigma);
			printf("\tn = %" PRIu64 "\n", n);
			printf("\trho = %" PRIu64 "\n", rho);
			printf("\tx = %" PRIx64 "\n", P.X);
			printf("\tz = %" PRIx64 "\n", P.Z);
			printf("\tb = %" PRIx64 "\n", work.s);
		}
		uecm_stage1(rho, &work, &P);
		result = ucheck_factor(P.Z, n, &tmp1);
			
        if (verbose)
        {
            printf("after stage1: P = %" PRIx64 ", %" PRIx64 "\n", P.X, P.Z);
        }

		if (result == 1)
		{
			if (verbose)
				printf("\nfound factor %" PRIx64 " in stage 1 with sigma = %u\n",
					tmp1, work.sigma);

			*f = tmp1;
			break;
		}

		if (B2 > B1)
		{
			uecm_stage2(&P, rho, &work);

            if (verbose)
            {
                printf("after stage2: A = %" PRIx64 "\n", work.stg2acc);
            }

			if (verbose)
				printf("performed %d pair-multiplies for %" PRIu64 " primes in stage 2\n",
					work.paired, work.numprimes);
			if (verbose)
				printf("performed %u point-additions and %u point-doubles in stage 2\n",
					work.ptadds + work.stg1Add, work.stg1Doub);

			result = ucheck_factor(work.stg2acc, n, &tmp1);

			if (result == 1)
			{
				if (verbose)
					printf("\nfound factor %" PRIx64 " in stage 2 with sigma = %u\n",
						tmp1, work.sigma);

				*f = tmp1;
				break;
			}
		}
	}

	return;
}

void uecm_stage1(uint64_t rho, uecm_work *work, uecm_pt *P)
{
	int i;
	uint64_t q;
	uint64_t stg1 = (uint64_t)work->stg1_max;
	 
	// handle the only even case 
	q = 2;
	while (q < stg1)
	{
		work->diff1 = submod(P->X, P->Z, work->n);
		work->sum1 = addmod(P->X, P->Z, work->n);
		udup(rho, work, work->sum1, work->diff1, P);
		q *= 2;
	}

    if (stg1 == 47)
    {
        uprac(rho, work, P, 3, 0.618033988749894903);
        uprac(rho, work, P, 3, 0.618033988749894903);
        uprac(rho, work, P, 3, 0.618033988749894903);
        uprac(rho, work, P, 5, 0.618033988749894903);
        uprac(rho, work, P, 5, 0.618033988749894903);
        uprac(rho, work, P, 7, 0.618033988749894903);
        uprac(rho, work, P, 11, 0.580178728295464130);
        uprac(rho, work, P, 13, 0.618033988749894903);
        uprac(rho, work, P, 17, 0.618033988749894903);
        uprac(rho, work, P, 19, 0.618033988749894903);
        uprac(rho, work, P, 23, 0.522786351415446049);
        uprac(rho, work, P, 29, 0.548409048446403258);
        uprac(rho, work, P, 31, 0.618033988749894903);
        uprac(rho, work, P, 37, 0.580178728295464130);
        uprac(rho, work, P, 41, 0.548409048446403258);
        uprac(rho, work, P, 43, 0.618033988749894903);
        uprac(rho, work, P, 47, 0.548409048446403258);
        i = 15;
    }
    else if (stg1 == 59)
    {
        uprac(rho, work, P, 3, 0.61803398874989485);
        uprac(rho, work, P, 3, 0.61803398874989485);
        uprac(rho, work, P, 3, 0.61803398874989485);
        uprac(rho, work, P, 5, 0.618033988749894903);
        uprac(rho, work, P, 5, 0.618033988749894903);
        uprac(rho, work, P, 7, 0.618033988749894903);
        uprac(rho, work, P, 7, 0.618033988749894903);
        uprac(rho, work, P, 11, 0.580178728295464130);
        uprac(rho, work, P, 13, 0.618033988749894903);
        uprac(rho, work, P, 17, 0.618033988749894903);
        uprac(rho, work, P, 19, 0.618033988749894903);
        uprac(rho, work, P, 23, 0.522786351415446049);
        uprac(rho, work, P, 29, 0.548409048446403258);
        uprac(rho, work, P, 31, 0.618033988749894903);
        uprac(rho, work, P, 1961, 0.552936068843375);	// 37 * 53
        uprac(rho, work, P, 41, 0.548409048446403258);
        uprac(rho, work, P, 43, 0.618033988749894903);
        uprac(rho, work, P, 47, 0.548409048446403258);
        uprac(rho, work, P, 59, 0.548409048446403258);
        i = 17;
    }
    else if (stg1 == 70)
	{
		// call prac with best ratio found in deep search.
		// some composites are cheaper than their 
		// constituent primes.
        uprac70(rho, work, P);
		i = 19;
	}
	else if (stg1 >= 85)
	{
        uprac85(rho, work, P);
		if (stg1 < 100)
		{
			// paired into a composite for larger bounds
			uprac(rho, work, P, 61, 0.522786351415446049);
		}
		i = 23;

		if (stg1 >= 125)
		{
			uprac(rho, work, P, 5, 0.618033988749894903);
			uprac(rho, work, P, 11, 0.580178728295464130);
			uprac(rho, work, P, 61, 0.522786351415446049);
			uprac(rho, work, P, 89, 0.618033988749894903);
			uprac(rho, work, P, 97, 0.723606797749978936);
			uprac(rho, work, P, 101, 0.556250337855490828);
			uprac(rho, work, P, 107, 0.580178728295464130);
			uprac(rho, work, P, 109, 0.548409048446403258);
			uprac(rho, work, P, 113, 0.618033988749894903);

			i = 30;
		}
		if (stg1 < 130)
		{
			uprac(rho, work, P, 103, 0.632839806088706269);
		}

		if (stg1 >= 165)
		{
			uprac(rho, work, P, 7747, 0.552188778811121); // 61 x 127
			uprac(rho, work, P, 131, 0.618033988749894903);
			uprac(rho, work, P, 14111, 0.632839806088706);	// 103 x 137
			uprac(rho, work, P, 20989, 0.620181980807415);	// 139 x 151
			uprac(rho, work, P, 157, 0.640157392785047019);
			uprac(rho, work, P, 163, 0.551390822543526449);

			i = 38;
		}
		if (stg1 < 200)
		{
			uprac(rho, work, P, 149, 0.580178728295464130);
		}

		if (stg1 >= 205)
		{
			uprac(rho, work, P, 13, 0.618033988749894903);
			uprac(rho, work, P, 167, 0.580178728295464130);
			uprac(rho, work, P, 173, 0.612429949509495031);
			uprac(rho, work, P, 179, 0.618033988749894903);
			uprac(rho, work, P, 181, 0.551390822543526449);
			uprac(rho, work, P, 191, 0.618033988749894903);
			uprac(rho, work, P, 193, 0.618033988749894903);
			uprac(rho, work, P, 29353, 0.580178728295464);	// 149 x 197
			uprac(rho, work, P, 199, 0.551390822543526449);

			i = 46;
		}
	}

	work->last_pid = i;
	return;
}

// pre-paired sequences for various B1 and B2 = 25*B1
static const int numb1_70 = 186;
static uint8_t b1_70[186] = {
    53,49,47,43,41,37,23,19,13,11,1,7,17,29,31,0,59,47,43,41,37,31,29,19,13,7,1,11,23,0,59,53,43,41,37,
    31,23,17,11,7,1,19,29,49,0,53,49,47,43,31,23,19,11,7,1,13,37,59,0,59,53,43,37,31,29,23,17,13,11,1,47,
    0,59,49,41,31,23,17,11,7,1,19,37,47,0,59,49,47,43,41,31,17,13,11,7,37,0,53,49,43,37,23,19,13,7,1,29,
    31,41,59,0,59,49,47,41,23,19,17,13,7,1,43,53,0,59,49,43,37,29,17,13,7,1,19,47,53,0,59,53,49,47,43,31,
    29,23,11,17,0,47,43,41,37,31,23,19,17,11,1,13,29,53,0,59,47,41,37,31,23,19,11,7,17,29,0,53,47,43,41,
    17,13,11,1,23,31,37,49 };

static const int numb1_85 = 225;
static uint8_t b1_85[225] = {
    1,53,49,47,43,41,37,23,19,13,11,1,7,17,29,31,0,59,47,43,41,37,31,29,19,13,7,1,11,23,0,59,53,43,41,37,
    31,23,17,11,7,1,19,29,49,0,53,49,47,43,31,23,19,11,7,1,13,37,59,0,59,53,43,37,31,29,23,17,13,11,1,47,
    0,59,49,41,31,23,17,11,7,1,19,37,47,0,59,49,47,43,41,31,17,13,11,7,37,0,53,49,43,37,23,19,13,7,1,29,
    31,41,59,0,59,49,47,41,23,19,17,13,7,1,43,53,0,59,49,43,37,29,17,13,7,1,19,47,53,0,59,53,49,47,43,31,
    29,23,11,17,0,47,43,41,37,31,23,19,17,11,1,13,29,53,0,59,47,41,37,31,23,19,11,7,17,29,0,53,47,43,41,
    17,13,11,1,23,31,37,49,0,53,47,43,41,29,19,7,1,17,31,37,49,59,0,49,43,37,19,17,1,23,29,47,53,0,59,53,
    43,41,31,17,7,1,11,13,19,29 };

static const int numb1_125 = 319;
static uint8_t b1_125[319] = {
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
static uint8_t b1_165[425] = {
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
static uint8_t b1_205[511] = {
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


void uecm_stage2(uecm_pt *P, uint64_t rho, uecm_work *work)
{
	int b;
	int i, j, k, m;
	uecm_pt *Pa = &work->Pa;
	uecm_pt *Pb = work->Pb;
	uecm_pt *Pd;
	uint64_t acc = work->stg2acc;
    uint8_t *barray;
    uint32_t numb;

	work->paired = 0;
	work->numprimes = 0;
	work->ptadds = 0;
	work->stg1Add = 0;
	work->stg1Doub = 0;

	//stage 2 init
	//Q = P = result of stage 1
	//compute [d]Q for 0 < d <= D
	Pd = &Pb[map[D]];

	// [1]Q
	Pb[1].Z = P->Z;
	Pb[1].X = P->X;
	work->Pbprod[1] = mulredc(Pb[1].X, Pb[1].Z, work->n, rho);

	// [2]Q
	Pb[2].Z = P->Z;
	Pb[2].X = P->X;
	work->diff1 = submod(P->X, P->Z, work->n);
	work->sum1 = addmod(P->X, P->Z, work->n);
	udup(rho, work, work->sum1, work->diff1, &Pb[2]);
	work->Pbprod[2] = mulredc(Pb[2].X, Pb[2].Z, work->n, rho);

	/*
	D is small in tinyecm, so it is straightforward to just enumerate the needed
	points.  We can do it efficiently with two progressions mod 6.
	Pb[0] = scratch
	Pb[1] = [1]Q;
	Pb[2] = [2]Q;
	Pb[3] = [7]Q;	prog2
	Pb[4] = [11]Q;	prog1
	Pb[5] = [13]Q;	prog2
	Pb[6] = [17]Q;	prog1
	Pb[7] = [19]Q;	prog2
	Pb[8] = [23]Q;	prog1
	Pb[9] = [29]Q;	prog1
	Pb[10] = [30 == D]Q;
	Pb[11] = [31]Q;	prog2
	Pb[12] = [37]Q;	prog2
	Pb[13] = [41]Q;	prog1
	Pb[14] = [43]Q;	prog2
	Pb[15] = [47]Q;	prog1
	Pb[16] = [49]Q;	prog2
	Pb[17] = [53]Q;	prog1
	Pb[18] = [59]Q;	prog1

	two progressions with total of 17 adds to get 15 values of Pb.
	6 + 5(1) -> 11 + 6(5) -> 17 + 6(11) -> 23 + 6(17) -> 29 + 6(23) -> 35 + 6(29) -> 41 + 6(35) -> 47 + 6(41) -> 53 + 6(47) -> 59
	6 + 1(5) -> 7 + 6(1) -> 13 + 6(7) -> 19 + 6(13) -> 25 + 6(19) -> 31 + 6(25) -> 37 + 6(31) -> 43 + 6(37) -> 49

	we also need [2D]Q = [60]Q
	to get [60]Q we just need one more add:
	compute [60]Q from [31]Q + [29]Q([2]Q), all of which we
	have after the above progressions are computed.

	we also need [A]Q = [((B1 + D) / (2D) * 2D]Q
	which is equal to the following for various common B1:
	B1		[x]Q
	65		[120]Q
	85		[120]Q
	125		[180]Q
	165		[180]Q
	205		[240]Q

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

	// Calculate all Pb: the following is specialized for D=60
	// [2]Q + [1]Q([1]Q) = [3]Q
	uadd(rho, work, &Pb[1], &Pb[2], &Pb[1], &Pb[3]);		// <-- temporary

	// 2*[3]Q = [6]Q
	work->diff1 = submod(Pb[3].X, Pb[3].Z, work->n);
	work->sum1 = addmod(Pb[3].X, Pb[3].Z, work->n);
	udup(rho, work, work->sum1, work->diff1, &work->pt3);	// pt3 = [6]Q

	// [3]Q + [2]Q([1]Q) = [5]Q
	uadd(rho, work, &Pb[3], &Pb[2], &Pb[1], &work->pt1);	// <-- pt1 = [5]Q
	Pb[3].X = work->pt1.X;
	Pb[3].Z = work->pt1.Z;

	// [6]Q + [5]Q([1]Q) = [11]Q
	uadd(rho, work, &work->pt3, &work->pt1, &Pb[1], &Pb[4]);	// <-- [11]Q

	i = 3;
	k = 4;
	j = 5;
	while ((j + 12) < D)
	{
		// [j+6]Q + [6]Q([j]Q) = [j+12]Q
		uadd(rho, work, &work->pt3, &Pb[k], &Pb[i], &Pb[map[j + 12]]);
		i = k;
		k = map[j + 12];
		j += 6;
	}

	// [6]Q + [1]Q([5]Q) = [7]Q
	uadd(rho, work, &work->pt3, &Pb[1], &work->pt1, &Pb[3]);	// <-- [7]Q
	i = 1;
	k = 3;
	j = 1;
	while ((j + 12) < D)
	{
		// [j+6]Q + [6]Q([j]Q) = [j+12]Q
		uadd(rho, work, &work->pt3, &Pb[k], &Pb[i], &Pb[map[j + 12]]);
		i = k;
		k = map[j + 12];
		j += 6;
	}

	// Pd = [2w]Q
	// [31]Q + [29]Q([2]Q) = [60]Q
	uadd(rho, work, &Pb[9], &Pb[10], &Pb[2], Pd);	// <-- [60]Q
	work->ptadds++;

	// make all of the Pbprod's
	for (i = 3; i < 19; i++)
	{
		work->Pbprod[i] = mulredc(Pb[i].X, Pb[i].Z, work->n, rho);
	}

#if 1

    //initialize info needed for giant step
    // temporary - make [4]Q
    work->diff1 = submod(Pb[2].X, Pb[2].Z, work->n);
    work->sum1 = addmod(Pb[2].X, Pb[2].Z, work->n);
    udup(rho, work, work->sum1, work->diff1, &work->pt3);	// pt3 = [4]Q

    // Pd = [w]Q
    // [17]Q + [13]Q([4]Q) = [30]Q
    uadd(rho, work, &Pb[map[17]], &Pb[map[13]], &work->pt3, &work->Pad);	// <-- [30]Q

    // [60]Q + [30]Q([30]Q) = [90]Q
    uadd(rho, work, Pd, &work->Pad, &work->Pad, Pa);
    work->pt1.X = Pa->X;
    work->pt1.Z = Pa->Z;

    // [90]Q + [30]Q([60]Q) = [120]Q
    uadd(rho, work, Pa, &work->Pad, Pd, Pa);
    Pd->X = Pa->X;
    Pd->Z = Pa->Z;

    // [120]Q + [30]Q([90]Q) = [150]Q
    uadd(rho, work, Pa, &work->Pad, &work->pt1, Pa);

    // adjustment of Pa and Pad for larger B1.
    // Currently we have Pa=150, Pd=120, Pad=30
    if (work->stg1_max == 165)
    {
        // need Pa = 180, Pad = 60
        // [150]Q + [30]Q([120]Q) = [180]Q
        uadd(rho, work, Pa, &work->Pad, Pd, Pa);

        work->diff1 = submod(work->Pad.X, work->Pad.Z, work->n);
        work->sum1 = addmod(work->Pad.X, work->Pad.Z, work->n);
        udup(rho, work, work->sum1, work->diff1, &work->Pad);	// Pad = [60]Q
    }
    else if (work->stg1_max == 205)
    {
        // need Pa = 210, Pad = 90.
        // have pt1 = 90

        work->diff1 = submod(work->Pad.X, work->Pad.Z, work->n);
        work->sum1 = addmod(work->Pad.X, work->Pad.Z, work->n);
        udup(rho, work, work->sum1, work->diff1, &work->Pad);	// Pad = [60]Q

        // [150]Q + [60]Q([90]Q) = [210]Q
        uadd(rho, work, Pa, &work->Pad, &work->pt1, Pa);
        work->Pad.X = work->pt1.X;
        work->Pad.Z = work->pt1.Z;
    }

    //initialize accumulator and Paprod
    acc = u64div(1, work->n);
    work->Paprod = mulredc(Pa->X, Pa->Z, work->n, rho);

    if (work->stg1_max <= 70)
    {
        barray = b1_70;
        numb = numb1_70;
    }
    else if (work->stg1_max == 85)
    {
        barray = b1_85;
        numb = numb1_85;
    }
    else if (work->stg1_max == 125)
    {
        barray = b1_125;
        numb = numb1_125;
    }
    else if (work->stg1_max == 165)
    {
        barray = b1_165;
        numb = numb1_165;
    }
    else if (work->stg1_max == 205)
    {
        barray = b1_205;
        numb = numb1_205;
    }

    for (i = 0; i < numb; i++)
    {
        uint64_t tmp = acc;

        if (barray[i] == 0)
        {
            //giant step - use the addition formula for ECM
            work->pt1.X = Pa->X;
            work->pt1.Z = Pa->Z;

            //Pa + Pd
            uadd(rho, work, Pa, Pd, &work->Pad, Pa);

            //Pad holds the previous Pa
            work->Pad.X = work->pt1.X;
            work->Pad.Z = work->pt1.Z;

            //and Paprod
            work->Paprod = mulredc(Pa->X, Pa->Z, work->n, rho);

            i++;
        }

        //we accumulate XrZd - XdZr = (Xr - Xd) * (Zr + Zd) + XdZd - XrZr
        //in CP notation, Pa -> (Xr,Zr), Pb -> (Xd,Zd)

        b = barray[i];
        // accumulate the cross product  (zimmerman syntax).
        // page 342 in C&P
        work->tt1 = submod(Pa->X, Pb[map[b]].X, work->n);
        work->tt2 = addmod(Pa->Z, Pb[map[b]].Z, work->n);
        work->tt3 = mulredc(work->tt1, work->tt2, work->n, rho);
        work->tt1 = addmod(work->tt3, work->Pbprod[map[b]], work->n);
        work->tt2 = submod(work->tt1, work->Paprod, work->n);
        acc = mulredc(acc, work->tt2, work->n, rho);
        if (acc == 0)
        {
            acc = tmp;
            break;
        }
    }

    work->stg2acc = acc;
#else

    //first A value: first multiple of D greater than B1
    work->A = work->stg1_max / D + (work->stg1_max % D != 0);
    work->A *= D;

    //initialize info needed for giant step
    Pa->Z = P->Z;
    Pa->X = P->X;

    // currently we have Pd = 60.  Need to set Pa and Pad
    // appropriately for the stage 1 bound.
    if (work->stg1_max < 60)
    {
        // set Pa = [120]Q.  This means we may skip a few 
        // primes between stg1 and 60.
        work->diff1 = submod(Pd->X, Pd->Z, work->n);
        work->sum1 = addmod(Pd->X, Pd->Z, work->n);
        udup(rho, work, work->sum1, work->diff1, Pa);
        work->Pad.X = Pd->X;
        work->Pad.Z = Pd->Z;

        work->A = 120;
        work->last_pid = 17;
    }
    else
    {
        // set Pa = [120]Q and Pad = 120-60 = [60]Q
        work->diff1 = submod(Pd->X, Pd->Z, work->n);
        work->sum1 = addmod(Pd->X, Pd->Z, work->n);
        udup(rho, work, work->sum1, work->diff1, Pa);
        work->Pad.X = Pd->X;
        work->Pad.Z = Pd->Z;
    }
    
    if (work->stg1_max > 180)
    {
        // set Pa = 240... double again.
        // But first, Pad = 120+60[60]=180.
        uadd(rho, work, Pa, Pd, Pd, &work->Pad);
        work->diff1 = submod(Pa->X, Pa->Z, work->n);
        work->sum1 = addmod(Pa->X, Pa->Z, work->n);
        udup(rho, work, work->sum1, work->diff1, Pa);
    }
    else if (work->stg1_max > 120)
    {
        // set Pad = 120
        work->Pad.X = Pa->X;
        work->Pad.Z = Pa->Z;
        // now Pa = 120+60[60]=180.
        uadd(rho, work, Pa, Pd, Pd, Pa);
    }

    //and Paprod
    work->Paprod = mulredc(Pa->X, Pa->Z, work->n, rho);

    //initialize accumulator
    acc = u64div(1, work->n);

    //printf("starting accumulator is %" PRIx64 "\n", acc);
    //printf("Pad = %" PRIx64 ", %" PRIx64 "\n", work->Pad.X, work->Pad.Z);
    //printf("Pa = %" PRIx64 ", %" PRIx64 "\n", Pa->X, Pa->Z);
    //printf("Pd = %" PRIx64 ", %" PRIx64 "\n", Pd->X, Pd->Z);

    memset(marks, 0, D * sizeof(uint8_t));
    memset(nmarks, 0, D * sizeof(uint8_t));

    //begin stage 2
    i = work->last_pid;
    work->numprimes = 0;

    while (primes[i] > work->A)
    {
        //giant step - use the addition formula for ECM
        work->pt1.X = Pa->X;
        work->pt1.Z = Pa->Z;

        //Pa + Pd
        uadd(rho, work, Pa, Pd, &work->Pad, Pa);

        //Pad holds the previous Pa
        work->Pad.X = work->pt1.X;
        work->Pad.Z = work->pt1.Z;

        //and Paprod
        work->Paprod = mulredc(Pa->X, Pa->Z, work->n, rho);

        //next range
        work->A += D;
    }

    // compiler bug?  If I remove this conditional, the code no
    // longer works correctly, even if verbose is always 0.
    if (0)
    {
        printf("beginning work at prime %lu\n", primes[i]);
    }

    while ((primes[i] < work->stg2_max) && (i < 801))
    {
        b = work->A - primes[i];
        work->numprimes++;

        //printf("prime %u, index %d, acc is now %" PRIx64 "\n", primes[i], i, acc);

        if (!marks[b])
        {
            uint64 tmp = acc;
            //not marked, so doesn't have a match on the other side of the previous 'a'.
            //accumulate it, and mark the next range of 'a'.
            //we accumulate XrZd - XdZr = (Xr - Xd) * (Zr + Zd) + XdZd - XrZr
            //in CP notation, Pa -> (Xr,Zr), Pb -> (Xd,Zd)

            // accumulate the cross product  (zimmerman syntax).
            // page 342 in C&P
            work->tt1 = submod(Pa->X, Pb[map[b]].X, work->n);
            work->tt2 = addmod(Pa->Z, Pb[map[b]].Z, work->n);
            work->tt3 = mulredc(work->tt1, work->tt2, work->n, rho);
            work->tt1 = addmod(work->tt3, work->Pbprod[map[b]], work->n);
            work->tt2 = submod(work->tt1, work->Paprod, work->n);

            //printf("test: %" PRIx64 " * %" PRIx64 " mod %" PRIx64 " = %" PRIx64 "\n",
            //    acc, work->tt2, work->n, mulredc(acc, work->tt2, work->n, rho));

            acc = mulredc(acc, work->tt2, work->n, rho);

            if (acc == 0)
            {
                work->stg2acc = tmp;
                return;
            }

            // mark it so that if we hit a prime in a progression on the other side of D we
            // can skip it. 
            nmarks[D - b] = 1;
            work->paired++;
        }

        i++;

        while (primes[i] > work->A)
        {
            //set marks = nextmarks, then clear nextmarks
            memcpy(marks, nmarks, D * sizeof(uint8_t));
            memset(nmarks, 0, D * sizeof(uint8_t));

            //giant step - use the addition formula for ECM
            work->pt1.X = Pa->X;
            work->pt1.Z = Pa->Z;

            //Pa + Pd
            uadd(rho, work, Pa, Pd, &work->Pad, Pa);

            //Pad holds the previous Pa
            work->Pad.X = work->pt1.X;
            work->Pad.Z = work->pt1.Z;

            //and Paprod
            work->Paprod = mulredc(Pa->X, Pa->Z, work->n, rho);

            //next range
            work->A += D;
        }
    }

    work->last_pid = i;
    work->stg2acc = acc;

    return;
#endif
	
}

int ucheck_factor(uint64 Z, uint64 n, uint64 * f)
{
	int status;

	*f = bingcd64(n, Z);

	status = 0;
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



uint64 do_uecm(uint64 q64)
{
    int B1, curves, targetBits;
    uint64 f64;

    targetBits = spBits(q64) / 2;
    if (targetBits <= 24)
    {
        // multi-thread issue here...
        //f64 = LehmanFactor(q64, 0, 0, 0);
        B1 = 70;
        curves = 24;
        microecm(q64, &f64, B1, 25 * B1, curves, 0);
    }
    else if (targetBits <= 26)
    {
        B1 = 85;
        curves = 24;
        microecm(q64, &f64, B1, 25 * B1, curves, 0);
    }
    else if (targetBits <= 29)
    {
        B1 = 125;
        curves = 24;
        microecm(q64, &f64, B1, 25 * B1, curves, 0);
    }
    else if (targetBits <= 31)
    {
        B1 = 165;
        curves = 32;
        microecm(q64, &f64, B1, 25 * B1, curves, 0);
    }
    else if (targetBits <= 32)
    {
        B1 = 205;
        curves = 32;
        microecm(q64, &f64, B1, 25 * B1, curves, 0);
    }

    return f64;
}
