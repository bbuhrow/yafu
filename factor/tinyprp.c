// MIT License
// 
// Copyright (c) 2024 Ben Buhrow
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// MIT License
// 
// Copyright (c) 2024 Pierre
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <stdint.h>
#include <immintrin.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef GMP_CHECK
#include "gmp.h"
#endif

#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)

typedef __uint128_t uint128_t;
#include <sys/time.h>	//for gettimeofday using gcc

#elif defined(__GNUC__)

typedef __uint128_t uint128_t;
#include <sys/time.h>	//for gettimeofday using gcc

#elif defined(_MSC_VER)


#endif


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

/* ============= Begin routines borrowed or adapted from Perig ==============
See https://github.com/Boutoukoat/Euler-Sprp-fast-primality-checks
and https://www.mersenneforum.org/node/22163/page10
*/


static uint128_t my_random(void)
{
	// based on linear congruential generator, period = 2^128
	static uint128_t seed = ((uint128_t)0x123456789ull << 92) + ((uint128_t)0xabcdef << 36) + 0x987654321ull;
	seed = seed * 137 + 13;
	// shuffle
	uint128_t x = seed ^ (seed >> 17) ^ (seed << 13);
	return x;
}

// count trailing zeroes in binary representation 
static __inline
uint64_t my_ctz52(uint64_t n)
{
#if (INLINE_ASM && defined(__x86_64__))
#if defined(__BMI1__)
	uint64_t t;
	asm(" tzcntq %1, %0\n": "=r"(t) : "r"(n) : "flags");
	return t;
#else
	if (n)
		return __builtin_ctzll(n);
	return 52;
#endif
#else
#if defined(__GNUC__)
	if (n)
		return __builtin_ctzll(n);
	return 52;
#else
	if (n == 0)
		return 52;
	uint64_t r = 0;
	if ((n & 0xFFFFFFFFull) == 0)
		r += 32, n >>= 32;
	if ((n & 0xFFFFull) == 0)
		r += 16, n >>= 16;
	if ((n & 0xFFull) == 0)
		r += 8, n >>= 8;
	if ((n & 0xFull) == 0)
		r += 4, n >>= 4;
	if ((n & 0x3ull) == 0)
		r += 2, n >>= 2;
	if ((n & 0x1ull) == 0)
		r += 1;
	return r;
#endif
#endif
}

// count trailing zeroes in binary representation 
static __inline
uint64_t my_ctz104(uint64_t n_lo, uint64_t n_hi)
{
	if (n_lo) {
		return my_ctz52(n_lo);
	}
	return 52 + my_ctz52(n_hi);
}

// count leading zeroes in binary representation
static __inline
uint64_t my_clz52(uint64_t n)
{
#if (INLINE_ASM && defined(__x86_64__))
#ifdef __BMI1__
	uint64_t t;
	asm(" lzcntq %1, %0\n": "=r"(t) : "r"(n) : "flags");
	return t;
#else
	if (n)
		return __builtin_clzll(n);
	return 52;
#endif
#else
#if defined(__GNUC__)
	if (n)
		return __builtin_clzll(n);
	return 52;
#else
	if (n == 0)
		return 52;
	uint64_t r = 0;
	if ((n & (0xFFFFFFFFull << 32)) == 0)
		r += 32, n <<= 32;
	if ((n & (0xFFFFull << 48)) == 0)
		r += 16, n <<= 16;
	if ((n & (0xFFull << 56)) == 0)
		r += 8, n <<= 8;
	if ((n & (0xFull << 60)) == 0)
		r += 4, n <<= 4;
	if ((n & (0x3ull << 62)) == 0)
		r += 2, n <<= 2;
	if ((n & (0x1ull << 63)) == 0)
		r += 1;
	return r;
#endif
#endif
}

// count leading zeroes in binary representation 
static __inline
uint64_t my_clz104(uint64_t n_lo, uint64_t n_hi)
{
	if (n_hi) {
		return my_clz52(n_hi);
	}
	return 52 + my_clz52(n_lo);
}

// count leading zeroes in binary representation
static __inline
uint64_t my_clz64(uint64_t n)
{
#if (INLINE_ASM && defined(__x86_64__))
#ifdef __BMI1__
	uint64_t t;
	asm(" lzcntq %1, %0\n": "=r"(t) : "r"(n) : "flags");
	return t;
#else
	if (n)
		return __builtin_clzll(n);
	return 64;
#endif
#else
#if defined(__GNUC__)
	if (n)
		return __builtin_clzll(n);
	return 64;
#else
	if (n == 0)
		return 52;
	uint64_t r = 0;
	if ((n & (0xFFFFFFFFull << 32)) == 0)
		r += 32, n <<= 32;
	if ((n & (0xFFFFull << 48)) == 0)
		r += 16, n <<= 16;
	if ((n & (0xFFull << 56)) == 0)
		r += 8, n <<= 8;
	if ((n & (0xFull << 60)) == 0)
		r += 4, n <<= 4;
	if ((n & (0x3ull << 62)) == 0)
		r += 2, n <<= 2;
	if ((n & (0x1ull << 63)) == 0)
		r += 1;
	return r;
#endif
#endif
}


static __inline
uint64_t my_clz128(uint64_t n_lo, uint64_t n_hi)
{
	if (n_hi) {
		return my_clz64(n_hi);
	}
	return 64 + my_clz64(n_lo);
}

static inline uint64_t my_rdtsc(void)
{
#if defined(__x86_64__)
	// supported by GCC and Clang for x86 platform
	return _rdtsc();
#elif INLINE_ASM && defined(__aarch64__)
	// should be a 64 bits wallclock counter
	// document for old/recent architecture and/or BMC chipsets mention it
	// could be a 56 bit counter.
	uint64_t val;

	asm volatile ("mrs %0, cntvct_el0":"=r" (val));

	// I am not sure what the clock unit is, it depends on pre-scaler setup
	// A multiplication by 32 might be needed on my platform 
	return val * 32;	// aarch64 emulation on x86_64 ?
	return ((val / 3) * 25) << 4;	// maybe for ARM M1 ?
	return val;
#else
#error "todo : unsupported _rdtsc implementation\n"
	return 0;
#endif
}

/* ============= End routines borrowed or adapted from Perig ============== */

static __m512i lo52mask;

#ifdef IFMA

#define FORCE_INLINE __inline

FORCE_INLINE static __m512i mul52hi(__m512i b, __m512i c)
{
	return _mm512_madd52hi_epu64(_mm512_set1_epi64(0), c, b);
}

FORCE_INLINE static __m512i mul52lo(__m512i b, __m512i c)
{
	return _mm512_madd52lo_epu64(_mm512_set1_epi64(0), c, b);
}

FORCE_INLINE static void mul52lohi(__m512i b, __m512i c, __m512i* l, __m512i* h)
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

__inline static __m512i mul52lo(__m512i b, __m512i c)
{
	return _mm512_and_si512(_mm512_mullo_epi64(b, c), _mm512_set1_epi64(0x000fffffffffffffull));
}
__inline static __m512i mul52hi(__m512i b, __m512i c)
{
	__m512d prod1_ld = _mm512_cvtepu64_pd(b);
	__m512d prod2_ld = _mm512_cvtepu64_pd(c);
	prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	return _mm512_sub_epi64(castepu(prod1_ld), vbias1);
}
__inline static void mul52lohi(__m512i b, __m512i c, __m512i* l, __m512i* h)
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

void printvec(char* msg, __m512i v)
{
	uint64_t m[8];
	storeu64(m, v);
	int i;
	printf("%s: ", msg);
	for (i = 0; i < 8; i++)
		printf("%016lx ", m[i]);
	printf("\n");
	return;
}


double _difftime(struct timeval* start, struct timeval* end)
{
	double secs;
	double usecs;

	if (start->tv_sec == end->tv_sec) {
		secs = 0;
		usecs = end->tv_usec - start->tv_usec;
	}
	else {
		usecs = 1000000 - start->tv_usec;
		secs = end->tv_sec - (start->tv_sec + 1);
		usecs += end->tv_usec;
		if (usecs >= 1000000) {
			usecs -= 1000000;
			secs += 1;
		}
	}

	return secs + usecs / 1000000.;
}

#define carryprop(lo, hi, mask) \
	{ __m512i carry = _mm512_srli_epi64(lo, 52);	\
	hi = _mm512_add_epi64(hi, carry);		\
	lo = _mm512_and_epi64(mask, lo); }

__m512i __inline _mm512_addsetc_epi52(__m512i a, __m512i b, __mmask8* cout)
{
	__m512i t = _mm512_add_epi64(a, b);
	*cout = _mm512_cmpgt_epu64_mask(t, _mm512_set1_epi64(0xfffffffffffffULL));
	t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
	return t;
}
__m512i __inline _mm512_mask_addsetc_epi52(__m512i c, __mmask8 mask, __m512i a, __m512i b, __mmask8* cout)
{
	__m512i t = _mm512_add_epi64(a, b);
	*cout = _mm512_mask_cmpgt_epu64_mask(mask, t, _mm512_set1_epi64(0xfffffffffffffULL));
	t = _mm512_mask_and_epi64(c, mask, t, _mm512_set1_epi64(0xfffffffffffffULL));
	return t;
}
__m512i __inline _mm512_subsetc_epi52(__m512i a, __m512i b, __mmask8* cout)
{
	__m512i t = _mm512_sub_epi64(a, b);
	*cout = _mm512_cmpgt_epu64_mask(b, a);
	t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
	return t;
}
__m512i __inline _mm512_mask_subsetc_epi52(__m512i c, __mmask8 mask, __m512i a, __m512i b, __mmask8* cout)
{
	__m512i t = _mm512_sub_epi64(a, b);
	*cout = _mm512_mask_cmpgt_epu64_mask(mask, b, a);
	t = _mm512_mask_and_epi64(c, mask, t, _mm512_set1_epi64(0xfffffffffffffULL));
	return t;
}
__m512i __inline _mm512_adc_epi52(__m512i a, __mmask8 c, __m512i b, __mmask8* cout)
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


__inline static void mulredc52_mask_add_vec(__m512i* c0, __mmask8 addmsk, __m512i a0, __m512i b0, __m512i n0, __m512i vrho)
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
__inline static void mask_mulredc104_vec(__m512i* c1, __m512i* c0, __mmask8 mulmsk,
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
__inline static void sqrredc104_vec(__m512i* c1, __m512i* c0,
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
__inline static void mask_sqrredc104_vec(__m512i* c1, __m512i* c0, __mmask8 mulmsk,
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
__inline static void mask_sqrredc104_vec_pos(__m512i* c1, __m512i* c0, __mmask8 mulmsk,
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
__inline static void mask_sqrredc104_exact_vec(__m512i* c1, __m512i* c0, __mmask8 mulmsk,
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
__inline static void addmod104_x8(__m512i* c1, __m512i* c0, __m512i a1, __m512i a0,
	__m512i b1, __m512i b0, __m512i n1, __m512i n0)
{
	// add
	__mmask8 bmsk;
	//a0 = _mm512_addsetc_epi52(a0, b0, &bmsk);
	//a1 = _mm512_adc_epi52(a1, bmsk, b1, &bmsk);
	a0 = _mm512_add_epi64(a0, b0);
	a1 = _mm512_add_epi64(a1, b1);
	a1 = _mm512_add_epi64(a1, _mm512_srli_epi64(a0, 52));
	a0 = _mm512_and_epi64(a0, lo52mask);

	// compare
	//__mmask8 msk = bmsk | _mm512_cmpgt_epu64_mask(a1, n1);
	__mmask8 msk = _mm512_cmpgt_epu64_mask(a1, n1);
	msk |= (_mm512_cmpeq_epu64_mask(a1, n1) & _mm512_cmpge_epu64_mask(a0, n0));

	// conditionally subtract N
	*c0 = _mm512_mask_subsetc_epi52(a0, msk, a0, n0, &bmsk);
	*c1 = _mm512_mask_sbb_epi52(a1, msk, bmsk, n1, &bmsk);
	// *c0 = _mm512_mask_sub_epi64(a0, msk, a0, n0);
	// *c1 = _mm512_mask_sub_epi64(a1, msk, a1, n1);
	// *c1 = _mm512_mask_sub_epi64(*c1, msk, *c1, _mm512_srli_epi64(*c0, 63));
	return;
}
__inline static void mask_addmod104_x8(__m512i* c1, __m512i* c0, __mmask8 addmsk,
	__m512i a1, __m512i a0, __m512i b1, __m512i b0, __m512i n1, __m512i n0)
{
	// add
	__mmask8 bmsk;
	a0 = _mm512_mask_addsetc_epi52(a0, addmsk, a0, b0, &bmsk);
	a1 = _mm512_mask_adc_epi52(a1, addmsk, bmsk, b1, &bmsk);

	// compare
	__mmask8 msk = bmsk | _mm512_cmpgt_epu64_mask(a1, n1);
	msk |= (_mm512_cmpeq_epu64_mask(a1, n1) & _mm512_cmpge_epu64_mask(a0, n0));

	// conditionally subtract N
	*c0 = _mm512_mask_subsetc_epi52(a0, addmsk & msk, a0, n0, &bmsk);
	*c1 = _mm512_mask_sbb_epi52(a1, addmsk & msk, bmsk, n1, &bmsk);
	return;
}
__inline static void mask_dblmod104_x8(__m512i* c1, __m512i* c0, __mmask8 addmsk,
	__m512i a1, __m512i a0, __m512i n1, __m512i n0)
{
	// add
	__mmask8 bmsk;
	//a0 = _mm512_mask_addsetc_epi52(a0, addmsk, a0, b0, &bmsk);
	//a1 = _mm512_mask_adc_epi52(a1, addmsk, bmsk, b1, &bmsk);

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
	*c0 = _mm512_mask_subsetc_epi52(a0, addmsk & msk, a0, n0, &bmsk);
	*c1 = _mm512_mask_sbb_epi52(a1, addmsk & msk, bmsk, n1, &bmsk);
	return;
}
__inline static void mask_redsub104_x8(__m512i* c1, __m512i* c0, __mmask8 addmsk,
	__m512i a1, __m512i a0, __m512i n1, __m512i n0)
{
	__mmask8 bmsk;

	// compare
	__mmask8 msk = _mm512_cmpgt_epu64_mask(a1, n1);
	msk |= (_mm512_cmpeq_epu64_mask(a1, n1) & _mm512_cmpge_epu64_mask(a0, n0));

	// conditionally subtract N
	*c0 = _mm512_mask_subsetc_epi52(a0, addmsk & msk, a0, n0, &bmsk);
	*c1 = _mm512_mask_sbb_epi52(a1, addmsk & msk, bmsk, n1, &bmsk);
	return;
}
__inline static void redsub104_x8(__m512i* c1, __m512i* c0,
	__m512i a1, __m512i a0, __m512i n1, __m512i n0)
{
	__mmask8 bmsk;

	// compare
	__mmask8 msk = _mm512_cmpgt_epu64_mask(a1, n1);
	msk |= (_mm512_cmpeq_epu64_mask(a1, n1) & _mm512_cmpge_epu64_mask(a0, n0));

	// conditionally subtract N
	*c0 = _mm512_mask_subsetc_epi52(a0, msk, a0, n0, &bmsk);
	*c1 = _mm512_mask_sbb_epi52(a1, msk, bmsk, n1, &bmsk);

	return;
}
__inline static void submod104_x8(__m512i* c1, __m512i* c0, __m512i a1, __m512i a0,
	__m512i b1, __m512i b0, __m512i n1, __m512i n0)
{
	// compare
	__mmask8 msk = _mm512_cmplt_epu64_mask(a1, b1);
	msk |= _mm512_cmpeq_epu64_mask(a1, b1) & _mm512_cmplt_epu64_mask(a0, b0);

	// subtract
	__mmask8 bmsk;
	a0 = _mm512_subsetc_epi52(a0, b0, &bmsk);
	a1 = _mm512_sbb_epi52(a1, bmsk, b1, &bmsk);

	// conditionally add N
	*c0 = _mm512_mask_addsetc_epi52(a0, msk, a0, n0, &bmsk);
	*c1 = _mm512_mask_adc_epi52(a1, msk, bmsk, n1, &bmsk);
	return;
}

static __inline uint64_t multiplicative_inverse(uint64_t a)
{
	// compute the 64-bit inverse of a mod 2^64
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
	return x4;
}

static __inline __m512i multiplicative_inverse104_x8(uint64_t* a)
{
	//    assert(a%2 == 1);  // the inverse (mod 2<<64) only exists for odd values
	__m512i x0, x1, x2, x3, x4, x5, y, n, i0, i1;
	__m512i three = _mm512_set1_epi64(3), two = _mm512_set1_epi64(2), one = _mm512_set1_epi64(1);
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);

	n = loadu64(a);
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

#define USE_INLINE_ASM_X86

#if defined(USE_INLINE_ASM_X86) && !defined(_MSC_VER)

FORCE_INLINE uint64_t submod64(uint64_t a, uint64_t b, uint64_t n)
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

FORCE_INLINE uint64_t addmod64(uint64_t x, uint64_t y, uint64_t n)
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

FORCE_INLINE uint64_t submod64(uint64_t a, uint64_t b, uint64_t n)
{
	uint64_t r0;
	if (_subborrow_u64(0, a, b, &r0))
		r0 += n;
	return r0;
}

FORCE_INLINE uint64_t addmod64(uint64_t x, uint64_t y, uint64_t n)
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
FORCE_INLINE static uint64_t mulredc_alt(uint64_t x, uint64_t y, uint64_t N, uint64_t invN)
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
#if defined(USE_INLINE_ASM_X86) && !defined(_MSC_VER)
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
	return result;
}

/* --- end Hurchalla functions --- */


// full strength mul/sqr redc
FORCE_INLINE static uint64_t mulredc64(uint64_t x, uint64_t y, uint64_t n, uint64_t inv)
{
#if defined(__unix__) && defined(__x86_64__)
	// On Intel Skylake: 9 cycles latency, 7 fused uops.
	__uint128_t prod = (__uint128_t)x * y;
	uint64_t Thi = (uint64_t)(prod >> 64);
	uint64_t rrax = (uint64_t)(prod);
	__asm__(
		"imulq %[inv], %%rax \n\t"        /* m = T_lo * invN */
		"mulq %[n] \n\t"                  /* mN = m * N */
		"leaq (%[Thi], %[n]), %%rax \n\t" /* rax = T_hi + N */
		"subq %%rdx, %%rax \n\t"          /* rax = rax - mN_hi */
		"subq %%rdx, %[Thi] \n\t"         /* t_hi = T_hi - mN_hi */
		"cmovbq %%rax, %[Thi] \n\t"       /* t_hi = (T_hi<mN_hi) ? rax : t_hi */
		: [Thi] "+&bcSD"(Thi), "+&a"(rrax)
		: [n] "r"(n), [inv]"r"(inv)
		: "rdx", "cc");
	uint64_t result = Thi;
	return result;

#else
	return mulredc_alt(x, y, n, inv);
#endif
}
FORCE_INLINE static uint64_t sqrredc64(uint64_t x, uint64_t n, uint64_t inv)
{
#if defined(__unix__) && defined(__x86_64__)
	// On Intel Skylake: 9 cycles latency, 7 fused uops.
	__uint128_t prod = (__uint128_t)x * x;
	uint64_t Thi = (uint64_t)(prod >> 64);
	uint64_t rrax = (uint64_t)(prod);
	__asm__(
		"imulq %[inv], %%rax \n\t"        /* m = T_lo * invN */
		"mulq %[n] \n\t"                  /* mN = m * N */
		"leaq (%[Thi], %[n]), %%rax \n\t" /* rax = T_hi + N */
		"subq %%rdx, %%rax \n\t"          /* rax = rax - mN_hi */
		"subq %%rdx, %[Thi] \n\t"         /* t_hi = T_hi - mN_hi */
		"cmovbq %%rax, %[Thi] \n\t"       /* t_hi = (T_hi<mN_hi) ? rax : t_hi */
		: [Thi] "+&bcSD"(Thi), "+&a"(rrax)
		: [n] "r"(n), [inv]"r"(inv)
		: "rdx", "cc");
	uint64_t result = Thi;
	return result;

#else
	return mulredc_alt(x, x, n, inv);
#endif
}


#define USE_PERIG_128BIT

/********************* begin Perig's 128-bit code  **********************/

#ifdef USE_PERIG_128BIT
typedef __uint128_t uint128_t;

static void ciosModMul128(uint64_t* res_lo, uint64_t* res_hi, uint64_t b_lo, uint64_t b_hi, uint64_t mod_lo, uint64_t mod_hi,
	uint64_t mmagic)
{
	uint64_t a_lo = *res_lo, a_hi = *res_hi;
	uint128_t cs, cc;
	uint64_t t0, t1, t2, t3, m, ignore;

	cc = (uint128_t)a_lo * b_lo;	// #1
	t0 = (uint64_t)cc;
	cc = cc >> 64;
	cc += (uint128_t)a_lo * b_hi;	// #2
	t1 = (uint64_t)cc;
	cc = cc >> 64;
	t2 = (uint64_t)cc;
#if PARANOID
	assert(cc >> 64 == 0);
#endif

	m = t0 * mmagic;	// #3
	cs = (uint128_t)m * mod_lo;	// #4
	cs += t0;
	cs = cs >> 64;
	cs += (uint128_t)m * mod_hi;	// #5
	cs += t1;
	t0 = (uint64_t)cs;
	cs = cs >> 64;
	cs += t2;
	t1 = (uint64_t)cs;
	cs = cs >> 64;
	t2 = (uint64_t)cs;
#if PARANOID
	assert(cs >> 64 == 0);
#endif

	cc = (uint128_t)a_hi * b_lo;	// #6
	cc += t0;
	t0 = (uint64_t)cc;
	cc = cc >> 64;
	cc += (uint128_t)a_hi * b_hi;	// #7
	cc += t1;
	t1 = (uint64_t)cc;
	cc = cc >> 64;
	cc += t2;
	t2 = (uint64_t)cc;
	cc = cc >> 64;
	t3 = (uint64_t)cc;
#if PARANOID
	assert(cc >> 64 == 0);
#endif

	m = t0 * mmagic;	// #8
	cs = (uint128_t)m * mod_lo;	// #9
	cs += t0;
	cs = cs >> 64;
	cs += (uint128_t)m * mod_hi;	// #10
	cs += t1;
	t0 = (uint64_t)cs;
	cs = cs >> 64;
	cs += t2;
	t1 = (uint64_t)cs;
	cs = cs >> 64;
	cs += t3;
	t2 = (uint64_t)cs;
	if (t2) {
		unsigned char carry = _subborrow_u64(0, t0, mod_lo, &t0);
		_subborrow_u64(carry, t1, mod_hi, &t1);
		//ciosSubtract128(&t0, &t1, t2, mod_lo, mod_hi);
	}

	*res_lo = t0;
	*res_hi = t1;
}

/********************* end of Perig's 128-bit code **********************/

// modSqr version I wrote based on the modMul
static void ciosModSqr128(uint64_t* res_lo, uint64_t* res_hi, uint64_t b_lo, uint64_t b_hi, uint64_t mod_lo, uint64_t mod_hi,
	uint64_t mmagic)
{
	uint128_t cs, cc, b_lohi;
	uint64_t t0, t1, t2, t3, m, ignore;

	cc = (uint128_t)b_lo * b_lo;	// #1
	t0 = (uint64_t)cc;
	cc = cc >> 64;
	b_lohi = (uint128_t)b_lo * b_hi;	// #2
	cc += b_lohi;
	t1 = (uint64_t)cc;
	cc = cc >> 64;
	t2 = (uint64_t)cc;
#if PARANOID
	assert(cc >> 64 == 0);
#endif

	m = t0 * mmagic;	// #3
	cs = (uint128_t)m * mod_lo;	// #4
	cs += t0;
	cs = cs >> 64;
	cs += (uint128_t)m * mod_hi;	// #5
	cs += t1;
	t0 = (uint64_t)cs;
	cs = cs >> 64;
	cs += t2;
	t1 = (uint64_t)cs;
	cs = cs >> 64;
	t2 = (uint64_t)cs;
#if PARANOID
	assert(cs >> 64 == 0);
#endif

	cc = b_lohi + t0;
	t0 = (uint64_t)cc;
	cc = cc >> 64;
	cc += (uint128_t)b_hi * b_hi;	// #6
	cc += t1;
	t1 = (uint64_t)cc;
	cc = cc >> 64;
	cc += t2;
	t2 = (uint64_t)cc;
	cc = cc >> 64;
	t3 = (uint64_t)cc;
#if PARANOID
	assert(cc >> 64 == 0);
#endif

	m = t0 * mmagic;	// #8
	cs = (uint128_t)m * mod_lo;	// #9
	cs += t0;
	cs = cs >> 64;
	cs += (uint128_t)m * mod_hi;	// #10
	cs += t1;
	t0 = (uint64_t)cs;
	cs = cs >> 64;
	cs += t2;
	t1 = (uint64_t)cs;
	cs = cs >> 64;
	cs += t3;
	t2 = (uint64_t)cs;
	if (t2) {
		unsigned char carry = _subborrow_u64(0, t0, mod_lo, &t0);
		_subborrow_u64(carry, t1, mod_hi, &t1);
		//ciosSubtract128(&t0, &t1, t2, mod_lo, mod_hi);
	}

	*res_lo = t0;
	*res_hi = t1;
}

// modSqr version I wrote based on the modMul
static void ciosModSqr128_pos_debug(uint64_t* res_lo, uint64_t* res_hi,
	uint64_t b_lo, uint64_t b_hi, uint64_t mod_lo, uint64_t mod_hi, uint64_t mmagic, int verbose)
{
	uint128_t cs, cc, b_lohi, tt;
	uint64_t t0, t1, t2, t3, m, ignore;

	cc = (uint128_t)b_lo * b_lo;	// #1
	t0 = (uint64_t)cc;
	cc = cc >> 64;
	b_lohi = (uint128_t)b_lo * b_hi;	// #2
	cc += b_lohi;
	t1 = (uint64_t)cc;
	cc = cc >> 64;
	t2 = (uint64_t)cc;
#if PARANOID
	assert(cc >> 64 == 0);
#endif
	tt = ((uint128_t)t2 << 64) | t1;

	if (verbose)
	{
		printf("x  = 0x%016lx%016lx\n", b_hi, b_lo);
		printf("n  = 0x%016lx%016lx\n", mod_hi, mod_lo);
		printf("n' = 0x%016lx\n", mmagic);
		printf("t  = 0x%016lx%016lx%016lx\n", t2, t1, t0);
	}

	// in the pos variant, we subtract mN instead of add mN.
	// The advantage is that the lo part of mN never generates
	// a borrow when subtracted from T, so we can ignore
	// a step of the carry propagation.
	m = t0 * mmagic;	// #3
	if (verbose)
	{
		printf("m  = 0x%016lx\n", m);
	}
	cs = (uint128_t)m * mod_lo;	// #4
	if (verbose)
	{
		printf("mNlo  = 0x%016lx%016lx\n", (uint64_t)(cs >> 64), (uint64_t)cs);
	}
	//cs += t0;		// equals 0, and never generates a borrow
	cs = cs >> 64;
	tt -= cs;
	if (verbose)
	{
		printf("t - mn0hi = 0x%016lx%016lx\n", (uint64_t)(tt >> 64), (uint64_t)tt);
	}
	cs = (uint128_t)m * mod_hi;	// #5
	tt -= cs;
	int sign;
	if (verbose)
	{
		printf("t - mn1   = 0x%016lx%016lx\n", (uint64_t)(tt >> 64), (uint64_t)tt);
	}

#if PARANOID
	assert(cs >> 64 == 0);
#endif


	tt += b_lohi;
	cc = (uint128_t)b_hi * b_hi;	// #6

#if PARANOID
	assert(cc >> 64 == 0);
#endif

	sign = tt >> 127;
	if (verbose)
	{
		printf("t     = 0x%016lx%016lx (sign %d)\n", (uint64_t)(tt >> 64), (uint64_t)tt, sign);
	}

	// in the pos variant, we subtract mN instead of add mN.
	// The advantage is that the lo part of mN never generates
	// a borrow when subtracted from T, so we can ignore
	// a step of the carry propagation.
	m = (uint64_t)tt * mmagic;	// #3
	if (verbose)
	{
		printf("m  = 0x%016lx\n", m);
	}
	cs = (uint128_t)m * mod_lo;	// #4
	if (verbose)
	{
		printf("mNlo  = 0x%016lx%016lx\n", (uint64_t)(cs >> 64), (uint64_t)cs);
	}
	//cs += t0;		// equals 0, and never generates a borrow
	cs = cs >> 64;
	tt = (tt >> 64);
	if (sign)
		tt |= ((uint128_t)0xffffffffffffffffull << 64);
	tt -= cs;
	if (verbose)
	{
		printf("t - mn0hi = 0x%016lx%016lx\n", (uint64_t)(tt >> 64), (uint64_t)tt);
	}
	cs = (uint128_t)m * mod_hi;	// #5
	tt -= cs;
	if (tt < cs) tt -= ((uint128_t)1 << 64);
	if (verbose)
	{
		printf("t - mn1   = 0x%016lx%016lx\n", (uint64_t)(tt >> 64), (uint64_t)tt);
	}

	tt += cc;
	if (verbose)
	{
		printf("t + x1^2  = 0x%016lx%016lx\n", (uint64_t)(tt >> 64), (uint64_t)tt);
	}

	t1 = (uint64_t)(tt >> 64);
	t0 = (uint64_t)tt;
	if (t1 >> 63) {
		unsigned char carry = _addcarry_u64(0, t0, mod_lo, &t0);
		_addcarry_u64(carry, t1, mod_hi, &t1);
		//ciosSubtract128(&t0, &t1, t2, mod_lo, mod_hi);
	}

	if (verbose)
	{
		printf("t = 0x%016lx%016lx\n", t1, t0);
		exit(1);
	}
	*res_lo = t0;
	*res_hi = t1;
}

static void ciosModSqr128_pos(uint64_t* res_lo, uint64_t* res_hi,
	uint64_t b_lo, uint64_t b_hi, uint64_t mod_lo, uint64_t mod_hi, uint64_t mmagic)
{
	uint128_t cs, cc, b_lohi, tt;
	uint64_t t0, t1, m;
	int sign;

	// in the pos variant, we subtract mN instead of add mN.
	// The advantage is that the lo part of mN never generates
	// a borrow when subtracted from T, so we can ignore
	// the low word of T
	tt = (uint128_t)b_lo * b_lo;	// #1
	t0 = (uint64_t)tt;
	tt = tt >> 64;
	b_lohi = (uint128_t)b_lo * b_hi;	// #2
	tt += b_lohi;

	// and we can skip the first part of the carry propagation
	// when subtracting m*N
	m = t0 * mmagic;	// #3
	cs = (uint128_t)m * mod_lo;	// #4
	cs = cs >> 64;
	tt -= cs;
	cs = (uint128_t)m * mod_hi;	// #5
	tt -= cs;
	tt += b_lohi;
	cc = (uint128_t)b_hi * b_hi;	// #6

	// did the subtract generate a borrow?  will need 
	// to generate a signed shift if so.
	sign = tt >> 127;

	// round 2.  again we ignore the low word and first part
	// of carry propagation.
	m = (uint64_t)tt * mmagic;	// #3
	cs = (uint128_t)m * mod_lo;	// #4
	cs = cs >> 64;
	tt = tt >> 64;

	// that's a signed right shift, so pull in the f's as needed
	tt = sign ? ((uint128_t)0xffffffffffffffffull << 64) | tt : tt;
	tt -= cs;

	//unsigned char b;
	//b = _subborrow_u64(0, (uint64_t)tt, (uint64_t)cs, &t0);
	//_subborrow_u64(0, sign ? 0xffffffffffffffffull : 0, 0, &t1);
	//tt = ((uint128_t)t1 << 64) | t0;

	cs = (uint128_t)m * mod_hi;	// #5
	tt -= cs;
	tt += cc;

	t1 = (uint64_t)(tt >> 64);
	t0 = (uint64_t)tt;

	// comparison is also easier in the positive variant since we
	// only have to check for a borrow (not >= n);
	if (t1 >> 63) {
		unsigned char carry = _addcarry_u64(0, t0, mod_lo, &t0);
		_addcarry_u64(carry, t1, mod_hi, &t1);
	}

	*res_lo = t0;
	*res_hi = t1;
}
#endif

#define POS_VARIANT
#define LO(x) ((uint64_t)((x) & 0xffffffffffffffffull))
#define HI(x) ((uint64_t)((x) >> 64))

__inline void chkmod128(uint64_t* a, uint64_t* n)
{

#if defined(__x86_64__) && defined(__unix__)

	// not really faster than the branching version below
	// but probably should be checked on more cpus.
	uint64_t t0, t1;
	t0 = a[0];
	t1 = a[1];
	__asm__ volatile (
		"xorq %%r8, %%r8 \n\t"
		"xorq %%r9, %%r9 \n\t"
		"subq %2, %%r8 \n\t"		/* t = 0 - n */
		"sbbq %3, %%r9 \n\t"
		"addq %0, %%r8 \n\t"		/* t += x */
		"adcq %1, %%r9 \n\t"
		"cmovc %%r8, %0 \n\t"
		"cmovc %%r9, %1 \n\t"
		: "+&r"(t0), "+&r"(t1)
		: "r"(n[0]), "r"(n[1])
		: "r8", "r9", "cc", "memory");

	a[0] = t0;
	a[1] = t1;
#else

	if ((a[1] > n[1]) || ((a[1] == n[1]) && (a[0] >= n[0])))
	{
		uint8_t c1 = _subborrow_u64(0, a[0], n[0], &a[0]);
		_subborrow_u64(c1, a[1], n[1], &a[1]);
	}

#endif
	return;
}

__inline void dblmod128(uint64_t* a, uint64_t* n)
{

#if defined(__x86_64__) && defined(__unix__)
	// requires GCC_ASM64 syntax
	uint64_t t1, t0;

	// do the 2-word variant of this:
	//uint64_t r;
	//uint64_t tmp = x - n;
	//uint8_t c = _addcarry_u64(0, tmp, y, &r);
	//return (c) ? r : x + y;
	t1 = a[1];
	t0 = a[0];

	__asm__ volatile (
		"movq %%rax, %%r8 \n\t"
		"movq %%rdx, %%r9 \n\t"
		"subq %4, %%r8 \n\t"		/* t = x - n */
		"sbbq %5, %%r9 \n\t"
		"addq %%rax, %0 \n\t"		/* x += x */
		"adcq %%rdx, %1 \n\t"
		"addq %%rax, %%r8 \n\t"		/* t = t + x */
		"adcq %%rdx, %%r9 \n\t"
		"cmovc %%r8, %0 \n\t"
		"cmovc %%r9, %1 \n\t"
		: "+&r"(a[0]), "+&r"(a[1])
		: "a"(t0), "d"(t1), "r"(n[0]), "r"(n[1])
		: "r8", "r9", "cc", "memory");

#else

	uint128_t x = ((uint128_t)a[1] << 64) | a[0];
	uint128_t n128 = ((uint128_t)n[1] << 64) | n[0];
	uint128_t r = (x >= n128 - x) ? x - (n128 - x) : x + x;
	a[1] = (uint64_t)(r >> 64);
	a[0] = (uint64_t)r;

#endif
	return;
}

__inline void submod128(uint64_t* a, uint64_t* b, uint64_t* w, uint64_t* n)
{
#if 1 //defined(__x86_64__) && defined(__unix__)
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

static void mulmod128(uint64_t* u, uint64_t* v, uint64_t* w, uint64_t* n, uint64_t rho)
{
	w[0] = v[0];
	w[1] = v[1];
	ciosModMul128(&w[0], &w[1], u[0], u[1], n[0], n[1], rho);
	return;
}

static void sqrmod128(uint64_t* u, uint64_t* w, uint64_t* n, uint64_t rho)
{
#ifdef POS_VARIANT

	ciosModSqr128_pos(&w[0], &w[1], u[0], u[1], n[0], n[1], rho);

	//ciosModSqr128(&res0, &res1, test0, test1, n[0], n[1], 0ULL - rho);
	//
	//if ((res0 != w[0]) || (res1 != w[1]))
	//{
	//	printf("problem!\n");
	//	printf("neg variant: 0x%016lx%016lx\n", w[1], w[0]);
	//	printf("pos variant: 0x%016lx%016lx\n", res1, res0);
	//	ciosModSqr128_pos(&w[0], &w[1], test0, test1, n[0], n[1], rho, 1);
	//}

#else
	ciosModSqr128(&w[0], &w[1], u[0], u[1], n[0], n[1], rho);
#endif
	return;
}

int fermat_prp_64x1(uint64_t n)
{
	// assumes has no small factors.
	// do a base-2 fermat prp test using LR binexp.

	uint64_t rho = multiplicative_inverse(n);
	uint64_t unityval = ((uint64_t)0 - n) % n;  // unityval == R  (mod n)
	uint64_t m = 1ULL << (62 - __lzcnt64(n));   // set a mask at the leading bit - 2
	uint64_t r = unityval;
	uint64_t e = n - 1;

	// we know the first bit is set and the first squaring is of unity,
	// so we can do the first iteration manually with no squaring.
	r = addmod64(r, r, n);

	while (m > 0)
	{
		r = sqrredc64(r, n, rho);
		if (e & m) r = addmod64(r, r, n);
		m >>= 1;
	}
	return (r == unityval);
}

int fermat_prp_128x1(uint64_t* n)
{
	// assumes has no small factors.  
	// assumes n has two 64-bit words, n0 and n1.
	// do a base-2 fermat prp test using LR binexp.

	uint64_t rho = multiplicative_inverse(n[0]);
	uint64_t unityval[2]; // = ((uint64_t)0 - n) % n;  // unityval == R  (mod n)
	uint64_t m;
	uint64_t r[2];
	uint64_t e[2];
	uint64_t one[2] = { 1ULL, 0ULL };
	uint64_t zero[2] = { 0ULL, 0ULL };

#ifdef POS_VARIANT

#else
	rho = 0ULL - rho;
#endif

	uint128_t n128 = ((uint128_t)n[1] << 64) | n[0];
	uint128_t unity128 = (uint128_t)0 - n128;
	unity128 = unity128 % n128;
	unityval[1] = (uint64_t)(unity128 >> 64);
	unityval[0] = (uint64_t)unity128;

	e[1] = n[1];
	e[0] = n[0] - 1;	// n odd: won't overflow

	r[1] = unityval[1];
	r[0] = unityval[0];

	int lzcnt = my_clz64(e[1]);

#ifndef POS_VARIANT
	int protect = (lzcnt < 3) ? 1 : 0;
#endif

	// we know the first bit is set and the first squaring is of unity,
	// so we can do the first iteration manually with no squaring.
	dblmod128(r, n);

#ifndef POS_VARIANT
	if (protect)
	{
		m = (e[1] <= 1) ? 0 : 1ULL << (62 - lzcnt);   // set a mask at the leading bit - 2
		while (m > 0)
		{
			sqrmod128(r, r, n, rho);
			chkmod128(r, n);

			if (e[1] & m) dblmod128(r, n);
			m >>= 1;
		}

		m = (e[1] >= 1) ? 1ULL << 63 : 1ULL << (62 - my_clz64(e[0]));   // set a mask at the leading bit - 2
		while (m > 0)
		{
			sqrmod128(r, r, n, rho);
			chkmod128(r, n);

			if (e[0] & m) dblmod128(r, n);
			m >>= 1;
		}
	}
	else
#endif
	{
		m = (e[1] <= 1) ? 0 : 1ULL << (62 - lzcnt);   // set a mask at the leading bit - 2
		while (m > 0)
		{
			sqrmod128(r, r, n, rho);

			if (e[1] & m) dblmod128(r, n);
			m >>= 1;
		}

		m = (e[1] >= 1) ? 1ULL << 63 : 1ULL << (62 - my_clz64(e[0]));   // set a mask at the leading bit - 2
		while (m > 0)
		{
			sqrmod128(r, r, n, rho);

			if (e[0] & m) dblmod128(r, n);
			m >>= 1;
		}
	}

	chkmod128(r, n);
	return ((r[0] == unityval[0]) && (r[1] == unityval[1]));
}

__m512i rem_epu64_x8(__m512i n, __m512i d)
{
	// DANGER: I haven't proven this works for every possible input.
	__m512d d1pd = _mm512_cvtepu64_pd(d);
	__m512d n1pd = _mm512_cvtepu64_pd(n);
	__m512i q, q2, r;

	//n1pd = _mm512_div_pd(n1pd, d1pd);
	//q = _mm512_cvt_roundpd_epu64(n1pd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

	n1pd = _mm512_div_round_pd(n1pd, d1pd, ROUNDING_MODE);
	q = _mm512_cvttpd_epu64(n1pd);

	__m512i qd = _mm512_mullox_epi64(q, d);
	r = _mm512_sub_epi64(n, qd);

	// fix q too big by a little, with special check for
	// numerators close to 2^64 and denominators close to 1
	// DANGER: the special check is unused for 64-bits, only for 32-bits.
	// This routine is only used in modmul32 and input numerators
	// shouldn't get that large in normal cases.  The factor base
	// would need to be close to 2^32...
	__mmask8 err = _mm512_cmpgt_epu64_mask(r, n); // |
		//(_mm512_cmpgt_epu64_mask(r, d) & _mm512_cmplt_epu64_mask(
		//	_mm512_sub_epi64(_mm512_set1_epi64(0), r), _mm512_set1_epi64(1024)));
	if (err)
	{
		n1pd = _mm512_cvtepu64_pd(_mm512_sub_epi64(_mm512_set1_epi64(0), r));

		//n1pd = _mm512_div_pd(n1pd, d1pd);
		//q2 = _mm512_add_epi64(_mm512_set1_epi64(1), _mm512_cvt_roundpd_epu64(n1pd,
		//	(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));

		n1pd = _mm512_div_round_pd(n1pd, d1pd, ROUNDING_MODE);
		q2 = _mm512_add_epi64(_mm512_set1_epi64(1), _mm512_cvttpd_epu64(n1pd));

		q = _mm512_mask_sub_epi64(q, err, q, q2);
		r = _mm512_mask_add_epi64(r, err, r, _mm512_mullox_epi64(q2, d));
	}

	// fix q too small by a little bit
	err = _mm512_cmpge_epu64_mask(r, d);
	if (err)
	{
		n1pd = _mm512_cvtepu64_pd(r);

		//n1pd = _mm512_div_pd(n1pd, d1pd);
		//q2 = _mm512_cvt_roundpd_epu64(n1pd,
		//	(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

		n1pd = _mm512_div_round_pd(n1pd, d1pd, ROUNDING_MODE);
		q2 = _mm512_cvttpd_epu64(n1pd);

		q = _mm512_mask_add_epi64(q, err, q, q2);
		r = _mm512_mask_sub_epi64(r, err, r, _mm512_mullox_epi64(q2, d));
	}

	return r;
}

// a Fermat PRP test on 8x 52-bit inputs
uint8_t fermat_prp_52x8(uint64_t* n)
{
	// assumes has no small factors.  assumes n <= 52 bits.
	// assumes n is a list of 8 52-bit integers
	// do a base-2 fermat prp test on each using LR binexp.
	__m512i vrho = multiplicative_inverse104_x8(n);
	__m512i unity;
	__m512i r;
	__m512i nvec;
	__m512i evec;
	__m512i m;
	__m512i zero = _mm512_setzero_si512();
	__m512i one = _mm512_set1_epi64(1);

	vrho = _mm512_and_epi64(_mm512_sub_epi64(zero, vrho), lo52mask);
	nvec = loadu64(n);
	evec = _mm512_sub_epi64(nvec, one);

#if defined(INTEL_COMPILER) || defined(INTEL_LLVM_COMPILER)
	r = _mm512_rem_epu64(_mm512_set1_epi64(1ULL << 52), nvec);
#else
	r = rem_epu64_x8(_mm512_set1_epi64(1ULL << 52), nvec);
#endif

	// penultimate-hi-bit mask
	m = _mm512_sub_epi64(_mm512_set1_epi64(62), _mm512_lzcnt_epi64(evec));
	m = _mm512_sllv_epi64(_mm512_set1_epi64(1), m);

	// we know the first bit is set and the first squaring is of unity,
	// so we can do the first iteration manually with no squaring.
	unity = r;

	r = _mm512_add_epi64(r, r);
	__mmask8 ge = _mm512_cmpge_epi64_mask(r, nvec);
	r = _mm512_mask_sub_epi64(r, ge, r, nvec);

	while (_mm512_cmpgt_epu64_mask(m, zero))
	{
		__mmask8 bitcmp = _mm512_test_epi64_mask(m, evec);
		mulredc52_mask_add_vec(&r, bitcmp, r, r, nvec, vrho);
		m = _mm512_srli_epi64(m, 1);
	}

	// AMM possibly needs a final correction by n
	ge = _mm512_cmpge_epi64_mask(r, nvec);
	r = _mm512_mask_sub_epi64(r, ge, r, nvec);

	return _mm512_cmpeq_epu64_mask(unity, r);
}

// a Fermat PRP test on 8x 104-bit inputs
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

// a Miller-Rabin SPRP test on 8x 52-bit inputs using base 2
uint8_t MR_2sprp_52x8(uint64_t* n)
{
	// assumes has no small factors.  assumes n <= 52 bits.
	// assumes n is a list of 8 52-bit integers
	// do a base-2 MR sprp test on each using LR binexp.
	__m512i vrho = multiplicative_inverse104_x8(n);
	__m512i unity;
	__m512i r;
	__m512i nvec;
	__m512i evec;
	__m512i m;
	__m512i zero = _mm512_setzero_si512();
	__m512i one = _mm512_set1_epi64(1);

	vrho = _mm512_and_epi64(_mm512_sub_epi64(zero, vrho), lo52mask);
	nvec = loadu64(n);
	evec = _mm512_sub_epi64(nvec, one);

#if defined(INTEL_COMPILER) || defined(INTEL_LLVM_COMPILER)
	r = _mm512_rem_epu64(_mm512_set1_epi64(1ULL << 52), nvec);
#else
	r = rem_epu64_x8(_mm512_set1_epi64(1ULL << 52), nvec);
#endif

	// penultimate-hi-bit mask
	m = _mm512_sub_epi64(_mm512_set1_epi64(62), _mm512_lzcnt_epi64(evec));
	m = _mm512_sllv_epi64(_mm512_set1_epi64(1), m);

	// we know the first bit is set and the first squaring is of unity,
	// so we can do the first iteration manually with no squaring.
	unity = r;

	r = _mm512_add_epi64(r, r);
	__mmask8 ge = _mm512_cmpge_epi64_mask(r, nvec);
	r = _mm512_mask_sub_epi64(r, ge, r, nvec);

	while (_mm512_cmpgt_epu64_mask(m, zero))
	{
		__mmask8 bitcmp = _mm512_test_epi64_mask(m, evec);
		mulredc52_mask_add_vec(&r, bitcmp, r, r, nvec, vrho);
		m = _mm512_srli_epi64(m, 1);
	}

	// AMM possibly needs a final correction by n
	ge = _mm512_cmpge_epi64_mask(r, nvec);
	r = _mm512_mask_sub_epi64(r, ge, r, nvec);

	return _mm512_cmpeq_epu64_mask(unity, r);
}

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

	nv[0] = loadu64(&n[0]);
	nv[1] = loadu64(&n[8]);

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

	mone[0] = loadu64(unity.data[0]);
	mone[1] = loadu64(unity.data[1]);

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

// a Miller-Rabin SPRP test on 8x 104-bit inputs using an
// independent arbitrary 52-bit base on each
uint8_t MR_sprp_104x8(uint64_t* n, uint64_t* bases)
{
	// assumes has no small factors.  assumes n >= 54 bits.
	// assumes n is a list of 8 104-bit integers (16 52-bit words)
	// in the format: 8 lo-words, 8 hi-words.
	// assume bases is a list of 8 small (single-word) bases, one for each input n.
	// do a Miller-Rabin sprp test on each using the supplied bases.
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

	nv[0] = loadu64(&n[0]);
	nv[1] = loadu64(&n[8]);

	// the 128-bit division we do one at a time
	int i;
	for (i = 0; i < 8; i++)
	{
		uint128_t mod = ((uint128_t)n[i + 8] << 52) + n[i];
		uint128_t one = (uint128_t)1 << 104;
		one %= mod;

		unity.data[0][i] = (uint64_t)one & 0xfffffffffffffULL;
		unity.data[1][i] = (uint64_t)(one >> 52);
	}

	mone[0] = loadu64(unity.data[0]);
	mone[1] = loadu64(unity.data[1]);

	// get bases into Monty rep
	bv[0] = loadu64(bases);
	bv[1] = zerov;

	__m512i mpow[2];
	mpow[0] = mone[0];
	mpow[1] = mone[1];

	rv[0] = mone[0];
	rv[1] = mone[1];

	bv[0] = _mm512_srli_epi64(bv[0], 1);
	__mmask8 done = _mm512_cmpeq_epi64_mask(bv[0], zerov);
	while (done != 0xff)
	{
		addmod104_x8(&mpow[1], &mpow[0], mpow[1], mpow[0], mpow[1], mpow[0], nv[1], nv[0]);
		__mmask8 bitcmp = _mm512_test_epi64_mask(onev, bv[0]);
		mask_addmod104_x8(&rv[1], &rv[0], (~done) & bitcmp, rv[1], rv[0], mpow[1], mpow[0], nv[1], nv[0]);

		bv[0] = _mm512_srli_epi64(bv[0], 1);
		done = _mm512_cmpeq_epi64_mask(bv[0], zerov);
	}

	bv[0] = rv[0];
	bv[1] = rv[1];

	// compute d and tzcnt
	submod104_x8(&n1v[1], &n1v[0], nv[1], nv[0], zerov, onev, nv[1], nv[0]);

	done = 0;
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
	rv[0] = bv[0];
	rv[1] = bv[1];

	// compute b^d
	done = _mm512_cmpeq_epu64_mask(m, zerov);
	while (done != 0xff)
	{
		__mmask8 bitcmp = _mm512_test_epi64_mask(m, dv[1]);

		mask_sqrredc104_vec(&rv[1], &rv[0], ~done, rv[1], rv[0], nv[1], nv[0], vrho);
		mask_mulredc104_vec(&rv[1], &rv[0], (~done) & bitcmp, rv[1], rv[0], bv[1], bv[0], nv[1], nv[0], vrho);

		m = _mm512_srli_epi64(m, 1);
		done = _mm512_cmpeq_epu64_mask(m, zerov);
	}

	m = _mm512_sub_epi64(_mm512_set1_epi64(62), _mm512_lzcnt_epi64(dv[0]));
	m = _mm512_sllv_epi64(_mm512_set1_epi64(1), m);
	m = _mm512_mask_set1_epi64(m, _mm512_cmpge_epi64_mask(dv[1], onev), 1ULL << 51);

	done = _mm512_cmpeq_epu64_mask(m, zerov);
	while (done != 0xff)
	{
		__mmask8 bitcmp = _mm512_test_epi64_mask(m, dv[0]);

		mask_sqrredc104_vec(&rv[1], &rv[0], ~done, rv[1], rv[0], nv[1], nv[0], vrho);
		mask_mulredc104_vec(&rv[1], &rv[0], (~done) & bitcmp, rv[1], rv[0], bv[1], bv[0], nv[1], nv[0], vrho);

		m = _mm512_srli_epi64(m, 1);
		done = _mm512_cmpeq_epu64_mask(m, zerov);
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
		addmod104_x8(&rv[1], &rv[0], rv[1], rv[0], zerov, zerov, nv[1], nv[0]);

		// definitely not prp by 1 check, stop checking
		done |= (_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[1], zerov) &
			_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[0], onev));

		done |= _mm512_mask_cmple_epu64_mask(~done, tzcntv, onev);
	}

	addmod104_x8(&rv[1], &rv[0], rv[1], rv[0], zerov, zerov, nv[1], nv[0]);

	// check current result == m-1
	ism1prp |= (_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[1], n1v[1]) &
		_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[0], n1v[0]));

	return (is1prp | ism1prp);
}

// a Miller-Rabin SPRP test on 1 104-bit input using 8x
// different bases
uint8_t MR_sprp_104x8base(uint64_t* n, uint64_t* one, uint64_t* bases)
{
	// assumes has no small factors and is odd.  assumes n >= 54 bits.
	// assumes n is a 104-bit integer with two 52 bit words: [lo,hi].
	// assumes one is a 104-bit integer equal to (1 << 104) mod n.
	// assume bases is a list of 8 small (single-word) bases.
	// do a Miller-Rabin sprp test on the input to each supplied base.
	// uint128_t n128 = ((uint128_t)n[1] << 52) + (uint128_t)n[0];
	__m512i vrho = _mm512_set1_epi64(multiplicative_inverse(n[0]));
	__m512i mone[2];
	__m512i nv[2];
	__m512i dv[2];
	__m512i rv[2];
	__m512i bv[2];
	__m512i n1v[2];
	__m512i tv[2];
	__m512i zerov = _mm512_setzero_si512();
	__m512i onev = _mm512_set1_epi64(1);
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
	uint64_t tmp = 0;

	vrho = _mm512_and_epi64(_mm512_sub_epi64(zerov, vrho), lo52mask);

	nv[0] = _mm512_set1_epi64(n[0]);
	nv[1] = _mm512_set1_epi64(n[1]);

	mone[0] = _mm512_set1_epi64(one[0]);
	mone[1] = _mm512_set1_epi64(one[1]);

	// get bases into Monty rep
	bv[0] = loadu64(bases);
	bv[1] = zerov;

	__m512i mpow[2];
	mpow[0] = mone[0];
	mpow[1] = mone[1];

	rv[0] = mone[0];
	rv[1] = mone[1];

	bv[0] = _mm512_srli_epi64(bv[0], 1);

	__mmask8 done = _mm512_cmpeq_epi64_mask(bv[0], zerov);
	while (done != 0xff)
	{
		addmod104_x8(&mpow[1], &mpow[0], mpow[1], mpow[0], mpow[1], mpow[0], nv[1], nv[0]);
		__mmask8 bitcmp = _mm512_test_epi64_mask(onev, bv[0]);
		mask_addmod104_x8(&rv[1], &rv[0], (~done) & bitcmp, rv[1], rv[0], mpow[1], mpow[0], nv[1], nv[0]);

		bv[0] = _mm512_srli_epi64(bv[0], 1);
		done = _mm512_cmpeq_epi64_mask(bv[0], zerov);
	}

	bv[0] = rv[0];
	bv[1] = rv[1];

	// compute d and tzcnt
	uint64_t d[2];
	d[1] = n[1];
	d[0] = n[0] - 1;			// n odd, so this won't carry
	int ntz = my_ctz104(d[0], d[1]);
	n1v[0] = _mm512_set1_epi64(d[0]);
	n1v[1] = _mm512_set1_epi64(d[1]);
	if (ntz < 52)
	{
		uint64_t shift = d[1] & ((1ULL << ntz) - 1);
		d[0] = (d[0] >> ntz) + (shift << (52 - ntz));
		d[1] >>= ntz;
	}
	else
	{
		d[0] = d[1];
		d[1] = 0;
		d[0] >>= (ntz - 52);
	}


#if 0
	// LR
	uint64_t m;
	if (d[1] <= 1) m = 0;
	else m = (1ull << (62 - my_clz52(d[1])));

	while (m > 0)
	{
		sqrredc104_vec(&rv[1], &rv[0], rv[1], rv[0], nv[1], nv[0], vrho);
		if (m & d[1])
			mask_mulredc104_vec(&rv[1], &rv[0], 0xff, rv[1], rv[0], bv[1], bv[0], nv[1], nv[0], vrho);
		m >>= 1;
	}

	if (d[1] == 0) m = (1ull << (62 - my_clz52(d[0])));
	else m = (1ull << 51);

	while (m > 0)
	{
		sqrredc104_vec(&rv[1], &rv[0], rv[1], rv[0], nv[1], nv[0], vrho);
		if (m & d[0])
			mask_mulredc104_vec(&rv[1], &rv[0], 0xff, rv[1], rv[0], bv[1], bv[0], nv[1], nv[0], vrho);
		m >>= 1;
	}
#elif 0
	// RL-104
	uint128_t d128 = ((uint128_t)d[1] << 52) + (uint128_t)d[0];
	rv[0] = mone[0];
	rv[1] = mone[1];
	while (d128 > 0)
	{
		if (d128 & 1)
			mask_mulredc104_vec(&rv[1], &rv[0], 0xff,
				rv[1], rv[0], bv[1], bv[0], nv[1], nv[0], vrho);

		d128 >>= 1;

		if (d128)
			sqrredc104_vec(&bv[1], &bv[0], bv[1], bv[0], nv[1], nv[0], vrho);
	}
#elif 0
	// LR-kary
	uint64_t g[256];
	// 
	// precomputation
	rv[0] = bv[0];
	rv[1] = bv[1];

	int i;
	storeu64(&g[0 * 8], mone[0]);		// g0 = 1
	storeu64(&g[1 * 8], mone[1]);		// g0 = 1
	storeu64(&g[2 * 8], rv[0]);			// g1 = g 
	storeu64(&g[3 * 8], rv[1]);			// g1 = g 
	sqrredc104_vec(&rv[1], &rv[0], rv[1], rv[0], nv[1], nv[0], vrho);
	storeu64(&g[4 * 8], rv[0]);			// g2 = g^2
	storeu64(&g[5 * 8], rv[1]);			// g2 = g^2
	for (i = 3; i < 16; i++)
	{
		mask_mulredc104_vec(&rv[1], &rv[0], 0xff, rv[1], rv[0], bv[1], bv[0], nv[1], nv[0], vrho);
		storeu64(&g[(i * 2) * 8], rv[0]);			// gi = g^i
		storeu64(&g[(i * 2 + 1) * 8], rv[1]);		// gi = g^i
	}

	uint128_t d128 = ((uint128_t)d[1] << 52) + (uint128_t)d[0];
	rv[0] = mone[0];
	rv[1] = mone[1];
	int lz = my_clz104(d[0], d[1]);
	int msb = 112 - lz;
	int m;

	m = (d128 >> msb) & 0xf;
	rv[0] = loadu64(&g[(2 * m) * 8]);
	rv[1] = loadu64(&g[(2 * m + 1) * 8]);
	msb -= 4;

	while (msb > 0)
	{
		sqrredc104_vec(&rv[1], &rv[0], rv[1], rv[0], nv[1], nv[0], vrho);
		sqrredc104_vec(&rv[1], &rv[0], rv[1], rv[0], nv[1], nv[0], vrho);
		sqrredc104_vec(&rv[1], &rv[0], rv[1], rv[0], nv[1], nv[0], vrho);
		sqrredc104_vec(&rv[1], &rv[0], rv[1], rv[0], nv[1], nv[0], vrho);
		m = (d128 >> msb) & 0xf;

		if (m > 0)
			mask_mulredc104_vec(&rv[1], &rv[0], 0xff,
				loadu64(&g[(2 * m + 1) * 8]), loadu64(&g[(2 * m) * 8]),
				rv[1], rv[0], nv[1], nv[0], vrho);

		msb -= 4;
	}

	msb += 4;
	m = (int)d128 & ((1 << msb) - 1);
	while (msb > 0)
	{
		sqrredc104_vec(&rv[1], &rv[0], rv[1], rv[0], nv[1], nv[0], vrho);
		msb--;
	}

	mask_mulredc104_vec(&rv[1], &rv[0], 0xff,
		loadu64(&g[(2 * m + 1) * 8]), loadu64(&g[(2 * m) * 8]),
		rv[1], rv[0], nv[1], nv[0], vrho);
#else
	// RL-52x2
	rv[0] = mone[0];
	rv[1] = mone[1];
	int i = 0;
	while (d[0] > 0)
	{
		if (d[0] & 1)
			mask_mulredc104_vec(&rv[1], &rv[0], 0xff,
				rv[1], rv[0], bv[1], bv[0], nv[1], nv[0], vrho);

		d[0] >>= 1;
		i++;

		sqrredc104_vec(&bv[1], &bv[0], bv[1], bv[0], nv[1], nv[0], vrho);
	}

	for (; (i < 52) && (d[1] > 0); i++)
	{
		sqrredc104_vec(&bv[1], &bv[0], bv[1], bv[0], nv[1], nv[0], vrho);
	}

	while (d[1] > 0)
	{
		if (d[1] & 1)
			mask_mulredc104_vec(&rv[1], &rv[0], 0xff,
				rv[1], rv[0], bv[1], bv[0], nv[1], nv[0], vrho);

		d[1] >>= 1;

		if (d[1])
			sqrredc104_vec(&bv[1], &bv[0], bv[1], bv[0], nv[1], nv[0], vrho);
	}

#endif

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

	while ((done != 0xff) && (ntz > 0))
	{
		ntz--;

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
	}

	addmod104_x8(&rv[1], &rv[0], zerov, zerov, rv[1], rv[0], nv[1], nv[0]);

	// check current result == m-1
	ism1prp |= (_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[1], n1v[1]) &
		_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[0], n1v[0]));

	return (is1prp | ism1prp);
}

int main(int argc, char** argv)
{
	uint64_t prp[16];
	int correct = 0;
	int i, k;

#ifdef GMP_CHECK
	mpz_t gmp2, gmpn, gmpn1;
	mpz_init(gmp2);
	mpz_init(gmpn);
	mpz_init(gmpn1);
#endif

#ifndef IFMA
	dbias = _mm512_castsi512_pd(set64(0x4670000000000000ULL));
	vbias1 = set64(0x4670000000000000ULL);
	vbias2 = set64(0x4670000000000001ULL);
	vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);
#endif

	lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);

	struct timeval start, stop;
	int bits;
	uint64_t elapsed;

	// test of fermat_prp_64x1 on random 6k+1 inputs
	for (bits = 20; bits <= 0; bits += 1)
	{
		uint32_t numprp = 0;
		uint64_t ticks1 = my_rdtsc();
		uint64_t ticks2;
		uint32_t num = 1000000;
		double telapsed = 0;
		int k;

		numprp = 0;
		k = 0;
		elapsed = 0;
		telapsed = 0;
		do {

			uint128_t x;
			do {
				x = my_random();
				uint128_t maskAnd = ((uint128_t)1 << (bits - 1)) - 1;	// clear msbits
				uint128_t maskOr = ((uint128_t)1 << (bits - 1)) | ((uint128_t)1 << (bits / 2));	// force msb, force another bit
				x &= maskAnd;
				x |= maskOr;
				x /= 6;
				x *= 6;	// now a multiple of 6
				x += 1;	// number like 6*k + 1
			} while (x >> (bits - 1) != 1);

			prp[0] = (uint64_t)x;

			ticks1 = my_rdtsc();
			gettimeofday(&start, NULL);

			uint64_t inc = 4;
			int j;

			for (j = 0; j < num; j++)
			{
				numprp += fermat_prp_64x1(prp[0]);
				prp[0] += inc;
				inc = 6 - inc;
			}

			k++;
			ticks2 = my_rdtsc();
			elapsed += (ticks2 - ticks1);
			gettimeofday(&stop, NULL);
			telapsed += _difftime(&start, &stop);

		} while (elapsed < (1ull << 30));

		printf("total ticks = %lu, ticks per %d-bit input = %lu\n",
			elapsed, bits, (elapsed) / (k * num));
		printf("found %d fermat-prp out of %u %d-bit inputs: %1.2f%%\n",
			numprp, k * num, bits, 100. * (double)numprp / (double)(k * num));
		printf("elapsed time: %1.4f sec, %1.4f us / input\n", telapsed, 1000000. * telapsed / (double)(k * num));
	}
	printf("\n");

	// test of fermat_prp_128x1 on random 6k+1 inputs
	for (bits = 118; bits <= 0; bits += 1)
	{
		uint32_t numprp = 0;
		uint64_t ticks1 = my_rdtsc();
		uint64_t ticks2;
		uint32_t num = 1000000;
		double telapsed = 0;
		int k;

		numprp = 0;
		k = 0;
		elapsed = 0;
		telapsed = 0;
		do {

			uint128_t x;
			do {
				x = my_random();
				uint128_t maskAnd = ((uint128_t)1 << (bits - 1)) - 1;	// clear msbits
				uint128_t maskOr = ((uint128_t)1 << (bits - 1)) | ((uint128_t)1 << (bits / 2));	// force msb, force another bit
				x &= maskAnd;
				x |= maskOr;
				x /= 6;
				x *= 6;	// now a multiple of 6
				x += 1;	// number like 6*k + 1
			} while (x >> (bits - 1) != 1);

			prp[0] = (uint64_t)x;
			prp[1] = (uint64_t)(x >> 64);

			ticks1 = my_rdtsc();
			gettimeofday(&start, NULL);

			uint64_t inc = 4;
			int j;

			for (j = 0; j < num; j++)
			{
				numprp += fermat_prp_128x1(prp);
				prp[0] += inc;
				inc = 6 - inc;
			}

			k++;
			ticks2 = my_rdtsc();
			elapsed += (ticks2 - ticks1);
			gettimeofday(&stop, NULL);
			telapsed += _difftime(&start, &stop);

		} while (elapsed < (1ull << 30));

		printf("total ticks = %lu, ticks per %d-bit input = %lu\n",
			elapsed, bits, (elapsed) / (k * num));
		printf("found %d fermat-prp out of %u %d-bit inputs: %1.2f%%\n",
			numprp, k * num, bits, 100. * (double)numprp / (double)(k * num));
		printf("elapsed time: %1.4f sec, %1.4f us / input\n", telapsed, 1000000. * telapsed / (double)(k * num));
	}
	printf("\n");

	// test of fermat_prp_52x8 on random 6k+1 inputs
	for (bits = 20; bits <= 0; bits += 1)
	{
		uint32_t numprp = 0;
		uint64_t ticks1 = my_rdtsc();
		uint64_t ticks2;
		uint32_t num = 1000000;
		double telapsed = 0;
		int k;
		if (bits > 52) bits = 52;

		numprp = 0;
		k = 0;
		elapsed = 0;
		telapsed = 0;
		do {
			for (i = 0; i < 8; i++)
			{
				uint64_t x;
				do {
					x = my_random();
					uint128_t maskAnd = ((uint128_t)1 << (bits - 1)) - 1;	// clear msbits
					uint128_t maskOr = ((uint128_t)1 << (bits - 1)) | ((uint128_t)1 << (bits / 2));	// force msb, force another bit
					x &= maskAnd;
					x |= maskOr;
					x /= 6;
					x *= 6;	// now a multiple of 6
					x += 1;	// number like 6*k + 1
				} while (x >> (bits - 1) != 1);
				prp[i] = (uint64_t)x & 0xfffffffffffffull;
			}

			ticks1 = my_rdtsc();
			gettimeofday(&start, NULL);

			uint64_t inc = 4;
			int j;

			for (j = 0; j < num; j++)
			{
				numprp += _mm_popcnt_u32(fermat_prp_52x8(prp));
				for (i = 0; i < 8; i++)
				{
					prp[i] += inc;
				}
				inc = 6 - inc;
			}

			k++;
			ticks2 = my_rdtsc();
			elapsed += (ticks2 - ticks1);
			gettimeofday(&stop, NULL);
			telapsed += _difftime(&start, &stop);

		} while (elapsed < (1ull << 30));

		printf("total ticks = %lu, ticks per %d-bit input = %lu\n",
			elapsed, bits, (elapsed) / (k * num * 8));
		printf("found %d fermat-prp out of %u %d-bit inputs: %1.2f%%\n",
			numprp, k * num * 8, bits, 100. * (double)numprp / (double)(k * num * 8));
		printf("elapsed time: %1.4f sec, %1.4f us / input\n", telapsed, 1000000. * telapsed / (double)(k * num * 8));
	}
	printf("\n");

	// test of fermat_prp_104x8 on random 6k+1 inputs
	for (bits = 80; bits <= 104; bits += 1)
	{
		uint32_t numprp = 0;
		uint64_t ticks1 = my_rdtsc();
		uint64_t ticks2;
		uint32_t num = 100000;
		double telapsed = 0;
		int k;
		if (bits > 104) bits = 104;

		k = 0;
		int fail = 0;
		elapsed = 0;
		do {
			for (i = 0; i < 8; i++)
			{
				uint128_t x;
				do {
					x = my_random();
					uint128_t maskAnd = ((uint128_t)1 << (bits - 1)) - 1;	// clear msbits
					uint128_t maskOr = ((uint128_t)1 << (bits - 1)) | ((uint128_t)1 << (bits / 2));	// force msb, force another bit
					x &= maskAnd;
					x |= maskOr;
					x /= 6;
					x *= 6;	// now a multiple of 6
					x += 1;	// number like 6*k + 1
				} while (x >> (bits - 1) != 1);
				prp[i] = (uint64_t)x & 0xfffffffffffffull;
				prp[i + 8] = (uint64_t)(x >> 52) & 0xfffffffffffffull;
			}

			uint64_t inc = 4;
			int j;

			ticks1 = my_rdtsc();
			gettimeofday(&start, NULL);

			for (j = 0; j < num; j++)
			{
				uint8_t prpmask = fermat_prp_104x8(prp);
				numprp += _mm_popcnt_u32(prpmask);
				for (i = 0; i < 8; i++)
				{
#ifdef GMP_CHECK
					if ((prpmask & (1 << i)) == 0)
					{
						mpz_set_ui(gmp2, 2);
						mpz_set_ui(gmpn, prp[8 + i]);
						mpz_mul_2exp(gmpn, gmpn, 52);
						mpz_add_ui(gmpn, gmpn, prp[i]);
						mpz_sub_ui(gmpn1, gmpn, 1);
						mpz_powm(gmp2, gmp2, gmpn1, gmpn);
						if (mpz_cmp_ui(gmp2, 1) == 0)
						{
							//printf("prp %016lx%016lx failed in lane %d\n", prp[i+8], prp[i], i);
							//gmp_printf("mpz result = %Zx\n", gmp2);
							//exit(1);
							fail++;
						}
					}
#endif
					prp[i] += inc;
				}
				inc = 6 - inc;
			}

			k++;
			ticks2 = my_rdtsc();
			elapsed += (ticks2 - ticks1);
			gettimeofday(&stop, NULL);
			telapsed = +_difftime(&start, &stop);
		} while (elapsed < (1ull << 30));

		printf("total ticks = %lu, ticks per %d-bit input = %lu\n",
			elapsed, bits, (elapsed) / (k * num * 8));
		printf("found %d fermat-prp out of %u %d-bit inputs: %1.2f%%\n",
			numprp, k * num * 8, bits, 100. * (double)numprp / (double)(k * num * 8));
#ifdef GMP_CHECK
		printf("GMP checked number of failures: %d\n", fail);
#endif
		printf("elapsed time: %1.4f sec, %1.4f us / input\n", telapsed, 1000000. * telapsed / (double)(k * num * 8));
	}
	printf("\n");

	// test of fermat_prp_52x8 on random 6k+1 inputs
	for (bits = 20; bits <= 0; bits += 1)
	{
		uint32_t numprp = 0;
		uint64_t ticks1 = my_rdtsc();
		uint64_t ticks2;
		uint32_t num = 1000000;
		double telapsed = 0;
		int k;
		if (bits > 52) bits = 52;

		numprp = 0;
		k = 0;
		elapsed = 0;
		telapsed = 0;
		do {
			for (i = 0; i < 8; i++)
			{
				uint64_t x;
				do {
					x = my_random();
					uint128_t maskAnd = ((uint128_t)1 << (bits - 1)) - 1;	// clear msbits
					uint128_t maskOr = ((uint128_t)1 << (bits - 1)) | ((uint128_t)1 << (bits / 2));	// force msb, force another bit
					x &= maskAnd;
					x |= maskOr;
					x /= 6;
					x *= 6;	// now a multiple of 6
					x += 1;	// number like 6*k + 1
				} while (x >> (bits - 1) != 1);
				prp[i] = (uint64_t)x & 0xfffffffffffffull;
			}

			ticks1 = my_rdtsc();
			gettimeofday(&start, NULL);

			uint64_t inc = 4;
			int j;

			for (j = 0; j < num; j++)
			{
				numprp += _mm_popcnt_u32(MR_2sprp_52x8(prp));
				for (i = 0; i < 8; i++)
				{
					prp[i] += inc;
				}
				inc = 6 - inc;
			}

			k++;
			ticks2 = my_rdtsc();
			elapsed += (ticks2 - ticks1);
			gettimeofday(&stop, NULL);
			telapsed += _difftime(&start, &stop);

		} while (elapsed < (1ull << 30));

		printf("total ticks = %lu, ticks per %d-bit input = %lu\n",
			elapsed, bits, (elapsed) / (k * num * 8));
		printf("found %d MR-2sprp out of %u %d-bit inputs: %1.2f%%\n",
			numprp, k * num * 8, bits, 100. * (double)numprp / (double)(k * num * 8));
		printf("elapsed time: %1.4f sec, %1.4f us / input\n", telapsed, 1000000. * telapsed / (double)(k * num * 8));
	}
	printf("\n");

	// test of MR_2sprp_104x8 on random 6k+1 inputs
	for (bits = 50; bits <= 0; bits += 1)
	{
		uint32_t numprp = 0;
		uint64_t ticks1 = my_rdtsc();
		uint64_t ticks2;
		uint32_t num = 100000;
		double telapsed = 0;
		int k;
		if (bits > 104) bits = 104;

		k = 0;
		elapsed = 0;
		int fail = 0;
		do {
			for (i = 0; i < 8; i++)
			{
				uint128_t x;
				do {
					x = my_random();
					uint128_t maskAnd = ((uint128_t)1 << (bits - 1)) - 1;	// clear msbits
					uint128_t maskOr = ((uint128_t)1 << (bits - 1)) | ((uint128_t)1 << (bits / 2));	// force msb, force another bit
					x &= maskAnd;
					x |= maskOr;
					x /= 6;
					x *= 6;	// now a multiple of 6
					x += 1;	// number like 6*k + 1
				} while (x >> (bits - 1) != 1);
				prp[i] = (uint64_t)x & 0xfffffffffffffull;
				prp[i + 8] = (uint64_t)(x >> 52) & 0xfffffffffffffull;
			}

			uint64_t inc = 4;
			int j;

			ticks1 = my_rdtsc();
			gettimeofday(&start, NULL);

			for (j = 0; j < num; j++)
			{
				uint8_t prpmask = MR_2sprp_104x8(prp);
				numprp += _mm_popcnt_u32(prpmask);
				for (i = 0; i < 8; i++)
				{
#ifdef GMP_CHECK
					if ((prpmask & (1 << i)) == 0)
					{
						mpz_set_ui(gmp2, 2);
						mpz_set_ui(gmpn, prp[8 + i]);
						mpz_mul_2exp(gmpn, gmpn, 52);
						mpz_add_ui(gmpn, gmpn, prp[i]);
						mpz_sub_ui(gmpn1, gmpn, 1);
						mpz_powm(gmp2, gmp2, gmpn1, gmpn);
						if (mpz_cmp_ui(gmp2, 1) == 0)
						{
							//printf("prp %016lx%016lx failed in lane %d\n", prp[i + 8], prp[i], i);
							//gmp_printf("mpz result = %Zx\n", gmp2);
							//exit(1);
							fail++;
						}
					}
#endif
					prp[i] += inc;
				}
				inc = 6 - inc;
			}

			k++;
			ticks2 = my_rdtsc();
			elapsed += (ticks2 - ticks1);
			gettimeofday(&stop, NULL);
			telapsed = +_difftime(&start, &stop);
		} while (elapsed < (1ull << 30));

		printf("total ticks = %lu, ticks per %d-bit input = %lu\n",
			elapsed, bits, (elapsed) / (k * num * 8));
		printf("found %d MR-2sprp out of %u %d-bit inputs: %1.2f%%\n",
			numprp, k * num * 8, bits, 100. * (double)numprp / (double)(k * num * 8));
#ifdef GMP_CHECK
		printf("GMP checked number of failures: %d\n", fail);
#endif
		printf("elapsed time: %1.4f sec, %1.4f us / input\n", telapsed, 1000000. * telapsed / (double)(k * num * 8));
	}
	printf("\n");

	// bases for MR-sprp check:
	uint64_t bases[24] = { 3, 5, 7, 11,
		13, 17, 19, 23,
		29, 31, 37, 41,
		43, 47, 53, 59,
		61, 67, 71, 73,
		79, 83, 89, 97 };

	// test of MR_sprp_104x8 on PRP 6k+1 inputs
	for (bits = 50; bits <= 0; bits += 1)
	{
		uint32_t numprp = 0;
		uint64_t ticks1 = my_rdtsc();
		uint64_t ticks2;
		uint32_t num = 100000;
		double telapsed = 0;

		if (bits > 104) bits = 104;

		//printf("commencing test of random 6k+1 %d-bit inputs\n", bits);
		elapsed = 0;

		uint64_t inc[8] = { 4, 4, 4, 4, 4, 4, 4, 4 };

		for (i = 0; i < 8; i++)
		{
			uint128_t x;
			do {
				x = my_random();
				uint128_t maskAnd = ((uint128_t)1 << (bits - 1)) - 1;	// clear msbits
				uint128_t maskOr = ((uint128_t)1 << (bits - 1)) | ((uint128_t)1 << (bits / 2));	// force msb, force another bit
				x &= maskAnd;
				x |= maskOr;
				x /= 6;
				x *= 6;	// now a multiple of 6
				x += 1;	// number like 6*k + 1
			} while (x >> (bits - 1) != 1);
			prp[i] = (uint64_t)x & 0xfffffffffffffull;
			prp[i + 8] = (uint64_t)(x >> 52);
		}

		uint8_t isprp = 0;
		while (isprp != 0xff)
		{
			isprp = fermat_prp_104x8(prp);

			for (i = 0; i < 8; i++)
			{
				if ((isprp & (1 << i)) == 0)
				{
					prp[i] += inc[i];
					inc[i] = 6 - inc[i];
				}
			}
		}

		ticks1 = my_rdtsc();

		uint64_t basecount[8];
		uint64_t maxcount[8];
		uint64_t currentbase[8];
		for (i = 0; i < 8; i++)
		{
			//printf("prp%d = %016lx%016lx\n", i, prp[i+8], prp[i]); fflush(stdout);
			basecount[i] = 0;
			currentbase[i] = bases[0];
			if (bits <= 62)
				maxcount[i] = 8;
			else if (bits <= 82)
				maxcount[i] = 16;	// was 12, but need a multiple of 8 (same cost anyway with avx512)
			else if (bits <= 112)
				maxcount[i] = 16;
			else
				maxcount[i] = 24;	// was 20, but need a multiple of 8 (same cost anyway with avx512)
		}

		uint32_t tested = 0;
		gettimeofday(&start, NULL);

		elapsed = 0;
		int fail = 0;
		while ((elapsed < (1ull << 30)))
		{
			uint8_t prpmask = MR_sprp_104x8(prp, currentbase);
			for (i = 0; i < 8; i++)
			{
#ifdef GMP_CHECK
				if ((prpmask & (1 << i)) == 0)
				{
					mpz_set_ui(gmpn, prp[8 + i]);
					mpz_mul_2exp(gmpn, gmpn, 52);
					mpz_add_ui(gmpn, gmpn, prp[i]);
					if (mpz_probab_prime_p(gmpn, 1) > 0)
					{
						//printf("prp %016lx%016lx failed in lane %d\n", prp[i + 8], prp[i], i);
						//gmp_printf("mpz result = %Zx\n", gmp2);
						//exit(1);
						fail++;
					}
				}
#endif

				if (prpmask & (1 << i))
				{
					// the input in position i could be prp, increment the
					// base if there are more of them
					if (basecount[i] < maxcount[i])
					{
						basecount[i]++;
						currentbase[i] = bases[basecount[i]];
					}
					else
					{
						// we've tested enough bases to know this is prime
						tested++;
						numprp++;
						//prp[i] += inc[i];
						//inc[i] = 6 - inc[i];
						currentbase[i] = bases[0];
						basecount[i] = 0;
					}
				}
				else
				{
					// the input in position i is definitely not prime,
					// replace it and increment num tested
					tested++;
					//prp[i] += inc[i];
					//inc[i] = 6 - inc[i];
					currentbase[i] = bases[0];
					basecount[i] = 0;
				}
			}

			ticks2 = my_rdtsc();
			elapsed = (ticks2 - ticks1);
		}

		gettimeofday(&stop, NULL);
		telapsed = +_difftime(&start, &stop);

		printf("total ticks = %lu, ticks per %d-bit input = %lu\n",
			elapsed, bits, (elapsed) / tested);
		printf("found %d MR-sprp out of %u %d-bit inputs: %1.2f%%\n",
			numprp, tested, bits, 100. * (double)numprp / (double)tested);
#ifdef GMP_CHECK
		printf("GMP checked number of failures: %d\n", fail);
#endif
		printf("elapsed time: %1.4f sec, %1.4f us / input\n", telapsed, 1000000. * telapsed / (double)tested);
	}
	printf("\n");

	// test of MR_sprp_104x8base on PRP 6k+1 inputs
	for (bits = 50; bits <= 0; bits += 1)
	{
		uint32_t numprp = 0;
		uint64_t ticks1;
		uint64_t ticks2;
		uint32_t num = 100000;
		double telapsed = 0;
		uint64_t one[16];

		//printf("commencing test of random 6k+1 %d-bit inputs\n", bits);
		elapsed = 0;

		if (bits > 104) bits = 104;

		uint64_t inc[8] = { 4, 4, 4, 4, 4, 4, 4, 4 };

		for (i = 0; i < 8; i++)
		{
			uint128_t x;
			do {
				x = my_random();
				uint128_t maskAnd = ((uint128_t)1 << (bits - 1)) - 1;	// clear msbits
				uint128_t maskOr = ((uint128_t)1 << (bits - 1)) | ((uint128_t)1 << (bits / 2));	// force msb, force another bit
				x &= maskAnd;
				x |= maskOr;
				x /= 6;
				x *= 6;	// now a multiple of 6
				x += 1;	// number like 6*k + 1
			} while (x >> (bits - 1) != 1);
			prp[i] = (uint64_t)x & 0xfffffffffffffull;
			prp[i + 8] = (uint64_t)(x >> 52);
		}

		uint8_t isprp = 0;
		while (isprp != 0xff)
		{
			isprp = fermat_prp_104x8(prp);

			for (i = 0; i < 8; i++)
			{
				//printf("%016lx%016lx : %u (%u)\n", prp[i+8], prp[i], isprp & (1 << i), isprp);
				if ((isprp & (1 << i)) == 0)
				{
					prp[i] += inc[i];
					inc[i] = 6 - inc[i];
				}
			}
		}

		uint128_t o128;
		uint128_t n128;
		for (i = 0; i < 8; i++)
		{
			//printf("prp%d = %016lx%016lx\n", i, prp[i+8], prp[i]); fflush(stdout);
			n128 = ((uint128_t)prp[i + 8] << 52) + prp[i];
			o128 = (uint128_t)1 << 104;
			o128 = o128 % n128;
			one[i] = (uint64_t)o128 & 0xfffffffffffffull;
			one[i + 8] = (uint64_t)(o128 >> 52);
			//printf("one%d = %016lx%016lx\n", i, one[i + 8], one[i]); fflush(stdout);
		}

		uint64_t basecount = 0;
		uint64_t maxcount;
		uint64_t currentbase[8];

		if (bits <= 62)
			maxcount = 8;
		else if (bits <= 82)
			maxcount = 16;	// was 12, but need a multiple of 8 (same cost anyway with avx512)
		else if (bits <= 112)
			maxcount = 16;
		else
			maxcount = 24;	// was 20, but need a multiple of 8 (same cost anyway with avx512)

		for (i = 0; i < 8; i++)
		{
			currentbase[i] = bases[i];
		}

		uint32_t tested = 0;
		gettimeofday(&start, NULL);

		ticks1 = my_rdtsc();

		elapsed = 0;
		int tnum = 0;
		while ((elapsed < (1ull << 30)))
		{
			uint64_t ntest[2], otest[2];
			ntest[1] = prp[tnum + 8];
			ntest[0] = prp[tnum];
			otest[1] = one[tnum + 8];
			otest[0] = one[tnum];

			// so far, simple RL is faster than kary-LR
			uint8_t prpmask = MR_sprp_104x8base(ntest, otest, currentbase);

			if (prpmask == 0xff)
			{
				// the input is prp to all current bases, increment the
				// base if there are more of them
				if ((basecount + 8) <= maxcount)
				{
					for (i = 0; i < 8; i++)
					{
						currentbase[i] = bases[basecount + i];
					}
					basecount += 8;
				}
				else
				{
					// we've tested enough bases to know this is prime
					tested++;
					numprp++;
					basecount = 0;
					for (i = 0; i < 8; i++)
					{
						currentbase[i] = bases[basecount + i];
					}
					tnum = (tnum + 1) & 7;
				}
			}
			else
			{
				// the input in position i is definitely not prime,
				// replace it and increment num tested
				tested++;
				basecount = 0;
				for (i = 0; i < 8; i++)
				{
					currentbase[i] = bases[basecount + i];
				}
				tnum = (tnum + 1) & 7;
			}

			ticks2 = my_rdtsc();
			elapsed = (ticks2 - ticks1);
		}

		gettimeofday(&stop, NULL);
		telapsed = +_difftime(&start, &stop);

		printf("total ticks = %lu, ticks per %d-bit input = %lu\n",
			elapsed, bits, (elapsed) / tested);
		printf("found %d MR-sprp out of %u %d-bit inputs: %1.2f%%\n",
			numprp, tested, bits, 100. * (double)numprp / (double)tested);
		printf("elapsed time: %1.4f sec, %1.4f us / input\n", telapsed, 1000000. * telapsed / (double)tested);
	}
	printf("\n");

#ifdef GMP_CHECK
	mpz_clear(gmpn);
	mpz_clear(gmpn1);
	mpz_clear(gmp2);
#endif
	return 0;
}