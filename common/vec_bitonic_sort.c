// MIT License
// 
// Copyright (c) 2025 Ben Buhrow
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


// fast bitonic sort using AVX-512 and test code to find
// collisions of N-bit keys in many contiguous medium-
// sized lists of random keys.
// Written by Ben Buhrow, July 2025
// example compile line:
// clang -O2 -g -march=icelake-client vec_bitonic_sort.c -o vecsort
// REQUIRES AVX-512 or will fail with illegal instruction 


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <immintrin.h>
#include <x86intrin.h>
#include <ytools.h>

#if defined(WIN32) || defined(_WIN64)
	#define WIN32_LEAN_AND_MEAN

	#include <windows.h>
#else
	#include <fcntl.h>
	#include <unistd.h>
	#include <sys/resource.h>
#endif

typedef uint8_t uint8;
typedef uint16_t uint16;
typedef uint32_t uint32;
typedef uint64_t uint64;

/* #define HAVE_PROF */
#ifdef HAVE_PROF
#define SHOW_PROF __attribute__((noinline))
#else
#define SHOW_PROF /* nothing */
#endif

static inline uint64_t my_rdtsc(void)
{
#if defined(__x86_64__)
	// supported by GCC and Clang for x86 platform
	return __rdtsc();
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

static uint32 
get_rand(uint32 *rand_seed, uint32 *rand_carry) {
   
	/* A multiply-with-carry generator by George Marsaglia.
	   The period is about 2^63. */

	#define RAND_MULT 2131995753

	uint64 temp;

	temp = (uint64)(*rand_seed) * 
		       (uint64)RAND_MULT + 
		       (uint64)(*rand_carry);
	*rand_seed = (uint32)temp;
	*rand_carry = (uint32)(temp >> 32);
	return (uint32)temp;
}

uint32_t my_clz32(uint64_t n)
{
#if (INLINE_ASM && defined(__x86_64__))
#ifdef __BMI1__
	uint32_t t;
	asm(" lzcnt %1, %0\n": "=r"(t) : "r"(n) : "flags");
	return t;
#else
	if (n)
		return __builtin_clz(n);
	return 32;
#endif
#else
#if defined(__GNUC__)
	if (n)
		return __builtin_clz(n);
	return 32;
#else
	if (n == 0)
		return 32;
	uint32_t r = 0;
	if ((n & (0xFFFFull << 16)) == 0)
		r += 16, n <<= 16;
	if ((n & (0xFFull << 24)) == 0)
		r += 8, n <<= 8;
	if ((n & (0xFull << 28)) == 0)
		r += 4, n <<= 4;
	if ((n & (0x3ull << 30)) == 0)
		r += 2, n <<= 2;
	if ((n & (0x1ull << 31)) == 0)
		r += 1;
	return r;
#endif
#endif
}

uint64_t my_ctz32(uint32_t n)
{
#if (INLINE_ASM && defined(__x86_64__))
#if defined(__BMI1__)
	uint32_t t;
	asm(" tzcnt %1, %0\n": "=r"(t) : "r"(n) : "flags");
	return t;
#else
	if (n)
		return __builtin_ctz(n);
	return 32;
#endif
#else
#if defined(__GNUC__)
	if (n)
		return __builtin_ctz(n);
	return 32;
#else
	if (n == 0)
		return 32;
	uint32_t r = 0;
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

uint32_t next_power_2(uint32_t sz)
{
	uint32_t new_sz = my_clz32(sz);
	if (new_sz == 0)
	{
		printf("buffer too big, sz must be <= 2^31 in sort()\n");
		exit(0);
	}
	return (1 << (32 - new_sz));
}

#if defined(USE_AVX512F)

// intrinsics for swapping N-bit chunks of data within a 512-bit vector
// that use immediates (faster and fewer registers than needing to load index vectors)
#define SWAP16(x) _mm512_rol_epi32((x), 16)
#define SWAP32(x) _mm512_shuffle_epi32((x), 0xB1)
#define SWAP64(x) _mm512_shuffle_epi32((x), 0x4E)
#define SWAP128(x) _mm512_permutex_epi64((x), 0x4E)
#define SWAP256(x) _mm512_shuffle_i64x2((x), (x), 0x4E)

// intrinsics for swapping N-bit chunks of data within a 512-bit vector, with N < 32.
// these require index vectors and the additional extension AVX-512BW
static __m512i swap8bit_idx;

#define SWAP8(x) _mm512_shuffle_epi8((x), swap8bit_idx)

#define SWAP16x2(v1, v2) \
	dv1_swap = SWAP16(v1);\
	dv2_swap = SWAP16(v2);
	
#define SWAP32x2(v1, v2) \
	dv1_swap = SWAP32(v1);\
	dv2_swap = SWAP32(v2);
	
#define SWAP64x2(v1, v2) \
	dv1_swap = SWAP64(v1);\
	dv2_swap = SWAP64(v2);

#define SWAP128x2(v1, v2) \
	dv1_swap = SWAP128(v1);\
	dv2_swap = SWAP128(v2);

#define SWAP256x2(v1, v2) \
	dv1_swap = SWAP256(v1);\
	dv2_swap = SWAP256(v2);
	
#define SWAP16x4(v1, v2, v3, v4) \
	dv1_swap = SWAP16(v1);\
	dv2_swap = SWAP16(v2);\
	dv3_swap = SWAP16(v3);\
	dv4_swap = SWAP16(v4);
	
#define SWAP32x4(v1, v2, v3, v4) \
	dv1_swap = SWAP32(v1);\
	dv2_swap = SWAP32(v2);\
	dv3_swap = SWAP32(v3);\
	dv4_swap = SWAP32(v4);
	
#define SWAP64x4(v1, v2, v3, v4) \
	dv1_swap = SWAP64(v1);\
	dv2_swap = SWAP64(v2);\
	dv3_swap = SWAP64(v3);\
	dv4_swap = SWAP64(v4);

#define SWAP128x4(v1, v2, v3, v4) \
	dv1_swap = SWAP128(v1);\
	dv2_swap = SWAP128(v2);\
	dv3_swap = SWAP128(v3);\
	dv4_swap = SWAP128(v4);

#define SWAP256x4(v1, v2, v3, v4) \
	dv1_swap = SWAP256(v1);\
	dv2_swap = SWAP256(v2);\
	dv3_swap = SWAP256(v3);\
	dv4_swap = SWAP256(v4);
	
#define SWAP16x8(v1, v2, v3, v4, v5, v6, v7, v8) \
	dv1_swap = SWAP16(v1);\
	dv2_swap = SWAP16(v2);\
	dv3_swap = SWAP16(v3);\
	dv4_swap = SWAP16(v4);\
	dv5_swap = SWAP16(v5);\
	dv6_swap = SWAP16(v6);\
	dv7_swap = SWAP16(v7);\
	dv8_swap = SWAP16(v8);
	
#define SWAP32x8(v1, v2, v3, v4, v5, v6, v7, v8) \
	dv1_swap = SWAP32(v1);\
	dv2_swap = SWAP32(v2);\
	dv3_swap = SWAP32(v3);\
	dv4_swap = SWAP32(v4);\
	dv5_swap = SWAP32(v5);\
	dv6_swap = SWAP32(v6);\
	dv7_swap = SWAP32(v7);\
	dv8_swap = SWAP32(v8);
	
#define SWAP64x8(v1, v2, v3, v4, v5, v6, v7, v8) \
	dv1_swap = SWAP64(v1);\
	dv2_swap = SWAP64(v2);\
	dv3_swap = SWAP64(v3);\
	dv4_swap = SWAP64(v4);\
	dv5_swap = SWAP64(v5);\
	dv6_swap = SWAP64(v6);\
	dv7_swap = SWAP64(v7);\
	dv8_swap = SWAP64(v8);

#define SWAP128x8(v1, v2, v3, v4, v5, v6, v7, v8) \
	dv1_swap = SWAP128(v1);\
	dv2_swap = SWAP128(v2);\
	dv3_swap = SWAP128(v3);\
	dv4_swap = SWAP128(v4);\
	dv5_swap = SWAP128(v5);\
	dv6_swap = SWAP128(v6);\
	dv7_swap = SWAP128(v7);\
	dv8_swap = SWAP128(v8);

#define SWAP256x8(v1, v2, v3, v4, v5, v6, v7, v8) \
	dv1_swap = SWAP256(v1);\
	dv2_swap = SWAP256(v2);\
	dv3_swap = SWAP256(v3);\
	dv4_swap = SWAP256(v4);\
	dv5_swap = SWAP256(v5);\
	dv6_swap = SWAP256(v6);\
	dv7_swap = SWAP256(v7);\
	dv8_swap = SWAP256(v8);

#define MINMAX16x2(m1, m2, v1, v2) \
	t1 = _mm512_mask_max_epu16(v1, m1, v1, dv1_swap);		\
	t2 = _mm512_mask_max_epu16(v2, m1, v2, dv2_swap);     \
	v1 = _mm512_mask_min_epu16(t1, m2, v1, dv1_swap);     \
	v2 = _mm512_mask_min_epu16(t2, m2, v2, dv2_swap);
	
#define MINMAX16x2_alt1(m1, m2, v1, v2) \
	t1 = _mm512_mask_max_epu16(v1, m1, v1, dv1_swap);		\
	t2 = _mm512_mask_min_epu16(v2, m1, v2, dv2_swap);     \
	v1 = _mm512_mask_min_epu16(t1, m2, v1, dv1_swap);     \
	v2 = _mm512_mask_max_epu16(t2, m2, v2, dv2_swap);
	
#define MINMAX16x2_alt2(m1, m2, v1, v2) \
	t1 = _mm512_mask_max_epu16(v1, m1, v1, dv1_swap);		\
	t2 = _mm512_mask_max_epu16(v2, m1, v2, dv2_swap);     \
	v1 = _mm512_mask_min_epu16(t1, m2, v1, dv1_swap);     \
	v2 = _mm512_mask_min_epu16(t2, m2, v2, dv2_swap);
	
#define MINMAX16x4(m1, m2, v1, v2, v3, v4) \
	t1 = _mm512_mask_max_epu16(v1, m1, v1, dv1_swap);		\
	t2 = _mm512_mask_max_epu16(v2, m1, v2, dv2_swap);     \
	t3 = _mm512_mask_max_epu16(v3, m1, v3, dv3_swap);     \
	t4 = _mm512_mask_max_epu16(v4, m1, v4, dv4_swap);     \
	v1 = _mm512_mask_min_epu16(t1, m2, v1, dv1_swap);     \
	v2 = _mm512_mask_min_epu16(t2, m2, v2, dv2_swap);     \
	v3 = _mm512_mask_min_epu16(t3, m2, v3, dv3_swap);     \
	v4 = _mm512_mask_min_epu16(t4, m2, v4, dv4_swap);
	
#define MINMAX16x4_alt1(m1, m2, v1, v2, v3, v4) \
	t1 = _mm512_mask_max_epu16(v1, m1, v1, dv1_swap);		\
	t2 = _mm512_mask_min_epu16(v2, m1, v2, dv2_swap);     \
	t3 = _mm512_mask_max_epu16(v3, m1, v3, dv3_swap);     \
	t4 = _mm512_mask_min_epu16(v4, m1, v4, dv4_swap);     \
	v1 = _mm512_mask_min_epu16(t1, m2, v1, dv1_swap);     \
	v2 = _mm512_mask_max_epu16(t2, m2, v2, dv2_swap);     \
	v3 = _mm512_mask_min_epu16(t3, m2, v3, dv3_swap);     \
	v4 = _mm512_mask_max_epu16(t4, m2, v4, dv4_swap);
	
#define MINMAX16x4_alt2(m1, m2, v1, v2, v3, v4) \
	t1 = _mm512_mask_max_epu16(v1, m1, v1, dv1_swap);		\
	t2 = _mm512_mask_max_epu16(v2, m1, v2, dv2_swap);     \
	t3 = _mm512_mask_min_epu16(v3, m1, v3, dv3_swap);     \
	t4 = _mm512_mask_min_epu16(v4, m1, v4, dv4_swap);     \
	v1 = _mm512_mask_min_epu16(t1, m2, v1, dv1_swap);     \
	v2 = _mm512_mask_min_epu16(t2, m2, v2, dv2_swap);     \
	v3 = _mm512_mask_max_epu16(t3, m2, v3, dv3_swap);     \
	v4 = _mm512_mask_max_epu16(t4, m2, v4, dv4_swap);
	
#define MINMAX16x8(m1, m2, v1, v2, v3, v4, v5, v6, v7, v8) \
	t1 = _mm512_mask_max_epu16(v1, m1, v1, dv1_swap);		\
	t2 = _mm512_mask_max_epu16(v2, m1, v2, dv2_swap);     \
	t3 = _mm512_mask_max_epu16(v3, m1, v3, dv3_swap);     \
	t4 = _mm512_mask_max_epu16(v4, m1, v4, dv4_swap);     \
	t5 = _mm512_mask_max_epu16(v5, m1, v5, dv5_swap);     \
	t6 = _mm512_mask_max_epu16(v6, m1, v6, dv6_swap);     \
	t7 = _mm512_mask_max_epu16(v7, m1, v7, dv7_swap);     \
	t8 = _mm512_mask_max_epu16(v8, m1, v8, dv8_swap);     \
	v1 = _mm512_mask_min_epu16(t1, m2, v1, dv1_swap);     \
	v2 = _mm512_mask_min_epu16(t2, m2, v2, dv2_swap);     \
	v3 = _mm512_mask_min_epu16(t3, m2, v3, dv3_swap);     \
	v4 = _mm512_mask_min_epu16(t4, m2, v4, dv4_swap);     \
	v5 = _mm512_mask_min_epu16(t5, m2, v5, dv5_swap);     \
	v6 = _mm512_mask_min_epu16(t6, m2, v6, dv6_swap);     \
	v7 = _mm512_mask_min_epu16(t7, m2, v7, dv7_swap);     \
	v8 = _mm512_mask_min_epu16(t8, m2, v8, dv8_swap); 
	
#define MINMAX16x8_alt1(m1, m2, v1, v2, v3, v4, v5, v6, v7, v8) \
	t1 = _mm512_mask_max_epu16(v1, m1, v1, dv1_swap);		\
	t2 = _mm512_mask_min_epu16(v2, m1, v2, dv2_swap);     \
	t3 = _mm512_mask_max_epu16(v3, m1, v3, dv3_swap);     \
	t4 = _mm512_mask_min_epu16(v4, m1, v4, dv4_swap);     \
	t5 = _mm512_mask_max_epu16(v5, m1, v5, dv5_swap);     \
	t6 = _mm512_mask_min_epu16(v6, m1, v6, dv6_swap);     \
	t7 = _mm512_mask_max_epu16(v7, m1, v7, dv7_swap);     \
	t8 = _mm512_mask_min_epu16(v8, m1, v8, dv8_swap);     \
	v1 = _mm512_mask_min_epu16(t1, m2, v1, dv1_swap);     \
	v2 = _mm512_mask_max_epu16(t2, m2, v2, dv2_swap);     \
	v3 = _mm512_mask_min_epu16(t3, m2, v3, dv3_swap);     \
	v4 = _mm512_mask_max_epu16(t4, m2, v4, dv4_swap);     \
	v5 = _mm512_mask_min_epu16(t5, m2, v5, dv5_swap);     \
	v6 = _mm512_mask_max_epu16(t6, m2, v6, dv6_swap);     \
	v7 = _mm512_mask_min_epu16(t7, m2, v7, dv7_swap);     \
	v8 = _mm512_mask_max_epu16(t8, m2, v8, dv8_swap); 
	
#define MINMAX16x8_alt2(m1, m2, v1, v2, v3, v4, v5, v6, v7, v8) \
	t1 = _mm512_mask_max_epu16(v1, m1, v1, dv1_swap);		\
	t2 = _mm512_mask_max_epu16(v2, m1, v2, dv2_swap);     \
	t3 = _mm512_mask_min_epu16(v3, m1, v3, dv3_swap);     \
	t4 = _mm512_mask_min_epu16(v4, m1, v4, dv4_swap);     \
	t5 = _mm512_mask_max_epu16(v5, m1, v5, dv5_swap);     \
	t6 = _mm512_mask_max_epu16(v6, m1, v6, dv6_swap);     \
	t7 = _mm512_mask_min_epu16(v7, m1, v7, dv7_swap);     \
	t8 = _mm512_mask_min_epu16(v8, m1, v8, dv8_swap);     \
	v1 = _mm512_mask_min_epu16(t1, m2, v1, dv1_swap);     \
	v2 = _mm512_mask_min_epu16(t2, m2, v2, dv2_swap);     \
	v3 = _mm512_mask_max_epu16(t3, m2, v3, dv3_swap);     \
	v4 = _mm512_mask_max_epu16(t4, m2, v4, dv4_swap);     \
	v5 = _mm512_mask_min_epu16(t5, m2, v5, dv5_swap);     \
	v6 = _mm512_mask_min_epu16(t6, m2, v6, dv6_swap);     \
	v7 = _mm512_mask_max_epu16(t7, m2, v7, dv7_swap);     \
	v8 = _mm512_mask_max_epu16(t8, m2, v8, dv8_swap); 
	
#define MINMAX16x8_alt4(m1, m2, v1, v2, v3, v4, v5, v6, v7, v8) \
	t1 = _mm512_mask_max_epu16(v1, m1, v1, dv1_swap);		\
	t2 = _mm512_mask_max_epu16(v2, m1, v2, dv2_swap);     \
	t3 = _mm512_mask_max_epu16(v3, m1, v3, dv3_swap);     \
	t4 = _mm512_mask_max_epu16(v4, m1, v4, dv4_swap);     \
	t5 = _mm512_mask_min_epu16(v5, m1, v5, dv5_swap);     \
	t6 = _mm512_mask_min_epu16(v6, m1, v6, dv6_swap);     \
	t7 = _mm512_mask_min_epu16(v7, m1, v7, dv7_swap);     \
	t8 = _mm512_mask_min_epu16(v8, m1, v8, dv8_swap);     \
	v1 = _mm512_mask_min_epu16(t1, m2, v1, dv1_swap);     \
	v2 = _mm512_mask_min_epu16(t2, m2, v2, dv2_swap);     \
	v3 = _mm512_mask_min_epu16(t3, m2, v3, dv3_swap);     \
	v4 = _mm512_mask_min_epu16(t4, m2, v4, dv4_swap);     \
	v5 = _mm512_mask_max_epu16(t5, m2, v5, dv5_swap);     \
	v6 = _mm512_mask_max_epu16(t6, m2, v6, dv6_swap);     \
	v7 = _mm512_mask_max_epu16(t7, m2, v7, dv7_swap);     \
	v8 = _mm512_mask_max_epu16(t8, m2, v8, dv8_swap); 
	
#define CMPGT(v1, v2, v3, v4, v5, v6, v7, v8) \
	m1 = _mm512_cmpgt_epu32_mask(v1, dv1_swap);  \
	m2 = _mm512_cmpgt_epu32_mask(v2, dv2_swap);  \
	m3 = _mm512_cmpgt_epu32_mask(v3, dv3_swap);  \
	m4 = _mm512_cmpgt_epu32_mask(v4, dv4_swap);  \
	m5 = _mm512_cmpgt_epu32_mask(v5, dv5_swap);  \
	m6 = _mm512_cmpgt_epu32_mask(v6, dv6_swap);  \
	m7 = _mm512_cmpgt_epu32_mask(v7, dv7_swap);  \
	m8 = _mm512_cmpgt_epu32_mask(v8, dv8_swap);
	
#define CMPLT(v1, v2, v3, v4, v5, v6, v7, v8) \
	m1 = _mm512_cmplt_epu32_mask(v1, dv1_swap);  \
	m2 = _mm512_cmplt_epu32_mask(v2, dv2_swap);  \
	m3 = _mm512_cmplt_epu32_mask(v3, dv3_swap);  \
	m4 = _mm512_cmplt_epu32_mask(v4, dv4_swap);  \
	m5 = _mm512_cmplt_epu32_mask(v5, dv5_swap);  \
	m6 = _mm512_cmplt_epu32_mask(v6, dv6_swap);  \
	m7 = _mm512_cmplt_epu32_mask(v7, dv7_swap);  \
	m8 = _mm512_cmplt_epu32_mask(v8, dv8_swap);
		
#define CMP(v1, v2, v3, v4, v5, v6, v7, v8, c1, c2, c3, c4, c5, c6, c7, c8) \
	m1 = _mm512_cmp_epu32_mask(v1, dv1_swap, c1);  \
	m2 = _mm512_cmp_epu32_mask(v2, dv2_swap, c2);  \
	m3 = _mm512_cmp_epu32_mask(v3, dv3_swap, c3);  \
	m4 = _mm512_cmp_epu32_mask(v4, dv4_swap, c4);  \
	m5 = _mm512_cmp_epu32_mask(v5, dv5_swap, c5);  \
	m6 = _mm512_cmp_epu32_mask(v6, dv6_swap, c6);  \
	m7 = _mm512_cmp_epu32_mask(v7, dv7_swap, c7);  \
	m8 = _mm512_cmp_epu32_mask(v8, dv8_swap, c8);
	
#define BLENDMASK(M, v1, v2, v3, v4, v5, v6, v7, v8) \
	v1 = _mm512_mask_blend_epi32(m1 ^ (M), v1, dv1_swap);  \
	v2 = _mm512_mask_blend_epi32(m2 ^ (M), v2, dv2_swap);  \
	v3 = _mm512_mask_blend_epi32(m3 ^ (M), v3, dv3_swap);  \
	v4 = _mm512_mask_blend_epi32(m4 ^ (M), v4, dv4_swap);  \
	v5 = _mm512_mask_blend_epi32(m5 ^ (M), v5, dv5_swap);  \
	v6 = _mm512_mask_blend_epi32(m6 ^ (M), v6, dv6_swap);  \
	v7 = _mm512_mask_blend_epi32(m7 ^ (M), v7, dv7_swap);  \
	v8 = _mm512_mask_blend_epi32(m8 ^ (M), v8, dv8_swap);

#define MINMAX32x4(m1, m2, v1, v2, v3, v4) \
	t1 = _mm512_mask_max_epu32(v1, m1, v1, dv1_swap);		\
	t2 = _mm512_mask_max_epu32(v2, m1, v2, dv2_swap);     \
	t3 = _mm512_mask_max_epu32(v3, m1, v3, dv3_swap);     \
	t4 = _mm512_mask_max_epu32(v4, m1, v4, dv4_swap);     \
	v1 = _mm512_mask_min_epu32(t1, m2, v1, dv1_swap);     \
	v2 = _mm512_mask_min_epu32(t2, m2, v2, dv2_swap);     \
	v3 = _mm512_mask_min_epu32(t3, m2, v3, dv3_swap);     \
	v4 = _mm512_mask_min_epu32(t4, m2, v4, dv4_swap);
	
#define MINMAX32x4_alt1(m1, m2, v1, v2, v3, v4) \
	t1 = _mm512_mask_max_epu32(v1, m1, v1, dv1_swap);		\
	t2 = _mm512_mask_min_epu32(v2, m1, v2, dv2_swap);     \
	t3 = _mm512_mask_max_epu32(v3, m1, v3, dv3_swap);     \
	t4 = _mm512_mask_min_epu32(v4, m1, v4, dv4_swap);     \
	v1 = _mm512_mask_min_epu32(t1, m2, v1, dv1_swap);     \
	v2 = _mm512_mask_max_epu32(t2, m2, v2, dv2_swap);     \
	v3 = _mm512_mask_min_epu32(t3, m2, v3, dv3_swap);     \
	v4 = _mm512_mask_max_epu32(t4, m2, v4, dv4_swap);
	
#define MINMAX32x4_alt2(m1, m2, v1, v2, v3, v4) \
	t1 = _mm512_mask_max_epu32(v1, m1, v1, dv1_swap);		\
	t2 = _mm512_mask_max_epu32(v2, m1, v2, dv2_swap);     \
	t3 = _mm512_mask_min_epu32(v3, m1, v3, dv3_swap);     \
	t4 = _mm512_mask_min_epu32(v4, m1, v4, dv4_swap);     \
	v1 = _mm512_mask_min_epu32(t1, m2, v1, dv1_swap);     \
	v2 = _mm512_mask_min_epu32(t2, m2, v2, dv2_swap);     \
	v3 = _mm512_mask_max_epu32(t3, m2, v3, dv3_swap);     \
	v4 = _mm512_mask_max_epu32(t4, m2, v4, dv4_swap);
	
#define MINMAX(m1, m2, v1, v2, v3, v4, v5, v6, v7, v8) \
	t1 = _mm512_mask_max_epu32(v1, m1, v1, dv1_swap);		\
	t2 = _mm512_mask_max_epu32(v2, m1, v2, dv2_swap);     \
	t3 = _mm512_mask_max_epu32(v3, m1, v3, dv3_swap);     \
	t4 = _mm512_mask_max_epu32(v4, m1, v4, dv4_swap);     \
	t5 = _mm512_mask_max_epu32(v5, m1, v5, dv5_swap);     \
	t6 = _mm512_mask_max_epu32(v6, m1, v6, dv6_swap);     \
	t7 = _mm512_mask_max_epu32(v7, m1, v7, dv7_swap);     \
	t8 = _mm512_mask_max_epu32(v8, m1, v8, dv8_swap);     \
	v1 = _mm512_mask_min_epu32(t1, m2, v1, dv1_swap);     \
	v2 = _mm512_mask_min_epu32(t2, m2, v2, dv2_swap);     \
	v3 = _mm512_mask_min_epu32(t3, m2, v3, dv3_swap);     \
	v4 = _mm512_mask_min_epu32(t4, m2, v4, dv4_swap);     \
	v5 = _mm512_mask_min_epu32(t5, m2, v5, dv5_swap);     \
	v6 = _mm512_mask_min_epu32(t6, m2, v6, dv6_swap);     \
	v7 = _mm512_mask_min_epu32(t7, m2, v7, dv7_swap);     \
	v8 = _mm512_mask_min_epu32(t8, m2, v8, dv8_swap); 
	
#define MINMAX_alt1(m1, m2, v1, v2, v3, v4, v5, v6, v7, v8) \
	t1 = _mm512_mask_max_epu32(v1, m1, v1, dv1_swap);		\
	t2 = _mm512_mask_min_epu32(v2, m1, v2, dv2_swap);     \
	t3 = _mm512_mask_max_epu32(v3, m1, v3, dv3_swap);     \
	t4 = _mm512_mask_min_epu32(v4, m1, v4, dv4_swap);     \
	t5 = _mm512_mask_max_epu32(v5, m1, v5, dv5_swap);     \
	t6 = _mm512_mask_min_epu32(v6, m1, v6, dv6_swap);     \
	t7 = _mm512_mask_max_epu32(v7, m1, v7, dv7_swap);     \
	t8 = _mm512_mask_min_epu32(v8, m1, v8, dv8_swap);     \
	v1 = _mm512_mask_min_epu32(t1, m2, v1, dv1_swap);     \
	v2 = _mm512_mask_max_epu32(t2, m2, v2, dv2_swap);     \
	v3 = _mm512_mask_min_epu32(t3, m2, v3, dv3_swap);     \
	v4 = _mm512_mask_max_epu32(t4, m2, v4, dv4_swap);     \
	v5 = _mm512_mask_min_epu32(t5, m2, v5, dv5_swap);     \
	v6 = _mm512_mask_max_epu32(t6, m2, v6, dv6_swap);     \
	v7 = _mm512_mask_min_epu32(t7, m2, v7, dv7_swap);     \
	v8 = _mm512_mask_max_epu32(t8, m2, v8, dv8_swap); 
	
#define MINMAX_alt2(m1, m2, v1, v2, v3, v4, v5, v6, v7, v8) \
	t1 = _mm512_mask_max_epu32(v1, m1, v1, dv1_swap);		\
	t2 = _mm512_mask_max_epu32(v2, m1, v2, dv2_swap);     \
	t3 = _mm512_mask_min_epu32(v3, m1, v3, dv3_swap);     \
	t4 = _mm512_mask_min_epu32(v4, m1, v4, dv4_swap);     \
	t5 = _mm512_mask_max_epu32(v5, m1, v5, dv5_swap);     \
	t6 = _mm512_mask_max_epu32(v6, m1, v6, dv6_swap);     \
	t7 = _mm512_mask_min_epu32(v7, m1, v7, dv7_swap);     \
	t8 = _mm512_mask_min_epu32(v8, m1, v8, dv8_swap);     \
	v1 = _mm512_mask_min_epu32(t1, m2, v1, dv1_swap);     \
	v2 = _mm512_mask_min_epu32(t2, m2, v2, dv2_swap);     \
	v3 = _mm512_mask_max_epu32(t3, m2, v3, dv3_swap);     \
	v4 = _mm512_mask_max_epu32(t4, m2, v4, dv4_swap);     \
	v5 = _mm512_mask_min_epu32(t5, m2, v5, dv5_swap);     \
	v6 = _mm512_mask_min_epu32(t6, m2, v6, dv6_swap);     \
	v7 = _mm512_mask_max_epu32(t7, m2, v7, dv7_swap);     \
	v8 = _mm512_mask_max_epu32(t8, m2, v8, dv8_swap); 
	
#define MINMAX_alt4(m1, m2, v1, v2, v3, v4, v5, v6, v7, v8) \
	t1 = _mm512_mask_max_epu32(v1, m1, v1, dv1_swap);		\
	t2 = _mm512_mask_max_epu32(v2, m1, v2, dv2_swap);     \
	t3 = _mm512_mask_max_epu32(v3, m1, v3, dv3_swap);     \
	t4 = _mm512_mask_max_epu32(v4, m1, v4, dv4_swap);     \
	t5 = _mm512_mask_min_epu32(v5, m1, v5, dv5_swap);     \
	t6 = _mm512_mask_min_epu32(v6, m1, v6, dv6_swap);     \
	t7 = _mm512_mask_min_epu32(v7, m1, v7, dv7_swap);     \
	t8 = _mm512_mask_min_epu32(v8, m1, v8, dv8_swap);     \
	v1 = _mm512_mask_min_epu32(t1, m2, v1, dv1_swap);     \
	v2 = _mm512_mask_min_epu32(t2, m2, v2, dv2_swap);     \
	v3 = _mm512_mask_min_epu32(t3, m2, v3, dv3_swap);     \
	v4 = _mm512_mask_min_epu32(t4, m2, v4, dv4_swap);     \
	v5 = _mm512_mask_max_epu32(t5, m2, v5, dv5_swap);     \
	v6 = _mm512_mask_max_epu32(t6, m2, v6, dv6_swap);     \
	v7 = _mm512_mask_max_epu32(t7, m2, v7, dv7_swap);     \
	v8 = _mm512_mask_max_epu32(t8, m2, v8, dv8_swap); 

#define CMPGT64(v1, v2, v3, v4, v5, v6, v7, v8) \
	m1 = _mm512_cmpgt_epu64_mask(v1, dv1_swap);  \
	m2 = _mm512_cmpgt_epu64_mask(v2, dv2_swap);  \
	m3 = _mm512_cmpgt_epu64_mask(v3, dv3_swap);  \
	m4 = _mm512_cmpgt_epu64_mask(v4, dv4_swap);  \
	m5 = _mm512_cmpgt_epu64_mask(v5, dv5_swap);  \
	m6 = _mm512_cmpgt_epu64_mask(v6, dv6_swap);  \
	m7 = _mm512_cmpgt_epu64_mask(v7, dv7_swap);  \
	m8 = _mm512_cmpgt_epu64_mask(v8, dv8_swap);
	
#define CMPLT64(v1, v2, v3, v4, v5, v6, v7, v8) \
	m1 = _mm512_cmplt_epu64_mask(v1, dv1_swap);  \
	m2 = _mm512_cmplt_epu64_mask(v2, dv2_swap);  \
	m3 = _mm512_cmplt_epu64_mask(v3, dv3_swap);  \
	m4 = _mm512_cmplt_epu64_mask(v4, dv4_swap);  \
	m5 = _mm512_cmplt_epu64_mask(v5, dv5_swap);  \
	m6 = _mm512_cmplt_epu64_mask(v6, dv6_swap);  \
	m7 = _mm512_cmplt_epu64_mask(v7, dv7_swap);  \
	m8 = _mm512_cmplt_epu64_mask(v8, dv8_swap);
		
#define CMP64(v1, v2, v3, v4, v5, v6, v7, v8, c1, c2, c3, c4, c5, c6, c7, c8) \
	m1 = _mm512_cmp_epu64_mask(v1, dv1_swap, c1);  \
	m2 = _mm512_cmp_epu64_mask(v2, dv2_swap, c2);  \
	m3 = _mm512_cmp_epu64_mask(v3, dv3_swap, c3);  \
	m4 = _mm512_cmp_epu64_mask(v4, dv4_swap, c4);  \
	m5 = _mm512_cmp_epu64_mask(v5, dv5_swap, c5);  \
	m6 = _mm512_cmp_epu64_mask(v6, dv6_swap, c6);  \
	m7 = _mm512_cmp_epu64_mask(v7, dv7_swap, c7);  \
	m8 = _mm512_cmp_epu64_mask(v8, dv8_swap, c8);
	
#define BLENDMASK64(M, v1, v2, v3, v4, v5, v6, v7, v8) \
	v1 = _mm512_mask_blend_epi64(m1 ^ (M), v1, dv1_swap);  \
	v2 = _mm512_mask_blend_epi64(m2 ^ (M), v2, dv2_swap);  \
	v3 = _mm512_mask_blend_epi64(m3 ^ (M), v3, dv3_swap);  \
	v4 = _mm512_mask_blend_epi64(m4 ^ (M), v4, dv4_swap);  \
	v5 = _mm512_mask_blend_epi64(m5 ^ (M), v5, dv5_swap);  \
	v6 = _mm512_mask_blend_epi64(m6 ^ (M), v6, dv6_swap);  \
	v7 = _mm512_mask_blend_epi64(m7 ^ (M), v7, dv7_swap);  \
	v8 = _mm512_mask_blend_epi64(m8 ^ (M), v8, dv8_swap);

#define MINMAX64(m1, m2, v1, v2, v3, v4, v5, v6, v7, v8) \
	t1 = _mm512_mask_max_epu64(v1, m1, v1, dv1_swap);		\
	t2 = _mm512_mask_max_epu64(v2, m1, v2, dv2_swap);     \
	t3 = _mm512_mask_max_epu64(v3, m1, v3, dv3_swap);     \
	t4 = _mm512_mask_max_epu64(v4, m1, v4, dv4_swap);     \
	t5 = _mm512_mask_max_epu64(v5, m1, v5, dv5_swap);     \
	t6 = _mm512_mask_max_epu64(v6, m1, v6, dv6_swap);     \
	t7 = _mm512_mask_max_epu64(v7, m1, v7, dv7_swap);     \
	t8 = _mm512_mask_max_epu64(v8, m1, v8, dv8_swap);     \
	v1 = _mm512_mask_min_epu64(t1, m2, v1, dv1_swap);     \
	v2 = _mm512_mask_min_epu64(t2, m2, v2, dv2_swap);     \
	v3 = _mm512_mask_min_epu64(t3, m2, v3, dv3_swap);     \
	v4 = _mm512_mask_min_epu64(t4, m2, v4, dv4_swap);     \
	v5 = _mm512_mask_min_epu64(t5, m2, v5, dv5_swap);     \
	v6 = _mm512_mask_min_epu64(t6, m2, v6, dv6_swap);     \
	v7 = _mm512_mask_min_epu64(t7, m2, v7, dv7_swap);     \
	v8 = _mm512_mask_min_epu64(t8, m2, v8, dv8_swap); 
	
#define MINMAX64_alt1(m1, m2, v1, v2, v3, v4, v5, v6, v7, v8) \
	t1 = _mm512_mask_max_epu64(v1, m1, v1, dv1_swap);		\
	t2 = _mm512_mask_min_epu64(v2, m1, v2, dv2_swap);     \
	t3 = _mm512_mask_max_epu64(v3, m1, v3, dv3_swap);     \
	t4 = _mm512_mask_min_epu64(v4, m1, v4, dv4_swap);     \
	t5 = _mm512_mask_max_epu64(v5, m1, v5, dv5_swap);     \
	t6 = _mm512_mask_min_epu64(v6, m1, v6, dv6_swap);     \
	t7 = _mm512_mask_max_epu64(v7, m1, v7, dv7_swap);     \
	t8 = _mm512_mask_min_epu64(v8, m1, v8, dv8_swap);     \
	v1 = _mm512_mask_min_epu64(t1, m2, v1, dv1_swap);     \
	v2 = _mm512_mask_max_epu64(t2, m2, v2, dv2_swap);     \
	v3 = _mm512_mask_min_epu64(t3, m2, v3, dv3_swap);     \
	v4 = _mm512_mask_max_epu64(t4, m2, v4, dv4_swap);     \
	v5 = _mm512_mask_min_epu64(t5, m2, v5, dv5_swap);     \
	v6 = _mm512_mask_max_epu64(t6, m2, v6, dv6_swap);     \
	v7 = _mm512_mask_min_epu64(t7, m2, v7, dv7_swap);     \
	v8 = _mm512_mask_max_epu64(t8, m2, v8, dv8_swap); 
	
#define MINMAX64_alt2(m1, m2, v1, v2, v3, v4, v5, v6, v7, v8) \
	t1 = _mm512_mask_max_epu64(v1, m1, v1, dv1_swap);		\
	t2 = _mm512_mask_max_epu64(v2, m1, v2, dv2_swap);     \
	t3 = _mm512_mask_min_epu64(v3, m1, v3, dv3_swap);     \
	t4 = _mm512_mask_min_epu64(v4, m1, v4, dv4_swap);     \
	t5 = _mm512_mask_max_epu64(v5, m1, v5, dv5_swap);     \
	t6 = _mm512_mask_max_epu64(v6, m1, v6, dv6_swap);     \
	t7 = _mm512_mask_min_epu64(v7, m1, v7, dv7_swap);     \
	t8 = _mm512_mask_min_epu64(v8, m1, v8, dv8_swap);     \
	v1 = _mm512_mask_min_epu64(t1, m2, v1, dv1_swap);     \
	v2 = _mm512_mask_min_epu64(t2, m2, v2, dv2_swap);     \
	v3 = _mm512_mask_max_epu64(t3, m2, v3, dv3_swap);     \
	v4 = _mm512_mask_max_epu64(t4, m2, v4, dv4_swap);     \
	v5 = _mm512_mask_min_epu64(t5, m2, v5, dv5_swap);     \
	v6 = _mm512_mask_min_epu64(t6, m2, v6, dv6_swap);     \
	v7 = _mm512_mask_max_epu64(t7, m2, v7, dv7_swap);     \
	v8 = _mm512_mask_max_epu64(t8, m2, v8, dv8_swap); 
	
#define MINMAX64_alt4(m1, m2, v1, v2, v3, v4, v5, v6, v7, v8) \
	t1 = _mm512_mask_max_epu64(v1, m1, v1, dv1_swap);		\
	t2 = _mm512_mask_max_epu64(v2, m1, v2, dv2_swap);     \
	t3 = _mm512_mask_max_epu64(v3, m1, v3, dv3_swap);     \
	t4 = _mm512_mask_max_epu64(v4, m1, v4, dv4_swap);     \
	t5 = _mm512_mask_min_epu64(v5, m1, v5, dv5_swap);     \
	t6 = _mm512_mask_min_epu64(v6, m1, v6, dv6_swap);     \
	t7 = _mm512_mask_min_epu64(v7, m1, v7, dv7_swap);     \
	t8 = _mm512_mask_min_epu64(v8, m1, v8, dv8_swap);     \
	v1 = _mm512_mask_min_epu64(t1, m2, v1, dv1_swap);     \
	v2 = _mm512_mask_min_epu64(t2, m2, v2, dv2_swap);     \
	v3 = _mm512_mask_min_epu64(t3, m2, v3, dv3_swap);     \
	v4 = _mm512_mask_min_epu64(t4, m2, v4, dv4_swap);     \
	v5 = _mm512_mask_max_epu64(t5, m2, v5, dv5_swap);     \
	v6 = _mm512_mask_max_epu64(t6, m2, v6, dv6_swap);     \
	v7 = _mm512_mask_max_epu64(t7, m2, v7, dv7_swap);     \
	v8 = _mm512_mask_max_epu64(t8, m2, v8, dv8_swap); 


#if !defined(__clang__)
#define _mm512_storeu_epi32 _mm512_storeu_si512
#define _mm512_loadu_epi32 _mm512_loadu_si512
#define _mm512_storeu_epi64 _mm512_storeu_si512
#define _mm512_loadu_epi64 _mm512_loadu_si512
#endif

void print128(__m512i i1, __m512i i2, __m512i i3, __m512i i4,
	__m512i i5, __m512i i6, __m512i i7, __m512i i8)
{
	uint32_t t[16];
	int i;
	
	_mm512_storeu_si512(t, i1);
	for (i = 0; i < 16; i++)
	{
		printf("%08u ", t[i]);
	}
	printf("\n");
	
	_mm512_storeu_si512(t, i2);
	for (i = 0; i < 16; i++)
	{
		printf("%08u ", t[i]);
	}
	printf("\n");
	
	_mm512_storeu_si512(t, i3);
	for (i = 0; i < 16; i++)
	{
		printf("%08u ", t[i]);
	}
	printf("\n");
	
	_mm512_storeu_si512(t, i4);
	for (i = 0; i < 16; i++)
	{
		printf("%08u ", t[i]);
	}
	printf("\n");
	
	_mm512_storeu_si512(t, i5);
	for (i = 0; i < 16; i++)
	{
		printf("%08u ", t[i]);
	}
	printf("\n");
	
	_mm512_storeu_si512(t, i6);
	for (i = 0; i < 16; i++)
	{
		printf("%08u ", t[i]);
	}
	printf("\n");
	
	_mm512_storeu_si512(t, i7);
	for (i = 0; i < 16; i++)
	{
		printf("%08u ", t[i]);
	}
	printf("\n");
	
	_mm512_storeu_si512(t, i8);
	for (i = 0; i < 16; i++)
	{
		printf("%08u ", t[i]);
	}
	printf("\n");
}

void bitonic_merge16_dir_64(uint16_t* data, int dir) 
{
	// merge 64 16-bit elements
	__m512i t1;
	__m512i t2;
	__m512i dv1;
	__m512i dv2;
	__m512i dv1_swap;
	__m512i dv2_swap;

	dv1 = _mm512_load_si512(data);
	dv2 = _mm512_load_si512(data + 32);

	// phase 5: merge (all same compare 'dir')
	if (dir == 1)
	{
		// distance 32 swaps
		t1 =  _mm512_max_epu16(dv1, dv2);
		dv2 = _mm512_min_epu16(dv1, dv2);
		dv1 = t1;
		
		// distance 16 swaps
		SWAP256x2(dv1, dv2);
		MINMAX16x2(0x0000ffff, 0xFFff0000, dv1, dv2);
		
		// distance 8 swaps
		SWAP128x2(dv1, dv2);
		MINMAX16x2(0x00ff00ff, 0xFF00ff00, dv1, dv2);

		// distance 4 swaps
		SWAP64x2(dv1, dv2);
		MINMAX16x2(0x0f0f0f0f, 0xF0F0f0f0, dv1, dv2);

		// distance 2 swaps
		SWAP32x2(dv1, dv2);
		MINMAX16x2(0x33333333, 0xCCCCcccc, dv1, dv2);

		// adjacent swaps
		SWAP16x2(dv1, dv2);
		MINMAX16x2(0x55555555, 0xAAAAaaaa, dv1, dv2);

	}
	else
	{	
		// distance 32 swaps
		t1 =  _mm512_min_epu16(dv1, dv2);
		dv2 = _mm512_max_epu16(dv1, dv2);
		dv1 = t1;
		
		// distance 16 swaps
		SWAP256x2(dv1, dv2);
		MINMAX16x2(0xFFff0000, 0x0000ffff, dv1, dv2);
		
		// distance 8 swaps
		SWAP128x2(dv1, dv2);
		MINMAX16x2(0xFF00ff00, 0x00ff00ff, dv1, dv2);

		// distance 4 swaps
		SWAP64x2(dv1, dv2);
		MINMAX16x2(0xF0F0f0f0, 0x0f0f0f0f, dv1, dv2);

		// distance 2 swaps
		SWAP32x2(dv1, dv2);
		MINMAX16x2(0xCCCCcccc, 0x33333333, dv1, dv2);

		// adjacent swaps
		SWAP16x2(dv1, dv2);
		MINMAX16x2(0xAAAAaaaa, 0x55555555, dv1, dv2);

	}

	_mm512_store_si512(data, dv1);
	_mm512_store_si512(data + 32, dv2);
	
	return;
}

void bitonic_sort16_dir_64(uint16_t* data, int dir) 
{
	// sort 128 16-bit elements
	__m512i t1;
	__m512i t2;
	__m512i dv1;
	__m512i dv2;
	__m512i dv1_swap;
	__m512i dv2_swap;

	dv1 = _mm512_load_si512(data);
	dv2 = _mm512_load_si512(data + 32);

	// phase 0: dist-2 alternating compares ('CC')

	// adjacent swaps (AA)
	SWAP16x2(dv1, dv2);
	MINMAX16x2(0x66666666, 0x99999999, dv1, dv2);
	
	// phase 1: dist-4 alternating compares ('F0')
	
	// distance 2 swaps (CC)
	SWAP32x2(dv1, dv2);
	MINMAX16x2(0x3C3C3C3C, 0xC3C3C3C3, dv1, dv2)

	// adjacent swaps (AA)
	SWAP16x2(dv1, dv2);
	MINMAX16x2(0x5A5A5A5A, 0xA5A5A5A5, dv1, dv2);

	// phase 2: dist-8 alternating compares ('FF00')

	// distance 4 swaps (F0F0)
	SWAP64x2(dv1, dv2);
	MINMAX16x2(0x0FF00FF0, 0xf00ff00f, dv1, dv2);

	// distance 2 swaps (CC)
	SWAP32x2(dv1, dv2);
	MINMAX16x2(0x33CC33CC, 0xcc33cc33, dv1, dv2);

	// adjacent swaps (AA)
	SWAP16x2(dv1, dv2);
	MINMAX16x2(0x55AA55AA, 0xaa55aa55, dv1, dv2);

	// phase 2: dist-16 alternating compares ('FFFF0000')
	
	// distance 8 swaps (FF00FF00)
	SWAP128x2(dv1, dv2);
	MINMAX16x2(0x00FFFF00, 0xFF0000FF, dv1, dv2);

	// distance 4 swaps (F0F0)
	SWAP64x2(dv1, dv2);
	MINMAX16x2(0x0F0FF0F0, 0xf0f00f0f, dv1, dv2);

	// distance 2 swaps (CC)
	SWAP32x2(dv1, dv2);
	MINMAX16x2(0x3333CCCC, 0xCCCC3333, dv1, dv2);

	// adjacent swaps (AA)
	SWAP16x2(dv1, dv2);
	MINMAX16x2(0x5555AAAA, 0xaaaa5555, dv1, dv2);

	// phase 3: dist-32 alternating compares (alternating gt/lt)

	// distance 16 swaps
	SWAP256x2(dv1, dv2);
	MINMAX16x2_alt1(0xFFff0000, 0x0000ffff, dv1, dv2);
	
	// distance 8 swaps
	SWAP128x2(dv1, dv2);
	MINMAX16x2_alt1(0xFF00ff00, 0x00ff00ff, dv1, dv2);

	// distance 4 swaps
	SWAP64x2(dv1, dv2);
	MINMAX16x2_alt1(0xF0F0f0f0, 0x0f0f0f0f, dv1, dv2);

	// distance 2 swaps
	SWAP32x2(dv1, dv2);
	MINMAX16x2_alt1(0xCCCCcccc, 0x33333333, dv1, dv2);

	// adjacent swaps
	SWAP16x2(dv1, dv2);
	MINMAX16x2_alt1(0xAAAAaaaa, 0x55555555, dv1, dv2);

	// phase 5: merge (all same compare 'dir')
	if (dir == 1)
	{
		// distance 32 swaps
		t1 =  _mm512_max_epu16(dv1, dv2);
		dv2 = _mm512_min_epu16(dv1, dv2);
		dv1 = t1;
		
		// distance 16 swaps
		SWAP256x2(dv1, dv2);
		MINMAX16x2(0x0000ffff, 0xFFff0000, dv1, dv2);
		
		// distance 8 swaps
		SWAP128x2(dv1, dv2);
		MINMAX16x2(0x00ff00ff, 0xFF00ff00, dv1, dv2);

		// distance 4 swaps
		SWAP64x2(dv1, dv2);
		MINMAX16x2(0x0f0f0f0f, 0xF0F0f0f0, dv1, dv2);

		// distance 2 swaps
		SWAP32x2(dv1, dv2);
		MINMAX16x2(0x33333333, 0xCCCCcccc, dv1, dv2);

		// adjacent swaps
		SWAP16x2(dv1, dv2);
		MINMAX16x2(0x55555555, 0xAAAAaaaa, dv1, dv2);
	}
	else
	{	
		// distance 32 swaps
		t1 =  _mm512_min_epu16(dv1, dv2);
		dv2 = _mm512_max_epu16(dv1, dv2);
		dv1 = t1;
		
		// distance 16 swaps
		SWAP256x2(dv1, dv2);
		MINMAX16x2(0xFFff0000, 0x0000ffff, dv1, dv2);
		
		// distance 8 swaps
		SWAP128x2(dv1, dv2);
		MINMAX16x2(0xFF00ff00, 0x00ff00ff, dv1, dv2);

		// distance 4 swaps
		SWAP64x2(dv1, dv2);
		MINMAX16x2(0xF0F0f0f0, 0x0f0f0f0f, dv1, dv2);

		// distance 2 swaps
		SWAP32x2(dv1, dv2);
		MINMAX16x2(0xCCCCcccc, 0x33333333, dv1, dv2);

		// adjacent swaps
		SWAP16x2(dv1, dv2);
		MINMAX16x2(0xAAAAaaaa, 0x55555555, dv1, dv2);
	}

	_mm512_store_si512(data, dv1);
	_mm512_store_si512(data + 32, dv2);
	
	return;
}

void bitonic_merge16_dir_128(uint16_t* data, int dir) 
{
	// merge 128 16-bit elements
	__m512i t1;
	__m512i t2;
	__m512i t3;
	__m512i t4;
	__m512i dv1;
	__m512i dv2;
	__m512i dv3;
	__m512i dv4;
	__m512i dv1_swap;
	__m512i dv2_swap;
	__m512i dv3_swap;
	__m512i dv4_swap;

	dv1 = _mm512_load_si512(data);
	dv2 = _mm512_load_si512(data + 32);
	dv3 = _mm512_load_si512(data + 64);
	dv4 = _mm512_load_si512(data + 96);

	// phase 5: merge (all same compare 'dir')
	if (dir == 1)
	{
		// dist-64 swaps
		t1 =  _mm512_max_epu16(dv1, dv3);
		t2 =  _mm512_max_epu16(dv2, dv4);
		dv3 = _mm512_min_epu16(dv1, dv3);
		dv4 = _mm512_min_epu16(dv2, dv4);
		dv1 = t1;
		dv2 = t2;
		
		// distance 32 swaps
		t1 =  _mm512_max_epu16(dv1, dv2);
		t2 =  _mm512_max_epu16(dv3, dv4);
		dv2 = _mm512_min_epu16(dv1, dv2);
		dv4 = _mm512_min_epu16(dv3, dv4);
		dv1 = t1;
		dv3 = t2;
		
		// distance 16 swaps
		SWAP256x4(dv1, dv2, dv3, dv4);
		MINMAX16x4(0x0000ffff, 0xFFff0000, dv1, dv2, dv3, dv4);
		
		// distance 8 swaps
		SWAP128x4(dv1, dv2, dv3, dv4);
		MINMAX16x4(0x00ff00ff, 0xFF00ff00, dv1, dv2, dv3, dv4);

		// distance 4 swaps
		SWAP64x4(dv1, dv2, dv3, dv4);
		MINMAX16x4(0x0f0f0f0f, 0xF0F0f0f0, dv1, dv2, dv3, dv4);

		// distance 2 swaps
		SWAP32x4(dv1, dv2, dv3, dv4);
		MINMAX16x4(0x33333333, 0xCCCCcccc, dv1, dv2, dv3, dv4);

		// adjacent swaps
		SWAP16x4(dv1, dv2, dv3, dv4);
		MINMAX16x4(0x55555555, 0xAAAAaaaa, dv1, dv2, dv3, dv4);

	}
	else
	{	
		// dist-64 swaps
		t1 =  _mm512_min_epu16(dv1, dv3);
		t2 =  _mm512_min_epu16(dv2, dv4);
		dv3 = _mm512_max_epu16(dv1, dv3);
		dv4 = _mm512_max_epu16(dv2, dv4);
		dv1 = t1;
		dv2 = t2;
		
		// distance 32 swaps
		t1 =  _mm512_min_epu16(dv1, dv2);
		t2 =  _mm512_min_epu16(dv3, dv4);
		dv2 = _mm512_max_epu16(dv1, dv2);
		dv4 = _mm512_max_epu16(dv3, dv4);
		dv1 = t1;
		dv3 = t2;
		
		// distance 16 swaps
		SWAP256x4(dv1, dv2, dv3, dv4);
		MINMAX16x4(0xFFff0000, 0x0000ffff, dv1, dv2, dv3, dv4);
		
		// distance 8 swaps
		SWAP128x4(dv1, dv2, dv3, dv4);
		MINMAX16x4(0xFF00ff00, 0x00ff00ff, dv1, dv2, dv3, dv4);

		// distance 4 swaps
		SWAP64x4(dv1, dv2, dv3, dv4);
		MINMAX16x4(0xF0F0f0f0, 0x0f0f0f0f, dv1, dv2, dv3, dv4);

		// distance 2 swaps
		SWAP32x4(dv1, dv2, dv3, dv4);
		MINMAX16x4(0xCCCCcccc, 0x33333333, dv1, dv2, dv3, dv4);

		// adjacent swaps
		SWAP16x4(dv1, dv2, dv3, dv4);
		MINMAX16x4(0xAAAAaaaa, 0x55555555, dv1, dv2, dv3, dv4);

	}

	_mm512_store_si512(data, dv1);
	_mm512_store_si512(data + 32, dv2);
	_mm512_store_si512(data + 64, dv3);
	_mm512_store_si512(data + 96, dv4);
	
	return;
}

void bitonic_sort16_dir_128(uint16_t* data, int dir) 
{
	// sort 128 16-bit elements
	__m512i t1;
	__m512i t2;
	__m512i t3;
	__m512i t4;
	__m512i dv1;
	__m512i dv2;
	__m512i dv3;
	__m512i dv4;
	__m512i dv1_swap;
	__m512i dv2_swap;
	__m512i dv3_swap;
	__m512i dv4_swap;

	dv1 = _mm512_load_si512(data);
	dv2 = _mm512_load_si512(data + 32);
	dv3 = _mm512_load_si512(data + 64);
	dv4 = _mm512_load_si512(data + 96);

	// phase 0: dist-2 alternating compares ('CC')

	// adjacent swaps (AA)
	SWAP16x4(dv1, dv2, dv3, dv4);
	MINMAX16x4(0x66666666, 0x99999999, dv1, dv2, dv3, dv4);
	
	// phase 1: dist-4 alternating compares ('F0')
	
	// distance 2 swaps (CC)
	SWAP32x4(dv1, dv2, dv3, dv4);
	MINMAX16x4(0x3C3C3C3C, 0xC3C3C3C3, dv1, dv2, dv3, dv4)

	// adjacent swaps (AA)
	SWAP16x4(dv1, dv2, dv3, dv4);
	MINMAX16x4(0x5A5A5A5A, 0xA5A5A5A5, dv1, dv2, dv3, dv4);

	// phase 2: dist-8 alternating compares ('FF00')

	// distance 4 swaps (F0F0)
	SWAP64x4(dv1, dv2, dv3, dv4);
	MINMAX16x4(0x0FF00FF0, 0xf00ff00f, dv1, dv2, dv3, dv4);

	// distance 2 swaps (CC)
	SWAP32x4(dv1, dv2, dv3, dv4);
	MINMAX16x4(0x33CC33CC, 0xcc33cc33, dv1, dv2, dv3, dv4);

	// adjacent swaps (AA)
	SWAP16x4(dv1, dv2, dv3, dv4);
	MINMAX16x4(0x55AA55AA, 0xaa55aa55, dv1, dv2, dv3, dv4);

	// phase 2: dist-16 alternating compares ('FFFF0000')
	
	// distance 8 swaps (FF00FF00)
	SWAP128x4(dv1, dv2, dv3, dv4);
	MINMAX16x4(0x00FFFF00, 0xFF0000FF, dv1, dv2, dv3, dv4);

	// distance 4 swaps (F0F0)
	SWAP64x4(dv1, dv2, dv3, dv4);
	MINMAX16x4(0x0F0FF0F0, 0xf0f00f0f, dv1, dv2, dv3, dv4);

	// distance 2 swaps (CC)
	SWAP32x4(dv1, dv2, dv3, dv4);
	MINMAX16x4(0x3333CCCC, 0xCCCC3333, dv1, dv2, dv3, dv4);

	// adjacent swaps (AA)
	SWAP16x4(dv1, dv2, dv3, dv4);
	MINMAX16x4(0x5555AAAA, 0xaaaa5555, dv1, dv2, dv3, dv4);

	// phase 3: dist-32 alternating compares (alternating gt/lt)

	// distance 16 swaps
	SWAP256x4(dv1, dv2, dv3, dv4);
	MINMAX16x4_alt1(0xFFff0000, 0x0000ffff, dv1, dv2, dv3, dv4);
	
	// distance 8 swaps
	SWAP128x4(dv1, dv2, dv3, dv4);
	MINMAX16x4_alt1(0xFF00ff00, 0x00ff00ff, dv1, dv2, dv3, dv4);

	// distance 4 swaps
	SWAP64x4(dv1, dv2, dv3, dv4);
	MINMAX16x4_alt1(0xF0F0f0f0, 0x0f0f0f0f, dv1, dv2, dv3, dv4);

	// distance 2 swaps
	SWAP32x4(dv1, dv2, dv3, dv4);
	MINMAX16x4_alt1(0xCCCCcccc, 0x33333333, dv1, dv2, dv3, dv4);

	// adjacent swaps
	SWAP16x4(dv1, dv2, dv3, dv4);
	MINMAX16x4_alt1(0xAAAAaaaa, 0x55555555, dv1, dv2, dv3, dv4);
	
	// phase 4: dist-64 alternating compares (alternating gtgt/ltlt)

	// distance 32 swaps - just compare the vecs
	t1 =  _mm512_min_epu16(dv1, dv2);
	t2 =  _mm512_max_epu16(dv3, dv4);
	dv2 = _mm512_max_epu16(dv1, dv2);
	dv4 = _mm512_min_epu16(dv3, dv4);
	dv1 = t1;
	dv3 = t2;
	
	// distance 16 swaps
	SWAP256x4(dv1, dv2, dv3, dv4);
	MINMAX16x4_alt2(0xFFff0000, 0x0000ffff, dv1, dv2, dv3, dv4);
	
	// distance 8 swaps
	SWAP128x4(dv1, dv2, dv3, dv4);
	MINMAX16x4_alt2(0xFF00ff00, 0x00ff00ff, dv1, dv2, dv3, dv4);

	// distance 4 swaps
	SWAP64x4(dv1, dv2, dv3, dv4);
	MINMAX16x4_alt2(0xF0F0f0f0, 0x0f0f0f0f, dv1, dv2, dv3, dv4);

	// distance 2 swaps
	SWAP32x4(dv1, dv2, dv3, dv4);
	MINMAX16x4_alt2(0xCCCCcccc, 0x33333333, dv1, dv2, dv3, dv4);

	// adjacent swaps
	SWAP16x4(dv1, dv2, dv3, dv4);
	MINMAX16x4_alt2(0xAAAAaaaa, 0x55555555, dv1, dv2, dv3, dv4);


	// phase 5: merge (all same compare 'dir')
	if (dir == 1)
	{
		// dist-64 swaps
		t1 =  _mm512_max_epu16(dv1, dv3);
		t2 =  _mm512_max_epu16(dv2, dv4);
		dv3 = _mm512_min_epu16(dv1, dv3);
		dv4 = _mm512_min_epu16(dv2, dv4);
		dv1 = t1;
		dv2 = t2;
		
		// distance 32 swaps
		t1 =  _mm512_max_epu16(dv1, dv2);
		t2 =  _mm512_max_epu16(dv3, dv4);
		dv2 = _mm512_min_epu16(dv1, dv2);
		dv4 = _mm512_min_epu16(dv3, dv4);
		dv1 = t1;
		dv3 = t2;
		
		// distance 16 swaps
		SWAP256x4(dv1, dv2, dv3, dv4);
		MINMAX16x4(0x0000ffff, 0xFFff0000, dv1, dv2, dv3, dv4);
		
		// distance 8 swaps
		SWAP128x4(dv1, dv2, dv3, dv4);
		MINMAX16x4(0x00ff00ff, 0xFF00ff00, dv1, dv2, dv3, dv4);

		// distance 4 swaps
		SWAP64x4(dv1, dv2, dv3, dv4);
		MINMAX16x4(0x0f0f0f0f, 0xF0F0f0f0, dv1, dv2, dv3, dv4);

		// distance 2 swaps
		SWAP32x4(dv1, dv2, dv3, dv4);
		MINMAX16x4(0x33333333, 0xCCCCcccc, dv1, dv2, dv3, dv4);

		// adjacent swaps
		SWAP16x4(dv1, dv2, dv3, dv4);
		MINMAX16x4(0x55555555, 0xAAAAaaaa, dv1, dv2, dv3, dv4);

	}
	else
	{	
		// dist-64 swaps
		t1 =  _mm512_min_epu16(dv1, dv3);
		t2 =  _mm512_min_epu16(dv2, dv4);
		dv3 = _mm512_max_epu16(dv1, dv3);
		dv4 = _mm512_max_epu16(dv2, dv4);
		dv1 = t1;
		dv2 = t2;
		
		// distance 32 swaps
		t1 =  _mm512_min_epu16(dv1, dv2);
		t2 =  _mm512_min_epu16(dv3, dv4);
		dv2 = _mm512_max_epu16(dv1, dv2);
		dv4 = _mm512_max_epu16(dv3, dv4);
		dv1 = t1;
		dv3 = t2;
		
		// distance 16 swaps
		SWAP256x4(dv1, dv2, dv3, dv4);
		MINMAX16x4(0xFFff0000, 0x0000ffff, dv1, dv2, dv3, dv4);
		
		// distance 8 swaps
		SWAP128x4(dv1, dv2, dv3, dv4);
		MINMAX16x4(0xFF00ff00, 0x00ff00ff, dv1, dv2, dv3, dv4);

		// distance 4 swaps
		SWAP64x4(dv1, dv2, dv3, dv4);
		MINMAX16x4(0xF0F0f0f0, 0x0f0f0f0f, dv1, dv2, dv3, dv4);

		// distance 2 swaps
		SWAP32x4(dv1, dv2, dv3, dv4);
		MINMAX16x4(0xCCCCcccc, 0x33333333, dv1, dv2, dv3, dv4);

		// adjacent swaps
		SWAP16x4(dv1, dv2, dv3, dv4);
		MINMAX16x4(0xAAAAaaaa, 0x55555555, dv1, dv2, dv3, dv4);

	}

	_mm512_store_si512(data, dv1);
	_mm512_store_si512(data + 32, dv2);
	_mm512_store_si512(data + 64, dv3);
	_mm512_store_si512(data + 96, dv4);
	
	return;
}

void bitonic_merge16_dir_256(uint16_t* data, int dir) 
{
	// merge phase: 256 16-bit elements
	int i, j;

	__m512i t1;
	__m512i t2;
	__m512i t3;
	__m512i t4;
	__m512i t5;
	__m512i t6;
	__m512i t7;
	__m512i t8;
	__m512i dv1;
	__m512i dv2;
	__m512i dv3;
	__m512i dv4;
	__m512i dv5;
	__m512i dv6;
	__m512i dv7;
	__m512i dv8;
	__m512i dv1_swap;
	__m512i dv2_swap;
	__m512i dv3_swap;
	__m512i dv4_swap;
	__m512i dv5_swap;
	__m512i dv6_swap;
	__m512i dv7_swap;
	__m512i dv8_swap;
	__mmask32 m1;
	__mmask32 m2;
	__mmask32 m3;
	__mmask32 m4;
	__mmask32 m5;
	__mmask32 m6;
	__mmask32 m7;
	__mmask32 m8;

	dv1 = _mm512_load_si512(data);
	dv2 = _mm512_load_si512(data + 32);
	dv3 = _mm512_load_si512(data + 64);
	dv4 = _mm512_load_si512(data + 96);
	dv5 = _mm512_load_si512(data + 128);
	dv6 = _mm512_load_si512(data + 160);
	dv7 = _mm512_load_si512(data + 192);
	dv8 = _mm512_load_si512(data + 224);

	// phase 6: merge (all same compare 'dir')
	if (dir == 1)
	{
		// distance 128 swaps - just compare the vecs
		t1 =  _mm512_max_epu16(dv1, dv5);
		t2 =  _mm512_max_epu16(dv2, dv6);
		t3 =  _mm512_max_epu16(dv3, dv7);
		t4 =  _mm512_max_epu16(dv4, dv8);
		dv5 = _mm512_min_epu16(dv1, dv5);
		dv6 = _mm512_min_epu16(dv2, dv6);
		dv7 = _mm512_min_epu16(dv3, dv7);
		dv8 = _mm512_min_epu16(dv4, dv8);
		dv1 = t1;
		dv2 = t2;
		dv3 = t3;
		dv4 = t4;
		
		// distance 64 swaps - just compare the vecs
		t1 =  _mm512_max_epu16(dv1, dv3);
		t2 =  _mm512_max_epu16(dv2, dv4);
		t3 =  _mm512_max_epu16(dv5, dv7);
		t4 =  _mm512_max_epu16(dv6, dv8);
		dv3 = _mm512_min_epu16(dv1, dv3);
		dv4 = _mm512_min_epu16(dv2, dv4);
		dv7 = _mm512_min_epu16(dv5, dv7);
		dv8 = _mm512_min_epu16(dv6, dv8);
		dv1 = t1;
		dv2 = t2;
		dv5 = t3;
		dv6 = t4;
		
		// distance 32 swaps - just compare the vecs
		t1 =  _mm512_max_epu16(dv1, dv2);
		t2 =  _mm512_max_epu16(dv3, dv4);
		t3 =  _mm512_max_epu16(dv5, dv6);
		t4 =  _mm512_max_epu16(dv7, dv8);
		dv2 = _mm512_min_epu16(dv1, dv2);
		dv4 = _mm512_min_epu16(dv3, dv4);
		dv6 = _mm512_min_epu16(dv5, dv6);
		dv8 = _mm512_min_epu16(dv7, dv8);
		dv1 = t1;
		dv3 = t2;
		dv5 = t3;
		dv7 = t4;

		// distance 16 swaps
		SWAP256x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX16x8(0x0000ffff, 0xFFff0000, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		// distance 8 swaps
		SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX16x8(0x00ff00ff, 0xFF00ff00, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

		// distance 4 swaps
		SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX16x8(0x0f0f0f0f, 0xF0F0f0f0, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

		// distance 2 swaps
		SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX16x8(0x33333333, 0xCCCCcccc, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

		// adjacent swaps
		SWAP16x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX16x8(0x55555555, 0xAAAAaaaa, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	}
	else
	{	
		// distance 128 swaps - just compare the vecs
		t1 =  _mm512_min_epu16(dv1, dv5);
		t2 =  _mm512_min_epu16(dv2, dv6);
		t3 =  _mm512_min_epu16(dv3, dv7);
		t4 =  _mm512_min_epu16(dv4, dv8);
		dv5 = _mm512_max_epu16(dv1, dv5);
		dv6 = _mm512_max_epu16(dv2, dv6);
		dv7 = _mm512_max_epu16(dv3, dv7);
		dv8 = _mm512_max_epu16(dv4, dv8);
		dv1 = t1;
		dv2 = t2;
		dv3 = t3;
		dv4 = t4;
		
		// distance 64 swaps - just compare the vecs
		t1 =  _mm512_min_epu16(dv1, dv3);
		t2 =  _mm512_min_epu16(dv2, dv4);
		t3 =  _mm512_min_epu16(dv5, dv7);
		t4 =  _mm512_min_epu16(dv6, dv8);
		dv3 = _mm512_max_epu16(dv1, dv3);
		dv4 = _mm512_max_epu16(dv2, dv4);
		dv7 = _mm512_max_epu16(dv5, dv7);
		dv8 = _mm512_max_epu16(dv6, dv8);
		dv1 = t1;
		dv2 = t2;
		dv5 = t3;
		dv6 = t4;
		
		// distance 32 swaps - just compare the vecs
		t1 =  _mm512_min_epu16(dv1, dv2);
		t2 =  _mm512_min_epu16(dv3, dv4);
		t3 =  _mm512_min_epu16(dv5, dv6);
		t4 =  _mm512_min_epu16(dv7, dv8);
		dv2 = _mm512_max_epu16(dv1, dv2);
		dv4 = _mm512_max_epu16(dv3, dv4);
		dv6 = _mm512_max_epu16(dv5, dv6);
		dv8 = _mm512_max_epu16(dv7, dv8);
		dv1 = t1;
		dv3 = t2;
		dv5 = t3;
		dv7 = t4;
		
		// distance 16 swaps
		SWAP256x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX16x8(0xFFff0000, 0x0000ffff, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		// distance 8 swaps
		SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX16x8(0xFF00ff00, 0x00ff00ff, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

		// distance 4 swaps
		SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX16x8(0xF0F0f0f0, 0x0f0f0f0f, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

		// distance 2 swaps
		SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX16x8(0xCCCCcccc, 0x33333333, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

		// adjacent swaps
		SWAP16x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX16x8(0xAAAAaaaa, 0x55555555, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	}

	_mm512_store_si512(data + 0  , dv1);
	_mm512_store_si512(data + 32 , dv2);
	_mm512_store_si512(data + 64 , dv3);
	_mm512_store_si512(data + 96 , dv4);
	_mm512_store_si512(data + 128, dv5);
	_mm512_store_si512(data + 160, dv6);
	_mm512_store_si512(data + 192, dv7);
	_mm512_store_si512(data + 224, dv8);
	
	return;
}

void bitonic_sort16_dir_256(uint16_t* data, int dir) 
{
	// sort 256 16-bit elements
	__m512i t1;
	__m512i t2;
	__m512i t3;
	__m512i t4;
	__m512i t5;
	__m512i t6;
	__m512i t7;
	__m512i t8;
	__m512i dv1;
	__m512i dv2;
	__m512i dv3;
	__m512i dv4;
	__m512i dv5;
	__m512i dv6;
	__m512i dv7;
	__m512i dv8;
	__m512i dv1_swap;
	__m512i dv2_swap;
	__m512i dv3_swap;
	__m512i dv4_swap;
	__m512i dv5_swap;
	__m512i dv6_swap;
	__m512i dv7_swap;
	__m512i dv8_swap;
	__mmask32 m1;
	__mmask32 m2;

	dv1 = _mm512_load_si512(data);
	dv2 = _mm512_load_si512(data + 32);
	dv3 = _mm512_load_si512(data + 64);
	dv4 = _mm512_load_si512(data + 96);
	dv5 = _mm512_load_si512(data + 128);
	dv6 = _mm512_load_si512(data + 160);
	dv7 = _mm512_load_si512(data + 192);
	dv8 = _mm512_load_si512(data + 224);

	// phase 0: dist-2 alternating compares ('CC')
	
	// adjacent swaps
	SWAP16x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8(0x66666666, 0x99999999, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	// phase 1: dist-4 alternating compares ('F0')
	
	// distance 2 swaps
	SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8(0x3C3C3c3c, 0xC3C3c3c3, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8)

	// adjacent swaps
	SWAP16x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8(0x5A5A5a5a, 0xA5A5a5a5, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// phase 2: dist-8 alternating compares ('FF00')

	// distance 4 swaps
	SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8(0x0FF00ff0, 0xf00ff00f, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// distance 2 swaps
	SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8(0x33CC33cc, 0xcc33cc33, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// adjacent swaps
	SWAP16x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8(0x55AA55aa, 0xaa55aa55, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// phase 3: dist-16 alternating compares
	
	SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8(0x00FFFF00, 0xFF0000FF, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// distance 4 swaps (F0F0)
	SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8(0x0F0FF0F0, 0xf0f00f0f, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// distance 2 swaps (CC)
	SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8(0x3333CCCC, 0xCCCC3333, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// adjacent swaps (AA)
	SWAP16x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8(0x5555AAAA, 0xaaaa5555, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// phase 4: dist-32 alternating compares

	// distance 16 swaps
	SWAP256x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8_alt1(0xFFff0000, 0x0000ffff, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	// distance 8 swaps
	SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8_alt1(0xFF00ff00, 0x00ff00ff, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// distance 4 swaps
	SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8_alt1(0xF0F0f0f0, 0x0f0f0f0f, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// distance 2 swaps
	SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8_alt1(0xCCCCcccc, 0x33333333, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// adjacent swaps
	SWAP16x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8_alt1(0xAAAAaaaa, 0x55555555, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	// phase 5: dist-64 alternating compares

	// distance 32 swaps - just compare the vecs
	t1 =  _mm512_min_epu16(dv1, dv2);
	t2 =  _mm512_max_epu16(dv3, dv4);
	t3 =  _mm512_min_epu16(dv5, dv6);
	t4 =  _mm512_max_epu16(dv7, dv8);
	dv2 = _mm512_max_epu16(dv1, dv2);
	dv4 = _mm512_min_epu16(dv3, dv4);
	dv6 = _mm512_max_epu16(dv5, dv6);
	dv8 = _mm512_min_epu16(dv7, dv8);
	dv1 = t1;
	dv3 = t2;
	dv5 = t3;
	dv7 = t4;
	
	// distance 16 swaps
	SWAP256x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8_alt2(0xFFff0000, 0x0000ffff, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	// distance 8 swaps
	SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8_alt2(0xFF00ff00, 0x00ff00ff, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// distance 4 swaps
	SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8_alt2(0xF0F0f0f0, 0x0f0f0f0f, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// distance 2 swaps
	SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8_alt2(0xCCCCcccc, 0x33333333, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// adjacent swaps
	SWAP16x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8_alt2(0xAAAAaaaa, 0x55555555, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	// phase 6: dist-128 alternating compares

	// distance 64 swaps - just compare the vecs
	t1 =  _mm512_min_epu16(dv1, dv3);
	t2 =  _mm512_min_epu16(dv2, dv4);
	t3 =  _mm512_max_epu16(dv5, dv7);
	t4 =  _mm512_max_epu16(dv6, dv8);
	dv3 = _mm512_max_epu16(dv1, dv3);
	dv4 = _mm512_max_epu16(dv2, dv4);
	dv7 = _mm512_min_epu16(dv5, dv7);
	dv8 = _mm512_min_epu16(dv6, dv8);
	dv1 = t1;
	dv2 = t2;
	dv5 = t3;
	dv6 = t4;
	
	// distance 32 swaps - just compare the vecs
	t1 =  _mm512_min_epu16(dv1, dv2);
	t2 =  _mm512_min_epu16(dv3, dv4);
	t3 =  _mm512_max_epu16(dv5, dv6);
	t4 =  _mm512_max_epu16(dv7, dv8);
	dv2 = _mm512_max_epu16(dv1, dv2);
	dv4 = _mm512_max_epu16(dv3, dv4);
	dv6 = _mm512_min_epu16(dv5, dv6);
	dv8 = _mm512_min_epu16(dv7, dv8);
	dv1 = t1;
	dv3 = t2;
	dv5 = t3;
	dv7 = t4;
	
	// distance 16 swaps
	SWAP256x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8_alt4(0xFFff0000, 0x0000ffff, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	// distance 8 swaps
	SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8_alt4(0xFF00ff00, 0x00ff00ff, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// distance 4 swaps
	SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8_alt4(0xF0F0f0f0, 0x0f0f0f0f, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// distance 2 swaps
	SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8_alt4(0xCCCCcccc, 0x33333333, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// adjacent swaps
	SWAP16x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX16x8_alt4(0xAAAAaaaa, 0x55555555, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);


	// phase 7: merge (all same compare 'dir')
	if (dir == 1)
	{
		// distance 128 swaps - just compare the vecs
		t1 =  _mm512_max_epu16(dv1, dv5);
		t2 =  _mm512_max_epu16(dv2, dv6);
		t3 =  _mm512_max_epu16(dv3, dv7);
		t4 =  _mm512_max_epu16(dv4, dv8);
		dv5 = _mm512_min_epu16(dv1, dv5);
		dv6 = _mm512_min_epu16(dv2, dv6);
		dv7 = _mm512_min_epu16(dv3, dv7);
		dv8 = _mm512_min_epu16(dv4, dv8);
		dv1 = t1;
		dv2 = t2;
		dv3 = t3;
		dv4 = t4;
		
		// distance 64 swaps - just compare the vecs
		t1 =  _mm512_max_epu16(dv1, dv3);
		t2 =  _mm512_max_epu16(dv2, dv4);
		t3 =  _mm512_max_epu16(dv5, dv7);
		t4 =  _mm512_max_epu16(dv6, dv8);
		dv3 = _mm512_min_epu16(dv1, dv3);
		dv4 = _mm512_min_epu16(dv2, dv4);
		dv7 = _mm512_min_epu16(dv5, dv7);
		dv8 = _mm512_min_epu16(dv6, dv8);
		dv1 = t1;
		dv2 = t2;
		dv5 = t3;
		dv6 = t4;
		
		// distance 32 swaps - just compare the vecs
		t1 =  _mm512_max_epu16(dv1, dv2);
		t2 =  _mm512_max_epu16(dv3, dv4);
		t3 =  _mm512_max_epu16(dv5, dv6);
		t4 =  _mm512_max_epu16(dv7, dv8);
		dv2 = _mm512_min_epu16(dv1, dv2);
		dv4 = _mm512_min_epu16(dv3, dv4);
		dv6 = _mm512_min_epu16(dv5, dv6);
		dv8 = _mm512_min_epu16(dv7, dv8);
		dv1 = t1;
		dv3 = t2;
		dv5 = t3;
		dv7 = t4;

		// distance 16 swaps
		SWAP256x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX16x8(0x0000ffff, 0xFFff0000, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		// distance 8 swaps
		SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX16x8(0x00ff00ff, 0xFF00ff00, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

		// distance 4 swaps
		SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX16x8(0x0f0f0f0f, 0xF0F0f0f0, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

		// distance 2 swaps
		SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX16x8(0x33333333, 0xCCCCcccc, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

		// adjacent swaps
		SWAP16x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX16x8(0x55555555, 0xAAAAaaaa, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	}
	else
	{	
		// distance 128 swaps - just compare the vecs
		t1 =  _mm512_min_epu16(dv1, dv5);
		t2 =  _mm512_min_epu16(dv2, dv6);
		t3 =  _mm512_min_epu16(dv3, dv7);
		t4 =  _mm512_min_epu16(dv4, dv8);
		dv5 = _mm512_max_epu16(dv1, dv5);
		dv6 = _mm512_max_epu16(dv2, dv6);
		dv7 = _mm512_max_epu16(dv3, dv7);
		dv8 = _mm512_max_epu16(dv4, dv8);
		dv1 = t1;
		dv2 = t2;
		dv3 = t3;
		dv4 = t4;
		
		// distance 64 swaps - just compare the vecs
		t1 =  _mm512_min_epu16(dv1, dv3);
		t2 =  _mm512_min_epu16(dv2, dv4);
		t3 =  _mm512_min_epu16(dv5, dv7);
		t4 =  _mm512_min_epu16(dv6, dv8);
		dv3 = _mm512_max_epu16(dv1, dv3);
		dv4 = _mm512_max_epu16(dv2, dv4);
		dv7 = _mm512_max_epu16(dv5, dv7);
		dv8 = _mm512_max_epu16(dv6, dv8);
		dv1 = t1;
		dv2 = t2;
		dv5 = t3;
		dv6 = t4;
		
		// distance 32 swaps - just compare the vecs
		t1 =  _mm512_min_epu16(dv1, dv2);
		t2 =  _mm512_min_epu16(dv3, dv4);
		t3 =  _mm512_min_epu16(dv5, dv6);
		t4 =  _mm512_min_epu16(dv7, dv8);
		dv2 = _mm512_max_epu16(dv1, dv2);
		dv4 = _mm512_max_epu16(dv3, dv4);
		dv6 = _mm512_max_epu16(dv5, dv6);
		dv8 = _mm512_max_epu16(dv7, dv8);
		dv1 = t1;
		dv3 = t2;
		dv5 = t3;
		dv7 = t4;
		
		// distance 16 swaps
		SWAP256x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX16x8(0xFFff0000, 0x0000ffff, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		// distance 8 swaps
		SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX16x8(0xFF00ff00, 0x00ff00ff, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

		// distance 4 swaps
		SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX16x8(0xF0F0f0f0, 0x0f0f0f0f, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

		// distance 2 swaps
		SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX16x8(0xCCCCcccc, 0x33333333, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

		// adjacent swaps
		SWAP16x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX16x8(0xAAAAaaaa, 0x55555555, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	}

	_mm512_store_si512(data + 0  , dv1);
	_mm512_store_si512(data + 32 , dv2);
	_mm512_store_si512(data + 64 , dv3);
	_mm512_store_si512(data + 96 , dv4);
	_mm512_store_si512(data + 128, dv5);
	_mm512_store_si512(data + 160, dv6);
	_mm512_store_si512(data + 192, dv7);
	_mm512_store_si512(data + 224, dv8);
	
	return;
}

void bitonic_merge32_dir_64(uint32_t* data, int dir) 
{
	// perform the final merge phase on 64 32-bit elements
	int i, j;

	__m512i t1;
	__m512i t2;
	__m512i t3;
	__m512i t4;
	__m512i dv1;
	__m512i dv2;
	__m512i dv3;
	__m512i dv4;
	__m512i dv1_swap;
	__m512i dv2_swap;
	__m512i dv3_swap;
	__m512i dv4_swap;

	dv1 = _mm512_load_si512(data);
	dv2 = _mm512_load_si512(data + 16);
	dv3 = _mm512_load_si512(data + 32);
	dv4 = _mm512_load_si512(data + 48);

	if (dir == 1)
	{	
		// distance 32 swaps - just compare the vecs
		t1 =  _mm512_max_epu32(dv1, dv3);
		t2 =  _mm512_max_epu32(dv2, dv4);
		dv3 = _mm512_min_epu32(dv1, dv3);
		dv4 = _mm512_min_epu32(dv2, dv4);
		dv1 = t1;
		dv2 = t2;

		// distance 16 swaps - just compare the vecs
		t1 =  _mm512_max_epu32(dv1, dv2);
		t2 =  _mm512_max_epu32(dv3, dv4);
		dv2 = _mm512_min_epu32(dv1, dv2);
		dv4 = _mm512_min_epu32(dv3, dv4);
		dv1 = t1;
		dv3 = t2;
		
		// distance 8 swaps - 256 bits
		SWAP256x4(dv1, dv2, dv3, dv4);
		MINMAX32x4(0x00ff, 0xFF00, dv1, dv2, dv3, dv4);
		
		// distance 4 swaps
		SWAP128x4(dv1, dv2, dv3, dv4);
		MINMAX32x4(0x0f0f, 0xF0F0, dv1, dv2, dv3, dv4);
		
		// distance 2 swaps
		SWAP64x4(dv1, dv2, dv3, dv4);
		MINMAX32x4(0x3333, 0xCCCC, dv1, dv2, dv3, dv4);
		
		// adjacent swaps
		SWAP32x4(dv1, dv2, dv3, dv4);
		MINMAX32x4(0x5555, 0xaaaa, dv1, dv2, dv3, dv4);

	}
	else
	{		
		// distance 32 swaps - just compare the vecs
		t1 =  _mm512_min_epu32(dv1, dv3);
		t2 =  _mm512_min_epu32(dv2, dv4);
		dv3 = _mm512_max_epu32(dv1, dv3);
		dv4 = _mm512_max_epu32(dv2, dv4);
		dv1 = t1;
		dv2 = t2;

		// distance 16 swaps - just compare the vecs
		t1 =  _mm512_min_epu32(dv1, dv2);
		t2 =  _mm512_min_epu32(dv3, dv4);
		dv2 = _mm512_max_epu32(dv1, dv2);
		dv4 = _mm512_max_epu32(dv3, dv4);
		dv1 = t1;
		dv3 = t2;

		// distance 8 swaps - 256 bits
		SWAP256x4(dv1, dv2, dv3, dv4);
		MINMAX32x4(0xFF00, 0x00ff, dv1, dv2, dv3, dv4);
		
		// distance 4 swaps
		SWAP128x4(dv1, dv2, dv3, dv4);
		MINMAX32x4(0xf0f0, 0x0f0f, dv1, dv2, dv3, dv4);

		// distance 2 swaps
		SWAP64x4(dv1, dv2, dv3, dv4);
		MINMAX32x4(0xCCCC, 0x3333, dv1, dv2, dv3, dv4);
		
		// adjacent swaps
		SWAP32x4(dv1, dv2, dv3, dv4);
		MINMAX32x4(0xaaaa, 0x5555, dv1, dv2, dv3, dv4);

	}

	_mm512_store_si512(data, dv1);
	_mm512_store_si512(data + 16, dv2);
	_mm512_store_si512(data + 32, dv3);
	_mm512_store_si512(data + 48, dv4);

	return;
}

void bitonic_sort32_dir_64(uint32_t* data, int dir) 
{
	// sort 64 32-bit elements
	__m512i t1;
	__m512i t2;
	__m512i t3;
	__m512i t4;
	__m512i dv1;
	__m512i dv2;
	__m512i dv3;
	__m512i dv4;
	__m512i dv1_swap;
	__m512i dv2_swap;
	__m512i dv3_swap;
	__m512i dv4_swap;

	dv1 = _mm512_load_si512(data);
	dv2 = _mm512_load_si512(data + 16);
	dv3 = _mm512_load_si512(data + 32);
	dv4 = _mm512_load_si512(data + 48);

	// phase 0: dist-2 alternating compares ('CC')
	
	// adjacent swaps, alternating compares
	SWAP32x4(dv1, dv2, dv3, dv4);
	MINMAX32x4(0x6666, 0x9999, dv1, dv2, dv3, dv4);
	
	// phase 1: dist-4 alternating compares ('F0')
	
	// distance 2 swaps
	SWAP64x4(dv1, dv2, dv3, dv4);
	MINMAX32x4(0x3C3C, 0xC3C3, dv1, dv2, dv3, dv4)

	// adjacent swaps
	SWAP32x4(dv1, dv2, dv3, dv4);
	MINMAX32x4(0x5A5A, 0xA5A5, dv1, dv2, dv3, dv4);

	// phase 2: dist-8 alternating compares ('FF00')

	// distance 4 swaps
	SWAP128x4(dv1, dv2, dv3, dv4);
	MINMAX32x4(0x0FF0, 0xf00f, dv1, dv2, dv3, dv4);

	// distance 2 swaps
	SWAP64x4(dv1, dv2, dv3, dv4);
	MINMAX32x4(0x33CC, 0xcc33, dv1, dv2, dv3, dv4);

	// adjacent swaps
	SWAP32x4(dv1, dv2, dv3, dv4);
	MINMAX32x4(0x55AA, 0xaa55, dv1, dv2, dv3, dv4);

	// phase 2: dist-16 alternating compares (alternating gt/lt)
	
	// distance 8 swaps (256 bits)
	SWAP256x4(dv1, dv2, dv3, dv4);
	MINMAX32x4_alt1(0xFF00, 0x00ff, dv1, dv2, dv3, dv4);

	// distance 4 swaps
	SWAP128x4(dv1, dv2, dv3, dv4);
	MINMAX32x4_alt1(0xF0F0, 0x0f0f, dv1, dv2, dv3, dv4);

	// distance 2 swaps
	SWAP64x4(dv1, dv2, dv3, dv4);
	MINMAX32x4_alt1(0xCCCC, 0x3333, dv1, dv2, dv3, dv4);

	// adjacent swaps
	SWAP32x4(dv1, dv2, dv3, dv4);
	MINMAX32x4_alt1(0xAAAA, 0x5555, dv1, dv2, dv3, dv4);

	// phase 3: dist-32 alternating compares (alternating gtgt/ltlt)

	// distance 16 swaps - just compare the vecs
	t1 =  _mm512_min_epu32(dv1, dv2);
	t2 =  _mm512_max_epu32(dv3, dv4);
	dv2 = _mm512_max_epu32(dv1, dv2);
	dv4 = _mm512_min_epu32(dv3, dv4);
	dv1 = t1;
	dv3 = t2;

	// distance 8 swaps (256 bits)
	SWAP256x4(dv1, dv2, dv3, dv4);
	MINMAX32x4_alt2(0xFF00, 0x00ff, dv1, dv2, dv3, dv4);

	// distance 4 swaps
	SWAP128x4(dv1, dv2, dv3, dv4);
	MINMAX32x4_alt2(0xF0F0, 0x0f0f, dv1, dv2, dv3, dv4);

	// distance 2 swaps
	SWAP64x4(dv1, dv2, dv3, dv4);
	MINMAX32x4_alt2(0xCCCC, 0x3333, dv1, dv2, dv3, dv4);

	// adjacent swaps
	SWAP32x4(dv1, dv2, dv3, dv4);
	MINMAX32x4_alt2(0xAAAA, 0x5555, dv1, dv2, dv3, dv4);

	// phase 5: merge (all same compare 'dir')
	if (dir == 1)
	{	
		// distance 32 swaps - just compare the vecs
		t1 =  _mm512_max_epu32(dv1, dv3);
		t2 =  _mm512_max_epu32(dv2, dv4);
		dv3 = _mm512_min_epu32(dv1, dv3);
		dv4 = _mm512_min_epu32(dv2, dv4);
		dv1 = t1;
		dv2 = t2;

		// distance 16 swaps - just compare the vecs
		t1 =  _mm512_max_epu32(dv1, dv2);
		t2 =  _mm512_max_epu32(dv3, dv4);
		dv2 = _mm512_min_epu32(dv1, dv2);
		dv4 = _mm512_min_epu32(dv3, dv4);
		dv1 = t1;
		dv3 = t2;
		
		// distance 8 swaps - 256 bits
		SWAP256x4(dv1, dv2, dv3, dv4);
		MINMAX32x4(0x00ff, 0xFF00, dv1, dv2, dv3, dv4);
		
		// distance 4 swaps
		SWAP128x4(dv1, dv2, dv3, dv4);
		MINMAX32x4(0x0f0f, 0xF0F0, dv1, dv2, dv3, dv4);
		
		// distance 2 swaps
		SWAP64x4(dv1, dv2, dv3, dv4);
		MINMAX32x4(0x3333, 0xCCCC, dv1, dv2, dv3, dv4);
		
		// adjacent swaps
		SWAP32x4(dv1, dv2, dv3, dv4);
		MINMAX32x4(0x5555, 0xaaaa, dv1, dv2, dv3, dv4);

	}
	else
	{		
		// distance 32 swaps - just compare the vecs
		t1 =  _mm512_min_epu32(dv1, dv3);
		t2 =  _mm512_min_epu32(dv2, dv4);
		dv3 = _mm512_max_epu32(dv1, dv3);
		dv4 = _mm512_max_epu32(dv2, dv4);
		dv1 = t1;
		dv2 = t2;

		// distance 16 swaps - just compare the vecs
		t1 =  _mm512_min_epu32(dv1, dv2);
		t2 =  _mm512_min_epu32(dv3, dv4);
		dv2 = _mm512_max_epu32(dv1, dv2);
		dv4 = _mm512_max_epu32(dv3, dv4);
		dv1 = t1;
		dv3 = t2;

		// distance 8 swaps - 256 bits
		SWAP256x4(dv1, dv2, dv3, dv4);
		MINMAX32x4(0xFF00, 0x00ff, dv1, dv2, dv3, dv4);
		
		// distance 4 swaps
		SWAP128x4(dv1, dv2, dv3, dv4);
		MINMAX32x4(0xf0f0, 0x0f0f, dv1, dv2, dv3, dv4);

		// distance 2 swaps
		SWAP64x4(dv1, dv2, dv3, dv4);
		MINMAX32x4(0xCCCC, 0x3333, dv1, dv2, dv3, dv4);
		
		// adjacent swaps
		SWAP32x4(dv1, dv2, dv3, dv4);
		MINMAX32x4(0xaaaa, 0x5555, dv1, dv2, dv3, dv4);

	}

	_mm512_store_si512(data, dv1);
	_mm512_store_si512(data + 16, dv2);
	_mm512_store_si512(data + 32, dv3);
	_mm512_store_si512(data + 48, dv4);
	
	return;
}

void bitonic_merge32_dir_128(uint32_t* data, int dir) 
{
	// perform the final merge phase on 128 32-bit elements
	int i, j;

	__m512i t1;
	__m512i t2;
	__m512i t3;
	__m512i t4;
	__m512i t5;
	__m512i t6;
	__m512i t7;
	__m512i t8;
	__m512i dv1;
	__m512i dv2;
	__m512i dv3;
	__m512i dv4;
	__m512i dv5;
	__m512i dv6;
	__m512i dv7;
	__m512i dv8;
	__m512i dv1_swap;
	__m512i dv2_swap;
	__m512i dv3_swap;
	__m512i dv4_swap;
	__m512i dv5_swap;
	__m512i dv6_swap;
	__m512i dv7_swap;
	__m512i dv8_swap;
	__mmask16 m1;
	__mmask16 m2;
	__mmask16 m3;
	__mmask16 m4;
	__mmask16 m5;
	__mmask16 m6;
	__mmask16 m7;
	__mmask16 m8;

	dv1 = _mm512_load_si512(data);
	dv2 = _mm512_load_si512(data + 16);
	dv3 = _mm512_load_si512(data + 32);
	dv4 = _mm512_load_si512(data + 48);
	dv5 = _mm512_load_si512(data + 64);
	dv6 = _mm512_load_si512(data + 80);
	dv7 = _mm512_load_si512(data + 96);
	dv8 = _mm512_load_si512(data + 112);

	if (dir == 1)
	{
		// dist-64 swaps
		t1 =  _mm512_max_epu32(dv1, dv5);
		t2 =  _mm512_max_epu32(dv2, dv6);
		t3 =  _mm512_max_epu32(dv3, dv7);
		t4 =  _mm512_max_epu32(dv4, dv8);
		dv5 = _mm512_min_epu32(dv1, dv5);
		dv6 = _mm512_min_epu32(dv2, dv6);
		dv7 = _mm512_min_epu32(dv3, dv7);
		dv8 = _mm512_min_epu32(dv4, dv8);
		dv1 = t1;
		dv2 = t2;
		dv3 = t3;
		dv4 = t4;
		
		// distance 32 swaps - just compare the vecs
		t1 =  _mm512_max_epu32(dv1, dv3);
		t2 =  _mm512_max_epu32(dv2, dv4);
		t3 =  _mm512_max_epu32(dv5, dv7);
		t4 =  _mm512_max_epu32(dv6, dv8);
		dv3 = _mm512_min_epu32(dv1, dv3);
		dv4 = _mm512_min_epu32(dv2, dv4);
		dv7 = _mm512_min_epu32(dv5, dv7);
		dv8 = _mm512_min_epu32(dv6, dv8);
		dv1 = t1;
		dv2 = t2;
		dv5 = t3;
		dv6 = t4;

		// distance 16 swaps - just compare the vecs
		t1 =  _mm512_max_epu32(dv1, dv2);
		t2 =  _mm512_max_epu32(dv3, dv4);
		t3 =  _mm512_max_epu32(dv5, dv6);
		t4 =  _mm512_max_epu32(dv7, dv8);
		dv2 = _mm512_min_epu32(dv1, dv2);
		dv4 = _mm512_min_epu32(dv3, dv4);
		dv6 = _mm512_min_epu32(dv5, dv6);
		dv8 = _mm512_min_epu32(dv7, dv8);
		dv1 = t1;
		dv3 = t2;
		dv5 = t3;
		dv7 = t4;
		
		// distance 8 swaps - 256 bits
		SWAP256x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0x00ff, 0xFF00, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		// distance 4 swaps
		SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0x0f0f, 0xF0F0, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		// distance 2 swaps
		SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0x3333, 0xCCCC, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		// adjacent swaps
		SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0x5555, 0xaaaa, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	}
	else
	{	
		// dist-64 swaps
		t1 =  _mm512_min_epu32(dv1, dv5);
		t2 =  _mm512_min_epu32(dv2, dv6);
		t3 =  _mm512_min_epu32(dv3, dv7);
		t4 =  _mm512_min_epu32(dv4, dv8);
		dv5 = _mm512_max_epu32(dv1, dv5);
		dv6 = _mm512_max_epu32(dv2, dv6);
		dv7 = _mm512_max_epu32(dv3, dv7);
		dv8 = _mm512_max_epu32(dv4, dv8);
		dv1 = t1;
		dv2 = t2;
		dv3 = t3;
		dv4 = t4;
		
		// distance 32 swaps - just compare the vecs
		t1 =  _mm512_min_epu32(dv1, dv3);
		t2 =  _mm512_min_epu32(dv2, dv4);
		t3 =  _mm512_min_epu32(dv5, dv7);
		t4 =  _mm512_min_epu32(dv6, dv8);
		dv3 = _mm512_max_epu32(dv1, dv3);
		dv4 = _mm512_max_epu32(dv2, dv4);
		dv7 = _mm512_max_epu32(dv5, dv7);
		dv8 = _mm512_max_epu32(dv6, dv8);
		dv1 = t1;
		dv2 = t2;
		dv5 = t3;
		dv6 = t4;

		// distance 16 swaps - just compare the vecs
		t1 =  _mm512_min_epu32(dv1, dv2);
		t2 =  _mm512_min_epu32(dv3, dv4);
		t3 =  _mm512_min_epu32(dv5, dv6);
		t4 =  _mm512_min_epu32(dv7, dv8);
		dv2 = _mm512_max_epu32(dv1, dv2);
		dv4 = _mm512_max_epu32(dv3, dv4);
		dv6 = _mm512_max_epu32(dv5, dv6);
		dv8 = _mm512_max_epu32(dv7, dv8);
		dv1 = t1;
		dv3 = t2;
		dv5 = t3;
		dv7 = t4;

		// distance 8 swaps - 256 bits
		SWAP256x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0xFF00, 0x00ff, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		// distance 4 swaps
		SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0xf0f0, 0x0f0f, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

		// distance 2 swaps
		SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0xCCCC, 0x3333, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		// adjacent swaps
		SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0xaaaa, 0x5555, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	}

	_mm512_store_si512(data, dv1);
	_mm512_store_si512(data + 16, dv2);
	_mm512_store_si512(data + 32, dv3);
	_mm512_store_si512(data + 48, dv4);
	_mm512_store_si512(data + 64, dv5);
	_mm512_store_si512(data + 80, dv6);
	_mm512_store_si512(data + 96, dv7);
	_mm512_store_si512(data + 112, dv8);

	return;
}

void bitonic_sort32_dir_128(uint32_t* data, int dir) 
{
	// sort 128 32-bit elements
	int i, j;

	__m512i t1;
	__m512i t2;
	__m512i t3;
	__m512i t4;
	__m512i t5;
	__m512i t6;
	__m512i t7;
	__m512i t8;
	__m512i dv1;
	__m512i dv2;
	__m512i dv3;
	__m512i dv4;
	__m512i dv5;
	__m512i dv6;
	__m512i dv7;
	__m512i dv8;
	__m512i dv1_swap;
	__m512i dv2_swap;
	__m512i dv3_swap;
	__m512i dv4_swap;
	__m512i dv5_swap;
	__m512i dv6_swap;
	__m512i dv7_swap;
	__m512i dv8_swap;
	__mmask16 m1;
	__mmask16 m2;
	__mmask16 m3;
	__mmask16 m4;
	__mmask16 m5;
	__mmask16 m6;
	__mmask16 m7;
	__mmask16 m8;

	dv1 = _mm512_load_si512(data);
	dv2 = _mm512_load_si512(data + 16);
	dv3 = _mm512_load_si512(data + 32);
	dv4 = _mm512_load_si512(data + 48);
	dv5 = _mm512_load_si512(data + 64);
	dv6 = _mm512_load_si512(data + 80);
	dv7 = _mm512_load_si512(data + 96);
	dv8 = _mm512_load_si512(data + 112);

	// phase 0: dist-2 alternating compares ('CC')
	
	// adjacent swaps, alternating compares
	SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX(0x6666, 0x9999, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	// phase 1: dist-4 alternating compares ('F0')
	
	// distance 2 swaps
	SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX(0x3C3C, 0xC3C3, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8)

	// adjacent swaps
	SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX(0x5A5A, 0xA5A5, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// phase 2: dist-8 alternating compares ('FF00')

	// distance 4 swaps
	SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX(0x0FF0, 0xf00f, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// distance 2 swaps
	SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX(0x33CC, 0xcc33, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// adjacent swaps
	SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX(0x55AA, 0xaa55, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// phase 2: dist-16 alternating compares (alternating gt/lt)
	
	// distance 8 swaps (256 bits)
	SWAP256x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX_alt1(0xFF00, 0x00ff, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// distance 4 swaps
	SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX_alt1(0xF0F0, 0x0f0f, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// distance 2 swaps
	SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX_alt1(0xCCCC, 0x3333, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// adjacent swaps
	SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX_alt1(0xAAAA, 0x5555, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// phase 3: dist-32 alternating compares (alternating gtgt/ltlt)

	// distance 16 swaps - just compare the vecs
	t1 =  _mm512_min_epu32(dv1, dv2);
	t2 =  _mm512_max_epu32(dv3, dv4);
	t3 =  _mm512_min_epu32(dv5, dv6);
	t4 =  _mm512_max_epu32(dv7, dv8);
	dv2 = _mm512_max_epu32(dv1, dv2);
	dv4 = _mm512_min_epu32(dv3, dv4);
	dv6 = _mm512_max_epu32(dv5, dv6);
	dv8 = _mm512_min_epu32(dv7, dv8);
	dv1 = t1;
	dv3 = t2;
	dv5 = t3;
	dv7 = t4;

	// distance 8 swaps (256 bits)
	SWAP256x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX_alt2(0xFF00, 0x00ff, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// distance 4 swaps
	SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX_alt2(0xF0F0, 0x0f0f, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// distance 2 swaps
	SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX_alt2(0xCCCC, 0x3333, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// adjacent swaps
	SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX_alt2(0xAAAA, 0x5555, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	// phase 4: dist-64 alternating compares (alternating gtgtgtgt/ltltltlt)

	// distance 32 swaps - just compare the vecs
	t1 =  _mm512_min_epu32(dv1, dv3);
	t2 =  _mm512_min_epu32(dv2, dv4);
	t3 =  _mm512_max_epu32(dv5, dv7);
	t4 =  _mm512_max_epu32(dv6, dv8);
	dv3 = _mm512_max_epu32(dv1, dv3);
	dv4 = _mm512_max_epu32(dv2, dv4);
	dv7 = _mm512_min_epu32(dv5, dv7);
	dv8 = _mm512_min_epu32(dv6, dv8);
	dv1 = t1;
	dv2 = t2;
	dv5 = t3;
	dv6 = t4;
	
	// distance 16 swaps - just compare the vecs
	t1 =  _mm512_min_epu32(dv1, dv2);
	t2 =  _mm512_min_epu32(dv3, dv4);
	t3 =  _mm512_max_epu32(dv5, dv6);
	t4 =  _mm512_max_epu32(dv7, dv8);
	dv2 = _mm512_max_epu32(dv1, dv2);
	dv4 = _mm512_max_epu32(dv3, dv4);
	dv6 = _mm512_min_epu32(dv5, dv6);
	dv8 = _mm512_min_epu32(dv7, dv8);
	dv1 = t1;
	dv3 = t2;
	dv5 = t3;
	dv7 = t4;
	
	// distance 8 swaps (256 bits)
	SWAP256x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX_alt4(0xFF00, 0x00ff, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	// distance 4 swaps
	SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX_alt4(0xF0F0, 0x0f0f, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	// distance 2 swaps
	SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX_alt4(0xCCCC, 0x3333, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	// adjacent swaps
	SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX_alt4(0xAAAA, 0x5555, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);


	// phase 5: merge (all same compare 'dir')
	if (dir == 1)
	{
		// dist-64 swaps
		t1 =  _mm512_max_epu32(dv1, dv5);
		t2 =  _mm512_max_epu32(dv2, dv6);
		t3 =  _mm512_max_epu32(dv3, dv7);
		t4 =  _mm512_max_epu32(dv4, dv8);
		dv5 = _mm512_min_epu32(dv1, dv5);
		dv6 = _mm512_min_epu32(dv2, dv6);
		dv7 = _mm512_min_epu32(dv3, dv7);
		dv8 = _mm512_min_epu32(dv4, dv8);
		dv1 = t1;
		dv2 = t2;
		dv3 = t3;
		dv4 = t4;
		
		// distance 32 swaps - just compare the vecs
		t1 =  _mm512_max_epu32(dv1, dv3);
		t2 =  _mm512_max_epu32(dv2, dv4);
		t3 =  _mm512_max_epu32(dv5, dv7);
		t4 =  _mm512_max_epu32(dv6, dv8);
		dv3 = _mm512_min_epu32(dv1, dv3);
		dv4 = _mm512_min_epu32(dv2, dv4);
		dv7 = _mm512_min_epu32(dv5, dv7);
		dv8 = _mm512_min_epu32(dv6, dv8);
		dv1 = t1;
		dv2 = t2;
		dv5 = t3;
		dv6 = t4;

		// distance 16 swaps - just compare the vecs
		t1 =  _mm512_max_epu32(dv1, dv2);
		t2 =  _mm512_max_epu32(dv3, dv4);
		t3 =  _mm512_max_epu32(dv5, dv6);
		t4 =  _mm512_max_epu32(dv7, dv8);
		dv2 = _mm512_min_epu32(dv1, dv2);
		dv4 = _mm512_min_epu32(dv3, dv4);
		dv6 = _mm512_min_epu32(dv5, dv6);
		dv8 = _mm512_min_epu32(dv7, dv8);
		dv1 = t1;
		dv3 = t2;
		dv5 = t3;
		dv7 = t4;
		
		// distance 8 swaps - 256 bits
		SWAP256x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0x00ff, 0xFF00, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		// distance 4 swaps
		SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0x0f0f, 0xF0F0, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		// distance 2 swaps
		SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0x3333, 0xCCCC, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		// adjacent swaps
		SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0x5555, 0xaaaa, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	}
	else
	{	
		// dist-64 swaps
		t1 =  _mm512_min_epu32(dv1, dv5);
		t2 =  _mm512_min_epu32(dv2, dv6);
		t3 =  _mm512_min_epu32(dv3, dv7);
		t4 =  _mm512_min_epu32(dv4, dv8);
		dv5 = _mm512_max_epu32(dv1, dv5);
		dv6 = _mm512_max_epu32(dv2, dv6);
		dv7 = _mm512_max_epu32(dv3, dv7);
		dv8 = _mm512_max_epu32(dv4, dv8);
		dv1 = t1;
		dv2 = t2;
		dv3 = t3;
		dv4 = t4;
		
		// distance 32 swaps - just compare the vecs
		t1 =  _mm512_min_epu32(dv1, dv3);
		t2 =  _mm512_min_epu32(dv2, dv4);
		t3 =  _mm512_min_epu32(dv5, dv7);
		t4 =  _mm512_min_epu32(dv6, dv8);
		dv3 = _mm512_max_epu32(dv1, dv3);
		dv4 = _mm512_max_epu32(dv2, dv4);
		dv7 = _mm512_max_epu32(dv5, dv7);
		dv8 = _mm512_max_epu32(dv6, dv8);
		dv1 = t1;
		dv2 = t2;
		dv5 = t3;
		dv6 = t4;

		// distance 16 swaps - just compare the vecs
		t1 =  _mm512_min_epu32(dv1, dv2);
		t2 =  _mm512_min_epu32(dv3, dv4);
		t3 =  _mm512_min_epu32(dv5, dv6);
		t4 =  _mm512_min_epu32(dv7, dv8);
		dv2 = _mm512_max_epu32(dv1, dv2);
		dv4 = _mm512_max_epu32(dv3, dv4);
		dv6 = _mm512_max_epu32(dv5, dv6);
		dv8 = _mm512_max_epu32(dv7, dv8);
		dv1 = t1;
		dv3 = t2;
		dv5 = t3;
		dv7 = t4;

		// distance 8 swaps - 256 bits
		SWAP256x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0xFF00, 0x00ff, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		// distance 4 swaps
		SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0xf0f0, 0x0f0f, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

		// distance 2 swaps
		SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0xCCCC, 0x3333, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		// adjacent swaps
		SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0xaaaa, 0x5555, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	}

	_mm512_store_si512(data, dv1);
	_mm512_store_si512(data + 16, dv2);
	_mm512_store_si512(data + 32, dv3);
	_mm512_store_si512(data + 48, dv4);
	_mm512_store_si512(data + 64, dv5);
	_mm512_store_si512(data + 80, dv6);
	_mm512_store_si512(data + 96, dv7);
	_mm512_store_si512(data + 112, dv8);
	
	return;
}

void bitonic_merge32_dir_256(uint32_t* data, int dir) 
{
	// perform the final merge phase on 128 32-bit elements
	int i, j;

	__m512i dv1;
	__m512i dv2;
	__m512i dv3;
	__m512i dv4;
	__m512i dv5;
	__m512i dv6;
	__m512i dv7;
	__m512i dv8;
	__m512i dv9;
	__m512i dv10;
	__m512i dv11;
	__m512i dv12;
	__m512i dv13;
	__m512i dv14;
	__m512i dv15;
	__m512i dv16;
	__m512i dv1_swap;
	__m512i dv2_swap;
	__m512i dv3_swap;
	__m512i dv4_swap;
	__m512i dv5_swap;
	__m512i dv6_swap;
	__m512i dv7_swap;
	__m512i dv8_swap;
	__m512i t1;
	__m512i t2;
	__m512i t3;
	__m512i t4;
	__m512i t5;
	__m512i t6;
	__m512i t7;
	__m512i t8;
	__mmask16 m1;
	__mmask16 m2;
	__mmask16 m3;
	__mmask16 m4;
	__mmask16 m5;
	__mmask16 m6;
	__mmask16 m7;
	__mmask16 m8;

	dv1 = _mm512_load_si512(data);
	dv2 = _mm512_load_si512(data + 16);
	dv3 = _mm512_load_si512(data + 32);
	dv4 = _mm512_load_si512(data + 48);
	dv5 = _mm512_load_si512(data + 64);
	dv6 = _mm512_load_si512(data + 80);
	dv7 = _mm512_load_si512(data + 96);
	dv8 = _mm512_load_si512(data + 112);
	dv9 = _mm512_load_si512(data + 128);
	dv10 = _mm512_load_si512(data + 9 * 16);
	dv11 = _mm512_load_si512(data + 10 * 16);
	dv12 = _mm512_load_si512(data + 11 * 16);
	dv13 = _mm512_load_si512(data + 12 * 16);
	dv14 = _mm512_load_si512(data + 13 * 16);
	dv15 = _mm512_load_si512(data + 14 * 16);
	dv16 = _mm512_load_si512(data + 15 * 16);

	if (dir == 1)
	{
		// distance 128 swaps - just compare the vecs
		t1 =   _mm512_max_epu32(dv1, dv9 );
		t2 =   _mm512_max_epu32(dv2, dv10);
		t3 =   _mm512_max_epu32(dv3, dv11);
		t4 =   _mm512_max_epu32(dv4, dv12);
		dv9  = _mm512_min_epu32(dv1, dv9 );
		dv10 = _mm512_min_epu32(dv2, dv10);
		dv11 = _mm512_min_epu32(dv3, dv11);
		dv12 = _mm512_min_epu32(dv4, dv12);
		dv1 = t1;
		dv2 = t2;
		dv3 = t3;
		dv4 = t4;

		t1  =  _mm512_max_epu32(dv5, dv13);
		t2  =  _mm512_max_epu32(dv6, dv14);
		t3  =  _mm512_max_epu32(dv7, dv15);
		t4  =  _mm512_max_epu32(dv8, dv16);
		dv13 = _mm512_min_epu32(dv5, dv13);
		dv14 = _mm512_min_epu32(dv6, dv14);
		dv15 = _mm512_min_epu32(dv7, dv15);
		dv16 = _mm512_min_epu32(dv8, dv16);
		dv5 = t1;
		dv6 = t2;
		dv7 = t3;
		dv8 = t4;
		
		// dist-64 swaps
		t1 =  _mm512_max_epu32(dv1, dv5);
		t2 =  _mm512_max_epu32(dv2, dv6);
		t3 =  _mm512_max_epu32(dv3, dv7);
		t4 =  _mm512_max_epu32(dv4, dv8);
		dv5 = _mm512_min_epu32(dv1, dv5);
		dv6 = _mm512_min_epu32(dv2, dv6);
		dv7 = _mm512_min_epu32(dv3, dv7);
		dv8 = _mm512_min_epu32(dv4, dv8);
		dv1 = t1;
		dv2 = t2;
		dv3 = t3;
		dv4 = t4;

		t1  =  _mm512_max_epu32(dv9 , dv13);
		t2  =  _mm512_max_epu32(dv10, dv14);
		t3  =  _mm512_max_epu32(dv11, dv15);
		t4  =  _mm512_max_epu32(dv12, dv16);
		dv13 = _mm512_min_epu32(dv9 , dv13);
		dv14 = _mm512_min_epu32(dv10, dv14);
		dv15 = _mm512_min_epu32(dv11, dv15);
		dv16 = _mm512_min_epu32(dv12, dv16);
		dv9  = t1;
		dv10 = t2;
		dv11 = t3;
		dv12 = t4;
		
		
		// distance 32 swaps - just compare the vecs
		t1 =  _mm512_max_epu32(dv1, dv3);
		t2 =  _mm512_max_epu32(dv2, dv4);
		t3 =  _mm512_max_epu32(dv5, dv7);
		t4 =  _mm512_max_epu32(dv6, dv8);
		dv3 = _mm512_min_epu32(dv1, dv3);
		dv4 = _mm512_min_epu32(dv2, dv4);
		dv7 = _mm512_min_epu32(dv5, dv7);
		dv8 = _mm512_min_epu32(dv6, dv8);
		dv1 = t1;
		dv2 = t2;
		dv5 = t3;
		dv6 = t4;

		t1  =  _mm512_max_epu32(dv9 , dv11);
		t2  =  _mm512_max_epu32(dv10, dv12);
		t3  =  _mm512_max_epu32(dv13, dv15);
		t4  =  _mm512_max_epu32(dv14, dv16);
		dv11 = _mm512_min_epu32(dv9 , dv11);
		dv12 = _mm512_min_epu32(dv10, dv12);
		dv15 = _mm512_min_epu32(dv13, dv15);
		dv16 = _mm512_min_epu32(dv14, dv16);
		dv9  = t1;
		dv10 = t2;
		dv13 = t3;
		dv14 = t4;
		
		// distance 16 swaps - just compare the vecs
		t1 =  _mm512_max_epu32(dv1, dv2);
		t2 =  _mm512_max_epu32(dv3, dv4);
		t3 =  _mm512_max_epu32(dv5, dv6);
		t4 =  _mm512_max_epu32(dv7, dv8);
		dv2 = _mm512_min_epu32(dv1, dv2);
		dv4 = _mm512_min_epu32(dv3, dv4);
		dv6 = _mm512_min_epu32(dv5, dv6);
		dv8 = _mm512_min_epu32(dv7, dv8);
		dv1 = t1;
		dv3 = t2;
		dv5 = t3;
		dv7 = t4;

		t1  =  _mm512_max_epu32(dv9 , dv10);
		t2  =  _mm512_max_epu32(dv11, dv12);
		t3  =  _mm512_max_epu32(dv13, dv14);
		t4  =  _mm512_max_epu32(dv15, dv16);
		dv10 = _mm512_min_epu32(dv9 , dv10);
		dv12 = _mm512_min_epu32(dv11, dv12);
		dv14 = _mm512_min_epu32(dv13, dv14);
		dv16 = _mm512_min_epu32(dv15, dv16);
		dv9  = t1;
		dv11 = t2;
		dv13 = t3;
		dv15 = t4;
		
		// distance 8 swaps - 256 bits
		SWAP256x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0x00ff, 0xFF00, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		SWAP256x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		MINMAX(0x00ff, 0xFF00, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		
		// distance 4 swaps
		SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0x0f0f, 0xF0F0, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		SWAP128x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		MINMAX(0x0f0f, 0xF0F0, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		
		// distance 2 swaps
		SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0x3333, 0xCCCC, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		SWAP64x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		MINMAX(0x3333, 0xCCCC, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		
		// adjacent swaps
		SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0x5555, 0xaaaa, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		SWAP32x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		MINMAX(0x5555, 0xaaaa, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	}
	else
	{
		// distance 128 swaps - just compare the vecs
		t1 =  _mm512_min_epu32(dv1, dv9 );
		t2 =  _mm512_min_epu32(dv2, dv10);
		t3 =  _mm512_min_epu32(dv3, dv11);
		t4 =  _mm512_min_epu32(dv4, dv12);
		dv9  = _mm512_max_epu32(dv1, dv9 );
		dv10 = _mm512_max_epu32(dv2, dv10);
		dv11 = _mm512_max_epu32(dv3, dv11);
		dv12 = _mm512_max_epu32(dv4, dv12);
		dv1 = t1;
		dv2 = t2;
		dv3 = t3;
		dv4 = t4;

		t1  =  _mm512_min_epu32(dv5, dv13);
		t2  =  _mm512_min_epu32(dv6, dv14);
		t3  =  _mm512_min_epu32(dv7, dv15);
		t4  =  _mm512_min_epu32(dv8, dv16);
		dv13 = _mm512_max_epu32(dv5, dv13);
		dv14 = _mm512_max_epu32(dv6, dv14);
		dv15 = _mm512_max_epu32(dv7, dv15);
		dv16 = _mm512_max_epu32(dv8, dv16);
		dv5 = t1;
		dv6 = t2;
		dv7 = t3;
		dv8 = t4;
		
		// dist-64 swaps
		t1 =  _mm512_min_epu32(dv1, dv5);
		t2 =  _mm512_min_epu32(dv2, dv6);
		t3 =  _mm512_min_epu32(dv3, dv7);
		t4 =  _mm512_min_epu32(dv4, dv8);
		dv5 = _mm512_max_epu32(dv1, dv5);
		dv6 = _mm512_max_epu32(dv2, dv6);
		dv7 = _mm512_max_epu32(dv3, dv7);
		dv8 = _mm512_max_epu32(dv4, dv8);
		dv1 = t1;
		dv2 = t2;
		dv3 = t3;
		dv4 = t4;

		t1  =  _mm512_min_epu32(dv9 , dv13);
		t2  =  _mm512_min_epu32(dv10, dv14);
		t3  =  _mm512_min_epu32(dv11, dv15);
		t4  =  _mm512_min_epu32(dv12, dv16);
		dv13 = _mm512_max_epu32(dv9 , dv13);
		dv14 = _mm512_max_epu32(dv10, dv14);
		dv15 = _mm512_max_epu32(dv11, dv15);
		dv16 = _mm512_max_epu32(dv12, dv16);
		dv9  = t1;
		dv10 = t2;
		dv11 = t3;
		dv12 = t4;
		
		
		// distance 32 swaps - just compare the vecs
		t1 =  _mm512_min_epu32(dv1, dv3);
		t2 =  _mm512_min_epu32(dv2, dv4);
		t3 =  _mm512_min_epu32(dv5, dv7);
		t4 =  _mm512_min_epu32(dv6, dv8);
		dv3 = _mm512_max_epu32(dv1, dv3);
		dv4 = _mm512_max_epu32(dv2, dv4);
		dv7 = _mm512_max_epu32(dv5, dv7);
		dv8 = _mm512_max_epu32(dv6, dv8);
		dv1 = t1;
		dv2 = t2;
		dv5 = t3;
		dv6 = t4;

		t1  =  _mm512_min_epu32(dv9 , dv11);
		t2  =  _mm512_min_epu32(dv10, dv12);
		t3  =  _mm512_min_epu32(dv13, dv15);
		t4  =  _mm512_min_epu32(dv14, dv16);
		dv11 = _mm512_max_epu32(dv9 , dv11);
		dv12 = _mm512_max_epu32(dv10, dv12);
		dv15 = _mm512_max_epu32(dv13, dv15);
		dv16 = _mm512_max_epu32(dv14, dv16);
		dv9  = t1;
		dv10 = t2;
		dv13 = t3;
		dv14 = t4;
		
		// distance 16 swaps - just compare the vecs
		t1 =  _mm512_min_epu32(dv1, dv2);
		t2 =  _mm512_min_epu32(dv3, dv4);
		t3 =  _mm512_min_epu32(dv5, dv6);
		t4 =  _mm512_min_epu32(dv7, dv8);
		dv2 = _mm512_max_epu32(dv1, dv2);
		dv4 = _mm512_max_epu32(dv3, dv4);
		dv6 = _mm512_max_epu32(dv5, dv6);
		dv8 = _mm512_max_epu32(dv7, dv8);
		dv1 = t1;
		dv3 = t2;
		dv5 = t3;
		dv7 = t4;

		t1  =  _mm512_min_epu32(dv9 , dv10);
		t2  =  _mm512_min_epu32(dv11, dv12);
		t3  =  _mm512_min_epu32(dv13, dv14);
		t4  =  _mm512_min_epu32(dv15, dv16);
		dv10 = _mm512_max_epu32(dv9 , dv10);
		dv12 = _mm512_max_epu32(dv11, dv12);
		dv14 = _mm512_max_epu32(dv13, dv14);
		dv16 = _mm512_max_epu32(dv15, dv16);
		dv9  = t1;
		dv11 = t2;
		dv13 = t3;
		dv15 = t4;
		
		// distance 8 swaps - 256 bits
		SWAP256x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0xFF00, 0x00ff, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		SWAP256x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		MINMAX(0xFF00, 0x00ff, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		
		// distance 4 swaps
		SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0xf0f0, 0x0f0f, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		SWAP128x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		MINMAX(0xf0f0, 0x0f0f, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);

		// distance 2 swaps
		SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0xCCCC, 0x3333, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		SWAP64x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		MINMAX(0xCCCC, 0x3333, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		
		// adjacent swaps
		SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0xaaaa, 0x5555, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		SWAP32x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		MINMAX(0xaaaa, 0x5555, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);

	}

	_mm512_store_si512(data, dv1);
	_mm512_store_si512(data + 16, dv2);
	_mm512_store_si512(data + 32, dv3);
	_mm512_store_si512(data + 48, dv4);
	_mm512_store_si512(data + 64, dv5);
	_mm512_store_si512(data + 80, dv6);
	_mm512_store_si512(data + 96, dv7);
	_mm512_store_si512(data + 112, dv8);
	_mm512_store_si512(data + 128    , dv9);
	_mm512_store_si512(data + 9 * 16 , dv10);
	_mm512_store_si512(data + 10 * 16, dv11);
	_mm512_store_si512(data + 11 * 16, dv12);
	_mm512_store_si512(data + 12 * 16, dv13);
	_mm512_store_si512(data + 13 * 16, dv14);
	_mm512_store_si512(data + 14 * 16, dv15);
	_mm512_store_si512(data + 15 * 16, dv16);
	
	return;
}

void bitonic_sort32_dir_256(uint32_t* data, int dir) 
{
	// sort 256 32-bit elements
	int i, j;

	__m512i dv1;
	__m512i dv2;
	__m512i dv3;
	__m512i dv4;
	__m512i dv5;
	__m512i dv6;
	__m512i dv7;
	__m512i dv8;
	__m512i dv9;
	__m512i dv10;
	__m512i dv11;
	__m512i dv12;
	__m512i dv13;
	__m512i dv14;
	__m512i dv15;
	__m512i dv16;
	__m512i dv1_swap;
	__m512i dv2_swap;
	__m512i dv3_swap;
	__m512i dv4_swap;
	__m512i dv5_swap;
	__m512i dv6_swap;
	__m512i dv7_swap;
	__m512i dv8_swap;
	__m512i t1;
	__m512i t2;
	__m512i t3;
	__m512i t4;
	__m512i t5;
	__m512i t6;
	__m512i t7;
	__m512i t8;
	__mmask16 m1;
	__mmask16 m2;
	__mmask16 m3;
	__mmask16 m4;
	__mmask16 m5;
	__mmask16 m6;
	__mmask16 m7;
	__mmask16 m8;

	dv1 = _mm512_load_si512(data);
	dv2 = _mm512_load_si512(data + 16);
	dv3 = _mm512_load_si512(data + 32);
	dv4 = _mm512_load_si512(data + 48);
	dv5 = _mm512_load_si512(data + 64);
	dv6 = _mm512_load_si512(data + 80);
	dv7 = _mm512_load_si512(data + 96);
	dv8 = _mm512_load_si512(data + 112);
	dv9 = _mm512_load_si512(data + 128);
	dv10 = _mm512_load_si512(data + 9 * 16);
	dv11 = _mm512_load_si512(data + 10 * 16);
	dv12 = _mm512_load_si512(data + 11 * 16);
	dv13 = _mm512_load_si512(data + 12 * 16);
	dv14 = _mm512_load_si512(data + 13 * 16);
	dv15 = _mm512_load_si512(data + 14 * 16);
	dv16 = _mm512_load_si512(data + 15 * 16);

	// phase 0: dist-2 alternating compares ('CC')
	
	// adjacent swaps, alternating compares
	SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX(0x6666, 0x9999, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8)

	SWAP32x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	MINMAX(0x6666, 0x9999, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16)
	
	// phase 1: dist-4 alternating compares ('F0')
	
	// distance 2 swaps
	SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX(0x3C3C, 0xC3C3, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8)
	
	SWAP64x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	MINMAX(0x3C3C, 0xC3C3, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);

	// adjacent swaps
	SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX(0x5A5A, 0xA5A5, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	SWAP32x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	MINMAX(0x5A5A, 0xA5A5, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);

	// phase 2: dist-8 alternating compares ('FF00')

	// distance 4 swaps
	SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX(0x0FF0, 0xf00f, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	SWAP128x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	MINMAX(0x0FF0, 0xf00f, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);

	// distance 2 swaps
	SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX(0x33CC, 0xcc33, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	SWAP64x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	MINMAX(0x33CC, 0xcc33, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);

	// adjacent swaps
	SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX(0x55AA, 0xaa55, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	SWAP32x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	MINMAX(0x55AA, 0xaa55, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);

	// phase 2: dist-16 alternating compares (alternating gt/lt)
	
	// distance 8 swaps (256 bits)
	SWAP256x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX_alt1(0xFF00, 0x00ff, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	SWAP256x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	MINMAX_alt1(0xFF00, 0x00ff, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);

	// distance 4 swaps
	SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX_alt1(0xF0F0, 0x0f0f, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	SWAP128x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	MINMAX_alt1(0xF0F0, 0x0f0f, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);

	// distance 2 swaps
	SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX_alt1(0xCCCC, 0x3333, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	SWAP64x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	MINMAX_alt1(0xCCCC, 0x3333, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);

	// adjacent swaps
	SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX_alt1(0xAAAA, 0x5555, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	SWAP32x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	MINMAX_alt1(0xAAAA, 0x5555, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);

	// phase 3: dist-32 alternating compares (alternating gtgt/ltlt)

	// distance 16 swaps - just compare the vecs
	t1 =  _mm512_min_epu32(dv1, dv2);
	t2 =  _mm512_max_epu32(dv3, dv4);
	t3 =  _mm512_min_epu32(dv5, dv6);
	t4 =  _mm512_max_epu32(dv7, dv8);
	dv2 = _mm512_max_epu32(dv1, dv2);
	dv4 = _mm512_min_epu32(dv3, dv4);
	dv6 = _mm512_max_epu32(dv5, dv6);
	dv8 = _mm512_min_epu32(dv7, dv8);
	dv1 = t1;
	dv3 = t2;
	dv5 = t3;
	dv7 = t4;

	t1  =  _mm512_min_epu32(dv9 , dv10);
	t2  =  _mm512_max_epu32(dv11, dv12);
	t3  =  _mm512_min_epu32(dv13, dv14);
	t4  =  _mm512_max_epu32(dv15, dv16);
	dv10 = _mm512_max_epu32(dv9 , dv10);
	dv12 = _mm512_min_epu32(dv11, dv12);
	dv14 = _mm512_max_epu32(dv13, dv14);
	dv16 = _mm512_min_epu32(dv15, dv16);
	dv9  = t1;
	dv11 = t2;
	dv13 = t3;
	dv15 = t4;

	// distance 8 swaps (256 bits)
	SWAP256x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX_alt2(0xFF00, 0x00ff, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	SWAP256x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	MINMAX_alt2(0xFF00, 0x00ff, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);

	// distance 4 swaps
	SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX_alt2(0xF0F0, 0x0f0f, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	SWAP128x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	MINMAX_alt2(0xF0F0, 0x0f0f, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);

	// distance 2 swaps
	SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX_alt2(0xCCCC, 0x3333, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	SWAP64x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	MINMAX_alt2(0xCCCC, 0x3333, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);

	// adjacent swaps
	SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX_alt2(0xAAAA, 0x5555, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	SWAP32x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	MINMAX_alt2(0xAAAA, 0x5555, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	
	// phase 4: dist-64 alternating compares (alternating gtgtgtgt/ltltltlt)

	// distance 32 swaps - just compare the vecs
	t1 =  _mm512_min_epu32(dv1, dv3);
	t2 =  _mm512_min_epu32(dv2, dv4);
	t3 =  _mm512_max_epu32(dv5, dv7);
	t4 =  _mm512_max_epu32(dv6, dv8);
	dv3 = _mm512_max_epu32(dv1, dv3);
	dv4 = _mm512_max_epu32(dv2, dv4);
	dv7 = _mm512_min_epu32(dv5, dv7);
	dv8 = _mm512_min_epu32(dv6, dv8);
	dv1 = t1;
	dv2 = t2;
	dv5 = t3;
	dv6 = t4;

	t1  =  _mm512_min_epu32(dv9 , dv11);
	t2  =  _mm512_min_epu32(dv10, dv12);
	t3  =  _mm512_max_epu32(dv13, dv15);
	t4  =  _mm512_max_epu32(dv14, dv16);
	dv11 = _mm512_max_epu32(dv9 , dv11);
	dv12 = _mm512_max_epu32(dv10, dv12);
	dv15 = _mm512_min_epu32(dv13, dv15);
	dv16 = _mm512_min_epu32(dv14, dv16);
	dv9  = t1;
	dv10 = t2;
	dv13 = t3;
	dv14 = t4;
	
	// distance 16 swaps - just compare the vecs
	t1 =  _mm512_min_epu32(dv1, dv2);
	t2 =  _mm512_min_epu32(dv3, dv4);
	t3 =  _mm512_max_epu32(dv5, dv6);
	t4 =  _mm512_max_epu32(dv7, dv8);
	dv2 = _mm512_max_epu32(dv1, dv2);
	dv4 = _mm512_max_epu32(dv3, dv4);
	dv6 = _mm512_min_epu32(dv5, dv6);
	dv8 = _mm512_min_epu32(dv7, dv8);
	dv1 = t1;
	dv3 = t2;
	dv5 = t3;
	dv7 = t4;

	t1  =  _mm512_min_epu32(dv9 , dv10);
	t2  =  _mm512_min_epu32(dv11, dv12);
	t3  =  _mm512_max_epu32(dv13, dv14);
	t4  =  _mm512_max_epu32(dv15, dv16);
	dv10 = _mm512_max_epu32(dv9 , dv10);
	dv12 = _mm512_max_epu32(dv11, dv12);
	dv14 = _mm512_min_epu32(dv13, dv14);
	dv16 = _mm512_min_epu32(dv15, dv16);
	dv9  = t1;
	dv11 = t2;
	dv13 = t3;
	dv15 = t4;
	
	// distance 8 swaps (256 bits)
	SWAP256x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX_alt4(0xFF00, 0x00ff, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	SWAP256x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	MINMAX_alt4(0xFF00, 0x00ff, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	
	// distance 4 swaps
	SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX_alt4(0xF0F0, 0x0f0f, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	SWAP128x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	MINMAX_alt4(0xF0F0, 0x0f0f, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	
	// distance 2 swaps
	SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX_alt4(0xCCCC, 0x3333, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	SWAP64x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	MINMAX_alt4(0xCCCC, 0x3333, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	
	// adjacent swaps
	SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX_alt4(0xAAAA, 0x5555, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	SWAP32x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	MINMAX_alt4(0xAAAA, 0x5555, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	
	// phase 5: dist-128 alternating compares (alternating gtgtgtgtgtgtgtgt/ltltltltltltltlt)
	
	// dist-64 swaps
	t1 =  _mm512_min_epu32(dv1, dv5);
	t2 =  _mm512_min_epu32(dv2, dv6);
	t3 =  _mm512_min_epu32(dv3, dv7);
	t4 =  _mm512_min_epu32(dv4, dv8);
	dv5 = _mm512_max_epu32(dv1, dv5);
	dv6 = _mm512_max_epu32(dv2, dv6);
	dv7 = _mm512_max_epu32(dv3, dv7);
	dv8 = _mm512_max_epu32(dv4, dv8);
	dv1 = t1;
	dv2 = t2;
	dv3 = t3;
	dv4 = t4;

	t1  =  _mm512_max_epu32(dv9 , dv13);
	t2  =  _mm512_max_epu32(dv10, dv14);
	t3  =  _mm512_max_epu32(dv11, dv15);
	t4  =  _mm512_max_epu32(dv12, dv16);
	dv13 = _mm512_min_epu32(dv9 , dv13);
	dv14 = _mm512_min_epu32(dv10, dv14);
	dv15 = _mm512_min_epu32(dv11, dv15);
	dv16 = _mm512_min_epu32(dv12, dv16);
	dv9  = t1;
	dv10 = t2;
	dv11 = t3;
	dv12 = t4;
	
	
	// distance 32 swaps - just compare the vecs
	t1 =  _mm512_min_epu32(dv1, dv3);
	t2 =  _mm512_min_epu32(dv2, dv4);
	t3 =  _mm512_min_epu32(dv5, dv7);
	t4 =  _mm512_min_epu32(dv6, dv8);
	dv3 = _mm512_max_epu32(dv1, dv3);
	dv4 = _mm512_max_epu32(dv2, dv4);
	dv7 = _mm512_max_epu32(dv5, dv7);
	dv8 = _mm512_max_epu32(dv6, dv8);
	dv1 = t1;
	dv2 = t2;
	dv5 = t3;
	dv6 = t4;

	t1  =  _mm512_max_epu32(dv9 , dv11);
	t2  =  _mm512_max_epu32(dv10, dv12);
	t3  =  _mm512_max_epu32(dv13, dv15);
	t4  =  _mm512_max_epu32(dv14, dv16);
	dv11 = _mm512_min_epu32(dv9 , dv11);
	dv12 = _mm512_min_epu32(dv10, dv12);
	dv15 = _mm512_min_epu32(dv13, dv15);
	dv16 = _mm512_min_epu32(dv14, dv16);
	dv9  = t1;
	dv10 = t2;
	dv13 = t3;
	dv14 = t4;
	
	// distance 16 swaps - just compare the vecs
	t1 =  _mm512_min_epu32(dv1, dv2);
	t2 =  _mm512_min_epu32(dv3, dv4);
	t3 =  _mm512_min_epu32(dv5, dv6);
	t4 =  _mm512_min_epu32(dv7, dv8);
	dv2 = _mm512_max_epu32(dv1, dv2);
	dv4 = _mm512_max_epu32(dv3, dv4);
	dv6 = _mm512_max_epu32(dv5, dv6);
	dv8 = _mm512_max_epu32(dv7, dv8);
	dv1 = t1;
	dv3 = t2;
	dv5 = t3;
	dv7 = t4;

	t1  =  _mm512_max_epu32(dv9 , dv10);
	t2  =  _mm512_max_epu32(dv11, dv12);
	t3  =  _mm512_max_epu32(dv13, dv14);
	t4  =  _mm512_max_epu32(dv15, dv16);
	dv10 = _mm512_min_epu32(dv9 , dv10);
	dv12 = _mm512_min_epu32(dv11, dv12);
	dv14 = _mm512_min_epu32(dv13, dv14);
	dv16 = _mm512_min_epu32(dv15, dv16);
	dv9  = t1;
	dv11 = t2;
	dv13 = t3;
	dv15 = t4;
	
	// distance 8 swaps (256 bits)
	SWAP256x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX(0xff00, 0x00ff, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	SWAP256x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	MINMAX(0x00ff, 0xff00, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	
	// distance 4 swaps
	SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX(0xF0F0, 0x0f0f, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	SWAP128x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	MINMAX(0x0f0f, 0xF0F0, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	
	// distance 2 swaps
	SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX(0xcccc, 0x3333, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	SWAP64x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	MINMAX(0x3333, 0xcccc, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	
	// adjacent swaps
	SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	MINMAX(0xAAAA, 0x5555, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	
	SWAP32x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	MINMAX(0x5555, 0xAAAA, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	

	// phase 6: merge (all same compare 'dir')
	if (dir == 1)
	{
		// distance 128 swaps - just compare the vecs
		t1 =   _mm512_max_epu32(dv1, dv9 );
		t2 =   _mm512_max_epu32(dv2, dv10);
		t3 =   _mm512_max_epu32(dv3, dv11);
		t4 =   _mm512_max_epu32(dv4, dv12);
		dv9  = _mm512_min_epu32(dv1, dv9 );
		dv10 = _mm512_min_epu32(dv2, dv10);
		dv11 = _mm512_min_epu32(dv3, dv11);
		dv12 = _mm512_min_epu32(dv4, dv12);
		dv1 = t1;
		dv2 = t2;
		dv3 = t3;
		dv4 = t4;

		t1  =  _mm512_max_epu32(dv5, dv13);
		t2  =  _mm512_max_epu32(dv6, dv14);
		t3  =  _mm512_max_epu32(dv7, dv15);
		t4  =  _mm512_max_epu32(dv8, dv16);
		dv13 = _mm512_min_epu32(dv5, dv13);
		dv14 = _mm512_min_epu32(dv6, dv14);
		dv15 = _mm512_min_epu32(dv7, dv15);
		dv16 = _mm512_min_epu32(dv8, dv16);
		dv5 = t1;
		dv6 = t2;
		dv7 = t3;
		dv8 = t4;
		
		// dist-64 swaps
		t1 =  _mm512_max_epu32(dv1, dv5);
		t2 =  _mm512_max_epu32(dv2, dv6);
		t3 =  _mm512_max_epu32(dv3, dv7);
		t4 =  _mm512_max_epu32(dv4, dv8);
		dv5 = _mm512_min_epu32(dv1, dv5);
		dv6 = _mm512_min_epu32(dv2, dv6);
		dv7 = _mm512_min_epu32(dv3, dv7);
		dv8 = _mm512_min_epu32(dv4, dv8);
		dv1 = t1;
		dv2 = t2;
		dv3 = t3;
		dv4 = t4;

		t1  =  _mm512_max_epu32(dv9 , dv13);
		t2  =  _mm512_max_epu32(dv10, dv14);
		t3  =  _mm512_max_epu32(dv11, dv15);
		t4  =  _mm512_max_epu32(dv12, dv16);
		dv13 = _mm512_min_epu32(dv9 , dv13);
		dv14 = _mm512_min_epu32(dv10, dv14);
		dv15 = _mm512_min_epu32(dv11, dv15);
		dv16 = _mm512_min_epu32(dv12, dv16);
		dv9  = t1;
		dv10 = t2;
		dv11 = t3;
		dv12 = t4;
		
		
		// distance 32 swaps - just compare the vecs
		t1 =  _mm512_max_epu32(dv1, dv3);
		t2 =  _mm512_max_epu32(dv2, dv4);
		t3 =  _mm512_max_epu32(dv5, dv7);
		t4 =  _mm512_max_epu32(dv6, dv8);
		dv3 = _mm512_min_epu32(dv1, dv3);
		dv4 = _mm512_min_epu32(dv2, dv4);
		dv7 = _mm512_min_epu32(dv5, dv7);
		dv8 = _mm512_min_epu32(dv6, dv8);
		dv1 = t1;
		dv2 = t2;
		dv5 = t3;
		dv6 = t4;

		t1  =  _mm512_max_epu32(dv9 , dv11);
		t2  =  _mm512_max_epu32(dv10, dv12);
		t3  =  _mm512_max_epu32(dv13, dv15);
		t4  =  _mm512_max_epu32(dv14, dv16);
		dv11 = _mm512_min_epu32(dv9 , dv11);
		dv12 = _mm512_min_epu32(dv10, dv12);
		dv15 = _mm512_min_epu32(dv13, dv15);
		dv16 = _mm512_min_epu32(dv14, dv16);
		dv9  = t1;
		dv10 = t2;
		dv13 = t3;
		dv14 = t4;
		
		// distance 16 swaps - just compare the vecs
		t1 =  _mm512_max_epu32(dv1, dv2);
		t2 =  _mm512_max_epu32(dv3, dv4);
		t3 =  _mm512_max_epu32(dv5, dv6);
		t4 =  _mm512_max_epu32(dv7, dv8);
		dv2 = _mm512_min_epu32(dv1, dv2);
		dv4 = _mm512_min_epu32(dv3, dv4);
		dv6 = _mm512_min_epu32(dv5, dv6);
		dv8 = _mm512_min_epu32(dv7, dv8);
		dv1 = t1;
		dv3 = t2;
		dv5 = t3;
		dv7 = t4;

		t1  =  _mm512_max_epu32(dv9 , dv10);
		t2  =  _mm512_max_epu32(dv11, dv12);
		t3  =  _mm512_max_epu32(dv13, dv14);
		t4  =  _mm512_max_epu32(dv15, dv16);
		dv10 = _mm512_min_epu32(dv9 , dv10);
		dv12 = _mm512_min_epu32(dv11, dv12);
		dv14 = _mm512_min_epu32(dv13, dv14);
		dv16 = _mm512_min_epu32(dv15, dv16);
		dv9  = t1;
		dv11 = t2;
		dv13 = t3;
		dv15 = t4;
		
		// distance 8 swaps - 256 bits
		SWAP256x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0x00ff, 0xFF00, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		SWAP256x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		MINMAX(0x00ff, 0xFF00, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		
		// distance 4 swaps
		SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0x0f0f, 0xF0F0, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		SWAP128x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		MINMAX(0x0f0f, 0xF0F0, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		
		// distance 2 swaps
		SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0x3333, 0xCCCC, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		SWAP64x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		MINMAX(0x3333, 0xCCCC, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		
		// adjacent swaps
		SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0x5555, 0xaaaa, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		SWAP32x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		MINMAX(0x5555, 0xaaaa, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
	}
	else
	{
		// distance 128 swaps - just compare the vecs
		t1 =  _mm512_min_epu32(dv1, dv9 );
		t2 =  _mm512_min_epu32(dv2, dv10);
		t3 =  _mm512_min_epu32(dv3, dv11);
		t4 =  _mm512_min_epu32(dv4, dv12);
		dv9  = _mm512_max_epu32(dv1, dv9 );
		dv10 = _mm512_max_epu32(dv2, dv10);
		dv11 = _mm512_max_epu32(dv3, dv11);
		dv12 = _mm512_max_epu32(dv4, dv12);
		dv1 = t1;
		dv2 = t2;
		dv3 = t3;
		dv4 = t4;

		t1  =  _mm512_min_epu32(dv5, dv13);
		t2  =  _mm512_min_epu32(dv6, dv14);
		t3  =  _mm512_min_epu32(dv7, dv15);
		t4  =  _mm512_min_epu32(dv8, dv16);
		dv13 = _mm512_max_epu32(dv5, dv13);
		dv14 = _mm512_max_epu32(dv6, dv14);
		dv15 = _mm512_max_epu32(dv7, dv15);
		dv16 = _mm512_max_epu32(dv8, dv16);
		dv5 = t1;
		dv6 = t2;
		dv7 = t3;
		dv8 = t4;
		
		// dist-64 swaps
		t1 =  _mm512_min_epu32(dv1, dv5);
		t2 =  _mm512_min_epu32(dv2, dv6);
		t3 =  _mm512_min_epu32(dv3, dv7);
		t4 =  _mm512_min_epu32(dv4, dv8);
		dv5 = _mm512_max_epu32(dv1, dv5);
		dv6 = _mm512_max_epu32(dv2, dv6);
		dv7 = _mm512_max_epu32(dv3, dv7);
		dv8 = _mm512_max_epu32(dv4, dv8);
		dv1 = t1;
		dv2 = t2;
		dv3 = t3;
		dv4 = t4;

		t1  =  _mm512_min_epu32(dv9 , dv13);
		t2  =  _mm512_min_epu32(dv10, dv14);
		t3  =  _mm512_min_epu32(dv11, dv15);
		t4  =  _mm512_min_epu32(dv12, dv16);
		dv13 = _mm512_max_epu32(dv9 , dv13);
		dv14 = _mm512_max_epu32(dv10, dv14);
		dv15 = _mm512_max_epu32(dv11, dv15);
		dv16 = _mm512_max_epu32(dv12, dv16);
		dv9  = t1;
		dv10 = t2;
		dv11 = t3;
		dv12 = t4;
		
		
		// distance 32 swaps - just compare the vecs
		t1 =  _mm512_min_epu32(dv1, dv3);
		t2 =  _mm512_min_epu32(dv2, dv4);
		t3 =  _mm512_min_epu32(dv5, dv7);
		t4 =  _mm512_min_epu32(dv6, dv8);
		dv3 = _mm512_max_epu32(dv1, dv3);
		dv4 = _mm512_max_epu32(dv2, dv4);
		dv7 = _mm512_max_epu32(dv5, dv7);
		dv8 = _mm512_max_epu32(dv6, dv8);
		dv1 = t1;
		dv2 = t2;
		dv5 = t3;
		dv6 = t4;

		t1  =  _mm512_min_epu32(dv9 , dv11);
		t2  =  _mm512_min_epu32(dv10, dv12);
		t3  =  _mm512_min_epu32(dv13, dv15);
		t4  =  _mm512_min_epu32(dv14, dv16);
		dv11 = _mm512_max_epu32(dv9 , dv11);
		dv12 = _mm512_max_epu32(dv10, dv12);
		dv15 = _mm512_max_epu32(dv13, dv15);
		dv16 = _mm512_max_epu32(dv14, dv16);
		dv9  = t1;
		dv10 = t2;
		dv13 = t3;
		dv14 = t4;
		
		// distance 16 swaps - just compare the vecs
		t1 =  _mm512_min_epu32(dv1, dv2);
		t2 =  _mm512_min_epu32(dv3, dv4);
		t3 =  _mm512_min_epu32(dv5, dv6);
		t4 =  _mm512_min_epu32(dv7, dv8);
		dv2 = _mm512_max_epu32(dv1, dv2);
		dv4 = _mm512_max_epu32(dv3, dv4);
		dv6 = _mm512_max_epu32(dv5, dv6);
		dv8 = _mm512_max_epu32(dv7, dv8);
		dv1 = t1;
		dv3 = t2;
		dv5 = t3;
		dv7 = t4;

		t1  =  _mm512_min_epu32(dv9 , dv10);
		t2  =  _mm512_min_epu32(dv11, dv12);
		t3  =  _mm512_min_epu32(dv13, dv14);
		t4  =  _mm512_min_epu32(dv15, dv16);
		dv10 = _mm512_max_epu32(dv9 , dv10);
		dv12 = _mm512_max_epu32(dv11, dv12);
		dv14 = _mm512_max_epu32(dv13, dv14);
		dv16 = _mm512_max_epu32(dv15, dv16);
		dv9  = t1;
		dv11 = t2;
		dv13 = t3;
		dv15 = t4;
		
		// distance 8 swaps - 256 bits
		SWAP256x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0xFF00, 0x00ff, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		SWAP256x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		MINMAX(0xFF00, 0x00ff, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		
		// distance 4 swaps
		SWAP128x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0xf0f0, 0x0f0f, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		SWAP128x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		MINMAX(0xf0f0, 0x0f0f, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);

		// distance 2 swaps
		SWAP64x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0xCCCC, 0x3333, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		SWAP64x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		MINMAX(0xCCCC, 0x3333, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		
		// adjacent swaps
		SWAP32x8(dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		MINMAX(0xaaaa, 0x5555, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		SWAP32x8(dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);
		MINMAX(0xaaaa, 0x5555, dv9, dv10, dv11, dv12, dv13, dv14, dv15, dv16);

	}

	_mm512_store_si512(data, dv1);
	_mm512_store_si512(data + 16, dv2);
	_mm512_store_si512(data + 32, dv3);
	_mm512_store_si512(data + 48, dv4);
	_mm512_store_si512(data + 64, dv5);
	_mm512_store_si512(data + 80, dv6);
	_mm512_store_si512(data + 96, dv7);
	_mm512_store_si512(data + 112, dv8);
	_mm512_store_si512(data + 128    , dv9);
	_mm512_store_si512(data + 9 * 16 , dv10);
	_mm512_store_si512(data + 10 * 16, dv11);
	_mm512_store_si512(data + 11 * 16, dv12);
	_mm512_store_si512(data + 12 * 16, dv13);
	_mm512_store_si512(data + 13 * 16, dv14);
	_mm512_store_si512(data + 14 * 16, dv15);
	_mm512_store_si512(data + 15 * 16, dv16);
	
	return;
}

void bitonic_merge_dir_64(uint64_t* data, int dir)
{
	// perform the final merge pass on 64 elements
	__m512i t1;
	__m512i t2;
	__m512i t3;
	__m512i t4;
	__m512i t5;
	__m512i t6;
	__m512i t7;
	__m512i t8;
	__m512i dv1;
	__m512i dv2;
	__m512i dv3;
	__m512i dv4;
	__m512i dv5;
	__m512i dv6;
	__m512i dv7;
	__m512i dv8;
	__m512i dv1_swap;
	__m512i dv2_swap;
	__m512i dv3_swap;
	__m512i dv4_swap;
	__m512i dv5_swap;
	__m512i dv6_swap;
	__m512i dv7_swap;
	__m512i dv8_swap;
	__mmask8 m1;
	__mmask8 m2;
	__mmask8 m3;
	__mmask8 m4;
	__mmask8 m5;
	__mmask8 m6;
	__mmask8 m7;
	__mmask8 m8;

	dv1 = _mm512_load_si512(data);
	dv2 = _mm512_load_si512(data + 8);
	dv3 = _mm512_load_si512(data + 16);
	dv4 = _mm512_load_si512(data + 24);
	dv5 = _mm512_load_si512(data + 32);
	dv6 = _mm512_load_si512(data + 40);
	dv7 = _mm512_load_si512(data + 48);
	dv8 = _mm512_load_si512(data + 56);
	
	
	if (dir == 1)
	{
		// distance 32 swaps - just compare the vecs
		m1 = _mm512_cmplt_epu64_mask(dv1, dv5);
		m2 = _mm512_cmplt_epu64_mask(dv2, dv6);
		m3 = _mm512_cmplt_epu64_mask(dv3, dv7);
		m4 = _mm512_cmplt_epu64_mask(dv4, dv8);

		t1 = _mm512_mask_blend_epi64(m1, dv1, dv5);
		t2 = _mm512_mask_blend_epi64(m2, dv2, dv6);
		t3 = _mm512_mask_blend_epi64(m3, dv3, dv7);
		t4 = _mm512_mask_blend_epi64(m4, dv4, dv8);
		dv5 = _mm512_mask_blend_epi64(m1, dv5, dv1);
		dv6 = _mm512_mask_blend_epi64(m2, dv6, dv2);
		dv7 = _mm512_mask_blend_epi64(m3, dv7, dv3);
		dv8 = _mm512_mask_blend_epi64(m4, dv8, dv4);
		dv1 = t1;
		dv2 = t2;
		dv3 = t3;
		dv4 = t4;
		
		// distance 16 swaps - just compare the vecs
		m1 = _mm512_cmplt_epu64_mask(dv1, dv3);
		m2 = _mm512_cmplt_epu64_mask(dv2, dv4);
		m3 = _mm512_cmplt_epu64_mask(dv5, dv7);
		m4 = _mm512_cmplt_epu64_mask(dv6, dv8);

		t1 = _mm512_mask_blend_epi64(m1, dv1, dv3);
		t2 = _mm512_mask_blend_epi64(m2, dv2, dv4);
		t3 = _mm512_mask_blend_epi64(m3, dv5, dv7);
		t4 = _mm512_mask_blend_epi64(m4, dv6, dv8);

		dv3 = _mm512_mask_blend_epi64(m1, dv3, dv1);
		dv4 = _mm512_mask_blend_epi64(m2, dv4, dv2);
		dv7 = _mm512_mask_blend_epi64(m3, dv7, dv5);
		dv8 = _mm512_mask_blend_epi64(m4, dv8, dv6);

		dv1 = t1;
		dv2 = t2;
		dv5 = t3;
		dv6 = t4;

		// distance 8 swaps - just compare the vecs
		m1 = _mm512_cmplt_epu64_mask(dv1, dv2);
		m2 = _mm512_cmplt_epu64_mask(dv3, dv4);
		m3 = _mm512_cmplt_epu64_mask(dv5, dv6);
		m4 = _mm512_cmplt_epu64_mask(dv7, dv8);

		t1 = _mm512_mask_blend_epi64(m1, dv1, dv2);
		t2 = _mm512_mask_blend_epi64(m2, dv3, dv4);
		t3 = _mm512_mask_blend_epi64(m3, dv5, dv6);
		t4 = _mm512_mask_blend_epi64(m4, dv7, dv8);
		
		dv2 = _mm512_mask_blend_epi64(m1, dv2, dv1);
		dv4 = _mm512_mask_blend_epi64(m2, dv4, dv3);
		dv6 = _mm512_mask_blend_epi64(m3, dv6, dv5);
		dv8 = _mm512_mask_blend_epi64(m4, dv8, dv7);
		
		dv1 = t1;
		dv3 = t2;
		dv5 = t3;
		dv7 = t4;

		// distance 4 swaps
		dv1_swap = SWAP256(dv1);
		dv2_swap = SWAP256(dv2);
		dv3_swap = SWAP256(dv3);
		dv4_swap = SWAP256(dv4);
		dv5_swap = SWAP256(dv5);
		dv6_swap = SWAP256(dv6);
		dv7_swap = SWAP256(dv7);
		dv8_swap = SWAP256(dv8);
		
		m1 = _mm512_cmplt_epu64_mask(dv1, dv1_swap);
		m2 = _mm512_cmplt_epu64_mask(dv2, dv2_swap);
		m3 = _mm512_cmplt_epu64_mask(dv3, dv3_swap);
		m4 = _mm512_cmplt_epu64_mask(dv4, dv4_swap);
		m5 = _mm512_cmplt_epu64_mask(dv5, dv5_swap);
		m6 = _mm512_cmplt_epu64_mask(dv6, dv6_swap);
		m7 = _mm512_cmplt_epu64_mask(dv7, dv7_swap);
		m8 = _mm512_cmplt_epu64_mask(dv8, dv8_swap);

		// 'F0' for the swap between dist-4 lanes, non alternating
		dv1 = _mm512_mask_blend_epi64(m1 ^ 0xF0, dv1, dv1_swap);
		dv2 = _mm512_mask_blend_epi64(m2 ^ 0xF0, dv2, dv2_swap);
		dv3 = _mm512_mask_blend_epi64(m3 ^ 0xF0, dv3, dv3_swap);
		dv4 = _mm512_mask_blend_epi64(m4 ^ 0xF0, dv4, dv4_swap);
		dv5 = _mm512_mask_blend_epi64(m5 ^ 0xF0, dv5, dv5_swap);
		dv6 = _mm512_mask_blend_epi64(m6 ^ 0xF0, dv6, dv6_swap);
		dv7 = _mm512_mask_blend_epi64(m7 ^ 0xF0, dv7, dv7_swap);
		dv8 = _mm512_mask_blend_epi64(m8 ^ 0xF0, dv8, dv8_swap);

		// distance 2 swaps
		dv1_swap = SWAP128(dv1);
		dv2_swap = SWAP128(dv2);
		dv3_swap = SWAP128(dv3);
		dv4_swap = SWAP128(dv4);
		dv5_swap = SWAP128(dv5);
		dv6_swap = SWAP128(dv6);
		dv7_swap = SWAP128(dv7);
		dv8_swap = SWAP128(dv8);
		
		m1 = _mm512_cmplt_epu64_mask(dv1, dv1_swap);
		m2 = _mm512_cmplt_epu64_mask(dv2, dv2_swap);
		m3 = _mm512_cmplt_epu64_mask(dv3, dv3_swap);
		m4 = _mm512_cmplt_epu64_mask(dv4, dv4_swap);
		m5 = _mm512_cmplt_epu64_mask(dv5, dv5_swap);
		m6 = _mm512_cmplt_epu64_mask(dv6, dv6_swap);
		m7 = _mm512_cmplt_epu64_mask(dv7, dv7_swap);
		m8 = _mm512_cmplt_epu64_mask(dv8, dv8_swap);

		// 'CC' for the swap between dist-2 lanes, non alternating
		dv1 = _mm512_mask_blend_epi64(m1 ^ 0xCC, dv1, dv1_swap);
		dv2 = _mm512_mask_blend_epi64(m2 ^ 0xCC, dv2, dv2_swap);
		dv3 = _mm512_mask_blend_epi64(m3 ^ 0xCC, dv3, dv3_swap);
		dv4 = _mm512_mask_blend_epi64(m4 ^ 0xCC, dv4, dv4_swap);
		dv5 = _mm512_mask_blend_epi64(m5 ^ 0xCC, dv5, dv5_swap);
		dv6 = _mm512_mask_blend_epi64(m6 ^ 0xCC, dv6, dv6_swap);
		dv7 = _mm512_mask_blend_epi64(m7 ^ 0xCC, dv7, dv7_swap);
		dv8 = _mm512_mask_blend_epi64(m8 ^ 0xCC, dv8, dv8_swap);

		// adjacent swaps
		dv1_swap = SWAP64(dv1);
		dv2_swap = SWAP64(dv2);
		dv3_swap = SWAP64(dv3);
		dv4_swap = SWAP64(dv4);
		dv5_swap = SWAP64(dv5);
		dv6_swap = SWAP64(dv6);
		dv7_swap = SWAP64(dv7);
		dv8_swap = SWAP64(dv8);
		
		m1 = _mm512_cmplt_epu64_mask(dv1, dv1_swap);
		m2 = _mm512_cmplt_epu64_mask(dv2, dv2_swap);
		m3 = _mm512_cmplt_epu64_mask(dv3, dv3_swap);
		m4 = _mm512_cmplt_epu64_mask(dv4, dv4_swap);
		m5 = _mm512_cmplt_epu64_mask(dv5, dv5_swap);
		m6 = _mm512_cmplt_epu64_mask(dv6, dv6_swap);
		m7 = _mm512_cmplt_epu64_mask(dv7, dv7_swap);
		m8 = _mm512_cmplt_epu64_mask(dv8, dv8_swap);

		// 'AA' for the swap between adjacent lanes
		dv1 = _mm512_mask_blend_epi64(m1 ^ 0xAA, dv1, dv1_swap);
		dv2 = _mm512_mask_blend_epi64(m2 ^ 0xAA, dv2, dv2_swap);
		dv3 = _mm512_mask_blend_epi64(m3 ^ 0xAA, dv3, dv3_swap);
		dv4 = _mm512_mask_blend_epi64(m4 ^ 0xAA, dv4, dv4_swap);
		dv5 = _mm512_mask_blend_epi64(m5 ^ 0xAA, dv5, dv5_swap);
		dv6 = _mm512_mask_blend_epi64(m6 ^ 0xAA, dv6, dv6_swap);
		dv7 = _mm512_mask_blend_epi64(m7 ^ 0xAA, dv7, dv7_swap);
		dv8 = _mm512_mask_blend_epi64(m8 ^ 0xAA, dv8, dv8_swap);

	}
	else
	{

		// distance 32 swaps - just compare the vecs
		m1 = _mm512_cmpgt_epu64_mask(dv1, dv5);
		m2 = _mm512_cmpgt_epu64_mask(dv2, dv6);
		m3 = _mm512_cmpgt_epu64_mask(dv3, dv7);
		m4 = _mm512_cmpgt_epu64_mask(dv4, dv8);

		t1 = _mm512_mask_blend_epi64(m1, dv1, dv5);
		t2 = _mm512_mask_blend_epi64(m2, dv2, dv6);
		t3 = _mm512_mask_blend_epi64(m3, dv3, dv7);
		t4 = _mm512_mask_blend_epi64(m4, dv4, dv8);
		dv5 = _mm512_mask_blend_epi64(m1, dv5, dv1);
		dv6 = _mm512_mask_blend_epi64(m2, dv6, dv2);
		dv7 = _mm512_mask_blend_epi64(m3, dv7, dv3);
		dv8 = _mm512_mask_blend_epi64(m4, dv8, dv4);
		dv1 = t1;
		dv2 = t2;
		dv3 = t3;
		dv4 = t4;
		
		// distance 16 swaps - just compare the vecs
		m1 = _mm512_cmpgt_epu64_mask(dv1, dv3);
		m2 = _mm512_cmpgt_epu64_mask(dv2, dv4);
		m3 = _mm512_cmpgt_epu64_mask(dv5, dv7);
		m4 = _mm512_cmpgt_epu64_mask(dv6, dv8);

		t1 = _mm512_mask_blend_epi64(m1, dv1, dv3);
		t2 = _mm512_mask_blend_epi64(m2, dv2, dv4);
		t3 = _mm512_mask_blend_epi64(m3, dv5, dv7);
		t4 = _mm512_mask_blend_epi64(m4, dv6, dv8);

		dv3 = _mm512_mask_blend_epi64(m1, dv3, dv1);
		dv4 = _mm512_mask_blend_epi64(m2, dv4, dv2);
		dv7 = _mm512_mask_blend_epi64(m3, dv7, dv5);
		dv8 = _mm512_mask_blend_epi64(m4, dv8, dv6);

		dv1 = t1;
		dv2 = t2;
		dv5 = t3;
		dv6 = t4;

		// distance 8 swaps - just compare the vecs
		m1 = _mm512_cmpgt_epu64_mask(dv1, dv2);
		m2 = _mm512_cmpgt_epu64_mask(dv3, dv4);
		m3 = _mm512_cmpgt_epu64_mask(dv5, dv6);
		m4 = _mm512_cmpgt_epu64_mask(dv7, dv8);

		t1 = _mm512_mask_blend_epi64(m1, dv1, dv2);
		t2 = _mm512_mask_blend_epi64(m2, dv3, dv4);
		t3 = _mm512_mask_blend_epi64(m3, dv5, dv6);
		t4 = _mm512_mask_blend_epi64(m4, dv7, dv8);
		
		dv2 = _mm512_mask_blend_epi64(m1, dv2, dv1);
		dv4 = _mm512_mask_blend_epi64(m2, dv4, dv3);
		dv6 = _mm512_mask_blend_epi64(m3, dv6, dv5);
		dv8 = _mm512_mask_blend_epi64(m4, dv8, dv7);
		
		dv1 = t1;
		dv3 = t2;
		dv5 = t3;
		dv7 = t4;

		// distance 4 swaps
		dv1_swap = SWAP256(dv1);
		dv2_swap = SWAP256(dv2);
		dv3_swap = SWAP256(dv3);
		dv4_swap = SWAP256(dv4);
		dv5_swap = SWAP256(dv5);
		dv6_swap = SWAP256(dv6);
		dv7_swap = SWAP256(dv7);
		dv8_swap = SWAP256(dv8);
		
		m1 = _mm512_cmpgt_epu64_mask(dv1, dv1_swap);
		m2 = _mm512_cmpgt_epu64_mask(dv2, dv2_swap);
		m3 = _mm512_cmpgt_epu64_mask(dv3, dv3_swap);
		m4 = _mm512_cmpgt_epu64_mask(dv4, dv4_swap);
		m5 = _mm512_cmpgt_epu64_mask(dv5, dv5_swap);
		m6 = _mm512_cmpgt_epu64_mask(dv6, dv6_swap);
		m7 = _mm512_cmpgt_epu64_mask(dv7, dv7_swap);
		m8 = _mm512_cmpgt_epu64_mask(dv8, dv8_swap);

		// 'F0' for the swap between dist-4 lanes, non alternating
		dv1 = _mm512_mask_blend_epi64(m1 ^ 0xF0, dv1, dv1_swap);
		dv2 = _mm512_mask_blend_epi64(m2 ^ 0xF0, dv2, dv2_swap);
		dv3 = _mm512_mask_blend_epi64(m3 ^ 0xF0, dv3, dv3_swap);
		dv4 = _mm512_mask_blend_epi64(m4 ^ 0xF0, dv4, dv4_swap);
		dv5 = _mm512_mask_blend_epi64(m5 ^ 0xF0, dv5, dv5_swap);
		dv6 = _mm512_mask_blend_epi64(m6 ^ 0xF0, dv6, dv6_swap);
		dv7 = _mm512_mask_blend_epi64(m7 ^ 0xF0, dv7, dv7_swap);
		dv8 = _mm512_mask_blend_epi64(m8 ^ 0xF0, dv8, dv8_swap);

		// distance 2 swaps
		dv1_swap = SWAP128(dv1);
		dv2_swap = SWAP128(dv2);
		dv3_swap = SWAP128(dv3);
		dv4_swap = SWAP128(dv4);
		dv5_swap = SWAP128(dv5);
		dv6_swap = SWAP128(dv6);
		dv7_swap = SWAP128(dv7);
		dv8_swap = SWAP128(dv8);
		
		m1 = _mm512_cmpgt_epu64_mask(dv1, dv1_swap);
		m2 = _mm512_cmpgt_epu64_mask(dv2, dv2_swap);
		m3 = _mm512_cmpgt_epu64_mask(dv3, dv3_swap);
		m4 = _mm512_cmpgt_epu64_mask(dv4, dv4_swap);
		m5 = _mm512_cmpgt_epu64_mask(dv5, dv5_swap);
		m6 = _mm512_cmpgt_epu64_mask(dv6, dv6_swap);
		m7 = _mm512_cmpgt_epu64_mask(dv7, dv7_swap);
		m8 = _mm512_cmpgt_epu64_mask(dv8, dv8_swap);

		// 'CC' for the swap between dist-2 lanes, non alternating
		dv1 = _mm512_mask_blend_epi64(m1 ^ 0xCC, dv1, dv1_swap);
		dv2 = _mm512_mask_blend_epi64(m2 ^ 0xCC, dv2, dv2_swap);
		dv3 = _mm512_mask_blend_epi64(m3 ^ 0xCC, dv3, dv3_swap);
		dv4 = _mm512_mask_blend_epi64(m4 ^ 0xCC, dv4, dv4_swap);
		dv5 = _mm512_mask_blend_epi64(m5 ^ 0xCC, dv5, dv5_swap);
		dv6 = _mm512_mask_blend_epi64(m6 ^ 0xCC, dv6, dv6_swap);
		dv7 = _mm512_mask_blend_epi64(m7 ^ 0xCC, dv7, dv7_swap);
		dv8 = _mm512_mask_blend_epi64(m8 ^ 0xCC, dv8, dv8_swap);

		// adjacent swaps
		dv1_swap = SWAP64(dv1);
		dv2_swap = SWAP64(dv2);
		dv3_swap = SWAP64(dv3);
		dv4_swap = SWAP64(dv4);
		dv5_swap = SWAP64(dv5);
		dv6_swap = SWAP64(dv6);
		dv7_swap = SWAP64(dv7);
		dv8_swap = SWAP64(dv8);
		
		m1 = _mm512_cmpgt_epu64_mask(dv1, dv1_swap);
		m2 = _mm512_cmpgt_epu64_mask(dv2, dv2_swap);
		m3 = _mm512_cmpgt_epu64_mask(dv3, dv3_swap);
		m4 = _mm512_cmpgt_epu64_mask(dv4, dv4_swap);
		m5 = _mm512_cmpgt_epu64_mask(dv5, dv5_swap);
		m6 = _mm512_cmpgt_epu64_mask(dv6, dv6_swap);
		m7 = _mm512_cmpgt_epu64_mask(dv7, dv7_swap);
		m8 = _mm512_cmpgt_epu64_mask(dv8, dv8_swap);

		// 'AA' for the swap between adjacent lanes
		dv1 = _mm512_mask_blend_epi64(m1 ^ 0xAA, dv1, dv1_swap);
		dv2 = _mm512_mask_blend_epi64(m2 ^ 0xAA, dv2, dv2_swap);
		dv3 = _mm512_mask_blend_epi64(m3 ^ 0xAA, dv3, dv3_swap);
		dv4 = _mm512_mask_blend_epi64(m4 ^ 0xAA, dv4, dv4_swap);
		dv5 = _mm512_mask_blend_epi64(m5 ^ 0xAA, dv5, dv5_swap);
		dv6 = _mm512_mask_blend_epi64(m6 ^ 0xAA, dv6, dv6_swap);
		dv7 = _mm512_mask_blend_epi64(m7 ^ 0xAA, dv7, dv7_swap);
		dv8 = _mm512_mask_blend_epi64(m8 ^ 0xAA, dv8, dv8_swap);
	}

	_mm512_store_si512(data, dv1);
	_mm512_store_si512(data + 8, dv2);
	_mm512_store_si512(data + 16, dv3);
	_mm512_store_si512(data + 24, dv4);
	_mm512_store_si512(data + 32, dv5);
	_mm512_store_si512(data + 40, dv6);
	_mm512_store_si512(data + 48, dv7);
	_mm512_store_si512(data + 56, dv8);
	
	return;
}

void bitonic_sort_dir_64(uint64_t* data, int dir) 
{
	// sort 64 elements of 64-bit data
	int i, j;

	__m512i t1;
	__m512i t2;
	__m512i t3;
	__m512i t4;
	__m512i t5;
	__m512i t6;
	__m512i t7;
	__m512i t8;
	__m512i dv1;
	__m512i dv2;
	__m512i dv3;
	__m512i dv4;
	__m512i dv5;
	__m512i dv6;
	__m512i dv7;
	__m512i dv8;
	__m512i dv1_swap;
	__m512i dv2_swap;
	__m512i dv3_swap;
	__m512i dv4_swap;
	__m512i dv5_swap;
	__m512i dv6_swap;
	__m512i dv7_swap;
	__m512i dv8_swap;
	__mmask8 m1;
	__mmask8 m2;
	__mmask8 m3;
	__mmask8 m4;
	__mmask8 m5;
	__mmask8 m6;
	__mmask8 m7;
	__mmask8 m8;

	dv1 = _mm512_load_si512(data);
	dv2 = _mm512_load_si512(data + 8);
	dv3 = _mm512_load_si512(data + 16);
	dv4 = _mm512_load_si512(data + 24);
	dv5 = _mm512_load_si512(data + 32);
	dv6 = _mm512_load_si512(data + 40);
	dv7 = _mm512_load_si512(data + 48);
	dv8 = _mm512_load_si512(data + 56);

	// phase 1 : dist-2 alternating compares
	
	// adjacent swaps, alternating compares
	dv1_swap = SWAP64(dv1);
	dv2_swap = SWAP64(dv2);
	dv3_swap = SWAP64(dv3);
	dv4_swap = SWAP64(dv4);
	dv5_swap = SWAP64(dv5);
	dv6_swap = SWAP64(dv6);
	dv7_swap = SWAP64(dv7);
	dv8_swap = SWAP64(dv8);
	
	m1 = _mm512_cmpgt_epu64_mask(dv1, dv1_swap);
	m2 = _mm512_cmpgt_epu64_mask(dv2, dv2_swap);
	m3 = _mm512_cmpgt_epu64_mask(dv3, dv3_swap);
	m4 = _mm512_cmpgt_epu64_mask(dv4, dv4_swap);
	m5 = _mm512_cmpgt_epu64_mask(dv5, dv5_swap);
	m6 = _mm512_cmpgt_epu64_mask(dv6, dv6_swap);
	m7 = _mm512_cmpgt_epu64_mask(dv7, dv7_swap);
	m8 = _mm512_cmpgt_epu64_mask(dv8, dv8_swap);
	
	// 'AA' for the swap between adjacent lanes ^ 'CC' to alternate gt vs. lt --> 0x66
	dv1 = _mm512_mask_blend_epi64(m1 ^ 0x66, dv1, dv1_swap);
	dv2 = _mm512_mask_blend_epi64(m2 ^ 0x66, dv2, dv2_swap);
	dv3 = _mm512_mask_blend_epi64(m3 ^ 0x66, dv3, dv3_swap);
	dv4 = _mm512_mask_blend_epi64(m4 ^ 0x66, dv4, dv4_swap);
	dv5 = _mm512_mask_blend_epi64(m5 ^ 0x66, dv5, dv5_swap);
	dv6 = _mm512_mask_blend_epi64(m6 ^ 0x66, dv6, dv6_swap);
	dv7 = _mm512_mask_blend_epi64(m7 ^ 0x66, dv7, dv7_swap);
	dv8 = _mm512_mask_blend_epi64(m8 ^ 0x66, dv8, dv8_swap);
	
	// phase 2 : dist-4 alternating compares

	// distance 2 swaps
	dv1_swap = SWAP128(dv1);
	dv2_swap = SWAP128(dv2);
	dv3_swap = SWAP128(dv3);
	dv4_swap = SWAP128(dv4);
	dv5_swap = SWAP128(dv5);
	dv6_swap = SWAP128(dv6);
	dv7_swap = SWAP128(dv7);
	dv8_swap = SWAP128(dv8);
	
	m1 = _mm512_cmpgt_epu64_mask(dv1, dv1_swap);
	m2 = _mm512_cmpgt_epu64_mask(dv2, dv2_swap);
	m3 = _mm512_cmpgt_epu64_mask(dv3, dv3_swap);
	m4 = _mm512_cmpgt_epu64_mask(dv4, dv4_swap);
	m5 = _mm512_cmpgt_epu64_mask(dv5, dv5_swap);
	m6 = _mm512_cmpgt_epu64_mask(dv6, dv6_swap);
	m7 = _mm512_cmpgt_epu64_mask(dv7, dv7_swap);
	m8 = _mm512_cmpgt_epu64_mask(dv8, dv8_swap);

	// 'CC' for the swap between dist-2 lanes ^ 'F0' to alternate gt vs. lt --> 0xC3
	dv1 = _mm512_mask_blend_epi64(m1 ^ 0x3C, dv1, dv1_swap);
	dv2 = _mm512_mask_blend_epi64(m2 ^ 0x3C, dv2, dv2_swap);
	dv3 = _mm512_mask_blend_epi64(m3 ^ 0x3C, dv3, dv3_swap);
	dv4 = _mm512_mask_blend_epi64(m4 ^ 0x3C, dv4, dv4_swap);
	dv5 = _mm512_mask_blend_epi64(m5 ^ 0x3C, dv5, dv5_swap);
	dv6 = _mm512_mask_blend_epi64(m6 ^ 0x3C, dv6, dv6_swap);
	dv7 = _mm512_mask_blend_epi64(m7 ^ 0x3C, dv7, dv7_swap);
	dv8 = _mm512_mask_blend_epi64(m8 ^ 0x3C, dv8, dv8_swap);

	// adjacent swaps
	dv1_swap = SWAP64(dv1);
	dv2_swap = SWAP64(dv2);
	dv3_swap = SWAP64(dv3);
	dv4_swap = SWAP64(dv4);
	dv5_swap = SWAP64(dv5);
	dv6_swap = SWAP64(dv6);
	dv7_swap = SWAP64(dv7);
	dv8_swap = SWAP64(dv8);
	
	m1 = _mm512_cmpgt_epu64_mask(dv1, dv1_swap);
	m2 = _mm512_cmpgt_epu64_mask(dv2, dv2_swap);
	m3 = _mm512_cmpgt_epu64_mask(dv3, dv3_swap);
	m4 = _mm512_cmpgt_epu64_mask(dv4, dv4_swap);
	m5 = _mm512_cmpgt_epu64_mask(dv5, dv5_swap);
	m6 = _mm512_cmpgt_epu64_mask(dv6, dv6_swap);
	m7 = _mm512_cmpgt_epu64_mask(dv7, dv7_swap);
	m8 = _mm512_cmpgt_epu64_mask(dv8, dv8_swap);

	// 'AA' for the swap between adjacent lanes ^ 'F0' to alternate gt vs. lt --> 0xA5
	dv1 = _mm512_mask_blend_epi64(m1 ^ 0x5A, dv1, dv1_swap);
	dv2 = _mm512_mask_blend_epi64(m2 ^ 0x5A, dv2, dv2_swap);
	dv3 = _mm512_mask_blend_epi64(m3 ^ 0x5A, dv3, dv3_swap);
	dv4 = _mm512_mask_blend_epi64(m4 ^ 0x5A, dv4, dv4_swap);
	dv5 = _mm512_mask_blend_epi64(m5 ^ 0x5A, dv5, dv5_swap);
	dv6 = _mm512_mask_blend_epi64(m6 ^ 0x5A, dv6, dv6_swap);
	dv7 = _mm512_mask_blend_epi64(m7 ^ 0x5A, dv7, dv7_swap);
	dv8 = _mm512_mask_blend_epi64(m8 ^ 0x5A, dv8, dv8_swap);


	// phase 3 : dist-8 alternating compares

	// distance 4 swaps
	dv1_swap = SWAP256(dv1);
	dv2_swap = SWAP256(dv2);
	dv3_swap = SWAP256(dv3);
	dv4_swap = SWAP256(dv4);
	dv5_swap = SWAP256(dv5);
	dv6_swap = SWAP256(dv6);
	dv7_swap = SWAP256(dv7);
	dv8_swap = SWAP256(dv8);
	
	m1 = _mm512_cmpgt_epu64_mask(dv1, dv1_swap);
	m2 = _mm512_cmplt_epu64_mask(dv2, dv2_swap);
	m3 = _mm512_cmpgt_epu64_mask(dv3, dv3_swap);
	m4 = _mm512_cmplt_epu64_mask(dv4, dv4_swap);
	m5 = _mm512_cmpgt_epu64_mask(dv5, dv5_swap);
	m6 = _mm512_cmplt_epu64_mask(dv6, dv6_swap);
	m7 = _mm512_cmpgt_epu64_mask(dv7, dv7_swap);
	m8 = _mm512_cmplt_epu64_mask(dv8, dv8_swap);

	// '0F' for the swap between dist-4 lanes, non alternating
	dv1 = _mm512_mask_blend_epi64(m1 ^ 0xF0, dv1, dv1_swap);
	dv2 = _mm512_mask_blend_epi64(m2 ^ 0xF0, dv2, dv2_swap);
	dv3 = _mm512_mask_blend_epi64(m3 ^ 0xF0, dv3, dv3_swap);
	dv4 = _mm512_mask_blend_epi64(m4 ^ 0xF0, dv4, dv4_swap);
	dv5 = _mm512_mask_blend_epi64(m5 ^ 0xF0, dv5, dv5_swap);
	dv6 = _mm512_mask_blend_epi64(m6 ^ 0xF0, dv6, dv6_swap);
	dv7 = _mm512_mask_blend_epi64(m7 ^ 0xF0, dv7, dv7_swap);
	dv8 = _mm512_mask_blend_epi64(m8 ^ 0xF0, dv8, dv8_swap);

	// distance 2 swaps
	dv1_swap = SWAP128(dv1);
	dv2_swap = SWAP128(dv2);
	dv3_swap = SWAP128(dv3);
	dv4_swap = SWAP128(dv4);
	dv5_swap = SWAP128(dv5);
	dv6_swap = SWAP128(dv6);
	dv7_swap = SWAP128(dv7);
	dv8_swap = SWAP128(dv8);
	
	m1 = _mm512_cmpgt_epu64_mask(dv1, dv1_swap);
	m2 = _mm512_cmplt_epu64_mask(dv2, dv2_swap);
	m3 = _mm512_cmpgt_epu64_mask(dv3, dv3_swap);
	m4 = _mm512_cmplt_epu64_mask(dv4, dv4_swap);
	m5 = _mm512_cmpgt_epu64_mask(dv5, dv5_swap);
	m6 = _mm512_cmplt_epu64_mask(dv6, dv6_swap);
	m7 = _mm512_cmpgt_epu64_mask(dv7, dv7_swap);
	m8 = _mm512_cmplt_epu64_mask(dv8, dv8_swap);

	// 'CC' for the swap between dist-2 lanes, non alternating
	dv1 = _mm512_mask_blend_epi64(m1 ^ 0xCC, dv1, dv1_swap);
	dv2 = _mm512_mask_blend_epi64(m2 ^ 0xCC, dv2, dv2_swap);
	dv3 = _mm512_mask_blend_epi64(m3 ^ 0xCC, dv3, dv3_swap);
	dv4 = _mm512_mask_blend_epi64(m4 ^ 0xCC, dv4, dv4_swap);
	dv5 = _mm512_mask_blend_epi64(m5 ^ 0xCC, dv5, dv5_swap);
	dv6 = _mm512_mask_blend_epi64(m6 ^ 0xCC, dv6, dv6_swap);
	dv7 = _mm512_mask_blend_epi64(m7 ^ 0xCC, dv7, dv7_swap);
	dv8 = _mm512_mask_blend_epi64(m8 ^ 0xCC, dv8, dv8_swap);

	// adjacent swaps
	dv1_swap = SWAP64(dv1);
	dv2_swap = SWAP64(dv2);
	dv3_swap = SWAP64(dv3);
	dv4_swap = SWAP64(dv4);
	dv5_swap = SWAP64(dv5);
	dv6_swap = SWAP64(dv6);
	dv7_swap = SWAP64(dv7);
	dv8_swap = SWAP64(dv8);
	
	m1 = _mm512_cmpgt_epu64_mask(dv1, dv1_swap);
	m2 = _mm512_cmplt_epu64_mask(dv2, dv2_swap);
	m3 = _mm512_cmpgt_epu64_mask(dv3, dv3_swap);
	m4 = _mm512_cmplt_epu64_mask(dv4, dv4_swap);
	m5 = _mm512_cmpgt_epu64_mask(dv5, dv5_swap);
	m6 = _mm512_cmplt_epu64_mask(dv6, dv6_swap);
	m7 = _mm512_cmpgt_epu64_mask(dv7, dv7_swap);
	m8 = _mm512_cmplt_epu64_mask(dv8, dv8_swap);

	// 'AA' for the swap between adjacent lanes, non alternating
	dv1 = _mm512_mask_blend_epi64(m1 ^ 0xAA, dv1, dv1_swap);
	dv2 = _mm512_mask_blend_epi64(m2 ^ 0xAA, dv2, dv2_swap);
	dv3 = _mm512_mask_blend_epi64(m3 ^ 0xAA, dv3, dv3_swap);
	dv4 = _mm512_mask_blend_epi64(m4 ^ 0xAA, dv4, dv4_swap);
	dv5 = _mm512_mask_blend_epi64(m5 ^ 0xAA, dv5, dv5_swap);
	dv6 = _mm512_mask_blend_epi64(m6 ^ 0xAA, dv6, dv6_swap);
	dv7 = _mm512_mask_blend_epi64(m7 ^ 0xAA, dv7, dv7_swap);
	dv8 = _mm512_mask_blend_epi64(m8 ^ 0xAA, dv8, dv8_swap);


	// phase 4 : dist-16 alternating compares

	// distance 8 swaps - just compare the vecs
	m1 = _mm512_cmpgt_epu64_mask(dv1, dv2);
	m3 = _mm512_cmplt_epu64_mask(dv3, dv4);
	m5 = _mm512_cmpgt_epu64_mask(dv5, dv6);
	m7 = _mm512_cmplt_epu64_mask(dv7, dv8);
	t1 = _mm512_mask_blend_epi64(m1, dv1, dv2);
	t2 = _mm512_mask_blend_epi64(m3, dv3, dv4);
	t3 = _mm512_mask_blend_epi64(m5, dv5, dv6);
	t4 = _mm512_mask_blend_epi64(m7, dv7, dv8);
	dv2 = _mm512_mask_blend_epi64(m1, dv2, dv1);
	dv4 = _mm512_mask_blend_epi64(m3, dv4, dv3);
	dv6 = _mm512_mask_blend_epi64(m5, dv6, dv5);
	dv8 = _mm512_mask_blend_epi64(m7, dv8, dv7);
	dv1 = t1;
	dv3 = t2;
	dv5 = t3;
	dv7 = t4;

	// distance 4 swaps
	dv1_swap = SWAP256(dv1);
	dv2_swap = SWAP256(dv2);
	dv3_swap = SWAP256(dv3);
	dv4_swap = SWAP256(dv4);
	dv5_swap = SWAP256(dv5);
	dv6_swap = SWAP256(dv6);
	dv7_swap = SWAP256(dv7);
	dv8_swap = SWAP256(dv8);
	
	m1 = _mm512_cmpgt_epu64_mask(dv1, dv1_swap);
	m2 = _mm512_cmpgt_epu64_mask(dv2, dv2_swap);
	m3 = _mm512_cmplt_epu64_mask(dv3, dv3_swap);
	m4 = _mm512_cmplt_epu64_mask(dv4, dv4_swap);
	m5 = _mm512_cmpgt_epu64_mask(dv5, dv5_swap);
	m6 = _mm512_cmpgt_epu64_mask(dv6, dv6_swap);
	m7 = _mm512_cmplt_epu64_mask(dv7, dv7_swap);
	m8 = _mm512_cmplt_epu64_mask(dv8, dv8_swap);

	// '0F' for the swap between dist-4 lanes, non alternating
	dv1 = _mm512_mask_blend_epi64(m1 ^ 0xF0, dv1, dv1_swap);
	dv2 = _mm512_mask_blend_epi64(m2 ^ 0xF0, dv2, dv2_swap);
	dv3 = _mm512_mask_blend_epi64(m3 ^ 0xF0, dv3, dv3_swap);
	dv4 = _mm512_mask_blend_epi64(m4 ^ 0xF0, dv4, dv4_swap);
	dv5 = _mm512_mask_blend_epi64(m5 ^ 0xF0, dv5, dv5_swap);
	dv6 = _mm512_mask_blend_epi64(m6 ^ 0xF0, dv6, dv6_swap);
	dv7 = _mm512_mask_blend_epi64(m7 ^ 0xF0, dv7, dv7_swap);
	dv8 = _mm512_mask_blend_epi64(m8 ^ 0xF0, dv8, dv8_swap);

	// distance 2 swaps
	dv1_swap = SWAP128(dv1);
	dv2_swap = SWAP128(dv2);
	dv3_swap = SWAP128(dv3);
	dv4_swap = SWAP128(dv4);
	dv5_swap = SWAP128(dv5);
	dv6_swap = SWAP128(dv6);
	dv7_swap = SWAP128(dv7);
	dv8_swap = SWAP128(dv8);
	
	m1 = _mm512_cmpgt_epu64_mask(dv1, dv1_swap);
	m2 = _mm512_cmpgt_epu64_mask(dv2, dv2_swap);
	m3 = _mm512_cmplt_epu64_mask(dv3, dv3_swap);
	m4 = _mm512_cmplt_epu64_mask(dv4, dv4_swap);
	m5 = _mm512_cmpgt_epu64_mask(dv5, dv5_swap);
	m6 = _mm512_cmpgt_epu64_mask(dv6, dv6_swap);
	m7 = _mm512_cmplt_epu64_mask(dv7, dv7_swap);
	m8 = _mm512_cmplt_epu64_mask(dv8, dv8_swap);

	// 'CC' for the swap between dist-2 lanes, non alternating
	dv1 = _mm512_mask_blend_epi64(m1 ^ 0xCC, dv1, dv1_swap);
	dv2 = _mm512_mask_blend_epi64(m2 ^ 0xCC, dv2, dv2_swap);
	dv3 = _mm512_mask_blend_epi64(m3 ^ 0xCC, dv3, dv3_swap);
	dv4 = _mm512_mask_blend_epi64(m4 ^ 0xCC, dv4, dv4_swap);
	dv5 = _mm512_mask_blend_epi64(m5 ^ 0xCC, dv5, dv5_swap);
	dv6 = _mm512_mask_blend_epi64(m6 ^ 0xCC, dv6, dv6_swap);
	dv7 = _mm512_mask_blend_epi64(m7 ^ 0xCC, dv7, dv7_swap);
	dv8 = _mm512_mask_blend_epi64(m8 ^ 0xCC, dv8, dv8_swap);

	// adjacent swaps
	dv1_swap = SWAP64(dv1);
	dv2_swap = SWAP64(dv2);
	dv3_swap = SWAP64(dv3);
	dv4_swap = SWAP64(dv4);
	dv5_swap = SWAP64(dv5);
	dv6_swap = SWAP64(dv6);
	dv7_swap = SWAP64(dv7);
	dv8_swap = SWAP64(dv8);
	
	m1 = _mm512_cmpgt_epu64_mask(dv1, dv1_swap);
	m2 = _mm512_cmpgt_epu64_mask(dv2, dv2_swap);
	m3 = _mm512_cmplt_epu64_mask(dv3, dv3_swap);
	m4 = _mm512_cmplt_epu64_mask(dv4, dv4_swap);
	m5 = _mm512_cmpgt_epu64_mask(dv5, dv5_swap);
	m6 = _mm512_cmpgt_epu64_mask(dv6, dv6_swap);
	m7 = _mm512_cmplt_epu64_mask(dv7, dv7_swap);
	m8 = _mm512_cmplt_epu64_mask(dv8, dv8_swap);

	// 'AA' for the swap between adjacent lanes
	dv1 = _mm512_mask_blend_epi64(m1 ^ 0xAA, dv1, dv1_swap);
	dv2 = _mm512_mask_blend_epi64(m2 ^ 0xAA, dv2, dv2_swap);
	dv3 = _mm512_mask_blend_epi64(m3 ^ 0xAA, dv3, dv3_swap);
	dv4 = _mm512_mask_blend_epi64(m4 ^ 0xAA, dv4, dv4_swap);
	dv5 = _mm512_mask_blend_epi64(m5 ^ 0xAA, dv5, dv5_swap);
	dv6 = _mm512_mask_blend_epi64(m6 ^ 0xAA, dv6, dv6_swap);
	dv7 = _mm512_mask_blend_epi64(m7 ^ 0xAA, dv7, dv7_swap);
	dv8 = _mm512_mask_blend_epi64(m8 ^ 0xAA, dv8, dv8_swap);


	// phase 5 : dist-32 alternating compares

	// distance 16 swaps - just compare the vecs
	m1 = _mm512_cmpgt_epu64_mask(dv1, dv3);
	m2 = _mm512_cmpgt_epu64_mask(dv2, dv4);
	m3 = _mm512_cmplt_epu64_mask(dv5, dv7);
	m4 = _mm512_cmplt_epu64_mask(dv6, dv8);

	t1 = _mm512_mask_blend_epi64(m1, dv1, dv3);
	t2 = _mm512_mask_blend_epi64(m2, dv2, dv4);
	t3 = _mm512_mask_blend_epi64(m3, dv5, dv7);
	t4 = _mm512_mask_blend_epi64(m4, dv6, dv8);
		
	dv3 = _mm512_mask_blend_epi64(m1, dv3, dv1);
	dv4 = _mm512_mask_blend_epi64(m2, dv4, dv2);
	dv7 = _mm512_mask_blend_epi64(m3, dv7, dv5);
	dv8 = _mm512_mask_blend_epi64(m4, dv8, dv6);
	dv1 = t1;
	dv2 = t2;
	dv5 = t3;
	dv6 = t4;

	// distance 8 swaps - just compare the vecs
	m1 = _mm512_cmpgt_epu64_mask(dv1, dv2);
	m2 = _mm512_cmpgt_epu64_mask(dv3, dv4);
	m3 = _mm512_cmplt_epu64_mask(dv5, dv6);
	m4 = _mm512_cmplt_epu64_mask(dv7, dv8);

	t1 = _mm512_mask_blend_epi64(m1, dv1, dv2);
	t2 = _mm512_mask_blend_epi64(m2, dv3, dv4);
	t3 = _mm512_mask_blend_epi64(m3, dv5, dv6);
	t4 = _mm512_mask_blend_epi64(m4, dv7, dv8);
	
	dv2 = _mm512_mask_blend_epi64(m1, dv2, dv1);
	dv4 = _mm512_mask_blend_epi64(m2, dv4, dv3);
	dv6 = _mm512_mask_blend_epi64(m3, dv6, dv5);
	dv8 = _mm512_mask_blend_epi64(m4, dv8, dv7);
	
	dv1 = t1;
	dv3 = t2;
	dv5 = t3;
	dv7 = t4;

	// distance 4 swaps
	dv1_swap = SWAP256(dv1);
	dv2_swap = SWAP256(dv2);
	dv3_swap = SWAP256(dv3);
	dv4_swap = SWAP256(dv4);
	dv5_swap = SWAP256(dv5);
	dv6_swap = SWAP256(dv6);
	dv7_swap = SWAP256(dv7);
	dv8_swap = SWAP256(dv8);
	
	m1 = _mm512_cmpgt_epu64_mask(dv1, dv1_swap);
	m2 = _mm512_cmpgt_epu64_mask(dv2, dv2_swap);
	m3 = _mm512_cmpgt_epu64_mask(dv3, dv3_swap);
	m4 = _mm512_cmpgt_epu64_mask(dv4, dv4_swap);
	m5 = _mm512_cmplt_epu64_mask(dv5, dv5_swap);
	m6 = _mm512_cmplt_epu64_mask(dv6, dv6_swap);
	m7 = _mm512_cmplt_epu64_mask(dv7, dv7_swap);
	m8 = _mm512_cmplt_epu64_mask(dv8, dv8_swap);

	// '0F' for the swap between dist-4 lanes, non alternating
	dv1 = _mm512_mask_blend_epi64(m1 ^ 0xF0, dv1, dv1_swap);
	dv2 = _mm512_mask_blend_epi64(m2 ^ 0xF0, dv2, dv2_swap);
	dv3 = _mm512_mask_blend_epi64(m3 ^ 0xF0, dv3, dv3_swap);
	dv4 = _mm512_mask_blend_epi64(m4 ^ 0xF0, dv4, dv4_swap);
	dv5 = _mm512_mask_blend_epi64(m5 ^ 0xF0, dv5, dv5_swap);
	dv6 = _mm512_mask_blend_epi64(m6 ^ 0xF0, dv6, dv6_swap);
	dv7 = _mm512_mask_blend_epi64(m7 ^ 0xF0, dv7, dv7_swap);
	dv8 = _mm512_mask_blend_epi64(m8 ^ 0xF0, dv8, dv8_swap);

	// distance 2 swaps
	dv1_swap = SWAP128(dv1);
	dv2_swap = SWAP128(dv2);
	dv3_swap = SWAP128(dv3);
	dv4_swap = SWAP128(dv4);
	dv5_swap = SWAP128(dv5);
	dv6_swap = SWAP128(dv6);
	dv7_swap = SWAP128(dv7);
	dv8_swap = SWAP128(dv8);
	
	m1 = _mm512_cmpgt_epu64_mask(dv1, dv1_swap);
	m2 = _mm512_cmpgt_epu64_mask(dv2, dv2_swap);
	m3 = _mm512_cmpgt_epu64_mask(dv3, dv3_swap);
	m4 = _mm512_cmpgt_epu64_mask(dv4, dv4_swap);
	m5 = _mm512_cmplt_epu64_mask(dv5, dv5_swap);
	m6 = _mm512_cmplt_epu64_mask(dv6, dv6_swap);
	m7 = _mm512_cmplt_epu64_mask(dv7, dv7_swap);
	m8 = _mm512_cmplt_epu64_mask(dv8, dv8_swap);

	// 'CC' for the swap between dist-2 lanes, non alternating
	dv1 = _mm512_mask_blend_epi64(m1 ^ 0xCC, dv1, dv1_swap);
	dv2 = _mm512_mask_blend_epi64(m2 ^ 0xCC, dv2, dv2_swap);
	dv3 = _mm512_mask_blend_epi64(m3 ^ 0xCC, dv3, dv3_swap);
	dv4 = _mm512_mask_blend_epi64(m4 ^ 0xCC, dv4, dv4_swap);
	dv5 = _mm512_mask_blend_epi64(m5 ^ 0xCC, dv5, dv5_swap);
	dv6 = _mm512_mask_blend_epi64(m6 ^ 0xCC, dv6, dv6_swap);
	dv7 = _mm512_mask_blend_epi64(m7 ^ 0xCC, dv7, dv7_swap);
	dv8 = _mm512_mask_blend_epi64(m8 ^ 0xCC, dv8, dv8_swap);

	// adjacent swaps
	dv1_swap = SWAP64(dv1);
	dv2_swap = SWAP64(dv2);
	dv3_swap = SWAP64(dv3);
	dv4_swap = SWAP64(dv4);
	dv5_swap = SWAP64(dv5);
	dv6_swap = SWAP64(dv6);
	dv7_swap = SWAP64(dv7);
	dv8_swap = SWAP64(dv8);
	
	m1 = _mm512_cmpgt_epu64_mask(dv1, dv1_swap);
	m2 = _mm512_cmpgt_epu64_mask(dv2, dv2_swap);
	m3 = _mm512_cmpgt_epu64_mask(dv3, dv3_swap);
	m4 = _mm512_cmpgt_epu64_mask(dv4, dv4_swap);
	m5 = _mm512_cmplt_epu64_mask(dv5, dv5_swap);
	m6 = _mm512_cmplt_epu64_mask(dv6, dv6_swap);
	m7 = _mm512_cmplt_epu64_mask(dv7, dv7_swap);
	m8 = _mm512_cmplt_epu64_mask(dv8, dv8_swap);

	// 'AA' for the swap between adjacent lanes
	dv1 = _mm512_mask_blend_epi64(m1 ^ 0xAA, dv1, dv1_swap);
	dv2 = _mm512_mask_blend_epi64(m2 ^ 0xAA, dv2, dv2_swap);
	dv3 = _mm512_mask_blend_epi64(m3 ^ 0xAA, dv3, dv3_swap);
	dv4 = _mm512_mask_blend_epi64(m4 ^ 0xAA, dv4, dv4_swap);
	dv5 = _mm512_mask_blend_epi64(m5 ^ 0xAA, dv5, dv5_swap);
	dv6 = _mm512_mask_blend_epi64(m6 ^ 0xAA, dv6, dv6_swap);
	dv7 = _mm512_mask_blend_epi64(m7 ^ 0xAA, dv7, dv7_swap);
	dv8 = _mm512_mask_blend_epi64(m8 ^ 0xAA, dv8, dv8_swap);


	if (dir == 1)
	{
		// distance 32 swaps - just compare the vecs
		m1 = _mm512_cmplt_epu64_mask(dv1, dv5);
		m2 = _mm512_cmplt_epu64_mask(dv2, dv6);
		m3 = _mm512_cmplt_epu64_mask(dv3, dv7);
		m4 = _mm512_cmplt_epu64_mask(dv4, dv8);

		t1 = _mm512_mask_blend_epi64(m1, dv1, dv5);
		t2 = _mm512_mask_blend_epi64(m2, dv2, dv6);
		t3 = _mm512_mask_blend_epi64(m3, dv3, dv7);
		t4 = _mm512_mask_blend_epi64(m4, dv4, dv8);
		dv5 = _mm512_mask_blend_epi64(m1, dv5, dv1);
		dv6 = _mm512_mask_blend_epi64(m2, dv6, dv2);
		dv7 = _mm512_mask_blend_epi64(m3, dv7, dv3);
		dv8 = _mm512_mask_blend_epi64(m4, dv8, dv4);
		dv1 = t1;
		dv2 = t2;
		dv3 = t3;
		dv4 = t4;
		
		// distance 16 swaps - just compare the vecs
		m1 = _mm512_cmplt_epu64_mask(dv1, dv3);
		m2 = _mm512_cmplt_epu64_mask(dv2, dv4);
		m3 = _mm512_cmplt_epu64_mask(dv5, dv7);
		m4 = _mm512_cmplt_epu64_mask(dv6, dv8);

		t1 = _mm512_mask_blend_epi64(m1, dv1, dv3);
		t2 = _mm512_mask_blend_epi64(m2, dv2, dv4);
		t3 = _mm512_mask_blend_epi64(m3, dv5, dv7);
		t4 = _mm512_mask_blend_epi64(m4, dv6, dv8);

		dv3 = _mm512_mask_blend_epi64(m1, dv3, dv1);
		dv4 = _mm512_mask_blend_epi64(m2, dv4, dv2);
		dv7 = _mm512_mask_blend_epi64(m3, dv7, dv5);
		dv8 = _mm512_mask_blend_epi64(m4, dv8, dv6);

		dv1 = t1;
		dv2 = t2;
		dv5 = t3;
		dv6 = t4;

		// distance 8 swaps - just compare the vecs
		m1 = _mm512_cmplt_epu64_mask(dv1, dv2);
		m2 = _mm512_cmplt_epu64_mask(dv3, dv4);
		m3 = _mm512_cmplt_epu64_mask(dv5, dv6);
		m4 = _mm512_cmplt_epu64_mask(dv7, dv8);

		t1 = _mm512_mask_blend_epi64(m1, dv1, dv2);
		t2 = _mm512_mask_blend_epi64(m2, dv3, dv4);
		t3 = _mm512_mask_blend_epi64(m3, dv5, dv6);
		t4 = _mm512_mask_blend_epi64(m4, dv7, dv8);
		
		dv2 = _mm512_mask_blend_epi64(m1, dv2, dv1);
		dv4 = _mm512_mask_blend_epi64(m2, dv4, dv3);
		dv6 = _mm512_mask_blend_epi64(m3, dv6, dv5);
		dv8 = _mm512_mask_blend_epi64(m4, dv8, dv7);
		
		dv1 = t1;
		dv3 = t2;
		dv5 = t3;
		dv7 = t4;

		// distance 4 swaps
		dv1_swap = SWAP256(dv1);
		dv2_swap = SWAP256(dv2);
		dv3_swap = SWAP256(dv3);
		dv4_swap = SWAP256(dv4);
		dv5_swap = SWAP256(dv5);
		dv6_swap = SWAP256(dv6);
		dv7_swap = SWAP256(dv7);
		dv8_swap = SWAP256(dv8);
		
		m1 = _mm512_cmplt_epu64_mask(dv1, dv1_swap);
		m2 = _mm512_cmplt_epu64_mask(dv2, dv2_swap);
		m3 = _mm512_cmplt_epu64_mask(dv3, dv3_swap);
		m4 = _mm512_cmplt_epu64_mask(dv4, dv4_swap);
		m5 = _mm512_cmplt_epu64_mask(dv5, dv5_swap);
		m6 = _mm512_cmplt_epu64_mask(dv6, dv6_swap);
		m7 = _mm512_cmplt_epu64_mask(dv7, dv7_swap);
		m8 = _mm512_cmplt_epu64_mask(dv8, dv8_swap);

		// 'F0' for the swap between dist-4 lanes, non alternating
		dv1 = _mm512_mask_blend_epi64(m1 ^ 0xF0, dv1, dv1_swap);
		dv2 = _mm512_mask_blend_epi64(m2 ^ 0xF0, dv2, dv2_swap);
		dv3 = _mm512_mask_blend_epi64(m3 ^ 0xF0, dv3, dv3_swap);
		dv4 = _mm512_mask_blend_epi64(m4 ^ 0xF0, dv4, dv4_swap);
		dv5 = _mm512_mask_blend_epi64(m5 ^ 0xF0, dv5, dv5_swap);
		dv6 = _mm512_mask_blend_epi64(m6 ^ 0xF0, dv6, dv6_swap);
		dv7 = _mm512_mask_blend_epi64(m7 ^ 0xF0, dv7, dv7_swap);
		dv8 = _mm512_mask_blend_epi64(m8 ^ 0xF0, dv8, dv8_swap);

		// distance 2 swaps
		dv1_swap = SWAP128(dv1);
		dv2_swap = SWAP128(dv2);
		dv3_swap = SWAP128(dv3);
		dv4_swap = SWAP128(dv4);
		dv5_swap = SWAP128(dv5);
		dv6_swap = SWAP128(dv6);
		dv7_swap = SWAP128(dv7);
		dv8_swap = SWAP128(dv8);
		
		m1 = _mm512_cmplt_epu64_mask(dv1, dv1_swap);
		m2 = _mm512_cmplt_epu64_mask(dv2, dv2_swap);
		m3 = _mm512_cmplt_epu64_mask(dv3, dv3_swap);
		m4 = _mm512_cmplt_epu64_mask(dv4, dv4_swap);
		m5 = _mm512_cmplt_epu64_mask(dv5, dv5_swap);
		m6 = _mm512_cmplt_epu64_mask(dv6, dv6_swap);
		m7 = _mm512_cmplt_epu64_mask(dv7, dv7_swap);
		m8 = _mm512_cmplt_epu64_mask(dv8, dv8_swap);

		// 'CC' for the swap between dist-2 lanes, non alternating
		dv1 = _mm512_mask_blend_epi64(m1 ^ 0xCC, dv1, dv1_swap);
		dv2 = _mm512_mask_blend_epi64(m2 ^ 0xCC, dv2, dv2_swap);
		dv3 = _mm512_mask_blend_epi64(m3 ^ 0xCC, dv3, dv3_swap);
		dv4 = _mm512_mask_blend_epi64(m4 ^ 0xCC, dv4, dv4_swap);
		dv5 = _mm512_mask_blend_epi64(m5 ^ 0xCC, dv5, dv5_swap);
		dv6 = _mm512_mask_blend_epi64(m6 ^ 0xCC, dv6, dv6_swap);
		dv7 = _mm512_mask_blend_epi64(m7 ^ 0xCC, dv7, dv7_swap);
		dv8 = _mm512_mask_blend_epi64(m8 ^ 0xCC, dv8, dv8_swap);

		// adjacent swaps
		dv1_swap = SWAP64(dv1);
		dv2_swap = SWAP64(dv2);
		dv3_swap = SWAP64(dv3);
		dv4_swap = SWAP64(dv4);
		dv5_swap = SWAP64(dv5);
		dv6_swap = SWAP64(dv6);
		dv7_swap = SWAP64(dv7);
		dv8_swap = SWAP64(dv8);
		
		m1 = _mm512_cmplt_epu64_mask(dv1, dv1_swap);
		m2 = _mm512_cmplt_epu64_mask(dv2, dv2_swap);
		m3 = _mm512_cmplt_epu64_mask(dv3, dv3_swap);
		m4 = _mm512_cmplt_epu64_mask(dv4, dv4_swap);
		m5 = _mm512_cmplt_epu64_mask(dv5, dv5_swap);
		m6 = _mm512_cmplt_epu64_mask(dv6, dv6_swap);
		m7 = _mm512_cmplt_epu64_mask(dv7, dv7_swap);
		m8 = _mm512_cmplt_epu64_mask(dv8, dv8_swap);

		// 'AA' for the swap between adjacent lanes
		dv1 = _mm512_mask_blend_epi64(m1 ^ 0xAA, dv1, dv1_swap);
		dv2 = _mm512_mask_blend_epi64(m2 ^ 0xAA, dv2, dv2_swap);
		dv3 = _mm512_mask_blend_epi64(m3 ^ 0xAA, dv3, dv3_swap);
		dv4 = _mm512_mask_blend_epi64(m4 ^ 0xAA, dv4, dv4_swap);
		dv5 = _mm512_mask_blend_epi64(m5 ^ 0xAA, dv5, dv5_swap);
		dv6 = _mm512_mask_blend_epi64(m6 ^ 0xAA, dv6, dv6_swap);
		dv7 = _mm512_mask_blend_epi64(m7 ^ 0xAA, dv7, dv7_swap);
		dv8 = _mm512_mask_blend_epi64(m8 ^ 0xAA, dv8, dv8_swap);

	}
	else
	{

		// distance 32 swaps - just compare the vecs
		m1 = _mm512_cmpgt_epu64_mask(dv1, dv5);
		m2 = _mm512_cmpgt_epu64_mask(dv2, dv6);
		m3 = _mm512_cmpgt_epu64_mask(dv3, dv7);
		m4 = _mm512_cmpgt_epu64_mask(dv4, dv8);

		t1 = _mm512_mask_blend_epi64(m1, dv1, dv5);
		t2 = _mm512_mask_blend_epi64(m2, dv2, dv6);
		t3 = _mm512_mask_blend_epi64(m3, dv3, dv7);
		t4 = _mm512_mask_blend_epi64(m4, dv4, dv8);
		dv5 = _mm512_mask_blend_epi64(m1, dv5, dv1);
		dv6 = _mm512_mask_blend_epi64(m2, dv6, dv2);
		dv7 = _mm512_mask_blend_epi64(m3, dv7, dv3);
		dv8 = _mm512_mask_blend_epi64(m4, dv8, dv4);
		dv1 = t1;
		dv2 = t2;
		dv3 = t3;
		dv4 = t4;
		
		// distance 16 swaps - just compare the vecs
		m1 = _mm512_cmpgt_epu64_mask(dv1, dv3);
		m2 = _mm512_cmpgt_epu64_mask(dv2, dv4);
		m3 = _mm512_cmpgt_epu64_mask(dv5, dv7);
		m4 = _mm512_cmpgt_epu64_mask(dv6, dv8);

		t1 = _mm512_mask_blend_epi64(m1, dv1, dv3);
		t2 = _mm512_mask_blend_epi64(m2, dv2, dv4);
		t3 = _mm512_mask_blend_epi64(m3, dv5, dv7);
		t4 = _mm512_mask_blend_epi64(m4, dv6, dv8);

		dv3 = _mm512_mask_blend_epi64(m1, dv3, dv1);
		dv4 = _mm512_mask_blend_epi64(m2, dv4, dv2);
		dv7 = _mm512_mask_blend_epi64(m3, dv7, dv5);
		dv8 = _mm512_mask_blend_epi64(m4, dv8, dv6);

		dv1 = t1;
		dv2 = t2;
		dv5 = t3;
		dv6 = t4;

		// distance 8 swaps - just compare the vecs
		m1 = _mm512_cmpgt_epu64_mask(dv1, dv2);
		m2 = _mm512_cmpgt_epu64_mask(dv3, dv4);
		m3 = _mm512_cmpgt_epu64_mask(dv5, dv6);
		m4 = _mm512_cmpgt_epu64_mask(dv7, dv8);

		t1 = _mm512_mask_blend_epi64(m1, dv1, dv2);
		t2 = _mm512_mask_blend_epi64(m2, dv3, dv4);
		t3 = _mm512_mask_blend_epi64(m3, dv5, dv6);
		t4 = _mm512_mask_blend_epi64(m4, dv7, dv8);
		
		dv2 = _mm512_mask_blend_epi64(m1, dv2, dv1);
		dv4 = _mm512_mask_blend_epi64(m2, dv4, dv3);
		dv6 = _mm512_mask_blend_epi64(m3, dv6, dv5);
		dv8 = _mm512_mask_blend_epi64(m4, dv8, dv7);
		
		dv1 = t1;
		dv3 = t2;
		dv5 = t3;
		dv7 = t4;

		// distance 4 swaps
		dv1_swap = SWAP256(dv1);
		dv2_swap = SWAP256(dv2);
		dv3_swap = SWAP256(dv3);
		dv4_swap = SWAP256(dv4);
		dv5_swap = SWAP256(dv5);
		dv6_swap = SWAP256(dv6);
		dv7_swap = SWAP256(dv7);
		dv8_swap = SWAP256(dv8);
		
		m1 = _mm512_cmpgt_epu64_mask(dv1, dv1_swap);
		m2 = _mm512_cmpgt_epu64_mask(dv2, dv2_swap);
		m3 = _mm512_cmpgt_epu64_mask(dv3, dv3_swap);
		m4 = _mm512_cmpgt_epu64_mask(dv4, dv4_swap);
		m5 = _mm512_cmpgt_epu64_mask(dv5, dv5_swap);
		m6 = _mm512_cmpgt_epu64_mask(dv6, dv6_swap);
		m7 = _mm512_cmpgt_epu64_mask(dv7, dv7_swap);
		m8 = _mm512_cmpgt_epu64_mask(dv8, dv8_swap);

		// 'F0' for the swap between dist-4 lanes, non alternating
		dv1 = _mm512_mask_blend_epi64(m1 ^ 0xF0, dv1, dv1_swap);
		dv2 = _mm512_mask_blend_epi64(m2 ^ 0xF0, dv2, dv2_swap);
		dv3 = _mm512_mask_blend_epi64(m3 ^ 0xF0, dv3, dv3_swap);
		dv4 = _mm512_mask_blend_epi64(m4 ^ 0xF0, dv4, dv4_swap);
		dv5 = _mm512_mask_blend_epi64(m5 ^ 0xF0, dv5, dv5_swap);
		dv6 = _mm512_mask_blend_epi64(m6 ^ 0xF0, dv6, dv6_swap);
		dv7 = _mm512_mask_blend_epi64(m7 ^ 0xF0, dv7, dv7_swap);
		dv8 = _mm512_mask_blend_epi64(m8 ^ 0xF0, dv8, dv8_swap);

		// distance 2 swaps
		dv1_swap = SWAP128(dv1);
		dv2_swap = SWAP128(dv2);
		dv3_swap = SWAP128(dv3);
		dv4_swap = SWAP128(dv4);
		dv5_swap = SWAP128(dv5);
		dv6_swap = SWAP128(dv6);
		dv7_swap = SWAP128(dv7);
		dv8_swap = SWAP128(dv8);
		
		m1 = _mm512_cmpgt_epu64_mask(dv1, dv1_swap);
		m2 = _mm512_cmpgt_epu64_mask(dv2, dv2_swap);
		m3 = _mm512_cmpgt_epu64_mask(dv3, dv3_swap);
		m4 = _mm512_cmpgt_epu64_mask(dv4, dv4_swap);
		m5 = _mm512_cmpgt_epu64_mask(dv5, dv5_swap);
		m6 = _mm512_cmpgt_epu64_mask(dv6, dv6_swap);
		m7 = _mm512_cmpgt_epu64_mask(dv7, dv7_swap);
		m8 = _mm512_cmpgt_epu64_mask(dv8, dv8_swap);

		// 'CC' for the swap between dist-2 lanes, non alternating
		dv1 = _mm512_mask_blend_epi64(m1 ^ 0xCC, dv1, dv1_swap);
		dv2 = _mm512_mask_blend_epi64(m2 ^ 0xCC, dv2, dv2_swap);
		dv3 = _mm512_mask_blend_epi64(m3 ^ 0xCC, dv3, dv3_swap);
		dv4 = _mm512_mask_blend_epi64(m4 ^ 0xCC, dv4, dv4_swap);
		dv5 = _mm512_mask_blend_epi64(m5 ^ 0xCC, dv5, dv5_swap);
		dv6 = _mm512_mask_blend_epi64(m6 ^ 0xCC, dv6, dv6_swap);
		dv7 = _mm512_mask_blend_epi64(m7 ^ 0xCC, dv7, dv7_swap);
		dv8 = _mm512_mask_blend_epi64(m8 ^ 0xCC, dv8, dv8_swap);

		// adjacent swaps
		dv1_swap = SWAP64(dv1);
		dv2_swap = SWAP64(dv2);
		dv3_swap = SWAP64(dv3);
		dv4_swap = SWAP64(dv4);
		dv5_swap = SWAP64(dv5);
		dv6_swap = SWAP64(dv6);
		dv7_swap = SWAP64(dv7);
		dv8_swap = SWAP64(dv8);
		
		m1 = _mm512_cmpgt_epu64_mask(dv1, dv1_swap);
		m2 = _mm512_cmpgt_epu64_mask(dv2, dv2_swap);
		m3 = _mm512_cmpgt_epu64_mask(dv3, dv3_swap);
		m4 = _mm512_cmpgt_epu64_mask(dv4, dv4_swap);
		m5 = _mm512_cmpgt_epu64_mask(dv5, dv5_swap);
		m6 = _mm512_cmpgt_epu64_mask(dv6, dv6_swap);
		m7 = _mm512_cmpgt_epu64_mask(dv7, dv7_swap);
		m8 = _mm512_cmpgt_epu64_mask(dv8, dv8_swap);

		// 'AA' for the swap between adjacent lanes
		dv1 = _mm512_mask_blend_epi64(m1 ^ 0xAA, dv1, dv1_swap);
		dv2 = _mm512_mask_blend_epi64(m2 ^ 0xAA, dv2, dv2_swap);
		dv3 = _mm512_mask_blend_epi64(m3 ^ 0xAA, dv3, dv3_swap);
		dv4 = _mm512_mask_blend_epi64(m4 ^ 0xAA, dv4, dv4_swap);
		dv5 = _mm512_mask_blend_epi64(m5 ^ 0xAA, dv5, dv5_swap);
		dv6 = _mm512_mask_blend_epi64(m6 ^ 0xAA, dv6, dv6_swap);
		dv7 = _mm512_mask_blend_epi64(m7 ^ 0xAA, dv7, dv7_swap);
		dv8 = _mm512_mask_blend_epi64(m8 ^ 0xAA, dv8, dv8_swap);
	}

	_mm512_store_si512(data, dv1);
	_mm512_store_si512(data + 8, dv2);
	_mm512_store_si512(data + 16, dv3);
	_mm512_store_si512(data + 24, dv4);
	_mm512_store_si512(data + 32, dv5);
	_mm512_store_si512(data + 40, dv6);
	_mm512_store_si512(data + 48, dv7);
	_mm512_store_si512(data + 56, dv8);

	return;
}

// on Epyc, max/min_epu64 is slower: the instruction is 
// latency 3, throughput 1 instead of 1/0.5 for epu32
void bitonic_sort_dir_64_minmax(uint64_t* data, int dir) 
{
	// sort 64 64-bit elements
	int i, j;

	__m512i t1;
	__m512i t2;
	__m512i t3;
	__m512i t4;
	__m512i t5;
	__m512i t6;
	__m512i t7;
	__m512i t8;
	__m512i dv1;
	__m512i dv2;
	__m512i dv3;
	__m512i dv4;
	__m512i dv5;
	__m512i dv6;
	__m512i dv7;
	__m512i dv8;
	__m512i dv1_swap;
	__m512i dv2_swap;
	__m512i dv3_swap;
	__m512i dv4_swap;
	__m512i dv5_swap;
	__m512i dv6_swap;
	__m512i dv7_swap;
	__m512i dv8_swap;
	__mmask16 m1;
	__mmask16 m2;
	__mmask16 m3;
	__mmask16 m4;
	__mmask16 m5;
	__mmask16 m6;
	__mmask16 m7;
	__mmask16 m8;

	dv1 = _mm512_load_si512(data);
	dv2 = _mm512_load_si512(data + 8);
	dv3 = _mm512_load_si512(data + 16);
	dv4 = _mm512_load_si512(data + 24);
	dv5 = _mm512_load_si512(data + 32);
	dv6 = _mm512_load_si512(data + 40);
	dv7 = _mm512_load_si512(data + 48);
	dv8 = _mm512_load_si512(data + 56);

	// phase 1: dist-2 alternating compares ('CC')
	
	// adjacent swaps, alternating compares
	dv1_swap = SWAP64(dv1);
	dv2_swap = SWAP64(dv2);
	dv3_swap = SWAP64(dv3);
	dv4_swap = SWAP64(dv4);
	dv5_swap = SWAP64(dv5);
	dv6_swap = SWAP64(dv6);
	dv7_swap = SWAP64(dv7);
	dv8_swap = SWAP64(dv8);
	
	// 'AA' for the swap between adjacent lanes ^ 'CC' to alternate gt vs. le --> 0x66
	MINMAX64(0x66, 0x99, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8)
	
	// phase 2: dist-4 alternating compares ('F0')
	
	// distance 2 swaps
	dv1_swap = SWAP128(dv1);
	dv2_swap = SWAP128(dv2);
	dv3_swap = SWAP128(dv3);
	dv4_swap = SWAP128(dv4);
	dv5_swap = SWAP128(dv5);
	dv6_swap = SWAP128(dv6);
	dv7_swap = SWAP128(dv7);
	dv8_swap = SWAP128(dv8);

	// 'CC' for the swap between dist-2 lanes ^ 'F0' to alternate gt vs. le --> 0x3C
	MINMAX64(0x3C, 0xC3, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8)

	// adjacent swaps
	dv1_swap = SWAP64(dv1);
	dv2_swap = SWAP64(dv2);
	dv3_swap = SWAP64(dv3);
	dv4_swap = SWAP64(dv4);
	dv5_swap = SWAP64(dv5);
	dv6_swap = SWAP64(dv6);
	dv7_swap = SWAP64(dv7);
	dv8_swap = SWAP64(dv8);

	// 'AA' for the swap between adjacent lanes ^ 'F0' to alternate gt vs. le --> 0x5A
	MINMAX64(0x5A, 0xA5, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// phase 3: dist-8 alternating compares ('FF00')

	// distance 4 swaps
	dv1_swap = SWAP256(dv1);
	dv2_swap = SWAP256(dv2);
	dv3_swap = SWAP256(dv3);
	dv4_swap = SWAP256(dv4);
	dv5_swap = SWAP256(dv5);
	dv6_swap = SWAP256(dv6);
	dv7_swap = SWAP256(dv7);
	dv8_swap = SWAP256(dv8);

	// 'F0F0' for the swap between dist-4 lanes ^ 'FF00' to alternate gt vs. le --> 0x0FF0
	MINMAX64_alt1(0xF0, 0x0f, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// distance 2 swaps
	dv1_swap = SWAP128(dv1);
	dv2_swap = SWAP128(dv2);
	dv3_swap = SWAP128(dv3);
	dv4_swap = SWAP128(dv4);
	dv5_swap = SWAP128(dv5);
	dv6_swap = SWAP128(dv6);
	dv7_swap = SWAP128(dv7);
	dv8_swap = SWAP128(dv8);

	// 'CCCC' for the swap between dist-2 lanes ^ 'FF00' to alternate gt vs. le --> 0x33CC
	MINMAX64_alt1(0xCC, 0x33, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// adjacent swaps
	dv1_swap = SWAP64(dv1);
	dv2_swap = SWAP64(dv2);
	dv3_swap = SWAP64(dv3);
	dv4_swap = SWAP64(dv4);
	dv5_swap = SWAP64(dv5);
	dv6_swap = SWAP64(dv6);
	dv7_swap = SWAP64(dv7);
	dv8_swap = SWAP64(dv8);

	// 'AA' for the swap between adjacent lanes ^ 'FF00' to alternate gt vs. le --> 0x55AA
	MINMAX64_alt1(0xAA, 0x55, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// phase 4: dist-16 alternating compares (alternating gt/lt)
	
	// distance 8 swaps (256 bits)
	t1 =  _mm512_min_epu64(dv1, dv2);
	t2 =  _mm512_max_epu64(dv3, dv4);
	t3 =  _mm512_min_epu64(dv5, dv6);
	t4 =  _mm512_max_epu64(dv7, dv8);
	dv2 = _mm512_max_epu64(dv1, dv2);
	dv4 = _mm512_min_epu64(dv3, dv4);
	dv6 = _mm512_max_epu64(dv5, dv6);
	dv8 = _mm512_min_epu64(dv7, dv8);
	dv1 = t1;
	dv3 = t2;
	dv5 = t3;
	dv7 = t4;


	// distance 4 swaps
	dv1_swap = SWAP256(dv1);
	dv2_swap = SWAP256(dv2);
	dv3_swap = SWAP256(dv3);
	dv4_swap = SWAP256(dv4);
	dv5_swap = SWAP256(dv5);
	dv6_swap = SWAP256(dv6);
	dv7_swap = SWAP256(dv7);
	dv8_swap = SWAP256(dv8);

	// 'F0' for the swap between dist-4 lanes, non alternating
	MINMAX64_alt2(0xF0, 0x0f, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// distance 2 swaps
	dv1_swap = SWAP128(dv1);
	dv2_swap = SWAP128(dv2);
	dv3_swap = SWAP128(dv3);
	dv4_swap = SWAP128(dv4);
	dv5_swap = SWAP128(dv5);
	dv6_swap = SWAP128(dv6);
	dv7_swap = SWAP128(dv7);
	dv8_swap = SWAP128(dv8);

	// 'CC' for the swap between dist-2 lanes, non alternating
	MINMAX64_alt2(0xCC, 0x33, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// adjacent swaps
	dv1_swap = SWAP64(dv1);
	dv2_swap = SWAP64(dv2);
	dv3_swap = SWAP64(dv3);
	dv4_swap = SWAP64(dv4);
	dv5_swap = SWAP64(dv5);
	dv6_swap = SWAP64(dv6);
	dv7_swap = SWAP64(dv7);
	dv8_swap = SWAP64(dv8);

	// 'AA' for the swap between adjacent lanes, non alternating
	MINMAX64_alt2(0xAA, 0x55, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// phase 5: dist-32 alternating compares (alternating gtgt/ltlt)

	// distance 16 swaps - just compare the vecs
	t1 =  _mm512_min_epu64(dv1, dv3);
	t2 =  _mm512_min_epu64(dv2, dv4);
	t3 =  _mm512_max_epu64(dv5, dv7);
	t4 =  _mm512_max_epu64(dv6, dv8);
	dv3 = _mm512_max_epu64(dv1, dv3);
	dv4 = _mm512_max_epu64(dv2, dv4);
	dv7 = _mm512_min_epu64(dv5, dv7);
	dv8 = _mm512_min_epu64(dv6, dv8);
	dv1 = t1;
	dv2 = t2;
	dv5 = t3;
	dv6 = t4;

	// distance 8 swaps (256 bits)
	t1 =  _mm512_min_epu64(dv1, dv2);
	t2 =  _mm512_min_epu64(dv3, dv4);
	t3 =  _mm512_max_epu64(dv5, dv6);
	t4 =  _mm512_max_epu64(dv7, dv8);
	dv2 = _mm512_max_epu64(dv1, dv2);
	dv4 = _mm512_max_epu64(dv3, dv4);
	dv6 = _mm512_min_epu64(dv5, dv6);
	dv8 = _mm512_min_epu64(dv7, dv8);
	dv1 = t1;
	dv3 = t2;
	dv5 = t3;
	dv7 = t4;

	// distance 4 swaps
	dv1_swap = SWAP256(dv1);
	dv2_swap = SWAP256(dv2);
	dv3_swap = SWAP256(dv3);
	dv4_swap = SWAP256(dv4);
	dv5_swap = SWAP256(dv5);
	dv6_swap = SWAP256(dv6);
	dv7_swap = SWAP256(dv7);
	dv8_swap = SWAP256(dv8);

	// 'F0' for the swap between dist-4 lanes, non alternating
	MINMAX64_alt4(0xF0, 0x0f, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// distance 2 swaps
	dv1_swap = SWAP128(dv1);
	dv2_swap = SWAP128(dv2);
	dv3_swap = SWAP128(dv3);
	dv4_swap = SWAP128(dv4);
	dv5_swap = SWAP128(dv5);
	dv6_swap = SWAP128(dv6);
	dv7_swap = SWAP128(dv7);
	dv8_swap = SWAP128(dv8);

	// 'CC' for the swap between dist-2 lanes, non alternating
	MINMAX64_alt4(0xCC, 0x33, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	// adjacent swaps
	dv1_swap = SWAP64(dv1);
	dv2_swap = SWAP64(dv2);
	dv3_swap = SWAP64(dv3);
	dv4_swap = SWAP64(dv4);
	dv5_swap = SWAP64(dv5);
	dv6_swap = SWAP64(dv6);
	dv7_swap = SWAP64(dv7);
	dv8_swap = SWAP64(dv8);

	// 'AA' for the swap between adjacent lanes
	MINMAX64_alt4(0xAA, 0x55, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
	


	// phase 6: merge (all same compare 'dir')
	if (dir == 1)
	{	
		// distance 32 swaps - just compare the vecs
		t1 =  _mm512_max_epu64(dv1, dv5);
		t2 =  _mm512_max_epu64(dv2, dv6);
		t3 =  _mm512_max_epu64(dv3, dv7);
		t4 =  _mm512_max_epu64(dv4, dv8);
		dv5 = _mm512_min_epu64(dv1, dv5);
		dv6 = _mm512_min_epu64(dv2, dv6);
		dv7 = _mm512_min_epu64(dv3, dv7);
		dv8 = _mm512_min_epu64(dv4, dv8);
		dv1 = t1;
		dv2 = t2;
		dv3 = t3;
		dv4 = t4;

		// distance 16 swaps - just compare the vecs
		t1 =  _mm512_max_epu64(dv1, dv3);
		t2 =  _mm512_max_epu64(dv2, dv4);
		t3 =  _mm512_max_epu64(dv5, dv7);
		t4 =  _mm512_max_epu64(dv6, dv8);
		dv3 = _mm512_min_epu64(dv1, dv3);
		dv4 = _mm512_min_epu64(dv2, dv4);
		dv7 = _mm512_min_epu64(dv5, dv7);
		dv8 = _mm512_min_epu64(dv6, dv8);
		dv1 = t1;
		dv2 = t2;
		dv5 = t3;
		dv6 = t4;
		
		// distance 8 swaps - 256 bits
		t1 =  _mm512_max_epu64(dv1, dv2);
		t2 =  _mm512_max_epu64(dv3, dv4);
		t3 =  _mm512_max_epu64(dv5, dv6);
		t4 =  _mm512_max_epu64(dv7, dv8);
		dv2 = _mm512_min_epu64(dv1, dv2);
		dv4 = _mm512_min_epu64(dv3, dv4);
		dv6 = _mm512_min_epu64(dv5, dv6);
		dv8 = _mm512_min_epu64(dv7, dv8);
		dv1 = t1;
		dv3 = t2;
		dv5 = t3;
		dv7 = t4;
		
		// distance 4 swaps
		dv1_swap = SWAP256(dv1);
		dv2_swap = SWAP256(dv2);
		dv3_swap = SWAP256(dv3);
		dv4_swap = SWAP256(dv4);
		dv5_swap = SWAP256(dv5);
		dv6_swap = SWAP256(dv6);
		dv7_swap = SWAP256(dv7);
		dv8_swap = SWAP256(dv8);

		// 'F0' for the swap between dist-4 lanes, non alternating
		MINMAX64(0x0f, 0xF0, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		// distance 2 swaps
		dv1_swap = SWAP128(dv1);
		dv2_swap = SWAP128(dv2);
		dv3_swap = SWAP128(dv3);
		dv4_swap = SWAP128(dv4);
		dv5_swap = SWAP128(dv5);
		dv6_swap = SWAP128(dv6);
		dv7_swap = SWAP128(dv7);
		dv8_swap = SWAP128(dv8);

		// 'CC' for the swap between dist-2 lanes, non alternating
		MINMAX64(0x33, 0xCC, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		// adjacent swaps
		dv1_swap = SWAP64(dv1);
		dv2_swap = SWAP64(dv2);
		dv3_swap = SWAP64(dv3);
		dv4_swap = SWAP64(dv4);
		dv5_swap = SWAP64(dv5);
		dv6_swap = SWAP64(dv6);
		dv7_swap = SWAP64(dv7);
		dv8_swap = SWAP64(dv8);

		// 'AA' for the swap between adjacent lanes
		MINMAX64(0x55, 0xaa, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	}
	else
	{	
		// distance 32 swaps - just compare the vecs
		t1 =  _mm512_min_epu64(dv1, dv5);
		t2 =  _mm512_min_epu64(dv2, dv6);
		t3 =  _mm512_min_epu64(dv3, dv7);
		t4 =  _mm512_min_epu64(dv4, dv8);
		dv5 = _mm512_max_epu64(dv1, dv5);
		dv6 = _mm512_max_epu64(dv2, dv6);
		dv7 = _mm512_max_epu64(dv3, dv7);
		dv8 = _mm512_max_epu64(dv4, dv8);
		dv1 = t1;
		dv2 = t2;
		dv3 = t3;
		dv4 = t4;

		// distance 16 swaps - just compare the vecs
		t1 =  _mm512_min_epu64(dv1, dv3);
		t2 =  _mm512_min_epu64(dv2, dv4);
		t3 =  _mm512_min_epu64(dv5, dv7);
		t4 =  _mm512_min_epu64(dv6, dv8);
		dv3 = _mm512_max_epu64(dv1, dv3);
		dv4 = _mm512_max_epu64(dv2, dv4);
		dv7 = _mm512_max_epu64(dv5, dv7);
		dv8 = _mm512_max_epu64(dv6, dv8);
		dv1 = t1;
		dv2 = t2;
		dv5 = t3;
		dv6 = t4;

		// distance 8 swaps - 256 bits
		t1 =  _mm512_min_epu64(dv1, dv2);
		t2 =  _mm512_min_epu64(dv3, dv4);
		t3 =  _mm512_min_epu64(dv5, dv6);
		t4 =  _mm512_min_epu64(dv7, dv8);
		dv2 = _mm512_max_epu64(dv1, dv2);
		dv4 = _mm512_max_epu64(dv3, dv4);
		dv6 = _mm512_max_epu64(dv5, dv6);
		dv8 = _mm512_max_epu64(dv7, dv8);
		dv1 = t1;
		dv3 = t2;
		dv5 = t3;
		dv7 = t4;
		
		// distance 4 swaps
		dv1_swap = SWAP256(dv1);
		dv2_swap = SWAP256(dv2);
		dv3_swap = SWAP256(dv3);
		dv4_swap = SWAP256(dv4);
		dv5_swap = SWAP256(dv5);
		dv6_swap = SWAP256(dv6);
		dv7_swap = SWAP256(dv7);
		dv8_swap = SWAP256(dv8);

		// '0F' for the swap between dist-4 lanes, non alternating
		MINMAX64(0xf0, 0x0f, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

		// distance 2 swaps
		dv1_swap = SWAP128(dv1);
		dv2_swap = SWAP128(dv2);
		dv3_swap = SWAP128(dv3);
		dv4_swap = SWAP128(dv4);
		dv5_swap = SWAP128(dv5);
		dv6_swap = SWAP128(dv6);
		dv7_swap = SWAP128(dv7);
		dv8_swap = SWAP128(dv8);

		// 'CC' for the swap between dist-2 lanes, non alternating
		MINMAX64(0xCC, 0x33, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);
		
		// adjacent swaps
		dv1_swap = SWAP64(dv1);
		dv2_swap = SWAP64(dv2);
		dv3_swap = SWAP64(dv3);
		dv4_swap = SWAP64(dv4);
		dv5_swap = SWAP64(dv5);
		dv6_swap = SWAP64(dv6);
		dv7_swap = SWAP64(dv7);
		dv8_swap = SWAP64(dv8);

		// 'AA' for the swap between adjacent lanes
		MINMAX64(0xaa, 0x55, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8);

	}

	_mm512_store_si512(data, dv1);
	_mm512_store_si512(data + 8, dv2);
	_mm512_store_si512(data + 16, dv3);
	_mm512_store_si512(data + 24, dv4);
	_mm512_store_si512(data + 32, dv5);
	_mm512_store_si512(data + 40, dv6);
	_mm512_store_si512(data + 48, dv7);
	_mm512_store_si512(data + 56, dv8);
	
	return;
}

void bitonic_merge(uint64_t *data, uint32_t sz, int dir)
{
	if (sz <= 64)
	{
		// base case: do the hardcoded 64-element sort
		bitonic_merge_dir_64(data, dir);
		return;
	}
	
	// half-size cmp/swap
	int i;

	// we have sz/2 swaps to do at a stride of sz/2.
	// the number of swaps will be divisible by 64 because
	// sz is at least 128 (basecase is 64).
	// can batch up the swaps:

	if (dir == 1)
	{
		// 128-element merge passes at a stride of sz/2
		for (i = 0; i < sz / 128; i++)
		{
			__m512i t1;
			__m512i t2;
			__m512i t3;
			__m512i t4;
			__m512i t5;
			__m512i t6;
			__m512i t7;
			__m512i t8;
			__m512i dv1;
			__m512i dv2;
			__m512i dv3;
			__m512i dv4;
			__m512i dv5;
			__m512i dv6;
			__m512i dv7;
			__m512i dv8;
			__m512i dv9;
			__m512i dv10;
			__m512i dv11;
			__m512i dv12;
			__m512i dv13;
			__m512i dv14;
			__m512i dv15;
			__m512i dv16;

			__mmask8 m1;
			__mmask8 m2;
			__mmask8 m3;
			__mmask8 m4;
			__mmask8 m5;
			__mmask8 m6;
			__mmask8 m7;
			__mmask8 m8;

			// swap 64 elements at starting offset i * 64
			dv1 = _mm512_load_si512(data + i * 64 + 0);
			dv2 = _mm512_load_si512(data + i * 64 + 8);
			dv3 = _mm512_load_si512(data + i * 64 + 16);
			dv4 = _mm512_load_si512(data + i * 64 + 24);
			dv5 = _mm512_load_si512(data + i * 64 + 32);
			dv6 = _mm512_load_si512(data + i * 64 + 40);
			dv7 = _mm512_load_si512(data + i * 64 + 48);
			dv8 = _mm512_load_si512(data + i * 64 + 56);
			
			// with 64 elements at offset i * 64 + stride (sz/2)
			dv9 = _mm512_load_si512(data  + i * 64 + sz/2 + 0);
			dv10 = _mm512_load_si512(data + i * 64 + sz/2 + 8);
			dv11 = _mm512_load_si512(data + i * 64 + sz/2 + 16);
			dv12 = _mm512_load_si512(data + i * 64 + sz/2 + 24);
			dv13 = _mm512_load_si512(data + i * 64 + sz/2 + 32);
			dv14 = _mm512_load_si512(data + i * 64 + sz/2 + 40);
			dv15 = _mm512_load_si512(data + i * 64 + sz/2 + 48);
			dv16 = _mm512_load_si512(data + i * 64 + sz/2 + 56);

				m1 = _mm512_cmplt_epu64_mask(dv1, dv9);
				m2 = _mm512_cmplt_epu64_mask(dv2, dv10);
				m3 = _mm512_cmplt_epu64_mask(dv3, dv11);
				m4 = _mm512_cmplt_epu64_mask(dv4, dv12);
				m5 = _mm512_cmplt_epu64_mask(dv5, dv13);
				m6 = _mm512_cmplt_epu64_mask(dv6, dv14);
				m7 = _mm512_cmplt_epu64_mask(dv7, dv15);
				m8 = _mm512_cmplt_epu64_mask(dv8, dv16);


			t1 = _mm512_mask_blend_epi64(m1, dv1, dv9);
			t2 = _mm512_mask_blend_epi64(m2, dv2, dv10);
			t3 = _mm512_mask_blend_epi64(m3, dv3, dv11);
			t4 = _mm512_mask_blend_epi64(m4, dv4, dv12);
			t5 = _mm512_mask_blend_epi64(m5, dv5, dv13);
			t6 = _mm512_mask_blend_epi64(m6, dv6, dv14);
			t7 = _mm512_mask_blend_epi64(m7, dv7, dv15);
			t8 = _mm512_mask_blend_epi64(m8, dv8, dv16);
			
			dv9 = _mm512_mask_blend_epi64(m1, dv9, dv1);
			dv10 = _mm512_mask_blend_epi64(m2, dv10, dv2);
			dv11 = _mm512_mask_blend_epi64(m3, dv11, dv3);
			dv12 = _mm512_mask_blend_epi64(m4, dv12, dv4);
			dv13 = _mm512_mask_blend_epi64(m5, dv13, dv5);
			dv14 = _mm512_mask_blend_epi64(m6, dv14, dv6);
			dv15 = _mm512_mask_blend_epi64(m7, dv15, dv7);
			dv16 = _mm512_mask_blend_epi64(m8, dv16, dv8);

			dv1 = t1;
			dv2 = t2;
			dv3 = t3;
			dv4 = t4;
			dv5 = t5;
			dv6 = t6;
			dv7 = t7;
			dv8 = t8;
			
			_mm512_store_si512(data + i * 64 + 0, dv1);
			_mm512_store_si512(data + i * 64 + 8, dv2);
			_mm512_store_si512(data + i * 64 + 16, dv3);
			_mm512_store_si512(data + i * 64 + 24, dv4);
			_mm512_store_si512(data + i * 64 + 32, dv5);
			_mm512_store_si512(data + i * 64 + 40, dv6);
			_mm512_store_si512(data + i * 64 + 48, dv7);
			_mm512_store_si512(data + i * 64 + 56, dv8);
			
			_mm512_store_si512(data + i * 64 + sz/2 + 0,  dv9);
			_mm512_store_si512(data + i * 64 + sz/2 + 8,  dv10);
			_mm512_store_si512(data + i * 64 + sz/2 + 16, dv11);
			_mm512_store_si512(data + i * 64 + sz/2 + 24, dv12);
			_mm512_store_si512(data + i * 64 + sz/2 + 32, dv13);
			_mm512_store_si512(data + i * 64 + sz/2 + 40, dv14);
			_mm512_store_si512(data + i * 64 + sz/2 + 48, dv15);
			_mm512_store_si512(data + i * 64 + sz/2 + 56, dv16);

		}
	}
	else
	{
		// 128-element merge passes at a stride of sz/2
		for (i = 0; i < sz / 128; i++)
		{
			__m512i t1;
			__m512i t2;
			__m512i t3;
			__m512i t4;
			__m512i t5;
			__m512i t6;
			__m512i t7;
			__m512i t8;
			__m512i dv1;
			__m512i dv2;
			__m512i dv3;
			__m512i dv4;
			__m512i dv5;
			__m512i dv6;
			__m512i dv7;
			__m512i dv8;
			__m512i dv9;
			__m512i dv10;
			__m512i dv11;
			__m512i dv12;
			__m512i dv13;
			__m512i dv14;
			__m512i dv15;
			__m512i dv16;

			__mmask8 m1;
			__mmask8 m2;
			__mmask8 m3;
			__mmask8 m4;
			__mmask8 m5;
			__mmask8 m6;
			__mmask8 m7;
			__mmask8 m8;

			// swap 64 elements at starting offset i * 64
			dv1 = _mm512_load_si512(data + i * 64 + 0);
			dv2 = _mm512_load_si512(data + i * 64 + 8);
			dv3 = _mm512_load_si512(data + i * 64 + 16);
			dv4 = _mm512_load_si512(data + i * 64 + 24);
			dv5 = _mm512_load_si512(data + i * 64 + 32);
			dv6 = _mm512_load_si512(data + i * 64 + 40);
			dv7 = _mm512_load_si512(data + i * 64 + 48);
			dv8 = _mm512_load_si512(data + i * 64 + 56);
		
			// with 64 elements at offset i * 64 + stride (sz/2)
			dv9 = _mm512_load_si512(data  + i * 64 + sz/2 + 0);
			dv10 = _mm512_load_si512(data + i * 64 + sz/2 + 8);
			dv11 = _mm512_load_si512(data + i * 64 + sz/2 + 16);
			dv12 = _mm512_load_si512(data + i * 64 + sz/2 + 24);
			dv13 = _mm512_load_si512(data + i * 64 + sz/2 + 32);
			dv14 = _mm512_load_si512(data + i * 64 + sz/2 + 40);
			dv15 = _mm512_load_si512(data + i * 64 + sz/2 + 48);
			dv16 = _mm512_load_si512(data + i * 64 + sz/2 + 56);

				m1 = _mm512_cmpgt_epu64_mask(dv1, dv9);
				m2 = _mm512_cmpgt_epu64_mask(dv2, dv10);
				m3 = _mm512_cmpgt_epu64_mask(dv3, dv11);
				m4 = _mm512_cmpgt_epu64_mask(dv4, dv12);
				m5 = _mm512_cmpgt_epu64_mask(dv5, dv13);
				m6 = _mm512_cmpgt_epu64_mask(dv6, dv14);
				m7 = _mm512_cmpgt_epu64_mask(dv7, dv15);
				m8 = _mm512_cmpgt_epu64_mask(dv8, dv16);

			t1 = _mm512_mask_blend_epi64(m1, dv1, dv9);
			t2 = _mm512_mask_blend_epi64(m2, dv2, dv10);
			t3 = _mm512_mask_blend_epi64(m3, dv3, dv11);
			t4 = _mm512_mask_blend_epi64(m4, dv4, dv12);
			t5 = _mm512_mask_blend_epi64(m5, dv5, dv13);
			t6 = _mm512_mask_blend_epi64(m6, dv6, dv14);
			t7 = _mm512_mask_blend_epi64(m7, dv7, dv15);
			t8 = _mm512_mask_blend_epi64(m8, dv8, dv16);
			
			dv9 = _mm512_mask_blend_epi64(m1, dv9, dv1);
			dv10 = _mm512_mask_blend_epi64(m2, dv10, dv2);
			dv11 = _mm512_mask_blend_epi64(m3, dv11, dv3);
			dv12 = _mm512_mask_blend_epi64(m4, dv12, dv4);
			dv13 = _mm512_mask_blend_epi64(m5, dv13, dv5);
			dv14 = _mm512_mask_blend_epi64(m6, dv14, dv6);
			dv15 = _mm512_mask_blend_epi64(m7, dv15, dv7);
			dv16 = _mm512_mask_blend_epi64(m8, dv16, dv8);

			dv1 = t1;
			dv2 = t2;
			dv3 = t3;
			dv4 = t4;
			dv5 = t5;
			dv6 = t6;
			dv7 = t7;
			dv8 = t8;
			
			_mm512_store_si512(data + i * 64 + 0, dv1);
			_mm512_store_si512(data + i * 64 + 8, dv2);
			_mm512_store_si512(data + i * 64 + 16, dv3);
			_mm512_store_si512(data + i * 64 + 24, dv4);
			_mm512_store_si512(data + i * 64 + 32, dv5);
			_mm512_store_si512(data + i * 64 + 40, dv6);
			_mm512_store_si512(data + i * 64 + 48, dv7);
			_mm512_store_si512(data + i * 64 + 56, dv8);
			
			_mm512_store_si512(data + i * 64 + sz/2 + 0,  dv9);
			_mm512_store_si512(data + i * 64 + sz/2 + 8,  dv10);
			_mm512_store_si512(data + i * 64 + sz/2 + 16, dv11);
			_mm512_store_si512(data + i * 64 + sz/2 + 24, dv12);
			_mm512_store_si512(data + i * 64 + sz/2 + 32, dv13);
			_mm512_store_si512(data + i * 64 + sz/2 + 40, dv14);
			_mm512_store_si512(data + i * 64 + sz/2 + 48, dv15);
			_mm512_store_si512(data + i * 64 + sz/2 + 56, dv16);

		}
	}

	// two parallel half-size merges
	bitonic_merge(data, sz / 2, dir);
	bitonic_merge(data + sz / 2, sz / 2, dir);
}

void L1sort(uint64_t *data, int dir)
{
	// L1 is assumed 32k bytes. At 8 bytes per element
	// and 64 elements per sort there are 64 passes of
	// the base case over L1.
	int j;
	for (j = 0; j < 64; j++) {
		// alternating up/down sorts so we can finish using merge only
		bitonic_sort_dir_64(data + j * 64, j & 1);
	}
	
	uint32_t bitonic_sort_size = 128;
	while (bitonic_sort_size < 4096)
	{
		for (j = 0; j < 4096 / bitonic_sort_size; j++) {
			// up/down merges of the previous up/down sorts, output is up/down sorted
			bitonic_merge(&data[j * bitonic_sort_size], bitonic_sort_size, j & 1);
		}
		bitonic_sort_size *= 2;
	}
	
	// final merge in the specified direction
	bitonic_merge(data, 4096, dir);
	
	return;
}
		
//#define non_recursive
void bitonic_sort(uint64_t *data, uint32_t sz, int dir)
{
	if (sz == 64)
	{
		// base case: do the hardcoded 64-element sort
		bitonic_sort_dir_64(data, dir);
		return;
	}
	
#ifdef non_recursive

	if (sz < 4096)
	{	
		// two half-size bitonic sorts,
		// with opposite directions.
		bitonic_sort(data, sz / 2, 0);
		bitonic_sort(data + sz / 2, sz / 2, 1);
		
		// merge in the specified direction
		bitonic_merge(data, sz, dir);
	}
	else
	{
		// to make this more cache-friendly, get L1-sized chunks up/down sorted
		// and then merge those together.  L1 is typically 4096 uint64_t's (32k bytes).
		int j;
		for (j = 0; j < sz / 4096; j++) {
			// alternating up/down sorts so we can finish using merge only
			L1sort(data + j * 4096, j & 1);
		}
		
		uint32_t bitonic_sort_size = 8192;
		while (bitonic_sort_size < sz)
		{
			for (j = 0; j < sz / bitonic_sort_size; j++) {
				// up/down merges of the previous up/down sorts, output is up/down sorted
				bitonic_merge(&data[j * bitonic_sort_size], bitonic_sort_size, j & 1);
			}
			bitonic_sort_size *= 2;
		}
		
		// final merge in the specified direction
		bitonic_merge(data, sz, dir);
		return;
	}
#else
	
	// two half-size bitonic sorts,
	// with opposite directions.
	bitonic_sort(data, sz / 2, 0);
	bitonic_sort(data + sz / 2, sz / 2, 1);
	
	// merge in the specified direction
	bitonic_merge(data, sz, dir);
#endif

	return;
}

void sort(uint64_t *data, uint32_t sz, int dir)
{
	// top level sort dealing with two things:
	// 1) the bitonic sort function requires the data array
	// to be aligned and
	// 2) the bitonic sort works on a power-of-2 sized array
	// here we make sure we feed the bitonic sort function
	// an array that satifies those requirements.
	int is_aligned = (((uint64_t)data & 0x3full) == 0);
	if (is_aligned && ((sz & (sz - 1)) == 0))
	{
		// meets both requirements as-is
		bitonic_sort(data, sz, dir);
		return;
	}
	
	// otherwise we need to copy to an aligned buffer and/or 
	// change the buffer size
	uint32_t new_sz = sz;
	if ((sz & (sz - 1)) > 0)
	{
		new_sz = next_power_2(sz);
	}
		
	uint64_t *adata;
	
	if (!is_aligned)
	{
		adata = (uint64_t*)aligned_malloc(new_sz * sizeof(uint64_t), 64);
		memcpy(adata, data, sz * sizeof(uint64_t));
	}
	else
	{
		adata = data;
	}
	
	if ((new_sz - sz) > 0)
	{
		if (dir == 0)
			memset(adata + sz, 0xff, (new_sz - sz) * sizeof(uint64));
		else
			memset(adata + sz, 0, (new_sz - sz) * sizeof(uint64));
	}
	
	bitonic_sort(adata, new_sz, dir);
	
	if (!is_aligned)
	{
		memcpy(data, adata, sz * sizeof(uint64_t));
		aligned_free(adata);
	}
	
	return;
}

void bitonic_merge16(uint16_t *data, uint32_t sz, int dir)
{
	if (sz <= 64)
	{
		// base case: do the hardcoded 128-element sort
		bitonic_merge16_dir_64(data, dir);
		return;
	}
	else if (sz <= 128)
	{
		// base case: do the hardcoded 128-element sort
		bitonic_merge16_dir_128(data, dir);
		return;
	}
	else if (sz <= 256)
	{
		// base case: do the hardcoded 128-element sort
		bitonic_merge16_dir_256(data, dir);
		return;
	}
	
	// half-size cmp/swap
	int i;

	// we have sz/2 swaps to do at a stride of sz/2.
	// the number of swaps will be divisible by 256 because
	// sz is at least 512 (basecase is 256).
	// can batch up the swaps:
	__m512i t1;
	__m512i t2;
	__m512i t3;
	__m512i t4;
	__m512i t5;
	__m512i t6;
	__m512i t7;
	__m512i t8;
	__m512i dv1;
	__m512i dv2;
	__m512i dv3;
	__m512i dv4;
	__m512i dv5;
	__m512i dv6;
	__m512i dv7;
	__m512i dv8;
	__m512i dv9;
	__m512i dv10;
	__m512i dv11;
	__m512i dv12;
	__m512i dv13;
	__m512i dv14;
	__m512i dv15;
	__m512i dv16;

			
	if (dir == 1)
	{
		// 512-element merge passes at a stride of sz/2
		for (i = 0; i < sz / 512; i++)
		{
			// swap 128 elements at starting offset i * 128
			dv1 = _mm512_load_si512(data + i * 256 +  0);
			dv2 = _mm512_load_si512(data + i * 256 + 32);
			dv3 = _mm512_load_si512(data + i * 256 + 64);
			dv4 = _mm512_load_si512(data + i * 256 + 96);
			dv5 = _mm512_load_si512(data + i * 256 + 128);
			dv6 = _mm512_load_si512(data + i * 256 + 160);
			dv7 = _mm512_load_si512(data + i * 256 + 192);
			dv8 = _mm512_load_si512(data + i * 256 + 224);
			
			// with 128 elements at offset i * 128 + stride (sz/2)
			dv9  = _mm512_load_si512(data + i * 256 + sz/2 +  0);
			dv10 = _mm512_load_si512(data + i * 256 + sz/2 + 32);
			dv11 = _mm512_load_si512(data + i * 256 + sz/2 + 64);
			dv12 = _mm512_load_si512(data + i * 256 + sz/2 + 96);
			dv13 = _mm512_load_si512(data + i * 256 + sz/2 + 128);
			dv14 = _mm512_load_si512(data + i * 256 + sz/2 + 160);
			dv15 = _mm512_load_si512(data + i * 256 + sz/2 + 192);
			dv16 = _mm512_load_si512(data + i * 256 + sz/2 + 224);


			t1 =  _mm512_max_epu16(dv1, dv9 );
			t2 =  _mm512_max_epu16(dv2, dv10);
			t3 =  _mm512_max_epu16(dv3, dv11);
			t4 =  _mm512_max_epu16(dv4, dv12);
			dv9  = _mm512_min_epu16(dv1, dv9 );
			dv10 = _mm512_min_epu16(dv2, dv10);
			dv11 = _mm512_min_epu16(dv3, dv11);
			dv12 = _mm512_min_epu16(dv4, dv12);
			dv1 = t1;
			dv2 = t2;
			dv3 = t3;
			dv4 = t4;
			
			t1 =  _mm512_max_epu16(dv5, dv13);
			t2 =  _mm512_max_epu16(dv6, dv14);
			t3 =  _mm512_max_epu16(dv7, dv15);
			t4 =  _mm512_max_epu16(dv8, dv16);
			dv13 = _mm512_min_epu16(dv5, dv13);
			dv14 = _mm512_min_epu16(dv6, dv14);
			dv15 = _mm512_min_epu16(dv7, dv15);
			dv16 = _mm512_min_epu16(dv8, dv16);
			dv5 = t1;
			dv6 = t2;
			dv7 = t3;
			dv8 = t4;
			
			_mm512_store_si512(data + i * 256 +  0, dv1);
			_mm512_store_si512(data + i * 256 + 32, dv2);
			_mm512_store_si512(data + i * 256 + 64, dv3);
			_mm512_store_si512(data + i * 256 + 96, dv4);
			_mm512_store_si512(data + i * 256 + 128, dv5);
			_mm512_store_si512(data + i * 256 + 160, dv6);
			_mm512_store_si512(data + i * 256 + 192, dv7);
			_mm512_store_si512(data + i * 256 + 224, dv8);
			
			_mm512_store_si512(data + i * 256 + sz/2 +  0, dv9);
			_mm512_store_si512(data + i * 256 + sz/2 + 32, dv10);
			_mm512_store_si512(data + i * 256 + sz/2 + 64, dv11);
			_mm512_store_si512(data + i * 256 + sz/2 + 96, dv12);
			_mm512_store_si512(data + i * 256 + sz/2 + 128, dv13);
			_mm512_store_si512(data + i * 256 + sz/2 + 160, dv14);
			_mm512_store_si512(data + i * 256 + sz/2 + 192, dv15);
			_mm512_store_si512(data + i * 256 + sz/2 + 224, dv16);
		}
	}
	else
	{
		// 256-element merge passes at a stride of sz/2
		for (i = 0; i < sz / 512; i++)
		{
			// swap 128 elements at starting offset i * 128
			dv1 = _mm512_load_si512(data + i * 256 +  0);
			dv2 = _mm512_load_si512(data + i * 256 + 32);
			dv3 = _mm512_load_si512(data + i * 256 + 64);
			dv4 = _mm512_load_si512(data + i * 256 + 96);
			dv5 = _mm512_load_si512(data + i * 256 + 128);
			dv6 = _mm512_load_si512(data + i * 256 + 160);
			dv7 = _mm512_load_si512(data + i * 256 + 192);
			dv8 = _mm512_load_si512(data + i * 256 + 224);
			
			// with 128 elements at offset i * 128 + stride (sz/2)
			dv9  = _mm512_load_si512(data + i * 256 + sz/2 +  0);
			dv10 = _mm512_load_si512(data + i * 256 + sz/2 + 32);
			dv11 = _mm512_load_si512(data + i * 256 + sz/2 + 64);
			dv12 = _mm512_load_si512(data + i * 256 + sz/2 + 96);
			dv13 = _mm512_load_si512(data + i * 256 + sz/2 + 128);
			dv14 = _mm512_load_si512(data + i * 256 + sz/2 + 160);
			dv15 = _mm512_load_si512(data + i * 256 + sz/2 + 192);
			dv16 = _mm512_load_si512(data + i * 256 + sz/2 + 224);


			t1 =  _mm512_min_epu16(dv1, dv9 );
			t2 =  _mm512_min_epu16(dv2, dv10);
			t3 =  _mm512_min_epu16(dv3, dv11);
			t4 =  _mm512_min_epu16(dv4, dv12);
			dv9  = _mm512_max_epu16(dv1, dv9 );
			dv10 = _mm512_max_epu16(dv2, dv10);
			dv11 = _mm512_max_epu16(dv3, dv11);
			dv12 = _mm512_max_epu16(dv4, dv12);
			dv1 = t1;
			dv2 = t2;
			dv3 = t3;
			dv4 = t4;
			
			t1 =  _mm512_min_epu16(dv5, dv13);
			t2 =  _mm512_min_epu16(dv6, dv14);
			t3 =  _mm512_min_epu16(dv7, dv15);
			t4 =  _mm512_min_epu16(dv8, dv16);
			dv13 = _mm512_max_epu16(dv5, dv13);
			dv14 = _mm512_max_epu16(dv6, dv14);
			dv15 = _mm512_max_epu16(dv7, dv15);
			dv16 = _mm512_max_epu16(dv8, dv16);
			dv5 = t1;
			dv6 = t2;
			dv7 = t3;
			dv8 = t4;
			
			_mm512_store_si512(data + i * 256 +  0, dv1);
			_mm512_store_si512(data + i * 256 + 32, dv2);
			_mm512_store_si512(data + i * 256 + 64, dv3);
			_mm512_store_si512(data + i * 256 + 96, dv4);
			_mm512_store_si512(data + i * 256 + 128, dv5);
			_mm512_store_si512(data + i * 256 + 160, dv6);
			_mm512_store_si512(data + i * 256 + 192, dv7);
			_mm512_store_si512(data + i * 256 + 224, dv8);
			
			_mm512_store_si512(data + i * 256 + sz/2 +  0, dv9);
			_mm512_store_si512(data + i * 256 + sz/2 + 32, dv10);
			_mm512_store_si512(data + i * 256 + sz/2 + 64, dv11);
			_mm512_store_si512(data + i * 256 + sz/2 + 96, dv12);
			_mm512_store_si512(data + i * 256 + sz/2 + 128, dv13);
			_mm512_store_si512(data + i * 256 + sz/2 + 160, dv14);
			_mm512_store_si512(data + i * 256 + sz/2 + 192, dv15);
			_mm512_store_si512(data + i * 256 + sz/2 + 224, dv16);

		}
	}

	// two parallel half-size merges
	bitonic_merge16(data, sz / 2, dir);
	bitonic_merge16(data + sz / 2, sz / 2, dir);
}

void bitonic_sort16(uint16_t *data, uint32_t sz, int dir)
{
	if (sz == 64)
	{
		// base case: do the hardcoded 128-element sort
		bitonic_sort16_dir_64(data, dir);
		return;
	}
	else if (sz == 128)
	{
		// base case: do the hardcoded 128-element sort
		bitonic_sort16_dir_128(data, dir);
		return;
	}
	else if (sz == 256)
	{
		// base case: do the hardcoded 128-element sort
		bitonic_sort16_dir_256(data, dir);
		return;
	}

	// two half-size bitonic sorts,
	// with opposite directions.
	bitonic_sort16(data, sz / 2, 0);
	bitonic_sort16(data + sz / 2, sz / 2, 1);

	// merge in the specified direction
	bitonic_merge16(data, sz, dir);

	return;
}

void bitonic_merge32(uint32_t *data, uint32_t sz, int dir)
{
	if (sz <= 64)
	{
		// base case: do the hardcoded 128-element sort
		bitonic_merge32_dir_64(data, dir);
		return;
	}
	else if (sz <= 128)
	{
		// base case: do the hardcoded 128-element sort
		bitonic_merge32_dir_128(data, dir);
		return;
	}
	
	// half-size cmp/swap
	int i;

	// we have sz/2 swaps to do at a stride of sz/2.
	// the number of swaps will be divisible by 128 because
	// sz is at least 256 (basecase is 128).
	// can batch up the swaps:
	__m512i t1;
	__m512i t2;
	__m512i t3;
	__m512i t4;
	__m512i t5;
	__m512i t6;
	__m512i t7;
	__m512i t8;
	__m512i dv1;
	__m512i dv2;
	__m512i dv3;
	__m512i dv4;
	__m512i dv5;
	__m512i dv6;
	__m512i dv7;
	__m512i dv8;
	__m512i dv9;
	__m512i dv10;
	__m512i dv11;
	__m512i dv12;
	__m512i dv13;
	__m512i dv14;
	__m512i dv15;
	__m512i dv16;

	__mmask16 m1;
	__mmask16 m2;
	__mmask16 m3;
	__mmask16 m4;
	__mmask16 m5;
	__mmask16 m6;
	__mmask16 m7;
	__mmask16 m8;
			
	if (dir == 1)
	{
		// 256-element merge passes at a stride of sz/2
		for (i = 0; i < sz / 256; i++)
		{
			// swap 128 elements at starting offset i * 128
			dv1 = _mm512_load_si512(data + i * 128 +  0);
			dv2 = _mm512_load_si512(data + i * 128 + 16);
			dv3 = _mm512_load_si512(data + i * 128 + 32);
			dv4 = _mm512_load_si512(data + i * 128 + 48);
			dv5 = _mm512_load_si512(data + i * 128 + 64);
			dv6 = _mm512_load_si512(data + i * 128 + 80);
			dv7 = _mm512_load_si512(data + i * 128 + 96);
			dv8 = _mm512_load_si512(data + i * 128 + 112);
			
			// with 128 elements at offset i * 128 + stride (sz/2)
			dv9 = _mm512_load_si512(data  + i * 128 + sz/2 +  0);
			dv10 = _mm512_load_si512(data + i * 128 + sz/2 + 16);
			dv11 = _mm512_load_si512(data + i * 128 + sz/2 + 32);
			dv12 = _mm512_load_si512(data + i * 128 + sz/2 + 48);
			dv13 = _mm512_load_si512(data + i * 128 + sz/2 + 64);
			dv14 = _mm512_load_si512(data + i * 128 + sz/2 + 80);
			dv15 = _mm512_load_si512(data + i * 128 + sz/2 + 96);
			dv16 = _mm512_load_si512(data + i * 128 + sz/2 + 112);

			t1 =  _mm512_max_epu32(dv1, dv9 );
			t2 =  _mm512_max_epu32(dv2, dv10);
			t3 =  _mm512_max_epu32(dv3, dv11);
			t4 =  _mm512_max_epu32(dv4, dv12);
			t5 =  _mm512_max_epu32(dv5, dv13);
			t6 =  _mm512_max_epu32(dv6, dv14);
			t7 =  _mm512_max_epu32(dv7, dv15);
			t8 =  _mm512_max_epu32(dv8, dv16);
			dv9  = _mm512_min_epu32(dv1, dv9 );
			dv10 = _mm512_min_epu32(dv2, dv10);
			dv11 = _mm512_min_epu32(dv3, dv11);
			dv12 = _mm512_min_epu32(dv4, dv12);
			dv13 = _mm512_min_epu32(dv5, dv13);
			dv14 = _mm512_min_epu32(dv6, dv14);
			dv15 = _mm512_min_epu32(dv7, dv15);
			dv16 = _mm512_min_epu32(dv8, dv16);
			dv1 = t1;
			dv2 = t2;
			dv3 = t3;
			dv4 = t4;
			dv5 = t5;
			dv6 = t6;
			dv7 = t7;
			dv8 = t8;
			
			_mm512_store_si512(data + i * 128 +  0, dv1);
			_mm512_store_si512(data + i * 128 + 16, dv2);
			_mm512_store_si512(data + i * 128 + 32, dv3);
			_mm512_store_si512(data + i * 128 + 48, dv4);
			_mm512_store_si512(data + i * 128 + 64, dv5);
			_mm512_store_si512(data + i * 128 + 80, dv6);
			_mm512_store_si512(data + i * 128 + 96, dv7);
			_mm512_store_si512(data + i * 128 + 112, dv8);
			
			_mm512_store_si512(data + i * 128 + sz/2 +  0,  dv9);
			_mm512_store_si512(data + i * 128 + sz/2 + 16,  dv10);
			_mm512_store_si512(data + i * 128 + sz/2 + 32, dv11);
			_mm512_store_si512(data + i * 128 + sz/2 + 48, dv12);
			_mm512_store_si512(data + i * 128 + sz/2 + 64, dv13);
			_mm512_store_si512(data + i * 128 + sz/2 + 80, dv14);
			_mm512_store_si512(data + i * 128 + sz/2 + 96, dv15);
			_mm512_store_si512(data + i * 128 + sz/2 + 112, dv16);
			
			// swap 128 elements at starting offset i * 128
			//dv1 = _mm512_load_si512(data + i * 256 + 8*16);
			//dv2 = _mm512_load_si512(data + i * 256 + 9*16);
			//dv3 = _mm512_load_si512(data + i * 256 + 10*16);
			//dv4 = _mm512_load_si512(data + i * 256 + 11*16);
			//dv5 = _mm512_load_si512(data + i * 256 + 12*16);
			//dv6 = _mm512_load_si512(data + i * 256 + 13*16);
			//dv7 = _mm512_load_si512(data + i * 256 + 14*16);
			//dv8 = _mm512_load_si512(data + i * 256 + 15*16);
			//
			//// with 128 elements at offset i * 128 + stride (sz/2)
			//dv9 = _mm512_load_si512(data  + i * 256 + sz/2 + 8*16 );
			//dv10 = _mm512_load_si512(data + i * 256 + sz/2 + 9*16 );
			//dv11 = _mm512_load_si512(data + i * 256 + sz/2 + 10*16);
			//dv12 = _mm512_load_si512(data + i * 256 + sz/2 + 11*16);
			//dv13 = _mm512_load_si512(data + i * 256 + sz/2 + 12*16);
			//dv14 = _mm512_load_si512(data + i * 256 + sz/2 + 13*16);
			//dv15 = _mm512_load_si512(data + i * 256 + sz/2 + 14*16);
			//dv16 = _mm512_load_si512(data + i * 256 + sz/2 + 15*16);
			//
			//t1 =  _mm512_max_epu32(dv1, dv9 );
			//t2 =  _mm512_max_epu32(dv2, dv10);
			//t3 =  _mm512_max_epu32(dv3, dv11);
			//t4 =  _mm512_max_epu32(dv4, dv12);
			//t5 =  _mm512_max_epu32(dv5, dv13);
			//t6 =  _mm512_max_epu32(dv6, dv14);
			//t7 =  _mm512_max_epu32(dv7, dv15);
			//t8 =  _mm512_max_epu32(dv8, dv16);
			//dv9  = _mm512_min_epu32(dv1, dv9 );
			//dv10 = _mm512_min_epu32(dv2, dv10);
			//dv11 = _mm512_min_epu32(dv3, dv11);
			//dv12 = _mm512_min_epu32(dv4, dv12);
			//dv13 = _mm512_min_epu32(dv5, dv13);
			//dv14 = _mm512_min_epu32(dv6, dv14);
			//dv15 = _mm512_min_epu32(dv7, dv15);
			//dv16 = _mm512_min_epu32(dv8, dv16);
			//dv1 = t1;
			//dv2 = t2;
			//dv3 = t3;
			//dv4 = t4;
			//dv5 = t5;
			//dv6 = t6;
			//dv7 = t7;
			//dv8 = t8;
			//
			//_mm512_store_si512(data + i * 256 + 8*16 , dv1);
			//_mm512_store_si512(data + i * 256 + 9*16 , dv2);
			//_mm512_store_si512(data + i * 256 + 10*16, dv3);
			//_mm512_store_si512(data + i * 256 + 11*16, dv4);
			//_mm512_store_si512(data + i * 256 + 12*16, dv5);
			//_mm512_store_si512(data + i * 256 + 13*16, dv6);
			//_mm512_store_si512(data + i * 256 + 14*16, dv7);
			//_mm512_store_si512(data + i * 256 + 15*16, dv8);
			//
			//_mm512_store_si512(data + i * 256 + sz/2 + 8*16 , dv1);
			//_mm512_store_si512(data + i * 256 + sz/2 + 9*16 , dv2);
			//_mm512_store_si512(data + i * 256 + sz/2 + 10*16, dv3);
			//_mm512_store_si512(data + i * 256 + sz/2 + 11*16, dv4);
			//_mm512_store_si512(data + i * 256 + sz/2 + 12*16, dv5);
			//_mm512_store_si512(data + i * 256 + sz/2 + 13*16, dv6);
			//_mm512_store_si512(data + i * 256 + sz/2 + 14*16, dv7);
			//_mm512_store_si512(data + i * 256 + sz/2 + 15*16, dv8);

		}
	}
	else
	{
		// 256-element merge passes at a stride of sz/2
		for (i = 0; i < sz / 256; i++)
		{
			// swap 128 elements at starting offset i * 128
			dv1 = _mm512_load_si512(data + i * 128 +  0);
			dv2 = _mm512_load_si512(data + i * 128 + 16);
			dv3 = _mm512_load_si512(data + i * 128 + 32);
			dv4 = _mm512_load_si512(data + i * 128 + 48);
			dv5 = _mm512_load_si512(data + i * 128 + 64);
			dv6 = _mm512_load_si512(data + i * 128 + 80);
			dv7 = _mm512_load_si512(data + i * 128 + 96);
			dv8 = _mm512_load_si512(data + i * 128 + 112);
			
			// with 128 elements at offset i * 128 + stride (sz/2)
			dv9 = _mm512_load_si512(data  + i * 128 + sz/2 +  0);
			dv10 = _mm512_load_si512(data + i * 128 + sz/2 + 16);
			dv11 = _mm512_load_si512(data + i * 128 + sz/2 + 32);
			dv12 = _mm512_load_si512(data + i * 128 + sz/2 + 48);
			dv13 = _mm512_load_si512(data + i * 128 + sz/2 + 64);
			dv14 = _mm512_load_si512(data + i * 128 + sz/2 + 80);
			dv15 = _mm512_load_si512(data + i * 128 + sz/2 + 96);
			dv16 = _mm512_load_si512(data + i * 128 + sz/2 + 112);

			t1 =  _mm512_min_epu32(dv1, dv9 );
			t2 =  _mm512_min_epu32(dv2, dv10);
			t3 =  _mm512_min_epu32(dv3, dv11);
			t4 =  _mm512_min_epu32(dv4, dv12);
			t5 =  _mm512_min_epu32(dv5, dv13);
			t6 =  _mm512_min_epu32(dv6, dv14);
			t7 =  _mm512_min_epu32(dv7, dv15);
			t8 =  _mm512_min_epu32(dv8, dv16);
			dv9  = _mm512_max_epu32(dv1, dv9 );
			dv10 = _mm512_max_epu32(dv2, dv10);
			dv11 = _mm512_max_epu32(dv3, dv11);
			dv12 = _mm512_max_epu32(dv4, dv12);
			dv13 = _mm512_max_epu32(dv5, dv13);
			dv14 = _mm512_max_epu32(dv6, dv14);
			dv15 = _mm512_max_epu32(dv7, dv15);
			dv16 = _mm512_max_epu32(dv8, dv16);
			dv1 = t1;
			dv2 = t2;
			dv3 = t3;
			dv4 = t4;
			dv5 = t5;
			dv6 = t6;
			dv7 = t7;
			dv8 = t8;
			
			_mm512_store_si512(data + i * 128 +  0, dv1);
			_mm512_store_si512(data + i * 128 + 16, dv2);
			_mm512_store_si512(data + i * 128 + 32, dv3);
			_mm512_store_si512(data + i * 128 + 48, dv4);
			_mm512_store_si512(data + i * 128 + 64, dv5);
			_mm512_store_si512(data + i * 128 + 80, dv6);
			_mm512_store_si512(data + i * 128 + 96, dv7);
			_mm512_store_si512(data + i * 128 + 112, dv8);
			
			_mm512_store_si512(data + i * 128 + sz/2 +  0,  dv9);
			_mm512_store_si512(data + i * 128 + sz/2 + 16,  dv10);
			_mm512_store_si512(data + i * 128 + sz/2 + 32, dv11);
			_mm512_store_si512(data + i * 128 + sz/2 + 48, dv12);
			_mm512_store_si512(data + i * 128 + sz/2 + 64, dv13);
			_mm512_store_si512(data + i * 128 + sz/2 + 80, dv14);
			_mm512_store_si512(data + i * 128 + sz/2 + 96, dv15);
			_mm512_store_si512(data + i * 128 + sz/2 + 112, dv16);
			
			// swap 128 elements at starting offset i * 128
			//dv1 = _mm512_load_si512(data + i * 256 + 8*16);
			//dv2 = _mm512_load_si512(data + i * 256 + 9*16);
			//dv3 = _mm512_load_si512(data + i * 256 + 10*16);
			//dv4 = _mm512_load_si512(data + i * 256 + 11*16);
			//dv5 = _mm512_load_si512(data + i * 256 + 12*16);
			//dv6 = _mm512_load_si512(data + i * 256 + 13*16);
			//dv7 = _mm512_load_si512(data + i * 256 + 14*16);
			//dv8 = _mm512_load_si512(data + i * 256 + 15*16);
			//
			//// with 128 elements at offset i * 128 + stride (sz/2)
			//dv9 = _mm512_load_si512(data  + i * 256 + sz/2 + 8*16 );
			//dv10 = _mm512_load_si512(data + i * 256 + sz/2 + 9*16 );
			//dv11 = _mm512_load_si512(data + i * 256 + sz/2 + 10*16);
			//dv12 = _mm512_load_si512(data + i * 256 + sz/2 + 11*16);
			//dv13 = _mm512_load_si512(data + i * 256 + sz/2 + 12*16);
			//dv14 = _mm512_load_si512(data + i * 256 + sz/2 + 13*16);
			//dv15 = _mm512_load_si512(data + i * 256 + sz/2 + 14*16);
			//dv16 = _mm512_load_si512(data + i * 256 + sz/2 + 15*16);
			//
			//t1 =  _mm512_min_epu32(dv1, dv9 );
			//t2 =  _mm512_min_epu32(dv2, dv10);
			//t3 =  _mm512_min_epu32(dv3, dv11);
			//t4 =  _mm512_min_epu32(dv4, dv12);
			//t5 =  _mm512_min_epu32(dv5, dv13);
			//t6 =  _mm512_min_epu32(dv6, dv14);
			//t7 =  _mm512_min_epu32(dv7, dv15);
			//t8 =  _mm512_min_epu32(dv8, dv16);
			//dv9  = _mm512_max_epu32(dv1, dv9 );
			//dv10 = _mm512_max_epu32(dv2, dv10);
			//dv11 = _mm512_max_epu32(dv3, dv11);
			//dv12 = _mm512_max_epu32(dv4, dv12);
			//dv13 = _mm512_max_epu32(dv5, dv13);
			//dv14 = _mm512_max_epu32(dv6, dv14);
			//dv15 = _mm512_max_epu32(dv7, dv15);
			//dv16 = _mm512_max_epu32(dv8, dv16);
			//dv1 = t1;
			//dv2 = t2;
			//dv3 = t3;
			//dv4 = t4;
			//dv5 = t5;
			//dv6 = t6;
			//dv7 = t7;
			//dv8 = t8;
			//
			//_mm512_store_si512(data + i * 256 + 8*16 , dv1);
			//_mm512_store_si512(data + i * 256 + 9*16 , dv2);
			//_mm512_store_si512(data + i * 256 + 10*16, dv3);
			//_mm512_store_si512(data + i * 256 + 11*16, dv4);
			//_mm512_store_si512(data + i * 256 + 12*16, dv5);
			//_mm512_store_si512(data + i * 256 + 13*16, dv6);
			//_mm512_store_si512(data + i * 256 + 14*16, dv7);
			//_mm512_store_si512(data + i * 256 + 15*16, dv8);
			//
			//_mm512_store_si512(data + i * 256 + sz/2 + 8*16 , dv1);
			//_mm512_store_si512(data + i * 256 + sz/2 + 9*16 , dv2);
			//_mm512_store_si512(data + i * 256 + sz/2 + 10*16, dv3);
			//_mm512_store_si512(data + i * 256 + sz/2 + 11*16, dv4);
			//_mm512_store_si512(data + i * 256 + sz/2 + 12*16, dv5);
			//_mm512_store_si512(data + i * 256 + sz/2 + 13*16, dv6);
			//_mm512_store_si512(data + i * 256 + sz/2 + 14*16, dv7);
			//_mm512_store_si512(data + i * 256 + sz/2 + 15*16, dv8);

		}
	}

	// two parallel half-size merges
	bitonic_merge32(data, sz / 2, dir);
	bitonic_merge32(data + sz / 2, sz / 2, dir);
}

void bitonic_sort32(uint32_t *data, uint32_t sz, int dir)
{
	if (sz <= 64)
	{
		// base case: do the hardcoded 64-element sort
		bitonic_sort32_dir_64(data, dir);
		return;
	}
	else if (sz <= 128)
	{
		// base case: do the hardcoded 128-element sort
		bitonic_sort32_dir_128(data, dir);
		return;
	}
	else if (sz <= 256)
	{
		// base case: do the hardcoded 256-element sort
		bitonic_sort32_dir_256(data, dir);
		return;
	}

	// two half-size bitonic sorts,
	// with opposite directions.
	bitonic_sort32(data, sz / 2, 0);
	bitonic_sort32(data + sz / 2, sz / 2, 1);

	// merge in the specified direction
	bitonic_merge32(data, sz, dir);

	return;
}

void sort32(uint32_t *data, uint32_t sz, int dir)
{
	// top level sort dealing with two things:
	// 1) the bitonic sort function requires the data array
	// to be aligned and
	// 2) the bitonic sort works on a power-of-2 sized array
	// here we make sure we feed the bitonic sort function
	// an array that satifies those requirements.
	int is_aligned = (((uint64_t)data & 0x3full) == 0);
	if (is_aligned && ((sz & (sz - 1)) == 0))
	{
		// meets both requirements as-is
		bitonic_sort32(data, sz, dir);
		return;
	}
	
	// otherwise we need to copy to an aligned buffer and/or 
	// change the buffer size
	uint32_t new_sz = sz;
	if ((sz & (sz - 1)) > 0)
	{
		new_sz = my_clz32(sz);
		if (new_sz == 0)
		{
			printf("buffer too big, sz must be <= 2^31 in sort()\n");
			exit(0);
		}
		new_sz = 1 << (32 - new_sz + 1);
	}
		
	uint32_t *adata;
	
	if (!is_aligned)
	{
		adata = (uint32_t*)aligned_malloc(new_sz * sizeof(uint32_t), 64);
		memcpy(adata, data, sz * sizeof(uint32_t));
	}
	else
	{
		adata = data;
	}
	
	if ((new_sz - sz) > 0)
	{
		if (dir == 0)
			memset(adata + sz, 0xff, (new_sz - sz));
		else
			memset(adata + sz, 0, (new_sz - sz));
	}
	
	bitonic_sort32(adata, new_sz, dir);
	
	if (!is_aligned)
	{
		memcpy(data, adata, sz * sizeof(uint32_t));
		aligned_free(adata);
	}
	
	return;
}

// for 20M element problems
#ifndef SH
#define SH 8
#endif
#define NB (1<<SH)
#define NB1 ((NB)-1)
#define portion_sz 64

void bucket_sort16(uint16_t *buckets, uint32_t bucket_size, uint32_t *bucket_counts, 
	uint64_t *data, uint32_t key_bits, uint32_t sz)
{
	// sort data into buckets by most significant bits.
	// assume the buckets have been allocated of sufficient size for the random data.
	// the buckets are subdivided into cache-friendly portions of size portion_sz * NB.
	// once a particular portion is filled for a bucket, increment that portion's offset by
	// portion_sz * NB to begin filling the next portion.
	int i;
	uint16_t *current_bucket[NB];
	uint32_t portion_cnt[NB];

	for (i = 0; i < NB; i++) 
	{
		current_bucket[i] = buckets + i * portion_sz;
		bucket_counts[i] = 0;
		portion_cnt[i] = 0;
	}
	
	for (i = 0; i < sz; i++) {
		uint32 key = (uint32)((data[i] >> 16) & NB1);
		uint16_t *bptr = current_bucket[key];
		bptr[portion_cnt[key]] = data[i] & 0xffff;
		portion_cnt[key]++;

		if (portion_cnt[key] == portion_sz)
		{
			current_bucket[key] += portion_sz * NB;
			bucket_counts[key] += portion_sz;
			portion_cnt[key] = 0;
		}
	}
	
	for (i = 0; i < NB; i++)
	{
		bucket_counts[i] += portion_cnt[i];
	}
		
	return;
}

void bucket_sort32(uint32_t *buckets, uint32_t bucket_size, uint32_t *bucket_counts, 
	uint64_t *data, uint32_t key_bits, uint32_t sz)
{
	// sort data into buckets by most significant bits.
	// assume the buckets have been allocated of sufficient size for the random data.
	// the buckets are subdivided into cache-friendly portions of size portion_sz * NB.
	// once a particular portion is filled for a bucket, increment that portion's offset by
	// portion_sz * NB to begin filling the next portion.
	int i;
	uint32_t *current_bucket[NB];
	uint32_t portion_cnt[NB];
	uint32 shift, mask;
	
	if ((32 + SH) > key_bits)
	{
		shift = key_bits - SH;
		mask = (1 << shift)  - 1;
	}
	else
	{
		shift = 32;
		mask = 0xffffffff;
	}
	
	for (i = 0; i < NB; i++) 
	{
		current_bucket[i] = buckets + i * portion_sz;
		bucket_counts[i] = 0;
		portion_cnt[i] = 0;
	}
	
	for (i = 0; i < sz; i++) {
		uint32 key = (uint32)((data[i] >> shift) & NB1);
		uint32_t *bptr = current_bucket[key];
		bptr[portion_cnt[key]] = data[i] & mask;
		portion_cnt[key]++;

		if (portion_cnt[key] == portion_sz)
		{
			current_bucket[key] += portion_sz * NB;
			bucket_counts[key] += portion_sz;
			portion_cnt[key] = 0;
		}
	}
	
	for (i = 0; i < NB; i++)
	{
		bucket_counts[i] += portion_cnt[i];
	}
		
	return;
}

#endif

	