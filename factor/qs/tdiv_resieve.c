/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

Some parts of the code (and also this header), included in this 
distribution have been reused from other sources. In particular I 
have benefitted greatly from the work of Jason Papadopoulos's msieve @ 
www.boo.net/~jasonp, Scott Contini's mpqs implementation, and Tom St. 
Denis Tom's Fast Math library.  Many thanks to their kind donation of 
code to the public domain.
       				   --bbuhrow@gmail.com 11/24/09
----------------------------------------------------------------------*/

#include "yafu.h"
#include "qs.h"
#include "factor.h"
#include "util.h"
#include "common.h"

//#define SIQSDEBUG 1

/*
We are given an array of bytes that has been sieved.  The basic trial 
division strategy is as follows:

1) Scan through the array and 'mark' locations that meet criteria 
indicating they may factor completely over the factor base.  

2) 'Filter' the marked locations by trial dividing by small primes
that we did not sieve.  These primes are all less than 256.  If after
removing small primes the location does not meet another set of criteria,
remove it from the 'marked' list (do not subject it to further trial
division).

3) Divide out primes from the factor base between 256 and 2^13 or 2^14, 
depending on the version (2^13 for 32k version, 2^14 for 64k).  

4) Resieve primes between 2^{13|14} and 2^16, max.  

5) Primes larger than 2^16 will have been bucket sieved.  Remove these
by scanning the buckets for sieve hits equal to the current block location.

6) If applicable/appropriate, factor a remaining composite with squfof

this file contains code implementing 4)


*/


#if defined(GCC_ASM32X) || defined(GCC_ASM64X) || defined(__MINGW32__)
	
	#if defined(HAS_SSE2)
		#define SSE2_RESIEVING 1

		#define STEP_COMPARE_COMBINE \
			"psubw %%xmm1, %%xmm2 \n\t"		/* subtract primes from root1s */ \
			"psubw %%xmm1, %%xmm3 \n\t"		/* subtract primes from root2s */ \
			"pcmpeqw %%xmm2, %%xmm5 \n\t"	/* root1s ?= 0 */ \
			"pcmpeqw %%xmm3, %%xmm6 \n\t"	/* root2s ?= 0 */ \
			"por %%xmm5, %%xmm7 \n\t"		/* combine results */ \
			"por %%xmm6, %%xmm8 \n\t"		/* combine results */

		#define INIT_RESIEVE \
			"movdqa (%4), %%xmm4 \n\t"		/* bring in corrections to roots */				\
			"pxor %%xmm8, %%xmm8 \n\t"		/* zero xmm8 */ \
			"movdqa (%2), %%xmm2 \n\t"		/* bring in 8 root1s */ \
			"paddw %%xmm4, %%xmm2 \n\t"		/* correct root1s */ \
			"movdqa (%3), %%xmm3 \n\t"		/* bring in 8 root2s */ \
			"paddw %%xmm4, %%xmm3 \n\t"		/* correct root2s */ \
			"movdqa (%1), %%xmm1 \n\t"		/* bring in 8 primes */ \
			"pxor %%xmm7, %%xmm7 \n\t"		/* zero xmm7 */ \
			"pxor %%xmm5, %%xmm5 \n\t"		/* zero xmm5 */ \
			"pxor %%xmm6, %%xmm6 \n\t"		/* zero xmm6 */

		#ifdef YAFU_64K
			#define RESIEVE_8X_14BIT_MAX \
				asm ( \
					INIT_RESIEVE \
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					"por	%%xmm8, %%xmm7 \n\t" \
					"pmovmskb %%xmm7, %0 \n\t"		/* if one of these primes divides this location, this will be !0*/ \
					: "=r"(result) \
					: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections) \
					: "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "cc", "memory" \
					);

			#define RESIEVE_8X_15BIT_MAX \
				asm ( \
					INIT_RESIEVE \
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					"por	%%xmm8, %%xmm7 \n\t" \
					"pmovmskb %%xmm7, %0 \n\t"		/* if one of these primes divides this location, this will be !0*/ \
					: "=r"(result) \
					: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections) \
					: "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "cc", "memory" \
					);

			#define RESIEVE_8X_16BIT_MAX \
				asm ( \
					INIT_RESIEVE \
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					"por	%%xmm8, %%xmm7 \n\t" \
					"pmovmskb %%xmm7, %0 \n\t"		/* if one of these primes divides this location, this will be !0*/ \
					: "=r"(result) \
					: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections) \
					: "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "cc", "memory" \
					);
		#else
			#define RESIEVE_8X_14BIT_MAX \
				asm ( \
					INIT_RESIEVE \
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					"por	%%xmm8, %%xmm7 \n\t" \
					"pmovmskb %%xmm7, %0 \n\t"		/* if one of these primes divides this location, this will be !0*/ \
					: "=r"(result) \
					: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections) \
					: "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "cc", "memory" \
					);

			#define RESIEVE_8X_15BIT_MAX \
				asm ( \
					INIT_RESIEVE \
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					"por	%%xmm8, %%xmm7 \n\t" \
					"pmovmskb %%xmm7, %0 \n\t"		/* if one of these primes divides this location, this will be !0*/ \
					: "=r"(result) \
					: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections) \
					: "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "cc", "memory" \
					);

			#define RESIEVE_8X_16BIT_MAX \
				asm ( \
					INIT_RESIEVE \
					STEP_COMPARE_COMBINE	\
					"por	%%xmm8, %%xmm7 \n\t" \
					"pmovmskb %%xmm7, %0 \n\t"		/* if one of these primes divides this location, this will be !0*/ \
					: "=r"(result) \
					: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections) \
					: "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "cc", "memory" \
					);
		#endif

	#else
		#error SSE2 is required
	#endif


#elif defined(_MSC_VER)

	#define STEP_COMPARE_COMBINE \
		root1s = _mm_sub_epi16(root1s, primes); \
		root2s = _mm_sub_epi16(root2s, primes); \
		tmp1 = _mm_cmpeq_epi16(tmp1, root1s); \
		tmp2 = _mm_cmpeq_epi16(tmp2, root2s); \
		combine = _mm_xor_si128(combine, tmp1); \
		combine = _mm_xor_si128(combine, tmp2);

	#define INIT_RESIEVE \
		c = _mm_load_si128((__m128i *)corrections); \
		root1s = _mm_load_si128((__m128i *)(fbc->root1 + i)); \
		root1s = _mm_add_epi16(root1s, c); \
		root2s = _mm_load_si128((__m128i *)(fbc->root2 + i)); \
		root2s = _mm_add_epi16(root2s, c); \
		primes = _mm_load_si128((__m128i *)(fbc->prime + i)); \
		combine = _mm_xor_si128(combine, combine); \
		tmp1 = _mm_xor_si128(tmp1, tmp1); \
		tmp2 = _mm_xor_si128(tmp2, tmp2);

	#ifdef YAFU_64K

		#define RESIEVE_8X_14BIT_MAX \
			do { \
				__m128i tmp1;	\
				__m128i tmp2;	\
				__m128i root1s;	\
				__m128i root2s;	\
				__m128i primes;	\
				__m128i c;	\
				__m128i combine;	\
				INIT_RESIEVE \
				STEP_COMPARE_COMBINE	\
				STEP_COMPARE_COMBINE	\
				STEP_COMPARE_COMBINE	\
				STEP_COMPARE_COMBINE	\
				STEP_COMPARE_COMBINE	\
				STEP_COMPARE_COMBINE	\
				STEP_COMPARE_COMBINE	\
				STEP_COMPARE_COMBINE	\
				result = _mm_movemask_epi8(combine); \
			} while (0);

		#define RESIEVE_8X_15BIT_MAX \
			do { \
				__m128i tmp1;	\
				__m128i tmp2;	\
				__m128i root1s;	\
				__m128i root2s;	\
				__m128i primes;	\
				__m128i c;	\
				__m128i combine;	\
				INIT_RESIEVE \
				STEP_COMPARE_COMBINE	\
				STEP_COMPARE_COMBINE	\
				STEP_COMPARE_COMBINE	\
				STEP_COMPARE_COMBINE	\
				result = _mm_movemask_epi8(combine); \
			} while (0);

		#define RESIEVE_8X_16BIT_MAX \
			do { \
				__m128i tmp1;	\
				__m128i tmp2;	\
				__m128i root1s;	\
				__m128i root2s;	\
				__m128i primes;	\
				__m128i c;	\
				__m128i combine;	\
				INIT_RESIEVE \
				STEP_COMPARE_COMBINE	\
				STEP_COMPARE_COMBINE	\
				result = _mm_movemask_epi8(combine); \
			} while (0);

	#else

		#define RESIEVE_8X_14BIT_MAX \
			do { \
				__m128i tmp1;	\
				__m128i tmp2;	\
				__m128i root1s;	\
				__m128i root2s;	\
				__m128i primes;	\
				__m128i c;	\
				__m128i combine;	\
				INIT_RESIEVE \
				STEP_COMPARE_COMBINE	\
				STEP_COMPARE_COMBINE	\
				STEP_COMPARE_COMBINE	\
				STEP_COMPARE_COMBINE	\
				result = _mm_movemask_epi8(combine); \
			} while (0);

		#define RESIEVE_8X_15BIT_MAX \
			do { \
				__m128i tmp1;	\
				__m128i tmp2;	\
				__m128i root1s;	\
				__m128i root2s;	\
				__m128i primes;	\
				__m128i c;	\
				__m128i combine;	\
				INIT_RESIEVE \
				STEP_COMPARE_COMBINE	\
				STEP_COMPARE_COMBINE	\
				result = _mm_movemask_epi8(combine); \
			} while (0);

		#define RESIEVE_8X_16BIT_MAX \
			do { \
				__m128i tmp1;	\
				__m128i tmp2;	\
				__m128i root1s;	\
				__m128i root2s;	\
				__m128i primes;	\
				__m128i c;	\
				__m128i combine;	\
				INIT_RESIEVE \
				STEP_COMPARE_COMBINE	\
				result = _mm_movemask_epi8(combine); \
			} while (0);

	#endif

#else	/* compiler not recognized*/

	#define COMPARE_RESIEVE_VALS(x)	\
		if (r1 == 0) result |= x; \
		if (r2 == 0) result |= x;

	#define STEP_RESIEVE \
		r1 -= p;	\
		r2 -= p;

	#define RESIEVE_1X_14BIT_MAX(x)	\
		STEP_RESIEVE;						\
		COMPARE_RESIEVE_VALS(x);			\
		STEP_RESIEVE;						\
		COMPARE_RESIEVE_VALS(x);			\
		STEP_RESIEVE;						\
		COMPARE_RESIEVE_VALS(x);			\
		STEP_RESIEVE;						\
		COMPARE_RESIEVE_VALS(x);

	#define RESIEVE_1X_15BIT_MAX(x)	\
		STEP_RESIEVE;						\
		COMPARE_RESIEVE_VALS(x);			\
		STEP_RESIEVE;						\
		COMPARE_RESIEVE_VALS(x);

	#define RESIEVE_1X_16BIT_MAX(x)	\
		STEP_RESIEVE;						\
		COMPARE_RESIEVE_VALS(x);

	#ifdef YAFU_64K
		#define RESIEVE_8X_14BIT_MAX \
			do { \
				int p = (int)fbc->prime[i];										\
				int r1 = (int)fbc->root1[i] + BLOCKSIZE - block_loc;			\
				int r2 = (int)fbc->root2[i] + BLOCKSIZE - block_loc;			\
				RESIEVE_1X_14BIT_MAX(0x2);										\
				RESIEVE_1X_14BIT_MAX(0x2);										\
				\
				p = (int)fbc->prime[i+1];										\
				r1 = (int)fbc->root1[i+1] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+1] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_14BIT_MAX(0x8);										\
				RESIEVE_1X_14BIT_MAX(0x8);										\
				\
				p = (int)fbc->prime[i+2];										\
				r1 = (int)fbc->root1[i+2] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+2] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_14BIT_MAX(0x20);									\
				RESIEVE_1X_14BIT_MAX(0x20);									\
				\
				p = (int)fbc->prime[i+3];										\
				r1 = (int)fbc->root1[i+3] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+3] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_14BIT_MAX(0x80);									\
				RESIEVE_1X_14BIT_MAX(0x80);									\
				\
				p = (int)fbc->prime[i+4];										\
				r1 = (int)fbc->root1[i+4] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+4] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_14BIT_MAX(0x200);									\
				RESIEVE_1X_14BIT_MAX(0x200);									\
				\
				p = (int)fbc->prime[i+5];										\
				r1 = (int)fbc->root1[i+5] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+5] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_14BIT_MAX(0x800);									\
				RESIEVE_1X_14BIT_MAX(0x800);									\
				\
				p = (int)fbc->prime[i+6];										\
				r1 = (int)fbc->root1[i+6] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+6] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_14BIT_MAX(0x2000);									\
				RESIEVE_1X_14BIT_MAX(0x2000);									\
				\
				p = (int)fbc->prime[i+7];										\
				r1 = (int)fbc->root1[i+7] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+7] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_14BIT_MAX(0x8000);									\
				RESIEVE_1X_14BIT_MAX(0x8000);									\
			} while (0); 

		#define RESIEVE_8X_15BIT_MAX \
			do { \
				int p = (int)fbc->prime[i];										\
				int r1 = (int)fbc->root1[i] + BLOCKSIZE - block_loc;			\
				int r2 = (int)fbc->root2[i] + BLOCKSIZE - block_loc;			\
				RESIEVE_1X_15BIT_MAX(0x2);										\
				RESIEVE_1X_15BIT_MAX(0x2);										\
				\
				p = (int)fbc->prime[i+1];										\
				r1 = (int)fbc->root1[i+1] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+1] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_15BIT_MAX(0x8);										\
				RESIEVE_1X_15BIT_MAX(0x8);										\
				\
				p = (int)fbc->prime[i+2];										\
				r1 = (int)fbc->root1[i+2] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+2] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_15BIT_MAX(0x20);									\
				RESIEVE_1X_15BIT_MAX(0x20);									\
				\
				p = (int)fbc->prime[i+3];										\
				r1 = (int)fbc->root1[i+3] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+3] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_15BIT_MAX(0x80);									\
				RESIEVE_1X_15BIT_MAX(0x80);									\
				\
				p = (int)fbc->prime[i+4];										\
				r1 = (int)fbc->root1[i+4] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+4] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_15BIT_MAX(0x200);									\
				RESIEVE_1X_15BIT_MAX(0x200);									\
				\
				p = (int)fbc->prime[i+5];										\
				r1 = (int)fbc->root1[i+5] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+5] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_15BIT_MAX(0x800);									\
				RESIEVE_1X_15BIT_MAX(0x800);									\
				\
				p = (int)fbc->prime[i+6];										\
				r1 = (int)fbc->root1[i+6] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+6] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_15BIT_MAX(0x2000);									\
				RESIEVE_1X_15BIT_MAX(0x2000);									\
				\
				p = (int)fbc->prime[i+7];										\
				r1 = (int)fbc->root1[i+7] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+7] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_15BIT_MAX(0x8000);									\
				RESIEVE_1X_15BIT_MAX(0x8000);									\
			} while (0); 

		#define RESIEVE_8X_16BIT_MAX \
			do { \
				int p = (int)fbc->prime[i];										\
				int r1 = (int)fbc->root1[i] + BLOCKSIZE - block_loc;			\
				int r2 = (int)fbc->root2[i] + BLOCKSIZE - block_loc;			\
				RESIEVE_1X_16BIT_MAX(0x2);										\
				RESIEVE_1X_16BIT_MAX(0x2);										\
				\
				p = (int)fbc->prime[i+1];										\
				r1 = (int)fbc->root1[i+1] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+1] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_16BIT_MAX(0x8);										\
				RESIEVE_1X_16BIT_MAX(0x8);										\
				\
				p = (int)fbc->prime[i+2];										\
				r1 = (int)fbc->root1[i+2] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+2] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_16BIT_MAX(0x20);									\
				RESIEVE_1X_16BIT_MAX(0x20);									\
				\
				p = (int)fbc->prime[i+3];										\
				r1 = (int)fbc->root1[i+3] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+3] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_16BIT_MAX(0x80);									\
				RESIEVE_1X_16BIT_MAX(0x80);									\
				\
				p = (int)fbc->prime[i+4];										\
				r1 = (int)fbc->root1[i+4] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+4] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_16BIT_MAX(0x200);									\
				RESIEVE_1X_16BIT_MAX(0x200);									\
				\
				p = (int)fbc->prime[i+5];										\
				r1 = (int)fbc->root1[i+5] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+5] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_16BIT_MAX(0x800);									\
				RESIEVE_1X_16BIT_MAX(0x800);									\
				\
				p = (int)fbc->prime[i+6];										\
				r1 = (int)fbc->root1[i+6] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+6] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_16BIT_MAX(0x2000);									\
				RESIEVE_1X_16BIT_MAX(0x2000);									\
				\
				p = (int)fbc->prime[i+7];										\
				r1 = (int)fbc->root1[i+7] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+7] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_16BIT_MAX(0x8000);									\
				RESIEVE_1X_16BIT_MAX(0x8000);									\
			} while (0); 
	#else
		#define RESIEVE_8X_14BIT_MAX \
			do { \
				int p = (int)fbc->prime[i];										\
				int r1 = (int)fbc->root1[i] + BLOCKSIZE - block_loc;			\
				int r2 = (int)fbc->root2[i] + BLOCKSIZE - block_loc;			\
				RESIEVE_1X_14BIT_MAX(0x2);										\
				\
				p = (int)fbc->prime[i+1];										\
				r1 = (int)fbc->root1[i+1] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+1] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_14BIT_MAX(0x8);										\
				\
				p = (int)fbc->prime[i+2];										\
				r1 = (int)fbc->root1[i+2] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+2] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_14BIT_MAX(0x20);									\
				\
				p = (int)fbc->prime[i+3];										\
				r1 = (int)fbc->root1[i+3] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+3] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_14BIT_MAX(0x80);									\
				\
				p = (int)fbc->prime[i+4];										\
				r1 = (int)fbc->root1[i+4] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+4] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_14BIT_MAX(0x200);									\
				\
				p = (int)fbc->prime[i+5];										\
				r1 = (int)fbc->root1[i+5] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+5] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_14BIT_MAX(0x800);									\
				\
				p = (int)fbc->prime[i+6];										\
				r1 = (int)fbc->root1[i+6] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+6] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_14BIT_MAX(0x2000);									\
				\
				p = (int)fbc->prime[i+7];										\
				r1 = (int)fbc->root1[i+7] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+7] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_14BIT_MAX(0x8000);									\
			} while (0); 

		#define RESIEVE_8X_15BIT_MAX \
			do { \
				int p = (int)fbc->prime[i];										\
				int r1 = (int)fbc->root1[i] + BLOCKSIZE - block_loc;			\
				int r2 = (int)fbc->root2[i] + BLOCKSIZE - block_loc;			\
				RESIEVE_1X_15BIT_MAX(0x2);										\
				\
				p = (int)fbc->prime[i+1];										\
				r1 = (int)fbc->root1[i+1] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+1] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_15BIT_MAX(0x8);										\
				\
				p = (int)fbc->prime[i+2];										\
				r1 = (int)fbc->root1[i+2] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+2] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_15BIT_MAX(0x20);									\
				\
				p = (int)fbc->prime[i+3];										\
				r1 = (int)fbc->root1[i+3] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+3] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_15BIT_MAX(0x80);									\
				\
				p = (int)fbc->prime[i+4];										\
				r1 = (int)fbc->root1[i+4] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+4] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_15BIT_MAX(0x200);									\
				\
				p = (int)fbc->prime[i+5];										\
				r1 = (int)fbc->root1[i+5] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+5] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_15BIT_MAX(0x800);									\
				\
				p = (int)fbc->prime[i+6];										\
				r1 = (int)fbc->root1[i+6] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+6] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_15BIT_MAX(0x2000);									\
				\
				p = (int)fbc->prime[i+7];										\
				r1 = (int)fbc->root1[i+7] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+7] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_15BIT_MAX(0x8000);									\
			} while (0); 

		#define RESIEVE_8X_16BIT_MAX \
			do { \
				int p = (int)fbc->prime[i];										\
				int r1 = (int)fbc->root1[i] + BLOCKSIZE - block_loc;			\
				int r2 = (int)fbc->root2[i] + BLOCKSIZE - block_loc;			\
				RESIEVE_1X_16BIT_MAX(0x2);										\
				\
				p = (int)fbc->prime[i+1];										\
				r1 = (int)fbc->root1[i+1] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+1] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_16BIT_MAX(0x8);										\
				\
				p = (int)fbc->prime[i+2];										\
				r1 = (int)fbc->root1[i+2] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+2] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_16BIT_MAX(0x20);									\
				\
				p = (int)fbc->prime[i+3];										\
				r1 = (int)fbc->root1[i+3] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+3] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_16BIT_MAX(0x80);									\
				\
				p = (int)fbc->prime[i+4];										\
				r1 = (int)fbc->root1[i+4] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+4] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_16BIT_MAX(0x200);									\
				\
				p = (int)fbc->prime[i+5];										\
				r1 = (int)fbc->root1[i+5] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+5] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_16BIT_MAX(0x800);									\
				\
				p = (int)fbc->prime[i+6];										\
				r1 = (int)fbc->root1[i+6] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+6] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_16BIT_MAX(0x2000);									\
				\
				p = (int)fbc->prime[i+7];										\
				r1 = (int)fbc->root1[i+7] + BLOCKSIZE - block_loc;				\
				r2 = (int)fbc->root2[i+7] + BLOCKSIZE - block_loc;				\
				RESIEVE_1X_16BIT_MAX(0x8000);									\
			} while (0); 
		
	#endif
#endif

#ifdef USE_YAFU_TDIV
#define DIVIDE_ONE_PRIME \
	do	\
	{	\
		fb_offsets[++smooth_num] = i;	\
		zShortDiv32(tmp32, prime, tmp32);	\
	} while (zShortMod32(tmp32, prime) == 0);
#else
#define DIVIDE_ONE_PRIME \
	do \
	{						\
		fb_offsets[++smooth_num] = i;	\
		mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], prime); \
	} while (mpz_tdiv_ui(dconf->Qvals[report_num], prime) == 0); 
#endif

#ifdef USE_YAFU_TDIV
#define DIVIDE_RESIEVED_PRIME(j) \
	while (zShortMod32(tmp32, fbc->prime[i+j]) == 0)	\
	{	\
		fb_offsets[++smooth_num] = i+j;	\
		zShortDiv32(tmp32, fbc->prime[i+j], tmp32);	\
	}
#else
#define DIVIDE_RESIEVED_PRIME(j) \
	while (mpz_tdiv_ui(dconf->Qvals[report_num], fbc->prime[i+j]) == 0) \
	{						\
		fb_offsets[++smooth_num] = i+j;	\
		mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], fbc->prime[i+j]);		\
	}
#endif

void resieve_medprimes(uint8 parity, uint32 poly_id, uint32 bnum, 
						 static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//we have flagged this sieve offset as likely to produce a relation
	//nothing left to do now but check and see.
	int i;
	uint32 bound, report_num;
	int smooth_num;
	uint32 *fb_offsets;
	sieve_fb *fb;
	uint64 q64;
	sieve_fb_compressed *fbc;
	fb_element_siqs *fullfb_ptr, *fullfb = sconf->factor_base->list;
	uint32 block_loc;
	uint16 *corrections = dconf->corrections;

	fullfb_ptr = fullfb;
	if (parity)
	{
		fb = dconf->fb_sieve_n;
		fbc = dconf->comp_sieve_n;
	}
	else
	{
		fb = dconf->fb_sieve_p;
		fbc = dconf->comp_sieve_p;
	}		

#ifdef QS_TIMING
	gettimeofday(&qs_timing_start, NULL);
#endif

	for (report_num = 0; report_num < dconf->num_reports; report_num++)
	{
#ifdef USE_YAFU_TDIV
		z32 *tmp32 = &dconf->Qvals32[report_num];
#endif

		if (!dconf->valid_Qs[report_num])
			continue;

		// pull the details of this report to get started.
		fb_offsets = &dconf->fb_offsets[report_num][0];
		smooth_num = dconf->smooth_num[report_num];
		block_loc = dconf->reports[report_num];
		
		// where tdiv_medprimes left off
#if defined(YAFU_64K)
		i = sconf->factor_base->fb_14bit_B;
#else
		i = sconf->factor_base->fb_13bit_B;
#endif

#ifdef USE_RESIEVING

#if defined(SSE2_RESIEVING)
		corrections[0] = BLOCKSIZE - block_loc;
		corrections[1] = BLOCKSIZE - block_loc;
		corrections[2] = BLOCKSIZE - block_loc;
		corrections[3] = BLOCKSIZE - block_loc;		
		corrections[4] = BLOCKSIZE - block_loc;
		corrections[5] = BLOCKSIZE - block_loc;
		corrections[6] = BLOCKSIZE - block_loc;
		corrections[7] = BLOCKSIZE - block_loc;		
#endif
		
#ifndef YAFU_64K
		// since the blocksize is bigger for YAFU_64K, skip resieving of primes
		// up to 14 bits in size and instead use standard normal trial division.

		bound = sconf->factor_base->fb_14bit_B;
		while ((uint32)i < bound)
		{
			//minimum prime > blocksize / 2
			//maximum correction = blocksize
			//maximum starting value > blocksize * 3/2
			//max steps = 2
			//this misses all reports at block_loc = 0.  would need to check
			//for equality to blocksize in that case
			//printf("prime = %u, roots = %u,%u.  max steps = %u\n",
			//	fbc->prime[i],fbc->root1[i],fbc->root2[i],
			//	(fbc->prime[i] - 1 + BLOCKSIZE) / fbc->prime[i]);
			//printf("prime = %u, roots = %u,%u.  max steps = %u\n",
			//	fbc->prime[i+1],fbc->root1[i+1],fbc->root2[i+1],
			//	(fbc->prime[i+1] - 1 + BLOCKSIZE) / fbc->prime[i+1]);
			//printf("prime = %u, roots = %u,%u.  max steps = %u\n",
			//	fbc->prime[i+2],fbc->root1[i+2],fbc->root2[i+2],
			//	(fbc->prime[i+2] - 1 + BLOCKSIZE) / fbc->prime[i+2]);
			//printf("prime = %u, roots = %u,%u.  max steps = %u\n",
			//	fbc->prime[i+3],fbc->root1[i+3],fbc->root2[i+3],
			//	(fbc->prime[i+3] - 1 + BLOCKSIZE) / fbc->prime[i+3]);

			uint32 result = 0;
			RESIEVE_8X_14BIT_MAX;
			
			if (result == 0)
			{
				i += 8;
				continue;
			}

			if (result & 0x2)
			{
				DIVIDE_RESIEVED_PRIME(0);
			}

			if (result & 0x8)
			{
				DIVIDE_RESIEVED_PRIME(1);
			}

			if (result & 0x20)
			{
				DIVIDE_RESIEVED_PRIME(2);
			}

			if (result & 0x80)
			{
				DIVIDE_RESIEVED_PRIME(3);
			}

			if (result & 0x200)
			{
				DIVIDE_RESIEVED_PRIME(4);
			}

			if (result & 0x800)
			{
				DIVIDE_RESIEVED_PRIME(5);
			}

			if (result & 0x2000)
			{
				DIVIDE_RESIEVED_PRIME(6);
			}

			if (result & 0x8000)
			{
				DIVIDE_RESIEVED_PRIME(7);
			}

			i += 8;
		}
#endif
		bound = sconf->factor_base->fb_15bit_B;

		while ((uint32)i < bound)
		{
			uint32 result = 0;
			RESIEVE_8X_15BIT_MAX;

			if (result == 0)
			{
				i += 8;
				continue;
			}

			if (result & 0x2)
			{
				DIVIDE_RESIEVED_PRIME(0);
			}

			if (result & 0x8)
			{
				DIVIDE_RESIEVED_PRIME(1);
			}

			if (result & 0x20)
			{
				DIVIDE_RESIEVED_PRIME(2);
			}

			if (result & 0x80)
			{
				DIVIDE_RESIEVED_PRIME(3);
			}

			if (result & 0x200)
			{
				DIVIDE_RESIEVED_PRIME(4);
			}

			if (result & 0x800)
			{
				DIVIDE_RESIEVED_PRIME(5);
			}

			if (result & 0x2000)
			{
				DIVIDE_RESIEVED_PRIME(6);
			}

			if (result & 0x8000)
			{
				DIVIDE_RESIEVED_PRIME(7);
			}

			i += 8;
		}

		bound = sconf->factor_base->med_B;
		while ((uint32)i < bound)
		{

			uint32 result = 0;
			RESIEVE_8X_16BIT_MAX;

			if (result == 0)
			{
				i += 8;
				continue;
			}

			if (result & 0x2)
			{
				DIVIDE_RESIEVED_PRIME(0);
			}

			if (result & 0x8)
			{
				DIVIDE_RESIEVED_PRIME(1);
			}

			if (result & 0x20)
			{
				DIVIDE_RESIEVED_PRIME(2);
			}

			if (result & 0x80)
			{
				DIVIDE_RESIEVED_PRIME(3);
			}

			if (result & 0x200)
			{
				DIVIDE_RESIEVED_PRIME(4);
			}

			if (result & 0x800)
			{
				DIVIDE_RESIEVED_PRIME(5);
			}

			if (result & 0x2000)
			{
				DIVIDE_RESIEVED_PRIME(6);
			}

			if (result & 0x8000)
			{
				DIVIDE_RESIEVED_PRIME(7);
			}

			i += 8;
		}

#else
		// standard trial division methods, for when either the compiler or the OS or both
		// cause resieving to be slower.

		bound = sconf->factor_base->med_B;
		while ((uint32)i < bound)
		{
			uint32 tmp;
			uint32 prime = fbc->prime[i];
			uint32 root1 = fbc->root1[i];
			uint32 root2 = fbc->root2[i];

			//after sieving a block, the root is updated for the start of the next block
			//get it back on the current block's progression
			root1 = root1 + BLOCKSIZE - block_loc;		
			root2 = root2 + BLOCKSIZE - block_loc;	

			//there are faster methods if this is the case
			if (prime > BLOCKSIZE)
				break;

			//the difference root - blockoffset is less than prime about 20% of the time
			//and thus we don't have to divide at all in those cases.
			if (root2 >= prime)
			{
				//r2 is bigger than prime, it could be on the progression, check it.
				tmp = root2 + fullfb_ptr->correction[i];
				q64 = (uint64)tmp * (uint64)fullfb_ptr->small_inv[i];
				tmp = q64 >> 40; 
				tmp = root2 - tmp * prime;

				if (tmp == 0)
				{
					//it is, so it will divide Q(x).  do so as many times as we can.
					DIVIDE_ONE_PRIME;
				}
				else if (root1 >= prime)
				{			
					tmp = root1 + fullfb_ptr->correction[i];
					q64 = (uint64)tmp * (uint64)fullfb_ptr->small_inv[i];
					tmp = q64 >> 40; 
					tmp = root1 - tmp * prime;

					if (tmp == 0)
					{
						//r2 was a bust, but root1 met the criteria.  divide Q(x).	
						DIVIDE_ONE_PRIME;
					}
				}
			}
			i++;
		}

		//for primes bigger than the blocksize, we don't need to divide at all since
		//there can be at most one instance of the prime in the block.  thus the 
		//distance to the next one is equal to 'prime', rather than a multiple of it.
		bound = sconf->factor_base->med_B;
		while ((uint32)i < bound)
		{
			uint32 tmp;
			uint32 prime = fbc->prime[i];
			uint32 root1 = fbc->root1[i];
			uint32 root2 = fbc->root2[i];

			//the fbptr roots currently point to the next block
			//so adjust the current index
			tmp = block_loc + prime - BLOCKSIZE;
			if ((root1 == tmp) || (root2 == tmp))
			{
				DIVIDE_ONE_PRIME;
			}
			i++;
		}

#endif

		// after either resieving or standard trial division, record
		// how many factors we've found so far.
		dconf->smooth_num[report_num] = smooth_num;	

	}
			
#ifdef QS_TIMING
	gettimeofday (&qs_timing_stop, NULL);
	qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

	TF_STG4 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
	free(qs_timing_diff);
#endif

	return;
}
