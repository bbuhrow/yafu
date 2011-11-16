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

this file contains code implementing 3) and 4)


*/


#if defined(GCC_ASM32X) || defined(GCC_ASM64X) || defined(__MINGW32__)
	//these compilers support SIMD 
	#define SSE2_RESIEVING 1

	#if defined(HAS_SSE2)

		#define MOD_INIT_8X												\
			ASM_G (														\
				"movdqa (%0), %%xmm0 \n\t"		/* move in BLOCKSIZE */ \
				"movdqa (%1), %%xmm1 \n\t"		/* move in block_loc */ \
				:														\
				: "r" (bl_sizes), "r" (bl_locs)							\
				: "xmm0", "xmm1");

			// the original roots are the current roots + BLOCKSIZE
			// to test if this block_loc is on the root progression, we first
			// advance the current block_loc one step past blocksize and
			// then test if this value is equal to either root

			// to advance the block_loc, we compute 
			// steps = 1 + (BLOCKSIZE - block_loc) / prime

			// to do the div quickly using precomputed values, do:
			// tmp = ((BLOCKSIZE - block_loc) + correction) * inv >> shift
			// shift is either 24 or 26 bits, depending on the size of prime
			// in order to keep enough precision.  See Agner Fog's optimization
			// manuals.
			// also note that with 64k versions, the final addition will overflow,
			// but since the addition does not saturate and since the final 
			// subtraction always yeilds a number less than 2^16, the overflow
			// does not hurt.

#ifdef YAFU_64K
		#define MOD_CMP_8X																		\
			ASM_G (																				\
				"movdqa %%xmm0, %%xmm2 \n\t"	/* copy BLOCKSIZE */							\
				"pcmpeqw	%%xmm6, %%xmm6 \n\t"	/* create an array of 1's */				\
				"movdqa (%1), %%xmm3 \n\t"		/* move in primes */							\
				"psubw	%%xmm1, %%xmm2 \n\t"	/* BLOCKSIZE - block_loc */						\
				"movdqa (%2), %%xmm4 \n\t"		/* move in corrections */						\
				"psrlw	$15, %%xmm6 \n\t"		/* create an array of 1's */					\
				"movdqa %%xmm1, %%xmm7 \n\t"	/* copy block_loc */							\
				"paddw	%%xmm6, %%xmm2 \n\t"	/* add in 1's */								\
				"psubw	%%xmm6, %%xmm7 \n\t"	/* substract 1's */								\
				"movdqa (%3), %%xmm5 \n\t"		/* move in inverses */							\
				"paddw	%%xmm2, %%xmm4 \n\t"	/* apply corrections */							\
				"movdqa (%4), %%xmm6 \n\t"		/* move in root1s */							\
				"pmulhuw	%%xmm5, %%xmm4 \n\t"	/* (unsigned) multiply by inverses */		\
				"movdqa (%5), %%xmm2 \n\t"		/* move in root2s */							\
				"psrlw	$8, %%xmm4 \n\t"		/* to get to total shift of 24 bits */			\
				"paddw	%%xmm3, %%xmm7 \n\t"	/* add primes and block_loc */					\
				"pmullw	%%xmm3, %%xmm4 \n\t"	/* (signed) multiply by primes */				\
				"paddw	%%xmm7, %%xmm4 \n\t"	/* add in block_loc + primes */					\
				"psubw	%%xmm0, %%xmm4 \n\t"	/* substract blocksize */						\
				"pcmpeqw	%%xmm4, %%xmm6 \n\t"	/* compare to root1s */						\
				"pcmpeqw	%%xmm4, %%xmm2 \n\t"	/* compare to root2s */						\
				"por	%%xmm6, %%xmm2 \n\t"	/* combine compares */							\
				"pmovmskb %%xmm2, %0 \n\t"		/* export to result */							\
				: "=r" (tmp3)																		\
				: "r" (fbc->prime + i), "r" (fullfb_ptr->correction + i), "r" (fullfb_ptr->small_inv + i), "r" (fbc->root1 + i), "r" (fbc->root2 + i)	\
				: "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7");

		#define MOD_CMP_BIG_8X																		\
			ASM_G (																				\
				"movdqa %%xmm0, %%xmm2 \n\t"	/* copy BLOCKSIZE */							\
				"pcmpeqw	%%xmm6, %%xmm6 \n\t"	/* create an array of 1's */				\
				"movdqa (%1), %%xmm3 \n\t"		/* move in primes */							\
				"psubw	%%xmm1, %%xmm2 \n\t"	/* BLOCKSIZE - block_loc */						\
				"movdqa (%2), %%xmm4 \n\t"		/* move in corrections */						\
				"psrlw	$15, %%xmm6 \n\t"		/* create an array of 1's */					\
				"movdqa %%xmm1, %%xmm7 \n\t"	/* copy block_loc */							\
				"paddw	%%xmm6, %%xmm2 \n\t"	/* add in 1's */								\
				"psubw	%%xmm6, %%xmm7 \n\t"	/* substract 1's */								\
				"movdqa (%3), %%xmm5 \n\t"		/* move in inverses */							\
				"paddw	%%xmm2, %%xmm4 \n\t"	/* apply corrections */							\
				"movdqa (%4), %%xmm6 \n\t"		/* move in root1s */							\
				"pmulhuw	%%xmm5, %%xmm4 \n\t"	/* (unsigned) multiply by inverses */		\
				"movdqa (%5), %%xmm2 \n\t"		/* move in root2s */							\
				"psrlw	$10, %%xmm4 \n\t"		/* to get to total shift of 26 bits */			\
				"paddw	%%xmm3, %%xmm7 \n\t"	/* add primes and block_loc */					\
				"pmullw	%%xmm3, %%xmm4 \n\t"	/* (signed) multiply by primes */				\
				"paddw	%%xmm7, %%xmm4 \n\t"	/* add in block_loc + primes */					\
				"psubw	%%xmm0, %%xmm4 \n\t"	/* substract blocksize */						\
				"pcmpeqw	%%xmm4, %%xmm6 \n\t"	/* compare to root1s */						\
				"pcmpeqw	%%xmm4, %%xmm2 \n\t"	/* compare to root2s */						\
				"por	%%xmm6, %%xmm2 \n\t"	/* combine compares */							\
				"pmovmskb %%xmm2, %0 \n\t"		/* export to result */							\
				: "=r" (tmp3)																		\
				: "r" (fbc->prime + i), "r" (fullfb_ptr->correction + i), "r" (fullfb_ptr->small_inv + i), "r" (fbc->root1 + i), "r" (fbc->root2 + i)	\
				: "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7");

#else
		#define MOD_CMP_8X																		\
			ASM_G (																				\
				"movdqa (%6), %%xmm0 \n\t"		/* move in BLOCKSIZE */ \
				"movdqa (%7), %%xmm1 \n\t"		/* move in block_loc */ \
				"movdqa %%xmm0, %%xmm2 \n\t"	/* copy BLOCKSIZE */							\
				"movdqa (%1), %%xmm3 \n\t"		/* move in primes */							\
				"psubw	%%xmm1, %%xmm2 \n\t"	/* BLOCKSIZE - block_loc */						\
				"movdqa (%2), %%xmm4 \n\t"		/* move in corrections */							\
				"movdqa %%xmm1, %%xmm7 \n\t"	/* copy block_loc */							\
				"movdqa (%3), %%xmm5 \n\t"		/* move in inverses */							\
				"paddw	%%xmm2, %%xmm4 \n\t"	/* apply corrections */							\
				"movdqa (%4), %%xmm6 \n\t"		/* move in root1s */							\
				"pmulhuw	%%xmm5, %%xmm4 \n\t"	/* (unsigned) multiply by inverses */		\
				"movdqa (%5), %%xmm2 \n\t"		/* move in root2s */							\
				"psrlw	$8, %%xmm4 \n\t"		/* to get to total shift of 24 bits */			\
				"paddw	%%xmm3, %%xmm7 \n\t"	/* add primes and block_loc */					\
				"pmullw	%%xmm3, %%xmm4 \n\t"	/* (signed) multiply by primes */				\
				"paddw	%%xmm7, %%xmm4 \n\t"	/* add in block_loc + primes */					\
				"psubw	%%xmm0, %%xmm4 \n\t"	/* substract blocksize */						\
				"pcmpeqw	%%xmm4, %%xmm6 \n\t"	/* compare to root1s */						\
				"pcmpeqw	%%xmm4, %%xmm2 \n\t"	/* compare to root2s */						\
				"por	%%xmm6, %%xmm2 \n\t"	/* combine compares */								\
				"pmovmskb %%xmm2, %0 \n\t"		/* export to result */							\
				: "=r" (tmp3)																	\
				: "r" (fbc->prime + i), "r" (fullfb_ptr->correction + i), "r" (fullfb_ptr->small_inv + i), "r" (fbc->root1 + i), "r" (fbc->root2 + i), "r" (bl_sizes), "r" (bl_locs)	\
				: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7");

		#define MOD_CMP_BIG_8X																	\
			ASM_G (																				\
				"movdqa (%6), %%xmm0 \n\t"		/* move in BLOCKSIZE */ \
				"movdqa (%7), %%xmm1 \n\t"		/* move in block_loc */ \
				"movdqa %%xmm0, %%xmm2 \n\t"	/* copy BLOCKSIZE */							\
				"movdqa (%1), %%xmm3 \n\t"		/* move in primes */							\
				"psubw	%%xmm1, %%xmm2 \n\t"	/* BLOCKSIZE - block_loc */						\
				"movdqa (%2), %%xmm4 \n\t"		/* move in corrections */							\
				"movdqa %%xmm1, %%xmm7 \n\t"	/* copy block_loc */							\
				"movdqa (%3), %%xmm5 \n\t"		/* move in inverses */							\
				"paddw	%%xmm2, %%xmm4 \n\t"	/* apply corrections */							\
				"movdqa (%4), %%xmm6 \n\t"		/* move in root1s */							\
				"pmulhuw	%%xmm5, %%xmm4 \n\t"	/* (unsigned) multiply by inverses */		\
				"movdqa (%5), %%xmm2 \n\t"		/* move in root2s */							\
				"psrlw	$11, %%xmm4 \n\t"		/* to get to total shift of 26 bits */			\
				"paddw	%%xmm3, %%xmm7 \n\t"	/* add primes and block_loc */					\
				"pmullw	%%xmm3, %%xmm4 \n\t"	/* (signed) multiply by primes */				\
				"paddw	%%xmm7, %%xmm4 \n\t"	/* add in block_loc + primes */					\
				"psubw	%%xmm0, %%xmm4 \n\t"	/* substract blocksize */						\
				"pcmpeqw	%%xmm4, %%xmm6 \n\t"	/* compare to root1s */						\
				"pcmpeqw	%%xmm4, %%xmm2 \n\t"	/* compare to root2s */						\
				"por	%%xmm6, %%xmm2 \n\t"	/* combine compares */								\
				"pmovmskb %%xmm2, %0 \n\t"		/* export to result */							\
				: "=r" (tmp3)																	\
				: "r" (fbc->prime + i), "r" (fullfb_ptr->correction + i), "r" (fullfb_ptr->small_inv + i), "r" (fbc->root1 + i), "r" (fbc->root2 + i), "r" (bl_sizes), "r" (bl_locs)	\
				: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7");
#endif


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


#elif defined(MSC_ASM32A)
	#define SSE2_RESIEVING 1

	#if defined(HAS_SSE2)
		//top level sieve scanning with SSE2

		#define MOD_INIT_8X												\
			do { \
				ASM_M { \
				ASM_M movdqa xmm0, bl_sizes \
				ASM_M movdqa xmm1, bl_locs } \
			} while (0) ;

#ifdef YAFU_64K

		#define MOD_CMP_8X																		\
			do { \
					__m128i one;	\
					__m128i t1; \
					__m128i t2;	\
					__m128i lr1;	\
					__m128i lr2;	\
					__m128i lp;	\
					__m128i c;	\
					__m128i sinv; \
					__m128i blksz; \
					__m128i blkloc; \
				blksz = _mm_load_si128((__m128i *)(bl_sizes)); \
				blkloc = _mm_load_si128((__m128i *)(bl_locs)); \
				t1 = blksz; \
				one = _mm_cmpeq_epi16(one, one); \
				lp = _mm_load_si128((__m128i *)(fbc->prime + i)); \
				t1 = _mm_sub_epi16(t1, blkloc); \
				c = _mm_load_si128((__m128i *)(fullfb_ptr->correction + i));	\
				one = _mm_srli_epi16(one, 15); \
				t2 = blkloc; \
				t1 = _mm_add_epi16(t1, one); \
				t2 = _mm_sub_epi16(t2, one); \
				sinv = _mm_load_si128((__m128i *)(fullfb_ptr->small_inv + i));	\
				c = _mm_add_epi16(c, t1); \
				lr1 = _mm_load_si128((__m128i *)(fbc->root1 + i)); \
				c = _mm_mulhi_epu16(c, sinv); \
				lr2 = _mm_load_si128((__m128i *)(fbc->root2 + i)); \
				c = _mm_srli_epi16(c, 8); \
				t2 = _mm_add_epi16(t2, lp); \
				c = _mm_mullo_epi16(c, lp); \
				c = _mm_add_epi16(c, t2); \
				c = _mm_sub_epi16(c, blksz); \
				lr1 = _mm_cmpeq_epi16(lr1, c); \
				lr2 = _mm_cmpeq_epi16(lr2, c); \
				lr2 = _mm_or_si128(lr2, lr1); \
				tmp3 = _mm_movemask_epi8(lr2); \
			} while (0);


		#define MOD_CMP_BIG_8X																		\
			do { \
					__m128i one;	\
					__m128i t1; \
					__m128i t2;	\
					__m128i lr1;	\
					__m128i lr2;	\
					__m128i lp;	\
					__m128i c;	\
					__m128i sinv; \
					__m128i blksz; \
					__m128i blkloc; \
				blksz = _mm_load_si128((__m128i *)(bl_sizes)); \
				blkloc = _mm_load_si128((__m128i *)(bl_locs)); \
				t1 = blksz; \
				one = _mm_cmpeq_epi16(one, one); \
				lp = _mm_load_si128((__m128i *)(fbc->prime + i)); \
				t1 = _mm_sub_epi16(t1, blkloc); \
				c = _mm_load_si128((__m128i *)(fullfb_ptr->correction + i));	\
				one = _mm_srli_epi16(one, 15); \
				t2 = blkloc; \
				t1 = _mm_add_epi16(t1, one); \
				t2 = _mm_sub_epi16(t2, one); \
				sinv = _mm_load_si128((__m128i *)(fullfb_ptr->small_inv + i));	\
				c = _mm_add_epi16(c, t1); \
				lr1 = _mm_load_si128((__m128i *)(fbc->root1 + i)); \
				c = _mm_mulhi_epu16(c, sinv); \
				lr2 = _mm_load_si128((__m128i *)(fbc->root2 + i)); \
				c = _mm_srli_epi16(c, 8); \
				t2 = _mm_add_epi16(t2, lp); \
				c = _mm_mullo_epi16(c, lp); \
				c = _mm_add_epi16(c, t2); \
				c = _mm_sub_epi16(c, blksz); \
				lr1 = _mm_cmpeq_epi16(lr1, c); \
				lr2 = _mm_cmpeq_epi16(lr2, c); \
				lr2 = _mm_or_si128(lr2, lr1); \
				tmp3 = _mm_movemask_epi8(lr2); \
			} while (0);

#else
		#define MOD_CMP_8X																		\
			do { \
				uint16 *localprime = (uint16 *)(fbc->prime + i);	\
				uint16 *localsmall_inv = (uint16 *)(fullfb_ptr->small_inv + i);	\
				uint16 *localcorrections = (uint16 *)(fullfb_ptr->correction + i);	\
				uint16 *localroot1 = (uint16 *)(fbc->root1 + i);	\
				uint16 *localroot2 = (uint16 *)(fbc->root2 + i);	\
				uint16 *local_blksz = (uint16 *)bl_sizes; \
				uint16 *local_blkloc = (uint16 *)bl_locs; \
				ASM_M { \
				ASM_M mov eax, localprime \
				ASM_M mov ebx, local_blksz \
				ASM_M mov ecx, local_blkloc \
				ASM_M mov edx, localsmall_inv \
				ASM_M mov edi, localcorrections \
				ASM_M movdqa xmm0, XMMWORD PTR [ebx] \
				ASM_M movdqa xmm1, XMMWORD PTR [ecx] \
				ASM_M mov ebx, localroot1 \
				ASM_M mov ecx, localroot2 \
				ASM_M movdqa xmm2, xmm0 \
				ASM_M movdqa xmm3, XMMWORD PTR [eax] \
				ASM_M psubw xmm2, xmm1 \
				ASM_M movdqa xmm4, XMMWORD PTR [edi] \
				ASM_M movdqa xmm7, xmm1 \
				ASM_M movdqa xmm5, XMMWORD PTR [edx] \
				ASM_M paddw xmm4, xmm2 \
				ASM_M movdqa xmm6, XMMWORD PTR [ebx] \
				ASM_M pmulhuw xmm4, xmm5 \
				ASM_M movdqa xmm2, XMMWORD PTR [ecx] \
				ASM_M psrlw xmm4, 8 \
				ASM_M paddw xmm7, xmm3 \
				ASM_M pmullw xmm4, xmm3 \
				ASM_M paddw xmm4, xmm7 \
				ASM_M psubw xmm4, xmm0 \
				ASM_M pcmpeqw xmm6, xmm4 \
				ASM_M pcmpeqw xmm2, xmm4 \
				ASM_M por xmm2, xmm6 \
				ASM_M pmovmskb eax, xmm2 \
				ASM_M mov tmp3, eax } \
			} while (0);

		#define MOD_CMP_BIG_8X																		\
			do { \
				uint16 *localprime = (uint16 *)(fbc->prime + i);	\
				uint16 *localsmall_inv = (uint16 *)(fullfb_ptr->small_inv + i);	\
				uint16 *localcorrections = (uint16 *)(fullfb_ptr->correction + i);	\
				uint16 *localroot1 = (uint16 *)(fbc->root1 + i);	\
				uint16 *localroot2 = (uint16 *)(fbc->root2 + i);	\
				uint16 *local_blksz = (uint16 *)bl_sizes; \
				uint16 *local_blkloc = (uint16 *)bl_locs; \
				ASM_M { \
				ASM_M mov eax, localprime \
				ASM_M mov ebx, local_blksz \
				ASM_M mov ecx, local_blkloc \
				ASM_M mov edx, localsmall_inv \
				ASM_M mov edi, localcorrections \
				ASM_M movdqa xmm0, XMMWORD PTR [ebx] \
				ASM_M movdqa xmm1, XMMWORD PTR [ecx] \
				ASM_M mov ebx, localroot1 \
				ASM_M mov ecx, localroot2 \
				ASM_M movdqa xmm2, xmm0 \
				ASM_M movdqa xmm3, XMMWORD PTR [eax] \
				ASM_M psubw xmm2, xmm1 \
				ASM_M movdqa xmm4, XMMWORD PTR [edi] \
				ASM_M movdqa xmm7, xmm1 \
				ASM_M movdqa xmm5, XMMWORD PTR [edx] \
				ASM_M paddw xmm4, xmm2 \
				ASM_M movdqa xmm6, XMMWORD PTR [ebx] \
				ASM_M pmulhuw xmm4, xmm5 \
				ASM_M movdqa xmm2, XMMWORD PTR [ecx] \
				ASM_M psrlw xmm4, 10 \
				ASM_M paddw xmm7, xmm3 \
				ASM_M pmullw xmm4, xmm3 \
				ASM_M paddw xmm4, xmm7 \
				ASM_M psubw xmm4, xmm0 \
				ASM_M pcmpeqw xmm6, xmm4 \
				ASM_M pcmpeqw xmm2, xmm4 \
				ASM_M por xmm2, xmm6 \
				ASM_M pmovmskb eax, xmm2 \
				ASM_M mov tmp3, eax } \
			} while (0);
#endif

		#define STEP_COMPARE_COMBINE \
			ASM_M psubw xmm2, xmm1 \
			ASM_M psubw xmm3, xmm1 \
			ASM_M pcmpeqw xmm5, xmm2 \
			ASM_M pcmpeqw xmm6, xmm3 \
			ASM_M por xmm7, xmm5 \
			ASM_M por xmm8, xmm6


		#define INIT_RESIEVE \
			ASM_M movdqa xmm4, XMMWORD PTR [edx] \
			ASM_M pxor xmm8, xmm8 \
			ASM_M movdqa xmm2, XMMWORD PTR [ebx] \
			ASM_M paddw xmm2, xmm4 \
			ASM_M movdqa xmm3, XMMWORD PTR [ecx] \
			ASM_M paddw xmm3, xmm4 \
			ASM_M movdqa xmm1, XMMWORD PTR [eax] \
			ASM_M pxor xmm7, xmm7 \
			ASM_M pxor xmm6, xmm6 \
			ASM_M pxor xmm5, xmm5

		#ifdef YAFU_64K

			#define RESIEVE_8X_14BIT_MAX \
				do { \
					uint32 *localprime = (uint32 *)(fbc->prime + i);	\
					uint32 *localroot1 = (uint32 *)(fbc->root1 + i);	\
					uint32 *localroot2 = (uint32 *)(fbc->root2 + i);	\
					uint32 *localcorrect = (uint32 *)(corrections);	\
					ASM_M { \
					ASM_M mov eax, localprime \
					ASM_M mov ebx, localroot1 \
					ASM_M mov ecx, localroot2 \
					ASM_M mov edx, localcorrect \
					INIT_RESIEVE \
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					ASM_M por xmm7, xmm8 \
					ASM_M pmovmskb eax, xmm7	\
					ASM_M mov result, eax } \
				} while (0);

			#define RESIEVE_8X_15BIT_MAX \
				do { \
					uint32 *localprime = (uint32 *)(fbc->prime + i);	\
					uint32 *localroot1 = (uint32 *)(fbc->root1 + i);	\
					uint32 *localroot2 = (uint32 *)(fbc->root2 + i);	\
					uint32 *localcorrect = (uint32 *)(corrections);	\
					ASM_M { \
					ASM_M mov eax, localprime \
					ASM_M mov ebx, localroot1 \
					ASM_M mov ecx, localroot2 \
					ASM_M mov edx, localcorrect \
					INIT_RESIEVE \
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					ASM_M por xmm7, xmm8 \
					ASM_M pmovmskb eax, xmm7	\
					ASM_M mov result, eax } \
				} while (0);

			#define RESIEVE_8X_16BIT_MAX \
				do { \
					uint32 *localprime = (uint32 *)(fbc->prime + i);	\
					uint32 *localroot1 = (uint32 *)(fbc->root1 + i);	\
					uint32 *localroot2 = (uint32 *)(fbc->root2 + i);	\
					uint32 *localcorrect = (uint32 *)(corrections);	\
					ASM_M { \
					ASM_M mov eax, localprime \
					ASM_M mov ebx, localroot1 \
					ASM_M mov ecx, localroot2 \
					ASM_M mov edx, localcorrect \
					INIT_RESIEVE \
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					ASM_M por xmm7, xmm8 \
					ASM_M pmovmskb eax, xmm7	\
					ASM_M mov result, eax } \
				} while (0);

		#else

			#define RESIEVE_8X_14BIT_MAX \
				do { \
					uint16 *localprime = fbc->prime + i;	\
					uint16 *localroot1 = fbc->root1 + i;	\
					uint16 *localroot2 = fbc->root2 + i;	\
					ASM_M { \
					ASM_M mov eax, localprime \
					ASM_M mov ebx, localroot1 \
					ASM_M mov ecx, localroot2 \
					INIT_RESIEVE \
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					ASM_M por xmm7, xmm8 \
					ASM_M pmovmskb eax, xmm7	\
					ASM_M mov result, eax } \
				} while (0);

			#define RESIEVE_8X_15BIT_MAX \
				do { \
					uint16 *localprime = fbc->prime + i;	\
					uint16 *localroot1 = fbc->root1 + i;	\
					uint16 *localroot2 = fbc->root2 + i;	\
					ASM_M { \
					ASM_M mov eax, localprime \
					ASM_M mov ebx, localroot1 \
					ASM_M mov ecx, localroot2 \
					INIT_RESIEVE \
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					ASM_M por xmm7, xmm8 \
					ASM_M pmovmskb eax, xmm7	\
					ASM_M mov result, eax } \
				} while (0);

			#define RESIEVE_8X_16BIT_MAX \
				do { \
					uint16 *localprime = fbc->prime + i;	\
					uint16 *localroot1 = fbc->root1 + i;	\
					uint16 *localroot2 = fbc->root2 + i;	\
					ASM_M { \
					ASM_M mov eax, localprime \
					ASM_M mov ebx, localroot1 \
					ASM_M mov ecx, localroot2 \
					INIT_RESIEVE \
					STEP_COMPARE_COMBINE	\
					ASM_M por xmm7, xmm8 \
					ASM_M pmovmskb eax, xmm7	\
					ASM_M mov result, eax } \
				} while (0);

		#endif

	#elif defined(HAS_MMX)



	#else
		#error SSE2 is required
	#endif

#elif defined(_WIN64)

	#define SSE2_RESIEVING 1

	#if defined(HAS_SSE2)

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

void filter_medprimes(uint8 parity, uint32 poly_id, uint32 bnum, 
						 static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//we have flagged this sieve offset as likely to produce a relation
	//nothing left to do now but check and see.
	uint64 q64;
	int i;
	uint32 bound, tmp, prime, root1, root2, report_num;
	int smooth_num;
	uint32 *fb_offsets;
	sieve_fb *fb;
	sieve_fb_compressed *fbc;
	fb_element_siqs *fullfb_ptr, *fullfb = sconf->factor_base->list;
	uint32 block_loc;
#if defined(USE_8X_MOD) && !defined(USE_8X_MOD_ASM)
	uint16 tmp1, tmp2, tmp3, tmp4;
#else
	uint32 tmp1, tmp2, tmp3, tmp4;
#endif

	uint16 *corrections = dconf->corrections;
#ifdef USE_COMPRESSED_FB
	sieve_fb_compressed *fbptr;
#endif

#ifdef USE_8X_MOD_ASM
	uint16 *bl_sizes;
	uint16 *bl_locs;

	bl_sizes = (uint16 *)xmalloc_align(8 * sizeof(uint16));
	bl_locs = (uint16 *)xmalloc_align(8 * sizeof(uint16));

#endif

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

	for (report_num = 0; report_num < dconf->num_reports; report_num++)
	{
#ifdef USE_YAFU_TDIV
		z32 *tmp32 = &dconf->Qvals32[report_num];
#endif

		if (!dconf->valid_Qs[report_num])
			continue;

		fb_offsets = &dconf->fb_offsets[report_num][0];
		smooth_num = dconf->smooth_num[report_num];
		block_loc = dconf->reports[report_num];

#ifdef QS_TIMING
		gettimeofday(&qs_timing_start, NULL);
#endif

		//do the primes less than the blocksize.  primes bigger than the blocksize can be handled
		//even more efficiently.
		//a couple of observations from jasonp:
		//if a prime divides Q(x), then this index (j) and either
		//root1 or root2 are on the same arithmetic progression.  this we can
		//test with a single precision mod operation
		//do the first few until the rest can be done in batches of 4 that are aligned to 16 byte
		//boundaries.  this is necessary to use the SSE2 batch mod code, if it ever
		//becomes faster...

		i=sconf->sieve_small_fb_start;

#if !defined(USE_RESIEVING)
		bound = sconf->factor_base->small_B;
#else
		bound = sconf->factor_base->fb_11bit_B;
#endif
		
		// single-up test until i is a multiple of 8
		while ((uint32)i < bound && ((i & 7) != 0))
		{
#ifdef USE_COMPRESSED_FB
			fbptr = fbc + i;
			prime = fbptr->prime_and_logp & 0xFFFF;
			root1 = fbptr->roots & 0xFFFF;
			root2 = fbptr->roots >> 16;
#else
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
#endif
			//tmp = distance from this sieve block offset to the end of the block
			tmp = BLOCKSIZE - block_loc;
	
			//tmp = tmp/prime + 1 = number of steps to get past the end of the sieve
			//block, which is the state of the sieve now.
			tmp = 1+(uint32)(((uint64)(tmp + fullfb_ptr->correction[i])
					* (uint64)fullfb_ptr->small_inv[i]) >> FOGSHIFT); 
			tmp = block_loc + tmp*prime;
			tmp = tmp - BLOCKSIZE;

			//tmp = advance the offset to where it should be after the interval, and
			//check to see if that's where either of the roots are now.  if so, then
			//this offset is on the progression of the sieve for this prime
			if (tmp == root1 || tmp == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;
		}

#ifdef USE_8X_MOD_ASM

		tmp1 = (BLOCKSIZE << 16) | BLOCKSIZE;
		tmp2 = (block_loc << 16) | block_loc;

#ifdef YAFU_64K
		// the 64k blocksize will overflow the 16 bit word
		// so we use blocksize-1 and then compensate in 
		// the assembly by adding back 1 when needed
		bl_sizes[0] = BLOCKSIZEm1;
		bl_sizes[1] = BLOCKSIZEm1;
		bl_sizes[2] = BLOCKSIZEm1;
		bl_sizes[3] = BLOCKSIZEm1;
		bl_sizes[4] = BLOCKSIZEm1;
		bl_sizes[5] = BLOCKSIZEm1;
		bl_sizes[6] = BLOCKSIZEm1;
		bl_sizes[7] = BLOCKSIZEm1;
#else
		bl_sizes[0] = BLOCKSIZE;
		bl_sizes[1] = BLOCKSIZE;
		bl_sizes[2] = BLOCKSIZE;
		bl_sizes[3] = BLOCKSIZE;
		bl_sizes[4] = BLOCKSIZE;
		bl_sizes[5] = BLOCKSIZE;
		bl_sizes[6] = BLOCKSIZE;
		bl_sizes[7] = BLOCKSIZE;
#endif
		bl_locs[0] = block_loc;
		bl_locs[1] = block_loc;
		bl_locs[2] = block_loc;
		bl_locs[3] = block_loc;
		bl_locs[4] = block_loc;
		bl_locs[5] = block_loc;
		bl_locs[6] = block_loc;
		bl_locs[7] = block_loc;
		
		//MOD_INIT_8X;

		while ((uint32)i < bound)
		{
			tmp3 = 0;

			MOD_CMP_8X;

			if (tmp3 == 0)
			{
				i += 8;
				continue;
			}
			
			if (tmp3 & 0x2)
			{
				DIVIDE_RESIEVED_PRIME(0);
			}

			if (tmp3 & 0x8)
			{
				DIVIDE_RESIEVED_PRIME(1);
			}

			if (tmp3 & 0x20)
			{
				DIVIDE_RESIEVED_PRIME(2);
			}

			if (tmp3 & 0x80)
			{
				DIVIDE_RESIEVED_PRIME(3);
			}

			if (tmp3 & 0x200)
			{
				DIVIDE_RESIEVED_PRIME(4);
			}

			if (tmp3 & 0x800)
			{
				DIVIDE_RESIEVED_PRIME(5);
			}

			if (tmp3 & 0x2000)
			{
				DIVIDE_RESIEVED_PRIME(6);
			}

			if (tmp3 & 0x8000)
			{
				DIVIDE_RESIEVED_PRIME(7);
			}

			i += 8;			

		}

#if defined(GCC_ASM32X) || defined(GCC_ASM64X) || defined(__MINGW32__)
		asm volatile("emms");
#endif

#else
		//now do things in batches of 4 which are aligned on 16 byte boundaries.
		while ((uint32)i < bound)
		{

#if USE_8X_MOD
			tmp1 = BLOCKSIZE - block_loc;
			tmp2 = BLOCKSIZE - block_loc;
			tmp3 = BLOCKSIZE - block_loc;
			tmp4 = BLOCKSIZE - block_loc;

			tmp1 = tmp1 + fullfb_ptr->correction[i];
			q64 = (uint64)tmp1 * (uint64)fullfb_ptr->small_inv[i];
			tmp1 = q64 >> FOGSHIFT; 
			tmp1 = tmp1 + 1;
			tmp1 = block_loc + tmp1 * fullfb_ptr->prime[i];
			
			tmp2 = tmp2 + fullfb_ptr->correction[i+1];
			q64 = (uint64)tmp2 * (uint64)fullfb_ptr->small_inv[i+1];
			tmp2 = q64 >> FOGSHIFT; 
			tmp2 = tmp2 + 1;
			tmp2 = block_loc + tmp2 * fullfb_ptr->prime[i+1];

			tmp3 = tmp3 + fullfb_ptr->correction[i+2];
			q64 = (uint64)tmp3 * (uint64)fullfb_ptr->small_inv[i+2];
			tmp3 = q64 >> FOGSHIFT; 
			tmp3 = tmp3 + 1;
			tmp3 = block_loc + tmp3 * fullfb_ptr->prime[i+2];

			tmp4 = tmp4 + fullfb_ptr->correction[i+3];
			q64 = (uint64)tmp4 * (uint64)fullfb_ptr->small_inv[i+3];
			tmp4 = q64 >>  FOGSHIFT; 
			tmp4 = tmp4 + 1;
			tmp4 = block_loc + tmp4 * fullfb_ptr->prime[i+3];

			tmp1 = tmp1 - BLOCKSIZE;
			tmp2 = tmp2 - BLOCKSIZE;
			tmp3 = tmp3 - BLOCKSIZE;
			tmp4 = tmp4 - BLOCKSIZE;

#else

			tmp1 = BLOCKSIZE - block_loc;
			tmp2 = BLOCKSIZE - block_loc;
			tmp3 = BLOCKSIZE - block_loc;
			tmp4 = BLOCKSIZE - block_loc;

			tmp1 = tmp1 + fullfb_ptr->correction[i];
			q64 = (uint64)tmp1 * (uint64)fullfb_ptr->small_inv[i];
			tmp1 = q64 >> FOGSHIFT; 
			tmp1 = tmp1 + 1;
			tmp1 = block_loc + tmp1 * fullfb_ptr->prime[i];
			
			tmp2 = tmp2 + fullfb_ptr->correction[i+1];
			q64 = (uint64)tmp2 * (uint64)fullfb_ptr->small_inv[i+1];
			tmp2 = q64 >> FOGSHIFT; 
			tmp2 = tmp2 + 1;
			tmp2 = block_loc + tmp2 * fullfb_ptr->prime[i+1];

			tmp3 = tmp3 + fullfb_ptr->correction[i+2];
			q64 = (uint64)tmp3 * (uint64)fullfb_ptr->small_inv[i+2];
			tmp3 = q64 >> FOGSHIFT; 
			tmp3 = tmp3 + 1;
			tmp3 = block_loc + tmp3 * fullfb_ptr->prime[i+2];

			tmp4 = tmp4 + fullfb_ptr->correction[i+3];
			q64 = (uint64)tmp4 * (uint64)fullfb_ptr->small_inv[i+3];
			tmp4 = q64 >>  FOGSHIFT; 
			tmp4 = tmp4 + 1;
			tmp4 = block_loc + tmp4 * fullfb_ptr->prime[i+3];

			tmp1 = tmp1 - BLOCKSIZE;
			tmp2 = tmp2 - BLOCKSIZE;
			tmp3 = tmp3 - BLOCKSIZE;
			tmp4 = tmp4 - BLOCKSIZE;

#endif


#ifdef USE_COMPRESSED_FB
			fbptr = fbc + i;
			prime = fbptr->prime_and_logp & 0xFFFF;
			root1 = fbptr->roots & 0xFFFF;
			root2 = fbptr->roots >> 16;
#else
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
#endif

			if (tmp1 == root1 || tmp1 == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;

#ifdef USE_COMPRESSED_FB
			fbptr = fbc + i;
			prime = fbptr->prime_and_logp & 0xFFFF;
			root1 = fbptr->roots & 0xFFFF;
			root2 = fbptr->roots >> 16;
#else
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
#endif

			if (tmp2 == root1 || tmp2 == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;

#ifdef USE_COMPRESSED_FB
			fbptr = fbc + i;
			prime = fbptr->prime_and_logp & 0xFFFF;
			root1 = fbptr->roots & 0xFFFF;
			root2 = fbptr->roots >> 16;
#else
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
#endif

			if (tmp3 == root1 || tmp3 == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;

#ifdef USE_COMPRESSED_FB
			fbptr = fbc + i;
			prime = fbptr->prime_and_logp & 0xFFFF;
			root1 = fbptr->roots & 0xFFFF;
			root2 = fbptr->roots >> 16;
#else
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
#endif

			if (tmp4 == root1 || tmp4 == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;

		}

		//now cleanup any that don't fit in the last batch
		while ((uint32)i < bound)
		{
#ifdef USE_COMPRESSED_FB
			fbptr = fbc + i;
			prime = fbptr->prime_and_logp & 0xFFFF;
			root1 = fbptr->roots & 0xFFFF;
			root2 = fbptr->roots >> 16;
#else
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
#endif

			tmp = BLOCKSIZE - block_loc;
			tmp = 1+(uint32)(((uint64)(tmp + fullfb_ptr->correction[i])
					* (uint64)fullfb_ptr->small_inv[i]) >> FOGSHIFT); 
			tmp = block_loc + tmp*prime;
			tmp = tmp - BLOCKSIZE;

			if (tmp == root1 || tmp == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;
		}

#endif

		
#if !defined(USE_RESIEVING)
		bound = sconf->factor_base->small_B;
#elif defined(YAFU_64K)
		bound = sconf->factor_base->fb_14bit_B;
#else
		bound = sconf->factor_base->fb_13bit_B;
#endif				

#ifdef USE_8X_MOD_ASM

		tmp1 = (BLOCKSIZE << 16) | BLOCKSIZE;
		tmp2 = (block_loc << 16) | block_loc;
#ifdef YAFU_64K
		bl_sizes[0] = BLOCKSIZEm1;
		bl_sizes[1] = BLOCKSIZEm1;
		bl_sizes[2] = BLOCKSIZEm1;
		bl_sizes[3] = BLOCKSIZEm1;
		bl_sizes[4] = BLOCKSIZEm1;
		bl_sizes[5] = BLOCKSIZEm1;
		bl_sizes[6] = BLOCKSIZEm1;
		bl_sizes[7] = BLOCKSIZEm1;
#else
		bl_sizes[0] = BLOCKSIZE;
		bl_sizes[1] = BLOCKSIZE;
		bl_sizes[2] = BLOCKSIZE;
		bl_sizes[3] = BLOCKSIZE;
		bl_sizes[4] = BLOCKSIZE;
		bl_sizes[5] = BLOCKSIZE;
		bl_sizes[6] = BLOCKSIZE;
		bl_sizes[7] = BLOCKSIZE;
#endif
		bl_locs[0] = block_loc;
		bl_locs[1] = block_loc;
		bl_locs[2] = block_loc;
		bl_locs[3] = block_loc;
		bl_locs[4] = block_loc;
		bl_locs[5] = block_loc;
		bl_locs[6] = block_loc;
		bl_locs[7] = block_loc;

		//MOD_INIT_8X;

		while ((uint32)i < bound)
		{
			tmp3 = 0;

			MOD_CMP_BIG_8X;

			if (tmp3 == 0)
			{
				i += 8;
				continue;
			}
			
			if (tmp3 & 0x2)
			{
				DIVIDE_RESIEVED_PRIME(0);
			}

			if (tmp3 & 0x8)
			{
				DIVIDE_RESIEVED_PRIME(1);
			}

			if (tmp3 & 0x20)
			{
				DIVIDE_RESIEVED_PRIME(2);
			}

			if (tmp3 & 0x80)
			{
				DIVIDE_RESIEVED_PRIME(3);
			}

			if (tmp3 & 0x200)
			{
				DIVIDE_RESIEVED_PRIME(4);
			}

			if (tmp3 & 0x800)
			{
				DIVIDE_RESIEVED_PRIME(5);
			}

			if (tmp3 & 0x2000)
			{
				DIVIDE_RESIEVED_PRIME(6);
			}

			if (tmp3 & 0x8000)
			{
				DIVIDE_RESIEVED_PRIME(7);
			}

			i += 8;			

		}

#if defined(GCC_ASM32X) || defined(GCC_ASM64X) || defined(__MINGW32__)
		asm volatile("emms");
#endif

#else

		// single-up test until i is a multiple of 8
		while ((uint32)i < bound && ((i & 7) != 0))
		{
#ifdef USE_COMPRESSED_FB
			fbptr = fbc + i;
			prime = fbptr->prime_and_logp & 0xFFFF;
			root1 = fbptr->roots & 0xFFFF;
			root2 = fbptr->roots >> 16;
#else
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
#endif
			//tmp = distance from this sieve block offset to the end of the block
			tmp = BLOCKSIZE - block_loc;
	
			//tmp = tmp/prime + 1 = number of steps to get past the end of the sieve
			//block, which is the state of the sieve now.
			tmp = 1+(uint32)(((uint64)(tmp + fullfb_ptr->correction[i])
					* (uint64)fullfb_ptr->small_inv[i]) >> FOGSHIFT_2); 
			tmp = block_loc + tmp*prime;
			tmp = tmp - BLOCKSIZE;

			//tmp = advance the offset to where it should be after the interval, and
			//check to see if that's where either of the roots are now.  if so, then
			//this offset is on the progression of the sieve for this prime
			if (tmp == root1 || tmp == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;
		}

		//now do things in batches of 4 which are aligned on 16 byte boundaries.
		while ((uint32)i < bound)
		{

#if USE_8X_MOD

			tmp1 = BLOCKSIZE - block_loc;
			tmp2 = BLOCKSIZE - block_loc;
			tmp3 = BLOCKSIZE - block_loc;
			tmp4 = BLOCKSIZE - block_loc;

			tmp1 = tmp1 + fullfb_ptr->correction[i];
			q64 = (uint64)tmp1 * (uint64)fullfb_ptr->small_inv[i];
			tmp1 = q64 >> FOGSHIFT_2; 
			tmp1 = tmp1 + 1;
			tmp1 = block_loc + tmp1 * fullfb_ptr->prime[i];
			
			tmp2 = tmp2 + fullfb_ptr->correction[i+1];
			q64 = (uint64)tmp2 * (uint64)fullfb_ptr->small_inv[i+1];
			tmp2 = q64 >> FOGSHIFT_2; 
			tmp2 = tmp2 + 1;
			tmp2 = block_loc + tmp2 * fullfb_ptr->prime[i+1];

			tmp3 = tmp3 + fullfb_ptr->correction[i+2];
			q64 = (uint64)tmp3 * (uint64)fullfb_ptr->small_inv[i+2];
			tmp3 = q64 >> FOGSHIFT_2; 
			tmp3 = tmp3 + 1;
			tmp3 = block_loc + tmp3 * fullfb_ptr->prime[i+2];

			tmp4 = tmp4 + fullfb_ptr->correction[i+3];
			q64 = (uint64)tmp4 * (uint64)fullfb_ptr->small_inv[i+3];
			tmp4 = q64 >>  FOGSHIFT_2; 
			tmp4 = tmp4 + 1;
			tmp4 = block_loc + tmp4 * fullfb_ptr->prime[i+3];

			tmp1 = tmp1 - BLOCKSIZE;
			tmp2 = tmp2 - BLOCKSIZE;
			tmp3 = tmp3 - BLOCKSIZE;
			tmp4 = tmp4 - BLOCKSIZE;

#else

			tmp1 = BLOCKSIZE - block_loc;
			tmp2 = BLOCKSIZE - block_loc;
			tmp3 = BLOCKSIZE - block_loc;
			tmp4 = BLOCKSIZE - block_loc;

			tmp1 = tmp1 + fullfb_ptr->correction[i];
			q64 = (uint64)tmp1 * (uint64)fullfb_ptr->small_inv[i];
			tmp1 = q64 >> FOGSHIFT; 
			tmp1 = tmp1 + 1;
			tmp1 = block_loc + tmp1 * fullfb_ptr->prime[i];
			
			tmp2 = tmp2 + fullfb_ptr->correction[i+1];
			q64 = (uint64)tmp2 * (uint64)fullfb_ptr->small_inv[i+1];
			tmp2 = q64 >> FOGSHIFT; 
			tmp2 = tmp2 + 1;
			tmp2 = block_loc + tmp2 * fullfb_ptr->prime[i+1];

			tmp3 = tmp3 + fullfb_ptr->correction[i+2];
			q64 = (uint64)tmp3 * (uint64)fullfb_ptr->small_inv[i+2];
			tmp3 = q64 >> FOGSHIFT; 
			tmp3 = tmp3 + 1;
			tmp3 = block_loc + tmp3 * fullfb_ptr->prime[i+2];

			tmp4 = tmp4 + fullfb_ptr->correction[i+3];
			q64 = (uint64)tmp4 * (uint64)fullfb_ptr->small_inv[i+3];
			tmp4 = q64 >>  FOGSHIFT; 
			tmp4 = tmp4 + 1;
			tmp4 = block_loc + tmp4 * fullfb_ptr->prime[i+3];

			tmp1 = tmp1 - BLOCKSIZE;
			tmp2 = tmp2 - BLOCKSIZE;
			tmp3 = tmp3 - BLOCKSIZE;
			tmp4 = tmp4 - BLOCKSIZE;

#endif


#ifdef USE_COMPRESSED_FB
			fbptr = fbc + i;
			prime = fbptr->prime_and_logp & 0xFFFF;
			root1 = fbptr->roots & 0xFFFF;
			root2 = fbptr->roots >> 16;
#else
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
#endif

			if (tmp1 == root1 || tmp1 == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;

#ifdef USE_COMPRESSED_FB
			fbptr = fbc + i;
			prime = fbptr->prime_and_logp & 0xFFFF;
			root1 = fbptr->roots & 0xFFFF;
			root2 = fbptr->roots >> 16;
#else
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
#endif

			if (tmp2 == root1 || tmp2 == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;

#ifdef USE_COMPRESSED_FB
			fbptr = fbc + i;
			prime = fbptr->prime_and_logp & 0xFFFF;
			root1 = fbptr->roots & 0xFFFF;
			root2 = fbptr->roots >> 16;
#else
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
#endif

			if (tmp3 == root1 || tmp3 == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;

#ifdef USE_COMPRESSED_FB
			fbptr = fbc + i;
			prime = fbptr->prime_and_logp & 0xFFFF;
			root1 = fbptr->roots & 0xFFFF;
			root2 = fbptr->roots >> 16;
#else
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
#endif

			if (tmp4 == root1 || tmp4 == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;

		}
#endif

		//now cleanup any that don't fit in the last batch
		while ((uint32)i < bound)
		{
#ifdef USE_COMPRESSED_FB
			fbptr = fbc + i;
			prime = fbptr->prime_and_logp & 0xFFFF;
			root1 = fbptr->roots & 0xFFFF;
			root2 = fbptr->roots >> 16;
#else
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
#endif

			tmp = BLOCKSIZE - block_loc;
			tmp = 1+(uint32)(((uint64)(tmp + fullfb_ptr->correction[i])
					* (uint64)fullfb_ptr->small_inv[i]) >> FOGSHIFT_2); 
			tmp = block_loc + tmp*prime;
			tmp = tmp - BLOCKSIZE;

			if (tmp == root1 || tmp == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;
		}

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		TF_STG2 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);

		gettimeofday(&qs_timing_start, NULL);
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
		

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		TF_STG4 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);
#endif

		dconf->smooth_num[report_num] = smooth_num;	

#else
		// standard trial division methods, for when either the compiler or the OS or both
		// cause resieving to be slower.

		bound = sconf->factor_base->med_B;
		while ((uint32)i < bound)
		{
#ifdef USE_COMPRESSED_FB
			fbptr = fbc + i;
			prime = fbptr->prime_and_logp & 0xFFFF;
			root1 = fbptr->roots & 0xFFFF;
			root2 = fbptr->roots >> 16;
#else
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
#endif

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

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		TF_STG3 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);

		gettimeofday(&qs_timing_start, NULL);
#endif

		//for primes bigger than the blocksize, we don't need to divide at all since
		//there can be at most one instance of the prime in the block.  thus the 
		//distance to the next one is equal to 'prime', rather than a multiple of it.
		bound = sconf->factor_base->med_B;
		while ((uint32)i < bound)
		{
#ifdef USE_COMPRESSED_FB
			fbptr = fbc + i;
			prime = fbptr->prime_and_logp & 0xFFFF;
			root1 = fbptr->roots & 0xFFFF;
			root2 = fbptr->roots >> 16;
#else
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
#endif

			//the fbptr roots currently point to the next block
			//so adjust the current index
			tmp = block_loc + prime - BLOCKSIZE;
			if ((root1 == tmp) || (root2 == tmp))
			{
				DIVIDE_ONE_PRIME;
			}
			i++;
		}

		dconf->smooth_num[report_num] = smooth_num;	

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		TF_STG4 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);
#endif

#endif

	}

#ifdef USE_8X_MOD_ASM
	align_free(bl_sizes);
	align_free(bl_locs);
#endif

	return;
}


