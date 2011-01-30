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


#if defined(GCC_ASM32X) || defined(GCC_ASM64X) || defined(__MINGW32__)
	//these compilers support SIMD 
	#define SIMD_SIEVE_SCAN 1
	#define SCAN_CLEAN asm volatile("emms");	
	#define SSE2_RESIEVING 1

	#if defined(HAS_SSE2)
		//top level sieve scanning with SSE2
		#define SIEVE_SCAN_32	\
			asm volatile (		\
				"movdqa (%1), %%xmm0   \n\t"		\
				"orpd 16(%1), %%xmm0    \n\t"		\
				"pmovmskb %%xmm0, %0   \n\t"		\
				: "=r"(result)						\
				: "r"(sieveblock + j), "0"(result)	\
				: "%xmm0");

		#define SIEVE_SCAN_64		\
			asm volatile (							\
				"movdqa (%1), %%xmm0   \n\t"		\
				"orpd 16(%1), %%xmm0    \n\t"		\
				"orpd 32(%1), %%xmm0    \n\t"		\
				"orpd 48(%1), %%xmm0    \n\t"		\
				"pmovmskb %%xmm0, %0   \n\t"		\
				: "=r"(result)						\
				: "r"(sieveblock + j), "0"(result)	\
				: "%xmm0");

		#define SIEVE_SCAN_128		\
			asm volatile (			\
				"movdqa (%1), %%xmm0   \n\t"		\
				"orpd 16(%1), %%xmm0    \n\t"		\
				"orpd 32(%1), %%xmm0    \n\t"		\
				"orpd 48(%1), %%xmm0    \n\t"		\
				"orpd 64(%1), %%xmm0    \n\t"		\
				"orpd 80(%1), %%xmm0    \n\t"		\
				"orpd 96(%1), %%xmm0    \n\t"		\
				"orpd 112(%1), %%xmm0    \n\t"		\
				"pmovmskb %%xmm0, %0   \n\t"		\
				: "=r"(result)						\
				: "r"(sieveblock + j), "0"(result)	\
				: "%xmm0");

		#define SCAN_16X			\
			asm volatile (			\
				"movdqa (%2), %%xmm0	\n\t"	/*move mask into xmm0*/	\
				"movdqa (%1), %%xmm1	\n\t"	/*move 16 bptr locations into xmm regs*/	\
				"movdqa 16(%1), %%xmm2	\n\t"		\
				"movdqa 32(%1), %%xmm3	\n\t"		\
				"movdqa 48(%1), %%xmm4	\n\t"		\
				"pcmpeqw %%xmm0, %%xmm1	\n\t"	/*compare to mask*/	\
				"pcmpeqw %%xmm0, %%xmm2	\n\t"		\
				"pcmpeqw %%xmm0, %%xmm3	\n\t"		\
				"pcmpeqw %%xmm0, %%xmm4	\n\t"		\
				"por %%xmm1, %%xmm4		\n\t"	/*or the comparisons*/	\
				"por %%xmm2, %%xmm3		\n\t"		\
				"por %%xmm3, %%xmm4		\n\t"		\
				"pmovmskb %%xmm4, %0	\n\t"	/*if any are equal, this will be !0*/	\
				: "=r"(result)		\
				: "r"(bptr + j), "r"(mask)			\
				: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4");	

		#define STEP_COMPARE_COMBINE \
			"psubw %%xmm1, %%xmm2 \n\t"		/* subtract primes from root1s */ \
			"psubw %%xmm1, %%xmm3 \n\t"		/* subtract primes from root2s */ \
			"pcmpeqw %%xmm2, %%xmm5 \n\t"	/* root1s ?= 0 */ \
			"pcmpeqw %%xmm3, %%xmm6 \n\t"	/* root2s ?= 0 */ \
			"por %%xmm5, %%xmm7 \n\t"		/* combine results */ \
			"por %%xmm6, %%xmm7 \n\t"		/* combine results */

		#define INIT_RESIEVE \
			"movdqa (%4), %%xmm4 \n\t"		/* bring in corrections to roots */				\
			"movdqa (%2), %%xmm2 \n\t"		/* bring in 8 root1s */ \
			"paddw %%xmm4, %%xmm2 \n\t"		/* correct root1s */ \
			"movdqa (%3), %%xmm3 \n\t"		/* bring in 8 root2s */ \
			"paddw %%xmm4, %%xmm3 \n\t"		/* correct root2s */ \
			"movdqa (%1), %%xmm1 \n\t"		/* bring in 8 primes */ \
			"pxor %%xmm7, %%xmm7 \n\t"		/* zero xmm7 */ \
			"pxor %%xmm5, %%xmm5 \n\t"		/* zero xmm5 */ \
			"pxor %%xmm6, %%xmm6 \n\t"		/* zero xmm6 */

		#ifdef YAFU_64K
			#define RESIEVE_4X_14BIT_MAX \
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
					"pmovmskb %%xmm7, %0 \n\t"		/* if one of these primes divides this location, this will be !0*/ \
					: "=r"(result) \
					: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections) \
					: "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "cc", "memory" \
					);

			#define RESIEVE_4X_15BIT_MAX \
				asm ( \
					INIT_RESIEVE \
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					"pmovmskb %%xmm7, %0 \n\t"		/* if one of these primes divides this location, this will be !0*/ \
					: "=r"(result) \
					: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections) \
					: "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "cc", "memory" \
					);

			#define RESIEVE_4X_16BIT_MAX \
				asm ( \
					INIT_RESIEVE \
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					"pmovmskb %%xmm7, %0 \n\t"		/* if one of these primes divides this location, this will be !0*/ \
					: "=r"(result) \
					: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections) \
					: "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "cc", "memory" \
					);
		#else
			#define RESIEVE_4X_14BIT_MAX \
				asm ( \
					INIT_RESIEVE \
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					"pmovmskb %%xmm7, %0 \n\t"		/* if one of these primes divides this location, this will be !0*/ \
					: "=r"(result) \
					: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections) \
					: "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "cc", "memory" \
					);

			#define RESIEVE_4X_15BIT_MAX \
				asm ( \
					INIT_RESIEVE \
					STEP_COMPARE_COMBINE	\
					STEP_COMPARE_COMBINE	\
					"pmovmskb %%xmm7, %0 \n\t"		/* if one of these primes divides this location, this will be !0*/ \
					: "=r"(result) \
					: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections) \
					: "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "cc", "memory" \
					);

			#define RESIEVE_4X_16BIT_MAX \
				asm ( \
					INIT_RESIEVE \
					STEP_COMPARE_COMBINE	\
					"pmovmskb %%xmm7, %0 \n\t"		/* if one of these primes divides this location, this will be !0*/ \
					: "=r"(result) \
					: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections) \
					: "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "cc", "memory" \
					);
		#endif

	#elif defined(HAS_MMX)
		#define SIEVE_SCAN_32		\
			asm volatile (			\
				"movq (%1), %%mm0     \n\t"		\
				"por 8(%1), %%mm0     \n\t"		\
				"por 16(%1), %%mm0    \n\t"		\
				"por 24(%1), %%mm0    \n\t"		\
				"pmovmskb %%mm0, %0   \n\t"		\
				: "=r"(result)					\
				: "r"(sieveblock + j), "0"(result)	\
				: "%mm0");	

		#define SIEVE_SCAN_64		\
			asm volatile (			\
				"movq (%1), %%mm0     \n\t"		\
				"por 8(%1), %%mm0     \n\t"		\
				"por 16(%1), %%mm0    \n\t"		\
				"por 24(%1), %%mm0    \n\t"		\
				"por 32(%1), %%mm0    \n\t"		\
				"por 40(%1), %%mm0    \n\t"		\
				"por 48(%1), %%mm0    \n\t"		\
				"por 56(%1), %%mm0    \n\t"		\
				"pmovmskb %%mm0, %0   \n\t"		\
				: "=r"(result)					\
				: "r"(sieveblock + j), "0"(result)	\
				: "%mm0");

		#define SIEVE_SCAN_128		\
			asm volatile (				\
				"movq (%1), %%mm0     \n\t"		\
				"por 8(%1), %%mm0     \n\t"		\
				"por 16(%1), %%mm0    \n\t"		\
				"por 24(%1), %%mm0    \n\t"		\
				"por 32(%1), %%mm0    \n\t"		\
				"por 40(%1), %%mm0    \n\t"		\
				"por 48(%1), %%mm0    \n\t"		\
				"por 56(%1), %%mm0    \n\t"		\
				"por 64(%1), %%mm0     \n\t"	\
				"por 72(%1), %%mm0    \n\t"		\
				"por 80(%1), %%mm0    \n\t"		\
				"por 88(%1), %%mm0    \n\t"		\
				"por 96(%1), %%mm0    \n\t"		\
				"por 104(%1), %%mm0    \n\t"	\
				"por 112(%1), %%mm0    \n\t"	\
				"por 120(%1), %%mm0		\n\t"	\
				"pmovmskb %%mm0, %0   \n\t"		\
				: "=r"(result)					\
				: "r"(sieveblock + j), "0"(result)	\
				: "%mm0");

		#define SCAN_16X			\
			asm volatile (					/*this hasn't been tested yet...*/	\
				"movq (%2), %%mm0	\n\t"	/*move mask into xmm0*/	\
				"movq (%1), %%mm1	\n\t"	/*move 16 bptr locations into xmm regs*/	\
				"movq 8(%1), %%mm2	\n\t"		\
				"movq 16(%1), %%mm3	\n\t"		\
				"movq 24(%1), %%mm4	\n\t"		\
				"pcmpeqw %%mm0, %%mm1	\n\t"	/*compare to mask*/	\
				"pcmpeqw %%mm0, %%mm2	\n\t"		\
				"pcmpeqw %%mm0, %%mm3	\n\t"		\
				"pcmpeqw %%mm0, %%mm4	\n\t"		\
				"por %%mm1, %%mm4		\n\t"	/*or the comparisons*/	\
				"por %%mm2, %%mm3		\n\t"		\
				"por %%mm3, %%mm4		\n\t"		\
				"pmovmskb %%mm4, %0	\n\t"	/*if any are equal, this will be !0*/	\
				: "=r"(result)						\
				: "r"(bptr + j), "r"(mask)			\
				: "%mm0", "%mm1", "%mm2", "%mm3", "%mm4");		\
			asm volatile (			\
				"movl %0, %%ebx	\n\t"	/*remember result of first 8 comparisons*/	\
				"movq 8(%2), %%mm0	\n\t"	/*move mask into xmm0*/	\
				"movq 32(%1), %%mm1	\n\t"	/*move 16 bptr locations into xmm regs*/	\
				"movq 40(%1), %%mm2	\n\t"		\
				"movq 48(%1), %%mm3	\n\t"		\
				"movq 56(%1), %%mm4	\n\t"		\
				"pcmpeqw %%mm0, %%mm1	\n\t"	/*compare to mask*/	\
				"pcmpeqw %%mm0, %%mm2	\n\t"		\
				"pcmpeqw %%mm0, %%mm3	\n\t"		\
				"pcmpeqw %%mm0, %%mm4	\n\t"		\
				"por %%mm1, %%mm4		\n\t"	/*or the comparisons*/	\
				"por %%mm2, %%mm3		\n\t"		\
				"por %%mm3, %%mm4		\n\t"		\
				"pmovmskb %%mm4, %0	\n\t"	/*if any are equal, this will be !0*/	\
				"orl %%ebx, %0			\n\t"	/*combine with these 8 comparisons*/	\
				: "+r"(result)						\
				: "r"(bptr + j), "r"(mask)			\
				: "%mm0", "%mm1", "%mm2", "%mm3", "%mm4", "%ebx", "cc");	

	#else
		#define SCAN_16X	\
			result = 1;	/*dont know what compiler this is. force the normal method*/
		#undef SIMD_SIEVE_SCAN
	#endif

#elif defined(MSC_ASM32A)
	#define SIMD_SIEVE_SCAN 1
	#define SCAN_CLEAN ASM_M {emms};
	#define SSE2_RESIEVING 1

	#if defined(HAS_SSE2)
		//top level sieve scanning with SSE2
		#define SIEVE_SCAN_32	\
			do	{						\
				uint64 *localblock = sieveblock + j;	\
				ASM_M  {			\
					ASM_M mov edi, localblock			\
					ASM_M movdqa xmm0, XMMWORD PTR [edi]	\
					ASM_M por xmm0, XMMWORD PTR [edi + 16]	\
					ASM_M pmovmskb ecx, xmm0			\
					ASM_M mov result, ecx};			\
			} while (0);


		#define SIEVE_SCAN_64	\
			do	{						\
				uint64 *localblock = sieveblock + j;	\
				ASM_M  {			\
					ASM_M mov edi, localblock			\
					ASM_M movdqa xmm0, XMMWORD PTR [edi]	\
					ASM_M por xmm0, XMMWORD PTR [edi + 16]	\
					ASM_M por xmm0, XMMWORD PTR [edi + 32]	\
					ASM_M por xmm0, XMMWORD PTR [edi + 48]	\
					ASM_M pmovmskb ecx, xmm0			\
					ASM_M mov result, ecx};			\
			} while (0);

		#define SIEVE_SCAN_128	\
			do	{						\
				uint64 *localblock = sieveblock + j;	\
				ASM_M  {			\
					ASM_M mov edi, localblock			\
					ASM_M movdqa xmm0, XMMWORD PTR [edi]	\
					ASM_M por xmm0, XMMWORD PTR [edi + 16]	\
					ASM_M por xmm0, XMMWORD PTR [edi + 32]	\
					ASM_M por xmm0, XMMWORD PTR [edi + 48]	\
					ASM_M por xmm0, XMMWORD PTR [edi + 64]	\
					ASM_M por xmm0, XMMWORD PTR [edi + 80]	\
					ASM_M por xmm0, XMMWORD PTR [edi + 96]	\
					ASM_M por xmm0, XMMWORD PTR [edi + 112]	\
					ASM_M pmovmskb ecx, xmm0			\
					ASM_M mov result, ecx};			\
			} while (0);

		#define SCAN_16X			\
		do {							\
			/*bucket_element *local_bptr = bptr + j;*/	\
			uint32 *local_bptr = bptr + j; \
			uint16 *local_mask = &mask[0];	\
			ASM_M  {			\
				ASM_M mov edi, local_mask		\
				ASM_M movdqa xmm0, XMMWORD PTR [edi]		\
				ASM_M mov edi, local_bptr	\
				ASM_M movdqa xmm1, XMMWORD PTR [edi]	\
				ASM_M movdqa xmm2, XMMWORD PTR [edi + 16]	\
				ASM_M movdqa xmm3, XMMWORD PTR [edi + 32]	\
				ASM_M movdqa xmm4, XMMWORD PTR [edi + 48]	\
				ASM_M pcmpeqw xmm1, xmm0	\
				ASM_M pcmpeqw xmm2, xmm0	\
				ASM_M pcmpeqw xmm3, xmm0	\
				ASM_M pcmpeqw xmm4, xmm0	\
				ASM_M por xmm4, xmm1		\
				ASM_M por xmm2, xmm3		\
				ASM_M por xmm4, xmm2		\
				ASM_M pmovmskb eax, xmm4	\
				ASM_M mov result, eax		}	\
			} while (0);

		#define STEP_COMPARE_COMBINE \
			ASM_M psubw xmm2, xmm1 \
			ASM_M psubw xmm3, xmm1 \
			ASM_M pcmpeqw xmm5, xmm2 \
			ASM_M pcmpeqw xmm6, xmm3 \
			ASM_M por xmm7, xmm5 \
			ASM_M por xmm7, xmm6


		#define INIT_RESIEVE \
			ASM_M movdqa xmm4, XMMWORD PTR [edx] \
			ASM_M movdqa xmm2, XMMWORD PTR [ebx] \
			ASM_M paddw xmm2, xmm4 \
			ASM_M movdqa xmm3, XMMWORD PTR [ecx] \
			ASM_M paddw xmm3, xmm4 \
			ASM_M movdqa xmm1, XMMWORD PTR [eax] \
			ASM_M pxor xmm7, xmm7 \
			ASM_M pxor xmm6, xmm6 \
			ASM_M pxor xmm5, xmm5

		#ifdef YAFU_64K

			#define RESIEVE_4X_14BIT_MAX \
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
					ASM_M pmovmskb eax, xmm7	\
					ASM_M mov result, eax } \
				} while (0);

			#define RESIEVE_4X_15BIT_MAX \
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
					ASM_M pmovmskb eax, xmm7	\
					ASM_M mov result, eax } \
				} while (0);

			#define RESIEVE_4X_16BIT_MAX \
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
					ASM_M pmovmskb eax, xmm7	\
					ASM_M mov result, eax } \
				} while (0);

		#else

			#define RESIEVE_4X_14BIT_MAX \
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
					ASM_M pmovmskb eax, xmm7	\
					ASM_M mov result, eax } \
				} while (0);

			#define RESIEVE_4X_15BIT_MAX \
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
					ASM_M pmovmskb eax, xmm7	\
					ASM_M mov result, eax } \
				} while (0);

			#define RESIEVE_4X_16BIT_MAX \
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
					ASM_M pmovmskb eax, xmm7	\
					ASM_M mov result, eax } \
				} while (0);

		#endif

	#elif defined(HAS_MMX)

		#define SIEVE_SCAN_32	\
			do	{						\
				uint64 *localblock = sieveblock + j;	\
				ASM_M  {			\
					ASM_M mov edi, localblock			\
					ASM_M movq mm0, QWORD PTR [edi]	\
					ASM_M por mm0, QWORD PTR [edi + 8]	\
					ASM_M por mm0, QWORD PTR [edi + 16]	\
					ASM_M por mm0, QWORD PTR [edi + 24]	\
					ASM_M pmovmskb ecx, mm0			\
					ASM_M mov result, ecx};			\
			} while (0);

		#define SIEVE_SCAN_64	\
			do	{						\
				uint64 *localblock = sieveblock + j;	\
				ASM_M  {			\
					ASM_M mov edi, localblock			\
					ASM_M movq mm0, QWORD PTR [edi]	\
					ASM_M por mm0, QWORD PTR [edi + 8]	\
					ASM_M por mm0, QWORD PTR [edi + 16]	\
					ASM_M por mm0, QWORD PTR [edi + 24]	\
					ASM_M por mm0, QWORD PTR [edi + 32]	\
					ASM_M por mm0, QWORD PTR [edi + 40]	\
					ASM_M por mm0, QWORD PTR [edi + 48]	\
					ASM_M por mm0, QWORD PTR [edi + 56]	\
					ASM_M pmovmskb ecx, mm0			\
					ASM_M mov result, ecx};			\
			} while (0);

		#define SIEVE_SCAN_128	\
			do	{						\
				uint64 *localblock = sieveblock + j;	\
				ASM_M  {			\
					ASM_M mov edi, localblock			\
					ASM_M movq mm0, QWORD PTR [edi]	\
					ASM_M por mm0, QWORD PTR [edi + 8]	\
					ASM_M por mm0, QWORD PTR [edi + 16]	\
					ASM_M por mm0, QWORD PTR [edi + 24]	\
					ASM_M por mm0, QWORD PTR [edi + 32]	\
					ASM_M por mm0, QWORD PTR [edi + 40]	\
					ASM_M por mm0, QWORD PTR [edi + 48]	\
					ASM_M por mm0, QWORD PTR [edi + 56]	\
					ASM_M por mm0, QWORD PTR [edi + 64]	\
					ASM_M por mm0, QWORD PTR [edi + 72]	\
					ASM_M por mm0, QWORD PTR [edi + 80]	\
					ASM_M por mm0, QWORD PTR [edi + 88]	\
					ASM_M por mm0, QWORD PTR [edi + 96]	\
					ASM_M por mm0, QWORD PTR [edi + 104]	\
					ASM_M por mm0, QWORD PTR [edi + 112]	\
					ASM_M por mm0, QWORD PTR [edi + 120]	\
					ASM_M pmovmskb ecx, mm0			\
					ASM_M mov result, ecx};			\
			} while (0);


		#define SCAN_16X			\
		do {							\
			bucket_element *local_bptr = bptr + j;	\
			uint16 *local_mask = &mask[0];	\
			ASM_M {							\
			ASM_M mov edi, local_mask				\
			ASM_M movq mm0, QWORD PTR [edi]			\
			ASM_M mov edi, local_bptr				\
			ASM_M movq mm1, QWORD PTR [edi]			\
			ASM_M movq mm2, QWORD PTR [edi + 8]		\
			ASM_M movq mm3, QWORD PTR [edi + 16]	\
			ASM_M movq mm4, QWORD PTR [edi + 24]	\
			ASM_M pcmpeqw mm1, mm0	\
			ASM_M pcmpeqw mm2, mm0	\
			ASM_M pcmpeqw mm3, mm0	\
			ASM_M pcmpeqw mm4, mm0	\
			ASM_M por mm4, mm1		\
			ASM_M por mm2, mm3		\
			ASM_M por mm4, mm2		\
			ASM_M pmovmskb ecx, mm4	\
			ASM_M mov result, ecx	\
			}						\
		} while (0);				\
		do {							\
			bucket_element *local_bptr = bptr + j;	\
			uint16 *local_mask = &mask[4];			\
			ASM_M mov ebx, result					\
			ASM_M mov edi, local_mask				\
			ASM_M movq mm0, QWORD PTR [edi]			\
			ASM_M mov edi, local_bptr				\
			ASM_M movq mm1, QWORD PTR [edi + 32]	\
			ASM_M movq mm2, QWORD PTR [edi + 40]	\
			ASM_M movq mm3, QWORD PTR [edi + 48]	\
			ASM_M movq mm4, QWORD PTR [edi + 56]	\
			ASM_M pcmpeqw mm1, mm0		\
			ASM_M pcmpeqw mm2, mm0	\
			ASM_M pcmpeqw mm3, mm0	\
			ASM_M pcmpeqw mm4, mm0	\
			ASM_M por mm4, mm1		\
			ASM_M por mm2, mm3		\
			ASM_M por mm4, mm2		\
			ASM_M pmovmskb ecx, mm4	\
			ASM_M or ecx, ebx	\
			ASM_M mov result, ecx	\
		} while (0);

	#else
		#undef SIMD_SIEVE_SCAN
	#endif

#elif defined(_WIN64)

	#define SIMD_SIEVE_SCAN 1
	#define SSE2_RESIEVING 1

	#if defined(HAS_SSE2)
		//top level sieve scanning with SSE2
		#define SIEVE_SCAN_32	\
			do	{						\
				__m128i local_block;	\
				__m128i local_block2;	\
				local_block = _mm_load_si128(sieveblock + j); \
				local_block2 = _mm_load_si128(sieveblock + j + 2); \
				local_block = _mm_or_si128(local_block, local_block2); \
				result = _mm_movemask_epi8(local_block); \
			} while (0);


		#define SIEVE_SCAN_64	\
			do	{				  		\
				__m128i local_block;	\
				__m128i local_block2;	\
				__m128i local_block3;	\
				__m128i local_block4;	\
				local_block = _mm_load_si128(sieveblock + j); \
				local_block2 = _mm_load_si128(sieveblock + j + 2); \
				local_block3 = _mm_load_si128(sieveblock + j + 4); \
				local_block4 = _mm_load_si128(sieveblock + j + 6); \
				local_block = _mm_or_si128(local_block, local_block2); \
				local_block = _mm_or_si128(local_block, local_block3); \
				local_block = _mm_or_si128(local_block, local_block4); \
				result = _mm_movemask_epi8(local_block); \
			} while (0);

		#define SIEVE_SCAN_128	\
			do	{						\
				__m128i local_block;	\
				__m128i local_block2;	\
				__m128i local_block3;	\
				__m128i local_block4;	\
				__m128i local_block5;	\
				__m128i local_block6;	\
				__m128i local_block7;	\
				__m128i local_block8;	\
				local_block = _mm_load_si128(sieveblock + j); \
				local_block2 = _mm_load_si128(sieveblock + j + 2); \
				local_block3 = _mm_load_si128(sieveblock + j + 4); \
				local_block4 = _mm_load_si128(sieveblock + j + 6); \
				local_block5 = _mm_load_si128(sieveblock + j + 8); \
				local_block6 = _mm_load_si128(sieveblock + j + 10); \
				local_block7 = _mm_load_si128(sieveblock + j + 12); \
				local_block8 = _mm_load_si128(sieveblock + j + 14); \
				local_block = _mm_or_si128(local_block, local_block2); \
				local_block = _mm_or_si128(local_block, local_block3); \
				local_block = _mm_or_si128(local_block, local_block4); \
				local_block = _mm_or_si128(local_block, local_block5); \
				local_block = _mm_or_si128(local_block, local_block6); \
				local_block = _mm_or_si128(local_block, local_block7); \
				local_block = _mm_or_si128(local_block, local_block8); \
				result = _mm_movemask_epi8(local_block); \
			} while (0);

		#define SCAN_16X			\
		do {							\
			__m128i local_mask; \
			__m128i local_bptr; \
			__m128i local_bptr2; \
			__m128i local_bptr3; \
			__m128i local_bptr4; \
			local_mask = _mm_load_si128(&mask[0]); \
			local_bptr = _mm_load_si128(bptr + j); \
			local_bptr2 = _mm_load_si128(bptr + j + 4); \
			local_bptr3 = _mm_load_si128(bptr + j + 8); \
			local_bptr4 = _mm_load_si128(bptr + j + 12); \
			local_bptr = _mm_cmpeq_epi16(local_bptr, local_mask); \
			local_bptr2 = _mm_cmpeq_epi16(local_bptr2, local_mask); \
			local_bptr3 = _mm_cmpeq_epi16(local_bptr3, local_mask); \
			local_bptr4 = _mm_cmpeq_epi16(local_bptr4, local_mask); \
			local_bptr4 = _mm_or_si128(local_bptr4, local_bptr); \
			local_bptr2 = _mm_or_si128(local_bptr2, local_bptr3); \
			local_bptr4 = _mm_or_si128(local_bptr4, local_bptr2); \
			result = _mm_movemask_epi8(local_bptr4); \
			} while (0);

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

			#define RESIEVE_4X_14BIT_MAX \
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

			#define RESIEVE_4X_15BIT_MAX \
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

			#define RESIEVE_4X_16BIT_MAX \
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

			#define RESIEVE_4X_14BIT_MAX \
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

			#define RESIEVE_4X_15BIT_MAX \
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

			#define RESIEVE_4X_16BIT_MAX \
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

	#define SCAN_CLEAN /*nothing*/

#else	/* compiler not recognized*/
	#define SCAN_16X	\
			result = 1;	/*dont know what compiler this is. force the normal method*/
	#define SCAN_CLEAN /*nothing*/

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
		#define RESIEVE_4X_14BIT_MAX \
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

		#define RESIEVE_4X_15BIT_MAX \
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

		#define RESIEVE_4X_16BIT_MAX \
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
		#define RESIEVE_4X_14BIT_MAX \
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

		#define RESIEVE_4X_15BIT_MAX \
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

		#define RESIEVE_4X_16BIT_MAX \
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

#define DIVIDE_ONE_PRIME \
	do \
	{						\
		fb_offsets[++smooth_num] = i;	\
		zShortDiv32(Q,prime,Q);			\
	} while (zShortMod32(Q,prime) == 0);

#define PROTECTED_DIVIDE_ONE_PRIME \
	while (zShortMod32(Q,prime) == 0) \
	{						\
		fb_offsets[++smooth_num] = i;	\
		zShortDiv32(Q,prime,Q);			\
	}

//#define SCAN_MASK (( (uint64)0x80808080 << 32) | (uint64)0x80808080)
#define SCAN_MASK 0x8080808080808080


	//when we compress small primes into 16 bits of a 32 bit field, the
	//trick of fooling the sieve routine to not sieve those roots which
	//divide poly_a fails when the blocksize is 2^16, because we're doing this:
	//root1 = fbptr->roots & 0xFFFF;
	//root2 = fbptr->roots >> 16;
	//set_aprime_roots sets roots to all 1's, which then results in root1 and 
	//root2 being set to 65535 in the sieve routine.  this, of course, isn't right
	//so the sieve location 65535 is corrupted by many small prime hits when it
	//shouldn't be, and thus we might end up here more often then we should for
	//offset 65535.
	//
	//even if we do end up here when we shouldn't, often we'll fail to find many
	//small primes which actually divide this location, and we'll bail anyway.  this
	//is safe because we explicitly trial divide by these small primes.  
	//if we make it past the small prime test and go to check the progression
	//of a prime which divides poly_a then the roots we arrive at are false (65535 again)
	//but our computation of the progression will always be 65535 + prime - blocksize,
	//since we set the root to 65535 during the sieve step as well.  
	//65535 != 65535 + prime - blocksize, so we are safe here as well.
	//we may incur more trial division than necessary - is that better than always
	//throwing away block location 65535 - NO, empirically it is much better to just
	//always bail for location 65535 when the blocksize is 65536.

int check_relations_siqs_1(uint32 blocknum, uint8 parity, 
						   static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//not unrolled; for small inputs

	uint32 j,k,it=BLOCKSIZE>>3;
	uint32 thisloc;
	uint64 *sieveblock;
	uint64 mask = SCAN_MASK;

	sieveblock = (uint64 *)dconf->sieve;
	dconf->num_reports = 0;

	//check for relations
	for (j=0;j<it;j++)
	{
		//check 8 locations simultaneously
		if ((sieveblock[j] & mask) == (uint64)(0))
			continue;

		//at least one passed the check, find which one(s) and pass to 
		//trial division stage
		for (k=0;k<8;k++)
		{
			thisloc = (j<<3) + k;
			if ((dconf->sieve[thisloc] & 0x80) == 0)			
				continue;

#ifdef YAFU_64K
			//see discussion near line 377
			if (thisloc ==	65535)
				continue;
#endif
			// log this report
			dconf->reports[dconf->num_reports++] = thisloc;			
		}
	}

	if (dconf->num_reports >= MAX_SIEVE_REPORTS)
	{
		printf("error: too many sieve reports (found %d)\n",dconf->num_reports);
		exit(-1);
	}

	//remove small primes, and test if its worth continuing for each report
	filter_SPV(parity, dconf->sieve, dconf->numB-1,blocknum,sconf,dconf);
	filter_medprimes(parity, dconf->numB-1,blocknum,sconf,dconf);

	// factor all reports in this block
	for (j=0; j<dconf->num_reports; j++)
	{
		if (dconf->valid_Qs[j])
		{
			filter_LP(j, parity, blocknum, sconf, dconf);
			trial_divide_Q_siqs(j, parity, dconf->numB-1, blocknum,sconf,dconf);
		}
	}

	return 0;
}

int check_relations_siqs_4(uint32 blocknum, uint8 parity, 
						   static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//unrolled x32; for medium inputs

	uint32 i,j,k,it=BLOCKSIZE>>3;
	uint32 thisloc;
	uint64 *sieveblock;

	sieveblock = (uint64 *)dconf->sieve;
	dconf->num_reports = 0;

	//check for relations
	for (j=0;j<it;j+=4)	
	{

#ifdef SIMD_SIEVE_SCAN

		uint32 result = 0;

		SIEVE_SCAN_32;

		if (result == 0)
			continue;

#else
		uint64 mask = SCAN_MASK;

		if (((sieveblock[j] | sieveblock[j+1] | 
			sieveblock[j+2] | sieveblock[j+3]) & 
		      mask) == (uint64)(0))
			continue;
#endif
	
#if defined(SIMD_SIEVE_SCAN)
		// make it safe to perform floating point
		SCAN_CLEAN;
#endif

		//at least one passed the check, find which one(s) and pass to 
		//trial division stage
		for (i=0; i<4; i++)
		{
			//check 8 locations simultaneously
			if ((sieveblock[j + i] & 0x8080808080808080) == (uint64)(0))
				continue;

			//at least one passed the check, find which one(s) and pass to 
			//trial division stage
			for (k=0;k<8;k++)
			{
				thisloc = ((j+i)<<3) + k;
				if ((dconf->sieve[thisloc] & 0x80) == 0)
					continue;

#ifdef YAFU_64K
				//see discussion near line 377
				if (thisloc == 65535)
					continue;
#endif
				// log this report
				dconf->reports[dconf->num_reports++] = thisloc;
			}
		}
	}

#if defined(SIMD_SIEVE_SCAN)
	// make it safe to perform floating point
	SCAN_CLEAN;
#endif

	if (dconf->num_reports >= MAX_SIEVE_REPORTS)
	{
		printf("error: too many sieve reports (found %d)\n",dconf->num_reports);
		exit(-1);
	}

	//remove small primes, and test if its worth continuing for each report
	filter_SPV(parity, dconf->sieve,dconf->numB-1,blocknum,sconf,dconf);
	filter_medprimes(parity, dconf->numB-1,blocknum,sconf,dconf);

	// factor all reports in this block
	for (j=0; j<dconf->num_reports; j++)
	{
		if (dconf->valid_Qs[j])
		{
			filter_LP(j, parity, blocknum, sconf, dconf);
			trial_divide_Q_siqs(j, parity, dconf->numB-1, blocknum,sconf,dconf);
		}
	}

	return 0;
}

int check_relations_siqs_8(uint32 blocknum, uint8 parity, 
						   static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//unrolled x64; for large inputs
	uint32 i,j,k,it=BLOCKSIZE>>3;
	uint32 thisloc;
	uint64 *sieveblock;
	sieve_fb_compressed *fbptr, *fbc;
	int prime, root1, root2;

	sieveblock = (uint64 *)dconf->sieve;
	dconf->num_reports = 0;

	//check for relations
	for (j=0;j<it;j+=8)
	{

#ifdef SIMD_SIEVE_SCAN

		uint32 result = 0;

		SIEVE_SCAN_64;

		if (result == 0)
			continue;

#else
		uint64 mask = SCAN_MASK;

		if (((sieveblock[j] | sieveblock[j+1] | sieveblock[j+2] | sieveblock[j+3] |
		      sieveblock[j+4] | sieveblock[j+5] | sieveblock[j+6] | sieveblock[j+7]) 
			  & mask) == (uint64)(0))
			continue;
#endif
		
#if defined(SIMD_SIEVE_SCAN)
		// make it safe to perform floating point
		SCAN_CLEAN;
#endif

		//at least one passed the check, find which one(s) and pass to 
		//trial division stage
		for (i=0; i<8; i++)
		{
			//check 8 locations simultaneously
			if ((sieveblock[j + i] & 0x8080808080808080) == (uint64)(0))
				continue;

			for (k=0;k<8;k++)
			{
				thisloc = ((j+i)<<3) + k;
				if ((dconf->sieve[thisloc] & 0x80) == 0)
					continue;

#ifdef YAFU_64K
				//see discussion near line 377
				if (thisloc == 65535)
					continue;
#endif

				// log this report
				dconf->reports[dconf->num_reports++] = thisloc;
			}
		}
	}

#if defined(SIMD_SIEVE_SCAN)
	// make it safe to perform floating point
	SCAN_CLEAN;
#endif

	if (dconf->num_reports >= MAX_SIEVE_REPORTS)
	{
		printf("error: too many sieve reports (found %d)\n",dconf->num_reports);
		exit(-1);
	}

	//printf("block %d found %d reports\n", blocknum, dconf->num_reports);

	//remove small primes, and test if its worth continuing for each report
	filter_SPV(parity, dconf->sieve, dconf->numB-1, blocknum,sconf,dconf);
	filter_medprimes(parity, dconf->numB-1,blocknum,sconf,dconf);

	// factor all reports in this block
	for (j=0; j<dconf->num_reports; j++)
	{
		if (dconf->valid_Qs[j])
		{
			filter_LP(j, parity, blocknum, sconf, dconf);
			trial_divide_Q_siqs(j, parity, dconf->numB-1, blocknum,sconf,dconf);
		}
	}

	return 0;
}


int check_relations_siqs_16(uint32 blocknum, uint8 parity, 
						   static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//unrolled x128; for large inputs
	uint32 i,j,k,it=BLOCKSIZE>>3;
	uint32 thisloc;
	uint64 *sieveblock;

	sieveblock = (uint64 *)dconf->sieve;
	dconf->num_reports = 0;

	//check for relations
	for (j=0;j<it;j+=16)
	{

#ifdef SIMD_SIEVE_SCAN

		uint32 result = 0;

		SIEVE_SCAN_128;

		if (result == 0)
			continue;

#else
		uint64 mask = SCAN_MASK;

		if (((sieveblock[j] | sieveblock[j+1] | sieveblock[j+2] | sieveblock[j+3] |
		      sieveblock[j+4] | sieveblock[j+5] | sieveblock[j+6] | sieveblock[j+7] |
			  sieveblock[j+8] | sieveblock[j+9] | sieveblock[j+10] | sieveblock[j+11] |
		      sieveblock[j+12] | sieveblock[j+13] | sieveblock[j+14] | sieveblock[j+15]
				) & mask) == (uint64)(0))
			continue;
#endif
		

#if defined(SIMD_SIEVE_SCAN)
		// make it safe to perform floating point
		SCAN_CLEAN;
#endif

		//at least one passed the check, find which one(s) and pass to 
		//trial division stage
		for (i=0; i<16; i++)
		{
			//check 8 locations simultaneously
			if ((sieveblock[j + i] & 0x8080808080808080) == (uint64)(0))
				continue;

			for (k=0;k<8;k++)
			{
				thisloc = ((j+i)<<3) + k;
				if ((dconf->sieve[thisloc] & 0x80) == 0)
					continue;

#ifdef YAFU_64K
				//see discussion near line 377
				if (thisloc == 65535)
					continue;
#endif
				// log this report
				dconf->reports[dconf->num_reports++] = thisloc;
			}
		}
	}

#if defined(SIMD_SIEVE_SCAN)
	// make it safe to perform floating point
	SCAN_CLEAN;
#endif

	if (dconf->num_reports >= MAX_SIEVE_REPORTS)
	{
		printf("error: too many sieve reports (found %d)\n",dconf->num_reports);
		exit(-1);
	}

	//remove small primes, and test if its worth continuing for each report
	filter_SPV(parity, dconf->sieve, dconf->numB-1,blocknum,sconf,dconf);
	filter_medprimes(parity, dconf->numB-1,blocknum,sconf,dconf);

	// factor all reports in this block
	for (j=0; j<dconf->num_reports; j++)
	{
		if (dconf->valid_Qs[j])
		{
			filter_LP(j, parity, blocknum, sconf, dconf);
			trial_divide_Q_siqs(j, parity, dconf->numB-1, blocknum,sconf,dconf);
		}
	}

	return 0;
}

void filter_SPV(uint8 parity, uint8 *sieve, uint32 poly_id, uint32 bnum, 
				static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//we have flagged this sieve offset as likely to produce a relation
	//nothing left to do now but check and see.
	int i;
	uint32 bound, tmp, prime, root1, root2;
	int smooth_num;
	sieve_fb *fb;
	sieve_fb_compressed *fbptr, *fbc;
	fb_element_siqs *fullfb_ptr, *fullfb = sconf->factor_base->list;
	uint8 logp, bits;
	uint32 tmp1, tmp2, tmp3, tmp4, offset, report_num;
	z32 *Q;

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
		//this one qualifies to check further, log that fact.
		dconf->num++;

		smooth_num = -1;

		//this one is close enough, compute 
		//Q(x)/a = (ax + b)^2 - N, where x is the sieve index
		//Q(x)/a = (ax + 2b)x + c;	
		offset = (bnum << BLOCKBITS) + dconf->reports[report_num];

		//multiple precision arithmetic.  all the qstmp's are a global hack
		//but I don't want to Init/Free millions of times if I don't have to.
		zShiftLeft(&dconf->qstmp2,&dconf->curr_poly->poly_b,1);
		zShortMul(&dconf->curr_poly->poly_a,offset,&dconf->qstmp1);
		if (parity)
			zSub(&dconf->qstmp1,&dconf->qstmp2,&dconf->qstmp3);
		else
			zAdd(&dconf->qstmp1,&dconf->qstmp2,&dconf->qstmp3);

		zShortMul(&dconf->qstmp3,offset,&dconf->qstmp1);
		zAdd(&dconf->qstmp1,&dconf->curr_poly->poly_c,&dconf->qstmp4);
	
		//this is another hack because on most fast systems the multiple
		//precision arith is 64 bit based, but it turns out that the MP mod's
		//we have to do a lot of in trial division I've implemented faster
		//in 32 bit base.  The actual conversion here is just a cast.
	
#if BITS_PER_DIGIT == 32
		for (i=0; i<abs(dconf->qstmp4.size); i++)
			dconf->Qvals[report_num].val[i] = (uint32)dconf->qstmp4.val[i];
		dconf->Qvals[report_num].size = dconf->qstmp4.size;
		dconf->Qvals[report_num].type = dconf->qstmp4.type;

#else
		z64_to_z32(&dconf->qstmp4,&dconf->Qvals[report_num]);
#endif

		Q = &dconf->Qvals[report_num];

		//we have two signs to worry about.  the sign of the offset tells us how to calculate ax + b, while
		//the sign of Q(x) tells us how to factor Q(x) (with or without a factor of -1)
		//the square root phase will need to know both.  fboffset holds the sign of Q(x).  the sign of the 
		//offset is stored standalone in the relation structure.
		if (Q->size < 0)
		{
			dconf->fb_offsets[report_num][++smooth_num] = 0;
			Q->size = Q->size * -1;
		}
	
		//compute the bound for small primes.  if we can't find enough small
		//primes, then abort the trial division early because it is likely to fail to
		//produce even a partial relation.
		bits = sieve[dconf->reports[report_num]];
		bits = (255 - bits) + sconf->tf_closnuf + 1;

		//take care of powers of two
		while (!(Q->val[0] & 1))
		{
			zShiftRight32(Q,Q,1);
			dconf->fb_offsets[report_num][++smooth_num] = 1;
			bits++;
		}

		i=2;
		//explicitly trial divide by small primes which we have not
		//been sieving.  because we haven't been sieving, their progressions
		//have not been updated and thus we can't use the faster methods
		//seen below.  fortunately, there shouldn't be many of these to test
		//to speed things up, use multiplication by inverse rather than 
		//division, and do things in batches of 4 so we can use
		//the whole cache line at once (16 byte structure)

		//do the small primes in optimized batches of 4
		bound = (sconf->sieve_small_fb_start - 4);
		
		while ((uint32)i < bound)
		{
			uint64 q64;

			tmp1 = offset + fullfb_ptr->correction[i];
			q64 = (uint64)tmp1 * (uint64)fullfb_ptr->small_inv[i];
			tmp1 = q64 >> 32; 
			//at this point tmp1 is offset / prime
			tmp1 = offset - tmp1 * fullfb_ptr->prime[i];
			i++;

			tmp2 = offset + fullfb_ptr->correction[i];
			q64 = (uint64)tmp2 * (uint64)fullfb_ptr->small_inv[i];
			tmp2 = q64 >> 32; 
			tmp2 = offset - tmp2 * fullfb_ptr->prime[i];
			i++;

			tmp3 = offset + fullfb_ptr->correction[i];
			q64 = (uint64)tmp3 * (uint64)fullfb_ptr->small_inv[i];
			tmp3 = q64 >> 32;
			tmp3 = offset - tmp3 * fullfb_ptr->prime[i];
			i++;

			tmp4 = offset + fullfb_ptr->correction[i];
			q64 = (uint64)tmp4 * (uint64)fullfb_ptr->small_inv[i];
			tmp4 = q64 >> 32; 
			tmp4 = offset - tmp4 * fullfb_ptr->prime[i];
			
			i -= 3;

#ifdef USE_COMPRESSED_FB
			fbptr = fbc + i;
			prime = fbptr->prime_and_logp & 0xFFFF;
			root1 = fbptr->roots & 0xFFFF;
			root2 = fbptr->roots >> 16;
			logp = fbptr->prime_and_logp >> 16;
#else
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
			logp = fbc->logp[i];
#endif

			if ((tmp1 == root1 || tmp1 == root2) || 
				(root1 == prime && tmp1 == 0) || (root2 == prime && tmp1 == 0))
			{
				do
				{
					dconf->fb_offsets[report_num][++smooth_num] = i;
					zShortDiv32(Q,prime,Q);
					bits += logp;
				} while (zShortMod32(Q,prime) == 0);
			}

			i++;

#ifdef USE_COMPRESSED_FB
			fbptr = fbc + i;
			prime = fbptr->prime_and_logp & 0xFFFF;
			root1 = fbptr->roots & 0xFFFF;
			root2 = fbptr->roots >> 16;
			logp = fbptr->prime_and_logp >> 16;
#else
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
			logp = fbc->logp[i];
#endif

			if ((tmp2 == root1 || tmp2 == root2) || 
				(root1 == prime && tmp2 == 0) || (root2 == prime && tmp2 == 0))
			{
				do
				{
					dconf->fb_offsets[report_num][++smooth_num] = i;
					zShortDiv32(Q,prime,Q);
					bits += logp;
				} while (zShortMod32(Q,prime) == 0);
			}

			i++;

#ifdef USE_COMPRESSED_FB
			fbptr = fbc + i;
			prime = fbptr->prime_and_logp & 0xFFFF;
			root1 = fbptr->roots & 0xFFFF;
			root2 = fbptr->roots >> 16;
			logp = fbptr->prime_and_logp >> 16;
#else
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
			logp = fbc->logp[i];
#endif

			if ((tmp3 == root1 || tmp3 == root2) || 
				(root1 == prime && tmp3 == 0) || (root2 == prime && tmp3 == 0))
			{
				do
				{
					dconf->fb_offsets[report_num][++smooth_num] = i;
					zShortDiv32(Q,prime,Q);
					bits += logp;
				} while (zShortMod32(Q,prime) == 0);
			}

			i++;

#ifdef USE_COMPRESSED_FB
			fbptr = fbc + i;
			prime = fbptr->prime_and_logp & 0xFFFF;
			root1 = fbptr->roots & 0xFFFF;
			root2 = fbptr->roots >> 16;
			logp = fbptr->prime_and_logp >> 16;
#else
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
			logp = fbc->logp[i];
#endif

			if ((tmp4 == root1 || tmp4 == root2) || 
				(root1 == prime && tmp4 == 0) || (root2 == prime && tmp4 == 0))
			{
				do
				{
					dconf->fb_offsets[report_num][++smooth_num] = i;
					zShortDiv32(Q,prime,Q);
					bits += logp;
				} while (zShortMod32(Q,prime) == 0);
			}
			i++;
		}

		//finish up the rest of the small primes
		while ((uint32)i < sconf->sieve_small_fb_start)
		{
			uint64 q64;

#ifdef USE_COMPRESSED_FB
			fbptr = fbc + i;
			prime = fbptr->prime_and_logp & 0xFFFF;
			root1 = fbptr->roots & 0xFFFF;
			root2 = fbptr->roots >> 16;
			logp = fbptr->prime_and_logp >> 16;
#else
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
			logp = fbc->logp[i];
#endif
			
			//this is just offset % prime (but divisionless!)
			tmp = offset + fullfb_ptr->correction[i];
			q64 = (uint64)tmp * (uint64)fullfb_ptr->small_inv[i];
			tmp = q64 >>  32; 
			tmp = offset - tmp * prime;

			//if offset % prime == either root, it's on the progression.  also
			//need to check for the case if root1 or root2 == prime at the same
			//time as offset mod prime = 0.  for small primes, this happens fairly
			//often.  the simple offset % prime check will miss these cases.
			if (tmp == root1 || tmp == root2 || 
				(root1 == prime && tmp == 0) || (root2 == prime && tmp == 0))
			{
				do
				{
					dconf->fb_offsets[report_num][++smooth_num] = i;
					zShortDiv32(Q,prime,Q);
					bits += logp;
				} while (zShortMod32(Q,prime) == 0);
			}
			i++;
		}

		if (bits < (sconf->tf_closnuf + sconf->tf_small_cutoff))
			dconf->valid_Qs[report_num] = 0;
		else
			dconf->valid_Qs[report_num] = 1;

		dconf->smooth_num[report_num] = smooth_num;
	}

#ifdef QS_TIMING
	gettimeofday (&qs_timing_stop, NULL);
	qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

	TF_STG1 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
	free(qs_timing_diff);

	gettimeofday(&qs_timing_start, NULL);
#endif

	return;
}

void filter_LP(uint32 report_num,  uint8 parity, uint32 bnum, 
	static_conf_t *sconf, dynamic_conf_t *dconf)
{
	int i,j,k;
	uint32 basebucket, prime;
	int smooth_num;
	uint32 *fb_offsets;
	uint32 *bptr;
	sieve_fb *fb;
	uint32 block_loc;
	z32 *Q;
	uint16 *mask = dconf->mask;

#ifdef QS_TIMING
	gettimeofday(&qs_timing_start, NULL);
#endif

	fb_offsets = &dconf->fb_offsets[report_num][0];
	smooth_num = dconf->smooth_num[report_num];
	Q = &dconf->Qvals[report_num];
	block_loc = dconf->reports[report_num];

	mask[0] = block_loc;
	mask[2] = block_loc;
	mask[4] = block_loc;
	mask[6] = block_loc;
	
	if (parity)
		fb = dconf->fb_sieve_n;
	else
		fb = dconf->fb_sieve_p;

	//primes bigger than med_B are bucket sieved, so we need
	//only search through the bucket and see if any locations match the
	//current block index.
	bptr = dconf->buckets->list + (bnum << BUCKET_BITS);
	if (parity)
	{
		bptr += (sconf->num_blocks << BUCKET_BITS);
		basebucket = sconf->num_blocks;
	}
	else
		basebucket = 0;

	//times_checked++;
	for (k=0; (uint32)k < dconf->buckets->num_slices; k++)
	{
		uint32 lpnum = *(dconf->buckets->num + bnum + basebucket);

		uint32 fb_bound = *(dconf->buckets->fb_bounds + k);
		uint32 result;

		for (j=0; (uint32)j < (lpnum & (uint32)(~15)); j += 16)
		{
			SCAN_16X;

			if (result == 0)
				continue;

			//noticably faster to not put these in a loop!
			if ((bptr[j] & 0x0000ffff) == block_loc)
			{
				i = fb_bound + (bptr[j] >> 16);
				prime = fb[i].prime;
				DIVIDE_ONE_PRIME;
			}
			if ((bptr[j+1] & 0x0000ffff) == block_loc)
			{
				i = fb_bound + (bptr[j+1] >> 16);
				prime = fb[i].prime;
				DIVIDE_ONE_PRIME;
			}
			if ((bptr[j+2] & 0x0000ffff) == block_loc)
			{
				i = fb_bound + (bptr[j+2] >> 16);
				prime = fb[i].prime;
				DIVIDE_ONE_PRIME;
			}
			if ((bptr[j+3] & 0x0000ffff) == block_loc)
			{
				i = fb_bound + (bptr[j+3] >> 16);
				prime = fb[i].prime;
				DIVIDE_ONE_PRIME;
			}
			if ((bptr[j+4] & 0x0000ffff) == block_loc)
			{
				i = fb_bound + (bptr[j+4] >> 16);
				prime = fb[i].prime;
				DIVIDE_ONE_PRIME;
			}
			if ((bptr[j+5] & 0x0000ffff) == block_loc)
			{
				i = fb_bound + (bptr[j+5] >> 16);
				prime = fb[i].prime;
				DIVIDE_ONE_PRIME;
			}
			if ((bptr[j+6] & 0x0000ffff) == block_loc)
			{
				i = fb_bound + (bptr[j+6] >> 16);
				prime = fb[i].prime;
				DIVIDE_ONE_PRIME;
			}
			if ((bptr[j+7] & 0x0000ffff) == block_loc)
			{
				i = fb_bound + (bptr[j+7] >> 16);
				prime = fb[i].prime;
				DIVIDE_ONE_PRIME;
			}
			if ((bptr[j+8] & 0x0000ffff) == block_loc)
			{
				i = fb_bound + (bptr[j+8] >> 16);
				prime = fb[i].prime;
				DIVIDE_ONE_PRIME;
			}
			if ((bptr[j+9] & 0x0000ffff) == block_loc)
			{
				i = fb_bound + (bptr[j+9] >> 16);
				prime = fb[i].prime;
				DIVIDE_ONE_PRIME;
			}
			if ((bptr[j+10] & 0x0000ffff) == block_loc)
			{
				i = fb_bound + (bptr[j+10] >> 16);
				prime = fb[i].prime;
				DIVIDE_ONE_PRIME;
			}
			if ((bptr[j+11] & 0x0000ffff) == block_loc)
			{
				i = fb_bound + (bptr[j+11] >> 16);
				prime = fb[i].prime;
				DIVIDE_ONE_PRIME;
			}
			if ((bptr[j+12] & 0x0000ffff) == block_loc)
			{
				i = fb_bound + (bptr[j+12] >> 16);
				prime = fb[i].prime;
				DIVIDE_ONE_PRIME;
			}
			if ((bptr[j+13] & 0x0000ffff) == block_loc)
			{
				i = fb_bound + (bptr[j+13] >> 16);
				prime = fb[i].prime;
				DIVIDE_ONE_PRIME;
			}
			if ((bptr[j+14] & 0x0000ffff) == block_loc)
			{
				i = fb_bound + (bptr[j+14] >> 16);
				prime = fb[i].prime;
				DIVIDE_ONE_PRIME;
			}
			if ((bptr[j+15] & 0x0000ffff) == block_loc)
			{
				i = fb_bound + (bptr[j+15] >> 16);
				prime = fb[i].prime;
				DIVIDE_ONE_PRIME;
			}

			//if ((Q->size == 1) && (Q->val[0] < pmax ))
			//	goto early_abort;
		
		}
		
		for (; (uint32)j < lpnum; j++)
		{
			if ((bptr[j] & 0x0000ffff) == block_loc)
			{
				i = fb_bound + (bptr[j] >> 16);
				prime = fb[i].prime;
				//printf("block_loc = %u, bptr = %u, fb_bound = %u, fb_index = %u, prime = %u, Q mod prime = %u\n",
				//	block_loc, bptr[j].loc, fb_bound, bptr[j].fb_index, prime, zShortMod32(Q,prime));
				DIVIDE_ONE_PRIME;
			}
		}

		//point to the next slice of primes
		bptr += (sconf->num_blocks << (BUCKET_BITS + 1));
		basebucket += (sconf->num_blocks << 1);
	}

	SCAN_CLEAN;

#ifdef QS_TIMING
	gettimeofday (&qs_timing_stop, NULL);
	qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

	TF_STG5 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
	free(qs_timing_diff);
#endif

	dconf->smooth_num[report_num] = smooth_num;

	return;
}

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
	sieve_fb_compressed *fbptr, *fbc;
	fb_element_siqs *fullfb_ptr, *fullfb = sconf->factor_base->list;
	uint32 tmp1, tmp2, tmp3, tmp4, block_loc;
	uint16 *corrections;
	z32 *Q;	
	
#ifdef SSE2_RESIEVING
	#ifdef WIN32
		corrections = (uint16 *)_aligned_malloc(8 * sizeof(uint16),64);
	#else
		corrections = (uint16 *)memalign(64, 8 * sizeof(uint16));
	#endif
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
		if (!dconf->valid_Qs[report_num])
			continue;

		fb_offsets = &dconf->fb_offsets[report_num][0];
		smooth_num = dconf->smooth_num[report_num];
		Q = &dconf->Qvals[report_num];
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
#elif defined(YAFU_64K)
		bound = sconf->factor_base->fb_14bit_B;
#else
		bound = sconf->factor_base->fb_13bit_B;
#endif
		
		while ((uint32)i < bound && ((i & 3) != 0))
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
					* (uint64)fullfb_ptr->small_inv[i]) >> 40); 
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
			tmp1 = BLOCKSIZE - block_loc;
			tmp2 = BLOCKSIZE - block_loc;
			tmp3 = BLOCKSIZE - block_loc;
			tmp4 = BLOCKSIZE - block_loc;

			tmp1 = tmp1 + fullfb_ptr->correction[i];
			q64 = (uint64)tmp1 * (uint64)fullfb_ptr->small_inv[i];
			tmp1 = q64 >> 40; 
			tmp1 = tmp1 + 1;
			tmp1 = block_loc + tmp1 * fullfb_ptr->prime[i];
		
			i++;
		
			tmp2 = tmp2 + fullfb_ptr->correction[i];
			q64 = (uint64)tmp2 * (uint64)fullfb_ptr->small_inv[i];
			tmp2 = q64 >> 40; 
			tmp2 = tmp2 + 1;
			tmp2 = block_loc + tmp2 * fullfb_ptr->prime[i];

			i++;

			tmp3 = tmp3 + fullfb_ptr->correction[i];
			q64 = (uint64)tmp3 * (uint64)fullfb_ptr->small_inv[i];
			tmp3 = q64 >> 40; 
			tmp3 = tmp3 + 1;
			tmp3 = block_loc + tmp3 * fullfb_ptr->prime[i];

			i++;

			tmp4 = tmp4 + fullfb_ptr->correction[i];
			q64 = (uint64)tmp4 * (uint64)fullfb_ptr->small_inv[i];
			tmp4 = q64 >>  40; 
			tmp4 = tmp4 + 1;
			tmp4 = block_loc + tmp4 * fullfb_ptr->prime[i];

			tmp1 = tmp1 - BLOCKSIZE;
			tmp2 = tmp2 - BLOCKSIZE;
			tmp3 = tmp3 - BLOCKSIZE;
			tmp4 = tmp4 - BLOCKSIZE;

			i -= 3;

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

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		TF_STG2 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);

		gettimeofday(&qs_timing_start, NULL);
#endif

		//now cleanup any that don't fit in the last batch of 4
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
					* (uint64)fullfb_ptr->small_inv[i]) >> 40); 
			tmp = block_loc + tmp*prime;
			tmp = tmp - BLOCKSIZE;

			if (tmp == root1 || tmp == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;
		}

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
			RESIEVE_4X_14BIT_MAX;
			
			if (result == 0)
			{
				i += 8;
				continue;
			}

			if (result & 0x2)
			{
				while (zShortMod32(Q,fbc->prime[i]) == 0) 
				{						
					fb_offsets[++smooth_num] = i;	
					zShortDiv32(Q,fbc->prime[i],Q);			
				}
			}

			if (result & 0x8)
			{
				while (zShortMod32(Q,fbc->prime[i+1]) == 0) 
				{						
					fb_offsets[++smooth_num] = i+1;	
					zShortDiv32(Q,fbc->prime[i+1],Q);			
				}
			}

			if (result & 0x20)
			{
				while (zShortMod32(Q,fbc->prime[i+2]) == 0) 
				{						
					fb_offsets[++smooth_num] = i+2;	
					zShortDiv32(Q,fbc->prime[i+2],Q);			
				}
			}

			if (result & 0x80)
			{
				while (zShortMod32(Q,fbc->prime[i+3]) == 0) 
				{						
					fb_offsets[++smooth_num] = i+3;	
					zShortDiv32(Q,fbc->prime[i+3],Q);			
				}
			}

			if (result & 0x200)
			{
				while (zShortMod32(Q,fbc->prime[i+4]) == 0) 
				{						
					fb_offsets[++smooth_num] = i+4;	
					zShortDiv32(Q,fbc->prime[i+4],Q);			
				}
			}

			if (result & 0x800)
			{
				while (zShortMod32(Q,fbc->prime[i+5]) == 0) 
				{						
					fb_offsets[++smooth_num] = i+5;	
					zShortDiv32(Q,fbc->prime[i+5],Q);			
				}
			}

			if (result & 0x2000)
			{
				while (zShortMod32(Q,fbc->prime[i+6]) == 0) 
				{						
					fb_offsets[++smooth_num] = i+6;	
					zShortDiv32(Q,fbc->prime[i+6],Q);			
				}
			}

			if (result & 0x8000)
			{
				while (zShortMod32(Q,fbc->prime[i+7]) == 0) 
				{						
					fb_offsets[++smooth_num] = i+7;	
					zShortDiv32(Q,fbc->prime[i+7],Q);			
				}
			}

			i += 8;
		}
#endif
		bound = sconf->factor_base->fb_15bit_B;

		while ((uint32)i < bound)
		{
			uint32 result = 0;
			RESIEVE_4X_15BIT_MAX;

			if (result == 0)
			{
				i += 8;
				continue;
			}

			if (result & 0x2)
			{
				while (zShortMod32(Q,fbc->prime[i]) == 0) 
				{						
					fb_offsets[++smooth_num] = i;	
					zShortDiv32(Q,fbc->prime[i],Q);			
				}
			}

			if (result & 0x8)
			{
				while (zShortMod32(Q,fbc->prime[i+1]) == 0) 
				{						
					fb_offsets[++smooth_num] = i+1;	
					zShortDiv32(Q,fbc->prime[i+1],Q);			
				}
			}

			if (result & 0x20)
			{
				while (zShortMod32(Q,fbc->prime[i+2]) == 0) 
				{						
					fb_offsets[++smooth_num] = i+2;	
					zShortDiv32(Q,fbc->prime[i+2],Q);			
				}
			}

			if (result & 0x80)
			{
				while (zShortMod32(Q,fbc->prime[i+3]) == 0) 
				{						
					fb_offsets[++smooth_num] = i+3;	
					zShortDiv32(Q,fbc->prime[i+3],Q);			
				}
			}

			if (result & 0x200)
			{
				while (zShortMod32(Q,fbc->prime[i+4]) == 0) 
				{						
					fb_offsets[++smooth_num] = i+4;	
					zShortDiv32(Q,fbc->prime[i+4],Q);			
				}
			}

			if (result & 0x800)
			{
				while (zShortMod32(Q,fbc->prime[i+5]) == 0) 
				{						
					fb_offsets[++smooth_num] = i+5;	
					zShortDiv32(Q,fbc->prime[i+5],Q);			
				}
			}

			if (result & 0x2000)
			{
				while (zShortMod32(Q,fbc->prime[i+6]) == 0) 
				{						
					fb_offsets[++smooth_num] = i+6;	
					zShortDiv32(Q,fbc->prime[i+6],Q);			
				}
			}

			if (result & 0x8000)
			{
				while (zShortMod32(Q,fbc->prime[i+7]) == 0) 
				{						
					fb_offsets[++smooth_num] = i+7;	
					zShortDiv32(Q,fbc->prime[i+7],Q);			
				}
			}

			i += 8;
		}

		bound = sconf->factor_base->med_B;
		while ((uint32)i < bound)
		{

			uint32 result = 0;
			RESIEVE_4X_16BIT_MAX;

			if (result == 0)
			{
				i += 8;
				continue;
			}

			if (result & 0x2)
			{
				while (zShortMod32(Q,fbc->prime[i]) == 0) 
				{						
					fb_offsets[++smooth_num] = i;	
					zShortDiv32(Q,fbc->prime[i],Q);			
				}
			}

			if (result & 0x8)
			{
				while (zShortMod32(Q,fbc->prime[i+1]) == 0) 
				{						
					fb_offsets[++smooth_num] = i+1;	
					zShortDiv32(Q,fbc->prime[i+1],Q);			
				}
			}

			if (result & 0x20)
			{
				while (zShortMod32(Q,fbc->prime[i+2]) == 0) 
				{						
					fb_offsets[++smooth_num] = i+2;	
					zShortDiv32(Q,fbc->prime[i+2],Q);			
				}
			}

			if (result & 0x80)
			{
				while (zShortMod32(Q,fbc->prime[i+3]) == 0) 
				{						
					fb_offsets[++smooth_num] = i+3;	
					zShortDiv32(Q,fbc->prime[i+3],Q);			
				}
			}

			if (result & 0x200)
			{
				while (zShortMod32(Q,fbc->prime[i+4]) == 0) 
				{						
					fb_offsets[++smooth_num] = i+4;	
					zShortDiv32(Q,fbc->prime[i+4],Q);			
				}
			}

			if (result & 0x800)
			{
				while (zShortMod32(Q,fbc->prime[i+5]) == 0) 
				{						
					fb_offsets[++smooth_num] = i+5;	
					zShortDiv32(Q,fbc->prime[i+5],Q);			
				}
			}

			if (result & 0x2000)
			{
				while (zShortMod32(Q,fbc->prime[i+6]) == 0) 
				{						
					fb_offsets[++smooth_num] = i+6;	
					zShortDiv32(Q,fbc->prime[i+6],Q);			
				}
			}

			if (result & 0x8000)
			{
				while (zShortMod32(Q,fbc->prime[i+7]) == 0) 
				{						
					fb_offsets[++smooth_num] = i+7;	
					zShortDiv32(Q,fbc->prime[i+7],Q);			
				}
			}

			i += 8;

			/*
			if (result == 0)
			{
				i += 4;
				continue;
			}

			if (result & 0x8)
			{
				do 
				{						
					//if (zShortMod32(Q,fbc->prime[i]) != 0)
					//	printf("%u doesn't divide Q!\n",fbc->prime[i]);
					fb_offsets[++smooth_num] = i;	
					zShortDiv32(Q,fbc->prime[i],Q);			
				} while (zShortMod32(Q,fbc->prime[i]) == 0);
			}

			if (result & 0x80)
			{
				do 
				{						
					fb_offsets[++smooth_num] = i+1;	
					zShortDiv32(Q,fbc->prime[i+1],Q);			
				} while (zShortMod32(Q,fbc->prime[i+1]) == 0);
			}

			if (result & 0x800)
			{
				do 
				{					
					fb_offsets[++smooth_num] = i+2;	
					zShortDiv32(Q,fbc->prime[i+2],Q);			
				} while (zShortMod32(Q,fbc->prime[i+2]) == 0);
			}

			if (result & 0x8000)
			{
				do 
				{				
					fb_offsets[++smooth_num] = i+3;	
					zShortDiv32(Q,fbc->prime[i+3],Q);			
				} while (zShortMod32(Q,fbc->prime[i+3]) == 0);
			}

			i += 4;
			*/
		}
		

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		TF_STG4 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);
#endif

		dconf->smooth_num[report_num] = smooth_num;	

#else

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

		gettimeofday(&qs_timing_start, NULL);
#endif

#endif

	}

#if defined(SSE2_RESIEVING)
	#ifdef WIN32
		_aligned_free(corrections);
	#else
		free(corrections);
	#endif
#endif

	return;
}

void trial_divide_Q_siqs(uint32 report_num,  uint8 parity, 
						 uint32 poly_id, uint32 bnum, 
						 static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//we have flagged this sieve offset as likely to produce a relation
	//nothing left to do now but check and see.
	uint64 q64, f64;
	int j,it;
	uint32 prime;
	int smooth_num;
	uint32 *fb_offsets;
	uint32 polya_factors[20];
	sieve_fb *fb;
	uint32 offset, block_loc;
	z32 *Q;

	fb_offsets = &dconf->fb_offsets[report_num][0];
	smooth_num = dconf->smooth_num[report_num];
	Q = &dconf->Qvals[report_num];
	block_loc = dconf->reports[report_num];
	
#ifdef QS_TIMING
	gettimeofday(&qs_timing_start, NULL);
#endif

	offset = (bnum << BLOCKBITS) + block_loc;

	if (parity)
		fb = dconf->fb_sieve_n;
	else
		fb = dconf->fb_sieve_p;

	//check for additional factors of the a-poly factors
	//make a separate list then merge it with fb_offsets
	it=0;	//max 20 factors allocated for - should be overkill
	for (j = 0; (j < dconf->curr_poly->s) && (it < 20); j++)
	{
		//fbptr = fb + dconf->curr_poly->qlisort[j];
		//prime = fbptr->prime;
		prime = fb[dconf->curr_poly->qlisort[j]].prime;

		while ((zShortMod32(Q,prime) == 0) && (it < 20))
		{
			zShortDiv32(Q,prime,Q);
			polya_factors[it++] = dconf->curr_poly->qlisort[j];
		}
	}

	//check if it completely factored by looking at the unfactored portion in tmp
	if ((Q->size == 1) && (Q->val[0] < sconf->large_prime_max))
	{
		uint32 large_prime[2];
		
		large_prime[0] = Q->val[0];
		large_prime[1] = 1;

		//add this one
		buffer_relation(offset,large_prime,smooth_num+1,
			fb_offsets,poly_id,parity,dconf,polya_factors,it);

#ifdef QS_TIMING
	gettimeofday (&qs_timing_stop, NULL);
	qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

	TF_STG6 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
	free(qs_timing_diff);
#endif

		return;
	}

	if (sconf->use_dlp == 0)
		return;

	//quick check if Q is way too big for DLP (more than 64 bits)
	if (Q->size >= 3)
		return;

	q64 = (uint64)Q->val[0] + ((uint64)Q->val[1] << 32);
	sp642z(q64,&dconf->qstmp1);

	if ((zCompare(&dconf->qstmp1,&sconf->max_fb2) > 0) && 
		(zCompare(&dconf->qstmp1,&sconf->large_prime_max2) < 0))
	{	
		zShortSub(&dconf->qstmp1,1,&dconf->qstmp2);
		zModExp(&zTwo,&dconf->qstmp2,&dconf->qstmp1,&dconf->qstmp3);
		if (isOne(&dconf->qstmp3))
		{
#ifdef QS_TIMING
			gettimeofday (&qs_timing_stop, NULL);
			qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

			TF_STG6 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
			free(qs_timing_diff);
#endif
			return;
		}

		//zShortSub(&qstmp1,1,&qstmp2);
		//zModExp(&zThree,&qstmp2,&qstmp1,&qstmp3);
		//if (isOne(&qstmp3))
		//	return;
		
		//try to find a double large prime
		f64 = sp_shanks_loop(&dconf->qstmp1,sconf->obj);
		if (f64 > 1 && f64 != q64)
		{
			uint32 large_prime[2];

			large_prime[0] = (uint32)f64;
			large_prime[1] = (uint32)(q64 / (uint64)large_prime[0]);

			if (large_prime[0] < sconf->large_prime_max 
				&& large_prime[1] < sconf->large_prime_max)
			{
				//add this one
				buffer_relation(offset,large_prime,smooth_num+1,
					fb_offsets,poly_id,parity,dconf,polya_factors,it);
			}
		}
	}
#ifdef QS_TIMING
	gettimeofday (&qs_timing_stop, NULL);
	qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

	TF_STG6 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
	free(qs_timing_diff);
#endif
	
	return;
}

void buffer_relation(uint32 offset, uint32 *large_prime, uint32 num_factors, 
						  uint32 *fb_offsets, uint32 poly_id, uint32 parity,
						  dynamic_conf_t *conf, uint32 *polya_factors, 
						  uint32 num_polya_factors)
{
	//put this relations's info into a temporary buffer which
	//will get merged with other such buffers (if multi-threaded) and
	//dumped to file once the threads are joined.
	siqs_r *rel;
	uint32 i, j, k;

	//first check that this relation won't overflow the buffer
	if (conf->buffered_rels >= conf->buffered_rel_alloc)
	{
		printf("reallocating relation buffer\n");
		conf->relation_buf = (siqs_r *)realloc(conf->relation_buf, 
			conf->buffered_rel_alloc * 2 * sizeof(siqs_r));
		if (conf->relation_buf == NULL)
		{
			printf("error re-allocating temporary storage of relations\n");
			exit(-1);
		}
		conf->buffered_rel_alloc *= 2;
	}

	//then stick all the info in the buffer
	rel = conf->relation_buf + conf->buffered_rels;
	
	rel->sieve_offset = offset;
	rel->parity = parity;
	rel->poly_idx = poly_id;

	rel->fb_offsets = (uint32 *)malloc(
		(num_polya_factors + num_factors) * sizeof(uint32));

	//merge in extra factors of the apoly factors
	i = j = k = 0;
	while (k < num_factors && j < num_polya_factors) {
		if (fb_offsets[k] < polya_factors[j]) {
			rel->fb_offsets[i++] = fb_offsets[k++];
		}
		else if (fb_offsets[k] > polya_factors[j]) {
			rel->fb_offsets[i++] = polya_factors[j++];
		}
		else {
			rel->fb_offsets[i++] = fb_offsets[k++];
			rel->fb_offsets[i++] = polya_factors[j++];
		}
	}
	while (k < num_factors)
		rel->fb_offsets[i++] = fb_offsets[k++];
	while (j < num_polya_factors)
		rel->fb_offsets[i++] = polya_factors[j++];
		
	rel->num_factors = num_factors + num_polya_factors;
	rel->large_prime[0] = large_prime[0];
	rel->large_prime[1] = large_prime[1];

	conf->buffered_rels++;
	return;
}

void save_relation_siqs(uint32 offset, uint32 *large_prime, uint32 num_factors, 
						  uint32 *fb_offsets, uint32 poly_id, uint32 parity,
						  static_conf_t *conf)
{
	char buf[LINE_BUF_SIZE];
	fact_obj_t *obj = conf->obj;
	uint32 i, k;

	//store to file
	i = sprintf(buf, "R ");

	if (parity)
		i += sprintf(buf + i, "-%x ", offset);
	else
		i += sprintf(buf + i, "%x ", offset);
	
	i += sprintf(buf + i, "%x ", poly_id);
	
	k = 0;
	while (k < num_factors)
		i += sprintf(buf + i, "%x ", fb_offsets[k++]);

	if (large_prime[0] < large_prime[1])
		i += sprintf(buf + i, "L %x %x\n", large_prime[0], large_prime[1]);
	else
		i += sprintf(buf + i, "L %x %x\n", large_prime[1], large_prime[0]);

	qs_savefile_write_line(&obj->qs_obj.savefile, buf);

	/* for partial relations, also update the bookeeping for
	   tracking the number of fundamental cycles */

	if (large_prime[0] != large_prime[1]) {
		yafu_add_to_cycles(conf, obj->flags, large_prime[0], large_prime[1]);
		conf->num_cycles++;
	}
	else {
		conf->num_relations++;
	}

	return;
}

uint32 process_poly_a(static_conf_t *sconf)
{
	//given a poly a value, and some aux info about the factorization
	//generate all poly b values associated with that a and
	//store both a and all the b's in conf
	//also check that the generated polys are valid.

	//we will be reusing some routines here that are normally used
	//during sieving, and expect a dynamic_conf structure as input which
	//we don't have here.  So we need to create one for use in this 
	//routine only.
	dynamic_conf_t *dconf;
	int maxB,j,i;

	dconf = (dynamic_conf_t *)malloc(sizeof(dynamic_conf_t));

	//just initialize the stuff needed for processing poly's, namely, 
	//the curr_poly structure, the Bl array, and a few temp bigints.
	zInit(&dconf->qstmp1);
	zInit(&dconf->qstmp2);

	//this stuff changes with every new poly
	//allocate a polynomial structure which will hold the current
	//set of polynomial coefficients (a,b,c) and other info
	dconf->curr_poly = (siqs_poly *)malloc(sizeof(siqs_poly));
	zInit(&dconf->curr_poly->poly_a);
	zInit(&dconf->curr_poly->poly_b);
	zInit(&dconf->curr_poly->poly_c);
	dconf->curr_poly->qlisort = (int *)malloc(MAX_A_FACTORS*sizeof(int));
	dconf->curr_poly->gray = (char *) malloc( 65536 * sizeof(char));
	dconf->curr_poly->nu = (char *) malloc( 65536 * sizeof(char));

	//allocate the Bl array, space for MAX_Bl bigint numbers
	dconf->Bl = (z *)malloc(MAX_A_FACTORS * sizeof(z));
	for (i=0;i<MAX_A_FACTORS;i++)
		zInit(&dconf->Bl[i]);

	zCopy(&sconf->curr_a,&dconf->curr_poly->poly_a);
	
	//then compute all the 'b' poly's for this 'a'
	//and add them to the b-list

	//first, need to factorize this 'a'
	j = get_a_offsets(sconf->factor_base,dconf->curr_poly,&dconf->qstmp1);
	if (j)
	{
		//then this poly a might be corrupted - we couldn't factor it over the factor base
		free(dconf->curr_poly->gray);
		free(dconf->curr_poly->nu);
		free(dconf->curr_poly->qlisort);
		zFree(&dconf->curr_poly->poly_a);
		zFree(&dconf->curr_poly->poly_b);
		zFree(&dconf->curr_poly->poly_c);
		free(dconf->curr_poly);

		for (i=0;i<MAX_A_FACTORS;i++)
			zFree(&dconf->Bl[i]);
		free(dconf->Bl);

		//workspace bigints
		zFree(&dconf->qstmp1);
		zFree(&dconf->qstmp2);
		free(dconf);
		return 0;
	}

	//then initialize the gray code
	get_gray_code(dconf->curr_poly);

	//then compute all the Bl's
	//the first 'b' poly comes with computeBl
	computeBl(sconf,dconf);

	//compute how many 'b' values we can get from this 'a'
	maxB = 1 << (dconf->curr_poly->s - 1);

	//now we copy all the b coefficients over to sconf where they are
	//needed by the rest of the filtering routine.
	//make sure there is enough room for them.  curr_b could
	//currently be allocated for a number of poly smaller or bigger than
	//maxB

	//free any we won't be needing
	for (j = maxB; (uint32)j < sconf->bpoly_alloc; j++)
		zFree(&sconf->curr_b[j]);
	//reallocate the size of the array
	sconf->curr_b = (z *)realloc(sconf->curr_b, maxB * sizeof(z));
	//allocate any additional we need
	for (j = sconf->bpoly_alloc; j < maxB; j++)
		zInit(&sconf->curr_b[j]);
	sconf->bpoly_alloc = maxB;

	//generate all the b polys
	generate_bpolys(sconf,dconf,maxB);

	//we'll need to remember some things about the current poly,
	//so copy those over to sconf first...
	for (j=0; j<dconf->curr_poly->s; j++)
		sconf->curr_poly->qlisort[j] = dconf->curr_poly->qlisort[j];
	sconf->curr_poly->s = dconf->curr_poly->s;

	//then free the temp dynamic struct
	free(dconf->curr_poly->gray);
	free(dconf->curr_poly->nu);
	free(dconf->curr_poly->qlisort);
	zFree(&dconf->curr_poly->poly_a);
	zFree(&dconf->curr_poly->poly_b);
	zFree(&dconf->curr_poly->poly_c);
	free(dconf->curr_poly);

	for (i=0;i<MAX_A_FACTORS;i++)
		zFree(&dconf->Bl[i]);
	free(dconf->Bl);

	//workspace bigints
	zFree(&dconf->qstmp1);
	zFree(&dconf->qstmp2);
	free(dconf);

	return maxB;
}

int get_a_offsets(fb_list *fb, siqs_poly *poly, z *tmp)
{
	int j,k;
	fp_digit r;

	zCopy(&poly->poly_a,tmp);
	k=0;
	poly->s = 0;
	while (!(tmp->size == 1 && tmp->val[0] == 1) && (spSOEprimes[k] < 65536))
	{
		r = zShortMod(tmp,(fp_digit)spSOEprimes[k]);
		
		if (r != 0)
			k++;
		else
		{
			for (j=1;(uint32)j < fb->B;j++)
				if (fb->list->prime[j] == spSOEprimes[k])
					break;
			if ((uint32)j >= fb->B)
			{
				//then we didn't find this factor in the fb, thus the 
				//read in A is probably bad.
				printf("bad poly A encountered in savefile\n");
				return 1;
			}
			poly->qlisort[poly->s] = j;
			poly->s++;
			zShortDiv(tmp,(fp_digit)spSOEprimes[k],tmp);
		}
	}

	return 0;
}

void generate_bpolys(static_conf_t *sconf, dynamic_conf_t *dconf, int maxB)
{
	//given poly, which contains the first value of poly_b, and
	//info needed to generate the rest of the poly b's (namely,
	//the Bl vector), generate the rest of the Bl's and put
	//them in an output vector, which has been initialized 
	//by the callee.

	int numB=1;

	//iterate through all b's
	for ( ; numB < maxB; numB++)
	{
		zCopy(&dconf->curr_poly->poly_b,&sconf->curr_b[numB - 1]);
		dconf->numB = numB;
		nextB(dconf,sconf);
	}

	return;
}

int process_rel(char *substr, fb_list *fb, z *n,
				 static_conf_t *sconf, fact_obj_t *obj, siqs_r *rel)
{
	char *nextstr;
	uint32 lp[2];
	uint32 this_offset, this_id, this_num_factors, this_parity, this_val;
	uint32 fb_offsets[MAX_SMOOTH_PRIMES];
	int i,j,k, err_code = 0;

	rel->fb_offsets = (uint32 *)malloc(MAX_SMOOTH_PRIMES*sizeof(uint32));

	//read the offset
	this_offset = strtoul(substr,&nextstr,HEX);	//convert
	substr = nextstr;

	this_id = strtoul(substr,&nextstr,HEX);
	substr = nextstr;

	if (this_offset & 0x80000000)
	{
		this_parity = 1;
		this_offset = (~this_offset) + 1;
	}
	else
		this_parity = 0;

	j=0;
	do
	{
		this_val = strtoul(substr,&nextstr,HEX);
		substr = nextstr;
		fb_offsets[j] = this_val;
		j++;
	} while (substr[1] != 'L');
	this_num_factors = j;

	substr += 2;
	this_val = strtoul(substr,&nextstr,HEX);
	substr = nextstr;
	lp[0] = this_val;

	this_val = strtoul(substr,&nextstr,HEX);
	substr = nextstr;
	lp[1] = this_val;

	 //combine the factors of the sieve value with
	 //  the factors of the polynomial 'a' value; the 
	 //  linear algebra code has to know about both.
	 //  Because both lists are sorted, this is just
	 //  a merge operation 

	i = j = k = 0;
	while (i < (int)this_num_factors && j < sconf->curr_poly->s) {
		if (fb_offsets[i] < sconf->curr_poly->qlisort[j]) {
			rel->fb_offsets[k++] = fb_offsets[i++];
		}
		else if (fb_offsets[i] > sconf->curr_poly->qlisort[j]) {
			rel->fb_offsets[k++] = sconf->curr_poly->qlisort[j++];
		}
		else {
			rel->fb_offsets[k] = fb_offsets[i++];
			rel->fb_offsets[k+1] = sconf->curr_poly->qlisort[j++];
			k += 2;
		}
	}
	while (i < (int)this_num_factors)
		rel->fb_offsets[k++] = fb_offsets[i++];
	while (j < sconf->curr_poly->s)
		rel->fb_offsets[k++] = sconf->curr_poly->qlisort[j++];
	
	rel->sieve_offset = this_offset;
	rel->large_prime[0] = lp[0];
	rel->large_prime[1] = lp[1];
	rel->parity = this_parity;
	rel->num_factors = this_num_factors  + sconf->curr_poly->s;
	rel->poly_idx = this_id;

	if (!check_relation(&sconf->curr_a,
			&sconf->curr_b[rel->poly_idx],rel,fb,n))
	{
		if (lp[0] != lp[1]) {
			yafu_add_to_cycles(sconf, obj->flags, lp[0], lp[1]);
			sconf->num_cycles++;
		}
		else {
			sconf->num_relations++;
		}
	}
	else
	{
		err_code = 1;  
	}

	return err_code;	//error code, if there is one.
}

int restart_siqs(static_conf_t *sconf, dynamic_conf_t *dconf)
{
	int i,j;
	char *str, *substr;
	FILE *data;
	uint32 lp[2],pmax = sconf->large_prime_max / sconf->large_mult;
	//fact_obj_t *obj = sconf->obj;

	str = (char *)malloc(GSTR_MAXSIZE*sizeof(char));
	data = fopen(siqs_savefile,"r");
	i=0;
	j=0;
	
	if (data != NULL)
	{	
		fgets(str,1024,data);
		substr = str + 2;
		str2hexz(substr,&dconf->qstmp1);
		if (zCompare(&dconf->qstmp1,&sconf->n) == 0)
		{
			if (VFLAG > 1)
				printf("restarting siqs from saved data set\n");
			fflush(stdout);
			fflush(stderr);
			while (1)
			{
				//read a line
				if (feof(data))
					break;
				fgets(str,1024,data);
				substr = str + 2;

				if (str[0] == 'R')
				{	
					//process a relation
					//just trying to figure out how many relations we have
					//so read in the large primes and add to cycles
					substr = strchr(substr,'L');
					yafu_read_large_primes(substr,lp,lp+1);
					if (sconf->use_dlp)
					{
						if ((lp[0] > 1) && (lp[0] < pmax))
						{
							j++;
							continue;
						}
						if ((lp[1] > 1) && (lp[1] < pmax))
						{
							j++;
							continue;
						}
					}
					if (lp[0] != lp[1]) 
					{
						yafu_add_to_cycles(sconf, sconf->obj->flags, lp[0], lp[1]);
						sconf->num_cycles++;
					}
					else {
						sconf->num_relations++;
					}
				}
				else if (str[0] == 'A')
				{
					i++;
				}
			}

			sconf->num_r = sconf->num_relations + 
			sconf->num_cycles +
			sconf->components - sconf->vertices;

			if (VFLAG > 0)
			{
				printf("%d relations found: %d full + "
					"%d from %d partial\n",
					sconf->num_r,sconf->num_relations,
					sconf->num_cycles +
					sconf->components - sconf->vertices,
					sconf->num_cycles);
				printf("threw away %d relations with large primes too small\n",j);
				fflush(stdout);
				sconf->last_numfull = sconf->num_relations;
				sconf->last_numcycles = sconf->num_cycles;
			}

		}	
		fclose(data);
	}
	free(str);

	return 0;
}

#define QS_HASH_MULT ((uint32)(2654435761UL))
#define QS_HASH(a) (((a) * QS_HASH_MULT) >> (32 - QS_LOG2_CYCLE_HASH))

/**********************************************************
These 3 routines are used to add a relation to the cycle
tree every time one is found after sieving
(e.g., add_to_cycles is called in save_relation)
**********************************************************/
static qs_cycle_t *get_table_entry(qs_cycle_t *table, uint32 *hashtable,
				uint32 prime, uint32 new_entry_offset) {

	/* return a pointer to a unique qs_cycle_t specific
	   to 'prime'. The value of 'prime' is hashed and
	   the result used to index 'table'. If prime does
	   not appear in table, specify 'new_entry_offset'
	   as the table entry for prime.

	   Values of 'prime' which hash to the same offset
	   in hashtable are connected by a linked list of
	   offsets in 'table'. */

	uint32 offset, first_offset;
	qs_cycle_t *entry = NULL;

	first_offset = QS_HASH(prime);
	offset = hashtable[first_offset];
	
	/* follow the list of entries corresponding to
	   primes that hash to 'offset', stopping if either
	   the list runs out or we find an entry that's
	   dedicated to 'prime' already in the table */

	while (offset != 0) {
		entry = table + offset;
		if (entry->prime == prime)
			break;
		offset = entry->next;
	}

	/* if an entry was not found, initialize a new one */

	if (offset == 0) {
		entry = table + new_entry_offset;
		entry->next = hashtable[first_offset];
		entry->prime = prime;
		entry->data = 0;
		entry->count = 0;
		hashtable[first_offset] = new_entry_offset;
	}

	return entry;
}

/*--------------------------------------------------------------------*/
static uint32 add_to_hashtable(qs_cycle_t *table, uint32 *hashtable, 
			uint32 prime1, uint32 prime2, 
			uint32 default_table_entry, 
			uint32 *components, uint32 *vertices) {

	/* update the list of cycles to reflect the presence
	   of a partial relation with large primes 'prime1'
	   and 'prime2'.

	   There are three quantities to track, the number of
	   edges, components and vertices. The number of cycles 
	   in the graph is then e + c - v. There is one edge for
	   each partial relation, and one vertex for each prime
	   that appears in the graph (these are easy to count).

	   A connected component is a group of primes that are
	   all reachable from each other via edges in the graph.
	   All of the vertices that are in a cycle must belong
	   to the same connected component. Think of a component
	   as a tree of vertices; one prime is chosen arbitrarily
	   as the 'root' of that tree.
	   
	   The number of new primes added to the graph (0, 1, or 2)
	   is returned */

	uint32 root[2];
	uint32 root1, root2;
	uint32 i;
	uint32 num_new_entries = 0;

	/* for each prime */
	

	for (i = 0; i < 2; i++) {
		uint32 prime = ((i == 0) ? prime1 : prime2);
		uint32 offset; 
		qs_cycle_t *entry;

		/* retrieve the qs_cycle_t corresponding to that
		   prime from the graph (or allocate a new one) */

		entry = get_table_entry(table, hashtable,
					prime, default_table_entry);
		entry->count++;

		if (entry->data == 0) {

			/* if this prime has not occurred in the graph
			   before, increment the number of vertices and
			   the number of components, then make the table
			   entry point to itself. */
			//fprintf(globaltest,"new vertex: %u\n",prime);
			num_new_entries++;
			default_table_entry++;

			offset = entry - table;
			entry->data = offset;
			(*components)++;
			(*vertices)++;
		}
		else {
			/* the prime is already in the table, which
			   means we can follow a linked list of pointers
			   to other primes until we reach a qs_cycle_t
			   that points to itself. This last qs_cycle_t is
			   the 'root' of the connected component that
			   contains 'prime'. Save its value */
			//fprintf(globaltest,"data found: %u\n",prime);

			qs_cycle_t *first_entry, *next_entry;

			first_entry = entry;
			next_entry = table + entry->data;
			while (entry != next_entry) {
				entry = next_entry;
				next_entry = table + next_entry->data;
			}
				
			/* Also perform path compression: now that we
			   know the value of the root for this prime,
			   make all of the pointers in the primes we
			   visited along the way point back to this root.
			   This will speed up future root lookups */

			offset = entry->data;
			entry = first_entry;
			next_entry = table + entry->data;
			while (entry != next_entry) {
				entry->data = offset;
				entry = next_entry;
				next_entry = table + next_entry->data;
			}
		}

		root[i] = offset;
	}
				
	/* If the roots for prime1 and prime2 are different,
	   then they lie within separate connected components.
	   We're about to connect this edge to one of these
	   components, and the presence of the other prime
	   means that these two components are about to be
	   merged together. Hence the total number of components
	   in the graph goes down by one. */

	root1 = root[0];
	root2 = root[1];
	if (root1 != root2)
		(*components)--;
	
	/* This partial relation represents an edge in the
	   graph; we have to attach this edge to one or the
	   other of the connected components. Attach it to
	   the component whose representative prime is smallest;
	   since small primes are more common, this will give
	   the smaller root more edges, and will potentially
	   increase the number of cycles the graph contains */

	if (table[root1].prime < table[root2].prime)
		table[root2].data = root1;
	else
		table[root1].data = root2;
	
	return num_new_entries;
}

/*--------------------------------------------------------------------*/
void yafu_add_to_cycles(static_conf_t *conf, uint32 flags, uint32 prime1, uint32 prime2) {

	/* Top level routine for updating the graph of partial
	   relations */

	uint32 table_size = conf->cycle_table_size;
	uint32 table_alloc = conf->cycle_table_alloc;
	qs_cycle_t *table = conf->cycle_table;
	uint32 *hashtable = conf->cycle_hashtable;

	/* if we don't actually want to count cycles,
	   just increment the number of vertices. This is
	   equivalent to saying that the cycle-counting
	   code will never detect any cycles */
	
	if (flags & MSIEVE_FLAG_SKIP_QS_CYCLES) {
		conf->vertices++;
		return;
	}

	/* make sure there's room for new primes */

	if (table_size + 2 >= table_alloc) {
		table_alloc = conf->cycle_table_alloc = 2 * table_alloc;
		conf->cycle_table = (qs_cycle_t *)xrealloc(conf->cycle_table,
						table_alloc * sizeof(qs_cycle_t));
		table = conf->cycle_table;
	}

	conf->cycle_table_size += add_to_hashtable(table, hashtable, 
						prime1, prime2, 
						table_size, 
						&conf->components, 
						&conf->vertices);
}

/*******************************************************************************
These functions are used after sieving is complete to read in all
relations and find/optimize all the cycles
*******************************************************************************/
#define NUM_CYCLE_BINS 8

void yafu_qs_filter_relations(static_conf_t *sconf) {

	/* Perform all of the postprocessing on the list
	   of relations from the sieving phase. There are
	   two main jobs, reading in all the relations that
	   will be used and then determining the list of 
	   cycles in which partial relations appear. Care
	   should be taken to avoid wasting huge amounts of
	   memory */

	fact_obj_t *obj = sconf->obj;
	uint32 *hashtable = sconf->cycle_hashtable;
	qs_cycle_t *table = sconf->cycle_table;
	uint32 num_derived_poly;
	uint32 *final_poly_index;
	uint32 num_relations, num_cycles, num_poly;
	qs_la_col_t *cycle_list;
	siqs_r *relation_list;

	uint32 i, passes, start;
	uint32 curr_a_idx, curr_poly_idx, curr_rel; 
	uint32 curr_expected, curr_saved, curr_cycle; 
	uint32 total_poly_a;
	uint32 poly_saved;
	uint32 cycle_bins[NUM_CYCLE_BINS+1] = {0};
	char buf[LINE_BUF_SIZE];
	char *subbuf;
	uint32 last_id;
	int first;

 	/* Rather than reading all the relations in and 
	   then removing singletons, read only the large 
	   primes of each relation into an initial list,
	   remove the singletons, and then only read in
	   the relations that survive. This avoids reading
	   in useless relations (and usually the polynomials 
	   they would need) */

	i = 0;
	total_poly_a = 0;

	/* skip over the first line */
	qs_savefile_open(&obj->qs_obj.savefile, SAVEFILE_READ);
	qs_savefile_read_line(buf, sizeof(buf), &obj->qs_obj.savefile);

	//we don't know beforehand how many rels to expect, so start
	//with some amount and allow it to increase as we read them
	relation_list = (siqs_r *)xmalloc(10000 * sizeof(siqs_r));
	curr_rel = 10000;
	while (!qs_savefile_eof(&obj->qs_obj.savefile)) {
		char *start;

		switch (buf[0]) {
		case 'A':
			total_poly_a++;
			break;

		case 'R':
			start = strchr(buf, 'L');
			if (start != NULL) {
				uint32 prime1, prime2;
				yafu_read_large_primes(start, &prime1, &prime2);
				if (i == curr_rel) {
					curr_rel = 3 * curr_rel / 2;
					relation_list = (siqs_r *)xrealloc(
							relation_list,
							curr_rel *
							sizeof(siqs_r));
				}
				relation_list[i].poly_idx = i;
				relation_list[i].large_prime[0] = prime1;
				relation_list[i].large_prime[1] = prime2;
				i++;
			}
			break;
		}

		qs_savefile_read_line(buf, sizeof(buf), &obj->qs_obj.savefile);
	}
	
	num_relations = i;
	num_relations = qs_purge_singletons(obj, relation_list, num_relations,
					table, hashtable);

	relation_list = (siqs_r *)xrealloc(relation_list, num_relations * 
							sizeof(siqs_r));

	/* Now we know how many relations to expect. Also
	   initialize the lists of polynomial 'a' and 'b' values */

	num_poly = 10000;

	sconf->poly_list = (poly_t *)xmalloc(num_poly * sizeof(poly_t));
	sconf->total_poly_a = total_poly_a;
	sconf->poly_a_list = (z *)xmalloc(total_poly_a * sizeof(z));
	for (i=0; i<total_poly_a; i++)
		zInit(&sconf->poly_a_list[i]);

	final_poly_index = (uint32 *)xmalloc(1024 * 
						sizeof(uint32));
	
	/* initialize the running counts of relations and
	   polynomials */

	i = 0;
	last_id = 0;
	curr_expected = 0;
	curr_saved = 0;
	curr_rel = (uint32)(-1);
	curr_poly_idx = (uint32)(-1);
	curr_a_idx = (uint32)(-1);
	poly_saved = 0;
	sconf->poly_list_alloc = 0;
	if (VFLAG > 0)
		printf("attempting to read %u relations\n", num_relations);

	/* Read in the relations and the polynomials they use
	   at the same time. */

	qs_savefile_rewind(&obj->qs_obj.savefile);

	first = 1;
	while (curr_expected < num_relations) {
		char *tmp;
		uint32 bad_A_val = 0;
		siqs_r *r;

		/* read in the next entity */
		if (qs_savefile_eof(&obj->qs_obj.savefile))
			break;
		qs_savefile_read_line(buf, sizeof(buf), &obj->qs_obj.savefile);

		switch (buf[0]) {
		case 'A':
			/* Read in a new 'a' value */
			/* build all of the 'b' values associated with it */
			subbuf = buf + 2;	//skip the A and a space
			curr_a_idx++;
			str2hexz(subbuf,&sconf->curr_a);
			num_derived_poly = process_poly_a(sconf);
			if (num_derived_poly == 0)
			{
				//this is an error indicating a bad poly a.  skip all relations
				//until we see the next A
				bad_A_val = 1;
				continue;
			}
			else
				bad_A_val = 0;

			zCopy(&sconf->curr_a,sconf->poly_a_list + curr_a_idx);	

			/* all 'b' values start off unused */
			final_poly_index = (uint32 *)xrealloc(final_poly_index,
				num_derived_poly * sizeof(uint32));
			memset(final_poly_index, -1, num_derived_poly *
							sizeof(uint32));

			break;

		case 'R':
			/* handle a new relation. First find the 
			   large primes; these will determine
	     		   if a relation is full or partial */
			if (bad_A_val)
				continue;

			tmp = strchr(buf, 'L');
			if (tmp == NULL)
				break;

			/* Check if this relation is needed. If it
			   survived singleton removal then its 
			   ordinal ID will be in the next entry 
			   of relation_list. 
		   
			   First move up the index of relation_list 
			   until the relation index to check is >= 
			   the one we have (it may have gotten behind 
			   because relations were corrupted) */

			curr_rel++;

			while (curr_expected < num_relations &&
				relation_list[curr_expected].poly_idx <
						curr_rel) {
				curr_expected++;
			}

			/* now check if the relation should be saved */

			if (curr_expected >= num_relations ||
			    relation_list[curr_expected].poly_idx != curr_rel)
				break;

			curr_expected++;

			/* convert the ASCII text of the relation to a
		    relation_t, verifying correctness in the process */
			r = relation_list + curr_saved;
			subbuf = buf + 2;	//skip over the R and a space
			if (process_rel(subbuf, sconf->factor_base,
				&sconf->n, sconf, sconf->obj, r)) {
				logprint(obj->logfile, "failed to read relation %d\n", 
							curr_expected - 1);
				if (VFLAG > 0)
					printf("failed to read relation %d\n", 
							curr_expected - 1);
				break;
			}

			curr_saved++;

			/* if necessary, save the b value corresponding 
			   to this relation */

			if (final_poly_index[r->poly_idx] == (uint32)(-1)) {

				if (i == num_poly) {
					num_poly *= 2;
					sconf->poly_list = (poly_t *)xrealloc(
							sconf->poly_list,
							num_poly *
							sizeof(poly_t));
				}

				sconf->poly_list[i].a_idx = curr_a_idx;
				zInit(&(sconf->poly_list[i].b));
				zCopy(&(sconf->curr_b[r->poly_idx]), 
					&(sconf->poly_list[i].b));
				sconf->poly_list_alloc++;
				final_poly_index[r->poly_idx] = i;
				r->poly_idx = i++;
			}
			else {
				r->poly_idx = final_poly_index[r->poly_idx];
			}

			break;  /* done with this relation */
		}
	}

	/* update the structures with the counts of relations
	   and polynomials actually recovered */

	num_relations = curr_saved;
	logprint(obj->logfile, "recovered %u relations\n", num_relations);
	logprint(obj->logfile, "recovered %u polynomials\n", i);
	if (VFLAG > 0)
	{
		printf("recovered %u relations\n", num_relations);
		printf("recovered %u polynomials\n", i);
	}

	qs_savefile_close(&obj->qs_obj.savefile);
	free(final_poly_index);
	sconf->poly_list = (poly_t *)xrealloc(sconf->poly_list,
					   i * sizeof(poly_t));

	/* begin the cycle generation process by purging
	   duplicate relations. For the sake of consistency, 
	   always rebuild the graph afterwards */

	num_relations = qs_purge_duplicate_relations(obj, 
				relation_list, num_relations);

	memset(hashtable, 0, sizeof(uint32) << QS_LOG2_CYCLE_HASH);
	sconf->vertices = 0;
	sconf->components = 0;
	sconf->cycle_table_size = 1;

	for (i = 0; i < num_relations; i++) {
		siqs_r *r = relation_list + i;
		if (r->large_prime[0] != r->large_prime[1]) {
			yafu_add_to_cycles(sconf, sconf->obj->flags, r->large_prime[0], 
					r->large_prime[1]);
		}
	}
	
	/* compute the number of cycles to expect. Note that
	   this number includes cycles from both full and partial
	   relations (the cycle for a full relation is trivial) */

	num_cycles = num_relations + sconf->components - sconf->vertices;

	/* The idea behind the cycle-finding code is this: the 
	   graph is composed of a bunch of connected components, 
	   and each component contains one or more cycles. To 
	   find the cycles, you build the 'spanning tree' for 
	   each component.

	   Think of the spanning tree as a binary tree; there are
	   no cycles in it because leaves are only connected to a
	   common root and not to each other. Any time you connect 
	   together two leaves of the tree, though, a cycle is formed.
	   So, for a spanning tree like this:

	         1
	         o
		/ \
	    2  o   o  3
	      / \   \
	     o   o   o
	     4   5   6

	   if you connect leaves 4 and 5 you get a cycle (4-2-5). If
	   you connect leaves 4 and 6 you get another cycle (4-2-1-3-6)
	   that will reuse two of the nodes in the first cycle. It's
	   this reuse that makes double large primes so powerful.

	   For our purposes, every edge in the tree above represents
	   a partial relation. Every edge that would create a cycle
	   comes from another partial relation. So to find all the cycles,
	   you begin with the roots of all of the connected components,
	   and then iterate through the list of partial relations until 
	   all have been 'processed'. A partial relation is considered 
	   processed when one or both of its primes is in the tree. If 
	   one prime is present then the relation gets added to the tree; 
	   if both primes are present then the relation creates one cycle 
	   but is *not* added to the tree. 
	   
	   It's really great to see such simple ideas do something so
	   complicated as finding cycles (and doing it very quickly) */

	/* First traverse the entire graph and remove any vertices
	   that are not the roots of connected components (i.e.
	   remove any primes whose cycle_t entry does not point
	   to itself */

	for (i = 0; i < (1 << QS_LOG2_CYCLE_HASH); i++) {
		uint32 offset = hashtable[i];

		while (offset != 0) {
			qs_cycle_t *entry = table + offset;
			if (offset != entry->data)
				entry->data = 0;
			offset = entry->next;
		}
	}

	logprint(obj->logfile, "attempting to build %u cycles\n", num_cycles);
	if (VFLAG > 0)
		printf("attempting to build %u cycles\n", num_cycles);
	cycle_list = (qs_la_col_t *)xmalloc(num_cycles * sizeof(qs_la_col_t));

	/* keep going until either all cycles are found, all
	   relations are processed, or cycles stop arriving. 
	   Normally these conditions all occur at the same time */

	for (start = passes = curr_cycle = 0; start < num_relations && 
			curr_cycle < num_cycles; passes++) {

		/* The list of relations up to index 'start' is con-
		   sidered processed. For all relations past that... */

		uint32 start_cycles = curr_cycle;

		for (i = start; i < num_relations &&
				curr_cycle < num_cycles; i++) {

			qs_cycle_t *entry1, *entry2;
			siqs_r rtmp = relation_list[i];
			
			if (rtmp.large_prime[0] == rtmp.large_prime[1]) {

				/* this is a full relation, and forms a
				   cycle just by itself. Move it to position 
				   'start' of the relation list and increment 
				   'start'. The relation is now frozen at 
				   that position */

				qs_la_col_t *c = cycle_list + curr_cycle++;
				relation_list[i] = relation_list[start];
				relation_list[start] = rtmp;

				/* build a trivial cycle for the relation */

				c->cycle.num_relations = 1;
				c->cycle.list = (uint32 *)
						xmalloc(sizeof(uint32));
				c->cycle.list[0] = start++;
				continue;
			}

			/* retrieve the cycle_t entries associated
			   with the large primes in relation r. */

			entry1 = get_table_entry(table, hashtable, 
						rtmp.large_prime[0], 0);
			entry2 = get_table_entry(table, hashtable, 
						rtmp.large_prime[1], 0);

			/* if both vertices do not point to other
			   vertices, then neither prime has been added
			   to the graph yet, and r must remain unprocessed */

			if (entry1->data == 0 && entry2->data == 0)
				continue;

			/* if one or the other prime is part of the
			   graph, add r to the graph. The vertex not in
			   the graph points to the vertex that is, and
			   this entry also points to the relation that
			   is associated with rtmp.

			   If both primes are in the graph, recover the
			   cycle this generates */

			if (entry1->data == 0) {
				entry1->data = entry2 - table;
				entry1->count = start;
			}
			else if (entry2->data == 0) {
				entry2->data = entry1 - table;
				entry2->count = start;
			}
			else {
				qs_la_col_t *c = cycle_list + curr_cycle;
				c->cycle.list = NULL;
				qs_enumerate_cycle(obj, c, table, entry1,
						entry2, start);
				if (c->cycle.list)
					curr_cycle++;
			}

			/* whatever happened above, the relation is
			   processed now; move it to position 'start'
			   of the relation list and increment 'start'.
			   The relation is now frozen at that position */

			relation_list[i] = relation_list[start];
			relation_list[start++] = rtmp;
		}

		/* If this pass did not find any new cycles, then
		   we've reached steady state and are finished */

		if (curr_cycle == start_cycles)
			break;
	}
	num_cycles = curr_cycle;

	logprint(obj->logfile, "found %u cycles in %u passes\n", num_cycles, passes);
	if (VFLAG > 0)
		printf("found %u cycles in %u passes\n", num_cycles, passes);
	
	/* sort the list of cycles so that the cycles with
	   the largest number of relations will come last. 
	   If the linear algebra code skips any cycles it
	   can easily skip the most dense cycles */

	qsort(cycle_list, (size_t)num_cycles, sizeof(qs_la_col_t), yafu_sort_cycles);

	sconf->relation_list = relation_list;
	sconf->num_relations = num_relations;
	sconf->cycle_list = cycle_list;
	sconf->num_cycles = num_cycles;

	/* print out a histogram of cycle lengths for infor-
	   mational purposes */

	for (i = 0; i < num_cycles; i++) {
		num_relations = cycle_list[i].cycle.num_relations;

		if (num_relations >= NUM_CYCLE_BINS)
			cycle_bins[NUM_CYCLE_BINS]++;
		else
			cycle_bins[num_relations - 1]++;
	}
	logprint(obj->logfile, "distribution of cycle lengths:\n");
	if (VFLAG > 0)
		printf("distribution of cycle lengths:\n");
	for (i = 0; i < NUM_CYCLE_BINS; i++) {
		if (cycle_bins[i]) {
			logprint(obj->logfile, "   length %d : %d\n", 
					i + 1, cycle_bins[i]);
			if (VFLAG > 0)
				printf("   length %d : %d\n", 
					i + 1, cycle_bins[i]);
		}
	}
	if (cycle_bins[i])
	{
		logprint(obj->logfile, "   length %u+: %u\n", i + 1, cycle_bins[i]);
		if (VFLAG > 0)
			printf("   length %u+: %u\n", i + 1, cycle_bins[i]);
	}
	logprint(obj->logfile, "largest cycle: %u relations\n",
			cycle_list[num_cycles-1].cycle.num_relations);
	if (VFLAG > 0)
		printf("largest cycle: %u relations\n",
			cycle_list[num_cycles-1].cycle.num_relations);
}

void pull_large_primes()
{
	FILE *in, *out;
	char buf[1000];
	z tmp;

	zInit(&tmp);
	in = fopen(siqs_savefile,"r");
	out = fopen("dlp_composites.txt","a");

	while (!feof(in)) {
		char *start;
		fgets(buf,1000,in);

		switch (buf[0]) {
		case 'R':
			start = strchr(buf, 'L');
			if (start != NULL) {
				uint32 prime1, prime2;
				yafu_read_large_primes(start, &prime1, &prime2);

				if ((prime1 != 1) && (prime2 != 1))
				{
					sp2z(prime1,&tmp);
					zShortMul(&tmp,prime2,&tmp);
					fprintf(out,"%s,%u,%u\n",z2decstr(&tmp,&gstr1),prime1,prime2);
				}
			}
			break;
		}
	}

	fclose(in);
	fclose(out);
	zFree(&tmp);
	return;
}

static int compare_relations(const void *x, const void *y) {

	/* Callback used to sort a list of sieve relations.
	   Sorting is by size of large primes, then by number
	   of factors, then by factor values. Only the first
	   rule is needed for ordinary MPQS, but with self-
	   initialization we have to detect duplicate relations,
	   and this is easier if they are sorted as described */

	siqs_r *xx = (siqs_r *)x;
	siqs_r *yy = (siqs_r *)y;
	uint32 i;

	if (xx->large_prime[1] > yy->large_prime[1])
		return 1;
	if (xx->large_prime[1] < yy->large_prime[1])
		return -1;

	if (xx->large_prime[0] > yy->large_prime[0])
		return 1;
	if (xx->large_prime[0] < yy->large_prime[0])
		return -1;

	if (xx->num_factors > yy->num_factors)
		return 1;
	if (xx->num_factors < yy->num_factors)
		return -1;

	for (i = 0; i < xx->num_factors; i++) {
		if (xx->fb_offsets[i] > yy->fb_offsets[i])
			return 1;
		if (xx->fb_offsets[i] < yy->fb_offsets[i])
			return -1;
	}
	return 0;
}

/*--------------------------------------------------------------------*/
uint32 qs_purge_duplicate_relations(fact_obj_t *obj,
				siqs_r *rlist, 
				uint32 num_relations) {

	uint32 i, j;
	
	/* remove duplicates from rlist */

	if (num_relations < 2)
		return num_relations;

	qsort(rlist, (size_t)num_relations,
		sizeof(siqs_r), compare_relations);
	
	for (i = 1, j = 0; i < num_relations; i++) {
		if (compare_relations(rlist + j, rlist + i) == 0)
		{
			free(rlist[i].fb_offsets);
		}
		else
			rlist[++j] = rlist[i];
	}

	j++;
	if (j != num_relations)
	{
		logprint(obj->logfile, "freed %d duplicate relations\n", 
					num_relations - j);
		printf("freed %d duplicate relations\n", 
					num_relations - j);
	}
	return j;
}

void yafu_read_large_primes(char *buf, uint32 *prime1, uint32 *prime2) {

	char *next_field;
	uint32 p1, p2;

	*prime1 = p1 = 1;
	*prime2 = p2 = 2;
	if (*buf != 'L')
		return;

	buf++;
	while (isspace(*buf))
		buf++;
	if (isxdigit(*buf)) {
		p1 = strtoul(buf, &next_field, 16);
		buf = next_field;
	}
	else {
		return;
	}

	while (isspace(*buf))
		buf++;
	if (isxdigit(*buf))
		p2 = strtoul(buf, &next_field, 16);
	
	if (p1 < p2) {
		*prime1 = p1;
		*prime2 = p2;
	}
	else {
		*prime1 = p2;
		*prime2 = p1;
	}
}

uint32 qs_purge_singletons(fact_obj_t *obj, siqs_r *list, 
				uint32 num_relations,
				qs_cycle_t *table, uint32 *hashtable) {
	
	/* given a list of relations and the graph from the
	   sieving stage, remove any relation that contains
	   a prime that only occurs once in the graph. Because
	   removing a relation removes a second prime as well,
	   this process must be iterated until no more relations
	   are removed */

	uint32 num_left;
	uint32 i, j, k;
	uint32 passes = 0;

	if (VFLAG > 0)
		printf("begin with %u relations\n", num_relations);
	logprint(obj->logfile, "begin with %u relations\n", num_relations);

	do {
		num_left = num_relations;

		/* for each relation */

		for (i = j = 0; i < num_relations; i++) {
			siqs_r *r = list + i;
			uint32 prime;
			qs_cycle_t *entry;

			/* full relations always survive */

			if (r->large_prime[0] == r->large_prime[1]) {
				list[j++] = list[i];
				continue;
			}

			/* for each prime in that relation */

			for (k = 0; k < 2; k++) {
				prime = r->large_prime[k];
				entry = get_table_entry(table, hashtable,
							prime, 0);

				/* if the relation is due to be removed,
				   decrement the count of its other
				   primes in the graph. The following is
				   specialized for two primes */

				if (entry->count < 2) {
					prime = r->large_prime[k ^ 1];
					entry = get_table_entry(table, 
								hashtable, 
								prime, 0);
					entry->count--;
					break;
				}
			}

			if (k == 2)
				list[j++] = list[i];
		}
		num_relations = j;
		passes++;

	} while (num_left != num_relations);
				
	logprint(obj->logfile, "reduce to %u relations in %u passes\n", 
				num_left, passes);
	if (VFLAG > 0)
		printf("reduce to %u relations in %u passes\n", 
				num_left, passes);
	return num_left;
}

/*--------------------------------------------------------------------*/
void qs_enumerate_cycle(fact_obj_t *obj, 
			    qs_la_col_t *c, 
			    qs_cycle_t *table,
			    qs_cycle_t *entry1, qs_cycle_t *entry2,
			    uint32 final_relation) {

	/* given two entries out of the hashtable, corresponding
	   to two distinct primes, generate the list of relations
	   that participate in the cycle that these two primes
	   have just created. final_relation is the relation
	   to which the two primes belong, and the completed cycle
	   is packed into 'c' */

	uint32 traceback1[100];
	uint32 traceback2[100];
	uint32 num1, num2;
	uint32 i, j;

	/* Follow each cycle_t back up the graph until
	   the root component for this cycle is reached.
	   For each prime encountered along the way, save
	   the offset of the relation containing that prime */

	num1 = 0;
	while (entry1 != table + entry1->data) {
		if (num1 >= 100) {
			logprint(obj->logfile, "warning: cycle too long, "
					"skipping it\n");
			printf("warning: cycle too long, "
					"skipping it\n");
			return;
		}
		traceback1[num1++] = entry1->count;
		entry1 = table + entry1->data;
	}

	num2 = 0;
	while (entry2 != table + entry2->data) {
		if (num2 >= 100) {
			logprint(obj->logfile, "warning: cycle too long, "
					"skipping it\n");
			printf("warning: cycle too long, "
					"skipping it\n");
			return;
		}
		traceback2[num2++] = entry2->count;
		entry2 = table + entry2->data;
	}

	/* Now walk backwards through the lists, until
	   either one list runs out or a relation is
	   encountered that does not appear in both lists */

	while (num1 > 0 && num2 > 0) {
		if (traceback1[num1 - 1] != traceback2[num2 - 1])
			break;
		num1--; 
		num2--;
	}

	/* Now that we know how many relations are in the
	   cycle, allocate space to remember them */

	c->cycle.num_relations = num1 + num2 + 1;
	c->cycle.list = (uint32 *)xmalloc(c->cycle.num_relations * 
					sizeof(uint32));
	
	/* Combine the two lists of relations */
	for (i = 0; i < num1; i++)
		c->cycle.list[i] = traceback1[i];

	for (j = 0; j < num2; j++, i++)
		c->cycle.list[i] = traceback2[j];

	/* Add the relation that created the cycle in the
	   first place */

	c->cycle.list[i] = final_relation;
}

/*--------------------------------------------------------------------*/
int yafu_sort_cycles(const void *x, const void *y) {
	qs_la_col_t *xx = (qs_la_col_t *)x;
	qs_la_col_t *yy = (qs_la_col_t *)y;

	/* Callback for sorting a list of cycles by the
	   number of relations each contains */

	return xx->cycle.num_relations - yy->cycle.num_relations;
}
