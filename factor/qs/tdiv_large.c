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

this file contains code implementing 5)



*/


#if defined(GCC_ASM32X) || defined(GCC_ASM64X) || defined(__MINGW32__)
	//these compilers support SIMD 
	#define SCAN_CLEAN asm volatile("emms");	

    #if defined(USE_AVX2)

        #define SCAN_16X_VEC_b			\
			asm volatile (			\
				"vmovdqa (%2), %%xmm0	\n\t"	/*move mask into xmm0*/	\
				"vmovdqa (%1), %%xmm1	\n\t"	/*move 16 bptr locations into xmm regs*/	\
				"vmovdqa 16(%1), %%xmm2	\n\t"		\
				"vpcmpeqw %%xmm0, %%xmm1, %%xmm1	\n\t"	/*compare to mask*/	\
				"vmovdqa 32(%1), %%xmm3	\n\t"		\
				"vpcmpeqw %%xmm0, %%xmm2, %%xmm2	\n\t"		\
				"vmovdqa 48(%1), %%xmm4	\n\t"		\
				"vpcmpeqw %%xmm0, %%xmm3, %%xmm3	\n\t"		\
				"vpcmpeqw %%xmm0, %%xmm4, %%xmm4	\n\t"		\
                "vpmovmskb %%xmm1, %%r8   \n\t"		/* 1st 4 comparisons in 16 bits of r8  */		\
                "vpmovmskb %%xmm2, %%r9   \n\t"		/* 2nd 4 comparisons in 16 bits of r9  */		\
                "vpmovmskb %%xmm3, %%r10   \n\t"		/* 3rd 4 comparisons in 16 bits of r9  */		\
                "vpmovmskb %%xmm4, %%r11   \n\t"		/* 4th 4 comparisons in 16 bits of r9  */		\
                "salq $16, %%r9		\n\t"			/*  */ \
                "salq $32, %%r10		\n\t"			/*  */ \
                "salq $48, %%r11		\n\t"			/*  */ \
                "orq	%%r11,%%r10		\n\t"		/* r8 now holds 16 comparisons in 64 bits */ \
                "orq	%%r9,%%r8		\n\t"		/* r8 now holds 8 comparisons in 32 bits */ \                                
                "orq	%%r10,%%r8		\n\t"		/* r8 now holds 12 comparisons in 48 bits */ \                                
                "movq   $0x2222222222222222,%%r9    \n\t" /* clear the bytemask results we don't care about */ \
                "andq   %%r9,%%r8    \n\t"                /* clear the bytemask results we don't care about */ \                
                "movl	%0,%%r11d		\n\t"		/* initialize count of set bits */ \
                "xorq	%%r10,%%r10		\n\t"		/* initialize bit scan offset */ \
                "1:			\n\t"					/* top of bit scan loop */ \
                "bsfq	%%r8,%%rcx		\n\t"		/* put least significant set bit index into rcx */ \
                "jz 2f	\n\t"						/* jump out if zero (no hits).  high percentage. */ \
                "addq	%%rcx,%%r10	\n\t"			/* add in the offset of this index */ \
                "movq   %%r10,%%r9 \n\t" \
                "sarq   $2,%%r9 \n\t"               /* translate to offset within bptr */ \
                "addl   %4,%%r9d \n\t"   \
                "movw	%%r9w, (%3, %%r11, 2) \n\t"		/* put the bit index into the output buffer */ \
                "shrq	%%cl,%%r8	\n\t"			/* shift the bit scan register up to the bit we just processed */ \
                "incl	%%r11d		\n\t"			/* increment the count of set bits */ \
                "incq	%%r10		\n\t"			/* increment the index */ \
                "shrq	$1, %%r8 \n\t"				/* clear the bit */ \
                "jmp 1b		\n\t"					/* loop if so */ \
                "2:		\n\t"						/*  */ \
                "movl	%%r11d, %0 \n\t"			/* return the count of set bits */ \
                : "+r"(result)						\
                : "r"(bptr + j), "r"(mask), "r"(buffer), "r"(j)	\
                : "r8", "r9", "r10", "r11", "rcx", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "cc", "memory");	

        #define SCAN_16X_VEC			\
			asm volatile (			\
				"vmovdqa (%2), %%ymm0	\n\t"	/*move mask into xmm0*/	\
				"vmovdqa (%1), %%ymm1	\n\t"	/*move 16 bptr locations into xmm regs*/	\
				"vmovdqa 32(%1), %%ymm2	\n\t"		\
				"vpcmpeqw %%ymm0, %%ymm1, %%ymm1	\n\t"	/*compare to mask*/	\
				"vpcmpeqw %%ymm0, %%ymm2, %%ymm2	\n\t"		\
                "vpor   %%ymm1, %%ymm2, %%ymm3 \n\t" \
                "vpmovmskb %%ymm1, %%r8   \n\t"		/* 1st 4 comparisons in 16 bits of r8  */		\
                "vpmovmskb %%ymm2, %%r9   \n\t"		/* 2nd 4 comparisons in 16 bits of r9  */		\                
                "vpmovmskb %%ymm3, %%r10   \n\t"		/* 2nd 4 comparisons in 16 bits of r9  */		\
                "testq %%r10, %%r10 \n\t"			/* AND, and set ZF */ \
			    "jz 3f	\n\t"						/* jump out if zero (no hits).  high percentage. */ \
                "salq $32, %%r9		\n\t"			/*  */ \
                "orq	%%r9,%%r8		\n\t"		/* r8 now holds 8 comparisons in 32 bits */ \
                "movq   $0x2222222222222222,%%r9    \n\t" /* clear the bytemask results we don't care about */ \
                "movl   %0,%%r11d \n\t" \
                "andq   %%r9,%%r8    \n\t"                /* clear the bytemask results we don't care about */ \                
                "xorq	%%r10,%%r10		\n\t"		/* initialize bit scan offset */ \
                "1:			\n\t"					/* top of bit scan loop */ \
                "bsfq	%%r8,%%rcx		\n\t"		/* put least significant set bit index into rcx */ \
                "jz 2f	\n\t"						/* jump out if zero (no hits).  high percentage. */ \
                "addq	%%rcx,%%r10	\n\t"			/* add in the offset of this index */ \
                "movq   %%r10,%%r9 \n\t" \
                "shrq   $2,%%r9 \n\t"               /* translate to offset within bptr */ \
                "addl   %4,%%r9d \n\t"   \
                "movw	%%r9w, (%3, %%r11, 2) \n\t"		/* put the bit index into the output buffer */ \
                "shrq	%%cl,%%r8	\n\t"			/* shift the bit scan register up to the bit we just processed */ \
                "incl	%%r11d		\n\t"			/* increment the count of set bits */ \
                "incq	%%r10		\n\t"			/* increment the index */ \
                "shrq	$1, %%r8 \n\t"				/* clear the bit */ \
                "jmp 1b		\n\t"					/* repeat */ \
                "2:		\n\t"						/*  */ \
                "movl	%%r11d, %0 \n\t"			/* return the count of set bits */ \
                "3:     \n\t" \
                : "+r"(result)						\
                : "r"(bptr + j), "r"(mask), "r"(buffer), "r"(j)	\
                : "r8", "r9", "r10", "r11", "rcx", "xmm0", "xmm1", "xmm2", "xmm3", "cc", "memory");



	#elif defined(HAS_SSE2)
		// top level sieve scanning with SSE2
        // the block_loc that we are looking for is in the
        // bottom half of each 32-bit value of bptr, so each
        // 128-bit xmm register holds 4 locations to search.
        // we load 4 registers worth, 16 locations total, 
        // perform the test on each, and OR all of the results together.

        #define SCAN_16X_VEC			\
			asm volatile (			\
				"movdqa (%2), %%xmm0	\n\t"	/*move mask into xmm0*/	\
				"movdqa (%1), %%xmm1	\n\t"	/*move 16 bptr locations into xmm regs*/	\
				"movdqa 16(%1), %%xmm2	\n\t"		\
				"pcmpeqw %%xmm0, %%xmm1	\n\t"	/*compare to mask*/	\
				"movdqa 32(%1), %%xmm3	\n\t"		\
				"pcmpeqw %%xmm0, %%xmm2	\n\t"		\
				"movdqa 48(%1), %%xmm4	\n\t"		\
				"pcmpeqw %%xmm0, %%xmm3	\n\t"		\
				"pcmpeqw %%xmm0, %%xmm4	\n\t"		\
                "pmovmskb %%xmm1, %%r8   \n\t"		/* 1st 4 comparisons in 16 bits of r8  */		\
                "pmovmskb %%xmm2, %%r9   \n\t"		/* 2nd 4 comparisons in 16 bits of r9  */		\
                "pmovmskb %%xmm3, %%r10   \n\t"		/* 3rd 4 comparisons in 16 bits of r9  */		\
                "pmovmskb %%xmm4, %%r11   \n\t"		/* 4th 4 comparisons in 16 bits of r9  */		\
                "salq $16, %%r9		\n\t"			/*  */ \
                "salq $32, %%r10		\n\t"			/*  */ \
                "salq $48, %%r11		\n\t"			/*  */ \
                "orq	%%r11,%%r10		\n\t"		/* r8 now holds 16 comparisons in 64 bits */ \
                "orq	%%r9,%%r8		\n\t"		/* r8 now holds 8 comparisons in 32 bits */ \                                
                "orq	%%r10,%%r8		\n\t"		/* r8 now holds 12 comparisons in 48 bits */ \                                
                "movq   $0x2222222222222222,%%r9    \n\t" /* clear the bytemask results we don't care about */ \
                "andq   %%r9,%%r8    \n\t"                /* clear the bytemask results we don't care about */ \                
                "movl	%0,%%r11d		\n\t"		/* initialize count of set bits */ \
                "xorq	%%r10,%%r10		\n\t"		/* initialize bit scan offset */ \
                "1:			\n\t"					/* top of bit scan loop */ \
                "bsfq	%%r8,%%rcx		\n\t"		/* put least significant set bit index into rcx */ \
                "jz 2f	\n\t"						/* jump out if zero (no hits).  high percentage. */ \
                "addq	%%rcx,%%r10	\n\t"			/* add in the offset of this index */ \
                "movq   %%r10,%%r9 \n\t" \
                "sarq   $2,%%r9 \n\t"               /* translate to offset within bptr */ \
                "addl   %4,%%r9d \n\t"   \
                "movw	%%r9w, (%3, %%r11, 2) \n\t"		/* put the bit index into the output buffer */ \
                "shrq	%%cl,%%r8	\n\t"			/* shift the bit scan register up to the bit we just processed */ \
                "incl	%%r11d		\n\t"			/* increment the count of set bits */ \
                "incq	%%r10		\n\t"			/* increment the index */ \
                "shrq	$1, %%r8 \n\t"				/* clear the bit */ \
                "jmp 1b		\n\t"					/* loop if so */ \
                "2:		\n\t"						/*  */ \
                "movl	%%r11d, %0 \n\t"			/* return the count of set bits */ \
				: "+r"(result)						\
			    : "r"(bptr + j), "r"(mask), "r"(buffer), "r"(j)	\
                : "r8", "r9", "r10", "r11", "rcx", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "cc", "memory");	

    #endif

    #if defined(HAS_SSE2)

		#define SCAN_16X			\
			asm volatile (			\
				"movdqa (%2), %%xmm0	\n\t"	/*move mask into xmm0*/	\
				"movdqa (%1), %%xmm1	\n\t"	/*move 16 bptr locations into xmm regs*/	\
				"movdqa 16(%1), %%xmm2	\n\t"		\
				"pcmpeqw %%xmm0, %%xmm1	\n\t"	/*compare to mask*/	\
				"movdqa 32(%1), %%xmm3	\n\t"		\
				"pcmpeqw %%xmm0, %%xmm2	\n\t"		\
				"movdqa 48(%1), %%xmm4	\n\t"		\
				"pcmpeqw %%xmm0, %%xmm3	\n\t"		\
				"pcmpeqw %%xmm0, %%xmm4	\n\t"		\
				"por %%xmm1, %%xmm4		\n\t"	/*or the comparisons*/	\
				"por %%xmm2, %%xmm3		\n\t"		\
				"por %%xmm3, %%xmm4		\n\t"		\
				"pmovmskb %%xmm4, %0	\n\t"	/*if any are equal, this will be !0*/	\
				: "=r"(result)		\
				: "r"(bptr + j), "r"(mask)			\
				: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4");	

	#elif defined(HAS_MMX)
		#define SCAN_16X			\
			asm volatile (					/*this hasn't been tested yet...*/	\
				"movq (%2), %%mm0	\n\t"	/*move mask into xmm0*/	\
				"movq (%1), %%mm1	\n\t"	/*move 16 bptr locations into xmm regs*/	\
				"movq 8(%1), %%mm2	\n\t"		\
				"pcmpeqw %%mm0, %%mm1	\n\t"	/*compare to mask*/	\
				"movq 16(%1), %%mm3	\n\t"		\
				"pcmpeqw %%mm0, %%mm2	\n\t"		\
				"movq 24(%1), %%mm4	\n\t"		\
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
				"pcmpeqw %%mm0, %%mm1	\n\t"	/*compare to mask*/	\
				"movq 48(%1), %%mm3	\n\t"		\
				"pcmpeqw %%mm0, %%mm2	\n\t"		\
				"movq 56(%1), %%mm4	\n\t"		\
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
			result = 0xffff;	/*dont know what compiler this is. force the normal method*/
	#endif

#elif defined(MSC_ASM32A)

	#define SCAN_CLEAN ASM_M {emms};

	#if defined(HAS_SSE2)

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
				ASM_M pcmpeqw xmm1, xmm0	\
				ASM_M movdqa xmm3, XMMWORD PTR [edi + 32]	\
				ASM_M pcmpeqw xmm2, xmm0	\
				ASM_M movdqa xmm4, XMMWORD PTR [edi + 48]	\
				ASM_M pcmpeqw xmm3, xmm0	\
				ASM_M pcmpeqw xmm4, xmm0	\
				ASM_M por xmm4, xmm1		\
				ASM_M por xmm2, xmm3		\
				ASM_M por xmm4, xmm2		\
				ASM_M pmovmskb eax, xmm4	\
				ASM_M mov result, eax		}	\
			} while (0);



	#elif defined(HAS_MMX)

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
			ASM_M pcmpeqw mm1, mm0	\
			ASM_M movq mm3, QWORD PTR [edi + 16]	\
			ASM_M pcmpeqw mm2, mm0	\
			ASM_M movq mm4, QWORD PTR [edi + 24]	\
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
			ASM_M pcmpeqw mm1, mm0		\
			ASM_M movq mm3, QWORD PTR [edi + 48]	\
			ASM_M pcmpeqw mm2, mm0	\
			ASM_M movq mm4, QWORD PTR [edi + 56]	\
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
		#define SCAN_16X	\
			result = 1;	/*dont know what compiler this is. force the normal method*/
	#endif

#elif defined(_WIN64)

	#define SCAN_CLEAN /*nothing*/

	#if defined(HAS_SSE2)

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
			local_bptr = _mm_cmpeq_epi16(local_bptr, local_mask); \
			local_bptr3 = _mm_load_si128(bptr + j + 8); \
			local_bptr2 = _mm_cmpeq_epi16(local_bptr2, local_mask); \
			local_bptr4 = _mm_load_si128(bptr + j + 12); \
			local_bptr3 = _mm_cmpeq_epi16(local_bptr3, local_mask); \
			local_bptr4 = _mm_cmpeq_epi16(local_bptr4, local_mask); \
			local_bptr4 = _mm_or_si128(local_bptr4, local_bptr); \
			local_bptr2 = _mm_or_si128(local_bptr2, local_bptr3); \
			local_bptr4 = _mm_or_si128(local_bptr4, local_bptr2); \
			result = _mm_movemask_epi8(local_bptr4); \
			} while (0);

	#else

		#define SCAN_16X	\
			result = 1;	/* force the normal method*/

	#endif


#else	/* compiler not recognized*/

	#define SCAN_16X	\
		result = 1;	/* force the normal method */
	#define SCAN_CLEAN /*nothing*/

#endif

#ifdef USE_YAFU_TDIV
#define DIVIDE_ONE_PRIME \
	do	\
	{	\
		fb_offsets[++smooth_num] = i;	\
		zShortDiv32(tmp32, prime,tmp32);	\
	} while (zShortMod32(tmp32, prime) == 0);
#else
#define DIVIDE_ONE_PRIME \
	do \
	{						\
		fb_offsets[++smooth_num] = i;	\
		mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], prime); \
	} while (mpz_tdiv_ui(dconf->Qvals[report_num], prime) == 0); 

#define DIVIDE_RESIEVED_PRIME(j) \
            	while (mpz_tdiv_ui(dconf->Qvals[report_num], prime) == 0) \
                                    	{						\
		fb_offsets[++smooth_num] = j;	\
		mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], prime);		\
                                    	}
#endif

void tdiv_LP(uint32 report_num,  uint8 parity, uint32 bnum, 
	static_conf_t *sconf, dynamic_conf_t *dconf)
{
	int i,j,k;
	uint32 basebucket, prime;
	int smooth_num;
	uint32 *fb_offsets;
	uint32 *bptr;
    uint32 *fb = sconf->sieve_primes;
	uint32 block_loc;
	uint16 *mask = dconf->mask;
    uint16 buffer[32];

#ifdef USE_YAFU_TDIV
	z32 *tmp32 = &dconf->Qvals32[report_num];
#endif

#ifdef QS_TIMING
	gettimeofday(&qs_timing_start, NULL);
#endif

	fb_offsets = &dconf->fb_offsets[report_num][0];
	smooth_num = dconf->smooth_num[report_num];
	block_loc = dconf->reports[report_num];

	mask[0] = block_loc;
	mask[2] = block_loc;
	mask[4] = block_loc;
	mask[6] = block_loc;
    mask[8] = block_loc;
    mask[10] = block_loc;
    mask[12] = block_loc;
    mask[14] = block_loc;

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
    {
        basebucket = 0;
    }

	for (k=0; (uint32)k < dconf->buckets->num_slices; k++)
	{
		uint32 lpnum = *(dconf->buckets->num + bnum + basebucket);
        int r;
		uint32 fb_bound = *(dconf->buckets->fb_bounds + k);
		uint32 result = 0;

        //printf("tdiv_large for poly %d: slice %d bound = %u\n", 
        //    dconf->numB, k, fb_bound);

#if defined (_MSC_VER) || defined (FORCE_GENERIC)
        for (j=0; (uint32)j < (lpnum & (uint32)(~15)); j += 16)
        {
            SCAN_16X;

            if (result == 0)
                continue;

            //noticably faster to not put these in a loop!
            if (result & 0x2)
            {
                // could be j = 0, 4, 8, or 12
                if ((bptr[j] & 0x0000ffff) == block_loc)
                {
                    i = fb_bound + (bptr[j] >> 16);
                    prime = fb[i];
                    DIVIDE_ONE_PRIME;
                }
                if ((bptr[j + 4] & 0x0000ffff) == block_loc)
                {
                    i = fb_bound + (bptr[j + 4] >> 16);
                    prime = fb[i];
                    DIVIDE_ONE_PRIME;
                }
                if ((bptr[j + 8] & 0x0000ffff) == block_loc)
                {
                    i = fb_bound + (bptr[j + 8] >> 16);
                    prime = fb[i];
                    DIVIDE_ONE_PRIME;
                }
                if ((bptr[j + 12] & 0x0000ffff) == block_loc)
                {
                    i = fb_bound + (bptr[j + 12] >> 16);
                    prime = fb[i];
                    DIVIDE_ONE_PRIME;
                }
            }
            if (result & 0x20)
            {
                // could be j = 1, 5, 9, or 13
                if ((bptr[j + 1] & 0x0000ffff) == block_loc)
                {
                    i = fb_bound + (bptr[j + 1] >> 16);
                    prime = fb[i];
                    DIVIDE_ONE_PRIME;
                }
                if ((bptr[j + 5] & 0x0000ffff) == block_loc)
                {
                    i = fb_bound + (bptr[j + 5] >> 16);
                    prime = fb[i];
                    DIVIDE_ONE_PRIME;
                }
                if ((bptr[j + 9] & 0x0000ffff) == block_loc)
                {
                    i = fb_bound + (bptr[j + 9] >> 16);
                    prime = fb[i];
                    DIVIDE_ONE_PRIME;
                }
                if ((bptr[j + 13] & 0x0000ffff) == block_loc)
                {
                    i = fb_bound + (bptr[j + 13] >> 16);
                    prime = fb[i];
                    DIVIDE_ONE_PRIME;
                }
            }
            if (result & 0x200)
            {
                // could be j = 2, 6, 10, or 14
                if ((bptr[j + 2] & 0x0000ffff) == block_loc)
                {
                    i = fb_bound + (bptr[j + 2] >> 16);
                    prime = fb[i];
                    DIVIDE_ONE_PRIME;
                }
                if ((bptr[j + 6] & 0x0000ffff) == block_loc)
                {
                    i = fb_bound + (bptr[j + 6] >> 16);
                    prime = fb[i];
                    DIVIDE_ONE_PRIME;
                }
                if ((bptr[j + 10] & 0x0000ffff) == block_loc)
                {
                    i = fb_bound + (bptr[j + 10] >> 16);
                    prime = fb[i];
                    DIVIDE_ONE_PRIME;
                }
                if ((bptr[j + 14] & 0x0000ffff) == block_loc)
                {
                    i = fb_bound + (bptr[j + 14] >> 16);
                    prime = fb[i];
                    DIVIDE_ONE_PRIME;
                }
            }
            if (result & 0x2000)
            {
                // could be j= 3, 7, 11, or 15
                if ((bptr[j + 3] & 0x0000ffff) == block_loc)
                {
                    i = fb_bound + (bptr[j + 3] >> 16);
                    prime = fb[i];
                    DIVIDE_ONE_PRIME;
                }
                if ((bptr[j + 7] & 0x0000ffff) == block_loc)
                {
                    i = fb_bound + (bptr[j + 7] >> 16);
                    prime = fb[i];
                    DIVIDE_ONE_PRIME;
                }
                if ((bptr[j + 11] & 0x0000ffff) == block_loc)
                {
                    i = fb_bound + (bptr[j + 11] >> 16);
                    prime = fb[i];
                    DIVIDE_ONE_PRIME;
                }
                if ((bptr[j + 15] & 0x0000ffff) == block_loc)
                {
                    i = fb_bound + (bptr[j + 15] >> 16);
                    prime = fb[i];
                    DIVIDE_ONE_PRIME;
                }
            }
        }

        for (; (uint32)j < lpnum; j++)
        {
            if ((bptr[j] & 0x0000ffff) == block_loc)
            {
                i = fb_bound + (bptr[j] >> 16);
                prime = fb[i];
                //printf("block_loc = %u, bptr = %u, fb_bound = %u, fb_index = %u, prime = %u, Q mod prime = %u\n",
                //	block_loc, bptr[j].loc, fb_bound, bptr[j].fb_index, prime, zShortMod32(Q,prime));
                DIVIDE_ONE_PRIME;
            }
        }

#else



#if defined(USE_AVX2)
        CLEAN_AVX2;
#endif

        for (j = 0; (uint32)j < (lpnum & (uint32)(~15)); j += 16)
        {
            SCAN_16X_VEC;
        }

#if defined(USE_AVX2)
        CLEAN_AVX2;
#endif

        for (r = 0; r < result; r++)
        {
            i = fb_bound + (bptr[buffer[r]] >> 16);
            prime = fb[i];

            // Is this only necessary with AVX2, or with the new vector approach?
            if ((prime < 2) || (i >= sconf->factor_base->B))
            {
                dconf->lp_scan_failures++;
                continue;
            }

            DIVIDE_ONE_PRIME;

        }

		for (; (uint32)j < lpnum; j++)
		{
			if ((bptr[j] & 0x0000ffff) == block_loc)
			{
				i = fb_bound + (bptr[j] >> 16);
                prime = fb[i];

                if ((prime < 2) || (i >= sconf->factor_base->B))
                {
                    dconf->lp_scan_failures++;
                    continue;
                }

                DIVIDE_RESIEVED_PRIME(i);
			}
		}

#endif

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




