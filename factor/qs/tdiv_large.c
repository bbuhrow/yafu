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

#include "qs_impl.h"
#include "ytools.h"
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

#include <immintrin.h>

const uint32_t bitmask[16] = { 0x1, 0x2, 0x4, 0x8,
0x10, 0x20, 0x40, 0x80,
0x100, 0x200, 0x400, 0x800,
0x1000, 0x2000, 0x4000, 0x8000 };

#if (defined(GCC_ASM32X) || defined(GCC_ASM64X) || defined(__MINGW32__))
	
#ifdef _WIN32
#define ASM_ ASM_M
#else
#define ASM_ ASM_G
#endif



    #if defined(USE_AVX2)

        // these systems support SIMD 
        #define SCAN_CLEAN ASM_ volatile("emms");	


        #define SCAN_16X_VEC_b			\
			ASM_ volatile (			\
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
			ASM_ volatile (			\
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


	#elif defined(D_HAS_SSE2)

        // these systems support SIMD 
        #define SCAN_CLEAN ASM_ volatile("emms");	

		// top level sieve scanning with SSE2
        // the block_loc that we are looking for is in the
        // bottom half of each 32-bit value of bptr, so each
        // 128-bit xmm register holds 4 locations to search.
        // we load 4 registers worth, 16 locations total, 
        // perform the test on each, and OR all of the results together.

        #define SCAN_16X_VEC			\
			ASM_ volatile (			\
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


		#define SCAN_16X			\
			ASM_ volatile (			\
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
			ASM_ volatile (					/*this hasn't been tested yet...*/	\
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
			ASM_ volatile (			\
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

        #define SCAN_16X_VEC \
            for (i=0; i<16; i++) { \
                if ((bptr[j+i] & 0x0000ffff) == block_loc) { \
                    buffer[result++] = j+i; \
                } \
            }


        #define SCAN_CLEAN

		#define SCAN_16X	\
			result = 0xffff;	/*dont know what compiler this is. force the normal method*/
	#endif

#elif defined(_WIN64)

	#define SCAN_CLEAN /*nothing*/

	#if defined(D_HAS_SSE2)

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
			result = 0xffff;	/* force the normal method*/

	#endif


#else	/* compiler not recognized*/

	#define SCAN_16X	\
		result = 1;	/* force the normal method */
	#define SCAN_CLEAN /*nothing*/

#endif


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

#define DIVIDE_VLP_PRIME(j) \
    while (mpz_tdiv_ui(dconf->Qvals[report_num], prime) == 0) \
    {						\
        count++; \
        /*gmp_printf("vlp prime %u divides %Zu\n", prime, dconf->Qvals[report_num]); */ \
		fb_offsets[++smooth_num] = j;	\
		mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], prime);		\
    }

#if defined(USE_AVX2)
void tdiv_LP_avx2(uint32_t report_num,  uint8_t parity, uint32_t bnum, 
	static_conf_t *sconf, dynamic_conf_t *dconf)
{
	int i,j,k;
	uint32_t basebucket, prime;
	int smooth_num;
	uint32_t *fb_offsets;
	uint32_t *bptr;
    uint32_t *fb = sconf->sieve_primes;
	uint32_t block_loc;
	uint16_t *mask = dconf->mask;
    uint16_t buffer[32];    

#if defined(USE_BATCHPOLY) || defined(USE_BATCHPOLY_X2)
    int poly_offset = (dconf->numB % dconf->poly_batchsize) - 2;
    int pnum;

    if (dconf->numB == 1)
    {
        poly_offset = 0;
    }
    else if (poly_offset < 0)
    {
        poly_offset += dconf->poly_batchsize;
    }
    pnum = poly_offset;
    poly_offset = poly_offset * 2 * sconf->num_blocks * dconf->buckets->alloc_slices;

    //printf("begin large_tdiv on side %d with poly %d for location %u\n", 
    //    parity, pnum, dconf->reports[report_num]);

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
#if defined(USE_BATCHPOLY) || defined(USE_BATCHPOLY_X2)
    bptr = dconf->buckets->list + (bnum << BUCKET_BITS) + poly_offset * BUCKET_ALLOC;
#else
	bptr = dconf->buckets->list + (bnum << BUCKET_BITS);
#endif

    if (parity)
    {
        bptr += (sconf->num_blocks << BUCKET_BITS);
        basebucket = sconf->num_blocks;
    }
    else
    {
        basebucket = 0;
    }

#if defined(USE_BATCHPOLY_X2)
    for (k = 0; (uint32_t)k < dconf->buckets->num_slices_batch[pnum]; k++)
#else
	for (k=0; (uint32_t)k < dconf->buckets->num_slices; k++)
#endif
	{
        
#if defined(USE_BATCHPOLY) || defined(USE_BATCHPOLY_X2)
        uint32_t lpnum = *(dconf->buckets->num + bnum + basebucket + poly_offset);
#else
        uint32_t lpnum = *(dconf->buckets->num + bnum + basebucket);
        //uint32_t lpnum = bptr[0];
#endif

        int r, q;
#if defined(USE_BATCHPOLY_X2)
        uint32_t fb_bound = *(dconf->buckets->fb_bounds + k + pnum * dconf->buckets->alloc_slices);
#else
		uint32_t fb_bound = *(dconf->buckets->fb_bounds + k);
#endif
		uint32_t result = 0;

#if defined (_MSC_VER) && !defined(__INTEL_COMPILER) && !defined(__clang__)
        for (j = 0; (uint32_t)j < (lpnum & (uint32_t)(~15)); j += 16)
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

        // leftover bucket elements to check after doing 16x at a time
        for (; (uint32_t)j < lpnum; j++)
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

        CLEAN_AVX2;

        for (j = 0; (uint32_t)j < (lpnum & (uint32_t)(~15)); j += 16)
        {
            SCAN_16X_VEC;
        }

        CLEAN_AVX2;

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

		for (; (uint32_t)j < lpnum; j++)
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

	dconf->smooth_num[report_num] = smooth_num;

	return;
}

void tdiv_LP_sse2(uint32_t report_num, uint8_t parity, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf)
{
    return;
}

#elif defined(D_HAS_SSE2)
void tdiv_LP_sse2(uint32_t report_num, uint8_t parity, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf)
{
    int i, j, k;
    uint32_t basebucket, prime;
    int smooth_num;
    uint32_t* fb_offsets;
    uint32_t* bptr;
    uint32_t* fb = sconf->sieve_primes;
    uint32_t block_loc;
    uint16_t* mask = dconf->mask;
    uint16_t buffer[32];

#if defined(USE_BATCHPOLY) || defined(USE_BATCHPOLY_X2)
    int poly_offset = (dconf->numB % dconf->poly_batchsize) - 2;
    int pnum;

    if (dconf->numB == 1)
    {
        poly_offset = 0;
    }
    else if (poly_offset < 0)
    {
        poly_offset += dconf->poly_batchsize;
    }
    pnum = poly_offset;
    poly_offset = poly_offset * 2 * sconf->num_blocks * dconf->buckets->alloc_slices;

    //printf("begin large_tdiv on side %d with poly %d for location %u\n", 
    //    parity, pnum, dconf->reports[report_num]);

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
#if defined(USE_BATCHPOLY) || defined(USE_BATCHPOLY_X2)
    bptr = dconf->buckets->list + (bnum << BUCKET_BITS) + poly_offset * BUCKET_ALLOC;
#else
    bptr = dconf->buckets->list + (bnum << BUCKET_BITS);
#endif

    if (parity)
    {
        bptr += (sconf->num_blocks << BUCKET_BITS);
        basebucket = sconf->num_blocks;
    }
    else
    {
        basebucket = 0;
    }

#if defined(USE_BATCHPOLY_X2)
    for (k = 0; (uint32_t)k < dconf->buckets->num_slices_batch[pnum]; k++)
#else
    for (k = 0; (uint32_t)k < dconf->buckets->num_slices; k++)
#endif
    {

#if defined(USE_BATCHPOLY) || defined(USE_BATCHPOLY_X2)
        uint32_t lpnum = *(dconf->buckets->num + bnum + basebucket + poly_offset);
#else
        uint32_t lpnum = *(dconf->buckets->num + bnum + basebucket);
        //uint32_t lpnum = bptr[0];
#endif

        int r, q;
#if defined(USE_BATCHPOLY_X2)
        uint32_t fb_bound = *(dconf->buckets->fb_bounds + k + pnum * dconf->buckets->alloc_slices);
#else
        uint32_t fb_bound = *(dconf->buckets->fb_bounds + k);
#endif
        uint32_t result = 0;

        for (j = 0; (uint32_t)j < (lpnum & (uint32_t)(~15)); j += 16)
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

        // leftover bucket elements to check after doing 16x at a time
        for (; (uint32_t)j < lpnum; j++)
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


        //point to the next slice of primes
        bptr += (sconf->num_blocks << (BUCKET_BITS + 1));
        basebucket += (sconf->num_blocks << 1);
    }

    SCAN_CLEAN;

    dconf->smooth_num[report_num] = smooth_num;

    return;
}
#endif

#if defined(USE_AVX512F)
void tdiv_LP_avx512(uint32_t report_num, uint8_t parity, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf)
{
    int i, j, k;
    uint32_t basebucket, prime;
    int smooth_num;
    uint32_t* fb_offsets;
    uint32_t* bptr;
    uint32_t* fb = sconf->sieve_primes;
    uint32_t block_loc;
    uint16_t* mask = dconf->mask;
    uint16_t buffer[32];
    __m512i vmask, vblock;


#if defined(USE_BATCHPOLY) || defined(USE_BATCHPOLY_X2)
    int poly_offset = (dconf->numB % dconf->poly_batchsize) - 2;
    int pnum;

    if (dconf->numB == 1)
    {
        poly_offset = 0;
    }
    else if (poly_offset < 0)
    {
        poly_offset += dconf->poly_batchsize;
    }
    pnum = poly_offset;
    poly_offset = poly_offset * 2 * sconf->num_blocks * dconf->buckets->alloc_slices;

#ifdef DEBUGPRINT_BATCHPOLY
    printf("begin large_tdiv on side %d with poly %d for location %u... ", 
        parity, pnum, dconf->reports[report_num]);
#endif

#endif

    fb_offsets = &dconf->fb_offsets[report_num][0];
    smooth_num = dconf->smooth_num[report_num];
    block_loc = dconf->reports[report_num];

    // 16 copies of the 16-bit block_loc in the lower half of
    // each of the 32-bit vector elements.
    vblock = _mm512_set1_epi32(block_loc);
    vmask = _mm512_set1_epi32(0x0000ffff);

    //primes bigger than med_B are bucket sieved, so we need
    //only search through the bucket and see if any locations match the
    //current block index.
#if defined(USE_BATCHPOLY) || defined(USE_BATCHPOLY_X2)
    bptr = dconf->buckets->list + (bnum << BUCKET_BITS) + poly_offset * BUCKET_ALLOC;
#else
    bptr = dconf->buckets->list + (bnum << BUCKET_BITS);
#endif

    if (parity)
    {
        bptr += (sconf->num_blocks << BUCKET_BITS);
        basebucket = sconf->num_blocks;
    }
    else
    {
        basebucket = 0;
    }

#if defined(USE_BATCHPOLY_X2)
    for (k = 0; (uint32_t)k < dconf->buckets->num_slices_batch[pnum]; k++)
#else
    for (k = 0; (uint32_t)k < dconf->buckets->num_slices; k++)
#endif
    {

#if defined(USE_BATCHPOLY) || defined(USE_BATCHPOLY_X2)
        uint32_t lpnum = *(dconf->buckets->num + bnum + basebucket + poly_offset);
#else
        uint32_t lpnum = *(dconf->buckets->num + bnum + basebucket);
        //uint32_t lpnum = bptr[0];
#endif

        int r, q;
#if defined(USE_BATCHPOLY_X2)
        uint32_t fb_bound = *(dconf->buckets->fb_bounds + k + pnum * dconf->buckets->alloc_slices);
#else
        uint32_t fb_bound = *(dconf->buckets->fb_bounds + k);
#endif
        uint32_t result = 0;


        //#if defined(USE_BATCHPOLY_X2)
        //        printf("lp tdiv: checking %d primes from slice %d of %d for poly %d, bucket/block %d\n",
        //            lpnum, k, dconf->buckets->num_slices_batch[pnum], pnum, bnum);
        //#else
        //        printf("lp tdiv: checking %d primes from slice %d of %d for poly %d, bucket/block %d\n",
        //            lpnum, k, dconf->buckets->num_slices, 0, bnum);
        //#endif

        for (j = 0; (uint32_t)j < (lpnum & (uint32_t)(~15)); j += 16)
        {
            uint32_t idx;
            __m512i velements = _mm512_load_epi32(bptr + j);
            velements = _mm512_and_epi32(velements, vmask);
            result = _mm512_cmp_epu32_mask(velements, vblock, _MM_CMPINT_EQ);

            while (result > 0)
            {
                idx = _trail_zcnt(result);
                i = fb_bound + (bptr[j + idx] >> 16);
                prime = fb[i];
                DIVIDE_RESIEVED_PRIME(i);

                result = _reset_lsb(result);
            }

        }

        //point to the next slice of primes
        bptr += (sconf->num_blocks << (BUCKET_BITS + 1));
        basebucket += (sconf->num_blocks << 1);
    }

    SCAN_CLEAN;

    dconf->smooth_num[report_num] = smooth_num;

#ifdef DEBUGPRINT_BATCHPOLY
    printf("complete.\n"); fflush(stdout);
#endif

    return;
}


#endif

#ifdef USE_POLY_BUCKET_SS

void tdiv_SS(uint32_t report_num, uint8_t parity, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf)
{
    // bucket-sorted by polynomial
    int i;
    int smooth_num;
    uint32_t* fb_offsets;
    uint32_t block_loc;
    struct timeval start1, stop1;

    fb_offsets = &dconf->fb_offsets[report_num][0];
    smooth_num = dconf->smooth_num[report_num];
    block_loc = dconf->reports[report_num];
    int pidx = dconf->polymap[dconf->numB];
    int bucketalloc = dconf->ss_slices_p[0].alloc;

#ifdef SS_POLY_BUCKET_SMALL_GROUPS
    // if the mapped binary-encoded poly isn't in this block of 
    // poly buckets then just skip large prime sieving.
    //if ((pidx < dconf->ss_slices_p[0].curr_poly_idx) ||
    //    (pidx >= (dconf->ss_slices_p[0].curr_poly_idx + (1 << dconf->ss_set1.size))))
    //    return;

    //printf("mapping b-index %d to bucket %d\n", dconf->numB, pidx);

    //printf("commencing tdiv_ss on side %d on pidx %u (set2 instance %d), sizes %d,%d\n",
    //    parity, pidx, polymask / (1 << dconf->ss_set2.size),
    //    (1 << dconf->ss_set1.size), (1 << dconf->ss_set2.size));
#endif

    block_loc += bnum * 32768;
    
    uint32_t pid_offset = SS_SIGN_BIT + 1;

#ifdef USE_POLY_BUCKET_PN_COMBINED_VARIATION
    uint32_t signbit = SS_MAX_ROOT;
    uint32_t rootmask = SS_ROOT_MASK;
    
#ifdef USE_AVX512F
    __m512i vloc = _mm512_set1_epi32(block_loc);
    __m512i vrootmask = _mm512_set1_epi32(rootmask);
    __m512i vposmask = _mm512_set1_epi32(signbit);

#elif defined(USE_AVX2)
    __m256i vz = _mm256_setzero_si256();
    __m256i vloc = _mm256_set1_epi32(block_loc);
    __m256i vrootmask = _mm256_set1_epi32(rootmask);
    __m256i vposmask = _mm256_set1_epi32(signbit);

#endif

#else

#ifdef USE_AVX512F
    __m512i vrootmask = _mm512_set1_epi32(SS_ROOT_MASK);
#elif defined(USE_AVX2)
    __m256i vrootmask = _mm256_set1_epi32(SS_ROOT_MASK);
#endif

#endif

    if (parity == 0)
    {
        for (i = 0; i < dconf->num_ss_slices; i++)
        {
            uint32_t* bucketelements = dconf->ss_slices_p[i].elements + pidx * bucketalloc;
            uint32_t root;
            uint32_t fboffset = dconf->ss_slices_p[i].fboffset;

            int k = 0;

#ifdef USE_AVX512F
            for (k = 0; k < ((int)dconf->ss_slices_p[i].size[pidx] - 16); k += 16)
            {
                __m512i vr = _mm512_loadu_epi32(bucketelements + k);
#ifdef USE_POLY_BUCKET_PN_COMBINED_VARIATION
                __mmask16 mpos = ~_mm512_test_epi32_mask(vr, vposmask);
#else
                __mmask16 mpos = 0xffff;
#endif
                vr = _mm512_and_epi32(vr, vrootmask);

                mpos = _mm512_mask_cmpeq_epi32_mask(mpos, vr, vloc);

                while (mpos > 0)
                {
                    int idx = _trail_zcnt(mpos);

                    uint32_t pid = fboffset + (((bucketelements[k + idx]) >> pid_offset));
                    uint32_t prime = sconf->factor_base->list->prime[pid];

                    // these checks don't really impact performance at all, and
                    // they are useful for detecting bugs and problem configurations.
                    if ((mpz_tdiv_ui(dconf->Qvals[report_num], prime) != 0) && (dconf->numB > 1))
                    {
                        printf("tdiv invalid root %u for loc %u in slice %u, side %u, poly %u (%u) "
                            "pid = %u, fboffset = %u\n",
                            bucketelements[k + idx] & rootmask, block_loc, i, parity,
                            dconf->numB, pidx, pid, fboffset);
                        exit(1);
                    }

                    while (mpz_tdiv_ui(dconf->Qvals[report_num], prime) == 0)
                    {
                        fb_offsets[++smooth_num] = pid;
                        mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num],
                            prime);
                    }

                    mpos = _reset_lsb(mpos);
                }
            }

#elif defined(USE_AVX2)

            for (k = 0; k < ((int)dconf->ss_slices_p[i].size[pidx] - 8); k += 8)
            {
                __m256i vr = _mm256_load_si256((__m256i *)(&bucketelements[k]));
                //__mmask16 mpos = ~_mm512_test_epi32_mask(vr, vposmask);
                __m256i cmp = _mm256_and_si256(vr, vposmask);
                cmp = _mm256_cmpeq_epi32(cmp, vz);
                uint32_t mpos = _mm256_movemask_epi8(cmp) & 0x88888888;;

                //mpos = _mm512_mask_cmpeq_epi32_mask(mpos, vr, vloc);
                vr = _mm256_and_si256(vr, vrootmask);
                cmp = _mm256_cmpeq_epi32(vr, vloc);
                mpos &= _mm256_movemask_epi8(cmp);
                

                while (mpos > 0)
                {
                    int idx = _trail_zcnt(mpos) >> 2;

                    uint32_t pid = fboffset + (((bucketelements[k + idx]) >> pid_offset));
                    uint32_t prime = sconf->factor_base->list->prime[pid];

                    if ((mpz_tdiv_ui(dconf->Qvals[report_num], prime) != 0) && (dconf->numB > 1))
                    {
                        printf("tdiv invalid root %u for loc %u in slice %u, side %u, poly %u (%u) "
                            "pid = %u, fboffset = %u\n",
                            bucketelements[k + idx] & rootmask, block_loc, i, parity,
                            dconf->numB, pidx, pid, fboffset);
                        exit(1);
                    }

                    while (mpz_tdiv_ui(dconf->Qvals[report_num], prime) == 0)
                    {
                        fb_offsets[++smooth_num] = pid;
                        mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num],
                            prime);
                    }

                    mpos = _reset_lsb(mpos);
                }
            }
            
#endif

#if defined(USE_AVX512F) || (USE_AVX2)
            for ( ; k < dconf->ss_slices_p[i].size[pidx]; k++)
#else
            for (k = 0; k < dconf->ss_slices_p[i].size[pidx]; k++)
#endif
            {
#ifdef USE_POLY_BUCKET_PN_COMBINED_VARIATION
                root = (bucketelements[k] & rootmask);
                if ((bucketelements[k] & signbit)) continue;
#else
                root = (bucketelements[k] & 0x3ffff);
#endif

                if (block_loc == root)
                {
                    uint32_t pid = fboffset + (((bucketelements[k]) >> pid_offset));
                    uint32_t prime = sconf->factor_base->list->prime[pid];

                    if ((mpz_tdiv_ui(dconf->Qvals[report_num], prime) != 0) && (dconf->numB > 1))
                    {
                        printf("tdiv invalid root %u in slice %u, side %u, poly %u "
                            "pid = %u, fboffset = %u\n",
                            root, i, parity, dconf->numB, pid, fboffset);
                        exit(1);
                    }

                    while (mpz_tdiv_ui(dconf->Qvals[report_num], prime) == 0)
                    {
                        fb_offsets[++smooth_num] = pid;
                        mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num],
                            prime);
                    }

                }
            }
        }
    }
    else
    {
        for (i = 0; i < dconf->num_ss_slices; i++)
        {

#ifdef USE_POLY_BUCKET_PN_COMBINED_VARIATION


            uint32_t* bucketelements = dconf->ss_slices_p[i].elements + pidx * bucketalloc;
            uint32_t root;
            uint32_t fboffset = dconf->ss_slices_p[i].fboffset;

            int k = 0;
#ifdef USE_AVX512F
            for (k = 0; k < ((int)dconf->ss_slices_p[i].size[pidx] - 16); k += 16)
            {
                __m512i vr = _mm512_loadu_epi32(bucketelements + k);
                __mmask16 mpos = _mm512_test_epi32_mask(vr, vposmask);
                vr = _mm512_and_epi32(vr, vrootmask);

                mpos = _mm512_mask_cmpeq_epi32_mask(mpos, vr, vloc);

                while (mpos > 0)
                {
                    int idx = _trail_zcnt(mpos);

                    uint32_t pid = fboffset + (((bucketelements[k + idx]) >> pid_offset));
                    uint32_t prime = sconf->factor_base->list->prime[pid];

                    if ((mpz_tdiv_ui(dconf->Qvals[report_num], prime) != 0) && (dconf->numB > 1))
                    {
                        printf("tdiv invalid root %u for loc %u in slice %u, side %u, poly %u (%u) "
                            "pid = %u, fboffset = %u\n",
                            bucketelements[k + idx] & rootmask, block_loc, i, parity,
                            dconf->numB, pidx, pid, fboffset);
                        exit(1);
                    }

                    while (mpz_tdiv_ui(dconf->Qvals[report_num], prime) == 0)
                    {
                        fb_offsets[++smooth_num] = pid;
                        mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num],
                            prime);
                    }

                    mpos = _reset_lsb(mpos);
                }
            }

#elif defined(USE_AVX2)

            for (k = 0; k < ((int)dconf->ss_slices_p[i].size[pidx] - 8); k += 8)
            {
                __m256i vr = _mm256_load_si256((__m256i*)(&bucketelements[k]));
                //__mmask16 mpos = ~_mm512_test_epi32_mask(vr, vposmask);
                __m256i cmp = _mm256_and_si256(vr, vposmask);
                cmp = _mm256_cmpgt_epi32(cmp, vz);
                uint32_t mpos = _mm256_movemask_epi8(cmp) & 0x88888888;

                //mpos = _mm512_mask_cmpeq_epi32_mask(mpos, vr, vloc);
                vr = _mm256_and_si256(vr, vrootmask);
                cmp = _mm256_cmpeq_epi32(vr, vloc);
                mpos &= _mm256_movemask_epi8(cmp);


                while (mpos > 0)
                {
                    int idx = _trail_zcnt(mpos) >> 2;

                    uint32_t pid = fboffset + (((bucketelements[k + idx]) >> pid_offset));
                    uint32_t prime = sconf->factor_base->list->prime[pid];

                    if ((mpz_tdiv_ui(dconf->Qvals[report_num], prime) != 0) && (dconf->numB > 1))
                    {
                        printf("tdiv invalid root %u for loc %u in slice %u, side %u, poly %u (%u) "
                            "pid = %u, fboffset = %u\n",
                            bucketelements[k + idx] & rootmask, block_loc, i, parity,
                            dconf->numB, pidx, pid, fboffset);
                        exit(1);
                    }

                    while (mpz_tdiv_ui(dconf->Qvals[report_num], prime) == 0)
                    {
                        fb_offsets[++smooth_num] = pid;
                        mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num],
                            prime);
                    }

                    mpos = _reset_lsb(mpos);
                }
            }
#endif

#if defined(USE_AVX512F) || (USE_AVX2)
            for (; k < dconf->ss_slices_p[i].size[pidx]; k++)
#else
            for (k = 0; k < dconf->ss_slices_p[i].size[pidx]; k++)
#endif
            {
                root = (bucketelements[k] & rootmask);
                if ((bucketelements[k] & signbit) == 0) continue;

                if (block_loc == root)
                {
                    uint32_t pid = fboffset + (((bucketelements[k]) >> pid_offset));
                    uint32_t prime = sconf->factor_base->list->prime[pid];

                    if ((mpz_tdiv_ui(dconf->Qvals[report_num], prime) != 0) && (dconf->numB > 1))
                    {
                        printf("tdiv invalid root %u in slice %u, side %u, poly %u "
                            "pid = %u, fboffset = %u\n",
                            root, i, parity, dconf->numB, pid, fboffset);
                        exit(1);
                    }

                    while (mpz_tdiv_ui(dconf->Qvals[report_num], prime) == 0)
                    {
                        fb_offsets[++smooth_num] = pid;
                        mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num],
                            prime);
                    }

                }
            }
#else
            uint32_t* bucketelements = dconf->ss_slices_n[i].elements + pidx * bucketalloc;
            uint32_t root;
            uint32_t fboffset = dconf->ss_slices_n[i].fboffset;

            int k;
#ifdef USE_AVX512F
            for (k = 0; k < dconf->ss_slices_n[i].size[pidx] - 16; k += 16)
            {
                __m512i vr = _mm512_loadu_epi32(bucketelements + k);
                __mmask16 mpos = 0xffff; // _mm512_test_epi32_mask(vr, vposmask);
                vr = _mm512_and_epi32(vr, vrootmask);

                mpos = _mm512_mask_cmpeq_epi32_mask(mpos, vr, vloc);

                while (mpos > 0)
                {
                    int idx = _trail_zcnt(mpos);

                    uint32_t pid = fboffset + (((bucketelements[k + idx]) >> 18));
                    uint32_t prime = sconf->factor_base->list->prime[pid];

                    if ((mpz_tdiv_ui(dconf->Qvals[report_num], prime) != 0) && (dconf->numB > 1))
                    {
                        printf("tdiv invalid root %u for loc %u in slice %u, side %u, poly %u "
                            "pid = %u, fboffset = %u\n",
                            bucketelements[k + idx] & 0x1ffff, block_loc, i, parity, dconf->numB, pid, fboffset);
                        exit(1);
                    }

                    while (mpz_tdiv_ui(dconf->Qvals[report_num], prime) == 0)
                    {
                        fb_offsets[++smooth_num] = pid;
                        mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num],
                            prime);
                    }

                    mpos = _reset_lsb(mpos);
                }
            }

            for ( ; k < dconf->ss_slices_n[i].size[pidx]; k++)
#else
            for (k = 0; k < dconf->ss_slices_n[i].size[pidx]; k++)
#endif
            {
                root = (bucketelements[k] & 0x3ffff);

                if (block_loc == root)
                {
                    uint32_t pid = fboffset + (((bucketelements[k]) >> pid_offset));
                    uint32_t prime = sconf->factor_base->list->prime[pid];

                    if ((mpz_tdiv_ui(dconf->Qvals[report_num], prime) != 0) && (dconf->numB > 1))
                    {
                        printf("tdiv invalid root %u in slice %u, side %u, poly %u "
                            "pid = %u, fboffset = %u\n",
                            root, i, parity, dconf->numB, pid, fboffset);
                        exit(1);
                    }

                    while (mpz_tdiv_ui(dconf->Qvals[report_num], prime) == 0)
                    {
                        fb_offsets[++smooth_num] = pid;
                        mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num],
                            prime);
                    }

                }
            }

#endif
        }

    }

    dconf->smooth_num[report_num] = smooth_num;

    return;
}

#endif

