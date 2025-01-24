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

this file contains code implementing 1)


*/

#if defined(GCC_ASM64X) || defined(__MINGW64__)
	#define SCAN_CLEAN asm volatile("emms");	

#if defined(USE_AVX2)


    #define SIEVE_SCAN_32_VEC					\
		asm volatile (							\
			"vmovdqa (%1), %%ymm0   \n\t"		\
			"vpmovmskb %%ymm0, %%r11   \n\t"		/* output results to 64 bit register */		\
            "movq %%r11, %%r8   \n\t"		\
			                                    /* r8 now holds 64 byte mask results, in order, from sieveblock */ \
			"xorq	%%r11,%%r11		\n\t"		/* initialize count of set bits */ \
			"xorq	%%r10,%%r10		\n\t"		/* initialize bit scan offset */ \
			"1:			\n\t"					/* top of bit scan loop */ \
			"bsfq	%%r8,%%rcx		\n\t"		/* put least significant set bit index into rcx */ \
            "jz 2f	\n\t"						/* jump out if zero (no hits).  high percentage. */ \
			"addq	%%rcx,%%r10	\n\t"			/* add in the offset of this index */ \
			"movb	%%r10b, (%2, %%r11, 1) \n\t"		/* put the bit index into the output buffer */ \
			"shrq	%%cl,%%r8	\n\t"			/* shift the bit scan register up to the bit we just processed */ \
			"incq	%%r11		\n\t"			/* increment the count of set bits */ \
            "incq	%%r10		\n\t"			/* increment the index */ \
			"shrq	$1, %%r8 \n\t"				/* clear the bit */ \
			"jmp 1b		\n\t"					/* loop if so */ \
			"2:		\n\t"						/*  */ \
			"movl	%%r11d, %0 \n\t"			/* return the count of set bits */ \
			: "=r"(result)						\
			: "r"(sieveblock + j), "r"(buffer)	\
			: "xmm0", "r8", "r9", "r10", "r11", "rcx", "cc", "memory");

    #define SIEVE_SCAN_64_VEC					\
		asm volatile (							\
			"vmovdqa (%1), %%ymm0   \n\t"		\
			"vpor 32(%1), %%ymm0, %%ymm0    \n\t"		\
			"vpmovmskb %%ymm0, %%r11   \n\t"	/* output results to 64 bit register */		\
			"testq %%r11, %%r11 \n\t"			/* AND, and set ZF */ \
			"jz 2f	\n\t"						/* jump out if zero (no hits).  high percentage. */ \
			"vmovdqa (%1), %%ymm0   \n\t"		/* else, we had hits, move sections of sieveblock back in */ \
			"vmovdqa 32(%1), %%ymm1   \n\t"		/* extract high bit masks from each byte */ \
			"vpmovmskb %%ymm1, %%r9d   \n\t"		/*  */		\
			"salq $32, %%r9		\n\t"			/*  */ \
			"vpmovmskb %%ymm0, %%r8d   \n\t"		/*  */		\
			"orq	%%r9,%%r8		\n\t"		/* r8 now holds 64 byte mask results, in order, from sieveblock */ \
			"xorq	%%r11,%%r11		\n\t"		/* initialize count of set bits */ \
			"xorq	%%r10,%%r10		\n\t"		/* initialize bit scan offset */ \
			"1:			\n\t"					/* top of bit scan loop */ \
			"bsfq	%%r8,%%rcx		\n\t"		/* put least significant set bit index into rcx */ \
			"addq	%%rcx,%%r10	\n\t"			/* add in the offset of this index */ \
			"movb	%%r10b, (%2, %%r11, 1) \n\t"		/* put the bit index into the output buffer */ \
			"shrq	%%cl,%%r8	\n\t"			/* shift the bit scan register up to the bit we just processed */ \
			"incq	%%r11		\n\t"			/* increment the count of set bits */ \
            "incq	%%r10		\n\t"			/* increment the index */ \
			"shrq	$1, %%r8 \n\t"				/* clear the bit */ \
			"testq	%%r8,%%r8	\n\t"			/* check if there are any more set bits */ \
			"jnz 1b		\n\t"					/* loop if so */ \
			"2:		\n\t"						/*  */ \
			"movl	%%r11d, %0 \n\t"			/* return the count of set bits */ \
			: "=r"(result)						\
			: "r"(sieveblock + j), "r"(buffer)	\
			: "xmm0", "xmm1", "r8", "r9", "r10", "r11", "rcx", "cc", "memory");


#else

	//top level sieve scanning with SSE2
	#define SIEVE_SCAN_32_VEC					\
		asm volatile (							\
			"movdqa (%1), %%xmm0   \n\t"		\
			"por 16(%1), %%xmm0    \n\t"		\
			"pmovmskb %%xmm0, %%r11   \n\t"		/* output results to 64 bit register */		\
			"testq %%r11, %%r11 \n\t"			/* AND, and set ZF */ \
			"jz 2f	\n\t"						/* jump out if zero (no hits).  high percentage. */ \
			"movdqa (%1), %%xmm0   \n\t"		/* else, we had hits, move sections of sieveblock back in */ \
			"movdqa 16(%1), %%xmm1   \n\t"		/* there are 16 bytes in each section */ \
			"pmovmskb %%xmm1, %%r9d   \n\t"		/*  */		\
			"salq $16, %%r9		\n\t"			/*  */ \
			"pmovmskb %%xmm0, %%r8d   \n\t"		/*  */		\
			"orq	%%r9,%%r8		\n\t"		/* r8 now holds 64 byte mask results, in order, from sieveblock */ \
			"xorq	%%r11,%%r11		\n\t"		/* initialize count of set bits */ \
			"xorq	%%r10,%%r10		\n\t"		/* initialize bit scan offset */ \
			"1:			\n\t"					/* top of bit scan loop */ \
			"bsfq	%%r8,%%rcx		\n\t"		/* put least significant set bit index into rcx */ \
            "jz 2f	\n\t"						/* jump out if zero (no hits).  high percentage. */ \
			"addq	%%rcx,%%r10	\n\t"			/* add in the offset of this index */ \
			"movb	%%r10b, (%2, %%r11, 1) \n\t"		/* put the bit index into the output buffer */ \
			"shrq	%%cl,%%r8	\n\t"			/* shift the bit scan register up to the bit we just processed */ \
			"incq	%%r11		\n\t"			/* increment the count of set bits */ \
            "incq	%%r10		\n\t"			/* increment the index */ \
			"shrq	$1, %%r8 \n\t"				/* clear the bit */ \
			"jmp 1b		\n\t"					/* loop if so */ \
			"2:		\n\t"						/*  */ \
			"movl	%%r11d, %0 \n\t"			/* return the count of set bits */ \
			: "=r"(result)						\
			: "r"(sieveblock + j), "r"(buffer)	\
			: "xmm0", "xmm1", "xmm2", "xmm3", "r8", "r9", "r10", "r11", "rcx", "cc", "memory");

	#define SIEVE_SCAN_64_VEC					\
		asm volatile (							\
			"movdqa (%1), %%xmm0   \n\t"		\
			"por 16(%1), %%xmm0    \n\t"		\
			"por 32(%1), %%xmm0    \n\t"		\
			"por 48(%1), %%xmm0    \n\t"		\
			"pmovmskb %%xmm0, %%r11   \n\t"		/* output results to 64 bit register */		\
			"testq %%r11, %%r11 \n\t"			/* AND, and set ZF */ \
			"jz 2f	\n\t"						/* jump out if zero (no hits).  high percentage. */ \
			"movdqa (%1), %%xmm0   \n\t"		/* else, we had hits, move sections of sieveblock back in */ \
			"movdqa 16(%1), %%xmm1   \n\t"		/* there are 16 bytes in each section */ \
			"movdqa 32(%1), %%xmm2   \n\t"		/* extract high bit masks from each byte */ \
			"movdqa 48(%1), %%xmm3   \n\t"		/* and combine into one 64 bit register */ \
			"pmovmskb %%xmm1, %%r9d   \n\t"		/*  */		\
			"pmovmskb %%xmm3, %%r11d   \n\t"	/*  */		\
			"salq $16, %%r9		\n\t"			/*  */ \
			"pmovmskb %%xmm2, %%r10d   \n\t"	/*  */		\
			"salq $48, %%r11		\n\t"		/*  */ \
			"pmovmskb %%xmm0, %%r8d   \n\t"		/*  */		\
			"salq $32, %%r10		\n\t"		/*  */ \
			"orq	%%r11,%%r9		\n\t"		/*  */ \
			"orq	%%r10,%%r8		\n\t"		/*  */ \
			"xorq	%%r11,%%r11		\n\t"		/* initialize count of set bits */ \
			"orq	%%r9,%%r8		\n\t"		/* r8 now holds 64 byte mask results, in order, from sieveblock */ \
			"xorq	%%r10,%%r10		\n\t"		/* initialize bit scan offset */ \
			"1:			\n\t"					/* top of bit scan loop */ \
			"bsfq	%%r8,%%rcx		\n\t"		/* put least significant set bit index into rcx */ \
            "jz 2f	\n\t"						/* jump out if zero (no hits).  high percentage. */ \
			"addq	%%rcx,%%r10	\n\t"			/* add in the offset of this index */ \
			"movb	%%r10b, (%2, %%r11, 1) \n\t"		/* put the bit index into the output buffer */ \
			"shrq	%%cl,%%r8	\n\t"			/* shift the bit scan register up to the bit we just processed */ \
			"incq	%%r11		\n\t"			/* increment the count of set bits */ \
            "incq	%%r10		\n\t"			/* increment the index */ \
			"shrq	$1, %%r8 \n\t"				/* clear the bit */ \
			"jmp 1b		\n\t"					/* loop if so */ \
			"2:		\n\t"						/*  */ \
			"movl	%%r11d, %0 \n\t"			/* return the count of set bits */ \
			: "=r"(result)						\
			: "r"(sieveblock + j), "r"(buffer)	\
			: "xmm0", "xmm1", "xmm2", "xmm3", "r8", "r9", "r10", "r11", "rcx", "cc", "memory");

#endif

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


#elif defined(_WIN64) && defined(_MSC_VER)
	#define SCAN_CLEAN /*nothing*/

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
			local_block = _mm_or_si128(local_block, local_block2); \
			local_block = _mm_or_si128(local_block, local_block3); \
			local_block4 = _mm_load_si128(sieveblock + j + 6); \
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
			local_block = _mm_or_si128(local_block, local_block2); \
			local_block4 = _mm_load_si128(sieveblock + j + 6); \
			local_block = _mm_or_si128(local_block, local_block3); \
			local_block5 = _mm_load_si128(sieveblock + j + 8); \
			local_block = _mm_or_si128(local_block, local_block4); \
			local_block6 = _mm_load_si128(sieveblock + j + 10); \
			local_block = _mm_or_si128(local_block, local_block5); \
			local_block7 = _mm_load_si128(sieveblock + j + 12); \
			local_block = _mm_or_si128(local_block, local_block6); \
			local_block8 = _mm_load_si128(sieveblock + j + 14); \
			local_block = _mm_or_si128(local_block, local_block7); \
			local_block = _mm_or_si128(local_block, local_block8); \
			result = _mm_movemask_epi8(local_block); \
		} while (0);

#if defined(USE_AVX2)


    #define SIEVE_SCAN_32_VEC					\
        __m256i v_blk = _mm256_load_si256(sieveblock + j); \
        uint32_t pos, msk32 = _mm256_movemask_epi8(v_blk); \
        result = 0; \
        while (_BitScanForward(&pos, msk32)) { \
            buffer[result++] = pos; \
            _reset_lsb(msk32); \
        }

    #define SIEVE_SCAN_64_VEC				\
        __m256i v_blk = _mm256_or_si256(_mm256_load_si256(sieveblock + j + 4), _mm256_load_si256(sieveblock + j)); \
        uint32_t pos; \
        uint64_t msk64 = _mm256_movemask_epi8(v_blk); \
        result = 0; \
        if (msk64 > 0) { \
            v_blk = _mm256_load_si256(sieveblock + j + 4); \
            msk64 = ((uint64_t)(_mm256_movemask_epi8(v_blk)) << 32); \
            v_blk = _mm256_load_si256(sieveblock + j); \
            msk64 |= _mm256_movemask_epi8(v_blk); \
            while (_BitScanForward64(&pos, msk64)) { \
                        buffer[result++] = (uint8_t)pos; \
                        _reset_lsb64(msk64); \
            } }
#endif

#else	/* compiler not recognized*/

#define SCAN_CLEAN /*nothing*/
#undef SIMD_SIEVE_SCAN
#undef SIMD_SIEVE_SCAN_VEC

#endif

#define SCAN_MASK 0x8080808080808080ULL

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

	// also need to bail on location 65534, because otherwise the 8x sse2 asm division
	// can fail.  this is because we add 1 to the block loc and then add the correction
	// factor on top of that, which can overflow if blockloc >= 65534.

	// I think that throwing away 3/1000th of 1 percent of the sieve hits in exchange
	// for the speedups associated with 8x sse2 asm division and compression of
	// small primes is worth it (on 64k builds only).


int check_relations_siqs_4_sse2(uint32_t blocknum, uint8_t parity, 
	static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//unrolled x32; for medium inputs

	uint32_t i,j,it=sconf->qs_blocksize>>3;
	uint32_t thisloc;
	uint64_t *sieveblock;

	sieveblock = (uint64_t *)dconf->sieve;
	dconf->num_reports = 0;

#if defined(SIMD_SIEVE_SCAN)

	//check for relations
	for (j=0;j<it;j+=4)
	{
		uint32_t result = 0;

		SIEVE_SCAN_32;

		if (result == 0)
			continue;

		//at least one passed the check, find which one(s) and pass to 
		//trial division stage
		for (i=0; i<4; i++)
		{
			uint32_t k;

			//check 8 locations simultaneously
			if ((sieveblock[j + i] & SCAN_MASK) == (uint64_t)(0))
				continue;

			for (k=0;k<8;k++)
			{
				thisloc = ((j+i)<<3) + k;
				if ((dconf->sieve[thisloc] & 0x80) == 0)
					continue;

				//see discussion near line 323
				if ((thisloc >=	65534) || (thisloc == 0))
					continue;

				// log this report
				if (dconf->num_reports < MAX_SIEVE_REPORTS)
					dconf->reports[dconf->num_reports++] = thisloc;
			}
		}
	}

	// make it safe to perform floating point
	SCAN_CLEAN;

#else

	for (j=0;j<it;j+=4)	
	{
		uint32_t k;

		if (((sieveblock[j] | sieveblock[j+1] | sieveblock[j+2] | sieveblock[j+3]
			) & SCAN_MASK) == (uint64_t)(0))
			continue;

		//at least one passed the check, find which one(s) and pass to 
		//trial division stage
		for (i=0; i<4; i++)
		{
			//check 8 locations simultaneously
			if ((sieveblock[j + i] & SCAN_MASK) == (uint64_t)(0))
				continue;

			for (k=0;k<8;k++)
			{
				thisloc = ((j+i)<<3) + k;
				if ((dconf->sieve[thisloc] & 0x80) == 0)
					continue;

				//see discussion near line 323
				if ((thisloc >=	65534) || (thisloc == 0))
					continue;

				// log this report
				if (dconf->num_reports < MAX_SIEVE_REPORTS)
					dconf->reports[dconf->num_reports++] = thisloc;
			}
		}
	}


#endif


	if (dconf->num_reports >= MAX_SIEVE_REPORTS)
		dconf->num_reports = MAX_SIEVE_REPORTS-1;

    dconf->total_reports += dconf->num_reports;
    dconf->total_blocks++;

	//remove small primes, and test if its worth continuing for each report
	filter_SPV(parity, dconf->sieve,dconf->numB-1,blocknum,sconf,dconf);
	tdiv_med_ptr(parity, dconf->numB-1,blocknum,sconf,dconf);
	resieve_med_ptr(parity, dconf->numB-1,blocknum,sconf,dconf);

	// factor all reports in this block
	for (j=0; j<dconf->num_reports; j++)
	{
		if (dconf->valid_Qs[j])
		{
            dconf->total_surviving_reports++;
            tdiv_LP_ptr(j, parity, blocknum, sconf, dconf);
			trial_divide_Q_siqs(j, parity, dconf->numB-1, blocknum,sconf,dconf);
		}
	}

	return 0;
}

int check_relations_siqs_8_sse2(uint32_t blocknum, uint8_t parity, 
	static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//unrolled x64; for large inputs
	uint32_t i,j,it=sconf->qs_blocksize>>3;
	uint32_t thisloc;
	uint64_t *sieveblock;

	sieveblock = (uint64_t *)dconf->sieve;
	dconf->num_reports = 0;

#if defined(SIMD_SIEVE_SCAN)

	//check for relations
	for (j=0;j<it;j+=8)
	{
		uint32_t result = 0;

		SIEVE_SCAN_64;

		if (result == 0)
			continue;

		//at least one passed the check, find which one(s) and pass to 
		//trial division stage
		for (i=0; i<8; i++)
		{
			uint32_t k;

			//check 8 locations simultaneously
			if ((sieveblock[j + i] & SCAN_MASK) == (uint64_t)(0))
				continue;

			for (k=0;k<8;k++)
			{
				thisloc = ((j+i)<<3) + k;
				if ((dconf->sieve[thisloc] & 0x80) == 0)
					continue;

				//see discussion near line 323
				if ((thisloc >=	65534) || (thisloc == 0))
					continue;

				// log this report
				if (dconf->num_reports < MAX_SIEVE_REPORTS)
					dconf->reports[dconf->num_reports++] = thisloc;
			}
		}
	}

	// make it safe to perform floating point
	SCAN_CLEAN;

#else

	for (j=0;j<it;j+=8)	
	{
		uint32_t k;

		if (((sieveblock[j] | sieveblock[j+1] | sieveblock[j+2] | sieveblock[j+3] |
		      sieveblock[j+4] | sieveblock[j+5] | sieveblock[j+6] | sieveblock[j+7]
			) & SCAN_MASK) == (uint64_t)(0))
			continue;

		//at least one passed the check, find which one(s) and pass to 
		//trial division stage
		for (i=0; i<8; i++)
		{
			//check 8 locations simultaneously
			if ((sieveblock[j + i] & SCAN_MASK) == (uint64_t)(0))
				continue;

			for (k=0;k<8;k++)
			{
				thisloc = ((j+i)<<3) + k;
				if ((dconf->sieve[thisloc] & 0x80) == 0)
					continue;

				//see discussion near line 323
				if ((thisloc >=	65534) || (thisloc == 0))
					continue;

				// log this report
				if (dconf->num_reports < MAX_SIEVE_REPORTS)
					dconf->reports[dconf->num_reports++] = thisloc;
			}
		}
	}


#endif


	if (dconf->num_reports >= MAX_SIEVE_REPORTS)
		dconf->num_reports = MAX_SIEVE_REPORTS-1;

    dconf->total_reports += dconf->num_reports;
    dconf->total_blocks++;

	//remove small primes, and test if its worth continuing for each report
	filter_SPV(parity, dconf->sieve, dconf->numB-1, blocknum,sconf,dconf);
	tdiv_med_ptr(parity, dconf->numB-1,blocknum,sconf,dconf);
	resieve_med_ptr(parity, dconf->numB-1,blocknum,sconf,dconf);

	// factor all reports in this block
	for (j=0; j<dconf->num_reports; j++)
	{
		if (dconf->valid_Qs[j])
		{
            dconf->total_surviving_reports++;
            tdiv_LP_ptr(j, parity, blocknum, sconf, dconf);
			trial_divide_Q_siqs(j, parity, dconf->numB-1, blocknum,sconf,dconf);
		}
	}

	return 0;
}

int check_relations_siqs_16_sse2(uint32_t blocknum, uint8_t parity,
    static_conf_t *sconf, dynamic_conf_t *dconf)
{
    //unrolled x128; for large inputs
    uint32_t i, j, it = sconf->qs_blocksize >> 3;
    uint32_t thisloc;
    uint32_t num_reports = 0;
    uint64_t *sieveblock;

    dconf->num_reports = 0;
    sieveblock = (uint64_t *)dconf->sieve;


#if defined(SIMD_SIEVE_SCAN)

	//check for relations
	for (j=0;j<it;j+=16)
	{
		uint32_t result = 0;

		SIEVE_SCAN_128;

		if (result == 0)
			continue;

		//at least one passed the check, find which one(s) and pass to 
		//trial division stage
		for (i=0; i<16; i++)
		{
			uint32_t k;

			//check 8 locations simultaneously
			if ((sieveblock[j + i] & SCAN_MASK) == (uint64_t)(0))
				continue;

			for (k=0;k<8;k++)
			{
				thisloc = ((j+i)<<3) + k;
				if ((dconf->sieve[thisloc] & 0x80) == 0)
					continue;

				//see discussion near line 323
				if ((thisloc >=	65534) || (thisloc == 0))
					continue;

				// log this report
				if (dconf->num_reports < MAX_SIEVE_REPORTS)
					dconf->reports[dconf->num_reports++] = thisloc;
			}
		}
	}

	// make it safe to perform floating point
	SCAN_CLEAN;

#else

	for (j=0;j<it;j+=16)	
	{
		uint32_t k;

		if (((sieveblock[j] | sieveblock[j+1] | sieveblock[j+2] | sieveblock[j+3] |
		      sieveblock[j+4] | sieveblock[j+5] | sieveblock[j+6] | sieveblock[j+7] |
			  sieveblock[j+8] | sieveblock[j+9] | sieveblock[j+10] | sieveblock[j+11] |
		      sieveblock[j+12] | sieveblock[j+13] | sieveblock[j+14] | sieveblock[j+15]
				) & SCAN_MASK) == (uint64_t)(0))
			continue;

		//at least one passed the check, find which one(s) and pass to 
		//trial division stage
		for (i=0; i<16; i++)
		{
			//check 8 locations simultaneously
			if ((sieveblock[j + i] & SCAN_MASK) == (uint64_t)(0))
				continue;

			for (k=0;k<8;k++)
			{
				thisloc = ((j+i)<<3) + k;
				if ((dconf->sieve[thisloc] & 0x80) == 0)
					continue;

				//see discussion near line 323
				if ((thisloc >=	65534) || (thisloc == 0))
					continue;

				// log this report
				if (dconf->num_reports < MAX_SIEVE_REPORTS)
					dconf->reports[dconf->num_reports++] = thisloc;
			}
		}
	}


#endif	


	if (dconf->num_reports >= MAX_SIEVE_REPORTS)
		dconf->num_reports = MAX_SIEVE_REPORTS-1;

    dconf->total_reports += dconf->num_reports;
    dconf->total_blocks++;

	//remove small primes, and test if its worth continuing for each report
	filter_SPV(parity, dconf->sieve, dconf->numB-1,blocknum,sconf,dconf);
	tdiv_med_ptr(parity, dconf->numB-1,blocknum,sconf,dconf);
	resieve_med_ptr(parity, dconf->numB-1,blocknum,sconf,dconf);

	// factor all reports in this block
	for (j=0; j<dconf->num_reports; j++)
	{
		if (dconf->valid_Qs[j])
		{
            dconf->total_surviving_reports++;
            tdiv_LP_ptr(j, parity, blocknum, sconf, dconf);
			trial_divide_Q_siqs(j, parity, dconf->numB-1, blocknum,sconf,dconf);
		}
	}

	return 0;
}

#ifdef USE_AVX2
int check_relations_siqs_4_avx2(uint32_t blocknum, uint8_t parity,
    static_conf_t* sconf, dynamic_conf_t* dconf)
{
    //unrolled x32; for medium inputs

    uint32_t i, j, it = sconf->qs_blocksize >> 3;
    uint32_t thisloc;
    uint64_t* sieveblock;

    sieveblock = (uint64_t*)dconf->sieve;
    dconf->num_reports = 0;

    for (j = 0; j < it; j += 4)
    {
        uint32_t result;
        uint8_t buffer[32];

        // assumes BSF.
        // MSVC also assumes reset_lsb.
        // both need BMI1
        SIEVE_SCAN_32_VEC;

        if (result == 0)
            continue;

        if ((dconf->num_reports + result) < MAX_SIEVE_REPORTS)
        {
            for (i = 0; i < result; i++)
            {
                thisloc = (j << 3) + (uint32_t)buffer[i];
                dconf->reports[dconf->num_reports + i] = thisloc;
            }
            dconf->num_reports += result;
        }
        else
        {
            for (i = 0; i < result; i++)
            {
                thisloc = (j << 3) + (uint32_t)buffer[i];

                // log this report
                if (dconf->num_reports < MAX_SIEVE_REPORTS)
                    dconf->reports[dconf->num_reports++] = thisloc;
            }
        }
    }

    if (dconf->num_reports >= MAX_SIEVE_REPORTS)
        dconf->num_reports = MAX_SIEVE_REPORTS - 1;

    dconf->total_reports += dconf->num_reports;
    dconf->total_blocks++;

    //remove small primes, and test if its worth continuing for each report
    filter_SPV(parity, dconf->sieve, dconf->numB - 1, blocknum, sconf, dconf);
    tdiv_med_ptr(parity, dconf->numB - 1, blocknum, sconf, dconf);
    resieve_med_ptr(parity, dconf->numB - 1, blocknum, sconf, dconf);

    // factor all reports in this block
    for (j = 0; j < dconf->num_reports; j++)
    {
        if (dconf->valid_Qs[j])
        {
            dconf->total_surviving_reports++;
            tdiv_LP_ptr(j, parity, blocknum, sconf, dconf);
            trial_divide_Q_siqs(j, parity, dconf->numB - 1, blocknum, sconf, dconf);
        }
    }

    return 0;
}

int check_relations_siqs_8_avx2(uint32_t blocknum, uint8_t parity,
    static_conf_t* sconf, dynamic_conf_t* dconf)
{
    //unrolled x64; for large inputs
    uint32_t i, j, it = sconf->qs_blocksize >> 3;
    uint32_t thisloc;
    uint64_t* sieveblock;

    sieveblock = (uint64_t*)dconf->sieve;
    dconf->num_reports = 0;

    for (j = 0; j < it; j += 8)
    {
        uint32_t result;
        uint8_t buffer[64];

        // assumes BSF.
        // MSVC also assumes reset_lsb.
        // both need BMI1
        SIEVE_SCAN_64_VEC;

        if (result == 0)
            continue;

        if ((dconf->num_reports + result) < MAX_SIEVE_REPORTS)
        {
            for (i = 0; i < result; i++)
            {
                thisloc = (j << 3) + (uint32_t)buffer[i];
                dconf->reports[dconf->num_reports + i] = thisloc;
            }
            dconf->num_reports += result;
        }
        else
        {
            for (i = 0; i < result; i++)
            {
                thisloc = (j << 3) + (uint32_t)buffer[i];

                // log this report
                if (dconf->num_reports < MAX_SIEVE_REPORTS)
                    dconf->reports[dconf->num_reports++] = thisloc;
            }
        }
    }

    if (dconf->num_reports >= MAX_SIEVE_REPORTS)
    {
        dconf->num_reports = MAX_SIEVE_REPORTS - 1;
    }

    dconf->total_reports += dconf->num_reports;
    dconf->total_blocks++;

    //remove small primes, and test if its worth continuing for each report
    filter_SPV(parity, dconf->sieve, dconf->numB - 1, blocknum, sconf, dconf);
    tdiv_med_ptr(parity, dconf->numB - 1, blocknum, sconf, dconf);
    resieve_med_ptr(parity, dconf->numB - 1, blocknum, sconf, dconf);

    // factor all reports in this block
    for (j = 0; j < dconf->num_reports; j++)
    {
        if (dconf->valid_Qs[j])
        {
            dconf->total_surviving_reports++;
            tdiv_LP_ptr(j, parity, blocknum, sconf, dconf);
            trial_divide_Q_siqs(j, parity, dconf->numB - 1, blocknum, sconf, dconf);
        }
    }

    return 0;
}

int check_relations_siqs_16_avx2(uint32_t blocknum, uint8_t parity,
    static_conf_t* sconf, dynamic_conf_t* dconf)
{
    //unrolled x128; for large inputs
    uint32_t i, j, it = sconf->qs_blocksize >> 3;
    uint32_t thisloc;
    uint32_t num_reports = 0;
    uint64_t* sieveblock;

    dconf->num_reports = 0;
    sieveblock = (uint64_t*)dconf->sieve;

    for (j = 0; j < it; j += 8)
    {
        uint32_t result;
        uint8_t buffer[64];

        SIEVE_SCAN_64_VEC;

        if (result == 0)
            continue;

        if ((dconf->num_reports + result) < MAX_SIEVE_REPORTS)
        {
            for (i = 0; i < result; i++)
            {
                thisloc = (j << 3) + (uint32_t)buffer[i];
                dconf->reports[dconf->num_reports + i] = thisloc;
            }
            dconf->num_reports += result;
        }
        else
        {
            for (i = 0; i < result; i++)
            {
                thisloc = (j << 3) + (uint32_t)buffer[i];

                // log this report
                if (dconf->num_reports < MAX_SIEVE_REPORTS)
                    dconf->reports[dconf->num_reports++] = thisloc;
            }
        }
    }

    if (dconf->num_reports >= MAX_SIEVE_REPORTS)
        dconf->num_reports = MAX_SIEVE_REPORTS - 1;

    dconf->total_reports += dconf->num_reports;
    dconf->total_blocks++;

    //remove small primes, and test if its worth continuing for each report
    filter_SPV(parity, dconf->sieve, dconf->numB - 1, blocknum, sconf, dconf);
    tdiv_med_ptr(parity, dconf->numB - 1, blocknum, sconf, dconf);
    resieve_med_ptr(parity, dconf->numB - 1, blocknum, sconf, dconf);

    // factor all reports in this block
    for (j = 0; j < dconf->num_reports; j++)
    {
        if (dconf->valid_Qs[j])
        {
            dconf->total_surviving_reports++;
            tdiv_LP_ptr(j, parity, blocknum, sconf, dconf);

#if defined( USE_SS_SEARCH ) && defined(USE_POLY_BUCKET_SS)

			if ((sconf->factor_base->ss_start_B > 0) &&
				(sconf->factor_base->ss_start_B < sconf->factor_base->B))
			{
				tdiv_SS(j, parity, blocknum, sconf, dconf);
			}

#endif

            trial_divide_Q_siqs(j, parity, dconf->numB - 1, blocknum, sconf, dconf);
        }
    }

    return 0;
}

#endif

#ifdef USE_AVX512F
#include <immintrin.h>
int check_relations_siqs_16_avx512(uint32_t blocknum, uint8_t parity,
    static_conf_t* sconf, dynamic_conf_t* dconf)
{
    //unrolled x128; for large inputs
    uint32_t j;
    uint32_t thisloc;
    uint64_t* sieveblock;
    uint32_t* sieveblock32;

    if (sconf->use_dlp == 2)
    {
		// In the triple large prime variation we allow
		// the cutoff value to exceed 127, so that we can 
		// no longer use the hi-bit mask trick to detect 
		// likely smooth locations.  With AVX512BW it doesn't 
		// matter since we can do a direct byte compare on 64
		// simultaneous locations.
#ifdef USE_AVX512BW
        if (sconf->obj->HAS_AVX512BW)
        {
            uint32_t cutoff = sconf->tf_closnuf;
            __m512i vmask = _mm512_set1_epi32((cutoff << 24) + (cutoff << 16) + (cutoff << 8) + cutoff);

            dconf->num_reports = 0;
            sieveblock = (uint64_t*)dconf->sieve;
            sieveblock32 = (uint32_t*)dconf->sieve;

            for (j = 0; j < 4096; j += 8)
            {
                __mmask64 r_msk;
                int idx;

                __m512i vsieve = _mm512_load_epi32((__m512i*)(&sieveblock[j]));
                r_msk = _mm512_cmpgt_epu8_mask(vsieve, vmask);
                //r_msk = _mm512_test_epi32_mask(vsieve, vmask);

                thisloc = j * 8;

                while (r_msk > 0)
                {
                    idx = _trail_zcnt64(r_msk);
                    dconf->reports[dconf->num_reports++] = thisloc + idx;
                    r_msk = _blsr_u64(r_msk);
                }

            }
        }
        else
        {

            // I got nothin'!




        }
#else

#endif
    }
    else
    {
        __m512i vmask = _mm512_set1_epi32(0x80808080);

        dconf->num_reports = 0;
        sieveblock = (uint64_t*)dconf->sieve;
        sieveblock32 = (uint32_t*)dconf->sieve;

        for (j = 0; j < 4096; j += 8)
        {
            __mmask16 r_msk;
            int idx;
            int k;

            __m512i vsieve = _mm512_load_epi32((__m512i*)(&sieveblock[j]));
            r_msk = _mm512_test_epi32_mask(vsieve, vmask);

            thisloc = j * 8;

            while (r_msk > 0)
            {
                uint32_t a_msk;

                idx = _trail_zcnt(r_msk);

                // each lit bit identifies 4 possible bytes meeting our criteria
                // for a possible relation.
                a_msk = sieveblock32[(thisloc >> 2) + idx] & 0x80808080;

                do
                {
                    k = _trail_zcnt(a_msk) >> 3;
	                dconf->reports[dconf->num_reports++] = thisloc + k + idx * 4;
                    a_msk = _blsr_u32(a_msk);
                } while (a_msk > 0);

                r_msk = _blsr_u32(r_msk);
            }

        }
    }

    if (dconf->num_reports >= MAX_SIEVE_REPORTS)
        dconf->num_reports = MAX_SIEVE_REPORTS - 1;

    dconf->total_reports += dconf->num_reports;
    dconf->total_blocks++;

    //remove small primes, and test if its worth continuing for each report
    filter_SPV(parity, dconf->sieve, dconf->numB - 1, blocknum, sconf, dconf);
    tdiv_med_ptr(parity, dconf->numB - 1, blocknum, sconf, dconf);
    resieve_med_ptr(parity, dconf->numB - 1, blocknum, sconf, dconf);

    // factor all reports in this block
    for (j = 0; j < dconf->num_reports; j++)
    {
		if (dconf->valid_Qs[j])
		{
			dconf->total_surviving_reports++;
			tdiv_LP_ptr(j, parity, blocknum, sconf, dconf);

#if defined( USE_SS_SEARCH ) && defined(USE_POLY_BUCKET_SS)

			if ((sconf->factor_base->ss_start_B > 0) &&
				(sconf->factor_base->ss_start_B < sconf->factor_base->B))
			{
				tdiv_SS(j, parity, blocknum, sconf, dconf);
			}
		
#endif

            trial_divide_Q_siqs(j, parity, dconf->numB - 1, blocknum, sconf, dconf);
        }
    }

    return 0;
}
#endif
