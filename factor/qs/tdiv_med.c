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

this file contains code implementing 3)


*/


#if defined(GCC_ASM32X) || defined(GCC_ASM64X) || defined(__MINGW32__)

	#if defined(D_HAS_SSE2)

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
		#define MOD_CMP_8X(xtra_bits)															\
			ASM_G (																				\
				"movdqa %%xmm0, %%xmm4 \n\t"	/* copy BLOCKSIZE */							\
				"movdqa %%xmm1, %%xmm7 \n\t"	/* copy block_loc */							\
				"pcmpeqw	%%xmm6, %%xmm6 \n\t"	/* create an array of 1's */				\
				"movdqa (%1), %%xmm3 \n\t"		/* move in primes */							\
				"psubw	%%xmm1, %%xmm4 \n\t"	/* BLOCKSIZE - block_loc */						\
				"psrlw	$15, %%xmm6 \n\t"		/* create an array of 1's */					\
				"paddw	%%xmm6, %%xmm4 \n\t"	/* add in 1's */								\
				"psubw	%%xmm6, %%xmm7 \n\t"	/* substract 1's */								\
				"paddw	(%2), %%xmm4 \n\t"		/* apply corrections */							\
				"movdqa (%4), %%xmm6 \n\t"		/* move in root1s */							\
				"pmulhuw	(%3), %%xmm4 \n\t"	/* (unsigned) multiply by inverses */		\
				"movdqa (%5), %%xmm2 \n\t"		/* move in root2s */							\
				"psrlw	$" xtra_bits ", %%xmm4 \n\t"		/* to get to total shift of 24/26/28 bits */			\
				"paddw	%%xmm3, %%xmm7 \n\t"	/* add primes and block_loc */					\
				"pmullw	%%xmm3, %%xmm4 \n\t"	/* (signed) multiply by primes */				\
				"psubw	%%xmm0, %%xmm7 \n\t"	/* substract blocksize */						\
				"paddw	%%xmm7, %%xmm4 \n\t"	/* add in block_loc + primes - blocksize */		\
				"pcmpeqw	%%xmm4, %%xmm6 \n\t"	/* compare to root1s */						\
				"pcmpeqw	%%xmm4, %%xmm2 \n\t"	/* compare to root2s */						\
				"por	%%xmm6, %%xmm2 \n\t"	/* combine compares */							\
				"pmovmskb %%xmm2, %0 \n\t"		/* export to result */							\
				: "=r" (tmp3)																		\
				: "r" (fbc->prime + i), "r" (fullfb_ptr->correction + i), \
					"r" (fullfb_ptr->small_inv + i), "r" (fbc->root1 + i), \
						"r" (fbc->root2 + i)	\
				: "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7");

#else

		// will miss cases like this:
		//	index = 56, prime = 587, inv = 28581, corr = 1, root1 = 0, root2 = 322, 
		// loc = 1070, result = 0
		// which we can catch by comparing to prime as well as root1/2; but it probably isn't
		// worth it.

				/*"movdqa (%6), %%xmm0 \n\t"		 move in BLOCKSIZE */ \
				/*"movdqa (%7), %%xmm1 \n\t"		 move in block_loc */ \

		#define MOD_CMP_8X(xtra_bits)																		\
			ASM_G (																				\
				"movdqa %%xmm1, %%xmm7 \n\t"	/* copy block_loc */ \
				"movdqa	%%xmm0, %%xmm4 \n\t"	/* copy BLOCKSIZE */ \
				"movdqa (%1), %%xmm3 \n\t"		/* move in primes */							\
				"psubw	%%xmm1, %%xmm4 \n\t"	/* BLOCKSIZE - block_loc */						\
				"paddw	(%2), %%xmm4 \n\t"		/* apply corrections */							\
				"movdqa (%4), %%xmm6 \n\t"		/* move in root1s */							\
				"pmulhuw	(%3), %%xmm4 \n\t"	/* (unsigned) multiply by inverses */		\
				"movdqa (%5), %%xmm2 \n\t"		/* move in root2s */							\
				"psrlw	$" xtra_bits ", %%xmm4 \n\t"		/* to get to total shift of 24/26/28 bits */			\
				"paddw	%%xmm3, %%xmm7 \n\t"	/* add primes and block_loc */					\
				"pmullw	%%xmm3, %%xmm4 \n\t"	/* (signed) multiply by primes */				\
				"psubw	%%xmm0, %%xmm7 \n\t"	/* substract blocksize */						\
				"paddw	%%xmm7, %%xmm4 \n\t"	/* add in block_loc + primes - blocksize */		\
				"pcmpeqw	%%xmm4, %%xmm6 \n\t"	/* compare to root1s */						\
				"pcmpeqw	%%xmm4, %%xmm2 \n\t"	/* compare to root2s */						\
				"por	%%xmm6, %%xmm2 \n\t"	/* combine compares */							\
				"pmovmskb %%xmm2, %0 \n\t"		/* export to result */							\
				: "=r" (tmp3)																	\
				: "r" (fbc->prime + i), "r" (fullfb_ptr->correction + i), \
					"r" (fullfb_ptr->small_inv + i), "r" (fbc->root1 + i), \
						"r" (fbc->root2 + i) \
				: "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7");

#endif

	#endif

#elif defined(MSC_ASM32A)

	#if defined(D_HAS_SSE2)
		//top level sieve scanning with SSE2

		#define MOD_INIT_8X	

#ifdef YAFU_64K

		#define MOD_CMP_8X(xtra_bits)						\
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
				c = _mm_srli_epi16(c, xtra_bits); \
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

		#define MOD_CMP_8X(xtra_bits)						\
			do { \
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
				lp = _mm_load_si128((__m128i *)(fbc->prime + i)); \
				t1 = _mm_sub_epi16(t1, blkloc); \
				c = _mm_load_si128((__m128i *)(fullfb_ptr->correction + i));	\
				t2 = blkloc; \
				sinv = _mm_load_si128((__m128i *)(fullfb_ptr->small_inv + i));	\
				c = _mm_add_epi16(c, t1); \
				lr1 = _mm_load_si128((__m128i *)(fbc->root1 + i)); \
				c = _mm_mulhi_epu16(c, sinv); \
				lr2 = _mm_load_si128((__m128i *)(fbc->root2 + i)); \
				c = _mm_srli_epi16(c, xtra_bits); \
				t2 = _mm_add_epi16(t2, lp); \
				c = _mm_mullo_epi16(c, lp); \
				c = _mm_add_epi16(c, t2); \
				c = _mm_sub_epi16(c, blksz); \
				lr1 = _mm_cmpeq_epi16(lr1, c); \
				lr2 = _mm_cmpeq_epi16(lr2, c); \
				lr2 = _mm_or_si128(lr2, lr1); \
				tmp3 = _mm_movemask_epi8(lr2); \
			} while (0);

#endif

	#endif

#elif defined(_WIN64)

#if defined(D_HAS_SSE2)
		//top level sieve scanning with SSE2

		#define MOD_INIT_8X	

#ifdef YAFU_64K

		#define MOD_CMP_8X(xtra_bits)						\
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
				one = _mm_cmpeq_epi16(blksz, blksz); \
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
				c = _mm_srli_epi16(c, xtra_bits); \
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

		#define MOD_CMP_8X(xtra_bits)						\
			do { \
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
				lp = _mm_load_si128((__m128i *)(fbc->prime + i)); \
				t1 = _mm_sub_epi16(t1, blkloc); \
				c = _mm_load_si128((__m128i *)(fullfb_ptr->correction + i));	\
				t2 = blkloc; \
				sinv = _mm_load_si128((__m128i *)(fullfb_ptr->small_inv + i));	\
				c = _mm_add_epi16(c, t1); \
				lr1 = _mm_load_si128((__m128i *)(fbc->root1 + i)); \
				c = _mm_mulhi_epu16(c, sinv); \
				lr2 = _mm_load_si128((__m128i *)(fbc->root2 + i)); \
				c = _mm_srli_epi16(c, xtra_bits); \
				t2 = _mm_add_epi16(t2, lp); \
				c = _mm_mullo_epi16(c, lp); \
				c = _mm_add_epi16(c, t2); \
				c = _mm_sub_epi16(c, blksz); \
				lr1 = _mm_cmpeq_epi16(lr1, c); \
				lr2 = _mm_cmpeq_epi16(lr2, c); \
				lr2 = _mm_or_si128(lr2, lr1); \
				tmp3 = _mm_movemask_epi8(lr2); \
			} while (0);

#endif

#endif

#else	/* compiler not recognized*/


#endif


#define DIVIDE_ONE_PRIME \
	do \
	{						\
		fb_offsets[++smooth_num] = i;	\
		mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], prime); \
	} while (mpz_tdiv_ui(dconf->Qvals[report_num], prime) == 0); 


#define DIVIDE_RESIEVED_PRIME(j) \
	while (mpz_tdiv_ui(dconf->Qvals[report_num], fbc->prime[i+j]) == 0) \
	{						\
		fb_offsets[++smooth_num] = i+j;	\
		mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], fbc->prime[i+j]);		\
	}


void tdiv_medprimes(uint8_t parity, uint32_t poly_id, uint32_t bnum, 
						 static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//we have flagged this sieve offset as likely to produce a relation
	//nothing left to do now but check and see.
	int i;
	uint32_t bound, tmp, prime, root1, root2, report_num;
	int smooth_num;
	uint32_t *fb_offsets;
	sieve_fb_compressed *fbc;
	fb_element_siqs *fullfb_ptr, *fullfb = sconf->factor_base->list;
	uint32_t block_loc;

#ifdef USE_8X_MOD_ASM
	uint16_t *bl_sizes;
	uint16_t *bl_locs;

	bl_sizes = (uint16_t *)xmalloc_align(8 * sizeof(uint16_t));
	bl_locs = (uint16_t *)xmalloc_align(8 * sizeof(uint16_t));

#endif

	fullfb_ptr = fullfb;
	if (parity)
	{
		fbc = dconf->comp_sieve_n;
	}
	else
	{
		fbc = dconf->comp_sieve_p;
	}

#ifdef QS_TIMING
	gettimeofday(&qs_timing_start, NULL);
#endif

	for (report_num = 0; report_num < dconf->num_reports; report_num++)
	{

		if (!dconf->valid_Qs[report_num])
			continue;

		// for each report, we trial divide, and then either trial divide or
		// resieve.  for the first trial division step, we have either
		// unrolled C routines or SIMD assembly routines to choose from.  The
		// second trial division step only has a C routine - the more optimized
		// path is resieving.
		//
		// the basic idea of trial division is to test if the block location
		// in question lies on the arithmetic progression of a prime.  By the time
		// we get to this routine, the arithmetic progression has been reset for
		// the next sieve block, so we have to do a few manipulations to revert
		// to the "real" progression (add the blocksize back to the roots).  
		// there are various methods for doing the test depending on the size of
		// the primes (and therefore how many times it could have possibly landed
		// in the sieve block).  The most straightforward, and the method the SIMD
		// assembly uses, is to see if the roots (adjusted for the "real" progression)
		// minus the block location in question, divided by the prime, is zero.  if
		// so, this block location is on the arithmetic progression of that prime,
		// and we can proceed to trial divide the prime into Q(x) for this sieve hit.
		// 
		// the basic idea of resieving is to start from the end of the "real"
		// arithmetic progression, repeatedly subtract each prime, and test after
		// each subtraction if we've hit the sieve location in question.  if so,
		// we know this location is on the prime's arithmetic progression and can 
		// proceed to trial divide.  For "large enough" primes, this is very efficient
		// because we only need to do a few subtractions and tests instead of a 
		// division.  Since we are not really doing division tests, and instead are
		// doing multiplication by inverses, and futhermore since we might be doing those
		// multiplications 8 at a time using SIMD, resieving is only a win for the 
		// very largest primes less than 16 bits in size.  
		//
		// for OS/architecture/compilers where resieving isn't implemented, there are
		// further trial division steps instead.  These are more efficient than
		// the "check for exact division of a difference" method described above,
		// but are only implemented in portable C.  See code below for more detail.

		// pull the details of this report to get started.
		fb_offsets = &dconf->fb_offsets[report_num][0];
		smooth_num = dconf->smooth_num[report_num];
		block_loc = dconf->reports[report_num];

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

		// the bound before resieving takes over
#if defined(YAFU_64K)
		bound = sconf->factor_base->fb_14bit_B;
#else
		bound = sconf->factor_base->fb_13bit_B;
#endif

		// although if we are using 8x ASM division,
		// use a lower bound for a while first
#ifdef USE_8X_MOD_ASM
		bound = sconf->factor_base->fb_10bit_B;
#endif
		
		// single-up test until i is a multiple of 8
		while ((uint32_t)i < bound && ((i & 7) != 0))
		{
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];

			//tmp = distance from this sieve block offset to the end of the block
			tmp = BLOCKSIZE - block_loc;
	
			//tmp = tmp/prime + 1 = number of steps to get past the end of the sieve
			//block, which is the state of the sieve now.
			tmp = 1+(uint32_t)(((uint64_t)(tmp + fullfb_ptr->correction[i])
					* (uint64_t)fullfb_ptr->small_inv[i]) >> FOGSHIFT); 
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
		
		MOD_INIT_8X;

		while ((uint32_t)i < bound)
		{
			uint32_t tmp3 = 0;

#ifdef _MSC_VER
			MOD_CMP_8X(8);
#else
			MOD_CMP_8X("8");
#endif

			//if ((((((uint32_t)fbc->root1[i] + BLOCKSIZE - block_loc) % fbc->prime[i]) == 0) ||
			//	((((uint32_t)fbc->root2[i] + BLOCKSIZE - block_loc) % fbc->prime[i]) == 0)) &&
			//	((tmp3 & 0x2) == 0))
			//	printf("index = %d, prime = %u, inv = %u, corr = %u, root1 = %u, root2 = %u, "
			//	"loc = %u, result = %u\n",
			//		i, fbc->prime[i], fullfb_ptr->small_inv[i], 
			//		fullfb_ptr->correction[i], fbc->root1[i], fbc->root2[i], block_loc, tmp3);

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

		bound = sconf->factor_base->fb_12bit_B;
		while ((uint32_t)i < bound)
		{
			uint32_t tmp3 = 0;

#ifdef _MSC_VER
			MOD_CMP_8X(10);
#else
			MOD_CMP_8X("10");
#endif

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

#if defined(YAFU_64K)
		bound = sconf->factor_base->fb_14bit_B;
#else
		bound = sconf->factor_base->fb_13bit_B;
#endif				
		while ((uint32_t)i < bound)
		{
			uint32_t tmp3 = 0;

#ifdef _MSC_VER
			MOD_CMP_8X(12);
#else
			MOD_CMP_8X("12");
#endif

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
		while ((uint32_t)i < bound)
		{
			uint64_t q64;
			uint32_t tmp1 = BLOCKSIZE - block_loc;
			uint32_t tmp2 = BLOCKSIZE - block_loc;
			uint32_t tmp3 = BLOCKSIZE - block_loc;
			uint32_t tmp4 = BLOCKSIZE - block_loc;

			tmp1 = tmp1 + fullfb_ptr->correction[i];
			q64 = (uint64_t)tmp1 * (uint64_t)fullfb_ptr->small_inv[i];
			tmp1 = q64 >> FOGSHIFT; 
			tmp1 = tmp1 + 1;
			tmp1 = block_loc + tmp1 * fullfb_ptr->prime[i];
			
			tmp2 = tmp2 + fullfb_ptr->correction[i+1];
			q64 = (uint64_t)tmp2 * (uint64_t)fullfb_ptr->small_inv[i+1];
			tmp2 = q64 >> FOGSHIFT; 
			tmp2 = tmp2 + 1;
			tmp2 = block_loc + tmp2 * fullfb_ptr->prime[i+1];

			tmp3 = tmp3 + fullfb_ptr->correction[i+2];
			q64 = (uint64_t)tmp3 * (uint64_t)fullfb_ptr->small_inv[i+2];
			tmp3 = q64 >> FOGSHIFT; 
			tmp3 = tmp3 + 1;
			tmp3 = block_loc + tmp3 * fullfb_ptr->prime[i+2];

			tmp4 = tmp4 + fullfb_ptr->correction[i+3];
			q64 = (uint64_t)tmp4 * (uint64_t)fullfb_ptr->small_inv[i+3];
			tmp4 = q64 >>  FOGSHIFT; 
			tmp4 = tmp4 + 1;
			tmp4 = block_loc + tmp4 * fullfb_ptr->prime[i+3];

			tmp1 = tmp1 - BLOCKSIZE;
			tmp2 = tmp2 - BLOCKSIZE;
			tmp3 = tmp3 - BLOCKSIZE;
			tmp4 = tmp4 - BLOCKSIZE;

			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];

			if (tmp1 == root1 || tmp1 == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;

			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];

			if (tmp2 == root1 || tmp2 == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;

			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];

			if (tmp3 == root1 || tmp3 == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;

			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];

			if (tmp4 == root1 || tmp4 == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;

		}

		//now cleanup any that don't fit in the last batch
		while ((uint32_t)i < bound)
		{
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];

			tmp = BLOCKSIZE - block_loc;
			tmp = 1+(uint32_t)(((uint64_t)(tmp + fullfb_ptr->correction[i])
					* (uint64_t)fullfb_ptr->small_inv[i]) >> FOGSHIFT); 
			tmp = block_loc + tmp*prime;
			tmp = tmp - BLOCKSIZE;

			if (tmp == root1 || tmp == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;
		}

		// single-up test until i is a multiple of 8
		while ((uint32_t)i < bound && ((i & 7) != 0))
		{
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];

			//tmp = distance from this sieve block offset to the end of the block
			tmp = BLOCKSIZE - block_loc;
	
			//tmp = tmp/prime + 1 = number of steps to get past the end of the sieve
			//block, which is the state of the sieve now.
			tmp = 1+(uint32_t)(((uint64_t)(tmp + fullfb_ptr->correction[i])
					* (uint64_t)fullfb_ptr->small_inv[i]) >> FOGSHIFT_2); 
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
		while ((uint32_t)i < bound)
		{
			uint64_t q64;
			uint32_t tmp1 = BLOCKSIZE - block_loc;
			uint32_t tmp2 = BLOCKSIZE - block_loc;
			uint32_t tmp3 = BLOCKSIZE - block_loc;
			uint32_t tmp4 = BLOCKSIZE - block_loc;

			tmp1 = tmp1 + fullfb_ptr->correction[i];
			q64 = (uint64_t)tmp1 * (uint64_t)fullfb_ptr->small_inv[i];
			tmp1 = q64 >> FOGSHIFT; 
			tmp1 = tmp1 + 1;
			tmp1 = block_loc + tmp1 * fullfb_ptr->prime[i];
			
			tmp2 = tmp2 + fullfb_ptr->correction[i+1];
			q64 = (uint64_t)tmp2 * (uint64_t)fullfb_ptr->small_inv[i+1];
			tmp2 = q64 >> FOGSHIFT; 
			tmp2 = tmp2 + 1;
			tmp2 = block_loc + tmp2 * fullfb_ptr->prime[i+1];

			tmp3 = tmp3 + fullfb_ptr->correction[i+2];
			q64 = (uint64_t)tmp3 * (uint64_t)fullfb_ptr->small_inv[i+2];
			tmp3 = q64 >> FOGSHIFT; 
			tmp3 = tmp3 + 1;
			tmp3 = block_loc + tmp3 * fullfb_ptr->prime[i+2];

			tmp4 = tmp4 + fullfb_ptr->correction[i+3];
			q64 = (uint64_t)tmp4 * (uint64_t)fullfb_ptr->small_inv[i+3];
			tmp4 = q64 >>  FOGSHIFT; 
			tmp4 = tmp4 + 1;
			tmp4 = block_loc + tmp4 * fullfb_ptr->prime[i+3];

			tmp1 = tmp1 - BLOCKSIZE;
			tmp2 = tmp2 - BLOCKSIZE;
			tmp3 = tmp3 - BLOCKSIZE;
			tmp4 = tmp4 - BLOCKSIZE;

			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];

			if (tmp1 == root1 || tmp1 == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;

			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];

			if (tmp2 == root1 || tmp2 == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;

			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];

			if (tmp3 == root1 || tmp3 == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;

			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];

			if (tmp4 == root1 || tmp4 == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;

		}

		//now cleanup any that don't fit in the last batch
		while ((uint32_t)i < bound)
		{
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];

			tmp = BLOCKSIZE - block_loc;
			tmp = 1+(uint32_t)(((uint64_t)(tmp + fullfb_ptr->correction[i])
					* (uint64_t)fullfb_ptr->small_inv[i]) >> FOGSHIFT_2); 
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

		// either after 8x SSE2 ASM, or standard trial division, record
		// how many factors we've found so far
		dconf->smooth_num[report_num] = smooth_num;	

	}

#ifdef QS_TIMING
	gettimeofday (&qs_timing_stop, NULL);
    TF_STG2 += ytools_difftime (&qs_timing_start, &qs_timing_stop);
#endif

#ifdef USE_8X_MOD_ASM
	align_free(bl_sizes);
	align_free(bl_locs);
#endif

	return;
}


