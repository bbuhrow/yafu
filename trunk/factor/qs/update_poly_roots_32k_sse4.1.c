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
#include "util.h"
#include "common.h"
#include "poly_macros_32k.h"
#include "poly_macros_common.h"
#include "poly_macros_common_sse4.1.h"

// protect sse41 code under MSVC builds.  USE_SSE41 should be manually
// enabled at the top of qs.h for MSVC builds on supported hardware
#ifdef USE_SSE41

//this is in the poly library, even though the bulk of the time is spent
//bucketizing large primes, because it's where the roots of a poly are updated
void nextRoots_32k_sse41(static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//update the roots 
	sieve_fb_compressed *fb_p = dconf->comp_sieve_p;
	sieve_fb_compressed *fb_n = dconf->comp_sieve_n;
	int *rootupdates = dconf->rootupdates;

	update_t update_data = dconf->update_data;

	uint32 startprime = 2;
	uint32 bound = sconf->factor_base->B;

	char v = dconf->curr_poly->nu[dconf->numB];
	char sign = dconf->curr_poly->gray[dconf->numB];
	int *ptr;
	uint16 *sm_ptr;

	lp_bucket *lp_bucket_p = dconf->buckets;
	uint32 med_B = sconf->factor_base->med_B;
	uint32 large_B = sconf->factor_base->large_B;

	uint32 j, interval; //, fb_offset;
	int k,numblocks;
	uint32 root1, root2, prime;

	int bound_index=0;
	int check_bound = BUCKET_ALLOC/2 - 1;
	uint32 bound_val = med_B;
	uint32 *numptr_p, *numptr_n, *sliceptr_p,*sliceptr_n;
	
#if !defined(USE_POLY_SSE2_ASM) || defined(PROFILING)
	uint32 *bptr;
	int bnum, room;
#endif

	uint8 logp=0;
	polysieve_t helperstruct;

	numblocks = sconf->num_blocks;
	interval = numblocks << 15;
	
	if (lp_bucket_p->alloc_slices != 0) //NULL)
	{
		lp_bucket_p->fb_bounds[0] = med_B;

		sliceptr_p = lp_bucket_p->list;
		sliceptr_n = lp_bucket_p->list + (numblocks << BUCKET_BITS);

		numptr_p = lp_bucket_p->num;
		numptr_n = lp_bucket_p->num + numblocks;
		
		//reuse this for a sec...
		prime = 2*numblocks*lp_bucket_p->alloc_slices;

		//reset lp_buckets
		for (j=0;j<prime;j++)
			numptr_p[j] = 0;
	
		lp_bucket_p->num_slices = 0;

	}
	else
	{
		sliceptr_p = NULL;
		sliceptr_n = NULL;
		numptr_p = NULL;
		numptr_n = NULL;
	}

	k=0;
	ptr = &rootupdates[(v-1) * bound + startprime];	

	if (sign > 0)
	{
#ifdef QS_TIMING
		gettimeofday(&qs_timing_start, NULL);
#endif

		for (j=startprime;j<sconf->sieve_small_fb_start;j++,ptr++)
		{
			prime = update_data.prime[j];
			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];

			COMPUTE_NEXT_ROOTS_P;

			//we don't sieve these, so ordering doesn't matter
			update_data.firstroots1[j] = root1;
			update_data.firstroots2[j] = root2;

			fb_p->root1[j] = (uint16)root1;
			fb_p->root2[j] = (uint16)root2;
			fb_n->root1[j] = (uint16)(prime - root2);
			fb_n->root2[j] = (uint16)(prime - root1);
			if (fb_n->root1[j] == prime)
				fb_n->root1[j] = 0;
			if (fb_n->root2[j] == prime)
				fb_n->root2[j] = 0;

		}

		// do one at a time up to the 10bit boundary, where
		// we can start doing things 8 at a time and be
		// sure we can use aligned moves (static_data_init).		
		for (j=sconf->sieve_small_fb_start; 
			j < sconf->factor_base->fb_10bit_B; j++, ptr++)
		{
			// as soon as we are aligned, use more efficient sse2 based methods...
			if ((j&7) == 0)
				break;

			prime = update_data.prime[j];
			root1 = (uint32)update_data.sm_firstroots1[j];
			root2 = (uint32)update_data.sm_firstroots2[j];

			COMPUTE_NEXT_ROOTS_P;

			if (root2 < root1)
			{
				update_data.sm_firstroots1[j] = (uint16)root2;
				update_data.sm_firstroots2[j] = (uint16)root1;

				fb_p->root1[j] = (uint16)root2;
				fb_p->root2[j] = (uint16)root1;
				fb_n->root1[j] = (uint16)(prime - root1);
				fb_n->root2[j] = (uint16)(prime - root2);
			}
			else
			{
				update_data.sm_firstroots1[j] = (uint16)root1;
				update_data.sm_firstroots2[j] = (uint16)root2;

				fb_p->root1[j] = (uint16)root1;
				fb_p->root2[j] = (uint16)root2;
				fb_n->root1[j] = (uint16)(prime - root2);
				fb_n->root2[j] = (uint16)(prime - root1);
			}
		}
		
#if defined(GCC_ASM64X) || defined(_MSC_VER) //NOTDEF //GCC_ASM64X
		
		// update 8 at a time using SSE2 and no branching
		sm_ptr = &dconf->sm_rootupdates[(v-1) * bound];
		{
			small_update_t h;
			
			h.first_r1 = update_data.sm_firstroots1;		// 0
			h.first_r2 = update_data.sm_firstroots2;		// 8
			h.fbp1 = fb_p->root1;							// 16
			h.fbp2 = fb_p->root2;							// 24
			h.fbn1 = fb_n->root1;							// 32
			h.fbn2 = fb_n->root2;							// 40
			h.primes = fb_p->prime;							// 48
			h.updates = sm_ptr;								// 56
			h.start = j;									// 64
			h.stop = sconf->factor_base->fb_15bit_B;		// 68
			if ((h.stop - 8) > h.start)
				h.stop -= 8;

			COMPUTE_8X_SMALL_PROOTS;
			
			j = h.stop;
		}	

#else
		ptr = &dconf->rootupdates[(v-1) * bound + sconf->factor_base->fb_10bit_B];
		for (j=sconf->factor_base->fb_10bit_B; j < sconf->factor_base->fb_15bit_B; j++, ptr++)
		{
			prime = update_data.prime[j];
			root1 = update_data.sm_firstroots1[j];
			root2 = update_data.sm_firstroots2[j];

			COMPUTE_NEXT_ROOTS_P;

			if (root2 < root1)
			{
				update_data.sm_firstroots1[j] = (uint16)root2;
				update_data.sm_firstroots2[j] = (uint16)root1;

				fb_p->root1[j] = (uint16)root2;
				fb_p->root2[j] = (uint16)root1;
				fb_n->root1[j] = (uint16)(prime - root1);
				fb_n->root2[j] = (uint16)(prime - root2);
			}
			else
			{
				update_data.sm_firstroots1[j] = (uint16)root1;
				update_data.sm_firstroots2[j] = (uint16)root2;

				fb_p->root1[j] = (uint16)root1;
				fb_p->root2[j] = (uint16)root2;
				fb_n->root1[j] = (uint16)(prime - root2);
				fb_n->root2[j] = (uint16)(prime - root1);
			}
		}
#endif		

		// assembly code may not get all the way to 15 bits since we 
		// do things in blocks of 8 there.  Make sure we are at the 15 bit
		// boundary before we switch to using update_data.firstroots1/2.
		// this should only run a few iterations, if any.
		ptr = &dconf->rootupdates[(v-1) * bound + j];
		for ( ; j < med_B; j++, ptr++)
		{
			prime = update_data.prime[j];
			root1 = (uint16)update_data.sm_firstroots1[j];
			root2 = (uint16)update_data.sm_firstroots2[j];

			if ((prime > 32768) && ((j&7) == 0))
				break;

			COMPUTE_NEXT_ROOTS_P;

			if (root2 < root1)
			{
				update_data.sm_firstroots1[j] = (uint16)root2;
				update_data.sm_firstroots2[j] = (uint16)root1;

				fb_p->root1[j] = (uint16)root2;
				fb_p->root2[j] = (uint16)root1;
				fb_n->root1[j] = (uint16)(prime - root1);
				fb_n->root2[j] = (uint16)(prime - root2);
			}
			else
			{
				update_data.sm_firstroots1[j] = (uint16)root1;
				update_data.sm_firstroots2[j] = (uint16)root2;

				fb_p->root1[j] = (uint16)root1;
				fb_p->root2[j] = (uint16)root2;
				fb_n->root1[j] = (uint16)(prime - root2);
				fb_n->root2[j] = (uint16)(prime - root1);
			}
		}	

#if defined(GCC_ASM64X) || (defined(_MSC_VER) && defined (_WIN64))
		// update 8 at a time using SSE2 and no branching
		sm_ptr = &dconf->sm_rootupdates[(v-1) * bound];
		{
			small_update_t h;
			
			h.first_r1 = update_data.sm_firstroots1;		// 0
			h.first_r2 = update_data.sm_firstroots2;		// 8
			h.fbp1 = fb_p->root1;							// 16
			h.fbp2 = fb_p->root2;							// 24
			h.fbn1 = fb_n->root1;							// 32
			h.fbn2 = fb_n->root2;							// 40
			h.primes = fb_p->prime;							// 48
			h.updates = sm_ptr;								// 56
			h.start = j;									// 64
			h.stop = med_B;									// 68

			COMPUTE_8X_SMALL_PROOTS_SSE41;
			
			j = h.stop;
		}	
		sm_ptr = &dconf->sm_rootupdates[(v-1) * bound + j];

#else

		for ( ; j < med_B; j++, ptr++)
		{
			prime = update_data.prime[j];
			root1 = update_data.sm_firstroots1[j];
			root2 = update_data.sm_firstroots2[j];

			COMPUTE_NEXT_ROOTS_P;

			if (root2 < root1)
			{
				update_data.sm_firstroots1[j] = (uint16)root2;
				update_data.sm_firstroots2[j] = (uint16)root1;

				fb_p->root1[j] = (uint16)root2;
				fb_p->root2[j] = (uint16)root1;
				fb_n->root1[j] = (uint16)(prime - root1);
				fb_n->root2[j] = (uint16)(prime - root2);
			}
			else
			{
				update_data.sm_firstroots1[j] = (uint16)root1;
				update_data.sm_firstroots2[j] = (uint16)root2;

				fb_p->root1[j] = (uint16)root1;
				fb_p->root2[j] = (uint16)root2;
				fb_n->root1[j] = (uint16)(prime - root2);
				fb_n->root2[j] = (uint16)(prime - root1);
			}
		}


#endif


#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		POLY_STG2 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);

		gettimeofday(&qs_timing_start, NULL);
#endif

		bound_index = 0;
		bound_val = med_B;
		check_bound = med_B + BUCKET_ALLOC/2;
		

		
#if defined(USE_POLY_SSE2_ASM) && defined(GCC_ASM64X) && !defined(PROFILING)
		logp = update_data.logp[med_B-1];

		if (med_B % 16 != 0)
		{
			printf("med_B must be divisible by 16!\n");
			exit(-1);
		}
		if ((large_B - med_B) % 16 != 0)
		{
			printf("med range must be divisible by 16!\n");
			exit(-1);
		}

		//load up our helper struct, so we don't have a list
		//of a 100 things to stick into our asm block
		helperstruct.numptr_n = numptr_n;			//0
		helperstruct.numptr_p = numptr_p;			//8
		helperstruct.sliceptr_n = sliceptr_n;		//16
		helperstruct.sliceptr_p = sliceptr_p;		//24
		helperstruct.update_data_prime = update_data.prime;	//32
		helperstruct.update_data_root1 = update_data.firstroots1;	//40
		helperstruct.update_data_root2 = update_data.firstroots2;	//48
		helperstruct.update_data_logp = update_data.logp;	//56
		helperstruct.lp_bucket_p = lp_bucket_p;		//64
		helperstruct.ptr = &rootupdates[(v-1) * bound];						//72
		helperstruct.large_B = med_B;				//80
		helperstruct.B = large_B;					//84
		helperstruct.interval = interval;			//88
		helperstruct.numblocks = numblocks;			//92
		helperstruct.bound_val = bound_val;			//96
		helperstruct.bound_index = bound_index;		//100
		helperstruct.check_bound = check_bound;		//104
		helperstruct.logp = logp;					//108
		helperstruct.intervalm1 = interval-1;		//112

		ASM_G (		\
			"movq	%0,%%rsi \n\t"					/* move helperstruct into rsi */ \
			"movl   80(%%rsi,1),%%r15d \n\t"		/* large_B = j = r15d */ \
													/* do the loop comparison */ \
			"cmpl   84(%%rsi,1),%%r15d \n\t"		/* j >= bound ? */ \
			"jae    9f	\n\t"						/* jump to end of loop, if test fails */ \
			"8: \n\t"	\
				/* ================================================ */	\
				/* ========== BEGIN CHECK_NEW_SLICE BLOCK ========= */	\
				/* ================================================ */	\
			CHECK_NEW_SLICE_ASM	\
				/* ================================================ */	\
				/* ============ BEGIN - GET NEW ROOTS 1 =========== */	\
				/* ================================================ */	\
			"movq	72(%%rsi,1),%%rdi \n\t"			/* edi = ptr */ \
			"movq	40(%%rsi,1),%%r14 \n\t"			/* move updata_data.root1 pointer into r14 */ \
			"movdqa (%%rdi,%%r15,4), %%xmm3 \n\t"	/* xmm3 = next 4 values of rootupdates */ \
			"movq	48(%%rsi,1),%%r13 \n\t"			/* move updata_data.root2 pointer into r13 */ \
			"movdqa (%%r14,%%r15,4), %%xmm1 \n\t"	/* xmm1 = next 4 values of root1 */ \
			"psubd	%%xmm3, %%xmm1 \n\t"			/* root1 -= ptr */ \
			"movdqa (%%r13,%%r15,4), %%xmm2 \n\t"	/* xmm2 = next 4 values of root2 */ \
			"movq	32(%%rsi,1),%%r12 \n\t"			/* move updata_data.prime pointer into r12 */ \
			"psubd	%%xmm3, %%xmm2 \n\t"			/* root2 -= ptr */ \
			"pxor	%%xmm4, %%xmm4 \n\t"			/* zero xmm4 */ \
			"pxor	%%xmm5, %%xmm5 \n\t"			/* zero xmm5 */ \
			"movdqa (%%r12,%%r15,4), %%xmm0 \n\t"	/* xmm0 = next 4 primes */ \
			"pcmpgtd	%%xmm1, %%xmm4 \n\t"		/* signed comparison: 0 > root1? if so, set xmm4 dword to 1's */ \
			"pcmpgtd	%%xmm2, %%xmm5 \n\t"		/* signed comparison: 0 > root2? if so, set xmm5 dword to 1's */ \
			"movdqa %%xmm0, %%xmm6 \n\t"			/* copy of prime for neg root calculation */ \
			"movdqa %%xmm0, %%xmm7 \n\t"			/* copy of prime for neg root calculation */ \
			"movdqa %%xmm0, %%xmm8 \n\t"			/* copy of prime for neg root loops */ \
			"pand	%%xmm0, %%xmm4 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"pand	%%xmm0, %%xmm5 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"paddd	%%xmm4, %%xmm1 \n\t"			/* selectively add back prime (modular subtract) */ \
			"paddd	%%xmm5, %%xmm2 \n\t"			/* selectively add back prime (modular subtract) */ \
			"movdqa %%xmm1, %%xmm9 \n\t"					/* xmm5 = root1 copy */	\
			"pminud	%%xmm2, %%xmm1 \n\t"					/* xmm2 = root2 < root1 ? root2 : root1 */	\
			"pmaxud	%%xmm9, %%xmm2 \n\t"					/* xmm5 = root2 > root1 ? root2 : root1 */	\
			"movdqa %%xmm1, (%%r14,%%r15,4) \n\t"	/* save new root1 values */ \
			"movdqa %%xmm2, (%%r13,%%r15,4) \n\t"	/* save new root2 values */ \
			"psubd	%%xmm1, %%xmm6 \n\t"			/* form negative root1's; prime - root1 */ \
			"psubd	%%xmm2, %%xmm7 \n\t"			/* form negative root2's; prime - root2 */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
			"movl	88(%%rsi,1),%%r13d \n\t"		/* interval */ \
			"movl	96(%%rsi,1),%%r12d \n\t"		/* store bound_val in a register */ \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,1 ========== */	\
				/* ================================================ */	\
			"pextrd	$0,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$0,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"pextrd	$0,%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
			UPDATE_ROOT1_LOOP("0") \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,1 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2_LOOP("0") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,2 ========== */	\
				/* ================================================ */	\
			"pextrd	$1,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$1,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"pextrd	$1,%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
			UPDATE_ROOT1_LOOP("1") \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,2 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2_LOOP("1") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,3 ========== */	\
				/* ================================================ */	\
			"pextrd	$2,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$2,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"pextrd	$2,%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
			UPDATE_ROOT1_LOOP("2") \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,3 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2_LOOP("2") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,4 ========== */	\
				/* ================================================ */	\
			"pextrd	$3,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$3,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"pextrd	$3,%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
			UPDATE_ROOT1_LOOP("3") \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,4 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2_LOOP("3") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,1 ========== */	\
				/* ================================================ */	\
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"pextrd	$0,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$0,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"pextrd	$0,%%xmm8,%%edx \n\t"				/* else, extract prime from xmm0 */ \
			UPDATE_ROOT2_LOOP("0") \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,1 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1_LOOP("0") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,2 ========== */	\
				/* ================================================ */	\
			"pextrd	$1,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$1,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"pextrd	$1,%%xmm8,%%edx \n\t"				/* else, extract prime from xmm0 */ \
			UPDATE_ROOT2_LOOP("1") \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,2 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1_LOOP("1") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,3 ========== */	\
				/* ================================================ */	\
			"pextrd	$2,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$2,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"pextrd	$2,%%xmm8,%%edx \n\t"				/* else, extract prime from xmm0 */ \
			UPDATE_ROOT2_LOOP("2") \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,3 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1_LOOP("2") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,4 ========== */	\
				/* ================================================ */	\
			"pextrd	$3,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$3,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"pextrd	$3,%%xmm8,%%edx \n\t"				/* else, extract prime from xmm0 */ \
			UPDATE_ROOT2_LOOP("3") \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,4 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1_LOOP("3") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* ======== END OF LOOP - UPDATE AND CHECK ======== */	\
				/* ================================================ */	\
			"addl   $4,%%r15d \n\t"					/* increment j by 4*/ \
			"cmpl   84(%%rsi,1),%%r15d \n\t"		/* j < bound ? */ \
			"jb     8b \n\t"	\
			"9:		\n\t"				\
			"movl	%%r15d, %%eax \n\t" \
			:  \
			: "g"(&helperstruct) \
			: "rax", "rbx", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", 
				"xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm9", "memory", "cc");

		// refresh local pointers and constants before entering the next loop
		numptr_n = helperstruct.numptr_n;
		numptr_p = helperstruct.numptr_p;
		sliceptr_n = helperstruct.sliceptr_n;
		sliceptr_p = helperstruct.sliceptr_p;
		ptr = helperstruct.ptr;	
		bound_val = helperstruct.bound_val;	
		check_bound = helperstruct.check_bound;
		bound_index = helperstruct.bound_index;
		logp = helperstruct.logp;


#elif defined(HAS_SSE2)

		logp = update_data.logp[j-1];
		for (j=med_B;j<large_B; )
		{
			CHECK_NEW_SLICE(j);

			COMPUTE_4_PROOTS(j);

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);

			j++;
		}


#else
		logp = update_data.logp[j-1];
		for (j=med_B;j<large_B;j++,ptr++)
		{
			CHECK_NEW_SLICE(j);

			prime = update_data.prime[j];
			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];

			COMPUTE_NEXT_ROOTS_P;
			//fb_offset = (j - bound_val) << 16;

			update_data.firstroots1[j] = root1;
			update_data.firstroots2[j] = root2;

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);
		}

#endif

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		POLY_STG3 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);

		gettimeofday(&qs_timing_start, NULL);
#endif
			
		
#if defined(USE_POLY_SSE2_ASM) && defined(GCC_ASM64X) && !defined(PROFILING)
		logp = update_data.logp[large_B-1];

		if (large_B % 16 != 0)
		{
			printf("large_B must be divisible by 16!\n");
			exit(-1);
		}
		if ((bound - large_B) % 16 != 0)
		{
			printf("large range must be divisible by 16!\n");
			exit(-1);
		}

		//load up our helper struct, so we don't have a list
		//of a 100 things to stick into our asm block
		helperstruct.numptr_n = numptr_n;			//0
		helperstruct.numptr_p = numptr_p;			//8
		helperstruct.sliceptr_n = sliceptr_n;		//16
		helperstruct.sliceptr_p = sliceptr_p;		//24
		helperstruct.update_data_prime = update_data.prime;	//32
		helperstruct.update_data_root1 = update_data.firstroots1;	//40
		helperstruct.update_data_root2 = update_data.firstroots2;	//48
		helperstruct.update_data_logp = update_data.logp;	//56
		helperstruct.lp_bucket_p = lp_bucket_p;		//64
		helperstruct.ptr = &rootupdates[(v-1) * bound];  //72;						//72
		helperstruct.large_B = large_B;				//80
		helperstruct.B = bound;						//84
		helperstruct.interval = interval;			//88
		helperstruct.numblocks = numblocks;			//92
		helperstruct.bound_val = bound_val;			//96
		helperstruct.bound_index = bound_index;		//100
		helperstruct.check_bound = check_bound;		//104
		helperstruct.logp = logp;					//108
		helperstruct.intervalm1 = interval-1;		//112

		ASM_G (		\
			"movq	%0,%%rsi \n\t"					/* move helperstruct into rsi */ \
			"movl   80(%%rsi,1),%%r15d \n\t"		/* large_B = j = r15d */ \
													/* do the loop comparison */ \
			"cmpl   84(%%rsi,1),%%r15d \n\t"		/* j >= bound ? */ \
			"jae    9f	\n\t"						/* jump to end of loop, if test fails */ \
			"8: \n\t"	\
				/* ================================================ */	\
				/* ========== BEGIN CHECK_NEW_SLICE BLOCK ========= */	\
				/* ================================================ */	\
			CHECK_NEW_SLICE_ASM	\
				/* ================================================ */	\
				/* ============ BEGIN - GET NEW ROOTS 1 =========== */	\
				/* ================================================ */	\
			"movq	72(%%rsi,1),%%rdi \n\t"			/* edi = ptr */ \
			"movq	40(%%rsi,1),%%r14 \n\t"			/* move updata_data.root1 pointer into r14 */ \
			"movdqa (%%rdi,%%r15,4), %%xmm3 \n\t"	/* xmm3 = next 4 values of rootupdates */ \
			"movq	48(%%rsi,1),%%r13 \n\t"			/* move updata_data.root2 pointer into r13 */ \
			"movdqa (%%r14,%%r15,4), %%xmm1 \n\t"	/* xmm1 = next 4 values of root1 */ \
			"psubd	%%xmm3, %%xmm1 \n\t"			/* root1 -= ptr */ \
			"movdqa (%%r13,%%r15,4), %%xmm2 \n\t"	/* xmm2 = next 4 values of root2 */ \
			"movq	32(%%rsi,1),%%r12 \n\t"			/* move updata_data.prime pointer into r12 */ \
			"psubd	%%xmm3, %%xmm2 \n\t"			/* root2 -= ptr */ \
			"pxor	%%xmm4, %%xmm4 \n\t"			/* zero xmm4 */ \
			"pxor	%%xmm5, %%xmm5 \n\t"			/* zero xmm5 */ \
			"movdqa (%%r12,%%r15,4), %%xmm0 \n\t"	/* xmm0 = next 4 primes */ \
			"pcmpgtd	%%xmm1, %%xmm4 \n\t"		/* signed comparison: 0 > root1? if so, set xmm4 dword to 1's */ \
			"pcmpgtd	%%xmm2, %%xmm5 \n\t"		/* signed comparison: 0 > root2? if so, set xmm5 dword to 1's */ \
			"pand	%%xmm0, %%xmm4 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"movdqa %%xmm0, %%xmm6 \n\t"			/* copy prime to compute neg roots */ \
			"pand	%%xmm0, %%xmm5 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"movdqa %%xmm0, %%xmm7 \n\t"			/* copy prime to compute neg roots */ \
			"paddd	%%xmm4, %%xmm1 \n\t"			/* selectively add back prime (modular subtract) */ \
			"paddd	%%xmm5, %%xmm2 \n\t"			/* selectively add back prime (modular subtract) */ \
			"movdqa %%xmm1, %%xmm8 \n\t"					/* xmm5 = root1 copy */	\
			"pminud	%%xmm2, %%xmm1 \n\t"					/* xmm2 = root2 < root1 ? root2 : root1 */	\
			"pmaxud	%%xmm8, %%xmm2 \n\t"					/* xmm5 = root2 > root1 ? root2 : root1 */	\
			"movdqa %%xmm1, (%%r14,%%r15,4) \n\t"	/* save new root1 values */ \
			"movdqa %%xmm2, (%%r13,%%r15,4) \n\t"	/* save new root2 values */ \
			"psubd	%%xmm1, %%xmm6 \n\t"			/* form negative root1's; prime - root1 */ \
			"psubd	%%xmm2, %%xmm7 \n\t"			/* form negative root2's; prime - root2 */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
			"movl	88(%%rsi,1),%%r13d \n\t"		/* interval */ \
			"movl	96(%%rsi,1),%%r12d \n\t"		/* store bound_val in a register */ \
				/* edx is free at this point... can it be used? */ \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,1 ========== */	\
				/* ================================================ */	\
			"pextrd	$0,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$0,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			UPDATE_ROOT1("0") \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,1 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2("0") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,2 ========== */	\
				/* ================================================ */	\
			"pextrd	$1,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$1,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			UPDATE_ROOT1("1") \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,2 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2("1") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,3 ========== */	\
				/* ================================================ */	\
			"pextrd	$2,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$2,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			UPDATE_ROOT1("2") \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,3 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2("2") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,4 ========== */	\
				/* ================================================ */	\
			"pextrd	$3,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$3,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			UPDATE_ROOT1("3") \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,4 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2("3") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,1 ========== */	\
				/* ================================================ */	\
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"pextrd	$0,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$0,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			UPDATE_ROOT2("0") \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,1 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1("0") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,2 ========== */	\
				/* ================================================ */	\
			"pextrd	$1,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$1,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			UPDATE_ROOT2("1") \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,2 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1("1") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,3 ========== */	\
				/* ================================================ */	\
			"pextrd	$2,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$2,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			UPDATE_ROOT2("2") \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,3 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1("2") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,4 ========== */	\
				/* ================================================ */	\
			"pextrd	$3,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$3,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			UPDATE_ROOT2("3") \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,4 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1("3") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* ======== END OF LOOP - UPDATE AND CHECK ======== */	\
				/* ================================================ */	\
			"addl   $4,%%r15d \n\t"					/* increment j by 1*/ \
			"cmpl   84(%%rsi,1),%%r15d \n\t"		/* j < bound ? */ \
			"jb     8b \n\t"	\
			"9:		\n\t"				\
			"movl	%%r15d, %%eax \n\t" \
			:  \
			: "g"(&helperstruct) \
			: "rax", "rbx", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", 
				"xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "memory", "cc");

		bound_index = helperstruct.bound_index;
		logp = helperstruct.logp;


#elif defined(HAS_SSE2)

		logp = update_data.logp[j-1];
		for (j=large_B;j<bound; )
		{
			CHECK_NEW_SLICE(j);

			COMPUTE_4_PROOTS(j);

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);
			root2 = (prime - root2);

			FILL_ONE_PRIME_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);
			root2 = (prime - root2);

			FILL_ONE_PRIME_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);
			root2 = (prime - root2);

			FILL_ONE_PRIME_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);
			root2 = (prime - root2);

			FILL_ONE_PRIME_N(j);
			
			j++;
		}


#else
		logp = update_data.logp[j-1];
		for (j=large_B;j<bound;j++,ptr++)				
		{				
			CHECK_NEW_SLICE(j);

			prime = update_data.prime[j];			
			root1 = update_data.firstroots1[j];	
			root2 = update_data.firstroots2[j];	

			COMPUTE_NEXT_ROOTS_P;		

			update_data.firstroots1[j] = root1;	
			update_data.firstroots2[j] = root2;	

			FILL_ONE_PRIME_P(j);	

			root1 = (prime - root1);		
			root2 = (prime - root2);	
			
			FILL_ONE_PRIME_N(j);
		}

#endif

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		POLY_STG4 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);
#endif

	}
	else
	{

#ifdef QS_TIMING
		gettimeofday(&qs_timing_start, NULL);
#endif

		for (j=startprime;j<sconf->sieve_small_fb_start;j++,ptr++)
		{
			prime = update_data.prime[j];
			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];

			COMPUTE_NEXT_ROOTS_N;

			//we don't sieve these, so ordering doesn't matter
			update_data.firstroots1[j] = root1;
			update_data.firstroots2[j] = root2;

			fb_p->root1[j] = (uint16)root1;
			fb_p->root2[j] = (uint16)root2;
			fb_n->root1[j] = (uint16)(prime - root2);
			fb_n->root2[j] = (uint16)(prime - root1);
			if (fb_n->root1[j] == prime)
				fb_n->root1[j] = 0;
			if (fb_n->root2[j] == prime)
				fb_n->root2[j] = 0;

		}

		// do one at a time up to the 10bit boundary, where
		// we can start doing things 8 at a time and be
		// sure we can use aligned moves (static_data_init).	
		for (j=sconf->sieve_small_fb_start; 
			j < sconf->factor_base->fb_10bit_B; j++, ptr++)
		{

			// as soon as we are aligned, use more efficient sse2 based methods...
			if ((j&7) == 0)
				break;

			prime = update_data.prime[j];
			root1 = (uint32)update_data.sm_firstroots1[j];
			root2 = (uint32)update_data.sm_firstroots2[j];

			COMPUTE_NEXT_ROOTS_N;

			if (root2 < root1)
			{
				update_data.sm_firstroots1[j] = (uint16)root2;
				update_data.sm_firstroots2[j] = (uint16)root1;

				fb_p->root1[j] = (uint16)root2;
				fb_p->root2[j] = (uint16)root1;
				fb_n->root1[j] = (uint16)(prime - root1);
				fb_n->root2[j] = (uint16)(prime - root2);
			}
			else
			{
				update_data.sm_firstroots1[j] = (uint16)root1;
				update_data.sm_firstroots2[j] = (uint16)root2;

				fb_p->root1[j] = (uint16)root1;
				fb_p->root2[j] = (uint16)root2;
				fb_n->root1[j] = (uint16)(prime - root2);
				fb_n->root2[j] = (uint16)(prime - root1);
			}
		}
		
#if defined(GCC_ASM64X) || defined(_MSC_VER) //NOTDEF //GCC_ASM64X
		// update 8 at a time using SSE2 and no branching		
		sm_ptr = &dconf->sm_rootupdates[(v-1) * bound];
		{
			small_update_t h;
			
			h.first_r1 = update_data.sm_firstroots1;		// 0
			h.first_r2 = update_data.sm_firstroots2;		// 8
			h.fbp1 = fb_p->root1;							// 16
			h.fbp2 = fb_p->root2;							// 24
			h.fbn1 = fb_n->root1;							// 32
			h.fbn2 = fb_n->root2;							// 40
			h.primes = fb_p->prime;							// 48
			h.updates = sm_ptr;								// 56
			h.start = j;									// 64
			h.stop = sconf->factor_base->fb_15bit_B;		// 68
			if ((h.stop - 8) > h.start)
				h.stop -= 8;

			COMPUTE_8X_SMALL_NROOTS;
			
			j = h.stop;
		}
		sm_ptr = &dconf->sm_rootupdates[(v-1) * bound + j];
		

#else
		ptr = &dconf->rootupdates[(v-1) * bound + sconf->factor_base->fb_10bit_B];
		for (j=sconf->factor_base->fb_10bit_B; j < sconf->factor_base->fb_15bit_B; j++, ptr++)
		{
			prime = update_data.prime[j];
			root1 = update_data.sm_firstroots1[j];
			root2 = update_data.sm_firstroots2[j];

			COMPUTE_NEXT_ROOTS_N;

			if (root2 < root1)
			{
				update_data.sm_firstroots1[j] = (uint16)root2;
				update_data.sm_firstroots2[j] = (uint16)root1;

				fb_p->root1[j] = (uint16)root2;
				fb_p->root2[j] = (uint16)root1;
				fb_n->root1[j] = (uint16)(prime - root1);
				fb_n->root2[j] = (uint16)(prime - root2);
			}
			else
			{
				update_data.sm_firstroots1[j] = (uint16)root1;
				update_data.sm_firstroots2[j] = (uint16)root2;

				fb_p->root1[j] = (uint16)root1;
				fb_p->root2[j] = (uint16)root2;
				fb_n->root1[j] = (uint16)(prime - root2);
				fb_n->root2[j] = (uint16)(prime - root1);
			}
		}
#endif		

		// assembly code may not get all the way to 15 bits since we 
		// do things in blocks of 8 there.  Make sure we are at the 15 bit
		// boundary before we switch to using update_data.firstroots1/2.
		// this should only run a few iterations, if any.
		ptr = &dconf->rootupdates[(v-1) * bound + j];
		for ( ; j < med_B; j++, ptr++)
		{
			prime = update_data.prime[j];
			root1 = (uint16)update_data.sm_firstroots1[j];
			root2 = (uint16)update_data.sm_firstroots2[j];

			if ((prime > 32768) && ((j&7) == 0))
				break;

			COMPUTE_NEXT_ROOTS_N;

			if (root2 < root1)
			{
				update_data.sm_firstroots1[j] = (uint16)root2;
				update_data.sm_firstroots2[j] = (uint16)root1;

				fb_p->root1[j] = (uint16)root2;
				fb_p->root2[j] = (uint16)root1;
				fb_n->root1[j] = (uint16)(prime - root1);
				fb_n->root2[j] = (uint16)(prime - root2);
			}
			else
			{
				update_data.sm_firstroots1[j] = (uint16)root1;
				update_data.sm_firstroots2[j] = (uint16)root2;

				fb_p->root1[j] = (uint16)root1;
				fb_p->root2[j] = (uint16)root2;
				fb_n->root1[j] = (uint16)(prime - root2);
				fb_n->root2[j] = (uint16)(prime - root1);
			}
		}	

		// continue one at a time once we exceed 15 bits, because the 8x SSE2
		// code has a hard time with unsigned 16 bit comparisons
#if defined(GCC_ASM64X) || (defined(_MSC_VER) && defined(_WIN64))

		// update 8 at a time using SSE2 and no branching		
		sm_ptr = &dconf->sm_rootupdates[(v-1) * bound];
		{
			small_update_t h;
			
			h.first_r1 = update_data.sm_firstroots1;		// 0
			h.first_r2 = update_data.sm_firstroots2;		// 8
			h.fbp1 = fb_p->root1;							// 16
			h.fbp2 = fb_p->root2;							// 24
			h.fbn1 = fb_n->root1;							// 32
			h.fbn2 = fb_n->root2;							// 40
			h.primes = fb_p->prime;							// 48
			h.updates = sm_ptr;								// 56
			h.start = j;									// 64
			h.stop = med_B;									// 68

			COMPUTE_8X_SMALL_NROOTS_SSE41;

			j = h.stop;
		}
		sm_ptr = &dconf->sm_rootupdates[(v-1) * bound + j];
		
		
#else

		for ( ; j < med_B; j++, ptr++)
		{
			prime = update_data.prime[j];
			root1 = update_data.sm_firstroots1[j];
			root2 = update_data.sm_firstroots2[j];

			COMPUTE_NEXT_ROOTS_N;

			if (root2 < root1)
			{
				update_data.sm_firstroots1[j] = (uint16)root2;
				update_data.sm_firstroots2[j] = (uint16)root1;

				fb_p->root1[j] = (uint16)root2;
				fb_p->root2[j] = (uint16)root1;
				fb_n->root1[j] = (uint16)(prime - root1);
				fb_n->root2[j] = (uint16)(prime - root2);
			}
			else
			{
				update_data.sm_firstroots1[j] = (uint16)root1;
				update_data.sm_firstroots2[j] = (uint16)root2;

				fb_p->root1[j] = (uint16)root1;
				fb_p->root2[j] = (uint16)root2;
				fb_n->root1[j] = (uint16)(prime - root2);
				fb_n->root2[j] = (uint16)(prime - root1);
			}
		}	
		

#endif

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		POLY_STG2 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);

		gettimeofday(&qs_timing_start, NULL);
#endif

		bound_index = 0;
		bound_val = med_B;
		check_bound = med_B + BUCKET_ALLOC/2;
		
		
#if defined(USE_POLY_SSE2_ASM) && defined(GCC_ASM64X) && !defined(PROFILING)
		logp = update_data.logp[med_B-1];

		if (med_B % 16 != 0)
		{
			printf("med_B must be divisible by 16!\n");
			exit(-1);
		}
		if ((large_B - med_B) % 16 != 0)
		{
			printf("med range must be divisible by 16!\n");
			exit(-1);
		}

		//load up our helper struct, so we don't have a list
		//of a 100 things to stick into our asm block
		helperstruct.numptr_n = numptr_n;			//0
		helperstruct.numptr_p = numptr_p;			//8
		helperstruct.sliceptr_n = sliceptr_n;		//16
		helperstruct.sliceptr_p = sliceptr_p;		//24
		helperstruct.update_data_prime = update_data.prime;	//32
		helperstruct.update_data_root1 = update_data.firstroots1;	//40
		helperstruct.update_data_root2 = update_data.firstroots2;	//48
		helperstruct.update_data_logp = update_data.logp;	//56
		helperstruct.lp_bucket_p = lp_bucket_p;		//64
		helperstruct.ptr = &rootupdates[(v-1) * bound];  //72;						//72
		helperstruct.large_B = med_B;				//80
		helperstruct.B = large_B;					//84
		helperstruct.interval = interval;			//88
		helperstruct.numblocks = numblocks;			//92
		helperstruct.bound_val = bound_val;			//96
		helperstruct.bound_index = bound_index;		//100
		helperstruct.check_bound = check_bound;		//104
		helperstruct.logp = logp;					//108
		helperstruct.intervalm1 = interval-1;		//112

		ASM_G (		\
			"movq	%0,%%rsi \n\t"					/* move helperstruct into rsi */ \
			"movl   80(%%rsi,1),%%r15d \n\t"		/* large_B = j = r15d */ \
													/* do the loop comparison */ \
			"cmpl   84(%%rsi,1),%%r15d \n\t"		/* j >= bound ? */ \
			"jae    9f	\n\t"						/* jump to end of loop, if test fails */ \
			"8: \n\t"	\
				/* ================================================ */	\
				/* ========== BEGIN CHECK_NEW_SLICE BLOCK ========= */	\
				/* ================================================ */	\
			CHECK_NEW_SLICE_ASM	\
				/* ================================================ */	\
				/* ============ BEGIN - GET NEW ROOTS 1 =========== */	\
				/* ================================================ */	\
			"movq	72(%%rsi,1),%%rdi \n\t"			/* edi = ptr */ \
			"movq	40(%%rsi,1),%%r14 \n\t"			/* move updata_data.root1 pointer into r14 */ \
			"movdqa (%%rdi,%%r15,4), %%xmm3 \n\t"	/* xmm3 = next 4 values of rootupdates */ \
			"movq	48(%%rsi,1),%%r13 \n\t"			/* move updata_data.root2 pointer into r13 */ \
			"movdqa (%%r14,%%r15,4), %%xmm1 \n\t"	/* xmm1 = next 4 values of root1 */ \
			"paddd	%%xmm3, %%xmm1 \n\t"			/* root1 += ptr */ \
			"movdqa (%%r13,%%r15,4), %%xmm2 \n\t"	/* xmm2 = next 4 values of root2 */ \
			"movq	32(%%rsi,1),%%r12 \n\t"			/* move updata_data.prime pointer into r12 */ \
			"paddd	%%xmm3, %%xmm2 \n\t"			/* root2 += ptr */ \
			"movdqa	%%xmm1, %%xmm4 \n\t"			/* copy root1 to xmm4 */ \
			"movdqa (%%r12,%%r15,4), %%xmm0 \n\t"	/* xmm0 = next 4 primes */ \
			"movdqa	%%xmm2, %%xmm5 \n\t"			/* copy root2 to xmm5 */ \
			"pcmpgtd	%%xmm0, %%xmm4 \n\t"		/* signed comparison: root1 > p? if so, set xmm4 dword to 1's */ \
			"pcmpgtd	%%xmm0, %%xmm5 \n\t"		/* signed comparison: root2 > p? if so, set xmm5 dword to 1's */ \
			"movdqa %%xmm0, %%xmm6 \n\t"			/* copy of prime for neg root calculation */ \
			"movdqa %%xmm0, %%xmm7 \n\t"			/* copy of prime for neg root calculation */ \
			"movdqa %%xmm0, %%xmm8 \n\t"			/* copy of prime for neg root loops */ \
			"pand	%%xmm0, %%xmm4 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"pand	%%xmm0, %%xmm5 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"psubd	%%xmm4, %%xmm1 \n\t"			/* selectively sub back prime (modular addition) */ \
			"psubd	%%xmm5, %%xmm2 \n\t"			/* selectively sub back prime (modular addition) */ \
			"movdqa %%xmm1, %%xmm9 \n\t"					/* xmm5 = root1 copy */	\
			"pminud	%%xmm2, %%xmm1 \n\t"					/* xmm2 = root2 < root1 ? root2 : root1 */	\
			"pmaxud	%%xmm9, %%xmm2 \n\t"					/* xmm5 = root2 > root1 ? root2 : root1 */	\
			"movdqa %%xmm1, (%%r14,%%r15,4) \n\t"	/* save new root1 values */ \
			"movdqa %%xmm2, (%%r13,%%r15,4) \n\t"	/* save new root2 values */ \
			"psubd	%%xmm1, %%xmm6 \n\t"			/* form negative root1's; prime - root1 */ \
			"psubd	%%xmm2, %%xmm7 \n\t"			/* form negative root2's; prime - root2 */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
			"movl	88(%%rsi,1),%%r13d \n\t"		/* interval */ \
			"movl	96(%%rsi,1),%%r12d \n\t"		/* store bound_val in a register */ \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,1 ========== */	\
				/* ================================================ */	\
			"pextrd	$0,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$0,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"pextrd	$0,%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
			UPDATE_ROOT1_LOOP("0") \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,1 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2_LOOP("0") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,2 ========== */	\
				/* ================================================ */	\
			"pextrd	$1,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$1,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"pextrd	$1,%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
			UPDATE_ROOT1_LOOP("1") \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,2 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2_LOOP("1") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,3 ========== */	\
				/* ================================================ */	\
			"pextrd	$2,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$2,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"pextrd	$2,%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
			UPDATE_ROOT1_LOOP("2") \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,3 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2_LOOP("2") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,4 ========== */	\
				/* ================================================ */	\
			"pextrd	$3,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$3,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"pextrd	$3,%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
			UPDATE_ROOT1_LOOP("3") \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,4 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2_LOOP("3") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,1 ========== */	\
				/* ================================================ */	\
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"pextrd	$0,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$0,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"pextrd	$0,%%xmm8,%%edx \n\t"				/* else, extract prime from xmm0 */ \
			UPDATE_ROOT2_LOOP("0") \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,1 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1_LOOP("0") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,2 ========== */	\
				/* ================================================ */	\
			"pextrd	$1,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$1,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"pextrd	$1,%%xmm8,%%edx \n\t"				/* else, extract prime from xmm0 */ \
			UPDATE_ROOT2_LOOP("1") \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,2 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1_LOOP("1") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,3 ========== */	\
				/* ================================================ */	\
			"pextrd	$2,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$2,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"pextrd	$2,%%xmm8,%%edx \n\t"				/* else, extract prime from xmm0 */ \
			UPDATE_ROOT2_LOOP("2") \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,3 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1_LOOP("2") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,4 ========== */	\
				/* ================================================ */	\
			"pextrd	$3,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$3,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"pextrd	$3,%%xmm8,%%edx \n\t"				/* else, extract prime from xmm0 */ \
			UPDATE_ROOT2_LOOP("3") \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,4 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1_LOOP("3") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* ======== END OF LOOP - UPDATE AND CHECK ======== */	\
				/* ================================================ */	\
			"addl   $4,%%r15d \n\t"					/* increment j by 1*/ \
			"cmpl   84(%%rsi,1),%%r15d \n\t"		/* j < bound ? */ \
			"jb     8b \n\t"	\
			"9:		\n\t"				\
			"movl	%%r15d, %%eax \n\t" \
			:  \
			: "g"(&helperstruct) \
			: "rax", "rbx", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", 
				"xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm9", "memory", "cc");

		// refresh local pointers and constants before entering the next loop
		numptr_n = helperstruct.numptr_n;
		numptr_p = helperstruct.numptr_p;
		sliceptr_n = helperstruct.sliceptr_n;
		sliceptr_p = helperstruct.sliceptr_p;
		ptr = helperstruct.ptr;	
		bound_val = helperstruct.bound_val;	
		check_bound = helperstruct.check_bound;
		bound_index = helperstruct.bound_index;
		logp = helperstruct.logp;

#elif defined(HAS_SSE2)

		logp = update_data.logp[j-1];
		for (j=med_B;j<large_B; )
		{
			CHECK_NEW_SLICE(j);

			COMPUTE_4_NROOTS(j);

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);

			j++;
		}


#else

		logp = update_data.logp[j-1];
		for (j=med_B;j<large_B;j++,ptr++)
		{
			CHECK_NEW_SLICE(j);

			prime = update_data.prime[j];
			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];

			COMPUTE_NEXT_ROOTS_N;

			update_data.firstroots1[j] = root1;
			update_data.firstroots2[j] = root2;

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);
		}

#endif

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		POLY_STG3 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);

		gettimeofday(&qs_timing_start, NULL);
#endif

		
#if defined(USE_POLY_SSE2_ASM) && defined(GCC_ASM64X) && !defined(PROFILING)
		logp = update_data.logp[large_B-1];
		
		if (large_B % 16 != 0)
		{
			printf("large_B must be divisible by 16!\n");
			exit(-1);
		}
		if ((bound - large_B) % 16 != 0)
		{
			printf("large range must be divisible by 16!\n");
			exit(-1);
		}

		//load up our helper struct, so we don't have a list
		//of a 100 things to stick into our asm block
		helperstruct.numptr_n = numptr_n;			//0
		helperstruct.numptr_p = numptr_p;			//8
		helperstruct.sliceptr_n = sliceptr_n;		//16
		helperstruct.sliceptr_p = sliceptr_p;		//24
		helperstruct.update_data_prime = update_data.prime;	//32
		helperstruct.update_data_root1 = update_data.firstroots1;	//40
		helperstruct.update_data_root2 = update_data.firstroots2;	//48
		helperstruct.update_data_logp = update_data.logp;	//56
		helperstruct.lp_bucket_p = lp_bucket_p;		//64
		helperstruct.ptr = &rootupdates[(v-1) * bound];  //72;						//72
		helperstruct.large_B = large_B;				//80
		helperstruct.B = bound;						//84
		helperstruct.interval = interval;			//88
		helperstruct.numblocks = numblocks;			//92
		helperstruct.bound_val = bound_val;			//96
		helperstruct.bound_index = bound_index;		//100
		helperstruct.check_bound = check_bound;		//104
		helperstruct.logp = logp;					//108
		helperstruct.intervalm1 = interval-1;		//112

		ASM_G (		\
			"movq	%0,%%rsi \n\t"					/* move helperstruct into rsi */ \
			"movl   80(%%rsi,1),%%r15d \n\t"		/* large_B = j = r15d */ \
													/* do the loop comparison */ \
			"cmpl   84(%%rsi,1),%%r15d \n\t"		/* j >= bound ? */ \
			"jae    9f	\n\t"						/* jump to end of loop, if test fails */ \
			"8: \n\t"	\
				/* ================================================ */	\
				/* ========== BEGIN CHECK_NEW_SLICE BLOCK ========= */	\
				/* ================================================ */	\
			CHECK_NEW_SLICE_ASM	\
				/* ================================================ */	\
				/* ============ BEGIN - GET NEW ROOTS 1 =========== */	\
				/* ================================================ */	\
			"movq	72(%%rsi,1),%%rdi \n\t"			/* edi = ptr */ \
			"movq	40(%%rsi,1),%%r14 \n\t"			/* move updata_data.root1 pointer into r14 */ \
			"movdqa (%%rdi,%%r15,4), %%xmm3 \n\t"	/* xmm3 = next 4 values of rootupdates */ \
			"movq	48(%%rsi,1),%%r13 \n\t"			/* move updata_data.root2 pointer into r13 */ \
			"movdqa (%%r14,%%r15,4), %%xmm1 \n\t"	/* xmm1 = next 4 values of root1 */ \
			"paddd	%%xmm3, %%xmm1 \n\t"			/* root1 += ptr */ \
			"movdqa (%%r13,%%r15,4), %%xmm2 \n\t"	/* xmm2 = next 4 values of root2 */ \
			"movq	32(%%rsi,1),%%r12 \n\t"			/* move updata_data.prime pointer into r12 */ \
			"paddd	%%xmm3, %%xmm2 \n\t"			/* root2 += ptr */ \
			"movdqa	%%xmm1, %%xmm4 \n\t"			/* copy root1 to xmm4 */ \
			"movdqa (%%r12,%%r15,4), %%xmm0 \n\t"	/* xmm0 = next 4 primes */ \
			"movdqa	%%xmm2, %%xmm5 \n\t"			/* copy root2 to xmm5 */ \
			"pcmpgtd	%%xmm0, %%xmm4 \n\t"		/* signed comparison: root1 > p? if so, set xmm4 dword to 1's */ \
			"pcmpgtd	%%xmm0, %%xmm5 \n\t"		/* signed comparison: root2 > p? if so, set xmm5 dword to 1's */ \
			"movdqa %%xmm0, %%xmm6 \n\t"			/* copy prime to xmm6 */ \
			"pand	%%xmm0, %%xmm4 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"movdqa %%xmm0, %%xmm7 \n\t"			/* copy prime to xmm7 */ \
			"pand	%%xmm0, %%xmm5 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"psubd	%%xmm4, %%xmm1 \n\t"			/* selectively sub back prime (modular addition) */ \
			"psubd	%%xmm5, %%xmm2 \n\t"			/* selectively sub back prime (modular addition) */ \
			"movdqa %%xmm1, %%xmm8 \n\t"					/* xmm5 = root1 copy */	\
			"pminud	%%xmm2, %%xmm1 \n\t"					/* xmm2 = root2 < root1 ? root2 : root1 */	\
			"pmaxud	%%xmm8, %%xmm2 \n\t"					/* xmm5 = root2 > root1 ? root2 : root1 */	\
			"movdqa %%xmm1, (%%r14,%%r15,4) \n\t"	/* save new root1 values */ \
			"movdqa %%xmm2, (%%r13,%%r15,4) \n\t"	/* save new root2 values */ \
			"psubd	%%xmm1, %%xmm6 \n\t"			/* form negative root1's; prime - root1 */ \
			"psubd	%%xmm2, %%xmm7 \n\t"			/* form negative root2's; prime - root2 */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
			"movl	88(%%rsi,1),%%r13d \n\t"		/* interval */ \
			"movl	96(%%rsi,1),%%r12d \n\t"		/* store bound_val in a register */ \
				/* edx is free at this point... can it be used? */ \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,1 ========== */	\
				/* ================================================ */	\
			"pextrd	$0,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval?, if so then so is root2 */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$0,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			UPDATE_ROOT1("0") \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,1 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2("0") \
			"2: \n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,2 ========== */	\
				/* ================================================ */	\
			"pextrd	$1,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$1,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			UPDATE_ROOT1("1") \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,2 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2("1") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,3 ========== */	\
				/* ================================================ */	\
			"pextrd	$2,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$2,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			UPDATE_ROOT1("2") \
			
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,3 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2("2") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,4 ========== */	\
				/* ================================================ */	\
			"pextrd	$3,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$3,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			UPDATE_ROOT1("3") \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,4 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2("3") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,1 ========== */	\
				/* ================================================ */	\
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"pextrd	$0,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? , if so then so is root 1*/ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$0,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			UPDATE_ROOT2("0") \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,1 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1("0") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,2 ========== */	\
				/* ================================================ */	\
			"pextrd	$1,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? , if so then so is root 1*/ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$1,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			UPDATE_ROOT2("1") \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,2 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1("1") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,3 ========== */	\
				/* ================================================ */	\
			"pextrd	$2,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? , if so then so is root 1*/ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$2,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			UPDATE_ROOT2("2") \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,3 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1("2") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,4 ========== */	\
				/* ================================================ */	\
			"pextrd	$3,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? , if so then so is root 1*/ \
			"jae    2f \n\t" 						/* jump if CF = 1 */ \
			"pextrd	$3,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			UPDATE_ROOT2("3") \
			
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,4 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1("3") \
			"2:		\n\t" \
			"1:		\n\t" \
				/* ================================================ */	\
				/* ======== END OF LOOP - UPDATE AND CHECK ======== */	\
				/* ================================================ */	\
			"addl   $4,%%r15d \n\t"					/* increment j by 1*/ \
			"cmpl   84(%%rsi,1),%%r15d \n\t"		/* j < bound ? */ \
			"jb     8b \n\t"	\
			"9:		\n\t"				\
			"movl	%%r15d, %%eax \n\t" \
			:  \
			: "g"(&helperstruct) \
			: "rax", "rbx", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", 
				"xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "memory", "cc");

		bound_index = helperstruct.bound_index;
		logp = helperstruct.logp;

#elif defined(HAS_SSE2)

		logp = update_data.logp[j-1];
		for (j=large_B;j<bound; )
		{
			CHECK_NEW_SLICE(j);

			COMPUTE_4_NROOTS(j);

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);
			root2 = (prime - root2);

			FILL_ONE_PRIME_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);
			root2 = (prime - root2);

			FILL_ONE_PRIME_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);
			root2 = (prime - root2);

			FILL_ONE_PRIME_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);
			root2 = (prime - root2);

			FILL_ONE_PRIME_N(j);
			
			j++;
		}


#else

		logp = update_data.logp[j-1];
		for (j=large_B;j<bound;j++,ptr++)				
		{				
			CHECK_NEW_SLICE(j);

			prime = update_data.prime[j];			
			root1 = update_data.firstroots1[j];	
			root2 = update_data.firstroots2[j];	

			COMPUTE_NEXT_ROOTS_N;		

			update_data.firstroots1[j] = root1;	
			update_data.firstroots2[j] = root2;	

			FILL_ONE_PRIME_P(j);	

			root1 = (prime - root1);		
			root2 = (prime - root2);	
			
			FILL_ONE_PRIME_N(j);
		}

#endif

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		POLY_STG4 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);
#endif

	}

	if (lp_bucket_p->list != NULL)
	{
		lp_bucket_p->num_slices = bound_index + 1;
		lp_bucket_p->logp[bound_index] = logp;
	}

	return;
}

#endif // USE_SSE41
