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

#include "common.h"

#include "yafu.h"
#include "qs.h"
#include "util.h"
#include "poly_macros_32k.h"
#include "poly_macros_common.h"
#include "poly_macros_common_avx2.h"

// protect avx2 code under MSVC builds.  USE_AVX2 should be manually
// enabled at the top of qs.h for MSVC builds on supported hardware
#if defined( USE_AVX2 ) && defined(GCC_ASM64X)

// we have put (j - bound_val + i) | (root1 & 0x7fff) into ymm10.
#define UPDATE_ROOT1_NEW(it) \
		"movl   (%%r10,%%r8,4),%%r14d \n\t"	    /* numptr_p[bnum] */ \
        "incl   (%%r10,%%r8,4) \n\t"		    /* store new numptr to memory */ \
        "shll   $" BUCKET_BITStxt ",%%r8d \n\t"	/* bnum << BUCKET_BITS */ \
        "vpextrd $" it ",%%xmm10,%%edi \n\t"	/* load this lanes (j - bound_val + i) | (root1 & 0x7fff) */ \
		"addl   %%r14d,%%r8d \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
        "movl   %%edi,(%%r11,%%r8,4) \n\t"		/* store new fb_index/loc to memory */

#define UPDATE_ROOT2_NEW(it) \
		"movl   (%%r10,%%r9,4),%%r14d \n\t"	    /* numptr_p[bnum] */ \
        "incl   (%%r10,%%r9,4) \n\t"		    /* store new numptr to memory */ \
        "shll   $" BUCKET_BITStxt ",%%r9d \n\t"	/* bnum << BUCKET_BITS */ \
        "vpextrd $" it ",%%xmm11,%%edi \n\t"	/* load this lanes (j - bound_val + i) | (root1 & 0x7fff) */ \
		"addl   %%r14d,%%r9d \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
        "movl   %%edi,(%%r11,%%r9,4) \n\t"		/* store new fb_index/loc to memory */


// macro for iteratively adding an element (specified by r8d) to the end of a bucket list
#define UPDATE_ROOT1_LOOP_NEW(it) \
        "vpextrd $" it ",%%xmm10,%%edi \n\t"	/* load this lanes (j - bound_val + i) | (root1 & 0x7fff) */ \
		"1:		\n\t"	 						/* beginning of loop */ \
		"movl   %%r8d,%%ebx \n\t"				/* copy root1 to ebx */ \
        "movl   %%r8d,%%eax \n\t"				/* copy root1 to eax */ \
		"shrl   $15,%%ebx \n\t"	                /* right shift root1 by 15  = bnum */ \
        "andl   $32767,%%eax \n\t"	            /* root1 & BLOCKSIZEm1 */ \
		"movl   %%ebx,%%ecx \n\t"				/* copy bnum to start making address */ \
		"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
		"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
        "orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
		"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
		"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
		"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
		"addl	%%edx,%%r8d \n\t"				/* increment root by prime */ \
        "andl   $0xffff0000, %%edi \n\t"        /* clear low half to re-use high half */ \
		"cmpl   %%r13d,%%r8d \n\t"				/* root > interval? */ \
		"jb		1b \n\t"						/* repeat if necessary */

// macro for iteratively adding an element (specified by r9d) to the end of a bucket list
#define UPDATE_ROOT2_LOOP_NEW(it) \
		"vpextrd $" it ",%%xmm10,%%edi \n\t"	/* load this lanes (j - bound_val + i) | (root1 & 0x7fff) */ \
		"1:		\n\t"	 						/* beginning of loop */ \
		"movl   %%r9d,%%ebx \n\t"				/* copy root2 to ebx */ \
        "movl   %%r9d,%%eax \n\t"				/* copy root2 to eax */ \
		"shrl   $15,%%ebx \n\t"	                /* right shift root2 by 15  = bnum */ \
        "andl   $32767,%%eax \n\t"	            /* root2 & BLOCKSIZEm1 */ \
		"movl   %%ebx,%%ecx \n\t"				/* copy bnum to start making address */ \
		"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
		"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
        "orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
		"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
		"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
		"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
		"addl	%%edx,%%r9d \n\t"				/* increment root by prime */ \
        "andl   $0xffff0000, %%edi \n\t"        /* clear low half to re-use high half */ \
		"cmpl   %%r13d,%%r9d \n\t"				/* root > interval? */ \
		"jb		1b \n\t"						/* repeat if necessary */
   

#define CHECK_NEW_SLICE_ASM_NEW \
		"cmpl   104(%%rsi,1),%%r15d	\n\t"		/* compare j with check_bound */ \
			/* note this is the counter j, not the byte offset j */ \
		"jge     1f \n\t"						/* jump into "if" code if comparison works */ \
			/* else, this is the "else-if" check */ \
		"movl   %%r15d,%%ebx \n\t"				/* copy j into ebx */ \
		"subl   96(%%rsi,1),%%ebx \n\t"			/* ebx = j - bound_val */ \
		"cmpl   $0xffff,%%ebx \n\t"				/* compare to 2^16 */ \
		"jbe    2f \n\t"						/* exit CHECK_NEW_SLICE if this comparison fails too */ \
			/* now we are in the else-if block of CHECK_NEW_SLICE */ \
		"xorq	%%rdx, %%rdx \n\t"				/* clear rdx */ \
		"movl   100(%%rsi,1),%%edx \n\t"		/* move bound_index into rdx */ \
		"movq   64(%%rsi,1),%%r9 \n\t"			/* move lp_bucket_p ptr into r9 */ \
		"movq	16(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->logp ptr into r8 */ \
		"movq	56(%%rsi,1),%%r14 \n\t"			/* move updata_data.logp pointer into r14 */ \
		"movzbl (%%r14,%%r15,1),%%ebx \n\t"		/* bring in logp */ \
		"movb	%%bl, 108(%%rsi,1) \n\t"		/* shove logp into output */ \
		"movb   %%bl,(%%r8,%%rdx,1) \n\t"		/* mov logp into lp_bucket_p->logp[bound_index] */ \
		"incq   %%rdx \n\t"						/* increment bound_index locally */ \
		"movl   %%edx,100(%%rsi,1) \n\t"		/* copy bound_index back to structure */ \
		"movq	8(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->fb_bounds ptr into r8 */ \
		"movl   %%r15d,(%%r8,%%rdx,4) \n\t"		/* mov j into lp_bucket_p->fb_bounds[bound_index] */ \
			/* note this is the counter j, not the byte offset j */ \
		"movl   %%r15d,96(%%rsi,1) \n\t"		/* bound_val = j */ \
            "movq	120(%%rsi,1),%%rdi \n\t"			/* edi = polyscratch */ \
            "vmovdqa     (%%rdi), %%ymm15 \n\t" /* (j - bound_val) == 0, so just move in the lane offset */ \
            "vpslld	$16, %%ymm15, %%ymm15 \n\t"	/* put (j - bound_val + i) into high half of 32-bit words */ \
		"xorq	%%rbx, %%rbx \n\t"				/* clear rbx */ \
		"movl   92(%%rsi,1),%%ebx \n\t"			/* put numblocks into ebx */ \
		"shll	$2,%%ebx \n\t"					/* numblocks * 4 (translate to bytes) */ \
		"shll	$1,%%ebx \n\t"					/* numblocks << 1 (negative blocks are contiguous) */ \
		"addq   %%rbx,8(%%rsi,1) \n\t"			/* numptr_p += (numblocks << 1) */ \
		"addq   %%rbx,0(%%rsi,1) \n\t"			/* numptr_n += (numblocks << 1) */ \
		"shlq   $" BUCKET_BITStxt ",%%rbx \n\t"	/* numblocks << (BUCKET_BITS + 1) */ \
			/* note also, this works because we've already left shifted by 1 */ \
		"addq   %%rbx,24(%%rsi,1) \n\t"			/* sliceptr_p += (numblocks << 11) */ \
		"addq   %%rbx,16(%%rsi,1) \n\t"			/* sliceptr_n += (numblocks << 11) */ \
		"addl   $" HALFBUCKET_ALLOCtxt ",104(%%rsi,1) \n\t"		/* add 2^(BUCKET_BITS-1) to check_bound */ \
		"cmp	%%rax,%%rax \n\t"				/* force jump */ \
		"je		2f \n\t"						/* jump out of CHECK_NEW_SLICE */ \
		"1:		\n\t"									\
			/* now we are in the if block of CHECK_NEW_SLICE */ \
		"xorl   %%ecx,%%ecx \n\t"				/* ecx = room  = 0 */ \
		"xorq	%%rbx, %%rbx \n\t"				/* loop counter = 0 */ \
		"cmpl   92(%%rsi,1),%%ebx \n\t"			/* compare with numblocks */ \
		"jae    3f \n\t"						/* jump past loop if condition met */ \
			/* condition not met, put a couple things in registers */ \
		"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
		"movq	0(%%rsi,1),%%r11 \n\t"			/* numptr_n into r11 */ \
		"5:		\n\t"							\
			/* now we are in the room loop */ \
			/* room is in register ecx */ \
		"movl   (%%r10,%%rbx,4),%%eax \n\t"		/* value at numptr_p + k */ \
		"movl   (%%r11,%%rbx,4),%%edx \n\t"		/* value at numptr_n + k */ \
		"cmpl   %%ecx,%%eax \n\t"				/* *(numptr_p + k) > room ? */ \
		"cmova  %%eax,%%ecx \n\t"				/* new value of room if so */ \
		"cmpl   %%ecx,%%edx \n\t"				/* *(numptr_p + k) > room ? */ \
		"cmova  %%edx,%%ecx \n\t"				/* new value of room if so */ \
		"incq   %%rbx \n\t"						/* increment counter */ \
		"cmpl   92(%%rsi,1),%%ebx \n\t"			/* compare to numblocks */ \
		"jl     5b \n\t"						/* iterate loop if condition met */ \
		"3:		\n\t"							\
		"movl   $" BUCKET_ALLOCtxt ",%%ebx \n\t"	/* move bucket allocation into register for subtraction */ \
		"subl   %%ecx,%%ebx \n\t"				/* room = bucket_alloc - room */ \
		"cmpl   $31,%%ebx \n\t"					/* answer less than 32? */ \
		"movl   %%ebx,%%ecx \n\t"				/* copy answer back to room register */ \
		"jle    4f \n\t"						/* jump if less than */ \
		"sarl   %%ebx	\n\t"					/* room >> 1 (copy of room) */ \
		"addl   %%ebx,104(%%rsi,1) \n\t"		/* add (room >> 1) to check_bound */ \
		"cmpq	%%rax,%%rax \n\t"				/* force jump */ \
		"je     2f \n\t"						/* jump out of CHECK_NEW_SLICE */ \
		"4:		\n\t"							\
			/* now we are inside the (room < 2) block */ \
		"xorq	%%rax, %%rax \n\t" \
		"movl   %%r15d,%%eax \n\t"				/* copy j to scratch reg */ \
		"shll   $0x4,%%eax \n\t"				/* multiply by 16 bytes per j */ \
		"xorq	%%rdx, %%rdx \n\t" \
		"movl   100(%%rsi,1),%%edx \n\t"		/* move bound_index into rdx */ \
		"movq   64(%%rsi,1),%%r9 \n\t"			/* move lp_bucket_p ptr into r9 */ \
		"movq	16(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->logp ptr into r8 */ \
		"movq	56(%%rsi,1),%%r14 \n\t"			/* move updata_data.logp pointer into r14 */ \
		"movzbl (%%r14,%%r15,1),%%ebx \n\t"		/* bring in logp */ \
		"movb	%%bl, 108(%%rsi,1) \n\t"		/* shove logp into output */ \
		"movb   %%bl,(%%r8,%%rdx,1) \n\t"		/* mov logp into lp_bucket_p->logp[bound_index] */ \
		"incq   %%rdx \n\t"						/* increment bound_index locally */ \
		"movl   %%edx,100(%%rsi,1) \n\t"		/* copy bound_index back to structure */ \
		"movq	8(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->fb_bounds ptr into r8 */ \
		"movl   %%r15d,(%%r8,%%rdx,4) \n\t"		/* mov j into lp_bucket_p->fb_bounds[bound_index] */ \
			/* note this is the counter j, not the byte offset j */ \
		"movl   %%r15d,96(%%rsi,1) \n\t"		/* bound_val = j */ \
            "movq	120(%%rsi,1),%%rdi \n\t"			/* edi = polyscratch */ \
            "vmovdqa     (%%rdi), %%ymm15 \n\t" /* (j - bound_val) == 0, so just move in the lane offset */ \
            "vpslld	$16, %%ymm15, %%ymm15 \n\t"	/* put (j - bound_val + i) into high half of 32-bit words */ \
		"xorq	%%rbx, %%rbx \n\t" \
		"movl   92(%%rsi,1),%%ebx \n\t"			/* put numblocks into ebx */ \
		"shll	$2,%%ebx \n\t"					/* numblocks * 4 (bytes) */ \
		"shll	$1,%%ebx \n\t"					/* numblocks << 1 */ \
		"addq   %%rbx,8(%%rsi,1) \n\t"			/* numptr_p += (numblocks << 1) */ \
		"addq   %%rbx,0(%%rsi,1) \n\t"			/* numptr_n += (numblocks << 1) */ \
		"shll   $" BUCKET_BITStxt ",%%ebx \n\t"	/* numblocks << (BUCKET_BITS + 1) */ \
			/* note also, this works because we've already left shifted by 1 */ \
		"addq   %%rbx,24(%%rsi,1) \n\t"			/* sliceptr_p += (numblocks << 11) */ \
		"addq   %%rbx,16(%%rsi,1) \n\t"			/* sliceptr_n += (numblocks << 11) */ \
		"addl   $" HALFBUCKET_ALLOCtxt ",104(%%rsi,1) \n\t"		/* add 2^(BUCKET_BITS-1) to check_bound */ \
		"2:		\n\t"


//this is in the poly library, even though the bulk of the time is spent
//bucketizing large primes, because it's where the roots of a poly are updated
void nextRoots_32k_avx2(static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//update the roots 
	sieve_fb_compressed *fb_p = dconf->comp_sieve_p;
	sieve_fb_compressed *fb_n = dconf->comp_sieve_n;
	int *rootupdates = dconf->rootupdates;
    int it;

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

    CLEAN_AVX2;
	
	if (lp_bucket_p->alloc_slices != 0) //NULL)
	{
		lp_bucket_p->fb_bounds[0] = med_B;

		sliceptr_p = lp_bucket_p->list;
		sliceptr_n = lp_bucket_p->list + (numblocks << BUCKET_BITS);

		numptr_p = lp_bucket_p->num;
		numptr_n = lp_bucket_p->num + numblocks;
		
        // reset bucket counts
        for (j = 0; j < lp_bucket_p->list_size; j++)
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

    dconf->polyscratch[0] = 0;
    dconf->polyscratch[1] = 1;
    dconf->polyscratch[2] = 2;
    dconf->polyscratch[3] = 3;
    dconf->polyscratch[4] = 4;
    dconf->polyscratch[5] = 5;
    dconf->polyscratch[6] = 6;
    dconf->polyscratch[7] = 7;
    dconf->polyscratch[8] = 0x00007fff;
    dconf->polyscratch[9] = 0x00007fff;
    dconf->polyscratch[10] = 0x00007fff;
    dconf->polyscratch[11] = 0x00007fff;
    dconf->polyscratch[12] = 0x00007fff;
    dconf->polyscratch[13] = 0x00007fff;
    dconf->polyscratch[14] = 0x00007fff;
    dconf->polyscratch[15] = 0x00007fff;

	k=0;
	ptr = &rootupdates[(v-1) * bound + startprime];	

	if (sign > 0)
	{
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
		// we can start doing things 16 at a time and be
		// sure we can use aligned moves (static_data_init).		
		for (j=sconf->sieve_small_fb_start; 
			j < sconf->factor_base->fb_10bit_B; j++, ptr++)
		{
			// as soon as we are aligned, use more efficient AVX2 based methods...
			if ((j & 15) == 0)
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

		// update 16 at a time using AVX2
		sm_ptr = &dconf->sm_rootupdates[(v-1) * med_B];
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
			h.stop = sconf->factor_base->med_B;		// 68

			COMPUTE_16X_SMALL_PROOTS_AVX2;
			
			j = h.stop;
		}	

		bound_index = 0;
		bound_val = med_B;
		check_bound = med_B + BUCKET_ALLOC/2;
		
		logp = update_data.logp[med_B-1];

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
        helperstruct.filler = 42;
        helperstruct.scratch = dconf->polyscratch;        //120

		//if (VFLAG > 2)
		//{
		//	printf("commencing 1\n");
		//}

        ASM_G (		\
            "movq	%0,%%rsi \n\t"					/* move helperstruct into rsi */ \
            "movl   80(%%rsi,1),%%r15d \n\t"		/* large_B = j = r15d */ \
            "vmovd      %%r15d, %%xmm15 \n\t"           /* broadcast j to ymm15 */ \
            "vpshufd	    $0, %%xmm15, %%xmm15 \n\t"  /* broadcast j to ymm15 */ \
            "vinserti128    $1, %%xmm15, %%ymm15, %%ymm15 \n\t" \
            "movq	120(%%rsi,1),%%rdi \n\t"			/* edi = polyscratch */ \
            "vpaddd     (%%rdi), %%ymm15, %%ymm15 \n\t" /* add in lane offsets */ \
            "movl	    96(%%rsi,1),%%r13d \n\t"		/* store bound_val in a register */ \
            "vmovd      %%r13d, %%xmm14 \n\t"           /* broadcast bound_val to ymm14 */ \
            "vpshufd	    $0, %%xmm14, %%xmm14 \n\t"  /* broadcast bound_val to ymm14 */ \
            "vinserti128    $1, %%xmm14, %%ymm14, %%ymm14 \n\t" \
            "vpsubd         %%ymm14, %%ymm15, %%ymm15 \n\t" /* subtract bound_val */ \
            "vpslld	$16, %%ymm15, %%ymm15 \n\t"		/* put (j - bound_val + i) into high half of 32-bit words */ \
            "movl      $524288, %%edx \n\t"           /* make the shifted-j increment */ \
            "vmovd      %%edx, %%xmm13 \n\t"           /* make the shifted-j increment */ \
            "vpshufd	    $0, %%xmm13, %%xmm13 \n\t"  /* broadcast shifted-j to ymm13 */ \
            "vinserti128    $1, %%xmm13, %%ymm13, %%ymm13 \n\t" \
            /* do the loop comparison */ \
            "cmpl   84(%%rsi,1),%%r15d \n\t"		/* j >= bound ? */ \
            "jae    9f	\n\t"						/* jump to end of loop, if test fails */ \
            "8: \n\t"	\
            /* ================================================ */	\
            /* ========== BEGIN CHECK_NEW_SLICE BLOCK ========= */	\
            /* ================================================ */	\
            CHECK_NEW_SLICE_ASM_NEW	\
            /* ================================================ */	\
            /* ============ BEGIN - GET NEW ROOTS 1 =========== */	\
            /* ================================================ */	\
            "movq	    72(%%rsi,1),%%rdi \n\t"			/* edi = ptr */ \
            "movq	    40(%%rsi,1),%%r14 \n\t"			/* move updata_data.root1 pointer into r14 */ \
            "vmovdqa    (%%rdi,%%r15,4), %%ymm3 \n\t"	/* xmm3 = next 4 values of rootupdates */ \
            "movq	    48(%%rsi,1),%%r13 \n\t"			/* move updata_data.root2 pointer into r13 */ \
            "vmovdqa    (%%r14,%%r15,4), %%ymm1 \n\t"	/* xmm1 = next 4 values of root1 */ \
            "vpsubd	    %%ymm3, %%ymm1, %%ymm1 \n\t"			/* root1 -= ptr */ \
            "vmovdqa    (%%r13,%%r15,4), %%ymm2 \n\t"	/* xmm2 = next 4 values of root2 */ \
            "movq	    32(%%rsi,1),%%rdi \n\t"			/* move updata_data.prime pointer into r12 */ \
            "vpsubd	    %%ymm3, %%ymm2, %%ymm2 \n\t"			/* root2 -= ptr */ \
            "vpxor	    %%ymm4, %%ymm4, %%ymm4 \n\t"			/* zero xmm4 */ \
            "vpxor	    %%ymm5, %%ymm5, %%ymm5 \n\t"			/* zero xmm5 */ \
            "vmovdqa    (%%rdi,%%r15,4), %%ymm0 \n\t"	/* xmm0 = next 4 primes */ \
            "vpcmpgtd	%%ymm1, %%ymm4, %%ymm4 \n\t"		/* signed comparison: 0 > root1? if so, set xmm4 dword to 1's */ \
            "vpcmpgtd	%%ymm2, %%ymm5, %%ymm5 \n\t"		/* signed comparison: 0 > root2? if so, set xmm5 dword to 1's */ \
            "vmovdqa    %%ymm0, %%ymm6 \n\t"			/* copy of prime for neg root calculation */ \
            "vmovdqa    %%ymm0, %%ymm7 \n\t"			/* copy of prime for neg root calculation */ \
            "vmovdqa    %%ymm0, %%ymm8 \n\t"			/* copy of prime for neg root loops */ \
            "vpand	    %%ymm0, %%ymm4, %%ymm4 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
            "vpand	    %%ymm0, %%ymm5, %%ymm5 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
            "vpaddd	    %%ymm4, %%ymm1, %%ymm1 \n\t"			/* selectively add back prime (modular subtract) */ \
            "vpaddd	    %%ymm5, %%ymm2, %%ymm2 \n\t"			/* selectively add back prime (modular subtract) */ \
            "vmovdqa    %%ymm1, %%ymm9 \n\t"					/* xmm5 = root1 copy */	\
            "vpminud	%%ymm2, %%ymm1, %%ymm1 \n\t"					/* xmm2 = root2 < root1 ? root2 : root1 */	\
            "vpmaxud	%%ymm9, %%ymm2, %%ymm2 \n\t"					/* xmm5 = root2 > root1 ? root2 : root1 */	\
            "vmovdqa    %%ymm1, (%%r14,%%r15,4) \n\t"	/* save new root1 values */ \
            "vextracti128   $0,%%ymm15,%%xmm10 \n\t" \
            "vmovdqa    %%ymm2, (%%r13,%%r15,4) \n\t"	/* save new root2 values */ \
            "vpsubd	    %%ymm1, %%ymm6, %%ymm6 \n\t"			/* form negative root1's; prime - root1 */ \
            "vpsubd	    %%ymm2, %%ymm7, %%ymm7 \n\t"			/* form negative root2's; prime - root2 */ \
            "movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
            "movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
            "movl	88(%%rsi,1),%%r13d \n\t"		/* interval */ \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,1 ========== */	\
            /* ================================================ */	\
            "vpextrd	$0,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$0,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$0,%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT1_LOOP_NEW("0") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,1 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT2_LOOP_NEW("0") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,2 ========== */	\
            /* ================================================ */	\
            "vpextrd	$1,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$1,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$1,%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT1_LOOP_NEW("1") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,2 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT2_LOOP_NEW("1") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,3 ========== */	\
            /* ================================================ */	\
            "vpextrd	$2,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$2,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$2,%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT1_LOOP_NEW("2") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,3 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT2_LOOP_NEW("2") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,4 ========== */	\
            /* ================================================ */	\
            "vpextrd	$3,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$3,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$3,%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT1_LOOP_NEW("3") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,4 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT2_LOOP_NEW("3") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,5 ========== */	\
            /* ================================================ */	\
            "vextracti128   $1,%%ymm1,%%xmm1 \n\t" \
            "vextracti128   $1,%%ymm2,%%xmm2 \n\t" \
            "vextracti128   $1,%%ymm0,%%xmm0 \n\t" \
            "vextracti128   $1,%%ymm15,%%xmm10 \n\t" \
            "vpextrd	$0,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$0,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$0,%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT1_LOOP_NEW("4") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,5 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT2_LOOP_NEW("4") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,6 ========== */	\
            /* ================================================ */	\
            "vpextrd	$1,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$1,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$1,%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT1_LOOP_NEW("5") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,6 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT2_LOOP_NEW("5") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,7 ========== */	\
            /* ================================================ */	\
            "vpextrd	$2,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$2,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$2,%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "jae    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_LOOP_NEW("6") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,7 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT2_LOOP_NEW("6") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,8 ========== */	\
            /* ================================================ */	\
            "vpextrd	$3,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$3,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$3,%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT1_LOOP_NEW("7") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,8 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT2_LOOP_NEW("7") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,1 ========== */	\
            /* ================================================ */	\
            "movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
            "movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
            "vextracti128   $0,%%ymm15,%%xmm10 \n\t" \
            "vpextrd	$0,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$0,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$0,%%xmm8,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT2_LOOP_NEW("0") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,1 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT1_LOOP_NEW("0") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,2 ========== */	\
            /* ================================================ */	\
            "vpextrd	$1,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$1,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$1,%%xmm8,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT2_LOOP_NEW("1") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,2 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT1_LOOP_NEW("1") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,3 ========== */	\
            /* ================================================ */	\
            "vpextrd	$2,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$2,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$2,%%xmm8,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT2_LOOP_NEW("2") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,3 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT1_LOOP_NEW("2") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,4 ========== */	\
            /* ================================================ */	\
            "vpextrd	$3,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$3,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$3,%%xmm8,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT2_LOOP_NEW("3") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,4 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT1_LOOP_NEW("3") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,5 ========== */	\
            /* ================================================ */	\
            "vextracti128   $1,%%ymm7,%%xmm7 \n\t" \
            "vextracti128   $1,%%ymm6,%%xmm6 \n\t" \
            "vextracti128   $1,%%ymm8,%%xmm8 \n\t" \
            "vextracti128   $1,%%ymm15,%%xmm10 \n\t" \
            "vpextrd	$0,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$0,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$0,%%xmm8,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT2_LOOP_NEW("4") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,5 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT1_LOOP_NEW("4") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,6 ========== */	\
            /* ================================================ */	\
            "vpextrd	$1,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$1,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$1,%%xmm8,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT2_LOOP_NEW("5") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,6 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT1_LOOP_NEW("5") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,7 ========== */	\
            /* ================================================ */	\
            "vpextrd	$2,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$2,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$2,%%xmm8,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT2_LOOP_NEW("6") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,7 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT1_LOOP_NEW("6") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,8 ========== */	\
            /* ================================================ */	\
            "vpextrd	$3,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$3,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$3,%%xmm8,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT2_LOOP_NEW("7") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,8 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT1_LOOP_NEW("7") \
            /* ================================================ */	\
            /* ======== END OF LOOP - UPDATE AND CHECK ======== */	\
            /* ================================================ */	\
            "addl   $8,%%r15d \n\t"					/* increment j by 4*/ \
            "vpaddd  %%ymm13, %%ymm15, %%ymm15 \n\t" /* add in shifted-j: (j+8)<<16 */ \
            "cmpl   84(%%rsi,1),%%r15d \n\t"		/* j < bound ? */ \
            "jb     8b \n\t"	\
            "9:		\n\t"				\
            "movl	%%r15d, %%eax \n\t" \
            :  \
            : "g"(&helperstruct) \
            : "rax", "rbx", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11", "r13", "r14", "r15",
            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm10", "xmm11", "xmm13", "xmm14", "xmm15", "memory", "cc");


		// refresh local pointers and constants before entering the next loop
		numptr_n = helperstruct.numptr_n;
		numptr_p = helperstruct.numptr_p;
		sliceptr_n = helperstruct.sliceptr_n;
		sliceptr_p = helperstruct.sliceptr_p;
		ptr = helperstruct.ptr;	
		bound_val = helperstruct.bound_val;	
		check_bound = helperstruct.check_bound;
		bound_index = helperstruct.bound_index;
		logp = update_data.logp[large_B-1];

		//if (VFLAG > 2)
		//{
		//	printf("commencing 2\n");
		//}

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
        helperstruct.filler = 0;
        helperstruct.scratch = dconf->polyscratch;        //120

        ASM_G (		\
            "movq	%0,%%rsi \n\t"					/* move helperstruct into rsi */ \
            "movl   80(%%rsi,1),%%r15d \n\t"		/* large_B = j = r15d */ \
            "vmovd      %%r15d, %%xmm15 \n\t"           /* broadcast j to ymm15 */ \
            "vpshufd	    $0, %%xmm15, %%xmm15 \n\t"  /* broadcast j to ymm15 */ \
            "vinserti128    $1, %%xmm15, %%ymm15, %%ymm15 \n\t" \
            "movq	120(%%rsi,1),%%rdi \n\t"			/* edi = polyscratch */ \
            "vpaddd     (%%rdi), %%ymm15, %%ymm15 \n\t" /* add in lane offsets */ \
            "movl	    96(%%rsi,1),%%r13d \n\t"		/* store bound_val in a register */ \
            "vmovd      %%r13d, %%xmm14 \n\t"           /* broadcast bound_val to ymm14 */ \
            "vpshufd	    $0, %%xmm14, %%xmm14 \n\t"  /* broadcast bound_val to ymm14 */ \
            "vinserti128    $1, %%xmm14, %%ymm14, %%ymm14 \n\t" \
            "vpsubd         %%ymm14, %%ymm15, %%ymm15 \n\t" /* subtract bound_val */ \
            "vpslld	$16, %%ymm15, %%ymm15 \n\t"		/* put (j - bound_val + i) into high half of 32-bit words */ \
            "movl      $524288, %%edx \n\t"           /* make the shifted-j increment */ \
            "vmovd      %%edx, %%xmm13 \n\t"           /* make the shifted-j increment */ \
            "vpshufd	    $0, %%xmm13, %%xmm13 \n\t"  /* broadcast shifted-j to ymm13 */ \
            "vinserti128    $1, %%xmm13, %%ymm13, %%ymm13 \n\t" \
            "vmovdqa     32(%%rdi), %%ymm12 \n\t" /* put the root masks into ymm12 */ \
            /* do the loop comparison */ \
            "cmpl   84(%%rsi,1),%%r15d \n\t"		/* j >= bound ? */ \
            "jae    9f	\n\t"						/* jump to end of loop, if test fails */ \
            "8: \n\t"	\
            /* ================================================ */	\
            /* ========== BEGIN CHECK_NEW_SLICE BLOCK ========= */	\
            /* ================================================ */	\
            CHECK_NEW_SLICE_ASM_NEW	\
            /* ================================================ */	\
            /* ============ BEGIN - GET NEW ROOTS 1 =========== */	\
            /* ================================================ */	\
            "movq	72(%%rsi,1),%%rdi \n\t"			/* edi = ptr */ \
            "movq	40(%%rsi,1),%%r14 \n\t"			/* move updata_data.root1 pointer into r14 */ \
            "vmovdqa (%%rdi,%%r15,4), %%ymm3 \n\t"	/* xmm3 = next 4 values of rootupdates */ \
            "movq	48(%%rsi,1),%%r13 \n\t"			/* move updata_data.root2 pointer into r13 */ \
            "vmovdqa (%%r14,%%r15,4), %%ymm1 \n\t"	/* xmm1 = next 4 values of root1 */ \
            "vpxor	    %%ymm4, %%ymm4, %%ymm4 \n\t"			/* zero xmm4 */ \
            "vpxor	    %%ymm5, %%ymm5, %%ymm5 \n\t"			/* zero xmm5 */ \
            "vpsubd	    %%ymm3, %%ymm1, %%ymm1 \n\t"			/* root1 -= ptr */ \
            "vmovdqa (%%r13,%%r15,4), %%ymm2 \n\t"	/* xmm2 = next 4 values of root2 */ \
            "movq	32(%%rsi,1),%%rdi \n\t"			/* move updata_data.prime pointer into r12 */ \
            "vpsubd	    %%ymm3, %%ymm2, %%ymm2 \n\t"			/* root2 -= ptr */ \
            "vmovdqa    (%%rdi,%%r15,4), %%ymm0 \n\t"	/* xmm0 = next 4 primes */ \
            "vpcmpgtd	%%ymm1, %%ymm4, %%ymm4 \n\t"		/* signed comparison: 0 > root1? if so, set xmm4 dword to 1's */ \
            "vpcmpgtd	%%ymm2, %%ymm5, %%ymm5 \n\t"		/* signed comparison: 0 > root2? if so, set xmm5 dword to 1's */ \
            "vpand	    %%ymm0, %%ymm4, %%ymm4 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
            "vmovdqa    %%ymm0, %%ymm6 \n\t"			/* copy prime to compute neg roots */ \
            "vpand	    %%ymm0, %%ymm5, %%ymm5 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
            "vmovdqa    %%ymm0, %%ymm7 \n\t"			/* copy prime to compute neg roots */ \
            "vpaddd	    %%ymm4, %%ymm1, %%ymm1 \n\t"			/* selectively add back prime (modular subtract) */ \
            "vpaddd	    %%ymm5, %%ymm2, %%ymm2 \n\t"			/* selectively add back prime (modular subtract) */ \
            "vmovdqa    %%ymm1, %%ymm8 \n\t"					/* xmm5 = root1 copy */	\
            "vpminud	%%ymm2, %%ymm1, %%ymm1 \n\t"					/* xmm2 = root2 < root1 ? root2 : root1 */	\
            "vpmaxud	%%ymm8, %%ymm2, %%ymm2 \n\t"					/* xmm5 = root2 > root1 ? root2 : root1 */	\
            "vmovdqa    %%ymm1, (%%r14,%%r15,4) \n\t"	/* save new root1 values */ \
            "vmovdqa    %%ymm2, (%%r13,%%r15,4) \n\t"	/* save new root2 values */ \
            "vpsubd	    %%ymm1, %%ymm6, %%ymm6 \n\t"			/* form negative root1's; prime - root1 */ \
            "vpsubd	    %%ymm2, %%ymm7, %%ymm7 \n\t"			/* form negative root2's; prime - root2 */ \
            "vpand	    %%ymm1, %%ymm12, %%ymm10 \n\t"			/* mask the root1s */ \
            "vpand	    %%ymm2, %%ymm12, %%ymm11 \n\t"			/* mask the root2s */ \
            "vpsrld	    $15, %%ymm1, %%ymm1 \n\t"			 /* shift roots to compute blocknums */ \
            "vpsrld	    $15, %%ymm2, %%ymm2 \n\t"			 /* shift roots to compute blocknums */ \
            "vpor	    %%ymm10, %%ymm15, %%ymm10 \n\t"			/* combine with the prime index half of the 32-bit word */ \
            "vpor	    %%ymm11, %%ymm15, %%ymm11 \n\t"			/* combine with the prime index half of the 32-bit word */ \
            "movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
            "movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
            "movl	92(%%rsi,1),%%r13d \n\t"		/* 88 = interval, 92 = numblocks */ \
            /* do a vector compare and extract a bitmask, then the comparisons don't depend */ \
            /* on vpextrd and we can skip the vpextrd unless the test succeeds... */ \
            /* both eax and edx are available for use... */ \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,1 ========== */	\
            /* ================================================ */	\
            "vpextrd	$0,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$0,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "cmpl	%%r8d,%%r13d \n\t"				/* root1 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("0") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,1 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("0") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,2 ========== */	\
            /* ================================================ */	\
            "vpextrd	$1,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$1,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "cmpl	%%r8d,%%r13d \n\t"				/* root1 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("1") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,2 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("1") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,3 ========== */	\
            /* ================================================ */	\
            "vpextrd	$2,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$2,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "cmpl	%%r8d,%%r13d \n\t"				/* root1 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("2") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,3 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("2") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,4 ========== */	\
            /* ================================================ */	\
            "vpextrd	$3,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$3,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "cmpl	%%r8d,%%r13d \n\t"				/* root1 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("3") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,4 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("3") \
            "2:		\n\t" \
            "1:		\n\t" \
            "vextracti128   $1,%%ymm1,%%xmm1 \n\t" \
            "vextracti128   $1,%%ymm2,%%xmm2 \n\t" \
            "vextracti128   $1,%%ymm10,%%xmm10 \n\t" \
            "vextracti128   $1,%%ymm11,%%xmm11 \n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,5 ========== */	\
            /* ================================================ */	\
            "vpextrd	$0,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$0,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "cmpl	%%r8d,%%r13d \n\t"				/* root1 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("4") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,5 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("4") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,6 ========== */	\
            /* ================================================ */	\
            "vpextrd	$1,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$1,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "cmpl	%%r8d,%%r13d \n\t"				/* root1 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("5") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,6 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("5") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,7 ========== */	\
            /* ================================================ */	\
            "vpextrd	$2,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$2,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "cmpl	%%r8d,%%r13d \n\t"				/* root1 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("6") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,7 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("6") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,8 ========== */	\
            /* ================================================ */	\
            "vpextrd	$3,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$3,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "cmpl	%%r8d,%%r13d \n\t"				/* root1 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("7") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,8 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("7") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,1 ========== */	\
            /* ================================================ */	\
            "movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
            "movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
            "vpand	    %%ymm6, %%ymm12, %%ymm10 \n\t"			/* mask the root1s */ \
            "vpand	    %%ymm7, %%ymm12, %%ymm11 \n\t"			/* mask the root2s */ \
            "vpsrld	    $15, %%ymm6, %%ymm6 \n\t"			 /* shift roots to compute blocknums */ \
            "vpsrld	    $15, %%ymm7, %%ymm7 \n\t"			 /* shift roots to compute blocknums */ \
            "vpor	    %%ymm10, %%ymm15, %%ymm10 \n\t"			/* combine with the prime index half of the 32-bit word */ \
            "vpor	    %%ymm11, %%ymm15, %%ymm11 \n\t"			/* combine with the prime index half of the 32-bit word */ \
            "vpextrd	$0,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$0,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("0") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,1 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r8d,%%r13d \n\t"				/* root1 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("0") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,2 ========== */	\
            /* ================================================ */	\
            "vpextrd	$1,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$1,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("1") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,2 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r8d,%%r13d \n\t"				/* root1 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("1") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,3 ========== */	\
            /* ================================================ */	\
            "vpextrd	$2,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$2,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("2") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,3 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r8d,%%r13d \n\t"				/* root1 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("2") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,4 ========== */	\
            /* ================================================ */	\
            "vpextrd	$3,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$3,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("3") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,4 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r8d,%%r13d \n\t"				/* root1 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("3") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,5 ========== */	\
            /* ================================================ */	\
            "vextracti128   $1,%%ymm7,%%xmm7 \n\t" \
            "vextracti128   $1,%%ymm6,%%xmm6 \n\t" \
            "vextracti128   $1,%%ymm10,%%xmm10 \n\t" \
            "vextracti128   $1,%%ymm11,%%xmm11 \n\t" \
            "vpextrd	$0,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$0,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("4") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,5 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r8d,%%r13d \n\t"				/* root1 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("4") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,6 ========== */	\
            /* ================================================ */	\
            "vpextrd	$1,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$1,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("5") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,6 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r8d,%%r13d \n\t"				/* root1 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("5") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,7 ========== */	\
            /* ================================================ */	\
            "vpextrd	$2,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$2,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("6") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,7 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r8d,%%r13d \n\t"				/* root1 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("6") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,8 ========== */	\
            /* ================================================ */	\
            "vpextrd	$3,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$3,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("7") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,8 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r8d,%%r13d \n\t"				/* root1 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("7") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* ======== END OF LOOP - UPDATE AND CHECK ======== */	\
            /* ================================================ */	\
            "7: \n\t" \
            "addl   $8,%%r15d \n\t"					/* increment j by 1*/ \
            "vpaddd  %%ymm13, %%ymm15, %%ymm15 \n\t" /* add in shifted-j: (j+8)<<16 */ \
            "cmpl   84(%%rsi,1),%%r15d \n\t"		/* j < bound ? */ \
            "jb     8b \n\t"	\
            "9:		\n\t"				\
            "movl	%%r15d, %%eax \n\t" \
            :  \
            : "g"(&helperstruct) \
            : "rax", "rbx", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11", "r13", "r14", "r15",
            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm10", "xmm11", "xmm13", "xmm14", "xmm15", "memory", "cc");


            //printf("all positive roots skipped %d times out of %d\n",
            //    helperstruct.filler, (bound - large_B) / 8);
		bound_index = helperstruct.bound_index;
		logp = helperstruct.logp;

	}
	else
	{
        /////////////////////////////////////////////////////////////////////////////////
        // sign < 0
        /////////////////////////////////////////////////////////////////////////////////

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
		// we can start doing things 16 at a time and be
		// sure we can use aligned moves (static_data_init).	
		for (j=sconf->sieve_small_fb_start; 
			j < sconf->factor_base->fb_10bit_B; j++, ptr++)
		{

			// as soon as we are aligned, use more efficient avx2 based methods...
			if ((j & 15) == 0)
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
		
		// update 16 at a time using avx2 and no branching		
		sm_ptr = &dconf->sm_rootupdates[(v-1) * med_B];
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
			h.stop = sconf->factor_base->med_B;		// 68

            COMPUTE_16X_SMALL_NROOTS_AVX2;
			
			j = h.stop;
		}
		sm_ptr = &dconf->sm_rootupdates[(v-1) * med_B + j];

		bound_index = 0;
		bound_val = med_B;
		check_bound = med_B + BUCKET_ALLOC/2;
		logp = update_data.logp[med_B-1];

		//if (VFLAG > 2)
		//{
		//	printf("commencing 1n\n");
		//}

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
        helperstruct.filler = 42;
        helperstruct.scratch = dconf->polyscratch;        //120
		
        ASM_G (\
            "movq	%0,%%rsi \n\t"					/* move helperstruct into rsi */ \
            "movl   80(%%rsi,1),%%r15d \n\t"		/* large_B = j = r15d */ \
            "vmovd      %%r15d, %%xmm15 \n\t"           /* broadcast j to ymm15 */ \
            "vpshufd	    $0, %%xmm15, %%xmm15 \n\t"  /* broadcast j to ymm15 */ \
            "vinserti128    $1, %%xmm15, %%ymm15, %%ymm15 \n\t" \
            "movq	120(%%rsi,1),%%rdi \n\t"			/* edi = polyscratch */ \
            "vpaddd     (%%rdi), %%ymm15, %%ymm15 \n\t" /* add in lane offsets */ \
            "movl	    96(%%rsi,1),%%r13d \n\t"		/* store bound_val in a register */ \
            "vmovd      %%r13d, %%xmm14 \n\t"           /* broadcast bound_val to ymm14 */ \
            "vpshufd	    $0, %%xmm14, %%xmm14 \n\t"  /* broadcast bound_val to ymm14 */ \
            "vinserti128    $1, %%xmm14, %%ymm14, %%ymm14 \n\t" \
            "vpsubd         %%ymm14, %%ymm15, %%ymm15 \n\t" /* subtract bound_val */ \
            "vpslld	$16, %%ymm15, %%ymm15 \n\t"		/* put (j - bound_val + i) into high half of 32-bit words */ \
            "movl      $524288, %%edx \n\t"           /* make the shifted-j increment */ \
            "vmovd      %%edx, %%xmm13 \n\t"           /* make the shifted-j increment */ \
            "vpshufd	    $0, %%xmm13, %%xmm13 \n\t"  /* broadcast shifted-j to ymm13 */ \
            "vinserti128    $1, %%xmm13, %%ymm13, %%ymm13 \n\t" \
            /* do the loop comparison */ \
            "cmpl   84(%%rsi,1),%%r15d \n\t"		/* j >= bound ? */ \
            "jae    9f	\n\t"						/* jump to end of loop, if test fails */ \
            "8: \n\t"	\
            /* ================================================ */	\
            /* ========== BEGIN CHECK_NEW_SLICE BLOCK ========= */	\
            /* ================================================ */	\
            CHECK_NEW_SLICE_ASM_NEW	\
            /* ================================================ */	\
            /* ============ BEGIN - GET NEW ROOTS 1 =========== */	\
            /* ================================================ */	\
            "movq	72(%%rsi,1),%%rdi \n\t"			/* edi = ptr */ \
            "movq	40(%%rsi,1),%%r14 \n\t"			/* move updata_data.root1 pointer into r14 */ \
            "vmovdqa (%%rdi,%%r15,4), %%ymm3 \n\t"	/* xmm3 = next 4 values of rootupdates */ \
            "movq	48(%%rsi,1),%%r13 \n\t"			/* move updata_data.root2 pointer into r13 */ \
            "vmovdqa (%%r14,%%r15,4), %%ymm1 \n\t"	/* xmm1 = next 4 values of root1 */ \
            "vpaddd	%%ymm3, %%ymm1, %%ymm1 \n\t"			/* root1 += ptr */ \
            "vmovdqa (%%r13,%%r15,4), %%ymm2 \n\t"	/* xmm2 = next 4 values of root2 */ \
            "movq	32(%%rsi,1),%%rdi \n\t"			/* move updata_data.prime pointer into r12 */ \
            "vpaddd	    %%ymm3, %%ymm2, %%ymm2 \n\t"			/* root2 += ptr */ \
            "vmovdqa	%%ymm1, %%ymm4 \n\t"			/* copy root1 to xmm4 */ \
            "vmovdqa    (%%rdi,%%r15,4), %%ymm0 \n\t"	/* xmm0 = next 4 primes */ \
            "vmovdqa	%%ymm2, %%ymm5 \n\t"			/* copy root2 to xmm5 */ \
            "vpcmpgtd	%%ymm0, %%ymm4, %%ymm4 \n\t"		/* signed comparison: root1 > p? if so, set xmm4 dword to 1's */ \
            "vpcmpgtd	%%ymm0, %%ymm5, %%ymm5 \n\t"		/* signed comparison: root2 > p? if so, set xmm5 dword to 1's */ \
            "vmovdqa    %%ymm0, %%ymm6 \n\t"			/* copy of prime for neg root calculation */ \
            "vmovdqa    %%ymm0, %%ymm7 \n\t"			/* copy of prime for neg root calculation */ \
            "vmovdqa    %%ymm0, %%ymm8 \n\t"			/* copy of prime for neg root loops */ \
            "vpand	    %%ymm0, %%ymm4, %%ymm4 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
            "vpand	    %%ymm0, %%ymm5, %%ymm5 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
            "vpsubd	    %%ymm4, %%ymm1, %%ymm1 \n\t"			/* selectively sub back prime (modular addition) */ \
            "vpsubd	    %%ymm5, %%ymm2, %%ymm2 \n\t"			/* selectively sub back prime (modular addition) */ \
            "vmovdqa    %%ymm1, %%ymm9 \n\t"					/* xmm5 = root1 copy */	\
            "vpminud	%%ymm2, %%ymm1, %%ymm1 \n\t"					/* xmm2 = root2 < root1 ? root2 : root1 */	\
            "vpmaxud	%%ymm9, %%ymm2, %%ymm2 \n\t"					/* xmm5 = root2 > root1 ? root2 : root1 */	\
            "vmovdqa    %%ymm1, (%%r14,%%r15,4) \n\t"	/* save new root1 values */ \
            "vmovdqa    %%ymm2, (%%r13,%%r15,4) \n\t"	/* save new root2 values */ \
            "vpsubd	    %%ymm1, %%ymm6, %%ymm6 \n\t"			/* form negative root1's; prime - root1 */ \
            "vpsubd	    %%ymm2, %%ymm7, %%ymm7 \n\t"			/* form negative root2's; prime - root2 */ \
            "vextracti128   $0,%%ymm15,%%xmm10 \n\t" \
            "movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
            "movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
            "movl	88(%%rsi,1),%%r13d \n\t"		/* interval */ \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,1 ========== */	\
            /* ================================================ */	\
            "vpextrd	$0,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$0,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$0,%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT1_LOOP_NEW("0") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,1 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT2_LOOP_NEW("0") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,2 ========== */	\
            /* ================================================ */	\
            "vpextrd	$1,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$1,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$1,%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT1_LOOP_NEW("1") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,2 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT2_LOOP_NEW("1") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,3 ========== */	\
            /* ================================================ */	\
            "vpextrd	$2,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$2,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$2,%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT1_LOOP_NEW("2") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,3 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT2_LOOP_NEW("2") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,4 ========== */	\
            /* ================================================ */	\
            "vpextrd	$3,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$3,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$3,%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT1_LOOP_NEW("3") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,4 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT2_LOOP_NEW("3") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,5 ========== */	\
            /* ================================================ */	\
            "vextracti128   $1,%%ymm1,%%xmm1 \n\t" \
            "vextracti128   $1,%%ymm2,%%xmm2 \n\t" \
            "vextracti128   $1,%%ymm0,%%xmm0 \n\t" \
            "vextracti128   $1,%%ymm15,%%xmm10 \n\t" \
            "vpextrd	$0,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$0,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$0,%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT1_LOOP_NEW("4") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,5 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT2_LOOP_NEW("4") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,6 ========== */	\
            /* ================================================ */	\
            "vpextrd	$1,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$1,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$1,%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT1_LOOP_NEW("5") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,6 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT2_LOOP_NEW("5") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,7 ========== */	\
            /* ================================================ */	\
            "vpextrd	$2,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$2,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$2,%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT1_LOOP_NEW("6") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,7 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT2_LOOP_NEW("6") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,8 ========== */	\
            /* ================================================ */	\
            "vpextrd	$3,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$3,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$3,%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT1_LOOP_NEW("7") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,8 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT2_LOOP_NEW("7") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,1 ========== */	\
            /* ================================================ */	\
            "movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
            "movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
            "vextracti128   $0,%%ymm15,%%xmm10 \n\t" \
            "vpextrd	$0,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$0,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$0,%%xmm8,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT2_LOOP_NEW("0") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,1 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT1_LOOP_NEW("0") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,2 ========== */	\
            /* ================================================ */	\
            "vpextrd	$1,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$1,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$1,%%xmm8,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT2_LOOP_NEW("1") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,2 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT1_LOOP_NEW("1") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,3 ========== */	\
            /* ================================================ */	\
            "vpextrd	$2,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$2,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$2,%%xmm8,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT2_LOOP_NEW("2") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,3 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT1_LOOP_NEW("2") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,4 ========== */	\
            /* ================================================ */	\
            "vpextrd	$3,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$3,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$3,%%xmm8,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT2_LOOP_NEW("3") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,4 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT1_LOOP_NEW("3") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,5 ========== */	\
            /* ================================================ */	\
            "vextracti128   $1,%%ymm7,%%xmm7 \n\t" \
            "vextracti128   $1,%%ymm6,%%xmm6 \n\t" \
            "vextracti128   $1,%%ymm8,%%xmm8 \n\t" \
            "vextracti128   $1,%%ymm15,%%xmm10 \n\t" \
            "vpextrd	$0,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$0,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$0,%%xmm8,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT2_LOOP_NEW("4") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,5 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT1_LOOP_NEW("4") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,6 ========== */	\
            /* ================================================ */	\
            "vpextrd	$1,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$1,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$1,%%xmm8,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT2_LOOP_NEW("5") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,6 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT1_LOOP_NEW("5") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,7 ========== */	\
            /* ================================================ */	\
            "vpextrd	$2,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$2,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$2,%%xmm8,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT2_LOOP_NEW("6") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,7 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT1_LOOP_NEW("6") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,8 ========== */	\
            /* ================================================ */	\
            "vpextrd	$3,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$3,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$3,%%xmm8,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            UPDATE_ROOT2_LOOP_NEW("7") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,8 ========== */	\
            /* ================================================ */	\
            UPDATE_ROOT1_LOOP_NEW("7") \
            /* ================================================ */	\
            /* ======== END OF LOOP - UPDATE AND CHECK ======== */	\
            /* ================================================ */	\
            "addl   $8,%%r15d \n\t"					/* increment j by 1*/ \
            "vpaddd  %%ymm13, %%ymm15, %%ymm15 \n\t" /* add in shifted-j: (j+8)<<16 */ \
            "cmpl   84(%%rsi,1),%%r15d \n\t"		/* j < bound ? */ \
            "jb     8b \n\t"	\
            "9:		\n\t"				\
            "movl	%%r15d, %%eax \n\t" \
            :  \
            : "g"(&helperstruct) \
            : "rax", "rbx", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11", "r13", "r14", "r15", 
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
		logp = update_data.logp[large_B-1];

		//if (VFLAG > 2)
		//{
		//	printf("commencing 2n\n");
		//}

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
        helperstruct.filler = 42;
        helperstruct.scratch = dconf->polyscratch;        //120

        ASM_G ( \
            "movq	%0,%%rsi \n\t"					/* move helperstruct into rsi */ \
            "movl   80(%%rsi,1),%%r15d \n\t"		/* large_B = j = r15d */ \
            "vmovd      %%r15d, %%xmm15 \n\t"           /* broadcast j to ymm15 */ \
            "vpshufd	    $0, %%xmm15, %%xmm15 \n\t"  /* broadcast j to ymm15 */ \
            "vinserti128    $1, %%xmm15, %%ymm15, %%ymm15 \n\t" \
            "movq	120(%%rsi,1),%%rdi \n\t"			/* edi = polyscratch */ \
            "vpaddd     (%%rdi), %%ymm15, %%ymm15 \n\t" /* add in lane offsets */ \
            "movl	    96(%%rsi,1),%%r13d \n\t"		/* store bound_val in a register */ \
            "vmovd      %%r13d, %%xmm14 \n\t"           /* broadcast bound_val to ymm14 */ \
            "vpshufd	    $0, %%xmm14, %%xmm14 \n\t"  /* broadcast bound_val to ymm14 */ \
            "vinserti128    $1, %%xmm14, %%ymm14, %%ymm14 \n\t" \
            "vpsubd         %%ymm14, %%ymm15, %%ymm15 \n\t" /* subtract bound_val */ \
            "vpslld	$16, %%ymm15, %%ymm15 \n\t"		/* put (j - bound_val + i) into high half of 32-bit words */ \
            "movl      $524288, %%edx \n\t"           /* make the shifted-j increment */ \
            "vmovd      %%edx, %%xmm13 \n\t"           /* make the shifted-j increment */ \
            "vpshufd	    $0, %%xmm13, %%xmm13 \n\t"  /* broadcast shifted-j to ymm13 */ \
            "vinserti128    $1, %%xmm13, %%ymm13, %%ymm13 \n\t" \
            "vmovdqa     32(%%rdi), %%ymm12 \n\t" /* put the root masks into ymm12 */ \
            /* do the loop comparison */ \
            "cmpl   84(%%rsi,1),%%r15d \n\t"		/* j >= bound ? */ \
            "jae    9f	\n\t"						/* jump to end of loop, if test fails */ \
            "8: \n\t"	\
            /* ================================================ */	\
            /* ========== BEGIN CHECK_NEW_SLICE BLOCK ========= */	\
            /* ================================================ */	\
            CHECK_NEW_SLICE_ASM_NEW	\
            /* ================================================ */	\
            /* ============ BEGIN - GET NEW ROOTS 1 =========== */	\
            /* ================================================ */	\
            "movq	72(%%rsi,1),%%rdi \n\t"			/* edi = ptr */ \
            "movq	40(%%rsi,1),%%r14 \n\t"			/* move updata_data.root1 pointer into r14 */ \
            "vmovdqa (%%rdi,%%r15,4), %%ymm3 \n\t"	/* xmm3 = next 4 values of rootupdates */ \
            "movq	48(%%rsi,1),%%r13 \n\t"			/* move updata_data.root2 pointer into r13 */ \
            "vmovdqa (%%r14,%%r15,4), %%ymm1 \n\t"	/* xmm1 = next 4 values of root1 */ \
            "vpaddd	%%ymm3, %%ymm1, %%ymm1 \n\t"			/* root1 += ptr */ \
            "vmovdqa (%%r13,%%r15,4), %%ymm2 \n\t"	/* xmm2 = next 4 values of root2 */ \
            "movq	32(%%rsi,1),%%rdi \n\t"			/* move updata_data.prime pointer into r12 */ \
            "vpaddd	    %%ymm3, %%ymm2, %%ymm2 \n\t"			/* root2 += ptr */ \
            "vmovdqa	%%ymm1, %%ymm4 \n\t"			/* copy root1 to xmm4 */ \
            "vmovdqa    (%%rdi,%%r15,4), %%ymm0 \n\t"	/* xmm0 = next 4 primes */ \
            "vmovdqa	%%ymm2, %%ymm5 \n\t"			/* copy root2 to xmm5 */ \
            "vpcmpgtd	%%ymm0, %%ymm4, %%ymm4 \n\t"		/* signed comparison: root1 > p? if so, set xmm4 dword to 1's */ \
            "vpcmpgtd	%%ymm0, %%ymm5, %%ymm5 \n\t"		/* signed comparison: root2 > p? if so, set xmm5 dword to 1's */ \
            "vmovdqa    %%ymm0, %%ymm6 \n\t"			/* copy prime to xmm6 */ \
            "vpand	    %%ymm0, %%ymm4, %%ymm4 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
            "vmovdqa    %%ymm0, %%ymm7 \n\t"			/* copy prime to xmm7 */ \
            "vpand	    %%ymm0, %%ymm5, %%ymm5 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
            "vpsubd	    %%ymm4, %%ymm1, %%ymm1 \n\t"			/* selectively sub back prime (modular addition) */ \
            "vpsubd	    %%ymm5, %%ymm2, %%ymm2 \n\t"			/* selectively sub back prime (modular addition) */ \
            "vmovdqa    %%ymm1, %%ymm8 \n\t"					/* xmm5 = root1 copy */	\
            "vpminud	%%ymm2, %%ymm1, %%ymm1 \n\t"					/* xmm2 = root2 < root1 ? root2 : root1 */	\
            "vpmaxud	%%ymm8, %%ymm2, %%ymm2 \n\t"					/* xmm5 = root2 > root1 ? root2 : root1 */	\
            "vmovdqa    %%ymm1, (%%r14,%%r15,4) \n\t"	/* save new root1 values */ \
            "vmovdqa    %%ymm2, (%%r13,%%r15,4) \n\t"	/* save new root2 values */ \
            "vpsubd	    %%ymm1, %%ymm6, %%ymm6 \n\t"			/* form negative root1's; prime - root1 */ \
            "vpsubd	    %%ymm2, %%ymm7, %%ymm7 \n\t"			/* form negative root2's; prime - root2 */ \
            "vpand	    %%ymm1, %%ymm12, %%ymm10 \n\t"			/* mask the root1s */ \
            "vpand	    %%ymm2, %%ymm12, %%ymm11 \n\t"			/* mask the root2s */ \
            /* after this point the roots are only needed as block numbers */ \
            "vpsrld	    $15, %%ymm1, %%ymm1 \n\t"			 /* shift roots to compute blocknums */ \
            "vpsrld	    $15, %%ymm2, %%ymm2 \n\t"			 /* shift roots to compute blocknums */ \
            "vpor	    %%ymm10, %%ymm15, %%ymm10 \n\t"			/* combine with the prime index half of the 32-bit word */ \
            "vpor	    %%ymm11, %%ymm15, %%ymm11 \n\t"			/* combine with the prime index half of the 32-bit word */ \
            "movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
            "movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
            "movl	92(%%rsi,1),%%r13d \n\t"		/* 88 = interval, 92 = numblocks */ \
            /* edx is free at this point... can it be used? */ \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,1 ========== */	\
            /* ================================================ */	\
            "vpextrd	$0,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$0,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "cmpl	%%r8d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("0") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,1 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("0") \
            "2: \n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,2 ========== */	\
            /* ================================================ */	\
            "vpextrd	$1,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$1,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "cmpl	%%r8d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("1") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,2 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("1") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,3 ========== */	\
            /* ================================================ */	\
            "vpextrd	$2,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$2,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "cmpl	%%r8d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("2") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,3 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("2") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,4 ========== */	\
            /* ================================================ */	\
            "vpextrd	$3,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$3,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "cmpl	%%r8d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("3") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,4 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("3") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,5 ========== */	\
            /* ================================================ */	\
            "vextracti128   $1,%%ymm1,%%xmm1 \n\t" \
            "vextracti128   $1,%%ymm2,%%xmm2 \n\t" \
            "vextracti128   $1,%%ymm10,%%xmm10 \n\t" \
            "vextracti128   $1,%%ymm11,%%xmm11 \n\t" \
            "vpextrd	$0,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$0,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "cmpl	%%r8d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("4") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,5 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("4") \
            "2: \n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,6 ========== */	\
            /* ================================================ */	\
            "vpextrd	$1,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$1,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "cmpl	%%r8d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("5") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,6 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("5") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,7 ========== */	\
            /* ================================================ */	\
            "vpextrd	$2,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$2,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "cmpl	%%r8d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("6") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,7 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("6") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,8 ========== */	\
            /* ================================================ */	\
            "vpextrd	$3,%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "vpextrd	$3,%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "cmpl	%%r8d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("7") \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,8 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("7") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,1 ========== */	\
            /* ================================================ */	\
            "movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
            "movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
            "vpand	    %%ymm6, %%ymm12, %%ymm10 \n\t"			/* mask the root1s */ \
            "vpand	    %%ymm7, %%ymm12, %%ymm11 \n\t"			/* mask the root2s */ \
            "vpsrld	    $15, %%ymm6, %%ymm6 \n\t"			 /* shift roots to compute blocknums */ \
            "vpsrld	    $15, %%ymm7, %%ymm7 \n\t"			 /* shift roots to compute blocknums */ \
            "vpor	    %%ymm10, %%ymm15, %%ymm10 \n\t"			/* combine with the prime index half of the 32-bit word */ \
            "vpor	    %%ymm11, %%ymm15, %%ymm11 \n\t"			/* combine with the prime index half of the 32-bit word */ \
            "vpextrd	$0,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$0,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("0") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,1 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r8d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("0") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,2 ========== */	\
            /* ================================================ */	\
            "vpextrd	$1,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$1,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("1") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,2 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r8d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("1") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,3 ========== */	\
            /* ================================================ */	\
            "vpextrd	$2,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$2,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("2") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,3 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r8d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("2") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,4 ========== */	\
            /* ================================================ */	\
            "vpextrd	$3,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$3,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("3") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,4 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r8d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("3") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,5 ========== */	\
            /* ================================================ */	\
            "vextracti128   $1,%%ymm7,%%xmm7 \n\t" \
            "vextracti128   $1,%%ymm6,%%xmm6 \n\t" \
            "vextracti128   $1,%%ymm10,%%xmm10 \n\t" \
            "vextracti128   $1,%%ymm11,%%xmm11 \n\t" \
            "vpextrd	$0,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$0,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("4") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,5 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r8d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("4") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,6 ========== */	\
            /* ================================================ */	\
            "vpextrd	$1,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$1,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("5") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,6 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r8d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("5") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,7 ========== */	\
            /* ================================================ */	\
            "vpextrd	$2,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$2,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("6") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,7 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r8d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("6") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,8 ========== */	\
            /* ================================================ */	\
            "vpextrd	$3,%%xmm7,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "vpextrd	$3,%%xmm6,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "cmpl	%%r9d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    2f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_NEW("7") \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,8 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r8d,%%r13d \n\t"				/* root2 > interval? */ \
            "jbe    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_NEW("7") \
            "2:		\n\t" \
            "1:		\n\t" \
            /* ================================================ */	\
            /* ======== END OF LOOP - UPDATE AND CHECK ======== */	\
            /* ================================================ */	\
            "addl   $8,%%r15d \n\t"					/* increment j by 1*/ \
            "vpaddd  %%ymm13, %%ymm15, %%ymm15 \n\t" /* add in shifted-j: (j+8)<<16 */ \
            "cmpl   84(%%rsi,1),%%r15d \n\t"		/* j < bound ? */ \
            "jb     8b \n\t"	\
            "9:		\n\t"				\
            "movl	%%r15d, %%eax \n\t" \
            :  \
            : "g"(&helperstruct) \
            : "rax", "rbx", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11", "r13", "r14", "r15",
            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm10", "xmm11", "xmm13", "xmm14", "xmm15", "memory", "cc");



		bound_index = helperstruct.bound_index;
		logp = helperstruct.logp;

	}

    CLEAN_AVX2;

	if (lp_bucket_p->list != NULL)
	{
		lp_bucket_p->num_slices = bound_index + 1;
		lp_bucket_p->logp[bound_index] = logp;
	}

	return;
}

#endif

#if defined( USE_AVX2 )

void nextRoots_32k_avx2_small_test(static_conf_t* sconf, dynamic_conf_t* dconf)
{
    //update the roots 
    sieve_fb_compressed* fb_p = dconf->comp_sieve_p;
    sieve_fb_compressed* fb_n = dconf->comp_sieve_n;
    int* rootupdates = dconf->rootupdates;

    update_t update_data = dconf->update_data;

    uint32 startprime = 2;
    uint32 bound = sconf->factor_base->B;

    char v = dconf->curr_poly->nu[dconf->numB];
    char sign = dconf->curr_poly->gray[dconf->numB];
    int* ptr;
    uint16* sm_ptr;

    lp_bucket* lp_bucket_p = dconf->buckets;
    uint32 med_B = sconf->factor_base->med_B;
    uint32 large_B = sconf->factor_base->large_B;

    uint32 j, interval; //, fb_offset;
    int k, numblocks;
    uint32 root1, root2, nroot1, nroot2, prime;

    int bound_index = 0;
    int check_bound = BUCKET_ALLOC / 2 - 1;
    uint32 bound_val = med_B;
    uint32* numptr_p, * numptr_n, * sliceptr_p, * sliceptr_n;

#if !defined(USE_POLY_SSE2_ASM) || defined(PROFILING)
    uint32 * bptr;
    int bnum, room;
#endif

    uint8 logp = 0;
    polysieve_t helperstruct;

    numblocks = sconf->num_blocks;
    interval = numblocks << 15;

    if (lp_bucket_p->alloc_slices != 0) // != NULL)
    {
        lp_bucket_p->fb_bounds[0] = med_B;

        sliceptr_p = lp_bucket_p->list;
        sliceptr_n = lp_bucket_p->list + (numblocks << BUCKET_BITS);

        numptr_p = lp_bucket_p->num;
        numptr_n = lp_bucket_p->num + numblocks;

        // reset bucket counts
        for (j = 0; j < lp_bucket_p->list_size; j++)
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

    k = 0;
    ptr = &rootupdates[(v - 1) * bound + startprime];

    if (sign > 0)
    {
#ifdef QS_TIMING
        gettimeofday(&qs_timing_start, NULL);
#endif

        for (j = startprime; j < sconf->sieve_small_fb_start; j++, ptr++)
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
        for (j = sconf->sieve_small_fb_start;
            j < sconf->factor_base->fb_10bit_B; j++, ptr++)
        {
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

#if defined(HAS_SSE2) && (defined(GCC_ASM64X) || defined(_MSC_VER)) //NOTDEF //GCC_ASM64X

        // update 8 at a time using SSE2 and no branching
        sm_ptr = &dconf->sm_rootupdates[(v - 1) * med_B];
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
            h.start = sconf->factor_base->fb_10bit_B;		// 64
            h.stop = sconf->factor_base->fb_15bit_B;		// 68
            if ((h.stop - 8) > h.start)
                h.stop -= 8;

            COMPUTE_8X_SMALL_PROOTS;

            j = h.stop;
        }

#else


        ptr = &dconf->rootupdates[(v - 1) * bound + sconf->factor_base->fb_10bit_B];
        for (j = sconf->factor_base->fb_10bit_B; j < sconf->factor_base->fb_15bit_B; j++, ptr++)
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
        ptr = &dconf->rootupdates[(v - 1) * bound + j];
        for (; j < med_B; j++, ptr++)
        {
            prime = update_data.prime[j];
            root1 = (uint16)update_data.sm_firstroots1[j];
            root2 = (uint16)update_data.sm_firstroots2[j];

            if ((prime > 32768) && ((j & 7) == 0))
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

        // continue one at a time once we exceed 15 bits, because the 8x SSE2
        // code has a hard time with unsigned 16 bit comparisons
        for (; j < med_B; j++, ptr++)
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

#ifdef QS_TIMING
        gettimeofday(&qs_timing_stop, NULL);
        POLY_STG2 += yafu_difftime(&qs_timing_start, &qs_timing_stop);
        gettimeofday(&qs_timing_start, NULL);
#endif

        bound_index = 0;
        bound_val = med_B;
        check_bound = med_B + BUCKET_ALLOC / 2;



#if defined(USE_POLY_SSE2_ASM) && defined(GCC_ASM64X) && !defined(PROFILING)
        logp = update_data.logp[med_B - 1];

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
        helperstruct.ptr = &rootupdates[(v - 1) * bound];						//72
        helperstruct.large_B = med_B;				//80
        helperstruct.B = large_B;					//84
        helperstruct.interval = interval;			//88
        helperstruct.numblocks = numblocks;			//92
        helperstruct.bound_val = bound_val;			//96
        helperstruct.bound_index = bound_index;		//100
        helperstruct.check_bound = check_bound;		//104
        helperstruct.logp = logp;					//108

        ASM_G(\
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
            "movdqa %%xmm1, (%%r14,%%r15,4) \n\t"	/* save new root1 values */ \
            "paddd	%%xmm5, %%xmm2 \n\t"			/* selectively add back prime (modular subtract) */ \
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
            "movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "movd	%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
            "jae    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_LOOP("0") \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,1 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "psrldq $4,%%xmm2 \n\t" 				/* next root1 */ \
            "jae    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_LOOP("0") \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,2 ========== */	\
            /* ================================================ */	\
            "psrldq $4,%%xmm0 \n\t" 				/* next prime */ \
            "movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,2 from xmm1 */ \
            "movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,2 from xmm2 */ \
            "movd	%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
            "jae    3f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_LOOP("1") \
            "3:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,2 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "psrldq $4,%%xmm2 \n\t" 				/* next root1 */ \
            "jae    3f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_LOOP("1") \
            "3:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,3 ========== */	\
            /* ================================================ */	\
            "psrldq $4,%%xmm0 \n\t" 				/* next prime */ \
            "movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,3 from xmm1 */ \
            "movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,3 from xmm2 */ \
            "movd	%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
            "jae    5f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_LOOP("2") \
            "5:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,3 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "psrldq $4,%%xmm2 \n\t" 				/* next root1 */ \
            "jae    5f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_LOOP("2") \
            "5:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,4 ========== */	\
            /* ================================================ */	\
            "psrldq $4,%%xmm0 \n\t" 				/* next prime */ \
            "movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,4 from xmm1 */ \
            "movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,4 from xmm2 */ \
            "movd	%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "jae    7f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_LOOP("3") \
            "7:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,4 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "jae    7f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_LOOP("3") \
            "7:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,1 ========== */	\
            /* ================================================ */	\
            "movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
            "movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
            "movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,1 from xmm6 */ \
            "movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,1 from xmm7 */ \
            "movd	%%xmm8,%%edx \n\t"				/* else, extract prime from xmm8 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "psrldq $4,%%xmm6 \n\t" 				/* nextn root1 */ \
            "jae    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_LOOP("0") \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,1 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "psrldq $4,%%xmm7 \n\t" 				/* nextn nroot2 */ \
            "jae    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_LOOP("0") \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,2 ========== */	\
            /* ================================================ */	\
            "psrldq $4,%%xmm8 \n\t" 				/* next prime */ \
            "movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,2 from xmm6 */ \
            "movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,2 from xmm7 */ \
            "movd	%%xmm8,%%edx \n\t"				/* else, extract prime from xmm8 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "psrldq $4,%%xmm6 \n\t" 				/* nextn root1 */ \
            "jae    3f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_LOOP("1") \
            "3:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,2 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "psrldq $4,%%xmm7 \n\t" 				/* nextn nroot2 */ \
            "jae    3f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_LOOP("1") \
            "3:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,3 ========== */	\
            /* ================================================ */	\
            "psrldq $4,%%xmm8 \n\t" 				/* next prime */ \
            "movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,3 from xmm6 */ \
            "movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,3 from xmm7 */ \
            "movd	%%xmm8,%%edx \n\t"				/* else, extract prime from xmm8 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "psrldq $4,%%xmm6 \n\t" 				/* nextn root1 */ \
            "jae    5f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_LOOP("2") \
            "5:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,3 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "psrldq $4,%%xmm7 \n\t" 				/* nextn nroot2 */ \
            "jae    5f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_LOOP("2") \
            "5:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,4 ========== */	\
            /* ================================================ */	\
            "psrldq $4,%%xmm8 \n\t" 				/* next prime */ \
            "movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,4 from xmm6 */ \
            "movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,4 from xmm7 */ \
            "movd	%%xmm8,%%edx \n\t"				/* else, extract prime from xmm8 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "jae    7f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_LOOP("3") \
            "7:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,4 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "jae    7f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_LOOP("3") \
            "7:		\n\t" \
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
            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "memory", "cc");

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

        logp = update_data.logp[j - 1];
        for (j = med_B; j < large_B; )
        {
            CHECK_NEW_SLICE(j);

            COMPUTE_4_PROOTS(j);

            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];
            nroot1 = (prime - root1);
            nroot2 = (prime - root2);
            prime = update_data.prime[j];

            FILL_ONE_PRIME_LOOP_P(j);
            FILL_ONE_PRIME_LOOP_N(j);

            j++;

            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];
            nroot1 = (prime - root1);
            nroot2 = (prime - root2);
            prime = update_data.prime[j];

            FILL_ONE_PRIME_LOOP_P(j);
            FILL_ONE_PRIME_LOOP_N(j);

            j++;

            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];
            nroot1 = (prime - root1);
            nroot2 = (prime - root2);
            prime = update_data.prime[j];

            FILL_ONE_PRIME_LOOP_P(j);
            FILL_ONE_PRIME_LOOP_N(j);

            j++;

            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];
            nroot1 = (prime - root1);
            nroot2 = (prime - root2);
            prime = update_data.prime[j];

            FILL_ONE_PRIME_LOOP_P(j);
            FILL_ONE_PRIME_LOOP_N(j);

            j++;
        }


#else

        logp = update_data.logp[j - 1];
        for (j = med_B; j < large_B; j++, ptr++)
        {
            CHECK_NEW_SLICE(j);

            prime = update_data.prime[j];
            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];

            COMPUTE_NEXT_ROOTS_P;

            update_data.firstroots1[j] = root1;
            update_data.firstroots2[j] = root2;
            nroot1 = prime - root1;
            nroot2 = prime - root2;

            FILL_ONE_PRIME_LOOP_P(j);
            FILL_ONE_PRIME_LOOP_N(j);
        }

#endif

#ifdef QS_TIMING
        gettimeofday(&qs_timing_stop, NULL);
        POLY_STG3 += yafu_difftime(&qs_timing_start, &qs_timing_stop);
        gettimeofday(&qs_timing_start, NULL);
#endif


#if defined(USE_POLY_SSE2_ASM) && defined(GCC_ASM64X) && !defined(PROFILING)
        logp = update_data.logp[large_B - 1];

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
        helperstruct.ptr = &rootupdates[(v - 1) * bound];  //72;						//72
        helperstruct.large_B = large_B;				//80
        helperstruct.B = bound;						//84
        helperstruct.interval = interval;			//88
        helperstruct.numblocks = numblocks;			//92
        helperstruct.bound_val = bound_val;			//96
        helperstruct.bound_index = bound_index;		//100
        helperstruct.check_bound = check_bound;		//104
        helperstruct.logp = logp;					//108

        ASM_G(\
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
            "movdqa %%xmm1, (%%r14,%%r15,4) \n\t"	/* save new root1 values */ \
            "paddd	%%xmm5, %%xmm2 \n\t"			/* selectively add back prime (modular subtract) */ \
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
            "movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
            "jae    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1("0") \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,1 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "psrldq $4,%%xmm2 \n\t" 				/* next root1 */ \
            "jae    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2("0") \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,2 ========== */	\
            /* ================================================ */	\
            "movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,2 from xmm1 */ \
            "movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,2 from xmm2 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
            "jae    3f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1("1") \
            "3:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,2 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "psrldq $4,%%xmm2 \n\t" 				/* next root1 */ \
            "jae    3f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2("1") \
            "3:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,3 ========== */	\
            /* ================================================ */	\
            "movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,3 from xmm1 */ \
            "movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,3 from xmm2 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
            "jae    5f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1("2") \
            "5:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,3 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "psrldq $4,%%xmm2 \n\t" 				/* next root1 */ \
            "jae    5f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2("2") \
            "5:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,4 ========== */	\
            /* ================================================ */	\
            "movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,4 from xmm1 */ \
            "movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,4 from xmm2 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "jae    7f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1("3") \
            "7:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,4 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "jae    7f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2("3") \
            "7:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,1 ========== */	\
            /* ================================================ */	\
            "movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
            "movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
            "movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,1 from xmm6 */ \
            "movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,1 from xmm7 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "psrldq $4,%%xmm6 \n\t" 				/* nextn root1 */ \
            "jae    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1("0") \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,1 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "psrldq $4,%%xmm7 \n\t" 				/* nextn nroot2 */ \
            "jae    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2("0") \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,2 ========== */	\
            /* ================================================ */	\
            "movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,2 from xmm6 */ \
            "movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,2 from xmm7 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "psrldq $4,%%xmm6 \n\t" 				/* nextn root1 */ \
            "jae    3f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1("1") \
            "3:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,2 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "psrldq $4,%%xmm7 \n\t" 				/* nextn nroot2 */ \
            "jae    3f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2("1") \
            "3:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,3 ========== */	\
            /* ================================================ */	\
            "movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,3 from xmm6 */ \
            "movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,3 from xmm7 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "psrldq $4,%%xmm6 \n\t" 				/* nextn root1 */ \
            "jae    5f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1("2") \
            "5:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,3 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "psrldq $4,%%xmm7 \n\t" 				/* nextn nroot2 */ \
            "jae    5f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2("2") \
            "5:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,4 ========== */	\
            /* ================================================ */	\
            "movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,4 from xmm6 */ \
            "movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,4 from xmm7 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "jae    7f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1("3") \
            "7:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,4 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "jae    7f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2("3") \
            "7:		\n\t" \
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
            : "rax", "rbx", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "memory", "cc");

        bound_index = helperstruct.bound_index;
        logp = helperstruct.logp;


#elif defined(HAS_SSE2)

        logp = update_data.logp[j - 1];
        for (j = large_B; j < bound; )
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

        logp = update_data.logp[j - 1];
        for (j = large_B; j < bound; j++, ptr++)
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
        gettimeofday(&qs_timing_stop, NULL);
        POLY_STG4 += yafu_difftime(&qs_timing_start, &qs_timing_stop);
#endif

    }
    else
    {

#ifdef QS_TIMING
        gettimeofday(&qs_timing_start, NULL);
#endif

        for (j = startprime; j < sconf->sieve_small_fb_start; j++, ptr++)
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
        for (j = sconf->sieve_small_fb_start;
            j < sconf->factor_base->fb_10bit_B; j++, ptr++)
        {
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

#if defined(HAS_SSE2) && (defined(GCC_ASM64X) || defined(_MSC_VER)) //NOTDEF //GCC_ASM64X
        // update 8 at a time using SSE2 and no branching		
        sm_ptr = &dconf->sm_rootupdates[(v - 1) * med_B];
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
            h.start = sconf->factor_base->fb_10bit_B;		// 64
            h.stop = sconf->factor_base->fb_15bit_B;		// 68
            if ((h.stop - 8) > h.start)
                h.stop -= 8;

            COMPUTE_8X_SMALL_NROOTS;

            j = h.stop;
        }
        sm_ptr = &dconf->sm_rootupdates[(v - 1) * med_B + j];


#else

        ptr = &dconf->rootupdates[(v - 1) * bound + sconf->factor_base->fb_10bit_B];
        for (j = sconf->factor_base->fb_10bit_B; j < sconf->factor_base->fb_15bit_B; j++, ptr++)
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
        ptr = &dconf->rootupdates[(v - 1) * bound + j];
        for (; j < med_B; j++, ptr++)
        {
            prime = update_data.prime[j];
            root1 = (uint16)update_data.sm_firstroots1[j];
            root2 = (uint16)update_data.sm_firstroots2[j];

            if ((prime > 32768) && ((j & 7) == 0))
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
        for (; j < med_B; j++, ptr++)
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

#ifdef QS_TIMING
        gettimeofday(&qs_timing_stop, NULL);
        POLY_STG2 += yafu_difftime(&qs_timing_start, &qs_timing_stop);
        gettimeofday(&qs_timing_start, NULL);
#endif

        bound_index = 0;
        bound_val = med_B;
        check_bound = med_B + BUCKET_ALLOC / 2;


#if defined(USE_POLY_SSE2_ASM) && defined(GCC_ASM64X) && !defined(PROFILING)
        logp = update_data.logp[med_B - 1];

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
        helperstruct.ptr = &rootupdates[(v - 1) * bound];  //72;						//72
        helperstruct.large_B = med_B;				//80
        helperstruct.B = large_B;					//84
        helperstruct.interval = interval;			//88
        helperstruct.numblocks = numblocks;			//92
        helperstruct.bound_val = bound_val;			//96
        helperstruct.bound_index = bound_index;		//100
        helperstruct.check_bound = check_bound;		//104
        helperstruct.logp = logp;					//108

        ASM_G(\
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
            "movdqa %%xmm1, (%%r14,%%r15,4) \n\t"	/* save new root1 values */ \
            "psubd	%%xmm5, %%xmm2 \n\t"			/* selectively sub back prime (modular addition) */ \
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
            "movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "movd	%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
            "jae    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_LOOP("0") \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,1 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "psrldq $4,%%xmm2 \n\t" 				/* next root1 */ \
            "jae    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_LOOP("0") \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,2 ========== */	\
            /* ================================================ */	\
            "psrldq $4,%%xmm0 \n\t" 				/* next prime */ \
            "movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,2 from xmm1 */ \
            "movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,2 from xmm2 */ \
            "movd	%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
            "jae    3f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_LOOP("1") \
            "3:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,2 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "psrldq $4,%%xmm2 \n\t" 				/* next root1 */ \
            "jae    3f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_LOOP("1") \
            "3:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,3 ========== */	\
            /* ================================================ */	\
            "psrldq $4,%%xmm0 \n\t" 				/* next prime */ \
            "movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,3 from xmm1 */ \
            "movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,3 from xmm2 */ \
            "movd	%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
            "jae    5f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_LOOP("2") \
            "5:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,3 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "psrldq $4,%%xmm2 \n\t" 				/* next root1 */ \
            "jae    5f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_LOOP("2") \
            "5:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,4 ========== */	\
            /* ================================================ */	\
            "psrldq $4,%%xmm0 \n\t" 				/* next prime */ \
            "movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,4 from xmm1 */ \
            "movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,4 from xmm2 */ \
            "movd	%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "jae    7f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_LOOP("3") \
            "7:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,4 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "jae    7f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_LOOP("3") \
            "7:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,1 ========== */	\
            /* ================================================ */	\
            "movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
            "movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
            "movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,1 from xmm6 */ \
            "movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,1 from xmm7 */ \
            "movd	%%xmm8,%%edx \n\t"				/* else, extract prime from xmm8 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "psrldq $4,%%xmm6 \n\t" 				/* nextn root1 */ \
            "jae    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_LOOP("0") \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,1 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "psrldq $4,%%xmm7 \n\t" 				/* nextn nroot2 */ \
            "jae    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_LOOP("0") \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,2 ========== */	\
            /* ================================================ */	\
            "psrldq $4,%%xmm8 \n\t" 				/* next prime */ \
            "movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,2 from xmm6 */ \
            "movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,2 from xmm7 */ \
            "movd	%%xmm8,%%edx \n\t"				/* else, extract prime from xmm8 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "psrldq $4,%%xmm6 \n\t" 				/* nextn root1 */ \
            "jae    3f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_LOOP("1") \
            "3:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,2 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "psrldq $4,%%xmm7 \n\t" 				/* nextn nroot2 */ \
            "jae    3f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_LOOP("1") \
            "3:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,3 ========== */	\
            /* ================================================ */	\
            "psrldq $4,%%xmm8 \n\t" 				/* next prime */ \
            "movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,3 from xmm6 */ \
            "movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,3 from xmm7 */ \
            "movd	%%xmm8,%%edx \n\t"				/* else, extract prime from xmm8 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "psrldq $4,%%xmm6 \n\t" 				/* nextn root1 */ \
            "jae    5f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_LOOP("2") \
            "5:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,3 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "psrldq $4,%%xmm7 \n\t" 				/* nextn nroot2 */ \
            "jae    5f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_LOOP("2") \
            "5:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,4 ========== */	\
            /* ================================================ */	\
            "psrldq $4,%%xmm8 \n\t" 				/* next prime */ \
            "movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,4 from xmm6 */ \
            "movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,4 from xmm7 */ \
            "movd	%%xmm8,%%edx \n\t"				/* else, extract prime from xmm8 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "jae    7f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1_LOOP("3") \
            "7:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,4 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "jae    7f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2_LOOP("3") \
            "7:		\n\t" \
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
            : "rax", "rbx", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "memory", "cc");

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

        logp = update_data.logp[j - 1];
        for (j = med_B; j < large_B; )
        {
            CHECK_NEW_SLICE(j);

            COMPUTE_4_NROOTS(j);

            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];
            nroot1 = (prime - root1);
            nroot2 = (prime - root2);
            prime = update_data.prime[j];

            FILL_ONE_PRIME_LOOP_P(j);
            FILL_ONE_PRIME_LOOP_N(j);

            j++;

            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];
            nroot1 = (prime - root1);
            nroot2 = (prime - root2);
            prime = update_data.prime[j];

            FILL_ONE_PRIME_LOOP_P(j);
            FILL_ONE_PRIME_LOOP_N(j);

            j++;

            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];
            nroot1 = (prime - root1);
            nroot2 = (prime - root2);
            prime = update_data.prime[j];

            FILL_ONE_PRIME_LOOP_P(j);
            FILL_ONE_PRIME_LOOP_N(j);

            j++;

            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];
            nroot1 = (prime - root1);
            nroot2 = (prime - root2);
            prime = update_data.prime[j];

            FILL_ONE_PRIME_LOOP_P(j);
            FILL_ONE_PRIME_LOOP_N(j);

            j++;
        }


#else

        logp = update_data.logp[j - 1];
        for (j = med_B; j < large_B; j++, ptr++)
        {
            CHECK_NEW_SLICE(j);

            prime = update_data.prime[j];
            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];

            COMPUTE_NEXT_ROOTS_N;

            update_data.firstroots1[j] = root1;
            update_data.firstroots2[j] = root2;
            nroot1 = prime - root1;
            nroot2 = prime - root2;

            FILL_ONE_PRIME_LOOP_P(j);
            FILL_ONE_PRIME_LOOP_N(j);
        }

#endif

#ifdef QS_TIMING
        gettimeofday(&qs_timing_stop, NULL);
        POLY_STG3 += yafu_difftime(&qs_timing_start, &qs_timing_stop);
        gettimeofday(&qs_timing_start, NULL);
#endif


#if defined(USE_POLY_SSE2_ASM) && defined(GCC_ASM64X) && !defined(PROFILING)
        logp = update_data.logp[large_B - 1];

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
        helperstruct.ptr = &rootupdates[(v - 1) * bound];  //72;						//72
        helperstruct.large_B = large_B;				//80
        helperstruct.B = bound;						//84
        helperstruct.interval = interval;			//88
        helperstruct.numblocks = numblocks;			//92
        helperstruct.bound_val = bound_val;			//96
        helperstruct.bound_index = bound_index;		//100
        helperstruct.check_bound = check_bound;		//104
        helperstruct.logp = logp;					//108

        ASM_G(\
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
            "movdqa %%xmm1, (%%r14,%%r15,4) \n\t"	/* save new root1 values */ \
            "psubd	%%xmm5, %%xmm2 \n\t"			/* selectively sub back prime (modular addition) */ \
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
            "movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
            "movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
            "jae    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1("0") \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,1 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "psrldq $4,%%xmm2 \n\t" 				/* next root1 */ \
            "jae    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2("0") \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,2 ========== */	\
            /* ================================================ */	\
            "movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,2 from xmm1 */ \
            "movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,2 from xmm2 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
            "jae    3f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1("1") \
            "3:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,2 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "psrldq $4,%%xmm2 \n\t" 				/* next root1 */ \
            "jae    3f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2("1") \
            "3:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,3 ========== */	\
            /* ================================================ */	\
            "movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,3 from xmm1 */ \
            "movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,3 from xmm2 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
            "jae    5f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1("2") \
            "5:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,3 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "psrldq $4,%%xmm2 \n\t" 				/* next root1 */ \
            "jae    5f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2("2") \
            "5:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT1,4 ========== */	\
            /* ================================================ */	\
            "movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,4 from xmm1 */ \
            "movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,4 from xmm2 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "jae    7f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1("3") \
            "7:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE POS BUCKET - ROOT2,4 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "jae    7f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2("3") \
            "7:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,1 ========== */	\
            /* ================================================ */	\
            "movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
            "movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
            "movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,1 from xmm6 */ \
            "movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,1 from xmm7 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "psrldq $4,%%xmm6 \n\t" 				/* nextn root1 */ \
            "jae    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1("0") \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,1 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "psrldq $4,%%xmm7 \n\t" 				/* nextn nroot2 */ \
            "jae    1f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2("0") \
            "1:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,2 ========== */	\
            /* ================================================ */	\
            "movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,2 from xmm6 */ \
            "movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,2 from xmm7 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "psrldq $4,%%xmm6 \n\t" 				/* nextn root1 */ \
            "jae    3f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1("1") \
            "3:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,2 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "psrldq $4,%%xmm7 \n\t" 				/* nextn nroot2 */ \
            "jae    3f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2("1") \
            "3:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,3 ========== */	\
            /* ================================================ */	\
            "movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,3 from xmm6 */ \
            "movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,3 from xmm7 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "psrldq $4,%%xmm6 \n\t" 				/* nextn root1 */ \
            "jae    5f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1("2") \
            "5:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,3 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "psrldq $4,%%xmm7 \n\t" 				/* nextn nroot2 */ \
            "jae    5f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2("2") \
            "5:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT1,4 ========== */	\
            /* ================================================ */	\
            "movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,4 from xmm6 */ \
            "movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,4 from xmm7 */ \
            "cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
            "jae    7f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT1("3") \
            "7:		\n\t" \
            /* ================================================ */	\
            /* =========== UPDATE NEG BUCKET - ROOT2,4 ========== */	\
            /* ================================================ */	\
            "cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
            "jae    7f \n\t" 						/* jump if CF = 1 */ \
            UPDATE_ROOT2("3") \
            "7:		\n\t" \
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
            : "rax", "rbx", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "memory", "cc");

        bound_index = helperstruct.bound_index;
        logp = helperstruct.logp;

#elif defined(HAS_SSE2)

        logp = update_data.logp[j - 1];
        for (j = large_B; j < bound; )
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

        logp = update_data.logp[j - 1];
        for (j = large_B; j < bound; j++, ptr++)
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
        gettimeofday(&qs_timing_stop, NULL);
        POLY_STG4 += yafu_difftime(&qs_timing_start, &qs_timing_stop);
#endif

    }

    if (lp_bucket_p->list != NULL)
    {
        lp_bucket_p->num_slices = bound_index + 1;
        lp_bucket_p->logp[bound_index] = logp;
    }

    return;
}

void nextRoots_32k_avx2_small(static_conf_t *sconf, dynamic_conf_t *dconf)
{
    //update the roots 
    sieve_fb_compressed *fb_p = dconf->comp_sieve_p;
    sieve_fb_compressed *fb_n = dconf->comp_sieve_n;
    int *rootupdates = dconf->rootupdates;
    int it;

    update_t update_data = dconf->update_data;

    uint32 startprime = 2;
    uint32 bound = sconf->factor_base->B;

    char v = dconf->curr_poly->nu[dconf->numB];
    char sign = dconf->curr_poly->gray[dconf->numB];
    int *ptr;
    uint16 *sm_ptr;

    uint32 med_B = sconf->factor_base->med_B;

    uint32 j;
    int k, numblocks;
    uint32 root1, root2, prime;

    uint8 logp = 0;
    polysieve_t helperstruct;

    numblocks = sconf->num_blocks;

    CLEAN_AVX2;

    k = 0;
    ptr = &rootupdates[(v - 1) * bound + startprime];

    if (sign > 0)
    {

        for (j = startprime; j<sconf->sieve_small_fb_start; j++, ptr++)
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
        // we can start doing things 16 at a time and be
        // sure we can use aligned moves (static_data_init).		
        for (j = sconf->sieve_small_fb_start;
            j < sconf->factor_base->fb_10bit_B; j++, ptr++)
        {
            // as soon as we are aligned, use more efficient AVX2 based methods...
            if ((j & 15) == 0)
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

        // update 16 at a time using AVX2
        sm_ptr = &dconf->sm_rootupdates[(v - 1) * med_B];
        
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
            h.stop = sconf->factor_base->med_B;		// 68
        
            COMPUTE_16X_SMALL_PROOTS_AVX2;
        
            j = h.stop;
        }

        /*
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
            h.start = sconf->factor_base->fb_10bit_B;		// 64
            h.stop = sconf->factor_base->fb_15bit_B;		// 68
            if ((h.stop - 8) > h.start)
                h.stop -= 8;

            COMPUTE_8X_SMALL_PROOTS;

            j = h.stop;
        }
        sm_ptr = &dconf->sm_rootupdates[(v - 1) * med_B + j];

        // assembly code may not get all the way to 15 bits since we 
        // do things in blocks of 8 there.  Make sure we are at the 15 bit
        // boundary before we switch to using update_data.firstroots1/2.
        // this should only run a few iterations, if any.
        ptr = &dconf->rootupdates[(v - 1) * bound + j];
        for (; j < med_B; j++, ptr++)
        {
            prime = update_data.prime[j];
            root1 = (uint16)update_data.sm_firstroots1[j];
            root2 = (uint16)update_data.sm_firstroots2[j];

            if ((prime > 32768) && ((j & 7) == 0))
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

        // continue one at a time once we exceed 15 bits, because the 8x SSE2
        // code has a hard time with unsigned 16 bit comparisons
        for (; j < med_B; j++, ptr++)
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
        */

    }
    else
    {
        /////////////////////////////////////////////////////////////////////////////////
        // sign < 0
        /////////////////////////////////////////////////////////////////////////////////

        for (j = startprime; j<sconf->sieve_small_fb_start; j++, ptr++)
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
        // we can start doing things 16 at a time and be
        // sure we can use aligned moves (static_data_init).	
        for (j = sconf->sieve_small_fb_start;
            j < sconf->factor_base->fb_10bit_B; j++, ptr++)
        {

            // as soon as we are aligned, use more efficient AVX2 based methods...
            if ((j & 15) == 0)
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

        // update 16 at a time using AVX2 and no branching		
        sm_ptr = &dconf->sm_rootupdates[(v - 1) * med_B];
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
            h.stop = sconf->factor_base->med_B;		// 68
        
            COMPUTE_16X_SMALL_NROOTS_AVX2;
        
            j = h.stop;
        }

        /*
        // update 8 at a time using SSE2 and no branching		
        sm_ptr = &dconf->sm_rootupdates[(v - 1) * med_B];
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
            h.start = sconf->factor_base->fb_10bit_B;		// 64
            h.stop = sconf->factor_base->fb_15bit_B;		// 68
            if ((h.stop - 8) > h.start)
                h.stop -= 8;

            COMPUTE_8X_SMALL_NROOTS;

            j = h.stop;
        }
        sm_ptr = &dconf->sm_rootupdates[(v - 1) * med_B + j];

        // assembly code may not get all the way to 15 bits since we 
        // do things in blocks of 8 there.  Make sure we are at the 15 bit
        // boundary before we switch to using update_data.firstroots1/2.
        // this should only run a few iterations, if any.
        ptr = &dconf->rootupdates[(v - 1) * bound + j];
        for (; j < med_B; j++, ptr++)
        {
            prime = update_data.prime[j];
            root1 = (uint16)update_data.sm_firstroots1[j];
            root2 = (uint16)update_data.sm_firstroots2[j];

            if ((prime > 32768) && ((j & 7) == 0))
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
        for (; j < med_B; j++, ptr++)
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
        */
    }

    CLEAN_AVX2;

    return;
}

#endif // USE_AVX2
