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


#if defined( USE_AVX2 ) && defined (GCC_ASM64X)

#include "qs.h"

#ifdef USE_YAFU_TDIV
#define DIVIDE_ONE_PRIME \
	do	\
    	{	\
		fb_offsets[++smooth_num] = i;	\
		zShortDiv32(tmp32, prime, tmp32);	\
    	} while (zShortMod32(tmp32, prime) == 0);
#define DIVIDE_ONE_PRIME_2(j) \
	do	\
        	{	\
		fb_offsets[++smooth_num] = (j);	\
		zShortDiv32(tmp32, fbc->prime[j], tmp32);	\
            	} while (zShortMod32(tmp32, fbc->prime[j]) == 0);
#else
#define DIVIDE_ONE_PRIME \
	do \
    	{						\
		fb_offsets[++smooth_num] = i;	\
		mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], prime); \
    	} while (mpz_tdiv_ui(dconf->Qvals[report_num], prime) == 0); 
#define DIVIDE_ONE_PRIME_2(j) \
	do \
        	{						\
		fb_offsets[++smooth_num] = (j);	\
		mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], fbc->prime[j]); \
            } while (mpz_tdiv_ui(dconf->Qvals[report_num], fbc->prime[j]) == 0); 
#endif

#ifdef USE_YAFU_TDIV
#define DIVIDE_RESIEVED_PRIME(j) \
    	while (zShortMod32(tmp32, fbc->prime[i+j]) == 0)	\
        	{	\
		fb_offsets[++smooth_num] = i+j;	\
		zShortDiv32(tmp32, fbc->prime[i+j], tmp32);	\
        	}
#define DIVIDE_RESIEVED_PRIME_2(j) \
        	while (zShortMod32(tmp32, fbc->prime[j]) == 0)	\
                    	{	\
		fb_offsets[++smooth_num] = j;	\
		zShortDiv32(tmp32, fbc->prime[j], tmp32);	\
                    	}
#else
#define DIVIDE_RESIEVED_PRIME(j) \
    	while (mpz_tdiv_ui(dconf->Qvals[report_num], fbc->prime[i+j]) == 0) \
        	{						\
		fb_offsets[++smooth_num] = i+j;	\
		mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], fbc->prime[i+j]);		\
        	}
#define DIVIDE_RESIEVED_PRIME_2(j) \
        while (mpz_tdiv_ui(dconf->Qvals[report_num], fbc->prime[j]) == 0) \
            {						\
	    fb_offsets[++smooth_num] = j;	\
	    mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], fbc->prime[j]);		\
            }
#endif

#define TDIV_MED_CLEAN asm volatile("emms");


#define MOD_CMP_8X_vec(xtra_bits)																		\
		__asm__ (																				\
			"vmovdqa (%1), %%xmm3 \n\t"		/* move in primes */							\
			"vpsubw	%%xmm1, %%xmm0, %%xmm4 \n\t"	/* BLOCKSIZE - block_loc */						\
			"vpaddw	(%2), %%xmm4, %%xmm4 \n\t"		/* apply corrections */							\
			"vmovdqa (%4), %%xmm6 \n\t"		/* move in root1s */							\
			"vpmulhuw	(%3), %%xmm4, %%xmm4 \n\t"	/* (unsigned) multiply by inverses */		\
			"vmovdqa (%5), %%xmm2 \n\t"		/* move in root2s */							\
			"vpsrlw	$" xtra_bits ", %%xmm4, %%xmm4 \n\t"		/* to get to total shift of 24/26/28 bits */			\
			"vpaddw	%%xmm3, %%xmm1, %%xmm7 \n\t"	/* add primes and block_loc */					\
			"vpmullw	%%xmm3, %%xmm4, %%xmm4 \n\t"	/* (signed) multiply by primes */				\
			"vpsubw	%%xmm0, %%xmm7, %%xmm7 \n\t"	/* substract blocksize */						\
			"vpaddw	%%xmm7, %%xmm4, %%xmm4 \n\t"	/* add in block_loc + primes - blocksize */		\
			"vpcmpeqw	%%xmm4, %%xmm6, %%xmm6 \n\t"	/* compare to root1s */						\
			"vpcmpeqw	%%xmm4, %%xmm2, %%xmm2 \n\t"	/* compare to root2s */						\
			"vpor	%%xmm6, %%xmm2, %%xmm2 \n\t"	/* combine compares */							\
			"vpmovmskb %%xmm2, %%r8 \n\t"		/* export to result */							\
            "andl   $0xaaaaaaaa,%%r8d   \n\t"   /* mask the bits we don't care about */ \
            "movl	%0,%%r11d		\n\t"		/* initialize count of set bits */ \
            "xorq	%%r10,%%r10		\n\t"		/* initialize bit scan offset */ \
            "1:			\n\t"					/* top of bit scan loop */ \
            "bsfl	%%r8d,%%ecx		\n\t"		/* put least significant set bit index into rcx */ \
            "jz 2f	\n\t"						/* jump out if zero (no hits).  high percentage. */ \
            "addl	%%ecx,%%r10d	\n\t"			/* add in the offset of this index */ \
            "movl   %%r10d,%%r9d \n\t" \
            "shrl   $1,%%r9d \n\t"   \
            "addl   %7,%%r9d \n\t"   \
            "movw	%%r9w, (%6, %%r11, 2) \n\t"		/* put the bit index into the output buffer */ \
            "shrq	%%cl,%%r8	\n\t"			/* shift the bit scan register up to the bit we just processed */ \
            "incl	%%r11d		\n\t"			/* increment the count of set bits */ \
            "incq	%%r10		\n\t"			/* increment the index */ \
            "shrq	$1, %%r8 \n\t"				/* clear the bit */ \
            "jmp 1b		\n\t"					/* loop if so */ \
            "2:		\n\t"						/*  */ \
            "movl	%%r11d, %0 \n\t"			/* return the count of set bits */ \
			: "+r" (tmp3)																	\
			: "r" (fbc->prime + i), "r" (fullfb_ptr->correction + i), \
				"r" (fullfb_ptr->small_inv + i), "r" (fbc->root1 + i), \
					"r" (fbc->root2 + i), "r"(buffer), "r"(i) \
			: "r9", "r8", "r10", "r11", "rcx", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "memory", "cc");


#define MOD_CMP_16X_vec_b(xtra_bits)																		\
		__asm__ (																				\
            /* top half... load and work on last 8 roots */ \
            "vmovdqa    16(%1), %%xmm3 \n\t"		        /* move in primes */							\
			"vpsubw	    %%ymm1, %%ymm0, %%ymm4 \n\t"	    /* BLOCKSIZE - block_loc */						\
			"vpaddw	    16(%2), %%xmm4, %%xmm4 \n\t"		/* apply corrections */							\
			"vmovdqa    16(%4), %%xmm6 \n\t"		        /* move in root1s */							\
			"vpmulhuw	16(%3), %%xmm4, %%xmm4 \n\t"	    /* (unsigned) multiply by inverses */		\
			"vmovdqa    16(%5), %%xmm2 \n\t"		        /* move in root2s */							\
			"vpsrlw	    $" xtra_bits ", %%xmm4, %%xmm4 \n\t"		/* to get to total shift of 24/26/28 bits */			\
			"vpaddw	    %%ymm3, %%ymm1, %%ymm7 \n\t"	    /* add primes and block_loc */					\
			"vpmullw	%%xmm3, %%xmm4, %%xmm4 \n\t"	    /* (signed) multiply by primes */				\
			"vpsubw	    %%ymm0, %%ymm7, %%ymm7 \n\t"	    /* substract blocksize */						\
			"vpaddw	    %%ymm7, %%ymm4, %%ymm4 \n\t"	    /* add in block_loc + primes - blocksize */		\
			"vpcmpeqw	%%ymm4, %%ymm6, %%ymm6 \n\t"	    /* compare to root1s */						\
			"vpcmpeqw	%%ymm4, %%ymm2, %%ymm2 \n\t"	    /* compare to root2s */						\
			"vpor	    %%ymm6, %%ymm2, %%ymm2 \n\t"	    /* combine compares */							\
			"vpmovmskb  %%ymm2, %%r9 \n\t"		            /* export to result */							\
            /* bottom half... load and work on first 8 roots */ \
			"vmovdqa    (%1), %%ymm3 \n\t"		            /* move in primes */							\
			"vpsubw	    %%ymm1, %%ymm0, %%ymm4 \n\t"	    /* BLOCKSIZE - block_loc */						\
			"vpaddw	    (%2), %%ymm4, %%ymm4 \n\t"		    /* apply corrections */							\
			"vmovdqa    (%4), %%ymm6 \n\t"		            /* move in root1s */							\
			"vpmulhuw	(%3), %%xmm4, %%xmm4 \n\t"	        /* (unsigned) multiply by inverses */		\
            "salq       $16,%%r9    \n\t"                   /* move to top half of 32-bit word */ \
			"vmovdqa    (%5), %%ymm2 \n\t"		            /* move in root2s */							\
			"vpsrlw	$" xtra_bits ", %%ymm4, %%ymm4 \n\t"		/* to get to total shift of 24/26/28 bits */			\
			"vpaddw	    %%ymm3, %%ymm1, %%ymm7 \n\t"	    /* add primes and block_loc */					\
			"vpmullw	%%xmm3, %%xmm4, %%xmm4 \n\t"	    /* (signed) multiply by primes */				\
			"vpsubw	    %%ymm0, %%ymm7, %%ymm7 \n\t"	    /* substract blocksize */						\
			"vpaddw	    %%ymm7, %%ymm4, %%ymm4 \n\t"	    /* add in block_loc + primes - blocksize */		\
			"vpcmpeqw	%%ymm4, %%ymm6, %%ymm6 \n\t"	    /* compare to root1s */						\
			"vpcmpeqw	%%ymm4, %%ymm2, %%ymm2 \n\t"	    /* compare to root2s */						\
			"vpor	    %%ymm6, %%ymm2, %%ymm2 \n\t"	    /* combine compares */							\
			"vpmovmskb  %%ymm2, %%r8 \n\t"		            /* export to result */							\
            /* search... look for bits that are set in r8  */ \
            /* and move the indices that are set to a temp buffer */ \
            "orq    %%r9, %%r8  \n\t"                       /* now has 16 comparisons (taking 32 bits) */ \
            "movl	%0,%%r11d		\n\t"		            /* initialize count of set bits */ \
            "xorq	%%r10,%%r10		\n\t"		            /* initialize bit scan offset */ \
            "andl   $0xaaaaaaaa,%%r8d   \n\t"               /* mask the bits we don't care about */ \
            "1:			\n\t"					            /* top of bit scan loop */ \
            "bsfl	%%r8d,%%ecx		\n\t"		            /* put least significant set bit index into rcx */ \
            "jz 2f	\n\t"						            /* jump out if zero (no hits).  high percentage. */ \
            "addl	%%ecx,%%r10d	\n\t"			        /* add in the offset of this index */ \
            "movl   %%r10d,%%r9d \n\t"                      /* copy offset so we can modify it */    \
            "shrl   $1,%%r9d \n\t"                          /* divide by two  */    \
            "addl   %7,%%r9d \n\t"                          /* and add it loop counter */        \
            "movw	%%r9w, (%6, %%r11, 2) \n\t"		        /* put the bit index into the output buffer */ \
            "shrq	%%cl,%%r8	\n\t"			            /* shift the bit scan register up to the bit we just processed */ \
            "incl	%%r11d		\n\t"			            /* increment the count of set bits */ \
            "incq	%%r10		\n\t"			            /* increment the index */ \
            "shrq	$1, %%r8 \n\t"				            /* clear the bit */ \
            "jmp 1b		\n\t"					            /* loop if so */ \
            "2:		\n\t"						            /*  */ \
            "movl	%%r11d, %0 \n\t"			            /* return the count of set bits */ \
            : "+r" (tmp3)																	\
            : "r" (fbc->prime + i), "r" (fullfb_ptr->correction + i), \
            "r" (fullfb_ptr->small_inv + i), "r" (fbc->root1 + i), \
            "r" (fbc->root2 + i), "r"(buffer), "r"(i) \
            : "r9", "r8", "r10", "r11", "rcx", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "memory", "cc");


// for some reason, using ymm registers with either multiply instruction
// causes other parts of yafu to significantly slow down.  no idea why.
// so we workaround it by extract/inserting the high parts and doing
// 2x the operations with xmm regs (which don't impact the speed elsewhere).
#define MOD_CMP_16X_vec(xtra_bits)																		\
		__asm__ (																				\
			"vmovdqa (%1), %%ymm3 \n\t"		/* move in primes */							\
			"vpsubw	%%ymm1, %%ymm0, %%ymm4 \n\t"	/* BLOCKSIZE - block_loc */						\
			"vpaddw	(%2), %%ymm4, %%ymm4 \n\t"		/* apply corrections */							\
			"vmovdqa (%4), %%ymm6 \n\t"		/* move in root1s */							\
			"vpmulhuw	(%3), %%xmm4, %%xmm5 \n\t"	/* low 8 words, (unsigned) multiply by inverses */		\
            "vpmulhuw	16(%3), %%xmm4, %%xmm4 \n\t"	/* high 8 words, (unsigned) multiply by inverses */		\
            "vinserti128   $1, %%xmm4, %%ymm5, %%ymm4 \n\t" /* combine low and high parts */ \
			"vmovdqa (%5), %%ymm2 \n\t"		/* move in root2s */							\
			"vpsrlw	$" xtra_bits ", %%ymm4, %%ymm4 \n\t"		/* to get to total shift of 24/26/28 bits */			\
			"vpaddw	%%ymm3, %%ymm1, %%ymm7 \n\t"	/* add primes and block_loc */					\
            "vextracti128  $1, %%ymm4, %%xmm5 \n\t" /* put high part of op1 into xmm5 */ \
            "vextracti128  $1, %%ymm3, %%xmm8 \n\t" /* put high part of op2 into xmm8 */ \
            "vpmullw	%%xmm3, %%xmm4, %%xmm4 \n\t"	/* (signed) multiply by primes */				\
			"vpmullw	%%xmm8, %%xmm5, %%xmm5 \n\t"	/* (signed) multiply by primes */				\
            "vinserti128   $1, %%xmm5, %%ymm4, %%ymm4 \n\t" /* combine low and high parts of result */ \
			"vpsubw	%%ymm0, %%ymm7, %%ymm7 \n\t"	/* substract blocksize */						\
			"vpaddw	%%ymm7, %%ymm4, %%ymm4 \n\t"	/* add in block_loc + primes - blocksize */		\
			"vpcmpeqw	%%ymm4, %%ymm6, %%ymm6 \n\t"	/* compare to root1s */						\
			"vpcmpeqw	%%ymm4, %%ymm2, %%ymm2 \n\t"	/* compare to root2s */						\
			"vpor	%%ymm6, %%ymm2, %%ymm2 \n\t"	/* combine compares */							\
			"vpmovmskb %%ymm2, %%r8 \n\t"		/* export to result */							\
            "andl   $0xaaaaaaaa,%%r8d   \n\t"   /* mask the bits we don't care about */ \
            "movl	%0,%%r11d		\n\t"		/* initialize count of set bits */ \
            "xorl	%%r10d,%%r10d		\n\t"		/* initialize bit scan offset */ \
            "1:			\n\t"					/* top of bit scan loop */ \
            "bsfl	%%r8d,%%ecx		\n\t"		/* put least significant set bit index into rcx */ \
            "jz 2f	\n\t"						/* jump out if zero (no hits).  high percentage. */ \
            "addl	%%ecx,%%r10d	\n\t"			/* add in the offset of this index */ \
            "movl   %%r10d,%%r9d \n\t" \
            "shrl   $1,%%r9d \n\t"   \
            "addl   %7,%%r9d \n\t"   \
            "movw	%%r9w, (%6, %%r11, 2) \n\t"		/* put the bit index into the output buffer */ \
            "shrl	%%cl,%%r8d	\n\t"			/* shift the bit scan register up to the bit we just processed */ \
            "incl	%%r11d		\n\t"			/* increment the count of set bits */ \
            "incl	%%r10d		\n\t"			/* increment the index */ \
            "shrl	$1, %%r8d \n\t"				/* clear the bit */ \
            "jmp 1b		\n\t"					/* loop if so */ \
            "2:		\n\t"						/*  */ \
            "movl	%%r11d, %0 \n\t"			/* return the count of set bits */ \
            : "+r" (tmp3)																	\
            : "r" (fbc->prime + i), "r" (fullfb_ptr->correction + i), \
            "r" (fullfb_ptr->small_inv + i), "r" (fbc->root1 + i), \
            "r" (fbc->root2 + i), "r"(buffer), "r"(i)	 \
            : "r9", "r8", "r10", "r11", "rcx", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "memory", "cc");

#define MOD_CMP_16X_vec_c(xtra_bits)																		\
		__asm__ (																				\
			"vmovdqa (%1), %%ymm3 \n\t"		/* move in primes */							\
			"vpsubw	%%ymm1, %%ymm0, %%ymm4 \n\t"	/* BLOCKSIZE - block_loc */						\
			"vpaddw	(%2), %%ymm4, %%ymm4 \n\t"		/* apply corrections */							\
			"vmovdqa (%4), %%ymm6 \n\t"		/* move in root1s */							\
			"vpmulhuw	(%3), %%ymm4, %%ymm4 \n\t"	/* low 8 words, (unsigned) multiply by inverses */		\
			"vmovdqa (%5), %%ymm2 \n\t"		/* move in root2s */							\
			"vpsrlw	$" xtra_bits ", %%ymm4, %%ymm4 \n\t"		/* to get to total shift of 24/26/28 bits */			\
			"vpaddw	%%ymm3, %%ymm1, %%ymm7 \n\t"	/* add primes and block_loc */					\
            "vpmullw	%%xmm3, %%xmm4, %%xmm4 \n\t"	/* (signed) multiply by primes */				\
			"vpsubw	%%ymm0, %%ymm7, %%ymm7 \n\t"	/* substract blocksize */						\
			"vpaddw	%%ymm7, %%ymm4, %%ymm4 \n\t"	/* add in block_loc + primes - blocksize */		\
			"vpcmpeqw	%%ymm4, %%ymm6, %%ymm6 \n\t"	/* compare to root1s */						\
			"vpcmpeqw	%%ymm4, %%ymm2, %%ymm2 \n\t"	/* compare to root2s */						\
			"vpor	%%ymm6, %%ymm2, %%ymm2 \n\t"	/* combine compares */							\
			"vpmovmskb %%ymm2, %%r8 \n\t"		/* export to result */							\
            "andl   $0xaaaaaaaa,%%r8d   \n\t"   /* mask the bits we don't care about */ \
            "movl	%0,%%r11d		\n\t"		/* initialize count of set bits */ \
            "xorl	%%r10d,%%r10d		\n\t"		/* initialize bit scan offset */ \
            "1:			\n\t"					/* top of bit scan loop */ \
            "bsfl	%%r8d,%%ecx		\n\t"		/* put least significant set bit index into rcx */ \
            "jz 2f	\n\t"						/* jump out if zero (no hits).  high percentage. */ \
            "addl	%%ecx,%%r10d	\n\t"			/* add in the offset of this index */ \
            "movl   %%r10d,%%r9d \n\t" \
            "shrl   $1,%%r9d \n\t"   \
            "addl   %7,%%r9d \n\t"   \
            "movw	%%r9w, (%6, %%r11, 2) \n\t"		/* put the bit index into the output buffer */ \
            "shrl	%%cl,%%r8d	\n\t"			/* shift the bit scan register up to the bit we just processed */ \
            "incl	%%r11d		\n\t"			/* increment the count of set bits */ \
            "incl	%%r10d		\n\t"			/* increment the index */ \
            "shrl	$1, %%r8d \n\t"				/* clear the bit */ \
            "jmp 1b		\n\t"					/* loop if so */ \
            "2:		\n\t"						/*  */ \
            "movl	%%r11d, %0 \n\t"			/* return the count of set bits */ \
            : "+r" (tmp3)																	\
            : "r" (fbc->prime + i), "r" (fullfb_ptr->correction + i), \
            "r" (fullfb_ptr->small_inv + i), "r" (fbc->root1 + i), \
            "r" (fbc->root2 + i), "r"(buffer), "r"(i)	 \
            : "r9", "r8", "r10", "r11", "rcx", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "memory", "cc");

#define MOD_INIT_16X												\
		__asm__ (														\
			"vmovdqa (%0), %%ymm0 \n\t"		/* move in BLOCKSIZE */ \
            "vmovdqa (%1), %%ymm1 \n\t"		/* move in block_loc */ \
			:														\
			: "r" (bl_sizes), "r" (bl_locs)		\
			: "xmm0", "xmm1");



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

#define DO_16X

void tdiv_medprimes_32k_avx2(uint8 parity, uint32 poly_id, uint32 bnum,
    static_conf_t *sconf, dynamic_conf_t *dconf)
{
    //we have flagged this sieve offset as likely to produce a relation
    //nothing left to do now but check and see.
    uint32 i;
    uint32 tmp, prime, root1, root2, report_num;
    uint32 bound10;
    uint32 bound12;
    uint32 bound13;
    int smooth_num;
    uint32 *fb_offsets;
    sieve_fb *fb;
    sieve_fb_compressed *fbc;
    fb_element_siqs *fullfb_ptr, *fullfb = sconf->factor_base->list;
    uint32 block_loc;
    uint16 buffer[16];
    uint32 tmp3 = 0;
    uint32 r;

    uint16 *bl_sizes = dconf->bl_sizes;
    uint16 *bl_locs = dconf->bl_locs;

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

    bl_sizes[0] = 32768;
    bl_sizes[1] = 32768;
    bl_sizes[2] = 32768;
    bl_sizes[3] = 32768;
    bl_sizes[4] = 32768;
    bl_sizes[5] = 32768;
    bl_sizes[6] = 32768;
    bl_sizes[7] = 32768;

#ifdef DO_16X
    bl_sizes[8] = 32768;
    bl_sizes[9] = 32768;
    bl_sizes[10] = 32768;
    bl_sizes[11] = 32768;
    bl_sizes[12] = 32768;
    bl_sizes[13] = 32768;
    bl_sizes[14] = 32768;
    bl_sizes[15] = 32768;
#endif


#ifdef DO_16X

    // 16x trial division
    if ((sconf->factor_base->fb_10bit_B & 15) == 0)
    {
        bound10 = sconf->factor_base->fb_10bit_B;
    }
    else
    {
        bound10 = MAX(sconf->factor_base->fb_10bit_B - 8, sconf->sieve_small_fb_start);
    }

    // determine the next bound
    if ((sconf->factor_base->fb_12bit_B & 15) == 0)
    {
        bound12 = sconf->factor_base->fb_12bit_B;
    }
    else
    {
        bound12 = MAX(sconf->factor_base->fb_12bit_B - 8, sconf->factor_base->fb_10bit_B);
    }

    // determine the next bound
    if ((sconf->factor_base->fb_13bit_B & 15) == 0)
    {
        bound13 = sconf->factor_base->fb_13bit_B;
    }
    else
    {
        bound13 = MAX(sconf->factor_base->fb_13bit_B - 8, sconf->factor_base->fb_12bit_B);
    }

#else
    bound10 = sconf->factor_base->fb_10bit_B;
    bound12 = sconf->factor_base->fb_12bit_B;
    bound13 = sconf->factor_base->fb_13bit_B;
#endif



    for (report_num = 0; report_num < dconf->num_reports; report_num++)
    {

#ifdef USE_YAFU_TDIV
        z32 *tmp32 = &dconf->Qvals32[report_num];
#endif

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

        CLEAN_AVX2;

        // pull the details of this report to get started.
        fb_offsets = &dconf->fb_offsets[report_num][0];
        smooth_num = dconf->smooth_num[report_num];
        block_loc = dconf->reports[report_num];

        i = sconf->sieve_small_fb_start;

        //printf("start index = %u\n", i);
        //printf("10-bit bound = %u\n", sconf->factor_base->fb_10bit_B);
        //printf("12-bit bound = %u\n", sconf->factor_base->fb_12bit_B);
        //printf("13-bit bound = %u\n", sconf->factor_base->fb_13bit_B);

        // single-up test until i is a multiple of 8
#ifdef DO_16X
        while ((i < bound10) && ((i & 15) != 0))
#else
        while ((i < bound10) && ((i & 7) != 0))
#endif
        {
            prime = fbc->prime[i];
            root1 = fbc->root1[i];
            root2 = fbc->root2[i];

            //tmp = distance from this sieve block offset to the end of the block
            tmp = 32768 - block_loc;

            //tmp = tmp/prime + 1 = number of steps to get past the end of the sieve
            //block, which is the state of the sieve now.
            tmp = 1 + (uint32)(((uint64)(tmp + fullfb_ptr->correction[i])
                * (uint64)fullfb_ptr->small_inv[i]) >> 24);
            tmp = block_loc + tmp*prime;
            tmp = tmp - 32768;

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


        tmp3 = 0;

        CLEAN_AVX2;


        bl_locs[0] = block_loc;
        bl_locs[1] = block_loc;
        bl_locs[2] = block_loc;
        bl_locs[3] = block_loc;
        bl_locs[4] = block_loc;
        bl_locs[5] = block_loc;
        bl_locs[6] = block_loc;
        bl_locs[7] = block_loc;

#ifdef DO_16X
        bl_locs[8] = block_loc;
        bl_locs[9] = block_loc;
        bl_locs[10] = block_loc;
        bl_locs[11] = block_loc;
        bl_locs[12] = block_loc;
        bl_locs[13] = block_loc;
        bl_locs[14] = block_loc;
        bl_locs[15] = block_loc;

        MOD_INIT_16X;
#else
        MOD_INIT_8X;
#endif



        while (i < bound10)
        {
#ifdef DO_16X
            MOD_CMP_16X_vec("8");
            i += 16;
#else
            MOD_CMP_8X_vec("8");
            i += 8;
#endif
        }


        // transition to beginning of 12-bit loop
        if (i != sconf->factor_base->fb_10bit_B)
        {
            if ((i + 16) < sconf->factor_base->fb_12bit_B)
            {
                MOD_CMP_8X_vec("8");
                i += 8;
                MOD_CMP_8X_vec("10");
                i += 8;
            }
            else
            {
                while ((uint32)i < sconf->factor_base->fb_12bit_B)
                {
                    prime = fbc->prime[i];
                    root1 = fbc->root1[i];
                    root2 = fbc->root2[i];

                    //tmp = distance from this sieve block offset to the end of the block
                    tmp = 32768 - block_loc;

                    //tmp = tmp/prime + 1 = number of steps to get past the end of the sieve
                    //block, which is the state of the sieve now.
                    tmp = 1 + (uint32)(((uint64)(tmp + fullfb_ptr->correction[i])
                        * (uint64)fullfb_ptr->small_inv[i]) >> 24);
                    tmp = block_loc + tmp*prime;
                    tmp = tmp - 32768;

                    //tmp = advance the offset to where it should be after the interval, and
                    //check to see if that's where either of the roots are now.  if so, then
                    //this offset is on the progression of the sieve for this prime
                    if (tmp == root1 || tmp == root2)
                    {
                        buffer[tmp3++] = i;
                    }
                    i++;
                }
            }
        }

        while (i < bound12)
        {
#ifdef DO_16X
            MOD_CMP_16X_vec("10");
            i += 16;
#else
            MOD_CMP_8X_vec("10");
            i += 8;
#endif
        }

        // transition to beginning of 13-bit loop
        if (i != sconf->factor_base->fb_12bit_B)
        {
            if ((i + 16) < sconf->factor_base->fb_13bit_B)
            {
                MOD_CMP_8X_vec("10");
                i += 8;
                MOD_CMP_8X_vec("12");
                i += 8;
            }
            else
            {
                while (i < sconf->factor_base->fb_13bit_B)
                {
                    prime = fbc->prime[i];
                    root1 = fbc->root1[i];
                    root2 = fbc->root2[i];

                    //tmp = distance from this sieve block offset to the end of the block
                    tmp = 32768 - block_loc;

                    //tmp = tmp/prime + 1 = number of steps to get past the end of the sieve
                    //block, which is the state of the sieve now.
                    tmp = 1 + (uint32)(((uint64)(tmp + fullfb_ptr->correction[i])
                        * (uint64)fullfb_ptr->small_inv[i]) >> 26);
                    tmp = block_loc + tmp*prime;
                    tmp = tmp - 32768;

                    //tmp = advance the offset to where it should be after the interval, and
                    //check to see if that's where either of the roots are now.  if so, then
                    //this offset is on the progression of the sieve for this prime
                    if (tmp == root1 || tmp == root2)
                    {
                        //it will divide Q(x).  do so as many times as we can.
                        buffer[tmp3++] = i;
                    }
                    i++;
                }
            }
        }

        while (i < bound13)
        {
#ifdef DO_16X
            MOD_CMP_16X_vec("12");
            i += 16;
#else
            MOD_CMP_8X_vec("12");
            i += 8;
#endif
        }

        // transition to beginning of 14-bit loop
        if (i != sconf->factor_base->fb_13bit_B)
        {
            if ((i + 8) <= sconf->factor_base->fb_13bit_B)
            {
                MOD_CMP_8X_vec("12");
            }
        }

        CLEAN_AVX2;


        for (r = 0; r < tmp3; r++)
        {
            DIVIDE_RESIEVED_PRIME_2((buffer[r]));
        }

        CLEAN_AVX2;

        // either after 8x SSE2 ASM, or standard trial division, record
        // how many factors we've found so far
        dconf->smooth_num[report_num] = smooth_num;

    }

    TDIV_MED_CLEAN;

    return;
}

#endif // USE_AVX2
