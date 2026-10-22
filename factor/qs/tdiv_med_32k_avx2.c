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


#if defined( USE_AVX2 )

#include "tdiv_macros_common.h"
#include "qs_impl.h"
#include <immintrin.h>


#if defined(GCC_ASM64X)

#define TDIV_MED_CLEAN __asm__ volatile("emms");


#define MOD_CMP_8X_vec(xtra_bits)																		\
		ASM_G (																				\
			"vmovdqa (%1), %%xmm3 \n\t"		/* move in primes */							\
			"vpsubw	%%xmm1, %%xmm0, %%xmm4 \n\t"	/* BLOCKSIZE - block_loc */						\
			"vpaddw	(%2), %%xmm4, %%xmm4 \n\t"		/* apply corrections */							\
			"vmovdqa (%4), %%xmm6 \n\t"		/* move in root1s */							\
			"vpmulhuw	(%3), %%xmm4, %%xmm4 \n\t"	/* (unsigned) multiply by inverses */		\
			"vmovdqa (%5), %%xmm2 \n\t"		/* move in root2s */							\
			"vpsrlw	$" STRING(xtra_bits) ", %%xmm4, %%xmm4 \n\t"		/* to get to total shift of 24/26/28 bits */			\
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


// for some reason, using ymm registers with either multiply instruction
// causes other parts of yafu to significantly slow down.  no idea why.
// so we workaround it by extract/inserting the high parts and doing
// 2x the operations with xmm regs (which don't impact the speed elsewhere).
#define MOD_CMP_16X_vec(xtra_bits)																		\
		ASM_G (																				\
			"vmovdqa (%1), %%ymm3 \n\t"		/* move in primes */							\
			"vpsubw	%%ymm1, %%ymm0, %%ymm4 \n\t"	/* BLOCKSIZE - block_loc */						\
			"vpaddw	(%2), %%ymm4, %%ymm4 \n\t"		/* apply corrections */							\
			"vmovdqa (%4), %%ymm6 \n\t"		/* move in root1s */							\
			"vpmulhuw	(%3), %%xmm4, %%xmm5 \n\t"	/* low 8 words, (unsigned) multiply by inverses */		\
            "vpmulhuw	16(%3), %%xmm4, %%xmm4 \n\t"	/* high 8 words, (unsigned) multiply by inverses */		\
            "vinserti128   $1, %%xmm4, %%ymm5, %%ymm4 \n\t" /* combine low and high parts */ \
			"vmovdqa (%5), %%ymm2 \n\t"		/* move in root2s */							\
			"vpsrlw	$" STRING(xtra_bits) ", %%ymm4, %%ymm4 \n\t"		/* to get to total shift of 24/26/28 bits */			\
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


#define MOD_INIT_16X												\
		ASM_G (														\
			"vmovdqa (%0), %%ymm0 \n\t"		/* move in BLOCKSIZE */ \
            "vmovdqa (%1), %%ymm1 \n\t"		/* move in block_loc */ \
			:														\
			: "r" (bl_sizes), "r" (bl_locs)		\
			: "xmm0", "xmm1");
#else

#define TDIV_MED_CLEAN

#define MOD_CMP_8X_vec(xtra_bits)				{ \
        __m128i v_primes, v_y4, v_y7, v_r1, v_r2;       \
        v_primes = _mm_load_si128((__m128i *)(fbc->prime + i)); \
        v_y4 = _mm_sub_epi16(v_blksz128, v_blkloc128); \
        v_y4 = _mm_add_epi16(v_y4, _mm_load_si128((__m128i *)(fullfb_ptr->correction + i))); \
        v_r1 = _mm_load_si128((__m128i *)(fbc->root1 + i)); \
        v_y4 = _mm_mulhi_epu16(v_y4, _mm_load_si128((__m128i *)(fullfb_ptr->small_inv + i))); \
        v_r2 = _mm_load_si128((__m128i *)(fbc->root2 + i)); \
        v_y4 = _mm_srli_epi16(v_y4, xtra_bits); \
        v_y7 = _mm_add_epi16(v_blkloc128, v_primes); \
        v_y4 = _mm_mullo_epi16(v_y4, v_primes); \
        v_y7 = _mm_sub_epi16(v_y7, v_blksz128); \
        v_y4 = _mm_add_epi16(v_y7, v_y4); \
        v_r1 = _mm_cmpeq_epi16(v_y4, v_r1); \
        v_r2 = _mm_cmpeq_epi16(v_y4, v_r2); \
        v_y4 = _mm_or_si128(v_r1, v_r2); \
        msk32 = _mm_movemask_epi8(v_y4); \
        msk32 &= 0xaaaaaaaa; \
        while (_BitScanForward(&pos, msk32)) { \
            buffer[tmp3++] = (pos >> 1) + i; \
            _reset_lsb(msk32); \
        }}

#define MOD_CMP_16X_vec(xtra_bits)				{ \
        __m256i v_primes, v_y4, v_y7, v_r1, v_r2;       \
        v_primes = _mm256_load_si256((__m256i *)(fbc->prime + i)); \
        v_y4 = _mm256_sub_epi16(v_blksz, v_blkloc); \
        v_y4 = _mm256_add_epi16(v_y4, _mm256_load_si256((__m256i *)(fullfb_ptr->correction + i))); \
        v_r1 = _mm256_load_si256((__m256i *)(fbc->root1 + i)); \
        v_y4 = _mm256_mulhi_epu16(v_y4, _mm256_load_si256((__m256i *)(fullfb_ptr->small_inv + i))); \
        v_r2 = _mm256_load_si256((__m256i *)(fbc->root2 + i)); \
        v_y4 = _mm256_srli_epi16(v_y4, xtra_bits); \
        v_y7 = _mm256_add_epi16(v_blkloc, v_primes); \
        v_y4 = _mm256_mullo_epi16(v_y4, v_primes); \
        v_y7 = _mm256_sub_epi16(v_y7, v_blksz); \
        v_y4 = _mm256_add_epi16(v_y7, v_y4); \
        v_r1 = _mm256_cmpeq_epi16(v_y4, v_r1); \
        v_r2 = _mm256_cmpeq_epi16(v_y4, v_r2); \
        v_y4 = _mm256_or_si256(v_r1, v_r2); \
        msk32 = _mm256_movemask_epi8(v_y4); \
        msk32 &= 0xaaaaaaaa; \
        while (_BitScanForward(&pos, msk32)) { \
            buffer[tmp3++] = (pos >> 1) + i; \
            _reset_lsb(msk32); \
        }}

#define MOD_INIT_16X									\
        __m256i v_blksz = _mm256_load_si256((__m256i *)bl_sizes);  \
        __m256i v_blkloc = _mm256_load_si256((__m256i *)bl_locs);  \
        __m128i v_blksz128 = _mm_load_si128((__m128i *)bl_sizes);  \
        __m128i v_blkloc128 = _mm_load_si128((__m128i *)bl_locs);  \
        uint32_t msk32, pos;


#endif


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

// Fast Division utilizing the precomputed structure
uint32_t divide(uint32_t n, uint32_t m, uint32_t c, uint32_t s) 
{
    // Multiplies two 32-bit numbers, returning only the upper 32 bits of the 64-bit result
    uint32_t high = ((uint64_t)n * m) >> 32;

    if (c) {
        // Handles precision adjustments for edge case divisors
        uint32_t t = (n - high) >> 1;
        return (high + t) >> s;
    }
    else {
        return high >> s;
    }
}

void tdiv_medprimes_32k_avx2(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t *sconf, dynamic_conf_t *dconf)
{
    // we have flagged this sieve offset as likely to produce a relation
    // nothing left to do now but check and see.
    uint32_t i;
    uint32_t tmp, prime, root1, root2, report_num;
    uint32_t bound10;
    uint32_t bound12;
    uint32_t bound13;
    int smooth_num;
    uint32_t *fb_offsets;
    sieve_fb_compressed *fbc;
    fb_element_siqs *fullfb_ptr, *fullfb = sconf->factor_base->list;
    uint32_t block_loc;
    uint16_t buffer[32];
    uint32_t tmp3 = 0;
    uint32_t r;

    uint16_t *bl_sizes = dconf->bl_sizes;
    uint16_t *bl_locs = dconf->bl_locs;

    fullfb_ptr = fullfb;
    if (parity)
    {
        fbc = dconf->comp_sieve_n;
    }
    else
    {
        fbc = dconf->comp_sieve_p;
    }

#ifdef USE_AVX512BWno

    // this code path works, but on benchmarks it is slower.  I suspect
    // either due to switching back and forth from 512 to 128 bit vectors
    // or from clock penalties in using 512 bit for very short periods
    // of time.
    bl_sizes[0] = 32768;
    bl_sizes[1] = 32768;
    bl_sizes[2] = 32768;
    bl_sizes[3] = 32768;
    bl_sizes[4] = 32768;
    bl_sizes[5] = 32768;
    bl_sizes[6] = 32768;
    bl_sizes[7] = 32768;
    bl_sizes[8] = 32768;
    bl_sizes[9] = 32768;
    bl_sizes[10] = 32768;
    bl_sizes[11] = 32768;
    bl_sizes[12] = 32768;
    bl_sizes[13] = 32768;
    bl_sizes[14] = 32768;
    bl_sizes[15] = 32768;
    bl_sizes[17] = 32768;
    bl_sizes[18] = 32768;
    bl_sizes[19] = 32768;
    bl_sizes[20] = 32768;
    bl_sizes[21] = 32768;
    bl_sizes[22] = 32768;
    bl_sizes[23] = 32768;
    bl_sizes[24] = 32768;
    bl_sizes[25] = 32768;
    bl_sizes[26] = 32768;
    bl_sizes[27] = 32768;
    bl_sizes[28] = 32768;
    bl_sizes[29] = 32768;
    bl_sizes[30] = 32768;
    bl_sizes[31] = 32768;

    int j;

    ALIGNED_MEM uint16_t bits_10_12[32] = { 8, 8, 8, 8, 8, 8, 8, 8,
        8, 8, 8, 8, 8, 8, 8, 8,
        8, 8, 8, 8, 8, 8, 8, 8,
        8, 8, 8, 8, 8, 8, 8, 8 };
    ALIGNED_MEM uint16_t bits_12_13[32] = { 10, 10, 10, 10, 10, 10, 10, 10,
        10, 10, 10, 10, 10, 10, 10, 10,
        10, 10, 10, 10, 10, 10, 10, 10,
        10, 10, 10, 10, 10, 10, 10, 10 };

    __m512i v_bits_10 = _mm512_set1_epi16(8);
    __m512i v_bits_12 = _mm512_set1_epi16(10);
    __m512i v_bits_13 = _mm512_set1_epi16(12);
    __m512i v_bits_10_12 = _mm512_set1_epi16(8);
    __m512i v_bits_12_13 = _mm512_set1_epi16(10);

    // 32x trial division
    i = sconf->sieve_small_fb_start;
    while ((i < sconf->factor_base->fb_10bit_B) && ((i & 31) != 0))
    {
        i++;
    }
    //printf("start = %u\n", i);

    while (i < sconf->factor_base->fb_10bit_B)
    {
        i += 32;
    }
    if (i > sconf->factor_base->fb_10bit_B)
    {
        i -= 32;
        bound10 = i;
        j = 0;
        while (i < sconf->factor_base->fb_10bit_B)
        {
            i += 8;
            j += 8;
        }
        while ((i < sconf->factor_base->fb_12bit_B) && ((i & 31) != 0))
        {
            i += 8;
        }
        for (; j < 32; j++)
        {
            bits_10_12[j] = 10;
        }
        v_bits_10_12 = _mm512_load_si512(bits_10_12);
    }
    else
    {
        bound10 = i;
    }

    // determine the next bound
    while (i < sconf->factor_base->fb_12bit_B)
    {
        i += 32;
    }
    if (i > sconf->factor_base->fb_12bit_B)
    {
        i -= 32;
        bound12 = i;
        j = 0;
        while (i < sconf->factor_base->fb_12bit_B)
        {
            i += 8;
            j += 8;
        }
        while ((i < sconf->factor_base->fb_13bit_B) && ((i & 31) != 0))
        {
            i += 8;
        }
        for (; j < 32; j++)
        {
            bits_12_13[j] = 12;
        }
        v_bits_12_13 = _mm512_load_si512(bits_12_13);
    }
    else
    {
        bound12 = i;
    }

    // determine the next bound
    while (i < sconf->factor_base->fb_13bit_B)
    {
        i += 32;
    }
    if (i > sconf->factor_base->fb_13bit_B)
    {
        i -= 32;
        bound13 = i;
        while (i < sconf->factor_base->fb_13bit_B)
        {
            i += 8;
        }
    }
    else
    {
        bound13 = i;
    }

    
    //printf("bound10 = %u\n", bound10);
    //printf("bound12 = %u\n", bound12);
    //printf("bound13 = %u\n", bound13);
    //printf("bound10_8steps = %u\n", bound10_8steps);
    //printf("bound12_8steps = %u\n", bound12_8steps);
    //printf("bound13_8steps = %u\n", bound13_8steps);
    //
    //printf("ptrs:\n\t%lp\n\t%lp\n\t%lp\n\t%lp\n\t%lp\n",
    //    fbc->prime,
    //    fbc->root1,
    //    fbc->root2,
    //    fullfb_ptr->correction,
    //    fullfb_ptr->small_inv);

#else
    bl_sizes[0] = 32768;
    bl_sizes[1] = 32768;
    bl_sizes[2] = 32768;
    bl_sizes[3] = 32768;
    bl_sizes[4] = 32768;
    bl_sizes[5] = 32768;
    bl_sizes[6] = 32768;
    bl_sizes[7] = 32768;
    bl_sizes[8] = 32768;
    bl_sizes[9] = 32768;
    bl_sizes[10] = 32768;
    bl_sizes[11] = 32768;
    bl_sizes[12] = 32768;
    bl_sizes[13] = 32768;
    bl_sizes[14] = 32768;
    bl_sizes[15] = 32768;


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

        CLEAN_AVX2;

        // pull the details of this report to get started.
        fb_offsets = &dconf->fb_offsets[report_num][0];
        smooth_num = dconf->smooth_num[report_num];
        block_loc = dconf->reports[report_num];

        i = sconf->sieve_small_fb_start;

#if 0

        for (; i < sconf->factor_base->fb_13bit_B; i++)
        {
            prime = fbc->prime[i];
            root1 = fbc->root1[i];
            root2 = fbc->root2[i];

            // distance from this sieve block offset to the end of the block
            tmp = 32768 - block_loc;
            tmp = 1 + divide(tmp, fullfb_ptr->inv32[i], fullfb_ptr->corr32[i], fullfb_ptr->shift32[i]);
            tmp = block_loc + tmp * prime;
            tmp = tmp - 32768;

            // tmp = advance the offset to where it should be after the interval, and
            // check to see if that's where either of the roots are now.  if so, then
            // this offset is on the progression of the sieve for this prime
            if (tmp == root1 || tmp == root2)
            {
                //it will divide Q(x).  do so as many times as we can.
                DIVIDE_ONE_PRIME;
            }
        }

        dconf->smooth_num[report_num] = smooth_num;

#elif 1
        
#ifdef USE_AVX512BWno

        // single-up test until i is a multiple of 32
        while ((i < bound10) && ((i & 31) != 0))
        {
            prime = fbc->prime[i];
            root1 = fbc->root1[i];
            root2 = fbc->root2[i];

            //tmp = distance from this sieve block offset to the end of the block
            tmp = 32768 - block_loc;

            //tmp = tmp/prime + 1 = number of steps to get past the end of the sieve
            //block, which is the state of the sieve now.
            tmp = 1 + (uint32_t)(((uint64_t)(tmp + fullfb_ptr->correction[i])
                * (uint64_t)fullfb_ptr->small_inv[i]) >> 24);
            tmp = block_loc + tmp * prime;
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
        bl_locs[8] = block_loc;
        bl_locs[9] = block_loc;
        bl_locs[10] = block_loc;
        bl_locs[11] = block_loc;
        bl_locs[12] = block_loc;
        bl_locs[13] = block_loc;
        bl_locs[14] = block_loc;
        bl_locs[15] = block_loc;
        bl_locs[16] = block_loc;
        bl_locs[17] = block_loc;
        bl_locs[18] = block_loc;
        bl_locs[19] = block_loc;
        bl_locs[20] = block_loc;
        bl_locs[21] = block_loc;
        bl_locs[22] = block_loc;
        bl_locs[23] = block_loc;
        bl_locs[24] = block_loc;
        bl_locs[25] = block_loc;
        bl_locs[26] = block_loc;
        bl_locs[27] = block_loc;
        bl_locs[28] = block_loc;
        bl_locs[29] = block_loc;
        bl_locs[30] = block_loc;
        bl_locs[31] = block_loc;

        MOD_INIT_32X;

        while (i < bound10)
        {
            MOD_CMP_32X_vec(v_bits_10);
            i += 32;
        }

        if (bound10 != sconf->factor_base->fb_10bit_B)
        {
            // transition to beginning of 12-bit loop
            MOD_CMP_32X_vec(v_bits_10_12);
            i += 32;
        }

        while (i < bound12)
        {
            MOD_CMP_32X_vec(v_bits_12);
            i += 32;
        }

        if (bound12 != sconf->factor_base->fb_12bit_B)
        {
            // transition to beginning of 13-bit loop
            MOD_CMP_32X_vec(v_bits_12_13);
            i += 32;
        }

        while (i < bound13)
        {
            MOD_CMP_32X_vec(v_bits_13);
            i += 32;
        }

        // transition to beginning of 14-bit loop
        while (i < sconf->factor_base->fb_13bit_B)
        {
            MOD_CMP_8X_vec_bw(12);
            i += 8;
        }

#else
        // single-up test until i is a multiple of 16
        while ((i < bound10) && ((i & 15) != 0))
        {
            prime = fbc->prime[i];
            root1 = fbc->root1[i];
            root2 = fbc->root2[i];

            //tmp = distance from this sieve block offset to the end of the block
            tmp = 32768 - block_loc;

            //tmp = tmp/prime + 1 = number of steps to get past the end of the sieve
            //block, which is the state of the sieve now.
            tmp = 1 + (uint32_t)(((uint64_t)(tmp + fullfb_ptr->correction[i])
                * (uint64_t)fullfb_ptr->small_inv[i]) >> 24);
            tmp = block_loc + tmp * prime;
            tmp = tmp - 32768;

            // tmp = advance the offset to where it should be after the interval, and
            // check to see if that's where either of the roots are now.  if so, then
            // this offset is on the progression of the sieve for this prime
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
        bl_locs[8] = block_loc;
        bl_locs[9] = block_loc;
        bl_locs[10] = block_loc;
        bl_locs[11] = block_loc;
        bl_locs[12] = block_loc;
        bl_locs[13] = block_loc;
        bl_locs[14] = block_loc;
        bl_locs[15] = block_loc;

        MOD_INIT_16X;

        while (i < bound10)
        {
            MOD_CMP_16X_vec(8);
            i += 16;
        }

        // transition to beginning of 12-bit loop
        if (i != sconf->factor_base->fb_10bit_B)
        {
            if ((i + 16) < sconf->factor_base->fb_12bit_B)
            {
                MOD_CMP_8X_vec(8);
                i += 8;
                MOD_CMP_8X_vec(10);
                i += 8;
            }
            else
            {
                while ((uint32_t)i < sconf->factor_base->fb_12bit_B)
                {
                    prime = fbc->prime[i];
                    root1 = fbc->root1[i];
                    root2 = fbc->root2[i];

                    //tmp = distance from this sieve block offset to the end of the block
                    tmp = 32768 - block_loc;

                    //tmp = tmp/prime + 1 = number of steps to get past the end of the sieve
                    //block, which is the state of the sieve now.
                    tmp = 1 + (uint32_t)(((uint64_t)(tmp + fullfb_ptr->correction[i])
                        * (uint64_t)fullfb_ptr->small_inv[i]) >> 24);
                    tmp = block_loc + tmp * prime;
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
            MOD_CMP_16X_vec(10);
            i += 16;
        }

        // transition to beginning of 13-bit loop
        if (i != sconf->factor_base->fb_12bit_B)
        {
            if ((i + 16) < sconf->factor_base->fb_13bit_B)
            {
                MOD_CMP_8X_vec(10);
                i += 8;
                MOD_CMP_8X_vec(12);
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
                    tmp = 1 + (uint32_t)(((uint64_t)(tmp + fullfb_ptr->correction[i])
                        * (uint64_t)fullfb_ptr->small_inv[i]) >> 26);
                    tmp = block_loc + tmp * prime;
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
            MOD_CMP_16X_vec(12);
            i += 16;
        }

        // transition to beginning of 14-bit loop
        if (i != sconf->factor_base->fb_13bit_B)
        {
            if ((i + 8) <= sconf->factor_base->fb_13bit_B)
            {
                MOD_CMP_8X_vec(12);
            }
        }

        CLEAN_AVX2;
#endif


        for (r = 0; r < tmp3; r++)
        {
            DIVIDE_RESIEVED_PRIME_2((buffer[r]));
        }

        CLEAN_AVX2;

        // either after 8x SSE2 ASM, or standard trial division, record
        // how many factors we've found so far
        dconf->smooth_num[report_num] = smooth_num;

#endif

    }

    TDIV_MED_CLEAN;

    return;
}

void tdiv_L2primes_32k_avx2(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf)
{
    // we have flagged this sieve offset as likely to produce a relation
    // nothing left to do now but check and see.
    uint32_t i;
    uint32_t tmp, prime, root1, root2, report_num;
    int smooth_num;
    uint32_t* fb_offsets;
    sieve_fb* fbc;
    fb_element_siqs* fullfb_ptr, * fullfb = sconf->factor_base->list;
    uint32_t block_loc;

    fullfb_ptr = fullfb;
    if (parity)
    {
        fbc = dconf->sieve_n;
    }
    else
    {
        fbc = dconf->sieve_p;
    }

    for (report_num = 0; report_num < dconf->num_reports; report_num++)
    {
        if (!dconf->valid_Qs[report_num])
            continue;

        // pull the details of this report to get started.
        fb_offsets = &dconf->fb_offsets[report_num][0];
        smooth_num = dconf->smooth_num[report_num];
        block_loc = dconf->reports[report_num] + bnum * 32768;

        for (i = sconf->factor_base->med_B; i < sconf->factor_base->fb_blk8i; i++)
        {
            prime = fbc->prime[i];
            root1 = fbc->root1[i];
            root2 = fbc->root2[i];

            // distance from this sieve block offset to the end of the L2-block
            tmp = 32768 * 8 - block_loc;
            tmp = 1 + divide(tmp, fullfb_ptr->inv32[i], fullfb_ptr->corr32[i], fullfb_ptr->shift32[i]);
            tmp = block_loc + tmp * prime;
            tmp = tmp - 32768 * 8;

            // tmp = advance the offset to where it should be after the interval, and
            // check to see if that's where either of the roots are now.  if so, then
            // this offset is on the progression of the sieve for this prime
            if (tmp == root1 || tmp == root2)
            {
                //it will divide Q(x).  do so as many times as we can.
                //DIVIDE_RESIEVED_PRIME_2(i);
                DIVIDE_ONE_PRIME;
            }
        }

        dconf->smooth_num[report_num] = smooth_num;
    }

    return;
}

// ----------------------------------------------------------------------
// trial division for the experimental linked-list sieve (TRY_LIST_SIEVE).
//
// Runs the same block-by-block list procedure as
// list_sieveblock_32k_avx512bw (see med_sieve_32k_avx2.c), but with its
// own independent per-block lists (ctx->rhead/rtail/rnext) and running
// roots (ctx->rroot1/2), built by list_sieve_init_32k from the pristine
// roots the sieve saved in ctx->oroot1/oroot2.  Instead of decrementing
// sieve bytes, each visited hit checks whether its position is a known
// smooth location (one of this block's reports) and if so does the gmp
// divisions.
//
// Because the sieve and the resieve are called in lockstep, block by
// block, the resieve's list for block bnum contains exactly the primes
// that just hit the block, so in the normal path this routine only walks
// that one list: no rebuild and no fast-forward.  A dormant state machine
// (ctx->rblock) remains for robustness: if a block's call is ever
// skipped or comes out of order, the lists are rebuilt from the pristine
// roots and the procedure is fast-forwarded through the missing blocks
// (advance-only traversal), because later blocks' hits depend on earlier
// blocks being traversed.
//
// Divisions are guarded: divisibility is checked before the first
// division of Q by each prime (no unchecked progression matches).
//
// Notes for the caller:
//  - call once per block from the scan path, after list_sieve_init_32k
//    has run for the current poly (it builds both procedures' lists).
//  - ctx must be the same context (same side, same thread) that was
//    passed to the sieve for this poly.
//  - replaces tdiv_L2primes/tdiv_large and the bucket tdiv for the
//    prime range covered by the list sieve.
// ----------------------------------------------------------------------
void resieve_large_32k_avx2(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf, list_sieve_ctx_t* ctx)
{
    uint32_t i, b, id;
    uint32_t prime, root, pos, rn, rpt;
    uint32_t node, next;
    uint32_t start = ctx->start;
    uint32_t interval = ctx->numblocks << 15;
    uint32_t bloc = bnum << 15;
    uint32_t* fb_offsets;
    int smooth_num;
    fb_element_siqs* listp = sconf->factor_base->list;

    (void)parity;   // the side is implicit in ctx
    (void)poly_id;

    // dormant robustness path: if we're being asked for a block earlier
    // than the state has already advanced past (or the state was never
    // built), rebuild the resieve lists from the pristine roots.  In the
    // normal lockstep flow this never triggers: list_sieve_init_32k
    // builds the lists and sets rblock = 0.
    if ((!ctx->resieve_ready) || (bnum < ctx->rblock))
    {
        for (b = 0; b < ctx->numblocks; b++)
        {
            ctx->rhead[b] = LISTSIEVE_NULL;
            ctx->rtail[b] = LISTSIEVE_NULL;
        }

        for (i = start; i < ctx->stop; i++)
        {
            id = i - start;
            ctx->rroot1[id] = ctx->oroot1[id];
            ctx->rroot2[id] = ctx->oroot2[id];

            if (ctx->rroot1[id] < interval)
                LIST_SIEVE_APPEND(ctx->rhead, ctx->rtail, ctx->rnext,
                    ctx->rroot1[id] >> 15, 2 * id);
            if (ctx->rroot2[id] < interval)
                LIST_SIEVE_APPEND(ctx->rhead, ctx->rtail, ctx->rnext,
                    ctx->rroot2[id] >> 15, 2 * id + 1);
        }

        ctx->resieve_ready = 1;
        ctx->rblock = 0;
    }

    // fast-forward the procedure through any blocks that were skipped
    // since the last call (advance roots and re-append hits, no tdiv)
    while (ctx->rblock < bnum)
    {
        uint32_t fblk = ctx->rblock;

        node = ctx->rhead[fblk];
        ctx->rhead[fblk] = LISTSIEVE_NULL;
        ctx->rtail[fblk] = LISTSIEVE_NULL;

        while (node != LISTSIEVE_NULL)
        {
            next = ctx->rnext[node];
            id = node >> 1;
            i = start + id;
            prime = listp->prime[i];

            if ((node & 1) == 0)
            {
                root = ctx->rroot1[id] + prime;
                ctx->rroot1[id] = root;
                if (root < interval)
                    LIST_SIEVE_APPEND(ctx->rhead, ctx->rtail, ctx->rnext,
                        root >> 15, node);
            }
            else
            {
                root = ctx->rroot2[id] + prime;
                ctx->rroot2[id] = root;
                if (root < interval)
                    LIST_SIEVE_APPEND(ctx->rhead, ctx->rtail, ctx->rnext,
                        root >> 15, node);
            }

            node = next;
        }

        ctx->rblock++;
    }

    // mark the sieve positions of this block's valid reports; the stored
    // value is report_num + 1 so 0 can mean "no report here"
    for (rn = 0; rn < dconf->num_reports; rn++)
    {
        if (dconf->valid_Qs[rn])
            ctx->flags[dconf->reports[rn]] = (uint16_t)(rn + 1);
    }

    // walk this block's resieve list: each node is a (prime, root) that
    // just hit this block in the sieve
    node = ctx->rhead[bnum];
    ctx->rhead[bnum] = LISTSIEVE_NULL;
    ctx->rtail[bnum] = LISTSIEVE_NULL;

    while (node != LISTSIEVE_NULL)
    {
        // capture the link before processing: processing may re-append
        // this node to a later block's list
        next = ctx->rnext[node];
        id = node >> 1;
        i = start + id;
        prime = listp->prime[i];

        if ((node & 1) == 0)
        {
            root = ctx->rroot1[id];
            pos = root - bloc;
            root += prime;
            ctx->rroot1[id] = root;
            if (root < interval)
                LIST_SIEVE_APPEND(ctx->rhead, ctx->rtail, ctx->rnext,
                    root >> 15, node);
        }
        else
        {
            root = ctx->rroot2[id];
            pos = root - bloc;
            root += prime;
            ctx->rroot2[id] = root;
            if (root < interval)
                LIST_SIEVE_APPEND(ctx->rhead, ctx->rtail, ctx->rnext,
                    root >> 15, node);
        }

        // the list invariant says pos is in [0, 32768); guard the flags
        // lookup anyway so a broken invariant can't corrupt memory
        if ((pos < 32768) && (ctx->flags[pos] > 0))
        {
            rpt = ctx->flags[pos] - 1;
            fb_offsets = &dconf->fb_offsets[rpt][0];
            smooth_num = dconf->smooth_num[rpt];

            //it will divide Q(x).  do so as many times as we can.
            while (mpz_tdiv_ui(dconf->Qvals[rpt], prime) == 0)
            {
                fb_offsets[++smooth_num] = i;
                mpz_tdiv_q_ui(dconf->Qvals[rpt], dconf->Qvals[rpt], prime);
            }

            dconf->smooth_num[rpt] = smooth_num;
        }

        node = next;
    }

    ctx->rblock = bnum + 1;

    // clear the marks we set
    for (rn = 0; rn < dconf->num_reports; rn++)
    {
        if (dconf->valid_Qs[rn])
            ctx->flags[dconf->reports[rn]] = 0;
    }

    return;
}

#endif // USE_AVX2
