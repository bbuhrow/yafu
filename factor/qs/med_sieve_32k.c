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
#include "sieve_macros_32k.h"

#if defined(_MSC_VER)
	#include <mmintrin.h>
#endif

#ifdef TARGET_KNC
#include <immintrin.h>
#endif

#ifdef USE_AVX512F
#include <immintrin.h>
#endif


typedef struct
{
	uint8 *sieve;					//0
	uint16 *primeptr;				//8
	uint16 *root1ptr;				//16
	uint16 *root2ptr;				//24
	uint16 *logptr;					//32
	uint32 startprime;				//40
	uint32 med_B;					//44
} helperstruct_t;

void med_sieveblock_32k(uint8 *sieve, sieve_fb_compressed *fb, fb_list *full_fb,
    uint32 start_prime, uint8 s_init)
{
    uint32 i;
    uint32 med_B;

    uint32 prime, root1, root2, tmp, stop;
    uint8 logp;

    helperstruct_t asm_input;

#if defined( TARGET_KNC ) || defined(USE_AVX512F)
    __m512i vlomask = _mm512_set1_epi32(0x000000ff);
    __m512i vhimask = _mm512_set1_epi32(0xffffff00);
    __m512i vlosieve1, vhisieve1, vlosieve2, vhisieve2;
    __m512i vpmul = _mm512_setr_epi32(
        0, 1, 2, 3, 
        4, 5, 6, 7, 
        8, 9, 10, 11, 
        12, 13, 14, 15);

    __m512i vblock = _mm512_set1_epi32(32768);
    __m512i vzero = _mm512_setzero_epi32();

#endif

    med_B = full_fb->med_B;

#ifdef QS_TIMING
    gettimeofday(&qs_timing_start, NULL);
#endif

    //initialize the block
    BLOCK_INIT;

#if defined(USE_AVX512F)
    // beyond 11 bit we don't need to loop... 16 steps takes
    // us past blocksize all at once when prime > 2048
    for (i = start_prime; i < full_fb->fb_13bit_B; i++)
    {	
        __m512i vprime, vroot1, vroot2, vlogp, vidx1, vidx2, vbyte, v16p;
        __mmask16 mask1, mask2;
        int steps1 = 0, steps2 = 0;

        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];
        logp = fb->logp[i];

        vprime = _mm512_set1_epi32(prime);
        vroot1 = _mm512_set1_epi32(root1);
        vroot2 = _mm512_set1_epi32(root2);
        vlogp = _mm512_set1_epi32(logp);
        v16p = _mm512_slli_epi32(vprime, 4);
        vidx2 = _mm512_mullo_epi32(vprime, vpmul);
        vidx1 = _mm512_add_epi32(vroot1, vidx2);
        vidx2 = _mm512_add_epi32(vroot2, vidx2);

        mask2 = _mm512_cmp_epu32_mask(vidx2, vblock, _MM_CMPINT_LT);

        // while all 16 steps hit the block, loop using only mask2 (larger offset):
        do 
        {
            vhisieve2 = _mm512_mask_i32gather_epi32(vzero, mask2, vidx2, sieve, _MM_SCALE_1);
            vhisieve1 = _mm512_mask_i32gather_epi32(vzero, mask2, vidx1, sieve, _MM_SCALE_1);
            vlosieve2 = _mm512_and_epi32(vhisieve2, vlomask);
            vlosieve1 = _mm512_and_epi32(vhisieve1, vlomask);
            vhisieve2 = _mm512_and_epi32(vhisieve2, vhimask);
            vhisieve1 = _mm512_and_epi32(vhisieve1, vhimask);
            vlosieve2 = _mm512_sub_epi32(vlosieve2, vlogp);
            vlosieve1 = _mm512_sub_epi32(vlosieve1, vlogp);
            vlosieve2 = _mm512_or_epi32(vhisieve2, _mm512_and_epi32(vlosieve2, vlomask));
            vlosieve1 = _mm512_or_epi32(vhisieve1, _mm512_and_epi32(vlosieve1, vlomask));
            _mm512_mask_i32scatter_epi32(sieve, mask2, vidx2, vlosieve2, _MM_SCALE_1);
            _mm512_mask_i32scatter_epi32(sieve, mask2, vidx1, vlosieve1, _MM_SCALE_1);

            vidx1 = _mm512_mask_add_epi32(vidx1, mask2, vidx1, v16p);
            vidx2 = _mm512_mask_add_epi32(vidx2, mask2, vidx2, v16p);
            steps1 += 16;            
            steps2 += 16;

            mask2 = _mm512_cmp_epu32_mask(vidx2, vblock, _MM_CMPINT_LT);
        } while (mask2 == 0xffff);

        // last iteration using separate mask1 and mask2
        mask1 = _mm512_cmp_epu32_mask(vidx1, vblock, _MM_CMPINT_LT);

        vhisieve2 = _mm512_mask_i32gather_epi32(vzero, mask2, vidx2, sieve, _MM_SCALE_1);
        vhisieve1 = _mm512_mask_i32gather_epi32(vzero, mask1, vidx1, sieve, _MM_SCALE_1);
        vlosieve2 = _mm512_and_epi32(vhisieve2, vlomask);
        vlosieve1 = _mm512_and_epi32(vhisieve1, vlomask);
        vhisieve2 = _mm512_and_epi32(vhisieve2, vhimask);
        vhisieve1 = _mm512_and_epi32(vhisieve1, vhimask);
        vlosieve2 = _mm512_sub_epi32(vlosieve2, vlogp);
        vlosieve1 = _mm512_sub_epi32(vlosieve1, vlogp);
        vlosieve2 = _mm512_or_epi32(vhisieve2, _mm512_and_epi32(vlosieve2, vlomask));
        vlosieve1 = _mm512_or_epi32(vhisieve1, _mm512_and_epi32(vlosieve1, vlomask));
        _mm512_mask_i32scatter_epi32(sieve, mask2, vidx2, vlosieve2, _MM_SCALE_1);
        _mm512_mask_i32scatter_epi32(sieve, mask1, vidx1, vlosieve1, _MM_SCALE_1);

        // test to see if lower offset (mask1) took more steps
        // and determine the new roots.
        steps1 += _mm_popcnt_u32(mask1);
        steps2 += _mm_popcnt_u32(mask2);
        if (steps1 > steps2)
        {
            fb->root2[i] = (uint16)(steps1*prime - 32768);
            fb->root1[i] = (uint16)(steps2*prime - 32768);
        }
        else
        {
            fb->root1[i] = (uint16)(steps1*prime - 32768);
            fb->root2[i] = (uint16)(steps2*prime - 32768);
        }


    }

#elif defined(USE_ASM_SMALL_PRIME_SIEVING)

	asm_input.logptr = fb->logp;
	asm_input.primeptr = fb->prime;
	asm_input.root1ptr = fb->root1;
	asm_input.root2ptr = fb->root2;
	asm_input.sieve = sieve;
	asm_input.startprime = start_prime;
	asm_input.med_B = full_fb->fb_13bit_B-8;

	SIEVE_13b_ASM;

	i = asm_input.startprime;

#else
	for (i=start_prime;i< full_fb->fb_13bit_B-8;i++)
	{	
		uint8 *s2;		

		prime = fb->prime[i];
		root1 = fb->root1[i];
		root2 = fb->root2[i];
		logp = fb->logp[i];

		SIEVE_2X;
		SIEVE_1X;
		SIEVE_LAST;
		UPDATE_ROOTS;
	}
#endif

#if defined(HAS_SSE2) && !defined(TARGET_KNC)

	for (; i<med_B; i++)
	{	
		uint8 *s2;		

		prime = fb->prime[i];
		root1 = fb->root1[i];
		root2 = fb->root2[i];
		logp = fb->logp[i];

		// special exit condition: when prime > 8192 and i % 8 is 0;
		if ((prime > 8192) && ((i&7) == 0))
			break;

		// invalid root (part of poly->a)
		if (prime == 0) 
			continue;

		SIEVE_2X;
		SIEVE_1X;
		SIEVE_LAST;
		UPDATE_ROOTS;
	}

	_INIT_SSE2_SMALL_PRIME_SIEVE;
	_SSE2_SMALL_PRIME_SIEVE_32k_DIV3;

	// get past the 32k/3 boundary
	for (; i < full_fb->fb_32k_div3; i++)
	{	
		prime = fb->prime[i];
		root1 = fb->root1[i];
		root2 = fb->root2[i];
		logp = fb->logp[i];

		// invalid root (part of poly->a)
		if (prime == 0) 
			continue;

		SIEVE_1X;
		SIEVE_LAST;

		UPDATE_ROOTS;
	}

	_INIT_SSE2_SMALL_PRIME_SIEVE;
	_SSE2_SMALL_PRIME_SIEVE_14b;

	// get past the 14b boundary
	for (; i < full_fb->fb_14bit_B; i++)
	{	
		prime = fb->prime[i];
		root1 = fb->root1[i];
		root2 = fb->root2[i];
		logp = fb->logp[i];

		// invalid root (part of poly->a)
		if (prime == 0) 
			continue;

		SIEVE_1X;
		SIEVE_LAST;

		UPDATE_ROOTS;
	}

	_INIT_SSE2_SMALL_PRIME_SIEVE;
	_SSE2_SMALL_PRIME_SIEVE_15b;

	// get past the 15b boundary
	for (i=full_fb->fb_15bit_B-8;i<med_B;i++)
	{	
		prime = fb->prime[i];
		root1 = fb->root1[i];
		root2 = fb->root2[i];
		logp = fb->logp[i];

		// invalid root (part of poly->a)
		if (prime == 0) 
			continue;

		CHECK_1X_DONE;

		SIEVE_1X;
		SIEVE_LAST;

		UPDATE_ROOTS;
	}

#if defined(USE_ASM_SMALL_PRIME_SIEVING)

	asm_input.logptr = fb->logp;
	asm_input.primeptr = fb->prime;
	asm_input.root1ptr = fb->root1;
	asm_input.root2ptr = fb->root2;
	asm_input.sieve = sieve;
	asm_input.startprime = i;
	asm_input.med_B = med_B;

	SIEVE_GT_BLOCKSIZE_ASM;

#else

	//if there are primes left bigger than the blocksize, this will take
	//care of them.  if not, it doesn't run at all.
	for (;i<med_B;i++)
	{	
		prime = fb->prime[i];
		root1 = fb->root1[i];
		root2 = fb->root2[i];
		logp = fb->logp[i];

		SIEVE_BIG;
		UPDATE_ROOTS;
	}
#endif

#else


#if defined(USE_ASM_SMALL_PRIME_SIEVING)

#ifdef TARGET_KNC

    for (; i<full_fb->fb_15bit_B; i++)
    {	
        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];
        logp = fb->logp[i];

        if ((prime > 8192) && ((i & 15) == 0))
            break;

        // invalid root (part of poly->a)
        if (prime == 0)
            continue;

        SIEVE_1X;
        SIEVE_LAST;
        UPDATE_ROOTS;
    }

    logp = 13;
    if (((full_fb->fb_32k_div3 - 8 - i) & 15) != 0)
        stop = full_fb->fb_32k_div3 - 16;
    else
        stop = full_fb->fb_32k_div3 - 8;

    for (; i < stop; i += 16)
    {	
        __m512i vprime, vroot1, vroot2, vmax, vmin;
        __mmask16 mask1, mask2;
        int j, idx;

        vprime = _mm512_extload_epi32((__m512i *)(&fb->prime[i]),
            _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST32_NONE, _MM_HINT_NONE);
        vroot1 = _mm512_extload_epi32((__m512i *)(&fb->root1[i]),
            _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST32_NONE, _MM_HINT_NONE);
        vroot2 = _mm512_extload_epi32((__m512i *)(&fb->root2[i]),
            _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST32_NONE, _MM_HINT_NONE);

#pragma unroll(16)
        for (j = 0; j < 16; j++)
        {
            sieve[fb->root1[i+j]] -= logp;
            sieve[fb->root2[i+j]] -= logp;
        }

        vroot1 = _mm512_add_epi32(vroot1, vprime);
        vroot2 = _mm512_add_epi32(vroot2, vprime);

        _mm512_extstore_epi32((__m512i *)(&fb->root1[i]), vroot1,
            _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);
        _mm512_extstore_epi32((__m512i *)(&fb->root2[i]), vroot2,
            _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);

#pragma unroll(16)
        for (j = 0; j < 16; j++)
        {
            sieve[fb->root1[i+j]] -= logp;
            sieve[fb->root2[i+j]] -= logp;
        }

        vroot1 = _mm512_add_epi32(vroot1, vprime);
        vroot2 = _mm512_add_epi32(vroot2, vprime);

        _mm512_extstore_epi32((__m512i *)(&fb->root1[i]), vroot1,
            _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);
        _mm512_extstore_epi32((__m512i *)(&fb->root2[i]), vroot2,
            _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);

#pragma unroll(16)
        for (j = 0; j < 16; j++)
        {
            sieve[fb->root1[i+j]] -= logp;
            sieve[fb->root2[i+j]] -= logp;
        }

        vroot1 = _mm512_add_epi32(vroot1, vprime);
        vroot2 = _mm512_add_epi32(vroot2, vprime);

        mask1 = _mm512_cmp_epu32_mask(vroot1, vblock, _MM_CMPINT_LT);
        mask2 = _mm512_cmp_epu32_mask(vroot2, vblock, _MM_CMPINT_LT);

        _mm512_extstore_epi32((__m512i *)(&fb->root1[i]), vroot1,
            _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);
        _mm512_extstore_epi32((__m512i *)(&fb->root2[i]), vroot2,
            _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);

        vroot1 = _mm512_mask_add_epi32(vroot1, mask1, vroot1, vprime);
        vroot2 = _mm512_mask_add_epi32(vroot2, mask2, vroot2, vprime);

        while ((idx = _mm_tzcnt_32(mask1)) < 16)
        {
            sieve[fb->root1[i+idx]] -= logp;
            mask1 ^= (1 << idx);
            if (mask2 & (1 << idx))
                sieve[fb->root2[i+idx]] -= logp;
        }

        vroot1 = _mm512_sub_epi32(vroot1, vblock);
        vroot2 = _mm512_sub_epi32(vroot2, vblock);

        vmax = _mm512_max_epu32(vroot1, vroot2);
        vmin = _mm512_min_epu32(vroot1, vroot2);

        _mm512_extstore_epi32((__m512i *)(&fb->root1[i]), vmin,
            _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);
        _mm512_extstore_epi32((__m512i *)(&fb->root2[i]), vmax,
            _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);
    }


    for (; i<full_fb->fb_15bit_B; i++)
    {	
        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];
        logp = fb->logp[i];

        if ((prime > 10923) && ((i & 15) == 0))
            break;

        // invalid root (part of poly->a)
        if (prime == 0)
            continue;

        SIEVE_1X;
        SIEVE_LAST;
        UPDATE_ROOTS;
    }

    logp = 14;
    if (((full_fb->fb_14bit_B - 8 - i) & 15) != 0)
        stop = full_fb->fb_14bit_B - 16;
    else
        stop = full_fb->fb_14bit_B - 8;
   
    for (; i < stop; i += 16)
    {	
        __m512i vprime, vroot1, vroot2, vmax, vmin;
        __mmask16 mask1, mask2;
        int j, idx;

        vprime = _mm512_extload_epi32((__m512i *)(&fb->prime[i]),
            _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST32_NONE, _MM_HINT_NONE);
        vroot1 = _mm512_extload_epi32((__m512i *)(&fb->root1[i]),
            _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST32_NONE, _MM_HINT_NONE);
        vroot2 = _mm512_extload_epi32((__m512i *)(&fb->root2[i]),
            _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST32_NONE, _MM_HINT_NONE);

#pragma unroll(16)
        for (j = 0; j < 16; j++)
        {
            sieve[fb->root1[i+j]] -= logp;
            sieve[fb->root2[i+j]] -= logp;
        }

        vroot1 = _mm512_add_epi32(vroot1, vprime);
        vroot2 = _mm512_add_epi32(vroot2, vprime);

        _mm512_extstore_epi32((__m512i *)(&fb->root1[i]), vroot1,
            _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);
        _mm512_extstore_epi32((__m512i *)(&fb->root2[i]), vroot2,
            _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);

#pragma unroll(16)
        for (j = 0; j < 16; j++)
        {
            sieve[fb->root1[i+j]] -= logp;
            sieve[fb->root2[i+j]] -= logp;
        }

        vroot1 = _mm512_add_epi32(vroot1, vprime);
        vroot2 = _mm512_add_epi32(vroot2, vprime);

        mask1 = _mm512_cmp_epu32_mask(vroot1, vblock, _MM_CMPINT_LT);
        mask2 = _mm512_cmp_epu32_mask(vroot2, vblock, _MM_CMPINT_LT);

        _mm512_extstore_epi32((__m512i *)(&fb->root1[i]), vroot1,
            _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);
        _mm512_extstore_epi32((__m512i *)(&fb->root2[i]), vroot2,
            _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);

        vroot1 = _mm512_mask_add_epi32(vroot1, mask1, vroot1, vprime);
        vroot2 = _mm512_mask_add_epi32(vroot2, mask2, vroot2, vprime);

        while ((idx = _mm_tzcnt_32(mask1)) < 16)
        {
            sieve[fb->root1[i+idx]] -= logp;
            mask1 ^= (1 << idx);
            if (mask2 & (1 << idx))
                sieve[fb->root2[i+idx]] -= logp;
        }

        vroot1 = _mm512_sub_epi32(vroot1, vblock);
        vroot2 = _mm512_sub_epi32(vroot2, vblock);

        vmax = _mm512_max_epu32(vroot1, vroot2);
        vmin = _mm512_min_epu32(vroot1, vroot2);

        _mm512_extstore_epi32((__m512i *)(&fb->root1[i]), vmin,
            _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);
        _mm512_extstore_epi32((__m512i *)(&fb->root2[i]), vmax,
            _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);
    }

    for (; i<full_fb->fb_15bit_B; i++)
    {	
        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];

        if ((prime > 16384) && ((i & 15) == 0))
            break;

        // invalid root (part of poly->a)
        if (prime == 0)
            continue;

        SIEVE_1X;
        SIEVE_LAST;
        UPDATE_ROOTS;
    }

    logp = 15;
    for (; i<full_fb->fb_15bit_B; i += 16)
    {	
        __m512i vprime, vroot1, vroot2, vmax, vmin;
        __mmask16 mask1, mask2;
        int j, idx;

        vprime = _mm512_extload_epi32((__m512i *)(&fb->prime[i]),
            _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST32_NONE, _MM_HINT_NONE);
        vroot1 = _mm512_extload_epi32((__m512i *)(&fb->root1[i]),
            _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST32_NONE, _MM_HINT_NONE);
        vroot2 = _mm512_extload_epi32((__m512i *)(&fb->root2[i]),
            _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST32_NONE, _MM_HINT_NONE);

#pragma unroll(16)
        for (j = 0; j < 16; j++)
        {
            sieve[fb->root1[i+j]] -= logp;
            sieve[fb->root2[i+j]] -= logp;
        }

        vroot1 = _mm512_add_epi32(vroot1, vprime);
        vroot2 = _mm512_add_epi32(vroot2, vprime);

        mask1 = _mm512_cmp_epu32_mask(vroot1, vblock, _MM_CMPINT_LT);
        mask2 = _mm512_cmp_epu32_mask(vroot2, vblock, _MM_CMPINT_LT);

        _mm512_extstore_epi32((__m512i *)(&fb->root1[i]), vroot1,
            _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);
        _mm512_extstore_epi32((__m512i *)(&fb->root2[i]), vroot2,
            _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);

        vroot1 = _mm512_mask_add_epi32(vroot1, mask1, vroot1, vprime);
        vroot2 = _mm512_mask_add_epi32(vroot2, mask2, vroot2, vprime);

        while ((idx = _mm_tzcnt_32(mask1)) < 16)
        {
            sieve[fb->root1[i+idx]] -= logp;
            mask1 ^= (1 << idx);
            if (mask2 & (1 << idx))
                sieve[fb->root2[i+idx]] -= logp;
        }

        vroot1 = _mm512_sub_epi32(vroot1, vblock);
        vroot2 = _mm512_sub_epi32(vroot2, vblock);

        vmax = _mm512_max_epu32(vroot1, vroot2);
        vmin = _mm512_min_epu32(vroot1, vroot2);

        _mm512_extstore_epi32((__m512i *)(&fb->root1[i]), vmin,
            _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);
        _mm512_extstore_epi32((__m512i *)(&fb->root2[i]), vmax,
            _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);

    }

    logp = 15;
    for (; i<med_B; i += 16)
    {	
        __m512i vprime, vroot1, vroot2, vmax, vmin;
        __mmask16 mask1, mask2;
        int j, idx;

        vprime = _mm512_extload_epi32((__m512i *)(&fb->prime[i]),
            _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST32_NONE, _MM_HINT_NONE);
        vroot1 = _mm512_extload_epi32((__m512i *)(&fb->root1[i]),
            _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST32_NONE, _MM_HINT_NONE);
        vroot2 = _mm512_extload_epi32((__m512i *)(&fb->root2[i]),
            _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST32_NONE, _MM_HINT_NONE);

        mask1 = _mm512_cmp_epu32_mask(vroot1, vblock, _MM_CMPINT_LT);
        mask2 = _mm512_cmp_epu32_mask(vroot2, vblock, _MM_CMPINT_LT);

        vroot1 = _mm512_mask_add_epi32(vroot1, mask1, vroot1, vprime);
        vroot2 = _mm512_mask_add_epi32(vroot2, mask2, vroot2, vprime);

        while ((idx = _mm_tzcnt_32(mask1)) < 16)
        {
            sieve[fb->root1[i+idx]] -= logp;
            mask1 ^= (1 << idx);
            if (mask2 & (1 << idx))
                sieve[fb->root2[i+idx]] -= logp;
        }

        vroot1 = _mm512_sub_epi32(vroot1, vblock);
        vroot2 = _mm512_sub_epi32(vroot2, vblock);

        vmax = _mm512_max_epu32(vroot1, vroot2);
        vmin = _mm512_min_epu32(vroot1, vroot2);

        _mm512_extstore_epi32((__m512i *)(&fb->root1[i]), vmin,
            _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);
        _mm512_extstore_epi32((__m512i *)(&fb->root2[i]), vmax,
            _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);
    }

#else

    for (; i<full_fb->fb_15bit_B; i++)
    {	
        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];
        logp = fb->logp[i];

        if (prime > 8192)
            break;

        // invalid root (part of poly->a)
        if (prime == 0)
            continue;

        SIEVE_1X;
        SIEVE_LAST;
        UPDATE_ROOTS;
    }

    logp = 13;
    for (; i<full_fb->fb_32k_div3-8; i++)
    {	
        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];

        sieve[root1] -= logp;
        sieve[root2] -= logp;
        root1 += prime;
        root2 += prime;

        sieve[root1] -= logp;
        sieve[root2] -= logp;
        root1 += prime;
        root2 += prime;

        sieve[root1] -= logp;
        sieve[root2] -= logp;
        root1 += prime;
        root2 += prime;

        if (root1 < 32768) 
        {
            sieve[root1] -= logp;
            root1 += prime;
            if (root2 < 32768)
            {
                sieve[root2] -= logp;
                root2 += prime;
            }
            else
            {
                tmp = root2;
                root2 = root1;
                root1 = tmp;
            }
        }

        UPDATE_ROOTS;
    }

    for (; i<full_fb->fb_15bit_B; i++)
    {	
        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];
        logp = fb->logp[i];

        if (prime > 10923)
            break;

        // invalid root (part of poly->a)
        if (prime == 0)
            continue;

        SIEVE_1X;
        SIEVE_LAST;
        UPDATE_ROOTS;
    }

    logp = 14;
    for (; i<full_fb->fb_14bit_B - 8; i++)
    {	
        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];

        sieve[root1] -= logp;
        sieve[root2] -= logp;
        root1 += prime;
        root2 += prime;

        sieve[root1] -= logp;
        sieve[root2] -= logp;
        root1 += prime;
        root2 += prime;

        if (root1 < 32768) 
        {
            sieve[root1] -= logp;
            root1 += prime;
            if (root2 < 32768)
            {
                sieve[root2] -= logp;
                root2 += prime;
            }
            else
            {
                tmp = root2;
                root2 = root1;
                root1 = tmp;
            }
        }

        UPDATE_ROOTS;
    }

    for (; i<full_fb->fb_15bit_B; i++)
    {	
        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];

        if ((prime > 16384) && ((i & 15) == 0))
            break;

        // invalid root (part of poly->a)
        if (prime == 0)
            continue;

        SIEVE_1X;
        SIEVE_LAST;
        UPDATE_ROOTS;
    }

    logp = 15;
    for (; i<full_fb->fb_15bit_B; i++)
    {	
        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];

        sieve[root1] -= logp;
        sieve[root2] -= logp;
        root1 += prime;
        root2 += prime;

        if (root1 < 32768) 
        {
            sieve[root1] -= logp;
            root1 += prime;
            if (root2 < 32768)
            {
                sieve[root2] -= logp;
                root2 += prime;
            }
            else
            {
                tmp = root2;
                root2 = root1;
                root1 = tmp;
            }
        }

        UPDATE_ROOTS;
    }

    asm_input.logptr = fb->logp;
    asm_input.primeptr = fb->prime;
    asm_input.root1ptr = fb->root1;
    asm_input.root2ptr = fb->root2;
    asm_input.sieve = sieve;
    asm_input.startprime = i;
    asm_input.med_B = med_B;

    SIEVE_GT_BLOCKSIZE_ASM;

#endif


#else

    for (; i<full_fb->fb_15bit_B; i++)
    {	
        uint8 *s2;		

        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];
        logp = fb->logp[i];

        // invalid root (part of poly->a)
        if (prime == 0)
            continue;

        SIEVE_1X;
        SIEVE_LAST;
        UPDATE_ROOTS;
    }

    //if there are primes left bigger than the blocksize, this will take
    //care of them.  if not, it doesn't run at all.
    for (;i<med_B;i++)
    {	
        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];
        logp = fb->logp[i];

        SIEVE_BIG;
        UPDATE_ROOTS;
    }

    

#endif


    

#endif


#ifdef QS_TIMING
	gettimeofday (&qs_timing_stop, NULL);
    SIEVE_STG1 += yafu_difftime (&qs_timing_start, &qs_timing_stop);
	gettimeofday(&qs_timing_start, NULL);
#endif

	return;

}

