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

#ifdef TARGET_NOT_MIC
    __m512i vpmul = _mm512_setr_epi32(
        0, 1, 2, 3, 
        4, 5, 6, 7, 
        8, 9, 10, 11, 
        12, 13, 14, 15);

#endif

	med_B = full_fb->med_B;
	
#ifdef QS_TIMING
	gettimeofday(&qs_timing_start, NULL);
#endif

	//initialize the block
	BLOCK_INIT;

    /*
    // beyond 11 bit we don't need to loop... 16 steps takes
    // us past blocksize all at once when prime > 2048
    for (i = start_prime; i < full_fb->fb_13bit_B; i++)
    {	
        __m512i vprime, vroot1, vroot2, vlogp, vidx1, vidx2, vbyte, v16p;
        __mmask16 mask1, mask2;

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

        mask2 = _mm512_cmp_epu32_mask(vidx2, _mm512_set1_epi32(32768), _MM_CMPINT_LT);

        while (mask2 > 0)
        {
            vbyte = _mm512_mask_i32extgather_epi32(vzero, mask2, vidx1, sieve, 
                _MM_UPCONV_EPI32_UINT8, 1, _MM_HINT_NONE);

            vbyte = _mm512_sub_epi32(vbyte, vlogp);

            _mm512_mask_i32extscatter_epi32(sieve, mask2, vidx1, vbyte,
                _MM_DOWNCONV_EPI32_UINT8, 1, _MM_HINT_NONE);

            vbyte = _mm512_mask_i32extgather_epi32(vzero, mask2, vidx2, sieve,
                _MM_UPCONV_EPI32_UINT8, 1, _MM_HINT_NONE);

            vbyte = _mm512_sub_epi32(vbyte, vlogp);

            _mm512_mask_i32extscatter_epi32(sieve, mask2, vidx2, vbyte,
                _MM_DOWNCONV_EPI32_UINT8, 1, _MM_HINT_NONE);

            vidx1 = _mm512_mask_add_epi32(vidx1, mask2, vidx1, v16p);
            vidx2 = _mm512_mask_add_epi32(vidx2, mask2, vidx2, v16p);

            mask2 = _mm512_cmp_epu32_mask(vidx2, _mm512_set1_epi32(32768), _MM_CMPINT_LT);
        }

        mask1 = _mm512_cmp_epu32_mask(vidx1, _mm512_set1_epi32(32768), _MM_CMPINT_LT);

        if (mask1 > 0)
        {
            // some root1 step took an extra iteration so we need to swap
            vbyte = _mm512_mask_i32extgather_epi32(vzero, mask1, vidx1, sieve,
                _MM_UPCONV_EPI32_UINT8, 1, _MM_HINT_NONE);

            vbyte = _mm512_sub_epi32(vbyte, vlogp);

            _mm512_mask_i32extscatter_epi32(sieve, mask1, vidx1, vbyte,
                _MM_DOWNCONV_EPI32_UINT8, 1, _MM_HINT_NONE);

            // find first root >= blocksize
            // write idx vec, num trailing zeros in cmp mask + 1 id into idx vec
            // root update
        }
        else
        {
            // back up idx2 vec one step
            // find first root >= blocksize
            // extract, root update
        }

    }
    */

#if defined(USE_ASM_SMALL_PRIME_SIEVING)

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

    asm_input.logptr = fb->logp;
    asm_input.primeptr = fb->prime;
    asm_input.root1ptr = fb->root1;
    asm_input.root2ptr = fb->root2;
    asm_input.sieve = sieve;
    asm_input.startprime = i;
    asm_input.med_B = med_B;

    SIEVE_GT_BLOCKSIZE_ASM;

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
    SIEVE_STG1 += my_difftime (&qs_timing_start, &qs_timing_stop);
	gettimeofday(&qs_timing_start, NULL);
#endif

	return;

}

