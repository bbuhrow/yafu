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
#include "sieve_macros_32k_sse4.1.h"

// protect sse41 code under MSVC builds.  USE_SSE41 should be manually
// enabled at the top of qs.h for MSVC builds on supported hardware
#ifdef USE_SSE41

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

void med_sieveblock_32k_sse41(uint8 *sieve, sieve_fb_compressed *fb, fb_list *full_fb, 
		uint32 start_prime, uint8 s_init)
{
	uint32 i,j;
	uint32 med_B;
	
	uint32 prime, root1, root2, tmp, stop;
	uint8 logp;

	helperstruct_t asm_input;

	med_B = full_fb->med_B;
	
#ifdef QS_TIMING
	gettimeofday(&qs_timing_start, NULL);
#endif

	//initialize the block
	BLOCK_INIT;

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
	stop = MIN(8, full_fb->fb_32k_div3 - i);

	for (j=0;j<stop;i++, j++)
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
	stop = MIN(8, full_fb->fb_14bit_B - i);

	for (j=0; j<stop; i++, j++)
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
	for (i = full_fb->fb_15bit_B - 8; i<med_B; i++)
	{	
		prime = fb->prime[i];
		root1 = fb->root1[i];
		root2 = fb->root2[i];
		logp = fb->logp[i];

		if ((prime > 32768) && ((i&7) == 0))
			break;

		// invalid root (part of poly->a)
		if (prime == 0) 
			continue;

		SIEVE_1X;
		SIEVE_LAST;

		UPDATE_ROOTS;
	}

	_INIT_SSE2_SMALL_PRIME_SIEVE;
	_SSE41_SMALL_PRIME_SIEVE;


#ifdef QS_TIMING
	gettimeofday (&qs_timing_stop, NULL);
	qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

	SIEVE_STG1 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
	free(qs_timing_diff);

	gettimeofday(&qs_timing_start, NULL);
#endif

	return;

}

#endif // USE_SSE41

