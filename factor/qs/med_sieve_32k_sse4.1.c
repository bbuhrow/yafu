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
	uint8_t *sieve;					//0
	uint16_t *primeptr;				//8
	uint16_t *root1ptr;				//16
	uint16_t *root2ptr;				//24
	uint16_t *logptr;					//32
	uint32_t startprime;				//40
	uint32_t med_B;					//44
} helperstruct_t;

void med_sieveblock_32k_sse41(uint8_t *sieve, sieve_fb_compressed *fb, fb_list *full_fb, 
		uint32_t start_prime, uint8_t s_init)
{
	uint32_t i;
	uint32_t med_B;
	
	uint32_t prime, root1, root2, tmp, stop;
	uint8_t logp;

	helperstruct_t asm_input;

	med_B = full_fb->med_B;
	
#ifdef QS_TIMING
	gettimeofday(&qs_timing_start, NULL);
#endif

	//initialize the block
	BLOCK_INIT;

	// sieve primes less than 2^13 using optimized loops: it becomes
	// inefficient to do fully unrolled sse2 loops as the number of
	// steps through the block increases.  While I didn't specifically
	// test whether 2^13 is the best point to start using loops, it seemed
	// good enough.  any gains if it is not optmial will probably be minimal.
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
		uint8_t *s2;		

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

	// the small prime sieve stops just before prime exceeds 2^13
	// the next sse2 sieve assumes primes exceed 2^13.  since
	// some of the primes in the next set of 8 primes could be less
	// than the cutoff and some are greater than, we have to do this
	// small set of crossover primes manually, one at a time.
	for (; i<full_fb->fb_13bit_B; i++)
	{	
		uint8_t *s2;		

		prime = fb->prime[i];
		root1 = fb->root1[i];
		root2 = fb->root2[i];
		logp = fb->logp[i];

		//// special exit condition: when prime > 8192 and i % 8 is 0;
		//if ((prime > 8192) && ((i&7) == 0))
		//	break;

		// invalid root (part of poly->a)
		if (prime == 0) 
			continue;

		SIEVE_2X;
		SIEVE_1X;
		SIEVE_LAST;
		UPDATE_ROOTS;
	}

	// sieve primes 8 at a time, where 8192 < p < blocksize/3
	_INIT_SSE2_SMALL_PRIME_SIEVE;
	_SSE2_SMALL_PRIME_SIEVE_32k_DIV3;

	// the sse2 sieve stops just before prime exceeds blocksize/3
	// the next sse2 sieve assumes primes exceed blocksize/3.  since
	// some of the primes in the next set of 8 primes could be less
	// than the cutoff and some are greater than, we have to do this
	// small set of crossover primes manually, one at a time.
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

	// sieve primes 8 at a time, where blocksize/3 < p < 2^14
	_INIT_SSE2_SMALL_PRIME_SIEVE;
	_SSE2_SMALL_PRIME_SIEVE_14b;

	// do this small set of crossover primes manually, one at a time,
	// this time for the 14 bit crossover.
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

	// sieve primes 8 at a time, 2^14 < p < 2^15
	_INIT_SSE2_SMALL_PRIME_SIEVE;
	_SSE2_SMALL_PRIME_SIEVE_15b;

	// do this small set of crossover primes manually, one at a time,
	// this time for the 15 bit crossover.
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

	// sieve primes 8 at a time, 2^15 < p < med_B
	_INIT_SSE2_SMALL_PRIME_SIEVE;
	_SSE41_SMALL_PRIME_SIEVE;


#ifdef QS_TIMING
	gettimeofday (&qs_timing_stop, NULL);
    SIEVE_STG1 += ytools_difftime (&qs_timing_start, &qs_timing_stop);
	gettimeofday(&qs_timing_start, NULL);
#endif

	return;

}

#endif // USE_SSE41

