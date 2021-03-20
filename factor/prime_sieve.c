/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.

$Id: prime_sieve.c 228 2010-02-28 23:37:55Z jasonp_sf $
--------------------------------------------------------------------*/

#include <common.h>
#include "prime_sieve.h"

/* A somewhat basic prime sieve */

#define SIEVE_BLOCK_BITS 65536
#define SIEVE_BLOCK_BYTES (SIEVE_BLOCK_BITS/8)
#define MAX_SIEVE_PRIMES 6542		/* primes less than 2^16 */

static const uint32_t factors[] = {
	3,5,7,11,13,17,19,23,29,31,37,41,43,
	47,53,59,61,67,71,73,79,83,89,97,101,
	103,107,109,113,127,131,137,139,149,
	151,157,163,167,173,179,181,191,193,
	197,199,211,223,227,229,233,239,241,251
};

/*------------------------------------------------------------------*/
static void next_sieve(prime_sieve_t *s) { 

	/* Perform sieving for the next block */

	uint32_t i;
	uint32_t num_aux = s->num_aux;
	prime_aux_t *aux = s->aux;
	uint8_t *sieve = s->sieve;

	memset(sieve, 0, (size_t)SIEVE_BLOCK_BYTES);

	for (i = 0; i < num_aux; i++) {
		uint32_t p = aux[i].p;
		uint32_t r = aux[i].r;
		while (r < SIEVE_BLOCK_BITS) {
			sieve[r / 8] |= 1 << (r % 8);
			r += p;
		}
		aux[i].r = r - SIEVE_BLOCK_BITS;
	}
}

/*------------------------------------------------------------------*/
void init_prime_sieve(prime_sieve_t *s, 
			uint32_t min_prime, uint32_t max_prime) {

	uint32_t i, j;
	uint8_t *sieve;
	prime_aux_t *aux;
	uint32_t size;
	uint32_t block_start;

	sieve = s->sieve = (uint8_t *)xmalloc((size_t)SIEVE_BLOCK_BYTES);
	aux = s->aux = (prime_aux_t *)xmalloc(MAX_SIEVE_PRIMES * 
						sizeof(prime_aux_t));
	memset(sieve, 0, (size_t)SIEVE_BLOCK_BYTES);

	/* sieve with the odd primes < 256 */

	size = MIN(max_prime, 65536/2);
	for (i = 0; i < (sizeof(factors)/sizeof(uint32_t)); i++) {
		uint32_t p = factors[i];
		uint32_t r;

		/* skip the first sieve update (which is a prime),
		   and skip updates for even multiples of p. The
		   sieve is compressed to hold only odd entries,
		   so the doubling of p cancels the halving of the
		   sieve update. Thus, only the starting offset
		   needs modifying to account for the sieve being
		   compressed */

		r = 3 * p / 2;
		while (r < size) {
			sieve[r / 8] |= 1 << (r % 8);
			r += p;
		}
	}

	if (min_prime % 2)
		min_prime--;
	s->curr_block = min_prime / (2 * SIEVE_BLOCK_BITS);
	s->curr_off = (min_prime % (2 * SIEVE_BLOCK_BITS)) / 2;

	if (max_prime <= 65536)
		return;

	/* recover the odd primes less than sqrt(max_prime) */

	size = (uint32_t)(sqrt((double)max_prime) + 1);

	for (i = 1, j = 0; i < SIEVE_BLOCK_BITS; i++) {
		uint32_t p = 2 * i + 1;
		if (p > size)
			break;

		if (!(sieve[i / 8] & (1 << (i % 8)))) {
			aux[j++].p = p;
		}
	}
	s->num_aux = j;

	/* calculate the sieve offsets into the current block.
	   The sieve is compressed so that even multiples of
	   the sieving primes are skipped */

	size = s->num_aux;
	block_start = s->curr_block * 2 * SIEVE_BLOCK_BITS;
	if (block_start == 0) {

		/* if preparing block 0, also skip the first
		   sieve update (since it represents a prime) */

		for (i = 0; i < size; i++)
			aux[i].r = 3 * aux[i].p / 2;
	}
	else {
		for (i = 0; i < size; i++) {
			uint32_t p = aux[i].p;
			uint32_t r = p - (block_start % p);
			if (r % 2 == 0)
				r += p;
			aux[i].r = r / 2;
		}
	}
	next_sieve(s);
}

/*------------------------------------------------------------------*/
void free_prime_sieve(prime_sieve_t *s) { 
	free(s->aux);
	free(s->sieve);
	s->aux = NULL;
	s->sieve = NULL;
}

/*------------------------------------------------------------------*/
uint32_t get_next_prime(prime_sieve_t *s) { 
	
	uint32_t off = s->curr_off;
	uint8_t *sieve = s->sieve;

	if (off == 0 && s->curr_block == 0) {
		s->curr_off = 1;
		return 2;
	}

	while (1) {
		while (off < SIEVE_BLOCK_BITS) {
			if (!(sieve[off / 8] & (1 << (off % 8)))) {
				s->curr_off = off + 1;
				return 2 * (s->curr_block * 
						SIEVE_BLOCK_BITS + off) + 1;
			}
			off++;
		}
		s->curr_block++;
		off = 0;
		next_sieve(s);
	}

	return 0;	/* should never happen */
}

