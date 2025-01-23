/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: root_score.c 838 2012-12-17 03:41:49Z jasonp_sf $
--------------------------------------------------------------------*/

#include "poly.h"

#define ROOT_BUF_SIZE 200

/*------------------------------------------------------------------*/
static double get_root_freq(mpz_poly_t *poly, 
			uint32 initial_root, uint32 p) {
	
	/* find the number of roots of poly mod p^k for
	   all k where p^k <= PRIME_BOUND, then add up
	   the resulting contribution to random values of
	   poly */

	uint32 i, j;
	uint32 coeffs[MAX_POLY_DEGREE + 1];
	uint32 roots[ROOT_BUF_SIZE];
	uint32 new_roots[ROOT_BUF_SIZE];
	uint32 max_power, power, degree;
	uint32 num_roots, num_new_roots;
	double contrib;

	/* find the largest power of p that does not exceed 2^16 */

	max_power = 1;
	power = p;
	while (1) {
		uint32 next_power = power * p;
		if (next_power > 65000)
			break;

		power = next_power;
		max_power++;
	}

	/* compute poly mod power */

	degree = poly->degree;
	for (i = 0; i <= degree; i++)
		coeffs[i] = mpz_fdiv_ui(poly->coeff[i], power);
	while (degree > 0 && coeffs[degree] == 0)
		degree--;

	/* there is only one root mod p */

	num_roots = 1;
	roots[0] = initial_root;
	power = p;
	contrib = 1.0 / power;

	for (i = 2; i <= max_power; i++) {

		uint32 next_power = power * p;
		uint32 x;

		/* search for roots mod p^i using the roots
		   mod p^(i-1) */

		for (j = num_new_roots = 0; j < num_roots; j++) {

			/* all roots mod p^i must have a known value
			   mod p^(i-1). One root mod p^(i-1) can be a
			   multiple root, and can give rise to several 
			   roots mod p^i, or to no roots at all */

			for (x = roots[j]; x < next_power; x += power) {

				uint32 k;
				uint32 polyval = coeffs[degree] % next_power;
				for (k = degree; k; k--) {
					polyval = (polyval * x + 
						coeffs[k-1]) % next_power;
				}
				if (polyval == 0 && 
				    num_new_roots < ROOT_BUF_SIZE) {
					new_roots[num_new_roots++] = x;
				}
			}
		}

		if (num_new_roots == 0)
			break;

		/* add in the contribution due to roots mod p^i 
		   and make the new roots into the old roots */

		contrib += (double)num_new_roots / next_power;
		memcpy(roots, new_roots, num_new_roots * sizeof(uint32));
		power = next_power;
		num_roots = num_new_roots;
	}

	/* for ordinary roots, contrib would be 1/(p-1)
	   asymptotically */

	return contrib * p / (p + 1);
}

/*------------------------------------------------------------------*/
#define SMALL_PRIME_BOUND 100

uint32 analyze_poly_roots(mpz_poly_t *poly, uint32 prime_bound,
				double *result) {

	/* analyze a polynomial for root properties (i.e.
	   compute Murphy's 'alpha' value) */

	uint32 i, j;
	double root_score;
	uint32 prime;
	mpz_poly_t rev_poly;

	/* all linear polynomials have the same root properties;
	   note that we still must treat them as generating 
	   homogeneous polynomial values and *not* as random
	   numbers. A linear NFS polynomial will always have one
	   root modulo each prime p, leading to a fixed but
	   nonzero alpha value */

	if (poly->degree == 1) {
		*result = 0.569959993064325;
		return 0;
	}

	/* handling projective roots requires a reversed version
	   of poly (we find roots of poly(1/x) instead of poly(x)) */

	mpz_poly_init(&rev_poly);
	j = poly->degree;
	for (i = 0; i <= j; i++) {
		mpz_set(rev_poly.coeff[i], poly->coeff[j-i]);
	}
	while (j && mpz_cmp_ui(rev_poly.coeff[j], 0) == 0) {
		j--;
	}
	rev_poly.degree = j;

	/* roots of (poly mod p) can be either ordinary or special.
	   An ordinary root mod p corresponds to exactly one root
	   mod p^k for all k > 1. Murphy shows that for each ordinary 
	   root of (poly mod p) the contribution to the size of a
	   random sieve value due to all powers of p is log(p) * 
	   p/(p^2-1). 
	   
	   We compare this to the contribution of powers of p to 
	   the size of a random number. 1 in p random numbers are 
	   divisible by p, 1 in p^2 random numbers are divisible by
	   p^2, etc. Adding up all of these contributions and summing
	   the resulting geometric series yields an average contribution
	   of log(p) / (p-1). 
	   
	   Hence the root score for poly is the sum of the difference 
	   between the contribution of powers of p in sieve values 
	   compared to the contribution to a random number. When there 
	   are k ordinary roots mod p, root_score is the sum of 

	               log(p) * (1/(p-1) - k*p/(p^2-1))

	   for the first few p, meaning that dividing the small primes 
	   out of sieve values makes their size exp(root_score) times 
	   smaller than dividing the small primes out of a random number.

	   When the number of roots mod p^k changes with k, Murphy's 
	   formula for the contribution is inaccurate (usually it's too
	   small), so we have to take special measures to compute the 
	   contribution of the roots accurately. Roots become special 
	   if p divides the discriminant of poly, and basically the only 
	   effective way to compute the size of the contribution from 
	   special roots is to find the number of roots mod p^k for the
	   first few k and then explicitly sum the first few terms in 
	   the (formerly geometric) series. Many thanks to Kleinjung 
	   and Franke whose pol5 code handles special roots efficiently */

	root_score = 0;
	for (i = prime = 0; i < PRECOMPUTED_NUM_PRIMES; i++) {

		uint32 roots[MAX_POLY_DEGREE];
		uint32 mult[MAX_POLY_DEGREE];
		uint32 high_coeff;
		uint32 num_roots, num_ordinary_roots;
		double contrib;

		/* test if next prime will exceed the bound */

		if (prime + prime_delta[i] > SMALL_PRIME_BOUND || 
		    prime + prime_delta[i] > prime_bound)
			break;

		/* find the roots of poly mod prime and determine if
		   any are multiple roots (this is the telltale sign
		   that the root is special) */

		prime += prime_delta[i];
		num_roots = poly_get_zeros_and_mult(roots, mult, poly, 
					prime, &high_coeff);

		contrib = 1.0 / (prime - 1);

		/* test each root for specialness */

		for (j = num_ordinary_roots = 0; j < num_roots; j++) {
			if (mult[j] == 0)
				num_ordinary_roots++;
			else
				contrib -= get_root_freq(poly, roots[j], prime);
		}
		
		/* projective roots are always handled specially */

		if (high_coeff == 0)
			contrib -= get_root_freq(&rev_poly, 0, prime);

		/* add in the root score contribution for prime */

		root_score += log((double)prime) * (contrib - 
						(double)num_ordinary_roots * 
						prime / (prime * prime - 1));
	}

	/* for the larger primes, assume that all roots are
	   ordinary, so that only the number of roots modulo
	   prime i matter */

	for (; i < PRECOMPUTED_NUM_PRIMES; i++) {
		uint32 num_roots;
		uint32 dummy_roots[MAX_POLY_DEGREE];
		uint32 high_coeff;

		prime += prime_delta[i];
		if (prime > prime_bound)
			break;

		num_roots = poly_get_zeros(dummy_roots, poly, 
					prime, &high_coeff, 1);
		if (high_coeff == 0)
			num_roots++;

		root_score += (1.0 - (double)num_roots * prime / 
					(prime + 1)) * log((double)prime) / 
					(prime - 1);
	}

	*result = root_score;
	mpz_poly_free(&rev_poly);
	return 0;
}

/*------------------------------------------------------------------*/
uint32 analyze_poly_roots_projective(mpz_poly_t *poly, 
				uint32 prime_bound,
				double *result) {

	/* compute the contribution of projective roots
	   to random vaues of poly() */

	uint32 i, j;
	double root_score;
	uint32 prime;
	mpz_poly_t rev_poly;

	/* handling projective roots requires a reversed version
	   of poly (we find roots of poly(1/x) instead of poly(x)) */

	mpz_poly_init(&rev_poly);
	j = poly->degree;
	for (i = 0; i <= j; i++) {
		mpz_set(rev_poly.coeff[i], poly->coeff[j-i]);
	}
	while (j && mpz_cmp_ui(rev_poly.coeff[j], 0) == 0) {
		j--;
	}
	rev_poly.degree = j;

	root_score = 0;
	for (i = prime = 0; i < PRECOMPUTED_NUM_PRIMES; i++) {

		prime += prime_delta[i];
		if (prime > prime_bound)
			break;

		/* only projective roots contribute */

		if (mpz_fdiv_ui(poly->coeff[poly->degree], prime) == 0) {
			root_score -= log((double)prime) *
					get_root_freq(&rev_poly, 0, prime);
		}
	}

	*result = root_score;
	mpz_poly_free(&rev_poly);
	return 0;
}
