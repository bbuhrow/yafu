/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: mpqs.c 23 2009-07-20 02:59:07Z jasonp_sf $
--------------------------------------------------------------------*/

#include <ms_common.h>
#include "mpqs.h"

sieve_param_t prebuilt_params[] = {
	{64, 100, 40, 1 * 65536},
	{128, 450, 40, 1 * 65536},
	{183, 2000, 40, 1 * 65536},
	{200, 3000, 50, 1 * 65536},
	{212, 5400, 50, 3 * 65536},
	{233, 10000, 100, 3 * 65536},
	{249, 27000, 100, 3 * 65536},
	{266, 50000, 100, 3 * 65536},

	{283, 55000, 80, 3 * 65536},
	{298, 60000, 80, 9 * 65536},
	{315, 80000, 150, 9 * 65536},
	{332,100000, 150, 9 * 65536},
	{348,140000, 150, 9 * 65536},
	{363,210000, 150, 13 * 65536},
	{379,300000, 150, 17 * 65536},
	{395,400000, 150, 21 * 65536},
	{415,500000, 150, 25 * 65536}, /* beyond this point you're crazy */
	{440,700000, 150, 33 * 65536},
	{465,900000, 150, 50 * 65536},
	{490,1100000, 150, 75 * 65536},
	{512,1300000, 150, 100 * 65536},
};

static void get_sieve_params(uint32 bits, 
			     sieve_param_t *params);

static void build_factor_base(mp_t *n, 
				prime_list_t *prime_list,
				fb_t **out_fb,
				uint32 **out_modsqrt,
				uint32 *out_fb_size, 
				uint32 *out_multiplier);

static uint8 choose_multiplier(mp_t *n, 
				prime_list_t *prime_list,
				uint32 fb_size);

/*--------------------------------------------------------------------*/
uint32 factor_mpqs(msieve_obj *obj, mp_t *n, 
			factor_list_t *factor_list) {

	uint32 bits, fb_size;
	prime_list_t prime_list;
	fb_t *factor_base;
	uint32 *modsqrt_array;
	sieve_param_t params;
	relation_t *relation_list = NULL;
	uint32 num_relations = 0;
	la_col_t *cycle_list = NULL;
	uint32 num_cycles = 0;
	poly_t *poly_list = NULL;
	mp_t *poly_a_list = NULL;
	uint64 *bitfield = NULL;
	uint32 multiplier;
	uint32 factor_found = 0;

	/* Calculate the factor base bound */

	bits = mp_bits(n);
	get_sieve_params(bits, &params);
	if (params.fb_size < 100)
		params.fb_size = 100;

	logprintf(obj, "commencing quadratic sieve (%u-digit input)\n",
			strlen(mp_sprintf(n, 10, obj->mp_sprintf_buf)));

	/* Make a prime list, use it to build a factor base */

	fill_prime_list(&prime_list, 2 * params.fb_size + 100, MAX_FB_PRIME);
	build_factor_base(n, &prime_list, &factor_base, 
				&modsqrt_array, &params.fb_size, 
				&multiplier);
	fb_size = params.fb_size;
	logprintf(obj, "using multiplier of %u\n", multiplier);
	free(prime_list.list);
	
	/* Proceed with the algorithm */

	do_sieving(obj, n, &poly_a_list, &poly_list, 
		   factor_base, modsqrt_array, &params, multiplier, 
		   &relation_list, &num_relations,
		   &cycle_list, &num_cycles);

	free(modsqrt_array);
	if (relation_list == NULL || cycle_list == NULL ||
	    cycle_list == NULL || poly_a_list == NULL || poly_list == NULL) {

		free(factor_base);
		return 0;
	}

	solve_linear_system(obj, fb_size, &bitfield, 
			relation_list, cycle_list, &num_cycles);

	if (bitfield != NULL && num_cycles > 0) {
		factor_found = find_factors(obj, n, factor_base, fb_size, 
					cycle_list, num_cycles, 
					relation_list, bitfield, 
					multiplier, poly_a_list, poly_list,
					factor_list);
	}

	free(factor_base);
	free(bitfield);
	free_cycle_list(cycle_list, num_cycles);
	qs_free_relation_list(relation_list, num_relations);
	free(poly_list);
	free(poly_a_list);
	return factor_found;
}

/*--------------------------------------------------------------------*/

/* Implementation of the modified Knuth-Schroeppel multiplier
   algorithm. This borrows ideas from at least four different
   sources, and seems to choose multipliers that are better on
   average than many of the other methods available.
   
   There seem to be many misconceptions about what this algorithm
   is supposed to do. We want to multiply the input number n by a
   small squarefree constant k, chosen so that the factor base 
   for k * n contains as many small primes as possible. Since small primes
   occur more often than big ones, this makes sieve values smaller
   on average and so more likely to be smooth. We quantify this
   by measuring the average contribution of the first NUM_TEST_PRIMES
   primes to sieve values. There are two constraints: first, larger 
   multipliers mean a larger number to factor. Second, we can't spend 
   all day testing multipliers, so the set of multipliers to test should 
   be small. */

#define NUM_TEST_PRIMES 300

#define NUM_MULTIPLIERS (sizeof(mult_list)/sizeof(uint8))

static const uint8 mult_list[] =
	{1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 
	 21, 22, 23, 26, 29, 30, 31, 33, 34, 35, 37, 38,
	 39, 41, 42, 43, 46, 47, 51, 53, 55, 57, 58, 59, 
	 61, 62, 65, 66, 67, 69, 70, 71, 73};

static uint8 choose_multiplier(mp_t *n, prime_list_t *prime_list,
				uint32 fb_size) {
	uint32 i, j;
	uint32 num_primes = MIN(2 * fb_size, NUM_TEST_PRIMES);
	double best_score;
	uint8 best_mult;
	double scores[NUM_MULTIPLIERS];
	uint32 num_multipliers;
	double log2n = mp_log(n);

	num_primes = MIN(num_primes, prime_list->num_primes);

	/* measure the contribution of 2 as a factor of sieve
	   values. The multiplier itself must also be taken into
	   account in the score. scores[i] is the correction that
	   is implicitly applied to the size of sieve values for
	   multiplier i; a negative score makes sieve values 
	   smaller, and so is better */

	for (i = 0; i < NUM_MULTIPLIERS; i++) {
		uint8 curr_mult = mult_list[i];
		uint8 knmod8 = (curr_mult * n->val[0]) % 8;
		double logmult = log((double)curr_mult);

		/* only consider multipliers k such than
		   k*n will not overflow an mp_t */

		if (log2n + logmult > (32 * MAX_MP_WORDS - 2) * M_LN2)
			break;

		scores[i] = 0.5 * logmult;
		switch (knmod8) {
		case 1:
			scores[i] -= 2 * M_LN2;
			break;
		case 5:
			scores[i] -= M_LN2;
			break;
		case 3:
		case 7:
			scores[i] -= 0.5 * M_LN2;
			break;
		/* even multipliers start with a handicap */
		}
	}
	num_multipliers = i;

	/* for the rest of the small factor base primes */

	for (i = 1; i < num_primes; i++) {
		uint32 prime = prime_list->list[i];
		double contrib = log((double)prime) / (prime - 1);
		uint32 modp = mp_mod_1(n, prime);

		for (j = 0; j < num_multipliers; j++) {
			uint8 curr_mult = mult_list[j];
			uint32 knmodp = mp_modmul_1(modp, curr_mult, prime);

			/* if prime i is actually in the factor base
			   for k * n ... */

			if (knmodp == 0 || mp_legendre_1(knmodp, prime) == 1) {

				/* ...add its contribution. A prime p con-
				   tributes log(p) to 1 in p sieve values, plus
				   log(p) to 1 in p^2 sieve values, etc. The
				   average contribution of all multiples of p 
				   to a random sieve value is thus

				   log(p) * (1/p + 1/p^2 + 1/p^3 + ...)
				   = (log(p) / p) * 1 / (1 - (1/p)) 
				   = log(p) / (p-1)

				   This contribution occurs once for each
				   square root used for sieving. There are two
				   roots for each factor base prime, unless
				   the prime divides k*n. In that case there 
				   is only one root */

				if (knmodp == 0)
					scores[j] -= contrib;
				else
					scores[j] -= 2 * contrib;
			}
		}

	}

	/* use the multiplier that generates the best score */

	best_score = 1000.0;
	best_mult = 1;
	for (i = 0; i < num_multipliers; i++) {
		double score = scores[i];
		if (score < best_score) {
			best_score = score;
			best_mult = mult_list[i];
		}
	}
	return best_mult;
}

/*--------------------------------------------------------------------*/
static void build_factor_base(mp_t *n, 
				prime_list_t *prime_list,
				fb_t **out_fb,
				uint32 **out_modsqrt,
				uint32 *out_fb_size, 
				uint32 *out_multiplier) {

	/* Fills in all the factor base information */

	uint32 i, j;
	fb_t *fb_array, *fbptr;
	uint32 fb_size = *out_fb_size;
	uint32 num_primes = prime_list->num_primes;
	uint32 *modsqrt_array;

	*out_multiplier = choose_multiplier(n, prime_list, fb_size);
	mp_mul_1(n, *out_multiplier, n);

	fb_array = (fb_t *)xmalloc(fb_size * sizeof(fb_t));
	modsqrt_array = (uint32 *)xmalloc(fb_size * sizeof(uint32));
	fbptr = fb_array + MIN_FB_OFFSET;
	fbptr->prime = 2;
	fbptr++;

	for (i = 1, j = MIN_FB_OFFSET+1; i < num_primes; i++) {

		uint32 prime = prime_list->list[i];
		uint32 nmodp;

		if (j == fb_size)
			break;

		/* if n is a quadratic residue mod 'prime', then
		   'prime' belongs in the factor base */

		nmodp = mp_mod_1(n, prime);
		if (nmodp == 0 || mp_legendre_1(nmodp, prime) == 1) {

			fbptr->prime = prime;
			fbptr->logprime = (uint8)(log((double)prime) / 
							M_LN2 + .5);

			/* find x such that x^2 mod prime = n mod prime */

			if (nmodp != 0) {
				modsqrt_array[j] = mp_modsqrt_1(nmodp, prime);
			}
			fbptr++;
			j++;
		}
	}
	*out_fb_size = j;
	*out_fb = fb_array;
	*out_modsqrt = modsqrt_array;
}

/*--------------------------------------------------------------------*/
static void get_sieve_params(uint32 bits, sieve_param_t *params) {

	sieve_param_t *low, *high;
	uint32 max_size;
	uint32 i, j, dist;

	/* For small inputs, use the first set of precomputed
	   parameters */

	if (bits < prebuilt_params[0].bits) {
		*params = prebuilt_params[0];
		return;
	}

	/* bracket the input size between two table entries */

	max_size = sizeof(prebuilt_params) / sizeof(sieve_param_t);
	if (bits >= prebuilt_params[max_size - 1].bits) {
		*params = prebuilt_params[max_size - 1];
		return;
	}

	/* if the input is too large, just use the last table entry.
	   This means that the choice of parameters is increasingly
	   inappropriate as the input becomes larger, but there's no
	   guidance on what to do in this case anyway. */

	for (i = 0; i < max_size - 1; i++) {
		if (bits < prebuilt_params[i+1].bits)
			break;
	}

	/* Otherwise the parameters to use are a weighted average 
	   of the two table entries the input falls between */

	low = &prebuilt_params[i];
	high = &prebuilt_params[i+1];
	dist = high->bits - low->bits;
	i = bits - low->bits;
	j = high->bits - bits;

	params->bits = bits;
	params->fb_size = (uint32)(
			 ((double)low->fb_size * j +
			  (double)high->fb_size * i) / dist + 0.5);
	params->large_mult = (uint32)(
			 ((double)low->large_mult * j +
			  (double)high->large_mult * i) / dist + 0.5);
	params->sieve_size = (uint32)(
			 ((double)low->sieve_size * j +
			  (double)high->sieve_size * i) / dist + 0.5);
}
