/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: smallfact.c 23 2009-07-20 02:59:07Z jasonp_sf $
--------------------------------------------------------------------*/

#include <ms_common.h>

/*--------------------------------------------------------------------*/
uint32 trial_factor(msieve_obj *obj, mp_t *n, mp_t *reduced_n,
			factor_list_t *factor_list) {

	uint32 i, j;
	uint32 factor_found = 0;
	uint32 prime = 0;
	mp_t tmp_n, factor;

	/* Do trial division. If any factors
	   are found they'll get divided out of n */

	mp_copy(n, &tmp_n);
	mp_clear(&factor);
	factor.nwords = 1;

	for (i = 0; i < PRECOMPUTED_NUM_PRIMES; i++) {
		prime += prime_delta[i];

		j = mp_mod_1(&tmp_n, prime);
		if (j == 0) {
			factor.val[0] = prime;
			factor_list_add(obj, factor_list, &factor);
			factor_found = 1;

			/* divide out all instances of 'prime' from tmp_n */

			while (j == 0) {
				mp_divrem_1(&tmp_n, prime, &tmp_n);

				if (mp_is_one(&tmp_n))
					break;

				j = mp_mod_1(&tmp_n, prime);
			}

			/* check if tmp_n has been completely factored,
			   and exit early if so */

			if (tmp_n.nwords <= 2 && 
				((uint64)tmp_n.val[1] << 32 | 
				 tmp_n.val[0]) < ((uint64)prime * prime)) {

				if (!mp_is_one(&tmp_n)) {
					factor_list_add(obj, factor_list, 
								&tmp_n);
					mp_clear(&tmp_n);
					tmp_n.nwords = tmp_n.val[0] = 1;
				}
				break;
			}
		}
	}

	mp_copy(&tmp_n, reduced_n);
	return j;
}

/*--------------------------------------------------------------------*/
/* perform this many iterations between GCDs */
#define BATCH_SIZE 100

uint32 rho(msieve_obj *obj, mp_t *n, mp_t *reduced_n,
		factor_list_t *factor_list) {

	/* simplified implementation of Pollard Rho factoring
     	   (more specifially, Rho with the improvements due to Brent).
    	   The tuning assumes finding 6-8 digit factors */

	uint32 i, j, k, m;
	uint32 num_trials;
	uint32 num_batches; 
	mp_t tmp_n;
	uint32 factor_found = 0;
	uint32 bits = mp_bits(n);

	mp_copy(n, &tmp_n);

	/* Try num_trials pseudorandom functions */

	if (bits < 150)
		num_trials = 2;
	else if (bits < 250)
		num_trials = 3;
	else
		num_trials = 4;

	/* compute the number of iterations to use for
	   each pseudorandom function. This is a multiple
	   of the batch size so that GCD's only have to
	   occur in one place. We do not reduce this
	   bound if a factor is found; factors found
	   are likely to be the smallest in n, and we
	   don't want to hurt our chances of finding larger
	   factors by reducing the amount of work done */

	if (bits < 100)
		num_batches = 70;
	else if (bits < 200)
		num_batches = 200;
	else
		num_batches = 600;

	/* for each function */

	for (i = 0; i < num_trials; i++) {

		/* the function is x[n] = x[n-1]^2 + inc mod n,
		   with inc and x[0] chosen randomly. We can be
		   much more careful about the choice of inc and
		   x[0], and in particular Pari chooses inc from
		   a carefully constructed set of values that
		   ensure no duplicate iterations for x[0] = 2.
		   In practice I doubt that it makes much difference */

		mp_t x, y, acc, t, q, r;
		uint32 inc = get_rand(&obj->seed1, &obj->seed2);
		uint32 limit = 2;

		mp_clear(&x);
		x.nwords = 1;
		x.val[0] = get_rand(&obj->seed1, &obj->seed2);
		mp_copy(&x, &y);
		mp_clear(&acc);
		acc.nwords = acc.val[0] = 1;

		/* for each batch */

		for (j = 2, k = 0; k < num_batches; k++) {

			uint32 parity = 0;

			for (m = 0; m < BATCH_SIZE; m++, j++) {

				/* compute the next function value. Do not
				   reduce mod n after adding inc, since it
				   is so rare for x to exceed n that it's
				   faster to just try another function when
				   x gets corrupted */

				mp_modmul(&x, &x, &tmp_n, &x);
				mp_add_1(&x, inc, &x);

				/* multiply the product of gcd inputs by
				   (y - new_x). Remember the sign of the
				   product to avoid having to add n */

				if (mp_cmp(&y, &x) >= 0) {
					mp_sub(&y, &x, &t);
				}
				else {
					mp_sub(&x, &y, &t);
					parity ^= 1;
				}
				mp_modmul(&acc, &t, &tmp_n, &acc);

				/* possibly update y */

				if (j == limit) {
					mp_copy(&x, &y);
					limit *= 2;
				}
			}

			/* batch is done; compute a GCD, looking 
			   for a nontrivial factor of n */

			if (parity == 1)
				mp_sub(&tmp_n, &acc, &acc);
			mp_gcd(&acc, &tmp_n, &t);

			/* if the GCD is n, all future GCDs will
			   also be n. We could implement backtracking
			   to recover the factor found, but it's easier
			   to just give up and try another pseudorandom
			   function */

			if (mp_cmp(&t, &tmp_n) == 0)
				break;
			if (mp_is_one(&t))
				continue;

			/* factor found; remove all instances of the 
			   factor from n */

			while (1) {
				mp_divrem(&tmp_n, &t, &q, &r);
				if (mp_is_zero(&q) || !mp_is_zero(&r))
					break;
				mp_copy(&q, &tmp_n);
			}

			/* finish up if the new n is prime */

			factor_found = 1;
			if (factor_list_add(obj, factor_list, &t) == 0)
				goto clean_up;

			/* change all the iteration variables
			   to be modulo the new n */

			mp_mod(&x, &tmp_n, &t);
			mp_copy(&t, &x);
			mp_mod(&y, &tmp_n, &t);
			mp_copy(&t, &y);
			mp_mod(&acc, &tmp_n, &t);
			mp_copy(&t, &acc);
		}
	}

clean_up:
	mp_copy(&tmp_n, reduced_n);
	return factor_found;
}
