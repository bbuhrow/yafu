/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.

$Id: batch_factor.c 638 2011-09-11 15:31:19Z jasonp_sf $
--------------------------------------------------------------------*/

#include <ms_batch_factor.h>
#include "savefile.h"

/*------------------------------------------------------------------*/
#define BREAKOVER_WORDS 50

static void multiply_primes(uint32 first, uint32 last, 
				prime_sieve_t *sieve, mpz_t prod) {

	/* recursive routine to multiply all the elements of a
	   list of consecutive primes. The current invocation
	   multiplies elements first to last of the list, inclusive.
	   The primes are created on demand */

	mpz_t half_prod;
	uint32 mid = (last + first) / 2;

	/* base case; accumulate a few primes */

	if (last < first + BREAKOVER_WORDS) {
		mpz_set_ui(prod, (unsigned long)get_next_prime(sieve));
		while (++first <= last) {
			mpz_mul_ui(prod, prod, 
				 (unsigned long)get_next_prime(sieve));
		}
		return;
	}

	/* recursively handle the left and right halves of the list */

	mpz_init(half_prod);
	multiply_primes(first, mid, sieve, prod);
	multiply_primes(mid + 1, last, sieve, half_prod);

	/* multiply them together. We can take advantage of 
	   fast multiplication since in general the two half-
	   products are about the same size*/

	mpz_mul(prod, prod, half_prod);
	mpz_clear(half_prod);
}

/*------------------------------------------------------------------*/
static void relation_to_gmp(relation_batch_t *rb,
				uint32 index, mpz_t out) {

	/* multiply together the unfactored parts of an
	   NFS relation, then convert to an mpz_t */

	uint32 i, nwords;
	cofactor_t *c = rb->relations + index;
	uint32 *f = rb->factors + c->factor_list_word + 
			c->num_factors_r + c->num_factors_a;
	mp_t num1, num2, num3;

	if (c->lp_r_num_words > 1 && c->lp_a_num_words > 1) {

		/* rational and algebraic parts need to be
		   multiplied together first */

		mp_clear(&num1);
		mp_clear(&num2);
		num1.nwords = c->lp_r_num_words;
		for (i = 0; i < c->lp_r_num_words; i++)
			num1.val[i] = f[i];
		num2.nwords = c->lp_a_num_words;
		f += c->lp_r_num_words;
		for (i = 0; i < c->lp_a_num_words; i++)
			num2.val[i] = f[i];
		mp_mul(&num1, &num2, &num3);
	}
	else {
		/* no multiply necesary */

		mp_clear(&num3);
		nwords = c->lp_r_num_words;

		if (c->lp_a_num_words > 1) {
			nwords = c->lp_a_num_words;
			f += c->lp_r_num_words;
		}

		num3.nwords = nwords;
		for (i = 0; i < nwords; i++)
			num3.val[i] = f[i];
	}
	mp2gmp(&num3, out);
}

/*------------------------------------------------------------------*/
static void multiply_relations(uint32 first, uint32 last, 
				relation_batch_t *rb,
				mpz_t prod) {
	mpz_t half_prod;

	/* recursive routine to multiply together a list of
	   relations. The current invocation multiplies relations
	   first to last of the list, inclusive */

	/* base case of recursion */

	if (first == last) {
		relation_to_gmp(rb, first, prod);
		return;
	}

	/* recurse on the left and right halves */

	mpz_init(half_prod);
	if (last == first + 1) {
		relation_to_gmp(rb, first, half_prod);
		relation_to_gmp(rb, last, prod);
	}
	else {
		uint32 mid = (last + first) / 2;
		multiply_relations(first, mid, rb, prod);
		multiply_relations(mid + 1, last, rb, half_prod);
	}

	/* multiply the halves */

	mpz_mul(prod, prod, half_prod);
	mpz_clear(half_prod);
}

/*------------------------------------------------------------------*/
static mp_t two = {1, {2}};

/* the following recursion base-case is specialized for relations
   containing <= 3 rational and/or algebraic large primes. The rest
   of the batch factoring handles an arbitrary number of large primes,
   so only this routine needs to change when more advanced factoring
   code becomes available. */

static void check_relation(relation_batch_t *rb,
			uint32 index,
			mp_t *prime_product) {

	uint32 i;
	cofactor_t *c = rb->relations + index;
	uint32 *f = rb->factors + c->factor_list_word;
	uint32 *lp1 = f + c->num_factors_r + c->num_factors_a;
	uint32 *lp2 = lp1 + c->lp_r_num_words;
	mp_t n, f1r, f2r, f1a, f2a, t0, t1;
	uint32 lp_r[MAX_LARGE_PRIMES];
	uint32 lp_a[MAX_LARGE_PRIMES];
	uint32 num_r, num_a;
	mp_t *_small, *large;

	/* first compute gcd(prime_product, rational large cofactor).
	   The rational part will split into a cofactor with all 
	   factors <= the largest prime in rb->prime_product (stored
	   in f1r) and a cofactor with all factors larger (stored in f2r) */

	if (c->lp_r_num_words) {
		mp_clear(&n);
		n.nwords = c->lp_r_num_words;
		for (i = 0; i < c->lp_r_num_words; i++)
			n.val[i] = lp1[i];

		if (n.nwords == 1) {
			mp_copy(&n, &f1r);
			mp_clear(&f2r);
			f2r.nwords = f2r.val[0] = 1;
		}
		else {
			mp_gcd(prime_product, &n, &f1r);
			if (mp_is_one(&f1r))
				mp_copy(&n, &f2r);
			else
				mp_div(&n, &f1r, &f2r);
		}

		/* give up on this relation if
		     - f1r has a single factor, and that factor
		       exceeds the rational large prime cutoff
		     - f1r is more than 3 words long */

		if (f1r.nwords > 3 ||
		    (f1r.nwords == 1 && f1r.val[0] > rb->lp_cutoff_r))
			return;

		/* give up on this relation if
		     - f2r has a single factor, and that factor
		       exceeds the rational large prime cutoff
		     - f2r has two words and exceeds the square of the
		       cutoff (meaning at least one of the two factors in
		       f2r would exceed the large prime cutoff)
		     - f2r has two words and is smaller than the square
		       of the largest prime in rb->prime_product, meaning
		       that f2r is prime and way too large
		     - f2r has 3+ words. We don't know anything about
		       f2r in this case, except that all factors exceed the
		       largest prime in rb->prime_product, and it would be
		       much too expensive to find out anything more */

		if (f2r.nwords >= 3 ||
		    (f2r.nwords == 2 && 
		     		(mp_cmp(&f2r, &rb->max_prime2) <= 0 ||
		     		 mp_cmp(&f2r, &rb->lp_cutoff_r2) > 0)) ||
		    (f2r.nwords == 1 && f2r.val[0] > rb->lp_cutoff_r))
			return;
	}
	else {
		mp_clear(&f1r);
		mp_clear(&f2r);
	}

	/* repeat with the algebraic unfactored part, if any */

	if (c->lp_a_num_words) {
		mp_clear(&n);
		n.nwords = c->lp_a_num_words;
		for (i = 0; i < c->lp_a_num_words; i++)
			n.val[i] = lp2[i];

		if (n.nwords == 1) {
			mp_copy(&n, &f1a);
			mp_clear(&f2a);
			f2a.nwords = f2a.val[0] = 1;
		}
		else {
			mp_gcd(prime_product, &n, &f1a);
			if (mp_is_one(&f1a))
				mp_copy(&n, &f2a);
			else
				mp_div(&n, &f1a, &f2a);
		}

		if (f1a.nwords > 3 ||
		    (f1a.nwords == 1 && f1a.val[0] > rb->lp_cutoff_a))
			return;

		if (f2a.nwords >= 3 ||
		    (f2a.nwords == 2 && 
		     		(mp_cmp(&f2a, &rb->max_prime2) <= 0 ||
		     		 mp_cmp(&f2a, &rb->lp_cutoff_a2) > 0)) ||
		    (f2a.nwords == 1 && f2a.val[0] > rb->lp_cutoff_a))
			return;
	}
	else {
		mp_clear(&f1a);
		mp_clear(&f2a);
	}

	/* the relation isn't obviously bad; do more work
	   trying to factor everything. Note that when relations 
	   are expected to have three large primes then ~98% of 
	   relations do not make it to this point
	
	   Begin by performing compositeness tests on f2r and f2a,
	   which are necessary if they are two words in size and
	   f1r or f1a is not one (the latter being true means
	   that the sieving already performed the compositeness test) */

	for (i = num_r = num_a = 0; i < MAX_LARGE_PRIMES; i++)
		lp_r[i] = lp_a[i] = 1;

	if (f2r.nwords == 2 && !mp_is_one(&f1r)) {
		mp_t exponent, res;
		mp_sub_1(&f2r, 1, &exponent);
		mp_expo(&two, &exponent, &f2r, &res);
		if (mp_is_one(&res))
			return;
	}
	if (f2a.nwords == 2 && !mp_is_one(&f1a)) {
		mp_t exponent, res;
		mp_sub_1(&f2a, 1, &exponent);
		mp_expo(&two, &exponent, &f2a, &res);
		if (mp_is_one(&res))
			return;
	}

	/* now perform all the factorizations that
	   require SQUFOF, since it is much faster than the
	   QS code. We have to check all of f[1|2][r|a] but
	   for relations with three large primes then at most
	   two of the four choices need factoring */

	if (f1r.nwords == 1) {
		if (f1r.val[0] > 1)
			lp_r[num_r++] = f1r.val[0];
	}
	else if (f1r.nwords == 2) {
		i = squfof(&f1r);
		if (i <= 1 || i > rb->lp_cutoff_r)
			return;
		lp_r[num_r++] = i;
		mp_divrem_1(&f1r, i, &f1r);
		if (f1r.nwords > 1 || f1r.val[0] > rb->lp_cutoff_r)
			return;
		lp_r[num_r++] = f1r.val[0];
	}

	if (f2r.nwords == 1) {
		if (f2r.val[0] > 1)
			lp_r[num_r++] = f2r.val[0];
	}
	else if (f2r.nwords == 2) {
		i = squfof(&f2r);
		if (i <= 1 || i > rb->lp_cutoff_r)
			return;
		lp_r[num_r++] = i;
		mp_divrem_1(&f2r, i, &f2r);
		if (f2r.nwords > 1 || f2r.val[0] > rb->lp_cutoff_r)
			return;
		lp_r[num_r++] = f2r.val[0];
	}

	if (f1a.nwords == 1) {
		if (f1a.val[0] > 1)
			lp_a[num_a++] = f1a.val[0];
	}
	else if (f1a.nwords == 2) {
		i = squfof(&f1a);
		if (i <= 1 || i > rb->lp_cutoff_a)
			return;
		lp_a[num_a++] = i;
		mp_divrem_1(&f1a, i, &f1a);
		if (f1a.nwords > 1 || f1a.val[0] > rb->lp_cutoff_a)
			return;
		lp_a[num_a++] = f1a.val[0];
	}

	if (f2a.nwords == 1) {
		if (f2a.val[0] > 1)
			lp_a[num_a++] = f2a.val[0];
	}
	else if (f2a.nwords == 2) {
		i = squfof(&f2a);
		if (i <= 1 || i > rb->lp_cutoff_a)
			return;
		lp_a[num_a++] = i;
		mp_divrem_1(&f2a, i, &f2a);
		if (f2a.nwords > 1 || f2a.val[0] > rb->lp_cutoff_a)
			return;
		lp_a[num_a++] = f2a.val[0];
	}

	/* only use expensive factoring methods when we know 
	   f1r and/or f1a splits into three large primes and 
	   we know all three primes are smaller than the 
	   largest prime in rb->prime_product. When the latter 
	   is a good deal smaller than the large prime cutoff
	   this happens extremely rarely */

	if (f1r.nwords == 3) {
		if (tinyqs(&f1r, &t0, &t1) == 0)
			return;

		_small = &t0;
		large = &t1;
		if (mp_cmp(_small, large) > 0) {
			_small = &t1;
			large = &t0;
		}

		if (_small->nwords > 1 || _small->val[0] > rb->lp_cutoff_r)
			return;
		lp_r[num_r++] = _small->val[0];
		i = squfof(large);
		if (i <= 1 || i > rb->lp_cutoff_r)
			return;
		lp_r[num_r++] = i;
		mp_divrem_1(large, i, large);
		if (large->nwords > 1 || large->val[0] > rb->lp_cutoff_r)
			return;
		lp_r[num_r++] = large->val[0];
	}

	if (f1a.nwords == 3) {
		if (tinyqs(&f1a, &t0, &t1) == 0)
			return;

		_small = &t0;
		large = &t1;
		if (mp_cmp(_small, large) > 0) {
			_small = &t1;
			large = &t0;
		}

		if (_small->nwords > 1 || _small->val[0] > rb->lp_cutoff_a)
			return;
		lp_a[num_a++] = _small->val[0];
		i = squfof(large);
		if (i <= 1 || i > rb->lp_cutoff_a)
			return;
		lp_a[num_a++] = i;
		mp_divrem_1(large, i, large);
		if (large->nwords > 1 || large->val[0] > rb->lp_cutoff_a)
			return;
		lp_a[num_a++] = large->val[0];
	}

	/* yay! Another relation found */

	rb->num_success++;
	rb->print_relation(rb->savefile, c->a, c->b,
			f, c->num_factors_r, lp_r,
			f + c->num_factors_r, c->num_factors_a, lp_a);
}

/*------------------------------------------------------------------*/
static void compute_remainder_tree(uint32 first, uint32 last,
				relation_batch_t *rb,
				mpz_t numerator) {

	/* recursively compute numerator % (each relation in 
	   rb->relations) */

	uint32 mid = (first + last) / 2;
	mpz_t relation_prod, remainder;

	/* recursion base case: numerator already fits in
	   an mp_t, so manually compute the remainder and
	   postprocess each relation */

	if (mpz_size(numerator) * GMP_LIMB_BITS/32 <= MAX_MP_WORDS) {
		if (mpz_sgn(numerator) > 0) {
			mp_t num;
			gmp2mp(numerator, &num);
			while (first <= last)
				check_relation(rb, first++, &num);
		}
		return;
	}

	/* multiply together the unfactored parts of all the
	   relations from first to last */

	mpz_init(relation_prod);
	mpz_init(remainder);
	multiply_relations(first, last, rb, relation_prod);

	/* use the remainder to deal with the left and right
	   halves of the relation list */

	if (mpz_cmp(numerator, relation_prod) < 0) {
		mpz_clear(relation_prod);
		compute_remainder_tree(first, mid, rb, numerator);
		compute_remainder_tree(mid + 1, last, rb, numerator);
	}
	else {
		mpz_tdiv_r(remainder, numerator, relation_prod);
		mpz_clear(relation_prod);
		compute_remainder_tree(first, mid, rb, remainder);
		compute_remainder_tree(mid + 1, last, rb, remainder);
	}
	mpz_clear(remainder);
}

/*------------------------------------------------------------------*/
void relation_batch_init(msieve_obj *obj, relation_batch_t *rb,
			uint32 min_prime, uint32 max_prime,
			uint32 lp_cutoff_r, uint32 lp_cutoff_a, 
			savefile_t *savefile,
			print_relation_t print_relation) {

	prime_sieve_t sieve;
	uint32 num_primes, p;

	/* count the number of primes to multiply. Knowing this
	   in advance makes the recursion a lot easier, at the cost
	   of a small penalty in runtime */

	init_prime_sieve(&sieve, min_prime + 1, max_prime);
	p = min_prime;
	num_primes = 0;
	while (p < max_prime) {
		p = get_next_prime(&sieve);
		num_primes++;
	}
	free_prime_sieve(&sieve);

	/* compute the product of primes */

	logprintf(obj, "multiplying %u primes from %u to %u\n",
			num_primes, min_prime, max_prime);

	init_prime_sieve(&sieve, min_prime, max_prime);
	mpz_init(rb->prime_product);
	multiply_primes(0, num_primes - 2, &sieve, rb->prime_product);
	logprintf(obj, "multiply complete, product has %u bits\n", 
				(uint32)mpz_sizeinbase(rb->prime_product, 2));
					
	rb->savefile = savefile;
	rb->print_relation = print_relation;

	/* compute the cutoffs used by the recursion base-case. Large
	   primes have a maximum size specified as input arguments, 
	   but numbers that can be passed to the SQUFOF routine are
	   limited to size 2^62 */

	rb->lp_cutoff_r = lp_cutoff_r;
	lp_cutoff_r = MIN(lp_cutoff_r, 0x7fffffff);
	mp_clear(&rb->lp_cutoff_r2);
	rb->lp_cutoff_r2.nwords = 1;
	rb->lp_cutoff_r2.val[0] = lp_cutoff_r;
	mp_mul_1(&rb->lp_cutoff_r2, lp_cutoff_r, &rb->lp_cutoff_r2);

	rb->lp_cutoff_a = lp_cutoff_a;
	lp_cutoff_a = MIN(lp_cutoff_a, 0x7fffffff);
	mp_clear(&rb->lp_cutoff_a2);
	rb->lp_cutoff_a2.nwords = 1;
	rb->lp_cutoff_a2.val[0] = lp_cutoff_a;
	mp_mul_1(&rb->lp_cutoff_a2, lp_cutoff_a, &rb->lp_cutoff_a2);

	mp_clear(&rb->max_prime2);
	rb->max_prime2.nwords = 1;
	rb->max_prime2.val[0] = max_prime;
	mp_mul_1(&rb->max_prime2, max_prime, &rb->max_prime2);

	/* allocate lists for relations and their factors */

	rb->target_relations = 500000;
	rb->num_relations = 0;
	rb->num_relations_alloc = 1000;
	rb->relations = (cofactor_t *)xmalloc(rb->num_relations_alloc *
						sizeof(cofactor_t));

	rb->num_factors = 0;
	rb->num_factors_alloc = 10000;
	rb->factors = (uint32 *)xmalloc(rb->num_factors_alloc *
						sizeof(uint32));
}

/*------------------------------------------------------------------*/
void relation_batch_free(relation_batch_t *rb) {

	mpz_clear(rb->prime_product);
	free(rb->relations);
	free(rb->factors);
}

/*------------------------------------------------------------------*/
void relation_batch_add(int64 a, uint32 b, 
			uint32 *factors_r, uint32 num_factors_r, 
			mpz_t unfactored_r_in,
			uint32 *factors_a, uint32 num_factors_a, 
			mpz_t unfactored_a_in,
			relation_batch_t *rb) {

	uint32 i;
	uint32 *f;
	cofactor_t *c;
	mp_t tmp1, tmp2;
	mp_t *unfactored_r = &tmp1;
	mp_t *unfactored_a = &tmp2;

	/* add one relation to the batch */

	gmp2mp(unfactored_r_in, unfactored_r);
	gmp2mp(unfactored_a_in, unfactored_a);

	if (rb->num_relations == rb->num_relations_alloc) {
		rb->num_relations_alloc *= 2;
		rb->relations = (cofactor_t *)xrealloc(rb->relations,
					rb->num_relations_alloc *
					sizeof(cofactor_t));
	}
	c = rb->relations + rb->num_relations++;
	c->a = a;
	c->b = b;
	c->num_factors_r = num_factors_r;
	c->num_factors_a = num_factors_a;
	c->lp_r_num_words = 0;
	c->lp_a_num_words = 0;
	c->factor_list_word = rb->num_factors;

	/* add its small factors */

	if (rb->num_factors + num_factors_r + num_factors_a +
			2 * MAX_LARGE_PRIMES >= rb->num_factors_alloc) {
		rb->num_factors_alloc *= 2;
		rb->factors = (uint32 *)xrealloc(rb->factors,
					rb->num_factors_alloc *
					sizeof(uint32));
	}
	f = rb->factors + rb->num_factors;

	for (i = 0; i < num_factors_r; i++)
		f[i] = factors_r[i];
	f += i;
	rb->num_factors += i;

	for (i = 0; i < num_factors_a; i++)
		f[i] = factors_a[i];
	f += i;
	rb->num_factors += i;

	/* add its large factors (if any) */

	if (unfactored_r->nwords > 1 || 
	    unfactored_r->val[0] > 1) {
		for (i = 0; i < unfactored_r->nwords; i++)
			f[i] = unfactored_r->val[i];
		f += i;
		rb->num_factors += i;
		c->lp_r_num_words = i;
	}

	if (unfactored_a->nwords > 1 || 
	    unfactored_a->val[0] > 1) {
		for (i = 0; i < unfactored_a->nwords; i++)
			f[i] = unfactored_a->val[i];
		f += i;
		rb->num_factors += i;
		c->lp_a_num_words = i;
	}
}
	
/*------------------------------------------------------------------*/
uint32 relation_batch_run(relation_batch_t *rb) {

	rb->num_success = 0;
	if (rb->num_relations > 0) {
		compute_remainder_tree(0, rb->num_relations - 1,
					rb, rb->prime_product);
	}

	/* wipe out batched relations */

	rb->num_relations = 0;
	rb->num_factors = 0;
	return rb->num_success;
}
