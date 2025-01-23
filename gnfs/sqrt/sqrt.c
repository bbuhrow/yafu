/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: sqrt.c 904 2013-07-04 12:00:48Z jasonp_sf $
--------------------------------------------------------------------*/

#include "sqrt.h"

/* we will need to find primes q for which f(x) mod q
   is irreducible. This is the maximum number of q to try */

#define NUM_PRIME_RETRIES 100

/*--------------------------------------------------------------------*/
uint32 get_prime_for_sqrt(mpz_poly_t *alg_poly,
			  uint32 min_value,
			  uint32 *q_out) {

	uint32 i;
	uint32 status = 0;
	uint32 q = 0;
	prime_sieve_t prime_sieve;

	init_prime_sieve(&prime_sieve, min_value, (uint32)(-1));
	for (i = 0; i < NUM_PRIME_RETRIES; i++) {
		uint32 tmp_q = get_next_prime(&prime_sieve);
		if (is_irreducible(alg_poly, tmp_q)) {
			q = tmp_q;
			break;
		}
		else if (q == 0) {
			/* in rare cases, there is no q for which alg_poly
			   mod q is irreducible. Life becomes much more
			   difficult in this case, since square roots mod
			   q will not be unique. The only alternative when
			   this happens is to pick a q such that alg_poly 
			   mod q has no linear roots (so that all of the
			   relations mod q are relatively prime to alg_poly), 
			   then keep trying dependencies until by luck we 
			   start with the correct initial square root for 
			   the Newton iteration to succeed.

			   Buhler et. al. show that for polynomial degree d,
			   on average one in 2^(d/2) dependencies will lead to
			   a congruence of squares (and about half of those
			   will lead to a factor). Technically we also need to
			   check that alg_poly mod q is squarefree, but that
			   would require a full polynomial factoring routine;
			   I'm gambling that being squarefree is not rare. */

			uint32 roots[MAX_POLY_DEGREE];
			uint32 high_coeff;
			uint32 num_roots = poly_get_zeros(roots, alg_poly,
						tmp_q, &high_coeff, 1);
			if (high_coeff != 0 && num_roots == 0)
				q = tmp_q;
		}
	}
	free_prime_sieve(&prime_sieve);

	if (i == NUM_PRIME_RETRIES)
		status = 1;

	*q_out = q;
	return status;
}

/*--------------------------------------------------------------------*/
static void eval_poly_derivative(mpz_poly_t *poly, 
				mpz_t m1, mpz_t m0, 
				mpz_t n, mpz_t res) {
	uint32 i;
	mpz_t m1_pow;
	mpz_t tmp;

	mpz_init_set(m1_pow, m1);
	mpz_init(tmp);

	i = poly->degree;
	mpz_mul_ui(res, poly->coeff[i], i);

	mpz_neg(m0, m0);
	while (--i > 1) {
		mpz_mul(res, res, m0);
		mpz_mul_ui(tmp, poly->coeff[i], i);
		mpz_addmul(res, tmp, m1_pow);
		mpz_mul(m1_pow, m1_pow, m1);
	}
	mpz_mul(res, res, m0);
	mpz_addmul(res, poly->coeff[i], m1_pow);
	mpz_neg(m0, m0);

	mpz_clear(m1_pow);
	mpz_clear(tmp);
	mpz_mod(res, res, n);
}

/*--------------------------------------------------------------------*/
typedef struct {
	uint64 p;
	uint64 count;
} rat_prime_t;

static uint32 rat_square_root(relation_t *rlist, uint32 num_relations,
				mpz_t n, mpz_t sqrt_r) {
	uint32 i, j, num_primes;
	hashtable_t h;
	uint32 already_seen;
	uint32 array_size;
	mpz_t base, exponent, tmp;
	uint32 status = 0;
	rat_prime_t *curr;

	mpz_init(base);
	mpz_init(exponent);
	mpz_init(tmp);

	/* count up the number of times each prime factor in
	   rlist occurs */

	hashtable_init(&h, (uint32)WORDS_IN(rat_prime_t), 
				(uint32)WORDS_IN(uint64));

	for (i = 0; i < num_relations; i++) {
		relation_t *r = rlist + i;

		for (j = array_size = 0; j < r->num_factors_r; j++) {
			uint64 p = decompress_p(r->factors, &array_size);
			curr = (rat_prime_t *)hashtable_find(&h, &p, NULL,
							    &already_seen);
			if (!already_seen)
				curr->count = 1;
			else
				curr->count++;
		}
	}

	/* verify all such counts are even, and form the 
	   rational square root */

	mpz_set_ui(sqrt_r, 1);
	num_primes = hashtable_get_num(&h);
	curr = hashtable_get_first(&h);

	for (i = 0; i < num_primes; i++) {
		uint64 p = curr->p;
		uint64 count = curr->count;

		if (count % 2) {
			status = 1;
			break;
		}
		if (p > 0 && count > 0) {
			uint64_2gmp(p, base);
			uint64_2gmp(count / 2, exponent);
			mpz_powm(tmp, base, exponent, n);
			mpz_mul(sqrt_r, sqrt_r, tmp);
			mpz_tdiv_r(sqrt_r, sqrt_r, n);
		}
		curr = hashtable_get_next(&h, curr);
	}

	hashtable_free(&h);
	mpz_clear(base);
	mpz_clear(exponent);
	mpz_clear(tmp);
	return status;
}

/*--------------------------------------------------------------------*/
/* we will not do any computations involving the count,
   only verifying that it is even. Thus we can get away
   with storing only the low word of the count */

typedef struct {
	ideal_t ideal;
	uint32 count;
} alg_prime_t;

static uint32 verify_alg_ideal_powers(relation_t *rlist, 
					uint32 num_relations,
					uint32 *num_free_relations) {

	uint32 i, j, num_ideals;
	hashtable_t h;
	uint32 already_seen;
	alg_prime_t *curr;
	uint32 status = 0;

	/* count the multiplicity of each algebraic ideal (not
	   just the prime to which the ideal corresponds) in rlist */

	*num_free_relations = 0;

	hashtable_init(&h, (uint32)WORDS_IN(alg_prime_t),
			(uint32)WORDS_IN(ideal_t));

	for (i = 0; i < num_relations; i++) {
		relation_t *r = rlist + i;
		relation_lp_t rlp;

		find_large_ideals(r, &rlp, 0, 0);

		for (j = 0; j < rlp.ideal_count; j++) {
			ideal_t *curr_ideal = rlp.ideal_list + j;

			if (curr_ideal->rat_or_alg == RATIONAL_IDEAL)
				continue;

			curr = (alg_prime_t *)hashtable_find(&h, curr_ideal, 
						NULL, &already_seen);

			if (!already_seen)
				curr->count = 1;
			else
				curr->count++;
		}

		if (r->b == 0)
			(*num_free_relations)++;
	}

	/* verify each ideal occurs an even number of times */

	num_ideals = hashtable_get_num(&h);
	curr = hashtable_get_first(&h);

	for (i = 0; i < num_ideals; i++) {
		if (curr->count % 2) {
			status = 1;
			break;
		}
		curr = hashtable_get_next(&h, curr);
	}

	hashtable_free(&h);
	return status;
}

/*--------------------------------------------------------------------*/
uint32 nfs_find_factors(msieve_obj *obj, mpz_t n, 
			factor_list_t *factor_list) {

	/* external interface for the NFS square root */

	uint32 i, j;
	uint32 check_q;
	factor_base_t fb;
	mpz_poly_t monic_alg_poly;
	mpz_poly_t *rpoly;
	mpz_poly_t *apoly;
	mpz_t exponent, sqrt_r, sqrt_a;
	mpz_t c, tmp1, tmp2;
	uint32 dep_lower = 1;
	uint32 dep_upper = 64;
	uint32 factor_found = 0;
	time_t cpu_time;

	logprintf(obj, "\n");
	logprintf(obj, "commencing square root phase\n");

	memset(&fb, 0, sizeof(fb));
	apoly = &fb.afb.poly;
	rpoly = &fb.rfb.poly;
	mpz_poly_init(rpoly);
	mpz_poly_init(apoly);
	mpz_poly_init(&monic_alg_poly);
	mpz_init(exponent);
	mpz_init(sqrt_r);
	mpz_init(sqrt_a);
	mpz_init(c);
	mpz_init(tmp1);
	mpz_init(tmp2);

	/* read in the NFS polynomials */

	cpu_time = time(NULL);
	if (read_poly(obj, n, rpoly, apoly, NULL)) {
		logprintf(obj, "polynomials not found\n");
		goto finished;
	}

	/* find the values needed to convert the algebraic 
	   square root back to an integer */

	if (rpoly->degree != 1) {
		logprintf(obj, "cannot handle non-linear polynomials\n");
		goto finished;
	}

	/* construct a monic version of the algebraic poly,
	   saving off the leading coefficient separately */

	j = apoly->degree;
	if (mpz_cmp_ui(apoly->coeff[j], 0) < 0) {
		logprintf(obj, "cannot handle negative leading "
				"algebraic polynomial coefficient\n");
		goto finished;
	}

	mpz_set(c, apoly->coeff[j]);
	mpz_set(tmp1, c);
	mpz_set(monic_alg_poly.coeff[j-1], apoly->coeff[j-1]);
	monic_alg_poly.degree = j;
	mpz_set_ui(monic_alg_poly.coeff[j], 1);

	for (i = j - 2; (int32)i >= 0; i--) {
		mpz_mul(monic_alg_poly.coeff[i], apoly->coeff[i], tmp1);
		if (i > 0)
			mpz_mul(tmp1, tmp1, c);
	}
	get_prime_for_sqrt(&monic_alg_poly, (uint32)0x80000000, &check_q);

	/* determine the list of dependencies to compute */

	if (obj->nfs_args != NULL) {

		const char *tmp;
		const char *lower_limit;
		const char *upper_limit;

		tmp = strstr(obj->nfs_args, "dep_first=");
		if (tmp != NULL)
			dep_lower = strtoul(tmp + 10, NULL, 10);

		tmp = strstr(obj->nfs_args, "dep_last=");
		if (tmp != NULL)
			dep_upper = strtoul(tmp + 9, NULL, 10);

		/* old-style 'X,Y' format */

		upper_limit = strchr(obj->nfs_args, ',');
		if (upper_limit != NULL) {
			lower_limit = upper_limit - 1;
			while (lower_limit > obj->nfs_args &&
				isdigit(lower_limit[-1])) {
				lower_limit--;
			}
			upper_limit++;
			dep_lower = strtoul(lower_limit, NULL, 10);
			dep_upper = strtoul(upper_limit, NULL, 10);
		}

		dep_lower = MAX(dep_lower, 1);
		dep_upper = MAX(dep_upper, 1);
		dep_lower = MIN(dep_lower, 64);
		dep_upper = MIN(dep_upper, 64);
		logprintf(obj, "handling dependencies %u to %i\n",
				dep_lower, dep_upper);
	}

	/* for each dependency */

	for (i = dep_lower; i <= dep_upper; i++) {

		uint32 num_relations;
		uint32 num_free_relations;
		relation_t *rlist;
		abpair_t *abpairs;

		logprintf(obj, "reading relations for dependency %u\n", i);

		/* read in only the relations for dependency i */

		nfs_read_cycles(obj, &fb, NULL, NULL,
				&num_relations, &rlist, 0, i);

		if (num_relations == 0)
			continue;

		/* do some sanity checking, performing increasing
		   amounts of work as the dependency proves itself
		   to be valid */

		if (num_relations % 2) {
			/* the LA is supposed to force the number of 
			   relations in the dependency to be even. 
			   This isn't necessary if the leading coeff of
			   both NFS polynomials are squares, or if both
			   NFS polynomials are monic, since the 
			   corrections below that need the number of 
			   relations are avoided. But only a small 
			   minority of NFS jobs would satisfy this condition */

			logprintf(obj, "number of relations is not even\n");
			nfs_free_relation_list(rlist, num_relations);
			continue;
		}
		if (verify_alg_ideal_powers(rlist, 
				num_relations, &num_free_relations) != 0) {
			logprintf(obj, "algebraic side is not a square!\n");
			nfs_free_relation_list(rlist, num_relations);
			continue;
		}
		if (num_free_relations % 2) {
			logprintf(obj, "number of free relations (%u) is "
					"not even\n", num_free_relations);
			nfs_free_relation_list(rlist, num_relations);
			continue;
		}
		if (rat_square_root(rlist, num_relations, n, sqrt_r) != 0) {
			logprintf(obj, "rational side is not a square!\n");
			nfs_free_relation_list(rlist, num_relations);
			continue;
		}

		/* flatten the list of relations; each occurrence of
		   a relation gets its own abpair_t */

		abpairs = (abpair_t *)xmalloc(num_relations *
						sizeof(abpair_t));
		for (j = 0; j < num_relations; j++) {
			abpairs[j].a = rlist[j].a;
			abpairs[j].b = rlist[j].b;
		}
		nfs_free_relation_list(rlist, num_relations);

		/* perform the major work: the algebraic square root.
		   Note that to conserve memory, abpairs is freed in
		   the following call */

		mpz_set_ui(sqrt_a, 0);
		alg_square_root(obj, &monic_alg_poly, n, c, 
				rpoly->coeff[1], rpoly->coeff[0], 
				abpairs, num_relations, check_q, sqrt_a);
		if (mpz_sgn(sqrt_a) == 0) {
			logprintf(obj, "algebraic square root failed\n");
			continue;
		}

		/* an algebraic square root is available; move on
		   to the final congruence of squares. The arithmetic
		   is as given in Buhler et. al. with one exception:
		   when the rational poly is nonmonic there is a 
		   correction to the final square root value but the 
		   free relations *do not* figure into it. This latter
		   point is completely ignored in the literature! */

		eval_poly_derivative(apoly, rpoly->coeff[1], 
					rpoly->coeff[0], n, tmp1);
		mpz_mul(sqrt_r, sqrt_r, tmp1);
		mpz_mod(sqrt_r, sqrt_r, n);

		mpz_set_ui(exponent, 0);
		if (mpz_cmp_ui(c, 1) != 0) {
			mpz_set_ui(exponent, num_relations / 2 + 
						apoly->degree - 2);
			mpz_powm(tmp1, c, exponent, n);
			mpz_mul(sqrt_r, sqrt_r, tmp1);
			mpz_mod(sqrt_r, sqrt_r, n);
		}

		if (mpz_cmp_ui(rpoly->coeff[1], 1) != 0) {
			mpz_set_ui(exponent, (num_relations -
				       		num_free_relations) / 2);
			mpz_set(tmp1, rpoly->coeff[1]);
			if (mpz_sgn(tmp1) < 0)
				mpz_add(tmp1, tmp1, n);

			mpz_powm(tmp2, tmp1, exponent, n);
			mpz_mul(sqrt_a, sqrt_a, tmp2);
			mpz_mod(sqrt_a, sqrt_a, n);
		}

		/* a final sanity check: square the rational and algebraic 
		   square roots, expecting the same value modulo n */

		mpz_mul(tmp1, sqrt_r, sqrt_r);
		mpz_mul(tmp2, sqrt_a, sqrt_a);
		mpz_mod(tmp1, tmp1, n);
		mpz_mod(tmp2, tmp2, n);
		if (mpz_cmp(tmp1, tmp2) != 0) {
			logprintf(obj, "dependency does not form a "
					"congruence of squares!\n");
			continue;
		}

		/* look for a nontrivial factor of n */

		mpz_add(tmp1, sqrt_r, sqrt_a);
		mpz_gcd(tmp1, tmp1, n);
		if (mpz_cmp_ui(tmp1, 1) == 0) {
			logprintf(obj, "GCD is 1, no factor found\n");
		}
		else if (mpz_cmp(tmp1, n) == 0) {
			logprintf(obj, "GCD is N, no factor found\n");
		}
		else {
			/* factor found; add it to the list of factors. 
			   Stop trying dependencies if the remaining
			   composite is small enough that another method
			   will factor it faster.

			   Actually, we should be stopping when the remaining
			   composite is much larger (70-80 digits), but 
			   avoid doing this because the MPQS code will run
			   and wipe out all the NFS relations we've collected */

			uint32 composite_bits;
			mp_t junk;

			gmp2mp(tmp1, &junk);
			composite_bits = factor_list_add(obj, 
						factor_list, &junk);

			factor_found = 1;
			if (composite_bits < SMALL_COMPOSITE_CUTOFF_BITS) {
				break;
			}
			else {
				/* a single dependency could take hours,
				   and if N has more than two factors then
				   we'll need several dependencies to find
				   them all. So at least report the smallest
				   cofactor that we just found */

				mpz_divexact(tmp2, n, tmp1);
				gmp_sprintf(obj->mp_sprintf_buf, "%Zd",
						(mpz_cmp(tmp1, tmp2) < 0) ? 
						tmp1 : tmp2);
				logprintf(obj, "found factor: %s\n",
						obj->mp_sprintf_buf);
			}
		}
	}

finished:
	cpu_time = time(NULL) - cpu_time;
	logprintf(obj, "sqrtTime: %u\n", (uint32)cpu_time);

	mpz_poly_free(&fb.rfb.poly);
	mpz_poly_free(&fb.afb.poly);
	mpz_poly_free(&monic_alg_poly);
	mpz_clear(exponent);
	mpz_clear(sqrt_r);
	mpz_clear(sqrt_a);
	mpz_clear(c);
	mpz_clear(tmp1);
	mpz_clear(tmp2);

	return factor_found;
}
