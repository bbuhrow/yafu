/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: sqrt_a.c 887 2013-06-08 01:43:53Z jasonp_sf $
--------------------------------------------------------------------*/

#include "sqrt.h"

	/* This code computes the algebraic square root by
	   brute force. Given a collection of relations, each
	   of which is a polynomial p(x) = a*c-b*x for constant c
	   and various (a,b):

	   - Multiply all the relations together, modulo the
	     monic version of the algebraic polynomial f(x), 
	     of degree d. The result is a polynomial S(x) of degree
	     d-1 with arbitrary precision (very large) coefficients

	   - Find the polynomial that equals S(x) when squared
	     modulo f(x), using q-adic Newton iteration. This will
	     have coefficients about half the size of those in S(x)

	   - Convert the result to an integer modulo n by substituting
	     the root m of f(x) modulo n for x

	   None of the basic NFS papers consider this method to be
	   a viable option, but all of those papers were written
	   in the early to mid 1990s and computers were expensive and
	   limited then. On a modern machine, the memory consumed
	   by the direct method is not excessive (it's comparable to
	   the memory needed for the NFS linear algebra) and the 
	   runtime for optimized code using FFT-based multiplication
	   is quite reasonable. The big advantage of this method is
	   its simplicity; properly implementing the more standard
	   algorithms by Montgomery / Nguyen for the algebraic square
	   root would require approximately 100x as much code, plus
	   large amounts of very difficult algebraic number theory */
	     
/* When an mpz_poly_t is monic, the highest-order coefficient
   (i.e. one) is considered implicit, and is *not* reflected
   in the degree */

/* bag of quantities needed for computing S(x) */

typedef struct {
	mpz_poly_t *monic_poly;
	abpair_t *rlist;
	mpz_t c;
} relation_prod_t;

/*-------------------------------------------------------------------*/
static void mpz_poly_mod_q(mpz_poly_t *p, mpz_t q, mpz_poly_t *res) {

	uint32 i;
	uint64 pbits, resbits;

	/* also trim aggressively the memory use of the
	   computed remainders; the algebraic square root has
	   comparatively few arithmetic operations but the
	   memory use for large problems is a concern */

	for (i = 0; i <= p->degree; i++) {
		pbits = mpz_sizeinbase(p->coeff[i], 2);
		mpz_fdiv_r(res->coeff[i], p->coeff[i], q);
		resbits = mpz_sizeinbase(res->coeff[i], 2);

		if (pbits > resbits + 1000)
			mpz_realloc2(res->coeff[i], resbits + 500);
	}

	/* recalculate the degree */

	for (i = p->degree; i; i--) {
		if (mpz_sgn(res->coeff[i]) != 0)
			break;
	}
	res->degree = i;
}

/*-------------------------------------------------------------------*/
static void mpz_poly_bits(mpz_poly_t *p, 
			uint64 *total_bits, uint64 *max_bits) {

	uint32 i;
	uint64 bits1, bits2;

	for (i = bits1 = bits2 = 0; i <= p->degree; i++) {
		uint64 curr_bits = mpz_sizeinbase(p->coeff[i], 2);
		bits1 += curr_bits;
		bits2 = MAX(bits2, curr_bits);
	}

	*total_bits = bits1;
	*max_bits = bits2;
}

/*-------------------------------------------------------------------*/
static void mpz_poly_monic_derivative(mpz_poly_t *src, mpz_poly_t *dest) {
	
	uint32 i;

	/* compute the coefficients of the derivative of src,
	   assumed to be monic */

	for (i = 0; i < src->degree; i++) {
		mpz_mul_ui(dest->coeff[i], src->coeff[i+1], 
				(unsigned long)(i+1));
	}
	mpz_set_ui(dest->coeff[i], (unsigned long)(i+1));

	dest->degree = src->degree;
	for (i++; i <= MAX_POLY_DEGREE; i++)
		mpz_realloc2(dest->coeff[i], 1);
}

/*-------------------------------------------------------------------*/
static void mpz_poly_mul(mpz_poly_t *p1, mpz_poly_t *p2,
			mpz_poly_t *mod, uint32 free_p2) {

	/* multiply p1(x) by p2(x) modulo mod(x) (assumed monic)
	   If free_p2 is nonzero the coefficients of p2(x) are 
	   freed after being used */

	uint32 i, j;
	uint32 d = mod->degree;
	uint32 d1 = p1->degree;
	uint32 d2 = p2->degree;
	uint32 prod_degree;
	mpz_t tmp[MAX_POLY_DEGREE + 1];

	/* initialize */

	for (i = 0; i < MAX_POLY_DEGREE + 1; i++)
		mpz_init_set_ui(tmp[i], (unsigned long)0);

	/* multiply p1 by the leading coefficient of p2 */

	for (i = 0; i <= d1; i++) {
		mpz_mul(tmp[i], p1->coeff[i], p2->coeff[d2]);
	}
	prod_degree = d1;
	if (free_p2) {
		mpz_realloc2(p2->coeff[d2], 1);
	}

	/* for each of the other coefficients in p2 */

	for (i = d2 - 1; (int32)i >= 0; i--) {

		/* shift the accumulator up by one, bubbling
		   the highest-order coefficient to the lowest */

		for (j = prod_degree; (int32)j >= 0; j--)
			mpz_swap(tmp[j+1], tmp[j]);

		/* add in the product of p1(x) and coefficient
		   i of p2 */

		for (j = d1; j; j--) {
			mpz_addmul(tmp[j], p1->coeff[j], p2->coeff[i]);
		}
		mpz_mul(tmp[j], p1->coeff[j], p2->coeff[i]);
		if (free_p2) {
			mpz_realloc2(p2->coeff[i], 1);
		}

		/* recalculate the degree of the result */

		prod_degree = d + 1;
		while (prod_degree && mpz_sgn(tmp[prod_degree]) == 0)
			prod_degree--;

		/* if it exceeds the degree of mod(x), subtract
		   mod(x) * (leading accumulator coefficient) */

		if (prod_degree <= d)
			continue;

		for (j = d; (int32)j >= 0; j--) {
			mpz_submul(tmp[j], mod->coeff[j], tmp[prod_degree]);
		}
		prod_degree--;
	}

	/* move the result in the accumulator over to p1 */

	for (i = 0; i <= prod_degree; i++) {
		mpz_swap(p1->coeff[i], tmp[i]);
		mpz_clear(tmp[i]);
	}
	for (; i < MAX_POLY_DEGREE + 1; i++)
		mpz_clear(tmp[i]);

	/* recalculate the degree */

	i = prod_degree;
	while (i > 0 && mpz_sgn(p1->coeff[i]) == 0) {
		mpz_realloc2(p1->coeff[i], 1);
		i--;
	}
	p1->degree = i;
}

/*-------------------------------------------------------------------*/
static uint32 verify_product(mpz_poly_t *gmp_prod, abpair_t *abpairs, 
			uint32 num_relations, uint32 q, mpz_t c, 
			mpz_poly_t *alg_poly) {

	/* a sanity check on the computed value of S(x): for
	   a small prime q for which alg_poly is irreducible,
	   verify that gmp_prod mod q equals the product
	   mod q of the relations in abpairs[]. The latter can
	   be computed very quickly */

	uint32 i, j;
	uint32 c_mod_q = mpz_fdiv_ui(c, q);
	uint32 d = alg_poly->degree;
	uint32 ref_prod[MAX_POLY_DEGREE];
	uint32 prod[MAX_POLY_DEGREE];
	uint32 mod[MAX_POLY_DEGREE];
	uint32 accum[MAX_POLY_DEGREE + 1];

	/* compute the product mod q directly. First initialize
	   and reduce the coefficients of alg_poly and gmp_prod mod q */

	for (i = 0; i <= d; i++) {
		prod[i] = 0;
		ref_prod[i] = mpz_fdiv_ui(gmp_prod->coeff[i],
					(unsigned long)q);
		mod[i] = mpz_fdiv_ui(alg_poly->coeff[i],
					(unsigned long)q);
	}
	prod[0] = 1;

	/* multiply the product by each relation in
	   turn, modulo q */

	for (i = 0; i < num_relations; i++) {
		int64 a = abpairs[i].a;
		uint32 b = q - (abpairs[i].b % q);
		uint32 ac;

		a = a % (int64)q;
		if (a < 0)
			a += q;
		ac = mp_modmul_1((uint32)a, c_mod_q, q);

		for (j = accum[0] = 0; j <= d; j++) {
			accum[j+1] = mp_modmul_1(prod[j], b, q);
			accum[j] = mp_modadd_1(accum[j],
					mp_modmul_1(ac, prod[j], q), q);
		}

		for (j = 0; j <= d; j++) {
			prod[j] = mp_modsub_1(accum[j],
					mp_modmul_1(accum[d+1], 
						mod[j], q), q);
		}
	}

	/* do the polynomial compare */

	for (i = 0; i <= d; i++) {
		if (ref_prod[i] != prod[i])
			break;
	}
	if (i > d)
		return 1;
	return 0;
}

/*-------------------------------------------------------------------*/
static void relation_to_poly(abpair_t *abpair, mpz_t c,
				mpz_poly_t *poly) {

	/* given a and b in abpair, along with c, compute
	   poly(x) = a*c - b*x */

	mpz_set(poly->coeff[0], c);
	int64_2gmp(abpair->a, poly->coeff[1]);
	mpz_mul(poly->coeff[0], poly->coeff[0], poly->coeff[1]);
	mpz_set_ui(poly->coeff[1], (unsigned long)abpair->b);

	poly->degree = 0;
	if (abpair->b > 0) {
		poly->degree = 1;
		mpz_neg(poly->coeff[1], poly->coeff[1]);
	}
}

/*-------------------------------------------------------------------*/
static void multiply_relations(relation_prod_t *prodinfo, 
			uint32 index1, uint32 index2,
			mpz_poly_t *prod) {

	/* multiply together the relations from index1 
	   to index2, inclusive. We proceed recursively to
	   assure that polynomials with approximately equal-
	   size coefficients get multiplied, and also to
	   avoid wasting huge amounts of memory in the
	   beginning when all the polynomials are small
	   but the memory allocated for them is large */

	uint32 i;
	mpz_poly_t prod1, prod2;

	if (index1 == index2) {
		/* base case of recursion */

		relation_to_poly(prodinfo->rlist + index1,
				 prodinfo->c, prod);
		return;
	}

	mpz_poly_init(&prod1);
	mpz_poly_init(&prod2);
	
	if (index1 == index2 - 1) {
		/* base case of recursion */

		relation_to_poly(prodinfo->rlist + index1,
				 prodinfo->c, &prod1);
		relation_to_poly(prodinfo->rlist + index2,
				 prodinfo->c, &prod2);
	}
	else {
		/* recursively compute the product of the first
		   half and the last half of the relations */

		uint32 mid = (index1 + index2) / 2;
		multiply_relations(prodinfo, index1, mid, &prod1);
		multiply_relations(prodinfo, mid + 1, index2, &prod2);
	}

	/* multiply them together and save the result */
	mpz_poly_mul(&prod1, &prod2, prodinfo->monic_poly, 1);

	for (i = 0; i <= prod1.degree; i++)
		mpz_swap(prod->coeff[i], prod1.coeff[i]);
	prod->degree = prod1.degree;
	mpz_poly_free(&prod1);
	mpz_poly_free(&prod2);
}

/*-------------------------------------------------------------------*/
#define ISQRT_NUM_ATTEMPTS 10

static uint32 get_initial_inv_sqrt(msieve_obj *obj, mpz_poly_t *alg_poly,
				mpz_poly_t *prod, mpz_poly_t *isqrt_mod_q, 
				mpz_t q_out) {

	/* find the prime q_out and the initial value of the
	   reciprocal square root of prod(x) mod q_out to use 
	   for the Newton iteration */

	uint32 i;
	uint32 q, start_q;

	alg_poly->degree++;

	/* find a prime q for which mp_alg_poly mod q is
	   irreducible. The starting value to try was passed in */

	start_q = mpz_get_ui(q_out);

	for (i = 0; i < ISQRT_NUM_ATTEMPTS; i++) {
		if (get_prime_for_sqrt(alg_poly, start_q + 1, &q)) {
			if (start_q > 150) {
				logprintf(obj, "warning: no irreducible prime "
					"found, switching to small primes\n");
				start_q = 50;
				/* for octics, even mod 13 works and 
				   is resonably fast */
				if (alg_poly->degree > 6) 
					start_q = 12;
				continue;
			}
		}

		/* find the reciprocal square root mod q, or try
		   another q if this fails */

		if (inv_sqrt_mod_q(isqrt_mod_q, prod, 
				alg_poly, q, &obj->seed1, &obj->seed2)) {
			break;
		}
		start_q = q;
	}

	alg_poly->degree--;

	if (i == ISQRT_NUM_ATTEMPTS) {
		logprintf(obj, "error: cannot recover square root mod q\n");
		return 0;
	}

	logprintf(obj, "initial square root is modulo %u\n", q);
	mpz_set_ui(q_out, (unsigned long)q);
	return 1;
}

/*-------------------------------------------------------------------*/
static uint32 get_final_sqrt(msieve_obj *obj, mpz_poly_t *alg_poly,
			mpz_poly_t *prod, mpz_poly_t *isqrt_mod_q, 
			mpz_t q) {

	/* the main q-adic Newton iteration. On input, isqrt_mod_q
	   contains the starting value of the reciprocal square
	   root R[0](x) of the polynomial prod(x). The iteration is

	   R[k](x) = R[k-1](x) * (3 - prod(x)*R[k-1](x)^2) / 2 mod (q^(2^k))

	   and at the end of iteration k, prod(x)*R[k-1](x)^2 mod (q^(2^k))
	   is 1. We keep iterating until q^(2^k) is larger than the
	   size of the coefficients of the square root (i.e. about half
	   the size of the coefficients of prod(x)). Then the square
	   root to use is R[k](x) * prod(x) mod (q^(2^k)), which is
	   written to isqrt_mod_q */

	uint32 i, j;
	uint64 prod_bits, prod_max_bits;
	uint32 num_iter;

	/* initialize */

	mpz_poly_bits(prod, &prod_bits, &prod_max_bits);

	/* since prod(x) only matters mod q^(2^(final_k)), we can
	   cut the memory use in half by changing prod(x) to this.
	   Remember final_k as well */

	i = mpz_get_ui(q);
	for (num_iter = 0; mpz_sizeinbase(q, 2) < 
				prod_max_bits / 2 + 4000; num_iter++) {
		mpz_mul(q, q, q);
	}

	mpz_poly_mod_q(prod, q, prod);
	mpz_set_ui(q, (unsigned long)i);
	mpz_realloc2(q, 33);

	/* do the main iteration */

	for (i = 0; i < num_iter; i++) {

		mpz_poly_t tmp_poly;

		/* square the previous modulus */

		mpz_mul(q, q, q);

		/* compute prod(x) * (previous R)^2 */

		mpz_poly_init(&tmp_poly);
		mpz_poly_mod_q(prod, q, &tmp_poly);
		mpz_poly_mul(&tmp_poly, isqrt_mod_q, alg_poly, 0);
		mpz_poly_mod_q(&tmp_poly, q, &tmp_poly);
		mpz_poly_mul(&tmp_poly, isqrt_mod_q, alg_poly, 0);
		mpz_poly_mod_q(&tmp_poly, q, &tmp_poly);

		/* compute ( (3 - that) / 2 ) mod q */

		mpz_sub_ui(tmp_poly.coeff[0], tmp_poly.coeff[0], 
				(unsigned long)3);

		for (j = 0; j <= tmp_poly.degree; j++) {

			mpz_t *c = tmp_poly.coeff + j;

			if (mpz_sgn(*c) != 0) {
				mpz_neg(*c, *c);
				if (mpz_tstbit(*c, (unsigned long)0))
					mpz_add(*c, *c, q);
				mpz_tdiv_q_2exp(*c, *c, (unsigned long)1);
			}
		}

		/* finally, compute the new R(x) by multiplying the
		   result above by the old R(x) */

		mpz_poly_mul(&tmp_poly, isqrt_mod_q, alg_poly, 1);
		mpz_poly_mod_q(&tmp_poly, q, isqrt_mod_q);
		mpz_poly_free(&tmp_poly);
	}

	/* attempt to compute the square root. 
	   First multiply R(x) by prod(x), deleting prod(x) 
	   since we won't need it beyond this point */

	mpz_poly_mul(isqrt_mod_q, prod, alg_poly, 1);
	mpz_poly_mod_q(isqrt_mod_q, q, isqrt_mod_q);

	/* this is a little tricky. Up until now we've
	   been working modulo big numbers, but the coef-
	   ficients of the square root are just integers,
	   and may be negative. Negative numbers mod q
	   have a numerical value near that of +q, but we
	   want the square root to have a negative coef-
	   ficient in that case. Hence, if the top
	   few words of any coefficent of the square root
	   match the top few words of q, we assume this
	   coefficient is negative and subtract q from it.

	   Theoretically we could be wrong, and the 
	   coefficient really is supposed to be a big 
	   positive number near q in size. However, if
	   q is thousands of bits larger than the size we
	   expect for the square root coefficients, this
	   is so unlikely that it's not worth worrying about */

	for (i = 0; i <= isqrt_mod_q->degree; i++) {
		mpz_t *c = isqrt_mod_q->coeff + i;
		size_t limbs = mpz_size(*c);

		if (limbs == mpz_size(q) &&
		    mpz_getlimbn(*c, (mp_size_t)(limbs-1)) ==
			mpz_getlimbn(q, (mp_size_t)(limbs-1)) &&
		    mpz_getlimbn(*c, (mp_size_t)(limbs-2)) ==
			mpz_getlimbn(q, (mp_size_t)(limbs-2)) &&
		    mpz_getlimbn(*c, (mp_size_t)(limbs-3)) ==
			mpz_getlimbn(q, (mp_size_t)(limbs-3))) { 
			mpz_sub(*c, *c, q);
		}
	}

	/* another heuristic: we will assume the Newton
	   iteration has converged if, after applying the
	   correction above for negative square root
	   coefficients, the total number of bits in the 
	   coefficients of the resulting polynomial is
	   much smaller than we would expect from random
	   polynomials modulo q */

	mpz_poly_bits(isqrt_mod_q, &prod_bits, &prod_max_bits);
	if (prod_bits >= (isqrt_mod_q->degree + 1) * 
				mpz_sizeinbase(q, 2) - 100) {
		logprintf(obj, "Newton iteration failed to converge\n");
		return 0;
	}
	return 1;
}

/*-------------------------------------------------------------------*/
static void convert_to_integer(mpz_poly_t *alg_sqrt, mpz_t n,
				mpz_t c, mpz_t m1, 
				mpz_t m0, mpz_t res) {

	/* given the completed square root, apply the homomorphism
	   to convert the polynomial to an integer. We do this
	   by evaluating alg_sqrt at -c*m0/m1, with all calculations
	   performed mod n */

	uint32 i;
	mpz_t m1_pow;
	mpz_t m0c;

	mpz_init_set(m1_pow, m1);
	mpz_poly_mod_q(alg_sqrt, n, alg_sqrt);

	mpz_init_set_ui(m0c, 0);
	mpz_submul(m0c, m0, c);
	mpz_mod(m0c, m0c, n);

	i = alg_sqrt->degree;
	mpz_set(res, alg_sqrt->coeff[i]);
	while (--i) {
		mpz_mul(res, res, m0c);
		mpz_addmul(res, m1_pow, alg_sqrt->coeff[i]);
		mpz_mul(m1_pow, m1_pow, m1);
	}
	mpz_mul(res, res, m0c);
	mpz_addmul(res, m1_pow, alg_sqrt->coeff[i]);
	mpz_mod(res, res, n);

	mpz_clear(m1_pow);
	mpz_clear(m0c);
}

/*-------------------------------------------------------------------*/
void alg_square_root(msieve_obj *obj, mpz_poly_t *alg_poly, 
			mpz_t n, mpz_t c, mpz_t m1, 
			mpz_t m0, abpair_t *rlist, 
			uint32 num_relations, uint32 check_q,
			mpz_t sqrt_a) {
	
	/* external interface for computing the algebraic
	   square root */

	uint32 i;
	mpz_poly_t d_alg_poly;
	mpz_poly_t prod;
	mpz_poly_t alg_sqrt;
	relation_prod_t prodinfo;
	double log2_prodsize;
	mpz_t q;

	/* initialize */

	mpz_init(q);
	mpz_poly_init(&d_alg_poly);
	mpz_poly_init(&prod);
	mpz_poly_init(&alg_sqrt);

	if (mpz_cmp_ui(alg_poly->coeff[alg_poly->degree], 1) != 0) {
		printf("error: sqrt requires input poly to be monic\n");
		exit(-1);
	}
	alg_poly->degree--;

	/* multiply all the relations together */

	prodinfo.monic_poly = alg_poly;
	prodinfo.rlist = rlist;
	mpz_init_set(prodinfo.c, c);

	logprintf(obj, "multiplying %u relations\n", num_relations);
	multiply_relations(&prodinfo, 0, num_relations - 1, &prod);
	logprintf(obj, "multiply complete, coefficients have about "
			"%3.2lf million bits\n",
			(double)mpz_sizeinbase(prod.coeff[0], 2) / 1e6);

	/* perform a sanity check on the result */

	i = verify_product(&prod, rlist, num_relations, 
				check_q, c, alg_poly);
	free(rlist);
	mpz_clear(prodinfo.c);
	if (i == 0) {
		logprintf(obj, "error: relation product is incorrect\n");
		goto finished;
	}

	/* multiply by the square of the derivative of alg_poly;
	   this will guarantee that the square root of prod actually 
	   is an element of the number field defined by alg_poly.
	   If we didn't do this, we run the risk of the main Newton
	   iteration not converging */

	mpz_poly_monic_derivative(alg_poly, &d_alg_poly);
	mpz_poly_mul(&d_alg_poly, &d_alg_poly, alg_poly, 0);
	mpz_poly_mul(&prod, &d_alg_poly, alg_poly, 1);

	/* pick the initial small prime to start the Newton iteration.
	   To save both time and memory, choose an initial prime 
	   such that squaring it a large number of times will produce
	   a value just a little larger than we need to calculate
	   the square root.
	
	   Note that contrary to what some authors write, pretty much
	   any starting prime is okay. The Newton iteration has a division
	   by 2, so that 2 must be invertible mod the prime (this is
	   guaranteed for odd primes). Also, the Newton iteration will
	   fail if both square roots have the same value mod the prime;
	   however, even a 16-bit prime makes this very unlikely */

	i = mpz_size(prod.coeff[0]);
	log2_prodsize = (double)GMP_LIMB_BITS * (i - 2) +
			log(mpz_getlimbn(prod.coeff[0], (mp_size_t)(i-1)) *
				pow(2.0, (double)GMP_LIMB_BITS) +
			    mpz_getlimbn(prod.coeff[0], (mp_size_t)(i-2))) / 
			M_LN2 + 10000;

	while (log2_prodsize > 31.5)
		log2_prodsize *= 0.5;

	mpz_set_d(q, (uint32)pow(2.0, log2_prodsize) + 1);

	/* get the initial inverse square root */

	if (!get_initial_inv_sqrt(obj, alg_poly, 
				&prod, &alg_sqrt, q)) {
		goto finished;
	}

	/* compute the actual square root */

	if (get_final_sqrt(obj, alg_poly, &prod, &alg_sqrt, q))
		convert_to_integer(&alg_sqrt, n, c, m1, m0, sqrt_a);

finished:
	mpz_poly_free(&prod);
	mpz_poly_free(&alg_sqrt);
	mpz_poly_free(&d_alg_poly);
	mpz_clear(q);
	alg_poly->degree++;
}
