/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: polyutil.c 897 2013-06-22 13:16:18Z jasonp_sf $
--------------------------------------------------------------------*/

#include "poly.h"

/*------------------------------------------------------------------*/
static void get_bernstein_combined_score(poly_select_t *poly,
					double root_score) {

	/* An interesting problem: you have a root score and
	   a size score for each polynomial; how do you combine
	   them to get an overall rating? The value of root_score
	   implies that values of R(x)*A(x) effectively have size
	   R(x)*A(x)*exp(root_score) after dividing out the small
	   primes from sieve values. Technically this will
	   overestimate the size, because root_score includes a 
	   bias that accounts for random numbers; but this bias
	   is the same for all polynomials. The exponential term 
	   is constant, so pulling this expression out of the 
	   integrand used to calculate the size above gives: */

	poly->size_score *= pow(exp(root_score), -2.0 / 
				(poly->rpoly.degree + poly->apoly.degree));
}

/*------------------------------------------------------------------*/
static void poly_select_copy(poly_select_t *dest,
		  poly_select_t *src) {

	uint32 i;

	dest->size_score = src->size_score;
	dest->root_score = src->root_score;
	dest->combined_score = src->combined_score;
	dest->skewness = src->skewness;
	dest->num_real_roots = src->num_real_roots;

	dest->rpoly.degree = src->rpoly.degree;
	for (i = 0; i <= src->rpoly.degree; i++)
		mpz_set(dest->rpoly.coeff[i], src->rpoly.coeff[i]);

	dest->apoly.degree = src->apoly.degree;
	for (i = 0; i <= src->apoly.degree; i++)
		mpz_set(dest->apoly.coeff[i], src->apoly.coeff[i]);
}

/*------------------------------------------------------------------*/
static void poly_select_init(poly_select_t *s) {

	memset(s, 0, sizeof(poly_select_t));
	mpz_poly_init(&s->rpoly);
	mpz_poly_init(&s->apoly);
}

static void poly_select_free(poly_select_t *s) {

	mpz_poly_free(&s->rpoly);
	mpz_poly_free(&s->apoly);
}

/*------------------------------------------------------------------*/
void save_poly(poly_config_t *config, poly_select_t *poly) {

	/* save the polynomial if its combined score is the 
	   best seen so far. A really thorough implementation
	   would do this in two stages: the best X polynomials
	   would all be kept in stage 1, and in stage 2 the 
	   complete list (possibly combined from the output of
	   different machines) would get a more accurate 
	   root_score computed for each polynomial. The top
	   few such polynomials would then be subjected to 
	   sieving experiments */

	if (poly->combined_score > config->heap[0]->combined_score)
		poly_select_copy(config->heap[0], poly);
}

/*------------------------------------------------------------------*/
void analyze_poly(poly_config_t *config, poly_select_t *poly) {

	/* analyze a polynomial for sieving goodness
	  
	   The analysis routines are general enough so that
	   any polynomials can be tested, independent of
	   degree and skewness. The score assigned is
	   directly comparable to that of any other polynomials 
	   given to this routine */

	uint32 i;
	double root_score_r, root_score_a;
	ddpoly_t ddr, dda;
	mpz_poly_t *rpoly = &poly->rpoly;
	mpz_poly_t *apoly = &poly->apoly;

	poly->size_score = 0.0;
	poly->root_score = 0.0;
	poly->combined_score = 0.0;

	/* convert the polynomial coefficients from arbitrary
	   precision to dd_t floating point */

	ddr.degree = rpoly->degree;
	for (i = 0; i <= rpoly->degree; i++)
		ddr.coeff[i] = dd_gmp2dd(rpoly->coeff[i]);

	dda.degree = apoly->degree;
	for (i = 0; i <= apoly->degree; i++)
		dda.coeff[i] = dd_gmp2dd(apoly->coeff[i]);

	if (analyze_poly_roots(&poly->rpoly, PRIME_BOUND, &root_score_r))
		return;
	if (analyze_poly_roots(&poly->apoly, PRIME_BOUND, &root_score_a))
		return;

	if (analyze_poly_size(&config->integ_aux, 
				&ddr, &dda, &poly->size_score))
		return;

	get_bernstein_combined_score(poly, root_score_r + root_score_a);

	poly->root_score = root_score_a;

	analyze_poly_murphy(&config->integ_aux, &config->dickman_aux,
				&ddr, root_score_r, &dda, root_score_a, 
				poly->skewness, &poly->combined_score,
				&poly->num_real_roots);
}

/*------------------------------------------------------------------*/
void analyze_one_poly(msieve_obj *obj,
		      mpz_poly_t *rpoly, mpz_poly_t *apoly,
		      double skewness) {

	uint32 i;
	poly_config_t config;
	poly_select_t s;
	dd_precision_t prec = 0;
	uint32 prec_changed = 0;

	if (!dd_precision_is_ieee()) {
		prec_changed = 1;
		prec = dd_set_precision_ieee();
	}

	poly_config_init(&config);
	poly_select_init(&s);

	s.skewness = skewness;

	s.rpoly.degree = rpoly->degree;
	for (i = 0; i <= rpoly->degree; i++)
		mpz_set(s.rpoly.coeff[i], rpoly->coeff[i]);

	s.apoly.degree = apoly->degree;
	for (i = 0; i <= apoly->degree; i++)
		mpz_set(s.apoly.coeff[i], apoly->coeff[i]);

	analyze_poly(&config, &s);

	for (i = 0; i <= rpoly->degree; i++) {
		gmp_sprintf(obj->mp_sprintf_buf, 
				"R%u: %Zd\n", i, rpoly->coeff[i]);
		logprintf(obj, "%s", obj->mp_sprintf_buf);
	}
	for (i = 0; i <= apoly->degree; i++) {
		gmp_sprintf(obj->mp_sprintf_buf, 
				"A%u: %Zd\n", i, apoly->coeff[i]);
		logprintf(obj, "%s", obj->mp_sprintf_buf);
	}
	logprintf(obj, "skew %.2lf, size %.3e, "
			"alpha %.3lf, combined = %.3e rroots = %u\n",
			s.skewness, s.size_score, 
			s.root_score, s.combined_score,
			s.num_real_roots);
	poly_config_free(&config);
	poly_select_free(&s);

	if (prec_changed)
		dd_clear_precision(prec);
}

/*------------------------------------------------------------------*/
void poly_config_init(poly_config_t *config) {

	/* one-time initialization for polynomial search */

	uint32 i;

	config->heap_num_filled = 0;
	for (i = 0; i < POLY_HEAP_SIZE; i++) {
		config->heap[i] = (poly_select_t *)xcalloc((size_t)1, 
					sizeof(poly_select_t));
		poly_select_init(config->heap[i]);
	}

	integrate_init(&config->integ_aux, SIZE_EPS,
			double_exponential);

	dickman_init(&config->dickman_aux);
}

/*------------------------------------------------------------------*/
void poly_config_free(poly_config_t *config) {

	uint32 i;

	for (i = 0; i < POLY_HEAP_SIZE; i++) {
		poly_select_free(config->heap[i]);
		free(config->heap[i]);
	}
	integrate_free(&config->integ_aux);
	dickman_free(&config->dickman_aux);
}
