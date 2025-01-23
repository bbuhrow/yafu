/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: stage2.c 732 2012-08-04 02:32:46Z jasonp_sf $
--------------------------------------------------------------------*/

#include "stage2.h"

/*----------------------------------------------------------------------*/
void
assess_init(assess_t *a)
{
	integrate_init(&a->integ_aux, SIZE_EPS,
			double_exponential);
	dickman_init(&a->dickman_aux);
}

/*----------------------------------------------------------------------*/
void
assess_free(assess_t *a)
{
	integrate_free(&a->integ_aux);
	dickman_free(&a->dickman_aux);
}

/*-------------------------------------------------------------------------*/
static uint32
check_poly(curr_poly_t *c, mpz_t *coeffs, mpz_t lin0, 
		mpz_t gmp_N, uint32 degree) {

	uint32 i;

	mpz_set(c->gmp_help1, coeffs[degree]);
	mpz_set(c->gmp_help2, c->gmp_p);
	for (i = degree; i; i--) {
		mpz_mul(c->gmp_help1, c->gmp_help1, lin0);
		mpz_neg(c->gmp_help1, c->gmp_help1);
		mpz_addmul(c->gmp_help1, coeffs[i-1], c->gmp_help2);
		mpz_mul(c->gmp_help2, c->gmp_help2, c->gmp_p);
	}
	mpz_tdiv_r(c->gmp_help1, c->gmp_help1, gmp_N);
	if (mpz_cmp_ui(c->gmp_help1, (mp_limb_t)0) != 0) {
		printf("error: corrupt polynomial expand\n");
		return 0;
	}
	return 1;
}

/*-------------------------------------------------------------------------*/
#define MAX_CORRECT_STEPS 10

static int
pol_expand(curr_poly_t *c, mpz_t gmp_N, mpz_t high_coeff,
		mpz_t gmp_p, mpz_t gmp_d, 
		double coeff_bound, uint32 degree)
{
	uint32 i, j;

	if (mpz_cmp_ui(c->gmp_p, (mp_limb_t)1) == 0)
		mpz_set_ui(c->gmp_help1, (mp_limb_t)1);
	else {
		if (!mpz_invert(c->gmp_help1, gmp_d, gmp_p))
			return 0;
	}

	mpz_set(c->gmp_b[1], c->gmp_help1);
	for (i = 2; i < degree; i++)
		mpz_mul(c->gmp_b[i], c->gmp_b[i-1], c->gmp_help1);

	mpz_set(c->gmp_c[1], gmp_d);
	for (i = 2; i <= degree; i++)
		mpz_mul(c->gmp_c[i], c->gmp_c[i-1], gmp_d);

	mpz_set(c->gmp_a[degree], high_coeff);
	mpz_set(c->gmp_help2, gmp_N);

	for (i = degree - 1; (int32)i >= 0; i--) {

		mpz_mul(c->gmp_help3, c->gmp_a[i+1], c->gmp_c[i+1]);
		mpz_sub(c->gmp_help3, c->gmp_help2, c->gmp_help3);
		mpz_tdiv_q(c->gmp_help2, c->gmp_help3, gmp_p);

		if (i > 0) {
			mpz_tdiv_q(c->gmp_a[i], c->gmp_help2, c->gmp_c[i]);
			mpz_mul(c->gmp_help3, c->gmp_help2, c->gmp_b[i]);
			mpz_sub(c->gmp_help3, c->gmp_help3, c->gmp_a[i]);
			mpz_tdiv_r(c->gmp_help4, c->gmp_help3, gmp_p);

			if (mpz_sgn(c->gmp_help4) < 0)
				mpz_add(c->gmp_help4, c->gmp_help4, gmp_p);

			mpz_add(c->gmp_a[i], c->gmp_a[i], c->gmp_help4);
		}
	}
	mpz_set(c->gmp_a[0], c->gmp_help2);

	mpz_tdiv_q_2exp(c->gmp_help1, gmp_d, (mp_limb_t)1);
	for (i = 0; i < degree; i++) {
		for (j = 0; j < MAX_CORRECT_STEPS &&
			    mpz_cmpabs(c->gmp_a[i], c->gmp_help1) > 0; j++) {

			if (mpz_sgn(c->gmp_a[i]) < 0) {
				mpz_add(c->gmp_a[i], c->gmp_a[i], gmp_d);
				mpz_sub(c->gmp_a[i+1], c->gmp_a[i+1], gmp_p);
			}
			else {
				mpz_sub(c->gmp_a[i], c->gmp_a[i], gmp_d);
				mpz_add(c->gmp_a[i+1], c->gmp_a[i+1], gmp_p);
			}
		}

		if (j == MAX_CORRECT_STEPS)
			return 0;
	}

#if 0
	gmp_printf("%+Zd\n", c->gmp_lina[0]);
	gmp_printf("%+Zd\n", c->gmp_lina[1]);
	for (i = 0; i <= degree; i++)
		gmp_printf("%+Zd\n", c->gmp_a[i]);

	printf("coeff ratio = %.5lf\n",
		fabs(mpz_get_d(c->gmp_a[degree-2])) / coeff_bound);
#endif

	if (check_poly(c, c->gmp_a, 
			c->gmp_lina[0], gmp_N, degree) != 1) {
		return 0;
	}

	if (mpz_cmpabs_d(c->gmp_a[degree - 2], coeff_bound) > 0) {
		return 1;
	}
	return 2;
}

/*-------------------------------------------------------------------------*/
static void
curr_poly_init(curr_poly_t *c)
{
	int i;

	mpz_init(c->gmp_p);
	mpz_init(c->gmp_d);
	for (i = 0; i < 2; i++) {
		mpz_init(c->gmp_lina[i]);
		mpz_init(c->gmp_linb[i]);
		mpz_init(c->gmp_linc[i]);
	}
	for (i = 0; i < MAX_POLY_DEGREE + 1; i++) {
		mpz_init(c->gmp_a[i]);
		mpz_init(c->gmp_b[i]);
		mpz_init(c->gmp_c[i]);
	}
	mpz_init(c->gmp_help1);
	mpz_init(c->gmp_help2);
	mpz_init(c->gmp_help3);
	mpz_init(c->gmp_help4);
}

/*-------------------------------------------------------------------------*/
static void
curr_poly_free(curr_poly_t *c)
{
	int i;

	mpz_clear(c->gmp_p);
	mpz_clear(c->gmp_d);
	for (i = 0; i < 2; i++) {
		mpz_clear(c->gmp_lina[i]);
		mpz_clear(c->gmp_linb[i]);
		mpz_clear(c->gmp_linc[i]);
	}
	for (i = 0; i < MAX_POLY_DEGREE + 1; i++) {
		mpz_clear(c->gmp_a[i]);
		mpz_clear(c->gmp_b[i]);
		mpz_clear(c->gmp_c[i]);
	}
	mpz_clear(c->gmp_help1);
	mpz_clear(c->gmp_help2);
	mpz_clear(c->gmp_help3);
	mpz_clear(c->gmp_help4);
}

/*-------------------------------------------------------------------------*/
void
poly_sizeopt_init(poly_sizeopt_t *data, 
		 sizeopt_callback_t callback,
		 void *callback_data)
{
	memset(data, 0, sizeof(poly_sizeopt_t));
	mpz_init(data->gmp_N);

	data->internal = xmalloc(sizeof(curr_poly_t));
	curr_poly_init((curr_poly_t *)data->internal);

	data->callback = callback;
	data->callback_data = callback_data;
}

/*-------------------------------------------------------------------------*/
void
poly_sizeopt_free(poly_sizeopt_t *data)
{
	mpz_clear(data->gmp_N);
	curr_poly_free((curr_poly_t *)data->internal);
	free(data->internal);
}

/*-------------------------------------------------------------------------*/
void
poly_sizeopt_run(poly_sizeopt_t *data, mpz_t ad, mpz_t p, mpz_t d)
{
	double pol_norm;
	double alpha_proj;
	int status;
	curr_poly_t *c = (curr_poly_t *)(data->internal);

	mpz_set(c->gmp_d, d);
	mpz_set(c->gmp_p, p);
	mpz_neg(c->gmp_lina[0], c->gmp_d);
	mpz_set(c->gmp_lina[1], c->gmp_p);

	status = pol_expand(c, data->gmp_N, ad, p, d, 
			data->max_stage1_norm, data->degree);
	if (status != 2) {
		if (status == 0)
			fprintf(stderr, "expand failed\n");
		return;
	}

	optimize_initial(c, data->degree, &pol_norm, 0);

	stage2_root_score(data->degree, c->gmp_a, 100, &alpha_proj, 1);

	if (pol_norm * exp(alpha_proj) <= data->max_sizeopt_norm) {
		data->callback(data->degree, c->gmp_a, c->gmp_lina, 
				pol_norm, alpha_proj, 
				data->callback_data);
	}
}

/*-------------------------------------------------------------------------*/
void
poly_rootopt_init(poly_rootopt_t *data, 
		 msieve_obj *obj,
		 rootopt_callback_t callback,
		 void *callback_data)
{
	stage2_curr_data_t *s;

	memset(data, 0, sizeof(poly_rootopt_t));

	mpz_init(data->gmp_N);
	data->obj = obj;
	data->murphy_p_bound = PRIME_BOUND;
	data->callback = callback;
	data->callback_data = callback_data;

	s = (stage2_curr_data_t *)xmalloc(sizeof(stage2_curr_data_t));
	curr_poly_init(&s->curr_poly);
	root_sieve_init(&s->root_sieve);
	assess_init(&s->assess);
	data->internal = (void *)s;
}

/*-------------------------------------------------------------------------*/
void
poly_rootopt_free(poly_rootopt_t *data)
{
	stage2_curr_data_t *s = (stage2_curr_data_t *)(data->internal);

	mpz_clear(data->gmp_N);

	curr_poly_free(&s->curr_poly);
	root_sieve_free(&s->root_sieve);
	assess_free(&s->assess);

	free(data->internal);
}

/*-------------------------------------------------------------------------*/
void
poly_rootopt_run(poly_rootopt_t *data, mpz_t * alg_coeffs, 
		mpz_t * rat_coeffs, double sizeopt_norm, 
		double projective_alpha)
{
	uint32 i;
	stage2_curr_data_t *s = (stage2_curr_data_t *)(data->internal);
	curr_poly_t *c = &s->curr_poly;
	dd_precision_t precision = 0;
	uint32 precision_changed = 0;

	if (!dd_precision_is_ieee()) {
		precision_changed = 1;
		precision = dd_set_precision_ieee();
	}

	mpz_set(c->gmp_lina[0], rat_coeffs[0]);
	mpz_set(c->gmp_lina[1], rat_coeffs[1]);
	mpz_neg(c->gmp_d, rat_coeffs[0]);
	mpz_set(c->gmp_p, rat_coeffs[1]);

	for (i = 0; i <= data->degree; i++)
		mpz_set(c->gmp_a[i], alg_coeffs[i]);

	if (sizeopt_norm == 0) {

		/* size optimization possibly did not run
		   previously; check poly and get the norm */

		if (check_poly(c, c->gmp_a, c->gmp_lina[0], 
				data->gmp_N, data->degree) != 1) {
			goto finished;
		}

		optimize_initial(c, data->degree, &sizeopt_norm, 1);

		stage2_root_score(data->degree, c->gmp_a, 100, 
				&projective_alpha, 1);
	}

	if (sizeopt_norm * exp(projective_alpha) <= data->max_sizeopt_norm)
		root_sieve_run(data, sizeopt_norm, projective_alpha);

finished:
	if (precision_changed)
		dd_clear_precision(precision);
}
