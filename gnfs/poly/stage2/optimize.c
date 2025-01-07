/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: optimize.c 839 2013-01-13 14:31:02Z jasonp_sf $
--------------------------------------------------------------------*/

#include "stage2.h"

typedef double (*norm_t)(double *a, uint32 degree, double s);

/*-------------------------------------------------------------------------*/
static double
ifs_rectangular(double *a, uint32 degree, double s)
{
	double a0, a1, a2, a3, a4, a5, a6;
	double s2, s3, s4, s5, s6;
	double norm;

	if (s < 1)
		return 1e200;

	a0 = a[0];
	a1 = a[1] * s;
	s2 = s * s;
	a2 = a[2] * s2;
	s3 = s2 * s;
	a3 = a[3] * s3;
	s4 = s3 * s;
	a4 = a[4] * s4;

	switch (degree) {
	case 4:
		norm = 1.0 / 9.0 * (a4 * a4 + a0 * a0) + 
		       2.0 / 21.0 * (a4 * a2 + a2 * a0) + 
		       1.0 / 21.0 * (a3 * a3 + a1 * a1) + 
		       2.0 / 25.0 * (a4 * a0 + a3 * a1) + 
		       1.0 / 25.0 * a2 * a2;
		return norm / s4;

	case 5:
		s5 = s4 * s;
		a5 = a[5] * s5;
		norm = 1.0 / 11.0 * (a5 * a5 + a0 * a0) + 
		       2.0 / 27.0 * (a5 * a3 + a2 * a0) + 
		       1.0 / 27.0 * (a4 * a4 + a1 * a1) +
		       1.0 / 35.0 * (a3 * a3 + a2 * a2) +
		       2.0 / 35.0 * (a5 * a1 + a4 * a2 + a4 * a0 + a3 * a1);
		return norm / s5;

	case 6:
		s5 = s4 * s;
		a5 = a[5] * s5;
		s6 = s5 * s;
		a6 = a[6] * s6;
		norm = 1.0 / 13.0 * (a6 * a6 + a0 * a0) + 
		       2.0 / 33.0 * (a6 * a4 + a2 * a0) + 
		       1.0 / 33.0 * (a5 * a5 + a1 * a1) +
		       1.0 / 45.0 * (a4 * a4 + a2 * a2) + 
		       2.0 / 45.0 * (a6 * a2 + a5 * a3 + a4 * a0 + a3 * a1) + 
		       2.0 / 49.0 * (a6 * a0 + a5 * a1 + a4 * a2) + 
		       1.0 / 49.0 * a3 * a3;
		return norm / s6;
	}

	return 1e200;
}

static double
ifs_radial(double *a, uint32 degree, double s)
{
	double a0, a1, a2, a3, a4, a5, a6;
	double s2, s3, s4, s5, s6;
	double norm;

	if (s < 1)
		return 1e200;

	a0 = a[0];
	a1 = a[1] * s;
	s2 = s * s;
	a2 = a[2] * s2;
	s3 = s2 * s;
	a3 = a[3] * s3;
	s4 = s3 * s;
	a4 = a[4] * s4;

	switch (degree) {
	case 4:
		norm = 35.0 * (a4 * a4 + a0 * a0) + 
		       10.0 * (a4 * a2 + a2 * a0) + 
		        5.0 * (a3 * a3 + a1 * a1) + 
		        6.0 * (a4 * a0 + a3 * a1) + 
		        3.0 * a2 * a2;
		return norm / s4;

	case 5:
		s5 = s4 * s;
		a5 = a[5] * s5;
		norm = 63.0 * (a5 * a5 + a0 * a0) + 
		       14.0 * (a5 * a3 + a2 * a0) + 
		        7.0 * (a4 * a4 + a1 * a1) +
		        3.0 * (a3 * a3 + a2 * a2) +
		        6.0 * (a5 * a1 + a4 * a2 + a4 * a0 + a3 * a1);
		return norm / s5;

	case 6:
		s5 = s4 * s;
		a5 = a[5] * s5;
		s6 = s5 * s;
		a6 = a[6] * s6;
		norm = 231.0 * (a6 * a6 + a0 * a0) + 
		        42.0 * (a6 * a4 + a2 * a0) + 
		        21.0 * (a5 * a5 + a1 * a1) +
		         7.0 * (a4 * a4 + a2 * a2) + 
		        14.0 * (a6 * a2 + a5 * a3 + a4 * a0 + a3 * a1) + 
		        10.0 * (a6 * a0 + a5 * a1 + a4 * a2) + 
		         5.0 * a3 * a3;
		return norm / s6;
	}

	return 1e200;
}
 
/*----------------------------------------------------------------------*/
uint32
stage2_root_score(uint32 deg1, mpz_t *coeff1, 
		 uint32 prime_bound, double *score,
		 uint32 projective_only)
{
	uint32 i;
	uint32 status;
	mpz_poly_t apoly;

	mpz_poly_init(&apoly);

	apoly.degree = deg1;
	for (i = 0; i <= deg1; i++)
		mpz_set(apoly.coeff[i], coeff1[i]);

	if (projective_only)
		status = analyze_poly_roots_projective(&apoly, 
						prime_bound, score);
	else
		status = analyze_poly_roots(&apoly, prime_bound, score);

	mpz_poly_free(&apoly);
	return status;
}

static const double xlate_weights[6+1][6+1] = {
	{  1.,  1.,  1.,  1.,  1.,  1.,  1.},
	{  0.,  1.,  2.,  3.,  4.,  5.,  6.},
	{  0.,  0.,  1.,  3.,  6., 10., 15.},
	{  0.,  0.,  0.,  1.,  4., 10., 20.},
	{  0.,  0.,  0.,  0.,  1.,  5., 15.},
	{  0.,  0.,  0.,  0.,  0.,  1.,  6.},
	{  0.,  0.,  0.,  0.,  0.,  0.,  1.},
};

/*-------------------------------------------------------------------------*/
static void
translate_d(double *c, double *d, uint32 deg, double x)
{
	uint32 i, j;

	for (i = 0, x = -x; i < deg; i++) {

		const double *weights = xlate_weights[i];
		double accum = d[deg] * weights[deg];

		for (j = deg; j > i; j--) {
			accum = accum * x + d[j-1] * weights[j-1];
		}
		c[i] = accum;
	}
	c[i] = d[i];
}

/*-------------------------------------------------------------------------*/
static void
translate_dd(dd_t *c, dd_t *d, uint32 deg, double x)
{
	uint32 i, j;

	for (i = 0, x = -x; i < deg; i++) {

		const double *weights = xlate_weights[i];
		dd_t accum = dd_mul_d(d[deg], weights[deg]);

		for (j = deg; j > i; j--) {
			accum = dd_add_dd(dd_mul_d(accum, x),
					  dd_mul_d(d[j-1], weights[j-1]));
		}
		c[i] = accum;
	}
	c[i] = d[i];
}

/*-------------------------------------------------------------------------*/
static void
translate_gmp(curr_poly_t *c, mpz_t *gmp_c, uint32 deg,
		mpz_t *lin, int64 k)
{
	/* note that unlike the other translation functions, 
	   gmp_c[] and lin[] are overwritten */

	uint32 i, j;

	int64_2gmp(-k, c->gmp_help1);

	for (i = 0; i < deg; i++) {

		const double *weights = xlate_weights[i];

		mpz_set_d(c->gmp_help2, weights[deg]);
		mpz_mul(c->gmp_help2, c->gmp_help2, gmp_c[deg]);

		for (j = deg; j > i; j--) {
			mpz_set_d(c->gmp_help3, weights[j-1]);
			mpz_mul(c->gmp_help3, c->gmp_help3, gmp_c[j-1]);
			mpz_addmul(c->gmp_help3, c->gmp_help2, c->gmp_help1);
			mpz_swap(c->gmp_help2, c->gmp_help3);
		}
		mpz_swap(gmp_c[i], c->gmp_help2);
	}

	mpz_addmul(lin[0], lin[1], c->gmp_help1);
}

/*-------------------------------------------------------------------------*/
#define SKEWNESS 0
#define TRANSLATE_SIZE 1
#define ROTATE0 2
#define ROTATE1 3
#define ROTATE2 4

typedef struct {
	uint32 rotate_dim;
	dpoly_t *drpoly;
	dpoly_t *dapoly;

	ddpoly_t *rpoly;
	ddpoly_t *apoly;
	integrate_t *integ_aux;

	dickman_t *dickman_aux;
	double root_score_r;
	double root_score_a;

	uint32 num_real_roots;
	norm_t norm_callback;
} opt_data_t;

static double poly_xlate_callback(double *v, void *extra)
{
	opt_data_t *opt = (opt_data_t *)extra;
	dpoly_t *apoly = opt->dapoly;
	double translated[MAX_POLY_DEGREE + 1];
	double t = floor(v[TRANSLATE_SIZE] + 0.5);
	double s = v[SKEWNESS];

	translate_d(translated, apoly->coeff, apoly->degree, t);

	return ifs_rectangular(translated, apoly->degree, s);
}

static double poly_skew_callback(double *v, void *extra)
{
	opt_data_t *opt = (opt_data_t *)extra;
	dpoly_t *apoly = opt->dapoly;
	double s = v[SKEWNESS];

	return opt->norm_callback(apoly->coeff, apoly->degree, s);
}

static double poly_rotate_callback(double *v, void *extra)
{
	uint32 i;
	opt_data_t *opt = (opt_data_t *)extra;
	dpoly_t apoly = *(opt->dapoly);
	double translated[MAX_POLY_DEGREE + 1];
	double s = floor(v[SKEWNESS] + 0.5);
	double t = floor(v[TRANSLATE_SIZE] + 0.5);
	double r0 = opt->drpoly->coeff[0];
	double r1 = opt->drpoly->coeff[1];

	if (s < 1.0)
		return 1e200;

	for (i = 0; i <= opt->rotate_dim; i++) {
		double c = floor(v[ROTATE0 + i] + 0.5);
		apoly.coeff[i] += r0 * c;
		apoly.coeff[i+1] += r1 * c;
	}
	
	translate_d(translated, apoly.coeff, apoly.degree, t);
	return opt->norm_callback(translated, apoly.degree, s);
}

static double poly_murphy_callback(double *v, void *extra)
{
	opt_data_t *opt = (opt_data_t *)extra;
	ddpoly_t *apoly = opt->apoly;
	ddpoly_t *rpoly = opt->rpoly;
	ddpoly_t new_rpoly, new_apoly;
	double t = floor(v[TRANSLATE_SIZE] + 0.5);
	double s = v[SKEWNESS];
	double score = 0;

	if (s < 0.5)
		return 1e50;

	new_rpoly.degree = 1;
	new_rpoly.coeff[1] = rpoly->coeff[1];
	new_rpoly.coeff[0] = dd_sub_dd(rpoly->coeff[0],
					dd_mul_d(rpoly->coeff[1], t));

	new_apoly.degree = apoly->degree;
	translate_dd(new_apoly.coeff, apoly->coeff, apoly->degree, t);

	analyze_poly_murphy(opt->integ_aux, opt->dickman_aux,
				&new_rpoly, opt->root_score_r,
				&new_apoly, opt->root_score_a,
				s, &score, &opt->num_real_roots);

	return -score;
}

/*-------------------------------------------------------------------------*/
void
optimize_initial(curr_poly_t *c, uint32 deg, double *pol_norm, uint32 skew_only)
{
	uint32 rotate_dim = deg - 4;
	uint32 num_vars = rotate_dim + 3;
	opt_data_t opt_data;
	uint32 i, j;
	double best[MAX_VARS];
	double score, last_score, tol;
	dpoly_t rpoly, apoly;
	objective_func objective = poly_rotate_callback;

	opt_data.rotate_dim = rotate_dim;
	opt_data.drpoly = &rpoly;
	opt_data.dapoly = &apoly;
	opt_data.norm_callback = ifs_radial;

	best[SKEWNESS] = 1000;
	best[TRANSLATE_SIZE] = 0;
	best[ROTATE0] = 0;
	best[ROTATE1] = 0;
	best[ROTATE2] = 0;
	if (skew_only) {
		num_vars = 1;
		objective = poly_skew_callback;
	}
	else if (deg == 6) {
		score = optimize_initial_deg6(best, c, deg);
	}

	score = 1e200;
	tol = 1e-5;
	rpoly.degree = 1;
	for (i = 0; i <= 1; i++)
		rpoly.coeff[i] = mpz_get_d(c->gmp_lina[i]);
	apoly.degree = deg;
	for (i = 0; i <= deg; i++)
		apoly.coeff[i] = mpz_get_d(c->gmp_a[i]);

	for (i = 0; i < 2; i++) {

		uint32 num_minimize = 0;

		do {
			last_score = score;
			score = minimize(best, num_vars, tol, 40, 
					objective, &opt_data);

			for (j = 0; j <= rotate_dim; j++) {
				double cj = floor(best[ROTATE0 + j] + 0.5);
				mpz_set_d(c->gmp_help1, cj);
				mpz_addmul(c->gmp_a[j+1], c->gmp_help1, 
						c->gmp_p);
				mpz_submul(c->gmp_a[j], c->gmp_help1, 
						c->gmp_d);
			}
			translate_gmp(c, c->gmp_a, deg, c->gmp_lina,
					(int64)(best[TRANSLATE_SIZE] + 0.5));

			mpz_neg(c->gmp_d, c->gmp_lina[0]);
			for (j = 0; j <= 1; j++)
				rpoly.coeff[j] = mpz_get_d(c->gmp_lina[j]);
			for (j = 0; j <= deg; j++)
				apoly.coeff[j] = mpz_get_d(c->gmp_a[j]);
			best[TRANSLATE_SIZE] = 0;
			best[ROTATE0] = 0;
			best[ROTATE1] = 0;
			best[ROTATE2] = 0;

			if (++num_minimize >= 5)
				break;

		} while (fabs(score - last_score) > .001 * fabs(score));

		if (i == 0) {
			opt_data.norm_callback = ifs_rectangular;
			tol = 1e-5;
			score = ifs_rectangular(apoly.coeff, apoly.degree,
						best[SKEWNESS]);
		}
	}

	*pol_norm = sqrt(fabs(score));
#if 0
	printf("norm %.7e skew %lf\n", *pol_norm, best[SKEWNESS]);
	for (i = 0; i < 2; i++)
		gmp_printf("%+Zd\n", c->gmp_lina[i]);
	for (i = 0; i <= deg; i++)
		gmp_printf("%+Zd\n", c->gmp_a[i]);
#endif
}

/*-------------------------------------------------------------------------*/
double
optimize_basic(dpoly_t *apoly, double *best_skewness,
		double *best_translation)
{
	opt_data_t opt_data;
	double best[MAX_VARS];
	double score;

	opt_data.dapoly = apoly;
	best[TRANSLATE_SIZE] = 0;
	best[SKEWNESS] = 1000;

	score = minimize(best, 2, 1e-5, 40, poly_xlate_callback, &opt_data);

	*best_translation = floor(best[TRANSLATE_SIZE] + 0.5);
	*best_skewness  = best[SKEWNESS];
	return sqrt(score);
}

/*-------------------------------------------------------------------------*/
static void
optimize_final_core(curr_poly_t *c, assess_t *assess, uint32 deg,
			double root_score_r, double root_score_a,
			double *best_score_out,
			double *best_skewness_out, 
			uint32 *num_real_roots_out)
{
	uint32 i;
	opt_data_t opt_data;
	double best[MAX_VARS];
	double score;
	ddpoly_t rpoly, apoly;

	opt_data.rpoly = &rpoly;
	opt_data.apoly = &apoly;
	opt_data.integ_aux = &assess->integ_aux;
	opt_data.dickman_aux = &assess->dickman_aux;
	opt_data.root_score_r = root_score_r;
	opt_data.root_score_a = root_score_a;
	rpoly.degree = 1;
	apoly.degree = deg;

	best[TRANSLATE_SIZE] = 0;
	best[SKEWNESS] = 1000;

	for (i = 0; i <= 1; i++)
		rpoly.coeff[i] = dd_gmp2dd(c->gmp_linb[i]);
	for (i = 0; i <= deg; i++)
		apoly.coeff[i] = dd_gmp2dd(c->gmp_b[i]);

	score = minimize(best, 2, 1e-5, 40, 
			poly_murphy_callback, &opt_data);

	translate_gmp(c, c->gmp_b, deg, c->gmp_linb, 
			(int64)(best[TRANSLATE_SIZE] + 0.5));

	*best_score_out = fabs(score);
	*best_skewness_out = best[SKEWNESS];
	*num_real_roots_out = opt_data.num_real_roots;
}

/*-------------------------------------------------------------------------*/
static void
get_bernstein_score(curr_poly_t *c, assess_t *assess,
		uint32 deg, double root_score, double *eptr)
{
	uint32 i;
	double size_score;
	ddpoly_t rpoly, apoly;

	root_score = pow(exp(root_score), -2./(deg + 1));

	rpoly.degree = 1;
	for (i = 0; i <= 1; i++)
		rpoly.coeff[i] = dd_gmp2dd(c->gmp_linb[i]);

	apoly.degree = deg;
	for (i = 0; i <= deg; i++)
		apoly.coeff[i] = dd_gmp2dd(c->gmp_b[i]);

	analyze_poly_size(&assess->integ_aux, 
			&rpoly, &apoly, &size_score);
	if (size_score == 0.0)
		printf("error: size score computation failed\n");

	*eptr = size_score * root_score;
}

/*-------------------------------------------------------------------------*/
void
optimize_final(mpz_t x, mpz_t y, int64 z, poly_rootopt_t *data)
{
	uint32 i;
	uint32 deg = data->degree;
	uint32 num_real_roots;
	double alpha_r, alpha_a, skewness, bscore, combined_score;
	stage2_curr_data_t *s = (stage2_curr_data_t *)data->internal;
	curr_poly_t *c = &s->curr_poly;
	assess_t *assess = &s->assess;

	for (i = 0; i <= 1; i++)
		mpz_set(c->gmp_linb[i], c->gmp_lina[i]);

	for (i = 0; i <= deg; i++)
		mpz_set(c->gmp_b[i], c->gmp_a[i]);

	int64_2gmp(z, c->gmp_help1);
	mpz_addmul(c->gmp_b[3], c->gmp_p, c->gmp_help1);
	mpz_submul(c->gmp_b[2], c->gmp_d, c->gmp_help1);

	mpz_addmul(c->gmp_b[2], c->gmp_p, y);
	mpz_submul(c->gmp_b[1], c->gmp_d, y);

	mpz_addmul(c->gmp_b[1], c->gmp_p, x);
	mpz_submul(c->gmp_b[0], c->gmp_d, x);

	if (stage2_root_score(deg, c->gmp_b, 
			data->murphy_p_bound, &alpha_a, 0))
		return;

	if (stage2_root_score(1, c->gmp_linb, 
			data->murphy_p_bound, &alpha_r, 0))
		return;

	if (alpha_a > -4.5)
		return;

	get_bernstein_score(c, assess, deg, alpha_a, &bscore);

	if (bscore > data->min_e_bernstein) {

		optimize_final_core(c, assess, deg, alpha_r, alpha_a, 
				&combined_score, &skewness,
				&num_real_roots);

		if (combined_score > data->min_e) {
			data->callback(data->callback_data, deg, c->gmp_b, 
					c->gmp_linb, skewness, bscore,
					alpha_a, combined_score, 
					num_real_roots);
		}
	}
}
