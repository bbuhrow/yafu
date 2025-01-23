/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id$
--------------------------------------------------------------------*/

#include "stage2.h"

#define MAX_VARS 5
#define MAX_PARAMS 9
#define MAX_VARDEGREE 15
#define VAR_S 0
#define VAR_T 1
#define VAR_C0 2
#define VAR_C1 3
#define VAR_C2 4

typedef struct {
	double coeff;
	int32 params[MAX_PARAMS];
	int32 powers[MAX_VARS];
} term_t;

typedef struct {
	uint32 num_terms;
	uint32 num_terms_alloc;
	term_t *terms;
} multi_poly_t;

typedef struct {
	uint32 num_vars;
	multi_poly_t *objective;
	multi_poly_t grad[MAX_VARS];
	multi_poly_t hess[MAX_VARS][MAX_VARS];
	curr_poly_t *c;
	int32 max_powers[MAX_VARS];
	double params[MAX_PARAMS];
} opt_data_t;

static char *objective_deg6 = 
"231*a6^2*t^12/s^6-462*a5*a6*t^11/s^6+1386*a6^2*t^10/s^4+462*a4*a6*t^10/s^6+"
"231*a5^2*t^10/s^6-2310*a5*a6*t^9/s^4-462*a6*c2*r1*t^9/s^6-462*a3*a6*t^9/s^6-"
"462*a4*a5*t^9/s^6+3465*a6^2*t^8/s^2+1890*a4*a6*t^8/s^4+945*a5^2*t^8/s^4+"
"462*a5*c2*r1*t^8/s^6+462*a6*c1*r1*t^8/s^6+462*a6*c2*r0*t^8/s^6+"
"462*a2*a6*t^8/s^6+462*a3*a5*t^8/s^6+231*a4^2*t^8/s^6-4620*a5*a6*t^7/s^2-"
"1512*a6*c2*r1*t^7/s^4-1512*a3*a6*t^7/s^4-1512*a4*a5*t^7/s^4-"
"462*a4*c2*r1*t^7/s^6-462*a5*c1*r1*t^7/s^6-462*a6*c0*r1*t^7/s^6-"
"462*a5*c2*r0*t^7/s^6-462*a6*c1*r0*t^7/s^6-462*a1*a6*t^7/s^6-"
"462*a2*a5*t^7/s^6-462*a3*a4*t^7/s^6+2940*a4*a6*t^6/s^2+1470*a5^2*t^6/s^2+"
"1176*a5*c2*r1*t^6/s^4+1176*a6*c1*r1*t^6/s^4+1176*a6*c2*r0*t^6/s^4+"
"1176*a2*a6*t^6/s^4+1176*a3*a5*t^6/s^4+588*a4^2*t^6/s^4+231*c2^2*r1^2*t^6/s^6+"
"462*a3*c2*r1*t^6/s^6+462*a4*c1*r1*t^6/s^6+462*a5*c0*r1*t^6/s^6+"
"462*a4*c2*r0*t^6/s^6+462*a5*c1*r0*t^6/s^6+462*a6*c0*r0*t^6/s^6+"
"462*a0*a6*t^6/s^6+462*a1*a5*t^6/s^6+462*a2*a4*t^6/s^6+231*a3^2*t^6/s^6+"
"4620*a6^2*t^6-1764*a6*c2*r1*t^5/s^2-1764*a3*a6*t^5/s^2-1764*a4*a5*t^5/s^2-"
"882*a4*c2*r1*t^5/s^4-882*a5*c1*r1*t^5/s^4-882*a6*c0*r1*t^5/s^4-"
"882*a5*c2*r0*t^5/s^4-882*a6*c1*r0*t^5/s^4-882*a1*a6*t^5/s^4-882*a2*a5*t^5/s^4-"
"882*a3*a4*t^5/s^4-462*c1*c2*r1^2*t^5/s^6-462*c2^2*r0*r1*t^5/s^6-"
"462*a2*c2*r1*t^5/s^6-462*a3*c1*r1*t^5/s^6-462*a4*c0*r1*t^5/s^6-"
"462*a3*c2*r0*t^5/s^6-462*a4*c1*r0*t^5/s^6-462*a5*c0*r0*t^5/s^6-"
"462*a0*a5*t^5/s^6-462*a1*a4*t^5/s^6-462*a2*a3*t^5/s^6-4620*a5*a6*t^5+"
"3465*a6^2*s^2*t^4+980*a5*c2*r1*t^4/s^2+980*a6*c1*r1*t^4/s^2+"
"980*a6*c2*r0*t^4/s^2+980*a2*a6*t^4/s^2+980*a3*a5*t^4/s^2+490*a4^2*t^4/s^2+"
"315*c2^2*r1^2*t^4/s^4+630*a3*c2*r1*t^4/s^4+630*a4*c1*r1*t^4/s^4+"
"630*a5*c0*r1*t^4/s^4+630*a4*c2*r0*t^4/s^4+630*a5*c1*r0*t^4/s^4+"
"630*a6*c0*r0*t^4/s^4+630*a0*a6*t^4/s^4+630*a1*a5*t^4/s^4+630*a2*a4*t^4/s^4+"
"315*a3^2*t^4/s^4+462*c0*c2*r1^2*t^4/s^6+231*c1^2*r1^2*t^4/s^6+"
"924*c1*c2*r0*r1*t^4/s^6+462*a1*c2*r1*t^4/s^6+462*a2*c1*r1*t^4/s^6+"
"462*a3*c0*r1*t^4/s^6+231*c2^2*r0^2*t^4/s^6+462*a2*c2*r0*t^4/s^6"
"+462*a3*c1*r0*t^4/s^6+462*a4*c0*r0*t^4/s^6+462*a0*a4*t^4/s^6+"
"462*a1*a3*t^4/s^6+231*a2^2*t^4/s^6+2100*a4*a6*t^4+1050*a5^2*t^4-"
"2310*a5*a6*s^2*t^3-490*a4*c2*r1*t^3/s^2-490*a5*c1*r1*t^3/s^2-"
"490*a6*c0*r1*t^3/s^2-490*a5*c2*r0*t^3/s^2-490*a6*c1*r0*t^3/s^2-"
"490*a1*a6*t^3/s^2-490*a2*a5*t^3/s^2-490*a3*a4*t^3/s^2-420*c1*c2*r1^2*t^3/s^4-"
"420*c2^2*r0*r1*t^3/s^4-420*a2*c2*r1*t^3/s^4-420*a3*c1*r1*t^3/s^4-"
"420*a4*c0*r1*t^3/s^4-420*a3*c2*r0*t^3/s^4-420*a4*c1*r0*t^3/s^4-"
"420*a5*c0*r0*t^3/s^4-420*a0*a5*t^3/s^4-420*a1*a4*t^3/s^4-420*a2*a3*t^3/s^4-"
"462*c0*c1*r1^2*t^3/s^6-924*c0*c2*r0*r1*t^3/s^6-462*c1^2*r0*r1*t^3/s^6"
"-462*a0*c2*r1*t^3/s^6-462*a1*c1*r1*t^3/s^6-462*a2*c0*r1*t^3/s^6-"
"462*c1*c2*r0^2*t^3/s^6-462*a1*c2*r0*t^3/s^6-462*a2*c1*r0*t^3/s^6-"
"462*a3*c0*r0*t^3/s^6-462*a0*a3*t^3/s^6-462*a1*a2*t^3/s^6-840*a6*c2*r1*t^3-"
"840*a3*a6*t^3-840*a4*a5*t^3+1386*a6^2*s^4*t^2+630*a4*a6*s^2*t^2+"
"315*a5^2*s^2*t^2+105*c2^2*r1^2*t^2/s^2+210*a3*c2*r1*t^2/s^2+"
"210*a4*c1*r1*t^2/s^2+210*a5*c0*r1*t^2/s^2+210*a4*c2*r0*t^2/s^2+"
"210*a5*c1*r0*t^2/s^2+210*a6*c0*r0*t^2/s^2+210*a0*a6*t^2/s^2+"
"210*a1*a5*t^2/s^2+210*a2*a4*t^2/s^2+105*a3^2*t^2/s^2+252*c0*c2*r1^2*t^2/s^4+"
"126*c1^2*r1^2*t^2/s^4+504*c1*c2*r0*r1*t^2/s^4+252*a1*c2*r1*t^2/s^4+"
"252*a2*c1*r1*t^2/s^4+252*a3*c0*r1*t^2/s^4+126*c2^2*r0^2*t^2/s^4"
"+252*a2*c2*r0*t^2/s^4+252*a3*c1*r0*t^2/s^4+252*a4*c0*r0*t^2/s^4+"
"252*a0*a4*t^2/s^4+252*a1*a3*t^2/s^4+126*a2^2*t^2/s^4+231*c0^2*r1^2*t^2/s^6+"
"924*c0*c1*r0*r1*t^2/s^6+462*a0*c1*r1*t^2/s^6+462*a1*c0*r1*t^2/s^6+"
"462*c0*c2*r0^2*t^2/s^6+231*c1^2*r0^2*t^2/s^6+462*a0*c2*r0*t^2/s^6+"
"462*a1*c1*r0*t^2/s^6+462*a2*c0*r0*t^2/s^6+462*a0*a2*t^2/s^6+231*a1^2*t^2/s^6"
"+280*a5*c2*r1*t^2+280*a6*c1*r1*t^2+280*a6*c2*r0*t^2+280*a2*a6*t^2+"
"280*a3*a5*t^2+140*a4^2*t^2-462*a5*a6*s^4*t-126*a6*c2*r1*s^2*t-"
"126*a3*a6*s^2*t-126*a4*a5*s^2*t-70*c1*c2*r1^2*t/s^2-70*c2^2*r0*r1*t/s^2-"
"70*a2*c2*r1*t/s^2-70*a3*c1*r1*t/s^2-70*a4*c0*r1*t/s^2-70*a3*c2*r0*t/s^2-"
"70*a4*c1*r0*t/s^2-70*a5*c0*r0*t/s^2-70*a0*a5*t/s^2-70*a1*a4*t/s^2-"
"70*a2*a3*t/s^2-126*c0*c1*r1^2*t/s^4-252*c0*c2*r0*r1*t/s^4-"
"126*c1^2*r0*r1*t/s^4-126*a0*c2*r1*t/s^4-126*a1*c1*r1*t/s^4-"
"126*a2*c0*r1*t/s^4-126*c1*c2*r0^2*t/s^4-126*a1*c2*r0*t/s^4-"
"126*a2*c1*r0*t/s^4-126*a3*c0*r0*t/s^4-126*a0*a3*t/s^4-126*a1*a2*t/s^4-"
"462*c0^2*r0*r1*t/s^6-462*a0*c0*r1*t/s^6-462*c0*c1*r0^2*t/s^6-"
"462*a0*c1*r0*t/s^6-462*a1*c0*r0*t/s^6-462*a0*a1*t/s^6-70*a4*c2*r1*t-"
"70*a5*c1*r1*t-70*a6*c0*r1*t-70*a5*c2*r0*t-70*a6*c1*r0*t-70*a1*a6*t-"
"70*a2*a5*t-70*a3*a4*t+231*a6^2*s^6+42*a4*a6*s^4+21*a5^2*s^4+"
"14*a5*c2*r1*s^2+14*a6*c1*r1*s^2+14*a6*c2*r0*s^2+14*a2*a6*s^2+14*a3*a5*s^2+"
"7*a4^2*s^2+14*c0*c2*r1^2/s^2+7*c1^2*r1^2/s^2+28*c1*c2*r0*r1/s^2+"
"14*a1*c2*r1/s^2+14*a2*c1*r1/s^2+14*a3*c0*r1/s^2+7*c2^2*r0^2/s^2+"
"14*a2*c2*r0/s^2+14*a3*c1*r0/s^2+14*a4*c0*r0/s^2+14*a0*a4/s^2+14*a1*a3/s^2"
"+7*a2^2/s^2+21*c0^2*r1^2/s^4+84*c0*c1*r0*r1/s^4+42*a0*c1*r1/s^4+"
"42*a1*c0*r1/s^4+42*c0*c2*r0^2/s^4+21*c1^2*r0^2/s^4+42*a0*c2*r0/s^4+"
"42*a1*c1*r0/s^4+42*a2*c0*r0/s^4+42*a0*a2/s^4+21*a1^2/s^4+231*c0^2*r0^2/s^6+"
"462*a0*c0*r0/s^6+231*a0^2/s^6+5*c2^2*r1^2+10*a3*c2*r1+10*a4*c1*r1+"
"10*a5*c0*r1+10*a4*c2*r0+10*a5*c1*r0+10*a6*c0*r0+10*a0*a6+10*a1*a5+"
"10*a2*a4+5*a3^2";

/*--------------------------------------------------------------------*/
static void
poly_alloc(multi_poly_t *poly, uint32 num_terms)
{
	poly->terms = (term_t *)xcalloc(num_terms, sizeof(term_t));
	poly->num_terms_alloc = num_terms;
	poly->num_terms = 0;
}

/*--------------------------------------------------------------------*/
static void
poly_free(multi_poly_t *poly)
{
	free(poly->terms);
}

/*--------------------------------------------------------------------*/
static void 
poly_find_max_powers(multi_poly_t *poly, 
		int32 max_powers[MAX_VARS])
{
	uint32 i, j;

	for (i = 0; i < poly->num_terms; i++) {

		term_t *curr_term = poly->terms + i;

		for (j = 0; j < MAX_VARS; j++) {
			max_powers[j] = MAX(max_powers[j], 
					   abs(curr_term->powers[j]));
		}
	}
}

/*--------------------------------------------------------------------*/
static void
poly_diff(multi_poly_t *poly, uint32 var, multi_poly_t *res)
{
	uint32 i, j;

	for (i = j = 0; i < poly->num_terms; i++) {
		term_t *curr_term = poly->terms + i;
		term_t *diff_term = res->terms + j;
		int32 curr_pow = curr_term->powers[var];

		if (curr_pow != 0) {
			*diff_term = *curr_term;
			diff_term->coeff *= curr_pow;
			diff_term->powers[var]--;
			j++;
		}
	}

	res->num_terms = j;
}

/*--------------------------------------------------------------------*/
static double 
poly_eval(multi_poly_t *p,
		double params[MAX_PARAMS],
		double powers[MAX_VARS][MAX_VARDEGREE + 1])
{
	uint32 i, j, k;
	uint32 num_terms = p->num_terms;
	double f = 0;
	double t;

	for (i = 0; i < num_terms; i++) {
		term_t *term = p->terms + i;

		t = term->coeff;
		for (j = 0; j < MAX_VARS; j++) {
			int32 pow = term->powers[j];
			if (pow > 0)
				t *= powers[j][pow];
			else if (pow < 0)
				t /= powers[j][abs(pow)];
		}

		for (j = 0; j < MAX_PARAMS; j++) {
			int32 pow = term->params[j];
			if (pow > 0) {
				for (k = 0; k < pow; k++)
					t *= params[j];
			}
			else if (pow < 0) {
				for (k = 0; k < abs(pow); k++)
					t /= params[j];
			}
		}
		f += t;
	}

	return f;
}

/*--------------------------------------------------------------------*/
static void
parse_objective(multi_poly_t *poly, char *objective)
{
	uint32 i;
	uint32 num_terms = 1;
	char *tmp;

	tmp = objective;
	while (*tmp != 0) {
		if (*tmp == '+' || *tmp == '-')
			num_terms++;
		tmp++;
	}

	poly_alloc(poly, num_terms);

	tmp = objective;
	for (i = 0; i < num_terms; i++) {

		char *tmp1;
		int32 coeff;
		uint32 is_neg;
		uint32 power;
		term_t *curr_term = poly->terms + i;

		memset(curr_term, 0, sizeof(term_t));
		curr_term->coeff = 1;

		coeff = strtol(tmp, &tmp1, 10);
		if (coeff != 0)
			curr_term->coeff *= coeff;

		tmp = tmp1;
		is_neg = 0;
		power = 1;
		while (*tmp != 0 && *tmp != '+' && *tmp != '-') {

			uint32 which_var;

			switch (*tmp) {
			case '/':
				is_neg = 1;
				power = 1;
				tmp++;
				break;

			case '*':
				is_neg = 0;
				power = 1;
				tmp++;
				break;

			case 'a':
			case 'r':
				tmp1 = NULL;
				which_var = tmp[1] - '0';
				if (tmp[2] == '^')
					power = strtol(tmp + 3, &tmp1, 10);

				if (*tmp == 'r') {
					curr_term->params[which_var] = 
						(is_neg ? -1 : 1) * power;
				}
				else {
					curr_term->params[which_var+2] =
						(is_neg ? -1 : 1) * power;
				}

				tmp += 2;
				if (tmp1 != NULL)
					tmp = tmp1;
				break;

			case 'c':
				tmp1 = NULL;
				which_var = tmp[1] - '0';
				if (tmp[2] == '^')
					power = strtol(tmp + 3, &tmp1, 10);

				curr_term->powers[VAR_C0 + which_var] +=
					(is_neg ? -1 : +1) * power;

				tmp += 2;
				if (tmp1 != NULL)
					tmp = tmp1;
				break;

			case 's':
				tmp1 = NULL;
				if (tmp[1] == '^')
					power = strtol(tmp + 2, &tmp1, 10);

				curr_term->powers[VAR_S] +=
					(is_neg ? -1 : +1) * power;

				tmp++;
				if (tmp1 != NULL)
					tmp = tmp1;
				break;

			case 't':
				tmp1 = NULL;
				if (tmp[1] == '^')
					power = strtol(tmp + 2, &tmp1, 10);

				curr_term->powers[VAR_T] +=
					(is_neg ? -1 : +1) * power;

				tmp++;
				if (tmp1 != NULL)
					tmp = tmp1;
				break;

			default:
				printf("parse error: %u '%s'\n", i, tmp);
				exit(-1);
			}
		}
	}
	poly->num_terms = num_terms;
}

/*--------------------------------------------------------------------*/
static void
opt_data_alloc(multi_poly_t *objective, uint32 num_vars, uint32 deg,
		curr_poly_t *c, opt_data_t *s)
{
	uint32 i, j;

	memset(s, 0, sizeof(opt_data_t));

	s->num_vars = num_vars;
	s->objective = objective;
	s->c = c;

	for (i = 0; i <= 1; i++)
		s->params[i] = mpz_get_d(c->gmp_lina[i]);
	for (i = 0; i <= deg; i++)
		s->params[i+2] = mpz_get_d(c->gmp_a[i]);

	poly_find_max_powers(objective, s->max_powers);

	for (i = 0; i < num_vars; i++) {
		poly_alloc(&s->grad[i], objective->num_terms);
		poly_diff(objective, i, &s->grad[i]);
		poly_find_max_powers(&s->grad[i], s->max_powers);

		for (j = i; j < num_vars; j++) {
			poly_alloc(&s->hess[i][j], s->grad[i].num_terms);
			poly_diff(&s->grad[i], j, &s->hess[i][j]);
			poly_find_max_powers(&s->hess[i][j], s->max_powers);
		}
	}
}

/*--------------------------------------------------------------------*/
static void
opt_data_free(opt_data_t *s)
{
	uint32 i, j;

	for (i = 0; i < s->num_vars; i++)
		poly_free(s->grad + i);

	for (i = 0; i < s->num_vars; i++) {
		for (j = i; j < s->num_vars; j++) {
			poly_free(&s->hess[i][j]);
		}
	}
}

/*--------------------------------------------------------------------*/
static void
fill_powers(opt_data_t *s, double v[MAX_VARS], 
		double powers_in[MAX_VARS][MAX_VARDEGREE + 1])
{
	uint32 i, j;

	for (i = 0; i < MAX_VARS; i++) {
		double *powers = powers_in[i];
		int32 max_pow = s->max_powers[i];

		powers[1] = v[i];
		for (j = 2; j <= max_pow; j++)
			powers[j] = powers[j-1] * powers[1];
	}
}

/*--------------------------------------------------------------------*/
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
 
/*--------------------------------------------------------------------*/
static const double xlate_weights[6+1][6+1] = {
	{  1.,  1.,  1.,  1.,  1.,  1.,  1.},
	{  0.,  1.,  2.,  3.,  4.,  5.,  6.},
	{  0.,  0.,  1.,  3.,  6., 10., 15.},
	{  0.,  0.,  0.,  1.,  4., 10., 20.},
	{  0.,  0.,  0.,  0.,  1.,  5., 15.},
	{  0.,  0.,  0.,  0.,  0.,  1.,  6.},
	{  0.,  0.,  0.,  0.,  0.,  0.,  1.},
};

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

/*--------------------------------------------------------------------*/
static double 
callback(double v[MAX_VARS], void *extra)
{
	uint32 i;
	opt_data_t *data = (opt_data_t *)extra;
	double translated[MAX_POLY_DEGREE + 1];
	double apoly[MAX_POLY_DEGREE + 1];
	double s = v[0];
	double t = floor(v[1] + 0.5);
	double r0 = data->params[0];
	double r1 = data->params[1];

	if (s < 1.0)
		return 1e200;

	for (i = 0; i <= 6; i++)
		apoly[i] = data->params[i+2];

	for (i = 0; i <= 2; i++) {
		double c = floor(v[2 + i] + 0.5);
		apoly[i] += r0 * c;
		apoly[i+1] += r1 * c;
	}
	
	translate_d(translated, apoly, 6, t);
	return ifs_radial(translated, 6, s);
}

static double 
callback_hess(double v[MAX_VARS], 
		double grad[MAX_VARS], 
		double hess[MAX_VARS][MAX_VARS],
		void *extra)
{
	uint32 i, j;
	opt_data_t *s = (opt_data_t *)extra;
	double powers[MAX_VARS][MAX_VARDEGREE + 1];

	fill_powers(s, v, powers);

	for (i = 0; i < s->num_vars; i++)
		grad[i] = poly_eval(s->grad + i, s->params, powers);

	for (i = 0; i < s->num_vars; i++) {
		for (j = i; j < s->num_vars; j++) {
			hess[i][j] = hess[j][i] = 
				poly_eval(&s->hess[i][j], s->params, powers);
		}
	}

	return callback(v, extra);
}
				
/*--------------------------------------------------------------------*/
static void 
fixup(double v[MAX_VARS], void *extra)
{
	uint32 i;
	opt_data_t *s = (opt_data_t *)extra;
	curr_poly_t *c = s->c;

	for (i = 0; i <= 2; i++) {
		double cj = floor(v[i+2] + 0.5);
		mpz_set_d(c->gmp_help1, cj);
		mpz_addmul(c->gmp_a[i+1], c->gmp_help1, c->gmp_lina[1]);
		mpz_addmul(c->gmp_a[i], c->gmp_help1, c->gmp_lina[0]);
	}
	translate_gmp(c, c->gmp_a, 6, c->gmp_lina,
			(int64)(v[1] + 0.5));

	for (i = 0; i <= 1; i++)
		s->params[i] = mpz_get_d(c->gmp_lina[i]);
	for (i = 0; i <= 6; i++)
		s->params[i+2] = mpz_get_d(c->gmp_a[i]);

	v[1] = v[2] = v[3] = v[4] = 0;
}

/*--------------------------------------------------------------------*/
static void
find_start_point(opt_data_t *opt, double best[MAX_VARS])
{
	uint32 i, j, k;
	uint32 n;
	uint32 num_vars = opt->num_vars;
	double log_target = 0;
	multi_poly_t *objective = opt->objective;
	double mat[300][MAX_VARS];
	double rhs[300];
	double lstsqr[MAX_VARS][MAX_VARS];
	double lstsqr_rhs[MAX_VARS];

	for (i = 0; i < objective->num_terms; i++) {
		term_t *t = objective->terms + i;

		for (j = 0; j < MAX_VARS; j++) {
			if (t->powers[j])
				break;
		}
		if (j == MAX_VARS) {
			double curr_target = log(fabs(t->coeff));
			for (j = 0; j < MAX_PARAMS; j++) {
				if (t->params[j] == 0) 
					continue;
				if (opt->params[j] == 0) 
					break;
				curr_target += t->params[j] * log(
						fabs(opt->params[j]));
			}
			if (j == MAX_PARAMS)
				log_target = MAX(log_target, curr_target);
		}
	}

	memset(mat, 0, sizeof(mat));
	memset(rhs, 0, sizeof(rhs));
	memset(lstsqr, 0, sizeof(lstsqr));
	memset(lstsqr_rhs, 0, sizeof(lstsqr_rhs));

	for (i = n = 0; i < objective->num_terms; i++) {
		term_t *t = objective->terms + i;
		uint32 num_pow = 0;
		double curr_target = log(fabs(t->coeff));

		for (j = 0; j < MAX_PARAMS; j++) {
			if (t->params[j] == 0)
				continue;
			if (opt->params[j] == 0)
				break;
			curr_target += t->params[j] * 
					log(fabs(opt->params[j]));
		}
		if (j < MAX_PARAMS)
			continue;

		for (j = 0; j < MAX_VARS; j++) {
			if (t->powers[j] == 0)
				continue;
			mat[n][j] = t->powers[j];
			num_pow++;
		}
		if (num_pow > 0) {
			rhs[n++] = log_target - curr_target;
		}
	}

	for (i = 0; i < num_vars; i++) {
		double accum;
		for (j = 0; j < num_vars; j++) {
			for (k = 0, accum = 0; k < n; k++) {
				accum += mat[k][i] * mat[k][j];
			}
			lstsqr[i][j] = accum;
		}
		for (k = 0, accum = 0; k < n; k++)
			accum += mat[k][i] * rhs[k];
		lstsqr_rhs[i] = accum;
	}

	solve_dmatrix(lstsqr, best, lstsqr_rhs, num_vars);

	for (i = 0; i < num_vars; i++)
		best[i] = exp(best[i]);
}

/*--------------------------------------------------------------------*/
double
optimize_initial_deg6(double best[MAX_VARS], 
			curr_poly_t *c,
			uint32 degree)
{
	uint32 i, j, k;
	double best_score;
	multi_poly_t objective;
	opt_data_t opt_data;
	double lstsqr[MAX_VARS];
	double curr[MAX_VARS];
	uint32 num_vars = degree - 1;

	parse_objective(&objective, objective_deg6);

	opt_data_alloc(&objective, num_vars, 
			degree, c, &opt_data);

	find_start_point(&opt_data, lstsqr);

	for (i = 0; i <= degree; i++)
		mpz_set(c->gmp_c[i], c->gmp_a[i]);
	for (i = 0; i <= 1; i++)
		mpz_set(c->gmp_linc[i], c->gmp_lina[i]);

	best_score = 1e200;
	memset(best, 0, MAX_VARS * sizeof(double));

	for (i = 0; i < (1 << (num_vars - 1)); i++) {

		double curr_score = 1e200;
		double last_score;

		for (j = 0; j <= degree; j++)
			mpz_set(c->gmp_a[j], c->gmp_c[j]);
		for (j = 0; j <= 1; j++)
			mpz_set(c->gmp_lina[j], c->gmp_linc[j]);

		curr[0] = lstsqr[0];
		for (j = 1; j < num_vars; j++) {
			if (i & (1 << (j-1)))
				curr[j] = lstsqr[j];
			else
				curr[j] = -lstsqr[j];
		}

		for (j = 0; j < 20; j++) {
			last_score = curr_score;
			curr_score = minimize_hess(curr, num_vars, 1e-5, 20,
					   	callback, callback_hess,
						&opt_data);
			fixup(curr, &opt_data);
			curr_score = callback(curr, &opt_data);
			if (curr_score < best_score) {
				best_score = curr_score;
				best[0] = curr[0];

				for (k = 0; k <= degree; k++)
					mpz_set(c->gmp_b[k], c->gmp_a[k]);
				for (k = 0; k <= 1; k++)
					mpz_set(c->gmp_linb[k], c->gmp_lina[k]);
			}
			if (fabs(last_score - curr_score) < 
					1e-5 * fabs(last_score))
				break;
		}
#if 0
		printf("final: %le\n", curr_score);
#endif
	}

	for (i = 0; i <= degree; i++)
		mpz_set(c->gmp_a[i], c->gmp_b[i]);
	for (i = 0; i <= 1; i++)
		mpz_set(c->gmp_lina[i], c->gmp_linb[i]);
	mpz_neg(c->gmp_d, c->gmp_lina[0]);

	opt_data_free(&opt_data);
	poly_free(&objective);
	return best_score;
}
