/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: dickman.c 23 2009-07-20 02:59:07Z jasonp_sf $
--------------------------------------------------------------------*/

#include "ms_common.h"

#define MAX_LINE 30
#define MAX_COEFFS 55

/*------------------------------------------------------------------*/
void dickman_free(dickman_t *aux) {

	free(aux->lines);
	free(aux->coeffs);
}

/*------------------------------------------------------------------*/
void dickman_init(dickman_t *aux) {

	uint32 i, j, k;
	uint32 num_coeffs;
	uint32 num_coeffs_alloc;
	double last_coeffs[MAX_COEFFS];
	double *coeffs;
	double sum, s, t;

	/* initialize */

	num_coeffs = 0;
	num_coeffs_alloc = 1000;
	aux->coeffs = (double *)xmalloc(num_coeffs_alloc *
					sizeof(double));
	aux->lines = (dickman_line_t *)xmalloc((MAX_LINE + 1) *
					sizeof(dickman_line_t));

	/* dickman(x) is 1.0 for x <= 1; for x <= 2 the
	   value of dickman(x) is known explicitly */

	last_coeffs[0] = 1.0 - M_LN2;
	t = 0.5;
	for (i = 1; i < MAX_COEFFS; i++, t *= 0.5) {
		last_coeffs[i] = t / i;
	}

	/* save only enough coefficients for x=2 to 
	   achieve the desired accuracy, but use all the
	   rest for below to avoid roundoff error for
	   larger x 

	   We know the largest argument to the series
	   version of dickman(x) will be 1.0 and all terms
	   are positive, so this bounds the work needed */

	coeffs = aux->coeffs;
	sum = coeffs[0] = last_coeffs[0];
	num_coeffs = MAX_COEFFS;

	for (i = 1; i < MAX_COEFFS; i++) {
		coeffs[i] = last_coeffs[i];
		sum += coeffs[i];

		if (coeffs[i] < sum * 0.5 * DICKMAN_ACCURACY) {
			num_coeffs = i + 1;
			break;
		}
	}
	aux->lines[2].num_coeffs = num_coeffs;
	aux->lines[2].coeff_offset = 0;

	/* proceed with the rest of the integer aguments
	   to dickman(x) */

	for (i = 3; i <= MAX_LINE; i++) {
		dickman_line_t *line = aux->lines + i;
		double recip_di = 1.0 / i;

		if (num_coeffs + MAX_COEFFS >= num_coeffs_alloc) {
			num_coeffs_alloc *= 2;
			aux->coeffs = (double *)xrealloc(aux->coeffs,
							num_coeffs_alloc *
							sizeof(double));
		}

		/* derive the coefficients for dickman(x) from
		   those used in dickman(x-1) */

		sum = 0;
		for (j = MAX_COEFFS - 1; j; j--) {

			s = 0;
			t = recip_di / j;
			for (k = j - 1; (int32)k >= 0; k--) {
				s += t * last_coeffs[k];
				t *= recip_di;
			}

			last_coeffs[j] = s;
			sum += s / (j + 1);
		}
		last_coeffs[j] = sum / (i - 1);

		/* save enough coefficients to achieve the
		   desired accuracy for fractional arguments */

		coeffs = aux->coeffs + num_coeffs;
		sum = coeffs[0] = last_coeffs[0];
		line->coeff_offset = num_coeffs;
		line->num_coeffs = MAX_COEFFS;

		for (j = 1; j < MAX_COEFFS; j++) {
			coeffs[j] = last_coeffs[j];
			sum += coeffs[j];

			if (coeffs[j] < sum * 0.5 * DICKMAN_ACCURACY) {
				line->num_coeffs = j + 1;
				break;
			}
		}
		num_coeffs += line->num_coeffs;
	}
}

/*------------------------------------------------------------------*/
double dickman(dickman_t *aux, double arg) {

	uint32 i;
	double int_arg, frac_arg;
	double sum, term;
	dickman_line_t *line;
	uint32 num_coeffs;
	double *coeffs;

	if (arg <= 1.0)
		return 1.0;

	if (arg > MAX_LINE)
		return 0.0;

	int_arg = ceil(arg);
	frac_arg = int_arg - arg;
	line = aux->lines + (int32)int_arg;
	num_coeffs = line->num_coeffs;
	coeffs = aux->coeffs + line->coeff_offset;

	sum = coeffs[0];
	term = frac_arg;

	for (i = 1; i < num_coeffs; i++) {

		double iterm = term * coeffs[i];

		sum += iterm;

		if (iterm < sum * 0.5 * DICKMAN_ACCURACY)
			break;

		term *= frac_arg;
	}

	return sum;
}
