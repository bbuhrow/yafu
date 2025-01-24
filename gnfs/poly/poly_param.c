/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: poly_param.c 1023 2018-08-19 00:30:42Z jasonp_sf $
--------------------------------------------------------------------*/

#include "poly_skew.h"

static const poly_param_t prebuilt_params_deg4[] = {

	{ 80, 3.00E+013, 2.00E+013, 1.00E-007,  1 * 60},
	{ 85, 3.00E+014, 4.00E+013, 6.50E-008,  1 * 60},
	{ 90, 2.50E+015, 5.00E+014, 3.80E-008,  3 * 60},
	{ 92, 9.96E+015, 4.87E+013, 2.80E-008,  4 * 60},
	{ 95, 2.92E+016, 1.75E+014, 1.85E-008,  6 * 60},
	{100, 1.52E+017, 1.13E+015, 8.75E-009, 12 * 60},
	{105, 8.97E+017, 7.44E+015, 4.40E-009, 23 * 60},
	{110, 6.60E+018, 3.00E+016, 1.40E-009, 45 * 60},
	{115, 1.00E+019, 1.00E+017, 4.09E-010, 90 * 60},
	{120, 3.00E+020, 5.00E+017, 2.14E-010, 180 * 60},
	{125, 1.00E+022, 1.00E+099, 1.12E-010, 8 * 3600},
};

static const poly_param_t prebuilt_params_deg5[] = {

	{100, 9.89E+015, 4.00E+012, 2.95E-009,   12 * 60},
	{105, 4.60E+016, 1.56E+013, 1.60E-009,   23 * 60},
	{110, 2.63E+017, 9.80E+013, 8.76E-010,   45 * 60},
	{115, 1.50E+018, 5.37E+014, 4.87E-010,   90 * 60},
	{120, 8.75E+018, 2.96E+015, 2.61E-010,  180 * 60},
	{125, 4.99E+019, 1.57E+016, 1.32E-010,  360 * 60},
	{130, 2.62E+020, 8.24E+016, 6.57E-011,  720 * 60},
	{135, 1.04E+021, 3.90E+017, 3.59E-011, 1320 * 60},
	{140, 4.43E+021, 2.02E+018, 1.75E-011, 2520 * 60},
	{145, 1.77E+022, 1.03E+019, 8.70E-012, 4900 * 60},
	{150, 7.09E+022, 5.25E+019, 4.35E-012, 9000 * 60},
	{159, 2.00E+024, 2.00E+022, 1.00E-012, 300 * 3600},
	{165, 8.00E+024, 2.00E+023, 5.00E-013, 300 * 3600},
	{170, 5.00E+025, 1.58E+024, 1.50E-013, 300 * 3600},
	{175, 3.00E+026, 1.00E+025, 1.00E-013, 300 * 3600},
	{180, 1.80E+027, 5.36E+025, 7.00E-014, 300 * 3600},
	{185, 1.00E+028, 3.12E+026, 2.00E-014, 300 * 3600},
	{190, 6.00E+028, 1.82E+027, 4.00E-015, 300 * 3600},
	{197, 1.00E+030, 1.00E+029, 2.00E-015, 300 * 3600},
	{200, 3.10E+030, 1.10E+029, 1.50E-015, 300 * 3600},
	{205, 2.00E+031, 5.70E+029, 5.50E-016, 300 * 3600},
	{210, 1.00E+032, 3.00E+030, 1.90E-016, 300 * 3600},
	{215, 6.00E+032, 1.50E+031, 7.00E-017, 300 * 3600},
	{220, 2.40E+033, 7.70E+031, 3.00E-017, 300 * 3600},
};

static const poly_param_t prebuilt_params_deg6[] = {

 	{200, 1.00E+026, 1.00E+025, 8.0e-018, 300 * 3600},
 	{205, 1.00E+027, 1.00E+026, 6.0e-019, 300 * 3600},
 	{230, 3.00E+029, 3.00E+029, 6.0e-019, 300 * 3600},
 	{235, 1.50E+030, 8.00E+029, 6.0e-019, 300 * 3600},

	/* better than nothing (marginally) */

 	{305, 1.00E+040, 1.00E+041, 0, 300 * 3600},
 	{310, 1.00E+042, 1.00E+043, 0, 300 * 3600},
};

/*--------------------------------------------------------------------*/
static void get_default_params(double digits, poly_param_t *params,
				const poly_param_t *defaults, 
				uint32 num_default_entries) {

	uint32 i;
	const poly_param_t *low, *high;
	double j, k, dist;
	double max_digits;

	/* if the input is too small (large), give up */

	if (digits < defaults[0].digits)
		return;

	max_digits = defaults[num_default_entries - 1].digits;
	if (digits >= max_digits)
		return;

	/* Otherwise the parameters to use are a weighted average 
	   of the two table entries the input falls between */

	for (i = 0; i < num_default_entries - 1; i++) {
		if (digits < defaults[i+1].digits)
			break;
	}

	low = &defaults[i];
	high = &defaults[i+1];
	dist = high->digits - low->digits;
	j = digits - low->digits;
	k = high->digits - digits;

	/* use exponential interpolation */

	params->digits = digits;
	params->stage1_norm = exp((log(low->stage1_norm) * k +
			           log(high->stage1_norm) * j) / dist);
	params->stage2_norm = exp((log(low->stage2_norm) * k +
			           log(high->stage2_norm) * j) / dist);
	params->final_norm = exp((log(low->final_norm) * k +
			           log(high->final_norm) * j) / dist);
	params->deadline = exp((log(low->deadline) * k +
			           log(high->deadline) * j) / dist);
}

/*------------------------------------------------------------------*/
void get_poly_params(msieve_obj *obj, mpz_t n,
			uint32 *degree_out, 
			poly_param_t *params_out) {

	poly_param_t params;
	double digits = log(mpz_get_d(n)) / log(10.0);

	uint32 degree = 0;

	memset(&params, 0, sizeof(params));

	/* see if the degree is specified */

	if (obj->nfs_args != NULL) {
		const char *tmp = strstr(obj->nfs_args, "polydegree=");

		if (tmp != NULL) {
			degree = strtoul(tmp + 11, NULL, 10);
			logprintf(obj, "setting degree to %u\n", degree);
		}
	}

	/* if not, choose the degree automatically */

	if (degree == 0) {
		uint32 bits = mpz_sizeinbase(n, 2);

		if (digits < 108.0)
			degree = 4;
		else if (digits < 220.0)
			degree = 5;
		else
			degree = 6;
	}

	/* get default parameters, if any */

	switch (degree) {
	case 4:
		get_default_params(digits, &params, prebuilt_params_deg4,
				sizeof(prebuilt_params_deg4) / 
					sizeof(poly_param_t));
		break;

	case 5:
		get_default_params(digits, &params, prebuilt_params_deg5,
				sizeof(prebuilt_params_deg5) / 
					sizeof(poly_param_t));
		break;

	case 6:
		get_default_params(digits, &params, prebuilt_params_deg6,
				sizeof(prebuilt_params_deg6) / 
					sizeof(poly_param_t));
		break;

	default:
		printf("error: invalid degree %u\n", degree);
		exit(-1);
	}

	/* override with user-suplied params */

	if (obj->nfs_args != NULL) {
		const char *tmp;

		tmp = strstr(obj->nfs_args, "stage1_norm=");
		if (tmp != NULL)
			params.stage1_norm = strtod(tmp + 12, NULL);

		tmp = strstr(obj->nfs_args, "stage2_norm=");
		if (tmp != NULL)
			params.stage2_norm = strtod(tmp + 12, NULL);

		tmp = strstr(obj->nfs_args, "min_evalue=");
		if (tmp != NULL)
			params.final_norm = strtod(tmp + 11, NULL);

		tmp = strstr(obj->nfs_args, "poly_deadline=");
		if (tmp != NULL)
			params.deadline = strtoul(tmp + 14, NULL, 10);
	}

	logprintf(obj, "polynomial degree: %u\n", degree);
	logprintf(obj, "max stage 1 norm: %.2e\n", params.stage1_norm);
	logprintf(obj, "max stage 2 norm: %.2e\n", params.stage2_norm);
	logprintf(obj, "min E-value: %.2e\n", params.final_norm);
	logprintf(obj, "poly select deadline: %u\n", params.deadline);
	*degree_out = degree;
	*params_out = params;
}
