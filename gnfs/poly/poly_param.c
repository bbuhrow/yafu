/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: poly_param.c 817 2012-11-11 14:58:29Z jasonp_sf $
--------------------------------------------------------------------*/

#include "poly_skew.h"

static const poly_param_t prebuilt_params_deg4[] = {

	/* determined by experiment */

	{ 80, 3.00E+013, 2.00E+013, 1.00E-007, 4 * 60},
	{ 85, 3.00E+014, 4.00E+013, 6.50E-008, 8 * 60},
	{ 90, 2.50E+015, 5.00E+014, 3.80E-008, 12 * 60},
	{ 95, 1.50E+016, 1.00E+015, 1.50E-008, 15 * 60},
	{100, 2.10E+017, 4.00E+015, 8.30E-009, 23 * 60},
	{105, 1.00E+018, 1.50E+016, 3.00E-009, 30 * 60},
	{110, 6.00E+018, 5.00E+016, 1.00E-009, 60 * 60},
};

static const poly_param_t prebuilt_params_deg5[] = {

	/* shamelessly stolen from GGNFS */

	{100, 5.56E+014, 2.86E+013, 2.80E-009, 23 * 60},
	{101, 8.13E+014, 4.08E+013, 2.50E-009, 24 * 60},
	{102, 1.19E+015, 5.82E+013, 2.21E-009, 25 * 60},
	{103, 1.74E+015, 8.31E+013, 1.94E-009, 26 * 60},
	{104, 2.55E+015, 1.18E+014, 1.71E-009, 28 * 60},
	{105, 3.74E+015, 1.69E+014, 1.50E-009, 30 * 60},
	{106, 5.47E+015, 2.41E+014, 1.32E-009, 35 * 60},
	{107, 8.01E+015, 3.44E+014, 1.16E-009, 40 * 60},
	{108, 1.17E+016, 4.91E+014, 1.02E-009, 48 * 60},
	{109, 1.72E+016, 7.00E+014, 8.92E-010, 54 * 60},
	{110, 2.52E+016, 9.98E+014, 7.83E-010, 60 * 60},
	{111, 3.68E+016, 1.42E+015, 6.88E-010, 69 * 60},
	{112, 5.39E+016, 2.03E+015, 6.04E-010, 79 * 60},
	{113, 7.90E+016, 2.90E+015, 5.31E-010, 93 * 60},
	{114, 1.16E+017, 4.13E+015, 4.66E-010, 106 * 60},
	{115, 1.69E+017, 5.90E+015, 4.09E-010, 120 * 60},
	{116, 2.48E+017, 8.41E+015, 3.60E-010, 135 * 60},
	{117, 3.63E+017, 1.20E+016, 3.16E-010, 165 * 60},
	{118, 5.31E+017, 1.71E+016, 2.77E-010, 190 * 60},
	{119, 7.78E+017, 2.44E+016, 2.44E-010, 215 * 60},
	{120, 1.14E+018, 3.48E+016, 2.14E-010, 240 * 60},
	{121, 1.67E+018, 4.97E+016, 1.88E-010, 275 * 60},
	{122, 2.44E+018, 7.09E+016, 1.65E-010, 310 * 60},
	{123, 3.58E+018, 1.01E+017, 1.45E-010, 370 * 60},
	{124, 5.24E+018, 1.44E+017, 1.27E-010, 425 * 60},
	{125, 7.67E+018, 2.06E+017, 1.12E-010, 8 * 3600},
	{126, 1.12E+019, 2.94E+017, 9.82E-011, 9 * 3600},
	{127, 1.64E+019, 4.19E+017, 8.62E-011, 10 * 3600},
	{128, 2.41E+019, 5.98E+017, 7.57E-011, 12 * 3600},
	{129, 3.52E+019, 8.52E+017, 6.65E-011, 14 * 3600},
	{130, 5.16E+019, 1.22E+018, 5.84E-011, 16 * 3600},
	{131, 7.55E+019, 1.73E+018, 5.13E-011, 18 * 3600},
	{132, 1.11E+020, 2.47E+018, 4.51E-011, 20 * 3600},
	{133, 1.62E+020, 3.53E+018, 3.96E-011, 24 * 3600},
	{134, 2.37E+020, 5.04E+018, 3.48E-011, 28 * 3600},
	{135, 3.47E+020, 7.18E+018, 3.05E-011, 32 * 3600},
	{136, 5.08E+020, 1.02E+019, 2.68E-011, 36 * 3600},
	{137, 7.44E+020, 1.46E+019, 2.36E-011, 40 * 3600},
	{138, 1.09E+021, 2.09E+019, 2.07E-011, 48 * 3600},
	{139, 1.60E+021, 2.98E+019, 1.82E-011, 56 * 3600},
	{140, 2.34E+021, 4.24E+019, 1.60E-011, 64 * 3600},
	{141, 3.42E+021, 6.05E+019, 1.40E-011, 70 * 3600},
	{142, 5.01E+021, 8.64E+019, 1.23E-011, 76 * 3600},
	{143, 7.33E+021, 1.23E+020, 1.08E-011, 84 * 3600},
	{144, 1.07E+022, 1.76E+020, 9.50E-012, 92 * 3600},
	{145, 1.57E+022, 2.51E+020, 8.34E-012, 100 * 3600},
	{146, 2.30E+022, 3.58E+020, 7.33E-012, 120 * 3600},
	{147, 3.37E+022, 5.10E+020, 6.43E-012, 140 * 3600},
	{148, 4.94E+022, 7.28E+020, 5.65E-012, 160 * 3600},
	{149, 7.23E+022, 1.04E+021, 4.96E-012, 180 * 3600},
	{150, 1.06E+023, 1.48E+021, 4.36E-012, 200 * 3600},
	{151, 1.55E+023, 2.11E+021, 3.83E-012, 220 * 3600},
	{152, 2.27E+023, 3.01E+021, 3.36E-012, 240 * 3600},
	{153, 3.32E+023, 4.30E+021, 2.95E-012, 260 * 3600},
	{154, 8.00E+023, 4.90E+021, 2.40E-012, 280 * 3600},
	{155, 1.00E+024, 6.80E+021, 2.00E-012, 300 * 3600},

	/* contributed by Tom Womack */

	{159, 2.00E+024, 2.00E+022, 1.00E-012, 300 * 3600},
	{165, 8.00E+024, 2.00E+023, 5.00E-013, 300 * 3600},
	{170, 5.00E+025, 1.58E+024, 1.50E-013, 300 * 3600},

	/* contributed by Serge Batalov */

	{175, 3.00E+026, 1.00E+025, 1.00E-013, 300 * 3600},
	{180, 1.80E+027, 5.36E+025, 7.00E-014, 300 * 3600},
	{185, 1.00E+028, 3.12E+026, 2.00E-014, 300 * 3600},
	{190, 6.00E+028, 1.82E+027, 4.00E-015, 300 * 3600},

	/* from Tom Womack */

	{197, 1.00E+030, 1.00E+029, 2.00E-015, 300 * 3600},

	/* only a guess */

	{200, 3.10E+030, 1.10E+029, 1.50E-015, 300 * 3600},
	{205, 2.00E+031, 5.70E+029, 5.50E-016, 300 * 3600},
	{210, 1.00E+032, 3.00E+030, 1.90E-016, 300 * 3600},
	{215, 6.00E+032, 1.50E+031, 7.00E-017, 300 * 3600},
	{220, 2.40E+033, 7.70E+031, 3.00E-017, 300 * 3600},
};

static const poly_param_t prebuilt_params_deg6[] = {

	{140, 2.34E+017, 5.00E+018, 1.7e-012, 300 * 3600},
	{141, 2.34E+017, 5.00E+018, 1.7e-012, 300 * 3600},

	/* contributed by Paul Leyland */

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
	double digits;
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

		if (bits < 363) 		/* <= 110 digits */
			degree = 4;
		else if (bits < 726) 		/* 110-220 digits */
			degree = 5;
		else				/* 220+ digits */
			degree = 6;
	}

	/* get default parameters, if any */

	digits = log(mpz_get_d(n)) / log(10.0);

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
