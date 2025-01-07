/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: poly_skew.c 925 2013-07-20 17:01:02Z brgladman $
--------------------------------------------------------------------*/

#include "poly_skew.h"

/*------------------------------------------------------------------*/
/* callback for root optimization */

typedef struct {
	FILE *all_poly_file;
	poly_config_t *config;
} rootopt_callback_data_t;

static void
rootopt_callback(void *extra, uint32 degree, 
		mpz_t * coeff1, mpz_t * coeff2,
		double skewness, double size_score,
		double root_score, double combined_score,
		uint32 num_real_roots)
{
	uint32 i;
	poly_select_t poly;
	mpz_poly_t *rpoly;
	mpz_poly_t *apoly;
	rootopt_callback_data_t *data = (rootopt_callback_data_t *)extra;
	poly_config_t *config = data->config;

	memset(&poly, 0, sizeof(poly_select_t));
	rpoly = &poly.rpoly;
	apoly = &poly.apoly;
	mpz_poly_init(rpoly);
	mpz_poly_init(apoly);

	rpoly->degree = 1;
	for (i = 0; i <= 1; i++)
		mpz_set(rpoly->coeff[i], coeff2[i]);

	apoly->degree = degree;
	for (i = 0; i <= degree; i++)
		mpz_set(apoly->coeff[i], coeff1[i]);

	poly.root_score = root_score;
	poly.size_score = size_score;
	poly.combined_score = combined_score;
	poly.skewness = skewness;

	printf("save %le %.4lf %.2lf %le rroots %u\n", size_score,
			root_score, skewness, combined_score,
			num_real_roots);

	fprintf(data->all_poly_file, 
		"# norm %le alpha %lf e %.3le rroots %u\nskew: %.2lf\n", 
		size_score, root_score, combined_score, 
		num_real_roots, skewness);
	for (i = 0; i <= degree; i++)
		gmp_fprintf(data->all_poly_file, "c%u: %Zd\n", i, coeff1[i]);
	for (i = 0; i <= 1; i++)
		gmp_fprintf(data->all_poly_file, "Y%u: %Zd\n", i, coeff2[i]);
	fflush(data->all_poly_file);

	save_poly(config, &poly);
	mpz_poly_free(rpoly);
	mpz_poly_free(apoly);
}

/*------------------------------------------------------------------*/
/* callbacks for size optimization */

typedef struct {
	poly_rootopt_t *rootopt;
	rootopt_callback_data_t *rootopt_callback;
} sizeopt_callback_data_t;

static void sizeopt_callback(uint32 deg, mpz_t *alg_coeffs, mpz_t *rat_coeffs, 
				double sizeopt_norm, double projective_alpha, 
				void *extra)
{
	sizeopt_callback_data_t *callback = (sizeopt_callback_data_t *)extra;

	poly_rootopt_run(callback->rootopt, alg_coeffs, 
			rat_coeffs, sizeopt_norm, projective_alpha);
}

static void sizeopt_callback_log(uint32 deg, mpz_t *alg_coeffs, mpz_t *rat_coeffs, 
				double sizeopt_norm, double projective_alpha, 
				void *extra)
{
	uint32 i;
	FILE *mfile = (FILE *)extra;

	for (i = deg; (int32)i >= 0; i--)
		gmp_fprintf(mfile, "%Zd ", alg_coeffs[i]);

	gmp_fprintf(mfile, "%Zd %Zd %.2lf %le\n", rat_coeffs[1], 
				rat_coeffs[0], projective_alpha, 
				exp(projective_alpha) * sizeopt_norm);
	fflush(mfile);
}

/*------------------------------------------------------------------*/
/* callbacks for stage 1 */

static void stage1_callback(mpz_t ad, mpz_t p, mpz_t m, void *extra) {
	
	poly_sizeopt_t *data = (poly_sizeopt_t *)extra;

	//gmp_printf("%Zd %Zd %Zd\n", ad, p, m);
	poly_sizeopt_run(data, ad, p, m);
}

static void stage1_callback_log(mpz_t ad, mpz_t p, mpz_t m, void *extra) {
	
	FILE *mfile = (FILE *)extra;
	gmp_printf("%Zd %Zd %Zd\n", ad, p, m);
	gmp_fprintf(mfile, "%Zd %Zd %Zd\n", ad, p, m);
	fflush(mfile);
}

/*------------------------------------------------------------------*/
void find_poly_core(msieve_obj *obj, mpz_t n,
			poly_param_t *params,
			poly_config_t *config,
			uint32 degree) {

	poly_stage1_t stage1_data;
	poly_sizeopt_t sizeopt_data;
	poly_rootopt_t rootopt_data;
	sizeopt_callback_data_t sizeopt_callback_data;
	rootopt_callback_data_t rootopt_callback_data;
	char buf[2048];
	FILE *stage1_outfile = NULL;
	FILE *sizeopt_outfile = NULL;
	const char *lower_limit = NULL;
	const char *upper_limit = NULL;

	/* make sure the configured stages have the bounds that
	   they need. We only have to check maximum bounds, since
	   if not provided then no polynomials will be found */

	if ((obj->flags & MSIEVE_FLAG_NFS_POLY1) &&
	    params->stage1_norm == 0) {
		printf("error: stage 1 bound not provided\n");
		exit(-1);
	}
	if (((obj->flags & MSIEVE_FLAG_NFS_POLYSIZE) ||
	     (obj->flags & MSIEVE_FLAG_NFS_POLYROOT)) &&
	    params->stage2_norm == 0) {
		printf("error: stage 2 size bound not provided\n");
		exit(-1);
	}
	if ((obj->flags & MSIEVE_FLAG_NFS_POLY1) &&
	    (obj->flags & MSIEVE_FLAG_NFS_POLYROOT) &&
	    !(obj->flags & MSIEVE_FLAG_NFS_POLYSIZE)) {
		printf("error: middle polyselect stage missing\n");
		exit(-1);
	}

	/* parse arguments */

	if (obj->nfs_args != NULL) {
		const char *tmp;

		tmp = strstr(obj->nfs_args, "min_coeff=");
		if (tmp != NULL)
			lower_limit = tmp + 10;

		tmp = strstr(obj->nfs_args, "max_coeff=");
		if (tmp != NULL)
			upper_limit = tmp + 10;

		/* old-style 'X,Y' format */
		tmp = strchr(obj->nfs_args, ',');
		if (tmp != NULL) {
			lower_limit = tmp - 1;
			while (lower_limit > obj->nfs_args &&
				isdigit(lower_limit[-1])) {
				lower_limit--;
			}
			upper_limit = tmp + 1;
		}

	}

	/* set up stage 1 */

	if (obj->flags & MSIEVE_FLAG_NFS_POLY1) {

		double coeff_scale = 3.0;

		/* if we're doing both stage 1 and size optimization, 
		   every hit in stage 1 is immediately forwarded. For 
		   stage 1 alone, all the stage 1 hits are buffered
		   to file first */

		if (obj->flags & MSIEVE_FLAG_NFS_POLYSIZE) {
			poly_stage1_init(&stage1_data, 
					stage1_callback, &sizeopt_data);
		}
		else {
			sprintf(buf, "%s.m", obj->savefile.name);
			stage1_outfile = fopen(buf, "a");
			if (stage1_outfile == NULL) {
				printf("error: cannot open poly1 file\n");
				exit(-1);
			}
			poly_stage1_init(&stage1_data, stage1_callback_log, 
					stage1_outfile);
		}

		/* fill stage 1 data */

		mpz_set(stage1_data.gmp_N, n);
		stage1_data.degree = degree;
		stage1_data.norm_max = params->stage1_norm;
		stage1_data.deadline = params->deadline;

		if (lower_limit != NULL)
			gmp_sscanf(lower_limit, "%Zd",
					stage1_data.gmp_high_coeff_begin);
		else
			mpz_set_ui(stage1_data.gmp_high_coeff_begin, 
					(unsigned long)1);

		if (upper_limit != NULL)
			gmp_sscanf(upper_limit, "%Zd",
					stage1_data.gmp_high_coeff_end);
		else
			mpz_set_d(stage1_data.gmp_high_coeff_end,
				pow(mpz_get_d(stage1_data.gmp_N), 1.0 / 
					(double)(degree * (degree - 1))) / 
				coeff_scale );

		if (stage1_data.deadline != 0)
			logprintf(obj, "time limit set to %.2f CPU-hours\n",
				stage1_data.deadline / 3600.0);

		{ /* SB: tried L[1/3,c] fit; it is no better than this */
			double e0 = (params->digits >= 121) ? 
						(0.0607 * params->digits + 2.25):
				                (0.0526 * params->digits + 3.23);
			if (degree == 4)
				e0 = 0.0625 * params->digits + 1.69;
			e0 = exp(-log(10) * e0); 
#ifdef HAVE_CUDA
			e0 *= 1.15;
#endif
			logprintf(obj, "expecting poly E from %.2le to > %.2le\n",
				e0, 1.15 * e0);
			/* seen exceptional polys with +40% but that's */
			/* very rare. The fit is good for 88..232 digits */
		}
 
		logprintf(obj, "searching leading coefficients from "
				"%.0lf to %.0lf\n",
				mpz_get_d(stage1_data.gmp_high_coeff_begin),
				mpz_get_d(stage1_data.gmp_high_coeff_end));
	}

	/* set up size optimization */

	if (obj->flags & MSIEVE_FLAG_NFS_POLYSIZE) {

		/* if we're doing both size and root optimization, 
		   every size opt hit in is immediately forwarded. For 
		   size optimization alone, all the hits are buffered
		   to file first */

		if (obj->flags & MSIEVE_FLAG_NFS_POLYROOT) {
			sizeopt_callback_data.rootopt = &rootopt_data;
			sizeopt_callback_data.rootopt_callback = &rootopt_callback_data;

			poly_sizeopt_init(&sizeopt_data, 
					sizeopt_callback, &sizeopt_callback_data);
		}
		else {
			sprintf(buf, "%s.ms", obj->savefile.name);
			sizeopt_outfile = fopen(buf, "a");
			if (sizeopt_outfile == NULL) {
				printf("error: cannot open sizeopt file\n");
				exit(-1);
			}
			poly_sizeopt_init(&sizeopt_data, sizeopt_callback_log, 
					sizeopt_outfile);
		}

		/* fill size optimization data */

		mpz_set(sizeopt_data.gmp_N, n);
		sizeopt_data.degree = degree;
		sizeopt_data.max_stage1_norm = params->stage1_norm;
		sizeopt_data.max_sizeopt_norm = params->stage2_norm;
	}

	/* set up root optimization */

	if (obj->flags & MSIEVE_FLAG_NFS_POLYROOT) {

		poly_rootopt_init(&rootopt_data, obj, rootopt_callback, 
				&rootopt_callback_data);

		/* fill root optimization data */

		mpz_set(rootopt_data.gmp_N, n);
		rootopt_data.degree = degree;
		rootopt_data.max_sizeopt_norm = params->stage2_norm;
		rootopt_data.min_e = params->final_norm;
		rootopt_data.min_e_bernstein = 0;

		/* smaller problems (especially degree 5) run much faster 
		   when Bernstein's scoring function is used to weed out 
		   polynomials that probably cannot get their Murphy score 
		   optimized enough to exceed the E-value bound */

		if (degree == 4) {
			rootopt_data.min_e_bernstein = params->final_norm / 
					(30000 * pow(2.5, (params->digits - 90) / 5));
		}
		else if (degree == 5) {
			rootopt_data.min_e_bernstein = params->final_norm / 
					(2.2 * pow(3.0, (params->digits - 100) / 10));
		}

		sprintf(buf, "%s.p", obj->savefile.name);
		rootopt_callback_data.config = config;
		rootopt_callback_data.all_poly_file = fopen(buf, "a");
		if (rootopt_callback_data.all_poly_file == NULL) {
			printf("error: cannot open root opt file\n");
			exit(-1);
		}
	}


	if (obj->flags & MSIEVE_FLAG_NFS_POLY1) {

		poly_stage1_run(obj, &stage1_data);

		poly_stage1_free(&stage1_data);
		if (obj->flags & MSIEVE_FLAG_NFS_POLYROOT) {
			fclose(rootopt_callback_data.all_poly_file);
			poly_sizeopt_free(&sizeopt_data);
			poly_rootopt_free(&rootopt_data);
		}
		else if (obj->flags & MSIEVE_FLAG_NFS_POLYSIZE) {
			fclose(sizeopt_outfile);
			poly_sizeopt_free(&sizeopt_data);
		}
	}
	else if (obj->flags & MSIEVE_FLAG_NFS_POLYSIZE) {
		mpz_t ad, p, m;

		mpz_init(ad);
		mpz_init(p);
		mpz_init(m);

		sprintf(buf, "%s.m", obj->savefile.name);
		stage1_outfile = fopen(buf, "r");
		if (stage1_outfile == NULL) {
			printf("error: cannot open sizeopt input file\n");
			exit(-1);
		}

		while (1) {
			if (fgets(buf, sizeof(buf), stage1_outfile) == NULL)
				break;

			if (gmp_sscanf(buf, "%Zd %Zd %Zd", ad, p, m) != 3)
				continue;

			//gmp_printf("poly %Zd %Zd %Zd\n", ad, p, m);
			poly_sizeopt_run(&sizeopt_data, ad, p, m);
			if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)
				break;
		}

		mpz_clear(ad);
		mpz_clear(p);
		mpz_clear(m);
		fclose(stage1_outfile);
		poly_sizeopt_free(&sizeopt_data);
		if (obj->flags & MSIEVE_FLAG_NFS_POLYROOT) {
			fclose(rootopt_callback_data.all_poly_file);
			poly_rootopt_free(&rootopt_data);
		}
	}
	else if (obj->flags & MSIEVE_FLAG_NFS_POLYROOT) {
		uint32 i;
		mpz_t alg_coeffs[MAX_POLY_DEGREE + 1];
		mpz_t rat_coeffs[2];

		for (i = 0; i <= MAX_POLY_DEGREE; i++)
			mpz_init(alg_coeffs[i]);
		for (i = 0; i <= 1; i++)
			mpz_init(rat_coeffs[i]);

		sprintf(buf, "%s.ms", obj->savefile.name);
		sizeopt_outfile = fopen(buf, "r");
		if (sizeopt_outfile == NULL) {
			printf("error: cannot open rootopt input file\n");
			exit(-1);
		}

		while (1) {
			int c;
			char *tmp = buf;

			/* this reading scheme is unintuitive but more
			   resistant to corrupted input lines */

			if (fgets(tmp, sizeof(buf), sizeopt_outfile) == NULL)
				break;

			for (i = 0; i <= degree; i++) {
				if (gmp_sscanf(tmp, "%Zd%n", 
					    alg_coeffs[degree - i], &c) != 1)
						break;
				tmp += c;
			}

			if (gmp_sscanf(tmp, "%Zd %Zd", 
					rat_coeffs[1], rat_coeffs[0]) != 2)
				continue;

			/* make later code figure out the initial poly score */

			poly_rootopt_run(&rootopt_data, alg_coeffs, rat_coeffs, 0, 0);
			if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)
				break;
		}

		for (i = 0; i <= MAX_POLY_DEGREE; i++)
			mpz_clear(alg_coeffs[i]);
		for (i = 0; i <= 1; i++)
			mpz_clear(rat_coeffs[i]);
		fclose(sizeopt_outfile);
		fclose(rootopt_callback_data.all_poly_file);
		poly_rootopt_free(&rootopt_data);
	}
}
