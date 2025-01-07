/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: poly.c 937 2013-08-08 00:19:28Z jasonp_sf $
--------------------------------------------------------------------*/

#include "poly.h"

/*------------------------------------------------------------------*/
int32 read_poly(msieve_obj *obj, mpz_t n,
	       mpz_poly_t *rat_poly,
	       mpz_poly_t *alg_poly,
	       double *skewness) {
	
	uint32 i;
	FILE *fp;
	char buf[BIGNUM_BUF_SIZE];
	mpz_t read_n;
	mpz_t val, rpow;
	int32 status = 0;

	fp = fopen(obj->nfs_fbfile_name, "r");
	if (fp == NULL)
		return -1;
	
	buf[0] = 0;
	fgets(buf, (int)sizeof(buf), fp);
	if (buf[0] != 'N') {
		fclose(fp);
		logprintf(obj, "warning: factor base file uninitialized\n");
		return -1;
	}

	/* check that the polynomial is for the 
	   right number */

	mpz_init(read_n);
	mpz_init(val);
	mpz_init(rpow);

	gmp_sscanf(buf + 2, "%Zd", read_n);
	if (mpz_cmp(read_n, n) != 0) {
		fclose(fp);
		logprintf(obj, "warning: NFS input not found in "
				"factor base file\n");
		status = -1;
		goto finished;
	}

	/* read in skewness if present */

	fgets(buf, (int)sizeof(buf), fp);
	if (buf[0] == 'S') {
		if (skewness != NULL)
			*skewness = atof(buf + 5);
		fgets(buf, (int)sizeof(buf), fp);
	}
	else if (skewness != NULL) {
		*skewness = 1.0;
	}

	/* read one coefficient per line; 'R<number>' is
	   for rational coefficients, 'A<number>' for algebraic */

	while ((buf[0] == 'R' || buf[0] == 'A') && isdigit(buf[1])) {
		mpz_t *read_coeff;
		char *tmp;

		i = buf[1] - '0';
		if (i > MAX_POLY_DEGREE) {
			fclose(fp);
			logprintf(obj, "warning: polynomial degree exceeds "
					"%d\n", MAX_POLY_DEGREE);
			exit(-1);
		}

		if (buf[0] == 'R')
			read_coeff = rat_poly->coeff + i;
		else
			read_coeff = alg_poly->coeff + i;

		tmp = buf + 2;
		while (isspace(*tmp))
			tmp++;

		gmp_sscanf(tmp, "%Zd", *read_coeff);
		if (fgets(buf, (int)sizeof(buf), fp) == NULL)
			break;
	}

	for (i = MAX_POLY_DEGREE + 1; i; i--) {
		if (mpz_cmp_ui(rat_poly->coeff[i-1], 0) != 0)
			break;
	}
	if (i)
		rat_poly->degree = i - 1;

	for (i = MAX_POLY_DEGREE + 1; i; i--) {
		if (mpz_cmp_ui(alg_poly->coeff[i-1], 0) != 0)
			break;
	}
	if (i)
		alg_poly->degree = i - 1;

	fclose(fp);

	if (rat_poly->degree == 0 || alg_poly->degree == 0) {
		logprintf(obj, "error: polynomial is missing or corrupt\n");
		exit(-1);
	}
	if (rat_poly->degree != 1) {
		logprintf(obj, "error: no support for nonlinear "
				"rational polynomials\n");
		exit(-1);
	}
	
	/* plug the rational polynomial coefficients into the 
	   algebraic polynomial */

	i = alg_poly->degree;
	mpz_set(val, alg_poly->coeff[i]);
	mpz_set(rpow, rat_poly->coeff[1]);
	mpz_neg(rat_poly->coeff[0], rat_poly->coeff[0]);

	while (--i) {
		mpz_mul(val, val, rat_poly->coeff[0]);
		mpz_addmul(val, alg_poly->coeff[i], rpow);
		mpz_mul(rpow, rpow, rat_poly->coeff[1]);
	}
	mpz_mul(val, val, rat_poly->coeff[0]);
	mpz_addmul(val, alg_poly->coeff[i], rpow);
	mpz_neg(rat_poly->coeff[0], rat_poly->coeff[0]);

	/* verify that result % N == 0. 

	   The only place where we do any mod-N arithmetic is the 
	   NFS square root, which will not work if N has additional 
	   factors that are not reflected in the polynomials */

	mpz_mod(val, val, n);
	if (mpz_cmp_ui(val, 0) != 0) {
		logprintf(obj, "error: NFS input does not match polynomials\n");
		logprintf(obj, "check that input doesn't have small factors\n");
		exit(-1);
	}

finished:
	mpz_clear(read_n);
	mpz_clear(val);
	mpz_clear(rpow);
	return 0;
}

/*------------------------------------------------------------------*/
void write_poly(msieve_obj *obj, mpz_t n,
	       mpz_poly_t *rat_poly,
	       mpz_poly_t *alg_poly,
	       double skewness) {
	
	/* log a generated polynomial to the factor base file */

	uint32 i;
	FILE *fp;

	fp = fopen(obj->nfs_fbfile_name, "w");
	if (fp == NULL) {
		printf("error; cannot open factor base file '%s'\n",
					obj->nfs_fbfile_name);
		exit(-1);
	}

	gmp_fprintf(fp, "N %Zd\n", n);
	if (skewness > 0)
		fprintf(fp, "SKEW %.2lf\n", skewness);

	for (i = 0; i <= rat_poly->degree; i++)
		gmp_fprintf(fp, "R%u %Zd\n", i, rat_poly->coeff[i]);
	for (i = 0; i <= alg_poly->degree; i++)
		gmp_fprintf(fp, "A%u %Zd\n", i, alg_poly->coeff[i]);
	fclose(fp);
}

/*------------------------------------------------------------------*/
int32 find_poly(msieve_obj *obj, mpz_t n) {

	/* external entry point for NFS polynomial generation */

	poly_param_t params;
	poly_config_t config;
	uint32 degree;

	logprintf(obj, "commencing number field sieve polynomial selection\n");

	poly_config_init(&config);

	/* configure the selection process */

	get_poly_params(obj, n, &degree, &params);

	/* run the core polynomial finder */

	obj->flags |= MSIEVE_FLAG_SIEVING_IN_PROGRESS;
	find_poly_core(obj, n, &params, &config, degree);
	obj->flags &= ~MSIEVE_FLAG_SIEVING_IN_PROGRESS;

	/* save the best polynomial */

	logprintf(obj, "polynomial selection complete\n");
	if (config.heap[0]->rpoly.degree > 0) {
		write_poly(obj, n, &config.heap[0]->rpoly,
				&config.heap[0]->apoly,
				config.heap[0]->skewness);
	}
	poly_config_free(&config);

	return 0;
}
