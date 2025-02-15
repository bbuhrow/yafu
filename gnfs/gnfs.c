/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: gnfs.c 1025 2018-08-19 02:20:28Z jasonp_sf $
--------------------------------------------------------------------*/

#include <ms_common.h>
#include "gnfs.h"
#include "savefile.h"

/* The first few parameter sets (except the sieve length) are
   from GGNFS. That doesn't mean they're right for this
   implementation though. The line length is a total guess */

static sieve_param_t prebuilt_params[] = {
	{MIN_NFS_BITS, 1800000, 1800000, 1<<26, 1<<26, 4000000, 0, 0, 1},
	{342,          2300000, 2300000, 1<<26, 1<<26, 8000000, 0, 0, 1},
	{352,          2500000, 2500000, 1<<26, 1<<26, 12000000, 0, 0, 1},
	{365,          3200000, 3200000, 1<<27, 1<<27, 16000000, 0, 0, 1},
	{372,          3500000, 3500000, 1<<27, 1<<27, 20000000, 0, 0, 1},
	{392,          4500000, 4500000, 1<<27, 1<<27, 24000000, 0, 0, 1},
	{405,          6000000, 6000000, 1<<27, 1<<27, 28000000, 0, 0, 1}, 
	{416,          8000000, 8000000, 1<<27, 1<<27, 32000000, 0, 0, 1}, 
	{429,          11000000, 11000000, 1<<27, 1<<27, 36000000, 0, 0, 1}, 
	{443,          12000000, 12000000, 1<<28, 1<<28, 40000000, 0, 0, 1}, 
	{453,          15000000, 15000000, 1<<28, 1<<28, 44000000, 0, 0, 1}, 
	{463,          18000000, 18000000, 1<<28, 1<<28, 48000000, 0, 0, 1}, 
	{477,          21000000, 21000000, 1<<29, 1<<29, 52000000, 0, 0, 1},
	{492,          24000000, 24000000, 1<<29, 1<<29, 56000000, 0, 0, 1},
	{506,          27000000, 27000000, 1<<29, 1<<29, 60000000, 0, 0, 1},
	{520,          30000000, 30000000, 1<<29, 1<<29, 64000000, 0, 0, 1},
};

static void get_sieve_params(uint32 bits, sieve_param_t *params);

static uint32 nfs_init_savefile(msieve_obj *obj, mpz_t n);

/*--------------------------------------------------------------------*/
uint32 factor_gnfs(msieve_obj *obj, mp_t *input_n,
			factor_list_t *factor_list) {

	int32 status;
	uint32 bits;
	sieve_param_t params;
	mpz_t n;
	mpz_poly_t rat_poly;
	mpz_poly_t alg_poly;
	uint32 relations_found = 0;
	uint32 max_relations = 0;
	uint32 factor_found = 0;

	mpz_init(n);
	mp2gmp(input_n, n);

	/* Calculate the factor base bound */

	bits = mpz_sizeinbase(n, 2);
	get_sieve_params(bits, &params);

	gmp_sprintf(obj->mp_sprintf_buf, "%Zd", n);
	logprintf(obj, "commencing number field sieve (%d-digit input)\n",
			strlen(obj->mp_sprintf_buf));

	/* generate or read in the NFS polynomials */

	mpz_poly_init(&rat_poly);
	mpz_poly_init(&alg_poly);
	status = read_poly(obj, n, &rat_poly, &alg_poly, &params.skewness);
	if (status != 0 && 
	   (obj->flags & (MSIEVE_FLAG_NFS_POLY1 | 
			  MSIEVE_FLAG_NFS_POLYSIZE |
			  MSIEVE_FLAG_NFS_POLYROOT))) {
		status = find_poly(obj, n);
		status = read_poly(obj, n, &rat_poly, 
					&alg_poly, &params.skewness);
	}
	if (status != 0) {
		if ((obj->flags & MSIEVE_FLAG_NFS_POLY1) &&
			(obj->flags & MSIEVE_FLAG_NFS_POLYSIZE) &&
			(obj->flags & MSIEVE_FLAG_NFS_POLYROOT))
		{
			// it's potentially an error only if we are actually 
			// doing all of the polygen phases.  Even so, we may
			// have just not found any polynomials with the 
			// provided parameters. 
			printf("No NFS polynomials found or read\n");
		}
		goto finished;
	}
	analyze_one_poly(obj, &rat_poly, &alg_poly, params.skewness);

	if ((obj->flags & MSIEVE_FLAG_STOP_SIEVING) ||
	    !(obj->flags & (MSIEVE_FLAG_NFS_SIEVE |
	    		    MSIEVE_FLAG_NFS_FILTER |
			    MSIEVE_FLAG_NFS_LA |
			    MSIEVE_FLAG_NFS_SQRT))) {
		goto finished;
	}

	/* if we're supposed to be sieving, 
	   initialize the savefile */

	if (obj->flags & (MSIEVE_FLAG_NFS_SIEVE |
			  MSIEVE_FLAG_NFS_FILTER)) {

		/* figure out how many relations to look for, and 
		   quit early if that many have already been found */

		relations_found = nfs_init_savefile(obj, n);
		if (obj->max_relations > 0) {
			max_relations = relations_found + 
						obj->max_relations;
		}
		else {
			/* if no guidance on this is available, make a
			   wild guess: a fixed fraction of the total number 
			   of large primes possible */

			max_relations = 0.8 * (params.rfb_lp_size /
					(log((double)params.rfb_lp_size) - 1) +
					params.afb_lp_size /
					(log((double)params.afb_lp_size) - 1));
		}
	}

	/* this is a little tricky: perform only sieving (and
	   quit when enough relations are found), only filtering
	   (and quit if it fails, or continue the postprocessing
	   if it succeeds), or alternate between sieving and
	   filtering until the filtering succeeds */

	while (1) {
		if (!(obj->flags & (MSIEVE_FLAG_NFS_SIEVE |
				    MSIEVE_FLAG_NFS_FILTER))) {
			break;
		}

		if (obj->flags & MSIEVE_FLAG_NFS_SIEVE) {
			savefile_open(&obj->savefile, SAVEFILE_APPEND);
			relations_found = do_line_sieving(obj, &params, n, 
							relations_found,
							max_relations);
			savefile_close(&obj->savefile);
			if (relations_found == 0)
				break;
			if (!(obj->flags & MSIEVE_FLAG_NFS_FILTER))
				goto finished;
		}

		if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)
			goto finished;

		if (obj->flags & MSIEVE_FLAG_NFS_FILTER) {

			max_relations = nfs_filter_relations(obj, n);
			if (max_relations == 0)
				break;
			logprintf(obj, "filtering wants %u more relations\n",
							max_relations);
			if (!(obj->flags & MSIEVE_FLAG_NFS_SIEVE))
				goto finished;
			max_relations += relations_found;
		}
	}

	if (obj->flags & MSIEVE_FLAG_NFS_LA)
		nfs_solve_linear_system(obj, n);
		
	if (obj->flags & MSIEVE_FLAG_NFS_SQRT)
		factor_found = nfs_find_factors(obj, n, factor_list);

finished:
	mpz_poly_free(&rat_poly);
	mpz_poly_free(&alg_poly);
	mpz_clear(n);
	return factor_found;
}

/*--------------------------------------------------------------------*/
static void get_sieve_params(uint32 bits, sieve_param_t *params) {

	sieve_param_t *low, *high;
	uint32 max_size;
	uint32 i, j, dist;
	int64 sieve_size;

	/* For small inputs, use the first set of precomputed
	   parameters */

	if (bits < prebuilt_params[0].bits) {
		*params = prebuilt_params[0];
		return;
	}

	/* bracket the input size between two table entries */

	max_size = sizeof(prebuilt_params) / sizeof(sieve_param_t);
	if (bits >= prebuilt_params[max_size - 1].bits) {
		*params = prebuilt_params[max_size - 1];
		params->sieve_begin = -params->sieve_size;
		params->sieve_end = params->sieve_size;
		return;
	}

	/* if the input is too large, just use the last table entry.
	   This means that the choice of parameters is increasingly
	   inappropriate as the input becomes larger, but there's no
	   guidance on what to do in this case anyway. */

	for (i = 0; i < max_size - 1; i++) {
		if (bits < prebuilt_params[i+1].bits)
			break;
	}

	/* Otherwise the parameters to use are a weighted average 
	   of the two table entries the input falls between */

	low = &prebuilt_params[i];
	high = &prebuilt_params[i+1];
	dist = high->bits - low->bits;
	i = bits - low->bits;
	j = high->bits - bits;

	params->bits = bits;
	params->rfb_limit = (uint32)(
			 ((double)low->rfb_limit * j +
			  (double)high->rfb_limit * i) / dist + 0.5);
	params->afb_limit = (uint32)(
			 ((double)low->afb_limit * j +
			  (double)high->afb_limit * i) / dist + 0.5);
	params->rfb_lp_size = (uint32)(
			 ((double)low->rfb_lp_size * j +
			  (double)high->rfb_lp_size * i) / dist + 0.5);
	params->afb_lp_size = (uint32)(
			 ((double)low->afb_lp_size * j +
			  (double)high->afb_lp_size * i) / dist + 0.5);
	sieve_size = (uint64)(
			 ((double)low->sieve_size * j +
			  (double)high->sieve_size * i) / dist + 0.5);
	params->sieve_begin = -sieve_size;
	params->sieve_end = sieve_size;
	params->skewness = 1;		/* suboptimal but safe default */
}

/*--------------------------------------------------------------------*/
void eval_poly(mpz_t res, int64 a, uint32 b, mpz_poly_t *poly) {

	/* Evaluate one polynomial at 'a' and 'b' */

	uint32 d = poly->degree;

	int64_2gmp(a, poly->tmp1);
	mpz_set_ui(poly->tmp2, b);
	mpz_set(res, poly->coeff[d]);

	while (--d) {
		mpz_mul(res, res, poly->tmp1);
		mpz_addmul(res, poly->coeff[d], poly->tmp2);
		mpz_mul_ui(poly->tmp2, poly->tmp2, b);
	}
	mpz_mul(res, res, poly->tmp1);
	mpz_addmul(res, poly->coeff[d], poly->tmp2);
}

/*------------------------------------------------------------------*/
static uint32 nfs_init_savefile(msieve_obj *obj, mpz_t n) {

	char buf[LINE_BUF_SIZE];
	uint32 relations_found = 0;
	savefile_t *savefile = &obj->savefile;
	uint32 update = 1;

	/* open the savefile; if the file already
	   exists and the first line contains n,
	   then we are restarting from a previous factorization */

	if (savefile_exists(savefile)) {
		savefile_open(savefile, SAVEFILE_READ);
		buf[0] = 0;
		savefile_read_line(buf, sizeof(buf), savefile);
		if (buf[0] == 'N') {
			mpz_t read_n;

			mpz_init(read_n);
			gmp_sscanf(buf + 2, "%Zd", read_n);
			if (mpz_cmp(n, read_n) == 0)
				update = 0;
			mpz_clear(read_n);
		}
		savefile_close(savefile);
	}

	if (update && (obj->flags & MSIEVE_FLAG_NFS_SIEVE)) {
		/* If sieving is about to add new relations, truncate 
		   the file and write the present n. I hope you backed 
		   up savefiles you wanted! */

		savefile_open(savefile, SAVEFILE_WRITE);
		gmp_sprintf(buf, "N %Zd\n", n);
		savefile_write_line(savefile, buf);
		savefile_flush(savefile);
		savefile_close(savefile);
	}
	else {
		/* we don't care how many relations are present,
		   so don't waste time counting them */
		return 0;
	}


	/* count the number of relations in the savefile;
	   do not verify any of them */

	savefile_open(savefile, SAVEFILE_READ);

	while (!savefile_eof(savefile)) {

		/* count relations; No checking of relations 
		   happens here */

		savefile_read_line(buf, sizeof(buf), savefile);

		while (!savefile_eof(savefile)) {
			if (isdigit(buf[0]) || buf[0] == '-')
				relations_found++;
			savefile_read_line(buf, sizeof(buf), savefile);
		}

		if (relations_found)
			logprintf(obj, "restarting with %u relations\n",
							relations_found);
	}
	savefile_close(savefile);
	return relations_found;
}
