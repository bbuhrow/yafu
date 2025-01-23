/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.

$Id: fb.c 978 2015-03-14 01:40:32Z jasonp_sf $
--------------------------------------------------------------------*/

#include <ms_common.h>
#include "gnfs.h"

#define REPORT_INTERVAL 20000

/*------------------------------------------------------------------*/
void create_factor_base(msieve_obj *obj, 
			factor_base_t *fb, 
			uint32 report_progress) {

	uint32 i;
	uint32 next_report = REPORT_INTERVAL;
	prime_sieve_t prime_sieve;
	fb_side_t *rfb = &fb->rfb;
	fb_side_t *afb = &fb->afb;
	uint32 max_prime_r = rfb->max_prime;
	uint32 max_prime_a = afb->max_prime;
	uint32 num_found_r = 0;
	uint32 num_found_a = 0;

	/* the arrays for the factor bases will grow dynamically */

	rfb->num_alloc = 10000;
	rfb->entries = (fb_entry_t *)xmalloc(rfb->num_alloc*sizeof(fb_entry_t));
	afb->num_alloc = 10000;
	afb->entries = (fb_entry_t *)xmalloc(afb->num_alloc*sizeof(fb_entry_t));

	init_prime_sieve(&prime_sieve, 0, MAX(max_prime_r, max_prime_a));

	if (report_progress)
		logprintf(obj, "generating factor base\n");

	while(1) {
		uint32 p = get_next_prime(&prime_sieve);
		uint32 roots[MAX_POLY_DEGREE+1];
		uint32 num_roots, high_coeff;

		/* if p does not exceed the relevant factor base
		   bound, compute the roots mod p of the factor
		   base polynomial, determine if p is a projective
		   root, and write out the list of roots discovered. */

		if (p <= max_prime_r) {
			num_roots = poly_get_zeros(roots, &rfb->poly, p, 
							&high_coeff, 0);
			if (high_coeff == 0) {
				roots[num_roots++] = p;
			}
			if (num_found_r + num_roots >= rfb->num_alloc) {
				rfb->num_alloc *= 2;
				rfb->entries = (fb_entry_t *)xrealloc(
							rfb->entries,
							rfb->num_alloc * 
							sizeof(fb_entry_t));
			}
			for (i = 0; i < num_roots; i++) {
				rfb->entries[num_found_r+i].p = p;
				rfb->entries[num_found_r+i].r = roots[i];
			}
			num_found_r += num_roots;
		}

		if (p <= max_prime_a) {
			num_roots = poly_get_zeros(roots, &afb->poly, p, 
							&high_coeff, 0);
			if (high_coeff == 0) {
				roots[num_roots++] = p;
			}
			if (num_found_a + num_roots >= afb->num_alloc) {
				afb->num_alloc *= 2;
				afb->entries = (fb_entry_t *)xrealloc(
							afb->entries,
							afb->num_alloc * 
							sizeof(fb_entry_t));
			}
			for (i = 0; i < num_roots; i++) {
				afb->entries[num_found_a+i].p = p;
				afb->entries[num_found_a+i].r = roots[i];
			}
			num_found_a += num_roots;
		}

		if (p > max_prime_r && p > max_prime_a)
			break;

		if (report_progress && (num_found_r > next_report ||
					num_found_a > next_report)) {
			fprintf(stderr, "factor base: found "
					"%u rational and %u "
					"algebraic entries\r", 
					num_found_r, num_found_a);
			fflush(stderr);
			next_report += REPORT_INTERVAL;
		}
	}

	free_prime_sieve(&prime_sieve);
	rfb->entries = (fb_entry_t *)xrealloc(rfb->entries,
				num_found_r * sizeof(fb_entry_t));
	afb->entries = (fb_entry_t *)xrealloc(afb->entries,
				num_found_a * sizeof(fb_entry_t));

	/* update the factor base with the real numbers
	   used to create it */

	rfb->num_entries = num_found_r;
	afb->num_entries = num_found_a;
	rfb->max_prime = rfb->entries[num_found_r - 1].p;
	afb->max_prime = afb->entries[num_found_a - 1].p;
	if (report_progress) {
		fprintf(stderr, "\n");
		logprintf(obj, "factor base complete:\n");
		logprintf(obj, "%u rational roots (max prime = %u)\n",
					rfb->num_entries, rfb->max_prime);
		logprintf(obj, "%u algebraic roots (max prime = %u)\n",
					afb->num_entries, afb->max_prime);
	}
}

/*------------------------------------------------------------------*/
void free_factor_base(factor_base_t *fb) {
	free(fb->rfb.entries);
	free(fb->afb.entries);
}

/*------------------------------------------------------------------*/
int32 read_factor_base(msieve_obj *obj, mpz_t n,
			sieve_param_t *params,
			factor_base_t *fb) {
	
	int32 status = 0;
	uint32 i;
	FILE *fp;
	char buf[LINE_BUF_SIZE];
	uint32 rfb_size = 0;
	uint32 afb_size = 0;

	uint32 p, r;
	char *tmp;
	char *next_field;

	memset(fb, 0, sizeof(factor_base_t));
	fb->rfb.max_prime = params->rfb_limit;
	fb->afb.max_prime = params->afb_limit;

	fp = fopen(obj->nfs_fbfile_name, "r");
	if (fp == NULL)
		return -1;

	mpz_poly_init(&fb->rfb.poly);
	mpz_poly_init(&fb->afb.poly);
	if (read_poly(obj, n, &fb->rfb.poly, 
				&fb->afb.poly, &params->skewness)) {
		printf("error reading NFS polynomials\n");
		exit(-1);
	}

	/* Look for lines that start with 'F' or 'S' and
	   then have a number. The numbers give various NFS
	   parameters, which will override the parameters
	   that were automatially generated by the main NFS
	   driver. This gives a way to specify your own parameters
	   for a factorization; if they don't exactly match
	   the parameters already in the factor base file, this
	   call will fail and (presumably) the file will get 
	   regenerated */

	fgets(buf, (int)sizeof(buf), fp);
	while (!feof(fp)) {
		uint32 value;

		/* stop reading configuration info when a
		   number is encountered */

		if (isdigit(buf[0]))
			break;

		if (buf[0] != 'F' && buf[0] != 'S') {
			fgets(buf, (int)sizeof(buf), fp);
			continue;
		}

		tmp = buf;
		while (!isdigit(*tmp) && *tmp != '-')
			tmp++;
		value = strtoul(tmp, NULL, 10);

		if (strstr(buf, "FRNUM"))
			rfb_size = value;
		else if (strstr(buf, "FRMAX"))
			fb->rfb.max_prime = params->rfb_limit = value;
		else if (strstr(buf, "FANUM"))
			afb_size = value;
		else if (strstr(buf, "FAMAX"))
			fb->afb.max_prime = params->afb_limit = value;
		else if (strstr(buf, "SRLPMAX"))
			params->rfb_lp_size = value;
		else if (strstr(buf, "SALPMAX"))
			params->afb_lp_size = value;
		else if (strstr(buf, "SMIN"))
			params->sieve_begin = strtoll(tmp, NULL, 10);
		else if (strstr(buf, "SMAX"))
			params->sieve_end = strtoll(tmp, NULL, 10);
		else if (strstr(buf, "SLINE")) {
			int64 sieve_size = strtoll(tmp, NULL, 10);
			params->sieve_begin = -sieve_size;
			params->sieve_end = sieve_size;
		}

		fgets(buf, (int)sizeof(buf), fp);
	}

	if (rfb_size == 0 || afb_size == 0) {
		status = -2;
		goto cleanup;
	}

	/* allocate room for that many entries */

	fb->rfb.num_entries = rfb_size;
	fb->rfb.num_alloc = rfb_size;
	fb->rfb.entries = (fb_entry_t *)xmalloc((rfb_size + MAX_POLY_DEGREE) * 
					sizeof(fb_entry_t));
	fb->afb.num_entries = afb_size;
	fb->afb.num_alloc = afb_size;
	fb->afb.entries = (fb_entry_t *)xmalloc((afb_size + MAX_POLY_DEGREE) * 
					sizeof(fb_entry_t));
	if (fb->rfb.entries == NULL || fb->afb.entries == NULL) {
		printf("error: failed to allocate factor base arrays\n");
		exit(-1);
	}
	rewind(fp);
	fgets(buf, (int)sizeof(buf), fp);
	while (!feof(fp) && !isdigit(buf[0])) {
		fgets(buf, (int)sizeof(buf), fp);
	}

	/* read the rational factor base. Every prime has
	   all of the roots for that prime listed after it.
	   Keep reading until the file runs out, or all 
	   entries are filled */

	i = 0;
	while (i < rfb_size) {
		if (feof(fp))
			break;

		tmp = buf;
		p = strtoul(tmp, &next_field, 10);
		tmp = next_field + 1;
		while (1) {
			r = strtoul(tmp, &next_field, 10);
			fb->rfb.entries[i].p = p;
			fb->rfb.entries[i].r = r;
			i++;
			if (*next_field == '\n' || *next_field == '\r')
				break;
			tmp = next_field + 1;
		}
		fgets(buf, (int)sizeof(buf), fp);
	}
	if (i != rfb_size || fb->rfb.entries[i-1].p != params->rfb_limit) {
		status = -3;
		goto cleanup;
	}

	/* repeat for the algebraic factor base */

	i = 0;
	while (i < afb_size) {
		if (feof(fp))
			break;

		tmp = buf;
		p = strtoul(tmp, &next_field, 10);
		tmp = next_field + 1;
		while (1) {
			r = strtoul(tmp, &next_field, 10);
			fb->afb.entries[i].p = p;
			fb->afb.entries[i].r = r;
			i++;
			if (*next_field == '\n' || *next_field == '\r')
				break;
			tmp = next_field + 1;
		}
		fgets(buf, (int)sizeof(buf), fp);
	}
	if (i != afb_size || fb->afb.entries[i-1].p != params->afb_limit) {
		status = -4;
		goto cleanup;
	}

	fb->rfb.max_prime = fb->rfb.entries[rfb_size - 1].p;
	fb->afb.max_prime = fb->afb.entries[afb_size - 1].p;
	logprintf(obj, "factor base loaded:\n");
	logprintf(obj, "%u rational ideals (max prime = %u)\n",
				fb->rfb.num_entries, fb->rfb.max_prime);
	logprintf(obj, "%u algebraic ideals (max prime = %u)\n",
				fb->afb.num_entries, fb->afb.max_prime);

cleanup:
	if (status != 0)
		free_factor_base(fb);
	fclose(fp);
	return status;
}

/*------------------------------------------------------------------*/
void write_factor_base(msieve_obj *obj, mpz_t n,
			sieve_param_t *params,
			factor_base_t *fb) {
	
	/* completely replace the current factor base file
	   with a version corresponding to the current factor
	   base and current polynomials */

	uint32 i, j;
	FILE *fp;
	fb_side_t *side;

	write_poly(obj, n, &fb->rfb.poly, &fb->afb.poly, params->skewness);

	fp = fopen(obj->nfs_fbfile_name, "a");
	if (fp == NULL) {
		printf("error; cannot open factor base file '%s'\n",
					obj->nfs_fbfile_name);
		exit(-1);
	}
	fprintf(fp, "\n");
	fprintf(fp, "FRNUM %u\n", fb->rfb.num_entries);
	fprintf(fp, "FRMAX %u\n", fb->rfb.max_prime);
	fprintf(fp, "FANUM %u\n", fb->afb.num_entries);
	fprintf(fp, "FAMAX %u\n", fb->afb.max_prime);

	fprintf(fp, "SRLPMAX %u\n", params->rfb_lp_size);
	fprintf(fp, "SALPMAX %u\n", params->afb_lp_size);
	fprintf(fp, "SLINE %" PRIu64 "\n", 
			(uint64)(params->sieve_end - params->sieve_begin) / 2);
	fprintf(fp, "\n");

	i = 0;
	side = &fb->rfb;
	while (i < side->num_entries) {
		fprintf(fp, "%u %u", side->entries[i].p, side->entries[i].r);
		for (j = i + 1; j < side->num_entries &&
				side->entries[j].p == side->entries[i].p; j++) {
			fprintf(fp, " %u", side->entries[j].r);
		}
		fprintf(fp, "\n");
		i = j;
	}

	i = 0;
	side = &fb->afb;
	while (i < side->num_entries) {
		fprintf(fp, "%u %u", side->entries[i].p, side->entries[i].r);
		for (j = i + 1; j < side->num_entries &&
				side->entries[j].p == side->entries[i].p; j++) {
			fprintf(fp, " %u", side->entries[j].r);
		}
		fprintf(fp, "\n");
		i = j;
	}

	fclose(fp);
}
