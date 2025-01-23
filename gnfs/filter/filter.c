/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: filter.c 1023 2018-08-19 00:30:42Z jasonp_sf $
--------------------------------------------------------------------*/

#include "filter.h"

/*--------------------------------------------------------------------*/
static void find_fb_size(factor_base_t *fb, 
			uint32 limit_r, uint32 limit_a,
			uint32 *entries_r_out, uint32 *entries_a_out) {

	prime_sieve_t prime_sieve;
	uint32 entries_r = 0;
	uint32 entries_a = 0;

	/* If the filtering bounds are extremely large, just 
	   estimate the target as the number of primes less than 
	   the filtering bounds */

	if (limit_r > 20000000 && limit_a > 20000000) {
		*entries_r_out = 1.02 * limit_r / (log((double)limit_r) - 1);
		*entries_a_out = 1.02 * limit_a / (log((double)limit_a) - 1);
		return;
	}

	init_prime_sieve(&prime_sieve, 0, 
			MAX(limit_r, limit_a) + 1000);

	while (1) {
		uint32 p = get_next_prime(&prime_sieve);
		uint32 num_roots;
		uint32 high_coeff;
		uint32 roots[MAX_POLY_DEGREE + 1];

		if (p >= limit_r && p >= limit_a)
			break;

		if (p < limit_r) {
			num_roots = poly_get_zeros(roots, &fb->rfb.poly, p, 
							&high_coeff, 1);
			if (high_coeff == 0)
				num_roots++;
			entries_r += num_roots;
		}

		if (p < limit_a) {
			num_roots = poly_get_zeros(roots, &fb->afb.poly, p, 
							&high_coeff, 1);
			if (high_coeff == 0)
				num_roots++;
			entries_a += num_roots;
		}
	}

	free_prime_sieve(&prime_sieve);
	*entries_r_out = entries_r;
	*entries_a_out = entries_a;
}

/*--------------------------------------------------------------------*/
static uint32 check_excess(filter_t *filter) {

	/* give up if there is not enough excess 
	   to form a matrix */

	uint32 relations_needed = 0;

	if (filter->num_relations < filter->num_ideals ||
	    filter->num_relations - filter->num_ideals < 
	    				filter->target_excess) {
		relations_needed = 1000000;

		if (filter->num_relations > filter->num_ideals) {
			relations_needed = 3 * (filter->target_excess -
					(filter->num_relations - 
				 	 filter->num_ideals));
			relations_needed = MAX(relations_needed, 1000000);
		}
		free(filter->relation_array);
		filter->relation_array = NULL;
	}
	return relations_needed;
}

/*--------------------------------------------------------------------*/
static void dump_relation_numbers(msieve_obj *obj, filter_t *filter) {

	uint32 i;
	char buf[256];
	FILE *relation_fp;
	relation_ideal_t *r = filter->relation_array;

	sprintf(buf, "%s.d", obj->savefile.name);
	relation_fp = fopen(buf, "wb");
	if (relation_fp == NULL) {
		logprintf(obj, "error: reldump can't open out file\n");
		exit(-1);
	}

	for (i = 0; i < filter->num_relations; i++) {
		fwrite(&r->rel_index, (size_t)1, sizeof(uint32), relation_fp);
		r = next_relation_ptr(r);
	}
	fclose(relation_fp);
}

/*--------------------------------------------------------------------*/
static void set_filtering_bounds(msieve_obj *obj, factor_base_t *fb, 
			uint32 filtmin_r, uint32 filtmin_a, 
			uint32 *entries_r_out, uint32 *entries_a_out,
			uint32 num_relations, uint32 force_small, 
			filter_t *filter) {

	uint32 entries_r, entries_a;

	if (force_small) {
		if (num_relations < 2000000)
			filtmin_r = filtmin_a = 30000;
		else if (num_relations < 10000000) {
			filtmin_r = MIN(filtmin_r / 2, 100000);
			filtmin_a = MIN(filtmin_a / 2, 100000);
		}
		else {
			filtmin_r = MIN(filtmin_r / 2, 720000);
			filtmin_a = MIN(filtmin_a / 2, 720000);
		}
	}

	logprintf(obj, "reading ideals above %u\n", filtmin_r);
	find_fb_size(fb, filtmin_r, filtmin_a, &entries_r, &entries_a);
	filter->filtmin_r = filtmin_r;
	filter->filtmin_a = filtmin_a;
	filter->target_excess = entries_r + entries_a;

	*entries_r_out = entries_r;
	*entries_a_out = entries_a;
}

/*--------------------------------------------------------------------*/
/* the multiple of the amount of excess needed for
   merging to proceed */

#define FINAL_EXCESS_FRACTION 1.16

/* the default expected number of sparse nonzeros in the
   average matrix column (may be overriden if you know
   what you are doing) */

#define DEFAULT_TARGET_DENSITY 90.0

static uint32 do_merge(msieve_obj *obj, filter_t *filter, 
			merge_t *merge, double target_density) {

	uint32 relations_needed;
	uint32 extra_needed = filter->target_excess;

	/* make the clique removal more conservative by leaving 
	   some of the excess; this makes the merge phase easier. 
	   Note that the singleton removal probably threw away 
	   many large ideals that occur too often to be worth 
	   tracking, which forces the target matrix size to 
	   increase, so that target_excess is larger now */

	filter->target_excess *= FINAL_EXCESS_FRACTION;

	if ((relations_needed = check_excess(filter)) > 0)
		return relations_needed;

	/* build the matrix */

	merge->num_extra_relations = NUM_EXTRA_RELATIONS;

	merge->target_density = DEFAULT_TARGET_DENSITY;
	if (target_density != 0)
		merge->target_density = target_density;

	if (filter_make_relsets(obj, filter, merge, extra_needed) != 0)
		return 1000000;

	return 0;
}

/*--------------------------------------------------------------------*/
#define MAX_KEEP_WEIGHT 45

static uint32 do_partial_filtering(msieve_obj *obj, filter_t *filter,
				merge_t *merge, uint32 entries_r,
				uint32 entries_a, double target_density,
				uint32 max_weight) {

	uint32 relations_needed;
	uint32 num_relations = filter->num_relations;
	uint32 num_ideals = filter->num_ideals;

	if (filter->num_relations < 20000000 && max_weight < 25) {
		logprintf(obj, "raising initial max weight because there are few relations\n");
		max_weight = 25;
	}

	for (; max_weight < MAX_KEEP_WEIGHT; max_weight += 5) {

		filter->target_excess = entries_r + entries_a;
		filter->num_relations = num_relations;
		filter->num_ideals = num_ideals;

		filter_read_lp_file(obj, filter, max_weight);

		if ((relations_needed = do_merge(obj, filter, 
						merge, target_density)) > 0)
			return relations_needed;

		/* accept the collection of generated cycles 
		   if the matrix they form is dense enough or
		   max_weight has been incremented enough */

		if (merge->avg_cycle_weight > 63.0 ||
		    max_weight >= MAX_KEEP_WEIGHT - 5) 
			break;

		logprintf(obj, "matrix not dense enough, retrying\n");
		filter_free_relsets(merge);
	}

	return 0;
}

/*--------------------------------------------------------------------*/
uint32 nfs_filter_relations(msieve_obj *obj, mpz_t n) {

	filter_t filter;
	merge_t merge;
	uint32 filtmin_r, filtmin_a;
	uint32 entries_r, entries_a;
	uint32 num_relations;
	uint32 relations_needed = 0;
	factor_base_t fb;
	time_t wall_time = time(NULL);
	uint64 savefile_size = get_file_size(obj->savefile.name);
	uint64 ram_size = 0;
	uint32 max_relations = 0;
	uint32 filter_bound = 0;
	double target_density = 0;
	uint32 max_weight = 20;
	char lp_filename[256];

	logprintf(obj, "\n");
	logprintf(obj, "commencing relation filtering\n");

	/* parse arguments */

	if (obj->nfs_args != NULL) {

		const char *tmp;

		tmp = strstr(obj->nfs_args, "filter_mem_mb=");
		if (tmp != NULL) {
			ram_size = strtoull(tmp + 14, NULL, 10) << 20;
			logprintf(obj, "setting memory use to %.1f MB\n",
					(double)ram_size / 1048576);
		}

		tmp = strstr(obj->nfs_args, "filter_maxrels=");
		if (tmp != NULL) {
			max_relations = strtoul(tmp + 15, NULL, 10);
			logprintf(obj, "setting max relations to %u\n",
					max_relations);
		}

		tmp = strstr(obj->nfs_args, "filter_lpbound=");
		if (tmp != NULL) {
			filter_bound = strtoul(tmp + 15, NULL, 10);
			logprintf(obj, "setting large prime bound to %u\n",
					filter_bound);
		}

		tmp = strstr(obj->nfs_args, "target_density=");
		if (tmp != NULL) {
			target_density = strtod(tmp + 15, NULL);
			logprintf(obj, "setting target matrix density to %.1f\n",
					target_density);
		}

		tmp = strstr(obj->nfs_args, "max_weight=");
		if (tmp != NULL) {
			max_weight = strtoul(tmp + 11, NULL, 10);
			if (max_weight < MAX_KEEP_WEIGHT) {
				logprintf(obj, "setting initial max weight to %u\n",
						max_weight);
			}
			else {
#define str(s) #s
#define xstr(s) str(s)
				logprintf(obj, "initial max weight must be <= " xstr(MAX_KEEP_WEIGHT) "\n");
#undef xstr
#undef str
				exit(-1);
			}
		}

		/* old-style 'X,Y' format */

		tmp = strchr(obj->nfs_args, ',');
		if (tmp != NULL) {
			const char *tmp0 = tmp - 1;
			while (tmp0 > obj->nfs_args && isdigit(tmp0[-1]))
				tmp0--;
			max_relations = strtoul(tmp + 1, NULL, 10);
			filter_bound = strtoul(tmp0, NULL, 10);

			logprintf(obj, "setting max relations to %u\n",
					max_relations);
			logprintf(obj, "setting large prime bound to %u\n",
					filter_bound);
		}
	}

	if (ram_size == 0)
		ram_size = get_ram_size();

	memset(&filter, 0, sizeof(filter));
	memset(&merge, 0, sizeof(merge));
	memset(&fb, 0, sizeof(fb));
	mpz_poly_init(&fb.rfb.poly);
	mpz_poly_init(&fb.afb.poly);
	if (read_poly(obj, n, &fb.rfb.poly, &fb.afb.poly, NULL)) {
		printf("filtering failed to read polynomials\n");
		exit(-1);
	}
	logprintf(obj, "estimated available RAM is %.1lf MB\n", 
				(double)ram_size / 1048576);

	/* delete duplicate relations */

	filtmin_r = filtmin_a = nfs_purge_duplicates(obj, &fb, 
					max_relations, &num_relations);
	if (filter_bound > 0)
		filtmin_r = filtmin_a = filter_bound;

	/* set up the first disk-based pass; if the dataset is
	   "small", this will be the only such pass */

	set_filtering_bounds(obj, &fb, filtmin_r, filtmin_a,
				&entries_r, &entries_a, num_relations, 
				(uint32)(savefile_size < ram_size / 2), 
				&filter);

	/* separate out the large ideals and delete singletons
	   once they are all in memory. If the dataset is large,
	   first delete most of the singletons from the disk file */

	nfs_write_lp_file(obj, &fb, &filter, max_relations, 0);

	if (filter.lp_file_size > ram_size / 2) {
		filter_purge_lp_singletons(obj, &filter, ram_size);
#if 0
		/* also delete most of the cliques from the disk
		   file if it is still large, and there is a
		   great deal of excess relations */

		if (filter.lp_file_size > ram_size / 2) {
			filter_purge_lp_cliques(obj, &filter);
		}
#endif
	}
	filter_read_lp_file(obj, &filter, 0);

	if (savefile_size < ram_size / 2) {

		/* dataset is "small"; build the matrix immediately. 
		   Depending on how much memory the machine has, really
		   big datasets may get to do this */

		if ((relations_needed = do_merge(obj, &filter, 
						&merge, target_density)) > 0)
			goto finished;
	}
	else {  
		/* dataset is "large", perform multiple singleton passes.

		   The first pass used a large bound; the second 
		   filtering bound is much smaller. To allow reuse of 
		   previous results, the second bound is used 
		   during the rest of the filtering */

		dump_relation_numbers(obj, &filter);
		set_filtering_bounds(obj, &fb, filtmin_r, filtmin_a,
					&entries_r, &entries_a, 
					filter.num_relations, 1, &filter);

		free(filter.relation_array);
		filter.relation_array = NULL;

		nfs_write_lp_file(obj, &fb, &filter, max_relations, 1);

		if (filter.lp_file_size < ram_size / 2) {

			/* dataset is small enough for filtering to 
			   complete in one pass */

			filter_read_lp_file(obj, &filter, 0);
			if ((relations_needed = do_merge(obj, &filter, 
						&merge, target_density)) > 0) {
				goto finished;
			}
		}
		else {
			/* dataset is so large that even the pruned version
			   cannot fit comfortably in memory. We have to put
			   up with reading only the ideals that occur in the
			   fewest relations, forming the matrix, and then
			   determining whether the matrix incorporates enough
			   of the dataset so that the matrix will work */

			if ((relations_needed = do_partial_filtering(obj,
						&filter, &merge, entries_r,
						entries_a, target_density,
						max_weight)) > 0) {
				goto finished;
			}
		}
	}

	/* filtering succeeded; delete the LP file */

	sprintf(lp_filename, "%s.lp", obj->savefile.name);
	remove(lp_filename);

	/* optimize and then save the collection of relation-sets */

	filter_postproc_relsets(obj, &merge);
	filter_dump_relsets(obj, &merge);
	filter_free_relsets(&merge);
	wall_time = time(NULL) - wall_time;
	logprintf(obj, "RelProcTime: %u\n", (uint32)wall_time);
finished:
	mpz_poly_free(&fb.rfb.poly);
	mpz_poly_free(&fb.afb.poly);
	return relations_needed;
}
