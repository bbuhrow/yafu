/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.

$Id: sieve_line.c 982 2015-03-19 16:51:42Z jasonp_sf $
--------------------------------------------------------------------*/

#include "sieve.h"

/* sieving takes place in L1-cache-size blocks */
#define BLOCK_SIZE 65536
#define BLOCK_HASH(x) ((uint64)(x) / BLOCK_SIZE)

/* how often the log targets for the algebraic sieve are
   recalculated. Currently the rational sieve is only
   updated once per block */
#define A_SAMPLE_RATE 512

/* structure used to sieve with one factor base prime */
typedef struct {
	uint32 p;
	uint32 r;
	uint32 next;
	uint16 offset;
	uint8 logp;
	uint8 skip;
} fb_sieve_entry_t;

/* used for sieving by powers of factor base primes */
typedef struct {
	uint16 p;	/* the power of a factor base prime */
	uint16 r;	 /* sieve root (does not change during sieving) */
	uint16 offset; /* the offset of the next update */
	uint8 logp;	/* scaled logarithm of 'prime' */
	uint8 orig_p;	/* the original factor base prime */
	uint8 skip;	/* if nonzero, do not sieve with this entry */
} fb_power_t;

#define FB_POWER_LIMIT 40

/* packed version of fb_entry_t, used for large FB primes */
typedef struct {
	uint32 p;		/* factor base prime */
	uint16 offset;	 /* the offset of the current update */
	uint8 logp;	/* scaled logarithm of 'prime' */
} packed_fb_t;

/* structure for projective roots */
typedef struct {
	uint32 p;			 /* factor base prime */
} projective_t;
	
/* structure used for resieving. We buffer at most
   MAX_RESIEVE_ENTRIES sieve values at a time, and 
   each sieve value uses up two resieve_t structures */

#define MAX_RESIEVE_ENTRIES 255
#define RESIEVE_INVALID 0xff
#define NUM_RESIEVE_FACTORS 19

typedef struct {
	uint16 offset;
	uint16 num_factors;
	uint32 factors[NUM_RESIEVE_FACTORS];
} resieve_t;
	
/* main structure controlling rational or algebraic sieving */
typedef struct {
	mpz_poly_t *poly;	/* the sieve polynomial */

	uint8 *sieve_block;	/* piece of sieve interval (one block worth) */
	uint32 *bucket_list; /* head of linked list of updates (per bucket) */

	uint32 curr_num_lp;	/* the number of large primes for this side */
	uint32 LP1_max;		/* single large prime max bound */

	uint32 cutoff2;	 	/* fudge factor for trial factoring */
	uint32 scaled_cutoff2;	/* fudge factor for trial factoring */
	mpz_t LP2_min;		/* double large prime min bound */
	mpz_t LP2_max;		/* double large prime max bound */

	uint32 cutoff3;	 	/* fudge factor for trial factoring */
	uint32 scaled_cutoff3;	/* fudge factor for trial factoring */
	mpz_t LP3_min;		/* triple large prime min bound */
	mpz_t LP3_max;		/* triple large prime max bound */

	mpz_t res;     /* temporaries */
	mpz_t tmp2;
	mpz_t tmp3;

	double base;	/* base used for logs of sieve values and primes */
	double log_base;	/* logarithm of the above */

	uint32 med_fb_size; /* number of entries with primes below BLOCK_SIZE */
	uint32 small_fb_size;/* number of entries that are never resieved */
	uint32 fb_size;      /* total number of factor base entries */

	projective_t proj_entries[20]; /* projective roots */
	uint32 proj_fb_size;		 /* number of projective roots */
	uint8 proj_bias;	 /* sum of logs of projective roots that divide
					the current b value */ 

	fb_sieve_entry_t *fb_entries; /* the factor base */

	fb_power_t *fb_powers;   /* used for powers of factor base primes */
	uint32 powers_fb_size;   /* the number of such entries */
	uint32 powers_fb_alloc;  /* space allocated for list of powers */

	packed_fb_t *update_list; /* list of updates to current block */
	uint32 num_updates_alloc; /* number of allocated updates */
	uint32 num_updates;	  /* the current number of updates */
} sieve_t;
	
/* main sieving structure */
typedef struct {
	msieve_obj *obj;
	factor_base_t *fb;

	int64 min_a;
	int64 max_a;
	uint32 min_b;
	uint32 max_b;

	uint32 num_buckets;
	sieve_t sieve_rfb;
	sieve_t sieve_afb;

	resieve_t *resieve_array;

	relation_batch_t relation_batch;

} sieve_job_t;

static void init_one_fb(msieve_obj *obj, 
			fb_side_t *fb, sieve_t *out_fb, 
			uint32 num_buckets, uint32 lp_size,
			char *string);

static void free_one_sieve_fb(sieve_t *out_fb);

static uint32 do_one_line(sieve_job_t *job, uint32 b_offset);

static void init_one_sieve(sieve_t *out_fb,
			uint32 num_buckets, 
			int64 min_a, int64 max_a, 
			uint32 min_b, uint32 b_offset);

static void fill_one_block(sieve_t *sieve_fb, 
			uint32 block, uint32 num_buckets);

static uint32 do_factoring(sieve_job_t *job,
			int64 block_start, uint32 b,
			mpz_t log_scratch);

static uint32 do_resieve(sieve_job_t *job, uint32 low_offset,
			uint32 high_offset, uint32 num_resieve, 
			int64 block_start, uint32 b);

static void do_one_resieve(sieve_t *fb, resieve_t *resieve_array, 
			uint8 *hashtable, uint32 low_offset, 
			uint32 high_offset);

static uint32 do_one_factoring(sieve_job_t *job, 
				resieve_t *sieve_value,
				int64 block_start, uint32 b);

static uint32 do_one_tf(sieve_t *sieve_fb, resieve_t *resieve, 
			uint32 *factors, uint32 *num_factors_out, 
			uint32 offset, uint32 b);

/*------------------------------------------------------------------*/
uint32 do_line_sieving(msieve_obj *obj, sieve_param_t *params, mpz_t n,
			uint32 relations_found, uint32 max_relations) {

	uint32 i;
	sieve_job_t job;
	factor_base_t fb;
	const char *lower_limit = NULL;
	const char *upper_limit = NULL;

	if (relations_found >= max_relations)
		return relations_found;

	/* parse arguments */

	if (obj->nfs_args != NULL) {
		upper_limit = strchr(obj->nfs_args, ',');
		if (upper_limit != NULL) {
			lower_limit = upper_limit - 1;
			while (lower_limit > obj->nfs_args &&
				isdigit(lower_limit[-1])) {
				lower_limit--;
			}
			upper_limit++;
		}
	}

	/* generate or read the factor base */

	if (read_factor_base(obj, n, params, &fb)) {
		create_factor_base(obj, &fb, 1);
		write_factor_base(obj, n, params, &fb);
	}

	/* initialize sieve struct */

	memset(&job, 0, sizeof(job));
	job.obj = obj;
	job.fb = &fb;
	job.min_a = params->sieve_begin;
	job.max_a = params->sieve_end;
	job.min_b = 1;
	job.max_b = 0xffffffff;     /* default is to sieve forever */

	/* set user-specified limits, if any */

	if (lower_limit != NULL && upper_limit != NULL) {
		job.min_b = strtoul(lower_limit, NULL, 10);
		job.max_b = strtoul(upper_limit, NULL, 10);
		if (job.min_b > job.max_b) {
			printf("lower bound on b must be <= upper bound\n");
			return 0;
		}
	}
	else if ((lower_limit == NULL && upper_limit != NULL) ||
		 (lower_limit != NULL && upper_limit == NULL) ) {
		printf("lower/upper bounds on b must both be specified\n");
		return 0;
	}
	else {
		/* see if there is any guidance on the first b
		   value to use, i.e. if a run is being resumed */

		i = read_last_line(obj, n);
		if (i > 0)
			job.min_b = i;
	}

	/* the lower sieve limit is assumed to be even */
	if (job.min_a & 1)
		job.min_a--;

	/* perform all one-time initialization. Most of the
	   factor base will use a bucket sort for cache efficiency,
	   with the number of buckets chosen to cover the maximum
	   size of any factor base prime */

	job.num_buckets = MAX(BLOCK_HASH(fb.rfb.max_prime) + 1,
			      BLOCK_HASH(fb.afb.max_prime) + 1);
	
	logprintf(obj, "a range: [%" PRId64 ", %" PRId64 "]\n", 
					job.min_a, job.max_a);
	logprintf(obj, "b range: [%u, %u]\n", job.min_b, job.max_b);
	logprintf(obj, "number of hash buckets: %u\n", job.num_buckets);
	logprintf(obj, "sieve block size: %u\n", BLOCK_SIZE);
	logprintf(obj, "\n");

	init_one_fb(obj, &fb.rfb, &job.sieve_rfb, job.num_buckets, 
			params->rfb_lp_size, "RFB");
	init_one_fb(obj, &fb.afb, &job.sieve_afb, job.num_buckets, 
			params->afb_lp_size, "AFB");

	/* every sieve value needs two resieve_t entries, the
	   first for resieved algebraic factors and the second
	   for resieved rational factors */

	job.resieve_array = (resieve_t *)xmalloc(2 * MAX_RESIEVE_ENTRIES *
						sizeof(resieve_t));

	/* initialize the structures for batch factoring of relations.
	   We use batch factoring to split the parts of relations
	   containing large primes, and to do that we have to multiply
	   together all the primes from the factor base bound to 
	   somewhere below the large prime bound */

	i = MAX(job.sieve_rfb.LP1_max, job.sieve_afb.LP1_max);
	i = MIN(3 << 27, i / 4);
	
	relation_batch_init(job.obj, &job.relation_batch,
			MIN(fb.rfb.max_prime, fb.afb.max_prime),
			i,
			job.sieve_rfb.LP1_max,
			job.sieve_afb.LP1_max,
			&obj->savefile, print_relation);

	if (obj->flags & (MSIEVE_FLAG_USE_LOGFILE |
	    		   MSIEVE_FLAG_LOG_TO_STDOUT)) {
		fprintf(stderr, "\nsieving in progress "
				"(press Ctrl-C to pause)\n");
	}

	/* sieve one line at a time. An easy way to do this
	   in parallel with k machines available is to start 
	   with i being some number mod k and to do every k_th
	   sieve line. Each sieve line is initialized individually,
	   so this is easy to implement */

	obj->flags |= MSIEVE_FLAG_SIEVING_IN_PROGRESS;
	for (i = 0; i <= job.max_b - job.min_b; i++) {

		relations_found += do_one_line(&job, i);

		/* print a progress message if appropriate */

		if ((obj->flags & (MSIEVE_FLAG_USE_LOGFILE |
	    		   	  MSIEVE_FLAG_LOG_TO_STDOUT))) {
			fprintf(stderr, "b = %u, %u complete / "
				"%u batched relations (need %u)\r", 
				job.min_b + i, relations_found, 
				job.relation_batch.num_relations,
				max_relations);
			fflush(stderr);
		}

		/* quit if appropriate */

		if (relations_found >= max_relations ||
		    (obj->flags & MSIEVE_FLAG_STOP_SIEVING) ||
		    i == job.max_b - job.min_b) {

			/* finish up any batch factoring that's left */

			if (job.relation_batch.num_relations > 0) {
				relations_found += relation_batch_run(
							&job.relation_batch);
			}
						
			if (obj->flags & (MSIEVE_FLAG_USE_LOGFILE |
				     	MSIEVE_FLAG_LOG_TO_STDOUT))
				fprintf(stderr, "\n");

			logprintf(obj, "completed b = %u, "
				"found %u relations\n", 
				job.min_b + i, relations_found);
			break;
		}
	}
	obj->flags &= ~MSIEVE_FLAG_SIEVING_IN_PROGRESS;

	relation_batch_free(&job.relation_batch);
	savefile_flush(&obj->savefile);
	free_one_sieve_fb(&job.sieve_rfb);
	free_one_sieve_fb(&job.sieve_afb);
	free_factor_base(&fb);
	free(job.resieve_array);
	write_last_line(obj, n, job.min_b + i + 1);
	return relations_found;
}

/*------------------------------------------------------------------*/
static void init_one_fb(msieve_obj *obj, fb_side_t *fb, sieve_t *out_fb, 
			uint32 num_buckets, uint32 lp_size, char *string) {

	/* Set up all of the permanent parameters in 
	   one factor base */

	uint32 i, j, k;
	uint32 largest_p = fb->max_prime;
	uint32 high_bound;
	fb_entry_t *aliased;

	out_fb->poly = &fb->poly;

	/* make the version of the factor base used for sieving
	   alias right on top of the input factor base */

	out_fb->fb_size = fb->num_entries;
	out_fb->fb_entries = (fb_sieve_entry_t *)xrealloc(fb->entries,
						fb->num_entries * 
						sizeof(fb_sieve_entry_t));
	aliased = (fb_entry_t *)(out_fb->fb_entries);
	for (i = fb->num_entries - 1; (int32)i >= 0; i--) {
		out_fb->fb_entries[i].p = aliased[i].p;
		out_fb->fb_entries[i].r = aliased[i].r;
	}
	fb->entries = NULL;

	/* allocate the rest of the data structures */

	mpz_init(out_fb->res);
	mpz_init(out_fb->tmp2);
	mpz_init(out_fb->tmp3);
	out_fb->sieve_block = (uint8 *)xmalloc(BLOCK_SIZE * sizeof(uint8));
	out_fb->bucket_list = (uint32 *)xmalloc(num_buckets * sizeof(uint32));
	out_fb->num_updates_alloc = 10000;
	out_fb->update_list = (packed_fb_t *)xmalloc(
					out_fb->num_updates_alloc * 
					sizeof(packed_fb_t));

	/* Calculate the large prime cutoffs */

	out_fb->LP1_max = lp_size;

	mpz_init_set_ui(out_fb->LP2_min, largest_p);
	mpz_mul_ui(out_fb->LP2_min, out_fb->LP2_min, largest_p);

	mpz_init_set_ui(out_fb->LP2_max, lp_size);
	mpz_mul_ui(out_fb->LP2_max, out_fb->LP2_max, lp_size);
	mpz_tdiv_q_2exp(out_fb->LP2_max, out_fb->LP2_max, 3);

	/* we accept triple-large-prime cofactors if they are
	   somewhat smaller than lp_size^3 and somewhat larger
	   than (factor base limit)^3 */

	mpz_init_set_ui(out_fb->LP3_min, largest_p);
	mpz_pow_ui(out_fb->LP3_min, out_fb->LP3_min, 3);
	mpz_mul_ui(out_fb->LP3_min, out_fb->LP3_min, 5);

	mpz_init_set_ui(out_fb->LP3_max, lp_size);
	mpz_pow_ui(out_fb->LP3_max, out_fb->LP3_max, 3);
	mpz_tdiv_q_2exp(out_fb->LP3_max, out_fb->LP3_max, 3);

	/* Determine the cutoff for log values. Sieve offsets whose
	   unfactored log value is smaller than this will have trial 
	   factoring attempted. We want the cutoff to be a few bits larger
	   than the cutoff for large primes, to give relations a little
	   leeway for roundoff errors compared to the true log value */

	high_bound = (uint32)(2.6 * log((double)largest_p) / M_LN2 + 0.5);
	out_fb->cutoff2 = 4 + mpz_sizeinbase(out_fb->LP2_max, 2);
	out_fb->cutoff2 = MAX(out_fb->cutoff2, high_bound);

	high_bound = (uint32)(3.9 * log((double)largest_p) / M_LN2 + 0.5);
	out_fb->cutoff3 = 4 + mpz_sizeinbase(out_fb->LP3_max, 2);
	out_fb->cutoff3 = MAX(out_fb->cutoff3, high_bound);

	/* Make factor base entries out of the powers of all the factor
	   base primes up to FB_POWER_LIMIT */

	out_fb->powers_fb_size = 0;
	out_fb->powers_fb_alloc = 100;
	out_fb->fb_powers = (fb_power_t *)xmalloc(out_fb->powers_fb_alloc *
						sizeof(fb_power_t));

	for (i = 0; i < out_fb->fb_size; i++) {
		uint32 p = out_fb->fb_entries[i].p;
		uint32 power;
		if (p > FB_POWER_LIMIT)
			break;

		if (i > 0 && p == out_fb->fb_entries[i-1].p)
			continue;

		for (power = p * p; power < 1400; power *= p) {
			for (j = 0; j < power; j++) {
				eval_poly(out_fb->res, (int64)j, 
						1, out_fb->poly);
				if (mpz_fdiv_ui(out_fb->res, power) != 0)
					continue;

				if (out_fb->powers_fb_size >= 
					    out_fb->powers_fb_alloc) {
					out_fb->powers_fb_alloc *= 2;
					out_fb->fb_powers = (fb_power_t *)
						xrealloc(out_fb->fb_powers,
							out_fb->powers_fb_alloc 
							* sizeof(fb_power_t));
				}
				k = out_fb->powers_fb_size++;
				out_fb->fb_powers[k].p = power;
				out_fb->fb_powers[k].r = j;
				out_fb->fb_powers[k].orig_p = p;
			}
		}
		if (power == p * p) {	
		        /* all future powers will be too large */
			break;
		}
	}

	/* determine first factor base entry that is resieved */

	out_fb->small_fb_size = out_fb->fb_size;
	for (; i < out_fb->fb_size; i++) {
		if (out_fb->fb_entries[i].p > 1000) {
			out_fb->small_fb_size = i;
			break;
		}
	}

	/* determine the point at which the switch from ordinary
	   sieving to bucket sieving happens. The crossover point
	   is the factor base entry that is guaranteed to have only
	   one sieve update, at most, per block */

	out_fb->med_fb_size = out_fb->fb_size;
	for (; i < out_fb->fb_size; i++) {
		if (out_fb->fb_entries[i].p > BLOCK_SIZE) {
			out_fb->med_fb_size = i;
			break;
		}
	}

	/* separate out the projective roots */

	for (i = 0; i < out_fb->fb_size; i++) {
		fb_sieve_entry_t *entry = out_fb->fb_entries + i;
		projective_t *proj_entry;
		uint32 p = entry->p;

		if (entry->r != p)
			continue;

		k = out_fb->proj_fb_size++;
		proj_entry = out_fb->proj_entries + k;
		proj_entry->p = p;
	}

	/* log the choices above */

	logprintf(obj, "maximum %s prime: %u\n", string, largest_p);
	logprintf(obj, "%s entries: %u\n", string, out_fb->fb_size);
	logprintf(obj, "medium %s entries: %u\n", string, out_fb->med_fb_size);
	logprintf(obj, "resieved %s entries: %u\n", string, 
				out_fb->med_fb_size - out_fb->small_fb_size);
	logprintf(obj, "small %s prime powers: %u\n", string, 
					out_fb->powers_fb_size);
	logprintf(obj, "projective %s roots: %u\n", string, 
					out_fb->proj_fb_size);
	logprintf(obj, "%s trial factoring cutoff: %u or %u bits\n", string, 
					out_fb->cutoff2, out_fb->cutoff3);
	logprintf(obj, "single large prime %s range: %u - %u bits\n", string,
			(uint32)(log((double)largest_p)/M_LN2 + 0.5),
			(uint32)(log((double)out_fb->LP1_max)/M_LN2 + 0.5));
	logprintf(obj, "double large prime %s range: %u - %u bits\n", string,
			mpz_sizeinbase(out_fb->LP2_min, 2), 
			mpz_sizeinbase(out_fb->LP2_max, 2));
	logprintf(obj, "triple large prime %s range: %u - %u bits\n", string,
			mpz_sizeinbase(out_fb->LP3_min, 2), 
			mpz_sizeinbase(out_fb->LP3_max, 2));
	logprintf(obj, "\n");
}

/*------------------------------------------------------------------*/
static void free_one_sieve_fb(sieve_t *out_fb) {

	out_fb->fb_entries = NULL;
	mpz_clear(out_fb->res);
	mpz_clear(out_fb->tmp2);
	mpz_clear(out_fb->tmp3);
	mpz_clear(out_fb->LP2_min);
	mpz_clear(out_fb->LP2_max);
	mpz_clear(out_fb->LP3_min);
	mpz_clear(out_fb->LP3_max);
	free(out_fb->sieve_block);
	free(out_fb->bucket_list);
	free(out_fb->update_list);
	free(out_fb->fb_powers);
}

/*------------------------------------------------------------------*/
static uint32 do_one_line(sieve_job_t *job, uint32 b_offset)
{
	/* Perform all the sieving for one value of b */

	msieve_obj *obj = job->obj;
	uint32 rels = 0;
	uint32 i;
	int64 min_a = job->min_a;
	int64 max_a = job->max_a;
	uint32 min_b = job->min_b;
	uint32 b = min_b + b_offset;
	int64 block_base = min_a;
	uint32 num_buckets = job->num_buckets;
	mpz_t log_scratch;

	/* finish off the factor base initialization and fill
	   in the initial values of all the sieve updates */

	init_one_sieve(&job->sieve_rfb, num_buckets, 
			min_a, max_a, min_b, b_offset);
	init_one_sieve(&job->sieve_afb, num_buckets, 
			min_a, max_a, min_b, b_offset);
	mpz_init(log_scratch);

	while (min_a < max_a) {

		/* for each sieve block */

		for (i = 0; i < num_buckets; i++) {

			/* fill up the sieve block */

			fill_one_block(&job->sieve_rfb, i, num_buckets);
			fill_one_block(&job->sieve_afb, i, num_buckets);

			/* scan it for values to trial factor */

			rels += do_factoring(job, block_base, b, log_scratch);
			if (b % 2 == 0)
				block_base += 2 * BLOCK_SIZE;
			else
				block_base += BLOCK_SIZE;

			if (block_base >= max_a)
				break;
		}

		if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)
			break;
		min_a = block_base;
	}

	mpz_clear(log_scratch);
	return rels;
}

#define LOG_UPDATE_RATE 100

/*------------------------------------------------------------------*/
static void init_one_sieve(sieve_t *out_fb,
			uint32 num_buckets, 
			int64 min_a, int64 max_a, 
			uint32 min_b, uint32 b_offset) {

	/* initialize all of the sieve updates */

	uint32 i;
	uint32 b = min_b + b_offset;
	uint32 common;

	if (b_offset % LOG_UPDATE_RATE == 0) {
		/* determine smallest base of logarithms of sieve
		   values, such that the log in this base of the
		   largest possible sieve value fits in a byte */

		double log_of_base;

		out_fb->base = get_log_base(out_fb->poly, min_a, max_a, 
					b + LOG_UPDATE_RATE - 1);
		log_of_base = out_fb->log_base = log(out_fb->base);

		/* convert the trial factoring cutoffs to this base */

		out_fb->scaled_cutoff2 = (uint32)((double)out_fb->cutoff2 * 
						M_LN2 / log_of_base + 0.5);
		out_fb->scaled_cutoff3 = (uint32)((double)out_fb->cutoff3 * 
						M_LN2 / log_of_base + 0.5);

		/* compute logarithms of the factor base primes. Note that
		   for prime powers, the log value is always that of the
		   original factor base prime and not the actual power. This
		   is because sieve updates caused by powers of primes will
		   stack on top of each other */

		for (i = 0; i < out_fb->fb_size; i++) {
			fb_sieve_entry_t *entry = out_fb->fb_entries + i;
			entry->logp = (uint8)fplog(entry->p, log_of_base);
		}
		for (i = 0; i < out_fb->powers_fb_size; i++) {
			fb_power_t *entry = out_fb->fb_powers + i;
			entry->logp = (uint8)fplog(entry->orig_p, log_of_base);
		}
	}

	/* clear the linked lists of sieve updates */

	memset(out_fb->bucket_list, 0, num_buckets * sizeof(uint32));

	/* initialize the factor base, and put the sieve update
	   associated with each factor base prime into the correct
	   bucket. This amount to an in-place radix sort on the
	   collection of updates */

	for (i = 0; i < out_fb->fb_size; i++) {
		fb_sieve_entry_t *entry = out_fb->fb_entries + i;
		uint32 p = entry->p;
		uint32 br;
		int64 rem;

		/* compute the sieve update corresponding to
		   this factor base prime. The initial sieve value
		   is the first value of b*r + k*p that exceeds min_a.
		   
		   Do not sieve with this entry if p divides b or
		   if the entry is a projective root (i.e. p == r) */

		entry->skip = 0;
		br = mp_modmul_1(entry->r, b, p);
		if (br == 0 && (entry->r == p || b % p == 0)) {
			entry->skip = 1;
			continue;
		}

		rem = (int64)p - min_a % (int64)p;
		if (rem >= p)
			rem -= p;

		rem += br;
		if (rem >= p)
			rem -= p;

		/* if b is even, the sieve interval is compressed to
		   contain only odd values */

		if (b % 2 == 0) {
			if (rem % 2 == 0)
				rem += p;
			rem = rem / 2;
		}

		/* compute the sieve update corresponding to
		   this factor base prime */

		entry->offset = (uint16)(rem % BLOCK_SIZE);

		/* if p exceeds the block size, add it to the 
		   linked list corresponding to the bucket in which the 
		   update falls */

		if (p > BLOCK_SIZE) {
			uint32 bucket = BLOCK_HASH(rem);
			entry->next = out_fb->bucket_list[bucket];
			out_fb->bucket_list[bucket] = i;
		}
	}

	/* repeat for the powers of small factor base primes */

	for (i = 0; i < out_fb->powers_fb_size; i++) {
		fb_power_t *entry = out_fb->fb_powers + i;
		uint32 p = entry->p;
		uint32 br;
		int64 rem;

		entry->skip = 0;
		br = mp_modmul_1(entry->r, b, p);
		if (br == 0 && b % entry->orig_p == 0) {
			entry->skip = 1;
			continue;
		}

		rem = (int64)p - min_a % (int64)p;
		if (rem >= p)
			rem -= p;
		rem += br;
		if (rem >= p)
			rem -= p;

		if (b % 2 == 0) {
			if (rem % 2 == 0)
				rem += p;
			rem = rem / 2;
		}

		entry->offset = (uint16)(rem & (BLOCK_SIZE - 1));
	}

	/* find the projective roots that divide b, and compute
	   their combined log value. This is added to the cutoff
	   for the entire sieve later */

	common = mpz_fdiv_ui(out_fb->poly->coeff[out_fb->poly->degree], b);
	common = mp_gcd_1(common, b);
	out_fb->proj_bias = fplog(common, out_fb->log_base);
}

/*------------------------------------------------------------------*/
static void fill_one_block(sieve_t *sieve_fb, 
			uint32 block, uint32 num_buckets) {

	/* combine all of the sieve updates for one factor base */

	uint32 i;
	fb_sieve_entry_t *factor_base = sieve_fb->fb_entries;
	uint32 *bucket_list = sieve_fb->bucket_list;
	uint32 med_fb_size = sieve_fb->med_fb_size;
	uint32 powers_fb_size = sieve_fb->powers_fb_size;
	uint8 *sieve_block = sieve_fb->sieve_block;
	packed_fb_t *update_list = sieve_fb->update_list;
	uint32 num_updates_alloc = sieve_fb->num_updates_alloc;
	uint32 num_updates;
	uint32 fb_offset;

	/* first, pack together the sieve updates corresponding
	   to the large factor base primes. Begin by getting the
	   head of the linked list of updates. Then set the 
	   list of updates to empty, so updates that wrap around 
	   the end of the list of buckets have somewhere to go */

	fb_offset = bucket_list[block];
	bucket_list[block] = 0;
	i = 0;

	while (fb_offset != 0) {
		fb_sieve_entry_t *entry = factor_base + fb_offset;
		packed_fb_t update;
		uint32 tmp_offset;
		uint32 new_bucket;
		uint64 new_offset;

		PREFETCH(factor_base + entry->next);
		update.p = entry->p;
		update.offset = entry->offset;
		update.logp = entry->logp;

		if (i >= num_updates_alloc) {
			num_updates_alloc += 5000;
			sieve_fb->num_updates_alloc = num_updates_alloc;
			update_list = (packed_fb_t *)xrealloc(update_list, 
							num_updates_alloc * 
							sizeof(packed_fb_t));
			sieve_fb->update_list = update_list;
		}
		update_list[i++] = update;

		/* Since the act of packing things will pull the 
		   factor base entry into cache anyway, also fold 
		   in the computation of the next update for this 
		   factor base prime */

		new_offset = (uint64)block * BLOCK_SIZE + 
					entry->offset + entry->p;
		entry->offset = (uint16)(new_offset & (BLOCK_SIZE-1));
		new_bucket = BLOCK_HASH(new_offset);
		if (new_bucket >= num_buckets)
			new_bucket -= num_buckets;

		tmp_offset = entry->next;
		entry->next = bucket_list[new_bucket];
		bucket_list[new_bucket] = fb_offset;
		fb_offset = tmp_offset;
	}
	sieve_fb->num_updates = i;

	/* now do the sieving. First set the initial log value */

	memset(sieve_block, sieve_fb->proj_bias, (size_t)BLOCK_SIZE);

	/* add the updates from the small factor base primes */

	for (i = 0; i < med_fb_size; i++) {
		fb_sieve_entry_t *entry = factor_base + i;
		uint32 p = entry->p;
		uint32 r = entry->offset;
		uint8 logp = entry->logp;

		if (entry->skip)
			continue;

		while (r < BLOCK_SIZE) {
			sieve_block[r] += logp;
			r += p;
		}
		entry->offset = (uint16)(r - BLOCK_SIZE);
	}

	/* add the updates from the powers of small factor base primes */

	for (i = 0; i < powers_fb_size; i++) {
		fb_power_t *entry = sieve_fb->fb_powers + i;
		uint32 p = entry->p;
		uint32 r = entry->offset;
		uint8 logp = entry->logp;
	
		if (entry->skip)
			continue;

		while (r < BLOCK_SIZE) {
			sieve_block[r] += logp;
			r += p;
		}
		entry->offset = (uint16)(r - BLOCK_SIZE);
	}

	/* add the rest of the updates */

	num_updates = sieve_fb->num_updates;
	for (i = 0; i < (num_updates & (uint32)(~7)); i += 8) {

#ifdef MANUAL_PREFETCH
		PREFETCH(update_list + i + 8);
#endif
		sieve_block[update_list[i+0].offset] += update_list[i+0].logp;
		sieve_block[update_list[i+1].offset] += update_list[i+1].logp;
		sieve_block[update_list[i+2].offset] += update_list[i+2].logp;
		sieve_block[update_list[i+3].offset] += update_list[i+3].logp;
		sieve_block[update_list[i+4].offset] += update_list[i+4].logp;
		sieve_block[update_list[i+5].offset] += update_list[i+5].logp;
		sieve_block[update_list[i+6].offset] += update_list[i+6].logp;
		sieve_block[update_list[i+7].offset] += update_list[i+7].logp;
	}
	for (; i < num_updates; i++)
		sieve_block[update_list[i].offset] += update_list[i].logp;
}

/*------------------------------------------------------------------*/
static uint32 do_factoring(sieve_job_t *job, 
			int64 block_start, uint32 b,
			mpz_t log_scratch) {

	/* scan a sieve block for smooth numbers */

	int32 i, j;
	sieve_t *rfb = &job->sieve_rfb;
	sieve_t *afb = &job->sieve_afb;
	uint8 *sieve_r = rfb->sieve_block;
	uint8 *sieve_a = afb->sieve_block;
	int32 cutoff_r, cutoff_a;
	int32 cutoffL_r, cutoffR_r;
	int32 cutoffL_a, cutoffR_a;
	uint32 rbits, abits;
	int64 a, gcda; 
	int32 num_scanned, real_block_length;
	resieve_t *resieve_array = job->resieve_array;
	uint32 num_resieve = 0;
	uint32 resieve_base = 0;
	uint32 rels = 0;

	/* If b is even then the sieve interval is compressed;
	   offset i refers to an uncompressed index of 2*i+1,
	   and the sieve interval represents a region of size
	   2*BLOCK_SIZE so there are twice as many samples */

	if (b % 2 == 0) {
		real_block_length = 2 * BLOCK_SIZE;
		num_scanned = A_SAMPLE_RATE / 2;
	}
	else {
		real_block_length = BLOCK_SIZE;
		num_scanned = A_SAMPLE_RATE;
	}

	/* the log value for the rational sieve is only
	   computed once for this block */

	a = block_start;
	cutoffL_r = fplog_eval_poly(a, b, log_scratch,
				rfb->poly, rfb->log_base, &rbits);
	cutoffR_r = fplog_eval_poly(a + real_block_length, b, log_scratch, 
				rfb->poly, rfb->log_base, &rbits);

	cutoffL_a = fplog_eval_poly(a, b, log_scratch, 
				afb->poly, afb->log_base, &abits);

	for (i = 0; i < BLOCK_SIZE; i += num_scanned) {

		/* calculate the cutoffs for the next
		   num_scanned sieve values */

		a += A_SAMPLE_RATE;
		cutoffR_a = fplog_eval_poly(a, b, log_scratch, afb->poly, 
						afb->log_base, &abits);

		cutoff_r = (cutoffL_r + cutoffR_r) / 2;
		cutoff_a = (cutoffL_a + cutoffR_a) / 2;

		/* decide whether sieve values will have two
     		   large primes or three; this depends only on the size
     		   of the current norms. When the norm on one side
		   is large enough, use three large primes for both 
		   sides. This will noticeably boost the relation
		   discovery rate for only a little more runtime */

		rfb->curr_num_lp = 2;
		afb->curr_num_lp = 2;
		if (abits > 150) {
			afb->curr_num_lp = 3;
			cutoff_a -= afb->scaled_cutoff3;

			if (abits > 185) {
				rfb->curr_num_lp = 3;
				cutoff_r -= rfb->scaled_cutoff3;
			}
			else {
				cutoff_r -= rfb->scaled_cutoff2;
			}
		}
		else if (rbits > 150) {
			rfb->curr_num_lp = 3;
			cutoff_r -= rfb->scaled_cutoff3;

			if (rbits > 185) {
				afb->curr_num_lp = 3;
				cutoff_a -= afb->scaled_cutoff3;
			}
			else {
				cutoff_a -= afb->scaled_cutoff2;
			}
		}
		else {
			cutoff_r -= rfb->scaled_cutoff2;
			cutoff_a -= afb->scaled_cutoff2;
		}

		if (cutoff_r < 0)
			cutoff_r = 0;
		if (cutoff_a < 0)
			cutoff_a = 0;
		cutoffL_a = cutoffR_a;

		for (j = 0; j < num_scanned; j++) {
			int64 curr_a;
			resieve_t *r;

			/* we gradually increase the amount of effort
			   expended to find smooth relations. First check
			   that the logs of both sieve values meet the cutoff.
			   Check the algebraic value first, since it's less
			   likely to pass (this keeps the rational sieve
			   from polluting the cache most of the time) */

			if (sieve_a[i+j] <= cutoff_a ||
			    sieve_r[i+j] <= cutoff_r) {
				sieve_a[i+j] = RESIEVE_INVALID;
				continue;
			}

			if (b % 2 == 0)
				curr_a = block_start + 2 * (i+j) + 1;
			else
				curr_a = block_start + (i+j);

			/* do a little more work */

			gcda = curr_a % (int64)b;
			if (gcda < 0)
				gcda += b;

			if (mp_gcd_1((uint32)gcda, b) != 1) {
				sieve_a[i+j] = RESIEVE_INVALID;
				continue;
			}

			/* queue this value for resieving. There are a lot
			   of compromises in the resieving implementation;
			   the biggest are that we will not buffer too many
			   values at any given time, and only buffer one
			   sieve block worth of values. There are many
			   advantages to this approach: no bitfield is 
			   necessary to represent queued sieve values, and
			   no hashtable is needed to identify them. The
			   byte array for the sieve block can do both of
			   these functions, and initializing it is cheap.
			   A small resieve array also fits better in cache.

			   The only disadvantage is that there may not be
			   enough buffer slots for all the sieve values that
			   need resieving. Rather than ignoring potentially
			   good values in this case, we can just resieve a
			   given block more than once. In cases where it's
			   needed, the savings from resieving everything
			   dwarf the extra overhead of sieving more than once */

			r = resieve_array + 2 * num_resieve;
			r[0].offset = i+j;
			r[0].num_factors = 0;
			r[1].num_factors = 0;
			sieve_a[i+j] = num_resieve++;

			if (num_resieve == MAX_RESIEVE_ENTRIES) {
				rels += do_resieve(job, resieve_base, 
						(uint32)(i+j), num_resieve,
						block_start, b);
				num_resieve = 0;
				resieve_base = i+j+1;
			}
		}
	}

	if (num_resieve > 0) {
		rels += do_resieve(job, resieve_base, BLOCK_SIZE - 1,
					num_resieve, block_start, b);
	}
	return rels;
}

/*------------------------------------------------------------------*/
static uint32 do_resieve(sieve_job_t *job, uint32 low_offset,
			uint32 high_offset, uint32 num_resieve, 
			int64 block_start, uint32 b) {

	/* perform resieving for all the values between 
	   low_offset and high_offset (inclusive) in a 
	   sieve block */

	uint32 i;
	uint32 rels = 0;
	uint8 *hashtable = job->sieve_afb.sieve_block;
	resieve_t *resieve_array = job->resieve_array;

	/* if too few values are queued for resieving, just
	   attempt to trial factor them directly. The ability
	   to selectively turn resieving on and off in a sieve
	   block is very important: our blocks are small, and
	   even small values of 'b' have long stretches that
	   contain few worthwhile sieve values. We don't want 
	   the overhead of resieving unless the benefit is
	   overwhelming */

	if (num_resieve < 10) {
		for (i = 0; i < num_resieve; i++) {
			resieve_t *r = resieve_array + 2 * i;
			r[0].num_factors = RESIEVE_INVALID;
			r[1].num_factors = RESIEVE_INVALID;
			rels += do_one_factoring(job, r, block_start, b);
		}
		return rels;
	}

	/* do the resieving, first on the algebraic and then
	   on the rational parts of each queued value. Even-
	   indexed entries in resieve_array handle algebraic
	   factors, odd-indexed values handle the rational side */

	do_one_resieve(&job->sieve_afb, resieve_array,
			hashtable, low_offset, high_offset);

	do_one_resieve(&job->sieve_rfb, resieve_array + 1,
			hashtable, low_offset, high_offset);

	/* complete the trial factoring of each value 
	   individually */

	for (i = 0; i < num_resieve; i++)
		rels += do_one_factoring(job, resieve_array + 2 * i, 
					block_start, b);
	return rels;
}

/*------------------------------------------------------------------*/
static void do_one_resieve(sieve_t *fb, resieve_t *resieve_array, 
			uint8 *hashtable, uint32 low_offset, 
			uint32 high_offset) {

	/* handle resieving for one factor base */

	uint32 i;
	uint32 fb_start = fb->small_fb_size;
	uint32 fb_end = fb->med_fb_size;
	packed_fb_t *update_list = fb->update_list;
	uint32 num_updates = fb->num_updates;
	fb_sieve_entry_t *fb_entries = fb->fb_entries;

	/* for factor base primes smaller than the block size,
	   the current offset is the first into the *next*
	   sieve block. Hence we run the sieve *backwards*
	   for each such prime */

	for (i = fb_start; i < fb_end; i++) {
		uint32 p = fb_entries[i].p;
		uint32 r = fb_entries[i].offset + BLOCK_SIZE;

		if (fb_entries[i].skip)
			continue;

		/* back up (without touching memory) until 
		   sieve values below the upper limit of the 
		   region of interest are encountered */

		while ((int32)r > (int32)high_offset)
			r -= p;

		/* continue working backwards, and now add p 
		   to the list of factors for any sieve value 
		   that is encountered */

		while ((int32)r >= (int32)low_offset) {
			if (hashtable[r] != RESIEVE_INVALID) {
				resieve_t *entry = resieve_array + 
						2 * hashtable[r];

				/* if there are too many factors in the list,
				   turn resieving off for this sieve value
				   and this factor base */

				if (entry->num_factors >= NUM_RESIEVE_FACTORS)
					entry->num_factors = RESIEVE_INVALID;
				else
					entry->factors[
						entry->num_factors++] = p;
			}
			r -= p;
		}
	}

	/* repeat the process with the large factor base
	   primes. Each such prime only hits in the sieve
	   interval once, so no sieving is needed per se */

	for (i = 0; i < num_updates; i++) {
		uint32 p = update_list[i].p;
		uint32 r = update_list[i].offset;
		if (r >= low_offset && r <= high_offset) {
			if (hashtable[r] != RESIEVE_INVALID) {
				resieve_t *entry = resieve_array + 
						2 * hashtable[r];
				if (entry->num_factors >= NUM_RESIEVE_FACTORS)
					entry->num_factors = RESIEVE_INVALID;
				else
					entry->factors[
						entry->num_factors++] = p;
			}
		}
	}
}

/*------------------------------------------------------------------*/
static uint32 do_one_factoring(sieve_job_t *job, resieve_t *sieve_value,
				int64 block_start, uint32 b) {

	/* Determine the factorization of one sieve value */

	sieve_t *rfb = &job->sieve_rfb;
	sieve_t *afb = &job->sieve_afb;
	uint32 offset = sieve_value->offset;
	int64 a;
	uint32 num_factors_r;
	uint32 factors_r[100];
	uint32 num_factors_a;
	uint32 factors_a[100];
	sieve_t *small, *large;

	if (b % 2 == 0)
		a = block_start + 2 * offset + 1;
	else
		a = block_start + offset;

	/* do trial factoring for each sieve separately.
	   Start with the algebraic sieve value, since the
	   vast majority of such trial factoring attempts
	   will not succeed */

	eval_poly(afb->res, a, b, afb->poly);
	mpz_abs(afb->res, afb->res);

	if (!do_one_tf(afb, sieve_value, 
			factors_a, &num_factors_a, offset, b))
		return 0;

	eval_poly(rfb->res, a, b, rfb->poly);
	mpz_abs(rfb->res, rfb->res);

	if (!do_one_tf(rfb, sieve_value + 1, 
			factors_r, &num_factors_r, offset, b))
		return 0;

	/* continue only if the cofactors remaining after 
	   trial division are composite. A base-2 pseudoprime 
	   test is used to determine primality, since 
	   mp_is_prime does unnecessary trial division 
	   
	   With two large primes per side, the time needed for
	   primality tests is not significant. However, with 
	   three large primes the number of sieve reports is
	   very large, and resieving makes trial factoring of
	   sieve reports very cheap, so the exponentiations
	   below take a significant chunk of the total time.
	   In that case it's a big win to save the primality test
	   until we're sure it's needed, and to do the smallest
	   of the two tests first */
		 
	small = rfb;
	large = afb;
	if (mpz_cmp(small->res, large->res) > 0) {
		small = afb;
		large = rfb;
	}
	if (mpz_cmp_ui(small->res, (uint32)(-1)) > 0) {
		mpz_set_ui(small->tmp2, 2);
		mpz_sub_ui(small->tmp3, small->res, 1);
		mpz_powm(small->tmp2, small->tmp2, small->tmp3, small->res);
		if (mpz_cmp_ui(small->tmp2, 1) == 0)
			return 0;
	}
	if (mpz_cmp_ui(large->res, (uint32)(-1)) > 0) {
		mpz_set_ui(large->tmp2, 2);
		mpz_sub_ui(large->tmp3, large->res, 1);
		mpz_powm(large->tmp2, large->tmp2, large->tmp3, large->res);
		if (mpz_cmp_ui(large->tmp2, 1) == 0)
			return 0;
	}

	if (mpz_cmp_ui(rfb->res, (uint32)(-1)) < 0 &&
	    mpz_cmp_ui(afb->res, (uint32)(-1)) < 0) {

		/* No extra cofactors in this relation, so just
		   save it immediately */

		uint32 i;
		uint32 lp_r[MAX_LARGE_PRIMES];
		uint32 lp_a[MAX_LARGE_PRIMES];
		
		lp_r[0] = mpz_get_ui(rfb->res);
		lp_a[0] = mpz_get_ui(afb->res);
		for (i = 1; i < MAX_LARGE_PRIMES; i++)
			lp_r[i] = lp_a[i] = 1;

		print_relation(&job->obj->savefile, a, b, 
				factors_r, num_factors_r, lp_r,
				factors_a, num_factors_a, lp_a);
		return 1;
	}

	/* schedule the relation to be factored later */

	relation_batch_add(a, b, factors_r, num_factors_r, rfb->res,
				factors_a, num_factors_a, afb->res,
				&job->relation_batch);

	/* if enough unfactored relations have accumulated,
	   factor them */

	if (job->relation_batch.num_relations >= 
			job->relation_batch.target_relations) {
		return relation_batch_run(&job->relation_batch);
	}
	return 0;
}

/*------------------------------------------------------------------*/
static void divide_out_p(sieve_t *fb, uint32 p) {

	uint32 rem;

	while (1) {
		rem = mpz_tdiv_q_ui(fb->tmp2, fb->res, p);
		if (rem != 0)
			return;
		mpz_swap(fb->res, fb->tmp2);
	}
}

static uint32 do_one_tf(sieve_t *sieve_fb, resieve_t *resieve,
			uint32 *factors, uint32 *num_factors_out,
			uint32 offset, uint32 b) {

	/* Trial factor a sieve value (in sieve_fb->res) using 
	   the primes in one factor base */

	uint32 num_factors = 0;
	uint32 i;
	uint32 small_fb_size = sieve_fb->small_fb_size;
	uint32 med_fb_size = sieve_fb->med_fb_size;
	fb_sieve_entry_t *factor_base = sieve_fb->fb_entries;
	uint32 proj_fb_size = sieve_fb->proj_fb_size;
	projective_t *proj_factor_base = sieve_fb->proj_entries;
	packed_fb_t *update_list = sieve_fb->update_list;
	uint32 num_updates = sieve_fb->num_updates;
	uint32 num_resieve_factors = resieve->num_factors;

	/* Use trial division for the small factor base
	   primes. Cofactors that are divisible by a prime
	   p must lie on a particular arithmetic progression.
	   We have another point (offset+BLOCK_SIZE) that
	   also lies on that progression, so if the sieve
	   offset is a multiple of p away from this point the
	   sieve value is divisible by p. 

	   Note that the following will also handle powers
	   of factor base primes automatically, since the root
	   of a prime power lies on the same arithmetic 
	   progression as that of its underlying prime */

	for (i = 0; i < small_fb_size; i++) {
		fb_sieve_entry_t *entry = factor_base + i;
		uint32 p = entry->p;
		uint32 r = entry->r;

		if (entry->skip)
			continue;

		r = (uint32)entry->offset + BLOCK_SIZE - offset;
		if (r % p == 0) {
			if (p >= MAX_SKIPPED_FACTOR)
				factors[num_factors++] = p;
			divide_out_p(sieve_fb, p);
		}
	}

	if (num_resieve_factors != RESIEVE_INVALID) {

		/* the factors from the rest of the factor
		   base are available already; just determine 
		   their multiplicity */

		for (i = 0; i < num_resieve_factors; i++) {
			uint32 p = resieve->factors[i];
			factors[num_factors++] = p;
			divide_out_p(sieve_fb, p);
		}
	}
	else {
		/* if resieving failed or was not performed,
		   continue trial factoring */

		for (; i < med_fb_size; i++) {
			fb_sieve_entry_t *entry = factor_base + i;
			uint32 p = entry->p;
			uint32 r = entry->r;
	
			if (entry->skip)
				continue;
	
			r = (uint32)entry->offset + BLOCK_SIZE - offset;
			if (r % p == 0) {
				factors[num_factors++] = p;
				divide_out_p(sieve_fb, p);
			}
		}
	
		/* For the remainder of the factor base, reuse
		   the list of sieve updates to determine which of
		   them hit the current offset. Not only is this
		   division-free, but the size of the list is dozens
		   to hundreds of times smaller than the complete
		   factor base */
	
		for (i = 0; i < num_updates; i++) {
			uint32 list_offset = update_list[i].offset;
#ifdef MANUAL_PREFETCH
			if (i % 16 == 0)
				PREFETCH(update_list + i + 16);
#endif
	
			if (offset == list_offset) {
				uint32 p = update_list[i].p;
				factors[num_factors++] = p;
				divide_out_p(sieve_fb, p);
			}
		}
	}

	/* Now test the projective roots. All such 
	   primes p divide all sieve values if p 
	   divides b. Note that if b divides p then
	   any other factor base entries for p have
	   been turned off, so it's safe to start
	   pulling out factors of p here. We could 
	   move this up into the resieving phase to 
	   make it much cheaper, but it's probably 
	   not worth the effort */

	for (i = 0; i < proj_fb_size; i++) {
		uint32 p = proj_factor_base[i].p;

		if (b % p == 0) {
			if (p >= MAX_SKIPPED_FACTOR)
				factors[num_factors++] = p;
			divide_out_p(sieve_fb, p);
		}
	}

	/* continue only if the cofactor remaining after 
	   trial division is small, or the right size 
	   for large primes */
		 
	if (sieve_fb->curr_num_lp == 3) {
		if (mpz_cmp(sieve_fb->res, sieve_fb->LP3_max) > 0)
			return 0;

		if (mpz_cmp(sieve_fb->res, sieve_fb->LP3_min) < 0 &&
		    mpz_cmp(sieve_fb->res, sieve_fb->LP2_max) > 0) {
			return 0;
		}
	}
	else {
		if (mpz_cmp(sieve_fb->res, sieve_fb->LP2_max) > 0)
			return 0;
	}

	if (mpz_cmp_ui(sieve_fb->res, sieve_fb->LP1_max) > 0 &&
	    mpz_cmp(sieve_fb->res, sieve_fb->LP2_min) < 0) {
		return 0;
	}

	*num_factors_out = num_factors;
	return 1;
}
