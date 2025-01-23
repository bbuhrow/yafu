/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: sieve.c 820 2012-11-17 03:26:17Z jasonp_sf $
--------------------------------------------------------------------*/

#include <ms_common.h>
#include "mpqs.h"

static void collect_relations(sieve_conf_t *conf,
			      uint32 target_relations,
			      qs_core_sieve_fcn core_sieve_fcn);

static uint32 do_sieving_internal(sieve_conf_t *conf,
				  uint32 max_relations,
				  qs_core_sieve_fcn core_sieve_fcn);

#ifdef SIEVE_TIMING
#define PRINT_TIME(var) printf(#var ": %lf (%4.1f%%)\n",		\
		(double)(var) / 1.86e9,					\
		100.0 * (double)(var) / (double)total_time);

uint64 total_time;
uint64 base_poly_time;
uint64 next_poly_small_time;
uint64 next_poly_large_time;
uint64 plus_init_time;
uint64 bucket_time;
uint64 sieve_small_time;
uint64 sieve_large_time;
uint64 tf_plus_scan_time;
uint64 tf_total_time;
uint64 tf_small_time;
uint64 tf_large_time;

#else
#define PRINT_TIME(var) /* nothing */
#endif

/*--------------------------------------------------------------------*/
void do_sieving(msieve_obj *obj, mp_t *n, 
		mp_t **poly_a_list, poly_t **poly_list,
		fb_t *factor_base, uint32 *modsqrt_array,
		sieve_param_t *params, uint32 multiplier,
		relation_t **relation_list, uint32 *num_relations,
		la_col_t **cycle_list, uint32 *num_cycles) {

	sieve_conf_t conf;
	uint32 bound;
	uint32 i;
	uint32 bits;
	uint32 fb_size = params->fb_size;
	uint32 sieve_size = params->sieve_size;
	uint32 num_sieve_blocks;
	uint32 error_bits;
	uint32 max_relations, relations_found;
	uint32 sieve_block_size;
	uint32 recip_cutoff;
	qs_core_sieve_fcn core_sieve_fcn;

	/* fill in initial sieve parameters */

	memset(&conf, 0, sizeof(conf));
	conf.obj = obj;
	conf.n = n;
	conf.multiplier = multiplier;
	conf.factor_base = factor_base;
	conf.modsqrt_array = modsqrt_array;
	conf.fb_size = params->fb_size;
	bits = mp_bits(conf.n);

	/* decide on the size of one sieve block. If the L1
	   cache size is 32kB then this is also the sieve block
	   size, otherwise it is 64kB. The latter rule is needed
	   because Pentium 4 CPUs have very small L1 sizes, and
	   are designed to work out of L2 cache most of the
	   time anyway */

	if (obj->cache_size1 == 32768)
		conf.sieve_block_size = sieve_block_size = 32768;
	else 
		conf.sieve_block_size = sieve_block_size = 65536;
	conf.sieve_array = (uint8 *)aligned_malloc(
					(size_t)sieve_block_size, 64);

	/* decide on the core sieving routine to use */

	core_sieve_fcn = NULL;
#if defined(HAS_MSVC_SIEVE_CORE)
	if (sieve_block_size == 32768) {
		logprintf(obj, "using VC8 32kb sieve core\n");
		core_sieve_fcn = qs_core_sieve_vc8_32k;
	}
	else {
		logprintf(obj, "using VC8 64kb sieve core\n");
		core_sieve_fcn = qs_core_sieve_vc8_64k;
	}
#else
	if (sieve_block_size == 32768) {
		logprintf(obj, "using generic 32kb sieve core\n");
		core_sieve_fcn = qs_core_sieve_generic_32k;
	}
	else {
		logprintf(obj, "using generic 64kb sieve core\n");
		core_sieve_fcn = qs_core_sieve_generic_64k;
	}
#endif

	/* round the size of the sieve interval (positive plus
	   negative parts combined) up to an integral number 
	   of sieve blocks */

	num_sieve_blocks = (2 * sieve_size + sieve_block_size - 1) /
					sieve_block_size;

	logprintf(obj, "sieve interval: %u blocks of size %u\n", 
			num_sieve_blocks, sieve_block_size);

	/* choose the largest factor base prime that is not
	   sieved, and that uses small reciprocals during 
	   trial factoring. The small prime variation becomes
	   less important as the input size increases, and
	   we reduce its influence in that case */

	if (bits > 350)
		bound = 15;
	else if (bits > 290)
		bound = 20;
	else if (bits > 230)
		bound = 22;
	else if (bits > 200)
		bound = 22;
	else
		bound = 50;
	bound = MIN(bound, fb_size);

	/* to avoid using explicit remainder operations
	   when trial factoring by factor base primes,
	   we multiply by precomputed reciprocals of those
	   primes. See Agner Fog's brilliant optimization
	   manuals for why this works (www.agner.org) 
	 
	   The reciprocals must all fit in 32 bits, and be
	   given inputs not exceeding 2^r for r a compile-time
	   constant. Each reciprocal value for prime p is either 
	   2^r / p or 2^r / p + 1, so there is a range of factor 
	   base entries that have r = 32 (for small p) and
	   another range that has r = 40 (for larger p). 
	 
	   range 1: p < 256 or small prime bound, r = 32 */

	recip_cutoff = ((uint64)1 << 32) / ( 
				(uint64)num_sieve_blocks * sieve_block_size);	
	if (recip_cutoff < 256) {
		logprintf(obj, "error1: factor base and/or sieve interval "
				"is too large\n");
		exit(-1);
	}
	for (i = MIN_FB_OFFSET + 1; i < bound; i++) {
		fb_t *fbptr = conf.factor_base + i;
		uint32 prime = fbptr->prime;
		if (prime >= 1000 || prime >= recip_cutoff)
			break;

		fbptr->recip = ((uint64)1 << 32) / (uint64)prime;
		if (floor(MP_RADIX / (double)prime + 0.5) ==
						(double)fbptr->recip) {
			fbptr->rcorrect = 1;
		}
		else {
			fbptr->rcorrect = 0;
			fbptr->recip++;
		}
	}
	conf.tf_small_recip1_cutoff = i;

	/* range 2: p < (cutoff for unsieved factor base entries), r = 40.
	   This range is probably empty */

	for (; i < bound; i++) {
		fb_t *fbptr = conf.factor_base + i;
		uint32 prime = fbptr->prime;
		if (prime >= 1000)
			break;

		fbptr->recip = ((uint64)1 << 40) / (uint64)prime;
		if (floor(256 * MP_RADIX / (double)prime + 0.5) ==
						(double)fbptr->recip) {
			fbptr->rcorrect = 1;
		}
		else {
			fbptr->rcorrect = 0;
			fbptr->recip++;
		}
	}
	conf.tf_small_recip2_cutoff = i;
	conf.sieve_small_fb_start = i;

	/* range 3: p < block size, r = 32 */

	for (; i < fb_size; i++) {
		fb_t *fbptr = conf.factor_base + i;
		uint32 prime = fbptr->prime;
		if (prime >= sieve_block_size || prime >= recip_cutoff)
			break;

		fbptr->recip = ((uint64)1 << 32) / (uint64)prime;
		if (floor(MP_RADIX / (double)prime + 0.5) ==
						(double)fbptr->recip) {
			fbptr->rcorrect = 1;
		}
		else {
			fbptr->rcorrect = 0;
			fbptr->recip++;
		}
	}
	conf.tf_med_recip1_cutoff = i;

	/* repeat, but for large reciprocals. An implicit
	   limit here is that the trial division code will not
	   work if the product of the largest such factor base
	   prime and the largest sieve offset exceeds 2^40 

	   range 4: p < block size, r = 40 */

	recip_cutoff = ((uint64)1 << 40) / (
				(uint64)num_sieve_blocks * sieve_block_size);	
	if (recip_cutoff < sieve_block_size) {
		logprintf(obj, "error: factor base and/or sieve interval "
				"is too large\n");
		exit(-1);
	}
	for (; i < fb_size; i++) {
		fb_t *fbptr = conf.factor_base + i;
		uint32 prime = fbptr->prime;
		if (prime >= sieve_block_size)
			break;

		fbptr->recip = ((uint64)1 << 40) / (uint64)prime;
		if (floor(256 * MP_RADIX / (double)prime + 0.5) ==
						(double)fbptr->recip) {
			fbptr->rcorrect = 1;
		}
		else {
			fbptr->rcorrect = 0;
			fbptr->recip++;
		}
	}
	conf.tf_med_recip2_cutoff = i;

	/* choose the largest factor base prime that does not 
	   use hashtables or reciprocals. We could make this
	   range empty, but it's better to start using hashtables
	   with primes somewhat larger than the block size. 
	   This greatly reduces the amount of memory used by 
	   the hashtables, and trial factoring by primes in this 
	   range is extremely cheap */

	for (; i < fb_size && conf.factor_base[i].prime < 
				3 * sieve_block_size; i++) {
		/* nothing */
	}
	conf.tf_large_cutoff = i;
	conf.sieve_large_fb_start = i;
	conf.packed_fb = (packed_fb_t *)xmalloc(i * sizeof(packed_fb_t));

	/* The sieve code is optimized for sieving intervals that are
	   extremely small. To reduce the overhead of using a large
	   factor base, cache blocking is used for the sieve interval
	   *and* the factor base; this requires that sieving takes place
	   for several polynomials simultaneously. */
	
	/* the number of polynomials that are simultaneously
	   sieved depends on how long each sieve interval is;
	   more polynomials are more efficient but consume
	   more memory */

	conf.poly_block = (uint32)((double)65536 / sieve_block_size *
			100 / num_sieve_blocks + 1);

	logprintf(obj, "processing polynomials in batches of %u\n", 
			conf.poly_block);

	conf.fb_block = 200;

 	bound = conf.factor_base[conf.fb_size - 1].prime;
	logprintf(obj, "using a sieve bound of %u (%u primes)\n", 
				bound, conf.fb_size);

	/* The trial factoring code can only handle 32-bit
	   factors, so the single large prime bound must
	   fit in 32 bits. Note that the SQUFOF code can
	   only manage 31-bit factors, so that single large
	   prime relations can have a 32-bit factor but
	   double large prime relations are limited to 31-bit
	   factors */

	if ((uint32)(-1) / bound <= params->large_mult)
		bound = (uint32)(-1);
	else
		bound *= params->large_mult;

	logprintf(obj, "using large prime bound of %u (%d bits)\n", 
			bound, (int32)(log((double)bound) / M_LN2));

	/* compute max_fb2, the square of the largest factor base prime */

	i = conf.factor_base[fb_size - 1].prime;
	mp_clear(&conf.max_fb2);
	conf.max_fb2.nwords = 1;
	conf.max_fb2.val[0] = i;
	mp_mul_1(&conf.max_fb2, i, &conf.max_fb2);

	mp_clear(&conf.large_prime_max2);

	/* fill in sieving cutoffs that leverage the small
	   prime variation. Do not skip small primes if the
	   factor base is too small, and skip fewer small
	   primes when the factor base is large (i.e. when it's
	   more important to increase the relation rate) */

	error_bits = (uint32)(log((double)bound) / M_LN2 + 0.5);

	if (fb_size < 800) {
		conf.cutoff2 = (uint32)(1.5 * error_bits);
		conf.cutoff1 = conf.cutoff2;
		conf.sieve_small_fb_start = MIN_FB_OFFSET + 1;
		/* do not change the trial factoring ranges, since
		   they depend on the choice of reciprocal used */
	}
	else {
		/* Turn on double large primes if n is 85 digits or more */

		if (bits >= 282) {
			mp_t *large_max2 = &conf.large_prime_max2;

			/* the double-large-prime cutoff is equal to
			   the single-large-prime cutoff raised to the 1.8
			   power. We don't square the single cutoff because
			   relations with two really big large primes are
			   very rare, and will almost never survive the 
			   filtering step anyway */

			mp_clear(large_max2);
			large_max2->nwords = 1;
			large_max2->val[0] = bound;
			mp_mul_1(large_max2, 
				(uint32)pow((double)bound, 0.8), 
				large_max2);

			conf.cutoff2 = mp_bits(large_max2);
			logprintf(obj, "using double large prime bound of "
				"%s (%u-%u bits)\n",
				mp_sprintf(large_max2, 10, obj->mp_sprintf_buf),
				mp_bits(&conf.max_fb2), 
				conf.cutoff2);
		}
		else {
			conf.cutoff2 = error_bits;
		}
		conf.cutoff1 = conf.cutoff2 + 30;
		logprintf(obj, "using trial factoring cutoff "
				"of %u bits\n", conf.cutoff2);
	}

	/* set up the hashtable if it's going to be used; otherwise
	   the rest of the sieve code will silently ignore it.
	   
	   There is one hash bin for every sieve block in the sieve
	   interval; there are also twice as many such bins because
	   we collect positive and negative sieve offsets into separate
	   hashtables. Finally, we sieve over up to POLY_BLOCK_SIZE
	   polynomials at the same time, so the number of hash bins
	   is similarly multiplied.
	   
	   Needless to say, that's a *lot* of hash bins! The sieve
	   interval must be extremely small to keep the hashtable
	   size down; that's okay, very small sieve intervals are
	   actually more likely to contain smooth relations */

	conf.buckets = (bucket_t *)xcalloc((size_t)(conf.poly_block *
							num_sieve_blocks), 
							sizeof(bucket_t));
	if (fb_size > conf.sieve_large_fb_start) {
		for (i = 0; i < conf.poly_block * num_sieve_blocks; i++) {
			conf.buckets[i].num_alloc = 1000;
			conf.buckets[i].list = (bucket_entry_t *)
					xmalloc(1000 * sizeof(bucket_entry_t));
		}
	}

	/* fill in miscellaneous parameters */

	conf.num_sieve_blocks = num_sieve_blocks;
	conf.large_prime_max = bound;

	/* initialize the polynomial generation code. Note that
	   we do *not* use the sieve size that is input to specify
	   the sizing of polynomials, but instead use the rounded value 
	   derived from it */

	poly_init(&conf, num_sieve_blocks * sieve_block_size / 2);

	/* initialize the bookkeeping for tracking partial relations */

	conf.components = 0;
	conf.vertices = 0;
	if (!(obj->flags & MSIEVE_FLAG_SKIP_QS_CYCLES)) {
		conf.cycle_hashtable = (uint32 *)xcalloc(
					(size_t)(1 << LOG2_CYCLE_HASH),
					sizeof(uint32));
		conf.cycle_table_size = 1;
		conf.cycle_table_alloc = 10000;
		conf.cycle_table = (cycle_t *)xmalloc(conf.cycle_table_alloc * 
							sizeof(cycle_t));
	}

	/* stop when this many relations are found */

	max_relations = obj->max_relations;
	if (obj->max_relations == 0)
		max_relations = fb_size + 3 * NUM_EXTRA_RELATIONS / 2;

	/* do all the sieving, and save the relations 
	   that were generated */

	obj->flags |= MSIEVE_FLAG_SIEVING_IN_PROGRESS;

	TIME1(total_time)
	relations_found = do_sieving_internal(&conf, max_relations,
						core_sieve_fcn);
	TIME2(total_time)

	PRINT_TIME(total_time);
	PRINT_TIME(base_poly_time);
	PRINT_TIME(next_poly_small_time);
	PRINT_TIME(next_poly_large_time);
	PRINT_TIME(plus_init_time);
	PRINT_TIME(bucket_time);
	PRINT_TIME(sieve_small_time);
	PRINT_TIME(sieve_large_time);
	PRINT_TIME(tf_plus_scan_time);
	PRINT_TIME(tf_total_time);
	PRINT_TIME(tf_small_time);
	PRINT_TIME(tf_large_time);

	/* free all of the sieving structures first, to leave
	   more memory for the postprocessing step. Do *not* free
	   the polynomial subsystem yet. Sieving is only considered
	   finished when the savefile has received all pending output */

	savefile_flush(&obj->savefile);
	savefile_close(&obj->savefile);
	obj->flags &= ~MSIEVE_FLAG_SIEVING_IN_PROGRESS;

	for (i = 0; i < conf.poly_block * num_sieve_blocks; i++) {
		free(conf.buckets[i].list);
	}
	free(conf.buckets);
	free(conf.packed_fb);
	aligned_free(conf.sieve_array);

	/* if enough relations are available, do the postprocessing
	   and save the results where the rest of the program can
	   find them. Don't run filter_relations() if the cycle-
	   counting structures have not been initialized */

	if (relations_found >= max_relations &&
	    !(obj->flags & MSIEVE_FLAG_SKIP_QS_CYCLES)) {
	        if(obj->flags & (MSIEVE_FLAG_USE_LOGFILE |
	    		   MSIEVE_FLAG_LOG_TO_STDOUT)) {
			fprintf(stderr, "sieving complete, "
					"commencing postprocessing\n");
		}
		qs_filter_relations(&conf);
		*relation_list = conf.relation_list;
		*num_relations = conf.num_relations;
		*cycle_list = conf.cycle_list;
		*num_cycles = conf.num_cycles;
		*poly_a_list = conf.poly_a_list;
		*poly_list = conf.poly_list;
	}

	poly_free(&conf);
	free(conf.cycle_table);
	free(conf.cycle_hashtable);
}

/*--------------------------------------------------------------------*/
static uint32 do_sieving_internal(sieve_conf_t *conf, 
				uint32 max_relations,
				qs_core_sieve_fcn core_sieve_fcn) {

	uint32 num_relations = 0;
	uint32 update;
	msieve_obj *obj = conf->obj;
	savefile_t *savefile = &obj->savefile;
	char buf[300];

	/* open the savefile; if the file already
	   exists and the first line contains n in base 16,
	   then we are restarting from a previous factorization */

	update = 1;
	if (savefile_exists(savefile)) {
		savefile_open(savefile, SAVEFILE_READ);
		buf[0] = 0;
		savefile_read_line(buf, sizeof(buf), savefile);
		if (isxdigit(buf[0])) {
			mp_t read_n;
			mp_str2mp(buf, &read_n, 16);
			if (mp_cmp(conf->n, &read_n) == 0)
				update = 0;
		}
		savefile_close(savefile);
	}

	if (update) {
		/* Truncate the file and write the present n. 
		   I hope you backed up savefiles you wanted! */

		savefile_open(savefile, SAVEFILE_WRITE);
		savefile_write_line(savefile, mp_sprintf(
					conf->n, 16, obj->mp_sprintf_buf));
		savefile_write_line(savefile, "\n");
		savefile_flush(savefile);
		savefile_close(savefile);
	}

	/* Read in the large primes for all the relations
	   in the savefile, and count the number of
	   cycles that can be formed. Note that no check
	   for duplicate or corrupted relations is made here;
	   the cycle finder will rebuild everything from 
	   scratch when sieving finishes, and does all the 
	   verification at that point */

	savefile_open(savefile, SAVEFILE_READ);
	savefile_read_line(buf, sizeof(buf), savefile);

	while (!savefile_eof(savefile)) {
		uint32 prime1, prime2;
		char *tmp = strchr(buf, 'L');

		if (buf[0] != 'R' || tmp == NULL) {
			savefile_read_line(buf, sizeof(buf), savefile);
			continue;
		}

		read_large_primes(tmp, &prime1, &prime2);
		if (prime1 == prime2) {
			conf->num_relations++;
		}
		else {
			add_to_cycles(conf, prime1, prime2);
			conf->num_cycles++;
		}
		savefile_read_line(buf, sizeof(buf), savefile);
	}
	savefile_close(savefile);

	/* prepare the savefile for receiving more relations */

	savefile_open(savefile, SAVEFILE_APPEND);

	if (conf->num_relations || conf->num_cycles) {
		logprintf(obj, "restarting with %u full and %u "
				"partial relations\n", conf->num_relations, 
						conf->num_cycles);
	}

	num_relations = conf->num_relations + 
			conf->num_cycles +
			conf->components - conf->vertices;

	/* choose how many full relations to collect before
	   printing a progress update */

	if (max_relations >= 62000)
		update = 5;
	else if (max_relations >= 50000)
		update = 25;
	else if (max_relations >= 10000)
		update = 100;
	else
		update = 200;
	update = MIN(update, max_relations / 10);

	if (num_relations < max_relations &&
	    (obj->flags & (MSIEVE_FLAG_USE_LOGFILE |
	    		   MSIEVE_FLAG_LOG_TO_STDOUT))) {
		fprintf(stderr, "\nsieving in progress "
				"(press Ctrl-C to pause)\n");
	}

	/* sieve until at least that many relations have
	   been found, then update the number of fulls and
	   partials. This way we can declare sieving to be
	   finished the moment enough relations are available */

	while (!(obj->flags & MSIEVE_FLAG_STOP_SIEVING) && 
		num_relations < max_relations) {

		collect_relations(conf, update, core_sieve_fcn);

		num_relations = conf->num_relations + 
				conf->num_cycles +
				conf->components - conf->vertices;

	    	if (obj->flags & (MSIEVE_FLAG_USE_LOGFILE |
	    		   	  MSIEVE_FLAG_LOG_TO_STDOUT)) {
			fprintf(stderr, "%u relations (%u full + "
				"%u combined from %u partial), need %u\r",
					num_relations,
					conf->num_relations,
					conf->num_cycles +
					conf->components - conf->vertices,
					conf->num_cycles,
					max_relations);
			fflush(stderr);
		}
	}

	if (obj->flags & (MSIEVE_FLAG_USE_LOGFILE |
	    		   MSIEVE_FLAG_LOG_TO_STDOUT))
		fprintf(stderr, "\n");

	logprintf(obj, "%u relations (%u full + %u combined from "
			"%u partial), need %u\n",
				num_relations, conf->num_relations,
				conf->num_cycles +
				conf->components - conf->vertices,
				conf->num_cycles,
				max_relations);
	return num_relations;
}

/*--------------------------------------------------------------------*/
static void collect_relations(sieve_conf_t *conf, 
			      uint32 target_relations,
			      qs_core_sieve_fcn core_sieve_fcn) {
	
	uint32 i;
	uint32 relations_found = 0;
	uint32 num_poly = 1 << (conf->num_poly_factors - 1);
	uint32 poly_block = MIN(num_poly, conf->poly_block);

	/* top-level sieving loop: keep building polynomials
	   and sieving with them until at least target_relations
	   relations have been found */

	while (relations_found < target_relations) {
		
		/* build the next batch of polynomials. For
		   big factorizations there may be thousands
		   of them */

		TIME1(base_poly_time)
		build_base_poly(conf);
		TIME2(base_poly_time)

		/* Do the sieving for all polynomials, handling
		   batches of polynomials at a time. */

		i = 0;
		while (i < num_poly) {
			uint32 curr_num_poly = MIN(poly_block, num_poly - i);

			relations_found += core_sieve_fcn(conf, i, 
							curr_num_poly);
			i += curr_num_poly;

			/* see if anybody wants sieving to stop. If num_poly
	    		   is small we bail out after sieving for the entire 
			   current batch of polynomials has finished. If
			   num_poly is large we wait until a significant
			   number of polynomials have been sieved. Basically
			   you're entitled to a big pile of polynomials
			   whenever the loop runs, and you should get your
			   money's worth without having to wait a really
			   long time to finish up */

			if ((conf->obj->flags & MSIEVE_FLAG_STOP_SIEVING) &&
			    ((num_poly <= 1024 && i == num_poly) ||
			     (num_poly > 1024 && i > 2000))) {
				return;
			}
		}
	}
}

/*--------------------------------------------------------------------*/
static mp_t two = {1, {2}};

uint32 check_sieve_val(sieve_conf_t *conf, int32 sieve_offset,
			uint32 bits, mp_t *a, signed_mp_t *b, signed_mp_t *c,
			uint32 poly_index, bucket_t *hash_bucket) {

	/* check a single sieve value for smoothness. This
	   routine is called quite rarely but is very compu-
	   tationally intensive. Returns 1 if the input sieve
	   value is completely factored, 0 otherwise. 
	   
	   There are a lot of things in this routine that are
	   not explained very well in references. The quadratic 
	   sieve uses a polynomial p(x) = (sqrt(n) + x), and trial 
	   division is performed on p(x)^2 - n. MPQS instead uses
	   a polynomial p(x) = a * x + b, where a and b are several
	   digits smaller than sqrt(n). 
	   
	   Thus, for MPQS p(x)^2 - n is usually much larger than
	   sqrt(n). However, the way a and b were computed, p(x)^2 - n
	   is divisible by a. Rather than compute p(x)^2 - n and divide
	   manually by a, you can instead compute a * x^2 + 2 * b * x + c,
	   where c is precomputed to be (b*b-n)/a. This quadratic 
	   polynomial happens to be (p(x)^2-n)/a, and its value is a
	   little less than sqrt(n) so you can trial divide like in QS.
	   c and 2*b were both precomputed by calling code, although
	   c is not used after the sieving phase and need not be saved */

	uint32 i, j;
	uint32 num_factors = 0;
	uint32 sign_of_offset;
	uint32 fb_offsets[32 * MAX_MP_WORDS / 4];
	mp_t res;
	signed_mp_t polyval;
	fb_t *factor_base = conf->factor_base;
	packed_fb_t *packed_factor_base = conf->packed_fb;
	uint32 tf_small_recip1_cutoff = conf->tf_small_recip1_cutoff;
	uint32 tf_small_recip2_cutoff = conf->tf_small_recip2_cutoff;
	uint32 tf_med_recip1_cutoff = conf->tf_med_recip1_cutoff;
	uint32 tf_med_recip2_cutoff = conf->tf_med_recip2_cutoff;
	uint32 tf_large_cutoff = conf->tf_large_cutoff;
	bucket_entry_t *list;
	uint32 lowbits;
	mp_t exponent, ans;
	uint32 cutoff2;
	uint32 abs_offset;
	uint32 tf_offset;
	uint32 index;
	uint32 sieve_block_size = conf->sieve_block_size;
	int32 sieve_size = conf->num_sieve_blocks * sieve_block_size;

	/* Compute the polynomial index. We work with three numbers:

	     - the sieve array offset, a number in 
               [-sieve_size/2,+sieve_size/2] used to compute the
	       polynomial value

	     - the trial factoring offset, a number in [0, sieve_size] used
	       by the ordinary trial factoring code

	     - the sieve block offset, a number in [0,sieve_block_size] 
	       used by the hashtable-based trial factoring code
	*/

	tf_offset = (sieve_offset + sieve_size / 2);
	index = tf_offset & (sieve_block_size - 1);

	sign_of_offset = POSITIVE;
	abs_offset = (uint32)abs(sieve_offset);
	if (sieve_offset < 0)
		sign_of_offset = NEGATIVE;

	/* compute a*x^2 + 2*b*x + c, or using Horner's rule, 
	   (a*x+2*b)*x + c 

	   Note that the linear algebra code cares about the sign 
	   of the polynomial value, but the MPQS square root phase 
	   cares about the sign of sieve_offset. Thus, you have to 
	   track both these values for every relation */

	mp_mul_1(a, abs_offset, &polyval.num);
	polyval.sign = sign_of_offset;
	signed_mp_add(&polyval, b, &polyval);
	mp_mul_1(&polyval.num, abs_offset, &polyval.num);
	polyval.sign ^= sign_of_offset;
	signed_mp_add(&polyval, c, &polyval);
	mp_copy(&polyval.num, &res);

	if (polyval.sign == NEGATIVE)
		fb_offsets[num_factors++] = 0;

	/* now that the polynomial value has been computed, calculate
	   the exact number of bits that would be needed to indicate
	   trial factoring is worthwhile. If sieve_offset is near
	   a real root of the polynomial then res can be several digits
	   smaller than calling code expects, and we don't want to
	   throw it away by accident (especially since the odds are
	   better that it will be smooth) */

	cutoff2 = mp_bits(&res);
	if (cutoff2 >= conf->cutoff2)
		cutoff2 -= conf->cutoff2;
	else
		cutoff2 = 0;

	/* Do trial division; factor out all primes in
	   the factor base from 'res'. First pull out
	   factors of two */

	lowbits = mp_rjustify(&res, &res);
	bits += lowbits;
	for (i = 0; i < lowbits; i++)
		fb_offsets[num_factors++] = MIN_FB_OFFSET;

	/* Now deal with all the rest of the factor base
	   primes. Rather than perform multiple-precision
	   divides, this code uses the fact that if 'prime'
	   divides 'res', then 'sieve_offset' lies on one 
	   of two arithmetic progressions. In other words,
	   sieve_offset = root1 + i*prime or root2 + j*prime
	   for some i or j, if prime divides res. Because
	   sieve_offset is always less than 32 bits, and
	   root1 and root2 are available, a single remainder
	   operation is enough to determine divisibility.
	   
	   The first step is a small amount of trial division
	   before comparison of log2(sieve_value) to cutoff2.
	   Begin with factor base primes whose reciprocal assumes
	   a numerator up to 2^32 */

	for (i = MIN_FB_OFFSET + 1; i < tf_small_recip1_cutoff; i++) {
		fb_t *fbptr = factor_base + i;
		uint32 prime = fbptr->prime;
		uint32 root1 = fbptr->root1;
		uint32 root2 = fbptr->root2;
		uint32 logprime = fbptr->logprime;
		uint32 recip = fbptr->recip;
		uint32 rcorrect = fbptr->rcorrect;

		/* if the roots have not been computed, do
		   a multiple precision mod. The only value for
		   which this is true (in this loop) should be
		   the small multiplier */

		if (root1 == INVALID_ROOT) {
			if (mp_mod_1(&res, prime) == 0) {
				do {
					bits += logprime;
					fb_offsets[num_factors++] = i;
					mp_divrem_1(&res, prime, &res);
					j = mp_mod_1(&res, prime);
				} while (j == 0);
			}
			continue;
		}

		j = (uint32)(((uint64)(tf_offset + rcorrect) * 
					(uint64)recip) >> 32);
		j = tf_offset - j * prime;
		if (j == root1 || j == root2) {
			do {
				bits += logprime;
				fb_offsets[num_factors++] = i;
				mp_divrem_1(&res, prime, &res);
				j = mp_mod_1(&res, prime);
			} while (j == 0);
		}
	}

	/* continue with primes whose reciprocal assumes 2^40.
	   This range is probably empty */

	for (; i < tf_small_recip2_cutoff; i++) {
		fb_t *fbptr = factor_base + i;
		uint32 prime = fbptr->prime;
		uint32 root1 = fbptr->root1;
		uint32 root2 = fbptr->root2;
		uint32 logprime = fbptr->logprime;
		uint32 recip = fbptr->recip;
		uint32 rcorrect = fbptr->rcorrect;

		if (root1 == INVALID_ROOT) {
			if (mp_mod_1(&res, prime) == 0) {
				do {
					bits += logprime;
					fb_offsets[num_factors++] = i;
					mp_divrem_1(&res, prime, &res);
					j = mp_mod_1(&res, prime);
				} while (j == 0);
			}
			continue;
		}

		j = (uint32)(((uint64)(tf_offset + rcorrect) * 
					(uint64)recip) >> 40);
		j = tf_offset - j * prime;
		if (j == root1 || j == root2) {
			do {
				bits += logprime;
				fb_offsets[num_factors++] = i;
				mp_divrem_1(&res, prime, &res);
				j = mp_mod_1(&res, prime);
			} while (j == 0);
		}
	}

	if (bits <= cutoff2)
		return 0;

	/* Now perform trial division for the rest of the
	   "small" factor base primes. Begin with those whose
	   reciprocal assumes numerators up to 2^32 */

	TIME1(tf_small_time)
	for (; i < tf_med_recip1_cutoff; i++) {
		fb_t *fbptr = factor_base + i;
		uint32 prime = fbptr->prime;
		uint32 root1 = fbptr->root1;
		uint32 root2 = fbptr->root2;
		uint32 recip = fbptr->recip;
		uint32 rcorrect = fbptr->rcorrect;

#ifdef MANUAL_PREFETCH
		if (i % 4 == 0)
			PREFETCH(fbptr + 4);
#endif

		/* if the roots have not been computed, do
		   a multiple precision mod. The only values for
		   which this is true (in this loop) are the
		   factors of the polynomial 'a' value */

		if (root1 == INVALID_ROOT) {
			if (mp_mod_1(&res, prime) == 0) {
				do {
					fb_offsets[num_factors++] = i;
					mp_divrem_1(&res, prime, &res);
					j = mp_mod_1(&res, prime);
				} while (j == 0);
			}
			continue;
		}

		j = (uint32)(((uint64)(tf_offset + rcorrect) * 
					(uint64)recip) >> 32);
		j = tf_offset - j * prime;
		if (j == root1 || j == root2) {
			do {
				fb_offsets[num_factors++] = i;
				mp_divrem_1(&res, prime, &res);
				j = mp_mod_1(&res, prime);
			} while (j == 0);
		}
	}

	/* repeat for numerators up to 2^40 */

	for (; i < tf_med_recip2_cutoff; i++) {
		fb_t *fbptr = factor_base + i;
		uint32 prime = fbptr->prime;
		uint32 root1 = fbptr->root1;
		uint32 root2 = fbptr->root2;
		uint32 recip = fbptr->recip;
		uint32 rcorrect = fbptr->rcorrect;

#ifdef MANUAL_PREFETCH
		if (i % 4 == 0)
			PREFETCH(fbptr + 4);
#endif

		if (root1 == INVALID_ROOT) {
			if (mp_mod_1(&res, prime) == 0) {
				do {
					fb_offsets[num_factors++] = i;
					mp_divrem_1(&res, prime, &res);
					j = mp_mod_1(&res, prime);
				} while (j == 0);
			}
			continue;
		}

		j = (uint32)(((uint64)(tf_offset + rcorrect) * 
					(uint64)recip) >> 40);
		j = tf_offset - j * prime;
		if (j == root1 || j == root2) {
			do {
				fb_offsets[num_factors++] = i;
				mp_divrem_1(&res, prime, &res);
				j = mp_mod_1(&res, prime);
			} while (j == 0);
		}
	}

	/* handle the largest factor base primes that are
	   not hashed. Since by design all such primes p exceed 
	   the size of the sieve block, at most two offsets
	   per sieve block are divisible by p. We have the 
	   offsets of the *next* sieve values that are divisible
	   by p, so if the current offset is p away from either
	   of these then it is divisible by p as well. This test
	   is multiply-and-divide-free, and so is extremely fast */

	for (; i < tf_large_cutoff; i++) {
		packed_fb_t *pfbptr = packed_factor_base + i;
		uint32 prime = pfbptr->prime;
		uint32 root1 = pfbptr->next_loc1;
		uint32 root2 = pfbptr->next_loc2;

#ifdef MANUAL_PREFETCH
		if (i % 4 == 0)
			PREFETCH(pfbptr + 4);
#endif

		if (root1 == INVALID_ROOT) {
			if (mp_mod_1(&res, prime) == 0) {
				do {
					fb_offsets[num_factors++] = i;
					mp_divrem_1(&res, prime, &res);
					j = mp_mod_1(&res, prime);
				} while (j == 0);
			}
			continue;
		}

		j = index + prime - sieve_block_size;
		if (j == root1 || j == root2) {
			do {
				fb_offsets[num_factors++] = i;
				mp_divrem_1(&res, prime, &res);
				j = mp_mod_1(&res, prime);
			} while (j == 0);
		}
	}
	TIME2(tf_small_time)

	list = hash_bucket->list;

	/* This hashtable entry contains all the primes that
	   divide sieve offsets in this range, along with the
	   sieve offsets those primes divide. Hence trial
	   division for the entire factor base above small_fb_size
	   requires only iterating through this hashtable entry
	   and comparing sieve offsets

	   Not only does this not require any divisions, but
	   the number of entries in list[] is much smaller
	   than the full factor base (5-10x smaller) */

	TIME1(tf_large_time)
	for (i = 0; i < hash_bucket->num_used; i++) {

#ifdef MANUAL_PREFETCH
		if (i % 8 == 0)
			PREFETCH(list + i + 16);
#endif
	
		if (list[i].sieve_offset == index) {
			uint32 prime_index = list[i].prime_index;
			uint32 prime = factor_base[prime_index].prime;

			do {
				fb_offsets[num_factors++] = prime_index;
				mp_divrem_1(&res, prime, &res);
				j = mp_mod_1(&res, prime);
			} while (j == 0);
		}
	}
	TIME2(tf_large_time)

	/* encode the sign of sieve_offset into its top bit */

	abs_offset |= sign_of_offset << 31;
	
	/* If 'res' has been completely factored, save a full relation */

	if (res.nwords == 1 && res.val[0] == 1) {
		save_relation(conf, abs_offset, fb_offsets, 
				num_factors, poly_index, 1, 1);
		return 1;
	}

	/* if 'res' is smaller than the bound for partial 
	   relations, save a 1LP-relation */

	if (res.nwords == 1 && res.val[0] < conf->large_prime_max) {
		save_relation(conf, abs_offset, fb_offsets, 
				num_factors, poly_index, 1, res.val[0]);
		return 0;
	}

	/* if 'res' is smaller than the square of the 
	   largest factor base prime, then 'res' is
	   itself prime and is useless for our purposes. 
	   Note that single large prime relations will
	   always fail at this point */
	
	if (mp_cmp(&res, &conf->max_fb2) < 0)
		return 0;
	
	/* 'res' is not too small; see if it's too big */

	if (mp_cmp(&res, &conf->large_prime_max2) > 0)
		return 0;
	
	/* perform a base-2 pseudoprime test to make sure
	   'res' is composite */
	
	mp_sub_1(&res, 1, &exponent);
	mp_expo(&two, &exponent, &res, &ans);
	if (mp_is_one(&ans))
		return 0;
	
	/* *finally* attempt to factor 'res'; if successful,
	   and both factors are smaller than the single 
	   large prime bound, save 'res' as a partial-partial 
	   relation */
	
	i = squfof(&res);
	if (i > 1) {
		mp_divrem_1(&res, i, &res);
		if (i < conf->large_prime_max && res.nwords == 1 &&
				res.val[0] < conf->large_prime_max) {
		    
		        if (i == res.val[0]) {

				/* if 'res' is a perfect square, then this
				   is actually a full relation! */

				save_relation(conf, abs_offset, fb_offsets, 
						num_factors, poly_index,
						i, res.val[0]);
				return 1;
			}
			else {
				save_relation(conf, abs_offset, fb_offsets, 
						num_factors, poly_index,
						i, res.val[0]);
			}
		}
	}
	return 0;
}
