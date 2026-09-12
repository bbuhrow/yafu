/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: stage1.c 1088 2026-05-19 01:23:20Z jasonp_sf $
--------------------------------------------------------------------*/

#include <stage1.h>
#include <stage1_engine.h>

/* main driver for stage 1 */

/*------------------------------------------------------------------------*/
static void
stage1_bounds_update(poly_search_t *poly, poly_coeff_t *c)
{
	/* determine the parameters for the collision search,
	   given one leading algebraic coefficient a_d */

	uint32 degree = poly->degree;
	double N = mpz_get_d(poly->N);
	double high_coeff = mpz_get_d(c->high_coeff);
	double m0 = pow(N / high_coeff, 1./degree);
	double skewness_min, coeff_max;

	/* we don't know the optimal skewness for this polynomial
	   but at least can bound the skewness. The value of the
	   third-highest coefficient is from Kleinjung's 2006
	   poly selection algorithm as published in Math. Comp. */

	switch (degree) {
	case 4:
		skewness_min = sqrt(m0 / poly->norm_max);
		coeff_max = poly->norm_max;
		break;

	case 5:
		skewness_min = pow(m0 / poly->norm_max, 2./3.);
		coeff_max = poly->norm_max / sqrt(skewness_min);
		break;

	case 6:
		skewness_min = sqrt(m0 / poly->norm_max);
		coeff_max = poly->norm_max / skewness_min;
		break;

	default:
		printf("error: unexpected poly degree %d\n", degree);
		exit(-1);
	}

	c->degree = degree;
	c->m0 = m0;
	c->coeff_max = coeff_max;
	c->p_size_max = coeff_max / skewness_min;

	/* we perform the collision search on a transformed version
	   of N and the low-order rational coefficient m. In the
	   transformed coordinates, a_d is 1 and a_{d-1} is 0. When
	   a hit is found, we undo the transformation to recover
	   the correction to m that makes the new polynomial 'work' */

	mpz_mul_ui(c->trans_N, c->high_coeff, degree);
	mpz_pow_ui(c->trans_N, c->trans_N, degree - 1);
	mpz_mul_ui(c->trans_N, c->trans_N, degree);
	mpz_mul(c->trans_N, c->trans_N, poly->N);
	mpz_root(c->trans_m0, c->trans_N, degree);
}

/*------------------------------------------------------------------------*/
/* infrastructure for submitting stage 1 hits to the stage 2 thread pool */

typedef struct {
	stage1_callback_t callback;
	void *callback_data;
	poly_stage_stats_t* stats;
	stage1_sieve_data_t* d;

	mpz_t ad;
	mpz_t p;
	mpz_t m;
} stage1_hit_data_t;

static void
stage1_hit_free(void *data, int threadid)
{
	stage1_hit_data_t *hit_data = (stage1_hit_data_t *)data;

	mpz_clear(hit_data->ad);
	mpz_clear(hit_data->p);
	mpz_clear(hit_data->m);
	free(hit_data);
}

static void
stage1_hit_run(void* data, int threadid)
{
	stage1_hit_data_t* hit_data = (stage1_hit_data_t*)data;
	void* extra = hit_data->callback_data;   /* single-worker fallback (-np1, S=1) */

	if (hit_data->d && hit_data->d->num_stage2_workers > 0)
		extra = &hit_data->d->stage2_workers[threadid].sizeopt_data;

	hit_data->callback(hit_data->ad, hit_data->p, hit_data->m, extra);

	if (hit_data->stats)
		poly_stats_stage2_done(hit_data->stats);
}

void
handle_collision(task_data_t *task, 
		uint64 p, uint64 special_q,
		uint128 special_q_root, int64 res)
{
	/* the proposed rational coefficient is p*special_q;
	   p and special_q must be coprime. The 'trivial
	   special q' has special_q = 1 and special_q_root = 0 */

	poly_coeff_t *c = task->c;

	uint64_2gmp(p, c->p);
	uint64_2gmp(special_q, c->tmp1);
	mpz_gcd(c->tmp2, c->p, c->tmp1);
	if (mpz_cmp_ui(c->tmp2, 1))
		return;

	mpz_mul(c->p, c->p, c->tmp1);

	/* the corresponding correction to trans_m0 is 
	   special_q_root + res * special_q^2, and can be
	   positive or negative */

	mpz_import(c->tmp3, 4, -1, sizeof(uint32), 0, 0, special_q_root.w);
	int64_2gmp(res, c->tmp2);

	mpz_mul(c->tmp1, c->tmp1, c->tmp1);
	mpz_addmul(c->tmp3, c->tmp2, c->tmp1);
	if (fabs(mpz_get_d(c->tmp3)) >
	    c->coeff_max / c->m0 * mpz_get_d(c->p) * mpz_get_d(c->p)) {

		return;
	} 
	mpz_add(c->m, c->trans_m0, c->tmp3);

	/* a lot can go wrong before this function is called!
	   Check that Kleinjung's modular condition is satisfied */

	mpz_pow_ui(c->tmp1, c->m, c->degree);
	mpz_mul(c->tmp2, c->p, c->p);
	mpz_sub(c->tmp1, c->trans_N, c->tmp1);
	mpz_tdiv_r(c->tmp3, c->tmp1, c->tmp2);
	if (mpz_cmp_ui(c->tmp3, 0)) {
		gmp_printf("crap %Zd %Zd %Zd\n", c->high_coeff, c->p, c->m);
		return;
	}

	/* the pair works, now translate the computed m back
	   to the original polynomial. We have

	   computed_m = degree * high_coeff * real_m +
	   			(second_highest_coeff) * p

	   and need to solve for real_m and second_highest_coeff.
	   Per the CADO code: reducing the above modulo
	   degree*high_coeff causes the first term on the right
	   to disappear, so second_highest_coeff can be found
	   modulo degree*high_coeff and real_m then follows */

	mpz_mul_ui(c->tmp1, c->high_coeff, c->degree);
	mpz_tdiv_r(c->tmp2, c->m, c->tmp1);
	mpz_invert(c->tmp3, c->p, c->tmp1);
	mpz_mul(c->tmp2, c->tmp3, c->tmp2);
	mpz_tdiv_r(c->tmp2, c->tmp2, c->tmp1);

	/* make second_highest_coeff as small as possible in
	   absolute value */

	mpz_tdiv_q_2exp(c->tmp3, c->tmp1, 1);
	if (mpz_cmp(c->tmp2, c->tmp3) > 0) {
		mpz_sub(c->tmp2, c->tmp2, c->tmp1);
	}

	/* solve for real_m */
	mpz_submul(c->m, c->tmp2, c->p);
	mpz_tdiv_q(c->m, c->m, c->tmp1);

	//if (task->d->test_mode) {
	//	if (task->d->test_dump)
	//		gmp_fprintf(task->d->test_dump, "%Zd %Zd %Zd\n",
	//			c->high_coeff, c->p, c->m);
	//	//return;                           /* skip stage 2 entirely */
	//}

	{
		/* submit the hit to the stage 2 thread pool */

		task_control_t task_control;
		stage1_sieve_data_t *d = task->d;
		stage1_hit_data_t *hit_data = (stage1_hit_data_t *)
					xmalloc(sizeof(stage1_hit_data_t));

		hit_data->callback = d->poly->callback;
		hit_data->callback_data = d->poly->callback_data;
		hit_data->stats = d->stats;
		hit_data->d = d;
		mpz_init_set(hit_data->ad, c->high_coeff);
		mpz_init_set(hit_data->p, c->p);
		mpz_init_set(hit_data->m, c->m);

		if (d->stats)
			poly_stats_add_hit(d->stats);

		task_control.init = NULL;
		task_control.run = stage1_hit_run;
		task_control.shutdown = stage1_hit_free;
		task_control.data = hit_data;

		threadpool_add_task(d->stage2_threadpool,
					&task_control, 1);
	}
}

/*------------------------------------------------------------------------*/
static void
poly_search_init(poly_search_t *poly, poly_stage1_t *data)
{
	mpz_init_set(poly->N, data->gmp_N);

	mpz_init_set(poly->gmp_high_coeff_begin, 
			data->gmp_high_coeff_begin);
	mpz_init_set(poly->gmp_high_coeff_end, 
			data->gmp_high_coeff_end);
	mpz_init(poly->tmp1);

	poly->degree = data->degree;
	poly->norm_max = data->norm_max;
	poly->callback = data->callback;
	poly->callback_data = data->callback_data;
}

static void
poly_search_free(poly_search_t *poly)
{
	mpz_clear(poly->N);
	mpz_clear(poly->gmp_high_coeff_begin);
	mpz_clear(poly->gmp_high_coeff_end);
	mpz_clear(poly->tmp1);
}

/*------------------------------------------------------------------------*/
poly_coeff_t *
poly_coeff_init(void)
{
	poly_coeff_t *c = (poly_coeff_t *)xmalloc(sizeof(poly_coeff_t));

	mpz_init(c->high_coeff);
	mpz_init(c->trans_N);
	mpz_init(c->trans_m0);
	mpz_init(c->m);
	mpz_init(c->p);
	mpz_init(c->tmp1);
	mpz_init(c->tmp2);
	mpz_init(c->tmp3);
	return c;
}

void
poly_coeff_free(poly_coeff_t *c)
{
	mpz_clear(c->high_coeff);
	mpz_clear(c->trans_N);
	mpz_clear(c->trans_m0);
	mpz_clear(c->m);
	mpz_clear(c->p);
	mpz_clear(c->tmp1);
	mpz_clear(c->tmp2);
	mpz_clear(c->tmp3);
	free(c);
}

void
poly_coeff_copy(poly_coeff_t *dest, poly_coeff_t *src)
{
	dest->degree = src->degree;
	dest->coeff_max = src->coeff_max;
	dest->m0 = src->m0;
	dest->p_size_max = src->p_size_max;

	mpz_set(dest->high_coeff, src->high_coeff);
	mpz_set(dest->trans_N, src->trans_N);
	mpz_set(dest->trans_m0, src->trans_m0);
}

/*------------------------------------------------------------------------*/
typedef struct {
	uint32 p, r;
	uint8 log_val;
} sieve_prime_t;

#define SIEVE_ARRAY_SIZE 8192

typedef struct {
	uint8 *sieve_array;
	sieve_prime_t *primes;
	uint32 num_primes;
	uint32 num_primes_alloc;
	uint32 curr_offset;
	uint32 high_coeff_multiplier; /* divides all a_d */
	uint32 high_coeff_power_limit;
} sieve_t;

static void
sieve_ad_block(sieve_t *sieve, poly_search_t *poly)
{
	uint32 i;
	uint32 log_target;
	double target;

	target = mpz_get_d(poly->gmp_high_coeff_begin) /
				sieve->high_coeff_multiplier;
	target = MIN(target, HIGH_COEFF_SIEVE_LIMIT);

	log_target = floor(log(target) / M_LN2 + 0.5);
	memset(sieve->sieve_array, (int)(log_target - 4),
			SIEVE_ARRAY_SIZE);

	for (i = 0; i < sieve->num_primes; i++) {
		uint32 p = sieve->primes[i].p;
		uint32 r = sieve->primes[i].r;
		uint8 log_val = sieve->primes[i].log_val;

		while (r < SIEVE_ARRAY_SIZE) {
			sieve->sieve_array[r] -= log_val;
			r += p;
		}
		sieve->primes[i].r = r - SIEVE_ARRAY_SIZE;
	}
}

/*------------------------------------------------------------------------*/
static int
find_next_ad(sieve_t *sieve, poly_search_t *poly, mpz_t next_coeff)
{
	uint32 i, j, p, k;
	double td_test;
	uint8 *sieve_array = sieve->sieve_array;

	while (1) {

		for (i = sieve->curr_offset; i < SIEVE_ARRAY_SIZE; i++) {

			if (!(sieve_array[i] & 0x80))
				continue;

			mpz_divexact_ui(poly->tmp1, poly->gmp_high_coeff_begin,
					sieve->high_coeff_multiplier);
			mpz_add_ui(poly->tmp1, poly->tmp1, i);
			mpz_mul_ui(next_coeff, poly->tmp1,
					sieve->high_coeff_multiplier);

			if (mpz_cmp(next_coeff, poly->gmp_high_coeff_end) > 0)
				break;

			/* trial divide the a_d and skip it if it
			   does not have enough small factors */

			td_test = ceil(mpz_get_d(poly->tmp1) /
						HIGH_COEFF_SIEVE_LIMIT);

			for (j = p = 0; j < PRECOMPUTED_NUM_PRIMES; j++) {
				p += prime_delta[j];

				if (p > HIGH_COEFF_PRIME_LIMIT)
					break;

				for (k = 0; k < sieve->high_coeff_power_limit; k++) {
					if (mpz_divisible_ui_p(poly->tmp1, p))
						mpz_divexact_ui(poly->tmp1, 
							poly->tmp1, p);
					else
						break;
				}
			}
			if (mpz_get_d(poly->tmp1) > td_test)
				continue;

			/* a_d is okay, search it */

			sieve->curr_offset = i + 1;
			return 0;
		}

		/* update lower bound for next sieve block */

		mpz_set_ui(poly->tmp1, SIEVE_ARRAY_SIZE);
		mpz_mul_ui(poly->tmp1, poly->tmp1, sieve->high_coeff_multiplier);
		mpz_add(poly->gmp_high_coeff_begin,
				poly->gmp_high_coeff_begin, poly->tmp1);

		if (mpz_cmp(poly->gmp_high_coeff_begin,
					poly->gmp_high_coeff_end) > 0)
			break;

		sieve->curr_offset = 0;
		sieve_ad_block(sieve, poly);
	}

	return 1;
}

/*------------------------------------------------------------------------*/
static void
init_ad_sieve(sieve_t *sieve, poly_search_t *poly)
{
	uint32 i, j, p;
	uint32 digits = mpz_sizeinbase(poly->N, 10);

	if (poly->degree == 4) {
		sieve->high_coeff_multiplier = 420;
		sieve->high_coeff_power_limit = 4;
	}
	else if (digits > 200) {
		sieve->high_coeff_multiplier = 120120;
		sieve->high_coeff_power_limit = 4;
	}
	else if (digits > 120) {
		sieve->high_coeff_multiplier = 60;
		sieve->high_coeff_power_limit = 2;
	}
	else {
		sieve->high_coeff_multiplier = 12;
		sieve->high_coeff_power_limit = 2;
	}

	sieve->num_primes = 0;
	sieve->num_primes_alloc = 100;
	sieve->primes = (sieve_prime_t *)xmalloc(sizeof(sieve_prime_t) *
						sieve->num_primes_alloc);
	sieve->sieve_array = (uint8 *)xmalloc(sizeof(uint8) *
						SIEVE_ARRAY_SIZE);

	mpz_divexact_ui(poly->tmp1, poly->gmp_high_coeff_begin,
			(mp_limb_t)sieve->high_coeff_multiplier);
	for (i = p = 0; i < PRECOMPUTED_NUM_PRIMES; i++) {
		uint32 power;
		uint8 log_val;

		p += prime_delta[i];
		if (p > HIGH_COEFF_PRIME_LIMIT)
			break;

		log_val = floor(log(p) / M_LN2 + 0.5);
		power = p;
		for (j = 0; j < sieve->high_coeff_power_limit; j++) {
			uint32 r = mpz_cdiv_ui(poly->tmp1, (mp_limb_t)power);

			if (sieve->num_primes >= sieve->num_primes_alloc) {
				sieve->num_primes_alloc *= 2;
				sieve->primes = (sieve_prime_t *)xrealloc(
					sieve->primes,
					sieve->num_primes_alloc *
						sizeof(sieve_prime_t));
			}

			sieve->primes[sieve->num_primes].p = power;
			sieve->primes[sieve->num_primes].r = r;
			sieve->primes[sieve->num_primes].log_val = log_val;
			sieve->num_primes++;

			if ((uint32)(-1) / power < p)
				break;

			power *= p;
		}
	}

	sieve->curr_offset = 0;
	sieve_ad_block(sieve, poly);
}

/*------------------------------------------------------------------------*/
static void
free_ad_sieve(sieve_t *sieve)
{
	free(sieve->primes);
	free(sieve->sieve_array);
}

/*------------------------------------------------------------------------*/
static void
stage1_sieve_data_init(stage1_sieve_data_t *d, 
		msieve_obj *obj, poly_search_t *poly)
{
	uint32 i;
	uint32 num_threads;
	thread_control_t thread_control;

	/* pick the collision engine for this run (registry gates by what was
	   compiled in; HAVE_CUDA still decides whether the GPU engines exist) */
	const stage1_engine_vtable_t *v = stage1_engine_select(obj);

	d->obj = obj;
	d->poly = poly;
	d->engine = v;

	//d->test_mode = 0;
	//mpz_init(d->test_ad);
	//d->test_dump = NULL;
	//if (obj->nfs_args != NULL) {
	//	const char* tmp;
	//	if ((tmp = strstr(obj->nfs_args, "test_ad=")) != NULL) {
	//		uint64 ad = strtoull(tmp + 8, NULL, 10);   /* stops at the space */
	//		mpz_set_ui(d->test_ad, ad);                /* a_d always fits a word */
	//		d->test_mode = 1;
	//	}
	//	if ((tmp = strstr(obj->nfs_args, "test_pmin=")) != NULL) d->test_pmin = strtoul(tmp + 10, NULL, 10);
	//	if ((tmp = strstr(obj->nfs_args, "test_pmax=")) != NULL) d->test_pmax = strtoul(tmp + 10, NULL, 10);
	//	if ((tmp = strstr(obj->nfs_args, "test_qmin=")) != NULL) d->test_qmin = strtoull(tmp + 10, NULL, 10);
	//	if ((tmp = strstr(obj->nfs_args, "test_qmax=")) != NULL) d->test_qmax = strtoull(tmp + 10, NULL, 10);
	//}
	//if (d->test_mode) {
	//	char fn[256];
	//	num_threads = 1;                         /* single writer -> no dump lock */
	//	sprintf(fn, "test_%s.hits", v->name);
	//	d->test_dump = fopen(fn, "w");
	//	logprintf(obj, "TEST MODE %s: a_d=%" PRIu64
	//		"  p[%u,%u]  q[%" PRIu64 ",%" PRIu64 "] -> %s\n",
	//		v->name, (uint64)mpz_get_ui(d->test_ad),
	//		d->test_pmin, d->test_pmax,
	//		d->test_qmin, d->test_qmax, fn);
	//}

	/* account for multiple threads; we allocate a thread pool
	   with a number of threads requested, where each thread
	   deals with a single leading coefficient. We also allocate
	   another thread pool with a single thread, that runs stage
	   2. Eventually the latter can be made more concurrent. */

	num_threads = MAX(1, obj->num_threads);
	if (v->max_threads)                       /* GPU engines cap at 4 */
		num_threads = MIN(v->max_threads, num_threads);
	d->num_threads = num_threads;

	d->hw_data = NULL;                        /* per-run device ctx, if any */
	if (v->sieve_data_init)
		d->hw_data = v->sieve_data_init(obj, num_threads, v->id);   /* +v->id */

	thread_control.init = v->thread_data_init;
	thread_control.shutdown = v->thread_data_free;
	thread_control.data = d;

	d->threads = (stage1_sieve_thread_data_t *)xcalloc(
					num_threads,
					sizeof(stage1_sieve_thread_data_t));

	for (i = 0; i < num_threads; i++) {
		d->threads[i].sieve_p_fb = sieve_fb_alloc();
		d->threads[i].sieve_q_fb = sieve_fb_alloc();
	}

	d->stage1_threadpool = threadpool_init(num_threads,
					MAX(10, num_threads),
					&thread_control);

	thread_control.init = NULL;
	thread_control.shutdown = NULL;
	thread_control.data = NULL;
	d->stage2_threadpool = threadpool_init(
		MAX(1, d->num_stage2_workers), 1000, &thread_control);

}

/*------------------------------------------------------------------------*/
void stage1_sieve_data_free(stage1_sieve_data_t *d)
{
	uint32 i;

	if (!(d->obj->flags & MSIEVE_FLAG_STOP_SIEVING)) {
		/* we're allowed to try to shut down gracefully */

		threadpool_drain(d->stage1_threadpool, 1);
	}

	/* shut down the stage 1 threadpool first, since
	   we don't want it feeding the stage 2 threadpool
	   after it has been freed */

	threadpool_free(d->stage1_threadpool);
	threadpool_drain(d->stage2_threadpool, 1);   /* process the stragglers */
	threadpool_free(d->stage2_threadpool);

	//if (d->test_dump)
	//	fclose(d->test_dump);
	//mpz_clear(d->test_ad);

	for (i = 0; i < d->num_threads; i++) {
		sieve_fb_free(d->threads[i].sieve_p_fb);
		sieve_fb_free(d->threads[i].sieve_q_fb);
	}
	free(d->threads);

	if (d->engine->sieve_data_free && d->hw_data)
		d->engine->sieve_data_free(d->hw_data);
}

/*------------------------------------------------------------------------*/
static void
search_coeff_core(task_data_t * task, uint32 threadid)
{
	msieve_obj *obj = task->obj;
	poly_coeff_t *c = task->c;
	stage1_sieve_data_t *d = task->d;
	uint32 degree = d->poly->degree;
	uint32 num_pieces;
	uint32 p_min, p_max;
	uint64 special_q_min, special_q_max;
	uint64 special_q_min2, special_q_max2;
	uint32 special_q_fb_max;

	/* Kleinjung shows that the third-to-largest algebraic
	   polynomial coefficient is of size approximately

	             (correction to m0) * m0
		    --------------------------
		    (leading rational coeff)^2
	
	   We have a bound 'coeff_max' on what this number is 
	   supposed to be, and we know m0 and an upper bound on 
	   the size of the leading rational coefficient P. Let 
	   P = p1*p2*q, where p1 and p2 are drawn from a fixed
	   set of candidates, and q (the 'special-q') is arbitrary
	   except that gcd(q,p1,p2)=1. Then the correction to
	   m0 is < q0 + 0.5 * q^2 * max(p1,p2)^2 so that

	   coeff_max   0.5 * q^2 * max(p1,p2)^2
	   --------- < ------------------------ 
	      m0         (q * min(p1,p2)^2)^2

	   if p_max = P_SCALE * p_min then

	             0.5 * m0 * P_SCALE^4
	   p_max^2 < --------------------
	                  coeff_max
	*/

	p_max = MIN(MAX_P, sqrt(c->p_size_max));
	p_max = MIN(p_max, P_SCALE * P_SCALE *
			sqrt(0.5 * c->m0 / c->coeff_max));
	p_min = MAX(1, p_max / P_SCALE);

	special_q_max = MIN(MAX_SPECIAL_Q, 
			    c->p_size_max / p_min / p_min);
	special_q_max = MAX(special_q_max, 1);
	special_q_min = 1;

	//if (d->test_mode) {
	//	p_min = d->test_pmin;
	//	p_max = d->test_pmax;
	//	special_q_min = d->test_qmin;
	//	special_q_max = d->test_qmax;
	//}

	/* set up the special q factory; special-q may have 
	   arbitrary factors, but many small factors are 
	   preferred since that will allow for many more roots
	   per special q, so we choose the factors to be as 
	   small as possible */

	special_q_fb_max = MIN(200000, special_q_max);
	sieve_fb_init(d->threads[threadid].sieve_q_fb, c,
			2, special_q_fb_max,
			1, degree,
			1);

	/* because special-q can have any factors, we require that
	   the progressions we generate use p that have somewhat
	   large factors. This minimizes the chance that a given
	   special-q has factors in common with many progressions
	   in the set */

	sieve_fb_init(d->threads[threadid].sieve_p_fb, c, 
			100, 5000,
			1, degree,
		       	0);

	/* large search problems can be randomized so that
	   multiple runs over the same range of leading
	   a_d will likely generate different results */

	num_pieces = 1;
	if ((special_q_max - special_q_min > 500000)) //(!d->test_mode) && 
		num_pieces = MIN(200, (double)special_q_max * p_max
				/ log(special_q_max) / log(p_max)
				/ 3e10);

	if (num_pieces > 51) { /* randomize the special_q range */
		uint32 piece_length = (special_q_max - special_q_min)
				/ num_pieces;
		uint32 piece = get_rand(&obj->seed1, &obj->seed2)
				% num_pieces;
        piece = 100;
		special_q_min2 = special_q_min + piece * piece_length;
		special_q_max2 = special_q_min2 + piece_length;
	}
	else {
		special_q_min2 = special_q_min;
		special_q_max2 = special_q_max;
	}

	// gmp_printf("coeff %Zd specialq %" PRId64 " - %" PRId64 " p %u - %u\n",
	// 		c->high_coeff,
	// 		special_q_min2, special_q_max2,
	// 		p_min, p_max);

	/* dispatch to the selected engine. The registry replaces the old
	   compile-time CPU/GPU switch; the same call shape serves all four.
	   Route around anything outside the engine's envelope rather than
	   letting a capped engine crash: an over-p a_d is skipped, an over-q
	   window is clamped in place. (a_d-level routing ultimately belongs
	   in search_coeffs; skipping here is the safe interim.) */
	{
		const stage1_engine_vtable_t *v = d->engine;
		uint64 q_min = special_q_min2;
		uint64 q_max = special_q_max2;

		if (!stage1_engine_cell_fits(v, p_max, &q_max)) {
			logprintf(obj, "stage1 %s: skipping coeff, p_max %u "
				"exceeds engine cap %u\n",
				v->name, p_max, v->envelope.max_p);
			return;
		}

		v->specialq(task, threadid, q_min, q_max, p_min, p_max);
	}
}

/*------------------------------------------------------------------------*/
/* infrastructure for submitting new leading
   coeffs to the stage 1 thread pool */

static void
task_data_free(void *data, int threadid)
{
	task_data_t *task_data = (task_data_t *)data;

	poly_coeff_free(task_data->c);
	free(task_data);
}

static void
task_data_run(void *data, int threadid)
{
	task_data_t *task = (task_data_t *)data;

	search_coeff_core(task, threadid);
}

static double search_coeff_async(stage1_sieve_data_t * d,
			poly_coeff_t *c, double coeff_deadline)
{
	/* submit a leading coefficient asynchronously to the
	   thread pool; we copy the coefficient so the input
	   one can be overwritten by calling code */

	uint32 i;
	task_control_t task_control;
	task_data_t *task_data = (task_data_t *)xmalloc(sizeof(task_data_t));
	poly_coeff_t *c2 = poly_coeff_init();
	double cumulative_elapsed = 0;

	poly_coeff_copy(c2, c);

	task_data->obj = d->obj;
	task_data->c = c2;
	task_data->d = d;
	task_data->coeff_deadline = coeff_deadline;

	task_control.init = NULL;
	task_control.run = task_data_run;
	task_control.shutdown = task_data_free;
	task_control.data = task_data;

	threadpool_add_task(d->stage1_threadpool, &task_control, 1);

	for (i = 0; i < d->num_threads; i++)
		cumulative_elapsed += d->threads[i].cumulative_elapsed;
	return cumulative_elapsed;
}

/*------------------------------------------------------------------------*/
static void
search_coeffs(stage1_sieve_data_t *d, uint32 deadline)
{
	double deadline_per_coeff;
	double cumulative_time = 0;
	uint64 ad_index = 0;
	sieve_t ad_sieve;
	poly_search_t *poly = d->poly;
	poly_coeff_t *c = poly_coeff_init();

	deadline_per_coeff = 8640000;

	//if (d->test_mode) {
	//	mpz_set(c->high_coeff, d->test_ad);
	//	stage1_bounds_update(poly, c);
	//	search_coeff_async(d, c, deadline_per_coeff);
	//	threadpool_drain(d->stage1_threadpool, 1);      /* finish it + its dumps */
	//	d->obj->flags |= MSIEVE_FLAG_STOP_SIEVING;       /* stop the re-dispatch */
	//	poly_coeff_free(c);
	//	fflush(d->test_dump);
	//	fclose(d->test_dump);
	//	exit(0);
	//	return;
	//}

	/* set up lower limit on a_d */

	init_ad_sieve(&ad_sieve, poly);

	mpz_sub_ui(poly->tmp1, poly->gmp_high_coeff_begin, 1);
	mpz_fdiv_q_ui(poly->tmp1, poly->tmp1, 
			ad_sieve.high_coeff_multiplier);
	mpz_add_ui(poly->tmp1, poly->tmp1, 1);
	mpz_mul_ui(poly->gmp_high_coeff_begin, poly->tmp1, 
			ad_sieve.high_coeff_multiplier);

	/* count the a_d in range for the progress fraction */
	if (d->stats) {
		sieve_t count_sieve;
		poly_coeff_t* cc = poly_coeff_init();
		mpz_t save_begin;
		uint64 ad_count = 0;

		mpz_init_set(save_begin, poly->gmp_high_coeff_begin);
		init_ad_sieve(&count_sieve, poly);
		while (find_next_ad(&count_sieve, poly, cc->high_coeff) == 0)
			ad_count++;
		free_ad_sieve(&count_sieve);
		poly_coeff_free(cc);
		mpz_set(poly->gmp_high_coeff_begin, save_begin);
		mpz_clear(save_begin);
		d->stats->ad_total = ad_count;
	}

	while (1) {
		/* we only use a_d which are composed of
		   many small prime factors, in order to
		   have lots of projective roots going
		   into stage 2 */

		if (find_next_ad(&ad_sieve, poly, c->high_coeff))
			break;

		if (d->stats) {
			char adbuf[64];
			gmp_snprintf(adbuf, sizeof(adbuf), "%Zd",
				c->high_coeff);
			poly_stats_set_ad(d->stats, adbuf, ++ad_index);
		}

		/* recalculate internal parameters used
		   for search */

		stage1_bounds_update(poly, c);

		/* execute search */

		cumulative_time = search_coeff_async(
					d, c, deadline_per_coeff);

		if (d->obj->flags & MSIEVE_FLAG_STOP_SIEVING)
			break;

		if (deadline && cumulative_time > deadline)
			break;
	}

	free_ad_sieve(&ad_sieve);
	poly_coeff_free(c);
}

/*------------------------------------------------------------------------*/
void
poly_stage1_init(poly_stage1_t *data,
		 stage1_callback_t callback, void *callback_data)
{
	memset(data, 0, sizeof(poly_stage1_t));
	mpz_init_set_ui(data->gmp_N, (mp_limb_t)0);
	mpz_init_set_ui(data->gmp_high_coeff_begin, (mp_limb_t)0);
	mpz_init_set_ui(data->gmp_high_coeff_end, (mp_limb_t)0);
	data->callback = callback;
	data->callback_data = callback_data;
}

/*------------------------------------------------------------------------*/
void
poly_stage1_free(poly_stage1_t *data)
{
	mpz_clear(data->gmp_N);
	mpz_clear(data->gmp_high_coeff_begin);
	mpz_clear(data->gmp_high_coeff_end);
}

/*------------------------------------------------------------------------*/
void
poly_stage1_run(msieve_obj *obj, poly_stage1_t *data)
{
	/* pass external configuration in and run the search */

	poly_search_t poly;
	stage1_sieve_data_t sieve_data;

	poly_search_init(&poly, data);

	sieve_data.stats = data->stats;
	sieve_data.stage2_workers = data->stage2_workers;
	sieve_data.num_stage2_workers = data->num_stage2_workers;

	stage1_sieve_data_init(&sieve_data, obj, &poly);

	search_coeffs(&sieve_data, data->deadline);

	stage1_sieve_data_free(&sieve_data);
	poly_search_free(&poly);
}
