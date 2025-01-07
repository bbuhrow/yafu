/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: stage1.c 817 2012-11-11 14:58:29Z jasonp_sf $
--------------------------------------------------------------------*/

#include <stage1.h>

/* main driver for stage 1 */

/*------------------------------------------------------------------------*/
static void
stage1_bounds_update(poly_search_t *poly, poly_coeff_t *c)
{
	/* determine the parametrs for the collision search,
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
uint32
handle_collision(poly_coeff_t *c, uint64 p, uint32 special_q,
		uint64 special_q_root, int64 res)
{
	/* the proposed rational coefficient is p*special_q;
	   p and special_q must be coprime. The 'trivial
	   special q' has special_q = 1 and special_q_root = 0 */

	uint64_2gmp(p, c->p);
	mpz_gcd_ui(c->tmp1, c->p, special_q);
	if (mpz_cmp_ui(c->tmp1, 1))
		return 0;

	mpz_mul_ui(c->p, c->p, (unsigned long)special_q);

	/* the corresponding correction to trans_m0 is 
	   special_q_root + res * special_q^2, and can be
	   positive or negative */

	uint64_2gmp(special_q_root, c->tmp1);
	int64_2gmp(res, c->tmp2);
	mpz_set_ui(c->tmp3, special_q);

	mpz_mul(c->tmp3, c->tmp3, c->tmp3);
	mpz_addmul(c->tmp1, c->tmp2, c->tmp3);
	mpz_add(c->m, c->trans_m0, c->tmp1);

	/* a lot can go wrong before this function is called!
	   Check that Kleinjung's modular condition is satisfied */

	mpz_pow_ui(c->tmp1, c->m, c->degree);
	mpz_mul(c->tmp2, c->p, c->p);
	mpz_sub(c->tmp1, c->trans_N, c->tmp1);
	mpz_tdiv_r(c->tmp3, c->tmp1, c->tmp2);
	if (mpz_cmp_ui(c->tmp3, 0)) {
		gmp_printf("crap %Zd %Zd %Zd\n", c->high_coeff, c->p, c->m);
		return 0;
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
	return 1;
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
} sieve_t;

/*------------------------------------------------------------------------*/
static void
sieve_ad_block(sieve_t *sieve, poly_search_t *poly)
{
	uint32 i;
	uint32 log_target;
	double target;

	target = mpz_get_d(poly->gmp_high_coeff_begin) /
				HIGH_COEFF_MULTIPLIER;
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
					HIGH_COEFF_MULTIPLIER);
			mpz_add_ui(poly->tmp1, poly->tmp1, i);
			mpz_mul_ui(next_coeff, poly->tmp1,
					HIGH_COEFF_MULTIPLIER);

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

				for (k = 0; k < HIGH_COEFF_POWER_LIMIT; k++) {
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
		mpz_mul_ui(poly->tmp1, poly->tmp1, HIGH_COEFF_MULTIPLIER);
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

	sieve->num_primes = 0;
	sieve->num_primes_alloc = 100;
	sieve->primes = (sieve_prime_t *)xmalloc(sizeof(sieve_prime_t) *
						sieve->num_primes_alloc);
	sieve->sieve_array = (uint8 *)xmalloc(sizeof(uint8) *
						SIEVE_ARRAY_SIZE);

	mpz_divexact_ui(poly->tmp1, poly->gmp_high_coeff_begin,
			(mp_limb_t)HIGH_COEFF_MULTIPLIER);
	for (i = p = 0; i < PRECOMPUTED_NUM_PRIMES; i++) {
		uint32 power;
		uint8 log_val;

		p += prime_delta[i];
		if (p > HIGH_COEFF_PRIME_LIMIT)
			break;

		log_val = floor(log(p) / M_LN2 + 0.5);
		power = p;
		for (j = 0; j < HIGH_COEFF_POWER_LIMIT; j++) {
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
search_coeffs(msieve_obj *obj, poly_search_t *poly, uint32 deadline)
{
	uint32 digits = mpz_sizeinbase(poly->N, 10);
	double deadline_per_coeff;
	double cumulative_time = 0;
	sieve_t ad_sieve;
	poly_coeff_t *c = poly_coeff_init();
#ifdef HAVE_CUDA
	void *gpu_data = gpu_data_init(obj, poly);
#endif
	/* determine the CPU time limit; I have no idea if
	   the following is appropriate */

	if (digits <= 100)
		deadline_per_coeff = 5;
	else if (digits <= 105)
		deadline_per_coeff = 20;
	else if (digits <= 110)
		deadline_per_coeff = 30;
	else if (digits <= 120)
		deadline_per_coeff = 50;
	else if (digits <= 130)
		deadline_per_coeff = 100;
	else if (digits <= 140)
		deadline_per_coeff = 200;
	else if (digits <= 150)
		deadline_per_coeff = 400;
	else if (digits <= 175)
		deadline_per_coeff = 800;
	else if (digits <= 200)
		deadline_per_coeff = 1600;
	else
		deadline_per_coeff = 3200;

	printf("deadline: %.0lf CPU-seconds per coefficient\n",
					deadline_per_coeff);

	/* set up lower limit on a_d */

	mpz_sub_ui(poly->tmp1, poly->gmp_high_coeff_begin, 1);
	mpz_fdiv_q_ui(poly->tmp1, poly->tmp1, 
			HIGH_COEFF_MULTIPLIER);
	mpz_add_ui(poly->tmp1, poly->tmp1, 1);
	mpz_mul_ui(poly->gmp_high_coeff_begin, poly->tmp1, 
			HIGH_COEFF_MULTIPLIER);

	init_ad_sieve(&ad_sieve, poly);

	while (1) {
		double elapsed;

		/* we only use a_d which are composed of
		   many small prime factors, in order to
		   have lots of projective roots going
		   into stage 2 */

		if (find_next_ad(&ad_sieve, poly, c->high_coeff))
			break;

		/* recalculate internal parameters used
		   for search */

		stage1_bounds_update(poly, c);

		/* finally, sieve for polynomials using
		   Kleinjung's improved algorithm */

#ifdef HAVE_CUDA
		cumulative_time = sieve_lattice_gpu(obj, poly, c,
					gpu_data, deadline_per_coeff);
#else
		elapsed = sieve_lattice_cpu(obj, poly, c, deadline_per_coeff);
		cumulative_time += elapsed;
#endif

		if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)
			break;

		if (deadline && cumulative_time > deadline)
			break;
	}

	free_ad_sieve(&ad_sieve);
#ifdef HAVE_CUDA
	gpu_data_free(gpu_data);
#endif
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

	poly_search_init(&poly, data);

	search_coeffs(obj, &poly, data->deadline);

	poly_search_free(&poly);
}
