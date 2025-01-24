/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: stage1_roots.c 826 2012-12-01 21:30:10Z jaysonking $
--------------------------------------------------------------------*/

#include <stage1.h>

#define P_PRIME_LIMIT 0xfffff000
#define MAX_P_FACTORS 7
#define MAX_POWERS 4

/* structure for building arithmetic progressions */

typedef struct {
	uint32 power;
	uint32 roots[MAX_POLYSELECT_DEGREE];
} aprog_power_t;

typedef struct {
	uint32 p;

	/* number of 'degree-th roots' of N (mod p) */
	uint32 num_roots;

	/* largest e for which p^e < 2^32 */
	uint32 max_power;

	/* power p^e_j and its 'degree-th roots' of N (mod p^e_j) */
	aprog_power_t powers[MAX_POWERS];

	/* maximum product of other p which may be
	   combined with this p */
	uint32 cofactor_max;
} aprog_t;

typedef struct {
	aprog_t *aprogs;
	uint32 num_aprogs;
	uint32 num_aprogs_alloc;
} aprog_list_t;

/* structures for finding arithmetic progressions via
   explicit enumeration */

typedef struct {
	/* number of distinct primes p_i in the current enum product */
	uint32 num_factors;

	/* index into aprog list for each distinct prime p_i */
	uint32 curr_factor[MAX_P_FACTORS + 1];

	/* exponent e_i of each prime p_i */
	uint32 curr_power[MAX_P_FACTORS + 1];

	/* the 'running product' of p_1^e_1..p_{i-1}^e_{i-1} for
	   each i, or 1 if i==1 */
	uint32 curr_prod[MAX_P_FACTORS + 1];

	/* number of roots in each 'running product' */
	uint32 curr_num_roots[MAX_P_FACTORS + 1];
} p_enum_t;

#define ALGO_ENUM  0x1
#define ALGO_PRIME 0x2

/* a factory for building arithmetic progressions */

typedef struct {
	uint32 num_roots_min;
	uint32 num_roots_max;
	uint32 avail_algos;
	uint32 fb_only;
	uint32 degree;
	uint32 p_min, p_max;
	uint64 roots[MAX_ROOTS];
	mpz_poly_t tmp_poly;

	aprog_list_t aprog_data;

	prime_sieve_t p_prime;

	p_enum_t p_enum;

	mpz_t p, pp, nmodpp, tmp1, tmp2, tmp3, gmp_root;
} sieve_fb_t;

/*------------------------------------------------------------------------*/
static uint32
lift_root_32(uint32 n, uint32 r, uint32 old_power, 
		uint32 p, uint32 d)
{
	/* given r, a d_th root of n mod old_power, compute
	   the corresponding root mod (old_power*p) via Hensel lifting */

	uint32 q;
	uint32 p2 = old_power * p;
	uint64 rsave = r;

	q = mp_modsub_1(n % p2, mp_expo_1(r, d, p2), p2) / old_power;
	r = mp_modmul_1(d, mp_expo_1(r % p, d - 1, p), p);
	r = mp_modmul_1(q, mp_modinv_1(r, p), p);
	return rsave + old_power * r;
}

/*------------------------------------------------------------------------*/
void
sieve_fb_free(void *s_in)
{
	sieve_fb_t *s = (sieve_fb_t *)s_in;
	aprog_list_t *list = &s->aprog_data;

	free_prime_sieve(&s->p_prime);
	free(list->aprogs);

	mpz_clear(s->p);
	mpz_clear(s->pp);
	mpz_clear(s->nmodpp);
	mpz_clear(s->tmp1);
	mpz_clear(s->tmp2);
	mpz_clear(s->tmp3);
	mpz_clear(s->gmp_root);
	mpz_poly_free(&s->tmp_poly);

	free(s);
}

/*------------------------------------------------------------------------*/
static uint32
get_prime_roots(poly_coeff_t *c, uint32 p, uint32 *roots,
		mpz_poly_t *tmp_poly)
{
	/* find all nonzero roots of (N' - x^d) mod p, where 
	   d is the desired polynomial degree and N' is the 
	   transformed version of N modulo p. We throw out roots
	   of zero because a zero root implies p divides the
	   degree or the leading algebraic poly coefficient, and
	   neither of these is allowed in later stages */

	uint32 i, high_coeff;

	mpz_tdiv_r_ui(tmp_poly->coeff[0], c->trans_N, p);

	if (mpz_cmp_ui(tmp_poly->coeff[0], 0) == 0) {
		/* when p divides trans_N, only a root of zero
		   exists, so skip this p */
		return 0;
	}
	for (i = 1; i < c->degree; i++)
		mpz_set_ui(tmp_poly->coeff[i], 0);

	tmp_poly->degree = i;
	mpz_set_ui(tmp_poly->coeff[i], p - 1);

	return poly_get_zeros(roots, tmp_poly, p, &high_coeff, 0);
}

/*------------------------------------------------------------------------*/
static void
sieve_add_aprog(sieve_fb_t *s, poly_coeff_t *c, uint32 p, 
		uint32 fb_roots_min, uint32 fb_roots_max)
{
	uint32 i, j, nmodp, num_roots;
	uint32 degree = c->degree;
	uint32 power, power_limit;
	uint32 roots[MAX_POLYSELECT_DEGREE];
	aprog_list_t *list = &s->aprog_data;
	aprog_t *a;

	/* p will be able to generate arithmetic progressions;
	   add it to our list of them... */

	if (list->num_aprogs == list->num_aprogs_alloc) {
		list->num_aprogs_alloc *= 2;
		list->aprogs = (aprog_t *)xrealloc(list->aprogs,
						list->num_aprogs_alloc *
						sizeof(aprog_t));
	}

	/* ...if trans_N has any degree_th roots mod p */

	a = list->aprogs + list->num_aprogs;
	a->p = p;
	num_roots = get_prime_roots(c, p, roots, &s->tmp_poly);

	if (num_roots == 0 ||
	    num_roots < fb_roots_min ||
	    num_roots > fb_roots_max)
		return;

	list->num_aprogs++;
	a->num_roots = num_roots;

	power = p;
	power_limit = (uint32)(-1) / p;
	for (i = 1; i < MAX_POWERS && power < power_limit; i++)
		power *= p;

	a->max_power = i;
	a->powers[0].power = p;

	for (i = 0; i < num_roots; i++)
		a->powers[0].roots[i] = roots[i];

	/* add powers of p as well */

	power = p;
	for (i = 1; i < a->max_power; i++) {

		power *= p;
		a->powers[i].power = power;

		nmodp = mpz_tdiv_ui(c->trans_N, (mp_limb_t)power);

		for (j = 0; j < num_roots; j++)
			a->powers[i].roots[j] = lift_root_32(nmodp,
						a->powers[i - 1].roots[j],
						a->powers[i - 1].power,
						p, degree);
	}
}

/*------------------------------------------------------------------------*/
void *
sieve_fb_alloc(void)
{
	sieve_fb_t *s = (sieve_fb_t *)xcalloc(1, sizeof(sieve_fb_t));
	aprog_list_t *aprog = &s->aprog_data;

	mpz_init(s->p);
	mpz_init(s->pp);
	mpz_init(s->nmodpp);
	mpz_init(s->tmp1);
	mpz_init(s->tmp2);
	mpz_init(s->tmp3);
	mpz_init(s->gmp_root);
	mpz_poly_init(&s->tmp_poly);

	aprog->num_aprogs_alloc = 500;
	aprog->aprogs = (aprog_t *)xmalloc(aprog->num_aprogs_alloc *
						sizeof(aprog_t));
	return s;
}

/*------------------------------------------------------------------------*/
void 
sieve_fb_init(void *s_in, poly_coeff_t *c,
		uint32 factor_min, uint32 factor_max,
		uint32 fb_roots_min, uint32 fb_roots_max,
		uint32 fb_only)
{
	uint32 p;
	sieve_fb_t *s = (sieve_fb_t *)s_in;
	aprog_list_t *aprog = &s->aprog_data;
	prime_sieve_t prime_sieve;

	s->degree = c->degree;
	s->fb_only = fb_only;

	factor_max = MIN(factor_max, 1000000);

	if (factor_max <= factor_min)
		return;

	aprog->num_aprogs = 0;
	init_prime_sieve(&prime_sieve, factor_min, factor_max);

	while (1) {
		p = get_next_prime(&prime_sieve);

		if (p > factor_max)
			break;

		sieve_add_aprog(s, c, p, fb_roots_min, fb_roots_max);
	}

	free_prime_sieve(&prime_sieve);
}

/*------------------------------------------------------------------------*/
void 
sieve_fb_reset(void *s_in, uint32 p_min, uint32 p_max,
		uint32 num_roots_min, uint32 num_roots_max)
{
	uint32 i;
	sieve_fb_t *s = (sieve_fb_t *)s_in;
	aprog_t *aprogs = s->aprog_data.aprogs;
	uint32 num_aprogs = s->aprog_data.num_aprogs;

	if (num_roots_max > MAX_ROOTS) {
		printf("error: num_roots_max exceeds %d\n", MAX_ROOTS);
		exit(-1);
	}

	s->p_min = p_min;
	s->p_max = p_max;
	s->num_roots_min = num_roots_min;
	s->num_roots_max = num_roots_max;
	s->avail_algos = 0;

	if (num_roots_max == 0) {
		/* nothing to do */
		return;
	}

	/* set up for finding arithmetic progressions by
	   enumerating combinations of small factors p_i
	   from the aprog list */

	if (num_aprogs > 0) {
		p_enum_t *p_enum = &s->p_enum; 

		s->avail_algos |= ALGO_ENUM;

		/* clear the enum state */

		p_enum->curr_factor[0] = 0;
		p_enum->num_factors = 0;
		p_enum->curr_power[0] = 0;
		p_enum->curr_prod[0] = 1;
		p_enum->curr_num_roots[0] = 1;

		for (i = 0; i < num_aprogs; i++) {
			aprog_t *a = aprogs + i;

			a->cofactor_max = p_max / a->p;
		}
	}

	/* set up for finding arithmetic progressions where
	   p is a large prime exceeding factor_max; do so if
	   that's allowed and the range of such p is higher
	   than what's already in the aprog list */

	if (s->fb_only == 0 &&
	    p_min < P_PRIME_LIMIT &&
	    s->degree >= num_roots_min) {
		uint32 last_p;

		if (num_aprogs == 0)
			last_p = 0;
		else
			last_p = aprogs[num_aprogs - 1].p;

		if (p_max > last_p) {
			s->avail_algos |= ALGO_PRIME;

			free_prime_sieve(&s->p_prime);
			init_prime_sieve(&s->p_prime, MAX(p_min, last_p + 1),
					 MIN(p_max, P_PRIME_LIMIT));
		}
	}
}

/*------------------------------------------------------------------------*/
static void
lift_roots(sieve_fb_t *s, poly_coeff_t *c, uint32 p, uint32 num_roots)
{
	/* we have num_roots arithmetic progressions mod p;
	   convert the progressions to be mod p^2, using
	   Hensel lifting, and then move the origin of the
	   result trans_m0 units to the left. */

	uint32 i;
	unsigned long degree = s->degree;

	mpz_set_ui(s->p, (unsigned long)p);
	mpz_mul(s->pp, s->p, s->p);
	mpz_tdiv_r(s->nmodpp, c->trans_N, s->pp);
	mpz_tdiv_r(s->tmp3, c->trans_m0, s->pp);

	for (i = 0; i < num_roots; i++) {

		uint64_2gmp(s->roots[i], s->gmp_root);

		mpz_powm_ui(s->tmp1, s->gmp_root, degree, s->pp);
		mpz_sub(s->tmp1, s->nmodpp, s->tmp1);
		if (mpz_cmp_ui(s->tmp1, (mp_limb_t)0) < 0)
			mpz_add(s->tmp1, s->tmp1, s->pp);
		mpz_tdiv_q(s->tmp1, s->tmp1, s->p);

		mpz_powm_ui(s->tmp2, s->gmp_root, degree-1, s->p);
		mpz_mul_ui(s->tmp2, s->tmp2, degree);
		mpz_invert(s->tmp2, s->tmp2, s->p);

		mpz_mul(s->tmp1, s->tmp1, s->tmp2);
		mpz_tdiv_r(s->tmp1, s->tmp1, s->p);
		mpz_addmul(s->gmp_root, s->tmp1, s->p);
		mpz_sub(s->gmp_root, s->gmp_root, s->tmp3);
		if (mpz_cmp_ui(s->gmp_root, (unsigned long)0) < 0)
			mpz_add(s->gmp_root, s->gmp_root, s->pp);

		s->roots[i] = gmp2uint64(s->gmp_root);
	}
}

/*------------------------------------------------------------------------*/
static uint32
combine_roots(sieve_fb_t *s, uint32 p,
		uint32 num_factors, uint32 p_i[MAX_P_FACTORS],
		uint32 num_roots[MAX_P_FACTORS],
		uint32 roots[MAX_P_FACTORS][MAX_POLYSELECT_DEGREE])
{
	/* given a composite p and its factors p_i,
	   combine the roots mod p_i into roots mod
	   p using the Chinese Remainder Theorem */

	uint32 i, i0, i1, i2, i3, i4, i5, i6;
	uint32 prod[MAX_P_FACTORS];
	uint64 accum[MAX_P_FACTORS + 1];

	if (num_factors == 1) { 
		/* no CRT needed */

		for (i = 0; i < num_roots[0]; i++)
			s->roots[i] = roots[0][i];

		return MIN(num_roots[0], s->num_roots_max);
	}

	/* fill in auxiliary CRT quantities */

	for (i = 0; i < num_factors; i++) {
		prod[i] = p / p_i[i];
		prod[i] *= mp_modinv_1(prod[i] % p_i[i], p_i[i]);
	}
	accum[i] = 0;

#if MAX_P_FACTORS > 7
#error "MAX_P_FACTORS exceeds 7"
#endif
	/* loop over all combinations of roots, changing one
	   root at a time. The accumulator value in the innermost
	   loop will exceed p by a few bits, so we need the
	   accum array to have wide integers. */

	i0 = i1 = i2 = i3 = i4 = i5 = i6 = i = 0;
	switch (num_factors) {
	case 7:
		for (i6 = num_roots[6] - 1; (int32)i6 >= 0; i6--) {
			accum[6] = accum[7] + (uint64)prod[6] * roots[6][i6];
	case 6:
		for (i5 = num_roots[5] - 1; (int32)i5 >= 0; i5--) {
			accum[5] = accum[6] + (uint64)prod[5] * roots[5][i5];
	case 5:
		for (i4 = num_roots[4] - 1; (int32)i4 >= 0; i4--) {
			accum[4] = accum[5] + (uint64)prod[4] * roots[4][i4];
	case 4:
		for (i3 = num_roots[3] - 1; (int32)i3 >= 0; i3--) {
			accum[3] = accum[4] + (uint64)prod[3] * roots[3][i3];
	case 3:
		for (i2 = num_roots[2] - 1; (int32)i2 >= 0; i2--) {
			accum[2] = accum[3] + (uint64)prod[2] * roots[2][i2];
	case 2:
		for (i1 = num_roots[1] - 1; (int32)i1 >= 0; i1--) {
			accum[1] = accum[2] + (uint64)prod[1] * roots[1][i1];

		for (i0 = num_roots[0] - 1; (int32)i0 >= 0; i0--) {
			accum[0] = accum[1] + (uint64)prod[0] * roots[0][i0];
			s->roots[i++] = accum[0] % p;

			if (i == s->num_roots_max)
				goto finished;
		}}}}}}}
	}

finished:
	return i;
}

/*------------------------------------------------------------------------*/
static uint32
get_enum_roots(sieve_fb_t *s)
{
	uint32 i, j;
	aprog_t *aprogs = s->aprog_data.aprogs;
	p_enum_t *p_enum = &s->p_enum;
	uint32 *curr_factor = p_enum->curr_factor;
	uint32 *curr_power = p_enum->curr_power;
	uint32 num_factors = p_enum->num_factors;
	uint32 p = p_enum->curr_prod[num_factors];
	uint32 p_i[MAX_P_FACTORS];
	uint32 num_roots[MAX_P_FACTORS];
	uint32 roots[MAX_P_FACTORS][MAX_POLYSELECT_DEGREE];

	for (i = 0; i < num_factors; i++) {
		aprog_t *a = aprogs + curr_factor[i];

		/* p_i may be a power */

		p_i[i] = a->powers[curr_power[i]].power;
		num_roots[i] = a->num_roots;
		for (j = 0; j < num_roots[i]; j++)
			roots[i][j] = a->powers[curr_power[i]].roots[j];
	}

	return combine_roots(s, p, num_factors, p_i, num_roots, roots);
}

/*------------------------------------------------------------------------*/
static uint32
get_next_enum(sieve_fb_t *s)
{
	/* find the next p by enumerating combinations of
	   aprogs p_i. We do this by maintaining a "stack"
	   of unique chosen p_i which are strictly in
	   non-increasing order, when read first-to-last.
	   While maintaining this ordering, and respecting
	   the maximum size of p, we "push" a new p_i
	   whenever possible. Otherwise, we "pop" p_i and,
	   as long as ordering can be maintained, "push"
	   the successor of the p_i which was last popped.
	   Multiple copies of the same p_i are combined
	   together into a single p_i (to simplify the CRT
	   later).

	   enumeration ends when there are no more aprogs
	   p_i available to be pushed

	   this method allows the p_i to be recovered
	   directly by reading the stack; no trial division
	   is necessary */

	p_enum_t *p_enum = &s->p_enum;
	uint32 *curr_factor = p_enum->curr_factor;
	uint32 *curr_power = p_enum->curr_power;
	uint32 *curr_prod = p_enum->curr_prod;
	uint32 *curr_num_roots = p_enum->curr_num_roots;

	while (1) {

		uint32 i = p_enum->num_factors;
		uint32 power_up = (i && curr_factor[i] == curr_factor[i - 1]);
		uint32 num_roots = curr_num_roots[i];
		aprog_t *a = NULL;

		if (curr_factor[i] < s->aprog_data.num_aprogs)
			a = s->aprog_data.aprogs + curr_factor[i];

		if (a != NULL && curr_prod[i] <= a->cofactor_max
		    && !(power_up && ++curr_power[i - 1] >= a->max_power)
		    && (power_up || i < MAX_P_FACTORS)) {

			uint32 p = curr_prod[i] * a->p;

			if (!power_up) {
				p_enum->num_factors = ++i;
				num_roots *= a->num_roots;
				curr_num_roots[i] = num_roots;
			}

			curr_prod[i] = p;
			curr_factor[i] = 0;
			curr_power[i] = 0;

			if (p >= s->p_min && num_roots >= s->num_roots_min)
				return p;
		}
		else if (i) {

			i--;
			curr_factor[i]++;
			curr_power[i] = 0;
			p_enum->num_factors = i;
		}
		else {
			return P_SEARCH_DONE;
		}
	}
}

/*------------------------------------------------------------------------*/
uint32
sieve_fb_next(void *s_in, poly_coeff_t *c,
		root_callback callback, void *extra)
{
	/* main external interface */

	uint32 i, p, num_roots;
	sieve_fb_t *s = (sieve_fb_t *)s_in;

	while (1) {
		if (s->avail_algos & ALGO_ENUM) {

			/* first attempt to find a p by 
			   combining smaller suitable aprogs p_i */

			p = get_next_enum(s);

			if (p == P_SEARCH_DONE) {
				s->avail_algos &= ~ALGO_ENUM;
				continue;
			}

			num_roots = get_enum_roots(s);
		}
		else if (s->avail_algos & ALGO_PRIME) {

			/* then try to find a large prime p */

			uint32 roots[MAX_POLYSELECT_DEGREE];

			p = get_next_prime(&s->p_prime);

			if (p >= s->p_max || p >= P_PRIME_LIMIT) {
				s->avail_algos &= ~ALGO_PRIME;
				continue;
			}

			num_roots = get_prime_roots(c, p, roots,
						    &s->tmp_poly);
			num_roots = MIN(num_roots, s->num_roots_max);

			if (num_roots == 0 ||
			    num_roots < s->num_roots_min)
				continue;

			for (i = 0; i < num_roots; i++)
				s->roots[i] = roots[i];
		}
		else {
			return P_SEARCH_DONE;
		}

		/* p found; generate all the arithmetic
		   progressions it allows and postprocess the
		   whole batch */

		lift_roots(s, c, p, num_roots);
		callback(p, num_roots, s->roots, extra);

		return p;
	}
}
