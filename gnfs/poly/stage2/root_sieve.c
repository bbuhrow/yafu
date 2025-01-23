/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: root_sieve.c 817 2012-11-11 14:58:29Z jasonp_sf $
--------------------------------------------------------------------*/

#include "stage2.h"

#define MAX_SIEVE_PRIME_POWER 1000

/*-------------------------------------------------------------------------*/
static double
init_sieve_core(curr_poly_t *c, sieve_power_t *sp, uint32 deg, uint32 p)
{
	uint32 i, j;
	uint32 fi, mi, md, mp, step;
	uint32 co[MAX_POLY_DEGREE + 1];
	double contrib = 0;
	uint32 num_roots = 0;
	uint32 pk = sp->power;

	for (i = 0; i <= deg; i++)
		co[i] = mpz_fdiv_ui(c->gmp_a[i], (mp_limb_t)pk);
	md = mpz_fdiv_ui(c->gmp_d, (mp_limb_t)pk);
	mp = mpz_fdiv_ui(c->gmp_p, (mp_limb_t)pk);

	for (i = 0; i < pk; i++) {
		mi = mp_modsub_1(md, mp_modmul_1(i, mp, pk), pk);

		fi = co[deg];
		for (j = deg; j; j--)
			fi = (fi * i + co[j - 1]) % pk;

		if (mi != 0) {
			step = pk;
			while (mi % p == 0 && fi % p == 0) {
				mi /= p;
				fi /= p;
				step /= p;
			}
			if (mi % p) {
				sieve_root_t *r = sp->roots + num_roots++;
				r->resclass = i;
				r->step = step;
				r->start = mp_modmul_1(fi, 
						mp_modinv_1(mi, step), step);
			}	/* otherwise no solutions */
		}
		else if (fi == 0) {
			/* p^k | fi+j*(p*i-m) for all j */
			contrib += sp->root_contrib;
		}
	}

	sp->num_roots = num_roots;
	return contrib;
}

/*-------------------------------------------------------------------------*/
static void
init_sieve(curr_poly_t *c, root_sieve_t *rs, 
		uint32 deg, double sieve_bias)
{
	uint32 i, j;

	for (i = 0; i < rs->num_primes; i++) {
		sieve_prime_t *sp = rs->primes + i;

		for (j = 0; j < sp->num_powers; j++) {
			sieve_bias += init_sieve_core(c, sp->powers + j,
							deg, sp->prime);
		}
	}

	rs->sieve_bias = sieve_bias;
}


/*-------------------------------------------------------------------------*/
/* boilerplate code for managing heaps */

#define HEAP_SWAP(a,b) { tmp = a; a = b; b = tmp; }
#define HEAP_PARENT(i)  (((i)-1) >> 1)
#define HEAP_LEFT(i)    (2 * (i) + 1)
#define HEAP_RIGHT(i)   (2 * (i) + 2)

static void
heapify(rotation_t *h, uint32 index, uint32 size) {

	uint32 c;
	rotation_t tmp;
	for (c = HEAP_LEFT(index); c < (size-1); 
			index = c, c = HEAP_LEFT(index)) {

		if (h[c].score < h[c+1].score)
			c++;

		if (h[index].score < h[c].score) {
			HEAP_SWAP(h[index], h[c]);
		}
		else {
			return;
		}
	}
	if (c == (size-1) && h[index].score < h[c].score) {
		HEAP_SWAP(h[index], h[c]);
	}
}

static void
make_heap(rotation_t *h, uint32 size) {

	int32 i;
	for (i = HEAP_PARENT(size); i >= 0; i--)
		heapify(h, (uint32)i, size);
}

void
save_rotation(root_heap_t *heap, mpz_t x, mpz_t y,
		int64 z, float score)
{
	rotation_t *h = heap->entries;

	if (heap->num_entries <= heap->max_entries - 1) {
		rotation_t *r = h + heap->num_entries++;
		mpz_set(r->x, x);
		mpz_set(r->y, y);
		r->z = z;
		r->score = score;
		if (heap->num_entries == heap->max_entries)
			make_heap(h, heap->max_entries);
	}
	else if (h->score > score) {
		mpz_set(h->x, x);
		mpz_set(h->y, y);
		h->z = z;
		h->score = score;
		heapify(h, 0, heap->max_entries);
	}
}

/*-------------------------------------------------------------------------*/
typedef struct {
	double max_norm;
	double line_min;
	double line_max;
	uint64 lattice_size;
} bound_t;

#define MAX_BOUNDS 5

void
root_sieve_run_core(poly_rootopt_t *data, double initial_norm,
			double alpha_proj)
{
	uint32 i;
	msieve_obj *obj = data->obj;
	stage2_curr_data_t *s = (stage2_curr_data_t *)data->internal;
	root_sieve_t *rs = &s->root_sieve;
	curr_poly_t *c = &s->curr_poly;
	double alpha_bias = exp(-alpha_proj);
	uint32 num_bounds;
	bound_t bounds[MAX_BOUNDS];
	double norm_multiple = 1.01;
	double curr_norm;
	double line_min, line_max;
	double direction[3] = {0};
	double max_norm;

	rs->data = data;
	rs->root_heap.num_entries = 0;
	rs->dbl_p = mpz_get_d(c->gmp_p);
	rs->dbl_d = mpz_get_d(c->gmp_d);
	rs->apoly.degree = data->degree;
	for (i = 0; i <= data->degree; i++)
		rs->apoly.coeff[i] = mpz_get_d(c->gmp_a[i]);

	direction[data->degree - 4] = 1.0;
	rs->xyzdata.lattice_size = 1;
	mpz_set_ui(rs->xydata.mp_lattice_size, 1);
	mpz_set_ui(rs->xdata.mp_lattice_size, 1);
	line_min = -10;
	line_max = 10;
	max_norm = MIN(data->max_sizeopt_norm * alpha_bias, 100 * initial_norm);

	num_bounds = 0;
	for (curr_norm = 0; curr_norm < max_norm;
			norm_multiple *= 1.10) {

		uint64 lattice_size = 1;
		bound_t tmp_bound;

		curr_norm = initial_norm * norm_multiple;
		if (curr_norm > max_norm)
			curr_norm = max_norm;

		tmp_bound.max_norm = curr_norm;
		compute_line_size(tmp_bound.max_norm,
				&rs->apoly, rs->dbl_p, rs->dbl_d,
				direction, line_min, line_max,
				&line_min, &line_max);

		if (line_min >= line_max) {
			line_min = -10;
			line_max = 10;
			if (curr_norm == max_norm)
				break;
			else
				continue;
		}
		tmp_bound.line_min = line_min;
		tmp_bound.line_max = line_max;

		switch (data->degree) {
		case 4:
			lattice_size = 
				find_lattice_size_x(rs->xydata.mp_lattice_size,
						line_max - line_min);
			break;
		case 5:
			lattice_size = 
				find_lattice_size_y(line_max - line_min);
			break;
		case 6: 
			lattice_size = 
				find_lattice_size_z(line_max - line_min);
			break;
		}

		if (num_bounds == 0 ||
		    lattice_size > bounds[num_bounds - 1].lattice_size) {

			if (num_bounds < MAX_BOUNDS)
				num_bounds++;
		}

		tmp_bound.lattice_size = lattice_size;
		bounds[num_bounds - 1] = tmp_bound;
	}

	for (i = 0; i < num_bounds; i++) {
		bound_t *curr_bound = bounds + i;

		rs->max_norm = curr_bound->max_norm;
		switch (data->degree) {
		case 4:
			sieve_x_run_deg4(rs, curr_bound->lattice_size,
					curr_bound->line_min,
					curr_bound->line_max);
			break;
		case 5:
			sieve_xy_run_deg5(rs, curr_bound->lattice_size,
					curr_bound->line_min,
					curr_bound->line_max);
			break;
		case 6:
			sieve_xyz_run_deg6(rs, curr_bound->lattice_size,
					curr_bound->line_min,
					curr_bound->line_max);
			break;
		}

		if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)
			return;
	}

	for (i = 0; i < rs->root_heap.num_entries; i++) {
		rotation_t *r = rs->root_heap.entries + i;

		optimize_final(r->x, r->y, r->z, data);

		if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)
			return;
	}
}

/*-------------------------------------------------------------------------*/
void
root_sieve_run(poly_rootopt_t *data, double curr_norm, double alpha_proj)
{
	stage2_curr_data_t *s = (stage2_curr_data_t *)data->internal;
	curr_poly_t *c = &s->curr_poly;
	root_sieve_t *rs = &s->root_sieve;

	init_sieve(c, rs, data->degree, -alpha_proj);

	root_sieve_run_core(data, curr_norm, alpha_proj);
}

/*-------------------------------------------------------------------------*/
void root_sieve_init(root_sieve_t *rs)
{
	uint32 i, j;
	uint32 p;
	uint32 num_primes;

	memset(rs, 0, sizeof(root_sieve_t));

	rs->root_heap.max_entries = ROOT_HEAP_SIZE;
	rs->root_heap.entries = (rotation_t *)xmalloc(ROOT_HEAP_SIZE * 
							sizeof(rotation_t));
	for (i = 0; i < ROOT_HEAP_SIZE; i++) {
		mpz_init(rs->root_heap.entries[i].x);
		mpz_init(rs->root_heap.entries[i].y);
	}

	rs->sieve_block = (uint16 *)aligned_malloc(MAX(DEFAULT_BLOCK_SIZE,
						      MAX_SIEVE_PRIME_POWER) *
						   sizeof(uint16), 64);

	for (i = p = 0; i < PRECOMPUTED_NUM_PRIMES; i++) {
		p += prime_delta[i];
		if (p > MAX_SIEVE_PRIME)
			break;

		rs->random_root_score += log((double)p) / (p - 1);
	}

	num_primes = rs->num_primes = i;
	rs->primes = (sieve_prime_t *)xmalloc(num_primes * 
						sizeof(sieve_prime_t));

	for (i = p = 0; i < num_primes; i++) {

		uint32 num_powers;
		uint32 power;
		double dlog;
		sieve_prime_t *curr_prime = rs->primes + i;

		p += prime_delta[i];
		num_powers = 1;
		power = p;
		while (power * p < MAX_SIEVE_PRIME_POWER) {
			num_powers++;
			power *= p;
		}

		dlog = log((double)p) * p / (p + 1);
		curr_prime->prime = p;
		curr_prime->num_powers = num_powers;
		curr_prime->powers = (sieve_power_t *)xmalloc(num_powers *
						sizeof(sieve_power_t));
		curr_prime->contrib_array_size = power * 
					(MAX_SIEVE_PRIME_POWER / power);
		curr_prime->contrib_array = (uint16 *)aligned_malloc(
					UNROLL * sizeof(uint16) * 
					curr_prime->contrib_array_size, 64);

		for (j = 0, power = p; j < num_powers; j++, power *= p) {

			sieve_power_t *curr_power = curr_prime->powers + j;

			curr_power->power = power;
			curr_power->root_contrib = dlog / power;
			curr_power->sieve_contrib = (uint16)(0.5 +
						LOG_SCALE_FACTOR * 
						curr_power->root_contrib);
			curr_power->roots = (sieve_root_t *)xmalloc(power *
							sizeof(sieve_root_t));
		}
	}

	sieve_xyz_alloc(&rs->xyzdata);
	sieve_xy_alloc(&rs->xydata);
	sieve_x_alloc(&rs->xdata);
	mpz_init(rs->curr_x);
	mpz_init(rs->curr_y);
}

/*-------------------------------------------------------------------------*/
void root_sieve_free(root_sieve_t *rs)
{
	uint32 i, j;

	for (i = 0; i < ROOT_HEAP_SIZE; i++) {
		mpz_clear(rs->root_heap.entries[i].x);
		mpz_clear(rs->root_heap.entries[i].y);
	}

	for (i = 0; i < rs->num_primes; i++) {
		sieve_prime_t *prime = rs->primes + i;

		for (j = 0; j < prime->num_powers; j++) {
			free(prime->powers[j].roots);
		}
		aligned_free(prime->contrib_array);
	}
	free(rs->primes);
	free(rs->root_heap.entries);

	aligned_free(rs->sieve_block);

	sieve_xyz_free(&rs->xyzdata);
	sieve_xy_free(&rs->xydata);
	sieve_x_free(&rs->xdata);
	mpz_clear(rs->curr_x);
	mpz_clear(rs->curr_y);

	memset(rs, 0, sizeof(root_sieve_t));
}
