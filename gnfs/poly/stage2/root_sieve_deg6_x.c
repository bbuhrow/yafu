/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: root_sieve_deg6_x.c 529 2011-01-17 14:39:34Z jasonp_sf $
--------------------------------------------------------------------*/

#include "stage2.h"

typedef struct {
	uint16 start;
	uint16 stride_y;
} xprog_t;

typedef struct {
	uint32 p;
	uint32 latsize_mod_p;
	uint32 num_roots;
	uint16 contrib;
	xprog_t *roots;
	uint16 *sieve;
} xdata_t;

typedef struct {
	uint16 score;
	uint32 which_y_block;
	uint64 resclass;
} xline_t;

#define XLINE_HEAP_SIZE 100

typedef struct {
	uint32 num_entries;
	xline_t entries[XLINE_HEAP_SIZE];
} xline_heap_t;

/*-------------------------------------------------------------------------*/
/* boilerplate code for managing heaps */

#define HEAP_SWAP(a,b) { tmp = a; a = b; b = tmp; }
#define HEAP_PARENT(i)  (((i)-1) >> 1)
#define HEAP_LEFT(i)    (2 * (i) + 1)
#define HEAP_RIGHT(i)   (2 * (i) + 2)

static void
heapify(xline_t *h, uint32 index, uint32 size) {

	uint32 c;
	xline_t tmp;
	for (c = HEAP_LEFT(index); c < (size-1); 
			index = c, c = HEAP_LEFT(index)) {

		if (h[c].score > h[c+1].score)
			c++;

		if (h[index].score > h[c].score) {
			HEAP_SWAP(h[index], h[c]);
		}
		else {
			return;
		}
	}
	if (c == (size-1) && h[index].score > h[c].score) {
		HEAP_SWAP(h[index], h[c]);
	}
}

static void
make_heap(xline_t *h, uint32 size) {

	int32 i;
	for (i = HEAP_PARENT(size); i >= 0; i--)
		heapify(h, (uint32)i, size);
}

static void
save_xline(xline_heap_t *heap, uint32 score, 
		uint32 which_y_block, uint64 resclass)
{
	xline_t *h = heap->entries;

	if (heap->num_entries <= XLINE_HEAP_SIZE - 1) {
		xline_t *s = h + heap->num_entries++;
		s->score = score;
		s->which_y_block = which_y_block;
		s->resclass = resclass;
		if (heap->num_entries == XLINE_HEAP_SIZE)
			make_heap(h, XLINE_HEAP_SIZE);
	}
	else if (h->score < score) {
		h->score = score;
		h->which_y_block = which_y_block;
		h->resclass = resclass;
		heapify(h, 0, XLINE_HEAP_SIZE);
	}
}

/*-------------------------------------------------------------------------*/
static uint32
find_lattice_primes(sieve_prime_t *primes, uint32 num_primes,
			mpz_t mp_lattice_size, sieve_prime_t *lattice_primes,
			uint64 *lattice_size_x, double line_length)
{
	uint32 i;
	uint32 num_lattice_primes = 0;
	uint64 tmp = 1;
	double target_size = line_length / mpz_get_d(mp_lattice_size) / 1e6;
	double curr_size = 1.0;

	for (i = 0; i < num_primes; i++) {

		sieve_prime_t *curr_prime = primes + i;

		if (num_lattice_primes == MAX_CRT_FACTORS)
			break;

		if (curr_prime->powers[0].num_roots == 0 ||
		    mpz_tdiv_ui(mp_lattice_size, curr_prime->prime) == 0)
			continue;

		if (curr_size * curr_prime->prime >= target_size)
			break;

		tmp *= curr_prime->prime;
		lattice_primes[num_lattice_primes] = *curr_prime;
		lattice_primes[num_lattice_primes].num_powers = 1;
		curr_size *= curr_prime->prime;
		num_lattice_primes++;
	}

	*lattice_size_x = tmp;
	return num_lattice_primes;
}

/*-------------------------------------------------------------------------*/
#define MAX_X_LATTICES 60

static void 
root_sieve_x(root_sieve_t *rs, xdata_t *xdata, 
		uint32 num_lattice_primes, uint32 which_y_block,
		xline_heap_t *heap, uint32 cutoff)
{
	uint32 i, j, k;
	hit_t hitlist[MAX_CRT_FACTORS];
	lattice_t lattices_x[MAX_X_LATTICES];
	sieve_x_t *x = &rs->xdata;
	uint32 num_lattices = 1;

	for (i = 0; i < num_lattice_primes; i++) {

		xdata_t *curr_xdata = xdata + i;
		uint32 p = curr_xdata->p;
		uint16 *sieve = curr_xdata->sieve;
		uint16 max_sieve_val = 0;
		hit_t *hits = hitlist + i;

		for (j = k = 0; j < p; j++) {
			if (sieve[j] >= max_sieve_val) {

				if (sieve[j] > max_sieve_val) {
					max_sieve_val = sieve[j];
					k = 0;
				}
				hits->score[k] = sieve[j];
				hits->roots[k][0] = j;
				k++;
			}
		}

		hits->power = p;
		hits->num_roots = k;
		num_lattices *= k;
	}

	num_lattices = MIN(num_lattices, MAX_X_LATTICES);

	compute_lattices(hitlist, num_lattice_primes,
			lattices_x, x->lattice_size, num_lattices, 1);

	if (lattices_x[0].score < cutoff)
		return;

	if (heap->num_entries < XLINE_HEAP_SIZE ||
	    lattices_x[0].score > heap->entries[0].score) {

		for (i = 0; i < num_lattices; i++) {
			save_xline(heap, lattices_x[i].score, 
					which_y_block, lattices_x[i].x);
		}
	}

}

/*-------------------------------------------------------------------------*/
static uint32
xdata_alloc(sieve_prime_t *lattice_primes, 
		uint32 num_lattice_primes, 
		mpz_t mp_lattice_size,
		xdata_t *xdata)
{
	uint32 i;
	uint32 cutoff = 0;

	for (i = 0; i < num_lattice_primes; i++) {

		sieve_prime_t *curr_prime = lattice_primes + i;
		xdata_t *curr_xdata = xdata + i;
		uint32 p = curr_prime->prime;
		uint32 num_roots = curr_prime->powers[0].num_roots;

		curr_xdata->p = p;
		curr_xdata->latsize_mod_p = mpz_tdiv_ui(mp_lattice_size, p);
		curr_xdata->num_roots = num_roots;
		curr_xdata->contrib = curr_prime->powers[0].sieve_contrib;
		curr_xdata->roots = (xprog_t *)xmalloc(num_roots * 
							sizeof(xprog_t));
		curr_xdata->sieve = (uint16 *)xmalloc(p * sizeof(uint16));

		cutoff += curr_xdata->contrib * MIN(3, num_roots - 1);
	}

	return cutoff;
}

/*-------------------------------------------------------------------------*/
static void
xdata_free(xdata_t *xdata, uint32 num_lattice_primes)
{
	uint32 i;

	for (i = 0; i < num_lattice_primes; i++) {

		xdata_t *curr_xdata = xdata + i;

		free(curr_xdata->roots);
		free(curr_xdata->sieve);
	}
}

/*-------------------------------------------------------------------------*/
static void 
xdata_init(sieve_prime_t *lattice_primes, xdata_t *xdata,
		uint32 num_lattice_primes, 
		mpz_t y_base, mpz_t resclass_y, int64 curr_z)
{
	uint32 i, j;

	for (i = 0; i < num_lattice_primes; i++) {

		sieve_prime_t *curr_prime = lattice_primes + i;
		sieve_power_t *sp = curr_prime->powers + 0;
		xdata_t *curr_xdata = xdata + i;

		uint32 p = curr_xdata->p;
		uint32 num_roots = curr_xdata->num_roots;
		xprog_t *roots = curr_xdata->roots;

		uint32 latsize_mod_p = curr_xdata->latsize_mod_p;
		uint32 yres_mod_p = mpz_tdiv_ui(resclass_y, p);
		uint32 ybase_mod_p = mpz_fdiv_ui(y_base, p);
		uint32 y_mod_p = mp_modadd_1(yres_mod_p, ybase_mod_p, p);
		int32 z_mod = curr_z % p;
		uint32 z_mod_p = (z_mod < 0) ?  (z_mod + (int32)p) : z_mod;

		for (j = 0; j < num_roots; j++) {
			sieve_root_t *r = sp->roots + j;
			xprog_t *curr_xprog = roots + j;
			uint32 start = r->start;
			uint32 resclass = r->resclass;
			uint32 resclass2 = mp_modmul_1(resclass, resclass, p);

			curr_xprog->stride_y = mp_modmul_1(resclass, 
							latsize_mod_p, p);
			start = mp_modsub_1(start, 
					mp_modmul_1(resclass, y_mod_p, p), 
					p);
			curr_xprog->start = mp_modsub_1(start, 
					mp_modmul_1(resclass2, z_mod_p, p), 
					p);
		}
	}
}

/*-------------------------------------------------------------------------*/
static void 
find_hits(root_sieve_t *rs, xdata_t *xdata, 
		uint32 num_lattice_primes, uint32 y_blocks,
		xline_heap_t *heap, uint32 cutoff)
{
	uint32 i, j, k;

	for (i = 0; i < y_blocks; i++) {

		for(j = 0; j < num_lattice_primes; j++) {

			xdata_t *curr_xdata = xdata + j;
			uint32 p = curr_xdata->p;
			uint32 num_roots = curr_xdata->num_roots;
			uint16 contrib = curr_xdata->contrib;
			xprog_t *roots = curr_xdata->roots;
			uint16 *sieve = curr_xdata->sieve;

			memset(sieve, 0, p * sizeof(uint16));

			for (k = 0; k < num_roots; k++) {

				xprog_t *curr_prog = roots + k;
				uint32 start = curr_prog->start;

				sieve[start] += contrib;

				curr_prog->start = mp_modsub_1(start,
						curr_prog->stride_y, p);
			}
		}

		root_sieve_x(rs, xdata, num_lattice_primes, 
				i, heap, cutoff);
	}
}

/*-------------------------------------------------------------------------*/
static int 
compare_xlines(const void *x, const void *y)
{
	xline_t *xx = (xline_t *)x;
	xline_t *yy = (xline_t *)y;
	return (int)yy->score - (int)xx->score;
}

void
sieve_x_run_deg6(root_sieve_t *rs)
{
	uint32 i;
	sieve_xy_t *xy = &rs->xydata;
	sieve_x_t *x = &rs->xdata;
	sieve_prime_t *lattice_primes = x->lattice_primes;
	uint32 num_lattice_primes;

	double direction[3] = {1, 0, 0};
	double line_min, line_max;
	xdata_t xdata[MAX_CRT_FACTORS];
	xline_heap_t xline_heap;
	uint32 cutoff_score;

	compute_line_size(rs->max_norm, &xy->apoly,
			rs->dbl_p, rs->dbl_d, direction,
			-10000, 10000, &line_min, &line_max);
	if (line_min > line_max)
		return;

	num_lattice_primes = x->num_lattice_primes = 
				find_lattice_primes(rs->primes, 
					rs->num_primes, xy->mp_lattice_size, 
					lattice_primes, &x->lattice_size,
					line_max - line_min);
	x->last_line_min = line_min;
	x->last_line_max = line_max;

	uint64_2gmp(x->lattice_size, x->tmp1);
	mpz_invert(x->crt0, xy->mp_lattice_size, x->tmp1);
	mpz_invert(x->crt1, x->tmp1, xy->mp_lattice_size);
	mpz_mul(x->crt0, x->crt0, xy->mp_lattice_size);
	mpz_mul(x->crt1, x->crt1, x->tmp1);
	mpz_mul(x->mp_lattice_size, xy->mp_lattice_size, x->tmp1);
	x->dbl_lattice_size = mpz_get_d(x->mp_lattice_size);

	cutoff_score = xdata_alloc(lattice_primes, 
				num_lattice_primes, 
				xy->mp_lattice_size, xdata);

	xdata_init(lattice_primes, xdata, num_lattice_primes,
			xy->y_base, xy->resclass_y, 
			rs->curr_z);

	xline_heap.num_entries = 0;
	find_hits(rs, xdata, num_lattice_primes, 
			xy->y_blocks, &xline_heap, 
			cutoff_score);

	xdata_free(xdata, num_lattice_primes);

	qsort(xline_heap.entries, xline_heap.num_entries,
			sizeof(xline_t), compare_xlines);
	cutoff_score = 0.95 * xline_heap.entries[0].score;

	for (i = 0; i < xline_heap.num_entries; i++) {

		xline_t *curr_xline = xline_heap.entries + i;

		if (curr_xline->score < cutoff_score)
			break;

		mpz_add(rs->curr_y, xy->y_base, xy->resclass_y);
		mpz_addmul_ui(rs->curr_y, xy->mp_lattice_size, 
				curr_xline->which_y_block);

		x->apoly = xy->apoly;
		x->apoly.coeff[2] += mpz_get_d(rs->curr_y) * rs->dbl_p;
		x->apoly.coeff[1] -= mpz_get_d(rs->curr_y) * rs->dbl_d;
		compute_line_size(rs->max_norm, &x->apoly,
				rs->dbl_p, rs->dbl_d, direction,
				x->last_line_min, x->last_line_max,
				&line_min, &line_max);
		if (line_min > line_max ||
		    line_max - line_min < 10000 * x->dbl_lattice_size)
			continue;

		mpz_set_d(x->tmp1, line_min);
		mpz_tdiv_q(x->x_base, x->tmp1, x->mp_lattice_size);
		mpz_mul(x->x_base, x->x_base, x->mp_lattice_size);
		x->x_blocks = (line_max - line_min) / x->dbl_lattice_size;

		uint64_2gmp(curr_xline->resclass, x->tmp1);
		mpz_mul(x->resclass, x->tmp1, x->crt0);
		mpz_addmul(x->resclass, xy->resclass_x, x->crt1);

		mpz_tdiv_r(x->resclass, x->resclass, x->mp_lattice_size);

		x->curr_score = xy->curr_score + curr_xline->score;

		root_sieve_line(rs);
	}
}
