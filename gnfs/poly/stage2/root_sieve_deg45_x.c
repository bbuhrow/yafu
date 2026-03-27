/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: root_sieve_deg45_x.c 817 2012-11-11 14:58:29Z jasonp_sf $
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
} xpower_t;

typedef struct {
	uint32 num_powers;
	xpower_t *powers;
	uint32 sieve_size;
	uint16 *sieve;
} xdata_t;

typedef struct {
	uint16 score;
	uint32 which_y_block;
	uint32 which_lattice;
	uint64 resclass;
} xline_t;

#define XLINE_HEAP_SIZE 300

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
		uint32 which_y_block, uint32 which_lattice,
		uint64 resclass)
{
	xline_t *h = heap->entries;

	if (heap->num_entries <= XLINE_HEAP_SIZE - 1) {
		xline_t *s = h + heap->num_entries++;
		s->score = score;
		s->which_y_block = which_y_block;
		s->which_lattice = which_lattice;
		s->resclass = resclass;
		if (heap->num_entries == XLINE_HEAP_SIZE)
			make_heap(h, XLINE_HEAP_SIZE);
	}
	else if (h->score < score) {
		h->score = score;
		h->which_y_block = which_y_block;
		h->which_lattice = which_lattice;
		h->resclass = resclass;
		heapify(h, 0, XLINE_HEAP_SIZE);
	}
}

/*-------------------------------------------------------------------------*/
uint64 find_lattice_size_x(mpz_t prev_lattice_size,
				double line_length)
{
	uint32 i, p;
	uint64 curr_size = 1;
	double target_size = line_length / 
			mpz_get_d(prev_lattice_size) / 1e5;
	uint32 num_lattice_primes = 0;

	if (line_length < 2000000)
		return 1;

	for (i = p = 0; i < PRECOMPUTED_NUM_PRIMES; i++) {

		uint32 power;

		p += prime_delta[i];
		if (p >= MAX_SIEVE_PRIME)
			break;

		if (mpz_tdiv_ui(prev_lattice_size, p) == 0)
			continue;

		switch (p) {
		case 2: power = (target_size < 1000) ? 4 : 8; break;
		case 3: power = (target_size < 1000) ? 9 : 27; break;
		case 5: power = (target_size < 50000) ? 5 : 25; break;
		default: power = p; break;
		}

		if (curr_size * power >= target_size)
			break;

		curr_size *= power;
		if (++num_lattice_primes == MAX_CRT_FACTORS)
			break;
	}

	return curr_size;
}

/*-------------------------------------------------------------------------*/
static uint32
find_lattice_primes(sieve_prime_t *primes, uint32 num_primes,
			sieve_prime_t *lattice_primes, 
			uint64 lattice_size_x)
{
	uint32 i, num_lattice_primes;

	for (i = num_lattice_primes = 0; i < num_primes; i++) {

		sieve_prime_t *curr_prime = primes + i;
		uint32 p = curr_prime->prime;
		uint32 num_powers = 0;

		if (lattice_size_x == 1)
			return num_lattice_primes;

		if (lattice_size_x % p != 0) 
			continue;

		do {
			lattice_size_x /= p;
			num_powers++;
		} while (lattice_size_x % p == 0);

		lattice_primes[num_lattice_primes] = *curr_prime;
		lattice_primes[num_lattice_primes].num_powers = num_powers;
		num_lattice_primes++;
	}

	/* shouldn't get here */
	printf("error trial dividing lattice_size_x\n");
	return 0;
}

/*-------------------------------------------------------------------------*/
#define MAX_X_LATTICES 60

static void 
root_sieve_x(root_sieve_t *rs, xdata_t *xdata, 
		uint32 num_lattice_primes, uint32 which_y_block,
		xline_heap_t *heap, uint32 which_lattice,
		uint32 xy_score)
{
	uint32 i, j, k;
	hit_t hitlist[MAX_CRT_FACTORS];
	lattice_t lattices_x[MAX_X_LATTICES];
	sieve_x_t *x = &rs->xdata;
	uint32 num_lattices = 1;

	for (i = 0; i < num_lattice_primes; i++) {

		xdata_t *curr_xdata = xdata + i;
		uint32 sieve_size = curr_xdata->sieve_size;
		uint16 *sieve = curr_xdata->sieve;
		uint16 max_sieve_val = 0;
		hit_t *hits = hitlist + i;

		for (j = 0; j < sieve_size; j++) {
			if (sieve[j] > max_sieve_val) {
				max_sieve_val = sieve[j];
			}
		}
		max_sieve_val = 0.7 * max_sieve_val;

		for (j = k = 0; j < sieve_size; j++) {
			if (sieve[j] >= max_sieve_val) {
				if (k == MAX_ROOTS)
					break;

				hits->score[k] = sieve[j];
				hits->roots[k][0] = j;
				k++;
			}
		}

		hits->power = sieve_size;
		hits->num_roots = k;
		num_lattices *= k;
	}

	num_lattices = MIN(num_lattices, MAX_X_LATTICES);

	compute_lattices(hitlist, num_lattice_primes,
			lattices_x, x->lattice_size, num_lattices, 1);

	for (i = 0; i < num_lattices; i++) {
		uint32 score = xy_score + lattices_x[i].score;

		if (heap->num_entries < XLINE_HEAP_SIZE ||
		    score > heap->entries[0].score) {

			save_xline(heap, score, which_y_block, 
					which_lattice, lattices_x[i].x);
		}
	}
}

/*-------------------------------------------------------------------------*/
static void
xdata_alloc(sieve_prime_t *lattice_primes, 
		uint32 num_lattice_primes, 
		mpz_t mp_lattice_size,
		xdata_t *xdata)
{
	uint32 i, j;

	memset(xdata, 0, sizeof(xdata_t));

	for (i = 0; i < num_lattice_primes; i++) {

		sieve_prime_t *curr_prime = lattice_primes + i;
		xdata_t *curr_xdata = xdata + i;
		uint32 num_powers = curr_prime->num_powers;
		uint32 sieve_size = curr_prime->powers[num_powers - 1].power;

		curr_xdata->sieve_size = sieve_size;
		curr_xdata->sieve = (uint16 *)xmalloc(sieve_size *
						sizeof(uint16));
		curr_xdata->num_powers = num_powers;
		curr_xdata->powers = (xpower_t *)xmalloc(num_powers *
						sizeof(xpower_t));

		for (j = 0; j < num_powers; j++) {
			sieve_power_t *sp = curr_prime->powers + j;
			xpower_t *curr_xpower = curr_xdata->powers + j;
			uint32 power = sp->power;

			curr_xpower->p = power;
			curr_xpower->latsize_mod_p = mpz_tdiv_ui(
						mp_lattice_size, power);
			curr_xpower->num_roots = sp->num_roots;
			curr_xpower->contrib = sp->sieve_contrib;
			curr_xpower->roots = (xprog_t *)xmalloc(sp->num_roots * 
							sizeof(xprog_t));
		}
	}
}

/*-------------------------------------------------------------------------*/
static void
xdata_free(xdata_t *xdata, uint32 num_lattice_primes)
{
	uint32 i, j;

	for (i = 0; i < num_lattice_primes; i++) {

		xdata_t *curr_xdata = xdata + i;

		for (j = 0; j < curr_xdata->num_powers; j++) {
			xpower_t *curr_xpower = curr_xdata->powers + j;
	
			free(curr_xpower->roots);
		}
		free(curr_xdata->powers);
		free(curr_xdata->sieve);
	}
}

/*-------------------------------------------------------------------------*/
static void 
xdata_init(sieve_prime_t *lattice_primes, xdata_t *xdata,
		uint32 num_lattice_primes, 
		mpz_t y_base, uint64 resclass_y)
{
	uint32 i, j, k;

	for (i = 0; i < num_lattice_primes; i++) {

		sieve_prime_t *curr_prime = lattice_primes + i;
		xdata_t *curr_xdata = xdata + i;
		uint32 num_powers = curr_xdata->num_powers;

		for (j = 0; j < num_powers; j++) {

			sieve_power_t *sp = curr_prime->powers + j;
			xpower_t * curr_xpower = curr_xdata->powers + j;
			uint32 p = curr_xpower->p;
			uint32 num_roots = curr_xpower->num_roots;
			xprog_t *roots = curr_xpower->roots;

			uint32 latsize_mod_p = curr_xpower->latsize_mod_p;
			uint32 yres_mod_p = (uint32)(resclass_y % p);
			uint32 ybase_mod_p = mpz_fdiv_ui(y_base, p);
			uint32 y_mod_p = mp_modadd_1(yres_mod_p, 
							ybase_mod_p, p);

			for (k = 0; k < num_roots; k++) {
				sieve_root_t *r = sp->roots + k;

				roots[k].stride_y = mp_modmul_1(r->resclass, 
							latsize_mod_p, p);
				roots[k].start = mp_modsub_1(r->start, 
						mp_modmul_1(r->resclass,
							y_mod_p, p), p);
			}
		}
	}
}

/*-------------------------------------------------------------------------*/
static void
do_sieving(xdata_t *xdata)
{
	uint32 i, j, k;
	uint32 num_powers = xdata->num_powers;
	uint16 *sieve = xdata->sieve;
	uint32 sieve_size = xdata->sieve_size;

	memset(sieve, 0, sieve_size * sizeof(uint16));

	for (i = 0; i < num_powers; i++) {

		xpower_t * curr_xpower = xdata->powers + i;
		uint32 p = curr_xpower->p;
		uint32 num_roots = curr_xpower->num_roots;
		uint16 contrib = curr_xpower->contrib;
		xprog_t *roots = curr_xpower->roots;

		for (j = 0; j < num_roots; j++) {

			xprog_t *curr_prog = roots + j;
			uint32 start = curr_prog->start;

			for (k = start; k < sieve_size; k += p)
				sieve[k] += contrib;

			curr_prog->start = mp_modsub_1(start,
					curr_prog->stride_y, p);
		}
	}
}

/*-------------------------------------------------------------------------*/
static void 
find_hits(root_sieve_t *rs, xdata_t *xdata, 
		uint32 num_lattice_primes, uint32 y_blocks,
		xline_heap_t *heap, uint32 which_lattice,
		uint32 xy_score)
{
	uint32 i, j;

	for (i = 0; i < y_blocks; i++) {

		for(j = 0; j < num_lattice_primes; j++)
			do_sieving(xdata + j);

		root_sieve_x(rs, xdata, num_lattice_primes, 
				i, heap, which_lattice, xy_score);
	}
}

/*-------------------------------------------------------------------------*/
void
sieve_x_alloc(sieve_x_t *x)
{
	mpz_init(x->mp_lattice_size);
	mpz_init(x->resclass);
	mpz_init(x->crt0);
	mpz_init(x->crt1);
	mpz_init(x->tmp1);
}

/*-------------------------------------------------------------------------*/
void
sieve_x_free(sieve_x_t *x)
{
	mpz_clear(x->mp_lattice_size);
	mpz_clear(x->resclass);
	mpz_clear(x->crt0);
	mpz_clear(x->crt1);
	mpz_clear(x->tmp1);
}

/*-------------------------------------------------------------------------*/
static int 
compare_xlines(const void *x, const void *y)
{
	xline_t *xx = (xline_t *)x;
	xline_t *yy = (xline_t *)y;
	return (int)yy->score - (int)xx->score;
}

/*-------------------------------------------------------------------------*/
void
sieve_x_run_deg5(root_sieve_t *rs)
{
	uint32 i, j;
	sieve_xy_t *xy = &rs->xydata;
	sieve_x_t *x = &rs->xdata;
	sieve_prime_t *lattice_primes = x->lattice_primes;
	uint32 num_lattice_primes;
	msieve_obj *obj = rs->data->obj;

	double line_min, line_max;
	xdata_t xdata[MAX_CRT_FACTORS];
	xline_heap_t xline_heap;
	uint32 cutoff_score;

	line_min = line_max = 0;
	for (i = 0; i < xy->y_blocks; i++) {
		double curr_min = xy->x_line_min[i];
		double curr_max = xy->x_line_max[i];
		if (curr_max - curr_min > line_max - line_min) {
			line_min = curr_min;
			line_max = curr_max;
		}
	}
	if (line_min == line_max)
		return;

	x->lattice_size = find_lattice_size_x(xy->mp_lattice_size,
					line_max - line_min);

	num_lattice_primes = x->num_lattice_primes = 
				find_lattice_primes(rs->primes, 
					rs->num_primes, lattice_primes, 
					x->lattice_size);

	uint64_2gmp(x->lattice_size, x->tmp1);
	mpz_mul(x->mp_lattice_size, xy->mp_lattice_size, x->tmp1);
	x->dbl_lattice_size = mpz_get_d(x->mp_lattice_size);

	if (x->lattice_size == 1) {
		for (i = 0; i < xy->num_lattices; i++) {
			lattice_t * curr_lattice = xy->lattices + i;

			uint64_2gmp(curr_lattice->x, xy->resclass_x);
			uint64_2gmp(curr_lattice->y, xy->resclass_y);
			mpz_set(x->resclass, xy->resclass_x);
			mpz_add(rs->curr_y, xy->y_base, xy->resclass_y);

			x->curr_score = curr_lattice->score;

			for (j = 0; j < xy->y_blocks; j++) {
				line_min = xy->x_line_min[j];
				line_max = xy->x_line_max[j];

				x->apoly = rs->apoly;
				x->apoly.coeff[2] += mpz_get_d(rs->curr_y) * 
								rs->dbl_p;
				x->apoly.coeff[1] -= mpz_get_d(rs->curr_y) * 
								rs->dbl_d;

				mpz_set_d(x->x_base, line_min);
				mpz_tdiv_q(x->x_base, x->x_base, 
						x->mp_lattice_size);
				mpz_mul(x->x_base, x->x_base, 
						x->mp_lattice_size);
				x->x_blocks = (line_max - line_min) /
						x->dbl_lattice_size;

				root_sieve_line(rs);

				mpz_add(rs->curr_y, rs->curr_y, 
						xy->mp_lattice_size);
			}

			if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)
				break;
		}
		return;
	}

	xdata_alloc(lattice_primes, num_lattice_primes, 
			xy->mp_lattice_size, xdata);

	xline_heap.num_entries = 0;

	for (i = 0; i < xy->num_lattices; i++) {
		xdata_init(lattice_primes, xdata, num_lattice_primes,
				xy->y_base, xy->lattices[i].y);

		find_hits(rs, xdata, num_lattice_primes, 
				xy->y_blocks, &xline_heap, i,
				xy->lattices[i].score);
	}

	xdata_free(xdata, num_lattice_primes);

	qsort(xline_heap.entries, xline_heap.num_entries,
			sizeof(xline_t), compare_xlines);

	cutoff_score = 0.4 * xline_heap.entries[0].score;
	mpz_invert(x->crt0, xy->mp_lattice_size, x->tmp1);
	mpz_invert(x->crt1, x->tmp1, xy->mp_lattice_size);
	mpz_mul(x->crt0, x->crt0, xy->mp_lattice_size);
	mpz_mul(x->crt1, x->crt1, x->tmp1);

	for (i = 0; i < xline_heap.num_entries; i++) {

		xline_t *curr_xline = xline_heap.entries + i;
		lattice_t *curr_xy_lattice = xy->lattices + 
					curr_xline->which_lattice;

		if (curr_xline->score < cutoff_score)
			break;

		uint64_2gmp(curr_xy_lattice->x, xy->resclass_x);
		uint64_2gmp(curr_xy_lattice->y, xy->resclass_y);

		mpz_add(rs->curr_y, xy->y_base, xy->resclass_y);
		mpz_addmul_ui(rs->curr_y, xy->mp_lattice_size, 
				curr_xline->which_y_block);

		x->apoly = rs->apoly;
		x->apoly.coeff[2] += mpz_get_d(rs->curr_y) * rs->dbl_p;
		x->apoly.coeff[1] -= mpz_get_d(rs->curr_y) * rs->dbl_d;

		line_min = xy->x_line_min[curr_xline->which_y_block];
		line_max = xy->x_line_max[curr_xline->which_y_block];

		mpz_set_d(x->tmp1, line_min);
		mpz_tdiv_q(x->x_base, x->tmp1, x->mp_lattice_size);
		mpz_mul(x->x_base, x->x_base, x->mp_lattice_size);
		x->x_blocks = (line_max - line_min) / x->dbl_lattice_size;

		uint64_2gmp(curr_xline->resclass, x->tmp1);
		mpz_mul(x->resclass, x->tmp1, x->crt0);
		mpz_addmul(x->resclass, xy->resclass_x, x->crt1);

		mpz_tdiv_r(x->resclass, x->resclass, x->mp_lattice_size);

		x->curr_score = curr_xline->score;

		root_sieve_line(rs);

		if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)
			break;
	}
}

/*-------------------------------------------------------------------------*/
void
sieve_x_run_deg4(root_sieve_t *rs, uint64 lattice_size,
			double line_min, double line_max)
{
	uint32 i;
	sieve_xy_t *xy = &rs->xydata;
	sieve_x_t *x = &rs->xdata;
	sieve_prime_t *lattice_primes = x->lattice_primes;
	uint32 num_lattice_primes;
	msieve_obj *obj = rs->data->obj;

	xdata_t xdata[MAX_CRT_FACTORS];
	xline_heap_t xline_heap;
	uint32 cutoff_score;

	mpz_set_ui(xy->y_base, (unsigned long)0);
	xy->y_blocks = 0;
	mpz_set_ui(rs->curr_y, (unsigned long)0);
	rs->curr_z = 0;
	x->lattice_size = lattice_size;
	uint64_2gmp(x->lattice_size, x->mp_lattice_size);
	x->dbl_lattice_size = (double)x->lattice_size;
	x->apoly = rs->apoly;

	if (x->lattice_size == 1) {
		mpz_set_ui(x->resclass, (unsigned long)0);
		mpz_set_d(x->x_base, line_min);
		x->x_blocks = line_max - line_min;
		x->curr_score = 0;
		root_sieve_line(rs);
		return;
	}

	xline_heap.num_entries = 0;

	num_lattice_primes = x->num_lattice_primes = 
				find_lattice_primes(rs->primes, 
					rs->num_primes, lattice_primes, 
					x->lattice_size);

	xdata_alloc(lattice_primes, num_lattice_primes, 
			xy->mp_lattice_size, xdata);
	xdata_init(lattice_primes, xdata, num_lattice_primes, xy->y_base, 0);
	find_hits(rs, xdata, num_lattice_primes, 1, &xline_heap, 0, 0);
	xdata_free(xdata, num_lattice_primes);

	qsort(xline_heap.entries, xline_heap.num_entries,
			sizeof(xline_t), compare_xlines);

	cutoff_score = 0.4 * xline_heap.entries[0].score;

	for (i = 0; i < xline_heap.num_entries; i++) {

		xline_t *curr_xline = xline_heap.entries + i;

		if (curr_xline->score < cutoff_score)
			break;

		mpz_set_d(x->x_base, line_min);
		mpz_tdiv_q(x->x_base, x->x_base, x->mp_lattice_size);
		mpz_mul(x->x_base, x->x_base, x->mp_lattice_size);
		x->x_blocks = (line_max - line_min) / x->dbl_lattice_size;

		uint64_2gmp(curr_xline->resclass, x->resclass);

		x->curr_score = curr_xline->score;
		root_sieve_line(rs);

		if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)
			break;
	}
}
