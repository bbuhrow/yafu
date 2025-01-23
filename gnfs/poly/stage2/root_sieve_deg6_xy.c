/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: root_sieve_deg6_xy.c 585 2011-06-20 04:16:29Z jaysonking $
--------------------------------------------------------------------*/

#include "stage2.h"

typedef struct {
	uint8 base_start;
	uint8 resclass;
	uint8 start;
	uint8 stride_z;
	uint8 *invtable_y;
} xyprog_t;

typedef struct {
	uint32 power;
	uint32 latsize_mod;
	uint16 contrib;
	uint32 num_roots;
	xyprog_t *roots;
} xypower_t;

typedef struct {
	uint32 p;
	uint32 table_size;
	uint32 num_powers;
	xypower_t *powers;
	uint16 *sieve;
} xydata_t;

typedef struct {
	lattice_t plane;
	uint32 which_z_block;
	uint32 which_lattice_xyz;
} plane_t;

#define PLANE_HEAP_SIZE 50

typedef struct {
	uint32 num_entries;
	plane_t entries[PLANE_HEAP_SIZE];
} plane_heap_t;

/*-------------------------------------------------------------------------*/
/* boilerplate code for managing heaps */

#define HEAP_SWAP(a,b) { tmp = a; a = b; b = tmp; }
#define HEAP_PARENT(i)  (((i)-1) >> 1)
#define HEAP_LEFT(i)    (2 * (i) + 1)
#define HEAP_RIGHT(i)   (2 * (i) + 2)

static void
heapify(plane_t *h, uint32 index, uint32 size) {

	uint32 c;
	plane_t tmp;
	for (c = HEAP_LEFT(index); c < (size-1); 
			index = c, c = HEAP_LEFT(index)) {

		if (h[c].plane.score > h[c+1].plane.score)
			c++;

		if (h[index].plane.score > h[c].plane.score) {
			HEAP_SWAP(h[index], h[c]);
		}
		else {
			return;
		}
	}
	if (c == (size-1) && h[index].plane.score > h[c].plane.score) {
		HEAP_SWAP(h[index], h[c]);
	}
}

static void
make_heap(plane_t *h, uint32 size) {

	int32 i;
	for (i = HEAP_PARENT(size); i >= 0; i--)
		heapify(h, (uint32)i, size);
}

static void
save_plane(plane_heap_t *heap, lattice_t *lattice, 
		uint32 which_z_block, uint32 which_lattice_xyz)
{
	plane_t *h = heap->entries;

	if (heap->num_entries <= PLANE_HEAP_SIZE - 1) {
		plane_t *s = h + heap->num_entries++;
		s->plane = *lattice;
		s->which_z_block = which_z_block;
		s->which_lattice_xyz = which_lattice_xyz;
		if (heap->num_entries == PLANE_HEAP_SIZE)
			make_heap(h, PLANE_HEAP_SIZE);
	}
	else if (h->plane.score < lattice->score) {
		h->plane = *lattice;
		h->which_z_block = which_z_block;
		h->which_lattice_xyz = which_lattice_xyz;
		heapify(h, 0, PLANE_HEAP_SIZE);
	}
}

/*-------------------------------------------------------------------------*/
static uint32
find_lattice_primes(sieve_prime_t *primes, uint32 num_primes,
			uint64 lattice_size_xyz, sieve_prime_t *lattice_primes,
			uint64 *lattice_size_xy, double line_length)
{
	uint32 i;
	uint32 num_lattice_primes = 0;
	uint64 tmp = 1;
	double target_size = line_length / (double)lattice_size_xyz / 1e4;
	double curr_size = 1.0;

	for (i = 0; i < num_primes; i++) {

		sieve_prime_t *curr_prime = primes + i;
		uint32 p = curr_prime->prime;
		uint32 num_powers = 1;

		if (num_lattice_primes == MAX_CRT_FACTORS)
			break;

		if (curr_prime->powers[0].num_roots == 0 ||
		    lattice_size_xyz % p == 0)
			continue;

		if (curr_size * p >= target_size)
			break;

		switch (p) {
		case 2:
			if (target_size < 1000) {
				p = 16;
				num_powers = 4;
			}
			else {
				p = 32;
				num_powers = 5;
			}
			break;

		case 3:
			if (target_size < 1000) {
				p = 9;
				num_powers = 2;
			}
			else {
				p = 27;
				num_powers = 3;
			}
			break;

		case 5:
			if (target_size > 50000) {
				p = 25;
				num_powers = 2;
			}
			break;
		}

		tmp *= p;
		lattice_primes[num_lattice_primes] = *curr_prime;
		lattice_primes[num_lattice_primes].num_powers = num_powers;
		curr_size *= p;
		num_lattice_primes++;
	}

	*lattice_size_xy = tmp;
	return num_lattice_primes;
}

/*-------------------------------------------------------------------------*/
#define MAX_XY_LATTICES 32

static void 
root_sieve_xy(root_sieve_t *rs, xydata_t *xydata, 
		uint32 num_lattice_primes, 
		uint32 which_z_block, uint32 which_lattice_xyz,
		plane_heap_t *plane_heap)
{
	uint32 i, j, k;
	hit_t hitlist[MAX_CRT_FACTORS];
	lattice_t lattices_xy[MAX_XY_LATTICES];
	sieve_xy_t *xy = &rs->xydata;
	uint32 num_lattices = 1;

	for (i = 0; i < num_lattice_primes; i++) {

		xydata_t *curr_xydata = xydata + i;
		uint32 p = curr_xydata->p;
		uint32 table_size = curr_xydata->table_size;
		uint16 *sieve = curr_xydata->sieve;
		uint32 max_sieve_val = 0;
		hit_t *hits = hitlist + i;

		for (j = k = 0; j < table_size; j++) {
			if (sieve[j] >= max_sieve_val) {

				if (sieve[j] > max_sieve_val) {
					max_sieve_val = sieve[j];
					k = 0;
				}

				hits->score[k] = sieve[j];
				hits->roots[k][0] = j % p;
				hits->roots[k][1] = j / p;
				if (++k == MAX_ROOTS)
					break;
			}
		}

		hits->power = p;
		hits->num_roots = k;
		num_lattices *= k;
	}

	num_lattices = MIN(num_lattices, MAX_XY_LATTICES);

	compute_lattices(hitlist, num_lattice_primes,
			lattices_xy, xy->lattice_size, 
			num_lattices, 2);

	if (plane_heap->num_entries < PLANE_HEAP_SIZE ||
	    lattices_xy[0].score > plane_heap->entries[0].plane.score) {

		for (i = 0; i < num_lattices; i++) {
			save_plane(plane_heap, lattices_xy + i,
					which_z_block, 
					which_lattice_xyz);
		}
	}
}

/*-------------------------------------------------------------------------*/
static void
xydata_alloc(sieve_prime_t *lattice_primes, 
		uint32 num_lattice_primes, 
		uint64 lattice_size_xyz,
		xydata_t *xydata)
{
	uint32 i, j, k;

	for (i = 0; i < num_lattice_primes; i++) {

		sieve_prime_t *curr_prime = lattice_primes + i;
		xydata_t *curr_xydata = xydata + i;
		uint32 num_powers = curr_prime->num_powers;
		uint32 p = curr_prime->powers[num_powers-1].power;
		uint32 table_size = p * p;

		curr_xydata->p = p;
		curr_xydata->table_size = table_size;
		curr_xydata->num_powers = num_powers;
		curr_xydata->powers = (xypower_t *)xmalloc(num_powers *
						sizeof(xypower_t));
		curr_xydata->sieve = (uint16 *)xmalloc(table_size *
						sizeof(uint16));

		for (j = 0; j < num_powers; j++) {
			sieve_power_t *sp = curr_prime->powers + j;
			xypower_t *curr_xypower = curr_xydata->powers + j;
			uint32 power = sp->power;
			uint32 num_roots = sp->num_roots;
	
			curr_xypower->power = power;
			curr_xypower->latsize_mod = lattice_size_xyz % power;
			curr_xypower->num_roots = num_roots;
			curr_xypower->contrib = sp->sieve_contrib;
			curr_xypower->roots = (xyprog_t *)xmalloc(num_roots * 
							sizeof(xyprog_t));

			for (k = 0; k < num_roots; k++) {
				sieve_root_t *r = sp->roots + k;
				xyprog_t *curr_xyprog = curr_xypower->roots + k;

				curr_xyprog->base_start = r->start;
				curr_xyprog->resclass = r->resclass;
				curr_xyprog->invtable_y = 
					(uint8 *)xmalloc(power * sizeof(uint8));
			}
		}
	}
}

/*-------------------------------------------------------------------------*/
static void
xydata_free(xydata_t *xydata, uint32 num_lattice_primes)
{
	uint32 i, j, k;

	for (i = 0; i < num_lattice_primes; i++) {

		xydata_t *curr_xydata = xydata + i;

		for (j = 0; j < curr_xydata->num_powers; j++) {
			xypower_t *curr_xypower = curr_xydata->powers + j;
	
			for (k = 0; k < curr_xypower->num_roots; k++)
				free(curr_xypower->roots[k].invtable_y);

			free(curr_xypower->roots);
		}
		free(curr_xydata->powers);
		free(curr_xydata->sieve);
	}
}

/*-------------------------------------------------------------------------*/
static void 
xydata_init(xydata_t *xydata, uint32 num_lattice_primes, 
		lattice_t *lattice_xyz, int64 z_base)
{
	uint32 i, j, k, m, n;

	for (i = 0; i < num_lattice_primes; i++) {

		xydata_t *curr_xydata = xydata + i;
		uint32 num_powers = curr_xydata->num_powers;

		for (j = 0; j < num_powers; j++) {
			xypower_t *curr_xypower = curr_xydata->powers + j;
			uint32 num_roots = curr_xypower->num_roots;

			uint32 p = curr_xypower->power;
			uint32 latsize_mod = curr_xypower->latsize_mod;
			uint32 y_mod_p = lattice_xyz->y % p;
			int64 z_start = z_base + lattice_xyz->z;
			int32 z_start_mod = z_start % p;
			uint32 z_mod_p = (z_start_mod < 0) ? 
					(z_start_mod + (int32)p) : z_start_mod;

			for (k = 0; k < num_roots; k++) {
				xyprog_t *curr_xyprog = curr_xypower->roots + k;

				uint8 *invtable_y = curr_xyprog->invtable_y;
				uint32 start = curr_xyprog->base_start;
				uint32 resclass = curr_xyprog->resclass;
				uint32 resclass2 = mp_modmul_1(resclass, 
							resclass, p);
				uint32 ytmp = y_mod_p;
				uint32 stride_y = mp_modmul_1(resclass, 
							latsize_mod, p);

				curr_xyprog->stride_z = mp_modmul_1(resclass2, 
							latsize_mod, p);

				start = mp_modsub_1(start, 
						mp_modmul_1(resclass, 
							y_mod_p, p), p);
				curr_xyprog->start = mp_modsub_1(start, 
						mp_modmul_1(resclass2, 
							z_mod_p, p), p);

				for (m = n = 0; m < p; m++) {
					invtable_y[ytmp] = n;
					ytmp = mp_modadd_1(ytmp, 
							latsize_mod, p);
					n = mp_modadd_1(n, stride_y, p);
				}
			}
		}
	}
}

/*-------------------------------------------------------------------------*/
static void 
do_sieving(xydata_t *curr_xydata) 
{
	uint32 i, j;

	uint16 *sieve = curr_xydata->sieve;
	xypower_t *curr_xypower = curr_xydata->powers + 0;
	uint32 p = curr_xypower->power;
	uint16 contrib = curr_xypower->contrib;
	uint32 num_roots = curr_xypower->num_roots;
	xyprog_t *roots = curr_xypower->roots;

	for (i = 0; i < num_roots; i++) {

		uint16 *row = sieve;
		xyprog_t *curr_prog = roots + i;
		uint32 start = curr_prog->start;
		uint8 *invtable_y = curr_prog->invtable_y;

		for (j = 0; j < p; j++) {
			uint32 curr_start = mp_modsub_1(start, 
						invtable_y[j], p);
			row[curr_start] += contrib;
			row += p;
		}

		curr_prog->start = mp_modsub_1(start, 
					curr_prog->stride_z, p);
	}
}

/*-------------------------------------------------------------------------*/
static void 
do_sieving_powers(xydata_t *curr_xydata) 
{
	uint32 i, j, k, m;

	uint16 *sieve = curr_xydata->sieve;
	uint32 p = curr_xydata->p;
	uint32 num_powers = curr_xydata->num_powers;
	xypower_t *powers = curr_xydata->powers;

	for (i = 0; i < num_powers; i++) {
		xypower_t *curr_xypower = powers + i;
		uint32 power = curr_xypower->power;
		uint16 contrib = curr_xypower->contrib;
		uint32 num_roots = curr_xypower->num_roots;
		xyprog_t *roots = curr_xypower->roots;

		for (j = 0; j < num_roots; j++) {

			uint16 *row = sieve;
			xyprog_t *curr_prog = roots + j;
			uint32 start = curr_prog->start;
			uint8 *invtable_y = curr_prog->invtable_y;

			for (k = 0; k < p; k += power) {

				for (m = 0; m < power; m++) {
					uint32 curr_start = mp_modsub_1(start, 
								invtable_y[m], 
								power);
					do {
						row[curr_start] += contrib;
						curr_start += power;
					} while (curr_start < p);

					row += p;
				}
			}

			curr_prog->start = mp_modsub_1(start, 
						curr_prog->stride_z, p);
		}
	}
}

/*-------------------------------------------------------------------------*/
static void 
find_hits(root_sieve_t *rs, xydata_t *xydata, 
		uint32 num_lattice_primes, 
		uint32 which_lattice_xyz,
		plane_heap_t *plane_heap)
{
	uint32 i, j;
	sieve_xyz_t *xyz = &rs->xyzdata;
	uint32 z_blocks = xyz->z_blocks;

	for (i = 0; i < z_blocks; i++) {

		for(j = 0; j < num_lattice_primes; j++) {

			xydata_t *curr_xydata = xydata + j;
			uint32 table_size = curr_xydata->table_size;
			uint16 *sieve = curr_xydata->sieve;

			memset(sieve, 0, table_size * sizeof(uint16));

			if (curr_xydata->num_powers == 1)
				do_sieving(curr_xydata);
			else
				do_sieving_powers(curr_xydata);
		}

		if (xyz->y_line_min[i] < xyz->y_line_max[i]) {
			root_sieve_xy(rs, xydata, 
					num_lattice_primes,
					i, which_lattice_xyz,
					plane_heap);
		}
	}
}

/*-------------------------------------------------------------------------*/
static int 
compare_planes(const void *x, const void *y)
{
	plane_t *xx = (plane_t *)x;
	plane_t *yy = (plane_t *)y;
	return (int)yy->plane.score - (int)xx->plane.score;
}

void
sieve_xy_run_deg6(root_sieve_t *rs)
{
	uint32 i;

	sieve_xyz_t *xyz = &rs->xyzdata;
	int64 z_base = xyz->z_base;

	sieve_xy_t *xy = &rs->xydata;
	sieve_prime_t *lattice_primes = xy->lattice_primes;
	uint32 num_lattice_primes;
	msieve_obj *obj = rs->data->obj;

	double direction[3] = {0, 1, 0};
	double line_min, line_max;
	uint16 cutoff_score;
	xydata_t xydata[MAX_CRT_FACTORS];
	plane_heap_t plane_heap;

	uint64 inv_xy;
	uint64 inv_xyz;

	compute_line_size(rs->max_norm, &rs->apoly,
			rs->dbl_p, rs->dbl_d, direction,
			-10000, 10000, &line_min, &line_max);
	if (line_min > line_max)
		return;

	num_lattice_primes = xy->num_lattice_primes = 
				find_lattice_primes(rs->primes, 
					rs->num_primes, xyz->lattice_size, 
					lattice_primes, &xy->lattice_size,
					line_max - line_min);

	inv_xy = mp_modinv_2(xyz->lattice_size, xy->lattice_size);
	inv_xyz = mp_modinv_2(xy->lattice_size, xyz->lattice_size);
	uint64_2gmp(xy->lattice_size, xy->tmp1);
	uint64_2gmp(inv_xy, xy->tmp2);
	uint64_2gmp(xyz->lattice_size, xy->tmp3);
	uint64_2gmp(inv_xyz, xy->tmp4);
	mpz_mul(xy->mp_lattice_size, xy->tmp1, xy->tmp3);
	mpz_mul(xy->crt0, xy->tmp2, xy->tmp3);
	mpz_mul(xy->crt1, xy->tmp1, xy->tmp4);
	xy->dbl_lattice_size = mpz_get_d(xy->mp_lattice_size);

	xydata_alloc(lattice_primes, num_lattice_primes, 
			xyz->lattice_size, xydata);

	plane_heap.num_entries = 0;

	for (i = 0; i < xyz->num_lattices; i++) {

		lattice_t *curr_lattice_xyz = xyz->lattices + i;

		xydata_init(xydata, num_lattice_primes,
				curr_lattice_xyz, z_base);

		find_hits(rs, xydata, num_lattice_primes, 
				i, &plane_heap);
	}

	xydata_free(xydata, num_lattice_primes);

	qsort(plane_heap.entries, plane_heap.num_entries,
			sizeof(plane_t), compare_planes);
	cutoff_score = 0.9 * plane_heap.entries[0].plane.score;

	for (i = 0; i < plane_heap.num_entries; i++) {

		plane_t *curr_plane = plane_heap.entries + i;
		lattice_t *lattice_xy = &curr_plane->plane;
		lattice_t *lattice_xyz = xyz->lattices + 
					curr_plane->which_lattice_xyz;

		if (lattice_xy->score < cutoff_score)
			break;

		line_min = xyz->y_line_min[curr_plane->which_z_block];
		line_max = xyz->y_line_max[curr_plane->which_z_block];
		z_base = xyz->z_base + curr_plane->which_z_block *
					xyz->lattice_size;

		xy->apoly = rs->apoly;
		xy->apoly.coeff[3] += z_base * rs->dbl_p;
		xy->apoly.coeff[2] -= z_base * rs->dbl_d;

		mpz_set_d(xy->tmp1, line_min);
		mpz_tdiv_q(xy->y_base, xy->tmp1, xy->mp_lattice_size);
		mpz_mul(xy->y_base, xy->y_base, xy->mp_lattice_size);
		xy->y_blocks = (line_max - line_min) / xy->dbl_lattice_size;

		uint64_2gmp(lattice_xy->x, xy->tmp1);
		uint64_2gmp(lattice_xyz->x, xy->tmp2);
		mpz_mul(xy->resclass_x, xy->tmp1, xy->crt0);
		mpz_addmul(xy->resclass_x, xy->tmp2, xy->crt1);

		uint64_2gmp(lattice_xy->y, xy->tmp1);
		uint64_2gmp(lattice_xyz->y, xy->tmp2);
		mpz_mul(xy->resclass_y, xy->tmp1, xy->crt0);
		mpz_addmul(xy->resclass_y, xy->tmp2, xy->crt1);

		mpz_tdiv_r(xy->resclass_x, xy->resclass_x, 
				xy->mp_lattice_size);
		mpz_tdiv_r(xy->resclass_y, xy->resclass_y, 
				xy->mp_lattice_size);

		xy->curr_score = lattice_xyz->score + lattice_xy->score;
		rs->curr_z = z_base + lattice_xyz->z; 

		sieve_x_run_deg6(rs);

		if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)
			break;
	}
}
