/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: root_sieve_deg5_xy.c 817 2012-11-11 14:58:29Z jasonp_sf $
--------------------------------------------------------------------*/

#include "stage2.h"

/*-------------------------------------------------------------------------*/
uint64
find_lattice_size_y(double line_length)
{
	if (line_length < 1e3)
		return 1;
	else if (line_length < 1e4)
		return (uint64)2*2*3*5;
	else if (line_length < 1e5)
		return (uint64)2*2*2*3*3*5;
	return (uint64)2*2*2*3*3*5*7;
}

/*-------------------------------------------------------------------------*/
static uint32
find_lattice_primes(sieve_prime_t *primes, uint32 num_primes,
			uint64 lattice_size, sieve_prime_t *lattice_primes)
{
	uint32 i, j;
	uint32 num_lattice_primes = 0;

	for (i = 0; i < num_primes; i++) {

		sieve_prime_t *curr_prime = primes + i;
		uint32 num_powers = curr_prime->num_powers;
		uint32 num_lattice_powers = 0;

		if (lattice_size % curr_prime->prime)
			break;

		for (j = 0; j < num_powers; j++, num_lattice_powers++) {
			sieve_power_t *sp = curr_prime->powers + j;
			if (lattice_size % sp->power)
				break;
		}

		lattice_primes[num_lattice_primes] = *curr_prime;
		lattice_primes[num_lattice_primes].num_powers = 
						num_lattice_powers;
		num_lattice_primes++;
	}

	return num_lattice_primes;
}

/*-------------------------------------------------------------------------*/
static void 
do_sieving(sieve_root_t *r, uint16 *sieve,
		uint32 contrib, uint32 dim)
{
	uint32 i;
	uint32 start = r->start;
	uint32 step = r->step;
	uint32 resclass = r->resclass;

	if (resclass >= step)
		resclass %= step;

	for (i = 0; i < dim; i++) {

		uint32 ri = start;

		do {
			sieve[ri] += contrib;
			ri += step;
		} while (ri < dim);

		sieve += dim;
		start = mp_modsub_1(start, resclass, step);
	}
}

/*-------------------------------------------------------------------------*/
#define MAX_DIM 16

static void 
find_hits(sieve_prime_t *lattice_primes,
	uint32 num_lattice_primes, hit_t *hitlist)
{
	uint32 i, j, k;
	uint16 sieve[MAX_DIM * MAX_DIM];

	for (i = 0; i < num_lattice_primes; i++) {

		sieve_prime_t *curr_prime = lattice_primes + i;
		uint32 num_powers = curr_prime->num_powers;
		uint32 dim = curr_prime->powers[num_powers-1].power;
		uint32 size = dim * dim;
		uint32 max_sieve_val = 0;
		hit_t *hits = hitlist + i;

		memset(sieve, 0, size * sizeof(uint16));

		for (j = 0; j < num_powers; j++) {

			sieve_power_t *sp = curr_prime->powers + j;
			uint32 num_roots = sp->num_roots;
			uint32 contrib = sp->sieve_contrib;

			for (k = 0; k < num_roots; k++) {
				do_sieving(sp->roots + k, sieve, 
						contrib, dim);
			}
		}

		for (j = k = 0; j < size; j++) {
			uint32 curr_score = sieve[j];

			if (curr_score >= max_sieve_val) {
				if (curr_score > max_sieve_val) {
					max_sieve_val = curr_score;
					k =0;
				}

				if (k == MAX_ROOTS)
					break;

				hits->score[k] = curr_score;
				hits->roots[k][0] = j % dim;
				hits->roots[k][1] = (j / dim) % dim;
				hits->roots[k][2] = 0;
				k++;
			}
		}
		hits->power = dim;
		hits->num_roots = k;
	}
}

/*-------------------------------------------------------------------------*/
void
sieve_xy_alloc(sieve_xy_t *xy)
{
	memset(xy, 0, sizeof(sieve_xy_t));
	mpz_init(xy->y_base);
	mpz_init(xy->mp_lattice_size);
	mpz_init(xy->resclass_x);
	mpz_init(xy->resclass_y);
	mpz_init(xy->crt0);
	mpz_init(xy->crt1);
	mpz_init(xy->tmp1);
	mpz_init(xy->tmp2);
	mpz_init(xy->tmp3);
	mpz_init(xy->tmp4);
}

/*-------------------------------------------------------------------------*/
void
sieve_xy_free(sieve_xy_t *xy)
{
	mpz_clear(xy->y_base);
	mpz_clear(xy->mp_lattice_size);
	mpz_clear(xy->resclass_x);
	mpz_clear(xy->resclass_y);
	mpz_clear(xy->crt0);
	mpz_clear(xy->crt1);
	mpz_clear(xy->tmp1);
	mpz_clear(xy->tmp2);
	mpz_clear(xy->tmp3);
	mpz_clear(xy->tmp4);
	free(xy->x_line_min);
	free(xy->x_line_max);
}

/*-------------------------------------------------------------------------*/
void
sieve_xy_run_deg5(root_sieve_t *rs, uint64 lattice_size,
		double line_min, double line_max)
{
	uint32 i, j;
	sieve_xy_t *xy = &rs->xydata;
	hit_t hitlist[MAX_CRT_FACTORS];
	uint32 num_lattice_primes;
	uint32 num_lattices;
	uint32 y_blocks;
	int64 curr_y;
	double direction[3] = {0, 1, 0};

	xy->lattice_size = lattice_size;
	xy->dbl_lattice_size = (double)lattice_size;
	uint64_2gmp(lattice_size, xy->mp_lattice_size);

	mpz_set_d(xy->y_base, line_min / lattice_size - 1);
	mpz_mul(xy->y_base, xy->y_base, xy->mp_lattice_size);

	y_blocks = (line_max - line_min) / lattice_size + 1;
	if (y_blocks > xy->y_blocks) {
		xy->x_line_min = (double *)xrealloc(xy->x_line_min,
						y_blocks * sizeof(double));
		xy->x_line_max = (double *)xrealloc(xy->x_line_max,
						y_blocks * sizeof(double));
	}
	xy->y_blocks = y_blocks;

	xy->num_lattices = 0;
	if (lattice_size == 1) {
		num_lattice_primes = xy->num_lattice_primes = 0;
		num_lattices = 1;
		if (num_lattices > xy->num_lattices) {
			xy->lattices = (lattice_t *)xrealloc(xy->lattices,
					num_lattices * sizeof(lattice_t));
		}
		memset(xy->lattices, 0, sizeof(lattice_t));
	}
	else {
		num_lattice_primes = xy->num_lattice_primes = 
				find_lattice_primes(rs->primes, 
						rs->num_primes, lattice_size, 
						xy->lattice_primes);

		find_hits(xy->lattice_primes, num_lattice_primes, hitlist);

		for (i = 0, num_lattices = 1; i < num_lattice_primes; i++) {
			num_lattices *= hitlist[i].num_roots;
		}

		if (num_lattices > xy->num_lattices) {
			xy->lattices = (lattice_t *)xrealloc(xy->lattices,
					num_lattices * sizeof(lattice_t));
		}

		compute_lattices(hitlist, num_lattice_primes, xy->lattices,
				lattice_size, num_lattices, 2);

	}
	xy->num_lattices = num_lattices;

	line_min = -10000;
	line_max = 10000;
	direction[0] = 1;
	direction[1] = 0;
	direction[2] = 0;
	curr_y = gmp2int64(xy->y_base);
	for (i = 0; i < y_blocks; i++) {
		dpoly_t apoly = rs->apoly;

		apoly.coeff[2] += rs->dbl_p * curr_y;
		apoly.coeff[1] -= rs->dbl_d * curr_y;

		compute_line_size(rs->max_norm, &apoly,
			rs->dbl_p, rs->dbl_d, direction,
			line_min, line_max, &line_min, &line_max);

		if (line_min >= line_max) {
			xy->x_line_min[i] = 0;
			xy->x_line_max[i] = 0;
			line_min = -10000;
			line_max = 10000;
		}
		else {
			xy->x_line_min[i] = line_min;
			xy->x_line_max[i] = line_max;
		}

		curr_y += lattice_size;
	}

	for (i = y_blocks; i; i--) {
		if (xy->x_line_min[i-1] != xy->x_line_max[i-1])
			break;
	}
	y_blocks = i;
	for (i = 0; i < y_blocks; i++) {
		if (xy->x_line_min[i] != xy->x_line_max[i])
			break;
	}
	mpz_addmul_ui(xy->y_base, xy->mp_lattice_size, i);
	y_blocks -= i;
	if (i > 0) {
		for (j = 0; j < y_blocks; j++) {
			xy->x_line_min[j] = xy->x_line_min[j+i];
			xy->x_line_max[j] = xy->x_line_max[j+i];
		}
	}
	xy->y_blocks = y_blocks;

#if 0
	printf("\n%.0lf %u %u\n", (double)lattice_size, 
			y_blocks, num_lattices);
#endif

	sieve_x_run_deg5(rs);
}
