/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: tinyqs.c 590 2011-06-29 03:16:04Z jaysonking $
--------------------------------------------------------------------*/

#include <ms_common.h>

/* Basic MPQS implementation, intended for tiny inputs. Note
   that this code is very old and inefficient; it's only here
   to factor numbers that are small enough to be 'below the radar' 
   of the main SIQS routines. The practical upper limit on the
   input size is 85 bits, about 25-26 digits. */

#define LOGPRIME_SCALE_TINY 2

static const uint8 logprime_list[] = {
2,  3,  5,  6,  7,  7,  8,  8,  9,  10,
10, 10, 11, 11, 11, 11, 12, 12, 12, 12,
12, 13, 13, 13, 13, 13, 13, 13, 14, 14,
14, 14, 14, 14, 14, 14, 15, 15, 15, 15,
15, 15, 15, 15, 15, 15, 15, 16, 16, 16,
16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
16, 16, 17, 17, 17, 17, 17, 17, 17, 17,
17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
17, 17, 18, 18, 18, 18, 18, 18, 18, 18,
18, 18, 18, 18, 18, 18, 18, 18, 18, 18,
18, 18, 18, 18, 18, 18, 18, 18, 18, 18,
18, 19, 19, 19, 19, 19, 19, 19, 19, 19,
19, 19, 19, 19, 19, 19, 19, 19, 19, 19,
19, 19, 19, 19, 19, 19, 19, 19, 19, 19,
19, 19, 19, 19, 19, 19, 19, 19, 19, 20,
20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
20, 20, 20, 20, 20, 20, 20, 20, 20, 21,
};

#define BIT(x) ((uint64)(1) << (x))

static const uint64 bitmask[64] = {
	BIT( 0), BIT( 1), BIT( 2), BIT( 3), BIT( 4), BIT( 5), BIT( 6), BIT( 7),
	BIT( 8), BIT( 9), BIT(10), BIT(11), BIT(12), BIT(13), BIT(14), BIT(15),
	BIT(16), BIT(17), BIT(18), BIT(19), BIT(20), BIT(21), BIT(22), BIT(23),
	BIT(24), BIT(25), BIT(26), BIT(27), BIT(28), BIT(29), BIT(30), BIT(31),
	BIT(32), BIT(33), BIT(34), BIT(35), BIT(36), BIT(37), BIT(38), BIT(39),
	BIT(40), BIT(41), BIT(42), BIT(43), BIT(44), BIT(45), BIT(46), BIT(47),
	BIT(48), BIT(49), BIT(50), BIT(51), BIT(52), BIT(53), BIT(54), BIT(55),
	BIT(56), BIT(57), BIT(58), BIT(59), BIT(60), BIT(61), BIT(62), BIT(63),
};

#define MIN_FB_OFFSET 1
#define NUM_PRIMES_TINY (sizeof(logprime_list) / sizeof(uint8))
#define MAX_FB_SIZE_TINY 85
#define NUM_EXTRA_RELATIONS_TINY 16
#define MAX_RELATIONS_TINY (MAX_FB_SIZE_TINY + NUM_EXTRA_RELATIONS_TINY)
#define MIN_FB_OFFSET_TO_SIEVE_TINY 7
#define SMALL_PRIME_FUDGE_TINY (LOGPRIME_SCALE_TINY * 10)
#define MAX_POLY_TINY 32
#define SIEVE_BLOCK_SIZE_TINY 32768
#define MAX_FACTORS_TINY 20
#define LOG2_PARTIAL_TABLE_SIZE 10
#define LARGE_PRIME_HASH(x) (((uint32)(x) * ((uint32)40499 * 65543)) >> \
				(32 - LOG2_PARTIAL_TABLE_SIZE))

typedef struct {
	uint32 sieve_offset;
	uint16 large_prime;
	uint8 poly_num;
	uint8 num_factors;
	uint8 fb_offsets[MAX_FACTORS_TINY];
} tiny_relation;

typedef struct {
	uint16 prime;
	uint8 logprime;
	uint16 modsqrt;
	uint16 root1;
	uint16 root2;
	uint32 next_loc1;
	uint32 next_loc2;
} tiny_fb;
	
typedef struct {
	mp_t n;
	mp_t kn;
	uint32 multiplier;

	mp_t target_a;
	uint8 poly_num;
	uint32 curr_poly_offset;
	uint32 poly_offset[MAX_POLY_TINY];
	mp_t poly_b[MAX_POLY_TINY];

	double align_me;
	uint8 sieve_block[SIEVE_BLOCK_SIZE_TINY];
	uint8 num_sieve_blocks;

	uint32 fb_size;
	tiny_fb factor_base[MAX_FB_SIZE_TINY];

	uint32 num_relations;
	uint32 target_relations;
	uint16 large_prime_max;
	uint32 error_bits;
	tiny_relation relation_list[MAX_RELATIONS_TINY];
	tiny_relation partial_list[1 << LOG2_PARTIAL_TABLE_SIZE];

	uint64 matrix[MAX_FB_SIZE_TINY][(MAX_RELATIONS_TINY+63) / 64];
	uint16 null_vectors[MAX_RELATIONS_TINY];

	uint32 seed1;
	uint32 seed2;
} tiny_qs_params;

/*----------------------------------------------------------------------*/
static uint32 init_one_fb_tiny(uint32 fb_size, 
			mp_t *n, tiny_fb *factor_base) {

	uint32 i, j;
	uint32 prime;

	i = MIN_FB_OFFSET;
	prime = 2;
	factor_base[i].prime = prime;

	for (i++, j = 1; i < fb_size && j < NUM_PRIMES_TINY; j++) {

		tiny_fb *fbptr = factor_base + i;
		uint32 nmodp;

		prime += prime_delta[j];
		nmodp = mp_mod_1(n, prime);
		if (nmodp == 0) {
			fbptr->prime = (uint16)prime;
			fbptr->logprime = logprime_list[j];
			i++;
		}
		else if (mp_legendre_1(nmodp, prime) == 1) {
			fbptr->prime = (uint16)prime;
			fbptr->logprime = logprime_list[j];
			fbptr->modsqrt = mp_modsqrt_1(nmodp, prime);
			i++;
		}
	}

	return i;
}

/*----------------------------------------------------------------------*/
static void init_fb_tiny(tiny_qs_params *params) {

	uint32 i, j;
	tiny_fb next_fb_array[MAX_FB_SIZE_TINY];
	uint32 next_fb_size;
	double score, best_score;
	mp_t kn;
	uint16 mult_list[] = {1, 3, 5, 7, 11, 13, 15, 17, 19, 21, 23,
			 29, 31, 33, 35, 37, 39, 41, 43, 47, 51,
			 53, 55, 57, 59, 61, 65, 67, 69, 71, 73};

	best_score = 1000.0;
	score = 0;

	for (i = 0; i < sizeof(mult_list) / sizeof(uint16); i++) {

		uint32 curr_mult = mult_list[i];

		mp_mul_1(&params->n, curr_mult, &kn);

		next_fb_size = init_one_fb_tiny(MAX_FB_SIZE_TINY,
						&kn, next_fb_array);
		
		if ((params->n.val[0] & 7) == 1)
			score = 0.5 * log((double)curr_mult) - 2 * M_LN2;
		else if ((params->n.val[0] & 7) == 5)
			score = 0.5 * log((double)curr_mult) - M_LN2;
		else
			score = 0.5 * log((double)curr_mult) - 0.5 * M_LN2;

		for (j = MIN_FB_OFFSET + 1; j < next_fb_size; j++) {
			tiny_fb *fbptr = next_fb_array + j;
			uint32 prime = fbptr->prime;

			if (prime <= 73 && curr_mult % prime == 0)
				score -= log((double)prime) / (prime - 1);
			else
				score -= 2.0 * log((double)prime) / (prime - 1);
		}

		if (score < best_score) {
			memcpy(params->factor_base, next_fb_array,
					sizeof(params->factor_base));
			best_score = score;
			params->fb_size = next_fb_size;
			params->multiplier = curr_mult;
			mp_copy(&kn, &params->kn);
		}
	}
}

/*----------------------------------------------------------------------*/
static void fill_sieve_block_tiny(tiny_qs_params *params) {

	uint32 i;
	uint32 fb_size = params->fb_size;
	uint8 *sieve_block = params->sieve_block;
	tiny_fb *factor_base = params->factor_base;

	for (i = MIN_FB_OFFSET_TO_SIEVE_TINY; i < fb_size; i++) {
		tiny_fb *fbptr = factor_base + i;
		uint32 prime = fbptr->prime;
		uint8 logprime = fbptr->logprime;
		uint32 root1 = fbptr->next_loc1;
		uint32 root2 = fbptr->next_loc2;

		if (fbptr->root1 == (uint16)(-1))
			continue;

		while (root2 < SIEVE_BLOCK_SIZE_TINY) {
			sieve_block[root1] -= logprime;
			sieve_block[root2] -= logprime;
			root1 += prime;
			root2 += prime;
		}

		if (root1 < SIEVE_BLOCK_SIZE_TINY) {
			sieve_block[root1] -= logprime;
			root1 += prime;
			fbptr->next_loc1 = root2 - SIEVE_BLOCK_SIZE_TINY;
			fbptr->next_loc2 = root1 - SIEVE_BLOCK_SIZE_TINY;
		}
		else {
			fbptr->next_loc1 = root1 - SIEVE_BLOCK_SIZE_TINY;
			fbptr->next_loc2 = root2 - SIEVE_BLOCK_SIZE_TINY;
		}
	}
}

/*----------------------------------------------------------------------*/
static void check_sieve_val_tiny(tiny_qs_params *params, 
				mp_t *a, mp_t *b, mp_t *c, 
				uint32 sieve_offset,
				uint32 sign_of_index, 
				uint32 bits) {

	uint32 i, j;
	uint32 fb_size = params->fb_size;
	tiny_fb *factor_base = params->factor_base;
	uint32 num_factors = 0;
	mp_t res;
	uint32 cutoff2;
	tiny_relation *relation = params->relation_list +
					params->num_relations;
	uint8 *fb_offsets = relation->fb_offsets;

	mp_mul_1(a, sieve_offset, &res);
	if (sign_of_index == POSITIVE)
		mp_add(&res, b, &res);
	else
		mp_sub(&res, b, &res);

	mp_mul_1(&res, sieve_offset, &res);
	if (mp_cmp(&res, c) >= 0) {
		mp_sub(&res, c, &res);
	}
	else {
		mp_sub(c, &res, &res);
		fb_offsets[num_factors++] = 0;
	}

	cutoff2 = LOGPRIME_SCALE_TINY * mp_bits(&res) - params->error_bits;

	i = mp_rjustify(&res, &res);
	for (bits += i * LOGPRIME_SCALE_TINY; i; i--)
		fb_offsets[num_factors++] = MIN_FB_OFFSET;

	for (i = MIN_FB_OFFSET + 1; i < MIN_FB_OFFSET_TO_SIEVE_TINY; i++) {
		tiny_fb *fbptr = factor_base + i;
		uint16 prime = fbptr->prime;
		uint8 logprime = fbptr->logprime;
		uint16 root1 = fbptr->root1;
		uint16 root2 = fbptr->root2;

		j = sieve_offset % prime;

		if (root1 == (uint16)(-1)) {
			if (mp_mod_1(&res, prime) == 0) {
				do {
					if (num_factors >= MAX_FACTORS_TINY)
						return;
					bits += logprime;
					fb_offsets[num_factors++] = i;
					mp_divrem_1(&res, prime, &res);
					j = mp_mod_1(&res, prime);
				} while (j == 0);
			}
			continue;
		}

		if (sign_of_index == NEGATIVE) {
			root2 = prime - root2;
			if (root1)
				root1 = prime - root1;
		}

		if (j == root1 || j == root2) {
			do {
				if (num_factors >= MAX_FACTORS_TINY)
					return;
				bits += logprime;
				fb_offsets[num_factors++] = i;
				mp_divrem_1(&res, prime, &res);
				j = mp_mod_1(&res, prime);
			} while (j == 0);
		}
	}

	if (bits < cutoff2)
		return;

	for (; i < fb_size; i++) {
		tiny_fb *fbptr = factor_base + i;
		uint16 prime = fbptr->prime;
		uint16 root1 = fbptr->root1;
		uint16 root2 = fbptr->root2;

		j = sieve_offset % prime;

		if (root1 == (uint16)(-1)) {
			if (mp_mod_1(&res, prime) == 0) {
				do {
					if (num_factors >= MAX_FACTORS_TINY)
						return;
					fb_offsets[num_factors++] = i;
					mp_divrem_1(&res, prime, &res);
					j = mp_mod_1(&res, prime);
				} while (j == 0);
			}
			continue;
		}

		if (sign_of_index == NEGATIVE) {
			root2 = prime - root2;
			if (root1)
				root1 = prime - root1;
		}

		if (j == root1 || j == root2) {
			do {
				if (num_factors >= MAX_FACTORS_TINY)
					return;
				fb_offsets[num_factors++] = i;
				mp_divrem_1(&res, (uint32)prime, &res);
				j = mp_mod_1(&res, prime);
			} while (j == 0);
		}
	}

	if (res.nwords > 1)
		return;

	if (sign_of_index == NEGATIVE)
		sieve_offset |= 0x80000000;

	relation->sieve_offset = sieve_offset;
	relation->num_factors = num_factors;
	relation->poly_num = params->poly_num - 1;

	if (res.val[0] == 1) {
		relation->large_prime = 1;
		params->num_relations++;
		return;
	}

	if (res.val[0] < params->large_prime_max) {
		uint32 table_idx = LARGE_PRIME_HASH(res.val[0]);
		uint16 hash_entry = params->partial_list[table_idx].large_prime;

		relation->large_prime = res.val[0];
		if (hash_entry == res.val[0])
			params->num_relations++;
		else if (hash_entry == 0)
			memcpy(params->partial_list + table_idx, relation,
					sizeof(tiny_relation));
	}
}

/*----------------------------------------------------------------------*/
#define PACKED_SIEVE_MASK ((uint64)0x80808080 << 32 | 0x80808080)

static void sieve_next_poly_tiny(tiny_qs_params *params) {

	uint32 i, j, k;
	uint32 fb_size = params->fb_size;
	uint8 *sieve_block = params->sieve_block;
	uint64 *packed_sieve_block = (uint64 *)params->sieve_block;
	uint32 block_start;
	tiny_fb *factor_base = params->factor_base;
	mp_t a, b, c, tmp;
	mp_t *kn = &params->kn;
	uint32 cutoff1;
	uint8 poly_num = params->poly_num;
	uint8 num_sieve_blocks = params->num_sieve_blocks;
	uint32 target_relations = params->target_relations;

	i = params->curr_poly_offset;
	do {
		mp_add_1(&params->target_a, i, &a);
		i += mp_next_prime(&a, &tmp, &params->seed1,
					&params->seed2);
	} while(mp_legendre(kn, &tmp) != 1);

	mp_modsqrt2(kn, &tmp, &b, &params->seed1, &params->seed2);
	mp_mul(&tmp, &tmp, &a);

	params->curr_poly_offset = i;
	params->poly_offset[poly_num] = i;
	mp_copy(&b, &params->poly_b[poly_num]);
	params->poly_num++;

	mp_mul(&b, &b, &tmp);
	mp_sub(kn, &tmp, &tmp);
	mp_div(&tmp, &a, &c);

	for (i = MIN_FB_OFFSET + 1; i < fb_size; i++) {
		tiny_fb *fbptr = factor_base + i;
		uint32 prime = fbptr->prime;
		uint32 modsqrt = fbptr->modsqrt;
		uint32 amodp = mp_mod_1(&a, prime);
		uint32 bmodp = prime - mp_mod_1(&b, prime);
		uint32 root1, root2;

		if (amodp == 0 || 
		    (prime < 64 && params->multiplier % prime == 0)) {
			fbptr->root1 = (uint16)(-1);
			fbptr->root2 = (uint16)(-1);
			continue;
		}

		amodp = mp_modinv_1(amodp, prime);
		root1 = mp_modmul_1(amodp, bmodp + modsqrt, prime);
		root2 = mp_modmul_1(amodp, bmodp + prime - modsqrt, prime);

		if (root1 < root2) {
			fbptr->root1 = root1;
			fbptr->root2 = root2;
			fbptr->next_loc1 = root1;
			fbptr->next_loc2 = root2;
		}
		else {
			fbptr->root1 = root2;
			fbptr->root2 = root1;
			fbptr->next_loc1 = root2;
			fbptr->next_loc2 = root1;
		}
	}

	cutoff1 = LOGPRIME_SCALE_TINY * mp_bits(&c) - 
			params->error_bits - SMALL_PRIME_FUDGE_TINY;
	mp_add(&b, &b, &b);
	
	for (i = block_start = 0; i < num_sieve_blocks; i++) {

		memset(sieve_block, (int32)(cutoff1 - 1), 
					(size_t)SIEVE_BLOCK_SIZE_TINY);
		fill_sieve_block_tiny(params);

		for (j = 0; j < SIEVE_BLOCK_SIZE_TINY / 8; j += 4) {
			uint64 accum = packed_sieve_block[j] |
					packed_sieve_block[j+1] |
					packed_sieve_block[j+2] |
					packed_sieve_block[j+3];

			if ((accum & PACKED_SIEVE_MASK) == (uint64)(0))
				continue;

			for (k = 0; k < 32; k++) {
				uint32 bits = sieve_block[8 * j + k];
				if (bits <= cutoff1)
					continue;
				check_sieve_val_tiny(params, &a, &b, &c,
					block_start + 8 * j + k, 
					POSITIVE, cutoff1 + 257 - bits);
				if (params->num_relations == target_relations)
					return;
			}
		}
		block_start += SIEVE_BLOCK_SIZE_TINY;
	}

	for (i = MIN_FB_OFFSET + 1; i < fb_size; i++) {
		tiny_fb *fbptr = factor_base + i;
		uint32 prime = fbptr->prime;

		fbptr->next_loc1 = prime - fbptr->root2;
		fbptr->next_loc2 = prime - fbptr->root1;
	}
		
	for (i = block_start = 0; i < num_sieve_blocks; i++) {

		memset(sieve_block, (int32)(cutoff1 - 1), 
					(size_t)SIEVE_BLOCK_SIZE_TINY);
		fill_sieve_block_tiny(params);

		for (j = 0; j < SIEVE_BLOCK_SIZE_TINY / 8; j += 4) {
			uint64 accum = packed_sieve_block[j] |
					packed_sieve_block[j+1] |
					packed_sieve_block[j+2] |
					packed_sieve_block[j+3];

			if ((accum & PACKED_SIEVE_MASK) == (uint64)(0))
				continue;

			for (k = 0; k < 32; k++) {
				uint32 bits = sieve_block[8 * j + k];
				if (bits <= cutoff1)
					continue;
				check_sieve_val_tiny(params, &a, &b, &c,
					block_start + 8 * j + k, 
					NEGATIVE, cutoff1 + 257 - bits);
				if (params->num_relations == target_relations)
					return;
			}
		}
		block_start += SIEVE_BLOCK_SIZE_TINY;
	}
}

/*----------------------------------------------------------------------*/
static void solve_linear_system_tiny(tiny_qs_params *params) {

	uint32 i, j, k, start_row;
	uint32 nrows = params->fb_size;
	uint32 ncols = params->num_relations;
	uint8 rowperm[MAX_FB_SIZE_TINY];
	uint8 pivot[MAX_FB_SIZE_TINY];
	uint8 row = 0;

	memset(params->matrix, 0, sizeof(params->matrix));

	for (i = 0; i < ncols; i++) {
		tiny_relation *r = params->relation_list + i;
		for (j = 0; j < r->num_factors; j++) {
			row = r->fb_offsets[j];
			params->matrix[row][i / 64] ^= bitmask[i & 63];
		}
		if (r->large_prime > 1) {
			r = params->partial_list + 
				LARGE_PRIME_HASH(r->large_prime);
			for (j = 0; j < r->num_factors; j++) {
				row = r->fb_offsets[j];
				params->matrix[row][i / 64] ^= bitmask[i & 63];
			}
		}
	}
	for (i = 0; i < nrows; i++)
		rowperm[i] = i;
	for (i = 0; i < ncols; i++)
		params->null_vectors[i] = (uint16)get_rand(
						&params->seed1,
						&params->seed2);

	for (i = start_row = 0; start_row < nrows && i < ncols; i++) {
		
		for (j = start_row; j < nrows; j++) {
			row = rowperm[j];
			if (params->matrix[row][i / 64] & bitmask[i & 63])
				break;
		}
		if (j == nrows)
			continue;

		rowperm[j] = rowperm[start_row];
		rowperm[start_row] = row;
		pivot[start_row++] = i;

		for (j++; j < nrows; j++) {
			uint8 row2 = rowperm[j];
			if (params->matrix[row2][i / 64] & bitmask[i & 63]) {
				for (k = i / 64; k < (ncols + 63) / 64; k++) {
					params->matrix[row2][k] ^=
					params->matrix[row][k];
				}
			}
		}
	}

	for (i = start_row - 1; (int32)i >= 0; i--) {
		uint16 accum;
		row = rowperm[i];

		for (j = pivot[i] + 1, accum = 0; j < ncols; j++) {
			if (params->matrix[row][j / 64] & bitmask[j & 63])
				accum ^= params->null_vectors[j];
		}
		params->null_vectors[pivot[i]] = accum;
	}
}

/*----------------------------------------------------------------------*/
static uint32 find_factors_tiny(tiny_qs_params *params,
				mp_t *factor1, mp_t *factor2) {

	mp_t x, y, t0, t1;
	mp_t *n = &params->n;
	mp_t *kn = &params->kn;
	uint16 i, j, k;
	uint16 mask;
	uint16 fb_counts[MAX_FB_SIZE_TINY];

	for (mask = 1; mask; mask <<= 1) {

		memset(fb_counts, 0, sizeof(fb_counts));
		mp_clear(&x);
		mp_clear(&y);
		x.nwords = x.val[0] = 1;
		y.nwords = y.val[0] = 1;

		for (i = 0; i < params->num_relations; i++) {

			if (!(params->null_vectors[i] & mask))
				continue;

			for (j = 0; j < 2; j++) {
				tiny_relation *r = params->relation_list + i;
				uint32 poly_num;

				if (j == 1) {
					if (r->large_prime > 1) {
						r = params->partial_list +
					    		LARGE_PRIME_HASH(
								r->large_prime);
						mp_mul_1(&y, 
							r->large_prime, &t0);
						mp_mod(&t0, kn, &y);
					}
					else {
						break;
					}
				}
				poly_num = r->poly_num;
	
				mp_add_1(&params->target_a, 
					params->poly_offset[poly_num], &t0);
				mp_modmul(&y, &t0, kn, &y);
	
				for (k = 0; k < r->num_factors; k++)
					fb_counts[r->fb_offsets[k]]++;
	
				mp_mul(&t0, &t0, &t1);
				mp_mul_1(&t1, r->sieve_offset & 
						0x7fffffff, &t1);
				if (r->sieve_offset >> 31)
					mp_sub(&t1, &params->poly_b[poly_num],
								&t1);
				else
					mp_add(&t1, &params->poly_b[poly_num], 
								&t1);
				mp_modmul(&x, &t1, kn, &x);
			}
		}

		for (i = MIN_FB_OFFSET; i < params->fb_size; i++) {
			uint16 mask2 = 0x8000;
			uint16 exponent = fb_counts[i] / 2;
			uint32 prime = params->factor_base[i].prime;

			if (exponent == 0)
				continue;

			mp_clear(&t0);
			t0.nwords = 1;
			t0.val[0] = prime;
			mp_copy(&t0, &t1);

			while (!(exponent & mask2))
				mask2 >>= 1;

			for (mask2 >>= 1; mask2; mask2 >>= 1) {
				mp_modmul(&t0, &t0, kn, &t0);
				if (exponent & mask2) {
					mp_modmul(&t0, &t1, kn, &t0);
				}
			}
			mp_modmul(&t0, &y, kn, &y);
		}

		for (i = 0; i < 2; i++) {
			if (i == 0)
				mp_add(&x, &y, &t0);
			else if (mp_cmp(&x, &y) > 0)
				mp_sub(&x, &y, &t0);
			else
				mp_sub(&y, &x, &t0);

			mp_gcd(&t0, kn, &t1);
			if (!mp_is_one(&t1) && mp_cmp(&t1, kn)) {
				if (params->multiplier > 1) {
					mp_t m;
					mp_clear(&m);
					m.nwords = 1; 
					m.val[0] = params->multiplier;
					mp_gcd(&t1, &m, &t0);
					mp_divrem_1(&t1, t0.val[0], &t1);
				}
				if (!mp_is_one(&t1) && mp_cmp(&t1, n)) {
					mp_div(n, &t1, factor1);
					mp_copy(&t1, factor2);
					return 1;
				}
			}
		}
	}

	return 0;
}

/*----------------------------------------------------------------------*/
uint32 tinyqs(mp_t *n, mp_t *factor1, mp_t *factor2) {

	tiny_qs_params *params;
	uint32 i, bits;
	uint32 fb_size, status = 0;
	uint16 bound;
	uint16 large_prime_mult;
	mp_t tmp1, tmp2;

	params = (tiny_qs_params *)xmalloc(sizeof(tiny_qs_params));

	mp_copy(n, &params->n);
	params->num_relations = 0;
	params->seed1 = 287643287;
	params->seed2 = 833277363;
	params->curr_poly_offset = 0;
	params->poly_num = 0;

	init_fb_tiny(params);
	bits = mp_bits(&params->kn);

	if (bits < 70) {
		fb_size = MIN(params->fb_size, 40);
		params->num_sieve_blocks = 4;
		large_prime_mult = 1;
	}
	else if (bits < 75) {
		fb_size = MIN(params->fb_size, 50);
		params->num_sieve_blocks = 4;
		large_prime_mult = 5;
	}
	else if (bits < 80) {
		fb_size = MIN(params->fb_size, 60);
		params->num_sieve_blocks = 4;
		large_prime_mult = 10;
	}
	else if (bits < 85) {
		fb_size = MIN(params->fb_size, 70);
		params->num_sieve_blocks = 8;
		large_prime_mult = 10;
	}
	else {
		fb_size = MIN(params->fb_size, 85);
		params->num_sieve_blocks = 8;
		large_prime_mult = 15;
	}
	params->fb_size = fb_size;

	mp_add(&params->n, &params->n, &tmp1);
	mp_isqrt(&tmp1, &tmp2);
	mp_divrem_1(&tmp2, (uint32)(params->num_sieve_blocks * 
			SIEVE_BLOCK_SIZE_TINY), &tmp2);
	mp_isqrt(&tmp2, &params->target_a);

	bound = params->factor_base[params->fb_size - 1].prime;
	bound *= large_prime_mult;
	params->large_prime_max = bound;
	params->error_bits = (uint32)(LOGPRIME_SCALE_TINY * 
					(log((double)bound) / M_LN2 + 1));
	for (i = 0; i < (1 << LOG2_PARTIAL_TABLE_SIZE); i++)
		params->partial_list[i].large_prime = 0;

	params->target_relations = params->fb_size + NUM_EXTRA_RELATIONS_TINY;
	while (params->poly_num < MAX_POLY_TINY) {
		sieve_next_poly_tiny(params);

		if (params->num_relations == params->target_relations) {

			solve_linear_system_tiny(params);
			status = find_factors_tiny(params, factor1, factor2);

			if (!status) {
				/* Failed to find a nontrivial solution.
				   Throw away some relations and try again. */

				uint32 discard = params->num_relations / 10 + 1;

				params->num_relations -= discard;
			}
			else {
				break;
			}
		}
	}

	free(params);
	return status;
}
