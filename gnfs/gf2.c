/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: gf2.c 1025 2018-08-19 02:20:28Z jasonp_sf $
--------------------------------------------------------------------*/

#include <ms_common.h>
#include "gnfs.h"

/* the number of quadratic characters for each
   matrix column. A practical upper limit given
   in Buhler et. al. is 3 * log2(n), while 
   Bernstein writes that a QCB size of 50 'should 
   be enough for any possible NFS factorization'. 
   If each quadratic character reduces by half the
   odds that the linear algebra will not produce
   an algebraic square, then a very few characters
   (say 16) will be enough. The matrix-building code 
   is flexible enough so that this number may take 
   any positive value and everything else will 
   just work 
   
   The size of the QCB is limited to allow caching
   the quadratic characters */

#define QCB_SIZE 20

#if QCB_SIZE > 32
#error "QCB size must be 32 or less"
#endif

/*------------------------------------------------------------------*/
#define IDEAL_MINUS_ONE  ((uint64)0x7fff << 32 | (uint32)(-1))

static int compare_ideals(const void *x, const void *y) {
	
	/* used to determine the ordering of two ideals.
	   Ordering is by prime, then by root of prime,
	   then by rational or algebraic type. The ordering 
	   by prime is tricky, because -1 has a special value 
	   that must be explicitly accounted for. This ordering
	   is designed to put the most dense matrix rows first */

	ideal_t *k = (ideal_t *)x;
	ideal_t *t = (ideal_t *)y;
	uint64 p_k = (uint64)k->p_hi << 32 | k->p_lo;
	uint64 p_t = (uint64)t->p_hi << 32 | t->p_lo;
	uint64 r_k, r_t;

	if (p_k == IDEAL_MINUS_ONE || p_t == IDEAL_MINUS_ONE) {
		if (p_k == IDEAL_MINUS_ONE && p_t != IDEAL_MINUS_ONE)
			return -1;
		if (p_k != IDEAL_MINUS_ONE && p_t == IDEAL_MINUS_ONE)
			return 1;
		return 0;
	}

	if (p_k < p_t)
		return -1;
	if (p_k > p_t)
		return 1;
		
	r_k = (uint64)k->r_hi << 32 | k->r_lo;
	r_t = (uint64)t->r_hi << 32 | t->r_lo;

	if (r_k < r_t)
		return -1;
	if (r_k > r_t)
		return 1;

	if (k->rat_or_alg < t->rat_or_alg)
		return -1;
	if (k->rat_or_alg > t->rat_or_alg)
		return 1;
		
	return 0;
}

/*------------------------------------------------------------------*/
#define MAX_SMALL_IDEALS 200

static ideal_t *fill_small_ideals(factor_base_t *fb,
				uint32 *num_ideals_out,
				uint32 *max_small_ideal) {
	uint32 i, j;
	uint32 num_ideals;
	ideal_t *small_ideals;
	uint32 p;

	if (MAX_PACKED_PRIME == 0) {
		*num_ideals_out = 0;
		*max_small_ideal = 0;
		return NULL;
	}

	small_ideals = (ideal_t *)xmalloc(MAX_SMALL_IDEALS * sizeof(ideal_t));

	/* fill in the rational ideal of -1 */

	small_ideals[0].rat_or_alg = RATIONAL_IDEAL;
	small_ideals[0].p_lo = (uint32)(-1);
	small_ideals[0].p_hi = 0x7fff;
	small_ideals[0].r_lo = (uint32)(-1);
	small_ideals[0].r_hi = 0xffff;
	num_ideals = 1;

	/* for each small prime */

	for (i = p = 0; i < PRECOMPUTED_NUM_PRIMES; i++) {

		uint32 num_roots_r;
		uint32 num_roots_a;
		uint32 roots_r[MAX_POLY_DEGREE + 1];
		uint32 roots_a[MAX_POLY_DEGREE + 1];
		uint32 high_coeff;

		if (p + prime_delta[i] >= MAX_PACKED_PRIME)
			break;

		/* count the number of ideals; all rational ideals
		   are only counted once */

		p += prime_delta[i];
		num_roots_r = poly_get_zeros(roots_r, &fb->rfb.poly, p, 
					&high_coeff, 0);
		if (high_coeff == 0 || num_roots_r > 0)
			num_roots_r = 1;

		num_roots_a = poly_get_zeros(roots_a, &fb->afb.poly, p, 
					&high_coeff, 0);
		if (high_coeff == 0)
			roots_a[num_roots_a++] = p;

		/* if there's room in the array, save the ideals */

		if (num_ideals + num_roots_r + 
				num_roots_a >= MAX_SMALL_IDEALS)
			break;

		if (num_roots_r > 0) {
			small_ideals[num_ideals].rat_or_alg = RATIONAL_IDEAL;
			small_ideals[num_ideals].p_lo = (p - 1) / 2;
			small_ideals[num_ideals].p_hi = 0;
			small_ideals[num_ideals].r_lo = p;
			small_ideals[num_ideals].r_hi = 0;
			num_ideals++;
		}
		for (j = 0; j < num_roots_a; j++) {
			small_ideals[num_ideals].rat_or_alg = ALGEBRAIC_IDEAL;
			small_ideals[num_ideals].p_lo = (p - 1) / 2;
			small_ideals[num_ideals].p_hi = 0;
			small_ideals[num_ideals].r_lo = roots_a[j];
			small_ideals[num_ideals].r_hi = 0;
			num_ideals++;
		}
	}

	/* put the ideals in order of increasing size */

	qsort(small_ideals, (size_t)num_ideals, 
			sizeof(ideal_t), compare_ideals);

	*num_ideals_out = num_ideals;
	*max_small_ideal = p;
	return small_ideals;
}

/*------------------------------------------------------------------*/
#define QCB_VALS(r) ((r)->rel_index)
#define QCB_NUM_CHOICES 50000

static uint32 fill_qcb(msieve_obj *obj, mpz_poly_t *apoly, 
			relation_t *rlist, uint32 num_relations) {
	uint32 i, j;
	prime_sieve_t sieve;
	fb_entry_t qcb[QCB_SIZE];
	uint32 min_qcb_ideal = ((uint32)(-1) - QCB_NUM_CHOICES) &
					~(uint32)(0x1); /* must be even */
	uint8 bits[(QCB_NUM_CHOICES + 15) / 16] = {0};
	uint32 qcb_size;

	/* strike out all the algebraic factors of relations that
	   are between min_qcb_ideal and 2^32 */

	for (i = 0; i < num_relations; i++) {
		relation_t *r = rlist + i;
		uint32 num_r = r->num_factors_r;
		uint32 num_a = r->num_factors_a;
		uint32 array_size = 0;
		for (j = 0; j < num_r + num_a; j++) {
			uint64 p = decompress_p(r->factors, &array_size);
			if (j >= num_r && p > min_qcb_ideal &&
					p < ((uint64)1 << 32)) {
				uint32 pos = (p - min_qcb_ideal) / 2;
				bits[pos / 8] |= 1 << (pos % 8);
			}
		}
	}

	/* construct the quadratic character base, starting
	   with the primes above min_qcb_ideal */

	init_prime_sieve(&sieve, min_qcb_ideal + 1, 0xffffffff);

	i = 0;
	while (i < QCB_SIZE) {
		uint32 roots[MAX_POLY_DEGREE];
		uint32 num_roots, high_coeff;
		uint32 p = get_next_prime(&sieve);
		uint32 pos;

		if (p >= 0xffffff00)
			break;

		/* skip any p that appears in relations */

		pos = (p - min_qcb_ideal) / 2;
		if (bits[pos / 8] & (1 << (pos % 8)))
			continue;

		num_roots = poly_get_zeros(roots, apoly, p, 
					&high_coeff, 0);

		/* p cannot be a projective root of the algebraic poly */

		if (high_coeff == 0)
			continue;

		/* save the next batch of roots */

		for (j = 0; i < QCB_SIZE && j < num_roots; i++, j++) {
			qcb[i].p = p;
			qcb[i].r = roots[j];
		}
	}

	free_prime_sieve(&sieve);

	qcb_size = i;
	logprintf(obj, "using %u quadratic characters above %u\n",
				qcb_size, min_qcb_ideal + 1);

	/* cache each relation's quadratic characters for later use */

	for (i = 0; i < num_relations; i++) {
		relation_t *rel = rlist + i;
		int64 a = rel->a;
		uint32 b = rel->b;

		QCB_VALS(rel) = 0;
		for (j = 0; j < qcb_size; j++) {
			uint32 p = qcb[j].p;
			uint32 r = qcb[j].r;
			int64 res = a % (int64)p;
			int32 symbol;

			if (res < 0)
				res += (int64)p;

			symbol = mp_legendre_1(mp_modsub_1((uint32)res,
					mp_modmul_1(b, r, p), p), p);

			/* symbol must be 1 or -1; if it's 0,
			   there's something wrong with the choice
			   of primes in the QCB but this isn't
			   a fatal error */

			if (symbol == -1)
				QCB_VALS(rel) |= 1 << (j % 32);
			else if (symbol == 0)
				printf("warning: zero character\n");
		}
	}

	return qcb_size;
}

/*------------------------------------------------------------------*/
static uint32 combine_relations(la_col_t *col, relation_t *rlist,
				ideal_t *merged_ideals, uint32 *dense_rows,
				uint32 num_dense_rows, uint32 qcb_size) {

	uint32 i, j, k;
	uint32 num_merged = 0;
	ideal_t tmp_ideals[MAX_COL_IDEALS];
	uint32 num_tmp_ideals;

	/* form the matrix column corresponding to a 
	   collection of relations */

	for (i = 0; i < col->cycle.num_relations; i++) {
		relation_t *r = rlist + col->cycle.list[i];
		relation_lp_t new_ideals;

		/* fold in the quadratic characters for r */

		dense_rows[0] ^= QCB_VALS(r);

		/* if r is a free relation, modify the last dense row */

		if (r->b == 0) {
			dense_rows[(num_dense_rows - 1) / 32] ^=
					1 << ((num_dense_rows - 1) % 32);
		}

		/* get the ideal decomposition of relation i, sort
		   by size of prime */

		if (find_large_ideals(r, &new_ideals, 0, 0) > 
						TEMP_FACTOR_LIST_SIZE) {
			printf("error: overflow reading ideals\n");
			exit(-1);
		}
		if (num_merged + new_ideals.ideal_count >= MAX_COL_IDEALS) {
			printf("error: overflow merging ideals\n");
			exit(-1);
		}
		if (new_ideals.ideal_count > 1) {
			qsort(new_ideals.ideal_list, 
					(size_t)new_ideals.ideal_count,
					sizeof(ideal_t), compare_ideals);
		}

		/* merge it with the current list of ideals */

		j = k = 0;
		num_tmp_ideals = 0;
		while (j < num_merged && k < new_ideals.ideal_count) {
			int32 compare_result = compare_ideals(
						merged_ideals + j,
						new_ideals.ideal_list + k);
			if (compare_result < 0) {
				tmp_ideals[num_tmp_ideals++] = 
						merged_ideals[j++];
			}
			else if (compare_result > 0) {
				tmp_ideals[num_tmp_ideals++] = 
						new_ideals.ideal_list[k++];
			}
			else {
				j++; k++;
			}
		}
		while (j < num_merged) {
			tmp_ideals[num_tmp_ideals++] = merged_ideals[j++];
		}
		while (k < new_ideals.ideal_count) {
			tmp_ideals[num_tmp_ideals++] = 
						new_ideals.ideal_list[k++];
		}

		num_merged = num_tmp_ideals;
		memcpy(merged_ideals, tmp_ideals, 
				num_merged * sizeof(ideal_t));
	}

	/* fill in the parity row, and place at dense 
	   row position qcb_size */

	if (col->cycle.num_relations % 2)
		dense_rows[qcb_size / 32] |= 1 << (qcb_size % 32);

	return num_merged;
}

/*------------------------------------------------------------------*/
#define MAX_DENSE_ROW_WORDS 32

static void build_matrix_core(msieve_obj *obj, la_col_t *cycle_list, 
			uint32 num_cycles, relation_t *rlist, 
			uint32 num_relations, uint32 num_dense_rows, 
			ideal_t *small_ideals, uint32 num_small_ideals, 
			uint32 qcb_size, FILE *matrix_fp) {

	uint32 i, j, k;
	hashtable_t unique_ideals;
	uint32 max_small_ideal;
	uint32 dense_rows[MAX_DENSE_ROW_WORDS];
	uint32 dense_row_words;
	size_t mem_use;

	logprintf(obj, "building initial matrix\n");

	dense_row_words = (num_dense_rows + 31) / 32;
	if (dense_row_words > MAX_DENSE_ROW_WORDS) {
		printf("error: too many dense rows\n");
		exit(-1);
	}

	max_small_ideal = 0;
	if (num_small_ideals > 0) {
		max_small_ideal = small_ideals[
					num_small_ideals - 1].p_lo;
	}

	hashtable_init(&unique_ideals, (uint32)WORDS_IN(ideal_t), 0);

	fseek(matrix_fp, 3 * sizeof(uint32), SEEK_SET);

	/* for each cycle */

	for (i = 0; i < num_cycles; i++) {
		la_col_t *c = cycle_list + i;
		ideal_t merged_ideals[MAX_COL_IDEALS];
		uint32 mapped_ideals[MAX_COL_IDEALS];
		uint32 num_merged;

		/* dense rows start off empty */

		for (j = 0; j < dense_row_words; j++)
			dense_rows[j] = 0;

		/* merge the relations and quadratic characters
		   in the cycle */

		num_merged = combine_relations(c, rlist, merged_ideals, 
						dense_rows, num_dense_rows,
						qcb_size);

		/* assign a unique number to each ideal in 
		   the cycle. This will automatically ignore
		   empty rows in the matrix */

		for (j = k = 0; j < num_merged; j++) {
			ideal_t *ideal = merged_ideals + j;
			uint64 p = (uint64)ideal->p_hi << 32 | ideal->p_lo;

			if (max_small_ideal > 0 && (p == IDEAL_MINUS_ONE || 
					p <= max_small_ideal) ) {
				/* dense ideal; store in compressed format */
				ideal_t *loc = (ideal_t *)bsearch(ideal, 
						small_ideals,
						(size_t)num_small_ideals,
						sizeof(ideal_t),
						compare_ideals);
				uint32 idx = qcb_size + 1 +
						(loc - small_ideals);
				if (loc == NULL) {
					printf("error: unexpected dense "
						"ideal found\n");
					exit(-1);
				}
				dense_rows[idx / 32] |= 1 << (idx % 32);
			}
			else {
				uint32 idx;
				hashtable_find(&unique_ideals, 
						ideal, &idx, NULL);
				mapped_ideals[k++] = num_dense_rows + idx;
			}
		}

		/* save the matrix entries to disk */

		fwrite(&k, sizeof(uint32), (size_t)1, matrix_fp);
		fwrite(mapped_ideals, sizeof(uint32), (size_t)k, matrix_fp);
		fwrite(dense_rows, sizeof(uint32), 
				(size_t)dense_row_words, matrix_fp);
	}

	/* save the matrix dimensions to disk */

	i = num_dense_rows + hashtable_get_num(&unique_ideals);
	rewind(matrix_fp);
	fwrite(&i, sizeof(uint32), (size_t)1, matrix_fp);
	fwrite(&num_dense_rows, sizeof(uint32), (size_t)1, matrix_fp);
	fwrite(&num_cycles, sizeof(uint32), (size_t)1, matrix_fp);

	/* report memory use */

	mem_use = num_relations * sizeof(relation_t) +
			num_cycles * sizeof(la_col_t) +
			hashtable_sizeof(&unique_ideals);

	for (i = 0; i < num_cycles; i++) {
		la_col_t *c = cycle_list + i;
		mem_use += c->cycle.num_relations * sizeof(uint32);
	}
	for (i = 0; i < num_relations; i++) {
		relation_t *r = rlist + i;
		/* this is an upper bound */
		mem_use += (r->num_factors_r + r->num_factors_a) *
				sizeof(uint32);
	}
	logprintf(obj, "memory use: %.1f MB\n", (double)mem_use / 1048576);

	hashtable_free(&unique_ideals);
}

/*------------------------------------------------------------------*/
static void build_matrix(msieve_obj *obj, mpz_t n) {

	/* read in the relations that contribute to the
	   initial matrix, and form the quadratic characters
	   for each column */

	uint32 num_relations;
	relation_t *rlist;
	uint32 num_cycles;
	la_col_t *cycle_list;
	uint32 num_dense_rows;
	ideal_t *small_ideals;
	uint32 num_small_ideals;
	uint32 max_small_ideal;
	uint32 qcb_size;
	FILE *matrix_fp;
	char buf[256];
	factor_base_t fb;

	sprintf(buf, "%s.mat", obj->savefile.name);
	matrix_fp = fopen(buf, "w+b");
	if (matrix_fp == NULL) {
		logprintf(obj, "error: can't open matrix file '%s'\n", buf);
		exit(-1);
	}

	/* read in the NFS polynomials */

	memset(&fb, 0, sizeof(fb));
	mpz_poly_init(&fb.rfb.poly);
	mpz_poly_init(&fb.afb.poly);
	if (read_poly(obj, n, &fb.rfb.poly, &fb.afb.poly, NULL)) {
		printf("error: failed to read NFS polynomials\n");
		exit(-1);
	}

	/* fill in the small ideals */

	small_ideals = fill_small_ideals(&fb, &num_small_ideals,
					&max_small_ideal);

	/* read in the cycles that form the matrix columns,
	   and the relations they will need */

	nfs_read_cycles(obj, &fb, &num_cycles, &cycle_list, 
			&num_relations, &rlist, 1, 0);

	/* assign quadratic characters to each relation */

	qcb_size = fill_qcb(obj, &fb.afb.poly, rlist, num_relations);

	/* we need extra matrix rows to make sure that each
	   dependency has an even number of relations, and also an
	   even number of free relations. If the rational 
	   poly R(x) is monic, and we weren't using free relations,
	   the sign of R(x) is negative for all relations, meaning we 
	   already get the effect of the extra rows. However, in 
	   general we can't assume both of these are true */
	
	num_dense_rows = qcb_size + 2 + num_small_ideals;

	/* build the matrix columns, store to disk */

	build_matrix_core(obj, cycle_list, num_cycles, rlist, 
			num_relations, num_dense_rows, 
			small_ideals, num_small_ideals, 
			qcb_size, matrix_fp);

	nfs_free_relation_list(rlist, num_relations);
	free_cycle_list(cycle_list, num_cycles);
	free(small_ideals);
	fclose(matrix_fp);
	mpz_poly_free(&fb.rfb.poly);
	mpz_poly_free(&fb.afb.poly);
}

/*------------------------------------------------------------------*/
void nfs_solve_linear_system(msieve_obj *obj, mpz_t n) {

	/* convert the list of relations from the sieving 
	   stage into a matrix */

	uint32 i;
	la_col_t *cols;
	uint32 nrows; 
	uint32 max_nrows; 
	uint32 start_row;
	uint32 ncols; 
	uint32 max_ncols; 
	uint32 start_col;
	uint32 num_dense_rows;
	uint32 deps_found;
	uint64 *dependencies;
	uint32 skip_matbuild = 0;
	uint32 cado_filter = 0;
	time_t cpu_time = time(NULL);
#ifdef HAVE_MPI
	int32 grid_bools[2] = {0};
	int32 grid_dims[2];
	int32 mpi_nrows = 0;
	int32 mpi_ncols = 0;
#endif

	logprintf(obj, "\n");
	logprintf(obj, "commencing linear algebra\n");

	/* parse input arguments */

	if (obj->nfs_args != NULL) {
#ifdef HAVE_MPI
		const char *tmp0, *tmp1;

		tmp1 = strstr(obj->nfs_args, "mpi_nrows=");
		if (tmp1 != NULL)
			mpi_nrows = atoi(tmp1 + 10);

		tmp1 = strstr(obj->nfs_args, "mpi_ncols=");
		if (tmp1 != NULL)
			mpi_ncols = atoi(tmp1 + 10);

		/* old-style 'X,Y' format */

		tmp1 = strchr(obj->nfs_args, ',');
		if (tmp1 != NULL) {
			tmp0 = tmp1 - 1;
			while (tmp0 > obj->nfs_args && isdigit(tmp0[-1]))
				tmp0--;
			mpi_nrows = atoi(tmp0);
			mpi_ncols = atoi(tmp1 + 1);
		}
#endif
		if (strstr(obj->nfs_args, "skip_matbuild=1")) {
			logprintf(obj, "skipping matrix build\n");
			skip_matbuild = 1;
		}
		if (strstr(obj->nfs_args, "cado_filter=1")) {
			logprintf(obj, "assuming CADO-NFS filtering\n");
			cado_filter = 1;
		}
	}

#ifdef HAVE_MPI
	/* do not allow the square root to run if the
	   LA is running on a grid; this is a side-effect
	   of the LA not cleaning up very gracefully */

	if (obj->mpi_size > 1 && (obj->flags & MSIEVE_FLAG_NFS_SQRT)) {
		printf("error: square root is incompatible with "
			"multiple MPI processes\n");
		exit(-1);
	}

	/* create the grid */

	obj->mpi_nrows = grid_dims[0] = 1;
	obj->mpi_ncols = grid_dims[1] = obj->mpi_size;
	if (mpi_nrows && mpi_ncols) {
		if (obj->mpi_size != mpi_nrows * mpi_ncols) {
			printf("error: MPI size %u incompatible with "
				"%d x %d grid\n", obj->mpi_size, 
				mpi_nrows, mpi_ncols);
			MPI_Abort(MPI_COMM_WORLD, MPI_ERR_TOPOLOGY);
		}
		obj->mpi_nrows = grid_dims[0] = mpi_nrows;
		obj->mpi_ncols = grid_dims[1] = mpi_ncols;
	}
	if (obj->mpi_nrows > MAX_MPI_GRID_DIM ||
	    obj->mpi_ncols > MAX_MPI_GRID_DIM) {
		printf("error: MPI grid can be at most %u on a side\n",
			MAX_MPI_GRID_DIM);
		MPI_Abort(MPI_COMM_WORLD, MPI_ERR_TOPOLOGY);
	}

	MPI_TRY(MPI_Cart_create(MPI_COMM_WORLD, 2, grid_dims,
			grid_bools, 1, &obj->mpi_la_grid))

	/* get the rank of the current process in the grid */

	MPI_TRY(MPI_Comm_rank(obj->mpi_la_grid, (int *)&i))

	/* convert to grid coordinates */

	MPI_TRY(MPI_Cart_coords(obj->mpi_la_grid, i, 2, grid_dims))
	obj->mpi_la_row_rank = grid_dims[0];
	obj->mpi_la_col_rank = grid_dims[1];

	/* build communicators for the current row and column */

	grid_bools[0] = 1;
	grid_bools[1] = 0;
	MPI_TRY(MPI_Cart_sub(obj->mpi_la_grid, grid_bools, 
				&obj->mpi_la_col_grid))

	grid_bools[0] = 0;
	grid_bools[1] = 1;
	MPI_TRY(MPI_Cart_sub(obj->mpi_la_grid, grid_bools, 
				&obj->mpi_la_row_grid))

	logprintf(obj, "initialized process (%u,%u) of %u x %u grid\n", 
			obj->mpi_la_row_rank, obj->mpi_la_col_rank,
			obj->mpi_nrows, obj->mpi_ncols);
#endif

	if (!skip_matbuild && !(obj->flags & MSIEVE_FLAG_NFS_LA_RESTART)) {

		/* build the matrix; if using MPI, only process
		   0 does this, the rest are stalled. This isn't very
		   elegant, but avoiding it means either solving the
		   matrix twice, with the first pass having foreknowledge 
		   of the number of MPI processes that will eventually 
		   be used, or doing it in one pass with the matrix build
		   occurring in parallel. That actually is a nice idea
		   but would need a lot more more memory (a distributed
		   hashtable, or multiple copies of all the relations
		   involved) */

#ifdef HAVE_MPI
		if (obj->mpi_la_row_rank + obj->mpi_la_col_rank == 0) {
#endif
		uint64 sparse_weight;

		if (cado_filter)
			nfs_convert_cado_cycles(obj);

		/* build the initial matrix that is the output from
		   the filtering */

		build_matrix(obj, n);

		/* read the matrix and the list of cycles into memory
		   again, now that the underlying relations have been freed */

		read_matrix(obj, &nrows, NULL, NULL, &num_dense_rows, 
				&ncols, NULL, NULL, &cols, NULL, NULL);
		read_cycles(obj, &ncols, &cols, 0, NULL);

		count_matrix_nonzero(obj, nrows, num_dense_rows, ncols, cols);

		/* perform light filtering on the matrix */

		sparse_weight = reduce_matrix(obj, &nrows, num_dense_rows, 
				&ncols, cols, NUM_EXTRA_RELATIONS);
		if (ncols == 0) {
			logprintf(obj, "matrix is corrupt; skipping "
					"linear algebra\n");
			free(cols);
			return;
		}

		/* save the reduced matrix on disk; if MPI is configured,
		   also save the file offsets where each MPI process will
		   begin reading its own slab of matrix columns */

		dump_matrix(obj, nrows, num_dense_rows, 
				ncols, cols, sparse_weight);

		/* free the matrix */
		for (i = 0; i < ncols; i++) {
			free(cols[i].data);
			free(cols[i].cycle.list);
		}
		free(cols);
#if 0
		/* optimize the layout of large matrices */
		if (ncols > MIN_REORDER_SIZE) {

			uint32 *rowperm;
			uint32 *colperm;

			/* permute the rows and columns to concentrate
			   the nonzeros in specific places */

			reorder_matrix(obj, &rowperm, &colperm);

			/* read the matrix back into memory, applying
			   the permutation in the process */

			read_matrix(obj, &nrows, NULL, NULL, 
					&num_dense_rows, &ncols, 
					NULL, NULL, &cols, rowperm, colperm);
			read_cycles(obj, &ncols, &cols, 0, colperm);

			/* save the permuted matrix */

			dump_matrix(obj, nrows, num_dense_rows, 
					ncols, cols, sparse_weight);

			/* free everything */
			free(rowperm);
			free(colperm);
			for (i = 0; i < ncols; i++) {
				free(cols[i].data);
				free(cols[i].cycle.list);
			}
			free(cols);
		}
#endif

#ifdef HAVE_MPI
		}
		MPI_TRY(MPI_Barrier(obj->mpi_la_grid))
#endif
	}

	/* read the matrix in; if configured for MPI, this reads
	   in only the submatrix used by the current MPI process.
	   Without MPI, this reads the whole matrix, ncols = max_ncols, 
	   nrows = max_nrows, and start_row = start_col = 0.
	
	   Do not read in the relation numbers, the Lanczos code
	   doesn't need them */

	read_matrix(obj, &nrows, &max_nrows, &start_row,
			&num_dense_rows, 
			&ncols, &max_ncols, &start_col,
			&cols, NULL, NULL);
	logprintf(obj, "matrix starts at (%u, %u)\n", start_row, start_col);
	count_matrix_nonzero(obj, nrows, num_dense_rows, ncols, cols);

	/* solve the linear system */

	dependencies = block_lanczos(obj, 
				nrows, max_nrows, start_row,
				num_dense_rows,
				ncols, max_ncols, start_col,
				cols, &deps_found);
	if (deps_found)
		dump_dependencies(obj, dependencies, max_ncols);
	free(dependencies);
	free(cols);

#ifdef HAVE_MPI
	MPI_TRY(MPI_Comm_free(&obj->mpi_la_grid))
	MPI_TRY(MPI_Comm_free(&obj->mpi_la_row_grid))
	MPI_TRY(MPI_Comm_free(&obj->mpi_la_col_grid))
#endif
	cpu_time = time(NULL) - cpu_time;
	logprintf(obj, "BLanczosTime: %u\n", (uint32)cpu_time);
}
