/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: lanczos.c 1051 2024-05-07 14:42:24Z vegamink $
--------------------------------------------------------------------*/

#include "lanczos.h"

#define DEFAULT_DUMP_INTERVAL 2000

#ifdef HAVE_MPI
	#define MPI_NODE_0_START if (obj->mpi_la_row_rank + \
				obj->mpi_la_col_rank == 0) {

	#define MPI_NODE_0_END }
#else
	#define MPI_NODE_0_START /* nothing */
	#define MPI_NODE_0_END /* nothing */
#endif

#define BIT0(x) {{(uint64)(1) << (x)}}
#define BIT1(x) {{(uint64)0, (uint64)(1) << (x)}}
#define BIT2(x) {{(uint64)0, (uint64)0, (uint64)(1) << (x)}}
#define BIT3(x) {{(uint64)0, (uint64)0, (uint64)0, (uint64)(1) << (x)}}

static const v_t bitmask[VBITS] = {

BIT0( 0), BIT0( 1), BIT0( 2), BIT0( 3), BIT0( 4), BIT0( 5), BIT0( 6), BIT0( 7),
BIT0( 8), BIT0( 9), BIT0(10), BIT0(11), BIT0(12), BIT0(13), BIT0(14), BIT0(15),
BIT0(16), BIT0(17), BIT0(18), BIT0(19), BIT0(20), BIT0(21), BIT0(22), BIT0(23),
BIT0(24), BIT0(25), BIT0(26), BIT0(27), BIT0(28), BIT0(29), BIT0(30), BIT0(31),
BIT0(32), BIT0(33), BIT0(34), BIT0(35), BIT0(36), BIT0(37), BIT0(38), BIT0(39),
BIT0(40), BIT0(41), BIT0(42), BIT0(43), BIT0(44), BIT0(45), BIT0(46), BIT0(47),
BIT0(48), BIT0(49), BIT0(50), BIT0(51), BIT0(52), BIT0(53), BIT0(54), BIT0(55),
BIT0(56), BIT0(57), BIT0(58), BIT0(59), BIT0(60), BIT0(61), BIT0(62), BIT0(63),
#if VBITS > 64
BIT1( 0), BIT1( 1), BIT1( 2), BIT1( 3), BIT1( 4), BIT1( 5), BIT1( 6), BIT1( 7),
BIT1( 8), BIT1( 9), BIT1(10), BIT1(11), BIT1(12), BIT1(13), BIT1(14), BIT1(15),
BIT1(16), BIT1(17), BIT1(18), BIT1(19), BIT1(20), BIT1(21), BIT1(22), BIT1(23),
BIT1(24), BIT1(25), BIT1(26), BIT1(27), BIT1(28), BIT1(29), BIT1(30), BIT1(31),
BIT1(32), BIT1(33), BIT1(34), BIT1(35), BIT1(36), BIT1(37), BIT1(38), BIT1(39),
BIT1(40), BIT1(41), BIT1(42), BIT1(43), BIT1(44), BIT1(45), BIT1(46), BIT1(47),
BIT1(48), BIT1(49), BIT1(50), BIT1(51), BIT1(52), BIT1(53), BIT1(54), BIT1(55),
BIT1(56), BIT1(57), BIT1(58), BIT1(59), BIT1(60), BIT1(61), BIT1(62), BIT1(63),
#if VBITS > 128
BIT2( 0), BIT2( 1), BIT2( 2), BIT2( 3), BIT2( 4), BIT2( 5), BIT2( 6), BIT2( 7),
BIT2( 8), BIT2( 9), BIT2(10), BIT2(11), BIT2(12), BIT2(13), BIT2(14), BIT2(15),
BIT2(16), BIT2(17), BIT2(18), BIT2(19), BIT2(20), BIT2(21), BIT2(22), BIT2(23),
BIT2(24), BIT2(25), BIT2(26), BIT2(27), BIT2(28), BIT2(29), BIT2(30), BIT2(31),
BIT2(32), BIT2(33), BIT2(34), BIT2(35), BIT2(36), BIT2(37), BIT2(38), BIT2(39),
BIT2(40), BIT2(41), BIT2(42), BIT2(43), BIT2(44), BIT2(45), BIT2(46), BIT2(47),
BIT2(48), BIT2(49), BIT2(50), BIT2(51), BIT2(52), BIT2(53), BIT2(54), BIT2(55),
BIT2(56), BIT2(57), BIT2(58), BIT2(59), BIT2(60), BIT2(61), BIT2(62), BIT2(63),
#if VBITS > 192
BIT3( 0), BIT3( 1), BIT3( 2), BIT3( 3), BIT3( 4), BIT3( 5), BIT3( 6), BIT3( 7),
BIT3( 8), BIT3( 9), BIT3(10), BIT3(11), BIT3(12), BIT3(13), BIT3(14), BIT3(15),
BIT3(16), BIT3(17), BIT3(18), BIT3(19), BIT3(20), BIT3(21), BIT3(22), BIT3(23),
BIT3(24), BIT3(25), BIT3(26), BIT3(27), BIT3(28), BIT3(29), BIT3(30), BIT3(31),
BIT3(32), BIT3(33), BIT3(34), BIT3(35), BIT3(36), BIT3(37), BIT3(38), BIT3(39),
BIT3(40), BIT3(41), BIT3(42), BIT3(43), BIT3(44), BIT3(45), BIT3(46), BIT3(47),
BIT3(48), BIT3(49), BIT3(50), BIT3(51), BIT3(52), BIT3(53), BIT3(54), BIT3(55),
BIT3(56), BIT3(57), BIT3(58), BIT3(59), BIT3(60), BIT3(61), BIT3(62), BIT3(63),
#endif
#endif
#endif

};

/*-------------------------------------------------------------------*/
static uint32 form_post_lanczos_matrix(msieve_obj *obj, uint32 *nrows, 
				uint32 *dense_rows_out, 
				uint32 ncols, la_col_t *cols,
				v_t **post_lanczos_matrix) {

	uint32 i, j, k;
	uint32 num_dense_rows = *dense_rows_out;
	uint32 dense_row_words;
	uint32 new_dense_rows;
	uint32 new_dense_row_words;
	uint32 final_dense_row_words;
	v_t *submatrix;
	mp_t tmp;

	/* if the matrix is going to have cache blocking applied,
	   proceed but do not form a post-Lanczos matrix if one
	   is not desired. We have to do this because the block
	   matrix multiply expects the number of dense rows to be
	   a multiple of 64.

	   Otherwise, don't do anything if the Lanczos iteration 
	   would finish quickly */

	submatrix = NULL;

#ifdef HAVE_MPI
	/* for a 2-D grid of MPI processes, only the top row
	   of processes construct a post-lanczos matrix */

	if (obj->mpi_la_row_rank != 0)
		return 0;
#endif

	if (*nrows >= MIN_NROWS_TO_PACK ||
	    (POST_LANCZOS_ROWS > 0 && *nrows >= MIN_POST_LANCZOS_DIM)) {

		if (POST_LANCZOS_ROWS > 0) {
			logprintf(obj, "saving the first %u matrix rows "
					"for later\n", POST_LANCZOS_ROWS);
			submatrix = (v_t *)xmalloc(ncols * sizeof(v_t));
		}
	}
	else {
		return 0;
	}

	dense_row_words = (num_dense_rows + 31) / 32;
	mp_clear(&tmp);

	/* we will be removing the first POST_LANCZOS_ROWS rows
	   from the matrix entirely, and packing together the
	   next few rows. The matrix may have dense rows already, 
	   or these rows may be partially or completely sparse, 
	   in which case we'll have to pack them manually. After
	   the post-lanczos rows are removed, the number of dense 
	   rows remaining is a multiple of VBITS (minimum of VBITS) */

	new_dense_rows = MAX(num_dense_rows, POST_LANCZOS_ROWS);
	new_dense_rows += VBITS - (new_dense_rows - POST_LANCZOS_ROWS) % VBITS;
	new_dense_row_words = (new_dense_rows + 31) / 32;
	final_dense_row_words = (new_dense_rows - POST_LANCZOS_ROWS) / 32;

	for (i = 0; i < ncols; i++) {
		uint32 curr_weight = cols[i].weight;
		uint32 *curr_row = cols[i].data;

		/* build up a bitfield of the rows that will be
		   stored in packed format. Start with the rows
		   that are already packed */

		for (j = 0; j < dense_row_words; j++)
			tmp.val[j] = curr_row[curr_weight + j];

		/* add in the rows from the sparse part of the matrix.
		   Entries from these rows are either added to the
		   new dense bitfield, or moved to fill the holes
		   created by packing the first few sparse rows. In 
		   the latter case, the row index must be biased to 
		   reflect the removed rows */

		for (; j < new_dense_row_words; j++)
			tmp.val[j] = 0;

		for (j = k = 0; j < curr_weight; j++) {
			uint32 curr_index = curr_row[j];

			if (curr_index < new_dense_rows)
				tmp.val[curr_index / 32] |= 
						(uint32)1 << (curr_index % 32);
			else
				curr_row[k++] = curr_index - POST_LANCZOS_ROWS;
		}

		tmp.nwords = new_dense_row_words;
#if POST_LANCZOS_ROWS > 0
		/* remove the first POST_LANCZOS_ROWS bits from
		   the bitfield */
		for (j = 0; j < (POST_LANCZOS_ROWS + 63) / 64; j++) {
			submatrix[i].w[j] = ((uint64)tmp.val[2*j] |
						(uint64)tmp.val[2*j+1] << 32);
		}
#endif

		/* move the rest of the bitfield and repack the (hopefully
		   shorter) current column in the heap */
		cols[i].weight = k;
		if (k + final_dense_row_words > 0) {
			cols[i].data = (uint32 *)xrealloc(curr_row, (k + 
						final_dense_row_words) * 
						sizeof(uint32));
			mp_rshift(&tmp, POST_LANCZOS_ROWS, &tmp);
			memcpy(cols[i].data + k, tmp.val, 
						final_dense_row_words * 
						sizeof(uint32));
		}
		else {
			free(cols[i].data);
			cols[i].data = NULL;
		}
	}

	*nrows -= POST_LANCZOS_ROWS;
	*dense_rows_out = new_dense_rows - POST_LANCZOS_ROWS;
	*post_lanczos_matrix = submatrix;
	return (submatrix != NULL);
}

/*-------------------------------------------------------------------*/
static void mul_BxB_BxB(v_t *a, v_t *b, v_t *c ) {

	/* c[][] = x[][] * y[][], where all operands are 
	   VBITS x VBITS (i.e. contain VBITS vectors of
	   VBITS bits each). The result may overwrite a or b. */

	uint32 i, j;
	v_t tmp[VBITS];

	for (i = 0; i < VBITS; i++) {
		v_t accum = v_zero;

		for (j = 0; j < VWORDS; j++) {
			uint64 aij = a[i].w[j];
			uint32 k = 0;

			while (aij) {
				if (aij & 1)
					accum = v_xor(accum, b[64*j+k]);
				aij >>= 1;
				k++;
			}
		}

		tmp[i] = accum;
	}
	memcpy(c, tmp, sizeof(tmp));
}

/*-----------------------------------------------------------------------*/
static void transpose_BxB(v_t *a, v_t *b) {

	uint32 i, j;
	v_t tmp[VBITS];

	memset(tmp, 0, sizeof(tmp));

	for (i = 0; i < VBITS; i++) {
		v_t mask = bitmask[i];

		for (j = 0; j < VWORDS; j++) {
			uint64 word = a[i].w[j];
			uint32 k = 0;
			while (word) {
				if (word & 1)
					tmp[64*j+k] = v_or(tmp[64*j+k], mask);
				word >>= 1;
				k++;
			}
		}
	}
	memcpy(b, tmp, sizeof(tmp));
}

/*-------------------------------------------------------------------*/
static uint32 find_nonsingular_sub(msieve_obj *obj,
				v_t *t, uint32 *s, 
				uint32 *last_s, uint32 last_dim, 
				v_t *w) {

	/* given a VBITSxVBITS matrix t[][] and a list 
	   of 'last_dim' column indices enumerated in last_s[]: 
	   
	     - find a submatrix of t that is invertible 
	     - invert it and copy to w[][]
	     - enumerate in s[] the columns represented in w[][] */

	uint32 i, j;
	uint32 dim, curr_dim;
	uint32 cols[VBITS];
	v_t M[VBITS][2];
	v_t mask, *row_i, *row_j;
	v_t m0, m1;

	/* M = [t | I] for I the VBITS x VBITS identity matrix */

	for (i = 0; i < VBITS; i++) {
		M[i][0] = t[i]; 
		M[i][1] = bitmask[i];
	}

	/* put the column indices from last_s[] into the
	   back of cols[], and copy to the beginning of cols[]
	   any column indices not in last_s[] */

	mask = v_zero;
	for (i = 0; i < last_dim; i++) {
		cols[VBITS - 1 - i] = last_s[i];
		mask = v_or(mask, bitmask[last_s[i]]);
	}
	for (i = j = 0; i < VBITS; i++) {
		if (!v_bitset(mask, i))
			cols[j++] = i;
	}

	/* compute the inverse of t[][] */

	for (i = dim = 0; i < VBITS; i++) {
	
		/* find the next pivot row and put in row i */

		curr_dim = cols[i];
		row_i = M[curr_dim];

		for (j = i; j < VBITS; j++) {
			row_j = M[cols[j]];
			if (v_bitset(row_j[0], curr_dim)) {
				m0 = row_j[0];
				m1 = row_j[1];
				row_j[0] = row_i[0];
				row_j[1] = row_i[1];
				row_i[0] = m0; 
				row_i[1] = m1;
				break;
			}
		}
				
		/* if a pivot row was found, eliminate the pivot
		   column from all other rows */

		if (j < VBITS) {
			for (j = 0; j < VBITS; j++) {
				row_j = M[cols[j]];
				if (row_i != row_j && 
				    v_bitset(row_j[0], curr_dim)) {
					row_j[0] = v_xor(row_j[0], row_i[0]);
					row_j[1] = v_xor(row_j[1], row_i[1]);
				}
			}

			/* add the pivot column to the list of 
			   accepted columns */

			s[dim++] = curr_dim;
			continue;
		}

		/* otherwise, use the right-hand half of M[]
		   to compensate for the absence of a pivot column */

		for (j = i; j < VBITS; j++) {
			row_j = M[cols[j]];
			if (v_bitset(row_j[1], curr_dim)) {
				m0 = row_j[0];
				m1 = row_j[1];
				row_j[0] = row_i[0];
				row_j[1] = row_i[1];
				row_i[0] = m0; 
				row_i[1] = m1;
				break;
			}
		}
				
		if (j == VBITS) {
			logprintf(obj, "lanczos error: submatrix "
					"is not invertible\n");
			return 0;
		}
			
		/* eliminate the pivot column from the other rows
		   of the inverse */

		for (j = 0; j < VBITS; j++) {
			row_j = M[cols[j]];
			if (row_i != row_j && 
			    v_bitset(row_j[1], curr_dim)) {
				row_j[0] = v_xor(row_j[0], row_i[0]);
				row_j[1] = v_xor(row_j[1], row_i[1]);
			}
		}

		/* wipe out the pivot row */

		row_i[0] = row_i[1] = v_zero;
	}

	/* the right-hand half of M[] is the desired inverse */
	
	for (i = 0; i < VBITS; i++) 
		w[i] = M[i][1];

	return dim;
}

/*-----------------------------------------------------------------------*/
static void transpose_vector(uint32 ncols, v_t *v, uint64 **trans) {

	/* Hideously inefficient routine to transpose a
	   vector v[] of VBITS-size words into a 2-D array
	   trans[][] of 64-bit words */

	uint32 i, j, k;
	uint32 col;
	uint64 mask, word;

	for (i = 0; i < ncols; i++) {
		col = i / 64;
		mask = bitmask[i % 64].w[0];
		for (j = 0; j < VWORDS; j++) {
			word = v[i].w[j];
			k = 0;
			while (word) {
				if (word & 1) {
					trans[64*j+k][col] |= mask;
				}
				word = word >> 1;
				k++;
			}
		}
	}
}

/*-----------------------------------------------------------------------*/
static uint32 combine_cols(uint32 ncols, 
			v_t *x, v_t *v, v_t *ax, v_t *av) {

	/* Once the block Lanczos iteration has finished, 
	   x[] and v[] will contain mostly nullspace vectors
	   between them, as well as possibly some columns
	   that are linear combinations of nullspace vectors.
	   Given vectors ax[] and av[] that are the result of
	   multiplying x[] and v[] by the matrix, this routine 
	   will use Gauss elimination on the columns of [ax | av] 
	   to find all of the linearly dependent columns. The
	   column operations needed to accomplish this are mir-
	   rored in [x | v] and the columns that are independent
	   are skipped. Finally, the dependent columns are copied
	   back into x[] and represent the nullspace vector output
	   of the block Lanczos code. */

	uint32 i, j, k, bitpos, col, col_words;
	uint64 mask;
	uint64 *matrix[2*VBITS], *amatrix[2*VBITS], *tmp;

	col_words = (ncols + 63) / 64;

	for (i = 0; i < 2*VBITS; i++) {
		matrix[i] = (uint64 *)xcalloc((size_t)col_words, 
					     sizeof(uint64));
		amatrix[i] = (uint64 *)xcalloc((size_t)col_words, 
					      sizeof(uint64));
	}

	/* operations on columns can more conveniently become 
	   operations on rows if all the vectors are first
	   transposed */

	transpose_vector(ncols, x, matrix);
	transpose_vector(ncols, ax, amatrix);
	transpose_vector(ncols, v, matrix + VBITS);
	transpose_vector(ncols, av, amatrix + VBITS);

	/* Keep eliminating rows until the unprocessed part
	   of amatrix[][] is all zero. The rows where this
	   happens correspond to linearly dependent vectors
	   in the nullspace */

	for (i = bitpos = 0; i < 2*VBITS && bitpos < ncols; bitpos++) {

		/* find the next pivot row */

		mask = bitmask[bitpos % 64].w[0];
		col = bitpos / 64;
		for (j = i; j < 2*VBITS; j++) {
			if (amatrix[j][col] & mask) {
				tmp = matrix[i];
				matrix[i] = matrix[j];
				matrix[j] = tmp;
				tmp = amatrix[i];
				amatrix[i] = amatrix[j];
				amatrix[j] = tmp;
				break;
			}
		}
		if (j == 2*VBITS)
			continue;

		/* a pivot was found; eliminate it from the
		   remaining rows */

		for (j++; j < 2*VBITS; j++) {
			if (amatrix[j][col] & mask) {

				/* Note that the entire row, *not*
				   just the nonzero part of it, must
				   be eliminated; this is because the
				   corresponding (dense) row of matrix[][]
				   must have the same operation applied */

				for (k = 0; k < col_words; k++) {
					amatrix[j][k] ^= amatrix[i][k];
					matrix[j][k] ^= matrix[i][k];
				}
			}
		}
		i++;
	}

	/* transpose rows i to VBITS back into x[]. Pack the
	   dependencies into the low-order bits of x[] */

	for (j = 0; j < ncols; j++) {
		v_t word = v_zero;

		col = j / 64;
		mask = bitmask[j % 64].w[0];

		for (k = i; k < VBITS; k++) {
			if (matrix[k][col] & mask) {
				word = v_or(word, bitmask[k - i]);
			}
		}
		x[j] = word;
	}

	for (j = 0; j < 2*VBITS; j++) {
		free(matrix[j]);
		free(amatrix[j]);
	}

	if (i > VBITS)
		return 0;
	return VBITS - i;
}

/*-----------------------------------------------------------------------*/
#ifdef HAVE_MPI
static v_t * gather_ncols(msieve_obj *obj,
			packed_matrix_t *packed_matrix,
			void *v, void *scratch, 
			v_t *out) {

	MPI_NODE_0_START
	if (out == NULL)
		out = (v_t *)aligned_malloc(packed_matrix->max_ncols * 
						sizeof(v_t), 64);
	MPI_NODE_0_END

	/* gather v into MPI row 0 */

	MPI_TRY(MPI_Gatherv(v,
			VWORDS * packed_matrix->nsubcols, 
			MPI_LONG_LONG, scratch,
			packed_matrix->subcol_counts,
			packed_matrix->subcol_offsets,
			MPI_LONG_LONG, 0, 
			obj->mpi_la_col_grid))

	/* gather row 0 into the root node */

	if (obj->mpi_la_row_rank == 0) {
		MPI_TRY(MPI_Gatherv(scratch,
				VWORDS * packed_matrix->ncols, 
				MPI_LONG_LONG, out,
				packed_matrix->col_counts,
				packed_matrix->col_offsets,
				MPI_LONG_LONG, 0, 
				obj->mpi_la_row_grid))
	}
	
	return out;
}

/*-----------------------------------------------------------------------*/
static v_t * gather_nrows(msieve_obj *obj,
			packed_matrix_t *packed_matrix,
			void *scratch, v_t *out) {

	MPI_NODE_0_START
	if (out == NULL)
		out = (v_t *)aligned_malloc(packed_matrix->max_ncols * 
						sizeof(v_t), 64);
	MPI_NODE_0_END

	/* gather column 0 into the root node */

	if (obj->mpi_la_col_rank == 0) {
		MPI_TRY(MPI_Gatherv(scratch,
				VWORDS * packed_matrix->nrows, 
				MPI_LONG_LONG, out,
				packed_matrix->row_counts,
				packed_matrix->row_offsets,
				MPI_LONG_LONG, 0, 
				obj->mpi_la_col_grid))
	}
	
	return out;
}

/*-----------------------------------------------------------------------*/
static void scatter_ncols(msieve_obj *obj,
			packed_matrix_t *packed_matrix,
			void *out, void *scratch, 
			v_t *in) {

	/* push out to the top MPI row */

	if (obj->mpi_la_row_rank == 0)
		MPI_TRY(MPI_Scatterv(in, packed_matrix->col_counts,
				packed_matrix->col_offsets, 
				MPI_LONG_LONG, scratch,
				VWORDS * packed_matrix->ncols,
				MPI_LONG_LONG, 0, 
				obj->mpi_la_row_grid))

	/* push down each column */

	MPI_TRY(MPI_Scatterv(scratch, packed_matrix->subcol_counts,
	       			packed_matrix->subcol_offsets, 
	      			MPI_LONG_LONG, out,
				VWORDS * packed_matrix->ncols,
	   			MPI_LONG_LONG, 0, 
       				obj->mpi_la_col_grid))
}
#endif

/*-----------------------------------------------------------------------*/
static void dump_lanczos_state(msieve_obj *obj, 
			packed_matrix_t *packed_matrix,
			void *x, v_t **vt_v0, void **v, void *v0,
			v_t **vt_a_v, v_t **vt_a2_v, v_t **winv,
			uint32 n, uint32 max_n, uint32 dim_solved, uint32 iter,
			uint32 s[2][VBITS], uint32 dim1, void *scratch) {

	uint32 vbits = VBITS;
	char buf[256];
	char buf_old[256];
	char buf_bak[256];
	FILE *fp;
	uint32 status = 1;
	v_t *tmp = NULL;

	MPI_NODE_0_START

	sprintf(buf, "%s.chk0", obj->savefile.name);
	sprintf(buf_old, "%s.chk", obj->savefile.name);
	sprintf(buf_bak, "%s.bak.chk", obj->savefile.name);
	fp = fopen(buf, "wb");
	if (fp == NULL) {
		printf("error: cannot open matrix checkpoint file\n");
		exit(-1);
	}

	status &= (fwrite(&max_n, sizeof(uint32), (size_t)1, fp) == 1);
	status &= (fwrite(&dim_solved, sizeof(uint32), (size_t)1, fp) == 1);
	status &= (fwrite(&iter, sizeof(uint32), (size_t)1, fp) == 1);
	status &= (fwrite(&vbits, sizeof(uint32), (size_t)1, fp) == 1);

	status &= (fwrite(vt_a_v[1], sizeof(v_t), (size_t)VBITS, fp) == VBITS);
	status &= (fwrite(vt_a2_v[1], sizeof(v_t), (size_t)VBITS, fp) == VBITS);
	status &= (fwrite(winv[1], sizeof(v_t), (size_t)VBITS, fp) == VBITS);
	status &= (fwrite(winv[2], sizeof(v_t), (size_t)VBITS, fp) == VBITS);
	status &= (fwrite(vt_v0[0], sizeof(v_t), (size_t)VBITS, fp) == VBITS);
	status &= (fwrite(vt_v0[1], sizeof(v_t), (size_t)VBITS, fp) == VBITS);
	status &= (fwrite(vt_v0[2], sizeof(v_t), (size_t)VBITS, fp) == VBITS);
	status &= (fwrite(s[1], sizeof(uint32), (size_t)VBITS, fp) == VBITS);
	status &= (fwrite(&dim1, sizeof(uint32), (size_t)1, fp) == 1);
	tmp = (v_t *)aligned_malloc(max_n * sizeof(v_t), 64);
	MPI_NODE_0_END

#ifdef HAVE_MPI
	tmp = gather_ncols(obj, packed_matrix, x, scratch, tmp);
	MPI_NODE_0_START
	status &= (fwrite(tmp, sizeof(v_t), (size_t)max_n, fp) == max_n);
	MPI_NODE_0_END

	tmp = gather_ncols(obj, packed_matrix, v[0], scratch, tmp);
	MPI_NODE_0_START
	status &= (fwrite(tmp, sizeof(v_t), (size_t)max_n, fp) == max_n);
	MPI_NODE_0_END

	tmp = gather_ncols(obj, packed_matrix, v[1], scratch, tmp);
	MPI_NODE_0_START
	status &= (fwrite(tmp, sizeof(v_t), (size_t)max_n, fp) == max_n);
	MPI_NODE_0_END

	tmp = gather_ncols(obj, packed_matrix, v[2], scratch, tmp);
	MPI_NODE_0_START
	status &= (fwrite(tmp, sizeof(v_t), (size_t)max_n, fp) == max_n);
	MPI_NODE_0_END

	tmp = gather_ncols(obj, packed_matrix, v0, scratch, tmp);
	MPI_NODE_0_START
	status &= (fwrite(tmp, sizeof(v_t), (size_t)max_n, fp) == max_n);
	MPI_NODE_0_END
#else
	vv_copyout(tmp, x, max_n);
	status &= (fwrite(tmp, sizeof(v_t), (size_t)max_n, fp) == max_n);

	vv_copyout(tmp, v[0], max_n);
	status &= (fwrite(tmp, sizeof(v_t), (size_t)max_n, fp) == max_n);

	vv_copyout(tmp, v[1], max_n);
	status &= (fwrite(tmp, sizeof(v_t), (size_t)max_n, fp) == max_n);

	vv_copyout(tmp, v[2], max_n);
	status &= (fwrite(tmp, sizeof(v_t), (size_t)max_n, fp) == max_n);

	vv_copyout(tmp, v0, max_n);
	status &= (fwrite(tmp, sizeof(v_t), (size_t)max_n, fp) == max_n);
#endif

	MPI_NODE_0_START
	aligned_free(tmp);
	fclose(fp);

	/* only delete an old checkpoint file if the current 
	   checkpoint completed writing. More paranoid: compute a 
	   cryptographic hash of the file and then verify against 
	   the disk image */

	if (status == 0) {
		printf("error: cannot write new checkpoint file\n");
		printf("error: previous checkpoint file not overwritten\n");
		exit(-1);
	}
	remove(buf_bak);
	rename(buf_old, buf_bak);
	if (rename(buf, buf_old)) {
		printf("error: cannot update checkpoint file\n");
		exit(-1);
	}

	MPI_NODE_0_END
}

/*-----------------------------------------------------------------------*/
static void read_lanczos_state(msieve_obj *obj, 
			packed_matrix_t *packed_matrix,
			void *x, v_t **vt_v0, void **v, void *v0,
			v_t **vt_a_v, v_t **vt_a2_v, v_t **winv,
			uint32 n, uint32 max_n, uint32 *dim_solved, 
			uint32 *iter, uint32 s[2][VBITS], uint32 *dim1,
			void *scratch) {

	uint32 read_n;
	uint32 status;
	char buf[256];
	FILE *fp;
	v_t *tmp = NULL;
	uint32 vbits = 0;

	sprintf(buf, "%s.chk", obj->savefile.name);
	fp = fopen(buf, "rb");
	if (fp == NULL) {
		printf("error: cannot open matrix checkpoint file\n");
		exit(-1);
	}

	status = 1;
	fread(&read_n, sizeof(uint32), (size_t)1, fp);
	if (read_n != max_n) {
		printf("error: unexpected vector size\n");
		exit(-1);
	}
	status &= (fread(dim_solved, sizeof(uint32), (size_t)1, fp) == 1);
	status &= (fread(iter, sizeof(uint32), (size_t)1, fp) == 1);
	status &= (fread(&vbits, sizeof(uint32), (size_t)1, fp) == 1);
	if (vbits != VBITS) {
		printf("error: vector length mismatch\n");
		exit(-1);
	}

	status &= (fread(vt_a_v[1], sizeof(v_t), (size_t)VBITS, fp) == VBITS);
	status &= (fread(vt_a2_v[1], sizeof(v_t), (size_t)VBITS, fp)== VBITS);
	status &= (fread(winv[1], sizeof(v_t), (size_t)VBITS, fp) == VBITS);
	status &= (fread(winv[2], sizeof(v_t), (size_t)VBITS, fp) == VBITS);
	status &= (fread(vt_v0[0], sizeof(v_t), (size_t)VBITS, fp) == VBITS);
	status &= (fread(vt_v0[1], sizeof(v_t), (size_t)VBITS, fp) == VBITS);
	status &= (fread(vt_v0[2], sizeof(v_t), (size_t)VBITS, fp) == VBITS);
	status &= (fread(s[1], sizeof(uint32), (size_t)VBITS, fp) == VBITS);
	status &= (fread(dim1, sizeof(uint32), (size_t)1, fp) == 1);

#ifdef HAVE_MPI
	MPI_NODE_0_START
	tmp = (v_t *)xmalloc(max_n * sizeof(v_t));
	status &= (fread(tmp, sizeof(v_t), (size_t)max_n, fp) == max_n);
	MPI_NODE_0_END
	scatter_ncols(obj, packed_matrix, x, scratch, tmp);

	MPI_NODE_0_START
	status &= (fread(tmp, sizeof(v_t), (size_t)max_n, fp) == max_n);
	MPI_NODE_0_END
	scatter_ncols(obj, packed_matrix, v[0], scratch, tmp);

	MPI_NODE_0_START
	status &= (fread(tmp, sizeof(v_t), (size_t)max_n, fp) == max_n);
	MPI_NODE_0_END
	scatter_ncols(obj, packed_matrix, v[1], scratch, tmp);

	MPI_NODE_0_START
	status &= (fread(tmp, sizeof(v_t), (size_t)max_n, fp) == max_n);
	MPI_NODE_0_END
	scatter_ncols(obj, packed_matrix, v[2], scratch, tmp);

	MPI_NODE_0_START
	status &= (fread(tmp, sizeof(v_t), (size_t)max_n, fp) == max_n);
	MPI_NODE_0_END
	scatter_ncols(obj, packed_matrix, v0, scratch, tmp);
#else
	tmp = (v_t *)xmalloc(max_n * sizeof(v_t));

	status &= (fread(tmp, sizeof(v_t), (size_t)max_n, fp) == max_n);
	vv_copyin(x, tmp, max_n);

	status &= (fread(tmp, sizeof(v_t), (size_t)max_n, fp) == max_n);
	vv_copyin(v[0], tmp, max_n);

	status &= (fread(tmp, sizeof(v_t), (size_t)max_n, fp) == max_n);
	vv_copyin(v[1], tmp, max_n);

	status &= (fread(tmp, sizeof(v_t), (size_t)max_n, fp) == max_n);
	vv_copyin(v[2], tmp, max_n);

	status &= (fread(tmp, sizeof(v_t), (size_t)max_n, fp) == max_n);
	vv_copyin(v0, tmp, max_n);
#endif
	free(tmp);

	fclose(fp);
	if (status == 0) {
		printf("error: checkpoint recovery failed\n");
		exit(-1);
	}
}

/*-----------------------------------------------------------------------*/
static void init_lanczos_state(msieve_obj *obj, 
			packed_matrix_t *packed_matrix, void *scratch,
			void *x, void *v0, v_t **vt_v0, void **v, 
			v_t **vt_a_v, v_t **vt_a2_v, v_t **winv,
			uint32 n, uint32 s[2][VBITS], uint32 *dim1) {

	uint32 i;
	v_t *tmp;

	/* The computed solution 'x' starts off random,
	   and v[0] starts off as B*x. This initial copy
	   of v[0] must be saved off separately */

#ifdef HAVE_MPI
	/* all nodes work with vectors of size ncols/P */
	n = packed_matrix->nsubcols;
#endif
	tmp = (v_t *)xmalloc(n * sizeof(v_t));
	for (i = 0; i < n; i++)
		tmp[i] = v_random(&obj->seed1, &obj->seed2);

	vv_copyin(x, tmp, n);
	free(tmp);

	mul_sym_NxN_NxB(packed_matrix, x, v[0], scratch);
	vv_copy(v0, v[0], n);

	/* Subscripts larger than zero represent past versions of 
	   these quantities, which start off empty (except for the 
	   past version of s[], which contains all the column 
	   indices) */
	   
	vv_clear(v[1], n);
	vv_clear(v[2], n);
	for (i = 0; i < VBITS; i++) {
		s[1][i] = i;
		vt_a_v[1][i] = v_zero;
		vt_a2_v[1][i] = v_zero;
		winv[1][i] = v_zero;
		winv[2][i] = v_zero;
		vt_v0[0][i] = v_zero;
		vt_v0[1][i] = v_zero;
		vt_v0[2][i] = v_zero;
	}
	*dim1 = VBITS;
}

/*-----------------------------------------------------------------------*/
static v_t * block_lanczos_core(msieve_obj *obj, 
				packed_matrix_t *packed_matrix,
				uint32 *num_deps_found,
				v_t *post_lanczos_matrix,
				uint32 dump_interval) {
	
	/* Solve Bx = 0 for some nonzero x; the computed
	   solution, containing up to 64 of these nullspace
	   vectors, is returned */

	uint32 n = packed_matrix->ncols;
	uint32 max_n = packed_matrix->max_ncols;
	void *vnext, *v[3], *x, *v0, *scratch, *tmp;
	v_t *out0, *out1, *out2, *out3;
	v_t *winv[3], *vt_v0_next;
	v_t *vt_a_v[2], *vt_a2_v[2], *vt_v0[3];
	uint32 s[2][VBITS];
	v_t d[VBITS], e[VBITS], f[VBITS], f2[VBITS];
	uint32 i; 
	uint32 dim0, dim1;
	v_t mask0, mask1;

	uint32 iter = 0;
	uint32 dim_solved = 0;
	uint32 first_dim_solved = 0;
	uint32 report_interval = 0;
	uint32 check_interval = 0;
	uint32 next_report = 0;
	uint32 log_eta_once = 0;
	uint32 next_check = 0;
	uint32 next_dump = 0;
	time_t first_time;

	if (packed_matrix->num_threads > 1)
		logprintf(obj, "commencing Lanczos iteration (%u threads)\n",
					packed_matrix->num_threads);
	else
		logprintf(obj, "commencing Lanczos iteration\n");

#ifdef HAVE_MPI
	
	/* all nodes work with vectors of size ncols/P */
	n = packed_matrix->nsubcols;
	   
	/* we'll need 2 scratch vectors for the matrix multiply
	   and for scatter-gather operations */

	scratch = vv_alloc(2 * MAX(packed_matrix->nrows, 
				packed_matrix->ncols),
			  packed_matrix->extra);
#else	
	/* without MPI, all vectors are the maximum size */
	scratch = vv_alloc(n, packed_matrix->extra);
#endif    

	v[0] = vv_alloc(n, packed_matrix->extra);
	v[1] = vv_alloc(n, packed_matrix->extra);
	v[2] = vv_alloc(n, packed_matrix->extra);
	vnext = vv_alloc(n, packed_matrix->extra);
	x = vv_alloc(n, packed_matrix->extra);
	v0 = vv_alloc(n, packed_matrix->extra);
    
	/* VBITSxVBITS data */

	winv[0] = (v_t *)aligned_malloc(VBITS * sizeof(v_t), 64);
	winv[1] = (v_t *)aligned_malloc(VBITS * sizeof(v_t), 64);
	winv[2] = (v_t *)aligned_malloc(VBITS * sizeof(v_t), 64);
	vt_a_v[0] = (v_t *)aligned_malloc(VBITS * sizeof(v_t), 64);
	vt_a_v[1] = (v_t *)aligned_malloc(VBITS * sizeof(v_t), 64);
	vt_a2_v[0] = (v_t *)aligned_malloc(VBITS * sizeof(v_t), 64);
	vt_a2_v[1] = (v_t *)aligned_malloc(VBITS * sizeof(v_t), 64);
	vt_v0[0] = (v_t *)aligned_malloc(VBITS * sizeof(v_t), 64);
	vt_v0[1] = (v_t *)aligned_malloc(VBITS * sizeof(v_t), 64);
	vt_v0[2] = (v_t *)aligned_malloc(VBITS * sizeof(v_t), 64);
	vt_v0_next = (v_t *)aligned_malloc(VBITS * sizeof(v_t), 64);

	logprintf(obj, "memory use: %.1f MB\n", (double)
			(packed_matrix_sizeof(packed_matrix)) / 1048576);

	logprintf(obj, "VBITS = %d\n", VBITS);

	/* initialize */

	*num_deps_found = 0;
	iter = 0;
	dim0 = 0;

	if (obj->flags & MSIEVE_FLAG_NFS_LA_RESTART) {
		read_lanczos_state(obj, packed_matrix, 
				x, vt_v0, v, v0, vt_a_v, vt_a2_v,
				winv, packed_matrix->ncols, max_n, 
				&dim_solved, &iter, s, &dim1, scratch);
		logprintf(obj, "restarting at iteration %u (dim = %u)\n",
				iter, dim_solved);
	}
	else {
		init_lanczos_state(obj, packed_matrix, scratch, x, 
				v0, vt_v0, v, vt_a_v, vt_a2_v, 
				winv, packed_matrix->ncols, s, &dim1);
	}

	mask1 = v_zero;
	for (i = 0; i < dim1; i++)
		mask1 = v_or(mask1, bitmask[s[1][i]]);

	/* determine if the solver will run long enough that
	   it would be worthwhile to report progress */

	first_time = time(NULL);
	if (max_n > 60000 &&
	    obj->flags & (MSIEVE_FLAG_USE_LOGFILE |
	    		  MSIEVE_FLAG_LOG_TO_STDOUT)) {
		if (max_n > 1000000)
			report_interval = 200;
		else if (max_n > 500000)
			report_interval = 500;
		else if (max_n > 100000)
			report_interval = 2000;
		else
			report_interval = 8000;
		first_dim_solved = dim_solved;
		next_report = dim_solved + report_interval;
	}

	if (dump_interval) {
		/* avoid check (at dump) within the next few iterations */
		next_dump = ((dim_solved + 6 * VBITS) / dump_interval + 1) * 
					dump_interval;
		check_interval = 10000;
		/* avoid next_check within 4*64 dim + some cushion */
		next_check = ((dim_solved + 6 * VBITS) / check_interval + 1) * 
					check_interval;
	}

	/* perform the iteration */

	while (1) {
		iter++;

		/* multiply the current v[0] by the matrix and write
		   to vnext */
              
		mul_sym_NxN_NxB(packed_matrix, v[0], vnext, scratch);
                
		/* compute v0'*A*v0 and (A*v0)'(A*v0) */

		vv_mul_BxN_NxB(packed_matrix, v[0], vnext, vt_a_v[0], n);
		vv_mul_BxN_NxB(packed_matrix, vnext, vnext, vt_a2_v[0], n);

		/* if the former is orthogonal to itself, then
		   the iteration has finished */

		for (i = 0; i < VBITS; i++) {
			if (!v_is_all_zeros(vt_a_v[0][i]))
				break;
		}
		if (i == VBITS)
			break;

		/* Find the size-'dim0' nonsingular submatrix
		   of v0'*A*v0, invert it, and list the column
		   indices present in the submatrix */

		dim0 = find_nonsingular_sub(obj, vt_a_v[0], s[0], 
					    s[1], dim1, winv[0]);
		if (dim0 == 0)
			break;

		/* mask0 contains one set bit for every column
		   that participates in the inverted submatrix
		   computed above */

		mask0 = v_zero;
		for (i = 0; i < dim0; i++)
			mask0 = v_or(mask0, bitmask[s[0][i]]);

		/* The block Lanczos recurrence depends on all columns
		   of v'Av appearing in the current and/or previous iteration. 
		   Verify that condition here
		  
		   Note that the test only applies if this is not the
		   last Lanczos iteration. I'm not sure that this is right, 
		   but the test fails on the last iteration much more often 
		   than would be expected by chance alone, yet ignoring
		   the failure still produces good dependencies 
		   
		   Note that the last iteration typically has dim_solved
		   slightly less than the number of rows, not the number
		   of columns (=n) */
	
		if (dim_solved < packed_matrix->max_nrows - VBITS) {
			if (!v_is_all_ones(v_or(mask0, mask1))) {
				logprintf(obj, "lanczos error (dim = %u): "
						"not all columns used\n",
						dim_solved);
				dim0 = 0;
				break;
			}
		}

		/* begin the computation of the next v. First mask
		   off the vectors that are included in this iteration */

		dim_solved += dim0;
		if (!v_is_all_ones(mask0)) {
			vv_mask(vnext, mask0, n);
		}

		/* begin the computation of the next v' * v0. For 
		   the first three iterations, this requires a full 
		   inner product. For all succeeding iterations, the 
		   next v' * v0 is the sum of three 64x64 products 
		   and is stored in vt_v0_next. */

		if (iter < 4) {
			vv_mul_BxN_NxB(packed_matrix, v[0], v0, vt_v0[0], n);
		}
		else if (iter == 4) {
			/* v0 is not needed from now on; recycle it 
			   for use as a check vector */
			vv_clear(v0, n);
		}

		/* perform an integrity check on the iteration. This 
		   verifies that the current value of vnext is orthogonal 
		   to the vnext that was computed about check_interval 
		   dimensions ago
		
		   Checks happen on a fixed schedule, as well as 
		   right before a checkpoint file is written */

		if (check_interval && (dim_solved >= next_check ||
#ifndef HAVE_MPI
		    obj->flags & MSIEVE_FLAG_STOP_SIEVING ||
#endif
		    dim_solved >= next_dump)) {

			vv_mul_BxN_NxB(packed_matrix, v0, vnext, d, n);
			for (i = 0; i < VBITS; i++) {
				if (!v_is_all_zeros(d[i])) {
					logprintf(obj, "error: corrupt state, "
					       "please restart from "
					       "checkpoint\n");
					printf("\nerror: corrupt state, "
					       "please restart from "
					       "checkpoint\n");
#ifdef HAVE_MPI
					MPI_Abort(MPI_COMM_WORLD, 
							MPI_ERR_ASSERT);
#else
					exit(-1);
#endif
				}
			}
			/* check passed */
			next_check = ((dim_solved + 6 * VBITS) / 
					check_interval + 1) * check_interval;
			vv_copy(v0, vnext, n);
		}

		/* compute d, fold it into vnext and update v'*v0 */

		for (i = 0; i < VBITS; i++) {
			d[i] = v_xor(vt_a_v[0][i],
					v_and(vt_a2_v[0][i], mask0));
		}

		mul_BxB_BxB(winv[0], d, d);

		for (i = 0; i < VBITS; i++)
			d[i] = v_xor(d[i], bitmask[i]);

		vv_mul_NxB_BxB_acc(packed_matrix, v[0], d, vnext, n);

		transpose_BxB(d, d);
		mul_BxB_BxB(d, vt_v0[0], vt_v0_next);

		/* compute e, fold it into vnext and update v'*v0 */

		mul_BxB_BxB(winv[1], vt_a_v[0], e);

		for (i = 0; i < VBITS; i++)
			e[i] = v_and(e[i], mask0);

		vv_mul_NxB_BxB_acc(packed_matrix, v[1], e, vnext, n);

		transpose_BxB(e, e);
		mul_BxB_BxB(e, vt_v0[1], e);
		for (i = 0; i < VBITS; i++)
			vt_v0_next[i] = v_xor(vt_v0_next[i], e[i]);

		/* compute f, fold it in. Montgomery shows that 
		   this is unnecessary (f would be zero) if the 
		   previous value of v had full rank */
		if (!v_is_all_ones(mask1)) {
			mul_BxB_BxB(vt_a_v[1], winv[1], f);

			for (i = 0; i < VBITS; i++)
				f[i] = v_xor(f[i], bitmask[i]);

			mul_BxB_BxB(winv[2], f, f);

			for (i = 0; i < VBITS; i++) {
				f2[i] = v_and(mask0, 
					    v_xor(vt_a_v[1][i],
						  v_and(vt_a2_v[1][i], mask1)));
			}

			mul_BxB_BxB(f, f2, f);

			vv_mul_NxB_BxB_acc(packed_matrix, v[2], f, vnext, n);

			transpose_BxB(f, f);
			mul_BxB_BxB(f, vt_v0[2], f);
			for (i = 0; i < VBITS; i++)
				vt_v0_next[i] = v_xor(vt_v0_next[i], f[i]);
		}

		/* update the computed solution 'x' */

		mul_BxB_BxB(winv[0], vt_v0[0], d);
		vv_mul_NxB_BxB_acc(packed_matrix, v[0], d, x, n);

		/* rotate all the variables */

		tmp = v[2]; 
		v[2] = v[1]; 
		v[1] = v[0]; 
		v[0] = vnext; 
		vnext = tmp;
		
		tmp = winv[2]; 
		winv[2] = winv[1]; 
		winv[1] = winv[0]; 
		winv[0] = tmp;
		
		tmp = vt_v0[2]; 
		vt_v0[2] = vt_v0[1]; 
		vt_v0[1] = vt_v0[0]; 
		vt_v0[0] = vt_v0_next; 
		vt_v0_next = tmp;
		
		tmp = vt_a_v[1]; vt_a_v[1] = vt_a_v[0]; vt_a_v[0] = tmp;
		
		tmp = vt_a2_v[1]; vt_a2_v[1] = vt_a2_v[0]; vt_a2_v[0] = tmp;

		memcpy(s[1], s[0], VBITS * sizeof(uint32));
		mask1 = mask0;
		dim1 = dim0;

		MPI_NODE_0_START

		/* possibly print a status update */

		if (report_interval) {
			if (dim_solved >= next_report) {
				time_t curr_time = time(NULL);
				double elapsed = curr_time - first_time;
				uint32 eta = elapsed * (max_n - dim_solved) /
						(dim_solved - first_dim_solved);

				fprintf(stderr, "linear algebra completed %u "
					"of %u dimensions (%1.1f%%, ETA "
					"%dh%2dm)    \r",
					dim_solved, max_n, 100.0 * dim_solved / 
					max_n, eta / 3600, (eta % 3600) / 60);

				/* report the ETA to the logfile once 
				   (wait 6 intervals for a better ETA) */

				if (++log_eta_once == 6) {
					logprintf(obj, "linear algebra at "
						   "%1.1f%%, ETA %dh%2dm\n",
						100.0 * dim_solved / max_n,
						eta / 3600, 
						(eta % 3600) / 60);
				}
				next_report = dim_solved + report_interval;
				fflush(stderr);
			}
		}

		MPI_NODE_0_END

		/* possibly dump a checkpoint file, check for interrupt.

		   Note that MPI cannot reliably dump a checkpoint when
		   interrupted, because multiple MPI processes must 
		   participate in the dump process but there is no real
		   way to signal them all to enter the following at the
		   same time, short of a logical-or of all the obj->flags
		   fields on every iteration. That costs about 3% of the
		   total runtime, and grows with increasing grid size */

		if (dump_interval) {
			if (dim_solved >= next_dump &&
			    dump_interval == DEFAULT_DUMP_INTERVAL) {

				/* the dump interval is the initial one,
				   chosen to accumulate some timing information.
				   Now compute the real dump interval, 
				   calibrated to happen about once per hour.

				   For MPI, the root node computes the
				   dump interval used by everyone */

				MPI_NODE_0_START
				time_t curr_time = time(NULL);
				double elapsed = curr_time - first_time;

				dump_interval = (3600.0 / elapsed) *
					       (dim_solved - first_dim_solved); 
				dump_interval = MAX(dump_interval,
						   DEFAULT_DUMP_INTERVAL + 1);

				/* make the dump interval a multiple of
				   the check interval. If this is not done,
				   eventually we will perform a check and
				   then less than three iterations later
				   will get a dump, which performs another
				   check. The Lanczos recurrence only
				   guarantees that check vectors more than
				   three iterations back will be orthogonal
				   to the current x, so this will cause 
				   spurious failures */

				dump_interval += check_interval -
						dump_interval %
						check_interval;

				logprintf(obj, "checkpointing every %u "
					   "dimensions\n", dump_interval);
				MPI_NODE_0_END
#ifdef HAVE_MPI
				MPI_TRY(MPI_Bcast(&dump_interval, 1, 
						MPI_INT, 0, obj->mpi_la_grid))
#endif
				next_dump = ((dim_solved + 6 * VBITS) / 
						dump_interval + 1) * dump_interval;
				continue;
			}

			if (
#ifndef HAVE_MPI
			    obj->flags & MSIEVE_FLAG_STOP_SIEVING ||
#endif
			    dim_solved >= next_dump) {

				dump_lanczos_state(obj, packed_matrix, 
						   x, vt_v0, v, v0, 
						   vt_a_v, vt_a2_v, winv, 
						   n, max_n, dim_solved, 
						   iter, s, dim1, scratch);
				next_dump = ((dim_solved + 6 * VBITS) / dump_interval + 1) * 
							dump_interval;
			}
			if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)
				break;
		}

	}

	MPI_NODE_0_START
	if (report_interval)
		fprintf(stderr, "\n");
	MPI_NODE_0_END

	logprintf(obj, "lanczos halted after %u iterations (dim = %u)\n", 
					iter, dim_solved);

	/* free unneeded storage */

	vv_free(vnext);
	vv_free(v0);

	aligned_free(vt_a_v[0]);
	aligned_free(vt_a_v[1]);
	aligned_free(vt_a2_v[0]);
	aligned_free(vt_a2_v[1]);
	aligned_free(winv[0]);
	aligned_free(winv[1]);
	aligned_free(winv[2]);
	aligned_free(vt_v0_next);
	aligned_free(vt_v0[0]);
	aligned_free(vt_v0[1]);
	aligned_free(vt_v0[2]);

	MPI_NODE_0_START

	if (dim0 == 0 || (obj->flags & MSIEVE_FLAG_STOP_SIEVING)) {
		vv_free(x);
		vv_free(scratch);
		vv_free(v[0]);
		vv_free(v[1]);
		vv_free(v[2]);
		if (dim0 == 0)
			logprintf(obj, "linear algebra failed, aborting\n");
#ifdef HAVE_MPI
		/* MPI cannot shut down gracefully from an interrupt */
		MPI_Abort(MPI_COMM_WORLD, MPI_ERR_ASSERT);
#endif
		return NULL;
	}

	MPI_NODE_0_END


	/* convert the output of the iteration to an actual
	   collection of nullspace vectors. Begin by multiplying
	   the output from the iteration by B */

#ifdef HAVE_MPI
	/* pull the result vectors into rank 0 */

	vv_free(v[1]);
	vv_free(v[2]);

	mul_MxN_NxB(packed_matrix, x, scratch);

	out2 = gather_nrows(obj, packed_matrix, scratch, NULL);
	out0 = gather_ncols(obj, packed_matrix, x, scratch, NULL);
	vv_free(x);
    
	mul_MxN_NxB(packed_matrix, v[0], scratch);

	out3 = gather_nrows(obj, packed_matrix, scratch, NULL);
	out1 = gather_ncols(obj, packed_matrix, v[0], scratch, NULL);
	vv_free(v[0]);
    
	vv_free(scratch);
#else
	mul_MxN_NxB(packed_matrix, x, v[1]);
	mul_MxN_NxB(packed_matrix, v[0], v[2]);
	vv_free(scratch);

	out0 = (v_t *)aligned_malloc(max_n * sizeof(v_t), 64);
	vv_copyout(out0, x, max_n);
	vv_free(x);

	out1 = (v_t *)aligned_malloc(max_n * sizeof(v_t), 64);
	vv_copyout(out1, v[0], max_n);
	vv_free(v[0]);

	out2 = (v_t *)aligned_malloc(max_n * sizeof(v_t), 64);
	vv_copyout(out2, v[1], max_n);
	vv_free(v[1]);

	out3 = (v_t *)aligned_malloc(max_n * sizeof(v_t), 64);
	vv_copyout(out3, v[2], max_n);
	vv_free(v[2]);
#endif

	MPI_NODE_0_START
        
	/* make sure the last few words of the above matrix products
	   are zero, since the postprocessing will be using them */

	for (i = packed_matrix->max_nrows; 
			i < packed_matrix->max_ncols; i++) {
		out2[i] = out3[i] = v_zero;
	}

	/* if necessary, add in the contribution of the
	   first few rows that were originally in B. We 
	   expect there to be about VBITS - POST_LANCZOS_ROWS 
	   bit vectors that are in the nullspace of B and
	   post_lanczos_matrix simultaneously */

	if (post_lanczos_matrix) {
		for (i = 0; i < POST_LANCZOS_ROWS; i++) {
			v_t accum0 = v_zero;
			v_t accum1 = v_zero;
			uint32 j;
			for (j = 0; j < max_n; j++) {
				if (v_bitset(post_lanczos_matrix[j], i)) {
					accum0 = v_xor(accum0, out0[j]);
					accum1 = v_xor(accum1, out1[j]);
				}
			}
			out2[i] = v_xor(out2[i], accum0);
			out3[i] = v_xor(out3[i], accum1);
		}
	}

	*num_deps_found = combine_cols(max_n, out0, out1, out2, out3);

	MPI_NODE_0_END

	aligned_free(out1);
	aligned_free(out2);
	aligned_free(out3);

	if (*num_deps_found == 0)
		logprintf(obj, "lanczos error: only trivial "
				"dependencies found\n");
	else
		logprintf(obj, "recovered %u nontrivial dependencies\n", 
				*num_deps_found);
	return out0;
}

/*-----------------------------------------------------------------------*/
uint64 * block_lanczos(msieve_obj *obj, 
			uint32 nrows, uint32 max_nrows, uint32 start_row,
			uint32 num_dense_rows, 
			uint32 ncols, uint32 max_ncols, uint32 start_col,
			la_col_t *B, uint32 *num_deps_found) {
	
	/* External interface to the linear algebra */

	v_t *post_lanczos_matrix = NULL;
	uint64 *dependencies = NULL;
	v_t *lanczos_output = NULL;
	packed_matrix_t packed_matrix;
	uint32 dump_interval;
	uint32 have_post_lanczos;
#ifdef HAVE_MPI
	uint32 start_sub;
#endif

	if (max_ncols <= max_nrows) {
		logprintf(obj, "matrix needs more columns than rows; "
                 "try adding 2-3%% more relations\n");
		// this is a recoverable issue when run from yafu.
		//exit(-1);
	}

	/* optionally remove the densest rows of the matrix, and
	   optionally pack a few more rows into dense format */

	have_post_lanczos = form_post_lanczos_matrix(obj, &nrows,
					&num_dense_rows, ncols, B, 
					&post_lanczos_matrix);
	if (num_dense_rows) {
		logprintf(obj, "matrix includes %u packed rows\n", 
					num_dense_rows);
	}

	memset(&packed_matrix, 0, sizeof(packed_matrix_t));

#ifndef HAVE_MPI
	if (have_post_lanczos)
		max_nrows -= POST_LANCZOS_ROWS;

#else
	/* tell all the MPI processes whether a post lanczos matrix
	   was constructed */

	MPI_TRY(MPI_Bcast(&have_post_lanczos, 1, MPI_INT, 0,
			obj->mpi_la_col_grid))

	if (have_post_lanczos) {
		/* adjust the number of rows to reflect the fact
		   that the matrix is now POST_LANCZOS_ROWS smaller.
		   Fortunately, MPI processes below the top row of
		   the grid have their row numbers relative to offset
		   start_row, so we don't have to adjust all the data */
		
		max_nrows -= POST_LANCZOS_ROWS;
		if (obj->mpi_la_row_rank > 0)
			start_row -= POST_LANCZOS_ROWS;
	}

	/* give the bounds for scatter-gather operations 
	   to the MPI root node */

	if (obj->mpi_la_row_rank == 0) {
		MPI_TRY(MPI_Gather(&ncols, 1, MPI_INT, 
				packed_matrix.col_counts, 
				1, MPI_INT, 0, obj->mpi_la_row_grid))
		MPI_TRY(MPI_Gather(&start_col, 1, MPI_INT, 
				packed_matrix.col_offsets, 
				1, MPI_INT, 0, obj->mpi_la_row_grid))
	}

	if (obj->mpi_la_col_rank == 0) {
		MPI_TRY(MPI_Gather(&nrows, 1, MPI_INT, 
				packed_matrix.row_counts, 
				1, MPI_INT, 0, obj->mpi_la_col_grid))
		MPI_TRY(MPI_Gather(&start_row, 1, MPI_INT, 
				packed_matrix.row_offsets, 
				1, MPI_INT, 0, obj->mpi_la_col_grid))
	}
    
	/* figure out the bounds of scatter-gather operations
	   down each MPI column */

	global_chunk_info(ncols, obj->mpi_nrows, obj->mpi_la_row_rank,
			 &packed_matrix.nsubcols, &start_sub);   
	
	MPI_TRY(MPI_Allgather(&packed_matrix.nsubcols, 1, MPI_INT, 
			packed_matrix.subcol_counts, 
			1, MPI_INT, obj->mpi_la_col_grid))
	MPI_TRY(MPI_Allgather(&start_sub, 1, MPI_INT, 
			packed_matrix.subcol_offsets, 
			1, MPI_INT, obj->mpi_la_col_grid))

#if VWORDS > 1
	/* scatter-gather operations count 64-bit words and not 
	   VBITS-bit vectors, so scale the counts */
	{
		uint32 i;
		for (i = 0; i < obj->mpi_nrows; i++) {
			packed_matrix.row_counts[i] *= VWORDS;
			packed_matrix.row_offsets[i] *= VWORDS;
			packed_matrix.subcol_counts[i] *= VWORDS;
			packed_matrix.subcol_offsets[i] *= VWORDS;
		}
		for (i = 0; i < obj->mpi_ncols; i++) {
			packed_matrix.col_counts[i] *= VWORDS;
			packed_matrix.col_offsets[i] *= VWORDS;
		}
	}
#endif

	/* if using a post-lanczos matrix, gather the matrix elements
	   at the root node since all of them will be necessary at once */

	if (post_lanczos_matrix != NULL && obj->mpi_la_row_rank == 0) {
		if (obj->mpi_la_col_rank == 0) {
			post_lanczos_matrix = xrealloc(post_lanczos_matrix,
						max_ncols * sizeof(v_t));
		}

		MPI_TRY(MPI_Gatherv((obj->mpi_la_col_rank == 0) ?
					MPI_IN_PLACE : post_lanczos_matrix, 
				VWORDS * ncols, MPI_LONG_LONG, 
				post_lanczos_matrix,
				packed_matrix.col_counts,
				packed_matrix.col_offsets,
				MPI_LONG_LONG, 0, obj->mpi_la_row_grid))

		if (obj->mpi_la_col_rank != 0) {
			free(post_lanczos_matrix);
			post_lanczos_matrix = NULL;
		}
	}
#endif
	if (have_post_lanczos)
		count_matrix_nonzero(obj, nrows, num_dense_rows, ncols, B);

	packed_matrix_init(obj, &packed_matrix, B, 
			   nrows, max_nrows, start_row,
			   ncols, max_ncols, start_col, 
			   num_dense_rows,
#ifdef HAVE_MPI
			   NUM_MEDIUM_ROWS / obj->mpi_nrows
#else
			   NUM_MEDIUM_ROWS
#endif
			   );

	/* set up for writing checkpoint files. This only applies
	   to the largest matrices. The initial dump interval is
	   just to establish timing information */

	dump_interval = 0;
	if (max_nrows > 1000000) {
		dump_interval = DEFAULT_DUMP_INTERVAL;
		obj->flags |= MSIEVE_FLAG_SIEVING_IN_PROGRESS;
	}

	/* solve the matrix */

	lanczos_output = block_lanczos_core(obj, &packed_matrix,
						num_deps_found,
						post_lanczos_matrix,
						dump_interval);

	if (dump_interval)
		obj->flags &= ~MSIEVE_FLAG_SIEVING_IN_PROGRESS;

	if (*num_deps_found) {
		uint32 i;

		dependencies = (uint64 *)xmalloc(max_ncols * sizeof(uint64));

		if (*num_deps_found > 64)
			logprintf(obj, "saving only 64 dependencies\n");

		for (i = 0; i < max_ncols; i++)
			dependencies[i] = lanczos_output[i].w[0];
	}

	/* note that the following frees any auxiliary packed
	   matrix structures, and also frees the column entries from
	   the input matrix (whether packed or not) */

	packed_matrix_free(&packed_matrix);
	aligned_free(lanczos_output);
	return dependencies;
}
