/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: lanczos_matmul.c 1025 2018-08-19 02:20:28Z jasonp_sf $
--------------------------------------------------------------------*/

#include "lanczos.h"

/*-------------------------------------------------------------------*/
void mul_unpacked(packed_matrix_t *matrix,
			  v_t *x, v_t *b) 
{
	uint32 ncols = matrix->ncols;
	uint32 num_dense_rows = matrix->num_dense_rows;
	la_col_t *A = matrix->unpacked_cols;
	uint32 i, j, k;

	memset(b, 0, ncols * sizeof(v_t));
	
	for (i = 0; i < ncols; i++) {
		la_col_t *col = A + i;
		uint32 *row_entries = col->data;
		v_t tmp = x[i];

		for (j = 0; j < col->weight; j++) {
			k = row_entries[j]; 
			b[k] = v_xor(b[k], tmp);
		}
	}

	if (num_dense_rows) {
		for (i = 0; i < ncols; i++) {
			la_col_t *col = A + i;
			uint32 *row_entries = col->data + col->weight;
			v_t tmp = x[i];
	
			for (j = 0; j < num_dense_rows; j++) {
				if (row_entries[j / 32] & 
						((uint32)1 << (j % 32))) {
					b[j] = v_xor(b[j], tmp);
				}
			}
		}
	}
}

/*-------------------------------------------------------------------*/
void mul_trans_unpacked(packed_matrix_t *matrix,
				v_t *x, v_t *b) 
{
	uint32 ncols = matrix->ncols;
	uint32 num_dense_rows = matrix->num_dense_rows;
	la_col_t *A = matrix->unpacked_cols;
	uint32 i, j;

	for (i = 0; i < ncols; i++) {
		la_col_t *col = A + i;
		uint32 *row_entries = col->data;
		v_t accum = v_zero;

		for (j = 0; j < col->weight; j++) {
			accum = v_xor(accum, x[row_entries[j]]);
		}
		b[i] = accum;
	}

	if (num_dense_rows) {
		for (i = 0; i < ncols; i++) {
			la_col_t *col = A + i;
			uint32 *row_entries = col->data + col->weight;
			v_t accum = b[i];
	
			for (j = 0; j < num_dense_rows; j++) {
				if (row_entries[j / 32] &
						((uint32)1 << (j % 32))) {
					accum = v_xor(accum, x[j]);
				}
			}
			b[i] = accum;
		}
	}
}

/*-------------------------------------------------------------------*/
void packed_matrix_init(msieve_obj *obj,
			packed_matrix_t *p, la_col_t *A,
			uint32 nrows, uint32 max_nrows, uint32 start_row, 
			uint32 ncols, uint32 max_ncols, uint32 start_col, 
			uint32 num_dense_rows, uint32 first_block_size) {

	/* initialize */

	p->unpacked_cols = A;
	p->nrows = nrows;
	p->max_nrows = max_nrows;
	p->start_row = start_row;
	p->ncols = ncols;
	p->max_ncols = max_ncols;
	p->start_col = start_col;
	p->num_dense_rows = num_dense_rows;
	p->num_threads = 1;
#ifdef HAVE_MPI
	p->mpi_size = obj->mpi_size;
	p->mpi_nrows = obj->mpi_nrows;
	p->mpi_ncols = obj->mpi_ncols;
	p->mpi_la_row_rank = obj->mpi_la_row_rank;
	p->mpi_la_col_rank = obj->mpi_la_col_rank;
	p->mpi_la_row_grid = obj->mpi_la_row_grid;
	p->mpi_la_col_grid = obj->mpi_la_col_grid;
#endif

	matrix_extra_init(obj, p, first_block_size);
}

/*-------------------------------------------------------------------*/
void packed_matrix_free(packed_matrix_t *p) {

	matrix_extra_free(p);
}

/*-------------------------------------------------------------------*/
void mul_MxN_NxB(packed_matrix_t *A, void *x, 
			void *scratch) {
    
	/* Multiply the vector x[] by the matrix A and put the 
	   result in scratch[]. The MPI version needs an extra
	   scratch array because MPI reduction operations really
	   want to be out-of-place */

#ifdef HAVE_MPI
	v_t *scratch2 = (v_t *)scratch + MAX(A->ncols, A->nrows);

	if (A->mpi_size <= 1) {
#endif
		mul_core(A, x, scratch);
#ifdef HAVE_MPI
		return;
	}
    
	/* make each MPI column gather its own part of x */
	
	global_allgather(x, scratch, A->ncols, A->mpi_nrows, 
			A->mpi_la_row_rank, A->mpi_la_col_grid);
		
	mul_core(A, scratch, scratch2);
	
	/* make each MPI row combine all of its vectors. The
	   matrix-vector product is redundantly stored in each
	   MPI column, but this routine is called very rarely
	   so it's not worth removing the redundancy */
	
	global_xor(scratch2, scratch, A->nrows, A->mpi_ncols,
			   A->mpi_la_col_rank, A->mpi_la_row_grid);

#endif
}

/*-------------------------------------------------------------------*/
void mul_sym_NxN_NxB(packed_matrix_t *A, void *x, 
			void *b, void *scratch) {

	/* Multiply x by A and write to scratch, then
	   multiply scratch by the transpose of A and
	   write to b. x may alias b, but the two must
	   be distinct from scratch */

#ifdef HAVE_MPI
	v_t *scratch2 = (v_t *)scratch + MAX(A->ncols, A->nrows);
        
	if (A->mpi_size <= 1) {
#endif
		mul_core(A, x, scratch);
		mul_trans_core(A, scratch, b);
#ifdef HAVE_MPI
		return;
	}
    
	/* make each MPI column gather its own part of x */
	 
	global_allgather(x, scratch, A->ncols, A->mpi_nrows, 
			A->mpi_la_row_rank, A->mpi_la_col_grid);
	
	mul_core(A, scratch, scratch2);
		
	/* make each MPI row combine its own part of A*x */
	
	global_xor(scratch2, scratch, A->nrows, A->mpi_ncols,
			   A->mpi_la_col_rank, A->mpi_la_row_grid);
		
	mul_trans_core(A, scratch, scratch2);
		
	/* make each MPI row combine and scatter its own part of A^T * A*x */
		
	global_xor_scatter(scratch2, b, scratch,  A->ncols, A->mpi_nrows, 
			A->mpi_la_row_rank, A->mpi_la_col_grid);
#endif
}
