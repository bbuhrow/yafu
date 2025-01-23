/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id$
--------------------------------------------------------------------*/

#ifndef _COMMON_LANCZOS_LANCZOS_H_
#define _COMMON_LANCZOS_LANCZOS_H_

#include <ms_common.h>

#ifdef __cplusplus
extern "C" {
#endif

/* the number of dependencies that will be generated
   internally (a max of 64 dependencies will be exposed
   to calling code) */

#ifndef VBITS
#error "linear algebra vector length not specified"
#endif

#define VWORDS ((VBITS + 63) / 64)

#if VBITS!=64 && VBITS!=128 && VBITS!=192 && VBITS!=256 && VBITS!=320 && VBITS!=384 && VBITS!=448 && VBITS!=512
#error "unsupported vector size"
#endif

typedef struct {
	uint64 w[VWORDS];
} v_t;

static INLINE v_t v_and(v_t a, v_t b) {
	v_t res;
	int i;
	for (i = 0; i < VWORDS; i++) res.w[i] = a.w[i] & b.w[i];
	return res;
}

static INLINE v_t v_or(v_t a, v_t b) {
	v_t res;
	int i;
	for (i = 0; i < VWORDS; i++) res.w[i] = a.w[i] | b.w[i];
	return res;
}

static INLINE v_t v_xor(v_t a, v_t b) {
	v_t res;
	int i;
	for (i = 0; i < VWORDS; i++) res.w[i] = a.w[i] ^ b.w[i];
	return res;
}

static INLINE uint32 v_bitset(v_t a, uint32 bit) {
	if (a.w[bit >> 6] & ((uint64)1 << (bit & 63)))
		return 1;
	return 0;
}

static INLINE uint32 v_is_all_zeros(v_t a) {

	uint32 i;

	for (i = 0; i < VWORDS; i++)
		if (a.w[i])
			return 0;
	return 1;
}

static INLINE uint32 v_is_all_ones(v_t a) {

	uint32 i;

	for (i = 0; i < VWORDS; i++)
		if (a.w[i] != (uint64)(-1))
			return 0;
	return 1;
}

static INLINE v_t v_random(uint32 *seed1, uint32 *seed2) {
	v_t res;
	uint32 i;

	for (i = 0; i < VWORDS; i++)
		res.w[i] = (uint64)(get_rand(seed1, seed2)) << 32 |
			   (uint64)(get_rand(seed1, seed2));

	return res;
}

static const v_t v_zero = {{0}};

/* for matrices of dimension exceeding MIN_POST_LANCZOS_DIM,
   the first POST_LANCZOS_ROWS rows are handled in a separate
   Gauss elimination phase after the Lanczos iteration
   completes. This means the lanczos code will produce about
   VBITS - POST_LANCZOS_ROWS dependencies on average. 
   
   The code will still work if POST_LANCZOS_ROWS is 0, but I 
   don't know why you would want to do that. The first rows are 
   essentially completely dense, and removing them from the main 
   Lanczos iteration greatly reduces the amount of arithmetic 
   in a matrix multiply, as well as the memory footprint of 
   the matrix */

#define POST_LANCZOS_ROWS (VBITS - 16)
#define MIN_POST_LANCZOS_DIM 10000

/* the smallest matrix size that will be converted 
   to packed format */

#define MIN_NROWS_TO_PACK 30000

/* routines for cache-efficient multiplication of
   sparse matrices */

/* the number of moderately dense rows that are
   packed less tightly */

#define NUM_MEDIUM_ROWS 3000

/* struct representing a packed matrix */

typedef struct packed_matrix_t {
	uint32 nrows;
	uint32 max_nrows;
	uint32 start_row;

	uint32 ncols;
	uint32 max_ncols;
	uint32 start_col;

	uint32 num_dense_rows;
	uint32 num_threads;

	uint32 block_nnz; /* used by the CUDA code */

	la_col_t *unpacked_cols;  /* used if no packing takes place */

	void * extra; /* implementation-specific stuff */

#ifdef HAVE_MPI
	uint32 mpi_size;
	uint32 mpi_nrows;
	uint32 mpi_ncols;
	uint32 mpi_la_row_rank;
	uint32 mpi_la_col_rank;
	MPI_Datatype mpi_word;
	MPI_Comm mpi_la_grid;
	MPI_Comm mpi_la_row_grid;
	MPI_Comm mpi_la_col_grid;

	uint32 nsubcols;
	int32 subcol_counts[MAX_MPI_GRID_DIM];
	int32 subcol_offsets[MAX_MPI_GRID_DIM];    

	/* needed on root node only */
	int32 col_counts[MAX_MPI_GRID_DIM];
	int32 col_offsets[MAX_MPI_GRID_DIM]; 
	int32 row_counts[MAX_MPI_GRID_DIM];
	int32 row_offsets[MAX_MPI_GRID_DIM];

#endif

} packed_matrix_t;

void packed_matrix_init(msieve_obj *obj, 
			packed_matrix_t *packed_matrix,
			la_col_t *A, 
			uint32 nrows, uint32 max_nrows, uint32 start_row, 
			uint32 ncols, uint32 max_ncols, uint32 start_col, 
			uint32 num_dense_rows, uint32 first_block_size);

void packed_matrix_free(packed_matrix_t *packed_matrix);

size_t packed_matrix_sizeof(packed_matrix_t *packed_matrix);

void matrix_extra_init(msieve_obj *obj, 
			packed_matrix_t *packed_matrix,
			uint32 first_block_size);

void matrix_extra_free(packed_matrix_t *packed_matrix);

/* top-level calls for matrix multiplies */

void mul_MxN_NxB(packed_matrix_t *A, 
			void *x, void *scratch, void *scratch2);

void mul_sym_NxN_NxB(packed_matrix_t *A, void *x, 
			void *b, void *scratch, void *scratch2);

/* easy base cases for small problems */

void mul_unpacked(packed_matrix_t *matrix, v_t *x, v_t *b); 

void mul_trans_unpacked(packed_matrix_t *matrix, v_t *x, v_t *b);

/* implementation-specific matrix-vector product */

void mul_core(packed_matrix_t *A, void *x, void *b);
void mul_trans_core(packed_matrix_t *A, void *x, void *b);

#ifdef HAVE_MPI
void global_xor(void *send_buf, void *recv_buf, 
		uint32 bufsize, uint32 mpi_nodes, 
		uint32 mpi_rank, MPI_Datatype, MPI_Comm comm);

void global_chunk_info(uint32 total_size, uint32 num_nodes, 
		uint32 my_id, uint32 *chunk_size, uint32 *chunk_start);

void global_allgather(void *send_buf, void *recv_buf, 
                        uint32 bufsize, uint32 mpi_nodes, 
                        uint32 mpi_rank, MPI_Datatype, MPI_Comm comm);

void global_xor_scatter(void *send_buf, void *recv_buf, 
			void *scratch, uint32 bufsize, 
			uint32 mpi_nodes, uint32 mpi_rank, 
			MPI_Datatype, MPI_Comm comm);
			
v_t * gather_ncols(msieve_obj *obj,
			packed_matrix_t *packed_matrix,
			void *v_in, void *scratch_in, 
			v_t *out);
			
v_t * gather_nrows(msieve_obj *obj,
			packed_matrix_t *packed_matrix,
			void *scratch_in, v_t *out);
			
void scatter_ncols(msieve_obj *obj,
			packed_matrix_t *packed_matrix,
			void *out_in, void *scratch_in, 
			v_t *in);
#endif

/* top-level calls for vector-vector operations */

void *vv_alloc(uint32 n, void *extra);
void vv_free(void *v);
void vv_copyin(void *dest, v_t *src, uint32 n);
void vv_copy(void *dest, void *src, uint32 n);
void vv_copyout(v_t *dest, void *src, uint32 n);
void vv_clear(void *v, uint32 n);
void vv_xor(void *dest, void *src, uint32 n);
void vv_mask(void *v, v_t mask, uint32 n);

void vv_mul_NxB_BxB_acc(packed_matrix_t *A, void *v, v_t *x, 
			void *y, uint32 n);

void vv_mul_BxN_NxB(packed_matrix_t *A, void *x, void *y, 
			v_t *xy, uint32 n);

#ifdef __cplusplus
}
#endif

#endif /* !_COMMON_LANCZOS_LANCZOS_H_ */
