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

#ifndef _COMMON_LANCZOS_CPU_LANCZOS_CPU_H_
#define _COMMON_LANCZOS_CPU_LANCZOS_CPU_H_

#ifdef A64FX
#include <arm_sve.h>
#endif

#include "../lanczos.h"

#ifdef __cplusplus
extern "C" {
#endif

/* structure representing a nonzero element of
   the matrix after packing into block format. 
   The two fields are the row and column offsets
   from the top left corner of the block */

typedef struct {
	uint16 row_off;
	uint16 col_off;
} entry_idx_t;

/* struct representing one block */

typedef struct {
	uint32 num_entries;       /* number of nonzero matrix entries */
	union {
		entry_idx_t *entries;     /* nonzero entries */
		uint16 *med_entries;	  /* nonzero entries for medium dense rows */
	} d;
} packed_block_t;

#define MIN_NROWS_TO_THREAD 200000

/* implementation-specific structure */

typedef struct {

	/* used for block matrix multiplies */

	uint32 block_size;
	uint32 num_block_rows;
	uint32 num_block_cols;

	uint32 superblock_size;  /* in units of blocks */
	uint32 num_superblock_rows;
	uint32 num_superblock_cols;

	uint32 first_block_size;/* block size for the smallest row numbers */

	v_t **dense_blocks;  /* for holding dense matrix rows; 
				   dense_blocks[i] holds the i_th batch of
				   64 matrix rows */
	packed_block_t *blocks; /* sparse part of matrix, in block format */

} cpudata_t;

/* for big jobs, we use a multithreaded framework that calls
   these routines for the heavy lifting */

void mul_one_med_block(packed_block_t *curr_block,
			v_t *curr_col, v_t *curr_b);

void mul_one_block(const packed_block_t *curr_block,
			const v_t *curr_col, v_t * __restrict__ curr_b);

void mul_trans_one_med_block(packed_block_t *curr_block,
			v_t *curr_row, v_t *curr_b);

void mul_trans_one_block(const packed_block_t *curr_block,
				const v_t *curr_row, v_t * __restrict__ curr_b);

/* internal stuff for vector-vector operations within the
   matrix multiply */

void mul_NxB_BxB_acc(v_t *v, v_t *x, v_t *y, uint32 n);

void mul_BxN_NxB(v_t *x, v_t *y, v_t *xy, uint32 n);

void mul_BxN_NxB_2(v_t *x, v_t *y, v_t *xy, uint32 n);

#pragma omp declare reduction(^ : v_t : \
        omp_out = v_xor(omp_out, omp_in)) \
        initializer( omp_priv = {0})

#ifdef __cplusplus
}
#endif

#endif /* !_COMMON_LANCZOS_CPU_LANCZOS_CPU_H_ */
