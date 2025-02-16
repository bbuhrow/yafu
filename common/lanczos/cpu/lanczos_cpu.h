/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: lanczos_cpu.h 1012 2017-06-11 17:42:14Z jasonp_sf $
--------------------------------------------------------------------*/

#ifndef _COMMON_LANCZOS_CPU_LANCZOS_CPU_H_
#define _COMMON_LANCZOS_CPU_LANCZOS_CPU_H_

#include <thread.h>
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

#define MAX_THREADS 32
#define MIN_NROWS_TO_THREAD 60000

/* struct used by threads for computing partial
   matrix multiplies */

typedef struct {
	/* items for matrix-vector operations */

	v_t *tmp_b;

	/* items for vector-vector operations */

	v_t *x;
	v_t *b;
	v_t *y;
	uint32 vsize;

} thread_data_t;

typedef struct {
	struct packed_matrix_t *matrix;
	uint32 task_num;
	uint32 block_num;
} la_task_t;

/* implementation-specific structure */

typedef struct {

	v_t *x;
	v_t *b;

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

	/* threading stuff */

	struct threadpool *threadpool;
	thread_data_t thread_data[MAX_THREADS];
	la_task_t *tasks;
} cpudata_t;

/* for big jobs, we use a multithreaded framework that calls
   these routines for the heavy lifting */

void mul_packed_core(void *data, int thread_num);

void mul_packed_small_core(void *data, int thread_num);

void mul_trans_packed_core(void *data, int thread_num);

void mul_trans_packed_small_core(void *data, int thread_num);

/* internal stuff for vector-vector operations within the
   matrix multiply */

void mul_NxB_BxB_acc(v_t *v, v_t *x, v_t *y, uint32 n);

void mul_BxN_NxB(v_t *x, v_t *y, v_t *xy, uint32 n);

#ifdef __cplusplus
}
#endif

#endif /* !_COMMON_LANCZOS_CPU_LANCZOS_CPU_H_ */
