/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	
       				   --jasonp@boo.net 9/24/08

Modified:	Ben Buhrow
Date:		11/24/09
Purpose:	Port into Yafu-1.14.
--------------------------------------------------------------------*/

#ifndef _LANCZOS_H_
#define _LANCZOS_H_

#include "qs_impl.h"
#include "factor.h"

#if defined(__GNUC__)
#include <pthread.h>
#endif

#define HAS_MMX
#define HAS_CMOV

#ifdef __cplusplus
extern "C" {
#endif

/* routines for cache-efficient multiplication of
   sparse matrices */

/* the smallest number of columns that will be
   converted to packed format */

#define QS_MIN_NCOLS_TO_PACK 18000

/* the number of moderately dense rows that are
   packed less tightly */

#define QS_NUM_MEDIUM_ROWS 3000

/* structure representing a nonzero element of
   the matrix after packing into block format. 
   The two fields are the row and column offsets
   from the top left corner of the block */

typedef struct {
	uint16_t row_off;
	uint16_t col_off;
} qs_entry_idx_t;

/* struct representing one block */

typedef struct {
	uint32_t start_row;
	uint32_t start_col;         /* coordinates of top left corner */
	uint32_t num_rows;
	uint32_t num_entries;       /* number of nonzero matrix entries */
	uint32_t num_entries_alloc; /* nonzero matrix entries allocated */
	qs_entry_idx_t *entries;     /* nonzero entries */
	uint16_t *med_entries;      /* nonzero entries for medium dense rows */
} qs_packed_block_t;

enum thread_command {
	COMMAND_INIT,
	COMMAND_WAIT,
	COMMAND_RUN,
	COMMAND_RUN_TRANS,
	COMMAND_END
};

/* struct used by threads for computing partial
   matrix multiplies */

typedef struct {
	/* items used during initialization */

	uint32_t my_oid;		/* number assigned to this thread */
	qs_la_col_t *initial_cols; /* unpacked matrix columns */
	uint32_t col_min;
	uint32_t col_max;		/* range of column indices to handle */
	uint32_t nrows_in;	/* number of rows in the matrix */
	uint32_t ncols_in;	/* number of columns in the matrix */
	uint32_t block_size;	/* used to pack the column entries */

	/* items used during matrix multiplies */

	uint32_t ncols;		/* number of columns used by this thread */
	uint32_t num_dense_rows;  /* number of rows packed by dense_blocks */
	uint64_t **dense_blocks;  /* for holding dense matrix rows; 
				   dense_blocks[i] holds the i_th batch of
				   64 matrix rows */
	uint32_t num_blocks;
	uint64_t *x;
	uint64_t *b;
	qs_packed_block_t *blocks; /* sparse part of matrix, in block format */

	/* fields for thread pool synchronization */

	volatile enum thread_command command;

#if defined(WIN32) || defined(_WIN64)
	HANDLE thread_id;
	HANDLE run_event;
	HANDLE finish_event;
#else
	pthread_t thread_id;
	pthread_mutex_t run_lock;
	pthread_cond_t run_cond;
#endif

} qs_msieve_thread_data_t;

#define QS_MAX_THREADS 64
#define QS_MIN_NCOLS_TO_THREAD 200000

/* struct representing a packed matrix */

typedef struct {
	uint32_t nrows;
	uint32_t ncols;
	uint32_t num_dense_rows;
	uint32_t num_threads;

	qs_la_col_t *unpacked_cols;  /* used if no packing takes place */

	qs_msieve_thread_data_t thread_data[QS_MAX_THREADS];

} qs_packed_matrix_t;

void yafu_packed_matrix_init(fact_obj_t *obj, 
			qs_packed_matrix_t *packed_matrix,
			qs_la_col_t *A, uint32_t nrows, uint32_t ncols,
			uint32_t num_dense_rows);

void yafu_packed_matrix_free(qs_packed_matrix_t *packed_matrix);

size_t yafu_packed_matrix_sizeof(qs_packed_matrix_t *packed_matrix);

void yafu_mul_MxN_Nx64(qs_packed_matrix_t *A, uint64_t *x, uint64_t *b);

void yafu_mul_trans_MxN_Nx64(qs_packed_matrix_t *A, uint64_t *x, uint64_t *b);

void yafu_mul_Nx64_64x64_acc(uint64_t *v, uint64_t *x, uint64_t *y, uint32_t n);

void yafu_mul_64xN_Nx64(uint64_t *x, uint64_t *y, uint64_t *xy, uint32_t n);

/* for big jobs, we use a multithreaded framework that calls
   these two routines for the heavy lifting */

void yafu_mul_packed_core(qs_msieve_thread_data_t *t);

void yafu_mul_trans_packed_core(qs_msieve_thread_data_t *t);

#ifdef __cplusplus
}
#endif

#endif /* !_LANCZOS_H_ */

