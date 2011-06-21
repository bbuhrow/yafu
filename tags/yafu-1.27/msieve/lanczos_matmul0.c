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

#include "lanczos.h"
#include "util.h"

/*-------------------------------------------------------------------*/
static void yafu_mul_unpacked(qs_packed_matrix_t *matrix,
			  uint64 *x, uint64 *b) {

	uint32 ncols = matrix->ncols;
	uint32 num_dense_rows = matrix->num_dense_rows;
	qs_la_col_t *A = matrix->unpacked_cols;
	uint32 i, j;
	
	memset(b, 0, ncols * sizeof(uint64));
	
	for (i = 0; i < ncols; i++) {
		qs_la_col_t *col = A + i;
		uint32 *row_entries = col->data;
		uint64 tmp = x[i];

		for (j = 0; j < col->weight; j++) {
			b[row_entries[j]] ^= tmp;
		}
	}

	if (num_dense_rows) {
		for (i = 0; i < ncols; i++) {
			qs_la_col_t *col = A + i;
			uint32 *row_entries = col->data + col->weight;
			uint64 tmp = x[i];
	
			for (j = 0; j < num_dense_rows; j++) {
				if (row_entries[j / 32] & 
						((uint32)1 << (j % 32))) {
					b[j] ^= tmp;
				}
			}
		}
	}
}

/*-------------------------------------------------------------------*/
static void yafu_mul_trans_unpacked(qs_packed_matrix_t *matrix,
				uint64 *x, uint64 *b) {

	uint32 ncols = matrix->ncols;
	uint32 num_dense_rows = matrix->num_dense_rows;
	qs_la_col_t *A = matrix->unpacked_cols;
	uint32 i, j;
	
	for (i = 0; i < ncols; i++) {
		qs_la_col_t *col = A + i;
		uint32 *row_entries = col->data;
		uint64 accum = 0;

		for (j = 0; j < col->weight; j++) {
			accum ^= x[row_entries[j]];
		}
		b[i] = accum;
	}

	if (num_dense_rows) {
		for (i = 0; i < ncols; i++) {
			qs_la_col_t *col = A + i;
			uint32 *row_entries = col->data + col->weight;
			uint64 accum = b[i];
	
			for (j = 0; j < num_dense_rows; j++) {
				if (row_entries[j / 32] &
						((uint32)1 << (j % 32))) {
					accum ^= x[j];
				}
			}
			b[i] = accum;
		}
	}
}

/*-------------------------------------------------------------------*/
static void yafu_mul_packed(qs_packed_matrix_t *matrix, uint64 *x, uint64 *b) {

	uint32 i;
	uint32 ncols = matrix->ncols;
	
	for (i = 0; i < matrix->num_threads; i++) {
		qs_msieve_thread_data_t *t = matrix->thread_data + i;

		/* use each thread's scratch vector, except the
		   first thead, which has no scratch vector but
		   uses b instead */

		t->x = x;
		if (i == 0)
			t->b = b;
		memset(t->b, 0, ncols * sizeof(uint64));

		/* fire off each part of the matrix multiply
		   in a separate thread from the thread pool, 
		   except the last part. The current thread 
		   does the last partial multiply, and this 
		   saves one synchronize operation */

		if (i == matrix->num_threads - 1) {
			yafu_mul_packed_core(t);
		}
		else {
			t->command = COMMAND_RUN;
#if defined(WIN32) || defined(_WIN64)
			SetEvent(t->run_event);
#else
			pthread_cond_signal(&t->run_cond);
			pthread_mutex_unlock(&t->run_lock);
#endif
		}
	}

	/* wait for each thread to finish. All the scratch
	   vectors used by threads get xor-ed into the final b
	   vector */

	for (i = 0; i < matrix->num_threads; i++) {
		qs_msieve_thread_data_t *t = matrix->thread_data + i;

		if (i < matrix->num_threads - 1) {
#if defined(WIN32) || defined(_WIN64)
			WaitForSingleObject(t->finish_event, INFINITE);
#else
			pthread_mutex_lock(&t->run_lock);
			while (t->command != COMMAND_WAIT)
				pthread_cond_wait(&t->run_cond, &t->run_lock);
#endif
		}

		if (i > 0) {
			uint64 *curr_b = t->b;
			uint32 j;

			for (j = 0; j < ncols; j++)
				b[j] ^= curr_b[j];
		}
	}

#if (defined(__GNUC__) || defined(__ICL)) && \
	defined(__i386__) && defined(HAS_MMX)
	asm volatile ("emms");
#elif defined(_MSC_VER) && !defined(_WIN64) && \
	defined(NDEBUG) && defined(HAS_MMX)
	__asm emms
#endif
}

/*-------------------------------------------------------------------*/
void yafu_mul_trans_packed(qs_packed_matrix_t *matrix, uint64 *x, uint64 *b) {

	uint32 i;
	uint32 ncols = matrix->ncols;
	uint64 *tmp_b[QS_MAX_THREADS];

	memset(b, 0, ncols * sizeof(uint64));
	
	for (i = 0; i < matrix->num_threads; i++) {
		qs_msieve_thread_data_t *t = matrix->thread_data + i;

		/* separate threads fill up disjoint portions
		   of a single b vector, and do not need 
		   per-thread scratch space */

		tmp_b[i] = t->b;
		t->x = x;
		t->b = b;

		/* fire off each part of the matrix multiply
		   in a separate thread from the thread pool, 
		   except the last part. The current thread 
		   does the last partial multiply, and this 
		   saves one synchronize operation */

		if (i == matrix->num_threads - 1) {
			yafu_mul_trans_packed_core(t);
		}
		else {
			t->command = COMMAND_RUN_TRANS;
#if defined(WIN32) || defined(_WIN64)
			SetEvent(t->run_event);
#else
			pthread_cond_signal(&t->run_cond);
			pthread_mutex_unlock(&t->run_lock);
#endif
		}
	}

	/* wait for each thread to finish */

	for (i = 0; i < matrix->num_threads; i++) {
		qs_msieve_thread_data_t *t = matrix->thread_data + i;

		if (i < matrix->num_threads - 1) {
#if defined(WIN32) || defined(_WIN64)
			WaitForSingleObject(t->finish_event, INFINITE);
#else
			pthread_mutex_lock(&t->run_lock);
			while (t->command != COMMAND_WAIT)
				pthread_cond_wait(&t->run_cond, &t->run_lock);
#endif
		}
		t->b = tmp_b[i];
	}

#if (defined(__GNUC__) || defined(__ICL)) && \
	defined(__i386__) && defined(HAS_MMX)
	asm volatile ("emms");
#elif defined(_MSC_VER) && !defined(_WIN64) && \
	defined(NDEBUG) && defined(HAS_MMX)
	__asm emms
#endif
}

/*-------------------------------------------------------------------*/
int yafu_compare_row_off(const void *x, const void *y) {
	qs_entry_idx_t *xx = (qs_entry_idx_t *)x;
	qs_entry_idx_t *yy = (qs_entry_idx_t *)y;
	
	if (xx->row_off > yy->row_off)
		return 1;
	if (xx->row_off < yy->row_off)
		return -1;

	return (int)xx->col_off - (int)yy->col_off;
}

static void yafu_matrix_thread_init(qs_msieve_thread_data_t *t) {

	uint32 i, j, k, m;
	uint32 num_row_blocks;
	uint32 num_col_blocks;
	uint32 dense_row_blocks;
	qs_packed_block_t *curr_stripe;
	qs_entry_idx_t *e;

	qs_la_col_t *A = t->initial_cols;
	uint32 nrows = t->nrows_in;
	uint32 col_min = t->col_min;
	uint32 col_max = t->col_max;
	uint32 block_size = t->block_size;
	uint32 num_dense_rows = t->num_dense_rows;
	
	/* each thread needs scratch space to store
	   matrix products. The first thread doesn't need
	   scratch space, it's provided by calling code */

	if (t->my_oid > 0)
		t->b = (uint64 *)xmalloc(t->ncols_in * sizeof(uint64));

	/* pack the dense rows 64 at a time */

	t->ncols = col_max - col_min + 1;
	dense_row_blocks = (num_dense_rows + 63) / 64;
	if (dense_row_blocks) {
		t->dense_blocks = (uint64 **)xmalloc(dense_row_blocks *
						sizeof(uint64 *));
		for (i = 0; i < dense_row_blocks; i++) {
			t->dense_blocks[i] = (uint64 *)xmalloc(t->ncols *
							sizeof(uint64));
		}

		for (i = 0; i < t->ncols; i++) {
			qs_la_col_t *c = A + col_min + i;
			uint32 *src = c->data + c->weight;
			for (j = 0; j < dense_row_blocks; j++) {
				t->dense_blocks[j][i] = 
						(uint64)src[2 * j + 1] << 32 |
						(uint64)src[2 * j];
			}
		}
	}

	/* allocate blocks in row-major order; a 'stripe' is
	   a vertical column of blocks */

	num_row_blocks = (nrows - QS_NUM_MEDIUM_ROWS + 
				(block_size-1)) / block_size + 1;
	num_col_blocks = ((col_max - col_min + 1) + 
				(block_size-1)) / block_size;
	t->num_blocks = num_row_blocks * num_col_blocks;
	t->blocks = curr_stripe = (qs_packed_block_t *)xcalloc(
						(size_t)t->num_blocks,
						sizeof(qs_packed_block_t));

	/* we convert the sparse part of the matrix to packed
	   format one stripe at a time. This limits the worst-
	   case memory use of the packing process */

	for (i = 0; i < num_col_blocks; i++, curr_stripe++) {

		uint32 curr_cols = MIN(block_size, (col_max - col_min + 1) - 
						i * block_size);
		qs_packed_block_t *b;

		/* initialize the blocks in stripe i. The first
		   block has QS_NUM_MEDIUM_ROWS rows, and all the
		   rest have block_size rows */

		for (j = 0, b = curr_stripe; j < num_row_blocks; j++) {

			if (j == 0) {
				b->start_row = 0;
				b->num_rows = QS_NUM_MEDIUM_ROWS;
			}
			else {
				b->start_row = QS_NUM_MEDIUM_ROWS +
						(j - 1) * block_size;
				b->num_rows = block_size;
			}

			b->start_col = col_min + i * block_size;
			b += num_col_blocks;
		}

		/* count the number of nonzero entries in each block */

		for (j = 0; j < curr_cols; j++) {
			qs_la_col_t *c = A + col_min + i * block_size + j;

			for (k = 0, b = curr_stripe; k < c->weight; k++) {
				uint32 index = c->data[k];

				while (index >= b->start_row + b->num_rows)
					b += num_col_blocks;
				b->num_entries_alloc++;
			}
		}

		/* concatenate the nonzero elements of the matrix
		   columns corresponding to this stripe. Note that
		   the ordering of elements within a single block
		   is arbitrary, and it's possible that some non-
		   default ordering can enhance performance. However,
		   I've tried a few different ideas and they're much
		   slower than just storing the nonzeros in the order
		   in which they occur in the unpacked matrix 
		   
		   We technically can combine the previous pass through
		   the columns with this pass, but on some versions of
		   libc the number of reallocations causes an incredible
		   slowdown */

		for (j = 0, b = curr_stripe; j < num_row_blocks; j++) {
			b->entries = (qs_entry_idx_t *)xmalloc(
						b->num_entries_alloc *
						sizeof(qs_entry_idx_t));
			b += num_col_blocks;
		}

		for (j = 0; j < curr_cols; j++) {
			qs_la_col_t *c = A + col_min + i * block_size + j;

			for (k = 0, b = curr_stripe; k < c->weight; k++) {
				uint32 index = c->data[k];

				while (index >= b->start_row + b->num_rows)
					b += num_col_blocks;

				e = b->entries + b->num_entries++;
				e->row_off = index - b->start_row;
				e->col_off = j;
			}

			free(c->data);
			c->data = NULL;
		}

		/* now convert the first block in the stripe to
		   a somewhat-compressed format. Entries in this
		   first block are stored by row, and all rows
		   are concatenated into a single 16-bit array */

		b = curr_stripe;
		e = b->entries;
		qsort(e, (size_t)b->num_entries, 
				sizeof(qs_entry_idx_t), yafu_compare_row_off);
		for (j = k = 0; j < b->num_entries; j++) {
			if (j == 0 || e[j].row_off != e[j-1].row_off)
				k++;
		}

		/* we need a 16-bit word for each element and two more
		   16-bit words at the start of each of the k packed
		   arrays making up med_entries. The first extra word
		   gives the row number and the second gives the number
		   of entries in that row. We also need a few extra words 
		   at the array end because the multiply code uses a 
		   software pipeline and would fetch off the end of 
		   med_entries otherwise */

		b->med_entries = (uint16 *)xmalloc((b->num_entries + 
						2 * k + 8) * sizeof(uint16));
		j = k = 0;
		while (j < b->num_entries) {
			for (m = 0; j + m < b->num_entries; m++) {
				if (m > 0 && 
				    e[j+m].row_off != e[j+m-1].row_off)
					break;
				b->med_entries[k+m+2] = e[j+m].col_off;
			}
			b->med_entries[k] = e[j].row_off;
			b->med_entries[k+1] = m;
			j += m;
			k += m + 2;
		}
		b->med_entries[k] = b->med_entries[k+1] = 0;
		free(b->entries);
		b->entries = NULL;
	}
}

/*-------------------------------------------------------------------*/
static void yafu_matrix_thread_free(qs_msieve_thread_data_t *t) {

	uint32 i;
	
	for (i = 0; i < (t->num_dense_rows + 63) / 64; i++)
		free(t->dense_blocks[i]);
	free(t->dense_blocks);

	for (i = 0; i < t->num_blocks; i++) {
		free(t->blocks[i].entries);
		free(t->blocks[i].med_entries);
	}
	free(t->blocks);
	if (t->my_oid > 0)
		free(t->b);
}

/*-------------------------------------------------------------------*/
#if defined(WIN32) || defined(_WIN64)
static DWORD WINAPI yafu_worker_thread_main(LPVOID thread_data) {
#else
static void *yafu_worker_thread_main(void *thread_data) {
#endif
	qs_msieve_thread_data_t *t = (qs_msieve_thread_data_t *)thread_data;

	while(1) {

		/* wait forever for work to do */
#if defined(WIN32) || defined(_WIN64)
		WaitForSingleObject(t->run_event, INFINITE);
#else
		pthread_mutex_lock(&t->run_lock);
		while (t->command == COMMAND_WAIT) {
			pthread_cond_wait(&t->run_cond, &t->run_lock);
		}
#endif
		/* do work */
		
		if (t->command == COMMAND_RUN)
			yafu_mul_packed_core(t);
		else if (t->command == COMMAND_RUN_TRANS)
			yafu_mul_trans_packed_core(t);
		else if (t->command == COMMAND_INIT)
			yafu_matrix_thread_init(t);
		else if (t->command == COMMAND_END)
			break;

		/* signal completion */

		t->command = COMMAND_WAIT;
#if defined(WIN32) || defined(_WIN64)
		SetEvent(t->finish_event);
#else
		pthread_cond_signal(&t->run_cond);
		pthread_mutex_unlock(&t->run_lock);
#endif
	}

	yafu_matrix_thread_free(t);

#if defined(WIN32) || defined(_WIN64)
	return 0;
#else
	return NULL;
#endif
}

/*-------------------------------------------------------------------*/
static void yafu_start_worker_thread(qs_msieve_thread_data_t *t, 
				uint32 is_master_thread) {

	/* create a thread that will handle matrix multiplies 
	   for block k of the matrix. The last block does 
	   not get its own thread (the current thread handles it) */
					
	if (is_master_thread) {
		yafu_matrix_thread_init(t);
		return;
	}

	t->command = COMMAND_INIT;
#if defined(WIN32) || defined(_WIN64)
	t->run_event = CreateEvent(NULL, FALSE, TRUE, NULL);
	t->finish_event = CreateEvent(NULL, FALSE, FALSE, NULL);
	t->thread_id = CreateThread(NULL, 0, yafu_worker_thread_main, t, 0, NULL);

	WaitForSingleObject(t->finish_event, INFINITE); /* wait for ready */
#else
	pthread_mutex_init(&t->run_lock, NULL);
	pthread_cond_init(&t->run_cond, NULL);

	pthread_cond_signal(&t->run_cond);
	pthread_mutex_unlock(&t->run_lock);
	pthread_create(&t->thread_id, NULL, yafu_worker_thread_main, t);

	pthread_mutex_lock(&t->run_lock); /* wait for ready */
	while (t->command != COMMAND_WAIT)
		pthread_cond_wait(&t->run_cond, &t->run_lock);
#endif
}

/*-------------------------------------------------------------------*/
static void yafu_stop_worker_thread(qs_msieve_thread_data_t *t,
				uint32 is_master_thread)
{
	if (is_master_thread) {
		yafu_matrix_thread_free(t);
		return;
	}
	
	t->command = COMMAND_END;
#if defined(WIN32) || defined(_WIN64)
	SetEvent(t->run_event);
	WaitForSingleObject(t->thread_id, INFINITE);
	CloseHandle(t->thread_id);
	CloseHandle(t->run_event);
	CloseHandle(t->finish_event);
#else
	pthread_cond_signal(&t->run_cond);
	pthread_mutex_unlock(&t->run_lock);
	pthread_join(t->thread_id, NULL);
	pthread_cond_destroy(&t->run_cond);
	pthread_mutex_destroy(&t->run_lock);
#endif
}

/*-------------------------------------------------------------------*/
void yafu_packed_matrix_init(fact_obj_t *obj,
			qs_packed_matrix_t *p, qs_la_col_t *A,
			uint32 nrows, uint32 ncols,
			uint32 num_dense_rows) {

	uint32 i, j, k;
	uint32 block_size;
	uint32 num_threads;
	uint32 num_nonzero;
	uint32 num_nonzero_per_thread;

	/* initialize */
	
	memset(p, 0, sizeof(qs_packed_matrix_t));
	p->unpacked_cols = A;
	p->nrows = nrows;
	p->ncols = ncols;
	p->num_dense_rows = num_dense_rows;

	if (ncols <= QS_MIN_NCOLS_TO_PACK)
		return;

	p->unpacked_cols = NULL;

	/* decide on the block size. We want to split the
	   cache into thirds; one third for x, one third for
	   b and the last third for streaming access to 
	   each block. If the matrix is small, adjust the
	   block size so that the matrix is divided 2.5
	   ways in each dimension.
	  
	   Note that the following does not compensate for
	   having multiple threads. If each thread gets a 
	   separate processor, no compensation is needed.
	   This could conceivably cause problems if multiple 
	   cores share the same cache, but multicore processors 
	   typically have pretty big caches anyway */

	block_size = obj->cache_size2 / (3 * sizeof(uint64));
	block_size = MIN(block_size, ncols / 2.5);
	block_size = MIN(block_size, 65536);
	if (block_size == 0)
		block_size = 32768;

	//block_size = 32768;

	if (VFLAG > 0)
		printf("using block size %u for "
			"processor cache size %u kB\n", 
				block_size, obj->cache_size2 / 1024);
	logprint(obj->logfile, "using block size %u for "
			"processor cache size %u kB\n", 
				block_size, obj->cache_size2 / 1024);

	/* determine the number of threads to use */

	num_threads = obj->num_threads;
	if (num_threads < 2 || ncols < QS_MIN_NCOLS_TO_THREAD)
		num_threads = 1;
	p->num_threads = num_threads = MIN(num_threads, QS_MAX_THREADS);

	/* compute the number of nonzero elements in the submatrix
	   given to each thread; overestimate the number slightly
	   so that we don't have one thread with almost no columns */

	for (i = num_nonzero = 0; i < ncols; i++)
		num_nonzero += A[i].weight;
	num_nonzero_per_thread = num_nonzero / num_threads + 1000;

	/* divide the matrix into groups of columns, one group
	   per thread, and pack each group separately */

	for (i = j = k = num_nonzero = 0; i < ncols; i++) {

		num_nonzero += A[i].weight;

		if (i == ncols - 1 || num_nonzero >= num_nonzero_per_thread) {
			qs_msieve_thread_data_t *t = p->thread_data + k;

			t->my_oid = k++;
			t->initial_cols = A;
			t->col_min = j;
			t->col_max = i;
			t->nrows_in = nrows;
			t->ncols_in = ncols;
			t->block_size = block_size;
			t->num_dense_rows = num_dense_rows;
			j = i + 1;
			num_nonzero = 0;
		}
	}

	/* update the number of threads, in case it's 
	   not what we estimated */

	p->num_threads = k;

	/* activate the threads one at a time. The last is the
	   master thread (i.e. not a thread at all). Each thread
	   packs its own portion of the matrix, so that memory
	   allocation is local on NUMA architectures */

	for (i = 0; i < p->num_threads - 1; i++)
		yafu_start_worker_thread(p->thread_data + i, 0);

	yafu_start_worker_thread(p->thread_data + i, 1);
}

/*-------------------------------------------------------------------*/
void yafu_packed_matrix_free(qs_packed_matrix_t *p) {

	uint32 i;
	
	if (p->unpacked_cols) {
		qs_la_col_t *A = p->unpacked_cols;
		for (i = 0; i < p->ncols; i++) {
			free(A[i].data);
			A[i].data = NULL;
		}
	}
	else {
		/* stop the worker threads; each will free
		   its own memory */

		for (i = 0; i < p->num_threads - 1; i++)
			yafu_stop_worker_thread(p->thread_data + i, 0);

		yafu_stop_worker_thread(p->thread_data + i, 1);
	}
}

/*-------------------------------------------------------------------*/
size_t yafu_packed_matrix_sizeof(qs_packed_matrix_t *p) {

	uint32 i, j;
	size_t mem_use = 0;
	
	if (p->unpacked_cols) {
		qs_la_col_t *A = p->unpacked_cols;
		mem_use = p->ncols * sizeof(qs_la_col_t);
		for (i = 0; i < p->ncols; i++) {
			mem_use += A[i].weight * sizeof(uint32);
		}
	}
	else {
		for (i = 0; i < p->num_threads; i++) {
			qs_msieve_thread_data_t *t = p->thread_data + i;

			mem_use += p->ncols * sizeof(uint64) +
				   t->num_blocks * sizeof(qs_packed_block_t) +
				   t->ncols * sizeof(uint64) *
					((t->num_dense_rows + 63) / 64);

			for (j = 0; j < t->num_blocks; j++) {
				qs_packed_block_t *b = t->blocks + j;
				if (b->entries) {
					mem_use += b->num_entries * 
							sizeof(uint32);
				}
				else {
					mem_use += (b->num_entries + 
						    2 * QS_NUM_MEDIUM_ROWS) * 
							sizeof(uint16);
				}
			}
		}
	}
	return mem_use;
}

/*-------------------------------------------------------------------*/
void yafu_mul_MxN_Nx64(qs_packed_matrix_t *A, uint64 *x, uint64 *b) {

	/* Multiply the vector x[] by the matrix A (stored
	   columnwise) and put the result in b[]. */
	
	if (A->unpacked_cols)
		yafu_mul_unpacked(A, x, b);
	else
		yafu_mul_packed(A, x, b);
}

/*-------------------------------------------------------------------*/
void yafu_mul_trans_MxN_Nx64(qs_packed_matrix_t *A, uint64 *x, uint64 *b) {

	/* Multiply the vector x[] by the transpose of the
	   matrix A and put the result in b[]. Since A is stored
	   by columns, this is just a matrix-vector product */
	
	if (A->unpacked_cols)
		yafu_mul_trans_unpacked(A, x, b);
	else
		yafu_mul_trans_packed(A, x, b);
}
