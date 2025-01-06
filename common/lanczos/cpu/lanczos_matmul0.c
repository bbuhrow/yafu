/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: lanczos_matmul0.c 1014 2017-06-12 03:00:48Z jasonp_sf $
--------------------------------------------------------------------*/

#include "lanczos_cpu.h"

/*-------------------------------------------------------------------*/
static void mul_packed(packed_matrix_t *p, v_t *x, v_t *b) 
{
	uint32 i, j;
	task_control_t task = {NULL, NULL, NULL, NULL};
	cpudata_t *c = (cpudata_t *)p->extra;

	c->x = x;
	c->b = b;

	/* start accumulating the dense matrix multiply results;
	   each thread has scratch space for these, so we don't have
	   to wait for the tasks to finish */

	task.run = mul_packed_small_core;

	for (i = 0; i < p->num_threads - 1; i++) {
		task.data = c->tasks + i;
		threadpool_add_task(c->threadpool, &task, 1);
	}
	mul_packed_small_core(c->tasks + i, i);

	/* switch to the sparse blocks */

	task.run = mul_packed_core;

	for (i = 0; i < c->num_superblock_cols; i++) {

		la_task_t *t = c->tasks;
				
		for (j = 0; j < p->num_threads; j++)
			t[j].block_num = i;

		for (j = 0; j < p->num_threads - 1; j++) {
			task.data = t + j;
			threadpool_add_task(c->threadpool, &task, 1);
		}

		mul_packed_core(t + j, j);
		if (j) {
			threadpool_drain(c->threadpool, 1);
		}
	}

	/* xor the small vectors from each thread */

	vv_copy(b, c->thread_data[0].tmp_b, 
			MAX(c->first_block_size, VBITS * 
				((p->num_dense_rows + VBITS - 1) / VBITS)));

	for (i = 1; i < p->num_threads; i++)
		vv_xor(b, c->thread_data[i].tmp_b,
			MAX(c->first_block_size, VBITS * 
				((p->num_dense_rows + VBITS - 1) / VBITS)));

#if defined(GCC_ASM32A) && defined(HAS_MMX)
	ASM_G volatile ("emms");
#elif defined(MSC_ASM32A) && defined(HAS_MMX)
	ASM_M emms
#endif
}

/*-------------------------------------------------------------------*/
static void mul_trans_packed(packed_matrix_t *p, v_t *x, v_t *b) 
{
	uint32 i, j;
	task_control_t task = {NULL, NULL, NULL, NULL};
	cpudata_t *c = (cpudata_t *)p->extra;

	c->x = x;
	c->b = b;

	task.run = mul_trans_packed_core;

	for (i = 0; i < c->num_superblock_rows; i++) {

		la_task_t *t = c->tasks;
				
		for (j = 0; j < p->num_threads; j++)
			t[j].block_num = i;

		for (j = 0; j < p->num_threads - 1; j++) {
			task.data = t + j;
			threadpool_add_task(c->threadpool, &task, 1);
		}

		mul_trans_packed_core(t + j, j);
		if (j) {
			threadpool_drain(c->threadpool, 1);
		}
	}

	if (p->num_dense_rows) {
		/* add in the dense matrix multiply blocks; these don't 
		   use scratch space, but need all of b to accumulate 
		   results so we have to wait until all tasks finish */

		task.run = mul_trans_packed_small_core;

		for (i = 0; i < p->num_threads - 1; i++) {
			task.data = c->tasks + i;
			threadpool_add_task(c->threadpool, &task, 1);
		}

		mul_trans_packed_small_core(c->tasks + i, i);
		if (i) {
			threadpool_drain(c->threadpool, 1);
		}
	}

#if defined(GCC_ASM32A) && defined(HAS_MMX)
	ASM_G volatile ("emms");
#elif defined(MSC_ASM32A) && defined(HAS_MMX)
	ASM_M emms
#endif
}

/*--------------------------------------------------------------------*/
static void matrix_thread_init(void *data, int thread_num) {

	packed_matrix_t *p = (packed_matrix_t *)data;
	cpudata_t *c = (cpudata_t *)p->extra;
	thread_data_t *t = c->thread_data + thread_num;

	/* we use this scratch vector for both matrix multiplies
	   and vector-vector operations; it has to be large enough
	   to support both. Note that first_block_size is split across
	   MPI rows, so it is conceivable with enough MPI processes that
	   the MAX() is necessary */

	t->tmp_b = (v_t *)vv_alloc(MAX(c->first_block_size, VBITS *
			(1 + (p->num_dense_rows + VBITS - 1) / VBITS)), p->extra);
}

/*-------------------------------------------------------------------*/
static void matrix_thread_free(void *data, int thread_num) {

	packed_matrix_t *p = (packed_matrix_t *)data;
	cpudata_t *c = (cpudata_t *)p->extra;
	thread_data_t *t = c->thread_data + thread_num;

	vv_free(t->tmp_b);
}

/*-------------------------------------------------------------------*/
static int compare_row_off(const void *x, const void *y) {
	entry_idx_t *xx = (entry_idx_t *)x;
	entry_idx_t *yy = (entry_idx_t *)y;

	if (xx->row_off > yy->row_off)
		return 1;
	if (xx->row_off < yy->row_off)
		return -1;

	return (int)xx->col_off - (int)yy->col_off;
}

/*--------------------------------------------------------------------*/
static void pack_med_block(packed_block_t *b)
{
	uint32 j, k, m;
	uint16 *med_entries;
	entry_idx_t *e;

	/* convert the first block in the stripe to a somewhat-
	   compressed format. Entries in this first block are stored 
	   by row, and all rows are concatenated into a single 
	   16-bit array */

	e = b->d.entries;
	qsort(e, (size_t)b->num_entries, 
			sizeof(entry_idx_t), compare_row_off);
	for (j = k = 1; j < b->num_entries; j++) {
		if (e[j].row_off != e[j-1].row_off)
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

	med_entries = (uint16 *)xmalloc((b->num_entries + 
					2 * k + 8) * sizeof(uint16));
	j = k = 0;
	while (j < b->num_entries) {
		for (m = 0; j + m < b->num_entries; m++) {
			if (m > 0 && e[j+m].row_off != e[j+m-1].row_off)
				break;
			med_entries[k+m+2] = e[j+m].col_off;
		}
		med_entries[k] = e[j].row_off;
		med_entries[k+1] = m;
		j += m;
		k += m + 2;
	}
	med_entries[k] = med_entries[k+1] = 0;
	free(b->d.entries);
	b->d.med_entries = med_entries;
}

/*--------------------------------------------------------------------*/
static void pack_matrix_core(packed_matrix_t *p)
{
	uint32 i, j, k;
	la_col_t *A = p->unpacked_cols;
	cpudata_t *c = (cpudata_t *)p->extra;
	uint32 dense_row_blocks;
	packed_block_t *curr_stripe;

	uint32 ncols = p->ncols;
	uint32 block_size = c->block_size;
	uint32 num_block_rows = c->num_block_rows;
	uint32 num_block_cols = c->num_block_cols;
	uint32 first_block_size = c->first_block_size;

	/* pack the dense rows VBITS at a time */

	dense_row_blocks = (p->num_dense_rows + VBITS - 1) / VBITS;
	if (dense_row_blocks) {
		c->dense_blocks = (v_t **)xmalloc(dense_row_blocks *
						sizeof(v_t *));
		for (i = 0; i < dense_row_blocks; i++) {
			c->dense_blocks[i] = (v_t *)vv_alloc(ncols, p->extra);
		}

		for (i = 0; i < ncols; i++) {
			la_col_t *col = A + i;
			uint32 *src = col->data + col->weight;
			for (j = 0; j < dense_row_blocks; j++) {
				for (k = 0; k < VWORDS; k++) {
					uint32 t = j * VWORDS + k;
					c->dense_blocks[j][i].w[k] = 
						(uint64)src[2 * t + 1] << 32 |
						(uint64)src[2 * t];
				}
			}
		}
	}

	/* allocate blocks in row-major order; a 'stripe' is
	   a vertical column of blocks. The first block in each
	   column has first_block_size rows instead of block_size */

	c->blocks = curr_stripe = (packed_block_t *)xcalloc(
						(size_t)num_block_rows *
						        num_block_cols,
						sizeof(packed_block_t));

	/* we convert the sparse part of the matrix to packed
	   format one stripe at a time. This limits the worst-
	   case memory use of the packing process */

	for (i = 0; i < num_block_cols; i++, curr_stripe++) {

		uint32 curr_cols = MIN(block_size, ncols - i * block_size);
		packed_block_t *b;

		/* count the number of nonzero entries in each block */

		for (j = 0; j < curr_cols; j++) {
			la_col_t *col = A + i * block_size + j;

			for (k = 0, b = curr_stripe; k < col->weight; k++) {
				uint32 index = col->data[k];
				uint32 row = 0;

				if (index >= first_block_size) {
					row = 1 + (index - first_block_size) /
						block_size;
				}

				b[row * num_block_cols].num_entries++;
			}
		}

		/* concatenate the nonzero elements of the matrix
		   columns corresponding to this stripe.
		   
		   We technically can combine the previous pass through
		   the columns with this pass, but on some versions of
		   libc the number of reallocations causes an incredible
		   slowdown */

		for (j = 0, b = curr_stripe; j < num_block_rows; 
						j++, b += num_block_cols) {
			b->d.entries = (entry_idx_t *)xmalloc(
						b->num_entries *
						sizeof(entry_idx_t));
			b->num_entries = 0;
		}

		for (j = 0; j < curr_cols; j++) {
			la_col_t *col = A + i * block_size + j;

			for (k = 0; k < col->weight; k++) {
				entry_idx_t *e;
				uint32 index = col->data[k];
				uint32 row = 0;

				if (index >= first_block_size) {
					row = 1 + (index - first_block_size) /
						block_size;
					index -= first_block_size +
						(row - 1) * block_size;
				}

				b = curr_stripe + row * num_block_cols;
				e = b->d.entries + b->num_entries++;
				e->row_off = index;
				e->col_off = j;
			}

			free(col->data);
			col->data = NULL;
		}

		pack_med_block(curr_stripe);
	}

	p->unpacked_cols = NULL;
}

/*-------------------------------------------------------------------*/
void matrix_extra_init(msieve_obj *obj, packed_matrix_t *p,
			uint32 first_block_size) {

	uint32 i;
	uint32 num_threads;
	uint32 block_size;
	uint32 superblock_size;
	thread_control_t control;
	cpudata_t *c;

	p->extra = c = (cpudata_t *)xcalloc(1, sizeof(cpudata_t));

	/* determine the number of threads to use */

	num_threads = obj->num_threads;
	if (num_threads < 2 || p->max_nrows < MIN_NROWS_TO_THREAD)
		num_threads = 1;

	p->num_threads = num_threads = MIN(num_threads, MAX_THREADS);

	/* start the thread pool; note that even single-threaded
	   runs need these structures to be allocated */

	c->first_block_size = first_block_size;

	control.init = matrix_thread_init;
	control.shutdown = matrix_thread_free;
	control.data = p;

	if (num_threads > 1) {
		c->threadpool = threadpool_init(num_threads - 1, 
						200, &control);
	}
	matrix_thread_init(p, num_threads - 1);

	/* pre-generate the structures to drive the thread pool */

	c->tasks = (la_task_t *)xmalloc(sizeof(la_task_t) * num_threads);

	for (i = 0; i < num_threads; i++) {
		c->tasks[i].matrix = p;
		c->tasks[i].task_num = i;
	}

	if (p->max_nrows <= MIN_NROWS_TO_PACK)
		return;

	/* determine the block sizes. We assume that the largest
	   cache in the system is unified and shared across all
	   threads. When performing matrix multiplies 'A*x=b', 
	   we choose the block size small enough so that one block
	   of b fits in L1 cache, and choose the superblock size
	   to be small enough so that a superblock's worth of x
	   or b takes up 3/4 of the largest cache in the system.
	   
	   Making the block size too small increases memory use
	   and puts more pressure on the larger caches, while
	   making the superblock too small reduces the effectiveness
	   of L1 cache and increases the synchronization overhead
	   in multithreaded runs */

	block_size = 8192;
	superblock_size = 3 * obj->cache_size2 / (4 * sizeof(v_t));

	/* possibly override from the command line */

	if (obj->nfs_args != NULL) {

		const char *tmp;

		tmp = strstr(obj->nfs_args, "la_block=");
		if (tmp != NULL)
			block_size = atoi(tmp + 9);

		tmp = strstr(obj->nfs_args, "la_superblock=");
		if (tmp != NULL)
			superblock_size = atoi(tmp + 14);
	}

	logprintf(obj, "using block size %u and superblock size %u for "
			"processor cache size %u kB\n", 
				block_size, superblock_size,
				obj->cache_size2 / 1024);

	c->block_size = block_size;
	c->num_block_cols = (p->ncols + block_size - 1) / block_size;
	c->num_block_rows = 1 + (p->nrows - first_block_size + 
				block_size - 1) / block_size;

	c->superblock_size = (superblock_size + block_size - 1) / block_size;
	c->num_superblock_cols = (c->num_block_cols + c->superblock_size - 1) / 
					c->superblock_size;
	c->num_superblock_rows = (c->num_block_rows - 1 + 
				c->superblock_size - 1) / c->superblock_size;

	/* do the core work of packing the matrix */

	pack_matrix_core(p);
}

/*-------------------------------------------------------------------*/
void matrix_extra_free(packed_matrix_t *p) {

	uint32 i;
	cpudata_t *c = (cpudata_t *)p->extra;

	if (p->unpacked_cols) {
		la_col_t *A = p->unpacked_cols;
		for (i = 0; i < p->ncols; i++) {
			free(A[i].data);
			A[i].data = NULL;
		}
	}
	else {
		for (i = 0; i < (p->num_dense_rows + VBITS - 1) / VBITS; i++)
			vv_free(c->dense_blocks[i]);

		for (i = 0; i < c->num_block_rows * c->num_block_cols; i++) 
			free(c->blocks[i].d.entries);

		free(c->dense_blocks);
		free(c->blocks);
	}

	if (p->num_threads > 1) {
		threadpool_drain(c->threadpool, 1);
		threadpool_free(c->threadpool);
	}
	matrix_thread_free(p, p->num_threads - 1);

	free(c->tasks);
	free(c);
}

/*-------------------------------------------------------------------*/
size_t packed_matrix_sizeof(packed_matrix_t *p) {

	uint32 i, j;
	size_t mem_use;

	/* account for the vectors used in the lanczos iteration */

#ifdef HAVE_MPI
	mem_use = (6 * p->nsubcols + 2 * 
			MAX(p->nrows, p->ncols)) * sizeof(v_t);
#else
	mem_use = 7 * p->max_ncols * sizeof(v_t);
#endif

	/* and for the matrix */

	if (p->unpacked_cols) {
		la_col_t *A = p->unpacked_cols;
		mem_use += p->ncols * (sizeof(la_col_t) +
				(p->num_dense_rows + 31) / 32);
		for (i = 0; i < p->ncols; i++) {
			mem_use += A[i].weight * sizeof(uint32);
		}
	}
	else {
		cpudata_t *c = (cpudata_t *)p->extra;

		uint32 num_blocks = c->num_block_rows * 
					c->num_block_cols;

		mem_use += sizeof(v_t) * p->num_threads * 
				c->first_block_size;

		mem_use += sizeof(packed_block_t) * num_blocks;

		mem_use += p->ncols * sizeof(v_t) *
				((p->num_dense_rows + VBITS - 1) / VBITS);

		for (j = 0; j < num_blocks; j++) {
			packed_block_t *b = c->blocks + j;

			if (j < c->num_block_cols) {
				mem_use += (b->num_entries + 
					    2 * c->first_block_size) * 
						sizeof(uint16);
			}
			else {
				mem_use += b->num_entries *
						sizeof(entry_idx_t);
			}
		}
	}
	return mem_use;
}

/*-------------------------------------------------------------------*/
void mul_core(packed_matrix_t *A, void *x_in, void *b_in) {
    
	/* Multiply the vector x[] by the matrix A and put the 
	   result in b[]. x must not alias b */

	v_t *x = (v_t *)x_in;
	v_t *b = (v_t *)b_in;

	if (A->unpacked_cols)
		mul_unpacked(A, x, b);
	else
		mul_packed(A, x, b);
}

/*-------------------------------------------------------------------*/
void mul_trans_core(packed_matrix_t *A, void *x_in, void *b_in) {
    
	/* Multiply the vector x[] by the transpose of matrix A 
	   and put the result in b[]. x must not alias b */

	v_t *x = (v_t *)x_in;
	v_t *b = (v_t *)b_in;

	if (A->unpacked_cols)
		mul_trans_unpacked(A, x, b);
	else
		mul_trans_packed(A, x, b);
}
