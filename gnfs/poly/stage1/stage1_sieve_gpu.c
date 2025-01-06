/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: stage1_sieve_gpu.c 1056 2024-06-09 13:04:11Z brgladman $
--------------------------------------------------------------------*/

#include <sort_engine.h> /* interface to GPU sorting library */
#include <stage1.h>
#include <stage1_core_gpu/stage1_core.h>

/* GPU collision search; this code looks for self-collisions
   among arithmetic progressions, by finding k1 and k2 such that
   for two arithmetic progressions r1+k*p1^2 and r2+k*p2^2 we
   have

      r1 + k1*p1^2 = r2 + k2*p2^2

   such that
      - p1 and p2 are coprime and < 2^32
      - the value where they coincide is of size smaller
        than a fixed bound

   This code uses a sort routine to find collisions across all the
   p1 and p2 in the set simultaneously. We further use a 'special-q'
   formulation where all the inputs to the sort routine are
   constrained to fall on a third arithmetic progression r3 + k*q^2
   for some k. We choose a given q and for each of its roots run the
   complete sort. This is analogous to lattice sieving across the
   interval.
   
   This allows us to choose q so that the sort problem is of
   reasonable size but the collisions found are still over the
   original, impractically large range. */

enum {
	GPU_TRANS_PP32_R32 = 0,
	GPU_TRANS_PP32_R64,
	GPU_TRANS_PP64_R64,
	GPU_FINAL_32,
	GPU_FINAL_64,
	NUM_GPU_FUNCTIONS /* must be last */
};

static const char * gpu_kernel_names[] = 
{
	"sieve_kernel_trans_pp32_r32",
	"sieve_kernel_trans_pp32_r64",
	"sieve_kernel_trans_pp64_r64",
	"sieve_kernel_final_32",
	"sieve_kernel_final_64",
};

static const gpu_arg_type_list_t gpu_kernel_args[] = 
{
	/* sieve_kernel_trans_pp{32|64}_r{32|64} */
	{ 12,
		{
		  GPU_ARG_PTR,
		  GPU_ARG_UINT32,
		  GPU_ARG_PTR,
		  GPU_ARG_UINT32,
		  GPU_ARG_PTR,
		  GPU_ARG_PTR,
		  GPU_ARG_PTR,
		  GPU_ARG_UINT32,
		  GPU_ARG_UINT32,
		  GPU_ARG_UINT32,
		  GPU_ARG_UINT32,
		  GPU_ARG_UINT32,
		}
	},
	/* sieve_kernel_final_{32|64} */
	{ 6,
		{
		  GPU_ARG_PTR,
		  GPU_ARG_PTR,
		  GPU_ARG_UINT32,
		  GPU_ARG_PTR,
		  GPU_ARG_PTR,
		  GPU_ARG_UINT32,
		}
	},
};

/*------------------------------------------------------------------------*/

typedef struct {
	uint32 num_roots;
	uint32 num_p;
	uint32 num_p_alloc;

	uint32 *p;
	uint64 *start_roots;

	CUdeviceptr dev_p;
	CUdeviceptr dev_start_roots;

	union { 
		uint64 *roots64;
		uint32 *roots32;
		void *roots;
	} r;

} p_soa_var_t;

#define MAX_P_SOA_ARRAYS 5

typedef struct {
	uint32 num_arrays;
	uint32 max_arrays;
	uint32 max_p_roots;
	uint32 pp_is_64;

	p_soa_var_t start_soa[MAX_P_SOA_ARRAYS];
	p_soa_var_t *soa[MAX_P_SOA_ARRAYS];
} p_soa_array_t;

static p_soa_array_t *
p_soa_array_init(uint32 degree)
{
	uint32 i;
	p_soa_array_t *s = (p_soa_array_t *)xcalloc(1, 
					sizeof(p_soa_array_t));
	p_soa_var_t *start_soa = s->start_soa;

	switch (degree) {
	case 4:
		s->max_arrays = 3;
		start_soa[0].num_roots = 2;
		start_soa[1].num_roots = 4;
		start_soa[2].num_roots = 8;
		break;

	case 5:
		s->max_arrays = 3;
		start_soa[0].num_roots = 1;
		start_soa[1].num_roots = 5;
		start_soa[2].num_roots = 25;
		break;

	case 6:
		s->max_arrays = 5;
		start_soa[0].num_roots = 2;
		start_soa[1].num_roots = 4;
		start_soa[2].num_roots = 6;
		start_soa[3].num_roots = 12;
		start_soa[4].num_roots = 36;
		break;

	case 7: /* ;) */
		s->max_arrays = 3;
		start_soa[0].num_roots = 1;
		start_soa[1].num_roots = 7;
		start_soa[2].num_roots = 49;
		break;
	}
	s->max_p_roots = start_soa[s->max_arrays - 1].num_roots;

	for (i = 0; i < s->max_arrays; i++) {
		p_soa_var_t *soa = s->start_soa + i;

		soa->num_p = 0;
		soa->num_p_alloc = 256;
		soa->p = (uint32 *)xmalloc(soa->num_p_alloc * sizeof(uint32));
		soa->start_roots = (uint64 *)xmalloc(soa->num_roots *
					soa->num_p_alloc * sizeof(uint64));
		soa->r.roots = xmalloc(soa->num_roots *
					soa->num_p_alloc * sizeof(uint64));
		CUDA_TRY(cuMemAlloc(&soa->dev_p, 
					soa->num_p_alloc * 
					sizeof(uint32)))
		CUDA_TRY(cuMemAlloc(&soa->dev_start_roots, 
					soa->num_p_alloc *
					soa->num_roots * sizeof(uint64)))
	}

	return s;
}

static void
p_soa_array_free(p_soa_array_t *s)
{
	uint32 i;

	for (i = 0; i < s->max_arrays; i++) {
		p_soa_var_t *soa = s->start_soa + i;

		free(soa->p);
		free(soa->start_roots);
		free(soa->r.roots);
		CUDA_TRY(cuMemFree(soa->dev_p))
		CUDA_TRY(cuMemFree(soa->dev_start_roots))
	}

	free(s);
}

static void
p_soa_array_reset(p_soa_array_t *s)
{
	uint32 i;

	s->num_arrays = 0;
	s->pp_is_64 = 0;
	for (i = 0; i < s->max_arrays; i++)
		s->start_soa[i].num_p = 0;
}

static void
p_soa_array_start(p_soa_array_t *s, uint32 pp_is_64,
		  CUstream stream)
{
	uint32 i, j, k, m;
	uint32 root_bytes = pp_is_64 ? sizeof(uint64) : sizeof(uint32);

	for (i = j = 0; i < s->max_arrays; i++) {
		p_soa_var_t *soa = s->start_soa + i;
		uint32 num_p = soa->num_p;
		uint32 num_roots = soa->num_roots;
		uint64 *rs = soa->start_roots;
		void *root_array = soa->r.roots;

		if (num_p * num_roots < 50)
			continue;

		if (pp_is_64) {

			if (num_roots == 1) {
				root_array = rs;
			}
			else {
				uint64 *rd = soa->r.roots64;

				for (k = 0; k < num_p; k++) {
					uint64 *rd2 = rd;

					for (m = 0; m < num_roots; m++) {
						*rd2 = rs[m];
						rd2 += num_p;
					}
					rs += num_roots;
					rd++;
				}
			}
		}
		else {
			uint32 *rd = soa->r.roots32;

			for (k = 0; k < num_p; k++) {
				uint32 *rd2 = rd;

				for (m = 0; m < num_roots; m++) {
					*rd2 = (uint32)rs[m];
					rd2 += num_p;
				}
				rs += num_roots;
				rd++;
			}
		}

		CUDA_TRY(cuMemcpyHtoDAsync(soa->dev_p, 
				soa->p,
				num_p * sizeof(uint32),
				stream))
		CUDA_TRY(cuMemcpyHtoDAsync(soa->dev_start_roots, 
				root_array,
				num_p * num_roots * root_bytes,
				stream))
		s->soa[j++] = soa;
	}
	s->num_arrays = j;
	s->pp_is_64 = pp_is_64;
}

static void
store_p_soa(uint32 p, uint32 num_roots, uint64 *roots, void *extra)
{
	uint32 i, j;
	p_soa_array_t *s = (p_soa_array_t *)extra;

	for (i = 0; i < s->max_arrays; i++) {

		p_soa_var_t *soa = s->start_soa + i;
		uint32 num_p;

		if (soa->num_roots != num_roots)
			continue;

		num_p = soa->num_p;

		if (soa->num_p_alloc == num_p) {
			soa->num_p_alloc *= 2;
			soa->p = (uint32 *)xrealloc(soa->p, soa->num_p_alloc *
							sizeof(uint32));
			soa->start_roots = (uint64 *)xrealloc(soa->start_roots,
						soa->num_p_alloc *
						num_roots *
						sizeof(uint64));
			soa->r.roots = xrealloc(soa->r.roots,
						soa->num_p_alloc *
						num_roots *
						sizeof(uint64));

			CUDA_TRY(cuMemFree(soa->dev_p))
			CUDA_TRY(cuMemFree(soa->dev_start_roots))
			CUDA_TRY(cuMemAlloc(&soa->dev_p, 
					soa->num_p_alloc * 
					sizeof(uint32)))
			CUDA_TRY(cuMemAlloc(&soa->dev_start_roots, 
					soa->num_p_alloc *
					soa->num_roots * sizeof(uint64)))
		}

		soa->p[num_p] = p;
		for (j = 0; j < num_roots; j++)
			soa->start_roots[num_p * num_roots + j] = roots[j];
		soa->num_p++;
		return;
	}
}

/*------------------------------------------------------------------------*/

typedef struct {
	uint32 num_specialq;
	uint32 max_specialq;
	specialq_t *specialq;
	CUdeviceptr dev_q;
} specialq_array_t;

static specialq_array_t *
specialq_array_init(void)
{
	specialq_array_t *q_array = (specialq_array_t *)xcalloc(1,
					sizeof(specialq_array_t));

	q_array->max_specialq = 512;
	q_array->specialq = (specialq_t *)xmalloc(q_array->max_specialq *
						sizeof(specialq_t));
	CUDA_TRY(cuMemAlloc(&q_array->dev_q,
				q_array->max_specialq * 
				sizeof(specialq_t)))
	return q_array;
}

static void
specialq_array_free(specialq_array_t *q_array)
{
	CUDA_TRY(cuMemFree(q_array->dev_q))
	free(q_array->specialq);
	free(q_array);
}

static void
specialq_array_reset(specialq_array_t *q_array)
{
	q_array->num_specialq = 0;
}

static void
specialq_array_nextbatch(specialq_array_t *q_array, 
			uint32 num_removed)
{
	uint32 new_size;

	if (num_removed >= q_array->num_specialq) {
		q_array->num_specialq = 0;
		return;
	}

	new_size = q_array->num_specialq - num_removed;
	memmove(q_array->specialq, 
		q_array->specialq + num_removed,
		new_size * sizeof(specialq_t));
	q_array->num_specialq = new_size;
}

static void
specialq_array_start(specialq_array_t *q_array, 
			uint32 num_specialq,
			CUstream stream)
{
	CUDA_TRY(cuMemcpyHtoDAsync(q_array->dev_q,
			q_array->specialq,
			num_specialq * sizeof(specialq_t),
			stream))
}

static void
store_specialq(uint32 q, uint32 num_roots, uint64 *roots, void *extra)
{
	uint32 i;
	uint64 q2 = (uint64)q * q;
	specialq_array_t *q_array = (specialq_array_t *)extra;

	if (q_array->num_specialq + num_roots >= q_array->max_specialq) {

		q_array->max_specialq *= 2;
		q_array->specialq = (specialq_t *)xrealloc(
						q_array->specialq,
						q_array->max_specialq *
						sizeof(specialq_t));
		CUDA_TRY(cuMemFree(q_array->dev_q))
		CUDA_TRY(cuMemAlloc(&q_array->dev_q,
					q_array->max_specialq * 
					sizeof(specialq_t)))
	}

	for (i = 0; i < num_roots; i++) {
		specialq_t *s = q_array->specialq + 
				q_array->num_specialq + i;

		s->p = q;
		s->pp = q2;
		s->root = roots[i];
	}

	q_array->num_specialq += num_roots;
}

typedef struct {

	CUcontext gpu_context;
	CUmodule gpu_module;

	void *sieve_p_fb;
	void *sieve_q_fb;

	p_soa_array_t *p_array;
	specialq_array_t *q_array;

	CUstream stream;

	uint32 num_entries;

	gpu_launch_t *launch;

	CUdeviceptr gpu_p_array;
	CUdeviceptr gpu_p_array_scratch;

	CUdeviceptr gpu_root_array;
	CUdeviceptr gpu_root_array_scratch;

	CUdeviceptr gpu_found_array;
	found_t *found_array;

	void * sort_engine;

	CUevent start_event;
	CUevent end_event;
	double gpu_elapsed;
	double cumulative_elapsed;

} device_thread_data_t;

typedef struct {

	msieve_obj *obj;

	poly_search_t *poly;

	gpu_info_t *gpu_info;

	libhandle_t sort_engine_handle;
	sort_engine_init_func sort_engine_init;
	sort_engine_free_func sort_engine_free;
	sort_engine_run_func sort_engine_run;

	size_t max_sort_entries32;
	size_t max_sort_entries64;

	uint32 num_threads;
	device_thread_data_t *threads;

	struct threadpool *gpu_threadpool;
	struct threadpool *stage2_threadpool;

} device_data_t;

/*------------------------------------------------------------------------*/
/* infrastructure for submitting stage 1 hits to the stage 2 thread pool */

typedef struct {
	stage1_callback_t callback;
	void *callback_data;

	mpz_t ad;
	mpz_t p;
	mpz_t m;

} stage1_hit_data_t;

static void
stage1_hit_free(void *data, int threadid)
{
	stage1_hit_data_t *hit_data = (stage1_hit_data_t *)data;

	mpz_clear(hit_data->ad);
	mpz_clear(hit_data->p);
	mpz_clear(hit_data->m);
	free(hit_data);
}

static void
stage1_hit_run(void *data, int threadid)
{
	stage1_hit_data_t *hit_data = (stage1_hit_data_t *)data;

	hit_data->callback(hit_data->ad,
			   hit_data->p,
			   hit_data->m,
			   hit_data->callback_data);
}

static void
check_found_array(poly_coeff_t *c, device_data_t *d,
			device_thread_data_t *t)
{
	uint32 i;
	uint32 found_array_size;
	found_t *found_array = t->found_array;

	CUDA_TRY(cuMemcpyDtoHAsync(found_array, t->gpu_found_array,
			FOUND_ARRAY_SIZE * sizeof(found_t), t->stream))

	/* we have to synchronize now */

	CUDA_TRY(cuStreamSynchronize(t->stream))

	found_array_size = MIN(FOUND_ARRAY_SIZE - 1,
				found_array[0].p1);

	if (found_array_size == 0)
		return;

	/* clear only the first element */
	CUDA_TRY(cuMemsetD8Async(t->gpu_found_array, 0, 
				sizeof(found_t), t->stream))

	for (i = 1; i <= found_array_size; i++) {
		found_t *found = found_array + i;
		uint32 p1 = found->p1;
		uint32 p2 = found->p2;
		uint32 q = found->q;
		uint64 qroot = found->qroot;
		int64 offset = found->offset;

		double dp = (double)q * p1 * p2;
		double coeff = c->m0 * fabs((double)qroot + 
					(double)offset * q * q) /
					(dp * dp);

		if (coeff <= c->coeff_max &&
		    handle_collision(c, (uint64)p1 * p2, q,
					qroot, offset) != 0) {

			/* submit the hit to the stage 2 thread pool */

			task_control_t task_control;
			stage1_hit_data_t *hit_data = (stage1_hit_data_t *)
						xmalloc(sizeof(stage1_hit_data_t));

			hit_data->callback = d->poly->callback;
			hit_data->callback_data = d->poly->callback_data;
			mpz_init_set(hit_data->ad, c->high_coeff);
			mpz_init_set(hit_data->p, c->p);
			mpz_init_set(hit_data->m, c->m);

			task_control.init = NULL;
			task_control.run = stage1_hit_run;
			task_control.shutdown = stage1_hit_free;
			task_control.data = hit_data;

			threadpool_add_task(d->stage2_threadpool,
						&task_control, 1);
		}
	}
}

#define MAX_SPECIAL_Q ((uint32)(-1))
#define MAX_OTHER ((uint32)1 << 27)

/*------------------------------------------------------------------------*/
static uint32
handle_special_q_batch(msieve_obj *obj, device_data_t *d, 
			device_thread_data_t *t, uint32 num_specialq,
		       	uint32 shift, uint32 key_bits, uint32 num_aprog_vals)
{
	uint32 i, j;
	uint32 quit = 0;
	p_soa_array_t *p_array = t->p_array;
	specialq_array_t *q_array = t->q_array;
	uint32 num_blocks;
	gpu_arg_t gpu_args[GPU_MAX_KERNEL_ARGS];
	sort_data_t sort_data;
	gpu_launch_t *launch;
	uint32 num_q, curr_q;
	uint32 root_bytes = (key_bits > 32) ? sizeof(uint64) : sizeof(uint32);
	float elapsed_ms;

	CUDA_TRY(cuEventRecord(t->start_event, t->stream))

	specialq_array_start(q_array, num_specialq, t->stream);

	CUDA_TRY(cuMemsetD8Async(t->gpu_root_array, 0,
			num_specialq * t->num_entries * 
			num_aprog_vals * root_bytes,
			t->stream))

	for (i = num_q = curr_q = 0; i < num_specialq; i++) {
		if (q_array->specialq[i].p != curr_q) {
			num_q++;
			curr_q = q_array->specialq[i].p;
		}
	}
	
	for (i = j = 0; i < p_array->num_arrays; i++) {
		p_soa_var_t *soa = p_array->soa[i];
		uint32 num_p = soa->num_p;
		uint32 blocks_x, blocks_y;
		uint32 size_x, size_y;
		uint32 total_blocks;

		if (p_array->pp_is_64)
			launch = t->launch + GPU_TRANS_PP64_R64;
		else if (root_bytes == sizeof(uint64))
			launch = t->launch + GPU_TRANS_PP32_R64;
		else
			launch = t->launch + GPU_TRANS_PP32_R32;

		/* perform a block decomposition so that all the
		   soa's generate blocks with about the same amount
		   of arithmetic. There is a modular multiply for
		   each root and a modular inverse for each (p,q) pair, 
		   which we count as 3 multiplies */

		total_blocks = (3 * num_p * num_q +
			        num_p * soa->num_roots * num_specialq) /
				50000;
		total_blocks = MIN(total_blocks, 1000);
		total_blocks = MAX(total_blocks, 1);

		/* choose the number of threads per block to be
		   - a multiple of the warp size between 128 and 256
		   - that can generate the desired number of blocks
		     so the whole dataset is covered, while maximizing
		     the size of the blocks on the borders */

		size_x = MIN(256, launch->threads_per_block);
		while (1) {

			blocks_x = (num_p + size_x - 1) / size_x;
			blocks_y = (total_blocks + blocks_x - 1) / blocks_x;
			size_y = (num_specialq + blocks_y - 1) / blocks_y;

			if (size_x == 128 ||
			    blocks_x * size_x - num_p <= size_x / 3)
				break;

			size_x -= d->gpu_info->warp_size;
		}

		gpu_args[0].ptr_arg = (void *)(size_t)soa->dev_p;
		gpu_args[1].uint32_arg = num_p;
		gpu_args[2].ptr_arg = (void *)(size_t)soa->dev_start_roots;
		gpu_args[3].uint32_arg = soa->num_roots;
		gpu_args[4].ptr_arg = (void *)(size_t)(
				t->gpu_p_array + j * sizeof(uint32));
		gpu_args[5].ptr_arg = (void *)(size_t)(
				t->gpu_root_array + j * root_bytes);
		gpu_args[6].ptr_arg = (void *)(size_t)q_array->dev_q;
		gpu_args[7].uint32_arg = num_specialq;
		gpu_args[8].uint32_arg = size_y;
		gpu_args[9].uint32_arg = t->num_entries;
		gpu_args[10].uint32_arg = shift;
		gpu_args[11].uint32_arg = num_aprog_vals;
		gpu_launch_set(launch, gpu_args);

		CUDA_TRY(cuFuncSetBlockShape(launch->kernel_func, 
				size_x, 1, 1))

		CUDA_TRY(cuLaunchGridAsync(launch->kernel_func,
				blocks_x, blocks_y, t->stream))

		j += num_p * soa->num_roots;
	}

	sort_data.keys_in = t->gpu_root_array;
	sort_data.keys_in_scratch = t->gpu_root_array_scratch;
	sort_data.data_in = t->gpu_p_array;
	sort_data.data_in_scratch = t->gpu_p_array_scratch;
	sort_data.num_elements = num_specialq * t->num_entries * num_aprog_vals;
	sort_data.num_arrays = 1;
	sort_data.key_bits = key_bits;
	sort_data.stream = t->stream;
	d->sort_engine_run(t->sort_engine, &sort_data);

	/* the sort engine may have swapped the input arrays */
	t->gpu_p_array = sort_data.data_in;
	t->gpu_p_array_scratch = sort_data.data_in_scratch;
	t->gpu_root_array = sort_data.keys_in;
	t->gpu_root_array_scratch = sort_data.keys_in_scratch;

	if (root_bytes == sizeof(uint64))
		launch = t->launch + GPU_FINAL_64;
	else
		launch = t->launch + GPU_FINAL_32;

	gpu_args[0].ptr_arg = (void *)(size_t)(t->gpu_p_array);
	gpu_args[1].ptr_arg = (void *)(size_t)(t->gpu_root_array);
	gpu_args[2].uint32_arg = num_specialq * t->num_entries * num_aprog_vals;
	gpu_args[3].ptr_arg = (void *)(size_t)q_array->dev_q;
	gpu_args[4].ptr_arg = (void *)(size_t)(t->gpu_found_array);
	gpu_args[5].uint32_arg = shift;
	gpu_launch_set(launch, gpu_args);

	num_blocks = 1 + (num_specialq * t->num_entries * 
			num_aprog_vals - 1) /
			launch->threads_per_block;
	num_blocks = MIN(num_blocks, 1000);

	CUDA_TRY(cuLaunchGridAsync(launch->kernel_func, 
				num_blocks, 1, t->stream))

	CUDA_TRY(cuEventRecord(t->end_event, t->stream))
	CUDA_TRY(cuEventSynchronize(t->end_event))
	CUDA_TRY(cuEventElapsedTime(&elapsed_ms, 
			t->start_event, t->end_event))
	if (elapsed_ms < 60000 && elapsed_ms > 0) {
		/* this function should execute in under a second. If
		   it takes a very long time, assume that the system
		   was in hibernation and don't let it count. */
		t->gpu_elapsed += elapsed_ms / 1000;
	}

	if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)
		quit = 1;

	return quit;
}

/*------------------------------------------------------------------------*/
static uint32
sieve_specialq(msieve_obj *obj,
		poly_coeff_t *c, device_data_t *d,
		device_thread_data_t *t,
		uint32 special_q_min, uint32 special_q_max,
		uint32 p_min, uint32 p_max, 
		uint32 max_aprog_vals, double deadline)
{
	uint32 i, j;
	uint32 quit = 0;
	uint32 all_q_done = 0;
	uint32 degree = d->poly->degree;
	p_soa_array_t *p_array = t->p_array;
	void *p_fb = t->sieve_p_fb;
	specialq_array_t *q_array = t->q_array;
	void *q_fb = t->sieve_q_fb;
	double cpu_start_time = get_cpu_time();
	uint32 unused_bits;
	uint32 pp_is_64 = (p_max >= 65536);
	uint32 max_batch_specialq32;
	uint32 max_batch_specialq64;
	double elapsed = 0;

	t->gpu_elapsed = 0;

	/* build all the arithmetic progressions */

	p_soa_array_reset(p_array);
	sieve_fb_reset(p_fb, p_min, p_max, 1, p_array->max_p_roots);
	while (sieve_fb_next(p_fb, c, store_p_soa,
			p_array) != P_SEARCH_DONE) {
		;
	}

	p_soa_array_start(p_array, pp_is_64, t->stream);
	if (p_array->num_arrays == 0)
		return 0;

	for (i = j = 0; i < p_array->num_arrays; i++) {
		p_soa_var_t *soa = p_array->soa[i];

		j += soa->num_p * soa->num_roots;
	}

	t->num_entries = j;

	unused_bits = 1;
	while (!(p_max & (1 << (31 - unused_bits))))
		unused_bits++;

	max_batch_specialq32 = d->max_sort_entries32 / t->num_entries;
	max_batch_specialq64 = d->max_sort_entries64 / t->num_entries;

	/* account for 'trivial' special-q */

	specialq_array_reset(q_array);

	if (special_q_min == 1) {
		uint64 trivroots[1] = { 0 };

		store_specialq(1, 1, trivroots, q_array);
	}

	/* handle special-q in batches */

	sieve_fb_reset(q_fb, special_q_min, 
			special_q_max, degree, MAX_ROOTS);

	while (!quit && !all_q_done) {

		uint32 batch_size;
		uint32 max_batch_size;
		uint32 key_bits;
		uint32 num_aprog_vals = 1;

		max_batch_size = pp_is_64 ? 
				max_batch_specialq64 :
				max_batch_specialq32;
		max_batch_size = MIN(max_batch_size, 
				(uint32)1 << unused_bits);
		if (max_batch_size == 0) {
			printf("error: max_batch_size == 0\n");
			exit(-1);
		}

		while (q_array->num_specialq < max_batch_size) {
			if (sieve_fb_next(q_fb, c,
				store_specialq, q_array) == P_SEARCH_DONE) {

				all_q_done = 1;
				break;
			}
		}
		if (q_array->num_specialq == 0)
			continue;

		batch_size = MIN(max_batch_size, q_array->num_specialq);

		if (batch_size < max_batch_size / 3) {

			/* current batch of q is too small to utilize 
			   the card efficiently. Have each (p,q) pair
			   generate multiple offsets, centered about 0 */

			num_aprog_vals = max_batch_size / batch_size;
			num_aprog_vals = MIN(num_aprog_vals, max_aprog_vals);

			/* if we were set up for 32-bit sort keys but
			   now require 64-bit sort keys, make sure to 
			   respect the 64-bit special-q limit */

			key_bits = 1 + ceil(log((double)p_max * p_max *
					((num_aprog_vals + 1) / 2)) / M_LN2);

			if (key_bits > 32)
				num_aprog_vals = MIN(max_batch_size,
				        max_batch_specialq64) / batch_size;
		}

		key_bits = ceil(log((double)p_max * p_max *
				((num_aprog_vals + 1) / 2)) / M_LN2);
		if (num_aprog_vals > 1)
			key_bits++;

		quit = handle_special_q_batch(obj, d, t, batch_size, 
				32 - unused_bits, key_bits, num_aprog_vals);

		check_found_array(c, d, t);

		specialq_array_nextbatch(q_array, batch_size);

		elapsed = get_cpu_time() - cpu_start_time + t->gpu_elapsed;
		if (elapsed > deadline)
			quit = 1;
	}

	t->cumulative_elapsed += elapsed;
	return quit;
}

/*------------------------------------------------------------------------*/
static void
sieve_lattice_gpu_core(msieve_obj *obj,
		poly_coeff_t *c, device_data_t *d, 
		device_thread_data_t *t, double deadline)
{
	uint32 degree = d->poly->degree;
	uint32 num_pieces;
	uint32 p_min, p_max;
	uint32 special_q_min, special_q_max;
	uint32 special_q_min2, special_q_max2;
	uint32 special_q_fb_max;
	double target = c->coeff_max / c->m0;
	uint32 max_aprog_vals = ceil(2 * P_SCALE);

	/* Kleinjung shows that the third-to-largest algebraic
	   polynomial coefficient is of size approximately

	             (correction to m0) * m0
		    --------------------------
		    (leading rational coeff)^2
	
	   We have a bound 'coeff_max' on what this number is 
	   supposed to be, and we know m0 and an upper bound on 
	   the size of the leading rational coefficient P. Let 
	   P = p1*p2*q, where p1 and p2 are drawn from a fixed
	   set of candidates, and q (the 'special-q') is arbitrary
	   except that gcd(q,p1,p2)=1. Then the correction 'C' to 
	   m0 is < max(p1,p2)^2 */
	   
	p_max = MIN(MAX_OTHER, sqrt(c->p_size_max));
	p_max = MIN(p_max, sqrt(0.5 / target));

	special_q_max = MIN(MAX_SPECIAL_Q, 
			    c->p_size_max / p_max / p_max);
	special_q_max = MAX(special_q_max, 1);

	p_min = MAX(1, p_max / P_SCALE);
	special_q_min = 1;

	/* set up the special q factory; special-q may have 
	   arbitrary factors, but many small factors are 
	   preferred since that will allow for many more roots
	   per special q, so we choose the factors to be as 
	   small as possible */

	special_q_fb_max = MIN(200000, special_q_max);
	sieve_fb_init(t->sieve_q_fb, c,
			2, special_q_fb_max,
			1, degree,
			1);

	/* because special-q can have any factors, we require that
	   the progressions we generate use p that have somewhat
	   large factors. This minimizes the chance that a given
	   special-q has factors in common with many progressions
	   in the set */

	sieve_fb_init(t->sieve_p_fb, c, 
			100, 5000,
			1, degree,
		       	0);

	/* large search problems can be randomized so that
	   multiple runs over the same range of leading
	   a_d will likely generate different results */

	num_pieces = 1;
	if (special_q_max - special_q_min > 500000)
		num_pieces = MIN(50, (double)special_q_max * p_max
				/ log(special_q_max) / log(p_max)
				/ 3e10);

	if (num_pieces > 51) { /* randomize the special_q range */
		uint32 piece_length = (special_q_max - special_q_min)
				/ num_pieces;
		uint32 piece = get_rand(&obj->seed1, &obj->seed2)
				% num_pieces;
#if 1
		printf("randomizing rational coefficient: "
			"using piece #%u of %u\n",
			piece + 1, num_pieces);
#endif
		special_q_min2 = special_q_min + piece * piece_length;
		special_q_max2 = special_q_min2 + piece_length;
	}
	else {
		special_q_min2 = special_q_min;
		special_q_max2 = special_q_max;
	}
#if 1
	gmp_printf("coeff %Zd specialq %u - %u other %u - %u\n",
			c->high_coeff,
			special_q_min2, special_q_max2,
			p_min, p_max);
#endif
	sieve_specialq(obj, c, d, t,
			special_q_min2, special_q_max2, p_min, p_max,
			max_aprog_vals, deadline);
}

/*------------------------------------------------------------------------*/
static void
load_sort_engine(msieve_obj *obj, device_data_t *d)
{
	char libname[256];
	#if defined(WIN32) || defined(_WIN64)
	const char *suffix = ".dll";
	#else
	const char *suffix = ".so";
	#endif

	sprintf(libname, "cub/sort_engine%s", suffix);

	/* override from input args */

	if (obj->nfs_args != NULL) {
		char *tmp = strstr(obj->nfs_args, "sortlib=");

		if (tmp != NULL) {
			uint32 i;
			for (i = 0, tmp += 8; i < sizeof(libname) - 1; i++) {
				if (*tmp == 0 || isspace(*tmp))
					break;

				libname[i] = *tmp++;
			}
			libname[i] = 0;
		}
	}

#ifdef _MSC_VER && _MSC_VER >= 1900
/*  convert the sort engine file path to an absolute path using
    backslash directory separators  */
	char libpath[256], *p, *q;
	int len;

	_getcwd(libpath, 256);
	len = strlen(libpath);
	p = libpath + len;
	if(*(p - 1) != '\\')
	{
		*p++ = '\\';
		len++;
	}
	q = libname - 1;
	while(*++q && len < 256)
	{
		*p++ = (*q != '/' ? *q : '\\');
		len++;
	}
	if(len < 256) {
		*p++ = 0;
	}
	else {
		printf("error: buffer overflow in %s at line %d\n", __FILE__, __LINE__);
		exit(-1);
	}

	if(_access_s(libpath, 0)) {
		printf("error: GPU sort engine not found at \"%s\"\n", libpath);
		exit(-1);
	}

	d->sort_engine_handle = load_dynamic_lib(libpath);
	if(d->sort_engine_handle == NULL) {
		printf("error: cannot load GPU sort engine from \"%s\" (check dependencies)\n", libpath);
		exit(-1);
}
#else
	d->sort_engine_handle = load_dynamic_lib(libname);
	if(d->sort_engine_handle == NULL) {
		printf("error: failed to load GPU sorting engine from \"%s\"\n", libname);
		exit(-1);
	}
#endif
	/* the sort engine uses the same CUDA context */

	d->sort_engine_init = get_lib_symbol(
					d->sort_engine_handle,
					"sort_engine_init");
	d->sort_engine_free = get_lib_symbol(
					d->sort_engine_handle,
					"sort_engine_free");
	d->sort_engine_run = get_lib_symbol(
					d->sort_engine_handle,
					"sort_engine_run");
	if (d->sort_engine_init == NULL ||
	    d->sort_engine_free == NULL ||
	    d->sort_engine_run == NULL) {
		printf("error: cannot find GPU sorting function\n");
		exit(-1);
	}
}


/*------------------------------------------------------------------------*/
static void
gpu_thread_data_init(void *data, int threadid)
{
	uint32 i, j;
	device_data_t *d = (device_data_t *)data;
	device_thread_data_t *t = d->threads + threadid;

	/* every thread needs its own context; making all
	   threads share the same context causes problems
	   with the sort engine, because apparently it
	   changes the GPU cache size on the fly */

	CUDA_TRY(cuCtxCreate(&t->gpu_context, 
			CU_CTX_BLOCKING_SYNC,
			d->gpu_info->device_handle))

	/* load GPU kernels */

	if (d->gpu_info->compute_version_major == 2) {
		CUDA_TRY(cuModuleLoad(&t->gpu_module, "stage1_core_sm20.ptx"))
	}
	else if (d->gpu_info->compute_version_major == 3) {
		if (d->gpu_info->compute_version_minor < 5)
			CUDA_TRY(cuModuleLoad(&t->gpu_module, "stage1_core_sm30.ptx"))
		else
			CUDA_TRY(cuModuleLoad(&t->gpu_module, "stage1_core_sm35.ptx"))
	}
	else if (d->gpu_info->compute_version_major == 5) {
		if (d->gpu_info->compute_version_minor < 2)
			CUDA_TRY(cuModuleLoad(&t->gpu_module, "stage1_core_sm50.ptx"))
		else
			CUDA_TRY(cuModuleLoad(&t->gpu_module, "stage1_core_sm52.ptx"))
	}
	else if (d->gpu_info->compute_version_major == 6) {
		if (d->gpu_info->compute_version_minor < 1)
			CUDA_TRY(cuModuleLoad(&t->gpu_module, "stage1_core_sm60.ptx"))
		else
			CUDA_TRY(cuModuleLoad(&t->gpu_module, "stage1_core_sm61.ptx"))
	}
	else if (d->gpu_info->compute_version_major == 7) {
		if (d->gpu_info->compute_version_minor < 5)
			CUDA_TRY(cuModuleLoad(&t->gpu_module, "stage1_core_sm70.ptx"))
		else
			CUDA_TRY(cuModuleLoad(&t->gpu_module, "stage1_core_sm75.ptx"))
	}
	else if (d->gpu_info->compute_version_major >= 8) {
		if (d->gpu_info->compute_version_minor < 6)
			CUDA_TRY(cuModuleLoad(&t->gpu_module, "stage1_core_sm80.ptx"))
		else if (d->gpu_info->compute_version_minor < 9)
			CUDA_TRY(cuModuleLoad(&t->gpu_module, "stage1_core_sm86.ptx"))
		else
     		CUDA_TRY(cuModuleLoad(&t->gpu_module, "stage1_core_sm89.ptx"))
	}
	else
	{
	    printf("sorry, Nvidia doesn't want to support your card\n");
		exit(-1);
	}

	t->launch = (gpu_launch_t *)xmalloc(NUM_GPU_FUNCTIONS *
				sizeof(gpu_launch_t));

	for (i = 0; i < NUM_GPU_FUNCTIONS; i++) {
		gpu_launch_t *launch = t->launch + i;

		gpu_launch_init(t->gpu_module, gpu_kernel_names[i],
				gpu_kernel_args + (i / 3), launch);

		if (i == GPU_FINAL_32 || i == GPU_FINAL_64) {
			/* performance of the cleanup functions is not 
			   that sensitive to the block shape; set it 
			   once up front */

			launch->threads_per_block = 
					MIN(256, launch->threads_per_block);
			CUDA_TRY(cuFuncSetBlockShape(launch->kernel_func,
					launch->threads_per_block, 1, 1))
		}
	}

	/* threads each send a stream of kernel calls */

	CUDA_TRY(cuStreamCreate(&t->stream, 0))

	/* set up found array */

	CUDA_TRY(cuMemAlloc(&t->gpu_found_array, sizeof(found_t) *
			FOUND_ARRAY_SIZE))
	t->found_array = (found_t *)xmalloc(sizeof(found_t) * 
			FOUND_ARRAY_SIZE);

	CUDA_TRY(cuMemsetD8(t->gpu_found_array, 0, sizeof(found_t)))

	/* set up root generation arrays */

	t->sieve_p_fb = sieve_fb_alloc();
	t->sieve_q_fb = sieve_fb_alloc();

	t->p_array = p_soa_array_init(d->poly->degree);
	t->q_array = specialq_array_init();

	i = sizeof(uint32) * MAX(d->max_sort_entries32, d->max_sort_entries64);
	j = MAX(d->max_sort_entries32 * sizeof(uint32),
	        d->max_sort_entries64 * sizeof(uint64));

	CUDA_TRY(cuMemAlloc(&t->gpu_p_array, i))
	CUDA_TRY(cuMemAlloc(&t->gpu_p_array_scratch, i))
	CUDA_TRY(cuMemAlloc(&t->gpu_root_array, j))
	CUDA_TRY(cuMemAlloc(&t->gpu_root_array_scratch, j))

	t->sort_engine = d->sort_engine_init();

	CUDA_TRY(cuEventCreate(&t->start_event, CU_EVENT_BLOCKING_SYNC))
	CUDA_TRY(cuEventCreate(&t->end_event, CU_EVENT_BLOCKING_SYNC))
}


/*------------------------------------------------------------------------*/
static void
gpu_thread_data_free(void *data, int threadid)
{
	device_data_t *d = (device_data_t *)data;
	device_thread_data_t *t = d->threads + threadid;

	CUDA_TRY(cuEventDestroy(t->start_event))
	CUDA_TRY(cuEventDestroy(t->end_event))

	d->sort_engine_free(t->sort_engine);

	CUDA_TRY(cuMemFree(t->gpu_p_array))
	CUDA_TRY(cuMemFree(t->gpu_p_array_scratch))
	CUDA_TRY(cuMemFree(t->gpu_root_array))
	CUDA_TRY(cuMemFree(t->gpu_root_array_scratch))

	sieve_fb_free(t->sieve_p_fb);
	sieve_fb_free(t->sieve_q_fb);

	CUDA_TRY(cuStreamDestroy(t->stream))

	free(t->found_array);
	CUDA_TRY(cuMemFree(t->gpu_found_array))

	free(t->launch);

	p_soa_array_free(t->p_array);
	specialq_array_free(t->q_array);

	CUDA_TRY(cuCtxDestroy(t->gpu_context)) 
}

/*------------------------------------------------------------------------*/
void *
gpu_data_init(msieve_obj *obj, poly_search_t *poly)
{
	device_data_t *d;
	gpu_config_t gpu_config;
	gpu_info_t *gpu_info;
	size_t gpu_mem;

	uint32 num_threads;
	thread_control_t thread_control;

	gpu_init(&gpu_config);
	if (gpu_config.num_gpu == 0) {
		printf("error: no CUDA-enabled GPUs found\n");
		exit(-1);
	}
	if (obj->which_gpu >= (uint32)gpu_config.num_gpu) {
		printf("error: GPU %u does not exist "
			"or is not CUDA-enabled\n", obj->which_gpu);
		exit(-1);
	}

	d = (device_data_t *)xcalloc(1, sizeof(device_data_t));

	d->obj = obj;
	d->gpu_info = gpu_info = (gpu_info_t *)xmalloc(sizeof(gpu_info_t));
	memcpy(gpu_info, gpu_config.info + obj->which_gpu,
			sizeof(gpu_info_t)); 

	logprintf(obj, "using GPU %u (%s)\n", obj->which_gpu, gpu_info->name);
	logprintf(obj, "selected card has CUDA arch %d.%d\n",
			gpu_info->compute_version_major,
			gpu_info->compute_version_minor);

	load_sort_engine(obj, d);

	/* a single transformed array will not have enough
	   elements for the GPU to sort efficiently; instead,
	   we batch multiple arrays together and sort the
	   entire batch. Array elements get their array number
	   added into unused bits above each p, so the number
	   of such unused bits is one limit on the batch size.
	   Another limit is the amount of RAM on the card;
	   we try to use less than 30% of it. Finally, we top
	   out at a few ten millions of elements, which is far 
	   into asymptotically optimal throughput territory 
	   for all current GPUs
	
	   The sizing below assumes we can batch as many special-Q
	   as we want. However, for small degree-5 problems the 
	   number of special-Q available is too small to efficiently 
	   use the card. In that case, we will compensate later by 
	   making each arithmetic progression contribute multiple 
	   offsets */

	gpu_mem = 0.3 * d->gpu_info->global_mem_size;
	if (obj->nfs_args != NULL) {

		char *tmp = strstr(obj->nfs_args, "gpu_mem_mb=");
		if (tmp != NULL) {
			size_t m = strtoul(tmp + 11, NULL, 10);

			gpu_mem = MIN(d->gpu_info->global_mem_size,
					m * 1048576);
			gpu_mem = MAX(10 * 1048576, gpu_mem);
			logprintf(obj, "setting max GPU mem use to %.1lf MB\n",
					(double)gpu_mem / 1048576);
		}
	}

	d->max_sort_entries32 = MIN(50000000, gpu_mem /
				 (2 * (sizeof(uint32) + 
				       sizeof(uint32))));

	d->max_sort_entries64 = MIN(35000000, gpu_mem /
				 (2 * (sizeof(uint32) + 
				       sizeof(uint64))));

	/* account for multiple threads; we allocate a thread pool
	   with a number of threads requested, where each thread
	   deals with a single leading coefficient. We also allocate
	   another thread pool with a single thread, that runs stage
	   2. Eventually the latter can be made more concurrent.  
	 
	   Each GPU thread gets an equal share of the total memory
	   used on the card. The GPU thread pool should only accumulate
	   a few leading coefficients at a time, but the stage 2 thread
	   pool should have a deeper queue of work */

	d->poly = poly;
	num_threads = MAX(1, obj->num_threads);
	d->num_threads = num_threads = MIN(4, num_threads);
	d->max_sort_entries32 /= num_threads;
	d->max_sort_entries64 /= num_threads;
	d->threads = (device_thread_data_t *)xcalloc(num_threads,
					sizeof(device_thread_data_t));

	thread_control.init = gpu_thread_data_init;
	thread_control.shutdown = gpu_thread_data_free;
	thread_control.data = d;
	d->gpu_threadpool = threadpool_init(num_threads,
					MAX(10, num_threads),
					&thread_control);

	thread_control.init = NULL;
	thread_control.shutdown = NULL;
	thread_control.data = NULL;
	d->stage2_threadpool = threadpool_init(1, 1000, &thread_control);

	return d;
}

/*------------------------------------------------------------------------*/
void gpu_data_free(void *gpu_data)
{
	device_data_t *d = (device_data_t *)gpu_data;

	if (!(d->obj->flags & MSIEVE_FLAG_STOP_SIEVING)) {
		/* we're allowed to try to shut down gracefully */

		threadpool_drain(d->gpu_threadpool, 1);
	}

	/* shut down the GPU threadpool first, since
	   we don't want it feeding the stage 2 threadpool
	   after it's been freed */

	threadpool_free(d->gpu_threadpool);
	threadpool_free(d->stage2_threadpool);

	free(d->threads);

	unload_dynamic_lib(d->sort_engine_handle);

	free(d->gpu_info);
	free(d);
}

/*------------------------------------------------------------------------*/
/* infrastructure for submitting new leading
   coeffs to the GPU thread pool */

typedef struct {
	msieve_obj *obj;
	poly_coeff_t *c;
	device_data_t *d;
	uint32 deadline;
} task_data_t;

static void
task_data_free(void *data, int threadid)
{
	task_data_t *task_data = (task_data_t *)data;

	poly_coeff_free(task_data->c);
	free(task_data);
}

static void
task_data_run(void *data, int threadid)
{
	task_data_t *task_data = (task_data_t *)data;

	sieve_lattice_gpu_core(task_data->obj, task_data->c, task_data->d,
				task_data->d->threads + threadid,
				task_data->deadline);
}

/* external entry point */
double sieve_lattice_gpu(msieve_obj *obj, poly_search_t *poly,
			poly_coeff_t *c, void *gpu_data,
			double deadline)
{
	/* submit a leading coefficient asynchronously to the
	   GPU thread pool; we copy the coefficient so the input
	   one can be overwritten by calling code */

	uint32 i;
	task_control_t task_control;
	task_data_t *task_data = (task_data_t *)xmalloc(sizeof(task_data_t));
	device_data_t *d = (device_data_t *)gpu_data;
	poly_coeff_t *c2 = poly_coeff_init();
	double cumulative_elapsed = 0;

	poly_coeff_copy(c2, c);

	task_data->obj = obj;
	task_data->c = c2;
	task_data->d = d;
	task_data->deadline = deadline;

	task_control.init = NULL;
	task_control.run = task_data_run;
	task_control.shutdown = task_data_free;
	task_control.data = task_data;

	threadpool_add_task(d->gpu_threadpool, &task_control, 1);

	/* return total time spent by all stage 1 threads */

	for (i = 0; i < d->num_threads; i++)
		cumulative_elapsed += d->threads[i].cumulative_elapsed;

	return cumulative_elapsed;
}


