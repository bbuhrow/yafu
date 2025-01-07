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

#ifndef _COMMON_LANCZOS_GPU_LANCZOS_GPU_H_
#define _COMMON_LANCZOS_GPU_LANCZOS_GPU_H_

#include <cuda_xface_la.h>
#include <spmv_engine.h>
#include "../lanczos.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	uint32 num_rows;
	uint32 num_cols;
	uint32 num_col_entries;         /* but int32 in cub */
	uint32 blocksize;
	CUdeviceptr col_entries;        /* uint32 */
	CUdeviceptr row_entries;        /* uint32 */
} block_row_t;

/* implementation-specific structure */

typedef struct {

	gpu_info_la_t *gpu_info;

	CUcontext gpu_context;
	CUmodule gpu_module;

	gpu_launch_la_t *launch;

	/* gpu product data */

	CUdeviceptr gpu_scratch;

	/* matrix data */

	CUdeviceptr *dense_blocks;

	uint32 num_block_rows;
	block_row_t *block_rows;

	uint32 num_trans_block_rows;
	block_row_t *trans_block_rows;

	/* scan engine data */

	libhandle_t spmv_engine_handle;
	spmv_engine_init_func spmv_engine_init;
	spmv_engine_free_func spmv_engine_free;
	spmv_engine_run_func spmv_engine_run;
	void * spmv_engine;

	/* use managed memory to store the matrix data */
	uint32 use_cudamanaged;

} gpudata_t;


typedef struct {
	gpudata_t *gpudata;
	v_t *host_vec;
	CUdeviceptr gpu_vec;
} gpuvec_t;

/* #define LANCZOS_GPU_DEBUG */

/* ordinal list of GPU kernels */
enum {
	GPU_K_MASK = 0,
	GPU_K_XOR,
	GPU_K_INNER_PROD,
	GPU_K_OUTER_PROD,
	NUM_GPU_FUNCTIONS /* must be last */
};

void vv_xor_gpu(void *dest, void *src, uint32 n, gpudata_t *d);

void mul_BxN_NxB_gpu(packed_matrix_t *matrix,
		   CUdeviceptr x, CUdeviceptr y,
		   CUdeviceptr xy, uint32 n);

void mul_NxB_BxB_acc_gpu(packed_matrix_t *matrix, 
			CUdeviceptr v, CUdeviceptr x,
			CUdeviceptr y, uint32 n);

#ifdef __cplusplus
}
#endif

#endif /* !_COMMON_LANCZOS_GPU_LANCZOS_GPU_H_ */
