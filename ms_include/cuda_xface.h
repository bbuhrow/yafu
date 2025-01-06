/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: cuda_xface.h 817 2012-11-11 14:58:29Z jasonp_sf $
--------------------------------------------------------------------*/

#ifndef _GPU_XFACE_H
#define _GPU_XFACE_H

#if defined(HAVE_CUDA)

#include <util.h>
#include <cuda.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_GPU 4

typedef struct {
	char name[32];
	int32 compute_version_major;
	int32 compute_version_minor;
	int32 clock_speed; /* in kHz */
	int32 num_compute_units;
	int32 constant_mem_size;
	int32 shared_mem_size;
	size_t global_mem_size;
	int32 registers_per_block;
	int32 max_threads_per_block;
	int32 can_overlap;
	int32 warp_size;
	int32 max_thread_dim[3];
	int32 max_grid_size[3];
	int32 has_timeout;
	CUdevice device_handle;
} gpu_info_t;

typedef struct {
	int32 num_gpu;
	gpu_info_t info[MAX_GPU];
} gpu_config_t;

char * cuGetErrorMessage(CUresult result);

void gpu_init(gpu_config_t *config);

#define CUDA_TRY(func) \
	{ 			 				\
		CUresult status = func;				\
		if (status != CUDA_SUCCESS) {			\
			printf("error (line %d): %s\n", __LINE__,\
				cuGetErrorMessage(status));	\
			exit(-1);				\
		}						\
	}

#define CUDA_ALIGN_PARAM(offset, alignment) \
	(offset) = ((offset) + (alignment) - 1) & ~((alignment) - 1)

/* defines for streamlining the handling of arguments to GPU kernels */

typedef enum {
	GPU_ARG_NONE = 0,
	GPU_ARG_PTR,
	GPU_ARG_INT32,
	GPU_ARG_UINT32,
	GPU_ARG_INT64,
	GPU_ARG_UINT64
} gpu_arg_type_t;

#define GPU_MAX_KERNEL_ARGS 15

typedef struct {
	uint32 num_args;
	gpu_arg_type_t arg_type[GPU_MAX_KERNEL_ARGS];
} gpu_arg_type_list_t;

typedef union {
	void * ptr_arg;
	int32 int32_arg;
	uint32 uint32_arg;
	int64 int64_arg;
	uint64 uint64_arg;
} gpu_arg_t;

typedef struct {
	CUfunction kernel_func;
	int32 threads_per_block;
	int32 arg_offsets[GPU_MAX_KERNEL_ARGS];
	gpu_arg_type_list_t arg_desc;
} gpu_launch_t;

void gpu_launch_init(CUmodule gpu_module, const char *func_name,
			const gpu_arg_type_list_t *arg_desc,
			gpu_launch_t *launch);

void gpu_launch_set(gpu_launch_t *launch, gpu_arg_t *args);

#ifdef __cplusplus
}
#endif

#endif /* HAVE_CUDA */

#endif /* !_CUDA_XFACE_H_ */
