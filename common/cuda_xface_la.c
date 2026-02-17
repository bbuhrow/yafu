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

#include <cuda_xface_la.h>

#ifdef  HAVE_LA_CUDA

/*------------------------------------------------------------------------*/
void
cuGetErrorMessageLA(CUresult result, int line) 
{
	char * error_name = NULL;
	char * error_string = NULL;

	cuGetErrorName(result, &error_name);
	cuGetErrorString(result, &error_string);
	printf("error (line %d): %s, %s\n", line, error_name, error_string);
}

/*------------------------------------------------------------------------*/
void
gpu_init_la(gpu_config_la_t *config)
{
	int32 i;

	/* determine the specifics of CUDA GPUs using the 14 different
	   methods in the Nvidia documentation */

	memset(config, 0, sizeof(gpu_config_la_t));

	CUDA_TRY(cuInit(0))
	CUDA_TRY(cuDeviceGetCount(&config->num_gpu))
	if (config->num_gpu == 0)
		return;

	for (i = 0; i < (int32)config->num_gpu; i++) {
		CUdevice device;
		gpu_info_la_t *info = config->info + i;

		CUDA_TRY(cuDeviceGet(&device, i))

		info->device_handle = device;

		CUDA_TRY(cuDeviceGetName(info->name,
				sizeof(info->name), device))
		CUDA_TRY(cuDeviceGetAttribute(
				&info->compute_version_major,
				CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, device))
		CUDA_TRY(cuDeviceGetAttribute(
				&info->compute_version_minor,
				CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, device))
		CUDA_TRY(cuDeviceGetAttribute(
				&info->clock_speed,
				CU_DEVICE_ATTRIBUTE_CLOCK_RATE, device))
		CUDA_TRY(cuDeviceGetAttribute(
				&info->constant_mem_size,
				CU_DEVICE_ATTRIBUTE_TOTAL_CONSTANT_MEMORY, device))
		CUDA_TRY(cuDeviceGetAttribute(
				&info->shared_mem_size,
				CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK, device))
		CUDA_TRY(cuDeviceGetAttribute(
				&info->registers_per_block,
				CU_DEVICE_ATTRIBUTE_MAX_REGISTERS_PER_BLOCK, device))
		CUDA_TRY(cuDeviceGetAttribute(
				&info->max_threads_per_block,
				CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK, device))
		CUDA_TRY(cuDeviceGetAttribute(
				&info->warp_size,
				CU_DEVICE_ATTRIBUTE_WARP_SIZE, device))
		CUDA_TRY(cuDeviceGetAttribute(
				&info->max_thread_dim[0],
				CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_X, device))
		CUDA_TRY(cuDeviceGetAttribute(
				&info->max_thread_dim[1],
				CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_Y, device))
		CUDA_TRY(cuDeviceGetAttribute(
				&info->max_thread_dim[2],
				CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_Z, device))
		CUDA_TRY(cuDeviceGetAttribute(
				&info->max_grid_size[0],
				CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_X, device))
		CUDA_TRY(cuDeviceGetAttribute(
				&info->max_grid_size[1],
				CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_Y, device))
		CUDA_TRY(cuDeviceGetAttribute(
				&info->max_grid_size[2],
				CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_Z, device))
		CUDA_TRY(cuDeviceTotalMem(
			&info->global_mem_size, device))
		CUDA_TRY(cuDeviceGetAttribute(&info->can_overlap,
				CU_DEVICE_ATTRIBUTE_GPU_OVERLAP, device))
		CUDA_TRY(cuDeviceGetAttribute(&info->num_compute_units,
				CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT,
				device))
		CUDA_TRY(cuDeviceGetAttribute(&info->has_timeout,
				CU_DEVICE_ATTRIBUTE_KERNEL_EXEC_TIMEOUT,
				device))
		CUDA_TRY(cuDeviceGetAttribute(&info->concurrent_managed_access,
				CU_DEVICE_ATTRIBUTE_CONCURRENT_MANAGED_ACCESS,
				device))
	}
}

/*------------------------------------------------------------------------*/
void gpu_launch_init_la(CUmodule gpu_module, const char *func_name,
			gpu_launch_la_t *launch)
{
	memset(launch, 0, sizeof(gpu_launch_la_t));

	CUDA_TRY(cuModuleGetFunction(&launch->kernel_func,
			gpu_module, func_name))

	CUDA_TRY(cuFuncGetAttribute(&launch->threads_per_block,
			CU_FUNC_ATTRIBUTE_MAX_THREADS_PER_BLOCK,
			launch->kernel_func))
}
#endif
