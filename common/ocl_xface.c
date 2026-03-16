/*--------------------------------------------------------------------
 * ocl_xface.c  --  OpenCL replacement for cuda_xface.c
 *
 * Implements:
 *   clGetErrorString  -- replaces cuGetErrorMessage
 *   gpu_init          -- replaces CUDA gpu_init (device enumeration)
 *   gpu_launch_init   -- replaces CUDA version (no offset arithmetic)
 *   gpu_launch_set    -- replaces CUDA version (uses clSetKernelArg)
 *--------------------------------------------------------------------*/

#include "ocl_xface.h"
#include <stdint.h>

#if defined(HAVE_OCL_BATCH_FACTOR)

/* -----------------------------------------------------------------------
 * Error string table
 * Replaces cuGetErrorMessage().  OpenCL has clGetErrorString in some
 * implementations, but it is not part of the standard, so we provide our
 * own complete table.
 * --------------------------------------------------------------------- */
const char *
clGetErrorString(cl_int err)
{
    switch (err) {
    case CL_SUCCESS:                          return "CL_SUCCESS";
    case CL_DEVICE_NOT_FOUND:                 return "CL_DEVICE_NOT_FOUND";
    case CL_DEVICE_NOT_AVAILABLE:             return "CL_DEVICE_NOT_AVAILABLE";
    case CL_COMPILER_NOT_AVAILABLE:           return "CL_COMPILER_NOT_AVAILABLE";
    case CL_MEM_OBJECT_ALLOCATION_FAILURE:    return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
    case CL_OUT_OF_RESOURCES:                 return "CL_OUT_OF_RESOURCES";
    case CL_OUT_OF_HOST_MEMORY:               return "CL_OUT_OF_HOST_MEMORY";
    case CL_PROFILING_INFO_NOT_AVAILABLE:     return "CL_PROFILING_INFO_NOT_AVAILABLE";
    case CL_MEM_COPY_OVERLAP:                 return "CL_MEM_COPY_OVERLAP";
    case CL_IMAGE_FORMAT_MISMATCH:            return "CL_IMAGE_FORMAT_MISMATCH";
    case CL_IMAGE_FORMAT_NOT_SUPPORTED:       return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
    case CL_BUILD_PROGRAM_FAILURE:            return "CL_BUILD_PROGRAM_FAILURE";
    case CL_MAP_FAILURE:                      return "CL_MAP_FAILURE";
    case CL_INVALID_VALUE:                    return "CL_INVALID_VALUE";
    case CL_INVALID_DEVICE_TYPE:              return "CL_INVALID_DEVICE_TYPE";
    case CL_INVALID_PLATFORM:                 return "CL_INVALID_PLATFORM";
    case CL_INVALID_DEVICE:                   return "CL_INVALID_DEVICE";
    case CL_INVALID_CONTEXT:                  return "CL_INVALID_CONTEXT";
    case CL_INVALID_QUEUE_PROPERTIES:         return "CL_INVALID_QUEUE_PROPERTIES";
    case CL_INVALID_COMMAND_QUEUE:            return "CL_INVALID_COMMAND_QUEUE";
    case CL_INVALID_HOST_PTR:                 return "CL_INVALID_HOST_PTR";
    case CL_INVALID_MEM_OBJECT:               return "CL_INVALID_MEM_OBJECT";
    case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:  return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
    case CL_INVALID_IMAGE_SIZE:               return "CL_INVALID_IMAGE_SIZE";
    case CL_INVALID_SAMPLER:                  return "CL_INVALID_SAMPLER";
    case CL_INVALID_BINARY:                   return "CL_INVALID_BINARY";
    case CL_INVALID_BUILD_OPTIONS:            return "CL_INVALID_BUILD_OPTIONS";
    case CL_INVALID_PROGRAM:                  return "CL_INVALID_PROGRAM";
    case CL_INVALID_PROGRAM_EXECUTABLE:       return "CL_INVALID_PROGRAM_EXECUTABLE";
    case CL_INVALID_KERNEL_NAME:              return "CL_INVALID_KERNEL_NAME";
    case CL_INVALID_KERNEL_DEFINITION:        return "CL_INVALID_KERNEL_DEFINITION";
    case CL_INVALID_KERNEL:                   return "CL_INVALID_KERNEL";
    case CL_INVALID_ARG_INDEX:                return "CL_INVALID_ARG_INDEX";
    case CL_INVALID_ARG_VALUE:                return "CL_INVALID_ARG_VALUE";
    case CL_INVALID_ARG_SIZE:                 return "CL_INVALID_ARG_SIZE";
    case CL_INVALID_KERNEL_ARGS:              return "CL_INVALID_KERNEL_ARGS";
    case CL_INVALID_WORK_DIMENSION:           return "CL_INVALID_WORK_DIMENSION";
    case CL_INVALID_WORK_GROUP_SIZE:          return "CL_INVALID_WORK_GROUP_SIZE";
    case CL_INVALID_WORK_ITEM_SIZE:           return "CL_INVALID_WORK_ITEM_SIZE";
    case CL_INVALID_GLOBAL_OFFSET:            return "CL_INVALID_GLOBAL_OFFSET";
    case CL_INVALID_EVENT_WAIT_LIST:          return "CL_INVALID_EVENT_WAIT_LIST";
    case CL_INVALID_EVENT:                    return "CL_INVALID_EVENT";
    case CL_INVALID_OPERATION:                return "CL_INVALID_OPERATION";
    case CL_INVALID_GL_OBJECT:                return "CL_INVALID_GL_OBJECT";
    case CL_INVALID_BUFFER_SIZE:              return "CL_INVALID_BUFFER_SIZE";
    case CL_INVALID_GLOBAL_WORK_SIZE:         return "CL_INVALID_GLOBAL_WORK_SIZE";
    default:                                  return "unknown OpenCL error";
    }
}

/* -----------------------------------------------------------------------
 * gpu_init
 *
 * Replaces the CUDA version which called cuInit / cuDeviceGetCount /
 * cuDeviceComputeCapability / cuDeviceGetProperties etc.
 *
 * In OpenCL:
 *   - We enumerate cl_platform_id first, then cl_device_id inside each
 *     platform.  We collect up to MAX_GPU GPU devices total.
 *   - There is no single "compute capability" version; we report the
 *     OpenCL version string parsed into major/minor instead.
 *   - CL_DEVICE_LOCAL_MEM_SIZE replaces sharedMemPerBlock.
 *   - CL_DEVICE_PREFERRED_WORK_GROUP_SIZE_MULTIPLE replaces warp_size.
 *   - Registers per block have no OpenCL equivalent; we set 0.
 * --------------------------------------------------------------------- */
void
gpu_init(gpu_config_t *config)
{
    cl_uint         num_platforms = 0;
    cl_platform_id  platforms[8];
    cl_uint         num_devices;
    cl_device_id    devices[MAX_GPU];
    cl_int          err;
    int             i, j;

    memset(config, 0, sizeof(gpu_config_t));

    err = clGetPlatformIDs(8, platforms, &num_platforms);
    if (err != CL_SUCCESS || num_platforms == 0) {
        printf("gpu_init: no OpenCL platforms found\n");
        return;
    }

    for (i = 0; i < (int)num_platforms && config->num_gpu < MAX_GPU; i++) {

        num_devices = 0;
        err = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_GPU,
                             MAX_GPU - config->num_gpu,
                             devices, &num_devices);
        if (err != CL_SUCCESS || num_devices == 0)
            continue;

        for (j = 0; j < (int)num_devices && config->num_gpu < MAX_GPU; j++) {

            gpu_info_t *info = config->info + config->num_gpu;
            cl_device_id dev  = devices[j];
            char ver_str[64]  = {0};
            size_t wg_mult    = 0;
            size_t max_wg     = 0;
            size_t max_wi[3]  = {0, 0, 0};
            cl_ulong lmem     = 0, cmem = 0, gmem = 0;
            cl_uint  cu = 0, freq = 0;

            info->device_handle   = dev;
            info->platform_handle = platforms[i];

            clGetDeviceInfo(dev, CL_DEVICE_NAME,
                            sizeof(info->name), info->name, NULL);

            /* Parse "OpenCL X.Y" from CL_DEVICE_OPENCL_C_VERSION */
            clGetDeviceInfo(dev, CL_DEVICE_OPENCL_C_VERSION,
                            sizeof(ver_str), ver_str, NULL);
            /* ver_str is e.g. "OpenCL C 2.0" */
            {
                int maj = 1, min = 0;
                sscanf(ver_str, "OpenCL C %d.%d", &maj, &min);
                info->compute_version_major = maj;
                info->compute_version_minor = min;
            }

            clGetDeviceInfo(dev, CL_DEVICE_MAX_CLOCK_FREQUENCY,
                            sizeof(freq), &freq, NULL);
            info->clock_speed = (int32_t)freq;   /* MHz */

            clGetDeviceInfo(dev, CL_DEVICE_MAX_COMPUTE_UNITS,
                            sizeof(cu), &cu, NULL);
            info->num_compute_units = (int32_t)cu;

            clGetDeviceInfo(dev, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE,
                            sizeof(cmem), &cmem, NULL);
            info->constant_mem_size = (int32_t)(cmem > INT32_MAX ? INT32_MAX : cmem);

            clGetDeviceInfo(dev, CL_DEVICE_LOCAL_MEM_SIZE,
                            sizeof(lmem), &lmem, NULL);
            info->shared_mem_size = (int32_t)(lmem > INT32_MAX ? INT32_MAX : lmem);

            clGetDeviceInfo(dev, CL_DEVICE_GLOBAL_MEM_SIZE,
                            sizeof(gmem), &gmem, NULL);
            info->global_mem_size = (size_t)gmem;

            info->registers_per_block = 0;   /* not exposed in OpenCL */

            clGetDeviceInfo(dev, CL_DEVICE_MAX_WORK_GROUP_SIZE,
                            sizeof(max_wg), &max_wg, NULL);
            info->max_threads_per_block = (int32_t)max_wg;

            info->can_overlap = 1;   /* all OpenCL queues can overlap */

            clGetDeviceInfo(dev, CL_DEVICE_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
                            sizeof(wg_mult), &wg_mult, NULL);
            info->warp_size = (int32_t)wg_mult;

            clGetDeviceInfo(dev, CL_DEVICE_MAX_WORK_ITEM_SIZES,
                            sizeof(max_wi), max_wi, NULL);
            for (int k = 0; k < 3; k++) {
                info->max_thread_dim[k]  = (int32_t)(max_wi[k] > INT32_MAX
                                                     ? INT32_MAX : max_wi[k]);
                /* OpenCL max NDRange per dimension = max global work size.
                 * Use max_wg * max_cu as a conservative estimate since
                 * CL_DEVICE_MAX_WORK_ITEM_SIZES is per-dimension WG only. */
                info->max_grid_size[k]   = INT32_MAX;
            }

            info->has_timeout = 0;   /* compute devices don't time out */

            config->num_gpu++;
        }
    }
}

/* -----------------------------------------------------------------------
 * gpu_launch_init
 *
 * CUDA version:
 *   - Called cuModuleGetFunction to get a CUfunction.
 *   - Called cuFuncGetAttribute to get max threads per block.
 *   - Iterated over the arg-type list computing byte offsets for the
 *     manual cuParamSet* packing scheme.
 *   - Called cuParamSetSize to commit the total parameter buffer size.
 *
 * OpenCL version:
 *   - Calls clCreateKernel to get a cl_kernel.
 *   - Queries CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE for the
 *     preferred local work size (analogous to threads_per_block hint).
 *   - No offset arithmetic needed; clSetKernelArg handles everything.
 * --------------------------------------------------------------------- */
void
gpu_launch_init(cl_program program, const char *func_name,
                const gpu_arg_type_list_t *arg_desc,
                gpu_launch_t *launch,
                cl_device_id device)
{
    cl_int  err;
    size_t  wg_mult = 0;

    memset(launch, 0, sizeof(gpu_launch_t));

    launch->kernel_func = clCreateKernel(program, func_name, &err);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "gpu_launch_init: clCreateKernel(%s) failed: %s\n",
                func_name, clGetErrorString(err));
        exit(-1);
    }

    /* Query preferred work-group size multiple (≈ warp size).
     * This is the closest OpenCL analogue to
     * CU_FUNC_ATTRIBUTE_MAX_THREADS_PER_BLOCK. */
    clGetKernelWorkGroupInfo(launch->kernel_func, device,
                             CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
                             sizeof(wg_mult), &wg_mult, NULL);
    launch->threads_per_block = (int32_t)wg_mult;

    launch->arg_desc = *arg_desc;
    /* No offset table to fill -- clSetKernelArg is index-based. */
}

/* -----------------------------------------------------------------------
 * gpu_launch_set
 *
 * CUDA version used cuParamSeti / cuParamSetv with pre-computed byte
 * offsets.  OpenCL replaces that with one clSetKernelArg call per
 * argument, indexed by position.
 *
 * Argument type mapping:
 *   GPU_ARG_PTR    -> cl_mem  (the actual device buffer object)
 *   GPU_ARG_INT32  -> int32_t
 *   GPU_ARG_UINT32 -> uint32_t
 *   GPU_ARG_INT64  -> int64_t
 *   GPU_ARG_UINT64 -> uint64_t
 * --------------------------------------------------------------------- */
void
gpu_launch_set(gpu_launch_t *launch, gpu_arg_t *args)
{
    uint32_t i;
    cl_int   err;

    for (i = 0; i < launch->arg_desc.num_args; i++) {
        switch (launch->arg_desc.arg_type[i]) {

        case GPU_ARG_PTR:
            /* args[i].ptr_arg is now a cl_mem, not a raw pointer */
            err = clSetKernelArg(launch->kernel_func, i,
                                 sizeof(cl_mem), &args[i].ptr_arg);
            break;

        case GPU_ARG_INT32:
            err = clSetKernelArg(launch->kernel_func, i,
                                 sizeof(int32_t), &args[i].int32_arg);
            break;

        case GPU_ARG_UINT32:
            err = clSetKernelArg(launch->kernel_func, i,
                                 sizeof(uint32_t), &args[i].uint32_arg);
            break;

        case GPU_ARG_INT64:
            err = clSetKernelArg(launch->kernel_func, i,
                                 sizeof(int64_t), &args[i].int64_arg);
            break;

        case GPU_ARG_UINT64:
            err = clSetKernelArg(launch->kernel_func, i,
                                 sizeof(uint64_t), &args[i].uint64_arg);
            break;

        default:
            fprintf(stderr, "gpu_launch_set: unknown GPU argument type\n");
            exit(-1);
        }

        if (err != CL_SUCCESS) {
            fprintf(stderr, "gpu_launch_set: clSetKernelArg(%u) failed: %s\n",
                    i, clGetErrorString(err));
            exit(-1);
        }
    }
}

#endif /* HAVE_CUDA_BATCH_FACTOR */
