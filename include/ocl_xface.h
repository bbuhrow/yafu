/*--------------------------------------------------------------------
 * ocl_xface.h  --  OpenCL replacement for cuda_xface.h
 *
 * Translation notes
 * -----------------
 * CUDA Driver API concept        -> OpenCL equivalent
 * ---------------------------       --------------------------
 * CUdevice                       -> cl_device_id
 * CUcontext                      -> cl_context
 * CUmodule                       -> cl_program
 * CUfunction                     -> cl_kernel
 * CUstream                       -> cl_command_queue
 * CUevent                        -> cl_event  (for timing)
 * CUdeviceptr                    -> cl_mem
 * cuMemAlloc / cuMemFree         -> clCreateBuffer / clReleaseMemObject
 * cuMemcpyHtoDAsync              -> clEnqueueWriteBuffer
 * cuMemcpyDtoHAsync              -> clEnqueueReadBuffer
 * cuStreamSynchronize            -> clFinish
 * cuEventRecord/Synchronize      -> CL_QUEUE_PROFILING_ENABLE + cl_event
 * cuEventElapsedTime             -> clGetEventProfilingInfo
 * cuFuncSetBlockShape +          -> local_work_size / global_work_size
 *   cuLaunchGridAsync               passed to clEnqueueNDRangeKernel
 * cuParamSet* / cuParamSetSize   -> clSetKernelArg  (one call per arg)
 *
 * The whole gpu_launch_t / gpu_arg_type_list_t machinery in cuda_xface
 * existed to work around CUDA's old manual parameter-packing API.
 * In OpenCL each argument is set individually with clSetKernelArg, so
 * the arg-type list is simplified: we still keep the structure so that
 * the call sites in gpu_cofactorization.c can remain largely unchanged,
 * but it no longer needs offset arithmetic.
 *--------------------------------------------------------------------*/

#ifndef _OCL_XFACE_H
#define _OCL_XFACE_H

#ifdef HAVE_OCL_BATCH_FACTOR   /* reuse the existing guard flag */

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifdef __APPLE__
#  include <OpenCL/opencl.h>
#else
#  include <CL/cl.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* -----------------------------------------------------------------------
 * Error handling
 * --------------------------------------------------------------------- */

const char *clGetErrorString(cl_int err);

/* Drop-in replacement for CUDA_TRY.
 * Evaluates the OpenCL call, checks cl_int status, prints and exits on
 * error.  Wrap any cl* call that returns cl_int. */
#define OCL_TRY(func) \
    { \
        cl_int _status = (func); \
        if (_status != CL_SUCCESS) { \
            fprintf(stderr, "OpenCL error (line %d): %s\n", \
                    __LINE__, clGetErrorString(_status)); \
            exit(-1); \
        } \
    }

/* Variant for calls that return the object and take cl_int* as last arg
 * (e.g. clCreateBuffer, clCreateCommandQueue).  Use as:
 *   buf = OCL_CREATE(clCreateBuffer(ctx, flags, sz, NULL, &_err), _err); */
#define OCL_CREATE(expr, err_var) \
    ( (err_var) = CL_SUCCESS, (expr) )

#define MAX_GPU 4

/* -----------------------------------------------------------------------
 * Device / platform info
 * --------------------------------------------------------------------- */

typedef struct {
    char            name[128];
    int32_t         compute_version_major;  /* OpenCL major version   */
    int32_t         compute_version_minor;  /* OpenCL minor version   */
    int32_t         clock_speed;            /* in MHz                 */
    int32_t         num_compute_units;
    int32_t         constant_mem_size;      /* bytes                  */
    int32_t         shared_mem_size;        /* local mem per WG       */
    size_t          global_mem_size;
    int32_t         registers_per_block;    /* not exposed in OpenCL; set 0 */
    int32_t         max_threads_per_block;  /* max work-group size    */
    int32_t         can_overlap;            /* always 1 in OpenCL     */
    int32_t         warp_size;              /* preferred WG size multiple */
    int32_t         max_thread_dim[3];      /* max WG dims            */
    int32_t         max_grid_size[3];       /* max NDRange dims       */
    int32_t         has_timeout;            /* 0 for compute devices  */
    cl_device_id    device_handle;
    cl_platform_id  platform_handle;
} gpu_info_t;

typedef struct {
    int32_t     num_gpu;
    gpu_info_t  info[MAX_GPU];
} gpu_config_t;

void gpu_init(gpu_config_t *config);

/* -----------------------------------------------------------------------
 * Kernel argument descriptor
 * (kept structurally identical to the CUDA version so call-sites compile
 *  unchanged; the actual parameter-setting logic moves to gpu_launch_set)
 * --------------------------------------------------------------------- */

typedef enum {
    GPU_ARG_NONE   = 0,
    GPU_ARG_PTR,        /* cl_mem   */
    GPU_ARG_INT32,
    GPU_ARG_UINT32,
    GPU_ARG_INT64,
    GPU_ARG_UINT64
} gpu_arg_type_t;

#define GPU_MAX_KERNEL_ARGS 15

typedef struct {
    uint32_t        num_args;
    gpu_arg_type_t  arg_type[GPU_MAX_KERNEL_ARGS];
} gpu_arg_type_list_t;

typedef union {
    cl_mem      ptr_arg;    /* NOTE: was void* in CUDA; now cl_mem       */
    int32_t     int32_arg;
    uint32_t    uint32_arg;
    int64_t     int64_arg;
    uint64_t    uint64_arg;
} gpu_arg_t;

/* -----------------------------------------------------------------------
 * Launch descriptor
 *
 * CUDA stored pre-computed byte offsets for each argument and used
 * cuParamSetv/cuParamSeti to pack them into a buffer.  In OpenCL we call
 * clSetKernelArg(kernel, idx, size, &value) directly -- no offsets needed.
 * The arg_offsets array is removed; everything else stays the same shape.
 * --------------------------------------------------------------------- */

typedef struct {
    cl_kernel           kernel_func;
    int32_t             threads_per_block;  /* preferred local WG size    */
    gpu_arg_type_list_t arg_desc;
} gpu_launch_t;

/* Initialise a gpu_launch_t: look up the named kernel in the program,
 * query its preferred work-group size, and store the arg descriptor. */
void gpu_launch_init(cl_program program, const char *func_name,
                     const gpu_arg_type_list_t *arg_desc,
                     gpu_launch_t *launch,
                     cl_device_id device);

/* Set all kernel arguments from the args array. */
void gpu_launch_set(gpu_launch_t *launch, gpu_arg_t *args);

#ifdef __cplusplus
}
#endif

#endif /* HAVE_OCL_BATCH_FACTOR */
#endif /* _OCL_XFACE_H */
