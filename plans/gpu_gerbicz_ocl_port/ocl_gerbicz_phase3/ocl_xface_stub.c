/*--------------------------------------------------------------------
 * ocl_xface_stub.c
 *
 * NOT part of the Phase 1 deliverable proper. gpu_init(), gpu_launch_init()
 * and gpu_launch_set() are declared in ocl_xface.h and already called by
 * gpu_cofactorization_cl.c, so an implementation must exist somewhere in
 * the real repo already -- it simply wasn't among the files handed to
 * this phase. This is a straightforward, standard-OpenCL reference
 * implementation used ONLY so Phase 1's smoke test (task 1.7) can
 * actually build and run in this sandbox. See STATUS "Open issues":
 * reconcile/replace with the real one before Phase 3 relies on it.
 *--------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ocl_shared.h"   /* pulls in ocl_xface.h + clGetErrorString() decl */

void
gpu_init(gpu_config_t *config)
{
    cl_platform_id platforms[8];
    cl_uint num_platforms = 0;
    cl_uint p;

    memset(config, 0, sizeof(*config));

    if (clGetPlatformIDs(8, platforms, &num_platforms) != CL_SUCCESS)
        return;

    for (p = 0; p < num_platforms && config->num_gpu < MAX_GPU; p++) {
        cl_device_id devices[MAX_GPU];
        cl_uint num_devices = 0;
        cl_uint i;

        /* CL_DEVICE_TYPE_ALL so a CPU-only OpenCL runtime (e.g. PoCL, used
         * in this sandbox in place of the real AMD RX 6700 XT) is still
         * found; production code targeting a real GPU can narrow this. */
        if (clGetDeviceIDs(platforms[p], CL_DEVICE_TYPE_ALL,
                MAX_GPU - config->num_gpu, devices, &num_devices) != CL_SUCCESS)
            continue;

        for (i = 0; i < num_devices && config->num_gpu < MAX_GPU; i++) {
            gpu_info_t *gi = &config->info[config->num_gpu];
            cl_device_id d = devices[i];
            cl_uint cu = 0, clk = 0;
            cl_ulong cmem = 0, gmem = 0;
            cl_ulong lmem = 0;
            size_t mwgs = 0;
            size_t maxdim[3] = {0,0,0};
            cl_uint maxdims = 0;
            char ver[64];

            memset(gi, 0, sizeof(*gi));
            clGetDeviceInfo(d, CL_DEVICE_NAME, sizeof(gi->name), gi->name, NULL);
            clGetDeviceInfo(d, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cu), &cu, NULL);
            clGetDeviceInfo(d, CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(clk), &clk, NULL);
            clGetDeviceInfo(d, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(cmem), &cmem, NULL);
            clGetDeviceInfo(d, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(lmem), &lmem, NULL);
            clGetDeviceInfo(d, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(gmem), &gmem, NULL);
            clGetDeviceInfo(d, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(mwgs), &mwgs, NULL);
            clGetDeviceInfo(d, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(maxdims), &maxdims, NULL);
            clGetDeviceInfo(d, CL_DEVICE_VERSION, sizeof(ver), ver, NULL);
            {
                size_t wisz_ret = 0;
                size_t wisz[3] = {0,0,0};
                clGetDeviceInfo(d, CL_DEVICE_MAX_WORK_ITEM_SIZES,
                        sizeof(wisz), wisz, &wisz_ret);
                maxdim[0] = wisz[0]; maxdim[1] = wisz[1]; maxdim[2] = wisz[2];
            }

            gi->compute_version_major = 1;
            gi->compute_version_minor = 2;
            {
                int maj = 1, min = 2;
                sscanf(ver, "OpenCL %d.%d", &maj, &min);
                gi->compute_version_major = maj;
                gi->compute_version_minor = min;
            }
            gi->clock_speed = (int32_t)clk;
            gi->num_compute_units = (int32_t)cu;
            gi->constant_mem_size = (int32_t)cmem;
            gi->shared_mem_size = (int32_t)lmem;
            gi->global_mem_size = (size_t)gmem;
            gi->registers_per_block = 0;
            gi->max_threads_per_block = (int32_t)mwgs;
            gi->can_overlap = 1;
            {
                size_t pwgs = 0;
                clGetDeviceInfo(d, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(pwgs), &pwgs, NULL);
                gi->warp_size = 32; /* placeholder; real value needs a built kernel */
            }
            gi->max_thread_dim[0] = (int32_t)maxdim[0];
            gi->max_thread_dim[1] = (int32_t)maxdim[1];
            gi->max_thread_dim[2] = (int32_t)maxdim[2];
            gi->max_grid_size[0] = gi->max_thread_dim[0];
            gi->max_grid_size[1] = gi->max_thread_dim[1];
            gi->max_grid_size[2] = gi->max_thread_dim[2];
            gi->has_timeout = 0;
            gi->device_handle = d;
            gi->platform_handle = platforms[p];

            config->num_gpu++;
        }
    }
}

void
gpu_launch_init(cl_program program, const char *func_name,
                 const gpu_arg_type_list_t *arg_desc,
                 gpu_launch_t *launch, cl_device_id device)
{
    cl_int err;
    size_t pref = 0;

    memset(launch, 0, sizeof(*launch));
    launch->kernel_func = clCreateKernel(program, func_name, &err);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "gpu_launch_init: clCreateKernel('%s') failed: %s\n",
                func_name, clGetErrorString(err));
        exit(-1);
    }

    if (clGetKernelWorkGroupInfo(launch->kernel_func, device,
            CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
            sizeof(pref), &pref, NULL) == CL_SUCCESS && pref > 0) {
        launch->threads_per_block = (int32_t)pref;
    } else {
        launch->threads_per_block = 64;
    }

    launch->arg_desc = *arg_desc;
}

void
gpu_launch_set(gpu_launch_t *launch, gpu_arg_t *args)
{
    uint32_t i;
    for (i = 0; i < launch->arg_desc.num_args; i++) {
        cl_int err = CL_SUCCESS;
        switch (launch->arg_desc.arg_type[i]) {
        case GPU_ARG_PTR:
            err = clSetKernelArg(launch->kernel_func, i, sizeof(cl_mem),
                                  &args[i].ptr_arg);
            break;
        case GPU_ARG_INT32:
            err = clSetKernelArg(launch->kernel_func, i, sizeof(int32_t),
                                  &args[i].int32_arg);
            break;
        case GPU_ARG_UINT32:
            err = clSetKernelArg(launch->kernel_func, i, sizeof(uint32_t),
                                  &args[i].uint32_arg);
            break;
        case GPU_ARG_INT64:
            err = clSetKernelArg(launch->kernel_func, i, sizeof(int64_t),
                                  &args[i].int64_arg);
            break;
        case GPU_ARG_UINT64:
            err = clSetKernelArg(launch->kernel_func, i, sizeof(uint64_t),
                                  &args[i].uint64_arg);
            break;
        case GPU_ARG_LOCAL:
            /* Dynamically-sized __local argument (Phase 3 addition --
             * see ocl_xface_arg_local.h.patch). args[i].uint32_arg holds
             * the byte size; NULL host pointer means "allocate in local
             * memory, don't initialize from host". */
            err = clSetKernelArg(launch->kernel_func, i,
                                  (size_t)args[i].uint32_arg, NULL);
            break;
        default:
            continue;
        }
        if (err != CL_SUCCESS) {
            fprintf(stderr, "gpu_launch_set: clSetKernelArg(%u) failed: %s\n",
                    i, clGetErrorString(err));
            exit(-1);
        }
    }
}
