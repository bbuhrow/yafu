/*--------------------------------------------------------------------
 * gpu_cofactorization.h  --  OpenCL translation
 *
 * All CUxxx types in device_thread_ctx_t are replaced with their
 * OpenCL equivalents:
 *
 *   CUcontext    -> cl_context
 *   CUmodule     -> cl_program
 *   CUstream     -> cl_command_queue  (queues are the stream equivalent)
 *   CUevent      -> cl_event          (timing events)
 *   CUdeviceptr  -> cl_mem            (opaque device buffer handle)
 *
 * The include of cuda_xface.h is replaced with ocl_xface.h.
 *--------------------------------------------------------------------*/

#pragma once
#include <stdint.h>
#include "batch_factor.h"
#include "ocl_xface.h"
#include "cofactorize.h"

#ifdef HAVE_OCL_BATCH_FACTOR

typedef struct {
    int         gpunum;
    gpu_info_t *gpu_info;
} device_ctx_t;

typedef struct {
    /* reference to the configured device */
    device_ctx_t *dev;

    /* config from calling code */
    uint32_t lpba;          /* lpb on the a-side                    */
    uint32_t lpbr;          /* lpb on the r-side                    */
    uint32_t mfba;          /* to determine if 3LP or 2LP           */
    uint32_t mfbr;          /* to determine if 3LP or 2LP           */
    int      verbose;

    /* internal config */
    int      lpb_2lp;
    int      lpb_3lp;
    int      first_side;    /* first side is the 2LP side           */
    uint32_t b1_2lp;
    uint32_t b2_2lp;
    uint32_t curves_2lp;
    uint32_t b1_3lp;
    uint32_t b2_3lp;
    uint32_t curves_3lp;
    uint32_t stop_nofactor;
    tiny_qs_params* params;

    /* ----------------------------------------------------------------
     * OpenCL runtime objects
     * (replaces: CUcontext, CUmodule, CUstream, CUevent x2)
     * -------------------------------------------------------------- */
    cl_context       gpu_context;   /* was: CUcontext gpu_context    */
    cl_program       gpu_program;   /* was: CUmodule  gpu_module     */
    cl_command_queue queue;         /* was: CUstream  stream         */
    gpu_launch_t    *launch;
    cl_event         start_event;   /* was: CUevent   start_event    */
    cl_event         end_event;     /* was: CUevent   end_event      */

    /* the data structure to work on */
    relation_batch_t *rb;

    /* ----------------------------------------------------------------
     * Device (GPU) buffers
     * (replaces CUdeviceptr -- now cl_mem opaque handles)
     * -------------------------------------------------------------- */
    cl_mem  gpu_a_array;        /* was: CUdeviceptr gpu_a_array      */
    cl_mem  gpu_b_array;        /* was: CUdeviceptr gpu_b_array      */
    cl_mem  gpu_n_array;        /* was: CUdeviceptr gpu_n_array      */
    cl_mem  gpu_one_array;      /* was: CUdeviceptr gpu_one_array    */
    cl_mem  gpu_rsq_array;      /* was: CUdeviceptr gpu_rsq_array    */
    cl_mem  gpu_u32_array;      /* was: CUdeviceptr gpu_u32_array    */
    cl_mem  gpu_res32_array;    /* was: CUdeviceptr gpu_res32_array  */
    cl_mem  gpu_rho_array;      /* was: CUdeviceptr gpu_rho_array    */

    uint32_t  array_sz;
    uint32_t *u32_array;
    uint32_t  flags;

    /* host data */
    uint64_t *a;
    uint64_t *b;

    /* common */
    uint32_t *rho;

    /* for ecm64 */
    uint64_t *rsq;
    uint64_t *one;
    uint64_t *modulus_in;
    uint64_t *sigma_in;

    /* for ecm96 */
    uint32_t *rsq96;
    uint32_t *one96;
    uint32_t *modulus96_in;
    uint32_t *factors96;

    /* tracking factors to input residues */
    uint32_t *residues_r_in;
    uint32_t *residues_a_in;
    uint32_t  numres_r;
    uint32_t  numres_a;
    uint32_t *rb_idx_r;
    uint32_t *rb_idx_a;
    uint32_t *rb_idx_2lp;
    uint32_t *rb_idx_3lp;
    uint32_t  num_factors_2lp;
    uint32_t  num_factors_3lp;

    int mode_2lp;   /* 0: 2LP factorizations.  1: 3LP-cofactorizations */

} device_thread_ctx_t;

/* create and destroy gpu context info */
device_ctx_t        *gpu_device_init(int which_gpu, int verbose);
device_thread_ctx_t *gpu_ctx_init(device_ctx_t *d);
void                 gpu_dev_free(device_ctx_t *d);
void                 gpu_ctx_free(device_thread_ctx_t *t);

/* do gpu cofactorization work */
int do_gpu_cofactorization(device_thread_ctx_t *t, relation_batch_t* rb, uint64_t *lcg,
                           int b1_3lp_ovr, int b2_3lp_ovr,
                           int b1_2lp_ovr, int b2_2lp_ovr,
                           int curves_3lp_ovr, int curves_2lp_ovr);

#endif /* HAVE_CUDA_BATCH_FACTOR */
