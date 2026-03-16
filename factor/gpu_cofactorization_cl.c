/*--------------------------------------------------------------------
 * gpu_cofactorization.c  --  OpenCL translation
 *
 * Per-call translation summary
 * ----------------------------
 *
 * Memory transfers (all were Async; OpenCL uses blocking=CL_FALSE for
 * the same effect with an in-order queue):
 *   cuMemcpyHtoDAsync(dst, src, sz, stream)
 *     -> clEnqueueWriteBuffer(queue, dst, CL_FALSE, 0, sz, src, 0,NULL,NULL)
 *   cuMemcpyDtoHAsync(dst, src, sz, stream)
 *     -> clEnqueueReadBuffer (queue, src, CL_FALSE, 0, sz, dst, 0,NULL,NULL)
 *
 * Kernel launch (CUDA split this across two calls):
 *   cuFuncSetBlockShape(func, tpb, 1, 1)   -- sets local size
 *   cuLaunchGridAsync(func, blocks, 1, stream) -- launches
 *     -> clEnqueueNDRangeKernel(queue, kernel, 1, NULL,
 *                               &global_size, &local_size, 0,NULL,NULL)
 *   where global_size = num_blocks * threads_per_block
 *         local_size  = threads_per_block
 *
 * Synchronization:
 *   cuStreamSynchronize(stream) -> clFinish(queue)
 *
 * Timing:
 *   cuEventRecord(ev, stream)      \
 *   cuEventSynchronize(ev)          >  replaced with wall-clock timing
 *   cuEventElapsedTime(&ms, s, e)  /   via clock_gettime or a profiling
 *                                      event approach (see below).
 *   OpenCL timing uses CL_QUEUE_PROFILING_ENABLE on the queue, then
 *   clGetEventProfilingInfo after clFinish.  Because the existing code
 *   wraps an entire multi-iteration loop between start/end events, we
 *   use simple wall-clock time (clock_gettime) for the overall elapsed
 *   measurement -- it is equivalent and simpler.
 *
 * Allocation / free:
 *   cuMemAlloc(&ptr, sz)  -> ptr = clCreateBuffer(ctx, CL_MEM_READ_WRITE, sz, NULL, &err)
 *   cuMemFree(ptr)        -> clReleaseMemObject(ptr)
 *
 * Context / module / stream / events:
 *   cuCtxCreate  -> clCreateContext (see gpu_ctx_init)
 *   cuModuleLoad -> clCreateProgramWithSource + clBuildProgram
 *   cuStreamCreate -> clCreateCommandQueue (in-order, no special flags)
 *   cuEventCreate  -> not needed as objects; we use wall-clock timing
 *   cuEventDestroy -> not needed
 *   cuStreamDestroy -> clReleaseCommandQueue
 *   cuCtxDestroy    -> clReleaseContext + clReleaseProgram
 *
 * gpu_module rename:
 *   t->gpu_module  (CUmodule) -> t->gpu_program  (cl_program)
 *
 * ptr_arg in gpu_arg_t:
 *   Was void* (cast from CUdeviceptr).  Is now cl_mem directly.
 *   Call sites change from:
 *     gpu_args[N].ptr_arg = (void*)(t->gpu_xxx_array);
 *   to:
 *     gpu_args[N].ptr_arg = t->gpu_xxx_array;
 *--------------------------------------------------------------------*/

// mpqs builds in msys2 but crashes right away when run.
// something to figure out later.
#if defined(__GNUC__) && !defined(__MINGW32__)
#define HAVE_LASIEVE_MPQS
#endif

// posix feature-test macros to enable things in posix headers.
// this one is for clock_gettime and others in time.h
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <string.h>
#include <time.h>
#include "gmp.h"
#include "microecm.h"
#include "batch_factor.h"
#include "gpu_cofactorization_cl.h"
#include "gmp-aux.h"
#ifdef HAVE_LASIEVE_MPQS
#include "mpqs3/mpqs.h"
#include "mpqs3/mpqs3.h"
#include "mpqs3/if.h"
#ifdef ULL_NO_UL
#include "mpqs3/gmp-aux.h"
#endif
#endif
#include "cofactorize.h"

//#define HAVE_CUDA_BATCH_FACTOR

#ifdef HAVE_OCL_BATCH_FACTOR

#define MAX_RESIDUE_WORDS 3

enum test_flags {
    DEFAULT_FLAGS       = 0,
    FLAG_USE_LOGFILE    = 0x01,
    FLAG_LOG_TO_STDOUT  = 0x02,
    FLAG_STOP_GRACEFULLY= 0x04
};

enum {
    GPU_ECM_VEC   = 0,
    GPU_ECM96_VEC,
    GPU_PM196_VEC,
    NUM_GPU_FUNCTIONS
};

static const char *gpu_kernel_names[] = {
    "gbl_ecm",
    "gbl_ecm96",
    "gbl_pm196",
};

static const gpu_arg_type_list_t gpu_kernel_args[] = {
    /* ecm -- 9 args */
    { 9, { GPU_ARG_INT32, GPU_ARG_PTR, GPU_ARG_PTR, GPU_ARG_PTR,
           GPU_ARG_PTR,   GPU_ARG_PTR, GPU_ARG_PTR, GPU_ARG_UINT32,
           GPU_ARG_INT32 } },
    /* ecm96 -- 10 args */
    { 10, { GPU_ARG_INT32,  GPU_ARG_PTR, GPU_ARG_PTR, GPU_ARG_PTR,
            GPU_ARG_PTR,    GPU_ARG_PTR, GPU_ARG_PTR, GPU_ARG_UINT32,
            GPU_ARG_UINT32, GPU_ARG_INT32 } },
    /* pm196 -- 7 args */
    { 7, { GPU_ARG_INT32, GPU_ARG_PTR, GPU_ARG_PTR, GPU_ARG_PTR,
           GPU_ARG_PTR,   GPU_ARG_UINT32, GPU_ARG_UINT32 } },
};



/* -----------------------------------------------------------------------
 * Convenience: simple wall-clock elapsed time in ms between two
 * struct timespec values.  Replaces cuEventElapsedTime.
 * --------------------------------------------------------------------- */
static float
timespec_elapsed_ms(struct timespec *start, struct timespec *end)
{
    double s = (double)(end->tv_sec  - start->tv_sec)  * 1000.0
             + (double)(end->tv_nsec - start->tv_nsec) / 1e6;
    return (float)s;
}


#define BIGWORDS32 3

void mpz_to_bignum32(uint32_t* bignum, mpz_t gmp_in)
{
    int i;
    mpz_t t;
    mpz_init(t);
    mpz_set(t, gmp_in);

    for (i = 0; i < BIGWORDS32; i++)
    {
        bignum[i] = mpz_get_ui(t) & 0xffffffff;
        mpz_tdiv_q_2exp(t, t, 32);
    }

    mpz_clear(t);
    return;
}

void bignum32_to_mpz(mpz_t gmp_out, uint32_t* bignum)
{
    int i;

    mpz_set_ui(gmp_out, bignum[BIGWORDS32 - 1]);
    for (i = BIGWORDS32 - 2; i >= 0; i--)
    {
        mpz_mul_2exp(gmp_out, gmp_out, 32);
        mpz_add_ui(gmp_out, gmp_out, bignum[i]);
    }

    return;
}

static const double INV_2_POW_32 = 1.0 / (double)((uint64_t)(1) << 32);
static uint32_t uecm_lcg_rand_32B(uint32_t lower, uint32_t upper, uint64_t* ploc_lcg)
{
    *ploc_lcg = 6364136223846793005ULL * (*ploc_lcg) + 1442695040888963407ULL;
    return lower + (uint32_t)(
        (double)(upper - lower) * (double)((*ploc_lcg) >> 32) * INV_2_POW_32);
}

uint32_t multiplicative_neg_inverse32(uint64_t a)
{
    uint32_t res = 2 + a;
    res = res * (2 + a * res);
    res = res * (2 + a * res);
    res = res * (2 + a * res);
    return res * (2 + a * res);
}

int handle_96b_factorization(device_thread_ctx_t* t, int idx,
    mpz_t zf, mpz_t zc, mpz_t zn, mpz_t* flist,
    int num2lp_retest, int* mpqs_success, int* mpqs_err, int* num_mpqs)
{
    // here is a non-trivial factor that divides the modulus.
    // check it against LPB size constraint.
    int bits1 = mpz_sizeinbase(zf, 2);
    mpz_tdiv_q(zc, zn, zf);
    int bits2 = mpz_sizeinbase(zc, 2);

    if (bits1 <= t->lpb_3lp)
    {
        // the factor we found is good.
        // the cofactor needs further analysis
        // that we either assign to a list for
        // more gpu-ecm work, or tackle immediately
        // if too big for that.

        // check if cofactor is small enough to re-process with 64-bit kernel.
        if (bits2 <= 64)
        {
            cofactor_t* c = t->rb->relations + t->rb_idx_3lp[idx];
            uint64_t cofactor = mpz_get_ull(zc);

            // check if cofactor is composite and small enough to
            // possibly yield 2 correctly-sized primes
            //if (prp_uecm(cofactor) == 0)
            if ((bits2 <= 2 * t->lpb_3lp) && (mpz_probab_prime_p(zc, 1) == 0))
            {
                // record the factor we found
                if (t->first_side == 0)
                {
                    c->lp_a[0] = mpz_get_ull(zf);
                }
                else
                {
                    c->lp_r[0] = mpz_get_ull(zf);
                }
                // and load the cofactor for further factorization
                t->modulus_in[num2lp_retest] = cofactor;
                t->rb_idx_2lp[num2lp_retest] = t->rb_idx_3lp[idx];
                num2lp_retest++;
            }
            else if (bits2 > 2 * t->lpb_3lp)
            {
                // this cofactor is too big
                c->success = 0;
            }
            else if (bits2 <= t->lpb_3lp)
            {
                // we just factored a 2LP larger than 64 bits.
                if (t->first_side == 0)
                {
                    c->lp_a[0] = mpz_get_ull(zf);
                    c->lp_a[1] = cofactor;
                }
                else
                {
                    c->lp_r[0] = mpz_get_ull(zf);
                    c->lp_r[1] = cofactor;
                }
                c->success |= 0x0f;
            }
        }
        else
        {
            // the cofactor is larger than 64 bits. if it's 
            // not prime we could either try to factor it here
            // (slow) or do another 96-bit pass on it (more code).
            // for now try to do it here with mpqs.
            
            // check if cofactor is composite and small enough to
            // possibly yield 2 correctly-sized primes
            //if (prp_uecm(cofactor) == 0)
            if ((bits2 <= 2 * t->lpb_3lp) && (mpz_probab_prime_p(zc, 1) == 0))
            {
#ifdef HAVE_LASIEVE_MPQS
                int nf = mpqs_factor(zc, t->lpb_3lp, &flist);
                (*num_mpqs)++;
#else
                (*num_mpqs)++;
                int nf = tinysiqs(t->params, zc, flist[0], flist[1], flist[2], t->lpb_3lp);
#endif
                if (nf == 2)
                {
                    //gmp_printf("factored %Zd as %Zd * %Zd by tinysiqs\n", zc, flist[0], flist[1]);
                    (*mpqs_success)++;
                    cofactor_t* c = t->rb->relations + t->rb_idx_3lp[idx];
                    if (t->first_side == 0)
                    {
                        c->lp_a[0] = mpz_get_ull(zf);
                        c->lp_a[1] = mpz_get_ull(flist[0]);
                        c->lp_a[2] = mpz_get_ull(flist[1]);
                    }
                    else
                    {
                        c->lp_r[0] = mpz_get_ull(zf);
                        c->lp_r[1] = mpz_get_ull(flist[0]);
                        c->lp_r[2] = mpz_get_ull(flist[1]);
                    }
                    c->success |= 0x0f;
                }
                else if (nf == 0)
                {
                    //gmp_printf("tinysiqs failed to factor %Zd, found %d factors\n", zc, nf);
                    // do some rho?
                    // submit it back for more 96-bit ecm?
                }
                else if (nf > 5)
                {
                    (*mpqs_err)++;
                }
            }
            else
            {
                cofactor_t* c = t->rb->relations + t->rb_idx_3lp[idx];
                // if it's too big or prime, it's no good
                c->success = 0;
            }

        }
    }
    else if (bits2 <= t->lpb_3lp)
    {
        // we found either an improbably large prime factor
        // or two smaller factors simultaneously.
        // submit this for further gpu analysis.
        // build up a list on which we'll do 64-bit
        // factorizations as needed.  possible 
        // to maybe also do the prp checks on gpu
        // but these are extremely cheap on cpu.

        // check cofactor
        if (bits1 <= 64)
        {
            cofactor_t* c = t->rb->relations + t->rb_idx_3lp[idx];
            uint64_t cofactor = mpz_get_ull(zf);
            // check if cofactor is composite and small enough to
            // possibly yield 2 correctly-sized primes
            //if (prp_uecm(cofactor) == 0)
            if ((bits1 <= 2 * t->lpb_3lp) && (mpz_probab_prime_p(zf, 1) == 0))
            {
                // record the factor we found
                if (t->first_side == 0)
                {
                    c->lp_a[0] = mpz_get_ull(zc);
                }
                else
                {
                    c->lp_r[0] = mpz_get_ull(zc);
                }

                // and load the cofactor for further factorization
                t->modulus_in[num2lp_retest] = cofactor;
                t->rb_idx_2lp[num2lp_retest] = t->rb_idx_3lp[idx];
                num2lp_retest++;
            }
            else if (bits2 > 2 * t->lpb_3lp)
            {
                // this cofactor is too big
                c->success = 0;
            }
            else if (bits1 <= t->lpb_3lp)
            {
                // we just factored a 2LP larger than 64 bits.
                if (t->first_side == 0)
                {
                    c->lp_a[0] = mpz_get_ull(zc);
                    c->lp_a[1] = cofactor;
                }
                else
                {
                    c->lp_r[0] = mpz_get_ull(zc);
                    c->lp_r[1] = cofactor;
                }
                c->success |= 0x0f;
            }
        }
        else
        {
            // the cofactor is larger than 64 bits, if it's 
            // not prime we could either try to factor it here
            // (slow) or do another 96-bit pass on it (more code)
            // for now try to do it here with mpqs.
            // 
            // check if cofactor is composite and small enough to
            // possibly yield 2 correctly-sized primes
            //if (prp_uecm(cofactor) == 0)
            if ((bits1 <= 2 * t->lpb_3lp) && (mpz_probab_prime_p(zf, 1) == 0))
            {
#ifdef HAVE_LASIEVE_MPQS
                int nf = mpqs_factor(zf, t->lpb_3lp, &flist);
                (*num_mpqs)++;
#else
                (*num_mpqs)++;
                int nf = tinysiqs(t->params, zf, flist[0], flist[1], flist[2], t->lpb_3lp);
#endif

                if (nf == 2)
                {
                    //gmp_printf("factored %Zd as %Zd * %Zd by tinysiqs\n", zf, flist[0], flist[1]);
                    cofactor_t* c = t->rb->relations + t->rb_idx_3lp[idx];
                    (*mpqs_success)++;
                    if (t->first_side == 0)
                    {
                        c->lp_a[0] = mpz_get_ull(zc);
                        c->lp_a[1] = mpz_get_ull(flist[0]);
                        c->lp_a[2] = mpz_get_ull(flist[1]);
                    }
                    else
                    {
                        c->lp_r[0] = mpz_get_ull(zc);
                        c->lp_r[1] = mpz_get_ull(flist[0]);
                        c->lp_r[2] = mpz_get_ull(flist[1]);
                    }
                    c->success |= 0x0f;
                }
                else if (nf == 0)
                {
                    //gmp_printf("tinysiqs failed to factor %Zd, found %d factors\n", zf, nf);
                }
                else if (nf > 5)
                {
                    (*mpqs_err)++;
                }
            }
            else
            {
                cofactor_t* c = t->rb->relations + t->rb_idx_3lp[idx];
                // if it's too big or prime, it's no good
                c->success = 0;
            }
        }
    }
    return num2lp_retest;
}


/* -----------------------------------------------------------------------
 * do_gpu_ecm64
 * --------------------------------------------------------------------- */
uint32_t
do_gpu_ecm64(device_thread_ctx_t *t)
{
    uint32_t     quit = 0;
    gpu_arg_t    gpu_args[GPU_MAX_KERNEL_ARGS];
    gpu_launch_t *launch;
    struct timespec ts_start, ts_end;
    float        elapsed_ms;

    int threads_per_block = 128;
    int num_blocks = t->array_sz / threads_per_block
                   + ((t->array_sz % threads_per_block) > 0);

    printf("commencing gpu 64-bit ecm in mode %d on %d inputs\n",
           t->mode_2lp, t->array_sz);
    fflush(stdout);

    /* copy sigma into device memory */
    OCL_TRY(clEnqueueWriteBuffer(t->queue, t->gpu_u32_array, CL_FALSE,
        0, t->array_sz * sizeof(uint32_t), t->u32_array, 0, NULL, NULL))

    /* copy n into device memory */
    OCL_TRY(clEnqueueWriteBuffer(t->queue, t->gpu_n_array, CL_FALSE,
        0, t->array_sz * sizeof(uint64_t), t->modulus_in, 0, NULL, NULL))

    clock_gettime(CLOCK_MONOTONIC, &ts_start);  /* replaces cuEventRecord */

    int curve         = 0;
    int total_factors = 0;
    int i;
    int lf            = 0;

    /* initialize on cpu: compute rho, one, and Rsq */
    mpz_t rsq, zn;
    mpz_init(rsq);
    mpz_init(zn);
    for (i = 0; i < t->array_sz; i++) {
        t->rho[i] = multiplicative_neg_inverse32(t->modulus_in[i]);
        t->one[i] = ((uint64_t)0 - t->modulus_in[i]) % t->modulus_in[i];
        mpz_set_ui(rsq, 1);
        mpz_mul_2exp(rsq, rsq, 128);
        t->rsq[i] = mpz_tdiv_ui(rsq, t->modulus_in[i]);
    }

    OCL_TRY(clEnqueueWriteBuffer(t->queue, t->gpu_rsq_array, CL_FALSE,
        0, t->array_sz * sizeof(uint64_t), t->rsq, 0, NULL, NULL))
    OCL_TRY(clEnqueueWriteBuffer(t->queue, t->gpu_rho_array, CL_FALSE,
        0, t->array_sz * sizeof(uint32_t), t->rho, 0, NULL, NULL))
    OCL_TRY(clEnqueueWriteBuffer(t->queue, t->gpu_one_array, CL_FALSE,
        0, t->array_sz * sizeof(uint64_t), t->one, 0, NULL, NULL))

    launch = t->launch + GPU_ECM_VEC;

    int orig_size  = t->array_sz;
    int nofactors  = 0;

    while ((curve < t->curves_2lp) && (total_factors < orig_size)) {

        gpu_args[0].int32_arg  = t->array_sz;
        gpu_args[1].ptr_arg    = t->gpu_n_array;
        gpu_args[2].ptr_arg    = t->gpu_rho_array;
        gpu_args[3].ptr_arg    = t->gpu_one_array;
        gpu_args[4].ptr_arg    = t->gpu_rsq_array;
        gpu_args[5].ptr_arg    = t->gpu_u32_array;
        gpu_args[6].ptr_arg    = t->gpu_a_array;
        gpu_args[7].uint32_arg = 205;
        gpu_args[8].int32_arg  = curve;

        gpu_launch_set(launch, gpu_args);

        /* Launch kernel.
         * Replaces: cuFuncSetBlockShape + cuLaunchGridAsync.
         * global_work_size = total threads; local_work_size = block size. */
        size_t local_sz  = (size_t)threads_per_block;
        size_t global_sz = (size_t)(num_blocks * threads_per_block);
        OCL_TRY(clEnqueueNDRangeKernel(t->queue, launch->kernel_func,
            1, NULL, &global_sz, &local_sz, 0, NULL, NULL))

        /* copy factors back to host */
        OCL_TRY(clEnqueueReadBuffer(t->queue, t->gpu_a_array, CL_FALSE,
            0, t->array_sz * sizeof(uint64_t), t->a, 0, NULL, NULL))

        /* Must flush before touching host results */
        OCL_TRY(clFinish(t->queue))

        /* swap factored inputs to end of list */
        int n = t->array_sz;
        int c = 0;
        for (i = 0; i < n; i++) {
            mpz_set_ull(zn, t->modulus_in[i]);
            uint64_t factor = t->a[i];

            if ((factor > 1) && (factor < t->modulus_in[i])) {

                mpz_set_ull(rsq, factor);

                int bits1 = mpz_sizeinbase(rsq, 2);
                mpz_tdiv_q(rsq, zn, rsq);
                int bits2 = mpz_sizeinbase(rsq, 2);

                if (t->rb_idx_2lp[i] < t->rb->num_relations)
                {
                    if ((bits1 <= t->lpb_2lp) && (bits2 <= t->lpb_2lp))
                    {
                        // valid factorization, save it.
                        if (t->mode_2lp == 0)
                        {
                            cofactor_t* c = t->rb->relations + t->rb_idx_2lp[i];

                            if (t->first_side == 0)
                            {
                                c->lp_r[0] = factor;
                                c->lp_r[1] = t->modulus_in[i] / factor;
                            }
                            else
                            {
                                c->lp_a[0] = factor;
                                c->lp_a[1] = t->modulus_in[i] / factor;
                            }
                        }
                        else
                        {
                            // in mode 1 these are a-side LPs whose indices
                            // into the RB have been loaded into the r-side array.
                            cofactor_t* c = t->rb->relations + t->rb_idx_2lp[i];
                            uint8_t success = c->success;

                            // we end up here because the 3LP kernel already
                            // found one valid factor and put it in position 0.
                            // here we record the final two factors.
                            if (t->first_side == 0)
                            {
                                c->lp_a[1] = factor;
                                c->lp_a[2] = t->modulus_in[i] / factor;
                            }
                            else
                            {
                                c->lp_r[1] = factor;
                                c->lp_r[2] = t->modulus_in[i] / factor;
                            }
                            c->success |= 0x0f;
                        }
                        t->num_factors_2lp++;
                    }
                    else
                    {
                        if (t->mode_2lp == 0)
                        {
                            // gmp_printf("2LP invalid factor sizes: %Zx = %Zx (%d) * ",
                            //  	zn, rsq, mpz_sizeinbase(rsq, 2));
                            // mpz_tdiv_q(rsq, zn, rsq);
                            // gmp_printf("%Zx (%d)\n", rsq, mpz_sizeinbase(rsq, 2));
                            cofactor_t* c = t->rb->relations + t->rb_idx_2lp[i];
                            c->success = 0;
                        }
                    }
                }
                else
                {
                    printf("invalid relation index %d (of %d)\n", t->rb_idx_2lp[i],
                        t->rb->num_relations);
                }

                t->modulus_in[i]  = t->modulus_in[n - 1];
                t->rsq[i]         = t->rsq[n - 1];
                t->one[i]         = t->one[n - 1];
                t->rho[i]         = t->rho[n - 1];
                t->a[i]           = t->a[n - 1];
                t->rb_idx_2lp[i]  = t->rb_idx_2lp[n - 1];

                //// also sync the sigma
                //t->u32_array[i] = t->u32_array[n - 1];
                //
                //// advance the sigma every curve
                //uint64_t sigma64 = (uint64_t)t->u32_array[i];
                //t->u32_array[i] = uecm_lcg_rand_32B(7, 0xffffffff, &sigma64);

                // shrink the list
                n--;
                c++;

                // visit this index again
                i--;
            }
            //else
            //{
            //    // advance the sigma every curve
            //    uint64_t sigma64 = (uint64_t)t->u32_array[i];
            //    t->u32_array[i] = uecm_lcg_rand_32B(7, 0xffffffff, &sigma64);
            //}

        }

        int lastfactors = total_factors;
        total_factors  += c;
        t->array_sz     = n;

        if (lastfactors == total_factors)
            nofactors++;

        num_blocks = t->array_sz / threads_per_block
                   + ((t->array_sz % threads_per_block) > 0);

        if (t->array_sz == 0)
            break;

        // temporary: update and sync sigma locally for debug purposes
        //OCL_TRY(clEnqueueWriteBuffer(t->queue, t->gpu_u32_array, CL_FALSE,
        //    0, t->array_sz * sizeof(uint32_t), t->u32_array, 0, NULL, NULL))

        OCL_TRY(clEnqueueWriteBuffer(t->queue, t->gpu_n_array, CL_FALSE,
            0, t->array_sz * sizeof(uint64_t), t->modulus_in, 0, NULL, NULL))
        OCL_TRY(clEnqueueWriteBuffer(t->queue, t->gpu_rsq_array, CL_FALSE,
            0, t->array_sz * sizeof(uint64_t), t->rsq, 0, NULL, NULL))
        OCL_TRY(clEnqueueWriteBuffer(t->queue, t->gpu_rho_array, CL_FALSE,
            0, t->array_sz * sizeof(uint32_t), t->rho, 0, NULL, NULL))
        OCL_TRY(clEnqueueWriteBuffer(t->queue, t->gpu_one_array, CL_FALSE,
            0, t->array_sz * sizeof(uint64_t), t->one, 0, NULL, NULL))

        curve += 1;
    }

    // anything unfactored we mark as invalid
    nofactors = 0;
    if (t->mode_2lp == 0)
    {
        for (i = 0; i < t->array_sz; i++) {
            cofactor_t* c = t->rb->relations + t->rb_idx_2lp[i];
            c->success = 0;
            nofactors++;
        }

        printf("marked %d unfactored residues as invalid\n", nofactors);
    }

    /* end timing */
    OCL_TRY(clFinish(t->queue))
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    elapsed_ms = timespec_elapsed_ms(&ts_start, &ts_end);

    printf("found %d total factors (%d valid) in %1.4f ms\n",
           total_factors, t->num_factors_2lp, elapsed_ms);

    mpz_clear(rsq);
    mpz_clear(zn);

    /* final sync */
    OCL_TRY(clFinish(t->queue))

    return quit;
}

/* -----------------------------------------------------------------------
 * do_gpu_ecm_96b
 * --------------------------------------------------------------------- */
uint32_t
do_gpu_ecm_96b(device_thread_ctx_t *t)
{
    uint32_t     quit = 0;
    gpu_arg_t    gpu_args[GPU_MAX_KERNEL_ARGS];
    gpu_launch_t *launch;
    struct timespec ts_start, ts_end;
    float        elapsed_ms;

    int threads_per_block = 128;
    int num_blocks = t->array_sz / threads_per_block
                   + ((t->array_sz % threads_per_block) > 0);

    printf("commencing gpu 96-bit ecm on %d inputs (b1 = %d, b2 = %d, curves = %d)\n",
           t->array_sz, t->b1_3lp, t->b2_3lp * t->b1_3lp, t->curves_3lp);
    fflush(stdout);

    OCL_TRY(clEnqueueWriteBuffer(t->queue, t->gpu_u32_array, CL_FALSE,
        0, t->array_sz * sizeof(uint32_t), t->u32_array, 0, NULL, NULL))
    OCL_TRY(clEnqueueWriteBuffer(t->queue, t->gpu_n_array, CL_FALSE,
        0, t->array_sz * 3 * sizeof(uint32_t), t->modulus96_in, 0, NULL, NULL))

    clock_gettime(CLOCK_MONOTONIC, &ts_start);

    int total_factors = 0;
    int i;

    mpz_t rsq, zn, zf, zc;
    mpz_init(rsq); mpz_init(zn); mpz_init(zf); mpz_init(zc);
    mpz_t fac[3];
    mpz_init(fac[0]); mpz_init(fac[1]); mpz_init(fac[2]);
    mpz_t *flist = fac;

    for (i = 0; i < t->array_sz; i++) {
        bignum32_to_mpz(zn, &t->modulus96_in[3 * i]);
        uint32_t n32 = t->modulus96_in[3 * i];
        t->rho[i] = multiplicative_neg_inverse32(n32);

        mpz_set_ui(rsq, 1);
        mpz_mul_2exp(rsq, rsq, 96);
        mpz_sub(rsq, rsq, zn);
        mpz_tdiv_r(rsq, rsq, zn);
        mpz_to_bignum32(&t->one96[3 * i], rsq);

        mpz_set_ui(rsq, 1);
        mpz_mul_2exp(rsq, rsq, 192);
        mpz_tdiv_r(rsq, rsq, zn);
        mpz_to_bignum32(&t->rsq96[3 * i], rsq);
    }

    OCL_TRY(clEnqueueWriteBuffer(t->queue, t->gpu_rsq_array, CL_FALSE,
        0, t->array_sz * 3 * sizeof(uint32_t), t->rsq96, 0, NULL, NULL))
    OCL_TRY(clEnqueueWriteBuffer(t->queue, t->gpu_rho_array, CL_FALSE,
        0, t->array_sz * sizeof(uint32_t), t->rho, 0, NULL, NULL))
    OCL_TRY(clEnqueueWriteBuffer(t->queue, t->gpu_one_array, CL_FALSE,
        0, t->array_sz * 3 * sizeof(uint32_t), t->one96, 0, NULL, NULL))

    launch = t->launch + GPU_ECM96_VEC;

    uint64_t lcg        = 0xbaddecafbaddecafull;
    int total_curves    = 0;
    int num2lp_retest   = 0;
    int orig_size       = t->array_sz;
    int no_factors      = 0;
    int max_no_factors  = t->stop_nofactor;
    int max_curves      = t->curves_3lp;
    int curve           = 0;
    int num_mpqs        = 0;
    int mpqs_success    = 0;
    int mpqs_err        = 0;

    while ((curve < max_curves) && (total_factors < orig_size)) {

        gpu_args[0].int32_arg  = t->array_sz;
        gpu_args[1].ptr_arg    = t->gpu_n_array;
        gpu_args[2].ptr_arg    = t->gpu_rho_array;
        gpu_args[3].ptr_arg    = t->gpu_one_array;
        gpu_args[4].ptr_arg    = t->gpu_rsq_array;
        gpu_args[5].ptr_arg    = t->gpu_u32_array;
        gpu_args[6].ptr_arg    = t->gpu_res32_array;
        gpu_args[7].uint32_arg = t->b1_3lp;
        gpu_args[8].uint32_arg = t->b2_3lp * t->b1_3lp;
        gpu_args[9].int32_arg  = curve;

        gpu_launch_set(launch, gpu_args);

        int last_factors = num2lp_retest;

        //printf("launching curve %d on %u 96-bit inputs\n", curve, t->array_sz);

        size_t local_sz  = (size_t)threads_per_block;
        size_t global_sz = (size_t)(num_blocks * threads_per_block);

        //printf("kernel %s, inputs %d, size <%d,%d>, ", gpu_kernel_names[GPU_ECM96_VEC],
        //    t->array_sz, num_blocks, threads_per_block);

        OCL_TRY(clEnqueueNDRangeKernel(t->queue, launch->kernel_func,
            1, NULL, &global_sz, &local_sz, 0, NULL, NULL))

        OCL_TRY(clEnqueueReadBuffer(t->queue, t->gpu_res32_array, CL_FALSE,
            0, t->array_sz * 3 * sizeof(uint32_t), t->factors96, 0, NULL, NULL))

        OCL_TRY(clFinish(t->queue))

        total_curves += threads_per_block * num_blocks;

        int n = t->array_sz;
        int c = 0;
        for (i = 0; i < n; i++) {
            bignum32_to_mpz(zf, &t->factors96[3 * i]);
            bignum32_to_mpz(zn, &t->modulus96_in[3 * i]);

            if ((mpz_cmp_ui(zf, 1) > 0) && (mpz_cmp(zf, zn) < 0)) {
                mpz_tdiv_r(rsq, zn, zf);
                if (mpz_cmp_ui(rsq, 0) == 0) {

                    num2lp_retest = handle_96b_factorization(t, i, zf, zc, zn,
                        flist, num2lp_retest, &mpqs_success, &mpqs_err, &num_mpqs);

                    t->modulus96_in[3 * i + 0] = t->modulus96_in[3 * (n - 1) + 0];
                    t->modulus96_in[3 * i + 1] = t->modulus96_in[3 * (n - 1) + 1];
                    t->modulus96_in[3 * i + 2] = t->modulus96_in[3 * (n - 1) + 2];
                    t->rsq96[3 * i + 0] = t->rsq96[3 * (n - 1) + 0];
                    t->rsq96[3 * i + 1] = t->rsq96[3 * (n - 1) + 1];
                    t->rsq96[3 * i + 2] = t->rsq96[3 * (n - 1) + 2];
                    t->one96[3 * i + 0] = t->one96[3 * (n - 1) + 0];
                    t->one96[3 * i + 1] = t->one96[3 * (n - 1) + 1];
                    t->one96[3 * i + 2] = t->one96[3 * (n - 1) + 2];
                    t->rho[i] = t->rho[n - 1];
                    t->factors96[3 * i + 0] = t->factors96[3 * (n - 1) + 0];
                    t->factors96[3 * i + 1] = t->factors96[3 * (n - 1) + 1];
                    t->factors96[3 * i + 2] = t->factors96[3 * (n - 1) + 2];
                    t->rb_idx_3lp[i] = t->rb_idx_3lp[n - 1];

                    // shrink the list
                    n--;
                    c++;

                    // visit this index again
                    i--;
                }
            }
        }

        //printf("found %d factors, %d marked for 2lp retest\n", 
        //    c, num2lp_retest - last_factors);

        if (last_factors == num2lp_retest)
            no_factors++;

        total_factors  += c;
        t->array_sz     = n;

        if (no_factors >= max_no_factors) {
            printf("halting after %d curves (%d/%d/%d mpqs success/calls/err): "
                   "%d curves yielded no factors\n",
                   curve + 1, mpqs_success, num_mpqs, mpqs_err, max_no_factors);
            printf("giving up on %d likely unproductive 96-bit residues\n", n);
            curve = max_curves;
        } else if (curve == (max_curves - 1)) {
            curve = max_curves;
            printf("halting after running the maximum specified %d curves "
                   "(%d/%d mpqs calls)\n", max_curves, mpqs_success, num_mpqs);
        }

        num_blocks = t->array_sz / threads_per_block
                   + ((t->array_sz % threads_per_block) > 0);

        if (t->array_sz == 0)
            break;

        OCL_TRY(clEnqueueWriteBuffer(t->queue, t->gpu_n_array, CL_FALSE,
            0, t->array_sz * 3 * sizeof(uint32_t), t->modulus96_in, 0, NULL, NULL))
        OCL_TRY(clEnqueueWriteBuffer(t->queue, t->gpu_rsq_array, CL_FALSE,
            0, t->array_sz * 3 * sizeof(uint32_t), t->rsq96, 0, NULL, NULL))
        OCL_TRY(clEnqueueWriteBuffer(t->queue, t->gpu_rho_array, CL_FALSE,
            0, t->array_sz * sizeof(uint32_t), t->rho, 0, NULL, NULL))
        OCL_TRY(clEnqueueWriteBuffer(t->queue, t->gpu_one_array, CL_FALSE,
            0, t->array_sz * 3 * sizeof(uint32_t), t->one96, 0, NULL, NULL))

        curve += 1;
    }

    mpz_clear(rsq); mpz_clear(zn); mpz_clear(zf); mpz_clear(zc);
    mpz_clear(fac[0]); mpz_clear(fac[1]); mpz_clear(fac[2]);

    OCL_TRY(clFinish(t->queue))
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    elapsed_ms = timespec_elapsed_ms(&ts_start, &ts_end);

    printf("found %d total factors with %d total curves in %1.4f ms\n",
           total_factors, total_curves, elapsed_ms);
    t->array_sz = total_factors;

    OCL_TRY(clFinish(t->queue))

    if (num2lp_retest > 0)
    {
        printf("running 2LP kernel on %d 3LP-cofactors\n", num2lp_retest);
        t->array_sz = num2lp_retest;
        t->mode_2lp = 1;
        t->num_factors_2lp = 0;
        do_gpu_ecm64(t);

        printf("found %d valid factors\n", t->num_factors_2lp);
        t->num_factors_3lp = t->num_factors_2lp;
    }
    else
    {
        t->num_factors_3lp = 0;
    }

    return quit;
}

/* -----------------------------------------------------------------------
 * do_gpu_pm1_96b
 * --------------------------------------------------------------------- */
uint32_t
do_gpu_pm1_96b(device_thread_ctx_t *t)
{
    uint32_t     quit = 0;
    gpu_arg_t    gpu_args[GPU_MAX_KERNEL_ARGS];
    gpu_launch_t *launch;
    struct timespec ts_start, ts_end;
    float        elapsed_ms;

    int threads_per_block = 128;
    int num_blocks = t->array_sz / threads_per_block
                   + ((t->array_sz % threads_per_block) > 0);

    printf("commencing gpu 96-bit pm1 on %d inputs (b1 = %d, b2 = %d)\n",
           t->array_sz, 500, 500 * 50);
    fflush(stdout);

    OCL_TRY(clEnqueueWriteBuffer(t->queue, t->gpu_n_array, CL_FALSE,
        0, t->array_sz * 3 * sizeof(uint32_t), t->modulus96_in, 0, NULL, NULL))

    clock_gettime(CLOCK_MONOTONIC, &ts_start);

    int total_factors = 0;
    int i;

    mpz_t rsq, zn, zf, zc;
    mpz_init(rsq); mpz_init(zn); mpz_init(zf); mpz_init(zc);
    mpz_t fac[3];
    mpz_init(fac[0]); mpz_init(fac[1]); mpz_init(fac[2]);
    mpz_t *flist = fac;

    for (i = 0; i < t->array_sz; i++) {
        bignum32_to_mpz(zn, &t->modulus96_in[3 * i]);
        uint32_t n32 = t->modulus96_in[3 * i];
        t->rho[i] = multiplicative_neg_inverse32(n32);

        mpz_set_ui(rsq, 1);
        mpz_mul_2exp(rsq, rsq, 96);
        mpz_sub(rsq, rsq, zn);
        mpz_tdiv_r(rsq, rsq, zn);
        mpz_to_bignum32(&t->one96[3 * i], rsq);
    }

    OCL_TRY(clEnqueueWriteBuffer(t->queue, t->gpu_rho_array, CL_FALSE,
        0, t->array_sz * sizeof(uint32_t), t->rho, 0, NULL, NULL))
    OCL_TRY(clEnqueueWriteBuffer(t->queue, t->gpu_one_array, CL_FALSE,
        0, t->array_sz * 3 * sizeof(uint32_t), t->one96, 0, NULL, NULL))

    launch = t->launch + GPU_PM196_VEC;
    int orig_size     = t->array_sz;
    int num2lp_retest = 0;
    int num_mpqs      = 0;
    int mpqs_success  = 0;
    int mpqs_err      = 0;

    /* P-1 runs once */
    {
        gpu_args[0].int32_arg  = t->array_sz;
        gpu_args[1].ptr_arg    = t->gpu_n_array;
        gpu_args[2].ptr_arg    = t->gpu_rho_array;
        gpu_args[3].ptr_arg    = t->gpu_one_array;
        gpu_args[4].ptr_arg    = t->gpu_res32_array;
        gpu_args[5].uint32_arg = 500;
        gpu_args[6].uint32_arg = t->b2_3lp * t->b1_3lp;

        gpu_launch_set(launch, gpu_args);

        size_t local_sz  = (size_t)threads_per_block;
        size_t global_sz = (size_t)(num_blocks * threads_per_block);
        OCL_TRY(clEnqueueNDRangeKernel(t->queue, launch->kernel_func,
            1, NULL, &global_sz, &local_sz, 0, NULL, NULL))

        OCL_TRY(clEnqueueReadBuffer(t->queue, t->gpu_res32_array, CL_FALSE,
            0, t->array_sz * 3 * sizeof(uint32_t), t->factors96, 0, NULL, NULL))

        OCL_TRY(clFinish(t->queue))

        int n = t->array_sz;
        int c = 0;
        for (i = 0; i < n; i++) {
            bignum32_to_mpz(zf, &t->factors96[3 * i]);
            bignum32_to_mpz(zn, &t->modulus96_in[3 * i]);

            if ((mpz_cmp_ui(zf, 1) > 0) && (mpz_cmp(zf, zn) < 0))
            {
                mpz_tdiv_r(rsq, zn, zf);

                if (mpz_cmp_ui(rsq, 0) == 0)
                {
                    num2lp_retest = handle_96b_factorization(t, i, zf, zc, zn, flist,
                        num2lp_retest, &mpqs_success, &mpqs_err, &num_mpqs);

                    t->modulus96_in[3*i+0] = t->modulus96_in[3*(n-1)+0];
                    t->modulus96_in[3*i+1] = t->modulus96_in[3*(n-1)+1];
                    t->modulus96_in[3*i+2] = t->modulus96_in[3*(n-1)+2];
                    t->rsq96[3*i+0]        = t->rsq96[3*(n-1)+0];
                    t->rsq96[3*i+1]        = t->rsq96[3*(n-1)+1];
                    t->rsq96[3*i+2]        = t->rsq96[3*(n-1)+2];
                    t->one96[3*i+0]        = t->one96[3*(n-1)+0];
                    t->one96[3*i+1]        = t->one96[3*(n-1)+1];
                    t->one96[3*i+2]        = t->one96[3*(n-1)+2];
                    t->rho[i]              = t->rho[n-1];
                    t->factors96[3*i+0]    = t->factors96[3*(n-1)+0];
                    t->factors96[3*i+1]    = t->factors96[3*(n-1)+1];
                    t->factors96[3*i+2]    = t->factors96[3*(n-1)+2];
                    t->rb_idx_3lp[i]       = t->rb_idx_3lp[n-1];
                    n--;
                    c++;
                    i--;
                }
                else
                {
                    gmp_printf("found bogus factor %Zd, does not divide %Zd\n",
                        zf, zn);
                }
            }
        }
        total_factors = c;
        t->array_sz   = n;
    }

    mpz_clear(rsq); mpz_clear(zn); mpz_clear(zf); mpz_clear(zc);
    mpz_clear(fac[0]); mpz_clear(fac[1]); mpz_clear(fac[2]);

    OCL_TRY(clFinish(t->queue))
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    elapsed_ms = timespec_elapsed_ms(&ts_start, &ts_end);

    orig_size = t->array_sz;
    printf("found %d total factors (%d/%d/%d mpqs success/calls/err) with P-1 in %1.4f ms\n",
           total_factors, mpqs_success, num_mpqs, mpqs_err, elapsed_ms);
    t->array_sz = total_factors;

    OCL_TRY(clFinish(t->queue))

    if (total_factors > 0)
    {
        printf("running 2LP kernel on %d 3LP-cofactors\n", num2lp_retest);
        t->array_sz = num2lp_retest;
        t->mode_2lp = 1;
        t->num_factors_2lp = 0;
        do_gpu_ecm64(t);

        printf("found %d valid factors\n", t->num_factors_2lp);

        t->num_factors_3lp = t->num_factors_2lp;
        t->array_sz = orig_size;
    }
    else
    {
        t->num_factors_3lp = 0;
        t->array_sz = orig_size;
    }

    return quit;
}

uint32_t gpu_cofactorization(device_thread_ctx_t* t)
{
    uint32_t quit = 0;
    int i;
    int j;

    // which list is 2LP?
    // load the residues into the 64-bit input array.
    if (t->first_side == 0)
    {
        // 2LP side has 2 words per input
        printf("setting up gpu to factor %d r-side 2LPs\n", t->numres_r);
        for (i = 0; i < t->numres_r; i++) {
            t->modulus_in[i] = ((uint64_t)t->residues_r_in[i * 2 + 1] << 32) |
                (uint64_t)t->residues_r_in[i * 2 + 0];
        }

        // the 2LP factorization code is agnostic to side, so
        // point it to which side it should be tracking.
        t->array_sz = t->numres_r;

        memcpy(t->rb_idx_2lp, t->rb_idx_r, t->numres_r * sizeof(uint32_t));
        //t->rb_idx_2lp = t->rb_idx_r;		
    }
    else
    {
        // 2LP side has 2 words per input
        printf("setting up gpu to factor %d a-side 2LPs\n", t->numres_a);
        for (i = 0; i < t->numres_a; i++) {
            t->modulus_in[i] = ((uint64_t)t->residues_a_in[i * 2 + 1] << 32) |
                (uint64_t)t->residues_a_in[i * 2 + 0];
        }
        // the 2LP factorization code is agnostic to side, so
        // point it to which side it should be tracking.
        t->array_sz = t->numres_a;

        memcpy(t->rb_idx_2lp, t->rb_idx_a, t->numres_a * sizeof(uint32_t));
        //t->rb_idx_2lp = t->rb_idx_a;
    }

    // try to completely factor the 2LP list.  The last handful
    // of curves don't typically make sense to run on the gpu (only
    // a few inputs left) but it also doesn't take much time, so
    // to be lazy we just finish it all here.
    t->mode_2lp = 0;
    t->num_factors_2lp = 0;
    do_gpu_ecm64(t);

    // sometimes the factors of a 2LP are not correctly sized.
    // when that happens, we can ignore the corresponding 3LP side cofactor.
    // here we build up a list of 3lp candidates to try to factor
    // with 96-bit ecm code.
    j = 0;

    for (i = 0; i < t->rb->num_relations; i++)
    {
        if (t->rb->relations[i].success == 0)
        {
            // skip relations in the rb that didn't have a 
            // valid 2LP factorization.
        }
        else
        {
            // the 3lp factorization code needs to know the
            // moduli to factor and the indices of those moduli in
            // the rb structure.  Copy from whichever side has
            // the 3lps for this successful 2lp-side factorization.
            if (t->first_side == 0)
            {
                t->modulus96_in[3 * j + 0] = t->residues_a_in[3 * i + 0];
                t->modulus96_in[3 * j + 1] = t->residues_a_in[3 * i + 1];
                t->modulus96_in[3 * j + 2] = t->residues_a_in[3 * i + 2];
                t->rb_idx_3lp[j] = t->rb_idx_a[i];
            }
            else
            {
                t->modulus96_in[3 * j + 0] = t->residues_r_in[3 * i + 0];
                t->modulus96_in[3 * j + 1] = t->residues_r_in[3 * i + 1];
                t->modulus96_in[3 * j + 2] = t->residues_r_in[3 * i + 2];
                t->rb_idx_3lp[j] = t->rb_idx_r[i];
            }

            // 2LP factorization was good
            j++;
        }
    }

    t->array_sz = j;
    printf("ignoring %d 3LP-side cofactors due to invalid 2LP-side factorizations\n",
        t->rb->num_relations - j);

    // now run the 3LP kernels
    t->num_factors_3lp = 0;

    do_gpu_pm1_96b(t);
    do_gpu_ecm_96b(t);

    // any survivors have now survived both sides (had factors 
    // found on both R and A sides). 
    // double check the number that have success fully flagged
    t->num_factors_3lp = 0;
    t->rb->num_success = 0;
    for (i = 0; i < t->rb->num_relations; i++)
    {
        if (t->rb->relations[i].success == 0xff)
        {
            t->num_factors_3lp++;
            t->rb->relations[i].success = 1;
            t->rb->num_success++;
        }
        else
        {
            t->rb->relations[i].success = 0;
        }
    }
    printf("%d relations have been flagged as completely factored\n", t->rb->num_success);

    return quit;
}

/* -----------------------------------------------------------------------
 * gpu_device_init
 *
 * Replaces the CUDA version that called gpu_init() + cuCtxCreate etc.
 * We only fill the device_ctx_t here; the OpenCL context is created later
 * in gpu_ctx_init where we have both device and context together.
 * --------------------------------------------------------------------- */
device_ctx_t *
gpu_device_init(int which_gpu, int verbose)
{
    gpu_config_t  gpu_config;
    gpu_info_t   *gpu_info;

    device_ctx_t *d = (device_ctx_t *)xcalloc(1, sizeof(device_ctx_t));

    gpu_init(&gpu_config);
    if (gpu_config.num_gpu == 0) {
        printf("error: no OpenCL-capable GPUs found\n");
        exit(-1);
    }

    d->gpunum   = which_gpu;
    d->gpu_info = gpu_info = (gpu_info_t *)xmalloc(sizeof(gpu_info_t));
    memcpy(gpu_info, gpu_config.info + which_gpu, sizeof(gpu_info_t));

    printf("using GPU %u (%s)\n", which_gpu, gpu_info->name);
    printf("selected card has OpenCL version %d.%d\n",
           gpu_info->compute_version_major,
           gpu_info->compute_version_minor);
    if (verbose)
    {
        printf("more GPU info:\n");
        printf("\tmax_grid_size: %d x %d x %d\n",
            gpu_info->max_grid_size[0],
            gpu_info->max_grid_size[1],
            gpu_info->max_grid_size[2]);
        printf("\tglobal_mem_size: %zd\n", gpu_info->global_mem_size);
        printf("\tconstant_mem_size: %d\n", gpu_info->constant_mem_size);
        printf("\tmax_threads_per_block: %d\n", gpu_info->max_threads_per_block);
        printf("\tmax_thread_dim: %d x %d x %d\n",
            gpu_info->max_thread_dim[0],
            gpu_info->max_thread_dim[1],
            gpu_info->max_thread_dim[2]);
        printf("\tnum_compute_units: %d\n", gpu_info->num_compute_units);
        printf("\tregisters_per_block: %d\n", gpu_info->registers_per_block);
        printf("\tshared_mem_size: %d\n", gpu_info->shared_mem_size);
        printf("\twarp_size: %d\n", gpu_info->warp_size);
    }
    return d;
}

void
gpu_dev_free(device_ctx_t *d)
{
    free(d->gpu_info);
    free(d);
}

/* -----------------------------------------------------------------------
 * gpu_ctx_init
 *
 * CUDA version:
 *   cuCtxCreate     -- create a per-thread CUDA context
 *   cuModuleLoad    -- load a pre-compiled .ptx file
 *   gpu_launch_init -- look up CUfunction handles
 *   cuStreamCreate  -- create an async stream
 *   cuEventCreate x2
 *
 * OpenCL version:
 *   clCreateContext         -- one context per device
 *   clCreateCommandQueue    -- in-order queue (equivalent to CUDA stream
 *                              with CU_CTX_BLOCKING_SYNC)
 *   clCreateProgramWithSource + clBuildProgram
 *                           -- replaces cuModuleLoad from .ptx;
 *                              loads the .cl source file at runtime.
 *   gpu_launch_init x N     -- creates cl_kernel objects
 *
 * The .ptx file selection logic (sm_20 / sm_30 / sm_35 / sm_50 / sm_80
 * / sm_90) is replaced with a single .cl file -- "opencl_tinyecm.cl" --
 * that the AMD driver JIT-compiles for whatever GPU is present.
 *
 * Events (CUevent start_event / end_event) are replaced by wall-clock
 * timing in each do_gpu_* function, so we don't create them here.
 * --------------------------------------------------------------------- */

/* Helper: read an entire text file into a malloc'd buffer. */
static char *
read_file(const char *path)
{
    FILE  *f = fopen(path, "rb");
    char  *buf;
    long   sz;

    if (!f) {
        fprintf(stderr, "gpu_ctx_init: cannot open kernel file '%s'\n", path);
        exit(-1);
    }
    fseek(f, 0, SEEK_END);
    sz = ftell(f);
    rewind(f);
    buf = (char *)xmalloc((size_t)sz + 1);
    if (fread(buf, 1, (size_t)sz, f) != (size_t)sz) {
        fprintf(stderr, "gpu_ctx_init: read error on '%s'\n", path);
        exit(-1);
    }
    buf[sz] = '\0';
    fclose(f);
    return buf;
}

device_thread_ctx_t*
gpu_ctx_init(device_ctx_t* d)
{
    device_thread_ctx_t* t;
    cl_int               err;
    cl_device_id         dev = d->gpu_info->device_handle;
    cl_platform_id       plat = d->gpu_info->platform_handle;
    int                  i;

    t = (device_thread_ctx_t*)xcalloc(1, sizeof(device_thread_ctx_t));
    t->dev = d;

    /* Create OpenCL context.
     * Replaces: cuCtxCreate(&t->gpu_context, ..., device_handle)        */
    cl_context_properties props[] = {
        CL_CONTEXT_PLATFORM, (cl_context_properties)plat, 0
    };
    t->gpu_context = clCreateContext(props, 1, &dev, NULL, NULL, &err);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "gpu_ctx_init: clCreateContext failed: %s\n",
            clGetErrorString(err));
        exit(-1);
    }

    /* Load and build the OpenCL kernel source.
     * Replaces: cuModuleLoad(&t->gpu_module, ptxfile)
     *
     * AMD's comgr compiler runs a full LLVM pipeline on first build, which
     * can take several minutes for dense kernel code like this.  We cache
     * the compiled binary to disk (keyed on device name + build options)
     * so subsequent runs load in under a second.
     *
     * Cache logic:
     *   1. If <cache_file> exists on disk, load it and call
     *      clCreateProgramWithBinary -- no compilation needed.
     *   2. Otherwise compile from source with clBuildProgram, then
     *      extract the binary with clGetProgramInfo(CL_PROGRAM_BINARIES)
     *      and write it to <cache_file> for next time.
     *
     * The cache filename embeds the device name so it is automatically
     * invalidated when run on a different GPU.  Delete the file manually
     * if you change the .cl sources or build options.                    */


    const char* intrinsics_file = "factor/opencl_intrinsics.cl";
    const char* tinyecm_file = "factor/opencl_tinyecm.cl";
    const char* build_options = "-cl-std=CL2.0 -cl-mad-enable";

    /* Build a cache filename: "ocl_ecm_<devname>.bin"
     * Sanitise the device name -- replace spaces and slashes with '_'.  */
    char cache_file[256];
    {
        char safe_name[128];
        const char* src_n = d->gpu_info->name;
        char* dst_n = safe_name;
        size_t      lim = sizeof(safe_name) - 1;
        while (*src_n && (size_t)(dst_n - safe_name) < lim) {
            char c = *src_n++;
            *dst_n++ = (c == ' ' || c == '/' || c == '\\') ? '_' : c;
        }
        *dst_n = '\0';
        snprintf(cache_file, sizeof(cache_file),
            "ocl_ecm_%s.bin", safe_name);
    }

    /* --- Try loading from cache first -------------------------------- */
    int loaded_from_cache = 0;
    {
        FILE* f = fopen(cache_file, "rb");
        if (f) {
            fseek(f, 0, SEEK_END);
            long bin_sz = ftell(f);
            rewind(f);

            if (bin_sz > 0) {
                unsigned char* bin = (unsigned char*)xmalloc((size_t)bin_sz);
                if (fread(bin, 1, (size_t)bin_sz, f) == (size_t)bin_sz) {
                    const unsigned char* bin_ptr = bin;
                    size_t               bin_len = (size_t)bin_sz;
                    cl_int               bin_status = CL_SUCCESS;

                    t->gpu_program = clCreateProgramWithBinary(
                        t->gpu_context, 1, &dev,
                        &bin_len, &bin_ptr,
                        &bin_status, &err);

                    if (err == CL_SUCCESS && bin_status == CL_SUCCESS) {
                        err = clBuildProgram(t->gpu_program, 1, &dev,
                            build_options, NULL, NULL);
                        if (err == CL_SUCCESS) {
                            printf("loaded compiled kernels from cache: %s\n",
                                cache_file);
                            loaded_from_cache = 1;
                        }
                        else {
                            /* Stale / corrupt cache -- fall through to recompile */
                            printf("cache load failed (stale?), recompiling...\n");
                            clReleaseProgram(t->gpu_program);
                            t->gpu_program = NULL;
                        }
                    }
                    else {
                        clReleaseProgram(t->gpu_program);
                        t->gpu_program = NULL;
                    }
                }
                free(bin);
            }
            fclose(f);
        }
    }

    /* --- Compile from source if cache miss or stale cache ------------ */
    if (!loaded_from_cache) {
        printf("compiling kernels from %s + %s (this takes a while the first time)...\n",
            intrinsics_file, tinyecm_file);
        fflush(stdout);

        char* src_intrinsics = read_file(intrinsics_file);
        char* src_tinyecm = read_file(tinyecm_file);

        /* Skip the '#include "opencl_intrinsics.cl"' line in tinyecm so
         * the intrinsics are not compiled twice.                         */
        const char* tinyecm_body = src_tinyecm;
        {
            const char* inc = strstr(src_tinyecm,
                "#include \"opencl_intrinsics.cl\"");
            if (inc) {
                const char* nl = strchr(inc, '\n');
                if (nl) tinyecm_body = nl + 1;
            }
        }

        const char* sources[2] = { src_intrinsics, tinyecm_body };
        size_t      lengths[2] = { strlen(src_intrinsics),
                                    strlen(tinyecm_body) };

        t->gpu_program = clCreateProgramWithSource(t->gpu_context,
            2, sources, lengths, &err);
        free(src_intrinsics);
        free(src_tinyecm);
        if (err != CL_SUCCESS) {
            fprintf(stderr, "clCreateProgramWithSource failed: %s\n",
                clGetErrorString(err));
            exit(-1);
        }

        err = clBuildProgram(t->gpu_program, 1, &dev,
            build_options, NULL, NULL);
        if (err != CL_SUCCESS) {
            size_t log_sz = 0;
            clGetProgramBuildInfo(t->gpu_program, dev,
                CL_PROGRAM_BUILD_LOG, 0, NULL, &log_sz);
            char* log = (char*)xmalloc(log_sz + 1);
            clGetProgramBuildInfo(t->gpu_program, dev,
                CL_PROGRAM_BUILD_LOG, log_sz, log, NULL);
            log[log_sz] = '\0';
            fprintf(stderr, "clBuildProgram failed:\n%s\n", log);
            free(log);
            exit(-1);
        }

        printf("compilation successful, saving to cache: %s\n", cache_file);

        /* Extract the compiled binary and write it to the cache file.
         * clGetProgramInfo(CL_PROGRAM_BINARY_SIZES) returns one size per
         * device; we only have one device so we read index [0].          */
        size_t bin_sz = 0;
        clGetProgramInfo(t->gpu_program, CL_PROGRAM_BINARY_SIZES,
            sizeof(bin_sz), &bin_sz, NULL);
        if (bin_sz > 0) {
            unsigned char* bin = (unsigned char*)xmalloc(bin_sz);
            unsigned char* bins[1] = { bin };
            clGetProgramInfo(t->gpu_program, CL_PROGRAM_BINARIES,
                sizeof(bins), bins, NULL);

            FILE* f = fopen(cache_file, "wb");
            if (f) {
                fwrite(bin, 1, bin_sz, f);
                fclose(f);
            }
            else {
                fprintf(stderr, "warning: could not write cache file %s\n",
                    cache_file);
            }
            free(bin);
        }
    }

    printf("kernels ready\n");

    /* Create kernel objects.
     * Replaces: gpu_launch_init(t->gpu_module, name, args, launch)       */
    t->launch = (gpu_launch_t*)xmalloc(NUM_GPU_FUNCTIONS * sizeof(gpu_launch_t));

    printf("initializing kernels\n");
    for (i = 0; i < NUM_GPU_FUNCTIONS; i++) {
        gpu_launch_init(t->gpu_program, gpu_kernel_names[i],
            gpu_kernel_args + i,
            t->launch + i,
            dev);
    }

    /* Create in-order command queue.
     * Replaces: cuStreamCreate(&t->stream, 0)
     * CU_CTX_BLOCKING_SYNC -> 0 flags (in-order queue blocks on finish) */
    printf("creating command queue\n");
    t->queue = clCreateCommandQueue(t->gpu_context, dev, 0, &err);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "clCreateCommandQueue failed: %s\n",
            clGetErrorString(err));
        exit(-1);
    }

    /* No separate event objects needed; timing is done with clock_gettime
     * in each do_gpu_* function.
     * Replaces: cuEventCreate x2                                          */

    return t;
}


/* -----------------------------------------------------------------------
 * gpu_ctx_free
 *
 * Replaces: cuEventDestroy x2, cuStreamDestroy, cuCtxDestroy
 * --------------------------------------------------------------------- */
void
gpu_ctx_free(device_thread_ctx_t *d)
{
    int i;

    /* Release kernel objects */
    for (i = 0; i < NUM_GPU_FUNCTIONS; i++)
        clReleaseKernel(d->launch[i].kernel_func);
    free(d->launch);

    clReleaseCommandQueue(d->queue);    /* was: cuStreamDestroy    */
    clReleaseProgram(d->gpu_program);   /* was: part of cuCtxDestroy */
    clReleaseContext(d->gpu_context);   /* was: cuCtxDestroy       */
}

/* -----------------------------------------------------------------------
 * do_gpu_cofactorization  (external entry point)
 *
 * cuMemAlloc  -> clCreateBuffer(ctx, CL_MEM_READ_WRITE, size, NULL, &err)
 * cuMemFree   -> clReleaseMemObject
 *
 * Note on buffer sizes:
 *   gpu_n_array   -- 64-bit ECM uses uint64_t * array_sz
 *                    96-bit ECM uses uint32_t * 3 * array_sz
 *   The original allocated "sizeof(uint32_t) * 3 * array_sz" for gpu_n_array,
 *   gpu_one_array, and gpu_rsq_array (sufficient for both 64-bit and 96-bit
 *   usage since 3*uint32 >= uint64).  We keep the same conservative sizing.
 * --------------------------------------------------------------------- */
int
do_gpu_cofactorization(device_thread_ctx_t *t, relation_batch_t* rb, uint64_t *lcg,
                       int b1_3lp_ovr, int b2_3lp_ovr,
                       int b1_2lp_ovr, int b2_2lp_ovr,
                       int curves_3lp_ovr, int curves_2lp_ovr)
{
    cl_int err;
    int i;
    t->rb = rb;

    // the relation batch tracks all factors and metadata of a relation.
    // the gpu only cares about the large factors.  we need
    // to extract this info and put it in the data structures the gpu
    // code wants.
    t->array_sz = rb->num_relations;
    t->residues_r_in = (uint32_t*)xmalloc(sizeof(uint32_t) * MAX_RESIDUE_WORDS * rb->num_relations);
    t->residues_a_in = (uint32_t*)xmalloc(sizeof(uint32_t) * MAX_RESIDUE_WORDS * rb->num_relations);

    // the rb stores all factors, large and small, in one giant list.
    uint32_t* factors = rb->factors;

    // reference lists to relation_batch indices
    t->rb_idx_r = (uint32_t*)xmalloc(sizeof(uint32_t) * t->array_sz);
    t->rb_idx_a = (uint32_t*)xmalloc(sizeof(uint32_t) * t->array_sz);
    t->rb_idx_2lp = (uint32_t*)xmalloc(sizeof(uint32_t) * t->array_sz);
    t->rb_idx_3lp = (uint32_t*)xmalloc(sizeof(uint32_t) * t->array_sz);

    // get ready for siqs
    t->params = init_tinysiqs();

    uint32_t kr = 0;
    uint32_t ka = 0;
    t->numres_r = 0;
    t->numres_a = 0;
    t->first_side = -1;
    int max_words[2] = { 0,0 };
    int j;
    for (i = 0; i < rb->num_relations; i++)
    {
        // metadata for this relation
        cofactor_t* c = rb->relations + i;

        // advance past the small factors for this relation
        factors += c->num_factors_r;
        factors += c->num_factors_a;

        // initialize success as a successful 2LP factorization.
        // this only won't be true if we have an actual 2LP that
        // fails to factor into correctly sized factors, in which
        // case the 2LP ecm factorization code will zero it.  Many
        // 2LPs are already completely factored and the relation
        // just needs work on the 3LP side.
        c->success = 0xf0;

        // next in the list is the r-side residue
        if (c->lp_r_num_words > 0)
        {
            if (c->lp_r_num_words == 3)
            {
                if (t->first_side == 0)
                {
                    // we've already picked the r-side as first side,
                    // meaning that the a-side has 3lps as well.
                    printf("cuda cofactorization can't handle 3lp on both sides yet\n");
                    exit(1);
                }

                t->lpb_3lp = t->lpbr;		// 3lp's are on the r-side
                t->lpb_2lp = t->lpba;		// 2lp's are on the a-side
                t->first_side = 1;		// so do the a-side first.
                max_words[0] = MAX(c->lp_r_num_words, max_words[0]);
            }
            for (j = 0; j < c->lp_r_num_words; j++)
            {
                t->residues_r_in[kr++] = factors[j];
            }
            factors += c->lp_r_num_words;
            // assign an index back to this cofactor position in the relation_batch_t 
            t->rb_idx_r[t->numres_r] = i;
            t->numres_r++;
        }
        else if (c->lp_r[0] > 1)
        {
            max_words[0] = MAX(1, max_words[0]);
        }

        // and then the a-side
        if (c->lp_a_num_words > 0)
        {
            if (c->lp_a_num_words == 3)
            {
                if (t->first_side == 1)
                {
                    // we've already picked the a-side as first side,
                    // meaning that the r-side has 3lps as well.
                    printf("cuda cofactorization can't handle 3lp on both sides yet\n");
                    exit(1);
                }

                t->lpb_3lp = t->lpba;		// 3lp's are on the a-side
                t->lpb_2lp = t->lpbr;		// 2lp's are on the r-side
                t->first_side = 0;		// so do the r-side first.
                max_words[1] = MAX(c->lp_a_num_words, max_words[1]);
            }
            for (j = 0; j < c->lp_a_num_words; j++)
            {
                t->residues_a_in[ka++] = factors[j];
            }
            factors += c->lp_a_num_words;
            // assign an index back to this cofactor position in the relation_batch_t 
            t->rb_idx_a[t->numres_a] = i;
            t->numres_a++;
        }
        else if (c->lp_a[0] > 1)
        {
            max_words[1] = MAX(1, max_words[1]);
        }
    }

    if (t->first_side < 0)
    {
        printf("could not determine first side to factor\n");
        exit(1);
    }

    if (max_words[0] == 0)
    {
        printf("only one side has unfactored residues of max size %d\n", max_words[1]);
        if (max_words[1] == 2)
            t->first_side = -2;
        else if (max_words[1] == 3)
            t->first_side = -3;
        else
        {
            printf("could not find a list of unfactored residues to process\n");
            exit(1);
        }
    }

    if (max_words[1] == 0)
    {
        printf("only one side has unfactored residues of max size %d\n", max_words[0]);
        if (max_words[0] == 2)
            t->first_side = -4;
        else if (max_words[0] == 3)
            t->first_side = -5;
        else
        {
            printf("could not find a list of unfactored residues to process\n");
            exit(1);
        }
    }

    // determine ECM parameters from LPB sizes.
    if (t->lpb_2lp <= 26)
    {
        t->b1_2lp = 85;
        t->curves_2lp = 32;
    }
    else if (t->lpb_2lp <= 28)
    {
        t->b1_2lp = 125;
        t->curves_2lp = 32;
    }
    else if (t->lpb_2lp <= 30)
    {
        t->b1_2lp = 165;
        t->curves_2lp = 40;
    }
    else // <= 32
    {
        t->b1_2lp = 205;
        t->curves_2lp = 40;
    }

    if (t->lpb_3lp <= 26)
    {
        t->b1_3lp = 85;
        t->curves_3lp = 64;
    }
    else if (t->lpb_3lp <= 28)
    {
        t->b1_3lp = 125;
        t->curves_3lp = 64;
    }
    else if (t->lpb_3lp <= 30)
    {
        t->b1_3lp = 165;
        t->curves_3lp = 80;
    }
    else // <= 32
    {
        t->b1_3lp = 205;
        t->curves_3lp = 80;
    }

    if (b1_3lp_ovr > 0) t->b1_3lp = b1_3lp_ovr;
    if (b1_2lp_ovr > 0) t->b1_2lp = b1_2lp_ovr;
    if (b2_3lp_ovr > 0) t->b2_3lp = b2_3lp_ovr;
    if (b2_2lp_ovr > 0) t->b2_2lp = b2_2lp_ovr;
    if (curves_3lp_ovr > 0) t->curves_3lp = curves_3lp_ovr;
    if (curves_2lp_ovr > 0) t->curves_2lp = curves_2lp_ovr;

/* Helper macro so each allocation is one line and exits on failure */
#define GPU_ALLOC(field, nbytes) \
    do { \
        t->field = clCreateBuffer(t->gpu_context, CL_MEM_READ_WRITE, \
                                  (nbytes), NULL, &err); \
        if (err != CL_SUCCESS) { \
            fprintf(stderr, "clCreateBuffer(" #field ") failed: %s\n", \
                    clGetErrorString(err)); \
            exit(-1); \
        } \
    } while (0)

    /* Device buffer allocation -- mirrors the cuMemAlloc block exactly */
    GPU_ALLOC(gpu_u32_array,   sizeof(uint32_t) *     t->array_sz);
    GPU_ALLOC(gpu_rsq_array,   sizeof(uint32_t) * 3 * t->array_sz);
    GPU_ALLOC(gpu_a_array,     sizeof(uint64_t) *     t->array_sz);
    GPU_ALLOC(gpu_one_array,   sizeof(uint32_t) * 3 * t->array_sz);
    GPU_ALLOC(gpu_n_array,     sizeof(uint32_t) * 3 * t->array_sz);
    GPU_ALLOC(gpu_res32_array, sizeof(uint32_t) * 3 * t->array_sz);
    GPU_ALLOC(gpu_rho_array,   sizeof(uint32_t) *     t->array_sz);

#undef GPU_ALLOC

    /* Host buffer allocation -- unchanged */
    t->u32_array   = (uint32_t *)xmalloc(sizeof(uint32_t) * t->array_sz);
    t->rsq         = (uint64_t *)xmalloc(sizeof(uint64_t) * t->array_sz);
    t->modulus_in  = (uint64_t *)xmalloc(sizeof(uint64_t) * t->array_sz);
    t->one         = (uint64_t *)xmalloc(sizeof(uint64_t) * t->array_sz);
    t->a           = (uint64_t *)xmalloc(sizeof(uint64_t) * t->array_sz);

    t->rsq96        = (uint32_t *)xmalloc(sizeof(uint32_t) * 3 * t->array_sz);
    t->modulus96_in = (uint32_t *)xmalloc(sizeof(uint32_t) * 3 * t->array_sz);
    t->one96        = (uint32_t *)xmalloc(sizeof(uint32_t) * 3 * t->array_sz);
    t->factors96    = (uint32_t *)xmalloc(sizeof(uint32_t) * 3 * t->array_sz);

    t->rho = (uint32_t *)xmalloc(sizeof(uint32_t) * t->array_sz);

    /* generate a sigma for each input */
    for (i = 0; i < t->array_sz; i++)
    {
        t->u32_array[i] = uecm_lcg_rand_32B(7, 0xffffffff, lcg);
    }

    gpu_cofactorization(t);

    /* Host memory cleanup -- unchanged */
    free(t->a);
    free(t->u32_array);
    free(t->modulus_in);
    free(t->rho);
    free(t->rsq);
    free(t->one);

    free(t->residues_r_in);
    free(t->residues_a_in);
    free(t->rb_idx_r);
    free(t->rb_idx_a);
    free(t->rb_idx_2lp);
    free(t->rb_idx_3lp);

    free(t->modulus96_in);
    free(t->rsq96);
    free(t->one96);
    free(t->factors96);

    t->params = free_tinysiqs(t->params);

    /* Device buffer cleanup -- replaces cuMemFree */
    clReleaseMemObject(t->gpu_a_array);
    clReleaseMemObject(t->gpu_u32_array);
    clReleaseMemObject(t->gpu_n_array);
    clReleaseMemObject(t->gpu_rho_array);
    clReleaseMemObject(t->gpu_rsq_array);
    clReleaseMemObject(t->gpu_one_array);
    clReleaseMemObject(t->gpu_res32_array);

    return 0;
}

#endif /* HAVE_CUDA_BATCH_FACTOR */
