/* test_ocl_foundation.c -- Phase 1 task 1.7 smoke test.
 *
 * Runs the whole new pipeline end to end:
 *   gpu_init -> ocl_gerbicz_device_init -> ocl_build_program_cached
 *   -> ocl_thread_init -> gpu_launch_set -> clEnqueueNDRangeKernel
 *   -> verify -> re-run to prove the disk cache path also works.
 *
 * Also exercises the base-offset convention: the two input buffers are
 * allocated bigger than needed and read from a non-zero offset, and the
 * output buffer is written at a non-zero offset, so a convention bug
 * (e.g. forgetting to add the offset in the kernel) shows up as wrong
 * output rather than silently working.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ocl_gerbicz_ctx.h"

#define N 1024
#define A_OFF 7
#define B_OFF 13
#define OUT_OFF 3

static int
run_once(ocl_gerbicz_device_t *gd, const char *label)
{
    char build_opts[192];
    const char *src;
    char *src_owned;
    cl_int err;
    ocl_thread_t th;
    const char *kname = "smoke_vecadd";
    gpu_arg_type_list_t arg_desc = {
        7, { GPU_ARG_PTR, GPU_ARG_UINT64,
             GPU_ARG_PTR, GPU_ARG_UINT64,
             GPU_ARG_PTR, GPU_ARG_UINT64,
             GPU_ARG_UINT32 }
    };
    cl_mem buf_a, buf_b, buf_out;
    uint32_t *h_a, *h_b, *h_out;
    unsigned i;
    size_t global;
    gpu_arg_t args[7];
    int ok = 1;

    ocl_gerbicz_build_opts(gd, build_opts, sizeof(build_opts));
    printf("[%s] build options: -cl-std=%s %s\n", label,
           gd->dev.cl_c_std[0] ? gd->dev.cl_c_std : "CL2.0", build_opts);

    src_owned = ocl_read_text_file("ocl_gerbicz_smoke.cl");
    if (!src_owned) { fprintf(stderr, "cannot open ocl_gerbicz_smoke.cl\n"); exit(1); }
    src = src_owned;

    gd->program_collision = ocl_build_program_cached(&gd->dev, &src, 1,
            gd->dev.cl_c_std[0] ? gd->dev.cl_c_std : "CL2.0",
            build_opts, ".", "gerbicz_smoke");
    free(src_owned);
    if (!gd->program_collision) {
        fprintf(stderr, "[%s] program build failed\n", label);
        return 0;
    }

    if (ocl_thread_init(&th, &gd->dev, gd->program_collision,
                         &kname, &arg_desc, 1) != 0) {
        fprintf(stderr, "[%s] ocl_thread_init failed\n", label);
        return 0;
    }
    printf("[%s] preferred work-group multiple: %zu\n",
           label, gd->dev.pref_wg_multiple);

    h_a = malloc((N + A_OFF) * sizeof(uint32_t));
    h_b = malloc((N + B_OFF) * sizeof(uint32_t));
    h_out = calloc(N + OUT_OFF, sizeof(uint32_t));
    for (i = 0; i < N; i++) {
        h_a[A_OFF + i] = i;
        h_b[B_OFF + i] = 2 * i + 1;
    }

    buf_a = clCreateBuffer(gd->dev.context, CL_MEM_READ_ONLY,
            (N + A_OFF) * sizeof(uint32_t), NULL, &err);
    buf_b = clCreateBuffer(gd->dev.context, CL_MEM_READ_ONLY,
            (N + B_OFF) * sizeof(uint32_t), NULL, &err);
    buf_out = clCreateBuffer(gd->dev.context, CL_MEM_READ_WRITE,
            (N + OUT_OFF) * sizeof(uint32_t), NULL, &err);

    clEnqueueWriteBuffer(th.queue, buf_a, CL_TRUE, 0,
            (N + A_OFF) * sizeof(uint32_t), h_a, 0, NULL, NULL);
    clEnqueueWriteBuffer(th.queue, buf_b, CL_TRUE, 0,
            (N + B_OFF) * sizeof(uint32_t), h_b, 0, NULL, NULL);

    args[0].ptr_arg = buf_a;    args[1].uint64_arg = A_OFF;
    args[2].ptr_arg = buf_b;    args[3].uint64_arg = B_OFF;
    args[4].ptr_arg = buf_out;  args[5].uint64_arg = OUT_OFF;
    args[6].uint32_arg = N;
    gpu_launch_set(&th.launch[0], args);

    global = ((N + 63) / 64) * 64;
    err = clEnqueueNDRangeKernel(th.queue, th.launch[0].kernel_func,
            1, NULL, &global, NULL, 0, NULL, NULL);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "[%s] clEnqueueNDRangeKernel failed: %s\n",
                label, clGetErrorString(err));
        ok = 0;
    } else {
        clFinish(th.queue);
        clEnqueueReadBuffer(th.queue, buf_out, CL_TRUE, 0,
                (N + OUT_OFF) * sizeof(uint32_t), h_out, 0, NULL, NULL);
        for (i = 0; i < N; i++) {
            uint32_t expect = i + (2 * i + 1);
            if (h_out[OUT_OFF + i] != expect) {
                fprintf(stderr, "[%s] mismatch at %u: got %u want %u\n",
                        label, i, h_out[OUT_OFF + i], expect);
                ok = 0;
                break;
            }
        }
        if (ok)
            printf("[%s] %d elements verified OK (base-offset convention exercised: "
                   "a_off=%d b_off=%d out_off=%d)\n", label, N, A_OFF, B_OFF, OUT_OFF);
    }

    clReleaseMemObject(buf_a);
    clReleaseMemObject(buf_b);
    clReleaseMemObject(buf_out);
    free(h_a); free(h_b); free(h_out);
    ocl_thread_free(&th);
    return ok;
}

int
main(void)
{
    gpu_config_t config;
    ocl_gerbicz_device_t gd;
    int ok;

    gpu_init(&config);
    if (config.num_gpu == 0) {
        fprintf(stderr, "no OpenCL devices found -- build-only mode not implemented "
                        "in this driver; see STATUS\n");
        return 1;
    }

    printf("found %d OpenCL device(s); using device 0: %s (OpenCL %d.%d)\n",
           config.num_gpu, config.info[0].name,
           config.info[0].compute_version_major,
           config.info[0].compute_version_minor);

    if (ocl_gerbicz_device_init(&gd, &config.info[0]) != 0)
        return 1;

    printf("local_mem_size=%llu max_mem_alloc_size=%llu max_work_group_size=%zu "
           "has_subgroups=%d\n",
           (unsigned long long)gd.dev.local_mem_size,
           (unsigned long long)gd.dev.max_mem_alloc_size,
           gd.dev.max_work_group_size, gd.dev.has_subgroups);

    /* First run: compiles from source and writes the cache. */
    remove("gerbicz_smoke_*.bin"); /* best-effort; glob doesn't work in remove(),
                                       real cleanup is in the Makefile's `clean` */
    ok = run_once(&gd, "compile");
    if (gd.program_collision) { clReleaseProgram(gd.program_collision); gd.program_collision = NULL; }

    /* Second run: must hit the disk cache (prints "loaded cached binary"). */
    ok = ok && run_once(&gd, "cached");

    ocl_gerbicz_device_free(&gd);

    printf(ok ? "PHASE 1 SMOKE TEST: PASS\n" : "PHASE 1 SMOKE TEST: FAIL\n");
    return ok ? 0 : 1;
}
