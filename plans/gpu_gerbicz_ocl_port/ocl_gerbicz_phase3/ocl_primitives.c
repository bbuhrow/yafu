/*--------------------------------------------------------------------
 * ocl_primitives.c  --  see ocl_primitives.h for the design rationale.
 *--------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ocl_primitives.h"

#define NUM_PRIM_KERNELS 8

static uint32_t
ceil_div_u32(uint32_t a, uint32_t b)
{
    return (a + b - 1u) / b;
}

static uint32_t
round_up_u32(uint32_t a, uint32_t b)
{
    return ceil_div_u32(a, b) * b;
}

static int
ensure_u32_cap(ocl_gerbicz_device_t *gd, cl_mem *buf, size_t *cap, size_t need_elems)
{
    cl_int err;
    if (*cap >= need_elems)
        return 0;
    if (*buf) {
        clReleaseMemObject(*buf);
        *buf = NULL;
        *cap = 0;
    }
    *buf = clCreateBuffer(gd->dev.context, CL_MEM_READ_WRITE,
            need_elems * sizeof(uint32_t), NULL, &err);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "ensure_u32_cap: clCreateBuffer failed: %s\n", clGetErrorString(err));
        return -1;
    }
    *cap = need_elems;
    return 0;
}

static int
ensure_u64_cap(ocl_gerbicz_device_t *gd, cl_mem *buf, size_t *cap, size_t need_elems)
{
    cl_int err;
    if (*cap >= need_elems)
        return 0;
    if (*buf) {
        clReleaseMemObject(*buf);
        *buf = NULL;
        *cap = 0;
    }
    *buf = clCreateBuffer(gd->dev.context, CL_MEM_READ_WRITE,
            need_elems * sizeof(uint64_t), NULL, &err);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "ensure_u64_cap: clCreateBuffer failed: %s\n", clGetErrorString(err));
        return -1;
    }
    *cap = need_elems;
    return 0;
}

/* ---------------------------------------------------------------------
 * init / free
 * --------------------------------------------------------------------- */
int
ocl_primitives_init(ocl_primitives_t *p, ocl_gerbicz_device_t *gd,
                     const char *cl_source_path, const char *cache_dir)
{
    char *src;
    const char *src_ptr;
    char extra_opts[256];
    char subgroup_opt[32];
    static const char *kernel_names[NUM_PRIM_KERNELS] = {
        "ocl_fill_u32", "ocl_fill_u64", "ocl_reduce_max_u32",
        "ocl_scan_local_u32", "ocl_scan_single_wg_u32", "ocl_scan_addback_u32",
        "ocl_radix_histogram_u64", "ocl_radix_scatter_u64"
    };
    static const gpu_arg_type_list_t arg_descs[NUM_PRIM_KERNELS] = {
        /* ocl_fill_u32: buf, buf_off, n, value */
        { 4, { GPU_ARG_PTR, GPU_ARG_UINT64, GPU_ARG_UINT32, GPU_ARG_UINT32 } },
        /* ocl_fill_u64: buf, buf_off, n, value */
        { 4, { GPU_ARG_PTR, GPU_ARG_UINT64, GPU_ARG_UINT32, GPU_ARG_UINT64 } },
        /* ocl_reduce_max_u32: in, in_off, out, out_off, n */
        { 5, { GPU_ARG_PTR, GPU_ARG_UINT64, GPU_ARG_PTR, GPU_ARG_UINT64, GPU_ARG_UINT32 } },
        /* ocl_scan_local_u32: in, in_off, out, out_off, block_sums, block_sums_off, n */
        { 7, { GPU_ARG_PTR, GPU_ARG_UINT64, GPU_ARG_PTR, GPU_ARG_UINT64,
               GPU_ARG_PTR, GPU_ARG_UINT64, GPU_ARG_UINT32 } },
        /* ocl_scan_single_wg_u32: buf, buf_off, n */
        { 3, { GPU_ARG_PTR, GPU_ARG_UINT64, GPU_ARG_UINT32 } },
        /* ocl_scan_addback_u32: buf, buf_off, block_sums, block_sums_off, n */
        { 5, { GPU_ARG_PTR, GPU_ARG_UINT64, GPU_ARG_PTR, GPU_ARG_UINT64, GPU_ARG_UINT32 } },
        /* ocl_radix_histogram_u64: keys, keys_off, hist, hist_off, n, num_wg, shift */
        { 7, { GPU_ARG_PTR, GPU_ARG_UINT64, GPU_ARG_PTR, GPU_ARG_UINT64,
               GPU_ARG_UINT32, GPU_ARG_UINT32, GPU_ARG_UINT32 } },
        /* ocl_radix_scatter_u64: keys_in, keys_in_off, keys_out, keys_out_off,
           offsets, offsets_off, n, num_wg, shift */
        { 9, { GPU_ARG_PTR, GPU_ARG_UINT64, GPU_ARG_PTR, GPU_ARG_UINT64,
               GPU_ARG_PTR, GPU_ARG_UINT64, GPU_ARG_UINT32, GPU_ARG_UINT32, GPU_ARG_UINT32 } },
    };

    memset(p, 0, sizeof(*p));
    p->gd = gd;

    src = ocl_read_text_file(cl_source_path);
    if (!src) {
        fprintf(stderr, "ocl_primitives_init: cannot read '%s'\n", cl_source_path);
        return -1;
    }
    src_ptr = src;

    snprintf(subgroup_opt, sizeof(subgroup_opt), "-DHAVE_SUBGROUPS=%d",
              gd->dev.has_subgroups ? 1 : 0);
    snprintf(extra_opts, sizeof(extra_opts),
             "%s -DPRIM_LOCAL_SIZE=%uu -DPRIM_ELEMS_PER_THREAD=%uu "
             "-DPRIM_ELEMS_PER_WG=%uu -DPRIM_RADIX_BITS=%uu -DPRIM_RADIX_BINS=%uu",
             subgroup_opt,
             (unsigned)OCL_PRIM_LOCAL_SIZE, (unsigned)OCL_PRIM_ELEMS_PER_THREAD,
             (unsigned)OCL_PRIM_ELEMS_PER_WG, (unsigned)OCL_PRIM_RADIX_BITS,
             (unsigned)OCL_PRIM_RADIX_BINS);

    p->program = ocl_build_program_cached(&gd->dev, &src_ptr, 1,
            gd->dev.cl_c_std[0] ? gd->dev.cl_c_std : "CL2.0",
            extra_opts, cache_dir, "gerbicz_primitives");
    free(src);
    if (!p->program)
        return -1;

    if (ocl_thread_init(&p->th, &gd->dev, p->program,
                         kernel_names, arg_descs, NUM_PRIM_KERNELS) != 0) {
        clReleaseProgram(p->program);
        p->program = NULL;
        return -1;
    }

    p->k_fill_u32            = 0;
    p->k_fill_u64             = 1;
    p->k_reduce_max_u32       = 2;
    p->k_scan_local_u32       = 3;
    p->k_scan_single_wg_u32   = 4;
    p->k_scan_addback_u32     = 5;
    p->k_radix_histogram_u64  = 6;
    p->k_radix_scatter_u64    = 7;
    return 0;
}

void
ocl_primitives_free(ocl_primitives_t *p)
{
    if (p->scratch_a) clReleaseMemObject(p->scratch_a);
    if (p->scratch_b) clReleaseMemObject(p->scratch_b);
    if (p->scratch_c) clReleaseMemObject(p->scratch_c);
    if (p->scratch_keys) clReleaseMemObject(p->scratch_keys);
    ocl_thread_free(&p->th);
    if (p->program) clReleaseProgram(p->program);
    memset(p, 0, sizeof(*p));
}

/* ---------------------------------------------------------------------
 * fill
 * --------------------------------------------------------------------- */
int
ocl_fill_u32(ocl_primitives_t *p, cl_mem buf, cl_ulong buf_off, uint32_t n, uint32_t value)
{
    gpu_arg_t args[4];
    size_t global, local = OCL_PRIM_LOCAL_SIZE;
    cl_int err;

    if (n == 0) return 0;
    args[0].ptr_arg = buf;      args[1].uint64_arg = buf_off;
    args[2].uint32_arg = n;     args[3].uint32_arg = value;
    gpu_launch_set(&p->th.launch[p->k_fill_u32], args);

    global = round_up_u32(n, (uint32_t)local);
    err = clEnqueueNDRangeKernel(p->th.queue, p->th.launch[p->k_fill_u32].kernel_func,
            1, NULL, &global, &local, 0, NULL, NULL);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "ocl_fill_u32: launch failed: %s\n", clGetErrorString(err));
        return -1;
    }
    return 0;
}

int
ocl_fill_u64(ocl_primitives_t *p, cl_mem buf, cl_ulong buf_off, uint32_t n, uint64_t value)
{
    gpu_arg_t args[4];
    size_t global, local = OCL_PRIM_LOCAL_SIZE;
    cl_int err;

    if (n == 0) return 0;
    args[0].ptr_arg = buf;      args[1].uint64_arg = buf_off;
    args[2].uint32_arg = n;     args[3].uint64_arg = value;
    gpu_launch_set(&p->th.launch[p->k_fill_u64], args);

    global = round_up_u32(n, (uint32_t)local);
    err = clEnqueueNDRangeKernel(p->th.queue, p->th.launch[p->k_fill_u64].kernel_func,
            1, NULL, &global, &local, 0, NULL, NULL);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "ocl_fill_u64: launch failed: %s\n", clGetErrorString(err));
        return -1;
    }
    return 0;
}

/* ---------------------------------------------------------------------
 * reduce-max
 * --------------------------------------------------------------------- */
int
ocl_reduce_max_u32(ocl_primitives_t *p, cl_mem in, cl_ulong in_off,
                    uint32_t n, uint32_t *out_max)
{
    cl_mem cur_buf = in;
    cl_ulong cur_off = in_off;
    uint32_t cur_n = n;
    int use_a_as_dest = 1;
    cl_int err;

    if (n == 0) {
        fprintf(stderr, "ocl_reduce_max_u32: n == 0 has no defined result\n");
        return -1;
    }

    while (cur_n > 1) {
        uint32_t num_wg = ceil_div_u32(cur_n, OCL_PRIM_ELEMS_PER_WG);
        cl_mem dest;
        gpu_arg_t args[5];
        size_t global = (size_t)num_wg * OCL_PRIM_LOCAL_SIZE;
        size_t local = OCL_PRIM_LOCAL_SIZE;

        if (use_a_as_dest) {
            if (ensure_u32_cap(p->gd, &p->scratch_a, &p->scratch_a_cap, num_wg) != 0)
                return -1;
            dest = p->scratch_a;
        } else {
            if (ensure_u32_cap(p->gd, &p->scratch_b, &p->scratch_b_cap, num_wg) != 0)
                return -1;
            dest = p->scratch_b;
        }

        args[0].ptr_arg = cur_buf;  args[1].uint64_arg = cur_off;
        args[2].ptr_arg = dest;     args[3].uint64_arg = 0;
        args[4].uint32_arg = cur_n;
        gpu_launch_set(&p->th.launch[p->k_reduce_max_u32], args);

        err = clEnqueueNDRangeKernel(p->th.queue, p->th.launch[p->k_reduce_max_u32].kernel_func,
                1, NULL, &global, &local, 0, NULL, NULL);
        if (err != CL_SUCCESS) {
            fprintf(stderr, "ocl_reduce_max_u32: launch failed: %s\n", clGetErrorString(err));
            return -1;
        }

        cur_buf = dest;
        cur_off = 0;
        cur_n = num_wg;
        use_a_as_dest = !use_a_as_dest;
    }

    err = clEnqueueReadBuffer(p->th.queue, cur_buf, CL_TRUE,
            (size_t)cur_off * sizeof(uint32_t), sizeof(uint32_t), out_max, 0, NULL, NULL);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "ocl_reduce_max_u32: readback failed: %s\n", clGetErrorString(err));
        return -1;
    }
    return 0;
}

/* ---------------------------------------------------------------------
 * exclusive scan
 * --------------------------------------------------------------------- */
int
ocl_scan_exclusive_u32(ocl_primitives_t *p, cl_mem buf, cl_ulong buf_off, uint32_t n)
{
    uint32_t num_blocks;
    gpu_arg_t args[7];
    size_t global, local = OCL_PRIM_LOCAL_SIZE;
    cl_int err;

    if (n == 0)
        return 0;
    if (n > OCL_PRIM_SCAN_MAX_N) {
        fprintf(stderr, "ocl_scan_exclusive_u32: n=%u exceeds two-level scan capacity %u\n",
                n, (unsigned)OCL_PRIM_SCAN_MAX_N);
        return -1;
    }

    num_blocks = ceil_div_u32(n, OCL_PRIM_ELEMS_PER_WG);
    if (ensure_u32_cap(p->gd, &p->scratch_c, &p->scratch_c_cap, num_blocks) != 0)
        return -1;

    /* Level 1: local scan, in-place (in==out is safe -- see header comment
     * in ocl_primitives.h / the kernel comment: each thread only ever
     * reads the exact indices it later writes, no cross-thread reads). */
    args[0].ptr_arg = buf;          args[1].uint64_arg = buf_off;
    args[2].ptr_arg = buf;          args[3].uint64_arg = buf_off;
    args[4].ptr_arg = p->scratch_c; args[5].uint64_arg = 0;
    args[6].uint32_arg = n;
    gpu_launch_set(&p->th.launch[p->k_scan_local_u32], args);
    global = (size_t)num_blocks * OCL_PRIM_LOCAL_SIZE;
    err = clEnqueueNDRangeKernel(p->th.queue, p->th.launch[p->k_scan_local_u32].kernel_func,
            1, NULL, &global, &local, 0, NULL, NULL);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "ocl_scan_exclusive_u32: local scan launch failed: %s\n", clGetErrorString(err));
        return -1;
    }

    /* Level 2: scan the (small) block-sums array in place, one workgroup. */
    {
        gpu_arg_t args2[3];
        size_t g2 = OCL_PRIM_LOCAL_SIZE;
        args2[0].ptr_arg = p->scratch_c; args2[1].uint64_arg = 0;
        args2[2].uint32_arg = num_blocks;
        gpu_launch_set(&p->th.launch[p->k_scan_single_wg_u32], args2);
        err = clEnqueueNDRangeKernel(p->th.queue, p->th.launch[p->k_scan_single_wg_u32].kernel_func,
                1, NULL, &g2, &local, 0, NULL, NULL);
        if (err != CL_SUCCESS) {
            fprintf(stderr, "ocl_scan_exclusive_u32: block-sums scan launch failed: %s\n", clGetErrorString(err));
            return -1;
        }
    }

    /* Level 3: add block prefix back into every element. */
    {
        gpu_arg_t args3[5];
        size_t g3 = round_up_u32(n, (uint32_t)local);
        args3[0].ptr_arg = buf;          args3[1].uint64_arg = buf_off;
        args3[2].ptr_arg = p->scratch_c; args3[3].uint64_arg = 0;
        args3[4].uint32_arg = n;
        gpu_launch_set(&p->th.launch[p->k_scan_addback_u32], args3);
        err = clEnqueueNDRangeKernel(p->th.queue, p->th.launch[p->k_scan_addback_u32].kernel_func,
                1, NULL, &g3, &local, 0, NULL, NULL);
        if (err != CL_SUCCESS) {
            fprintf(stderr, "ocl_scan_exclusive_u32: add-back launch failed: %s\n", clGetErrorString(err));
            return -1;
        }
    }

    return 0;
}

/* ---------------------------------------------------------------------
 * LSD radix sort (ulong keys)
 * --------------------------------------------------------------------- */
int
ocl_radix_sort_u64(ocl_primitives_t *p,
                    cl_mem keys_in,  cl_ulong keys_in_off,
                    cl_mem keys_out, cl_ulong keys_out_off,
                    uint32_t n)
{
    uint32_t num_wg, hist_size, pass;
    cl_mem cur_src; cl_ulong cur_src_off;
    cl_mem bufs[2]; cl_ulong offs[2];
    int dst_idx;
    cl_int err;

    if (n == 0)
        return 0;

    num_wg = ceil_div_u32(n, OCL_PRIM_ELEMS_PER_WG);
    hist_size = (uint32_t)((uint64_t)num_wg * OCL_PRIM_RADIX_BINS);
    if (hist_size > OCL_PRIM_SCAN_MAX_N) {
        fprintf(stderr, "ocl_radix_sort_u64: n=%u needs a %u-entry histogram, "
                "exceeds scan capacity %u\n", n, hist_size, (unsigned)OCL_PRIM_SCAN_MAX_N);
        return -1;
    }
    if (ensure_u32_cap(p->gd, &p->scratch_a, &p->scratch_a_cap, hist_size) != 0)
        return -1;
    if (ensure_u64_cap(p->gd, &p->scratch_keys, &p->scratch_keys_cap, n) != 0)
        return -1;

    cur_src = keys_in;
    cur_src_off = keys_in_off;
    bufs[0] = keys_out;       offs[0] = keys_out_off;
    bufs[1] = p->scratch_keys; offs[1] = 0;
    dst_idx = 0;

    for (pass = 0; pass < OCL_PRIM_RADIX_PASSES_U64; pass++) {
        uint32_t shift = pass * OCL_PRIM_RADIX_BITS;
        cl_mem dst = bufs[dst_idx];
        cl_ulong dst_off = offs[dst_idx];
        size_t global = (size_t)num_wg * OCL_PRIM_LOCAL_SIZE;
        size_t local = OCL_PRIM_LOCAL_SIZE;
        gpu_arg_t hargs[7];
        gpu_arg_t sargs[9];

        hargs[0].ptr_arg = cur_src;      hargs[1].uint64_arg = cur_src_off;
        hargs[2].ptr_arg = p->scratch_a; hargs[3].uint64_arg = 0;
        hargs[4].uint32_arg = n;         hargs[5].uint32_arg = num_wg;
        hargs[6].uint32_arg = shift;
        gpu_launch_set(&p->th.launch[p->k_radix_histogram_u64], hargs);
        err = clEnqueueNDRangeKernel(p->th.queue, p->th.launch[p->k_radix_histogram_u64].kernel_func,
                1, NULL, &global, &local, 0, NULL, NULL);
        if (err != CL_SUCCESS) {
            fprintf(stderr, "ocl_radix_sort_u64: histogram launch failed (pass %u): %s\n",
                    pass, clGetErrorString(err));
            return -1;
        }

        if (ocl_scan_exclusive_u32(p, p->scratch_a, 0, hist_size) != 0)
            return -1;

        sargs[0].ptr_arg = cur_src;      sargs[1].uint64_arg = cur_src_off;
        sargs[2].ptr_arg = dst;          sargs[3].uint64_arg = dst_off;
        sargs[4].ptr_arg = p->scratch_a; sargs[5].uint64_arg = 0;
        sargs[6].uint32_arg = n;         sargs[7].uint32_arg = num_wg;
        sargs[8].uint32_arg = shift;
        gpu_launch_set(&p->th.launch[p->k_radix_scatter_u64], sargs);
        err = clEnqueueNDRangeKernel(p->th.queue, p->th.launch[p->k_radix_scatter_u64].kernel_func,
                1, NULL, &global, &local, 0, NULL, NULL);
        if (err != CL_SUCCESS) {
            fprintf(stderr, "ocl_radix_sort_u64: scatter launch failed (pass %u): %s\n",
                    pass, clGetErrorString(err));
            return -1;
        }

        cur_src = dst;
        cur_src_off = dst_off;
        dst_idx ^= 1;
    }

    if (cur_src != keys_out || cur_src_off != keys_out_off) {
        err = clEnqueueCopyBuffer(p->th.queue, cur_src, keys_out,
                (size_t)cur_src_off * sizeof(uint64_t),
                (size_t)keys_out_off * sizeof(uint64_t),
                (size_t)n * sizeof(uint64_t), 0, NULL, NULL);
        if (err != CL_SUCCESS) {
            fprintf(stderr, "ocl_radix_sort_u64: final copy failed: %s\n", clGetErrorString(err));
            return -1;
        }
    }
    clFinish(p->th.queue);
    return 0;
}
