/*--------------------------------------------------------------------
 * ocl_collision.c -- see ocl_collision.h for design notes.
 *--------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ocl_collision.h"
#include "collision_bucket.h"   /* NUM_BUCKETS, LOG2_NUM_BUCKETS, BUCKET_HASH_MIX */

#define NUM_COLL_KERNELS 10

/* Portable ilog2-ish helpers (collision_engine.cu's host_ilog2 uses
 * __builtin_clz directly, which is GCC/clang-only; kept portable here
 * since this file must also build under MSVC eventually). */
static uint32_t
msb_index(uint32_t v)
{
#if defined(__GNUC__) || defined(__clang__)
    return 31u - (uint32_t)__builtin_clz(v);
#else
    uint32_t r = 0;
    while (v >>= 1) r++;
    return r;
#endif
}

static uint32_t
host_ilog2(uint32_t cnt)
{
    uint32_t lg, il;
    if (cnt <= 1u)
        return 5u;
    lg = msb_index(cnt);
    il = lg + 5u;
    return il < 5u ? 5u : il;
}

static uint32_t
ceil_div_u32(uint32_t a, uint32_t b) { return (a + b - 1u) / b; }

static size_t
round_up_to(uint32_t n, uint32_t multiple)
{
    return (size_t)ceil_div_u32(n, multiple) * multiple;
}

static int
ensure_u32(ocl_gerbicz_device_t *gd, cl_mem *buf, size_t n_elems)
{
    cl_int err;
    *buf = clCreateBuffer(gd->dev.context, CL_MEM_READ_WRITE,
            n_elems * sizeof(uint32_t), NULL, &err);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "ocl_collision_init: alloc failed (%zu u32): %s\n",
                n_elems, clGetErrorString(err));
        return -1;
    }
    return 0;
}

static int
ensure_u64(ocl_gerbicz_device_t *gd, cl_mem *buf, size_t n_elems)
{
    cl_int err;
    *buf = clCreateBuffer(gd->dev.context, CL_MEM_READ_WRITE,
            n_elems * sizeof(uint64_t), NULL, &err);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "ocl_collision_init: alloc failed (%zu u64): %s\n",
                n_elems, clGetErrorString(err));
        return -1;
    }
    return 0;
}

/* ---------------------------------------------------------------------
 * init / free
 * --------------------------------------------------------------------- */
int
ocl_collision_init(ocl_collision_engine_t *e, ocl_gerbicz_device_t *gd,
                    ocl_primitives_t *prim,
                    const char *cl_source_path, const char *cache_dir)
{
    char *src; const char *src_ptr;
    char extra_opts[384];
    static const char *knames[NUM_COLL_KERNELS] = {
        "ocl_scatter_roots", "ocl_filter_per_bucket", "ocl_dedup",
        "ocl_count_secondary", "ocl_scatter_secondary",
        "ocl_count_matched_values", "ocl_scatter_matched_values",
        "ocl_count_and_store_matched_values",
        "ocl_emit_found", "ocl_emit_found_arena"
    };
    static const gpu_arg_type_list_t adescs[NUM_COLL_KERNELS] = {
        /* ocl_scatter_roots: roots,roots_off,root_bytes,n,max_per_bucket,
           hash_mode,bucket_count,overflow_flag,bucket_storage */
        { 9, { GPU_ARG_PTR, GPU_ARG_UINT64, GPU_ARG_UINT32, GPU_ARG_UINT32,
               GPU_ARG_UINT32, GPU_ARG_INT32, GPU_ARG_PTR, GPU_ARG_PTR, GPU_ARG_PTR } },
        /* ocl_filter_per_bucket: bucket_count_in,max_per_bucket,arr_a,arr_b,
           candidate_keys,candidate_cnt,candidate_overflow,candidate_cap,
           key_bits,max_tsize_words,iters_hist_out,s_hash(local) */
        { 12, { GPU_ARG_PTR, GPU_ARG_UINT32, GPU_ARG_PTR, GPU_ARG_PTR,
                GPU_ARG_PTR, GPU_ARG_PTR, GPU_ARG_PTR, GPU_ARG_UINT32,
                GPU_ARG_UINT32, GPU_ARG_UINT32, GPU_ARG_PTR, GPU_ARG_LOCAL } },
        /* ocl_dedup: sorted_keys,csize,dedup_keys,dedup_cnt */
        { 4, { GPU_ARG_PTR, GPU_ARG_UINT32, GPU_ARG_PTR, GPU_ARG_PTR } },
        /* ocl_count_secondary: dedup_keys,csize,hash_value,D,S */
        { 5, { GPU_ARG_PTR, GPU_ARG_UINT32, GPU_ARG_UINT32, GPU_ARG_PTR, GPU_ARG_PTR } },
        /* ocl_scatter_secondary: dedup_keys,csize,hash_value,D_pos,X */
        { 5, { GPU_ARG_PTR, GPU_ARG_UINT32, GPU_ARG_UINT32, GPU_ARG_PTR, GPU_ARG_PTR } },
        /* ocl_count_matched_values: roots,roots_off,root_bytes,values,values_off,
           n,hash_value,S,D,X,dedup_count,value_counts,match_total,value_overflow,value_cap */
        { 15, { GPU_ARG_PTR, GPU_ARG_UINT64, GPU_ARG_UINT32, GPU_ARG_PTR, GPU_ARG_UINT64,
                GPU_ARG_UINT32, GPU_ARG_UINT32, GPU_ARG_PTR, GPU_ARG_PTR, GPU_ARG_PTR,
                GPU_ARG_UINT32, GPU_ARG_PTR, GPU_ARG_PTR, GPU_ARG_PTR, GPU_ARG_UINT32 } },
        /* ocl_scatter_matched_values: roots,roots_off,root_bytes,values,values_off,
           n,hash_value,S,D,X,dedup_count,value_cursor,matched_values,value_cap */
        { 14, { GPU_ARG_PTR, GPU_ARG_UINT64, GPU_ARG_UINT32, GPU_ARG_PTR, GPU_ARG_UINT64,
                GPU_ARG_UINT32, GPU_ARG_UINT32, GPU_ARG_PTR, GPU_ARG_PTR, GPU_ARG_PTR,
                GPU_ARG_UINT32, GPU_ARG_PTR, GPU_ARG_PTR, GPU_ARG_UINT32 } },
        /* ocl_count_and_store_matched_values: roots,roots_off,root_bytes,values,values_off,
           n,hash_value,S,D,X,dedup_count,value_counts,match_total,value_overflow,
           matched_values,value_cap */
        { 16, { GPU_ARG_PTR, GPU_ARG_UINT64, GPU_ARG_UINT32, GPU_ARG_PTR, GPU_ARG_UINT64,
                GPU_ARG_UINT32, GPU_ARG_UINT32, GPU_ARG_PTR, GPU_ARG_PTR, GPU_ARG_PTR,
                GPU_ARG_UINT32, GPU_ARG_PTR, GPU_ARG_PTR, GPU_ARG_PTR,
                GPU_ARG_PTR, GPU_ARG_UINT32 } },
        /* ocl_emit_found: X,matched_values,value_offsets,dedup_count,pshift,
           root_bytes,q_batch,q_batch_off,found_array,found_array_off */
        { 10, { GPU_ARG_PTR, GPU_ARG_PTR, GPU_ARG_PTR, GPU_ARG_UINT32, GPU_ARG_UINT32,
                GPU_ARG_UINT32, GPU_ARG_PTR, GPU_ARG_UINT64, GPU_ARG_PTR, GPU_ARG_UINT64 } },
        /* ocl_emit_found_arena: X,matched_values,value_counts,dedup_count,pshift,
           root_bytes,q_batch,q_batch_off,found_array,found_array_off */
        { 10, { GPU_ARG_PTR, GPU_ARG_PTR, GPU_ARG_PTR, GPU_ARG_UINT32, GPU_ARG_UINT32,
                GPU_ARG_UINT32, GPU_ARG_PTR, GPU_ARG_UINT64, GPU_ARG_PTR, GPU_ARG_UINT64 } },
    };

    memset(e, 0, sizeof(*e));
    e->gd = gd;
    e->prim = prim;

    src = ocl_read_text_file(cl_source_path);
    if (!src) {
        fprintf(stderr, "ocl_collision_init: cannot read '%s'\n", cl_source_path);
        return -1;
    }
    src_ptr = src;

    snprintf(extra_opts, sizeof(extra_opts),
             "-DHAVE_SUBGROUPS=%d -DLOG2_NUM_BUCKETS=%uu -DBUCKET_HASH_MIX=0x%llxULL "
             "-DCOLL_BLOCK_THREADS=%uu -DCOLL_MAX_FILTER_ITERS=%d "
             "-DCOLL_MATCH_ARENA_WIDTH=%uu -DFOUND_ARRAY_SIZE=%uu",
             gd->dev.has_subgroups ? 1 : 0,
             (unsigned)LOG2_NUM_BUCKETS, (unsigned long long)BUCKET_HASH_MIX,
             (unsigned)OCL_COLL_BLOCK_THREADS, OCL_COLL_MAX_FILTER_ITERS,
             (unsigned)OCL_COLL_MATCH_ARENA_WIDTH, (unsigned)OCL_COLL_FOUND_ARRAY_SIZE);

    e->program = ocl_build_program_cached(&gd->dev, &src_ptr, 1,
            gd->dev.cl_c_std[0] ? gd->dev.cl_c_std : "CL2.0",
            extra_opts, cache_dir, "gerbicz_collision");
    free(src);
    if (!e->program)
        return -1;

    if (ocl_thread_init(&e->th, &gd->dev, e->program, knames, adescs, NUM_COLL_KERNELS) != 0) {
        clReleaseProgram(e->program);
        e->program = NULL;
        return -1;
    }
    /* Share ONE queue with the primitives engine. ocl_thread_init
     * (Phase 1) always creates a fresh queue, but this engine and
     * ocl_primitives_t constantly touch the SAME buffers back and
     * forth (fill before scatter, reduce-max after scatter, sort then
     * dedup, scan then scatter_secondary, ...) -- two independent
     * command queues give OpenCL no reason to order those operations
     * relative to each other, which is a real race, not just an
     * efficiency concern (found the hard way: see STATUS). Decision #5
     * says one queue per THREAD, not one per "kind of kernel", so
     * discard the queue ocl_thread_init just made and reuse prim's. */
    clReleaseCommandQueue(e->th.queue);
    e->th.queue = prim->th.queue;   /* borrowed, not owned -- see ocl_collision_free */
    e->k_scatter_roots = 0;
    e->k_filter_per_bucket = 1;
    e->k_dedup = 2;
    e->k_count_secondary = 3;
    e->k_scatter_secondary = 4;
    e->k_count_matched_values = 5;
    e->k_scatter_matched_values = 6;
    e->k_count_and_store_matched_values = 7;
    e->k_emit_found = 8;
    e->k_emit_found_arena = 9;

    /* Fixed-size buffers -- allocated once, matching CUDA's own
     * struct (sized off constants, not off n/key_bits). */
    if (ensure_u32(gd, &e->d_bucket_count, NUM_BUCKETS) ||
        ensure_u32(gd, &e->d_bucket_overflow, 1) ||
        ensure_u64(gd, &e->d_candidate_keys, OCL_COLL_CANDIDATE_CAP) ||
        ensure_u64(gd, &e->d_sorted_keys, OCL_COLL_CANDIDATE_CAP) ||
        ensure_u64(gd, &e->d_dedup_keys, OCL_COLL_CANDIDATE_CAP) ||
        ensure_u32(gd, &e->d_candidate_cnt, 1) ||
        ensure_u32(gd, &e->d_candidate_overflow, 1) ||
        ensure_u32(gd, &e->d_dedup_cnt, 1) ||
        ensure_u32(gd, &e->d_D, OCL_COLL_MAX_DSIZE) ||
        ensure_u32(gd, &e->d_D_pos, OCL_COLL_MAX_DSIZE) ||
        ensure_u32(gd, &e->d_S, OCL_COLL_MAX_SSIZE) ||
        ensure_u64(gd, &e->d_X, OCL_COLL_CANDIDATE_CAP) ||
        ensure_u32(gd, &e->d_value_counts, (size_t)OCL_COLL_VALUE_MATCH_CAP + 1) ||
        ensure_u32(gd, &e->d_value_cursor, (size_t)OCL_COLL_VALUE_MATCH_CAP + 1) ||
        ensure_u32(gd, &e->d_matched_values, OCL_COLL_VALUE_MATCH_CAP) ||
        ensure_u32(gd, &e->d_value_match_total, 1) ||
        ensure_u32(gd, &e->d_value_overflow, 1) ||
        ensure_u32(gd, &e->d_filter_iters_hist, 102)) {
        ocl_collision_free(e);
        return -1;
    }

    return 0;
}

static void
release_if(cl_mem *m)
{
    if (*m) { clReleaseMemObject(*m); *m = NULL; }
}

void
ocl_collision_free(ocl_collision_engine_t *e)
{
    release_if(&e->d_arr_a); release_if(&e->d_arr_b);
    release_if(&e->d_bucket_count); release_if(&e->d_bucket_overflow);
    release_if(&e->d_candidate_keys); release_if(&e->d_sorted_keys); release_if(&e->d_dedup_keys);
    release_if(&e->d_candidate_cnt); release_if(&e->d_candidate_overflow); release_if(&e->d_dedup_cnt);
    release_if(&e->d_D); release_if(&e->d_D_pos); release_if(&e->d_S); release_if(&e->d_X);
    release_if(&e->d_value_counts); release_if(&e->d_value_cursor); release_if(&e->d_matched_values);
    release_if(&e->d_value_match_total); release_if(&e->d_value_overflow);
    release_if(&e->d_filter_iters_hist);
    e->th.queue = NULL; /* borrowed from prim's ocl_thread_t -- do not release here */
    if (e->th.launch) ocl_thread_free(&e->th);
    if (e->program) clReleaseProgram(e->program);
    memset(e, 0, sizeof(*e));
}

/* Grows d_arr_a/d_arr_b only -- mirrors collision_engine::ensure_capacity
 * exactly (same mean/sigma/estimate formula). */
static int
ensure_capacity(ocl_collision_engine_t *e, uint32_t n, uint32_t key_bits,
                 uint32_t min_per_bucket)
{
    uint32_t alloc_n, alloc_key_bits, min_bucket;
    uint32_t mean_per_bucket, sigma, max_bucket_est;
    size_t bucket_words;

    if (n <= e->max_n && key_bits <= e->max_key_bits && min_per_bucket <= e->max_per_bucket)
        return 0;

    alloc_n = n > e->max_n ? n : e->max_n;
    alloc_key_bits = key_bits > e->max_key_bits ? key_bits : e->max_key_bits;
    min_bucket = min_per_bucket > e->max_per_bucket ? min_per_bucket : e->max_per_bucket;

    release_if(&e->d_arr_a);
    release_if(&e->d_arr_b);
    e->max_n = alloc_n;
    e->max_key_bits = alloc_key_bits;

    mean_per_bucket = ceil_div_u32(e->max_n, NUM_BUCKETS);
    sigma = (uint32_t)sqrt((double)mean_per_bucket + 1.0);
    max_bucket_est = mean_per_bucket + 6u * sigma + 32u;
    e->max_per_bucket = max_bucket_est > min_bucket ? max_bucket_est : min_bucket;
    bucket_words = (size_t)NUM_BUCKETS * e->max_per_bucket;

    if (ensure_u64(e->gd, &e->d_arr_a, bucket_words) ||
        ensure_u64(e->gd, &e->d_arr_b, bucket_words))
        return -1;
    return 0;
}

/* ---------------------------------------------------------------------
 * run
 * --------------------------------------------------------------------- */
ocl_collision_status_t
ocl_collision_run(ocl_collision_engine_t *e, ocl_collision_data_t *data)
{
    uint32_t n = data->num_elements;
    uint32_t key_bits = data->key_bits;
    cl_command_queue q = e->th.queue;
    uint32_t h_max_bucket = 0, h_bucket_overflow = 0;
    uint32_t max_bucket_for_hash, lg, ilog2, hash_bits_avail, max_tsize_words;
    size_t hash_bytes;
    uint32_t cand_cnt = 0, h_candidate_overflow = 0;
    uint32_t dedup_cnt = 0;
    uint32_t c_ilog2, hash_value, dsize, ssize;
    uint32_t h_value_match_total = 0, h_value_overflow = 0;
    uint32_t use_match_arena;

    if (key_bits < OCL_COLL_MIN_KEY_BITS || key_bits > 64 ||
        (data->root_bytes != 4 && data->root_bytes != 8)) {
        fprintf(stderr, "ocl_collision_run: invalid input\n");
        return OCL_COLLISION_INVALID_INPUT;
    }

    memset(data->filter_iters_hist, 0, sizeof(data->filter_iters_hist));
    data->bucket_grow_count = 0;
    data->hash_cap_count = 0;

    if (ensure_capacity(e, n, key_bits, 0) != 0)
        return OCL_COLLISION_CL_ERROR;

    /* --- 3.3a: scatter + reduce-max, grow-and-retry loop --- */
    for (;;) {
        gpu_arg_t sargs[9];
        size_t global = round_up_to(n, 256u);
        cl_int err;

        ocl_fill_u32(e->prim, e->d_bucket_count, 0, NUM_BUCKETS, 0u);
        ocl_fill_u32(e->prim, e->d_bucket_overflow, 0, 1, 0u);

        sargs[0].ptr_arg = data->keys_in;         sargs[1].uint64_arg = data->keys_in_off;
        sargs[2].uint32_arg = data->root_bytes;   sargs[3].uint32_arg = n;
        sargs[4].uint32_arg = e->max_per_bucket;  sargs[5].int32_arg = data->bucket_hash;
        sargs[6].ptr_arg = e->d_bucket_count;     sargs[7].ptr_arg = e->d_bucket_overflow;
        sargs[8].ptr_arg = e->d_arr_a;
        gpu_launch_set(&e->th.launch[e->k_scatter_roots], sargs);
        {
            size_t local = 256;
            err = clEnqueueNDRangeKernel(q, e->th.launch[e->k_scatter_roots].kernel_func,
                    1, NULL, &global, &local, 0, NULL, NULL);
        }
        if (err != CL_SUCCESS) {
            fprintf(stderr, "ocl_collision_run: scatter_roots failed: %s\n", clGetErrorString(err));
            return OCL_COLLISION_CL_ERROR;
        }

        if (ocl_reduce_max_u32(e->prim, e->d_bucket_count, 0, NUM_BUCKETS, &h_max_bucket) != 0)
            return OCL_COLLISION_CL_ERROR;
        clEnqueueReadBuffer(q, e->d_bucket_overflow, CL_TRUE, 0, sizeof(uint32_t),
                &h_bucket_overflow, 0, NULL, NULL);

        if (!h_bucket_overflow && h_max_bucket <= e->max_per_bucket)
            break;

        {
            uint32_t observed = h_max_bucket;
            uint32_t grown_cap;
            if (observed <= e->max_per_bucket)
                observed = e->max_per_bucket + 1u;
            grown_cap = observed + observed / 4u + 64u;
            data->bucket_grow_count++;
            if (data->debug)
                fprintf(stderr, "ocl_collision: growing bucket cap %u -> %u (observed %u)\n",
                        e->max_per_bucket, grown_cap, h_max_bucket);
            if (ensure_capacity(e, n, key_bits, grown_cap) != 0)
                return OCL_COLLISION_CL_ERROR;
        }
    }
    data->bucket_max = h_max_bucket;

    /* --- hash-table sizing (collision_engine_run 930-961) --- */
    max_bucket_for_hash = h_max_bucket == 0 ? 1u : h_max_bucket;
    lg = msb_index(max_bucket_for_hash);
    ilog2 = lg + 5u;
    if (ilog2 < 5u) ilog2 = 5u;
    hash_bits_avail = key_bits - LOG2_NUM_BUCKETS;
    if (ilog2 > hash_bits_avail) ilog2 = hash_bits_avail;
    max_tsize_words = 1u << (ilog2 - 5u);

    /* CUDA queries cudaDevAttrMaxSharedMemoryPerBlockOptin at this
     * point; we already have the equivalent from Phase 1's device
     * query (ocl_device_t.local_mem_size) -- no new query needed. */
    if (e->gd->dev.local_mem_size > 0) {
        uint32_t max_words = (uint32_t)(e->gd->dev.local_mem_size / (3u * sizeof(uint32_t)));
        uint32_t capped_words = 1u;
        while ((capped_words << 1) <= max_words)
            capped_words <<= 1;
        if (max_tsize_words > capped_words) {
            data->hash_cap_count++;
            if (data->debug)
                fprintf(stderr, "ocl_collision: capping hash words %u -> %u (bucket max %u)\n",
                        max_tsize_words, capped_words, max_bucket_for_hash);
            max_tsize_words = capped_words;
        }
    }

    /* --- 3.3b: per-bucket filter --- */
    ocl_fill_u32(e->prim, e->d_candidate_cnt, 0, 1, 0u);
    ocl_fill_u32(e->prim, e->d_candidate_overflow, 0, 1, 0u);
    if (data->collect_stats)
        ocl_fill_u32(e->prim, e->d_filter_iters_hist, 0, 102, 0u);

    hash_bytes = 3u * (size_t)max_tsize_words * sizeof(uint32_t);
    {
        gpu_arg_t fargs[12];
        size_t global = (size_t)NUM_BUCKETS * OCL_COLL_BLOCK_THREADS;
        size_t local = OCL_COLL_BLOCK_THREADS;
        cl_int err;

        fargs[0].ptr_arg = e->d_bucket_count;      fargs[1].uint32_arg = e->max_per_bucket;
        fargs[2].ptr_arg = e->d_arr_a;              fargs[3].ptr_arg = e->d_arr_b;
        fargs[4].ptr_arg = e->d_candidate_keys;     fargs[5].ptr_arg = e->d_candidate_cnt;
        fargs[6].ptr_arg = e->d_candidate_overflow; fargs[7].uint32_arg = OCL_COLL_CANDIDATE_CAP;
        fargs[8].uint32_arg = key_bits;             fargs[9].uint32_arg = max_tsize_words;
        fargs[10].ptr_arg = data->collect_stats ? e->d_filter_iters_hist : NULL;
        fargs[11].uint32_arg = (uint32_t)hash_bytes;  /* GPU_ARG_LOCAL: size, not value */
        gpu_launch_set(&e->th.launch[e->k_filter_per_bucket], fargs);

        err = clEnqueueNDRangeKernel(q, e->th.launch[e->k_filter_per_bucket].kernel_func,
                1, NULL, &global, &local, 0, NULL, NULL);
        if (err != CL_SUCCESS) {
            fprintf(stderr, "ocl_collision_run: filter_per_bucket failed: %s\n", clGetErrorString(err));
            return OCL_COLLISION_CL_ERROR;
        }
    }

    clEnqueueReadBuffer(q, e->d_candidate_cnt, CL_FALSE, 0, sizeof(uint32_t), &cand_cnt, 0, NULL, NULL);
    clEnqueueReadBuffer(q, e->d_candidate_overflow, CL_FALSE, 0, sizeof(uint32_t), &h_candidate_overflow, 0, NULL, NULL);
    if (data->collect_stats)
        clEnqueueReadBuffer(q, e->d_filter_iters_hist, CL_FALSE, 0, sizeof(data->filter_iters_hist),
                data->filter_iters_hist, 0, NULL, NULL);
    clFinish(q);
    data->candidate_count = cand_cnt;

    if (h_candidate_overflow || cand_cnt > OCL_COLL_CANDIDATE_CAP)
        return OCL_COLLISION_CANDIDATE_OVERFLOW;
    if (cand_cnt == 0)
        return OCL_COLLISION_OK;

    /* --- 3.3c: sort, dedup, secondary hash --- */
    if (ocl_radix_sort_u64(e->prim, e->d_candidate_keys, 0, e->d_sorted_keys, 0, cand_cnt) != 0)
        return OCL_COLLISION_CL_ERROR;

    ocl_fill_u32(e->prim, e->d_dedup_cnt, 0, 1, 0u);
    {
        gpu_arg_t dargs[4];
        size_t global = round_up_to(cand_cnt, 256u);
        size_t local = 256;
        cl_int err;
        dargs[0].ptr_arg = e->d_sorted_keys; dargs[1].uint32_arg = cand_cnt;
        dargs[2].ptr_arg = e->d_dedup_keys;  dargs[3].ptr_arg = e->d_dedup_cnt;
        gpu_launch_set(&e->th.launch[e->k_dedup], dargs);
        err = clEnqueueNDRangeKernel(q, e->th.launch[e->k_dedup].kernel_func,
                1, NULL, &global, &local, 0, NULL, NULL);
        if (err != CL_SUCCESS) { fprintf(stderr, "dedup failed: %s\n", clGetErrorString(err)); return OCL_COLLISION_CL_ERROR; }
    }
    clEnqueueReadBuffer(q, e->d_dedup_cnt, CL_TRUE, 0, sizeof(uint32_t), &dedup_cnt, 0, NULL, NULL);
    data->dedup_count = dedup_cnt;
    if (dedup_cnt == 0)
        return OCL_COLLISION_OK;

    c_ilog2 = host_ilog2(dedup_cnt + 1u);
    if (c_ilog2 > OCL_COLL_MAX_C_ILOG2) c_ilog2 = OCL_COLL_MAX_C_ILOG2;
    hash_value = (1u << c_ilog2) - 1u;
    dsize = (1u << c_ilog2) + 1u;
    ssize = (1u << c_ilog2) / 32u + 1u;

    ocl_fill_u32(e->prim, e->d_D, 0, dsize, 0u);
    ocl_fill_u32(e->prim, e->d_S, 0, ssize, 0u);
    {
        gpu_arg_t cargs[5];
        size_t global = round_up_to(dedup_cnt, 256u);
        size_t local = 256;
        cl_int err;
        cargs[0].ptr_arg = e->d_dedup_keys; cargs[1].uint32_arg = dedup_cnt;
        cargs[2].uint32_arg = hash_value;   cargs[3].ptr_arg = e->d_D; cargs[4].ptr_arg = e->d_S;
        gpu_launch_set(&e->th.launch[e->k_count_secondary], cargs);
        err = clEnqueueNDRangeKernel(q, e->th.launch[e->k_count_secondary].kernel_func,
                1, NULL, &global, &local, 0, NULL, NULL);
        if (err != CL_SUCCESS) { fprintf(stderr, "count_secondary failed: %s\n", clGetErrorString(err)); return OCL_COLLISION_CL_ERROR; }
    }

    if (ocl_scan_exclusive_u32(e->prim, e->d_D, 0, dsize) != 0)
        return OCL_COLLISION_CL_ERROR;
    clEnqueueCopyBuffer(q, e->d_D, e->d_D_pos, 0, 0, (size_t)dsize * sizeof(uint32_t), 0, NULL, NULL);
    {
        gpu_arg_t sargs[5];
        size_t global = round_up_to(dedup_cnt, 256u);
        size_t local = 256;
        cl_int err;
        sargs[0].ptr_arg = e->d_dedup_keys; sargs[1].uint32_arg = dedup_cnt;
        sargs[2].uint32_arg = hash_value;   sargs[3].ptr_arg = e->d_D_pos; sargs[4].ptr_arg = e->d_X;
        gpu_launch_set(&e->th.launch[e->k_scatter_secondary], sargs);
        err = clEnqueueNDRangeKernel(q, e->th.launch[e->k_scatter_secondary].kernel_func,
                1, NULL, &global, &local, 0, NULL, NULL);
        if (err != CL_SUCCESS) { fprintf(stderr, "scatter_secondary failed: %s\n", clGetErrorString(err)); return OCL_COLLISION_CL_ERROR; }
    }


    /* --- 3.3d: match + emit --- */
    ocl_fill_u32(e->prim, e->d_value_counts, 0, dedup_cnt + 1u, 0u);
    ocl_fill_u32(e->prim, e->d_value_match_total, 0, 1, 0u);
    ocl_fill_u32(e->prim, e->d_value_overflow, 0, 1, 0u);

    use_match_arena = (dedup_cnt <= OCL_COLL_VALUE_MATCH_CAP / OCL_COLL_MATCH_ARENA_WIDTH) ? 1u : 0u;
    data->match_arena_attempt_count = use_match_arena;
    data->match_arena_capacity_skip_count = use_match_arena ? 0u : 1u;
    {
        size_t global = round_up_to(n, 256u);
        size_t local = 256;
        cl_int err;
        if (use_match_arena) {
            gpu_arg_t margs[16];
            margs[0].ptr_arg = data->keys_in;       margs[1].uint64_arg = data->keys_in_off;
            margs[2].uint32_arg = data->root_bytes; margs[3].ptr_arg = data->data_in;
            margs[4].uint64_arg = data->data_in_off;margs[5].uint32_arg = n;
            margs[6].uint32_arg = hash_value;       margs[7].ptr_arg = e->d_S;
            margs[8].ptr_arg = e->d_D;              margs[9].ptr_arg = e->d_X;
            margs[10].uint32_arg = dedup_cnt;       margs[11].ptr_arg = e->d_value_counts;
            margs[12].ptr_arg = e->d_value_match_total; margs[13].ptr_arg = e->d_value_overflow;
            margs[14].ptr_arg = e->d_matched_values;    margs[15].uint32_arg = OCL_COLL_VALUE_MATCH_CAP;
            gpu_launch_set(&e->th.launch[e->k_count_and_store_matched_values], margs);
            err = clEnqueueNDRangeKernel(q, e->th.launch[e->k_count_and_store_matched_values].kernel_func,
                    1, NULL, &global, &local, 0, NULL, NULL);
        } else {
            gpu_arg_t margs[15];
            margs[0].ptr_arg = data->keys_in;       margs[1].uint64_arg = data->keys_in_off;
            margs[2].uint32_arg = data->root_bytes; margs[3].ptr_arg = data->data_in;
            margs[4].uint64_arg = data->data_in_off;margs[5].uint32_arg = n;
            margs[6].uint32_arg = hash_value;       margs[7].ptr_arg = e->d_S;
            margs[8].ptr_arg = e->d_D;              margs[9].ptr_arg = e->d_X;
            margs[10].uint32_arg = dedup_cnt;       margs[11].ptr_arg = e->d_value_counts;
            margs[12].ptr_arg = e->d_value_match_total; margs[13].ptr_arg = e->d_value_overflow;
            margs[14].uint32_arg = OCL_COLL_VALUE_MATCH_CAP;
            gpu_launch_set(&e->th.launch[e->k_count_matched_values], margs);
            err = clEnqueueNDRangeKernel(q, e->th.launch[e->k_count_matched_values].kernel_func,
                    1, NULL, &global, &local, 0, NULL, NULL);
        }
        if (err != CL_SUCCESS) {
            fprintf(stderr, "ocl_collision_run: match kernel failed: %s\n", clGetErrorString(err));
            return OCL_COLLISION_CL_ERROR;
        }
    }
    clEnqueueReadBuffer(q, e->d_value_match_total, CL_TRUE, 0, sizeof(uint32_t), &h_value_match_total, 0, NULL, NULL);
    clEnqueueReadBuffer(q, e->d_value_overflow, CL_TRUE, 0, sizeof(uint32_t), &h_value_overflow, 0, NULL, NULL);
    data->value_match_count = h_value_match_total;

    if ((!use_match_arena && h_value_overflow) || h_value_match_total > OCL_COLL_VALUE_MATCH_CAP)
        return OCL_COLLISION_VALUE_MATCH_OVERFLOW;
    if (h_value_match_total == 0)
        return OCL_COLLISION_OK;

    if (use_match_arena && !h_value_overflow) {
        gpu_arg_t eargs[10];
        size_t global = round_up_to(dedup_cnt, 256u);
        size_t local = 256;
        cl_int err;
        eargs[0].ptr_arg = e->d_X;             eargs[1].ptr_arg = e->d_matched_values;
        eargs[2].ptr_arg = e->d_value_counts;  eargs[3].uint32_arg = dedup_cnt;
        eargs[4].uint32_arg = data->shift;     eargs[5].uint32_arg = data->root_bytes;
        eargs[6].ptr_arg = data->q_batch;      eargs[7].uint64_arg = data->q_batch_off;
        eargs[8].ptr_arg = data->found_array;  eargs[9].uint64_arg = data->found_array_off;
        gpu_launch_set(&e->th.launch[e->k_emit_found_arena], eargs);
        err = clEnqueueNDRangeKernel(q, e->th.launch[e->k_emit_found_arena].kernel_func,
                1, NULL, &global, &local, 0, NULL, NULL);
        if (err != CL_SUCCESS) { fprintf(stderr, "emit_found_arena failed: %s\n", clGetErrorString(err)); return OCL_COLLISION_CL_ERROR; }
        clFinish(q);
        return OCL_COLLISION_OK;
    }

    if (data->debug && use_match_arena)
        fprintf(stderr, "ocl_collision: match arena overflow (dedup %u matched %u width %u), fallback\n",
                dedup_cnt, h_value_match_total, (unsigned)OCL_COLL_MATCH_ARENA_WIDTH);
    if (use_match_arena)
        data->match_arena_fallback_count = 1u;

    if (ocl_scan_exclusive_u32(e->prim, e->d_value_counts, 0, dedup_cnt + 1u) != 0)
        return OCL_COLLISION_CL_ERROR;
    clEnqueueCopyBuffer(q, e->d_value_counts, e->d_value_cursor, 0, 0,
            (size_t)dedup_cnt * sizeof(uint32_t), 0, NULL, NULL);
    {
        gpu_arg_t sargs[14];
        size_t global = round_up_to(n, 256u);
        size_t local = 256;
        cl_int err;
        sargs[0].ptr_arg = data->keys_in;       sargs[1].uint64_arg = data->keys_in_off;
        sargs[2].uint32_arg = data->root_bytes; sargs[3].ptr_arg = data->data_in;
        sargs[4].uint64_arg = data->data_in_off;sargs[5].uint32_arg = n;
        sargs[6].uint32_arg = hash_value;       sargs[7].ptr_arg = e->d_S;
        sargs[8].ptr_arg = e->d_D;              sargs[9].ptr_arg = e->d_X;
        sargs[10].uint32_arg = dedup_cnt;       sargs[11].ptr_arg = e->d_value_cursor;
        sargs[12].ptr_arg = e->d_matched_values;sargs[13].uint32_arg = OCL_COLL_VALUE_MATCH_CAP;
        gpu_launch_set(&e->th.launch[e->k_scatter_matched_values], sargs);
        err = clEnqueueNDRangeKernel(q, e->th.launch[e->k_scatter_matched_values].kernel_func,
                1, NULL, &global, &local, 0, NULL, NULL);
        if (err != CL_SUCCESS) { fprintf(stderr, "scatter_matched_values failed: %s\n", clGetErrorString(err)); return OCL_COLLISION_CL_ERROR; }
    }

    {
        gpu_arg_t eargs[10];
        size_t global = round_up_to(dedup_cnt, 256u);
        size_t local = 256;
        cl_int err;
        eargs[0].ptr_arg = e->d_X;              eargs[1].ptr_arg = e->d_matched_values;
        eargs[2].ptr_arg = e->d_value_counts;   eargs[3].uint32_arg = dedup_cnt;
        eargs[4].uint32_arg = data->shift;      eargs[5].uint32_arg = data->root_bytes;
        eargs[6].ptr_arg = data->q_batch;       eargs[7].uint64_arg = data->q_batch_off;
        eargs[8].ptr_arg = data->found_array;   eargs[9].uint64_arg = data->found_array_off;
        gpu_launch_set(&e->th.launch[e->k_emit_found], eargs);
        err = clEnqueueNDRangeKernel(q, e->th.launch[e->k_emit_found].kernel_func,
                1, NULL, &global, &local, 0, NULL, NULL);
        if (err != CL_SUCCESS) { fprintf(stderr, "emit_found failed: %s\n", clGetErrorString(err)); return OCL_COLLISION_CL_ERROR; }
    }
    clFinish(q);
    return OCL_COLLISION_OK;
}
