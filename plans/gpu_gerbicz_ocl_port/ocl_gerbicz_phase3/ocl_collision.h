/*--------------------------------------------------------------------
 * ocl_collision.h -- OpenCL port of collision_engine.cu's host struct
 * and orchestration (collision_engine / collision_engine_run).
 *
 * Built on Phase 1 (ocl_gerbicz_ctx.h) and Phase 2 (ocl_primitives.h,
 * reused directly for reduce-max/scan/sort rather than reimplemented).
 *
 * =======================================================================
 * WHAT CARRIES A BASE-OFFSET ARGUMENT AND WHAT DOESN'T
 * =======================================================================
 * Per decision #4 rule 5 ("a buffer used only in full never takes an
 * offset arg"): this engine's own scratch buffers (arr_a/arr_b,
 * candidate_keys, D/S/X, value_counts, matched_values, ...) are
 * allocated once per engine instance and used in full on every call --
 * they never carry an offset. The four buffers that mirror CUDA's
 * collision_data_t external device pointers -- keys_in, data_in
 * (values), q_batch, found_array -- are genuinely per-batch, externally
 * owned allocations that Phase 5's real driver will slice out of a
 * larger buffer, so all four carry a cl_ulong element offset, exactly
 * like collision_data_t's own CUdeviceptr fields do today (implicitly,
 * via whatever the CUDA driver's own pointer arithmetic already does
 * upstream of this engine).
 *
 * =======================================================================
 * GPU_ARG_LOCAL (new)
 * =======================================================================
 * filter_per_bucket_kernel's hash table is CUDA dynamic shared memory,
 * sized per-launch from a runtime-computed max_tsize_words. OpenCL's
 * equivalent is a dynamically-sized __local kernel PARAMETER, whose
 * size is set via clSetKernelArg(kernel, idx, size, NULL) -- but
 * ocl_xface.h's gpu_arg_type_t had no case for "set this arg's SIZE,
 * not its VALUE". Added GPU_ARG_LOCAL (see ocl_xface_arg_local.h.patch)
 * for exactly this; gpu_arg_t.uint32_arg is repurposed to carry the
 * byte size for that one argument slot.
 *
 * =======================================================================
 * OVERFLOW HANDLING (task 3.4)
 * =======================================================================
 * CUDA calls exit(-1) on candidate overflow and value-match overflow.
 * Per decision #7 and the plan's own open issue, this port returns an
 * ocl_collision_status_t instead. The registry already has an unused
 * STAGE1_OVERFLOW_SKIP (stage1_engine.h) for exactly this shape of
 * failure; Phase 5's driver is where that actually gets wired to a
 * "skip this cell" decision -- this phase only defines the contract
 * (the two overflow status codes below) so Phase 5 doesn't have to
 * invent one.
 * --------------------------------------------------------------------- */

#ifndef _OCL_COLLISION_H_
#define _OCL_COLLISION_H_

#include "ocl_gerbicz_ctx.h"
#include "ocl_primitives.h"

#ifdef __cplusplus
extern "C" {
#endif

/* collision_engine.cu constants (30-37), unchanged. */
#define OCL_COLL_BLOCK_THREADS     128u
#define OCL_COLL_MAX_FILTER_ITERS  20
#define OCL_COLL_CANDIDATE_CAP     (1u << 22)
#define OCL_COLL_VALUE_MATCH_CAP   OCL_COLL_CANDIDATE_CAP
#define OCL_COLL_MATCH_ARENA_WIDTH 8u
#define OCL_COLL_MAX_C_ILOG2       20u
#define OCL_COLL_MAX_DSIZE         ((1u << OCL_COLL_MAX_C_ILOG2) + 64u)
#define OCL_COLL_MAX_SSIZE         (((1u << OCL_COLL_MAX_C_ILOG2) >> 5) + 64u)
/* stage1_core.h */
#define OCL_COLL_FOUND_ARRAY_SIZE  1000u
#define OCL_COLL_MIN_KEY_BITS      20  /* COLLISION_ENGINE_MIN_KEY_BITS */

typedef enum {
    OCL_COLLISION_OK = 0,
    OCL_COLLISION_INVALID_INPUT,
    OCL_COLLISION_CANDIDATE_OVERFLOW,     /* -> STAGE1_OVERFLOW_SKIP, Phase 5 */
    OCL_COLLISION_VALUE_MATCH_OVERFLOW,   /* -> STAGE1_OVERFLOW_SKIP, Phase 5 */
    OCL_COLLISION_CL_ERROR                /* an OpenCL call itself failed */
} ocl_collision_status_t;

/* Mirrors collision_data_t (collision_engine.h), with cl_mem+offset
 * pairs replacing CUdeviceptr for the four external buffers. */
typedef struct {
    cl_mem   keys_in;      cl_ulong keys_in_off;      /* root_bytes-element units */
    cl_mem   data_in;      cl_ulong data_in_off;      /* uint32 element units */
    cl_mem   q_batch;      cl_ulong q_batch_off;      /* specialq_t element units */
    cl_mem   found_array;  cl_ulong found_array_off;  /* found_t element units,
                                                          OCL_COLL_FOUND_ARRAY_SIZE elements */

    uint32_t num_elements;
    uint32_t key_bits;
    uint32_t root_bytes;   /* 4 or 8 */
    uint32_t shift;
    int      bucket_hash;  /* 0 or 1 */
    int      debug;
    int      collect_stats;

    /* outputs -- filled in by ocl_collision_run() */
    uint32_t bucket_max;
    uint32_t candidate_count;
    uint32_t dedup_count;
    uint32_t value_match_count;
    uint32_t bucket_grow_count;
    uint32_t hash_cap_count;
    uint32_t match_arena_attempt_count;
    uint32_t match_arena_fallback_count;
    uint32_t match_arena_capacity_skip_count;
    uint32_t filter_iters_hist[102];
} ocl_collision_data_t;

typedef struct {
    ocl_gerbicz_device_t *gd;    /* not owned */
    ocl_primitives_t     *prim;  /* not owned -- reused for reduce-max/scan/sort */
    cl_program            program;
    ocl_thread_t          th;

    int k_scatter_roots;
    int k_filter_per_bucket;
    int k_dedup;
    int k_count_secondary;
    int k_scatter_secondary;
    int k_count_matched_values;
    int k_scatter_matched_values;
    int k_count_and_store_matched_values;
    int k_emit_found;
    int k_emit_found_arena;

    /* Sizing state, mirroring collision_engine's own struct. */
    uint32_t max_n;
    uint32_t max_key_bits;
    uint32_t max_per_bucket;

    /* Per-max_per_bucket buffers (grown via ensure_capacity). */
    cl_mem d_arr_a, d_arr_b;

    /* Fixed-size buffers, allocated once at init (mirrors CUDA: these
     * are sized off compile-time constants, not off n/key_bits).
     *
     * Simplification vs. the CUDA original: CUDA keeps d_D separate
     * from d_D_scan (and d_value_counts separate from d_value_offsets)
     * because cub::DeviceScan::ExclusiveSum always writes to a distinct
     * output buffer. Phase 2's ocl_scan_exclusive_u32 scans IN PLACE,
     * and in both cases the CUDA original never reads the pre-scan
     * buffer again once the scanned version exists -- so this port
     * scans d_D and d_value_counts in place and reuses the same buffer
     * for what CUDA calls d_D_scan / d_value_offsets, one buffer fewer
     * each. d_D_pos / d_value_cursor remain separate, mutable COPIES
     * of the scanned array (made via clEnqueueCopyBuffer), exactly as
     * CUDA's own d_D_pos / d_value_cursor are copies of d_D_scan /
     * d_value_offsets -- those two really do need to diverge from the
     * static scanned array during the scatter step. */
    cl_mem d_bucket_count;         /* uint, NUM_BUCKETS */
    cl_mem d_bucket_overflow;      /* uint, 1 */
    cl_mem d_candidate_keys;       /* ulong, CANDIDATE_CAP */
    cl_mem d_sorted_keys;          /* ulong, CANDIDATE_CAP */
    cl_mem d_dedup_keys;           /* ulong, CANDIDATE_CAP */
    cl_mem d_candidate_cnt;        /* uint, 1 */
    cl_mem d_candidate_overflow;   /* uint, 1 */
    cl_mem d_dedup_cnt;            /* uint, 1 */
    cl_mem d_D, d_D_pos;           /* uint, MAX_DSIZE -- d_D doubles as CUDA's d_D_scan once scanned in place */
    cl_mem d_S;                    /* uint, MAX_SSIZE */
    cl_mem d_X;                    /* ulong, CANDIDATE_CAP */
    cl_mem d_value_counts;         /* uint, VALUE_MATCH_CAP+1 -- doubles as CUDA's d_value_offsets once scanned in place */
    cl_mem d_value_cursor;         /* uint, VALUE_MATCH_CAP+1 */
    cl_mem d_matched_values;       /* uint, VALUE_MATCH_CAP */
    cl_mem d_value_match_total;    /* uint, 1 */
    cl_mem d_value_overflow;       /* uint, 1 */
    cl_mem d_filter_iters_hist;    /* uint, 102 */
} ocl_collision_engine_t;

int  ocl_collision_init(ocl_collision_engine_t *e, ocl_gerbicz_device_t *gd,
                         ocl_primitives_t *prim,
                         const char *cl_source_path, const char *cache_dir);
void ocl_collision_free(ocl_collision_engine_t *e);

/* Runs the full pipeline (scatter -> filter -> sort/dedup -> secondary
 * hash -> match -> emit) for one batch, matching collision_engine_run()
 * exactly except for the two overflow cases (see ocl_collision_status_t
 * above) and the no-sub-group compaction (decision #6). */
ocl_collision_status_t ocl_collision_run(ocl_collision_engine_t *e,
                                          ocl_collision_data_t *data);

#ifdef __cplusplus
}
#endif

#endif /* _OCL_COLLISION_H_ */
