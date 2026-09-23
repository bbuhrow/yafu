/* ocl_collision_kernels.cl -- OpenCL port of collision_engine.cu's
 * kernels (scatter_roots_kernel, filter_per_bucket_kernel, dedup_kernel,
 * count/scatter_secondary_kernel, the three match kernels,
 * emit_found_kernel/emit_found_arena_kernel).
 *
 * Locked decision #6 (no-sub-group baseline): every place the CUDA
 * source uses __ballot_sync/__shfl_sync/__match_any_sync for warp-scoped
 * stream compaction (scatter_roots_kernel's SCATTER_DIRECT_ATOMIC #else
 * branch, filter_per_bucket_kernel's survivor packing) is replaced here
 * with a single global or local atomicAdd per surviving thread -- the
 * CUDA source's own SCATTER_DIRECT_ATOMIC branch is explicitly named in
 * the plan as "the model" for this. Per the plan's own Fact ("filter
 * output is order-independent"), this changes nothing about
 * correctness, only which thread gets which output slot.
 *
 * -D flags baked in by the host (see ocl_collision.c):
 *   LOG2_NUM_BUCKETS, BUCKET_HASH_MIX  -- from collision_bucket.h, via
 *     ocl_gerbicz_build_opts() (Phase 1) so this translation unit can
 *     never drift from collision_bucket.h's own values.
 *   COLL_BLOCK_THREADS  -- matches CUDA's BLOCK_THREADS (128); kept as
 *     a named constant (not silently reusing Phase 2's PRIM_LOCAL_SIZE)
 *     because filter_per_bucket_kernel's local arrays are sized off it.
 *   COLL_MAX_FILTER_ITERS, COLL_MATCH_ARENA_WIDTH -- match
 *     collision_engine.cu's MAX_FILTER_ITERS / MATCH_ARENA_WIDTH.
 *
 * Base-offset convention: only buffers that are genuinely EXTERNAL,
 * shared, per-batch allocations in the real driver (Phase 5's
 * data->keys_in / data->data_in / data->q_batch / data->found_array,
 * matching collision_data_t) carry an offset argument, per decision #4
 * rule 5 ("a buffer used only in full never takes an offset arg") --
 * the engine's own internal scratch buffers (arr_a, D, S, X,
 * candidate_keys, ...) are never sliced within this phase, so they
 * don't carry one. See ocl_collision.h for the exact rationale.
 */

#ifndef HAVE_SUBGROUPS
#define HAVE_SUBGROUPS 0
#endif
#ifndef FOUND_ARRAY_SIZE
#define FOUND_ARRAY_SIZE 1000u   /* stage1_core.h; also passed as -D, this is just a safety net */
#endif

#define NUM_BUCKETS (1u << LOG2_NUM_BUCKETS)
#define BUCKET_MASK (NUM_BUCKETS - 1u)

/* ---- collision_bucket.h::compute_bucket(), hand-ported (OpenCL C has
 * no <stdint.h>/__host__ __device__; see file header). Must stay
 * byte-for-byte equivalent to the CUDA original. ---- */
inline uint
compute_bucket(ulong key, int hash_mode)
{
    if (hash_mode)
        return (uint)((key * (ulong)BUCKET_HASH_MIX) >> (64 - LOG2_NUM_BUCKETS));
    return (uint)(key & BUCKET_MASK);
}

/* ---- stage1_core.h mirrors (must stay byte-layout-identical to the
 * host's found_t/specialq_t -- both are plain C structs of fixed-width
 * scalar types with natural alignment, so field order alone is
 * sufficient; no packing pragmas needed on either side). ---- */
typedef struct {
    uint  p;
    uint  pad;
    ulong pp;
    ulong root;
} specialq_t;

typedef struct {
    uint  p1;
    uint  p2;
    uint  q;
    uint  pad;
    ulong qroot;
    long  offset;
} found_t;

/* ---- collision_engine.cu device helpers, hand-ported ---- */

inline uint
compute_ilog2(uint cnt, uint key_bits)
{
    uint lg, il, hash_bits_available;
    if (cnt == 0)
        return 5u;
    lg = 31u - clz(cnt);
    il = max(5u, lg + 5u);
    hash_bits_available = key_bits - LOG2_NUM_BUCKETS;
    if (il > hash_bits_available)
        il = hash_bits_available;
    return il;
}

inline uint
compute_capped_ilog2(uint cnt, uint key_bits, uint max_tsize_words)
{
    uint il = compute_ilog2(cnt, key_bits);
    uint max_il = 5u;
    if (max_tsize_words > 1u)
        max_il += 31u - clz(max_tsize_words);
    if (il > max_il)
        il = max_il;
    return il;
}

inline void
clear_table_parallel(__local uint *tbl, uint nwords, uint lid, uint lsize)
{
    uint j;
    for (j = lid; j < nwords; j += lsize)
        tbl[j] = 0u;
}

inline uint
gcd32_generic(uint a, uint b)
{
    while (b != 0u) {
        uint t = a % b;
        a = b;
        b = t;
    }
    return a;
}

inline void
store_hit_collision(__global found_t *found_array, uint found_array_size,
                     uint p1, uint p2, long root,
                     __global const specialq_t *q)
{
    uint index = atomic_add((volatile __global uint *)&found_array[0].p1, 1u);
    if (index < found_array_size - 1u) {
        __global found_t *f = found_array + index + 1u;
        f->p1 = p1;
        f->p2 = p2;
        f->q = q->p;
        f->qroot = q->root;
        f->offset = root;
    }
}

/* ===================================================================
 * scatter_roots_kernel
 *
 * `roots` arrives as a byte pointer (OpenCL C has no portable
 * kernel-parameter void*, unlike CUDA/C) plus root_bytes to reinterpret
 * it at either 4 or 8 bytes per element, exactly mirroring the CUDA
 * source's own `(const uint32*)roots` / `(const uint64*)roots` cast.
 * roots_off is in ELEMENTS of root_bytes (not raw bytes) to keep the
 * base-offset convention's "offset is in elements of the pointee type"
 * rule meaningful even though the pointee type is runtime-chosen.
 * =================================================================== */
__kernel void ocl_scatter_roots(
        __global const uchar *roots, ulong roots_off, uint root_bytes,
        uint n, uint max_per_bucket, int hash_mode,
        __global uint *bucket_count,
        __global uint *overflow_flag,
        __global ulong *bucket_storage)
{
    uint tid = (uint)get_global_id(0);
    ulong k;
    uint bucket, slot;

    if (tid >= n)
        return;

    roots += roots_off * root_bytes;
    if (root_bytes == 4u)
        k = ((__global const uint *)roots)[tid];
    else
        k = ((__global const ulong *)roots)[tid];
    if (k == 0)
        return;

    bucket = compute_bucket(k, hash_mode);
    /* No-sub-group baseline (decision #6): direct global atomic, the
     * CUDA source's own SCATTER_DIRECT_ATOMIC branch. */
    slot = atomic_add(bucket_count + bucket, 1u);
    if (slot >= max_per_bucket) {
        atomic_or(overflow_flag, 1u);
        return;
    }
    bucket_storage[(size_t)bucket * max_per_bucket + slot] = k;
}

/* ===================================================================
 * filter_per_bucket_kernel
 *
 * s_hash is a dynamically-sized __local argument (GPU_ARG_LOCAL, see
 * ocl_xface_arg_local.h.patch) of 3*max_tsize_words*4 bytes, the exact
 * OpenCL equivalent of CUDA's cudaFuncSetAttribute +
 * MaxDynamicSharedMemorySize dance -- partitioned into T/T2/T3 exactly
 * as the CUDA source partitions its own s_hash.
 *
 * iters_hist_out may be NULL (portable in OpenCL: a NULL __global
 * pointer argument, checked before use) when collect_stats==0, mirroring
 * the CUDA original's "zero atomic ops when off" behavior exactly --
 * host passes a null cl_mem via gpu_arg_t.ptr_arg = NULL in that case.
 * =================================================================== */
__kernel __attribute__((reqd_work_group_size(COLL_BLOCK_THREADS, 1, 1)))
void ocl_filter_per_bucket(
        __global const uint *bucket_count_in,
        uint max_per_bucket,
        __global ulong *arr_a,
        __global ulong *arr_b,
        __global ulong *candidate_keys,
        __global uint *candidate_cnt,
        __global uint *candidate_overflow,
        uint candidate_cap,
        uint key_bits,
        uint max_tsize_words,
        __global uint *iters_hist_out,
        __local uint *s_hash)
{
    uint bucket = get_group_id(0);
    uint lid = get_local_id(0);
    uint cnt = bucket_count_in[bucket];
    uint cnt0;
    size_t offset;
    __global ulong *arr_in, *arr_out;
    __local uint *T, *T2, *T3;
    uint ilog2, tsize, hash_value2, my_shift2;
    int it;
    bool emit = false;

    __local uint s_nsize[COLL_MAX_FILTER_ITERS];
    __local uint s_cnt_out;

    if (cnt > max_per_bucket)
        cnt = max_per_bucket;
    if (cnt == 0)
        return;
    cnt0 = cnt;

    offset = (size_t)bucket * max_per_bucket;
    arr_in = arr_a + offset;
    arr_out = arr_b + offset;

    T = s_hash;
    T2 = s_hash + max_tsize_words;
    T3 = s_hash + 2u * max_tsize_words;

    ilog2 = compute_capped_ilog2(cnt, key_bits, max_tsize_words);
    tsize = 1u << (ilog2 - 5u);
    hash_value2 = (ilog2 >= 32u) ? 0xFFFFFFFFu : ((1u << ilog2) - 1u);
    my_shift2 = LOG2_NUM_BUCKETS;

    clear_table_parallel(T, tsize, lid, COLL_BLOCK_THREADS);
    clear_table_parallel(T2, tsize, lid, COLL_BLOCK_THREADS);
    barrier(CLK_LOCAL_MEM_FENCE);

    {
        uint j;
        for (j = lid; j < cnt; j += COLL_BLOCK_THREADS) {
            ulong item = arr_in[j];
            uint hv = (uint)((item >> my_shift2) & hash_value2);
            uint bit = 1u << (hv & 31u);
            uint prev = atomic_or(&T2[hv >> 5], bit);
            if (prev & bit)
                atomic_or(&T[hv >> 5], bit);
        }
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    for (it = 0; it < COLL_MAX_FILTER_ITERS; it++) {
        uint my_shift = my_shift2;
        uint hash_value = hash_value2;
        uint j_base;
        uint cnt_out;
        bool stop_zero, stop_cap, stop_conv, stop;
        __local uint *U, *U2, *U3;

        ilog2 = compute_capped_ilog2(cnt, key_bits, max_tsize_words);
        tsize = 1u << (ilog2 - 5u);

        if ((it & 1) == 0) { U = T; U2 = T2; U3 = T3; }
        else { U = T3; U2 = T2; U3 = T; }

        my_shift2 = my_shift + 6u;
        if (my_shift2 + ilog2 > key_bits)
            my_shift2 = LOG2_NUM_BUCKETS;
        hash_value2 = (ilog2 >= 32u) ? 0xFFFFFFFFu : ((1u << ilog2) - 1u);

        clear_table_parallel(U2, tsize, lid, COLL_BLOCK_THREADS);
        clear_table_parallel(U3, tsize, lid, COLL_BLOCK_THREADS);
        if (lid == 0)
            s_cnt_out = 0u;
        barrier(CLK_LOCAL_MEM_FENCE);

        for (j_base = 0; j_base < cnt; j_base += COLL_BLOCK_THREADS) {
            uint j = j_base + lid;
            bool survives = false;
            ulong item = 0;

            if (j < cnt) {
                item = arr_in[j];
                {
                    uint hv = (uint)((item >> my_shift) & hash_value);
                    uint bit = 1u << (hv & 31u);
                    if (U[hv >> 5] & bit) {
                        uint hv2 = (uint)((item >> my_shift2) & hash_value2);
                        uint bit2 = 1u << (hv2 & 31u);
                        uint prev = atomic_or(&U2[hv2 >> 5], bit2);
                        if (prev & bit2)
                            atomic_or(&U3[hv2 >> 5], bit2);
                        survives = true;
                    }
                }
            }

            /* No-sub-group baseline (decision #6): one atomic per
             * surviving thread instead of CUDA's warp-ballot+shuffle
             * compaction. Output order differs from CUDA's; per the
             * plan's Fact, filter output is order-independent, so this
             * cannot change the result set. */
            if (survives) {
                uint slot = atomic_add(&s_cnt_out, 1u);
                arr_out[slot] = item;
            }
        }

        barrier(CLK_LOCAL_MEM_FENCE);
        cnt_out = s_cnt_out;
        if (lid == 0)
            s_nsize[it] = cnt_out;
        barrier(CLK_LOCAL_MEM_FENCE);

        cnt = cnt_out;
        { __global ulong *tmp = arr_in; arr_in = arr_out; arr_out = tmp; }

        stop_zero = (cnt == 0u);
        stop_cap = !stop_zero && (it == COLL_MAX_FILTER_ITERS - 1);
        stop_conv = !stop_zero && !stop_cap && (it >= 3 &&
                s_nsize[it - 3] == s_nsize[it]);
        stop = stop_zero || stop_cap || stop_conv;
        if (stop) {
            if (iters_hist_out && lid == 0) {
                uint idx = (uint)it;
                uint size_bin, cat_base;
                if (idx > 20u) idx = 20u;
                atomic_add(iters_hist_out + idx, 1u);
                if (stop_zero)
                    atomic_add(iters_hist_out + 21u + idx, 1u);
                else if (stop_cap)
                    atomic_add(iters_hist_out + 42u + idx, 1u);

                size_bin = (cnt0 == 0u) ? 0u : (31u - clz(cnt0));
                if (size_bin > 12u) size_bin = 12u;
                cat_base = (it <= 3) ? 63u : (it == 4) ? 76u : 89u;
                atomic_add(iters_hist_out + cat_base + size_bin, 1u);
            }
            emit = true;
            break;
        }
    }

    if (!emit || cnt == 0)
        return;

    if (lid == 0)
        s_cnt_out = atomic_add(candidate_cnt, cnt);
    barrier(CLK_LOCAL_MEM_FENCE);

    {
        uint base = s_cnt_out;
        uint j;
        if (base + cnt > candidate_cap) {
            if (lid == 0)
                atomic_or(candidate_overflow, 1u);
            return;
        }
        for (j = lid; j < cnt; j += COLL_BLOCK_THREADS)
            candidate_keys[base + j] = arr_in[j];
    }
}

/* ===================================================================
 * dedup / secondary hash
 * =================================================================== */
__kernel void ocl_dedup(
        __global const ulong *sorted_keys, uint csize,
        __global ulong *dedup_keys, __global uint *dedup_cnt)
{
    uint tid = (uint)get_global_id(0);
    ulong k;
    bool is_first, has_neighbor;

    if (tid >= csize)
        return;

    k = sorted_keys[tid];
    is_first = (tid == 0) || (sorted_keys[tid - 1] != k);
    has_neighbor = (tid + 1u < csize) && (sorted_keys[tid + 1u] == k);
    if (is_first && has_neighbor) {
        uint pos = atomic_add(dedup_cnt, 1u);
        dedup_keys[pos] = k;
    }
}

__kernel void ocl_count_secondary(
        __global const ulong *dedup_keys, uint csize, uint hash_value,
        __global uint *D, __global uint *S)
{
    uint tid = (uint)get_global_id(0);
    ulong k; uint hv;
    if (tid >= csize)
        return;
    k = dedup_keys[tid];
    hv = (uint)(k & hash_value);
    atomic_add(D + hv, 1u);
    atomic_or(S + (hv >> 5), 1u << (hv & 31u));
}

__kernel void ocl_scatter_secondary(
        __global const ulong *dedup_keys, uint csize, uint hash_value,
        __global uint *D_pos, __global ulong *X)
{
    uint tid = (uint)get_global_id(0);
    ulong k; uint hv, slot;
    if (tid >= csize)
        return;
    k = dedup_keys[tid];
    hv = (uint)(k & hash_value);
    slot = atomic_add(D_pos + hv, 1u);
    X[slot] = k;
}

inline int
find_candidate_slot(ulong k, uint hash_value,
                     __global const uint *S, __global const uint *D,
                     __global const ulong *X, uint *slot_out)
{
    uint hv = (uint)(k & hash_value);
    uint lo, hi, j;

    if ((S[hv >> 5] & (1u << (hv & 31u))) == 0u)
        return 0;

    lo = D[hv];
    hi = D[hv + 1u];
    for (j = lo; j < hi; j++) {
        if (X[j] == k) {
            *slot_out = j;
            return 1;
        }
    }
    return 0;
}

/* ===================================================================
 * Match kernels. `roots`/`values` are external, per-batch buffers
 * (collision_data_t's keys_in/data_in) -- both carry the base-offset
 * convention.
 * =================================================================== */
__kernel void ocl_count_matched_values(
        __global const uchar *roots, ulong roots_off, uint root_bytes,
        __global const uint *values, ulong values_off,
        uint n, uint hash_value,
        __global const uint *S, __global const uint *D, __global const ulong *X,
        uint dedup_count,
        __global uint *value_counts, __global uint *match_total,
        __global uint *value_overflow, uint value_cap)
{
    uint tid = (uint)get_global_id(0);
    ulong k; uint slot;

    if (tid >= n)
        return;
    roots += roots_off * root_bytes;
    k = (root_bytes == 4u) ? ((__global const uint *)roots)[tid]
                            : ((__global const ulong *)roots)[tid];
    if (k == 0)
        return;

    if (find_candidate_slot(k, hash_value, S, D, X, &slot) && slot < dedup_count) {
        uint total_pos = atomic_add(match_total, 1u);
        atomic_add(value_counts + slot, 1u);
        if (total_pos >= value_cap)
            atomic_or(value_overflow, 1u);
    }
    (void)values; (void)values_off;
}

__kernel void ocl_scatter_matched_values(
        __global const uchar *roots, ulong roots_off, uint root_bytes,
        __global const uint *values, ulong values_off,
        uint n, uint hash_value,
        __global const uint *S, __global const uint *D, __global const ulong *X,
        uint dedup_count,
        __global uint *value_cursor, __global uint *matched_values, uint value_cap)
{
    uint tid = (uint)get_global_id(0);
    ulong k; uint slot;

    if (tid >= n)
        return;
    roots += roots_off * root_bytes;
    values += values_off;
    k = (root_bytes == 4u) ? ((__global const uint *)roots)[tid]
                            : ((__global const ulong *)roots)[tid];
    if (k == 0)
        return;

    if (find_candidate_slot(k, hash_value, S, D, X, &slot) && slot < dedup_count) {
        uint pos = atomic_add(value_cursor + slot, 1u);
        if (pos < value_cap)
            matched_values[pos] = values[tid];
    }
}

__kernel void ocl_count_and_store_matched_values(
        __global const uchar *roots, ulong roots_off, uint root_bytes,
        __global const uint *values, ulong values_off,
        uint n, uint hash_value,
        __global const uint *S, __global const uint *D, __global const ulong *X,
        uint dedup_count,
        __global uint *value_counts, __global uint *match_total,
        __global uint *value_overflow, __global uint *matched_values, uint value_cap)
{
    uint tid = (uint)get_global_id(0);
    ulong k; uint slot;

    if (tid >= n)
        return;
    roots += roots_off * root_bytes;
    values += values_off;
    k = (root_bytes == 4u) ? ((__global const uint *)roots)[tid]
                            : ((__global const ulong *)roots)[tid];
    if (k == 0)
        return;

    if (find_candidate_slot(k, hash_value, S, D, X, &slot) && slot < dedup_count) {
        uint total_pos = atomic_add(match_total, 1u);
        uint pos = atomic_add(value_counts + slot, 1u);
        size_t arena_pos = (size_t)slot * COLL_MATCH_ARENA_WIDTH + pos;

        if (total_pos >= value_cap) {
            atomic_or(value_overflow, 1u);
        } else if (pos < COLL_MATCH_ARENA_WIDTH && arena_pos < value_cap) {
            matched_values[arena_pos] = values[tid];
        } else {
            atomic_or(value_overflow, 1u);
        }
    }
}

/* ===================================================================
 * Emit kernels. q_batch and found_array are external, per-batch
 * buffers -- base-offset convention applies to both.
 * =================================================================== */
__kernel void ocl_emit_found(
        __global const ulong *X, __global const uint *matched_values,
        __global const uint *value_offsets, uint dedup_count,
        uint pshift, uint root_bytes,
        __global const specialq_t *q_batch, ulong q_batch_off,
        __global found_t *found_array, ulong found_array_off)
{
    uint slot = (uint)get_global_id(0);
    uint lo, hi, mask, i;
    long root;

    if (slot >= dedup_count)
        return;
    q_batch += q_batch_off;
    found_array += found_array_off;

    lo = value_offsets[slot];
    hi = value_offsets[slot + 1u];
    mask = (1u << pshift) - 1u;
    root = (root_bytes == 4u) ? (long)(int)((uint)X[slot]) : (long)X[slot];

    for (i = lo; i + 1u < hi; i++) {
        uint v1 = matched_values[i];
        uint q1 = v1 >> pshift;
        uint p1 = v1 & mask;
        uint j;
        for (j = i + 1u; j < hi; j++) {
            uint v2 = matched_values[j];
            uint q2 = v2 >> pshift;
            uint p2 = v2 & mask;
            if (q1 == q2 && gcd32_generic(p1, p2) == 1u) {
                store_hit_collision(found_array, FOUND_ARRAY_SIZE, p1, p2, root, q_batch + q1);
            }
        }
    }
}

__kernel void ocl_emit_found_arena(
        __global const ulong *X, __global const uint *matched_values,
        __global const uint *value_counts, uint dedup_count,
        uint pshift, uint root_bytes,
        __global const specialq_t *q_batch, ulong q_batch_off,
        __global found_t *found_array, ulong found_array_off)
{
    uint slot = (uint)get_global_id(0);
    uint cnt, base, mask, i;
    long root;

    if (slot >= dedup_count)
        return;
    q_batch += q_batch_off;
    found_array += found_array_off;

    cnt = value_counts[slot];
    if (cnt > COLL_MATCH_ARENA_WIDTH)
        cnt = COLL_MATCH_ARENA_WIDTH;
    base = slot * COLL_MATCH_ARENA_WIDTH;
    mask = (1u << pshift) - 1u;
    root = (root_bytes == 4u) ? (long)(int)((uint)X[slot]) : (long)X[slot];

    for (i = 0; i + 1u < cnt; i++) {
        uint v1 = matched_values[base + i];
        uint q1 = v1 >> pshift;
        uint p1 = v1 & mask;
        uint j;
        for (j = i + 1u; j < cnt; j++) {
            uint v2 = matched_values[base + j];
            uint q2 = v2 >> pshift;
            uint p2 = v2 & mask;
            if (q1 == q2 && gcd32_generic(p1, p2) == 1u) {
                store_hit_collision(found_array, FOUND_ARRAY_SIZE, p1, p2, root, q_batch + q1);
            }
        }
    }
}
