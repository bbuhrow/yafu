/* ocl_primitives_kernels.cl -- Phase 2 primitives (fill, reduce-max,
 * exclusive scan, LSD radix sort over ulong keys), replacing the CUB
 * calls in collision_engine.cu 741-757.
 *
 * All chunk-size constants come in as -D flags matching
 * ocl_primitives.h's OCL_PRIM_* macros (see ocl_primitives.c) so the
 * host and device sides can never drift apart -- do not hardcode 256,
 * 16, 4096, 8 or 256 (radix bins) a second time in this file.
 *
 * Every kernel follows the base-offset convention from ocl_gerbicz_ctx.h:
 * a ulong element-offset argument immediately follows the buffer it
 * offsets, added to the pointer as the kernel's first statement.
 *
 * Decision #7 ("bit-exact first, then optimize"): several kernels below
 * are deliberately not maximally parallel -- a single-thread-per-workgroup
 * scan of a small local array, a single-thread-per-workgroup radix
 * scatter. Each is marked "Phase 7 target" where a real optimization
 * exists. Correctness came first; none of this is accidental.
 */

#ifndef HAVE_SUBGROUPS
#define HAVE_SUBGROUPS 0
#endif

/* ===================================================================
 * Fill
 * =================================================================== */
__kernel void ocl_fill_u32(__global uint *buf, ulong buf_off, uint n, uint value)
{
    buf += buf_off;
    uint gid = (uint)get_global_id(0);
    if (gid < n)
        buf[gid] = value;
}

__kernel void ocl_fill_u64(__global ulong *buf, ulong buf_off, uint n, ulong value)
{
    buf += buf_off;
    uint gid = (uint)get_global_id(0);
    if (gid < n)
        buf[gid] = value;
}

/* ===================================================================
 * Reduce-max (uint32). Host loops: partials = reduce(in, n); if
 * count(partials) > 1, reduce again on the partials, until one value
 * remains. See ocl_reduce_max_u32() in ocl_primitives.c for the loop.
 * =================================================================== */
__kernel void ocl_reduce_max_u32(__global const uint *in, ulong in_off,
                                  __global uint *out, ulong out_off,
                                  uint n)
{
    in += in_off;
    out += out_off;

    __local uint lmax[PRIM_LOCAL_SIZE];
    uint lid = get_local_id(0);
    uint wg  = get_group_id(0);
    uint base = wg * PRIM_ELEMS_PER_WG;

    uint m = 0;
    for (uint i = lid; i < PRIM_ELEMS_PER_WG; i += PRIM_LOCAL_SIZE) {
        uint idx = base + i;
        if (idx < n)
            m = max(m, in[idx]);
    }
    lmax[lid] = m;
    barrier(CLK_LOCAL_MEM_FENCE);

    for (uint stride = PRIM_LOCAL_SIZE / 2; stride > 0; stride >>= 1) {
        if (lid < stride)
            lmax[lid] = max(lmax[lid], lmax[lid + stride]);
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if (lid == 0)
        out[wg] = lmax[0];
}

/* ===================================================================
 * Exclusive scan (uint32), two-level / three-kernel, no decoupled
 * look-back (decision #7).
 *
 * Level 1 (this kernel): each workgroup covers PRIM_ELEMS_PER_WG
 * elements. Each thread serially exclusive-scans its own contiguous
 * PRIM_ELEMS_PER_THREAD-element slice (registers, not shared), leaving
 * its total in thread_sum[lid]. Thread 0 then serially scans the
 * PRIM_LOCAL_SIZE thread_sum values -- Phase 7 target: replace with a
 * parallel Hillis-Steele scan of thread_sum[]; PRIM_LOCAL_SIZE (256)
 * serial adds is cheap enough relative to everything else that it
 * wasn't worth the extra complexity for a first correct version. Each
 * workgroup writes its block total to block_sums[wg].
 * =================================================================== */
__kernel void ocl_scan_local_u32(__global const uint *in, ulong in_off,
                                  __global uint *out, ulong out_off,
                                  __global uint *block_sums, ulong block_sums_off,
                                  uint n)
{
    in += in_off;
    out += out_off;
    block_sums += block_sums_off;

    __local uint thread_sum[PRIM_LOCAL_SIZE];
    __local uint thread_off[PRIM_LOCAL_SIZE];

    uint lid = get_local_id(0);
    uint wg  = get_group_id(0);
    uint base = wg * PRIM_ELEMS_PER_WG + lid * PRIM_ELEMS_PER_THREAD;

    uint local_vals[PRIM_ELEMS_PER_THREAD];
    uint running = 0;
    for (uint i = 0; i < PRIM_ELEMS_PER_THREAD; i++) {
        uint idx = base + i;
        uint v = (idx < n) ? in[idx] : 0u;
        local_vals[i] = running;
        running += v;
    }
    thread_sum[lid] = running;
    barrier(CLK_LOCAL_MEM_FENCE);

    if (lid == 0) {
        uint acc = 0;
        for (uint t = 0; t < PRIM_LOCAL_SIZE; t++) {
            uint v = thread_sum[t];
            thread_off[t] = acc;
            acc += v;
        }
        block_sums[wg] = acc;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    uint off = thread_off[lid];
    for (uint i = 0; i < PRIM_ELEMS_PER_THREAD; i++) {
        uint idx = base + i;
        if (idx < n)
            out[idx] = local_vals[i] + off;
    }
}

/* Level 2: exclusive-scan block_sums[0..n) IN PLACE using exactly one
 * workgroup. Host guarantees n <= PRIM_ELEMS_PER_WG (see
 * OCL_PRIM_SCAN_MAX_N in ocl_primitives.h) -- this is the two-level
 * scan's documented capacity limit, not a bug. Same thread-serial
 * structure as ocl_scan_local_u32, minus the "write to a separate out[]"
 * step. */
__kernel void ocl_scan_single_wg_u32(__global uint *buf, ulong buf_off, uint n)
{
    buf += buf_off;

    __local uint thread_sum[PRIM_LOCAL_SIZE];
    __local uint thread_off[PRIM_LOCAL_SIZE];

    uint lid = get_local_id(0);
    uint base = lid * PRIM_ELEMS_PER_THREAD;

    uint local_vals[PRIM_ELEMS_PER_THREAD];
    uint running = 0;
    for (uint i = 0; i < PRIM_ELEMS_PER_THREAD; i++) {
        uint idx = base + i;
        uint v = (idx < n) ? buf[idx] : 0u;
        local_vals[i] = running;
        running += v;
    }
    thread_sum[lid] = running;
    barrier(CLK_LOCAL_MEM_FENCE);

    if (lid == 0) {
        uint acc = 0;
        for (uint t = 0; t < PRIM_LOCAL_SIZE; t++) {
            uint v = thread_sum[t];
            thread_off[t] = acc;
            acc += v;
        }
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    uint off = thread_off[lid];
    for (uint i = 0; i < PRIM_ELEMS_PER_THREAD; i++) {
        uint idx = base + i;
        if (idx < n)
            buf[idx] = local_vals[i] + off;
    }
}

/* Level 3: add each (now-exclusive) block sum to every element of its
 * block. buf[] holds ocl_scan_local_u32's per-block-local exclusive
 * scan; after this kernel it's the full exclusive scan of the original
 * array. */
__kernel void ocl_scan_addback_u32(__global uint *buf, ulong buf_off,
                                    __global const uint *block_sums, ulong block_sums_off,
                                    uint n)
{
    buf += buf_off;
    block_sums += block_sums_off;

    uint gid = (uint)get_global_id(0);
    if (gid >= n)
        return;
    uint wg = gid / PRIM_ELEMS_PER_WG;
    buf[gid] += block_sums[wg];
}

/* ===================================================================
 * LSD radix sort over ulong keys, one pass = 8 bits (PRIM_RADIX_BITS),
 * 8 passes total for a 64-bit key. Two kernels per pass:
 *   1. histogram: per-workgroup per-bin counts, written TRANSPOSED as
 *      hist[bin * num_wg + wg] so that scanning the whole flattened
 *      array (host calls ocl_scan_exclusive_u32 on it) gives, for each
 *      (bin, workgroup) pair, exactly the global output offset a
 *      standard stable counting-sort digit pass needs -- bin-major
 *      order means "all of bin 0 across every workgroup, in workgroup
 *      order" comes before "all of bin 1", which is the ordering a
 *      digit pass must produce.
 *   2. scatter: re-reads its own chunk and writes each key to its
 *      global offset, advancing a per-bin cursor seeded from the
 *      scanned histogram.
 *
 * Stability (required for LSD correctness across passes) comes from:
 *   - workgroups process disjoint, ORDER-PRESERVING chunks of the
 *     input array (workgroup 0's chunk precedes workgroup 1's), and
 *   - scatter processes its chunk with a SINGLE thread, sequentially in
 *     original order (Phase 7 target: parallelize with a proper
 *     per-thread local prefix sum instead of one thread doing the
 *     whole chunk -- correctness first per decision #7).
 * =================================================================== */
__kernel void ocl_radix_histogram_u64(__global const ulong *keys, ulong keys_off,
                                       __global uint *hist, ulong hist_off,
                                       uint n, uint num_wg, uint shift)
{
    keys += keys_off;
    hist += hist_off;

    __local uint local_hist[PRIM_RADIX_BINS];
    uint lid = get_local_id(0);
    uint wg  = get_group_id(0);

    for (uint b = lid; b < PRIM_RADIX_BINS; b += PRIM_LOCAL_SIZE)
        local_hist[b] = 0;
    barrier(CLK_LOCAL_MEM_FENCE);

    uint base = wg * PRIM_ELEMS_PER_WG;
    for (uint i = lid; i < PRIM_ELEMS_PER_WG; i += PRIM_LOCAL_SIZE) {
        uint idx = base + i;
        if (idx < n) {
            uint digit = (uint)((keys[idx] >> shift) & (PRIM_RADIX_BINS - 1));
            atomic_inc(&local_hist[digit]);
        }
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    for (uint b = lid; b < PRIM_RADIX_BINS; b += PRIM_LOCAL_SIZE)
        hist[(size_t)b * num_wg + wg] = local_hist[b];
}

__kernel void ocl_radix_scatter_u64(__global const ulong *keys_in, ulong keys_in_off,
                                     __global ulong *keys_out, ulong keys_out_off,
                                     __global const uint *offsets, ulong offsets_off,
                                     uint n, uint num_wg, uint shift)
{
    keys_in += keys_in_off;
    keys_out += keys_out_off;
    offsets += offsets_off;

    /* __local, not a private per-work-item array: only lid==0 ever
     * touches it, so one workgroup-shared 256-entry buffer is enough --
     * declaring PRIM_RADIX_BINS uints as a PRIVATE array here blew the
     * per-work-item stack on PoCL's CPU backend (256 work-items x 1KB
     * private each, even though 255 of them never touch it). */
    __local uint cursor[PRIM_RADIX_BINS];

    uint lid = get_local_id(0);
    uint wg  = get_group_id(0);
    if (lid != 0)
        return; /* single-thread-per-workgroup scatter -- see header comment;
                    Phase 7 target, not a correctness shortcut being hidden */

    for (uint b = 0; b < PRIM_RADIX_BINS; b++)
        cursor[b] = offsets[(size_t)b * num_wg + wg];

    uint base = wg * PRIM_ELEMS_PER_WG;
    uint end = base + PRIM_ELEMS_PER_WG;
    if (end > n)
        end = n;
    for (uint idx = base; idx < end; idx++) {
        ulong k = keys_in[idx];
        uint digit = (uint)((k >> shift) & (PRIM_RADIX_BINS - 1));
        keys_out[cursor[digit]] = k;
        cursor[digit]++;
    }
}
