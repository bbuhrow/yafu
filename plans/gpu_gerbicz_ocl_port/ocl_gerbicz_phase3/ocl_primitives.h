/*--------------------------------------------------------------------
 * ocl_primitives.h  --  Phase 2: fill / reduce-max / exclusive-scan /
 * LSD radix sort, replacing the CUB calls in collision_engine.cu
 * 741-757 (cub::DeviceReduce::Max, cub::DeviceRadixSort::SortKeys,
 * cub::DeviceScan::ExclusiveSum x2).
 *
 * Built on Phase 1's ocl_shared.h / ocl_gerbicz_ctx.h: one program is
 * built (via ocl_build_program_cached) holding all of this phase's
 * kernels, and each primitive gets its own gpu_launch_t via
 * ocl_thread_init(). Every kernel here follows ocl_gerbicz_ctx.h's
 * base-offset convention (ulong element offsets immediately after the
 * buffer they offset) even though Phase 2 itself doesn't need slicing
 * yet -- Phase 3 will call these primitives on sub-ranges of its own
 * larger buffers, so the convention has to be there from day one.
 *
 * =======================================================================
 * SCOPE DECISION: keys are always ulong (64-bit), never native uint.
 * =======================================================================
 * collision_engine.cu's scatter_roots_kernel (111-153) reads a 4-byte
 * root as `uint32` and assigns it to a `uint64 k` -- an implicit ZERO
 * extension, not sign extension (sign-extension only happens later, on
 * emit, per the plan's Section 3 fact). d_candidate_keys/d_sorted_keys
 * are declared `uint64 *` unconditionally (collision_engine.cu 715-719),
 * regardless of root_bytes. So the "sort 32- or 64-bit keys" requirement
 * in the plan is about never truncating to `key_bits`, NOT about needing
 * two separate sort implementations for two different storage widths --
 * there is only ever one storage width (64-bit) at the point where
 * sorting happens. This phase therefore implements ONE radix sort, over
 * ulong keys, and does not special-case root_bytes==4. See STATUS for
 * the discrepancy writeup.
 *
 * Because the stored value is a zero-extended raw bit pattern (not a
 * sign-extended one), a plain unsigned LSD radix sort (no sign-bit-flip
 * trick) is correct here -- matches cub::DeviceRadixSort::SortKeys's own
 * default behavior on a plain uint64_t* buffer.
 * --------------------------------------------------------------------- */

#ifndef _OCL_PRIMITIVES_H_
#define _OCL_PRIMITIVES_H_

#include "ocl_gerbicz_ctx.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Chunk sizes baked into the .cl source at build time (see
 * ocl_primitives_kernels.cl); exposed here so host code can compute
 * buffer sizes (number of workgroups, block-sums array size, etc.)
 * without hardcoding them a second time. */
#define OCL_PRIM_LOCAL_SIZE        256u
#define OCL_PRIM_ELEMS_PER_THREAD  16u
#define OCL_PRIM_ELEMS_PER_WG      (OCL_PRIM_LOCAL_SIZE * OCL_PRIM_ELEMS_PER_THREAD) /* 4096 */
#define OCL_PRIM_RADIX_BITS        8u
#define OCL_PRIM_RADIX_BINS        256u   /* 1u << OCL_PRIM_RADIX_BITS */
#define OCL_PRIM_RADIX_PASSES_U64  8u     /* 64 / OCL_PRIM_RADIX_BITS */

/* Two-level scan capacity: the block-sums array (one entry per
 * workgroup) must itself fit in a single OCL_PRIM_ELEMS_PER_WG-sized
 * workgroup scan (see ocl_scan_single_wg_u32 in the .cl source) --
 * decision #7 "bit-exact first" says a plain two-pass scan is fine for
 * now, so this is a real, documented limit rather than a bug: */
#define OCL_PRIM_SCAN_MAX_N        (OCL_PRIM_ELEMS_PER_WG * OCL_PRIM_ELEMS_PER_WG) /* ~16.7M */

typedef struct {
    ocl_gerbicz_device_t *gd;         /* not owned */
    cl_program             program;
    ocl_thread_t           th;        /* one queue + all this phase's kernels */

    /* Kernel indices into th.launch[], set by ocl_primitives_init(). */
    int k_fill_u32;
    int k_fill_u64;
    int k_reduce_max_u32;
    int k_scan_local_u32;
    int k_scan_single_wg_u32;
    int k_scan_addback_u32;
    int k_radix_histogram_u64;
    int k_radix_scatter_u64;

    /* Scratch device buffers, grown on demand (ensure_capacity style,
     * mirroring collision_engine's own struct -- see ocl_shared.h's
     * design note on why this shape was chosen in Phase 1). Callers
     * never touch these directly. */
    cl_mem scratch_a;   size_t scratch_a_cap;    /* generic uint32: reduce-max ping-pong, radix histogram/offsets */
    cl_mem scratch_b;   size_t scratch_b_cap;    /* generic uint32: reduce-max ping-pong (2nd buffer) */
    cl_mem scratch_c;   size_t scratch_c_cap;    /* generic uint32: ocl_scan_exclusive_u32's OWN block-sums array --
                                                     kept separate from scratch_a/b so radix sort can scan its
                                                     histogram (which lives in scratch_a) without the scan
                                                     primitive clobbering the very data it's scanning */
    cl_mem scratch_keys;size_t scratch_keys_cap; /* ulong, radix sort's second ping-pong key buffer */
} ocl_primitives_t;

/* Builds the program (from ocl_primitives_kernels.cl, read from disk)
 * and creates every kernel object. cl_source_path is the path to
 * ocl_primitives_kernels.cl (kept as a parameter rather than a fixed
 * relative path since Phase 3+ may relocate it). */
int  ocl_primitives_init(ocl_primitives_t *p, ocl_gerbicz_device_t *gd,
                          const char *cl_source_path, const char *cache_dir);
void ocl_primitives_free(ocl_primitives_t *p);

/* out[i] = value for i in [0, n). buf_off is in elements. */
int ocl_fill_u32(ocl_primitives_t *p, cl_mem buf, cl_ulong buf_off, uint32_t n, uint32_t value);
int ocl_fill_u64(ocl_primitives_t *p, cl_mem buf, cl_ulong buf_off, uint32_t n, uint64_t value);

/* *out_max = max(in[0..n)). n == 0 is an error (no identity element
 * defined here -- callers with a possibly-empty range should check n
 * themselves, matching collision_engine's own usage where NUM_BUCKETS
 * is always > 0). */
int ocl_reduce_max_u32(ocl_primitives_t *p, cl_mem in, cl_ulong in_off,
                        uint32_t n, uint32_t *out_max);

/* In-place exclusive scan of buf[0..n). n must be <= OCL_PRIM_SCAN_MAX_N
 * (see above) -- returns -1 and does nothing if it isn't. */
int ocl_scan_exclusive_u32(ocl_primitives_t *p, cl_mem buf, cl_ulong buf_off, uint32_t n);

/* Sorts keys_in[0..n) into keys_out[0..n) ascending, treating each key
 * as a plain unsigned 64-bit integer (see the header comment above for
 * why that's correct for this codebase's actual key representation).
 * keys_in and keys_out must be different allocations (unlike the scan
 * primitives above, radix sort needs true ping-pong buffers across its
 * 8 passes) and both must have at least n elements from their
 * respective offsets. n must be <= OCL_PRIM_SCAN_MAX_N / OCL_PRIM_RADIX_BINS
 * (the per-pass histogram, RADIX_BINS entries per workgroup, must itself
 * fit the scan primitive -- see ocl_primitives.c for the exact check). */
int ocl_radix_sort_u64(ocl_primitives_t *p,
                        cl_mem keys_in,  cl_ulong keys_in_off,
                        cl_mem keys_out, cl_ulong keys_out_off,
                        uint32_t n);

#ifdef __cplusplus
}
#endif

#endif /* _OCL_PRIMITIVES_H_ */
