/*--------------------------------------------------------------------
 * ocl_gerbicz_ctx.h  --  Gerbicz-engine-specific OpenCL plumbing.
 *
 * Built on top of ocl_shared.h. This is the header later phases include;
 * ocl_shared.h stays generic (also usable by gpu_cofactorization_cl.c).
 *
 * =======================================================================
 * BASE-OFFSET KERNEL-ARGUMENT CONVENTION  (locked decision #4)
 * =======================================================================
 * No sub-buffers, no SVM. A kernel that must operate on a sub-range of a
 * larger buffer (e.g. batch B's slice of d_candidate_keys) receives the
 * FULL buffer plus a separate offset scalar, and adds the offset to the
 * pointer itself as its first statement:
 *
 *   __kernel void foo(__global ulong *keys, ulong keys_off, ...) {
 *       keys += keys_off;
 *       ... use keys[i] as if it were keys[0] of the sub-range ...
 *   }
 *
 * Rules:
 *   1. Offsets are in ELEMENTS of the pointee type, never bytes.
 *   2. Every offset argument is cl_ulong (8 bytes) host-side and `ulong`
 *      device-side, even though today's buffers (CANDIDATE_CAP = 2^22,
 *      NUM_BUCKETS * max_per_bucket) fit in 32 bits -- so a later phase
 *      that grows a buffer past 4G elements does not need to touch the
 *      call sites, only the value passed.
 *   3. In the kernel's parameter list and in the gpu_arg_type_list_t /
 *      gpu_arg_t arrays that describe it, the offset argument
 *      IMMEDIATELY FOLLOWS the buffer it offsets: (..., GPU_ARG_PTR,
 *      GPU_ARG_UINT64, ...). This keeps host call sites and kernel
 *      signatures visually paired and lets a future arg-list generator
 *      (if one is ever written) infer pairing from position alone.
 *   4. A kernel with several base-offset buffers repeats the pattern for
 *      each one, in the buffer's own declaration order -- not grouped
 *      at the end of the arg list.
 *   5. A buffer used only in full (never sliced) takes no offset arg.
 *   6. Host side, the offset is always "index into the allocation this
 *      phase's collision_engine (later, ocl_collision_engine) owns" --
 *      i.e. it plays exactly the role CUDA's raw pointer arithmetic
 *      (`p_out += off`) played, just expressed as an explicit argument
 *      because OpenCL's cl_mem is an opaque handle, not a pointer value
 *      the host can add to.
 *
 * This must be settled now (Phase 1) because every Phase 2/3 kernel
 * signature depends on it; changing it later means re-touching every
 * kernel and every gpu_arg_type_list_t entry.
 *
 * =======================================================================
 * SUB-GROUP GATING  (locked decisions #1, #6)
 * =======================================================================
 * OpenCL 2.0 is the floor; sub-groups (cl_khr_subgroups) are optional.
 * Gating is a BUILD-TIME macro, decided at program-build time from the
 * already-queried ocl_device_t::has_subgroups (see ocl_shared.h), not a
 * runtime branch inside the kernel:
 *
 *   ocl_build_program_cached(..., extra_build_opts, ...)
 *   where extra_build_opts includes "-DHAVE_SUBGROUPS=1" or "=0"
 *   depending on dev->has_subgroups.
 *
 * .cl sources guard the sub-group path with:
 *   #if HAVE_SUBGROUPS
 *     ... sub_group_broadcast / sub_group_reduce_add path ...
 *   #else
 *     ... local-memory + barrier path (the baseline, decision #6) ...
 *   #endif
 *
 * Rationale for build-time over runtime: the two code paths use
 * different work-group-size assumptions and different local-memory
 * layouts; branching on it inside the kernel would mean allocating
 * local memory for both paths simultaneously. Phase 3's baseline
 * collision kernels are no-sub-group per decision #6 regardless of what
 * the device supports -- HAVE_SUBGROUPS is here for the *optional fast
 * path* mentioned in decision #1, not required for Phase 3 itself.
 * --------------------------------------------------------------------- */

#ifndef _OCL_GERBICZ_CTX_H_
#define _OCL_GERBICZ_CTX_H_

#include "ocl_shared.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Two separate programs (plan Phase 1 bullet "Separate programs for the
 * sieve kernels and the collision kernels"): they have unrelated build
 * options (the sieve kernels never need HAVE_SUBGROUPS, for instance)
 * and unrelated release cadences (Phase 3 vs Phase 4). One ocl_device_t
 * (one context) is still shared between them per decision #5. */
typedef struct {
    ocl_device_t  dev;
    cl_program    program_collision;   /* built in Phase 3 */
    cl_program    program_sieve;       /* built in Phase 4 */
} ocl_gerbicz_device_t;

/* Per host thread: one ocl_thread_t per program actually in use by that
 * thread (a thread that only launches collision kernels doesn't need
 * sieve kernel objects, and vice versa; Phase 5 decides which threads
 * need which). Phase 1 does not populate collision/sieve thread_t
 * bodies yet -- there are no real kernels to look up until Phase 3/4 --
 * this struct just fixes the shape so Phase 3/5 don't redesign it. */
typedef struct {
    ocl_thread_t  collision;   /* .launch == NULL / num_kernels == 0 until Phase 3 */
    ocl_thread_t  sieve;       /* .launch == NULL / num_kernels == 0 until Phase 4 */
} ocl_gerbicz_thread_t;

int  ocl_gerbicz_device_init(ocl_gerbicz_device_t *gd, gpu_info_t *info);
void ocl_gerbicz_device_free(ocl_gerbicz_device_t *gd);

/* Builds the -DHAVE_SUBGROUPS=0/1 flag string from gd->dev.has_subgroups
 * and appends collision_bucket.h's LOG2_NUM_BUCKETS/BUCKET_HASH_MIX as
 * further -D flags (plan Phase 1 bullet), so both are baked into the
 * program instead of duplicated as .cl-side #defines that could drift
 * from collision_bucket.h (the header explicitly warns against drift
 * between translation units). buf must be >= 160 bytes. */
void ocl_gerbicz_build_opts(const ocl_gerbicz_device_t *gd, char *buf, size_t buf_sz);

#ifdef __cplusplus
}
#endif

#endif /* _OCL_GERBICZ_CTX_H_ */
