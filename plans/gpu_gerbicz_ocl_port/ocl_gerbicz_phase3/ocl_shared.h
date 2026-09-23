/*--------------------------------------------------------------------
 * ocl_shared.h  --  shared OpenCL context/program/cache infrastructure
 *
 * Phase 1 of the gpu_gerbicz -> OpenCL port.
 *
 * This module factors out the device-selection / program-build / binary
 * cache / per-thread queue+kernel machinery that gpu_cofactorization_cl.c
 * already implements inline for its own ECM/P-1 kernels, so that the new
 * ocl_gerbicz collision engine (Phase 3+) can reuse exactly the same
 * lifecycle shape without duplicating it, and so gpu_cofactorization_cl.c
 * can (optionally, later) be migrated onto it to pick up the fixed cache
 * key (see ocl_cache_key_t below).
 *
 * Nothing here is gerbicz-specific. Gerbicz-specific pieces (the actual
 * kernel names, arg descriptors, and the sieve/collision program split)
 * live in ocl_gerbicz_ctx.h/.c.
 *
 * Locked decisions this module implements (see OCL_GERBICZ_PLAN.md §2):
 *   1. OpenCL >= 2.0 required; sub-groups optional and gated
 *      (ocl_check_min_version, ocl_query_subgroup_support).
 *   2. Static link, .cl sources loaded and built at runtime; binaries are
 *      cached to disk (ocl_build_program_cached).
 *   5. One cl_context + one built cl_program per device (ocl_device_t);
 *      per-thread cl_command_queue and cl_kernel objects (ocl_thread_t).
 *   4. Base-offset kernel-argument convention -- see the comment block
 *      in ocl_gerbicz_ctx.h. This module does not encode that convention
 *      itself (it is a per-kernel-signature convention, not a runtime
 *      object), but ocl_thread_set_args() below is the single place that
 *      calls clSetKernelArg, so it is the enforcement point.
 *--------------------------------------------------------------------*/

#ifndef _OCL_SHARED_H_
#define _OCL_SHARED_H_

#include <stddef.h>
#include <stdint.h>

/* Locked decision #1: OpenCL >= 2.0 is the floor. We pin the HOST API
 * surface (this macro) to 300, not 200: querying whether a device
 * actually accepts "-cl-std=CL2.0" needs the 3.0 host query
 * CL_DEVICE_OPENCL_C_ALL_VERSIONS (see ocl_device_query) even on a
 * 2.0-only target device -- clGetDeviceInfo simply returns
 * CL_INVALID_VALUE there and we fall back. This also silences the 1.x
 * deprecation warnings (e.g. clCreateCommandQueue) since we use the
 * 2.0 replacements directly. Define before this header if a
 * translation unit needs a different target for some reason. */
#ifndef CL_TARGET_OPENCL_VERSION
#define CL_TARGET_OPENCL_VERSION 300
#endif

#ifdef __APPLE__
#  include <OpenCL/opencl.h>
#else
#  include <CL/cl.h>
#endif

#include "ocl_xface.h"   /* gpu_info_t, gpu_arg_type_list_t, OCL_TRY, ... */

#ifdef __cplusplus
extern "C" {
#endif

/* -----------------------------------------------------------------------
 * clGetErrorString
 *
 * ocl_xface.h declares this (OCL_TRY expands to a call to it) but does
 * not define it anywhere in the files handed to this phase, and
 * gpu_cofactorization_cl.c calls it without defining it either -- so it
 * must already live in an ocl_xface.c we were not given, OR it simply
 * does not exist yet. Phase 1 cannot tell which from the sources on
 * hand. We provide our own definition here, guarded so it becomes a
 * silent no-op duplicate (link error) if one already exists elsewhere --
 * see STATUS "Open issues" for the flag and what to do if the real repo
 * already has one.
 * --------------------------------------------------------------------- */
const char *clGetErrorString(cl_int err);

/* -----------------------------------------------------------------------
 * Program identity / cache key
 *
 * Fixes the Phase 0 Fact: gpu_cofactorization_cl.c's existing cache file
 * is keyed on device name ONLY ("ocl_ecm_<devname>.bin"), so an edited
 * .cl source or changed build flags silently loads a stale binary. Here
 * the key is device name + driver version + build-options string + a
 * source hash (FNV-1a 64, order-sensitive over the concatenated source
 * texts) so any of those changing forces a rebuild.
 * --------------------------------------------------------------------- */
typedef struct {
    char     device_name[128];
    char     driver_version[64];
    char     build_options[256];
    uint64_t source_hash;      /* FNV-1a 64 over concatenated sources */
} ocl_cache_key_t;

/* FNV-1a 64. Exposed so callers can hash their own source text before
 * filling in ocl_cache_key_t.source_hash. */
uint64_t ocl_fnv1a64(const void *data, size_t len);

/* Reads a whole text file into a malloc'd, NUL-terminated buffer (caller
 * frees). Returns NULL on any I/O error. This is the "runtime-loaded .cl
 * files" half of locked decision #2 -- callers pass the returned strings
 * straight to ocl_build_program_cached(). */
char *ocl_read_text_file(const char *path);

/* Formats a filesystem-safe cache file name from a key, e.g.
 *   "oclcache_gfx1031_5f3a1c2b9e7d4410.bin"
 * (device name sanitised, driver+options+source folded into one hex
 * hash so the filename stays short). buf must be >= 128 bytes. */
void ocl_cache_key_filename(const ocl_cache_key_t *key,
                             const char *prefix,
                             char *buf, size_t buf_sz);

/* -----------------------------------------------------------------------
 * Device-level object: one per physical device, shared by every thread.
 * Callers embed/point to a gpu_info_t (from ocl_xface.h / gpu_init())
 * for device_handle/platform_handle and capability fields.
 * --------------------------------------------------------------------- */
typedef struct {
    gpu_info_t   *info;            /* not owned; caller retains ownership */
    cl_context    context;

    /* OpenCL >= 2.0 check result and sub-group support, cached here so
     * every thread and every program build can consult it without
     * re-querying the device. */
    int           version_major;
    int           version_minor;
    int           has_subgroups;   /* cl_khr_subgroups or >= CL 2.1 */

    /* Queried limits, needed from Phase 1 on so later phases don't have
     * to re-derive them (task 1.2 / plan Phase 1 bullet list). */
    cl_ulong      local_mem_size;       /* CL_DEVICE_LOCAL_MEM_SIZE */
    cl_ulong      max_mem_alloc_size;   /* CL_DEVICE_MAX_MEM_ALLOC_SIZE */
    size_t        max_work_group_size;  /* CL_DEVICE_MAX_WORK_GROUP_SIZE */
    size_t        pref_wg_multiple;     /* kernel-specific; 0 until a
                                            kernel is queried, see
                                            ocl_query_preferred_wg_multiple */

    /* The actual "-cl-std=CLx.y" string this device will accept, e.g.
     * "CL2.0". NOT the same thing as version_major.version_minor (the
     * device-slash-platform version): a device can report OpenCL 3.0 and
     * still reject "-cl-std=CL2.0" if its CL_DEVICE_OPENCL_C_ALL_VERSIONS
     * list happens to skip 2.0 (observed on PoCL in this sandbox -- see
     * STATUS Facts). ocl_build_program_cached() uses this instead of
     * trusting a hardcoded "CL2.0". */
    char          cl_c_std[32];
} ocl_device_t;

/* Fills version_major/minor, has_subgroups, local_mem_size,
 * max_mem_alloc_size, max_work_group_size from info->device_handle.
 * Does NOT create the context (call ocl_device_create_context after,
 * or use ocl_device_init which does both). Returns 0 on success. */
int ocl_device_query(ocl_device_t *dev, gpu_info_t *info);

/* Creates dev->context from dev->info->device_handle/platform_handle.
 * Returns 0 on success. */
int ocl_device_create_context(ocl_device_t *dev);

/* ocl_device_query + ocl_device_create_context in one call. */
int ocl_device_init(ocl_device_t *dev, gpu_info_t *info);

void ocl_device_free(ocl_device_t *dev);   /* releases context only */

/* Returns 1 if dev->version_{major,minor} >= (min_major, min_minor). */
int ocl_check_min_version(const ocl_device_t *dev,
                           int min_major, int min_minor);

/* Preferred kernel work-group-size multiple (wave/warp size equivalent).
 * Must be called after the kernel is built (needs a cl_kernel handle).
 * Do NOT assume 32 -- AMD RDNA2 wavefronts are 32 or 64 depending on
 * kernel occupancy and compile mode. */
size_t ocl_query_preferred_wg_multiple(cl_kernel kernel, cl_device_id device);

/* -----------------------------------------------------------------------
 * Program build with disk cache.
 *
 * num_sources source strings are concatenated (in order) to compute the
 * hash and are passed to clCreateProgramWithSource as separate strings
 * (so line numbers in compiler diagnostics stay meaningful).
 *
 * cache_dir may be NULL for the current directory. extra_build_opts is
 * appended after the module's own "-cl-std=CLx.y" (pass NULL for none);
 * this is where a caller adds e.g. "-DHAVE_SUBGROUPS=1" or the
 * collision_bucket.h-derived -D flags (task 1.2 bullet 5 / plan Phase 1
 * bullet "Generate -D options from collision_bucket.h").
 *
 * Returns the built cl_program, or NULL on failure (build log is
 * printed to stderr; this function does not exit(-1) itself so callers
 * can decide whether a failed build is fatal).
 * --------------------------------------------------------------------- */
cl_program ocl_build_program_cached(ocl_device_t *dev,
                                     const char **sources,
                                     int num_sources,
                                     const char *cl_std,          /* e.g. "CL2.0" */
                                     const char *extra_build_opts,/* may be NULL */
                                     const char *cache_dir,       /* may be NULL */
                                     const char *cache_prefix);   /* e.g. "gerbicz" */

/* -----------------------------------------------------------------------
 * Thread-level object: per-thread queue + kernel handles.
 * clSetKernelArg is not thread-safe on a shared kernel object, so every
 * thread gets its OWN cl_kernel for each entry point (created via
 * clCreateKernel from the shared, already-built cl_program), matching
 * ocl_xface.h's existing gpu_launch_t / gpu_launch_init pattern.
 * --------------------------------------------------------------------- */
typedef struct {
    ocl_device_t     *dev;          /* not owned */
    cl_command_queue  queue;
    gpu_launch_t     *launch;       /* array, num_kernels entries */
    int               num_kernels;
} ocl_thread_t;

/* Creates the queue and calls gpu_launch_init() (from ocl_xface.h) once
 * per (kernel_names[i], arg_descs[i]) pair against program. */
int ocl_thread_init(ocl_thread_t *th, ocl_device_t *dev, cl_program program,
                     const char **kernel_names,
                     const gpu_arg_type_list_t *arg_descs,
                     int num_kernels);

void ocl_thread_free(ocl_thread_t *th);

#ifdef __cplusplus
}
#endif

#endif /* _OCL_SHARED_H_ */
