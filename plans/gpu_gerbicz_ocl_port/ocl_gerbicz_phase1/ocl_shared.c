/*--------------------------------------------------------------------
 * ocl_shared.c  --  see ocl_shared.h for the design rationale.
 *--------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ocl_shared.h"

/* -----------------------------------------------------------------------
 * clGetErrorString
 *
 * See the note in ocl_shared.h. Only the codes that actually appear in
 * this module (and the ones gpu_cofactorization_cl.c already triggers)
 * are named individually; everything else falls back to the numeric
 * value so OCL_TRY never prints nothing.
 * --------------------------------------------------------------------- */
const char *
clGetErrorString(cl_int err)
{
    switch (err) {
    case CL_SUCCESS:                          return "CL_SUCCESS";
    case CL_DEVICE_NOT_FOUND:                 return "CL_DEVICE_NOT_FOUND";
    case CL_DEVICE_NOT_AVAILABLE:             return "CL_DEVICE_NOT_AVAILABLE";
    case CL_COMPILER_NOT_AVAILABLE:           return "CL_COMPILER_NOT_AVAILABLE";
    case CL_MEM_OBJECT_ALLOCATION_FAILURE:    return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
    case CL_OUT_OF_RESOURCES:                 return "CL_OUT_OF_RESOURCES";
    case CL_OUT_OF_HOST_MEMORY:               return "CL_OUT_OF_HOST_MEMORY";
    case CL_PROFILING_INFO_NOT_AVAILABLE:     return "CL_PROFILING_INFO_NOT_AVAILABLE";
    case CL_MEM_COPY_OVERLAP:                 return "CL_MEM_COPY_OVERLAP";
    case CL_IMAGE_FORMAT_MISMATCH:            return "CL_IMAGE_FORMAT_MISMATCH";
    case CL_IMAGE_FORMAT_NOT_SUPPORTED:       return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
    case CL_BUILD_PROGRAM_FAILURE:            return "CL_BUILD_PROGRAM_FAILURE";
    case CL_MAP_FAILURE:                      return "CL_MAP_FAILURE";
    case CL_MISALIGNED_SUB_BUFFER_OFFSET:     return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
    case CL_INVALID_VALUE:                    return "CL_INVALID_VALUE";
    case CL_INVALID_DEVICE_TYPE:              return "CL_INVALID_DEVICE_TYPE";
    case CL_INVALID_PLATFORM:                 return "CL_INVALID_PLATFORM";
    case CL_INVALID_DEVICE:                   return "CL_INVALID_DEVICE";
    case CL_INVALID_CONTEXT:                  return "CL_INVALID_CONTEXT";
    case CL_INVALID_QUEUE_PROPERTIES:         return "CL_INVALID_QUEUE_PROPERTIES";
    case CL_INVALID_COMMAND_QUEUE:            return "CL_INVALID_COMMAND_QUEUE";
    case CL_INVALID_HOST_PTR:                 return "CL_INVALID_HOST_PTR";
    case CL_INVALID_MEM_OBJECT:               return "CL_INVALID_MEM_OBJECT";
    case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:  return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
    case CL_INVALID_IMAGE_SIZE:               return "CL_INVALID_IMAGE_SIZE";
    case CL_INVALID_SAMPLER:                  return "CL_INVALID_SAMPLER";
    case CL_INVALID_BINARY:                   return "CL_INVALID_BINARY";
    case CL_INVALID_BUILD_OPTIONS:            return "CL_INVALID_BUILD_OPTIONS";
    case CL_INVALID_PROGRAM:                  return "CL_INVALID_PROGRAM";
    case CL_INVALID_PROGRAM_EXECUTABLE:       return "CL_INVALID_PROGRAM_EXECUTABLE";
    case CL_INVALID_KERNEL_NAME:              return "CL_INVALID_KERNEL_NAME";
    case CL_INVALID_KERNEL_DEFINITION:        return "CL_INVALID_KERNEL_DEFINITION";
    case CL_INVALID_KERNEL:                   return "CL_INVALID_KERNEL";
    case CL_INVALID_ARG_INDEX:                return "CL_INVALID_ARG_INDEX";
    case CL_INVALID_ARG_VALUE:                return "CL_INVALID_ARG_VALUE";
    case CL_INVALID_ARG_SIZE:                 return "CL_INVALID_ARG_SIZE";
    case CL_INVALID_KERNEL_ARGS:              return "CL_INVALID_KERNEL_ARGS";
    case CL_INVALID_WORK_DIMENSION:           return "CL_INVALID_WORK_DIMENSION";
    case CL_INVALID_WORK_GROUP_SIZE:          return "CL_INVALID_WORK_GROUP_SIZE";
    case CL_INVALID_WORK_ITEM_SIZE:           return "CL_INVALID_WORK_ITEM_SIZE";
    case CL_INVALID_GLOBAL_OFFSET:            return "CL_INVALID_GLOBAL_OFFSET";
    case CL_INVALID_EVENT_WAIT_LIST:          return "CL_INVALID_EVENT_WAIT_LIST";
    case CL_INVALID_EVENT:                    return "CL_INVALID_EVENT";
    case CL_INVALID_OPERATION:                return "CL_INVALID_OPERATION";
    case CL_INVALID_BUFFER_SIZE:              return "CL_INVALID_BUFFER_SIZE";
    case CL_INVALID_GLOBAL_WORK_SIZE:         return "CL_INVALID_GLOBAL_WORK_SIZE";
    default: {
        static char buf[32];   /* NOTE: not thread-safe; see STATUS. */
        snprintf(buf, sizeof(buf), "CL error %d", (int)err);
        return buf;
    }
    }
}

/* -----------------------------------------------------------------------
 * FNV-1a 64
 * --------------------------------------------------------------------- */
uint64_t
ocl_fnv1a64(const void *data, size_t len)
{
    const unsigned char *p = (const unsigned char *)data;
    uint64_t h = 0xcbf29ce484222325ULL;
    size_t i;
    for (i = 0; i < len; i++) {
        h ^= (uint64_t)p[i];
        h *= 0x100000001b3ULL;
    }
    return h;
}

static void
sanitize_into(const char *src, char *dst, size_t dst_sz)
{
    size_t lim = dst_sz - 1;
    size_t n = 0;
    while (*src && n < lim) {
        char c = *src++;
        int ok = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
                 (c >= '0' && c <= '9') || c == '_' || c == '-';
        dst[n++] = ok ? c : '_';
    }
    dst[n] = '\0';
}

void
ocl_cache_key_filename(const ocl_cache_key_t *key,
                        const char *prefix,
                        char *buf, size_t buf_sz)
{
    char safe_name[128];
    uint64_t combined;
    char combined_buf[128 + 64 + 256 + 32];
    int n;

    sanitize_into(key->device_name, safe_name, sizeof(safe_name));

    n = snprintf(combined_buf, sizeof(combined_buf), "%s|%s|%s|%016llx",
                 key->device_name, key->driver_version, key->build_options,
                 (unsigned long long)key->source_hash);
    combined = ocl_fnv1a64(combined_buf, (size_t)(n > 0 ? n : 0));

    snprintf(buf, buf_sz, "%s_%s_%016llx.bin",
             prefix ? prefix : "oclcache", safe_name,
             (unsigned long long)combined);
}

/* -----------------------------------------------------------------------
 * Device query / context
 * --------------------------------------------------------------------- */
int
ocl_device_query(ocl_device_t *dev, gpu_info_t *info)
{
    cl_device_id d = info->device_handle;
    char ext[4096];
    size_t ext_len = 0;

    memset(dev, 0, sizeof(*dev));
    dev->info = info;

    /* gpu_info_t is assumed already populated by gpu_init() (existing
     * cofactorization code path) -- reuse rather than re-query. */
    dev->version_major = info->compute_version_major;
    dev->version_minor = info->compute_version_minor;

    if (clGetDeviceInfo(d, CL_DEVICE_LOCAL_MEM_SIZE,
            sizeof(dev->local_mem_size), &dev->local_mem_size, NULL) != CL_SUCCESS)
        return -1;
    if (clGetDeviceInfo(d, CL_DEVICE_MAX_MEM_ALLOC_SIZE,
            sizeof(dev->max_mem_alloc_size), &dev->max_mem_alloc_size, NULL) != CL_SUCCESS)
        return -1;
    if (clGetDeviceInfo(d, CL_DEVICE_MAX_WORK_GROUP_SIZE,
            sizeof(dev->max_work_group_size), &dev->max_work_group_size, NULL) != CL_SUCCESS)
        return -1;

    /* Sub-group support: OpenCL 2.0 devices expose it (if at all) only
     * via the cl_khr_subgroups extension string; do not assume the 2.1
     * query mechanism (CL_DEVICE_MAX_NUM_SUB_GROUPS) exists. */
    if (clGetDeviceInfo(d, CL_DEVICE_EXTENSIONS, sizeof(ext), ext, &ext_len) == CL_SUCCESS) {
        dev->has_subgroups =
            (strstr(ext, "cl_khr_subgroups") != NULL) ||
            (strstr(ext, "cl_intel_subgroups") != NULL);
    } else {
        dev->has_subgroups = 0;
    }

    dev->pref_wg_multiple = 0; /* filled in lazily, needs a built kernel */

    /* Determine the actual "-cl-std=CLx.y" the device will accept.
     * CL_DEVICE_VERSION (already reflected in version_major/minor) is
     * the *platform/device* version and is NOT reliable proof that a
     * given -cl-std value builds -- on CL3.0 devices the accepted C
     * standards are whatever CL_DEVICE_OPENCL_C_ALL_VERSIONS lists, and
     * that list can skip versions (observed: PoCL 5.0 reports device
     * version 3.0 and lists OpenCL C 1.0/1.1/1.2/3.0 -- no 2.0 entry at
     * all -- so "-cl-std=CL2.0" is rejected even though "CL3.0" and
     * "CL1.2" both work). Prefer the highest version >= 2.0 that the
     * device actually lists; if CL_DEVICE_OPENCL_C_ALL_VERSIONS isn't
     * queryable (pre-3.0 ICDs), fall back to parsing the single
     * CL_DEVICE_OPENCL_C_VERSION string. */
    {
        cl_name_version versions[32];
        size_t ret_bytes = 0;
        int best_major = -1, best_minor = -1;

        if (clGetDeviceInfo(d, CL_DEVICE_OPENCL_C_ALL_VERSIONS,
                sizeof(versions), versions, &ret_bytes) == CL_SUCCESS && ret_bytes > 0) {
            size_t i, n = ret_bytes / sizeof(cl_name_version);
            for (i = 0; i < n; i++) {
                int maj = (int)((versions[i].version >> 22) & 0x3ff);
                int min = (int)((versions[i].version >> 12) & 0x3ff);
                if (maj < 2)
                    continue; /* below our floor, decision #1 */
                if (maj > best_major || (maj == best_major && min > best_minor)) {
                    best_major = maj;
                    best_minor = min;
                }
            }
        }
        if (best_major < 0) {
            char cver[128];
            int maj = 0, min = 0;
            if (clGetDeviceInfo(d, CL_DEVICE_OPENCL_C_VERSION, sizeof(cver), cver, NULL) == CL_SUCCESS &&
                sscanf(cver, "OpenCL C %d.%d", &maj, &min) == 2) {
                best_major = maj;
                best_minor = min;
            }
        }
        if (best_major >= 0)
            snprintf(dev->cl_c_std, sizeof(dev->cl_c_std), "CL%d.%d", best_major, best_minor);
        else
            dev->cl_c_std[0] = '\0'; /* couldn't determine; caller must supply one */
    }

    return 0;
}

int
ocl_device_create_context(ocl_device_t *dev)
{
    cl_int err;
    cl_device_id d = dev->info->device_handle;
    cl_context_properties props[] = {
        CL_CONTEXT_PLATFORM, (cl_context_properties)dev->info->platform_handle, 0
    };

    dev->context = clCreateContext(props, 1, &d, NULL, NULL, &err);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "ocl_device_create_context: clCreateContext failed: %s\n",
            clGetErrorString(err));
        return -1;
    }
    return 0;
}

int
ocl_device_init(ocl_device_t *dev, gpu_info_t *info)
{
    if (ocl_device_query(dev, info) != 0)
        return -1;
    return ocl_device_create_context(dev);
}

void
ocl_device_free(ocl_device_t *dev)
{
    if (dev->context)
        clReleaseContext(dev->context);
    dev->context = NULL;
}

int
ocl_check_min_version(const ocl_device_t *dev, int min_major, int min_minor)
{
    if (dev->version_major != min_major)
        return dev->version_major > min_major;
    return dev->version_minor >= min_minor;
}

size_t
ocl_query_preferred_wg_multiple(cl_kernel kernel, cl_device_id device)
{
    size_t multiple = 0;
    if (clGetKernelWorkGroupInfo(kernel, device,
            CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
            sizeof(multiple), &multiple, NULL) != CL_SUCCESS)
        return 0;
    return multiple;
}

/* -----------------------------------------------------------------------
 * Program build with disk cache
 * --------------------------------------------------------------------- */
char *
ocl_read_text_file(const char *path)
{
    FILE *f = fopen(path, "rb");
    long sz;
    char *buf;

    if (!f)
        return NULL;
    fseek(f, 0, SEEK_END);
    sz = ftell(f);
    rewind(f);
    buf = (char *)malloc((size_t)sz + 1);
    if (!buf) { fclose(f); return NULL; }
    if (fread(buf, 1, (size_t)sz, f) != (size_t)sz) {
        free(buf);
        fclose(f);
        return NULL;
    }
    buf[sz] = '\0';
    fclose(f);
    return buf;
}

static int
try_load_cached_binary(ocl_device_t *dev, const char *path,
                        const char *build_opts, cl_program *out_program)
{
    FILE *f = fopen(path, "rb");
    long bin_sz;
    unsigned char *bin;
    cl_int err, bin_status = CL_SUCCESS;
    cl_device_id d = dev->info->device_handle;
    cl_program prog;

    if (!f)
        return -1;
    fseek(f, 0, SEEK_END);
    bin_sz = ftell(f);
    rewind(f);
    if (bin_sz <= 0) { fclose(f); return -1; }

    bin = (unsigned char *)malloc((size_t)bin_sz);
    if (!bin || fread(bin, 1, (size_t)bin_sz, f) != (size_t)bin_sz) {
        free(bin);
        fclose(f);
        return -1;
    }
    fclose(f);

    {
        const unsigned char *bin_ptr = bin;
        size_t bin_len = (size_t)bin_sz;
        prog = clCreateProgramWithBinary(dev->context, 1, &d,
                &bin_len, &bin_ptr, &bin_status, &err);
    }
    free(bin);

    if (err != CL_SUCCESS || bin_status != CL_SUCCESS) {
        if (prog) clReleaseProgram(prog);
        return -1;
    }

    err = clBuildProgram(prog, 1, &d, build_opts, NULL, NULL);
    if (err != CL_SUCCESS) {
        clReleaseProgram(prog);
        return -1;   /* stale/incompatible cache -- caller recompiles */
    }

    *out_program = prog;
    return 0;
}

static void
save_binary_to_cache(cl_program program, const char *path)
{
    size_t bin_sz = 0;
    unsigned char *bin;
    unsigned char *bins[1];
    FILE *f;

    if (clGetProgramInfo(program, CL_PROGRAM_BINARY_SIZES,
            sizeof(bin_sz), &bin_sz, NULL) != CL_SUCCESS || bin_sz == 0)
        return;

    bin = (unsigned char *)malloc(bin_sz);
    if (!bin)
        return;
    bins[0] = bin;
    if (clGetProgramInfo(program, CL_PROGRAM_BINARIES,
            sizeof(bins), bins, NULL) != CL_SUCCESS) {
        free(bin);
        return;
    }

    f = fopen(path, "wb");
    if (f) {
        fwrite(bin, 1, bin_sz, f);
        fclose(f);
    } else {
        fprintf(stderr, "warning: could not write OpenCL binary cache '%s'\n", path);
    }
    free(bin);
}

cl_program
ocl_build_program_cached(ocl_device_t *dev,
                          const char **sources,
                          int num_sources,
                          const char *cl_std,
                          const char *extra_build_opts,
                          const char *cache_dir,
                          const char *cache_prefix)
{
    cl_device_id d = dev->info->device_handle;
    cl_int err;
    cl_program program = NULL;
    char build_opts[512];
    char driver_version[64];
    ocl_cache_key_t key;
    char cache_name[192];
    char cache_path[512];
    int i;
    size_t *lengths;
    uint64_t hash;

    snprintf(build_opts, sizeof(build_opts), "-cl-std=%s%s%s",
             cl_std ? cl_std : "CL2.0",
             extra_build_opts ? " " : "",
             extra_build_opts ? extra_build_opts : "");

    if (clGetDeviceInfo(d, CL_DRIVER_VERSION, sizeof(driver_version),
            driver_version, NULL) != CL_SUCCESS)
        strcpy(driver_version, "unknown");

    /* Hash the concatenated sources (order matters -- matches how they
     * are handed to clCreateProgramWithSource below). */
    hash = 0xcbf29ce484222325ULL;
    for (i = 0; i < num_sources; i++) {
        size_t l = strlen(sources[i]);
        /* fold each source in with FNV-1a chaining so hash depends on
         * source content AND source order/boundaries */
        uint64_t h2 = ocl_fnv1a64(sources[i], l);
        hash ^= h2 + 0x9e3779b97f4a7c15ULL + (hash << 6) + (hash >> 2);
    }

    memset(&key, 0, sizeof(key));
    strncpy(key.device_name, dev->info->name, sizeof(key.device_name) - 1);
    strncpy(key.driver_version, driver_version, sizeof(key.driver_version) - 1);
    strncpy(key.build_options, build_opts, sizeof(key.build_options) - 1);
    key.source_hash = hash;

    ocl_cache_key_filename(&key, cache_prefix, cache_name, sizeof(cache_name));
    if (cache_dir && cache_dir[0])
        snprintf(cache_path, sizeof(cache_path), "%s/%s", cache_dir, cache_name);
    else
        snprintf(cache_path, sizeof(cache_path), "%s", cache_name);

    if (try_load_cached_binary(dev, cache_path, build_opts, &program) == 0) {
        fprintf(stderr, "ocl: loaded cached binary '%s'\n", cache_path);
        return program;
    }

    lengths = (size_t *)malloc((size_t)num_sources * sizeof(size_t));
    for (i = 0; i < num_sources; i++)
        lengths[i] = strlen(sources[i]);

    program = clCreateProgramWithSource(dev->context, (cl_uint)num_sources,
            sources, lengths, &err);
    free(lengths);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "ocl_build_program_cached: clCreateProgramWithSource failed: %s\n",
                clGetErrorString(err));
        return NULL;
    }

    err = clBuildProgram(program, 1, &d, build_opts, NULL, NULL);
    if (err != CL_SUCCESS) {
        size_t log_sz = 0;
        char *log;
        clGetProgramBuildInfo(program, d, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_sz);
        log = (char *)malloc(log_sz + 1);
        clGetProgramBuildInfo(program, d, CL_PROGRAM_BUILD_LOG, log_sz, log, NULL);
        log[log_sz] = '\0';
        fprintf(stderr, "ocl_build_program_cached: clBuildProgram failed:\n%s\n", log);
        free(log);
        clReleaseProgram(program);
        return NULL;
    }

    fprintf(stderr, "ocl: compiled, saving cache '%s'\n", cache_path);
    save_binary_to_cache(program, cache_path);
    return program;
}

/* -----------------------------------------------------------------------
 * Per-thread queue + kernel objects
 * --------------------------------------------------------------------- */
int
ocl_thread_init(ocl_thread_t *th, ocl_device_t *dev, cl_program program,
                 const char **kernel_names,
                 const gpu_arg_type_list_t *arg_descs,
                 int num_kernels)
{
    cl_int err;
    int i;

    memset(th, 0, sizeof(*th));
    th->dev = dev;
    th->num_kernels = num_kernels;

    /* CL2.0 replacement for the deprecated clCreateCommandQueue (locked
     * decision #1 makes 2.0 the floor, so use the 2.0 API directly). */
    th->queue = clCreateCommandQueueWithProperties(dev->context,
            dev->info->device_handle, NULL, &err);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "ocl_thread_init: clCreateCommandQueue failed: %s\n",
                clGetErrorString(err));
        return -1;
    }

    th->launch = (gpu_launch_t *)calloc((size_t)num_kernels, sizeof(gpu_launch_t));
    if (!th->launch) {
        clReleaseCommandQueue(th->queue);
        th->queue = NULL;
        return -1;
    }

    for (i = 0; i < num_kernels; i++) {
        gpu_launch_init(program, kernel_names[i], &arg_descs[i],
                         &th->launch[i], dev->info->device_handle);
    }

    if (dev->pref_wg_multiple == 0 && num_kernels > 0) {
        dev->pref_wg_multiple = ocl_query_preferred_wg_multiple(
                th->launch[0].kernel_func, dev->info->device_handle);
    }

    return 0;
}

void
ocl_thread_free(ocl_thread_t *th)
{
    int i;
    for (i = 0; i < th->num_kernels; i++)
        if (th->launch[i].kernel_func)
            clReleaseKernel(th->launch[i].kernel_func);
    free(th->launch);
    th->launch = NULL;
    if (th->queue)
        clReleaseCommandQueue(th->queue);
    th->queue = NULL;
}
