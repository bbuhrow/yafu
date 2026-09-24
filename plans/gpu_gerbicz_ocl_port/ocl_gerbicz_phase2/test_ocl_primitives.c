/* test_ocl_primitives.c -- Phase 2 tasks 2.6 (tests) and 2.7 (benchmarks).
 *
 * CPU reference implementations are deliberately simple, not fast --
 * they exist only to check the GPU output, per the plan's own
 * instruction for Phase 0's harness ("small, deliberately simple").
 */
/* _POSIX_C_SOURCE for clock_gettime/CLOCK_MONOTONIC under strict -std=c99 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "ocl_primitives.h"

static int g_fail_count = 0;
static int g_pass_count = 0;

static double
now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* Simple deterministic PRNG (xorshift64) so runs are reproducible
 * without pulling in a dependency. */
static uint64_t g_rng_state = 0x9E3779B97F4A7C15ULL;
static uint64_t
rng_next(void)
{
    uint64_t x = g_rng_state;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    g_rng_state = x;
    return x;
}

/* ---------------------------------------------------------------------
 * CPU references
 * --------------------------------------------------------------------- */
static uint32_t
cpu_reduce_max_u32(const uint32_t *a, uint32_t n)
{
    uint32_t m = a[0];
    uint32_t i;
    for (i = 1; i < n; i++)
        if (a[i] > m) m = a[i];
    return m;
}

static void
cpu_scan_exclusive_u32(const uint32_t *in, uint32_t *out, uint32_t n)
{
    uint32_t acc = 0, i;
    for (i = 0; i < n; i++) {
        out[i] = acc;
        acc += in[i];
    }
}

static int
cmp_u64(const void *pa, const void *pb)
{
    uint64_t a = *(const uint64_t *)pa, b = *(const uint64_t *)pb;
    return (a > b) - (a < b);
}

static void
cpu_sort_u64(const uint64_t *in, uint64_t *out, uint32_t n)
{
    memcpy(out, in, (size_t)n * sizeof(uint64_t));
    qsort(out, n, sizeof(uint64_t), cmp_u64);
}

/* ---------------------------------------------------------------------
 * Test drivers
 * --------------------------------------------------------------------- */
static void
report(const char *name, uint32_t n, int ok, const char *detail)
{
    if (ok) {
        g_pass_count++;
        printf("  [PASS] %-24s n=%-9u %s\n", name, n, detail ? detail : "");
    } else {
        g_fail_count++;
        printf("  [FAIL] %-24s n=%-9u %s\n", name, n, detail ? detail : "");
    }
}

static void
test_fill(ocl_primitives_t *p, uint32_t n)
{
    cl_int err;
    cl_mem buf = clCreateBuffer(p->gd->dev.context, CL_MEM_READ_WRITE,
            (size_t)(n ? n : 1) * sizeof(uint32_t), NULL, &err);
    uint32_t *host = (uint32_t *)malloc((size_t)(n ? n : 1) * sizeof(uint32_t));
    uint32_t i;
    int ok = 1;
    char detail[64] = "";

    if (ocl_fill_u32(p, buf, 0, n, 0xDEADBEEFu) != 0) ok = 0;
    if (ok && n > 0) {
        clEnqueueReadBuffer(p->th.queue, buf, CL_TRUE, 0, (size_t)n * sizeof(uint32_t), host, 0, NULL, NULL);
        for (i = 0; i < n; i++) {
            if (host[i] != 0xDEADBEEFu) { ok = 0; snprintf(detail, sizeof(detail), "mismatch at %u", i); break; }
        }
    }
    report("fill_u32", n, ok, detail);
    free(host);
    clReleaseMemObject(buf);
}

static void
test_reduce_max(ocl_primitives_t *p, uint32_t n)
{
    uint32_t *host = (uint32_t *)malloc((size_t)n * sizeof(uint32_t));
    uint32_t expect, got = 0;
    cl_int err;
    cl_mem buf;
    uint32_t i;
    int ok;
    char detail[64];

    for (i = 0; i < n; i++)
        host[i] = (uint32_t)(rng_next() & 0xFFFFFFu);
    expect = cpu_reduce_max_u32(host, n);

    buf = clCreateBuffer(p->gd->dev.context, CL_MEM_READ_WRITE,
            (size_t)n * sizeof(uint32_t), NULL, &err);
    clEnqueueWriteBuffer(p->th.queue, buf, CL_TRUE, 0, (size_t)n * sizeof(uint32_t), host, 0, NULL, NULL);

    ok = (ocl_reduce_max_u32(p, buf, 0, n, &got) == 0) && (got == expect);
    snprintf(detail, sizeof(detail), "expect=%u got=%u", expect, got);
    report("reduce_max_u32", n, ok, detail);

    free(host);
    clReleaseMemObject(buf);
}

static void
test_scan(ocl_primitives_t *p, uint32_t n)
{
    uint32_t *host = (uint32_t *)malloc((size_t)n * sizeof(uint32_t));
    uint32_t *expect = (uint32_t *)malloc((size_t)n * sizeof(uint32_t));
    uint32_t *got = (uint32_t *)malloc((size_t)n * sizeof(uint32_t));
    cl_int err;
    cl_mem buf;
    uint32_t i;
    int ok = 1;
    char detail[64] = "";

    for (i = 0; i < n; i++)
        host[i] = (uint32_t)(rng_next() & 0xFFu); /* small values so sums don't overflow at 4M scale */
    cpu_scan_exclusive_u32(host, expect, n);

    buf = clCreateBuffer(p->gd->dev.context, CL_MEM_READ_WRITE,
            (size_t)n * sizeof(uint32_t), NULL, &err);
    clEnqueueWriteBuffer(p->th.queue, buf, CL_TRUE, 0, (size_t)n * sizeof(uint32_t), host, 0, NULL, NULL);

    if (ocl_scan_exclusive_u32(p, buf, 0, n) != 0) ok = 0;
    if (ok) {
        clEnqueueReadBuffer(p->th.queue, buf, CL_TRUE, 0, (size_t)n * sizeof(uint32_t), got, 0, NULL, NULL);
        for (i = 0; i < n; i++) {
            if (got[i] != expect[i]) {
                ok = 0;
                snprintf(detail, sizeof(detail), "mismatch at %u: got=%u want=%u", i, got[i], expect[i]);
                break;
            }
        }
    }
    report("scan_exclusive_u32", n, ok, detail);

    free(host); free(expect); free(got);
    clReleaseMemObject(buf);
}

static void
test_radix_sort(ocl_primitives_t *p, uint32_t n, int use_adversarial)
{
    uint64_t *host = (uint64_t *)malloc((size_t)n * sizeof(uint64_t));
    uint64_t *expect = (uint64_t *)malloc((size_t)n * sizeof(uint64_t));
    uint64_t *got = (uint64_t *)malloc((size_t)n * sizeof(uint64_t));
    cl_int err;
    cl_mem in_buf, out_buf;
    uint32_t i;
    int ok = 1;
    char detail[80] = "";
    char name[40];

    for (i = 0; i < n; i++) {
        if (use_adversarial) {
            /* Force lots of duplicate keys and top-bit-set values (the
             * "signed key, zero-extended storage" case from the plan's
             * Section 3 fact) to stress the sort beyond uniform random
             * data. */
            uint64_t r = rng_next();
            if ((r & 3) == 0)
                host[i] = 0xFFFFFFFFu; /* zero-extended -1 (int32) pattern, per the Fact */
            else if ((r & 3) == 1)
                host[i] = 0;
            else
                host[i] = r % 17; /* heavy duplication */
        } else {
            host[i] = rng_next();
        }
    }
    cpu_sort_u64(host, expect, n);

    in_buf = clCreateBuffer(p->gd->dev.context, CL_MEM_READ_ONLY,
            (size_t)n * sizeof(uint64_t), NULL, &err);
    out_buf = clCreateBuffer(p->gd->dev.context, CL_MEM_READ_WRITE,
            (size_t)n * sizeof(uint64_t), NULL, &err);
    clEnqueueWriteBuffer(p->th.queue, in_buf, CL_TRUE, 0, (size_t)n * sizeof(uint64_t), host, 0, NULL, NULL);

    if (ocl_radix_sort_u64(p, in_buf, 0, out_buf, 0, n) != 0) ok = 0;
    if (ok) {
        clEnqueueReadBuffer(p->th.queue, out_buf, CL_TRUE, 0, (size_t)n * sizeof(uint64_t), got, 0, NULL, NULL);
        for (i = 0; i < n; i++) {
            if (got[i] != expect[i]) {
                ok = 0;
                snprintf(detail, sizeof(detail), "mismatch at %u: got=%llu want=%llu",
                         i, (unsigned long long)got[i], (unsigned long long)expect[i]);
                break;
            }
        }
    }
    snprintf(name, sizeof(name), "radix_sort_u64%s", use_adversarial ? "(adv)" : "");
    report(name, n, ok, detail);

    free(host); free(expect); free(got);
    clReleaseMemObject(in_buf);
    clReleaseMemObject(out_buf);
}

static void
bench_radix_sort(ocl_primitives_t *p, uint32_t n)
{
    uint64_t *host = (uint64_t *)malloc((size_t)n * sizeof(uint64_t));
    cl_int err;
    cl_mem in_buf, out_buf;
    uint32_t i;
    double t0, t1;

    for (i = 0; i < n; i++)
        host[i] = rng_next();

    in_buf = clCreateBuffer(p->gd->dev.context, CL_MEM_READ_ONLY,
            (size_t)n * sizeof(uint64_t), NULL, &err);
    out_buf = clCreateBuffer(p->gd->dev.context, CL_MEM_READ_WRITE,
            (size_t)n * sizeof(uint64_t), NULL, &err);
    clEnqueueWriteBuffer(p->th.queue, in_buf, CL_TRUE, 0, (size_t)n * sizeof(uint64_t), host, 0, NULL, NULL);
    clFinish(p->th.queue);

    t0 = now_sec();
    ocl_radix_sort_u64(p, in_buf, 0, out_buf, 0, n);
    t1 = now_sec();

    printf("  [BENCH] radix_sort_u64 n=%-9u %.3f s  (%.1f M keys/s)\n",
           n, t1 - t0, (double)n / (t1 - t0) / 1e6);

    free(host);
    clReleaseMemObject(in_buf);
    clReleaseMemObject(out_buf);
}

static void
bench_scan(ocl_primitives_t *p, uint32_t n)
{
    uint32_t *host = (uint32_t *)malloc((size_t)n * sizeof(uint32_t));
    cl_int err;
    cl_mem buf;
    uint32_t i;
    double t0, t1;

    for (i = 0; i < n; i++)
        host[i] = (uint32_t)(rng_next() & 0xFFu);

    buf = clCreateBuffer(p->gd->dev.context, CL_MEM_READ_WRITE,
            (size_t)n * sizeof(uint32_t), NULL, &err);
    clEnqueueWriteBuffer(p->th.queue, buf, CL_TRUE, 0, (size_t)n * sizeof(uint32_t), host, 0, NULL, NULL);
    clFinish(p->th.queue);

    t0 = now_sec();
    ocl_scan_exclusive_u32(p, buf, 0, n);
    clFinish(p->th.queue);
    t1 = now_sec();

    printf("  [BENCH] scan_exclusive_u32 n=%-9u %.3f s  (%.1f M elems/s)\n",
           n, t1 - t0, (double)n / (t1 - t0) / 1e6);

    free(host);
    clReleaseMemObject(buf);
}

int
main(void)
{
    gpu_config_t config;
    ocl_gerbicz_device_t gd;
    ocl_primitives_t prim;
    uint32_t sizes[] = { 1, 17, 255, 256, 4096, 4097, 8193, 100000 };
    size_t i, num_sizes = sizeof(sizes) / sizeof(sizes[0]);

    gpu_init(&config);
    if (config.num_gpu == 0) {
        fprintf(stderr, "no OpenCL devices found\n");
        return 1;
    }
    printf("using device 0: %s (OpenCL %d.%d, cl_c_std pending query)\n",
           config.info[0].name, config.info[0].compute_version_major,
           config.info[0].compute_version_minor);

    if (ocl_gerbicz_device_init(&gd, &config.info[0]) != 0)
        return 1;
    printf("cl_c_std=%s has_subgroups=%d local_mem=%llu max_wg=%zu\n",
           gd.dev.cl_c_std, gd.dev.has_subgroups,
           (unsigned long long)gd.dev.local_mem_size, gd.dev.max_work_group_size);

    if (ocl_primitives_init(&prim, &gd, "ocl_primitives_kernels.cl", ".") != 0) {
        ocl_gerbicz_device_free(&gd);
        return 1;
    }

    printf("\n== fill ==\n");
    test_fill(&prim, 0);
    test_fill(&prim, 1);
    test_fill(&prim, 4096);
    test_fill(&prim, 4097);

    printf("\n== reduce_max ==\n");
    for (i = 0; i < num_sizes; i++)
        test_reduce_max(&prim, sizes[i]);
    test_reduce_max(&prim, 16384); /* NUM_BUCKETS, the actual collision_engine.cu use case */

    printf("\n== scan_exclusive ==\n");
    for (i = 0; i < num_sizes; i++)
        test_scan(&prim, sizes[i]);
    test_scan(&prim, 1048640);  /* MAX_DSIZE */
    test_scan(&prim, 4194305);  /* VALUE_MATCH_CAP + 1, the actual ~4M CUB call */

    printf("\n== radix_sort_u64 ==\n");
    for (i = 0; i < num_sizes; i++) {
        test_radix_sort(&prim, sizes[i], 0);
        test_radix_sort(&prim, sizes[i], 1);
    }
    test_radix_sort(&prim, 4194304, 0); /* CANDIDATE_CAP */
    test_radix_sort(&prim, 4194304, 1); /* CANDIDATE_CAP, adversarial */

    printf("\n== benchmarks (order-of-magnitude only; see README for device caveat) ==\n");
    bench_scan(&prim, 4194305);
    bench_radix_sort(&prim, 4194304);

    printf("\n%d passed, %d failed\n", g_pass_count, g_fail_count);

    ocl_primitives_free(&prim);
    ocl_gerbicz_device_free(&gd);
    return g_fail_count ? 1 : 0;
}
