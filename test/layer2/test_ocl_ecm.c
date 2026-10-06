/*----------------------------------------------------------------------
 Layer 2 -- OpenCL batch-factoring kernels (factor/opencl_tinyecm.cl).

 Builds the production OpenCL kernels the same way gpu_cofactorization_cl.c
 does (opencl_intrinsics.cl + opencl_tinyecm.cl, -cl-std=CL2.0
 -cl-mad-enable), feeds them composites of known structure, and checks every
 factor the GPU reports:

     gbl_ecm    64-bit ECM, semiprimes of 61..62 bits
     gbl_ecm96  96-bit ECM, semiprimes of ~82..87 bits (30..34-bit p)
     gbl_pm196  96-bit P-1, with p-1 built from distinct primes < 500

 Each input is classified like the CPU ECM tests in test_ecm.c:
     pass  -- a proper divisor (1 < f < N, f | N) was reported;
     miss  -- no factor was reported (ECM is probabilistic: a small miss rate
              is allowed, not a correctness failure);
     hard  -- f > 1, f < N but f does not divide N: a real bug.
 A kernel fails on any hard failure, or if its miss rate exceeds a budget.

 If no OpenCL GPU is available the tests print a note and pass, so a build
 made with WITH_OPENCL=1 still runs on a machine without a GPU.

 Needs a build with HAVE_OCL_BATCH_FACTOR (make WITH_OPENCL=1 test) and must
 be run from the yafu source root, or with YAFU_OCL_KERNEL_DIR pointing at the
 directory holding opencl_tinyecm.cl and opencl_intrinsics.cl.
 YAFU_OCL_PLATFORM (substring of the platform name) selects a platform.
 Public domain.
----------------------------------------------------------------------*/
#include "testkit.h"

#ifdef HAVE_OCL_BATCH_FACTOR

#include "test_data.h"
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef CL_TARGET_OPENCL_VERSION
#define CL_TARGET_OPENCL_VERSION 200
#endif
#ifdef __APPLE__
#  include <OpenCL/opencl.h>
#else
#  include <CL/cl.h>
#endif

#define OCL_NTRIALS_64   16384
#define OCL_NTRIALS_96   8192
#define OCL_CURVES_64    30
#define OCL_CURVES_96    64
#define OCL_MAX_MISS_PCT 2.0

/* ------------------------------------------------------------------ *
 * One shared context/program for the whole module (the first build of
 * the kernels takes tens of seconds on some drivers).
 * ------------------------------------------------------------------ */
static struct {
    int            tried, ok;
    cl_context     ctx;
    cl_command_queue q;
    cl_program     prog;
    cl_kernel      k_ecm, k_ecm96, k_pm196;
    char           devname[128];
} G;

static char *slurp(const char *path)
{
    FILE *f = fopen(path, "rb");
    long sz;
    char *b;
    if (!f) return NULL;
    fseek(f, 0, SEEK_END); sz = ftell(f); rewind(f);
    b = (char *)malloc((size_t)sz + 1);
    if (!b || fread(b, 1, (size_t)sz, f) != (size_t)sz) { free(b); fclose(f); return NULL; }
    b[sz] = 0; fclose(f);
    return b;
}

/* returns 1 if the kernels are ready; prints the reason and returns 0 if not */
static int ocl_setup(tk_ctx *tk)
{
    cl_platform_id plats[16];
    cl_uint np = 0, i;
    cl_device_id dev = 0;
    cl_int e;
    const char *want = getenv("YAFU_OCL_PLATFORM");
    const char *dir  = getenv("YAFU_OCL_KERNEL_DIR");
    char path[1024], *s1, *s2, *inc, *body;
    const char *srcs[2];
    double t0;

    if (G.tried) return G.ok;
    G.tried = 1;
    if (!dir) dir = "factor";

    if (clGetPlatformIDs(16, plats, &np) != CL_SUCCESS || np == 0) {
        printf("        (no OpenCL platform: skipping)\n"); return 0;
    }
    for (i = 0; i < np && !dev; i++) {
        char nm[128] = {0};
        clGetPlatformInfo(plats[i], CL_PLATFORM_NAME, sizeof nm, nm, NULL);
        if (want && !strstr(nm, want)) continue;
        if (clGetDeviceIDs(plats[i], CL_DEVICE_TYPE_GPU, 1, &dev, NULL) != CL_SUCCESS) dev = 0;
    }
    if (!dev) { printf("        (no OpenCL GPU device: skipping)\n"); return 0; }
    clGetDeviceInfo(dev, CL_DEVICE_NAME, sizeof G.devname, G.devname, NULL);

    snprintf(path, sizeof path, "%s/opencl_intrinsics.cl", dir); s1 = slurp(path);
    snprintf(path, sizeof path, "%s/opencl_tinyecm.cl", dir);    s2 = slurp(path);
    if (!s1 || !s2) {
        TK_CHECKF(tk, 0, "cannot read the kernel sources from '%s' "
                  "(run from the source root or set YAFU_OCL_KERNEL_DIR)", dir);
        free(s1); free(s2); return 0;
    }
    /* same assembly as the host code: tinyecm's own #include is dropped */
    body = s2;
    inc = strstr(s2, "#include \"opencl_intrinsics.cl\"");
    if (inc) { char *nl = strchr(inc, '\n'); if (nl) body = nl + 1; }
    srcs[0] = s1; srcs[1] = body;

    G.ctx = clCreateContext(NULL, 1, &dev, NULL, NULL, &e);
    if (e != CL_SUCCESS) { TK_CHECKF(tk, 0, "clCreateContext: %d", e); return 0; }
    G.q = clCreateCommandQueueWithProperties(G.ctx, dev, NULL, &e);
    if (e != CL_SUCCESS) { TK_CHECKF(tk, 0, "clCreateCommandQueue: %d", e); return 0; }

    printf("        device: %s  (building kernels...)\n", G.devname);
    fflush(stdout);
    t0 = tk_now_sec();
    G.prog = clCreateProgramWithSource(G.ctx, 2, srcs, NULL, &e);
    if (e == CL_SUCCESS)
        e = clBuildProgram(G.prog, 1, &dev, "-cl-std=CL2.0 -cl-mad-enable", NULL, NULL);
    if (e != CL_SUCCESS) {
        char log[4096] = {0};
        clGetProgramBuildInfo(G.prog, dev, CL_PROGRAM_BUILD_LOG, sizeof log - 1, log, NULL);
        TK_CHECKF(tk, 0, "kernel build failed (%d): %s", e, log);
        return 0;
    }
    if (tk_verbose(tk)) printf("        kernel build: %.1f s\n", tk_now_sec() - t0);
    G.k_ecm   = clCreateKernel(G.prog, "gbl_ecm",   &e);
    if (e == CL_SUCCESS) G.k_ecm96 = clCreateKernel(G.prog, "gbl_ecm96", &e);
    if (e == CL_SUCCESS) G.k_pm196 = clCreateKernel(G.prog, "gbl_pm196", &e);
    if (e != CL_SUCCESS) { TK_CHECKF(tk, 0, "clCreateKernel: %d", e); return 0; }
    G.ok = 1;
    return 1;
}

/* ------------------------------------------------------------------ *
 * helpers (same setup math as gpu_cofactorization_cl.c)
 * ------------------------------------------------------------------ */
static uint32_t neg_inverse32(uint64_t a)
{
    uint32_t r = 2 + (uint32_t)a;
    r = r * (2 + (uint32_t)a * r);
    r = r * (2 + (uint32_t)a * r);
    r = r * (2 + (uint32_t)a * r);
    return r * (2 + (uint32_t)a * r);
}

/* write z as 3 little-endian 32-bit words */
static void to_w96(uint32_t *w, const mpz_t z)
{
    size_t cnt = 0, i;
    uint32_t tmp[8] = {0};
    mpz_export(tmp, &cnt, -1, 4, 0, 0, z);
    for (i = 0; i < 3; i++) w[i] = (i < cnt) ? tmp[i] : 0;
}

static void from_w96(mpz_t z, const uint32_t *w)
{
    mpz_import(z, 3, -1, 4, 0, 0, w);
}

static cl_mem mkbuf(size_t bytes, const void *src)
{
    cl_int e;
    cl_mem m = clCreateBuffer(G.ctx, CL_MEM_READ_WRITE, bytes, NULL, &e);
    if (e == CL_SUCCESS && src)
        e = clEnqueueWriteBuffer(G.q, m, CL_TRUE, 0, bytes, src, 0, NULL, NULL);
    return (e == CL_SUCCESS) ? m : NULL;
}

static int run1d(cl_kernel k, size_t n)
{
    size_t local = 64, global = (n + local - 1) / local * local;
    if (clEnqueueNDRangeKernel(G.q, k, 1, NULL, &global, &local, 0, NULL, NULL) != CL_SUCCESS)
        return -1;
    return (clFinish(G.q) == CL_SUCCESS) ? 0 : -1;
}

#define SETARG(k, i, v) clSetKernelArg((k), (cl_uint)(i), sizeof(v), &(v))

static void report(tk_ctx *tk, const char *name, long pass, long miss, long hard, long n)
{
    double pct = (double)miss * 100.0 / (double)n;
    if (tk_verbose(tk))
        printf("        %s: %ld pass, %ld miss, %ld hard / %ld\n", name, pass, miss, hard, n);
    TK_CHECKF(tk, hard == 0, "%s: %ld wrong factor(s)", name, hard);
    TK_CHECKF(tk, pct <= OCL_MAX_MISS_PCT,
              "%s: miss rate %.2f%% (%ld/%ld) exceeds %.1f%%", name, pct, miss, n, OCL_MAX_MISS_PCT);
}

/* ================================================================== *
 * gbl_ecm: 64-bit moduli
 * ================================================================== */
static void t_ocl_ecm64(tk_ctx *tk)
{
    const int N = OCL_NTRIALS_64;
    uint64_t *n, *one, *rsq, *f;
    uint32_t *rho, *sg;
    char *found;
    mpz_t z, r, m;
    cl_mem dn, drho, done, drsq, dsg, dfo;
    long pass = 0, miss = 0, hard = 0;
    tk_rng *rng = tk_rng_of(tk);
    uint32_t stg1 = 205;
    int i, c;

    if (!ocl_setup(tk)) return;
    n = malloc(8 * (size_t)N); one = malloc(8 * (size_t)N); rsq = malloc(8 * (size_t)N);
    f = calloc((size_t)N, 8); rho = malloc(4 * (size_t)N); sg = malloc(4 * (size_t)N);
    found = calloc((size_t)N, 1);
    mpz_init(z); mpz_init(r); mpz_init(m);
    for (i = 0; i < N; i++) {
        uint64_t p, q;
        n[i] = tk_gen_semiprime_u64(tk, 61 + (i & 1), &p, &q);
        rho[i] = neg_inverse32(n[i]);
        one[i] = (uint64_t)0 - n[i];  one[i] %= n[i];       /* 2^64 mod n  */
        mpz_import(m, 1, -1, 8, 0, 0, &n[i]);
        mpz_set_ui(z, 1); mpz_mul_2exp(z, z, 128);           /* 2^128 mod n */
        mpz_tdiv_r(r, z, m);
        rsq[i] = 0;
        { size_t cnt = 0; mpz_export(&rsq[i], &cnt, -1, 8, 0, 0, r); }
        sg[i] = (uint32_t)tk_rng_u64(rng);
    }
    dn = mkbuf(8 * (size_t)N, n);    drho = mkbuf(4 * (size_t)N, rho);
    done = mkbuf(8 * (size_t)N, one); drsq = mkbuf(8 * (size_t)N, rsq);
    dsg = mkbuf(4 * (size_t)N, sg);  dfo = mkbuf(8 * (size_t)N, f);
    TK_REQUIRE(tk, dn && drho && done && drsq && dsg && dfo, "buffer allocation failed");

    for (c = 0; c < OCL_CURVES_64; c++) {
        int num = N, curve = c;
        SETARG(G.k_ecm, 0, num);  SETARG(G.k_ecm, 1, dn);   SETARG(G.k_ecm, 2, drho);
        SETARG(G.k_ecm, 3, done); SETARG(G.k_ecm, 4, drsq); SETARG(G.k_ecm, 5, dsg);
        SETARG(G.k_ecm, 6, dfo);  SETARG(G.k_ecm, 7, stg1); SETARG(G.k_ecm, 8, curve);
        TK_REQUIRE(tk, run1d(G.k_ecm, (size_t)N) == 0, "gbl_ecm launch failed");
        TK_REQUIRE(tk, clEnqueueReadBuffer(G.q, dfo, CL_TRUE, 0, 8 * (size_t)N, f, 0, NULL, NULL)
                       == CL_SUCCESS, "read-back failed");
        for (i = 0; i < N; i++) {
            if (found[i]) continue;
            if (f[i] > 1 && f[i] < n[i]) {
                if (n[i] % f[i] == 0) found[i] = 1;
                else { hard++; found[i] = 2;
                       TK_CHECKF(tk, 0, "gbl_ecm: n=%llu got non-divisor %llu",
                                 (unsigned long long)n[i], (unsigned long long)f[i]); }
            }
        }
    }
    for (i = 0; i < N; i++) { if (found[i] == 1) pass++; else if (!found[i]) miss++; }
    report(tk, "gbl_ecm   (64-bit)", pass, miss, hard, N);

    clReleaseMemObject(dn); clReleaseMemObject(drho); clReleaseMemObject(done);
    clReleaseMemObject(drsq); clReleaseMemObject(dsg); clReleaseMemObject(dfo);
    mpz_clear(z); mpz_clear(r); mpz_clear(m);
    free(n); free(one); free(rsq); free(f); free(rho); free(sg); free(found);
}

/* ================================================================== *
 * gbl_ecm96 / gbl_pm196: 96-bit moduli.  mode 0 = ECM, 1 = P-1.
 * ================================================================== */
static void run96(tk_ctx *tk, int pm1)
{
    const int N = OCL_NTRIALS_96;
    uint32_t *n, *one, *rsq, *f, *rho, *sg;
    char *found;
    mpz_t zn, zp, zq, zr, zf;
    cl_mem dn, drho, done, drsq, dsg, dfo;
    long pass = 0, miss = 0, hard = 0;
    tk_rng *rng = tk_rng_of(tk);
    cl_kernel k = pm1 ? G.k_pm196 : G.k_ecm96;
    int i, c, ncur = pm1 ? 1 : OCL_CURVES_96;

    if (!ocl_setup(tk)) return;
    n = calloc(3 * (size_t)N, 4); one = calloc(3 * (size_t)N, 4); rsq = calloc(3 * (size_t)N, 4);
    f = calloc(3 * (size_t)N, 4); rho = calloc((size_t)N, 4); sg = calloc((size_t)N, 4);
    found = calloc((size_t)N, 1);
    mpz_inits(zn, zp, zq, zr, zf, NULL);

    for (i = 0; i < N; i++) {
        uint64_t p, q;
        if (!pm1) {
            p = tk_gen_prime_u64(tk, 30 + (int)tk_rng_range(rng, 5));
        } else {   /* p-1 = 2 * distinct odd primes < 500, p of 34..40 bits */
            for (;;) {
                uint64_t m = 2, pr;
                int bits = 2;
                while (bits < 34) {
                    pr = 3 + 2 * tk_rng_range(rng, 248);
                    if (tk_is_prime_u64(pr) && m % pr) {
                        m *= pr; bits = 0; { uint64_t t = m; while (t) { bits++; t >>= 1; } }
                    }
                }
                if (bits <= 40 && tk_is_prime_u64(m + 1)) { p = m + 1; break; }
            }
        }
        q = tk_gen_prime_u64(tk, 52 + (int)tk_rng_range(rng, 2));
        mpz_set_ui(zp, (unsigned long)p); mpz_set_ui(zq, (unsigned long)q);
        mpz_mul(zn, zp, zq);                          /* < 2^96 by construction */
        to_w96(&n[3 * i], zn);
        rho[i] = neg_inverse32(n[3 * i]);
        mpz_set_ui(zr, 1); mpz_mul_2exp(zr, zr, 96); mpz_sub(zr, zr, zn); mpz_tdiv_r(zr, zr, zn);
        to_w96(&one[3 * i], zr);                      /* 2^96 mod n  */
        mpz_set_ui(zr, 1); mpz_mul_2exp(zr, zr, 192); mpz_tdiv_r(zr, zr, zn);
        to_w96(&rsq[3 * i], zr);                      /* 2^192 mod n */
        sg[i] = (uint32_t)tk_rng_u64(rng);
    }
    dn = mkbuf(12 * (size_t)N, n);   drho = mkbuf(4 * (size_t)N, rho);
    done = mkbuf(12 * (size_t)N, one); drsq = mkbuf(12 * (size_t)N, rsq);
    dsg = mkbuf(4 * (size_t)N, sg);  dfo = mkbuf(12 * (size_t)N, f);
    TK_REQUIRE(tk, dn && drho && done && drsq && dsg && dfo, "buffer allocation failed");

    for (c = 0; c < ncur; c++) {
        int num = N, curve = c, a = 0;
        uint32_t b1 = pm1 ? 500 : 205, b2 = b1 * 50;
        SETARG(k, a, num); a++;  SETARG(k, a, dn); a++;  SETARG(k, a, drho); a++;
        SETARG(k, a, done); a++;
        if (!pm1) { SETARG(k, a, drsq); a++; SETARG(k, a, dsg); a++; }
        SETARG(k, a, dfo); a++;  SETARG(k, a, b1); a++;  SETARG(k, a, b2); a++;
        if (!pm1) { SETARG(k, a, curve); a++; }
        TK_REQUIRE(tk, run1d(k, (size_t)N) == 0, "96-bit kernel launch failed");
        TK_REQUIRE(tk, clEnqueueReadBuffer(G.q, dfo, CL_TRUE, 0, 12 * (size_t)N, f, 0, NULL, NULL)
                       == CL_SUCCESS, "read-back failed");
        for (i = 0; i < N; i++) {
            if (found[i]) continue;
            from_w96(zf, &f[3 * i]); from_w96(zn, &n[3 * i]);
            if (mpz_cmp_ui(zf, 1) > 0 && mpz_cmp(zf, zn) < 0) {
                if (mpz_divisible_p(zn, zf)) found[i] = 1;
                else { hard++; found[i] = 2;
                       TK_CHECKF(tk, 0, "%s: reported a non-divisor of N",
                                 pm1 ? "gbl_pm196" : "gbl_ecm96"); }
            }
        }
    }
    for (i = 0; i < N; i++) { if (found[i] == 1) pass++; else if (!found[i]) miss++; }
    report(tk, pm1 ? "gbl_pm196  (96-bit P-1)" : "gbl_ecm96 (96-bit ECM)", pass, miss, hard, N);

    clReleaseMemObject(dn); clReleaseMemObject(drho); clReleaseMemObject(done);
    clReleaseMemObject(drsq); clReleaseMemObject(dsg); clReleaseMemObject(dfo);
    mpz_clears(zn, zp, zq, zr, zf, NULL);
    free(n); free(one); free(rsq); free(f); free(rho); free(sg); free(found);
}

static void t_ocl_ecm96(tk_ctx *tk) { run96(tk, 0); }
static void t_ocl_pm196(tk_ctx *tk) { run96(tk, 1); }

static const tk_test tk__ocl_ecm_tests[] = {
    { "gbl_ecm",   t_ocl_ecm64, "slow ecm ocl gpu" },
    { "gbl_ecm96", t_ocl_ecm96, "slow ecm ocl gpu" },
    { "gbl_pm196", t_ocl_pm196, "slow pm1 ocl gpu" }
};
const tk_module tk_module_ocl_ecm = {
    "ocl_ecm",
    "OpenCL batch-factoring kernels (gbl_ecm, gbl_ecm96, gbl_pm196)",
    tk__ocl_ecm_tests,
    (int)(sizeof tk__ocl_ecm_tests / sizeof tk__ocl_ecm_tests[0])
};

#endif /* HAVE_OCL_BATCH_FACTOR */
