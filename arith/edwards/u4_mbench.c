/* Throughput of the x4 FP52 Montgomery product vs the scalar 64-bit one:
   K independent product chains, ns per lane-product.
   gcc -O2 -mavx2 -mfma -mbmi2 -DUSE_AVX2 -I../../factor/gmp-ecm u4_mbench.c \
       ../../factor/gmp-ecm/microecm_edw.c */
#include <stdio.h>
#include <time.h>
#include "../../factor/gmp-ecm/microecm_edw_x4.c"

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

static uint64_t smul(uint64_t x, uint64_t y, uint64_t n, uint64_t inv)
{
    unsigned __int128 p = (unsigned __int128)x * y;
    uint64_t lo = (uint64_t)p, hi = (uint64_t)(p >> 64);
    uint64_t m = lo * inv;
    uint64_t mh = (uint64_t)(((unsigned __int128)m * n) >> 64);
    uint64_t r = hi - mh;
    return (hi < mh) ? r + n : r;
}

int main(void)
{
    const int iters = 20000000;
    uint64_t n = 4503599627370449ULL;           // < 2^52
    ALIGNED_MEM double nA[4], rA[4], o[4];
    int i, l;
    for (l = 0; l < 4; l++) {
        nA[l] = (double)n;
        rA[l] = (double)((0 - uedw4_inverse(n)) & 0x000fffffffffffffULL);
    }
    __m256d N = _mm256_load_pd(nA), rho = _mm256_load_pd(rA);
    __m256d a[8];
    for (l = 0; l < 8; l++)
        a[l] = _mm256_set_pd(12345.0 + l, 67890.0 + l, 31337.0 + l, 4242.0 + l);
    __m256d b = _mm256_set1_pd(987654321.0);
    double t0 = now();
    for (i = 0; i < iters; i++) {
        a[0] = uedw4_mulredc(a[0], b, N, rho);
        a[1] = uedw4_mulredc(a[1], b, N, rho);
        a[2] = uedw4_mulredc(a[2], b, N, rho);
        a[3] = uedw4_mulredc(a[3], b, N, rho);
        a[4] = uedw4_mulredc(a[4], b, N, rho);
        a[5] = uedw4_mulredc(a[5], b, N, rho);
        a[6] = uedw4_mulredc(a[6], b, N, rho);
        a[7] = uedw4_mulredc(a[7], b, N, rho);
    }
    double t1 = now();
    __m256d s = a[0];
    for (l = 1; l < 8; l++) s = _mm256_add_pd(s, a[l]);
    _mm256_store_pd(o, s);
    printf("x4 fp52 mulredc: %.3f ns per lane-product (%.0f)\n",
        (t1 - t0) * 1e9 / ((double)iters * 32), o[0]);

    uint64_t inv = 0 - uedw4_inverse(n);        // -1/n mod 2^64 (scalar R = 2^64)
    inv = uedw4_inverse(n);
    uint64_t x[8], y = 987654321;
    for (l = 0; l < 8; l++) x[l] = 12345 + l;
    t0 = now();
    for (i = 0; i < iters; i++) {
        for (l = 0; l < 8; l++)
            x[l] = smul(x[l], y, n, inv);
    }
    t1 = now();
    uint64_t sx = 0;
    for (l = 0; l < 8; l++) sx += x[l];
    printf("scalar mulredc:  %.3f ns per product (%llu)\n",
        (t1 - t0) * 1e9 / ((double)iters * 8), (unsigned long long)sx);
    return 0;
}
