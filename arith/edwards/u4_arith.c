/* Arithmetic check of microecm_edw_x4.c's FP52 Montgomery product, add and
   sub against __int128, with edge cases (n near 2^52, x, y near n, small n).
   gcc -O2 -mavx2 -mfma -DUSE_AVX2 -I../../factor/gmp-ecm u4_arith.c \
       ../../factor/gmp-ecm/microecm_edw.c */
#include <stdio.h>
#include "../../factor/gmp-ecm/microecm_edw_x4.c"

static uint64_t rs = 88172645463325252ULL;
static uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }

int main(void)
{
    uint64_t bad = 0, tests = 0, ties = 0;
    int it, l;
    for (it = 0; it < 4000000; it++) {
        ALIGNED_MEM double nA[4], rA[4], xA[4], yA[4], o1[4], o2[4], o3[4];
        uint64_t nU[4], xU[4], yU[4];
        for (l = 0; l < 4; l++) {
            uint64_t n;
            int kind = (int)(rnd() % 6);
            if (kind == 0) n = ((1ULL << 52) - 1 - 2 * (rnd() % 1000)) | 1;
            else if (kind == 1) n = (rnd() % 1000) | 1;
            else if (kind == 2) n = ((rnd() >> 12) | (1ULL << 51)) | 1;
            else n = (rnd() >> (12 + rnd() % 40)) | 1;
            if (n < 3) n = 3;
            uint64_t x = rnd() % n, y = rnd() % n;
            if ((rnd() & 7) == 0) x = n - 1 - (rnd() % 3) % n;
            if ((rnd() & 7) == 0) y = n - 1 - (rnd() % 3) % n;
            if (x >= n) x = n - 1;
            if (y >= n) y = n - 1;
            if ((rnd() & 3) == 0) {
                // a tie: x y = 2^51 mod 2^52 (x = odd 2^a, y = odd 2^(51 - a))
                int a = (int)(rnd() % 52);
                uint64_t xx = ((rnd() | 1) << a) & 0x000fffffffffffffULL;
                uint64_t yy = ((rnd() | 1) << (51 - a)) & 0x000fffffffffffffULL;
                if (xx < n && yy < n) { x = xx; y = yy; ties++; }
            }
            nU[l] = n; xU[l] = x; yU[l] = y;
            nA[l] = (double)n;
            rA[l] = (double)((0 - uedw4_inverse(n)) & 0x000fffffffffffffULL);
            xA[l] = (double)x; yA[l] = (double)y;
        }
        __m256d n = _mm256_load_pd(nA), rho = _mm256_load_pd(rA);
        __m256d x = _mm256_load_pd(xA), y = _mm256_load_pd(yA);
        _mm256_store_pd(o1, uedw4_mulredc(x, y, n, rho));
        _mm256_store_pd(o2, uedw4_addmod(x, y, n));
        _mm256_store_pd(o3, uedw4_submod(x, y, n));
        for (l = 0; l < 4; l++) {
            unsigned __int128 R = (unsigned __int128)1 << 52;
            unsigned __int128 inv = 1, e = 0;
            uint64_t n_ = nU[l];
            // expected x y / 2^52 mod n: (x y * (2^-52 mod n)) mod n
            // 2^-52 mod n by halving 52 times
            e = 1;
            int k;
            for (k = 0; k < 52; k++) e = (e & 1) ? (e + n_) / 2 : e / 2;
            uint64_t want = (uint64_t)(((unsigned __int128)xU[l] * yU[l] % n_) * e % n_);
            uint64_t got = (uint64_t)o1[l];
            uint64_t wa = (xU[l] + yU[l]) % n_;
            uint64_t ws = (xU[l] + n_ - yU[l]) % n_;
            (void)R; (void)inv;
            tests++;
            if (got != want || o1[l] < 0 || (uint64_t)o2[l] != wa || (uint64_t)o3[l] != ws ||
                o2[l] < 0 || o3[l] < 0 || o1[l] != (double)want) {
                if (bad < 10)
                    printf("bad n=%llu x=%llu y=%llu mul %.0f want %llu add %.0f sub %.0f\n",
                        (unsigned long long)n_, (unsigned long long)xU[l],
                        (unsigned long long)yU[l], o1[l], (unsigned long long)want, o2[l], o3[l]);
                bad++;
            }
        }
    }
    printf("%llu of %llu bad (%llu ties)\n", (unsigned long long)bad, (unsigned long long)tests,
        (unsigned long long)ties);
    return bad != 0;
}
