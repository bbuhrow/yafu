/*----------------------------------------------------------------------
 Layer 1 -- arith/limb checked against GMP (links limb1.c, limb2.c, GMP).

 The header-only kernels have self-contained tests in Layer 0
 (layer0/test_limb_inline.c).  This module covers the out-of-line library
 routines and re-checks the main kernels with an independent oracle:

   limb2   two-limb add/sub/shift/mul/sqr, uint128_div / mod_128, the
           GMP-free monty128_init / to_monty128, mulmod128[n] / sqrmod128[n]
           / addmod128 / submod128 / dblmod128 / chkmod128 with aliased
           operands and moduli pinned near 2^128, the
           104/128-bit scans, gcd128 / bin_gcd128, and (AVX-512) every
           104-bit x8 kernel
   limb1   positive-inverse REDC, spPRP2, the sp* helpers that
           test_sp_arith.c does not reach (spAdd/spAdd3/spSub/spSub3/
           spMulAdd/u64div/spBits/ndigits_1/modinv_1b), the scalar gcds,
           jacobi_1, ShanksTonelli_1, and (AVX-512) the 52-bit x8 kernels

 GMP is the oracle throughout (the library itself does not use it).
 Public domain.
----------------------------------------------------------------------*/
#include "testkit.h"
#include "limb1.h"
#include "limb2.h"
#include <gmp.h>
#include <string.h>

/* ---- helpers: limb arrays <-> mpz ---- */
static void set_limbs(mpz_t z, const uint64_t *x, int n)
{
    mpz_import(z, (size_t)n, -1, sizeof(uint64_t), 0, 0, x);
}

static void get_limbs(uint64_t *x, int n, const mpz_t z)
{
    size_t cnt = 0, i;
    uint64_t tmp[8] = { 0 };
    mpz_export(tmp, &cnt, -1, sizeof(uint64_t), 0, 0, z);
    for (i = 0; i < (size_t)n; i++)
        x[i] = (i < cnt) ? tmp[i] : 0;
}

/* random two-limb value with a random bit length (so short operands, zero
   high limbs and full-width values all show up) */
static void rand128(tk_rng *r, uint64_t *x)
{
    int bits = (int)tk_rng_range(r, 129);
    x[0] = tk_rng_u64(r);
    x[1] = tk_rng_u64(r);
    if (bits == 0) { x[0] = x[1] = 0; }
    else if (bits <= 64) { x[1] = 0; x[0] >>= (64 - bits); }
    else if (bits < 128) { x[1] >>= (128 - bits); }
}

/* ================================================================== */
static void t_addsub_shift(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    mpz_t a, b, c, m;
    int k;
    mpz_inits(a, b, c, m, NULL);
    mpz_setbit(m, 128);                         /* 2^128 */
    for (k = 0; k < 100000; k++) {
        uint64_t x[2], y[2], z[2], e[2];
        int s = (int)tk_rng_range(r, 128);
        mp_carry_t cy;
        rand128(r, x); rand128(r, y);
        set_limbs(a, x, 2); set_limbs(b, y, 2);

        cy = add128(x, y, z);
        mpz_add(c, a, b);
        TK_EQ_U64(tk, cy, mpz_tstbit(c, 128));
        mpz_fdiv_r_2exp(c, c, 128); get_limbs(e, 2, c);
        if (!TK_EQ_U64(tk, z[0], e[0]) || !TK_EQ_U64(tk, z[1], e[1])) break;

        cy = sub128(x, y, z);
        mpz_sub(c, a, b);
        TK_EQ_U64(tk, cy, mpz_sgn(c) < 0);
        mpz_fdiv_r_2exp(c, c, 128); get_limbs(e, 2, c);
        if (!TK_EQ_U64(tk, z[0], e[0]) || !TK_EQ_U64(tk, z[1], e[1])) break;

        TK_CHECK(tk, cmp128(x, y) == (mpz_cmp(a, b) > 0) - (mpz_cmp(a, b) < 0));

        shl128(x, s, z);
        mpz_mul_2exp(c, a, (mp_bitcnt_t)s); mpz_fdiv_r_2exp(c, c, 128); get_limbs(e, 2, c);
        if (!TK_EQ_U64(tk, z[0], e[0]) || !TK_EQ_U64(tk, z[1], e[1])) break;

        shr128(x, s, z);
        mpz_fdiv_q_2exp(c, a, (mp_bitcnt_t)s); get_limbs(e, 2, c);
        if (!TK_EQ_U64(tk, z[0], e[0]) || !TK_EQ_U64(tk, z[1], e[1])) break;
    }
    mpz_clears(a, b, c, m, NULL);
}

static void t_mul_sqr(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    mpz_t a, b, c;
    int k, i;
    mpz_inits(a, b, c, NULL);
    for (k = 0; k < 100000; k++) {
        uint64_t x[2], y[2], z[4], e[4];
        rand128(r, x); rand128(r, y);
        if ((k & 15) == 0) { x[0] = x[1] = ~(uint64_t)0; }
        if ((k & 31) == 0) { y[0] = y[1] = ~(uint64_t)0; }
        set_limbs(a, x, 2); set_limbs(b, y, 2);

        mul128(x, y, z);
        mpz_mul(c, a, b); get_limbs(e, 4, c);
        for (i = 0; i < 4; i++) if (!TK_EQ_U64(tk, z[i], e[i])) break;

        sqr128(x, z);
        mpz_mul(c, a, a); get_limbs(e, 4, c);
        for (i = 0; i < 4; i++) if (!TK_EQ_U64(tk, z[i], e[i])) break;

        mul128x64(x, y[0], z);
        mpz_set_ui(c, 0); set_limbs(c, y, 1); mpz_mul(c, a, c); get_limbs(e, 3, c);
        for (i = 0; i < 3; i++) if (!TK_EQ_U64(tk, z[i], e[i])) break;
    }
    if (tk_bench_enabled(tk)) {
        uint64_t x[2], y[2], z[4];
        rand128(r, x); rand128(r, y);
        TK_BENCH(tk, "mul128 (128x128->256)", 1, 0.0,
                 { x[0] += 0x9E37; mul128(x, y, z); tk_sink_u64 += z[0] ^ z[3]; });
        TK_BENCH(tk, "sqr128 (128->256)", 1, 0.0,
                 { x[0] += 0x9E37; sqr128(x, z); tk_sink_u64 += z[0] ^ z[3]; });
    }
    mpz_clears(a, b, c, NULL);
}

static void t_div128(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    mpz_t a, b, q, rem;
    int k;
    uint64_t z[2] = { 0, 0 }, one[2] = { 1, 0 };
    mpz_inits(a, b, q, rem, NULL);

    TK_CHECK(tk, uint128_div(one, z, NULL, NULL) == -1);    /* divide by zero */

    for (k = 0; k < 200000; k++) {
        uint64_t n[2], d[2], qq[2], rr[2], eq[2], er[2];
        rand128(r, n); rand128(r, d);
        /* stress the quotient-estimate corrections: divisors with a top limb
           of all ones / a low limb of all ones, dividends just above them */
        if ((k % 7) == 0)  { d[1] |= 0x8000000000000000ULL; d[0] = ~(uint64_t)0; }
        if ((k % 11) == 0) { n[1] = d[1]; n[0] = tk_rng_u64(r); }
        if ((k % 13) == 0) { n[0] = n[1] = ~(uint64_t)0; }
        if ((d[0] | d[1]) == 0) d[0] = 1;
        set_limbs(a, n, 2); set_limbs(b, d, 2);
        mpz_tdiv_qr(q, rem, a, b);
        get_limbs(eq, 2, q); get_limbs(er, 2, rem);

        TK_CHECK(tk, uint128_div(n, d, qq, rr) == 0);
        if (!TK_EQ_U64(tk, qq[0], eq[0]) || !TK_EQ_U64(tk, qq[1], eq[1])) break;
        if (!TK_EQ_U64(tk, rr[0], er[0]) || !TK_EQ_U64(tk, rr[1], er[1])) break;

        mod_128(n, d, rr);
        if (!TK_EQ_U64(tk, rr[0], er[0]) || !TK_EQ_U64(tk, rr[1], er[1])) break;
    }
    if (tk_bench_enabled(tk)) {
        uint64_t n[2], d[2], qq[2], rr[2];
        n[0] = tk_rng_u64(r); n[1] = tk_rng_u64(r);
        d[0] = tk_rng_u64(r); d[1] = tk_rng_u64(r) >> 20;
        TK_BENCH(tk, "uint128_div (2-limb divisor)", 1, 0.0,
                 { n[0] += 0x9E37; uint128_div(n, d, qq, rr); tk_sink_u64 += qq[0] ^ rr[0]; });
    }
    mpz_clears(a, b, q, rem, NULL);
}

/* 128-bit Montgomery setup without GMP, checked with GMP */
static void t_monty128_setup(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    mpz_t n, x, e, R;
    int k;
    mpz_inits(n, x, e, R, NULL);
    mpz_setbit(R, 128);
    for (k = 0; k < 50000; k++) {
        monty128_t m;
        uint64_t nn[2], xx[2], ee[2], w[2];
        rand128(r, nn); rand128(r, xx);
        nn[0] |= 1;
        if (nn[1] == 0 && nn[0] < 3) nn[0] = 3;
        set_limbs(n, nn, 2); set_limbs(x, xx, 2);

        monty128_init(&m, nn);
        TK_EQ_U64(tk, m.rho * nn[0], ~(uint64_t)0);            /* rho = -n^-1 */
        mpz_mod(e, R, n); get_limbs(ee, 2, e);
        if (!TK_EQ_U64(tk, m.one[0], ee[0]) || !TK_EQ_U64(tk, m.one[1], ee[1])) break;

        /* to_monty128 on an arbitrary (possibly unreduced) x */
        w[0] = xx[0]; w[1] = xx[1];
        to_monty128(&m, w);
        mpz_mul(e, x, R); mpz_mod(e, e, n); get_limbs(ee, 2, e);
        if (!TK_EQ_U64(tk, w[0], ee[0]) || !TK_EQ_U64(tk, w[1], ee[1])) break;

        /* and back out through a Montgomery multiply by 1: w * 1 * R^-1 = x mod n */
        {
            uint64_t one[2] = { 1, 0 }, back[2];
            mulmod128(w, one, back, &m);
            chkmod128(back, nn);
            mpz_mod(e, x, n); get_limbs(ee, 2, e);
            if (!TK_EQ_U64(tk, back[0], ee[0]) || !TK_EQ_U64(tk, back[1], ee[1])) break;
        }
    }
    mpz_clears(n, x, e, R, NULL);
}

/* positive-inverse 64-bit Montgomery kernels, addmod/submod, mfma64 */
static void t_redc64_pos(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    mpz_t n, a, b, c, e, Rinv, R;
    int k;
    mpz_inits(n, a, b, c, e, Rinv, R, NULL);
    mpz_setbit(R, 64);
    for (k = 0; k < 100000; k++) {
        int bits = 2 + (int)tk_rng_range(r, 63);
        uint64_t nn = tk_rng_nbits(r, bits) | 1, x, y, z, inv, got;
        if (nn < 3) nn = 3;
        x = tk_rng_u64(r) % nn; y = tk_rng_u64(r) % nn; z = tk_rng_u64(r) % nn;
        inv = multiplicative_inverse(nn);
        TK_EQ_U64(tk, inv * nn, 1);
        TK_EQ_U64(tk, (multiplicative_inverse52(nn) * nn) & 0x000fffffffffffffULL, 1);

        set_limbs(n, &nn, 1); set_limbs(a, &x, 1); set_limbs(b, &y, 1); set_limbs(c, &z, 1);
        mpz_invert(Rinv, R, n);

        mpz_mul(e, a, b); mpz_mul(e, e, Rinv); mpz_mod(e, e, n);
        get_limbs(&got, 1, e);
        if (!TK_EQ_U64(tk, mulredc_pos(x, y, nn, inv), got)) break;
        if (!TK_EQ_U64(tk, mulredc_pos_alt(x, y, nn, inv), got)) break;
        if (!TK_EQ_U64(tk, mulredc(x, y, nn, (uint64_t)0 - inv), got)) break;

        mpz_mul(e, a, a); mpz_mul(e, e, Rinv); mpz_mod(e, e, n);
        get_limbs(&got, 1, e);
        if (!TK_EQ_U64(tk, sqrredc_pos(x, nn, inv), got)) break;
        if (!TK_EQ_U64(tk, sqrredc(x, nn, (uint64_t)0 - inv), got)) break;

        /* mfma64(x, y, c) = (x*y + c*R) * R^-1 = x*y*R^-1 + c  (mod n) */
        mpz_mul(e, a, b); mpz_mul(e, e, Rinv); mpz_add(e, e, c); mpz_mod(e, e, n);
        get_limbs(&got, 1, e);
        if (!TK_EQ_U64(tk, mfma64(x, y, z, nn, inv), got)) break;

        mpz_add(e, a, b); mpz_mod(e, e, n); get_limbs(&got, 1, e);
        if (!TK_EQ_U64(tk, addmod(x, y, nn), got)) break;
        mpz_sub(e, a, b); mpz_mod(e, e, n); get_limbs(&got, 1, e);
        if (!TK_EQ_U64(tk, submod(x, y, nn), got)) break;

        if (nn <= 0xffffffffULL) {
            TK_EQ_U64(tk, addmod32((uint32_t)x, (uint32_t)y, (uint32_t)nn), (x + y) % nn);
            TK_EQ_U64(tk, submod32((uint32_t)x, (uint32_t)y, (uint32_t)nn), (x + nn - y) % nn);
        }

        /* gcd */
        {
            uint64_t u = tk_rng_u64(r) | 1, v = tk_rng_u64(r) >> tk_rng_range(r, 64);
            set_limbs(a, &u, 1); set_limbs(b, &v, 1);
            mpz_gcd(e, a, b); get_limbs(&got, 1, e);
            if (!TK_EQ_U64(tk, bingcd64(u, v), got)) break;
        }
    }
    if (tk_bench_enabled(tk)) {
        uint64_t nn = tk_rng_u64(r) | 0x8000000000000001ULL, x = tk_rng_u64(r) % nn;
        uint64_t inv = multiplicative_inverse(nn);
        TK_BENCH(tk, "mulredc_pos latency", 1, 0.0,
                 { x = mulredc_pos(x, x, nn, inv); tk_sink_u64 += x; });
        TK_BENCH(tk, "mulredc (neg inverse) latency", 1, 0.0,
                 { x = mulredc(x, x, nn, (uint64_t)0 - inv); tk_sink_u64 += x; });
    }
    mpz_clears(n, a, b, c, e, Rinv, R, NULL);
}

static void t_spPRP2(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    mpz_t n, e, two;
    int k;
    mpz_inits(n, e, two, NULL);
    mpz_set_ui(two, 2);
    for (k = 0; k < 20000; k++) {
        int bits = 33 + (int)tk_rng_range(r, 32);        /* spPRP2 wants p > 2^32 */
        uint64_t p = tk_rng_nbits(r, bits) | 1, got;
        if ((k & 3) == 0) {                              /* mix in real primes */
            set_limbs(n, &p, 1); mpz_nextprime(n, n);
            if (mpz_sizeinbase(n, 2) > 64) continue;
            get_limbs(&p, 1, n);
        }
        set_limbs(n, &p, 1);
        mpz_sub_ui(e, n, 1);
        mpz_powm(e, two, e, n);
        get_limbs(&got, 1, e);
        if (!TK_EQ_U64(tk, spPRP2(p), got)) break;
    }
    mpz_clears(n, e, two, NULL);
}

#if MP_USE_AVX512F && !defined(LIMB1_NO_VECTOR)
/* shared 8-lane 52-bit Montgomery kernels vs a scalar reference */
static void t_redc52_x8(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    const uint64_t M52 = 0x000fffffffffffffULL;
    int k, i;
    for (k = 0; k < 20000; k++) {
        uint64_t nn[8], x[8], y[8], inv[8], out[8];
        __m512i vn, vx, vy, vinv, v;
        for (i = 0; i < 8; i++) {
            int bits = 2 + (int)tk_rng_range(r, 51);
            nn[i] = tk_rng_nbits(r, bits) | 1;
            if (nn[i] < 3) nn[i] = 3;
            x[i] = tk_rng_u64(r) % nn[i];
            y[i] = tk_rng_u64(r) % nn[i];
            inv[i] = multiplicative_inverse52(nn[i]);
        }
        vn = _mm512_loadu_si512((const void *)nn);
        vx = _mm512_loadu_si512((const void *)x);
        vy = _mm512_loadu_si512((const void *)y);
        vinv = _mm512_loadu_si512((const void *)inv);

        v = mulredc52_x8(vx, vy, vn, vinv);
        _mm512_storeu_si512((void *)out, v);
        for (i = 0; i < 8; i++) {
            /* reference: x*y*2^-52 mod n, via the textbook REDC on 128-bit */
            uint64_t hi, lo = _umul128(x[i], y[i], &hi);
            uint64_t m = ((lo & M52) * inv[i]) & M52;       /* T - m*n = 0 mod 2^52 */
            uint64_t mh, ml = _umul128(m, nn[i], &mh);
            uint64_t t_hi = (hi << 12) | (lo >> 52);        /* T >> 52 */
            uint64_t mn_hi = (mh << 12) | (ml >> 52);       /* mN >> 52 */
            uint64_t e = (t_hi >= mn_hi) ? t_hi - mn_hi : t_hi + nn[i] - mn_hi;
            if (!TK_EQ_U64(tk, out[i], e)) break;
        }

        v = sqrredc52_x8(vx, vn, vinv);
        _mm512_storeu_si512((void *)out, v);
        v = mulredc52_x8(vx, vx, vn, vinv);
        {
            uint64_t o2[8];
            _mm512_storeu_si512((void *)o2, v);
            for (i = 0; i < 8; i++) TK_EQ_U64(tk, out[i], o2[i]);
        }

        v = addmod52_x8(vx, vy, vn);
        _mm512_storeu_si512((void *)out, v);
        for (i = 0; i < 8; i++)
            if (!TK_EQ_U64(tk, out[i], (x[i] + y[i]) % nn[i])) break;

        v = submod52_x8(vx, vy, vn);
        _mm512_storeu_si512((void *)out, v);
        for (i = 0; i < 8; i++)
            if (!TK_EQ_U64(tk, out[i], (x[i] + nn[i] - y[i]) % nn[i])) break;
    }
    if (tk_bench_enabled(tk)) {
        uint64_t nn[8], x[8], inv[8];
        __m512i vn, vx, vinv;
        for (i = 0; i < 8; i++) {
            nn[i] = tk_rng_nbits(r, 52) | 1;
            x[i] = tk_rng_u64(r) % nn[i];
            inv[i] = multiplicative_inverse52(nn[i]);
        }
        vn = _mm512_loadu_si512((const void *)nn);
        vx = _mm512_loadu_si512((const void *)x);
        vinv = _mm512_loadu_si512((const void *)inv);
        TK_BENCH(tk, "mulredc52_x8 latency (per lane)", 8, 0.0,
                 { vx = mulredc52_x8(vx, vx, vn, vinv); });
        _mm512_storeu_si512((void *)x, vx);
        tk_sink_u64 += x[0];
    }
}
#endif

/* ================================================================== *
 * limb1.c single-word helpers that test_sp_arith.c does not reach
 * ================================================================== */
static void t_sp_helpers(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    mpz_t a, b, c, d, e;
    int k;
    mpz_inits(a, b, c, d, e, NULL);
    for (k = 0; k < 100000; k++) {
        uint64_t u = tk_rng_u64(r), v = tk_rng_u64(r), w = tk_rng_u64(r), t = tk_rng_u64(r);
        uint64_t lo, hi, ex[3];
        if ((k % 5) == 0) u = ~(uint64_t)0;
        if ((k % 7) == 0) v = ~(uint64_t)0;
        if ((k % 11) == 0) w = ~(uint64_t)0;
        if ((k % 13) == 0) t = ~(uint64_t)0;
        if ((k % 17) == 0) u = 0;
        set_limbs(a, &u, 1); set_limbs(b, &v, 1); set_limbs(c, &w, 1); set_limbs(d, &t, 1);

        spAdd(u, v, &lo, &hi);
        mpz_add(e, a, b); get_limbs(ex, 2, e);
        if (!TK_EQ_U64(tk, lo, ex[0]) || !TK_EQ_U64(tk, hi, ex[1])) break;

        spAdd3(u, v, w, &lo, &hi);
        mpz_add(e, a, b); mpz_add(e, e, c); get_limbs(ex, 2, e);
        if (!TK_EQ_U64(tk, lo, ex[0]) || !TK_EQ_U64(tk, hi, ex[1])) break;

        /* differences: low word mod 2^64, borrow = how many times 2^64 was
           borrowed (0 or 1 for spSub, 0..2 for spSub3) */
        spSub(u, v, &lo, &hi);
        if (!TK_EQ_U64(tk, lo, u - v) || !TK_EQ_U64(tk, hi, u < v)) break;
        spSub3(u, v, w, &lo, &hi);
        if (!TK_EQ_U64(tk, lo, u - v - w)) break;
        if (!TK_EQ_U64(tk, hi, (uint64_t)(u < v) + (uint64_t)((u - v) < w))) break;

        /* u*v + w + t always fits in two words */
        spMulAdd(u, v, w, t, &lo, &hi);
        mpz_mul(e, a, b); mpz_add(e, e, c); mpz_add(e, e, d); get_limbs(ex, 3, e);
        if (!TK_EQ_U64(tk, lo, ex[0]) || !TK_EQ_U64(tk, hi, ex[1]) || !TK_EQ_U64(tk, ex[2], 0)) break;

        /* u64div(c, n) = c * 2^64 mod n, for c < n */
        {
            uint64_t n = tk_rng_nbits(r, 1 + (int)tk_rng_range(r, 64)) | 1, cc;
            if (n < 3) n = 3;
            cc = tk_rng_u64(r) % n;
            set_limbs(a, &cc, 1); set_limbs(b, &n, 1);
            mpz_mul_2exp(e, a, 64); mpz_mod(e, e, b); get_limbs(ex, 1, e);
            if (!TK_EQ_U64(tk, u64div(cc, n), ex[0])) break;
        }

        /* bit and digit counts */
        {
            uint64_t x = tk_rng_u64(r) >> tk_rng_range(r, 64);
            set_limbs(a, &x, 1);
            if (!TK_EQ_U64(tk, spBits(x), x ? mpz_sizeinbase(a, 2) : 0)) break;
            {
                char buf[32];           /* mpz_sizeinbase(.,10) may be one too big */
                gmp_sprintf(buf, "%Zd", a);
                if (!TK_EQ_U64(tk, (uint64_t)ndigits_1(x), (uint64_t)strlen(buf))) break;
            }
        }

        /* modinv_1b: inverse mod 2^32 of an odd a (p is ignored) */
        {
            uint32_t x = (uint32_t)tk_rng_u64(r) | 1u;
            if (!TK_EQ_U64(tk, (uint32_t)(modinv_1b(x, 0) * x), 1)) break;
        }
    }
    TK_EQ_U64(tk, (uint64_t)ndigits_1(0), 1);
    TK_EQ_U64(tk, (uint64_t)ndigits_1(9), 1);
    TK_EQ_U64(tk, (uint64_t)ndigits_1(10), 2);
    TK_EQ_U64(tk, (uint64_t)ndigits_1(~(uint64_t)0), 20);
    mpz_clears(a, b, c, d, e, NULL);
}

static void t_gcds(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    mpz_t a, b, g;
    int k;
    mpz_inits(a, b, g, NULL);
    for (k = 0; k < 100000; k++) {
        /* plant a common factor, vary sizes and powers of two */
        uint64_t f = tk_rng_nbits(r, 1 + (int)tk_rng_range(r, 24));
        uint64_t u = tk_rng_nbits(r, 1 + (int)tk_rng_range(r, 40)) * f;
        uint64_t v = tk_rng_nbits(r, 1 + (int)tk_rng_range(r, 40)) * f;
        uint64_t ref;
        double dw;
        set_limbs(a, &u, 1); set_limbs(b, &v, 1);
        mpz_gcd(g, a, b); get_limbs(&ref, 1, g);

        if (!TK_EQ_U64(tk, spGCD(u, v), ref)) break;
        if (!TK_EQ_U64(tk, gcd64(u, v), ref)) break;
        if (!TK_EQ_U64(tk, spBinGCD(u, v), ref)) break;         /* non-zero inputs */
        if (u & 1) {                                            /* u odd */
            if (!TK_EQ_U64(tk, spBinGCD_odd(u, v), ref)) break;
            if (!TK_EQ_U64(tk, bingcd64(u, v), ref)) break;
        }
        if (u < ((uint64_t)1 << 52) && v < ((uint64_t)1 << 52)) {   /* exact in doubles */
            dblGCD((double)u, (double)v, &dw);
            if (!TK_EQ_U64(tk, (uint64_t)dw, ref)) break;
        }
    }
    TK_EQ_U64(tk, spGCD(0, 7), 7);
    TK_EQ_U64(tk, spGCD(7, 0), 7);
    TK_EQ_U64(tk, gcd64(0, 7), 7);

    for (k = 0; k < 50000; k++) {
        uint64_t f[2] = { 0, 0 }, uu[2], vv[2], w[2], e[2], t[4], p[2], q[2];
        /* u = f * p, v = f * q with f up to 60 bits and p, q up to 64 */
        f[0] = tk_rng_nbits(r, 1 + (int)tk_rng_range(r, 60));
        p[0] = tk_rng_nbits(r, 1 + (int)tk_rng_range(r, 64)); p[1] = 0;
        q[0] = tk_rng_nbits(r, 1 + (int)tk_rng_range(r, 64)); q[1] = 0;
        mul128(f, p, t); uu[0] = t[0]; uu[1] = t[1];
        mul128(f, q, t); vv[0] = t[0]; vv[1] = t[1];
        set_limbs(a, uu, 2); set_limbs(b, vv, 2);
        mpz_gcd(g, a, b); get_limbs(e, 2, g);

        gcd128(uu, vv, w);
        if (!TK_CHECKF(tk, w[0] == e[0] && w[1] == e[1], "gcd128 (k=%d)", k)) break;

#if defined(__SIZEOF_INT128__) && (__SIZEOF_INT128__ == 16)
        /* bin_gcd128 wants u odd and works on a scratch copy of v */
        if (uu[0] & 1) {
            uint64_t vs[2] = { vv[0], vv[1] }, us[2] = { uu[0], uu[1] };
            bin_gcd128(us, vs, w);
            if (!TK_CHECKF(tk, w[0] == e[0] && w[1] == e[1],
                           "bin_gcd128 (k=%d) u=%016llx%016llx v=%016llx%016llx", k,
                           (unsigned long long)uu[1], (unsigned long long)uu[0],
                           (unsigned long long)vv[1], (unsigned long long)vv[0])) break;
        }
#endif
    }
    mpz_clears(a, b, g, NULL);
}

static void t_jacobi_sqrt(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    mpz_t a, p, s;
    int k;
    mpz_inits(a, p, s, NULL);

    /* jacobi_1(n, p) for any odd p, against GMP */
    for (k = 0; k < 100000; k++) {
        uint64_t pp = tk_rng_nbits(r, 2 + (int)tk_rng_range(r, 31)) | 1, nn;
        if (pp < 3) pp = 3;
        nn = tk_rng_u64(r) >> tk_rng_range(r, 64);
        if ((k % 9) == 0) nn = pp * (1 + tk_rng_range(r, 5));       /* shares a factor: 0 */
        set_limbs(a, &nn, 1); set_limbs(p, &pp, 1);
        if (!TK_CHECKF(tk, jacobi_1(nn, pp) == mpz_jacobi(a, p), "jacobi_1(%llu, %llu)",
                       (unsigned long long)nn, (unsigned long long)pp)) break;
    }
    TK_CHECK(tk, jacobi_1(5, 8) == -2);                              /* even p: error code */

    /* Tonelli-Shanks on real primes, including p = 1 mod 2^e for large e
       (the case that needs many squarings; 2^32 + 15*2^... style primes) */
    {
        static const uint64_t fixed[] = { 3, 5, 7, 13, 17, 97, 193, 257, 12289, 65537,
            40961, 786433, 5767169, 7340033, 469762049ULL, 998244353ULL, 2013265921ULL,
            3221225473ULL, 4179340454199820289ULL, 1945555039024054273ULL };
        int nf = (int)(sizeof fixed / sizeof fixed[0]), i, j;
        for (i = 0; i < nf + 300; i++) {
            uint64_t pp, sq = 0;
            if (i < nf) pp = fixed[i];
            else {
                pp = tk_rng_nbits(r, 3 + (int)tk_rng_range(r, 60));
                set_limbs(p, &pp, 1); mpz_nextprime(p, p);
                if (mpz_sizeinbase(p, 2) > 63) continue;
                get_limbs(&pp, 1, p);
            }
            set_limbs(p, &pp, 1);
            for (j = 0; j < 20; j++) {
                uint64_t x = 1 + tk_rng_u64(r) % (pp - 1), aa, back;
                set_limbs(s, &x, 1); mpz_mul(a, s, s); mpz_mod(a, a, p);
                get_limbs(&aa, 1, a);                               /* a residue */
                ShanksTonelli_1(aa, pp, &sq);
                set_limbs(s, &sq, 1); mpz_mul(s, s, s); mpz_mod(s, s, p);
                get_limbs(&back, 1, s);
                if (!TK_CHECKF(tk, sq < pp && back == aa, "ShanksTonelli_1(%llu mod %llu) = %llu",
                               (unsigned long long)aa, (unsigned long long)pp, (unsigned long long)sq)) {
                    i = nf + 300; break;
                }
            }
        }
    }
    mpz_clears(a, p, s, NULL);
}

/* ================================================================== *
 * 128-bit Montgomery: what test_modular.c's monty128 does not cover --
 * the ...n variants, dblmod / chkmod, operand aliasing, results that must
 * be fully reduced, and moduli pinned to the top of the range.
 * ================================================================== */
static void t_monty128_edges(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    mpz_t n, x, y, e, R, Rinv;
    int k;
    mpz_inits(n, x, y, e, R, Rinv, NULL);
    mpz_setbit(R, 128);
    for (k = 0; k < 40000; k++) {
        monty128_t m;
        uint64_t nn[2], xx[2], yy[2], ee[2], w[2];
        int kind = k % 5;

        /* moduli: top bit set, just below 2^128, 65-bit, one limb, random */
        nn[0] = tk_rng_u64(r) | 1; nn[1] = tk_rng_u64(r);
        if (kind == 0) nn[1] |= (uint64_t)1 << 63;
        if (kind == 1) { nn[1] = ~(uint64_t)0; nn[0] |= 0xffffffff00000000ULL; }
        if (kind == 2) nn[1] = 1;
        if (kind == 3) { nn[1] = 0; if (nn[0] < 3) nn[0] = 3; }
        set_limbs(n, nn, 2);

        /* residues, biased towards n - 1 and small values */
        xx[0] = tk_rng_u64(r); xx[1] = tk_rng_u64(r); set_limbs(x, xx, 2); mpz_mod(x, x, n);
        yy[0] = tk_rng_u64(r); yy[1] = tk_rng_u64(r); set_limbs(y, yy, 2); mpz_mod(y, y, n);
        if ((k % 7) == 0) mpz_sub_ui(x, n, 1);
        if ((k % 11) == 0) mpz_sub_ui(y, n, 1);
        if ((k % 13) == 0) mpz_set_ui(x, 0);
        get_limbs(xx, 2, x); get_limbs(yy, 2, y);

        monty128_init(&m, nn);
        mpz_invert(Rinv, R, n);

#define EXPECT(val, name)                                                        \
        do { get_limbs(ee, 2, e);                                                \
             if (!TK_CHECKF(tk, (val)[0] == ee[0] && (val)[1] == ee[1],          \
                 "%s (k=%d kind=%d) n=%016llx%016llx", name, k, kind,            \
                 (unsigned long long)nn[1], (unsigned long long)nn[0])) goto done; \
        } while (0)

        /* products: exact, fully reduced value x*y*R^-1 mod n, from every entry point */
        mpz_mul(e, x, y); mpz_mul(e, e, Rinv); mpz_mod(e, e, n);
        mulmod128(xx, yy, w, &m);               EXPECT(w, "mulmod128 reduced");
        mulmod128n(xx, yy, w, nn, m.rho);       EXPECT(w, "mulmod128n");
        w[0] = xx[0]; w[1] = xx[1];
        mulmod128(w, yy, w, &m);                EXPECT(w, "mulmod128 w aliases u");
        w[0] = yy[0]; w[1] = yy[1];
        mulmod128(xx, w, w, &m);                EXPECT(w, "mulmod128 w aliases v");

        mpz_mul(e, x, x); mpz_mul(e, e, Rinv); mpz_mod(e, e, n);
        sqrmod128(xx, w, &m);                   EXPECT(w, "sqrmod128 reduced");
        sqrmod128n(xx, w, nn, m.rho);           EXPECT(w, "sqrmod128n");
        w[0] = xx[0]; w[1] = xx[1];
        sqrmod128(w, w, &m);                    EXPECT(w, "sqrmod128 in place");
        mulmod128(xx, xx, w, &m);               EXPECT(w, "mulmod128 u == v");

        /* add / sub / double, including every aliasing pattern */
        mpz_add(e, x, y); mpz_mod(e, e, n);
        addmod128(xx, yy, w, nn);               EXPECT(w, "addmod128");
        w[0] = xx[0]; w[1] = xx[1];
        addmod128(w, yy, w, nn);                EXPECT(w, "addmod128 w aliases a");
        w[0] = yy[0]; w[1] = yy[1];
        addmod128(xx, w, w, nn);                EXPECT(w, "addmod128 w aliases b");

        mpz_add(e, x, x); mpz_mod(e, e, n);
        addmod128(xx, xx, w, nn);               EXPECT(w, "addmod128 a == b");
        w[0] = xx[0]; w[1] = xx[1];
        addmod128(w, w, w, nn);                 EXPECT(w, "addmod128 all aliased");
        w[0] = xx[0]; w[1] = xx[1];
        dblmod128(w, nn);                       EXPECT(w, "dblmod128");

        mpz_sub(e, x, y); mpz_mod(e, e, n);
        submod128(xx, yy, w, nn);               EXPECT(w, "submod128");
        w[0] = xx[0]; w[1] = xx[1];
        submod128(w, yy, w, nn);                EXPECT(w, "submod128 w aliases a");
        w[0] = yy[0]; w[1] = yy[1];
        submod128(xx, w, w, nn);                EXPECT(w, "submod128 w aliases b");
        mpz_set_ui(e, 0);
        submod128(xx, xx, w, nn);               EXPECT(w, "submod128 a == b");

        /* chkmod128: a in [0, 2n) -> a mod n (only when 2n fits in 128 bits) */
        if (!(nn[1] >> 63)) {
            mpz_add(e, x, n); get_limbs(w, 2, e);       /* x + n, in [n, 2n) */
            mpz_set(e, x);
            chkmod128(w, nn);                   EXPECT(w, "chkmod128 above n");
        }
        w[0] = xx[0]; w[1] = xx[1];
        mpz_set(e, x);
        chkmod128(w, nn);                       EXPECT(w, "chkmod128 below n");
        w[0] = nn[0]; w[1] = nn[1];
        mpz_set_ui(e, 0);
        chkmod128(w, nn);                       EXPECT(w, "chkmod128 at n");
#undef EXPECT
    }
done:
    if (tk_bench_enabled(tk)) {
        monty128_t m;
        uint64_t nn[2], xx[2], yy[2];
        nn[0] = tk_rng_u64(r) | 1; nn[1] = tk_rng_u64(r) | ((uint64_t)1 << 62);
        xx[0] = tk_rng_u64(r); xx[1] = tk_rng_u64(r) >> 3;
        yy[0] = tk_rng_u64(r); yy[1] = tk_rng_u64(r) >> 3;
        monty128_init(&m, nn);
        TK_BENCH(tk, "mulmod128 latency", 1, 0.0, { mulmod128(xx, yy, xx, &m); });
        TK_BENCH(tk, "sqrmod128 latency", 1, 0.0, { sqrmod128(xx, xx, &m); });
        TK_BENCH(tk, "addmod128 latency", 1, 0.0, { addmod128(xx, yy, xx, nn); });
        TK_BENCH(tk, "monty128_init", 1, 0.0, { nn[0] += 2; monty128_init(&m, nn); });
        TK_BENCH(tk, "to_monty128", 1, 0.0, { to_monty128(&m, yy); });
        tk_sink_u64 += xx[0] + yy[0];
    }
    mpz_clears(n, x, y, e, R, Rinv, NULL);
}

/* zero-safe 104/128-bit scans */
static void t_scan128(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    int k, b;
    TK_EQ_U64(tk, my_clz128(0, 0), 128);
    TK_EQ_U64(tk, my_ctz128(0, 0), 128);
    TK_EQ_U64(tk, my_ctz104(0, 0), 104);
    TK_EQ_U64(tk, my_clz104(0, 0), 104);
    for (b = 0; b < 128; b++) {
        uint64_t lo = (b < 64) ? (uint64_t)1 << b : 0, hi = (b >= 64) ? (uint64_t)1 << (b - 64) : 0;
        TK_EQ_U64(tk, my_ctz128(lo, hi), (uint64_t)b);
        TK_EQ_U64(tk, my_clz128(lo, hi), (uint64_t)(127 - b));
    }
    for (b = 0; b < 104; b++) {
        /* 104-bit values are two 52-bit words */
        uint64_t lo = (b < 52) ? (uint64_t)1 << b : 0, hi = (b >= 52) ? (uint64_t)1 << (b - 52) : 0;
        TK_EQ_U64(tk, my_ctz104(lo, hi), (uint64_t)b);
        /* my_clz104 counts from bit 63 of the word holding the top set bit
           (it does not discount the 12 padding bits) -- long-standing
           behaviour that modexp callers rely on; pinned here so a change is
           deliberate. */
        TK_EQ_U64(tk, my_clz104(lo, hi), (b >= 52) ? (uint64_t)(63 - (b - 52)) : (uint64_t)(52 + 63 - b));
    }
    for (k = 0; k < 50000; k++) {
        uint64_t lo = tk_rng_u64(r) << tk_rng_range(r, 64), hi = tk_rng_u64(r) >> tk_rng_range(r, 64);
        uint64_t tz, lz;
        if ((k & 3) == 0) hi = 0;
        if ((k & 7) == 1) lo = 0;
        if ((lo | hi) == 0) continue;
        for (tz = 0; tz < 128; tz++)
            if (((tz < 64) ? (lo >> tz) : (hi >> (tz - 64))) & 1) break;
        for (lz = 0; lz < 128; lz++)
            if (((lz < 64) ? (hi >> (63 - lz)) : (lo >> (127 - lz))) & 1) break;
        if (!TK_EQ_U64(tk, my_ctz128(lo, hi), tz)) break;
        if (!TK_EQ_U64(tk, my_clz128(lo, hi), lz)) break;
    }
}

#if MP_USE_AVX512F && !defined(LIMB2_NO_VECTOR) && !defined(LIMB1_SCALAR_ONLY)
/* ================================================================== *
 * 104-bit vector Montgomery kernels (8 lanes of 2 x 52 bits) and the
 * 52-bit conditional-double multiply, checked lane by lane against GMP.
 * Layout: x[0..7] low 52-bit words, x[8..15] high 52-bit words.
 * ================================================================== */
#define M52 0x000fffffffffffffULL

static void set104(mpz_t z, const uint64_t *x, int lane)
{
    mpz_set_ui(z, 0);
    {
        uint64_t hi = x[8 + lane], lo = x[lane];
        mpz_t t; mpz_init(t);
        set_limbs(z, &hi, 1); mpz_mul_2exp(z, z, 52);
        set_limbs(t, &lo, 1); mpz_add(z, z, t);
        mpz_clear(t);
    }
}
static void get104(uint64_t *x, int lane, const mpz_t z)
{
    uint64_t w[2];
    get_limbs(w, 2, z);
    x[lane] = w[0] & M52;
    x[8 + lane] = (w[0] >> 52) | (w[1] << 12);
}

static void t_redc104_x8(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    mpz_t n[8], x[8], y[8], e, R, Rinv, got;
    gmp_randstate_t gs;
    int k, i;
    for (i = 0; i < 8; i++) mpz_inits(n[i], x[i], y[i], NULL);
    mpz_inits(e, R, Rinv, got, NULL);
    mpz_setbit(R, 104);
    gmp_randinit_default(gs);
    gmp_randseed_ui(gs, (unsigned long)tk_rng_u64(r));

    for (k = 0; k < 20000; k++) {
        uint64_t nn[16], xx[16], yy[16], o[16];
        __m512i n0, n1, x0, x1, y0, y1, c0, c1, vrho, vrho_pos;
        const __m512i zero = _mm512_setzero_si512(), mask = _mm512_set1_epi64((long long)M52);
        int bits = 54 + (int)tk_rng_range(r, 51);           /* 54 .. 104 */
        __mmask8 msk = (__mmask8)tk_rng_u64(r);
        if ((k % 4) == 0) bits = 104;
        if ((k % 9) == 0) msk = 0xff;

        for (i = 0; i < 8; i++) {
            mpz_urandomb(n[i], gs, (mp_bitcnt_t)bits);
            mpz_setbit(n[i], (mp_bitcnt_t)(bits - 1)); mpz_setbit(n[i], 0);
            mpz_urandomm(x[i], gs, n[i]); mpz_urandomm(y[i], gs, n[i]);
            if ((k % 7) == 0) mpz_sub_ui(x[i], n[i], 1);
            if ((k % 11) == 0) mpz_sub_ui(y[i], n[i], 1);
            if ((k % 13) == 0 && (i & 1)) mpz_set_ui(x[i], 0);
            get104(nn, i, n[i]); get104(xx, i, x[i]); get104(yy, i, y[i]);
        }
        n0 = _mm512_loadu_si512((const void *)&nn[0]); n1 = _mm512_loadu_si512((const void *)&nn[8]);
        x0 = _mm512_loadu_si512((const void *)&xx[0]); x1 = _mm512_loadu_si512((const void *)&xx[8]);
        y0 = _mm512_loadu_si512((const void *)&yy[0]); y1 = _mm512_loadu_si512((const void *)&yy[8]);

        /* 52-bit inverse of the low word: positive, and its negation */
        vrho_pos = multiplicative_inverse104_x8(nn);
        _mm512_storeu_si512((void *)o, vrho_pos);
        for (i = 0; i < 8; i++)
            if (!TK_EQ_U64(tk, (o[i] * nn[i]) & M52, 1)) goto done;
        vrho_pos = _mm512_and_si512(vrho_pos, mask);
        vrho = _mm512_and_si512(_mm512_sub_epi64(zero, vrho_pos), mask);

#define STORE() do { _mm512_storeu_si512((void *)&o[0], c0); _mm512_storeu_si512((void *)&o[8], c1); } while (0)
/* lane i holds exactly `e` (a fully reduced value) */
#define LANE_IS(name) do { set104(got, o, i);                                        \
        if (!TK_CHECKF(tk, o[i] <= M52 && o[8 + i] <= M52 && mpz_cmp(got, e) == 0,   \
            "%s lane %d (k=%d bits=%d)%s%s", name, i, k, bits,                       \
            mpz_congruent_p(got, e, n[i]) ? " [congruent mod n]" : " [NOT congruent]", \
            (mpz_cmp(got, n[i]) >= 0) ? " [>= n]" : "")) goto done; } while (0)

/* lane i is congruent to `e` mod n and fits in 104 bits.  The plain
   multiply / square kernels reduce lazily: their result can be >= n (seen
   for moduli of 103-104 bits).  Only mask_sqrredc104_exact_vec, and the
   add / sub / chk kernels, return fully reduced values. */
#define LANE_CONG(name) do { set104(got, o, i);                                      \
        if (!TK_CHECKF(tk, o[i] <= M52 && o[8 + i] <= M52 &&                         \
            mpz_congruent_p(got, e, n[i]), "%s lane %d (k=%d bits=%d): not congruent", \
            name, i, k, bits)) goto done; } while (0)

        /* x * y * 2^-104 mod n */
        mulredc104_vec(&c1, &c0, x1, x0, y1, y0, n1, n0, vrho); STORE();
        for (i = 0; i < 8; i++) {
            mpz_invert(Rinv, R, n[i]);
            mpz_mul(e, x[i], y[i]); mpz_mul(e, e, Rinv); mpz_mod(e, e, n[i]);
            LANE_CONG("mulredc104_vec");
        }
        mask_mulredc104_vec(&c1, &c0, msk, x1, x0, y1, y0, n1, n0, vrho); STORE();
        for (i = 0; i < 8; i++) {
            if (!((msk >> i) & 1)) continue;
            mpz_invert(Rinv, R, n[i]);
            mpz_mul(e, x[i], y[i]); mpz_mul(e, e, Rinv); mpz_mod(e, e, n[i]);
            LANE_CONG("mask_mulredc104_vec");
        }

        /* x^2 * 2^-104 mod n, three entry points */
        sqrredc104_vec(&c1, &c0, x1, x0, n1, n0, vrho); STORE();
        for (i = 0; i < 8; i++) {
            mpz_invert(Rinv, R, n[i]);
            mpz_mul(e, x[i], x[i]); mpz_mul(e, e, Rinv); mpz_mod(e, e, n[i]);
            LANE_CONG("sqrredc104_vec");
        }
        mask_sqrredc104_vec(&c1, &c0, msk, x1, x0, n1, n0, vrho); STORE();
        for (i = 0; i < 8; i++) {
            if (!((msk >> i) & 1)) continue;
            mpz_invert(Rinv, R, n[i]);
            mpz_mul(e, x[i], x[i]); mpz_mul(e, e, Rinv); mpz_mod(e, e, n[i]);
            LANE_CONG("mask_sqrredc104_vec");
        }
        mask_sqrredc104_exact_vec(&c1, &c0, msk, x1, x0, n1, n0, vrho); STORE();
        for (i = 0; i < 8; i++) {
            if (!((msk >> i) & 1)) continue;
            mpz_invert(Rinv, R, n[i]);
            mpz_mul(e, x[i], x[i]); mpz_mul(e, e, Rinv); mpz_mod(e, e, n[i]);
            LANE_IS("mask_sqrredc104_exact_vec");
        }
        /* mask_sqrredc104_vec_pos (positive-inverse experiment, no callers) has
           its own conventions and is deliberately not pinned here. */
        (void)vrho_pos;

        /* modular add / sub / double */
        addmod104_x8(&c1, &c0, x1, x0, y1, y0, n1, n0); STORE();
        for (i = 0; i < 8; i++) { mpz_add(e, x[i], y[i]); mpz_mod(e, e, n[i]); LANE_IS("addmod104_x8"); }
        submod104_x8(&c1, &c0, x1, x0, y1, y0, n1, n0); STORE();
        for (i = 0; i < 8; i++) { mpz_sub(e, x[i], y[i]); mpz_mod(e, e, n[i]); LANE_IS("submod104_x8"); }
        mask_addmod104_x8(&c1, &c0, msk, x1, x0, y1, y0, n1, n0); STORE();
        for (i = 0; i < 8; i++) {
            if ((msk >> i) & 1) { mpz_add(e, x[i], y[i]); mpz_mod(e, e, n[i]); }
            else mpz_set(e, x[i]);
            LANE_IS("mask_addmod104_x8");
        }
        mask_dblmod104_x8(&c1, &c0, msk, x1, x0, n1, n0); STORE();
        for (i = 0; i < 8; i++) {
            if ((msk >> i) & 1) { mpz_add(e, x[i], x[i]); mpz_mod(e, e, n[i]); }
            else mpz_set(e, x[i]);
            LANE_IS("mask_dblmod104_x8");
        }

        /* conditional subtract of n: on x (already reduced: unchanged) and on n itself (-> 0) */
        chkmod104_x8(&c1, &c0, x1, x0, n1, n0); STORE();
        for (i = 0; i < 8; i++) { mpz_set(e, x[i]); LANE_IS("chkmod104_x8 below n"); }
        chkmod104_x8(&c1, &c0, n1, n0, n1, n0); STORE();
        for (i = 0; i < 8; i++) { mpz_set_ui(e, 0); LANE_IS("chkmod104_x8 at n"); }
#undef STORE
#undef LANE_IS
#undef LANE_CONG

        /* 52-bit multiply with conditional doubling: (x*y*2^-52 mod n) * (1 or 2) mod n */
        {
            uint64_t n52[8], a52[8], b52[8], o52[8], rho52[8];
            __m512i vn, va, vb, vr, vc;
            for (i = 0; i < 8; i++) {
                n52[i] = (tk_rng_nbits(r, 2 + (int)tk_rng_range(r, 50)) | 1);
                if (n52[i] < 3) n52[i] = 3;
                a52[i] = tk_rng_u64(r) % n52[i]; b52[i] = tk_rng_u64(r) % n52[i];
                rho52[i] = ((uint64_t)0 - multiplicative_inverse52(n52[i])) & M52;
            }
            vn = _mm512_loadu_si512((const void *)n52); va = _mm512_loadu_si512((const void *)a52);
            vb = _mm512_loadu_si512((const void *)b52); vr = _mm512_loadu_si512((const void *)rho52);
            mulredc52_mask_add_vec(&vc, msk, va, vb, vn, vr);
            _mm512_storeu_si512((void *)o52, vc);
            for (i = 0; i < 8; i++) {
                uint64_t ex;
                mpz_t nz, t; mpz_inits(nz, t, NULL);
                set_limbs(nz, &n52[i], 1);
                mpz_set_ui(R, 0); mpz_setbit(R, 52); mpz_invert(Rinv, R, nz);
                set_limbs(e, &a52[i], 1); set_limbs(t, &b52[i], 1);
                mpz_mul(e, e, t); mpz_mul(e, e, Rinv);
                if ((msk >> i) & 1) mpz_mul_2exp(e, e, 1);
                mpz_mod(e, e, nz); get_limbs(&ex, 1, e);
                mpz_clears(nz, t, NULL);
                if (!TK_CHECKF(tk, o52[i] == ex, "mulredc52_mask_add_vec lane %d n=%llx",
                               i, (unsigned long long)n52[i])) goto done;
            }
            mpz_set_ui(R, 0); mpz_setbit(R, 104);
        }
    }
done:
    if (tk_bench_enabled(tk)) {
        uint64_t nn[16], xx[16];
        __m512i n0, n1, x0, x1, vrho;
        for (i = 0; i < 8; i++) {
            nn[i] = (tk_rng_u64(r) & M52) | 1; nn[8 + i] = (tk_rng_u64(r) & M52) | ((uint64_t)1 << 51);
            xx[i] = tk_rng_u64(r) & M52; xx[8 + i] = tk_rng_u64(r) & (M52 >> 1);
        }
        n0 = _mm512_loadu_si512((const void *)&nn[0]); n1 = _mm512_loadu_si512((const void *)&nn[8]);
        x0 = _mm512_loadu_si512((const void *)&xx[0]); x1 = _mm512_loadu_si512((const void *)&xx[8]);
        vrho = _mm512_and_si512(_mm512_sub_epi64(_mm512_setzero_si512(), multiplicative_inverse104_x8(nn)),
                                _mm512_set1_epi64((long long)M52));
        TK_BENCH(tk, "sqrredc104_vec latency (per lane)", 8, 0.0,
                 { sqrredc104_vec(&x1, &x0, x1, x0, n1, n0, vrho); });
        TK_BENCH(tk, "mulredc104_vec latency (per lane)", 8, 0.0,
                 { mulredc104_vec(&x1, &x0, x1, x0, x1, x0, n1, n0, vrho); });
        _mm512_storeu_si512((void *)xx, x0); tk_sink_u64 += xx[0];
    }
    gmp_randclear(gs);
    for (i = 0; i < 8; i++) mpz_clears(n[i], x[i], y[i], NULL);
    mpz_clears(e, R, Rinv, got, NULL);
}
#endif

static const tk_test tk__limb_tests[] = {
    { "addsub_shift128", t_addsub_shift,   "fast limb2 arith gmp" },
    { "mul_sqr128",      t_mul_sqr,        "fast limb2 arith gmp" },
    { "div128",          t_div128,         "fast limb2 arith gmp" },
    { "monty128_setup",  t_monty128_setup, "fast limb2 modular gmp" },
    { "monty128_edges",  t_monty128_edges, "fast limb2 modular gmp" },
    { "scan128",         t_scan128,        "fast limb2 bitscan" },
    { "redc64_pos",      t_redc64_pos,     "fast limb1 modular gmp" },
    { "spPRP2",          t_spPRP2,         "fast limb1 modular gmp" },
    { "sp_helpers",      t_sp_helpers,     "fast limb1 arith sp gmp" },
    { "gcds",            t_gcds,           "fast limb1 limb2 gcd gmp" },
    { "jacobi_sqrt",     t_jacobi_sqrt,    "fast limb1 modular gmp" },
#if MP_USE_AVX512F && !defined(LIMB1_NO_VECTOR)
    { "redc52_x8",       t_redc52_x8,      "fast limb1 modular vector gmp" },
#endif
#if MP_USE_AVX512F && !defined(LIMB2_NO_VECTOR) && !defined(LIMB1_SCALAR_ONLY)
    { "redc104_x8",      t_redc104_x8,     "fast limb2 modular vector gmp" },
#endif
};

const tk_module tk_module_limb = {
    "limb",
    "arith/limb vs GMP: two-limb arithmetic, 128-bit Montgomery, sp helpers, gcd, jacobi, 104-bit x8",
    tk__limb_tests,
    (int)(sizeof tk__limb_tests / sizeof tk__limb_tests[0])
};
