/*----------------------------------------------------------------------
 Layer 0 -- arith/limb header-only kernels.

 Everything here is an inline function from mp_platform.h, limb1.h or
 limb2.h, so (like the rest of Layer 0) this module needs no GMP and links
 no YAFU object.  It fills the gaps left by mp_arith / mp_bitscan:

   mp_platform.h  32-bit carry ops, 32-bit bit scans, popcount, the
                  software 128/64 divide
   limb1.h        addmod / submod (64 and 32 bit), every 64-bit Montgomery
                  flavour (positive and negative inverse, 63- and 60-bit
                  variants, multiply-add), inverses mod 2^64 / 2^52, the
                  binary gcd, the zero-safe 52-bit scans
   limb2.h        two-limb add / sub / compare / shift / multiply / square
   AVX-512 only   the 8-lane 52-bit Montgomery kernels, the 52 x 52
                  product halves (both the shared and the legacy forms),
                  and the 52-bit carry primitives

 References are deliberately slow and intrinsic-free: bit-at-a-time shifts,
 adds and compares on limb arrays.  Public domain.
----------------------------------------------------------------------*/
#include "testkit.h"
#include "limb1.h"
#include "limb2.h"
#include <string.h>

/* ---- portable references on little-endian limb arrays ---- */
static int ref_cmp(const uint64_t *a, const uint64_t *b, int n)
{
    int i;
    for (i = n - 1; i >= 0; i--)
        if (a[i] != b[i]) return (a[i] < b[i]) ? -1 : 1;
    return 0;
}
static uint64_t ref_add(uint64_t *r, const uint64_t *a, const uint64_t *b, int n)
{
    uint64_t c = 0; int i;
    for (i = 0; i < n; i++) {
        uint64_t s = a[i] + c, c1 = (s < c);
        r[i] = s + b[i]; c = c1 | (r[i] < s);
    }
    return c;
}
static uint64_t ref_sub(uint64_t *r, const uint64_t *a, const uint64_t *b, int n)
{
    uint64_t c = 0; int i;
    for (i = 0; i < n; i++) {
        uint64_t d = a[i] - c, c1 = (a[i] < c);
        r[i] = d - b[i]; c = c1 | (d < b[i]);
    }
    return c;
}
static void ref_shl1(uint64_t *a, int n)
{
    uint64_t c = 0; int i;
    for (i = 0; i < n; i++) { uint64_t t = a[i] >> 63; a[i] = (a[i] << 1) | c; c = t; }
}
static void ref_shr1(uint64_t *a, int n)
{
    uint64_t c = 0; int i;
    for (i = n - 1; i >= 0; i--) { uint64_t t = a[i] & 1; a[i] = (a[i] >> 1) | (c << 63); c = t; }
}
static void ref_mul64(uint64_t a, uint64_t b, uint64_t *lo, uint64_t *hi)
{
    uint64_t al = (uint32_t)a, ah = a >> 32, bl = (uint32_t)b, bh = b >> 32;
    uint64_t ll = al * bl, lh = al * bh, hl = ah * bl, hh = ah * bh;
    uint64_t mid = (ll >> 32) + (uint32_t)lh + (uint32_t)hl;
    *lo = (mid << 32) | (uint32_t)ll;
    *hi = hh + (lh >> 32) + (hl >> 32) + (mid >> 32);
}
static void ref_mul(uint64_t *r, const uint64_t *a, int n, const uint64_t *b, int m)
{
    int i, j;
    memset(r, 0, (size_t)(n + m) * sizeof(uint64_t));
    for (i = 0; i < n; i++) {
        uint64_t carry = 0;
        for (j = 0; j < m; j++) {
            uint64_t lo, hi, s;
            ref_mul64(a[i], b[j], &lo, &hi);
            s = r[i + j] + lo; hi += (s < lo);
            r[i + j] = s + carry; hi += (r[i + j] < s);
            carry = hi;
        }
        r[i + m] = carry;
    }
}
/* (hi:lo) mod d, restoring division one bit at a time */
static uint64_t ref_mod(uint64_t hi, uint64_t lo, uint64_t d)
{
    uint64_t r = 0;
    int i;
    for (i = 127; i >= 0; i--) {
        uint64_t top = r >> 63;
        uint64_t bit = (i >= 64) ? (hi >> (i - 64)) & 1 : (lo >> i) & 1;
        r = (r << 1) | bit;
        if (top || r >= d) r -= d;
    }
    return r;
}
static uint64_t ref_gcd(uint64_t a, uint64_t b)
{
    while (b) { uint64_t t = a % b; a = b; b = t; }
    return a;
}
/* random odd modulus of a random bit length in [2, maxbits] */
static uint64_t rand_modulus(tk_rng *r, int maxbits)
{
    int bits = 2 + (int)tk_rng_range(r, (uint64_t)(maxbits - 1));
    uint64_t n = tk_rng_nbits(r, bits) | 1;
    return (n < 3) ? 3 : n;
}

/* ================================================================== */
static void t_carry32(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    static const uint32_t corner[] = { 0u, 1u, 0xffffffffu, 0x80000000u, 0x7fffffffu };
    int k, i, j, c;
    for (i = 0; i < 5; i++) for (j = 0; j < 5; j++) for (c = 0; c < 2; c++) {
        uint32_t x = corner[i], y = corner[j], z;
        mp_carry_t co = mp_addcarry_u32((mp_carry_t)c, x, y, &z);
        TK_EQ_U64(tk, z, (uint32_t)(x + y + (uint32_t)c));
        TK_EQ_U64(tk, co, ((uint64_t)x + y + (uint64_t)c) >> 32);
        co = mp_subborrow_u32((mp_carry_t)c, x, y, &z);
        TK_EQ_U64(tk, z, (uint32_t)(x - y - (uint32_t)c));
        TK_EQ_U64(tk, co, (uint64_t)x < (uint64_t)y + (uint64_t)c);
    }
    for (k = 0; k < 200000; k++) {
        uint32_t x = (uint32_t)tk_rng_u64(r), y = (uint32_t)tk_rng_u64(r), z;
        mp_carry_t cin = (mp_carry_t)(k & 1), co;
        co = mp_addcarry_u32(cin, x, y, &z);
        if (!TK_EQ_U64(tk, z, (uint32_t)(x + y + cin))) break;
        if (!TK_EQ_U64(tk, co, ((uint64_t)x + y + cin) >> 32)) break;
        co = mp_subborrow_u32(cin, x, y, &z);
        if (!TK_EQ_U64(tk, z, (uint32_t)(x - y - cin))) break;
        if (!TK_EQ_U64(tk, co, (uint64_t)x < (uint64_t)y + cin)) break;
    }
}

static void t_scan32_popcount(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    int k, b;
    for (b = 0; b < 32; b++) {
        TK_EQ_U64(tk, _trail_zcnt((uint32_t)1 << b), (uint64_t)b);
        TK_EQ_U64(tk, _trail_zcnt(0xffffffffu << b), (uint64_t)b);
    }
    TK_EQ_U64(tk, mp_popcount64(0), 0);
    TK_EQ_U64(tk, mp_popcount64(~(uint64_t)0), 64);
    TK_EQ_U64(tk, mp_popcount32(0), 0);
    TK_EQ_U64(tk, mp_popcount32(0xffffffffu), 32);
    /* the 52-bit zero-safe scans: full 64-bit scan for nonzero, 52 at zero */
    TK_EQ_U64(tk, _lead_full_zcnt52(0), 52);
    TK_EQ_U64(tk, _trail_full_zcnt52(0), 52);
    for (k = 0; k < 100000; k++) {
        uint64_t x = tk_rng_u64(r) >> tk_rng_range(r, 64), t;
        uint32_t x32 = (uint32_t)x, c = 0, c32 = 0;
        for (t = x; t; t &= t - 1) c++;
        for (t = x32; t; t &= t - 1) c32++;
        if (!TK_EQ_U64(tk, mp_popcount64(x), c)) break;
        if (!TK_EQ_U64(tk, mp_popcount32(x32), c32)) break;
        if (x32) {
            uint32_t tz = 0, y = x32;
            while (!((x32 >> tz) & 1)) tz++;
            if (!TK_EQ_U64(tk, _trail_zcnt(x32), tz)) break;
            y = _reset_lsb(y);
            if (!TK_EQ_U64(tk, y, x32 & (x32 - 1))) break;
        }
        if (x) {
            TK_EQ_U64(tk, _lead_full_zcnt52(x), _lead_zcnt64(x));
            TK_EQ_U64(tk, _trail_full_zcnt52(x), _trail_zcnt64(x));
        }
    }
    if (tk_bench_enabled(tk)) {
        uint64_t x = tk_rng_u64(r);
        TK_BENCH(tk, "popcount64", 1, 0.0, { tk_sink_u64 += mp_popcount64(x += 0x9E3779B9); });
    }
}

/* the software 128/64 divide (the _udiv128 of targets with no hardware or
   __int128 path) against the native one and the bitwise reference */
static void t_udiv128_soft(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    int k;
    for (k = 0; k < 300000; k++) {
        uint64_t d, hi, lo, q1, r1, q2, r2;
        int bits = 1 + (int)tk_rng_range(r, 64);
        d = tk_rng_nbits(r, bits);
        if ((k % 9) == 0) d = ~(uint64_t)0 >> (64 - bits);      /* all ones */
        if ((k % 10) == 0) d = (uint64_t)1 << (bits - 1);       /* power of 2 */
        hi = tk_rng_u64(r) % d;
        lo = tk_rng_u64(r);
        if ((k % 17) == 0) { hi = d - 1; lo = ~(uint64_t)0; }   /* largest quotient */
        q1 = _udiv128(hi, lo, d, &r1);
        q2 = mp_udiv128_soft(hi, lo, d, &r2);
        if (!TK_EQ_U64(tk, q2, q1) || !TK_EQ_U64(tk, r2, r1)) break;
        if ((k & 15) == 0 && !TK_EQ_U64(tk, r2, ref_mod(hi, lo, d))) break;
    }
    if (tk_bench_enabled(tk)) {
        uint64_t d = tk_rng_u64(r) | 1, hi = tk_rng_u64(r) % d, lo = tk_rng_u64(r), rr;
        TK_BENCH(tk, "udiv128 soft", 1, 0.0,
                 { tk_sink_u64 += mp_udiv128_soft(hi, lo += 0x9E37, d, &rr) ^ rr; });
    }
}

/* r is a correct Montgomery product of x and y with R = 2^64:
   r < n  and  r * 2^64 == x * y (mod n) */
static int redc_ok(uint64_t r, uint64_t x, uint64_t y, uint64_t n)
{
    uint64_t lo, hi;
    ref_mul64(x, y, &lo, &hi);
    return r < n && ref_mod(r, 0, n) == ref_mod(hi, lo, n);
}

static void t_modular64(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    int k;
    for (k = 0; k < 100000; k++) {
        uint64_t n = rand_modulus(r, 64), x, y, c, inv, nhat, s[2], got;
        if ((k % 13) == 0) n |= (uint64_t)1 << 63;          /* full-width moduli */
        if ((k % 97) == 0) n = ~(uint64_t)0;                /* 2^64 - 1 */
        x = tk_rng_u64(r) % n; y = tk_rng_u64(r) % n; c = tk_rng_u64(r) % n;
        if ((k % 11) == 0) x = n - 1;
        if ((k % 23) == 0) y = n - 1;
        if ((k % 29) == 0) x = 0;

        inv = multiplicative_inverse(n);
        nhat = (uint64_t)0 - inv;
        if (!TK_EQ_U64(tk, inv * n, 1)) break;
        if (!TK_EQ_U64(tk, (multiplicative_inverse52(n) * n) & 0x000fffffffffffffULL, 1)) break;

        /* add / sub */
        s[0] = x + y; s[1] = (s[0] < x);
        if (!TK_EQ_U64(tk, addmod(x, y, n), ref_mod(s[1], s[0], n))) break;
        if (!TK_EQ_U64(tk, submod(x, y, n), (x >= y) ? x - y : x + (n - y))) break;
        if (!TK_EQ_U64(tk, submod(addmod(x, y, n), y, n), x)) break;
        if (n <= 0xffffffffULL) {
            TK_EQ_U64(tk, addmod32((uint32_t)x, (uint32_t)y, (uint32_t)n), (x + y) % n);
            TK_EQ_U64(tk, submod32((uint32_t)x, (uint32_t)y, (uint32_t)n), (x + n - y) % n);
        }

        /* every Montgomery flavour computes the same reduced product */
        got = mulredc_pos(x, y, n, inv);
        if (!TK_CHECKF(tk, redc_ok(got, x, y, n), "mulredc_pos n=%llx", (unsigned long long)n)) break;
        if (!TK_EQ_U64(tk, mulredc_pos_alt(x, y, n, inv), got)) break;
        if (!TK_EQ_U64(tk, mulredc(x, y, n, nhat), got)) break;
        if (n < ((uint64_t)1 << 63) && !TK_EQ_U64(tk, mulredc63(x, y, n, nhat), got)) break;
        if (n < ((uint64_t)1 << 60)) {                      /* lazy: [0, 2n) */
            uint64_t lz = mulredc60(x, y, n, nhat);
            if (!TK_CHECK(tk, lz < 2 * n)) break;
            if (!TK_EQ_U64(tk, (lz >= n) ? lz - n : lz, got)) break;
        }

        got = sqrredc_pos(x, n, inv);
        if (!TK_CHECKF(tk, redc_ok(got, x, x, n), "sqrredc_pos n=%llx", (unsigned long long)n)) break;
        if (!TK_EQ_U64(tk, sqrredc(x, n, nhat), got)) break;
        if (n < ((uint64_t)1 << 63) && !TK_EQ_U64(tk, sqrredc63(x, n, nhat), got)) break;
        if (n < ((uint64_t)1 << 60)) {
            uint64_t lz = sqrredc60(x, n, nhat);
            if (!TK_CHECK(tk, lz < 2 * n)) break;
            if (!TK_EQ_U64(tk, (lz >= n) ? lz - n : lz, got)) break;
        }

        /* multiply-add: x*y/R + c */
        if (!TK_EQ_U64(tk, mfma64(x, y, c, n, inv), addmod(mulredc_pos(x, y, n, inv), c, n))) break;

        /* into and out of Montgomery form: (x*R) * 1 / R == x */
        {
            uint64_t xm = ref_mod(x, 0, n);
            if (!TK_EQ_U64(tk, mulredc_pos(xm, 1, n, inv), x)) break;
        }
    }
    if (tk_bench_enabled(tk)) {
        uint64_t n = tk_rng_u64(r) | 0x8000000000000001ULL, x = tk_rng_u64(r) % n;
        uint64_t n62 = (n >> 2) | 1, x62 = x % n62;
        uint64_t inv = multiplicative_inverse(n), inv62 = multiplicative_inverse(n62);
        TK_BENCH(tk, "mulredc_pos latency", 1, 0.0, { x = mulredc_pos(x, x, n, inv); });
        TK_BENCH(tk, "mulredc_pos_alt latency", 1, 0.0, { x = mulredc_pos_alt(x, x, n, inv); });
        TK_BENCH(tk, "mulredc latency", 1, 0.0, { x = mulredc(x, x, n, (uint64_t)0 - inv); });
        TK_BENCH(tk, "mulredc63 latency", 1, 0.0, { x62 = mulredc63(x62, x62, n62, (uint64_t)0 - inv62); });
        TK_BENCH(tk, "addmod+submod latency", 2, 0.0, { x = submod(addmod(x, inv, n), n >> 1, n); });
        tk_sink_u64 += x + x62;
    }
}

static void t_bingcd64(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    int k;
    TK_EQ_U64(tk, bingcd64(0, 12345), 12345);
    TK_EQ_U64(tk, bingcd64(12345, 0), 12345);
    TK_EQ_U64(tk, bingcd64(1, ~(uint64_t)0), 1);
    for (k = 0; k < 100000; k++) {
        /* u odd (the routine's contract); v anything, with planted common
           factors and heavy powers of two */
        uint64_t g = tk_rng_nbits(r, 1 + (int)tk_rng_range(r, 20)) | 1;
        uint64_t u = (tk_rng_nbits(r, 1 + (int)tk_rng_range(r, 40)) | 1) * g;
        uint64_t v = tk_rng_nbits(r, 1 + (int)tk_rng_range(r, 40)) * g;
        if (k & 1) v <<= tk_rng_range(r, 4);
        if ((k % 7) == 0) v = u;
        if (!TK_EQ_U64(tk, bingcd64(u, v), ref_gcd(u, v))) break;
        if (!TK_EQ_U64(tk, bin_gcd64(u, v), ref_gcd(u, v))) break;   /* legacy alias */
    }
    if (tk_bench_enabled(tk)) {
        uint64_t u = tk_rng_u64(r) | 1, v = tk_rng_u64(r);
        TK_BENCH(tk, "bingcd64 (64-bit operands)", 1, 0.0,
                 { tk_sink_u64 += bingcd64(u, v += 0x9E3779B97F4A7C15ULL); });
    }
}

static void rand128(tk_rng *r, uint64_t *x)
{
    int bits = (int)tk_rng_range(r, 129);
    x[0] = tk_rng_u64(r);
    x[1] = tk_rng_u64(r);
    if (bits == 0) { x[0] = x[1] = 0; }
    else if (bits <= 64) { x[1] = 0; x[0] >>= (64 - bits); }
    else if (bits < 128) { x[1] >>= (128 - bits); }
}

static void t_arith128(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    int k, i;
    for (k = 0; k < 100000; k++) {
        uint64_t a[2], b[2], z[4], e[4];
        int s = (int)tk_rng_range(r, 128);
        rand128(r, a); rand128(r, b);
        if ((k % 9) == 0) { a[0] = a[1] = ~(uint64_t)0; }
        if ((k % 11) == 0) { b[0] = b[1] = ~(uint64_t)0; }
        if ((k % 13) == 0) { b[0] = a[0]; b[1] = a[1]; }

        if (!TK_EQ_U64(tk, add128(a, b, z), ref_add(e, a, b, 2)) || !TK_CHECK(tk, !ref_cmp(z, e, 2))) break;
        if (!TK_EQ_U64(tk, sub128(a, b, z), ref_sub(e, a, b, 2)) || !TK_CHECK(tk, !ref_cmp(z, e, 2))) break;
        if (!TK_CHECK(tk, cmp128(a, b) == ref_cmp(a, b, 2))) break;

        /* in-place forms: the result may alias either operand */
        z[0] = a[0]; z[1] = a[1]; (void)add128(z, b, z); (void)ref_add(e, a, b, 2);
        if (!TK_CHECK(tk, !ref_cmp(z, e, 2))) break;
        z[0] = b[0]; z[1] = b[1]; (void)sub128(a, z, z); (void)ref_sub(e, a, b, 2);
        if (!TK_CHECK(tk, !ref_cmp(z, e, 2))) break;

        shl128(a, s, z);
        e[0] = a[0]; e[1] = a[1];
        for (i = 0; i < s; i++) ref_shl1(e, 2);
        if (!TK_CHECKF(tk, !ref_cmp(z, e, 2), "shl128 s=%d", s)) break;
        shr128(a, s, z);
        e[0] = a[0]; e[1] = a[1];
        for (i = 0; i < s; i++) ref_shr1(e, 2);
        if (!TK_CHECKF(tk, !ref_cmp(z, e, 2), "shr128 s=%d", s)) break;

        mul128(a, b, z);        ref_mul(e, a, 2, b, 2);
        if (!TK_CHECK(tk, !ref_cmp(z, e, 4))) break;
        sqr128(a, z);           ref_mul(e, a, 2, a, 2);
        if (!TK_CHECK(tk, !ref_cmp(z, e, 4))) break;
        mul128x64(a, b[0], z);  ref_mul(e, a, 2, b, 1);
        if (!TK_CHECK(tk, !ref_cmp(z, e, 3))) break;
    }
    /* shift edge cases */
    {
        uint64_t a[2] = { 0x0123456789abcdefULL, 0xfedcba9876543210ULL }, z[2];
        shl128(a, 0, z);   TK_CHECK(tk, z[0] == a[0] && z[1] == a[1]);
        shr128(a, 0, z);   TK_CHECK(tk, z[0] == a[0] && z[1] == a[1]);
        shl128(a, 64, z);  TK_CHECK(tk, z[0] == 0 && z[1] == a[0]);
        shr128(a, 64, z);  TK_CHECK(tk, z[0] == a[1] && z[1] == 0);
        shl128(a, 127, z); TK_CHECK(tk, z[0] == 0 && z[1] == (a[0] << 63));
        shr128(a, 127, z); TK_CHECK(tk, z[0] == (a[1] >> 63) && z[1] == 0);
    }
    if (tk_bench_enabled(tk)) {
        uint64_t a[2], b[2], z[4];
        a[0] = tk_rng_u64(r); a[1] = tk_rng_u64(r); b[0] = tk_rng_u64(r); b[1] = tk_rng_u64(r);
        TK_BENCH(tk, "mul128 (128x128->256)", 1, 0.0,
                 { a[0] += 0x9E37; mul128(a, b, z); tk_sink_u64 += z[0] ^ z[3]; });
        TK_BENCH(tk, "sqr128 (128->256)", 1, 0.0,
                 { a[0] += 0x9E37; sqr128(a, z); tk_sink_u64 += z[0] ^ z[3]; });
        TK_BENCH(tk, "add128+sub128", 2, 0.0,
                 { (void)add128(a, b, z); (void)sub128(z, b, a); tk_sink_u64 += a[0]; });
    }
}

#if MP_USE_AVX512F && !defined(LIMB1_NO_VECTOR)
#define M52 0x000fffffffffffffULL

static void ld8(__m512i *v, const uint64_t *x) { *v = _mm512_loadu_si512((const void *)x); }
static void st8(uint64_t *x, __m512i v) { _mm512_storeu_si512((void *)x, v); }

/* shared 8-lane 52-bit Montgomery kernels and product halves */
static void t_redc52_x8(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    int k, i;
    for (k = 0; k < 30000; k++) {
        uint64_t n[8], x[8], y[8], inv[8], o[8], lo[8], hi[8];
        __m512i vn, vx, vy, vinv, l, h;
        for (i = 0; i < 8; i++) {
            n[i] = rand_modulus(r, 52);
            if ((k % 7) == 0) n[i] |= (uint64_t)1 << 51;        /* full 52 bits */
            if ((k % 31) == 0) n[i] = M52;                      /* 2^52 - 1 */
            x[i] = tk_rng_u64(r) % n[i]; y[i] = tk_rng_u64(r) % n[i];
            if ((k % 5) == 0) x[i] = n[i] - 1;
            if ((k % 9) == 0) y[i] = (i & 1) ? 0 : n[i] - 1;
            inv[i] = multiplicative_inverse52(n[i]);
        }
        ld8(&vn, n); ld8(&vx, x); ld8(&vy, y); ld8(&vinv, inv);

        st8(o, mulredc52_x8(vx, vy, vn, vinv));
        for (i = 0; i < 8; i++) {
            /* o * 2^52 == x*y (mod n), o < n */
            uint64_t p0, p1;
            ref_mul64(x[i], y[i], &p0, &p1);
            if (!TK_CHECKF(tk, o[i] < n[i] && ref_mod(o[i] >> 12, o[i] << 52, n[i]) == ref_mod(p1, p0, n[i]),
                           "mulredc52_x8 lane %d n=%llx", i, (unsigned long long)n[i])) return;
        }
        st8(lo, sqrredc52_x8(vx, vn, vinv));
        st8(hi, mulredc52_x8(vx, vx, vn, vinv));
        for (i = 0; i < 8; i++) if (!TK_EQ_U64(tk, lo[i], hi[i])) return;

        st8(o, addmod52_x8(vx, vy, vn));
        for (i = 0; i < 8; i++) if (!TK_EQ_U64(tk, o[i], (x[i] + y[i]) % n[i])) return;
        st8(o, submod52_x8(vx, vy, vn));
        for (i = 0; i < 8; i++) if (!TK_EQ_U64(tk, o[i], (x[i] + n[i] - y[i]) % n[i])) return;

        /* 52 x 52 product halves on arbitrary 52-bit operands */
        for (i = 0; i < 8; i++) {
            x[i] = tk_rng_u64(r) & M52; y[i] = tk_rng_u64(r) & M52;
            if ((k % 3) == 0) { x[i] = M52; }
            if ((k % 6) == 0) { y[i] = M52; }
        }
        ld8(&vx, x); ld8(&vy, y);
        mp_mul52lohi_x8(vx, vy, &l, &h); st8(lo, l); st8(hi, h);
        for (i = 0; i < 8; i++) {
            uint64_t p0, p1;
            ref_mul64(x[i], y[i], &p0, &p1);
            /* the low half is only specified mod 2^52 without IFMA */
            if (!TK_EQ_U64(tk, lo[i] & M52, p0 & M52)) return;
            if (!TK_EQ_U64(tk, hi[i], (p1 << 12) | (p0 >> 52))) return;
        }
        st8(lo, mp_mul52lo_x8(vx, vy)); st8(hi, mp_mul52hi_x8(vx, vy));
        for (i = 0; i < 8; i++) {
            uint64_t p0, p1;
            ref_mul64(x[i], y[i], &p0, &p1);
            if (!TK_EQ_U64(tk, lo[i], p0 & M52)) return;
            if (!TK_EQ_U64(tk, hi[i], (p1 << 12) | (p0 >> 52))) return;
        }
    }
    if (tk_bench_enabled(tk)) {
        uint64_t n[8], x[8], inv[8];
        __m512i vn, vx, vinv;
        for (i = 0; i < 8; i++) {
            n[i] = tk_rng_nbits(r, 52) | 1; x[i] = tk_rng_u64(r) % n[i];
            inv[i] = multiplicative_inverse52(n[i]);
        }
        ld8(&vn, n); ld8(&vx, x); ld8(&vinv, inv);
        TK_BENCH(tk, "mulredc52_x8 latency (per lane)", 8, 0.0,
                 { vx = mulredc52_x8(vx, vx, vn, vinv); });
        st8(x, vx); tk_sink_u64 += x[0];
    }
}

#if !defined(LIMB1_SCALAR_ONLY)
/* the 52-bit carry primitives the 104-bit and AVX-ECM code is built from */
static void t_carry52_x8(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    int k, i;
    for (k = 0; k < 50000; k++) {
        uint64_t a[8], b[8], c[8], o[8];
        __m512i va, vb, vc;
        __mmask8 cin = (__mmask8)tk_rng_u64(r), msk = (__mmask8)tk_rng_u64(r), co;
        for (i = 0; i < 8; i++) {
            a[i] = tk_rng_u64(r) & M52; b[i] = tk_rng_u64(r) & M52; c[i] = tk_rng_u64(r) & M52;
            if ((k % 3) == 0 && (i & 1)) a[i] = M52;
            if ((k % 5) == 0 && (i & 2)) b[i] = M52;
            if ((k % 7) == 0 && (i & 4)) a[i] = 0;
            if ((k % 11) == 0) b[i] = a[i];
        }
        ld8(&va, a); ld8(&vb, b); ld8(&vc, c);

#define LANE_C(m, i) (((m) >> (i)) & 1)
        st8(o, _mm512_addsetc_epi52(va, vb, &co));
        for (i = 0; i < 8; i++) {
            uint64_t s = a[i] + b[i];
            if (!TK_EQ_U64(tk, o[i], s & M52) || !TK_EQ_U64(tk, LANE_C(co, i), s >> 52)) return;
        }
        st8(o, _mm512_subsetc_epi52(va, vb, &co));
        for (i = 0; i < 8; i++) {
            if (!TK_EQ_U64(tk, o[i], (a[i] - b[i]) & M52) || !TK_EQ_U64(tk, LANE_C(co, i), a[i] < b[i])) return;
        }
        st8(o, _mm512_adc_epi52(va, cin, vb, &co));
        for (i = 0; i < 8; i++) {
            uint64_t s = a[i] + b[i] + LANE_C(cin, i);
            if (!TK_EQ_U64(tk, o[i], s & M52) || !TK_EQ_U64(tk, LANE_C(co, i), s >> 52)) return;
        }
        st8(o, _mm512_sbb_epi52(va, cin, vb, &co));
        for (i = 0; i < 8; i++) {
            uint64_t bi = b[i] + LANE_C(cin, i);
            if (!TK_EQ_U64(tk, o[i], (a[i] - bi) & M52) || !TK_EQ_U64(tk, LANE_C(co, i), a[i] < bi)) return;
        }
        st8(o, _mm512_addcarry_epi52(va, cin, &co));
        for (i = 0; i < 8; i++) {
            uint64_t s = a[i] + LANE_C(cin, i);
            if (!TK_EQ_U64(tk, o[i], s & M52) || !TK_EQ_U64(tk, LANE_C(co, i), s >> 52)) return;
        }
        /* _mm512_subborrow_epi52: the difference is exact; its borrow-out is
           "a == 0" whatever the borrow-in, so it is only meaningful in lanes
           that had a borrow coming in -- which is how it is used.  Check the
           difference everywhere and the borrow in those lanes. */
        st8(o, _mm512_subborrow_epi52(va, cin, &co));
        for (i = 0; i < 8; i++) {
            if (!TK_EQ_U64(tk, o[i], (a[i] - LANE_C(cin, i)) & M52)) return;
            if (LANE_C(cin, i) && !TK_EQ_U64(tk, LANE_C(co, i), a[i] == 0)) return;
        }

        /* masked forms: selected lanes compute, the rest pass the first
           operand through and report no carry */
        st8(o, _mm512_mask_addsetc_epi52(vc, msk, va, vb, &co));
        for (i = 0; i < 8; i++) {
            uint64_t s = a[i] + b[i];
            if (LANE_C(msk, i)) {
                if (!TK_EQ_U64(tk, o[i], s & M52) || !TK_EQ_U64(tk, LANE_C(co, i), s >> 52)) return;
            } else {
                if (!TK_EQ_U64(tk, o[i], c[i]) || !TK_EQ_U64(tk, LANE_C(co, i), 0)) return;
            }
        }
        st8(o, _mm512_mask_subsetc_epi52(vc, msk, va, vb, &co));
        for (i = 0; i < 8; i++) {
            if (LANE_C(msk, i)) {
                if (!TK_EQ_U64(tk, o[i], (a[i] - b[i]) & M52) || !TK_EQ_U64(tk, LANE_C(co, i), a[i] < b[i])) return;
            } else {
                if (!TK_EQ_U64(tk, o[i], c[i]) || !TK_EQ_U64(tk, LANE_C(co, i), 0)) return;
            }
        }
        st8(o, _mm512_mask_adc_epi52(va, msk, cin, vb, &co));
        for (i = 0; i < 8; i++) {
            if (LANE_C(msk, i)) {
                uint64_t s = a[i] + b[i] + LANE_C(cin, i);
                if (!TK_EQ_U64(tk, o[i], s & M52) || !TK_EQ_U64(tk, LANE_C(co, i), s >> 52)) return;
            } else {
                if (!TK_EQ_U64(tk, o[i], a[i])) return;
            }
        }
        st8(o, _mm512_mask_sbb_epi52(va, msk, cin, vb, &co));
        for (i = 0; i < 8; i++) {
            if (LANE_C(msk, i)) {
                uint64_t bi = b[i] + LANE_C(cin, i);
                if (!TK_EQ_U64(tk, o[i], (a[i] - bi) & M52) || !TK_EQ_U64(tk, LANE_C(co, i), a[i] < bi)) return;
            } else {
                if (!TK_EQ_U64(tk, o[i], a[i]) || !TK_EQ_U64(tk, LANE_C(co, i), 0)) return;
            }
        }
#undef LANE_C

        /* carryprop: move bits >= 52 of lo into hi */
        {
            uint64_t l[8], h[8], ol[8], oh[8];
            __m512i vl, vh, vmask = _mm512_set1_epi64((long long)M52);
            for (i = 0; i < 8; i++) { l[i] = tk_rng_u64(r) >> 4; h[i] = tk_rng_u64(r) & M52; }
            ld8(&vl, l); ld8(&vh, h);
            carryprop(vl, vh, vmask);
            st8(ol, vl); st8(oh, vh);
            for (i = 0; i < 8; i++)
                if (!TK_EQ_U64(tk, ol[i], l[i] & M52) || !TK_EQ_U64(tk, oh[i], h[i] + (l[i] >> 52))) return;
        }
    }
}
#endif /* !LIMB1_SCALAR_ONLY */
#endif /* AVX-512 */

static const tk_test tk__limb_inline_tests[] = {
    { "carry32",        t_carry32,         "fast primitives arith" },
    { "scan32_popcount", t_scan32_popcount, "fast primitives bitscan" },
    { "udiv128_soft",   t_udiv128_soft,    "fast primitives arith" },
    { "modular64",      t_modular64,       "fast limb1 modular" },
    { "bingcd64",       t_bingcd64,        "fast limb1 gcd" },
    { "arith128",       t_arith128,        "fast limb2 arith" },
#if MP_USE_AVX512F && !defined(LIMB1_NO_VECTOR)
    { "redc52_x8",      t_redc52_x8,       "fast limb1 modular vector" },
#if !defined(LIMB1_SCALAR_ONLY)
    { "carry52_x8",     t_carry52_x8,      "fast limb1 vector" },
#endif
#endif
};

const tk_module tk_module_limb_inline = {
    "limb_inline",
    "arith/limb header-only kernels: 32-bit primitives, 64-bit modular, 128-bit arith, 52-bit x8",
    tk__limb_inline_tests,
    (int)(sizeof tk__limb_inline_tests / sizeof tk__limb_inline_tests[0])
};
