/*
tinyrho: Pollard rho / Brent for n < 2^64 and n < 2^128 (scalar) and lists
of n < 2^52 and n < 2^104 (AVX-512, one input per lane).  See tinyrho.h.
The arithmetic is that of microecm_edw.c (64-bit and x8 52-bit) and
tinyecm_edw.c (128-bit and x8 104-bit), copied so this file stands alone.
This file is dual-licensed like microecm.h (FreeBSD or MPL 2.0).
*/
#include "tinyrho.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#if defined(_MSC_VER)
#include <intrin.h>
#endif

#ifndef ALIGNED_MEM
#if defined(_MSC_VER)
#define ALIGNED_MEM __declspec(align(64))
#else
#define ALIGNED_MEM __attribute__((aligned(64)))
#endif
#endif

#if defined(_MSC_VER)
#define TRHO_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define TRHO_INLINE inline __attribute__((always_inline))
#else
#define TRHO_INLINE __inline
#endif

#if defined(__GNUC__) || defined(__clang__)
#define trho_ctz64 __builtin_ctzll
#define trho_clz64 __builtin_clzll
#else
static __inline int trho_ctz64(uint64_t x)
{
    unsigned long i;
    _BitScanForward64(&i, x);
    return (int)i;
}
static __inline int trho_clz64(uint64_t x)
{
    unsigned long i;
    _BitScanReverse64(&i, x);
    return 63 - (int)i;
}
#endif

/* ------------------------------------------------------------- common */
static int trho_bits64(uint64_t x)
{
    return x ? 64 - trho_clz64(x) : 0;
}

/* the default iteration limit: ~20x the mean for a balanced semiprime */
static uint64_t trho_default_maxit(int bits)
{
    if (bits > 120)
        bits = 120;
    return ((uint64_t)32 << (bits / 4)) + 4096;
}

/* steps between gcds (scalar): ~2^(bits/4) / 8, in [8, 1024] (tune.sh:
   within a few % of the best at 32..64 bits) */
#ifndef TRHO_SB_OFF
#define TRHO_SB_OFF 3
#endif
#ifndef TRHO_SB_MAX
#define TRHO_SB_MAX 10
#endif
static uint32_t trho_batch(int bits)
{
    int e = bits / 4 - TRHO_SB_OFF;
    if (e < 3) e = 3;
    if (e > TRHO_SB_MAX) e = TRHO_SB_MAX;
    return (uint32_t)1 << e;
}

/* a^-1 mod 2^64, a odd */
static uint64_t trho_inv64(uint64_t a)
{
    uint64_t x0 = (3 * a) ^ 2;
    uint64_t y = 1 - a * x0;
    uint64_t x1 = x0 * (1 + y);
    y *= y;
    uint64_t x2 = x1 * (1 + y);
    y *= y;
    uint64_t x3 = x2 * (1 + y);
    y *= y;
    return x3 * (1 + y);
}

/* gcd(u, v), u odd */
static uint64_t trho_gcd64(uint64_t u, uint64_t v)
{
    if (v == 0)
        return u;
    v >>= trho_ctz64(v);
    while (u != v) {
        if (u > v) {
            uint64_t t = u;
            u = v;
            v = t;
        }
        v -= u;
        v >>= trho_ctz64(v);
    }
    return u;
}

/* ============================================================ scalar 64 */
/* Montgomery multiplication with the positive inverse, R = 2^64, inputs
   and outputs in [0, n) (microecm_edw.c's uedw_mulredc, C version) */
static TRHO_INLINE uint64_t trho_mulredc(uint64_t x, uint64_t y, uint64_t n, uint64_t inv)
{
#if defined(_MSC_VER)
    uint64_t T_hi;
    uint64_t T_lo = _umul128(x, y, &T_hi);
    uint64_t m = T_lo * inv;
    uint64_t mN_hi = __umulh(m, n);
#else
    __uint128_t prod = (__uint128_t)x * y;
    uint64_t T_hi = (uint64_t)(prod >> 64);
    uint64_t T_lo = (uint64_t)(prod);
    uint64_t m = T_lo * inv;
    uint64_t mN_hi = (uint64_t)(((__uint128_t)m * n) >> 64);
#endif
    uint64_t tmp = T_hi + n - mN_hi;
    uint64_t result = T_hi - mN_hi;
    return (T_hi < mN_hi) ? tmp : result;
}

static TRHO_INLINE uint64_t trho_addmod(uint64_t x, uint64_t y, uint64_t n)
{
    return (x >= n - y) ? x - (n - y) : x + y;
}

static TRHO_INLINE uint64_t trho_submod(uint64_t a, uint64_t b, uint64_t n)
{
    uint64_t r = a - b;
    return (a < b) ? r + n : r;
}

#define F64(v) trho_addmod(trho_mulredc((v), (v), n, inv), c, n)

uint64_t tinyrho64(uint64_t n, uint32_t c0, uint64_t maxit)
{
    uint64_t inv, it = 0;
    uint32_t att, m;
    int bits;
    if ((n & 1) == 0)
        return (n > 2) ? 2 : 1;
    if (n < 9)
        return 1;
    bits = trho_bits64(n);
    if (maxit == 0)
        maxit = trho_default_maxit(bits);
    m = trho_batch(bits);
    inv = trho_inv64(n);
    for (att = 0; it < maxit; att++) {
        uint64_t c = (uint64_t)(c0 + att) % n;
        uint64_t x, y = 2, ys, q = 1, g = 1, r = 1, k, i;
        if (c == 0)
            c = 1;
        do {
            x = y;
            for (i = 0; i < r; i++)
                y = F64(y);
            it += r;
            k = 0;
            do {
                uint64_t lim = (r - k < m) ? r - k : m;
                ys = y;
                for (i = 0; i < lim; i++) {
                    y = F64(y);
                    q = trho_mulredc(q, trho_submod(x, y, n), n, inv);
                }
                it += lim;
                k += lim;
                g = trho_gcd64(n, q);
            } while (k < r && g == 1);
            r *= 2;
        } while (g == 1 && it < maxit);
        if (g == 1)
            break;
        if (g == n) {
            // backtrack: the first term of the last batch with a common factor
            for (i = 0; i < m; i++) {
                ys = F64(ys);
                g = trho_gcd64(n, trho_submod(x, ys, n));
                if (g != 1)
                    break;
            }
            if (g == 1 || g == n)
                continue;           // both factors at once: next c
        }
        return g;
    }
    return 1;
}
#undef F64

/* =========================================================== scalar 128 */
typedef struct { uint64_t lo, hi; } trho_u128;

/* hi:lo = a b + c + d */
static TRHO_INLINE uint64_t trho_mac(uint64_t a, uint64_t b, uint64_t c, uint64_t d, uint64_t *hi)
{
#if defined(__SIZEOF_INT128__)
    unsigned __int128 t = (unsigned __int128)a * b + c + d;
    *hi = (uint64_t)(t >> 64);
    return (uint64_t)t;
#else
    uint64_t h;
    uint64_t l = _umul128(a, b, &h);
    l += c;
    h += (l < c);
    l += d;
    h += (l < d);
    *hi = h;
    return l;
#endif
}

static TRHO_INLINE trho_u128 trho_sub_b(trho_u128 a, trho_u128 b, uint64_t *borrow)
{
    trho_u128 r;
    uint64_t b1 = (uint64_t)(a.lo < b.lo);
    r.lo = a.lo - b.lo;
    r.hi = a.hi - b.hi - b1;
    *borrow = (uint64_t)(a.hi < b.hi) | ((uint64_t)(a.hi == b.hi) & b1);
    return r;
}

static TRHO_INLINE trho_u128 trho_add_c(trho_u128 a, trho_u128 b, uint64_t *carry)
{
    trho_u128 r;
    r.lo = a.lo + b.lo;
    uint64_t c1 = (uint64_t)(r.lo < a.lo);
    uint64_t s = a.hi + b.hi;
    uint64_t c2 = (uint64_t)(s < a.hi);
    r.hi = s + c1;
    *carry = c2 | (uint64_t)(r.hi < s);
    return r;
}

static TRHO_INLINE trho_u128 trho_sel(uint64_t mask, trho_u128 a, trho_u128 b)
{
    trho_u128 r;
    r.lo = (a.lo & mask) | (b.lo & ~mask);
    r.hi = (a.hi & mask) | (b.hi & ~mask);
    return r;
}

/* CIOS Montgomery product, R = 2^128, rho = -1/n mod 2^64, inputs < n
   (tinyecm_edw.c's tedw_mulredc; branch-free corrections) */
static TRHO_INLINE trho_u128 trho_mulredc128(trho_u128 u, trho_u128 v, trho_u128 n, uint64_t rho)
{
    uint64_t t0, t1, t2, cy, m, h, bw;
    t0 = trho_mac(u.lo, v.lo, 0, 0, &cy);
    t1 = trho_mac(u.lo, v.hi, cy, 0, &t2);
    m = t0 * rho;
    (void)trho_mac(m, n.lo, t0, 0, &cy);
    t0 = trho_mac(m, n.hi, t1, cy, &cy);
    t1 = t2 + cy;
    t2 = (uint64_t)(t1 < cy);
    t0 = trho_mac(u.hi, v.lo, t0, 0, &cy);
    t1 = trho_mac(u.hi, v.hi, t1, cy, &h);
    t2 += h;
    m = t0 * rho;
    (void)trho_mac(m, n.lo, t0, 0, &cy);
    t0 = trho_mac(m, n.hi, t1, cy, &cy);
    t1 = t2 + cy;
    t2 = (uint64_t)(t1 < cy);
    trho_u128 r;
    r.lo = t0;
    r.hi = t1;
    trho_u128 d = trho_sub_b(r, n, &bw);
    return trho_sel((uint64_t)0 - (t2 | (bw ^ 1)), d, r);
}

static TRHO_INLINE trho_u128 trho_addmod128(trho_u128 x, trho_u128 y, trho_u128 n)
{
    uint64_t carry, bw;
    trho_u128 s = trho_add_c(x, y, &carry);
    trho_u128 d = trho_sub_b(s, n, &bw);
    return trho_sel((uint64_t)0 - (carry | (bw ^ 1)), d, s);
}

static TRHO_INLINE trho_u128 trho_submod128(trho_u128 x, trho_u128 y, trho_u128 n)
{
    uint64_t bw, carry;
    trho_u128 r = trho_sub_b(x, y, &bw);
    trho_u128 nm;
    nm.lo = n.lo & ((uint64_t)0 - bw);
    nm.hi = n.hi & ((uint64_t)0 - bw);
    return trho_add_c(r, nm, &carry);
}

static TRHO_INLINE int trho_eq128(trho_u128 a, trho_u128 b)
{
    return a.lo == b.lo && a.hi == b.hi;
}

static TRHO_INLINE int trho_gt128(trho_u128 a, trho_u128 b)
{
    return (a.hi > b.hi) || (a.hi == b.hi && a.lo > b.lo);
}

static TRHO_INLINE trho_u128 trho_shr128(trho_u128 a, int s)
{
    if (s >= 64) {
        a.lo = a.hi >> (s - 64);
        a.hi = 0;
    } else if (s > 0) {
        a.lo = (a.lo >> s) | (a.hi << (64 - s));
        a.hi >>= s;
    }
    return a;
}

static TRHO_INLINE int trho_ctz128(trho_u128 a)
{
    return a.lo ? trho_ctz64(a.lo) : 64 + trho_ctz64(a.hi);
}

/* gcd(u, v), u odd */
static trho_u128 trho_gcd128(trho_u128 u, trho_u128 v)
{
    if ((v.lo | v.hi) == 0)
        return u;
    v = trho_shr128(v, trho_ctz128(v));
    while (!trho_eq128(u, v)) {
        if ((u.hi | v.hi) == 0) {
            u.lo = trho_gcd64(u.lo, v.lo);
            return u;
        }
        if (trho_gt128(u, v)) {
            trho_u128 t = u;
            u = v;
            v = t;
        }
        uint64_t bw;
        v = trho_sub_b(v, u, &bw);
        v = trho_shr128(v, trho_ctz128(v));
    }
    return u;
}

static int trho_bits128(trho_u128 a)
{
    return a.hi ? 64 + trho_bits64(a.hi) : trho_bits64(a.lo);
}

#define F128(v) trho_addmod128(trho_mulredc128((v), (v), N, rho), c, N)

int tinyrho128(const uint64_t n[2], uint64_t f[2], uint32_t c0, uint64_t maxit)
{
    trho_u128 N, one;
    uint64_t rho, it = 0;
    uint32_t att, m;
    int bits;
    N.lo = n[0];
    N.hi = n[1];
    f[0] = 1;
    f[1] = 0;
    one.lo = 1;
    one.hi = 0;
    if (N.hi == 0) {
        uint64_t g = tinyrho64(N.lo, c0, maxit);
        f[0] = g;
        return g > 1;
    }
    if ((N.lo & 1) == 0) {
        f[0] = 2;
        return 1;
    }
    bits = trho_bits128(N);
    if (maxit == 0)
        maxit = trho_default_maxit(bits);
    m = trho_batch(bits);
    rho = (uint64_t)0 - trho_inv64(N.lo);
    for (att = 0; it < maxit; att++) {
        trho_u128 c, x, y, ys, q = one, g = one;
        uint64_t r = 1, k, i;
        c.lo = c0 + att;
        c.hi = 0;
        if (c.lo == 0)
            c.lo = 1;
        y.lo = 2;
        y.hi = 0;
        do {
            x = y;
            for (i = 0; i < r; i++)
                y = F128(y);
            it += r;
            k = 0;
            do {
                uint64_t lim = (r - k < m) ? r - k : m;
                ys = y;
                for (i = 0; i < lim; i++) {
                    y = F128(y);
                    q = trho_mulredc128(q, trho_submod128(x, y, N), N, rho);
                }
                it += lim;
                k += lim;
                g = trho_gcd128(N, q);
            } while (k < r && trho_eq128(g, one));
            r *= 2;
        } while (trho_eq128(g, one) && it < maxit);
        if (trho_eq128(g, one))
            break;
        if (trho_eq128(g, N)) {
            for (i = 0; i < m; i++) {
                ys = F128(ys);
                g = trho_gcd128(N, trho_submod128(x, ys, N));
                if (!trho_eq128(g, one))
                    break;
            }
            if (trho_eq128(g, one) || trho_eq128(g, N))
                continue;
        }
        f[0] = g.lo;
        f[1] = g.hi;
        return 1;
    }
    return 0;
}
#undef F128

/* =========================================================== x8 lists */
static uint64_t trho_list_its;

uint64_t tinyrho_list_iterations(void)
{
    uint64_t r = trho_list_its;
    trho_list_its = 0;
    return r;
}

#ifdef USE_AVX512F
#include <immintrin.h>

/* tuning knobs (compile-time overridable for experiments) */
#ifndef TRHO_VB_OFF
#define TRHO_VB_OFF 1
#endif
#ifndef TRHO_VB_MIN
#define TRHO_VB_MIN 6
#endif
#ifndef TRHO_VB_MAX
#define TRHO_VB_MAX 10
#endif

/* steps between gcds (vector lists; also the first x = y point): the gcds
   are per lane and scalar (8 U of them per batch, branchy), so they are
   spaced much further apart than in the scalar code: 2^(bits/4) / 2, in
   [64, 1024] (arith/tinyrho/tune.sh: 1.3-1.8x faster than 2^(bits/4) / 8
   at 36-48 bits, equal at 52) */
static uint32_t trho_vbatch(int bits)
{
    int e = bits / 4 - TRHO_VB_OFF;
    if (e < TRHO_VB_MIN) e = TRHO_VB_MIN;
    if (e > TRHO_VB_MAX) e = TRHO_VB_MAX;
    return (uint32_t)1 << e;
}

#define TRHO_M52 0x000fffffffffffffULL

#ifndef IFMA
static TRHO_INLINE void trho8_lohi(__m512i a, __m512i b, __m512i *l, __m512i *h)
{
    const __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    const __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);
    const __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);
    __m512d a_ld = _mm512_cvtepu64_pd(a);
    __m512d b_ld = _mm512_cvtepu64_pd(b);
    __m512d hd = _mm512_fmadd_round_pd(a_ld, b_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
    *h = _mm512_sub_epi64(_mm512_castpd_si512(hd), vbias1);
    hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), hd);
    *l = _mm512_and_si512(_mm512_castpd_si512(_mm512_fmadd_round_pd(a_ld, b_ld, hd,
        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC))), _mm512_set1_epi64(TRHO_M52));
}
#endif

/* ------------------------------------------------------------- 52-bit */
/* Montgomery multiplication, R = 2^52, positive inverse; inputs in [0, n)
   (microecm_edw.c's uedw8_mulredc) */
static TRHO_INLINE __m512i trho8_mulredc52(__m512i x, __m512i y, __m512i N, __m512i invN)
{
    __m512i T_hi, T_lo, m, mN_hi;
#ifdef IFMA
    const __m512i z = _mm512_setzero_si512();
    T_lo = _mm512_madd52lo_epu64(z, x, y);
    T_hi = _mm512_madd52hi_epu64(z, x, y);
    m = _mm512_madd52lo_epu64(z, T_lo, invN);
    mN_hi = _mm512_madd52hi_epu64(z, m, N);
#else
    __m512i l;
    trho8_lohi(x, y, &T_lo, &T_hi);
    m = _mm512_and_si512(_mm512_mullo_epi64(T_lo, invN), _mm512_set1_epi64(TRHO_M52));
    trho8_lohi(m, N, &l, &mN_hi);
    (void)l;
#endif
    __m512i tmp = _mm512_sub_epi64(_mm512_add_epi64(T_hi, N), mN_hi);
    __m512i result = _mm512_sub_epi64(T_hi, mN_hi);
    __mmask8 msk = _mm512_cmplt_epu64_mask(T_hi, mN_hi);
    return _mm512_mask_mov_epi64(result, msk, tmp);
}

static TRHO_INLINE __m512i trho8_addmod52(__m512i x, __m512i y, __m512i N)
{
    __m512i t = _mm512_add_epi64(x, y);
    __mmask8 m = _mm512_cmpge_epu64_mask(t, N);
    return _mm512_mask_sub_epi64(t, m, t, N);
}

static TRHO_INLINE __m512i trho8_submod52(__m512i x, __m512i y, __m512i N)
{
    __m512i t = _mm512_sub_epi64(x, y);
    __mmask8 m = _mm512_cmpgt_epu64_mask(y, x);
    return _mm512_mask_add_epi64(t, m, t, N);
}

#ifndef TRHO52_U
#define TRHO52_U 4                      /* vectors in flight */
#endif
#define TRHO52_L (8 * TRHO52_U)         /* lanes */

typedef struct {
    ALIGNED_MEM uint64_t n[TRHO52_L], inv[TRHO52_L], c[TRHO52_L];
    ALIGNED_MEM uint64_t y[TRHO52_L], x[TRHO52_L], q[TRHO52_L];
    uint64_t it[TRHO52_L];      // iterations of this input, all tries
    uint64_t s[TRHO52_L];       // iterations of this try
    uint64_t r[TRHO52_L];       // next x = y point (s value)
    int64_t idx[TRHO52_L];      // input index, -1 = idle
    uint32_t att[TRHO52_L];
} trho52_lanes;

static void trho52_try(trho52_lanes *L, int k)
{
    L->c[k] = (uint64_t)(1 + L->att[k]) % L->n[k];
    if (L->c[k] == 0)
        L->c[k] = 1;
    L->y[k] = L->x[k] = 2;
    L->q[k] = 1;
    L->s[k] = 0;
    L->r[k] = 0;                // set by the caller (the batch)
}

/* next input into lane k (inputs not for the vector code are done here) */
static void trho52_fill(trho52_lanes *L, int k, const uint64_t *q64, uint64_t *f64,
    uint32_t num, uint32_t *next, uint32_t G, uint64_t maxit)
{
    while (*next < num) {
        uint32_t i = (*next)++;
        uint64_t n = q64[i];
        if ((n & 1) == 0 || n < 9 || n >= ((uint64_t)1 << 52)) {
            f64[i] = tinyrho64(n, 1, maxit);
            continue;
        }
        L->idx[k] = i;
        L->n[k] = n;
        L->inv[k] = trho_inv64(n) & TRHO_M52;
        L->it[k] = 0;
        L->att[k] = 0;
        trho52_try(L, k);
        L->r[k] = G;
        return;
    }
    // idle lane: a dummy modulus, results ignored
    L->idx[k] = -1;
    L->n[k] = 9;
    L->inv[k] = trho_inv64(9) & TRHO_M52;
    L->c[k] = 1;
    L->y[k] = L->x[k] = 2;
    L->q[k] = 1;
}

void tinyrho52_x8_list(const uint64_t *q64, uint64_t *f64, uint32_t num, uint64_t maxit)
{
    trho52_lanes L;
    uint32_t next = 0, G, i;
    int k, active, j;
    if (num == 0)
        return;
    G = trho_vbatch(trho_bits64(q64[0]));
    for (k = 0; k < TRHO52_L; k++)
        trho52_fill(&L, k, q64, f64, num, &next,
            G, maxit ? maxit : trho_default_maxit(trho_bits64(q64[0])));
    while (1) {
        active = 0;
        for (k = 0; k < TRHO52_L; k++)
            active += (L.idx[k] >= 0);
        if (!active)
            break;
        __m512i n[TRHO52_U], inv[TRHO52_U], c[TRHO52_U], y[TRHO52_U], x[TRHO52_U], q[TRHO52_U];
        for (j = 0; j < TRHO52_U; j++) {
            n[j] = _mm512_load_si512(L.n + 8 * j);
            inv[j] = _mm512_load_si512(L.inv + 8 * j);
            c[j] = _mm512_load_si512(L.c + 8 * j);
            y[j] = _mm512_load_si512(L.y + 8 * j);
            x[j] = _mm512_load_si512(L.x + 8 * j);
            q[j] = _mm512_load_si512(L.q + 8 * j);
        }
        for (i = 0; i < G; i++) {
            for (j = 0; j < TRHO52_U; j++) {
                y[j] = trho8_addmod52(trho8_mulredc52(y[j], y[j], n[j], inv[j]), c[j], n[j]);
                q[j] = trho8_mulredc52(q[j], trho8_submod52(x[j], y[j], n[j]), n[j], inv[j]);
            }
        }
        for (j = 0; j < TRHO52_U; j++) {
            _mm512_store_si512(L.y + 8 * j, y[j]);
            _mm512_store_si512(L.q + 8 * j, q[j]);
        }
        trho_list_its += (uint64_t)G * (uint64_t)active;
        for (k = 0; k < TRHO52_L; k++) {
            int64_t ix = L.idx[k];
            uint64_t lim;
            if (ix < 0)
                continue;
            L.it[k] += G;
            L.s[k] += G;
            uint64_t g = trho_gcd64(L.n[k], L.q[k]);
            lim = maxit ? maxit : trho_default_maxit(trho_bits64(L.n[k]));
            if (g == 1) {
                if (L.s[k] == L.r[k]) {
                    L.x[k] = L.y[k];
                    L.r[k] *= 2;
                }
                if (L.it[k] < lim)
                    continue;
                f64[ix] = 1;
            } else if (g == L.n[k]) {
                // both factors in one batch (or q = 0): the next c
                if (L.it[k] < lim) {
                    L.att[k]++;
                    trho52_try(&L, k);
                    L.r[k] = G;
                    continue;
                }
                f64[ix] = 1;
            } else {
                f64[ix] = g;
            }
            trho52_fill(&L, k, q64, f64, num, &next, G, lim);
        }
    }
}

/* ------------------------------------------------------------ 104-bit */
typedef struct { __m512i d0, d1; } trho_v;

/* 104-bit Montgomery product (R = 2^104, positive inverse invN mod 2^104),
   inputs < N (tinyecm_edw.c's tedw8_mulredc, own final correction) */
static TRHO_INLINE trho_v trho8_mulredc104(trho_v x, trho_v y, trho_v N, trho_v invN)
{
    const __m512i mask = _mm512_set1_epi64(TRHO_M52);
    __m512i T0, T1, T2, T3, m0, m1, mN1, mN2, mN3;
#define CARRY52(lo, hi) { hi = _mm512_add_epi64(hi, _mm512_srli_epi64(lo, 52)); \
    lo = _mm512_and_si512(lo, mask); }
#ifdef IFMA
    const __m512i z = _mm512_setzero_si512();
#define TL(a, b, cc) _mm512_madd52lo_epu64((a), (b), (cc))
#define TH(a, b, cc) _mm512_madd52hi_epu64((a), (b), (cc))
    T0 = TL(z, x.d0, y.d0);
    T1 = TL(TL(TH(z, x.d0, y.d0), x.d0, y.d1), x.d1, y.d0);
    T2 = TL(TH(TH(z, x.d0, y.d1), x.d1, y.d0), x.d1, y.d1);
    T3 = TH(z, x.d1, y.d1);
    CARRY52(T1, T2);
    CARRY52(T2, T3);
    m0 = TL(z, T0, invN.d0);
    m1 = _mm512_and_si512(TL(TL(TH(z, T0, invN.d0), T0, invN.d1), T1, invN.d0), mask);
    mN1 = TL(TL(TH(z, m0, N.d0), m0, N.d1), m1, N.d0);
    mN2 = TL(TH(TH(z, m0, N.d1), m1, N.d0), m1, N.d1);
    mN3 = TH(z, m1, N.d1);
#undef TL
#undef TH
#else
    __m512i l, h, l2, h2, l3, h3, l4, h4, mN0;
    trho8_lohi(x.d0, y.d0, &T0, &T1);
    trho8_lohi(x.d0, y.d1, &l2, &h2);
    trho8_lohi(x.d1, y.d0, &l3, &h3);
    trho8_lohi(x.d1, y.d1, &l4, &h4);
    T1 = _mm512_add_epi64(_mm512_add_epi64(T1, l2), l3);
    T2 = _mm512_add_epi64(_mm512_add_epi64(h2, h3), l4);
    T3 = h4;
    CARRY52(T1, T2);
    CARRY52(T2, T3);
    trho8_lohi(T0, invN.d0, &m0, &m1);
    trho8_lohi(T0, invN.d1, &l, &h);
    m1 = _mm512_add_epi64(m1, l);
    trho8_lohi(T1, invN.d0, &l, &h);
    m1 = _mm512_and_si512(_mm512_add_epi64(m1, l), mask);
    trho8_lohi(m0, N.d0, &mN0, &mN1);
    trho8_lohi(m0, N.d1, &l2, &h2);
    trho8_lohi(m1, N.d0, &l3, &h3);
    trho8_lohi(m1, N.d1, &l4, &h4);
    mN1 = _mm512_add_epi64(_mm512_add_epi64(mN1, l2), l3);
    mN2 = _mm512_add_epi64(_mm512_add_epi64(h2, h3), l4);
    mN3 = h4;
    (void)mN0;
#endif
    CARRY52(mN1, mN2);
    CARRY52(mN2, mN3);
    // (T - mN) / 2^104 = (T3:T2) - (mN3:mN2) (the low halves cancel),
    // + N if negative
    __m512i d2 = _mm512_sub_epi64(T2, mN2);
    __m512i d3 = _mm512_sub_epi64(_mm512_sub_epi64(T3, mN3), _mm512_srli_epi64(d2, 63));
    d2 = _mm512_and_si512(d2, mask);
    __mmask8 neg = _mm512_cmplt_epi64_mask(d3, _mm512_setzero_si512());
    d2 = _mm512_mask_add_epi64(d2, neg, d2, N.d0);
    d3 = _mm512_mask_add_epi64(d3, neg, d3, _mm512_add_epi64(N.d1, _mm512_srli_epi64(d2, 52)));
    trho_v r;
    r.d0 = _mm512_and_si512(d2, mask);
    r.d1 = _mm512_and_si512(d3, mask);
#undef CARRY52
    return r;
}

/* x + y mod N, inputs < N < 2^104 */
static TRHO_INLINE trho_v trho8_addmod104(trho_v x, trho_v y, trho_v N)
{
    const __m512i mask = _mm512_set1_epi64(TRHO_M52);
    __m512i s0 = _mm512_add_epi64(x.d0, y.d0);
    __m512i s1 = _mm512_add_epi64(_mm512_add_epi64(x.d1, y.d1), _mm512_srli_epi64(s0, 52));
    s0 = _mm512_and_si512(s0, mask);
    // s - N, keep it if not negative
    __m512i d0 = _mm512_sub_epi64(s0, N.d0);
    __m512i d1 = _mm512_sub_epi64(_mm512_sub_epi64(s1, N.d1), _mm512_srli_epi64(d0, 63));
    __mmask8 ok = _mm512_cmpge_epi64_mask(d1, _mm512_setzero_si512());
    trho_v r;
    r.d0 = _mm512_mask_mov_epi64(s0, ok, _mm512_and_si512(d0, mask));
    r.d1 = _mm512_mask_mov_epi64(s1, ok, d1);
    return r;
}

/* x - y mod N, inputs < N < 2^104 */
static TRHO_INLINE trho_v trho8_submod104(trho_v x, trho_v y, trho_v N)
{
    const __m512i mask = _mm512_set1_epi64(TRHO_M52);
    __m512i d0 = _mm512_sub_epi64(x.d0, y.d0);
    __m512i d1 = _mm512_sub_epi64(_mm512_sub_epi64(x.d1, y.d1), _mm512_srli_epi64(d0, 63));
    d0 = _mm512_and_si512(d0, mask);
    __mmask8 neg = _mm512_cmplt_epi64_mask(d1, _mm512_setzero_si512());
    d0 = _mm512_mask_add_epi64(d0, neg, d0, N.d0);
    d1 = _mm512_mask_add_epi64(d1, neg, d1, _mm512_add_epi64(N.d1, _mm512_srli_epi64(d0, 52)));
    trho_v r;
    r.d0 = _mm512_and_si512(d0, mask);
    r.d1 = _mm512_and_si512(d1, mask);
    return r;
}

#ifndef TRHO104_U
#define TRHO104_U 2
#endif
#define TRHO104_L (8 * TRHO104_U)

typedef struct {
    ALIGNED_MEM uint64_t n[2][TRHO104_L], inv[2][TRHO104_L], c[2][TRHO104_L];
    ALIGNED_MEM uint64_t y[2][TRHO104_L], x[2][TRHO104_L], q[2][TRHO104_L];
    uint64_t it[TRHO104_L], s[TRHO104_L], r[TRHO104_L];
    int64_t idx[TRHO104_L];
    uint32_t att[TRHO104_L];
} trho104_lanes;

static trho_u128 trho_lane128(uint64_t (*a)[TRHO104_L], int k)
{
    trho_u128 v;
    v.lo = a[0][k] | (a[1][k] << 52);
    v.hi = a[1][k] >> 12;
    return v;
}

/* n^-1 mod 2^104, n odd, as two 52-bit digits */
static void trho_inv104(trho_u128 n, uint64_t *i0, uint64_t *i1)
{
    // x = n^-1 mod 2^64, one Newton step to 2^128: x (2 - n x)
    uint64_t x = trho_inv64(n.lo), hi, t_hi;
    uint64_t t_lo = trho_mac(n.lo, x, 0, 0, &t_hi);     // n.lo x
    t_hi += n.hi * x;                                   // n x mod 2^128 = 1 + 2^64 t_hi
    // x (2 - n x) = x (1 - 2^64 t_hi) mod 2^128 = x - 2^64 (x t_hi)
    (void)t_lo;
    hi = (uint64_t)0 - x * t_hi;
    *i0 = x & TRHO_M52;
    *i1 = ((x >> 52) | (hi << 12)) & TRHO_M52;
}

static void trho104_try(trho104_lanes *L, int k)
{
    L->c[0][k] = 1 + L->att[k];
    L->c[1][k] = 0;
    L->y[0][k] = L->x[0][k] = 2;
    L->y[1][k] = L->x[1][k] = 0;
    L->q[0][k] = 1;
    L->q[1][k] = 0;
    L->s[k] = 0;
}

static void trho104_fill(trho104_lanes *L, int k, const uint64_t *q64, uint64_t *f64,
    uint32_t num, uint32_t *next, uint32_t G, uint64_t maxit)
{
    while (*next < num) {
        uint32_t i = (*next)++;
        uint64_t n0 = q64[2 * i], n1 = q64[2 * i + 1];
        trho_u128 n;
        n.lo = n0 | (n1 << 52);
        n.hi = n1 >> 12;
        if (n.hi == 0 && n.lo < ((uint64_t)1 << 60)) {
            // small (or even) input: the scalar code
            uint64_t nn[2] = { n.lo, n.hi };
            tinyrho128(nn, f64 + 2 * i, 1, maxit);
            continue;
        }
        if ((n0 & 1) == 0 || (n1 >> 52) != 0) {
            uint64_t nn[2] = { n.lo, n.hi };
            tinyrho128(nn, f64 + 2 * i, 1, maxit);
            continue;
        }
        L->idx[k] = i;
        L->n[0][k] = n0;
        L->n[1][k] = n1;
        trho_inv104(n, &L->inv[0][k], &L->inv[1][k]);
        L->it[k] = 0;
        L->att[k] = 0;
        trho104_try(L, k);
        L->r[k] = G;
        return;
    }
    L->idx[k] = -1;
    L->n[0][k] = 9;
    L->n[1][k] = 0;
    {
        trho_u128 nine = { 9, 0 };
        trho_inv104(nine, &L->inv[0][k], &L->inv[1][k]);
    }
    L->c[0][k] = 1;
    L->c[1][k] = 0;
    L->y[0][k] = L->x[0][k] = 2;
    L->y[1][k] = L->x[1][k] = 0;
    L->q[0][k] = 1;
    L->q[1][k] = 0;
}

static trho_v trho104_load(uint64_t (*a)[TRHO104_L], int j)
{
    trho_v v;
    v.d0 = _mm512_load_si512(a[0] + 8 * j);
    v.d1 = _mm512_load_si512(a[1] + 8 * j);
    return v;
}

static void trho104_store(uint64_t (*a)[TRHO104_L], int j, trho_v v)
{
    _mm512_store_si512(a[0] + 8 * j, v.d0);
    _mm512_store_si512(a[1] + 8 * j, v.d1);
}

void tinyrho104_x8_list(const uint64_t *q64, uint64_t *f64, uint32_t num, uint64_t maxit)
{
    trho104_lanes L;
    uint32_t next = 0, G, i;
    int k, active, j, bits;
    if (num == 0)
        return;
    for (i = 0; i < num; i++) {
        f64[2 * i] = 1;
        f64[2 * i + 1] = 0;
    }
    {
        trho_u128 n0;
        n0.lo = q64[0] | (q64[1] << 52);
        n0.hi = q64[1] >> 12;
        bits = trho_bits128(n0);
    }
    G = trho_vbatch(bits);
    for (k = 0; k < TRHO104_L; k++)
        trho104_fill(&L, k, q64, f64, num, &next, G, maxit);
    while (1) {
        active = 0;
        for (k = 0; k < TRHO104_L; k++)
            active += (L.idx[k] >= 0);
        if (!active)
            break;
        trho_v n[TRHO104_U], inv[TRHO104_U], c[TRHO104_U], y[TRHO104_U], x[TRHO104_U], q[TRHO104_U];
        for (j = 0; j < TRHO104_U; j++) {
            n[j] = trho104_load(L.n, j);
            inv[j] = trho104_load(L.inv, j);
            c[j] = trho104_load(L.c, j);
            y[j] = trho104_load(L.y, j);
            x[j] = trho104_load(L.x, j);
            q[j] = trho104_load(L.q, j);
        }
        for (i = 0; i < G; i++) {
            for (j = 0; j < TRHO104_U; j++) {
                y[j] = trho8_addmod104(trho8_mulredc104(y[j], y[j], n[j], inv[j]), c[j], n[j]);
                q[j] = trho8_mulredc104(q[j], trho8_submod104(x[j], y[j], n[j]), n[j], inv[j]);
            }
        }
        for (j = 0; j < TRHO104_U; j++) {
            trho104_store(L.y, j, y[j]);
            trho104_store(L.q, j, q[j]);
        }
        trho_list_its += (uint64_t)G * (uint64_t)active;
        for (k = 0; k < TRHO104_L; k++) {
            int64_t ix = L.idx[k];
            if (ix < 0)
                continue;
            trho_u128 nk = trho_lane128(L.n, k);
            uint64_t lim = maxit ? maxit : trho_default_maxit(trho_bits128(nk));
            L.it[k] += G;
            L.s[k] += G;
            trho_u128 g = trho_gcd128(nk, trho_lane128(L.q, k));
            if (g.lo == 1 && g.hi == 0) {
                if (L.s[k] == L.r[k]) {
                    L.x[0][k] = L.y[0][k];
                    L.x[1][k] = L.y[1][k];
                    L.r[k] *= 2;
                }
                if (L.it[k] < lim)
                    continue;
            } else if (trho_eq128(g, nk)) {
                if (L.it[k] < lim) {
                    L.att[k]++;
                    trho104_try(&L, k);
                    L.r[k] = G;
                    continue;
                }
            } else {
                f64[2 * ix] = g.lo;
                f64[2 * ix + 1] = g.hi;
            }
            trho104_fill(&L, k, q64, f64, num, &next, G, maxit);
        }
    }
}

#else /* !USE_AVX512F: the scalar code per input */

void tinyrho52_x8_list(const uint64_t *q64, uint64_t *f64, uint32_t num, uint64_t maxit)
{
    uint32_t i;
    for (i = 0; i < num; i++)
        f64[i] = tinyrho64(q64[i], 1, maxit);
}

void tinyrho104_x8_list(const uint64_t *q64, uint64_t *f64, uint32_t num, uint64_t maxit)
{
    uint32_t i;
    for (i = 0; i < num; i++) {
        uint64_t n[2];
        n[0] = q64[2 * i] | (q64[2 * i + 1] << 52);
        n[1] = q64[2 * i + 1] >> 12;
        tinyrho128(n, f64 + 2 * i, 1, maxit);
    }
}

#endif /* USE_AVX512F */
