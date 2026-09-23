/*----------------------------------------------------------------------
 Layer 0 -- mp_platform.h arithmetic primitives.

 Tests _umul128, _udiv128, mp_addcarry_u64, mp_subborrow_u64 against
 self-contained portable-C reference oracles (no intrinsics, no __int128
 in the references), plus contract/reconstruction invariants. Optional
 micro-benchmarks report ns/op vs the reference.

 This is the mp_arith_test harness re-homed as a testkit module. It
 #includes the project's real mp_platform.h; nothing here changes.
 Public domain.
----------------------------------------------------------------------*/
#include "testkit.h"
#include "mp_platform.h"

/* ---- portable reference oracles (deliberately intrinsic-free) ---- */

static uint64_t ref_umul128(uint64_t a, uint64_t b, uint64_t *hi)
{
    uint64_t al = a & 0xFFFFFFFFULL, ah = a >> 32;
    uint64_t bl = b & 0xFFFFFFFFULL, bh = b >> 32;
    uint64_t ll = al * bl;
    uint64_t lh = al * bh;
    uint64_t hl = ah * bl;
    uint64_t hh = ah * bh;
    uint64_t cross = (ll >> 32) + (lh & 0xFFFFFFFFULL) + (hl & 0xFFFFFFFFULL);
    uint64_t lo = (ll & 0xFFFFFFFFULL) | (cross << 32);
    *hi = hh + (lh >> 32) + (hl >> 32) + (cross >> 32);
    return lo;
}

/* 128/64 -> 64 restoring division. Requires hi < d (quotient fits in 64). */
static uint64_t ref_udiv128(uint64_t hi, uint64_t lo, uint64_t d, uint64_t *rem)
{
    uint64_t q = 0, r = hi;
    int i;
    for (i = 63; i >= 0; i--) {
        /* r = (r << 1) | bit i of lo ; watch for r's top bit */
        uint64_t top = r >> 63;
        r = (r << 1) | ((lo >> i) & 1ULL);
        if (top || r >= d) { r -= d; q |= (1ULL << i); }
    }
    *rem = r;
    return q;
}

static unsigned char ref_addcarry(unsigned char cin, uint64_t a, uint64_t b, uint64_t *out)
{
    uint64_t s = a + b;
    unsigned char c = (s < a);
    s += cin;
    c = (unsigned char)(c | (s < (uint64_t)cin));
    *out = s;
    return c;
}

static unsigned char ref_subborrow(unsigned char bin, uint64_t a, uint64_t b, uint64_t *out)
{
    uint64_t d = a - b;
    unsigned char bo = (a < b);
    uint64_t d2 = d - bin;
    bo = (unsigned char)(bo | (d < (uint64_t)bin));
    *out = d2;
    return bo;
}

/* ---- interesting corner inputs ---- */
static const uint64_t corner[] = {
    0ULL, 1ULL, 2ULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFULL, 0x100000000ULL,
    0x8000000000000000ULL, 0x7FFFFFFFFFFFFFFFULL, 0xDEADBEEFCAFEBABEULL,
    0xFFFFFFFF00000000ULL, 0x00000000FFFFFFFFULL
};
#define NCORNER (int)(sizeof corner / sizeof corner[0])

/* ================================================================== */
static void t_umul128(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    int i, j, k;
    for (i = 0; i < NCORNER; i++)
        for (j = 0; j < NCORNER; j++) {
            uint64_t rh, ph, rl, pl;
            rl = ref_umul128(corner[i], corner[j], &rh);
            pl = _umul128(corner[i], corner[j], &ph);
            TK_EQ_U64(tk, pl, rl);
            TK_EQ_U64(tk, ph, rh);
        }
    for (k = 0; k < 200000; k++) {
        uint64_t a = tk_rng_u64(r), b = tk_rng_u64(r), rh, ph, rl, pl;
        rl = ref_umul128(a, b, &rh);
        pl = _umul128(a, b, &ph);
        if (!TK_EQ_U64(tk, pl, rl)) break;
        if (!TK_EQ_U64(tk, ph, rh)) break;
    }

    if (tk_bench_enabled(tk)) {
        uint64_t a = tk_rng_u64(r) | 1, b = tk_rng_u64(r) | 1;
        TK_BENCH(tk, "umul128 ref", 1, 0.0,
                 { uint64_t h; tk_sink_u64 += ref_umul128(a += 0x9E37, b, &h) ^ h; });
        TK_BENCH(tk, "umul128 native", 1, 0.0,
                 { uint64_t h; tk_sink_u64 += _umul128(a += 0x9E37, b, &h) ^ h; });
    }
}

static void t_udiv128(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    int k;
    for (k = 0; k < 200000; k++) {
        uint64_t d = tk_rng_u64(r); uint64_t hi, lo, rq, rr, pq, pr;
        if (d == 0) d = 1;
        hi = tk_rng_u64(r) % d;        /* precondition: hi < d */
        lo = tk_rng_u64(r);
        rq = ref_udiv128(hi, lo, d, &rr);
        pq = _udiv128(hi, lo, d, &pr);
        if (!TK_EQ_U64(tk, pq, rq)) break;
        if (!TK_EQ_U64(tk, pr, rr)) break;
        /* independent reconstruction: hi:lo == q*d + r */
        {
            uint64_t mh, ml = _umul128(pq, d, &mh);
            uint64_t sl, carry;
            carry = mp_addcarry_u64(0, ml, pr, &sl);
            TK_EQ_U64(tk, sl, lo);
            TK_EQ_U64(tk, mh + carry, hi);
        }
    }
    if (tk_bench_enabled(tk)) {
        uint64_t d = (tk_rng_u64(r) | 1);
        uint64_t hi = tk_rng_u64(r) % d, lo = tk_rng_u64(r), rr;
        TK_BENCH(tk, "udiv128 ref", 1, 0.0,
                 { tk_sink_u64 += ref_udiv128(hi, lo += 0x9E37, d, &rr) ^ rr; });
        TK_BENCH(tk, "udiv128 native", 1, 0.0,
                 { tk_sink_u64 += _udiv128(hi, lo += 0x9E37, d, &rr) ^ rr; });
    }
}

#define NCHAIN 64
static void t_addcarry(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    int k, i;
    for (k = 0; k < 100000; k++) {
        unsigned char cin = (unsigned char)(tk_rng_u64(r) & 1);
        uint64_t a = tk_rng_u64(r), b = tk_rng_u64(r), ro, po;
        unsigned char rc = ref_addcarry(cin, a, b, &ro);
        unsigned char pc = mp_addcarry_u64(cin, a, b, &po);
        if (!TK_EQ_U64(tk, po, ro)) break;
        if (!TK_EQ_U64(tk, pc, rc)) break;
    }
    /* multi-limb chain equivalence */
    for (k = 0; k < 20000; k++) {
        uint64_t x[NCHAIN], y[NCHAIN], zr[NCHAIN], zp[NCHAIN];
        unsigned char cr = 0, cp = 0;
        for (i = 0; i < NCHAIN; i++) { x[i] = tk_rng_u64(r); y[i] = tk_rng_u64(r); }
        for (i = 0; i < NCHAIN; i++) cr = ref_addcarry(cr, x[i], y[i], &zr[i]);
        for (i = 0; i < NCHAIN; i++) cp = mp_addcarry_u64(cp, x[i], y[i], &zp[i]);
        TK_EQ_U64(tk, cp, cr);
        for (i = 0; i < NCHAIN; i++)
            if (!TK_EQ_U64(tk, zp[i], zr[i])) break;
    }
    if (tk_bench_enabled(tk)) {
        uint64_t x[NCHAIN], y[NCHAIN], z[NCHAIN];
        for (i = 0; i < NCHAIN; i++) { x[i] = tk_rng_u64(r); y[i] = tk_rng_u64(r); }
        TK_BENCH(tk, "addcarry ref (per-limb)", NCHAIN, 0.0,
                 { unsigned char c = 0; for (i = 0; i < NCHAIN; i++)
                   c = ref_addcarry(c, x[i], y[i], &z[i]); tk_sink_u64 += z[0] + c; });
        TK_BENCH(tk, "addcarry native (per-limb)", NCHAIN, 0.0,
                 { unsigned char c = 0; for (i = 0; i < NCHAIN; i++)
                   c = mp_addcarry_u64(c, x[i], y[i], &z[i]); tk_sink_u64 += z[0] + c; });
    }
}

static void t_subborrow(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    int k, i;
    for (k = 0; k < 100000; k++) {
        unsigned char bin = (unsigned char)(tk_rng_u64(r) & 1);
        uint64_t a = tk_rng_u64(r), b = tk_rng_u64(r), ro, po;
        unsigned char rb = ref_subborrow(bin, a, b, &ro);
        unsigned char pb = mp_subborrow_u64(bin, a, b, &po);
        if (!TK_EQ_U64(tk, po, ro)) break;
        if (!TK_EQ_U64(tk, pb, rb)) break;
    }
    for (k = 0; k < 20000; k++) {
        uint64_t x[NCHAIN], y[NCHAIN], zr[NCHAIN], zp[NCHAIN];
        unsigned char br = 0, bp = 0;
        for (i = 0; i < NCHAIN; i++) { x[i] = tk_rng_u64(r); y[i] = tk_rng_u64(r); }
        for (i = 0; i < NCHAIN; i++) br = ref_subborrow(br, x[i], y[i], &zr[i]);
        for (i = 0; i < NCHAIN; i++) bp = mp_subborrow_u64(bp, x[i], y[i], &zp[i]);
        TK_EQ_U64(tk, bp, br);
        for (i = 0; i < NCHAIN; i++)
            if (!TK_EQ_U64(tk, zp[i], zr[i])) break;
    }
    if (tk_bench_enabled(tk)) {
        uint64_t x[NCHAIN], y[NCHAIN], z[NCHAIN];
        for (i = 0; i < NCHAIN; i++) { x[i] = tk_rng_u64(r); y[i] = tk_rng_u64(r); }
        TK_BENCH(tk, "subborrow ref (per-limb)", NCHAIN, 0.0,
                 { unsigned char c = 0; for (i = 0; i < NCHAIN; i++)
                   c = ref_subborrow(c, x[i], y[i], &z[i]); tk_sink_u64 += z[0] + c; });
        TK_BENCH(tk, "subborrow native (per-limb)", NCHAIN, 0.0,
                 { unsigned char c = 0; for (i = 0; i < NCHAIN; i++)
                   c = mp_subborrow_u64(c, x[i], y[i], &z[i]); tk_sink_u64 += z[0] + c; });
    }
}

static const tk_test tk__mp_arith_tests[] = {
    { "umul128",   t_umul128,   "fast primitives arith" },
    { "udiv128",   t_udiv128,   "fast primitives arith" },
    { "addcarry",  t_addcarry,  "fast primitives arith" },
    { "subborrow", t_subborrow, "fast primitives arith" }
};

const tk_module tk_module_mp_arith = {
    "mp_arith",
    "mp_platform 128-bit mul/div and add/sub carry chains",
    tk__mp_arith_tests,
    (int)(sizeof tk__mp_arith_tests / sizeof tk__mp_arith_tests[0])
};
