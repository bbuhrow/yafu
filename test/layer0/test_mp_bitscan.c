/*----------------------------------------------------------------------
 Layer 0 -- mp_platform.h bit-scan primitives.

 Tests _lead_zcnt64 / _trail_zcnt64 (leading/trailing zero COUNT for
 nonzero inputs), the zero-safe wrappers _lead_full_zcnt / _trail_full_zcnt
 (return 64 at zero), and _reset_lsb64 (== x & (x-1)).

 References are bit-serial loops (obviously correct). Inputs mix a
 systematic single-bit / mask sweep with random values, and are kept
 NONZERO for the raw count entry points (zero is undefined there on the
 builtin paths). The reset-lsb check uses a scratch copy and verifies both
 exact x&(x-1) equality and the structural invariant that iterating to zero
 takes exactly popcount(x) steps, clearing bits low-to-high.

 Re-homed mp_bitscan_test harness. Public domain.

 NOTE: _reset_lsb64 has historically had an inconsistent contract across
 build paths -- value-returning on the BMI2 _blsr_u64 branch but a mutating
 `((x)&=((x)-1))` macro on the fallback (evaluates its argument twice and
 needs a modifiable lvalue). The recommended fix is to make every path a
 value-returning inline; until/unless that lands, the tk_reset_lsb64()
 wrapper below is written to be correct under BOTH conventions, so this
 module compiles and passes either way.
----------------------------------------------------------------------*/
#include "testkit.h"
#include "mp_platform.h"

/* Adapt _reset_lsb64 to a value-returning helper, independent of whether the
 * header defines it as the mutating fallback macro `((x)&=((x)-1))`, a
 * value-returning macro, or the BMI2 _blsr_u64 intrinsic: `x` here is a local
 * parameter, so `return _reset_lsb64(x);` yields x&(x-1) under all three (the
 * mutating form writes the local once, so no lvalue error, no double-write). */
static uint64_t tk_reset_lsb64(uint64_t x) { return _reset_lsb64(x); }

/* ---- references ---- */
static int ref_clz64(uint64_t x)          /* count leading zeros, x != 0 */
{
    int n = 0;
    while ((x & 0x8000000000000000ULL) == 0) { n++; x <<= 1; }
    return n;
}
static int ref_ctz64(uint64_t x)          /* count trailing zeros, x != 0 */
{
    int n = 0;
    while ((x & 1ULL) == 0) { n++; x >>= 1; }
    return n;
}
static int ref_popcount64(uint64_t x)
{
    int n = 0;
    while (x) { n += (int)(x & 1); x >>= 1; }
    return n;
}

/* nonzero test vectors: single bits, high masks, low masks, two-bit, random */
static void gen_nonzero(tk_ctx *tk, uint64_t *v, int *nv, int cap)
{
    tk_rng *r = tk_rng_of(tk);
    int b, n = 0;
    for (b = 0; b < 64 && n < cap; b++) v[n++] = (1ULL << b);          /* single bit  */
    for (b = 0; b < 64 && n < cap; b++) v[n++] = (~0ULL << b);         /* high run    */
    for (b = 0; b < 63 && n < cap; b++) v[n++] = (~0ULL >> b);         /* low run     */
    while (n < cap) {
        uint64_t x = tk_rng_u64(r);
        if (x == 0) continue;
        x >>= (tk_rng_u64(r) % 64);   /* spread the top bit position    */
        if (x == 0) x = 1;
        v[n++] = x;
    }
    *nv = n;
}

/* ================================================================== */
static void t_lead_zcnt(tk_ctx *tk)
{
    uint64_t v[512]; int nv, i;
    gen_nonzero(tk, v, &nv, 512);
    for (i = 0; i < nv; i++)
        TK_EQ_U64(tk, _lead_zcnt64(v[i]), (uint64_t)ref_clz64(v[i]));

    if (tk_bench_enabled(tk)) {
        tk_rng *r = tk_rng_of(tk);
        uint64_t x = tk_rng_u64(r) | 1;
        TK_BENCH(tk, "lead_zcnt64 ref", 1, 0.0,
                 { tk_sink_u64 += (uint64_t)ref_clz64((x += 0x9E37) | 1); });
        TK_BENCH(tk, "lead_zcnt64 native", 1, 0.0,
                 { tk_sink_u64 += _lead_zcnt64((x += 0x9E37) | 1); });
    }
}

static void t_trail_zcnt(tk_ctx *tk)
{
    uint64_t v[512]; int nv, i;
    gen_nonzero(tk, v, &nv, 512);
    for (i = 0; i < nv; i++)
        TK_EQ_U64(tk, _trail_zcnt64(v[i]), (uint64_t)ref_ctz64(v[i]));

    if (tk_bench_enabled(tk)) {
        tk_rng *r = tk_rng_of(tk);
        uint64_t x = tk_rng_u64(r) | 1;
        TK_BENCH(tk, "trail_zcnt64 ref", 1, 0.0,
                 { tk_sink_u64 += (uint64_t)ref_ctz64((x += 0x9E37) | 1); });
        TK_BENCH(tk, "trail_zcnt64 native", 1, 0.0,
                 { tk_sink_u64 += _trail_zcnt64((x += 0x9E37) | 1); });
    }
}

static void t_full_zcnt_zero_safe(tk_ctx *tk)
{
    /* the zero-safe wrappers must agree with the raw ones for nonzero
       inputs and return 64 at zero. */
    uint64_t v[256]; int nv, i;
    gen_nonzero(tk, v, &nv, 256);
    for (i = 0; i < nv; i++) {
        TK_EQ_U64(tk, _lead_full_zcnt(v[i]),  _lead_zcnt64(v[i]));
        TK_EQ_U64(tk, _trail_full_zcnt(v[i]), _trail_zcnt64(v[i]));
    }
    TK_EQ_U64(tk, _lead_full_zcnt(0),  64);
    TK_EQ_U64(tk, _trail_full_zcnt(0), 64);
}

static void t_reset_lsb(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    int k;
    /* exact equality vs x & (x-1) */
    for (k = 0; k < 100000; k++) {
        uint64_t x = tk_rng_u64(r);
        TK_EQ_U64(tk, tk_reset_lsb64(x), x & (x - 1));
    }
    /* structural: iterating reset_lsb to zero clears exactly popcount(x)
       bits, lowest first. */
    for (k = 0; k < 20000; k++) {
        uint64_t x = tk_rng_u64(r);
        int steps = 0, want = ref_popcount64(x);
        uint64_t cur = x, prev_low = 0;
        int first = 1;
        while (cur != 0) {
            uint64_t low = cur & (~cur + 1);            /* lowest set bit */
            if (!first) TK_CHECK(tk, low > prev_low);   /* ascending */
            prev_low = low; first = 0;
            cur = tk_reset_lsb64(cur);
            steps++;
            if (steps > 64) { TK_CHECK(tk, 0); break; }
        }
        TK_EQ_U64(tk, (uint64_t)steps, (uint64_t)want);
    }

    if (tk_bench_enabled(tk)) {
        uint64_t x = tk_rng_u64(r);
        /* latency-shaped: dependent chain, like real bit-iteration loops.
           Re-set the top bit each step so the chain never drains to zero. */
        TK_BENCH(tk, "reset_lsb64 latency", 1, 0.0,
                 { x = tk_reset_lsb64(x | 0x8000000000000000ULL); tk_sink_u64 += x; });
    }
}

static const tk_test tk__mp_bitscan_tests[] = {
    { "lead_zcnt",       t_lead_zcnt,           "fast primitives bitscan" },
    { "trail_zcnt",      t_trail_zcnt,          "fast primitives bitscan" },
    { "full_zcnt_zero",  t_full_zcnt_zero_safe, "fast primitives bitscan" },
    { "reset_lsb",       t_reset_lsb,           "fast primitives bitscan" }
};

const tk_module tk_module_mp_bitscan = {
    "mp_bitscan",
    "mp_platform leading/trailing zero-count and reset-lsb",
    tk__mp_bitscan_tests,
    (int)(sizeof tk__mp_bitscan_tests / sizeof tk__mp_bitscan_tests[0])
};
