/*----------------------------------------------------------------------
 YAFU modular test system -- entry point.

 Lists the modules to run and hands off to the testkit runner. Layer 1
 (which links against YAFU + GMP) is compiled in unless TK_NO_LAYER1 is
 defined, so the framework/Layer-0 subset can be built and run stand-alone.
 Public domain.
----------------------------------------------------------------------*/
#include "testkit.h"
#include "test_data.h"

/* ---- a tiny framework self-test (all checks pass; keeps exit status 0) ---- */
static void t_meta_rng(tk_ctx *tk)
{
    tk_rng a, b;
    int i;
    tk_rng_seed(&a, 42);
    tk_rng_seed(&b, 42);
    for (i = 0; i < 1000; i++)
        TK_EQ_U64(tk, tk_rng_u64(&a), tk_rng_u64(&b));   /* determinism */
    /* range bound */
    for (i = 0; i < 10000; i++) {
        uint64_t n = 1 + tk_rng_u64(&a) % 1000;
        TK_CHECK(tk, tk_rng_range(&a, n) < n);
    }
    /* exact bit length */
    for (i = 1; i <= 64; i++) {
        uint64_t x = tk_rng_nbits(&a, i);
        int hb = 0; uint64_t t = x;
        while (t) { hb++; t >>= 1; }
        TK_CHECKF(tk, hb == i, "nbits(%d) produced %d-bit value", i, hb);
    }
}

static void t_meta_assert(tk_ctx *tk)
{
    /* exercise the assertion API on values known to pass */
    TK_CHECK(tk, 1);
    TK_CHECK(tk, 2 + 2 == 4);
    TK_EQ_U64(tk, 0xFFFFFFFFFFFFFFFFULL, ~0ULL);
    TK_CHECKF(tk, sizeof(uint64_t) == 8, "uint64_t is %zu bytes", sizeof(uint64_t));
    TK_REQUIRE(tk, 1, "this REQUIRE must not fire");
}

static const tk_test tk__meta_tests[] = {
    { "rng",    t_meta_rng,    "fast meta" },
    { "assert", t_meta_assert, "fast meta" }
};
static const tk_module tk_module_meta = {
    "meta", "framework self-test (RNG determinism, assertion mechanics)",
    tk__meta_tests, (int)(sizeof tk__meta_tests / sizeof tk__meta_tests[0])
};

/* ---- Layer 0 ---- */
extern const tk_module tk_module_mp_arith;
extern const tk_module tk_module_mp_bitscan;

/* ---- Layer 1 (YAFU + GMP) ---- */
#ifndef TK_NO_LAYER1
extern const tk_module tk_module_sp_arith;
extern const tk_module tk_module_modular;
extern const tk_module tk_module_primality;
extern const tk_module tk_module_sieve;
extern const tk_module tk_module_microecm;
extern const tk_module tk_module_tinyecm;
#ifdef TK_WITH_LAYER3
extern const tk_module tk_module_siqs;   /* Layer 3: links the full factoring archives */
#endif
#endif

static const tk_module *const modules[] = {
    &tk_module_meta,
    &tk_module_selfcheck,
    &tk_module_mp_arith,
    &tk_module_mp_bitscan,
#ifndef TK_NO_LAYER1
    &tk_module_sp_arith,
    &tk_module_modular,
    &tk_module_primality,
    &tk_module_sieve,
    &tk_module_microecm,
    &tk_module_tinyecm,
#ifdef TK_WITH_LAYER3
    &tk_module_siqs,
#endif
#endif
};

int main(int argc, char **argv)
{
    return tk_run(modules, (int)(sizeof modules / sizeof modules[0]), argc, argv);
}
