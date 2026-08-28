/*----------------------------------------------------------------------
 YAFU modular test system -- shared test data & self-contained oracles.
 See test_data.h. Public domain.
----------------------------------------------------------------------*/
#include "test_data.h"
#include <string.h>

/* ================================================================== *
 * 64-bit modular arithmetic helpers (no GMP, no compiler intrinsics
 * beyond __int128 where available).
 * ================================================================== */
static uint64_t tk__mulmod64(uint64_t a, uint64_t b, uint64_t m)
{
#if defined(__SIZEOF_INT128__)
    return (uint64_t)(((__uint128_t)a * b) % m);
#else
    /* portable double-and-add; m < 2^63 assumed on this path is NOT safe,
     * so do the general binary version. */
    uint64_t r = 0;
    a %= m;
    while (b) {
        if (b & 1) { r += a; if (r < a || r >= m) r -= m; }
        a <<= 1; if (a < (a >> 1) /*ovf*/ || a >= m) a -= m;
        b >>= 1;
    }
    return r;
#endif
}

static uint64_t tk__powmod64(uint64_t a, uint64_t e, uint64_t m)
{
    uint64_t r = 1 % m;
    a %= m;
    while (e) {
        if (e & 1) r = tk__mulmod64(r, a, m);
        a = tk__mulmod64(a, a, m);
        e >>= 1;
    }
    return r;
}

/* one strong-probable-prime (Miller-Rabin) round to base a */
static int tk__sprp(uint64_t n, uint64_t a)
{
    uint64_t d = n - 1, x;
    int r = 0, i;
    if ((a % n) == 0) return 1;
    while ((d & 1) == 0) { d >>= 1; r++; }
    x = tk__powmod64(a, d, n);
    if (x == 1 || x == n - 1) return 1;
    for (i = 1; i < r; i++) {
        x = tk__mulmod64(x, x, n);
        if (x == n - 1) return 1;
    }
    return 0;
}

int tk_is_prime_u64(uint64_t n)
{
    /* deterministic base set: sufficient for all n < 2^64
       (Sinclair / best-known set). */
    static const uint64_t bases[] = {
        2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL,
        23ULL, 29ULL, 31ULL, 37ULL
    };
    size_t i;
    if (n < 2) return 0;
    for (i = 0; i < sizeof bases / sizeof bases[0]; i++) {
        if (n == bases[i]) return 1;
        if (n % bases[i] == 0) return 0;
    }
    for (i = 0; i < sizeof bases / sizeof bases[0]; i++)
        if (!tk__sprp(n, bases[i])) return 0;
    return 1;
}

/* ================================================================== *
 * Generators
 * ================================================================== */
uint64_t tk_gen_prime_u64(tk_ctx *tk, int bits)
{
    tk_rng *r = tk_rng_of(tk);
    for (;;) {
        uint64_t x = tk_rng_nbits(r, bits) | 1ULL;   /* odd, exact length */
        if (tk_is_prime_u64(x)) return x;
    }
}

uint64_t tk_gen_semiprime_u64(tk_ctx *tk, int bits, uint64_t *p, uint64_t *q)
{
    int pb = bits / 2;
    int qb = bits - pb;
    uint64_t pp, qq, t;
    if (pb < 2) pb = 2;
    if (qb < 2) qb = 2;
    for (;;) {
        pp = tk_gen_prime_u64(tk, pb);
        qq = tk_gen_prime_u64(tk, qb);
        if (pp == qq) continue;
        /* guard against overflow of the product past 64 bits */
        if (qq != 0 && pp > (uint64_t)0xFFFFFFFFFFFFFFFFULL / qq) continue;
        break;
    }
    if (pp > qq) { t = pp; pp = qq; qq = t; }
    if (p) *p = pp;
    if (q) *q = qq;
    return pp * qq;
}

/* ================================================================== *
 * Corpora (verified at authoring time with sympy; re-checked at runtime).
 * ================================================================== */
const char *const tk_primes_dec[] = {
    "97",
    "1000003",
    "2147483647",                                       /* 2^31-1  */
    "2305843009213693951",                              /* 2^61-1  */
    "618970019642690137449562111",                      /* 2^89-1  */
    "170141183460469231731687303715884105727",          /* 2^127-1 */
    "10000000000000000000000000000000000043",           /* ~10^37  */
    "1606938044258990275541962092341162602522202993782792835301611", /* ~2^200 */
    NULL
};

const char *const tk_composites_dec[] = {
    "18446744073709551616",                             /* 2^64            */
    "147573952589676412927",                            /* 2^67-1          */
    "3215031751",                                        /* spsp 2,3,5,7    */
    "170141183460469231731687303715884105728",          /* 2^127           */
    "10000000000000000000000000000603000000000000000000000000001881", /* RSA-like */
    NULL
};

const uint64_t tk_spsp2_u64[] = {
    2047ULL, 3277ULL, 4033ULL, 4681ULL, 8321ULL, 15841ULL, 29341ULL,
    42799ULL, 49141ULL, 52633ULL, 65281ULL, 74665ULL, 80581ULL,
    85489ULL, 88357ULL, 90751ULL
};
const int tk_spsp2_u64_count = (int)(sizeof tk_spsp2_u64 / sizeof tk_spsp2_u64[0]);

const uint64_t tk_carmichael_u64[] = {
    561ULL, 1105ULL, 1729ULL, 2465ULL, 2821ULL, 6601ULL, 8911ULL,
    10585ULL, 15841ULL, 29341ULL, 41041ULL, 46657ULL, 52633ULL,
    62745ULL, 63973ULL, 75361ULL
};
const int tk_carmichael_u64_count =
    (int)(sizeof tk_carmichael_u64 / sizeof tk_carmichael_u64[0]);

const tk_semiprime64 tk_semiprimes_u64[] = {
    { 227761879448843ULL,      14569531ULL,   15632753ULL   },  /* 48-bit */
    { 31703017636450211ULL,    135274121ULL,  234361291ULL  },  /* 55-bit */
    { 2437069920619449697ULL,  1240906309ULL, 1963943533ULL },  /* 62-bit */
    { 8720161814859239021ULL,  2320360421ULL, 3758106601ULL }   /* 63-bit */
};
const int tk_semiprimes_u64_count =
    (int)(sizeof tk_semiprimes_u64 / sizeof tk_semiprimes_u64[0]);

/* ================================================================== *
 * Layer -1: corpus self-check. Validates the data itself with the
 * oracle so a mistaken edit is caught before it can mask a real bug.
 * ================================================================== */
static void t_spsp2(tk_ctx *tk)
{
    int i;
    for (i = 0; i < tk_spsp2_u64_count; i++) {
        uint64_t n = tk_spsp2_u64[i];
        /* must be composite ... */
        TK_CHECKF(tk, !tk_is_prime_u64(n), "spsp2[%d]=%llu should be composite",
                  i, (unsigned long long)n);
        /* ... yet pass a Fermat/MR test to base 2 (that's what makes it a trap) */
        TK_CHECKF(tk, tk__sprp(n, 2), "spsp2[%d]=%llu should be sprp base 2",
                  i, (unsigned long long)n);
    }
}

static void t_carmichael(tk_ctx *tk)
{
    int i;
    for (i = 0; i < tk_carmichael_u64_count; i++) {
        uint64_t n = tk_carmichael_u64[i];
        TK_CHECKF(tk, !tk_is_prime_u64(n), "carmichael[%d]=%llu should be composite",
                  i, (unsigned long long)n);
    }
}

static void t_semiprimes(tk_ctx *tk)
{
    int i;
    for (i = 0; i < tk_semiprimes_u64_count; i++) {
        const tk_semiprime64 *s = &tk_semiprimes_u64[i];
        TK_CHECKF(tk, s->p <= s->q, "semiprime[%d]: p<=q ordering", i);
        TK_CHECKF(tk, tk_is_prime_u64(s->p), "semiprime[%d]: p=%llu prime",
                  i, (unsigned long long)s->p);
        TK_CHECKF(tk, tk_is_prime_u64(s->q), "semiprime[%d]: q=%llu prime",
                  i, (unsigned long long)s->q);
        TK_CHECKF(tk, s->p * s->q == s->n, "semiprime[%d]: p*q == N", i);
    }
}

static void t_oracle_sanity(tk_ctx *tk)
{
    /* spot-check the oracle against tiny hand-known values */
    static const uint64_t small_primes[] = {2,3,5,7,11,13,9973,104729};
    static const uint64_t small_comps[]  = {1,4,6,9,15,21,100,104730,999999};
    size_t i;
    for (i = 0; i < sizeof small_primes/sizeof small_primes[0]; i++)
        TK_CHECKF(tk, tk_is_prime_u64(small_primes[i]),
                  "%llu should be prime", (unsigned long long)small_primes[i]);
    for (i = 0; i < sizeof small_comps/sizeof small_comps[0]; i++)
        TK_CHECKF(tk, !tk_is_prime_u64(small_comps[i]),
                  "%llu should be composite", (unsigned long long)small_comps[i]);
    TK_CHECK(tk, !tk_is_prime_u64(0));
    /* the largest 64-bit prime */
    TK_CHECK(tk, tk_is_prime_u64(18446744073709551557ULL));
}

static void t_generators(tk_ctx *tk)
{
    int i;
    for (i = 0; i < 200; i++) {
        uint64_t p, q, n;
        int bits = 20 + (int)tk_rng_range(tk_rng_of(tk), 40); /* 20..59 */
        n = tk_gen_semiprime_u64(tk, bits, &p, &q);
        TK_CHECKF(tk, tk_is_prime_u64(p) && tk_is_prime_u64(q),
                  "generated factors prime (n=%llu)", (unsigned long long)n);
        TK_CHECKF(tk, p != q, "generated p != q");
        TK_CHECKF(tk, p * q == n, "generated p*q == n");
        TK_CHECKF(tk, !tk_is_prime_u64(n), "generated semiprime is composite");
    }
}

static const tk_test tk__selfcheck_tests[] = {
    { "oracle_sanity", t_oracle_sanity, "fast selfcheck" },
    { "spsp2_corpus",  t_spsp2,         "fast selfcheck" },
    { "carmichael",    t_carmichael,    "fast selfcheck" },
    { "semiprimes",    t_semiprimes,    "fast selfcheck" },
    { "generators",    t_generators,    "fast selfcheck" }
};

const tk_module tk_module_selfcheck = {
    "selfcheck",
    "validate the KAT corpus and generators with the built-in oracle",
    tk__selfcheck_tests,
    (int)(sizeof tk__selfcheck_tests / sizeof tk__selfcheck_tests[0])
};
