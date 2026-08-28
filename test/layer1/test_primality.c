/*----------------------------------------------------------------------
 Layer 1 -- primality: APRCL / BPSW / MR (mpz_aprcl.h) and the small-int
 tinyprp kernels (tinyprp.h), checked against the known-answer corpus in
 test_data and the built-in exact uint64 oracle.

 What each layer of the check buys us:
   * mpz_aprcl(N) is a proof -> must return PRP_PRIME on primes,
     PRP_COMPOSITE on composites, with no probabilistic slack.
   * mpz_bpsw_prp / mpz_strongbpsw_prp have no known counterexample, so
     they must agree with truth on the whole corpus -- including the
     strong-pseudoprime-base-2 and Carmichael traps, which are exactly the
     inputs a naive Fermat/MR-base-2 test gets wrong.
   * mpz_sprp(n,2) SHOULD be fooled by the spsp2 corpus (returns PRP_PRP on
     composites) -- we assert that trap behaviour so a regression that
     accidentally "fixes" it by weakening the base handling is caught.
   * tinyprp's uint64/128 kernels are checked against tk_is_prime_u64
     (exact < 2^64): no false negatives ever, and bpsw_prp_128x1 must be
     exact in the 64-bit range.

 APRCL return codes (mpz_aprcl.h): PRP_COMPOSITE 0, PRP_PRP 1, PRP_PRIME 2.

 NOTE (verify on first build): fermat_prp_64x1 is assumed to be a base-2
 Fermat test returning 1 for probable-prime, 0 for composite. bpsw_prp_128x1
 takes a uint64[2] (little-endian lo,hi) and returns 1/0. Adjust the two
 spots marked [PRP-CONV] if your build differs.
 Public domain.
----------------------------------------------------------------------*/
#include "testkit.h"
#include "test_data.h"
#include <gmp.h>
#include "mpz_aprcl.h"
#include "tinyprp.h"

static void set_mpz_u128(mpz_t out, const uint64_t v[2])
{
    mpz_set_ui(out, v[1]);
    mpz_mul_2exp(out, out, 64);
    mpz_add_ui(out, out, v[0]);
}
static void get_mpz_u128(uint64_t v[2], const mpz_t in)
{
    mpz_t t; mpz_init(t);
    v[0] = mpz_get_ui(in) & 0xFFFFFFFFULL;               /* build 64b safely */
    mpz_tdiv_q_2exp(t, in, 32); v[0] |= (mpz_get_ui(t) & 0xFFFFFFFFULL) << 32;
    mpz_tdiv_q_2exp(t, in, 64);
    v[1] = mpz_get_ui(t) & 0xFFFFFFFFULL;
    mpz_tdiv_q_2exp(t, in, 96); v[1] |= (mpz_get_ui(t) & 0xFFFFFFFFULL) << 32;
    mpz_clear(t);
}

/* ================================================================== *
 * mpz APRCL / BPSW / MR against the decimal-string corpus
 * ================================================================== */
static void t_aprcl_known(tk_ctx *tk)
{
    mpz_t n;
    int i;
    mpz_init(n);

    for (i = 0; tk_primes_dec[i]; i++) {
        TK_REQUIRE(tk, mpz_set_str(n, tk_primes_dec[i], 10) == 0,
                   "bad decimal literal: %s", tk_primes_dec[i]);
        TK_CHECKF(tk, mpz_aprcl(n) == PRP_PRIME,
                  "APRCL should prove prime: %s", tk_primes_dec[i]);
        TK_CHECKF(tk, mpz_bpsw_prp(n) != PRP_COMPOSITE,
                  "BPSW should accept prime: %s", tk_primes_dec[i]);
    }
    for (i = 0; tk_composites_dec[i]; i++) {
        TK_REQUIRE(tk, mpz_set_str(n, tk_composites_dec[i], 10) == 0,
                   "bad decimal literal: %s", tk_composites_dec[i]);
        TK_CHECKF(tk, mpz_aprcl(n) == PRP_COMPOSITE,
                  "APRCL should prove composite: %s", tk_composites_dec[i]);
        TK_CHECKF(tk, mpz_bpsw_prp(n) == PRP_COMPOSITE,
                  "BPSW should reject composite: %s", tk_composites_dec[i]);
    }
    mpz_clear(n);
}

static void t_bpsw_traps(tk_ctx *tk)
{
    /* strong pseudoprimes to base 2: MR-base-2 is fooled, BPSW is not */
    mpz_t n, two;
    int i;
    mpz_inits(n, two, NULL);
    mpz_set_ui(two, 2);

    for (i = 0; i < tk_spsp2_u64_count; i++) {
        mpz_set_ui(n, (unsigned long)(tk_spsp2_u64[i] >> 32));
        mpz_mul_2exp(n, n, 32);
        mpz_add_ui(n, n, (unsigned long)(tk_spsp2_u64[i] & 0xFFFFFFFFUL));
        TK_CHECKF(tk, mpz_sprp(n, two) == PRP_PRP,
                  "spsp2[%d] should fool MR base 2", i);          /* the trap  */
        TK_CHECKF(tk, mpz_bpsw_prp(n) == PRP_COMPOSITE,
                  "spsp2[%d] must not fool BPSW", i);             /* the point */
        TK_CHECKF(tk, mpz_aprcl(n) == PRP_COMPOSITE,
                  "spsp2[%d] must not fool APRCL", i);
    }
    /* Carmichael numbers: Fermat to many bases is fooled; BPSW is not */
    for (i = 0; i < tk_carmichael_u64_count; i++) {
        mpz_set_ui(n, (unsigned long)(tk_carmichael_u64[i] >> 32));
        mpz_mul_2exp(n, n, 32);
        mpz_add_ui(n, n, (unsigned long)(tk_carmichael_u64[i] & 0xFFFFFFFFUL));
        TK_CHECKF(tk, mpz_bpsw_prp(n) == PRP_COMPOSITE,
                  "carmichael[%d] must not fool BPSW", i);
    }
    mpz_clears(n, two, NULL);
}

/* ================================================================== *
 * Fermat base 2 uint64/128 kernels vs the exact oracle
 * ================================================================== */
static void t_fermat_u64(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    int i, k;

    /* random values: Fermat has no false negatives, and agrees with the
       oracle except on rare pseudoprimes (which the corpus covers). */
    for (k = 0; k < 100000; k++) {
        uint64_t n = tk_rng_u64(r) | 1ULL;         /* odd, > 1 almost surely */
        int fp, truth;
        if (n < 3) n = 3;
        fp = fermat_prp_64x1(n);                    /* [PRP-CONV] base-2 Fermat */
        truth = tk_is_prime_u64(n);
        if (truth) TK_CHECKF(tk, fp, "%llu prime but Fermat says composite",
                             (unsigned long long)n);
        if (!fp)   TK_CHECKF(tk, !truth, "%llu Fermat-composite but oracle prime",
                             (unsigned long long)n);
    }
    /* the spsp2 corpus is exactly where base-2 Fermat is fooled */
    for (i = 0; i < tk_spsp2_u64_count; i++)
        TK_CHECKF(tk, fermat_prp_64x1(tk_spsp2_u64[i]),
                  "spsp2[%d] should fool 64-bit Fermat base 2", i);
}

static void t_fermat_u128(tk_ctx* tk)
{
    tk_rng* r = tk_rng_of(tk);
    mpz_t n;
    uint64_t n128[2];
    gmp_randstate_t gmp_rng_state;

    mpz_init(n);
    /* GMP's own randomness, seeded from this test's tk stream so runs remain
       reproducible (for a given GMP version) from the printed seed. */
    gmp_randinit_default(gmp_rng_state);
    gmp_randseed_ui(gmp_rng_state, (unsigned long)tk_rng_u64(r));

    int i, k;

    /* random values: Fermat has no false negatives, and agrees with the
       oracle except on rare pseudoprimes (which the corpus covers). */
    for (k = 0; k < 100000; k++) {
        mpz_urandomb(n, gmp_rng_state, 128);

        get_mpz_u128(n128, n);

        int fp, truth;
        fp = fermat_prp_128x1(n128);
        truth = mpz_probab_prime_p(n, 20);

        if (truth) TK_CHECKF(tk, fp, "%llu prime but Fermat says composite",
            (unsigned long long)n);
        if (!fp)   TK_CHECKF(tk, !truth, "%llu Fermat-composite but oracle prime",
            (unsigned long long)n);
    }
    /* the spsp2 corpus should fool Fermat */
    for (i = 0; i < tk_spsp2_u64_count; i++)
    {
        n128[0] = tk_spsp2_u64[i];
        n128[1] = 0;
        TK_CHECKF(tk, fermat_prp_128x1(n128),
            "spsp2[%d] should fool 128-bit Fermat base 2", i);
    }

    gmp_randclear(gmp_rng_state);
    mpz_clear(n);
}

#ifdef USE_AVX512F
static void t_fermat_prp_52x8(tk_ctx* tk)
{
    tk_rng* r = tk_rng_of(tk);
    int i, k;

    /* random values: Fermat has no false negatives, and agrees with the
       oracle except on rare pseudoprimes (which the corpus covers). */
    for (k = 0; k < 100000; k+=8) {
        uint64_t n[8];
        
        for (i = 0; i < 8; i++)
        {
            n[i] = (tk_rng_u64(r) >> 12) | 1ULL;         /* odd, > 1 almost surely */
            if (n[i] < 3) n[i] = 3;
        }

        uint8_t fpmsk;
        int truth;
        
        fpmsk = fermat_prp_52x8(n);                    /* [PRP-CONV] base-2 Fermat */

        for (i = 0; i < 8; i++)
        {
            truth = tk_is_prime_u64(n[i]);
            if (truth) TK_CHECKF(tk, (fpmsk & (1 << i)), "%llu prime but Fermat says composite",
                (unsigned long long)n[i]);
            if (!(fpmsk & (1 << i)))   TK_CHECKF(tk, !truth, "%llu Fermat-composite but oracle prime",
                (unsigned long long)n[i]);
        }
    }
}

static void t_fermat_prp_104x8(tk_ctx* tk)
{
    tk_rng* r = tk_rng_of(tk);
    mpz_t n[8];
    uint64_t n128[16];
    gmp_randstate_t gmp_rng_state;
    int i, k;

    for (i = 0; i < 8; i++)
        mpz_init(n[i]);

    /* GMP's own randomness, seeded from this test's tk stream so runs remain
       reproducible (for a given GMP version) from the printed seed. */
    gmp_randinit_default(gmp_rng_state);
    gmp_randseed_ui(gmp_rng_state, (unsigned long)tk_rng_u64(r));

    /* random values: Fermat has no false negatives, and agrees with the
       oracle except on rare pseudoprimes (which the corpus covers). */
    for (k = 0; k < 100000; k += 8) {

        for (i = 0; i < 8; i++)
        {
            uint64_t ntmp[2];

            mpz_urandomb(n[i], gmp_rng_state, 104);
            mpz_setbit(n[i], 0);
            get_mpz_u128(ntmp, n[i]);

            n128[i    ] = (ntmp[0] & 0xFFFFFFFFFFFFFull);
            n128[i + 8] = (ntmp[0] >> 52) | ((ntmp[1] & 0xFFFFFFFFFFull) << 12);
        }

        uint8_t fpmsk;
        int truth;

        fpmsk = fermat_prp_104x8(n128);

        for (i = 0; i < 8; i++)
        {
            truth = mpz_probab_prime_p(n[i], 20);
            if (truth) TK_CHECKF(tk, (fpmsk & (1 << i)), "%013llx%013llx prime but Fermat says composite",
                (unsigned long long)n128[i + 8], (unsigned long long)n128[i]);
            if (!(fpmsk & (1 << i)))   TK_CHECKF(tk, !truth, "%013llx%013llx Fermat-composite but oracle prime",
                (unsigned long long)n128[i + 8], (unsigned long long)n128[i]);
        }
    }

    gmp_randclear(gmp_rng_state);
    for (i = 0; i < 8; i++)
        mpz_clear(n[i]);
}

#endif

/* ================================================================== *
 * MR base 2 strong PRP uint64/128 kernels vs the exact oracle
 * ================================================================== */
static void t_MR_2sprp_u64(tk_ctx* tk)
{
    tk_rng* r = tk_rng_of(tk);
    int i, k;

    /* random values: Fermat has no false negatives, and agrees with the
       oracle except on rare pseudoprimes (which the corpus covers). */
    for (k = 0; k < 100000; k++) {
        uint64_t n = tk_rng_u64(r) | 1ULL;         /* odd, > 1 almost surely */
        int fp, truth;
        if (n < 3) n = 3;
        fp = MR_2sprp_64x1(n);                    /* [PRP-CONV] base-2 Fermat */
        truth = tk_is_prime_u64(n);
        if (truth) TK_CHECKF(tk, fp, "%llu prime but Fermat says composite",
            (unsigned long long)n);
        if (!fp)   TK_CHECKF(tk, !truth, "%llu Fermat-composite but oracle prime",
            (unsigned long long)n);
    }
    /* the spsp2 corpus should not fool MR */
    for (i = 0; i < tk_spsp2_u64_count; i++)
    {
        TK_CHECKF(tk, MR_2sprp_64x1(tk_spsp2_u64[i]), 
            "%llu should not fool 64-bit MR sprp base 2", tk_spsp2_u64[i]);
    }
}

static void t_MR_2sprp_u128(tk_ctx* tk)
{
    tk_rng* r = tk_rng_of(tk);
    mpz_t n;
    uint64_t n128[2];
    gmp_randstate_t gmp_rng_state;

    mpz_init(n);
    /* GMP's own randomness, seeded from this test's tk stream so runs remain
       reproducible (for a given GMP version) from the printed seed. */
    gmp_randinit_default(gmp_rng_state);
    gmp_randseed_ui(gmp_rng_state, (unsigned long)tk_rng_u64(r));

    int i, k;

    /* random values: Fermat has no false negatives, and agrees with the
       oracle except on rare pseudoprimes (which the corpus covers). */
    for (k = 0; k < 100000; k++) {
        mpz_urandomb(n, gmp_rng_state, 128);

        get_mpz_u128(n128, n);

        int fp, truth;
        fp = MR_2sprp_128x1(n128);
        truth = mpz_probab_prime_p(n, 20);

        if (truth) TK_CHECKF(tk, fp, "%llu prime but MR says composite",
            (unsigned long long)n);
        if (!fp)   TK_CHECKF(tk, !truth, "%llu MR-composite but oracle prime",
            (unsigned long long)n);
    }
    /* the spsp2 corpus should not fool MR */
    for (i = 0; i < tk_spsp2_u64_count; i++)
    {
        n128[0] = tk_spsp2_u64[i];
        n128[1] = 0;
        TK_CHECKF(tk, MR_2sprp_128x1(n128),
            "%llu should not fool 128-bit MR sprp base 2", tk_spsp2_u64[i]);
    }

    gmp_randclear(gmp_rng_state);
    mpz_clear(n);
}

#ifdef USE_AVX512F
static void t_MR_2sprp_52x8(tk_ctx* tk)
{
    tk_rng* r = tk_rng_of(tk);
    int i, k;

    /* random values: Fermat has no false negatives, and agrees with the
       oracle except on rare pseudoprimes (which the corpus covers). */
    for (k = 0; k < 100000; k += 8) {
        uint64_t n[8];

        for (i = 0; i < 8; i++)
        {
            n[i] = (tk_rng_u64(r) >> 12) | 1ULL;         /* odd, > 1 almost surely */
            if (n[i] < 3) n[i] = 3;
        }

        uint8_t fpmsk;
        int truth;

        fpmsk = MR_2sprp_52x8(n);                    /* [PRP-CONV] base-2 Fermat */

        for (i = 0; i < 8; i++)
        {
            truth = tk_is_prime_u64(n[i]);
            if (truth) TK_CHECKF(tk, (fpmsk & (1 << i)), "%llu prime but MR says composite",
                (unsigned long long)n[i]);
            if (!(fpmsk & (1 << i)))   TK_CHECKF(tk, !truth, "%llu MR-composite but oracle prime",
                (unsigned long long)n[i]);
        }
    }
}

static void t_MR_2sprp_104x8(tk_ctx* tk)
{
    tk_rng* r = tk_rng_of(tk);
    mpz_t n[8];
    uint64_t n128[16];
    gmp_randstate_t gmp_rng_state;
    int i, k;

    for (i = 0; i < 8; i++)
        mpz_init(n[i]);

    /* GMP's own randomness, seeded from this test's tk stream so runs remain
       reproducible (for a given GMP version) from the printed seed. */
    gmp_randinit_default(gmp_rng_state);
    gmp_randseed_ui(gmp_rng_state, (unsigned long)tk_rng_u64(r));

    /* random values: Fermat has no false negatives, and agrees with the
       oracle except on rare pseudoprimes (which the corpus covers). */
    for (k = 0; k < 100000; k += 8) {

        for (i = 0; i < 8; i++)
        {
            uint64_t ntmp[2];

            mpz_urandomb(n[i], gmp_rng_state, 104);
            mpz_setbit(n[i], 0);
            get_mpz_u128(ntmp, n[i]);

            n128[i  ] = (ntmp[0] & 0xFFFFFFFFFFFFFull);
            n128[i+8] = (ntmp[0] >> 52) | ((ntmp[1] & 0xFFFFFFFFFFull) << 12);
        }

        uint8_t fpmsk;
        int truth;

        fpmsk = MR_2sprp_104x8(n128);

        for (i = 0; i < 8; i++)
        {
            truth = mpz_probab_prime_p(n[i], 20);
            if (truth) TK_CHECKF(tk, (fpmsk & (1 << i)), "%013llx%013llx prime but MR says composite",
                (unsigned long long)n128[i + 8], (unsigned long long)n128[i]);
            if (!(fpmsk & (1 << i)))   TK_CHECKF(tk, !truth, "%013llx%013llx MR-composite but oracle prime",
                (unsigned long long)n128[i + 8], (unsigned long long)n128[i]);
        }
    }

    gmp_randclear(gmp_rng_state);
    for (i = 0; i < 8; i++)
        mpz_clear(n[i]);
}

#endif

/* ================================================================== *
 * declare module test routines and tags: a list of what gets run.
 * ================================================================== */
static const tk_test tk__primality_tests[] = {
    { "aprcl_known",  t_aprcl_known,  "slow gmp primality" },
    { "bpsw_traps",   t_bpsw_traps,   "slow gmp primality" },
    { "fermat_u64",   t_fermat_u64,   "fast primality fermat_u64" },
    { "MR_2sprp_u64", t_MR_2sprp_u64, "fast primality MR_2sprp_u64" },
    { "fermat_u128",   t_fermat_u128,   "fast primality fermat_u128" },
    { "MR_2sprp_u128", t_MR_2sprp_u128, "fast primality MR_2sprp_u128" }
#ifdef USE_AVX512F
    ,
    { "fermat_u52",   t_fermat_prp_52x8,   "fast primality fermat_u52" },
    { "MR_2sprp_u52", t_MR_2sprp_52x8, "fast primality MR_2sprp_u52" },
    { "fermat_u104",   t_fermat_prp_104x8,   "fast primality fermat_u104" },
    { "MR_2sprp_u104", t_MR_2sprp_104x8, "fast primality MR_2sprp_u104" }
#endif
};

const tk_module tk_module_primality = {
    "primality",
    "APRCL/BPSW/MR and tinyprp kernels",
    tk__primality_tests,
    (int)(sizeof tk__primality_tests / sizeof tk__primality_tests[0])
};
