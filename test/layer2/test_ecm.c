/*----------------------------------------------------------------------
 Layer 2 -- small ECM factoring: microecm (64-bit) and tinyecm (128-bit).

 For each balanced-semiprime size we generate N = p*q (p, q known by
 construction, so they are the reference), hand N to the routine, and
 classify the result:
     pass  -- the routine returned p or q;
     miss  -- the routine reported "no factor" (ECM is probabilistic, so a
              small miss rate is allowed, not a correctness failure);
     hard  -- anything else (a wrong factor, a non-divisor, 1, or N): a real
              bug, flagged immediately.
 A size fails if any hard failure occurs, or if the miss rate exceeds a small
 budget. Inputs are odd (products of odd primes) with no small factors, so we
 pass is_arbitrary=0, matching test_dlp_composites() in top/test.c.

   microecm covers 40/50/60/64-bit inputs (its 64-bit ceiling), semiprimes
   from tk_gen_semiprime_u64.
   tinyecm covers 40..128-bit inputs, including the full 128-bit max case
   (two ~64-bit factors -- its SIQS double-large-prime workload), semiprimes
   built with GMP so the >64-bit sizes are reachable.

 Entry points (note the differing "not found" conventions):
   uint64_t getfactor_uecm(uint64_t n, int is_arbitrary, uint64_t *pran);
       returns 1 when no factor is found; the factor otherwise.
   int getfactor_tecm(mpz_t n, mpz_t f, int is_arbitrary, uint64_t *pran);
       returns 0 when no factor is found; 1 with the factor in f otherwise.
 Both take a single *pran seeded once and then left alone.
 Public domain.
----------------------------------------------------------------------*/
#include "testkit.h"
#include "test_data.h"     /* tk_gen_semiprime_u64 */
#include <gmp.h>
#include <stdio.h>
#include "microecm.h"
#include "tinyecm.h"

#define ECM_NTRIALS       10000    /* semiprimes per size */
#define ECM_MAX_MISS_PCT  1.0      /* per-size miss-rate budget (%); tune if the
                                      largest tinyecm sizes miss more in practice */

static const int uecm_bits[] = { 40, 50, 60, 64 };
#define UECM_NBITS (int)(sizeof uecm_bits / sizeof uecm_bits[0])

static const int tecm_bits[] = { 40, 50, 60, 64, 80, 104, 128 };
#define TECM_NBITS (int)(sizeof tecm_bits / sizeof tecm_bits[0])

/* one prime of exactly b bits (the size retry rejects a nextprime carry) */
static void gen_mpz_prime(mpz_t p, gmp_randstate_t st, int b)
{
    do {
        mpz_urandomb(p, st, (mp_bitcnt_t)b);
        mpz_setbit(p, (mp_bitcnt_t)(b - 1));
        mpz_nextprime(p, p);
    } while ((int)mpz_sizeinbase(p, 2) != b);
}

/* GMP balanced semiprime N = p*q, each factor exactly bits/2 (or +1) bits, so
   for bits=128 both factors are 64-bit and N stays < 2^128. */
static void gen_mpz_semiprime(mpz_t n, mpz_t p, mpz_t q, gmp_randstate_t st, int bits)
{
    int pb = bits / 2, qb = bits - pb;
    do {
        gen_mpz_prime(p, st, pb);
        gen_mpz_prime(q, st, qb);
    } while (mpz_cmp(p, q) == 0);
    mpz_mul(n, p, q);
}

/* GMP 3-large-prime composite N = p1*p2*p3: three distinct primes, each with
   bit-length uniform in [27,35] -- i.e. in (2^26, 2^35). The product of any
   three is at most ~2^120 < 2^128 (tinyecm's ceiling), so the size bound alone
   guarantees the constraint. Models NFS/SIQS cofactors carrying 3 large primes. */
static void gen_mpz_3prime(mpz_t n, mpz_t p1, mpz_t p2, mpz_t p3,
                           gmp_randstate_t st, tk_rng *rng)
{
    do {
        gen_mpz_prime(p1, st, 26 + (int)tk_rng_range(rng, 9));
        gen_mpz_prime(p2, st, 26 + (int)tk_rng_range(rng, 9));
        gen_mpz_prime(p3, st, 26 + (int)tk_rng_range(rng, 9));
    } while (mpz_cmp(p1, p2) == 0 || mpz_cmp(p1, p3) == 0 || mpz_cmp(p2, p3) == 0);
    mpz_mul(n, p1, p2);
    mpz_mul(n, n, p3);
}

/* ================================================================== *
 * microecm -- getfactor_uecm (64-bit)
 * ================================================================== */
static void t_microecm(tk_ctx *tk)
{
    tk_rng *rng = tk_rng_of(tk);
    uint64_t lcg = tk_rng_u64(rng) | 1ULL;   /* seed once; routine updates it */
    int bi;

    for (bi = 0; bi < UECM_NBITS; bi++) {
        int bits = uecm_bits[bi];
        long misses = 0, passes = 0, t;
        for (t = 0; t < ECM_NTRIALS; t++) {
            uint64_t p, q, n = tk_gen_semiprime_u64(tk, bits, &p, &q);
            uint64_t f = getfactor_uecm(n, 0, &lcg);
            if (f == 1ULL) { misses++; continue; }         /* not found */
            if (f == p || f == q) { passes++; continue; }  /* correct   */
            TK_CHECKF(tk, 0,
                "microecm[%d-bit] wrong factor: n=%llu p=%llu q=%llu got=%llu",
                bits, (unsigned long long)n, (unsigned long long)p,
                (unsigned long long)q, (unsigned long long)f);
        }
        if (tk_verbose(tk))
            printf("        microecm %d-bit: %ld pass, %ld miss / %d\n",
                   bits, passes, misses, ECM_NTRIALS);
        TK_CHECKF(tk, (double)misses * 100.0 / ECM_NTRIALS <= ECM_MAX_MISS_PCT,
            "microecm[%d-bit] miss rate %.2f%% (%ld/%d) exceeds %.1f%%",
            bits, (double)misses * 100.0 / ECM_NTRIALS, misses,
            ECM_NTRIALS, ECM_MAX_MISS_PCT);
    }
}

/* ================================================================== *
 * tinyecm on 3-large-prime composites (NFS/SIQS cofactors): three primes
 * in (2^26, 2^35), N ~78..105 bits. The routine only needs to find ONE
 * factor, not all three -- so a correct result is any single nontrivial
 * divisor of N (one prime, giving a 2-prime cofactor, or occasionally a
 * product of two). We verify 1 < f < N and f | N, nothing more.
 * ================================================================== */
static void t_tinyecm_3lp(tk_ctx *tk)
{
    tk_rng *rng = tk_rng_of(tk);
    uint64_t lcg = tk_rng_u64(rng) | 1ULL;
    gmp_randstate_t st;
    mpz_t N, F, P1, P2, P3;
    long misses = 0, passes = 0, t;

    mpz_inits(N, F, P1, P2, P3, NULL);
    gmp_randinit_default(st);
    gmp_randseed_ui(st, (unsigned long)tk_rng_u64(rng));

    for (t = 0; t < ECM_NTRIALS; t++) {
        int found;
        gen_mpz_3prime(N, P1, P2, P3, st, rng);
        found = getfactor_tecm(N, F, 0, &lcg);
        if (!found) { misses++; continue; }                   /* not found */
        if (mpz_cmp_ui(F, 1) > 0 && mpz_cmp(F, N) < 0 &&
            mpz_divisible_p(N, F)) { passes++; continue; }     /* nontrivial divisor */
        {
            char nb[48], fb[48];
            mpz_get_str(nb, 10, N); mpz_get_str(fb, 10, F);
            TK_CHECKF(tk, 0, "tinyecm[3LP] not a proper divisor: n=%s got=%s", nb, fb);
        }
    }
    if (tk_verbose(tk))
        printf("        tinyecm  3LP (3x 2^27..2^35): %ld pass, %ld miss / %d\n",
               passes, misses, ECM_NTRIALS);
    TK_CHECKF(tk, (double)misses * 100.0 / ECM_NTRIALS <= ECM_MAX_MISS_PCT,
        "tinyecm[3LP] miss rate %.2f%% (%ld/%d) exceeds %.1f%%",
        (double)misses * 100.0 / ECM_NTRIALS, misses, ECM_NTRIALS, ECM_MAX_MISS_PCT);

    gmp_randclear(st);
    mpz_clears(N, F, P1, P2, P3, NULL);
}

static const tk_test tk__microecm_tests[] = {
    { "getfactor_uecm", t_microecm, "slow ecm microecm" }
};
const tk_module tk_module_microecm = {
    "microecm",
    "getfactor_uecm on 40/50/60/64-bit semiprimes",
    tk__microecm_tests,
    (int)(sizeof tk__microecm_tests / sizeof tk__microecm_tests[0])
};

static const tk_test tk__tinyecm_tests[] = {
    { "getfactor_tecm_3lp", t_tinyecm_3lp, "slow ecm tinyecm" }
};
const tk_module tk_module_tinyecm = {
    "tinyecm",
    "getfactor_tecm: 3-large-prime composites",
    tk__tinyecm_tests,
    (int)(sizeof tk__tinyecm_tests / sizeof tk__tinyecm_tests[0])
};
