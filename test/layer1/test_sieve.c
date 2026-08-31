/*----------------------------------------------------------------------
 Layer 1 -- ysieve (sieve of Eratosthenes) via the wrapper.c public API.

 Two brief checks, each run single- and multi-threaded:
   * prime counts in ~10^6-wide windows near powers of 2 (up to 2^40) vs
     known scalars (sympy-verified here; publicly checkable);
   * the actual primes in a few small (10^4) windows at offsets below 10^14
     vs GMP's mpz_nextprime enumeration.

 Every window uses EVEN endpoints, so the count/enumeration is identical
 whether soe_wrapper treats the range as inclusive or half-open (there is no
 prime at an even bound > 2) -- that keeps the checks independent of that
 convention. Each check runs under THREADS=1 and THREADS=2 and asserts the
 two agree, covering the threaded path.

 API (ysieve/soe.h, defined in ysieve/wrapper.c):
   soe_staticdata_t *soe_init(int vflag, int threads, int blocksize);
   uint64_t *soe_wrapper(sdata, lowlimit, highlimit, count, *num_p,
                         PRIMES_TO_FILE, PRIMES_TO_SCREEN);
       count != 0 -> count only, returns NULL, sets *num_p to the count;
       count == 0 -> returns a malloc'd array of the primes, *num_p = length.
   void soe_finalize(soe_staticdata_t *sdata);
 Public domain.
----------------------------------------------------------------------*/
#include "testkit.h"
#include <gmp.h>
#include <stdlib.h>
#include "soe.h"

#define SOE_VFLAG      0        /* quiet */
#define SOE_BLOCKSIZE  32768    /* 32 KB (>1024 -> used as bytes) */

/* portable mpz <- uint64 (mpz_set_ui takes only the low 32 bits on LLP64) */
static void mpz_set_u64(mpz_t r, uint64_t v)
{
    mpz_set_ui(r, (unsigned long)(v >> 32));
    mpz_mul_2exp(r, r, 32);
    mpz_add_ui(r, r, (unsigned long)(v & 0xFFFFFFFFUL));
}

/* ---- prime-count windows near powers of 2 (endpoints even) ---- */
typedef struct { uint64_t lo, hi, count; } sieve_window;
static const sieve_window sieve_windows[] = {
    { 1048576ULL,        2048576ULL,        70265ULL },   /* [2^20,      +10^6) */
    { 1073741824ULL,     1074741824ULL,     48129ULL },   /* [2^30,      +10^6) */
    { 68719476736ULL,    68720476736ULL,    40178ULL },   /* [2^36,      +10^6) */
    { 1099510627776ULL,  1099511627776ULL,  36000ULL }    /* [2^40-10^6, 2^40)  */
};
#define N_WINDOWS (int)(sizeof sieve_windows / sizeof sieve_windows[0])

/* ---- small value windows (base even, width 10^4), offsets < 10^14 ---- */
typedef struct { uint64_t base, width; } sieve_vrange;
static const sieve_vrange sieve_vranges[] = {
    { 1000000000000ULL,  10000ULL },   /* 10^12 */
    { 10000000000000ULL, 10000ULL },   /* 10^13 */
    { 90000000000000ULL, 10000ULL }    /* 9*10^13 */
};
#define N_VRANGES (int)(sizeof sieve_vranges / sizeof sieve_vranges[0])

/* ================================================================== *
 * Known-scalar counts, single- and multi-threaded.
 * ================================================================== */
static void t_prime_counts(tk_ctx *tk)
{
    soe_staticdata_t *s1 = soe_init(SOE_VFLAG, 1, SOE_BLOCKSIZE);
    soe_staticdata_t *s2 = soe_init(SOE_VFLAG, 2, SOE_BLOCKSIZE);
    int i;
    TK_REQUIRE(tk, s1 && s2, "soe_init failed");

    for (i = 0; i < N_WINDOWS; i++) {
        const sieve_window *w = &sieve_windows[i];
        uint64_t n1 = 0, n2 = 0;
        soe_wrapper(s1, w->lo, w->hi, 1, &n1, 0, 0);   /* count only, 1 thread  */
        soe_wrapper(s2, w->lo, w->hi, 1, &n2, 0, 0);   /* count only, 2 threads */
        TK_CHECKF(tk, n1 == w->count,
                  "pi[%llu,%llu) = %llu, expected %llu (1 thread)",
                  (unsigned long long)w->lo, (unsigned long long)w->hi,
                  (unsigned long long)n1, (unsigned long long)w->count);
        TK_CHECKF(tk, n2 == w->count,
                  "pi[%llu,%llu) = %llu, expected %llu (2 threads)",
                  (unsigned long long)w->lo, (unsigned long long)w->hi,
                  (unsigned long long)n2, (unsigned long long)w->count);
    }

    soe_finalize(s1);
    soe_finalize(s2);
}

/* ================================================================== *
 * Actual prime values vs GMP mpz_nextprime, single- and multi-threaded.
 * Returns the number of primes the sieve reported (for a 1-vs-2 cross-check).
 * ================================================================== */
static uint64_t check_values(tk_ctx *tk, soe_staticdata_t *s, int threads,
                             const sieve_vrange *v)
{
    uint64_t lo = v->base, hi = v->base + v->width;
    uint64_t np = 0, i;
    uint64_t *primes = soe_wrapper(s, lo, hi, 0, &np, 0, 0);   /* return the primes */
    mpz_t p, sp, HI;
    if (!TK_CHECKF(tk, primes != NULL, "soe_wrapper returned no array (%d thr)", threads))
        return 0;

    mpz_inits(p, sp, HI, NULL);
    mpz_set_u64(HI, hi);
    mpz_set_u64(p, lo);          /* lo is even (composite): nextprime walks up from here */

    for (i = 0; i < np; i++) {
        mpz_nextprime(p, p);                 /* GMP's i-th prime above lo */
        mpz_set_u64(sp, primes[i]);          /* sieve's i-th prime        */
        if (!TK_CHECKF(tk, mpz_cmp(sp, p) == 0,
                       "value[%llu] sieve/gmp mismatch (%d thr, base=%llu)",
                       (unsigned long long)i, threads, (unsigned long long)lo))
            break;
    }
    /* nothing skipped: the next prime past the last reported one is >= hi */
    if (i == np) {
        mpz_nextprime(p, p);
        TK_CHECKF(tk, mpz_cmp(p, HI) >= 0,
                  "sieve missed a prime below hi (%d thr, base=%llu)",
                  threads, (unsigned long long)lo);
    }

    mpz_clears(p, sp, HI, NULL);
    free(primes);
    return np;
}

static void t_prime_values(tk_ctx *tk)
{
    soe_staticdata_t *s1 = soe_init(SOE_VFLAG, 1, SOE_BLOCKSIZE);
    soe_staticdata_t *s2 = soe_init(SOE_VFLAG, 2, SOE_BLOCKSIZE);
    int i;
    TK_REQUIRE(tk, s1 && s2, "soe_init failed");

    for (i = 0; i < N_VRANGES; i++) {
        uint64_t n1 = check_values(tk, s1, 1, &sieve_vranges[i]);
        uint64_t n2 = check_values(tk, s2, 2, &sieve_vranges[i]);
        TK_CHECKF(tk, n1 == n2,
                  "thread-count disagreement: 1thr=%llu 2thr=%llu (base=%llu)",
                  (unsigned long long)n1, (unsigned long long)n2,
                  (unsigned long long)sieve_vranges[i].base);
    }

    soe_finalize(s1);
    soe_finalize(s2);
}

static const tk_test tk__sieve_tests[] = {
    { "prime_counts", t_prime_counts, "slow sieve" },
    { "prime_values", t_prime_values, "slow sieve gmp" }
};

const tk_module tk_module_sieve = {
    "sieve",
    "ysieve prime counts vs known scalars and primes vs GMP nextprime",
    tk__sieve_tests,
    (int)(sizeof tk__sieve_tests / sizeof tk__sieve_tests[0])
};
