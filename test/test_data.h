/*----------------------------------------------------------------------
 YAFU modular test system -- shared test data & self-contained oracles.

 This layer is GMP-free and depends only on <stdint.h>, so it compiles and
 runs everywhere the framework does. It provides:

   * A deterministic 64-bit primality oracle (tk_is_prime_u64), correct for
     all n < 2^64 via a fixed Miller-Rabin base set. Used as the reference
     for the uint64 primality paths and by the generators.
   * Known-answer corpora: primes / composites (as decimal strings for the
     GMP-backed Layer-1 tests), strong pseudoprimes base 2, Carmichael
     numbers, and 64-bit semiprimes {N,p,q}.
   * Generators for random n-bit primes and semiprimes.

 Every entry is verified at authoring time; test_selfcheck re-verifies the
 whole corpus at runtime with the oracle, so a bad edit fails loudly.

 Public domain.
----------------------------------------------------------------------*/
#ifndef YAFU_TEST_DATA_H
#define YAFU_TEST_DATA_H

#include <stdint.h>
#include "testkit.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Deterministic, exact for all n < 2^64. Returns 1 if prime, else 0. */
int tk_is_prime_u64(uint64_t n);

/* Random prime with exactly `bits` bits (1 < bits <= 62 recommended so a
 * product of two still fits in 64 bits). Uses the ctx RNG (reproducible). */
uint64_t tk_gen_prime_u64(tk_ctx *tk, int bits);

/* Random semiprime N = p*q with p<=q, each ~bits/2 bits, N ~bits bits
 * (bits <= 64). p and q are returned; N is the product. p != q guaranteed. */
uint64_t tk_gen_semiprime_u64(tk_ctx *tk, int bits, uint64_t *p, uint64_t *q);

/* ---- decimal-string corpora (for GMP-backed tests) ---- */
extern const char *const tk_primes_dec[];      /* NULL-terminated */
extern const char *const tk_composites_dec[];  /* NULL-terminated */

/* ---- uint64 corpora ---- */
typedef struct { uint64_t n, p, q; } tk_semiprime64;

extern const uint64_t        tk_spsp2_u64[];       /* strong pseudoprimes, base 2 */
extern const int             tk_spsp2_u64_count;
extern const uint64_t        tk_carmichael_u64[];  /* Carmichael numbers          */
extern const int             tk_carmichael_u64_count;
extern const tk_semiprime64  tk_semiprimes_u64[];  /* {N,p,q}, p<=q               */
extern const int             tk_semiprimes_u64_count;

/* The module descriptor for the corpus self-check (Layer -1). */
extern const tk_module tk_module_selfcheck;

#ifdef __cplusplus
}
#endif
#endif /* YAFU_TEST_DATA_H */
