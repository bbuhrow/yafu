/*
tinyrho: Pollard rho with Brent's cycle finding for small inputs, in the
style of microecm / tinyecm:
  - tinyrho64:   one n < 2^64, scalar (Montgomery, R = 2^64);
  - tinyrho128:  one n < 2^128, scalar (two-word Montgomery, R = 2^128);
  - tinyrho52_x8_list:  a list of n < 2^52, AVX-512 (IFMA or AVX-512F),
                        one input per lane, lanes refilled as inputs finish;
  - tinyrho104_x8_list: a list of n < 2^104, likewise, two 52-bit digits.
Without USE_AVX512F the list routines run the scalar code per input.
Iteration f(x) = x^2 + c on Montgomery residues (the constant c is used as
a Montgomery residue directly, so no conversions are needed); a product of
(x - y) is accumulated and a gcd taken every few hundred steps; a gcd equal
to n is backtracked (scalar) or retried with the next c (vector).
maxit = the iteration limit per input (f evaluations, over all retries);
0 picks a default, 32 * 2^(bits(n) / 4) + 4096: about 20x the mean for a
balanced semiprime.
This file is dual-licensed like microecm.h (FreeBSD or MPL 2.0).
*/
#ifndef TINYRHO_H_INCLUDED
#define TINYRHO_H_INCLUDED
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// A nontrivial factor of n, or 1 if none was found (n even: 2).
uint64_t tinyrho64(uint64_t n, uint32_t c, uint64_t maxit);

// n = n[0] + 2^64 n[1].  Returns 1 and a nontrivial factor in f[0], f[1],
// else 0 and f = 1.
int tinyrho128(const uint64_t n[2], uint64_t f[2], uint32_t c, uint64_t maxit);

// f64[i] = a factor of q64[i] (< 2^52), or 1.
void tinyrho52_x8_list(const uint64_t *q64, uint64_t *f64, uint32_t num, uint64_t maxit);

// q64[2i], q64[2i+1] = low / high 52 bits of input i (< 2^104, as
// getfactor_tecm_x8_list); f64[2i], f64[2i+1] = low / high 64 bits of a
// factor, or 1, 0.
void tinyrho104_x8_list(const uint64_t *q64, uint64_t *f64, uint32_t num, uint64_t maxit);

// iterations (f evaluations, all lanes) run by the list routines since the
// last call, for benchmarks (not thread-safe)
uint64_t tinyrho_list_iterations(void);

#ifdef __cplusplus
}
#endif

#endif
