/* prng.h -- deterministic, portable PRNG.
 *
 * Never use libc rand()/srand() for test-case generation: its
 * sequence is implementation-defined, so the same seed produces
 * different data under glibc, MSVC's CRT, etc, which would break
 * reproducibility across the gcc/clang/MSVC builds this harness
 * must support. splitmix64 is a tiny, well-known, fully specified
 * generator -- same seed, same bits, everywhere.
 */
#ifndef COLLHARNESS_PRNG_H
#define COLLHARNESS_PRNG_H

#include <stdint.h>

typedef struct {
	uint64_t state;
} cc_rng_t;

static inline void
cc_rng_seed(cc_rng_t *r, uint64_t seed)
{
	r->state = seed;
}

static inline uint64_t
cc_rng_next(cc_rng_t *r)
{
	uint64_t z = (r->state += 0x9E3779B97F4A7C15ULL);
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
	return z ^ (z >> 31);
}

/* Uniform in [0, bound). bound must be > 0. Rejection-free (slight
 * modulo bias for non-power-of-two bounds, fine for test data). */
static inline uint64_t
cc_rng_below(cc_rng_t *r, uint64_t bound)
{
	if (bound == 0)
		return 0;
	return cc_rng_next(r) % bound;
}

/* Uniform in [lo, hi] inclusive. static inline (rather than plain
 * static) so unused-in-this-TU is not a warning: not every file that
 * includes this header needs every helper. */
static inline uint64_t
cc_rng_range(cc_rng_t *r, uint64_t lo, uint64_t hi)
{
	if (hi <= lo)
		return lo;
	return lo + cc_rng_below(r, hi - lo + 1);
}

#endif /* COLLHARNESS_PRNG_H */
