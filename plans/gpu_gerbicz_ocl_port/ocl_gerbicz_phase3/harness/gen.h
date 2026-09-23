/* gen.h -- deterministic synthetic collcase generator (plan task 0.3). */
#ifndef COLLHARNESS_GEN_H
#define COLLHARNESS_GEN_H

#include <stdint.h>
#include "collcase.h"

typedef struct {
	uint32_t n;              /* total element count (best-effort; edge
	                             cases + skew are added on top and may
	                             push the actual count slightly over) */
	uint32_t root_bytes;      /* 4 or 8 */
	uint32_t key_bits;        /* signed key range is +/-2^(key_bits-1) */
	uint32_t shift;           /* bits allocated to p in each value */
	uint32_t bucket_hash;     /* 0 or 1 -- which compute_bucket mode */
	uint32_t hash_word_cap;   /* 0 = default (4096) */
	uint64_t seed;            /* splitmix64 seed -- fully determines
	                             every random choice below */
	double collision_density; /* fraction of filler slots that form a
	                             random duplicate pair instead of a
	                             singleton, in [0,1] */
	uint32_t bucket_skew_count; /* 0 disables; else this many extra
	                             distinct keys are forced into one
	                             bucket under the case's bucket_hash
	                             mode, to exercise bucket growth */
	int include_edge_cases;  /* 1 = prepend the fixed 0.3 edge suite */
	int with_expected;       /* 1 = compute and embed cpuref_exact +
	                             cpuref_stats output */
} gen_params_t;

void gen_params_defaults(gen_params_t *p);

/* Builds *out from scratch (out must not already hold live pointers --
 * call collcase_init() first, or collcase_free() an old one). Returns
 * 0 on success. */
int gen_generate(const gen_params_t *p, collcase_t *out);

#endif /* COLLHARNESS_GEN_H */
