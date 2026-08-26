/*--------------------------------------------------------------------
Shared bucket-hash definitions for the Gerbicz collision engine.

This header is the single source of truth for the bucket count, mask,
hash mix constant, and bucket-index function used by both:

  - cub/collision_engine.cu (host-side engine + scatter_roots_kernel +
    filter_per_bucket_kernel)
  - gnfs/poly/stage1/stage1_core_gpu/stage1_core.cu (sieve_kernel_trans_fused_*
    when the collfuse path is enabled)

Both translation units must agree on these values exactly. If they
drift, the fused trans kernel will write keys into different buckets
than the engine's filter expects to read them from, and parity will
fail in confusing ways.
--------------------------------------------------------------------*/

#ifndef _COLLISION_BUCKET_H_
#define _COLLISION_BUCKET_H_

#include <stdint.h>

#define LOG2_NUM_BUCKETS  14u
#define NUM_BUCKETS       (1u << LOG2_NUM_BUCKETS)
#define BUCKET_MASK       (NUM_BUCKETS - 1u)
#define BUCKET_HASH_MIX   0x9E3779B97F4A7C15ULL

#ifdef __CUDACC__
__host__ __device__
#endif
static inline uint32_t
compute_bucket(uint64_t key, int hash_mode) {
	if (hash_mode) {
		return (uint32_t)((key * BUCKET_HASH_MIX) >>
				(64 - LOG2_NUM_BUCKETS));
	}
	return (uint32_t)(key & BUCKET_MASK);
}

#endif /* !_COLLISION_BUCKET_H_ */
