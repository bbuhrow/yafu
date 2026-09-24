/* bucket_hash.h -- exact mirror of the project's collision_bucket.h.
 * Both translation units (the real engine and this harness) must
 * agree on these constants bit-for-bit, or the CPU stats model will
 * not reproduce the GPU engine's bucket assignment. */
#ifndef COLLHARNESS_BUCKET_HASH_H
#define COLLHARNESS_BUCKET_HASH_H

#include <stdint.h>

#define LOG2_NUM_BUCKETS  14u
#define NUM_BUCKETS       (1u << LOG2_NUM_BUCKETS)
#define BUCKET_MASK       (NUM_BUCKETS - 1u)
#define BUCKET_HASH_MIX   0x9E3779B97F4A7C15ULL

static inline uint32_t
compute_bucket(uint64_t key, int hash_mode)
{
	if (hash_mode) {
		return (uint32_t)((key * BUCKET_HASH_MIX) >>
				(64 - LOG2_NUM_BUCKETS));
	}
	return (uint32_t)(key & BUCKET_MASK);
}

#endif /* COLLHARNESS_BUCKET_HASH_H */
