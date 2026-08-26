/* DSO interface for Gerbicz-style GPU collision search in NFS stage 1. */

#ifndef _COLLISION_ENGINE_H_
#define _COLLISION_ENGINE_H_

#include <stdlib.h>
#include <cuda.h>

#ifdef __cplusplus
extern "C"
{
#endif

#define COLLISION_ENGINE_MIN_KEY_BITS 20

typedef struct {
	CUdeviceptr keys_in;
	CUdeviceptr data_in;
	CUdeviceptr q_batch;
	CUdeviceptr found_array;
	CUstream stream;

	size_t num_elements;
	int key_bits;
	int root_bytes;
	unsigned int shift;
	int bucket_hash;
	int debug;
	int collect_stats;

	unsigned int bucket_max;
	unsigned int candidate_count;
	unsigned int dedup_count;
	unsigned int value_match_count;
	unsigned int bucket_grow_count;
	unsigned int hash_cap_count;
	unsigned int match_arena_attempt_count;
	unsigned int match_arena_fallback_count;
	unsigned int match_arena_capacity_skip_count;
	/* filter_per_bucket_kernel stop-iteration histogram + per-category
	   bucket-size histograms. Populated only when collect_stats != 0;
	   otherwise the engine passes NULL and zero atomic ops happen on
	   the GPU side, so non-stats runs pay no overhead.

	   Stop-iter segments (21 iters each):
	     [0..20]   = total stops at each iter
	     [21..41]  = stops because cnt==0 (already filtered out)
	     [42..62]  = stops because hit MAX_FILTER_ITERS-1 (cap)
	   Converged stops = total - zero - cap (host-side derived).

	   Bucket-size histograms by stop-iter category (13 log2 bins each:
	   1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096+):
	     [63..75]  = initial cnt for fast stops (it <= 3)
	     [76..88]  = initial cnt for medium stops (it == 4)
	     [89..101] = initial cnt for slow stops (it >= 5)
	   Tells whether the slow tail buckets are the densest ones (which
	   would make "split oversized buckets" the natural filter lever). */
	unsigned int filter_iters_hist[102];
	float elapsed_ms;
} collision_data_t;

typedef void * (*collision_engine_init_func)(void);

typedef void (*collision_engine_free_func)(void * engine);

typedef void (*collision_engine_run_func)(void * engine,
				collision_data_t * collision_data);

#ifdef __cplusplus
}
#endif

#endif /* !_COLLISION_ENGINE_H_ */
