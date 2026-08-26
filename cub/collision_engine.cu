#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <cuda.h>
#include <cub/cub.cuh>

#include "collision_engine.h"
#include "../gnfs/poly/stage1/stage1_core_gpu/stage1_core.h"

typedef unsigned int uint32;
typedef unsigned long long uint64;
typedef long long int64;

#if defined(_WIN32) || defined (_WIN64)
	#define COLLISION_ENGINE_DECL __declspec(dllexport)
#else
	#define COLLISION_ENGINE_DECL __attribute__((visibility("default")))
#endif

#include "collision_bucket.h"

/* LOG2_NUM_BUCKETS, NUM_BUCKETS, BUCKET_MASK, BUCKET_HASH_MIX and
   compute_bucket() live in collision_bucket.h so any future kernel
   that wants to bucket-scatter (e.g., a fused trans + bucket-scatter
   attempt — see FUSED_TRANS_SCATTER_PLAN.md) can include the same
   constants. The fused-trans spike itself was retired 2026-05-25
   (register pressure pushed trans 5 → 4 blocks/SM, net regression). */

#define MAX_FILTER_ITERS  20
#define BLOCK_THREADS     128
#define CANDIDATE_CAP     (1u << 22)
#define VALUE_MATCH_CAP   CANDIDATE_CAP
#define MATCH_ARENA_WIDTH 8u
#define MAX_C_ILOG2       20u
#define MAX_DSIZE         ((1u << MAX_C_ILOG2) + 64u)
#define MAX_SSIZE         (((1u << MAX_C_ILOG2) >> 5) + 64u)

#define CUDA_TRY(func) \
	{ 			 					\
		cudaError_t status = func;				\
		if (status != cudaSuccess) {				\
			const char * str = cudaGetErrorString(status);	\
			if (!str)					\
				str = "Unknown";			\
			printf("error (%s:%d): %s\n", __FILE__, __LINE__, str);\
			exit(-1);					\
		}							\
	}

#define CUDA_CHECK_LAST() CUDA_TRY(cudaGetLastError())

__device__ static inline uint32
compute_ilog2(uint32 cnt, uint32 key_bits) {
	if (cnt == 0)
		return 5;
	uint32 lg = 31u - __clz(cnt);
	uint32 il = max(5u, lg + 5u);
	uint32 hash_bits_available = key_bits - LOG2_NUM_BUCKETS;
	if (il > hash_bits_available)
		il = hash_bits_available;
	return il;
}

__device__ static inline uint32
compute_capped_ilog2(uint32 cnt, uint32 key_bits, uint32 max_tsize_words) {
	uint32 il = compute_ilog2(cnt, key_bits);
	uint32 max_il = 5u;

	if (max_tsize_words > 1u)
		max_il += 31u - __clz(max_tsize_words);

	if (il > max_il)
		il = max_il;
	return il;
}

__device__ static void
clear_table_parallel(uint32 *tbl, uint32 nwords) {
	for (uint32 j = threadIdx.x; j < nwords; j += BLOCK_THREADS)
		tbl[j] = 0;
}

__device__ static uint32
gcd32_generic(uint32 a, uint32 b) {
	while (b != 0) {
		uint32 t = a % b;
		a = b;
		b = t;
	}
	return a;
}

__device__ static void
store_hit_collision(found_t *found_array, uint32 found_array_size,
		uint32 p1, uint32 p2, int64 root, specialq_t *q) {

	uint32 index = atomicAdd(&found_array[0].p1, 1);

	if (index < found_array_size - 1) {
		found_t *f = found_array + index + 1;

		f->p1 = p1;
		f->p2 = p2;
		f->q = q->p;
		f->qroot = q->root;
		f->offset = root;
	}
}

__global__ void
scatter_roots_kernel(const void * __restrict__ roots,
		uint32 root_bytes, uint32 N, uint32 max_per_bucket,
		int hash_mode, uint32 * __restrict__ bucket_count,
		uint32 * __restrict__ overflow_flag,
		uint64 * __restrict__ bucket_storage) {

	uint32 tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= N)
		return;

	uint64 k;
	if (root_bytes == sizeof(uint32))
		k = ((const uint32 *)roots)[tid];
	else
		k = ((const uint64 *)roots)[tid];
	if (k == 0)
		return;

	uint32 bucket = compute_bucket(k, hash_mode);
#ifdef SCATTER_DIRECT_ATOMIC
	uint32 slot = atomicAdd(&bucket_count[bucket], 1u);
#else
	uint32 lane = threadIdx.x & 31u;
	uint32 mask = __match_any_sync(0xFFFFFFFFu, bucket);
	uint32 leader = __ffs(mask) - 1u;
	uint32 count = __popc(mask);

	uint32 base = 0;
	if (lane == leader)
		base = atomicAdd(&bucket_count[bucket], count);
	base = __shfl_sync(0xFFFFFFFFu, base, leader);

	uint32 within = __popc(mask & ((1u << lane) - 1u));
	uint32 slot = base + within;
#endif
	if (slot >= max_per_bucket) {
		atomicOr(overflow_flag, 1u);
		return;
	}

	bucket_storage[(size_t)bucket * max_per_bucket + slot] = k;
}

extern __shared__ uint32 s_hash[];

__global__ void
filter_per_bucket_kernel(const uint32 * __restrict__ bucket_count_in,
		uint32 max_per_bucket,
		uint64 * __restrict__ arr_a,
		uint64 * __restrict__ arr_b,
		uint64 * __restrict__ candidate_keys,
		uint32 * __restrict__ candidate_cnt,
		uint32 * __restrict__ candidate_overflow,
		uint32 candidate_cap,
		uint32 key_bits,
		uint32 max_tsize_words,
		uint32 * iters_hist_out) {

	uint32 bucket = blockIdx.x;
	uint32 cnt = bucket_count_in[bucket];
	if (cnt > max_per_bucket)
		cnt = max_per_bucket;
	if (cnt == 0)
		return;

	/* Save initial bucket size for tail diagnostics. */
	uint32 cnt0 = cnt;

	size_t offset = (size_t)bucket * max_per_bucket;
	uint64 *arr_in = arr_a + offset;
	uint64 *arr_out = arr_b + offset;

	uint32 *T = s_hash;
	uint32 *T2 = s_hash + max_tsize_words;
	uint32 *T3 = s_hash + 2u * max_tsize_words;

	uint32 ilog2 = compute_capped_ilog2(cnt, key_bits, max_tsize_words);
	uint32 tsize = 1u << (ilog2 - 5u);
	uint32 hash_value2 = (ilog2 >= 32u) ?
			0xFFFFFFFFu : ((1u << ilog2) - 1u);
	uint32 my_shift2 = LOG2_NUM_BUCKETS;

	clear_table_parallel(T, tsize);
	clear_table_parallel(T2, tsize);
	__syncthreads();

	for (uint32 j = threadIdx.x; j < cnt; j += BLOCK_THREADS) {
		uint64 item = arr_in[j];
		uint32 hv = (uint32)((item >> my_shift2) & hash_value2);
		uint32 bit = 1u << (hv & 31u);
		uint32 prev = atomicOr(&T2[hv >> 5], bit);
		if (prev & bit)
			atomicOr(&T[hv >> 5], bit);
	}
	__syncthreads();

	__shared__ uint32 s_nsize[MAX_FILTER_ITERS];
	__shared__ uint32 s_cnt_out;

	bool emit = false;
	for (int it = 0; it < MAX_FILTER_ITERS; it++) {
		uint32 my_shift = my_shift2;
		uint32 hash_value = hash_value2;

		ilog2 = compute_capped_ilog2(cnt, key_bits, max_tsize_words);
		tsize = 1u << (ilog2 - 5u);

		uint32 *U, *U2, *U3;
		if ((it & 1) == 0) { U = T; U2 = T2; U3 = T3; }
		else { U = T3; U2 = T2; U3 = T; }

		my_shift2 = my_shift + 6u;
		if (my_shift2 + ilog2 > key_bits)
			my_shift2 = LOG2_NUM_BUCKETS;
		hash_value2 = (ilog2 >= 32u) ?
				0xFFFFFFFFu : ((1u << ilog2) - 1u);

		clear_table_parallel(U2, tsize);
		clear_table_parallel(U3, tsize);
		if (threadIdx.x == 0)
			s_cnt_out = 0;
		__syncthreads();

		for (uint32 j_base = 0; j_base < cnt; j_base += BLOCK_THREADS) {
			uint32 j = j_base + threadIdx.x;
			bool survives = false;
			uint64 item = 0;

			if (j < cnt) {
				item = arr_in[j];
				uint32 hv = (uint32)((item >> my_shift) & hash_value);
				uint32 bit = 1u << (hv & 31u);
				if (U[hv >> 5] & bit) {
					uint32 hv2 = (uint32)((item >> my_shift2) & hash_value2);
					uint32 bit2 = 1u << (hv2 & 31u);
					uint32 prev = atomicOr(&U2[hv2 >> 5], bit2);
					if (prev & bit2)
						atomicOr(&U3[hv2 >> 5], bit2);
					survives = true;
				}
			}

			uint32 vote = __ballot_sync(0xFFFFFFFFu, survives);
			uint32 lane = threadIdx.x & 31u;
			uint32 leader_slot = 0;
			if (vote) {
				uint32 wcount = __popc(vote);
				if (lane == 0)
					leader_slot = atomicAdd(&s_cnt_out, wcount);
				leader_slot = __shfl_sync(0xFFFFFFFFu, leader_slot, 0);
				if (survives) {
					uint32 within = __popc(vote & ((1u << lane) - 1u));
					arr_out[leader_slot + within] = item;
				}
			}
		}

		__syncthreads();
		uint32 cnt_out = s_cnt_out;
		if (threadIdx.x == 0)
			s_nsize[it] = cnt_out;
		__syncthreads();

		cnt = cnt_out;
		uint64 *tmp = arr_in;
		arr_in = arr_out;
		arr_out = tmp;

		bool stop_zero = (cnt == 0);
		bool stop_cap = !stop_zero && (it == MAX_FILTER_ITERS - 1);
		bool stop_conv = !stop_zero && !stop_cap && (it >= 3 &&
				s_nsize[it - 3] == s_nsize[it]);
		bool stop = stop_zero || stop_cap || stop_conv;
		if (stop) {
			if (iters_hist_out && threadIdx.x == 0) {
				/* Layout (102 slots):
				   [0..20]   total stops per iter
				   [21..41]  stops because cnt==0
				   [42..62]  stops because cap hit
				   [63..75]  initial-cnt log2 hist, fast (it<=3)
				   [76..88]  initial-cnt log2 hist, medium (it==4)
				   [89..101] initial-cnt log2 hist, slow (it>=5) */
				uint32 idx = (uint32)it;
				if (idx > 20u)
					idx = 20u;
				atomicAdd(&iters_hist_out[idx], 1u);
				if (stop_zero)
					atomicAdd(&iters_hist_out[21u + idx], 1u);
				else if (stop_cap)
					atomicAdd(&iters_hist_out[42u + idx], 1u);

				/* Bucket-size log2 hist: 13 bins (1, 2, 4, ...,
				   4096+). __clz(0)==undefined so guard cnt0==0 even
				   though we returned early above. */
				uint32 size_bin = (cnt0 == 0u) ? 0u :
						(31u - (uint32)__clz(cnt0));
				if (size_bin > 12u)
					size_bin = 12u;
				uint32 cat_base = (it <= 3) ? 63u :
						(it == 4) ? 76u : 89u;
				atomicAdd(&iters_hist_out[cat_base + size_bin], 1u);
			}
			emit = true;
			break;
		}
	}

	if (!emit || cnt == 0)
		return;

	__shared__ uint32 s_emit_base;
	if (threadIdx.x == 0)
		s_emit_base = atomicAdd(candidate_cnt, cnt);
	__syncthreads();

	uint32 base = s_emit_base;
	if (base + cnt > candidate_cap) {
		if (threadIdx.x == 0)
			atomicOr(candidate_overflow, 1u);
		return;
	}

	for (uint32 j = threadIdx.x; j < cnt; j += BLOCK_THREADS)
		candidate_keys[base + j] = arr_in[j];
}

__global__ void
dedup_kernel(const uint64 * __restrict__ sorted_keys, uint32 csize,
		uint64 * __restrict__ dedup_keys,
		uint32 * __restrict__ dedup_cnt) {

	uint32 tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= csize)
		return;

	uint64 k = sorted_keys[tid];
	bool is_first = (tid == 0) || (sorted_keys[tid - 1] != k);
	bool has_neighbor = (tid + 1 < csize) && (sorted_keys[tid + 1] == k);
	if (is_first && has_neighbor) {
		uint32 pos = atomicAdd(dedup_cnt, 1u);
		dedup_keys[pos] = k;
	}
}

__global__ void
count_secondary_kernel(const uint64 * __restrict__ dedup_keys,
		uint32 csize, uint32 hash_value,
		uint32 * __restrict__ D, uint32 * __restrict__ S) {

	uint32 tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= csize)
		return;

	uint64 k = dedup_keys[tid];
	uint32 hv = (uint32)(k & hash_value);
	atomicAdd(&D[hv], 1u);
	atomicOr(&S[hv >> 5], 1u << (hv & 31u));
}

__global__ void
scatter_secondary_kernel(const uint64 * __restrict__ dedup_keys,
		uint32 csize, uint32 hash_value,
		uint32 * __restrict__ D_pos,
		uint64 * __restrict__ X) {

	uint32 tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= csize)
		return;

	uint64 k = dedup_keys[tid];
	uint32 hv = (uint32)(k & hash_value);
	uint32 slot = atomicAdd(&D_pos[hv], 1u);
	X[slot] = k;
}

__device__ static inline int
find_candidate_slot(uint64 k, uint32 hash_value,
		const uint32 * __restrict__ S,
		const uint32 * __restrict__ D,
		const uint64 * __restrict__ X,
		uint32 *slot_out) {

	uint32 hv = (uint32)(k & hash_value);
	if ((S[hv >> 5] & (1u << (hv & 31u))) == 0)
		return 0;

	uint32 lo = D[hv];
	uint32 hi = D[hv + 1u];
	for (uint32 j = lo; j < hi; j++) {
		if (X[j] == k) {
			*slot_out = j;
			return 1;
		}
	}
	return 0;
}

__global__ void
count_matched_values_kernel(const void * __restrict__ roots,
		const uint32 * __restrict__ values,
		uint32 root_bytes, uint32 N,
		uint32 hash_value,
		const uint32 * __restrict__ S,
		const uint32 * __restrict__ D,
		const uint64 * __restrict__ X,
		uint32 dedup_count,
		uint32 * __restrict__ value_counts,
		uint32 * __restrict__ match_total,
		uint32 * __restrict__ value_overflow,
		uint32 value_cap) {

	uint32 tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= N)
		return;

	uint64 k;
	if (root_bytes == sizeof(uint32))
		k = ((const uint32 *)roots)[tid];
	else
		k = ((const uint64 *)roots)[tid];
	if (k == 0)
		return;

	uint32 slot;
	if (find_candidate_slot(k, hash_value, S, D, X, &slot) &&
	    slot < dedup_count) {
		uint32 total_pos = atomicAdd(match_total, 1u);
		atomicAdd(&value_counts[slot], 1u);
		if (total_pos >= value_cap)
			atomicOr(value_overflow, 1u);
	}
	(void)values;
}

__global__ void
scatter_matched_values_kernel(const void * __restrict__ roots,
		const uint32 * __restrict__ values,
		uint32 root_bytes, uint32 N,
		uint32 hash_value,
		const uint32 * __restrict__ S,
		const uint32 * __restrict__ D,
		const uint64 * __restrict__ X,
		uint32 dedup_count,
		uint32 * __restrict__ value_cursor,
		uint32 * __restrict__ matched_values,
		uint32 value_cap) {

	uint32 tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= N)
		return;

	uint64 k;
	if (root_bytes == sizeof(uint32))
		k = ((const uint32 *)roots)[tid];
	else
		k = ((const uint64 *)roots)[tid];
	if (k == 0)
		return;

	uint32 slot;
	if (find_candidate_slot(k, hash_value, S, D, X, &slot) &&
	    slot < dedup_count) {
		uint32 pos = atomicAdd(&value_cursor[slot], 1u);
		if (pos < value_cap)
			matched_values[pos] = values[tid];
	}
}

__global__ void
count_and_store_matched_values_kernel(const void * __restrict__ roots,
		const uint32 * __restrict__ values,
		uint32 root_bytes, uint32 N,
		uint32 hash_value,
		const uint32 * __restrict__ S,
		const uint32 * __restrict__ D,
		const uint64 * __restrict__ X,
		uint32 dedup_count,
		uint32 * __restrict__ value_counts,
		uint32 * __restrict__ match_total,
		uint32 * __restrict__ value_overflow,
		uint32 * __restrict__ matched_values,
		uint32 value_cap) {

	uint32 tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= N)
		return;

	uint64 k;
	if (root_bytes == sizeof(uint32))
		k = ((const uint32 *)roots)[tid];
	else
		k = ((const uint64 *)roots)[tid];
	if (k == 0)
		return;

	uint32 slot;
	if (find_candidate_slot(k, hash_value, S, D, X, &slot) &&
	    slot < dedup_count) {
		uint32 total_pos = atomicAdd(match_total, 1u);
		uint32 pos = atomicAdd(&value_counts[slot], 1u);
		size_t arena_pos = (size_t)slot * MATCH_ARENA_WIDTH + pos;

		if (total_pos >= value_cap) {
			atomicOr(value_overflow, 1u);
		}
		else if (pos < MATCH_ARENA_WIDTH && arena_pos < value_cap) {
			matched_values[arena_pos] = values[tid];
		}
		else {
			atomicOr(value_overflow, 1u);
		}
	}
}

__global__ void
emit_found_kernel(const uint64 * __restrict__ X,
		const uint32 * __restrict__ matched_values,
		const uint32 * __restrict__ value_offsets,
		uint32 dedup_count, uint32 pshift, uint32 root_bytes,
		specialq_t * __restrict__ q_batch,
		found_t * __restrict__ found_array) {

	uint32 slot = blockIdx.x * blockDim.x + threadIdx.x;
	if (slot >= dedup_count)
		return;

	uint32 lo = value_offsets[slot];
	uint32 hi = value_offsets[slot + 1u];
	uint32 mask = (1u << pshift) - 1u;
	int64 root = (root_bytes == sizeof(uint32)) ?
			(int64)(int)((uint32)X[slot]) : (int64)X[slot];

	for (uint32 i = lo; i + 1u < hi; i++) {
		uint32 v1 = matched_values[i];
		uint32 q1 = v1 >> pshift;
		uint32 p1 = v1 & mask;

		for (uint32 j = i + 1u; j < hi; j++) {
			uint32 v2 = matched_values[j];
			uint32 q2 = v2 >> pshift;
			uint32 p2 = v2 & mask;

			if (q1 == q2 && gcd32_generic(p1, p2) == 1u) {
				store_hit_collision(found_array, FOUND_ARRAY_SIZE,
						p1, p2, root, q_batch + q1);
			}
		}
	}
}

__global__ void
emit_found_arena_kernel(const uint64 * __restrict__ X,
		const uint32 * __restrict__ matched_values,
		const uint32 * __restrict__ value_counts,
		uint32 dedup_count, uint32 pshift, uint32 root_bytes,
		specialq_t * __restrict__ q_batch,
		found_t * __restrict__ found_array) {

	uint32 slot = blockIdx.x * blockDim.x + threadIdx.x;
	if (slot >= dedup_count)
		return;

	uint32 cnt = value_counts[slot];
	if (cnt > MATCH_ARENA_WIDTH)
		cnt = MATCH_ARENA_WIDTH;
	uint32 base = slot * MATCH_ARENA_WIDTH;
	uint32 mask = (1u << pshift) - 1u;
	int64 root = (root_bytes == sizeof(uint32)) ?
			(int64)(int)((uint32)X[slot]) : (int64)X[slot];

	for (uint32 i = 0; i + 1u < cnt; i++) {
		uint32 v1 = matched_values[base + i];
		uint32 q1 = v1 >> pshift;
		uint32 p1 = v1 & mask;

		for (uint32 j = i + 1u; j < cnt; j++) {
			uint32 v2 = matched_values[base + j];
			uint32 q2 = v2 >> pshift;
			uint32 p2 = v2 & mask;

			if (q1 == q2 && gcd32_generic(p1, p2) == 1u) {
				store_hit_collision(found_array, FOUND_ARRAY_SIZE,
						p1, p2, root, q_batch + q1);
			}
		}
	}
}

struct collision_engine {
	collision_engine()
	  : max_n(0), max_key_bits(0), max_per_bucket(0),
	    d_bucket_count(0), d_bucket_overflow(0),
	    d_candidate_overflow(0), d_arr_a(0), d_arr_b(0),
	    d_candidate_keys(0), d_candidate_cnt(0),
	    d_sorted_keys(0), d_dedup_keys(0), d_dedup_cnt(0),
	    d_D(0), d_D_scan(0), d_D_pos(0), d_S(0), d_X(0),
	    d_max_bucket(0), d_value_counts(0), d_value_offsets(0),
	    d_value_cursor(0), d_matched_values(0),
	    d_value_match_total(0), d_value_overflow(0),
	    d_filter_iters_hist(0),
	    d_cub_max_buckets(0), max_buckets_bytes(0),
	    d_cub_sort_candidates(0), sort_candidates_bytes(0),
	    d_cub_scan_D(0), scan_D_bytes(0), start_event(0), end_event(0),
	    d_cub_scan_values(0), scan_values_bytes(0) {
		CUDA_TRY(cudaMallocHost((void **)&h_max_bucket, sizeof(uint32)))
		CUDA_TRY(cudaMallocHost((void **)&h_bucket_overflow, sizeof(uint32)))
		CUDA_TRY(cudaMallocHost((void **)&h_candidate_overflow, sizeof(uint32)))
		CUDA_TRY(cudaMallocHost((void **)&h_dedup_cnt, sizeof(uint32)))
		CUDA_TRY(cudaMallocHost((void **)&h_value_match_total, sizeof(uint32)))
		CUDA_TRY(cudaMallocHost((void **)&h_value_overflow, sizeof(uint32)))
		CUDA_TRY(cudaMalloc((void **)&d_filter_iters_hist,
				102u * sizeof(uint32)))
		CUDA_TRY(cudaEventCreateWithFlags(&start_event,
				cudaEventBlockingSync))
		CUDA_TRY(cudaEventCreateWithFlags(&end_event,
				cudaEventBlockingSync))
	}

	~collision_engine() {
		free_device();
		cudaFree(d_filter_iters_hist);
		cudaEventDestroy(start_event);
		cudaEventDestroy(end_event);
		cudaFreeHost(h_max_bucket);
		cudaFreeHost(h_bucket_overflow);
		cudaFreeHost(h_candidate_overflow);
		cudaFreeHost(h_dedup_cnt);
		cudaFreeHost(h_value_match_total);
		cudaFreeHost(h_value_overflow);
	}

	void free_device() {
		cudaFree(d_bucket_count);
		cudaFree(d_bucket_overflow);
		cudaFree(d_candidate_overflow);
		cudaFree(d_arr_a);
		cudaFree(d_arr_b);
		cudaFree(d_candidate_keys);
		cudaFree(d_candidate_cnt);
		cudaFree(d_sorted_keys);
		cudaFree(d_dedup_keys);
		cudaFree(d_dedup_cnt);
		cudaFree(d_D);
		cudaFree(d_D_scan);
		cudaFree(d_D_pos);
		cudaFree(d_S);
		cudaFree(d_X);
		cudaFree(d_max_bucket);
		cudaFree(d_value_counts);
		cudaFree(d_value_offsets);
		cudaFree(d_value_cursor);
		cudaFree(d_matched_values);
		cudaFree(d_value_match_total);
		cudaFree(d_value_overflow);
		cudaFree(d_cub_max_buckets);
		cudaFree(d_cub_sort_candidates);
		cudaFree(d_cub_scan_D);
		cudaFree(d_cub_scan_values);

		d_bucket_count = d_bucket_overflow = d_candidate_overflow = 0;
		d_arr_a = d_arr_b = d_candidate_keys = d_sorted_keys = 0;
		d_dedup_keys = d_X = 0;
		d_candidate_cnt = d_dedup_cnt = d_D = d_D_scan = d_D_pos = 0;
		d_S = d_max_bucket = d_value_counts = d_value_offsets = 0;
		d_value_cursor = d_matched_values = d_value_match_total = 0;
		d_value_overflow = 0;
		d_cub_max_buckets = d_cub_sort_candidates = 0;
		d_cub_scan_D = d_cub_scan_values = 0;
		max_n = max_key_bits = max_per_bucket = 0;
		max_buckets_bytes = sort_candidates_bytes = 0;
		scan_D_bytes = scan_values_bytes = 0;
	}

	void ensure_capacity(uint32 n, uint32 key_bits, uint32 min_per_bucket) {
		if (n <= max_n && key_bits <= max_key_bits &&
		    min_per_bucket <= max_per_bucket)
			return;

		uint32 alloc_n = n > max_n ? n : max_n;
		uint32 alloc_key_bits = key_bits > max_key_bits ?
				key_bits : max_key_bits;
		uint32 min_bucket = min_per_bucket > max_per_bucket ?
				min_per_bucket : max_per_bucket;

		free_device();
		max_n = alloc_n;
		max_key_bits = alloc_key_bits;

		uint32 mean_per_bucket = (max_n + NUM_BUCKETS - 1u) / NUM_BUCKETS;
		uint32 sigma = (uint32)sqrt((double)mean_per_bucket + 1.0);
		uint32 max_bucket_est = mean_per_bucket + 6u * sigma + 32u;
		max_per_bucket = max_bucket_est > min_bucket ?
				max_bucket_est : min_bucket;
		size_t bucket_words = (size_t)NUM_BUCKETS * max_per_bucket;

		CUDA_TRY(cudaMalloc((void **)&d_bucket_count,
				NUM_BUCKETS * sizeof(uint32)))
		CUDA_TRY(cudaMalloc((void **)&d_bucket_overflow, sizeof(uint32)))
		CUDA_TRY(cudaMalloc((void **)&d_candidate_overflow, sizeof(uint32)))
		CUDA_TRY(cudaMalloc((void **)&d_arr_a,
				bucket_words * sizeof(uint64)))
		CUDA_TRY(cudaMalloc((void **)&d_arr_b,
				bucket_words * sizeof(uint64)))
		CUDA_TRY(cudaMalloc((void **)&d_candidate_keys,
				CANDIDATE_CAP * sizeof(uint64)))
		CUDA_TRY(cudaMalloc((void **)&d_candidate_cnt, sizeof(uint32)))
		CUDA_TRY(cudaMalloc((void **)&d_sorted_keys,
				CANDIDATE_CAP * sizeof(uint64)))
		CUDA_TRY(cudaMalloc((void **)&d_dedup_keys,
				CANDIDATE_CAP * sizeof(uint64)))
		CUDA_TRY(cudaMalloc((void **)&d_dedup_cnt, sizeof(uint32)))
		CUDA_TRY(cudaMalloc((void **)&d_D, MAX_DSIZE * sizeof(uint32)))
		CUDA_TRY(cudaMalloc((void **)&d_D_scan, MAX_DSIZE * sizeof(uint32)))
		CUDA_TRY(cudaMalloc((void **)&d_D_pos, MAX_DSIZE * sizeof(uint32)))
		CUDA_TRY(cudaMalloc((void **)&d_S, MAX_SSIZE * sizeof(uint32)))
		CUDA_TRY(cudaMalloc((void **)&d_X,
				CANDIDATE_CAP * sizeof(uint64)))
		CUDA_TRY(cudaMalloc((void **)&d_max_bucket, sizeof(uint32)))
		CUDA_TRY(cudaMalloc((void **)&d_value_counts,
				(VALUE_MATCH_CAP + 1u) * sizeof(uint32)))
		CUDA_TRY(cudaMalloc((void **)&d_value_offsets,
				(VALUE_MATCH_CAP + 1u) * sizeof(uint32)))
		CUDA_TRY(cudaMalloc((void **)&d_value_cursor,
				(VALUE_MATCH_CAP + 1u) * sizeof(uint32)))
		CUDA_TRY(cudaMalloc((void **)&d_matched_values,
				VALUE_MATCH_CAP * sizeof(uint32)))
		CUDA_TRY(cudaMalloc((void **)&d_value_match_total, sizeof(uint32)))
		CUDA_TRY(cudaMalloc((void **)&d_value_overflow, sizeof(uint32)))

		cub::DeviceReduce::Max(0, max_buckets_bytes,
				d_bucket_count, d_max_bucket, NUM_BUCKETS);
		CUDA_TRY(cudaMalloc(&d_cub_max_buckets, max_buckets_bytes))

		cub::DeviceRadixSort::SortKeys(0, sort_candidates_bytes,
				d_candidate_keys, d_sorted_keys, CANDIDATE_CAP);
		CUDA_TRY(cudaMalloc(&d_cub_sort_candidates,
				sort_candidates_bytes))

		cub::DeviceScan::ExclusiveSum(0, scan_D_bytes,
				d_D, d_D_scan, MAX_DSIZE);
		CUDA_TRY(cudaMalloc(&d_cub_scan_D, scan_D_bytes))

		cub::DeviceScan::ExclusiveSum(0, scan_values_bytes,
				d_value_counts, d_value_offsets,
				VALUE_MATCH_CAP + 1u);
		CUDA_TRY(cudaMalloc(&d_cub_scan_values, scan_values_bytes))
	}

	uint32 max_n;
	uint32 max_key_bits;
	uint32 max_per_bucket;

	uint32 *d_bucket_count;
	uint32 *d_bucket_overflow;
	uint32 *d_candidate_overflow;
	uint64 *d_arr_a;
	uint64 *d_arr_b;
	uint64 *d_candidate_keys;
	uint32 *d_candidate_cnt;
	uint64 *d_sorted_keys;
	uint64 *d_dedup_keys;
	uint32 *d_dedup_cnt;
	uint32 *d_D;
	uint32 *d_D_scan;
	uint32 *d_D_pos;
	uint32 *d_S;
	uint64 *d_X;
	uint32 *d_max_bucket;
	uint32 *d_value_counts;
	uint32 *d_value_offsets;
	uint32 *d_value_cursor;
	uint32 *d_matched_values;
	uint32 *d_value_match_total;
	uint32 *d_value_overflow;
	uint32 *d_filter_iters_hist;  /* 102 slots: 3×21 stop-iter + 3×13 size */

	void *d_cub_max_buckets;
	size_t max_buckets_bytes;
	void *d_cub_sort_candidates;
	size_t sort_candidates_bytes;
	void *d_cub_scan_D;
	size_t scan_D_bytes;
	cudaEvent_t start_event;
	cudaEvent_t end_event;
	void *d_cub_scan_values;
	size_t scan_values_bytes;

	uint32 *h_max_bucket;
	uint32 *h_bucket_overflow;
	uint32 *h_candidate_overflow;
	uint32 *h_dedup_cnt;
	uint32 *h_value_match_total;
	uint32 *h_value_overflow;
};

static uint32
host_ilog2(uint32 cnt) {
	if (cnt <= 1)
		return 5;
	uint32 lg = 31u - __builtin_clz(cnt);
	uint32 il = lg + 5u;
	return il < 5u ? 5u : il;
}

static void
collision_stats_clear(collision_data_t *data) {
	data->bucket_max = 0;
	data->candidate_count = 0;
	data->dedup_count = 0;
	data->value_match_count = 0;
	data->bucket_grow_count = 0;
	data->hash_cap_count = 0;
	data->match_arena_attempt_count = 0;
	data->match_arena_fallback_count = 0;
	data->match_arena_capacity_skip_count = 0;
	for (int hi = 0; hi < 102; hi++)
		data->filter_iters_hist[hi] = 0;
	data->elapsed_ms = 0.0f;
}

static void
collision_stats_begin(collision_engine *engine, collision_data_t *data,
		cudaStream_t stream) {

	collision_stats_clear(data);
	if (data->collect_stats)
		CUDA_TRY(cudaEventRecord(engine->start_event, stream))
}

static void
collision_stats_finish(collision_engine *engine, collision_data_t *data,
		cudaStream_t stream) {

	if (data->collect_stats) {
		CUDA_TRY(cudaEventRecord(engine->end_event, stream))
		CUDA_TRY(cudaEventSynchronize(engine->end_event))
		CUDA_TRY(cudaEventElapsedTime(&data->elapsed_ms,
				engine->start_event, engine->end_event))
	}
}

extern "C"
{

COLLISION_ENGINE_DECL void *
collision_engine_init(void) {
	return new collision_engine;
}

COLLISION_ENGINE_DECL void
collision_engine_free(void *e) {
	delete (collision_engine *)e;
}

COLLISION_ENGINE_DECL void
collision_engine_run(void *e, collision_data_t *data) {
	collision_engine *engine = (collision_engine *)e;
	uint32 n = (uint32)data->num_elements;
	uint32 key_bits = (uint32)data->key_bits;
	cudaStream_t stream = (cudaStream_t)data->stream;

	if (data->num_elements > 0xFFFFFFFFu ||
	    key_bits < COLLISION_ENGINE_MIN_KEY_BITS ||
	    key_bits > 64 ||
		    (data->root_bytes != sizeof(uint32) &&
		     data->root_bytes != sizeof(uint64))) {
		printf("collision_engine: invalid input\n");
		exit(-1);
	}
	collision_stats_begin(engine, data, stream);

	dim3 scatter_block(256);
	dim3 scatter_grid((n + scatter_block.x - 1u) / scatter_block.x);

	engine->ensure_capacity(n, key_bits, 0);
	for (;;) {
		CUDA_TRY(cudaMemsetAsync(engine->d_bucket_count, 0,
				NUM_BUCKETS * sizeof(uint32), stream))
		CUDA_TRY(cudaMemsetAsync(engine->d_bucket_overflow, 0,
				sizeof(uint32), stream))

		scatter_roots_kernel<<<scatter_grid, scatter_block, 0, stream>>>(
				(const void *)(size_t)data->keys_in,
				(uint32)data->root_bytes, n, engine->max_per_bucket,
				data->bucket_hash, engine->d_bucket_count,
				engine->d_bucket_overflow, engine->d_arr_a);
		CUDA_CHECK_LAST();

		cub::DeviceReduce::Max(engine->d_cub_max_buckets,
				engine->max_buckets_bytes, engine->d_bucket_count,
				engine->d_max_bucket, NUM_BUCKETS, stream);
		CUDA_TRY(cudaMemcpyAsync(engine->h_max_bucket,
				engine->d_max_bucket, sizeof(uint32),
				cudaMemcpyDeviceToHost, stream))
		CUDA_TRY(cudaMemcpyAsync(engine->h_bucket_overflow,
				engine->d_bucket_overflow, sizeof(uint32),
				cudaMemcpyDeviceToHost, stream))
		CUDA_TRY(cudaStreamSynchronize(stream))

		if (!*engine->h_bucket_overflow &&
		    *engine->h_max_bucket <= engine->max_per_bucket)
			break;

		uint32 observed = *engine->h_max_bucket;
		if (observed <= engine->max_per_bucket)
			observed = engine->max_per_bucket + 1u;
		uint32 grown_cap = observed + observed / 4u + 64u;
		data->bucket_grow_count++;
		if (data->debug) {
			printf("collision_engine: growing bucket cap from %u to %u "
			       "(observed max %u)\n",
			       engine->max_per_bucket, grown_cap,
			       *engine->h_max_bucket);
		}
		engine->ensure_capacity(n, key_bits, grown_cap);
	}
	data->bucket_max = *engine->h_max_bucket;

	uint32 max_bucket_for_hash = *engine->h_max_bucket;
	if (max_bucket_for_hash == 0)
		max_bucket_for_hash = 1;
	uint32 lg = 31u - __builtin_clz(max_bucket_for_hash);
	uint32 ilog2 = lg + 5u;
	if (ilog2 < 5u)
		ilog2 = 5u;
	uint32 hash_bits_avail = key_bits - LOG2_NUM_BUCKETS;
	if (ilog2 > hash_bits_avail)
		ilog2 = hash_bits_avail;
	uint32 max_tsize_words = 1u << (ilog2 - 5u);
	int max_optin_smem = 0;
	CUDA_TRY(cudaDeviceGetAttribute(&max_optin_smem,
			cudaDevAttrMaxSharedMemoryPerBlockOptin, 0))
	if (max_optin_smem > 0) {
		uint32 max_words = (uint32)max_optin_smem /
				(3u * (uint32)sizeof(uint32));
		uint32 capped_words = 1u;

		while ((capped_words << 1) <= max_words)
			capped_words <<= 1;
		if (max_tsize_words > capped_words) {
			data->hash_cap_count++;
			if (data->debug) {
				printf("collision_engine: capping hash words "
				       "from %u to %u (bucket max %u)\n",
				       max_tsize_words, capped_words,
				       max_bucket_for_hash);
			}
			max_tsize_words = capped_words;
		}
	}

	CUDA_TRY(cudaMemsetAsync(engine->d_candidate_cnt, 0,
			sizeof(uint32), stream))
	CUDA_TRY(cudaMemsetAsync(engine->d_candidate_overflow, 0,
			sizeof(uint32), stream))

	size_t hash_bytes = 3u * max_tsize_words * sizeof(uint32);
	CUDA_TRY(cudaFuncSetAttribute(filter_per_bucket_kernel,
			cudaFuncAttributeMaxDynamicSharedMemorySize,
			(int)hash_bytes))

	uint32 *hist_dev = NULL;
	if (data->collect_stats) {
		hist_dev = engine->d_filter_iters_hist;
		CUDA_TRY(cudaMemsetAsync(hist_dev, 0,
				102u * sizeof(uint32), stream))
	}

	filter_per_bucket_kernel<<<NUM_BUCKETS, BLOCK_THREADS,
			hash_bytes, stream>>>(engine->d_bucket_count,
			engine->max_per_bucket, engine->d_arr_a, engine->d_arr_b,
			engine->d_candidate_keys, engine->d_candidate_cnt,
			engine->d_candidate_overflow, CANDIDATE_CAP,
			key_bits, max_tsize_words, hist_dev);
	CUDA_CHECK_LAST();

	uint32 cand_cnt = 0;
	CUDA_TRY(cudaMemcpyAsync(&cand_cnt, engine->d_candidate_cnt,
			sizeof(uint32), cudaMemcpyDeviceToHost, stream))
	CUDA_TRY(cudaMemcpyAsync(engine->h_candidate_overflow,
			engine->d_candidate_overflow, sizeof(uint32),
			cudaMemcpyDeviceToHost, stream))
	if (hist_dev) {
		/* Single contiguous DtoH for all 102 slots — 3×21 stop-iter
		   segments followed by 3×13 size segments. Layout matches
		   the kernel's index scheme. */
		CUDA_TRY(cudaMemcpyAsync(data->filter_iters_hist, hist_dev,
				102u * sizeof(uint32),
				cudaMemcpyDeviceToHost, stream))
	}
	CUDA_TRY(cudaStreamSynchronize(stream))
	data->candidate_count = cand_cnt;

	if (*engine->h_candidate_overflow || cand_cnt > CANDIDATE_CAP) {
		printf("collision_engine: candidate overflow %u > %u\n",
				cand_cnt, CANDIDATE_CAP);
		exit(-1);
	}
	if (cand_cnt == 0) {
		collision_stats_finish(engine, data, stream);
		return;
	}

	cub::DeviceRadixSort::SortKeys(engine->d_cub_sort_candidates,
			engine->sort_candidates_bytes, engine->d_candidate_keys,
			engine->d_sorted_keys, cand_cnt, 0,
			sizeof(uint64) * 8, stream);

	CUDA_TRY(cudaMemsetAsync(engine->d_dedup_cnt, 0,
			sizeof(uint32), stream))
	dim3 block256(256);
	dim3 dedup_grid((cand_cnt + block256.x - 1u) / block256.x);
	dedup_kernel<<<dedup_grid, block256, 0, stream>>>(
			engine->d_sorted_keys, cand_cnt,
			engine->d_dedup_keys, engine->d_dedup_cnt);
	CUDA_CHECK_LAST();
	CUDA_TRY(cudaMemcpyAsync(engine->h_dedup_cnt, engine->d_dedup_cnt,
			sizeof(uint32), cudaMemcpyDeviceToHost, stream))
	CUDA_TRY(cudaStreamSynchronize(stream))

	uint32 dedup_cnt = *engine->h_dedup_cnt;
	data->dedup_count = dedup_cnt;
	if (dedup_cnt == 0) {
		collision_stats_finish(engine, data, stream);
		return;
	}

	uint32 c_ilog2 = host_ilog2(dedup_cnt + 1u);
	if (c_ilog2 > MAX_C_ILOG2)
		c_ilog2 = MAX_C_ILOG2;
	uint32 hash_value = (1u << c_ilog2) - 1u;
	uint32 dsize = (1u << c_ilog2) + 1u;
	uint32 ssize = (1u << c_ilog2) / 32u + 1u;

	CUDA_TRY(cudaMemsetAsync(engine->d_D, 0,
			dsize * sizeof(uint32), stream))
	CUDA_TRY(cudaMemsetAsync(engine->d_S, 0,
			ssize * sizeof(uint32), stream))
	dim3 dedup_count_grid((dedup_cnt + block256.x - 1u) / block256.x);
	count_secondary_kernel<<<dedup_count_grid, block256, 0, stream>>>(
			engine->d_dedup_keys, dedup_cnt, hash_value,
			engine->d_D, engine->d_S);
	CUDA_CHECK_LAST();

	cub::DeviceScan::ExclusiveSum(engine->d_cub_scan_D,
			engine->scan_D_bytes, engine->d_D, engine->d_D_scan,
			dsize, stream);
	CUDA_TRY(cudaMemcpyAsync(engine->d_D_pos, engine->d_D_scan,
			dsize * sizeof(uint32), cudaMemcpyDeviceToDevice, stream))
	scatter_secondary_kernel<<<dedup_count_grid, block256, 0, stream>>>(
			engine->d_dedup_keys, dedup_cnt, hash_value,
			engine->d_D_pos, engine->d_X);
	CUDA_CHECK_LAST();

	CUDA_TRY(cudaMemsetAsync(engine->d_value_counts, 0,
			(dedup_cnt + 1u) * sizeof(uint32), stream))
	CUDA_TRY(cudaMemsetAsync(engine->d_value_match_total, 0,
			sizeof(uint32), stream))
	CUDA_TRY(cudaMemsetAsync(engine->d_value_overflow, 0,
			sizeof(uint32), stream))
	dim3 scan_grid((n + block256.x - 1u) / block256.x);
	uint32 use_match_arena = dedup_cnt <=
			VALUE_MATCH_CAP / MATCH_ARENA_WIDTH;
	data->match_arena_attempt_count = use_match_arena ? 1u : 0u;
	data->match_arena_capacity_skip_count = use_match_arena ? 0u : 1u;
	if (use_match_arena) {
		count_and_store_matched_values_kernel<<<scan_grid, block256,
				0, stream>>>(
				(const void *)(size_t)data->keys_in,
				(const uint32 *)(size_t)data->data_in,
				(uint32)data->root_bytes, n, hash_value,
				engine->d_S, engine->d_D_scan, engine->d_X,
				dedup_cnt, engine->d_value_counts,
				engine->d_value_match_total,
				engine->d_value_overflow,
				engine->d_matched_values, VALUE_MATCH_CAP);
	}
	else {
		count_matched_values_kernel<<<scan_grid, block256, 0, stream>>>(
				(const void *)(size_t)data->keys_in,
				(const uint32 *)(size_t)data->data_in,
				(uint32)data->root_bytes, n, hash_value,
				engine->d_S, engine->d_D_scan, engine->d_X,
				dedup_cnt, engine->d_value_counts,
				engine->d_value_match_total,
				engine->d_value_overflow, VALUE_MATCH_CAP);
	}
	CUDA_CHECK_LAST();
	CUDA_TRY(cudaMemcpyAsync(engine->h_value_match_total,
			engine->d_value_match_total, sizeof(uint32),
			cudaMemcpyDeviceToHost, stream))
	CUDA_TRY(cudaMemcpyAsync(engine->h_value_overflow,
			engine->d_value_overflow, sizeof(uint32),
			cudaMemcpyDeviceToHost, stream))
	CUDA_TRY(cudaStreamSynchronize(stream))
	data->value_match_count = *engine->h_value_match_total;
	if ((!use_match_arena && *engine->h_value_overflow) ||
	    *engine->h_value_match_total > VALUE_MATCH_CAP) {
		printf("collision_engine: value-match overflow %u > %u\n",
				*engine->h_value_match_total, VALUE_MATCH_CAP);
		exit(-1);
	}
	if (*engine->h_value_match_total == 0) {
		collision_stats_finish(engine, data, stream);
		return;
	}
	if (use_match_arena && !*engine->h_value_overflow) {
		emit_found_arena_kernel<<<dedup_count_grid, block256, 0,
				stream>>>(
				engine->d_X, engine->d_matched_values,
				engine->d_value_counts, dedup_cnt, data->shift,
				(uint32)data->root_bytes,
				(specialq_t *)(size_t)data->q_batch,
				(found_t *)(size_t)data->found_array);
		CUDA_CHECK_LAST();
		collision_stats_finish(engine, data, stream);
		return;
	}
	if (data->debug && use_match_arena) {
		printf("collision_engine: match arena overflow "
		       "(dedup %u matched %u width %u), using fallback\n",
		       dedup_cnt, *engine->h_value_match_total,
		       MATCH_ARENA_WIDTH);
	}
	if (use_match_arena)
		data->match_arena_fallback_count = 1u;

	cub::DeviceScan::ExclusiveSum(engine->d_cub_scan_values,
			engine->scan_values_bytes, engine->d_value_counts,
			engine->d_value_offsets, dedup_cnt + 1u, stream);
	CUDA_TRY(cudaMemcpyAsync(engine->d_value_cursor,
			engine->d_value_offsets, dedup_cnt * sizeof(uint32),
			cudaMemcpyDeviceToDevice, stream))
	scatter_matched_values_kernel<<<scan_grid, block256, 0, stream>>>(
			(const void *)(size_t)data->keys_in,
			(const uint32 *)(size_t)data->data_in,
			(uint32)data->root_bytes, n, hash_value,
			engine->d_S, engine->d_D_scan, engine->d_X,
			dedup_cnt, engine->d_value_cursor,
			engine->d_matched_values, VALUE_MATCH_CAP);
	CUDA_CHECK_LAST();

	emit_found_kernel<<<dedup_count_grid, block256, 0, stream>>>(
			engine->d_X, engine->d_matched_values,
			engine->d_value_offsets, dedup_cnt, data->shift,
			(uint32)data->root_bytes,
			(specialq_t *)(size_t)data->q_batch,
			(found_t *)(size_t)data->found_array);
	CUDA_CHECK_LAST();
	collision_stats_finish(engine, data, stream);
}

} // extern "C"
