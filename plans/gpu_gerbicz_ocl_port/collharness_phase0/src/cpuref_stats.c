#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "collcase.h"
#include "cpuref.h"
#include "bucket_hash.h"

#define MAX_FILTER_ITERS 20u
#define DEFAULT_HASH_WORD_CAP 4096u

static uint32_t
cc_clz32(uint32_t x)
{
	uint32_t n = 0;
	if (x == 0)
		return 32;
	while (!(x & 0x80000000u)) {
		x <<= 1;
		n++;
	}
	return n;
}

static uint32_t
compute_ilog2(uint32_t cnt, uint32_t key_bits)
{
	uint32_t lg, il, hash_bits_available;
	if (cnt == 0)
		return 5;
	lg = 31u - cc_clz32(cnt);
	il = lg + 5u;
	if (il < 5u)
		il = 5u;
	hash_bits_available = key_bits - LOG2_NUM_BUCKETS;
	if (il > hash_bits_available)
		il = hash_bits_available;
	return il;
}

static uint32_t
compute_capped_ilog2(uint32_t cnt, uint32_t key_bits, uint32_t max_tsize_words)
{
	uint32_t il = compute_ilog2(cnt, key_bits);
	uint32_t max_il = 5u;
	if (max_tsize_words > 1u)
		max_il += 31u - cc_clz32(max_tsize_words);
	if (il > max_il)
		il = max_il;
	return il;
}

static uint32_t
bits_or_get_prev(uint32_t *tbl, uint32_t hv)
{
	uint32_t word = hv >> 5;
	uint32_t bit = 1u << (hv & 31u);
	uint32_t prev = tbl[word];
	tbl[word] = prev | bit;
	return prev;
}

static void
bit_set(uint32_t *tbl, uint32_t hv)
{
	tbl[hv >> 5] |= (1u << (hv & 31u));
}

static int
bit_test(const uint32_t *tbl, uint32_t hv)
{
	return (tbl[hv >> 5] & (1u << (hv & 31u))) != 0;
}

typedef struct {
	uint64_t *data;
	size_t cap, len;
} u64_vec_t;

static int
u64_vec_push(u64_vec_t *v, uint64_t x)
{
	if (v->len == v->cap) {
		size_t ncap = v->cap ? v->cap * 2 : 4096;
		uint64_t *nd = (uint64_t *)realloc(v->data,
				ncap * sizeof(uint64_t));
		if (!nd)
			return 1;
		v->data = nd;
		v->cap = ncap;
	}
	v->data[v->len++] = x;
	return 0;
}

/* Faithful sequential replay of filter_per_bucket_kernel for one
 * bucket. bucket_items[0..cnt0) is the bucket's raw-key population
 * (order irrelevant -- the filter's output multiset does not depend
 * on insertion order, only on the multiset itself; every round's
 * hash tables are built from an OR over the *whole* current item set
 * before being read, so no ordering dependency is introduced).
 * Surviving keys are appended to *candidates_out; hist (may be NULL)
 * accumulates the 102-slot histogram exactly as collision_engine.cu
 * does. */
static int
filter_bucket(const uint64_t *bucket_items, uint32_t cnt0,
		uint32_t key_bits, uint32_t max_tsize_words,
		u64_vec_t *candidates_out, uint32_t *hist)
{
	uint64_t *arr_in, *arr_out;
	uint32_t *T, *T2, *T3;
	uint32_t cnt = cnt0;
	uint32_t ilog2, tsize, hash_value2, my_shift2;
	uint32_t s_nsize[MAX_FILTER_ITERS];
	int emit = 0;
	uint32_t it, j;

	if (cnt0 == 0)
		return 0;

	arr_in = (uint64_t *)malloc((size_t)cnt0 * sizeof(uint64_t));
	arr_out = (uint64_t *)malloc((size_t)cnt0 * sizeof(uint64_t));
	T = (uint32_t *)calloc(max_tsize_words, sizeof(uint32_t));
	T2 = (uint32_t *)calloc(max_tsize_words, sizeof(uint32_t));
	T3 = (uint32_t *)calloc(max_tsize_words, sizeof(uint32_t));
	if (!arr_in || !arr_out || !T || !T2 || !T3) {
		free(arr_in); free(arr_out); free(T); free(T2); free(T3);
		fprintf(stderr, "filter_bucket: out of memory\n");
		return 1;
	}
	memcpy(arr_in, bucket_items, (size_t)cnt0 * sizeof(uint64_t));

	ilog2 = compute_capped_ilog2(cnt, key_bits, max_tsize_words);
	tsize = 1u << (ilog2 - 5u);
	hash_value2 = (ilog2 >= 32u) ? 0xFFFFFFFFu : ((1u << ilog2) - 1u);
	my_shift2 = LOG2_NUM_BUCKETS;

	/* T, T2 already zero from calloc -- pre-pass builds them. */
	for (j = 0; j < cnt; j++) {
		uint64_t item = arr_in[j];
		uint32_t hv = (uint32_t)((item >> my_shift2) & hash_value2);
		uint32_t bit = 1u << (hv & 31u);
		uint32_t prev = bits_or_get_prev(T2, hv);
		if (prev & bit)
			bit_set(T, hv);
	}

	for (it = 0; it < MAX_FILTER_ITERS; it++) {
		uint32_t my_shift = my_shift2;
		uint32_t hash_value = hash_value2;
		uint32_t *U, *U2, *U3;
		uint32_t cnt_out;
		int stop_zero, stop_cap, stop_conv, stop;

		ilog2 = compute_capped_ilog2(cnt, key_bits, max_tsize_words);
		tsize = 1u << (ilog2 - 5u);

		if ((it & 1u) == 0u) { U = T;  U2 = T2; U3 = T3; }
		else                 { U = T3; U2 = T2; U3 = T;  }

		my_shift2 = my_shift + 6u;
		if (my_shift2 + ilog2 > key_bits)
			my_shift2 = LOG2_NUM_BUCKETS;
		hash_value2 = (ilog2 >= 32u) ? 0xFFFFFFFFu : ((1u << ilog2) - 1u);

		memset(U2, 0, (size_t)tsize * sizeof(uint32_t));
		memset(U3, 0, (size_t)tsize * sizeof(uint32_t));

		cnt_out = 0;
		for (j = 0; j < cnt; j++) {
			uint64_t item = arr_in[j];
			uint32_t hv = (uint32_t)((item >> my_shift) & hash_value);
			if (!bit_test(U, hv))
				continue;
			{
				uint32_t hv2 = (uint32_t)((item >> my_shift2) &
						hash_value2);
				uint32_t bit2 = 1u << (hv2 & 31u);
				uint32_t prev2 = bits_or_get_prev(U2, hv2);
				if (prev2 & bit2)
					bit_set(U3, hv2);
			}
			arr_out[cnt_out++] = item;
		}

		s_nsize[it] = cnt_out;
		cnt = cnt_out;
		{ uint64_t *tmp = arr_in; arr_in = arr_out; arr_out = tmp; }

		stop_zero = (cnt == 0);
		stop_cap = !stop_zero && (it == MAX_FILTER_ITERS - 1);
		stop_conv = !stop_zero && !stop_cap && (it >= 3 &&
				s_nsize[it - 3] == s_nsize[it]);
		stop = stop_zero || stop_cap || stop_conv;

		if (stop) {
			if (hist) {
				uint32_t idx = it;
				if (idx > 20u)
					idx = 20u;
				hist[idx]++;
				if (stop_zero)
					hist[21u + idx]++;
				else if (stop_cap)
					hist[42u + idx]++;
				{
					uint32_t size_bin = (cnt0 == 0u) ? 0u :
							(31u - cc_clz32(cnt0));
					uint32_t cat_base;
					if (size_bin > 12u)
						size_bin = 12u;
					cat_base = (it <= 3) ? 63u :
							(it == 4) ? 76u : 89u;
					hist[cat_base + size_bin]++;
				}
			}
			emit = 1;
			break;
		}
	}

	if (emit && cnt > 0) {
		for (j = 0; j < cnt; j++) {
			if (u64_vec_push(candidates_out, arr_in[j])) {
				free(arr_in); free(arr_out);
				free(T); free(T2); free(T3);
				fprintf(stderr, "filter_bucket: out of memory\n");
				return 1;
			}
		}
	}

	free(arr_in); free(arr_out); free(T); free(T2); free(T3);
	return 0;
}

static int
cmp_u64(const void *a, const void *b)
{
	uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
	if (x < y) return -1;
	if (x > y) return 1;
	return 0;
}

int
cpuref_stats(const collcase_t *c, cc_stats_t *out)
{
	uint32_t *bucket_count = NULL;
	uint32_t *bucket_offset = NULL;
	uint64_t *bucket_storage = NULL;
	uint32_t max_tsize_words;
	uint32_t b, i, n = c->n;
	uint64_t match_total = 0;
	u64_vec_t candidates;
	int rc = 0;

	memset(out, 0, sizeof(*out));
	memset(&candidates, 0, sizeof(candidates));

	max_tsize_words = c->hash_word_cap ? c->hash_word_cap :
			DEFAULT_HASH_WORD_CAP;

	bucket_count = (uint32_t *)calloc(NUM_BUCKETS, sizeof(uint32_t));
	bucket_offset = (uint32_t *)calloc((size_t)NUM_BUCKETS + 1,
			sizeof(uint32_t));
	if (!bucket_count || !bucket_offset) {
		fprintf(stderr, "cpuref_stats: out of memory\n");
		rc = 1;
		goto done;
	}

	for (i = 0; i < n; i++) {
		uint64_t k = c->keys[i];
		if (k == 0)
			continue;
		bucket_count[compute_bucket(k, (int)c->bucket_hash)]++;
	}

	for (b = 0; b < NUM_BUCKETS; b++) {
		bucket_offset[b + 1] = bucket_offset[b] + bucket_count[b];
		if (bucket_count[b] > out->bucket_max)
			out->bucket_max = bucket_count[b];
	}

	bucket_storage = (uint64_t *)malloc(
			(size_t)bucket_offset[NUM_BUCKETS] * sizeof(uint64_t));
	if (bucket_offset[NUM_BUCKETS] > 0 && !bucket_storage) {
		fprintf(stderr, "cpuref_stats: out of memory\n");
		rc = 1;
		goto done;
	}
	{
		uint32_t *cursor = (uint32_t *)malloc(
				(size_t)NUM_BUCKETS * sizeof(uint32_t));
		if (!cursor) {
			fprintf(stderr, "cpuref_stats: out of memory\n");
			rc = 1;
			goto done;
		}
		memcpy(cursor, bucket_offset, (size_t)NUM_BUCKETS * sizeof(uint32_t));
		for (i = 0; i < n; i++) {
			uint64_t k = c->keys[i];
			uint32_t bkt;
			if (k == 0)
				continue;
			bkt = compute_bucket(k, (int)c->bucket_hash);
			bucket_storage[cursor[bkt]++] = k;
		}
		free(cursor);
	}

	for (b = 0; b < NUM_BUCKETS; b++) {
		uint32_t pop = bucket_offset[b + 1] - bucket_offset[b];
		if (pop == 0)
			continue;
		if (filter_bucket(bucket_storage + bucket_offset[b], pop,
				c->key_bits, max_tsize_words, &candidates,
				out->filter_iters_hist)) {
			rc = 1;
			goto done;
		}
	}
	out->candidate_count = (uint32_t)candidates.len;

	if (candidates.len == 0)
		goto done;

	qsort(candidates.data, candidates.len, sizeof(uint64_t), cmp_u64);

	/* dedup_count: number of distinct values appearing >=2 times in
	 * the candidate multiset (mirrors dedup_kernel's is_first &&
	 * has_neighbor exact-duplicate test). value_match_count: total
	 * original (pre-filter) occurrences of every such deduped value
	 * -- computed directly against the full original array, which is
	 * mathematically identical to what the engine's exact D/S/X
	 * secondary hash counts (see cpuref.h for why). */
	{
		size_t idx = 0;
		/* index original array by key for O(n log n) multiplicity
		 * lookup, reusing the same sort-by-key approach as 0.4. */
		uint64_t *okeys = (uint64_t *)malloc((size_t)n * sizeof(uint64_t));
		uint32_t nz = 0;
		if (n > 0 && !okeys) {
			fprintf(stderr, "cpuref_stats: out of memory\n");
			rc = 1;
			goto done;
		}
		for (i = 0; i < n; i++) {
			if (c->keys[i] != 0)
				okeys[nz++] = c->keys[i];
		}
		qsort(okeys, nz, sizeof(uint64_t), cmp_u64);

		while (idx < candidates.len) {
			size_t run_start = idx, run_len;
			while (idx < candidates.len &&
					candidates.data[idx] == candidates.data[run_start])
				idx++;
			run_len = idx - run_start;
			if (run_len >= 2) {
				uint64_t val = candidates.data[run_start];
				/* count true multiplicity in the full original
				 * array via binary search + linear scan of the
				 * matching run. */
				size_t lo = 0, hi = nz, mid;
				size_t first, last, mult;
				while (lo < hi) {
					mid = lo + (hi - lo) / 2;
					if (okeys[mid] < val) lo = mid + 1;
					else hi = mid;
				}
				first = lo;
				lo = first; hi = nz;
				while (lo < hi) {
					mid = lo + (hi - lo) / 2;
					if (okeys[mid] <= val) lo = mid + 1;
					else hi = mid;
				}
				last = lo; /* [first,last) all equal val */
				mult = last - first;

				out->dedup_count++;
				match_total += mult;
			}
		}
		free(okeys);
	}
	out->value_match_count = (uint32_t)match_total;

done:
	free(bucket_count);
	free(bucket_offset);
	free(bucket_storage);
	free(candidates.data);
	return rc;
}
