#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "collcase.h"
#include "gen.h"
#include "cpuref.h"
#include "prng.h"
#include "bucket_hash.h"

void
gen_params_defaults(gen_params_t *p)
{
	memset(p, 0, sizeof(*p));
	p->n = 10000;
	p->root_bytes = 4;
	p->key_bits = 24;
	p->shift = 20;
	p->bucket_hash = 0;
	p->hash_word_cap = 0;
	p->seed = 1;
	p->collision_density = 0.01;
	p->bucket_skew_count = 0;
	p->include_edge_cases = 1;
	p->with_expected = 1;
}

/* ---- dynamic vectors ---------------------------------------------- */

typedef struct {
	uint64_t *keys;
	uint32_t *values;
	size_t cap, len;
} kvvec_t;

static int
kv_reserve(kvvec_t *v, size_t more)
{
	if (v->len + more <= v->cap)
		return 0;
	{
		size_t ncap = v->cap ? v->cap * 2 : 1024;
		uint64_t *nk;
		uint32_t *nv;
		while (ncap < v->len + more)
			ncap *= 2;
		nk = (uint64_t *)realloc(v->keys, ncap * sizeof(uint64_t));
		nv = (uint32_t *)realloc(v->values, ncap * sizeof(uint32_t));
		if (!nk || !nv)
			return 1;
		v->keys = nk;
		v->values = nv;
		v->cap = ncap;
	}
	return 0;
}

static int
kv_push(kvvec_t *v, uint64_t k, uint32_t val)
{
	if (kv_reserve(v, 1))
		return 1;
	v->keys[v->len] = k;
	v->values[v->len] = val;
	v->len++;
	return 0;
}

typedef struct {
	cc_specialq_t *data;
	size_t cap, len;
} qvec_t;

static int
q_push(qvec_t *v, uint32_t p, uint64_t root, uint32_t *out_idx)
{
	if (v->len == v->cap) {
		size_t ncap = v->cap ? v->cap * 2 : 256;
		cc_specialq_t *nd = (cc_specialq_t *)realloc(v->data,
				ncap * sizeof(cc_specialq_t));
		if (!nd)
			return 1;
		v->data = nd;
		v->cap = ncap;
	}
	v->data[v->len].p = p;
	v->data[v->len].pad = 0;
	v->data[v->len].pp = (uint64_t)p * (uint64_t)p;
	v->data[v->len].root = root;
	*out_idx = (uint32_t)v->len;
	v->len++;
	return 0;
}

/* ---- small helpers -------------------------------------------------- */

static int
is_prime(uint32_t x)
{
	uint32_t i;
	if (x < 2)
		return 0;
	if (x % 2 == 0)
		return x == 2;
	for (i = 3; (uint64_t)i * i <= x; i += 2) {
		if (x % i == 0)
			return 0;
	}
	return 1;
}

/* Fills out[0..count) with `count` distinct primes >= start, each
 * strictly less than limit (0 = no limit). Any two entries are
 * therefore pairwise coprime by construction. Returns 0 on success. */
static int
fill_prime_pool(uint32_t *out, uint32_t count, uint32_t start, uint32_t limit)
{
	uint32_t found = 0, x = start | 1u; /* start odd search near start */
	if (start <= 2) {
		if (count > 0)
			out[found++] = 2;
		x = 3;
	}
	for (; found < count; x += 2) {
		if (limit && x >= limit) {
			fprintf(stderr, "fill_prime_pool: ran out of primes "
					"below %u\n", limit);
			return 1;
		}
		if (is_prime(x))
			out[found++] = x;
	}
	return 0;
}

static int64_t
key_min(uint32_t key_bits)
{
	if (key_bits >= 64)
		return (int64_t)((uint64_t)1 << 63);
	return -(int64_t)((uint64_t)1 << (key_bits - 1));
}

static int64_t
key_max(uint32_t key_bits)
{
	if (key_bits >= 64)
		return (int64_t)(((uint64_t)1 << 63) - 1);
	return (int64_t)(((uint64_t)1 << (key_bits - 1)) - 1);
}

static uint64_t
encode_key(int64_t v, uint32_t root_bytes)
{
	if (root_bytes == 4)
		return (uint64_t)(uint32_t)(int32_t)v;
	return (uint64_t)v;
}

static int64_t
rand_key(cc_rng_t *r, uint32_t key_bits)
{
	int64_t v;
	do {
		if (key_bits >= 64) {
			v = (int64_t)cc_rng_next(r);
		} else {
			uint64_t span = (uint64_t)1 << key_bits;
			uint64_t off = cc_rng_next(r) % span;
			v = key_min(key_bits) + (int64_t)off;
		}
	} while (v == 0);
	return v;
}

static uint32_t
make_value(uint32_t q_index, uint32_t p, uint32_t shift)
{
	return (q_index << shift) | p;
}

/* ---- edge-case scenarios (plan task 0.3) ---------------------------- */

static int
add_edge_cases(kvvec_t *kv, qvec_t *qv, cc_rng_t *r, const gen_params_t *p,
		const uint32_t *pool_lo, const uint32_t *pool_hi)
{
	uint32_t qi;
	int64_t kmin = key_min(p->key_bits), kmax = key_max(p->key_bits);

	/* 1: basic positive pair -- same key, same q, coprime p -> emits
	 * exactly one pair. */
	if (q_push(qv, 90001u, cc_rng_next(r), &qi))
		return 1;
	{
		int64_t k = -12345;
		if (k < kmin || k > kmax)
			k = kmin + 1;
		if (kv_push(kv, encode_key(k, p->root_bytes),
				make_value(qi, pool_lo[0], p->shift)))
			return 1;
		if (kv_push(kv, encode_key(k, p->root_bytes),
				make_value(qi, pool_lo[1], p->shift)))
			return 1;
	}

	/* 2: same key, different q_index -> must NOT emit. */
	{
		uint32_t qa, qb;
		int64_t k = 424242;
		if (k < kmin || k > kmax)
			k = kmax - 1;
		if (q_push(qv, 90007u, cc_rng_next(r), &qa))
			return 1;
		if (q_push(qv, 90011u, cc_rng_next(r), &qb))
			return 1;
		if (kv_push(kv, encode_key(k, p->root_bytes),
				make_value(qa, pool_lo[2], p->shift)))
			return 1;
		if (kv_push(kv, encode_key(k, p->root_bytes),
				make_value(qb, pool_lo[2], p->shift)))
			return 1;
	}

	/* 3: same key, same q, gcd(p1,p2) > 1 -> must NOT emit. */
	{
		uint32_t q3;
		int64_t k = -999999;
		if (k < kmin || k > kmax)
			k = kmin + 2;
		if (q_push(qv, 90013u, cc_rng_next(r), &q3))
			return 1;
		if (kv_push(kv, encode_key(k, p->root_bytes),
				make_value(q3, 6u, p->shift)))
			return 1;
		if (kv_push(kv, encode_key(k, p->root_bytes),
				make_value(q3, 9u, p->shift)))
			return 1;
	}

	/* 4: multiplicity 10, one shared key/q, pairwise-coprime p's --
	 * exercises the >MATCH_ARENA_WIDTH(8) fallback path. Emits
	 * C(10,2) = 45 pairs. */
	{
		uint32_t q4, i;
		int64_t k = 555555;
		if (k < kmin || k > kmax)
			k = kmax - 2;
		if (q_push(qv, 90017u, cc_rng_next(r), &q4))
			return 1;
		for (i = 0; i < 10; i++) {
			if (kv_push(kv, encode_key(k, p->root_bytes),
					make_value(q4, pool_lo[i], p->shift)))
				return 1;
		}
	}

	/* 5: zero keys -- always skipped, whatever the value. */
	{
		uint32_t i;
		for (i = 0; i < 5; i++) {
			if (kv_push(kv, 0u, make_value(qi, pool_lo[3], p->shift)))
				return 1;
		}
	}

	/* 6: sign / boundary coverage -- min and max representable keys,
	 * each with a coprime duplicate pair. */
	{
		uint32_t q6a, q6b;
		if (q_push(qv, 90019u, cc_rng_next(r), &q6a))
			return 1;
		if (q_push(qv, 90023u, cc_rng_next(r), &q6b))
			return 1;
		if (kv_push(kv, encode_key(kmin, p->root_bytes),
				make_value(q6a, pool_lo[4], p->shift)))
			return 1;
		if (kv_push(kv, encode_key(kmin, p->root_bytes),
				make_value(q6a, pool_lo[5], p->shift)))
			return 1;
		if (kv_push(kv, encode_key(kmax, p->root_bytes),
				make_value(q6b, pool_lo[6], p->shift)))
			return 1;
		if (kv_push(kv, encode_key(kmax, p->root_bytes),
				make_value(q6b, pool_lo[7], p->shift)))
			return 1;
	}

	/* 7: p below and above 65536, each in its own duplicate group. */
	{
		uint32_t q7a, q7b;
		int64_t ka = 111111, kb = -222222;
		if (ka < kmin || ka > kmax) ka = kmin + 3;
		if (kb < kmin || kb > kmax) kb = kmax - 3;
		if (q_push(qv, 90029u, cc_rng_next(r), &q7a))
			return 1;
		if (q_push(qv, 90031u, cc_rng_next(r), &q7b))
			return 1;
		if (kv_push(kv, encode_key(ka, p->root_bytes),
				make_value(q7a, pool_lo[8], p->shift)))
			return 1;
		if (kv_push(kv, encode_key(ka, p->root_bytes),
				make_value(q7a, pool_lo[9], p->shift)))
			return 1;
		if (pool_hi) {
			if (kv_push(kv, encode_key(kb, p->root_bytes),
					make_value(q7b, pool_hi[0], p->shift)))
				return 1;
			if (kv_push(kv, encode_key(kb, p->root_bytes),
					make_value(q7b, pool_hi[1], p->shift)))
				return 1;
		}
	}

	return 0;
}

/* Forces `count` distinct, singleton keys into one target bucket under
 * the case's hash mode, to exercise bucket-population growth in later
 * phases. For hash_mode 0 this is a direct construction; for
 * hash_mode 1 (multiplicative mix) it is a bounded random search. */
static int
add_bucket_skew(kvvec_t *kv, qvec_t *qv, cc_rng_t *r, const gen_params_t *p,
		uint32_t target_bucket, uint32_t count)
{
	uint32_t made = 0, qi;
	uint64_t attempts = 0;
	uint64_t max_attempts = (uint64_t)count * NUM_BUCKETS * 64u + 1000000u;

	if (count == 0)
		return 0;
	if (q_push(qv, 90097u, cc_rng_next(r), &qi))
		return 1;

	if (p->bucket_hash == 0) {
		uint64_t upper = 1;
		while (made < count) {
			uint64_t cand = (upper << LOG2_NUM_BUCKETS) | target_bucket;
			/* keep candidate within the case's key range and
			 * within root_bytes width */
			if (p->root_bytes == 4 && cand > 0xFFFFFFFFu)
				break;
			if (cand == 0) {
				upper++;
				continue;
			}
			if (kv_push(kv, cand, make_value(qi,
					100003u + (made % 97u), p->shift))) {
				return 1;
			}
			made++;
			upper++;
		}
	} else {
		while (made < count && attempts < max_attempts) {
			uint64_t cand;
			attempts++;
			if (p->root_bytes == 4)
				cand = (uint64_t)(uint32_t)cc_rng_next(r);
			else
				cand = cc_rng_next(r);
			if (cand == 0)
				continue;
			if (compute_bucket(cand, 1) != target_bucket)
				continue;
			if (kv_push(kv, cand, make_value(qi,
					100003u + (made % 97u), p->shift))) {
				return 1;
			}
			made++;
		}
		if (made < count) {
			fprintf(stderr, "add_bucket_skew: only found %u/%u "
					"keys for bucket %u within %llu "
					"attempts\n", made, count,
					target_bucket,
					(unsigned long long)max_attempts);
		}
	}
	return 0;
}

/* ---- top level -------------------------------------------------------- */

int
gen_generate(const gen_params_t *p, collcase_t *out)
{
	kvvec_t kv;
	qvec_t qv;
	cc_rng_t r;
	uint32_t pool_lo[16];
	uint32_t pool_hi[16];
	int have_hi = 1;
	uint32_t i;

	memset(&kv, 0, sizeof(kv));
	memset(&qv, 0, sizeof(qv));
	cc_rng_seed(&r, p->seed);

	if (p->root_bytes != 4 && p->root_bytes != 8) {
		fprintf(stderr, "gen_generate: root_bytes must be 4 or 8\n");
		return 1;
	}
	if (p->root_bytes == 4 && p->key_bits > 32) {
		fprintf(stderr, "gen_generate: key_bits > 32 needs root_bytes=8\n");
		return 1;
	}
	if (p->key_bits < 2 || p->key_bits > 64) {
		fprintf(stderr, "gen_generate: key_bits out of range\n");
		return 1;
	}
	if (p->shift == 0 || p->shift >= 32) {
		fprintf(stderr, "gen_generate: shift must be in [1,31]\n");
		return 1;
	}

	if (fill_prime_pool(pool_lo, 16, 2, 65536u)) {
		return 1;
	}
	{
		uint32_t limit = (p->shift >= 32) ? 0u : (1u << p->shift);
		if (limit == 0 || limit > 65536u) {
			if (fill_prime_pool(pool_hi, 16, 65537u, limit))
				have_hi = 0;
		} else {
			have_hi = 0;
		}
	}

	if (p->include_edge_cases) {
		if (add_edge_cases(&kv, &qv, &r, p, pool_lo,
				have_hi ? pool_hi : NULL)) {
			goto fail;
		}
	}

	if (p->bucket_skew_count > 0) {
		uint32_t target_bucket = (uint32_t)cc_rng_below(&r, NUM_BUCKETS);
		if (add_bucket_skew(&kv, &qv, &r, p, target_bucket,
				p->bucket_skew_count)) {
			goto fail;
		}
	}

	/* filler q pool, shared by all random filler entries */
	{
		uint32_t filler_q_count = p->n / 8 + 4;
		uint32_t *filler_q_idx;
		uint32_t remaining;

		if (filler_q_count > 100000u)
			filler_q_count = 100000u;
		filler_q_idx = (uint32_t *)malloc(filler_q_count * sizeof(uint32_t));
		if (!filler_q_idx)
			goto fail;
		for (i = 0; i < filler_q_count; i++) {
			uint32_t qi;
			uint32_t synthetic_p = 200003u + i * 2u;
			if (q_push(&qv, synthetic_p, cc_rng_next(&r), &qi)) {
				free(filler_q_idx);
				goto fail;
			}
			filler_q_idx[i] = qi;
		}

		remaining = (kv.len >= p->n) ? 0 : (p->n - (uint32_t)kv.len);
		while (remaining > 0) {
			uint32_t qidx = filler_q_idx[cc_rng_below(&r, filler_q_count)];
			int64_t k = rand_key(&r, p->key_bits);
			uint64_t ek = encode_key(k, p->root_bytes);
			int make_pair = (remaining >= 2) &&
					(cc_rng_below(&r, 1000000u) <
					 (uint64_t)(p->collision_density * 1000000.0));

			if (make_pair) {
				uint32_t a = cc_rng_below(&r, 16);
				uint32_t b = (a + 1 + cc_rng_below(&r, 15)) % 16;
				if (kv_push(&kv, ek, make_value(qidx, pool_lo[a], p->shift)) ||
				    kv_push(&kv, ek, make_value(qidx, pool_lo[b], p->shift))) {
					free(filler_q_idx);
					goto fail;
				}
				remaining -= 2;
			} else {
				uint32_t a = cc_rng_below(&r, 16);
				if (kv_push(&kv, ek, make_value(qidx, pool_lo[a], p->shift))) {
					free(filler_q_idx);
					goto fail;
				}
				remaining -= 1;
			}
		}
		free(filler_q_idx);
	}

	collcase_init(out);
	out->n = (uint32_t)kv.len;
	out->root_bytes = p->root_bytes;
	out->key_bits = p->key_bits;
	out->shift = p->shift;
	out->bucket_hash = p->bucket_hash;
	out->num_q = (uint32_t)qv.len;
	out->hash_word_cap = p->hash_word_cap;
	out->keys = kv.keys;
	out->values = kv.values;
	out->q_batch = qv.data;

	if (p->with_expected) {
		cc_found_t *entries = NULL;
		uint32_t entry_count = 0, found_count = 0;
		if (cpuref_exact(out, &entries, &entry_count, &found_count)) {
			collcase_free(out);
			return 1;
		}
		if (cpuref_stats(out, &out->stats)) {
			free(entries);
			collcase_free(out);
			return 1;
		}
		out->has_expected = 1;
		out->found_count = found_count;
		out->entry_count = entry_count;
		out->entries = entries;
	}

	return 0;

fail:
	free(kv.keys);
	free(kv.values);
	free(qv.data);
	fprintf(stderr, "gen_generate: out of memory or generation error\n");
	return 1;
}
