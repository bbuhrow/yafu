#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "collcase.h"
#include "cpuref.h"

typedef struct {
	uint64_t key;
	uint32_t idx;
} kv_t;

static int
cmp_kv(const void *a, const void *b)
{
	const kv_t *x = (const kv_t *)a;
	const kv_t *y = (const kv_t *)b;
	if (x->key < y->key) return -1;
	if (x->key > y->key) return 1;
	if (x->idx < y->idx) return -1;
	if (x->idx > y->idx) return 1;
	return 0;
}

static uint32_t
gcd32(uint32_t a, uint32_t b)
{
	while (b != 0) {
		uint32_t t = a % b;
		a = b;
		b = t;
	}
	return a;
}

static int64_t
sign_extend_key(uint64_t raw, uint32_t root_bytes)
{
	if (root_bytes == 4)
		return (int64_t)(int32_t)(uint32_t)raw;
	return (int64_t)raw;
}

typedef struct {
	cc_found_t *data;
	size_t cap, len;
} found_vec_t;

static int
found_vec_push(found_vec_t *v, const cc_found_t *e)
{
	if (v->len == v->cap) {
		size_t ncap = v->cap ? v->cap * 2 : 1024;
		cc_found_t *nd = (cc_found_t *)realloc(v->data,
				ncap * sizeof(cc_found_t));
		if (!nd)
			return 1;
		v->data = nd;
		v->cap = ncap;
	}
	v->data[v->len++] = *e;
	return 0;
}

int
cpuref_exact(const collcase_t *c, cc_found_t **out_entries,
		uint32_t *out_entry_count, uint32_t *out_found_count)
{
	kv_t *kv;
	uint32_t n = c->n;
	uint32_t nz = 0, i;
	uint32_t mask = (c->shift >= 32) ? 0xFFFFFFFFu :
			((1u << c->shift) - 1u);
	found_vec_t fv;

	memset(&fv, 0, sizeof(fv));

	if (c->shift == 0 || c->shift >= 32) {
		fprintf(stderr, "cpuref_exact: shift must be in [1,31] (got %u)\n",
				c->shift);
		return 1;
	}

	kv = (kv_t *)malloc((size_t)n * sizeof(kv_t));
	if (n > 0 && !kv) {
		fprintf(stderr, "cpuref_exact: out of memory\n");
		return 1;
	}

	for (i = 0; i < n; i++) {
		if (c->keys[i] == 0)
			continue;
		kv[nz].key = c->keys[i];
		kv[nz].idx = i;
		nz++;
	}

	qsort(kv, nz, sizeof(kv_t), cmp_kv);

	i = 0;
	while (i < nz) {
		uint32_t run_start = i, run_end;
		while (i < nz && kv[i].key == kv[run_start].key)
			i++;
		run_end = i; /* [run_start, run_end) all share the same key */

		if (run_end - run_start >= 2) {
			int64_t offset = sign_extend_key(kv[run_start].key,
					c->root_bytes);
			uint32_t a;
			for (a = run_start; a < run_end; a++) {
				uint32_t b;
				uint32_t v1 = c->values[kv[a].idx];
				uint32_t q1 = v1 >> c->shift;
				uint32_t p1 = v1 & mask;
				for (b = a + 1; b < run_end; b++) {
					uint32_t v2 = c->values[kv[b].idx];
					uint32_t q2 = v2 >> c->shift;
					uint32_t p2 = v2 & mask;

					if (q1 != q2 || gcd32(p1, p2) != 1u)
						continue;
					if (q1 >= c->num_q) {
						fprintf(stderr, "cpuref_exact: "
							"q_index %u out of range "
							"(num_q=%u)\n",
							q1, c->num_q);
						free(kv);
						free(fv.data);
						return 1;
					}
					{
						cc_found_t e;
						e.p1 = p1;
						e.p2 = p2;
						e.q = c->q_batch[q1].p;
						e.pad = 0;
						e.qroot = c->q_batch[q1].root;
						e.offset = offset;
						if (found_vec_push(&fv, &e)) {
							fprintf(stderr,
								"cpuref_exact: out of memory\n");
							free(kv);
							free(fv.data);
							return 1;
						}
					}
				}
			}
		}
	}

	free(kv);
	*out_entries = fv.data;
	*out_entry_count = (uint32_t)fv.len;
	*out_found_count = (uint32_t)fv.len;
	return 0;
}
