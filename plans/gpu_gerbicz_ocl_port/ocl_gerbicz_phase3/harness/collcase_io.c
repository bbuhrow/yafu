#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "collcase.h"
#include "byteio.h"

void
collcase_init(collcase_t *c)
{
	memset(c, 0, sizeof(*c));
}

void
collcase_free(collcase_t *c)
{
	free(c->keys);
	free(c->values);
	free(c->q_batch);
	free(c->entries);
	memset(c, 0, sizeof(*c));
}

int
collcase_write(const char *path, const collcase_t *c)
{
	FILE *f;
	uint32_t i;

	if (c->root_bytes != 4 && c->root_bytes != 8) {
		fprintf(stderr, "collcase_write: bad root_bytes %u\n",
				c->root_bytes);
		return 1;
	}

	f = fopen(path, "wb");
	if (!f) {
		fprintf(stderr, "collcase_write: cannot open %s\n", path);
		return 1;
	}

	wr_bytes(f, COLLCASE_MAGIC, COLLCASE_MAGIC_LEN);
	wr_u32(f, COLLCASE_VERSION);
	wr_u32(f, c->n);
	wr_u32(f, c->root_bytes);
	wr_u32(f, c->key_bits);
	wr_u32(f, c->shift);
	wr_u32(f, c->bucket_hash);
	wr_u32(f, c->num_q);
	wr_u32(f, c->hash_word_cap);
	wr_u32(f, c->has_expected ? 1u : 0u);
	for (i = 0; i < 5; i++)
		wr_u32(f, 0u);

	for (i = 0; i < c->n; i++) {
		if (c->root_bytes == 4)
			wr_u32(f, (uint32_t)c->keys[i]);
		else
			wr_u64(f, c->keys[i]);
	}
	for (i = 0; i < c->n; i++)
		wr_u32(f, c->values[i]);
	for (i = 0; i < c->num_q; i++) {
		wr_u32(f, c->q_batch[i].p);
		wr_u32(f, c->q_batch[i].pad);
		wr_u64(f, c->q_batch[i].pp);
		wr_u64(f, c->q_batch[i].root);
	}

	if (c->has_expected) {
		wr_u32(f, c->found_count);
		wr_u32(f, c->entry_count);
		for (i = 0; i < c->entry_count; i++) {
			const cc_found_t *e = &c->entries[i];
			wr_u32(f, e->p1);
			wr_u32(f, e->p2);
			wr_u32(f, e->q);
			wr_u32(f, e->pad);
			wr_u64(f, e->qroot);
			wr_i64(f, e->offset);
		}
		wr_u32(f, c->stats.candidate_count);
		wr_u32(f, c->stats.dedup_count);
		wr_u32(f, c->stats.value_match_count);
		wr_u32(f, c->stats.bucket_max);
		for (i = 0; i < COLLCASE_HIST_SLOTS; i++)
			wr_u32(f, c->stats.filter_iters_hist[i]);
	}

	if (ferror(f)) {
		fprintf(stderr, "collcase_write: write error on %s\n", path);
		fclose(f);
		return 1;
	}
	fclose(f);
	return 0;
}

#define CC_FAIL(msg) do { \
		fprintf(stderr, "collcase_read: %s (%s)\n", msg, path); \
		fclose(f); \
		return 1; \
	} while (0)

int
collcase_read(const char *path, collcase_t *c)
{
	FILE *f;
	char magic[COLLCASE_MAGIC_LEN];
	uint32_t version, has_expected, reserved;
	uint32_t i;

	collcase_init(c);

	f = fopen(path, "rb");
	if (!f) {
		fprintf(stderr, "collcase_read: cannot open %s\n", path);
		return 1;
	}

	if (!rd_bytes(f, magic, COLLCASE_MAGIC_LEN))
		CC_FAIL("truncated magic");
	if (memcmp(magic, COLLCASE_MAGIC, COLLCASE_MAGIC_LEN) != 0)
		CC_FAIL("bad magic (not a collcase v1 file)");
	if (!rd_u32(f, &version))
		CC_FAIL("truncated header");
	if (version != COLLCASE_VERSION)
		CC_FAIL("unsupported version");

	if (!rd_u32(f, &c->n) || !rd_u32(f, &c->root_bytes) ||
	    !rd_u32(f, &c->key_bits) || !rd_u32(f, &c->shift) ||
	    !rd_u32(f, &c->bucket_hash) || !rd_u32(f, &c->num_q) ||
	    !rd_u32(f, &c->hash_word_cap) || !rd_u32(f, &has_expected))
		CC_FAIL("truncated header");
	for (i = 0; i < 5; i++) {
		if (!rd_u32(f, &reserved))
			CC_FAIL("truncated header reserved words");
	}
	c->has_expected = has_expected ? 1 : 0;

	if (c->root_bytes != 4 && c->root_bytes != 8)
		CC_FAIL("bad root_bytes");

	if (c->n > 0) {
		c->keys = (uint64_t *)malloc((size_t)c->n * sizeof(uint64_t));
		c->values = (uint32_t *)malloc((size_t)c->n * sizeof(uint32_t));
		if (!c->keys || !c->values)
			CC_FAIL("out of memory");
	}
	for (i = 0; i < c->n; i++) {
		if (c->root_bytes == 4) {
			uint32_t k32;
			if (!rd_u32(f, &k32))
				CC_FAIL("truncated keys[]");
			c->keys[i] = (uint64_t)k32;
		} else {
			if (!rd_u64(f, &c->keys[i]))
				CC_FAIL("truncated keys[]");
		}
	}
	for (i = 0; i < c->n; i++) {
		if (!rd_u32(f, &c->values[i]))
			CC_FAIL("truncated values[]");
	}

	if (c->num_q > 0) {
		c->q_batch = (cc_specialq_t *)malloc(
				(size_t)c->num_q * sizeof(cc_specialq_t));
		if (!c->q_batch)
			CC_FAIL("out of memory");
	}
	for (i = 0; i < c->num_q; i++) {
		if (!rd_u32(f, &c->q_batch[i].p) ||
		    !rd_u32(f, &c->q_batch[i].pad) ||
		    !rd_u64(f, &c->q_batch[i].pp) ||
		    !rd_u64(f, &c->q_batch[i].root))
			CC_FAIL("truncated q_batch[]");
	}

	if (c->has_expected) {
		if (!rd_u32(f, &c->found_count) || !rd_u32(f, &c->entry_count))
			CC_FAIL("truncated expected-output header");
		if (c->entry_count > 0) {
			c->entries = (cc_found_t *)malloc(
					(size_t)c->entry_count *
					sizeof(cc_found_t));
			if (!c->entries)
				CC_FAIL("out of memory");
		}
		for (i = 0; i < c->entry_count; i++) {
			cc_found_t *e = &c->entries[i];
			if (!rd_u32(f, &e->p1) || !rd_u32(f, &e->p2) ||
			    !rd_u32(f, &e->q) || !rd_u32(f, &e->pad) ||
			    !rd_u64(f, &e->qroot) || !rd_i64(f, &e->offset))
				CC_FAIL("truncated entries[]");
		}
		if (!rd_u32(f, &c->stats.candidate_count) ||
		    !rd_u32(f, &c->stats.dedup_count) ||
		    !rd_u32(f, &c->stats.value_match_count) ||
		    !rd_u32(f, &c->stats.bucket_max))
			CC_FAIL("truncated stats block");
		for (i = 0; i < COLLCASE_HIST_SLOTS; i++) {
			if (!rd_u32(f, &c->stats.filter_iters_hist[i]))
				CC_FAIL("truncated stats histogram");
		}
	}

	fclose(f);
	return 0;
}
