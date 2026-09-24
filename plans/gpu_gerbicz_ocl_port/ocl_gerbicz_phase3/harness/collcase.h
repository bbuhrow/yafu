/* collcase.h -- engine-neutral test-case format for the gpu_gerbicz
 * collision engine (CUDA) and its planned OpenCL port.
 *
 * ============================================================
 * collcase v1 binary layout (all integers little-endian)
 * ============================================================
 *
 * HEADER (64 bytes, fixed size, written/read field-by-field --
 * never as a raw struct -- so it never depends on compiler padding):
 *
 *   offset  size  field
 *   ------  ----  -----
 *   0       8     magic[8]        = "COLLCAS1" (ASCII, no NUL)
 *   8       4     version         = 1
 *   12      4     n               number of (key,value) elements
 *   16      4     root_bytes      4 or 8 -- width of each key
 *   20      4     key_bits        valid signed range is +/-2^(key_bits-1)
 *   24      4     shift           value = (q_index << shift) | p
 *   28      4     bucket_hash     0 = mask hash, 1 = multiplicative hash
 *   32      4     num_q           number of specialq_t entries
 *   36      4     hash_word_cap   0 = derive from a 64 KiB local-mem
 *                                 limit (4096 words); else explicit
 *                                 override, for phase 5/7 stats-model
 *                                 experiments
 *   40      4     has_expected    0 or 1 -- whether the optional
 *                                 expected-output section follows
 *   44      20    reserved        five u32 zero words, for v2 growth
 *
 * ARRAYS (immediately after the header, no padding between
 * sections):
 *
 *   keys[n]        each root_bytes bytes (u32 or u64 LE). This is
 *                  the RAW bit pattern the engine treats as the sort/
 *                  hash key -- for root_bytes==4 it is the two's-
 *                  complement encoding of a signed value in a 32-bit
 *                  word (NOT sign-extended to 64 bits at this stage:
 *                  the engine zero-extends the 32-bit pattern when it
 *                  builds the internal uint64 sort key, and only
 *                  sign-extends back to a signed 64-bit "offset" at
 *                  emit time). A key of all-zero bytes is the engine's
 *                  reserved "empty" sentinel and is always skipped.
 *
 *   values[n]      each 4 bytes, u32 LE. Packed as
 *                  (q_index << shift) | p, where p occupies the low
 *                  `shift` bits and q_index the rest.
 *
 *   q_batch[num_q] each 24 bytes, matching specialq_t exactly:
 *                    u32 p; u32 pad(=0); u64 pp; u64 root;
 *
 * OPTIONAL EXPECTED-OUTPUT SECTION (present iff has_expected == 1):
 *
 *   found_count     u32  total attempted-store count -- the value
 *                        found_array[0].p1 would reach on the real
 *                        engine, i.e. every emitted (p1,p2) pair,
 *                        with NO 999-entry cap applied.
 *   entry_count     u32  number of entries that follow (== found_count
 *                        for a CPU reference file, since nothing here
 *                        enforces FOUND_ARRAY_SIZE - 1 = 999; a
 *                        saturated comparison against a real engine's
 *                        999-capped dump is the comparator's job, not
 *                        this file's).
 *   entries[entry_count]  each 32 bytes, matching found_t exactly:
 *                    u32 p1; u32 p2; u32 q; u32 pad(=0);
 *                    u64 qroot; i64 offset;
 *   stats block (fixed 424 bytes):
 *     candidate_count      u32
 *     dedup_count          u32
 *     value_match_count    u32
 *     bucket_max           u32
 *     filter_iters_hist[102]  u32 each (408 bytes) -- see
 *                             collision_engine.h for the 3x21 + 3x13
 *                             segment layout.
 *
 * Canonical entry order is NOT fixed by this file format -- the
 * comparator canonicalizes before comparing (see cmp.c).
 */
#ifndef COLLHARNESS_COLLCASE_H
#define COLLHARNESS_COLLCASE_H

#include <stdint.h>
#include <stddef.h>

#define COLLCASE_MAGIC       "COLLCAS1"
#define COLLCASE_MAGIC_LEN   8
#define COLLCASE_VERSION     1u
#define COLLCASE_HEADER_SIZE 64u
#define COLLCASE_HIST_SLOTS  102u
#define COLLCASE_STATS_SIZE  (4u * 4u + COLLCASE_HIST_SLOTS * 4u) /* 424 */

/* Mirrors of the project's stage1_core.h structs, kept local so this
 * harness has zero dependency on CUDA/OpenCL headers. Field order and
 * types must track stage1_core.h exactly. */
typedef struct {
	uint32_t p;
	uint32_t pad;
	uint64_t pp;
	uint64_t root;
} cc_specialq_t;

typedef struct {
	uint32_t p1;
	uint32_t p2;
	uint32_t q;
	uint32_t pad;
	uint64_t qroot;
	int64_t  offset;
} cc_found_t;

/* Static asserts (C99-compatible: no _Static_assert). A negative
 * array size is a compile error on every C99/C++ compiler including
 * MSVC's cl. */
typedef char cc_assert_specialq_is_24_bytes[
		(sizeof(cc_specialq_t) == 24) ? 1 : -1];
typedef char cc_assert_found_is_32_bytes[
		(sizeof(cc_found_t) == 32) ? 1 : -1];

typedef struct {
	uint32_t candidate_count;
	uint32_t dedup_count;
	uint32_t value_match_count;
	uint32_t bucket_max;
	uint32_t filter_iters_hist[COLLCASE_HIST_SLOTS];
} cc_stats_t;

typedef struct {
	uint32_t n;
	uint32_t root_bytes;   /* 4 or 8 */
	uint32_t key_bits;
	uint32_t shift;
	uint32_t bucket_hash;  /* 0 or 1 */
	uint32_t num_q;
	uint32_t hash_word_cap; /* 0 = derive default */

	/* keys, as raw zero-extended 64-bit bit patterns regardless of
	 * root_bytes (top bits are always 0 for root_bytes==4 entries in
	 * memory; only the on-disk encoding trims to 4 bytes). */
	uint64_t *keys;
	uint32_t *values;
	cc_specialq_t *q_batch;

	int has_expected;
	uint32_t found_count;
	uint32_t entry_count;
	cc_found_t *entries;    /* length entry_count, may be NULL */
	cc_stats_t stats;       /* valid iff has_expected */
} collcase_t;

void collcase_init(collcase_t *c);
void collcase_free(collcase_t *c);

/* Returns 0 on success, nonzero (with a message on stderr) on error. */
int collcase_write(const char *path, const collcase_t *c);
int collcase_read(const char *path, collcase_t *c);

#endif /* COLLHARNESS_COLLCASE_H */
