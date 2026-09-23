/* cpuref.h -- CPU-only, sequential reference implementations of the
 * gpu_gerbicz collision engine's semantics.
 *
 * Two independent pieces, matching plan tasks 0.4 and 0.5:
 *
 *  - cpuref_exact(): the found SET at result level, computed by a
 *    simple sort/group -- independent of the hash filter entirely.
 *    This is the ground truth every engine (CUDA, OpenCL, this
 *    harness) must agree with.
 *
 *  - cpuref_stats(): a faithful, order-independent sequential replay
 *    of filter_per_bucket_kernel's multi-round hash-collision filter,
 *    to reproduce candidate_count and filter_iters_hist bit-for-bit.
 *    dedup_count and value_match_count are NOT re-derived via a
 *    second hash-table simulation (the D/S/X secondary-hash tables in
 *    collision_engine.cu are an exact, collision-free lookup with no
 *    false positives -- an implementation detail for GPU parallelism,
 *    not a source of approximation) -- they are computed directly
 *    from the exact-duplicate structure of the candidate set the
 *    filter produces, which is mathematically identical to what the
 *    secondary hash computes as long as no capacity/overflow path is
 *    hit (assumed throughout this harness; see STATUS "Open issues").
 */
#ifndef COLLHARNESS_CPUREF_H
#define COLLHARNESS_CPUREF_H

#include "collcase.h"

/* 0.4: exact found-set reference, independent of any hash filter.
 * On success (0), *out_entries is a malloc'd array of *out_entry_count
 * cc_found_t (caller frees), and *out_found_count is the total
 * attempted-store count (== *out_entry_count here, since this
 * reference applies no FOUND_ARRAY_SIZE cap). Returns nonzero and
 * leaves outputs untouched on error (e.g. a value's q_index out of
 * range for num_q). */
int cpuref_exact(const collcase_t *c, cc_found_t **out_entries,
		uint32_t *out_entry_count, uint32_t *out_found_count);

/* 0.5: sequential exact mirror of the filter pipeline. Fills *out
 * with candidate_count, dedup_count, value_match_count, bucket_max,
 * and the 102-slot filter_iters_hist. Uses c->hash_word_cap if
 * nonzero, else the default 4096-word cap (64 KiB / (3 tables *
 * 4 bytes), rounded down to a power of two). Returns nonzero on
 * error. */
int cpuref_stats(const collcase_t *c, cc_stats_t *out);

#endif /* COLLHARNESS_CPUREF_H */
