/* cmp.h -- canonicalize and compare two collcase result files. */
#ifndef COLLHARNESS_CMP_H
#define COLLHARNESS_CMP_H

/* Returns 0 if the files match under the saturation-aware rule below,
 * 1 on a genuine mismatch, 2 on a usage/IO error. Prints a short
 * human-readable report to stdout either way.
 *
 * observed_path is the file under test (e.g. a real engine's dump,
 * possibly capped at 999 stored entries per FOUND_ARRAY_SIZE - 1).
 * reference_path is assumed untruncated (e.g. produced by
 * `collharness ref`).
 *
 * Rule: if reference.found_count < 999, the two files' entry sets
 * must match exactly after canonicalization (order p1 < p2 within an
 * entry; sort by (q, qroot, offset, p1, p2)), and found_count must be
 * equal. If reference.found_count >= 999 (the true total is at or
 * past the FOUND_ARRAY_SIZE - 1 cap), found_count must still match
 * exactly between the two files, but only a subset check is applied:
 * every one of observed's stored entries must appear in reference's
 * full set. */
int cmp_run(const char *observed_path, const char *reference_path);

#endif /* COLLHARNESS_CMP_H */
