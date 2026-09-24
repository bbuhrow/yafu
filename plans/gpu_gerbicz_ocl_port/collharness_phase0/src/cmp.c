#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "collcase.h"
#include "cmp.h"

#define SATURATION_THRESHOLD 999u

static void
canon_entry(cc_found_t *e)
{
	if (e->p1 > e->p2) {
		uint32_t t = e->p1;
		e->p1 = e->p2;
		e->p2 = t;
	}
}

static int
cmp_found(const void *a, const void *b)
{
	const cc_found_t *x = (const cc_found_t *)a;
	const cc_found_t *y = (const cc_found_t *)b;
	if (x->q != y->q) return (x->q < y->q) ? -1 : 1;
	if (x->qroot != y->qroot) return (x->qroot < y->qroot) ? -1 : 1;
	if (x->offset != y->offset) return (x->offset < y->offset) ? -1 : 1;
	if (x->p1 != y->p1) return (x->p1 < y->p1) ? -1 : 1;
	if (x->p2 != y->p2) return (x->p2 < y->p2) ? -1 : 1;
	return 0;
}

static cc_found_t *
canonicalize(const cc_found_t *src, uint32_t n)
{
	cc_found_t *out = (cc_found_t *)malloc((size_t)n * sizeof(cc_found_t));
	uint32_t i;
	if (n > 0 && !out)
		return NULL;
	for (i = 0; i < n; i++) {
		out[i] = src[i];
		canon_entry(&out[i]);
	}
	qsort(out, n, sizeof(cc_found_t), cmp_found);
	return out;
}

/* Binary search for `key` in a sorted canonical array. */
static int
canon_contains(const cc_found_t *arr, uint32_t n, const cc_found_t *key)
{
	uint32_t lo = 0, hi = n;
	while (lo < hi) {
		uint32_t mid = lo + (hi - lo) / 2;
		int c = cmp_found(&arr[mid], key);
		if (c == 0)
			return 1;
		if (c < 0)
			lo = mid + 1;
		else
			hi = mid;
	}
	return 0;
}

int
cmp_run(const char *observed_path, const char *reference_path)
{
	collcase_t obs, ref;
	cc_found_t *obs_c = NULL, *ref_c = NULL;
	int rc = 0;

	collcase_init(&obs);
	collcase_init(&ref);

	if (collcase_read(observed_path, &obs)) {
		fprintf(stderr, "cmp: failed to read %s\n", observed_path);
		return 2;
	}
	if (collcase_read(reference_path, &ref)) {
		fprintf(stderr, "cmp: failed to read %s\n", reference_path);
		collcase_free(&obs);
		return 2;
	}
	if (!obs.has_expected || !ref.has_expected) {
		fprintf(stderr, "cmp: both files must carry an expected-output "
				"section\n");
		collcase_free(&obs);
		collcase_free(&ref);
		return 2;
	}

	printf("observed: found_count=%u entry_count=%u (%s)\n",
			obs.found_count, obs.entry_count, observed_path);
	printf("reference: found_count=%u entry_count=%u (%s)\n",
			ref.found_count, ref.entry_count, reference_path);

	if (obs.found_count != ref.found_count) {
		printf("MISMATCH: found_count differs (%u vs %u)\n",
				obs.found_count, ref.found_count);
		rc = 1;
		goto done;
	}

	obs_c = canonicalize(obs.entries, obs.entry_count);
	ref_c = canonicalize(ref.entries, ref.entry_count);
	if ((obs.entry_count && !obs_c) || (ref.entry_count && !ref_c)) {
		fprintf(stderr, "cmp: out of memory\n");
		rc = 2;
		goto done;
	}

	if (ref.found_count < SATURATION_THRESHOLD) {
		uint32_t i;
		int ok = (obs.entry_count == ref.entry_count);
		if (ok) {
			for (i = 0; i < ref.entry_count; i++) {
				if (cmp_found(&obs_c[i], &ref_c[i]) != 0) {
					ok = 0;
					break;
				}
			}
		}
		if (ok) {
			printf("MATCH: exact entry sets agree (%u entries)\n",
					ref.entry_count);
		} else {
			printf("MISMATCH: entry sets differ (unsaturated -- "
					"exact match required)\n");
			rc = 1;
		}
	} else {
		uint32_t i, missing = 0;
		printf("(saturated at/above %u -- subset check only)\n",
				SATURATION_THRESHOLD);
		for (i = 0; i < obs.entry_count; i++) {
			if (!canon_contains(ref_c, ref.entry_count, &obs_c[i]))
				missing++;
		}
		if (missing == 0) {
			printf("MATCH: all %u observed entries are present in "
					"the reference set\n", obs.entry_count);
		} else {
			printf("MISMATCH: %u of %u observed entries are NOT in "
					"the reference set\n", missing,
					obs.entry_count);
			rc = 1;
		}
	}

done:
	free(obs_c);
	free(ref_c);
	collcase_free(&obs);
	collcase_free(&ref);
	return rc;
}
