#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "collcase.h"
#include "gen.h"
#include "cpuref.h"
#include "cmp.h"

static void
usage(const char *prog)
{
	fprintf(stderr,
		"usage: %s <command> [options]\n"
		"\n"
		"commands:\n"
		"  gen     --out FILE [options]      generate a synthetic collcase\n"
		"  ref     --in FILE --out FILE      (re)compute the expected-output\n"
		"                                    section for an existing collcase\n"
		"  cmp     OBSERVED REFERENCE        compare two collcase result files\n"
		"  suite   --out-dir DIR [--large]   generate the named test suite\n"
		"  info    FILE                      print a collcase file's header\n"
		"\n"
		"gen options:\n"
		"  --out FILE            output path (required)\n"
		"  --n N                 element count (default 10000)\n"
		"  --root-bytes 4|8      key width (default 4)\n"
		"  --key-bits K          signed key range +/-2^(K-1) (default 24)\n"
		"  --shift S             bits allocated to p in each value (default 20)\n"
		"  --bucket-hash 0|1     0 = mask hash, 1 = multiplicative (default 0)\n"
		"  --hash-word-cap W     override filter hash-table cap (default: 4096)\n"
		"  --seed S              PRNG seed (default 1)\n"
		"  --density D           filler collision density in [0,1] (default 0.01)\n"
		"  --skew N              force N extra keys into one bucket (default 0)\n"
		"  --no-edge-cases       omit the fixed 0.3 edge-case suite\n"
		"  --no-expected         skip computing the expected-output section\n",
		prog);
}

static int
streq(const char *a, const char *b) { return strcmp(a, b) == 0; }

static int
cmd_gen(int argc, char **argv)
{
	gen_params_t p;
	const char *out = NULL;
	int i;
	collcase_t c;

	gen_params_defaults(&p);

	for (i = 0; i < argc; i++) {
		if (streq(argv[i], "--out") && i + 1 < argc) out = argv[++i];
		else if (streq(argv[i], "--n") && i + 1 < argc) p.n = (uint32_t)strtoul(argv[++i], NULL, 10);
		else if (streq(argv[i], "--root-bytes") && i + 1 < argc) p.root_bytes = (uint32_t)strtoul(argv[++i], NULL, 10);
		else if (streq(argv[i], "--key-bits") && i + 1 < argc) p.key_bits = (uint32_t)strtoul(argv[++i], NULL, 10);
		else if (streq(argv[i], "--shift") && i + 1 < argc) p.shift = (uint32_t)strtoul(argv[++i], NULL, 10);
		else if (streq(argv[i], "--bucket-hash") && i + 1 < argc) p.bucket_hash = (uint32_t)strtoul(argv[++i], NULL, 10);
		else if (streq(argv[i], "--hash-word-cap") && i + 1 < argc) p.hash_word_cap = (uint32_t)strtoul(argv[++i], NULL, 10);
		else if (streq(argv[i], "--seed") && i + 1 < argc) p.seed = (uint64_t)strtoull(argv[++i], NULL, 10);
		else if (streq(argv[i], "--density") && i + 1 < argc) p.collision_density = strtod(argv[++i], NULL);
		else if (streq(argv[i], "--skew") && i + 1 < argc) p.bucket_skew_count = (uint32_t)strtoul(argv[++i], NULL, 10);
		else if (streq(argv[i], "--no-edge-cases")) p.include_edge_cases = 0;
		else if (streq(argv[i], "--no-expected")) p.with_expected = 0;
		else { fprintf(stderr, "gen: unknown option '%s'\n", argv[i]); return 2; }
	}
	if (!out) { fprintf(stderr, "gen: --out is required\n"); return 2; }

	if (gen_generate(&p, &c))
		return 1;

	printf("generated %s: n=%u root_bytes=%u key_bits=%u shift=%u "
			"bucket_hash=%u num_q=%u\n", out, c.n, c.root_bytes,
			c.key_bits, c.shift, c.bucket_hash, c.num_q);
	if (c.has_expected) {
		printf("  expected: found_count=%u candidate_count=%u "
				"dedup_count=%u value_match_count=%u bucket_max=%u\n",
				c.found_count, c.stats.candidate_count,
				c.stats.dedup_count, c.stats.value_match_count,
				c.stats.bucket_max);
	}

	if (collcase_write(out, &c)) {
		collcase_free(&c);
		return 1;
	}
	collcase_free(&c);
	return 0;
}

static int
cmd_ref(int argc, char **argv)
{
	const char *in = NULL, *out = NULL;
	int i;
	collcase_t c;
	cc_found_t *entries = NULL;
	uint32_t entry_count = 0, found_count = 0;

	for (i = 0; i < argc; i++) {
		if (streq(argv[i], "--in") && i + 1 < argc) in = argv[++i];
		else if (streq(argv[i], "--out") && i + 1 < argc) out = argv[++i];
		else { fprintf(stderr, "ref: unknown option '%s'\n", argv[i]); return 2; }
	}
	if (!in || !out) { fprintf(stderr, "ref: --in and --out are required\n"); return 2; }

	if (collcase_read(in, &c))
		return 1;

	if (cpuref_exact(&c, &entries, &entry_count, &found_count)) {
		collcase_free(&c);
		return 1;
	}
	free(c.entries);
	c.entries = entries;
	c.entry_count = entry_count;
	c.found_count = found_count;

	if (cpuref_stats(&c, &c.stats)) {
		collcase_free(&c);
		return 1;
	}
	c.has_expected = 1;

	printf("ref %s: found_count=%u candidate_count=%u dedup_count=%u "
			"value_match_count=%u bucket_max=%u\n",
			in, c.found_count, c.stats.candidate_count,
			c.stats.dedup_count, c.stats.value_match_count,
			c.stats.bucket_max);

	if (collcase_write(out, &c)) {
		collcase_free(&c);
		return 1;
	}
	collcase_free(&c);
	return 0;
}

static int
cmd_cmp(int argc, char **argv)
{
	if (argc != 2) {
		fprintf(stderr, "cmp: usage: cmp OBSERVED REFERENCE\n");
		return 2;
	}
	return cmp_run(argv[0], argv[1]);
}

static int
cmd_info(int argc, char **argv)
{
	collcase_t c;
	if (argc != 1) {
		fprintf(stderr, "info: usage: info FILE\n");
		return 2;
	}
	if (collcase_read(argv[0], &c))
		return 1;
	printf("%s\n", argv[0]);
	printf("  n=%u root_bytes=%u key_bits=%u shift=%u bucket_hash=%u "
			"num_q=%u hash_word_cap=%u\n",
			c.n, c.root_bytes, c.key_bits, c.shift, c.bucket_hash,
			c.num_q, c.hash_word_cap);
	if (c.has_expected) {
		printf("  expected: found_count=%u entry_count=%u "
				"candidate_count=%u dedup_count=%u "
				"value_match_count=%u bucket_max=%u\n",
				c.found_count, c.entry_count,
				c.stats.candidate_count, c.stats.dedup_count,
				c.stats.value_match_count, c.stats.bucket_max);
	} else {
		printf("  (no expected-output section)\n");
	}
	collcase_free(&c);
	return 0;
}

/* ---- suite (plan task 0.7) ------------------------------------------- */

typedef struct {
	const char *name;
	uint32_t n;
	uint32_t root_bytes;
	uint32_t key_bits;
	uint32_t shift;
	uint32_t bucket_hash;
	uint32_t skew;
	int large;
} suite_case_t;

static uint64_t
name_seed(const char *name)
{
	/* FNV-1a -- deterministic, trivial, good enough to turn a suite
	 * case's name into a distinct reproducible seed. */
	uint64_t h = 1469598103934665603ULL;
	const unsigned char *s = (const unsigned char *)name;
	while (*s) {
		h ^= (uint64_t)*s++;
		h *= 1099511628211ULL;
	}
	return h;
}

static int
cmd_suite(int argc, char **argv)
{
	const char *out_dir = NULL;
	int with_large = 0;
	int i;
	static const suite_case_t cases[] = {
		{ "tiny_basic4",      64,      4, 24, 20, 0,   0, 0 },
		{ "tiny_basic8",      64,      8, 40, 20, 0,   0, 0 },
		{ "tiny_hashmode1",   64,      4, 24, 20, 1,   0, 0 },
		{ "tiny_skew_mask",   200,     4, 24, 20, 0, 300, 0 },
		{ "tiny_skew_mix",    200,     4, 24, 20, 1, 300, 0 },
		{ "med_default4",     2000000, 4, 24, 20, 0,   0, 0 },
		{ "med_default8",     2000000, 8, 40, 20, 0,   0, 0 },
		{ "med_hashmode1",    2000000, 4, 24, 20, 1,   0, 0 },
		{ "large_default4",  30000000, 4, 27, 20, 0,   0, 1 },
	};
	size_t ncases = sizeof(cases) / sizeof(cases[0]);
	char path[4096];

	for (i = 0; i < argc; i++) {
		if (streq(argv[i], "--out-dir") && i + 1 < argc) out_dir = argv[++i];
		else if (streq(argv[i], "--large")) with_large = 1;
		else { fprintf(stderr, "suite: unknown option '%s'\n", argv[i]); return 2; }
	}
	if (!out_dir) { fprintf(stderr, "suite: --out-dir is required\n"); return 2; }

	for (i = 0; i < (int)ncases; i++) {
		const suite_case_t *sc = &cases[i];
		gen_params_t p;
		collcase_t c;

		if (sc->large && !with_large) {
			printf("skipping %s (pass --large to include it)\n", sc->name);
			continue;
		}

		gen_params_defaults(&p);
		p.n = sc->n;
		p.root_bytes = sc->root_bytes;
		p.key_bits = sc->key_bits;
		p.shift = sc->shift;
		p.bucket_hash = sc->bucket_hash;
		p.bucket_skew_count = sc->skew;
		p.seed = name_seed(sc->name);
		p.with_expected = 1;
		p.include_edge_cases = 1;

		printf("generating %s (n=%u root_bytes=%u key_bits=%u "
				"bucket_hash=%u skew=%u)...\n", sc->name, sc->n,
				sc->root_bytes, sc->key_bits, sc->bucket_hash,
				sc->skew);

		if (gen_generate(&p, &c))
			return 1;

		printf("  -> n=%u num_q=%u found_count=%u candidate_count=%u "
				"dedup_count=%u value_match_count=%u bucket_max=%u\n",
				c.n, c.num_q, c.found_count,
				c.stats.candidate_count, c.stats.dedup_count,
				c.stats.value_match_count, c.stats.bucket_max);

#if defined(_WIN32)
		snprintf(path, sizeof(path), "%s\\%s.collcase", out_dir, sc->name);
#else
		snprintf(path, sizeof(path), "%s/%s.collcase", out_dir, sc->name);
#endif
		if (collcase_write(path, &c)) {
			collcase_free(&c);
			return 1;
		}
		collcase_free(&c);
	}

	return 0;
}

int
main(int argc, char **argv)
{
	if (argc < 2) {
		usage(argv[0]);
		return 2;
	}
	if (streq(argv[1], "gen"))   return cmd_gen(argc - 2, argv + 2);
	if (streq(argv[1], "ref"))   return cmd_ref(argc - 2, argv + 2);
	if (streq(argv[1], "cmp"))   return cmd_cmp(argc - 2, argv + 2);
	if (streq(argv[1], "info"))  return cmd_info(argc - 2, argv + 2);
	if (streq(argv[1], "suite")) return cmd_suite(argc - 2, argv + 2);
	if (streq(argv[1], "-h") || streq(argv[1], "--help")) {
		usage(argv[0]);
		return 0;
	}
	fprintf(stderr, "unknown command '%s'\n", argv[1]);
	usage(argv[0]);
	return 2;
}
