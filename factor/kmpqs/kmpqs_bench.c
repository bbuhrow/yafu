/*
kmpqs_bench: benchmark and validation driver for factor/kmpqs.

Part of the kmpqs module (GPL, see kmpqs.h).

Usage:
  kmpqs_bench [options] file      numbers, one per line (text after the
                                  first ',' or space is ignored)
  kmpqs_bench [options] -gen B C  C random semiprimes of B bits (two primes
                                  of B/2 bits), fixed seed

Options:
  -1 | -3 | -a   use kmpqs1_factor (default), kmpqs3_factor, or the
                 size-dispatching kmpqs_factor
  -maxbits M     max_bits argument (default 80 for -1, 200 otherwise; these
                 match lasieve's mpqstest / mpqs3test)
  -seed S        seed for -gen (default 42)
  -dump          write the generated numbers to stdout and exit
  -v             kmpqs_verbose = 1
  -threads T     also factor the inputs with T threads (thread-local
                 contexts) and compare with a single-threaded run

When built with -DMPQS_STAT -DMPQS3_STAT the "Stat:" lines match the output
of lasieve's mpqstest / mpqs3test, so the two can be diffed.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <time.h>
#endif
#include "gmp.h"
#include "kmpqs.h"
#if !defined(MPQS_STAT) && !defined(MPQS3_STAT)
#include <pthread.h>
#endif

#if defined(MPQS_STAT) || defined(MPQS3_STAT)
uint64_t stat_asm_eval = 0, stat_asm_td = 0;
uint64_t stat_td_cand = 0, stat_td_surv = 0;
uint64_t stat_ff = 0, stat_pf = 0, stat_comb = 0, stat_pp = 0;
uint32_t stat_asm_div = 0, stat_final_mulmod = 0;
uint32_t stat_counter0 = 0, stat_counter1 = 0, stat_retry = 0;
uint32_t stat_size[15] = { 0 };
#endif

static double now_sec(void)
{
#ifdef _WIN32
	LARGE_INTEGER f, c;

	QueryPerformanceFrequency(&f);
	QueryPerformanceCounter(&c);
	return (double)c.QuadPart / (double)f.QuadPart;
#else
	struct timespec ts;

	clock_gettime(CLOCK_MONOTONIC, &ts);
	return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
#endif
}

/* Product of the factors equals N and each factor is a probable prime. */
static int check_factors(mpz_t N, mpz_t *f, int n)
{
	mpz_t prod;
	int i, ok = 1;

	mpz_init_set_ui(prod, 1);
	for (i = 0; i < n; i++) {
		mpz_mul(prod, prod, f[i]);
		if (!mpz_probab_prime_p(f[i], 10))
			ok = 0;
	}
	if (mpz_cmp(prod, N) != 0)
		ok = 0;
	mpz_clear(prod);
	return ok;
}

typedef struct {
	mpz_t *inputs;
	int count, variant;
	size_t maxbits;
	int *nf;		/* per input: return value */
	uint64_t *fhash;	/* per input: sum of the factors mod 2^64 */
} jobs_t;

typedef struct {
	jobs_t *jobs;
	int start, stride;
	int suc, fail, wrong;
} worker_arg_t;

/* factor inputs start, start+stride, ... with the calling thread's context */
static void bench_worker(worker_arg_t *w)
{
	jobs_t *jb = w->jobs;
	mpz_t N, *f;
	int i, j, nf;

	mpz_init(N);
	for (i = w->start; i < jb->count; i += w->stride) {
		mpz_set(N, jb->inputs[i]);
		if (jb->variant == 1) nf = kmpqs1_factor(N, jb->maxbits, &f);
		else if (jb->variant == 3) nf = kmpqs3_factor(N, jb->maxbits, &f);
		else nf = kmpqs_factor(N, jb->maxbits, &f);
		jb->nf[i] = nf;
		jb->fhash[i] = 0;
		if (nf >= 1) {
			w->suc++;
			for (j = 0; j < nf; j++)
				jb->fhash[i] += (uint64_t)mpz_getlimbn(f[j], 0);
			if (!check_factors(jb->inputs[i], f, nf)) {
				w->wrong++;
				gmp_fprintf(stderr, "wrong factorization of %Zd\n", jb->inputs[i]);
			}
		}
		else {
			w->fail++;
			if (kmpqs_verbose)
				gmp_fprintf(stderr, "failed (%d): %Zd\n", nf, jb->inputs[i]);
		}
	}
	mpz_clear(N);
}

#if !defined(MPQS_STAT) && !defined(MPQS3_STAT)
static void *bench_worker_thread(void *arg)
{
	bench_worker((worker_arg_t *)arg);
	kmpqs_thread_cleanup();
	return NULL;
}
#endif

int main(int argc, char **argv)
{
	int nthreads = 1;
	jobs_t jobs;
	int variant = 1, i, dump = 0;
	int suc = 0, fail = 0, wrong = 0, count = 0;
	long maxbits = -1;
	unsigned long seed = 42;
	int gen_bits = 0, gen_count = 0;
	const char *fname = NULL;
	mpz_t N, p, q;
	mpz_t *inputs = NULL;
	gmp_randstate_t rs;
	double t0, t1;

	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-1") == 0) variant = 1;
		else if (strcmp(argv[i], "-3") == 0) variant = 3;
		else if (strcmp(argv[i], "-a") == 0) variant = 0;
		else if (strcmp(argv[i], "-v") == 0) kmpqs_verbose = 1;
		else if (strcmp(argv[i], "-dump") == 0) dump = 1;
		else if (strcmp(argv[i], "-maxbits") == 0 && i + 1 < argc) maxbits = atol(argv[++i]);
		else if (strcmp(argv[i], "-threads") == 0 && i + 1 < argc) nthreads = atoi(argv[++i]);
		else if (strcmp(argv[i], "-seed") == 0 && i + 1 < argc) seed = strtoul(argv[++i], NULL, 10);
		else if (strcmp(argv[i], "-gen") == 0 && i + 2 < argc) {
			gen_bits = atoi(argv[++i]);
			gen_count = atoi(argv[++i]);
		}
		else fname = argv[i];
	}
	if (maxbits < 0)
		maxbits = (variant == 1) ? 80 : 200;

	mpz_init(N);
	mpz_init(p);
	mpz_init(q);

	if (gen_bits > 0) {
		gmp_randinit_default(rs);
		gmp_randseed_ui(rs, seed);
		inputs = (mpz_t *)malloc(gen_count * sizeof(mpz_t));
		for (i = 0; i < gen_count; i++) {
			int b1 = gen_bits / 2, b2 = gen_bits - gen_bits / 2;

			mpz_init(inputs[i]);
			do {
				mpz_urandomb(p, rs, b1);
				mpz_setbit(p, b1 - 1);
				mpz_nextprime(p, p);
				mpz_urandomb(q, rs, b2);
				mpz_setbit(q, b2 - 1);
				mpz_nextprime(q, q);
				mpz_mul(inputs[i], p, q);
			} while (mpz_sizeinbase(inputs[i], 2) != (size_t)gen_bits ||
				mpz_cmp(p, q) == 0);
		}
		count = gen_count;
		gmp_randclear(rs);
	}
	else {
		FILE *fi;
		char line[1024];
		int alloc = 1024;

		if (fname == NULL) {
			fprintf(stderr, "usage: kmpqs_bench [-1|-3|-a] [-maxbits M] [-v] (file | -gen bits count)\n");
			return 1;
		}
		fi = (strcmp(fname, "-") == 0) ? stdin : fopen(fname, "r");
		if (fi == NULL) {
			fprintf(stderr, "cannot open %s\n", fname);
			return 1;
		}
		inputs = (mpz_t *)malloc(alloc * sizeof(mpz_t));
		while (fgets(line, sizeof(line), fi) != NULL) {
			line[strcspn(line, ", \t\r\n")] = '\0';
			if (line[0] == '\0')
				continue;
			if (count == alloc) {
				alloc *= 2;
				inputs = (mpz_t *)realloc(inputs, alloc * sizeof(mpz_t));
			}
			mpz_init(inputs[count]);
			if (mpz_set_str(inputs[count], line, 10)) {
				mpz_clear(inputs[count]);
				break;
			}
			count++;
		}
		if (fi != stdin)
			fclose(fi);
	}

	if (dump) {
		for (i = 0; i < count; i++)
			gmp_printf("%Zd\n", inputs[i]);
		return 0;
	}

	jobs.inputs = inputs;
	jobs.count = count;
	jobs.variant = variant;
	jobs.maxbits = (size_t)maxbits;
	jobs.nf = (int *)calloc(count, sizeof(int));
	jobs.fhash = (uint64_t *)calloc(count, sizeof(uint64_t));

	if (nthreads > 1) {
#if defined(MPQS_STAT) || defined(MPQS3_STAT)
		fprintf(stderr, "-threads needs a build without MPQS_STAT (the counters are global)\n");
		return 1;
#else
		/* reference results, single-threaded */
		int *nf1 = (int *)calloc(count, sizeof(int));
		uint64_t *h1 = (uint64_t *)calloc(count, sizeof(uint64_t));
		worker_arg_t w0 = { &jobs, 0, 1, 0, 0, 0 };
		worker_arg_t *wa = (worker_arg_t *)calloc(nthreads, sizeof(worker_arg_t));
		pthread_t *tid = (pthread_t *)calloc(nthreads, sizeof(pthread_t));
		int mismatch = 0;

		bench_worker(&w0);
		memcpy(nf1, jobs.nf, count * sizeof(int));
		memcpy(h1, jobs.fhash, count * sizeof(uint64_t));

		t0 = now_sec();
		for (i = 0; i < nthreads; i++) {
			wa[i].jobs = &jobs;
			wa[i].start = i;
			wa[i].stride = nthreads;
			pthread_create(&tid[i], NULL, bench_worker_thread, &wa[i]);
		}
		for (i = 0; i < nthreads; i++) {
			pthread_join(tid[i], NULL);
			suc += wa[i].suc;
			fail += wa[i].fail;
			wrong += wa[i].wrong;
		}
		t1 = now_sec();
		for (i = 0; i < count; i++)
			if (nf1[i] != jobs.nf[i] || h1[i] != jobs.fhash[i])
				mismatch++;
		printf("threads: %d, results differing from the single-threaded run: %d\n",
			nthreads, mismatch);
		if (mismatch)
			wrong += mismatch;
#endif
	}
	else {
		worker_arg_t w = { &jobs, 0, 1, 0, 0, 0 };

		t0 = now_sec();
		bench_worker(&w);
		t1 = now_sec();
		suc = w.suc;
		fail = w.fail;
		wrong = w.wrong;
	}

	printf("success: %d  failure: %d\n", suc, fail);
#if defined(MPQS_STAT) || defined(MPQS3_STAT)
	if (stat_retry) printf("Warning: %u retries\n", stat_retry);
	printf("Stat: sieves %llu, td %llu->%llu->%llu\n",
		(unsigned long long)stat_asm_eval, (unsigned long long)stat_td_cand,
		(unsigned long long)stat_asm_td, (unsigned long long)stat_td_surv);
	printf("  ff %llu, pf %llu->%llu\n", (unsigned long long)stat_ff,
		(unsigned long long)stat_pf, (unsigned long long)stat_comb);
#endif
	printf("wrong: %d\n", wrong);
	printf("time: %.3f s (%.2f us per input)\n", t1 - t0,
		count ? 1e6 * (t1 - t0) / count : 0.0);

	return wrong ? 2 : 0;
}
