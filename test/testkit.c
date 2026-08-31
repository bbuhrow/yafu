/*----------------------------------------------------------------------
 YAFU modular test system -- framework implementation ("testkit").
 See testkit.h for the design and public API. Public domain.
----------------------------------------------------------------------*/
#include "testkit.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <setjmp.h>

#if defined(_WIN32)
  #define WIN32_LEAN_AND_MEAN
  #include <windows.h>
#else
  #include <time.h>
#endif

/* ================================================================== *
 * Timer
 * ================================================================== */
double tk_now_sec(void)
{
#if defined(_WIN32)
    static LARGE_INTEGER freq;
    static int have_freq = 0;
    LARGE_INTEGER now;
    if (!have_freq) { QueryPerformanceFrequency(&freq); have_freq = 1; }
    QueryPerformanceCounter(&now);
    return (double)now.QuadPart / (double)freq.QuadPart;
#elif defined(CLOCK_MONOTONIC)
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
#else
    return (double)clock() / (double)CLOCKS_PER_SEC;
#endif
}

/* ================================================================== *
 * Compiler id
 * ================================================================== */
const char *tk_compiler_id(void)
{
    static char buf[64];
#if defined(__INTEL_LLVM_COMPILER)
    snprintf(buf, sizeof buf, "icx %d", __INTEL_LLVM_COMPILER);
#elif defined(__INTEL_COMPILER)
    snprintf(buf, sizeof buf, "icc %d", __INTEL_COMPILER);
#elif defined(__clang__)
    snprintf(buf, sizeof buf, "clang %d.%d.%d",
             __clang_major__, __clang_minor__, __clang_patchlevel__);
#elif defined(_MSC_VER)
    snprintf(buf, sizeof buf, "msvc %d", _MSC_VER);
#elif defined(__GNUC__)
    snprintf(buf, sizeof buf, "gcc %d.%d.%d",
             __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#else
    snprintf(buf, sizeof buf, "unknown");
#endif
    return buf;
}

/* ================================================================== *
 * RNG -- splitmix64. Small, fast, good enough for test data; the point
 * is determinism, not cryptographic quality.
 * ================================================================== */
void tk_rng_seed(tk_rng *r, uint64_t seed) { r->s = seed; }

uint64_t tk_rng_u64(tk_rng *r)
{
    uint64_t z = (r->s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

uint64_t tk_rng_range(tk_rng *r, uint64_t n)
{
    /* Lemire-style unbiased reduction. n > 0 assumed. */
    uint64_t x, m_hi, m_lo, t;
    if (n == 0) return 0;
    x = tk_rng_u64(r);
#if defined(__SIZEOF_INT128__)
    {
        __uint128_t m = (__uint128_t)x * (__uint128_t)n;
        m_lo = (uint64_t)m;
        m_hi = (uint64_t)(m >> 64);
    }
#else
    /* 64x64->128 without __int128 */
    {
        uint64_t xl = x & 0xFFFFFFFFULL, xh = x >> 32;
        uint64_t nl = n & 0xFFFFFFFFULL, nh = n >> 32;
        uint64_t ll = xl * nl, lh = xl * nh, hl = xh * nl, hh = xh * nh;
        uint64_t cross = (ll >> 32) + (lh & 0xFFFFFFFFULL) + (hl & 0xFFFFFFFFULL);
        m_lo = (ll & 0xFFFFFFFFULL) | (cross << 32);
        m_hi = hh + (lh >> 32) + (hl >> 32) + (cross >> 32);
    }
#endif
    if (m_lo < n) {
        t = (uint64_t)(-(int64_t)n) % n;       /* 2^64 mod n */
        while (m_lo < t) {
            x = tk_rng_u64(r);
#if defined(__SIZEOF_INT128__)
            {
                __uint128_t m = (__uint128_t)x * (__uint128_t)n;
                m_lo = (uint64_t)m; m_hi = (uint64_t)(m >> 64);
            }
#else
            {
                uint64_t xl = x & 0xFFFFFFFFULL, xh = x >> 32;
                uint64_t nl = n & 0xFFFFFFFFULL, nh = n >> 32;
                uint64_t ll = xl*nl, lh = xl*nh, hl = xh*nl, hh = xh*nh;
                uint64_t cross = (ll >> 32) + (lh & 0xFFFFFFFFULL) + (hl & 0xFFFFFFFFULL);
                m_lo = (ll & 0xFFFFFFFFULL) | (cross << 32);
                m_hi = hh + (lh >> 32) + (hl >> 32) + (cross >> 32);
            }
#endif
        }
    }
    return m_hi;
}

uint64_t tk_rng_nbits(tk_rng *r, int bits)
{
    uint64_t x;
    if (bits <= 0)  return 0;
    if (bits >= 64) return tk_rng_u64(r) | (1ULL << 63);
    x = tk_rng_u64(r) & ((1ULL << bits) - 1);
    x |= (1ULL << (bits - 1));       /* force exact bit-length */
    return x;
}

/* ================================================================== *
 * Context
 * ================================================================== */
#define TK_MAX_BENCH 16

typedef struct { char label[48]; double ns; double speedup; } tk_bench_row;

struct tk_ctx {
    const char *module;
    const char *test;
    long        checks;
    long        fails;
    int         verbose;
    int         do_bench;
    tk_rng      rng;
    jmp_buf     abort_env;
    int         have_abort_env;
    tk_bench_row bench[TK_MAX_BENCH];
    int         nbench;
};

tk_rng *tk_rng_of(tk_ctx *tk)        { return &tk->rng; }
int     tk_bench_enabled(tk_ctx *tk) { return tk->do_bench; }
int     tk_verbose(tk_ctx *tk)       { return tk->verbose; }

/* ================================================================== *
 * Assertions
 * ================================================================== */
static void tk__fail_banner(tk_ctx *tk, const char *file, int line)
{
    fprintf(stderr, "    FAIL %s::%s  (%s:%d)\n",
            tk->module, tk->test, file, line);
}

int tk_check(tk_ctx *tk, int cond, const char *expr, const char *file, int line)
{
    tk->checks++;
    if (!cond) {
        tk->fails++;
        tk__fail_banner(tk, file, line);
        fprintf(stderr, "        expected true: %s\n", expr);
    }
    return cond ? 1 : 0;
}

int tk_checkf(tk_ctx *tk, int cond, const char *file, int line, const char *fmt, ...)
{
    tk->checks++;
    if (!cond) {
        va_list ap;
        tk->fails++;
        tk__fail_banner(tk, file, line);
        fprintf(stderr, "        ");
        va_start(ap, fmt);
        vfprintf(stderr, fmt, ap);
        va_end(ap);
        fprintf(stderr, "\n");
    }
    return cond ? 1 : 0;
}

int tk_eq_u64(tk_ctx *tk, uint64_t got, uint64_t exp,
              const char *expr, const char *file, int line)
{
    tk->checks++;
    if (got != exp) {
        tk->fails++;
        tk__fail_banner(tk, file, line);
        fprintf(stderr, "        %s\n", expr);
        fprintf(stderr, "        got 0x%016llx (%llu), expected 0x%016llx (%llu)\n",
                (unsigned long long)got, (unsigned long long)got,
                (unsigned long long)exp, (unsigned long long)exp);
        return 0;
    }
    return 1;
}

void tk_abort(tk_ctx *tk, const char *file, int line, const char *fmt, ...)
{
    va_list ap;
    tk->checks++;
    tk->fails++;
    tk__fail_banner(tk, file, line);
    fprintf(stderr, "        ABORT: ");
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fprintf(stderr, "\n");
    if (tk->have_abort_env) longjmp(tk->abort_env, 1);
    /* no jmp context (shouldn't happen inside the runner): give up hard */
    fprintf(stderr, "        (no abort context; exiting)\n");
    exit(2);
}

/* ================================================================== *
 * Benchmark reporting
 * ================================================================== */
volatile uint64_t tk_sink_u64 = 0;

void tk_report_bench(tk_ctx *tk, const char *label, double ns_per_op, double speedup)
{
    tk_bench_row *row;
    if (tk->nbench >= TK_MAX_BENCH) return;
    row = &tk->bench[tk->nbench++];
    strncpy(row->label, label, sizeof row->label - 1);
    row->label[sizeof row->label - 1] = '\0';
    row->ns = ns_per_op;
    row->speedup = speedup;
}

/* ================================================================== *
 * Runner
 * ================================================================== */
typedef struct {
    const char **names;   int nnames;   /* --module filters (NULL = all)  */
    const char **tags;    int ntags;    /* --tag filters    (NULL = all)  */
    uint64_t seed;
    int do_bench;
    int verbose;
    int stop_on_fail;     /* stop launching further tests after 1st fail  */
    int list_only;
} tk_config;

static int tk__name_selected(const tk_config *cfg, const char *name)
{
    int i;
    if (cfg->nnames == 0) return 1;
    for (i = 0; i < cfg->nnames; i++)
        if (strcmp(cfg->names[i], name) == 0) return 1;
    return 0;
}

/* a tag matches if it appears as a whitespace-delimited token in `tags` */
static int tk__has_tag(const char *tags, const char *want)
{
    const char *p = tags;
    size_t wl = strlen(want);
    if (!tags) return 0;
    while (*p) {
        const char *start;
        size_t len;
        while (*p == ' ' || *p == '\t') p++;
        start = p;
        while (*p && *p != ' ' && *p != '\t') p++;
        len = (size_t)(p - start);
        if (len == wl && strncmp(start, want, wl) == 0) return 1;
    }
    return 0;
}

static int tk__tag_selected(const tk_config *cfg, const char *tags)
{
    int i;
    if (cfg->ntags == 0) return 1;
    for (i = 0; i < cfg->ntags; i++)
        if (tk__has_tag(tags, cfg->tags[i])) return 1;
    return 0;
}

/* Run one test with an abort context. Keeping setjmp in its own function
 * means the runner's loop variables are not in the setjmp scope, so they
 * cannot be clobbered by a longjmp out of a TK_REQUIRE. */
static void tk__run_one(tk_ctx *tk, const tk_test *t)
{
    tk->have_abort_env = 1;
    if (setjmp(tk->abort_env) == 0)
        t->fn(tk);
    /* else: test aborted via tk_abort(); counters already updated */
}

static void tk__usage(const char *argv0)
{
    printf(
"Usage: %s [options] [module ...]\n"
"  --module NAME   run only this module (repeatable; bare args also work)\n"
"  --tag TAG       run only tests carrying TAG (repeatable)\n"
"  --seed N        RNG seed (default 1); printed so failures reproduce\n"
"  --bench         run micro-benchmarks (off by default)\n"
"  --stop          stop after the first failing test\n"
"  --list          list modules and tests, then exit\n"
"  -v, --verbose   extra output\n"
"  -h, --help      this help\n",
        argv0);
}

static void tk__list(const tk_module *const *m, int n)
{
    int i, j;
    for (i = 0; i < n; i++) {
        printf("%-14s %s\n", m[i]->name, m[i]->desc ? m[i]->desc : "");
        for (j = 0; j < m[i]->ntests; j++)
            printf("    %-24s [%s]\n", m[i]->tests[j].name,
                   m[i]->tests[j].tags ? m[i]->tests[j].tags : "");
    }
}

int tk_run(const tk_module *const *modules, int nmodules, int argc, char **argv)
{
    tk_config cfg;
    const char *name_buf[64];
    const char *tag_buf[64];
    int i, j;
    long total_checks = 0, total_fails = 0, total_tests = 0, total_run = 0;
    int  any_module_ran = 0;

    memset(&cfg, 0, sizeof cfg);
    cfg.names = name_buf; cfg.tags = tag_buf;
    cfg.seed = 1;

    for (i = 1; i < argc; i++) {
        const char *a = argv[i];
        if      (!strcmp(a, "-h") || !strcmp(a, "--help")) { tk__usage(argv[0]); return 0; }
        else if (!strcmp(a, "--list"))    cfg.list_only = 1;
        else if (!strcmp(a, "--bench"))   cfg.do_bench = 1;
        else if (!strcmp(a, "--stop"))    cfg.stop_on_fail = 1;
        else if (!strcmp(a, "-v") || !strcmp(a, "--verbose")) cfg.verbose = 1;
        else if (!strcmp(a, "--seed") && i + 1 < argc)
            cfg.seed = strtoull(argv[++i], NULL, 0);
        else if (!strcmp(a, "--module") && i + 1 < argc) {
            if (cfg.nnames < 64) cfg.names[cfg.nnames++] = argv[++i];
        }
        else if (!strcmp(a, "--tag") && i + 1 < argc) {
            if (cfg.ntags < 64) cfg.tags[cfg.ntags++] = argv[++i];
        }
        else if (a[0] == '-') {
            fprintf(stderr, "unknown option: %s\n", a);
            tk__usage(argv[0]); return 2;
        }
        else { /* bare positional = module name */
            if (cfg.nnames < 64) cfg.names[cfg.nnames++] = a;
        }
    }

    if (cfg.list_only) { tk__list(modules, nmodules); return 0; }

    printf("=========================================================\n");
    printf(" YAFU testkit %s   compiler: %s   seed: %llu\n",
           TK_VERSION, tk_compiler_id(), (unsigned long long)cfg.seed);
    printf("=========================================================\n");

    for (i = 0; i < nmodules; i++) {
        const tk_module *mod = modules[i];
        int module_header_printed = 0;

        if (!tk__name_selected(&cfg, mod->name)) continue;

        for (j = 0; j < mod->ntests; j++) {
            const tk_test *t = &mod->tests[j];
            tk_ctx tk;
            int k;
            double tk__elapsed = 0.0;

            total_tests++;
            if (!tk__tag_selected(&cfg, t->tags)) continue;

            if (!module_header_printed) {
                printf("\n--- %s : %s ---\n", mod->name, mod->desc ? mod->desc : "");
                module_header_printed = 1;
                any_module_ran = 1;
            }

            memset(&tk, 0, sizeof tk);
            tk.module = mod->name;
            tk.test   = t->name;
            tk.verbose = cfg.verbose;
            tk.do_bench = cfg.do_bench;
            /* independent, reproducible stream per test */
            tk_rng_seed(&tk.rng, cfg.seed + 0x1000ULL * (uint64_t)(i + 1)
                                          + (uint64_t)(j + 1));

            {
                double tk__t0 = tk_now_sec();
                tk__run_one(&tk, t);
                tk__elapsed = tk_now_sec() - tk__t0;
            }

            total_run++;
            total_checks += tk.checks;
            total_fails  += tk.fails;

            printf("  %-24s %6ld checks  %-4s  %9.3fs\n",
                   t->name, tk.checks,
                   tk.fails == 0 ? "ok" : "FAIL", tk__elapsed);
            for (k = 0; k < tk.nbench; k++) {
                if (tk.bench[k].speedup > 0.0)
                    printf("        bench %-28s %10.3f ns/op  %6.2fx\n",
                           tk.bench[k].label, tk.bench[k].ns, tk.bench[k].speedup);
                else
                    printf("        bench %-28s %10.3f ns/op\n",
                           tk.bench[k].label, tk.bench[k].ns);
            }
            fflush(stdout);   /* live progress for long-running tests */

            if (cfg.stop_on_fail && tk.fails > 0) {
                printf("\n[--stop] halting after first failing test\n");
                goto done;
            }
        }
    }

done:
    printf("\n=========================================================\n");
    printf(" tests run: %ld/%ld   checks: %ld   failures: %ld   -> %s\n",
           total_run, total_tests, total_checks, total_fails,
           total_fails == 0 ? "PASS" : "FAIL");
    printf("=========================================================\n");

    if (!any_module_ran)
        fprintf(stderr, "warning: no tests matched the given filters\n");

    return total_fails == 0 ? 0 : 1;
}
