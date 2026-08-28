/*----------------------------------------------------------------------
 YAFU modular test system -- framework header ("testkit").

 A tiny, dependency-free unit-test + micro-benchmark harness tailored to
 YAFU's portability matrix (gcc/clang/icc/msvc on linux/wsl/msys2/windows).

 Design notes:
   * Explicit registration -- no linker-section auto-registration (that is
     an MSVC portability minefield). Each test module exposes a single
     `const tk_module` descriptor; test_main.c lists the ones to run.
   * Assertions are non-aborting by default (TK_CHECK*): a failing check is
     recorded and the test continues, so one test reports every failure it
     finds in a single run. TK_REQUIRE aborts just the current test (via
     setjmp/longjmp in the runner) when continuing makes no sense.
   * Benchmarking (TK_BENCH) is opt-in (--bench) and never affects pass/fail.
   * Deterministic, seedable RNG so failures reproduce from the printed seed.
   * Exit status is 0 iff every selected check passed -- CI-wireable.

 Public domain, following the rest of YAFU.
----------------------------------------------------------------------*/
#ifndef YAFU_TESTKIT_H
#define YAFU_TESTKIT_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define TK_VERSION "0.1"

/* ------------------------------------------------------------------ *
 * Deterministic RNG (splitmix64). Seeded once per run; each test gets
 * an independent, reproducible stream derived from (seed, test index).
 * ------------------------------------------------------------------ */
typedef struct { uint64_t s; } tk_rng;

void     tk_rng_seed(tk_rng *r, uint64_t seed);
uint64_t tk_rng_u64(tk_rng *r);
uint64_t tk_rng_range(tk_rng *r, uint64_t n);   /* uniform in [0, n), n > 0 */
uint64_t tk_rng_nbits(tk_rng *r, int bits);     /* exact bit-length (top bit set), 1..64 */

/* ------------------------------------------------------------------ *
 * Per-test context. Opaque: tests touch it only through the API below.
 * ------------------------------------------------------------------ */
typedef struct tk_ctx tk_ctx;

/* Accessors a test may need. */
tk_rng *tk_rng_of(tk_ctx *tk);        /* this test's private RNG stream    */
int     tk_bench_enabled(tk_ctx *tk); /* nonzero if --bench was given      */
int     tk_verbose(tk_ctx *tk);       /* nonzero if -v/--verbose was given */

/* Monotonic wall-clock seconds (QPC / clock_gettime / clock fallback). */
double  tk_now_sec(void);

/* Short identifier for the building compiler, e.g. "gcc 13.3.0". */
const char *tk_compiler_id(void);

/* ------------------------------------------------------------------ *
 * Assertions. All return 1 on pass, 0 on fail (so a test can branch on
 * the result). Failures are recorded against the current test.
 * ------------------------------------------------------------------ */
int tk_check (tk_ctx *tk, int cond, const char *expr, const char *file, int line);
int tk_checkf(tk_ctx *tk, int cond, const char *file, int line, const char *fmt, ...);
int tk_eq_u64(tk_ctx *tk, uint64_t got, uint64_t exp,
              const char *expr, const char *file, int line);
/* Abort the current test immediately (does not return). */
void tk_abort(tk_ctx *tk, const char *file, int line, const char *fmt, ...);

#define TK_CHECK(tk, cond)        tk_check((tk), (cond) ? 1 : 0, #cond, __FILE__, __LINE__)
#define TK_CHECKF(tk, cond, ...)  tk_checkf((tk), (cond) ? 1 : 0, __FILE__, __LINE__, __VA_ARGS__)
#define TK_EQ_U64(tk, got, exp)   tk_eq_u64((tk), (uint64_t)(got), (uint64_t)(exp), \
                                            #got " == " #exp, __FILE__, __LINE__)
#define TK_REQUIRE(tk, cond, ...) \
    do { if (!(cond)) tk_abort((tk), __FILE__, __LINE__, __VA_ARGS__); } while (0)

/* ------------------------------------------------------------------ *
 * Benchmarking. TK_BENCH runs BODY, auto-calibrating the iteration count
 * to at least TK_BENCH_MINSEC, then reports the best of TK_BENCH_REPS.
 *   ops_per_iter : number of measured operations one BODY performs
 *                  (result is normalised to ns/op).
 *   speedup      : ratio vs a reference for the report, or 0.0 if none.
 * Only runs when benchmarking is enabled; never affects pass/fail.
 * ------------------------------------------------------------------ */
void tk_report_bench(tk_ctx *tk, const char *label, double ns_per_op, double speedup);

/* volatile sink to defeat dead-code elimination in benchmark bodies. */
extern volatile uint64_t tk_sink_u64;

#define TK_BENCH_MINSEC 0.10
#define TK_BENCH_REPS   5
#define TK_BENCH(tk, label, ops_per_iter, speedup, BODY)                        \
    do {                                                                        \
        if (tk_bench_enabled(tk)) {                                             \
            long   tk__iters = 1;                                               \
            double tk__best  = 1e30;                                            \
            long   tk__i;                                                       \
            int    tk__r;                                                       \
            for (;;) {                                                          \
                double tk__t0 = tk_now_sec();                                   \
                for (tk__i = 0; tk__i < tk__iters; ++tk__i) { BODY; }           \
                double tk__el = tk_now_sec() - tk__t0;                          \
                if (tk__el >= TK_BENCH_MINSEC) break;                           \
                if (tk__el <= 0.0) tk__iters *= 4;                              \
                else {                                                          \
                    double tk__f = (TK_BENCH_MINSEC / tk__el) * 1.3;            \
                    if (tk__f < 2.0) tk__f = 2.0;                               \
                    tk__iters = (long)((double)tk__iters * tk__f);             \
                }                                                               \
                if (tk__iters <= 0) tk__iters = (1L << 20);                     \
            }                                                                   \
            for (tk__r = 0; tk__r < TK_BENCH_REPS; ++tk__r) {                   \
                double tk__t0 = tk_now_sec();                                   \
                for (tk__i = 0; tk__i < tk__iters; ++tk__i) { BODY; }           \
                double tk__el = tk_now_sec() - tk__t0;                          \
                double tk__ns = (tk__el * 1e9) /                               \
                                ((double)tk__iters * (double)(ops_per_iter));  \
                if (tk__ns < tk__best) tk__best = tk__ns;                       \
            }                                                                   \
            tk_report_bench((tk), (label), tk__best, (double)(speedup));        \
        }                                                                       \
    } while (0)

/* ------------------------------------------------------------------ *
 * Module / test descriptors. A module is a named array of tests; each
 * test carries space-separated tags used for --tag filtering
 * (e.g. "fast primitives arith", "slow gmp").
 * ------------------------------------------------------------------ */
typedef void (*tk_test_fn)(tk_ctx *tk);

typedef struct {
    const char *name;
    tk_test_fn  fn;
    const char *tags;
} tk_test;

typedef struct {
    const char   *name;
    const char   *desc;
    const tk_test *tests;
    int           ntests;
} tk_module;

/* ------------------------------------------------------------------ *
 * Runner. Parses argv (see --help), runs the selected tests, prints a
 * per-module + grand-total table, returns 0 iff all checks passed.
 * ------------------------------------------------------------------ */
int tk_run(const tk_module *const *modules, int nmodules, int argc, char **argv);

#ifdef __cplusplus
}
#endif
#endif /* YAFU_TESTKIT_H */
