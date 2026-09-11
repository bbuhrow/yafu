/*--------------------------------------------------------------------
Shared, thread-safe progress statistics for polynomial selection.

One instance per poly-select run, reachable from every engine (via the
driver's handle_collision) and from the stage-2 callbacks. Because all
engines and both stage-2 phases funnel through the same instance, the
progress output is identical regardless of which stage-1 engine ran.

The funnel it tracks:
    stage-1 hits  ->  survive size-opt  ->  survive root-opt (saved)
plus the best Murphy-E seen so far and coefficient progress.

Reporting is a single job-wide roll-up line (verbosity 1), refreshed in
place on stderr, printed by whichever thread wins a throttled check while
holding the lock -- so there is exactly one writer at a time and the line
never interleaves. Verbosity 2 suppresses the roll-up (the driver's
per-thread detail owns the terminal instead); verbosity 0 is silent.

The a_d is passed in as a pre-formatted string so this module carries no
GMP dependency; callers gmp_snprintf the mpz at the (infrequent) hooks.
--------------------------------------------------------------------*/

#ifndef _GNFS_POLY_POLY_STATS_H_
#define _GNFS_POLY_POLY_STATS_H_

#include "poly.h"      /* uint32/uint64, msieve_gettimeofday, struct timeval */
#include "thread.h"    /* mutex_t, mutex_init/lock/unlock/free              */

#ifdef __cplusplus
extern "C" {
#endif

#define POLY_STATS_ADSTRLEN 64

typedef struct {
	mutex_t lock;

	/* the funnel (monotone counters) */
	uint64 hits;          /* stage-1 collisions submitted to stage 2 */
	uint64 sizeopt_pass;  /* survivors of size optimization          */
	uint64 rootopt_pass;  /* survivors of root optimization (saved)  */
	uint64 stage2_done;   /* stage-2 tasks completed (for queue depth) */

	/* best quality so far */
	double best_e;
	char   best_ad[POLY_STATS_ADSTRLEN];
	int    have_best;

	/* coefficient progress */
	char   cur_ad[POLY_STATS_ADSTRLEN];
	uint64 ad_done;
	uint64 ad_total;      /* 0 if unknown */

	/* reporting */
	int    vflag;         /* 0 silent, 1 roll-up, 2 detail (roll-up off) */
	double interval;      /* seconds between roll-up refreshes           */
	struct timeval start;
	struct timeval last_print;
	int    line_pending;  /* a \r line is on screen awaiting a newline    */

	/* early-abort control (0/NULL = disabled); checked on each saved poly */
	msieve_obj *obj;      /* raise MSIEVE_FLAG_STOP_SIEVING here          */
	double e_threshold;   /* abort once a saved poly's E >= this          */
	uint64 max_polys;     /* abort once rootopt_pass >= this              */
} poly_stage_stats_t;

/* ad_total may be 0 if the count of leading coefficients isn't known up
   front; the progress fraction is then omitted. */
void poly_stats_init(poly_stage_stats_t *s, int vflag, uint64 ad_total);
void poly_stats_free(poly_stage_stats_t *s);

/* driver marks the coefficient currently being searched */
void poly_stats_set_ad(poly_stage_stats_t *s, const char *ad_str,
			uint64 ad_index);

/* the three funnel hooks; safe to call from any thread */
void poly_stats_add_hit(poly_stage_stats_t *s);
void poly_stats_add_sizeopt(poly_stage_stats_t *s);
void poly_stats_add_rootopt(poly_stage_stats_t *s, double combined_e,
			const char *ad_str);

/* stage-2 worker calls this when it finishes a hit's task; the live queue
   depth is (hits submitted - tasks completed). */
void poly_stats_stage2_done(poly_stage_stats_t *s);

/* force!=0 prints unconditionally and terminates the line (use at end of
   run); otherwise the roll-up is throttled to `interval`. */
void poly_stats_report(poly_stage_stats_t *s, int force);

/* arm early abort: when a saved poly reaches e_threshold, or rootopt_pass
   reaches max_polys, MSIEVE_FLAG_STOP_SIEVING is raised on obj and every
   engine drains. 0 for either threshold disables it. */
void poly_stats_set_abort(poly_stage_stats_t *s, msieve_obj *obj,
			double e_threshold, uint64 max_polys);

#ifdef __cplusplus
}
#endif

#endif /* !_GNFS_POLY_POLY_STATS_H_ */
