/*--------------------------------------------------------------------
Implementation of the shared poly-select progress statistics.
See poly_stats.h.
--------------------------------------------------------------------*/

#include "poly_stats.h"
#include <stdio.h>
#include <string.h>

/* how often (in hits) to even check whether a roll-up is due; keeps the
   common hit path from calling the clock on every collision */
#define HIT_REPORT_MASK 0x0FFFu

static double
stats_elapsed(struct timeval *a, struct timeval *b)
{
	return (double)(b->tv_sec - a->tv_sec) +
	       (double)(b->tv_usec - a->tv_usec) * 1e-6;
}

/* 1234567 -> "1.23M", 47100 -> "47.1k", 312 -> "312" */
static void
fmt_count(char *buf, size_t n, uint64 v)
{
	if (v < 1000)
		snprintf(buf, n, "%llu", (unsigned long long)v);
	else if (v < 1000000)
		snprintf(buf, n, "%.1fk", (double)v / 1e3);
	else if (v < 1000000000ULL)
		snprintf(buf, n, "%.2fM", (double)v / 1e6);
	else
		snprintf(buf, n, "%.2fG", (double)v / 1e9);
}

/* assumes the lock is held */
static void
report_locked(poly_stage_stats_t *s, int force)
{
	struct timeval now;
	double t;
	char hbuf[16], sbuf[16], rbuf[16], adbuf[96], bestbuf[96];
	uint64 depth;

	if (s->vflag < 1)
		return;                       /* silent */

	msieve_gettimeofday(&now, NULL);
	if (!force && stats_elapsed(&s->last_print, &now) < s->interval)
		return;                       /* throttled */
	s->last_print = now;
	t = stats_elapsed(&s->start, &now);

	fmt_count(hbuf, sizeof(hbuf), s->hits);
	fmt_count(sbuf, sizeof(sbuf), s->sizeopt_pass);
	fmt_count(rbuf, sizeof(rbuf), s->rootopt_pass);
	depth = (s->hits >= s->stage2_done) ? (s->hits - s->stage2_done) : 0;

	if (s->ad_total)
		snprintf(adbuf, sizeof(adbuf), "a_d %s (%llu/%llu)",
			s->cur_ad[0] ? s->cur_ad : "-",
			(unsigned long long)s->ad_done,
			(unsigned long long)s->ad_total);
	else
		snprintf(adbuf, sizeof(adbuf), "a_d %s",
			s->cur_ad[0] ? s->cur_ad : "-");

	if (s->have_best)
		snprintf(bestbuf, sizeof(bestbuf), "best E %.3le @ a_d %s",
			s->best_e, s->best_ad);
	else
		snprintf(bestbuf, sizeof(bestbuf), "best E -");

	/* pad with trailing spaces so a shrinking line fully overwrites the
	   previous one; \r returns to column 0 without a newline */
	fprintf(stderr, "\r%s | hits %s  cand %s  saved %s | q %llu | %s | %.0fs      ",
		adbuf, hbuf, sbuf, rbuf, (unsigned long long)depth, bestbuf, t);
	if (force || s->vflag >= 2)
		fprintf(stderr, "\n");   /* -vv scrolls; final line terminates */
	else
		fflush(stderr);          /* -v overwrites in place */
	s->line_pending = !(force || s->vflag >= 2);
}

void
poly_stats_init(poly_stage_stats_t *s, int vflag, uint64 ad_total)
{
	memset(s, 0, sizeof(*s));
	mutex_init(&s->lock);
	s->vflag = vflag;
	s->ad_total = ad_total;
	s->interval = 0.5;
	s->best_e = 0.0;
	msieve_gettimeofday(&s->start, NULL);
	s->last_print = s->start;
}

void
poly_stats_free(poly_stage_stats_t *s)
{
	/* leave the final line terminated if one is still pending */
	if (s->line_pending)
		fprintf(stderr, "\n");
	mutex_free(&s->lock);
}

void
poly_stats_set_ad(poly_stage_stats_t *s, const char *ad_str, uint64 ad_index)
{
	mutex_lock(&s->lock);
	if (ad_str) {
		strncpy(s->cur_ad, ad_str, POLY_STATS_ADSTRLEN - 1);
		s->cur_ad[POLY_STATS_ADSTRLEN - 1] = 0;
	}
	s->ad_done = ad_index;
	report_locked(s, 0);
	mutex_unlock(&s->lock);
}

void
poly_stats_add_hit(poly_stage_stats_t *s)
{
	mutex_lock(&s->lock);
	s->hits++;
	if ((s->hits & HIT_REPORT_MASK) == 0)
		report_locked(s, 0);
	mutex_unlock(&s->lock);
}

void
poly_stats_add_sizeopt(poly_stage_stats_t *s)
{
	mutex_lock(&s->lock);
	s->sizeopt_pass++;
	report_locked(s, 0);
	mutex_unlock(&s->lock);
}

void
poly_stats_add_rootopt(poly_stage_stats_t *s, double combined_e,
			const char *ad_str)
{
	mutex_lock(&s->lock);
	s->rootopt_pass++;
	if (!s->have_best || combined_e > s->best_e) {
		s->best_e = combined_e;
		s->have_best = 1;
		if (ad_str) {
			strncpy(s->best_ad, ad_str, POLY_STATS_ADSTRLEN - 1);
			s->best_ad[POLY_STATS_ADSTRLEN - 1] = 0;
		}
	}

	/* early abort: a good-enough poly, or enough of them. Monotone flag
	   set; every engine reads it locklessly in its inner loop and drains. */
	if (s->obj != NULL &&
	    ((s->max_polys   != 0   && s->rootopt_pass >= s->max_polys) ||
	     (s->e_threshold > 0.0  && combined_e     >= s->e_threshold)))
		s->obj->flags |= MSIEVE_FLAG_STOP_SIEVING;

	report_locked(s, 0);
	mutex_unlock(&s->lock);
}

void
poly_stats_stage2_done(poly_stage_stats_t *s)
{
	mutex_lock(&s->lock);
	s->stage2_done++;
	report_locked(s, 0);
	mutex_unlock(&s->lock);
}

void
poly_stats_report(poly_stage_stats_t *s, int force)
{
	mutex_lock(&s->lock);
	report_locked(s, force);
	mutex_unlock(&s->lock);
}

void
poly_stats_set_abort(poly_stage_stats_t *s, msieve_obj *obj,
			double e_threshold, uint64 max_polys)
{
	mutex_lock(&s->lock);
	s->obj = obj;
	s->e_threshold = e_threshold;
	s->max_polys = max_polys;
	mutex_unlock(&s->lock);
}
