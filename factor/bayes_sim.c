/*
Copyright (c) 2026, Ben Buhrow
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.
*/

/*
 * ========== Note: created by Claude with prompting by bbuhrow ==========
 * bayessim(): replay pretest plans on inputs whose factorizations are known
 * (docs/ecm_bayes.md section 6).
 *
 * The real scheduler (schedule_work(), with its curve planning, level maxima,
 * P-1 states, trivial-input rule and -plan bayesian) picks every step.  Only
 * the work itself is simulated:
 *   - trial division and rho remove the primes below rho_bound;
 *   - P-1 finds p when p - 1 is B1-smooth but for one prime <= B2;
 *   - ECM finds p in the step where the accumulated hazard
 *     sum(curves * -ln(1 - P(curve finds p))) first passes an Exp(1) threshold
 *     drawn per (input, replicate, prime).  The same thresholds are used for
 *     every plan (common random numbers), so the plans are compared on the
 *     same luck;
 *   - the sieve finishes the composite left.
 * Costs come from a fitted model of the machine the inputs were run on
 * (docs/bayes_sim_prep.py writes it from a factor.json), so a whole benchmark
 * can be replayed in seconds and without timing noise.
 *
 * -plan bayesian can be given the model's true costs (scaled, to test how
 * wrong cost estimates hurt) through bayes_cost_hook, or use the tune data.
 *
 * Files (working directory):
 *   bayes_sim_in.txt     one input per line: <class> <prime> <prime> ...
 *                        (all prime factors, with multiplicity)
 *   bayes_sim_cost.txt   cost model, replicate count and plans (see
 *                        docs/ecm_bayes.md)
 *   bayes_sim_out.csv    per input and plan: mean costs over replicates
 *   bayes_sim_paths.csv  per input and plan: the steps when ECM finds nothing
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#if !defined(_MSC_VER)
#include <sys/time.h>
#endif
#include "factor.h"
#include "yafu_ecm.h"
#include "ecm_prob.h"
#include "ytools.h"
#include "gmp.h"
#include "autofactor_priv.h"

#define SIM_IN_FILE   "bayes_sim_in.txt"
#define SIM_COST_FILE "bayes_sim_cost.txt"
#define SIM_OUT_FILE  "bayes_sim_out.csv"
#define SIM_PATH_FILE "bayes_sim_paths.csv"

#define SIM_MAXP      48      // prime factors per input
#define SIM_MAXSIQS   512     // siqs table points
#define SIM_MAXPOL    32      // plans
#define SIM_MAXCLS    64      // input classes
#define SIM_HCACHE    64      // curve kinds per input in the hazard cache
#define SIM_PATHLEN   1024

typedef struct
{
	// ECM: per call overhead + curves * c (B1/1e5)^e (digits/70)^q seconds,
	// with exponent ecm_e_hi above ecm_b1max (the largest fitted B1)
	double ecm_c, ecm_e, ecm_q, ecm_step, ecm_b1max, ecm_e_hi;
	// P-1 levels: B1, B2, seconds = a + b * digits
	double pm1_B1[3], pm1_B2[3], pm1_a[3], pm1_b[3];
	int have_pm1[3];
	// sieve seconds vs digits (log-linear interpolation)
	double sq_d[SIM_MAXSIQS], sq_t[SIM_MAXSIQS];
	int nsq;
	double rho_bound;		// primes below this are found by tdiv / rho
	int lanes;				// AVX-ECM curves per thread are rounded up to this
	int reps;
	uint64_t seed;
} sim_cost_t;

enum { POL_PLAN, POL_BAYES_HOOK, POL_BAYES_TUNE };

typedef struct
{
	char name[64];
	int kind;				// POL_*
	int plan;				// PRETEST_* (POL_PLAN)
	double ratio;			// PRETEST_CUSTOM ratio
	double siqs_scale;		// POL_BAYES_HOOK: believed sieve time = scale * true
	double ecm_scale;		// POL_BAYES_HOOK: believed ECM time = scale * true
} sim_policy_t;

typedef struct
{
	int cls;
	int np;
	mpz_t p[SIM_MAXP];
	double pd[SIM_MAXP];
	int first[SIM_MAXP];			// index of the first equal prime
	signed char pm1[SIM_MAXP][3];	// P-1 level finds it: -1 not computed
} sim_input_t;

typedef struct
{
	ecm_curve_param_t cp;
	double h[SIM_MAXP];
} sim_hazard_t;

typedef struct
{
	sim_hazard_t e[SIM_HCACHE];
	int n, next;
} sim_hcache_t;

typedef struct
{
	double ecm, pm1, sieve, tlevel;
	char path[SIM_PATHLEN];
} sim_result_t;

typedef struct
{
	const sim_cost_t *C;
	double siqs_scale, ecm_scale;
	int threads;
} sim_hook_ctx_t;

// ---------------------------------------------------------------------------
// cost model
// ---------------------------------------------------------------------------
static double sim_sieve_sec(const sim_cost_t *C, double d)
{
	int i, n = C->nsq;
	double l0, l1, s;

	if (n == 0)
		return 0.0;
	if (n == 1)
		return C->sq_t[0];
	if (d <= C->sq_d[0])
		i = 0;
	else if (d >= C->sq_d[n - 1])
		i = n - 2;
	else
		for (i = 0; (i < n - 2) && (d > C->sq_d[i + 1]); i++);
	l0 = log(C->sq_t[i]);
	l1 = log(C->sq_t[i + 1]);
	s = (l1 - l0) / (C->sq_d[i + 1] - C->sq_d[i]);
	return exp(l0 + s * (d - C->sq_d[i]));
}

static double sim_ecm_spc(const sim_cost_t *C, double d, double B1)
{
	double s;
	if ((C->ecm_b1max > 0) && (B1 > C->ecm_b1max))
		s = pow(C->ecm_b1max / 1e5, C->ecm_e) * pow(B1 / C->ecm_b1max, C->ecm_e_hi);
	else
		s = pow(B1 / 1e5, C->ecm_e);
	return C->ecm_c * s * pow(d / 70.0, C->ecm_q);
}

static double sim_pm1_sec(const sim_cost_t *C, int lvl, double d)
{
	return C->pm1_a[lvl] + C->pm1_b[lvl] * d;
}

// believed costs for -plan bayesian (bayes_cost_hook)
static double hook_sieve(void *vctx, double digits)
{
	sim_hook_ctx_t *h = (sim_hook_ctx_t *)vctx;
	return h->siqs_scale * sim_sieve_sec(h->C, digits);
}

static double hook_ecm(void *vctx, double digits, double B1)
{
	// the per-call overhead spread over a minimum-size step
	sim_hook_ctx_t *h = (sim_hook_ctx_t *)vctx;
	double minc = (double)(h->C->lanes * (h->threads > 0 ? h->threads : 1));
	return h->ecm_scale * (sim_ecm_spc(h->C, digits, B1) + h->C->ecm_step / minc);
}

// ---------------------------------------------------------------------------
// P-1: p - 1 smooth to B1 but for one prime <= B2
// ---------------------------------------------------------------------------
static uint32_t *sim_primes = NULL;
static uint32_t sim_nprimes = 0, sim_plim = 0;

static void sim_primes_to(uint32_t lim)
{
	uint8_t *s;
	uint32_t i, j, n = 0;

	if (lim <= sim_plim)
		return;
	free(sim_primes);
	s = (uint8_t *)calloc((size_t)lim + 1, 1);
	for (i = 2; (uint64_t)i * i <= lim; i++)
		if (!s[i])
			for (j = i * i; j <= lim; j += i)
				s[j] = 1;
	for (i = 2; i <= lim; i++)
		n += !s[i];
	sim_primes = (uint32_t *)malloc((size_t)n * sizeof(uint32_t));
	n = 0;
	for (i = 2; i <= lim; i++)
		if (!s[i])
			sim_primes[n++] = i;
	free(s);
	sim_nprimes = n;
	sim_plim = lim;
}

static int sim_pm1_finds(mpz_t p, double B1, double B2)
{
	mpz_t m;
	uint32_t i;
	int found;

	sim_primes_to((uint32_t)B1);
	mpz_init(m);
	mpz_sub_ui(m, p, 1);
	for (i = 0; (i < sim_nprimes) && (sim_primes[i] <= B1); i++)
	{
		while (mpz_divisible_ui_p(m, sim_primes[i]))
			mpz_divexact_ui(m, m, sim_primes[i]);
		if (mpz_cmp_ui(m, 1) == 0)
			break;
	}
	found = (mpz_cmp_ui(m, 1) == 0) ||
		((mpz_cmp_d(m, B2) <= 0) && mpz_probab_prime_p(m, 10));
	mpz_clear(m);
	return found;
}

// ---------------------------------------------------------------------------
// ECM hazards and thresholds
// ---------------------------------------------------------------------------
static int sim_same_curve(const ecm_curve_param_t *a, const ecm_curve_param_t *b)
{
	return (a->B1 == b->B1) && (a->B2 == b->B2) && (a->nr == b->nr) &&
		(a->S == b->S) && (a->bsmodel == b->bsmodel) && (a->family == b->family);
}

// -ln(1 - P(one curve of kind cp finds p)) for each prime of the input
static const double *sim_hazards(sim_hcache_t *hc, const sim_input_t *in,
	const ecm_curve_param_t *cp)
{
	sim_hazard_t *e;
	int i;

	for (i = 0; i < hc->n; i++)
		if (sim_same_curve(&hc->e[i].cp, cp))
			return hc->e[i].h;
	if (hc->n < SIM_HCACHE)
		e = &hc->e[hc->n++];
	else
	{
		e = &hc->e[hc->next];
		hc->next = (hc->next + 1) % SIM_HCACHE;
	}
	e->cp = *cp;
	for (i = 0; i < in->np; i++)
	{
		double pr = ecmprob_prob(cp, in->pd[i]);
		e->h[i] = (pr >= 1.0) ? 1e300 : -log1p(-pr);
	}
	return e->h;
}

static uint64_t splitmix64(uint64_t *x)
{
	uint64_t z = (*x += 0x9E3779B97F4A7C15ull);
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
	return z ^ (z >> 31);
}

// Exp(1) threshold of prime k of input idx in replicate rep
static double sim_threshold(uint64_t seed, int idx, int rep, int k)
{
	uint64_t x = seed;
	double u;
	x ^= splitmix64(&x) + (uint64_t)(uint32_t)idx * 0x100000001B3ull;
	x ^= splitmix64(&x) + (uint64_t)(uint32_t)rep * 0xC2B2AE3D27D4EB4Full;
	x ^= splitmix64(&x) + (uint64_t)(uint32_t)k * 0x165667B19E3779F9ull;
	u = ((double)(splitmix64(&x) >> 11) + 0.5) * (1.0 / 9007199254740992.0);
	return -log(u);
}

// ---------------------------------------------------------------------------
// one factorization
// ---------------------------------------------------------------------------
static int is_ecm_state(enum factorization_state s)
{
	return (s >= state_ecm_15digit) && (s <= state_ecm_65digit);
}

static void sim_reset_ecm(fact_obj_t *fobj)
{
	int i;
	if (fobj->ecm_obj.num_records > 0)
	{
		free(fobj->ecm_obj.curve_b1_rec);
		free(fobj->ecm_obj.curve_b2_rec);
		free(fobj->ecm_obj.param_rec);
		free(fobj->ecm_obj.num_rec);
	}
	fobj->ecm_obj.num_records = 0;
	for (i = 0; i < ECMPROB_NDIG; i++)
		fobj->ecm_obj.twork[i] = 0.0;
	for (i = 0; i < NUM_ECM_LEVELS; i++)
		fobj->ecm_obj.tlevels[i] = 0.0;
	fobj->ecm_obj.total_work = 0.0;
}

static uint32_t sim_round_curves(fact_obj_t *fobj, const sim_cost_t *C,
	uint32_t n, uint64_t B1)
{
	// as AVX-ECM runs them (vec_ecm_main): at least one per thread, the same
	// number per thread, a whole vector per thread.  GMP-ECM runs exactly n.
	uint32_t t = (uint32_t)(fobj->THREADS > 0 ? fobj->THREADS : 1), per;
	if (fobj->ecm_obj.prefer_gmpecm || (B1 > fobj->ecm_obj.ecm_ext_xover) || (C->lanes <= 1))
		return n;
	if (n < t)
		n = t;
	per = n / t + ((n % t) != 0);
	per = ((per + (uint32_t)C->lanes - 1) / (uint32_t)C->lanes) * (uint32_t)C->lanes;
	return per * t;
}

static void sim_one(fact_obj_t *fobj, const sim_cost_t *C, sim_input_t *in,
	int idx, int rep, int allfail, sim_hcache_t *hc, sim_result_t *R)
{
	factor_work_t fwork;
	int alive[SIM_MAXP];
	double H[SIM_MAXP], E[SIM_MAXP];
	double tbound;
	mpz_t b;
	int k, iter, plen = 0;

	memset(R, 0, sizeof(*R));
	mpz_init(b);
	sim_reset_ecm(fobj);
	memset(&fwork, 0, sizeof(fwork));
	fwork.initial_work = 0.0;
	init_factor_work(&fwork, fobj);

	// trial division, fermat and rho are done
	fwork.tdiv_limit = fwork.tdiv_max_limit;
	fwork.fermat_iterations = fwork.fermat_max_iterations;
	fwork.rho_bases = fwork.rho_max_bases;
	tbound = MAX(C->rho_bound, (double)fwork.tdiv_max_limit);

	for (k = 0; k < in->np; k++)
	{
		alive[k] = (in->pd[k] >= tbound);
		H[k] = 0.0;
		E[k] = allfail ? 1e300 : sim_threshold(C->seed, idx, rep, in->first[k]);
	}

	for (iter = 0; iter < 10000; iter++)
	{
		enum factorization_state st;
		int nalive = 0, digits, found = 0;

		mpz_set_ui(b, 1);
		for (k = 0; k < in->np; k++)
			if (alive[k])
			{
				mpz_mul(b, b, in->p[k]);
				nalive++;
			}
		if (nalive <= 1)
			break;

		digits = gmp_base10(b);
		mpz_set(fobj->N, b);
		mpz_set(fobj->ecm_obj.gmp_n, b);
		// the job type, as determine_job_type() would find it for a general
		// input (skips its snfs search)
		if (fwork.target_job_type == job_unknown)
			fwork.target_job_type = ((double)digits < fobj->autofact_obj.qs_gnfs_xover) ?
				job_siqs : job_gnfs;

		st = schedule_work(&fwork, b, fobj);

		if (is_ecm_state(st))
		{
			uint64_t B1 = fwork.B1;
			uint32_t n = sim_round_curves(fobj, C, fwork.curves, B1);
			ecm_curve_param_t cp;
			const double *h;

			if (n == 0)
				break;
			ecm_plan_curve_params(fobj, B1, &cp);
			record_curves_completed(&fobj->ecm_obj, (int)n, &cp);
			set_ecm_curves_done(&fwork, st, get_ecm_curves_done(&fwork, st) + n);
			R->ecm += C->ecm_step + (double)n * sim_ecm_spc(C, (double)digits, (double)B1);
			if (plen < SIM_PATHLEN - 32)
				plen += sprintf(R->path + plen, "%s%" PRIu64 ":%u", plen ? ";" : "", B1, n);

			h = sim_hazards(hc, in, &cp);
			for (k = 0; k < in->np; k++)
			{
				if (!alive[k])
					continue;
				H[k] += (double)n * h[k];
				if (H[k] >= E[k])
				{
					alive[k] = 0;
					found = 1;
				}
			}
		}
		else if ((st == state_pm1_lvl1) || (st == state_pm1_lvl2) || (st == state_pm1_lvl3))
		{
			int lvl = (st == state_pm1_lvl1) ? 0 : (st == state_pm1_lvl2) ? 1 : 2;
			if (lvl == 0) fwork.pm1_lvl1_curves = fwork.pm1_max_lvl1_curves;
			else if (lvl == 1) fwork.pm1_lvl2_curves = fwork.pm1_max_lvl2_curves;
			else fwork.pm1_lvl3_curves = fwork.pm1_max_lvl3_curves;
			if (!C->have_pm1[lvl])
				continue;
			R->pm1 += sim_pm1_sec(C, lvl, (double)digits);
			if (plen < SIM_PATHLEN - 32)
				plen += sprintf(R->path + plen, "%spm1_%.0f", plen ? ";" : "", C->pm1_B1[lvl]);
			for (k = 0; k < in->np; k++)
			{
				int f = in->first[k];
				if (!alive[k])
					continue;
				if (in->pm1[f][lvl] < 0)
					in->pm1[f][lvl] = (signed char)sim_pm1_finds(in->p[f],
						C->pm1_B1[lvl], C->pm1_B2[lvl]);
				if (in->pm1[f][lvl])
				{
					alive[k] = 0;
					found = 1;
				}
			}
		}
		else if ((st == state_pp1_lvl1) || (st == state_pp1_lvl2) || (st == state_pp1_lvl3))
		{
			// not simulated (none by default): mark done
			fwork.pp1_lvl1_curves = fwork.pp1_max_lvl1_curves;
			if (st != state_pp1_lvl1) fwork.pp1_lvl2_curves = fwork.pp1_max_lvl2_curves;
			if (st == state_pp1_lvl3) fwork.pp1_lvl3_curves = fwork.pp1_max_lvl3_curves;
		}
		else if ((st == state_qs) || (st == state_nfs))
		{
			R->sieve += sim_sieve_sec(C, (double)digits);
			if (plen < SIM_PATHLEN - 32)
				plen += sprintf(R->path + plen, "%ssieve_c%d", plen ? ";" : "", digits);
			break;
		}
		else
		{
			// done / idle / unexpected: stop
			break;
		}

		if (found)
		{
			// equal primes go together
			for (k = 0; k < in->np; k++)
				if (alive[k] && !alive[in->first[k]])
					alive[k] = 0;
			fwork.target_job_type = job_unknown;
		}
	}

	R->tlevel = fobj->ecm_obj.total_work;
	mpz_clear(b);
}

// ---------------------------------------------------------------------------
// files
// ---------------------------------------------------------------------------
static int sim_read_cost(const char *fname, sim_cost_t *C, sim_policy_t *pol, int *npol)
{
	FILE *f = fopen(fname, "r");
	char line[1024];
	int i, j;

	memset(C, 0, sizeof(*C));
	C->rho_bound = 1e6;
	C->lanes = 8;
	C->reps = 16;
	C->seed = 12345;
	C->ecm_e = 1.0;
	C->ecm_e_hi = 1.0;
	*npol = 0;

	if (f == NULL)
	{
		printf("bayessim: could not open %s\n", fname);
		return 0;
	}
	while (fgets(line, sizeof(line), f) != NULL)
	{
		char key[64], a1[64];
		double v[4];
		char *h = strchr(line, '#');
		if (h) *h = '\0';
		if (sscanf(line, "%63s", key) != 1)
			continue;

		if (!strcmp(key, "rho_bound") && (sscanf(line, "%*s %lf", &v[0]) == 1)) C->rho_bound = v[0];
		else if (!strcmp(key, "lanes") && (sscanf(line, "%*s %lf", &v[0]) == 1)) C->lanes = (int)v[0];
		else if (!strcmp(key, "reps") && (sscanf(line, "%*s %lf", &v[0]) == 1)) C->reps = (int)v[0];
		else if (!strcmp(key, "seed") && (sscanf(line, "%*s %lf", &v[0]) == 1)) C->seed = (uint64_t)v[0];
		else if (!strcmp(key, "ecm_b1max") && (sscanf(line, "%*s %lf", &v[0]) == 1)) C->ecm_b1max = v[0];
		else if (!strcmp(key, "ecm_e_hi") && (sscanf(line, "%*s %lf", &v[0]) == 1)) C->ecm_e_hi = v[0];
		else if (!strcmp(key, "ecm_step") && (sscanf(line, "%*s %lf", &v[0]) == 1)) C->ecm_step = v[0];
		else if (!strcmp(key, "ecm_curve") &&
			(sscanf(line, "%*s %lf %lf %lf", &v[0], &v[1], &v[2]) == 3))
		{
			C->ecm_c = v[0]; C->ecm_e = v[1]; C->ecm_q = v[2];
		}
		else if (!strcmp(key, "pm1") &&
			(sscanf(line, "%*s %lf %lf %lf %lf", &v[0], &v[1], &v[2], &v[3]) == 4))
		{
			int lvl = (v[0] < 1e6) ? 0 : (v[0] < 1e7) ? 1 : 2;
			C->pm1_B1[lvl] = v[0]; C->pm1_B2[lvl] = v[1];
			C->pm1_a[lvl] = v[2]; C->pm1_b[lvl] = v[3];
			C->have_pm1[lvl] = 1;
		}
		else if (!strcmp(key, "siqs") && (sscanf(line, "%*s %lf %lf", &v[0], &v[1]) == 2))
		{
			if ((C->nsq < SIM_MAXSIQS) && (v[1] > 0))
			{
				C->sq_d[C->nsq] = v[0];
				C->sq_t[C->nsq] = v[1];
				C->nsq++;
			}
		}
		else if (!strcmp(key, "policy") && (*npol < SIM_MAXPOL))
		{
			// policy <name> normal|light|deep|custom <ratio>|bayes <siqs_scale> <ecm_scale>|bayes_tune
			sim_policy_t *p = &pol[*npol];
			char kind[64];
			memset(p, 0, sizeof(*p));
			if (sscanf(line, "%*s %63s %63s", p->name, kind) != 2)
				continue;
			p->siqs_scale = 1.0;
			p->ecm_scale = 1.0;
			if (!strcmp(kind, "normal")) { p->kind = POL_PLAN; p->plan = PRETEST_NORMAL; }
			else if (!strcmp(kind, "light")) { p->kind = POL_PLAN; p->plan = PRETEST_LIGHT; }
			else if (!strcmp(kind, "deep")) { p->kind = POL_PLAN; p->plan = PRETEST_DEEP; }
			else if (!strcmp(kind, "noecm")) { p->kind = POL_PLAN; p->plan = PRETEST_NOECM; }
			else if (!strcmp(kind, "custom") &&
				(sscanf(line, "%*s %*s %*s %lf", &v[0]) == 1))
			{
				p->kind = POL_PLAN; p->plan = PRETEST_CUSTOM; p->ratio = v[0];
			}
			else if (!strcmp(kind, "bayes"))
			{
				int n = sscanf(line, "%*s %*s %*s %lf %lf", &v[0], &v[1]);
				p->kind = POL_BAYES_HOOK;
				if (n >= 1) p->siqs_scale = v[0];
				if (n >= 2) p->ecm_scale = v[1];
			}
			else if (!strcmp(kind, "bayes_tune")) { p->kind = POL_BAYES_TUNE; }
			else
			{
				printf("bayessim: unknown policy kind '%s'\n", kind);
				continue;
			}
			(*npol)++;
		}
		else
		{
			sscanf(line, "%63s", a1);
			printf("bayessim: ignoring cost line '%s ...'\n", a1);
		}
	}
	fclose(f);

	// sort the sieve table by size
	for (i = 1; i < C->nsq; i++)
		for (j = i; (j > 0) && (C->sq_d[j - 1] > C->sq_d[j]); j--)
		{
			double t = C->sq_d[j]; C->sq_d[j] = C->sq_d[j - 1]; C->sq_d[j - 1] = t;
			t = C->sq_t[j]; C->sq_t[j] = C->sq_t[j - 1]; C->sq_t[j - 1] = t;
		}

	if (*npol == 0)
	{
		strcpy(pol[0].name, "normal"); pol[0].kind = POL_PLAN; pol[0].plan = PRETEST_NORMAL;
		strcpy(pol[1].name, "bayes"); pol[1].kind = POL_BAYES_HOOK;
		pol[1].siqs_scale = 1.0; pol[1].ecm_scale = 1.0;
		*npol = 2;
	}
	if ((C->nsq == 0) || (C->ecm_c <= 0))
	{
		printf("bayessim: %s needs siqs and ecm_curve lines\n", fname);
		return 0;
	}
	return 1;
}

static sim_input_t *sim_read_inputs(const char *fname, int *nin)
{
	FILE *f = fopen(fname, "r");
	char *line;
	int alloc = 1024, n = 0;
	sim_input_t *in;
	const int linelen = 1 << 16;

	*nin = 0;
	if (f == NULL)
	{
		printf("bayessim: could not open %s\n", fname);
		return NULL;
	}
	line = (char *)malloc((size_t)linelen);
	in = (sim_input_t *)malloc((size_t)alloc * sizeof(sim_input_t));
	while (fgets(line, linelen, f) != NULL)
	{
		char *tok;
		sim_input_t *x;
		int k, j;

		tok = strtok(line, " \t\r\n,");
		if ((tok == NULL) || (tok[0] == '#'))
			continue;
		if (n == alloc)
		{
			alloc *= 2;
			in = (sim_input_t *)realloc(in, (size_t)alloc * sizeof(sim_input_t));
		}
		x = &in[n];
		x->cls = atoi(tok);
		x->np = 0;
		while (((tok = strtok(NULL, " \t\r\n,")) != NULL) && (x->np < SIM_MAXP))
		{
			mpz_init(x->p[x->np]);
			if (mpz_set_str(x->p[x->np], tok, 10) != 0)
			{
				mpz_clear(x->p[x->np]);
				continue;
			}
			x->pd[x->np] = mpz_get_d(x->p[x->np]);
			x->np++;
		}
		for (k = 0; k < x->np; k++)
		{
			x->first[k] = k;
			for (j = 0; j < k; j++)
				if (mpz_cmp(x->p[j], x->p[k]) == 0)
				{
					x->first[k] = j;
					break;
				}
			x->pm1[k][0] = x->pm1[k][1] = x->pm1[k][2] = -1;
		}
		if (x->np > 0)
			n++;
	}
	free(line);
	fclose(f);
	*nin = n;
	return in;
}

// ---------------------------------------------------------------------------
// bayessim()
// ---------------------------------------------------------------------------
void bayes_sim(fact_obj_t *fobj)
{
	sim_cost_t C;
	sim_policy_t pol[SIM_MAXPOL];
	sim_input_t *in;
	int npol, nin, i, j, r, c;
	double *tot;		// [npol][nin] mean total per input
	double *parts;		// [npol][nin][4] mean ecm, pm1, sieve, tlevel
	int cls[SIM_MAXCLS], ncls = 0;
	FILE *out, *paths;
	sim_hcache_t *hc;
	sim_hook_ctx_t hctx;
	bayes_cost_hook_t hook;
	struct timeval t0, t1;

	// what the replay changes in fobj
	int save_v = fobj->VFLAG, save_log = fobj->LOGFLAG;
	int save_plan = fobj->autofact_obj.yafu_pretest_plan;
	double save_ratio = fobj->autofact_obj.target_pretest_ratio;
	char save_logname[1024];

	if (!sim_read_cost(SIM_COST_FILE, &C, pol, &npol))
		return;
	in = sim_read_inputs(SIM_IN_FILE, &nin);
	if ((in == NULL) || (nin == 0))
	{
		printf("bayessim: no inputs in %s\n", SIM_IN_FILE);
		free(in);
		return;
	}

	for (i = 0; i < nin; i++)
	{
		for (c = 0; c < ncls; c++)
			if (cls[c] == in[i].cls)
				break;
		if ((c == ncls) && (ncls < SIM_MAXCLS))
			cls[ncls++] = in[i].cls;
	}
	for (i = 1; i < ncls; i++)
		for (j = i; (j > 0) && (cls[j - 1] > cls[j]); j--)
		{
			int t = cls[j]; cls[j] = cls[j - 1]; cls[j - 1] = t;
		}

	printf("bayessim: %d inputs in %d classes, %d plans, %d replicates, %d threads\n",
		nin, ncls, npol, C.reps, fobj->THREADS);

	tot = (double *)calloc((size_t)npol * (size_t)nin, sizeof(double));
	parts = (double *)calloc((size_t)npol * (size_t)nin * 4, sizeof(double));
	hc = (sim_hcache_t *)malloc(sizeof(sim_hcache_t));
	out = fopen(SIM_OUT_FILE, "w");
	paths = fopen(SIM_PATH_FILE, "w");
	if (out) fprintf(out, "plan,class,index,ecm,pm1,sieve,total,tlevel\n");
	if (paths) fprintf(paths, "plan,class,index,path\n");

	strcpy(save_logname, fobj->flogname);
	fobj->flogname[0] = '\0';
	fobj->VFLAG = -1;
	fobj->LOGFLAG = 0;
	ecmprob_init();

	hctx.C = &C;
	hctx.threads = fobj->THREADS;
	hook.sieve_sec = hook_sieve;
	hook.ecm_sec_per_curve = hook_ecm;
	hook.ctx = &hctx;

	gettimeofday(&t0, NULL);
	for (j = 0; j < npol; j++)
	{
		sim_policy_t *p = &pol[j];

		bayes_cost_hook = NULL;
		if (p->kind == POL_PLAN)
		{
			fobj->autofact_obj.yafu_pretest_plan = p->plan;
			if (p->plan == PRETEST_CUSTOM)
				fobj->autofact_obj.target_pretest_ratio = p->ratio;
		}
		else
		{
			fobj->autofact_obj.yafu_pretest_plan = PRETEST_BAYES;
			if (p->kind == POL_BAYES_HOOK)
			{
				hctx.siqs_scale = p->siqs_scale;
				hctx.ecm_scale = p->ecm_scale;
				hook.gen = 1000 + j;
				bayes_cost_hook = &hook;
			}
		}

		for (i = 0; i < nin; i++)
		{
			sim_result_t R;
			double s[4] = { 0, 0, 0, 0 };

			hc->n = hc->next = 0;
			for (r = 0; r < C.reps; r++)
			{
				sim_one(fobj, &C, &in[i], i, r, 0, hc, &R);
				s[0] += R.ecm; s[1] += R.pm1; s[2] += R.sieve; s[3] += R.tlevel;
			}
			for (c = 0; c < 4; c++)
				parts[((size_t)j * (size_t)nin + (size_t)i) * 4 + (size_t)c] = s[c] / C.reps;
			tot[(size_t)j * (size_t)nin + (size_t)i] = (s[0] + s[1] + s[2]) / C.reps;
			if (out)
				fprintf(out, "%s,%d,%d,%.6g,%.6g,%.6g,%.6g,%.3f\n", p->name, in[i].cls, i,
					s[0] / C.reps, s[1] / C.reps, s[2] / C.reps,
					(s[0] + s[1] + s[2]) / C.reps, s[3] / C.reps);

			// the path when ECM finds nothing
			sim_one(fobj, &C, &in[i], i, 0, 1, hc, &R);
			if (paths)
				fprintf(paths, "%s,%d,%d,%s\n", p->name, in[i].cls, i, R.path);
		}
		gettimeofday(&t1, NULL);
		printf("bayessim: plan %s done (%1.1f s)\n", p->name, ytools_difftime(&t0, &t1));
	}
	bayes_cost_hook = NULL;

	// summary: per class, each plan's mean time per input and the paired
	// difference from the first plan
	printf("\nmean seconds per input (simulated: ECM + P-1 + sieve; tdiv/rho excluded)\n");
	printf("%5s %-14s %9s %9s %9s %9s %7s %10s %8s\n", "class", "plan", "total",
		"ecm", "pm1", "sieve", "t-lvl", "vs first", "se");
	for (c = 0; c < ncls; c++)
	{
		for (j = 0; j < npol; j++)
		{
			double sum = 0, sp[4] = { 0, 0, 0, 0 }, sd = 0, sd2 = 0, base = 0;
			int n = 0, k;
			for (i = 0; i < nin; i++)
			{
				double d;
				if (in[i].cls != cls[c])
					continue;
				n++;
				sum += tot[(size_t)j * (size_t)nin + (size_t)i];
				base += tot[i];
				for (k = 0; k < 4; k++)
					sp[k] += parts[((size_t)j * (size_t)nin + (size_t)i) * 4 + (size_t)k];
				d = tot[(size_t)j * (size_t)nin + (size_t)i] - tot[i];
				sd += d;
				sd2 += d * d;
			}
			if (n == 0)
				continue;
			{
				double md = sd / n;
				double se = (n > 1) ? sqrt((sd2 / n - md * md) / (n - 1)) : 0.0;
				printf("C%-4d %-14s %9.4f %9.4f %9.4f %9.4f %7.2f %+9.2f%% %7.2f%%\n",
					cls[c], pol[j].name, sum / n, sp[0] / n, sp[1] / n, sp[2] / n, sp[3] / n,
					(base > 0) ? 100.0 * sd / base : 0.0,
					(base > 0) ? 100.0 * se * n / base : 0.0);
			}
		}
	}
	printf("wrote %s and %s\n", SIM_OUT_FILE, SIM_PATH_FILE);

	if (out) fclose(out);
	if (paths) fclose(paths);
	strcpy(fobj->flogname, save_logname);
	fobj->VFLAG = save_v;
	fobj->LOGFLAG = save_log;
	fobj->autofact_obj.yafu_pretest_plan = save_plan;
	fobj->autofact_obj.target_pretest_ratio = save_ratio;
	sim_reset_ecm(fobj);
	mpz_set_ui(fobj->ecm_obj.gmp_n, 0);

	for (i = 0; i < nin; i++)
		for (r = 0; r < in[i].np; r++)
			mpz_clear(in[i].p[r]);
	free(in);
	free(tot);
	free(parts);
	free(hc);
	free(sim_primes);
	sim_primes = NULL;
	sim_nprimes = sim_plim = 0;
	return;
}
