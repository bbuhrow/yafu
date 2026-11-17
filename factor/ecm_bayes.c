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
 * Bayesian ECM-vs-sieve stopping rule.  See include/ecm_bayes.h for the
 * model and docs/ecm_bayes.md for the derivation and the test driver.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ecm_bayes.h"

static double band_mass(double lo, double hi)
{
    // expected number of primes p with lo <= log10 p < hi in a random
    // integer (Mertens): ln(hi / lo)
    if (hi <= lo || lo <= 0)
        return 0.0;
    return log(hi / lo);
}

void ecmbayes_prior(double *lambda, int type, double digits, double dmin,
    double margin)
{
    int i;
    double half = digits / 2.0;     // the smallest factor is <= sqrt(N)
    double lo0 = dmin;

    if (lo0 < (double)(ECMPROB_DMIN - 1))
        lo0 = (double)(ECMPROB_DMIN - 1);
    if (type == ECMBAYES_PRIOR_HARD)
    {
        double h = half - (margin > 0 ? margin : 1.0);
        if (h > lo0)
            lo0 = h;
    }

    for (i = 0; i < ECMPROB_NDIG; i++)
    {
        double d = (double)(ECMPROB_DMIN + i);
        double lo = d - 1.0, hi = d;
        if (lo < lo0) lo = lo0;
        if (hi > half) hi = half;
        lambda[i] = band_mass(lo, hi);
    }

    // above the grid
    {
        double lo = (double)ECMPROB_DMAX;
        if (lo < lo0) lo = lo0;
        lambda[ECMPROB_NDIG] = band_mass(lo, half);
    }
}

// remaining (undetected) mass per band and its total
static double post_mass(const double *lambda, const double *W, double *m)
{
    int i;
    double L = lambda[ECMPROB_NDIG];
    for (i = 0; i < ECMPROB_NDIG; i++)
    {
        m[i] = lambda[i] * exp(-W[i]);
        L += m[i];
    }
    return L;
}

double ecmbayes_pfind(const double *lambda, const double *W, const double *w)
{
    int i;
    double m[ECMPROB_NDIG];
    double L = post_mass(lambda, W, m);
    double Ls = 0.0;

    for (i = 0; i < ECMPROB_NDIG; i++)
        Ls += m[i] * -expm1(-w[i]);

    if (L <= 0.0)
        return 0.0;
    // P(the step finds >= 1 | >= 1 undetected factor exists)
    return -expm1(-Ls) / -expm1(-L);
}

// expected cost of a c-digit cofactor: P(composite) * T(c), where
// P(composite) comes from the posterior mass at or below c/2 digits
static double cofactor_cost(const double *m, double tail_lambda, double c,
    ecmbayes_time_fn T, void *ctx)
{
    int i;
    double half = c / 2.0;
    double L = 0.0;

    if (c < 2.0)
        return 0.0;

    for (i = 0; i < ECMPROB_NDIG; i++)
    {
        double d = (double)(ECMPROB_DMIN + i);
        if (d - 1.0 >= half)
            break;
        if (d <= half)
            L += m[i];
        else
            L += m[i] * (half - (d - 1.0));     // partial band
    }
    if (half > (double)ECMPROB_DMAX)
        L += tail_lambda;   // (approximate: whole tail)

    return -expm1(-L) * T(c, ctx);
}

double ecmbayes_saving(const double *lambda, const double *W, const double *w,
    double digits, ecmbayes_time_fn T, void *ctx)
{
    int i;
    double m[ECMPROB_NDIG];
    double TD = T(digits, ctx);
    double num = 0.0, den = 0.0;

    post_mass(lambda, W, m);
    for (i = 0; i < ECMPROB_NDIG; i++)
    {
        double pi = m[i] * -expm1(-w[i]);
        double d = (double)(ECMPROB_DMIN + i) - 0.5;
        if (pi <= 0.0)
            continue;
        num += pi * (TD - cofactor_cost(m, lambda[ECMPROB_NDIG],
            digits - d, T, ctx));
        den += pi;
    }
    return (den > 0.0) ? num / den : 0.0;
}

void ecmbayes_option_init(ecmbayes_option_t *opt, const ecm_curve_param_t *cp,
    double sec_per_curve, double min_curves)
{
    opt->cp = *cp;
    opt->sec_per_curve = sec_per_curve;
    opt->min_curves = min_curves;
    memset(opt->w1, 0, sizeof(opt->w1));
    ecmprob_work_add(opt->w1, cp, 1.0);
}

int ecmbayes_decide(const double *lambda, const double *W, double digits,
    ecmbayes_time_fn T, void *ctx, const ecmbayes_option_t *opts, int nopt,
    double quantum, ecmbayes_decision_t *dec)
{
    int j, i;
    double m[ECMPROB_NDIG];
    double L;

    memset(dec, 0, sizeof(*dec));
    dec->best = -1;
    dec->stop = 1;
    dec->tsieve = T(digits, ctx);
    if (quantum <= 0.0)
        quantum = dec->tsieve / 1000.0;

    L = post_mass(lambda, W, m);
    dec->pcomposite_left = -expm1(-L);

    for (j = 0; j < nopt; j++)
    {
        double w[ECMPROB_NDIG];
        double n, cost, pf, sv, ben, r;

        if (opts[j].sec_per_curve <= 0.0)
            continue;
        n = quantum / opts[j].sec_per_curve;
        if (n < opts[j].min_curves)
            n = opts[j].min_curves;
        n = ceil(n);
        cost = n * opts[j].sec_per_curve;

        for (i = 0; i < ECMPROB_NDIG; i++)
            w[i] = n * opts[j].w1[i];

        pf = ecmbayes_pfind(lambda, W, w);
        sv = ecmbayes_saving(lambda, W, w, digits, T, ctx);
        ben = pf * sv;
        r = ben / cost;

        if ((dec->best < 0) || (r > dec->ratio))
        {
            dec->best = j;
            dec->curves = n;
            dec->cost = cost;
            dec->pfind = pf;
            dec->benefit = ben;
            dec->ratio = r;
        }
    }

    dec->stop = (dec->best < 0) || (dec->ratio < 1.0);
    return dec->stop;
}

#ifdef ECMBAYES_MAIN

/*
 * Test driver: simulate the policy on a model machine, assuming every ECM
 * step fails (the path whose length matters), and compare the expected total
 * time with the fixed-ratio plans.
 *
 * Model machine (all wall clock, all threads):
 *   T(D)        = a exp(b D) for GNFS (D >= xover), qa exp(qb D) for SIQS
 *   ECM curve   = k * (B1 / 1e6) * ((D + 10) / 160)^2 seconds
 * Defaults are a rough 6-thread desktop: C120 GNFS in 3960 s (measured,
 * docs/nfs_job_data/c120), doubling every 7 digits.  They are placeholders
 * until tune data exist.
 */

typedef struct
{
    double a, b;        // gnfs
    double qa, qb;      // siqs
    double xover;
} model_t;

static double model_T(double D, void *ctx)
{
    model_t *M = (model_t *)ctx;
    double g = M->a * exp(M->b * D);
    double q = M->qa * exp(M->qb * D);
    if (D < 1.0)
        return 0.0;
    return (D < M->xover) ? q : (g < q ? g : q);
}

static const double std_b1[] = { 11e3, 5e4, 25e4, 1e6, 3e6, 11e6, 43e6,
    11e7, 26e7, 85e7 };
#define NSTD ((int)(sizeof(std_b1) / sizeof(std_b1[0])))

typedef struct
{
    double ecm_time;    // ECM time on the all-fail path
    double tlevel;      // t-level reached
    double etotal;      // expected total time (ECM + sieve on N or cofactor)
    int steps;
} path_t;

// run a policy: rule 0 = Bayesian stop; rule 1 = stop at t-level >= target
// (steps still chosen by the Bayesian benefit/cost ranking)
static path_t run_policy(int rule, double target, double D, model_t *M,
    ecmbayes_option_t *opts, int nopt, double quantum, int verbose)
{
    double lambda[ECMBAYES_NLAM], W[ECMPROB_NDIG];
    double S = 1.0;     // P(no factor found yet)
    double TD = model_T(D, M);
    path_t P;
    int j;

    memset(&P, 0, sizeof(P));
    memset(W, 0, sizeof(W));
    ecmbayes_prior(lambda, ECMBAYES_PRIOR_RANDOM, D, 12.0, 0);

    for (;;)
    {
        ecmbayes_decision_t dec;
        double t = ecmprob_work_tlevel(W);
        int stop = ecmbayes_decide(lambda, W, D, model_T, M, opts, nopt,
            quantum, &dec);

        if (rule == 1)
            stop = (t >= target) || (dec.best < 0);

        if (stop || P.steps > 100000)
        {
            P.etotal += S * TD;
            P.tlevel = t;
            break;
        }

        // expected time: run the step (if still running), then on a find
        // the cofactor costs TD - saving
        P.etotal += S * dec.cost;
        P.etotal += S * dec.pfind * (TD - dec.benefit / dec.pfind);
        S *= (1.0 - dec.pfind);
        P.ecm_time += dec.cost;
        P.steps++;
        for (j = 0; j < ECMPROB_NDIG; j++)
            W[j] += dec.curves * opts[dec.best].w1[j];

        if (verbose > 1)
            printf("  step %d: B1=%.0f curves=%.0f cost=%.1f pfind=%.4g "
                "ratio=%.3f t=%.2f\n", P.steps, opts[dec.best].cp.B1,
                dec.curves, dec.cost, dec.pfind, dec.ratio,
                ecmprob_work_tlevel(W));
    }
    return P;
}

int main(int argc, char **argv)
{
    model_t M;
    double k = 0.1;             // s per curve at B1=1e6, C150, all threads
    double quantum = 0;         // default T/1000
    int d0 = 90, d1 = 200, dd = 10, verbose = 1, i, j;
    ecmbayes_option_t opts[NSTD];

    // GNFS: 3960 s at C120, x2 per 7 digits
    M.b = log(2.0) / 7.0;
    M.a = 3960.0 / exp(M.b * 120.0);
    // SIQS: x2 per 4 digits, crossover ~C95 (placeholder)
    M.qb = log(2.0) / 4.0;
    M.xover = 95.0;
    M.qa = (M.a * exp(M.b * M.xover)) / exp(M.qb * M.xover);

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-k") && i + 1 < argc) k = atof(argv[++i]);
        else if (!strcmp(argv[i], "-q") && i + 1 < argc) quantum = atof(argv[++i]);
        else if (!strcmp(argv[i], "-d") && i + 3 < argc)
        {
            d0 = atoi(argv[++i]); d1 = atoi(argv[++i]); dd = atoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "-tscale") && i + 1 < argc)
        {
            double s = atof(argv[++i]);
            M.a *= s; M.qa *= s;
        }
        else if (!strcmp(argv[i], "-v")) verbose++;
        else
        {
            fprintf(stderr, "usage: %s [-k sec_per_1M_curve_at_C150] "
                "[-q quantum_s] [-d from to step] [-tscale x] [-v]\n", argv[0]);
            return 1;
        }
    }

    ecmprob_init();

    printf("model: GNFS T(120) = %.0f s, x2 per %.1f digits; SIQS x2 per %.1f "
        "digits below C%.0f; ECM k = %.3g s/curve (B1=1e6, C150)\n",
        model_T(120, &M), log(2.0) / M.b, log(2.0) / M.qb, M.xover, k);
    printf("%5s %10s | %7s %9s %6s %10s | %7s %10s | %7s %10s | %10s\n",
        "D", "T(D) s", "bayes t", "ECM s", "ECM/T", "E[total]",
        "4/13 t", "E[total]", "1/3 t", "E[total]", "noECM");

    for (int D = d0; D <= d1; D += dd)
    {
        path_t pb, p413, p13;
        double TD = model_T(D, &M);
        double scale = ((D + 10.0) / 160.0) * ((D + 10.0) / 160.0);

        for (j = 0; j < NSTD; j++)
        {
            ecm_curve_param_t cp;
            ecmprob_gmpecm_stage2(&cp, std_b1[j], 0, 0,
                ECMPROB_FAM_EDW_IV, NULL, NULL);
            ecmbayes_option_init(&opts[j], &cp,
                k * (std_b1[j] / 1e6) * scale, 8);
        }

        pb = run_policy(0, 0, D, &M, opts, NSTD, quantum, verbose);
        p413 = run_policy(1, 4.0 / 13.0 * D, D, &M, opts, NSTD, quantum, 0);
        p13 = run_policy(1, D / 3.0, D, &M, opts, NSTD, quantum, 0);

        printf("%5d %10.4g | %7.2f %9.4g %6.3f %10.5g | %7.2f %10.5g | "
            "%7.2f %10.5g | %10.5g\n", D, TD, pb.tlevel, pb.ecm_time,
            pb.ecm_time / TD, pb.etotal, p413.tlevel, p413.etotal,
            p13.tlevel, p13.etotal, TD);
    }
    return 0;
}

#endif
