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
 * Bayesian ECM-vs-sieve stopping rule (factor/ecm_bayes.c, docs/ecm_bayes.md).
 *
 * Prime factors of N in each digit band d are modelled as independent
 * Poisson counts with mean lambda[d] (the prior).  ECM work W[d] (the
 * expected finds of a d-digit factor, ecm_obj.twork[] from ecm_prob.c)
 * thins the process: the undetected factors have mean lambda[d] exp(-W[d]).
 * Conditioning on N being composite (at least one factor at or below
 * sqrt(N) remains), a step that adds work w[d] finds a factor with
 * probability
 *
 *     P = (1 - exp(-L_s)) / (1 - exp(-L)),
 *     L   = sum_d lambda[d] exp(-W[d]),
 *     L_s = sum_d lambda[d] exp(-W[d]) (1 - exp(-w[d])).
 *
 * The value of a find of a d-digit factor is the sieve time it saves,
 * T(D) - V(D - d), with V the expected cost of the remaining cofactor.
 * The stopping problem is monotone (each failure lowers every future
 * success probability), so the one-step look-ahead rule is optimal: run
 * the step with the best benefit per second while its expected benefit
 * exceeds its cost, and otherwise start the sieve.
 *
 * Self-contained apart from ecm_prob.h, so it can be built and tested on
 * its own:
 *     gcc -O2 -DECMBAYES_MAIN -Iinclude factor/ecm_bayes.c \
 *         factor/gmp-ecm/ecm_prob.c -lm
 */

#ifndef ECM_BAYES_H
#define ECM_BAYES_H

#include "ecm_prob.h"

#ifdef __cplusplus
extern "C" {
#endif

// prior arrays: lambda[i] for the band of d = ECMPROB_DMIN + i digits
// (10^(d-1) <= p < 10^d, the band ecm_prob's work grid entry d stands for),
// i < ECMPROB_NDIG, plus lambda[ECMPROB_NDIG] = the mass above the grid
// (factors ECM never finds, but which count when conditioning on N being
// composite).  Work arrays W[], w[] are ecm_prob work grids (ECMPROB_NDIG).
#define ECMBAYES_NLAM   (ECMPROB_NDIG + 1)

// prior types
#define ECMBAYES_PRIOR_RANDOM   0   // random integer: lambda[d] = ln(d / (d-1))
#define ECMBAYES_PRIOR_HARD     1   // no factor expected below D/2 - margin

// expected sieve time (seconds, wall clock) for a composite of the given
// number of decimal digits: min(SIQS, GNFS) or the SNFS time for this input
typedef double (*ecmbayes_time_fn)(double digits, void *ctx);

// one candidate ECM step
typedef struct
{
    ecm_curve_param_t cp;   // curve parameters at this B1 (ecm_plan_curve_params)
    double sec_per_curve;   // wall-clock seconds per curve with all threads
    double min_curves;      // smallest sensible step (e.g. one vector batch)
    double w1[ECMPROB_NDIG];// work of one curve (ecmbayes_option_init)
} ecmbayes_option_t;

// fill an option (computes w1 from cp: call once per option, not per step)
void ecmbayes_option_init(ecmbayes_option_t *opt, const ecm_curve_param_t *cp,
    double sec_per_curve, double min_curves);

// the decision for one step
typedef struct
{
    int stop;               // 1: start the sieve; 0: run option 'best'
    int best;               // index into the options (-1 if none)
    double curves;          // curves to run at that option
    double cost;            // their expected time (s)
    double pfind;           // probability that they find a factor
    double benefit;         // expected saving (s) = pfind * E[saving | find]
    double ratio;           // benefit / cost of the best option
    double tsieve;          // T(D)
    double pcomposite_left; // P(a factor remains below sqrt N | failures so
                            //   far), unconditioned on compositeness; for
                            //   reporting only
} ecmbayes_decision_t;

// fill lambda[ECMBAYES_NLAM] (indexed d - ECMPROB_DMIN) with the prior for a
// composite of 'digits' digits with no prime factor below 'dmin' digits
// (trial division / rho / known factors; clamped to >= ECMPROB_DMIN - 1).  'margin' is used by the HARD
// prior only (digits below D/2 with no factor expected).
void ecmbayes_prior(double *lambda, int type, double digits, double dmin,
    double margin);

// probability that adding work w[] (expected finds per band) finds a factor,
// given prior lambda[], completed work W[] and that N is composite with its
// smallest factor in the grid
double ecmbayes_pfind(const double *lambda, const double *W, const double *w);

// expected saving (s) if a factor is found by work w[]: averages
// T(D) - V(D - d) over the bands the step reaches, V(c) the expected cost of
// a c-digit cofactor (P(composite) * T(c))
double ecmbayes_saving(const double *lambda, const double *W, const double *w,
    double digits, ecmbayes_time_fn T, void *ctx);

// the one-step look-ahead decision.  Each option is evaluated for a step of
// max(min_curves, quantum / sec_per_curve) curves; quantum <= 0 uses
// T(D) / 1000.  Returns dec->stop.
int ecmbayes_decide(const double *lambda, const double *W, double digits,
    ecmbayes_time_fn T, void *ctx, const ecmbayes_option_t *opts, int nopt,
    double quantum, ecmbayes_decision_t *dec);

#ifdef __cplusplus
}
#endif

#endif
