/*
Copyright (c) 2019, Ben Buhrow
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
 * ========== Note: createed by Claude with prompting by bbuhrow ==========
 * ECM success probabilities and t-level accounting (factor/gmp-ecm/ecm_prob.c,
 * docs/ecm_prob.md).
 *
 * The model is GMP-ECM's (Kruppa's thesis, ch. 5): the group order of a
 * curve mod p behaves like a random integer of size p / exp(delta), delta
 * the extra smoothness of the curve family; a curve finds p when that
 * number is B1-smooth, or B1-smooth but for one prime in (B1, B2], or (with
 * a Brent-Suyama polynomial) one prime above B2 that a stage 2 pair happens
 * to catch.  This file is self-contained (no yafu headers) so the model can
 * be built and tested on its own:
 *     gcc -O2 -DECMPROB_MAIN -Iinclude factor/gmp-ecm/ecm_prob.c -lm
 */

#ifndef ECM_PROB_H
#define ECM_PROB_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// curve families (their delta: ecmprob_delta)
#define ECMPROB_FAM_SUYAMA  0   // Suyama sigma: GMP-ECM param 0, -avxecm_curve 0
#define ECMPROB_FAM_P1      1   // GMP-ECM param 1 (A = 4d - 2, d a square)
#define ECMPROB_FAM_P2      2   // GMP-ECM param 2
#define ECMPROB_FAM_P3      3   // GMP-ECM param 3 (A = 4d - 2, d 32 bits; GPU)
#define ECMPROB_FAM_EDW_IV  4   // twisted Edwards a = -1, Gelin-Kleinjung-Lenstra
                                // family (iv): -avxecm_curve 1 and 2
#define ECMPROB_NFAM        5

// Brent-Suyama catch model
#define ECMPROB_BS_GMPECM   0   // GMP-ECM's root classes (rho.c)
#define ECMPROB_BS_AVXECM   1   // AVX-ECM FFT stage 2 (measured r histograms,
                                // arith/vecpoly/ecm_tests/bs_model.py)

// the parameters of one kind of curve, as run
typedef struct
{
    double B1;          // stage 1 bound
    double B2;          // largest prime stage 2 covers (<= B1: no stage 2)
    double nr;          // stage 2 pairs (baby x giant steps), for Brent-Suyama
    int S;              // Brent-Suyama polynomial: |S| <= 1 none,
                        //   S > 1 x^S, S < -1 Dickson D_|S|
    int bsmodel;        // ECMPROB_BS_*
    int family;         // ECMPROB_FAM_*
} ecm_curve_param_t;

// t-level work grid: expected finds of a d-digit factor for
// d = ECMPROB_DMIN .. ECMPROB_DMAX
#define ECMPROB_DMIN        10
#define ECMPROB_DMAX        90
#define ECMPROB_NDIG        (ECMPROB_DMAX - ECMPROB_DMIN + 1)

// build the Dickman rho and prime tables (idempotent; call once from a
// single thread before the other functions run in parallel)
void ecmprob_init(void);

// extra smoothness (natural log) of a curve family
double ecmprob_delta(int family);

// probability that one curve finds a prime factor near p
double ecmprob_prob(const ecm_curve_param_t *c, double p);

// the same for a factor of the given (fractional) number of decimal digits,
// taken as p = 10^(digits - 0.5) like GMP-ECM
double ecmprob_prob_digits(const ecm_curve_param_t *c, double digits);

// expected curves to find a factor of that many digits (1 / probability;
// 1e300 if the probability is 0)
double ecmprob_curves(const ecm_curve_param_t *c, double digits);

// Dickman rho(u)
double ecmprob_rho(double u);

// GMP-ECM's stage 2 for B1 and requested B2 (B2 <= 0: its default
// (11/6 B1)^1.43): fills c->B2 (the covered bound), c->nr = dF^2 k and
// c->S (its default Dickson degree, or S if S != 0), bsmodel GMP-ECM.
// Assumes the NTT stage 2 (power-of-2 dF), GMP-ECM's default for inputs
// below several thousand bits.  Returns dF and k through the pointers if
// they are not NULL.
void ecmprob_gmpecm_stage2(ecm_curve_param_t *c, double B1, double B2, int S,
    int family, uint32_t *dF, uint32_t *k);

// t-level accounting on the work grid w[ECMPROB_NDIG]:
// add curves of kind c
void ecmprob_work_add(double *w, const ecm_curve_param_t *c, double curves);
// the t-level: the digit size at which the expected finds are 1
// (log-linear interpolation in the grid; 0 if under ECMPROB_DMIN)
double ecmprob_work_tlevel(const double *w);
// expected finds of a factor of that many digits (interpolated)
double ecmprob_work_at(const double *w, double digits);
// curves of kind c that bring the t-level to target (0 if already there;
// UINT32_MAX if it can't be reached)
uint32_t ecmprob_work_curves_to(const double *w, const ecm_curve_param_t *c,
    double target);

// short description of a family / curve kind ("Edwards (iv)", ...)
const char *ecmprob_family_name(int family);

#ifdef __cplusplus
}
#endif

#endif
