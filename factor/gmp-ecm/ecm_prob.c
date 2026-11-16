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
 * ========== Note: createed by Claude with prompting by bbuhrow ==========
 * ECM success probabilities and t-level accounting.  See include/ecm_prob.h
 * for the interface and docs/ecm_prob.md for the model, its validation
 * against GMP-ECM and the measured extra smoothness of the curve families.
 *
 * An independent implementation of the model GMP-ECM prints with -v
 * (Kruppa, "Speeding up integer multiplication and factorization", PhD
 * thesis 2010, ch. 5): the order of the group mod p is taken to be as
 * likely smooth as a random integer near p / exp(delta), and
 *
 *   P = rho_loc(u) + sum over primes q in (B1, B2] of rho_loc(u_q) / q
 *       + sum over primes q > B2 of rho_loc(u_q) / q * P(Brent-Suyama hit)
 *
 * with u = ln(p e^-delta) / ln B1, u_q = ln(p e^-delta / q) / ln B1 and
 * the "local" Dickman function rho_loc(u) = rho(u) - gamma rho(u-1) / ln x
 * for numbers near x (thesis eq. 5.6).  The prime sums are exact below
 * ECMPROB_QSUM and integrals (dq / (q ln q)) above.  Like GMP-ECM the
 * integrals take the local correction at x = p e^-delta rather than at the
 * cofactor x / q (which gives ~1.7% fewer finds; t-levels are GMP-ECM's
 * convention, so we keep its choice).  Unlike GMP-ECM, for B1 < ECMPROB_QSUM
 * the integral covers exactly (ECMPROB_QSUM, B2] (GMP-ECM 7 mixes log
 * bases there and misses part of stage 2, so its B1 = 2e3 and 11e3 tables
 * need up to 2x too many curves).
 *
 * Standalone test driver: build with -DECMPROB_MAIN, see main() below.
 */

#include "ecm_prob.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define RHO_INV     256             // table step 1/256 in u
#define RHO_UMAX    80              // rho(u) = 0 beyond (rho(80) ~ 1e-150)
#define ECMPROB_QSUM 20000          // exact prime sum below this q
#define EULER_GAMMA 0.57721566490153286061

static double *rho_tab = NULL;
static int rho_n = 0;
static uint32_t *qprimes = NULL;    // primes < ECMPROB_QSUM
static int nqprimes = 0;

// ---------------------------------------------------------------------------
// Dickman rho: rho = 1 on [0,1], 1 - ln u on [1,2], then from the identity
// u rho(u) = int_{u-1}^{u} rho(t) dt by Simpson's rule over the RHO_INV
// (even) table steps of [u-1, u], solved for the rho(u) end point.  Every
// term is positive, so the relative accuracy holds for tiny rho (a
// recurrence on rho' = -rho(u-1)/u loses it by cancellation).
// ---------------------------------------------------------------------------
static void rho_build(void)
{
    const double h = 1.0 / RHO_INV;
    int n = RHO_UMAX * RHO_INV + 1, i, j, R = RHO_INV;
    double *r = (double *)malloc((size_t)n * sizeof(double));

    for (i = 0; i <= R; i++)
        r[i] = 1.0;
    for (; i <= 2 * R; i++)
        r[i] = 1.0 - log(i * h);
    for (; i < n; i++)
    {
        // Simpson weights 1 4 2 4 ... 2 4 1 on r[i-R .. i]
        double s = r[i - R];
        for (j = 1; j < R; j++)
            s += ((j & 1) ? 4.0 : 2.0) * r[i - R + j];
        r[i] = (h / 3.0) * s / (i * h - h / 3.0);
    }
    rho_tab = r;
    rho_n = n;
}

static void primes_build(void)
{
    char *c = (char *)calloc(ECMPROB_QSUM, 1);
    int i, j, k = 0;
    for (i = 2; i < ECMPROB_QSUM; i++)
        if (!c[i])
            for (j = 2 * i; j < ECMPROB_QSUM; j += i)
                c[j] = 1;
    for (i = 2; i < ECMPROB_QSUM; i++)
        k += !c[i];
    qprimes = (uint32_t *)malloc((size_t)k * sizeof(uint32_t));
    for (i = 2, k = 0; i < ECMPROB_QSUM; i++)
        if (!c[i])
            qprimes[k++] = (uint32_t)i;
    nqprimes = k;
    free(c);
}

void ecmprob_init(void)
{
    if (rho_tab == NULL)
        rho_build();
    if (qprimes == NULL)
        primes_build();
}

double ecmprob_rho(double u)
{
    double x, f;
    int i;
    if (u <= 1.0)
        return (u < 0) ? 0.0 : 1.0;
    if (u <= 2.0)
        return 1.0 - log(u);
    if (rho_tab == NULL)
        ecmprob_init();
    x = u * RHO_INV;
    i = (int)x;
    if (i + 1 >= rho_n)
        return 0.0;
    f = x - i;
    return rho_tab[i] * (1.0 - f) + rho_tab[i + 1] * f;
}

// local rho for numbers near x (ln x = lnx): rho(u) - gamma rho(u-1) / ln x
static double rho_loc(double u, double lnx)
{
    double v;
    if (u <= 1.0)
        return (u < 0) ? 0.0 : 1.0;
    v = ecmprob_rho(u) - EULER_GAMMA * ecmprob_rho(u - 1.0) / lnx;
    return v > 0 ? v : 0.0;
}

// ---------------------------------------------------------------------------
// curve families
// ---------------------------------------------------------------------------

// GMP-ECM's ECM_EXTRA_SMOOTHNESS for Suyama's curves (2^3.323 3^1.687 and
// the primes >= 5, against a random integer)
#define DELTA_SUYAMA    3.134

double ecmprob_delta(int family)
{
    switch (family)
    {
    case ECMPROB_FAM_SUYAMA:
    case ECMPROB_FAM_P2:
        return DELTA_SUYAMA;
    case ECMPROB_FAM_P1:
        // GMP-ECM ecm.c EXTRA_SMOOTHNESS_SQUARE (Bouvier, Barbulescu):
        // average torsion of Suyama times 2^(1/3) / (3 * 3^(1/128))
        return DELTA_SUYAMA + log(0.416384512396064);
    case ECMPROB_FAM_P3:
        // EXTRA_SMOOTHNESS_32BITS_D: 1 / (3 * 3^(1/128))
        return DELTA_SUYAMA + log(0.330484606500389);
    case ECMPROB_FAM_EDW_IV:
        // measured (arith/edwards/fam_delta.py, docs/ecm_prob.md):
        // E[v2(#E)] = 5.335, E[v3(#E)] = 0.676 against Suyama's 3.333 and
        // 1.684 on the same primes; delta(iv) - delta(Suyama) = +0.279
        // +- 0.020 over 2 x 40k paired group orders (32- and 40-bit p)
        return DELTA_SUYAMA + 0.279;
    default:
        return DELTA_SUYAMA;
    }
}

const char *ecmprob_family_name(int family)
{
    switch (family)
    {
    case ECMPROB_FAM_SUYAMA: return "Suyama (param 0)";
    case ECMPROB_FAM_P1: return "Montgomery param 1";
    case ECMPROB_FAM_P2: return "Montgomery param 2";
    case ECMPROB_FAM_P3: return "Montgomery param 3";
    case ECMPROB_FAM_EDW_IV: return "twisted Edwards family (iv)";
    default: return "unknown";
    }
}

// ---------------------------------------------------------------------------
// Brent-Suyama: a prime q > B2 is caught when some stage 2 pair has
// f(a) = +-f(b) mod q with a != +-b; for a q in a residue class with r such
// non-trivial roots the expected number of hits is nr r / q.  The classes
// (r, weight) per polynomial:
// ---------------------------------------------------------------------------
#define BS_MAXCLS 64
typedef struct { int n; double r[BS_MAXCLS], w[BS_MAXCLS]; } bs_cls_t;

static int gcd_i(int a, int b)
{
    while (b) { int t = a % b; a = b; b = t; }
    return a < 0 ? -a : a;
}

// AVX-ECM FFT stage 2 (avxecm_stage2fft.c s2_rhist, counted over random
// primes with arith/vecpoly/ecm_tests/bs_model.py); the weights don't sum
// to 1: the rest has r = 0
typedef struct { int S, cnt; double r[6], w[6]; } bs_hist_t;
static const bs_hist_t bs_hist_avx[] = {
    {  2, 1, { 2 }, { 0.500 } },
    {  3, 1, { 4 }, { 0.467 } },
    {  6, 3, { 2, 4, 10 }, { 0.273, 0.240, 0.227 } },
    { -3, 1, { 4 }, { 0.485 } },
    { -6, 3, { 2, 4, 10 }, { 0.251, 0.262, 0.223 } },
    { -12, 5, { 2, 4, 6, 10, 22 }, { 0.147, 0.262, 0.104, 0.102, 0.122 } },
};

// GMP-ECM's classes (rho.c): with m = 2|S|, x^S: i in [1, m) coprime to m,
// r = gcd(i-1, m) - 2; Dickson: i in [1, m/2] coprime to m,
// r = (gcd(i-1, m) + gcd(i+1, m) - 4) / 2; equal weights
static void bs_classes(int S, int model, bs_cls_t *c)
{
    int i, m, k;
    c->n = 0;
    if ((S >= -1) && (S <= 1))
        return;
    if (model == ECMPROB_BS_AVXECM)
    {
        for (i = 0; i < (int)(sizeof(bs_hist_avx) / sizeof(bs_hist_avx[0])); i++)
        {
            if (bs_hist_avx[i].S == S)
            {
                for (k = 0; k < bs_hist_avx[i].cnt; k++)
                {
                    c->r[k] = bs_hist_avx[i].r[k];
                    c->w[k] = bs_hist_avx[i].w[k];
                }
                c->n = bs_hist_avx[i].cnt;
                return;
            }
        }
        // not tabulated: GMP-ECM's classes below
    }
    m = 2 * (S > 0 ? S : -S);
    if (m > 2 * 120)
        m = 2 * 120;
    k = 0;
    if (S > 0)
    {
        for (i = 1; i < m; i++)
            if (gcd_i(i, m) == 1)
                k++;
        for (i = 1; i < m; i++)
            if (gcd_i(i, m) == 1)
            {
                double r = (double)(gcd_i(i - 1, m) - 2);
                if (c->n < BS_MAXCLS)
                {
                    c->r[c->n] = r;
                    c->w[c->n] = 1.0 / k;
                    c->n++;
                }
            }
    }
    else
    {
        for (i = 1; i <= m / 2; i++)
            if (gcd_i(i, m) == 1)
                k++;
        c->n = 0;
        for (i = 1; i <= m / 2; i++)
            if (gcd_i(i, m) == 1)
            {
                double r = (double)(gcd_i(i - 1, m) + gcd_i(i + 1, m) - 4) / 2.0;
                if (c->n < BS_MAXCLS)
                {
                    c->r[c->n] = r;
                    c->w[c->n] = 1.0 / k;
                    c->n++;
                }
            }
    }
}

// ---------------------------------------------------------------------------
// probability
// ---------------------------------------------------------------------------

// Simpson intervals (even) for a log-q range of width wl, in units of
// a = ln B1: about RHO_INV / 4 per unit of u
static int simpson_n(double wl, double a)
{
    int n = (int)ceil(wl / a * (RHO_INV / 4)) + 16;
    return n + (n & 1);
}

double ecmprob_prob(const ecm_curve_param_t *c, double p)
{
    double B1 = c->B1, B2 = c->B2;
    double L, a, P, lo, hi, h;
    int i, n;

    if ((B1 < 2.0) || (p < 2.0))
        return 0.0;
    if (rho_tab == NULL)
        ecmprob_init();

    L = log(p) - ecmprob_delta(c->family);      // ln of the effective size
    a = log(B1);
    if (L <= a)
        return 1.0;

    // stage 1: the whole effective order is B1-smooth
    P = rho_loc(L / a, L);

    // stage 2: B1-smooth but for one prime q in (B1, B2]
    if (B2 > B1)
    {
        for (i = 0; i < nqprimes; i++)
        {
            double q = (double)qprimes[i], l;
            if (q <= B1)
                continue;
            if ((q > B2) || (log(q) >= L))
                break;
            l = log(q);
            P += rho_loc((L - l) / a, L - l) / q;
        }
        lo = log(B1 > ECMPROB_QSUM ? B1 : (double)ECMPROB_QSUM);
        hi = log(B2);
        if (hi > L) hi = L;
        if (hi > lo)
        {
            double s = 0;
            n = simpson_n(hi - lo, a);
            h = (hi - lo) / n;
            for (i = 0; i <= n; i++)
            {
                double l = lo + i * h;
                double wt = (i == 0 || i == n) ? 1.0 : ((i & 1) ? 4.0 : 2.0);
                s += wt * rho_loc((L - l) / a, L) / l;
            }
            P += s * h / 3.0;
        }
    }

    // Brent-Suyama: one prime q in (B2, effective size] caught by a pair
    if ((c->S < -1 || c->S > 1) && (c->nr > 0) && (B2 > B1) && (L > log(B2)))
    {
        bs_cls_t cls;
        double s = 0;
        int j;
        memset(&cls, 0, sizeof(cls));
        bs_classes(c->S, c->bsmodel, &cls);
        if (cls.n > 0)
        {
            lo = log(B2);
            hi = L;
            n = simpson_n(hi - lo, a);
            h = (hi - lo) / n;
            for (i = 0; i <= n; i++)
            {
                double l = lo + i * h, E = c->nr * exp(-l), hit = 0;
                double wt = (i == 0 || i == n) ? 1.0 : ((i & 1) ? 4.0 : 2.0);
                for (j = 0; j < cls.n; j++)
                    hit += cls.w[j] * -expm1(-E * cls.r[j]);
                s += wt * rho_loc((L - l) / a, L) / l * hit;
            }
            P += s * h / 3.0;
        }
    }

    return P < 1.0 ? P : 1.0;
}

double ecmprob_prob_digits(const ecm_curve_param_t *c, double digits)
{
    return ecmprob_prob(c, pow(10.0, digits - 0.5));
}

double ecmprob_curves(const ecm_curve_param_t *c, double digits)
{
    double P = ecmprob_prob_digits(c, digits);
    return (P > 1e-300) ? 1.0 / P : 1e300;
}

// ---------------------------------------------------------------------------
// GMP-ECM's stage 2 parameters (its bestD with the NTT, power-of-2 dF, and
// default k; the Dickson degree by stage 2 length).  Values checked against
// ecm -v for B1 = 2e3 .. 2.9e9 (docs/ecm_prob.md).
// ---------------------------------------------------------------------------
static const uint32_t gmp_d_po2[] = { 12, 30, 60, 120, 240, 510, 1020, 2310,
    4620, 9240, 19110, 39270, 79170, 158340, 324870, 690690, 1345890, 2852850,
    5705700, 11741730, 23130030, 48498450, 96996900 };

static uint64_t eulerphi_u(uint64_t n)
{
    uint64_t r = n, q;
    for (q = 2; q * q <= n; q++)
    {
        if (n % q == 0)
        {
            while (n % q == 0) n /= q;
            r -= r / q;
        }
    }
    if (n > 1) r -= r / n;
    return r;
}

// count of i in [i0, x] not divisible by d2 (d2 > 1)
static double count_nondiv(double i0, double x, double d2)
{
    return (x - i0 + 1) - (floor(x / d2) - floor((i0 - 1) / d2));
}

void ecmprob_gmpecm_stage2(ecm_curve_param_t *c, double B1, double B2, int S,
    int family, uint32_t *pdF, uint32_t *pk)
{
    int i, nd = (int)(sizeof(gmp_d_po2) / sizeof(gmp_d_po2[0]));
    double d1 = 0, d2 = 1, dF = 1, k = 0, i0 = 0, i1, len;

    c->B1 = B1;
    c->family = family;
    c->bsmodel = ECMPROB_BS_GMPECM;
    c->nr = 0;
    c->S = 0;
    if (B2 <= 0)
        B2 = pow(11.0 / 6.0 * B1, 1.43);
    if (B2 <= B1)
    {
        c->B2 = B1;
        if (pdF) *pdF = 0;
        if (pk) *pk = 0;
        return;
    }

    for (i = 0; i < nd; i++)
    {
        double phid, j, t;
        int lg = 0;
        d1 = (double)gmp_d_po2[i];
        phid = (double)(eulerphi_u(gmp_d_po2[i]) / 2);
        while ((double)(1u << lg) < phid) lg++;
        dF = (double)(1u << lg);
        // smallest prime 5 <= d2 < 25 (not 3) not dividing d1
        for (d2 = 5; d2 < 25; d2 += 2)
        {
            if (((int)d2 % 3) == 0)
                continue;
            if (fmod(d1, d2) != 0)
                break;
        }
        if ((d2 >= 25) || (d2 - 1 > dF))
            d2 = 1;
        i0 = ceil((B1 - (d1 - 1) * d2) / d1);
        t = floor((B2 + (d1 - 1) * d2) / d1);
        j = t - i0 + 1;
        if (d2 > 1)
            j -= floor(t / d2) - floor(i0 / d2);
        k = ceil(j / dF);
        if (k <= 6)
            break;
    }
    if (k < 1) k = 1;

    // the last root: k dF values of i >= i0 coprime to d2
    if (d2 > 1)
    {
        double want = k * dF;
        i1 = i0 - 1 + floor(want * d2 / (d2 - 1));
        while (count_nondiv(i0, i1, d2) < want) i1 += 1;
        while ((i1 > i0) && (count_nondiv(i0, i1 - 1, d2) >= want)) i1 -= 1;
    }
    else
    {
        i1 = i0 + k * dF - 1;
    }
    c->B2 = d1 * (i1 - d2 + 1) + 5 * d2 - 1;
    if (c->B2 < B2)
        c->B2 = B2;
    c->nr = dF * dF * k;

    if (S != 0)
        c->S = S;
    else
    {
        len = c->B2 - B1;
        if (len < 1e7) c->S = 1;
        else if (len < 1e8) c->S = 2;
        else if (len < 1e9) c->S = -3;
        else if (len < 1e10) c->S = -6;
        else if (len < 3e11) c->S = -12;
        else c->S = -30;
    }
    if (pdF) *pdF = (uint32_t)dF;
    if (pk) *pk = (uint32_t)k;
}

// ---------------------------------------------------------------------------
// t-level accounting on the digit grid
// ---------------------------------------------------------------------------
void ecmprob_work_add(double *w, const ecm_curve_param_t *c, double curves)
{
    int i;
    if (curves <= 0)
        return;
    for (i = 0; i < ECMPROB_NDIG; i++)
        w[i] += curves * ecmprob_prob_digits(c, (double)(ECMPROB_DMIN + i));
}

double ecmprob_work_tlevel(const double *w)
{
    int i;
    double l0, l1;
    if (w[0] < 1.0)
        return 0.0;
    for (i = 1; i < ECMPROB_NDIG; i++)
        if (w[i] < 1.0)
            break;
    if (i == ECMPROB_NDIG)
        return (double)ECMPROB_DMAX;
    // ln w is close to linear in the digits between grid points
    l0 = log(w[i - 1]);
    l1 = log(w[i] > 1e-300 ? w[i] : 1e-300);
    return (double)(ECMPROB_DMIN + i - 1) + l0 / (l0 - l1);
}

double ecmprob_work_at(const double *w, double digits)
{
    double x = digits - ECMPROB_DMIN, f, a, b;
    int i;
    if (x <= 0)
        return w[0];
    if (x >= ECMPROB_NDIG - 1)
        return w[ECMPROB_NDIG - 1];
    i = (int)x;
    f = x - i;
    a = w[i];
    b = w[i + 1];
    if ((a <= 0) || (b <= 0))
        return a * (1 - f) + b * f;
    return exp(log(a) * (1 - f) + log(b) * f);
}

uint32_t ecmprob_work_curves_to(const double *w, const ecm_curve_param_t *c,
    double target)
{
    double pv[ECMPROB_NDIG], t[ECMPROB_NDIG];
    uint64_t lo = 0, hi = 1, mid;
    int i;

    if (ecmprob_work_tlevel(w) >= target)
        return 0;
    for (i = 0; i < ECMPROB_NDIG; i++)
        pv[i] = ecmprob_prob_digits(c, (double)(ECMPROB_DMIN + i));

    // smallest n with tlevel(w + n pv) >= target: grow, then bisect
    while (1)
    {
        for (i = 0; i < ECMPROB_NDIG; i++)
            t[i] = w[i] + (double)hi * pv[i];
        if (ecmprob_work_tlevel(t) >= target)
            break;
        lo = hi;
        hi *= 2;
        if (hi > 0xffffffffull)
            return 0xffffffffu;
    }
    while (hi - lo > 1)
    {
        mid = lo + (hi - lo) / 2;
        for (i = 0; i < ECMPROB_NDIG; i++)
            t[i] = w[i] + (double)mid * pv[i];
        if (ecmprob_work_tlevel(t) >= target)
            hi = mid;
        else
            lo = mid;
    }
    return (uint32_t)hi;
}

// ---------------------------------------------------------------------------
// standalone test driver
// ---------------------------------------------------------------------------
#ifdef ECMPROB_MAIN

static const double std_b1[] = { 2e3, 11e3, 5e4, 25e4, 1e6, 3e6, 11e6, 43e6,
    11e7, 26e7, 85e7, 29e8 };
#define NSTD ((int)(sizeof(std_b1) / sizeof(std_b1[0])))

static int fam_arg(const char *s)
{
    if (!strcmp(s, "suyama") || !strcmp(s, "p0")) return ECMPROB_FAM_SUYAMA;
    if (!strcmp(s, "p1")) return ECMPROB_FAM_P1;
    if (!strcmp(s, "p2")) return ECMPROB_FAM_P2;
    if (!strcmp(s, "p3")) return ECMPROB_FAM_P3;
    if (!strcmp(s, "iv") || !strcmp(s, "edw")) return ECMPROB_FAM_EDW_IV;
    fprintf(stderr, "family: suyama|p0|p1|p2|p3|iv\n");
    exit(1);
}

static void print_curves_row(const ecm_curve_param_t *c, int d0, int d1, int dd)
{
    int d;
    for (d = d0; d <= d1; d += dd)
    {
        double n = ecmprob_curves(c, d);
        if (n < 1e7) printf("%.0f\t", floor(n + 0.5));
        else if (n < 1e299) printf("%.2g\t", n);
        else printf("Inf\t");
    }
    printf("\n");
}

static void usage(void)
{
    printf("ecm_prob test driver\n"
        "  rho                         Dickman rho checks\n"
        "  gmp fam B1 [B2]             GMP-ECM stage 2 params and its -v table (35..80)\n"
        "  avx fam B1 B2 [S nr]        expected curves, AVX-ECM stage 2 (15..80)\n"
        "  std fam gmp|avx100|edwfft   table at the standard B1 levels (rows: digits)\n"
        "  prob fam B1 B2 S nr bsmodel digits\n"
        "  fam: suyama|p0|p1|p2|p3|iv\n");
}

int main(int argc, char **argv)
{
    ecm_curve_param_t c;
    ecmprob_init();
    if (argc < 2) { usage(); return 0; }

    if (!strcmp(argv[1], "rho"))
    {
        // reference values (Marsaglia-Zaman-Marsaglia / tables)
        static const double ref[][2] = { { 2.5, 0.130319035 }, { 3, 0.0486083882911 },
            { 4, 0.00491092564776 }, { 5, 0.000354724700 }, { 6, 1.96496963539e-5 },
            { 8, 3.23206930422e-8 }, { 10, 2.77017183772596e-11 },
            { 20, 2.46178622e-29 } };
        int i;
        for (i = 0; i < (int)(sizeof(ref) / sizeof(ref[0])); i++)
            printf("rho(%g) = %.12g  ref %.12g  rel err %.2e\n", ref[i][0],
                ecmprob_rho(ref[i][0]), ref[i][1],
                ecmprob_rho(ref[i][0]) / ref[i][1] - 1.0);
        return 0;
    }
    if (!strcmp(argv[1], "gmp") && argc >= 4)
    {
        uint32_t dF, k;
        double B2 = (argc > 4) ? atof(argv[4]) : 0;
        ecmprob_gmpecm_stage2(&c, atof(argv[3]), B2, 0, fam_arg(argv[2]), &dF, &k);
        printf("B1=%.0f, B2=%.0f, S=%d, dF=%u, k=%u, nr=%.0f\n", c.B1, c.B2, c.S, dF, k, c.nr);
        print_curves_row(&c, 35, 80, 5);
        return 0;
    }
    if (!strcmp(argv[1], "avx") && argc >= 5)
    {
        c.family = fam_arg(argv[2]);
        c.B1 = atof(argv[3]);
        c.B2 = atof(argv[4]);
        c.S = (argc > 5) ? atoi(argv[5]) : 0;
        c.nr = (argc > 6) ? atof(argv[6]) : 0;
        c.bsmodel = ECMPROB_BS_AVXECM;
        print_curves_row(&c, 15, 80, 5);
        return 0;
    }
    if (!strcmp(argv[1], "prob") && argc >= 9)
    {
        c.family = fam_arg(argv[2]);
        c.B1 = atof(argv[3]);
        c.B2 = atof(argv[4]);
        c.S = atoi(argv[5]);
        c.nr = atof(argv[6]);
        c.bsmodel = atoi(argv[7]);
        printf("%.10g\n", ecmprob_prob_digits(&c, atof(argv[8])));
        return 0;
    }
    if (!strcmp(argv[1], "std") && argc >= 4)
    {
        int fam = fam_arg(argv[2]), i, d;
        ecm_curve_param_t cs[NSTD];
        for (i = 0; i < NSTD; i++)
        {
            if (!strcmp(argv[3], "gmp"))
                ecmprob_gmpecm_stage2(&cs[i], std_b1[i], 0, 0, fam, NULL, NULL);
            else
            {
                cs[i].family = fam;
                cs[i].B1 = std_b1[i];
                cs[i].B2 = 100 * std_b1[i];
                cs[i].S = 0;
                cs[i].nr = 0;
                cs[i].bsmodel = ECMPROB_BS_AVXECM;
                if (!strcmp(argv[3], "edwfft"))
                    cs[i].B2 = pow(11.0 / 6.0 * std_b1[i], 1.43);
            }
        }
        printf("/*        ");
        for (i = 0; i < NSTD; i++) printf("%-9.3g", std_b1[i]);
        printf("*/\n");
        for (d = 15; d <= 70; d += 5)
        {
            printf("/*t%d*/ {", d);
            for (i = 0; i < NSTD; i++)
            {
                double n = ecmprob_curves(&cs[i], d);
                if (n < 1e7) printf("%.0f%s", floor(n + 0.5), i < NSTD - 1 ? ", " : "");
                else printf("%.2g%s", n > 9e99 ? 9e99 : n, i < NSTD - 1 ? ", " : "");
            }
            printf("},\n");
        }
        return 0;
    }
    usage();
    return 1;
}
#endif
