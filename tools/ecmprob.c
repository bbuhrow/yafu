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
 * ecmprob: t-level and missed-factor probabilities for a set of ECM curves,
 * using yafu's probability model (factor/gmp-ecm/ecm_prob.c,
 * docs/ecm_prob.md).  Curve strings follow the t-level utility
 * (github.com/brubsby/t-level) with yafu's curve families added:
 *
 *     [done/]curves@B1[,B2][,param][,S=s][,nr=n][,s2=gmp|avx|pair|none]
 *         [,mem=MB][,threads=t][,bits=b | digits=d]
 *
 * `make ecmprob` builds it against libyecm with the build's ISA flags, so
 * AVX-ECM's FFT stage 2 is planned exactly as yafu plans it
 * (avxecm_s2fft_plan: automatic B2 and Brent-Suyama polynomial from the
 * input size, threads and memory limit, with this ISA's cost model).
 * Without AVX2 / AVX-512 (or with -DECMPROB_NO_PLANNER) it needs only the
 * model, and the FFT stage 2 falls back to the fixed B2 (11/6 B1)^1.43:
 *     gcc -O2 -Iinclude tools/ecmprob.c factor/gmp-ecm/ecm_prob.c -lm -o ecmprob
 */

#include "ecm_prob.h"
#if (defined(USE_AVX512F) || defined(USE_AVX2)) && !defined(ECMPROB_NO_PLANNER)
#define ECMPROB_PLANNER 1
#include "../factor/avx-ecm/avxecm_stage2fft.h"
#endif
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <io.h>
#define ISATTY_STDIN() _isatty(_fileno(stdin))
#else
#include <unistd.h>
#define ISATTY_STDIN() isatty(fileno(stdin))
#endif

#define S2_DEFAULT  0       // the family's usual stage 2
#define S2_GMP      1       // GMP-ECM's stage 2 (NTT, default k, Dickson poly)
#define S2_AVX      2       // AVX-ECM FFT stage 2
#define S2_NONE     3       // stage 1 only
#define S2_PAIR     4       // AVX-ECM pairing stage 2 (-avxecm_stg2 0)

#define S_UNSET     INT_MIN

// what the AVX-ECM FFT stage 2 plan depends on (yafu's defaults:
// -avxecm_stg2mem 256, -threads 1)
typedef struct
{
    double mem_mb;          // memory limit per thread, MB
    int threads;
    int bits;               // input size
} plan_env_t;

typedef struct
{
    char text[256];         // the curve string as given
    double done, total;     // curves with stage 2, curves in all
    ecm_curve_param_t c;    // the curves with stage 2
    ecm_curve_param_t c1;   // the stage-1-only remainder (total - done)
    int s2;                 // S2_* as resolved
    uint32_t dF, k;         // GMP-ECM stage 2 shape (S2_GMP)
    int nr_est;             // nr estimated, not given
    int planned;            // S2_AVX: 1 planned by avxecm_s2fft_plan,
                            //   0 fixed fallback
    int auto_b2;            // S2_AVX: B2 chosen by the planner
    plan_env_t env;
} spec_t;

static const double std_b1[] = { 2e3, 11e3, 5e4, 25e4, 1e6, 3e6, 11e6, 43e6,
    11e7, 26e7, 85e7, 29e8, 76e8, 25e9 };
#define NSTD ((int)(sizeof(std_b1) / sizeof(std_b1[0])))

static const char *s2_name[] = { "default", "gmp", "avx", "none", "pair" };

// ---------------------------------------------------------------------------
// parsing
// ---------------------------------------------------------------------------

static char *trim(char *s)
{
    char *e;
    while (isspace((unsigned char)*s)) s++;
    e = s + strlen(s);
    while ((e > s) && isspace((unsigned char)e[-1])) *--e = '\0';
    return s;
}

static int lower_eq(const char *a, const char *b)
{
    for (; *a && *b; a++, b++)
        if (tolower((unsigned char)*a) != tolower((unsigned char)*b))
            return 0;
    return *a == *b;
}

// a number: 11e6, 35133391030, 1.5e9
static int parse_num(const char *s, double *v)
{
    char *end;
    if (*s == '\0')
        return 0;
    *v = strtod(s, &end);
    return (*trim(end) == '\0') && (*v >= 0);
}

// a curve family: 0-4, p0-p3, param0-3, or a name; -1 if none
static int parse_family(const char *s)
{
    double v;
    if (!strncmp(s, "param", 5) && isdigit((unsigned char)s[5])) s += 5;
    else if (((s[0] == 'p') || (s[0] == 'P')) && isdigit((unsigned char)s[1])) s += 1;
    if (parse_num(s, &v) && (v == floor(v)) && (v < ECMPROB_NFAM))
        return (int)v;
    if (lower_eq(s, "suyama") || lower_eq(s, "sigma")) return ECMPROB_FAM_SUYAMA;
    if (lower_eq(s, "edw") || lower_eq(s, "edwards") || lower_eq(s, "iv") ||
        lower_eq(s, "tedw") || lower_eq(s, "edwiv"))
        return ECMPROB_FAM_EDW_IV;
    return -1;
}

// Brent-Suyama polynomial: D6 / d6 / -6 Dickson, x6 / x^6 / 6 power, 0 / 1 none
static int parse_bs(const char *s, int *S)
{
    double v;
    int sign = 1;
    if ((*s == 'D') || (*s == 'd')) { sign = -1; s++; }
    else if ((*s == 'x') || (*s == 'X')) { s++; if (*s == '^') s++; }
    else if (*s == '-') { sign = -1; s++; }
    if (!parse_num(s, &v) || (v != floor(v)) || (v > 120))
        return 0;
    *S = sign * (int)v;
    return 1;
}

// memory per thread: MB, or with a K / M / G suffix
static int parse_mem(const char *s, double *mb)
{
    char tmp[64], *e;
    double scale = 1.0;
    size_t n;
    snprintf(tmp, sizeof(tmp), "%s", s);
    n = strlen(tmp);
    if ((n > 1) && ((tmp[n - 1] == 'b') || (tmp[n - 1] == 'B')))
        tmp[--n] = '\0';
    if (n > 1)
    {
        switch (tmp[n - 1])
        {
        case 'k': case 'K': scale = 1.0 / 1024; tmp[--n] = '\0'; break;
        case 'm': case 'M': scale = 1.0; tmp[--n] = '\0'; break;
        case 'g': case 'G': scale = 1024.0; tmp[--n] = '\0'; break;
        default: break;
        }
    }
    e = trim(tmp);
    if (!parse_num(e, mb) || (*mb <= 0))
        return 0;
    *mb *= scale;
    return 1;
}

static int parse_s2(const char *s)
{
    if (lower_eq(s, "gmp") || lower_eq(s, "gmpecm")) return S2_GMP;
    if (lower_eq(s, "avx") || lower_eq(s, "avxecm") || lower_eq(s, "fft")) return S2_AVX;
    if (lower_eq(s, "none") || lower_eq(s, "0")) return S2_NONE;
    if (lower_eq(s, "pair") || lower_eq(s, "pairing")) return S2_PAIR;
    return -1;
}

// the parameters the curves run with (see usage for the defaults)
static void make_params(spec_t *sp, double B1, double B2, int fam, int s2,
    int S, double nr, const plan_env_t *env)
{
    ecm_curve_param_t *c = &sp->c;

    if (s2 == S2_DEFAULT)
        s2 = (fam == ECMPROB_FAM_EDW_IV) ? S2_AVX : S2_GMP;
    if ((B2 == 0) || ((B2 > 0) && (B2 <= B1)))
        s2 = S2_NONE;
    sp->s2 = s2;
    sp->dF = sp->k = 0;
    sp->nr_est = 0;
    sp->planned = 0;
    sp->auto_b2 = 0;
    sp->env = *env;

    memset(c, 0, sizeof(*c));
    c->B1 = B1;
    c->B2 = B1;
    c->family = fam;
    c->bsmodel = ECMPROB_BS_GMPECM;

    if (s2 == S2_GMP)
    {
        ecmprob_gmpecm_stage2(c, B1, B2 < 0 ? 0 : B2, (S == S_UNSET) ? 0 : S,
            fam, &sp->dF, &sp->k);
        if (nr >= 0)
            c->nr = nr;
    }
    else if (s2 == S2_PAIR)
    {
        // AVX-ECM's pairing stage 2: default B2 = 100 B1, no polynomial
        c->B2 = (B2 < 0) ? 100.0 * B1 : B2;
        c->bsmodel = ECMPROB_BS_AVXECM;
    }
    else if (s2 == S2_AVX)
    {
        // as yafu's ecm_avxecm_curve_params: with the default B2 the
        // planner picks B2 and the polynomial for the input size, threads
        // and memory limit; with a given B2, only the polynomial if S= isn't
        // given (none).  If the FFT stage 2 is unavailable (B1 below the
        // smallest D) the requested B2 runs without a polynomial.
        c->B2 = (B2 < 0) ? floor(pow(11.0 / 6.0 * B1, 1.43)) : B2;
        c->bsmodel = ECMPROB_BS_AVXECM;
        c->S = (S == S_UNSET) ? 0 : S;
#ifdef ECMPROB_PLANNER
        {
            uint64_t cB2;
            int bs, auto_b2 = (B2 < 0);
            double pairs;
            int req = (S == S_UNSET) ? AVXECM_BS_DEFAULT : (((S >= -1) && (S <= 1)) ? 0 : S);
            if (avxecm_s2fft_plan((uint32_t)env->bits, (uint64_t)B1, (uint64_t)c->B2,
                env->threads, env->mem_mb, auto_b2, req, fam == ECMPROB_FAM_EDW_IV,
                &cB2, &bs, &pairs))
            {
                c->B2 = (double)cB2;
                c->S = bs;
                c->nr = (nr >= 0) ? nr : pairs;
                sp->planned = 1;
                sp->auto_b2 = auto_b2;
            }
            else
                c->S = 0;
        }
        if (!sp->planned && ((c->S < -1) || (c->S > 1)))
#else
        if ((c->S < -1) || (c->S > 1))
#endif
        {
            if (nr >= 0)
                c->nr = nr;
            else
            {
                // pairs ~ (B2 - B1) phi(D) / 2D with D = 2*3*5*7*11*...
                c->nr = 0.1 * (c->B2 - B1);
                sp->nr_est = 1;
            }
        }
    }

    sp->c1 = *c;
    sp->c1.B2 = B1;
    sp->c1.S = 0;
    sp->c1.nr = 0;
}

// one curve string; returns 0 with a message in err on a bad string
static int parse_spec(const char *in, int deffam, int defs2, const plan_env_t *defenv,
    spec_t *sp, char *err, size_t errlen)
{
    char buf[256], *at, *f, *next;
    double B1 = 0, B2 = -1, nr = -1;
    int fam = -1, s2 = defs2, S = S_UNSET, npos = 0;
    plan_env_t env = *defenv;

    memset(sp, 0, sizeof(*sp));
    snprintf(sp->text, sizeof(sp->text), "%s", in);
    snprintf(buf, sizeof(buf), "%s", in);

    at = strchr(buf, '@');
    if (at == NULL)
    {
        snprintf(err, errlen, "missing '@' (curves@B1)");
        return 0;
    }
    *at = '\0';

    // curves, or done/total
    {
        char *slash = strchr(buf, '/');
        if (slash != NULL)
        {
            *slash = '\0';
            if (!parse_num(trim(buf), &sp->done) || !parse_num(trim(slash + 1), &sp->total) ||
                (sp->done > sp->total))
            {
                snprintf(err, errlen, "bad curve count (curves or done/total)");
                return 0;
            }
        }
        else
        {
            if (!parse_num(trim(buf), &sp->total))
            {
                snprintf(err, errlen, "bad curve count");
                return 0;
            }
            sp->done = sp->total;
        }
    }

    // comma-separated fields: B1 first, then positional B2 and param, or
    // key=value in any order
    for (f = at + 1; f != NULL; f = next)
    {
        char *eq, *key = NULL, *val;
        next = strchr(f, ',');
        if (next != NULL)
            *next++ = '\0';
        val = trim(f);
        eq = strchr(val, '=');
        if (eq != NULL)
        {
            *eq = '\0';
            key = trim(val);
            val = trim(eq + 1);
        }

        if (key == NULL)
        {
            // a non-numeric name is a family wherever it is
            if ((npos > 0) && !isdigit((unsigned char)val[0]) && (val[0] != '\0') &&
                (val[0] != '-') && (val[0] != '.') && (parse_family(val) >= 0))
            {
                fam = parse_family(val);
                continue;
            }
            switch (npos++)
            {
            case 0: key = "B1"; break;
            case 1: key = "B2"; break;
            case 2: key = "param"; break;
            default:
                snprintf(err, errlen, "too many fields at '%s'", val);
                return 0;
            }
        }

        if (lower_eq(key, "B1"))
        {
            if (!parse_num(val, &B1) || (B1 < 2))
            {
                snprintf(err, errlen, "bad B1 '%s'", val);
                return 0;
            }
        }
        else if (lower_eq(key, "B2"))
        {
            // empty or '-' is the default B2
            if ((val[0] == '\0') || !strcmp(val, "-") || lower_eq(val, "default"))
                B2 = -1;
            else if (!parse_num(val, &B2))
            {
                snprintf(err, errlen, "bad B2 '%s'", val);
                return 0;
            }
        }
        else if (lower_eq(key, "param") || lower_eq(key, "p") ||
            lower_eq(key, "type") || lower_eq(key, "curve"))
        {
            if ((fam = parse_family(val)) < 0)
            {
                snprintf(err, errlen, "bad curve type '%s' (0-4, p0-p3, suyama, edw)", val);
                return 0;
            }
        }
        else if (lower_eq(key, "S") || lower_eq(key, "bs"))
        {
            if (!parse_bs(val, &S))
            {
                snprintf(err, errlen, "bad Brent-Suyama polynomial '%s' (D6, x2, -6, 1)", val);
                return 0;
            }
        }
        else if (lower_eq(key, "nr") || lower_eq(key, "pairs"))
        {
            if (!parse_num(val, &nr))
            {
                snprintf(err, errlen, "bad pair count '%s'", val);
                return 0;
            }
        }
        else if (lower_eq(key, "s2"))
        {
            if ((s2 = parse_s2(val)) < 0)
            {
                snprintf(err, errlen, "bad stage 2 '%s' (gmp, avx, pair, none)", val);
                return 0;
            }
        }
        else if (lower_eq(key, "mem"))
        {
            if (!parse_mem(val, &env.mem_mb))
            {
                snprintf(err, errlen, "bad memory limit '%s' (MB per thread, or 512M, 2G)", val);
                return 0;
            }
        }
        else if (lower_eq(key, "threads") || lower_eq(key, "t"))
        {
            double v;
            if (!parse_num(val, &v) || (v < 1) || (v != floor(v)))
            {
                snprintf(err, errlen, "bad thread count '%s'", val);
                return 0;
            }
            env.threads = (int)v;
        }
        else if (lower_eq(key, "bits") || lower_eq(key, "digits") || lower_eq(key, "c"))
        {
            double v;
            if (!parse_num(val, &v) || (v < 1))
            {
                snprintf(err, errlen, "bad input size '%s'", val);
                return 0;
            }
            env.bits = lower_eq(key, "bits") ? (int)v : (int)ceil(v * 3.321928094887362);
        }
        else
        {
            snprintf(err, errlen, "unknown field '%s'", key);
            return 0;
        }
    }

    if (B1 < 2)
    {
        snprintf(err, errlen, "missing B1");
        return 0;
    }
    if (fam < 0)
        fam = deffam;
    make_params(sp, B1, B2, fam, s2, S, nr, &env);
    return 1;
}

// ---------------------------------------------------------------------------
// input: specs separated by ';' or newlines, '#' starts a comment
// ---------------------------------------------------------------------------

static spec_t *specs = NULL;
static int nspecs = 0, aspecs = 0;
static int nerrors = 0;
static int g_fam = ECMPROB_FAM_P1, g_s2 = S2_DEFAULT;
static plan_env_t g_env = { 256.0, 1, 500 };

static void add_chunk(char *s, const char *where)
{
    char *tok, *next;
    for (tok = s; tok != NULL; tok = next)
    {
        char err[160], *t;
        next = strpbrk(tok, ";\n");
        if (next != NULL)
            *next++ = '\0';
        if ((t = strchr(tok, '#')) != NULL)
            *t = '\0';
        t = trim(tok);
        if (*t == '\0')
            continue;
        if (nspecs == aspecs)
        {
            aspecs = aspecs ? 2 * aspecs : 16;
            specs = (spec_t *)realloc(specs, (size_t)aspecs * sizeof(spec_t));
        }
        if (parse_spec(t, g_fam, g_s2, &g_env, &specs[nspecs], err, sizeof(err)))
            nspecs++;
        else
        {
            fprintf(stderr, "ecmprob: %s: '%s': %s\n", where, t, err);
            nerrors++;
        }
    }
}

static void add_file(FILE *fp, const char *where)
{
    char line[1024];
    while (fgets(line, sizeof(line), fp) != NULL)
        add_chunk(line, where);
}

// ---------------------------------------------------------------------------
// output
// ---------------------------------------------------------------------------

static void fmt_bound(char *out, size_t len, double x)
{
    // compact exact form when there is one (11e6, 25e4, 35e9), else digits
    double m = floor(x + 0.5);
    int e = 0;
    if ((m != x) || (m < 1e3) || (m >= 1e18))
    {
        snprintf(out, len, (x < 1e18) ? "%.0f" : "%.4g", x);
        return;
    }
    while ((m >= 10) && (fmod(m, 10) == 0)) { m /= 10; e++; }
    if ((e >= 3) && (m < 1000))
        snprintf(out, len, "%.0fe%d", m, e);
    else
        snprintf(out, len, "%.0f", x);
}

static void fmt_poly(char *out, size_t len, int S)
{
    if (S > 1) snprintf(out, len, "x^%d", S);
    else if (S < -1) snprintf(out, len, "D%d", -S);
    else snprintf(out, len, "-");
}

static void fmt_curves(char *out, size_t len, double n)
{
    if (n >= 1e299) snprintf(out, len, "inf");
    else if (n < 1e7) snprintf(out, len, "%.0f", ceil(n));
    else snprintf(out, len, "%.3g", n);
}

static void print_spec(const spec_t *sp)
{
    char b1[32], b2[32], poly[16], t[32];
    double w[ECMPROB_NDIG] = { 0 }, tl;

    fmt_bound(b1, sizeof(b1), sp->c.B1);
    fmt_bound(b2, sizeof(b2), sp->c.B2);
    fmt_poly(poly, sizeof(poly), sp->c.S);
    ecmprob_work_add(w, &sp->c, sp->done);
    ecmprob_work_add(w, &sp->c1, sp->total - sp->done);
    tl = ecmprob_work_tlevel(w);

    printf("  %s\n", sp->text);
    printf("      %s, delta %.3f; ", ecmprob_family_name(sp->c.family),
        ecmprob_delta(sp->c.family));
    if (sp->total > sp->done)
        printf("%.0f of %.0f curves with stage 2; ", sp->done, sp->total);
    else
        printf("%.0f curves; ", sp->total);
    printf("t%.2f alone\n", tl);
    if (sp->s2 == S2_NONE)
        printf("      B1=%s, no stage 2\n", b1);
    else if (sp->s2 == S2_GMP)
        printf("      B1=%s, B2=%s (GMP-ECM stage 2, dF=%u, k=%u), poly %s, %.4g pairs\n",
            b1, b2, sp->dF, sp->k, poly, sp->c.nr);
    else if (sp->s2 == S2_PAIR)
        printf("      B1=%s, B2=%s (AVX-ECM pairing stage 2)\n", b1, b2);
    else if (sp->planned)
        printf("      B1=%s, B2=%s (AVX-ECM FFT stage 2, %s B2; %d bits, %d thread%s, %.0f MB/thread),"
            " poly %s, %.4g pairs\n", b1, b2, sp->auto_b2 ? "automatic" : "requested",
            sp->env.bits, sp->env.threads, sp->env.threads > 1 ? "s" : "", sp->env.mem_mb,
            poly, sp->c.nr);
    else
        printf("      B1=%s, B2=%s (AVX-ECM stage 2, not planned%s), poly %s%s%s\n",
#ifdef ECMPROB_PLANNER
            b1, b2, ": FFT unavailable at this B1",
#else
            b1, b2, ": built without the planner, fixed default B2",
#endif
            poly, ((sp->c.S < -1) || (sp->c.S > 1)) ? ", pairs " : "",
            ((sp->c.S < -1) || (sp->c.S > 1)) ? (sp->nr_est ? "~0.1(B2-B1) (estimated)" : "given") : "");

    // expected curves per t-level over this B1's useful range: from the
    // first size that takes more than a few curves
    {
        int d, d0 = 15;
        while ((d0 < ECMPROB_DMAX - 20) && (ecmprob_curves(&sp->c, (double)d0) < 10))
            d0 += 5;
        printf("      expected curves:");
        for (d = d0; d <= d0 + 20; d += 5)
        {
            fmt_curves(t, sizeof(t), ecmprob_curves(&sp->c, (double)d));
            printf(" t%d %s%s", d, t, d < d0 + 20 ? "," : "\n");
        }
    }
}

static void usage(void)
{
    printf(
"ecmprob: ECM t-level and missed-factor probabilities (yafu's model)\n"
"\n"
"usage: ecmprob [options] curve-string ...\n"
"       ecmprob [options] -i file          (or curve strings on stdin)\n"
"\n"
"curve string:  [done/]curves@B1[,B2][,param][,key=value ...]\n"
"               keys: B1 B2 p S nr s2 mem threads bits digits\n"
"  curves      curves run; done/total: all 'total' ran stage 1, 'done' of them stage 2\n"
"  B1, B2      plain or e-notation (11e6, 35e9); B1=, B2= prefixes optional.\n"
"              B2 omitted, empty or '-': the stage 2's default; B2=0: stage 1 only\n"
"  param       curve type (p=, param= prefix optional), number or name:\n"
"                0, p0, suyama  Suyama sigma (GMP-ECM param 0, AVX-ECM -avxecm_curve 0)\n"
"                1, p1          GMP-ECM param 1 (GMP-ECM's default on the CPU)\n"
"                2, p2          GMP-ECM param 2\n"
"                3, p3          GMP-ECM param 3 (GPU)\n"
"                4, edw, iv     twisted Edwards family (iv), AVX-ECM's default\n"
"  s2=         stage 2 the curves ran (default: avx for edw, gmp otherwise)\n"
"                gmp   GMP-ECM's NTT stage 2: B2 is rounded up to what it\n"
"                      covers, default B2 and Dickson polynomial as GMP-ECM\n"
"                avx   AVX-ECM FFT stage 2, planned as yafu plans it%s:\n"
"                      default B2 and polynomial chosen for the input size,\n"
"                      threads and memory limit (mem=, threads=, bits=)\n"
"                pair  AVX-ECM pairing stage 2 (-avxecm_stg2 0), default B2 100 B1\n"
"                none  stage 1 only (same as B2=0)\n"
"  S=          Brent-Suyama polynomial: D6 (or -6) Dickson, x2 (or 2) power, 1 none\n"
"  nr=         stage 2 pairs for Brent-Suyama (default: from the plan)\n"
"  mem=        s2=avx memory limit per thread: MB, or 512M, 2G (default 256,\n"
"              yafu's -avxecm_stg2mem)\n"
"  threads=    s2=avx threads (default 1, yafu's -threads)\n"
"  bits=, digits=  s2=avx input size (default 500 bits)\n"
"  Strings are separated by ';', newlines or arguments; '#' starts a comment.\n"
"\n"
"options:\n"
"  -i file     read curve strings from file ('-' for stdin)\n"
"  -p param    default curve type for strings without one (default 1)\n"
"  -s2 kind    default stage 2 (gmp, avx, pair, none)\n"
"  -mem MB     default s2=avx memory limit per thread (512M, 2G also work)\n"
"  -threads n  default s2=avx threads\n"
"  -bits n, -digits n   default s2=avx input size\n"
"  -r lo,hi    table rows from lo to hi digits (default: around the t-level)\n"
"  -d step     table row step in digits (default 5)\n"
"  -t level    suggest curves to reach this t-level (at -b B1, or each standard B1)\n"
"  -b B1       B1 for -t\n"
"  -q          print only the t-level (tNN.NNN)\n"
"  -v          show each string's parameters and expected curves\n"
"  -h          this help\n"
"\n"
"The t-level is the factor size at which the expected number of finds is 1\n"
"(miss chance exp(-1) = 36.8%%); 't-mult' is the expected finds at each size,\n"
"the multiple of that t-level covered, and P(miss) = exp(-t-mult).\n"
"\n"
"examples:\n"
"  ecmprob 5208@11e6                     GMP-ECM param 1, default B2\n"
"  ecmprob \"7550@43e6;5000@11e7,p=3\"\n"
"  ecmprob -v 2000@3e6,edw,mem=1G,digits=150   yafu's default stage 2 plan\n"
"  ecmprob 2000@3e6,5.2e9,edw,S=D6,nr=5.3e8    a stage 2 given exactly\n"
"  ecmprob 500/1000@11e7,p=3            GPU stage 1, half had stage 2\n"
"  ecmprob -t 55 -p edw -mem 1024 7550@43e6\n",
#ifdef ECMPROB_PLANNER
        ""
#else
        " (this build has no planner:\n"
        "                      fixed B2 (11/6 B1)^1.43; build with make ecmprob)"
#endif
        );
}

int main(int argc, char **argv)
{
    int i, quiet = 0, verbose = 0, any_input = 0, step = 5, lo = -1, hi = -1;
    double target = 0, sugB1 = 0;
    double w[ECMPROB_NDIG] = { 0 }, tl;

    ecmprob_init();

    // options first (they set the defaults the strings are parsed with)
    for (i = 1; i < argc; i++)
    {
        const char *a = argv[i];
        if (a[0] != '-' || a[1] == '\0')
            continue;
        if (!strcmp(a, "-h") || !strcmp(a, "--help")) { usage(); return 0; }
        else if (!strcmp(a, "-q")) quiet = 1;
        else if (!strcmp(a, "-v")) verbose = 1;
        else if ((i + 1 < argc) && (!strcmp(a, "-p") || !strcmp(a, "-s2") ||
            !strcmp(a, "-r") || !strcmp(a, "-d") || !strcmp(a, "-t") ||
            !strcmp(a, "-b") || !strcmp(a, "-i") || !strcmp(a, "-mem") ||
            !strcmp(a, "-threads") || !strcmp(a, "-bits") || !strcmp(a, "-digits")))
        {
            const char *v = argv[++i];
            if (!strcmp(a, "-mem"))
            {
                if (!parse_mem(v, &g_env.mem_mb))
                {
                    fprintf(stderr, "ecmprob: bad memory limit '%s'\n", v);
                    return 1;
                }
            }
            else if (!strcmp(a, "-threads") || !strcmp(a, "-bits") || !strcmp(a, "-digits"))
            {
                int n = atoi(v);
                if (n < 1)
                {
                    fprintf(stderr, "ecmprob: bad value for %s: '%s'\n", a, v);
                    return 1;
                }
                if (!strcmp(a, "-threads")) g_env.threads = n;
                else if (!strcmp(a, "-bits")) g_env.bits = n;
                else g_env.bits = (int)ceil(n * 3.321928094887362);
            }
            else if (!strcmp(a, "-p"))
            {
                if ((g_fam = parse_family(v)) < 0)
                {
                    fprintf(stderr, "ecmprob: bad curve type '%s'\n", v);
                    return 1;
                }
            }
            else if (!strcmp(a, "-s2"))
            {
                if ((g_s2 = parse_s2(v)) < 0)
                {
                    fprintf(stderr, "ecmprob: bad stage 2 '%s'\n", v);
                    return 1;
                }
            }
            else if (!strcmp(a, "-r"))
            {
                if ((sscanf(v, "%d,%d", &lo, &hi) != 2) || (lo > hi))
                {
                    fprintf(stderr, "ecmprob: -r lo,hi\n");
                    return 1;
                }
            }
            else if (!strcmp(a, "-d"))
            {
                step = atoi(v);
                if (step < 1) step = 1;
            }
            else if (!strcmp(a, "-t")) target = atof(v);
            else if (!strcmp(a, "-b"))
            {
                if (!parse_num(v, &sugB1) || (sugB1 < 2))
                {
                    fprintf(stderr, "ecmprob: bad B1 '%s'\n", v);
                    return 1;
                }
            }
        }
        else
        {
            fprintf(stderr, "ecmprob: unknown option or missing value: %s (-h for help)\n", a);
            return 1;
        }
    }

    // then the curve strings and files, in order
    for (i = 1; i < argc; i++)
    {
        const char *a = argv[i];
        if ((a[0] == '-') && (a[1] != '\0'))
        {
            if (!strcmp(a, "-i"))
            {
                const char *fn = argv[++i];
                any_input = 1;
                if (!strcmp(fn, "-"))
                    add_file(stdin, "stdin");
                else
                {
                    FILE *fp = fopen(fn, "r");
                    if (fp == NULL)
                    {
                        fprintf(stderr, "ecmprob: can't open %s\n", fn);
                        return 1;
                    }
                    add_file(fp, fn);
                    fclose(fp);
                }
            }
            else if (strcmp(a, "-q") && strcmp(a, "-v") && strcmp(a, "-h"))
                i++;        // an option with a value, handled above
            continue;
        }
        {
            char *s = (char *)malloc(strlen(a) + 1);
            strcpy(s, a);
            any_input = 1;
            add_chunk(s, "argument");
            free(s);
        }
    }
    if (!any_input)
    {
        if (ISATTY_STDIN() && (target <= 0))
        {
            usage();
            return 0;
        }
        if (!ISATTY_STDIN())
            add_file(stdin, "stdin");
    }
    if (nerrors > 0)
        return 1;

    for (i = 0; i < nspecs; i++)
    {
        ecmprob_work_add(w, &specs[i].c, specs[i].done);
        ecmprob_work_add(w, &specs[i].c1, specs[i].total - specs[i].done);
    }
    tl = ecmprob_work_tlevel(w);

    if (quiet)
    {
        printf("t%.3f\n", tl);
    }
    else
    {
        int d;
        if (verbose && (nspecs > 0))
        {
            printf("curves:\n");
            for (i = 0; i < nspecs; i++)
                print_spec(&specs[i]);
            printf("\n");
        }

        if (nspecs > 0)
        {
            // rows: from two steps below the t-level row until a factor would
            // almost surely have been missed
            int fixed = (lo >= 0);
            if (!fixed)
            {
                lo = step * (int)(tl / step) - 2 * step;
                if (lo < 15) lo = 15;
                hi = ECMPROB_DMAX;
            }
            if (lo < ECMPROB_DMIN) lo = ECMPROB_DMIN;
            if (hi > ECMPROB_DMAX) hi = ECMPROB_DMAX;

            printf("  digits    t-mult     P(miss)\n");
            for (d = lo; d <= hi; d += step)
            {
                double x = ecmprob_work_at(w, (double)d), m = exp(-x);
                if (m >= 0.0001)
                    printf("  %5d  %9.3f  %9.2f%%\n", d, x, 100.0 * m);
                else
                    printf("  %5d  %9.3f  %10s\n", d, x, "<0.01%");
                if (!fixed && (m > 0.99) && (d >= tl))
                    break;
            }
            printf("t-level: t%.2f  (expected finds 1, P(miss) 36.8%%)\n", tl);
        }
    }

    // suggestions to reach a target t-level
    if (target > 0)
    {
        spec_t sp;
        int n1 = (sugB1 > 0) ? 1 : NSTD;
        if (tl >= target)
        {
            printf("already at t%.2f >= t%.2f\n", tl, target);
            return 0;
        }
        printf("to reach t%.2f from t%.2f, any one of:\n", target, tl);
        for (i = 0; i < n1; i++)
        {
            double B1 = (sugB1 > 0) ? sugB1 : std_b1[i];
            char b1[32], fam[8];
            uint32_t n;
            make_params(&sp, B1, -1, g_fam, g_s2, S_UNSET, -1, &g_env);
            n = ecmprob_work_curves_to(w, &sp.c, target);
            if (n == 0xffffffffu)
            {
                if (sugB1 > 0)
                    printf("  t%.2f is out of reach at B1 = %.4g\n", target, B1);
                continue;
            }
            // without -b, skip B1 far below (millions of curves) or far
            // above (a handful of curves) the target's useful range
            if ((sugB1 <= 0) && ((n > 2000000) || (n < 4)))
                continue;
            fmt_bound(b1, sizeof(b1), B1);
            if (g_fam == ECMPROB_FAM_EDW_IV)
                snprintf(fam, sizeof(fam), "edw");
            else
                snprintf(fam, sizeof(fam), "%d", g_fam);
            printf("  %u@%s,p=%s%s%s", n, b1, fam,
                (g_s2 != S2_DEFAULT) ? ",s2=" : "",
                (g_s2 != S2_DEFAULT) ? s2_name[g_s2] : "");
            // the plan's inputs, when not the defaults, so the string
            // reproduces this plan
            if (sp.planned)
            {
                if (g_env.mem_mb != 256.0) printf(",mem=%.0f", g_env.mem_mb);
                if (g_env.threads != 1) printf(",threads=%d", g_env.threads);
                if (g_env.bits != 500) printf(",bits=%d", g_env.bits);
            }
            printf("\n");
        }
    }
    return 0;
}
