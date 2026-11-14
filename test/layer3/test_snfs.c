/*----------------------------------------------------------------------
 Layer 3 -- SNFS special-form detection, algebraic factor removal and
 polynomial construction, plus the factor() decisions that depend on them.

 Two kinds of tests:

 1) In-process (always run).  Mirrors what factor()'s determine_job_type()
    does: the full form is the input (fobj->input_N), the number being
    worked on (fobj->nfs_obj.gmp_n) is that input with primes < 10^5 removed
    (or an explicit cofactor), then

        poly   = snfs_find_form(fobj, input_N);
        status = snfs_choose_poly(fobj, &job, poly, 0);

    and checks:
      - the expected special form is detected (brent, homogeneous
        cunningham, xyyxf, direct);
      - nothing is lost: the product of the factor list is still the input,
        and the reduced input is one of its entries (algebraic, Aurifeuillian
        and difference-of-powers splits only refine the list);
      - the expected splits happened (pinned: size of the reduced input,
        number of new factor-list entries);
      - the polynomial is for the reduced input (n == gmp_n, n > 1) and
        f, g share a root mod n: sum c_i (-Y0)^i Y1^(d-i) == 0 (mod n);
      - degree and difficulty (pinned).
    Constructions covered: brent with prime / composite exponent (one to four
    distinct odd primes, even exponents with mixed algebraic halves), with
    coefficients, composite and large bases; homogeneous cunninghams;
    xyyxf; direct forms; Aurifeuillian forms (brent and homogeneous); sums /
    differences of k-th powers with non-unit coefficients; prime exponents
    with no algebraic factor at degrees 4, 5 and 6.

 2) Black-box (run when the yafu binary is found).  Runs the built yafu
    in a scratch directory (no yafu.ini there, so defaults apply) and checks
    the command line options: snfs / gnfs / siqs / ecm method choice with
    -gnfs, -xover, -snfs_xover, -pretest; the SNFS ECM depth; -np job files
    from factor() and snfs(); and an SNFS input where ECM with a fixed
    -sigma finds a factor and the job degenerates to SIQS (-max_siqs stops it
    before sieving).  The binary is $YAFU_BIN, else ./yafu.exe / ./yafu in
    the current directory (run the test from the yafu source root).  Without
    a binary these tests print a note and pass.  Scratch files go to
    ./yafu_test_snfs.tmp/<test>/.

 The in-process tests leave YAFU_get_poly_score.out (msieve's poly scoring
 scratch file) in the working directory.
 Public domain.
----------------------------------------------------------------------*/
#include "testkit.h"
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "factor.h"
#include "nfs_impl.h"

#if defined(_WIN32)
#include <direct.h>
#define TS_MKDIR(p) _mkdir(p)
#define TS_GETCWD(b, n) _getcwd(b, n)
#define TS_EXE ".exe"
#else
#include <sys/stat.h>
#include <unistd.h>
#define TS_MKDIR(p) mkdir(p, 0755)
#define TS_GETCWD(b, n) getcwd(b, n)
#define TS_EXE ""
#endif

/* ------------------------------------------------------------------ *
 * a tiny integer expression evaluator: + - * / ^ ( ) and decimals.
 * '/' must be exact.  Returns 0 on success.
 * ------------------------------------------------------------------ */
typedef struct { const char *p; int err; } ts_lex;

static void ts_expr(ts_lex *L, mpz_t r);

static void ts_skip(ts_lex *L) { while (isspace((unsigned char)*L->p)) L->p++; }

static void ts_atom(ts_lex *L, mpz_t r)
{
    ts_skip(L);
    if (*L->p == '(') {
        L->p++;
        ts_expr(L, r);
        ts_skip(L);
        if (*L->p == ')') L->p++; else L->err = 1;
    } else if (*L->p == '-') {
        L->p++;
        ts_atom(L, r);
        mpz_neg(r, r);
    } else if (isdigit((unsigned char)*L->p)) {
        const char *s = L->p;
        char *buf;
        while (isdigit((unsigned char)*L->p)) L->p++;
        buf = (char *)malloc((size_t)(L->p - s) + 1);
        memcpy(buf, s, (size_t)(L->p - s));
        buf[L->p - s] = 0;
        mpz_set_str(r, buf, 10);
        free(buf);
    } else {
        L->err = 1;
        mpz_set_ui(r, 0);
    }
}

static void ts_power(ts_lex *L, mpz_t r)
{
    ts_atom(L, r);
    ts_skip(L);
    if (*L->p == '^') {
        mpz_t e;
        L->p++;
        mpz_init(e);
        ts_power(L, e);                 /* right associative */
        if (!mpz_fits_ulong_p(e)) L->err = 1;
        else mpz_pow_ui(r, r, mpz_get_ui(e));
        mpz_clear(e);
    }
}

static void ts_term(ts_lex *L, mpz_t r)
{
    mpz_t t;
    mpz_init(t);
    ts_power(L, r);
    for (;;) {
        char op;
        ts_skip(L);
        op = *L->p;
        if (op != '*' && op != '/') break;
        L->p++;
        ts_power(L, t);
        if (op == '*') mpz_mul(r, r, t);
        else if (mpz_sgn(t) == 0 || !mpz_divisible_p(r, t)) L->err = 1;
        else mpz_divexact(r, r, t);
    }
    mpz_clear(t);
}

static void ts_expr(ts_lex *L, mpz_t r)
{
    mpz_t t;
    mpz_init(t);
    ts_term(L, r);
    for (;;) {
        char op;
        ts_skip(L);
        op = *L->p;
        if (op != '+' && op != '-') break;
        L->p++;
        ts_term(L, t);
        if (op == '+') mpz_add(r, r, t); else mpz_sub(r, r, t);
    }
    mpz_clear(t);
}

static int ts_eval(mpz_t r, const char *s)
{
    ts_lex L;
    L.p = s;
    L.err = 0;
    ts_expr(&L, r);
    ts_skip(&L);
    return (L.err || *L.p != 0) ? -1 : 0;
}

/* f(m) == 0 mod n for the homogeneous pair: sum c_i (-Y0)^i Y1^(d-i) */
static int ts_common_root(mpz_t n, mpz_t *c, int deg, mpz_t Y0, mpz_t Y1)
{
    mpz_t acc, t, a, b;
    int i, ok;
    mpz_inits(acc, t, a, b, NULL);
    mpz_set_ui(acc, 0);
    for (i = 0; i <= deg; i++) {
        mpz_neg(a, Y0);
        mpz_powm_ui(a, a, (unsigned long)i, n);
        mpz_powm_ui(b, Y1, (unsigned long)(deg - i), n);
        mpz_mul(t, a, b);
        mpz_mul(t, t, c[i]);
        mpz_add(acc, acc, t);
        mpz_mod(acc, acc, n);
    }
    ok = (mpz_sgn(acc) == 0);
    mpz_clears(acc, t, a, b, NULL);
    return ok;
}

/* ------------------------------------------------------------------ *
 * In-process cases
 * ------------------------------------------------------------------ */
#define ST_POLY  1      /* snfs_choose_poly returned a polynomial        */
#define ST_DONE -2      /* splits finished the input (rest is prp)        */
#define ST_NONE  0      /* no usable polynomial                           */
#define ANY     -1      /* don't check                                    */

typedef struct {
    const char *name;
    const char *expr;     /* full form: fobj->input_N                       */
    const char *cof;      /* working cofactor, or NULL: expr less p < 1e5   */
    int form;             /* expected snfs_form_e                           */
    int status;           /* ST_*                                           */
    int ndigits;          /* digits of the reduced input / poly n, or ANY   */
    int npieces;          /* new factor-list entries from splits, or ANY    */
    int degree;           /* algebraic degree, or ANY                       */
    double difficulty;    /* SNFS difficulty (+-0.01), or < 0 for ANY        */
} snfs_case;

/* the Phi_63(10) * Phi_126(10) product: divides 10^126 - 1, draws one factor
   from each of 10^63 - 1 and 10^63 + 1 (the even-exponent halving bug).    */
#define MIXED_63_126 "999999000000000000999999000000000000999999999999000000999999999999000001"

/* cofactor of 30^74 - 50625 from mersenneforum (P52 * P40) */
#define C92_DIFF_SQ "11122701950483278721937968112660782095113141752223238370030932602445600600874494559062373921"

static const snfs_case cases[] = {
    /* --- brent, prime exponent, no algebraic factor: degrees 4, 5, 7 --- */
    { "brent_2_401m",   "2^401-1",   NULL, SNFS_BRENT, ST_POLY, 121, 0, 4, 120.71 },
    { "brent_2_613m",   "2^613-1",   NULL, SNFS_BRENT, ST_POLY, 185, 0, 5, 184.53 },
    { "brent_2_1061m",  "2^1061-1",  NULL, SNFS_BRENT, ST_POLY, 320, 0, 7, 319.39 },
    /* --- brent with leading coefficient, composite base --- */
    { "brent_coeff",    "275*2^335-1",   NULL, SNFS_BRENT, ST_POLY, 102, 0, 4, 103.28 },
    { "brent_coeff_b5", "374*5^374+1",   NULL, SNFS_BRENT, ST_POLY, 261, 0, 6, 263.99 },
    { "brent_cbase_c2", "35^120+7",      NULL, SNFS_BRENT, ST_POLY, 180, 0, 5, 185.29 },
    /* --- large bases (x - 1 removed; quartic via Phi_5, sextic) --- */
    { "brent_big_b5",   "5444131015388735415835538025644022987013813^5-1", NULL, SNFS_BRENT, ST_POLY, 168, 1, 4, 170.94 },
    { "brent_big_b7",   "70332287049806043638180632371100261^7-1",         NULL, SNFS_BRENT, ST_POLY, 210, 1, 6, 209.08 },
    /* --- composite exponents: cyclotomic splits --- */
    { "brent_e1155p",   "2^1155+1",  NULL, SNFS_BRENT, ST_POLY, 137, 6, 4, 185.43 },
    { "brent_e1296m",   "2^1296-1",  NULL, SNFS_BRENT, ST_POLY, 127, 9, 4, 130.04 },
    { "brent_e1208m",   "2^1208-1",  NULL, SNFS_BRENT, ST_POLY, 178, 3, 5, 181.82 },
    { "brent_10_222m",  "10^222-1",  NULL, SNFS_BRENT, ST_POLY, 71, 3, 4, 74.00 },
    { "brent_10_126m",  "10^126-1",  NULL, SNFS_BRENT, ST_POLY, 37, 5, 6, 36.00 },
    { "brent_mixed",    MIXED_63_126, NULL, SNFS_BRENT, ST_POLY, 37, 1, 6, 36.00 },
    /* --- homogeneous cunningham --- */
    { "hcunn_13_10",    "13^199-10^199",  NULL, SNFS_H_CUNNINGHAM, ST_POLY, 219, 0, 6, 221.67 },
    { "hcunn_50_49e",   "50^124+49^124",  NULL, SNFS_H_CUNNINGHAM, ST_POLY, 204, 1, 5, 213.06 },
    /* --- Aurifeuillian (brent and homogeneous) --- */
    { "aurif_2_398p",   "2^398+1",        NULL, SNFS_BRENT, ST_DONE, 0, 1, ANY, -1 },
    { "aurif_2_990p",   "2^990+1",        NULL, SNFS_BRENT, ST_POLY, 73, 6, 8, 79.47 },
    { "aurif_2_1086p",  "2^1086+1",       NULL, SNFS_BRENT, ST_DONE, 0, 2, ANY, -1 },
    { "aurif_h_50_94",  "50^94+49^94",    NULL, SNFS_H_CUNNINGHAM, ST_POLY, 79, 1, 6, 83.24 },
    { "aurif_h_50_142", "50^142+49^142",  NULL, SNFS_H_CUNNINGHAM, ST_POLY, 119, 1, 6, 124.02 },
    { "aurif_h_28_9",   "28^161+9^161",   NULL, SNFS_H_CUNNINGHAM, ST_POLY, 96, 2, 6, 99.85 },
    /* --- xyyxf and direct --- */
    { "xyyxf_84_83",    "84^83+83^84",    NULL, SNFS_XYYXF, ST_POLY, 159, 0, 5, 163.18 },
    { "direct_16_19",   "(2^152+9)*2^152+1", NULL, SNFS_DIRECT, ST_POLY, 92, 0, 4, 91.51 },
    /* --- sums / differences of k-th powers with non-unit coefficients --- */
    { "pow_diff_sq",    "30^74-50625", C92_DIFF_SQ, SNFS_BRENT, ST_DONE, 0, 1, ANY, -1 },
    { "pow_diff_cube",  "8*7^99-27",   NULL,        SNFS_BRENT, ST_DONE, 0, 1, ANY, -1 },
};
#define NCASES ((int)(sizeof cases / sizeof cases[0]))

static int ts_report = 0;    /* set YAFU_TEST_SNFS_REPORT=1 to print actuals */

static void snfs_case_run(tk_ctx *tk, const snfs_case *c)
{
    mpz_t N, w, p, prod, Y0, Y1;
    fact_obj_t *fobj;
    snfs_t *poly;
    nfs_job_t job;
    int status, list0, deg = 0, nd = 0, form;
    double diff = 0.0;

    {
        const char *rep = getenv("YAFU_TEST_SNFS_REPORT");
        ts_report = (rep != NULL && rep[0] != 0 && rep[0] != '0');
    }

    mpz_inits(N, w, p, prod, Y0, Y1, NULL);
    TK_REQUIRE(tk, ts_eval(N, c->expr) == 0, "bad expression %s", c->expr);

    fobj = new_default_factorization(N);
    TK_REQUIRE(tk, fobj != NULL, "new_default_factorization returned NULL");
    fobj->VFLAG = -1;
    fobj->LOGFLAG = 0;
    fobj->autofact_obj.autofact_active = 1;

    /* the working cofactor: like factor() after trial division */
    if (c->cof != NULL) {
        TK_REQUIRE(tk, ts_eval(w, c->cof) == 0, "bad cofactor %s", c->cof);
        TK_REQUIRE(tk, mpz_divisible_p(N, w), "cofactor does not divide the input");
        mpz_divexact(p, N, w);
        add_to_factor_list(fobj->factors, p, fobj->VFLAG, fobj->NUM_WITNESSES, 0);
    } else {
        mpz_set(w, N);
        mpz_set_ui(p, 2);
        while (mpz_cmp_ui(p, 100000) < 0) {
            mpz_t pe;
            unsigned long e;
            mpz_init(pe);
            e = mpz_remove(w, w, p);
            if (e > 0) {
                mpz_pow_ui(pe, p, e);
                add_to_factor_list(fobj->factors, pe, fobj->VFLAG, fobj->NUM_WITNESSES, 0);
            }
            mpz_clear(pe);
            mpz_nextprime(p, p);
        }
    }
    mpz_set(fobj->N, w);
    mpz_set(fobj->nfs_obj.gmp_n, w);
    list0 = fobj->factors->num_factors;

    poly = snfs_find_form(fobj, fobj->input_N);
    TK_REQUIRE(tk, poly != NULL, "no special form detected for %s", c->expr);
    form = (int)poly->form_type;
    TK_CHECKF(tk, form == c->form, "%s: form type %d, expected %d", c->expr, form, c->form);

    memset(&job, 0, sizeof(job));
    status = snfs_choose_poly(fobj, &job, poly, 0);
    snfs_clear(poly);
    free(poly);

    /* nothing lost: the factor list still multiplies to the input */
    get_prod_of_factors(fobj->factors, prod);
    TK_CHECKF(tk, mpz_cmp(prod, N) == 0, "%s: factor list product != input", c->expr);

    /* the reduced input is what's left to factor: an entry of the list */
    if (mpz_cmp_ui(fobj->nfs_obj.gmp_n, 1) > 0)
        TK_CHECKF(tk, find_in_factor_list(fobj->factors, fobj->nfs_obj.gmp_n) >= 0,
            "%s: reduced input is not a factor-list entry", c->expr);

    if (status > 0) {
        snfs_t *s = job.snfs;
        TK_REQUIRE(tk, s != NULL, "status %d but no polynomial", status);
        deg = (int)s->poly->alg.degree;
        diff = s->difficulty;
        nd = (int)mpz_sizeinbase(s->n, 10);

        TK_CHECKF(tk, mpz_cmp_ui(s->n, 1) > 0, "%s: polynomial n <= 1", c->expr);
        TK_CHECKF(tk, mpz_cmp(s->n, fobj->nfs_obj.gmp_n) == 0,
            "%s: polynomial n differs from the reduced input", c->expr);
        TK_CHECKF(tk, mpz_divisible_p(N, s->n), "%s: polynomial n does not divide the input", c->expr);
        TK_CHECKF(tk, s->poly->rat.degree == 1, "%s: rational poly degree %u", c->expr,
            s->poly->rat.degree);
        TK_CHECKF(tk, deg >= 2 && deg <= 8, "%s: algebraic degree %d", c->expr, deg);

        mpz_set(Y0, s->poly->rat.coeff[0]);
        mpz_set(Y1, s->poly->rat.coeff[1]);
        TK_CHECKF(tk, mpz_sgn(Y1) != 0, "%s: Y1 == 0", c->expr);
        if (mpz_cmp_ui(s->n, 1) > 0)
            TK_CHECKF(tk, ts_common_root(s->n, s->poly->alg.coeff, deg, Y0, Y1),
                "%s: f and g have no common root mod n", c->expr);

        snfs_clear(job.snfs);
        free(job.snfs);
    } else {
        nd = (mpz_cmp_ui(fobj->nfs_obj.gmp_n, 1) > 0) ?
            (int)mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10) : 0;
    }

    if (c->status != ANY)
        TK_CHECKF(tk, (status > 0 ? ST_POLY : status) == c->status,
            "%s: snfs_choose_poly status %d, expected %d", c->expr, status, c->status);
    if (status == ST_DONE)
        TK_CHECKF(tk, mpz_cmp_ui(fobj->nfs_obj.gmp_n, 1) == 0,
            "%s: finished (-2) but input not reduced to 1", c->expr);
    if (c->ndigits != ANY)
        TK_CHECKF(tk, nd == c->ndigits, "%s: reduced input has %d digits, expected %d",
            c->expr, nd, c->ndigits);
    if (c->npieces != ANY)
        TK_CHECKF(tk, fobj->factors->num_factors - list0 == c->npieces,
            "%s: %d new factor-list entries, expected %d", c->expr,
            fobj->factors->num_factors - list0, c->npieces);
    if (c->degree != ANY && status > 0)
        TK_CHECKF(tk, deg == c->degree, "%s: degree %d, expected %d", c->expr, deg, c->degree);
    if (c->difficulty >= 0 && status > 0)
        TK_CHECKF(tk, diff > c->difficulty - 0.01 && diff < c->difficulty + 0.01,
            "%s: difficulty %1.2f, expected %1.2f", c->expr, diff, c->difficulty);

    if (ts_report || tk_verbose(tk)) {
        printf("        %-15s form %d status %2d  n c%-3d  +%d pieces  deg %d  diff %1.2f\n",
            c->name, form, status > 0 ? 1 : status, nd,
            fobj->factors->num_factors - list0, deg, diff);
        fflush(stdout);
    }

    free_factobj(fobj);
    mpz_clears(N, w, p, prod, Y0, Y1, NULL);
}

#define CASE_FN(k) static void t_case_##k(tk_ctx *tk) { \
    if (k < NCASES) snfs_case_run(tk, &cases[k]); }
CASE_FN(0)  CASE_FN(1)  CASE_FN(2)  CASE_FN(3)  CASE_FN(4)  CASE_FN(5)
CASE_FN(6)  CASE_FN(7)  CASE_FN(8)  CASE_FN(9)  CASE_FN(10) CASE_FN(11)
CASE_FN(12) CASE_FN(13) CASE_FN(14) CASE_FN(15) CASE_FN(16) CASE_FN(17)
CASE_FN(18) CASE_FN(19) CASE_FN(20) CASE_FN(21) CASE_FN(22) CASE_FN(23)
CASE_FN(24) CASE_FN(25)

/* ------------------------------------------------------------------ *
 * Black-box: run the yafu binary
 * ------------------------------------------------------------------ */
static char ts_bin[4096];
static int  ts_bin_state = 0;   /* 0 unknown, 1 found, -1 missing */

static int ts_find_bin(void)
{
    const char *env;
    char cwd[3072];
    FILE *f;

    if (ts_bin_state != 0)
        return ts_bin_state > 0;

    env = getenv("YAFU_BIN");
    if (env != NULL && env[0] != 0) {
        snprintf(ts_bin, sizeof ts_bin, "%s", env);
    } else if (TS_GETCWD(cwd, sizeof cwd) != NULL) {
        snprintf(ts_bin, sizeof ts_bin, "%s/yafu%s", cwd, TS_EXE);
    } else {
        ts_bin_state = -1;
        return 0;
    }

    f = fopen(ts_bin, "rb");
    if (f == NULL) {
        printf("        note: yafu binary not found at %s (set YAFU_BIN); "
            "skipping command-line tests\n", ts_bin);
        ts_bin_state = -1;
        return 0;
    }
    fclose(f);
    ts_bin_state = 1;
    return 1;
}

static char *ts_slurp(const char *path)
{
    FILE *f = fopen(path, "rb");
    long n;
    char *buf;
    if (f == NULL) return NULL;
    fseek(f, 0, SEEK_END);
    n = ftell(f);
    fseek(f, 0, SEEK_SET);
    buf = (char *)malloc((size_t)n + 1);
    n = (long)fread(buf, 1, (size_t)n, f);
    buf[n] = 0;
    fclose(f);
    return buf;
}

/* run yafu with 'args' in ./yafu_test_snfs.tmp/<tag>/; returns its output
   (malloc'd, caller frees) and leaves the directory path in dir. */
static char *ts_run(const char *tag, const char *args, char *dir, size_t dirlen)
{
    static const char *junk[] = { "nfs.job", "nfs.fb", "nfs.log", "nfs.dat",
        "nfs.dat.p", "nfs.dat.ranges", "factor.log", "session.log",
        "factor.json", "siqs.dat", "YAFU_get_poly_score.out", "out.txt" };
    char path[4600], cmd[9000];
    size_t i;

    TS_MKDIR("yafu_test_snfs.tmp");
    snprintf(dir, dirlen, "yafu_test_snfs.tmp/%s", tag);
    TS_MKDIR(dir);

    /* stale nfs.job / siqs.dat files would make yafu resume old jobs */
    for (i = 0; i < sizeof junk / sizeof junk[0]; i++) {
        snprintf(path, sizeof path, "%s/%s", dir, junk[i]);
        remove(path);
    }

#if defined(_WIN32)
    /* system() runs "cmd /c <cmd>"; cmd strips the first and last quote of
       the line when it starts with one, so wrap the whole command in an
       extra pair of quotes. */
    for (i = 0; dir[i]; i++) if (dir[i] == '/') dir[i] = '\\';
    snprintf(cmd, sizeof cmd, "\"cd /d \"%s\" && \"%s\" %s > out.txt 2>&1\"", dir, ts_bin, args);
#else
    snprintf(cmd, sizeof cmd, "cd '%s' && '%s' %s > out.txt 2>&1", dir, ts_bin, args);
#endif
    fflush(stdout);
    if (system(cmd) == -1)
        return NULL;

    snprintf(path, sizeof path, "%s/out.txt", dir);
    return ts_slurp(path);
}

/* the method factor() settled on last ("fac: job type determined to be X") */
static int ts_job_type(const char *out, const char *which, int nth)
{
    const char *key = "job type determined to be ";
    const char *p = out;
    int k = 0;
    while ((p = strstr(p, key)) != NULL) {
        p += strlen(key);
        if (strncmp(p, which, strlen(which)) == 0) {
            if (++k >= nth)
                return 1;
        }
    }
    return 0;
}

static const char *ts_first_job_type(const char *out)
{
    const char *p = strstr(out, "job type determined to be ");
    return p ? p + strlen("job type determined to be ") : "";
}

/* the ECM target factor() set for an snfs job */
static double ts_snfs_ecm_target(const char *out, int *equiv)
{
    const char *p = strstr(out, "input has usable snfs form");
    const char *q;
    double from = 0, to = -1;
    if (p == NULL) return -1.0;
    q = p;
    while (q > out && strncmp(q, "ecm effort set from", 19) != 0) q--;
    if (sscanf(q, "ecm effort set from %lf to %lf", &from, &to) != 2)
        return -1.0;
    q = strstr(p, "equivalent gnfs size ");
    if (q != NULL && equiv != NULL)
        *equiv = atoi(q + strlen("equivalent gnfs size "));
    return to;
}

/* parse and validate dir/nfs.job against the full input */
static void ts_check_job(tk_ctx *tk, const char *dir, const char *expr, int ndigits,
    const char *type)
{
    char path[4600], *txt, *line;
    mpz_t N, n, Y0, Y1, c[9];
    int i, deg = -1, have_n = 0, have_y0 = 0, have_y1 = 0;
    char jtype[32] = "";

    snprintf(path, sizeof path, "%s/nfs.job", dir);
    txt = ts_slurp(path);
    TK_REQUIRE(tk, txt != NULL, "%s: no nfs.job written", expr);

    mpz_inits(N, n, Y0, Y1, NULL);
    for (i = 0; i < 9; i++) mpz_init_set_ui(c[i], 0);
    ts_eval(N, expr);

    for (line = strtok(txt, "\r\n"); line != NULL; line = strtok(NULL, "\r\n")) {
        if (!have_n && strncmp(line, "n: ", 3) == 0) { mpz_set_str(n, line + 3, 10); have_n = 1; }
        else if (strncmp(line, "Y0: ", 4) == 0) { mpz_set_str(Y0, line + 4, 10); have_y0 = 1; }
        else if (strncmp(line, "Y1: ", 4) == 0) { mpz_set_str(Y1, line + 4, 10); have_y1 = 1; }
        else if (strncmp(line, "type: ", 6) == 0) sscanf(line + 6, "%31s", jtype);
        else if (line[0] == 'c' && isdigit((unsigned char)line[1]) && line[2] == ':') {
            int k = line[1] - '0';
            if (k <= 8) {
                mpz_set_str(c[k], line + 4, 10);
                if (k > deg) deg = k;
            }
        }
    }

    TK_CHECKF(tk, have_n && mpz_cmp_ui(n, 1) > 0, "%s: nfs.job n missing or <= 1", expr);
    if (have_n && mpz_cmp_ui(n, 1) > 0) {
        TK_CHECKF(tk, mpz_divisible_p(N, n), "%s: nfs.job n does not divide the input", expr);
        if (ndigits != ANY)
            TK_CHECKF(tk, (int)mpz_sizeinbase(n, 10) == ndigits,
                "%s: nfs.job n has %d digits, expected %d", expr,
                (int)mpz_sizeinbase(n, 10), ndigits);
        TK_CHECKF(tk, have_y0 && have_y1 && deg >= 2, "%s: nfs.job is missing the polynomial", expr);
        if (have_y0 && have_y1 && deg >= 2)
            TK_CHECKF(tk, ts_common_root(n, c, deg, Y0, Y1),
                "%s: nfs.job f and g have no common root mod n", expr);
    }
    if (type != NULL)
        TK_CHECKF(tk, strcmp(jtype, type) == 0, "%s: nfs.job type '%s', expected '%s'",
            expr, jtype, type);

    for (i = 0; i < 9; i++) mpz_clear(c[i]);
    mpz_clears(N, n, Y0, Y1, NULL);
    free(txt);
}

#define BB_SNFS_C190 "2^688+1"
#define BB_C85 "4026441132641338206935076539055337224157349068413338560230768350534670023276796041989"
#define BB_C100 "6583474337077286035901994243868706192275025152834104495868770200652039120096298882866287304200021941"

/* the sieve method factor() actually started: -max_nfs 1 / -max_siqs 1 make
   it report the choice ("... job size larger than specified maximum") and
   stop instead of sieving.  Returns "nfs", "siqs" or "". */
static const char *ts_sieve(const char *out)
{
    const char *a = strstr(out, "NFS job size larger than specified maximum");
    const char *b = strstr(out, "SIQS job size larger than specified maximum");
    if (a && (!b || a < b)) return "nfs";
    if (b) return "siqs";
    return "";
}

/* the method decision only: -plan none jumps straight to the sieve choice,
   then -max_nfs 1 / -max_siqs 1 report it and stop (no ecm, no sieving).
   With pretest != NULL, "-pretest <pretest>" is used instead (ecm-only jobs
   stop before any sieve choice). */
static char *ts_decide(tk_ctx *tk, const char *tag, const char *input, const char *opts,
    const char *pretest, char *dir, size_t dirlen)
{
    char args[2048];
    char *out;
    if (pretest != NULL)
        snprintf(args, sizeof args, "\"factor(%s)\" -plan none -pretest %s %s", input, pretest, opts);
    else
        snprintf(args, sizeof args, "\"factor(%s)\" -plan none -max_nfs 1 -max_siqs 1 -v %s",
            input, opts);
    out = ts_run(tag, args, dir, dirlen);
    TK_REQUIRE(tk, out != NULL, "could not run yafu");
    if (tk_verbose(tk))
        printf("        %s: first job type '%.4s', sieve '%s'\n", tag, ts_first_job_type(out),
            ts_sieve(out));
    return out;
}

/* job type and sieve method both as expected */
static void ts_expect(tk_ctx *tk, const char *out, const char *what, const char *job,
    const char *notjob, const char *sieve)
{
    TK_CHECKF(tk, ts_job_type(out, job, 1), "%s: expected job type %s", what, job);
    if (notjob != NULL)
        TK_CHECKF(tk, !ts_job_type(out, notjob, 1), "%s: unexpected job type %s", what, notjob);
    if (sieve != NULL)
        TK_CHECKF(tk, strcmp(ts_sieve(out), sieve) == 0, "%s: started '%s', expected %s",
            what, ts_sieve(out), sieve);
}

static void t_cli_snfs_default(tk_ctx *tk)
{
    char dir[1024];
    char *out;
    int equiv = 0;
    double tgt;
    if (!ts_find_bin()) return;
    out = ts_decide(tk, "snfs_default", BB_SNFS_C190, "", NULL, dir, sizeof dir);
    ts_expect(tk, out, "c190 snfs", "snfs", "gnfs", "nfs");
    tgt = ts_snfs_ecm_target(out, &equiv);
    /* plan ratio 4/13 times the equivalent gnfs size of the snfs job */
    TK_CHECKF(tk, equiv > 0 && tgt > 0 && tgt > 4.0 / 13.0 * equiv - 0.01 &&
        tgt < 4.0 / 13.0 * equiv + 0.01,
        "snfs ecm target %1.2f does not match 4/13 * %d", tgt, equiv);
    free(out);
}

static void t_cli_gnfs_flag(tk_ctx *tk)
{
    char dir[1024];
    char *out;
    if (!ts_find_bin()) return;
    /* -gnfs on an snfs input: gnfs */
    out = ts_decide(tk, "gnfs_flag", BB_SNFS_C190, "-gnfs", NULL, dir, sizeof dir);
    ts_expect(tk, out, "c190 -gnfs", "gnfs", "snfs", "nfs");
    free(out);
    /* -gnfs below the qs/gnfs crossover: still gnfs (it used to run siqs) */
    out = ts_decide(tk, "gnfs_flag85", BB_C85, "-gnfs", NULL, dir, sizeof dir);
    ts_expect(tk, out, "c85 -gnfs", "gnfs", "siqs", "nfs");
    free(out);
    /* and without -gnfs that c85 is a siqs job */
    out = ts_decide(tk, "c85_dflt", BB_C85, "", NULL, dir, sizeof dir);
    ts_expect(tk, out, "c85", "siqs", "gnfs", "siqs");
    free(out);
}

static void t_cli_snfs_xover(tk_ctx *tk)
{
    char dir[1024];
    char *out;
    if (!ts_find_bin()) return;
    /* the c190 is below a 200-digit qs/snfs crossover: siqs, and it must
       actually start siqs (the scheduler used to start nfs regardless) */
    out = ts_decide(tk, "snfs_xover", BB_SNFS_C190, "-snfs_xover 200", NULL, dir, sizeof dir);
    ts_expect(tk, out, "c190 -snfs_xover 200", "siqs", "snfs", "siqs");
    free(out);
}

static void t_cli_gnfs_xover(tk_ctx *tk)
{
    char dir[1024];
    char *out;
    if (!ts_find_bin()) return;
    out = ts_decide(tk, "xover_dflt", BB_C100, "", NULL, dir, sizeof dir);
    ts_expect(tk, out, "c100", "gnfs", "siqs", "nfs");
    free(out);
    out = ts_decide(tk, "xover_105", BB_C100, "-xover 105", NULL, dir, sizeof dir);
    ts_expect(tk, out, "c100 -xover 105", "siqs", "gnfs", "siqs");
    free(out);
}

static void t_cli_pretest(tk_ctx *tk)
{
    char dir[1024];
    char *out;
    if (!ts_find_bin()) return;
    out = ts_decide(tk, "pretest30", BB_SNFS_C190, "", "30", dir, sizeof dir);
    TK_CHECKF(tk, ts_job_type(out, "ecm", 1) && !ts_job_type(out, "snfs", 1),
        "-pretest 30: expected job type ecm");
    TK_CHECKF(tk, strstr(out, "input has usable snfs form") == NULL ||
        strstr(out, "ecm effort maintained at 30.00") != NULL,
        "-pretest 30: ecm target was changed by the snfs form");
    free(out);
}

static void t_cli_factor_np(tk_ctx *tk)
{
    char dir[1024];
    char *out;
    if (!ts_find_bin()) return;
    out = ts_run("factor_np", "\"factor(" BB_SNFS_C190 ")\" -plan none -np -testsieve 99999",
        dir, sizeof dir);
    TK_REQUIRE(tk, out != NULL, "could not run yafu");
    TK_CHECKF(tk, ts_job_type(out, "snfs", 1), "expected job type snfs");
    ts_check_job(tk, dir, BB_SNFS_C190, 190, "snfs");
    free(out);
}

static void t_cli_snfs_np(tk_ctx *tk)
{
    char dir[1024];
    char *out;
    if (!ts_find_bin()) return;
    out = ts_run("snfs_np", "\"snfs(13^199-10^199,13^199-10^199)\" -np -testsieve 99999",
        dir, sizeof dir);
    TK_REQUIRE(tk, out != NULL, "could not run yafu");
    ts_check_job(tk, dir, "13^199-10^199", ANY, "snfs");
    free(out);
}

static void t_cli_ecm_then_siqs(tk_ctx *tk)
{
    char dir[1024], path[1100];
    char *out, *log;
    if (!ts_find_bin()) return;
    /* 2^401-1: the c110 left after trial division / P-1 is an snfs job.  ECM
       curve 0 at B1=11000 with sigma 938848200 finds the prp22; the c88 left is
       a siqs job, which -max_siqs 85 declines, so nothing is sieved. */
    out = ts_run("ecm_sigma", "\"factor(2^401-1)\" -sigma 938848200 -threads 1 "
        "-plan light -max_siqs 85", dir, sizeof dir);
    TK_REQUIRE(tk, out != NULL, "could not run yafu");
    TK_CHECKF(tk, ts_job_type(out, "snfs", 1), "expected an snfs job before ecm");
    TK_CHECKF(tk, strstr(out, "found prp22 factor = 2136958965524920285681") != NULL,
        "ecm did not find the prp22");
    TK_CHECKF(tk, ts_job_type(out, "siqs", 1), "expected siqs for the remaining c88");
    TK_CHECKF(tk, strstr(out, "C88 = 2820107517997538057112522872095070422330922339778101711"
        "779609977566650816255226156185129") != NULL, "expected the c88 cofactor to remain");
    snprintf(path, sizeof path, "%s/factor.log", dir);
    log = ts_slurp(path);
    TK_CHECKF(tk, log != NULL && strstr(log, "sigma=938848200") != NULL,
        "factor.log does not record the ecm find with sigma 938848200");
    free(log);
    free(out);
}

static const tk_test tk__snfs_tests[] = {
    { "brent_2_401m",   t_case_0,  "fast snfs integration" },
    { "brent_2_613m",   t_case_1,  "fast snfs integration" },
    { "brent_2_1061m",  t_case_2,  "fast snfs integration" },
    { "brent_coeff",    t_case_3,  "fast snfs integration" },
    { "brent_coeff_b5", t_case_4,  "fast snfs integration" },
    { "brent_cbase_c2", t_case_5,  "fast snfs integration" },
    { "brent_big_b5",   t_case_6,  "fast snfs integration" },
    { "brent_big_b7",   t_case_7,  "fast snfs integration" },
    { "brent_e1155p",   t_case_8,  "fast snfs integration" },
    { "brent_e1296m",   t_case_9,  "fast snfs integration" },
    { "brent_e1208m",   t_case_10, "fast snfs integration" },
    { "brent_10_222m",  t_case_11, "fast snfs integration" },
    { "brent_10_126m",  t_case_12, "fast snfs integration" },
    { "brent_mixed",    t_case_13, "fast snfs integration" },
    { "hcunn_13_10",    t_case_14, "fast snfs integration" },
    { "hcunn_50_49e",   t_case_15, "fast snfs integration" },
    { "aurif_2_398p",   t_case_16, "fast snfs integration" },
    { "aurif_2_990p",   t_case_17, "fast snfs integration" },
    { "aurif_2_1086p",  t_case_18, "fast snfs integration" },
    { "aurif_h_50_94",  t_case_19, "fast snfs integration" },
    { "aurif_h_50_142", t_case_20, "fast snfs integration" },
    { "aurif_h_28_9",   t_case_21, "fast snfs integration" },
    { "xyyxf_84_83",    t_case_22, "fast snfs integration" },
    { "direct_16_19",   t_case_23, "fast snfs integration" },
    { "pow_diff_sq",    t_case_24, "fast snfs integration" },
    { "pow_diff_cube",  t_case_25, "fast snfs integration" },
    { "cli_snfs_default", t_cli_snfs_default, "snfs integration cli" },
    { "cli_gnfs_flag",    t_cli_gnfs_flag,    "snfs integration cli" },
    { "cli_snfs_xover",   t_cli_snfs_xover,   "snfs integration cli" },
    { "cli_gnfs_xover",   t_cli_gnfs_xover,   "snfs integration cli" },
    { "cli_pretest",      t_cli_pretest,      "snfs integration cli" },
    { "cli_factor_np",    t_cli_factor_np,    "snfs integration cli" },
    { "cli_snfs_np",      t_cli_snfs_np,      "snfs integration cli" },
    { "cli_ecm_then_siqs", t_cli_ecm_then_siqs, "snfs integration cli" },
};

const tk_module tk_module_snfs = {
    "snfs",
    "SNFS form detection, algebraic factors, polys, and factor() method choice",
    tk__snfs_tests,
    (int)(sizeof tk__snfs_tests / sizeof tk__snfs_tests[0])
};
