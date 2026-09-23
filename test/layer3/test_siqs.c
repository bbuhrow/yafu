/*----------------------------------------------------------------------
 Layer 3 -- individual factoring routine: SIQS (self-initializing QS).

 Runs SIQS in-process on the 50-, 60-, and 70-digit composites from
 siqsbench and verifies the result is a complete, correct factorization:
 the product of the returned factors^count equals N and every returned
 factor is prime (GMP-verified). SIQS is a complete algorithm, so this is
 strict pass/fail -- no probabilistic miss budget.

 Each size is its own test, so the runner's per-test wall-clock column
 reports the time for that number; the test also prints a one-line
 proof-of-life detail (the split it found).

 In-process recipe (from calc.c's siqs case): SIQS reads its input from
 fobj->qs_obj.gmp_n (NOT fobj->N), so we set that before calling.
   fact_obj_t *fobj = new_default_factorization(N);
   fobj->VFLAG = -1;  fobj->LOGFLAG = 0;
   mpz_set(fobj->qs_obj.gmp_n, N);
   SIQS(fobj);
   // factors in fobj->factors->factors[i].{factor,count}

 SIQS expects an input with no small factors and that is not prime/perfect
 power (these balanced semiprimes satisfy that). It may drop a relations
 savefile (siqs.dat) in the working directory -- a side effect, not a
 failure. Threads default to whatever new_default_factorization sets (1);
 fobj->THREADS could be raised for a multi-threaded timing.
 Public domain.
----------------------------------------------------------------------*/
#include "testkit.h"
#include <gmp.h>
#include <stdio.h>
#include "factor.h"
#include "qs.h"

static void siqs_case(tk_ctx *tk, const char *decimal)
{
    mpz_t N, prod, pe;
    fact_obj_t *fobj;
    yfactor_list_t *fl;
    int i;

    mpz_inits(N, prod, pe, NULL);
    TK_REQUIRE(tk, mpz_set_str(N, decimal, 10) == 0, "bad decimal literal");

    fobj = new_default_factorization(N);
    TK_REQUIRE(tk, fobj != NULL, "new_default_factorization returned NULL");
    fobj->VFLAG = -1;      /* silent */
    fobj->LOGFLAG = 0;     /* no logfile */
    mpz_set(fobj->qs_obj.gmp_n, N);    /* SIQS reads its input from here */

    SIQS(fobj);
    fl = fobj->factors;

    /* completeness + correctness: product of factors^count == N, each prime */
    mpz_set_ui(prod, 1);
    for (i = 0; i < fl->num_factors; i++) {
        int cnt = fl->factors[i].count;
        TK_CHECKF(tk, mpz_probab_prime_p(fl->factors[i].factor, 30) > 0,
                  "SIQS returned a non-prime factor");
        TK_CHECKF(tk, cnt >= 1, "factor with non-positive multiplicity");
        mpz_pow_ui(pe, fl->factors[i].factor, (unsigned long)cnt);
        mpz_mul(prod, prod, pe);
    }
    TK_CHECKF(tk, fl->num_factors >= 2, "SIQS did not split the composite");
    TK_CHECKF(tk, mpz_cmp(prod, N) == 0, "product of returned factors != N");

    /* proof-of-life detail (the split found) */
    printf("        %zu-digit N ->", mpz_sizeinbase(N, 10));
    for (i = 0; i < fl->num_factors; i++)
        printf(" %s%zu-digit", i ? "x " : "",
               mpz_sizeinbase(fl->factors[i].factor, 10));
    printf("  (%d factor%s)\n", fl->num_factors, fl->num_factors == 1 ? "" : "s");
    fflush(stdout);

    free_factobj(fobj);
    mpz_clears(N, prod, pe, NULL);
}

static void t_siqs_c50(tk_ctx *tk)
{ siqs_case(tk, "29660734457033883936073030405220515257819037444591"); }

static void t_siqs_c60(tk_ctx *tk)
{ siqs_case(tk, "349594255864176572614071853194924838158088864370890996447417"); }

static void t_siqs_c70(tk_ctx *tk)
{ siqs_case(tk, "6470287906463336878241474855987746904297564226439499503918586590778209"); }

static const tk_test tk__siqs_tests[] = {
    { "siqs_c50", t_siqs_c50, "slow siqs integration" },
    { "siqs_c60", t_siqs_c60, "slow siqs integration" },
    { "siqs_c70", t_siqs_c70, "slow siqs integration" }
};

const tk_module tk_module_siqs = {
    "siqs",
    "SIQS on the 50/60/70-digit siqsbench composites",
    tk__siqs_tests,
    (int)(sizeof tk__siqs_tests / sizeof tk__siqs_tests[0])
};
