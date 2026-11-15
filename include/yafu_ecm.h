/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

Some parts of the code (and also this header), included in this 
distribution have been reused from other sources. In particular I 
have benefitted greatly from the work of Jason Papadopoulos's msieve @ 
www.boo.net/~jasonp, Scott Contini's mpqs implementation, and Tom St. 
Denis Tom's Fast Math library.  Many thanks to their kind donation of 
code to the public domain.
       				   --bbuhrow@gmail.com 11/24/09
----------------------------------------------------------------------*/
#ifndef _YECM_H_
#define _YECM_H_

#include "factor.h"
#include "ytools.h"
#include "arith.h"
#include <stdio.h>
#include <gmp_u64_xface.h>


/* The standard t-levels and their B1 values.  The expected curves per
   t-level are no longer tabulated: they come from the probability model in
   ecm_prob.h for the parameters the curves actually run with (GMP-ECM
   param 0/1/3 with its stage 2, AVX-ECM Suyama or Edwards curves with the
   FFT or pairing stage 2 and Brent-Suyama); docs/ecm_prob.md has the
   tables it gives and the comparison with GMP-ECM's.

                                       2k   11k   50k   250k   1M   3M   11M   43M   110M  260M  850M  2900M */

UNUSED_VAR static int ecm_levels[NUM_ECM_LEVELS] = {
    15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70 };

UNUSED_VAR static uint64_t ecm_std_b1[NUM_ECM_LEVELS] = { 2000, 11000, 50000, 250000, 1000000,
    3000000, 11000000, 43000000, 110000000, 260000000, 850000000, 2900000000 };



/* ============================ interface to gmpecm/pm1/pp1 ============================ */
extern void pollard_loop(fact_obj_t* fobj);
extern void williams_loop(fact_obj_t* fobj);
extern int ecm_loop(fact_obj_t* fobj);


/* ============================ interface to avxecm ============================ */
extern void vec_ecm_main(fact_obj_t* fobj, uint32_t numcurves, uint64_t B1,
    uint64_t B2, int threads, int* numfactors, int verbose,
    int save_b1, uint32_t* curves_run);
extern void vecPP1(fact_obj_t* fobj);
// 1 if the AVX-ECM RNS backend is built in (factor/avx-ecm/avxecm_rns.c)
extern int avxecm_rns_available(void);
// RNS size-class tables for N (NULL if RNS cannot handle it); free with
// avxecm_rns_free_fixed.  *words = vector words the RNS layout needs.
extern void *avxecm_rns_fixed(mpz_t N, uint32_t *words);
extern void avxecm_rns_free_fixed(void *fx);
// the FFT stage 2 avxecm_s2fft_setup would build for an nbits-bit input
// (avxecm_stage2fft.h): 0 if unavailable, else the covered B2, the
// Brent-Suyama polynomial and the stage 2 pairs
extern int avxecm_s2fft_plan(uint32_t nbits, uint64_t B1, uint64_t B2, int threads,
    double mem_mb, int auto_b2, int bs, int edw_stage1,
    uint64_t *B2out, int *bsout, double *pairs);
extern void vecPM1(fact_obj_t* fobj);

/* ============================ interface to utilities and table data ============================ */
// the curve parameters (family, B2, Brent-Suyama, pairs) that curves at
// b1 would run with under the current options and input
extern void ecm_plan_curve_params(fact_obj_t* fobj, uint64_t b1, ecm_curve_param_t* cp);
// curves of kind cp for one full t-level at ecm_levels[tlevel] digits
extern uint32_t get_curves_for_tlevel(int tlevel, const ecm_curve_param_t* cp);
// curves of kind cp that bring the completed work to target_tlevel
extern uint32_t get_curves_required(ecm_obj_t* ecm_obj, double target_tlevel, const ecm_curve_param_t* cp);
extern void print_std_ecm_work_done(ecm_obj_t* ecm_obj, int disp_levels, FILE* log, int VFLAG, int LOGFLAG);
// add completed curves of kind cp to the work record and t-levels
extern void record_curves_completed(ecm_obj_t* ecm_obj, int curves, const ecm_curve_param_t* cp);
int print_B1B2(fact_obj_t* fobj, FILE* fid);

#endif // #ifndef _YECM_H_
