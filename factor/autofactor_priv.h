/*
 * ========== Note: created by Claude with prompting by bbuhrow ==========
 * factor()'s work-state machine, shared by autofactor.c and bayes_sim.c
 * (the -plan replay simulator, docs/ecm_bayes.md section 6).  Not a public
 * interface.
 */

#ifndef AUTOFACTOR_PRIV_H
#define AUTOFACTOR_PRIV_H

#include <stdint.h>
#include "factor.h"

enum job_type_e {
	job_snfs,
	job_gnfs,
	job_siqs,
	job_ecm,
	job_unknown
};

typedef struct
{
	// total effort so far
	double total_time;
	double qs_time;
	double nfs_time;
	double trialdiv_time;
	double fermat_time;
	double rho_time;
	double pp1_time;
	double pm1_time;
	double ecm_time;
	double pp1_lvl1_time_per_curve;
	double pp1_lvl2_time_per_curve;
	double pp1_lvl3_time_per_curve;
	double pm1_lvl1_time_per_curve;
	double pm1_lvl2_time_per_curve;
	double pm1_lvl3_time_per_curve;
	double ecm_15digit_time_per_curve;
	double ecm_20digit_time_per_curve;
	double ecm_25digit_time_per_curve;
	double ecm_30digit_time_per_curve;
	double ecm_35digit_time_per_curve;
	double ecm_40digit_time_per_curve;
	double ecm_45digit_time_per_curve;
	double ecm_50digit_time_per_curve;
	double ecm_55digit_time_per_curve;
	double ecm_60digit_time_per_curve;
	double ecm_65digit_time_per_curve;
    double initial_work;

	// amount of work we've done in various areas
	uint32_t  tdiv_limit;
	uint32_t  fermat_iterations;
	uint32_t  rho_iterations;
	uint32_t  rho_bases;
	uint32_t  pp1_lvl1_curves;
	uint32_t  pm1_lvl1_curves;
	uint32_t  pp1_lvl2_curves;
	uint32_t  pm1_lvl2_curves;
	uint32_t  pp1_lvl3_curves;
	uint32_t  pm1_lvl3_curves;
	uint32_t  ecm_15digit_curves;
	uint32_t  ecm_20digit_curves;
	uint32_t  ecm_25digit_curves;
	uint32_t  ecm_30digit_curves;
	uint32_t  ecm_35digit_curves;
	uint32_t  ecm_40digit_curves;
	uint32_t  ecm_45digit_curves;
	uint32_t  ecm_50digit_curves;
	uint32_t  ecm_55digit_curves;
	uint32_t  ecm_60digit_curves;
	uint32_t  ecm_65digit_curves;
	int min_pretest_done;

	// max amount of work we'll allow in various areas.
	// to be filled in during init, or overriden by user
	uint32_t  tdiv_max_limit;
	uint32_t  fermat_max_iterations;
	uint32_t  rho_max_iterations;
	uint32_t  rho_max_bases;
	uint32_t  pp1_max_lvl1_curves;
	uint32_t  pm1_max_lvl1_curves;
	uint32_t  pp1_max_lvl2_curves;
	uint32_t  pm1_max_lvl2_curves;
	uint32_t  pp1_max_lvl3_curves;
	uint32_t  pm1_max_lvl3_curves;
	uint32_t  ecm_max_15digit_curves;
	uint32_t  ecm_max_20digit_curves;
	uint32_t  ecm_max_25digit_curves;
	uint32_t  ecm_max_30digit_curves;
	uint32_t  ecm_max_35digit_curves;
	uint32_t  ecm_max_40digit_curves;
	uint32_t  ecm_max_45digit_curves;
	uint32_t  ecm_max_50digit_curves;
	uint32_t  ecm_max_55digit_curves;
	uint32_t  ecm_max_60digit_curves;
	uint32_t  ecm_max_65digit_curves;

	// current parameters
	uint32_t  B1;
	uint64_t  B2;
	uint32_t  curves;

	// target job type
	enum job_type_e target_job_type;

	// equivalent gnfs size of the snfs job, when target_job_type == job_snfs
	int snfs_equiv_size;

} factor_work_t;

enum factorization_state {
	state_idle,
	state_trialdiv,
	state_fermat,
	state_rho,
	state_pp1_lvl1,
	state_pm1_lvl1,
	state_pp1_lvl2,
	state_pm1_lvl2,
	state_pp1_lvl3,
	state_pm1_lvl3,
	state_ecm_15digit,
	state_ecm_20digit,
	state_ecm_25digit,
	state_ecm_30digit,
	state_ecm_35digit,
	state_ecm_40digit,
	state_ecm_45digit,
	state_ecm_50digit,
	state_ecm_55digit,
	state_ecm_60digit,
	state_ecm_65digit,
	state_qs,
	state_nfs,
	state_done
};

enum factorization_state schedule_work(factor_work_t *fwork, mpz_t b, fact_obj_t *fobj);
uint32_t  get_ecm_curves_done(factor_work_t *fwork, enum factorization_state state);
uint32_t  set_ecm_curves_done(factor_work_t *fwork, enum factorization_state state, uint32_t  curves_done);
void init_factor_work(factor_work_t *fwork, fact_obj_t *fobj);

// Replacement cost estimates for -plan bayesian.  When bayes_cost_hook is
// not NULL, bayes_schedule() takes the sieve time and the ECM seconds per
// curve from it instead of the tune data (used by the bayessim() replay to
// give the rule the true costs, or deliberately wrong ones).  Change gen
// whenever the costs change: it invalidates the cached step options.
typedef struct
{
	double (*sieve_sec)(void *ctx, double digits);
	double (*ecm_sec_per_curve)(void *ctx, double digits, double B1);
	void *ctx;
	int gen;
} bayes_cost_hook_t;

extern bayes_cost_hook_t *bayes_cost_hook;

#endif
