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
       				   --bbuhrow@gmail.com 3/26/10
----------------------------------------------------------------------*/

#include "yafu.h"
#include "soe.h"
#include "factor.h"
#include "qs.h"
#include "nfs.h"
#include "yafu_ecm.h"
#include "util.h"
#include "yafu_string.h"
#include "mpz_aprcl.h"

/* produced using ecm -v -v -v for the various B1 bounds (default B2).
/	Thanks A. Schindel !
/
/					2k			11k			50k			250k		1M			3M			11M			43M			110M	260M	850M */
#define NUM_ECM_LEVELS 11
static int ecm_levels[NUM_ECM_LEVELS] = {
	15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65};
double ecm_data[NUM_ECM_LEVELS][NUM_ECM_LEVELS] = {	
	/*t15, 2000,	*/	{30,		12,			7,			5,			3,			2,			2,			2,			2,		1,		1},
	/*t20, 11000,	*/	{844,		74,			21,			8,			5,			3,			2,			2,			2,		2,		1},
	/*t25, 50000,	*/	{58129,		1539,		214,		50,			20,			11,			7,			5,			4,		3,		3},
	/*t30, 250000,	*/	{6711967,	49962,		3288,		430,		118,		54,			26,			14,			10,		8,		6},
	/*t35, 1E+06,	*/	{1.20E+09,	2292278,	68422,		4914,		904,		322,		122,		54,			34,		23,		15},
	/*t40, 3E+06,	*/	{2.90E+12,	1.40E+08,	1849287,	70293,		8613,		2350,		681,		242,		135,	82,		47},
	/*t45, 11E+06,	*/	{9.00E+99,	1.10E+10,	6.10E+07,	1214949,	97057,		20265,		4480,		1263,		613,	333,	168},
	/*t50, 44E+06,	*/	{9.00E+99,	9.00E+99,	2.50E+09,	2.50E+07,	1270662,	199745,		33652,		7404,		3133,	1512,	661},
	/*t55, 110E+06,	*/	{9.00E+99,	9.00E+99,	1.30E+11,	5.90E+08,	1.90E+07,	2246256,	283939,		48714,		17769,	7643,	2865},
	/*t60, 260E+06,	*/	{9.00E+99,	9.00E+99,	5.80E+16,	1.60E+10,	3.20E+08,	2.80E+07,	2655154,	350439,		111196,	42017,	13611},
	/*t65, 850E+06,	*/	{9.00E+99,	9.00E+99,	8.20E+21,	2.70E+13,	6.10E+09,	4.00E+08,	2.70E+07,	2768535,	751771,	250214,	69408}};



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
	uint32 tdiv_limit;
	uint32 fermat_iterations;
	uint32 rho_iterations;
	uint32 rho_bases;
	uint32 pp1_lvl1_curves;
	uint32 pm1_lvl1_curves;
	uint32 pp1_lvl2_curves;
	uint32 pm1_lvl2_curves;
	uint32 pp1_lvl3_curves;
	uint32 pm1_lvl3_curves;
	uint32 ecm_15digit_curves;
	uint32 ecm_20digit_curves;
	uint32 ecm_25digit_curves;
	uint32 ecm_30digit_curves;
	uint32 ecm_35digit_curves;
	uint32 ecm_40digit_curves;
	uint32 ecm_45digit_curves;
	uint32 ecm_50digit_curves;
	uint32 ecm_55digit_curves;
	uint32 ecm_60digit_curves;
	uint32 ecm_65digit_curves;
	int min_pretest_done;

	// max amount of work we'll allow in various areas.
	// to be filled in during init, or overriden by user
	uint32 tdiv_max_limit;
	uint32 fermat_max_iterations;
	uint32 rho_max_iterations;
	uint32 rho_max_bases;
	uint32 pp1_max_lvl1_curves;
	uint32 pm1_max_lvl1_curves;
	uint32 pp1_max_lvl2_curves;
	uint32 pm1_max_lvl2_curves;
	uint32 pp1_max_lvl3_curves;
	uint32 pm1_max_lvl3_curves;
	uint32 ecm_max_15digit_curves;
	uint32 ecm_max_20digit_curves;
	uint32 ecm_max_25digit_curves;
	uint32 ecm_max_30digit_curves;
	uint32 ecm_max_35digit_curves;
	uint32 ecm_max_40digit_curves;
	uint32 ecm_max_45digit_curves;
	uint32 ecm_max_50digit_curves;	
	uint32 ecm_max_55digit_curves;	
	uint32 ecm_max_60digit_curves;	
	uint32 ecm_max_65digit_curves;	

	// current parameters
	uint32 B1;
	uint64 B2;
	uint32 curves;

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

// local functions to do state based factorization
double get_qs_time_estimate(fact_obj_t *fobj, mpz_t b);
double get_gnfs_time_estimate(fact_obj_t *fobj, mpz_t b);
void do_work(enum factorization_state method, factor_work_t *fwork, mpz_t b, fact_obj_t *fobj);
enum factorization_state schedule_work(factor_work_t *fwork, mpz_t b, fact_obj_t *fobj);
int check_if_done(fact_obj_t *fobj, mpz_t N);
uint32 get_ecm_curves_done(factor_work_t *fwork, enum factorization_state state);
uint32 set_ecm_curves_done(factor_work_t *fwork, enum factorization_state state, uint32 curves_done);
uint32 get_max_ecm_curves(factor_work_t *fwork, enum factorization_state state);
void set_work_params(factor_work_t *fwork, enum factorization_state state);
int check_tune_params(fact_obj_t *fobj);
enum factorization_state get_next_state(factor_work_t *fwork, fact_obj_t *fobj);
double compute_ecm_work_done(factor_work_t *fwork, int disp, FILE *log);
void init_factor_work(factor_work_t *fwork, fact_obj_t *fobj);
void interp_and_set_curves(factor_work_t *fwork, fact_obj_t *fobj, 
	enum factorization_state state, double work_done,
	double target_digits, int log_results);


void init_factobj(fact_obj_t *fobj)
{
	// get space for everything
	alloc_factobj(fobj);

	// initialize global stuff in fobj
	fobj->seed1 = g_rand.low;
	fobj->seed2 = g_rand.hi;
	yafu_get_cache_sizes(&fobj->cache_size1,&fobj->cache_size2);
	fobj->flags = 0;
	fobj->num_threads = 1;		//read from input arguments	
	strcpy(fobj->flogname,"factor.log");	
	fobj->do_logging = 1;

	// initialize stuff for rho	
	fobj->rho_obj.iterations = 1000;		
	fobj->rho_obj.curr_poly = 0;	

	// initialize stuff for pm1	
	fobj->pm1_obj.B1 = 100000;
	fobj->pm1_obj.B2 = 10000000;
	fobj->pm1_obj.stg2_is_default = 1;		
	fobj->pm1_obj.pm1_exponent = 0;
	fobj->pm1_obj.pm1_multiplier = 0;
	fobj->pm1_obj.pm1_tune_freq = 0;

	// initialize stuff for pp1	
	fobj->pp1_obj.B1 = 20000;
	fobj->pp1_obj.B2 = 1000000;
	fobj->pp1_obj.stg2_is_default = 1;	
	fobj->pp1_obj.pp1_exponent = 0;
	fobj->pp1_obj.pp1_multiplier = 0;
	fobj->pp1_obj.pp1_tune_freq = 0;

	// initialize stuff for ecm	
	fobj->ecm_obj.B1 = 11000;
	fobj->ecm_obj.B2 = 1100000;
	fobj->ecm_obj.stg2_is_default = 1;
	fobj->ecm_obj.sigma = 0;
	fobj->ecm_obj.num_curves = 90;
	fobj->ecm_obj.curves_run = 0;
	fobj->ecm_obj.ecm_exponent = 0;
	fobj->ecm_obj.ecm_multiplier = 0;
	fobj->ecm_obj.ecm_tune_freq = 0;
    fobj->ecm_obj.bail_on_factor = 1;
    fobj->ecm_obj.save_b1 = 0;
    

	// unlike ggnfs, ecm does not *require* external binaries.  
	// an empty string indicates the use of the built-in GMP-ECM hooks, while
	// a non-empty string (filled in by the user) will indicate the use of
	// an external binary
	strcpy(fobj->ecm_obj.ecm_path,"");
	fobj->ecm_obj.use_external = 0;
#ifdef USE_AVX512F
    fobj->ecm_obj.prefer_gmpecm = 0;
    fobj->ecm_obj.ecm_ext_xover = 40000000;
#else
    fobj->ecm_obj.prefer_gmpecm = 1;
	fobj->ecm_obj.ecm_ext_xover = 48000;
#endif

	// initialize stuff for squfof
	fobj->squfof_obj.num_factors = 0;	

	// initialize stuff for qs	
	fobj->qs_obj.gbl_override_B_flag = 0;
	fobj->qs_obj.gbl_override_B = 0;
    fobj->qs_obj.gbl_override_small_cutoff_flag = 0;
    fobj->qs_obj.gbl_override_small_cutoff = 0;
	fobj->qs_obj.gbl_override_blocks_flag = 0;
	fobj->qs_obj.gbl_override_blocks = 0 ;
	fobj->qs_obj.gbl_override_lpmult_flag = 0;
	fobj->qs_obj.gbl_override_lpmult = 0;
	fobj->qs_obj.gbl_override_rel_flag = 0;
	fobj->qs_obj.gbl_override_rel = 0;
	fobj->qs_obj.gbl_override_tf_flag = 0;
	fobj->qs_obj.gbl_override_tf = 0;
	fobj->qs_obj.gbl_override_time_flag = 0;
	fobj->qs_obj.gbl_override_time = 0;
	fobj->qs_obj.gbl_override_mfbd = 0.;
	fobj->qs_obj.gbl_override_mfbt = 0.;
	fobj->qs_obj.gbl_override_lpb = 0;
    fobj->qs_obj.gbl_override_bdiv_flag = 0;
    fobj->qs_obj.gbl_override_bdiv = 3;
    fobj->qs_obj.gbl_btarget = 500000;
	fobj->qs_obj.flags = 0;
	fobj->qs_obj.gbl_force_DLP = 0;
	fobj->qs_obj.gbl_force_TLP = 0;
	fobj->qs_obj.qs_exponent = 0;
	fobj->qs_obj.qs_multiplier = 0;
	fobj->qs_obj.qs_tune_freq = 0;
	fobj->qs_obj.no_small_cutoff_opt = 0;
	strcpy(fobj->qs_obj.siqs_savefile,"siqs.dat");
	init_lehman();

	// initialize stuff for trial division	
	fobj->div_obj.print = 0;	
	fobj->div_obj.limit = 10000;
	fobj->div_obj.fmtlimit = 1000000;
	
	//initialize stuff for nfs
	fobj->nfs_obj.snfs = 0;
	fobj->nfs_obj.gnfs = 0;
	fobj->nfs_obj.gnfs_exponent = 0;
	fobj->nfs_obj.gnfs_multiplier = 0;
	fobj->nfs_obj.gnfs_tune_freq = 0;
	fobj->nfs_obj.min_digits = 85;
	fobj->nfs_obj.filter_min_rels_nudge = 1.05;	// raise min_rels bounds by 5%
													// on unsuccessful filtering
	fobj->nfs_obj.siever = 0;					//default, use automatic selection
	fobj->nfs_obj.startq = 0;					//default, not used
	fobj->nfs_obj.rangeq = 0;					//default, not used
	fobj->nfs_obj.polystart = 0;					//default, not used
	fobj->nfs_obj.polyrange = 0;					//default, not used
	strcpy(fobj->nfs_obj.outputfile,"nfs.dat");			//default
	strcpy(fobj->nfs_obj.logfile,"nfs.log");			//default
	strcpy(fobj->nfs_obj.fbfile,"nfs.fb");				//default
	fobj->nfs_obj.sq_side = 0;					//default = algebraic
	fobj->nfs_obj.timeout = 0;					//default, not used
	strcpy(fobj->nfs_obj.job_infile,"nfs.job");			//default
	fobj->nfs_obj.poly_option = 0;					//default = fast search
										//1 = wide
									//2 = deep
	fobj->nfs_obj.restart_flag = 0;					//default = not a restart
	fobj->nfs_obj.nfs_phases = NFS_DEFAULT_PHASES;
	fobj->nfs_obj.snfs_testsieve_threshold = 160;
	*fobj->nfs_obj.filearg = '\0';

	fobj->nfs_obj.polybatch = 250;						//default	
#if defined(_WIN64)
	strcpy(fobj->nfs_obj.ggnfs_dir,".\\");
#elif defined(WIN32)
	strcpy(fobj->nfs_obj.ggnfs_dir,".\\");
#else
	strcpy(fobj->nfs_obj.ggnfs_dir,"./");
#endif

	//initialize autofactor object
	//whether we want to output certain info to their own files...
	fobj->autofact_obj.want_output_primes = 0;
	fobj->autofact_obj.want_output_factors = 0;
	fobj->autofact_obj.want_output_unfactored = 0;
	fobj->autofact_obj.want_output_expressions = 1;
	fobj->autofact_obj.qs_gnfs_xover = 95;
    fobj->autofact_obj.qs_snfs_xover = 75;
	// use xover even when timing info is available
	fobj->autofact_obj.prefer_xover = 0;			
	fobj->autofact_obj.want_only_1_factor = 0;
	fobj->autofact_obj.no_ecm = 0;
	fobj->autofact_obj.target_pretest_ratio = 4.0 / 13.0;
	fobj->autofact_obj.initial_work = 0.0;
	fobj->autofact_obj.has_snfs_form = -1;		// not checked yet

	//pretesting plan used by factor()
	fobj->autofact_obj.yafu_pretest_plan = PRETEST_NORMAL;
	strcpy(fobj->autofact_obj.plan_str,"normal");
	fobj->autofact_obj.only_pretest = 0;
	fobj->autofact_obj.autofact_active = 0;

	// if a number is <= aprcl_prove_cutoff, we will prove it prime or composite
	fobj->aprcl_prove_cutoff = 500;
	// if a number is >= aprcl_display_cutoff, we will show the APRCL progress
	fobj->aprcl_display_cutoff = 200;

	return;
}

void free_factobj(fact_obj_t *fobj)
{
	uint32 i;

	// free stuff in rho
	free(fobj->rho_obj.polynomials);
	mpz_clear(fobj->rho_obj.gmp_n);
	mpz_clear(fobj->rho_obj.gmp_f);

	// free stuff in pm1
	mpz_clear(fobj->pm1_obj.gmp_n);
	mpz_clear(fobj->pm1_obj.gmp_f);

	// free stuff in pp1
	mpz_clear(fobj->pp1_obj.gmp_n);
	mpz_clear(fobj->pp1_obj.gmp_f);

	// free any factors found in ecm
	for (i=0; i<fobj->ecm_obj.num_factors; i++)
		zFree(&fobj->ecm_obj.factors[i]);
	free(fobj->ecm_obj.factors);

	// then free other stuff in ecm
	mpz_clear(fobj->ecm_obj.gmp_n);
	mpz_clear(fobj->ecm_obj.gmp_f);

	// free any factors found in squfof
	for (i=0; i<fobj->squfof_obj.num_factors; i++)
		zFree(&fobj->squfof_obj.factors[i]);
	free(fobj->squfof_obj.factors);

	// then free other stuff in squfof
	mpz_clear(fobj->squfof_obj.gmp_n);
	mpz_clear(fobj->squfof_obj.gmp_f);

	// free any factors found in qs
	for (i=0; i<fobj->qs_obj.num_factors; i++)
		zFree(&fobj->qs_obj.factors[i]);
	free(fobj->qs_obj.factors);

	// then free other stuff in qs
	mpz_clear(fobj->qs_obj.gmp_n);

	// free any factors found in div
	for (i=0; i<fobj->div_obj.num_factors; i++)
		zFree(&fobj->div_obj.factors[i]);
	free(fobj->div_obj.factors);

	// then free other stuff in div
	mpz_clear(fobj->div_obj.gmp_n);
	mpz_clear(fobj->div_obj.gmp_f);

	// free any factors found in nfs
	for (i=0; i<fobj->nfs_obj.num_factors; i++)
		zFree(&fobj->nfs_obj.factors[i]);
	free(fobj->nfs_obj.factors);

	// then free other stuff in nfs
	mpz_clear(fobj->nfs_obj.gmp_n);
	mpz_clear(fobj->nfs_obj.snfs_cofactor);

	//free general fobj stuff
	mpz_clear(fobj->N);
	sFree(&fobj->str_N);

	clear_factor_list(fobj);
	free(fobj->fobj_factors);

	return;
}

void alloc_factobj(fact_obj_t *fobj)
{
	int i;
	
	mpz_init(fobj->N);
	sInit(&fobj->str_N);	

	fobj->rho_obj.num_poly = 3;
	fobj->rho_obj.polynomials = (uint32 *)malloc(fobj->rho_obj.num_poly * sizeof(uint32));
	fobj->rho_obj.polynomials[0] = 3;
	fobj->rho_obj.polynomials[1] = 2;
	fobj->rho_obj.polynomials[2] = 1;
	mpz_init(fobj->rho_obj.gmp_n);
	mpz_init(fobj->rho_obj.gmp_f);

	mpz_init(fobj->pm1_obj.gmp_n);
	mpz_init(fobj->pm1_obj.gmp_f);
	
	mpz_init(fobj->pp1_obj.gmp_n);
	mpz_init(fobj->pp1_obj.gmp_f);

	fobj->ecm_obj.factors = (z *)malloc(sizeof(z));
	mpz_init(fobj->ecm_obj.gmp_n);
	mpz_init(fobj->ecm_obj.gmp_f);

	fobj->squfof_obj.factors = (z *)malloc(sizeof(z));
	mpz_init(fobj->squfof_obj.gmp_n);
	mpz_init(fobj->squfof_obj.gmp_f);

	fobj->qs_obj.factors = (z *)malloc(sizeof(z));
	mpz_init(fobj->qs_obj.gmp_n);

	fobj->div_obj.factors = (z *)malloc(sizeof(z));
	mpz_init(fobj->div_obj.gmp_n);
	mpz_init(fobj->div_obj.gmp_f);

	fobj->nfs_obj.factors = (z *)malloc(sizeof(z));
	mpz_init(fobj->nfs_obj.gmp_n);
	mpz_init(fobj->nfs_obj.snfs_cofactor);

	fobj->allocated_factors = 256;
	fobj->fobj_factors = (factor_t *)malloc(256 * sizeof(factor_t));
	for (i=0; i < fobj->allocated_factors; i++)
	{
		fobj->fobj_factors[i].type = UNKNOWN;
		fobj->fobj_factors[i].count = 0;
	}

	fobj->ecm_obj.num_factors = 0;	
	fobj->qs_obj.num_factors = 0;	
	fobj->div_obj.num_factors = 0;	
	fobj->nfs_obj.num_factors = 0;	
	fobj->num_factors = 0;

	return;
}

void reset_factobj(fact_obj_t *fobj)
{
	// keep all of the settings in fobj, but do an init/free cycle on all
	// allocated structures
	free_factobj(fobj);
	alloc_factobj(fobj);

	return;
}

void add_to_factor_list(fact_obj_t *fobj, mpz_t n)
{
	//stick the number n into the global factor list
	uint32 i;
	int found = 0, v = 0;

	if (fobj->num_factors >= fobj->allocated_factors)
	{
		fobj->allocated_factors *= 2;
		fobj->fobj_factors = (factor_t *)realloc(fobj->fobj_factors,
			fobj->allocated_factors * sizeof(factor_t));
	}

	//look to see if this factor is already in the list
	for (i=0;i<fobj->num_factors && !found; i++)
	{
		if (mpz_cmp(n, fobj->fobj_factors[i].factor) == 0)
		{
			found = 1;
			fobj->fobj_factors[i].count++;
			return;
		}
	}

	//else, put it in the list
	mpz_init(fobj->fobj_factors[fobj->num_factors].factor);
	mpz_set(fobj->fobj_factors[fobj->num_factors].factor, n);
	fobj->fobj_factors[fobj->num_factors].count = 1;
	if (gmp_base10(n) <= fobj->aprcl_prove_cutoff) /* prove primality of numbers <= aprcl_prove_cutoff digits */
	{
		int ret = 0;

		if (VFLAG > 0)
			v = (gmp_base10(n) < fobj->aprcl_display_cutoff) ? APRTCLE_VERBOSE0 : APRTCLE_VERBOSE1;
		else
			v = APRTCLE_VERBOSE0;

		if (v == APRTCLE_VERBOSE1)
			printf("\n");
		ret = mpz_aprtcle(n, v);
		if (v == APRTCLE_VERBOSE1)
			printf("\n");

		if (ret == APRTCLE_PRIME)
			fobj->fobj_factors[fobj->num_factors].type = PRIME;
		else
		{
			if (mpz_bpsw_prp(n) != PRP_COMPOSITE)
			{
				// if BPSW doesn't think this composite number is actually composite, then ask the user
				// to report this fact to the YAFU sub-forum at mersenneforum.org
				printf(" *** ATTENTION: BPSW issue found.  Please report the following number to:\n");
				printf(" *** ATTENTION: http://www.mersenneforum.org/forumdisplay.php?f=96\n");
				gmp_printf(" *** ATTENTION: n = %Zd\n", n);
			}
			fobj->fobj_factors[fobj->num_factors].type = COMPOSITE;
		}
	}
	else if (is_mpz_prp(n))
	{
		if (mpz_cmp_ui(n, 100000000) < 0)
			fobj->fobj_factors[fobj->num_factors].type = PRIME;
		else
			fobj->fobj_factors[fobj->num_factors].type = PRP;
	}
	else
		fobj->fobj_factors[fobj->num_factors].type = COMPOSITE;

	fobj->num_factors++;

	return;
}

void delete_from_factor_list(fact_obj_t *fobj, mpz_t n)
{
	//remove the number n from the global factor list
	uint32 i;

	//find the factor
	for (i=0;i<fobj->num_factors; i++)
	{
		if (mpz_cmp(n,fobj->fobj_factors[i].factor) == 0)
		{
			int j;
			// copy everything above this in the list back one position
			for (j=i; j<fobj->num_factors-1; j++)
			{
				mpz_set(fobj->fobj_factors[j].factor, fobj->fobj_factors[j+1].factor);
				fobj->fobj_factors[j].count = fobj->fobj_factors[j+1].count;
			}
			// remove the last one in the list
			fobj->fobj_factors[j].count = 0;
			mpz_clear(fobj->fobj_factors[j].factor);

			fobj->num_factors--;
			break;
		}
	}

	return;
}

void clear_factor_list(fact_obj_t *fobj)
{
	uint32 i;

	//clear this info
	for (i=0; i<fobj->num_factors; i++)
	{
		fobj->fobj_factors[i].count = 0;
		mpz_clear(fobj->fobj_factors[i].factor);
	}
	fobj->num_factors = 0;

	return;
}

int resume_check_input_match(mpz_t file_n, mpz_t input_n, mpz_t common_fact)
{
	// see if the value we
	// read from the file is really the same as the input number.
	// if the gcd of the two is a divisor of the file, then
	// this is a resume.  if the gcd is non-trival,
	// set the input to the value found in the file.

	mpz_t r;
	int retval;

	mpz_init(r);
	if (VFLAG > 1)
	{
		gmp_printf("input from file = %Zd\n", file_n);
		gmp_printf("input to yafu = %Zd\n", input_n);
	}

	//mpz_gcd(common_fact, file_n, input_n);	
	if (mpz_cmp_ui(file_n, 0) == 0)
		return 0;

	mpz_tdiv_qr(common_fact, r, input_n, file_n);	

	if (mpz_cmp_ui(r, 0) == 0)
	{
		retval = 1;
		if (VFLAG > 1)
		{
			gmp_printf("input matches with multiple of %Zd\n", common_fact);
		}
	}
	else 
	{
		retval = 0;
		mpz_set_ui(common_fact, 1);
	}
		
	mpz_clear(r);
	return retval;
}

void print_factors(fact_obj_t *fobj)
{
	uint32 i;
	int j, v;
	mpz_t tmp, tmp2;

	//always print factors unless complete silence is requested
	if (VFLAG >= 0)
	{
		printf("\n\n***factors found***\n\n");
		mpz_init(tmp);
		mpz_set_ui(tmp, 1);

		for (i=0; i<fobj->num_factors; i++)
		{

			//if (fobj->fobj_factors[i].type == COMPOSITE)
			//{
			//	for (j=0;j<fobj->fobj_factors[i].count;j++)
			//	{
			//		mpz_mul(tmp, tmp, fobj->fobj_factors[i].factor);
			//		gmp_printf("C%d = %Zd\n", gmp_base10(fobj->fobj_factors[i].factor),
			//			fobj->fobj_factors[i].factor);
			//	}
			//}
			//else if (fobj->fobj_factors[i].type == PRP)
			//{
			//	for (j=0;j<fobj->fobj_factors[i].count;j++)
			//	{
			//		mpz_mul(tmp, tmp, fobj->fobj_factors[i].factor);
			//		gmp_printf("PRP%d = %Zd\n", gmp_base10(fobj->fobj_factors[i].factor),
			//			fobj->fobj_factors[i].factor);
			//	}
			//}
			if (fobj->fobj_factors[i].type == PRIME)
			{
				// don't redo APR-CL calculations already performed by add_to_factor_list
				for (j=0;j<fobj->fobj_factors[i].count;j++)
				{
					mpz_mul(tmp, tmp, fobj->fobj_factors[i].factor);
					gmp_printf("P%d = %Zd\n", gmp_base10(fobj->fobj_factors[i].factor),
						fobj->fobj_factors[i].factor);
				}
			}
			else
			{
				//type not set, determine it now
				/* prove primality of numbers <= aprcl_prove_cutoff digits */
				if (gmp_base10(fobj->fobj_factors[i].factor) <= fobj->aprcl_prove_cutoff)
				{
					int ret = 0;
					if (VFLAG > 0)
						v = (gmp_base10(fobj->fobj_factors[i].factor) < fobj->aprcl_display_cutoff) ? APRTCLE_VERBOSE0 : APRTCLE_VERBOSE1;
					else
						v = APRTCLE_VERBOSE0;

					if (v == APRTCLE_VERBOSE1)
						printf("\n");
					ret = mpz_aprtcle(fobj->fobj_factors[i].factor, v);
					if (v == APRTCLE_VERBOSE1)
						printf("\n");

					if (ret == APRTCLE_PRIME)
					{
						for (j=0;j<fobj->fobj_factors[i].count;j++)
						{
							mpz_mul(tmp, tmp, fobj->fobj_factors[i].factor);
							gmp_printf("P%d = %Zd\n", gmp_base10(fobj->fobj_factors[i].factor),
								fobj->fobj_factors[i].factor);
						}
					}
					else
					{
						if (mpz_bpsw_prp(fobj->fobj_factors[i].factor) != PRP_COMPOSITE)
						{
							// if BPSW doesn't think this composite number is actually composite, then ask the user
							// to report this fact to the YAFU sub-forum at mersenneforum.org
							printf(" *** ATTENTION: BPSW issue found.  Please report the following number to:\n");
							printf(" *** ATTENTION: http://www.mersenneforum.org/forumdisplay.php?f=96\n");
							gmp_printf(" *** ATTENTION: n = %Zd\n", fobj->fobj_factors[i].factor);
						}
						for (j=0;j<fobj->fobj_factors[i].count;j++)
						{
							mpz_mul(tmp, tmp, fobj->fobj_factors[i].factor);
							gmp_printf("C%d = %Zd\n", gmp_base10(fobj->fobj_factors[i].factor),
								fobj->fobj_factors[i].factor);
						}
					}
				}
				else if (is_mpz_prp(fobj->fobj_factors[i].factor))
				{
					for (j=0;j<fobj->fobj_factors[i].count;j++)
					{
						mpz_mul(tmp, tmp, fobj->fobj_factors[i].factor);
						if (mpz_cmp_ui(fobj->fobj_factors[i].factor, 100000000) < 0)
							gmp_printf("P%d = %Zd\n", gmp_base10(fobj->fobj_factors[i].factor),
								fobj->fobj_factors[i].factor);
						else
							gmp_printf("PRP%d = %Zd\n", gmp_base10(fobj->fobj_factors[i].factor),
								fobj->fobj_factors[i].factor);
					}
				}
				else
				{
					for (j=0;j<fobj->fobj_factors[i].count;j++)
					{
						mpz_mul(tmp, tmp, fobj->fobj_factors[i].factor);
						gmp_printf("C%d = %Zd\n", gmp_base10(fobj->fobj_factors[i].factor),
							fobj->fobj_factors[i].factor);
					}
				}
			}
		}

		mpz_init(tmp2);
		mpz_set(tmp2, fobj->N);
		if (mpz_cmp_ui(tmp2, 1) > 0)
		{
			// non-trivial N remaining... compute and display the known co-factor
			mpz_tdiv_q(tmp2, tmp2, tmp);
			if (mpz_cmp_ui(tmp2, 1) > 0)
			{
				/* prove primality of numbers <= aprcl_prove_cutoff digits */
				if (gmp_base10(tmp2) <= fobj->aprcl_prove_cutoff)
				{
					int ret = 0;
					if (VFLAG > 0)
						v = (gmp_base10(tmp2) < fobj->aprcl_display_cutoff) ? APRTCLE_VERBOSE0 : APRTCLE_VERBOSE1;
					else
						v = APRTCLE_VERBOSE0;

					if (v == APRTCLE_VERBOSE1)
						printf("\n");
					ret = mpz_aprtcle(tmp2, v);
					if (v == APRTCLE_VERBOSE1)
						printf("\n");

					if (ret == APRTCLE_PRIME)
						gmp_printf("\n***co-factor***\nP%d = %Zd\n",
							gmp_base10(tmp2), tmp2);
					else
					{
						if (mpz_bpsw_prp(tmp2) != PRP_COMPOSITE)
						{
							// if BPSW doesn't think this composite number is actually composite, then ask the user
							// to report this fact to the YAFU sub-forum at mersenneforum.org
							printf(" *** ATTENTION: BPSW issue found.  Please report the following number to:\n");
							printf(" *** ATTENTION: http://www.mersenneforum.org/forumdisplay.php?f=96\n");
							gmp_printf(" *** ATTENTION: n = %Zd\n", tmp2);
						}
						gmp_printf("\n***co-factor***\nC%d = %Zd\n",
							gmp_base10(tmp2), tmp2);
					}
				}
				else if (is_mpz_prp(tmp2))
				{
					gmp_printf("\n***co-factor***\nPRP%d = %Zd\n",
						gmp_base10(tmp2), tmp2);
				}
				else
				{
					gmp_printf("\n***co-factor***\nC%d = %Zd\n",
						gmp_base10(tmp2), tmp2);
				}
			}			
		}
		mpz_clear(tmp);
		mpz_clear(tmp2);
	}

	return;
}

double get_qs_time_estimate(fact_obj_t *fobj, mpz_t b)
{
	//using rough empirical scaling equations, number size, information
	//on cpu type, architecture, speed, and compilation options, 
	//compute how long we think siqs would take to finish a factorization
	enum cpu_type cpu;
	double estimate;
	double freq = MEAS_CPU_FREQUENCY;
	int digits = gmp_base10(b);

	cpu = yafu_get_cpu_type();

	estimate = fobj->qs_obj.qs_multiplier * exp(fobj->qs_obj.qs_exponent * digits);
	estimate = estimate * fobj->qs_obj.qs_tune_freq / freq; 	

	//adjust for multi-threaded qs
	//if we assume threading is perfect, we'll get a smaller estimate for
	//qs than we can really achieve, resulting in less ECM, so fudge it a bit
	if (THREADS > 1)
	{
		switch (cpu)
		{
		case 0:
		case 1:
		case 2:	
		case 3:
		case 4:
		case 5:			
		case 6:
		case 7:
		case 8:
			estimate = estimate / ((double)THREADS * 0.75);
			break;
		case 9:
		case 10:
			estimate = estimate / ((double)THREADS * 0.90);
			break;

		default:
			estimate = estimate / ((double)THREADS * 0.75);
			break;
		}
	}

	if (VFLAG >= 2)
		printf("fac: QS time estimation from tune data = %1.2f sec\n", estimate);

	return estimate;
}

double get_gnfs_time_estimate(fact_obj_t *fobj, mpz_t b)
{
	//using rough empirical scaling equations, number size, information
	//on cpu type, architecture, speed, and compilation options, 
	//compute how long we think gnfs would take to finish a factorization
	enum cpu_type cpu;
	double estimate;
	double freq = MEAS_CPU_FREQUENCY;
	int digits = gmp_base10(b);

	cpu = yafu_get_cpu_type();
	
	estimate = fobj->nfs_obj.gnfs_multiplier * exp(fobj->nfs_obj.gnfs_exponent * digits);
	estimate = estimate * fobj->nfs_obj.gnfs_tune_freq / freq; 

	//adjust for multi-threaded nfs
	//if we assume threading is perfect, we'll get a smaller estimate for
	//nfs than we can really achieve, resulting in less ECM, so fudge it a bit
	if (THREADS > 1)
	{
		switch (cpu)
		{
		case 0:
		case 1:
		case 2:	
		case 3:
		case 4:
		case 5:			
		case 6:
		case 7:
		case 8:
			estimate = estimate / ((double)THREADS * 0.75);
			break;
		case 9:
		case 10:
			estimate = estimate / ((double)THREADS * 0.90);
			break;

		default:
			estimate = estimate / ((double)THREADS * 0.75);
			break;
		}
	}

	if (VFLAG >= 2)
		printf("fac: GNFS time estimation from tune data = %1.2f sec\n", estimate);

	return estimate;
}

void do_work(enum factorization_state method, factor_work_t *fwork, 
	mpz_t b, fact_obj_t *fobj)
{
	uint32 tmp1;
	uint64 tmp2;	
	struct timeval tstart, tstop;
	double t_time;
	uint32 curves_done;

	gettimeofday(&tstart, NULL);

	switch (method)
	{
	case state_trialdiv:

        // if larger than a small bound, do a perfect power check
        fobj->prime_threshold = fwork->tdiv_max_limit * fwork->tdiv_max_limit;
        
        if ((mpz_cmp_ui(b, fobj->prime_threshold) > 1) && 
            mpz_perfect_power_p(b))
        {
            FILE *flog;

            if (VFLAG > 0)
            {
                printf("fac: input is a perfect power\n");
            }

            flog = fopen(fobj->flogname, "a");
            if (flog != NULL)
            {
                logprint(flog, "input is a perfect power\n");
                fclose(flog);
            }

            factor_perfect_power(fobj, b);

            mpz_set(fobj->N, b);
            break;
        }

        // then do all of the tdiv work requested
        if (VFLAG >= 0)
            printf("div: primes less than %d\n", fwork->tdiv_max_limit);
        
        mpz_set(fobj->div_obj.gmp_n, b);
        fobj->div_obj.print = 1;
        fobj->div_obj.limit = fwork->tdiv_max_limit;
        zTrial(fobj);
        mpz_set(b, fobj->div_obj.gmp_n);

        // record the work done
        fwork->tdiv_limit = fwork->tdiv_max_limit;

        // measure time for this completed work
        gettimeofday(&tstop, NULL);
        t_time = yafu_difftime(&tstart, &tstop);

        fwork->trialdiv_time = t_time;
        fwork->total_time += t_time;

		break;

	case state_rho:
		// do all of the rho work requested
		mpz_set(fobj->rho_obj.gmp_n,b);
		brent_loop(fobj);
		mpz_set(b,fobj->rho_obj.gmp_n);

		// record the work done
		fwork->rho_bases = fwork->rho_max_bases;
		fwork->rho_iterations = fwork->rho_max_iterations;

		// measure time for this completed work
		gettimeofday (&tstop, NULL);
        t_time = yafu_difftime(&tstart, &tstop);

		fwork->rho_time = t_time;
		fwork->total_time += t_time;
		break;

	case state_fermat:
		// do all of the fermat work requested
		if (VFLAG >= 0)
			printf("fmt: %d iterations\n", fwork->fermat_max_iterations);
		mpz_set(fobj->div_obj.gmp_n,b);
		zFermat(fwork->fermat_max_iterations, 1, fobj);
		mpz_set(b,fobj->div_obj.gmp_n);

		// record the work done
		fwork->fermat_iterations = fwork->fermat_max_iterations;

		// measure time for this completed work
		gettimeofday (&tstop, NULL);
        t_time = yafu_difftime(&tstart, &tstop);

		fwork->fermat_time = t_time;
		fwork->total_time += t_time;

		break;

	case state_ecm_15digit:
	case state_ecm_20digit:
	case state_ecm_25digit:
	case state_ecm_30digit:
	case state_ecm_35digit:
	case state_ecm_40digit:
	case state_ecm_45digit:
	case state_ecm_50digit:
	case state_ecm_55digit:
	case state_ecm_60digit:
	case state_ecm_65digit:
		tmp1 = fobj->ecm_obj.B1;
		tmp2 = fobj->ecm_obj.B2;
		fobj->ecm_obj.B1 = fwork->B1;
		fobj->ecm_obj.B2 = fwork->B2;
		fobj->ecm_obj.num_curves = fwork->curves;
		mpz_set(fobj->ecm_obj.gmp_n, b);
		curves_done = ecm_loop(fobj);
		mpz_set(b, fobj->ecm_obj.gmp_n);
		fobj->ecm_obj.B1 = tmp1;
		fobj->ecm_obj.B2 = tmp2;

		// record the work done
		set_ecm_curves_done(fwork, method, 
			get_ecm_curves_done(fwork, method) + curves_done);
		
		// measure time for this completed work
		gettimeofday (&tstop, NULL);
        t_time = yafu_difftime(&tstart, &tstop);

		fwork->ecm_time += t_time;
		fwork->total_time += t_time;
		break;

	case state_pp1_lvl1:
	case state_pp1_lvl2:
	case state_pp1_lvl3:
		tmp1 = fobj->pp1_obj.B1;
		tmp2 = fobj->pp1_obj.B2;
		fobj->pp1_obj.B1 = fwork->B1;
		fobj->pp1_obj.B2 = fwork->B2;
		mpz_set(fobj->pp1_obj.gmp_n,b);
		fobj->pp1_obj.numbases = fwork->curves;
		williams_loop(fobj);
		mpz_set(b,fobj->pp1_obj.gmp_n);
		fobj->pp1_obj.B1 = tmp1;
		fobj->pp1_obj.B2 = tmp2;

		// record the work done
		if (method == state_pp1_lvl1)
			fwork->pp1_lvl1_curves = fwork->curves;
		else if (method == state_pp1_lvl2)
			fwork->pp1_lvl2_curves = fwork->curves;
		else if (method == state_pp1_lvl3)
			fwork->pp1_lvl3_curves = fwork->curves;

		// measure time for this completed work
		gettimeofday (&tstop, NULL);
        t_time = yafu_difftime(&tstart, &tstop);

		fwork->pp1_time += t_time;
		fwork->total_time += t_time;
		break;

	case state_pm1_lvl1:
	case state_pm1_lvl2:
	case state_pm1_lvl3:
		tmp1 = fobj->pm1_obj.B1;
		tmp2 = fobj->pm1_obj.B2;
		fobj->pm1_obj.B1 = fwork->B1;
		fobj->pm1_obj.B2 = fwork->B2;
		mpz_set(fobj->pm1_obj.gmp_n,b);
		pollard_loop(fobj);
		mpz_set(b,fobj->pm1_obj.gmp_n);
		fobj->pm1_obj.B1 = tmp1;
		fobj->pm1_obj.B2 = tmp2;

		// record the work done
		if (method == state_pm1_lvl1)
			fwork->pm1_lvl1_curves = fwork->curves;
		else if (method == state_pm1_lvl2)
			fwork->pm1_lvl2_curves = fwork->curves;
		else if (method == state_pm1_lvl3)
			fwork->pm1_lvl3_curves = fwork->curves;
		
		// measure time for this completed work
		gettimeofday (&tstop, NULL);
        t_time = yafu_difftime(&tstart, &tstop);

		fwork->pm1_time += t_time;
		fwork->total_time += t_time;
		break;

	case state_qs:
		mpz_set(fobj->qs_obj.gmp_n,b);
		SIQS(fobj);
		mpz_set(b,fobj->qs_obj.gmp_n);

		// measure time for this completed work
		gettimeofday (&tstop, NULL);
        t_time = yafu_difftime(&tstart, &tstop);

		fwork->qs_time = t_time;
		if (VFLAG > 0)
			printf("pretesting / qs ratio was %1.2f\n", 
				fwork->total_time / t_time); 
		break;

	case state_nfs:
		mpz_set(fobj->nfs_obj.gmp_n,b);
		nfs(fobj);
		mpz_set(b,fobj->nfs_obj.gmp_n);

		// measure time for this completed work
		gettimeofday (&tstop, NULL);
        t_time = yafu_difftime(&tstart, &tstop);

		fwork->nfs_time = t_time;
		if (VFLAG > 0)
			printf("pretesting / nfs ratio was %1.2f\n", 
				fwork->total_time / t_time); 
		break;

	default:
		printf("nothing to do for method %d\n", method);
		break;
	}

	return;
}

int check_if_done(fact_obj_t *fobj, mpz_t N)
{
	int i, done = 0;
	mpz_t tmp;

	mpz_init(tmp);
	mpz_set_ui(tmp, 1);

	/* if the user only wants to find one factor, check for that here... */
	if (fobj->autofact_obj.want_only_1_factor && (fobj->num_factors >= 1))
	{
		done = 1;
		mpz_clear(tmp);
		return done;
	}

	// check if the number is completely factorized
	for (i=0; i<fobj->num_factors; i++)
	{		
		int j;
		for (j=0; j<fobj->fobj_factors[i].count; j++)
			mpz_mul(tmp, tmp, fobj->fobj_factors[i].factor);
	}

	if (mpz_cmp(N,tmp) == 0)
	{		
		// yes, they are equal.  make sure everything is prp or prime.
		done = 0;
		while (!done)
		{
			done = 1;
			for (i=0; i<fobj->num_factors; i++)
			{
				if (!is_mpz_prp(fobj->fobj_factors[i].factor))
				{
					fobj->refactor_depth++;
					if (fobj->refactor_depth > 3)
					{
						printf("too many refactorization attempts, aborting\n");
						done = 1;
						break;
					}
					else
					{
						//printf("ignoring refactorization of composite factor\n");
						fact_obj_t *fobj_refactor;
						int j;

						if (VFLAG > 0)
							printf("\nComposite result found, starting re-factorization\n");

						// load the new fobj with this number
						fobj_refactor = (fact_obj_t *)malloc(sizeof(fact_obj_t));
						init_factobj(fobj_refactor);
						mpz_set(fobj_refactor->N, fobj->fobj_factors[i].factor);

						// recurse on factor
						factor(fobj_refactor);

						// remove the factor from the original list
						delete_from_factor_list(fobj, fobj->fobj_factors[i].factor);

						// add all factors found during the refactorization
						for (j=0; j< fobj_refactor->num_factors; j++)
						{
							int k;
							for (k=0; k < fobj_refactor->fobj_factors[j].count; k++)
								add_to_factor_list(fobj, fobj_refactor->fobj_factors[j].factor);
						}

						// free temps
						free_factobj(fobj_refactor);
						free(fobj_refactor);

						// check again, since this factorization could have added new
						// composite factors
						done = 0;
					}
				}
			}
		}
	}

	mpz_clear(tmp);
	return done;
}

enum factorization_state get_next_state(factor_work_t *fwork, fact_obj_t *fobj)
{
	enum factorization_state next_state;

	// check each state's completed work against the maximum.
	// return the first one not complete.
	if (fwork->tdiv_limit < fwork->tdiv_max_limit)
		return state_trialdiv;	// always do this state if not done
	else if (fwork->fermat_iterations < fwork->fermat_max_iterations)
		return state_fermat;	// always do this state if not done
	else if (fwork->rho_bases < fwork->rho_max_bases)
		return state_rho;		// always do this state if not done
	else if (fwork->pp1_lvl1_curves < fwork->pp1_max_lvl1_curves)
		next_state = state_pp1_lvl1;
	else if (fwork->pm1_lvl1_curves < fwork->pm1_max_lvl1_curves)
		next_state = state_pm1_lvl1;
	else if (fwork->ecm_15digit_curves < fwork->ecm_max_15digit_curves)
		next_state = state_ecm_15digit;
	else if (fwork->ecm_20digit_curves < fwork->ecm_max_20digit_curves)
		next_state = state_ecm_20digit;
	else if (fwork->ecm_25digit_curves < fwork->ecm_max_25digit_curves)
		next_state = state_ecm_25digit;
	else if (fwork->pp1_lvl2_curves < fwork->pp1_max_lvl2_curves)
		next_state = state_pp1_lvl2;
	else if (fwork->pm1_lvl2_curves < fwork->pm1_max_lvl2_curves)
		next_state = state_pm1_lvl2;
	else if (fwork->ecm_30digit_curves < fwork->ecm_max_30digit_curves)
		next_state = state_ecm_30digit;
	else if (fwork->pp1_lvl3_curves < fwork->pp1_max_lvl3_curves)
		next_state = state_pp1_lvl3;
	else if (fwork->pm1_lvl3_curves < fwork->pm1_max_lvl3_curves)
		next_state = state_pm1_lvl3;
	else if (fwork->ecm_35digit_curves < fwork->ecm_max_35digit_curves)
		next_state = state_ecm_35digit;
	else if (fwork->ecm_40digit_curves < fwork->ecm_max_40digit_curves)
		next_state = state_ecm_40digit;
	else if (fwork->ecm_45digit_curves < fwork->ecm_max_45digit_curves)
		next_state = state_ecm_45digit;
	else if (fwork->ecm_50digit_curves < fwork->ecm_max_50digit_curves)
		next_state = state_ecm_50digit;
    else if (fwork->ecm_55digit_curves < fwork->ecm_max_55digit_curves)
        next_state = state_ecm_55digit;
    else if (fwork->ecm_60digit_curves < fwork->ecm_max_60digit_curves)
        next_state = state_ecm_60digit;
    else if (fwork->ecm_65digit_curves < fwork->ecm_max_65digit_curves)
        next_state = state_ecm_65digit;
	else
		next_state = state_nfs;

	// modify according to user preferences if necessary
	switch (next_state)
	{
		case state_ecm_15digit:
		case state_ecm_20digit:
		case state_ecm_25digit:
		case state_ecm_30digit:
		case state_ecm_35digit:
		case state_ecm_40digit:
		case state_ecm_45digit:
		case state_ecm_50digit:
		case state_ecm_55digit:
		case state_ecm_60digit:
		case state_ecm_65digit:
			if (fobj->autofact_obj.yafu_pretest_plan == PRETEST_NOECM)
				next_state = state_nfs;
			break;
			
		default:

			break;
	}

	if (fobj->autofact_obj.yafu_pretest_plan == PRETEST_NONE)
		next_state = state_nfs;

	return next_state;
}

int check_tune_params(fact_obj_t *fobj)
{
	if (fobj->qs_obj.qs_multiplier == 0 || 
		fobj->qs_obj.qs_exponent == 0 || 
		fobj->qs_obj.qs_tune_freq == 0 ||
		fobj->nfs_obj.gnfs_multiplier == 0 || 
		fobj->nfs_obj.gnfs_exponent == 0 || 
		fobj->nfs_obj.gnfs_tune_freq == 0)
	{
        if (VFLAG >= 0)
        {
            printf("check tune params contained invalid parameter(s), ignoring tune info.\n");
        }

        if (VFLAG > 0)
        {
            printf("\tqs_mult = %e\n", fobj->qs_obj.qs_multiplier);
            printf("\tqs_exp = %e\n", fobj->qs_obj.qs_exponent);
            printf("\tqs_freq = %e\n", fobj->qs_obj.qs_tune_freq);
            printf("\tnfs_mult = %e\n", fobj->nfs_obj.gnfs_multiplier);
            printf("\tnfs_exp = %e\n", fobj->nfs_obj.gnfs_exponent);
            printf("\tnfs_freq = %e\n", fobj->nfs_obj.gnfs_tune_freq);
        }
		return 0;
	}

	return 1;
}

void set_work_params(factor_work_t *fwork, enum factorization_state state)
{
	switch (state)
	{
	case state_pp1_lvl1:
		fwork->B1 = 20000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = 3;
		break;

	case state_pp1_lvl2:
		fwork->B1 = 1250000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = 3;
		break;

	case state_pp1_lvl3:
		fwork->B1 = 5000000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = 3;
		break;

	case state_pm1_lvl1:
		fwork->B1 = 150000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = 1;
		break;

	case state_pm1_lvl2:
		fwork->B1 = 3750000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = 1;
		break;

	case state_pm1_lvl3:
		fwork->B1 = 15000000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = 1;
		break;

	case state_ecm_15digit:
		fwork->B1 = 2000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = fwork->ecm_max_15digit_curves;
		break;

	case state_ecm_20digit:
		fwork->B1 = 11000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = fwork->ecm_max_20digit_curves;
		break;

	case state_ecm_25digit:
		fwork->B1 = 50000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = fwork->ecm_max_25digit_curves;
		break;

	case state_ecm_30digit:
		fwork->B1 = 250000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = fwork->ecm_max_30digit_curves;
		break;

	case state_ecm_35digit:
		fwork->B1 = 1000000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = fwork->ecm_max_35digit_curves;
		break;

	case state_ecm_40digit:
		fwork->B1 = 3000000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = fwork->ecm_max_40digit_curves;
		break;

	case state_ecm_45digit:
		fwork->B1 = 11000000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = fwork->ecm_max_45digit_curves;
		break;

	case state_ecm_50digit:
		fwork->B1 = 43000000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = fwork->ecm_max_50digit_curves;
		break;

	case state_ecm_55digit:
		fwork->B1 = 110000000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = fwork->ecm_max_55digit_curves;
		break;

	case state_ecm_60digit:
		fwork->B1 = 260000000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = fwork->ecm_max_60digit_curves;
		break;

	case state_ecm_65digit:
		fwork->B1 = 850000000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = fwork->ecm_max_65digit_curves;
		break;

	default:
		fwork->B1 = 0;	//error condition
		fwork->B2 = 0;
		fwork->curves = 0;
		break;

	}

	return;
}

uint32 get_ecm_curves_done(factor_work_t *fwork, enum factorization_state state)
{
	uint32 curves_done;

	switch (state)
	{
	case state_ecm_15digit:
		curves_done = fwork->ecm_15digit_curves;
		break;
	case state_ecm_20digit:
		curves_done = fwork->ecm_20digit_curves;
		break;
	case state_ecm_25digit:
		curves_done = fwork->ecm_25digit_curves;
		break;
	case state_ecm_30digit:
		curves_done = fwork->ecm_30digit_curves;
		break;
	case state_ecm_35digit:
		curves_done = fwork->ecm_35digit_curves;
		break;
	case state_ecm_40digit:
		curves_done = fwork->ecm_40digit_curves;
		break;
	case state_ecm_45digit:
		curves_done = fwork->ecm_45digit_curves;
		break;
	case state_ecm_50digit:
		curves_done = fwork->ecm_50digit_curves;
		break;
	case state_ecm_55digit:
		curves_done = fwork->ecm_55digit_curves;
		break;
	case state_ecm_60digit:
		curves_done = fwork->ecm_60digit_curves;
		break;
	case state_ecm_65digit:
		curves_done = fwork->ecm_65digit_curves;
		break;
	default:
		curves_done = 0;
		break;
	}

	return curves_done;
}

uint32 set_ecm_curves_done(factor_work_t *fwork, enum factorization_state state, uint32 curves_done)
{
	switch (state)
	{
	case state_ecm_15digit:
		fwork->ecm_15digit_curves = curves_done;
		break;
	case state_ecm_20digit:
		fwork->ecm_20digit_curves = curves_done;
		break;
	case state_ecm_25digit:
		fwork->ecm_25digit_curves = curves_done;
		break;
	case state_ecm_30digit:
		fwork->ecm_30digit_curves = curves_done;
		break;
	case state_ecm_35digit:
		fwork->ecm_35digit_curves = curves_done;
		break;
	case state_ecm_40digit:
		fwork->ecm_40digit_curves = curves_done;
		break;
	case state_ecm_45digit:
		fwork->ecm_45digit_curves = curves_done;
		break;
	case state_ecm_50digit:
		fwork->ecm_50digit_curves = curves_done;
		break;
	case state_ecm_55digit:
		fwork->ecm_55digit_curves = curves_done;
		break;
	case state_ecm_60digit:
		fwork->ecm_60digit_curves = curves_done;
		break;
	case state_ecm_65digit:
		fwork->ecm_65digit_curves = curves_done;
		break;
	default:
		printf("don't know how to set curves for state %d\n", state);
		exit(1);
		break;
	}

	return curves_done;
}

uint32 get_max_ecm_curves(factor_work_t *fwork, enum factorization_state state)
{
	uint32 max_curves;

	switch (state)
	{
	case state_ecm_15digit:
		max_curves = fwork->ecm_max_15digit_curves;
		break;
	case state_ecm_20digit:
		max_curves = fwork->ecm_max_20digit_curves;
		break;
	case state_ecm_25digit:
		max_curves = fwork->ecm_max_25digit_curves;
		break;
	case state_ecm_30digit:
		max_curves = fwork->ecm_max_30digit_curves;
		break;
	case state_ecm_35digit:
		max_curves = fwork->ecm_max_35digit_curves;
		break;
	case state_ecm_40digit:
		max_curves = fwork->ecm_max_40digit_curves;
		break;
	case state_ecm_45digit:
		max_curves = fwork->ecm_max_45digit_curves;
		break;
	case state_ecm_50digit:
		max_curves = fwork->ecm_max_50digit_curves;
		break;
	case state_ecm_55digit:
		max_curves = fwork->ecm_max_55digit_curves;
		break;
	case state_ecm_60digit:
		max_curves = fwork->ecm_max_60digit_curves;
		break;
	case state_ecm_65digit:
		max_curves = fwork->ecm_max_65digit_curves;
		break;
	default:
		max_curves = 0;
		break;
	}

	return max_curves;
}

double compute_ecm_work_done(factor_work_t *fwork, int disp_levels, FILE *log)
{
	// there is probably a more elegant way to do this involving dickman's function
	// or something, but we can get a reasonable estimate using empirical data
	// for our fixed set of B1/B2 values.
	double tlevels[NUM_ECM_LEVELS];
	uint32 curves_done;
	int i, j;

    if (log != NULL)
    {
        logprint(log, "ecm work completed:\n");
    }

	// compute the %done of each tlevel
	for (i=0; i < NUM_ECM_LEVELS; i++)
	{
		enum factorization_state k;
		tlevels[i] = 0;
		
		for (k=state_ecm_15digit, j=0; k <= state_ecm_65digit; k++, j++)
		{            
			curves_done = get_ecm_curves_done(fwork, k);			
			tlevels[i] += (double)curves_done / ecm_data[i][j];
		}

        if ((VFLAG >= 1) && disp_levels && (tlevels[i] > 0.01))
        {
            printf("\tt%d: %1.2f\n", ecm_levels[i], tlevels[i]);
        }

        if ((log != NULL) && (tlevels[i] > 0.01))
        {
            logprint(log, "\tt%d: %1.2f\n", ecm_levels[i], tlevels[i]);
        }
	}

	// find the first one less than 1
    for (i = 0; i < NUM_ECM_LEVELS; i++)
    {
        if (tlevels[i] < 1)
        {
            break;
        }
    }

	// estimate the t level done by extrapolating between this and the previous one
	// assuming they are all spaced 5 digits apart.
    if (i == 0)
    {
        return 0;
    }
    else
    {
        return ecm_levels[i - 1] + 5 * tlevels[i];
    }

}

enum factorization_state schedule_work(factor_work_t *fwork, mpz_t b, fact_obj_t *fobj)
{
	int have_tune;
    int i;
	enum factorization_state next_state;
	int numdigits = gmp_base10(b);
	double target_digits;
	double work_done;
	FILE *flog;
	snfs_t *poly;

	// get the next factorization state that hasn't been completed
	next_state = get_next_state(fwork, fobj);

	// make sure a minimum amount of work is done
	if (next_state == state_trialdiv ||
		next_state == state_fermat ||
		next_state == state_rho)
	{
		return next_state;
	}

	// check to see if 'tune' has been run or not
	have_tune = check_tune_params(fobj);		

	// set target pretesting depth, depending on user selection and whether or not
	// the input is both big enough and snfsable...
    if ((numdigits >= fobj->autofact_obj.qs_snfs_xover) && (fobj->autofact_obj.has_snfs_form < 0))
    {
        mpz_set(fobj->nfs_obj.gmp_n, b);
#ifdef USE_NFS
        poly = snfs_find_form(fobj);

        if (poly != NULL)
        {
            fobj->autofact_obj.has_snfs_form = 1;
            // The actual poly is not needed now, so just get rid of it.
            snfs_clear(poly);
            free(poly);
        }
        else
        {
            fobj->autofact_obj.has_snfs_form = 0;
        }
#else
		fobj->autofact_obj.has_snfs_form = 0;
#endif
	}
	
    if (fobj->autofact_obj.only_pretest > 1)
    {
        target_digits = fobj->autofact_obj.only_pretest;
    }
    else if (fobj->autofact_obj.yafu_pretest_plan == PRETEST_DEEP)
    {
        target_digits = 1. * (double)numdigits / 3.;
    }
    else if (fobj->autofact_obj.yafu_pretest_plan == PRETEST_LIGHT)
    {
        target_digits = 2. * (double)numdigits / 9.;
    }
    else if (fobj->autofact_obj.yafu_pretest_plan == PRETEST_CUSTOM)
    {
        target_digits = (double)numdigits * fobj->autofact_obj.target_pretest_ratio;
    }
    else
    {
        target_digits = 4. * (double)numdigits / 13.;
    }

#ifdef USE_NFS
	if ((fobj->autofact_obj.has_snfs_form == 1) && (fobj->nfs_obj.gnfs == 0) && 
        (strcmp(fobj->autofact_obj.plan_str,"normal") == 0) &&
        (fobj->autofact_obj.only_pretest <= 1))
	{
		// 1) we found a snfs polynomial for the input.
		// 2) user has not specifically chosen gnfs.
		// 3) user has not already adjusted the ECM plan.
        // 4) user has not specified an only-pretest bound.
		// So it's looking like this will be an SNFS job and we should scale back the ECM effort .
		// however we also need to check if easier by gnfs... this involves a little more work.
		// This work will be duplicated if we actually get to SNFS (i.e., we don't find an
		// ECM factor), but this goes fast.
		// temporarily set verbosity high so we don't spam the screen.  
		int tmpV = VFLAG;
		snfs_t *polys = NULL;
		int npoly;
		int gnfs_size;

		VFLAG = -1;

		mpz_set(fobj->nfs_obj.gmp_n, b);
		poly = snfs_find_form(fobj);

		// with the form detected, create a good polynomial
        if (poly->form_type == SNFS_XYYXF)
        {
            polys = gen_xyyxf_poly(fobj, poly, &npoly);
        }
        else
        {
            polys = gen_brent_poly(fobj, poly, &npoly);
        }

		// then scale and rank them
		snfs_scale_difficulty(polys, npoly);
		npoly = snfs_rank_polys(fobj, polys, npoly);

		VFLAG = tmpV;

		// and test the best one, compared to gnfs or qs, depending on 
        // which one will run.
		gnfs_size = est_gnfs_size_via_poly(&polys[0]);

		if (gnfs_size <= (mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10) + 3))
		{
			// Finally - to the best of our knowledge this will be a SNFS job.
			// Since we are in factor(), we'll proceed with any ecm required, but adjust 
			// the plan ratio in accord with the snfs job.
			if (VFLAG >= 0) 
			{
				printf("fac: ecm effort reduced from %1.2f to %1.2f: input has snfs form\n",
					target_digits, target_digits / 1.2857);
			}
			target_digits /= 1.2857;
		}
        else
        {
            if (mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10) < fobj->autofact_obj.qs_gnfs_xover)
            {
                // don't consider the qs/snfs cutoff any more
                fobj->autofact_obj.has_snfs_form = 0;

                if (VFLAG >= 0)
                {
                    printf("fac: ecm effort maintained at %1.2f: input better by qs\n",
                        target_digits);
                }
            }
            else
            {
                if (VFLAG >= 0)
                {
                    printf("fac: ecm effort maintained at %1.2f: input better by gnfs\n",
                        target_digits);
                }
            }
        }

        // don't need the poly anymore
        snfs_clear(poly);
        free(poly);

        // or the list of best polys
        for (i = 0; i<npoly; i++)
        {
            snfs_clear(&polys[i]);
        }
        free(polys);
	}
#endif

	// get the current amount of work done - only print status prior to 
	// ecm steps
	switch (next_state)
	{
		case state_ecm_15digit:
		case state_ecm_20digit:
		case state_ecm_25digit:
		case state_ecm_30digit:
		case state_ecm_35digit:
		case state_ecm_40digit:
		case state_ecm_45digit:
		case state_ecm_50digit:
		case state_ecm_55digit:
		case state_ecm_60digit:
		case state_ecm_65digit:
			if (VFLAG >= 1)
				printf("fac: setting target pretesting digits to %1.2f\n", target_digits);
			
			work_done = compute_ecm_work_done(fwork, 1, NULL);
			
			if (VFLAG >= 1)
				printf("fac: estimated sum of completed work is t%1.2f\n", work_done);

			break;

		default:
			work_done = compute_ecm_work_done(fwork, 0, NULL);
			break;
	}

	// if there is a trivial amount of ecm to do, skip directly to a sieve method
	if ((target_digits < 15) && (numdigits <= 45))
	{
		if (VFLAG > 0)
			printf("fac: trivial ECM work to do... skipping to sieve method\n");
		next_state = state_nfs;
	}

	// handle the case where the next state is a sieve method
	if ((next_state == state_nfs) || (work_done > target_digits) ||
		((work_done > fobj->autofact_obj.only_pretest) && 
		(fobj->autofact_obj.only_pretest > 1)))
	{
		flog = fopen(fobj->flogname,"a");
		if (flog == NULL)
		{
			printf("fopen error: %s\n", strerror(errno));
			printf("could not open %s for writing\n",fobj->flogname);
			flog = stderr;
		}
		logprint(flog,"final ECM pretested depth: %1.2f\n", work_done);		

		// if the user specified -pretest, with or without arguments,
		// we should stop factoring now that ecm is done.  this covers the
		// case where the user specified a pretest work amount that was
		// too large as determined by factor
		if (fobj->autofact_obj.only_pretest)
		{
			fclose(flog);
			return state_done;
		}

		logprint(flog,"scheduler: switching to sieve method\n");
		if (flog != stderr)
			fclose(flog);

		if (!have_tune || fobj->autofact_obj.prefer_xover)
		{
			// use a hard cutoff - within reason
            if ((((numdigits > fobj->autofact_obj.qs_snfs_xover) && 
                (fobj->autofact_obj.has_snfs_form)) ||
                (numdigits > fobj->autofact_obj.qs_gnfs_xover)) &&
				(numdigits >= 75))
				return state_nfs;
			else
				return state_qs;			
		}
		else
		{
			double qs_time_est, gnfs_time_est;
		
			// compute the time to factor using estimates derived during 'tune'.
			qs_time_est = get_qs_time_estimate(fobj, b);
			gnfs_time_est = get_gnfs_time_estimate(fobj, b);

			if (qs_time_est < gnfs_time_est)
				return state_qs;
			else
				return state_nfs;
		}
	}

	// set the work parameters for the current state
	set_work_params(fwork, next_state);

	switch (next_state)
	{
		case state_ecm_15digit:
		case state_ecm_20digit:
		case state_ecm_25digit:
		case state_ecm_30digit:
		case state_ecm_35digit:
		case state_ecm_40digit:
		case state_ecm_45digit:
		case state_ecm_50digit:
		case state_ecm_55digit:
		case state_ecm_60digit:
		case state_ecm_65digit:
			// figure out how many curves at this level need to be done 
			// to get to the target level
			interp_and_set_curves(fwork, fobj, next_state, work_done,
				target_digits, 1);

			break;

		default:
			// non-ecm curves are set with set_work_params above
			break;
	}

	return next_state;
}

void interp_and_set_curves(factor_work_t *fwork, fact_obj_t *fobj, 
	enum factorization_state state, double work_done,
	double target_digits, int log_results)
{
	// do a binary search on the target state's amount of work.
	// probably there is a more elegant way to compute this, but this seems
	// to work.
	double work_low, work_high, work;
	uint32 tmp_curves;

	// if there is a user specified pretest value, use it, regardless if it
    // means over or under ecm'ing something.
    if (fobj->autofact_obj.only_pretest > 1)
    {
        //target_digits = MIN(target_digits, fobj->autofact_obj.only_pretest);
        target_digits = fobj->autofact_obj.only_pretest;
    }


	work_low = get_ecm_curves_done(fwork, state);
	work_high = get_max_ecm_curves(fwork, state);		
	work = (work_low + work_high) / 2;

    if ((VFLAG >= 1) && log_results)
    {
        printf("fac: work done at B1=%u: %1.0f curves, max work = %1.0f curves\n",
            fwork->B1, work_low, work_high);
    }

	tmp_curves = work_low;		
	while ((work_high - work_low) > 1)
	{
        double compute;

		set_ecm_curves_done(fwork, state, (uint32)work);       
        compute = compute_ecm_work_done(fwork, 0, NULL);

		if (compute > target_digits)
		{
			work_high = work;
			work = (work_high + work_low) / 2;							
		}
		else					
		{
			work_low = work;
			work = (work_high + work_low) / 2;							
		}
	}

	
    set_ecm_curves_done(fwork, state, tmp_curves);    
	fwork->curves = (uint32)ceil(work);

    if ((tmp_curves + fwork->curves) > get_max_ecm_curves(fwork, state))
    {
        fwork->curves = get_max_ecm_curves(fwork, state) - tmp_curves;
    }

    if ((VFLAG >= 1) && log_results)
    {
        printf("fac: %u more curves at B1=%u needed to get to t%1.2f\n",
            fwork->curves, fwork->B1, target_digits);
    }

    if (log_results)
	{
		FILE *flog;
		flog = fopen(fobj->flogname,"a");
		logprint(flog,"current ECM pretesting depth: %1.2f\n", work_done);
		logprint(flog,"scheduled %u curves at B1=%u toward target "
			"pretesting depth of %1.2f\n", fwork->curves, fwork->B1, target_digits);
		fclose(flog);
	}

	return;
}

void init_factor_work(factor_work_t *fwork, fact_obj_t *fobj)
{
	enum factorization_state interp_state = state_idle;
	// initialize max allowed work fields (note: maybe this structure should
	// be visible to the top level driver so that the user can edit values in it).
	// default values taken from gmp-ecm README, version 6.3
	fwork->ecm_max_15digit_curves = 30;		//2k
	fwork->ecm_max_20digit_curves = 74;		//11k
	fwork->ecm_max_25digit_curves = 214;	//50k
	fwork->ecm_max_30digit_curves = 430;	//250k
	fwork->ecm_max_35digit_curves = 904;	//1M
	fwork->ecm_max_40digit_curves = 2350;	//3M
	fwork->ecm_max_45digit_curves = 4480;	//11M
	fwork->ecm_max_50digit_curves = 7553;	//43M
	fwork->ecm_max_55digit_curves = 17769;	//110M
	fwork->ecm_max_60digit_curves = 42017;	//260M
	fwork->ecm_max_65digit_curves = 69408;	//850M
	fwork->tdiv_max_limit = fobj->div_obj.limit;
    fwork->fermat_iterations = 0;
	fwork->fermat_max_iterations = fobj->div_obj.fmtlimit;
	fwork->rho_max_bases = 3;
	fwork->rho_max_iterations = fobj->rho_obj.iterations;
	fwork->pp1_max_lvl1_curves = 0;
	fwork->pp1_max_lvl2_curves = 0;
	fwork->pp1_max_lvl3_curves = 0;
	fwork->pm1_max_lvl1_curves = 1;
	fwork->pm1_max_lvl2_curves = 1;
	fwork->pm1_max_lvl3_curves = 1;	
	fwork->total_time = 0;	
	fwork->trialdiv_time = 0;
	fwork->rho_time = 0;
	fwork->pp1_time = 0;
	fwork->pm1_time = 0;
	fwork->ecm_time = 0;
	fwork->qs_time = 0;
	fwork->nfs_time = 0;

	fwork->rho_bases = 0;

	fwork->pp1_lvl1_curves = 0;
	fwork->pp1_lvl2_curves = 0;
	fwork->pp1_lvl3_curves = 0;

	fwork->pm1_lvl1_curves = 0;
	fwork->pm1_lvl2_curves = 0;
	fwork->pm1_lvl3_curves = 0;

	fwork->ecm_15digit_curves = 0;
	fwork->ecm_20digit_curves = 0;
	fwork->ecm_25digit_curves = 0;
	fwork->ecm_30digit_curves = 0;
	fwork->ecm_35digit_curves = 0;
	fwork->ecm_40digit_curves = 0;
	fwork->ecm_45digit_curves = 0;
	fwork->ecm_50digit_curves = 0;
	fwork->ecm_55digit_curves = 0;
	fwork->ecm_60digit_curves = 0;
	fwork->ecm_65digit_curves = 0;

    
	// preload work structure with curves appropriate to the amount
	// of specified initial work
	if (fwork->initial_work >= 60.0)
	{
		fwork->pp1_max_lvl1_curves = 0;
		fwork->pp1_max_lvl2_curves = 0;
		fwork->pp1_max_lvl3_curves = 0;
		fwork->pm1_max_lvl1_curves = 0;
		fwork->pm1_max_lvl2_curves = 0;
		fwork->pm1_max_lvl3_curves = 0;
		fwork->ecm_max_15digit_curves = 0;
		fwork->ecm_max_20digit_curves = 0;
		fwork->ecm_max_25digit_curves = 0;		
		fwork->ecm_max_30digit_curves = 0;
		fwork->ecm_max_35digit_curves = 0;
		fwork->ecm_max_40digit_curves = 0;
		fwork->ecm_max_45digit_curves = 0;
		fwork->ecm_max_50digit_curves = 0;
		fwork->ecm_max_55digit_curves = 0;
		fwork->ecm_60digit_curves = get_max_ecm_curves(fwork, state_ecm_60digit);				
		interp_state = state_ecm_65digit;
	}
    else if (fwork->initial_work >= 55.0)
	{
		fwork->pp1_max_lvl1_curves = 0;
		fwork->pp1_max_lvl2_curves = 0;
		fwork->pp1_max_lvl3_curves = 0;
		fwork->pm1_max_lvl1_curves = 0;
		fwork->pm1_max_lvl2_curves = 0;
		fwork->pm1_max_lvl3_curves = 0;
		fwork->ecm_max_15digit_curves = 0;
		fwork->ecm_max_20digit_curves = 0;
		fwork->ecm_max_25digit_curves = 0;		
		fwork->ecm_max_30digit_curves = 0;
		fwork->ecm_max_35digit_curves = 0;
		fwork->ecm_max_40digit_curves = 0;
		fwork->ecm_max_45digit_curves = 0;
		fwork->ecm_max_50digit_curves = 0;
		fwork->ecm_55digit_curves = get_max_ecm_curves(fwork, state_ecm_55digit);				
		interp_state = state_ecm_60digit;
	}
    else if (fwork->initial_work >= 50.0)
	{
		fwork->pp1_max_lvl1_curves = 0;
		fwork->pp1_max_lvl2_curves = 0;
		fwork->pp1_max_lvl3_curves = 0;
		fwork->pm1_max_lvl1_curves = 0;
		fwork->pm1_max_lvl2_curves = 0;
		fwork->pm1_max_lvl3_curves = 0;
		fwork->ecm_max_15digit_curves = 0;
		fwork->ecm_max_20digit_curves = 0;
		fwork->ecm_max_25digit_curves = 0;		
		fwork->ecm_max_30digit_curves = 0;
		fwork->ecm_max_35digit_curves = 0;
		fwork->ecm_max_40digit_curves = 0;
		fwork->ecm_max_45digit_curves = 0;
		fwork->ecm_50digit_curves = get_max_ecm_curves(fwork, state_ecm_50digit);			
		interp_state = state_ecm_55digit;
	}
    else if (fwork->initial_work >= 45.0)
	{
		fwork->pp1_max_lvl1_curves = 0;
		fwork->pp1_max_lvl2_curves = 0;
		fwork->pp1_max_lvl3_curves = 0;
		fwork->pm1_max_lvl1_curves = 0;
		fwork->pm1_max_lvl2_curves = 0;
		fwork->pm1_max_lvl3_curves = 0;
		fwork->ecm_max_15digit_curves = 0;
		fwork->ecm_max_20digit_curves = 0;
		fwork->ecm_max_25digit_curves = 0;		
		fwork->ecm_max_30digit_curves = 0;
		fwork->ecm_max_35digit_curves = 0;
		fwork->ecm_max_40digit_curves = 0;
		fwork->ecm_45digit_curves = get_max_ecm_curves(fwork, state_ecm_45digit);				
		interp_state = state_ecm_50digit;
	}
    else if (fwork->initial_work >= 40.0)
	{
		fwork->pp1_max_lvl1_curves = 0;
		fwork->pp1_max_lvl2_curves = 0;
		fwork->pp1_max_lvl3_curves = 0;
		fwork->pm1_max_lvl1_curves = 0;
		fwork->pm1_max_lvl2_curves = 0;
		fwork->pm1_max_lvl3_curves = 0;
		fwork->ecm_max_15digit_curves = 0;
		fwork->ecm_max_20digit_curves = 0;
		fwork->ecm_max_25digit_curves = 0;		
		fwork->ecm_max_30digit_curves = 0;
		fwork->ecm_max_35digit_curves = 0;
		fwork->ecm_40digit_curves = get_max_ecm_curves(fwork, state_ecm_40digit);		
		interp_state = state_ecm_45digit;
	}
    else if (fwork->initial_work >= 35.0)
	{
		fwork->pp1_max_lvl1_curves = 0;
		fwork->pp1_max_lvl2_curves = 0;
		fwork->pp1_max_lvl3_curves = 0;
		fwork->pm1_max_lvl1_curves = 0;
		fwork->pm1_max_lvl2_curves = 0;
		fwork->pm1_max_lvl3_curves = 0;
		fwork->ecm_max_15digit_curves = 0;
		fwork->ecm_max_20digit_curves = 0;
		fwork->ecm_max_25digit_curves = 0;		
		fwork->ecm_max_30digit_curves = 0;
		fwork->ecm_35digit_curves = get_max_ecm_curves(fwork, state_ecm_35digit);
		interp_state = state_ecm_40digit;
	}
    else if (fwork->initial_work >= 30.0)
	{
		fwork->pp1_max_lvl1_curves = 0;
		fwork->pp1_max_lvl2_curves = 0;
		fwork->pm1_max_lvl1_curves = 0;
		fwork->pm1_max_lvl2_curves = 0;
		fwork->ecm_max_15digit_curves = 0;
		fwork->ecm_max_20digit_curves = 0;
		fwork->ecm_max_25digit_curves = 0;
		fwork->ecm_30digit_curves = get_max_ecm_curves(fwork, state_ecm_30digit);
		interp_state = state_ecm_35digit;
	}
    else if (fwork->initial_work >= 25.0)
	{
		fwork->pp1_max_lvl1_curves = 0;
		fwork->pm1_max_lvl1_curves = 0;
		fwork->ecm_max_15digit_curves = 0;
		fwork->ecm_max_20digit_curves = 0;
		fwork->ecm_25digit_curves = get_max_ecm_curves(fwork, state_ecm_25digit);
		interp_state = state_ecm_30digit;
	}
    else if (fwork->initial_work >= 20.0)
	{
		fwork->pp1_max_lvl1_curves = 0;
		fwork->pm1_max_lvl1_curves = 0;
		fwork->ecm_max_15digit_curves = 0;
		fwork->ecm_20digit_curves = get_max_ecm_curves(fwork, state_ecm_20digit);
		interp_state = state_ecm_25digit;
	}
    else
    {
        interp_state = state_idle;
    }

    if (interp_state != state_idle)
    {

        // initializing with an indicated amount of work.  we are using this function
        // to try to figure out how many curves at the current level this is
        // equivalent too.  But the function is normally used to compute the number
        // of curves left toward the target, which could be either a digit level or
        // a pretest value.  Conclusion: temporarily set any pretest value to zero,
        // so that this function finds the equivalent curves to the indicated work.
        double tmp_pretest = fobj->autofact_obj.only_pretest;

        fobj->autofact_obj.only_pretest = 0;
        interp_and_set_curves(fwork, fobj, interp_state, fwork->initial_work,
            fwork->initial_work, 0);
            
        // restore any pretest value
        fobj->autofact_obj.only_pretest = tmp_pretest;
            
        // then fill in the equivalent curves to the indicated work amount.
        switch (interp_state)
        {
        case state_ecm_25digit: fwork->ecm_25digit_curves = fwork->curves; break;
        case state_ecm_30digit: fwork->ecm_30digit_curves = fwork->curves; break;
        case state_ecm_35digit: fwork->ecm_35digit_curves = fwork->curves; break;
        case state_ecm_40digit: fwork->ecm_40digit_curves = fwork->curves; break;
        case state_ecm_45digit: fwork->ecm_45digit_curves = fwork->curves; break;
        case state_ecm_50digit: fwork->ecm_50digit_curves = fwork->curves; break;
        case state_ecm_55digit: fwork->ecm_55digit_curves = fwork->curves; break;
        case state_ecm_60digit: fwork->ecm_60digit_curves = fwork->curves; break;
        case state_ecm_65digit: fwork->ecm_65digit_curves = fwork->curves; break;
        }

    }
	
	return;
}

void factor(fact_obj_t *fobj)
{
	//run a varity of factoring algorithms on b.
	//return any composite number left over.
	//the factoring routines will build up a list of factors.

	mpz_t b, origN, copyN;
	enum factorization_state fact_state;
	factor_work_t fwork;
	FILE *flog;
	struct timeval start, stop;
	double t_time;
	int user_defined_ecm_b2 = fobj->ecm_obj.stg2_is_default;
	int user_defined_pp1_b2 = fobj->pp1_obj.stg2_is_default;
	int user_defined_pm1_b2 = fobj->pm1_obj.stg2_is_default;
	FILE *data;
	char tmpstr[GSTR_MAXSIZE];
	int quit_after_sieve_method = 0;
    double initial_work = fobj->autofact_obj.initial_work;

	//factor() always ignores user specified B2 values
	fobj->ecm_obj.stg2_is_default = 1;
	fobj->pp1_obj.stg2_is_default = 1;
	fobj->pm1_obj.stg2_is_default = 1;

	mpz_init(origN);
	mpz_init(copyN);
	mpz_init(b);

	mpz_set(origN, fobj->N);
	mpz_set(copyN, origN);
	mpz_set(b, origN);

	if (mpz_cmp_ui(b,1) <= 0)
	{
		mpz_clear(copyN);
		mpz_clear(origN);
		mpz_clear(b);
		return;
	}	
	
	gettimeofday(&start, NULL);

	flog = fopen(fobj->flogname,"a");
	logprint(flog,"\n");
	logprint(flog,"****************************\n");
	{
		char *s;
		s = mpz_get_str(NULL, 10, b);
		// use ... when we have very big numbers?
		logprint(flog,"Starting factorization of %s\n", s);
		free(s);
	}
	logprint(flog,"using pretesting plan: %s\n",fobj->autofact_obj.plan_str);
    if (fobj->autofact_obj.yafu_pretest_plan == PRETEST_CUSTOM)
    {
        logprint(flog, "custom pretest ratio is: %1.4f\n",
            fobj->autofact_obj.target_pretest_ratio);
    }
    if (fobj->autofact_obj.only_pretest > 1)
    {
        logprint(flog, "custom pretesting limit is: %d\n",
            fobj->autofact_obj.only_pretest);
    }

	if (check_tune_params(fobj))
	{
        if (fobj->autofact_obj.prefer_xover)
        {
            logprint(flog, "using specified qs/gnfs crossover of %1.0f digits\n",
                fobj->autofact_obj.qs_gnfs_xover);
            logprint(flog, "using specified qs/snfs crossover of %1.0f digits\n",
                fobj->autofact_obj.qs_snfs_xover);
        }
        else
        {
            logprint(flog, "using tune info for qs/gnfs crossover\n");
        }
	}
    else
    {
        logprint(flog, "no tune info: using qs/gnfs crossover of %1.0f digits\n",
            fobj->autofact_obj.qs_gnfs_xover);
        logprint(flog, "no tune info: using qs/snfs crossover of %1.0f digits\n",
            fobj->autofact_obj.qs_snfs_xover);
    }

    // if the user input a scaling factor rather than a digit level
    // then compute the effective digit number for this input.
    if (fobj->autofact_obj.initial_work < 1.0)
    {
        initial_work = fobj->autofact_obj.initial_work * mpz_sizeinbase(fobj->N, 10);
    }

    // put the initial work done into the work structure.  It's important
    // to not modify the autofact_obj.initial_work element because if it was
    // input as a scaling factor then modifying it would destroy the scaling
    // factor for other inputs on this run (potentially a batch job).
    fwork.initial_work = initial_work;

    if (initial_work > 0.0)
    {       
        logprint(flog, "input indicated to have been pretested to t%1.2f\n",
            initial_work);
    }

	logprint(flog,"****************************\n");
	fclose(flog);

	fobj->autofact_obj.autofact_active = 1;

	if (VFLAG >= 0)
	{
		gmp_printf("fac: factoring %Zd\n",b);
		printf("fac: using pretesting plan: %s\n",fobj->autofact_obj.plan_str);
		if (fobj->autofact_obj.yafu_pretest_plan == PRETEST_CUSTOM)
			printf("fac: custom pretest ratio is: %1.4f\n",fobj->autofact_obj.target_pretest_ratio);
		if (fobj->autofact_obj.only_pretest > 1)
			printf("fac: custom pretesting limit is: %d\n",fobj->autofact_obj.only_pretest);
		if (check_tune_params(fobj))
		{
            if (fobj->autofact_obj.prefer_xover)
            {
                printf("fac: using specified qs/gnfs crossover of %1.0f digits\n",
                    fobj->autofact_obj.qs_gnfs_xover);
                printf("fac: using specified qs/snfs crossover of %1.0f digits\n",
                    fobj->autofact_obj.qs_snfs_xover);
            }
            else
            {
                printf("fac: using tune info for qs/gnfs crossover\n");
            }
		}
        else
        {
            printf("fac: no tune info: using qs/gnfs crossover of %1.0f digits\n",
                fobj->autofact_obj.qs_gnfs_xover);
            printf("fac: no tune info: using qs/snfs crossover of %1.0f digits\n",
                fobj->autofact_obj.qs_snfs_xover);
        }

        if (initial_work > 0.0)
        {
            printf("fac: input indicated to have been pretested to t%1.2f\n",
                initial_work);
        }

	}	

	init_factor_work(&fwork, fobj);

	//starting point of factorization effort
	fact_state = state_trialdiv;

	//check to see if a siqs savefile exists for this input	
	data = fopen(fobj->qs_obj.siqs_savefile,"r");

	if (data != NULL)
	{	
		char *substr;
		mpz_t tmpz;
		mpz_t g;

		//read in the number from the savefile
		mpz_init(tmpz);
		mpz_init(g);

		fgets(tmpstr,1024,data);
		substr = tmpstr + 2;
		mpz_set_str(tmpz, substr, 0);	//auto detect the base

		if (resume_check_input_match(tmpz, b, g))
		{
			if (VFLAG > 0)
				printf("fac: found siqs savefile, resuming siqs\n");

			// remove any common factor so the input exactly matches
			// the file
			// mpz_tdiv_q(b, b, g);
			// mpz_set(fobj->N, b);
			// mpz_set(origN, b);
			// mpz_set(copyN, b);

			//override default choice
			fact_state = state_qs;

			// if for some reason qs doesn't find factors (such as
			// a user specified time out), don't continue ecm-ing, etc.
			quit_after_sieve_method = 1;
		}
		mpz_clear(tmpz);
		mpz_clear(g);
		fclose(data);
	}

	//check to see if a nfs job file exists for this input	
	data = fopen(fobj->nfs_obj.job_infile,"r");

	if (data != NULL)
	{	
		char *substr;
		mpz_t tmpz;
		mpz_t g;

		//read in the number from the job file
		mpz_init(tmpz);
		mpz_init(g);

		// may not be sufficient, in extreme cases...
		fgets(tmpstr,1024,data);
		substr = tmpstr + 2;
		mpz_set_str(tmpz, substr, 0);	//auto detect the base

		if (resume_check_input_match(tmpz, b, g))
		{
            // check if this is a snfsable number.  If the input is
            // really small and we don't check this, the resume
            // may default back to siqs.
            if ((mpz_sizeinbase(b,10) >= fobj->autofact_obj.qs_snfs_xover) && 
                (fobj->autofact_obj.has_snfs_form < 0))
            {
                snfs_t *poly;
                mpz_set(fobj->nfs_obj.gmp_n, b);
#ifdef USE_NFS
                poly = snfs_find_form(fobj);

                if (poly != NULL)
                {
                    fobj->autofact_obj.has_snfs_form = 1;
                    // The actual poly is not needed now, so just get rid of it.
                    snfs_clear(poly);
                    free(poly);

                    if (VFLAG > 0)
                        printf("fac: found nfs job file and snfs form, resuming snfs\n");
                }
                else
                {
                    fobj->autofact_obj.has_snfs_form = 0;
                    if (VFLAG > 0)
                        printf("fac: found nfs job file, resuming nfs\n");
                }
#else
                fobj->autofact_obj.has_snfs_form = 0;
#endif
            }
            else
            {
                if (VFLAG > 0)
                    printf("fac: found nfs job file, resuming nfs\n");
            }		

			// remove any common factor so the input exactly matches
			// the file
			mpz_tdiv_q(b, b, g);
			mpz_set(fobj->N, b);
			mpz_set(origN, b);
			mpz_set(copyN, b);

			//override default choice
			fact_state = state_nfs;

			// if for some reason nfs doesn't find factors (such as
			// a user specified time out or -ns, -nc, etc.), 
			// don't continue ecm-ing, etc.
			quit_after_sieve_method = 1;
		}
		mpz_clear(tmpz);
		mpz_clear(g);
		fclose(data);
	}

	// state machine to factor the number using a variety of methods
	while (fact_state != state_done)
	{	
        // do the next item of work
		do_work(fact_state, &fwork, b, fobj);

        // check if we are done:
        // * number is completely factored
        // * sieve method was performed and either finished or was interrupted.
        if (check_if_done(fobj, origN) ||
            (quit_after_sieve_method &&
            ((fact_state == state_qs) ||
            (fact_state == state_nfs))) ||
            ((fact_state == state_nfs) &&
            (fobj->flags == FACTOR_INTERRUPT)))
        {
            fact_state = state_done;
        }
        else if ((fact_state >= state_ecm_15digit) && (fact_state <= state_ecm_65digit))
        {
            // if we ran ecm, check the ecm exit code and
            // handle appropriately.
            if (fobj->ecm_obj.exit_cond == ECM_EXIT_ABORT)
            {
                FILE *flog;
                double work_done;

                flog = fopen(fobj->flogname, "a");
                work_done = compute_ecm_work_done(&fwork, 1, flog);
                logprint(flog, "\testimated sum of completed work is t%1.2f\n", work_done);
                fclose(flog);
                fact_state = state_done;
            }
        }

        // if not done, figure out the next item of work.
        if (fact_state != state_done)
        {
            fact_state = schedule_work(&fwork, b, fobj);
        }
	}

	// optionally record output in one or more file formats
	if (fobj->num_factors >= 1) 
	{
		//If the only factor in our array == N, then N is prime or prp...
		if (fobj->autofact_obj.want_output_primes && (mpz_cmp(fobj->fobj_factors[0].factor,origN) == 0))
		{
			if ((fobj->autofact_obj.op_file = fopen(fobj->autofact_obj.op_str, "a")) == NULL)
				printf(" ***Error: unable to open %s\n", fobj->autofact_obj.op_str);
			else
			{
				if (fobj->autofact_obj.want_output_expressions)
					fprintf(fobj->autofact_obj.op_file, "%s\n", fobj->str_N.s);
				else
					gmp_fprintf(fobj->autofact_obj.op_file, "%Zd\n", origN);
				if (fclose(fobj->autofact_obj.op_file) != 0)
					printf(" ***Error: problem closing file %s\n", fobj->autofact_obj.op_str);
			}
		}

		//If the first factor in the array != N, then is composite and we have factors...
		if (fobj->autofact_obj.want_output_factors && (mpz_cmp(fobj->fobj_factors[0].factor,origN) != 0))
		{
			if ((fobj->autofact_obj.of_file = fopen(fobj->autofact_obj.of_str, "a")) == NULL)
				printf(" ***Error: unable to open %s\n", fobj->autofact_obj.of_str);
			else
			{
				int i;
				//fprintf(fobj->autofact_obj.of_file, "%s\n", z2decstr(&origN,&gstr1));
				if (fobj->autofact_obj.want_output_expressions)
					fprintf(fobj->autofact_obj.of_file, "(%s)", fobj->str_N.s);
				else
					gmp_fprintf(fobj->autofact_obj.of_file, "%Zd", origN);
				for (i=0; i<fobj->num_factors; i++)
				{
					gmp_fprintf(fobj->autofact_obj.of_file, "/%Zd", fobj->fobj_factors[i].factor);
					if (fobj->fobj_factors[i].count > 1)
						fprintf(fobj->autofact_obj.of_file, "^%d", fobj->fobj_factors[i].count);
					//fprintf(fobj->autofact_obj.of_file, "\n");
				}
				fprintf(fobj->autofact_obj.of_file,"\n");
				if (fclose(fobj->autofact_obj.of_file) != 0)
					printf(" ***Error: problem closing file %s\n", fobj->autofact_obj.of_str);
			}
		}
	}
	else //assume: composite with no known factors... (need to clarify)
	{
		if (fobj->autofact_obj.want_output_unfactored)
		{
			if ((fobj->autofact_obj.ou_file = fopen(fobj->autofact_obj.ou_str, "a")) == NULL)
				printf(" ***Error: unable to open %s\n", fobj->autofact_obj.ou_str);
			else
			{
				if (fobj->autofact_obj.want_output_expressions)
					fprintf(fobj->autofact_obj.ou_file, "%s\n", fobj->str_N.s);
				else
					gmp_fprintf(fobj->autofact_obj.ou_file, "%s\n", origN);
				if (fclose(fobj->autofact_obj.ou_file) != 0)
					printf(" ***Error: problem closing file %s\n", fobj->autofact_obj.ou_str);
			}
		}
	}

		
	if (mpz_cmp_ui(b, 1) != 0)
	{
		if (is_mpz_prp(b))
		{
			flog = fopen(fobj->flogname,"a");
			logprint(flog,"prp%d cofactor = %s\n",gmp_base10(b),
				mpz_conv2str(&gstr1.s, 10, b));
			fclose(flog);
			add_to_factor_list(fobj,b);
		}
		else
		{
			flog = fopen(fobj->flogname,"a");
			logprint(flog,"c%d cofactor = %s\n",gmp_base10(b),
				mpz_conv2str(&gstr1.s, 10, b));
			fclose(flog);
			add_to_factor_list(fobj,b);
		}
	}
    
	mpz_set(fobj->N, b);

	gettimeofday (&stop, NULL);
    t_time = yafu_difftime(&start, &stop);

	if (VFLAG >= 0)
		printf("Total factoring time = %6.4f seconds\n",t_time);

	flog = fopen(fobj->flogname,"a");
	if (flog == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("Could not open %s for appending\n",fobj->flogname);
	}
	else
	{
		logprint(flog,"Total factoring time = %6.4f seconds\n",t_time);
		fclose(flog);
	}

	fobj->autofact_obj.autofact_active=0;

	//restore flags
	fobj->ecm_obj.stg2_is_default = user_defined_ecm_b2;
	fobj->pp1_obj.stg2_is_default = user_defined_pp1_b2;
	fobj->pm1_obj.stg2_is_default = user_defined_pm1_b2;    

	mpz_clear(origN);
	mpz_clear(copyN);
	mpz_clear(b);
	return;
}

void spfactorlist(uint64 nstart, uint64 nrange)
{
	int i, c, s, e, m;
	double t;
	struct timeval tstart, tstop;
	mpz_t tmp;
	fact_obj_t f;
	uint64 sf;
	uint64 tflimit = (uint64)sqrt(sqrt(nstart + nrange));

	mpz_init(tmp);
	init_factobj(&f);

	gettimeofday (&tstart, NULL);
	c = 0;
	s = 0;
	e = 0;
	m = 0;

	for (i=0; i<nrange; i++)
	{
		uint64 n = nstart + i;
		int j;
			
		// factors of 2
		while (!(n & 0x1))
			n >>= 1;

		j = 1;
		while ((n > 1) && (j < NUM_P) && (PRIMES[j] < tflimit))
		{
			if (n % PRIMES[j] == 0)
				n /= PRIMES[j];
			else
				j++;
		}

		if (n == 1)
		{ 
			c++;
			continue;
		}
			
		mpz_set_64(tmp, n);
		if (mpz_probab_prime_p(tmp, 3))
		{
			c++;
			continue;
		}

		sf = sp_shanks_loop(tmp, &f);
			
		if (sf > 1)
		{
			// squfof found a factor, divide it out
			n /= sf;

			if (n > 1)
			{
				// check to see if residue is prime
				mpz_set_64(tmp, n);
				if (mpz_probab_prime_p(tmp, 3))
				{
					c++;
					s++;
					continue;
				}
				else
				{
					// n is still composite.  try squfof again.
					mpz_set_64(tmp, n);
					n /= sp_shanks_loop(tmp, &f);
					if (n > 1)
					{
						mpz_set_64(tmp, n);
						if (mpz_probab_prime_p(tmp, 3))
						{
							c++;
							s++;
							continue;
						}
					}
				}
			}			
		}
		// else, try fermat
		sf = spfermat(n, tflimit);

		if (sf > 1)
		{ 
			n /= sf;

			if (n > 1)
			{
				// check to see if residue is prime
				mpz_set_64(tmp, n);
				if (mpz_probab_prime_p(tmp, 3))
				{
					m++;
					c++;
					continue;
				}
			}			
		}

		// else, try rho
		mpz_set_64(f.rho_obj.gmp_n, n);
		brent_loop(&f);
		e++;
		c++;
		clear_factor_list(&f);
	}

	free_factobj(&f);
	mpz_clear(tmp);
	gettimeofday (&tstop, NULL);
	t = yafu_difftime (&tstart, &tstop);

	printf("completed %d factorizations (%d SQUFOF, %d fermat, %d rho) "
		"in %6.4f seconds\n", c, s, m, e, t);

	return;
}

