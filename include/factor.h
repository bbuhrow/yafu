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

#ifndef _FACTOR_H_
#define _FACTOR_H_

#include "yafu.h"
#include "arith.h"
#include "ytools.h"
#include "monty.h"
#include "cmdOptions.h"
#include "gmp.h"
#include "msieve_common.h"

//#define NO_ZLIB



#if defined( _MSC_VER )
#define USE_NFS
#endif

/* declarations of factorization routines not grouped in with QS methods
or group theoretic methods, as well as higher level factorization
routines and bookkeeping */

// factorization objects //

typedef struct
{
	mpz_t gmp_n;
	mpz_t gmp_f;
	uint32 B1;
	uint64 B2;
	int stg2_is_default;
	double ttime;
	uint32 base;				//we compute base^B1

	// fit parameters to compute time_per_curve as a function of B1
	double pm1_exponent;
	double pm1_multiplier;
	double pm1_tune_freq;
} pm1_obj_t;

typedef struct
{
	mpz_t gmp_n;
	mpz_t gmp_f;
	uint32 B1;
	uint64 B2;
	int stg2_is_default;
	double ttime;
	uint32 base;				//we compute base^B1
	uint32 numbases;

	// fit parameters to compute time_per_curve as a function of B1
	double pp1_exponent;
	double pp1_multiplier;
	double pp1_tune_freq;

    // RNG state for picking bases
    uint64 lcg_state;

} pp1_obj_t;

typedef struct
{
	mpz_t gmp_n;
	mpz_t gmp_f;

    int save_b1;
    int prefer_gmpecm;
	char ecm_path[1024];
	int use_external;
	uint64 B1;
	uint64 B2;
	int stg2_is_default;
	int curves_run;
	uint32 num_curves;
	uint32 sigma;				//sigma value of successful curve
	uint32 num_factors;			//number of factors found in this method
	double ttime;
	uint64 ecm_ext_xover;
    int bail_on_factor;
    enum ecm_exit_cond_e exit_cond;              // exit condition

	// fit parameters to compute time_per_curve as a function of B1
	double ecm_exponent;
	double ecm_multiplier;
    double ecm_tune_freq;

    uint32 rand_seed1;
    uint32 rand_seed2;

    // RNG state for each thread
    uint64 *lcg_state;

} ecm_obj_t;

typedef struct
{
	mpz_t gmp_n;
	mpz_t gmp_f;
	uint32 iterations;
	uint32 num_poly;
	uint32 *polynomials;
	uint32 curr_poly;			//current polynomial in the list of polynomials
	double ttime;

} rho_obj_t;

typedef struct
{
	mpz_t gmp_n;
	mpz_t gmp_f;
	uint32 limit;				//trial div limit
	uint32 fmtlimit;			//fermat max iterations
	uint32 num_factors;			//number of factors found in this method
	double ttime;
	int print;

} div_obj_t;

typedef struct
{
	mpz_t gmp_n;
	mpz_t gmp_f;
	uint32 num_factors;			//number of factors found in this method
	double ttime;

} squfof_obj_t;

typedef struct
{
	mpz_t gmp_n;
	mpz_t N; // full form used for snfs

	int snfs; // if this is a snfs job
	int gnfs; // user wants gnfs
	mpz_t snfs_cofactor;
	int pref_degree;
	int alt_degree;

	char ggnfs_dir[GSTR_MAXSIZE];
	uint32 siever;
	int sq_side;
	uint32 startq;
	uint32 rangeq;
	uint32 polystart;
	uint32 polyrange;
	double filter_min_rels_nudge;
	char outputfile[GSTR_MAXSIZE];
	char logfile[GSTR_MAXSIZE];
	char fbfile[GSTR_MAXSIZE];
	uint32 timeout;
	char job_infile[GSTR_MAXSIZE];
	int poly_option;
	int restart_flag;
	uint32 polybatch;
	uint32 nfs_phases;
	uint32 snfs_testsieve_threshold;

	double gnfs_exponent;
	double gnfs_multiplier;
	double gnfs_tune_freq;

	uint32 min_digits;

	uint32 num_factors;			//number of factors found in this method
	double ttime;
	
	// an object used to carry around information needed by the msieve library
	msieve_obj *mobj;

	char filearg[GSTR_MAXSIZE]; // used to facilitate external trial sieving

} nfs_obj_t;


//globals for implementing the "plan" and "pretest" switches
enum pretest_plan {
	PRETEST_NONE = 0,
	PRETEST_NOECM = 1,
	PRETEST_LIGHT = 2,
	PRETEST_NORMAL = 3,
	PRETEST_DEEP = 4,
	PRETEST_CUSTOM = 5,
};

typedef struct
{
	//crossover between qs and gnfs
	double qs_gnfs_xover;
    //crossover between qs and snfs
    double qs_snfs_xover;
	int prefer_xover;

	//balance of ecm and various sieve methods
	double target_pretest_ratio;
	int has_snfs_form;			// 1 = input has snfs form
								// 0 = input does not have snfs form
								// -1 = input has not been checked yet

	//double target_ecm_gnfs_ratio;
	//double target_ecm_snfs_ratio;

	int no_ecm;

	int want_only_1_factor;
	int want_output_primes;
	int want_output_factors;
	int want_output_unfactored;
	int want_output_expressions;
	FILE *op_file;
	FILE *of_file;
	FILE *ou_file;
	char op_str[1024];
	char of_str[1024];
	char ou_str[1024];

	enum pretest_plan yafu_pretest_plan;
	char plan_str[1024];
	int only_pretest;
	int autofact_active;

	// user supplied value indicating prior pretesting work
	double initial_work;

	double ttime;

} autofact_obj_t;

typedef struct
{
	mpz_t N;					//numerical representation of input
	uint32 digits;				//number of digits in input
	uint32 bits;				//number of bits in input
	char flogname[1024];		//name of the factorization logfile to use
	FILE *logfile;				//the logfile
	char savefile_name[80];		//data savefile name
	uint32 flags;				//state flags

	// info for work done in various places during this factorization
	div_obj_t div_obj;			//info for any trial division work 
	squfof_obj_t squfof_obj;	//info for any squfof work 
	rho_obj_t rho_obj;			//info for any rho work 
	pm1_obj_t pm1_obj;			//info for any pm1 work 
	pp1_obj_t pp1_obj;			//info for any pp1 work 
	ecm_obj_t ecm_obj;			//info for any ecm work 
	qs_obj_t qs_obj;			//info for any qs work 
	nfs_obj_t nfs_obj;			//info for any nfs work
	autofact_obj_t autofact_obj;

	uint32 cache_size1;
	uint32 cache_size2;
	int num_threads;
	uint32 seed1;
	uint32 seed2;

	// threshold at which we know number are prime, as determined by trial division
	uint64 prime_threshold;

	// manage a list of factors
	yfactor_list_t *factors;
	int do_logging;
	int refactor_depth;
    options_t* options;

    // configurable options
    int VFLAG;
    int THREADS;
    int LOGFLAG;
    int LATHREADS;
    int NUM_WITNESSES;

    // computer info
    double MEAS_CPU_FREQUENCY;
    char CPU_ID_STR[80];
    int HAS_SSE41;
    int HAS_AVX;
    int HAS_AVX2;
    int L1CACHE;
    int L2CACHE;
    int L3CACHE;
    
    // RNG state
    uint64_t lcg_state;

} fact_obj_t;

void init_factobj(fact_obj_t *fobj, options_t *options);
void free_factobj(fact_obj_t *fobj);
void reset_factobj(fact_obj_t *fobj);
void alloc_factobj(fact_obj_t *fobj);

//#if defined(WIN32)
// windows machines also need these declarations for functions located
// within common.lib or gnfs.lib, for NFS factorizations using msieve.
void savefile_init(savefile_t *s, char *filename);
void savefile_free(savefile_t *s);
void savefile_open(savefile_t *s, uint32 flags);
void savefile_close(savefile_t *s);
uint32 savefile_eof(savefile_t *s);
uint32 savefile_exists(savefile_t *s);
void savefile_read_line(char *buf, size_t max_len, savefile_t *s);
void savefile_write_line(savefile_t *s, char *buf);
void savefile_flush(savefile_t *s);
void savefile_rewind(savefile_t *s);

//#endif




/*-----------TOP LEVEL ENTRY POINT FOR ALL FACTORING ROUTINES ----------*/

uint32 test_qn_res[128];

uint64 spbrent(uint64 N, uint64 c, int imax);
uint64 spbrent64(uint64 N, int imax);
int mbrent(fact_obj_t *fobj);
int montybrent(monty_t *mdata, mpz_t n, mpz_t f, uint32 a, uint32 imax);
void brent_loop(fact_obj_t *fobj);
void pollard_loop(fact_obj_t *fobj);
void williams_loop(fact_obj_t *fobj);
int ecm_loop(fact_obj_t *fobj);
yfactor_t * vec_ecm_main(mpz_t N, uint32 numcurves, uint64 B1, 
    uint64 B2, int threads, int *numfactors, int verbose, 
    int save_b1, uint32 *curves_run);
void tinyecm(mpz_t n, mpz_t f, uint32 B1, uint32 B2, uint32 curves,
    uint64* lcg_state, int verbose);
void microecm(uint64 n, uint64 *f, uint32 B1, uint32 B2, uint32 curves, int verbose);
uint64 do_uecm(uint64 q);
uint64 sp_shanks_loop(mpz_t N, fact_obj_t *fobj);
uint64 LehmanFactor(uint64 N, double Tune, int DoTrial, double CutFrac);
void init_lehman();
void zTrial(fact_obj_t *fobj);
void zFermat(uint64 limit, uint32 mult, fact_obj_t *fobj);
void factor_perfect_power(fact_obj_t *fobj, mpz_t b);
void nfs(fact_obj_t *fobj);
int par_shanks_loop(uint64 *N, uint64 *f, int num_in);

int sptestsqr(uint64 n);
uint64 spfermat(uint64 n, uint64 limit);

//auto factor routine
void factor(fact_obj_t *fobj);

// factoring related utility
int resume_check_input_match(mpz_t file_n, mpz_t input_n, mpz_t common_fact, int VFLAG);

/* Factor a number using GNFS. Returns
   1 if any factors were found and 0 if not */

uint32 factor_gnfs(msieve_obj *obj, mp_t *n, factor_list_t *factor_list);

void factor_tune(fact_obj_t *fobj);

#endif //_FACTOR_H
