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

#include "arith.h"
#include "ytools.h"
#include "monty.h"
#include "gmp.h"
#include "msieve_common.h"

#ifdef __MINGW32__
#include <Windows.h>
#endif

#define GSTR_MAXSIZE 1024

//#define NO_ZLIB

#if defined( _MSC_VER )

#define USE_NFS
#define MySleep(x) Sleep((x))

#else

#include <unistd.h> // usleep
//sleep in milliseconds
#define MySleep(x) usleep((x)*1000)	

#endif

// these are similar to things msieve defines.  Differences:
// the factor type contains more info about how the factor
// was found and the list is not fixed size.

typedef struct
{
    mpz_t factor;
    int count;
    int type;
    // new parameters to support returning factors
    // from avx-ecm.  In general I think it makes sense
    // to have this structure contain information not
    // only about the factor itself but how it was found.
    int method;             // factorization method used
    uint32_t curve_num;     // curve found
    uint64_t sigma;         // sigma value
    int tid;                // thread found
    int vid;                // vector position found
} yfactor_t;

typedef struct
{
    yfactor_t* factors;
    int num_factors;
    int alloc_factors;

    // the primitude of factors can be proven by aprcl: these
    // determine when and what is displayed while proving.
    int aprcl_prove_cutoff;
    int aprcl_display_cutoff;
} yfactor_list_t;


/* structure encapsulating the savefile used in a factorization */
typedef struct {

#if defined(WIN32) || defined(_WIN64)
    HANDLE file_handle;
    uint32_t read_size;
    uint32_t eof;
#else
    FILE* fp;
#endif
    char* name;
    char* buf;
    uint32_t buf_off;
} qs_savefile_t;


/*--------------DECLARATIONS FOR MANAGING FACTORS FOUND  and factor I/O -----------------*/

int add_to_factor_list(yfactor_list_t* factors, mpz_t n,
    int VFLAG, int NUM_WITNESSES);
void print_factors(yfactor_list_t* fobj, mpz_t N, int VFLAG, int NUM_WITNESSES);
void clear_factor_list(yfactor_list_t* fobj);
void delete_from_factor_list(yfactor_list_t* fobj, mpz_t n);

typedef struct
{
	mpz_t gmp_n;
	mpz_t gmp_f;
	uint64_t B1;
	uint64_t B2;
	int stg2_is_default;
	double ttime;
	uint32_t base;				// we compute base^B1
    char vpm1_work_file[256];
    mpz_t vecn[8];
    int vecnum;
    char resume_file[256];

	// fit parameters to compute time_per_curve as a function of B1
	double pm1_exponent;
	double pm1_multiplier;
	double pm1_tune_freq;
} pm1_obj_t;

typedef struct
{
	mpz_t gmp_n;
	mpz_t gmp_f;
	uint64_t B1;
	uint64_t B2;
	int stg2_is_default;
	double ttime;
	uint32_t base;				// we compute base^B1
	uint32_t numbases;
    char vpp1_work_file[256];
    mpz_t vecn[8];
    int vecnum;
    char resume_file[256];

	// fit parameters to compute time_per_curve as a function of B1
	double pp1_exponent;
	double pp1_multiplier;
	double pp1_tune_freq;

    // RNG state for picking bases
    uint64_t lcg_state;

} pp1_obj_t;

typedef struct
{
	mpz_t gmp_n;
	mpz_t gmp_f;

    int save_b1;
    int prefer_gmpecm;
    int prefer_gmpecm_stg2;
    int prefer_avxecm_stg2;
	char ecm_path[1024];
	int use_external;
	uint64_t B1;
	uint64_t B2;
	int stg2_is_default;
	int curves_run;
	uint32_t num_curves;
	uint64_t sigma;				//sigma value of successful curve
	uint32_t num_factors;			//number of factors found in this method
	double ttime;
	uint64_t ecm_ext_xover;
    int bail_on_factor;
    enum ecm_exit_cond_e exit_cond;              // exit condition
    int gpucurves;
    int use_gpudev;
    int use_gpuecm;
    int use_cgbn;

	// fit parameters to compute time_per_curve as a function of B1
	double ecm_exponent;
	double ecm_multiplier;
    double ecm_tune_freq;

    uint32_t rand_seed1;
    uint32_t rand_seed2;

    // RNG state for each thread
    uint64_t *lcg_state;

} ecm_obj_t;

typedef struct
{
	mpz_t gmp_n;
	mpz_t gmp_f;
	uint32_t iterations;
	uint32_t num_poly;
	uint32_t *polynomials;
	uint32_t curr_poly;			//current polynomial in the list of polynomials
	double ttime;

} rho_obj_t;

typedef struct
{
	mpz_t gmp_n;
	mpz_t gmp_f;
	uint32_t limit;				//trial div limit
	uint32_t fmtlimit;			//fermat max iterations
	uint32_t num_factors;			//number of factors found in this method
	double ttime;
	int print;

} div_obj_t;

typedef struct
{
	mpz_t gmp_n;
	mpz_t gmp_f;
	uint32_t num_factors;			//number of factors found in this method
	double ttime;

} squfof_obj_t;

typedef struct
{
    mpz_t gmp_n;
    char* savefile_name;
    char* flogname;
    qs_savefile_t savefile;		//savefile object
    FILE* logfile;
    char siqs_savefile[1024];

    double qs_exponent;
    double qs_multiplier;
    double qs_tune_freq;
    int no_small_cutoff_opt;	//1 is true - perform no optimization.  0 = optimize.

    // siqs parameters to override defaults
    uint32_t gbl_override_B;			//override the # of factor base primes
    uint32_t gbl_override_small_cutoff;			//override the tf_small_cutoff value
    uint32_t gbl_override_tf;			//extra reduction of the TF bound by X bits
    uint32_t gbl_override_time;		//stop after this many seconds
    uint32_t gbl_override_rel;		//stop after collecting this many relations
    uint32_t gbl_override_blocks;		//override the # of blocks used
    uint32_t gbl_override_lpmult;		//override the large prime multiplier
    float gbl_override_bdiv;        // override the lpmax divider for batch GCD
    uint32_t gbl_btarget;             // the target number of batch relations
    uint32_t gbl_override_lpb;		// override the large prime bound (specified in bits)
    double gbl_override_mfbt;		// override the mfbt exponent
    double gbl_override_mfbd;		// override the mfbd exponent
    uint32_t gbl_override_3lp_bat;    // don't do 3lp batch factoring (default is do_batch)
    uint32_t inmem_cutoff;          // below X digits, don't use savefile
    int gbl_force_DLP;
    int gbl_force_TLP;

    int gbl_override_B_flag;
    int gbl_override_small_cutoff_flag;
    int gbl_override_tf_flag;
    int gbl_override_time_flag;
    int gbl_override_rel_flag;
    int gbl_override_blocks_flag;
    int gbl_override_lpmult_flag;
    int gbl_override_bdiv_flag;

    uint32_t num_factors;			//number of factors found in this method
    uint32_t flags;				//each bit corresponds to a location in the 
                                //flags enum
    double rels_per_sec;
    double qs_time;
    double total_time;

    // configurable options
    int VFLAG;
    int THREADS;
    int LOGFLAG;
    int LATHREADS;
    int NUM_WITNESSES;

    // computer info
    double MEAS_CPU_FREQUENCY;
    char CPU_ID_STR[80];
    int HAS_SSE2;
    int HAS_SSE41;
    int HAS_AVX;
    int HAS_AVX2;
    int L1CACHE;
    int L2CACHE;
    int L3CACHE;
    uint32_t cache_size2;
    uint32_t num_threads;

    uint32_t seed1;
    uint32_t seed2;
    uint64_t lcg_state;

    uint32_t bits;
    uint32_t digits;

} qs_obj_t;

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
	uint32_t siever;
	int sq_side;
	uint32_t startq;
	uint32_t rangeq;
	uint32_t polystart;
	uint32_t polyrange;
	double filter_min_rels_nudge;
	char outputfile[GSTR_MAXSIZE];
	char logfile[GSTR_MAXSIZE];
	char fbfile[GSTR_MAXSIZE];
	uint32_t timeout;
	char job_infile[GSTR_MAXSIZE];
	int poly_option;
	int restart_flag;
	uint32_t polybatch;
	uint32_t nfs_phases;
	uint32_t snfs_testsieve_threshold;

	double gnfs_exponent;
	double gnfs_multiplier;
	double gnfs_tune_freq;
    double poly_time;
    double sieve_time;
    double filter_time;
    double la_time;
    double sqrt_time;
	uint32_t min_digits;

	uint32_t num_factors;			//number of factors found in this method
	double ttime;
	
	// an object used to carry around information needed by the msieve library
	msieve_obj *mobj;

	char filearg[GSTR_MAXSIZE]; // used to facilitate external trial sieving

    uint32_t cadoMsieve;
    char cado_dir[GSTR_MAXSIZE];
    char convert_poly_path[GSTR_MAXSIZE];
} nfs_obj_t;

// enum for implementing the "plan" and "pretest" switches
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
    int json_pretty;
    int want_only_1_factor;
    int want_output_primes;
    int want_output_factors;
    int want_output_unfactored;
    int want_output_expressions;
    FILE* op_file;
    FILE* of_file;
    FILE* ou_file;
    char op_str[1024];
    char of_str[1024];
    char ou_str[1024];
    int stople;
    int stoplt;
    int stopgt;
    int stopge;
    int stopbase;
    int stopeq;
    int stopprime;
    int check_stop_conditions;

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
    mpz_t input_N;					// numerical representation of input
	mpz_t N;					    // numerical representation of input, can change while processing
	uint32_t digits;			    // number of digits in input
	uint32_t bits;				    // number of bits in input
	char flogname[1024];		    // name of the factorization logfile to use
	FILE *logfile;				    // the logfile
	char savefile_name[80];		    // data savefile name
    char factor_json_name[80];      // json file name
	uint32_t flags;				    // state flags
    char *input_str;                // a copy of the input string
    int input_str_alloc;
    char argc;                      // number of input arguments to yafu
    char** argv;                    // pointer to the input arguments to yafu (read only)

	// info for work done in various places during this factorization
	div_obj_t div_obj;			    // info for any trial division work 
	squfof_obj_t squfof_obj;	    // info for any squfof work 
	rho_obj_t rho_obj;			    // info for any rho work 
	pm1_obj_t pm1_obj;			    // info for any pm1 work 
	pp1_obj_t pp1_obj;			    // info for any pp1 work 
	ecm_obj_t ecm_obj;			    // info for any ecm work 
	qs_obj_t qs_obj;			    // info for any qs work 
	nfs_obj_t nfs_obj;			    // info for any nfs work
	autofact_obj_t autofact_obj;

	// list of primes.  To be filled in once if multiple
    // algorithms will need it.
    uint64_t* primes;
    uint64_t num_p;
    uint64_t max_p;
    uint64_t min_p;

	// threshold at which we know number are prime, as determined by trial division
	uint64_t prime_threshold;

	// manage a list of factors
	yfactor_list_t *factors;
	int do_logging;
	int refactor_depth;

    // configurable options
    int VFLAG;
    int THREADS;
    int LOGFLAG;
    int LATHREADS;
    int NUM_WITNESSES;
    int num_threads;

    // computer info
    double MEAS_CPU_FREQUENCY;
    char CPU_ID_STR[80];
    int HAS_SSE2;
    int HAS_SSE41;
    int HAS_AVX;
    int HAS_AVX2;
    int HAS_BMI2;
    int HAS_AVX512;
    int HAS_AVX512F;
    int HAS_AVX512BW;
    int HAS_AVX512VL;
    int HAS_AVX512PF;
    int HAS_AVX512IFMA52;
    int L1CACHE;
    int L2CACHE;
    int L3CACHE;
    uint32_t cache_size1;
    uint32_t cache_size2;
    
    // RNG state
    uint64_t lcg_state;
    uint32_t seed1;
    uint32_t seed2;

} fact_obj_t;

void init_factobj(fact_obj_t *fobj);
void free_factobj(fact_obj_t *fobj);
void reset_factobj(fact_obj_t *fobj);
void copy_factobj(fact_obj_t* dest, fact_obj_t* src);
void alloc_factobj(fact_obj_t *fobj);


//#if defined(WIN32)
// windows machines also need these declarations for functions located
// within common.lib or gnfs.lib, for NFS factorizations using msieve.
void savefile_init(savefile_t *s, char *filename);
void savefile_free(savefile_t *s);
void savefile_open(savefile_t *s, uint32_t flags);
void savefile_close(savefile_t *s);
uint32_t savefile_eof(savefile_t *s);
uint32_t savefile_exists(savefile_t *s);
void savefile_read_line(char *buf, size_t max_len, savefile_t *s);
void savefile_write_line(savefile_t *s, char *buf);
void savefile_flush(savefile_t *s);
void savefile_rewind(savefile_t *s);

//#endif


/* ============================ interface to tinyecm ============================ */
extern void tinyecm(mpz_t n, mpz_t f, uint32_t B1, uint32_t B2, uint32_t curves,
    uint64_t* lcg_state, int verbose);

// getfactor_tecm() returns 0 if unable to find a factor of n,
// Otherwise it returns 1 and a factor of n in argument f.
// 
// if the input is known to have no small factors, set is_arbitrary=0, 
// otherwise, set is_arbitrary=1 and a few curves targetting small factors
// will be run prior to the standard sequence of curves for the input size.
//  
// Prior to your first call of getfactor_tecm(), set *pran = 0  (or set it to
// some other arbitrary value); after that, don't change *pran.
// FYI: *pran is used within this file by a random number generator, and it
// holds the current value of a pseudo random sequence.  Your first assigment
// to *pran seeds the sequence, and after seeding it you don't want to
// change *pran, since that would restart the sequence.
int getfactor_tecm(mpz_t n, mpz_t f, int is_arbitrary, uint64_t* pran);
void getfactor_tecm_x8_list(uint64_t* n, uint64_t* f, int target_bits, uint32_t num_in, uint64_t* pran);
int getfactor_tecm_x8(mpz_t n, mpz_t f, int target_bits, uint64_t* pran);
int getfactor_tpm1(mpz_t n, mpz_t f, uint32_t b1);

/* =============== interface to various small-factor-finding routines =========== */
uint64_t spbrent(uint64_t N, uint64_t c, int imax);
uint64_t spbrent64(uint64_t N, int imax);
int mbrent(fact_obj_t *fobj);
int montybrent(monty_t *mdata, mpz_t n, mpz_t f, uint32_t a, uint32_t imax);
void brent_loop(fact_obj_t *fobj);
uint64_t sp_shanks_loop(mpz_t N, fact_obj_t *fobj);
uint64_t LehmanFactor(uint64_t N, double Tune, int DoTrial, double CutFrac);
void init_lehman();
void zTrial(fact_obj_t *fobj);
void zFermat(uint64_t limit, uint32_t mult, fact_obj_t *fobj);
void factor_perfect_power(fact_obj_t *fobj, mpz_t b);
int par_shanks_loop(uint64_t*N, uint64_t*f, int num_in);

int sptestsqr(uint64_t n);
uint64_t spfermat(uint64_t limit, uint32_t mult, uint64_t n);


// factoring related utility
int resume_check_input_match(mpz_t file_n, mpz_t input_n, mpz_t common_fact, int VFLAG);

// find a crossover between QS and NFS
void factor_tune(fact_obj_t* inobj);

#endif //_FACTOR_H
