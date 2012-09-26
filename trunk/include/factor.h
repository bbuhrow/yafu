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
#include "util.h"

//#define NO_ZLIB

#if !defined(NO_ZLIB) && !defined(__MINGW32__) 
#include <zlib.h>
#else
#define NO_ZLIB
#endif

#ifdef _MSC_VER
#define USE_NFS
#endif

/* declarations of factorization routines not grouped in with QS methods
or group theoretic methods, as well as higher level factorization
routines and bookkeeping */

#define MAX_FACTORS 10

// stuff that needs to be visible to the msieve routines and 
// the yafu sieve routines

/* structure encapsulating the savefile used in a factorization */
typedef struct {

#if defined(WIN32) || defined(_WIN64)
	HANDLE file_handle;
	uint32 read_size;
	uint32 eof;
#else
	FILE *fp;
#endif
	char *name;
	char *buf;
	uint32 buf_off;
} qs_savefile_t;

typedef struct {

#if defined(NO_ZLIB) && (defined(WIN32) || defined(_WIN64))
	HANDLE file_handle;
	uint32 read_size;
	uint32 eof;
#else
	gzFile *fp;
	char isCompressed;
	char is_a_FILE;
#endif
	char *name;
	char *buf;
	uint32 buf_off;
} savefile_t;

enum msieve_flags {
	MSIEVE_DEFAULT_FLAGS = 0,		/* just a placeholder */
	MSIEVE_FLAG_USE_LOGFILE = 0x01,	    /* append log info to a logfile */
	MSIEVE_FLAG_LOG_TO_STDOUT = 0x02,   /* print log info to the screen */
	MSIEVE_FLAG_STOP_SIEVING = 0x04,    /* tell library to stop sieving
					       when it is safe to do so */
	MSIEVE_FLAG_FACTORIZATION_DONE = 0x08,  /* set by the library if a
						   factorization completed */
	MSIEVE_FLAG_SIEVING_IN_PROGRESS = 0x10, /* set by the library when 
						   any sieving operations are 
						   in progress */
	MSIEVE_FLAG_SKIP_QS_CYCLES = 0x20,  /* do not perform exact tracking of
	                                    the number of cycles while sieving
					    is in progress; for distributed
					    sieving where exact progress info
					    is not needed, sieving clients can
					    save a lot of memory with this */
	MSIEVE_FLAG_NFS_POLY1 = 0x40,     /* if input is large enough, perform
	                                    stage 1 polynomial selection for NFS */
	MSIEVE_FLAG_NFS_POLY2 = 0x80,     /* if input is large enough, perform
	                                    stage 2 polynomial selection for NFS */
	MSIEVE_FLAG_NFS_SIEVE = 0x100,   /* if input is large enough, perform
	                                    sieving for NFS */
	MSIEVE_FLAG_NFS_FILTER = 0x200,  /* if input is large enough, perform
	                                    filtering phase for NFS */
	MSIEVE_FLAG_NFS_LA = 0x400,      /* if input is large enough, perform
	                                    linear algebra phase for NFS */
	MSIEVE_FLAG_NFS_SQRT = 0x800,    /* if input is large enough, perform
	                                    square root phase for NFS */
	MSIEVE_FLAG_NFS_LA_RESTART = 0x1000,/* restart the NFS linear algbra */
	MSIEVE_FLAG_DEEP_ECM = 0x2000    /* perform nontrivial-size ECM */
};
	

enum msieve_factor_type {
	MSIEVE_COMPOSITE,
	MSIEVE_PRIME,
	MSIEVE_PROBABLE_PRIME
};

typedef struct msieve_factor {
	enum msieve_factor_type factor_type;
	char *number;
	struct msieve_factor *next;
} msieve_factor;


//* One factorization is represented by a msieve_obj
//   structure. This contains all the static information
//   that gets passed from one stage of the factorization
//   to another. If this was C++ it would be a simple object */
// this must be declared when calling msieve routines through
// the msieve libraries
typedef struct {
	char *input;		  /* pointer to string version of the 
				     integer to be factored */
	msieve_factor *factors;   /* linked list of factors found (in
				     ascending order */
	volatile uint32 flags;	  /* input/output flags */
	savefile_t savefile;      /* data for savefile */
	char *logfile_name;       /* name of the logfile that will be
				     used for this factorization */
	uint32 seed1, seed2;      /* current state of random number generator
				     (updated as random numbers are created) */
	char *nfs_fbfile_name;    /* name of factor base file */
	uint32 max_relations;      /* the number of relations that the sieving
	                              stage will try to find. The default (0)
				      is to keep sieving until all necessary 
				      relations are found. */
	uint64 nfs_lower;         /* lower bound for NFS-related thing to do */
	uint64 nfs_upper;         /* upper bound for NFS-related thing to do */

	uint32 cache_size1;       /* bytes in level 1 cache */
	uint32 cache_size2;       /* bytes in level 2 cache */
	enum cpu_type cpu;
	uint32 num_threads;

#ifdef HAVE_MPI
	uint32 mpi_size;          /* number of MPI processes, each with
                                     num_threads threads */
	uint32 mpi_rank;          /* from 0 to mpi_size - 1 */

	uint32 mpi_nrows;         /* a 2-D MPI lanczos grid */
	uint32 mpi_ncols;
	MPI_Comm mpi_la_grid;
	MPI_Comm mpi_la_row_grid; /* communicator for the current MPI row */
	MPI_Comm mpi_la_col_grid; /* communicator for the current MPI col */
	uint32 mpi_la_row_rank;
	uint32 mpi_la_col_rank;
#endif


	uint32 mem_mb;            /* megabytes usable for NFS filtering */

	uint32 which_gpu;         /* ordinal ID of GPU to use */

	double target_density;   /* the target number of nonzeros per matrix
				    column, used in NFS filtering */
	char *mp_sprintf_buf;    /* scratch space for printing big integers */
} msieve_obj;

// these must be declared when calling msieve routines through
// the msieve libraries.  they are defined in the msieve library.
msieve_obj * msieve_obj_new(char *input_integer,
			    uint32 flags,
			    char *savefile_name,
			    char *logfile_name,
			    char *nfs_fbfile_name,
			    uint32 seed1,
			    uint32 seed2,
			    uint32 max_relations,
			    uint64 nfs_lower,
			    uint64 nfs_upper,
			    enum cpu_type cpu,
			    uint32 cache_size1,
			    uint32 cache_size2,
			    uint32 num_threads,
			    uint32 mem_mb,
			    uint32 which_gpu,
			    double target_density);

msieve_obj * msieve_obj_free(msieve_obj *obj);


typedef struct {
	mp_t factor;
	enum msieve_factor_type type;
} final_factor_t;

typedef struct {
	uint32 num_factors;
	final_factor_t *final_factors[256];
} factor_list_t;


/*--------------LINEAR ALGEBRA RELATED DECLARATIONS ---------------------*/

/* Used to represent a list of relations */

typedef struct {
	uint32 num_relations;  /* number of relations in the cycle */
	uint32 *list;          /* list of offsets into an array of relations */
} qs_la_cycle_t;

/* A column of the matrix */

typedef struct {
	uint32 *data;		/* The list of occupied rows in this column */
	uint32 weight;		/* Number of nonzero entries in this column */
	qs_la_cycle_t cycle;       /* list of relations comprising this column */
} qs_la_col_t;

/*---------------- SAVEFILE RELATED DECLARATIONS ---------------------*/

#define LINE_BUF_SIZE 300
#define SAVEFILE_READ 0x01
#define SAVEFILE_WRITE 0x02
#define SAVEFILE_APPEND 0x04

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

} pp1_obj_t;

typedef struct
{
	mpz_t gmp_n;
	mpz_t gmp_f;

	char ecm_path[1024];
	int use_external;
	uint32 B1;
	uint64 B2;
	int stg2_is_default;
	int curves_run;
	uint32 num_curves;
	uint32 sigma;				//sigma value of successful curve
	uint32 num_factors;			//number of factors found in this method
	z *factors;					//array of bigint factors found in this method
	double ttime;

	// fit parameters to compute time_per_curve as a function of B1
	double ecm_exponent;
	double ecm_multiplier;
	double ecm_tune_freq;		

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
	z *factors;					//array of bigint factors found in this method
	double ttime;
	int print;

} div_obj_t;

typedef struct
{
	mpz_t gmp_n;
	mpz_t gmp_f;
	uint32 num_factors;			//number of factors found in this method
	z *factors;					//array of bigint factors found in this method
	double ttime;

} squfof_obj_t;

typedef struct
{
	mpz_t gmp_n;
	qs_savefile_t savefile;		//savefile object
	char siqs_savefile[1024];

	double qs_exponent;
	double qs_multiplier;
	double qs_tune_freq;
	int no_small_cutoff_opt;	//1 is true - perform no optimization.  0 is false.

	int gbl_override_B_flag;
	uint32 gbl_override_B;			//override the # of factor base primes
	int gbl_override_tf_flag;
	uint32 gbl_override_tf;			//extra reduction of the TF bound by X bits
	int gbl_override_time_flag;
	uint32 gbl_override_time;		//stop after this many seconds
	int gbl_override_rel_flag;
	uint32 gbl_override_rel;		//stop after collecting this many relations
	int gbl_override_blocks_flag;
	uint32 gbl_override_blocks;		//override the # of blocks used
	int gbl_override_lpmult_flag;
	uint32 gbl_override_lpmult;		//override the large prime multiplier
	int gbl_force_DLP;

	uint32 num_factors;			//number of factors found in this method
	z *factors;					//array of bigint factors found in this method
	uint32 flags;				//each bit corresponds to a location in the 
								//flags enum
	double rels_per_sec;
	double qs_time;
	double total_time;

} qs_obj_t;

typedef struct
{
	mpz_t gmp_n;
	char ggnfs_dir[1024];
	uint32 siever;
	uint32 startq;
	uint32 rangeq;
	uint32 polystart;
	uint32 polyrange;
	char outputfile[1024];
	char logfile[1024];
	char fbfile[1024];
	int sq_side;
	uint32 timeout;
	char job_infile[1024];
	int sieve_only;
	int poly_only;
	int post_only;
	int poly_option;
	int restart_flag;
	int la_restart;
	int user_job;
	uint32 polybatch;

	double gnfs_exponent;
	double gnfs_multiplier;
	double gnfs_tune_freq;

	uint32 min_digits;

	uint32 num_factors;			//number of factors found in this method
	z *factors;					//array of bigint factors found in this method
	double ttime;

	msieve_obj *mobj;

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
	int prefer_xover;

	//balance of ecm and various sieve methods
	double target_pretest_ratio;
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
	z N;						//numerical representation of input
	str_t str_N;				//string representation of input (expression form)
	uint32 digits;				//number of digits in input
	uint32 bits;				//number of bits in input
	char flogname[1024];		//name of the factorization logfile to use
	FILE *logfile;				//the logfile
	char savefile_name[80];		//data savefile name
	uint32 flags;				//state flags

	//info for work done in various places during this factorization
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

	//global storage for a list of factors
	factor_t *fobj_factors;
	uint32 num_factors;
	uint32 allocated_factors;

} fact_obj_t;

void init_factobj(fact_obj_t *fobj);
void free_factobj(fact_obj_t *fobj);
void reset_factobj(fact_obj_t *fobj);
void alloc_factobj(fact_obj_t *fobj);

/* ---------------- DECLARATIONS FOR SAVEFILE MANIPULATION -------------- */
// copied from msieve source code
void qs_savefile_init(qs_savefile_t *s, char *filename);
void qs_savefile_free(qs_savefile_t *s);
void qs_savefile_open(qs_savefile_t *s, uint32 flags);
void qs_savefile_close(qs_savefile_t *s);
uint32 qs_savefile_eof(qs_savefile_t *s);
uint32 qs_savefile_exists(qs_savefile_t *s);
void qs_savefile_rewind(qs_savefile_t *s);
void qs_savefile_read_line(char *buf, size_t max_len, qs_savefile_t *s);
void qs_savefile_write_line(qs_savefile_t *s, char *buf);
void qs_savefile_flush(qs_savefile_t *s);

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


/*--------------DECLARATIONS FOR MANAGING FACTORS FOUND -----------------*/

//msieve
void factor_list_init(factor_list_t *list);
uint32 factor_list_max_composite(factor_list_t *list);
void factor_list_free(z *n, factor_list_t *list, fact_obj_t *obj);

//yafu
void add_to_factor_list(fact_obj_t *fobj, mpz_t n);
void print_factors(fact_obj_t *fobj);
void free_factor_list(fact_obj_t *fobj);
void clear_factor_list(fact_obj_t *fobj);
void delete_from_factor_list(fact_obj_t *fobj, mpz_t n);

/*-----------TOP LEVEL ENTRY POINT FOR ALL FACTORING ROUTINES ----------*/

uint32 test_qn_res[128];

void brent_loop(fact_obj_t *fobj);
void pollard_loop(fact_obj_t *fobj);
void williams_loop(fact_obj_t *fobj);
int ecm_loop(fact_obj_t *fobj);
uint64 sp_shanks_loop(mpz_t N, fact_obj_t *fobj);
uint64 LehmanFactor(uint64 N, double Tune, int DoTrial, double CutFrac);
void init_lehman();
void zTrial(fact_obj_t *fobj);
void zFermat(uint64 limit, uint32 mult, fact_obj_t *fobj);
void factor_perfect_power(fact_obj_t *fobj, mpz_t b);
void nfs(fact_obj_t *fobj);
void SIQS(fact_obj_t *fobj);
void smallmpqs(fact_obj_t *fobj);
//void tinySIQS(fact_obj_t *fobj);
void tinySIQS(mpz_t n, mpz_t *factors, uint32 *num_factors);

int sptestsqr(uint64 n);
uint64 spfermat(uint64 n, uint64 limit);
void spfactorlist(uint64 nstart, uint64 nrange);

//auto factor routine
void factor(fact_obj_t *fobj);

// factoring related utility
int resume_check_input_match(mpz_t file_n, mpz_t input_n, mpz_t common_fact);

/* Factor a number using GNFS. Returns
   1 if any factors were found and 0 if not */

uint32 factor_gnfs(msieve_obj *obj, mp_t *n, factor_list_t *factor_list);

void factor_tune(fact_obj_t *fobj);

#endif //_FACTOR_H
