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

//definitions
#ifndef HEAD_DEF
#define HEAD_DEF

//#define DEBUG

#define _CRT_SECURE_NO_WARNINGS 

#define VERSION_STRING "1.26.5"

//basics
#define POSITIVE 0
#define NEGATIVE 1

//bases
#define DEC 10
#define HEX 16
#define OCT 8
#define BIN 2

//max words for fixed precision msieve bignum
#define MAX_MP_WORDS 32

//default maximum size in chars for a str_t
#define GSTR_MAXSIZE 1024

//support libraries
#include "types.h"

//global typedefs
typedef struct
{
	fp_digit *val;
	int alloc;
	int size;
	int type;
} z;

typedef struct
{
	uint32 *val;
	int alloc;
	int size;
	int type;
} z32;

typedef struct {
	uint32 nwords;		/* number of nonzero words in val[] */
	uint32 val[MAX_MP_WORDS];
} mp_t;

typedef struct {
	uint32 sign;	/* POSITIVE or NEGATIVE */
	mp_t num;
} signed_mp_t;

typedef struct
{
	char *s;		//pointer to beginning of s
	int nchars;		//number of valid characters in s (including \0)
	int alloc;		//bytes allocated to s
} str_t;

typedef struct
{
	char name[40];
	z data;
} uvar_t;

typedef struct
{
	uvar_t *vars;
	int num;
	int alloc;
} uvars_t;

typedef struct
{
	z r;
	z rhat;
	z nhat;
	z one;
} monty;

typedef struct
{
	z factor;
	int count;
} factor_t;

typedef struct
{
	str_t **elements;	//an array of pointers to elements
	int num;			//number of elements
	int size;			//allocated number of elements in stack
	int top;			//the top element
	int type;			//is this a stack (0) or a queue (1)?
} bstack_t;

typedef struct
{
	uint32 hi;
	uint32 low;
} rand_t;

typedef struct
{
	uint32 stg1;
	uint64 stg2;
	uint32 base;
	z *factors;
	uint32 *exp;
	double stg1time;
	double stg2time;
} pm1_data_t;

typedef struct
{
	uint32 stg1;
	uint64 stg2;
	uint32 numbases;
	uint32 base;
	z *factors;
	uint32 *exp;
	double *stg1time;
	double *stg2time;
} pp1_data_t;

typedef struct
{
	uint32 stg1;
	uint64 stg2;
	uint32 numcurves;
	uint32 sigma;
	z *factors;
	uint32 *exp;
	double *stg1time;
	double *stg2time;
} ecm_data_t;

typedef struct
{
	uint32 iterations;
	uint32 numpoly;
	z *factors;
	uint32 *exp;
	double time;
} rho_data_t;

typedef struct
{
	uint32 limit;
	z *factors;
	uint32 *exp;
	double time;
} trial_data_t;

//user dimis:
//http://cboard.cprogramming.com/cplusplus-programming/
//101085-how-measure-time-multi-core-machines-pthreads.html
//
typedef struct {
	long		secs;
	long		usecs;
} TIME_DIFF;

//global variables
int IBASE;
int OBASE;
uint32 BRENT_MAX_IT;
uint32 FMTMAX;
uint32 WILL_STG1_MAX;
uint64 WILL_STG2_MAX;
uint32 ECM_STG1_MAX;
uint64 ECM_STG2_MAX;
uint32 POLLARD_STG1_MAX;
uint64 POLLARD_STG2_MAX;
int ECM_STG2_ISDEFAULT;
int PP1_STG2_ISDEFAULT;
int PM1_STG2_ISDEFAULT;
int VFLAG, LOGFLAG;
uint32 LOWER_POLYPOOL_INDEX;
uint32 UPPER_POLYPOOL_INDEX;
uint32 QS_DUMP_CUTOFF;
uint32 MIN_NFS_DIGITS;

// used to determine an estimate of QS runtime
double QS_EXPONENT;
double QS_MULTIPLIER;
double QS_TUNE_FREQ;

//used to determine an estimate of GNFS runtime
double GNFS_EXPONENT;
double GNFS_MULTIPLIER;
double GNFS_TUNE_FREQ;

//crossover between qs and gnfs
double QS_GNFS_XOVER;

//balance of ecm and various sieve methods
double TARGET_ECM_QS_RATIO;
double TARGET_ECM_GNFS_RATIO;
double TARGET_ECM_SNFS_RATIO;

uint32 NUM_WITNESSES;
uint32 SIGMA;
int PRIMES_TO_FILE;
int PRIMES_TO_SCREEN;
int THREADS;
int AUTO_FACTOR;
int SCALE;
int NO_ECM;
int USEBATCHFILE;
int USERSEED;
uint32 L1CACHE, L2CACHE;
int CLSIZE;
uint64 PRIME_THRESHOLD;
uint64 gcounts[10];
uint32 maxbn;
int NO_SIQS_OPT;
double MEAS_CPU_FREQUENCY;
int VERBOSE_PROC_INFO;
char CPU_ID_STR[80];
int WANT_ONLY_1_FACTOR;
int WANT_OUTPUT_PRIMES;
int WANT_OUTPUT_FACTORS;
int WANT_OUTPUT_UNFACTORED;
int WANT_OUTPUT_EXPRESSIONS;
FILE *op_file;
FILE *of_file;
FILE *ou_file;
char op_str[1024];
char of_str[1024];
char ou_str[1024];

//NFS options
uint32 GGNFS_STARTQ;
uint32 GGNFS_RANGEQ;
uint32 GGNFS_POLYSTART;
uint32 GGNFS_POLYRANGE;
char GGNFS_OUTPUTFILE[1024];
char GGNFS_LOGFILE[1024];
char GGNFS_FBFILE[1024];
int GGNFS_SQ_SIDE;
uint32 GGNFS_TIMEOUT;
char GGNFS_JOB_INFILE[1024];
int GGNFS_SIEVE_ONLY;
int GGNFS_POLY_ONLY;
int GGNFS_POST_ONLY;
int GGNFS_POLY_OPTION;
int GGNFS_RESTART_FLAG;
uint32 GGNFS_POLYBATCH;


//globals for implementing the "plan" and "pretest" switches
enum pretest_plan {
	PRETEST_NONE = 0,
	PRETEST_NOECM = 1,
	PRETEST_LIGHT = 2,
	PRETEST_NORMAL = 3,
	PRETEST_DEEP = 4,
	PRETEST_CUSTOM = 5
};
enum pretest_plan yafu_pretest_plan;
char plan_str[1024];
int ONLY_PRETEST;


//globals for testing siqs
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

//global workspace variables
z zZero, zOne, zTwo, zThree;

//this array holds a global store of prime numbers
uint64 *spSOEprimes;	//the primes	
uint64 szSOEp;			//count of primes

//this array holds NUM_P primes in the range P_MIN to P_MAX, and
//can change as needed - always check the range and size to see
//if the primes you need are in there before using it
uint64 *PRIMES;
uint64 NUM_P;
uint64 P_MIN;
uint64 P_MAX;

//a few strings - mostly for logprint stuff
str_t gstr1, gstr2, gstr3;

//factorization log file
char flogname[1024];

//input batch file
char batchfilename[1024];

//siqs save file name
char siqs_savefile[1024];

//session name
char sessionname[1024];

//ggnfs binaries directory
char ggnfs_dir[1024];

//random seeds and the state of the LCG
rand_t g_rand;
uint64 LCGSTATE;

#if defined(WIN32)
//system info
char sysname[MAX_COMPUTERNAME_LENGTH + 1];
int sysname_sz;

#else
//system info
char sysname[80];
int sysname_sz;

#endif

//option flag processing
void applyOpt(char *opt, char *arg);
unsigned process_flags(int argc, char **argv);

#endif //ifndef HEAD_DEF



