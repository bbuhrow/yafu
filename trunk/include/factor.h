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

/*
declarations of factorization routines not grouped in with QS methods
or group theoretic methods, as well as higher level factorization
routines and bookkeeping
*/

#define MAX_FACTORS 10
#define POLLARD_METHOD 0
#define WILLIAMS_METHOD 1
#define STG1 1
#define STG2 2

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
	MSIEVE_FLAG_NFS_POLY = 0x40,     /* if input is large enough, perform
	                                    polynomial selection for NFS */
	MSIEVE_FLAG_NFS_SIEVE = 0x80,    /* if input is large enough, perform
	                                    sieving for NFS */
	MSIEVE_FLAG_NFS_FILTER = 0x100,  /* if input is large enough, perform
	                                    filtering phase for NFS */
	MSIEVE_FLAG_NFS_LA = 0x200,      /* if input is large enough, perform
	                                    linear algebra phase for NFS */
	MSIEVE_FLAG_NFS_SQRT = 0x400,    /* if input is large enough, perform
	                                    square root phase for NFS */
	MSIEVE_FLAG_NFS_LA_RESTART = 0x800,/* restart the NFS linear algbra */
	MSIEVE_FLAG_DEEP_ECM = 0x1000    /* perform nontrivial-size ECM */
};

enum msieve_factor_type {
	MSIEVE_COMPOSITE,
	MSIEVE_PRIME,
	MSIEVE_PROBABLE_PRIME
};

typedef struct {
	z factor;
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
} la_cycle_t;

/* A column of the matrix */

typedef struct {
	uint32 *data;		/* The list of occupied rows in this column */
	uint32 weight;		/* Number of nonzero entries in this column */
	la_cycle_t cycle;       /* list of relations comprising this column */
} la_col_t;

/*---------------- SAVEFILE RELATED DECLARATIONS ---------------------*/

#define LINE_BUF_SIZE 300
#define SAVEFILE_READ 0x01
#define SAVEFILE_WRITE 0x02
#define SAVEFILE_APPEND 0x04

// factorization objects //

typedef struct
{
	z n;
	str_t in, out;
	uint32 B1;
	uint64 B2;
	uint32 num_factors;			//number of factors found in this method
	z *factors;					//array of bigint factors found in this method

} pm1_obj_t;

typedef struct
{
	z n;
	str_t in, out;
	uint32 B1;
	uint64 B2;
	uint32 *bases;
	uint32 num_factors;			//number of factors found in this method
	z *factors;					//array of bigint factors found in this method

} pp1_obj_t;

typedef struct
{
	z n;
	str_t in, out;
	uint32 B1;
	uint64 B2;
	uint32 num_curves;
	uint32 sigma;				//sigma value of successful curve
	uint32 num_factors;			//number of factors found in this method
	z *factors;					//array of bigint factors found in this method

} ecm_obj_t;

typedef struct
{
	z n;
	str_t in, out;
	uint32 iterations;
	uint32 num_poly;
	uint32 *polynomials;
	uint32 num_factors;			//number of factors found in this method
	z *factors;					//array of bigint factors found in this method

} rho_obj_t;

typedef struct
{
	z n;
	str_t in, out;
	uint32 limit;
	uint32 num_factors;			//number of factors found in this method
	z *factors;					//array of bigint factors found in this method

} div_obj_t;

typedef struct
{
	z n;
	str_t in, out;
	uint32 num_factors;			//number of factors found in this method
	z *factors;					//array of bigint factors found in this method

} squfof_obj_t;

typedef struct
{
	z n;
	str_t in, out;
	savefile_t savefile;		//savefile object
	uint32 override1;			//override of default parameters
	uint32 override2;			//override of default parameters
	uint32 override3;			//override of default parameters
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
	z N;						//numerical representation of input
	str_t str_N;				//string representation of input (expression form)
	uint32 digits;				//number of digits in input
	uint32 bits;				//number of bits in input
	char logname[80];			//name of the factorization logfile to use
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

	uint32 cache_size1;
	uint32 cache_size2;
	int num_threads;
	uint32 seed1;
	uint32 seed2;

	//global storage for a list of factors
	factor_t *fobj_factors;
	uint32 num_factors;
	uint32 allocated_factors;

} fact_obj_t;

void init_factobj(fact_obj_t *fobj);
void record_new_factor(fact_obj_t *fobj,char *method, z *n);
void free_factobj(fact_obj_t *fobj);

void savefile_init(savefile_t *s, char *filename);
void savefile_free(savefile_t *s);
void savefile_open(savefile_t *s, uint32 flags);
void savefile_close(savefile_t *s);
uint32 savefile_eof(savefile_t *s);
uint32 savefile_exists(savefile_t *s);
void savefile_rewind(savefile_t *s);
void savefile_read_line(char *buf, size_t max_len, savefile_t *s);
void savefile_write_line(savefile_t *s, char *buf);
void savefile_flush(savefile_t *s);


/*--------------DECLARATIONS FOR MANAGING FACTORS FOUND -----------------*/

void factor_list_init(factor_list_t *list);

uint32 factor_list_max_composite(factor_list_t *list);

void factor_list_free(z *n, factor_list_t *list, fact_obj_t *obj);

//common to the QS variants
uint8 choose_multiplier(z *n, uint32 fbsize);

// factoring routines not declared inside their own header //

//shanks SQUFOF
void shanks_mult(uint64 N, uint64 *f);
uint64 sp_shanks_loop(z *N, fact_obj_t *fobj);
uint32 squfof(z *n);
//jasonp's squfof, for comparison
uint32 squfof_jp(z *n);
int SQUFOF_alpertron(int64 N, int64 queue[]);
void enqu(int q,int *iter);
int squfof_rds(int64 n, int *fac1, int *fac2);

//trial division
void zTrial(fp_digit limit, int print, fact_obj_t *fobj);
void zFermat(fp_digit limit, fact_obj_t *fobj);
void Trial64(int64 n, int print);
void Trial32(int32 n, int print);

//auto factoring routines
double get_qs_time_estimate(double freq, int bits);
void factor(fact_obj_t *fobj);
void factor_test(void);
void factor_big(z *b);
void factor_range_pp1(int kmin, int kmax, int n, int sign);
void factor_range_pm1(int kmin, int kmax, int n, int sign);
void add_to_factor_list(fact_obj_t *fobj, z *n);
void print_factors(fact_obj_t *fobj);
void free_factor_list(fact_obj_t *fobj);
void clear_factor_list(fact_obj_t *fobj);
void delete_from_factor_list(fact_obj_t *fobj, z *n);
void aliquot(z *input, fact_obj_t *fobj);

//algebraic factoring
int algebraic(str_t *str);

#endif //_FACTOR_H
