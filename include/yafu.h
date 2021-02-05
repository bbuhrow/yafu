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
#ifndef _YAFU_HEAD_DEF
#define _YAFU_HEAD_DEF

//#define DEBUG

#define _CRT_SECURE_NO_WARNINGS 

#define VERSION_STRING "2.0"

//max words for fixed precision msieve bignum
#define MAX_MP_WORDS 64

//default maximum size in chars for a str_t
#define GSTR_MAXSIZE 1024

//#define FORCE_GENERIC 1

//support libraries
#include "types.h"
#include <gmp.h>

//global typedefs
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
	mpz_t factor;
	int count;
	int type;
    // new parameters to support returning factors
    // from avx-ecm.  In general I think it makes sense
    // to have this structure contain information not
    // only about the factor itself but how it was found.
    int method;     // factorization method used
    uint32 curve_num;   // curve found
    uint64 sigma;   // sigma value
    int tid;        // thread found
    int vid;        // vector position found
} factor_t;

// a structure to hold a bunch of configuration info
// for yafu, instead of declaring a bunch of globals.
typedef struct
{
    // isprime
    uint32 NUM_WITNESSES;

    // output behavior - used everywhere
    int VFLAG, LOGFLAG;

    // threading - used many places (factoring, SoE)
    int THREADS;
    int LATHREADS;

    // input options
    int USEBATCHFILE;
    int USERSEED;
    int CMD_LINE_REPEAT;
    char batchfilename[1024];
    char sessionname[1024];
    char scriptname[1024];
    int NO_CLK_TEST;

    // machine info
    double MEAS_CPU_FREQUENCY;
    int VERBOSE_PROC_INFO;
    char CPU_ID_STR[80];
    uint32 L1CACHE, L2CACHE;
    int CLSIZE;
    char HAS_SSE41;
    char HAS_AVX;
    char HAS_AVX2;
#if defined(WIN32)
    char sysname[MAX_COMPUTERNAME_LENGTH + 1];
    int sysname_sz;
#else
    char sysname[256];
    int sysname_sz;
#endif

    

} yafu_obj_t;

//this array holds a global store of prime numbers
uint32* spSOEprimes;	//the primes	
uint32 szSOEp;			//count of primes

//this array holds NUM_P primes in the range P_MIN to P_MAX, and
//can change as needed - always check the range and size to see
//if the primes you need are in there before using it
uint64* PRIMES;
uint64 NUM_P;
uint64 P_MIN;
uint64 P_MAX;

extern void yafu_init(yafu_obj_t* yobj);
extern void yafu_finalize(yafu_obj_t* yobj);

#endif //ifndef HEAD_DEF
