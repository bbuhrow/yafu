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

#define VERSION_STRING "1.35-beta"

//basics
#define POSITIVE 0
#define NEGATIVE 1

//bases
#define DEC 10
#define HEX 16
#define OCT 8
#define BIN 2

//max words for fixed precision msieve bignum
#define MAX_MP_WORDS 64

//default maximum size in chars for a str_t
#define GSTR_MAXSIZE 1024

//#define FORCE_GENERIC 1

//support libraries
#include "types.h"
#include <gmp.h>

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
	mpz_t data;
} uvar_t;

typedef struct
{
    char name[40];
    char *data;
} strvar_t;

typedef struct
{
	uvar_t *vars;
	int num;
	int alloc;
} uvars_t;

typedef struct
{
    strvar_t *vars;
    int num;
    int alloc;
} strvars_t;

typedef struct
{
	mpz_t factor;
	int count;
	int type;
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

//global variables

// calculator
int IBASE;
int OBASE;

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

// random numbers - used everywhere
rand_t g_rand;
gmp_randstate_t gmp_randstate;
uint64 LCGSTATE;

//SoE
int PRIMES_TO_FILE;
int PRIMES_TO_SCREEN;

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

//global workspace variables
z zZero, zOne, zTwo, zThree, zFive;

//this array holds a global store of prime numbers
uint32 *spSOEprimes;	//the primes	
uint32 szSOEp;			//count of primes

//this array holds NUM_P primes in the range P_MIN to P_MAX, and
//can change as needed - always check the range and size to see
//if the primes you need are in there before using it
uint64 *PRIMES;
uint64 NUM_P;
uint64 P_MIN;
uint64 P_MAX;

//a few strings - mostly for logprint stuff
str_t gstr1, gstr2, gstr3;

#endif //ifndef HEAD_DEF
