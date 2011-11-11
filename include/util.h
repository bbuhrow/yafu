/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	
       				   --jasonp@boo.net 9/24/08

Modified:	Ben Buhrow
Date:		11/24/09
Purpose:	Port into Yafu-1.14.
--------------------------------------------------------------------*/

#ifndef _UTIL_H_
#define _UTIL_H_

#include "yafu.h"

/* system-specific stuff ---------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

/* basic types  -------------------------------------------------------*/

	//moved to types.h

/* useful functions ---------------------------------------------------*/

#define INLINE __inline
#if defined(_MSC_VER)
	#define getpid _getpid
#endif

#if defined(__GNUC__) && __GNUC__ >= 3
	#define PREFETCH(addr) __builtin_prefetch(addr) 
#elif defined(_MSC_VER) && _MSC_VER >= 1400
	#define PREFETCH(addr) PreFetchCacheLine(PF_TEMPORAL_LEVEL_1, addr)
#else
	#define PREFETCH(addr) /* nothing */
#endif

#define MIN(a,b) ((a) < (b)? (a) : (b))
#define MAX(a,b) ((a) > (b)? (a) : (b))

static INLINE void * xmalloc_align(size_t len)
{
#if defined (_MSC_VER) || defined(__MINGW32__)
	void *ptr = _aligned_malloc(len, 64);

#elif defined (__APPLE__)
	void *ptr = malloc(len);

#elif defined (__GNUC__)
	void *ptr = memalign(64, len);

#else
	void *ptr = malloc(len);

#endif

	return ptr;
}

static INLINE void * xmalloc(size_t len) {
	void *ptr = malloc(len);
	if (ptr == NULL) {
		printf("failed to allocate %u bytes\n", (uint32)len);
		exit(-1);
	}
	return ptr;
}

static INLINE void * xcalloc(size_t num, size_t len) {
	void *ptr = calloc(num, len);
	if (ptr == NULL) {
		printf("failed to calloc %u bytes\n", (uint32)(num * len));
		exit(-1);
	}
	return ptr;
}

static INLINE void * xrealloc(void *iptr, size_t len) {
	void *ptr = realloc(iptr, len);
	if (ptr == NULL) {
		printf("failed to reallocate %u bytes\n", (uint32)len);
		exit(-1);
	}
	return ptr;
}

void get_random_seeds(rand_t *r);

static INLINE uint32 
get_rand(uint32 *rand_seed, uint32 *rand_carry) {
   
	/* A multiply-with-carry generator by George Marsaglia.
	   The period is about 2^63. */

	#define RAND_MULT 2131995753

	uint64 temp;

	temp = (uint64)(*rand_seed) * 
		       (uint64)RAND_MULT + 
		       (uint64)(*rand_carry);
	*rand_seed = (uint32)temp;
	*rand_carry = (uint32)(temp >> 32);
	return (uint32)temp;
}


//user dimis:
//http://cboard.cprogramming.com/cplusplus-programming/
//101085-how-measure-time-multi-core-machines-pthreads.html
//
TIME_DIFF * my_difftime (struct timeval *, struct timeval *);

//http://www.openasthra.com/c-tidbits/gettimeofday-function-for-windows/
#if defined (_MSC_VER)
	int gettimeofday(struct timeval *tv, struct timezone *tz);
#endif

//routines used all over
void logprint(FILE *infile, char *args, ...);
void logprint_oc(const char *name, const char *method, char *args, ...);
char *gettimever(char *s);
void dbl2z(double n, z *a);
//char *strrev(char *str);
fp_digit spRand(fp_digit lower, fp_digit upper);
void zRand(z *n, uint32 ndigits);
void zRandb(z *n, int bits);
void build_RSA(int bits, z *n);
void gordon(int bits, z *p);
void zNextPrime(z *n, z *p, int dir);
void zNextPrime_1(fp_digit n, fp_digit *p, z *work, int dir);
void helpfunc(char *s);
int qcomp_int(const void *x, const void *y);
void * aligned_malloc(size_t len, uint32 align);
void aligned_free(void *newptr);
uint64 measure_processor_speed(void);
int lock_thread_to_core(void);
int unlock_thread_from_core(void);
void set_idle_priority(void);
int qcomp_int(const void *x, const void *y);
int qcomp_uint16(const void *x, const void *y);
int qcomp_uint32(const void *x, const void *y);
int qcomp_uint64(const void *x, const void *y);
void generate_pseudoprime_list(int num, int bits);

//routines for testing various aspects of code
void test_dlp_composites(void);
void modtest(int it);
void test_qsort(void);
void arith_timing(int num);

/* for turning on CPU-specific code */

//enum cpu_type {
//	cpu_generic,
//	cpu_pentium,
//	cpu_pentium2,
//	cpu_pentium3,
//	cpu_pentium4,
//	cpu_pentium_m,
//	cpu_core,
//	cpu_athlon,
//	cpu_athlon_xp,
//	cpu_opteron
//};

enum cpu_type {
	cpu_generic,
	cpu_pentium,
	cpu_pentium2,
	cpu_pentium3,
	cpu_pentium4,
	cpu_pentium_m,
	cpu_core,
	cpu_athlon,
	cpu_athlon_xp,
	cpu_opteron,
	//cpu_nehalem
};

uint64 yafu_read_clock(void);
void yafu_get_cache_sizes(uint32 *level1_cache, uint32 *level2_cache);
enum cpu_type yafu_get_cpu_type(void);
//enum cpu_type get_cpu_type(void);

int extended_cpuid(char *idstr, int *cachelinesize, int do_print);

/* CPU-specific capabilities */

/* assume for all CPUs, even non-x86 CPUs. These guard
   assembly language that has other guards anyway, and
   the only CPU that doesn't have these instructions is
   the classic Pentium */

#define HAS_CMOV
#define FORCE_MODERN
//#define HAS_SSE2
//#define CACHE_LINE_64

#if defined(CPU_GENERIC)
	#define MANUAL_PREFETCH
	#if !defined(WIN32) && !defined(__i386__)
		#define HAS_MANY_REGISTERS
	#endif
	#define CACHE_LINE_32

#elif defined(CPU_PENTIUM2) 
	#define MANUAL_PREFETCH
	#define CACHE_LINE_32

#elif defined(CPU_ATHLON)
	#define MANUAL_PREFETCH
	#define HAS_AMD_MMX
	#define CACHE_LINE_32

#elif defined(CPU_PENTIUM3) 
	#define MANUAL_PREFETCH
	#define HAS_SSE
	#define CACHE_LINE_32

#elif defined(CPU_ATHLON_XP)
	#define HAS_SSE

#elif defined(CPU_PENTIUM4) || defined(CPU_PENTIUM_M) || \
	defined(CPU_CORE) || defined(CPU_OPTERON)
	#define HAS_SSE
	#define HAS_SSE2
	#if !defined(WIN32) && !defined(__i386__)
		#define HAS_MANY_REGISTERS
	#endif
	#define CACHE_LINE_64

#elif defined(FORCE_MODERN)
	#define HAS_MMX
	#define HAS_SSE
	#define HAS_SSE2
	#define CACHE_LINE_64
#endif


#ifdef __cplusplus
}
#endif

#endif /* _UTIL_H_ */
