/*--------------------------------------------------------------------
MIT License

Copyright (c) 2021 bbuhrow

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
--------------------------------------------------------------------*/

#ifndef _YTOOLS_UTIL_H_
#define _YTOOLS_UTIL_H_


#ifdef __cplusplus
extern "C" {
#endif


#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stddef.h>
#include <time.h>
#include <sys/types.h>
#include <stdarg.h>     // va_start, va_end, va_list, ...
#include <errno.h>      // strerror, errno ...

#if defined(WIN32) || defined(_WIN64) 
#define WIN32_LEAN_AND_MEAN

#if defined(__clang__)
#include <time.h>
#endif
#include <windows.h>
#include <process.h>
#include <winsock.h>

#else
#include <sys/time.h>	//for gettimeofday using gcc
#include <unistd.h>
#endif

#ifdef __MINGW32__
#include <Windows.h>
#endif


// ============================================================================
// useful definitions
// ============================================================================

#define INLINE __inline
#if defined(_MSC_VER)
#define getpid _getpid
#endif

#if defined(__GNUC__) && __GNUC__ >= 3
#define PREFETCH(addr) __builtin_prefetch(addr) 
#elif defined(_MSC_VER) && (_MSC_VER >= 1400)
#define PREFETCH(addr) PreFetchCacheLine(PF_TEMPORAL_LEVEL_1, addr)
#else
#define PREFETCH(addr) /* nothing */
#endif

#define MIN(a,b) ((a) < (b)? (a) : (b))
#define MAX(a,b) ((a) > (b)? (a) : (b))

#ifdef _MSC_VER
#define strto_uint64 _strtoui64
#else
#define strto_uint64 strtoull
#endif

#ifndef PRId64
// portable 64-bit formatting
#if defined(_MSC_VER) || defined(__MINGW32__)
#define PRId64 "I64d"
#define PRIu64 "I64u"
#define PRIx64 "I64x"
#elif defined(__x86_64__)
#define PRId64 "ld"
#define PRIu64 "lu"
#define PRIx64 "lx"
#elif defined(__i386__)
#define PRId64 "lld"
#define PRIu64 "llu"
#define PRIx64 "llx"
#endif
#endif

// aligned memory allocation
#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
    
#define ALIGNED_MEM __declspec(align(64))
#define align_free _mm_free

#elif defined(_MSC_VER)

#define align_free _aligned_free	
#define ALIGNED_MEM __declspec(align(64))     

#elif defined(__GNUC__)

#if defined(__MINGW64__) || defined(__MINGW32__) || defined(__MSYS__)
#define align_free _aligned_free //_mm_free
#else
#define align_free free
#endif

#define ALIGNED_MEM __attribute__((aligned(64)))

#endif

// some constants
#define LN2		0.69314718055994530942
#define INV_2_POW_48 3.5527136788005009293556213378906e-15
#define INV_2_POW_52 2.2204460492503130808472633361816e-16

#define DEFAULT_L1_CACHE_SIZE (32 * 1024)
#define DEFAULT_L2_CACHE_SIZE (512 * 1024)

// ============================================================================
// precision time
// ============================================================================

#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif


#ifdef _MSC_VER
    struct timezone
    {
        int  tz_minuteswest; /* minutes W of Greenwich */
        int  tz_dsttime;     /* type of dst correction */
    };
#endif


double ytools_difftime(struct timeval* start, struct timeval* end);


#if defined(_MSC_VER)
    int gettimeofday(struct timeval* tv, struct timezone* tz);
#endif

    extern int lock_thread_to_core(void);
    extern int unlock_thread_from_core(void);
    extern char* time_from_secs(char* str, unsigned long time);

    extern int portable_sleep(int sleep_time_ms);

// ============================================================================
// randomness
// ============================================================================

    extern void get_random_seeds(uint32_t* seed1, uint32_t* seed2);
    extern uint32_t lcg_rand_32(uint64_t* state);
    extern uint64_t lcg_rand_64(uint64_t* state);
    extern uint32_t lcg_rand_32_range(uint32_t lower, uint32_t upper, uint64_t* state);
    extern uint64_t lcg_rand_64_range(uint64_t lower, uint64_t upper, uint64_t* state);
    extern double lcg_rand_d(uint64_t* state);

// ============================================================================
// 64-bit hashing
// ============================================================================

    typedef struct
    {
        uint8_t** hashBins;
        uint64_t** hashKey;
        uint32_t* binSize;
        uint32_t numBins;
        uint32_t numBinsPow2;
        uint32_t numStored;
        uint32_t elementSizeB;
    } hash_t;

    extern hash_t* initHash(uint32_t elementSizeB, uint32_t pow2numElements);
    extern void deleteHash(hash_t* hash);
    extern void hashPut(hash_t* hash, uint8_t* element, uint64_t key);
    extern void hashGet(hash_t* hash, uint64_t key, uint8_t* element);

    extern uint64_t hash64(uint64_t in);

// ============================================================================
// allocation
// ============================================================================
    extern void aligned_free(void* ptr);
    extern void* xmalloc_align(size_t len);
    extern void* xmalloc(size_t len);
    extern void* xcalloc(size_t num, size_t len);
    extern void* xrealloc(void* iptr, size_t len);

// ============================================================================
// computer info
// ============================================================================

    typedef struct
    {
        uint32_t L1cache;
        uint32_t L2cache;

#if defined(WIN32)
        char sysname[MAX_COMPUTERNAME_LENGTH + 1];
        int sysname_sz;
#else
        char sysname[256];
        int sysname_sz;
#endif
        char idstr[256];
        int cachelinesize;
        char bSSE41Extensions;
        char BMI1;
        char AVX;
        char AVX2;
        char BMI2;
        char AVX512BW;
        char AVX512DQ;
        char AVX512ER;
        char AVX512PF;
        char AVX512CD;
        char AVX512VL;
        char AVX512IFMA;
        char AVX512F;


    } info_t;

    typedef union {
        uint32_t data;

        struct {
            uint32_t cache_type : 5;
            uint32_t cache_level : 3;
            uint32_t i_dont_care : 24;
        } s;
    } cache_type_t;

    typedef union {
        uint32_t data;

        struct {
            uint32_t line_size : 12;
            uint32_t num_lines : 10;
            uint32_t ways : 10;
        } s;
    } cache_size_t;

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
    };

    extern enum cpu_type ytools_get_cpu_type(void);
    extern void ytools_get_cache_sizes(uint32_t* level1_size_out, uint32_t* level2_size_out);
    extern void ytools_get_computer_info(info_t* info, int do_print);
    extern void ytools_get_cache_sizes(uint32_t* level1_cache, uint32_t* level2_cache);
    extern int ytools_extended_cpuid(char* idstr, int* cachelinesize, char* bSSE41Extensions,
        char* BMI1, char* AVX, char* AVX2, char* BMI2, char* AVX512F, char* AVX512BW, char* AVX512ER,
        char* AVX512PF, char* AVX512CD, char* AVX512VL, char* AVX512IFMA, char* AVX512DQ, int do_print);

// ============================================================================
// sorting
// ============================================================================

    extern int qcomp_int(const void* x, const void* y);
    extern int qcomp_uint16(const void* x, const void* y);
    extern int qcomp_uint32(const void* x, const void* y);
    extern int qcomp_uint64(const void* x, const void* y);
    extern int qcomp_double(const void* x, const void* y);
    extern uint32_t* mergesort(uint32_t* a, uint32_t* b, int sz_a, int sz_b);

// ============================================================================
// searching
// ============================================================================

    extern int bin_search_uint32(int idp, int idm, uint32_t q, uint32_t* input);
    extern int bin_search_uint64(int idp, int idm, uint64_t q, uint64_t* input);

// ============================================================================
// queue/stack
// ============================================================================

    typedef struct
    {
        uint32_t* Q;
        uint32_t len;
        uint32_t sz;
        uint32_t head;
        uint32_t tail;
        int isStack;
    } Queue_t;

    extern void clearQueue(Queue_t* Q);
    extern uint32_t peekqueue(Queue_t* Q);
    extern uint32_t dequeue(Queue_t* Q);
    extern void enqueue(Queue_t* Q, uint32_t e);
    extern Queue_t* newQueue(uint32_t sz, int isStack);

// ============================================================================
// logging and file i/o
// ============================================================================

    extern void logprint(FILE* infile, char* args, ...);
    extern void logprint_oc(const char* name, const char* method, char* args, ...);
    extern char* get_full_line(char* line, int* sz, FILE* fid);

#ifdef __cplusplus
}
#endif



#endif /* _YTOOLS_UTIL_H_ */
