/*----------------------------------------------------------------------
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
----------------------------------------------------------------------*/

#if defined(_MSC_VER) && !defined(__clang__)
#include <intrin.h>
#pragma intrinsic(__rdtsc)
#elif defined(_MSC_VER) && defined(__clang__)
#include <x86intrin.h>
#pragma intrinsic(__rdtsc)
#endif


#if (defined(__unix__) || defined(__MINGW32__) || defined(__clang__))
#define asm __asm__
#endif

#if defined(WIN32)
#include <Windows.h>
#endif

#include "ytools.h"
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

// ============================================================================
// precision time
// ============================================================================



#if defined(_MSC_VER) // && !defined(__clang__)

    /* Core aware timing on Windows, courtesy of Brian Gladman */

#if defined( _WIN64 )

#define current_processor_number GetCurrentProcessorNumber

#else

unsigned long current_processor_number(void)
{
    __asm
    {
        mov     eax, 1
        cpuid
        shr     ebx, 24
        mov     eax, ebx
    }
}

#endif

int lock_thread_to_core(void)
{
    DWORD_PTR afp, afs;

    if (GetProcessAffinityMask(GetCurrentProcess(), &afp, &afs))
    {
        afp &= (DWORD_PTR)(1 << current_processor_number());
        if (SetThreadAffinityMask(GetCurrentThread(), afp))
            return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
}

int unlock_thread_from_core(void)
{
    DWORD_PTR afp, afs;

    if (GetProcessAffinityMask(GetCurrentProcess(), &afp, &afs))
    {
        if (SetThreadAffinityMask(GetCurrentThread(), afp))
            return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
}

#endif


#if defined(_MSC_VER)

#if 0 // defined(__clang__)
int gettimeofday(struct timeval* tv, struct timezone* tz)
{
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);

    //printf("timespec_get returned sec = %lu, nsec = %lu\n", ts.tv_sec, ts.tv_nsec);

    tv->tv_sec = ts.tv_sec;
    tv->tv_usec = ts.tv_nsec / 1000;

    return 0;
}
#else
int gettimeofday(struct timeval* tv, struct timezone* tz)
{
    FILETIME ft;
    unsigned __int64 tmpres = 0;
    static int tzflag;

    if (NULL != tv)
    {
        GetSystemTimeAsFileTime(&ft);

        tmpres |= ft.dwHighDateTime;
        tmpres <<= 32;
        tmpres |= ft.dwLowDateTime;

        /*converting file time to unix epoch*/
        tmpres /= 10;  /*convert into microseconds*/
        tmpres -= DELTA_EPOCH_IN_MICROSECS;
        tv->tv_sec = (long)(tmpres / 1000000UL);
        tv->tv_usec = (long)(tmpres % 1000000UL);
    }

    if (NULL != tz)
    {
        if (!tzflag)
        {
            _tzset();
            tzflag++;
        }
        tz->tz_minuteswest = _timezone / 60;
        tz->tz_dsttime = _daylight;
    }

    return 0;
}
#endif
#endif

double ytools_difftime(struct timeval* start, struct timeval* end)
{
    double secs;
    double usecs;

    if (start->tv_sec == end->tv_sec) {
        secs = 0;
        usecs = end->tv_usec - start->tv_usec;
    }
    else {
        usecs = 1000000 - start->tv_usec;
        secs = end->tv_sec - (start->tv_sec + 1);
        usecs += end->tv_usec;
        if (usecs >= 1000000) {
            usecs -= 1000000;
            secs += 1;
        }
    }

    return secs + usecs / 1000000.;
}

char* time_from_secs(char* str, unsigned long time)
{
    // use an input scratch string
    unsigned long d;

    strcpy(str, "");
    if (time > 3600 * 24)
    {
        d = time / (3600 * 24);
        time %= (3600 * 24);
        sprintf(str, "%lu day%s ", d, d > 1 ? "s" : "");
    }
    if (time > 3600)
    {
        d = time / 3600;
        time %= 3600;
        sprintf(str, "%s%luh ", str, d);
    }
    if (time > 60)
    {
        d = time / 60;
        time %= 60;
        sprintf(str, "%s%lum ", str, d);
    }
    sprintf(str, "%s%lus", str, time);
    return str;
}


// ============================================================================
// randomness
// ============================================================================
void get_random_seeds(uint32_t *seed1, uint32_t *seed2) {

    uint32_t tmp_seed1, tmp_seed2;

#ifndef WIN32

    FILE* rand_device = fopen("/dev/urandom", "r");

    if (rand_device != NULL) {

        /* Yay! Cryptographic-quality nondeterministic randomness! */

        fread(&tmp_seed1, sizeof(uint32_t), (size_t)1, rand_device);
        fread(&tmp_seed2, sizeof(uint32_t), (size_t)1, rand_device);
        fclose(rand_device);
    }
    else

#endif
    {
        /* <Shrug> For everyone else, sample the current time,
           the high-res timer (hopefully not correlated to the
           current time), and the process ID. Multithreaded
           applications should fold in the thread ID too */
        struct timeval start;
        gettimeofday(&start, NULL);

        uint64_t high_res_time = (uint64_t)start.tv_sec * 1000000 + (uint64_t)start.tv_usec;
        tmp_seed1 = ((uint32_t)(high_res_time >> 32) ^
            (uint32_t)time(NULL)) *
#ifdef _MSC_VER
            (uint32_t)_getpid();
#else
            (uint32_t)getpid();
#endif
        tmp_seed2 = (uint32_t)high_res_time;
    }

    /* The final seeds are the result of a multiplicative
       hash of the initial seeds */

    *seed1 = tmp_seed1 * ((uint32_t)40499 * 65543);
    *seed2 = tmp_seed2 * ((uint32_t)40499 * 65543);
}

// Knuth's 64 bit MMIX LCG, using a global 64 bit state variable.
uint32_t lcg_rand_32(uint64_t* state)
{
    // advance the state of the LCG and return the appropriate result
    *state = 6364136223846793005ULL * (*state) + 1442695040888963407ULL;
    return (uint32_t)(*state);
}

uint32_t lcg_rand_32_range(uint32_t lower, uint32_t upper, uint64_t* state)
{
    // advance the state of the LCG and return the appropriate result
    *state = 6364136223846793005ULL * (*state) + 1442695040888963407ULL;
    return lower + (uint32_t)(
        (double)(upper - lower) * (double)((*state) >> 12) * INV_2_POW_52);
}

uint64_t lcg_rand_64(uint64_t* state)
{
    // advance the state of the LCG and return the appropriate result
    *state = 6364136223846793005ULL * (*state) + 1442695040888963407ULL;
    return *state;
}

uint64_t lcg_rand_64_range(uint64_t lower, uint64_t upper, uint64_t* state)
{
    // advance the state of the LCG and return the appropriate result
    *state = 6364136223846793005ULL * (*state) + 1442695040888963407ULL;
    return lower + (uint64_t)(
        (double)(upper - lower) * (double)((*state) >> 12) * INV_2_POW_52);

}

double lcg_rand_d(uint64_t* state)
{
    // Knuth's MMIX LCG converted to a floating point random number
    // between EPS and 1 by sampling the top 48 bits and dividing by 2^48.
    // values less than floating point precision (i.e., 0) are coerced to EPS.
    // http://en.wikipedia.org/wiki/Linear_congruential_generator
    double d;

    *state = 6364136223846793005ULL * (*state) + 1442695040888963407ULL;
    d = (double)(*state >> 16) * INV_2_POW_48;
    d = d <= 1e-14 ? 1e-14 : d;
    return d;
}

// ============================================================================
// 64-bit hashing
// ============================================================================

// FNV-1 hash algorithm:
// http://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function
uint64_t hash64(uint64_t in)
{
    uint64_t hash = 14695981039346656037ULL;
    uint64_t prime = 1099511628211ULL;
    uint64_t hash_mask;
    uint64_t xor;

    hash = hash * prime;
    hash_mask = 0xffffffffffffff00ULL;
    xor = hash ^ in;
    hash = (hash & hash_mask) | (xor &(~hash_mask));

    hash = hash * prime;
    hash_mask = 0xffffffffffff00ffULL;
    xor = hash ^ in;
    hash = (hash & hash_mask) | (xor &(~hash_mask));

    hash = hash * prime;
    hash_mask = 0xffffffffff00ffffULL;
    xor = hash ^ in;
    hash = (hash & hash_mask) | (xor &(~hash_mask));

    hash = hash * prime;
    hash_mask = 0xffffffff00ffffffULL;
    xor = hash ^ in;
    hash = (hash & hash_mask) | (xor &(~hash_mask));

    hash = hash * prime;
    hash_mask = 0xffffff00ffffffffULL;
    xor = hash ^ in;
    hash = (hash & hash_mask) | (xor &(~hash_mask));

    hash = hash * prime;
    hash_mask = 0xffff00ffffffffffULL;
    xor = hash ^ in;
    hash = (hash & hash_mask) | (xor &(~hash_mask));

    hash = hash * prime;
    hash_mask = 0xff00ffffffffffffULL;
    xor = hash ^ in;
    hash = (hash & hash_mask) | (xor &(~hash_mask));

    hash = hash * prime;
    hash_mask = 0x00ffffffffffffffULL;
    xor = hash ^ in;
    hash = (hash & hash_mask) | (xor &(~hash_mask));

    return hash;
}

hash_t* initHash(uint32_t elementSizeB, uint32_t pow2numElements)
{
    hash_t* hashTable;
    uint32_t i;

    hashTable = (hash_t*)xmalloc(sizeof(hash_t));
    hashTable->hashKey = (uint64_t * *)xmalloc((1ULL << pow2numElements) * sizeof(uint64_t*));
    hashTable->hashBins = (uint8_t * *)xmalloc((1ULL << pow2numElements) * sizeof(uint8_t*));
    hashTable->binSize = (uint32_t*)xcalloc((1ULL << pow2numElements), sizeof(uint32_t));
    hashTable->elementSizeB = elementSizeB;
    hashTable->numStored = 0;
    hashTable->numBinsPow2 = pow2numElements;
    hashTable->numBins = 1 << pow2numElements;

    printf("initialized hash array of size %u with elements of size %u\n",
        hashTable->numBins, hashTable->elementSizeB);

    for (i = 0; i < hashTable->numBins; i++)
    {
        hashTable->hashBins[i] = (uint8_t*)xmalloc(elementSizeB);
        hashTable->hashKey[i] = (uint64_t*)xmalloc(sizeof(uint64_t));
    }

    return hashTable;
}

void deleteHash(hash_t* hash)
{
    uint32_t i;

    for (i = 0; i < hash->numBins; i++)
    {
        free(hash->hashBins[i]);
        free(hash->hashKey[i]);
    }

    free(hash->hashBins);
    free(hash->hashKey);
    free(hash->binSize);
    hash->numStored = 0;
    hash->numBins = 0;
    free(hash);
}

void hashPut(hash_t* hash, uint8_t* element, uint64_t key)
{
    uint32_t binNum = (uint32_t)((((key)+18932479UL) * 2654435761UL) >> (64 - hash->numBinsPow2));
    //printf("hashPut into bin %u with key %u\n", binNum, key);

    if (hash->binSize[binNum] > 0)
    {
        //printf("growing bin %u size to %u\n", binNum + 1, hash->binSize[binNum] + 1);
        hash->hashBins[binNum] = (uint8_t*)xrealloc(hash->hashBins[binNum],
            (hash->binSize[binNum] + 1) * hash->elementSizeB);
        hash->hashKey[binNum] = (uint64_t*)xrealloc(hash->hashKey[binNum],
            (hash->binSize[binNum] + 1) * sizeof(uint64_t));
        memcpy(&hash->hashBins[binNum][hash->elementSizeB * hash->binSize[binNum]],
            element, hash->elementSizeB);
        hash->hashKey[binNum][hash->binSize[binNum]] = key;
        hash->binSize[binNum]++;
    }
    else
    {
        memcpy(hash->hashBins[binNum], element, hash->elementSizeB);
        hash->hashKey[binNum][hash->binSize[binNum]] = key;
        hash->binSize[binNum]++;
    }

    hash->numStored++;

    return;
}

void hashGet(hash_t* hash, uint64_t key, uint8_t* element)
{
    uint32_t binNum = (uint32_t)((((key)+18932479UL) * 2654435761UL) >> (64 - hash->numBinsPow2));
    uint32_t i;

    for (i = 0; i < hash->binSize[binNum]; i++)
    {
        if (hash->hashKey[binNum][i] == key)
        {
            memcpy(element, &hash->hashBins[binNum][hash->elementSizeB * i],
                hash->elementSizeB);
            return;
        }
    }
    return;
}

// ============================================================================
// allocation
// ============================================================================

void* xmalloc_align(size_t len)
{

#if (defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER))
    void* ptr = _mm_malloc(len, 64);

#elif (defined (_MSC_VER) || defined(__MINGW32__))
    void* ptr = _aligned_malloc(len, 64);

#elif defined (__APPLE__)
    void* ptr = malloc(len);

#elif defined (__GNUC__) // || defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
    //void* ptr = memalign(64, len);
    void* ptr;
    if ((len % 64) == 0)
        ptr = aligned_alloc(64, len);
    else
        ptr = aligned_alloc(64, len + 64 - (len % 64));

    //void* ptr;
    //
    //if ((len % sizeof(void *)) == 0)
    //    ptr = posix_memalign(&ptr, 64, len);
    //else
    //    ptr = posix_memalign(&ptr, 64, len + 64 - (len % sizeof(void*)));
    
#else
    void* ptr = malloc(len);

#endif

    if (ptr == NULL)
    {
        if (ptr == NULL) {
            printf("failed to allocate %u bytes in xmalloc_align\n", (uint32_t)len);
            exit(-1);
        }
    }
    return ptr;
}

void* xmalloc(size_t len) {
    void* ptr = malloc(len);
    if (ptr == NULL) {
        printf("failed to allocate %u bytes\n", (uint32_t)len);
        exit(-1);
    }
    return ptr;
}

void* xcalloc(size_t num, size_t len) {
    void* ptr = calloc(num, len);
    if (ptr == NULL) {
        printf("failed to calloc %u bytes\n", (uint32_t)(num * len));
        exit(-1);
    }
    return ptr;
}

void* xrealloc(void* iptr, size_t len) {
    void* ptr = realloc(iptr, len);
    if (ptr == NULL) {
        printf("failed to reallocate %u bytes\n", (uint32_t)len);
        exit(-1);
    }
    return ptr;
}

// ============================================================================
// computer info
// ============================================================================

void ytools_get_computer_info(info_t* info, int do_print)
{

#ifdef __APPLE__
    // something in extended cpuid causes a segfault on mac builds.
    // just disable it for now - this information is typically not critical for
    // program operation.
    strcpy(idstr, "N/A");
    info->L1cache = DEFAULT_L1_CACHE_SIZE;
    info->L2cache = DEFAULT_L2_CACHE_SIZE;

#else
    //read cache sizes
    ytools_get_cache_sizes(&info->L1cache, &info->L2cache);

//#if defined(WIN32)
//
//    info->sysname_sz = MAX_COMPUTERNAME_LENGTH + 1;
//    GetComputerName((LPWSTR)info->sysname, (LPDWORD)& info->sysname_sz);
//
//#else
//
//    int ret = gethostname(info->sysname, sizeof(info->sysname) / sizeof(*info->sysname));
//    info->sysname[(sizeof(info->sysname) - 1) / sizeof(*info->sysname)] = 0;	// null terminate
//    if (ret != 0)
//    {
//        printf("error occured when getting host name\n");
//        strcpy(info->sysname, "N/A");
//    }
//    info->sysname_sz = strlen(info->sysname);
//
//#endif

    ytools_extended_cpuid(info->idstr, &info->cachelinesize, &info->bSSE41Extensions,
        &info->BMI1, &info->AVX, &info->AVX2, &info->BMI2, &info->AVX512F, &info->AVX512BW, &info->AVX512ER,
        &info->AVX512PF, &info->AVX512CD, &info->AVX512VL, &info->AVX512IFMA, &info->AVX512DQ, do_print);

#endif
    return;
}


/* macro to execute the x86 CPUID instruction. Note that
   this is more verbose than it needs to be; Intel Macs reserve
   the EBX or RBX register for the PIC base address, and so
   this register cannot get clobbered by inline assembly */



#if (defined(__unix__) || defined(__MINGW32__)) && defined(__x86_64__)
#define HAS_CPUID
#define CPUID(code, a, b, c, d) 			\
		asm volatile(					\
			"movq %%rbx, %%rsi   \n\t"		\
			"cpuid               \n\t"		\
			"movl %%ebx, %1      \n\t"		\
			"movq %%rsi, %%rbx   \n\t"		\
			:"=a"(a), "=m"(b), "=c"(c), "=d"(d) 	\
			:"0"(code) : "%rsi")
#define CPUID2(code1, code2, a, b, c, d)		\
		asm volatile(					\
			"movq %%rbx, %%rsi   \n\t"		\
			"cpuid               \n\t"		\
			"movl %%ebx, %1      \n\t"		\
			"movq %%rsi, %%rbx   \n\t"		\
			:"=a"(a), "=m"(b), "=c"(c), "=d"(d) 	\
			:"0"(code1), "2"(code2) : "%rsi")

#elif defined(__unix__) && defined(__i386__)
#define HAS_CPUID
#define CPUID(code, a, b, c, d) 			\
		asm volatile(					\
			"movl %%ebx, %%esi   \n\t"		\
			"cpuid               \n\t"		\
			"movl %%ebx, %1      \n\t"		\
			"movl %%esi, %%ebx   \n\t"		\
			:"=a"(a), "=m"(b), "=c"(c), "=d"(d) 	\
			:"0"(code) : "%esi")
#define CPUID2(code1, code2, a, b, c, d) 			\
		asm volatile(					\
			"movl %%ebx, %%esi   \n\t"		\
			"cpuid               \n\t"		\
			"movl %%ebx, %1      \n\t"		\
			"movl %%esi, %%ebx   \n\t"		\
			:"=a"(a), "=m"(b), "=c"(c), "=d"(d) 	\
			:"0"(code1), "2"(code2) : "%esi")

//#elif defined(_MSC_VER) && defined(__clang__)
//#include <x86intrin.h>
//#define HAS_CPUID
//#define CPUID(__leaf, __eax, __ebx, __ecx, __edx) \
//    __asm("cpuid" : "=a"(__eax), "=b" (__ebx), "=c"(__ecx), "=d"(__edx) \
//                  : "0"(__leaf))
//#define CPUID2(code1, code2, a, b, c, d) \
//	__asm("cpuid" : "=a"(a), "=b" (b), "=c"(c), "=d"(d) \
//                  : "0"(code1), "2"(code2))
//
#elif defined(_MSC_VER) 
#include <intrin.h>
#define HAS_CPUID
#define CPUID(code, a, b, c, d)	\
	{	uint32_t _z[4]; \
		__cpuid(_z, code); \
		a = _z[0]; \
		b = _z[1]; \
		c = _z[2]; \
		d = _z[3]; \
	}
#define CPUID2(code1, code2, a, b, c, d) \
	{	uint32_t _z[4]; \
		__cpuidex(_z, code1, code2); \
		a = _z[0]; \
		b = _z[1]; \
		c = _z[2]; \
		d = _z[3]; \
	}
#endif

void ytools_get_cache_sizes(uint32_t* level1_size_out,
    uint32_t* level2_size_out)
{

    /* attempt to automatically detect the size of
       the L2 cache; this helps tune the choice of
       parameters or algorithms used in the sieve-
       based methods. It should guess right for most
       PCs and Macs when using gcc.

       Otherwise, you have the source so just fill in
       the correct number. */

    uint32_t cache_size1 = DEFAULT_L1_CACHE_SIZE;
    uint32_t cache_size2 = DEFAULT_L2_CACHE_SIZE;

#if defined(HAS_CPUID) && !defined(__APPLE__)

    /* reading the CPU-specific features of x86
       processors is a simple 57-step process.
       The following should be able to retrieve
       the L1/L2/L3 cache size of any Intel or AMD
       processor made after ~1995 */

    uint32_t a, b, c, d;
    uint8_t is_intel, is_amd;

    CPUID(0, a, b, c, d);
    is_intel = ((b & 0xff) == 'G');		/* "GenuineIntel" */
    is_amd = ((b & 0xff) == 'A');		/* "AuthenticAMD" */

    if (is_intel && a >= 2) {

        uint32_t i;
        uint8_t features[15];
        uint32_t max_special;
        uint32_t j1 = 0;
        uint32_t j2 = 0;

        /* handle newer Intel */

        if (a >= 4) {
            for (i = 0; i < 100; i++) {
                uint32_t num_sets;
                cache_type_t type;
                cache_size_t size;

                CPUID2(4, i, type.data, size.data, num_sets, d);

                /* must be data cache or unified cache */

                if (type.s.cache_type == 0)
                    break;
                else if (type.s.cache_type != 1 &&
                    type.s.cache_type != 3)
                    continue;

                d = (size.s.line_size + 1) *
                    (size.s.num_lines + 1) *
                    (size.s.ways + 1) *
                    (num_sets + 1);

                if (type.s.cache_level == 1)
                    j1 = MAX(j1, d);
                else
                    j2 = MAX(j2, d);
            }
        }

        CPUID(0x80000000, max_special, b, c, d);
        if (max_special >= 0x80000006) {
            CPUID(0x80000006, a, b, c, d);
            j2 = MAX(j2, 1024 * (c >> 16));
        }

        /* handle older Intel, possibly overriding the above */

        CPUID(2, a, b, c, d);

        features[0] = (a >> 8);
        features[1] = (a >> 16);
        features[2] = (a >> 24);
        features[3] = b;
        features[4] = (b >> 8);
        features[5] = (b >> 16);
        features[6] = (b >> 24);
        features[7] = c;
        features[8] = (c >> 8);
        features[9] = (c >> 16);
        features[10] = (c >> 24);
        features[11] = d;
        features[12] = (d >> 8);
        features[13] = (d >> 16);
        features[14] = (d >> 24);

        /* use the maximum of the (known) L2 and L3 cache sizes */

        for (i = 0; i < sizeof(features); i++) {
            switch (features[i]) {
                /* level 1 cache codes */
            case 0x06:
            case 0x0a:
            case 0x66:
                j1 = MAX(j1, 8 * 1024); break;
            case 0x08:
            case 0x0c:
            case 0x0d:
            case 0x60:
            case 0x67:
                j1 = MAX(j1, 16 * 1024); break;
            case 0x0e:
                j1 = MAX(j1, 24 * 1024); break;
            case 0x09:
            case 0x2c:
            case 0x30:
            case 0x68:
                j1 = MAX(j1, 32 * 1024); break;

                /* level 2 and level 3 cache codes */
            case 0x41:
            case 0x79:
                j2 = MAX(j2, 128 * 1024); break;
            case 0x21:
            case 0x42:
            case 0x7a:
            case 0x82:
                j2 = MAX(j2, 256 * 1024); break;
            case 0x22:
            case 0x43:
            case 0x7b:
            case 0x7f:
            case 0x80:
            case 0x83:
            case 0x86:
                j2 = MAX(j2, 512 * 1024); break;
            case 0x23:
            case 0x44:
            case 0x78:
            case 0x7c:
            case 0x84:
            case 0x87:
                j2 = MAX(j2, 1 * 1024 * 1024); break;
            case 0x25:
            case 0x45:
            case 0x7d:
            case 0x85:
                j2 = MAX(j2, 2 * 1024 * 1024); break;
            case 0x48:
                j2 = MAX(j2, 3 * 1024 * 1024); break;
            case 0x29:
            case 0x46:
            case 0x49:
                j2 = MAX(j2, 4 * 1024 * 1024); break;
            case 0x4a:
            case 0x4e:
                j2 = MAX(j2, 6 * 1024 * 1024); break;
            case 0x47:
            case 0x4b:
            case 0xe4:
                j2 = MAX(j2, 8 * 1024 * 1024); break;
            case 0x4c:
            case 0xea:
                j2 = MAX(j2, 12 * 1024 * 1024); break;
            case 0x4d:
                j2 = MAX(j2, 16 * 1024 * 1024); break;
            case 0xeb:
                j2 = MAX(j2, 18 * 1024 * 1024); break;
            case 0xec:
                j2 = MAX(j2, 24 * 1024 * 1024); break;
            }
        }
        if (j1 > 0)
            cache_size1 = j1;
        if (j2 > 0)
            cache_size2 = j2;
    }
    else if (is_amd) {

        uint32_t max_special;
        CPUID(0x80000000, max_special, b, c, d);

        if (max_special >= 0x80000005) {
            CPUID(0x80000005, a, b, c, d);
            cache_size1 = 1024 * (c >> 24);

            if (max_special >= 0x80000006) {
                CPUID(0x80000006, a, b, c, d);
                cache_size2 = MAX(1024 * (c >> 16),
                    512 * 1024 * (d >> 18));
            }
        }
    }
#endif

    * level1_size_out = cache_size1;
    *level2_size_out = cache_size2;
}

enum cpu_type ytools_get_cpu_type(void) {

    enum cpu_type cpu = cpu_generic;

#if defined(HAS_CPUID)
    uint32_t a, b, c, d;

    CPUID(0, a, b, c, d);
    if ((b & 0xff) == 'G') {	/* "GenuineIntel" */

        uint8_t family, model;

        switch (a) {
        case 1:
            cpu = cpu_pentium;
            break;
        case 2:
            CPUID(1, a, b, c, d);
            family = (a >> 8) & 0xf;
            model = (a >> 4) & 0xf;
            if (family == 6) {
                if (model == 9 || model == 13)
                    cpu = cpu_pentium_m;
                else
                    cpu = cpu_pentium2;
            }
            else if (family == 15) {
                cpu = cpu_pentium4;
            }
            break;
        case 3:
            cpu = cpu_pentium3;
            break;
        case 5:
        case 6:
            cpu = cpu_pentium4;
            break;
        default:
            /* a = 10+; some subspecies of core or core2 */
            CPUID(1, a, b, c, d);
            family = (a >> 8) & 0xf;
            model = (a >> 4) & 0xf;

            if (model >= 10)
                cpu = cpu_core; //cpu_nehalem;
            else
                cpu = cpu_core;

            break;
        }
    }
    else if ((b & 0xff) == 'A') {		/* "AuthenticAMD" */

        uint8_t family, model;

        CPUID(1, a, b, c, d);
        family = (a >> 8) & 0xf;
        model = (a >> 4) & 0xf;
        if (family == 15)
            cpu = cpu_opteron;
        else if (family == 6) {
            CPUID(0x80000001, a, b, c, d);
            if (d & 0x1000000)		/* full SSE */
                cpu = cpu_athlon_xp;
            else				/* partial SSE */
                cpu = cpu_athlon;
        }
    }
#endif

    return cpu;
}


// http://msdn.microsoft.com/en-us/library/hskdteyh.aspx
// cpuid.cpp 
// processor: x86, x64
// Use the __cpuid intrinsic to get information about a CPU
// modified for c compliers and use of CPUID macros 
//		- brb, 10/26/10



int ytools_extended_cpuid(char* idstr, int* cachelinesize, char* bSSE41Extensions,
    char* BMI1, char* AVX, char* AVX2, char* BMI2, char* AVX512F, char* AVX512BW, char* AVX512ER,
    char* AVX512PF, char* AVX512CD, char* AVX512VL, char* AVX512IFMA, char* AVX512DQ, int do_print)
{
    char CPUString[0x20];
    char CPUBrandString[0x40];
    int CPUInfo[4] = { -1 };
    int nSteppingID = 0;
    int nModel = 0;
    int nFamily = 0;
    int nProcessorType = 0;
    int nExtendedmodel = 0;
    int nExtendedfamily = 0;
    int nBrandIndex = 0;
    int nCLFLUSHcachelinesize = 0;
    int nLogicalProcessors = 0;
    int nAPICPhysicalID = 0;
    int nFeatureInfo = 0;
    int nCacheLineSize = 0;
    int nL2Associativity = 0;
    int nCacheSizeK = 0;
    int nPhysicalAddress = 0;
    int nVirtualAddress = 0;
    int nRet = 0;

    int nCores = 0;
    int nCacheType = 0;
    int nCacheLevel = 0;
    int nMaxThread = 0;
    int nSysLineSize = 0;
    int nPhysicalLinePartitions = 0;
    int nWaysAssociativity = 0;
    int nNumberSets = 0;

    unsigned    nIds, nExIds, i;

    char    bSSE3Instructions = 0;
    char    bMONITOR_MWAIT = 0;
    char    bCPLQualifiedDebugStore = 0;
    char    bVirtualMachineExtensions = 0;
    char    bEnhancedIntelSpeedStepTechnology = 0;
    char    bThermalMonitor2 = 0;
    char    bSupplementalSSE3 = 0;
    char    bL1ContextID = 0;
    char    bCMPXCHG16B = 0;
    char    bxTPRUpdateControl = 0;
    char    bPerfDebugCapabilityMSR = 0;
    //char    bSSE41Extensions = 0;
    char    bSSE42Extensions = 0;
    char    bPOPCNT = 0;

    char    bMultithreading = 0;

    char    bLAHF_SAHFAvailable = 0;
    char    bCmpLegacy = 0;
    char    bSVM = 0;
    char    bExtApicSpace = 0;
    char    bAltMovCr8 = 0;
    char    bLZCNT = 0;
    char    bSSE4A = 0;
    char    bMisalignedSSE = 0;
    char    bPREFETCH = 0;
    char    bSKINITandDEV = 0;
    char    bSYSCALL_SYSRETAvailable = 0;
    char    bExecuteDisableBitAvailable = 0;
    char    bMMXExtensions = 0;
    char    bFFXSR = 0;
    char    b1GBSupport = 0;
    char    bRDTSCP = 0;
    char    b64Available = 0;
    char    b3DNowExt = 0;
    char    b3DNow = 0;
    char    bNestedPaging = 0;
    char    bLBRVisualization = 0;
    char    bFP128 = 0;
    char    bMOVOptimization = 0;

    char    bSelfInit = 0;
    char    bFullyAssociative = 0;

    const char* szFeatures[] =
    {
        "x87 FPU On Chip",
        "Virtual-8086 Mode Enhancement",
        "Debugging Extensions",
        "Page Size Extensions",
        "Time Stamp Counter",
        "RDMSR and WRMSR Support",
        "Physical Address Extensions",
        "Machine Check Exception",
        "CMPXCHG8B Instruction",
        "APIC On Chip",
        "Unknown1",
        "SYSENTER and SYSEXIT",
        "Memory Type Range Registers",
        "PTE Global Bit",
        "Machine Check Architecture",
        "Conditional Move/Compare Instruction",
        "Page Attribute Table",
        "36-bit Page Size Extension",
        "Processor Serial Number",
        "CFLUSH Extension",
        "Unknown2",
        "Debug Store",
        "Thermal Monitor and Clock Ctrl",
        "MMX Technology",
        "FXSAVE/FXRSTOR",
        "SSE Extensions",
        "SSE2 Extensions",
        "Self Snoop",
        "Multithreading Technology",
        "Thermal Monitor",
        "Unknown4",
        "Pending Break Enable"
    };


    *bSSE41Extensions = 0;

    // __cpuid with an InfoType argument of 0 returns the number of
    // valid Ids in CPUInfo[0] and the CPU identification string in
    // the other three array elements. The CPU identification string is
    // not in linear order. The code below arranges the information 
    // in a human readable form.
    CPUID(0, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
    //__cpuid(CPUInfo, 0);
    nIds = CPUInfo[0];
    memset(CPUString, 0, sizeof(CPUString));
    *((int*)CPUString) = CPUInfo[1];
    *((int*)(CPUString + 4)) = CPUInfo[3];
    *((int*)(CPUString + 8)) = CPUInfo[2];

    // Get the information associated with each valid Id
    for (i = 0; i <= nIds; ++i)
    {
        CPUID(i, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
        //__cpuid(CPUInfo, i);

        if (do_print)
        {
            printf("\nFor InfoType %d\n", i);
            printf("CPUInfo[0] = 0x%x\n", CPUInfo[0]);
            printf("CPUInfo[1] = 0x%x\n", CPUInfo[1]);
            printf("CPUInfo[2] = 0x%x\n", CPUInfo[2]);
            printf("CPUInfo[3] = 0x%x\n", CPUInfo[3]);
        }

        // Interpret CPU feature information.
        if (i == 1)
        {
            nSteppingID = CPUInfo[0] & 0xf;
            nModel = (CPUInfo[0] >> 4) & 0xf;
            nFamily = (CPUInfo[0] >> 8) & 0xf;
            nProcessorType = (CPUInfo[0] >> 12) & 0x3;
            nExtendedmodel = (CPUInfo[0] >> 16) & 0xf;
            nExtendedfamily = (CPUInfo[0] >> 20) & 0xff;
            nBrandIndex = CPUInfo[1] & 0xff;
            *cachelinesize = nCLFLUSHcachelinesize = ((CPUInfo[1] >> 8) & 0xff) * 8;
            nLogicalProcessors = ((CPUInfo[1] >> 16) & 0xff);
            nAPICPhysicalID = (CPUInfo[1] >> 24) & 0xff;
            bSSE3Instructions = (CPUInfo[2] & 0x1) || 0;
            bMONITOR_MWAIT = (CPUInfo[2] & 0x8) || 0;
            bCPLQualifiedDebugStore = (CPUInfo[2] & 0x10) || 0;
            bVirtualMachineExtensions = (CPUInfo[2] & 0x20) || 0;
            bEnhancedIntelSpeedStepTechnology = (CPUInfo[2] & 0x80) || 0;
            bThermalMonitor2 = (CPUInfo[2] & 0x100) || 0;
            bSupplementalSSE3 = (CPUInfo[2] & 0x200) || 0;
            bL1ContextID = (CPUInfo[2] & 0x300) || 0;
            bCMPXCHG16B = (CPUInfo[2] & 0x2000) || 0;
            bxTPRUpdateControl = (CPUInfo[2] & 0x4000) || 0;
            bPerfDebugCapabilityMSR = (CPUInfo[2] & 0x8000) || 0;
            *bSSE41Extensions = (CPUInfo[2] & 0x80000) || 0;
            bSSE42Extensions = (CPUInfo[2] & 0x100000) || 0;
            bPOPCNT = (CPUInfo[2] & 0x800000) || 0;
            *AVX = (CPUInfo[2] & 0x10000000) || 0;
            nFeatureInfo = CPUInfo[3];
            bMultithreading = (nFeatureInfo & (1 << 28)) || 0;
        }
    }

    // Calling __cpuid with 0x80000000 as the InfoType argument
    // gets the number of valid extended IDs.
    CPUID(0x80000000, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
    //__cpuid(CPUInfo, 0x80000000);
    nExIds = CPUInfo[0];
    memset(CPUBrandString, 0, sizeof(CPUBrandString));

    // Get the information associated with each extended ID.
    for (i = 0x80000000; i <= nExIds; ++i)
    {
        CPUID(i, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
        //__cpuid(CPUInfo, i);
        if (do_print)
        {
            printf("\nFor InfoType %x\n", i);
            printf("CPUInfo[0] = 0x%x\n", CPUInfo[0]);
            printf("CPUInfo[1] = 0x%x\n", CPUInfo[1]);
            printf("CPUInfo[2] = 0x%x\n", CPUInfo[2]);
            printf("CPUInfo[3] = 0x%x\n", CPUInfo[3]);
        }

        if (i == 0x80000001)
        {
            bLAHF_SAHFAvailable = (CPUInfo[2] & 0x1) || 0;
            bCmpLegacy = (CPUInfo[2] & 0x2) || 0;
            bSVM = (CPUInfo[2] & 0x4) || 0;
            bExtApicSpace = (CPUInfo[2] & 0x8) || 0;
            bAltMovCr8 = (CPUInfo[2] & 0x10) || 0;
            bLZCNT = (CPUInfo[2] & 0x20) || 0;
            bSSE4A = (CPUInfo[2] & 0x40) || 0;
            bMisalignedSSE = (CPUInfo[2] & 0x80) || 0;
            bPREFETCH = (CPUInfo[2] & 0x100) || 0;
            bSKINITandDEV = (CPUInfo[2] & 0x1000) || 0;
            bSYSCALL_SYSRETAvailable = (CPUInfo[3] & 0x800) || 0;
            bExecuteDisableBitAvailable = (CPUInfo[3] & 0x10000) || 0;
            bMMXExtensions = (CPUInfo[3] & 0x40000) || 0;
            bFFXSR = (CPUInfo[3] & 0x200000) || 0;
            b1GBSupport = (CPUInfo[3] & 0x400000) || 0;
            bRDTSCP = (CPUInfo[3] & 0x8000000) || 0;
            b64Available = (CPUInfo[3] & 0x20000000) || 0;
            b3DNowExt = (CPUInfo[3] & 0x40000000) || 0;
            b3DNow = (CPUInfo[3] & 0x80000000) || 0;
        }

        // Interpret CPU brand string and cache information.
        if (i == 0x80000002)
            memcpy(CPUBrandString, CPUInfo, sizeof(CPUInfo));
        else if (i == 0x80000003)
            memcpy(CPUBrandString + 16, CPUInfo, sizeof(CPUInfo));
        else if (i == 0x80000004)
            memcpy(CPUBrandString + 32, CPUInfo, sizeof(CPUInfo));
        else if (i == 0x80000006)
        {
            nCacheLineSize = CPUInfo[2] & 0xff;
            nL2Associativity = (CPUInfo[2] >> 12) & 0xf;
            nCacheSizeK = (CPUInfo[2] >> 16) & 0xffff;
        }
        else if (i == 0x80000008)
        {
            nPhysicalAddress = CPUInfo[0] & 0xff;
            nVirtualAddress = (CPUInfo[0] >> 8) & 0xff;
        }
        else if (i == 0x8000000A)
        {
            bNestedPaging = (CPUInfo[3] & 0x1) || 0;
            bLBRVisualization = (CPUInfo[3] & 0x2) || 0;
        }
        else if (i == 0x8000001A)
        {
            bFP128 = (CPUInfo[0] & 0x1) || 0;
            bMOVOptimization = (CPUInfo[0] & 0x2) || 0;
        }
    }
    strncpy(idstr, CPUBrandString, 255);
    // Display all the information in user-friendly format.
    if (do_print)
        printf("\n\nCPU String: %s\n", CPUString);

    if (nIds >= 1)
    {
        if (do_print)
        {
            if (nSteppingID)
                printf("Stepping ID = %d\n", nSteppingID);
            if (nModel)
                printf("Model = %d\n", nModel);
            if (nFamily)
                printf("Family = %d\n", nFamily);
            if (nProcessorType)
                printf("Processor Type = %d\n", nProcessorType);
            if (nExtendedmodel)
                printf("Extended model = %d\n", nExtendedmodel);
            if (nExtendedfamily)
                printf("Extended family = %d\n", nExtendedfamily);
            if (nBrandIndex)
                printf("Brand Index = %d\n", nBrandIndex);
            if (nCLFLUSHcachelinesize)
                printf("CLFLUSH cache line size = %d\n",
                    nCLFLUSHcachelinesize);
            if (bMultithreading && (nLogicalProcessors > 0))
                printf("Logical Processor Count = %d\n", nLogicalProcessors);
            if (nAPICPhysicalID)
                printf("APIC Physical ID = %d\n", nAPICPhysicalID);

            if (nFeatureInfo || bSSE3Instructions ||
                bMONITOR_MWAIT || bCPLQualifiedDebugStore ||
                bVirtualMachineExtensions || bEnhancedIntelSpeedStepTechnology ||
                bThermalMonitor2 || bSupplementalSSE3 || bL1ContextID ||
                bCMPXCHG16B || bxTPRUpdateControl || bPerfDebugCapabilityMSR ||
                *bSSE41Extensions || bSSE42Extensions || bPOPCNT ||
                bLAHF_SAHFAvailable || bCmpLegacy || bSVM ||
                bExtApicSpace || bAltMovCr8 ||
                bLZCNT || bSSE4A || bMisalignedSSE ||
                bPREFETCH || bSKINITandDEV || bSYSCALL_SYSRETAvailable ||
                bExecuteDisableBitAvailable || bMMXExtensions || bFFXSR || b1GBSupport ||
                bRDTSCP || b64Available || b3DNowExt || b3DNow || bNestedPaging ||
                bLBRVisualization || bFP128 || bMOVOptimization)
            {
                printf("\nThe following features are supported:\n");

                if (bSSE3Instructions)
                    printf("\tSSE3\n");
                if (bMONITOR_MWAIT)
                    printf("\tMONITOR/MWAIT\n");
                if (bCPLQualifiedDebugStore)
                    printf("\tCPL Qualified Debug Store\n");
                if (bVirtualMachineExtensions)
                    printf("\tVirtual Machine Extensions\n");
                if (bEnhancedIntelSpeedStepTechnology)
                    printf("\tEnhanced Intel SpeedStep Technology\n");
                if (bThermalMonitor2)
                    printf("\tThermal Monitor 2\n");
                if (bSupplementalSSE3)
                    printf("\tSupplemental Streaming SIMD Extensions 3\n");
                if (bL1ContextID)
                    printf("\tL1 Context ID\n");
                if (bCMPXCHG16B)
                    printf("\tCMPXCHG16B Instruction\n");
                if (bxTPRUpdateControl)
                    printf("\txTPR Update Control\n");
                if (bPerfDebugCapabilityMSR)
                    printf("\tPerf\\Debug Capability MSR\n");
                if (*bSSE41Extensions)
                    printf("\tSSE4.1 Extensions\n");
                if (bSSE42Extensions)
                    printf("\tSSE4.2 Extensions\n");
                if (*AVX)
                    printf("\tAVX Extensions\n");
                if (bPOPCNT)
                    printf("\tPPOPCNT Instruction\n");

                i = 0;
                nIds = 1;
                while (i < (sizeof(szFeatures) / sizeof(const char*)))
                {
                    if (nFeatureInfo & nIds)
                    {
                        printf("\t");
                        printf("%s", szFeatures[i]);
                        printf("\n");
                    }

                    nIds <<= 1;
                    ++i;
                }
                if (bLAHF_SAHFAvailable)
                    printf("\tLAHF/SAHF in 64-bit mode\n");
                if (bCmpLegacy)
                    printf("\tCore multi-processing legacy mode\n");
                if (bSVM)
                    printf("\tSecure Virtual Machine\n");
                if (bExtApicSpace)
                    printf("\tExtended APIC Register Space\n");
                if (bAltMovCr8)
                    printf("\tAltMovCr8\n");
                if (bLZCNT)
                    printf("\tLZCNT instruction\n");
                if (bSSE4A)
                    printf("\tSSE4A (EXTRQ, INSERTQ, MOVNTSD, MOVNTSS)\n");
                if (bMisalignedSSE)
                    printf("\tMisaligned SSE mode\n");
                if (bPREFETCH)
                    printf("\tPREFETCH and PREFETCHW Instructions\n");
                if (bSKINITandDEV)
                    printf("\tSKINIT and DEV support\n");
                if (bSYSCALL_SYSRETAvailable)
                    printf("\tSYSCALL/SYSRET in 64-bit mode\n");
                if (bExecuteDisableBitAvailable)
                    printf("\tExecute Disable Bit\n");
                if (bMMXExtensions)
                    printf("\tExtensions to MMX Instructions\n");
                if (bFFXSR)
                    printf("\tFFXSR\n");
                if (b1GBSupport)
                    printf("\t1GB page support\n");
                if (bRDTSCP)
                    printf("\tRDTSCP instruction\n");
                if (b64Available)
                    printf("\t64 bit Technology\n");
                if (b3DNowExt)
                    printf("\t3Dnow Ext\n");
                if (b3DNow)
                    printf("\t3Dnow! instructions\n");
                if (bNestedPaging)
                    printf("\tNested Paging\n");
                if (bLBRVisualization)
                    printf("\tLBR Visualization\n");
                if (bFP128)
                    printf("\tFP128 optimization\n");
                if (bMOVOptimization)
                    printf("\tMOVU Optimization\n");
            }
        }
    }

    if (nExIds >= 0x80000004 && do_print)
        printf("\nCPU Brand String: %s\n", CPUBrandString);

    if (nExIds >= 0x80000006 && do_print)
    {
        printf("Cache Line Size = %d\n", nCacheLineSize);
        printf("L2 Associativity = %d\n", nL2Associativity);
        printf("Cache Size = %dK\n", nCacheSizeK);
    }


    for (i = 0;; i++)
    {
        CPUID2(0x4, i, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
        //__cpuidex(CPUInfo, 0x4, i);
        if (!(CPUInfo[0] & 0xf0)) break;

        if (i == 0)
        {
            nCores = CPUInfo[0] >> 26;
            if (do_print)
                printf("\n\nNumber of Cores = %d\n", nCores + 1);
        }

        nCacheType = (CPUInfo[0] & 0x1f);
        nCacheLevel = (CPUInfo[0] & 0xe0) >> 5;
        bSelfInit = (CPUInfo[0] & 0x100) >> 8;
        bFullyAssociative = (CPUInfo[0] & 0x200) >> 9;
        nMaxThread = (CPUInfo[0] & 0x03ffc000) >> 14;
        nSysLineSize = (CPUInfo[1] & 0x0fff);
        nPhysicalLinePartitions = (CPUInfo[1] & 0x03ff000) >> 12;
        nWaysAssociativity = (CPUInfo[1]) >> 22;
        nNumberSets = CPUInfo[2];

        if (do_print)
        {
            printf("\n");

            printf("ECX Index %d\n", i);
            switch (nCacheType)
            {
            case 0:
                printf("   Type: Null\n");
                break;
            case 1:
                printf("   Type: Data Cache\n");
                break;
            case 2:
                printf("   Type: Instruction Cache\n");
                break;
            case 3:
                printf("   Type: Unified Cache\n");
                break;
            default:
                printf("   Type: Unknown\n");
            }

            printf("   Level = %d\n", nCacheLevel + 1);
            if (bSelfInit)
            {
                printf("   Self Initializing\n");
            }
            else
            {
                printf("   Not Self Initializing\n");
            }
            if (bFullyAssociative)
            {
                printf("   Is Fully Associatve\n");
            }
            else
            {
                printf("   Is Not Fully Associatve\n");
            }
            printf("   Max Threads = %d\n",
                nMaxThread + 1);
            printf("   System Line Size = %d\n",
                nSysLineSize + 1);
            printf("   Physical Line Partions = %d\n",
                nPhysicalLinePartitions + 1);
            printf("   Ways of Associativity = %d\n",
                nWaysAssociativity + 1);
            printf("   Number of Sets = %d\n",
                nNumberSets + 1);
        }
    }

    CPUID2(0x7, 0, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);

    if (do_print)
    {
        printf("EAX=7 CPUID feature bits:\n");
        printf("EAX=%08x\n", CPUInfo[0]);
        printf("EBX=%08x\n", CPUInfo[1]);
        printf("ECX=%08x\n", CPUInfo[2]);
        printf("EDX=%08x\n", CPUInfo[3]);
    }

    *BMI1 = (CPUInfo[1] & (1 << 3)) || 0;
    *AVX2 = (CPUInfo[1] & (1 << 5)) || 0;
    *BMI2 = (CPUInfo[1] & (1 << 8)) || 0;
    *AVX512F = (CPUInfo[1] & (1 << 16)) || 0;
    *AVX512DQ = (CPUInfo[1] & (1 << 17)) || 0;
    *AVX512IFMA = (CPUInfo[1] & (1 << 21)) || 0;
    *AVX512PF = (CPUInfo[1] & (1 << 26)) || 0;
    *AVX512ER = (CPUInfo[1] & (1 << 27)) || 0;
    *AVX512CD = (CPUInfo[1] & (1 << 28)) || 0;
    *AVX512BW = (CPUInfo[1] & (1 << 30)) || 0;
    *AVX512VL = (CPUInfo[1] & (1 << 31)) || 0;

    if ((*BMI1) && do_print)
        printf("\n\n\tBMI1 Extensions\n");

    if ((*AVX2) && do_print)
        printf("\n\n\tAVX2 Extensions\n");

    if ((*BMI2) && do_print)
        printf("\n\n\tBMI2 Extensions\n");

    if ((*AVX512F) && do_print)
        printf("\n\n\tAVX512F Extensions\n");

    if ((*AVX512DQ) && do_print)
        printf("\n\n\tAVX512DQ Extensions\n");

    if ((*AVX512IFMA) && do_print)
        printf("\n\n\tAVX512IFMA Extensions\n");

    if ((*AVX512PF) && do_print)
        printf("\n\n\tAVX512PF Extensions\n");

    if ((*AVX512ER) && do_print)
        printf("\n\n\tAVX512ER Extensions\n");

    if ((*AVX512CD) && do_print)
        printf("\n\n\tAVX512CD Extensions\n");

    if ((*AVX512BW) && do_print)
        printf("\n\n\tAVX512BW Extensions\n");

    if ((*AVX512VL) && do_print)
        printf("\n\n\tAVX512VL Extensions\n");

    return  nRet;
}



// ============================================================================
// sorting
// ============================================================================

int qcomp_uint16(const void *x, const void *y)
{
	uint16_t *xx = (uint16_t *)x;
	uint16_t *yy = (uint16_t *)y;
	
	if (*xx > *yy)
		return 1;
	else if (*xx == *yy)
		return 0;
	else
		return -1;
}

int qcomp_uint32(const void *x, const void *y)
{
	uint32_t *xx = (uint32_t *)x;
	uint32_t *yy = (uint32_t *)y;
	
	if (*xx > *yy)
		return 1;
	else if (*xx == *yy)
		return 0;
	else
		return -1;
}

int qcomp_uint64(const void *x, const void *y)
{
	uint64_t *xx = (uint64_t *)x;
	uint64_t *yy = (uint64_t *)y;
	
	if (*xx > *yy)
		return 1;
	else if (*xx == *yy)
		return 0;
	else
		return -1;
}

int qcomp_int(const void *x, const void *y)
{
	int *xx = (int *)x;
	int *yy = (int *)y;
	
	if (*xx > *yy)
		return 1;
	else if (*xx == *yy)
		return 0;
	else
		return -1;
}

int qcomp_double(const void *x, const void *y)
{
	double *xx = (double *)x;
	double *yy = (double *)y;

	if (*xx > *yy)
		return 1;
	else if (*xx == *yy)
		return 0;
	else
		return -1;
}



uint32_t* mergesort(uint32_t* a, uint32_t* b, int sz_a, int sz_b)
{
    uint32_t* c = (uint32_t*)malloc((sz_a + sz_b) * sizeof(uint32_t));
    int i = 0, j = 0, k = 0;

    while ((i < sz_a) && (j < sz_b)) {
        if (a[i] < b[j]) {
            c[k++] = a[i++];
        }
        else if (a[i] > b[j]) {
            c[k++] = b[j++];
        }
        else {
            c[k++] = a[i++];
            c[k++] = b[j++];
        }
    }

    while (i < sz_a)
        c[k++] = a[i++];

    while (j < sz_b)
        c[k++] = b[j++];

    return c;
}

// ============================================================================
// searching
// ============================================================================


int bin_search_uint32(int idp, int idm, uint32_t q, uint32_t *input)
{
	int next = (idp + idm) / 2;
	
	while ((idp - idm) > 10)
	{
		if (input[next] > q)
		{
			idp = next;
			next = (next + idm) / 2;							
		}
		else					
		{
			idm = next;
			next = (idp + next) / 2;							
		}
	}

	for (next = idm; next < idm + 10; next++)
		if (input[next] == q)
			return next;

	if (input[next] != q)
		next = -1;

	return next;
}

int bin_search_uint64(int idp, int idm, uint64_t q, uint64_t* input)
{
    int next = (idp + idm) / 2;

    while ((idp - idm) > 10)
    {
        if (input[next] > q)
        {
            idp = next;
            next = (next + idm) / 2;
        }
        else
        {
            idm = next;
            next = (idp + next) / 2;
        }
    }

    for (next = idm; next < idm + 10; next++)
        if (input[next] == q)
            return next;

    if (input[next] != q)
        next = -1;

    return next;
}

// ============================================================================
// queue/stack
// ============================================================================

Queue_t* newQueue(uint32_t sz, int isStack)
{
    Queue_t* Q = (Queue_t*)malloc(sizeof(Queue_t));
    Q->Q = (uint32_t*)malloc(sz * sizeof(uint32_t));
    Q->sz = sz;
    Q->head = 0;
    Q->tail = 0;
    Q->len = 0;
    Q->isStack = isStack;
    return Q;
}

void enqueue(Queue_t* Q, uint32_t e)
{
    Q->Q[Q->tail++] = e;
    Q->len++;

    if (Q->tail == Q->sz)
    {
        Q->tail = 0;
    }

    if (Q->len >= Q->sz)
    {
        printf("warning: Q overflowed\n");
    }
    return;
}

uint32_t dequeue(Queue_t* Q)
{
    uint32_t e = -1;

    if (Q->len > 0)
    {
        if (Q->isStack)
        {
            e = Q->Q[Q->tail];

            if (Q->tail > 0)
            {
                Q->tail--;
            }
        }
        else
        {
            e = Q->Q[Q->head];
            Q->head++;

            if (Q->head == Q->sz)
            {
                Q->head = 0;
            }
        }

        Q->len--;
    }
    else
    {
        printf("warning: attempted to dequeue from an empty queue\n");
    }

    return e;
}

uint32_t peekqueue(Queue_t* Q)
{
    uint32_t e = -1;
    if (Q->len > 0)
    {
        if (Q->isStack)
            e = Q->Q[Q->tail];
        else
            e = Q->Q[Q->head];
    }
    return e;
}

void clearQueue(Queue_t* Q)
{
    free(Q->Q);
    Q->len = 0;
    Q->sz = 0;
    Q->head = 0;
    Q->tail = 0;

    return;
}

// ============================================================================
// logging and file i/o
// ============================================================================

void logprint(FILE* infile, char* args, ...)
{
    va_list ap;
    time_t curtime;
    struct tm* loctime;
    char s[256];

    if (infile == NULL)
        return;

    curtime = time(NULL);
    loctime = localtime(&curtime);

    if (infile == NULL)
        infile = stdout;

    //args has at least one argument, which may be NULL
    va_start(ap, args);

    //print the date and version stamp
    strftime(s, 256, "%m/%d/%y %H:%M:%S", loctime);
    //fprintf(infile,"%s v%d.%02d @ %s, ",s,MAJOR_VER,MINOR_VER,sysname);
    //fprintf(infile,"%s v%s @ %s, ",s,VERSION_STRING,sysname);
    //fprintf(infile, "%s v%s, ", s, VERSION_STRING);
    fprintf(infile, "%s, ", s);
    vfprintf(infile, args, ap);

    va_end(ap);
    return;
}

void logprint_oc(const char* name, const char* method, char* args, ...)
{
    //print formatted string to logfile, given the file name
    //and access method
    va_list ap;
    time_t curtime;
    struct tm* loctime;
    char s[256];
    FILE* logfile;

    if (strlen(name) == 0)
        return;

    logfile = fopen(name, method);
    if (logfile == NULL)
    {
        printf("fopen error: %s\n", strerror(errno));
        printf("could not open logfile for method %c\n", method[0]);
        logfile = stdout;
    }

    curtime = time(NULL);
    loctime = localtime(&curtime);

    //args has at least one argument, which may be NULL
    va_start(ap, args);

    //print the date and version stamp
    strftime(s, 256, "%m/%d/%y %H:%M:%S", loctime);
    //fprintf(infile,"%s v%d.%02d @ %s, ",s,MAJOR_VER,MINOR_VER,sysname);
    //fprintf(logfile,"%s v%s @ %s, ",s,VERSION_STRING,sysname);
    //fprintf(logfile, "%s v%s, ", s, VERSION_STRING);
    fprintf(logfile, "%s, ", s);
    vfprintf(logfile, args, ap);

    va_end(ap);

    if (logfile != stdout)
        fclose(logfile);

    return;
}

char* get_full_line(char* line, int *sz, FILE* fid)
{
    // fid is open to a valid stream.
    // if line is allocated (non-null), then use it, otherwise allocate it.
    // return the (possibly modified) pointer to the next full line in the file.
    // also return the updated size in the *sz pointer.
    char* ptr;
    char tmpline[1024];

    if (line == NULL)
    {
        *sz = 1024;
        line = (char*)malloc(1024 * sizeof(char));
    }

    strcpy(line, "");
    do
    {
        int j;
        while (1)
        {
            ptr = fgets(tmpline, 1024, fid);
            strcpy(line + strlen(line), tmpline);

            // stop if we didn't read anything
            if ((feof(fid)) || (ptr == NULL))
            {
                free(line);
                return NULL;
            }

            // if we got the end of the line, stop reading
            if ((line[strlen(line) - 1] == 0xa) ||
                (line[strlen(line) - 1] == 0xd))
                break;

            // else reallocate the buffer and get some more
            *sz += 1024;
            line = (char*)realloc(line, (strlen(line) + 1024) * sizeof(char));
        }

        // remove trailing LF and CRs from line
        for (j = (int)strlen(line) - 1; j > 0; j--)
        {
            switch (line[j])
            {
            case 13:
            case 10:
                line[j] = '\0';
                break;
            default:
                j = 0;
                break;
            }
        }
    } while (strlen(line) == 0);

    return line;
}



// ============================================================================
// combinatorics
// ============================================================================

uint64_t choose(int n, int k)
{
    // compute binomial: n choose k
    // multiplication method:
    // prod(i=1,k,(n+1-i)/i)
    uint64_t p = 1;
    int i;
    for (i = 1; i <= k; i++)
    {
        p = p * (n + 1 - i) / i;
    }
    return p;
}

/** [combination c n p x]
 * get the [x]th lexicographically ordered set of [p] elements in [n]
 * output is in [c], and should be sizeof(int)*[p] */
void combination(int* c, int n, int p, int x) {
    int i, r, k = 0;
    for (i = 0; i < p - 1; i++) {
        c[i] = (i != 0) ? c[i - 1] : 0;
        do {
            c[i]++;
            r = choose(n - c[i], p - (i + 1));
            k = k + r;
        } while (k < x);
        k = k - r;
    }
    c[p - 1] = c[p - 2] + x - k;
}

/** [combination c n p x]
 * get all lexicographically ordered sets of [p] elements in [n]
 * output is in [c], and will be a matrix of nChooseK x p integers */
uint64_t combinations(int** c, int n, int p, int x) {
    int j;
    uint64_t k = choose(n, p);
    for (j = 0; j < k; j++)
    {
        combination(&c[j], n, p, j);
    }
    return k;
}

// helper function for permute
void swap(int* v, const int i, const int j)
{
    int t;
    t = v[i];
    v[i] = v[j];
    v[j] = t;
}

// helper function for permute
void rotateLeft(int* v, const int start, const int n)
{
    int tmp = v[start];
    for (int i = start; i < n - 1; i++) {
        v[i] = v[i + 1];
    }
    v[n - 1] = tmp;
} // rotateLeft

// in-place recursive permutation generator.
// adapted from:
// https://stackoverflow.com/questions/3862191/permutation-generator-on-c
// provide a list of integers, a start value (typically 0) and size
// of the list, and a function pointer that performs an action on
// each evaluation.  The eval function will get the list, list size,
// and an optional user-supplied void-pointer to arbitrary data as arguments.
void permute(int* v, const int start, const int n,
    void (*eval_fcn)(int*, const int, void*), void* eval_params)
{
    eval_fcn(v, n, eval_params);
    if (start < n) {
        int i, j;
        for (i = n - 2; i >= start; i--) {
            for (j = i + 1; j < n; j++) {
                swap(v, i, j);
                permute(v, i + 1, n, eval_fcn, eval_params);
            } // for j
            rotateLeft(v, i, n);
        } // for i
    }
} // permute

// example eval function - to print each permutation
void print(const int* v, const int size, void* custom_data)
{
    if (v != 0) {
        for (int i = 0; i < size; i++) {
            printf("%4d", v[i]);
        }
        printf("\n");
    }
} // print


// recursive combination generator, with repetitions
// adapted from: https://www.geeksforgeeks.org/combinations-with-repetitions/
/* arr[]  ---> Input Array
  chosen[] ---> Temporary array to store indices of
                   current combination
   start & end ---> Starting and Ending indexes in arr[]
   r ---> Size of a combination to be printed */
void CombinationRepetitionUtil(int chosen[], int arr[],
    int index, int r, int start, int end, int* num)
{
    // Since index has become r, current combination is
    // ready to be printed, print
    if (index == r)
    {
        printf("%05d: ", (*num)++);
        for (int i = 0; i < r; i++)
            printf("%d ", arr[chosen[i]]);
        printf("\n");
        return;
    }

    // One by one choose all elements (without considering
    // the fact whether element is already chosen or not)
    // and recur
    for (int i = start; i <= end; i++)
    {
        chosen[index] = i;
        CombinationRepetitionUtil(chosen, arr, index + 1,
            r, i, end, num);
    }
    return;
}

// The main function that prints all combinations of size r
// in arr[] of size n with repetitions. This function mainly
// uses CombinationRepetitionUtil()
void CombinationRepetition(int arr[], int n, int r)
{
    // Allocate memory
    int* chosen = (int*)xmalloc((r + 1) * sizeof(int));
    int num;

    // Call the recursive function
    CombinationRepetitionUtil(chosen, arr, 0, r, 0, n - 1, &num);

    free(chosen);
}



