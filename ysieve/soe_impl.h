/*
MIT License

Copyright (c) 2021 Ben Buhrow

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
*/

#ifndef SOE_IMPL_H
#define SOE_IMPL_H


#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#if defined(_MSC_VER) && defined(__clang__)
#include <x86intrin.h>
#else
#include <immintrin.h>
#endif
#include "gmp.h"
#include "ytools.h"
#include "soe.h"

#define BITSINBYTE 8
#define MAXSIEVEPRIMECOUNT 100000000	//# primes less than ~2e9: limit of 2e9^2 = 4e18


#ifdef __INTEL_COMPILER
// leading and trailing zero count are ABM instructions
// that require haswell or later on Intel or ABM on AMD.
// The same basic functionality exists with the 
// bsf and bsr instructions that are standard x86, if
// those requirements are not met.
#if defined( USE_BMI2 ) || defined (TARGET_KNL) || defined( USE_AVX512F )
#define _reset_lsb(x) _blsr_u32(x)
#define _reset_lsb64(x) _blsr_u64(x)
#define _lead_zcnt64 __lzcnt64
#define _trail_zcnt _tzcnt_u32
#define _trail_zcnt64 _tzcnt_u64
#else
#define _reset_lsb(x) ((x) &= ((x) - 1))
#define _reset_lsb64(x) ((x) &= ((x) - 1))
__inline uint32_t _trail_zcnt(uint32_t x)
{
    uint32_t pos;
    if (_BitScanForward(&pos, x))
        return pos;
    else
        return 32;
}
__inline uint64_t _trail_zcnt64(uint64_t x)
{
    uint32_t pos;
    if (_BitScanForward64(&pos, x))
        return pos;
    else
        return 64;
}
__inline uint64_t _lead_zcnt64(uint64_t x)
{
    uint32_t pos;
    if (_BitScanReverse64(&pos, x))
        return pos;
    else
        return 64;
}
#endif
#elif defined(__GNUC__)
#if defined( USE_BMI2 ) || defined (TARGET_KNL) || defined( USE_AVX512F )
#define _reset_lsb(x) _blsr_u32(x)
#define _reset_lsb64(x) _blsr_u64(x)
#define _lead_zcnt64 __builtin_clzll
#define _trail_zcnt __builtin_ctzl
#define _trail_zcnt64 __builtin_ctzll
#else
#define _reset_lsb(x) ((x) &= ((x) - 1))
#define _reset_lsb64(x) ((x) &= ((x) - 1))
#define _lead_zcnt64 __builtin_clzll
#define _trail_zcnt __builtin_ctzl
#define _trail_zcnt64 __builtin_ctzll

#endif
#elif defined(_MSC_VER)

#include <intrin.h>

#ifdef USE_BMI2
#define _reset_lsb(x) _blsr_u32(x)
#define _reset_lsb64(x) _blsr_u64(x)
#else
#define _reset_lsb(x) ((x) &= ((x) - 1))
#define _reset_lsb64(x) ((x) &= ((x) - 1))
#endif

#ifdef __clang__
    __inline uint32_t _trail_zcnt(uint32_t x)
    {
        unsigned long pos;
        if (_BitScanForward(&pos, x))
            return (uint32_t)pos;
        else
            return 32;
    }
    __inline uint32_t _trail_zcnt64(uint64_t x)
    {
        unsigned long  pos;
        if (_BitScanForward64(&pos, x))
            return (uint32_t)pos;
        else
            return 64;
    }
    __inline uint32_t _lead_zcnt64(uint64_t x)
    {
        unsigned long  pos;
        if (_BitScanReverse64(&pos, x))
            return (uint32_t)pos;
        else
            return 64;
    }
#else
__inline uint32_t _trail_zcnt(uint32_t x)
{
    uint32_t pos;
    if (_BitScanForward(&pos, x))
        return pos;
    else
        return 32;
}
__inline uint32_t _trail_zcnt64(uint64_t x)
{
    uint32_t pos;
    if (_BitScanForward64(&pos, x))
        return pos;
    else
        return 64;
}
__inline uint32_t _lead_zcnt64(uint64_t x)
{
    uint32_t pos;
    if (_BitScanReverse64(&pos, x))
        return pos;
    else
        return 64;
}
#endif


//#endif

#else

__inline uint64_t _lead_zcnt64(uint64_t x)
{
    uint64_t pos;
    if (x)
    {
        pos = 0;
        for (pos = 0; ; pos++)
        {
            if (x & (1ULL << (63 - pos)))
                break;
        }
    }
    else
    {
#ifdef CHAR_BIT
        pos = CHAR_BIT * sizeof(x);
#else
        pos = 8 * sizeof(x);
#endif
    }
    return pos;
}

__inline uint32_t _trail_zcnt(uint32_t x)
{
    uint32_t pos;
    if (x)
    {
        x = (x ^ (x - 1)) >> 1;  // Set x's trailing 0s to 1s and zero rest
        for (pos = 0; x; pos++)
        {
            x >>= 1;
        }
    }
    else
    {
#ifdef CHAR_BIT
        pos = CHAR_BIT * sizeof(x);
#else
        pos = 8 * sizeof(x);
#endif
    }
    return pos;
}

__inline uint64_t _trail_zcnt64(uint64_t x)
{
    uint64_t pos;
    if (x)
    {
        x = (x ^ (x - 1)) >> 1;  // Set x's trailing 0s to 1s and zero rest
        for (pos = 0; x; pos++)
        {
            x >>= 1;
        }
    }
    else
    {
#ifdef CHAR_BIT
        pos = CHAR_BIT * sizeof(x);
#else
        pos = 8 * sizeof(x);
#endif
    }
    return pos;
}
#define _reset_lsb(x) ((x) &= ((x) - 1))
#define _reset_lsb64(x) ((x) &= ((x) - 1))

#endif

#ifdef USE_AVX2

#ifdef USE_AVX512F
extern ALIGNED_MEM uint64_t presieve_largemasks[16][173][8];
extern ALIGNED_MEM uint32_t presieve_steps[32];
extern ALIGNED_MEM uint32_t presieve_primes[32];
extern ALIGNED_MEM uint32_t presieve_p1[32];

#else
// for storage of presieving lists from prime index 24 to 40 (97 to 173 inclusive)
extern ALIGNED_MEM uint64_t presieve_largemasks[16][173][4];
extern ALIGNED_MEM uint32_t presieve_steps[32];
extern ALIGNED_MEM uint32_t presieve_primes[32];
extern ALIGNED_MEM uint32_t presieve_p1[32];
#endif

#endif

// thread ready sieving functions
void sieve_line(thread_soedata_t* thread_data);
void sieve_line_avx2_32k(thread_soedata_t* thread_data);
void sieve_line_avx2_64k(thread_soedata_t* thread_data);
void sieve_line_avx2_128k(thread_soedata_t* thread_data);
void sieve_line_avx2_256k(thread_soedata_t* thread_data);
void sieve_line_avx2_512k(thread_soedata_t* thread_data);
void sieve_line_avx512_32k(thread_soedata_t* thread_data);
void sieve_line_avx512_64k(thread_soedata_t* thread_data);
void sieve_line_avx512_128k(thread_soedata_t* thread_data);
void sieve_line_avx512_256k(thread_soedata_t* thread_data);
void sieve_line_avx512_512k(thread_soedata_t* thread_data);


uint64_t count_line(soe_staticdata_t* sdata, uint32_t current_line);
uint64_t count_twins(soe_staticdata_t* sdata, thread_soedata_t* thread_data);
void count_line_special(thread_soedata_t* thread_data);
uint32_t compute_32_bytes(soe_staticdata_t* sdata,
    uint32_t pcount, uint64_t* primes, uint64_t byte_offset);
uint64_t primes_from_lineflags(soe_staticdata_t* sdata, thread_soedata_t* thread_data,
    uint32_t start_count, uint64_t* primes);
void get_offsets(thread_soedata_t* thread_data);
void getRoots(soe_staticdata_t* sdata, thread_soedata_t* thread_data);

//void stop_soe_worker_thread(thread_soedata_t* t);
//void start_soe_worker_thread(thread_soedata_t* t);
//#if defined(WIN32) || defined(_WIN64)
//DWORD WINAPI soe_worker_thread_main(LPVOID thread_data);
//#else
//void* soe_worker_thread_main(void* thread_data);
//#endif

// routines for finding small numbers of primes; seed primes for main SOE
uint32_t tiny_soe(uint32_t limit, uint32_t* primes);

// top level sieving routines
uint64_t* GetPRIMESRange(soe_staticdata_t* sdata,
    mpz_t* offset, uint64_t lowlimit, uint64_t highlimit, uint64_t* num_p);
uint64_t spSOE(soe_staticdata_t* sdata, mpz_t* offset,
    uint64_t lowlimit, uint64_t* highlimit, int count, uint64_t* primes);

// misc and helper functions
uint64_t estimate_primes_in_range(uint64_t lowlimit, uint64_t highlimit);
uint64_t mpz_estimate_primes_in_range(mpz_t lowlimit, mpz_t highlimit);
void get_numclasses(uint64_t highlimit, uint64_t lowlimit, soe_staticdata_t* sdata);
int check_input(uint64_t highlimit, uint64_t lowlimit, uint32_t num_sp, uint32_t* sieve_p,
    soe_staticdata_t* sdata, mpz_t offset);
uint64_t init_sieve(soe_staticdata_t* sdata);
void set_bucket_depth(soe_staticdata_t* sdata);
uint64_t alloc_threaddata(soe_staticdata_t* sdata, thread_soedata_t* thread_data);
void do_soe_sieving(soe_staticdata_t* sdata, thread_soedata_t* thread_data, int count);
void finalize_sieve(soe_staticdata_t* sdata,
    thread_soedata_t* thread_data, int count, uint64_t* primes);
void trim_line(soe_staticdata_t* sdata, int current_line);

uint32_t modinv1(uint32_t a, uint32_t p);
uint32_t modinv2(uint32_t a, uint32_t p);
uint32_t modinv3(uint32_t a, uint32_t p);
uint64_t gcd_1(uint64_t x, uint64_t y);

void pre_sieve(soe_dynamicdata_t* ddata, soe_staticdata_t* sdata, uint8_t* flagblock);
void pre_sieve_avx2(soe_dynamicdata_t* ddata, soe_staticdata_t* sdata, uint8_t* flagblock);
void pre_sieve_avx512(soe_dynamicdata_t* ddata, soe_staticdata_t* sdata, uint8_t* flagblock);


uint32_t compute_8_bytes(soe_staticdata_t* sdata,
    uint32_t pcount, uint64_t* primes, uint64_t byte_offset);
uint32_t compute_8_bytes_bmi2(soe_staticdata_t* sdata,
    uint32_t pcount, uint64_t* primes, uint64_t byte_offset);

// declare the fat-binary function pointers	
extern uint32_t(*compute_8_bytes_ptr)(soe_staticdata_t*, uint32_t, uint64_t*, uint64_t);
extern void (*pre_sieve_ptr)(soe_dynamicdata_t*, soe_staticdata_t*, uint8_t*);
extern void(*sieve_line_ptr)(thread_soedata_t*);



#ifdef __cplusplus
}
#endif

#endif /* #ifndef SOE_IMPL_H */

