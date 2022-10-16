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

#ifndef ARITH_H
#define ARITH_H

#include <stdint.h>
#include <sys/types.h>
#include "gmp.h"

#define LIMB_BLKSZ 10	
#define MAX_DIGITS 100
#define MP_RADIX 4294967296.0
#define LN2		0.69314718055994530942

// types of numbers
#define PRIME 0
#define PRP 1
#define COMPOSITE 2
#define UNKNOWN 3


#if defined( __INTEL_COMPILER)
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
    uint64_t pos;
    if (_BitScanForward64(&pos, x))
        return pos;
    else
        return 64;
}
__inline uint64_t _lead_zcnt64(uint64_t x)
{
    uint64_t pos;
    if (_BitScanReverse64(&pos, x))
        return pos;
    else
        return 64;
}
#endif
#elif defined(__GNUC__) || defined(__INTEL_LLVM_COMPILER)
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
#include <immintrin.h>
//#ifdef USE_BMI2

// not safe to assume these just because this is defined.
// instead always use the inline functions below.

//#define _lead_zcnt64 __lzcnt64
//#define _trail_zcnt _tzcnt_u32
//#define _trail_zcnt64 _tzcnt_u64
//// MSVC pages say these are available, but using them
//// leads to immediate crashes in msvc-190 (v160)
//// https://docs.microsoft.com/en-us/cpp/intrinsics/x64-amd64-intrinsics-list?view=msvc-160
////#define _reset_lsb(x) _blsr_u32(x)
////#define _reset_lsb64(x) _blsr_u64(x)
//#define _reset_lsb(x) ((x) &= ((x) - 1))
//#define _reset_lsb64(x) ((x) &= ((x) - 1))
//
//#else
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
#define _reset_lsb(x) ((x) &= ((x) - 1))
#define _reset_lsb64(x) ((x) &= ((x) - 1))
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

// arbitrary precision arith routines
/********************* single precision arith **********************/
void spAdd(uint64_t u, uint64_t v, uint64_t *sum, uint64_t *carry);
void spAdd3(uint64_t u, uint64_t v, uint64_t w, uint64_t *sum, uint64_t *carry);
void spSub3(uint64_t u, uint64_t v, uint64_t w, uint64_t *sub, uint64_t *borrow);
void spSub(uint64_t u, uint64_t v, uint64_t *sub, uint64_t *borrow);
void spMultiply(uint64_t u, uint64_t v, uint64_t *product, uint64_t *carry);
void spMulAdd(uint64_t u, uint64_t v, uint64_t w, uint64_t t, uint64_t *lower, uint64_t *carry);
void spMulMod(uint64_t u, uint64_t v, uint64_t m, uint64_t *w);
void spModExp(uint64_t a, uint64_t b, uint64_t m, uint64_t *u);
uint64_t spDivide(uint64_t *q, uint64_t *r, uint64_t u[2], uint64_t v);
uint64_t spBits(uint64_t n);
int bits64(uint64_t n);
uint32_t modinv_1(uint32_t a, uint32_t p);
uint32_t modinv_1b(uint32_t a, uint32_t p);
uint32_t modinv_1c(uint32_t a, uint32_t p);
uint64_t spPRP2(uint64_t p);
void ShanksTonelli_1(uint64_t a, uint64_t p, uint64_t *sq);
uint64_t spGCD(uint64_t x, uint64_t y);
uint64_t spBinGCD(uint64_t x, uint64_t y);
uint64_t spBinGCD_odd(uint64_t u, uint64_t v);
uint64_t gcd64(uint64_t x, uint64_t y);
uint64_t bingcd64(uint64_t x, uint64_t y);
void dblGCD(double x, double y, double* w);
int jacobi_1(uint64_t n, uint64_t p);
int ndigits_1(uint64_t n);

/********************* a few gmp-based utilities **********************/
double zlog(mpz_t x);
int llt(uint32_t exp, int vflag);
int gmp_base10(mpz_t x);
int is_mpz_prp(mpz_t n, int num_witnesses);
uint64_t mpz_get_64(mpz_t src);
void mpz_set_64(mpz_t dest, uint64_t src);
void build_RSA(int bits, mpz_t n, gmp_randstate_t gmp_randstate);
void gordon(int bits, mpz_t n, gmp_randstate_t gmp_randstate);


#endif
