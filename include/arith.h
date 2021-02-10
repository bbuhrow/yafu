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

#include "yafu.h"
#include <sys/types.h>

#define LIMB_BLKSZ 10	
#define MAX_DIGITS 100
#define MP_RADIX 4294967296.0
#define LN2		0.69314718055994530942

// types of numbers
#define PRIME 0
#define PRP 1
#define COMPOSITE 2
#define UNKNOWN 3


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
__inline uint32 _trail_zcnt(uint32 x)
{
    uint32 pos;
    if (_BitScanForward(&pos, x))
        return pos;
    else
        return 32;
}
__inline uint64 _trail_zcnt64(uint64 x)
{
    uint64 pos;
    if (_BitScanForward64(&pos, x))
        return pos;
    else
        return 64;
}
__inline uint64 _lead_zcnt64(uint64 x)
{
    uint64 pos;
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
#include <immintrin.h>
#ifdef USE_BMI2
#define _lead_zcnt64 __lzcnt64
#define _trail_zcnt _tzcnt_u32
#define _trail_zcnt64 _tzcnt_u64
// MSVC pages say these are available, but using them
// leads to immediate crashes in msvc-190 (v160)
// https://docs.microsoft.com/en-us/cpp/intrinsics/x64-amd64-intrinsics-list?view=msvc-160
//#define _reset_lsb(x) _blsr_u32(x)
//#define _reset_lsb64(x) _blsr_u64(x)
#define _reset_lsb(x) ((x) &= ((x) - 1))
#define _reset_lsb64(x) ((x) &= ((x) - 1))

#else
__inline uint32 _trail_zcnt(uint32 x)
{
    uint32 pos;
    if (_BitScanForward(&pos, x))
        return pos;
    else
        return 32;
}
__inline uint64 _trail_zcnt64(uint64 x)
{
    uint64 pos;
    if (_BitScanForward64(&pos, x))
        return pos;
    else
        return 64;
}
__inline uint64 _lead_zcnt64(uint64 x)
{
    uint64 pos;
    if (_BitScanReverse64(&pos, x))
        return pos;
    else
        return 64;
}
#define _reset_lsb(x) ((x) &= ((x) - 1))
#define _reset_lsb64(x) ((x) &= ((x) - 1))
#endif

#else

__inline uint64 _lead_zcnt64(uint64 x)
{
    uint64 pos;
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

__inline uint32 _trail_zcnt(uint32 x)
{
    uint32 pos;
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

__inline uint64 _trail_zcnt64(uint64 x)
{
    uint64 pos;
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
void spAdd(fp_digit u, fp_digit v, fp_digit *sum, fp_digit *carry);
void spAdd3(fp_digit u, fp_digit v, fp_digit w, fp_digit *sum, fp_digit *carry);
void spSub3(fp_digit u, fp_digit v, fp_digit w, fp_digit *sub, fp_digit *borrow);
void spSub(fp_digit u, fp_digit v, fp_digit *sub, fp_digit *borrow);
void spMultiply(fp_digit u, fp_digit v, fp_digit *product, fp_digit *carry);
void spMulAdd(fp_digit u, fp_digit v, fp_digit w, fp_digit t, fp_digit *lower, fp_digit *carry);
void spMulMod(fp_digit u, fp_digit v, fp_digit m, fp_digit *w);
void spModExp(fp_digit a, fp_digit b, fp_digit m, fp_digit *u);
fp_digit spDivide(fp_digit *q, fp_digit *r, fp_digit u[2], fp_digit v);
fp_digit spBits(fp_digit n);
int bits64(uint64 n);
uint32 modinv_1(uint32 a, uint32 p);
uint32 modinv_1b(uint32 a, uint32 p);
uint32 modinv_1c(uint32 a, uint32 p);
uint64 spPRP2(uint64 p);
void ShanksTonelli_1(fp_digit a, fp_digit p, fp_digit *sq);
fp_digit spGCD(fp_digit x, fp_digit y);
uint64 spBinGCD(uint64 x, uint64 y);
uint64 spBinGCD_odd(uint64 u, uint64 v);
uint64 gcd64(uint64 x, uint64 y);
uint64 bingcd64(uint64 x, uint64 y);
void dblGCD(double x, double y, double* w);
int jacobi_1(fp_digit n, fp_digit p);
int ndigits_1(fp_digit n);

/********************* a few gmp-based utilities **********************/
double zlog(mpz_t x);
int llt(uint32 exp, int vflag);
int gmp_base10(mpz_t x);
int is_mpz_prp(mpz_t n, int num_witnesses);
uint64 mpz_get_64(mpz_t src);
void mpz_set_64(mpz_t dest, uint64 src);
void build_RSA(int bits, mpz_t n, gmp_randstate_t gmp_randstate);
void gordon(int bits, mpz_t n, gmp_randstate_t gmp_randstate);


#endif
