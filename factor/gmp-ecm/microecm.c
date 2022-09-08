/*
Copyright (c) 2014, Ben Buhrow and (c) 2022, Jeff Hurchalla.
All rights reserved.
This software is provided "as is," without warranty of any kind,
express or implied.  In no event shall the author or contributors
be held liable for any damages arising in any way from the use of
this software.
The contents of this file are DUAL-LICENSED.  You may modify and/or
redistribute this software according to the terms of one of the
following two licenses (at your option):


LICENSE 1 ("FreeBSD")
---------------------
Copyright (c) 2014, Ben Buhrow and (c) 2022, Jeff Hurchalla.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.


LICENSE 2 (MPL 2.0)
-------------------
Copyright (c) 2014, Ben Buhrow and (c) 2022, Jeff Hurchalla.
All rights reserved.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
*/


#include "microecm.h"
#include "monty.h"

#include <stddef.h>
#include <stdint.h>
#if defined(_MSC_VER)
#  ifndef _WIN64
#    error "64 bit compilation mode is required for MSVC"
#  endif
#  include <immintrin.h>
#  include <intrin.h>
#endif



// If we define (or predefine) MICRO_ECM_EXPECT_ARBITRARY_SIZE_FACTORS, then
// the micro-ECM factoring will first try a few curves that target mid-size to
// semi-large factors for the given number to factor.  If no factor is found,
// then factoring will proceed in the normal way - it will try curves that are
// always targetted toward finding large factors.

//#define MICRO_ECM_EXPECT_ARBITRARY_SIZE_FACTORS

// Using the inline asm in this file can increase performance by ~20-25%
// (surprisingly).  Hence these macros are defined by default.
#if defined(__x86_64__) || defined(_M_X64)
#  define MICRO_ECM_ALT_MULREDC_USE_INLINE_ASM_X86
#endif

// We rarely will ever want to use debugging printfs
//#define MICRO_ECM_VERBOSE_PRINTF




#ifdef _MSC_VER
#  define MICRO_ECM_FORCE_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
#  define MICRO_ECM_FORCE_INLINE inline __attribute__((always_inline))
#else
#  define MICRO_ECM_FORCE_INLINE __inline
#endif


#ifdef MICRO_ECM_VERBOSE_PRINTF
#  include <stdio.h>
// these globals will cause data races, if we run multithreaded
    uint32_t stg1Doub;
    uint32_t stg1Add;
    uint32_t ptadds;
#endif


typedef struct
{
    uint64_t X;
    uint64_t Z;
} uecm_pt;


static const uint32_t map[61] = {
    0, 1, 2, 0, 0, 0, 0, 3, 0, 0,
    0, 4, 0, 5, 0, 0, 0, 6, 0, 7,
    0, 0, 0, 8, 0, 0, 0, 0, 0, 9,
    0, 10, 0, 0, 0, 0, 0, 11, 0, 0,
    0, 12, 0, 13, 0, 0, 0, 14, 0, 15,
    0, 0, 0, 16, 0, 0, 0, 0, 0, 17,
    18 };

static const double INV_2_POW_32 = 1.0 / (double)((uint64_t)(1) << 32);


static uint32_t uecm_lcg_rand_32B(uint32_t lower, uint32_t upper, uint64_t *ploc_lcg)
{
    *ploc_lcg = 6364136223846793005ULL * (*ploc_lcg) + 1442695040888963407ULL;
    return lower + (uint32_t)(
        (double)(upper - lower) * (double)((*ploc_lcg) >> 32) * INV_2_POW_32);
}



/* --- The following two functions are written by Jeff Hurchalla, Copyright 2022 --- */

// for this algorithm, see https://jeffhurchalla.com/2022/04/28/montgomery-redc-using-the-positive-inverse-mod-r/
MICRO_ECM_FORCE_INLINE static uint64_t uecm_mulredc_alt(uint64_t x, uint64_t y, uint64_t N, uint64_t invN)
{
#if defined(_MSC_VER)
    uint64_t T_hi;
    uint64_t T_lo = _umul128(x, y, &T_hi);
    uint64_t m = T_lo * invN;
    uint64_t mN_hi = __umulh(m, N);
#else
    __uint128_t prod = (__uint128_t)x * y;
    uint64_t T_hi = (uint64_t)(prod >> 64);
    uint64_t T_lo = (uint64_t)(prod);
    uint64_t m = T_lo * invN;
    __uint128_t mN = (__uint128_t)m * N;
    uint64_t mN_hi = (uint64_t)(mN >> 64);
#endif
    uint64_t tmp = T_hi + N;
#if defined(MICRO_ECM_ALT_MULREDC_USE_INLINE_ASM_X86) && !defined(_MSC_VER)
    __asm__ (
        "subq %[mN_hi], %[tmp] \n\t"    /* tmp = T_hi + N - mN_hi */
        "subq %[mN_hi], %[T_hi] \n\t"   /* T_hi = T_hi - mN_hi */
        "cmovaeq %[T_hi], %[tmp] \n\t"  /* tmp = (T_hi >= mN_hi) ? T_hi : tmp */
        : [tmp]"+&r"(tmp), [T_hi]"+&r"(T_hi)
        : [mN_hi]"r"(mN_hi)
        : "cc");
    uint64_t result = tmp;
#else
    tmp = tmp - mN_hi;
    uint64_t result = T_hi - mN_hi;
    result = (T_hi < mN_hi) ? tmp : result;
#endif
   return result;
}

// for this algorithm, see https://jeffhurchalla.com/2022/04/25/a-faster-multiplicative-inverse-mod-a-power-of-2/
static uint64_t uecm_multiplicative_inverse(uint64_t a)
{
//    assert(a%2 == 1);  // the inverse (mod 2<<64) only exists for odd values
    uint64_t x0 = (3*a)^2;
    uint64_t y = 1 - a*x0;
    uint64_t x1 = x0*(1 + y);
    y *= y;
    uint64_t x2 = x1*(1 + y);
    y *= y;
    uint64_t x3 = x2*(1 + y);
    y *= y;
    uint64_t x4 = x3*(1 + y);
    return x4;
}

/* --- end Hurchalla functions --- */





// full strength mul/sqr redc
MICRO_ECM_FORCE_INLINE static uint64_t uecm_mulredc(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
{
    return uecm_mulredc_alt(x, y, n, 0 - nhat);
}
MICRO_ECM_FORCE_INLINE static uint64_t uecm_sqrredc(uint64_t x, uint64_t n, uint64_t nhat)
{
    return uecm_mulredc_alt(x, x, n, 0 - nhat);
}





MICRO_ECM_FORCE_INLINE static void uecm_uadd(uint64_t rho, uint64_t n, const uecm_pt P1, const uecm_pt P2,
    const uecm_pt Pin, uecm_pt *Pout)
{
    // compute:
    //x+ = z- * [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
    //z+ = x- * [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
    // where:
    //x- = original x
    //z- = original z
    // given the sums and differences of the original points
    uint64_t diff1 = submod(P1.X, P1.Z, n);
    uint64_t sum1 = addmod(P1.X, P1.Z, n);
    uint64_t diff2 = submod(P2.X, P2.Z, n);
    uint64_t sum2 = addmod(P2.X, P2.Z, n);

    uint64_t tt1 = uecm_mulredc(diff1, sum2, n, rho); //U
    uint64_t tt2 = uecm_mulredc(sum1, diff2, n, rho); //V

    uint64_t tt3 = addmod(tt1, tt2, n);
    uint64_t tt4 = submod(tt1, tt2, n);
    tt1 = uecm_sqrredc(tt3, n, rho);   //(U + V)^2
    tt2 = uecm_sqrredc(tt4, n, rho);   //(U - V)^2

    uint64_t tmpx = uecm_mulredc(tt1, Pin.Z, n, rho);     //Z * (U + V)^2
    uint64_t tmpz = uecm_mulredc(tt2, Pin.X, n, rho);     //x * (U - V)^2
    Pout->X = tmpx;
    Pout->Z = tmpz;

#ifdef MICRO_ECM_VERBOSE_PRINTF
    stg1Add++;
#endif
    return;
}

MICRO_ECM_FORCE_INLINE static void uecm_udup(uint64_t s, uint64_t rho, uint64_t n,
    uint64_t insum, uint64_t indiff, uecm_pt *P)
{
    uint64_t tt1 = uecm_sqrredc(indiff, n, rho);          // U=(x1 - z1)^2
    uint64_t tt2 = uecm_sqrredc(insum, n, rho);           // V=(x1 + z1)^2
    P->X = uecm_mulredc(tt1, tt2, n, rho);         // x=U*V

    uint64_t tt3 = submod(tt2, tt1, n);          // w = V-U
    tt2 = uecm_mulredc(tt3, s, n, rho);      // w = (A+2)/4 * w
    tt2 = addmod(tt2, tt1, n);          // w = w + U
    P->Z = uecm_mulredc(tt2, tt3, n, rho);         // Z = w*(V-U)
#ifdef MICRO_ECM_VERBOSE_PRINTF
    stg1Doub++;
#endif
    return;
}


static void uecm_uprac70(uint64_t rho, uint64_t n, uecm_pt *P, uint64_t s)
{
    uint64_t s1, s2, d1, d2;
    uint64_t swp;
    int i;
    static const uint8_t steps[116] = {
        0,6,0,6,0,6,0,4,6,0,4,6,0,4,4,6,
        0,4,4,6,0,5,4,6,0,3,3,4,6,0,3,5,
        4,6,0,3,4,3,4,6,0,5,5,4,6,0,5,3,
        3,4,6,0,3,3,4,3,4,6,0,5,3,3,3,3,
        3,3,3,3,4,3,3,4,6,0,5,4,3,3,4,6,
        0,3,4,3,5,4,6,0,5,3,3,3,4,6,0,5,
        4,3,5,4,6,0,5,5,3,3,4,6,0,4,3,3,
        3,5,4,6 };

    uecm_pt pt1, pt2, pt3;
    for (i = 0; i < 116; i++)
    {
        if (steps[i] == 0)
        {
            pt1.X = pt2.X = pt3.X = P->X;
            pt1.Z = pt2.Z = pt3.Z = P->Z;

            d1 = submod(pt1.X, pt1.Z, n);
            s1 = addmod(pt1.X, pt1.Z, n);
            uecm_udup(s, rho, n, s1, d1, &pt1);
        }
        else if (steps[i] == 3)
        {
            // integrate step 4 followed by swap(1,2)
            uecm_pt pt4;
            uecm_uadd(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

            swp = pt1.X;
            pt1.X = pt4.X;
            pt4.X = pt3.X;
            pt3.X = pt2.X;
            pt2.X = swp;
            swp = pt1.Z;
            pt1.Z = pt4.Z;
            pt4.Z = pt3.Z;
            pt3.Z = pt2.Z;
            pt2.Z = swp;
        }
        else if (steps[i] == 4)
        {
            uecm_pt pt4;
            uecm_uadd(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

            swp = pt2.X;
            pt2.X = pt4.X;
            pt4.X = pt3.X;
            pt3.X = swp;
            swp = pt2.Z;
            pt2.Z = pt4.Z;
            pt4.Z = pt3.Z;
            pt3.Z = swp;
        }
        else if (steps[i] == 5)
        {
            d2 = submod(pt1.X, pt1.Z, n);
            s2 = addmod(pt1.X, pt1.Z, n);

            uecm_uadd(rho, n, pt2, pt1, pt3, &pt2);        // B = B + A (C)
            uecm_udup(s, rho, n, s2, d2, &pt1);        // A = 2A
        }
        else if (steps[i] == 6)
        {
            uecm_uadd(rho, n, pt1, pt2, pt3, P);     // A = A + B (C)
        }
    }
    return;
}

static void uecm_uprac85(uint64_t rho, uint64_t n, uecm_pt *P, uint64_t s)
{
    uint64_t s1, s2, d1, d2;
    uint64_t swp;
    int i;
    static const uint8_t steps[146] = {
        0,6,0,6,0,6,0,6,0,4,
        6,0,4,6,0,4,4,6,0,4,
        4,6,0,5,4,6,0,3,3,4,
        6,0,3,5,4,6,0,3,4,3,
        4,6,0,5,5,4,6,0,5,3,
        3,4,6,0,3,3,4,3,4,6,
        0,4,3,4,3,5,3,3,3,3,
        3,3,3,3,4,6,0,3,3,3,
        3,3,3,3,3,3,4,3,4,3,
        4,6,0,3,4,3,5,4,6,0,
        5,3,3,3,4,6,0,5,4,3,
        5,4,6,0,4,3,3,3,5,4,
        6,0,4,3,5,3,3,4,6,0,
        3,3,3,3,5,4,6,0,3,3,
        3,4,3,3,4,6 };

    uecm_pt pt1, pt2, pt3;
    for (i = 0; i < 146; i++)
    {
        if (steps[i] == 0)
        {
            pt1.X = pt2.X = pt3.X = P->X;
            pt1.Z = pt2.Z = pt3.Z = P->Z;

            d1 = submod(pt1.X, pt1.Z, n);
            s1 = addmod(pt1.X, pt1.Z, n);
            uecm_udup(s, rho, n, s1, d1, &pt1);
        }
        else if (steps[i] == 3)
        {
            // integrate step 4 followed by swap(1,2)
            uecm_pt pt4;
            uecm_uadd(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

            swp = pt1.X;
            pt1.X = pt4.X;
            pt4.X = pt3.X;
            pt3.X = pt2.X;
            pt2.X = swp;
            swp = pt1.Z;
            pt1.Z = pt4.Z;
            pt4.Z = pt3.Z;
            pt3.Z = pt2.Z;
            pt2.Z = swp;
        }
        else if (steps[i] == 4)
        {
            uecm_pt pt4;
            uecm_uadd(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

            swp = pt2.X;
            pt2.X = pt4.X;
            pt4.X = pt3.X;
            pt3.X = swp;
            swp = pt2.Z;
            pt2.Z = pt4.Z;
            pt4.Z = pt3.Z;
            pt3.Z = swp;
        }
        else if (steps[i] == 5)
        {
            d2 = submod(pt1.X, pt1.Z, n);
            s2 = addmod(pt1.X, pt1.Z, n);

            uecm_uadd(rho, n, pt2, pt1, pt3, &pt2);        // B = B + A (C)
            uecm_udup(s, rho, n, s2, d2, &pt1);        // A = 2A
        }
        else if (steps[i] == 6)
        {
            uecm_uadd(rho, n, pt1, pt2, pt3, P);     // A = A + B (C)
        }
    }
    return;
}

static void uecm_uprac(uint64_t rho, uint64_t n, uecm_pt *P, uint64_t c, double v, uint64_t s)
{
    uint64_t d, e, r;
    int i;
    uint64_t s1, s2, d1, d2;
    uint64_t swp;

    // we require c != 0
    int shift = _trail_zcnt64(c);
    c = c >> shift;

    d = c;
    r = (uint64_t)((double)d * v + 0.5);

    d = c - r;
    e = 2 * r - c;

    uecm_pt pt1, pt2, pt3;

    // the first one is always a doubling
    // point1 is [1]P
    pt1.X = pt2.X = pt3.X = P->X;
    pt1.Z = pt2.Z = pt3.Z = P->Z;

    d1 = submod(pt1.X, pt1.Z, n);
    s1 = addmod(pt1.X, pt1.Z, n);

    // point2 is [2]P
    uecm_udup(s, rho, n, s1, d1, &pt1);

    while (d != e)
    {
        if (d < e)
        {
            r = d;
            d = e;
            e = r;
            swp = pt1.X;
            pt1.X = pt2.X;
            pt2.X = swp;
            swp = pt1.Z;
            pt1.Z = pt2.Z;
            pt2.Z = swp;
        }
        if (d - e <= e / 4 && ((d + e) % 3) == 0)
        {
            d = (2 * d - e) / 3;
            e = (e - d) / 2;

            uecm_pt pt4;
            uecm_uadd(rho, n, pt1, pt2, pt3, &pt4); // T = A + B (C)
            uecm_pt pt5;
            uecm_uadd(rho, n, pt4, pt1, pt2, &pt5); // T2 = T + A (B)
            uecm_uadd(rho, n, pt2, pt4, pt1, &pt2); // B = B + T (A)

            swp = pt1.X;
            pt1.X = pt5.X;
            pt5.X = swp;
            swp = pt1.Z;
            pt1.Z = pt5.Z;
            pt5.Z = swp;
        }
        else if (d - e <= e / 4 && (d - e) % 6 == 0)
        {
            d = (d - e) / 2;

            d1 = submod(pt1.X, pt1.Z, n);
            s1 = addmod(pt1.X, pt1.Z, n);

            uecm_uadd(rho, n, pt1, pt2, pt3, &pt2);        // B = A + B (C)
            uecm_udup(s, rho, n, s1, d1, &pt1);        // A = 2A
        }
        else if ((d + 3) / 4 <= e)
        {
            d -= e;

            uecm_pt pt4;
            uecm_uadd(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

            swp = pt2.X;
            pt2.X = pt4.X;
            pt4.X = pt3.X;
            pt3.X = swp;
            swp = pt2.Z;
            pt2.Z = pt4.Z;
            pt4.Z = pt3.Z;
            pt3.Z = swp;
        }
        else if ((d + e) % 2 == 0)
        {
            d = (d - e) / 2;

            d2 = submod(pt1.X, pt1.Z, n);
            s2 = addmod(pt1.X, pt1.Z, n);

            uecm_uadd(rho, n, pt2, pt1, pt3, &pt2);        // B = B + A (C)
            uecm_udup(s, rho, n, s2, d2, &pt1);        // A = 2A
        }
        else if (d % 2 == 0)
        {
            d /= 2;

            d2 = submod(pt1.X, pt1.Z, n);
            s2 = addmod(pt1.X, pt1.Z, n);

            uecm_uadd(rho, n, pt3, pt1, pt2, &pt3);        // C = C + A (B)
            uecm_udup(s, rho, n, s2, d2, &pt1);        // A = 2A
        }
        else if (d % 3 == 0)
        {
            d = d / 3 - e;

            d1 = submod(pt1.X, pt1.Z, n);
            s1 = addmod(pt1.X, pt1.Z, n);

            uecm_pt pt4;
            uecm_udup(s, rho, n, s1, d1, &pt4);        // T = 2A
            uecm_pt pt5;
            uecm_uadd(rho, n, pt1, pt2, pt3, &pt5);        // T2 = A + B (C)
            uecm_uadd(rho, n, pt4, pt1, pt1, &pt1);        // A = T + A (A)
            uecm_uadd(rho, n, pt4, pt5, pt3, &pt4);        // T = T + T2 (C)

            swp = pt3.X;
            pt3.X = pt2.X;
            pt2.X = pt4.X;
            pt4.X = swp;
            swp = pt3.Z;
            pt3.Z = pt2.Z;
            pt2.Z = pt4.Z;
            pt4.Z = swp;

        }
        else if ((d + e) % 3 == 0)
        {
            d = (d - 2 * e) / 3;

            uecm_pt pt4;
            uecm_uadd(rho, n, pt1, pt2, pt3, &pt4);        // T = A + B (C)


            d2 = submod(pt1.X, pt1.Z, n);
            s2 = addmod(pt1.X, pt1.Z, n);
            uecm_uadd(rho, n, pt4, pt1, pt2, &pt2);        // B = T + A (B)
            uecm_udup(s, rho, n, s2, d2, &pt4);        // T = 2A
            uecm_uadd(rho, n, pt1, pt4, pt1, &pt1);        // A = A + T (A) = 3A
        }
        else if ((d - e) % 3 == 0)
        {
            d = (d - e) / 3;

            uecm_pt pt4;
            uecm_uadd(rho, n, pt1, pt2, pt3, &pt4);        // T = A + B (C)

            d2 = submod(pt1.X, pt1.Z, n);
            s2 = addmod(pt1.X, pt1.Z, n);
            uecm_uadd(rho, n, pt3, pt1, pt2, &pt3);        // C = C + A (B)

            swp = pt2.X;
            pt2.X = pt4.X;
            pt4.X = swp;
            swp = pt2.Z;
            pt2.Z = pt4.Z;
            pt4.Z = swp;

            uecm_udup(s, rho, n, s2, d2, &pt4);        // T = 2A
            uecm_uadd(rho, n, pt1, pt4, pt1, &pt1);        // A = A + T (A) = 3A
        }
        else
        {
            e /= 2;

            d2 = submod(pt2.X, pt2.Z, n);
            s2 = addmod(pt2.X, pt2.Z, n);

            uecm_uadd(rho, n, pt3, pt2, pt1, &pt3);        // C = C + B (A)
            uecm_udup(s, rho, n, s2, d2, &pt2);        // B = 2B
        }
    }
    uecm_uadd(rho, n, pt1, pt2, pt3, P);     // A = A + B (C)

    for (i = 0; i < shift; i++)
    {
        d1 = submod(P->X, P->Z, n);
        s1 = addmod(P->X, P->Z, n);
        uecm_udup(s, rho, n, s1, d1, P);     // P = 2P
    }
    return;
}


// jeff: "likely_gcd" is probably always the correct gcd, but I didn't add this
// parameter by using any proof; it's conceivable it might be wrong sometimes.
// See comments within the function.
inline static uint64_t uecm_modinv_64(uint64_t a, uint64_t p, uint64_t* plikely_gcd)
{
    /* thanks to the folks at www.mersenneforum.org */

    uint64_t ps1, ps2, parity, dividend, divisor, rem, q, t;

    q = 1;
    rem = a;
    dividend = p;
    divisor = a;
    ps1 = 1;
    ps2 = 0;
    parity = 0;

    while (divisor > 1) {
        rem = dividend - divisor;
        t = rem - divisor;
        if (rem >= divisor) {
            q += ps1; rem = t; t -= divisor;
            if (rem >= divisor) {
                q += ps1; rem = t; t -= divisor;
                if (rem >= divisor) {
                    q += ps1; rem = t; t -= divisor;
                    if (rem >= divisor) {
                        q += ps1; rem = t; t -= divisor;
                        if (rem >= divisor) {
                            q += ps1; rem = t; t -= divisor;
                            if (rem >= divisor) {
                                q += ps1; rem = t; t -= divisor;
                                if (rem >= divisor) {
                                    q += ps1; rem = t; t -= divisor;
                                    if (rem >= divisor) {
                                        q += ps1; rem = t;
                                        if (rem >= divisor) {
                                            q = dividend / divisor;
                                            rem = dividend % divisor;
                                            q *= ps1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        q += ps2;
        parity = ~parity;
        dividend = divisor;
        divisor = rem;
        ps2 = ps1;
        ps1 = q;
    }

    // jeff: added "likely_gcd".  this function seems to be a variant on the
    // extended euclidean algorithm, and thus the gcd likely equals the dividend
    // as calculated below.  However I'm doing this by analogy and educated
    // guess, not by proof.  It appears to work in all tests, so I suspect it is
    // correct, but it's possible this could be wrong in some cases.
    if (divisor == 1)
        dividend = divisor;
    *plikely_gcd = dividend;


    if (parity == 0)
        return ps1;
    else
        return p - ps1;
}


static uint64_t uecm_build(uecm_pt *P, uint64_t rho, uint64_t n, uint64_t *ploc_lcg, uint64_t* ps, uint64_t five, uint64_t Rsqr)
{
    uint64_t t1, t2, t3, t4;
    uint64_t u, v;

    uint32_t sigma = uecm_lcg_rand_32B(7, (uint32_t)-1, ploc_lcg);

    u = uecm_mulredc((uint64_t)sigma, Rsqr, n, rho);  // to_monty(sigma)

    //printf("sigma = %" PRIu64 ", u = %" PRIu64 ", n = %" PRIu64 "\n", sigma, u, n);

    v = addmod(u, u, n);
    v = addmod(v, v, n);            // 4*sigma

    //printf("v = %" PRIu64 "\n", v);

    u = uecm_sqrredc(u, n, rho);
    t1 = five;

    //printf("monty(5) = %" PRIu64 "\n", t1);

    u = submod(u, t1, n);           // sigma^2 - 5

    //printf("u = %" PRIu64 "\n", u);

    t1 = uecm_mulredc(u, u, n, rho);
    uint64_t tmpx = uecm_mulredc(t1, u, n, rho);  // u^3

    uint64_t v2 = addmod(v, v, n);             // 2*v
    uint64_t v4 = addmod(v2, v2, n);           // 4*v
    uint64_t v8 = addmod(v4, v4, n);           // 8*v
    uint64_t v16 = addmod(v8, v8, n);          // 16*v
    uint64_t t5 = uecm_mulredc(v16, tmpx, n, rho);    // 16*u^3*v

    t1 = uecm_mulredc(v, v, n, rho);
    uint64_t tmpz = uecm_mulredc(t1, v, n, rho);  // v^3

    //compute parameter A
    t1 = submod(v, u, n);           // (v - u)
    t2 = uecm_sqrredc(t1, n, rho);
    t4 = uecm_mulredc(t2, t1, n, rho);   // (v - u)^3

    t1 = addmod(u, u, n);           // 2u
    t2 = addmod(u, v, n);           // u + v
    t3 = addmod(t1, t2, n);         // 3u + v

    t1 = uecm_mulredc(t3, t4, n, rho);   // a = (v-u)^3 * (3u + v)

    // u holds the denom (jeff note: isn't it t5 that has the denom?)
    // t1 holds the numer
    // accomplish the division by multiplying by the modular inverse
    t2 = 1;
    t5 = uecm_mulredc(t5, t2, n, rho);   // take t5 out of monty rep

    uint64_t likely_gcd;
    t3 = uecm_modinv_64(t5, n, &likely_gcd);

    t3 = uecm_mulredc(t3, Rsqr, n, rho); // to_monty(t3)
    *ps = uecm_mulredc(t3, t1, n, rho);

    P->X = tmpx;
    P->Z = tmpz;

    return likely_gcd;
}


static int uecm_check_factor(uint64_t Z, uint64_t n, uint64_t* f)
{
    int status = 0;
    *f = bingcd64(n, Z);

    if (*f > 1)
    {
        if (*f == n)
        {
            *f = 0;
            status = 0;
        }
        else
        {
            status = 1;
        }
    }
    return status;
}


static void uecm_stage1(uint64_t rho, uint64_t n, uecm_pt *P, uint64_t stg1, uint64_t s)
{
    uint64_t q;

    // handle the only even case
    q = 2;
    while (q < stg1 * 4)  // jeff: multiplying by 4 improves perf ~1%
    {
        uint64_t diff1 = submod(P->X, P->Z, n);
        uint64_t sum1 = addmod(P->X, P->Z, n);
        uecm_udup(s, rho, n, sum1, diff1, P);
        q *= 2;
    }

    if (stg1 == 47)
    {
        // jeff: improved perf slightly by using one more uprac for 3,
        // and removing uprac for 47.
        uecm_uprac(rho, n, P, 3, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 3, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 3, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 3, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 5, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 5, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 7, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 11, 0.580178728295464130, s);
        uecm_uprac(rho, n, P, 13, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 17, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 19, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 23, 0.522786351415446049, s);
        uecm_uprac(rho, n, P, 29, 0.548409048446403258, s);
        uecm_uprac(rho, n, P, 31, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 37, 0.580178728295464130, s);
        uecm_uprac(rho, n, P, 41, 0.548409048446403258, s);
        uecm_uprac(rho, n, P, 43, 0.618033988749894903, s);
//        uecm_uprac(rho, n, P, 47, 0.548409048446403258, s);
    }
    else if (stg1 == 59)
    {   // jeff: probably stg1 of 59 would benefit from similar changes
        // as stg1 of 47 above, but I didn't bother. Stg1 of 59 seems to
        // always perform worse than stg1 of 47, so there doesn't seem
        // to be any reason to ever use stg1 of 59.
        uecm_uprac(rho, n, P, 3, 0.61803398874989485, s);
        uecm_uprac(rho, n, P, 3, 0.61803398874989485, s);
        uecm_uprac(rho, n, P, 3, 0.61803398874989485, s);
        uecm_uprac(rho, n, P, 5, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 5, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 7, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 7, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 11, 0.580178728295464130, s);
        uecm_uprac(rho, n, P, 13, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 17, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 19, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 23, 0.522786351415446049, s);
        uecm_uprac(rho, n, P, 29, 0.548409048446403258, s);
        uecm_uprac(rho, n, P, 31, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 1961, 0.552936068843375, s);   // 37 * 53
        uecm_uprac(rho, n, P, 41, 0.548409048446403258, s);
        uecm_uprac(rho, n, P, 43, 0.618033988749894903, s);
        uecm_uprac(rho, n, P, 47, 0.548409048446403258, s);
        uecm_uprac(rho, n, P, 59, 0.548409048446403258, s);
    }
    else if (stg1 == 70)
    {
        // call prac with best ratio found in deep search.
        // some composites are cheaper than their
        // constituent primes.
        uecm_uprac70(rho, n, P, s);
    }
    else // if (stg1 >= 85)
    {
        uecm_uprac85(rho, n, P, s);

        if (stg1 == 85)
        {
            uecm_uprac(rho, n, P, 61, 0.522786351415446049, s);
        }
        else
        {
            uecm_uprac(rho, n, P, 5, 0.618033988749894903, s);
            uecm_uprac(rho, n, P, 11, 0.580178728295464130, s);
//            uecm_uprac(rho, n, P, 61, 0.522786351415446049, s);
            uecm_uprac(rho, n, P, 89, 0.618033988749894903, s);
            uecm_uprac(rho, n, P, 97, 0.723606797749978936, s);
            uecm_uprac(rho, n, P, 101, 0.556250337855490828, s);
            uecm_uprac(rho, n, P, 107, 0.580178728295464130, s);
            uecm_uprac(rho, n, P, 109, 0.548409048446403258, s);
            uecm_uprac(rho, n, P, 113, 0.618033988749894903, s);

            if (stg1 == 125)
            {
                // jeff: moved 61 to here
                uecm_uprac(rho, n, P, 61, 0.522786351415446049, s);
                uecm_uprac(rho, n, P, 103, 0.632839806088706269, s);
            }
            else
            {
                uecm_uprac(rho, n, P, 7747, 0.552188778811121, s); // 61 x 127
                uecm_uprac(rho, n, P, 131, 0.618033988749894903, s);
                uecm_uprac(rho, n, P, 14111, 0.632839806088706, s);  // 103 x 137
                uecm_uprac(rho, n, P, 20989, 0.620181980807415, s);  // 139 x 151
                uecm_uprac(rho, n, P, 157, 0.640157392785047019, s);
                uecm_uprac(rho, n, P, 163, 0.551390822543526449, s);

                if (stg1 == 165)
                {
                    uecm_uprac(rho, n, P, 149, 0.580178728295464130, s);
                }
                else
                {
                    uecm_uprac(rho, n, P, 13, 0.618033988749894903, s);
                    uecm_uprac(rho, n, P, 167, 0.580178728295464130, s);
                    uecm_uprac(rho, n, P, 173, 0.612429949509495031, s);
                    uecm_uprac(rho, n, P, 179, 0.618033988749894903, s);
                    uecm_uprac(rho, n, P, 181, 0.551390822543526449, s);
                    uecm_uprac(rho, n, P, 191, 0.618033988749894903, s);
                    uecm_uprac(rho, n, P, 193, 0.618033988749894903, s);
                    uecm_uprac(rho, n, P, 29353, 0.580178728295464, s);  // 149 x 197
                    uecm_uprac(rho, n, P, 199, 0.551390822543526449, s);
                }
            }
        }
    }
    return;
}



// pre-paired sequences for various B1 and B2 = 25*B1
static const int numb1_70 = 186;
static const uint8_t b1_70[186] = {
    53,49,47,43,41,37,23,19,13,11,1,7,17,29,31,0,59,47,43,41,37,31,29,19,13,7,1,11,23,0,59,53,43,41,37,
    31,23,17,11,7,1,19,29,49,0,53,49,47,43,31,23,19,11,7,1,13,37,59,0,59,53,43,37,31,29,23,17,13,11,1,47,
    0,59,49,41,31,23,17,11,7,1,19,37,47,0,59,49,47,43,41,31,17,13,11,7,37,0,53,49,43,37,23,19,13,7,1,29,
    31,41,59,0,59,49,47,41,23,19,17,13,7,1,43,53,0,59,49,43,37,29,17,13,7,1,19,47,53,0,59,53,49,47,43,31,
    29,23,11,17,0,47,43,41,37,31,23,19,17,11,1,13,29,53,0,59,47,41,37,31,23,19,11,7,17,29,0,53,47,43,41,
    17,13,11,1,23,31,37,49 };

static const int numb1_85 = 225;
static const uint8_t b1_85[225] = {
    1,53,49,47,43,41,37,23,19,13,11,1,7,17,29,31,0,59,47,43,41,37,31,29,19,13,7,1,11,23,0,59,53,43,41,37,
    31,23,17,11,7,1,19,29,49,0,53,49,47,43,31,23,19,11,7,1,13,37,59,0,59,53,43,37,31,29,23,17,13,11,1,47,
    0,59,49,41,31,23,17,11,7,1,19,37,47,0,59,49,47,43,41,31,17,13,11,7,37,0,53,49,43,37,23,19,13,7,1,29,
    31,41,59,0,59,49,47,41,23,19,17,13,7,1,43,53,0,59,49,43,37,29,17,13,7,1,19,47,53,0,59,53,49,47,43,31,
    29,23,11,17,0,47,43,41,37,31,23,19,17,11,1,13,29,53,0,59,47,41,37,31,23,19,11,7,17,29,0,53,47,43,41,
    17,13,11,1,23,31,37,49,0,53,47,43,41,29,19,7,1,17,31,37,49,59,0,49,43,37,19,17,1,23,29,47,53,0,59,53,
    43,41,31,17,7,1,11,13,19,29 };

static const int numb1_125 = 319;
static const uint8_t b1_125[319] = {
    23,19,13,11,1,7,17,29,31,0,59,47,43,41,37,31,29,19,13,7,1,11,23,0,59,53,43,41,37,31,23,17,11,7,1,19,
    29,49,0,53,49,47,43,31,23,19,11,7,1,13,37,59,0,59,53,43,37,31,29,23,17,13,11,1,47,0,59,49,41,31,23,
    17,11,7,1,19,37,47,0,59,49,47,43,41,31,17,13,11,7,37,0,53,49,43,37,23,19,13,7,1,29,31,41,59,0,59,49,
    47,41,23,19,17,13,7,1,43,53,0,59,49,43,37,29,17,13,7,1,19,47,53,0,59,53,49,47,43,31,29,23,11,17,0,47,
    43,41,37,31,23,19,17,11,1,13,29,53,0,59,47,41,37,31,23,19,11,7,17,29,0,53,47,43,41,17,13,11,1,23,31,
    37,49,0,53,47,43,41,29,19,7,1,17,31,37,49,59,0,49,43,37,19,17,1,23,29,47,53,0,59,53,43,41,31,17,7,1,
    11,13,19,29,0,59,53,49,47,37,29,11,13,17,23,31,0,59,43,41,37,29,23,17,13,1,31,47,0,59,53,49,47,41,37,
    31,19,13,7,11,17,29,43,0,47,29,19,11,7,1,41,43,59,0,53,49,37,23,13,11,7,1,17,19,29,41,43,59,0,59,49,
    41,37,23,13,1,7,11,29,43,47,53,0,59,53,49,31,23,13,7,1,17,29,43,47,0,59,31,29,19,11,7,37,49,53 };

static const int numb1_165 = 425;
static const uint8_t b1_165[425] = {
    13,7,1,11,19,47,59,0,59,49,43,37,31,29,23,19,17,7,11,13,47,53,0,53,47,41,37,31,23,19,11,1,13,29,43,
    59,0,53,49,41,37,31,19,17,1,7,23,29,47,59,0,59,53,47,43,41,29,19,17,13,7,1,23,31,49,0,53,47,41,37,29,
    23,19,11,7,17,31,43,49,59,0,47,43,41,37,23,19,17,13,7,11,29,53,0,53,49,43,37,29,23,11,7,1,13,19,31,41,
    0,53,49,47,43,37,31,23,17,11,13,41,0,59,47,43,37,31,29,23,11,1,17,19,41,0,59,53,19,13,7,1,29,43,47,49,
    0,53,49,47,41,29,19,17,13,11,7,1,23,31,43,59,0,53,49,41,37,23,19,13,11,7,1,17,43,47,0,47,43,41,31,19,
    17,7,1,13,37,49,0,59,49,37,29,13,1,7,11,17,19,41,47,53,0,49,47,31,29,7,1,13,17,19,23,37,59,0,47,37,31,
    19,17,13,11,1,29,41,43,53,0,59,41,17,13,7,1,19,23,31,47,49,53,0,59,53,47,43,31,29,7,1,11,17,37,41,49,
    0,49,43,37,23,19,13,1,7,17,0,59,49,41,37,31,29,23,1,11,13,53,0,53,43,41,37,29,23,17,13,11,7,1,19,31,49,
    0,53,43,31,29,23,19,17,1,13,37,41,59,0,53,43,37,31,23,13,1,17,29,59,0,59,49,41,37,23,19,11,1,7,29,0,59,
    43,17,13,11,1,7,23,29,37,41,49,0,49,47,43,41,29,1,7,13,19,23,31,59,0,59,49,47,31,29,13,7,37,41,43,0,49,
    41,29,23,13,11,7,1,17,19,31,43,53,0,53,47,43,37,29,23,17,1,11,13,31,41,49,59,0,53,47,41,19,13,11,1,17,
    23,43,0,53,49,47,37,23,19,11,7,17,29,31,43,0,53,31,19,17,13,7,1,29,37,59 };

static const int numb1_205 = 511;
static const uint8_t b1_205[511] = {
    1,23,41,0,59,53,49,47,37,23,19,17,13,1,7,29,43,0,53,49,41,31,29,19,17,11,7,1,13,37,59,0,49,47,29,23,
    13,7,1,17,31,37,43,0,59,49,47,43,37,31,29,17,13,7,1,11,19,53,0,59,53,49,41,37,23,13,1,11,17,19,29,43,
    47,0,53,49,47,43,23,19,11,1,7,17,37,41,0,59,53,41,37,31,29,19,17,11,1,13,43,47,0,53,47,41,19,17,7,1,
    11,23,31,43,59,0,59,53,41,31,13,11,7,1,17,29,37,0,49,43,37,29,11,1,13,17,19,23,41,0,59,49,47,43,41,37,
    31,19,7,1,13,23,29,53,0,53,49,43,41,37,31,29,23,13,7,17,19,47,59,0,49,47,37,29,23,17,11,7,13,19,31,41,
    53,0,59,43,29,23,19,17,13,11,1,41,0,59,37,31,23,17,13,11,7,1,19,29,43,53,0,49,47,43,41,31,19,17,1,7,11,
    13,23,0,47,43,37,29,13,11,7,1,17,19,23,31,59,0,59,37,31,29,23,19,13,1,7,11,41,47,53,0,53,49,43,31,23,
    17,13,41,59,0,59,53,31,19,17,1,7,11,23,37,47,49,0,59,53,47,43,41,37,31,23,19,17,11,1,0,59,53,49,47,31,
    17,13,7,1,11,29,37,0,53,43,31,17,13,7,1,29,41,49,0,53,49,41,29,23,11,7,1,19,31,47,0,47,43,41,29,23,19,
    7,1,11,49,0,59,31,29,23,17,11,7,1,13,41,43,0,59,43,37,17,1,7,11,13,19,41,49,0,59,53,43,41,37,31,29,23,
    13,11,1,47,0,59,53,47,31,19,17,13,1,7,11,29,37,43,49,0,49,43,41,31,17,13,7,11,23,37,53,0,53,49,41,23,
    19,13,11,7,1,17,37,59,0,49,47,43,37,31,29,23,1,7,41,0,59,43,41,37,31,17,13,11,7,47,49,0,59,49,47,37,31,
    29,19,17,7,1,0,53,47,37,19,13,1,11,31,41,0,49,47,37,23,17,13,11,7,19,31,53,0,59,53,47,29,13,11,7,1,23,
    41,0,49,47,41,37,19,11,13,17,23,29,31,43,0,59,29,19,13,1,41,43,47,53,0,59,53,43,41,37,23,17,11,7,1,13,
    29,49 };


static uint64_t uecm_stage2(uecm_pt *P, uint64_t rho, uint64_t n, uint32_t stg1_max, uint64_t s, uint64_t unityval)
{
    int b;
    int i, j, k;
    uecm_pt Pa1;
    uecm_pt *Pa = &Pa1;
    uecm_pt Pb[20];
    uecm_pt *Pd;
    const uint8_t *barray = 0;
    int numb;

#ifdef MICRO_ECM_VERBOSE_PRINTF
    ptadds = 0;
    stg1Doub = 0;
    stg1Add = 0;
#endif

    // this function has been written for MICRO_ECM_PARAM_D of 60, so you
    // probably don't want to change it.
    const int MICRO_ECM_PARAM_D = 60;

    //stage 2 init
    //Q = P = result of stage 1
    //compute [d]Q for 0 < d <= MICRO_ECM_PARAM_D
    Pd = &Pb[map[MICRO_ECM_PARAM_D]];

    uint64_t Pbprod[20];

    // [1]Q
    Pb[1].Z = P->Z;
    Pb[1].X = P->X;
    Pbprod[1] = uecm_mulredc(Pb[1].X, Pb[1].Z, n, rho);

    // [2]Q
    Pb[2].Z = P->Z;
    Pb[2].X = P->X;
    uint64_t diff1 = submod(P->X, P->Z, n);
    uint64_t sum1 = addmod(P->X, P->Z, n);
    uecm_udup(s, rho, n, sum1, diff1, &Pb[2]);
    Pbprod[2] = uecm_mulredc(Pb[2].X, Pb[2].Z, n, rho);

    /*
    MICRO_ECM_PARAM_D is small in tinyecm, so it is straightforward to just enumerate the needed
    points.  We can do it efficiently with two progressions mod 6.
    Pb[0] = scratch
    Pb[1] = [1]Q;
    Pb[2] = [2]Q;
    Pb[3] = [7]Q;   prog2
    Pb[4] = [11]Q;  prog1
    Pb[5] = [13]Q;  prog2
    Pb[6] = [17]Q;  prog1
    Pb[7] = [19]Q;  prog2
    Pb[8] = [23]Q;  prog1
    Pb[9] = [29]Q;  prog1
    Pb[10] = [30 == MICRO_ECM_PARAM_D]Q;
    Pb[11] = [31]Q; prog2
    Pb[12] = [37]Q; prog2
    Pb[13] = [41]Q; prog1
    Pb[14] = [43]Q; prog2
    Pb[15] = [47]Q; prog1
    Pb[16] = [49]Q; prog2
    Pb[17] = [53]Q; prog1
    Pb[18] = [59]Q; prog1

    two progressions with total of 17 adds to get 15 values of Pb.
    6 + 5(1) -> 11 + 6(5) -> 17 + 6(11) -> 23 + 6(17) -> 29 + 6(23) -> 35 + 6(29) -> 41 + 6(35) -> 47 + 6(41) -> 53 + 6(47) -> 59
    6 + 1(5) -> 7 + 6(1) -> 13 + 6(7) -> 19 + 6(13) -> 25 + 6(19) -> 31 + 6(25) -> 37 + 6(31) -> 43 + 6(37) -> 49

    we also need [2 MICRO_ECM_PARAM_D]Q = [60]Q
    to get [60]Q we just need one more add:
    compute [60]Q from [31]Q + [29]Q([2]Q), all of which we
    have after the above progressions are computed.

    we also need [A]Q = [((B1 + MICRO_ECM_PARAM_D) / (2 MICRO_ECM_PARAM_D) * 2 MICRO_ECM_PARAM_D]Q
    which is equal to the following for various common B1:
    B1      [x]Q
    65      [120]Q
    85      [120]Q
    125     [180]Q
    165     [180]Q
    205     [240]Q

    and we need [x-MICRO_ECM_PARAM_D]Q as well, for the above [x]Q.
    So far we are getting [x]Q and [x-2 MICRO_ECM_PARAM_D]Q each from prac(x,Q).
    There is a better way using progressions of [2 MICRO_ECM_PARAM_D]Q
    [120]Q = 2*[60]Q
    [180]Q = [120]Q + [60]Q([60]Q)
    [240]Q = 2*[120]Q
    [300]Q = [240]Q + [60]Q([180]Q)
    ...
    etc.

    */

    uecm_pt pt1, pt3;

    // Calculate all Pb: the following is specialized for MICRO_ECM_PARAM_D=60
    // [2]Q + [1]Q([1]Q) = [3]Q
    uecm_uadd(rho, n, Pb[1], Pb[2], Pb[1], &Pb[3]);        // <-- temporary

    // 2*[3]Q = [6]Q
    diff1 = submod(Pb[3].X, Pb[3].Z, n);
    sum1 = addmod(Pb[3].X, Pb[3].Z, n);
    uecm_udup(s, rho, n, sum1, diff1, &pt3);   // pt3 = [6]Q

    // [3]Q + [2]Q([1]Q) = [5]Q
    uecm_uadd(rho, n, Pb[3], Pb[2], Pb[1], &pt1);    // <-- pt1 = [5]Q
    Pb[3].X = pt1.X;
    Pb[3].Z = pt1.Z;

    // [6]Q + [5]Q([1]Q) = [11]Q
    uecm_uadd(rho, n, pt3, pt1, Pb[1], &Pb[4]);    // <-- [11]Q

    i = 3;
    k = 4;
    j = 5;
    while ((j + 12) < MICRO_ECM_PARAM_D)
    {
        // [j+6]Q + [6]Q([j]Q) = [j+12]Q
        uecm_uadd(rho, n, pt3, Pb[k], Pb[i], &Pb[map[j + 12]]);
        i = k;
        k = map[j + 12];
        j += 6;
    }

    // [6]Q + [1]Q([5]Q) = [7]Q
    uecm_uadd(rho, n, pt3, Pb[1], pt1, &Pb[3]);    // <-- [7]Q
    i = 1;
    k = 3;
    j = 1;
    while ((j + 12) < MICRO_ECM_PARAM_D)
    {
        // [j+6]Q + [6]Q([j]Q) = [j+12]Q
        uecm_uadd(rho, n, pt3, Pb[k], Pb[i], &Pb[map[j + 12]]);
        i = k;
        k = map[j + 12];
        j += 6;
    }

    // Pd = [2w]Q
    // [31]Q + [29]Q([2]Q) = [60]Q
    uecm_uadd(rho, n, Pb[9], Pb[10], Pb[2], Pd);   // <-- [60]Q

#ifdef MICRO_ECM_VERBOSE_PRINTF
    ptadds++;
#endif

    // make all of the Pbprod's
    for (i = 3; i < 19; i++)
    {
        Pbprod[i] = uecm_mulredc(Pb[i].X, Pb[i].Z, n, rho);
    }


    //initialize info needed for giant step
    // temporary - make [4]Q
    diff1 = submod(Pb[2].X, Pb[2].Z, n);
    sum1 = addmod(Pb[2].X, Pb[2].Z, n);
    uecm_udup(s, rho, n, sum1, diff1, &pt3);   // pt3 = [4]Q

    uecm_pt Pad;

    // Pd = [w]Q
    // [17]Q + [13]Q([4]Q) = [30]Q
    uecm_uadd(rho, n, Pb[map[17]], Pb[map[13]], pt3, &Pad);    // <-- [30]Q

    // [60]Q + [30]Q([30]Q) = [90]Q
    uecm_uadd(rho, n, *Pd, Pad, Pad, Pa);
    pt1.X = Pa->X;
    pt1.Z = Pa->Z;

    // [90]Q + [30]Q([60]Q) = [120]Q
    uecm_uadd(rho, n, *Pa, Pad, *Pd, Pa);
    Pd->X = Pa->X;
    Pd->Z = Pa->Z;

    // [120]Q + [30]Q([90]Q) = [150]Q
    uecm_uadd(rho, n, *Pa, Pad, pt1, Pa);

    // adjustment of Pa and Pad for larger B1.
    // Currently we have Pa=150, Pd=120, Pad=30
    if (stg1_max == 165)
    {
        // need Pa = 180, Pad = 60
        // [150]Q + [30]Q([120]Q) = [180]Q
        uecm_uadd(rho, n, *Pa, Pad, *Pd, Pa);

        diff1 = submod(Pad.X, Pad.Z, n);
        sum1 = addmod(Pad.X, Pad.Z, n);
        uecm_udup(s, rho, n, sum1, diff1, &Pad);   // Pad = [60]Q
    }
    else if (stg1_max == 205)
    {
        // need Pa = 210, Pad = 90.
        // have pt1 = 90

        diff1 = submod(Pad.X, Pad.Z, n);
        sum1 = addmod(Pad.X, Pad.Z, n);
        uecm_udup(s, rho, n, sum1, diff1, &Pad);   // Pad = [60]Q

        // [150]Q + [60]Q([90]Q) = [210]Q
        uecm_uadd(rho, n, *Pa, Pad, pt1, Pa);
        Pad.X = pt1.X;
        Pad.Z = pt1.Z;
    }

    //initialize accumulator and Paprod
    uint64_t acc = unityval;
    uint64_t Paprod = uecm_mulredc(Pa->X, Pa->Z, n, rho);

    if (stg1_max <= 70)
    {
        barray = b1_70;
        numb = numb1_70;
    }
    else if (stg1_max == 85)
    {
        barray = b1_85;
        numb = numb1_85;
    }
    else if (stg1_max == 125)
    {
        barray = b1_125;
        numb = numb1_125;
    }
    else if (stg1_max == 165)
    {
        barray = b1_165;
        numb = numb1_165;
    }
    else if (stg1_max == 205)
    {
        barray = b1_205;
        numb = numb1_205;
    }

    for (i = 0; i < numb; i++)
    {
        if (barray[i] == 0)
        {
            //giant step - use the addition formula for ECM
            pt1.X = Pa->X;
            pt1.Z = Pa->Z;

            //Pa + Pd
            uecm_uadd(rho, n, *Pa, *Pd, Pad, Pa);

            //Pad holds the previous Pa
            Pad.X = pt1.X;
            Pad.Z = pt1.Z;

            //and Paprod
            Paprod = uecm_mulredc(Pa->X, Pa->Z, n, rho);

            i++;
        }

        //we accumulate XrZd - XdZr = (Xr - Xd) * (Zr + Zd) + XdZd - XrZr
        //in CP notation, Pa -> (Xr,Zr), Pb -> (Xd,Zd)

        b = barray[i];
        // accumulate the cross product  (zimmerman syntax).
        // page 342 in C&P
        uint64_t tt1 = submod(Pa->X, Pb[map[b]].X, n);
        uint64_t tt2 = addmod(Pa->Z, Pb[map[b]].Z, n);
        uint64_t tt3 = uecm_mulredc(tt1, tt2, n, rho);
        tt1 = addmod(tt3, Pbprod[map[b]], n);
        tt2 = submod(tt1, Paprod, n);

        uint64_t tmp = uecm_mulredc(acc, tt2, n, rho);
        if (tmp == 0)
            break;
        acc = tmp;
    }

    return acc;
}


static void microecm(uint64_t n, uint64_t *f, uint32_t B1, uint32_t B2,
                     uint32_t curves, uint64_t *ploc_lcg)
{
    //attempt to factor n with the elliptic curve method
    //following brent and montgomery's papers, and CP's book
    int curve;
    int found = 0;
    int result;
    uecm_pt P;
    uint64_t tmp1;

    uint64_t rho = (uint64_t)0 - uecm_multiplicative_inverse(n);

    uint32_t stg1_max = B1;
//    uint32_t stg2_max = B2;

//    uint64_t unityval = u64div(1, n);
    // Let R = 2^64.  We can see R%n ≡ (R-n)%n  (mod n)
    uint64_t unityval = ((uint64_t)0 - n) % n;   // unityval ≡ R  (mod n)

    uint64_t two = addmod(unityval, unityval, n);
    uint64_t four = addmod(two, two, n);
    uint64_t five = addmod(unityval, four, n);
    uint64_t eight = addmod(four, four, n);
    uint64_t sixteen = addmod(eight, eight, n);
    uint64_t two_8 = uecm_sqrredc(sixteen, n, rho);   // R*2^8         (mod n)
    uint64_t two_16 = uecm_sqrredc(two_8, n, rho);    // R*2^16        (mod n)
    uint64_t two_32 = uecm_sqrredc(two_16, n, rho);   // R*2^32        (mod n)
    uint64_t Rsqr = uecm_sqrredc(two_32, n, rho);     // R*2^64 ≡ R*R  (mod n)

    *f = 1;
    for (curve = 0; (uint32_t)curve < curves; curve++)
    {
#ifdef MICRO_ECM_VERBOSE_PRINTF
        stg1Add = 0;
        stg1Doub = 0;
            printf("commencing curve %d of %u\n", curve, curves);
#endif
        uint64_t s;
        uint64_t likely_gcd = uecm_build(&P, rho, n, ploc_lcg, &s, five, Rsqr);
        if (likely_gcd > 1)
        {
            // If the gcd gave us a factor, we're done.  If not, since gcd != 1
            // the inverse calculated in uecm_build would have bogus, and so this
            // curve is probably set up for failure (hence we continue).
            if (likely_gcd == n || n % likely_gcd != 0)
                continue;
            *f = likely_gcd;
            break;
        }

#ifdef MICRO_ECM_VERBOSE_PRINTF
        {
            printf("curve parameters:\n");
            printf("\tn = %" PRIu64 "\n", n);
            printf("\trho = %" PRIu64 "\n", rho);
            printf("\tx = %" PRIx64 "\n", P.X);
            printf("\tz = %" PRIx64 "\n", P.Z);
            printf("\tb = %" PRIx64 "\n", s);
        }
#endif
        uecm_stage1(rho, n, &P, (uint64_t)stg1_max, s);
        result = uecm_check_factor(P.Z, n, &tmp1);

#ifdef MICRO_ECM_VERBOSE_PRINTF
        {
            printf("after stage1: P = %" PRIx64 ", %" PRIx64 "\n", P.X, P.Z);
        }
#endif
        if (result == 1)
        {
#ifdef MICRO_ECM_VERBOSE_PRINTF
                printf("\nfound factor %" PRIx64 " in stage 1\n", tmp1);
#endif
            *f = tmp1;
            break;
        }

        if (B2 > B1)
        {
            uint64_t stg2acc = uecm_stage2(&P, rho, n, stg1_max, s, unityval);

#ifdef MICRO_ECM_VERBOSE_PRINTF
            {
                printf("after stage2: A = %" PRIx64 "\n", stg2acc);
            }
                uint32_t paired = 0;
                uint64_t numprimes = 0;
                printf("performed %d pair-multiplies for %" PRIu64 " primes in stage 2\n",
                    paired, numprimes);
                printf("performed %u point-additions and %u point-doubles in stage 2\n",
                    ptadds + stg1Add, stg1Doub);
#endif
            result = uecm_check_factor(stg2acc, n, &tmp1);

            if (result == 1)
            {
#ifdef MICRO_ECM_VERBOSE_PRINTF
                    printf("\nfound factor %" PRIx64 " in stage 2\n", tmp1);
#endif
                *f = tmp1;
                break;
            }
        }
    }

    return;
}


static uint64_t uecm_dispatch(uint64_t n, int targetBits, uint64_t *ploc_lcg)
{
    int B1, curves;
    uint64_t f64 = 1;

#ifdef MICRO_ECM_EXPECT_ARBITRARY_SIZE_FACTORS
    // try fast attempts to find possible small factors.
    {
        B1 = 47;
        curves = 1;
        microecm(n, &f64, B1, 25 * B1, curves, ploc_lcg);
        if (f64 > 1)
            return f64;
    }
    {
        B1 = 70;
        curves = 1;
        microecm(n, &f64, B1, 25 * B1, curves, ploc_lcg);
        if (f64 > 1)
            return f64;
    }
    if (targetBits > 58)
    {
        B1 = 125;
        curves = 1;
        microecm(n, &f64, B1, 25 * B1, curves, ploc_lcg);
        if (f64 > 1)
            return f64;
    }
#endif

    if (targetBits <= 48)
    {
        // multi-thread issue here...
        //f64 = LehmanFactor(n, 0, 0, 0);
        B1 = 70;
        curves = 32;
        microecm(n, &f64, B1, 25 * B1, curves, ploc_lcg);
    }
    else if (targetBits <= 52)
    {
        B1 = 85;
        curves = 32;
        microecm(n, &f64, B1, 25 * B1, curves, ploc_lcg);
    }
    else if (targetBits <= 58)
    {
        B1 = 125;
        curves = 32;
        microecm(n, &f64, B1, 25 * B1, curves, ploc_lcg);
    }
    else if (targetBits <= 62)
    {
        B1 = 165;
        curves = 42;
        microecm(n, &f64, B1, 25 * B1, curves, ploc_lcg);
    }
    else if (targetBits <= 64)
    {
        B1 = 205;
        curves = 42;
        microecm(n, &f64, B1, 25 * B1, curves, ploc_lcg);
    }

    return f64;
}


static int uecm_get_bits(uint64_t n)
{
    int i = 0;
    while (n != 0)
    {
        n >>= 1;
        i++;
    }
    return i;
}




/* ----- PREFERRED API ------- */

// getfactor_uecm() returns 1 if unable to find a factor of q64,
// Otherwise it returns a factor of q64.
// Prior to your first call of getfactor_uecm(), set *pran = 0  (or set it to
// some other arbitrary value); after that, don't change *pran.
// FYI: *pran is used within this file by a random number generator, and it
// holds the current value of a pseudo random sequence.  Your first assigment
// to *pran seeds the sequence, and after seeding it you don't want to
// change *pran, since that would restart the sequence.
uint64_t getfactor_uecm(uint64_t q64, uint64_t *pran)
{
    if (q64 % 2 == 0)
        return 2;
    int bits = uecm_get_bits(q64);
    return uecm_dispatch(q64, bits, pran);
}



/* ----- DEPRECATED API ------- */
// This api is deprecated because it is not thread safe.
static uint64_t LOC_LCG = 0;
void init_uecm(uint64_t lcg)  // This is a deprecated function
{
    LOC_LCG = lcg;
    return;
}
uint64_t do_uecm(uint64_t q64)
{
    return getfactor_uecm(q64, &LOC_LCG);
}
/* ----- END OF DEPRECATED API ------- */
