/*
Copyright (c) 2019, Ben Buhrow
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
*/

/* Elliptic Curve Method: toplevel and stage 1 routines.

Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011,
2012, 2016 Paul Zimmermann, Alexander Kruppa, Cyril Bouvier, David Cleaver.

This file is part of the ECM Library.

The ECM Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The ECM Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the ECM Library; see the file COPYING.LIB.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#include "yafu_ecm.h"
#include "avx_ecm.h"
#include "soe.h"
#include "arith.h"
#include "factor.h"
#include "ytools.h"
#include "nfs_impl.h"
#include <stdio.h>
#include <ecm.h>

// define function pointers to the type of reduction needed
void(*vecmulmod_ptr)(vec_bignum_t*, vec_bignum_t*, vec_bignum_t*, vec_bignum_t*, vec_bignum_t*, vec_monty_t*);
void(*vecsqrmod_ptr)(vec_bignum_t*, vec_bignum_t*, vec_bignum_t*, vec_bignum_t*, vec_monty_t*);
void(*vecaddmod_ptr)(vec_bignum_t*, vec_bignum_t*, vec_bignum_t*, vec_monty_t*);
void(*vecsubmod_ptr)(vec_bignum_t*, vec_bignum_t*, vec_bignum_t*, vec_monty_t*);
void(*vecaddsubmod_ptr)(vec_bignum_t*, vec_bignum_t*, vec_bignum_t*, vec_bignum_t*, vec_monty_t*);

void vecLucasV(vec_monty_t* mdata, uint64_t j, vec_bignum_t *v)
{
    //compute v_j mod p using the duplication and addition formula's
    //for lucas V sequences.  Algorithm described by Mersenne Wiki entry.
    //montgomery ladder powering algorithm...

    uint8_t t, bi[64];
    int i;
    vec_bignum_t *x = mdata->mtmp1;
    vec_bignum_t* y = mdata->mtmp2;
    vec_bignum_t* n = mdata->n;
    vec_bignum_t* s = mdata->mtmp4;
    vec_bignum_t* tmp = mdata->mtmp3;
    vec_bignum_t* two;

    //x=h=Vn 
    vecCopy(v, x);

    two = vecInit(mdata->NWORDS);
    vecaddmod_ptr(mdata->one, mdata->one, two, mdata);

    //y=(B^2-2) mod N
    vecsqrmod_ptr(v, tmp, n, s, mdata);
    vecsubmod_ptr(tmp, two, y, mdata);

    //get j in binary form
    i = 0;
    while (j > 0)
    {
        i++;
        bi[i] = j % 2;
        j = j / 2;
    }
    t = i;

    //compute loop
    //for each bit of M to the right of the most significant bit
    for (i = t - 1; i >= 1; i--)
    {
        //if the bit is 1
        if (bi[i])
        {
            //x=(x*y-B) mod N 
            vecmulmod_ptr(x, y, tmp, n, s, mdata);
            vecsubmod_ptr(tmp, v, x, mdata);
            //y=(y^2-2) mod N 
            vecsqrmod_ptr(y, tmp, n, s, mdata);
            vecsubmod_ptr(tmp, two, y, mdata);
        }
        else
        {
            //y=(x*y-B) mod N 
            vecmulmod_ptr(x, y, tmp, n, s, mdata);
            vecsubmod_ptr(tmp, v, y, mdata);
            //x=(x^2-2) mod N 
            vecsqrmod_ptr(x, tmp, n, s, mdata);
            vecsubmod_ptr(tmp, two, x, mdata);
        }
    }
    vecCopy(x, v);

    vecFree(two);
    return;
}

// P <- V_2(Q)
static void pp1_dup(vec_bignum_t* P, vec_bignum_t* Q, vec_monty_t* mdata)
{
    vecsqrmod_ptr(Q, P, mdata->n, mdata->mtmp1, mdata);
    vecsubmod_ptr(P, mdata->mtmp2, P, mdata);
}

// P <- V_{m+n} where Q = V_m, R = V_n, S = V_{m-n}.
static void pp1_add(vec_bignum_t* P, vec_bignum_t* Q, vec_bignum_t* R, 
    vec_bignum_t* S, vec_monty_t* mdata)
{
    vecmulmod_ptr(Q, R, mdata->mtmp1, mdata->n, mdata->mtmp3, mdata);
    vecsubmod_ptr(mdata->mtmp1, S, P, mdata);
}

// computes V_k(P) from P=A and puts the result in P=A. Assumes k>2.
void pp1_prac(vec_bignum_t* A, uint64_t k, vec_bignum_t* B,
    vec_bignum_t* C, vec_bignum_t* T, vec_bignum_t* T2, vec_monty_t* mdata)
{
    // same conditions as ecm prac, but using Lucas chain duplication and
    // addition rules instead of elliptic curve duplication and addition.
    uint64_t d, e, r;
    static double val = 0.61803398874989485; /* 1/(golden ratio) */
    base_t* swapptr;

    d = k;
    r = (uint64_t)((double)d * val + 0.5);

    /* first iteration always begins by Condition 3, then a swap */
    d = k - r;
    e = 2 * r - k;

    // mtmp2 holds the value 2 throughout (used in duplication) - don't reuse it.
    vecaddmod_ptr(mdata->one, mdata->one, mdata->mtmp2, mdata);
    vecCopy(A, B);
    vecCopy(A, C);
    pp1_dup(A, A, mdata);

    while (d != e)
    {
        if (d < e)
        {
            r = d;
            d = e;
            e = r;
            swapptr = A->data;
            A->data = B->data;
            B->data = swapptr;
        }
        /* do the first line of Table 4 whose condition qualifies */
        if (d - e <= e / 4 && ((d + e) % 3) == 0)
        {
            d = (2 * d - e) / 3;
            e = (e - d) / 2;
            pp1_add(T, A, B, C, mdata);
            pp1_add(T2, T, A, B, mdata);
            pp1_add(B, B, T, A, mdata);
            swapptr = A->data;
            A->data = T2->data;
            T2->data = swapptr;
        }
        else if (d - e <= e / 4 && (d - e) % 6 == 0)
        {
            d = (d - e) / 2;
            pp1_add(B, A, B, C, mdata);
            pp1_dup(A, A, mdata);
        }
        else if ((d + 3) / 4 <= e) /* <==>  (d <= 4 * e) */
        {
            d -= e;
            pp1_add(C, B, A, C, mdata);
            swapptr = B->data;
            B->data = C->data;
            C->data = swapptr;
        }
        else if ((d + e) % 2 == 0)
        {
            d = (d - e) / 2;
            pp1_add(B, B, A, C, mdata);
            pp1_dup(A, A, mdata);
        }
        /* d+e is now odd */
        else if (d % 2 == 0)
        {
            d /= 2;
            pp1_add(C, C, A, B, mdata);
            pp1_dup(A, A, mdata);
        }
        /* d is odd, e even */
        else if (d % 3 == 0)
        {
            d = d / 3 - e;
            pp1_dup(T, A, mdata);
            pp1_add(T2, A, B, C, mdata);
            pp1_add(A, T, A, A, mdata);
            pp1_add(C, T, T2, C, mdata);
            swapptr = B->data;
            B->data = C->data;
            C->data = swapptr;
        }
        else if ((d + e) % 3 == 0) /* d+e <= val[i]*k < k < 2^32 */
        {
            d = (d - 2 * e) / 3;
            pp1_add(T, A, B, C, mdata);
            pp1_add(B, T, A, B, mdata);
            pp1_dup(T, A, mdata);
            pp1_add(A, A, T, A, mdata);
        }
        else if ((d - e) % 3 == 0)
        {
            d = (d - e) / 3;
            pp1_add(T, A, B, C, mdata);
            pp1_add(C, C, A, B, mdata);
            swapptr = B->data;
            B->data = T->data;
            T->data = swapptr;
            pp1_dup(T, A, mdata);
            pp1_add(A, A, T, A, mdata);
        }
        else /* necessarily e is even */
        {
            e /= 2;
            pp1_add(C, C, B, A, mdata);
            pp1_dup(B, B, mdata);
        }
    }

    pp1_add(A, A, B, C, mdata);

    if (d != 1)
    {
        printf("problem in pp1_prac, d = %lu\n", d);
    }
    return;
}


void vecmodexp_1(vec_bignum_t* A, uint64_t k, vec_monty_t* mdata)
{
    // just a stupid RL binary method
    vec_bignum_t* r = mdata->mtmp1;

    if (k == 1)
        return;

    vecCopy(mdata->one, r);

    while (k > 0)
    {
        if (k & 1)
        {
            vecmulmod_ptr(r, A, r, mdata->n, mdata->mtmp2, mdata);
        }
        k >>= 1;
        vecsqrmod_ptr(A, A, mdata->n, mdata->mtmp2, mdata);
    }

    vecCopy(r, A);
    return;
}

#define MAX_HEIGHT 32

unsigned int compute_s(mpz_t s, uint64_t * primes, uint64_t nump, uint64_t B1)
{
    mpz_t acc[MAX_HEIGHT]; /* To accumulate products of prime powers */
    mpz_t ppz;
    unsigned int i, j, it;
    uint64_t pi = primes[0], pp, maxpp, qi;

    for (j = 0; j < MAX_HEIGHT; j++)
        mpz_init(acc[j]); /* sets acc[j] to 0 */
    mpz_init(ppz);

    i = 0;
    while (pi <= B1)
    {
        pp = qi = pi;
        maxpp = B1 / qi;

        while (pp <= maxpp)
            pp *= qi;

        mpz_set_ui(ppz, pp);

        if ((i & 1) == 0)
            mpz_set(acc[0], ppz);
        else
            mpz_mul(acc[0], acc[0], ppz);

        j = 0;
        /* We have accumulated i+1 products so far. If bits 0..j of i are all
           set, then i+1 is a multiple of 2^(j+1). */
        while ((i & (1 << j)) != 0)
        {
            /* we use acc[MAX_HEIGHT-1] as 0-sentinel below, thus we need
               j+1 < MAX_HEIGHT-1 */
            if ((i & (1 << (j + 1))) == 0) /* i+1 is not multiple of 2^(j+2),
                                              thus add[j+1] is "empty" */
                mpz_swap(acc[j + 1], acc[j]); /* avoid a copy with mpz_set */
            else
                mpz_mul(acc[j + 1], acc[j + 1], acc[j]); /* accumulate in acc[j+1] */
            mpz_set_ui(acc[j], 1);
            j++;
        }

        i++;
        if (i >= nump)
            break;
        pi = primes[i];
    }
    it = i;

    for (mpz_set(s, acc[0]), j = 1; mpz_cmp_ui(acc[j], 0) != 0; j++)
        mpz_mul(s, s, acc[j]);

    for (i = 0; i < MAX_HEIGHT; i++)
        mpz_clear(acc[i]);
    mpz_clear(ppz);

    return it;
}

int get_winsize(uint32_t bits)
{
    // the window size is based on minimizing the total number of multiplications
    // in the windowed exponentiation.  experiments show that this is best;
    // the growing size of the table doesn't change the calculus, at least
    // on the KNL.
    int size;
    int muls;
    int minmuls = 99999999;
    int minsize = 4;

    for (size = 2; size <= MAX_WINSIZE; size++)
    {
        muls = (bits / size) + (1 << size);
        if (muls < minmuls)
        {
            minmuls = muls;
            minsize = size;
        }
    }

    return minsize;
}

void vecslidingmodexp(vec_bignum_t* d, vec_bignum_t* b, mpz_t e, vec_bignum_t* m,
    vec_bignum_t* s, vec_bignum_t* one, vec_monty_t* mdata)
{
    // d = b^e mod m
    // all b and e vector elements can be different.
    // all m elements can be different
    // proceed in a left to right fashion.
    int i, j;
    int bit = mpz_sizeinbase(e, 2) - 1;
    // the bit string length.  needs to divide the word size (e.g., 32);
    int k = get_winsize(bit);
    int mask;
    int nsqr = 0, nmul = 0;
    int bstr;
    uint8_t done[1 << MAX_WINSIZE];

    vec_bignum_t* t = mdata->mtmp3;
    vec_bignum_t** g = mdata->g;
    vec_bignum_t* w = mdata->mtmp4;

    mask = 0;
    for (j = 0; j < k; j++)
    {
        mask = (mask << 1) | 1;
    }

    vecCopy(b, g[1]);
    vecsqrmod_ptr(b, g[2], m, s, mdata);
    nsqr++;
    for (i = 1; i < (1 << (k - 1)); i++)
    {
        vecmulmod_ptr(g[2 * i - 1], g[2], g[2 * i + 1], m, s, mdata);
        nmul++;
    }

    vecCopyn(one, d, mdata->NWORDS);

    while (bit >= 0)
    {
        if (mpz_tstbit(e, bit) == 0)
        {
            vecsqrmod_ptr(d, d, m, s, mdata);
            nsqr++;
            bit--;
        }
        else
        {
            int thisk;

            // find the longest bitstring ending in 1 that we can make
            int a = bit - k + 1;
            if (a < 0) a = 0;

            a = mpz_scan1(e, a);
            thisk = bit - a + 1;

            bstr = 1;
            for (j = 1; j < thisk; j++)
            {
                bstr <<= 1;
                bstr |= mpz_tstbit(e, bit - j);
            }

            for (j = 0; j < thisk; j++)
            {
                vecsqrmod_ptr(d, d, m, s, mdata);
                nsqr++;
            }

            vecmulmod_ptr(d, g[bstr], d, m, s, mdata);
            nmul++;
            bit -= thisk;
        }
    }

    printf("pm1: modexp performed %d muls and %d sqrs\n", nmul, nsqr);

    return;
}

void vecmodexp(vec_bignum_t* d, vec_bignum_t* b, mpz_t e, vec_bignum_t* m,
    vec_bignum_t* s, vec_bignum_t* one, vec_monty_t* mdata)
{
    // d = b^e mod m
    // all b and e vector elements can be different.
    // all m elements can be different
    // proceed in a left to right fashion.
    int i, j;
    int bit = mpz_sizeinbase(e, 2) - 1;
    // the bit string length.  needs to divide the word size (e.g., 32);
    int k = get_winsize(bit);
    int mask;
    int nsqr = 0, nmul = 0;
    int bstr;
    uint8_t done[1<<MAX_WINSIZE];

    vec_bignum_t* t = mdata->mtmp3;
    vec_bignum_t** g = mdata->g;
    vec_bignum_t* w = mdata->mtmp4;

    mask = 0;
    for (j = 0; j < k; j++)
    {
        mask = (mask << 1) | 1;
    }

    for (j = 16; j < (1 << k); j++)
    {
        done[j] = 0;
    }

    // precomputation.  smallest window is 4 bits, so these
    // computations can be pulled out of the loop
    vecCopy(b, g[1]);
    vecsqrmod_ptr(g[1], g[2], m, s, mdata);
    vecsqrmod_ptr(g[2], g[4], m, s, mdata);
    vecmulmod_ptr(g[2], b, g[3], m, s, mdata);
    vecsqrmod_ptr(g[3], g[6], m, s, mdata);
    vecmulmod_ptr(g[6], b, g[7], m, s, mdata);
    vecsqrmod_ptr(g[6], g[12], m, s, mdata);
    vecmulmod_ptr(g[4], b, g[5], m, s, mdata);
    vecsqrmod_ptr(g[4], g[8], m, s, mdata);
    vecsqrmod_ptr(g[5], g[10], m, s, mdata);
    vecsqrmod_ptr(g[7], g[14], m, s, mdata);
    vecmulmod_ptr(g[8], b, g[9], m, s, mdata);
    vecmulmod_ptr(g[10], b, g[11], m, s, mdata);
    vecmulmod_ptr(g[12], b, g[13], m, s, mdata);
    vecmulmod_ptr(g[14], b, g[15], m, s, mdata);

#ifdef DEBUGLANE

    print_vechex(g[2]->data, DEBUGLANE, NWORDS, "g[2]: ");
    print_vechex(g[4]->data, DEBUGLANE, NWORDS, "g[4]: ");
    print_vechex(g[3]->data, DEBUGLANE, NWORDS, "g[3]: ");
    print_vechex(g[6]->data, DEBUGLANE, NWORDS, "g[6]: ");
    print_vechex(g[7]->data, DEBUGLANE, NWORDS, "g[7]: ");
    print_vechex(g[12]->data, DEBUGLANE, NWORDS, "g[12]: ");
    print_vechex(g[5]->data, DEBUGLANE, NWORDS, "g[5]: ");
    print_vechex(g[8]->data, DEBUGLANE, NWORDS, "g[8]: ");
    print_vechex(g[10]->data, DEBUGLANE, NWORDS, "g[10]: ");
    print_vechex(g[14]->data, DEBUGLANE, NWORDS, "g[14]: ");
    print_vechex(g[9]->data, DEBUGLANE, NWORDS, "g[9]: ");
    print_vechex(g[11]->data, DEBUGLANE, NWORDS, "g[11]: ");
    print_vechex(g[13]->data, DEBUGLANE, NWORDS, "g[13]: ");
    print_vechex(g[15]->data, DEBUGLANE, NWORDS, "g[15]: ");

#endif

    nsqr = 7;
    nmul = 7;

    for (i = 16; i < (1 << k); i++)
    {
        int q = i;

        if (done[i])
            continue;

        vecmulmod_ptr(g[i - 1], b, g[i], m, s, mdata);
        nmul++;

        while ((q * 2) < (1 << k))
        {
            if (!done[2 * q])
            {
                vecsqrmod_ptr(g[q], g[2 * q], m, s, mdata);
                nsqr++;
                done[2 * q] = 1;
            }
            q *= 2;
        }
    }

    //printf("init performed %d sqrs and %d muls for window size %d\n", nsqr, nmul, k);

    vecCopyn(one, d, mdata->NWORDS);

    while (bit >= 0)
    {
        if (bit < k)
        {
            mask = 0x0;
            for (j = 0; j < (bit + 1); j++)
            {
                vecsqrmod_ptr(d, d, m, s, mdata);
                nsqr++;
                mask = (mask << 1) | 1;
            }

            bstr = mpz_get_ui(e) & mask;
        }
        else
        {
            bstr = 0;
            for (j = 0; j < k; j++)
            {
                bstr <<= 1;
                bstr |= mpz_tstbit(e, bit - j);
            }

            for (j = 0; j < k; j++)
            {
                vecsqrmod_ptr(d, d, m, s, mdata);
                nsqr++;
            }
        }
        
        if (bstr > 0)
        {
            nmul++;
            vecmulmod_ptr(d, g[bstr], d, m, s, mdata);
        }

        bit -= k;
    }

    printf("pm1: modexp performed %d muls and %d sqrs\n", nmul, nsqr);
    return;
}

// test case for B1=20000
// 138104364401754474985334425082633970310479697396601597658245297623349837820784076171
// has factor 4874944133233 that rho doesn't find and which factors over B1 by either
// p+1 or p-1

// test case for B1>21k
// 136416892096153201112953487577103662921779668537279516591234800748153097
// has factor 13166819699867 with smooth B1 for p+1 but not for p-1.

// should find factor 8678835016483 
// of 2539450738391012463035957619341325187941394170998736241
// in stage 2 with B1>5k, <54k

void vecPP1(fact_obj_t* fobj)
{
    // run William's P+1 factoring method on 8 inputs simultaneously.
    // The inputs should be about the same size so that all can use 
    // the same multplier lengths.  We use the longest one necessary
    // for all inputs.
    int i;
    vec_bignum_t* vecV;
    mpz_t N[VECLEN], r, g;
    int verbose = fobj->VFLAG;
    uint32_t maxbits, nwords, nblocks, nmaxwords = 0, maxinbits = 0;
    vec_monty_t* montyconst;
    // these track the range over which we currently have a prime list.
    uint64_t rangemin;
    uint64_t rangemax;
    uint64_t* ppm1_primes;
    uint64_t ppm1_nump;
    uint64_t ppm1_minp;
    uint64_t ppm1_maxp;
    uint64_t STAGE1_MAX = fobj->pp1_obj.B1;
    uint64_t STAGE2_MAX = 100 * fobj->pp1_obj.B1;
    uint64_t PRIME_RANGE = 100000000;
    uint32_t size_n;
    int isMersenne = 0, allMersenne = 1, forceNoMersenne = 0;
    uint64_t vecbase[VECLEN];
    
    rangemin = 0;
    rangemax = MIN(STAGE1_MAX + 1000, (uint64_t)PRIME_RANGE);

    soe_staticdata_t* sdata = soe_init(0, 1, 32768);
    ppm1_primes = soe_wrapper(sdata, rangemin, rangemax, 0, &ppm1_nump, 0, 0);
    ppm1_minp = ppm1_primes[0];
    ppm1_maxp = ppm1_primes[ppm1_nump - 1];

    if (verbose > 1)
    {
        printf("pp1: found %lu primes in range [%lu : %lu]\n", ppm1_nump, rangemin, rangemax);
    }

    mpz_init(r);
    mpz_init(g);

    int k = 0;
    while (k < fobj->pp1_obj.vecnum)
    {
        mpz_init(N[k]);
        mpz_set(N[k], fobj->pp1_obj.vecn[k]);

        // check for Mersenne inputs
        size_n = mpz_sizeinbase(N[k], 2);

        if (size_n > maxinbits)
            maxinbits = size_n;

        /*
        for (i = size_n; i < 2048; i++)
        {
            mpz_set_ui(r, 1);
            mpz_mul_2exp(r, r, i);
            mpz_sub_ui(r, r, 1);
            mpz_mod(g, r, N[k]);
            if (mpz_cmp_ui(g, 0) == 0)
            {
                size_n = i;
                isMersenne = 1;
                break;
            }

            mpz_set_ui(r, 1);
            mpz_mul_2exp(r, r, i);
            mpz_add_ui(r, r, 1);
            mpz_mod(g, r, N[k]);
            if (mpz_cmp_ui(g, 0) == 0)
            {
                size_n = i;
                isMersenne = -1;
                break;
            }

            // detect pseudo-Mersennes
            mpz_set_ui(r, 1);
            mpz_mul_2exp(r, r, i);
            mpz_mod(g, r, N[k]);
            if (mpz_sizeinbase(g, 2) < DIGITBITS)
            {
                size_n = i;
                isMersenne = mpz_get_ui(g);
                break;
            }
        }

        // if the input is Mersenne and still contains algebraic factors, remove them.
        if (abs(isMersenne) == 1)
        {
#ifdef USE_NFS
            snfs_t poly;

            snfs_init(&poly);
            poly.form_type = SNFS_BRENT;
            mpz_set_ui(poly.base1, 2);
            poly.coeff1 = 1;
            poly.coeff2 = -isMersenne;
            poly.exp1 = size_n;
            find_primitive_factor(&poly, fobj->primes, fobj->num_p, verbose);

            mpz_tdiv_q(g, N[k], poly.primitive);
            mpz_gcd(N[k], N[k], poly.primitive);

            if (mpz_cmp_ui(g, 1) > 0)
            {
                add_to_factor_list(fobj->factors, g,
                    fobj->VFLAG, fobj->NUM_WITNESSES);

                int fid = fobj->factors->num_factors - 1;

                fobj->factors->factors[fid].tid = 0;
                fobj->factors->factors[fid].curve_num = 0;
                fobj->factors->factors[fid].sigma = 0;
                fobj->factors->factors[fid].vid = 0;
                fobj->factors->factors[fid].method = 0;
            }

            snfs_clear(&poly);
#endif
        }
        */

        // compute NBLOCKS if using the actual size of the input (non-Mersenne)
        if (DIGITBITS == 52)
        {
            maxbits = 208;
            while (maxbits <= mpz_sizeinbase(N[k], 2))
            {
                maxbits += 208;
            }
        }
        else
        {
            maxbits = 128;
            while (maxbits <= mpz_sizeinbase(N[k], 2))
            {
                maxbits += 128;
            }
        }

        nwords = maxbits / DIGITBITS;
        nblocks = nwords / BLOCKWORDS;

        // and compute NBLOCKS if using Mersenne mod
        if (DIGITBITS == 52)
        {
            maxbits = 208;
            while (maxbits <= size_n)
            {
                maxbits += 208;
            }
        }
        else
        {
            maxbits = 128;
            while (maxbits <= size_n)
            {
                maxbits += 128;
            }
        }

        if (verbose > 0)
        {
            gmp_printf("pp1: commencing parallel P+1 to B1 = %lu on %Zd\n", STAGE1_MAX, N[k]);
        }

        if ((double)nwords / ((double)maxbits / (double)DIGITBITS) < 0.7)
        {
            if (verbose > 1)
            {
                printf("Mersenne input 2^%d - 1 determined to be faster by REDC\n",
                    size_n);
            }
            forceNoMersenne = 1;
        }
        else
        {
            nwords = maxbits / DIGITBITS;
            nblocks = nwords / BLOCKWORDS;
        }

        if (forceNoMersenne)
        {
            isMersenne = 0;
            size_n = mpz_sizeinbase(N[k], 2);
        }

        if (nwords > nmaxwords)
        {
            nmaxwords = nwords;
        }

        // we can only use Mersenne reduction if they are all the same form
        allMersenne &= isMersenne;

        k++;
    }

    montyconst = vec_monty_alloc(nmaxwords);

    montyconst->NBLOCKS = nmaxwords / BLOCKWORDS;
    montyconst->NWORDS = nmaxwords;
    montyconst->MAXBITS = nmaxwords * DIGITBITS;

    // random base
    for (i = 0; i < VECLEN; i++)
    {
        vecbase[i] = lcg_rand_32_range(3, 0xffffffff, &fobj->pp1_obj.lcg_state);
    }

    for (i = 0; i < k; i++)
    {
        if (allMersenne != 0)
        {
            montyconst->isMersenne = isMersenne;
            montyconst->nbits = size_n;
            mpz_set(montyconst->nhat, N[i]);           // remember input N
            // do all math w.r.t the Mersenne number
            mpz_set_ui(N[i], 1);
            mpz_mul_2exp(N[i], N[i], size_n);
            if (isMersenne > 0)
            {
                mpz_sub_ui(N[i], N[i], isMersenne);
            }
            else
            {
                mpz_add_ui(N[i], N[i], 1);
            }
            broadcast_mpz_to_vec(montyconst->n, N[i]);
            broadcast_mpz_to_vec(montyconst->vnhat, montyconst->nhat);
            mpz_set_ui(r, 1);
            broadcast_mpz_to_vec(montyconst->one, r);
        }
        else
        {
            montyconst->isMersenne = allMersenne;
            montyconst->nbits = maxinbits;
            mpz_set_ui(r, 1);
            mpz_mul_2exp(r, r, montyconst->MAXBITS);
            mpz_invert(montyconst->nhat, N[i], r);
            mpz_sub(montyconst->nhat, r, montyconst->nhat);
            mpz_invert(montyconst->rhat, r, N[i]);
            insert_mpz_to_vec(montyconst->n, N[i], i);
            insert_mpz_to_vec(montyconst->r, r, i);
            insert_mpz_to_vec(montyconst->vrhat, montyconst->rhat, i);
            insert_mpz_to_vec(montyconst->vnhat, montyconst->nhat, i);
            mpz_tdiv_r(r, r, N[i]);
            insert_mpz_to_vec(montyconst->one, r, i);
        }

        montyconst->vrho[i] = mpz_get_ui(montyconst->nhat) & VEC_MAXDIGIT;
    }

    if (verbose > 1)
    {
        printf("P+1 has been configured with DIGITBITS = %u, VECLEN = %d, GMP_LIMB_BITS = %d\n",
            DIGITBITS, VECLEN, GMP_LIMB_BITS);

        printf("Choosing MAXBITS = %u, NWORDS = %d, NBLOCKS = %d based on max input size %d\n",
            montyconst->MAXBITS, montyconst->NWORDS, montyconst->NBLOCKS, maxinbits);
    }

    if (DIGITBITS == 52)
    {
        if (montyconst->isMersenne > 1)
        {
            vecmulmod_ptr = &vecmulmod52_mersenne;
            vecsqrmod_ptr = &vecsqrmod52_mersenne;
            vecaddmod_ptr = &vecaddmod52_mersenne;
            vecsubmod_ptr = &vecsubmod52_mersenne;
            vecaddsubmod_ptr = &vec_simul_addsub52_mersenne;
            if (verbose > 1)
            {
                printf("Using special pseudo-Mersenne mod for factor of: 2^%d-%d\n",
                    montyconst->nbits, montyconst->isMersenne);
            }
        }
        else if (montyconst->isMersenne > 0)
        {
            vecmulmod_ptr = &vecmulmod52_mersenne;
            vecsqrmod_ptr = &vecsqrmod52_mersenne;
            vecaddmod_ptr = &vecaddmod52_mersenne;
            vecsubmod_ptr = &vecsubmod52_mersenne;
            vecaddsubmod_ptr = &vec_simul_addsub52_mersenne;
            if (verbose > 1)
            {
                printf("Using special Mersenne mod for factor of: 2^%d-1\n", montyconst->nbits);
            }
        }
        else if (montyconst->isMersenne < 0)
        {
            vecmulmod_ptr = &vecmulmod52_mersenne;
            vecsqrmod_ptr = &vecsqrmod52_mersenne;
            vecaddmod_ptr = &vecaddmod52_mersenne;
            vecsubmod_ptr = &vecsubmod52_mersenne;
            vecaddsubmod_ptr = &vec_simul_addsub52_mersenne;
            if (verbose > 1)
            {
                printf("Using special Mersenne mod for factor of: 2^%d+1\n", montyconst->nbits);
            }
        }
        else
        {
            vecmulmod_ptr = &vecmulmod52;
            vecsqrmod_ptr = &vecsqrmod52;
            vecaddmod_ptr = &vecaddmod52;
            vecsubmod_ptr = &vecsubmod52;
            vecaddsubmod_ptr = &vec_simul_addsub52;
        }
    }
    else
    {
        if (montyconst->isMersenne)
        {
            vecmulmod_ptr = &vecmulmod_mersenne;
            vecsqrmod_ptr = &vecsqrmod_mersenne;
            vecaddmod_ptr = &vecaddmod_mersenne;
            vecsubmod_ptr = &vecsubmod_mersenne;
            vecaddsubmod_ptr = &vec_simul_addsub_mersenne;
            if (verbose > 1)
            {
                printf("Using special Mersenne mod for factor of: 2^%d-1\n", montyconst->nbits);
            }
        }
        else
        {
            vecmulmod_ptr = &vecmulmod;
            vecsqrmod_ptr = &vecsqrmod;
            vecaddmod_ptr = &vecaddmod;
            vecsubmod_ptr = &vecsubmod;
            vecaddsubmod_ptr = &vec_simul_addsub;
        }
    }


    // timing variables
    struct timeval stopt;	// stop time of this job
    struct timeval startt;	// start time of this job
    double t_time;

    gettimeofday(&startt, NULL);

    // put the random starting points in a vector and convert them.
    vecV = vecInit(nmaxwords);
    for (i = 0; i < k; i++)
    {
        mpz_set_ui(r, vecbase[i]);
        if (!allMersenne)
        {
            // into Monty rep
            mpz_mul_2exp(r, r, montyconst->MAXBITS);
            mpz_tdiv_r(r, r, N[i]);
        }
        else
        {
            mpz_tdiv_r(r, r, N[i]);
        }

        insert_mpz_to_vec(vecV, r, i);
    }

    // stage 1
    i = 0;
    while (ppm1_primes[i] < STAGE1_MAX)
    {
        //compute more primes if we need to.  it is more efficient to 
        //get a large chunk, even if we won't use them all.
        if ((uint64_t)i >= ppm1_nump)
        {
            rangemin += PRIME_RANGE;
            rangemax = MIN(STAGE1_MAX + 1000, rangemin + PRIME_RANGE);
            ppm1_primes = soe_wrapper(sdata, rangemin, rangemax, 0, &ppm1_nump, 0, 0);
            ppm1_minp = ppm1_primes[0];
            ppm1_maxp = ppm1_primes[ppm1_nump - 1];

            if (verbose > 1)
            {
                printf("pp1: found %lu primes in range [%lu : %lu]\n", ppm1_nump, rangemin, rangemax);
            }
            i = 0;
        }

        uint64_t q;
        int e;
        q = ppm1_primes[i];
        e = floor(log(STAGE1_MAX) / log(q));

        if (e == 1)
            break;

        q = (uint64_t)pow(q, e);
        vecLucasV(montyconst, q, vecV);
        i++;
    }

    vec_bignum_t* vecT1 = vecInit(nmaxwords);
    vec_bignum_t* vecT2 = vecInit(nmaxwords);
    vec_bignum_t* vecB = vecInit(nmaxwords);
    vec_bignum_t* vecC = vecInit(nmaxwords);

    while (ppm1_primes[i] < STAGE1_MAX)
    {
        //compute more primes if we need to.  it is more efficient to 
        //get a large chunk, even if we won't use them all.
        if ((uint64_t)i >= ppm1_nump)
        {
            rangemin += PRIME_RANGE;
            rangemax = MIN(STAGE1_MAX + 1000, rangemin + PRIME_RANGE);
            ppm1_primes = soe_wrapper(sdata, rangemin, rangemax, 0, &ppm1_nump, 0, 0);
            ppm1_minp = ppm1_primes[0];
            ppm1_maxp = ppm1_primes[ppm1_nump - 1];

            if (verbose > 1)
            {
                printf("pp1: found %lu primes in range [%lu : %lu]\n", ppm1_nump, rangemin, rangemax);
            }
            i = 0;
        }

        uint64_t q;
        q = ppm1_primes[i];

        pp1_prac(vecV, q, vecB, vecC, vecT1, vecT2, montyconst);
        i++;

        if ((fobj->VFLAG > 0) && ((i & 131071) == 0))
        {
            printf("pp1: @ p = %lu\r", q);
        }
    }
    if ((fobj->VFLAG > 0) && (i > 131071))
    {
        printf("\n");
    }

    vecFree(vecB);
    vecFree(vecC);
    vecFree(vecT1);
    vecFree(vecT2);

    if (fobj->ecm_obj.save_b1)
    {
        FILE *save = fopen("save_vpp1_b1.txt", "a");

        if (save != NULL)
        {
            mpz_t tmpV;
            mpz_init(tmpV);
            for (i = 0; i < VECLEN; i++)
            {
                extract_bignum_from_vec_to_mpz(tmpV, vecV, i, nmaxwords);

                // take out of Monty
                extract_bignum_from_vec_to_mpz(r, montyconst->vnhat, i, nmaxwords);
                mpz_mul(r, r, tmpV);
                mpz_tdiv_r_2exp(r, r, montyconst->MAXBITS);
                mpz_mul(r, r, N[i]);
                mpz_add(r, r, tmpV);
                mpz_tdiv_q_2exp(r, r, montyconst->MAXBITS);
                if (mpz_cmp(r, N[i]) > 0)
                {
                    mpz_sub(r, r, N[i]);
                }

                fprintf(save, "METHOD=PP1; B1=%"PRIu64"; ", STAGE1_MAX);
                gmp_fprintf(save, "N=0x%Zx; ", N[i]);
                gmp_fprintf(save, "X=0x%Zx; PROGRAM=AVX-PP1;\n", r);
            }
            mpz_clear(tmpV);
            fclose(save);
        }
    }

    // extract and check for factors
    vecCopy(vecV, montyconst->mtmp1);
    vecsubmod_ptr(montyconst->mtmp1, montyconst->one, montyconst->mtmp1, montyconst);
    vecsubmod_ptr(montyconst->mtmp1, montyconst->one, montyconst->mtmp1, montyconst);
    int found[VECLEN];
    memset(found, 0, VECLEN * sizeof(int));

    for (i = 0; i < k; i++)
    {
        extract_bignum_from_vec_to_mpz(r, montyconst->mtmp1, i, nmaxwords);
        mpz_gcd(g, r, N[i]);

        if ((mpz_cmp_ui(g, 1) > 0) && (mpz_cmp(g, N[i]) < 0))
        {
            char* s = mpz_get_str(NULL, 10, N[i]);
            logprint_oc(fobj->flogname, "a", "AVX-PP1 with B1 = %lu on %s\n", STAGE1_MAX, s);
            free(s);

            //check if the factor is prime
            if (is_mpz_prp(g, fobj->NUM_WITNESSES))
            {
                gmp_printf("pp1: found prp%d factor = %Zd in p+1 stage 1, lane %d\n",
                    gmp_base10(g), g, i);

                s = mpz_get_str(NULL, 10, g);
                logprint_oc(fobj->flogname, "a", "found prp%d = %s\n",
                    gmp_base10(g), s);
                free(s);
            }
            else
            {
                gmp_printf("pp1: found c%d factor = %Zd in p+1 stage 1, lane %d\n",
                    gmp_base10(g), g, i);

                s = mpz_get_str(NULL, 10, g);
                logprint_oc(fobj->flogname, "a", "found c%d = %s\n",
                    gmp_base10(g), s);
                free(s);
            }

            mpz_tdiv_q(N[i], N[i], g);
            s = mpz_get_str(NULL, 10, N[i]);
            if (is_mpz_prp(N[i], fobj->NUM_WITNESSES))
            {
                logprint_oc(fobj->flogname, "a", "cofactor is prp%d = %s\n",
                    gmp_base10(N[i]), s);
            }
            else
            {
                logprint_oc(fobj->flogname, "a", "cofactor is c%d = %s\n",
                    gmp_base10(N[i]), s);
            }
            free(s);
            found[i] = 1;
        }
    }

    gettimeofday(&stopt, NULL);
    t_time = ytools_difftime(&startt, &stopt);
    if (verbose > 0)
    {
        printf("pp1: Stage1 took %1.4f seconds.\n", t_time);
    }

    // use gmp-ecm internal stage 2.
    // P+/-1 is not like ecm where we can just run more curves
    // to increase the probability of success.  So a vastly inferior
    // stage 2 (standard prime pairing continuation) doesn't make
    // much sense when we have gmp-ecm's stage 2 available.  The
    // greatly increased B2 makes up for not running vectorized.
    int status;
    ecm_params params;
    ecm_init(params);
    params->method = ECM_PP1;

    for (i = 0; i < k; i++)
    {
        if (found[i])
            continue;

        if ((STAGE1_MAX > PRIME_RANGE) && (fobj->VFLAG > 0))
        {
            printf("pp1: gmp-ecm stage2 @ lane %d\n", i);
        }

        params->B1done = STAGE1_MAX;
        if (fobj->VFLAG >= 3)
            params->verbose = fobj->VFLAG - 2;

        // grab this lane and take out of Montgomery rep.
        extract_bignum_from_vec_to_mpz(params->x, vecV, i, nmaxwords);
        extract_bignum_from_vec_to_mpz(g, montyconst->vnhat, i, nmaxwords);

        mpz_mul(r, params->x, g);
        mpz_tdiv_r_2exp(r, r, montyconst->MAXBITS);
        mpz_mul(r, r, N[i]);
        mpz_add(r, r, params->x);
        mpz_tdiv_q_2exp(r, r, montyconst->MAXBITS);
        if (mpz_cmp(r, N[i]) > 0)
        {
            mpz_sub(r, r, N[i]);
        }
        mpz_set(params->x, r);

        if (fobj->pp1_obj.stg2_is_default == 0)
        {
            //not default, tell gmp-ecm to use the requested B2
            //printf("using requested B2 value\n");
            uint64_2gmp(fobj->pp1_obj.B2, params->B2);
        }

        status = ecm_factor(g, N[i], STAGE1_MAX, params);

        if ((mpz_cmp_ui(g, 1) > 0) && (mpz_cmp(g, N[i]) < 0))
        {
            char* s = mpz_get_str(NULL, 10, N[i]);
            logprint_oc(fobj->flogname, "a", "AVX-PP1 with B1 = %lu on %s\n", STAGE1_MAX, s);
            free(s);

            //check if the factor is prime
            if (is_mpz_prp(g, fobj->NUM_WITNESSES))
            {
                gmp_printf("pp1: found prp%d factor = %Zd in p+1 stage 2, lane %d\n",
                    gmp_base10(g), g, i);

                s = mpz_get_str(NULL, 10, g);
                logprint_oc(fobj->flogname, "a", "found prp%d = %s\n",
                    gmp_base10(g), s);
                free(s);
            }
            else
            {
                gmp_printf("pp1: found c%d factor = %Zd in p+1 stage 2, lane %d\n",
                    gmp_base10(g), g, i);

                s = mpz_get_str(NULL, 10, g);
                logprint_oc(fobj->flogname, "a", "found c%d = %s\n",
                    gmp_base10(g), s);
                free(s);
            }

            mpz_tdiv_q(N[i], N[i], g);
            s = mpz_get_str(NULL, 10, N[i]);
            if (is_mpz_prp(N[i], fobj->NUM_WITNESSES))
            {
                logprint_oc(fobj->flogname, "a", "cofactor is prp%d = %s\n",
                    gmp_base10(N[i]), s);
            }
            else
            {
                logprint_oc(fobj->flogname, "a", "cofactor is c%d = %s\n",
                    gmp_base10(N[i]), s);
            }
            free(s);
        }
    }
    ecm_clear(params);

    gettimeofday(&stopt, NULL);
    t_time = ytools_difftime(&startt, &stopt);

    if (verbose > 0)
    {
        printf("pp1: Process took %1.4f seconds.\n", t_time);
    }

    soe_finalize(sdata);
    vec_monty_free(montyconst);
    for (i = 0; i < k; i++)
    {
        mpz_clear(N[i]);
    }
    mpz_clear(r);
    mpz_clear(g);
    free(ppm1_primes);
    return;
}

void vecPM1(fact_obj_t* fobj)
{
    // run Pollard's P-1 factoring method on 8 inputs simultaneously.
    // The inputs should be about the same size so that all can use 
    // the same multplier lengths.  We use the longest one necessary
    // for all inputs.
    int i;
    vec_bignum_t* vecV;
    mpz_t N[VECLEN], r, g, e;
    int verbose = fobj->VFLAG;
    uint32_t maxbits, nwords, nblocks, nmaxwords = 0, maxinbits = 0;
    vec_monty_t* montyconst;
    // these track the range over which we currently have a prime list.
    uint64_t rangemin;
    uint64_t rangemax;
    uint64_t* ppm1_primes;
    uint64_t ppm1_nump;
    uint64_t ppm1_minp;
    uint64_t ppm1_maxp;
    uint64_t STAGE1_MAX = fobj->pm1_obj.B1;
    uint64_t STAGE2_MAX = 100 * fobj->pm1_obj.B1;
    uint64_t PRIME_RANGE = 10000000;
    uint32_t size_n;
    int isMersenne = 0, allMersenne = 1, forceNoMersenne = 0;
    uint64_t vecbase[VECLEN];

    rangemin = 0;
    rangemax = MIN(STAGE1_MAX, (uint64_t)PRIME_RANGE);

    soe_staticdata_t* sdata = soe_init(0, 1, 32768);
    ppm1_primes = soe_wrapper(sdata, rangemin, rangemax, 0, &ppm1_nump, 0, 0);
    ppm1_minp = ppm1_primes[0];
    ppm1_maxp = ppm1_primes[ppm1_nump - 1];

    if (verbose > 1)
    {
        printf("pm1: found %lu primes in range [%lu : %lu]\n", ppm1_nump, rangemin, rangemax);
    }

    mpz_init(r);
    mpz_init(g);
    mpz_init(e);

    int k = 0;
    while (k < fobj->pm1_obj.vecnum)
    {
        mpz_init(N[k]);
        mpz_set(N[k], fobj->pm1_obj.vecn[k]);

        // check for Mersenne inputs
        size_n = mpz_sizeinbase(N[k], 2);

        if (size_n > maxinbits)
        {
            maxinbits = size_n;
        }
        
        /*
        // how do we treat cases where some inputs have special form and others don't?
        // or if some special-form numbers are cheaper to test without the special form
        // and others aren't?
        // we could treat all as non-special.
        // we could treat all as special (if supported).
        // we could abort.  
        // I think it's best to warn of the odd-case and continue as all-special 
        // if possible.  if mixed inputs, treat all as non-special.
        for (i = size_n; i < 2048; i++)
        {
            mpz_set_ui(r, 1);
            mpz_mul_2exp(r, r, i);
            mpz_sub_ui(r, r, 1);
            mpz_mod(g, r, N[k]);
            if (mpz_cmp_ui(g, 0) == 0)
            {
                size_n = i;
                isMersenne = 1;
                break;
            }

            mpz_set_ui(r, 1);
            mpz_mul_2exp(r, r, i);
            mpz_add_ui(r, r, 1);
            mpz_mod(g, r, N[k]);
            if (mpz_cmp_ui(g, 0) == 0)
            {
                size_n = i;
                isMersenne = -1;
                break;
            }

            // detect pseudo-Mersennes
            mpz_set_ui(r, 1);
            mpz_mul_2exp(r, r, i);
            mpz_mod(g, r, N[k]);
            if (mpz_sizeinbase(g, 2) < DIGITBITS)
            {
                size_n = i;
                isMersenne = mpz_get_ui(g);
                break;
            }
        }
        */

        // if the input is Mersenne and still contains algebraic factors, remove them.
        if (0) //(abs(isMersenne) == 1)
        {
#ifdef USE_NFS
            snfs_t poly;

            snfs_init(&poly);
            poly.form_type = SNFS_BRENT;
            mpz_set_ui(poly.base1, 2);
            poly.coeff1 = 1;
            poly.coeff2 = -isMersenne;
            poly.exp1 = size_n;
            find_primitive_factor(&poly, fobj->primes, fobj->num_p, verbose);

            mpz_tdiv_q(g, N[k], poly.primitive);
            mpz_gcd(N[k], N[k], poly.primitive);

            if (mpz_cmp_ui(g, 1) > 0)
            {
                add_to_factor_list(fobj->factors, g,
                    fobj->VFLAG, fobj->NUM_WITNESSES);

                int fid = fobj->factors->num_factors - 1;

                fobj->factors->factors[fid].tid = 0;
                fobj->factors->factors[fid].curve_num = 0;
                fobj->factors->factors[fid].sigma = 0;
                fobj->factors->factors[fid].vid = 0;
                fobj->factors->factors[fid].method = 0;
            }

            snfs_clear(&poly);
#endif
        }
        

        // compute NBLOCKS if using the actual size of the input (non-Mersenne)
        if (DIGITBITS == 52)
        {
            maxbits = 208;
            while (maxbits <= mpz_sizeinbase(N[k], 2))
            {
                maxbits += 208;
            }
        }
        else
        {
            maxbits = 128;
            while (maxbits <= mpz_sizeinbase(N[k], 2))
            {
                maxbits += 128;
            }
        }

        nwords = maxbits / DIGITBITS;
        nblocks = nwords / BLOCKWORDS;

        // and compute NBLOCKS if using Mersenne mod
        if (DIGITBITS == 52)
        {
            maxbits = 208;
            while (maxbits <= size_n)
            {
                maxbits += 208;
            }
        }
        else
        {
            maxbits = 128;
            while (maxbits <= size_n)
            {
                maxbits += 128;
            }
        }

        if (verbose > 0)
        {
            gmp_printf("pm1: commencing parallel P-1 @ B1=%lu on %Zd\n", STAGE1_MAX, N[k]);
        }

        if ((double)nwords / ((double)maxbits / (double)DIGITBITS) < 0.7)
        {
            if (verbose > 1)
            {
                printf("Mersenne input 2^%d - 1 determined to be faster by REDC\n",
                    size_n);
            }
            forceNoMersenne = 1;
        }
        else
        {
            nwords = maxbits / DIGITBITS;
            nblocks = nwords / BLOCKWORDS;
        }

        if (forceNoMersenne)
        {
            isMersenne = 0;
            size_n = mpz_sizeinbase(N[k], 2);
        }

        if (nwords > nmaxwords)
        {
            nmaxwords = nwords;
        }

        // we can only use Mersenne reduction if they are all the same form
        allMersenne &= isMersenne;

        k++;
    }

    montyconst = vec_monty_alloc(nmaxwords);

    //montyconst->NBLOCKS = nblocks;
    //montyconst->NWORDS = nwords;
    //montyconst->MAXBITS = maxbits;

    montyconst->NBLOCKS = nmaxwords / BLOCKWORDS;
    montyconst->NWORDS = nmaxwords;
    montyconst->MAXBITS = nmaxwords * DIGITBITS;

    // pm1 base
    for (i = 0; i < VECLEN; i++)
    {
        vecbase[i] = 3;
    }

    for (i = 0; i < k; i++)
    {
        if (allMersenne != 0)
        {
            montyconst->isMersenne = isMersenne;
            montyconst->nbits = size_n;
            mpz_set(montyconst->nhat, N[i]);           // remember input N
            // do all math w.r.t the Mersenne number
            mpz_set_ui(N[i], 1);
            mpz_mul_2exp(N[i], N[i], size_n);
            if (isMersenne > 0)
            {
                mpz_sub_ui(N[i], N[i], isMersenne);
            }
            else
            {
                mpz_add_ui(N[i], N[i], 1);
            }
            broadcast_mpz_to_vec(montyconst->n, N[i]);
            broadcast_mpz_to_vec(montyconst->vnhat, montyconst->nhat);
            mpz_set_ui(r, 1);
            broadcast_mpz_to_vec(montyconst->one, r);
        }
        else
        {
            montyconst->isMersenne = allMersenne;
            montyconst->nbits = maxinbits; // mpz_sizeinbase(N, 2);
            mpz_set_ui(r, 1);
            mpz_mul_2exp(r, r, montyconst->MAXBITS);
            mpz_invert(montyconst->nhat, N[i], r);
            mpz_sub(montyconst->nhat, r, montyconst->nhat);
            mpz_invert(montyconst->rhat, r, N[i]);
            insert_mpz_to_vec(montyconst->n, N[i], i);
            insert_mpz_to_vec(montyconst->r, r, i);
            insert_mpz_to_vec(montyconst->vrhat, montyconst->rhat, i);
            insert_mpz_to_vec(montyconst->vnhat, montyconst->nhat, i);
            mpz_tdiv_r(r, r, N[i]);
            insert_mpz_to_vec(montyconst->one, r, i);
        }

        montyconst->vrho[i] = mpz_get_ui(montyconst->nhat) & VEC_MAXDIGIT;
    }

    if (verbose > 1)
    {
        printf("P-1 has been configured with DIGITBITS = %u, VECLEN = %d, GMP_LIMB_BITS = %d\n",
            DIGITBITS, VECLEN, GMP_LIMB_BITS);

        printf("Choosing MAXBITS = %u, NWORDS = %d, NBLOCKS = %d based on max input size %d\n",
            montyconst->MAXBITS, montyconst->NWORDS, montyconst->NBLOCKS, maxinbits);
    }

    if (DIGITBITS == 52)
    {
        if (montyconst->isMersenne > 1)
        {
            vecmulmod_ptr = &vecmulmod52_mersenne;
            vecsqrmod_ptr = &vecsqrmod52_mersenne;
            vecaddmod_ptr = &vecaddmod52_mersenne;
            vecsubmod_ptr = &vecsubmod52_mersenne;
            vecaddsubmod_ptr = &vec_simul_addsub52_mersenne;
            if (verbose > 1)
            {
                printf("Using special pseudo-Mersenne mod for factor of: 2^%d-%d\n",
                    montyconst->nbits, montyconst->isMersenne);
            }
        }
        else if (montyconst->isMersenne > 0)
        {
            vecmulmod_ptr = &vecmulmod52_mersenne;
            vecsqrmod_ptr = &vecsqrmod52_mersenne;
            vecaddmod_ptr = &vecaddmod52_mersenne;
            vecsubmod_ptr = &vecsubmod52_mersenne;
            vecaddsubmod_ptr = &vec_simul_addsub52_mersenne;
            if (verbose > 1)
            {
                printf("Using special Mersenne mod for factor of: 2^%d-1\n", montyconst->nbits);
            }
        }
        else if (montyconst->isMersenne < 0)
        {
            vecmulmod_ptr = &vecmulmod52_mersenne;
            vecsqrmod_ptr = &vecsqrmod52_mersenne;
            vecaddmod_ptr = &vecaddmod52_mersenne;
            vecsubmod_ptr = &vecsubmod52_mersenne;
            vecaddsubmod_ptr = &vec_simul_addsub52_mersenne;
            if (verbose > 1)
            {
                printf("Using special Mersenne mod for factor of: 2^%d+1\n", montyconst->nbits);
            }
        }
        else
        {
            vecmulmod_ptr = &vecmulmod52;
            vecsqrmod_ptr = &vecsqrmod52;
            vecaddmod_ptr = &vecaddmod52;
            vecsubmod_ptr = &vecsubmod52;
            vecaddsubmod_ptr = &vec_simul_addsub52;
        }
    }
    else
    {
        if (montyconst->isMersenne)
        {
            vecmulmod_ptr = &vecmulmod_mersenne;
            vecsqrmod_ptr = &vecsqrmod_mersenne;
            vecaddmod_ptr = &vecaddmod_mersenne;
            vecsubmod_ptr = &vecsubmod_mersenne;
            vecaddsubmod_ptr = &vec_simul_addsub_mersenne;
            if (verbose > 1)
            {
                printf("Using special Mersenne mod for factor of: 2^%d-1\n", montyconst->nbits);
            }
        }
        else
        {
            vecmulmod_ptr = &vecmulmod;
            vecsqrmod_ptr = &vecsqrmod;
            vecaddmod_ptr = &vecaddmod;
            vecsubmod_ptr = &vecsubmod;
            vecaddsubmod_ptr = &vec_simul_addsub;
        }
    }


    // timing variables
    struct timeval stopt;	// stop time of this job
    struct timeval startt;	// start time of this job
    double t_time;

    gettimeofday(&startt, NULL);

    // put the starting points in a vector and convert them.
    vecV = vecInit(nmaxwords);
    for (i = 0; i < k; i++)
    {
        mpz_set_ui(r, vecbase[i]);
        if (!allMersenne)
        {
            // into Monty rep
            mpz_mul_2exp(r, r, montyconst->MAXBITS);
            mpz_tdiv_r(r, r, N[i]);
        }
        else
        {
            mpz_tdiv_r(r, r, N[i]);
        }

        insert_mpz_to_vec(vecV, r, i);
    }

    // stage 1
    vec_bignum_t* A = vecInit(montyconst->NWORDS);
    vec_bignum_t* T = vecInit(montyconst->NWORDS);

    while (rangemin < STAGE1_MAX)
    {
        struct timeval start1;	// start time of this job
        gettimeofday(&start1, NULL);

        compute_s(e, ppm1_primes, ppm1_nump, rangemax);

        gettimeofday(&stopt, NULL);
        t_time = ytools_difftime(&start1, &stopt);
        printf("pm1: primes product took %1.2f sec, exponent has %u bits\n", 
            t_time, mpz_sizeinbase(e, 2));

        //vecmodexp(A, vecV, e, montyconst->n, T, montyconst->one, montyconst);
        vecslidingmodexp(A, vecV, e, montyconst->n, T, montyconst->one, montyconst);
        vecCopy(A, vecV);

        rangemin = rangemax;
        rangemax = MIN(STAGE1_MAX, rangemin + PRIME_RANGE);

        if (rangemin < STAGE1_MAX)
        {
            ppm1_primes = soe_wrapper(sdata, rangemin, rangemax, 0, &ppm1_nump, 0, 0);
            ppm1_minp = ppm1_primes[0];
            ppm1_maxp = ppm1_primes[ppm1_nump - 1];

            if (verbose > 1)
            {
                printf("pm1: found %lu primes in range [%lu : %lu]\n", ppm1_nump, rangemin, rangemax);
            }
        }
    }

    vecFree(A);
    vecFree(T);

    // extract and check for factors
    vecsubmod_ptr(vecV, montyconst->one, montyconst->mtmp1, montyconst);
    int found[VECLEN];

    for (i = 0; i < k; i++)
    {
        found[i] = 0;
        extract_bignum_from_vec_to_mpz(r, montyconst->mtmp1, i, nmaxwords);
        mpz_gcd(g, r, N[i]);

        if ((mpz_cmp_ui(g, 1) > 0) && (mpz_cmp(g, N[i]) < 0))
        {
            char* s = mpz_get_str(NULL, 10, N[i]);
            logprint_oc(fobj->flogname, "a", "AVX-PM1 with B1 = %lu on %s\n", STAGE1_MAX, s);
            free(s);

            //check if the factor is prime
            if (is_mpz_prp(g, fobj->NUM_WITNESSES))
            {
                gmp_printf("pm1: found prp%d factor = %Zd in p-1 stage 1, lane %d\n",
                    gmp_base10(g), g, i);

                s = mpz_get_str(NULL, 10, g);
                logprint_oc(fobj->flogname, "a", "found prp%d = %s\n",
                    gmp_base10(g), s);
                free(s);
            }
            else
            {
                gmp_printf("pm1: found c%d factor = %Zd in p-1 stage 1, lane %d\n",
                    gmp_base10(g), g, i);

                s = mpz_get_str(NULL, 10, g);
                logprint_oc(fobj->flogname, "a", "found c%d = %s\n",
                    gmp_base10(g), s);
                free(s);
            }

            mpz_tdiv_q(N[i], N[i], g);
            s = mpz_get_str(NULL, 10, N[i]);
            if (is_mpz_prp(N[i], fobj->NUM_WITNESSES))
            {
                logprint_oc(fobj->flogname, "a", "cofactor is prp%d = %s\n",
                    gmp_base10(N[i]), s);
            }
            else
            {
                logprint_oc(fobj->flogname, "a", "cofactor is c%d = %s\n",
                    gmp_base10(N[i]), s);
            }
            free(s);
            found[i] = 1;
        }
    }

    gettimeofday(&stopt, NULL);
    t_time = ytools_difftime(&startt, &stopt);
    if (verbose > 0)
    {
        printf("pm1: Stage1 took %1.4f seconds.\n", t_time);
    }

    if (fobj->ecm_obj.save_b1)
    {
        FILE* save = fopen("save_vpm1_b1.txt", "a");

        if (save != NULL)
        {
            mpz_t tmpV;
            mpz_init(tmpV);
            for (i = 0; i < VECLEN; i++)
            {
                extract_bignum_from_vec_to_mpz(tmpV, vecV, i, nmaxwords);

                // take out of Monty
                extract_bignum_from_vec_to_mpz(r, montyconst->vnhat, i, nmaxwords);
                mpz_mul(r, r, tmpV);
                mpz_tdiv_r_2exp(r, r, montyconst->MAXBITS);
                mpz_mul(r, r, N[i]);
                mpz_add(r, r, tmpV);
                mpz_tdiv_q_2exp(r, r, montyconst->MAXBITS);
                if (mpz_cmp(r, N[i]) > 0)
                {
                    mpz_sub(r, r, N[i]);
                }

                fprintf(save, "METHOD=PM1; B1=%"PRIu64"; ", STAGE1_MAX);
                gmp_fprintf(save, "N=0x%Zx; ", N[i]);
                gmp_fprintf(save, "X=0x%Zx; PROGRAM=AVX-PP1;\n", r);
            }
            mpz_clear(tmpV);
            fclose(save);
        }
    }

    // use gmp-ecm internal stage 2.
    // P+/-1 is not like ecm where we can just run more curves
    // to increase the probability of success.  So a vastly inferior
    // stage 2 (standard prime pairing continuation) doesn't make
    // much sense when we have gmp-ecm's stage 2 available.  The
    // greatly increased B2 makes up for not running vectorized.
    int status;
    ecm_params params;
    ecm_init(params);
    params->method = ECM_PM1;

    for (i = 0; i < k; i++)
    {
        if (found[i])
            continue;

        if ((STAGE1_MAX > PRIME_RANGE) && (fobj->VFLAG > 0))
        {
            printf("pm1: gmp-ecm stage2 @ lane %d\n", i);
        }

        params->B1done = STAGE1_MAX;
        if (fobj->VFLAG >= 3)
            params->verbose = fobj->VFLAG - 2;

        // grab this lane and take out of Montgomery rep.
        extract_bignum_from_vec_to_mpz(params->x, vecV, i, montyconst->NWORDS);
        extract_bignum_from_vec_to_mpz(g, montyconst->vnhat, i, montyconst->NWORDS);

        mpz_mul(r, params->x, g);
        mpz_tdiv_r_2exp(r, r, montyconst->MAXBITS);
        mpz_mul(r, r, N[i]);
        mpz_add(r, r, params->x);
        mpz_tdiv_q_2exp(r, r, montyconst->MAXBITS);
        if (mpz_cmp(r, N[i]) > 0)
        {
            mpz_sub(r, r, N[i]);
        }
        mpz_set(params->x, r);

        if (fobj->pm1_obj.stg2_is_default == 0)
        {
            //not default, tell gmp-ecm to use the requested B2
            //printf("using requested B2 value\n");
            uint64_2gmp(fobj->pm1_obj.B2, params->B2);
        }

        status = ecm_factor(g, N[i], STAGE1_MAX, params);

        if ((mpz_cmp_ui(g, 1) > 0) && (mpz_cmp(g, N[i]) < 0))
        {
            char* s = mpz_get_str(NULL, 10, N[i]);
            logprint_oc(fobj->flogname, "a", "AVX-PM1 with B1 = %lu on %s\n", STAGE1_MAX, s);
            free(s);

            //check if the factor is prime
            if (is_mpz_prp(g, fobj->NUM_WITNESSES))
            {
                gmp_printf("pm1: found prp%d factor = %Zd in p-1 stage 2, lane %d\n",
                    gmp_base10(g), g, i);

                s = mpz_get_str(NULL, 10, g);
                logprint_oc(fobj->flogname, "a", "found prp%d = %s\n",
                    gmp_base10(g), s);
                free(s);
            }
            else
            {
                gmp_printf("pm1: found c%d factor = %Zd in p-1 stage 2, lane %d\n",
                    gmp_base10(g), g, i);

                s = mpz_get_str(NULL, 10, g);
                logprint_oc(fobj->flogname, "a", "found c%d = %s\n",
                    gmp_base10(g), s);
                free(s);
            }

            mpz_tdiv_q(N[i], N[i], g);
            s = mpz_get_str(NULL, 10, N[i]);
            if (is_mpz_prp(N[i], fobj->NUM_WITNESSES))
            {
                logprint_oc(fobj->flogname, "a", "cofactor is prp%d = %s\n",
                    gmp_base10(N[i]), s);
            }
            else
            {
                logprint_oc(fobj->flogname, "a", "cofactor is c%d = %s\n",
                    gmp_base10(N[i]), s);
            }
            free(s);
        }
    }
    ecm_clear(params);

    gettimeofday(&stopt, NULL);
    t_time = ytools_difftime(&startt, &stopt);

    if (verbose > 0)
    {
        printf("pm1: Process took %1.4f seconds.\n", t_time);
    }

    soe_finalize(sdata);
    vec_monty_free(montyconst);
    for (i = 0; i < k; i++)
    {
        mpz_clear(N[i]);
    }
    mpz_clear(r);
    mpz_clear(g);
    mpz_clear(e);
    free(ppm1_primes);
    return;
}

void vecPPM1_1N(fact_obj_t *fobj)
{
    // run William's P+1 factoring method with 8 random bases.
    // It is hoped that at least one of these will give a P+1 test
    // and the rest will give a P-1 test.  So we cover twice
    // the factor-space at the cost of duplicating several P-1/P+1 tests
    // and running slower than either P-1 or P+1 would take individually
    // using GMP-ECM.  Is this worth it?  Shrug.  Probably not.  But
    // the routines already existed; porting to vec arithmetic wasn't too
    // much work.
    int i, stg1i;
    vec_bignum_t *vecV;
    mpz_t N, r, g;
    int verbose = fobj->VFLAG;
    uint32_t maxbits, nwords, nblocks;
    vec_monty_t* montyconst;
    // these track the range over which we currently have a prime list.
    uint64_t rangemin;
    uint64_t rangemax;
    uint64_t* ppm1_primes;
    uint64_t ppm1_nump;
    uint64_t ppm1_minp;
    uint64_t ppm1_maxp;
    uint64_t STAGE1_MAX = fobj->pp1_obj.B1;
    uint64_t STAGE2_MAX = 100 * fobj->pp1_obj.B1;
    uint64_t PRIME_RANGE = 100000000;
    uint32_t size_n;
    int isMersenne = 0, forceNoMersenne = 0;
    uint64_t vecbase[VECLEN];

    rangemin = 0;
    rangemax = MIN(STAGE1_MAX + 1000, (uint64_t)PRIME_RANGE);

    soe_staticdata_t* sdata = soe_init(0, 1, 32768);
    ppm1_primes = soe_wrapper(sdata, rangemin, rangemax, 0, &ppm1_nump, 0, 0);
    ppm1_minp = ppm1_primes[0];
    ppm1_maxp = ppm1_primes[ppm1_nump - 1];

    

    if (verbose > 1)
    {
        printf("pp1: found %lu primes in range [%lu : %lu]\n", ppm1_nump, rangemin, rangemax);
    }

    mpz_init(N);
    mpz_init(r);
    mpz_init(g);

    mpz_set(N, fobj->pp1_obj.gmp_n);

    // check for Mersenne inputs
    size_n = mpz_sizeinbase(N, 2);

    for (i = size_n; i < 2048; i++)
    {
        mpz_set_ui(r, 1);
        mpz_mul_2exp(r, r, i);
        mpz_sub_ui(r, r, 1);
        mpz_mod(g, r, N);
        if (mpz_cmp_ui(g, 0) == 0)
        {
            size_n = i;
            isMersenne = 1;
            break;
        }

        mpz_set_ui(r, 1);
        mpz_mul_2exp(r, r, i);
        mpz_add_ui(r, r, 1);
        mpz_mod(g, r, N);
        if (mpz_cmp_ui(g, 0) == 0)
        {
            size_n = i;
            isMersenne = -1;
            break;
        }

        // detect pseudo-Mersennes
        mpz_set_ui(r, 1);
        mpz_mul_2exp(r, r, i);
        mpz_mod(g, r, N);
        if (mpz_sizeinbase(g, 2) < DIGITBITS)
        {
            size_n = i;
            isMersenne = mpz_get_ui(g);
            break;
        }
    }

    // if the input is Mersenne and still contains algebraic factors, remove them.
    if (abs(isMersenne) == 1)
    {
#ifdef USE_NFS
        snfs_t poly;

        snfs_init(&poly);
        poly.form_type = SNFS_BRENT;
        mpz_set_ui(poly.base1, 2);
        poly.coeff1 = 1;
        poly.coeff2 = -isMersenne;
        poly.exp1 = size_n;
        find_primitive_factor(&poly, fobj->primes, fobj->num_p, verbose);

        mpz_tdiv_q(g, N, poly.primitive);
        mpz_gcd(N, N, poly.primitive);

        if (mpz_cmp_ui(g, 1) > 0)
        {
            add_to_factor_list(fobj->factors, g,
                fobj->VFLAG, fobj->NUM_WITNESSES);

            int fid = fobj->factors->num_factors - 1;

            fobj->factors->factors[fid].tid = 0;
            fobj->factors->factors[fid].curve_num = 0;
            fobj->factors->factors[fid].sigma = 0;
            fobj->factors->factors[fid].vid = 0;
            fobj->factors->factors[fid].method = 0;
        }

        snfs_clear(&poly);
#endif
    }

    // force no Mersenne Mod, either specified by user or if the actual
    // input is smaller than the base Mersenne such that the arithmetic
    // becomes faster using fewer-digit Montgomery mod.
    // rough data for DIGITBITS 52 (Gold 6248 CPU @ 2.5 GHz)
    // NBLOCKS 6 MERSENNEMOD is faster than NBLOCKS 5 REDC ratio 0.83
    // NBLOCKS 5 MERSENNEMOD is faster than NBLOCKS 4 REDC ratio 0.8
    // NBLOCKS 4 MERSENNEMOD is faster than NBLOCKS 3 REDC ratio 0.75
    // NBLOCKS 3 MERSENNEMOD is slower than NBLOCKS 2 REDC ratio 0.66 (barely slower)

    // compute NBLOCKS if using the actual size of the input (non-Mersenne)
    if (DIGITBITS == 52)
    {
        maxbits = 208;
        while (maxbits <= mpz_sizeinbase(N, 2))
        {
            maxbits += 208;
        }
    }
    else
    {
        maxbits = 128;
        while (maxbits <= mpz_sizeinbase(N, 2))
        {
            maxbits += 128;
        }
    }

    nwords = maxbits / DIGITBITS;
    nblocks = nwords / BLOCKWORDS;

    // and compute NBLOCKS if using Mersenne mod
    if (DIGITBITS == 52)
    {
        maxbits = 208;
        while (maxbits <= size_n)
        {
            maxbits += 208;
        }
    }
    else
    {
        maxbits = 128;
        while (maxbits <= size_n)
        {
            maxbits += 128;
        }
    }

    if (verbose > 1)
    {
        gmp_printf("commencing parallel P+/-1 on %Zd\n", N);
    }

    if ((double)nwords / ((double)maxbits / (double)DIGITBITS) < 0.7)
    {
        if (verbose > 1)
        {
            printf("Mersenne input 2^%d - 1 determined to be faster by REDC\n",
                size_n);
        }
        forceNoMersenne = 1;
    }
    else
    {
        nwords = maxbits / DIGITBITS;
        nblocks = nwords / BLOCKWORDS;
    }

    if (forceNoMersenne)
    {
        isMersenne = 0;
        size_n = mpz_sizeinbase(N, 2);
    }

    montyconst = vec_monty_alloc(nwords);

    // random base
    for (i = 0; i < VECLEN; i++)
    {
        vecbase[i] = lcg_rand_32_range(3, 0xffffffff, &fobj->pp1_obj.lcg_state);
    }

    if (isMersenne != 0)
    {
        montyconst->isMersenne = isMersenne;
        montyconst->nbits = size_n;
        mpz_set(montyconst->nhat, N);           // remember input N
        // do all math w.r.t the Mersenne number
        mpz_set_ui(N, 1);
        mpz_mul_2exp(N, N, size_n);
        if (isMersenne > 0)
        {
            mpz_sub_ui(N, N, isMersenne);
        }
        else
        {
            mpz_add_ui(N, N, 1);
        }
        broadcast_mpz_to_vec(montyconst->n, N);
        broadcast_mpz_to_vec(montyconst->vnhat, montyconst->nhat);
        mpz_set_ui(r, 1);
        broadcast_mpz_to_vec(montyconst->one, r);
    }
    else
    {
        montyconst->isMersenne = 0;
        montyconst->nbits = mpz_sizeinbase(N, 2);
        mpz_set_ui(r, 1);
        mpz_mul_2exp(r, r, DIGITBITS* nwords);
        //gmp_printf("r = (1 << %d) = %Zd\n", DIGITBITS * NWORDS, r);
        mpz_invert(montyconst->nhat, N, r);
        mpz_sub(montyconst->nhat, r, montyconst->nhat);
        mpz_invert(montyconst->rhat, r, N);
        broadcast_mpz_to_vec(montyconst->n, N);
        broadcast_mpz_to_vec(montyconst->r, r);
        broadcast_mpz_to_vec(montyconst->vrhat, montyconst->rhat);
        broadcast_mpz_to_vec(montyconst->vnhat, montyconst->nhat);
        mpz_tdiv_r(r, r, N);
        broadcast_mpz_to_vec(montyconst->one, r);
    }

    montyconst->NBLOCKS = nblocks;
    montyconst->NWORDS = nwords;
    montyconst->MAXBITS = maxbits;

    for (i = 0; i < VECLEN; i++)
    {
        montyconst->vrho[i] = mpz_get_ui(montyconst->nhat) & VEC_MAXDIGIT;
    }

    if (verbose > 1)
    {
        printf("P+/-1 has been configured with DIGITBITS = %u, VECLEN = %d, GMP_LIMB_BITS = %d\n",
            DIGITBITS, VECLEN, GMP_LIMB_BITS);

        printf("Choosing MAXBITS = %u, NWORDS = %d, NBLOCKS = %d based on input size %d\n",
            maxbits, nwords, nblocks, montyconst->nbits);
    }

    if (DIGITBITS == 52)
    {
        if (montyconst->isMersenne > 1)
        {
            vecmulmod_ptr = &vecmulmod52_mersenne;
            vecsqrmod_ptr = &vecsqrmod52_mersenne;
            vecaddmod_ptr = &vecaddmod52_mersenne;
            vecsubmod_ptr = &vecsubmod52_mersenne;
            vecaddsubmod_ptr = &vec_simul_addsub52_mersenne;
            if (verbose > 1)
            {
                printf("Using special pseudo-Mersenne mod for factor of: 2^%d-%d\n",
                    montyconst->nbits, montyconst->isMersenne);
            }
        }
        else if (montyconst->isMersenne > 0)
        {
            vecmulmod_ptr = &vecmulmod52_mersenne;
            vecsqrmod_ptr = &vecsqrmod52_mersenne;
            vecaddmod_ptr = &vecaddmod52_mersenne;
            vecsubmod_ptr = &vecsubmod52_mersenne;
            vecaddsubmod_ptr = &vec_simul_addsub52_mersenne;
            if (verbose > 1)
            {
                printf("Using special Mersenne mod for factor of: 2^%d-1\n", montyconst->nbits);
            }
        }
        else if (montyconst->isMersenne < 0)
        {
            vecmulmod_ptr = &vecmulmod52_mersenne;
            vecsqrmod_ptr = &vecsqrmod52_mersenne;
            vecaddmod_ptr = &vecaddmod52_mersenne;
            vecsubmod_ptr = &vecsubmod52_mersenne;
            vecaddsubmod_ptr = &vec_simul_addsub52_mersenne;
            if (verbose > 1)
            {
                printf("Using special Mersenne mod for factor of: 2^%d+1\n", montyconst->nbits);
            }
        }
        else
        {
            vecmulmod_ptr = &vecmulmod52;
            vecsqrmod_ptr = &vecsqrmod52;
            vecaddmod_ptr = &vecaddmod52;
            vecsubmod_ptr = &vecsubmod52;
            vecaddsubmod_ptr = &vec_simul_addsub52;
        }
    }
    else
    {
        if (montyconst->isMersenne)
        {
            vecmulmod_ptr = &vecmulmod_mersenne;
            vecsqrmod_ptr = &vecsqrmod_mersenne;
            vecaddmod_ptr = &vecaddmod_mersenne;
            vecsubmod_ptr = &vecsubmod_mersenne;
            vecaddsubmod_ptr = &vec_simul_addsub_mersenne;
            if (verbose > 1)
            {
                printf("Using special Mersenne mod for factor of: 2^%d-1\n", montyconst->nbits);
            }
        }
        else
        {
            vecmulmod_ptr = &vecmulmod;
            vecsqrmod_ptr = &vecsqrmod;
            vecaddmod_ptr = &vecaddmod;
            vecsubmod_ptr = &vecsubmod;
            vecaddsubmod_ptr = &vec_simul_addsub;
        }
    }


    // timing variables
    struct timeval stopt;	// stop time of this job
    struct timeval startt;	// start time of this job
    double t_time;

    gettimeofday(&startt, NULL);

    // put the random starting points in a vector and convert them.
    vecV = vecInit(nwords);
    for (i = 0; i < VECLEN; i++)
    {
        mpz_set_ui(r, vecbase[i]);
        if (!isMersenne)
        {
            // into Monty rep
            mpz_mul_2exp(r, r, DIGITBITS * montyconst->NWORDS);
            mpz_tdiv_r(r, r, N);
        }
        else
        {
            mpz_tdiv_r(r, r, N);
        }

        insert_mpz_to_vec(vecV, r, i);
    }

    // stage 1
    i = 0;
    while (ppm1_primes[i] < STAGE1_MAX)
    {
        //compute more primes if we need to.  it is more efficient to 
        //get a large chunk, even if we won't use them all.
        if ((uint64_t)i >= ppm1_nump)
        {
            rangemin += PRIME_RANGE;
            rangemax = MIN(STAGE1_MAX + 1000, rangemin + PRIME_RANGE);
            ppm1_primes = soe_wrapper(sdata, rangemin, rangemax, 0, &ppm1_nump, 0, 0);
            ppm1_minp = ppm1_primes[0];
            ppm1_maxp = ppm1_primes[ppm1_nump - 1];

            if (verbose > 1)
            {
                printf("pp1: found %lu primes in range [%lu : %lu]\n", ppm1_nump, rangemin, rangemax);
            }
            i = 0;
        }

        uint64_t q;
        int e;
        q = ppm1_primes[i];
        e = floor(log(STAGE1_MAX) / log(q));

        if (e == 1)
            break;

        q = (uint64_t)pow(q, e);
        vecLucasV(montyconst, q, vecV);
        i++;
    }

    vec_bignum_t* vecT1 = vecInit(nwords);
    vec_bignum_t* vecT2 = vecInit(nwords);
    vec_bignum_t* vecB = vecInit(nwords);
    vec_bignum_t* vecC = vecInit(nwords);

    while (ppm1_primes[i] < STAGE1_MAX)
    {
        //compute more primes if we need to.  it is more efficient to 
        //get a large chunk, even if we won't use them all.
        if ((uint64_t)i >= ppm1_nump)
        {
            rangemin += PRIME_RANGE;
            rangemax = MIN(STAGE1_MAX + 1000, rangemin + PRIME_RANGE);
            ppm1_primes = soe_wrapper(sdata, rangemin, rangemax, 0, &ppm1_nump, 0, 0);
            ppm1_minp = ppm1_primes[0];
            ppm1_maxp = ppm1_primes[ppm1_nump - 1];

            if (verbose > 1)
            {
                printf("pp1: found %lu primes in range [%lu : %lu]\n", ppm1_nump, rangemin, rangemax);
            }
            i = 0;
        }

        uint64_t q;
        q = ppm1_primes[i];

        pp1_prac(vecV, q, vecB, vecC, vecT1, vecT2, montyconst);
        i++;

        if ((fobj->VFLAG > 0) && ((i & 1048575) == 0))
        {
            printf("pp1: @ p = %lu\r", q);
        }
    }
    if ((fobj->VFLAG > 0) && (i > 1048575))
    {
        printf("\n");
    }

    vecFree(vecB);
    vecFree(vecC);
    vecFree(vecT1);
    vecFree(vecT2);

    // extract and check for factors
    vecCopy(vecV, montyconst->mtmp1);
    vecsubmod_ptr(montyconst->mtmp1, montyconst->one, montyconst->mtmp1, montyconst);
    vecsubmod_ptr(montyconst->mtmp1, montyconst->one, montyconst->mtmp1, montyconst);
    int found = 0;
    for (i = 0; i < VECLEN; i++)
    {
        if ((STAGE1_MAX > PRIME_RANGE) && (fobj->VFLAG > 0))
        {
            printf("pp1: gmp-ecm stage2 @ lane %d\n", i);
        }

        extract_bignum_from_vec_to_mpz(r, montyconst->mtmp1, i, nwords);
        mpz_gcd(g, r, fobj->pp1_obj.gmp_n);

        if ((mpz_cmp_ui(g, 1) > 0) && (mpz_cmp(g, N) < 0))
        {
            mpz_set(fobj->pp1_obj.gmp_f, g);

            //check if the factor is prime
            if (is_mpz_prp(fobj->pp1_obj.gmp_f, fobj->NUM_WITNESSES))
            {
                add_to_factor_list(fobj->factors, fobj->pp1_obj.gmp_f,
                    fobj->VFLAG, fobj->NUM_WITNESSES);

                if (fobj->VFLAG > 0)
                    gmp_printf("pp1: found prp%d factor = %Zd in p+/-1 stage 1, lane %d\n",
                        gmp_base10(fobj->pp1_obj.gmp_f), fobj->pp1_obj.gmp_f, i);

                char* s = mpz_get_str(NULL, 10, fobj->pp1_obj.gmp_f);
                logprint_oc(fobj->flogname, "a", "prp%d = %s\n",
                    gmp_base10(fobj->pp1_obj.gmp_f), s);
                free(s);
            }
            else
            {
                add_to_factor_list(fobj->factors, fobj->pp1_obj.gmp_f,
                    fobj->VFLAG, fobj->NUM_WITNESSES);

                if (fobj->VFLAG > 0)
                    gmp_printf("pp1: found c%d factor = %Zd in p+/-1 stage 1, lane %d\n",
                        gmp_base10(fobj->pp1_obj.gmp_f), fobj->pp1_obj.gmp_f, i);

                char* s = mpz_get_str(NULL, 10, fobj->pp1_obj.gmp_f);
                logprint_oc(fobj->flogname, "a", "c%d = %s\n",
                    gmp_base10(fobj->pp1_obj.gmp_f), s);
                free(s);
            }

            found = 1;
            mpz_tdiv_q(fobj->pp1_obj.gmp_n, fobj->pp1_obj.gmp_n, fobj->pp1_obj.gmp_f);
        }
    }

    gettimeofday(&stopt, NULL);
    t_time = ytools_difftime(&startt, &stopt);
    if (verbose > 0)
    {
        printf("pp1: Stage1 took %1.4f seconds.\n", t_time);
    }

    if (found)
    {
        
        vecFree(vecV);
        vec_monty_free(montyconst);
        mpz_clear(N);
        mpz_clear(r);
        mpz_clear(g);
        soe_finalize(sdata);
        free(ppm1_primes);
        return;
    }

    // use gmp-ecm internal stage 2.
    // P+/-1 is not like ecm where we can just run more curves
    // to increase the probability of success.  So a vastly inferior
    // stage 2 (standard prime pairing continuation) doesn't make
    // much sense when we have gmp-ecm's stage 2 available.  The
    // greatly increased B2 makes up for not running vectorized.
    int status;
    ecm_params params;
    ecm_init(params);
    params->method = ECM_PP1;

    for (i = 0; i < VECLEN; i++)
    {
        params->B1done = STAGE1_MAX;
        if (fobj->VFLAG >= 3)
            params->verbose = fobj->VFLAG - 2;

        // grab this lane and take out of Montgomery rep.
        extract_bignum_from_vec_to_mpz(params->x, vecV, i, nwords);
        mpz_mul(r, params->x, montyconst->nhat);
        mpz_tdiv_r_2exp(r, r, maxbits);
        mpz_mul(r, r, N);
        mpz_add(r, r, params->x);
        mpz_tdiv_q_2exp(r, r, maxbits);
        if (mpz_cmp(r, N) > 0)
        {
            mpz_sub(r, r, N);
        }
        mpz_set(params->x, r);

        if (fobj->pp1_obj.stg2_is_default == 0)
        {
            //not default, tell gmp-ecm to use the requested B2
            //printf("using requested B2 value\n");
            uint64_2gmp(fobj->pp1_obj.B2, params->B2);
        }

        status = ecm_factor(g, N, STAGE1_MAX, params);

        if ((mpz_cmp_ui(g, 1) > 0) && (mpz_cmp(g, N) < 0))
        {
            //gmp_printf("found nontrivial factor %Zd in p+/-1 stage 2, lane %d\n", g, i);
            mpz_set(fobj->pp1_obj.gmp_f, g);

            //check if the factor is prime
            if (is_mpz_prp(fobj->pp1_obj.gmp_f, fobj->NUM_WITNESSES))
            {
                add_to_factor_list(fobj->factors, fobj->pp1_obj.gmp_f,
                    fobj->VFLAG, fobj->NUM_WITNESSES);

                if (fobj->VFLAG > 0)
                    gmp_printf("pp1: found prp%d factor = %Zd in p+/-1 stage 2, lane %d\n",
                        gmp_base10(fobj->pp1_obj.gmp_f), fobj->pp1_obj.gmp_f, i);

                char* s = mpz_get_str(NULL, 10, fobj->pp1_obj.gmp_f);
                logprint_oc(fobj->flogname, "a", "prp%d = %s\n",
                    gmp_base10(fobj->pp1_obj.gmp_f), s);
                free(s);
            }
            else
            {
                add_to_factor_list(fobj->factors, fobj->pp1_obj.gmp_f,
                    fobj->VFLAG, fobj->NUM_WITNESSES);

                if (fobj->VFLAG > 0)
                    gmp_printf("pp1: found c%d factor = %Zd in p+/-1 stage 2, lane %d\n",
                        gmp_base10(fobj->pp1_obj.gmp_f), fobj->pp1_obj.gmp_f, i);

                char* s = mpz_get_str(NULL, 10, fobj->pp1_obj.gmp_f);
                logprint_oc(fobj->flogname, "a", "c%d = %s\n",
                    gmp_base10(fobj->pp1_obj.gmp_f), s);
                free(s);
            }

            mpz_tdiv_q(fobj->pp1_obj.gmp_n, fobj->pp1_obj.gmp_n, fobj->pp1_obj.gmp_f);
            break;
        }
    }
    ecm_clear(params);

    gettimeofday(&stopt, NULL);
    t_time = ytools_difftime(&startt, &stopt);

    if (verbose > 0)
    {
        printf("pp1: Process took %1.4f seconds.\n", t_time);
    }

    soe_finalize(sdata);
    vec_monty_free(montyconst);
    mpz_clear(N);
    mpz_clear(r);
    mpz_clear(g);
    free(ppm1_primes);
    return;
}





