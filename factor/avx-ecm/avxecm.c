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

#include <stdio.h>
#include "avx_ecm.h"
#include "soe.h"
#include "arith.h"


//#define D 1155
//#define U 8
//#define R 483
//#define L 2 * U

//#define DEBUG 1

// local functions
void vec_add(vec_monty_t *mdata, ecm_work *work, ecm_pt *Pin, ecm_pt *Pout);
void vec_duplicate(vec_monty_t *mdata, ecm_work *work, vec_bignum_t *insum, vec_bignum_t *indiff, ecm_pt *P);
void next_pt_vec(vec_monty_t *mdata, ecm_work *work, ecm_pt *P, uint64_t c);
void euclid(vec_monty_t *mdata, ecm_work *work, ecm_pt *P, uint64_t c);
void vec_prac(vec_monty_t *mdata, ecm_work *work, ecm_pt *P, uint64_t c);
void vec_build_one_curve(thread_data_t *tdata, mpz_t X, mpz_t Z, mpz_t A, uint64_t sigma);
void build_one_curve_param1(thread_data_t *tdata, mpz_t X, mpz_t Z, mpz_t X2, mpz_t Z2, 
    mpz_t A, uint64_t sigma);
void vec_ecm_stage1(vec_monty_t *mdata, ecm_work *work, ecm_pt *P, 
    uint64_t b1, uint64_t*primes, uint64_t nump, int verbose);
int vec_ecm_stage2_init(ecm_pt *P, vec_monty_t *mdata, ecm_work *work, int verbose);
//void vec_ecm_stage2_pair(ecm_pt *P, vec_monty_t *mdata, ecm_work *work, base_t *primes, int verbose);
void vec_ecm_stage2_pair(uint32_t pairmap_steps, uint32_t* pairmap_v, uint32_t* pairmap_u,
    ecm_pt* P, vec_monty_t* mdata, ecm_work* work, int verbose);
uint32_t pair(uint32_t* pairmap_v, uint32_t* pairmap_u,
    ecm_work* work, Queue_t** Q, uint32_t* Qrmap, uint32_t* Qmap,
    uint64_t* primes, uint64_t nump, uint64_t B1, uint64_t B2, int verbose);

int vec_check_factor(mpz_t Z, mpz_t n, mpz_t f);

// GMP-ECM batch stage1
void array_mul(uint64_t* primes, uint64_t b1, int num_p, mpz_t piprimes);
unsigned int compute_s(mpz_t s, uint64_t *primes, uint64_t B1);
int vec_ecm_stage1_batch(mpz_t f, ecm_work *work, ecm_pt *P, vec_bignum_t * A, 
    vec_bignum_t * n, uint64_t B1, mpz_t s, vec_monty_t *mdata);

// support parallel MMIX by requiring the lcg's state to be passed in
uint64_t par_lcg_rand(uint64_t* lcg_state)
{
    // Knuth's MMIX LCG
    *lcg_state = 6364136223846793005ULL * *lcg_state + 1442695040888963407ULL;
    return *lcg_state;
}

#define INV_2_POW_32 2.3283064365386962890625e-10
// Knuth's 64 bit MMIX LCG, using a global 64 bit state variable.
uint32_t spRandp(uint64_t * lcg_state, uint32_t lower, uint32_t upper)
{
    // advance the state of the LCG and return the appropriate result
    *lcg_state = par_lcg_rand(lcg_state);
    return lower + (uint32_t)(
        (double)(upper - lower) * (double)(*lcg_state >> 32) * INV_2_POW_32);
}

int mulcnt[256];
int sqrcnt[256];

#include "threadpool.h"

void vec_ecm_sync(void *vptr);
void vec_ecm_dispatch(void *vptr);
void vec_ecm_build_curve_work_fcn(void *vptr);
void vec_ecm_stage1_work_fcn(void *vptr);
void vec_ecm_stage2_work_fcn(void *vptr);


void
mpres_get_z(mpz_t R, const mpz_t S, mpz_t modulus, vec_monty_t *mdata, uint32_t MAXBITS)
{
    // not fast, but gets the input S out of Mrep.
    mpz_set(mdata->gmp_t1, S);
    //ecm_redc_basecase(R, mdata->gmp_t1, modulus);
    mpz_mul(mdata->gmp_t2, mdata->nhat, S);
    mpz_tdiv_r_2exp(mdata->gmp_t2, mdata->gmp_t2, MAXBITS);
    mpz_mul(mdata->gmp_t2, modulus, mdata->gmp_t2);
    mpz_add(mdata->gmp_t2, S, mdata->gmp_t2);
    mpz_tdiv_q_2exp(mdata->gmp_t2, mdata->gmp_t2, MAXBITS);
    mpz_mod(R, mdata->gmp_t2, modulus);
}

/* R <- S / 2^n mod modulus. Does not need to be fast. */
void mpres_div_2exp(mpz_t R, mpz_t S, const unsigned int n, mpz_t modulus)
{
    int i;

    if (n == 0)
    {
        mpz_set(R, S);
        return;
    }

    if (mpz_odd_p(S))
    {
        mpz_add(R, S, modulus);
        mpz_tdiv_q_2exp(R, R, 1);
    }
    else
        mpz_tdiv_q_2exp(R, S, 1);

    for (i = n; i > 1; i--)
    {
        if (mpz_odd_p(R))
        {
            mpz_add(R, R, modulus);
        }
        mpz_tdiv_q_2exp(R, R, 1);
    }

}

/* R <- S * 2^k mod modulus */
void mpres_mul_2exp(mpz_t R, mpz_t S, const unsigned long k, mpz_t modulus)
{
    mpz_mul_2exp(R, S, k);
    /* This is the same for all methods: just reduce with original modulus */
    mpz_mod(R, R, modulus);
}

void mpres_add_ui(mpz_t R, mpz_t S, const unsigned long n, mpz_t modulus, uint32_t MAXBITS)
{
    mpz_set_ui(R, n);
    mpz_mul_2exp(R, R, MAXBITS);
    mpz_add(R, R, S);
    mpz_mod(R, R, modulus);
}


void vec_ecm_sync(void *vptr)
{
    tpool_t *tpdata = (tpool_t *)vptr;
    thread_data_t *udata = (thread_data_t *)tpdata->user_data;
    uint32_t tid = tpdata->tindex;

    if (udata[tid].ecm_phase == 0)
    {
        udata[tid].phase_done = 1;
    }
    else if (udata->ecm_phase == 1)
    {
        udata[tid].phase_done = 1;
    }
    else if (udata->ecm_phase == 2)
    {
        udata[tid].phase_done = 1;
    }
    else if (udata->ecm_phase == 3)
    {
        udata[tid].phase_done = 1;
    }
    
    return;
}

void vec_ecm_dispatch(void *vptr)
{
    tpool_t *tpdata = (tpool_t *)vptr;
    thread_data_t *udata = (thread_data_t *)tpdata->user_data;
    uint32_t tid = tpdata->tindex;
    
    if (udata[tid].ecm_phase == 0)
    {        
        if (udata[tid].phase_done == 0)
        {
            tpdata->work_fcn_id = 0;
        }
        else
        {
            tpdata->work_fcn_id = tpdata->num_work_fcn;
        }
    }
    else if (udata[tid].ecm_phase == 1)
    {
        if (udata[tid].phase_done == 0)
        {
            tpdata->work_fcn_id = 1;
        }
        else
        {
            tpdata->work_fcn_id = tpdata->num_work_fcn;
        }
    }
    else if (udata[tid].ecm_phase == 2)
    {
        if (udata[tid].phase_done == 0)
        {
            tpdata->work_fcn_id = 2;
        }
        else
        {
            tpdata->work_fcn_id = tpdata->num_work_fcn;
        }
    }
    else if (udata[tid].ecm_phase == 3)
    {
        if (udata[tid].phase_done == 0)
        {
            tpdata->work_fcn_id = 3;
        }
        else
        {
            tpdata->work_fcn_id = tpdata->num_work_fcn;
        }
    }

    return;
}

void vec_ecm_stage1_work_fcn(void *vptr)
{
    tpool_t *tpdata = (tpool_t *)vptr;
    thread_data_t *udata = (thread_data_t *)tpdata->user_data;
    uint32_t tid = tpdata->tindex;
    int verbose = (tid == 0) ? udata[tid].verbose : 0;

    vec_ecm_stage1(udata[tid].mdata, udata[tid].work, udata[tid].P, 
        udata[tid].STAGE1_MAX, udata[tid].primes, udata[tid].nump, 
        verbose);

    return;
}

void vec_ecm_stage2_init_work_fcn(void *vptr)
{
    tpool_t *tpdata = (tpool_t *)vptr;
    thread_data_t *udata = (thread_data_t *)tpdata->user_data;
    uint32_t tid = tpdata->tindex;
    int verbose = (tid == 0) ? udata[tid].verbose : 0;

    vec_ecm_stage2_init(udata[tid].P, udata[tid].mdata, udata[tid].work, verbose);

    return;
}

void vec_ecm_stage2_work_fcn(void *vptr)
{
    tpool_t *tpdata = (tpool_t *)vptr;
    thread_data_t *udata = (thread_data_t *)tpdata->user_data;
    uint32_t tid = tpdata->tindex;
    int verbose = (tid == 0) ? udata[tid].verbose : 0;

    vec_ecm_stage2_pair(udata[tid].pairmap_steps, udata[tid].pairmap_v, udata[tid].pairmap_u,
        udata[tid].P, udata[tid].mdata, udata[tid].work, verbose);
    udata[tid].work->last_pid = udata[tid].nump;

    return;
}

//#define PARAM1

void vec_ecm_build_curve_work_fcn(void *vptr)
{
    tpool_t *tpdata = (tpool_t *)vptr;
    thread_data_t *tdata = (thread_data_t *)tpdata->user_data;
    uint32_t tid = tpdata->tindex;
    mpz_t X, Z, A; // , t, n, one, r;
#ifdef PARAM1
    mpz_t X2, Z2;
#endif
    int i;

    mpz_init(X);
    mpz_init(Z);
    mpz_init(A);
#ifdef PARAM1
    mpz_init(X2);
    mpz_init(Z2);
#endif

    vecClear(tdata[tid].P->X);
    vecClear(tdata[tid].P->Z);
    vecClear(tdata[tid].work->s);

	for (i = 0; i < VECLEN; i++)
    {
        int j;
#ifdef PARAM1
        build_one_curve_param1(&tdata[tid], X, Z, X2, Z2, A, 0); // 65953669); // 36875098);
        insert_mpz_to_vec(tdata[tid].P->X, X, i);
        insert_mpz_to_vec(tdata[tid].P->Z, Z, i);
        insert_mpz_to_vec(tdata[tid].work->pt2.X, X2, i);
        insert_mpz_to_vec(tdata[tid].work->pt2.Z, Z2, i);
        insert_mpz_to_vec(tdata[tid].work->s, A, i);
#else
        // 8689346476060549ULL
        vec_build_one_curve(&tdata[tid], X, Z, A, tdata[tid].sigma[i]); // 6710676431370252287 + i);

        insert_mpz_to_vec(tdata[tid].P->X, X, i);
        insert_mpz_to_vec(tdata[tid].P->Z, Z, i);
        insert_mpz_to_vec(tdata[tid].work->s, A, i);
#endif
        tdata[tid].sigma[i] = tdata[tid].work->sigma;
    }

    //print_vechexbignum(tdata[tid].P->X, "POINT x: ");
    //print_vechexbignum(tdata[tid].P->Z, "POINT z: ");
    //print_vechexbignum(tdata[tid].work->s, "POINT a: ");

    tdata[tid].work->diff1->size = tdata[tid].work->n->size;
    tdata[tid].work->diff2->size = tdata[tid].work->n->size;
    tdata[tid].work->sum1->size  = tdata[tid].work->n->size;
    tdata[tid].work->sum2->size  = tdata[tid].work->n->size;
    tdata[tid].work->pt1.X->size = tdata[tid].work->n->size;
    tdata[tid].work->pt1.Z->size = tdata[tid].work->n->size;
    tdata[tid].work->pt2.X->size = tdata[tid].work->n->size;
    tdata[tid].work->pt2.Z->size = tdata[tid].work->n->size;
    tdata[tid].work->tt1->size   = tdata[tid].work->n->size;
    tdata[tid].work->tt2->size   = tdata[tid].work->n->size;
    tdata[tid].work->tt3->size   = tdata[tid].work->n->size;
    tdata[tid].work->tt4->size   = tdata[tid].work->n->size;

    mpz_clear(X);
    mpz_clear(Z);
    mpz_clear(A);

#ifdef PARAM1
    mpz_clear(X2);
    mpz_clear(Z2);
#endif

    return;
}

void vec_ecm_work_init(ecm_work* work)
{
    int i, j, m;
    uint32_t U = work->U;
    uint32_t L = work->L;
    uint32_t D = work->D;
    uint32_t R = work->R;
    uint32_t words = work->NWORDS;

    work->ptadds = 0;
    work->ptdups = 0;
    work->numinv = 0;
    work->numprimes = 0;
    work->paired = 0;

    work->diff1 = vecInit(words);
    work->diff2 = vecInit(words);
    work->sum1 = vecInit(words);
    work->sum2 = vecInit(words);
    work->tt1 = vecInit(words);
    work->tt2 = vecInit(words);
    work->tt3 = vecInit(words);
    work->tt4 = vecInit(words);
    work->tt5 = vecInit(words);
    work->s = vecInit(words);
    work->n = vecInit(words);

    work->Paprod = (vec_bignum_t **)malloc((2 * L) * sizeof(vec_bignum_t*));
    work->Pa_inv = (vec_bignum_t **)malloc((2 * L) * sizeof(vec_bignum_t*));
    work->Pa = (ecm_pt*)malloc((2 * L) * sizeof(ecm_pt));
    for (j = 0; j < (2 * L); j++)
    {
        work->Paprod[j] = vecInit(words);
        work->Pa_inv[j] = vecInit(words);
        vec_ecm_pt_init(&work->Pa[j], words);
    }

    work->Pad = (ecm_pt*)malloc(sizeof(ecm_pt));
    vec_ecm_pt_init(work->Pad, words);

    work->Pdnorm = (ecm_pt*)malloc(sizeof(ecm_pt));
    vec_ecm_pt_init(work->Pdnorm, words);

    // build an array to hold values of f(b).
    // will need to be of size U*R, allowed values will be
    // up to a multiple of U of the residues mod D
    work->Pb = (ecm_pt*)malloc(U * (R + 1) * sizeof(ecm_pt));
    work->Pbprod = (vec_bignum_t **)malloc(U * (R + 1) * sizeof(vec_bignum_t*));
    for (j = 0; j < U * (R + 1); j++)
    {
        vec_ecm_pt_init(&work->Pb[j], words);
        work->Pbprod[j] = vecInit(words);
    }

    work->map = (uint32_t*)calloc(U * (D + 1) + 3, sizeof(uint32_t));

    work->map[0] = 0;
    work->map[1] = 1;
    work->map[2] = 2;
    m = 3;
    for (i = 0; i < U; i++)
    {
        if (i == 0)
            j = 3;
        else
            j = 1;

        for (; j < D; j++)
        {
            if (spGCD(j, D) == 1)
            {
                work->map[i * D + j] = m++;
            }
            else
            {
                work->map[i * D + j] = 0;
            }
        }

        if (i == 0)
            work->map[i * D + j] = m++;

    }

    vec_ecm_pt_init(&work->pt1, words);
    vec_ecm_pt_init(&work->pt2, words);
    vec_ecm_pt_init(&work->pt3, words);
    vec_ecm_pt_init(&work->pt4, words);
    vec_ecm_pt_init(&work->pt5, words);

    work->stg2acc = vecInit(words);

    return;
}

void vec_ecm_work_free(ecm_work* work)
{
    int i, j;
    uint32_t U = work->U;
    uint32_t L = work->L;
    uint32_t D = work->D;
    uint32_t R = work->R;

    vecFree(work->diff1);
    vecFree(work->diff2);
    vecFree(work->sum1);
    vecFree(work->sum2);
    vecFree(work->tt1);
    vecFree(work->tt2);
    vecFree(work->tt3);
    vecFree(work->tt4);
    vecFree(work->tt5);
    vecFree(work->s);
    vecFree(work->n);
    vec_ecm_pt_free(&work->pt1);
    vec_ecm_pt_free(&work->pt2);
    vec_ecm_pt_free(&work->pt3);
    vec_ecm_pt_free(&work->pt4);
    vec_ecm_pt_free(&work->pt5);

    for (i = 0; i < (2 * L); i++)
    {
        vec_ecm_pt_free(&work->Pa[i]);
        vecFree(work->Paprod[i]);
        vecFree(work->Pa_inv[i]);
    }

    vec_ecm_pt_free(work->Pad);
    vec_ecm_pt_free(work->Pdnorm);
    free(work->Pa);
    free(work->Pad);
    free(work->Pdnorm);
    free(work->map);

    for (i = 0; i < U * (R + 1); i++)
    {
        vec_ecm_pt_free(&work->Pb[i]);
        vecFree(work->Pbprod[i]);
    }

    free(work->Pbprod);
    free(work->Paprod);
    free(work->Pb);
    vecFree(work->stg2acc);

    return;
}

void vec_ecm_pt_init(ecm_pt *pt, uint32_t words)
{
	pt->X = vecInit(words);
	pt->Z = vecInit(words);
}

void vec_ecm_pt_free(ecm_pt *pt)
{
	vecFree(pt->X);
	vecFree(pt->Z);
}

void mpz_mulmod(mpz_t out, mpz_t in1, mpz_t in2, mpz_t n, mpz_t nhat, mpz_t s, int bits)
{
    mpz_mul(out, in1, in2);

    mpz_mul(s, out, nhat);
    mpz_tdiv_r_2exp(s, s, bits);

    mpz_mul(s, s, n);
    mpz_add(out, s, out);
    mpz_tdiv_q_2exp(out, out, bits);
    mpz_tdiv_q_2exp(s, out, bits);

    if (mpz_cmp_ui(s, 0) > 0)
    {
        mpz_sub(out, out, n);
    }
    return;
}

void mpz_sqrmod(mpz_t out, mpz_t in, mpz_t n, mpz_t nhat, mpz_t s, int bits)
{
    mpz_mulmod(out, in, in, n, nhat, s, bits);
    return;
}

void mpz_submod(mpz_t out, mpz_t a, mpz_t b, mpz_t n)
{
    mpz_sub(out, a, b);

    if (mpz_cmp_ui(out, 0) < 0)
    {
        mpz_add(out, out, n);
    }

    if (mpz_cmp_ui(out, 0) < 0)
    {
        mpz_add(out, out, n);
    }

    return;
}

void mpz_addmod(mpz_t out, mpz_t a, mpz_t b, mpz_t n, mpz_t s, int bits)
{
    mpz_add(out, a, b);
    mpz_tdiv_q_2exp(s, out, bits);
    if (mpz_cmp_ui(s, 0) > 0)
    {
        mpz_sub(out, out, n);
    }
    return;
}

void vec_add(vec_monty_t *mdata, ecm_work *work, ecm_pt *Pin, ecm_pt *Pout)
{
    // compute:
    //x+ = z- * [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
    //z+ = x- * [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
    // where:
    //x- = original x
    //z- = original z
    // given the sums and differences of the original points (stored in work structure).
    
    vecmulmod_ptr(work->diff1, work->sum2, work->tt1, work->n, work->tt4, mdata);	//U
    vecmulmod_ptr(work->sum1, work->diff2, work->tt2, work->n, work->tt4, mdata);	//V

    vecsubmod_ptr(work->tt1, work->tt2, work->tt4, mdata);
    vecaddmod_ptr(work->tt1, work->tt2, work->tt3, mdata);

    vecsqrmod_ptr(work->tt3, work->tt1, work->n, work->tt5, mdata);					//(U + V)^2
    vecsqrmod_ptr(work->tt4, work->tt2, work->n, work->tt5, mdata);					//(U - V)^2

    // choosing the initial point Pz0 = 1 means that z_p-q = 1 and this mul isn't necessary...
    // but that involves a different way to initialize curves, so for now
    // we can't assume Z=1.
    // with streamlined prac, the if case is never executed...
	if (Pin->X == Pout->X)
	{
		base_t *swap;
		vecmulmod_ptr(work->tt1, Pin->Z, Pout->Z, work->n, work->tt4, mdata);		//Z * (U + V)^2
		vecmulmod_ptr(work->tt2, Pin->X, Pout->X, work->n, work->tt4, mdata);		//x * (U - V)^2
		swap = Pout->Z->data;
		Pout->Z->data = Pout->X->data;
		Pout->X->data = swap;
	}
	else
	{
		vecmulmod_ptr(work->tt1, Pin->Z, Pout->X, work->n, work->tt4, mdata);		//Z * (U + V)^2
		vecmulmod_ptr(work->tt2, Pin->X, Pout->Z, work->n, work->tt4, mdata);		//x * (U - V)^2
	}
	work->stg1Add++;
    return;
}

void vec_duplicate(vec_monty_t *mdata, ecm_work *work, vec_bignum_t *insum, vec_bignum_t *indiff, ecm_pt *P)
{
    vecsqrmod_ptr(indiff, work->tt1, work->n, work->tt4, mdata);			    // V=(x1 - z1)^2
    vecsqrmod_ptr(insum, work->tt2, work->n, work->tt4, mdata);			        // U=(x1 + z1)^2
    vecmulmod_ptr(work->tt1, work->tt2, P->X, work->n, work->tt4, mdata);       // x=U*V

    vecsubmod_ptr(work->tt2, work->tt1, work->tt3, mdata);	                // w = U-V
    vecmulmod_ptr(work->tt3, work->s, work->tt2, work->n, work->tt4, mdata);    // t = (A+2)/4 * w
    vecaddmod_ptr(work->tt2, work->tt1, work->tt2, mdata);                    // t = t + V
    vecmulmod_ptr(work->tt2, work->tt3, P->Z, work->n, work->tt4, mdata);       // Z = t*w
	work->stg1Doub++;
    return;
}

void mpz_elladd(mpz_t Px_out, mpz_t Pz_out, 
    mpz_t Px1, mpz_t Pz1, mpz_t Px2, mpz_t Pz2, mpz_t Px3, mpz_t Pz3, 
    mpz_t n, mpz_t nhat, mpz_t sum, mpz_t diff,
    mpz_t tmp, mpz_t tmp2, int bits)
{
    mpz_submod(diff, Px1, Pz1, n);
    mpz_addmod(sum, Px2, Pz2, n, tmp, bits);
    mpz_mulmod(tmp2, diff, sum, n, nhat, tmp, bits);    // U

    mpz_submod(diff, Px2, Pz2, n);
    mpz_addmod(sum, Px1, Pz1, n, tmp, bits);
    mpz_mulmod(sum, diff, sum, n, nhat, tmp, bits);     // V

    mpz_submod(diff, tmp2, sum, n);                     // U - V
    mpz_addmod(sum, tmp2, sum, n, tmp, bits);                // U + V

    mpz_sqrmod(diff, diff, n, nhat, tmp, bits);
    mpz_sqrmod(sum, sum, n, nhat, tmp, bits);

    mpz_mulmod(Px_out, Pz3, sum, n, nhat, tmp, bits);
    mpz_mulmod(Pz_out, Px3, diff, n, nhat, tmp, bits);

    return;
}

void mpz_elldup(mpz_t Px, mpz_t Pz, mpz_t n, mpz_t nhat, mpz_t sum, mpz_t diff, 
    mpz_t s, mpz_t tmp, mpz_t tmp2, int bits)
{
    mpz_submod(diff, Px, Pz, n);
    mpz_sqrmod(diff, diff, n, nhat, tmp, bits);     // V=(x1 - z1)^2
    mpz_addmod(sum, Px, Pz, n, tmp, bits);
    mpz_sqrmod(sum, sum, n, nhat, tmp, bits);       // U=(x1 + z1)^2
    mpz_mulmod(Px, sum, diff, n, nhat, tmp, bits);  // x=U*V
    mpz_submod(tmp2, sum, diff, n);                 // w = U-V
    mpz_mulmod(Pz, tmp2, s, n, nhat, tmp, bits);    // t = (A+2)/4 * w
    mpz_addmod(Pz, Pz, diff, n, tmp, bits);         // t = t + V
    mpz_mulmod(Pz, Pz, tmp2, n, nhat, tmp, bits);   // Z = t*w
    return;
}

double vec_getEcost(uint64_t d, uint64_t e)
{
	int doub = 0, add = 0;

	while (d > 0)
	{
		if ((e / 2) < d)
		{
			d = e - d;
		}
		else if ((d < (e / 4)) && ((e & 1) == 0))
		{
			e = e / 2;
			doub++;
			add++;
		}
		else
		{
			e = e - d;
			add++;
		}

	}
	return (doub + add) * 2 * 0.75 + add * 4 + doub * 3;
}

int * getEseq(uint64_t d, uint64_t e)
{
	uint64_t target = e;
	int doub = 0, add = 0;
	int seqLen = 64;
	int *seq = (int *)malloc(seqLen * sizeof(int));
	int it = 1;
	seq[0] = 0;

	d = e - d;

	while (d > 0)
	{
		if ((e / 2) < d)
		{
#ifdef DEBUG
			printf("[3]: d = %u\n", e - d);
#endif
			seq[it] = 3;
			d = e - d;
		}
		else if ((d < (e / 4)) && ((e & 1) == 0))
		{
#ifdef DEBUG
			printf("[1]: d = %u, e = %u\n", d, e);
			printf("[1]: d = %u, e = %u\n", d, e - d);
#endif
			e = e / 2;
			doub++;
			add++;
			seq[it] = 1;
		}
		else
		{
#ifdef DEBUG
			printf("[2]: d = %u, e = %u\n", d, e);
#endif
			e = e - d;
			add++;
			seq[it] = 2;
		}

		it++;
		if (it > seqLen)
		{
			seqLen *= 2;
			seq = (int *)realloc(seq, seqLen * sizeof(int));
		}
	}

	// reverse the sequence so it's constructive starting from 0,1
	seqLen = it;
	for (it = 0; it < seqLen / 2; it++)
	{
		uint32_t tmp = seq[it];
		seq[it] = seq[seqLen - it - 1];
		seq[seqLen - it - 1] = tmp;
	}

	// for fun.
#ifdef DEBUG
	if (1)
#else
	if (0)
#endif
	{
		uint32_t d = 0;
		uint32_t x = 1;
		it = 0;

		printf("target is %u, sequence length is %u\n", target, seqLen);
		while (seq[it] != 0)
		{
			if (seq[it] == 1)
			{
				printf("[1]: double to get %u and add %d to get %u\n", x * 2, x - d, x + x - d);
				x *= 2;
			}
			else if (seq[it] == 2)
			{
				x += d;
				printf("[2]: add %d to get %u\n", d, x);
			}
			else if (seq[it] == 3)
			{
				d = x - d;
				printf("[3]: d is now %u\n", d);
			}
			it++;
		}
	}

	return seq;
}

#define ADD 5.5     // accounts for sqr ~= 0.75 mul
#define DUP 4.5     // accounts for sqr ~= 0.75 mul

static double
lucas_cost(uint64_t n, double v)
{
	uint64_t d, e, r;
	double c; /* cost */

	d = n;
	r = (uint64_t)((double)d * v + 0.5);
	if (r >= n)
		return (ADD * (double)n);
	d = n - r;
	e = 2 * r - n;
	c = DUP + ADD; /* initial duplicate and final addition */
	while (d != e)
	{
		if (d < e)
		{
			r = d;
			d = e;
			e = r;
		}
#ifdef ORIG_PRAC
		if (d - e <= e / 4 && ((d + e) % 3) == 0)
		{ /* condition 1 */
			d = (2 * d - e) / 3;
			e = (e - d) / 2;
			c += 3.0 * ADD; /* 3 additions */
		}
		else if (d - e <= e / 4 && (d - e) % 6 == 0)
		{ /* condition 2 */
			d = (d - e) / 2;
			c += ADD + DUP; /* one addition, one duplicate */
		}
		else if ((d + 3) / 4 <= e)
#else
        if ((d + 3) / 4 <= e)
#endif
		{ /* condition 3 */
			d -= e;
			c += ADD; /* one addition */
		}
		else if ((d + e) % 2 == 0)
		{ /* condition 4 */
			d = (d - e) / 2;
			c += ADD + DUP; /* one addition, one duplicate */
		}
		/* now d+e is odd */
		else if (d % 2 == 0)
		{ /* condition 5 */
			d /= 2;
			c += ADD + DUP; /* one addition, one duplicate */
		}
		/* now d is odd and e is even */
#ifdef ORIG_PRAC
		else if (d % 3 == 0)
		{ /* condition 6 */
			d = d / 3 - e;
			c += 3.0 * ADD + DUP; /* three additions, one duplicate */
		}
		else if ((d + e) % 3 == 0)
		{ /* condition 7 */
			d = (d - 2 * e) / 3;
			c += 3.0 * ADD + DUP; /* three additions, one duplicate */
		}
		else if ((d - e) % 3 == 0)
		{ /* condition 8 */
			d = (d - e) / 3;
			c += 3.0 * ADD + DUP; /* three additions, one duplicate */
		}
#endif
		else /* necessarily e is even: catches all cases */
		{ /* condition 9 */
			e /= 2;
			c += ADD + DUP; /* one addition, one duplicate */
		}
	}

	return c;
}

//#define MPZ_VERIFY

void vec_prac(vec_monty_t *mdata, ecm_work *work, ecm_pt *P, uint64_t c)
{
	uint64_t d, e, r;
	double cmin, cost;
	int i;
	vec_bignum_t *s1, *s2, *d1, *d2;
	base_t *sw_x, *sw_z;

#ifdef MPZ_VERIFY
    mpz_t px1, px2, px3, px4, gtmp, gtmp2, gn, gnhat;
    mpz_t pz1, pz2, pz3, pz4, gs, gsum, gdiff;
#endif

#define NV 10  
	/* 1/val[0] = the golden ratio (1+sqrt(5))/2, and 1/val[i] for i>0
	   is the real number whose continued fraction expansion is all 1s
	   except for a 2 in i+1-st place */
	static double val[NV] =
	{ 0.61803398874989485, 0.72360679774997897, 0.58017872829546410,
	  0.63283980608870629, 0.61242994950949500, 0.62018198080741576,
	  0.61721461653440386, 0.61834711965622806, 0.61791440652881789,
	  0.61807966846989581 };

	/* chooses the best value of v */
    cmin = 999999999.0; // ADD* (double)c;
	for (i = d = 0; d < NV; d++)
	{
		cost = lucas_cost(c, val[d]);
		if (cost < cmin)
		{
			cmin = cost;
			i = d;
		}
	}

#ifdef MPZ_VERIFY
    mpz_init(px1);
    mpz_init(px2);
    mpz_init(px3);
    mpz_init(px4);
    mpz_init(pz1);
    mpz_init(pz2);
    mpz_init(pz3);
    mpz_init(pz4);
    mpz_init(gn);
    mpz_init(gnhat);
    mpz_init(gtmp);
    mpz_init(gtmp2);
    mpz_init(gs);
    mpz_init(gsum);
    mpz_init(gdiff);
#endif

	d = c;
	r = (uint64_t)((double)d * val[i] + 0.5);

	s1 = work->sum1;
	s2 = work->sum2;
	d1 = work->diff1;
	d2 = work->diff2;

	/* first iteration always begins by Condition 3, then a swap */
	d = c - r;
	e = 2 * r - c;

    //printf("c = %" PRIu64 ", d = %" PRIu64 " e = %" PRIu64 " r = %" PRIu64 "\n",
    //    c, d, e, r);

	// mpres_set(xB, xA, n);
	// mpres_set(zB, zA, n); /* B=A */
	// mpres_set(xC, xA, n);
	// mpres_set(zC, zA, n); /* C=A */
	// duplicate(xA, zA, xA, zA, n, b, u, v, w); /* A = 2*A */

	// the first one is always a doubling
	// point1 is [1]P
	vecCopy(P->X, work->pt1.X);
	vecCopy(P->Z, work->pt1.Z);
	vecCopy(P->X, work->pt2.X);
	vecCopy(P->Z, work->pt2.Z);
	vecCopy(P->X, work->pt3.X);
	vecCopy(P->Z, work->pt3.Z);
	vecsubmod_ptr(work->pt1.X, work->pt1.Z, d1, mdata);
	vecaddmod_ptr(work->pt1.X, work->pt1.Z, s1, mdata);

#ifdef MPZ_VERIFY

    mpz_set(px1, mdata->gmp_t1);
    mpz_set(pz1, mdata->gmp_t2);
    mpz_set(px2, px1);
    mpz_set(px3, px1);
    mpz_set(pz2, pz1);
    mpz_set(pz3, pz1);
    extract_bignum_from_vec_to_mpz(gs, work->s, 0, mdata->NWORDS);
    extract_bignum_from_vec_to_mpz(gn, work->n, 0, mdata->NWORDS);
    extract_bignum_from_vec_to_mpz(gnhat, mdata->vnhat, 0, mdata->NWORDS);

    int bits = mdata->NWORDS * 52;
    mpz_elldup(px1, pz1, gn, gnhat, gsum, gdiff, gs, gtmp, gtmp2, bits);
#endif

	// point2 is [2]P
	vec_duplicate(mdata, work, s1, d1, &work->pt1);

	while (d != e)
	{
		if (d < e)
		{
			r = d;
			d = e;
			e = r;
			//mpres_swap(xA, xB, n);
			//mpres_swap(zA, zB, n);
			sw_x = work->pt1.X->data;
			sw_z = work->pt1.Z->data;
			work->pt1.X->data = work->pt2.X->data;
			work->pt1.Z->data = work->pt2.Z->data;
			work->pt2.X->data = sw_x;
			work->pt2.Z->data = sw_z;

#ifdef MPZ_VERIFY
            mpz_set(gtmp, px1);
            mpz_set(gtmp2, pz1);
            mpz_set(px1, px2);
            mpz_set(pz1, pz2);
            mpz_set(px2, gtmp);
            mpz_set(pz2, gtmp2);
#endif
		}
		/* do the first line of Table 4 whose condition qualifies */
#ifdef ORIG_PRAC
        if (d - e <= e / 4 && ((d + e) % 3) == 0)
        { /* condition 1 */
            d = (2 * d - e) / 3;
            e = (e - d) / 2;

            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s1, d1, mdata);
            vecaddsubmod_ptr(work->pt2.X, work->pt2.Z, s2, d2, mdata);

            vec_add(mdata, work, &work->pt3, &work->pt4); // T = A + B (C)

            vecaddsubmod_ptr(work->pt4.X, work->pt4.Z, s1, d1, mdata);
            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s2, d2, mdata);

            vec_add(mdata, work, &work->pt2, &work->pt5); // T2 = T + A (B)

            vecaddsubmod_ptr(work->pt2.X, work->pt2.Z, s1, d1, mdata);
            vecaddsubmod_ptr(work->pt4.X, work->pt4.Z, s2, d2, mdata);

            vec_add(mdata, work, &work->pt1, &work->pt2); // B = B + T (A)

            //add3(xT, zT, xA, zA, xB, zB, xC, zC, n, u, v, w); /* T = f(A,B,C) */
            //add3(xT2, zT2, xT, zT, xA, zA, xB, zB, n, u, v, w); /* T2 = f(T,A,B) */
            //add3(xB, zB, xB, zB, xT, zT, xA, zA, n, u, v, w); /* B = f(B,T,A) */
            //mpres_swap(xA, xT2, n);
            //mpres_swap(zA, zT2, n); /* swap A and T2 */

            sw_x = work->pt1.X->data;
            sw_z = work->pt1.Z->data;
            work->pt1.X->data = work->pt5.X->data;
            work->pt1.Z->data = work->pt5.Z->data;
            work->pt5.X->data = sw_x;
            work->pt5.Z->data = sw_z;

        }
        else if (d - e <= e / 4 && (d - e) % 6 == 0)
        { /* condition 2 */
            d = (d - e) / 2;

            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s1, d1, mdata);
            vecaddsubmod_ptr(work->pt2.X, work->pt2.Z, s2, d2, mdata);

            vec_add(mdata, work, &work->pt3, &work->pt2);		// B = A + B (C)
            vec_duplicate(mdata, work, s1, d1, &work->pt1);		// A = 2A

            //add3(xB, zB, xA, zA, xB, zB, xC, zC, n, u, v, w); /* B = f(A,B,C) */
            //duplicate(xA, zA, xA, zA, n, b, u, v, w); /* A = 2*A */

        }
        else if ((d + 3) / 4 <= e)
#else
        if ((d + 3) / 4 <= e)
#endif
		{ /* condition 3 */
			d -= e;
			
            //vecaddsubmod_ptr(work->pt2.X, work->pt2.Z, s1, d1, mdata);
            //vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s2, d2, mdata);
            vecsubmod_ptr(work->pt1.X, work->pt1.Z, d2, mdata);
            vecaddmod_ptr(work->pt1.X, work->pt1.Z, s2, mdata);
            vecsubmod_ptr(work->pt2.X, work->pt2.Z, d1, mdata);
            vecaddmod_ptr(work->pt2.X, work->pt2.Z, s1, mdata);

			vec_add(mdata, work, &work->pt3, &work->pt4);		// T = B + A (C)
			//add3(xT, zT, xB, zB, xA, zA, xC, zC, n, u, v, w); /* T = f(B,A,C) */
			
			/* circular permutation (B,T,C) */
			//tmp = xB;
			//xB = xT;
			//xT = xC;
			//xC = tmp;
			//tmp = zB;
			//zB = zT;
			//zT = zC;
			//zC = tmp;

			sw_x = work->pt2.X->data;
			sw_z = work->pt2.Z->data;
			work->pt2.X->data = work->pt4.X->data;
			work->pt2.Z->data = work->pt4.Z->data;
			work->pt4.X->data = work->pt3.X->data;
			work->pt4.Z->data = work->pt3.Z->data;
			work->pt3.X->data = sw_x;
			work->pt3.Z->data = sw_z;

#ifdef MPZ_VERIFY
            mpz_elladd(px4, pz4, px2, pz2, px1, pz1, px3, pz3,
                gn, gnhat, gsum, gdiff, gtmp, gtmp2, bits);

            mpz_set(gtmp, px2);
            mpz_set(gtmp2, pz2);
            mpz_set(px2, px4);
            mpz_set(pz2, pz4);
            mpz_set(px4, px3);
            mpz_set(pz4, pz3);
            mpz_set(px3, gtmp);
            mpz_set(pz3, gtmp2);
#endif

		}
		else if ((d + e) % 2 == 0)
		{ /* condition 4 */
			d = (d - e) / 2;

            //vecaddsubmod_ptr(work->pt2.X, work->pt2.Z, s1, d1, mdata);
            //vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s2, d2, mdata);
            vecsubmod_ptr(work->pt1.X, work->pt1.Z, d2, mdata);
            vecaddmod_ptr(work->pt1.X, work->pt1.Z, s2, mdata);
            vecsubmod_ptr(work->pt2.X, work->pt2.Z, d1, mdata);
            vecaddmod_ptr(work->pt2.X, work->pt2.Z, s1, mdata);

			vec_add(mdata, work, &work->pt3, &work->pt2);		// B = B + A (C)
			vec_duplicate(mdata, work, s2, d2, &work->pt1);		// A = 2A
			
#ifdef MPZ_VERIFY
            mpz_elladd(px2, pz2, px2, pz2, px1, pz1, px3, pz3,
                gn, gnhat, gsum, gdiff, gtmp, gtmp2, bits);
            mpz_elldup(px1, pz1, gn, gnhat, gsum, gdiff, gs, gtmp, gtmp2, bits);
#endif
			//add3(xB, zB, xB, zB, xA, zA, xC, zC, n, u, v, w); /* B = f(B,A,C) */
			//duplicate(xA, zA, xA, zA, n, b, u, v, w); /* A = 2*A */
		}
		/* now d+e is odd */
		else if (d % 2 == 0)
		{ /* condition 5 */
			d /= 2;
			
            //vecaddsubmod_ptr(work->pt3.X, work->pt3.Z, s1, d1, mdata);
            //vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s2, d2, mdata);
            vecsubmod_ptr(work->pt3.X, work->pt3.Z, d1, mdata);
            vecaddmod_ptr(work->pt3.X, work->pt3.Z, s1, mdata);
            vecsubmod_ptr(work->pt1.X, work->pt1.Z, d2, mdata);
            vecaddmod_ptr(work->pt1.X, work->pt1.Z, s2, mdata);

			vec_add(mdata, work, &work->pt2, &work->pt3);		// C = C + A (B)
			vec_duplicate(mdata, work, s2, d2, &work->pt1);		// A = 2A

#ifdef MPZ_VERIFY
            mpz_elladd(px3, pz3, px3, pz3, px1, pz1, px2, pz2,
                gn, gnhat, gsum, gdiff, gtmp, gtmp2, bits);
            mpz_elldup(px1, pz1, gn, gnhat, gsum, gdiff, gs, gtmp, gtmp2, bits);
#endif
			//add3(xC, zC, xC, zC, xA, zA, xB, zB, n, u, v, w); /* C = f(C,A,B) */
			//duplicate(xA, zA, xA, zA, n, b, u, v, w); /* A = 2*A */
		}
		/* now d is odd, e is even */
#ifdef ORIG_PRAC
		else if (d % 3 == 0)
		{ /* condition 6 */
			d = d / 3 - e;

            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s1, d1, mdata);

			vec_duplicate(mdata, work, s1, d1, &work->pt4);		// T = 2A

            vecaddsubmod_ptr(work->pt2.X, work->pt2.Z, s2, d2, mdata);

			vec_add(mdata, work, &work->pt3, &work->pt5);		// T2 = A + B (C)
			
            vecaddsubmod_ptr(work->pt4.X, work->pt4.Z, s1, d1, mdata);
            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s2, d2, mdata);

			vec_add(mdata, work, &work->pt1, &work->pt1);		// A = T + A (A)

            vecaddsubmod_ptr(work->pt5.X, work->pt5.Z, s2, d2, mdata);

			vec_add(mdata, work, &work->pt3, &work->pt4);		// T = T + T2 (C)

			//duplicate(xT, zT, xA, zA, n, b, u, v, w); /* T = 2*A */
			//add3(xT2, zT2, xA, zA, xB, zB, xC, zC, n, u, v, w); /* T2 = f(A,B,C) */
			//add3(xA, zA, xT, zT, xA, zA, xA, zA, n, u, v, w); /* A = f(T,A,A) */
			//add3(xT, zT, xT, zT, xT2, zT2, xC, zC, n, u, v, w); /* T = f(T,T2,C) */

			/* circular permutation (C,B,T) */
			//tmp = xC;
			//xC = xB;
			//xB = xT;
			//xT = tmp;
			//tmp = zC;
			//zC = zB;
			//zB = zT;
			//zT = tmp;

			sw_x = work->pt3.X->data;
			sw_z = work->pt3.Z->data;
			work->pt3.X->data = work->pt2.X->data;
			work->pt3.Z->data = work->pt2.Z->data;
			work->pt2.X->data = work->pt4.X->data;
			work->pt2.Z->data = work->pt4.Z->data;
			work->pt4.X->data = sw_x;
			work->pt4.Z->data = sw_z;

		}
		else if ((d + e) % 3 == 0)
		{ /* condition 7 */
			d = (d - 2 * e) / 3;

            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s1, d1, mdata);
            vecaddsubmod_ptr(work->pt2.X, work->pt2.Z, s2, d2, mdata);

			vec_add(mdata, work, &work->pt3, &work->pt4);		// T = A + B (C)

            vecaddsubmod_ptr(work->pt4.X, work->pt4.Z, s1, d1, mdata);
            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s2, d2, mdata);

			vec_add(mdata, work, &work->pt2, &work->pt2);		// B = T + A (B)

			vec_duplicate(mdata, work, s2, d2, &work->pt4);		// T = 2A

            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s1, d1, mdata);
            vecaddsubmod_ptr(work->pt4.X, work->pt4.Z, s2, d2, mdata);

			vec_add(mdata, work, &work->pt1, &work->pt1);		// A = A + T (A) = 3A

			//add3(xT, zT, xA, zA, xB, zB, xC, zC, n, u, v, w); /* T = f(A,B,C) */
			//add3(xB, zB, xT, zT, xA, zA, xB, zB, n, u, v, w); /* B = f(T,A,B) */
			//duplicate(xT, zT, xA, zA, n, b, u, v, w);
			//add3(xA, zA, xA, zA, xT, zT, xA, zA, n, u, v, w); /* A = 3*A */
		}
		else if ((d - e) % 3 == 0)
		{ /* condition 8 */
			d = (d - e) / 3;
            // in a run with B1=100M, this condition only executed a few times.

            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s1, d1, mdata);
            vecaddsubmod_ptr(work->pt2.X, work->pt2.Z, s2, d2, mdata);

			vec_add(mdata, work, &work->pt3, &work->pt4);		// T = A + B (C)

            vecaddsubmod_ptr(work->pt3.X, work->pt3.Z, s1, d1, mdata);
            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s2, d2, mdata);

			vec_add(mdata, work, &work->pt2, &work->pt3);		// C = C + A (B)

			//add3(xT, zT, xA, zA, xB, zB, xC, zC, n, u, v, w); /* T = f(A,B,C) */
			//add3(xC, zC, xC, zC, xA, zA, xB, zB, n, u, v, w); /* C = f(A,C,B) */
			//mpres_swap(xB, xT, n);
			//mpres_swap(zB, zT, n); /* swap B and T */
			sw_x = work->pt2.X->data;
			sw_z = work->pt2.Z->data;
			work->pt2.X->data = work->pt4.X->data;
			work->pt2.Z->data = work->pt4.Z->data;
			work->pt4.X->data = sw_x;
			work->pt4.Z->data = sw_z;

            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s2, d2, mdata);

			vec_duplicate(mdata, work, s2, d2, &work->pt4);		// T = 2A

            vecaddsubmod_ptr(work->pt1.X, work->pt1.Z, s1, d1, mdata);
            vecaddsubmod_ptr(work->pt4.X, work->pt4.Z, s2, d2, mdata);

			vec_add(mdata, work, &work->pt1, &work->pt1);		// A = A + T (A) = 3A

			//duplicate(xT, zT, xA, zA, n, b, u, v, w);
			//add3(xA, zA, xA, zA, xT, zT, xA, zA, n, u, v, w); /* A = 3*A */
		}
#endif
		else /* necessarily e is even here */
		{ /* condition 9 */
			e /= 2;
            // in a run with B1=100M, this condition never executed.

            //vecaddsubmod_ptr(work->pt3.X, work->pt3.Z, s1, d1, mdata);
            //vecaddsubmod_ptr(work->pt2.X, work->pt2.Z, s2, d2, mdata);
            vecsubmod_ptr(work->pt2.X, work->pt2.Z, d2, mdata);
            vecaddmod_ptr(work->pt2.X, work->pt2.Z, s2, mdata);
            vecsubmod_ptr(work->pt3.X, work->pt3.Z, d1, mdata);
            vecaddmod_ptr(work->pt3.X, work->pt3.Z, s1, mdata);

			vec_add(mdata, work, &work->pt1, &work->pt3);		// C = C + B (A)
			vec_duplicate(mdata, work, s2, d2, &work->pt2);		// B = 2B

#ifdef MPZ_VERIFY
            mpz_elladd(px3, pz3, px3, pz3, px2, pz2, px1, pz1,
                gn, gnhat, gsum, gdiff, gtmp, gtmp2, bits);
            mpz_elldup(px2, pz2, gn, gnhat, gsum, gdiff, gs, gtmp, gtmp2, bits);
#endif
			//add3(xC, zC, xC, zC, xB, zB, xA, zA, n, u, v, w); /* C = f(C,B,A) */
			//duplicate(xB, zB, xB, zB, n, b, u, v, w); /* B = 2*B */
		}
	}

	vecsubmod_ptr(work->pt1.X, work->pt1.Z, d1, mdata);
	vecaddmod_ptr(work->pt1.X, work->pt1.Z, s1, mdata);
	vecsubmod_ptr(work->pt2.X, work->pt2.Z, d2, mdata);
	vecaddmod_ptr(work->pt2.X, work->pt2.Z, s2, mdata);

	vec_add(mdata, work, &work->pt3, P);		// A = A + B (C)

#ifdef MPZ_VERIFY
    mpz_elladd(px1, pz1, px1, pz1, px2, pz2, px3, pz3,
        gn, gnhat, gsum, gdiff, gtmp, gtmp2, bits);
#endif

	//add3(xA, zA, xA, zA, xB, zB, xC, zC, n, u, v, w);

	if (d != 1)
	{
		printf("problem: d != 1\n");
        exit(1);
	}

#ifdef MPZ_VERIFY

    mpz_set(mdata->gmp_t1, px1);
    mpz_set(mdata->gmp_t2, pz1);

    extract_bignum_from_vec_to_mpz(px1, P->X, 0, mdata->NWORDS);
    extract_bignum_from_vec_to_mpz(pz1, P->Z, 0, mdata->NWORDS);
    
    int result = 0;

    if (mpz_cmp(px1, mdata->gmp_t1) != 0)
    {
        mpz_set_ui(gtmp, 1);
        mpz_mul_2exp(gtmp, gtmp, bits);
        gmp_printf("X coords don't match!\n%Zx\n%Zx\n%Zx\n", px1, mdata->gmp_t1, gtmp);
        result = 1;
    }
    
    if (mpz_cmp(pz1, mdata->gmp_t2) != 0)
    {
        mpz_set_ui(gtmp, 1);
        mpz_mul_2exp(gtmp, gtmp, bits);
        gmp_printf("Z coords don't match!\n%Zx\n%Zx\n%Zx\n", pz1, mdata->gmp_t2, gtmp);
        result = 1;
    }

    mpz_clear(px1);
    mpz_clear(px2);
    mpz_clear(px3);
    mpz_clear(px4);
    mpz_clear(pz1);
    mpz_clear(pz2);
    mpz_clear(pz3);
    mpz_clear(pz4);
    mpz_clear(gn);
    mpz_clear(gnhat);
    mpz_clear(gtmp);
    mpz_clear(gtmp2);
    mpz_clear(gs);
    mpz_clear(gsum);
    mpz_clear(gdiff);

    if (result)
        exit(0);

#endif

	return;

}

void euclid(vec_monty_t *mdata, ecm_work *work, ecm_pt *P, uint64_t c)
{
	uint64_t startd;
	uint64_t numd = 10;
	uint64_t bestd = 0;
	int *seq;
	double best;
	double cost;
	int i;
	uint64_t d, x;
	vec_bignum_t *x1, *z1, *x2, *z2, *x3, *z3, *x4, *z4, *s1, *s2, *d1, *d2;
	base_t *sw_x, *sw_z;

	if (c == 1)
	{
		return;
	}
	
	// thank you gmp-ecm
	static double val[10] =
	{ 0.61803398874989485, 0.72360679774997897, 0.58017872829546410,
	  0.63283980608870629, 0.61242994950949500, 0.62018198080741576,
	  0.61721461653440386, 0.61834711965622806, 0.61791440652881789,
	  0.61807966846989581 };

	i = 0;
	best = 999999.0;

	//while (numd > 0)
	while (i < numd)
	{
		//uint64_t d = startd - i;
		uint64_t d = (uint64_t)((double)c * val[i]);
		uint64_t e = c;

		if (spGCD(d, e) != 1)
		{
			i++;
			continue;
		}

		cost = vec_getEcost(d, e);
		if (cost < best)
		{
			best = cost;
			bestd = d;
		}
		i++;
		numd--;
	}

	if ((c & 1) == 0)
	{
		printf("input should not be even\n");
	}

	seq = getEseq(bestd, c);

	// now follow the sequence
	x1 = work->pt1.X;
	z1 = work->pt1.Z;
	x2 = work->pt2.X;
	z2 = work->pt2.Z;
	x3 = work->pt3.X;
	z3 = work->pt3.Z;
	x4 = work->pt4.X;
	z4 = work->pt4.Z;
	s1 = work->sum1;
	s2 = work->sum2;
	d1 = work->diff1;
	d2 = work->diff2;

	// the first one is always a doubling
	// point1 is [1]P
	vecCopy(P->X, x1);
	vecCopy(P->Z, z1);
	vecCopy(P->X, x3);
	vecCopy(P->Z, z3);
	vecsubmod_ptr(P->X, P->Z, d1, mdata);
	vecaddmod_ptr(P->X, P->Z, s1, mdata);

	d = 1;
	x = 2;

	// point2 is [2]P
	vec_duplicate(mdata, work, s1, d1, &work->pt2);

#ifdef DEBUG
	printf("target is %lu\n", c);
	printf("pt1 holds [1]P\n");
	print_vechex(work->pt1.X->data, 0, NWORDS, "");
	print_vechex(work->pt1.Z->data, 0, NWORDS, "");
	printf("pt2 holds [2]P\n");
	print_vechex(work->pt2.X->data, 0, NWORDS, "");
	print_vechex(work->pt2.Z->data, 0, NWORDS, "");
#endif

	if (c == 2)
	{
		free(seq);
		vecCopy(x2, P->X);
		vecCopy(z2, P->Z);
		return;
	}
	
	startd = 0;
	i = 2;
	while (seq[i])
	{
		// both add and dup require the sum and difference 
		// of the X and Z coords of the two points we are tracking.
		vecaddsubmod_ptr(x2, z2, s2, d2, mdata);
		vecaddsubmod_ptr(x1, z1, s1, d1, mdata);

		if (seq[i] == 1)
		{
			// add point1 to point2, store in point1
			// double point2
			vec_add(mdata, work, &work->pt3, &work->pt1);
			vec_duplicate(mdata, work, s2, d2, &work->pt2);

			if (seq[i + 1] == 3)
			{
				d = x + d;
				x *= 2;
				i++;
			}
			else if (seq[i + 1] == 2)
			{
				d = x - d;
				x *= 2;

				// swap the "initial" point held in pt3 with the 
				// output point held in pt1
				sw_x = x3->data;
				sw_z = z3->data;
				x3->data = x1->data;
				z3->data = z1->data;
				x1->data = sw_x;
				z1->data = sw_z;
			}
			else
			{
				d = x + d;
				x *= 2;
			}

#ifdef DEBUG
			printf("[1]: pt1 holds [%d]P = ", d);
			print_vechex(work->pt1.X->data, 0, NWORDS, "");
			print_vechex(work->pt1.Z->data, 0, NWORDS, "");
			printf("[1]: pt2 holds [%d]P = ", x);
			print_vechex(work->pt2.X->data, 0, NWORDS, "");
			print_vechex(work->pt2.Z->data, 0, NWORDS, "");
#endif

		}
		else if(seq[i] == 2)
		{
			// add point1 (implicit in add/sub inputs) to point2, store in point2.
			// what is now in point 2 will either be copied to 
			// the new initial point (pt3) or to point 1.
			vecCopy(x2, x4);
			vecCopy(z2, z4);
			
#ifdef DEBUG
			printf("[2]: tmp holds [%d]P = ", x);
			print_vechex(work->pt4.X->data, 0, NWORDS, "");
			print_vechex(work->pt4.Z->data, 0, NWORDS, "");
#endif

			vec_add(mdata, work, &work->pt3, &work->pt2);
			x = x + d;

#ifdef DEBUG
			printf("[2]: pt2 holds [%d]P = ", x);
			print_vechex(work->pt2.X->data, 0, NWORDS, "");
			print_vechex(work->pt2.Z->data, 0, NWORDS, "");
#endif

			if (seq[i + 1] == 2)
			{
				vecCopy(x4, x3);
				vecCopy(z4, z3);

#ifdef DEBUG
				printf("no swap\n");
				printf("[2]: pt3 holds [%d]P = ", x - d);
				print_vechex(work->pt3.X->data, 0, NWORDS, "");
				print_vechex(work->pt3.Z->data, 0, NWORDS, "");
				printf("[2]: pt1 holds [%d]P = ", d);
				print_vechex(work->pt1.X->data, 0, NWORDS, "");
				print_vechex(work->pt1.Z->data, 0, NWORDS, "");
#endif
				
				
			}
			else
			{
				// change point1
				vecCopy(x1, x3);
				vecCopy(z1, z3);
				vecCopy(x4, x1);
				vecCopy(z4, z1);

#ifdef DEBUG
				printf("next is a swap\n");
				printf("[2]: pt3 holds [%d]P = ", d);
				print_vechex(work->pt3.X->data, 0, NWORDS, "");
				print_vechex(work->pt3.Z->data, 0, NWORDS, "");
				printf("[2]: pt1 holds [%d]P = ", x - d);
				print_vechex(work->pt1.X->data, 0, NWORDS, "");
				print_vechex(work->pt1.Z->data, 0, NWORDS, "");
#endif

				d = x - d;

				if (seq[i + 1] == 3)
					i++;
			}


		}

		i++;

#ifdef DEBUG
		printf("==================\n");
#endif
	}

	if (x != c)
	{
		printf("expected %lu, euclid returned %lu\n", c, x);
		exit(1);
	}

	// copy out answer
	vecCopy(x2, P->X);
	vecCopy(z2, P->Z);

#ifdef DEBUG
	//exit(1);
#endif

	free(seq);
	return;
}

void next_pt_vec(vec_monty_t *mdata, ecm_work *work, ecm_pt *P, uint64_t c)
{
	uint64_t mask;
	vec_bignum_t *x1, *z1, *x2, *z2, *s1, *s2, *d1, *d2;
	uint64_t e, d;

	x1 = work->pt1.X;
	z1 = work->pt1.Z;
	x2 = work->pt2.X;
	z2 = work->pt2.Z;
	s1 = work->sum1;
	s2 = work->sum2;
	d1 = work->diff1;
	d2 = work->diff2;

	if (c == 1)
	{
		return;
	}

	vecCopy(P->X, x1);
	vecCopy(P->Z, z1);
	vecsubmod_ptr(P->X, P->Z, d1, mdata);
	vecaddmod_ptr(P->X, P->Z, s1, mdata);
	vec_duplicate(mdata, work, s1, d1, &work->pt2);

	//mulcnt[tid] += 3;
	//sqrcnt[tid] += 2;

	if (c == 2)
	{
		vecCopy(x2, P->X);
		vecCopy(z2, P->Z);
		return;
	}

	// find the first '1' bit then skip it
#if defined( __INTEL_COMPILER) || defined (__INTEL_LLVM_COMPILER)
	mask = 1ULL << (64 - _lzcnt_u64((uint64_t)c) - 2);
#elif defined(__GNUC__)
	mask = 1ULL << (64 - __builtin_clzll((uint64_t)c) - 2);
#elif defined _MSC_VER
    mask = 1ULL << (64 - _lzcnt_u64((uint64_t)c) - 2);
#endif

	//goal is to compute x_c, z_c using montgomery's addition
	//and multiplication formula's for x and z
	//the procedure will be similar to p+1 lucas chain formula's
	//but this time we are simultaneously incrementing both x and z
	//rather than just Vn.  In each bit of the binary expansion of
	//c, we need just one addition and one duplication for both x and z.
	//compute loop
	//for each bit of M to the right of the most significant bit
	d = 1;
	e = 2;
	while (mask > 0)
	{
		vecaddsubmod_ptr(x2, z2, s2, d2, mdata);
		vecaddsubmod_ptr(x1, z1, s1, d1, mdata);
        
		//if the bit is 1
		if (c & mask)
		{
			//add x1,z1, duplicate x2,z2
			vec_add(mdata, work, P, &work->pt1);
			vec_duplicate(mdata, work, s2, d2, &work->pt2);
			d = d + e;
			e *= 2;
		}
		else
		{
			//add x2,z2, duplicate x1,z1
			vec_add(mdata, work, P, &work->pt2);
			vec_duplicate(mdata, work, s1, d1, &work->pt1);
			e = e + d;
			d *= 2;
		}

		//mulcnt[tid] += 7;
		//sqrcnt[tid] += 4;

		mask >>= 1;
	}
	if (d != c)
	{
		printf("expected %lu, ladder returned %lu\n", c, d);
		exit(1);
	}
	vecCopy(x1, P->X);
	vecCopy(z1, P->Z);

	return;
}

void array_mul(uint64_t *primes, uint64_t b1, int num_p, mpz_t piprimes)
{
	mpz_t *p;
	int i;
	int alloc;
	uint64_t stg1 = b1;

	if (num_p & 1)
	{
		p = (mpz_t *)malloc((num_p + 1) * sizeof(mpz_t));
		alloc = num_p + 1;
	}
	else
	{
		p = (mpz_t *)malloc((num_p + 0) * sizeof(mpz_t));
		alloc = num_p;
	}

	for (i = 0; i < num_p; i++)
	{
		uint64_t c = 1;
		uint64_t q;

		q = primes[i];
		do {
			c *= q;
		} while ((c * q) < stg1);

		mpz_init(p[i]);
		mpz_set_ui(p[i], q);
	}

	if (num_p & 1)
	{
		mpz_init(p[i]);
		mpz_set_ui(p[i], 1);
		num_p++;
	}

	while (num_p != 1)
	{
		for (i = 0; i < num_p / 2; i++)
		{
			mpz_mul(p[i], p[i], p[i + num_p / 2]);
		}

		num_p /= 2;

		if ((num_p > 1) && (num_p & 1))
		{
			mpz_set_ui(p[num_p], 1);
			num_p++;
		}
	}

	mpz_set(piprimes, p[0]);

	printf("piprimes has %lu bits\n\n", mpz_sizeinbase(piprimes, 2));

	for (i = 0; i < alloc; i++)
	{
		mpz_clear(p[i]);
	}
	free(p);

	return;
}

#if 0 
void vececm_old(thread_data_t *tdata)
{
	//attempt to factor n with the elliptic curve method
	//following brent and montgomery's papers, and CP's book
    tpool_t *tpool_data;
    uint32_t threads = tdata[0].total_threads;
    base_t retval;
	base_t i, j;
	int curve;
	int tid;
	FILE *save = NULL;
	char fname[80];
	char *wstr;
	int found = 0;
    int result;
    uint64_t num_found;
	vec_bignum_t *one = vecInit(words);
    mpz_t gmpt, gmpn;
    int verbose = tdata[0].verbose;
    mpz_t tmp_factor;
    uint32_t curves_run = 0;

	// timing variables
	struct timeval stopt;	// stop time of this job
	struct timeval startt;	// start time of this job
	struct timeval fullstopt;	// stop time of this job
	struct timeval fullstartt;	// start time of this job
	double t_time;

	gettimeofday(&fullstartt, NULL);
    mpz_init(gmpt);
    mpz_init(gmpn);
    mpz_init(tmp_factor);

    // in this function, gmpn is only used to for screen output and to 
    // check factors.  So if this is a Mersenne input, use the original
    // input number.  (all math is still done relative to the base Mersenne.)
    if (tdata[0].mdata->isMersenne)
    {
        extract_bignum_from_vec_to_mpz(gmpn, tdata[0].mdata->vnhat, 0, NWORDS);
    }
    else
    {
        extract_bignum_from_vec_to_mpz(gmpn, tdata[0].mdata->n, 0, NWORDS);
    }

	for (i = 0; i < VECLEN; i++)
	{
		one->data[i] = 1;
	}
    one->size = 1;

    if (!ecm_primes_initialized)
    {
        soe_staticdata_t* sdata = soe_init(0, 1, 32768);
        ecm_primes = soe_wrapper(sdata, 0, 100000000, 0, &ecm_nump, 0, 0);
        ecm_minp = ecm_primes[0];
        ecm_maxp = ecm_primes[ecm_nump - 1];
        soe_finalize(sdata);
        ecm_primes_initialized = 1;
    }
    
    tpool_data = tpool_setup(tdata[0].total_threads, NULL, NULL, &vec_ecm_sync,
        &vec_ecm_dispatch, tdata);

    for (i = 0; i < threads; i++)
        tpool_data[i].debug = 0;

    tpool_add_work_fcn(tpool_data, &vec_ecm_build_curve_work_fcn);
    tpool_add_work_fcn(tpool_data, &vec_ecm_stage1_work_fcn);
    tpool_add_work_fcn(tpool_data, &vec_ecm_stage2_init_work_fcn);
    tpool_add_work_fcn(tpool_data, &vec_ecm_stage2_work_fcn);

	for (curve = 0; curve < tdata[0].curves; curve += VECLEN)
	{
        uint64_t p;

        if (verbose > 1)
            printf("\n");

        if (verbose >= 0)
        {
            printf("ecm: %d/%d curves on C%d @ B1=%lu, B2=100*B1\r",
                (curve + VECLEN) * threads, tdata[0].curves * threads,
                (int)gmp_base10(gmpn), STAGE1_MAX);
        }
        else if (verbose > 1)
        {
            printf("ecm: %d/%d curves on C%d @ B1=%lu, B2=100*B1\n",
                (curve + VECLEN) * threads, tdata[0].curves * threads,
                (int)gmp_base10(gmpn), STAGE1_MAX);
        }

		gettimeofday(&startt, NULL);
        
        // parallel curve building        
        for (i = 0; i < threads; i++)
        {
			tdata[i].work->stg1Add = 0;
			tdata[i].work->stg1Doub = 0;
			tdata[i].work->last_pid = 0;
            tdata[i].phase_done = 0;
            tdata[i].ecm_phase = 0;
        }
        tpool_go(tpool_data);

		gettimeofday(&stopt, NULL);
		t_time = ytools_difftime (&startt, &stopt);
		
        if (verbose > 1)
		    printf("Building curves took %1.4f seconds.\n",t_time);

		// parallel stage 1
		gettimeofday(&startt, NULL);

        for (p = 0; p < STAGE1_MAX; p += PRIME_RANGE)
        {
			// first condition will fetch more primes if B1 is large
			// and thus stage 1 takes several iterations of PRIME_RANGE.
			// second condition resets primes array, if it is modified
			// by stage 2, before starting a new batch of curves.
			if ((tdata[0].work->last_pid == ecm_nump) || ((p == 0) && (ecm_minp != 2)))
			{
				if (ecm_primes != NULL) { free(ecm_primes); ecm_primes = NULL; };
                
                soe_staticdata_t* sdata = soe_init(0, 1, 32768);
                ecm_primes = soe_wrapper(sdata, p, 
                    MIN(STAGE2_MAX + 1000, p + (uint64_t)PRIME_RANGE), 
                    0, &ecm_nump, 0, 0);
                ecm_maxp = ecm_primes[ecm_nump - 1];
                soe_finalize(sdata);

                ecm_minp = ecm_primes[0];
                ecm_maxp = ecm_primes[ecm_nump - 1];

                if (verbose > 1)
				    printf("found %lu primes in range [%lu : %lu]\n", 
                        ecm_nump, ecm_minp, ecm_maxp);
			}

            for (i = 0; i < threads; i++)
            {
                tdata[i].phase_done = 0;
                tdata[i].ecm_phase = 1;
            }

            if (verbose > 1)
            {
                printf("commencing Stage 1 @ prime %lu\n", ecm_minp);
            }

            // start the stage 1 thread pool
            tpool_go(tpool_data);

            if (p < STAGE1_MAX)
            {
                // record results for this batch of primes
                if (tdata[0].save_b1) // && (save_intermediate)
                {
                    sprintf(fname, "save_b1_intermediate.txt");
                    save = fopen(fname, "a");
                }

                gettimeofday(&stopt, NULL);
                t_time = ytools_difftime(&startt, &stopt);
                if (verbose > 1)
                    printf("Stage 1 current elapsed time is %1.4f seconds\n", t_time);

                for (j = 0; j < threads; j++)
                {

                    // GMP-ECM wants X/Z.
                    // or, equivalently, X and Z listed separately.
                    vecmulmod_ptr(tdata[j].P->X, one, tdata[j].work->tt4, tdata[j].work->n,
                        tdata[j].work->tt2, tdata[j].mdata);

                    vecmulmod_ptr(tdata[j].P->Z, one, tdata[j].work->tt3, tdata[j].work->n,
                        tdata[j].work->tt2, tdata[j].mdata);

                    //print_vechexbignum(tdata[j].P->Z, "Z point\n");

                    for (i = 0; i < VECLEN; i++)
                    {
                        extract_bignum_from_vec_to_mpz(gmpt, tdata[j].P->Z, i, NWORDS);

                        if (mpz_cmp_ui(gmpt, 0) == 0)
                        {
                            printf("something failed: tid = %d, vec = %d has zero result\n", j, i);
                        }

                        result = vec_check_factor(gmpt, gmpn, tmp_factor);
                        if (result == 1)
                        {

                            //FILE *out = fopen("ecm_results.txt", "a");

                            if (verbose > 1)
                            {
                                gmp_printf("\nfound factor %Zd in stage 1 in thread %d, vec position %d, with sigma = ",
                                    tmp_factor, j, i);
                                printf("%"PRIu64"\n", tdata[j].sigma[i]);
                            }

                            //if (out != NULL)
                            //{
                            //	gmp_fprintf(out, "\nfound factor %Zd in stage 1 at curve %d, "
                            //		"in thread %d, vec position %d, with sigma = ",
                            //        tmp_factor, threads * curve + j * VECLEN + i, j, i);
                            //    fprintf(out, "%"PRIu64"\n", tdata[j].sigma[i]);
                            //	fclose(out);
                            //}
                            //
                            //fflush(stdout);
                            found = 1;

                            if (tdata[j].numfactors == 0)
                            {
                                tdata[j].factors = (avx_ecm_factor_t*)malloc(1 * sizeof(avx_ecm_factor_t));
                            }
                            else
                            {
                                tdata[j].factors = (avx_ecm_factor_t*)realloc(tdata[j].factors,
                                    (tdata[j].numfactors + 1) * sizeof(avx_ecm_factor_t));
                            }

                            tdata[j].factors[tdata[j].numfactors].thread_id = j;
                            tdata[j].factors[tdata[j].numfactors].stg_id = 1;
                            tdata[j].factors[tdata[j].numfactors].sigma = tdata[j].sigma[i];
                            tdata[j].factors[tdata[j].numfactors].vec_id = i;
                            tdata[j].factors[tdata[j].numfactors].curve_id =
                                (curve * threads) + (j * VECLEN) + i;
                            mpz_init(tdata[j].factors[tdata[j].numfactors].factor);
                            mpz_set(tdata[j].factors[tdata[j].numfactors++].factor, tmp_factor);

                            //gmp_printf("thread %d now has %d factors in stg1: %Zd\n",
                            //    j, tdata[j].numfactors, tmp_factor);
                        }

                        if (tdata[0].save_b1)
                        {
                            fprintf(save, "METHOD=ECM; SIGMA=%"PRIu64"; B1=%"PRIu64"; ",
                                tdata[j].sigma[i], ecm_primes[tdata[j].work->last_pid - 1]);
                            gmp_fprintf(save, "N=0x%Zx; ", gmpn);

                            extract_bignum_from_vec_to_mpz(gmpt, tdata[j].work->tt4, i, NWORDS);
                            gmp_fprintf(save, "X=0x%Zx; ", gmpt);

                            extract_bignum_from_vec_to_mpz(gmpt, tdata[j].work->tt3, i, NWORDS);
                            gmp_fprintf(save, "Z=0x%Zx; PROGRAM=AVX-ECM;\n", gmpt);
                        }
                    }

                }

                if (tdata[0].save_b1)
                {
                    fclose(save);
                }
            }
        }

        // record results for this batch of primes
        if (tdata[0].save_b1) // && (save_intermediate)
        {
            sprintf(fname, "save_b1.txt");
            save = fopen(fname, "a");
        }

        for (j = 0; j < threads; j++)
        {

            // GMP-ECM wants X/Z.
            // or, equivalently, X and Z listed separately.
            vecmulmod_ptr(tdata[j].P->X, one, tdata[j].work->tt4, tdata[j].work->n,
                tdata[j].work->tt2, tdata[j].mdata);

            vecmulmod_ptr(tdata[j].P->Z, one, tdata[j].work->tt3, tdata[j].work->n,
                tdata[j].work->tt2, tdata[j].mdata);

            //print_vechexbignum(tdata[j].P->Z, "Z point\n");

            for (i = 0; i < VECLEN; i++)
            {
                extract_bignum_from_vec_to_mpz(gmpt, tdata[j].P->Z, i, NWORDS);

                if (mpz_cmp_ui(gmpt, 0) == 0)
                {
                    printf("something failed: tid = %d, vec = %d has zero result\n", j, i);
                }

                result = vec_check_factor(gmpt, gmpn, tmp_factor);
                if (result == 1)
                {

                    //FILE *out = fopen("ecm_results.txt", "a");

                    if (verbose > 1)
                    {
                        gmp_printf("\nfound factor %Zd in stage 1 in thread %d, vec position %d, with sigma = ",
                            tmp_factor, j, i);
                        printf("%"PRIu64"\n", tdata[j].sigma[i]);
                    }

                    //if (out != NULL)
                    //{
                    //	gmp_fprintf(out, "\nfound factor %Zd in stage 1 at curve %d, "
                    //		"in thread %d, vec position %d, with sigma = ",
                    //        tmp_factor, threads * curve + j * VECLEN + i, j, i);
                    //    fprintf(out, "%"PRIu64"\n", tdata[j].sigma[i]);
                    //	fclose(out);
                    //}
                    //
                    //fflush(stdout);
                    found = 1;

                    if (tdata[j].numfactors == 0)
                    {
                        tdata[j].factors = (avx_ecm_factor_t*)malloc(1 * sizeof(avx_ecm_factor_t));
                    }
                    else
                    {
                        tdata[j].factors = (avx_ecm_factor_t*)realloc(tdata[j].factors,
                            (tdata[j].numfactors + 1) * sizeof(avx_ecm_factor_t));
                    }

                    tdata[j].factors[tdata[j].numfactors].thread_id = j;
                    tdata[j].factors[tdata[j].numfactors].stg_id = 1;
                    tdata[j].factors[tdata[j].numfactors].sigma = tdata[j].sigma[i];
                    tdata[j].factors[tdata[j].numfactors].vec_id = i;
                    tdata[j].factors[tdata[j].numfactors].curve_id =
                        (curve * threads) + (j * VECLEN) + i;
                    mpz_init(tdata[j].factors[tdata[j].numfactors].factor);
                    mpz_set(tdata[j].factors[tdata[j].numfactors++].factor, tmp_factor);

                    //gmp_printf("thread %d now has %d factors in stg1: %Zd\n",
                    //    j, tdata[j].numfactors, tmp_factor);
                }

                if (tdata[0].save_b1)
                {
                    fprintf(save, "METHOD=ECM; SIGMA=%"PRIu64"; B1=%"PRIu64"; ",
                        tdata[j].sigma[i], ecm_primes[tdata[j].work->last_pid - 1]);
                    gmp_fprintf(save, "N=0x%Zx; ", gmpn);

                    extract_bignum_from_vec_to_mpz(gmpt, tdata[j].work->tt4, i, NWORDS);
                    gmp_fprintf(save, "X=0x%Zx; ", gmpt);

                    extract_bignum_from_vec_to_mpz(gmpt, tdata[j].work->tt3, i, NWORDS);
                    gmp_fprintf(save, "Z=0x%Zx; PROGRAM=AVX-ECM;\n", gmpt);
                }
            }

        }

        if (tdata[0].save_b1)
        {
            fclose(save);
        }

        gettimeofday(&stopt, NULL);
        t_time = ytools_difftime(&startt, &stopt);
        if (verbose > 1)
            printf("Stage 1 took %1.4f seconds\n", t_time);
	
        curves_run += VECLEN * threads;

        // always stop when a factor is found
        //if (found)
        //	break;


        if (DO_STAGE2)
        {
            uint64_t last_p = ecm_primes[tdata[0].work->last_pid];

            // parallel stage 2
            gettimeofday(&startt, NULL);

			// stage 2 parallel init
            for (i = 0; i < threads; i++)
            {
                tdata[i].phase_done = 0;
                tdata[i].ecm_phase = 2;
            }
            tpool_go(tpool_data);

            for (; last_p < STAGE2_MAX; )
            {
				if (last_p == ecm_maxp)
				{
                    if (ecm_primes != NULL) { free(ecm_primes); ecm_primes = NULL; };

                    soe_staticdata_t* sdata = soe_init(0, 1, 32768);
                    ecm_primes = soe_wrapper(sdata, last_p, 
                        MIN(last_p + (uint64_t)PRIME_RANGE, STAGE2_MAX + 1000),
                        0, &ecm_nump, 0, 0);
                    ecm_maxp = ecm_primes[ecm_nump - 1];
                    soe_finalize(sdata);

                    ecm_minp = ecm_primes[0];
                    ecm_maxp = ecm_primes[ecm_nump - 1];

					for (i = 0; i < threads; i++)
					{
						tdata[i].work->last_pid = 1;
					}

                    if (verbose > 1)
                    {
                        printf("found %lu primes in range [%lu : %lu]\n",
                            ecm_nump, ecm_minp, ecm_maxp);
                    }
				}

                for (i = 0; i < threads; i++)
                {
                    tdata[i].phase_done = 0;
                    tdata[i].ecm_phase = 3;
                }
                tpool_go(tpool_data);

                //printf("stage 2 finished at P=%lu (maxP = %lu) (id %u of %lu)\n",
                //    PRIMES[tdata[0].work->last_pid], P_MAX, tdata[0].work->last_pid, NUM_P);

				if (tdata[0].work->last_pid == ecm_nump)
				{
					// we ended at the last prime we cached.  Check if we
					// need to do more.  
					if (STAGE2_MAX == PRIME_RANGE)
						break;
					else
						last_p = ecm_maxp;
				}
				else
				{
					last_p = ecm_primes[tdata[0].work->last_pid];
				}
            }

            gettimeofday(&stopt, NULL);
            t_time = ytools_difftime(&startt, &stopt);
            if (verbose > 1)
            {
                printf("Stage 2 took %1.4f seconds\n", t_time);
                printf("performed %d pair-multiplies for %lu primes in stage 2\n",
                    tdata[0].work->paired, tdata[0].work->numprimes);
                printf("performed %u point-additions and %u point-doubles in stage 2\n",
                    tdata[0].work->ptadds + tdata[0].work->stg1Add, tdata[0].work->stg1Doub);
            }

            for (j = 0; j < threads; j++)
            {
				for (i = 0; i < VECLEN; i++)
                {
                    extract_bignum_from_vec_to_mpz(gmpt, tdata[j].work->stg2acc, i, NWORDS);
                    result = vec_check_factor(gmpt, gmpn, tmp_factor);

                    if (mpz_cmp_ui(gmpt, 0) == 0)
                    {
                        printf("something failed: tid = %d, vec = %d has zero result\n", j, i);
                    }

                    if (result == 1)
                    {
						//FILE *out = fopen("ecm_results.txt", "a");

                        if (verbose > 1)
                        {
                            gmp_printf("\nfound factor %Zd in stage 2 in thread %d, vec position %d, with sigma = ",
                                tmp_factor, j, i);
                            printf("%"PRIu64"\n", tdata[j].sigma[i]);
                        }

						//if (out != NULL)
						//{
						//	gmp_fprintf(out, "\nfound factor %Zd in stage 1 at curve %d, "
						//		"in thread %d, vec position %d, with sigma = ",
                        //        tmp_factor, threads * curve + j * VECLEN + i, j, i);
                        //    fprintf(out, "%"PRIu64"\n", tdata[j].sigma[i]);
						//	fclose(out);
						//}
                        //
                        //fflush(stdout);
                        found = 1;

                        if (tdata[j].numfactors == 0)
                        {
                            tdata[j].factors = (avx_ecm_factor_t*)malloc(1 * sizeof(avx_ecm_factor_t));
                        }
                        else
                        {
                            tdata[j].factors = (avx_ecm_factor_t*)realloc(tdata[j].factors,
                                (tdata[j].numfactors + 1) * sizeof(avx_ecm_factor_t));
                        }

                        tdata[j].factors[tdata[j].numfactors].thread_id = j;
                        tdata[j].factors[tdata[j].numfactors].stg_id = 2;
                        tdata[j].factors[tdata[j].numfactors].sigma = tdata[j].sigma[i];
                        tdata[j].factors[tdata[j].numfactors].vec_id = i;
                        tdata[j].factors[tdata[j].numfactors].curve_id =
                            (curve * threads) + (j * VECLEN) + i;
                        mpz_init(tdata[j].factors[tdata[j].numfactors].factor);
                        mpz_set(tdata[j].factors[tdata[j].numfactors++].factor, tmp_factor);

                        //gmp_printf("thread %d now has %d factors in stg2: %Zd\n", 
                        //    j, tdata[j].numfactors, tmp_factor);
                    }
                }
            }
        }

        if (verbose >= 0)
        {
            printf("ecm: %d/%d curves on C%d @ B1=%lu, B2=100*B1\r",
                (curve + VECLEN) * threads, tdata[0].curves * threads,
                (int)gmp_base10(gmpn), STAGE1_MAX);
        }

		if (found)
			break;

	}

    if (verbose > 1)
    {
        printf("ecm: %d/%d curves on C%d @ B1=%lu, B2=100*B1\n",
            curves_run, tdata[0].curves * threads,
            (int)gmp_base10(gmpn), STAGE1_MAX);
    }
    else
    {
        printf("\n");
    }

    tdata[0].curves = curves_run;

	gettimeofday(&fullstopt, NULL);
	t_time = ytools_difftime(&fullstartt, &fullstopt);
    if (verbose > 1)
	    printf("Process took %1.4f seconds.\n", t_time);

	vecClear(one);
    mpz_clear(gmpn);
    mpz_clear(gmpt);
    mpz_clear(tmp_factor);
	return;
}
#endif

void vececm(thread_data_t* tdata)
{
    // attempt to factor n with the elliptic curve method
    // following brent and montgomery's papers, and CP's book
    tpool_t* tpool_data;
    uint32_t threads = tdata[0].total_threads;
    base_t i, j;
    int curve;
    FILE* save = NULL;
    char fname[80];
    int found = 0;
    int result;
    uint64_t num_found;
    mpz_t gmpt, gmpn, tmp_factor;
    // these track the range over which we currently have a prime list.
    uint64_t rangemin;
    uint64_t rangemax;
    uint64_t sigma_in[VECLEN];
    int verbose = tdata[0].verbose;
    uint32_t NWORDS = tdata[0].NWORDS;
    uint32_t NBLOCKS = tdata[0].NBLOCKS;
    uint32_t MAXBITS = tdata[0].MAXBITS;
    uint64_t STAGE1_MAX = tdata[0].STAGE1_MAX;
    uint64_t STAGE2_MAX = tdata[0].STAGE2_MAX;
    int DO_STAGE2 = tdata[0].DO_STAGE2;
    uint32_t PRIME_RANGE = tdata[0].PRIME_RANGE;
    vec_bignum_t* one = vecInit(NWORDS);
    uint32_t curves_to_run = tdata[0].curves;

    // timing variables
    struct timeval stopt;	// stop time of this job
    struct timeval startt;	// start time of this job
    struct timeval fullstopt;	// stop time of this job
    struct timeval fullstartt;	// start time of this job
    double t_time;

    // primes, allocated per-process
    uint64_t* ecm_primes;
    uint64_t ecm_nump;
    uint64_t ecm_minp;
    uint64_t ecm_maxp;

    gettimeofday(&fullstartt, NULL);
    mpz_init(gmpt);
    mpz_init(gmpn);
    mpz_init(tmp_factor);

    // in this function, gmpn is only used to for screen output and to 
    // check factors.  So if this is a Mersenne input, use the original
    // input number.  (all math is still done relative to the base Mersenne.)
    if (tdata[0].mdata->isMersenne == 0)
    {
        extract_bignum_from_vec_to_mpz(gmpn, tdata[0].mdata->n, 0, NWORDS);
    }
    else
    {
        // for Mersenne's, the original input is stored here.
        extract_bignum_from_vec_to_mpz(gmpn, tdata[0].mdata->vnhat, 0, NWORDS);
    }

    for (i = 0; i < VECLEN; i++)
    {
        one->data[i] = 1;
    }
    one->size = 1;


    rangemin = 0;
    rangemax = MIN(STAGE2_MAX + 1000, (uint64_t)PRIME_RANGE);

    soe_staticdata_t* sdata = soe_init(0, 1, 32768);
    ecm_primes = soe_wrapper(sdata, rangemin, rangemax, 0, &ecm_nump, 0, 0);
    ecm_minp = ecm_primes[0];
    ecm_maxp = ecm_primes[ecm_nump - 1];

    soe_finalize(sdata);

    if (verbose > 1)
    {
        printf("found %lu primes in range [%lu : %lu]\n", ecm_nump, rangemin, rangemax);
    }

    tpool_data = tpool_setup(tdata[0].total_threads, NULL, NULL, &vec_ecm_sync,
        &vec_ecm_dispatch, tdata);

    tpool_add_work_fcn(tpool_data, &vec_ecm_build_curve_work_fcn);
    tpool_add_work_fcn(tpool_data, &vec_ecm_stage1_work_fcn);
    tpool_add_work_fcn(tpool_data, &vec_ecm_stage2_init_work_fcn);
    tpool_add_work_fcn(tpool_data, &vec_ecm_stage2_work_fcn);

    for (j = 0; j < VECLEN; j++)
    {
        sigma_in[j] = tdata[0].sigma[j];
    }

    // set up recording results for this batch of primes if requested.
    if (tdata[0].save_b1) // && (save_intermediate)
    {
        if (tdata[0].save_b1 == 2)
        {
            remove("avx_ecm_resume.txt");
            sprintf(fname, "avx_ecm_resume.txt");
        }
        else
        {
            sprintf(fname, "save_b1.txt");
        }
        save = fopen(fname, "a");
    }

    for (curve = 0; curve < curves_to_run; curve += VECLEN)
    {
        uint64_t p;

        if (verbose >= 0)
        {
            printf("ecm: %d/%d curves on C%d @ B1=%lu, B2=100*B1\r",
                (curve + VECLEN) * threads, tdata[0].curves * threads,
                (int)gmp_base10(gmpn), STAGE1_MAX); fflush(stdout);
        }

        // get a new batch of primes if:
        // - the current range starts after B1
        // - the current range ends after B1 AND the range isn't maximal length
        if ((rangemin > STAGE1_MAX) ||
            ((rangemax < STAGE1_MAX) && ((rangemax - rangemin) < PRIME_RANGE)))
        {
            rangemin = 0;
            rangemax = MIN(STAGE2_MAX + 1000, (uint64_t)PRIME_RANGE);

            soe_staticdata_t* sdata = soe_init(0, 1, 32768);
            if (ecm_primes != NULL) {
                free(ecm_primes);
            }
            sdata->VFLAG = verbose;
            ecm_primes = soe_wrapper(sdata, rangemin, rangemax, 0, &ecm_nump, 0, 0);
            ecm_minp = ecm_primes[0];
            ecm_maxp = ecm_primes[ecm_nump - 1];
            soe_finalize(sdata);

            if (verbose > 1)
            {
                printf("Found %lu primes in range [%lu : %lu]\n", ecm_nump, rangemin, rangemax);
            }
        }

        gettimeofday(&startt, NULL);

        // parallel curve building        
        for (i = 0; i < threads; i++)
        {
            tdata[i].work->ptadds = 0;
            tdata[i].work->ptdups = 0;
            tdata[i].work->last_pid = 0;
            tdata[i].phase_done = 0;
            tdata[i].ecm_phase = 0;
            tdata[i].primes = ecm_primes;
            tdata[i].nump = ecm_nump;
            tdata[i].minp = ecm_minp;
            tdata[i].maxp = ecm_maxp;

            for (j = 0; j < VECLEN; j++)
            {
                if (sigma_in[j] > 0)
                {
                    tdata[i].sigma[j] = sigma_in[j] + curve;
                }
                else
                {
                    tdata[i].sigma[j] = 0;
                }
            }
        }
        tpool_go(tpool_data);

        gettimeofday(&stopt, NULL);
        t_time = ytools_difftime(&startt, &stopt);

        if (verbose > 1)
        {
            printf("\n");

            printf("Commencing curves %d-%d of %u\n", threads * curve,
                threads * (curve + VECLEN) - 1, threads * tdata[0].curves);

            printf("Building curves took %1.4f seconds.\n", t_time);
        }

        // parallel stage 1
        gettimeofday(&startt, NULL);

        if (STAGE1_MAX > PRIME_RANGE)
        {
            // this is a large B1:
            // prepare to record intermediate results for this batch of primes
            remove("save_b1_intermediate.txt");
        }

        for (p = 0; p < STAGE1_MAX; p += PRIME_RANGE)
        {
            // get a new batch of primes if the current range ends
            // before this new one starts
            if (p >= rangemax)
            {
                rangemin = rangemax;
                rangemax = MIN(STAGE2_MAX + 1000, rangemin + (uint64_t)PRIME_RANGE);

                soe_staticdata_t* sdata = soe_init(0, 1, 32768);
                if (ecm_primes != NULL) {
                    free(ecm_primes);
                }
                ecm_primes = soe_wrapper(sdata, rangemin, rangemax, 0, &ecm_nump, 0, 0);
                ecm_minp = ecm_primes[0];
                ecm_maxp = ecm_primes[ecm_nump - 1];
                soe_finalize(sdata);

                if (verbose > 1)
                {
                    printf("Found %lu primes in range [%lu : %lu]\n", ecm_nump, rangemin, rangemax);
                }
            }

            for (i = 0; i < threads; i++)
            {
                tdata[i].phase_done = 0;
                tdata[i].ecm_phase = 1;
                tdata[i].primes = ecm_primes;
                tdata[i].nump = ecm_nump;
                tdata[i].minp = ecm_minp;
                tdata[i].maxp = ecm_maxp;
            }

            if (verbose > 1)
            {
                printf("Commencing Stage 1 @ prime %lu\n", ecm_minp);
            }

            tpool_go(tpool_data);

#if 1
            if (STAGE1_MAX > PRIME_RANGE)
            {
                // this is a large B1:
                // record intermediate results for this batch of primes
                sprintf(fname, "save_b1_intermediate.txt");
                FILE *intsave = fopen(fname, "a");

                gettimeofday(&stopt, NULL);
                t_time = ytools_difftime(&startt, &stopt);
                if (verbose > 1)
                {
                    printf("Stage 1 current elapsed time is %1.4f seconds\n", t_time);
                }

                for (j = 0; j < threads; j++)
                {

                    // GMP-ECM wants X/Z.
                    // or, equivalently, X and Z listed separately.
                    vecmulmod_ptr(tdata[j].P->X, one, tdata[j].work->tt4, tdata[j].work->n,
                        tdata[j].work->tt2, tdata[j].mdata);

                    vecmulmod_ptr(tdata[j].P->Z, one, tdata[j].work->tt3, tdata[j].work->n,
                        tdata[j].work->tt2, tdata[j].mdata);

                    //print_vechexbignum(tdata[j].P->Z, "Z point\n");

                    for (i = 0; i < VECLEN; i++)
                    {
                        extract_bignum_from_vec_to_mpz(gmpt, tdata[j].P->Z, i, NWORDS);

                        if (mpz_cmp_ui(gmpt, 0) == 0)
                        {
                            printf("something failed: tid = %d, vec = %d has zero result\n", j, i);
                        }

                        result = vec_check_factor(gmpt, gmpn, tmp_factor);
                        if (result == 1)
                        {

                            //FILE *out = fopen("ecm_results.txt", "a");

                            if (verbose > 1)
                            {
                                gmp_printf("\nfound factor %Zd in stage 1 in thread %d, vec position %d, with sigma = ",
                                    tmp_factor, j, i);
                                printf("%"PRIu64"\n", tdata[j].sigma[i]);
                            }

                            //if (out != NULL)
                            //{
                            //	gmp_fprintf(out, "\nfound factor %Zd in stage 1 at curve %d, "
                            //		"in thread %d, vec position %d, with sigma = ",
                            //        tmp_factor, threads * curve + j * VECLEN + i, j, i);
                            //    fprintf(out, "%"PRIu64"\n", tdata[j].sigma[i]);
                            //	fclose(out);
                            //}
                            //
                            //fflush(stdout);
                            found = 1;

                            if (tdata[j].numfactors == 0)
                            {
                                tdata[j].factors = (avx_ecm_factor_t*)malloc(1 * sizeof(avx_ecm_factor_t));
                            }
                            else
                            {
                                tdata[j].factors = (avx_ecm_factor_t*)realloc(tdata[j].factors,
                                    (tdata[j].numfactors + 1) * sizeof(avx_ecm_factor_t));
                            }

                            tdata[j].factors[tdata[j].numfactors].thread_id = j;
                            tdata[j].factors[tdata[j].numfactors].stg_id = 1;
                            tdata[j].factors[tdata[j].numfactors].sigma = tdata[j].sigma[i];
                            tdata[j].factors[tdata[j].numfactors].vec_id = i;
                            tdata[j].factors[tdata[j].numfactors].curve_id =
                                (curve * threads) + (j * VECLEN) + i;
                            mpz_init(tdata[j].factors[tdata[j].numfactors].factor);
                            mpz_set(tdata[j].factors[tdata[j].numfactors++].factor, tmp_factor);

                            //gmp_printf("thread %d now has %d factors in stg1: %Zd\n",
                            //    j, tdata[j].numfactors, tmp_factor);
                        }

                        if (tdata[0].save_b1)
                        {
                            mpz_t tmp1, tmp2;
                            mpz_init(tmp2);
                            mpz_init(tmp1);

                            fprintf(intsave, "METHOD=ECM; SIGMA=%"PRIu64"; B1=%"PRIu64"; ",
                                tdata[j].sigma[i], ecm_primes[tdata[j].work->last_pid - 1]);
                            gmp_fprintf(intsave, "N=0x%Zx; ", gmpn);

                            extract_bignum_from_vec_to_mpz(tmp1, tdata[j].work->tt4, i, NWORDS);
                            extract_bignum_from_vec_to_mpz(tmp2, tdata[j].work->tt3, i, NWORDS);

                            mpz_invert(tmp2, tmp2, gmpn);
                            mpz_mul(gmpt, tmp2, tmp1);
                            mpz_mod(gmpt, gmpt, gmpn);

                            gmp_fprintf(intsave, "X=0x%Zx; PROGRAM=AVX-ECM;\n", gmpt);

                            mpz_clear(tmp1);
                            mpz_clear(tmp2);
                        }
                    }

                }

                fclose(intsave);

            }
#endif
        }

        gettimeofday(&stopt, NULL);
        t_time = ytools_difftime(&startt, &stopt);

        if (verbose > 1)
        {
            printf("Stage 1 took %1.4f seconds\n", t_time);
        }

        

        for (j = 0; j < threads; j++)
        {

            // GMP-ECM wants X/Z.
            // or, equivalently, X and Z listed separately.
            vecmulmod_ptr(tdata[j].P->X, one, tdata[j].work->tt4, tdata[j].work->n,
                tdata[j].work->tt2, tdata[j].mdata);

            vecmulmod_ptr(tdata[j].P->Z, one, tdata[j].work->tt3, tdata[j].work->n,
                tdata[j].work->tt2, tdata[j].mdata);
            //print_vechexbignum(tdata[j].P->Z, "Z point\n");

            for (i = 0; i < VECLEN; i++)
            {
                //uint32_t msk = vec_gte52(tdata[j].work->tt3, tdata[j].work->n);
                //vec_bignum52_mask_sub(tdata[j].work->tt3, tdata[j].work->n, tdata[j].work->tt2, msk);
                //extract_bignum_from_vec_to_mpz(gmpt, tdata[j].work->tt2, i, NWORDS);
                extract_bignum_from_vec_to_mpz(gmpt, tdata[j].P->Z, i, NWORDS);

                if (mpz_cmp_ui(gmpt, 0) == 0)
                {
                    printf("something failed: tid = %d, vec = %d has zero result\n", j, i);
                }

#ifdef MPZ_VERIFY
                gmp_printf("final Z coord is %Zx\n", tdata[j].mdata->gmp_t2);
                mpz_gcd(gmpt, tdata[j].mdata->gmp_t2, gmpn);
                gmp_printf("gcd is %Zd\n", gmpt);
#endif

                result = vec_check_factor(gmpt, gmpn, tmp_factor);
                if (result == 1)
                {

                    //FILE *out = fopen("ecm_results.txt", "a");

                    if (verbose > 1)
                    {
                        gmp_printf("\nfound factor %Zd in stage 1 in thread %d, vec position %d, with sigma = ",
                            tmp_factor, j, i);
                        printf("%"PRIu64"\n", tdata[j].sigma[i]);
                    }

                    //if (out != NULL)
                    //{
                    //	gmp_fprintf(out, "\nfound factor %Zd in stage 1 at curve %d, "
                    //		"in thread %d, vec position %d, with sigma = ",
                    //        tmp_factor, threads * curve + j * VECLEN + i, j, i);
                    //    fprintf(out, "%"PRIu64"\n", tdata[j].sigma[i]);
                    //	fclose(out);
                    //}
                    //
                    //fflush(stdout);
                    found = 1;

                    if (tdata[j].numfactors == 0)
                    {
                        tdata[j].factors = (avx_ecm_factor_t*)malloc(1 * sizeof(avx_ecm_factor_t));
                    }
                    else
                    {
                        tdata[j].factors = (avx_ecm_factor_t*)realloc(tdata[j].factors,
                            (tdata[j].numfactors + 1) * sizeof(avx_ecm_factor_t));
                    }

                    tdata[j].factors[tdata[j].numfactors].thread_id = j;
                    tdata[j].factors[tdata[j].numfactors].stg_id = 1;
                    tdata[j].factors[tdata[j].numfactors].sigma = tdata[j].sigma[i];
                    tdata[j].factors[tdata[j].numfactors].vec_id = i;
                    tdata[j].factors[tdata[j].numfactors].curve_id =
                        (curve * threads) + (j * VECLEN) + i;
                    mpz_init(tdata[j].factors[tdata[j].numfactors].factor);
                    mpz_set(tdata[j].factors[tdata[j].numfactors].factor, tmp_factor);

                    tdata[j].numfactors++;
                    //gmp_printf("thread %d now has %d factors in stg1: %Zd\n",
                    //    j, tdata[j].numfactors, tmp_factor);
                }

                if (tdata[0].save_b1)
                {
                    mpz_t tmp1, tmp2;
                    mpz_init(tmp2);
                    mpz_init(tmp1);

                    fprintf(save, "METHOD=ECM; SIGMA=%"PRIu64"; B1=%"PRIu64"; ",
                        tdata[j].sigma[i], STAGE1_MAX);
                    gmp_fprintf(save, "N=0x%Zx; ", gmpn);

                    extract_bignum_from_vec_to_mpz(tmp1, tdata[j].work->tt4, i, NWORDS);
                    extract_bignum_from_vec_to_mpz(tmp2, tdata[j].work->tt3, i, NWORDS);

                    mpz_invert(tmp2, tmp2, gmpn);
                    mpz_mul(gmpt, tmp2, tmp1);
                    mpz_mod(gmpt, gmpt, gmpn);

                    gmp_fprintf(save, "X=0x%Zx; PROGRAM=AVX-ECM;\n", gmpt);
                    //gmp_fprintf(save, "Z=0x%Zx; PROGRAM=AVX-ECM;\n", gmpt);

                    mpz_clear(tmp1);
                    mpz_clear(tmp2);
                }
            }

        }

        // always stop when a factor is found
        //if (found)
        //	break;

        if (DO_STAGE2)
        {
            uint64_t last_p = ecm_primes[tdata[0].work->last_pid];

            // parallel stage 2
            gettimeofday(&startt, NULL);

            // stage 2 parallel init
            for (i = 0; i < threads; i++)
            {
                tdata[i].phase_done = 0;
                tdata[i].ecm_phase = 2;
            }
            tpool_go(tpool_data);

            for (i = 0; i < threads; i++)
            {
                if (tdata[i].work->last_pid == (uint32_t)(-1))
                {
                    // found a factor while initializing stage 2
                    printf("received factor from stage 2 init\n");
                    //print_vechexbignum52(tdata[i].work->stg2acc, "stg2acc: ");
                    last_p = STAGE2_MAX;
                }
            }

            gettimeofday(&stopt, NULL);
            t_time = ytools_difftime(&startt, &stopt);

            if (verbose > 1)
            {
                printf("Stage 2 Init took %1.4f seconds\n", t_time);
            }

            for (p = STAGE1_MAX; p < STAGE2_MAX; p += PRIME_RANGE)
            {
                // get a new batch of primes if the current range ends
                // before this new one starts or if the range isn't large
                // enough to cover the current step
                // TODO: make the primes range and the pair generator line up
                // in such a way that new ranges can pick up where the last
                // left off, without recomputing 2 * L Pa's and re-inverting.
                // it is not a large, or even a medium expense on typical B1/B2's,
                // but is still wasteful.
                if ((p >= rangemax) ||
                    (rangemax < MIN(STAGE2_MAX, p + (uint64_t)PRIME_RANGE)))
                {
                    rangemin = p;
                    rangemax = MIN(STAGE2_MAX + 1000, p + (uint64_t)PRIME_RANGE);

                    soe_staticdata_t* sdata = soe_init(0, 1, 32768);
                    if (ecm_primes != NULL) {
                        free(ecm_primes);
                    }
                    sdata->VFLAG = verbose;
                    ecm_primes = soe_wrapper(sdata, rangemin, rangemax, 0, &ecm_nump, 0, 0);
                    ecm_minp = ecm_primes[0];
                    ecm_maxp = ecm_primes[ecm_nump - 1];
                    soe_finalize(sdata);

                    for (i = 0; i < threads; i++)
                    {
                        tdata[i].work->last_pid = 1;
                        tdata[i].primes = ecm_primes;
                        tdata[i].nump = ecm_nump;
                        tdata[i].minp = ecm_minp;
                        tdata[i].maxp = ecm_maxp;
                    }

                    if (verbose > 1)
                    {
                        printf("found %lu primes in range [%lu : %lu]\n", ecm_nump, rangemin, rangemax);
                    }
                }

                tdata[0].pairmap_steps = pair(tdata[0].pairmap_v, tdata[0].pairmap_u,
                    tdata[0].work, tdata[0].Q, tdata[0].Qrmap, tdata[0].Qmap,
                    ecm_primes, ecm_nump, p, MIN(p + (uint64_t)PRIME_RANGE, STAGE2_MAX), verbose);

                for (i = 0; i < threads; i++)
                {
                    tdata[i].pairmap_steps = tdata[0].pairmap_steps;
                    tdata[i].work->amin = tdata[0].work->amin;
                    tdata[i].phase_done = 0;
                    tdata[i].ecm_phase = 3;
                }
                tpool_go(tpool_data);

                //printf("\nlast amin: %u\n", tdata[0].work->amin);

                if (tdata[0].work->last_pid == ecm_nump)
                {
                    // we ended at the last prime we cached.  Check if we
                    // need to do more.  
                    if (STAGE2_MAX == PRIME_RANGE)
                        break;
                    else
                        last_p = ecm_maxp;
                }
                else
                {
                    last_p = ecm_primes[tdata[0].work->last_pid];
                }
            }

            gettimeofday(&stopt, NULL);
            t_time = ytools_difftime(&startt, &stopt);
            if (verbose > 1)
            {
                printf("Stage 2 took %1.4f seconds\n", t_time);
                printf("performed %u pt-adds, %u inversions, and %u pair-muls in stage 2\n",
                    tdata[0].work->ptadds, tdata[0].work->numinv, tdata[0].work->paired);
            }

            for (j = 0; j < threads; j++)
            {
                for (i = 0; i < VECLEN; i++)
                {
                    extract_bignum_from_vec_to_mpz(gmpt, tdata[j].work->stg2acc, i, NWORDS);
                    result = vec_check_factor(gmpt, gmpn, tmp_factor);

                    if (mpz_cmp_ui(gmpt, 0) == 0)
                    {
                        printf("something failed: tid = %d, vec = %d has zero result\n", j, i);
                    }

                    if (result == 1)
                    {
                        //FILE *out = fopen("ecm_results.txt", "a");

                        if (verbose > 1)
                        {
                            gmp_printf("\nfound factor %Zd in stage 2 in thread %d, vec position %d, with sigma = ",
                                tmp_factor, j, i);
                            printf("%"PRIu64"\n", tdata[j].sigma[i]);
                        }

                        //if (out != NULL)
                        //{
                        //	gmp_fprintf(out, "\nfound factor %Zd in stage 1 at curve %d, "
                        //		"in thread %d, vec position %d, with sigma = ",
                        //        tmp_factor, threads * curve + j * VECLEN + i, j, i);
                        //    fprintf(out, "%"PRIu64"\n", tdata[j].sigma[i]);
                        //	fclose(out);
                        //}
                        //
                        //fflush(stdout);
                        found = 1;

                        if (tdata[j].numfactors == 0)
                        {
                            tdata[j].factors = (avx_ecm_factor_t*)malloc(1 * sizeof(avx_ecm_factor_t));
                        }
                        else
                        {
                            tdata[j].factors = (avx_ecm_factor_t*)realloc(tdata[j].factors,
                                (tdata[j].numfactors + 1) * sizeof(avx_ecm_factor_t));
                        }

                        tdata[j].factors[tdata[j].numfactors].thread_id = j;
                        tdata[j].factors[tdata[j].numfactors].stg_id = 2;
                        tdata[j].factors[tdata[j].numfactors].sigma = tdata[j].sigma[i];
                        tdata[j].factors[tdata[j].numfactors].vec_id = i;
                        tdata[j].factors[tdata[j].numfactors].curve_id =
                            (curve * threads) + (j * VECLEN) + i;
                        mpz_init(tdata[j].factors[tdata[j].numfactors].factor);
                        mpz_set(tdata[j].factors[tdata[j].numfactors].factor, tmp_factor);

                        tdata[j].numfactors++;
                        //gmp_printf("thread %d now has %d factors in stg2: %Zd\n", 
                        //    j, tdata[j].numfactors, tmp_factor);
                    }
                }
            }
        }

        if (verbose >= 0)
        {
            printf("ecm: %d/%d curves on C%d @ B1=%lu, B2=100*B1\r",
                (curve + VECLEN) * threads, tdata[0].curves * threads,
                (int)gmp_base10(gmpn), STAGE1_MAX);  fflush(stdout);
        }

        if (found)
            break;

    }

    if (tdata[0].save_b1)
    {
        fclose(save);
    }

    tdata[0].curves = curve * threads;

    gettimeofday(&fullstopt, NULL);
    t_time = ytools_difftime(&fullstartt, &fullstopt);
    if (verbose >= 0)
    {
        printf("\necm: process took %1.4f seconds.\n", t_time);
    }

    vecClear(one);
    mpz_clear(gmpn);
    mpz_clear(gmpt);
    mpz_clear(tmp_factor);
    return;
}

//#define PRINT_DEBUG

void vec_build_one_curve(thread_data_t *tdata, mpz_t X, mpz_t Z, mpz_t A, uint64_t sigma)
{
    vec_monty_t *mdata = tdata->mdata;
    ecm_work *work = tdata->work;
    uint32_t tid = tdata->tid;
    ecm_pt *P = &work->pt1;

    mpz_t n, u, v, t1, t2, t3, t4;
    mpz_init(n);
    mpz_init(u);
    mpz_init(v);
    mpz_init(t1);
    mpz_init(t2);
    mpz_init(t3);
    mpz_init(t4);

    extract_bignum_from_vec_to_mpz(n, mdata->n, 0, mdata->NWORDS);

    if (sigma == 0)
    {
        do
        {
            work->sigma = spRandp(&tdata->lcg_state, 6, VEC_MAXDIGIT);
        } while (work->sigma < 6);
    }
    else
    {
        work->sigma = sigma;
    }


#ifdef PRINT_DEBUG
    printf("sigma = %lu\n", work->sigma);
#endif

    // v = 4*sigma
    mpz_set_ui(v, work->sigma);
    mpz_mul_2exp(v, v, 2);

#ifdef PRINT_DEBUG
    gmp_printf("v = %Zx\n", v);
#endif

    // u = sigma^2 - 5
    mpz_set_ui(u, work->sigma);
    mpz_mul(u, u, u);
    mpz_sub_ui(u, u, 5);

#ifdef PRINT_DEBUG
    gmp_printf("u = sigma^2 - 5 = %Zx\n", u);
#endif

    // x = u^3
    mpz_mul(X, u, u);
    mpz_mul(X, X, u);
    mpz_tdiv_r(X, X, n);

#ifdef PRINT_DEBUG
    gmp_printf("u^3 = %Zx\n", X);
#endif

    // z = v^3
    mpz_mul(Z, v, v);
    mpz_mul(Z, Z, v);
    mpz_tdiv_r(Z, Z, n);

#ifdef PRINT_DEBUG
    gmp_printf("v^3 = %Zx\n", Z);
#endif


    // compute parameter A
    // (v-u)
    if (mpz_cmp(u, v) > 0)
    {
        mpz_sub(t1, v, u);
        mpz_add(t1, t1, n);
    }
    else
    {
        mpz_sub(t1, v, u);
    }

    // (v-u)^3
    mpz_mul(t2, t1, t1);
    mpz_tdiv_r(t2, t2, n);
    mpz_mul(t4, t2, t1);
    mpz_tdiv_r(t4, t4, n);

    // 3u + v
    mpz_mul_ui(t1, u, 3);
    mpz_add(t3, t1, v);
    mpz_tdiv_r(t3, t3, n);

    // a = (v-u)^3 * (3u + v)
    mpz_mul(t1, t3, t4);
    mpz_tdiv_r(t1, t1, n);

#ifdef PRINT_DEBUG
    gmp_printf("(v-u)^3 = %Zx\n", t4);
    gmp_printf("(3u + v) = %Zx\n", t3);
    gmp_printf("a = %Zx\n", t1);
#endif

    if (0)
    {
        // This is how gmp-ecm does things since sometimes they want
        // the Montgomery parameter A.  We always use Suyama's parameterization
        // so we just go ahead and build (A+2)/4.
        // 4*x*v
        mpz_mul_ui(t2, X, 4);
        mpz_mul(t4, t2, v);
        mpz_tdiv_r(t4, t4, n);

        // 4*x*v * z
        mpz_mul(t3, t4, Z);
        mpz_tdiv_r(t3, t3, n);

        /* u = 1 / (v^3 * 4*u^3*v) */
        mpz_invert(t2, t3, n);

        /* v = z^(-1) (mod n)  = 1 / v^3   */
        mpz_mul(v, t2, t4);   
        mpz_tdiv_r(v, v, n);

        /* x = x * z^(-1)      = u^3 / v^3 */
        mpz_mul(X, X, v);
        mpz_tdiv_r(X, X, n);

        /* v = b^(-1) (mod n)  = 1 / 4*u^3*v */
        mpz_mul(v, t2, Z);
        mpz_tdiv_r(v, v, n);

        /* t = ((v-u)^3 * (3*u+v)) / 4*u^3*v */
        mpz_mul(t1, t1, v);       
        mpz_tdiv_r(t1, t1, n);

        /* A = ((v-u)^3 * (3*u+v)) / 4*u^3*v - 2*/
        mpz_sub_ui(A, t1, 2);
        if (mpz_sgn(A) < 0)
        {
            mpz_add(A, A, n);
        }

        mpz_mul_2exp(X, X, DIGITBITS * mdata->NWORDS);
        mpz_tdiv_r(X, X, n);
        mpz_mul_2exp(Z, Z, DIGITBITS * mdata->NWORDS);
        mpz_tdiv_r(Z, Z, n);
        mpz_mul_2exp(A, A, DIGITBITS * mdata->NWORDS);
        mpz_tdiv_r(A, A, n);

        //gmp_printf("X/Z = %Zx\n", X);
        //gmp_printf("Z = %Zx\n", Z);
        //gmp_printf("A = %Zx\n", A);
        //
        //printf("(A+2)*B/4\n");

        mpz_set_ui(t1, 2);
        mpz_mul_2exp(t1, t1, DIGITBITS * mdata->NWORDS);
        mpz_add(A, A, t1);
        mpz_tdiv_r(A, A, n);

        if (mpz_odd_p(A))
        {
            mpz_add(A, A, n);
        }
        mpz_tdiv_q_2exp(A, A, 1);

        if (mpz_odd_p(A))
        {
            mpz_add(A, A, n);
        }
        mpz_tdiv_q_2exp(A, A, 1);

        //gmp_printf("A = %Zx\n", A);
        //exit(1);
    }
    else
    {
        // We always use Suyama's parameterization
        // so we just go ahead and build (A+2)/4.

        // 16*u^3*v
        mpz_mul_ui(t2, X, 16);
        mpz_mul(t4, t2, v);
        mpz_tdiv_r(t4, t4, n);

#ifdef PRINT_DEBUG
        gmp_printf("16*u^3*v = %Zx\n", t4);
#endif

        // accomplish the division by multiplying by the modular inverse
        // of the denom
        mpz_invert(t2, t4, n);

#ifdef PRINT_DEBUG
        gmp_printf("inv = %Zx\n", t2);
#endif

        // t1 = b = (v - u)^3 * (3*u + v) / 16u^3v)
        mpz_mul(A, t1, t2);
        mpz_tdiv_r(A, A, n);

#ifdef PRINT_DEBUG
        gmp_printf("b = %Zx\n", A);
#endif

        mpz_invert(t1, Z, n);
        mpz_mul(X, X, t1);
        mpz_set_ui(Z, 1);

        if (!mdata->isMersenne)
        {
            // into Monty rep
            mpz_mul_2exp(X, X, DIGITBITS * mdata->NWORDS);
            mpz_tdiv_r(X, X, n);
            mpz_mul_2exp(Z, Z, DIGITBITS * mdata->NWORDS);
            mpz_tdiv_r(Z, Z, n);
            mpz_mul_2exp(A, A, DIGITBITS * mdata->NWORDS);
            mpz_tdiv_r(A, A, n);
        }
        else
        {
            // Mersenne's use the special properties of the 
            // representation for reduction, no Monty needed.
            mpz_tdiv_r(X, X, n);
            mpz_tdiv_r(Z, Z, n);
            mpz_tdiv_r(A, A, n);
        }
        

        //gmp_printf("X = %Zx\n", X);
        //gmp_printf("Z = %Zx\n", Z);
        //gmp_printf("A = %Zx\n", A);
        //exit(1);
    }
	
    mpz_clear(n);
    mpz_clear(u);
    mpz_clear(v);
    mpz_clear(t1);
    mpz_clear(t2);
    mpz_clear(t3);
    mpz_clear(t4);

	return;
}

//#define PRINT_DEBUG

void build_one_curve_param1(thread_data_t *tdata, mpz_t X, mpz_t Z, 
    mpz_t X2, mpz_t Z2, mpz_t A, uint64_t sigma)
{
    vec_monty_t *mdata = tdata->mdata;
    ecm_work *work = tdata->work;
    uint32_t tid = tdata->tid;
    ecm_pt *P = &work->pt1;
    int i;
    mpz_t n, v, d, dorig;
    uint32_t NWORDS = tdata[0].NWORDS;
    uint32_t NBLOCKS = tdata[0].NBLOCKS;
    uint32_t MAXBITS = tdata[0].MAXBITS;

    mpz_init(n);
    mpz_init(v);
    mpz_init(d);
    mpz_init(dorig);

    extract_bignum_from_vec_to_mpz(n, mdata->n, 0, NWORDS);

    if (sigma == 0)
    {
        do
        {
            work->sigma = spRandp(&tdata->lcg_state, 6, VEC_MAXDIGIT) & 0X3FFFFFF;  
        } while (work->sigma < 6);
    }
    else
    {
        work->sigma = sigma;
    }
    
#ifdef PRINT_DEBUG
    printf("sigma = %lu\n", work->sigma);
#endif

    // v = sigma^2
    mpz_set_ui(v, work->sigma);
    mpz_mul(v, v, v);

    /* A=4*d-2 with d = sigma^2/2^GMP_NUMB_BITS*/
      /* Compute d = sigma^2/2^GMP_NUMB_BITS */
    for (i = 0; i < 52; i++)
    {
        if (mpz_tstbit(v, 0) == 1)
            mpz_add(v, v, n);
        mpz_div_2exp(v, v, 1);
    }

    mpz_mod(v, v, n);

#ifdef PRINT_DEBUG
    gmp_printf("n = %Zd\n", n);
#endif
#ifdef PRINT_DEBUG
    gmp_printf("d = %Zd\n", v);
#endif
    mpz_set(dorig, v);

    /* TODO add d!=-1/8*/
    if (mpz_sgn(v) == 0 || mpz_cmp_ui(v, 1) == 0)
    {
        printf("parameter cannot be 0 or 1\n");
        exit(1);
    }

    mpz_mul_2exp(v, v, 2);           /* 4d */
    mpz_sub_ui(v, v, 2);             /* 4d-2 */

    // copy and convert to Montgomery representation
#ifdef PRINT_DEBUG
    gmp_printf("4d-2 = %Zd\n", v);
#endif
    //mpres_set_z(A, v, n);
    mpz_mul_2exp(A, v, MAXBITS);
    mpz_mod(A, A, n);
#ifdef PRINT_DEBUG
    gmp_printf("A = %Zd\n", A);
#endif

    /* Compute d=(A+2)/4 from A and d'=B*d thus d' = 2^(GMP_NUMB_BITS-2)*(A+2) */
    mpres_get_z(d, A, n, mdata, MAXBITS);
#ifdef PRINT_DEBUG
    gmp_printf("d (mpres_get_z) = %Zd\n", d);
#endif
    mpz_add_ui(d, d, 2);
#ifdef PRINT_DEBUG
    gmp_printf("d (mpz_add_ui) = %Zd\n", d);
#endif
    mpz_mul_2exp(d, d, 52 - 2);
#ifdef PRINT_DEBUG
    gmp_printf("d (mpz_mul_2exp) = %Zd\n", d);
#endif
    mpz_mod(d, d, n);
#ifdef PRINT_DEBUG
    gmp_printf("d (mpz_mod) = %Zd\n", d);
#endif

    mpz_set(A, d);

#ifdef PRINT_DEBUG
    gmp_printf("d' = (A+2)*2^(GMP_NUMB_BITS-2) = %Zd\n", d);
#endif
#ifdef PRINT_DEBUG
    gmp_printf("dorig = %Zd\n", dorig);
#endif

    //mpres_set_ui(X, 2, n);
    mpz_set_ui(X, 2);
    mpz_mul_2exp(X, X, MAXBITS);
    mpz_mod(X, X, n);

#ifdef PRINT_DEBUG
    gmp_printf("X = %Zd\n", X);
#endif

    //mpres_set_ui(Z, 1, n);
    mpz_set_ui(Z, 1);
    mpz_mul_2exp(Z, Z, MAXBITS);
    mpz_tdiv_r(Z, Z, n);

#ifdef PRINT_DEBUG
    gmp_printf("Z = %Zd\n", Z);
#endif

    /* Compute 2P : no need to duplicate P, the coordinates are simple. */
    mpz_set_ui(X2, 9);
    mpz_mul_2exp(X2, X2, MAXBITS);
    mpz_tdiv_r(X2, X2, n);

#ifdef PRINT_DEBUG
    gmp_printf("X2 = %Zd\n", X2);
#endif

    mpz_set(Z2, d);
    mpz_mul_2exp(Z2, Z2, MAXBITS);
    mpz_tdiv_r(Z2, Z2, n);
    mpres_div_2exp(Z2, Z2, GMP_NUMB_BITS, n);

    mpres_mul_2exp(Z2, Z2, 6, n);
    mpres_add_ui(Z2, Z2, 8, n, MAXBITS); /* P2 <- 2P = (9 : : 64d+8) */

#ifdef PRINT_DEBUG
    gmp_printf("Z2 = %Zd\n", Z2);
#endif

    mpz_clear(n);
    mpz_clear(v);
    mpz_clear(d);

    return;
}

//#define TESTMUL
void vec_ecm_stage1(vec_monty_t *mdata, ecm_work *work, ecm_pt *P, 
    uint64_t b1, uint64_t*primes, uint64_t nump, int verbose)
{
	int i;
	uint64_t q;

#ifdef PARAM1
    {
        mpz_t s, f;
        // timing variables
        struct timeval stopt;	// stop time of this job
        struct timeval startt;	// start time of this job
        double t_time;

        gettimeofday(&startt, NULL);
        mpz_init(s);
        mpz_init(f);
        work->last_pid = compute_s(s, PRIMES, stg1);
        gettimeofday(&stopt, NULL);
        t_time = ytools_difftime(&startt, &stopt);
        printf("Built product of %lu bits in %1.0f ms\n", mpz_sizeinbase(s, 2), t_time * 1000);
        ecm_stage1_batch(f, work, P, work->s, work->n, stg1, s, mdata);
        mpz_clear(s);
        mpz_clear(f);
        return;
    }
#endif

#ifdef MPZ_VERIFY
    mpz_t gtmp, gtmp2, gn, gnhat;
    mpz_t gs, gsum, gdiff;

    extract_bignum_from_vec_to_mpz(mdata->gmp_t1, P->X, 0, mdata->NWORDS);
    extract_bignum_from_vec_to_mpz(mdata->gmp_t2, P->Z, 0, mdata->NWORDS);

    mpz_init(gn);
    mpz_init(gnhat);
    mpz_init(gtmp);
    mpz_init(gtmp2);
    mpz_init(gs);
    mpz_init(gsum);
    mpz_init(gdiff);

    extract_bignum_from_vec_to_mpz(gn, work->n, 0, mdata->NWORDS);
    extract_bignum_from_vec_to_mpz(gs, work->s, 0, mdata->NWORDS);
    extract_bignum_from_vec_to_mpz(gnhat, mdata->vnhat, 0, mdata->NWORDS);

#endif

	// handle the only even case 
    q = 2;
	while (q < b1)
	{
		vecsubmod_ptr(P->X, P->Z, work->diff1, mdata);
		vecaddmod_ptr(P->X, P->Z, work->sum1, mdata);
		vec_duplicate(mdata, work, work->sum1, work->diff1, P);

#ifdef MPZ_VERIFY
        mpz_elldup(mdata->gmp_t1, mdata->gmp_t2, gn, gnhat,
            gsum, gdiff, gs, gtmp, gtmp2, mdata->NWORDS * 52);
#endif
		q *= 2;
	}

	for (i = 1; (i < nump) && ((uint32_t)primes[i] < b1); i++)
	{
		uint64_t c = 1;
	
		q = primes[i];

		do {
			vec_prac(mdata, work, P, q);
			c *= q;
		} while ((c * q) < b1);
	
#ifdef SKYLAKEX
		if ((verbose > 1) && ((i & 8191) == 0))
#else
        if ((verbose > 1) && ((i & 511) == 0))
#endif
		{
			printf("accumulating prime %lu\r", q);
			fflush(stdout);
		}
	}

	work->last_pid = i;

    if (verbose > 1)
	{
		printf("\nStage 1 completed at prime %lu with %u point-adds and %u point-doubles\n", 
			primes[i-1], work->stg1Add, work->stg1Doub);
		fflush(stdout);
	}

#ifdef MPZ_VERIFY
    mpz_clear(gn);
    mpz_clear(gnhat);
    mpz_clear(gtmp);
    mpz_clear(gtmp2);
    mpz_clear(gs);
    mpz_clear(gsum);
    mpz_clear(gdiff);
#endif

	return;
}

#if 0
void vec_ecm_stage2_init_old(ecm_pt *P, vec_monty_t *mdata, ecm_work *work, base_t *primes, int verbose)
{
	// run Montgomery's PAIR algorithm.  
	uint32_t D = work->D;
	uint32_t R = work->R;
	uint32_t w = D;
	uint32_t U = work->U;
	uint32_t L = work->L;
	uint32_t umax = U * w;
	int i, j, k, pid;
	
	uint32_t ainc = 2 * D;
	uint32_t ascale = 2 * D;
	uint32_t amin = work->amin = (STAGE1_MAX + w) / ascale;
	uint32_t s;
	uint32_t a;
	
	uint32_t numR;
	uint32_t u, ap;
	int q, mq;
	int debug = 0;

	uint32_t *rprime_map_U = work->map;
	ecm_pt *Pa = work->Pa;
	vec_bignum_t **Paprod = work->Paprod;
	vec_bignum_t **Pbprod = work->Pbprod;
	ecm_pt *Pb = work->Pb;
	ecm_pt *Pd;
	vec_bignum_t *acc = work->stg2acc;


    if (verbose > 1)
		printf("\n");

	work->paired = 0;
	work->numprimes = 0;
	work->ptadds = 0;
	work->stg1Add = 0;
	work->stg1Doub = 0;

	//stage 2 init
	//Q = P = result of stage 1
	//compute [d]Q for 0 < d <= D
	Pd = &Pb[rprime_map_U[w]];

	// [1]Q
	vecCopy(P->Z, Pb[1].Z);
	vecCopy(P->X, Pb[1].X);
	vecmulmod_ptr(Pb[1].X, Pb[1].Z, Pbprod[1], work->n, work->tt4, mdata);

	// [2]Q
	vecCopy(P->Z, Pb[2].Z);
	vecCopy(P->X, Pb[2].X);
    vecaddsubmod_ptr(P->X, P->Z, work->sum1, work->diff1, mdata);
	//vecsubmod_ptr(P->X, P->Z, work->diff1, mdata);
	//vecaddmod_ptr(P->X, P->Z, work->sum1, mdata);
	vec_duplicate(mdata, work, work->sum1, work->diff1, &Pb[2]);
	vecmulmod_ptr(Pb[2].X, Pb[2].Z, Pbprod[2], work->n, work->tt4, mdata);

	// [3]Q, [4]Q, ... [D]Q
	// because of the way we pick 'a' and 'D', we'll only need the 479 points that 
	// are relatively prime to D.  We need to keep a few points during the running 
	// computation i.e., ... j, j-1, j-2. The lookup table rprime_map maps the 
	// integer j into the relatively prime indices that we store. We also store
	// points 1, 2, and D, and storage point 0 is used as scratch for a total of
	// 483 stored points.

	vecCopy(Pb[1].X, work->pt2.X);
	vecCopy(Pb[1].Z, work->pt2.Z);
	vecCopy(Pb[2].X, work->pt1.X);
	vecCopy(Pb[2].Z, work->pt1.Z);

	for (j = 3; j <= U * D; j++)
	{
		ecm_pt *P1 = &work->pt1;			// Sd - 1
		ecm_pt *P2 = &Pb[1];				// S1
		ecm_pt *P3 = &work->pt2;			// Sd - 2
		ecm_pt *Pout = &Pb[rprime_map_U[j]];	// Sd

		// vecAdd:
		//x+ = z- * [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
		//z+ = x- * [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
		//x- = original x
		//z- = original z

		// compute Sd from Sd-1 + S1, requiring Sd-1 - S1 = Sd-2
        vecaddsubmod_ptr(P1->X, P1->Z, work->sum1, work->diff1, mdata);
        //vecsubmod_ptr(P1->X, P1->Z, work->diff1, mdata);
		//vecaddmod_ptr(P1->X, P1->Z, work->sum1, mdata);
        vecaddsubmod_ptr(P2->X, P2->Z, work->sum2, work->diff2, mdata);
		//vecsubmod_ptr(P2->X, P2->Z, work->diff2, mdata);
		//vecaddmod_ptr(P2->X, P2->Z, work->sum2, mdata);

		vecmulmod_ptr(work->diff1, work->sum2, work->tt1, work->n, work->tt4, mdata);	//U
		vecmulmod_ptr(work->sum1, work->diff2, work->tt2, work->n, work->tt4, mdata);	//V

        vecaddsubmod_ptr(work->tt1, work->tt2, Pout->X, Pout->Z, mdata);
		//vecaddmod_ptr(work->tt1, work->tt2, Pout->X, mdata);		//U + V
		//vecsubmod_ptr(work->tt1, work->tt2, Pout->Z, mdata);		//U - V
		vecsqrmod_ptr(Pout->X, work->tt1, work->n, work->tt4, mdata);					//(U + V)^2
		vecsqrmod_ptr(Pout->Z, work->tt2, work->n, work->tt4, mdata);					//(U - V)^2

		// if gcd(j,D) != 1, Pout maps to scratch space (Pb[0])
		vecmulmod_ptr(work->tt1, P3->Z, Pout->X, work->n, work->tt4, mdata);			//Z * (U + V)^2
		vecmulmod_ptr(work->tt2, P3->X, Pout->Z, work->n, work->tt4, mdata);			//x * (U - V)^2

		//store Pb[j].X * Pb[j].Z as well
		vecmulmod_ptr(Pout->X, Pout->Z, Pbprod[rprime_map_U[j]],
			work->n, work->tt4, mdata);

		work->ptadds++;

		// advance
		vecCopy(P1->X, P3->X);
		vecCopy(P1->Z, P3->Z);
		vecCopy(Pout->X, P1->X);
		vecCopy(Pout->Z, P1->Z);
	}

	//printf("B table generated to umax = %d\n", U * D);

	// Pd = [2w]Q
	vecCopy(P->Z, Pd->Z);
	vecCopy(P->X, Pd->X);
	next_pt_vec(mdata, work, Pd, ainc);
	//vec_prac(mdata, work, Pd, ainc);

	//first a value: first multiple of D greater than B1
	work->A = amin * ascale;

	//initialize info needed for giant step
	vecCopy(P->Z, Pa[0].Z);
	vecCopy(P->X, Pa[0].X);
	next_pt_vec(mdata, work, &Pa[0], work->A);
	//vec_prac(mdata, work, &Pa[0], work->A);

	//and Paprod
	vecmulmod_ptr(Pa[0].X, Pa[0].Z, Paprod[0], work->n, work->tt4, mdata);
	if (verbose & (debug == 2))
		printf("Pa[0] = [%lu]Q\n", work->A);

	vecCopy(P->Z, work->Pad->Z);
	vecCopy(P->X, work->Pad->X);
	next_pt_vec(mdata, work, work->Pad, work->A - ainc);
	//vec_prac(mdata, work, work->Pad, work->A - ainc);

	if (verbose & (debug == 2))
		printf("Pad = [%lu]Q\n", work->A - ainc);

	vecaddmod_ptr(Pa[0].X, Pa[0].Z, work->sum1, mdata);
	vecaddmod_ptr(Pd->X, Pd->Z, work->sum2, mdata);
	vecsubmod_ptr(Pa[0].X, Pa[0].Z, work->diff1, mdata);
	vecsubmod_ptr(Pd->X, Pd->Z, work->diff2, mdata);
	vec_add(mdata, work, work->Pad, &Pa[1]);
	vecmulmod_ptr(Pa[1].X, Pa[1].Z, Paprod[1], work->n, work->tt4, mdata);

	work->A += ainc;
	if (verbose & (debug == 2))
		printf("Pa[1] = [%lu]Q\n", work->A + ainc);

	for (i = 2; i < 2 * L; i++)
	{
		//giant step - use the addition formula for ECM
		//Pa + Pd
		//x+ = z- * [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
		//z+ = x- * [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
		//x- = [a-d]x
		//z- = [a-d]z
		//vecaddmod_ptr(Pa[i - 1].X, Pa[i - 1].Z, work->sum1, mdata);
		//vecaddmod_ptr(Pd->X, Pd->Z, work->sum2, mdata);
		//vecsubmod_ptr(Pa[i - 1].X, Pa[i - 1].Z, work->diff1, mdata);
		//vecsubmod_ptr(Pd->X, Pd->Z, work->diff2, mdata);
        vecaddsubmod_ptr(Pa[i - 1].X, Pa[i - 1].Z, work->sum1, work->diff1, mdata);
        vecaddsubmod_ptr(Pd->X, Pd->Z, work->sum2, work->diff2, mdata);
		vec_add(mdata, work, &Pa[i - 2], &Pa[i]);

		work->A += ainc;
		work->ptadds++;
		if (verbose & (debug == 2))
			printf("Pa[%d] = [%lu]Q\n", i, work->A + i * ainc);

		//and Paprod
		vecmulmod_ptr(Pa[i].X, Pa[i].Z, Paprod[i], work->n, work->tt4, mdata);
	}

	if (verbose & (debug == 2))
		printf("A table generated to 2 * L = %d\n", 2 * L);

	// initialize accumulator
    vecCopy(mdata->one, acc);

	return;
}
#endif

// batch inversions:
// From: H. Shacham and D. Boneh, "Improving SSL Handshake Performance via Batching"
// A general batched - inversion algorithm proceeds, in three phases, as follows.
// first, set A1 = x1 and Ai = xi * A(i-1) so that Ai = prod(j=1,i,xj).
// then, invert An and store in Bn and set Bi = x(i+1) * B(i+1) for i < n.
// Now we have Bi = prod(j=1,i,xj^-1).
// finally, set C1 = B1 and Ci = A(i-1) * B(i) for i > 1.
// Then Ci = xi^-1 for i > 1.
// each phase takes n-1 multiplications so we have 3n-3 total multiplications
// and one inversion mod N.

// Now, how do we use this to speed up stage 2?
// In projective coordinates (Px:Pz) = (Qx:Qz) does not
// imply that Px=Qx.  We need to cancel the Z-coordinates.
// to cancel a z-coord we can multiply a point by z^-1, then
// we have (Px/Pz:1).
// the cross product we compute and accumulate is:
// (Ax - Bx) * (Az + Bz) + BxBz - AxAz  == AxBz - BxAz
// In the baby-step giant-step algorithm for stage 2 we can normalize
// all pre-computed Pb's by the batch-inversion process.
// Likewise if we precompute all Pa's we can normalize in the same way.
// (If there are too many Pa's then maybe in batches at a time?)
// Then we simply work with the normalized values and directly 
// accumulate the cross products (Xr/Zr*1 - Xd/Zd*1)

#define CROSS_PRODUCT_INV \
    vecsubmod_ptr(work->Pa_inv[pa], Pb[rprime_map_U[pb]].X, work->tt1, mdata);          \
    vecmulmod_ptr(acc, work->tt1, acc, work->n, work->tt4, mdata);        

#define CROSS_PRODUCT \
    vecsubmod_ptr(Pa[pa].X, Pb[rprime_map_U[pb]].X, work->tt1, mdata);          \
    vecaddmod_ptr(Pa[pa].Z, Pb[rprime_map_U[pb]].Z, work->tt2, mdata);          \
    vecmulmod_ptr(work->tt1, work->tt2, work->tt3, work->n, work->tt4, mdata);    \
    vecaddmod_ptr(work->tt3, Pbprod[rprime_map_U[pb]], work->tt1, mdata);       \
    vecsubmod_ptr(work->tt1, Paprod[pa], work->tt2, mdata);                     \
    vecmulmod_ptr(acc, work->tt2, acc, work->n, work->tt4, mdata);      

int batch_invert_pt_inplace(ecm_pt* pts_to_Zinvert,
    vec_bignum_t** tmp_vec, vec_monty_t* mdata, ecm_work* work, int num)
{
    vec_bignum_t** B;
    vec_bignum_t** A = tmp_vec;
    int i;
    int j;
    int inverr;
    int foundDuringInv = 0;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t MAXBITS = mdata->MAXBITS;
    mpz_t gmptmp, gmptmp2, gmpn;

    mpz_init(gmptmp);
    mpz_init(gmptmp2);
    mpz_init(gmpn);

    work->numinv++;

    // here, we have temporary space for B, A is put into the unused Pbprod, and C is Pb.Z.
    // faster batch inversion in three phases, as follows:
    // first, set A1 = z1 and Ai = zi * A(i-1) so that Ai = prod(j=1,i,zj).
    vecCopy(pts_to_Zinvert[1].Z, A[1]);
    for (i = 2; i < num; i++)
    {
        vecmulmod_ptr(pts_to_Zinvert[i].Z, A[i - 1], A[i], work->n, work->tt4, mdata);
    }

    B = (vec_bignum_t **)malloc(num * sizeof(vec_bignum_t*));

    for (j = 0; j < num; j++)
    {
        B[j] = vecInit(NWORDS);
    }

    // now we have to take An out of monty rep so we can invert it.
    if (mdata->isMersenne == 0)
    {
        vecClear(work->tt1);
        for (j = 0; j < VECLEN; j++)
        {
            work->tt1->data[j] = 1;
        }
        work->tt1->size = 1;
        vecmulmod_ptr(A[num - 1], work->tt1, B[num - 1], work->n, work->tt4, mdata);
    }
    else
    {
        vecCopy(A[num - 1], B[num - 1]);
    }

    extract_bignum_from_vec_to_mpz(gmpn, mdata->n, 0, NWORDS);
    for (j = 0; j < VECLEN; j++)
    {
        // extract this vec position so we can use mpz_invert.
        extract_bignum_from_vec_to_mpz(gmptmp, B[num - 1], j, NWORDS);

        // invert it
        inverr = mpz_invert(gmptmp2, gmptmp, gmpn);

        if (inverr == 0)
        {
            //extract_bignum_from_vec_to_mpz(gmptmp, work->tt2, j, NWORDS);
            //printf("inversion error\n");
            //gmp_printf("tried to invert %Zd mod %Zd in stage2init Pb\n", gmptmp, gmpn);
            mpz_gcd(gmptmp, gmptmp, gmpn);
            //gmp_printf("the GCD is %Zd\n", gmptmp);
            int k;
            for (k = 0; k < NWORDS; k++)
                work->stg2acc->data[k * VECLEN + j] = 0;
            insert_mpz_to_vec(work->stg2acc, gmptmp, j);
            foundDuringInv = 1;
        }

        if (mdata->isMersenne == 0)
        {
            // now put it back into Monty rep.
            mpz_mul_2exp(gmptmp2, gmptmp2, MAXBITS);
            mpz_tdiv_r(gmptmp2, gmptmp2, gmpn);
        }

        // and stuff it back in the vector.
        insert_mpz_to_vec(B[num - 1], gmptmp2, j);
    }

    //if (doneIfFoundDuringInv && foundDuringInv)
    //{
    //    work->last_pid = -1;
    //
    //    mpz_clear(gmptmp);
    //    mpz_clear(gmptmp2);
    //    mpz_clear(gmpn);
    //
    //    return;
    //}


    // and continue.
    for (i = num - 2; i >= 0; i--)
    {
        vecmulmod_ptr(pts_to_Zinvert[i + 1].Z, B[i + 1], B[i], work->n, work->tt4, mdata);
    }

    // Now we have Bi = prod(j=1,i,zj^-1).
    // finally, set C1 = B1 and Ci = A(i-1) * B(i) for i > 1.
    // Then Ci = zi^-1 for i > 1.
    vecCopy(B[1], pts_to_Zinvert[1].Z);

    for (i = 2; i < num; i++)
    {
        vecmulmod_ptr(B[i], A[i - 1], pts_to_Zinvert[i].Z, work->n, work->tt4, mdata);
    }

    // each phase takes n-1 multiplications so we have 3n-3 total multiplications
    // and one inversion mod N.
    // but we still have to combine with the X coord.
    for (i = 1; i < num; i++)
    {
        vecmulmod_ptr(pts_to_Zinvert[i].X, pts_to_Zinvert[i].Z, pts_to_Zinvert[i].X,
            work->n, work->tt4, mdata);
    }

    for (j = 0; j < num; j++)
    {
        vecFree(B[j]);
    }
    free(B);

    mpz_clear(gmptmp);
    mpz_clear(gmptmp2);
    mpz_clear(gmpn);

    return foundDuringInv;

}

int batch_invert_pt_to_bignum(ecm_pt* pts_to_Zinvert, vec_bignum_t** out,
    vec_bignum_t** tmp_vec, vec_monty_t* mdata, ecm_work* work, int startid, int stopid)
{
    vec_bignum_t** B;
    vec_bignum_t** A = tmp_vec;
    int i;
    int j;
    int inverr = 0;
    int foundDuringInv = 0;
    int num = stopid - startid;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t MAXBITS = mdata->MAXBITS;
    mpz_t gmptmp, gmptmp2, gmpn;

    mpz_init(gmptmp);
    mpz_init(gmptmp2);
    mpz_init(gmpn);

    work->numinv++;

    // here, we have temporary space for B, A is put into the unused Pbprod, and C is Pb.Z.
    // faster batch inversion in three phases, as follows:
    // first, set A1 = z1 and Ai = zi * A(i-1) so that Ai = prod(j=1,i,zj).
    vecCopy(pts_to_Zinvert[startid].Z, A[0]);
    for (i = 1; i < num; i++)
    {
        vecmulmod_ptr(pts_to_Zinvert[startid + i].Z, A[i - 1], A[i], work->n, work->tt4, mdata);
    }

    B = (vec_bignum_t **)malloc(num * sizeof(vec_bignum_t*));

    for (j = 0; j < num; j++)
    {
        B[j] = vecInit(NWORDS);
    }

    // now we have to take An out of monty rep so we can invert it.
    if (mdata->isMersenne == 0)
    {
        vecClear(work->tt1);
        for (j = 0; j < VECLEN; j++)
        {
            work->tt1->data[j] = 1;
        }
        work->tt1->size = 1;
        vecmulmod_ptr(A[num - 1], work->tt1, B[num - 1], work->n, work->tt4, mdata);
    }
    else
    {
        vecCopy(A[num - 1], B[num - 1]);
    }

    extract_bignum_from_vec_to_mpz(gmpn, mdata->n, 0, NWORDS);
    for (j = 0; j < VECLEN; j++)
    {
        // extract this vec position so we can use mpz_invert.
        extract_bignum_from_vec_to_mpz(gmptmp, B[num - 1], j, NWORDS);

        // invert it
        inverr = mpz_invert(gmptmp2, gmptmp, gmpn);

        if (inverr == 0)
        {
            //extract_bignum_from_vec_to_mpz(gmptmp, work->tt2, j, NWORDS);
            //printf("inversion error\n");
            //gmp_printf("tried to invert %Zd mod %Zd in stage2init Pb\n", gmptmp, gmpn);
            mpz_gcd(gmptmp, gmptmp, gmpn);
            //gmp_printf("the GCD is %Zd\n", gmptmp);
            int k;
            for (k = 0; k < NWORDS; k++)
                work->stg2acc->data[k * VECLEN + j] = 0;
            insert_mpz_to_vec(work->stg2acc, gmptmp, j);
            foundDuringInv = 1;
        }

        if (mdata->isMersenne == 0)
        {
            // now put it back into Monty rep.
            mpz_mul_2exp(gmptmp2, gmptmp2, MAXBITS);
            mpz_tdiv_r(gmptmp2, gmptmp2, gmpn);
        }

        // and stuff it back in the vector.
        insert_mpz_to_vec(B[num - 1], gmptmp2, j);
    }

    //if (doneIfFoundDuringInv && foundDuringInv)
    //{
    //    work->last_pid = -1;
    //
    //    mpz_clear(gmptmp);
    //    mpz_clear(gmptmp2);
    //    mpz_clear(gmpn);
    //
    //    return;
    //}


    // and continue.
    for (i = num - 2; i >= 0; i--)
    {
        vecmulmod_ptr(pts_to_Zinvert[startid + i + 1].Z, B[i + 1], B[i], work->n, work->tt4, mdata);
    }

    // Now we have Bi = prod(j=1,i,zj^-1).
    // finally, set C1 = B1 and Ci = A(i-1) * B(i) for i > 1.
    // Then Ci = zi^-1 for i > 1.
    vecCopy(B[0], out[startid + 0]);

    for (i = 1; i < num; i++)
    {
        vecmulmod_ptr(B[i], A[i - 1], out[startid + i], work->n, work->tt4, mdata);
    }

    // each phase takes n-1 multiplications so we have 3n-3 total multiplications
    // and one inversion mod N.
    // but we still have to combine with the X coord.
    for (i = 0; i < num; i++)
    {
        vecmulmod_ptr(pts_to_Zinvert[startid + i].X, out[startid + i], out[startid + i],
            work->n, work->tt4, mdata);
    }

    for (j = 0; j < num; j++)
    {
        vecFree(B[j]);
    }
    free(B);

    mpz_clear(gmptmp);
    mpz_clear(gmptmp2);
    mpz_clear(gmpn);

    return foundDuringInv;

}

// test cases for generic input at B1=1e6, B2=100B1:
// n = 142946323174762557214361604817789197531833590620956958433836799929503392464892596183803921
//11919771003873180376
//827341355533811391
//6409678826612327146
//13778091190526084667
//10019108749973911965 *
//10593445070074576128
//16327347202299112611
//13768494887674349585
//17303758977955016383
//2123812563661387803
//2330438305415445111
//12942218412106273630
//5427613898610684157
//13727269399001077418
//3087408422684406072
//8338236510647016635
//18232185847183255223
//5070879816975737551
//9793972958987869750
//1683842010542383008
//16668736769625151751
//11148653366342049109
//6736437364141805734
//8860111571919296085
//15708855786729755459
//4263089024287634346
//10705409183485702771
//5104801995378138195
//9551766994217130412
//17824508581606173922
//4444245868135963544
//14755844915853888743
//4749513976499976002
//3933740986814285076
//2498288573977543008
//18051693002182940438
//421313926042840093
//1659254194582388863
//13762123388521706810
//1318769405167840394
//14979751960240161797
//4989253092822783329
//14628970911725975539
//4759771957864370849
//17870405635651283010
//472060146
//3776270672
//3954243165
//2576580518
//416265588


void addflag(uint8_t* flags, uint64_t loc, uint64_t lobound, uint64_t hibound)
{
    //printf("adding flag to %u\n", loc); fflush(stdout);
    //if (flags[loc] == 1)
    //    printf("duplicate location %u\n", loc);
    if ((loc < lobound) || (loc >= hibound))
    {
        //printf("attempted to add a flag to location %lu outside of bounds %lu:%lu\n",
        //    loc, lobound, hibound);
        return;
    }
    flags[loc] = 1;
    return;
}

int vec_ecm_stage2_init(ecm_pt* P, vec_monty_t* mdata, ecm_work* work, int verbose)
{
    // compute points used during the stage 2 pair algorithm.
    uint32_t w = work->D;
    uint32_t U = work->U;
    uint32_t L = work->L;
    int i, j;
    uint32_t amin = work->amin = (work->STAGE1_MAX + w) / (2 * w);
    int wscale = 1;

    int debug = 0;
    int inverr;
    int foundDuringInv = 0;
    int doneIfFoundDuringInv = 0;

    uint32_t* rprime_map_U = work->map;
    ecm_pt* Pa = work->Pa;
    vec_bignum_t** Paprod = work->Paprod;
    vec_bignum_t** Pbprod = work->Pbprod;
    ecm_pt* Pb = work->Pb;
    ecm_pt* Pd = work->Pdnorm;
    vec_bignum_t* acc = work->stg2acc;
    int lastMapID;

    work->paired = 0;
    work->numprimes = 0;
    work->ptadds = 0;
    work->ptdups = 0;
    work->numinv = 0;

    //stage 2 init
    //Q = P = result of stage 1
    //compute [d]Q for 0 < d <= D

    // [1]Q
    vecCopy(P->Z, Pb[1].Z);
    vecCopy(P->X, Pb[1].X);

    // [2]Q
    vecCopy(P->Z, Pb[2].Z);
    vecCopy(P->X, Pb[2].X);
    vecaddsubmod_ptr(P->X, P->Z, work->sum1, work->diff1, mdata);
    vec_duplicate(mdata, work, work->sum1, work->diff1, &Pb[2]);

    //printf("init: D = %d, ainc = %d, ascale = %d, U = %d, L = %d\n", D, ainc, ascale, U, L);

    // [3]Q, [4]Q, ... [D]Q
    // because of the way we pick 'a' and 'D', we'll only need the 479 points that 
    // are relatively prime to D.  We need to keep a few points during the running 
    // computation i.e., ... j, j-1, j-2. The lookup table rprime_map maps the 
    // integer j into the relatively prime indices that we store. We also store
    // points 1, 2, and D, and storage point 0 is used as scratch for a total of
    // 483 stored points.

    vecCopy(Pb[1].X, work->pt2.X);
    vecCopy(Pb[1].Z, work->pt2.Z);
    vecCopy(Pb[2].X, work->pt1.X);
    vecCopy(Pb[2].Z, work->pt1.Z);

    lastMapID = 0;
    for (j = 3; j <= U * w; j++)
    {
        ecm_pt* P1 = &work->pt1;			// Sd - 1
        ecm_pt* P2 = &Pb[1];				// S1
        ecm_pt* P3 = &work->pt2;			// Sd - 2
        ecm_pt* Pout = &Pb[rprime_map_U[j]];	// Sd

        if (rprime_map_U[j] > 0)
            lastMapID = rprime_map_U[j];

        // vecAdd:
        //x+ = z- * [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
        //z+ = x- * [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
        //x- = original x
        //z- = original z
        // compute Sd from Sd-1 + S1, requiring Sd-1 - S1 = Sd-2
        vecaddsubmod_ptr(P1->X, P1->Z, work->sum1, work->diff1, mdata);
        vecaddsubmod_ptr(P2->X, P2->Z, work->sum2, work->diff2, mdata);

        vecmulmod_ptr(work->diff1, work->sum2, work->tt1, work->n, work->tt4, mdata);	//U
        vecmulmod_ptr(work->sum1, work->diff2, work->tt2, work->n, work->tt4, mdata);	//V

        vecaddsubmod_ptr(work->tt1, work->tt2, Pout->X, Pout->Z, mdata);		        //U +/- V
        vecsqrmod_ptr(Pout->X, work->tt1, work->n, work->tt4, mdata);					//(U + V)^2
        vecsqrmod_ptr(Pout->Z, work->tt2, work->n, work->tt4, mdata);					//(U - V)^2

        // if gcd(j,D) != 1, Pout maps to scratch space (Pb[0])
        vecmulmod_ptr(work->tt1, P3->Z, Pout->X, work->n, work->tt4, mdata);			//Z * (U + V)^2
        vecmulmod_ptr(work->tt2, P3->X, Pout->Z, work->n, work->tt4, mdata);			//x * (U - V)^2

#ifndef DO_STAGE2_INV
        //store Pb[j].X * Pb[j].Z as well
        vecmulmod_ptr(Pout->X, Pout->Z, Pbprod[rprime_map_U[j]],
            work->n, work->tt4, mdata);
#endif

        work->ptadds++;

        // advance
        vecCopy(P1->X, P3->X);
        vecCopy(P1->Z, P3->Z);
        vecCopy(Pout->X, P1->X);
        vecCopy(Pout->Z, P1->Z);

        //sprintf(str, "Pb[%d].Z: ", rprime_map_U[j]);
        //print_vechex(Pout->Z->data, 0, NWORDS, str);
        if (verbose & (debug == 2))
            printf("rprime_map_U[%d] = %u\n", j, rprime_map_U[j]);
    }

    //printf("B table generated to umax = %d\n", U * D);

    // initialize accumulator
    vecCopy(mdata->one, acc);

#ifdef DO_STAGE2_INV
    // invert all of the Pb's
    foundDuringInv = batch_invert_pt_inplace(Pb, Pbprod, mdata, work, lastMapID + 1);

    if (doneIfFoundDuringInv && foundDuringInv)
    {
        work->last_pid = -1;
        return foundDuringInv;
    }
#endif

    // Pd = [w]Q
    vecCopy(P->Z, Pd->Z);
    vecCopy(P->X, Pd->X);
    next_pt_vec(mdata, work, Pd, wscale * w);

    if (verbose & (debug == 2))
        printf("Pd = [%u]Q\n", wscale * w);

    return foundDuringInv;
}

void vec_ecm_stage2_pair(uint32_t pairmap_steps, uint32_t* pairmap_v, uint32_t* pairmap_u,
    ecm_pt* P, vec_monty_t* mdata, ecm_work* work, int verbose)
{
    // use the output of the PAIR algorithm to perform stage 2.
    uint32_t w = work->D;
    uint32_t U = work->U;
    uint32_t L = work->L;
    uint32_t umax = U * w;
    int i, pid;
    uint32_t amin = work->amin;
    int foundDuringInv = 0;
    int doneIfFoundDuringInv = 0;
    int mapid;
    uint8_t* flags;
    int wscale = 1;
    int debug = 0;

    uint32_t* rprime_map_U = work->map;
    ecm_pt* Pa = work->Pa;      // non-inverted
    ecm_pt* Pb = work->Pb;      // inverted
    ecm_pt* Pd = work->Pdnorm;  // non-inverted Pd
    vec_bignum_t** Paprod = work->Paprod;
#ifndef DO_STAGE2_INV
    vec_bignum_t** Pbprod = work->Pbprod;
#endif
    vec_bignum_t* acc = work->stg2acc;


    if (1)
    {
        //first a value: first multiple of D greater than B1
        work->A = (uint64_t)amin * (uint64_t)w * 2;

        //initialize info needed for giant step
        vecCopy(P->Z, Pa[0].Z);
        vecCopy(P->X, Pa[0].X);

        next_pt_vec(mdata, work, &Pa[0], work->A);

        if (verbose & (debug == 2))
            printf("Pa[0] = [%lu]Q\n", work->A);

        vecCopy(P->Z, work->Pad->Z);
        vecCopy(P->X, work->Pad->X);
        next_pt_vec(mdata, work, work->Pad, work->A - wscale * w);

        if (verbose & (debug == 2))
            printf("Pad = [%lu]Q\n", work->A - wscale * w);

        vecaddmod_ptr(Pa[0].X, Pa[0].Z, work->sum1, mdata);
        vecaddmod_ptr(Pd->X, Pd->Z, work->sum2, mdata);
        vecsubmod_ptr(Pa[0].X, Pa[0].Z, work->diff1, mdata);
        vecsubmod_ptr(Pd->X, Pd->Z, work->diff2, mdata);
        vec_add(mdata, work, work->Pad, &Pa[1]);

        work->A += wscale * w;
        if (verbose & (debug == 2))
            printf("Pa[1] = [%lu]Q\n", work->A);

        for (i = 2; i < 2 * L; i++)
        {
            //giant step - use the addition formula for ECM
            //Pa + Pd
            //x+ = z- * [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
            //z+ = x- * [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
            //x- = [a-d]x
            //z- = [a-d]z
            vecaddsubmod_ptr(Pa[i - 1].X, Pa[i - 1].Z, work->sum1, work->diff1, mdata);
            vecaddsubmod_ptr(Pd->X, Pd->Z, work->sum2, work->diff2, mdata);
            vec_add(mdata, work, &Pa[i - 2], &Pa[i]);

#ifndef DO_STAGE2_INV
            vecmulmod_ptr(Pa[i].X, Pa[i].Z, work->Paprod[i], work->n, work->tt4, mdata);
#endif

            work->A += wscale * w;
            if (verbose & (debug == 2))
                printf("Pa[%d] = [%lu]Q\n", i, work->A);
        }

#ifdef DO_STAGE2_INV
        // and invert all of the Pa's into a separate vector
        foundDuringInv |= batch_invert_pt_to_bignum(Pa, work->Pa_inv, Paprod, mdata, work, 0, 2 * L);
        work->numinv++;
        if (doneIfFoundDuringInv && foundDuringInv)
        {
            work->last_pid = -1;
            return; // foundDuringInv;
        }
#endif

        if (verbose & (debug == 2))
            printf("A table generated to L = %d\n", 2 * L);
    }



    if (verbose > 1)
    {
        printf("commencing stage 2 at A=%lu\n"
            "w = %u, R = %u, L = %u, U = %d, umax = %u, amin = %u\n",
            2 * (uint64_t)amin * (uint64_t)w, w, work->R - 3, L, U, umax, amin);
    }

    for (mapid = 0; mapid < pairmap_steps; mapid++)
    {
        int pa, pb;

        if ((verbose > 1) && ((mapid & 65535) == 0))
        {
            printf("pairmap step %u of %u\r", mapid, pairmap_steps);
            fflush(stdout);
        }

        if ((pairmap_u[mapid] == 0) && (pairmap_v[mapid] == 0))
        {
            int shiftdist = 2;

            // shift out uneeded A's.
            // update amin by U * 2 * w.
            // each point-add increments by w, so shift 2 * U times;
            for (i = 0; i < 2 * L - shiftdist * U; i++)
            {
                vecCopy(Pa[i + shiftdist * U].X, Pa[i].X);
                vecCopy(Pa[i + shiftdist * U].Z, Pa[i].Z);
                vecCopy(work->Pa_inv[i + shiftdist * U], work->Pa_inv[i]);
            }

            // make new A's: need at least two previous points;
            // therefore we can't have U = 1
            for (i = 2 * L - shiftdist * U; i < 2 * L; i++)
            {
                //giant step - use the addition formula for ECM
                //Pa + Pd
                //x+ = z- * [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
                //z+ = x- * [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
                //x- = [a-d]x
                //z- = [a-d]z
                vecaddsubmod_ptr(Pa[i - 1].X, Pa[i - 1].Z, work->sum1, work->diff1, mdata);
                vecaddsubmod_ptr(Pd->X, Pd->Z, work->sum2, work->diff2, mdata);
                vec_add(mdata, work, &Pa[i - 2], &Pa[i]);

#ifndef DO_STAGE2_INV
                vecmulmod_ptr(Pa[i].X, Pa[i].Z, work->Paprod[i], work->n, work->tt4, mdata);
#endif

                work->A += wscale * w;
            }

            // amin tracks the Pa[0] position in units of 2 * w, so
            // shifting by w, 2 * U times, is equivalent to shifting
            // by 2 * w, U times.
            amin += U;

#ifdef DO_STAGE2_INV
            foundDuringInv = batch_invert_pt_to_bignum(Pa, work->Pa_inv,
                work->Paprod, mdata, work, 2 * L - shiftdist * U, 2 * L);
#endif
        }
        else
        {
            pa = pairmap_v[mapid] - amin;
            pb = pairmap_u[mapid];

            if (pa >= 2 * L)
            {
                printf("error: invalid A offset: %d,%d,%u\n", pa, pb, amin);
                exit(1);
            }

            if (rprime_map_U[pb] == 0)
            {
                printf("pb=%d doesn't exist\n", pb);
            }

            //if ((((2 * (uint64_t)amin + (uint64_t)pa) * (uint64_t)w - (uint64_t)pb) == 10078477ULL) ||
            //    (((2 * (uint64_t)amin + (uint64_t)pa) * (uint64_t)w + (uint64_t)pb) == 10078477ULL))
            //{
            //    printf("\naccumulated %lu @ amin = %u, pa = %d, pb = %d\n", 
            //        10078477ULL, amin, pa, pb);
            //}

#ifdef DO_STAGE2_INV
            CROSS_PRODUCT_INV;
#else
            CROSS_PRODUCT;
#endif
            work->paired++;
        }
    }

    work->amin = amin;

    return;
}

uint32_t pair(uint32_t* pairmap_v, uint32_t* pairmap_u,
    ecm_work* work, Queue_t** Q, uint32_t* Qrmap, uint32_t* Qmap,
    uint64_t* primes, uint64_t ecm_nump, uint64_t B1, uint64_t B2, int verbose)
{
    int i, j, pid = 0;
    int w = work->D;
    int U = work->U;
    int L = work->L;
    int R = work->R - 3;
    int umax = w * U;
    int64_t q, mq;
    uint64_t amin = work->amin = (B1 + w) / (2 * w);
    uint64_t a, s, ap, u;
    uint32_t pairs = 0;
    uint32_t nump = 0;
    uint32_t mapid = 0;
    uint8_t* flags = NULL;
    int printpairs = 0;
    int testcoverage = 0;
    int printpairmap = 0;

    // gives an index of a queue given a residue mod w
    //printf("Qmap: \n");
    // contains the value of q given an index
    //printf("Qrmap: \n");

    if (testcoverage)
    {
        flags = (uint8_t*)xcalloc((10000 + B2 - B1), sizeof(uint8_t));
    }

    if (verbose > 1)
    {
        printf("commencing pair on range %lu:%lu\n", B1, B2);
    }

    if (printpairmap || printpairs)
    {
        printf("commencing pair at A=%lu\n"
            "w = %u, R = %u, L = %u, U = %d, umax = %u, amin = %u\n",
            2 * (uint64_t)amin * (uint64_t)w, w, R, L, U, umax, amin);
    }

    while (primes[pid] < B1) { pid++; }

    while ((pid < ecm_nump) && (primes[pid] < B2))
    {
        s = primes[pid];
        a = (s + w) / (2 * w);
        nump++;

        //printf("s, a: %lu, %lu\n", s, a);

        while (a >= (amin + L))
        {
            int oldmin = amin;
            amin = amin + L - U;
            //printf("amin now %u\n", amin);

            for (i = 0; i < R; i++)
            {
                int len = Q[i]->len;

                if (Qrmap[i] > w)
                {
                    q = 2 * w - Qrmap[i];
                    for (j = 0; j < len; j++)
                    {
                        ap = dequeue(Q[i]);
                        if ((uint32_t)ap < amin)
                        {
                            pairmap_v[mapid] = 2 * ap - oldmin; // 2 * ap; //2 * (ap - oldmin);
                            pairmap_u[mapid] = q;
                            mapid++;

                            if (testcoverage)
                            {
                                addflag(flags, (2 * ap) * w + q - B1, 0, B2 - B1);
                                addflag(flags, (2 * ap) * w - q - B1, 0, B2 - B1);
                            }
                            if (printpairs)
                            {
                                printf("pair (ap,q):(%lu,%ld)  %lu:%lu\n",
                                    ap, q,
                                    2 * ap * w - q,
                                    2 * ap * w + q);
                            }
                            pairs++;
                        }
                        else
                        {
                            enqueue(Q[i], ap);
                        }
                    }
                }
                else
                {
                    for (j = 0; j < len; j++)
                    {
                        ap = dequeue(Q[i]);
                        if ((uint32_t)ap < amin)
                        {
                            pairmap_v[mapid] = 2 * ap - oldmin; //2 * ap; //2 * (ap - oldmin);
                            pairmap_u[mapid] = Qrmap[i];
                            mapid++;

                            if (testcoverage)
                            {
                                addflag(flags, (2 * ap) * w + Qrmap[i] - B1, 0, B2 - B1);
                                addflag(flags, (2 * ap) * w - Qrmap[i] - B1, 0, B2 - B1);
                            }
                            if (printpairs)
                            {
                                printf("pair (ap,q):(%lu,%u)  %lu:%lu\n",
                                    ap, Qrmap[i],
                                    2 * ap * w - Qrmap[i],
                                    2 * ap * w + Qrmap[i]);
                            }
                            pairs++;
                        }
                        else
                        {
                            enqueue(Q[i], ap);
                        }
                    }
                }
            }
            pairmap_u[mapid] = 0;
            pairmap_v[mapid] = 0;
            mapid++;
        }

        q = s - 2 * a * w;
        if (q < 0)
            mq = abs(q);
        else
            mq = 2 * w - q;

        //printf("q, mq: %d, %d\n", q, mq);

        do
        {
            if (Q[Qmap[mq]]->len > 0)
            {
                ap = dequeue(Q[Qmap[mq]]);
                if (q < 0)
                    u = w * (a - ap) - abs(q);
                else
                    u = w * (a - ap) + q;

                if (u > umax)
                {
                    if (q < 0)
                    {
                        int qq = abs(q);

                        pairmap_v[mapid] = 2 * ap - amin; //2 * ap; //2 * (ap - amin);
                        pairmap_u[mapid] = qq;
                        mapid++;

                        if (testcoverage)
                        {
                            addflag(flags, (2 * ap) * w + qq - B1, 0, B2 - B1);
                            addflag(flags, (2 * ap) * w - qq - B1, 0, B2 - B1);
                        }
                        if (printpairs)
                        {
                            printf("pair (ap,q):(%lu,%ld)  %lu:%lu\n",
                                ap, q,
                                2 * ap * w - qq,
                                2 * ap * w + qq);
                        }
                    }
                    else
                    {
                        int qq = q;

                        if (qq >= w)
                            qq = 2 * w - qq;

                        pairmap_v[mapid] = 2 * ap - amin; //2 * ap; //2 * (ap - amin);
                        pairmap_u[mapid] = qq;
                        mapid++;

                        if (testcoverage)
                        {
                            addflag(flags, (2 * ap) * w + qq - B1, 0, B2 - B1);
                            addflag(flags, (2 * ap) * w - qq - B1, 0, B2 - B1);
                        }
                        if (printpairs)
                        {
                            printf("pair (ap,q):(%lu,%ld)  %lu:%lu\n",
                                ap, q,
                                2 * ap * w - qq,
                                2 * ap * w + qq);
                        }
                    }
                    pairs++;
                }
                else
                {
                    pairmap_v[mapid] = a + ap - amin;
                    pairmap_u[mapid] = u;
                    mapid++;

                    if (testcoverage)
                    {
                        addflag(flags, (a + ap) * w + u - B1, 0, B2 - B1);
                        addflag(flags, (a + ap) * w - u - B1, 0, B2 - B1);
                    }
                    if (printpairs)
                    {
                        printf("pair (a,ap,u):(%lu,%lu,%lu)  %lu:%lu\n",
                            a, ap, u,
                            (a + ap) * w - u,
                            (a + ap) * w + u);
                    }
                    pairs++;
                }
            }
            else
            {
                //printf("queueing a=%lu in Q[%d]\n", a, abs(q));
                if (q < 0)
                {
                    //printf("queueing a=%lu in Q[%u](%u)\n", a, 2 * w + q, Qmap[2 * w + q]);
                    enqueue(Q[Qmap[2 * w + q]], a);
                }
                else
                {
                    //printf("queueing a=%lu in Q[%d]\n", a, q);
                    enqueue(Q[Qmap[q]], a);
                }
                u = 0;
            }
        } while (u > umax);

        pid++;
    }

    //printf("dumping leftovers in queues\n");
    // empty queues
    for (i = 0; i < R; i++)
    {
        //printf("queue %d (%u) has %d elements\n", i, Qrmap[i], Q[i]->len);
        int len = Q[i]->len;
        for (j = 0; j < len; j++)
        {
            ap = dequeue(Q[i]);
            if (Qrmap[i] > w)
            {
                q = 2 * w - Qrmap[i];

                pairmap_v[mapid] = 2 * ap - amin; //2 * ap; //2 * (ap - amin);
                pairmap_u[mapid] = q;
                mapid++;

                if (printpairs)
                {
                    printf("pair (ap,q):(%lu,%ld)  %lu:%lu\n",
                        ap, q,
                        2 * ap * w - q,
                        2 * ap * w + q);
                }
                if (testcoverage)
                {
                    addflag(flags, (2 * ap) * w + q - B1, 0, B2 - B1);
                    addflag(flags, (2 * ap) * w - q - B1, 0, B2 - B1);
                }
            }
            else
            {
                pairmap_v[mapid] = 2 * ap - amin; //2 * ap; // 2 * (ap - amin);
                pairmap_u[mapid] = Qrmap[i];
                mapid++;

                if (printpairs)
                {
                    printf("pair (ap,q):(%lu,%u)  %lu:%lu\n",
                        ap, Qrmap[i],
                        2 * ap * w - Qrmap[i],
                        2 * ap * w + Qrmap[i]);
                }
                if (testcoverage)
                {
                    addflag(flags, (2 * ap) * w + Qrmap[i] - B1, 0, B2 - B1);
                    addflag(flags, (2 * ap) * w - Qrmap[i] - B1, 0, B2 - B1);
                }
            }
            pairs++;
        }
    }

    if (printpairmap)
    {
        printf("%u pairing steps generated\n", mapid);
        amin = (B1 + w) / (2 * w);
        printf("amin is now %lu (A = %lu)\n", amin, 2 * amin * w);
        for (i = 0; i < mapid; i++)
        {
            printf("pair: %uw+/-%u => %lu:%lu\n", pairmap_v[i]+amin, pairmap_u[i],
                (pairmap_v[i]+amin) * w - pairmap_u[i],
                (pairmap_v[i]+amin) * w + pairmap_u[i]);
            if (pairmap_u[i] == 0)
            {
                amin = amin + L - U;
                printf("amin is now %u (A = %u)\n", amin, 2 * amin * w);
            }
        }

        printf("pairmap_v:\n{");
        for (i = 0; i < mapid; i++)
        {
            printf("%u, ", pairmap_v[i]);
        }
        printf("}\n");
        printf("pairmap_u:\n{");
        for (i = 0; i < mapid; i++)
        {
            printf("%u, ", pairmap_u[i]);
        }
        printf("}\n");
    }

    if (testcoverage)
    {
        pid = 0;
        while (primes[pid] < B1) { pid++; }

        int notcovered = 0;
        while ((pid < ecm_nump) && (primes[pid] < B2))
        {
            if (flags[primes[pid] - B1] != 1)
            {
                printf("prime %lu not covered!\n", primes[pid]);
                notcovered++;
            }
            pid++;
        }
        if (notcovered > 0)
        {
            printf("**** %d primes not covered during pairing!\n", notcovered);
            exit(1);
        }
        free(flags);
    }

    if (verbose > 1)
    {
        printf("%u pairs found from %u primes (ratio = %1.2f)\n",
            pairs, nump, (double)pairs / (double)nump);
    }

    if (printpairmap)
    {
        exit(0);
    }

    return mapid;
}

int vec_check_factor(mpz_t Z, mpz_t n, mpz_t f)
{
    //gmp_printf("checking point Z = %Zx against input N = %Zx\n", Z, n);
    mpz_gcd(f, Z, n);

	if (mpz_cmp_ui(f, 1) > 0)
	{
		if (mpz_cmp(f, n) == 0)
		{
            mpz_set_ui(f, 0);
			return 0;
		}
		return 1;
	}
	return 0;
}


#ifdef PARAM1


/* ECM stage 1 in batch mode, for initial point (x:z) with small coordinates,
   such that x and z fit into a mp_limb_t.
   For example we can start with (x=2:y=1) with the curve by^2 = x^3 + ax^2 + x
   with a = 4d-2 and b=16d+2, then we have to multiply by d=(a+2)/4 in the
   duplicates.
   With the change of variable x=b*X, y=b*Y, this curve becomes:
   Y^2 = X^3 + a/b*X^2 + 1/b^2*X.
*/

#define MAX_HEIGHT 32

#if defined(_WIN32)
/* Due to a limitation in GMP on 64-bit Windows, should also
    affect 32-bit Windows, sufficient memory cannot be allocated
    for the batch product s when using primes larger than the following */
#define MAX_B1_BATCH 3124253146UL
#else
/* nth_prime(2^(MAX_HEIGHT-1)) */
#define MAX_B1_BATCH 50685770167ULL
#endif

unsigned int compute_s(mpz_t s, uint64_t * primes, uint64_t B1)
{
    mpz_t acc[MAX_HEIGHT]; /* To accumulate products of prime powers */
    mpz_t ppz;
    unsigned int i, j, it;
    uint64_t pi = 2, pp, maxpp, qi;

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



/* R <- S*m/B mod modulus where m fits in a mp_limb_t.
   Here S (w in dup_add_batch1) is the result of a subtraction,
   thus with the notations from http://www.loria.fr/~zimmerma/papers/norm.pdf
   we have S < 2 \alpha N.
   Then R < (2 \alpha N \beta + \beta N) = (2 \alpha + 1) N.
   This result R is used in an addition with u being the result of a squaring
   thus u < \alpha N, which gives a result < (3 \alpha + 1) N.
   Finally this result is used in a multiplication with another operand less
   than 2 \alpha N, thus we want:
   ((2 \alpha) (3 \alpha + 1) N^2 + \beta N)/\beta \leq \alpha N, i.e.,
   2 \alpha (3 \alpha + 1) \varepsilon + 1 \leq \alpha
   This implies \varepsilon \leq 7/2 - sqrt(3)/2 ~ 0.0359, in which case
   we can take \alpha = 2/3*sqrt(3)+1 ~ 2.1547.
   In that case no adjustment is needed in mpresn_mul_1.
   However we prefer to keep the adjustment here, to allow a larger set of
   inputs (\varepsilon \leq 1/16 = 0.0625 instead of 0.0359).
*/
void
vecmulmod_1(vec_bignum_t* S, uint64_t* m, vec_bignum_t* R, vec_bignum_t* modulus, vec_bignum_t* t, vec_monty_t* mdata)
{
    //mp_ptr t1 = PTR(modulus->temp1);
    //mp_ptr t2 = PTR(modulus->temp2);
    //mp_size_t n = ABSIZ(modulus->orig_modulus);
    //mp_limb_t q;
    //
    //{
    //    t1[n] = mpn_mul_1(t1, PTR(S), n, m);
    //    q = t1[0] * modulus->Nprim[0];
    //    t2[n] = mpn_mul_1(t2, PTR(modulus->orig_modulus), n, q);
    //
    //    q = mpn_add_n(PTR(R), t1 + 1, t2 + 1, n);
    //    q += mpn_add_1(PTR(R), PTR(R), n, t1[0] != 0);
    //
    //    while (q != 0)
    //        q -= mpn_sub_n(PTR(R), PTR(R), PTR(modulus->orig_modulus), n);
    //}
    vecmulmod52_1(S, m, R, modulus, t, mdata);
}


/* (x1:z1) <- 2(x1:z1)
   (x2:z2) <- (x1:z1) + (x2:z2)
   assume (x2:z2) - (x1:z1) = (2:1)
   Uses 4 modular multiplies and 4 modular squarings.
   Inputs are x1, z1, x2, z2, d, n.
   Use two auxiliary variables: t, w (it seems using one only is not possible
   if all mpresn_mul and mpresn_sqr calls don't overlap input and output).

   In the batch 1 mode, we pass d_prime such that the actual d is d_prime/beta.
   Since beta is a square, if d_prime is a square (on 64-bit machines),
   so is d.
   In mpresn_mul_1, we multiply by d_prime = beta*d and divide by beta.
*/
static void
dup_add_batch1(vec_bignum_t* x1, vec_bignum_t* z1, vec_bignum_t* x2, vec_bignum_t* z2,
    vec_bignum_t* t, vec_bignum_t* w, vec_bignum_t* s, uint64_t* d_prime, vec_bignum_t* n, vec_monty_t* mdata)
{
    //vec_bignum_t *t2;
    //t2 = vecInit();
    //memcpy(t2->data, d_prime, VECLEN * sizeof(base_t));

    /* active: x1 z1 x2 z2 */
    //mpresn_addsub(w, z1, x1, z1, n); /* w = x1+z1, z1 = x1-z1 */
    vecaddsubmod_ptr(x1, z1, w, z1, mdata);

    /* active: w z1 x2 z2 */
    //mpresn_addsub(x1, x2, x2, z2, n); /* x1 = x2+z2, x2 = x2-z2 */
    vecaddsubmod_ptr(x2, z2, x1, x2, mdata);
    /* active: w z1 x1 x2 */

    //mpresn_mul(z2, w, x2, n); /* w = (x1+z1)(x2-z2) */
    vecmulmod_ptr(w, x2, z2, n, s, mdata);
    /* active: w z1 x1 z2 */
    //mpresn_mul(x2, z1, x1, n); /* x2 = (x1-z1)(x2+z2) */
    vecmulmod_ptr(z1, x1, x2, n, s, mdata);
    /* active: w z1 x2 z2 */
    //mpresn_sqr(t, z1, n);    /* t = (x1-z1)^2 */
    vecsqrmod_ptr(z1, t, n, s, mdata);
    /* active: w t x2 z2 */
    //mpresn_sqr(z1, w, n);    /* z1 = (x1+z1)^2 */
    vecsqrmod_ptr(w, z1, n, s, mdata);
    /* active: z1 t x2 z2 */

    //mpresn_mul(x1, z1, t, n); /* xdup = (x1+z1)^2 * (x1-z1)^2 */
    vecmulmod_ptr(z1, t, x1, n, s, mdata);
    /* active: x1 z1 t x2 z2 */

    //mpresn_sub(w, z1, t, n);   /* w = (x1+z1)^2 - (x1-z1)^2 */
    vecsubmod_ptr(z1, t, w, mdata);
    /* active: x1 w t x2 z2 */

    //mpresn_mul_1(z1, w, d_prime, n); /* z1 = d * ((x1+z1)^2 - (x1-z1)^2) */
    vecmulmod_1(w, d_prime, z1, n, s, mdata);
    //vecmulmod_ptr(w, t2, z1, n, s, mdata);
    /* active: x1 z1 w t x2 z2 */

    //mpresn_add(t, t, z1, n);  /* t = (x1-z1)^2 - d* ((x1+z1)^2 - (x1-z1)^2) */
    vecaddmod_ptr(z1, t, t, mdata);
    //vecsubmod_ptr(z1, t, t, mdata);

    /* active: x1 w t x2 z2 */
    //mpresn_mul(z1, w, t, n); /* zdup = w * [(x1-z1)^2 - d* ((x1+z1)^2 - (x1-z1)^2)] */
    vecmulmod_ptr(w, t, z1, n, s, mdata);
    /* active: x1 z1 x2 z2 */

    //mpresn_addsub(w, z2, x2, z2, n);
    vecaddsubmod_ptr(x2, z2, w, z2, mdata);
    /* active: x1 z1 w z2 */

    //mpresn_sqr(x2, w, n);
    vecsqrmod_ptr(w, x2, n, s, mdata);
    /* active: x1 z1 x2 z2 */
    //mpresn_sqr(w, z2, n);
    vecsqrmod_ptr(z2, w, n, s, mdata);
    /* active: x1 z1 x2 w */
    //mpresn_add(z2, w, w, n);
    vecaddmod_ptr(w, w, z2, mdata);

    //vecFree(t2);
}



/* Input: x is initial point
          A is curve parameter in Montgomery's form:
          g*y^2*z = x^3 + a*x^2*z + x*z^2
          n is the number to factor
      B1 is the stage 1 bound
   Output: If a factor is found, it is returned in x.
           Otherwise, x contains the x-coordinate of the point computed
           in stage 1 (with z coordinate normalized to 1).
       B1done is set to B1 if stage 1 completed normally,
       or to the largest prime processed if interrupted, but never
       to a smaller value than B1done was upon function entry.
   Return value: ECM_FACTOR_FOUND_STEP1 if a factor, otherwise
           ECM_NO_FACTOR_FOUND
*/
/*
For now we don't take into account go stop_asap and chkfilename
*/
int vec_ecm_stage1_batch(mpz_t f, ecm_work* work, ecm_pt* P, vec_bignum_t* A,
    vec_bignum_t* n, uint64_t B1, mpz_t s, vec_monty_t* mdata)
{
    uint64_t* d_1;
    mpz_t d_2;

    vec_bignum_t* x1, * z1, * x2, * z2;
    uint64_t i;
    vec_bignum_t* t, * u;
    int ret = 0;

    x1 = vecInit(words);
    z1 = vecInit(words);
    x2 = vecInit(words);
    z2 = vecInit(words);
    t = vecInit(words);
    u = vecInit(words);

    d_1 = (uint64_t*)xmalloc_align(VECLEN * sizeof(uint64_t));

    vecCopy(P->X, x1);
    vecCopy(P->Z, z1);
    vecCopy(work->pt2.X, x2);
    vecCopy(work->pt2.Z, z2);
    memcpy(d_1, work->s->data, VECLEN * sizeof(base_t));

    //printf("d_1 = ");
    //for (i = 0; i < VECLEN; i++)
    //    printf("%lu, ", d_1[i]);
    //printf("\n");

    /* initialize P */
    //mpres_set(x1, x, n);
    //mpres_set_ui(z1, 1, n); /* P1 <- 1P */

    /* Compute d=(A+2)/4 from A and d'=B*d thus d' = 2^(GMP_NUMB_BITS-2)*(A+2) */
    //mpres_get_z(u, A, n);
    //mpz_add_ui(u, u, 2);
    //mpz_mul_2exp(u, u, GMP_NUMB_BITS - 2);
    //mpres_set_z_for_gcd(u, u, n); /* reduces u mod n */
    //if (mpz_size(u) > 1)
    //{
    //    mpres_get_z(u, A, n);
    //    gmp_fprintf(stderr,
    //        "Error, 2^%d*(A+2) should fit in a mp_limb_t, A=%Zd\n",
    //        GMP_NUMB_BITS - 2, u);
    //    return -1;
    //}
    //d_1 = mpz_getlimbn(u, 0);

    ///* Compute 2P : no need to duplicate P, the coordinates are simple. */
    //mpres_set_ui(x2, 9, n);
    ///* here d = d_1 / GMP_NUMB_BITS */
    ///* warning: mpres_set_ui takes an unsigned long which has only 32 bits
    //    on Windows, while d_1 might have 64 bits */
    //if (mpz_size(u) != 1 || mpz_getlimbn(u, 0) != d_1_1)
    //{
    //    printf("problem with mpres_set_ui\n");
    //    exit(-2);
    //}
    //mpres_set_z(z2, u, n);
    //mpres_div_2exp(z2, z2, GMP_NUMB_BITS, n);
    //
    //mpres_mul_2exp(z2, z2, 6, n);
    //mpres_add_ui(z2, z2, 8, n); /* P2 <- 2P = (9 : : 64d+8) */

    /* invariant: if j represents the upper bits of s,
       then P1 = j*P and P2=(j+1)*P */

       //mpresn_pad(x1, n);
       //mpresn_pad(z1, n);
       //mpresn_pad(x2, n);
       //mpresn_pad(z2, n);

       /* now perform the double-and-add ladder */
    for (i = mpz_sizeinbase(s, 2) - 1; i-- > 0;)
    {
        if (mpz_tstbit(s, i) == 0) /* (j,j+1) -> (2j,2j+1) */
            /* P2 <- P1+P2    P1 <- 2*P1 */
            dup_add_batch1(x1, z1, x2, z2, t, u, mdata->mtmp1, d_1, n, mdata);
        else /* (j,j+1) -> (2j+1,2j+2) */
            /* P1 <- P1+P2     P2 <- 2*P2 */
            dup_add_batch1(x2, z2, x1, z1, t, u, mdata->mtmp1, d_1, n, mdata);
    }

    //mpresn_unpad(x1);
    //mpresn_unpad(z1);

    // this sets x = x1/z1, but our vececm expects
    // them to be separate.
    //if (!mpres_invert(u, z1, n)) /* Factor found? */
    //{
    //    mpres_gcd(f, z1, n);
    //    ret = 1;
    //}
    //mpres_mul(x, x1, u, n);
    //vecmulmod_ptr(x1, u, x, mdata->mtmp1, n, mdata);

    vecCopy(x1, P->X);
    vecCopy(z1, P->Z);

    vecFree(x1);
    vecFree(z1);
    vecFree(x2);
    vecFree(z2);
    vecFree(t);
    vecFree(u);
    align_free(d_1);

    return ret;
}

#endif


