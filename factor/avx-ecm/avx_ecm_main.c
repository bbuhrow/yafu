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

#include "avx_ecm.h"
#include "soe.h"
#include "gmp.h"
#include "ytools.h"
#include "nfs_impl.h"
#include <math.h>

// performance comparison
// http://www.mersenneforum.org/showthread.php?t=16480&page=20
// http://www.mersenneforum.org/showthread.php?t=5722&page=122

// t-level estimate
// http://mersenneforum.org/showpost.php?p=427989&postcount=2429


static int debugctr = 0;
void thread_init(thread_data_t* tdata, vec_monty_t* mdata, uint64_t B1, uint64_t B2);

void extract_bignum_from_vec_to_mpz(mpz_t dest, vec_bignum_t *vec_src, int num, int sz)
{
    int j;

    if (dest == NULL)
    {
        printf("invalid dest address in extract_vec_bignum_from_vec_to_mpz\n");
    }

    mpz_set_ui(dest, 0);
    for (j = sz - 1; j >= 0; j--)
    {
        
#if defined(WIN64) && (DIGITBITS == 52)
        mpz_mul_2exp(dest, dest, 20);
        mpz_add_ui(dest, dest, vec_src->data[num + j * VECLEN] >> 32);
        mpz_mul_2exp(dest, dest, 32);
        mpz_add_ui(dest, dest, (vec_src->data[num + j * VECLEN]) & 0xffffffff);
#else
        mpz_mul_2exp(dest, dest, DIGITBITS);
        mpz_add_ui(dest, dest, vec_src->data[num + j * VECLEN]);
#endif
        //gmp_printf("word is %016llx, dest is now %Zx\n", vec_src->data[num + j * VECLEN], dest);
    }

    return;
}

void broadcast_mpz_to_vec(vec_bignum_t *vec_dest, mpz_t src)
{
    mpz_t src_cp;
    int i, j;

    mpz_init(src_cp);
    mpz_set(src_cp, src);

    i = 0;
    vec_dest->size = 0;
    while (mpz_cmp_ui(src_cp, 0) > 0)
    {
        base_t thisword = mpz_get_ui(src_cp) & VEC_MAXDIGIT;
        for (j = 0; j < VECLEN; j++)
        {
            vec_dest->data[j + i * VECLEN] = thisword;
        }
        vec_dest->size++;
        i++;
        mpz_tdiv_q_2exp(src_cp, src_cp, DIGITBITS);
    }

    vec_dest->size = vec_dest->WORDS_ALLOC;
    mpz_clear(src_cp);
    return;
}

void insert_mpz_to_vec(vec_bignum_t *vec_dest, mpz_t src, int lane)
{
    mpz_t src_cp;
    int i;

    mpz_init(src_cp);
    mpz_set(src_cp, src);

    i = 0;
    //vec_dest->size = 0;
    while (mpz_cmp_ui(src_cp, 0) > 0)
    {
        base_t thisword = mpz_get_ui(src_cp) & VEC_MAXDIGIT;
        vec_dest->data[lane + i * VECLEN] = thisword;
        i++;
        mpz_tdiv_q_2exp(src_cp, src_cp, DIGITBITS);
    }

    vec_dest->size = vec_dest->WORDS_ALLOC; // MAX(vec_dest->size, i);
    mpz_clear(src_cp);
    return;
}

void vec_ecm_main(fact_obj_t* fobj, uint32_t numcurves, uint64_t B1,
    uint64_t B2, int threads, int* numfactors, int verbose, 
    int save_b1, uint32_t *curves_run)
{
    thread_data_t *tdata;
	vec_bignum_t **f, *n;
    mpz_t g, r, N;
	uint32_t *siglist;
	uint32_t numcurves_per_thread;
	uint64_t b1 = B1;
	uint32_t i, j;
	char **nextptr;
	vec_monty_t *montyconst;
	int pid = getpid();
    uint64_t limit;
    int size_n, isMersenne = 0, forceNoMersenne = 0;
    uint64_t sigma = fobj->ecm_obj.sigma;
    uint32_t maxbits, nwords, nblocks;

	// timing variables
	struct timeval stopt;	// stop time of this job
	struct timeval startt;	// start time of this job
	double t_time;

	gettimeofday(&startt, NULL);

    if (pid <= 0)
    {
        pid = startt.tv_usec;
    }

    if (verbose > 1)
    {
        printf("process id is %d\n", pid);
    }

    *numfactors = 0;
    mpz_init(g);
    mpz_init(r);
    mpz_init(N);

    mpz_set(N, fobj->ecm_obj.gmp_n);

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
        gmp_printf("commencing parallel ecm on %Zd with %d threads\n", N, threads);
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

    // now that we know NWORDS, can allocate the monty structure.
    montyconst = vec_monty_alloc(nwords);
    tdata = (thread_data_t*)malloc(threads * sizeof(thread_data_t));

    tdata[0].MAXBITS = montyconst->MAXBITS = maxbits;
    tdata[0].NBLOCKS = montyconst->NBLOCKS = nblocks;
    tdata[0].NWORDS = montyconst->NWORDS = nwords;

    if (verbose > 1)
    {
        printf("ECM has been configured with DIGITBITS = %u, VECLEN = %d, GMP_LIMB_BITS = %d\n",
            DIGITBITS, VECLEN, GMP_LIMB_BITS);

        printf("Choosing MAXBITS = %u, NWORDS = %d, NBLOCKS = %d based on input size %d\n",
            tdata[0].MAXBITS, tdata[0].NWORDS, tdata[0].NBLOCKS, size_n);

        if (sigma > 0)
        {
            printf("starting with sigma = %lu\n", sigma);
        }
    }

    if (isMersenne != 0)
    {
        montyconst->isMersenne = isMersenne;
        montyconst->nbits = size_n;
        montyconst->use_vnbits = 0;             // all lanes use the same N.
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
        montyconst->use_vnbits = 0;             // all lanes use the same N.
        mpz_set_ui(r, 1);
        mpz_mul_2exp(r, r, DIGITBITS * tdata[0].NWORDS);
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

    

    tdata[0].STAGE1_MAX = b1;
    tdata[0].DO_STAGE2 = 1;
    if (B2 == (uint64_t)-1)
	{
        tdata[0].STAGE2_MAX = tdata[0].STAGE1_MAX;
        tdata[0].DO_STAGE2 = 0;
	}
    else if (B2 == 0)
    {
        tdata[0].STAGE2_MAX = 100ULL * (uint64_t)b1;
    }
    else
    {
        tdata[0].STAGE2_MAX = B2;
    }

    // not currently used, but the mechanisms exist to pass in an initial
    // list of primes to use here.
    tdata[0].PRIME_RANGE = 100000000;

    if (numcurves < threads)
    {
        numcurves = threads;
    }
	
    if (verbose > 1)
    {
        printf("configuring avx-ecm with %d threads\n", threads); fflush(stdout);
    }

	numcurves_per_thread = numcurves / threads + (numcurves % threads != 0);
	numcurves = numcurves_per_thread * threads;

    if (verbose > 1)
    {
        printf("Input has %d bits, using %d threads (%d curves/thread)\n",
            mpz_sizeinbase(N, 2), threads, numcurves_per_thread);
        printf("Processing in batches of %u primes\n", tdata[0].PRIME_RANGE);
    }

    //gmp_printf("n = %Zx\n", gmpn);
    //gmp_printf("rhat = %Zx\n", montyconst->rhat);
    //gmp_printf("nhat = %Zx\n", montyconst->nhat);
    //gmp_printf("one = %Zx\n", r);
    //printf("rho = %016llx\n", mpz_get_ui(montyconst->nhat) & MAXDIGIT);
    for (i = 0; i < VECLEN; i++)
    {
        montyconst->vrho[i] = mpz_get_ui(montyconst->nhat) & VEC_MAXDIGIT;
    }

    // test karatsuba
    if (0)
    {
        mpz_t a, b, c;
        gmp_randstate_t grnd;
        uint64_t lcg_state;
        uint32_t sz;

        mpz_init(a);
        mpz_init(b);
        mpz_init(c);
        gmp_randinit_default(grnd);

        vec_bignum_t *va, *vb, *vc, * vkc, *vn, *vt;

        va = vecInit(montyconst->NWORDS);
        vb = vecInit(montyconst->NWORDS);
        vc = vecInit(montyconst->NWORDS);
        vkc = vecInit(montyconst->NWORDS);
        vn = vecInit(montyconst->NWORDS);
        vt = vecInit(montyconst->NWORDS);

        mpz_urandomb(a, grnd, montyconst->MAXBITS);
        mpz_nextprime(a, a);
        lcg_state = mpz_get_ui(a);
        sz = lcg_rand_32_range(montyconst->MAXBITS/2, montyconst->MAXBITS, &lcg_state);
        if ((sz & 1) == 0)
            sz--;

        printf("\nsize is %u\n", sz);

        montyconst->isMersenne = 1;
        montyconst->nbits = sz;
        montyconst->use_vnbits = 0;             // all lanes use the same N.
        

        // do a non-kara mul for reference
        vecmulmod_ptr = &vecmulmod52_mersenne;
        vecsqrmod_ptr = &vecsqrmod52_mersenne;
        //vecmulmod_ptr = &veckmul_mersenne;
        //vecsqrmod_ptr = &vecksqr_mersenne;
        vecaddmod_ptr = &vecaddmod52_mersenne;
        vecsubmod_ptr = &vecsubmod52_mersenne;
        vecaddsubmod_ptr = &vec_simul_addsub52_mersenne;

        montyconst->isMersenne = 1;
        mpz_set_ui(c, 1);
        mpz_mul_2exp(c, c, sz);
        mpz_sub_ui(c, c, 1);
        insert_mpz_to_vec(vn, c, i);
        mpz_set(montyconst->nhat, c);           // remember input N
        // do all math w.r.t the Mersenne number
        broadcast_mpz_to_vec(montyconst->n, c);
        broadcast_mpz_to_vec(montyconst->vnhat, montyconst->nhat);
        mpz_set_ui(r, 1);
        broadcast_mpz_to_vec(montyconst->one, r);

        for (i = 0; i < VECLEN; i++)
        {
            mpz_urandomb(a, grnd, sz);
            mpz_nextprime(a, a);
            insert_mpz_to_vec(va, a, i);
            mpz_urandomb(b, grnd, sz);
            mpz_nextprime(b, b);
            insert_mpz_to_vec(vb, b, i);

            insert_mpz_to_vec(vn, c, i);

            gmp_printf("a%d: %Zx\n", i, a);
            gmp_printf("b%d: %Zx\n", i, b);
            gmp_printf("n%d: %Zx\n", i, c);
        }

        vecmulmod_ptr(va, vb, vc, vn, vt, montyconst);

        vecmulmod_ptr = &veckmul_mersenne;
        vecsqrmod_ptr = &vecksqr_mersenne;

        vecmulmod_ptr(va, vb, vkc, vn, vt, montyconst);

        for (i = 0; i < VECLEN; i++)
        {
            extract_bignum_from_vec_to_mpz(a, vc, i, montyconst->NWORDS);
            extract_bignum_from_vec_to_mpz(b, vkc, i, 2*montyconst->NWORDS);

            gmp_printf("sch%d: %Zx\n", i, a);
            gmp_printf("kar%d: %Zx\n", i, b);
        }

        exit(1);
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
            //vecmulmod_ptr = &veckmul_mersenne;
            //vecsqrmod_ptr = &vecksqr_mersenne;
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
            if (tdata[0].MAXBITS == 416)
            {
                vecmulmod_ptr = &vecmulmod52_fixed416_bfips;
                vecsqrmod_ptr = &vecsqrmod52_fixed416_bfips;
                vecaddsubmod_ptr = &vec_simul_addsub52_fixed1040;
                if (verbose > 1)
                {
                    printf("Using 416-bit mul/sqr core\n");
                }
            }
            else if (tdata[0].MAXBITS == 624)
            {
                vecmulmod_ptr = &vecmulmod52_fixed624_bfips;
                vecsqrmod_ptr = &vecsqrmod52_fixed624_bfips;
                vecaddsubmod_ptr = &vec_simul_addsub52_fixed1040;
                if (verbose > 1)
                {
                    printf("Using 624-bit mul/sqr core\n");
                }
            }
            else if (tdata[0].MAXBITS == 832)
            {
                vecmulmod_ptr = &vecmulmod52_fixed832_bfips;
                vecsqrmod_ptr = &vecsqrmod52_fixed832_bfips;
                vecaddsubmod_ptr = &vec_simul_addsub52_fixed1040;
                if (verbose > 1)
                {
                    printf("Using 832-bit mul/sqr core\n");
                }
            }
            else if (tdata[0].MAXBITS == 1040)
            {
                vecmulmod_ptr = &vecmulmod52_fixed1040_bfips;
                vecsqrmod_ptr = &vecsqrmod52_fixed1040_bfips;
                vecaddsubmod_ptr = &vec_simul_addsub52_fixed1040;
                if (verbose > 1)
                {
                    printf("Using 1040-bit mul/sqr core\n");
                }
            }
            else
            {
                vecmulmod_ptr = &vecmulmod52;
                vecsqrmod_ptr = &vecsqrmod52;
                vecaddsubmod_ptr = &vec_simul_addsub52;
            }
            vecaddmod_ptr = &vecaddmod52;
            vecsubmod_ptr = &vecsubmod52;
            
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

	gettimeofday(&stopt, NULL);

    uint32_t seed1, seed2;
    get_random_seeds(&seed1, &seed2);
    for (i = 0; i < threads; i++)
    {
        int j;
        thread_init(&tdata[i], montyconst, tdata[0].STAGE1_MAX, tdata[0].STAGE2_MAX);
        tdata[i].curves = numcurves_per_thread;
        tdata[i].tid = i;
        tdata[i].numfactors = 0;
        // each thread gets a random LCG state that folds in the current pid,
        // this session's rand_t, and the current microsecond timer.  Plus the
        // thread id to make it unique per-thread.
		tdata[i].lcg_state = hash64(stopt.tv_usec) + hash64(pid) + hash64(i+1) + 
            seed1 + (uint64_t)seed2 << 32;
        tdata[i].total_threads = threads;
        tdata[i].verbose = verbose;
        tdata[i].save_b1 = save_b1;
        tdata[i].STAGE1_MAX = tdata[0].STAGE1_MAX;
        tdata[i].STAGE2_MAX = tdata[0].STAGE2_MAX;
        tdata[i].DO_STAGE2 = tdata[0].DO_STAGE2;
        tdata[i].PRIME_RANGE = tdata[0].PRIME_RANGE;
        tdata[i].MAXBITS = maxbits;
        tdata[i].NBLOCKS = nblocks;
        tdata[i].NWORDS = nwords;
        
        if (i == 0)
        {
            uint32_t D = tdata[i].work->D;
            int k;

            tdata[i].pairmap_v = (uint32_t*)calloc(tdata[0].PRIME_RANGE, sizeof(uint32_t));
            tdata[i].pairmap_u = (uint32_t*)calloc(tdata[0].PRIME_RANGE, sizeof(uint32_t));

            tdata[i].Qmap = (uint32_t*)malloc(2 * D * sizeof(uint32_t));
            tdata[i].Qrmap = (uint32_t*)malloc(2 * D * sizeof(uint32_t));

            for (j = 0, k = 0; k < 2 * D; k++)
            {
                if (spGCD(k, 2 * D) == 1)
                {
                    tdata[i].Qmap[k] = j;
                    tdata[i].Qrmap[j++] = k;
                }
                else
                {
                    tdata[i].Qmap[k] = (uint32_t)-1;
                }
            }

            for (k = j; k < 2 * D; k++)
            {
                tdata[i].Qrmap[k] = (uint32_t)-1;
            }

            tdata[i].Q = (Queue_t * *)malloc(j * sizeof(Queue_t*));
            for (k = 0; k < j; k++)
            {
                tdata[i].Q[k] = newQueue(D, 0);
            }
        }
        else
        {
            tdata[i].pairmap_v = tdata[0].pairmap_v;
            tdata[i].pairmap_u = tdata[0].pairmap_u;
        }


        if (sigma > 0)
        {
            for (j = 0; j < VECLEN; j++)
            {
                tdata[i].sigma[j] = sigma + VECLEN * i + j;
            }
        }
        else
        {
            for (j = 0; j < VECLEN; j++)
            {
                tdata[i].sigma[j] = 0;
            }
        }
    }

	gettimeofday(&stopt, NULL);
    t_time = ytools_difftime(&startt, &stopt);

    if (verbose > 1)
    {
        printf("Initialization took %1.4f seconds.\n", t_time);
    }

    FILE* flog = NULL;
    if (fobj->LOGFLAG && (strcmp(fobj->flogname, "") != 0))
    {
        flog = fopen(fobj->flogname, "a");
        if (flog == NULL)
        {
            printf("fopen error: %s\n", strerror(errno));
            printf("could not open %s for appending\n", fobj->flogname);
            return;
        }

        char* s;
        s = mpz_get_str(NULL, 10, fobj->ecm_obj.gmp_n);
        logprint(flog, "ecm: commencing %u curves using AVX-ECM method on %s, ",
            tdata[0].curves * threads, s);
        free(s);

        int tmp = fobj->ecm_obj.stg2_is_default;
        uint64_t tmp2 = fobj->ecm_obj.B2;
        fobj->ecm_obj.stg2_is_default = 0;
        fobj->ecm_obj.B2 = B2;

        print_B1B2(fobj, flog);

        fobj->ecm_obj.stg2_is_default = tmp;
        fobj->ecm_obj.B2 = tmp2;
        fprintf(flog, "\n");
        fclose(flog);
    }

    // top level ECM 
    vececm(tdata);    

    if (verbose > 1)
	    printf("\n");

    *curves_run = tdata[0].curves;

    // record factors and clean up thread data
    if (fobj->LOGFLAG && (strcmp(fobj->flogname, "") != 0))
    {
        flog = fopen(fobj->flogname, "a");
        if (flog == NULL)
        {
            printf("fopen error: %s\n", strerror(errno));
            printf("could not open %s for appending\n", fobj->flogname);
            return;
        }

        logprint(flog, "ecm: finished %u curves using AVX-ECM method on C%d input, ",
            *curves_run, gmp_base10(fobj->ecm_obj.gmp_n));

        int tmp = fobj->ecm_obj.stg2_is_default;
        uint64_t tmp2 = fobj->ecm_obj.B2;
        fobj->ecm_obj.stg2_is_default = 0;
        fobj->ecm_obj.B2 = B2;

        print_B1B2(fobj, flog);

        fobj->ecm_obj.stg2_is_default = tmp;
        fobj->ecm_obj.B2 = tmp2;
        fprintf(flog, "\n");
    }

    int total_factors = 0;
	for (i = 0; i < threads; i++)
	{
        //printf("thread %d found %d factors\n", i, tdata[i].numfactors);
        for (j = 0; j < tdata[i].numfactors; j++)
        {
            //factors = ecm_add_factor(factors, numfactors, tdata[i].factors[j].factor,
            //    tdata[i].factors[j].stg_id, tdata[i].factors[j].thread_id,
            //    tdata[i].factors[j].vec_id, tdata[i].factors[j].curve_id, 
            //    tdata[i].factors[j].sigma);
            int fid;

            // AVX-ECM can find the same factor multiple times simultaneously, so
            // check to make sure this factor still divides the input.
            // the factor can also be a composite where we have found part of it
            // already but not another.  So a divisibility check isn't enough...
            // use gcd. But now it depends on the order in which we find factors.
            // if we find a large composite factor first and remove it, then 
            // the smaller prime factors won't be added after.  To solve that 
            // dilemma, we scroll through the rest of the found factors list and
            // if any divide this one, we remove it from this one.
            mpz_gcd(g, fobj->ecm_obj.gmp_n, tdata[i].factors[j].factor);

            int k;
            for (k = j + 1; k < tdata[i].numfactors; k++)
            {
                if (mpz_divisible_p(g, tdata[i].factors[k].factor))
                {
                    mpz_tdiv_q(g, g, tdata[i].factors[k].factor);
                }
            }

            if (mpz_cmp_ui(g, 1) == 0)
            {
                mpz_clear(tdata[i].factors[j].factor);
                continue;
            }

            //tdata[i].factors[j].factor, 
            fid = add_to_factor_list(fobj->factors, g, 
                fobj->VFLAG, fobj->NUM_WITNESSES);

            total_factors++;
            //fid = fobj->factors->num_factors - 1;

            fobj->factors->factors[fid].tid = tdata[i].factors[j].thread_id;
            fobj->factors->factors[fid].curve_num = tdata[i].factors[j].curve_id;
            fobj->factors->factors[fid].sigma = tdata[i].factors[j].sigma;
            fobj->factors->factors[fid].vid = tdata[i].factors[j].vec_id;
            fobj->factors->factors[fid].method = tdata[i].factors[j].stg_id;

            if (is_mpz_prp(g, fobj->NUM_WITNESSES))
            {
                if (fobj->VFLAG > 0)
                    gmp_printf("\necm: found prp%d factor = %Zd\n",
                        gmp_base10(g), g);

                if (fobj->LOGFLAG && (strcmp(fobj->flogname, "") != 0))
                {
                    char* s = mpz_get_str(NULL, 10, g);
                    logprint(flog, "prp%d = %s (curve=%d stg=%d B1=%u B2=%lu sigma=%lu thread=%d vecpos=%d)\n",
                        gmp_base10(g),
                        s,
                        fobj->factors->factors[fid].curve_num, 
                        fobj->factors->factors[fid].method,
                        fobj->ecm_obj.B1, B2, 
                        fobj->factors->factors[fid].sigma,
                        fobj->factors->factors[fid].tid, 
                        fobj->factors->factors[fid].vid);
                    free(s);
                }
            }
            else
            {
                if (fobj->VFLAG > 0)
                    gmp_printf("\necm: found c%d factor = %Zd\n",
                        gmp_base10(g), g);

                if (fobj->LOGFLAG && (strcmp(fobj->flogname, "") != 0))
                {
                    char* s = mpz_get_str(NULL, 10, tdata[i].factors[j].factor);
                    logprint(flog, "c%d = %s (curve=%d stg=%d B1=%u B2=%lu sigma=%lu thread=%d vecpos=%d)\n",
                        gmp_base10(g),
                        s,
                        fobj->factors->factors[fid].curve_num,
                        fobj->factors->factors[fid].method,
                        fobj->ecm_obj.B1, B2,
                        fobj->factors->factors[fid].sigma,
                        fobj->factors->factors[fid].tid,
                        fobj->factors->factors[fid].vid);
                    free(s);
                }
            }

            mpz_tdiv_q(fobj->ecm_obj.gmp_n, fobj->ecm_obj.gmp_n, g); // tdata[i].factors[j].factor);

            mpz_clear(tdata[i].factors[j].factor);
        }

        // if no factors are found the structure will not
        // get allocated at all.
        if (tdata[i].numfactors > 0)
            free(tdata[i].factors);
       
        vec_ecm_work_free(tdata[i].work);
        vec_ecm_pt_free(tdata[i].P);
        vec_monty_free(tdata[i].mdata);
        free(tdata[i].mdata);
        free(tdata[i].work);
        free(tdata[i].P);
        free(tdata[i].sigma);
	}

    // if we found a factor but the input is not 1, check to see if it's prime
    if ((total_factors > 0) && mpz_cmp_ui(fobj->ecm_obj.gmp_n, 1) > 0)
    {
        if (is_mpz_prp(fobj->ecm_obj.gmp_n, fobj->NUM_WITNESSES))
        {
            add_to_factor_list(fobj->factors, fobj->ecm_obj.gmp_n,
                fobj->VFLAG, fobj->NUM_WITNESSES);

            if (fobj->VFLAG > 0)
                gmp_printf("\necm: found prp%d (co)factor = %Zd\n",
                    gmp_base10(fobj->ecm_obj.gmp_n), fobj->ecm_obj.gmp_n);

            if (fobj->LOGFLAG && (strcmp(fobj->flogname, "") != 0))
            {
                char* s = mpz_get_str(NULL, 10, fobj->ecm_obj.gmp_n);
                logprint(flog, "prp%d = %s (cofactor)\n",
                    gmp_base10(fobj->ecm_obj.gmp_n), s);
                free(s);
            }

            mpz_set_ui(fobj->ecm_obj.gmp_n, 1);
        }
    }

    if ((strcmp(fobj->flogname, "") != 0) && (flog != NULL))
    {
        fclose(flog);
    }

    // clean up local/global data
    mpz_clear(g);
    mpz_clear(r);
    vec_monty_free(montyconst);
	free(montyconst);
    free(tdata);

	return;
}

void thread_init(thread_data_t *tdata, vec_monty_t *mdata, uint64_t B1, uint64_t B2)
{    
    int i, j;

    tdata->work = (ecm_work *)malloc(sizeof(ecm_work));
    tdata->P = (ecm_pt *)malloc(sizeof(ecm_pt));
    tdata->sigma = (uint64_t *)malloc(VECLEN * sizeof(uint64_t));
    uint32_t D = tdata->work->D = 2310;
    tdata->work->STAGE1_MAX = B1;
    tdata->work->STAGE2_MAX = B2;

	// decide on the stage 2 parameters.  Larger U means
	// more memory and more setup overhead, but more prime pairs.
	// Smaller U means the opposite.  find a good balance.
    if (B1 <= 4096)
    {
        D = tdata->work->D = 1155;
    }

    if (B1 <= 2048)
    {
        D = tdata->work->D = 385;
    }

    if (B1 <= 512)
    {
        D = tdata->work->D = 210;
    }

    if (B1 <= 256)
    {
        D = tdata->work->D = 120;
    }

    if (B1 <= 128)
    {
        D = tdata->work->D = 60;
    }

    if (B1 <= 60)
    {
        D = tdata->work->D = 30;
    }

    for (j = 0, i = 0; i < 2 * D; i++)
    {
        if (spGCD(i, 2 * D) == 1)
        {
            j++;
        }
    }

    tdata->work->R = j + 3;

    // decide on the stage 2 parameters.  Larger U means
    // more memory and more setup overhead, but more prime pairs.
    // Smaller U means the opposite.  find a good balance.
    // these pairing ratios are estimates based on Montgomery's
    // Pair algorithm Table w assuming w = 1155, for various umax.
    //static double pairing[8] = { 0.7,  0.6446, 0.6043, 0.5794, 0.5535, 0.5401, 0.5266, 0.5015};
    // these ratios are based on observations with larger (more realistic) B1/B2, 
    // specifically, B2=1e9
    static double pairing[8] = { 0.8,  0.72, 0.67, 0.63, 0.59, 0.57, 0.55, 0.54 };
    int U[8] = { 1, 2, 3, 4, 6, 8, 12, 16 };
    double adds[8];
    int numadds;
    int numinv;
    double addcost, invcost;
    double best = 99999999999.;
    int bestU = 4;


#ifdef TARGET_KNL
    // for smaller B1, lower U can be better because initialization time
    // is more significant.  At B1=3M or above the max U is probably best.
    // If memory is an issue beyond the timings, then U will have to
    // be smaller.  (e.g., if we do a -maxmem option).
    for (i = 1; i < 6; i++)
#else
    for (i = 1; i < 8; i++)
#endif
    {
        double paircost;

        // setup cost of Pb.
        // should maybe take memory into account too, since this array can
        // get pretty big (especially with multiple threads)?
        adds[i] = (double)D * (double)U[i];;

        // runtime cost of stage 2
        numadds = (B2 - B1) / (D);

        // total adds cost
        addcost = 6.0 * ((double)numadds + adds[i]);

        // inversion cost - we invert in batches of 2U at a time, taking 3*U muls.
        // estimating the cost of the inversion as equal to an add, but we do the 
        // 8 vector elements one at a time.
        numinv = ((double)numadds / (double)U[i] / 2.0) + 2;
        invcost = (double)numinv * (VECLEN * 6.0) + (double)numinv * 3.0;

        //// estimate number of prime pairs times 1 (1 mul per pair with batch inversion)
        paircost = ((double)B2 / log((double)B2) -
        	(double)B1 / log((double)B1)) * pairing[i] * 1.0;
        //
        //printf("estimating %u primes to pair\n", (uint32_t)((double)STAGE2_MAX / log((double)STAGE2_MAX) -
        //	(double)STAGE1_MAX / log((double)STAGE1_MAX)));
        //printf("%d stg2 adds + %d setup adds + %d inversions\n", numadds, (int)adds[i], numinv);
        //printf("addcost = %f\n", addcost);
        //printf("paircost = %f\n", paircost);
        //printf("invcost = %f\n", invcost);
        //printf("totalcost = %f\n", addcost + paircost + invcost);

        if ((addcost + paircost + invcost) < best)
        {
            best = addcost + paircost + invcost;
            bestU = U[i];
        }
    }

    tdata->work->U = bestU;
    tdata->work->L = tdata->work->U * 2;

    tdata->work->NWORDS = mdata->NWORDS;
    tdata->work->MAXBITS = mdata->MAXBITS;
    tdata->work->NBLOCKS = mdata->NBLOCKS;
    vec_ecm_work_init(tdata->work);
    vec_ecm_pt_init(tdata->P, mdata->NWORDS);

    // allocate and then copy some constants over to this thread's mdata structure.
    tdata->mdata = vec_monty_alloc(mdata->NWORDS);
    tdata->mdata->NWORDS = mdata->NWORDS;
    tdata->mdata->MAXBITS = mdata->MAXBITS;
    tdata->mdata->NBLOCKS = mdata->NBLOCKS;
    mpz_set(tdata->mdata->nhat, mdata->nhat);
    mpz_set(tdata->mdata->rhat, mdata->rhat);
    vecCopy(mdata->n, tdata->mdata->n);
    vecCopy(mdata->n, tdata->work->n);
    vecCopy(mdata->one, tdata->mdata->one);
    vecCopy(mdata->vnhat, tdata->mdata->vnhat);
    vecCopy(mdata->vrhat, tdata->mdata->vrhat);
    tdata->mdata->nbits = mdata->nbits;
    tdata->mdata->use_vnbits = 0;             // all lanes use the same N.
    tdata->mdata->isMersenne = mdata->isMersenne;
    memcpy(tdata->mdata->vrho, mdata->vrho, VECLEN * sizeof(base_t));

    return;
}
