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

#include "yafu.h"
#include "avx_ecm.h"
#include "soe.h"
#include "queue.h"
#include "gmp.h"
#include "util.h"
#include "nfs.h"

// performance comparison
// http://www.mersenneforum.org/showthread.php?t=16480&page=20
// http://www.mersenneforum.org/showthread.php?t=5722&page=122

// t-level estimate
// http://mersenneforum.org/showpost.php?p=427989&postcount=2429


static int debugctr = 0;
void thread_init(thread_data_t* tdata, vec_monty_t* mdata);

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
    vec_dest->size = 0;
    while (mpz_cmp_ui(src_cp, 0) > 0)
    {
        base_t thisword = mpz_get_ui(src_cp) & VEC_MAXDIGIT;
        vec_dest->data[lane + i * VECLEN] = thisword;
        i++;
        mpz_tdiv_q_2exp(src_cp, src_cp, DIGITBITS);
    }

    vec_dest->size = MAX(vec_dest->size, i);
    mpz_clear(src_cp);
    return;
}




// =============== 64-bit hashing ================ //
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

factor_t* ecm_add_factor(factor_t* factors, int *numfactors, 
    mpz_t gmp_factor, int stg_id, int thread_id, int vec_id, 
    int curve_id, uint64_t sigma)
{
    if (*numfactors == 0)
    {
        factors = (factor_t*)malloc(1 * sizeof(factor_t));
    }
    else
    {
        factors = (factor_t*)realloc(factors, (*numfactors + 1) * sizeof(factor_t));
    }
    mpz_init(factors[*numfactors].factor);
    mpz_set(factors[*numfactors].factor, gmp_factor);
    factors[*numfactors].method = stg_id;
    factors[*numfactors].tid = thread_id;
    factors[*numfactors].vid = vec_id;
    factors[*numfactors].sigma = sigma;
    factors[*numfactors].curve_num = curve_id;
    (*numfactors)++;
    return factors;
}

factor_t * vec_ecm_main(mpz_t N, uint32 numcurves, uint64 B1, 
    uint64 B2, int threads, int* numfactors, int verbose, 
    int save_b1, uint32 *curves_run)
{
    thread_data_t *tdata;
	vec_bignum_t **f, *n;
    mpz_t g, r;
	uint32_t *siglist;
	uint32_t numcurves_per_thread;
	uint64_t b1 = B1;
	uint32_t i, j;
	char **nextptr;
	vec_monty_t *montyconst;
	int pid = getpid();
    uint64_t limit;
    int size_n, isMersenne = 0, forceNoMersenne = 0;
    factor_t * factors;
    
    // primes
    uint32_t seed_p[6542];
    uint32_t numSOEp;

	// timing variables
	struct timeval stopt;	// stop time of this job
	struct timeval startt;	// start time of this job
	double t_time;

	gettimeofday(&startt, NULL);

	if (pid <= 0)
		pid = startt.tv_usec;

    if (verbose > 1)
	    printf("process id is %d\n", pid);

    *numfactors = 0;
    mpz_init(g);
    mpz_init(r);

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
    }

    // if the input is Mersenne and still contains algebraic factors, remove them.
    if (isMersenne)
    {
#ifdef USE_NFS
        snfs_t poly;

        snfs_init(&poly);
        poly.form_type = SNFS_BRENT;
        mpz_set_ui(poly.base1, 2);
        poly.coeff1 = 1;
        poly.coeff2 = -1;
        poly.exp1 = size_n;
        find_primitive_factor(&poly, verbose);

        mpz_tdiv_q(g, N, poly.primitive);
        mpz_gcd(N, N, poly.primitive);
        factors = ecm_add_factor(factors, numfactors, g,
            0, 0, -1, -1, 0);

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
        MAXBITS = 208;
        while (MAXBITS <= mpz_sizeinbase(N, 2))
        {
            MAXBITS += 208;
        }
    }
    else
    {
        MAXBITS = 128;
        while (MAXBITS <= mpz_sizeinbase(N, 2))
        {
            MAXBITS += 128;
        }
    }

    NWORDS = MAXBITS / DIGITBITS;
    NBLOCKS = NWORDS / BLOCKWORDS;

    // and compute NBLOCKS if using Mersenne mod
    if (DIGITBITS == 52)
    {
        MAXBITS = 208;
        while (MAXBITS <= size_n)
        {
            MAXBITS += 208;
        }
    }
    else
    {
        MAXBITS = 128;
        while (MAXBITS <= size_n)
        {
            MAXBITS += 128;
        }
    }

    if (verbose > 0)
        gmp_printf("commencing parallel ecm on %Zd\n", N);

    if ((double)NWORDS / ((double)MAXBITS / (double)DIGITBITS) < 0.7)
    {
        if (verbose > 0)
            printf("Mersenne input 2^%d - 1 determined to be faster by REDC\n", size_n);
        forceNoMersenne = 1;
    }
    else
    {
        NWORDS = MAXBITS / DIGITBITS;
        NBLOCKS = NWORDS / BLOCKWORDS;
    }

    if (forceNoMersenne)
    {
        isMersenne = 0;
        size_n = mpz_sizeinbase(N, 2);
    }

    STAGE1_MAX = b1;
    STAGE2_MAX = 100ULL * (uint64_t)b1;

    // now that we know NWORDS, can allocate the monty structure.
    montyconst = vec_monty_alloc();
    tdata = (thread_data_t*)malloc(threads * sizeof(thread_data_t));

    if (isMersenne)
    {
        // TBD: may want to check if the input number is small
        // enough compared to the base Mersenne that normal Montgomery
        // arithmetic would be faster.
        montyconst->isMersenne = 1;
        montyconst->nbits = size_n;
        mpz_set(montyconst->nhat, N);           // remember input N
        // do all math w.r.t the Mersenne number
        mpz_set_ui(N, 1);
        mpz_mul_2exp(N, N, size_n);
        mpz_sub_ui(N, N, 1);
        broadcast_mpz_to_vec(montyconst->n, N);

        mpz_set_ui(r, 1);
        mpz_mul_2exp(r, r, DIGITBITS * NWORDS);
        //gmp_printf("r = (1 << %d) = %Zd\n", DIGITBITS * NWORDS, r);
        //mpz_invert(montyconst->nhat, N, r);
        //mpz_sub(montyconst->nhat, r, montyconst->nhat);
        mpz_invert(montyconst->rhat, r, N);
        broadcast_mpz_to_vec(montyconst->n, N);
        broadcast_mpz_to_vec(montyconst->r, r);
        broadcast_mpz_to_vec(montyconst->vrhat, montyconst->rhat);
        broadcast_mpz_to_vec(montyconst->vnhat, montyconst->nhat);
        //mpz_tdiv_r(r, r, N);
        mpz_set_ui(r, 1);
        broadcast_mpz_to_vec(montyconst->one, r);
    }
    else
    {
        montyconst->isMersenne = 0;
        montyconst->nbits = mpz_sizeinbase(N, 2);
        mpz_set_ui(r, 1);
        mpz_mul_2exp(r, r, DIGITBITS * NWORDS);
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

    if (verbose > 1)
        printf("ECM has been configured with DIGITBITS = %u, VECLEN = %d, GMP_LIMB_BITS = %d\n",
            DIGITBITS, VECLEN, GMP_LIMB_BITS);

    if (verbose > 1)
        printf("Choosing MAXBITS = %u, NWORDS = %d, NBLOCKS = %d based on input size %d\n",
            MAXBITS, NWORDS, NBLOCKS, size_n);

    DO_STAGE2 = 1;
    if (B2 == (uint64_t)-1)
	{
        STAGE2_MAX = STAGE1_MAX;
		DO_STAGE2 = 0;
	}

    if (STAGE1_MAX < 1000)
    {
        printf("stage 1 too small, resetting to 1000\n");
        STAGE1_MAX = 1000;
    }

    numSOEp = tiny_soe(65537, seed_p);
    if ((PRIMES[0] != 2) || PRIMES[NUM_P - 1] != 99999989)
    {
        PRIMES = soe_wrapper(seed_p, numSOEp, 0, 100000000, 0, &NUM_P);
        P_MAX = PRIMES[NUM_P - 1];
        P_MIN = PRIMES[0];
        PRIME_RANGE = 100000000;

        if (verbose > 1)
            printf("cached %u primes < %u\n", NUM_P, PRIMES[NUM_P - 1]);
    }

	if (numcurves < threads)
		numcurves = threads;
	
	numcurves_per_thread = numcurves / threads + (numcurves % threads != 0);
	numcurves = numcurves_per_thread * threads;

    if (verbose > 1)
    {
        printf("Input has %d bits, using %d threads (%d curves/thread)\n",
            mpz_sizeinbase(N, 2), threads, numcurves_per_thread);
        printf("Processing in batches of %u primes\n", PRIME_RANGE);
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
    
    if (DIGITBITS == 52)
    {
        if (montyconst->isMersenne)
        {
            vecmulmod_ptr = &vecmulmod52_mersenne;
            vecsqrmod_ptr = &vecsqrmod52_mersenne;
            vecaddmod_ptr = &vecaddmod52_mersenne;
            vecaddsubmod_ptr = &vec_simul_addsub52_mersenne;
            vecsubmod_ptr = &vecsubmod52_mersenne;
            printf("Using special Mersenne mod for factor of: 2^%d-1\n", montyconst->nbits);
        }
        else
        {
            vecmulmod_ptr = &vecmulmod52;
            vecsqrmod_ptr = &vecsqrmod52;
            vecaddmod_ptr = &vecaddmod52;
            vecaddsubmod_ptr = &vec_simul_addsub52;
            vecsubmod_ptr = &vecsubmod52;
        }
    }
    else
    {
        if (montyconst->isMersenne)
        {
            vecmulmod_ptr = &vecmulmod_mersenne;
            vecsqrmod_ptr = &vecsqrmod_mersenne;
            printf("Using special Mersenne mod for factor of: 2^%d-1\n", montyconst->nbits);
        }
        else
        {
            vecmulmod_ptr = &vecmulmod;
            vecsqrmod_ptr = &vecsqrmod;
        }

        vecaddmod_ptr = &vecaddmod;
        vecsubmod_ptr = &vecsubmod;
        vecaddsubmod_ptr = &vec_simul_addsub;
    }

	gettimeofday(&stopt, NULL);

    uint32_t seed1, seed2;
    get_random_seeds(&seed1, &seed2);
    for (i = 0; i < threads; i++)
    {
        int j;
        thread_init(&tdata[i], montyconst);
        tdata[i].curves = numcurves_per_thread;
        tdata[i].tid = i;
        tdata[i].numfactors = 0;
        // each thread gets a random LCG state that folds in the current pid,
        // this session's rand_t, and the current microsecond timer.  Plus the
        // thread id to make it unique per-thread.
		tdata[i].lcg_state = hash64(stopt.tv_usec) + hash64(pid) + hash64(i+1) + 
            seed1 + (uint64)seed2 << 32;
        tdata[i].total_threads = threads;
        tdata[i].verbose = verbose;
        tdata[i].save_b1 = save_b1;
    }

	gettimeofday(&stopt, NULL);
    t_time = yafu_difftime(&startt, &stopt);

    if (verbose > 1)
        printf("Initialization took %1.4f seconds.\n", t_time);

    // top level ECM 
    vececm(tdata);    

    if (verbose > 1)
	    printf("\n");

    // clean up thread data
	for (i = 0; i < threads; i++)
	{
        //printf("thread %d found %d factors\n", i, tdata[i].numfactors);
        for (j = 0; j < tdata[i].numfactors; j++)
        {
            factors = ecm_add_factor(factors, numfactors, tdata[i].factors[j].factor,
                tdata[i].factors[j].stg_id, tdata[i].factors[j].thread_id,
                tdata[i].factors[j].vec_id, tdata[i].factors[j].curve_id, 
                tdata[i].factors[j].sigma);
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

    *curves_run = tdata[0].curves;

    // clean up local/global data
    mpz_clear(g);
    mpz_clear(r);
    vec_monty_free(montyconst);
	free(montyconst);
    free(tdata);

	return factors;
}

void thread_init(thread_data_t *tdata, vec_monty_t *mdata)
{    
    //mpz_init(tdata->factor);
    tdata->work = (ecm_work *)malloc(sizeof(ecm_work));
    tdata->P = (ecm_pt *)malloc(sizeof(ecm_pt));
    tdata->sigma = (uint64_t *)malloc(VECLEN * sizeof(uint64_t));
	uint32_t D = tdata->work->D = 1155;
	tdata->work->R = 480 + 3;

	// decide on the stage 2 parameters.  Larger U means
	// more memory and more setup overhead, but more prime pairs.
	// Smaller U means the opposite.  find a good balance.
	static double pairing[4] = { 0.7283, 0.6446, 0.5794, 0.5401 };
	static double adds[4];
	double best = 99999999999.;
	int bestU = 4;

	if (STAGE1_MAX <= 8192)
	{
		D = tdata->work->D = 210;
		tdata->work->R = 48 + 3;
	}

	adds[0] = (double)D * 1.0;
	adds[1] = (double)D * 2.0;
	adds[2] = (double)D * 4.0;
	adds[3] = (double)D * 8.0;

	int i;
	for (i = 1; i < 4; i++)
	{
		int numadds = (STAGE2_MAX - STAGE1_MAX) / (2 * D);
		double addcost = 6.0 * ((double)numadds + adds[i]);
		double paircost = ((double)STAGE2_MAX / log((double)STAGE2_MAX) -
			(double)STAGE1_MAX / log((double)STAGE1_MAX)) * pairing[i] * 2.0;

		//printf("estimating %u primes paired\n", (uint32_t)((double)STAGE2_MAX / log((double)STAGE2_MAX) -
		//	(double)STAGE1_MAX / log((double)STAGE1_MAX)));
		//printf("%d adds + %d setup adds\n", numadds, (int)adds[i]);
		//printf("addcost = %f\n", addcost);
		//printf("paircost = %f\n", paircost);
		//printf("totalcost = %f\n", addcost + paircost);

		if ((addcost + paircost) < best)
		{
			best = addcost + paircost;
			bestU = 1 << i;
		}
	}

	tdata->work->U = bestU;
	tdata->work->L = bestU * 2;

    vec_ecm_work_init(tdata->work);
    vec_ecm_pt_init(tdata->P);

    // allocate and then copy some constants over to this thread's mdata structure.
    tdata->mdata = vec_monty_alloc();
    mpz_set(tdata->mdata->nhat, mdata->nhat);
    mpz_set(tdata->mdata->rhat, mdata->rhat);
    vecCopy(mdata->n, tdata->mdata->n);
    vecCopy(mdata->n, tdata->work->n);
    vecCopy(mdata->one, tdata->mdata->one);
    vecCopy(mdata->vnhat, tdata->mdata->vnhat);
    vecCopy(mdata->vrhat, tdata->mdata->vrhat);
    tdata->mdata->nbits = mdata->nbits;
    tdata->mdata->isMersenne = mdata->isMersenne;
    memcpy(tdata->mdata->vrho, mdata->vrho, VECLEN * sizeof(base_t));

    return;
}
