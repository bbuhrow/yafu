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

#include <inttypes.h>
#include "yafu_ecm.h"
#include "avx_ecm.h"
#include "avxecm_stage2fft.h"
#include "soe.h"
#include "gmp.h"
#include "ytools.h"
#include "nfs_impl.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef __MINGW32__
#include <sys/time.h>
#endif

// performance comparison
// http://www.mersenneforum.org/showthread.php?t=16480&page=20
// http://www.mersenneforum.org/showthread.php?t=5722&page=122

// t-level estimate
// http://mersenneforum.org/showpost.php?p=427989&postcount=2429


static int debugctr = 0;
void thread_init(thread_data_t* tdata, vec_monty_t* mdata, uint64_t B1, uint64_t B2);
void test_vecarith(int numbits, int iterations, int validate);
void vec_bench(int lo_words, int hi_words, int step_words, int reps, int validate);

/* Toom-3 supersedes Karatsuba from the Karatsuba band entry (~16400 bits =
   ~317 words; pointwise products ~107 words vs Karatsuba's ~159-word halves).
   Initial guess only -- the planned benchmark harness will establish the true
   Toom/Karatsuba crossover; lower this constant (keeping it > 16400) to move
   the boundary down inside the Karatsuba band.  (plans/toom3_finish.md) */
#define TOOM3_MIN_BITS 26400

void assign_mul_ptrs(vec_monty_t* montyconst, int maxbits, int verbose)
{
    if (montyconst->rns != NULL)
    {
        // RNS backend (avxecm_rns.c): arith/rns batch kernels
        int rn, rs;
        avxecm_rns_fixed_info(avxecm_rns_get_fixed(montyconst), &rn, &rs);
        vecmulmod_ptr = &avxecm_rns_mul;
        vecsqrmod_ptr = &avxecm_rns_sqr;
        vecaddmod_ptr = &avxecm_rns_add;
        vecsubmod_ptr = &avxecm_rns_sub;
        vecaddsubmod_ptr = &avxecm_rns_addsub;
        if (verbose > 1)
        {
            printf("Using RNS mul/sqr core (2 x %d channels of %d words)\n", rn, rs);
        }
        return;
    }

#if defined(USE_AVX2) && !defined(USE_AVX512F)
    // AVX2 build (no AVX-512): the AVX2 template cores at every size (loop
    // form above 32 words) and the AVX2 add/sub (vecarith52_avx2.c); the
    // other cores and the Mersenne forms are AVX-512 only
    if (!vecmulmod52_tmpl_get((int)montyconst->NWORDS, &vecmulmod_ptr, &vecsqrmod_ptr))
    {
        printf("no AVX2 mul/sqr core for %u words\n", montyconst->NWORDS);
        exit(1);
    }
    vecaddmod_ptr = &vecaddmod52_avx2;
    vecsubmod_ptr = &vecsubmod52_avx2;
    vecaddsubmod_ptr = &vec_simul_addsub52_avx2;
    if (verbose > 1)
    {
        printf("Using AVX2 constant-S template mul/sqr core (%u words)\n", montyconst->NWORDS);
    }
    return;
#endif

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
            // initial tests show the crossover to karatsuba is at
            // around 7k bits, but there are issues still with finding
            // factors.
            vecmulmod_ptr = &vecmulmod52_mersenne;
            vecsqrmod_ptr = &vecsqrmod52_mersenne;
            
            //if (maxbits > 1000)
            //{
            //    // faster starting around 16000 bits on Xeon 6254
            //    vecmulmod_ptr = &veckmul_mersenne;
            //    vecsqrmod_ptr = &vecksqr_mersenne;
            //    if (verbose > 1)
            //    {
            //        printf("Using karatsuba mul/sqr core\n");
            //    }
            //}

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
#if defined(USE_AVX512F) && !defined(AVXECM_TMPL_KERNELS)
#error "AVX-ECM needs the constant-S template cores (the old fixed / BFIPS cores are no longer built)"
#endif
            // constant-S template cores (arith/vecmont/vecarith52_tmpl.c) up to
            // AVXECM_KARATSUBA_BITS, karatsuba REDC above (inputs too large for
            // the RNS backend).  The old fixed416..1040 / BFIPS cores are in
            // vecarith52_oldcores.c, outside the builds.
            if (maxbits > AVXECM_KARATSUBA_BITS)
            {
                vecmulmod_ptr = &veckmul_redc;
                vecsqrmod_ptr = &vecksqr_redc;
                if (verbose > 1)
                {
                    printf("Using karatsuba mul/sqr core\n");
                }
            }
            else
            {
#ifdef AVXECM_TMPL_KERNELS
                if (!vecmulmod52_tmpl_get(maxbits / DIGITBITS, &vecmulmod_ptr, &vecsqrmod_ptr))
                {
                    printf("no mul/sqr core for %d words\n", maxbits / DIGITBITS);
                    exit(1);
                }
                if (verbose > 1)
                {
                    printf("Using constant-S template mul/sqr core (%d words)\n",
                        maxbits / DIGITBITS);
                }
#else
                printf("AVX-ECM needs an AVX-512 or AVX2 build\n");
                exit(1);
#endif
            }
            vecaddsubmod_ptr = &vec_simul_addsub52;
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
}

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

    // clear this lane's remaining words: callers reuse vectors (stage-2
    // batch inversion, curve building), and a shorter value must not keep
    // stale high words from the previous one.
    for (; i < (int)vec_dest->WORDS_ALLOC; i++)
    {
        vec_dest->data[lane + i * VECLEN] = 0;
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
    int size_n, isMersenne = 0, forceNoMersenne = 1;
    uint64_t sigma = fobj->ecm_obj.sigma;
    uint32_t maxbits, mmaxbits, nwords, nblocks;

	// timing variables
	struct timeval stopt;	// stop time of this job
	struct timeval startt;	// start time of this job
	double t_time;

    // Kara / Toom / School benchmarking, but using arith functions that no longer exist.
#if 0
	/* ---- microbenchmark mode (plans/toom3_bench_tune.md, Phase 1) ----
	   YAFU_VECBENCH=lo,hi,step,reps[,validate]
	   Widths in 52-bit words (step kept a multiple of 4).  For each width
	   and each REDC core (base, kara, toom) times tight mul/sqr/addsub
	   loops and prints one CSV row; then exits.  No ECM is run.
	   Example:  YAFU_VECBENCH=320,2048,32,3,1 ./yafu "ecm(<big N>)"  */
	if (getenv("YAFU_VECBENCH") != NULL)
	{
		int vblo = 320, vbhi = 2048, vbstep = 32, vbreps = 3, vbvalidate = 1;
		const char* s = getenv("YAFU_VECBENCH");
		char* endp;
		long v;
		vblo = (int)strtol(s, &endp, 10);
		if (*endp == ',') { v = strtol(endp + 1, &endp, 10); vbhi = (int)v; }
		if (*endp == ',') { v = strtol(endp + 1, &endp, 10); vbstep = (int)v; }
		if (*endp == ',') { v = strtol(endp + 1, &endp, 10); vbreps = (int)v; }
		if (*endp == ',') { v = strtol(endp + 1, &endp, 10); vbvalidate = (int)v; }
		vec_bench(vblo, vbhi, vbstep, vbreps, vbvalidate);
		exit(0);
	}
#endif

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

    // set N equal to the original input, so we can
    // detect Mersenne inputs correctly.
    mpz_set(N, fobj->N);

    // check for Mersenne inputs
    size_n = mpz_sizeinbase(N, 2);

    for (i = 31; i <= size_n; i++)
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
        if (mpz_sizeinbase(g, 2) < (DIGITBITS/2))
        {
            size_n = i;
            isMersenne = mpz_get_ui(g);
            break;
        }
    }
    //printf("found isMersenne = 2^%d %d\n", size_n, isMersenne);

    // now set N equal to the actual input, which may have had factors removed
    // by previous factoring routines.
    mpz_set(N, fobj->ecm_obj.gmp_n);

    // if the input is Mersenne and still contains algebraic factors, remove them.
    if ((abs(isMersenne) == 1) && (forceNoMersenne == 0))
    {
#ifdef USE_NFS
        snfs_t poly;

        snfs_init(&poly);
        poly.form_type = SNFS_BRENT;
        mpz_set_ui(poly.base1, 2);
        poly.coeff1 = 1;
        poly.coeff2 = -isMersenne;
        poly.exp1 = size_n;
        //find_primitive_factor(fobj, &poly, fobj->primes, fobj->num_p, verbose);
        remove_algebraic_factors(fobj, &poly, fobj->primes, fobj->num_p, verbose);

        mpz_tdiv_q(g, N, poly.primitive);
        mpz_gcd(N, N, poly.primitive);

        if (mpz_cmp_ui(g, 1) > 0)
        {
            add_to_factor_list(fobj->factors, g,
                fobj->VFLAG, fobj->NUM_WITNESSES, 0);

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
    // exact word count when the constant-S template cores are available,
    // else a multiple of 4 words; >= 3 spare bits either way.
    nwords = avxecm_redc_nwords((uint32_t)mpz_sizeinbase(N, 2));
    maxbits = nwords * DIGITBITS;
    // only the blocked (BFIPS/fixed/karatsuba) cores use nblocks; they
    // are never selected when nwords is not a multiple of BLOCKWORDS.
    nblocks = nwords / BLOCKWORDS;

    // Mersenne mod keeps the old (multiple-of-4 word) sizing; MAXBITS in
    // that case stays as before too.
    uint32_t mers_maxbits = 208;
    while (mers_maxbits <= mpz_sizeinbase(N, 2))
    {
        mers_maxbits += 208;
    }
    //printf("input has %d bits, %d words, %d blocks, mersenne size is %d\n", 
    //    mpz_sizeinbase(N, 2), nwords, nblocks, size_n);

    // and compute NBLOCKS if using Mersenne mod
    if (DIGITBITS == 52)
    {
        mmaxbits = 208;
        while (mmaxbits <= size_n)
        {
            mmaxbits += 208;
        }
    }
    else
    {
        mmaxbits = 128;
        while (mmaxbits <= size_n)
        {
            mmaxbits += 128;
        }
    }

    if (verbose > 1)
    {
        if (mpz_sizeinbase(N, 2) > 1024)
        {
            int ndigits = gmp_base10(N);
            mpz_set_ui(r, 10);
            mpz_pow_ui(r, r, ndigits - 10);
            mpz_tdiv_q(g, N, r);
            mpz_set_ui(r, 10);
            mpz_pow_ui(r, r, 10);
            mpz_tdiv_r(r, N, r);
            gmp_printf("commencing parallel ecm on %Zd...%Zd [%d digits] with %d threads\n", 
                g, r, ndigits, threads);
        }
        else
        {
            gmp_printf("commencing parallel ecm on %Zd with %d threads\n", N, threads);
        }
    }

    if ((isMersenne != 0) && (forceNoMersenne == 0) && 
        ((double)nwords / ((double)mmaxbits / (double)DIGITBITS) < 0.7))
    {
        if (verbose > 1)
        {
            printf("Mersenne input 2^%d - 1 determined to be faster by REDC\n", 
                size_n);
        }
        forceNoMersenne = 1;
        size_n = mpz_sizeinbase(N, 2);
    }
    else if ((isMersenne != 0) && (forceNoMersenne == 0))
    {
        maxbits = mers_maxbits;
        nwords = mmaxbits / DIGITBITS;
        nblocks = nwords / BLOCKWORDS;
    }
    else
    {
        size_n = mpz_sizeinbase(N, 2);
        nwords = maxbits / DIGITBITS;
        nblocks = nwords / BLOCKWORDS;
    }

    if (forceNoMersenne)
    {
        isMersenne = 0;
        size_n = mpz_sizeinbase(N, 2);
    }

    // RNS backend for large enough non-Mersenne inputs (avx_rns_bits,
    // option avxecm_rns_bits; 0 = never).  The vectors must hold the RNS
    // layout: n*s words.
    void *rnsfx = NULL;
    if ((isMersenne == 0) && (fobj->ecm_obj.avx_rns_bits > 0) &&
        (size_n >= fobj->ecm_obj.avx_rns_bits) && avxecm_rns_available())
    {
        uint32_t rnswords = 0;
        rnsfx = avxecm_rns_fixed(N, &rnswords);
        if (rnsfx != NULL)
        {
            if (rnswords > nwords)
                nwords = rnswords;
            maxbits = nwords * DIGITBITS;
            nblocks = nwords / BLOCKWORDS;
        }
        else if (verbose > 1)
        {
            printf("RNS arithmetic unavailable for this input, using radix 2^52\n");
        }
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
            printf("starting with sigma = %"PRIu64"\n", sigma);
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

        if ((rnsfx != NULL) && avxecm_rns_attach(montyconst, rnsfx, 1, N))
        {
            // values are RNS residues from here on (mdata->n stays radix
            // 2^52: it is only read back as an integer)
            vec_set_monty_one(montyconst, montyconst->one, N);
        }
        else if (rnsfx != NULL)
        {
            printf("RNS setup failed\n");
            exit(1);
        }
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

    assign_mul_ptrs(montyconst, tdata[0].MAXBITS, verbose);

	gettimeofday(&stopt, NULL);

    if (0)
    {
        size_n = mpz_sizeinbase(fobj->N, 2);
        test_vecarith(size_n, numcurves, 1);
        exit(0);

        for (size_n = 200, maxbits = 208, numcurves = 8000; size_n <= 2000; 
            size_n += 200, maxbits += 208, numcurves /= 2)
        {
            gettimeofday(&startt, NULL);
            assign_mul_ptrs(montyconst, maxbits, verbose);
            test_vecarith(size_n, numcurves, 1);
            gettimeofday(&stopt, NULL);
            t_time = ytools_difftime(&startt, &stopt);
            printf("Test took %1.4f seconds.\n", t_time);
        }

        numcurves = 1;
        for (size_n = 10000; size_n <= 100000; size_n += 4000)
        {
            gettimeofday(&startt, NULL);
            assign_mul_ptrs(montyconst, maxbits, verbose);
            test_vecarith(size_n, numcurves, 1);
            gettimeofday(&stopt, NULL);
            t_time = ytools_difftime(&startt, &stopt);
            printf("Test took %1.4f seconds.\n", t_time);
        }

        for (size_n = 200, maxbits = 208, numcurves = 8000; size_n <= 2000;
            size_n += 200, maxbits += 208, numcurves /= 2)
        {
            gettimeofday(&startt, NULL);
            assign_mul_ptrs(montyconst, maxbits, verbose);
            test_vecarith(size_n, numcurves, 0);
            gettimeofday(&stopt, NULL);
            t_time = ytools_difftime(&startt, &stopt);
            printf("Test took %1.4f seconds.\n", t_time);
        }

        numcurves = 1;
        for (size_n = 10000; size_n <= 100000; size_n += 4000)
        {
            gettimeofday(&startt, NULL);
            assign_mul_ptrs(montyconst, maxbits, verbose);
            test_vecarith(size_n, numcurves, 0);
            gettimeofday(&stopt, NULL);
            t_time = ytools_difftime(&startt, &stopt);
            printf("Test took %1.4f seconds.\n", t_time);
        }

        exit(0);
    }

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
        tdata[i].curve_type = fobj->ecm_obj.avx_curve;
        
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

    // FFT stage 2 (option avxecm_stg2): shared setup for this modulus and
    // range; the arithmetic modulus is montyconst->n (radix 2^52 in every
    // backend, the base Mersenne number for Mersenne inputs)
    {
        struct avxecm_s2fft_s *s2 = NULL;
        if (tdata[0].DO_STAGE2 && fobj->ecm_obj.avx_stg2_fft &&
            (tdata[0].STAGE2_MAX > tdata[0].STAGE1_MAX))
        {
            mpz_t M;
            mpz_init(M);
            extract_bignum_from_vec_to_mpz(M, montyconst->n, 0, (int)montyconst->NWORDS);
            // avx_stg2_fft = 2 and no explicit B2: B2 from the memory limit;
            // the Brent-Suyama polynomial is then chosen too unless given
            int auto_b2 = (fobj->ecm_obj.avx_stg2_fft == 2) && fobj->ecm_obj.stg2_is_default;
            int bs = fobj->ecm_obj.avx_stg2_bs;
            if (bs == AVXECM_BS_DEFAULT)
                bs = auto_b2 ? AVXECM_BS_AUTO : 0;
            s2 = avxecm_s2fft_setup(M, tdata[0].STAGE1_MAX, tdata[0].STAGE2_MAX,
                (int)threads, (double)fobj->ecm_obj.avx_stg2_mem, auto_b2, bs,
                fobj->ecm_obj.avx_curve == AVXECM_CURVE_EDW4X, verbose);
            mpz_clear(M);
            if (s2 != NULL)
            {
                B2 = avxecm_s2fft_B2(s2);
                for (i = 0; i < threads; i++)
                    tdata[i].STAGE2_MAX = B2;
            }
            else if (verbose > 0)
            {
                printf("ecm: FFT stage 2 unavailable for these parameters, using pairing\n");
            }
        }
        for (i = 0; i < threads; i++)
            tdata[i].s2fft = s2;
    }

	gettimeofday(&stopt, NULL);
    t_time = ytools_difftime(&startt, &stopt);

    if ((verbose > 0) && (fobj->ecm_obj.avx_curve != AVXECM_CURVE_SUYAMA))
    {
        printf("ecm: twisted Edwards curve [sigma](5,8) on a = -1 family (iv)\n"
            "ecm: using %s arithmetic in stage 1\n",
            (fobj->ecm_obj.avx_curve == AVXECM_CURVE_EDW4X) ?
            "Edwards extended-coordinate" : "X-coord-only Montgomery differential");
    }

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
                fobj->VFLAG, fobj->NUM_WITNESSES, 0);

            total_factors++;
            //fid = fobj->factors->num_factors - 1;

            fobj->factors->factors[fid].tid = tdata[i].factors[j].thread_id;
            fobj->factors->factors[fid].curve_num = tdata[i].factors[j].curve_id;
            fobj->factors->factors[fid].sigma = tdata[i].factors[j].sigma;
            fobj->factors->factors[fid].vid = tdata[i].factors[j].vec_id;
            fobj->factors->factors[fid].method = tdata[i].factors[j].stg_id;

            if (is_mpz_prp(g, fobj->NUM_WITNESSES))
            {
                if (fobj->VFLAG >= 0)
                    gmp_printf("\necm: found prp%d factor = %Zd\n",
                        gmp_base10(g), g);

                if (fobj->LOGFLAG && (strcmp(fobj->flogname, "") != 0))
                {
                    char* s = mpz_get_str(NULL, 10, g);
                    logprint(flog, "prp%d = %s (curve=%d stg=%d B1=%u B2=%"PRIu64" sigma=%"PRIu64" thread=%d vecpos=%d)\n",
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
                if (fobj->VFLAG >= 0)
                    gmp_printf("\necm: found c%d factor = %Zd\n",
                        gmp_base10(g), g);

                if (fobj->LOGFLAG && (strcmp(fobj->flogname, "") != 0))
                {
                    char* s = mpz_get_str(NULL, 10, tdata[i].factors[j].factor);
                    logprint(flog, "c%d = %s (curve=%d stg=%d B1=%u B2=%"PRIu64" sigma=%"PRIu64" thread=%d vecpos=%d)\n",
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
        avxecm_rns_detach(tdata[i].mdata);
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
                fobj->VFLAG, fobj->NUM_WITNESSES, 0);

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
    avxecm_rns_detach(montyconst);      // after the threads: owns the fixed tables
    vec_monty_free(montyconst);
	free(montyconst);
    avxecm_s2fft_free(tdata[0].s2fft);
    free(tdata);

	return;
}

void thread_init(thread_data_t *tdata, vec_monty_t *mdata, uint64_t B1, uint64_t B2)
{    
    int i, j;

    tdata->work = (ecm_work *)xmalloc(sizeof(ecm_work));
    tdata->P = (ecm_pt *)xmalloc(sizeof(ecm_pt));
    tdata->sigma = (uint64_t *)xmalloc(VECLEN * sizeof(uint64_t));
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

    if (mdata->rns != NULL)
    {
        // RNS: shared fixed tables, a context (and scratch) per thread
        mpz_t N;
        mpz_init(N);
        extract_bignum_from_vec_to_mpz(N, mdata->n, 0, (int)mdata->NWORDS);
        if (!avxecm_rns_attach(tdata->mdata, avxecm_rns_get_fixed(mdata), 0, N))
        {
            printf("RNS thread setup failed\n");
            exit(1);
        }
        mpz_clear(N);
    }

    return;
}

void test_vecarith(int numbits, int iterations, int validate)
{
    int numwords;
    vec_monty_t* mdata;
    mpz_t n, p, r, a, e, t1, t2, t3, t4;
    int i, j, numsuccess = 0;
    gmp_randstate_t state;
    vec_bignum_t* vec_a;
    vec_bignum_t* vec_p;
    vec_bignum_t* vec_n;
    vec_bignum_t* vec_e;
    vec_bignum_t* vec_s;
    vec_bignum_t* vec_b;
    vec_bignum_t* vec_one;
    vec_bignum_t* vec_am;

    numwords = (int)avxecm_redc_nwords((uint32_t)numbits);

    printf("commencing test with %d bits (%d words)\n", numbits, numwords);
    mdata = vec_monty_alloc(numwords);

    mdata->NBLOCKS = numwords / BLOCKWORDS;
    mdata->MAXBITS = numwords * DIGITBITS;
    mdata->NWORDS = numwords;

    mpz_init(n);
    mpz_init(r);
    mpz_init(a);
    mpz_init(e);
    mpz_init(p);
    mpz_init(t1);
    mpz_init(t2);
    mpz_init(t3);
    mpz_init(t4);
    vec_n = vecInit(numwords);
    vec_p = vecInit(numwords);
    vec_a = vecInit(numwords);
    vec_e = vecInit(numwords);
    vec_s = vecInit(numwords);
    vec_b = vecInit(numwords);
    vec_am = vecInit(numwords);
    vec_one = vecInit(numwords);
    gmp_randinit_default(state);

    vecClear(vec_one);
    for (j = 0; j < 8; j++)
    {
        vec_one->data[j] = 1;
    }

    for (i = 0; i < iterations; i++)
    {
        // get a random odd modulus
        mpz_urandomb(n, state, numbits);
        if (mpz_even_p(n))
        {
            mpz_add_ui(n, n, 1);
        }

        //gmp_printf("n = %Zd\n", n);
        vecClear(vec_n);
        broadcast_mpz_to_vec(vec_n, n);

        // get a random exponent
        mpz_urandomb(e, state, numbits);
        vecClear(vec_e);
        broadcast_mpz_to_vec(vec_e, e);

        printf("test %d of %d\r", i, iterations);

        // initialize monty arithmetic
        mdata->isMersenne = 0;
        mdata->nbits = mpz_sizeinbase(n, 2);
        mdata->use_vnbits = 0;             // all lanes use the same N.
        mpz_set_ui(r, 1);
        mpz_mul_2exp(r, r, numwords * DIGITBITS);
        mpz_invert(mdata->nhat, n, r);
        mpz_sub(mdata->nhat, r, mdata->nhat);
        //mpz_invert(mdata->rhat, r, n);
        mpz_mul(mdata->rhat, r, r);
        mpz_tdiv_r(mdata->rhat, mdata->rhat, n);
        vecClear(mdata->n);
        vecClear(mdata->r);
        vecClear(mdata->vrhat);
        vecClear(mdata->vnhat);
        vecClear(mdata->one);
        broadcast_mpz_to_vec(mdata->n, n);
        broadcast_mpz_to_vec(mdata->r, r);
        broadcast_mpz_to_vec(mdata->vrhat, mdata->rhat);
        broadcast_mpz_to_vec(mdata->vnhat, mdata->nhat);
        mpz_tdiv_r(mdata->rhat, r, n);
        broadcast_mpz_to_vec(mdata->one, mdata->rhat);

        for (j = 0; j < 8; j++)
        {
            mdata->vrho[j] = mpz_get_ui(mdata->nhat) & VEC_MAXDIGIT;
        }

        vecClear(vec_a);
        vecClear(vec_b);
        vecClear(vec_am);
        vecClear(vec_p);

        // get 8 random numbers and load them into a vector
        for (j = 0; j < 8; j++)
        {
            mpz_urandomb(a, state, numbits);
            mpz_tdiv_r(a, a, n);
            //if (j == 0) gmp_printf("a = %Zd\ne = %Zd\nn = %Zd\n", a, e, n);
            insert_mpz_to_vec(vec_a, a, j);
        }

        // back up input base
        vecCopyn(vec_a, vec_b, numwords);

        // monty rep
        vecmulmod_ptr(vec_a, mdata->vrhat, vec_am, vec_n, vec_s, mdata);        

        extract_bignum_from_vec_to_mpz(t1, vec_a, 0, numwords);
        extract_bignum_from_vec_to_mpz(t2, mdata->vrhat, 0, numwords);
        {

            // arguments are:
            // 1: T - input to reduce
            // 2: n - reduction modulus
            // 3: r - bits in R
            mpz_mul(t4, t1, t2);

            mpz_set_ui(t1, 1);
            mpz_mul_2exp(t1, t1, numwords * DIGITBITS);
            mpz_invert(t2, n, t1);
            mpz_sub(t2, t1, t2);
            mpz_mul(t3, t4, t2);
            mpz_tdiv_r_2exp(t3, t3, numwords * DIGITBITS);
            mpz_mul(t3, t3, n);
            mpz_add(t3, t3, t4);
            mpz_tdiv_q_2exp(t4, t3, numwords * DIGITBITS);
            if (mpz_cmp(t4, n) >= 0)
            {
                mpz_sub(t4, t4, n);
            }
        }
        extract_bignum_from_vec_to_mpz(t3, vec_am, 0, numwords);
        if (mpz_cmp(t3, t4) != 0)
        {
            printf("error during monty convert\n");
            extract_bignum_from_vec_to_mpz(t1, vec_a, 0, numwords);
            extract_bignum_from_vec_to_mpz(t2, mdata->vrhat, 0, numwords);
            gmp_printf("in1 : %Zd\n", t1);
            gmp_printf("in2 : %Zd\n", t2);
            gmp_printf("n   : %Zd\n", n);
            gmp_printf("out1: %Zd\n", t3);
            gmp_printf("out2: %Zd\n", t4);
        }

        // vec powm: tests mul and sqr
        // TODO: NEEDS THIS FROM THE VECPM1 CODE
        //vecslidingmodexp(vec_p, vec_am, e, vec_n, vec_s, mdata->one, mdata);
        //vecmodexp(vec_p, vec_am, e, vec_n, vec_s, mdata->one, mdata);

        // un-monty rep
        vecmulmod_ptr(vec_p, vec_one, vec_a, vec_n, vec_s, mdata);      

        // final reduction mod N
        uint32_t msk = vec_gte52(vec_p, vec_n);
        vec_bignum52_mask_sub(vec_p, vec_n, vec_p, msk);

        extract_bignum_from_vec_to_mpz(t1, vec_p, 0, numwords);
        extract_bignum_from_vec_to_mpz(t2, vec_one, 0, numwords);
        {

            // arguments are:
            // 1: T - input to reduce
            // 2: n - reduction modulus
            // 3: r - bits in R
            mpz_mul(t4, t1, t2);

            mpz_set_ui(t1, 1);
            mpz_mul_2exp(t1, t1, numwords * DIGITBITS);
            mpz_invert(t2, n, t1);
            mpz_sub(t2, t1, t2);
            mpz_mul(t3, t4, t2);
            mpz_tdiv_r_2exp(t3, t3, numwords * DIGITBITS);
            mpz_mul(t3, t3, n);
            mpz_add(t3, t3, t4);
            mpz_tdiv_q_2exp(t4, t3, numwords * DIGITBITS);
            if (mpz_cmp(t4, n) >= 0)
            {
                mpz_sub(t4, t4, n);
            }
        }
        extract_bignum_from_vec_to_mpz(t3, vec_a, 0, numwords);
        if (mpz_cmp(t3, t4) != 0)
        {
            printf("error during un-monty convert\n");
            gmp_printf("in1 : %Zd\n", t1);
            gmp_printf("in2 : %Zd\n", t2);
            gmp_printf("out1: %Zd\n", t3);
            gmp_printf("out2: %Zd\n", t4);
        }

        // compute reference powm and compare
        if (validate)
        {
            int failed = 0;
            for (j = 0; j < 8; j++)
            {
                extract_bignum_from_vec_to_mpz(a, vec_b, j, numwords);
                mpz_powm(p, a, e, n);
                extract_bignum_from_vec_to_mpz(a, vec_a, j, numwords);
                if (mpz_cmp(p, a) == 0)
                {
                    numsuccess++;
                }
                else
                {
                    extract_bignum_from_vec_to_mpz(a, vec_b, j, numwords);
                    gmp_printf("test failed in lane %d!\na = %Zd\ne = %Zd\nn = %Zd\nr = %Zd\n", 
                        j, a, e, n, r);
                    extract_bignum_from_vec_to_mpz(a, vec_a, j, numwords);
                    gmp_printf("gmp = %Zd\nvec = %Zd\n", p, a);
                    failed = 1;
                }
            }
            if (failed)
            {
                exit(1);
            }
        }

    }

    if (validate)
    {
        printf("successes: %d of %d\n", numsuccess, iterations * VECLEN);
    }

    mpz_clear(n);
    mpz_clear(r);
    mpz_clear(a);
    mpz_clear(e);
    mpz_clear(p);
    mpz_clear(t1);
    mpz_clear(t2);
    mpz_clear(t3);
    mpz_clear(t4);
    vecFree(vec_n);
    vecFree(vec_p);
    vecFree(vec_a);
    vecFree(vec_e);
    vecFree(vec_s);
    vecFree(vec_b);
    vecFree(vec_am);
    vecFree(vec_one);
    vec_monty_free(mdata);

    return;
}

// Kara / Toom / School benchmarking, but using arith functions that no longer exist.
#if 0

/* ---- microbenchmark of the REDC cores (plans/toom3_bench_tune.md) ----
   For each operand width and each core (base, kara, toom) prints one CSV
   row:  nwords,core,mul_ns,sqr_ns,addsub_ns,iters
   Times are the minimum over `reps` runs of a ping-pong loop calibrated to
   ~0.3 s.  Operands are random values in Montgomery form; each call's
   output feeds the next call (ping-pong, no aliasing), so the loop cannot
   be optimized away and cache state is realistic.  If validate is nonzero,
   one mul/sqr/addsub per (width, core) is checked against GMP first and
   the run aborts on any mismatch. */
void vec_bench(int lo_words, int hi_words, int step_words, int reps, int validate)
{
    int numwords, core, rep, it, iters, j;
    gmp_randstate_t state;
    mpz_t n, r, pa, pb, ma, mb, t1, t2, t3, t4;
    vec_bignum_t *va, *vb, *vc, *vd, *vn, *vs;
    vec_monty_t* mdata;
    struct timeval tv0, tv1;
    const char* corenames[3] = { "base", "kara", "toom" };
    void (*mulptrs[3])(vec_bignum_t *, vec_bignum_t *, vec_bignum_t *, vec_bignum_t *, vec_bignum_t *, vec_monty_t *);
    void (*sqrptrs[3])(vec_bignum_t *, vec_bignum_t *, vec_bignum_t *, vec_bignum_t *, vec_monty_t *);

    mulptrs[0] = &vecmulmod52;   sqrptrs[0] = &vecsqrmod52;
    mulptrs[1] = &veckmul_redc;  sqrptrs[1] = &vecksqr_redc;
    mulptrs[2] = &vectmul_redc;  sqrptrs[2] = &vectsqr_redc;

    if (step_words < 4)
        step_words = 4;
    step_words -= step_words % 4;
    if (lo_words < 8)
        lo_words = 8;
    if (reps < 1)
        reps = 1;

    gmp_randinit_default(state);
    mpz_init(n);
    mpz_init(r);
    mpz_init(pa);
    mpz_init(pb);
    mpz_init(ma);
    mpz_init(mb);
    mpz_init(t1);
    mpz_init(t2);
    mpz_init(t3);
    mpz_init(t4);

    printf("# vecbench: words %d..%d step %d reps %d validate %d\n",
        lo_words, hi_words, step_words, reps, validate);
    printf("nwords,core,mul_ns,sqr_ns,addsub_ns,iters\n");

    for (numwords = lo_words; numwords <= hi_words; numwords += step_words)
    {
        mdata = vec_monty_alloc(numwords);
        mdata->NBLOCKS = numwords / BLOCKWORDS;
        mdata->MAXBITS = numwords * DIGITBITS;
        mdata->NWORDS = numwords;

        va = vecInit(numwords);
        vb = vecInit(numwords);
        vc = vecInit(numwords);
        vd = vecInit(numwords);
        vn = vecInit(numwords);
        vs = vecInit(numwords);

        /* random full-size odd modulus just under R = 2^(52*numwords) */
        mpz_urandomb(n, state, numwords * DIGITBITS);
        mpz_setbit(n, numwords * DIGITBITS - 1);
        mpz_setbit(n, 0);

        /* monty constants - same setup as test_vecarith */
        mdata->isMersenne = 0;
        mdata->nbits = mpz_sizeinbase(n, 2);
        mdata->use_vnbits = 0;
        mpz_set_ui(r, 1);
        mpz_mul_2exp(r, r, numwords * DIGITBITS);
        mpz_invert(mdata->nhat, n, r);
        mpz_sub(mdata->nhat, r, mdata->nhat);
        mpz_mul(mdata->rhat, r, r);
        mpz_tdiv_r(mdata->rhat, mdata->rhat, n);
        vecClear(mdata->n);
        vecClear(mdata->r);
        vecClear(mdata->vrhat);
        vecClear(mdata->vnhat);
        vecClear(mdata->one);
        broadcast_mpz_to_vec(mdata->n, n);
        broadcast_mpz_to_vec(mdata->r, r);
        broadcast_mpz_to_vec(mdata->vrhat, mdata->rhat);
        broadcast_mpz_to_vec(mdata->vnhat, mdata->nhat);
        mpz_tdiv_r(mdata->rhat, r, n);
        broadcast_mpz_to_vec(mdata->one, mdata->rhat);
        for (j = 0; j < VECLEN; j++)
        {
            mdata->vrho[j] = mpz_get_ui(mdata->nhat) & VEC_MAXDIGIT;
        }

        broadcast_mpz_to_vec(vn, n);

        for (core = 0; core < 3; core++)
        {
            double mul_ns = -1.0, sqr_ns = -1.0, addsub_ns = -1.0;

            vecmulmod_ptr = mulptrs[core];
            vecsqrmod_ptr = sqrptrs[core];
            vecaddsubmod_ptr = &vec_simul_addsub52;

            /* fresh random operands, reduced and converted to monty form;
               keep lane-0 plain values for the GMP reference */
            for (j = 0; j < VECLEN; j++)
            {
                mpz_urandomb(t1, state, numwords * DIGITBITS);
                mpz_tdiv_r(t1, t1, n);
                mpz_urandomb(t2, state, numwords * DIGITBITS);
                mpz_tdiv_r(t2, t2, n);
                if (j == 0)
                {
                    mpz_set(pa, t1);
                    mpz_set(pb, t2);
                }
                insert_mpz_to_vec(va, t1, j);
                insert_mpz_to_vec(vb, t2, j);
            }
            vecmulmod_ptr(va, mdata->vrhat, va, vn, vs, mdata);
            vecmulmod_ptr(vb, mdata->vrhat, vb, vn, vs, mdata);

            if (validate)
            {
                extract_bignum_from_vec_to_mpz(ma, va, 0, numwords);
                extract_bignum_from_vec_to_mpz(mb, vb, 0, numwords);

                /* mul: monty(pa*pb) = (pa*pb mod n)*R mod n.
                   Tolerant mod-n comparison: the REDC wrappers are allowed
                   to return a result in [n, 2n) (final conditional
                   subtraction is their business, not ours). */
                vecmulmod_ptr(va, vb, vc, vn, vs, mdata);
                extract_bignum_from_vec_to_mpz(t1, vc, 0, numwords);
                mpz_mul(t2, pa, pb);
                mpz_mod(t2, t2, n);
                mpz_mul_2exp(t2, t2, numwords * DIGITBITS);
                mpz_mod(t2, t2, n);
                mpz_sub(t3, t1, t2);
                mpz_mod(t3, t3, n);
                if (mpz_cmp_ui(t3, 0) != 0)
                {
                    printf("vecbench VALIDATION FAILED: mul words=%d core=%s\n",
                        numwords, corenames[core]);
                    exit(1);
                }

                /* sqr: monty(pa*pa) */
                vecsqrmod_ptr(va, vc, vn, vs, mdata);
                extract_bignum_from_vec_to_mpz(t1, vc, 0, numwords);
                mpz_mul(t2, pa, pa);
                mpz_mod(t2, t2, n);
                mpz_mul_2exp(t2, t2, numwords * DIGITBITS);
                mpz_mod(t2, t2, n);
                mpz_sub(t3, t1, t2);
                mpz_mod(t3, t3, n);
                if (mpz_cmp_ui(t3, 0) != 0)
                {
                    printf("vecbench VALIDATION FAILED: sqr words=%d core=%s\n",
                        numwords, corenames[core]);
                    exit(1);
                }

                /* addsub: tolerant mod-n comparison of both outputs.
                   ma, mb hold the lane-0 monty-domain operands; the
                   reference sum/diff are reduced mod n and the vec
                   outputs are compared mod n so either output
                   convention (reduced or not) passes. */
                {
                    int bad = 0;
                    vecaddsubmod_ptr(va, vb, vc, vd, mdata);
                    extract_bignum_from_vec_to_mpz(t1, vc, 0, numwords);
                    extract_bignum_from_vec_to_mpz(t2, vd, 0, numwords);
                    /* sum reference: (ma + mb) mod n */
                    mpz_add(t3, ma, mb);
                    mpz_mod(t3, t3, n);
                    mpz_sub(t3, t1, t3);
                    mpz_mod(t3, t3, n);
                    if (mpz_cmp_ui(t3, 0) != 0) bad = 1;
                    /* diff reference: (ma - mb) mod n */
                    mpz_sub(t3, ma, mb);
                    mpz_mod(t3, t3, n);
                    mpz_sub(t3, t2, t3);
                    mpz_mod(t3, t3, n);
                    if (mpz_cmp_ui(t3, 0) != 0) bad = 1;
                    if (bad)
                    {
                        printf("vecbench VALIDATION FAILED: addsub words=%d core=%s\n",
                            numwords, corenames[core]);
                        exit(1);
                    }
                }
            }

            /* ---- time mul: ping-pong va,vb -> vc -> va ---- */
            iters = 32;
            gettimeofday(&tv0, NULL);
            for (it = 0; it < 32; it++)
            {
                vecmulmod_ptr(va, vb, vc, vn, vs, mdata);
                vecmulmod_ptr(vc, vb, va, vn, vs, mdata);
            }
            gettimeofday(&tv1, NULL);
            {
                double per = ((double)(tv1.tv_sec - tv0.tv_sec) +
                    1e-6 * (double)(tv1.tv_usec - tv0.tv_usec)) / 64.0;
                if (per <= 0.0)
                    per = 1e-9;
                iters = (int)(0.3 / per);
                if (iters < 32) iters = 32;
                if (iters > 20000000) iters = 20000000;
            }
            {
                double best = 1e30;
                for (rep = 0; rep < reps; rep++)
                {
                    gettimeofday(&tv0, NULL);
                    for (it = 0; it < iters; it++)
                    {
                        vecmulmod_ptr(va, vb, vc, vn, vs, mdata);
                        vecmulmod_ptr(vc, vb, va, vn, vs, mdata);
                    }
                    gettimeofday(&tv1, NULL);
                    {
                        double el = ((double)(tv1.tv_sec - tv0.tv_sec) +
                            1e-6 * (double)(tv1.tv_usec - tv0.tv_usec));
                        double ns = el * 1e9 / (2.0 * (double)iters);
                        if (ns < best) best = ns;
                    }
                }
                mul_ns = best;
            }

            /* ---- time sqr: ping-pong va -> vc -> va ---- */
            gettimeofday(&tv0, NULL);
            for (it = 0; it < 32; it++)
            {
                vecsqrmod_ptr(va, vc, vn, vs, mdata);
                vecsqrmod_ptr(vc, va, vn, vs, mdata);
            }
            gettimeofday(&tv1, NULL);
            {
                double per = ((double)(tv1.tv_sec - tv0.tv_sec) +
                    1e-6 * (double)(tv1.tv_usec - tv0.tv_usec)) / 64.0;
                if (per <= 0.0)
                    per = 1e-9;
                iters = (int)(0.3 / per);
                if (iters < 32) iters = 32;
                if (iters > 20000000) iters = 20000000;
            }
            {
                double best = 1e30;
                for (rep = 0; rep < reps; rep++)
                {
                    gettimeofday(&tv0, NULL);
                    for (it = 0; it < iters; it++)
                    {
                        vecsqrmod_ptr(va, vc, vn, vs, mdata);
                        vecsqrmod_ptr(vc, va, vn, vs, mdata);
                    }
                    gettimeofday(&tv1, NULL);
                    {
                        double el = ((double)(tv1.tv_sec - tv0.tv_sec) +
                            1e-6 * (double)(tv1.tv_usec - tv0.tv_usec));
                        double ns = el * 1e9 / (2.0 * (double)iters);
                        if (ns < best) best = ns;
                    }
                }
                sqr_ns = best;
            }

            /* ---- time addsub: fixed operands ---- */
            gettimeofday(&tv0, NULL);
            for (it = 0; it < 32; it++)
            {
                vecaddsubmod_ptr(va, vb, vc, vd, mdata);
            }
            gettimeofday(&tv1, NULL);
            {
                double per = ((double)(tv1.tv_sec - tv0.tv_sec) +
                    1e-6 * (double)(tv1.tv_usec - tv0.tv_usec)) / 32.0;
                if (per <= 0.0)
                    per = 1e-9;
                iters = (int)(0.3 / per);
                if (iters < 32) iters = 32;
                if (iters > 20000000) iters = 20000000;
            }
            {
                double best = 1e30;
                for (rep = 0; rep < reps; rep++)
                {
                    gettimeofday(&tv0, NULL);
                    for (it = 0; it < iters; it++)
                    {
                        vecaddsubmod_ptr(va, vb, vc, vd, mdata);
                    }
                    gettimeofday(&tv1, NULL);
                    {
                        double el = ((double)(tv1.tv_sec - tv0.tv_sec) +
                            1e-6 * (double)(tv1.tv_usec - tv0.tv_usec));
                        double ns = el * 1e9 / (double)iters;
                        if (ns < best) best = ns;
                    }
                }
                addsub_ns = best;
            }

            printf("%d,%s,%.1f,%.1f,%.1f,%d\n", numwords, corenames[core],
                mul_ns, sqr_ns, addsub_ns, iters);
            fflush(stdout);
        }

        vecFree(va);
        vecFree(vb);
        vecFree(vc);
        vecFree(vd);
        vecFree(vn);
        vecFree(vs);
        vec_monty_free(mdata);
    }

    gmp_randclear(state);
    mpz_clear(n);
    mpz_clear(r);
    mpz_clear(pa);
    mpz_clear(pb);
    mpz_clear(ma);
    mpz_clear(mb);
    mpz_clear(t1);
    mpz_clear(t2);
    return;
}


#endif