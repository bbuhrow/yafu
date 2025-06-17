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

#include <stdio.h>
#include "soe.h"
#include "soe_impl.h"
#include "ytools.h"
#include "threadpool.h"
#include <string.h>

// http://www.mersenneforum.org/showthread.php?t=21611
// http://www.mersenneforum.org/showthread.php?t=11900


typedef struct
{
    int *res_table;
    int *res_steps;
    uint32_t **res_nearest;
    uint32_t **res_nearest_class;
    int **res_classtab;

    soe_staticdata_t *sdata;
    thread_soedata_t *ddata;
} bitmap_userdata_t;

// define the fat-binary function pointers
uint32_t(*compute_8_bytes_ptr)(soe_staticdata_t*, uint32_t, uint64_t*, uint64_t);
void (*pre_sieve_ptr)(soe_dynamicdata_t*, soe_staticdata_t*, uint8_t*);
void(*sieve_line_ptr)(thread_soedata_t*);

void sieve_sync(void *vptr)
{
    tpool_t *tdata = (tpool_t *)vptr;
    soe_userdata_t *udata = (soe_userdata_t *)tdata->user_data;
    soe_staticdata_t *sdata = udata->sdata;
    thread_soedata_t *t = &udata->ddata[tdata->tindex];

    if (sdata->VFLAG > 1)
    {
        //don't print status if computing primes, because lots of routines within
        //yafu do this and they don't want this side effect
        printf("sieving: %d%%\r", 
            (int)((double)sdata->sync_count / (double)(sdata->numclasses)* 100.0));
        fflush(stdout);
    }

    sdata->num_found += t->linecount;

    if (t->ddata.min_sieved_val < sdata->min_sieved_val)
    {
        sdata->min_sieved_val = t->ddata.min_sieved_val;
    }

    return;
}

void sieve_dispatch(void *vptr)
{
    tpool_t *tdata = (tpool_t *)vptr;
    soe_userdata_t *udata = (soe_userdata_t *)tdata->user_data;
    soe_staticdata_t *sdata = udata->sdata;
    thread_soedata_t *t = &udata->ddata[tdata->tindex];

    // if not done, dispatch another line for sieving
    if (sdata->sync_count < sdata->numclasses)
    {
        t->current_line = (uint32_t)sdata->sync_count;
        tdata->work_fcn_id = 0;
        sdata->sync_count++;
    }
    else
    {
        tdata->work_fcn_id = tdata->num_work_fcn;
    }
    
    return;
}

void sieve_work_fcn(void *vptr)
{
    tpool_t *tdata = (tpool_t *)vptr;
    soe_userdata_t *udata = (soe_userdata_t *)tdata->user_data;
    soe_staticdata_t *sdata = udata->sdata;
    thread_soedata_t *t = &udata->ddata[tdata->tindex];

    if ((sdata->only_count == 0) || (sdata->analysis > 1) || (sdata->num_bitmap_primes > 0))
    {
        // if we are computing primes (not just counting) or if
        // bitmap sieving is enabled or if we're doing an analysis
        // that requires ordering, then we need to keep all lines
        // in memory at once.

        if ((sdata->analysis == 2) && (sdata->is_main_sieve))
        {
            // if counting/computing twins, then some residue classes can be
            // skipped entirely as they are never involved in a twin.
            if ((sdata->numclasses == 8) && ((t->current_line == 1) || (t->current_line == 6)))
            {
                t->linecount = 0;
                t->ddata.min_sieved_val = (1 << 31) - 1;
                return;
            }

            if (sdata->numclasses == 48)
            {
                // can skip classes: 5 8 11 12 15 18 19 20 21 26 27 28 29 32 35 36 39 42 
                switch (t->current_line)
                {
                case 5:
                case 8:
                case 11:
                case 12:
                case 15:
                case 18:
                case 19:
                case 20:
                case 21:
                case 26:
                case 27:
                case 28:
                case 29:
                case 32:
                case 35:
                case 36:
                case 39:
                case 42:
                    t->linecount = 0;
                    t->ddata.min_sieved_val = (1 << 31) - 1;
                    return;
                }
            }
        }

        sieve_line_ptr(t);
        trim_line(sdata, t->current_line);
        t->linecount = count_line(&t->sdata, t->current_line);
    }   
    else
    {
        t->sdata.lines[t->current_line] =
            (uint8_t *)xmalloc_align(t->sdata.numlinebytes * sizeof(uint8_t));
        sieve_line_ptr(t);
        trim_line(sdata, t->current_line);
        t->linecount = count_line(&t->sdata, t->current_line);
        align_free(t->sdata.lines[t->current_line]);        
    }

    return;
}

void bitmap_sync(void *vptr)
{
    tpool_t *tdata = (tpool_t *)vptr;
    bitmap_userdata_t *udata = (bitmap_userdata_t *)tdata->user_data;
    thread_soedata_t *t = &udata->ddata[tdata->tindex];

    t->startid = t->stopid;

    return;
}

void bitmap_dispatch(void *vptr)
{
    tpool_t *tdata = (tpool_t *)vptr;
    bitmap_userdata_t *udata = (bitmap_userdata_t *)tdata->user_data;
    thread_soedata_t *t = &udata->ddata[tdata->tindex];

    if (t->startid != t->stopid)
    {
        tdata->work_fcn_id = 0;
    }
    else
    {
        tdata->work_fcn_id = tdata->num_work_fcn;
    }

    return;
}

void bitmap_2class_work_fcn(void *vptr)
{
    tpool_t *tdata = (tpool_t *)vptr;
    bitmap_userdata_t *udata = (bitmap_userdata_t *)tdata->user_data;
    soe_staticdata_t *sdata = udata->sdata;
    thread_soedata_t *t = &udata->ddata[tdata->tindex];
    int i;

    int *res_table = udata->res_table;
    int *res_steps = udata->res_steps;
    uint32_t **res_nearest = udata->res_nearest;
    uint32_t **res_nearest_class = udata->res_nearest_class;
    int **res_classtab = udata->res_classtab;

    uint64_t llimit = sdata->lowlimit + (uint64_t)t->startid * (uint64_t)sdata->prodN * sdata->FLAGSIZE;
    uint64_t hlimit = sdata->lowlimit + (uint64_t)t->stopid * (uint64_t)sdata->prodN * sdata->FLAGSIZE;

    if (hlimit > sdata->orig_hlimit)
        hlimit = sdata->orig_hlimit;

    if (sdata->VFLAG > 2)
    {
        printf("thread %d working on range %lu - %lu\n", tdata->tindex,
            llimit, hlimit);
    }

    // do bitmap sieve work
    for (i = sdata->bitmap_start_id; i < sdata->pboundi; i++)
    {
        uint64_t interval = hlimit - llimit;
        uint64_t bloc;
        uint32_t prime = sdata->sieve_p[i];
        int pclass = res_table[prime % 6];    // residue class of current prime
        int bclass;                 // residue class of current bit
        int bclassnum;              // column in residue_pattern_mod30 of bclass

        // this is only worth doing for 2 classes, because in that
        // case is it possible that many of the sieve primes will 
        // be less than the sieve interval (which is 10^8 or less).
        if (prime > interval)
            break;

        // first bit offset above lowlimit
        bloc = prime - llimit % prime;

        // get it on a residue class
        bclassnum = res_nearest_class[pclass][bloc % 6];
        bloc = bloc + (uint64_t)res_nearest[pclass][bloc % 6] * (uint64_t)prime;
        bclass = bloc % 6;

        while (bloc < interval)
        {
            uint32_t brloc = t->startid * sdata->FLAGSIZE + (bloc - bclass) / 6;
            sdata->lines[res_table[bclass]][brloc >> 3] &= sdata->masks[brloc & 7];

            bloc += (uint64_t)res_steps[bclassnum] * (uint64_t)prime;

            bclassnum ^= 1;
            bclass = res_classtab[pclass][bclassnum];
        }
    }

    for (; i < sdata->pboundi; i++)
    {
        uint64_t interval = hlimit - llimit;
        uint64_t bloc;
        uint32_t prime = sdata->sieve_p[i];
        int pclass = res_table[prime % 6];    // residue class of current prime
        int bclass;                 // residue class of current bit
        //int bclassnum;              // column in residue_pattern_mod30 of bclass

        // first bit offset above lowlimit
        bloc = prime - llimit % prime;

        // get it on a residue class
        //bclassnum = res_nearest_class[pclass][bloc % 6];
        bloc = bloc + (uint64_t)res_nearest[pclass][bloc % 6] * (uint64_t)prime;
        bclass = bloc % 6;

        if (bloc < interval)
        {
            uint32_t brloc = t->startid * sdata->FLAGSIZE + (bloc - bclass) / 6;
            sdata->lines[res_table[bclass]][brloc >> 3] &= sdata->masks[brloc & 7];
        }
    }

    return;
}

void bitmap_8class_work_fcn(void *vptr)
{
    tpool_t *tdata = (tpool_t *)vptr;
    bitmap_userdata_t *udata = (bitmap_userdata_t *)tdata->user_data;
    soe_staticdata_t *sdata = udata->sdata;
    thread_soedata_t *t = &udata->ddata[tdata->tindex];
    int i;

    int *res_table = udata->res_table;
    int *res_steps = udata->res_steps;

    uint64_t llimit = sdata->lowlimit + (uint64_t)t->startid * (uint64_t)sdata->prodN * sdata->FLAGSIZE;
    uint64_t hlimit = sdata->lowlimit + (uint64_t)t->stopid * (uint64_t)sdata->prodN * sdata->FLAGSIZE;

    if (hlimit > sdata->orig_hlimit)
        hlimit = sdata->orig_hlimit;

    if (sdata->VFLAG > 2)
    {
        printf("thread %d working on range %lu - %lu\n", tdata->tindex,
            llimit, hlimit);
    }

    // do bitmap sieve work
    for (i = sdata->bitmap_start_id; i < sdata->pboundi; i++)
    {
        uint64_t interval = hlimit - llimit;
        uint64_t bloc;
        uint32_t prime = sdata->sieve_p[i];
        int pclass = res_table[prime % 30];    // residue class of current prime
        int bclass;                 // residue class of current bit
        int bclassnum;              // column in residue_pattern_mod30 of bclass

        // once prime exceeds the interval size we can use faster methods.
        if (prime > interval)
            break;

        // first bit offset above lowlimit
        bloc = prime - llimit % prime;

        // get it on a residue class
        bclassnum = udata->res_nearest_class[pclass][bloc % 30];
        bloc = bloc + (uint64_t)udata->res_nearest[pclass][bloc % 30] * (uint64_t)prime;
        bclass = bloc % 30;

        // sieve the entire space, jumping around residue classes
        while (bloc < interval)
        {
            // we need to compute (bloc - bclass) / prodN
            // to get the bit location in residue space, along with the 
            // bclassnum to select the correct line.  The classnum is easy
            // to track, but the division in the bit location calculation 
            // could be a problem.
            uint32_t brloc = t->startid * sdata->FLAGSIZE + (bloc - bclass) / 30;
            sdata->lines[res_table[bclass]][brloc >> 3] &= sdata->masks[brloc & 7];

            // increment the bit offset to the next valid residue class
            bloc += (uint64_t)res_steps[bclassnum] * (uint64_t)prime;

            // increment the residue class of the bit offset
            bclassnum++;
            bclassnum &= 7;

            bclass = udata->res_classtab[pclass][bclassnum];
        }
    }

    for (; i < sdata->pboundi; i++)
    {
        uint64_t interval = hlimit - llimit;
        uint64_t bloc;
        uint32_t prime = sdata->sieve_p[i];
        int pclass = res_table[prime % 30];    // residue class of current prime
        int bclass;                 // residue class of current bit
        //int bclassnum;              // column in residue_pattern_mod30 of bclass

        // first bit offset above lowlimit
        bloc = prime - llimit % prime;

        // get it on a residue class
        //bclassnum = res_nearest_class[pclass][bloc % 30];
        bloc = bloc + (uint64_t)udata->res_nearest[pclass][bloc % 30] * (uint64_t)prime;
        bclass = bloc % 30;

        if (bloc < interval)
        {
            uint32_t brloc = t->startid * sdata->FLAGSIZE + (bloc - bclass) / 30;
            sdata->lines[res_table[bclass]][brloc >> 3] &= sdata->masks[brloc & 7];
        }
    }

    return;
}

void bitmap_48class_work_fcn(void *vptr)
{
    tpool_t *tdata = (tpool_t *)vptr;
    bitmap_userdata_t *udata = (bitmap_userdata_t *)tdata->user_data;
    soe_staticdata_t *sdata = udata->sdata;
    thread_soedata_t *t = &udata->ddata[tdata->tindex];
    int i;

    int *res_table = udata->res_table;
    int *res_steps = udata->res_steps;
    uint32_t **res_nearest = udata->res_nearest;
    uint32_t **res_nearest_class = udata->res_nearest_class;
    int **res_classtab = udata->res_classtab;

    uint64_t llimit = sdata->lowlimit + (uint64_t)t->startid * sdata->prodN * sdata->FLAGSIZE;
    uint64_t hlimit = sdata->lowlimit + (uint64_t)t->stopid * sdata->prodN * sdata->FLAGSIZE;

    if (hlimit > sdata->orig_hlimit)
        hlimit = sdata->orig_hlimit;

    if (sdata->VFLAG > 2)
    {
        printf("thread %d working on range %lu - %lu\n", tdata->tindex, llimit, hlimit);
    }

    // do bitmap sieve work
    for (i = sdata->bitmap_start_id; i < sdata->pboundi; i++)
    {
        uint64_t interval = hlimit - llimit;
        uint64_t bloc;
        uint32_t prime = sdata->sieve_p[i];
        int pclass = res_table[prime % 210];    // residue class of current prime
        int bclass;                 // residue class of current bit
        uint32_t bclassnum;              // column in residue_pattern_mod30 of bclass

        if (prime > interval)
            break;

        // compute the bit location for this prime relative to the
        // starting lowlimit for this thread.  Also the class 
        // associated with that location and this prime.
        bloc = (uint64_t)prime - llimit % (uint64_t)prime;
        bclassnum = res_nearest_class[pclass][bloc % 210ULL];
        bloc = bloc + (uint64_t)res_nearest[pclass][bloc % 210ULL] * (uint64_t)prime;
        bclass = bloc % 210ULL;

        while (bloc < interval)
        {
            uint64_t brloc = (uint64_t)t->startid * sdata->FLAGSIZE + (bloc - (uint64_t)bclass) / 210ULL;
            sdata->lines[res_table[bclass]][brloc >> 3ULL] &= sdata->masks[brloc & 7ULL];

            bloc += ((uint64_t)res_steps[bclassnum] * (uint64_t)prime);

            bclassnum++;
            if (bclassnum == 48) bclassnum = 0;
            bclass = res_classtab[pclass][bclassnum];
        }
    }

    for (; i < sdata->pboundi; i++)
    {
        uint64_t interval = hlimit - llimit;
        uint64_t bloc;
        uint32_t prime = sdata->sieve_p[i];
        int pclass = res_table[prime % 210];    // residue class of current prime
        int bclass;                 // residue class of current bit
        //int bclassnum;              // column in residue_pattern_mod30 of bclass

        // first bit offset above lowlimit
        bloc = prime - llimit % prime;

        // get it on a residue class
        //bclassnum = res_nearest_class[pclass][bloc % 210];
        bloc = bloc + (uint64_t)res_nearest[pclass][bloc % 210] * (uint64_t)prime;
        bclass = bloc % 210;

        if (bloc < interval)
        {
            uint64_t brloc = (uint64_t)t->startid * sdata->FLAGSIZE + (bloc - (uint64_t)bclass) / 210ULL;
            sdata->lines[res_table[bclass]][brloc >> 3] &= sdata->masks[brloc & 7];
        }
    }

    return;
}


uint64_t spSOE(soe_staticdata_t *sdata, mpz_t offset,
	uint64_t lowlimit, uint64_t *highlimit, int count, uint64_t *primes)
{
	/*
	if count == 1, then the primes are simply counted, and not 
	explicitly calculated and saved in *primes.

	otherwise, store primes in the provided *primes array

	in either case, return the number of primes found
	*/
    // timing
    double t;
    struct timeval tstart, tstop;

	// keep track of how much memory we've used
	uint64_t allocated_bytes = 0;

	// structure of static info
    uint32_t* sieve_p = sdata->sieve_p;
    uint32_t num_sp = sdata->num_sp;
    int VFLAG = sdata->VFLAG;
    int THREADS = sdata->THREADS;
    info_t info;

	//thread data holds all data needed during sieving
	thread_soedata_t *thread_data;		//an array of thread data objects

	//*********************** BEGIN ******************************//
    sdata->only_count = count;

    if (VFLAG > 1)
    {
        gettimeofday(&tstart, NULL);
    }

    ytools_get_computer_info(&info, 0);
    sdata->has_avx2 = info.AVX2;
    sdata->has_avx512f = info.AVX512F;
    sdata->has_bmi1 = info.BMI1;
    sdata->has_bmi2 = info.BMI2;


	// sanity check the input
	if (check_input(*highlimit, lowlimit, num_sp, sieve_p, sdata, offset))
		return 0;

	// determine what kind of sieve to use based on the input
	get_numclasses(*highlimit, lowlimit, sdata);

    if (mpz_cmp_ui(offset, 0) > 0)
    {
        sdata->use_monty = 0;
    }

	// allocate and initialize some stuff
	allocated_bytes += init_sieve(sdata);
	*highlimit = sdata->highlimit;

	// allocate thread data structure
	thread_data = (thread_soedata_t *)malloc(THREADS * sizeof(thread_soedata_t));

	// find all roots of prime with prodN.  These are used when finding offsets.
	getRoots(sdata, thread_data);

	// init bucket sieving
	set_bucket_depth(sdata);

	// initialize stuff used in thread structures.
	// this is necessary even if THREADS = 1;	
	allocated_bytes += alloc_threaddata(sdata, thread_data);

    mpz_set(offset, sdata->offset);

	if (VFLAG > 2)
	{	
        char strfeatures[80];

        if (sdata->sieve_range == 0)
        {
            // normal sieve
            printf("finding requested range %" PRIu64 " to %" PRIu64 "\n",
                sdata->orig_llimit, sdata->orig_hlimit);
            printf("sieving range %" PRIu64 " to %" PRIu64 "\n",
                sdata->lowlimit, *highlimit);
            printf("using %" PRIu64 " primes, max prime = %" PRIu64 "  \n",
                sdata->pboundi, sdata->sieve_p[sdata->pboundi - 1]);
        }
        else
        {
            // sieve to depth - not guarenteed to find all primes
            gmp_printf("sieving range 0 to %" PRIu64 " from offset %Zd\n", 
                *highlimit, sdata->offset);
            gmp_printf("requested range is %Zd + %lu:%lu\n",
                sdata->offset, sdata->orig_llimit, sdata->orig_hlimit);
            printf("using %" PRIu64 " primes, max prime = %" PRIu64 "  \n",
                sdata->pboundi, sdata->sieve_p[sdata->pboundi - 1]);
        }

		printf("using %u residue classes\n",sdata->numclasses);

        if (sdata->use_monty)
            printf("using Montgomery enabled offset computations\n");

        strcpy(strfeatures, "");

#if defined(USE_AVX2)
        if (sdata->has_avx2)
        {
            sprintf(strfeatures, "AVX2");
        }
#endif

#if defined(USE_BMI2)
        if (sdata->has_bmi2)
        {
            if (strlen(strfeatures) > 0)
            {
                sprintf(strfeatures, "%s, BMI2", strfeatures);
            }
            else
            {
                sprintf(strfeatures, "BMI2");
            }
        }
#endif

#ifdef USE_AVX512F
        if (sdata->has_avx512f)
        {
            if (strlen(strfeatures) > 0)
            {
                sprintf(strfeatures, "%s, AVX512", strfeatures);
            }
            else
            {
                sprintf(strfeatures, "AVX512");
            }
        }
#endif

        if (strlen(strfeatures) > 0)
        {
            printf("Using cpu features: %s\n", strfeatures);
        }
        else
        {
            printf("Using cpu features: x86-64\n");
        }

		printf("lines have %" PRIu64 " bytes and %" PRIu64 " flags\n",
            sdata->numlinebytes,sdata->numlinebytes * 8);
		printf("lines broken into %" PRIu64 " blocks of size %u\n",
            sdata->blocks, sdata->SOEBLOCKSIZE);
		printf("blocks contain %u flags and span %" PRIu64 " integers\n", 
            sdata->FLAGSIZE, sdata->blk_r);
        
        if (sdata->num_bucket_primes > 0)
		{
			printf("bucket sieving %u primes > %u\n",
				sdata->num_bucket_primes, sdata->sieve_p[sdata->bucket_start_id]);
			printf("allocating space for %u hits per bucket\n", sdata->bucket_alloc);
			printf("allocating space for %u hits per large bucket\n", sdata->large_bucket_alloc);
		}
		if (sdata->num_bitmap_primes > 0)
		{
			printf("bitmap sieving %u primes > %u\n",
				sdata->num_bitmap_primes, sdata->sieve_p[sdata->bitmap_start_id]);
		}
		printf("using %" PRIu64 " bytes for sieving storage\n",allocated_bytes);
	}

    if (VFLAG > 1)
    {
        gettimeofday(&tstop, NULL);
        t = ytools_difftime(&tstart, &tstop);

        if (VFLAG > 2)
        {
            printf("setup took %1.6f seconds\n", t);
        }
    }

    if (sdata->num_bitmap_primes > 0)
    {
        uint8_t *line = sdata->lines[0];
        uint32_t linesize = sdata->FLAGSIZE * sdata->blocks;
        int i;

        int *res_table;
        int *res_steps;
        uint32_t **res_nearest;
        uint32_t **res_nearest_class;
        int **res_classtab;

        tpool_t *tpool_data;
        bitmap_userdata_t udata;        

        int threads = MIN(sdata->blocks, THREADS);
        int blocks_per_thread = sdata->blocks / threads;
        uint32_t b = 0;

		if (VFLAG > 1)
		{
			gettimeofday(&tstart, NULL);
            printf("commencing bitmap sieve\n");
		}

        res_table = (int *)malloc(sdata->prodN * sizeof(int));
        memset(res_table, -1, sdata->prodN * sizeof(int));
        for (i = 0; i < sdata->numclasses; i++)
            res_table[sdata->rclass[i]] = i;

		if (VFLAG > 3)
		{
			printf("res_table: ");
			for (i = 0; i < sdata->prodN; i++)
			{
				printf("%d ", res_table[i]);
				if ((i & 7) == 0)
					printf("\n");
			}
			printf("\n");
		}
		
        res_steps = (int *)malloc(sdata->numclasses * sizeof(int));
        for (i = 1; i < sdata->numclasses; i++)
            res_steps[i - 1] = sdata->rclass[i] - sdata->rclass[i - 1];
        res_steps[i - 1] = sdata->prodN - sdata->rclass[i - 1] + 1;

		if (VFLAG > 3)
		{
			printf("res_steps: ");
			for (i = 0; i < sdata->numclasses; i++)
			{
				printf("%d ", sdata->rclass[i]);
				if ((i & 7) == 0)
					printf("\n");
			}
			printf("\n");
		}

        res_classtab = (int **)malloc(sdata->numclasses * sizeof(int *));
        for (i = 0; i < sdata->numclasses; i++)
        {
            int j;
            res_classtab[i] = (int *)malloc(sdata->numclasses * sizeof(int));
            res_classtab[i][0] = sdata->rclass[i];
            for (j = 1; j < sdata->numclasses; j++)
                res_classtab[i][j] = (res_classtab[i][j-1] + 
                sdata->rclass[i] * res_steps[j-1]) % sdata->prodN;

			if (VFLAG > 3)
			{
				printf("res_classtab[%d]: ", i);
				for (j = 0; j < sdata->numclasses; j++)
				{
					printf("%d ", res_classtab[i][j]);
				}
				printf("\n");
			}
        }

        res_nearest = (uint32_t **)malloc(sdata->numclasses * sizeof(uint32_t *));
        res_nearest_class = (uint32_t **)malloc(sdata->numclasses * sizeof(uint32_t *));
        for (i = 0; i < sdata->numclasses; i++)
        {
            int j;
            res_nearest[i] = (uint32_t *)malloc(sdata->prodN * sizeof(uint32_t));
            res_nearest_class[i] = (uint32_t *)malloc(sdata->prodN * sizeof(uint32_t));

            for (j = 0; j < sdata->prodN; j++)
            {
                int k;
                uint32_t x = j;
                res_nearest[i][j] = 0;
                while (res_table[x] < 0)
                {
                    res_nearest[i][j]++;
                    x += sdata->rclass[i];
                    if (x >= sdata->prodN) x -= sdata->prodN;
                }

                for (k = 0; k < sdata->numclasses; k++)
                    if (res_classtab[i][k] == x)
                        res_nearest_class[i][j] = k;
            }
        }

		if (VFLAG > 3)
		{
			printf("res_nearest: ");
			for (i = 0; i < sdata->numclasses; i++)
			{
				int j;
				for (j = 0; j < sdata->prodN; j++)
				{
					printf("%d ", res_nearest[i][j]);
				}
				printf("\n");
			}
            printf("res_nearest_class: ");
            for (i = 0; i < sdata->numclasses; i++)
            {
                int j;
                for (j = 0; j < sdata->prodN; j++)
                {
                    printf("%d ", res_nearest_class[i][j]);
                }
                printf("\n");
            }
		}

        memset(line, 255, linesize * sdata->numclasses / 8);

        if (VFLAG > 1)
        {
            printf("commencing bitmap sieve starting at sieve prime index %u\n", 
                sdata->bitmap_start_id);
        }

        /*
        How to make this multi-threaded...
        Splitting up by primes index is out because any prime can hit
        all blocks/classes.
        Splitting up by class is out for the same reason.
        The threads must sieve separate columns of the bitmap.
        This means each thread will use a different lowlimit... I think that
        is all that needs to change (and interval, of course).
        */

        if (sdata->numclasses == 480)
        {
            printf("bitmap sieving should only be enabled with 2,8,or 48 classes!\n");
            exit(1);
        }

        udata.sdata = sdata;
        udata.ddata = thread_data;
        
        udata.res_classtab = (int**)malloc(sdata->numclasses * sizeof(int*));
        for (i = 0; i < sdata->numclasses; i++)
        {
            int j;
            udata.res_classtab[i] = (int*)malloc(sdata->numclasses * sizeof(int));
            udata.res_classtab[i][0] = sdata->rclass[i];
            for (j = 1; j < sdata->numclasses; j++)
                udata.res_classtab[i][j] = (udata.res_classtab[i][j - 1] +
                    sdata->rclass[i] * res_steps[j - 1]) % sdata->prodN;
        }

        udata.res_nearest = (uint32_t**)malloc(sdata->numclasses * sizeof(uint32_t*));
        udata.res_nearest_class = (uint32_t**)malloc(sdata->numclasses * sizeof(uint32_t*));
        for (i = 0; i < sdata->numclasses; i++)
        {
            int j;
            udata.res_nearest[i] = (uint32_t*)malloc(sdata->prodN * sizeof(uint32_t));
            udata.res_nearest_class[i] = (uint32_t*)malloc(sdata->prodN * sizeof(uint32_t));

            for (j = 0; j < sdata->prodN; j++)
            {
                int k;
                uint32_t x = j;
                udata.res_nearest[i][j] = 0;
                while (res_table[x] < 0)
                {
                    udata.res_nearest[i][j]++;
                    x += sdata->rclass[i];
                    if (x >= sdata->prodN) x -= sdata->prodN;
                }

                for (k = 0; k < sdata->numclasses; k++)
                    if (udata.res_classtab[i][k] == x)
                        udata.res_nearest_class[i][j] = k;
            }
        }

        udata.res_steps = res_steps;
        udata.res_table = res_table;

        tpool_data = tpool_setup(threads, NULL, NULL, &bitmap_sync,
            &bitmap_dispatch, &udata);

        for (i = 0; i < threads; i++, b += blocks_per_thread)
        {
            thread_data[i].startid = b;
            thread_data[i].stopid = b + blocks_per_thread;
        }
        thread_data[i - 1].stopid = sdata->blocks;

        if (sdata->numclasses == 2)
        {
            if (threads > 1)
            {
                tpool_add_work_fcn(tpool_data, &bitmap_2class_work_fcn);
                tpool_go(tpool_data);
            }
            else
                bitmap_2class_work_fcn(tpool_data);
        }
        else if (sdata->numclasses == 8)
        {
            if (threads > 1)
            {
                tpool_add_work_fcn(tpool_data, &bitmap_8class_work_fcn);
                tpool_go(tpool_data);
            }
            else
                bitmap_8class_work_fcn(tpool_data);
        }
        else if (sdata->numclasses == 48)
        {            
            if (threads > 1)
            {
                tpool_add_work_fcn(tpool_data, &bitmap_48class_work_fcn);
                tpool_go(tpool_data);
            }
            else
                bitmap_48class_work_fcn(tpool_data);
        }

        free(tpool_data);

        free(res_table);
        free(res_steps);
        for (i = 0; i < sdata->numclasses; i++)
        {
            free(res_classtab[i]);
            free(res_nearest[i]);
            free(res_nearest_class[i]);
        }
        free(res_classtab);
        free(res_nearest);
        free(res_nearest_class);

        if (VFLAG > 1)
        {
            gettimeofday(&tstop, NULL);
            t = ytools_difftime(&tstart, &tstop);

            if (VFLAG > 2)
            {
                printf("bitmap sieve took %1.6f seconds\n", t);
            }
        }
    }

    //printf("commencing sieve from %lu - %lu (originally %lu - %lu)\n",
    //    sdata->lowlimit, sdata->highlimit, sdata->orig_llimit, sdata->orig_hlimit);

	//get 'r done.
	do_soe_sieving(sdata, thread_data, count);

	//finish up
	finalize_sieve(sdata, thread_data, count, primes);

	return sdata->num_found;
}

void do_soe_sieving(soe_staticdata_t *sdata, thread_soedata_t *thread_data, int count)
{
	uint64_t i;
	//timing
	double t;
	struct timeval tstart, tstop;

	if (sdata->VFLAG > 1)
	{
		gettimeofday(&tstart, NULL);
	}

    // threading structures
    tpool_t *tpool_data;
    soe_userdata_t udata;

	//main sieve, line by line
    sdata->num_found = 0;
    sdata->only_count = count;

    udata.sdata = sdata;
    udata.ddata = thread_data;
    tpool_data = tpool_setup(sdata->THREADS, NULL, NULL, &sieve_sync,
        &sieve_dispatch, &udata);

    if (sdata->THREADS == 1)
    {
        thread_soedata_t *t = &thread_data[0];
        sdata->sync_count = 0;
        for (i = 0; i < sdata->numclasses; i++)
        {
            t->current_line = i;
            sieve_work_fcn(tpool_data);
            sieve_sync(tpool_data);
            sdata->sync_count++;
        }
    }
    else
    {
        sdata->sync_count = 0;
        tpool_add_work_fcn(tpool_data, &sieve_work_fcn);
        tpool_go(tpool_data);
    }

	if (sdata->VFLAG > 1)
	{
		gettimeofday(&tstop, NULL);
		t = ytools_difftime(&tstart, &tstop);
		printf("linesieve took %1.6f seconds\n", t);
	}
    
    free(tpool_data);

    // to test: make this a stop fcn
    for (i = 0; i < sdata->THREADS; i++)
    {
        align_free(thread_data[i].ddata.offsets);
    }

	return;
}

void finalize_sieve(soe_staticdata_t *sdata, 
	thread_soedata_t *thread_data, int count, uint64_t *primes)
{
	uint64_t i, j = 0, num_p = sdata->num_found;

	if (count)
	{
		//add in relevant sieving primes not captured in the flag arrays
		uint64_t ui_offset;

		if (sdata->sieve_range)
		{
            if (mpz_size(sdata->offset) == 1)
            {
                ui_offset = mpz_get_ui(sdata->offset);

                if (ui_offset > sdata->pbound)
                {
                    // huge offset, we don't need to add any primes.
                    ui_offset = 0;
                    sdata->min_sieved_val = 0;
                }
            }
			else
			{
				// huge offset, we don't need to add any primes.
				ui_offset = 0;
				sdata->min_sieved_val = 0;
			}
		}
        else
        {
            ui_offset = 0;
        }
		
        if (sdata->sieve_range)
        {
            sdata->min_sieved_val += ui_offset;
        }

        if (sdata->VFLAG > 2)
        {
            printf("num_found by main sieve = %lu\n", num_p);
        }

		// PRIMES is already sized appropriately by the wrapper
		// load in the sieve primes that we need
        //printf("min_sieved_val = %lu\n", sdata->min_sieved_val);
        //printf("bucket_start_id = %u\n", sdata->bucket_start_id);
        if ((sdata->analysis == 2) && (sdata->is_main_sieve == 1))
        {
            // count twins from the stored sieve lines
            num_p = count_twins(sdata, thread_data);

            //printf("main sieve found %lu twins, adding twins within sieve primes\n", num_p);
            i = 0;
            while (((uint64_t)sdata->sieve_p[i] < sdata->min_sieved_val) && (i < sdata->bucket_start_id))
            {
                // if we are doing twin prime or prime gap analysis 
                // then do that analysis here for the sieving primes
                // that are within the requested range.
                if (sdata->sieve_p[i] >= (sdata->orig_llimit + ui_offset))
                {
                    if ((sdata->sieve_p[i + 1] - sdata->sieve_p[i]) == 2)
                    {
                        num_p++;
                    }
                }
                i++;
            }
        }
        else
        {
            // num_p holds the primes counted in the main sieve.
            // add in sieve primes that are within the requested interval.        
            i = 0;
            while (((uint64_t)sdata->sieve_p[i] < sdata->min_sieved_val) && (i < sdata->bucket_start_id))
            {
                if (sdata->sieve_p[i] >= (sdata->orig_llimit + ui_offset))
                {
                    //printf("%u ", sdata->sieve_p[i]);
                    num_p++;
                }
                i++;
            }
        }
	}
	else
	{
		// now we need to raster vertically down the lines and horizontally
		// across the lines in order to compute the primes in order.

		// first put in any sieve primes if necessary.
		// if we are in this loop, and we are sieving a range, then offset
		// is a single precision number and we need to increment the 'prime'
		// we found above by it.
		uint64_t ui_offset;
			
		if (sdata->sieve_range)
		{
            if (mpz_size(sdata->offset) == 1)
            {
                ui_offset = mpz_get_ui(sdata->offset);
                if (ui_offset > sdata->pbound)
                {
                    // huge offset, we don't need to add any primes.
                    ui_offset = 0;
                    sdata->min_sieved_val = 0;
                }
            }
			else
			{
				// huge offset, we don't need to add any primes.
				ui_offset = 0;
				sdata->min_sieved_val = 0;
			}
		}
        else
        {
            ui_offset = 0;
        }
			
        if (sdata->sieve_range)
        {
            sdata->min_sieved_val += ui_offset;
        }

        if (sdata->VFLAG > 2)
        {
            printf("num_found by main sieve = %lu\n", num_p);
        }

		// PRIMES is already sized appropriately by the wrapper
		// load in the sieve primes that we need
		j = 0;
		i = 0;
        //printf("min_sieved_val = %lu\n", sdata->min_sieved_val);
        //printf("bucket_start_id = %u\n", sdata->bucket_start_id);

		while (((uint64_t)sdata->sieve_p[i] < sdata->min_sieved_val) && (i < sdata->bucket_start_id))
		{
            // if we are doing twin prime or prime gap analysis 
            // then do that analysis here for the sieving primes
            // that are within the requested range.
            if (sdata->sieve_p[i] >= (sdata->orig_llimit + ui_offset))
            {
                if ((sdata->analysis == 2) && (sdata->is_main_sieve == 1))
                {
                    if ((sdata->sieve_p[i + 1] - sdata->sieve_p[i]) == 2)
                    {
                        primes[j++] = (uint64_t)sdata->sieve_p[i];
                    }
                }
                else
                {
                    primes[j++] = (uint64_t)sdata->sieve_p[i];
                }
            }
            i++;
		}

		//and then the primes in the lines
		num_p = primes_from_lineflags(sdata, thread_data, j, primes);

	}

	// update count of found primes
	sdata->num_found = num_p;

    for (i = 0; i < sdata->THREADS; i++)
    {
        thread_soedata_t* thread = thread_data + i;
        free(thread->ddata.pbounds);
        if ((sdata->analysis > 1) || (sdata->gapmin > 0))
        {
            free(thread->ddata.analysis_carry_data);
        }
        align_free(thread->ddata.presieve_scratch);
    }

    if (sdata->num_bucket_primes > 0)
    {
        for (i = 0; i < sdata->THREADS; i++)
        {
            thread_soedata_t* thread = thread_data + i;

            free(thread->ddata.bucket_hits);
            if (thread->ddata.large_bucket_hits != NULL)
                free(thread->ddata.large_bucket_hits);
            for (j = 0; j < thread->sdata.blocks; j++)
            {
                free(thread->ddata.sieve_buckets[j]);
                if (thread->ddata.large_sieve_buckets != NULL)
                    free(thread->ddata.large_sieve_buckets[j]);
            }
            free(thread->ddata.sieve_buckets);
            if (thread->ddata.large_sieve_buckets != NULL)
                free(thread->ddata.large_sieve_buckets);
        }
    }

	//if (!sdata->only_count)
    if ((sdata->only_count == 0) || (sdata->num_bitmap_primes > 0))
	{
		//align_free(sdata->lines[0]);
        for (i = 0; i < sdata->numclasses; i++)
        {
            //sdata->lines[i] = sdata->lines[0] + i * numlinebytes;
            align_free(sdata->lines[i]);
        }
        align_free(sdata->lines);
	}
    else
    {
        align_free(sdata->lines);
    }

    align_free(sdata->r2modp);
    align_free(sdata->pinv);
    align_free(sdata->root);
    align_free(sdata->lower_mod_prime);
	free(thread_data);
	free(sdata->rclass);

	return;
}


