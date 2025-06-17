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

#include "soe.h"
#include "soe_impl.h"
#include "ytools.h"
#include "threadpool.h"
#if defined(_MSC_VER) && defined(__clang__)
#include <x86intrin.h>
#else
#include <immintrin.h>
#endif
#include <stdint.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <math.h>


//for testing one of 8 bits in a byte in one of 8 lines.
//bit num picks the row, lines num picks the col.	
const uint64_t nmasks64[8][8] = {
    { 1ULL, 256ULL, 65536ULL, 16777216ULL, 4294967296ULL, 1099511627776ULL, 281474976710656ULL, 72057594037927936ULL },
    { 2ULL, 512ULL, 131072ULL, 33554432ULL, 8589934592ULL, 2199023255552ULL, 562949953421312ULL, 144115188075855872ULL },
    { 4ULL, 1024ULL, 262144ULL, 67108864ULL, 17179869184ULL, 4398046511104ULL, 1125899906842624ULL, 288230376151711744ULL },
    { 8ULL, 2048ULL, 524288ULL, 134217728ULL, 34359738368ULL, 8796093022208ULL, 2251799813685248ULL, 576460752303423488ULL },
    { 16ULL, 4096ULL, 1048576ULL, 268435456ULL, 68719476736ULL, 17592186044416ULL, 4503599627370496ULL, 1152921504606846976ULL },
    { 32ULL, 8192ULL, 2097152ULL, 536870912ULL, 137438953472ULL, 35184372088832ULL, 9007199254740992ULL, 2305843009213693952ULL },
    { 64ULL, 16384ULL, 4194304ULL, 1073741824ULL, 274877906944ULL, 70368744177664ULL, 18014398509481984ULL, 4611686018427387904ULL },
    { 128ULL, 32768ULL, 8388608ULL, 2147483648ULL, 549755813888ULL, 140737488355328ULL, 36028797018963968ULL, 9223372036854775808ULL } };

void compute_primes_dispatch(void *vptr)
{
    tpool_t *tdata = (tpool_t *)vptr;
    soe_userdata_t *t = (soe_userdata_t *)tdata->user_data;
    soe_staticdata_t *sdata = t->sdata;

    // launch one range of computation for each thread.  don't really
    // need a threadpool for this, but the infrastructure is there...
    if (sdata->sync_count < sdata->THREADS)
    {
        tdata->work_fcn_id = 0;
        sdata->sync_count++;
    }
    else
    {
        tdata->work_fcn_id = tdata->num_work_fcn;
    }

    return;
}

void compute_primes_work_fcn(void *vptr)
{
    tpool_t *tdata = (tpool_t *)vptr;
    soe_userdata_t *udata = (soe_userdata_t *)tdata->user_data;
    soe_staticdata_t *sdata = udata->sdata;
    thread_soedata_t *t = &udata->ddata[tdata->tindex];
    int i;

    if (sdata->THREADS > 1)
    {
        t->linecount = 0;
    }

#if defined(USE_BMI2) || defined(USE_AVX512F)
    if ((sdata->has_bmi2) && (sdata->numclasses <= 48)) // && (!(sdata->analysis > 1)))
    {
        for (i = t->startid; i < t->stopid; i += 8)
        {
            t->linecount = compute_8_bytes_bmi2(sdata, t->linecount, t->ddata.primes, i);

            // if searching for prime constellations, then we need to look at the last 
            // flags of this block of 8 bytes and the first flags of the next one.
            // we do this by loading the trailing bits of this block into a carry
            // register.  The rest is handled by the block analysis function.
            if ((sdata->analysis == 2) && (sdata->is_main_sieve == 1))
            {
                uint8_t* lastline = sdata->lines[sdata->numclasses - 1];
                uint8_t* firstline = sdata->lines[0];

                if ((i + 8) < t->stopid)
                {
                    uint8_t lastflag = lastline[i + 7] & 0x80;
                    uint8_t firstflag = firstline[i + 8] & 0x1;

                    if (lastflag && firstflag)
                    {
                        uint64_t lowlimit = sdata->lowlimit + i * 8 * sdata->prodN;
                        uint64_t prime = lowlimit + 63 * sdata->prodN + sdata->rclass[sdata->numclasses - 1];

                        lowlimit = sdata->lowlimit + (i + 8) * 8 * sdata->prodN;
                        uint64_t p2 = lowlimit + 0 * sdata->prodN + sdata->rclass[0];

                        if ((prime >= sdata->orig_llimit) && (prime <= sdata->orig_hlimit))
                        {
                            t->ddata.primes[sdata->GLOBAL_OFFSET + t->linecount++] = prime;
                        }
                    }
                }
                else if ((i + 8) < sdata->numlinebytes)
                {
                    // this thread is done but if there is more data after
                    // this thread's chunk then check between thread boundaries.
                    uint8_t lastflag = lastline[i + 7] & 0x80;
                    uint8_t firstflag = firstline[i + 8] & 0x1;

                    if (lastflag && firstflag)
                    {
                        uint64_t lowlimit = sdata->lowlimit + i * 8 * sdata->prodN;
                        uint64_t prime = lowlimit + 63 * sdata->prodN + sdata->rclass[sdata->numclasses - 1];

                        lowlimit = sdata->lowlimit + (i + 8) * 8 * sdata->prodN;
                        uint64_t p2 = lowlimit + 0 * sdata->prodN + sdata->rclass[0];

                        if ((prime >= sdata->orig_llimit) && (prime <= sdata->orig_hlimit))
                        {
                            //printf("found twin %lu between threads\n", prime);
                            t->ddata.primes[sdata->GLOBAL_OFFSET + t->linecount++] = prime;
                        }
                    }
                }
            }
        }
    }
    else
    {
        // if we don't have BMI2, or numclasses > 48, or we're doing a more
        // complicated analysis, then we should end up here.
        for (i = t->startid; i < t->stopid; i += 8)
        {
            t->linecount = compute_8_bytes(sdata, t->linecount, t->ddata.primes, i);

            // if searching for prime constellations, then we need to look at the last 
            // flags of this block of 8 bytes and the first flags of the next one.
            // we do this by loading the trailing bits of this block into a carry
            // register.  The rest is handled by the block analysis function.
            if ((sdata->analysis == 2) && (sdata->is_main_sieve == 1))
            {
                uint8_t* lastline = sdata->lines[sdata->numclasses - 1];
                uint8_t* firstline = sdata->lines[0];

                if ((i + 8) < t->stopid)
                {
                    uint8_t lastflag = lastline[i + 7] & 0x80;
                    uint8_t firstflag = firstline[i + 8] & 0x1;

                    if (lastflag && firstflag)
                    {
                        uint64_t lowlimit = sdata->lowlimit + i * 8 * sdata->prodN;
                        uint64_t prime = lowlimit + 63 * sdata->prodN + sdata->rclass[sdata->numclasses - 1];

                        lowlimit = sdata->lowlimit + (i + 8) * 8 * sdata->prodN;
                        uint64_t p2 = lowlimit + 0 * sdata->prodN + sdata->rclass[0];
                        
                        if ((prime >= sdata->orig_llimit) && (prime <= sdata->orig_hlimit))
                        {
                            t->ddata.primes[sdata->GLOBAL_OFFSET + t->linecount++] = prime;
                        }
                    }
                }
                else if ((i + 8) < sdata->numlinebytes)
                {
                    // this thread is done but if there is more data after
                    // this thread's chunk then check between thread boundaries.
                    uint8_t lastflag = lastline[i + 7] & 0x80;
                    uint8_t firstflag = firstline[i + 8] & 0x1;

                    if (lastflag && firstflag)
                    {
                        uint64_t lowlimit = sdata->lowlimit + i * 8 * sdata->prodN;
                        uint64_t prime = lowlimit + 63 * sdata->prodN + sdata->rclass[sdata->numclasses - 1];

                        lowlimit = sdata->lowlimit + (i + 8) * 8 * sdata->prodN;
                        uint64_t p2 = lowlimit + 0 * sdata->prodN + sdata->rclass[0];

                        if ((prime >= sdata->orig_llimit) && (prime <= sdata->orig_hlimit))
                        {
                            //printf("found twin %lu between threads\n", prime);
                            t->ddata.primes[sdata->GLOBAL_OFFSET + t->linecount++] = prime;
                        }
                    }
                }
            }
        }
    }
#else
    for (i = t->startid; i < t->stopid; i += 8)
    {
        t->linecount = compute_8_bytes(sdata, t->linecount, t->ddata.primes, i);
    }
#endif

    return;
}

uint64_t primes_from_lineflags(soe_staticdata_t *sdata, thread_soedata_t *thread_data,
	uint32_t start_count, uint64_t *primes)
{
	// compute primes using all of the sieved lines we have stored
	uint64_t pcount = start_count;	
	uint64_t i;
	int j;
	uint32_t range, lastid;
    int GLOBAL_OFFSET = sdata->GLOBAL_OFFSET;

    // timing
    double t;
    struct timeval tstart, tstop;

    // threading structures
    tpool_t *tpool_data;
    soe_userdata_t udata;

    if (sdata->VFLAG > 1)
    {
        gettimeofday(&tstart, NULL);
    }

	// each thread needs to work on a number of bytes that is divisible by 32
	range = sdata->numlinebytes / sdata->THREADS;
	range -= (range % 32);
	lastid = 0;

    if (sdata->VFLAG > 1)
    {
        printf("computing primes from %lu to %lu, start_count is %u\n", 
            sdata->orig_llimit, sdata->orig_hlimit, start_count);
    }

    // divvy up the line bytes.  Unlike when counting primes,
    // here threading is by groups of bytes, over all classes.
    // this is necessary to order the primes.
    for (i = 0; i < sdata->THREADS; i++)
    {
        thread_soedata_t *t = thread_data + i;

        t->sdata = *sdata;
        t->startid = lastid;
        t->stopid = t->startid + range;
        lastid = t->stopid;

        if (sdata->VFLAG > 2)
        {
            printf("thread %d finding primes from byte offset %u to %u\n",
                (int)i, t->startid, t->stopid);
        }
    }

    // allocate a temporary array for each thread's primes
    if (sdata->THREADS > 1)
    {
        uint64_t memchunk;

        if (sdata->sieve_range)
        {
            // then just split the overall range into equal parts
            memchunk = (uint64_t)((double)(sdata->num_found / sdata->THREADS) * 1.1); 

            if (sdata->VFLAG > 2)
            {
                printf("allocating temporary space for %" PRIu64 " primes per thread\n",
                    memchunk);
            }

            for (i = 0; i < sdata->THREADS; i++)
            {
                thread_soedata_t *t = thread_data + i;

                t->ddata.primes = (uint64_t *)malloc(memchunk * sizeof(uint64_t));
            }
        }
        else
        {
            // then estimate the number of primes we'll find in each chunk.
            // it's important to do this chunk by chunk, because in some cases
            // the number of primes changes rapidly as a function of offset
            // from the start of the range (i.e., when start is 0)
            uint64_t hi_est, lo_est;
            uint64_t tmplo = sdata->orig_llimit;
            uint64_t tmphi;
            uint64_t chunk = 8 * sdata->numlinebytes / sdata->THREADS;

            chunk *= sdata->prodN;
            tmphi = tmplo + chunk;
            for (i = 0; i < sdata->THREADS; i++)
            {
                thread_soedata_t *t = thread_data + i;

                hi_est = (uint64_t)(tmphi / log((double)tmphi));
                if (tmplo > 1)
                    lo_est = (uint64_t)(tmplo / log((double)tmplo));
                else
                    lo_est = 0;

                memchunk = (uint64_t)((double)(hi_est - lo_est) * 1.25);

                if (sdata->VFLAG > 2)
                {
                    printf("allocating temporary space for %" PRIu64 " primes between %" PRIu64 " and %" PRIu64 "\n",
                        memchunk, tmplo, tmphi);
                }

                t->ddata.primes = (uint64_t *)malloc(memchunk * sizeof(uint64_t));

                tmplo += chunk;
                tmphi += chunk;
            }
        }
    }
    else
    {
        // with just one thread, don't bother with creating a temporary array
        thread_data[0].ddata.primes = primes;
    }

    udata.sdata = sdata;
    udata.ddata = thread_data;
    tpool_data = tpool_setup(sdata->THREADS, NULL, NULL, NULL,
        &compute_primes_dispatch, &udata);

    if (sdata->THREADS == 1)
    {
        thread_data->linecount = pcount;
        compute_primes_work_fcn(tpool_data);
    }
    else
    {
        sdata->sync_count = 0;
        tpool_add_work_fcn(tpool_data, &compute_primes_work_fcn);
        tpool_go(tpool_data);
    }
    free(tpool_data);    

	// now combine all of the temporary arrays, if necessary
	if (sdata->THREADS > 1)
	{
		pcount = start_count;
		for (j = 0; j < sdata->THREADS; j++)
		{
			thread_soedata_t *t = thread_data + j;

			if (t->linecount == 0)
			{
				free(t->ddata.primes);
				continue;
			}

            if (sdata->VFLAG > 2)
            {
                printf("adding %" PRIu64 " primes found in thread %d\n", t->linecount, j);
            }

			memcpy(primes + GLOBAL_OFFSET + pcount, t->ddata.primes, t->linecount * sizeof(uint64_t));

			pcount += t->linecount;
			free(t->ddata.primes);
		}
	}
    else
    {
        pcount = thread_data[0].linecount;
    }

	// and finally, get primes from any residual portion of the line arrays
	// using a direct method
	if (lastid != sdata->numlinebytes)
	{
        if (sdata->VFLAG > 2)
        {
            printf("adding primes from byte offset %u to %u\n",
                lastid, (uint32_t)sdata->numlinebytes);
        }

		for (i = lastid; i < sdata->numlinebytes; i+=8)
		{
			pcount = compute_8_bytes(sdata, pcount, primes, i);		
		}
	}

    if (sdata->VFLAG > 1)
    {
        gettimeofday(&tstop, NULL);

        t = ytools_difftime(&tstart, &tstop);

        if (sdata->VFLAG > 2)
        {
            printf("time to compute primes = %1.4f\n", t);
        }
    }

	return pcount;
}

uint32_t compute_8_bytes(soe_staticdata_t *sdata, 
	uint32_t pcount, uint64_t *primes, uint64_t byte_offset)
{
	uint32_t current_line;
	// re-ordering queues supporting arbitrary residue classes.
    uint64_t **pqueues; // [64][48]
    uint32_t pcounts[64];
    int i, j;
	uint32_t nc = sdata->numclasses;
	uint64_t lowlimit = sdata->lowlimit;
	uint64_t prodN = sdata->prodN;
	uint8_t **lines = sdata->lines;
	uint64_t olow = sdata->orig_llimit;
	uint64_t ohigh = sdata->orig_hlimit;
    int GLOBAL_OFFSET = sdata->GLOBAL_OFFSET;
		
	if ((byte_offset & 32767) == 0)
	{
		if (sdata->VFLAG > 1)
		{
			printf("computing: %d%%\r",(int)
				((double)byte_offset / (double)(sdata->numlinebytes) * 100.0));
			fflush(stdout);
		}
	}

    pqueues = (uint64_t**)xmalloc(64 * sizeof(uint64_t*));
    for (i = 0; i < 64; i++)
    {
        pqueues[i] = (uint64_t*)xmalloc(sdata->numclasses * sizeof(uint64_t));
    }

    // Compute the primes using ctz on the 64-bit words but push the results
    // into 64 different queues depending on the bit position.  Then
    // we pull from the queues in order while storing into the primes array.
    // This time the bottleneck is mostly in the queue-based sorting
    // and associated memory operations, so we don't bother with
    // switching between branch-free inner loops or not.
    memset(pcounts, 0, 64 * sizeof(uint32_t));

    lowlimit += byte_offset * 8 * prodN;
    for (current_line = 0; current_line < nc; current_line++)
    {
        uint64_t *line64 = (uint64_t *)lines[current_line];
        uint64_t flags64 = line64[byte_offset / 8];

        //if ((sdata->analysis == 2) && (sdata->is_main_sieve == 1))
        //    printf("sorting bits in twin residue %u\n", sdata->rclass[current_line]);

        while (flags64 > 0)
        {
            uint64_t pos = _trail_zcnt64(flags64);
            uint64_t prime = lowlimit + pos * prodN + sdata->rclass[current_line];

            if ((prime >= olow) && (prime <= ohigh))
            {
                pqueues[pos][pcounts[pos]] = prime;
                pcounts[pos]++;
            }
            flags64 ^= (1ULL << pos);
        }
    }
    //if ((sdata->analysis == 2) && (sdata->is_main_sieve == 1))
    //    exit(1);

    for (i = 0; i < 64; i++)
    {
        if ((sdata->analysis == 1) || (sdata->is_main_sieve == 0))
        {
            // now we can traverse the queues in order and
            // copy into the main prime array.
            for (j = 0; j < pcounts[i] / 2; j++)
            {
                __m128i t = _mm_loadu_si128((__m128i*)(&pqueues[i][j * 2]));
                _mm_storeu_si128((__m128i*)(&primes[GLOBAL_OFFSET + pcount]), t);
                pcount += 2;
            }
            for (j *= 2; j < pcounts[i]; j++)
            {
                primes[GLOBAL_OFFSET + pcount++] = pqueues[i][j];
            }
        }
        else if (sdata->analysis == 2)
        {
            // search for twins and only load the leading element/prime.
            // if depth-based sieving then these are candidate twins.
            if (pcounts[i] > 0)
            {
                for (j = 0; j < pcounts[i] - 1; j++)
                {
                    if ((pqueues[i][j + 1] - pqueues[i][j]) == 2)
                    {
                        primes[GLOBAL_OFFSET + pcount++] = pqueues[i][j];
                    }
                }
                if (i < 63)
                {
                    if (pcounts[i + 1] > 0)
                    {
                        if ((pqueues[i + 1][0] - pqueues[i][j]) == 2)
                        {
                            primes[GLOBAL_OFFSET + pcount++] = pqueues[i][j];
                        }
                    }
                }
            }
        }
        else if (sdata->analysis == 3)
        {
            // search for gaps and only load the leading element/prime.
            // if depth-based sieving then these are candidate gaps.
            for (j = 0; j < pcounts[i] / 2; j++)
            {
                __m128i t = _mm_loadu_si128((__m128i*)(&pqueues[i][j * 2]));
                _mm_storeu_si128((__m128i*)(&primes[GLOBAL_OFFSET + pcount]), t);
                pcount += 2;
            }
            for (j *= 2; j < pcounts[i]; j++)
            {
                primes[GLOBAL_OFFSET + pcount++] = pqueues[i][j];
            }
        }
    }

    for (i = 0; i < 64; i++)
    {
        free(pqueues[i]);
    }
    free(pqueues);

	return pcount;
}


#if defined(USE_BMI2) || defined(USE_AVX512F)

__inline uint64_t interleave_avx2_bmi2_pdep2x32(uint32_t x1, uint32_t x2)
{
    return _pdep_u64(x1, 0x5555555555555555) 
        | _pdep_u64(x2, 0xaaaaaaaaaaaaaaaa);
}

__inline uint64_t interleave_avx2_bmi2_pdep(uint8_t x1,
    uint8_t x2,
    uint8_t x3,
    uint8_t x4,
    uint8_t x5,
    uint8_t x6,
    uint8_t x7,
    uint8_t x8)
{
    return _pdep_u64(x1, 0x0101010101010101ull) |
        _pdep_u64(x2, 0x0202020202020202ull) |
        _pdep_u64(x3, 0x0404040404040404ull) |
        _pdep_u64(x4, 0x0808080808080808ull) |
        _pdep_u64(x5, 0x1010101010101010ull) |
        _pdep_u64(x6, 0x2020202020202020ull) |
        _pdep_u64(x7, 0x4040404040404040ull) |
        _pdep_u64(x8, 0x8080808080808080ull);
}

uint32_t compute_8_bytes_bmi2(soe_staticdata_t *sdata,
    uint32_t pcount, uint64_t *primes, uint64_t byte_offset)
{    
    uint32_t nc = sdata->numclasses;
    uint64_t lowlimit = sdata->lowlimit;
    uint8_t **lines = sdata->lines;
    uint64_t olow = sdata->orig_llimit;
    uint64_t ohigh = sdata->orig_hlimit;
    int GLOBAL_OFFSET = sdata->GLOBAL_OFFSET;

    if ((byte_offset & 32767) == 0)
    {
        if (sdata->VFLAG > 1)
        {
            printf("computing: %d%%\r", (int)
                ((double)byte_offset / (double)(sdata->numlinebytes) * 100.0));
            fflush(stdout);
        }
    }

    // BMI2 version, new instructions help quite a bit:
    // use _pdep_u64 to align/interleave bits from multiple bytes, 
    // _blsr_u64 to clear the last set bit, and depending on the 
    // number of residue classes, AVX2 vector load/store operations.

    // compute the minimum/maximum prime we could encounter in this range
        // and execute either a branch-free innermost loop or not.
    uint64_t plow = (byte_offset + 0) * 8 * sdata->prodN + 0 * sdata->prodN +
        sdata->rclass[0] + lowlimit;

    if (plow > ohigh)
        return pcount;

    // here is the 2 line version
    if (nc == 2)
    {
        int i;
        uint32_t last_bit = 0;
        uint32_t *lines32a = (uint32_t *)lines[0];
        uint32_t *lines32b = (uint32_t *)lines[1];

        // align the current bytes in next 2 residue classes
        lowlimit += byte_offset * 8 * sdata->prodN;
        for (i = 0; i < 2; i++)
        {
            uint64_t aligned_flags;

            aligned_flags = interleave_avx2_bmi2_pdep2x32(
                lines32a[byte_offset/4+i],
                lines32b[byte_offset/4+i]);

            if ((sdata->analysis == 2) && (sdata->is_main_sieve))
            {
                // alternate bits encode potential primes
                // in residue classes 1 and 5.  So twins can 
                // only exist with flags in class 5 followed by 1.
                uint64_t twins = aligned_flags & (aligned_flags >> 1);

                twins &= 0xaaaaaaaaaaaaaaaaULL;

                //pcount += _mm_popcnt_u64(twins);
                if (last_bit & aligned_flags)
                {
                    uint64_t prime = lowlimit - 1;
                    primes[GLOBAL_OFFSET + pcount++] = prime;
                }
                while (twins > 0)
                {
                    uint64_t pos = _trail_zcnt64(twins);
                    uint64_t prime = lowlimit + (pos / 2) * 6 + sdata->rclass[pos % 2];

                    primes[GLOBAL_OFFSET + pcount++] = prime;
                    twins = _reset_lsb64(twins);
                }
                lowlimit += 32 * sdata->prodN;
                last_bit = (aligned_flags >> 63);
            }
            else
            {
                // then compute primes in order for flags that are set.
                while (aligned_flags > 0)
                {
                    uint64_t pos = _trail_zcnt64(aligned_flags);
                    uint64_t prime = lowlimit + (pos / 2) * 6 + sdata->rclass[pos % 2];

                    primes[GLOBAL_OFFSET + pcount++] = prime;
                    aligned_flags = _reset_lsb64(aligned_flags);
                }
                lowlimit += 32 * sdata->prodN;
            }
        }

    }
    else if (nc == 8)
    {
        int i;
        uint32_t last_bit = 0;

        // align the current bytes in next 8 residue classes
        lowlimit += byte_offset * 8 * sdata->prodN;
        for (i = 0; i < 8; i++)
        {
            uint64_t aligned_flags;

            aligned_flags = interleave_avx2_bmi2_pdep(lines[0][byte_offset+i],
                lines[1][byte_offset+i],
                lines[2][byte_offset+i],
                lines[3][byte_offset+i],
                lines[4][byte_offset+i],
                lines[5][byte_offset+i],
                lines[6][byte_offset+i],
                lines[7][byte_offset+i]);

            if ((sdata->analysis == 2) && (sdata->is_main_sieve))
            {
                uint64_t twins = aligned_flags & (aligned_flags >> 1);
                twins &= 0x9494949494949494ULL;

                //pcount += _mm_popcnt_u64(twins);
                if (last_bit & aligned_flags)
                {
                    uint64_t prime = lowlimit - 1;
                    primes[GLOBAL_OFFSET + pcount++] = prime;
                }
                while (twins > 0)
                {
                    uint64_t pos = _trail_zcnt64(twins);
                    uint64_t prime = lowlimit + (pos / 8) * 30 + sdata->rclass[pos % 8];

                    primes[GLOBAL_OFFSET + pcount++] = prime;
                    twins = _reset_lsb64(twins);
                }
                lowlimit += 8 * sdata->prodN;
                last_bit = (aligned_flags >> 63);
            }
            else
            {
                // then compute primes in order for flags that are set.
                while (aligned_flags > 0)
                {
                    uint64_t pos = _trail_zcnt64(aligned_flags);
                    uint64_t prime = lowlimit + (pos / 8) * 30 + sdata->rclass[pos % 8];

                    primes[GLOBAL_OFFSET + pcount++] = prime;
                    aligned_flags = _reset_lsb64(aligned_flags);
                }
                lowlimit += 8 * sdata->prodN;
            }
        }

    }
    else if (nc == 30)
    {
        // twins residues for mod 210
        int i;
        uint32_t last_bit = 0;

        lowlimit += byte_offset * 8 * sdata->prodN;
        for (i = 0; i < 8; i++)
        {
            uint64_t aligned_flags1;
            uint64_t aligned_flags2;
            uint64_t aligned_flags3;
            uint64_t aligned_flags4;

            // first we partially order the lines such
            // that each 64-bit flags contains 8 ordered
            // bytes for a set of 8 classes.
            aligned_flags1 = interleave_avx2_bmi2_pdep(
                lines[0][byte_offset + i],
                lines[1][byte_offset + i],
                lines[2][byte_offset + i],
                lines[3][byte_offset + i],
                lines[4][byte_offset + i],
                lines[5][byte_offset + i],
                lines[6][byte_offset + i],
                lines[7][byte_offset + i]);

            aligned_flags2 = interleave_avx2_bmi2_pdep(
                lines[8 + 0][byte_offset + i],
                lines[8 + 1][byte_offset + i],
                lines[8 + 2][byte_offset + i],
                lines[8 + 3][byte_offset + i],
                lines[8 + 4][byte_offset + i],
                lines[8 + 5][byte_offset + i],
                lines[8 + 6][byte_offset + i],
                lines[8 + 7][byte_offset + i]);

            aligned_flags3 = interleave_avx2_bmi2_pdep(
                lines[16 + 0][byte_offset + i],
                lines[16 + 1][byte_offset + i],
                lines[16 + 2][byte_offset + i],
                lines[16 + 3][byte_offset + i],
                lines[16 + 4][byte_offset + i],
                lines[16 + 5][byte_offset + i],
                lines[16 + 6][byte_offset + i],
                lines[16 + 7][byte_offset + i]);

            aligned_flags4 = interleave_avx2_bmi2_pdep(
                lines[24 + 0][byte_offset + i],
                lines[24 + 1][byte_offset + i],
                lines[24 + 2][byte_offset + i],
                lines[24 + 3][byte_offset + i],
                lines[24 + 4][byte_offset + i],
                lines[24 + 5][byte_offset + i],
                0,
                0);


            // aligned_flags1 contains: b0(c0-7), b1(c0-7), b2(c0-7), ... b7(c0-7)
            // aligned_flags2 contains: b0(c8-15), b1(c8-15), b2(c8-15), ... b7(c8-15)
            // aligned_flags3 contains: b0(c16-23), b1(c16-23), b2(c16-23), ... b7(c16-23)
            // aligned_flags4 contains: b0(c24-29), b1(c24-29), b2(c24-29), ... b7(c24-29)

            // now, shuffle the bytes within the partially ordered chunks.
            uint8_t unordered_bytes[48];
            uint8_t ordered_bytes[64];
            uint64_t* unordered64 = (uint64_t*)unordered_bytes;
            uint64_t* ordered64 = (uint64_t*)ordered_bytes;
            uint32_t* ordered32 = (uint32_t*)ordered_bytes;
            unordered64[0] = aligned_flags1;
            unordered64[1] = aligned_flags2;
            unordered64[2] = aligned_flags3;
            unordered64[3] = aligned_flags4;

            ordered_bytes[0] = unordered_bytes[0];
            ordered_bytes[1] = unordered_bytes[8];
            ordered_bytes[2] = unordered_bytes[16];
            ordered_bytes[3] = unordered_bytes[24];
            ordered_bytes[4] = unordered_bytes[1];
            ordered_bytes[5] = unordered_bytes[9];
            ordered_bytes[6] = unordered_bytes[17];
            ordered_bytes[7] = unordered_bytes[25];

            ordered_bytes[8] = unordered_bytes[2];
            ordered_bytes[9] = unordered_bytes[10];
            ordered_bytes[10] = unordered_bytes[18];
            ordered_bytes[11] = unordered_bytes[26];
            ordered_bytes[12] = unordered_bytes[3];
            ordered_bytes[13] = unordered_bytes[11];
            ordered_bytes[14] = unordered_bytes[19];
            ordered_bytes[15] = unordered_bytes[27];

            ordered_bytes[16] = unordered_bytes[4];
            ordered_bytes[17] = unordered_bytes[12];
            ordered_bytes[18] = unordered_bytes[20];
            ordered_bytes[19] = unordered_bytes[28];
            ordered_bytes[20] = unordered_bytes[5];
            ordered_bytes[21] = unordered_bytes[13];
            ordered_bytes[22] = unordered_bytes[21];
            ordered_bytes[23] = unordered_bytes[29];

            ordered_bytes[24] = unordered_bytes[6];
            ordered_bytes[25] = unordered_bytes[14];
            ordered_bytes[26] = unordered_bytes[22];
            ordered_bytes[27] = unordered_bytes[30];
            ordered_bytes[28] = unordered_bytes[7];
            ordered_bytes[29] = unordered_bytes[15];
            ordered_bytes[30] = unordered_bytes[23];
            ordered_bytes[31] = unordered_bytes[31];

            // rotate 2 bits for continuity between
            // the two 30-bit sections.
            ordered32[0] |= ((ordered32[1] & 0x3) << 30);
            ordered32[2] |= ((ordered32[3] & 0x3) << 30);
            ordered32[4] |= ((ordered32[5] & 0x3) << 30);
            ordered32[6] |= ((ordered32[7] & 0x3) << 30);
            ordered32[1] >>= 2;
            ordered32[3] >>= 2;
            ordered32[5] >>= 2;
            ordered32[7] >>= 2;

            aligned_flags1 = ordered64[0];
            aligned_flags2 = ordered64[1];
            aligned_flags3 = ordered64[2];
            aligned_flags4 = ordered64[3];

            // twin residues mod 210:
            // 1 11 13 17 19 29 31 41 43 59 61 71 73 101 103 107 109 137 139 149 151 167 169 179 181 191 193 197 199 209
            // binary mask: 0101 0101 0101 0101 0101 0101 0101 01|01 0101 0101 0101 0101 0101 0101 0100
            // hex mask: AAAAAAAAAAAAAA20

            // then compute twins in order for flags that are set.
            // word 1
            uint64_t twins = aligned_flags1 & (aligned_flags1 >> 1);
            twins &= 0x02AAAAAAAAAAAAAAull;
            if (last_bit & aligned_flags1)
            {
                uint64_t prime = lowlimit - 1;

                primes[GLOBAL_OFFSET + pcount++] = prime;
                twins = twins & 0xfffffffffffffffeull;
            }
            while ((twins & 0x3fffffffull) > 0)
            {
                uint64_t pos = _trail_zcnt64((twins & 0x3fffffffull));
                uint64_t prime = lowlimit + sdata->rclass[pos];

                primes[GLOBAL_OFFSET + pcount++] = prime;
                twins = _reset_lsb64(twins);
            }
            lowlimit += sdata->prodN;
            while ((twins >> 30) > 0)
            {
                uint64_t pos = _trail_zcnt64((twins >> 30));
                uint64_t prime = lowlimit + sdata->rclass[pos];

                primes[GLOBAL_OFFSET + pcount++] = prime;
                twins = _reset_lsb64(twins);
            }
            lowlimit += sdata->prodN;
            last_bit = (aligned_flags1 >> 59);      // generate carry

            // 2
            twins = aligned_flags2 & (aligned_flags2 >> 1);
            twins &= 0x02AAAAAAAAAAAAAAull;
            if (last_bit & aligned_flags2)
            {
                uint64_t prime = lowlimit - 1;

                primes[GLOBAL_OFFSET + pcount++] = prime;
                twins = twins & 0xfffffffffffffffeull;
            }
            while ((twins & 0x3fffffffull) > 0)
            {
                uint64_t pos = _trail_zcnt64((twins & 0x3fffffffull));
                uint64_t prime = lowlimit + sdata->rclass[pos];

                primes[GLOBAL_OFFSET + pcount++] = prime;
                twins = _reset_lsb64(twins);
            }
            lowlimit += sdata->prodN;
            while ((twins >> 30) > 0)
            {
                uint64_t pos = _trail_zcnt64((twins >> 30));
                uint64_t prime = lowlimit + sdata->rclass[pos];

                primes[GLOBAL_OFFSET + pcount++] = prime;
                twins = _reset_lsb64(twins);
            }
            lowlimit += sdata->prodN;
            last_bit = (aligned_flags2 >> 59);      // generate carry

            // 3
            twins = aligned_flags3 & (aligned_flags3 >> 1);
            twins &= 0x02AAAAAAAAAAAAAAull;
            if (last_bit & aligned_flags3)
            {
                uint64_t prime = lowlimit - 1;

                primes[GLOBAL_OFFSET + pcount++] = prime;
                twins = twins & 0xfffffffffffffffeull;
            }
            while ((twins & 0x3fffffffull) > 0)
            {
                uint64_t pos = _trail_zcnt64((twins & 0x3fffffffull));
                uint64_t prime = lowlimit + sdata->rclass[pos];

                primes[GLOBAL_OFFSET + pcount++] = prime;
                twins = _reset_lsb64(twins);
            }
            lowlimit += sdata->prodN;
            while ((twins >> 30) > 0)
            {
                uint64_t pos = _trail_zcnt64((twins >> 30));
                uint64_t prime = lowlimit + sdata->rclass[pos];

                primes[GLOBAL_OFFSET + pcount++] = prime;
                twins = _reset_lsb64(twins);
            }
            lowlimit += sdata->prodN;
            last_bit = (aligned_flags3 >> 59);      // generate carry

            // 4
            twins = aligned_flags4 & (aligned_flags4 >> 1);
            twins &= 0x02AAAAAAAAAAAAAAull;
            if (last_bit & aligned_flags4)
            {
                uint64_t prime = lowlimit - 1;

                primes[GLOBAL_OFFSET + pcount++] = prime;
                twins = twins & 0xfffffffffffffffeull;
            }
            while ((twins & 0x3fffffffull) > 0)
            {
                uint64_t pos = _trail_zcnt64((twins & 0x3fffffffull));
                uint64_t prime = lowlimit + sdata->rclass[pos];

                primes[GLOBAL_OFFSET + pcount++] = prime;
                twins = _reset_lsb64(twins);
            }
            lowlimit += sdata->prodN;
            while ((twins >> 30) > 0)
            {
                uint64_t pos = _trail_zcnt64((twins >> 30));
                uint64_t prime = lowlimit + sdata->rclass[pos];

                primes[GLOBAL_OFFSET + pcount++] = prime;
                twins = _reset_lsb64(twins);
            }
            lowlimit += sdata->prodN;
            last_bit = (aligned_flags4 >> 59);      // generate carry
        }
    }
    else if (nc == 48)
    {
        int i;

        lowlimit += byte_offset * 8 * sdata->prodN;
        for (i = 0; i < 8; i++)
        {
            uint64_t aligned_flags1;
            uint64_t aligned_flags2;
            uint64_t aligned_flags3;
            uint64_t aligned_flags4;
            uint64_t aligned_flags5;
            uint64_t aligned_flags6;

            aligned_flags1 = interleave_avx2_bmi2_pdep(
                lines[0][byte_offset + i],
                lines[1][byte_offset + i],
                lines[2][byte_offset + i],
                lines[3][byte_offset + i],
                lines[4][byte_offset + i],
                lines[5][byte_offset + i],
                lines[6][byte_offset + i],
                lines[7][byte_offset + i]);

            aligned_flags2 = interleave_avx2_bmi2_pdep(
                lines[8+0][byte_offset + i],
                lines[8+1][byte_offset + i],
                lines[8+2][byte_offset + i],
                lines[8+3][byte_offset + i],
                lines[8+4][byte_offset + i],
                lines[8+5][byte_offset + i],
                lines[8+6][byte_offset + i],
                lines[8+7][byte_offset + i]);

            aligned_flags3 = interleave_avx2_bmi2_pdep(
                lines[16+0][byte_offset + i],
                lines[16+1][byte_offset + i],
                lines[16+2][byte_offset + i],
                lines[16+3][byte_offset + i],
                lines[16+4][byte_offset + i],
                lines[16+5][byte_offset + i],
                lines[16+6][byte_offset + i],
                lines[16+7][byte_offset + i]);

            aligned_flags4 = interleave_avx2_bmi2_pdep(
                lines[24+0][byte_offset + i],
                lines[24+1][byte_offset + i],
                lines[24+2][byte_offset + i],
                lines[24+3][byte_offset + i],
                lines[24+4][byte_offset + i],
                lines[24+5][byte_offset + i],
                lines[24+6][byte_offset + i],
                lines[24+7][byte_offset + i]);

            aligned_flags5 = interleave_avx2_bmi2_pdep(
                lines[32+0][byte_offset + i],
                lines[32+1][byte_offset + i],
                lines[32+2][byte_offset + i],
                lines[32+3][byte_offset + i],
                lines[32+4][byte_offset + i],
                lines[32+5][byte_offset + i],
                lines[32+6][byte_offset + i],
                lines[32+7][byte_offset + i]);

            aligned_flags6 = interleave_avx2_bmi2_pdep(
                lines[40+0][byte_offset + i],
                lines[40+1][byte_offset + i],
                lines[40+2][byte_offset + i],
                lines[40+3][byte_offset + i],
                lines[40+4][byte_offset + i],
                lines[40+5][byte_offset + i],
                lines[40+6][byte_offset + i],
                lines[40+7][byte_offset + i]);

            // aligned_flags1 contains: b0(c0-7), b1(c0-7), b2(c0-7), ... b7(c0-7)
            // aligned_flags2 contains: b0(c8-15), b1(c8-15), b2(c8-15), ... b7(c8-15)
            // so, each byte is ordered and we just need to rearrage the bytes.
            uint8_t unordered_bytes[48];
            uint8_t ordered_bytes[48];
            uint64_t* unordered64 = (uint64_t*)unordered_bytes;
            uint64_t* ordered64 = (uint64_t*)ordered_bytes;
            unordered64[0] = aligned_flags1;
            unordered64[1] = aligned_flags2;
            unordered64[2] = aligned_flags3;
            unordered64[3] = aligned_flags4;
            unordered64[4] = aligned_flags5;
            unordered64[5] = aligned_flags6;

            ordered_bytes[00] = unordered_bytes[00];
            ordered_bytes[01] = unordered_bytes[8];
            ordered_bytes[02] = unordered_bytes[16];
            ordered_bytes[03] = unordered_bytes[24];
            ordered_bytes[04] = unordered_bytes[32];
            ordered_bytes[05] = unordered_bytes[40];
            ordered_bytes[06] = unordered_bytes[01];
            ordered_bytes[07] = unordered_bytes[9];
            ordered_bytes[8] = unordered_bytes[17];
            ordered_bytes[9] = unordered_bytes[25];
            ordered_bytes[10] = unordered_bytes[33];
            ordered_bytes[11] = unordered_bytes[41];
            ordered_bytes[12] = unordered_bytes[02];
            ordered_bytes[13] = unordered_bytes[10];
            ordered_bytes[14] = unordered_bytes[18];
            ordered_bytes[15] = unordered_bytes[26];
            ordered_bytes[16] = unordered_bytes[34];
            ordered_bytes[17] = unordered_bytes[42];
            ordered_bytes[18] = unordered_bytes[03];
            ordered_bytes[19] = unordered_bytes[11];
            ordered_bytes[20] = unordered_bytes[19];
            ordered_bytes[21] = unordered_bytes[27];
            ordered_bytes[22] = unordered_bytes[35];
            ordered_bytes[23] = unordered_bytes[43];
            ordered_bytes[24] = unordered_bytes[04];
            ordered_bytes[25] = unordered_bytes[12];
            ordered_bytes[26] = unordered_bytes[20];
            ordered_bytes[27] = unordered_bytes[28];
            ordered_bytes[28] = unordered_bytes[36];
            ordered_bytes[29] = unordered_bytes[44];
            ordered_bytes[30] = unordered_bytes[05];
            ordered_bytes[31] = unordered_bytes[13];
            ordered_bytes[32] = unordered_bytes[21];
            ordered_bytes[33] = unordered_bytes[29];
            ordered_bytes[34] = unordered_bytes[37];
            ordered_bytes[35] = unordered_bytes[45];
            ordered_bytes[36] = unordered_bytes[06];
            ordered_bytes[37] = unordered_bytes[14];
            ordered_bytes[38] = unordered_bytes[22];
            ordered_bytes[39] = unordered_bytes[30];
            ordered_bytes[40] = unordered_bytes[38];
            ordered_bytes[41] = unordered_bytes[46];
            ordered_bytes[42] = unordered_bytes[07];
            ordered_bytes[43] = unordered_bytes[15];
            ordered_bytes[44] = unordered_bytes[23];
            ordered_bytes[45] = unordered_bytes[31];
            ordered_bytes[46] = unordered_bytes[39];
            ordered_bytes[47] = unordered_bytes[47];

            aligned_flags1 = ordered64[0];
            aligned_flags2 = ordered64[1];
            aligned_flags3 = ordered64[2];
            aligned_flags4 = ordered64[3];
            aligned_flags5 = ordered64[4];
            aligned_flags6 = ordered64[5];

            // then compute primes in order for flags that are set.
            while (aligned_flags1 > 0)
            {
                uint64_t pos = _trail_zcnt64(aligned_flags1);
                uint64_t prime = lowlimit + (pos / 48) * 210 + sdata->rclass[pos % 48];

                primes[GLOBAL_OFFSET + pcount++] = prime;
                aligned_flags1 = _reset_lsb64(aligned_flags1);
            }
            // then compute primes in order for flags that are set.
            while (aligned_flags2 > 0)
            {
                uint64_t pos = _trail_zcnt64(aligned_flags2);
                uint64_t prime = lowlimit + ((pos + 64) / 48) * 210 + sdata->rclass[(pos + 64) % 48];

                primes[GLOBAL_OFFSET + pcount++] = prime;
                aligned_flags2 = _reset_lsb64(aligned_flags2);
            }
            // then compute primes in order for flags that are set.
            while (aligned_flags3 > 0)
            {
                uint64_t pos = _trail_zcnt64(aligned_flags3);
                uint64_t prime = lowlimit + ((pos + 128) / 48) * 210 + sdata->rclass[(pos + 128) % 48];

                primes[GLOBAL_OFFSET + pcount++] = prime;
                aligned_flags3 = _reset_lsb64(aligned_flags3);
            }
            // then compute primes in order for flags that are set.
            while (aligned_flags4 > 0)
            {
                uint64_t pos = _trail_zcnt64(aligned_flags4);
                uint64_t prime = lowlimit + ((pos + 192) / 48) * 210 + sdata->rclass[(pos + 192) % 48];

                primes[GLOBAL_OFFSET + pcount++] = prime;
                aligned_flags4 = _reset_lsb64(aligned_flags4);
            }
            // then compute primes in order for flags that are set.
            while (aligned_flags5 > 0)
            {
                uint64_t pos = _trail_zcnt64(aligned_flags5);
                uint64_t prime = lowlimit + ((pos + 256) / 48) * 210 + sdata->rclass[(pos + 256) % 48];

                primes[GLOBAL_OFFSET + pcount++] = prime;
                aligned_flags5 = _reset_lsb64(aligned_flags5);
            }
            // then compute primes in order for flags that are set.
            while (aligned_flags6 > 0)
            {
                uint64_t pos = _trail_zcnt64(aligned_flags6);
                uint64_t prime = lowlimit + ((pos + 320) / 48) * 210 + sdata->rclass[(pos + 320) % 48];

                primes[GLOBAL_OFFSET + pcount++] = prime;
                aligned_flags6 = _reset_lsb64(aligned_flags6);
            }
            lowlimit += 8 * sdata->prodN;
        }
    }
    else if (nc == 270)
    {
        // twins residues for mod 2310

        /*

        twin residues mod 2310 (separated into groups of 64):
        1 17 19 29 31 41 43 59 61 71 73 101 103 107 109 137 139 149 151 167 169 179 181 191 193
        197 199 221 223 227 229 239 241 269 271 281 283 311 313 347 349 359 361 377 379 389 391
        401 403 419 421 431 433 437 439 461 463 479 481 491 493 521 523 527

        529 557 559 569 571 587 589 599 601 611 613 617 619 629 631 641 643 659 661 689 691 701
        703 731 733 767 769 797 799 809 811 821 823 827 829 839 841 851 853 857 859 881 883 899
        901 941 943 947 949 989 991 1007 1009 1019 1021 1031 1033 1037 1039 1049 1051 1061 1063 1079

        1081 1091 1093 1121 1123 1151 1153 1157 1159 1187 1189 1217 1219 1229 1231 1247 1249 1259
        1261 1271 1273 1277 1279 1289 1291 1301 1303 1319 1321 1361 1363 1367 1369 1409 1411 1427
        1429 1451 1453 1457 1459 1469 1471 1481 1483 1487 1489 1499 1501 1511 1513 1541 1543 1577
        1579 1607 1609 1619 1621 1649 1651 1667 1669 1679

        1681 1691 1693 1697 1699 1709 1711 1721 1723 1739 1741 1751 1753 1781 1783 1787 1789 1817
        1819 1829 1831 1847 1849 1871 1873 1877 1879 1889 1891 1907 1909 1919 1921 1931 1933 1949
        1951 1961 1963 1997 1999 2027 2029 2039 2041 2069 2071 2081 2083 2087 2089 2111 2113 2117
        2119 2129 2131 2141 2143 2159 2161 2171 2173 2201

        2203 2207 2209 2237 2239 2249 2251 2267 2269 2279 2281 2291 2293 2309

        */

        int i;
        uint32_t last_bit = 0;

        // sort this column of 270x64 bits
        for (i = 0; i < 8; i++)
        {
            // take one byte-column at a time (270x8 bits)
            uint64_t aligned_flags[34];
            int j;

            // first we partially order the lines such
            // that each 64-bit flags contains 8 ordered
            // bytes for a set of 8 classes.
            for (j = 0; j < 33; j++)
            {
                aligned_flags[j] = interleave_avx2_bmi2_pdep(
                    lines[(j * 8) + 0][byte_offset + i],
                    lines[(j * 8) + 1][byte_offset + i],
                    lines[(j * 8) + 2][byte_offset + i],
                    lines[(j * 8) + 3][byte_offset + i],
                    lines[(j * 8) + 4][byte_offset + i],
                    lines[(j * 8) + 5][byte_offset + i],
                    lines[(j * 8) + 6][byte_offset + i],
                    lines[(j * 8) + 7][byte_offset + i]);
            }

            aligned_flags[j] = interleave_avx2_bmi2_pdep(
                lines[(j * 8) + 0][byte_offset + i],
                lines[(j * 8) + 1][byte_offset + i],
                lines[(j * 8) + 2][byte_offset + i],
                lines[(j * 8) + 3][byte_offset + i],
                lines[(j * 8) + 4][byte_offset + i],
                lines[(j * 8) + 5][byte_offset + i],
                0,
                0);

            // now, shuffle the bytes within the partially ordered chunks.
            // every 270 bits there will be a discontinuity between columns
            // of ordered bits.  Rather than do the somewhat expensive
            // rotation like in the 30-class case we just group each 270-bit
            // column into its own set of 5 64-bit variables and then
            // carry propagate between the sets for each of the 8 columns 
            // in this byte_offset.  5 * 8 = 40 bytes per set. 
            uint8_t* unordered_bytes = (uint8_t*)aligned_flags;
            uint8_t ordered_bytes[320];
            uint64_t* ordered64 = (uint64_t*)ordered_bytes;

            for (j = 0; j < 8; j++)
            {
                ordered_bytes[(j * 40) + 0] = unordered_bytes[j + (0 * 8)];
                ordered_bytes[(j * 40) + 1] = unordered_bytes[j + (1 * 8)];
                ordered_bytes[(j * 40) + 2] = unordered_bytes[j + (2 * 8)];
                ordered_bytes[(j * 40) + 3] = unordered_bytes[j + (3 * 8)];
                ordered_bytes[(j * 40) + 4] = unordered_bytes[j + (4 * 8)];
                ordered_bytes[(j * 40) + 5] = unordered_bytes[j + (5 * 8)];
                ordered_bytes[(j * 40) + 6] = unordered_bytes[j + (6 * 8)];
                ordered_bytes[(j * 40) + 7] = unordered_bytes[j + (7 * 8)];
                ordered_bytes[(j * 40) + 8] = unordered_bytes[j + (8 * 8)];
                ordered_bytes[(j * 40) + 9] = unordered_bytes[j + (9 * 8)];
                ordered_bytes[(j * 40) + 10] = unordered_bytes[j + (10 * 8)];
                ordered_bytes[(j * 40) + 11] = unordered_bytes[j + (11 * 8)];
                ordered_bytes[(j * 40) + 12] = unordered_bytes[j + (12 * 8)];
                ordered_bytes[(j * 40) + 13] = unordered_bytes[j + (13 * 8)];
                ordered_bytes[(j * 40) + 14] = unordered_bytes[j + (14 * 8)];
                ordered_bytes[(j * 40) + 15] = unordered_bytes[j + (15 * 8)];
                ordered_bytes[(j * 40) + 16] = unordered_bytes[j + (16 * 8)];
                ordered_bytes[(j * 40) + 17] = unordered_bytes[j + (17 * 8)];
                ordered_bytes[(j * 40) + 18] = unordered_bytes[j + (18 * 8)];
                ordered_bytes[(j * 40) + 19] = unordered_bytes[j + (19 * 8)];
                ordered_bytes[(j * 40) + 20] = unordered_bytes[j + (20 * 8)];
                ordered_bytes[(j * 40) + 21] = unordered_bytes[j + (21 * 8)];
                ordered_bytes[(j * 40) + 22] = unordered_bytes[j + (22 * 8)];
                ordered_bytes[(j * 40) + 23] = unordered_bytes[j + (23 * 8)];
                ordered_bytes[(j * 40) + 24] = unordered_bytes[j + (24 * 8)];
                ordered_bytes[(j * 40) + 25] = unordered_bytes[j + (25 * 8)];
                ordered_bytes[(j * 40) + 26] = unordered_bytes[j + (26 * 8)];
                ordered_bytes[(j * 40) + 27] = unordered_bytes[j + (27 * 8)];
                ordered_bytes[(j * 40) + 28] = unordered_bytes[j + (28 * 8)];
                ordered_bytes[(j * 40) + 29] = unordered_bytes[j + (29 * 8)];
                ordered_bytes[(j * 40) + 30] = unordered_bytes[j + (30 * 8)];
                ordered_bytes[(j * 40) + 31] = unordered_bytes[j + (31 * 8)];
                ordered_bytes[(j * 40) + 32] = unordered_bytes[j + (32 * 8)];
                ordered_bytes[(j * 40) + 33] = unordered_bytes[j + (33 * 8)];
            }

            for (j = 0; j < 8; j++)
            {
                // then compute twins in order for flags that are set.
                while (ordered64[(j * 5) + 0] > 0)
                {
                    uint64_t pos = _trail_zcnt64(ordered64[(j * 5) + 0]);
                    uint64_t prime = lowlimit + (pos / 48) * 210 + sdata->rclass[pos % 48];

                    primes[GLOBAL_OFFSET + pcount++] = prime;
                    ordered64[(j * 5) + 0] = _reset_lsb64(ordered64[(j * 5) + 0]);
                }
            }
        }

    }
    else
    {
        // ordering the bits becomes inefficient with 48 lines because
        // they would need to be dispersed over too great a distance.
        // instead we compute the primes as before but push the results
        // into 64 different queues depending on the bit position.  Then
        // we pull from the queues in order while storing into the primes array.
        // This time the bottleneck is mostly in the queue-based sorting
        // and associated memory operations, so we don't bother with
        // switching between branch-free inner loops or not.
        uint64_t pqueues[64][48];
        uint8_t pcounts[64];
        int i,j;
        uint32_t current_line;

        memset(pcounts, 0, 64);
           
        lowlimit += byte_offset * 8 * 210;
        for (current_line = 0; current_line < nc; current_line++)
        {
            uint64_t *line64 = (uint64_t *)lines[current_line];
            uint64_t flags64 = line64[byte_offset/8];

            while (flags64 > 0)
            {
                uint64_t pos = _trail_zcnt64(flags64);
                uint64_t prime = lowlimit + pos * 210 + sdata->rclass[current_line];

                if ((prime >= olow) && (prime <= ohigh))
                {
                    pqueues[pos][pcounts[pos]] = prime;
                    pcounts[pos]++;
                }
                flags64 = _reset_lsb64(flags64);
            }
        }

        for (i = 0; i < 64; i++)
        {
            for (j = 0; j < pcounts[i] / 4; j++)
            {
                __m256i t = _mm256_loadu_si256((__m256i *)(&pqueues[i][j*4]));
                _mm256_storeu_si256((__m256i *)(&primes[GLOBAL_OFFSET + pcount]), t);
                pcount += 4;
            }
            for (j *= 4; j < pcounts[i]; j++)
            {
                primes[GLOBAL_OFFSET + pcount++] = pqueues[i][j];
            }
        }
    }

    return pcount;
}

#endif

