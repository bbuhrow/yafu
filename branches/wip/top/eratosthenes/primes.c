/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

       				   --bbuhrow@gmail.com 7/1/10
----------------------------------------------------------------------*/

#include "soe.h"
#include "threadpool.h"

void compute_primes_dispatch(void *vptr)
{
    tpool_t *tdata = (tpool_t *)vptr;
    soe_userdata_t *t = (soe_userdata_t *)tdata->user_data;
    soe_staticdata_t *sdata = t->sdata;

    if (sdata->sync_count < THREADS)
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

    if (THREADS > 1)
        t->linecount = 0;

    for (i = t->startid; i < t->stopid; i += 8)
    {
        t->linecount = compute_8_bytes(sdata, t->linecount, t->ddata.primes, i);
    }

    return;
}





uint64 primes_from_lineflags(soe_staticdata_t *sdata, thread_soedata_t *thread_data,
	uint32 start_count, uint64 *primes)
{
	//compute primes using all of the sieved lines we have stored
	uint32 pcount = start_count;	
	uint64 i;
	int j;
    int k;
	uint32 range, lastid;

    // threading structures
    tpool_t *tpool_data;
    soe_userdata_t udata;

	// each thread needs to work on a number of bytes that is divisible by 8
	range = sdata->numlinebytes / THREADS;
	range -= (range % 8);
	lastid = 0;

    // divvy up the line bytes
    for (i = 0; i < THREADS; i++)
    {
        thread_soedata_t *t = thread_data + i;

        t->sdata = *sdata;
        t->startid = lastid;
        t->stopid = t->startid + range;
        lastid = t->stopid;

        if (VFLAG > 2)
        {
            printf("thread %d finding primes from byte offset %u to %u\n",
                (int)i, t->startid, t->stopid);
        }
    }

    // allocate a temporary array for each thread's primes
    if (THREADS > 1)
    {
        uint64 memchunk;

        if (sdata->sieve_range)
        {
            // then just split the overall range into equal parts
            memchunk = (sdata->orig_hlimit - sdata->orig_llimit) / THREADS + THREADS;

            for (i = 0; i < THREADS; i++)
            {
                thread_soedata_t *t = thread_data + i;

                t->ddata.primes = (uint64 *)malloc(memchunk * sizeof(uint64));
            }
        }
        else
        {
            // then estimate the number of primes we'll find in each chunk.
            // it's important to do this chunk by chunk, because in some cases
            // the number of primes changes rapidly as a function of offset
            // from the start of the range (i.e., when start is 0)
            uint64 hi_est, lo_est;
            uint64 tmplo = sdata->orig_llimit;
            uint64 tmphi;
            uint64 chunk = 8 * sdata->numlinebytes / THREADS;

            chunk *= sdata->prodN;
            tmphi = tmplo + chunk;
            for (i = 0; i < THREADS; i++)
            {
                thread_soedata_t *t = thread_data + i;

                hi_est = (uint64)(tmphi / log((double)tmphi));
                if (tmplo > 1)
                    lo_est = (uint64)(tmplo / log((double)tmplo));
                else
                    lo_est = 0;

                memchunk = (uint64)((double)(hi_est - lo_est) * 1.25);

                if (VFLAG > 2)
                {
                    printf("allocating temporary space for %" PRIu64 " primes between %" PRIu64 " and %" PRIu64 "\n",
                        memchunk, tmplo, tmphi);
                }

                t->ddata.primes = (uint64 *)malloc(memchunk * sizeof(uint64));

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
    tpool_data = tpool_setup(THREADS, NULL, NULL, NULL,
        &compute_primes_dispatch, &udata);

    if (THREADS == 1)
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
	if (THREADS > 1)
	{
		pcount = start_count;
		for (j = 0; j < THREADS; j++)
		{
			thread_soedata_t *t = thread_data + j;

			if (t->linecount == 0)
			{
				free(t->ddata.primes);
				continue;
			}

			if (VFLAG > 2)
				printf("adding %" PRIu64 " primes found in thread %d\n", t->linecount, j);

			memcpy(primes + GLOBAL_OFFSET + pcount, t->ddata.primes, t->linecount * sizeof(uint64));

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
		if (VFLAG > 2)
			printf("adding primes from byte offset %u to %u\n", 
				lastid, (uint32)sdata->numlinebytes);

		for (i = lastid; i < sdata->numlinebytes; i+=8)
		{
			pcount = compute_8_bytes(sdata, pcount, primes, i);		
		}
	}

	return pcount;
}

uint32 compute_8_bytes(soe_staticdata_t *sdata, 
	uint32 pcount, uint64 *primes, uint64 byte_offset)
{
	int b;	
	uint32 current_line;
	//8 bytes from each of up to 48 sieve lines are packed into these words
	uint64 cache_word[64];
	//for testing one of 8 bits in a byte in one of 8 lines.
	//bit num picks the row, lines num picks the col.
	uint64 nmasks64[8][8] = {
		{1ULL,256ULL,65536ULL,16777216ULL,4294967296ULL,1099511627776ULL,281474976710656ULL,72057594037927936ULL},
		{2ULL,512ULL,131072ULL,33554432ULL,8589934592ULL,2199023255552ULL,562949953421312ULL,144115188075855872ULL},
		{4ULL,1024ULL,262144ULL,67108864ULL,17179869184ULL,4398046511104ULL,1125899906842624ULL,288230376151711744ULL},
		{8ULL,2048ULL,524288ULL,134217728ULL,34359738368ULL,8796093022208ULL,2251799813685248ULL,576460752303423488ULL},
		{16ULL,4096ULL,1048576ULL,268435456ULL,68719476736ULL,17592186044416ULL,4503599627370496ULL,1152921504606846976ULL},
		{32ULL,8192ULL,2097152ULL,536870912ULL,137438953472ULL,35184372088832ULL,9007199254740992ULL,2305843009213693952ULL},
		{64ULL,16384ULL,4194304ULL,1073741824ULL,274877906944ULL,70368744177664ULL,18014398509481984ULL,4611686018427387904ULL},
		{128ULL,32768ULL,8388608ULL,2147483648ULL,549755813888ULL,140737488355328ULL,36028797018963968ULL,9223372036854775808ULL}};
	
	uint64 prime;
	uint32 nc = sdata->numclasses;
	uint32 *rclass = sdata->rclass;
	uint64 lowlimit = sdata->lowlimit;
	uint64 prodN = sdata->prodN;
	uint8 **lines = sdata->lines;
	uint64 olow = sdata->orig_llimit;
	uint64 ohigh = sdata->orig_hlimit;
		
	if ((byte_offset & 32767) == 0)
	{
		if (VFLAG > 1)
		{
			int k;
			printf("computing: %d%%\r",(int)
				((double)byte_offset / (double)(sdata->numlinebytes) * 100.0));
			fflush(stdout);
		}
	}

	//get 8 bytes from each residue class and pack into a series of 64 bit words.
	//then we can raster across those 64 bits much more efficiently
	memset(cache_word, 0, 64 * sizeof(uint64));
	for (current_line = 0; current_line < nc; current_line++)
	{
		//put 8 bytes from the current line into each of 8 different 64 bit words.
		//shift the byte left according to the current line mod 8 so that
		//each 64 bit word will eventually hold bytes from up to 8 different lines.
		//the bytes from the current line are spaced 8 words apart, so that there is
		//room to store up to 64 lines (capacity enough for the 48 line case mod 210).
		uint32 line_div8 = current_line >> 3;
		uint32 line_mod8 = current_line & 7;
		cache_word[line_div8] |= ((uint64)lines[current_line][byte_offset] << (line_mod8 << 3));
		cache_word[8 + line_div8] |= ((uint64)lines[current_line][byte_offset+1] << (line_mod8 << 3));
		cache_word[16 + line_div8] |= ((uint64)lines[current_line][byte_offset+2] << (line_mod8 << 3));
		cache_word[24 + line_div8] |= ((uint64)lines[current_line][byte_offset+3] << (line_mod8 << 3));
		cache_word[32 + line_div8] |= ((uint64)lines[current_line][byte_offset+4] << (line_mod8 << 3));
		cache_word[40 + line_div8] |= ((uint64)lines[current_line][byte_offset+5] << (line_mod8 << 3));
		cache_word[48 + line_div8] |= ((uint64)lines[current_line][byte_offset+6] << (line_mod8 << 3));
		cache_word[56 + line_div8] |= ((uint64)lines[current_line][byte_offset+7] << (line_mod8 << 3));
	}

	//for each bit
	//for (b = 0; b < 8; b++)
	for (b = 0; b < 64; b++)
	{
		for (current_line = 0; current_line < nc; current_line++)
		{
			//compute the prime at this location if it is flagged and 
			//within our original boundaries.  
			//if (lines[current_line][i] & nmasks[b])
			//if (cache_word & nmasks64[b])
			//if (cache_word & ((uint64)nmasks[b] << (current_line << 3)))
			//if (cache_word[current_line >> 3] & nmasks64[b][current_line & 7])
			//select the appropriate word according to the bit and line.
			//then 'and' it with the appropriate mask according to the bit and line.
			//all these bit operations are cheaper than continually fetching new
			//bytes from many different lines (cache optimization)
			if (cache_word[((b >> 3) << 3) + (current_line >> 3)] & 
				nmasks64[b & 7][current_line & 7])
			{
				prime = prodN * ((byte_offset << 3) + b) + rclass[current_line] + lowlimit;

				if ((prime >= olow) && (prime <= ohigh))
				{
					if (NO_STORE)
						pcount++;
					else
						primes[GLOBAL_OFFSET + pcount++] = prime;
				}
			}
		}
	}

	return pcount;
}

