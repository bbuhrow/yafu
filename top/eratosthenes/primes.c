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

uint64 primes_from_lineflags(soe_staticdata_t *sdata, 
	uint32 start_count, uint64 *primes)
{
	//compute primes using all of the sieved lines we have stored
	uint32 current_line, pcount = start_count;
	uint64 prime;
	uint64 i;
	int pchar = 0;
	uint32 nc = sdata->numclasses;
	uint32 *rclass = sdata->rclass;
	uint64 lowlimit = sdata->lowlimit;
	uint64 prodN = sdata->prodN;
	uint8 **lines = sdata->lines;
	uint64 olow = sdata->orig_llimit;
	uint64 ohigh = sdata->orig_hlimit;
	
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

	for (i=0; i < sdata->numlinebytes; i+=8)
	{
		int b;	
		
		if ((i & 32767) == 0)
		{
			if (VFLAG > 0)
			{
				int k;
				for (k = 0; k<pchar; k++)
					printf("\b");
				pchar = printf("computing: %d%%",(int)((double)i / (double)(sdata->numlinebytes) * 100.0));
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
			cache_word[line_div8] |= ((uint64)lines[current_line][i] << (line_mod8 << 3));
			cache_word[8 + line_div8] |= ((uint64)lines[current_line][i+1] << (line_mod8 << 3));
			cache_word[16 + line_div8] |= ((uint64)lines[current_line][i+2] << (line_mod8 << 3));
			cache_word[24 + line_div8] |= ((uint64)lines[current_line][i+3] << (line_mod8 << 3));
			cache_word[32 + line_div8] |= ((uint64)lines[current_line][i+4] << (line_mod8 << 3));
			cache_word[40 + line_div8] |= ((uint64)lines[current_line][i+5] << (line_mod8 << 3));
			cache_word[48 + line_div8] |= ((uint64)lines[current_line][i+6] << (line_mod8 << 3));
			cache_word[56 + line_div8] |= ((uint64)lines[current_line][i+7] << (line_mod8 << 3));
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
				if (cache_word[((b >> 3) << 3) + (current_line >> 3)] & nmasks64[b & 7][current_line & 7])
				{
					prime = prodN * ((i << 3) + b) + rclass[current_line] + lowlimit;

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
	}

	if (VFLAG > 0)
	{
		//don't print status if computing primes, because lots of routines within
		//yafu do this and they don't want this side effect
		for (i = 0; i<pchar; i++)
			printf("\b");
	}

	return pcount;
}

