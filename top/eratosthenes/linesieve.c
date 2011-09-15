/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

       				   --bbuhrow@gmail.com 7/28/10
----------------------------------------------------------------------*/

#include "soe.h"

//masks for sieving multiple locations at once: small primes
//that hit a 64 bit interval more than once.  
int _64_MOD_P[9] = {4,1,9,12,13,7,18,6,2};

uint64 _5_MASKS[5] = {
	0xef7bdef7bdef7bdeULL,
	0xdef7bdef7bdef7bdULL,
	0xbdef7bdef7bdef7bULL,
	0x7bdef7bdef7bdef7ULL,
	0xf7bdef7bdef7bdefULL};

uint64 _7_MASKS[7] = {
	0x7efdfbf7efdfbf7eULL,
	0xfdfbf7efdfbf7efdULL,
	0xfbf7efdfbf7efdfbULL,
	0xf7efdfbf7efdfbf7ULL,
	0xefdfbf7efdfbf7efULL,
	0xdfbf7efdfbf7efdfULL,
	0xbf7efdfbf7efdfbfULL};

uint64 _11_MASKS[11] = {
	0xff7feffdffbff7feULL,
	0xfeffdffbff7feffdULL,
	0xfdffbff7feffdffbULL,
	0xfbff7feffdffbff7ULL,
	0xf7feffdffbff7fefULL,
	0xeffdffbff7feffdfULL,
	0xdffbff7feffdffbfULL,
	0xbff7feffdffbff7fULL,
	0x7feffdffbff7feffULL,
	0xffdffbff7feffdffULL,
	0xffbff7feffdffbffULL};

uint64 _13_MASKS[13] = {
	0xffefff7ffbffdffeULL,
	0xffdffefff7ffbffdULL,
	0xffbffdffefff7ffbULL,
	0xff7ffbffdffefff7ULL,
	0xfefff7ffbffdffefULL,
	0xfdffefff7ffbffdfULL,
	0xfbffdffefff7ffbfULL,
	0xf7ffbffdffefff7fULL,
	0xefff7ffbffdffeffULL,
	0xdffefff7ffbffdffULL,
	0xbffdffefff7ffbffULL,
	0x7ffbffdffefff7ffULL,
	0xfff7ffbffdffefffULL};

uint64 _17_MASKS[17] = {
	0xfff7fffbfffdfffeULL,
	0xffeffff7fffbfffdULL,
	0xffdfffeffff7fffbULL,
	0xffbfffdfffeffff7ULL,
	0xff7fffbfffdfffefULL,
	0xfeffff7fffbfffdfULL,
	0xfdfffeffff7fffbfULL,
	0xfbfffdfffeffff7fULL,
	0xf7fffbfffdfffeffULL,
	0xeffff7fffbfffdffULL,
	0xdfffeffff7fffbffULL,
	0xbfffdfffeffff7ffULL,
	0x7fffbfffdfffefffULL,
	0xffff7fffbfffdfffULL,
	0xfffeffff7fffbfffULL,
	0xfffdfffeffff7fffULL,
	0xfffbfffdfffeffffULL};

uint64 _19_MASKS[19] = {
	0xfdffffbffff7fffeULL,
	0xfbffff7fffeffffdULL,
	0xf7fffeffffdffffbULL,
	0xeffffdffffbffff7ULL,
	0xdffffbffff7fffefULL,
	0xbffff7fffeffffdfULL,
	0x7fffeffffdffffbfULL,
	0xffffdffffbffff7fULL,
	0xffffbffff7fffeffULL,
	0xffff7fffeffffdffULL,
	0xfffeffffdffffbffULL,
	0xfffdffffbffff7ffULL,
	0xfffbffff7fffefffULL,
	0xfff7fffeffffdfffULL,
	0xffeffffdffffbfffULL,
	0xffdffffbffff7fffULL,
	0xffbffff7fffeffffULL,
	0xff7fffeffffdffffULL,
	0xfeffffdffffbffffULL};

uint64 _23_MASKS[23] = {
	0xffffbfffff7ffffeULL,
	0xffff7ffffefffffdULL,
	0xfffefffffdfffffbULL,
	0xfffdfffffbfffff7ULL,
	0xfffbfffff7ffffefULL,
	0xfff7ffffefffffdfULL,
	0xffefffffdfffffbfULL,
	0xffdfffffbfffff7fULL,
	0xffbfffff7ffffeffULL,
	0xff7ffffefffffdffULL,
	0xfefffffdfffffbffULL,
	0xfdfffffbfffff7ffULL,
	0xfbfffff7ffffefffULL,
	0xf7ffffefffffdfffULL,
	0xefffffdfffffbfffULL,
	0xdfffffbfffff7fffULL,
	0xbfffff7ffffeffffULL,
	0x7ffffefffffdffffULL,
	0xfffffdfffffbffffULL,
	0xfffffbfffff7ffffULL,
	0xfffff7ffffefffffULL,
	0xffffefffffdfffffULL,
	0xffffdfffffbfffffULL};

uint64 _29_MASKS[29] = {
	0xfbffffffdffffffeULL,
	0xf7ffffffbffffffdULL,
	0xefffffff7ffffffbULL,
	0xdffffffefffffff7ULL,
	0xbffffffdffffffefULL,
	0x7ffffffbffffffdfULL,
	0xfffffff7ffffffbfULL,
	0xffffffefffffff7fULL,
	0xffffffdffffffeffULL,
	0xffffffbffffffdffULL,
	0xffffff7ffffffbffULL,
	0xfffffefffffff7ffULL,
	0xfffffdffffffefffULL,
	0xfffffbffffffdfffULL,
	0xfffff7ffffffbfffULL,
	0xffffefffffff7fffULL,
	0xffffdffffffeffffULL,
	0xffffbffffffdffffULL,
	0xffff7ffffffbffffULL,
	0xfffefffffff7ffffULL,
	0xfffdffffffefffffULL,
	0xfffbffffffdfffffULL,
	0xfff7ffffffbfffffULL,
	0xffefffffff7fffffULL,
	0xffdffffffeffffffULL,
	0xffbffffffdffffffULL,
	0xff7ffffffbffffffULL,
	0xfefffffff7ffffffULL,
	0xfdffffffefffffffULL};

void sieve_line(thread_soedata_t *thread_data)
{
	//extract stuff from the thread data structure
	soe_dynamicdata_t *ddata = &thread_data->ddata;
	soe_staticdata_t *sdata = &thread_data->sdata;
	uint8 *line = thread_data->ddata.line;
	uint32 current_line = thread_data->current_line;
	
	//stuff for bucket sieving
	soe_bucket_t *bptr;
	soe_bucket_t **buckets;
	uint32 *nptr;
	uint32 linesize = FLAGSIZE * sdata->blocks, bnum;
	//uint32 *bptr;

	uint8 *flagblock;
	uint64 startprime = sdata->startprime;
	uint64 i,j,k;
	uint32 prime;
	uint32 maxP;

	ddata->lblk_b = sdata->lowlimit + sdata->rclass[current_line];
	ddata->ublk_b = sdata->blk_r + ddata->lblk_b - sdata->prodN;
	ddata->blk_b_sqrt = (uint32)(sqrt(ddata->ublk_b + sdata->prodN)) + 1;

	//for the current line, find the offsets past the low limit
	get_offsets(thread_data);

	flagblock = line;
	for (i=0;i<sdata->blocks;i++)
	{
		uint64 *flagblock64;
		int mask_step, mask_step2;
		int mask_num, mask_num2;

		//sieve the block with each effective prime
		//set all flags for this block, which also puts it into cache for the sieving
		//to follow
		memset(flagblock,255,BLOCKSIZE);

		flagblock64 = (uint64 *)flagblock;		
		
		//do the smallest primes in predetermined 64 bit batches
		if (startprime == 2)
		{
			for (k=0, mask_step = _64_MOD_P[0], mask_num = ddata->offsets[2],
				mask_step2 = _64_MOD_P[1], mask_num2 = ddata->offsets[3]; 
				k<FLAGSIZE >> 6; k++)
			{
				flagblock64[k] &= (_5_MASKS[mask_num] & _7_MASKS[mask_num2]);
				mask_num -= mask_step;
				if (mask_num < 0) mask_num = 5 + mask_num;
				mask_num2 -= mask_step2;
				if (mask_num2 < 0) mask_num2 = 7 + mask_num2;
			}
			ddata->offsets[2]= (uint32)mask_num;
			ddata->offsets[3]= (uint32)mask_num2;

			for (k=0, mask_step = _64_MOD_P[2], mask_num = ddata->offsets[4],
				mask_step2 = _64_MOD_P[3], mask_num2 = ddata->offsets[5]; 
				k<FLAGSIZE >> 6; k++)
			{
				flagblock64[k] &= (_11_MASKS[mask_num] & _13_MASKS[mask_num2]);
				mask_num -= mask_step;
				if (mask_num < 0) mask_num = 11 + mask_num;
				mask_num2 -= mask_step2;
				if (mask_num2 < 0) mask_num2 = 13 + mask_num2;
			}
			ddata->offsets[4]= (uint32)mask_num;
			ddata->offsets[5]= (uint32)mask_num2;
		}
		else if (startprime == 3)
		{
			for (k=0, mask_step = _64_MOD_P[1], mask_num = ddata->offsets[3]; 
				k<FLAGSIZE >> 6; k++)
			{
				flagblock64[k] &= _7_MASKS[mask_num];
				mask_num -= mask_step;
				if (mask_num < 0) mask_num = 7 + mask_num;
			}
			ddata->offsets[3]= (uint32)mask_num;

			for (k=0, mask_step = _64_MOD_P[2], mask_num = ddata->offsets[4],
				mask_step2 = _64_MOD_P[3], mask_num2 = ddata->offsets[5]; 
				k<FLAGSIZE >> 6; k++)
			{
				flagblock64[k] &= (_11_MASKS[mask_num] & _13_MASKS[mask_num2]);
				mask_num -= mask_step;
				if (mask_num < 0) mask_num = 11 + mask_num;
				mask_num2 -= mask_step2;
				if (mask_num2 < 0) mask_num2 = 13 + mask_num2;
			}
			ddata->offsets[4]= (uint32)mask_num;
			ddata->offsets[5]= (uint32)mask_num2;
		}
		else if (startprime == 4)
		{

			for (k=0, mask_step = _64_MOD_P[2], mask_num = ddata->offsets[4],
				mask_step2 = _64_MOD_P[3], mask_num2 = ddata->offsets[5]; 
				k<FLAGSIZE >> 6; k++)
			{
				flagblock64[k] &= (_11_MASKS[mask_num] & _13_MASKS[mask_num2]);
				mask_num -= mask_step;
				if (mask_num < 0) mask_num = 11 + mask_num;
				mask_num2 -= mask_step2;
				if (mask_num2 < 0) mask_num2 = 13 + mask_num2;
			}
			ddata->offsets[4]= (uint32)mask_num;
			ddata->offsets[5]= (uint32)mask_num2;
		}
		else if (startprime == 5)
		{
			for (k=0, mask_step = _64_MOD_P[3], mask_num = ddata->offsets[5];
				k<FLAGSIZE >> 6; k++)
			{
				flagblock64[k] &= _13_MASKS[mask_num];
				mask_num -= mask_step;
				if (mask_num < 0) mask_num = 13 + mask_num;
			}
			ddata->offsets[5]= (uint32)mask_num;
		}	

		for (k=0, mask_step = _64_MOD_P[4], mask_num = ddata->offsets[6],
			mask_step2 = _64_MOD_P[5], mask_num2 = ddata->offsets[7]; 
			k<FLAGSIZE >> 6; k++)
		{
			flagblock64[k] &= (_17_MASKS[mask_num] & _19_MASKS[mask_num2]);
			mask_num -= mask_step;
			if (mask_num < 0) mask_num = 17 + mask_num;
			mask_num2 -= mask_step2;
			if (mask_num2 < 0) mask_num2 = 19 + mask_num2;
		}
		ddata->offsets[6]= (uint32)mask_num;
		ddata->offsets[7]= (uint32)mask_num2;

		for (k=0, mask_step = _64_MOD_P[6], mask_num = ddata->offsets[8],
			mask_step2 = _64_MOD_P[7], mask_num2 = ddata->offsets[9]; 
			k<FLAGSIZE >> 6; k++)
		{
			flagblock64[k] &= (_23_MASKS[mask_num] & _29_MASKS[mask_num2]);
			mask_num -= mask_step;
			if (mask_num < 0) mask_num = 23 + mask_num;
			mask_num2 -= mask_step2;
			if (mask_num2 < 0) mask_num2 = 29 + mask_num2;
		}
		ddata->offsets[8]= (uint32)mask_num;
		ddata->offsets[9]= (uint32)mask_num2;	
		
		//one is not a prime
		if ((sdata->rclass[current_line] == 1) &&
			(sdata->lowlimit <= 1) && (i == 0))
			flagblock[0] &= 0xfe;

		//unroll the loop: all primes less than this max hit the interval at least 16 times
		maxP = FLAGSIZE >> 4;

		for (j=10;j<ddata->pbounds[i];j++)
		{
			uint32 tmpP;
			uint64 stop;
			uint64 p1,p2,p3;

			prime = sdata->sieve_p[j];
			if (prime > maxP)
				break;

			tmpP = prime << 4;
			stop = FLAGSIZE - tmpP + prime;
			k=ddata->offsets[j];
			p1 = prime;
			p2 = p1 + prime;
			p3 = p2 + prime;
			while (k < stop)
			{
				flagblock[k>>3] &= masks[k&7];
				flagblock[(k+p1)>>3] &= masks[(k+p1)&7];
				flagblock[(k+p2)>>3] &= masks[(k+p2)&7];
				flagblock[(k+p3)>>3] &= masks[(k+p3)&7];
				k += (prime << 2);
				flagblock[k>>3] &= masks[k&7];
				flagblock[(k+p1)>>3] &= masks[(k+p1)&7];
				flagblock[(k+p2)>>3] &= masks[(k+p2)&7];
				flagblock[(k+p3)>>3] &= masks[(k+p3)&7];
				k += (prime << 2);
				flagblock[k>>3] &= masks[k&7];
				flagblock[(k+p1)>>3] &= masks[(k+p1)&7];
				flagblock[(k+p2)>>3] &= masks[(k+p2)&7];
				flagblock[(k+p3)>>3] &= masks[(k+p3)&7];
				k += (prime << 2);
				flagblock[k>>3] &= masks[k&7];
				flagblock[(k+p1)>>3] &= masks[(k+p1)&7];
				flagblock[(k+p2)>>3] &= masks[(k+p2)&7];
				flagblock[(k+p3)>>3] &= masks[(k+p3)&7];
				k += (prime << 2);
			}

			for (;k<FLAGSIZE;k+=prime)
				flagblock[k>>3] &= masks[k&7];


			
			//if ((j >= 2) && (j <= 10))
			//{
			//	printf("actual = %lx; next offset = %u\n",flagblock64[0], k - FLAGSIZE);
			//}
			ddata->offsets[j]= (uint32)(k - FLAGSIZE);
			
		}

		//unroll the loop: all primes less than this max hit the interval at least 8 times
		maxP = FLAGSIZE >> 3;

		for (;j<ddata->pbounds[i];j++)
		{
			uint32 tmpP;
			uint64 stop;
			uint64 p2, p4;

			prime = sdata->sieve_p[j];
			if (prime > maxP)
				break;

			tmpP = prime << 3;
			stop = FLAGSIZE - tmpP + prime;
			k=ddata->offsets[j];
			p2 = prime<<1;
			p4 = prime<<2;

			while (k < stop)
			{
				flagblock[k>>3] &= masks[k&7];								//0 * prime
				flagblock[(k+prime)>>3] &= masks[(k+prime)&7];				//1 * prime
				flagblock[(k+p2)>>3] &= masks[(k+p2)&7];					//2 * prime
				flagblock[(k+prime+p2)>>3] &= masks[(k+prime+p2)&7];		//3 * prime
				flagblock[(k+p4)>>3] &= masks[(k+p4)&7];					//4 * prime
				flagblock[(k+prime+p4)>>3] &= masks[(k+prime+p4)&7];		//5 * prime
				flagblock[(k+p2+p4)>>3] &= masks[(k+p2+p4)&7];				//6 * prime
				flagblock[(k+prime+p2+p4)>>3] &= masks[(k+prime+p2+p4)&7];	//7 * prime
				k += (prime << 3);											//advance
			}

			for (;k<FLAGSIZE;k+=prime)								//finish
				flagblock[k>>3] &= masks[k&7];

			ddata->offsets[j]= (uint32)(k - FLAGSIZE);
		}

		//unroll the loop: all primes less than this max hit the interval at least 4 times
		maxP = FLAGSIZE >> 2;

		for (;j<ddata->pbounds[i];j++)
		{
			uint32 tmpP;
			uint64 stop;
			uint64 p2;

			prime = sdata->sieve_p[j];
			if (prime > maxP)
				break;

			tmpP = prime << 2;
			stop = FLAGSIZE - tmpP + prime;
			k=ddata->offsets[j];
			p2 = prime<<1;
			while (k < stop)
			{
				flagblock[k>>3] &= masks[k&7];								//0 * prime
				flagblock[(k+prime)>>3] &= masks[(k+prime)&7];				//1 * prime
				flagblock[(k+p2)>>3] &= masks[(k+p2)&7];					//2 * prime
				flagblock[(k+prime+p2)>>3] &= masks[(k+prime+p2)&7];		//3 * prime
				k += (prime << 2);											//advance
			}

			for (;k<FLAGSIZE;k+=prime)								//finish
				flagblock[k>>3] &= masks[k&7];

			ddata->offsets[j]= (uint32)(k - FLAGSIZE);
		}
		
		if (ddata->bucket_depth > BUCKET_BUFFER)
		{
			for (;j<ddata->pbounds[i];j++)
			{
				prime = sdata->sieve_p[j];
				if (prime > BUCKETSTARTP)
					break;
				for (k=ddata->offsets[j];k<FLAGSIZE;k+=prime)
					flagblock[k>>3] &= masks[k&7];

				ddata->offsets[j]= (uint32)(k - FLAGSIZE);
			}

			//finally, fill any primes in this block's bucket
			bptr = ddata->sieve_buckets[i];	
			buckets = ddata->sieve_buckets;
			nptr = ddata->bucket_hits;

			//printf("unloading %d hits in block %d of line %d\n",nptr[i],i,thread_data->current_line);
			for (j=0; j < (nptr[i] & (uint32)(~7)); j+=8)
			{				
				//unload 8 hits
				flagblock[(bptr[j + 0].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 0].root & FLAGSIZEm1) & 7];
				flagblock[(bptr[j + 1].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 1].root & FLAGSIZEm1) & 7];
				flagblock[(bptr[j + 2].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 2].root & FLAGSIZEm1) & 7];
				flagblock[(bptr[j + 3].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 3].root & FLAGSIZEm1) & 7];
				flagblock[(bptr[j + 4].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 4].root & FLAGSIZEm1) & 7];
				flagblock[(bptr[j + 5].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 5].root & FLAGSIZEm1) & 7];
				flagblock[(bptr[j + 6].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 6].root & FLAGSIZEm1) & 7];
				flagblock[(bptr[j + 7].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 7].root & FLAGSIZEm1) & 7];
				
				//then compute their next hit and update the roots while they are
				//still fresh in the cache
				bptr[j + 0].root += bptr[j + 0].prime;	
				bptr[j + 1].root += bptr[j + 1].prime;		
				bptr[j + 2].root += bptr[j + 2].prime;	

				if (bptr[j + 0].root < linesize)			
				{	
					bnum = (bptr[j + 0].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j + 0].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j + 0].prime;
					nptr[bnum]++;						
				}	

				if (bptr[j + 1].root < linesize)			
				{	
					bnum = (bptr[j + 1].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j + 1].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j + 1].prime;
					nptr[bnum]++;						
				}	

				if (bptr[j + 2].root < linesize)			
				{	
					bnum = (bptr[j + 2].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j + 2].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j + 2].prime;
					nptr[bnum]++;						
				}	

				bptr[j + 3].root += bptr[j + 3].prime;		
				bptr[j + 4].root += bptr[j + 4].prime;	
				bptr[j + 5].root += bptr[j + 5].prime;		
				if (bptr[j + 3].root < linesize)			
				{	
					bnum = (bptr[j + 3].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j + 3].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j + 3].prime;
					nptr[bnum]++;						
				}	

				if (bptr[j + 4].root < linesize)			
				{	
					bnum = (bptr[j + 4].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j + 4].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j + 4].prime;
					nptr[bnum]++;						
				}	

				if (bptr[j + 5].root < linesize)			
				{	
					bnum = (bptr[j + 5].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j + 5].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j + 5].prime;
					nptr[bnum]++;						
				}	

				bptr[j + 6].root += bptr[j + 6].prime;		
				bptr[j + 7].root += bptr[j + 7].prime;		

				if (bptr[j + 6].root < linesize)			
				{	
					bnum = (bptr[j + 6].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j + 6].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j + 6].prime;
					nptr[bnum]++;						
				}	

				if (bptr[j + 7].root < linesize)			
				{	
					bnum = (bptr[j + 7].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j + 7].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j + 7].prime;
					nptr[bnum]++;						
				}	

			}

			//finish up those that didn't fit into a group of 8 hits
			for (;j < nptr[i]; j++)
			{
				flagblock[(bptr[j].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j].root & FLAGSIZEm1) & 7];

				bptr[j].root += bptr[j].prime;		
				if (bptr[j].root < linesize)			
				{	
					bnum = (bptr[j].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j].prime;
					nptr[bnum]++;						
				}
				
			}

			// repeat the dumping of bucket primes, this time with very large primes
			// that only hit the interval once.  thus, we don't need to update the root
			// with the next hit, and we can do more at once because each bucket hit is smaller
			if (ddata->largep_offset > 0)
			{
				uint32 *large_bptr = ddata->large_sieve_buckets[i];	
				uint32 *large_nptr = ddata->large_bucket_hits;

				for (j=0; j < (large_nptr[i] & (uint32)(~15)); j+=16)
				{				
					//unload 8 hits
					flagblock[(large_bptr[j + 0] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 0] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 1] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 1] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 2] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 2] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 3] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 3] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 4] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 4] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 5] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 5] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 6] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 6] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 7] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 7] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 8] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 8] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 9] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 9] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 10] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 10] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 11] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 11] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 12] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 12] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 13] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 13] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 14] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 14] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 15] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 15] & FLAGSIZEm1) & 7];
				}

				for (;j < large_nptr[i]; j++)
					flagblock[(large_bptr[j] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j] & FLAGSIZEm1) & 7];

			}

		}
		else
		{
			//didn't need to use a bucket sieve
			//finish with primes greater than (flagblocklimit >> 2) that we
			//didn't unroll.
			for (;j<ddata->pbounds[i];j++)
			{
				prime = sdata->sieve_p[j];
				for (k=ddata->offsets[j];k<FLAGSIZE;k+=prime)
					flagblock[k>>3] &= masks[k&7];

				ddata->offsets[j]= (uint32)(k - FLAGSIZE);
			}
		}

		flagblock += BLOCKSIZE;
	}

	return;
}
