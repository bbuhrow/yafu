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

int _64_MOD_P[9] = {4,1,9,12,13,7,18,6,2};

uint64 _5_MASKS[5] = {
	0xef7bdef7bdef7bde,
	0xdef7bdef7bdef7bd,
	0xbdef7bdef7bdef7b,
	0x7bdef7bdef7bdef7,
	0xf7bdef7bdef7bdef};

uint64 _7_MASKS[7] = {
	0x7efdfbf7efdfbf7e,
	0xfdfbf7efdfbf7efd,
	0xfbf7efdfbf7efdfb,
	0xf7efdfbf7efdfbf7,
	0xefdfbf7efdfbf7ef,
	0xdfbf7efdfbf7efdf,
	0xbf7efdfbf7efdfbf};

uint64 _11_MASKS[11] = {
	0xff7feffdffbff7fe,
	0xfeffdffbff7feffd,
	0xfdffbff7feffdffb,
	0xfbff7feffdffbff7,
	0xf7feffdffbff7fef,
	0xeffdffbff7feffdf,
	0xdffbff7feffdffbf,
	0xbff7feffdffbff7f,
	0x7feffdffbff7feff,
	0xffdffbff7feffdff,
	0xffbff7feffdffbff};

uint64 _13_MASKS[13] = {
	0xffefff7ffbffdffe,
	0xffdffefff7ffbffd,
	0xffbffdffefff7ffb,
	0xff7ffbffdffefff7,
	0xfefff7ffbffdffef,
	0xfdffefff7ffbffdf,
	0xfbffdffefff7ffbf,
	0xf7ffbffdffefff7f,
	0xefff7ffbffdffeff,
	0xdffefff7ffbffdff,
	0xbffdffefff7ffbff,
	0x7ffbffdffefff7ff,
	0xfff7ffbffdffefff};

uint64 _17_MASKS[17] = {
	0xfff7fffbfffdfffe,
	0xffeffff7fffbfffd,
	0xffdfffeffff7fffb,
	0xffbfffdfffeffff7,
	0xff7fffbfffdfffef,
	0xfeffff7fffbfffdf,
	0xfdfffeffff7fffbf,
	0xfbfffdfffeffff7f,
	0xf7fffbfffdfffeff,
	0xeffff7fffbfffdff,
	0xdfffeffff7fffbff,
	0xbfffdfffeffff7ff,
	0x7fffbfffdfffefff,
	0xffff7fffbfffdfff,
	0xfffeffff7fffbfff,
	0xfffdfffeffff7fff,
	0xfffbfffdfffeffff};

uint64 _19_MASKS[19] = {
	0xfdffffbffff7fffe,
	0xfbffff7fffeffffd,
	0xf7fffeffffdffffb,
	0xeffffdffffbffff7,
	0xdffffbffff7fffef,
	0xbffff7fffeffffdf,
	0x7fffeffffdffffbf,
	0xffffdffffbffff7f,
	0xffffbffff7fffeff,
	0xffff7fffeffffdff,
	0xfffeffffdffffbff,
	0xfffdffffbffff7ff,
	0xfffbffff7fffefff,
	0xfff7fffeffffdfff,
	0xffeffffdffffbfff,
	0xffdffffbffff7fff,
	0xffbffff7fffeffff,
	0xff7fffeffffdffff,
	0xfeffffdffffbffff};

uint64 _23_MASKS[23] = {
	0xffffbfffff7ffffe,
	0xffff7ffffefffffd,
	0xfffefffffdfffffb,
	0xfffdfffffbfffff7,
	0xfffbfffff7ffffef,
	0xfff7ffffefffffdf,
	0xffefffffdfffffbf,
	0xffdfffffbfffff7f,
	0xffbfffff7ffffeff,
	0xff7ffffefffffdff,
	0xfefffffdfffffbff,
	0xfdfffffbfffff7ff,
	0xfbfffff7ffffefff,
	0xf7ffffefffffdfff,
	0xefffffdfffffbfff,
	0xdfffffbfffff7fff,
	0xbfffff7ffffeffff,
	0x7ffffefffffdffff,
	0xfffffdfffffbffff,
	0xfffffbfffff7ffff,
	0xfffff7ffffefffff,
	0xffffefffffdfffff,
	0xffffdfffffbfffff};

uint64 _29_MASKS[29] = {
	0xfbffffffdffffffe,
	0xf7ffffffbffffffd,
	0xefffffff7ffffffb,
	0xdffffffefffffff7,
	0xbffffffdffffffef,
	0x7ffffffbffffffdf,
	0xfffffff7ffffffbf,
	0xffffffefffffff7f,
	0xffffffdffffffeff,
	0xffffffbffffffdff,
	0xffffff7ffffffbff,
	0xfffffefffffff7ff,
	0xfffffdffffffefff,
	0xfffffbffffffdfff,
	0xfffff7ffffffbfff,
	0xffffefffffff7fff,
	0xffffdffffffeffff,
	0xffffbffffffdffff,
	0xffff7ffffffbffff,
	0xfffefffffff7ffff,
	0xfffdffffffefffff,
	0xfffbffffffdfffff,
	0xfff7ffffffbfffff,
	0xffefffffff7fffff,
	0xffdffffffeffffff,
	0xffbffffffdffffff,
	0xff7ffffffbffffff,
	0xfefffffff7ffffff,
	0xfdffffffefffffff};

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
				
				//then compute their next hit
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
			//finish with primes greater than (flagblocklimit >> 2)
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
