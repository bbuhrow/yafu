/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

Some parts of the code (and also this header), included in this 
distribution have been reused from other sources. In particular I 
have benefitted greatly from the work of Jason Papadopoulos's msieve @ 
www.boo.net/~jasonp, Scott Contini's mpqs implementation, and Tom St. 
Denis Tom's Fast Math library.  Many thanks to their kind donation of 
code to the public domain.
       				   --bbuhrow@gmail.com 11/24/09
----------------------------------------------------------------------*/

#include "yafu.h"
#include "qs.h"
#include "ytools.h"
#include "common.h"
#include "poly_macros_common.h"
#include "poly_macros_32k.h"


#define FILL_ONE_PRIME_LOOP_P_test(i)				\
	bnum = root1 >> 15;					\
	while (bnum < numblocks)					\
	{											\
		bptr = sliceptr_p +						\
			(bnum << BUCKET_BITS) +				\
			numptr_p[bnum];					\
		*bptr = (((i) - bound_val) << 16 | (root1 & 32767)); \
		numptr_p[bnum]++;					\
		root1 += prime;							\
		bnum = root1 >> 15;				\
	}											\

void testfirstRoots_32k(static_conf_t *sconf, dynamic_conf_t *dconf)
{
	// the roots are computed using a and b as follows:
	// (+/-t - b)(a)^-1 mod p
	// where the t values are the roots to t^2 = N mod p, found by shanks_tonelli
	// when constructing the factor base.
	// assume b > t
       
	// compute the roots as if we were actually going to use this, but don't save
	// anything.  We are just trying to determine the size needed for each large 
	// prime bucket by sieving over just the first bucket

	uint32 i,logp;
	int root1, root2, prime, amodp, bmodp, inv, bnum,numblocks;
	int lpnum,last_bound;

	// unpack stuff from the job data
	siqs_poly *poly = dconf->curr_poly;
	fb_list *fb = sconf->factor_base;
	lp_bucket *lp_bucket_p = dconf->buckets;
	uint32 *modsqrt = sconf->modsqrt_array;

	numblocks = sconf->num_blocks;

	lpnum = 0;
	dconf->buckets->alloc_slices = 1;

	// extreme estimate for number of slices
	i = (sconf->factor_base->B - sconf->factor_base->med_B) / 512;

	last_bound = fb->med_B;
	for (i=fb->med_B;i<fb->B;i++)
	{
		prime = fb->list->prime[i];
		root1 = modsqrt[i]; 
		root2 = prime - root1; 
		logp = fb->list->logprime[i];

		amodp = (int)mpz_tdiv_ui(poly->mpz_poly_a,prime);
		bmodp = (int)mpz_tdiv_ui(poly->mpz_poly_b,prime);

		//find a^-1 mod p = inv(a mod p) mod p
		inv = modinv_1(amodp,prime);

		root1 = (int)root1 - bmodp;
		if (root1 < 0) root1 += prime;

		root2 = (int)root2 - bmodp;
		if (root2 < 0) root2 += prime;
	
		root1 = (uint32)((uint64)inv * (uint64)root1 % (uint64)prime);
		root2 = (uint32)((uint64)inv * (uint64)root2 % (uint64)prime);

		// just need to do this once, because the next step of prime will be 
		// into a different bucket
		bnum = root1 >> 15;
		if (bnum == 0)
			lpnum++;

		// repeat for the other root
		bnum = root2 >> 15;
		if (bnum == 0)
			lpnum++;

		if ((uint32)lpnum > (double)BUCKET_ALLOC * 0.75)
		{
			// we want to allocate more slices than we will probably need
			// assume alloc/2 is a safe amount of slack
			lp_bucket_p->alloc_slices++;
			lpnum = 0;
		}

		if (i - last_bound == 65536)
		{
			// when prime are really big, we may cross this boundary
			// before the buckets fill up
			lp_bucket_p->alloc_slices++;
			lpnum = 0;
			last_bound = i;
		}
	}

	// extra cushion - may increase the memory usage a bit, but in very
	// rare circumstances not enough slices allocated causes crashes.
	lp_bucket_p->alloc_slices++;

#ifdef DO_VLP_OPT
    // one more slice because we force a new one at the VLP boundary.
    lp_bucket_p->alloc_slices++;
#endif

	return;
}

void firstRoots_32k(static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//the roots are computed using a and b as follows:
	//(+/-t - b)(a)^-1 mod p
	//where the t values are the roots to t^2 = N mod p, found by shanks_tonelli
	//when constructing the factor base.
	//assume b > t

	//unpack stuff from the job data structures
	siqs_poly *poly = dconf->curr_poly;
	fb_list *fb = sconf->factor_base;
	uint32 start_prime = 2;
	int *rootupdates = dconf->rootupdates;
	update_t update_data = dconf->update_data;
	sieve_fb_compressed *fb_p = dconf->comp_sieve_p;
	sieve_fb_compressed *fb_n = dconf->comp_sieve_n;
	lp_bucket *lp_bucket_p = dconf->buckets;
	uint32 *modsqrt = sconf->modsqrt_array;

	//locals
	uint32 i, interval;
	uint8 logp;
    int root1, root2, nroot1, nroot2, prime, amodp, bmodp, inv, x, bnum, j, numblocks;
	int s = poly->s;
	int bound_index = 0, k;
	uint32 bound_val = fb->med_B;
	uint32 *bptr, *sliceptr_p, *sliceptr_n;
    uint64* bptr64, * sliceptr64_p, * sliceptr64_n;
	uint32 *numptr_p, *numptr_n;
	int check_bound = BUCKET_ALLOC/2 - 1, room;
    FILE *out;
    uint32 shift = 24;


	numblocks = sconf->num_blocks;
	interval = numblocks << 15;

	if (lp_bucket_p->list != NULL)
	{
		lp_bucket_p->fb_bounds[0] = fb->med_B;

		sliceptr_p = lp_bucket_p->list;
		sliceptr_n = lp_bucket_p->list + (numblocks << BUCKET_BITS);

		numptr_p = lp_bucket_p->num;
		numptr_n = lp_bucket_p->num + numblocks;
		//reset lp_buckets
		for (i=0;i< (2*numblocks*lp_bucket_p->alloc_slices) ;i++)
			numptr_p[i] = 0;

		lp_bucket_p->num_slices = 0;
	}
	else
	{
		sliceptr_p = NULL;
		sliceptr_n = NULL;
		numptr_p = NULL;
		numptr_n = NULL;
	}

	for (i = start_prime; i < sconf->sieve_small_fb_start; i++)
	{
		uint64 q64, tmp, t2;

		prime = fb->tinylist->prime[i];
		root1 = modsqrt[i]; 
		root2 = prime - root1; 

		amodp = (int)mpz_tdiv_ui(poly->mpz_poly_a,prime);
		bmodp = (int)mpz_tdiv_ui(poly->mpz_poly_b,prime);

		//find a^-1 mod p = inv(a mod p) mod p
		inv = modinv_1(amodp,prime);

		COMPUTE_FIRST_ROOTS
	
		// reuse integer inverse of prime that we've calculated for use
		// in trial division stage
		// inv * root1 % prime
		t2 = (uint64)inv * (uint64)root1;
		tmp = t2 + (uint64)fb->tinylist->correction[i];
		q64 = tmp * (uint64)fb->tinylist->small_inv[i];
		tmp = q64 >> 32; 
		root1 = t2 - tmp * prime;

		// inv * root2 % prime
		t2 = (uint64)inv * (uint64)root2;
		tmp = t2 + (uint64)fb->tinylist->correction[i];
		q64 = tmp * (uint64)fb->tinylist->small_inv[i];
		tmp = q64 >> 32; 
		root2 = t2 - tmp * prime;
	
		//we don't sieve these primes, so ordering doesn't matter
		update_data.firstroots1[i] = root1;
		update_data.firstroots2[i] = root2;

		fb_p->root1[i] = (uint16)root1;
		fb_p->root2[i] = (uint16)root2;
		fb_n->root1[i] = (uint16)(prime - root2);
		fb_n->root2[i] = (uint16)(prime - root1);
		//if we were sieving, this would double count the location on the 
		//positive side.  but since we're not, its easier to check for inclusion
		//on the progression if we reset the negative root to zero if it is == prime
		if (fb_n->root1[i] == prime)
			fb_n->root1[i] = 0;
		if (fb_n->root2[i] == prime)
			fb_n->root2[i] = 0;

		//for this factor base prime, compute the rootupdate value for all s
		//Bl values.  amodp holds a^-1 mod p
		//the rootupdate value is given by 2*Bj*amodp
		//Bl[j] now holds 2*Bl
		for (j=0;j<s;j++)
		{
			x = (int)mpz_tdiv_ui(dconf->Bl[j],prime);
			
			// x * inv % prime
			t2 = (uint64)inv * (uint64)x;
			tmp = t2 + (uint64)fb->tinylist->correction[i];
			q64 = tmp * (uint64)fb->tinylist->small_inv[i];
			tmp = q64 >> 32; 
			x = t2 - tmp * prime;

			rootupdates[(j)*fb->B+i] = x;
		}
	}

	for (i=sconf->sieve_small_fb_start;i<fb->fb_15bit_B;i++)
	{
		uint64 tmp, t2;

		prime = fb->list->prime[i];
		root1 = modsqrt[i]; 
		root2 = prime - root1; 

		amodp = (int)mpz_tdiv_ui(poly->mpz_poly_a,prime);
		bmodp = (int)mpz_tdiv_ui(poly->mpz_poly_b,prime);

		//find a^-1 mod p = inv(a mod p) mod p
		inv = modinv_1(amodp,prime);

		COMPUTE_FIRST_ROOTS

        // Barrett reduction for inv * root
        t2 = (uint64)inv * (uint64)root1;
        tmp = (t2 * fb->list->binv[i]) >> 32;
        t2 -= tmp * prime;
        root1 = (t2 >= prime) ? t2 - prime : t2;

        t2 = (uint64)inv * (uint64)root2;
        tmp = (t2 * fb->list->binv[i]) >> 32;
        t2 -= tmp * prime;
        root2 = (t2 >= prime) ? t2 - prime : t2;

		if (root2 < root1)
		{
			update_data.sm_firstroots1[i] = (uint16)root2;
			update_data.sm_firstroots2[i] = (uint16)root1;

			fb_p->root1[i] = (uint16)root2;
			fb_p->root2[i] = (uint16)root1;
			fb_n->root1[i] = (uint16)(prime - root1);
			fb_n->root2[i] = (uint16)(prime - root2);
		}
		else
		{
			update_data.sm_firstroots1[i] = (uint16)root1;
			update_data.sm_firstroots2[i] = (uint16)root2;

			fb_p->root1[i] = (uint16)root1;
			fb_p->root2[i] = (uint16)root2;
			fb_n->root1[i] = (uint16)(prime - root2);
			fb_n->root2[i] = (uint16)(prime - root1);
		}

		//for this factor base prime, compute the rootupdate value for all s
		//Bl values.  amodp holds a^-1 mod p
		//the rootupdate value is given by 2*Bj*amodp
		//Bl[j] now holds 2*Bl
		for (j=0;j<s;j++)
		{
			x = (int)mpz_tdiv_ui(dconf->Bl[j],prime);

			// x * inv % prime
            t2 = (uint64)inv * (uint64)x;
            tmp = (t2 * fb->list->binv[i]) >> 32;
            t2 -= tmp * prime;
            x = (t2 >= prime) ? t2 - prime : t2;

			rootupdates[(j)*fb->B+i] = x;
			dconf->sm_rootupdates[(j)*fb->med_B+i] = (uint16)x;
		}
	}

	//printf("prime[15bit-1] = %u\n", fb_p->prime[fb->fb_15bit_B-1]);
	for (i=fb->fb_15bit_B;i<fb->med_B;i++)
	{
		uint64 tmp, t2;

		prime = fb->list->prime[i];
		root1 = modsqrt[i]; 
		root2 = prime - root1; 

		amodp = (int)mpz_tdiv_ui(poly->mpz_poly_a,prime);
		bmodp = (int)mpz_tdiv_ui(poly->mpz_poly_b,prime);

		//find a^-1 mod p = inv(a mod p) mod p
		inv = modinv_1(amodp,prime);

		COMPUTE_FIRST_ROOTS

        // Barrett reduction for inv * root
        t2 = (uint64)inv * (uint64)root1;
        tmp = (t2 * fb->list->binv[i]) >> 32;
        t2 -= tmp * prime;
        root1 = (t2 >= prime) ? t2 - prime : t2;

        t2 = (uint64)inv * (uint64)root2;
        tmp = (t2 * fb->list->binv[i]) >> 32;
        t2 -= tmp * prime;
        root2 = (t2 >= prime) ? t2 - prime : t2;

		if (root2 < root1)
		{
			update_data.sm_firstroots1[i] = (uint16)root2;
			update_data.sm_firstroots2[i] = (uint16)root1;

			fb_p->root1[i] = (uint16)root2;
			fb_p->root2[i] = (uint16)root1;
			fb_n->root1[i] = (uint16)(prime - root1);
			fb_n->root2[i] = (uint16)(prime - root2);
		}
		else
		{
			update_data.sm_firstroots1[i] = (uint16)root1;
			update_data.sm_firstroots2[i] = (uint16)root2;

			fb_p->root1[i] = (uint16)root1;
			fb_p->root2[i] = (uint16)root2;
			fb_n->root1[i] = (uint16)(prime - root2);
			fb_n->root2[i] = (uint16)(prime - root1);
		}

		//for this factor base prime, compute the rootupdate value for all s
		//Bl values.  amodp holds a^-1 mod p
		//the rootupdate value is given by 2*Bj*amodp
		//Bl[j] now holds 2*Bl
		for (j=0;j<s;j++)
		{
			x = (int)mpz_tdiv_ui(dconf->Bl[j],prime);

			// x * inv % prime
            t2 = (uint64)inv * (uint64)x;
            tmp = (t2 * fb->list->binv[i]) >> 32;
            t2 -= tmp * prime;
            x = (t2 >= prime) ? t2 - prime : t2;

			rootupdates[(j)*fb->B+i] = x;
			dconf->sm_rootupdates[(j)*fb->med_B+i] = (uint16)x;
		}
	}

	check_bound = fb->med_B + BUCKET_ALLOC/2;
	logp = fb->list->logprime[fb->med_B-1];
	for (i=fb->med_B;i<fb->large_B;i++)
	{

#ifdef DO_VLP_OPT
        if (i >= check_bound)
        {
            room = 0;
            /* find the most filled bucket */
            for (k = 0; k < numblocks; k++)
            {
                if (*(numptr_p + k) > room)
                    room = *(numptr_p + k);
                if (*(numptr_n + k) > room)
                    room = *(numptr_n + k);
            }
            room = BUCKET_ALLOC - room;

            /* if it is filled close to the allocation, start recording in a new set of buckets */
            if (room < 32)
            {
                //printf("firstroots: bucket full, now at fb index %d, starting new slice %d\n",
                //    i, bound_index + 1);
                //uint32 ii;
                //uint32 bb;
                //for (bb = 0; bb < numblocks; bb++)
                //{
                //    printf("%u lp p-roots in slice %d, block %d:\n", numptr_p[bb], bound_index, bb);
                //    for (ii = 0; ii < numptr_p[bb]; ii++)
                //    {
                //        //bptr = sliceptr_p + ((uint64)bb << BLOCKBITS) + (uint64)ii;
                //        //printf("%08x ", *bptr);
                //        printf("%08x ", sliceptr_p[bb * BUCKET_ALLOC + ii]);
                //    }
                //    printf("\n");
                //}
                //
                //for (bb = 0; bb < numblocks; bb++)
                //{
                //    printf("%u lp n-roots in slice %d, block %d:\n", numptr_n[bb], bound_index, bb);
                //    for (ii = 0; ii < numptr_n[bb]; ii++)
                //    {
                //        //bptr = sliceptr_n + ((uint64)bb << BLOCKBITS) + (uint64)ii;
                //        //printf("%08x ", *bptr);
                //        printf("%08x ", sliceptr_n[bb * BUCKET_ALLOC + ii]);
                //    }
                //    printf("\n");
                //}
                logp = update_data.logp[i];
                lp_bucket_p->logp[bound_index] = logp;
                bound_index++;
                lp_bucket_p->fb_bounds[bound_index] = i;
                bound_val = i;
                sliceptr_p += (numblocks << (BUCKET_BITS + 1));
                sliceptr_n += (numblocks << (BUCKET_BITS + 1));
                numptr_p += (numblocks << 1);
                numptr_n += (numblocks << 1);
                check_bound += BUCKET_ALLOC >> 1;
            }
            else
            {
                check_bound += room >> 1;
            }
        }
        else if ((i - bound_val) >= 65536)
        {
            //printf("firstroots: prime slice limit, starting new slice %d\n",
            //    bound_index + 1);
            //int ii;
            //int bb;
            //for (bb = 0; bb < numblocks; bb++)
            //{
            //    printf("lp p-roots in slice %d, block %d:\n", bound_index, bb);
            //    for (ii = 0; ii < numptr_p[bb]; ii++)
            //    {
            //        printf("%u ", sliceptr_p[bb * BUCKET_ALLOC + ii]);
            //    }
            //    printf("\n");
            //}
            //
            //for (bb = 0; bb < numblocks; bb++)
            //{
            //    printf("lp n-roots in slice %d, block %d:\n", bound_index, bb);
            //    for (ii = 0; ii < numptr_n[bb]; ii++)
            //    {
            //        printf("%u ", sliceptr_n[bb * BUCKET_ALLOC + ii]);
            //    }
            //    printf("\n");
            //}
            lp_bucket_p->logp[bound_index] = logp;
            bound_index++;
            lp_bucket_p->fb_bounds[bound_index] = i;
            bound_val = i;
            sliceptr_p += (numblocks << (BUCKET_BITS + 1));
            sliceptr_n += (numblocks << (BUCKET_BITS + 1));
            numptr_p += (numblocks << 1);
            numptr_n += (numblocks << 1);
            check_bound += BUCKET_ALLOC >> 1;
        }
#else

        CHECK_NEW_SLICE(i);

#endif

		prime = fb->list->prime[i];
		root1 = modsqrt[i];
		root2 = prime - root1; 

		amodp = (int)mpz_tdiv_ui(poly->mpz_poly_a,prime);
		bmodp = (int)mpz_tdiv_ui(poly->mpz_poly_b,prime);

		//find a^-1 mod p = inv(a mod p) mod p
		inv = modinv_1(amodp,prime);

        COMPUTE_FIRST_ROOTS

        root1 = (uint32)((uint64)inv * (uint64)root1 % (uint64)prime);
        root2 = (uint32)((uint64)inv * (uint64)root2 % (uint64)prime);
		
		update_data.firstroots1[i] = root1;
		update_data.firstroots2[i] = root2;
        nroot1 = (prime - root1);
        nroot2 = (prime - root2);

		FILL_ONE_PRIME_LOOP_P(i);
		FILL_ONE_PRIME_LOOP_N(i);

		//for this factor base prime, compute the rootupdate value for all s
		//Bl values.  amodp holds a^-1 mod p
		//the rootupdate value is given by 2*Bj*amodp
		//Bl[j] now holds 2*Bl
		for (j=0;j<s;j++)
		{
			x = (int)mpz_tdiv_ui(dconf->Bl[j], prime);
			x = (int)((int64)x * (int64)inv % (int64)prime);

			rootupdates[(j)*fb->B+i] = x;
		}

	}

	logp = fb->list->logprime[fb->large_B-1];

#ifdef DO_VLP_OPT

    if (lp_bucket_p->list != NULL)
    {
        //printf("starting vlp with new slice %d\n", bound_index + 1);
        //for (i = 0, root1 = 0; i < 2 * numblocks; i++)
        //{
        //    root1 += numptr_p[i];
        //}
        //printf("slice %d had %d total entries (p and n)\n", bound_index, root1);
        // force a new slice
        i = fb->large_B;
        logp = update_data.logp[i];
        lp_bucket_p->logp[bound_index] = logp;
        bound_index++;
        // signal to large_sieve that we are now doing VLP_OPT.
        // we are not storing prime indices anyway, so this value 
        // is useless for lp_tdiv.
        lp_bucket_p->fb_bounds[bound_index] = 0x80000000 | i;
        bound_val = i;
        sliceptr_p += (numblocks << (BUCKET_BITS + 1));
        sliceptr_n += (numblocks << (BUCKET_BITS + 1));
        numptr_p += (numblocks << 1);
        numptr_n += (numblocks << 1);
        check_bound = bound_val + ((BUCKET_ALLOC / 2 * numblocks) >> 1);

        sliceptr64_p = (uint64*)sliceptr_p;
        sliceptr64_n = (uint64*)sliceptr_n;
        //printf("initial numptr_p[0] = %u, numptr_n[0] = %u\n", numptr_p[0], numptr_n[0]);
        //printf("bound_val = %u, check_bound = %u\n", bound_val, check_bound);
    }

    for (i = fb->large_B; i < fb->B; i++)
    {
        //CHECK_NEW_SLICE(i);

        if (i >= check_bound)
        {
            //printf("checking fullness of slice, num_p,num_n = %u,%u\n", numptr_p[0], numptr_n[0]);
            room = (BUCKET_ALLOC / 2 * numblocks) - MAX(numptr_p[0], numptr_n[0]);
            /* if it is filled close to the allocation, start recording in a new set of buckets */
            if (room < 32)
            {
                //printf("firstroots: bucket full, now at fb index %d, added %u,%u elements, starting new slice %d\n",
                //    i, numptr_p[0], numptr_n[0], bound_index + 1);
                //int ii;
                //printf("vlp p-roots in slice %d:\n", bound_index);
                //for (ii = 0; ii < numptr_p[0]; ii++)
                //{
                //    bptr = sliceptr_p + (uint64)ii;
                //    printf("%08x ", *bptr);
                //}
                //printf("\n");
                //printf("vlp n-roots in slice %d:\n", bound_index);
                //for (ii = 0; ii < numptr_n[0]; ii++)
                //{
                //    bptr = sliceptr_n + (uint64)ii;
                //    printf("%08x ", *bptr);
                //}
                //printf("\n");
                logp = update_data.logp[i];
                lp_bucket_p->logp[bound_index] = logp;
                bound_index++;
                lp_bucket_p->fb_bounds[bound_index] = i;
                bound_val = i;
                sliceptr_p += (numblocks << (BUCKET_BITS + 1));
                sliceptr_n += (numblocks << (BUCKET_BITS + 1));
                numptr_p += (numblocks << 1);
                numptr_n += (numblocks << 1);
                check_bound += ((BUCKET_ALLOC / 2 * numblocks) >> 1);

                sliceptr64_p = (uint64*)sliceptr_p;
                sliceptr64_n = (uint64*)sliceptr_n;
            }
            else
            {
                check_bound += room >> 1;
                //printf("have room for %u more elements, setting check_bound to %u\n", room, check_bound);
            }
        }
        else if ((i - bound_val) >= 65536)
        {
            //printf("firstroots: prime slice limit, index %d, added %u,%u elements, starting new slice %d\n",
            //    i, numptr_p[0], numptr_n[0], bound_index + 1);
            //int ii;
            //printf("vlp p-roots in slice %d:\n", bound_index);
            //for (ii = 0; ii < numptr_p[0]; ii++)
            //{
            //    printf("%u ", sliceptr_p[numptr_p[0] + ii]);
            //}
            //printf("\n");
            //printf("vlp n-roots in slice %d:\n", bound_index);
            //for (ii = 0; ii < numptr_n[0]; ii++)
            //{
            //    printf("%u ", sliceptr_n[numptr_n[0] + ii]);
            //}
            //printf("\n");
            lp_bucket_p->logp[bound_index] = logp;
            bound_index++;
            lp_bucket_p->fb_bounds[bound_index] = i;
            bound_val = i;
            sliceptr_p += (numblocks << (BUCKET_BITS + 1));
            sliceptr_n += (numblocks << (BUCKET_BITS + 1));
            numptr_p += (numblocks << 1);
            numptr_n += (numblocks << 1);
            check_bound += ((BUCKET_ALLOC / 2 * numblocks) >> 1);

            sliceptr64_p = (uint64*)sliceptr_p;
            sliceptr64_n = (uint64*)sliceptr_n;
        }

        prime = fb->list->prime[i];
        root1 = modsqrt[i];
        root2 = prime - root1;

        amodp = (int)mpz_tdiv_ui(poly->mpz_poly_a, prime);
        bmodp = (int)mpz_tdiv_ui(poly->mpz_poly_b, prime);

        //find a^-1 mod p = inv(a mod p) mod p
        inv = modinv_1(amodp, prime);

        COMPUTE_FIRST_ROOTS

        root1 = (uint32)((uint64)inv * (uint64)root1 % (uint64)prime);
        root2 = (uint32)((uint64)inv * (uint64)root2 % (uint64)prime);

        update_data.firstroots1[i] = root1;
        update_data.firstroots2[i] = root2;

        //FILL_ONE_PRIME_P(i);
        if (root1 < interval) 
            sliceptr64_p[numptr_p[0]++] = ((uint64)(i - bound_val) << 32) | root1;
        if (root2 < interval) 
            sliceptr64_p[numptr_p[0]++] = ((uint64)(i - bound_val) << 32) | root2;

        root1 = (prime - root1);
        root2 = (prime - root2);

        FILL_ONE_PRIME_N(i);
        if (root1 < interval) 
            sliceptr64_n[numptr_n[0]++] = ((uint64)(i - bound_val) << 32) | root1;
        if (root2 < interval) 
            sliceptr64_n[numptr_n[0]++] = ((uint64)(i - bound_val) << 32) | root2;

        //for this factor base prime, compute the rootupdate value for all s
        //Bl values.  amodp holds a^-1 mod p
        //the rootupdate value is given by 2*Bj*amodp
        //Bl[j] now holds 2*Bl
        //s is the number of primes in 'a'
        //fprintf(out, "%u,%u,%u", prime, root1, root2);
        for (j = 0; j < s; j++)
        {
            x = (int)mpz_tdiv_ui(dconf->Bl[j], prime);
            x = (int)((int64)x * (int64)inv % (int64)prime);
            rootupdates[(j)* fb->B + i] = x;
            //fprintf(out, ",%u", x);
        }
        //fprintf(out, "\n");
    }

    //printf("firstroots done\n");

#else

    for (i = fb->large_B; i < fb->B; i++)
    {
        CHECK_NEW_SLICE(i);

        prime = fb->list->prime[i];
        root1 = modsqrt[i];
        root2 = prime - root1;

        amodp = (int)mpz_tdiv_ui(poly->mpz_poly_a, prime);
        bmodp = (int)mpz_tdiv_ui(poly->mpz_poly_b, prime);

        //find a^-1 mod p = inv(a mod p) mod p
        inv = modinv_1(amodp, prime);

        COMPUTE_FIRST_ROOTS;

        root1 = (uint32)((uint64)inv * (uint64)root1 % (uint64)prime);
        root2 = (uint32)((uint64)inv * (uint64)root2 % (uint64)prime);

        update_data.firstroots1[i] = root1;
        update_data.firstroots2[i] = root2;

        FILL_ONE_PRIME_P(i);

        root1 = (prime - root1);
        root2 = (prime - root2);

        FILL_ONE_PRIME_N(i);

        //for this factor base prime, compute the rootupdate value for all s
        //Bl values.  amodp holds a^-1 mod p
        //the rootupdate value is given by 2*Bj*amodp
        //Bl[j] now holds 2*Bl
        //s is the number of primes in 'a'
        //fprintf(out, "%u,%u,%u", prime, root1, root2);
        for (j = 0; j < s; j++)
        {
            x = (int)mpz_tdiv_ui(dconf->Bl[j], prime);
            x = (int)((int64)x * (int64)inv % (int64)prime);
            rootupdates[(j)* fb->B + i] = x;
            //fprintf(out, ",%u", x);
        }
        //fprintf(out, "\n");
    }

#endif
	

    if (0)
    {
        char *v = dconf->curr_poly->nu;
        char *sign = dconf->curr_poly->gray;
        double avg = 0.;
        int hit;
        uint32 pstarti;

        printf("checking large prime hits in the interval of size %u\n", interval);

        pstarti = fb->large_B;
        for (i = fb->large_B; i < fb->B; i++)
        {
            int j;

            prime = fb->list->prime[i];
            root1 = update_data.firstroots1[i];
            root2 = update_data.firstroots2[i];

            for (j = 1; j < dconf->maxB; j++)
            {
                hit = 0;
                if (sign[j] > 0)
                {
                    root1 = (int)root1 - rootupdates[(v[j] - 1) * fb->B + i];
                    root2 = (int)root2 - rootupdates[(v[j] - 1) * fb->B + i];
                    root1 += (root1 < 0) * prime;
                    root1 += (root2 < 0) * prime;

                    if (root1 < interval)
                    {
                        hit = 1;
                    }

                    if (root2 < interval)
                    {
                        hit = 1;
                    }

                    if ((prime - root1) < interval)
                    {
                        hit = 1;
                    }

                    if ((prime - root2) < interval)
                    {
                        hit = 1;
                    }
                }
                else
                {
                    root1 = (int)root1 + rootupdates[(v[j] - 1) * fb->B + i];
                    root2 = (int)root2 + rootupdates[(v[j] - 1) * fb->B + i];
                    root1 -= ((root1 >= prime) * prime);
                    root2 -= ((root2 >= prime) * prime);

                    if (root1 < interval)
                    {
                        hit = 1;
                    }

                    if (root2 < interval)
                    {
                        hit = 1;
                    }

                    if ((prime - root1) < interval)
                    {
                        hit = 1;
                    }

                    if ((prime - root2) < interval)
                    {
                        hit = 1;
                    }
                }

                avg += (double)hit;
            }

            if ((i - pstarti) >= 1000)
            {
                printf("primes from %u(%d) to %u(%d) hit %1.0f times in all polys",
                    fb->list->prime[pstarti], pstarti, fb->list->prime[i], i, avg);

                avg /= (i - pstarti);
                avg /= dconf->maxB;

                printf(" (%1.2f%% hit rate of %u primes in %u polys (any root or side)\n",
                    avg * 100, (i - pstarti), dconf->maxB);
                
                avg = 0.;
                pstarti = i;
            }
        }

        exit(1);
    }

    

	if (lp_bucket_p->list != NULL)
		lp_bucket_p->num_slices = bound_index + 1;
	

	return;
}

