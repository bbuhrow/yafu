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

#include "qs_impl.h"
#include "ytools.h"
#include "common.h"
#include "poly_macros_common.h"
#include "poly_macros_32k.h"

//#define DEBUG_SS_SEARCH

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

//#define SS_TIMING

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

	uint32_t i,logp;
	int root1, root2, prime, amodp, bmodp, inv, bnum,numblocks;
	int lpnum,last_bound;

	// unpack stuff from the job data
	siqs_poly *poly = dconf->curr_poly;
	fb_list *fb = sconf->factor_base;
	lp_bucket *lp_bucket_p = dconf->buckets;
	uint32_t *modsqrt = sconf->modsqrt_array;

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
	
		root1 = (uint32_t)((uint64_t)inv * (uint64_t)root1 % (uint64_t)prime);
		root2 = (uint32_t)((uint64_t)inv * (uint64_t)root2 % (uint64_t)prime);

		// just need to do this once, because the next step of prime will be 
		// into a different bucket
		bnum = root1 >> 15;
		if (bnum == 0)
			lpnum++;

		// repeat for the other root
		bnum = root2 >> 15;
		if (bnum == 0)
			lpnum++;

		if ((uint32_t)lpnum > (double)BUCKET_ALLOC * 0.75)
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

// Subset-sum notes:
// siqs initialization amounts to fast ways to compute 
// x == a^-1 * (+/- t - (+/- B1 +/- B2 +/- B3 ... + Bs) mod p) mod p
// for all of the combinations of the B's, for each p.
// jasonp's subset-sum search method uses a time/space tradeoff to
// the sum.  With 's' B-values there are 2^(s-1) combinations of sums of B's.
// Instead of enumerating all of the sum combinations we split the list of B's
// and enumerate each half mod p (including the a^-1 and +/- t on one side).
// (sqrt factor of the number of sum enumerations.)
// These are set1 and set2.
// Now hash each set into a group of bins along the interval 0:prime, so
// that each enumeration mod p gets grouped with others about the same size.
// (sorting, essentially.  perhaps a sort would even be faster than a hash.)
// (... it doesn't matter much, this part of the algorithm is not the bottlenectk.)
// For each bin in set 1, compute the center of the set, x, and add each
// member of that bin to the bin in set 2 that contains prime - x.  So that
// the x's mostly cancel and we are left with something near in value to prime, 
// which taken mod p is near 0, i.e., likely to be very small and thus produce a 
// sieve hit for the polynomial represented by the enumeration in set1+set2.
// (Note that there are also hits when combining bin x with (p-x)+c and
// when combining bin x with (p-x)-c, for some value of 'c', depending on the size
// of the bins and simply because some binned values lie near the edge of the bin.)
// finally there is a choice of what to do with the resulting sieve hits.
// We could 1) sort the list of all sieve hits by their associated enumeration #,
// so that we can sieve "normally" by all small primes for each polynomial, and
// afterwards dump in the large prime sieve hits for that polynomial.  
// Sorting could occur using 
// 1a) a literal sort: put the sieve hits into a giant linear array as we
// find them and then sort it.
// 1b) a bucket sort: insert each sieve hit into its appropriate polynomial-bin
// as we find them.
// 2) put the sieve hits into linked lists (per polynomial) and traverse the lists
// when sieving.
// Note that we can still use all of the traditional sieve algorithms that
// yafu uses, including bucket sieving.  This new subset-sum technique would
// take over at some large factor base offset, replacing the bucket sieve for
// the largest set of fb primes.  Then we just have to find some parameterization
// (fb size, interval size, td cutoffs, SS start offset, etc.) that hopefully
// results in an overall faster algorithm for usefully sized N (i.e., smaller
// than the current qs-gnfs crossover).


// function for qsort to sort ss-elements by polynum
int cmp_ss_element(const void* x, const void* y)
{
	uint64_t* xx = (uint64_t*)x;
	uint64_t* yy = (uint64_t*)y;

	uint16_t polynumx = ((*xx) >> 40) & 0xffffff;
	uint16_t polynumy = ((*yy) >> 40) & 0xffffff;

	if (polynumx > polynumy)
		return 1;
	else if (polynumx == polynumy)
		return 0;
	else
		return -1;
}

#ifdef USE_LINKED_LIST_SS
void update_ll(int side, ss_slice_t *buckets, int slice, int idx, dynamic_conf_t* dconf)
{
	// add poly idx to the linked list slice on the indicated side.
	uint64_t lnext;
	uint64_t* poly_ll_ptr;

	if (side == 0)
	{
		poly_ll_ptr = dconf->poly_ll_ptr_p;
	}
	else
	{
		poly_ll_ptr = dconf->poly_ll_ptr_n;
	}

	if (poly_ll_ptr[idx] == (uint64_t)(-1))
	{
		// first member of the LL for this poly.
		// record the current slice and location within
		// the slice of the sieve hit for this poly.
		poly_ll_ptr[idx] = ((uint64_t)slice << 32) |
			(uint64_t)buckets[slice].size;

		if (side == 0)
		{
			dconf->poly_ll_first_p[idx] = buckets[slice].size;
		}
		else
		{
			dconf->poly_ll_first_n[idx] = buckets[slice].size;
		}
	}
	else
	{
		// unpack the slice and location within the slice of
		// the previous hit for this poly.
		uint32_t s = poly_ll_ptr[idx] >> 32;
		uint32_t loc = poly_ll_ptr[idx] & 0xffffffff;

		if (slice == s)
		{
			// if the current slice is the same record the distance
			// in the previous entry to point to this one.
			lnext = buckets[slice].size - loc;
			if (lnext >= (1 << 24))
			{
				printf("error: prevous entry for slice %d poly %d too large (%d)\n",
					slice, idx, (int)lnext);
				exit(1);
			}
			lnext <<= 40;
			buckets[s].element[loc] &= 0x000000ffffffffff;
			buckets[s].element[loc] |= lnext;
		}
		else
		{
			// if the current slice is not the same then
			// 1) check that it is at most 1 slice different.
			// 2) record the distance to the end of the 
			// previous slice plus the distance into this one.
			if ((slice - s) > 1)
			{
				printf("error: slice distance too far (%d)\n", slice - s);
				exit(1);
			}
			lnext = buckets[s].size - loc;
			lnext += buckets[slice].size;
			if (lnext >= (1 << 24))
			{
				printf("error: prevous entry for slice %d poly %d too large (%d)\n",
					slice, idx, (int)lnext);
				exit(1);
			}
			lnext <<= 40;
			buckets[s].element[loc] &= 0x000000ffffffffff;
			buckets[s].element[loc] |= lnext;
		}

		// record the current slice and location within
		// the slice of the sieve hit for this poly.
		poly_ll_ptr[idx] = ((uint64_t)slice << 32) |
			(uint64_t)buckets[slice].size;
	}
	return;
}

void ss_search_linked_lists(static_conf_t* sconf, dynamic_conf_t* dconf)
{
	// the subset-sum search algorithm using linked lists to
	// organize the sieve-hits per poly.  
	// Linked lists makes this portion of the algorithm very 
	// fast because it is cheap to add to a list.  But sieving
	// then gets very expensive because the list must be traversed
	// and there is no good way to accelerate that.
	siqs_poly* poly = dconf->curr_poly;
	fb_list* fb = sconf->factor_base;
	uint32_t numblocks = sconf->num_blocks;
	uint32_t interval;

	numblocks = sconf->num_blocks;
	interval = numblocks << 15;

#ifdef USE_SS_SEARCH

	int i, j, ii;
	int ss_num_poly_terms = poly->s / 2;
	ss_set_t ss_set1a, ss_set1b, ss_set2a, ss_set2b;

	mpz_t polyb1, polyb2;
	mpz_init(polyb1);
	mpz_init(polyb2);



#ifdef DEBUG_SS_SEARCH
	printf("allocating %u bytes for set1a/b of size %d\n",
		(1 << ss_num_poly_terms) * sizeof(int) * 2 * 2, ss_num_poly_terms);
#endif
	ss_set1a.root = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
	ss_set1a.polynum = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
	ss_set1a.size = ss_num_poly_terms;

	ss_set1b.root = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
	ss_set1b.polynum = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
	ss_set1b.size = ss_num_poly_terms;

	//gmp_printf("%Zd\n", dconf->Bl[0]);
	mpz_set(polyb1, dconf->Bl[0]);

	for (ii = 1; ii < ss_num_poly_terms; ii++) {
		//gmp_printf("%Zd\n", dconf->Bl[ii]);
		mpz_add(polyb1, polyb1, dconf->Bl[ii]);
	}
	mpz_tdiv_q_2exp(polyb1, polyb1, 1);

	ss_num_poly_terms = (poly->s - 1) - ss_num_poly_terms;

#ifdef DEBUG_SS_SEARCH
	printf("allocating %u bytes for set2a/b of size %d\n",
		(1 << ss_num_poly_terms) * sizeof(int) * 2 * 2, ss_num_poly_terms);
#endif
	ss_set2a.root = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
	ss_set2a.polynum = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
	ss_set2a.size = ss_num_poly_terms;

	ss_set2b.root = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
	ss_set2b.polynum = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
	ss_set2b.size = ss_num_poly_terms;

	// first poly b: sum of first n positive Bl
	//gmp_printf("%Zd\n", dconf->Bl[size1]);
	mpz_set(polyb2, dconf->Bl[ii]);

	for (ii++; ii < poly->s; ii++) {
		//gmp_printf("%Zd\n", dconf->Bl[size1 + ii]);
		mpz_add(polyb2, polyb2, dconf->Bl[ii]);
	}
	mpz_tdiv_q_2exp(polyb2, polyb2, 1);

	int numbins = 2 * fb->list->prime[fb->B - 1] / (interval)+1;
	int bindepth = 128;
	int binsize; // = p / numbins + 1;

	ss_set_t* bins1;
	ss_set_t* bins2;

	// init bins
	bins1 = (ss_set_t*)xmalloc(numbins * sizeof(ss_set_t));
	bins2 = (ss_set_t*)xmalloc(numbins * sizeof(ss_set_t));

	for (ii = 0; ii < numbins; ii++)
	{
		bins1[ii].root = (int*)xmalloc(bindepth * sizeof(int));
		bins2[ii].root = (int*)xmalloc(bindepth * sizeof(int));
		bins1[ii].polynum = (int*)xmalloc(bindepth * sizeof(int));
		bins2[ii].polynum = (int*)xmalloc(bindepth * sizeof(int));
		bins1[ii].alloc = bindepth;
		bins2[ii].alloc = bindepth;
		bins1[ii].size = 0;
		bins2[ii].size = 0;
	}

	//printf("allocated %d bins of depth %d\n", numbins, bindepth);
	//exit(0);

#ifdef DEBUG_SS_SEARCH
	mpz_add(dconf->gmptmp1, polyb2, polyb1);

	if (sconf->knmod8 == 1)
	{
		if (sconf->charcount == 1)
		{
			printf("including polya in polyb sum\n");
			mpz_add(dconf->gmptmp1, dconf->gmptmp1, poly->mpz_poly_a);
		}
	}

	if (mpz_cmp(dconf->gmptmp1, poly->mpz_poly_b) == 0)
	{
		prime = fb->list->prime[fb->B - 1];
		gmp_printf("polyb = %Zd + %Zd = %Zd\npolyb matches! (%u, %u, %u, %u)\n",
			polyb1, polyb2, dconf->gmptmp1,
			prime,
			mpz_tdiv_ui(polyb1, prime),
			mpz_tdiv_ui(polyb2, prime),
			mpz_tdiv_ui(poly->mpz_poly_a, prime));
	}
	else
	{
		gmp_printf("polyb = %Zd\npolyb doesn't match :(\n", dconf->gmptmp1);
	}
#endif

#if defined( DEBUG_SS_SEARCH) || defined(DEBUG_SS_SEARCH2)
	ss_set_t orighits;
	orighits.root = (int*)xmalloc(512 * sizeof(int));
	orighits.polynum = (int*)xmalloc(512 * sizeof(int));
#endif

	// create a map from poly number back to 
	// gray code enumeration order.
	int polynum = 0;
	int* polymap = (int*)xmalloc((1 << (poly->s - 1)) * sizeof(int));
	int* polynums = (int*)xmalloc((1 << (poly->s - 1)) * sizeof(int));
	int* polyv = (int*)xmalloc((1 << (poly->s - 1)) * sizeof(int));
	int* polysign = (int*)xmalloc((1 << (poly->s - 1)) * sizeof(int));

	polymap[0] = 1;
	polynums[0] = 0;

	//printf("%d <- %d (%d)\n", 1, 0, polynums[0]);

	for (ii = 1, polynum = 0; ii < (1 << (poly->s - 1)); ii++)
	{
		int v;
		int tmp;
		int sign;

		// next poly enumeration
		v = 1;
		j = ii;
		while ((j & 1) == 0)
			v++, j >>= 1;
		tmp = ii + (1 << v) - 1;
		tmp = (tmp >> v);
		polyv[ii] = v;
		if (tmp & 1)
			polysign[ii] = -1;
		else
			polysign[ii] = 1;

		// next polynum
		polynum ^= (1 << (v - 1));
		polynums[ii] = polynum;
		polymap[polynum] = ii + 1;

		//printf("%d <- %d (%d)\n", ii + 1, polynum, polynums[ii]);
	}

	//exit(0);

	uint64_t* poly_ll_ptr_p = dconf->poly_ll_ptr_p;
	uint64_t* poly_ll_ptr_n = dconf->poly_ll_ptr_n;
	uint32_t num_bpoly = (1 << (dconf->curr_poly->s - 1));

	for (i = 0; i < dconf->num_ss_slices; i++)
	{
		dconf->ss_slices_p[i].size = 0;
		dconf->ss_slices_n[i].size = 0;
		dconf->ss_slices_p[i].curr_poly_idx = 0;
		dconf->ss_slices_n[i].curr_poly_idx = 0;
		dconf->ss_slices_p[i].curr_poly_num = 0;
		dconf->ss_slices_n[i].curr_poly_num = 0;
	}

	double t_enum_roots;
	double t_sort_roots;
	double t_match_roots;
	double t_sort_buckets;
	struct timeval start, stop;


	// todo
	// allocate buckets for polysum hits, one for each b-poly and per sieve side.
	// X stop bucket sieving at the ss_start_B index.
	// X set the ss_start_B index in static poly init.
	// need a bucket dumping routine somewhere in tdiv.  maybe tdiv_large.
	// testing!

	// optimizations:
	// can still use vector root updates when enumerating roots
	// for each side of the sum.  just need to be careful to
	// do root matching per prime in that case.
	// probably can do vector match search since we take an element
	// on one side then sum and compare to many elements on the other side.
	// could either do this many primes at a time or maybe many elements per
	// prime at a time.
	// how many "other side" bins we need to search for matches needs to
	// be more exact than hardcoded +/- 1 bins from the center bin.
	int* rootupdates = dconf->rootupdates;
	update_t update_data = dconf->update_data;
	uint32_t* modsqrt = sconf->modsqrt_array;

	for (i = fb->x2_large_B; i < fb->B; i++)
	{
		int numB;
		int maxB = 1 << (dconf->curr_poly->s - 1);
		int r1 = update_data.firstroots1[i], r2 = update_data.firstroots2[i];
		uint32_t bound = sconf->factor_base->B;
		uint32_t med_B = sconf->factor_base->med_B;
		int polynum = 0;
		int slice = (int)((i - sconf->factor_base->x2_large_B) / 65536);

		int prime = fb->list->prime[i];
		int p = prime;
		int root1 = modsqrt[i];
		int root2 = prime - root1;
		uint8_t logp = fb->list->logprime[i];

		int amodp = (int)mpz_tdiv_ui(poly->mpz_poly_a, prime);
		int bmodp = (int)mpz_tdiv_ui(poly->mpz_poly_b, prime);
		int inv;

		//find a^-1 mod p = inv(a mod p) mod p
		if (sconf->knmod8 == 1)
		{
			inv = modinv_1(2 * amodp, prime);
		}
		else
		{
			inv = modinv_1(amodp, prime);
		}

		//if (i == sconf->factor_base->x2_large_B)
		//{
		//	printf("new slice %d starts at fboffset %u\n", slice, i);
		//	dconf->ss_slices_p[slice].fboffset = i;
		//	dconf->ss_slices_n[slice].fboffset = i;
		//}
		//
		//if (slice > 0)
		//{
		//	int previdx_slice = (int)((i - 1 - sconf->factor_base->x2_large_B) / 65536);
		//	if (previdx_slice != slice)
		//	{
		//		printf("new slice %d starts at fboffset %u\n", slice, i);
		//		dconf->ss_slices_p[slice].fboffset = i;
		//		dconf->ss_slices_n[slice].fboffset = i;
		//	}
		//}

#if defined( DEBUG_SS_SEARCH) || defined(DEBUG_SS_SEARCH2)
		orighits.size = 0;

		if ((r1 < interval) || (r2 < interval) || ((prime - r1) < interval) || ((prime - r2) < interval))
		{
			//printf("%04d,%04x: initial roots = %d,%d (%d,%d)\n", 0, polynum, r1, r2, prime - r1, prime - r2);
			if (r1 < interval) {
				orighits.root[orighits.size] = r1;
				orighits.polynum[orighits.size] = polynum;
				orighits.size++;
			}
			if (r2 < interval) {
				orighits.root[orighits.size] = r2;
				orighits.polynum[orighits.size] = polynum;
				orighits.size++;
			}
			if ((prime - r1) < interval) {
				orighits.root[orighits.size] = prime - r1;
				orighits.polynum[orighits.size] = polynum;
				orighits.size++;
			}
			if ((prime - r2) < interval) {
				orighits.root[orighits.size] = prime - r2;
				orighits.polynum[orighits.size] = polynum;
				orighits.size++;
			}
		}

		for (numB = 1; numB < maxB; numB++)
		{
			char v = dconf->curr_poly->nu[numB];
			char sign = dconf->curr_poly->gray[numB];
			int* ptr = &rootupdates[(v - 1) * bound];
			char hit1, hit2, hit3, hit4;

			if (sign > 0)
			{
				r1 -= ptr[i];
				r2 -= ptr[i];
				r1 = (r1 < 0) ? r1 + prime : r1;
				r2 = (r2 < 0) ? r2 + prime : r2;
			}
			else
			{
				r1 += ptr[i];
				r2 += ptr[i];
				r1 = (r1 >= prime) ? r1 - prime : r1;
				r2 = (r2 >= prime) ? r2 - prime : r2;
			}

			polynum ^= (1 << (v - 1));

			hit1 = r1 < interval ? '*' : ' ';
			hit2 = r2 < interval ? '*' : ' ';
			hit3 = (prime - r1) < interval ? '*' : ' ';
			hit4 = (prime - r2) < interval ? '*' : ' ';

			if ((r1 < interval) || (r2 < interval) || ((prime - r1) < interval) || ((prime - r2) < interval))
			{
				//printf("%04d, %04x: ptr = %09d, v = %03d, s = %d: roots = %d%c,%d%c : %d%c,%d%c\n",
				//	numB, polynum, ptr[i], v, sign > 0 ? 1 : 0, r1, hit1, r2, hit2, prime - r1, hit3, prime - r2, hit4);

				if (r1 < interval) {
					orighits.root[orighits.size] = r1;
					orighits.polynum[orighits.size] = polynum;
					orighits.size++;
				}
				if (r2 < interval) {
					orighits.root[orighits.size] = r2;
					orighits.polynum[orighits.size] = polynum;
					orighits.size++;
				}
				if ((prime - r1) < interval) {
					orighits.root[orighits.size] = prime - r1;
					orighits.polynum[orighits.size] = polynum;
					orighits.size++;
				}
				if ((prime - r2) < interval) {
					orighits.root[orighits.size] = prime - r2;
					orighits.polynum[orighits.size] = polynum;
					orighits.size++;
				}

	}
			}

#endif

		// create set1, an enumeration of (+/-t - b)(a)^-1 mod p
		// for b in the set [+/-B1 +/- B2 +/- ... Bs/2]
		int ii, v, sign, n = poly->s / 2;
		int tmp;

		gettimeofday(&start, NULL);


		prime = fb->list->prime[i];
		root1 = modsqrt[i];
		root2 = prime - root1;

		//printf("new prime %d, r1 = %d, r2 = %d\n", prime, root1, root2);

#ifdef DEBUG_SS_SEARCH
			//printf("p,t1,t2 = %d,%d,%d\n", prime, root1, root2);
		printf("building set1 B\n");
#endif

		// first poly b: sum of first n positive Bl
		bmodp = (int)mpz_tdiv_ui(polyb1, prime);

		// under these conditions polya is added to polyb.
		// here we account for that by adding amodp into bmodp.
		// only here in set1 enumeration.
		if (sconf->knmod8 == 1)
		{
			if (sconf->charcount == 1)
			{
				bmodp += amodp;
				if (bmodp >= prime)
				{
					bmodp -= prime;
				}
			}
		}

#ifdef DEBUG_SS_SEARCH
		int sumbmodp = bmodp;
		mpz_set(dconf->gmptmp1, polyb1);
#endif

		//COMPUTE_FIRST_ROOTS;
		r1 = (int)root1 - bmodp;
		r2 = (int)root2 - bmodp;
		if (r1 < 0) r1 += prime;
		if (r2 < 0) r2 += prime;

		r1 = (int)((uint64_t)inv * (uint64_t)r1 % (uint64_t)prime);
		r2 = (int)((uint64_t)inv * (uint64_t)r2 % (uint64_t)prime);

#ifdef DEBUG_SS_SEARCH
		printf("\np = %d, a^-1 mod p = %d, roots = %d,%d, s = %d\n", prime, inv, r1, r2, n);
		gmp_printf("b = %Zd, amodp = %d, bmodp = %d\n", polyb1, amodp, bmodp);
#endif

		int* ptr;

#ifdef DEBUG_SS_SEARCH
		printf("enumerating %d elements of set 1 with %d polyb terms\n", (1 << n), n);
#endif
		// enumerate set1 roots
		for (ii = 1, polynum = 0; ii < (1 << n); ii++) {
			// these roots go into the set
			ss_set1a.root[ii - 1] = r1;
			ss_set1b.root[ii - 1] = r2;
			ss_set1a.polynum[ii - 1] = polynum;
			ss_set1b.polynum[ii - 1] = polynum;

			// next polynum
			polynum = polynums[ii];
			sign = polysign[ii];
			v = polyv[ii];

			// next roots
			ptr = &rootupdates[(v - 1) * bound + i];
			if (sign > 0)
			{
				r1 = (int)r1 - *ptr;
				r2 = (int)r2 - *ptr;
				r1 = (r1 < 0) ? r1 + prime : r1;
				r2 = (r2 < 0) ? r2 + prime : r2;
			}
			else
			{
				r1 = (int)r1 + *ptr;
				r2 = (int)r2 + *ptr;
				r1 = (r1 >= prime) ? r1 - prime : r1;
				r2 = (r2 >= prime) ? r2 - prime : r2;
			}
		}
		// these roots go into the set
		ss_set1a.root[ii - 1] = r1;
		ss_set1a.polynum[ii - 1] = polynum;
		ss_set1b.root[ii - 1] = r2;
		ss_set1b.polynum[ii - 1] = polynum;

		int size1 = n;

		// create set2, an enumeration of (+/-t - b)(a)^-1 mod p
		// for b in the set [+/-Bs/2+1 +/- Bs/2+2 ... + Bs]
		n = (poly->s - 1) - n;

#ifdef DEBUG_SS_SEARCH
		printf("building set2 B\n");
#endif
		bmodp = (int)mpz_tdiv_ui(polyb2, prime);

#ifdef DEBUG_SS_SEARCH
		sumbmodp += bmodp;
		if (sumbmodp >= prime) sumbmodp -= prime;
#endif

		// - (B1 + B2 + B3 + ... Bs-1)
		bmodp = prime - bmodp;

		r1 = (uint32_t)((uint64_t)inv * (uint64_t)bmodp % (uint64_t)prime);

#ifdef DEBUG_SS_SEARCH
		printf("\np = %d, a^-1 mod p = %d, roots = %d,%d, s = %d\n", prime, inv, r1, r2, n);
		gmp_printf("b = %Zd, amodp = %d, bmodp = %d\n", polyb2, amodp, bmodp);

		gmp_printf("orig polya: %Zd\n", poly->mpz_poly_a);
		gmp_printf("orig polyb: %Zd\n", poly->mpz_poly_b);
		int oamodp = mpz_tdiv_ui(poly->mpz_poly_a, prime);
		int obmodp = mpz_tdiv_ui(poly->mpz_poly_b, prime);
		printf("orig amodp = %d, bmodp = %d, root1 = %d, root2 = %d\n",
			oamodp, obmodp, update_data.firstroots1[i], update_data.firstroots2[i]);

		if (obmodp == sumbmodp)
		{
			printf("bmodp matches!\n");
		}
		else
		{
			printf("bmodp doesn't match :(\n");
		}


		printf("enumerating %d elements of set 2 with %d polyb terms\n", (1 << n), n);
#endif

		// enumerate set2 roots
		for (ii = 1, polynum = 0; ii < (1 << n); ii++) {
			// these roots go into the set
			ss_set2a.root[ii - 1] = r1;
			ss_set2a.polynum[ii - 1] = polynum;

			polynum = polynums[ii] << size1;
			sign = polysign[ii];
			v = polyv[ii];

			// next roots
			ptr = &rootupdates[(size1 + v - 1) * bound + i];
			if (sign > 0)
			{
				r1 = (int)r1 - *ptr;
				r1 = (r1 < 0) ? r1 + prime : r1;
			}
			else
			{
				r1 = (int)r1 + *ptr;
				r1 = (r1 >= prime) ? r1 - prime : r1;
			}
		}

		// these roots go into the set
		ss_set2a.root[ii - 1] = r1;
		ss_set2a.polynum[ii - 1] = polynum;


		gettimeofday(&stop, NULL);
		t_enum_roots += ytools_difftime(&start, &stop);

		gettimeofday(&start, NULL);

		// now sort the sets into a moderate number of bins over the range 0:p
		numbins = p / (2 * interval) + 1;
		bindepth = 128;
		binsize = p / numbins + 1;

#ifdef DEBUG_SS_SEARCH
		printf("sorting rootset1 of size %d into %d bins of size %d\n", (1 << ss_set1a.size), numbins, binsize);
#endif

		for (ii = 0; ii < numbins; ii++)
		{
			bins1[ii].size = 0;
			bins2[ii].size = 0;
		}

		// sort root 1 into set 1
		for (ii = 0; ii < (1 << ss_set1a.size); ii++)
		{
			int binnum = ss_set1a.root[ii] / binsize;
			if (binnum < numbins)
			{
				//printf("bin %d <-- set1a root %d (L polynum %d)\n", binnum, ss_set1a.root[ii], ss_set1a.polynum[ii]);
				bins1[binnum].root[bins1[binnum].size] = ss_set1a.root[ii];
				bins1[binnum].polynum[bins1[binnum].size] = ss_set1a.polynum[ii];
				bins1[binnum].size++;
				if (bins1[binnum].size > bindepth)
				{
					printf("bin overflow\n");
					exit(1);
				}
			}
			else
			{
				printf("element %d of set 1, root %d, invalid bin %d [of %d]\n",
					ii, ss_set1a.root[ii], binnum, numbins);
			}
		}

#ifdef DEBUG_SS_SEARCH
		printf("sorting rootset2 of size %d into %d bins of size %d\n", (1 << ss_set2a.size), numbins, binsize);
#endif

		// sort set 2
		for (ii = 0; ii < (1 << ss_set2a.size); ii++)
		{
			int binnum = ss_set2a.root[ii] / binsize;
			if (binnum < numbins)
			{
				//printf("bin %d <-- set2a root %d (R polynum %d)\n", binnum, ss_set2a.root[ii], ss_set2a.polynum[ii]);
				bins2[binnum].root[bins2[binnum].size] = ss_set2a.root[ii];
				bins2[binnum].polynum[bins2[binnum].size] = ss_set2a.polynum[ii];
				bins2[binnum].size++;
				if (bins2[binnum].size > bindepth)
				{
					printf("bin overflow\n");
					exit(1);
				}
			}
			else
			{
				printf("element %d of set 2, root %d, invalid bin %d [of %d]\n",
					ii, ss_set2a.root[ii], binnum, numbins);
			}
		}

		// now match all numbers in bin x with the other set's bin that
		// contains p - x
#ifdef DEBUG_SS_SEARCH
		printf("matching root1\n");
#endif


		gettimeofday(&stop, NULL);
		t_sort_roots += ytools_difftime(&start, &stop);

		gettimeofday(&start, NULL);

		// match set1 holding root1 to set2
		int nummatch = 0;
		int matchp1 = 0;
		int matchm1 = 0;
		for (ii = 0; ii < numbins; ii++)
		{
			int x = binsize * ii + binsize / 2;
			int px = p - x;
			int b = px / binsize; // px / binsize;

			nummatch = 0;
			matchp1 = 0;
			matchm1 = 0;

			//printf("bin %d (cent %d) with %d elements matching with bin %d with %d elements\n", 
			//	ii, x, bins1[ii].size, b, bins2[b].size);
			//
			//if ((b + 1) < numbins)
			//{
			//	printf("bin %d (cent %d) with %d elements matching with bin %d with %d elements\n",
			//		ii, x, bins1[ii].size, b + 1, bins2[b + 1].size);
			//}
			//
			//if ((b - 1) >= 0)
			//{
			//	printf("bin %d (cent %d) with %d elements matching with bin %d with %d elements\n",
			//		ii, x, bins1[ii].size, b - 1, bins2[b - 1].size);
			//}

			{
				for (j = 0; j < bins1[ii].size; j++)
				{
					int k;
					for (k = 0; k < bins2[b].size; k++)
					{
						int sum = bins1[ii].root[j] + bins2[b].root[k];
						int polysum = bins1[ii].polynum[j] + bins2[b].polynum[k];

						if (polysum >= (1 << (poly->s - 1)))
						{
							printf("invalid polysum %d (%d + %d)\n",
								polysum, bins1[ii].polynum[j], bins2[b].polynum[k]);
							continue;
						}

						if ((sum >= prime) && ((sum - prime) < interval))
						{
							int idx = polymap[polysum];

							if (idx == num_bpoly) continue;

							uint64_t pid = (uint64_t)(i - dconf->ss_slices_p[slice].fboffset) << 24;
							uint64_t sid = (uint64_t)(sum - prime);
							uint64_t lnext;

							update_ll(0, dconf->ss_slices_p, slice, idx, dconf);

							// the current entry gets this stop marker.
							lnext = 0xffffff0000000000ULL;

							dconf->ss_slices_p[slice].element[dconf->ss_slices_p[slice].size++] =
								(lnext | pid | sid);
							//(lnext | pid | sid | 0); // ((uint64_t)idx));

						//dconf->ss_slices_p[slice].element[dconf->ss_slices_p[slice].size++] =
						//	((uint64_t)i << 32) | ((uint64_t)(sum - prime) << 16) | ((uint64_t)idx);

							if (dconf->ss_slices_p[slice].size >= dconf->ss_slices_p[slice].alloc)
							{
								dconf->ss_slices_p[slice].alloc *= 2;
								dconf->ss_slices_p[slice].element = (uint64_t*)xrealloc(
									dconf->ss_slices_p[slice].element, dconf->ss_slices_p[slice].alloc *
									sizeof(uint64_t));
								//printf("ss bucket allocation now %d elements\n", dconf->ss_slices_p[slice].alloc);
							}

							//dconf->ss_slices_p[idx].root[dconf->ss_size_p[idx]] = sum - prime;
							//dconf->ss_slices_p[idx].fbid[dconf->ss_size_p[idx]] = i;
							//dconf->ss_slices_p[idx].logp[dconf->ss_size_p[idx]] = logp;
							//dconf->ss_size_p[idx]++;

							nummatch++;

#if defined( DEBUG_SS_SEARCH)
							int result = check_poly_at_loc(sum - prime, 0, prime, polysum, idx, polyv, polysign, dconf, sconf);
							if (result == 0)
							{
								printf("more info on failure:\n");
								printf("Bindex %d <- polymap[%04x]\n", idx, polysum);
								printf("bin1 root = %d, bin1 polynum = %08x\n", bins1[ii].root[j], bins1[ii].polynum[j]);
								printf("bin2 root = %d, bin2 polynum = %08x\n", bins2[b].root[k], bins2[b].polynum[k]);
								printf("bin1: %d, bin2: %d, center match pos side\n", ii, b);
								if (polysum == 0) exit(0);
						}
#endif

							//if (dconf->ss_size_p[idx] > dconf->ss_alloc_p[idx])
							//{
							//	dconf->ss_alloc_p[idx] *= 2;
							//	printf("reallocating ss_slices_p[%d] to %d\n",
							//		idx, dconf->ss_alloc_p[idx]);
							//	dconf->ss_slices_p[idx].root = (uint32_t*)xrealloc(
							//		dconf->ss_slices_p[idx].root, dconf->ss_alloc_p[idx] * 
							//		sizeof(uint32_t));
							//	dconf->ss_slices_p[idx].fbid = (uint32_t*)xrealloc(
							//		dconf->ss_slices_p[idx].fbid, dconf->ss_alloc_p[idx] *
							//		sizeof(uint32_t));
							//	dconf->ss_slices_p[idx].logp = (uint8_t*)xrealloc(
							//		dconf->ss_slices_p[idx].logp, dconf->ss_alloc_p[idx] *
							//		sizeof(uint8_t));
							//}

#if defined( DEBUG_SS_SEARCH) || defined(DEBUG_SS_SEARCH2)

							printf("+hit %d + %d = %d, %d, poly (%d,%d) %d <- %04x\n",
								bins1[ii].root[j], bins2[b].root[k], sum, sum - prime,
								bins1[ii].polynum[j], bins2[b].polynum[k],
								polymap[polysum], polysum);
							//int y;
							//for (y = 0; y < orighits.size; y++)
							//{
							//	if ((orighits.root[y] == sum - prime) && (orighits.polynum[y] == polysum))
							//	{
							//		//printf("*");
							//		nummatch++;
							//		break;
							//	}
							//}
							//printf("\n");
#endif
					}
						else if ((sum < prime) && ((0 - (sum - prime)) < interval))
						{
							int idx = polymap[polysum];

							if (idx == num_bpoly) continue;

							uint64_t pid = (uint64_t)(i - dconf->ss_slices_n[slice].fboffset) << 24;
							uint64_t sid = (uint64_t)(0 - (sum - prime));
							uint64_t lnext;

							update_ll(1, dconf->ss_slices_n, slice, idx, dconf);

							// the current entry gets this stop marker.
							lnext = 0xffffff0000000000ULL;

							dconf->ss_slices_n[slice].element[dconf->ss_slices_n[slice].size++] =
								(lnext | pid | sid);

							//dconf->ss_slices_n[slice].element[dconf->ss_slices_n[slice].size++] =
							//	((uint64_t)i << 32) | ((uint64_t)(0 - (sum - prime)) << 16) | ((uint64_t)idx);

							if (dconf->ss_slices_n[slice].size >= dconf->ss_slices_n[slice].alloc)
							{
								dconf->ss_slices_n[slice].alloc *= 2;
								dconf->ss_slices_n[slice].element = (uint64_t*)xrealloc(
									dconf->ss_slices_n[slice].element, dconf->ss_slices_n[slice].alloc *
									sizeof(uint64_t));
								//printf("ss bucket allocation now %d elements\n", dconf->ss_slices_p[slice].alloc);
							}

							//dconf->ss_slices_n[idx].root[dconf->ss_size_n[idx]] = 0 - (sum - prime);
							//dconf->ss_slices_n[idx].fbid[dconf->ss_size_n[idx]] = i;
							//dconf->ss_slices_n[idx].logp[dconf->ss_size_n[idx]] = logp;
							//dconf->ss_size_n[idx]++;

							nummatch++;

#if defined( DEBUG_SS_SEARCH)
							int result = check_poly_at_loc(0 - (sum - prime), 1, prime, polysum, idx, polyv, polysign, dconf, sconf);
							if (result == 0)
							{
								printf("more info on failure:\n");
								printf("Bindex %d <- polymap[%04x]\n", idx, polysum);
								printf("bin1 root = %d, bin1 polynum = %08x\n", bins1[ii].root[j], bins1[ii].polynum[j]);
								printf("bin2 root = %d, bin2 polynum = %08x\n", bins2[b].root[k], bins2[b].polynum[k]);
								printf("bin1: %d, bin2: %d, center match, neg side\n", ii, b);
								if (polysum == 0) exit(0);
						}
#endif

							//if (dconf->ss_size_n[idx] > dconf->ss_alloc_n[idx])
							//{
							//	dconf->ss_alloc_n[idx] *= 2;
							//	printf("reallocating ss_slices_n[%d] to %d\n",
							//		idx, dconf->ss_alloc_n[idx]);
							//	dconf->ss_slices_n[idx].root = (uint32_t*)xrealloc(
							//		dconf->ss_slices_n[idx].root, dconf->ss_alloc_n[idx] * 
							//		sizeof(uint32_t));
							//	dconf->ss_slices_n[idx].fbid = (uint32_t*)xrealloc(
							//		dconf->ss_slices_n[idx].fbid, dconf->ss_alloc_n[idx] *
							//		sizeof(uint32_t));
							//	dconf->ss_slices_n[idx].logp = (uint8_t*)xrealloc(
							//		dconf->ss_slices_n[idx].logp, dconf->ss_alloc_n[idx] *
							//		sizeof(uint8_t));
							//}

#if defined( DEBUG_SS_SEARCH) || defined(DEBUG_SS_SEARCH2)
							printf("-hit %d + %d = %d, %d, poly (%d,%d) %d <- %04x\n",
								bins1[ii].root[j], bins2[b].root[k], sum,
								0 - (sum - prime), bins1[ii].polynum[j], bins2[b].polynum[k],
								polymap[polysum], polysum);
							//int y;
							//for (y = 0; y < orighits.size; y++)
							//{
							//	if ((orighits.root[y] == (0 - (sum - prime))) && (orighits.polynum[y] == polysum))
							//	{
							//		//printf("*");
							//		nummatch++;
							//		break;
							//	}
							//}
							//printf("\n");
#endif
				}
			}
#if 1
					if ((b + 1) < numbins)
					{
						for (k = 0; k < bins2[b + 1].size; k++)
						{
							int sum = bins1[ii].root[j] + bins2[b + 1].root[k];
							int polysum = bins1[ii].polynum[j] + bins2[b + 1].polynum[k];
							if ((sum >= prime) && ((sum - prime) < interval)) {

								int idx = polymap[polysum];

								if (idx == num_bpoly) continue;

								uint64_t pid = (uint64_t)(i - dconf->ss_slices_p[slice].fboffset) << 24;
								uint64_t sid = (uint64_t)(sum - prime);
								uint64_t lnext;

								update_ll(0, dconf->ss_slices_p, slice, idx, dconf);

								// the current entry gets this stop marker.
								lnext = 0xffffff0000000000ULL;

								dconf->ss_slices_p[slice].element[dconf->ss_slices_p[slice].size++] =
									(lnext | pid | sid);

								//dconf->ss_slices_p[slice].element[dconf->ss_slices_p[slice].size++] =
								//	((uint64_t)i << 32) | ((uint64_t)(sum - prime) << 16) | ((uint64_t)idx);

								if (dconf->ss_slices_p[slice].size >= dconf->ss_slices_p[slice].alloc)
								{
									dconf->ss_slices_p[slice].alloc *= 2;
									dconf->ss_slices_p[slice].element = (uint64_t*)xrealloc(
										dconf->ss_slices_p[slice].element, dconf->ss_slices_p[slice].alloc *
										sizeof(uint64_t));
									//printf("ss bucket allocation now %d elements\n", dconf->ss_slices_p[slice].alloc);
								}

								//dconf->ss_slices_p[idx].root[dconf->ss_size_p[idx]] = sum - prime;
								//dconf->ss_slices_p[idx].fbid[dconf->ss_size_p[idx]] = i;
								//dconf->ss_slices_p[idx].logp[dconf->ss_size_p[idx]] = logp;
								//dconf->ss_size_p[idx]++;

								matchp1++;

#if defined( DEBUG_SS_SEARCH)
								int result = check_poly_at_loc(sum - prime, 0, prime, polysum, idx, polyv, polysign, dconf, sconf);
								if (result == 0)
								{
									printf("more info on failure:\n");
									printf("Bindex %d <- polymap[%04x]\n", idx, polysum);
									printf("bin1 root = %d, bin1 polynum = %04x\n", bins1[ii].root[j], bins1[ii].polynum[j]);
									printf("bin2 root = %d, bin2 polynum = %04x\n", bins2[b + 1].root[k], bins2[b + 1].polynum[k]);
									printf("bin1: %d, bin2: %d, right match pos side\n", ii, b + 1);
									if (polysum == 0) exit(0);
							}
#endif

								//if (dconf->ss_size_p[idx] > dconf->ss_alloc_p[idx])
								//{
								//	dconf->ss_alloc_p[idx] *= 2;
								//	printf("reallocating ss_slices_p[%d] to %d\n",
								//		idx, dconf->ss_alloc_p[idx]);
								//	dconf->ss_slices_p[idx].root = (uint32_t*)xrealloc(
								//		dconf->ss_slices_p[idx].root, dconf->ss_alloc_p[idx] *
								//		sizeof(uint32_t));
								//	dconf->ss_slices_p[idx].fbid = (uint32_t*)xrealloc(
								//		dconf->ss_slices_p[idx].fbid, dconf->ss_alloc_p[idx] *
								//		sizeof(uint32_t));
								//	dconf->ss_slices_p[idx].logp = (uint8_t*)xrealloc(
								//		dconf->ss_slices_p[idx].logp, dconf->ss_alloc_p[idx] *
								//		sizeof(uint8_t));
								//}

#if defined( DEBUG_SS_SEARCH) || defined(DEBUG_SS_SEARCH2)
								printf("+hit %d + %d = %d, %d, poly (%d,%d) %04x",
									bins1[ii].root[j], bins2[b + 1].root[k], sum, sum - prime, bins1[ii].polynum[j], bins2[b + 1].polynum[k], polysum);
								int y;
								for (y = 0; y < orighits.size; y++)
								{
									if ((orighits.root[y] == sum - prime) && (orighits.polynum[y] == polysum))
									{
										printf("*");
										nummatch++;
										matchp1++;
										break;
						}
					}
								printf("\n");
#endif
		}
							else if ((0 - (sum - prime)) < interval)
							{

								int idx = polymap[polysum];

								if (idx == num_bpoly) continue;

								uint64_t pid = (uint64_t)(i - dconf->ss_slices_n[slice].fboffset) << 24;
								uint64_t sid = (uint64_t)(0 - (sum - prime));
								uint64_t lnext;

								update_ll(1, dconf->ss_slices_n, slice, idx, dconf);

								// the current entry gets this stop marker.
								lnext = 0xffffff0000000000ULL;

								dconf->ss_slices_n[slice].element[dconf->ss_slices_n[slice].size++] =
									(lnext | pid | sid);

								//dconf->ss_slices_n[slice].element[dconf->ss_slices_n[slice].size++] =
								//	((uint64_t)i << 32) | ((uint64_t)(0 - (sum - prime)) << 16) | ((uint64_t)idx);

								if (dconf->ss_slices_n[slice].size >= dconf->ss_slices_n[slice].alloc)
								{
									dconf->ss_slices_n[slice].alloc *= 2;
									dconf->ss_slices_n[slice].element = (uint64_t*)xrealloc(
										dconf->ss_slices_n[slice].element, dconf->ss_slices_n[slice].alloc *
										sizeof(uint64_t));
									//printf("ss bucket allocation now %d elements\n", dconf->ss_slices_p[slice].alloc);
								}

								//dconf->ss_slices_n[idx].root[dconf->ss_size_n[idx]] = 0 - (sum - prime);
								//dconf->ss_slices_n[idx].fbid[dconf->ss_size_n[idx]] = i;
								//dconf->ss_slices_n[idx].logp[dconf->ss_size_n[idx]] = logp;
								//dconf->ss_size_n[idx]++;

								matchp1++;

#if defined( DEBUG_SS_SEARCH)
								int result = check_poly_at_loc(0 - (sum - prime), 1, prime, polysum, idx, polyv, polysign, dconf, sconf);
								if (result == 0)
								{
									printf("more info on failure:\n");
									printf("Bindex %d <- polymap[%04x]\n", idx, polysum);
									printf("bin1 root = %d, bin1 polynum = %04x\n", bins1[ii].root[j], bins1[ii].polynum[j]);
									printf("bin2 root = %d, bin2 polynum = %04x\n", bins2[b + 1].root[k], bins2[b + 1].polynum[k]);
									printf("bin1: %d, bin2: %d, right match, neg side\n", ii, b + 1);
									if (polysum == 0) exit(0);
							}
#endif

								//if (dconf->ss_size_n[idx] > dconf->ss_alloc_n[idx])
								//{
								//	dconf->ss_alloc_n[idx] *= 2;
								//	printf("reallocating ss_slices_n[%d] to %d\n",
								//		idx, dconf->ss_alloc_n[idx]);
								//	dconf->ss_slices_n[idx].root = (uint32_t*)xrealloc(
								//		dconf->ss_slices_n[idx].root, dconf->ss_alloc_n[idx] *
								//		sizeof(uint32_t));
								//	dconf->ss_slices_n[idx].fbid = (uint32_t*)xrealloc(
								//		dconf->ss_slices_n[idx].fbid, dconf->ss_alloc_n[idx] *
								//		sizeof(uint32_t));
								//	dconf->ss_slices_n[idx].logp = (uint8_t*)xrealloc(
								//		dconf->ss_slices_n[idx].logp, dconf->ss_alloc_n[idx] *
								//		sizeof(uint8_t));
								//}

#if defined( DEBUG_SS_SEARCH) || defined(DEBUG_SS_SEARCH2)
								printf("-hit %d + %d = %d, %d, poly (%d,%d) %04x",
									bins1[ii].root[j], bins2[b + 1].root[k], sum, 0 - (sum - prime), bins1[ii].polynum[j], bins2[b + 1].polynum[k], polysum);
								int y;
								for (y = 0; y < orighits.size; y++)
								{
									if ((orighits.root[y] == (0 - (sum - prime))) && (orighits.polynum[y] == polysum))
									{
										printf("*");
										nummatch++;
										matchp1++;
										break;
			}
		}
								printf("\n");
#endif
			}
		}
						}
					if ((b - 1) >= 0)
					{
						for (k = 0; k < bins2[b - 1].size; k++)
						{
							int sum = bins1[ii].root[j] + bins2[b - 1].root[k];
							int polysum = bins1[ii].polynum[j] + bins2[b - 1].polynum[k];
							if ((sum >= prime) && ((sum - prime) < interval)) {

								int idx = polymap[polysum];

								if (idx == num_bpoly) continue;

								uint64_t pid = (uint64_t)(i - dconf->ss_slices_p[slice].fboffset) << 24;
								uint64_t sid = (uint64_t)(sum - prime);
								uint64_t lnext;

								update_ll(0, dconf->ss_slices_p, slice, idx, dconf);

								// the current entry gets this stop marker.
								lnext = 0xffffff0000000000ULL;

								dconf->ss_slices_p[slice].element[dconf->ss_slices_p[slice].size++] =
									(lnext | pid | sid);

								//dconf->ss_slices_p[slice].element[dconf->ss_slices_p[slice].size++] =
								//	((uint64_t)i << 32) | ((uint64_t)(sum - prime) << 16) | ((uint64_t)idx);

								if (dconf->ss_slices_p[slice].size >= dconf->ss_slices_p[slice].alloc)
								{
									dconf->ss_slices_p[slice].alloc *= 2;
									dconf->ss_slices_p[slice].element = (uint64_t*)xrealloc(
										dconf->ss_slices_p[slice].element, dconf->ss_slices_p[slice].alloc *
										sizeof(uint64_t));
									//printf("ss bucket allocation now %d elements\n", dconf->ss_slices_p[slice].alloc);
								}

								//dconf->ss_slices_p[idx].root[dconf->ss_size_p[idx]] = sum - prime;
								//dconf->ss_slices_p[idx].fbid[dconf->ss_size_p[idx]] = i;
								//dconf->ss_slices_p[idx].logp[dconf->ss_size_p[idx]] = logp;
								//dconf->ss_size_p[idx]++;

								matchm1++;

#if defined( DEBUG_SS_SEARCH)
								int result = check_poly_at_loc(sum - prime, 0, prime, polysum, idx, polyv, polysign, dconf, sconf);
								if (result == 0)
								{
									printf("more info on failure:\n");
									printf("Bindex %d <- polymap[%04x]\n", idx, polysum);
									printf("bin1 root = %d, bin1 polynum = %04x\n", bins1[ii].root[j], bins1[ii].polynum[j]);
									printf("bin2 root = %d, bin2 polynum = %04x\n", bins2[b - 1].root[k], bins2[b - 1].polynum[k]);
									printf("bin1: %d, bin2: %d, left match pos side\n", ii, b - 1);
									if (polysum == 0) exit(0);
							}
#endif

								//if (dconf->ss_size_p[idx] > dconf->ss_alloc_p[idx])
								//{
								//	dconf->ss_alloc_p[idx] *= 2;
								//	printf("reallocating ss_slices_p[%d] to %d\n",
								//		idx, dconf->ss_alloc_p[idx]);
								//	dconf->ss_slices_p[idx].root = (uint32_t*)xrealloc(
								//		dconf->ss_slices_p[idx].root, dconf->ss_alloc_p[idx] *
								//		sizeof(uint32_t));
								//	dconf->ss_slices_p[idx].fbid = (uint32_t*)xrealloc(
								//		dconf->ss_slices_p[idx].fbid, dconf->ss_alloc_p[idx] *
								//		sizeof(uint32_t));
								//	dconf->ss_slices_p[idx].logp = (uint8_t*)xrealloc(
								//		dconf->ss_slices_p[idx].logp, dconf->ss_alloc_p[idx] *
								//		sizeof(uint8_t));
								//}

#if defined( DEBUG_SS_SEARCH) || defined(DEBUG_SS_SEARCH2)
								printf("+hit %d + %d = %d, %d, poly (%d,%d) %04x",
									bins1[ii].root[j], bins2[b - 1].root[k], sum, sum - prime, bins1[ii].polynum[j], bins2[b - 1].polynum[k], polysum);
								int y;
								for (y = 0; y < orighits.size; y++)
								{
									if ((orighits.root[y] == sum - prime) && (orighits.polynum[y] == polysum))
									{
										printf("*");
										nummatch++;
										matchm1++;
										break;
						}
					}
								printf("\n");
#endif
								}
							else if ((0 - (sum - prime)) < interval)
							{

								int idx = polymap[polysum];

								if (idx == num_bpoly) continue;

								uint64_t pid = (uint64_t)(i - dconf->ss_slices_n[slice].fboffset) << 24;
								uint64_t sid = (uint64_t)(0 - (sum - prime));
								uint64_t lnext;

								update_ll(1, dconf->ss_slices_n, slice, idx, dconf);

								// the current entry gets this stop marker.
								lnext = 0xffffff0000000000ULL;

								dconf->ss_slices_n[slice].element[dconf->ss_slices_n[slice].size++] =
									(lnext | pid | sid);

								//dconf->ss_slices_n[slice].element[dconf->ss_slices_n[slice].size++] =
								//	((uint64_t)i << 32) | ((uint64_t)(0 - (sum - prime)) << 16) | ((uint64_t)idx);

								if (dconf->ss_slices_n[slice].size >= dconf->ss_slices_n[slice].alloc)
								{
									dconf->ss_slices_n[slice].alloc *= 2;
									dconf->ss_slices_n[slice].element = (uint64_t*)xrealloc(
										dconf->ss_slices_n[slice].element, dconf->ss_slices_n[slice].alloc *
										sizeof(uint64_t));
									//printf("ss bucket allocation now %d elements\n", dconf->ss_slices_p[slice].alloc);
								}

								//dconf->ss_slices_n[idx].root[dconf->ss_size_n[idx]] = 0 - (sum - prime);
								//dconf->ss_slices_n[idx].fbid[dconf->ss_size_n[idx]] = i;
								//dconf->ss_slices_n[idx].logp[dconf->ss_size_n[idx]] = logp;
								//dconf->ss_size_n[idx]++;

								matchm1++;

#if defined( DEBUG_SS_SEARCH)
								int result = check_poly_at_loc(0 - (sum - prime), 1, prime, polysum, idx, polyv, polysign, dconf, sconf);
								if (result == 0)
								{
									printf("more info on failure:\n");
									printf("Bindex %d <- polymap[%04x]\n", idx, polysum);
									printf("bin1 root = %d, bin1 polynum = %04x\n", bins1[ii].root[j], bins1[ii].polynum[j]);
									printf("bin2 root = %d, bin2 polynum = %04x\n", bins2[b - 1].root[k], bins2[b - 1].polynum[k]);
									printf("bin1: %d, bin2: %d, left match, neg side\n", ii, b - 1);
									if (polysum == 0) exit(0);
							}
#endif

								//if (dconf->ss_size_n[idx] > dconf->ss_alloc_n[idx])
								//{
								//	dconf->ss_alloc_n[idx] *= 2;
								//	printf("reallocating ss_slices_n[%d] to %d\n",
								//		idx, dconf->ss_alloc_n[idx]);
								//	dconf->ss_slices_n[idx].root = (uint32_t*)xrealloc(
								//		dconf->ss_slices_n[idx].root, dconf->ss_alloc_n[idx] *
								//		sizeof(uint32_t));
								//	dconf->ss_slices_n[idx].fbid = (uint32_t*)xrealloc(
								//		dconf->ss_slices_n[idx].fbid, dconf->ss_alloc_n[idx] *
								//		sizeof(uint32_t));
								//	dconf->ss_slices_n[idx].logp = (uint8_t*)xrealloc(
								//		dconf->ss_slices_n[idx].logp, dconf->ss_alloc_n[idx] *
								//		sizeof(uint8_t));
								//}


#if defined( DEBUG_SS_SEARCH) || defined(DEBUG_SS_SEARCH2)
								printf("-hit %d + %d = %d, %d, poly (%d,%d) %04x",
									bins1[ii].root[j], bins2[b - 1].root[k], sum, 0 - (sum - prime), bins1[ii].polynum[j], bins2[b - 1].polynum[k], polysum);
								int y;
								for (y = 0; y < orighits.size; y++)
								{
									if ((orighits.root[y] == (0 - (sum - prime))) && (orighits.polynum[y] == polysum))
									{
										printf("*");
										nummatch++;
										matchm1++;
										break;
										}
									}
								printf("\n");
#endif
								}
							}
						}
#endif
					}
				}

			//printf("found %d sieve hits matching bin %d to bin %d\n", nummatch, ii, b);
			//printf("found %d sieve hits matching bin %d to bin %d\n", matchp1, ii, b + 1);
			//printf("found %d sieve hits matching bin %d to bin %d\n", matchm1, ii, b - 1);

			}

		// reset set1 so we can put root2 into it
		for (ii = 0; ii < numbins; ii++)
		{
			bins1[ii].size = 0;
		}

#ifdef DEBUG_SS_SEARCH
		printf("sorting rootset1b of size %d into %d bins of size %d\n", (1 << ss_set1b.size), numbins, binsize);
#endif

		// sort root2 into set1
		for (ii = 0; ii < (1 << ss_set1b.size); ii++)
		{
			int binnum = ss_set1b.root[ii] / binsize;
			if (binnum < numbins)
			{
				// printf("bin %d <-- set1b root %d (polynum 0)\n", binnum, ss_set1a.root[ii]);
				bins1[binnum].root[bins1[binnum].size] = ss_set1b.root[ii];
				bins1[binnum].polynum[bins1[binnum].size] = ss_set1b.polynum[ii];
				bins1[binnum].size++;
				if (bins1[binnum].size > bindepth)
				{
					printf("bin overflow\n");
					exit(1);
				}
			}
			else
			{
				printf("element %d of set 1, root %d, invalid bin %d [of %d]\n",
					ii, ss_set1b.root[ii], binnum, numbins);
			}
		}


#ifdef DEBUG_SS_SEARCH
		printf("sorting rootset2b of size %d into %d bins of size %d\n", (1 << ss_set2b.size), numbins, binsize);
#endif

		// now match all numbers in bin x with the other set's bin that
		// contains p - x
#ifdef DEBUG_SS_SEARCH
		printf("matching root2\n");
#endif

		// match set1 holding root2 to set2
		for (ii = 0; ii < numbins; ii++)
		{
			int x = binsize * ii + binsize / 2;
			int px = p - x;
			int b = px / binsize; // px / binsize;

			nummatch = 0;
			matchp1 = 0;
			matchm1 = 0;

			//printf("bin %d (cent %d) with %d elements matching with bin %d with %d elements\n",
			//	ii, x, bins1[ii].size, b, bins2[b].size);
			//
			//if ((b + 1) < numbins)
			//{
			//	printf("bin %d (cent %d) with %d elements matching with bin %d with %d elements\n",
			//		ii, x, bins1[ii].size, b + 1, bins2[b + 1].size);
			//}
			//
			//if ((b - 1) >= 0)
			//{
			//	printf("bin %d (cent %d) with %d elements matching with bin %d with %d elements\n",
			//		ii, x, bins1[ii].size, b - 1, bins2[b - 1].size);
			//}
			//
			//for (b = minb; b < maxb; b++)

			// todo:
			// make root and fbid both 16 bit values in the ss_slice_t type.
			// remove logp.
			// create array of size values for each bucket that track when to
			// increment fb offset of fbids by 65536.  also record new logp at that point.
			// similar concept to fb slices in bucket sort, applies to each ss_bucket.
			// vector compare one entry of bin1 to 16 entries of bin2.  
			// scatter into ss_buckets.
			{
				for (j = 0; j < bins1[ii].size; j++)
				{
					int k;
					for (k = 0; k < bins2[b].size; k++)
					{
						int sum = bins1[ii].root[j] + bins2[b].root[k];
						int polysum = bins1[ii].polynum[j] + bins2[b].polynum[k];

						if (polysum >= (1 << (poly->s - 1)))
						{
							printf("invalid polysum %d (%d + %d)\n",
								polysum, bins1[ii].polynum[j], bins2[b].polynum[k]);
							continue;
						}

						if ((sum >= prime) && ((sum - prime) < interval))
						{
							int idx = polymap[polysum];

							if (idx == num_bpoly) continue;

							uint64_t pid = (uint64_t)(i - dconf->ss_slices_p[slice].fboffset) << 24;
							uint64_t sid = (uint64_t)(sum - prime);
							uint64_t lnext;

							update_ll(0, dconf->ss_slices_p, slice, idx, dconf);

							// the current entry gets this stop marker.
							lnext = 0xffffff0000000000ULL;

							dconf->ss_slices_p[slice].element[dconf->ss_slices_p[slice].size++] =
								(lnext | pid | sid);

							//dconf->ss_slices_p[slice].element[dconf->ss_slices_p[slice].size++] =
							//	((uint64_t)i << 32) | ((uint64_t)(sum - prime) << 16) | ((uint64_t)idx);

							if (dconf->ss_slices_p[slice].size >= dconf->ss_slices_p[slice].alloc)
							{
								dconf->ss_slices_p[slice].alloc *= 2;
								dconf->ss_slices_p[slice].element = (uint64_t*)xrealloc(
									dconf->ss_slices_p[slice].element, dconf->ss_slices_p[slice].alloc *
									sizeof(uint64_t));
								//printf("ss bucket allocation now %d elements\n", dconf->ss_slices_p[slice].alloc);
							}

							//dconf->ss_slices_p[idx].root[dconf->ss_size_p[idx]] = sum - prime;
							//dconf->ss_slices_p[idx].fbid[dconf->ss_size_p[idx]] = i;
							//dconf->ss_slices_p[idx].logp[dconf->ss_size_p[idx]] = logp;
							//dconf->ss_size_p[idx]++;

							nummatch++;

#if defined( DEBUG_SS_SEARCH)
							int result = check_poly_at_loc(sum - prime, 0, prime, polysum, idx, polyv, polysign, dconf, sconf);
							if (result == 0)
							{
								printf("more info on failure:\n");
								printf("Bindex %d <- polymap[%04x]\n", idx, polysum);
								printf("bin1 root = %d, bin1 polynum = %08x\n", bins1[ii].root[j], bins1[ii].polynum[j]);
								printf("bin2 root = %d, bin2 polynum = %08x\n", bins2[b].root[k], bins2[b].polynum[k]);
								printf("bin1: %d, bin2: %d, center match pos side\n", ii, b);
								if (polysum == 0) exit(0);
						}
#endif



							//if (dconf->ss_size_p[idx] > dconf->ss_alloc_p[idx])
							//{
							//	dconf->ss_alloc_p[idx] *= 2;
							//	printf("reallocating ss_slices_p[%d] to %d\n",
							//		idx, dconf->ss_alloc_p[idx]);
							//	dconf->ss_slices_p[idx].root = (uint32_t*)xrealloc(
							//		dconf->ss_slices_p[idx].root, dconf->ss_alloc_p[idx] *
							//		sizeof(uint32_t));
							//	dconf->ss_slices_p[idx].fbid = (uint32_t*)xrealloc(
							//		dconf->ss_slices_p[idx].fbid, dconf->ss_alloc_p[idx] *
							//		sizeof(uint32_t));
							//	dconf->ss_slices_p[idx].logp = (uint8_t*)xrealloc(
							//		dconf->ss_slices_p[idx].logp, dconf->ss_alloc_p[idx] *
							//		sizeof(uint8_t));
							//}

#if defined( DEBUG_SS_SEARCH) || defined(DEBUG_SS_SEARCH2)
							if (prime == 1389097)
								printf("+hit %d + %d = %d, %d, poly (%d,%d) %d <- %04x\n",
									bins1[ii].root[j], bins2[b].root[k], sum, sum - prime,
									bins1[ii].polynum[j], bins2[b].polynum[k],
									polymap[polysum], polysum);//int y;
							//for (y = 0; y < orighits.size; y++)
							//{
							//	if ((orighits.root[y] == sum - prime) && (orighits.polynum[y] == polysum))
							//	{
							//		//printf("*");
							//		nummatch++;
							//		break;
							//	}
							//}
							//printf("\n");
#endif
					}
						else if ((sum < prime) && ((0 - (sum - prime)) < interval))
						{
							int idx = polymap[polysum];

							if (idx == num_bpoly) continue;

							uint64_t pid = (uint64_t)(i - dconf->ss_slices_n[slice].fboffset) << 24;
							uint64_t sid = (uint64_t)(0 - (sum - prime));
							uint64_t lnext;

							update_ll(1, dconf->ss_slices_n, slice, idx, dconf);

							// the current entry gets this stop marker.
							lnext = 0xffffff0000000000ULL;

							dconf->ss_slices_n[slice].element[dconf->ss_slices_n[slice].size++] =
								(lnext | pid | sid);

							//dconf->ss_slices_n[slice].element[dconf->ss_slices_n[slice].size++] =
							//	((uint64_t)i << 32) | ((uint64_t)(0 - (sum - prime)) << 16) | ((uint64_t)idx);

							if (dconf->ss_slices_n[slice].size >= dconf->ss_slices_n[slice].alloc)
							{
								dconf->ss_slices_n[slice].alloc *= 2;
								dconf->ss_slices_n[slice].element = (uint64_t*)xrealloc(
									dconf->ss_slices_n[slice].element, dconf->ss_slices_n[slice].alloc *
									sizeof(uint64_t));
								//printf("ss bucket allocation now %d elements\n", dconf->ss_slices_p[slice].alloc);
							}

							//dconf->ss_slices_n[idx].root[dconf->ss_size_n[idx]] = 0 - (sum - prime);
							//dconf->ss_slices_n[idx].fbid[dconf->ss_size_n[idx]] = i;
							//dconf->ss_slices_n[idx].logp[dconf->ss_size_n[idx]] = logp;
							//dconf->ss_size_n[idx]++;

							nummatch++;

#if defined( DEBUG_SS_SEARCH)
							int result = check_poly_at_loc(0 - (sum - prime), 1, prime, polysum, idx, polyv, polysign, dconf, sconf);
							if (result == 0)
							{
								printf("more info on failure:\n");
								printf("Bindex %d <- polymap[%04x]\n", idx, polysum);
								printf("bin1 root = %d, bin1 polynum = %04x\n", bins1[ii].root[j], bins1[ii].polynum[j]);
								printf("bin2 root = %d, bin2 polynum = %04x\n", bins2[b].root[k], bins2[b].polynum[k]);
								printf("bin1: %d, bin2: %d, left match, neg side\n", ii, b);
								if (polysum == 0) exit(0);
						}
#endif

							//if (dconf->ss_size_n[idx] > dconf->ss_alloc_n[idx])
							//{
							//	dconf->ss_alloc_n[idx] *= 2;
							//	printf("reallocating ss_slices_n[%d] to %d\n",
							//		idx, dconf->ss_alloc_n[idx]);
							//	dconf->ss_slices_n[idx].root = (uint32_t*)xrealloc(
							//		dconf->ss_slices_n[idx].root, dconf->ss_alloc_n[idx] *
							//		sizeof(uint32_t));
							//	dconf->ss_slices_n[idx].fbid = (uint32_t*)xrealloc(
							//		dconf->ss_slices_n[idx].fbid, dconf->ss_alloc_n[idx] *
							//		sizeof(uint32_t));
							//	dconf->ss_slices_n[idx].logp = (uint8_t*)xrealloc(
							//		dconf->ss_slices_n[idx].logp, dconf->ss_alloc_n[idx] *
							//		sizeof(uint8_t));
							//}

#if defined( DEBUG_SS_SEARCH) || defined(DEBUG_SS_SEARCH2)
							if (prime == 1389097)
								printf("-hit %d + %d = %d, %d, poly (%d,%d) %d <- %04x\n",
									bins1[ii].root[j], bins2[b].root[k], sum,
									0 - (sum - prime), bins1[ii].polynum[j], bins2[b].polynum[k],
									polymap[polysum], polysum);
							//int y;
							//for (y = 0; y < orighits.size; y++)
							//{
							//	if ((orighits.root[y] == (0 - (sum - prime))) && (orighits.polynum[y] == polysum))
							//	{
							//		//printf("*");
							//		nummatch++;
							//		break;
							//	}
							//}
							//printf("\n");
#endif
				}
			}
#if 1
					if ((b + 1) < numbins)
					{
						for (k = 0; k < bins2[b + 1].size; k++)
						{
							int sum = bins1[ii].root[j] + bins2[b + 1].root[k];
							int polysum = bins1[ii].polynum[j] + bins2[b + 1].polynum[k];
							if ((sum >= prime) && ((sum - prime) < interval)) {

								int idx = polymap[polysum];

								if (idx == num_bpoly) continue;

								uint64_t pid = (uint64_t)(i - dconf->ss_slices_p[slice].fboffset) << 24;
								uint64_t sid = (uint64_t)(sum - prime);
								uint64_t lnext;

								update_ll(0, dconf->ss_slices_p, slice, idx, dconf);

								// the current entry gets this stop marker.
								lnext = 0xffffff0000000000ULL;

								dconf->ss_slices_p[slice].element[dconf->ss_slices_p[slice].size++] =
									(lnext | pid | sid);

								//dconf->ss_slices_p[slice].element[dconf->ss_slices_p[slice].size++] =
								//	((uint64_t)i << 32) | ((uint64_t)(sum - prime) << 16) | ((uint64_t)idx);

								if (dconf->ss_slices_p[slice].size >= dconf->ss_slices_p[slice].alloc)
								{
									dconf->ss_slices_p[slice].alloc *= 2;
									dconf->ss_slices_p[slice].element = (uint64_t*)xrealloc(
										dconf->ss_slices_p[slice].element, dconf->ss_slices_p[slice].alloc *
										sizeof(uint64_t));
									//printf("ss bucket allocation now %d elements\n", dconf->ss_slices_p[slice].alloc);
								}

								//dconf->ss_slices_p[idx].root[dconf->ss_size_p[idx]] = sum - prime;
								//dconf->ss_slices_p[idx].fbid[dconf->ss_size_p[idx]] = i;
								//dconf->ss_slices_p[idx].logp[dconf->ss_size_p[idx]] = logp;
								//dconf->ss_size_p[idx]++;

								matchp1++;

#if defined( DEBUG_SS_SEARCH)
								int result = check_poly_at_loc(sum - prime, 0, prime, polysum, idx, polyv, polysign, dconf, sconf);
								if (result == 0)
								{
									printf("more info on failure:\n");
									printf("Bindex %d <- polymap[%04x]\n", idx, polysum);
									printf("bin1 root = %d, bin1 polynum = %08x\n", bins1[ii].root[j], bins1[ii].polynum[j]);
									printf("bin2 root = %d, bin2 polynum = %08x\n", bins2[b + 1].root[k], bins2[b + 1].polynum[k]);
									printf("bin1: %d, bin2: %d, right match pos side\n", ii, b + 1);
									if (polysum == 0) exit(0);
							}
#endif

								//if (dconf->ss_size_p[idx] > dconf->ss_alloc_p[idx])
								//{
								//	dconf->ss_alloc_p[idx] *= 2;
								//	printf("reallocating ss_slices_p[%d] to %d\n",
								//		idx, dconf->ss_alloc_p[idx]);
								//	dconf->ss_slices_p[idx].root = (uint32_t*)xrealloc(
								//		dconf->ss_slices_p[idx].root, dconf->ss_alloc_p[idx] *
								//		sizeof(uint32_t));
								//	dconf->ss_slices_p[idx].fbid = (uint32_t*)xrealloc(
								//		dconf->ss_slices_p[idx].fbid, dconf->ss_alloc_p[idx] *
								//		sizeof(uint32_t));
								//	dconf->ss_slices_p[idx].logp = (uint8_t*)xrealloc(
								//		dconf->ss_slices_p[idx].logp, dconf->ss_alloc_p[idx] *
								//		sizeof(uint8_t));
								//}


#if defined( DEBUG_SS_SEARCH) || defined(DEBUG_SS_SEARCH2)
								printf("+hit %d + %d = %d, %d, poly (%d,%d) %04x",
									bins1[ii].root[j], bins2[b + 1].root[k], sum, sum - prime, bins1[ii].polynum[j], bins2[b + 1].polynum[k], polysum);
								int y;
								for (y = 0; y < orighits.size; y++)
								{
									if ((orighits.root[y] == sum - prime) && (orighits.polynum[y] == polysum))
									{
										printf("*");
										nummatch++;
										matchp1++;
										break;
						}
					}
								printf("\n");
#endif
		}
							else if ((0 - (sum - prime)) < interval)
							{

								int idx = polymap[polysum];

								if (idx == num_bpoly) continue;

								uint64_t pid = (uint64_t)(i - dconf->ss_slices_n[slice].fboffset) << 24;
								uint64_t sid = (uint64_t)(0 - (sum - prime));
								uint64_t lnext;

								update_ll(1, dconf->ss_slices_n, slice, idx, dconf);

								// the current entry gets this stop marker.
								lnext = 0xffffff0000000000ULL;

								dconf->ss_slices_n[slice].element[dconf->ss_slices_n[slice].size++] =
									(lnext | pid | sid);

								//dconf->ss_slices_n[slice].element[dconf->ss_slices_n[slice].size++] =
								//	((uint64_t)i << 32) | ((uint64_t)(0 - (sum - prime)) << 16) | ((uint64_t)idx);

								if (dconf->ss_slices_n[slice].size >= dconf->ss_slices_n[slice].alloc)
								{
									dconf->ss_slices_n[slice].alloc *= 2;
									dconf->ss_slices_n[slice].element = (uint64_t*)xrealloc(
										dconf->ss_slices_n[slice].element, dconf->ss_slices_n[slice].alloc *
										sizeof(uint64_t));
									//printf("ss bucket allocation now %d elements\n", dconf->ss_slices_p[slice].alloc);
								}

								//dconf->ss_slices_n[idx].root[dconf->ss_size_n[idx]] = 0 - (sum - prime);
								//dconf->ss_slices_n[idx].fbid[dconf->ss_size_n[idx]] = i;
								//dconf->ss_slices_n[idx].logp[dconf->ss_size_n[idx]] = logp;
								//dconf->ss_size_n[idx]++;

								matchp1++;

#if defined( DEBUG_SS_SEARCH)
								int result = check_poly_at_loc(0 - (sum - prime), 1, prime, polysum, idx, polyv, polysign, dconf, sconf);
								if (result == 0)
								{
									printf("more info on failure:\n");
									printf("Bindex %d <- polymap[%04x]\n", idx, polysum);
									printf("bin1 root = %d, bin1 polynum = %04x\n", bins1[ii].root[j], bins1[ii].polynum[j]);
									printf("bin2 root = %d, bin2 polynum = %04x\n", bins2[b + 1].root[k], bins2[b + 1].polynum[k]);
									printf("bin1: %d, bin2: %d, left match, neg side\n", ii, b + 1);
									if (polysum == 0) exit(0);
							}
#endif

								//if (dconf->ss_size_n[idx] > dconf->ss_alloc_n[idx])
								//{
								//	dconf->ss_alloc_n[idx] *= 2;
								//	printf("reallocating ss_slices_n[%d] to %d\n",
								//		idx, dconf->ss_alloc_n[idx]);
								//	dconf->ss_slices_n[idx].root = (uint32_t*)xrealloc(
								//		dconf->ss_slices_n[idx].root, dconf->ss_alloc_n[idx] *
								//		sizeof(uint32_t));
								//	dconf->ss_slices_n[idx].fbid = (uint32_t*)xrealloc(
								//		dconf->ss_slices_n[idx].fbid, dconf->ss_alloc_n[idx] *
								//		sizeof(uint32_t));
								//	dconf->ss_slices_n[idx].logp = (uint8_t*)xrealloc(
								//		dconf->ss_slices_n[idx].logp, dconf->ss_alloc_n[idx] *
								//		sizeof(uint8_t));
								//}

#if defined( DEBUG_SS_SEARCH) || defined(DEBUG_SS_SEARCH2)
								printf("-hit %d + %d = %d, %d, poly (%d,%d) %04x",
									bins1[ii].root[j], bins2[b + 1].root[k], sum, 0 - (sum - prime), bins1[ii].polynum[j], bins2[b + 1].polynum[k], polysum);
								int y;
								for (y = 0; y < orighits.size; y++)
								{
									if ((orighits.root[y] == (0 - (sum - prime))) && (orighits.polynum[y] == polysum))
									{
										printf("*");
										nummatch++;
										matchp1++;
										break;
		}
		}
								printf("\n");
#endif
								}
							}
						}
					if ((b - 1) >= 0)
					{
						for (k = 0; k < bins2[b - 1].size; k++)
						{
							int sum = bins1[ii].root[j] + bins2[b - 1].root[k];
							int polysum = bins1[ii].polynum[j] + bins2[b - 1].polynum[k];
							if ((sum >= prime) && ((sum - prime) < interval)) {

								int idx = polymap[polysum];

								if (idx == num_bpoly) continue;

								uint64_t pid = (uint64_t)(i - dconf->ss_slices_p[slice].fboffset) << 24;
								uint64_t sid = (uint64_t)(sum - prime);
								uint64_t lnext;

								update_ll(0, dconf->ss_slices_p, slice, idx, dconf);

								// the current entry gets this stop marker.
								lnext = 0xffffff0000000000ULL;

								dconf->ss_slices_p[slice].element[dconf->ss_slices_p[slice].size++] =
									(lnext | pid | sid);

								//dconf->ss_slices_p[slice].element[dconf->ss_slices_p[slice].size++] =
								//	((uint64_t)i << 32) | ((uint64_t)(sum - prime) << 16) | ((uint64_t)idx);

								if (dconf->ss_slices_p[slice].size >= dconf->ss_slices_p[slice].alloc)
								{
									dconf->ss_slices_p[slice].alloc *= 2;
									dconf->ss_slices_p[slice].element = (uint64_t*)xrealloc(
										dconf->ss_slices_p[slice].element, dconf->ss_slices_p[slice].alloc *
										sizeof(uint64_t));
									//printf("ss bucket allocation now %d elements\n", dconf->ss_slices_p[slice].alloc);
								}

								//dconf->ss_slices_p[idx].root[dconf->ss_size_p[idx]] = sum - prime;
								//dconf->ss_slices_p[idx].fbid[dconf->ss_size_p[idx]] = i;
								//dconf->ss_slices_p[idx].logp[dconf->ss_size_p[idx]] = logp;
								//dconf->ss_size_p[idx]++;

								matchm1++;

#if defined( DEBUG_SS_SEARCH)
								int result = check_poly_at_loc(sum - prime, 0, prime, polysum, idx, polyv, polysign, dconf, sconf);
								if (result == 0)
								{
									printf("more info on failure:\n");
									printf("Bindex %d <- polymap[%04x]\n", idx, polysum);
									printf("bin1 root = %d, bin1 polynum = %08x\n", bins1[ii].root[j], bins1[ii].polynum[j]);
									printf("bin2 root = %d, bin2 polynum = %08x\n", bins2[b - 1].root[k], bins2[b - 1].polynum[k]);
									printf("bin1: %d, bin2: %d, center match pos side\n", ii, b - 1);
									if (polysum == 0) exit(0);
							}
#endif

								//if (dconf->ss_size_p[idx] > dconf->ss_alloc_p[idx])
								//{
								//	dconf->ss_alloc_p[idx] *= 2;
								//	printf("reallocating ss_slices_p[%d] to %d\n",
								//		idx, dconf->ss_alloc_p[idx]);
								//	dconf->ss_slices_p[idx].root = (uint32_t*)xrealloc(
								//		dconf->ss_slices_p[idx].root, dconf->ss_alloc_p[idx] *
								//		sizeof(uint32_t));
								//	dconf->ss_slices_p[idx].fbid = (uint32_t*)xrealloc(
								//		dconf->ss_slices_p[idx].fbid, dconf->ss_alloc_p[idx] *
								//		sizeof(uint32_t));
								//	dconf->ss_slices_p[idx].logp = (uint8_t*)xrealloc(
								//		dconf->ss_slices_p[idx].logp, dconf->ss_alloc_p[idx] *
								//		sizeof(uint8_t));
								//}

#if defined( DEBUG_SS_SEARCH) || defined(DEBUG_SS_SEARCH2)
								printf("+hit %d + %d = %d, %d, poly (%d,%d) %04x",
									bins1[ii].root[j], bins2[b - 1].root[k], sum, sum - prime, bins1[ii].polynum[j], bins2[b - 1].polynum[k], polysum);
								int y;
								for (y = 0; y < orighits.size; y++)
								{
									if ((orighits.root[y] == sum - prime) && (orighits.polynum[y] == polysum))
									{
										printf("*");
										nummatch++;
										matchm1++;
										break;
						}
					}
								printf("\n");
#endif
								}
							else if ((0 - (sum - prime)) < interval)
							{

								int idx = polymap[polysum];

								if (idx == num_bpoly) continue;

								uint64_t pid = (uint64_t)(i - dconf->ss_slices_n[slice].fboffset) << 24;
								uint64_t sid = (uint64_t)(0 - (sum - prime));
								uint64_t lnext;

								update_ll(1, dconf->ss_slices_n, slice, idx, dconf);

								// the current entry gets this stop marker.
								lnext = 0xffffff0000000000ULL;

								dconf->ss_slices_n[slice].element[dconf->ss_slices_n[slice].size++] =
									(lnext | pid | sid);

								//dconf->ss_slices_n[slice].element[dconf->ss_slices_n[slice].size++] =
								//	((uint64_t)i << 32) | ((uint64_t)(0 - (sum - prime)) << 16) | ((uint64_t)idx);

								if (dconf->ss_slices_n[slice].size >= dconf->ss_slices_n[slice].alloc)
								{
									dconf->ss_slices_n[slice].alloc *= 2;
									dconf->ss_slices_n[slice].element = (uint64_t*)xrealloc(
										dconf->ss_slices_n[slice].element, dconf->ss_slices_n[slice].alloc *
										sizeof(uint64_t));
									//printf("ss bucket allocation now %d elements\n", dconf->ss_slices_p[slice].alloc);
								}

								//dconf->ss_slices_n[idx].root[dconf->ss_size_n[idx]] = 0 - (sum - prime);
								//dconf->ss_slices_n[idx].fbid[dconf->ss_size_n[idx]] = i;
								//dconf->ss_slices_n[idx].logp[dconf->ss_size_n[idx]] = logp;
								//dconf->ss_size_n[idx]++;

								matchm1++;

#if defined( DEBUG_SS_SEARCH)
								int result = check_poly_at_loc(0 - (sum - prime), 1, prime, polysum, idx, polyv, polysign, dconf, sconf);
								if (result == 0)
								{
									printf("more info on failure:\n");
									printf("Bindex %d <- polymap[%04x]\n", idx, polysum);
									printf("bin1 root = %d, bin1 polynum = %04x\n", bins1[ii].root[j], bins1[ii].polynum[j]);
									printf("bin2 root = %d, bin2 polynum = %04x\n", bins2[b - 1].root[k], bins2[b - 1].polynum[k]);
									printf("bin1: %d, bin2: %d, left match, neg side\n", ii, b - 1);
									if (polysum == 0) exit(0);
							}
#endif

								//if (dconf->ss_size_n[idx] > dconf->ss_alloc_n[idx])
								//{
								//	dconf->ss_alloc_n[idx] *= 2;
								//	printf("reallocating ss_slices_n[%d] to %d\n",
								//		idx, dconf->ss_alloc_n[idx]);
								//	dconf->ss_slices_n[idx].root = (uint32_t*)xrealloc(
								//		dconf->ss_slices_n[idx].root, dconf->ss_alloc_n[idx] *
								//		sizeof(uint32_t));
								//	dconf->ss_slices_n[idx].fbid = (uint32_t*)xrealloc(
								//		dconf->ss_slices_n[idx].fbid, dconf->ss_alloc_n[idx] *
								//		sizeof(uint32_t));
								//	dconf->ss_slices_n[idx].logp = (uint8_t*)xrealloc(
								//		dconf->ss_slices_n[idx].logp, dconf->ss_alloc_n[idx] *
								//		sizeof(uint8_t));
								//}

#if defined( DEBUG_SS_SEARCH) || defined(DEBUG_SS_SEARCH2)
								printf("-hit %d + %d = %d, %d, poly (%d,%d) %04x",
									bins1[ii].root[j], bins2[b - 1].root[k], sum, 0 - (sum - prime), bins1[ii].polynum[j], bins2[b - 1].polynum[k], polysum);
								int y;
								for (y = 0; y < orighits.size; y++)
								{
									if ((orighits.root[y] == (0 - (sum - prime))) && (orighits.polynum[y] == polysum))
									{
										printf("*");
										nummatch++;
										matchm1++;
										break;
										}
									}
								printf("\n");
#endif
								}
							}
						}
#endif
					}
				}

			//printf("found %d sieve hits matching bin %d to bin %d\n", nummatch, ii, b);
			//printf("found %d sieve hits matching bin %d to bin %d\n", matchp1, ii, b + 1);
			//printf("found %d sieve hits matching bin %d to bin %d\n", matchm1, ii, b - 1);
			}

		gettimeofday(&stop, NULL);
		t_match_roots += ytools_difftime(&start, &stop);

#if defined( DEBUG_SS_SEARCH) || defined(DEBUG_SS_SEARCH2)
		printf("found %d matches out of %d original hits for prime %d (%u)\n",
			nummatch, orighits.size, i, prime);

		//exit(1);
#endif

	}
#endif

	// now sort the bucket elements by polynum
	gettimeofday(&start, NULL);

	//for (i = 0; i < sconf->factor_base->num_ss_slices; i++)
	//{
	//	 printf("sorting %d elements in +side bucket %d\n", dconf->ss_slices_p[i].size, i);
	//	 qsort(dconf->ss_slices_p[i].element, dconf->ss_slices_p[i].size, sizeof(uint64_t), &cmp_ss_element);
	//	 printf("sorting %d elements in -side bucket %d\n", dconf->ss_slices_n[i].size, i);
	//	 qsort(dconf->ss_slices_n[i].element, dconf->ss_slices_n[i].size, sizeof(uint64_t), &cmp_ss_element);
	//}

	gettimeofday(&stop, NULL);
	t_sort_buckets += ytools_difftime(&start, &stop);

#ifdef USE_SS_SEARCH


	//printf("enumerating roots: %1.4f seconds\n", t_enum_roots);
	//printf("sorting roots: %1.4f seconds\n", t_sort_roots);
	//printf("matching roots: %1.4f seconds\n", t_match_roots);
	//printf("sorting buckets: %1.4f seconds\n", t_sort_buckets);
	//printf("%d max bins\n", numbins);

	//exit(0);

	mpz_clear(polyb1);
	mpz_clear(polyb2);
	free(polymap);
	free(polyv);
	free(polysign);
	free(polynums);

#ifdef DEBUG_SS_SEARCH
	free(orighits.root);
	free(orighits.polynum);
#endif
#endif

	return;
}
#endif

#ifdef USE_SORTED_LIST_SS
void ss_search_sorted_lists(static_conf_t* sconf, dynamic_conf_t* dconf)
{
	// the subset-sum search algorithm using sorted lists to
	// organize the sieve-hits per poly.  
	// Sorted lists makes sieving very fast because it is just
	// a linear dump from a sorted array into the sieve, and
	// SIMD can even help a little with that.
	// But sorting itself is very slow because the lists are so huge.
	// An idea is to not fully sort the whole huge list; instead
	// sort chunks of it.  Then we have several smaller linear
	// dumps into the sieve array and sorting the chunks can
	// stay in cache and hopefully be faster overall.
	siqs_poly* poly = dconf->curr_poly;
	fb_list* fb = sconf->factor_base;
	uint32_t numblocks = sconf->num_blocks;
	uint32_t interval;

	numblocks = sconf->num_blocks;
	interval = numblocks << 15;

#ifdef USE_SS_SEARCH

	int i, j, ii;
	int ss_num_poly_terms = poly->s / 2;
	ss_set_t ss_set1a, ss_set1b, ss_set2a, ss_set2b;

	mpz_t polyb1, polyb2;
	mpz_init(polyb1);
	mpz_init(polyb2);



#ifdef DEBUG_SS_SEARCH
	printf("allocating %u bytes for set1a/b of size %d\n",
		(1 << ss_num_poly_terms) * sizeof(int) * 2 * 2, ss_num_poly_terms);
#endif
	ss_set1a.root = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
	ss_set1a.polynum = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
	ss_set1a.size = ss_num_poly_terms;

	ss_set1b.root = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
	ss_set1b.polynum = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
	ss_set1b.size = ss_num_poly_terms;

	//gmp_printf("%Zd\n", dconf->Bl[0]);
	mpz_set(polyb1, dconf->Bl[0]);

	for (ii = 1; ii < ss_num_poly_terms; ii++) {
		//gmp_printf("%Zd\n", dconf->Bl[ii]);
		mpz_add(polyb1, polyb1, dconf->Bl[ii]);
	}
	mpz_tdiv_q_2exp(polyb1, polyb1, 1);

	ss_num_poly_terms = (poly->s - 1) - ss_num_poly_terms;

#ifdef DEBUG_SS_SEARCH
	printf("allocating %u bytes for set2a/b of size %d\n",
		(1 << ss_num_poly_terms) * sizeof(int) * 2 * 2, ss_num_poly_terms);
#endif
	ss_set2a.root = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
	ss_set2a.polynum = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
	ss_set2a.size = ss_num_poly_terms;

	ss_set2b.root = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
	ss_set2b.polynum = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
	ss_set2b.size = ss_num_poly_terms;

	// first poly b: sum of first n positive Bl
	//gmp_printf("%Zd\n", dconf->Bl[size1]);
	mpz_set(polyb2, dconf->Bl[ii]);

	for (ii++; ii < poly->s; ii++) {
		//gmp_printf("%Zd\n", dconf->Bl[size1 + ii]);
		mpz_add(polyb2, polyb2, dconf->Bl[ii]);
	}
	mpz_tdiv_q_2exp(polyb2, polyb2, 1);

	int numbins = 2 * fb->list->prime[fb->B - 1] / (interval)+1;
	int bindepth = 128;
	int binsize; // = p / numbins + 1;

	ss_set_t* bins1;
	ss_set_t* bins2;

	// init bins
	bins1 = (ss_set_t*)xmalloc(numbins * sizeof(ss_set_t));
	bins2 = (ss_set_t*)xmalloc(numbins * sizeof(ss_set_t));

	for (ii = 0; ii < numbins; ii++)
	{
		bins1[ii].root = (int*)xmalloc(bindepth * sizeof(int));
		bins2[ii].root = (int*)xmalloc(bindepth * sizeof(int));
		bins1[ii].polynum = (int*)xmalloc(bindepth * sizeof(int));
		bins2[ii].polynum = (int*)xmalloc(bindepth * sizeof(int));
		bins1[ii].alloc = bindepth;
		bins2[ii].alloc = bindepth;
		bins1[ii].size = 0;
		bins2[ii].size = 0;
	}

	//printf("allocated %d bins of depth %d\n", numbins, bindepth);
	//exit(0);

#ifdef DEBUG_SS_SEARCH
	mpz_add(dconf->gmptmp1, polyb2, polyb1);

	if (sconf->knmod8 == 1)
	{
		if (sconf->charcount == 1)
		{
			printf("including polya in polyb sum\n");
			mpz_add(dconf->gmptmp1, dconf->gmptmp1, poly->mpz_poly_a);
		}
	}

	if (mpz_cmp(dconf->gmptmp1, poly->mpz_poly_b) == 0)
	{
		prime = fb->list->prime[fb->B - 1];
		gmp_printf("polyb = %Zd + %Zd = %Zd\npolyb matches! (%u, %u, %u, %u)\n",
			polyb1, polyb2, dconf->gmptmp1,
			prime,
			mpz_tdiv_ui(polyb1, prime),
			mpz_tdiv_ui(polyb2, prime),
			mpz_tdiv_ui(poly->mpz_poly_a, prime));
	}
	else
	{
		gmp_printf("polyb = %Zd\npolyb doesn't match :(\n", dconf->gmptmp1);
	}
#endif

#if defined( DEBUG_SS_SEARCH) || defined(DEBUG_SS_SEARCH2)
	ss_set_t orighits;
	orighits.root = (int*)xmalloc(512 * sizeof(int));
	orighits.polynum = (int*)xmalloc(512 * sizeof(int));
#endif

	// create a map from poly number back to 
	// gray code enumeration order.
	int polynum = 0;
	int* polymap = (int*)xmalloc((1 << (poly->s - 1)) * sizeof(int));
	int* polynums = (int*)xmalloc((1 << (poly->s - 1)) * sizeof(int));
	int* polyv = (int*)xmalloc((1 << (poly->s - 1)) * sizeof(int));
	int* polysign = (int*)xmalloc((1 << (poly->s - 1)) * sizeof(int));

	polymap[0] = 1;
	polynums[0] = 0;

	//printf("%d <- %d (%d)\n", 1, 0, polynums[0]);

	for (ii = 1, polynum = 0; ii < (1 << (poly->s - 1)); ii++)
	{
		int v;
		int tmp;
		int sign;

		// next poly enumeration
		v = 1;
		j = ii;
		while ((j & 1) == 0)
			v++, j >>= 1;
		tmp = ii + (1 << v) - 1;
		tmp = (tmp >> v);
		polyv[ii] = v;
		if (tmp & 1)
			polysign[ii] = -1;
		else
			polysign[ii] = 1;

		// next polynum
		polynum ^= (1 << (v - 1));
		polynums[ii] = polynum;
		polymap[polynum] = ii + 1;

		//printf("%d <- %d (%d)\n", ii + 1, polynum, polynums[ii]);
	}

	//exit(0);
	int num_bpoly = 1 << (dconf->curr_poly->s - 1);

	for (i = 0; i < dconf->num_ss_slices; i++)
	{
		dconf->ss_slices_p[i].size = 0;
		dconf->ss_slices_n[i].size = 0;
		dconf->ss_slices_p[i].curr_poly_idx = 0;
		dconf->ss_slices_n[i].curr_poly_idx = 0;
		dconf->ss_slices_p[i].curr_poly_num = 0;
		dconf->ss_slices_n[i].curr_poly_num = 0;
	}

	double t_enum_roots;
	double t_sort_roots;
	double t_match_roots;
	double t_sort_buckets;
	struct timeval start, stop;


	// todo
	// allocate buckets for polysum hits, one for each b-poly and per sieve side.
	// X stop bucket sieving at the ss_start_B index.
	// X set the ss_start_B index in static poly init.
	// need a bucket dumping routine somewhere in tdiv.  maybe tdiv_large.
	// testing!

	// optimizations:
	// can still use vector root updates when enumerating roots
	// for each side of the sum.  just need to be careful to
	// do root matching per prime in that case.
	// probably can do vector match search since we take an element
	// on one side then sum and compare to many elements on the other side.
	// could either do this many primes at a time or maybe many elements per
	// prime at a time.
	// how many "other side" bins we need to search for matches needs to
	// be more exact than hardcoded +/- 1 bins from the center bin.
	int* rootupdates = dconf->rootupdates;
	update_t update_data = dconf->update_data;
	uint32_t* modsqrt = sconf->modsqrt_array;

	for (i = fb->x2_large_B; i < fb->B; i++)
	{
		int numB;
		int maxB = 1 << (dconf->curr_poly->s - 1);
		int r1 = update_data.firstroots1[i], r2 = update_data.firstroots2[i];
		uint32_t bound = sconf->factor_base->B;
		uint32_t med_B = sconf->factor_base->med_B;
		int polynum = 0;
		int slice = (int)((i - sconf->factor_base->x2_large_B) / 65536);

		int prime = fb->list->prime[i];
		int p = prime;
		int root1 = modsqrt[i];
		int root2 = prime - root1;
		uint8_t logp = fb->list->logprime[i];

		int amodp = (int)mpz_tdiv_ui(poly->mpz_poly_a, prime);
		int bmodp = (int)mpz_tdiv_ui(poly->mpz_poly_b, prime);
		int inv;

		//find a^-1 mod p = inv(a mod p) mod p
		if (sconf->knmod8 == 1)
		{
			inv = modinv_1(2 * amodp, prime);
		}
		else
		{
			inv = modinv_1(amodp, prime);
		}

#if defined( DEBUG_SS_SEARCH) || defined(DEBUG_SS_SEARCH2)
		orighits.size = 0;

		if ((r1 < interval) || (r2 < interval) || ((prime - r1) < interval) || ((prime - r2) < interval))
		{
			//printf("%04d,%04x: initial roots = %d,%d (%d,%d)\n", 0, polynum, r1, r2, prime - r1, prime - r2);
			if (r1 < interval) {
				orighits.root[orighits.size] = r1;
				orighits.polynum[orighits.size] = polynum;
				orighits.size++;
			}
			if (r2 < interval) {
				orighits.root[orighits.size] = r2;
				orighits.polynum[orighits.size] = polynum;
				orighits.size++;
			}
			if ((prime - r1) < interval) {
				orighits.root[orighits.size] = prime - r1;
				orighits.polynum[orighits.size] = polynum;
				orighits.size++;
			}
			if ((prime - r2) < interval) {
				orighits.root[orighits.size] = prime - r2;
				orighits.polynum[orighits.size] = polynum;
				orighits.size++;
			}
		}

		for (numB = 1; numB < maxB; numB++)
		{
			char v = dconf->curr_poly->nu[numB];
			char sign = dconf->curr_poly->gray[numB];
			int* ptr = &rootupdates[(v - 1) * bound];
			char hit1, hit2, hit3, hit4;

			if (sign > 0)
			{
				r1 -= ptr[i];
				r2 -= ptr[i];
				r1 = (r1 < 0) ? r1 + prime : r1;
				r2 = (r2 < 0) ? r2 + prime : r2;
			}
			else
			{
				r1 += ptr[i];
				r2 += ptr[i];
				r1 = (r1 >= prime) ? r1 - prime : r1;
				r2 = (r2 >= prime) ? r2 - prime : r2;
			}

			polynum ^= (1 << (v - 1));

			hit1 = r1 < interval ? '*' : ' ';
			hit2 = r2 < interval ? '*' : ' ';
			hit3 = (prime - r1) < interval ? '*' : ' ';
			hit4 = (prime - r2) < interval ? '*' : ' ';

			if ((r1 < interval) || (r2 < interval) || ((prime - r1) < interval) || ((prime - r2) < interval))
			{
				//printf("%04d, %04x: ptr = %09d, v = %03d, s = %d: roots = %d%c,%d%c : %d%c,%d%c\n",
				//	numB, polynum, ptr[i], v, sign > 0 ? 1 : 0, r1, hit1, r2, hit2, prime - r1, hit3, prime - r2, hit4);

				if (r1 < interval) {
					orighits.root[orighits.size] = r1;
					orighits.polynum[orighits.size] = polynum;
					orighits.size++;
				}
				if (r2 < interval) {
					orighits.root[orighits.size] = r2;
					orighits.polynum[orighits.size] = polynum;
					orighits.size++;
				}
				if ((prime - r1) < interval) {
					orighits.root[orighits.size] = prime - r1;
					orighits.polynum[orighits.size] = polynum;
					orighits.size++;
				}
				if ((prime - r2) < interval) {
					orighits.root[orighits.size] = prime - r2;
					orighits.polynum[orighits.size] = polynum;
					orighits.size++;
				}

			}
		}

#endif

		// create set1, an enumeration of (+/-t - b)(a)^-1 mod p
		// for b in the set [+/-B1 +/- B2 +/- ... Bs/2]
		int ii, v, sign, n = poly->s / 2;
		int tmp;

		gettimeofday(&start, NULL);


		prime = fb->list->prime[i];
		root1 = modsqrt[i];
		root2 = prime - root1;

		//printf("new prime %d, r1 = %d, r2 = %d\n", prime, root1, root2);

#ifdef DEBUG_SS_SEARCH
			//printf("p,t1,t2 = %d,%d,%d\n", prime, root1, root2);
		printf("building set1 B\n");
#endif

		// first poly b: sum of first n positive Bl
		bmodp = (int)mpz_tdiv_ui(polyb1, prime);

		// under these conditions polya is added to polyb.
		// here we account for that by adding amodp into bmodp.
		// only here in set1 enumeration.
		if (sconf->knmod8 == 1)
		{
			if (sconf->charcount == 1)
			{
				bmodp += amodp;
				if (bmodp >= prime)
				{
					bmodp -= prime;
				}
			}
		}

#ifdef DEBUG_SS_SEARCH
		int sumbmodp = bmodp;
		mpz_set(dconf->gmptmp1, polyb1);
#endif

		//COMPUTE_FIRST_ROOTS;
		r1 = (int)root1 - bmodp;
		r2 = (int)root2 - bmodp;
		if (r1 < 0) r1 += prime;
		if (r2 < 0) r2 += prime;

		r1 = (int)((uint64_t)inv * (uint64_t)r1 % (uint64_t)prime);
		r2 = (int)((uint64_t)inv * (uint64_t)r2 % (uint64_t)prime);

#ifdef DEBUG_SS_SEARCH
		printf("\np = %d, a^-1 mod p = %d, roots = %d,%d, s = %d\n", prime, inv, r1, r2, n);
		gmp_printf("b = %Zd, amodp = %d, bmodp = %d\n", polyb1, amodp, bmodp);
#endif

		int* ptr;

#ifdef DEBUG_SS_SEARCH
		printf("enumerating %d elements of set 1 with %d polyb terms\n", (1 << n), n);
#endif

		// enumerate set1 roots
		for (ii = 1, polynum = 0; ii < (1 << n); ii++) {
			// these roots go into the set
			ss_set1a.root[ii - 1] = r1;
			ss_set1b.root[ii - 1] = r2;
			ss_set1a.polynum[ii - 1] = polynum;
			ss_set1b.polynum[ii - 1] = polynum;

			// next polynum
			polynum = polynums[ii];
			sign = polysign[ii];
			v = polyv[ii];

			// next roots
			ptr = &rootupdates[(v - 1) * bound + i];
			if (sign > 0)
			{
				r1 = (int)r1 - *ptr;
				r2 = (int)r2 - *ptr;
				r1 = (r1 < 0) ? r1 + prime : r1;
				r2 = (r2 < 0) ? r2 + prime : r2;
			}
			else
			{
				r1 = (int)r1 + *ptr;
				r2 = (int)r2 + *ptr;
				r1 = (r1 >= prime) ? r1 - prime : r1;
				r2 = (r2 >= prime) ? r2 - prime : r2;
			}
		}
		// these roots go into the set
		ss_set1a.root[ii - 1] = r1;
		ss_set1a.polynum[ii - 1] = polynum;
		ss_set1b.root[ii - 1] = r2;
		ss_set1b.polynum[ii - 1] = polynum;

		int size1 = n;

		// create set2, an enumeration of (+/-t - b)(a)^-1 mod p
		// for b in the set [+/-Bs/2+1 +/- Bs/2+2 ... + Bs]
		n = (poly->s - 1) - n;

#ifdef DEBUG_SS_SEARCH
		printf("building set2 B\n");
#endif
		bmodp = (int)mpz_tdiv_ui(polyb2, prime);

#ifdef DEBUG_SS_SEARCH
		sumbmodp += bmodp;
		if (sumbmodp >= prime) sumbmodp -= prime;
#endif

		// - (B1 + B2 + B3 + ... Bs-1)
		bmodp = prime - bmodp;

		r1 = (uint32_t)((uint64_t)inv * (uint64_t)bmodp % (uint64_t)prime);

#ifdef DEBUG_SS_SEARCH
		printf("\np = %d, a^-1 mod p = %d, roots = %d,%d, s = %d\n", prime, inv, r1, r2, n);
		gmp_printf("b = %Zd, amodp = %d, bmodp = %d\n", polyb2, amodp, bmodp);

		gmp_printf("orig polya: %Zd\n", poly->mpz_poly_a);
		gmp_printf("orig polyb: %Zd\n", poly->mpz_poly_b);
		int oamodp = mpz_tdiv_ui(poly->mpz_poly_a, prime);
		int obmodp = mpz_tdiv_ui(poly->mpz_poly_b, prime);
		printf("orig amodp = %d, bmodp = %d, root1 = %d, root2 = %d\n",
			oamodp, obmodp, update_data.firstroots1[i], update_data.firstroots2[i]);

		if (obmodp == sumbmodp)
		{
			printf("bmodp matches!\n");
		}
		else
		{
			printf("bmodp doesn't match :(\n");
		}


		printf("enumerating %d elements of set 2 with %d polyb terms\n", (1 << n), n);
#endif

		// enumerate set2 roots
		for (ii = 1, polynum = 0; ii < (1 << n); ii++) {
			// these roots go into the set
			ss_set2a.root[ii - 1] = r1;
			ss_set2a.polynum[ii - 1] = polynum;

			polynum = polynums[ii] << size1;
			sign = polysign[ii];
			v = polyv[ii];

			// next roots
			ptr = &rootupdates[(size1 + v - 1) * bound + i];
			if (sign > 0)
			{
				r1 = (int)r1 - *ptr;
				r1 = (r1 < 0) ? r1 + prime : r1;
			}
			else
			{
				r1 = (int)r1 + *ptr;
				r1 = (r1 >= prime) ? r1 - prime : r1;
			}
		}

		// these roots go into the set
		ss_set2a.root[ii - 1] = r1;
		ss_set2a.polynum[ii - 1] = polynum;


		gettimeofday(&stop, NULL);
		t_enum_roots += ytools_difftime(&start, &stop);

		gettimeofday(&start, NULL);

		// now sort the sets into a moderate number of bins over the range 0:p
		numbins = p / (2 * interval) + 1;
		bindepth = 128;
		binsize = p / numbins + 1;

#ifdef DEBUG_SS_SEARCH
		printf("sorting rootset1 of size %d into %d bins of size %d\n", (1 << ss_set1a.size), numbins, binsize);
#endif

		for (ii = 0; ii < numbins; ii++)
		{
			bins1[ii].size = 0;
			bins2[ii].size = 0;
		}

		// sort root 1 into set 1
		for (ii = 0; ii < (1 << ss_set1a.size); ii++)
		{
			int binnum = ss_set1a.root[ii] / binsize;
			if (binnum < numbins)
			{
				//printf("bin %d <-- set1a root %d (L polynum %d)\n", binnum, ss_set1a.root[ii], ss_set1a.polynum[ii]);
				bins1[binnum].root[bins1[binnum].size] = ss_set1a.root[ii];
				bins1[binnum].polynum[bins1[binnum].size] = ss_set1a.polynum[ii];
				bins1[binnum].size++;
				if (bins1[binnum].size > bindepth)
				{
					printf("bin overflow\n");
					exit(1);
				}
			}
			else
			{
				printf("element %d of set 1, root %d, invalid bin %d [of %d]\n",
					ii, ss_set1a.root[ii], binnum, numbins);
			}
		}

#ifdef DEBUG_SS_SEARCH
		printf("sorting rootset2 of size %d into %d bins of size %d\n", (1 << ss_set2a.size), numbins, binsize);
#endif

		// sort set 2
		for (ii = 0; ii < (1 << ss_set2a.size); ii++)
		{
			int binnum = ss_set2a.root[ii] / binsize;
			if (binnum < numbins)
			{
				//printf("bin %d <-- set2a root %d (R polynum %d)\n", binnum, ss_set2a.root[ii], ss_set2a.polynum[ii]);
				bins2[binnum].root[bins2[binnum].size] = ss_set2a.root[ii];
				bins2[binnum].polynum[bins2[binnum].size] = ss_set2a.polynum[ii];
				bins2[binnum].size++;
				if (bins2[binnum].size > bindepth)
				{
					printf("bin overflow\n");
					exit(1);
				}
			}
			else
			{
				printf("element %d of set 2, root %d, invalid bin %d [of %d]\n",
					ii, ss_set2a.root[ii], binnum, numbins);
			}
		}

		// now match all numbers in bin x with the other set's bin that
		// contains p - x
#ifdef DEBUG_SS_SEARCH
		printf("matching root1\n");
#endif


		gettimeofday(&stop, NULL);
		t_sort_roots += ytools_difftime(&start, &stop);

		gettimeofday(&start, NULL);

		// match set1 holding root1 to set2
		int nummatch = 0;
		int matchp1 = 0;
		int matchm1 = 0;
		for (ii = 0; ii < numbins; ii++)
		{
			int x = binsize * ii + binsize / 2;
			int px = p - x;
			int b = px / binsize; // px / binsize;

			nummatch = 0;
			matchp1 = 0;
			matchm1 = 0;

			{
				for (j = 0; j < bins1[ii].size; j++)
				{
					int k;
					for (k = 0; k < bins2[b].size; k++)
					{
						int sum = bins1[ii].root[j] + bins2[b].root[k];
						int polysum = bins1[ii].polynum[j] + bins2[b].polynum[k];

						if (polysum >= (1 << (poly->s - 1)))
						{
							printf("invalid polysum %d (%d + %d)\n",
								polysum, bins1[ii].polynum[j], bins2[b].polynum[k]);
							continue;
						}

						if ((sum >= prime) && ((sum - prime) < interval))
						{
							int idx = polymap[polysum];

							if (idx == num_bpoly) continue;

							uint64_t pid = (uint64_t)(i - dconf->ss_slices_p[slice].fboffset) << 24;
							uint64_t sid = (uint64_t)(sum - prime);
							uint64_t polyid = (uint64_t)idx << 40;

							dconf->ss_slices_p[slice].element[dconf->ss_slices_p[slice].size++] =
								(polyid | pid | sid);

							if (dconf->ss_slices_p[slice].size >= dconf->ss_slices_p[slice].alloc)
							{
								dconf->ss_slices_p[slice].alloc *= 2;
								dconf->ss_slices_p[slice].element = (uint64_t*)xrealloc(
									dconf->ss_slices_p[slice].element, dconf->ss_slices_p[slice].alloc *
									sizeof(uint64_t));
							}

							nummatch++;

						}
						else if ((sum < prime) && ((0 - (sum - prime)) < interval))
						{
							int idx = polymap[polysum];

							if (idx == num_bpoly) continue;

							uint64_t pid = (uint64_t)(i - dconf->ss_slices_n[slice].fboffset) << 24;
							uint64_t sid = (uint64_t)(0 - (sum - prime));
							uint64_t polyid = (uint64_t)idx << 40;

							dconf->ss_slices_n[slice].element[dconf->ss_slices_n[slice].size++] =
								(polyid | pid | sid);

							if (dconf->ss_slices_n[slice].size >= dconf->ss_slices_n[slice].alloc)
							{
								dconf->ss_slices_n[slice].alloc *= 2;
								dconf->ss_slices_n[slice].element = (uint64_t*)xrealloc(
									dconf->ss_slices_n[slice].element, dconf->ss_slices_n[slice].alloc *
									sizeof(uint64_t));
								//printf("ss bucket allocation now %d elements\n", dconf->ss_slices_p[slice].alloc);
							}

							nummatch++;
						}
					}
#if 1
					if ((b + 1) < numbins)
					{
						for (k = 0; k < bins2[b + 1].size; k++)
						{
							int sum = bins1[ii].root[j] + bins2[b + 1].root[k];
							int polysum = bins1[ii].polynum[j] + bins2[b + 1].polynum[k];
							if ((sum >= prime) && ((sum - prime) < interval)) {

								int idx = polymap[polysum];

								if (idx == num_bpoly) continue;

								uint64_t pid = (uint64_t)(i - dconf->ss_slices_p[slice].fboffset) << 24;
								uint64_t sid = (uint64_t)(sum - prime);
								uint64_t polyid = (uint64_t)idx << 40;

								dconf->ss_slices_p[slice].element[dconf->ss_slices_p[slice].size++] =
									(polyid | pid | sid);

								if (dconf->ss_slices_p[slice].size >= dconf->ss_slices_p[slice].alloc)
								{
									dconf->ss_slices_p[slice].alloc *= 2;
									dconf->ss_slices_p[slice].element = (uint64_t*)xrealloc(
										dconf->ss_slices_p[slice].element, dconf->ss_slices_p[slice].alloc *
										sizeof(uint64_t));
									//printf("ss bucket allocation now %d elements\n", dconf->ss_slices_p[slice].alloc);
								}

								matchp1++;
							}
							else if ((0 - (sum - prime)) < interval)
							{

								int idx = polymap[polysum];

								if (idx == num_bpoly) continue;

								uint64_t pid = (uint64_t)(i - dconf->ss_slices_n[slice].fboffset) << 24;
								uint64_t sid = (uint64_t)(0 - (sum - prime));
								uint64_t polyid = (uint64_t)idx << 40;

								dconf->ss_slices_n[slice].element[dconf->ss_slices_n[slice].size++] =
									(polyid | pid | sid);

								if (dconf->ss_slices_n[slice].size >= dconf->ss_slices_n[slice].alloc)
								{
									dconf->ss_slices_n[slice].alloc *= 2;
									dconf->ss_slices_n[slice].element = (uint64_t*)xrealloc(
										dconf->ss_slices_n[slice].element, dconf->ss_slices_n[slice].alloc *
										sizeof(uint64_t));
									//printf("ss bucket allocation now %d elements\n", dconf->ss_slices_p[slice].alloc);
								}

								matchp1++;
							}
						}
					}
					if ((b - 1) >= 0)
					{
						for (k = 0; k < bins2[b - 1].size; k++)
						{
							int sum = bins1[ii].root[j] + bins2[b - 1].root[k];
							int polysum = bins1[ii].polynum[j] + bins2[b - 1].polynum[k];
							if ((sum >= prime) && ((sum - prime) < interval)) {

								int idx = polymap[polysum];

								if (idx == num_bpoly) continue;

								uint64_t pid = (uint64_t)(i - dconf->ss_slices_p[slice].fboffset) << 24;
								uint64_t sid = (uint64_t)(sum - prime);
								uint64_t polyid = (uint64_t)idx << 40;

								dconf->ss_slices_p[slice].element[dconf->ss_slices_p[slice].size++] =
									(polyid | pid | sid);

								if (dconf->ss_slices_p[slice].size >= dconf->ss_slices_p[slice].alloc)
								{
									dconf->ss_slices_p[slice].alloc *= 2;
									dconf->ss_slices_p[slice].element = (uint64_t*)xrealloc(
										dconf->ss_slices_p[slice].element, dconf->ss_slices_p[slice].alloc *
										sizeof(uint64_t));
									//printf("ss bucket allocation now %d elements\n", dconf->ss_slices_p[slice].alloc);
								}

								matchm1++;
							}
							else if ((0 - (sum - prime)) < interval)
							{

								int idx = polymap[polysum];

								if (idx == num_bpoly) continue;

								uint64_t pid = (uint64_t)(i - dconf->ss_slices_n[slice].fboffset) << 24;
								uint64_t sid = (uint64_t)(0 - (sum - prime));
								uint64_t polyid = (uint64_t)idx << 40;

								dconf->ss_slices_n[slice].element[dconf->ss_slices_n[slice].size++] =
									(polyid | pid | sid);

								if (dconf->ss_slices_n[slice].size >= dconf->ss_slices_n[slice].alloc)
								{
									dconf->ss_slices_n[slice].alloc *= 2;
									dconf->ss_slices_n[slice].element = (uint64_t*)xrealloc(
										dconf->ss_slices_n[slice].element, dconf->ss_slices_n[slice].alloc *
										sizeof(uint64_t));
									//printf("ss bucket allocation now %d elements\n", dconf->ss_slices_p[slice].alloc);
								}

								matchm1++;
							}
						}
					}
#endif
				}
			}

			//printf("found %d sieve hits matching bin %d to bin %d\n", nummatch, ii, b);
			//printf("found %d sieve hits matching bin %d to bin %d\n", matchp1, ii, b + 1);
			//printf("found %d sieve hits matching bin %d to bin %d\n", matchm1, ii, b - 1);

		}

		// reset set1 so we can put root2 into it
		for (ii = 0; ii < numbins; ii++)
		{
			bins1[ii].size = 0;
		}

#ifdef DEBUG_SS_SEARCH
		printf("sorting rootset1b of size %d into %d bins of size %d\n", (1 << ss_set1b.size), numbins, binsize);
#endif

		// sort root2 into set1
		for (ii = 0; ii < (1 << ss_set1b.size); ii++)
		{
			int binnum = ss_set1b.root[ii] / binsize;
			if (binnum < numbins)
			{
				// printf("bin %d <-- set1b root %d (polynum 0)\n", binnum, ss_set1a.root[ii]);
				bins1[binnum].root[bins1[binnum].size] = ss_set1b.root[ii];
				bins1[binnum].polynum[bins1[binnum].size] = ss_set1b.polynum[ii];
				bins1[binnum].size++;
				if (bins1[binnum].size > bindepth)
				{
					printf("bin overflow\n");
					exit(1);
				}
			}
			else
			{
				printf("element %d of set 1, root %d, invalid bin %d [of %d]\n",
					ii, ss_set1b.root[ii], binnum, numbins);
			}
		}


#ifdef DEBUG_SS_SEARCH
		printf("sorting rootset2b of size %d into %d bins of size %d\n", (1 << ss_set2b.size), numbins, binsize);
#endif

		// now match all numbers in bin x with the other set's bin that
		// contains p - x
#ifdef DEBUG_SS_SEARCH
		printf("matching root2\n");
#endif

		// match set1 holding root2 to set2
		for (ii = 0; ii < numbins; ii++)
		{
			int x = binsize * ii + binsize / 2;
			int px = p - x;
			int b = px / binsize; // px / binsize;

			nummatch = 0;
			matchp1 = 0;
			matchm1 = 0;

			{
				for (j = 0; j < bins1[ii].size; j++)
				{
					int k;
					for (k = 0; k < bins2[b].size; k++)
					{
						int sum = bins1[ii].root[j] + bins2[b].root[k];
						int polysum = bins1[ii].polynum[j] + bins2[b].polynum[k];

						if (polysum >= (1 << (poly->s - 1)))
						{
							printf("invalid polysum %d (%d + %d)\n",
								polysum, bins1[ii].polynum[j], bins2[b].polynum[k]);
							continue;
						}

						if ((sum >= prime) && ((sum - prime) < interval))
						{
							int idx = polymap[polysum];

							if (idx == num_bpoly) continue;

							uint64_t pid = (uint64_t)(i - dconf->ss_slices_p[slice].fboffset) << 24;
							uint64_t sid = (uint64_t)(sum - prime);
							uint64_t polyid = (uint64_t)idx << 40;

							dconf->ss_slices_p[slice].element[dconf->ss_slices_p[slice].size++] =
								(polyid | pid | sid);

							if (dconf->ss_slices_p[slice].size >= dconf->ss_slices_p[slice].alloc)
							{
								dconf->ss_slices_p[slice].alloc *= 2;
								dconf->ss_slices_p[slice].element = (uint64_t*)xrealloc(
									dconf->ss_slices_p[slice].element, dconf->ss_slices_p[slice].alloc *
									sizeof(uint64_t));
								//printf("ss bucket allocation now %d elements\n", dconf->ss_slices_p[slice].alloc);
							}

							nummatch++;

						}
						else if ((sum < prime) && ((0 - (sum - prime)) < interval))
						{
							int idx = polymap[polysum];

							if (idx == num_bpoly) continue;

							uint64_t pid = (uint64_t)(i - dconf->ss_slices_n[slice].fboffset) << 24;
							uint64_t sid = (uint64_t)(0 - (sum - prime));
							uint64_t polyid = (uint64_t)idx << 40;

							dconf->ss_slices_n[slice].element[dconf->ss_slices_n[slice].size++] =
								(polyid | pid | sid);

							if (dconf->ss_slices_n[slice].size >= dconf->ss_slices_n[slice].alloc)
							{
								dconf->ss_slices_n[slice].alloc *= 2;
								dconf->ss_slices_n[slice].element = (uint64_t*)xrealloc(
									dconf->ss_slices_n[slice].element, dconf->ss_slices_n[slice].alloc *
									sizeof(uint64_t));
								//printf("ss bucket allocation now %d elements\n", dconf->ss_slices_p[slice].alloc);
							}

							nummatch++;
						}
					}
#if 1
					if ((b + 1) < numbins)
					{
						for (k = 0; k < bins2[b + 1].size; k++)
						{
							int sum = bins1[ii].root[j] + bins2[b + 1].root[k];
							int polysum = bins1[ii].polynum[j] + bins2[b + 1].polynum[k];
							if ((sum >= prime) && ((sum - prime) < interval)) {

								int idx = polymap[polysum];

								if (idx == num_bpoly) continue;

								uint64_t pid = (uint64_t)(i - dconf->ss_slices_p[slice].fboffset) << 24;
								uint64_t sid = (uint64_t)(sum - prime);
								uint64_t polyid = (uint64_t)idx << 40;

								dconf->ss_slices_p[slice].element[dconf->ss_slices_p[slice].size++] =
									(polyid | pid | sid);

								if (dconf->ss_slices_p[slice].size >= dconf->ss_slices_p[slice].alloc)
								{
									dconf->ss_slices_p[slice].alloc *= 2;
									dconf->ss_slices_p[slice].element = (uint64_t*)xrealloc(
										dconf->ss_slices_p[slice].element, dconf->ss_slices_p[slice].alloc *
										sizeof(uint64_t));
									//printf("ss bucket allocation now %d elements\n", dconf->ss_slices_p[slice].alloc);
								}

								matchp1++;

							}
							else if ((0 - (sum - prime)) < interval)
							{

								int idx = polymap[polysum];

								if (idx == num_bpoly) continue;

								uint64_t pid = (uint64_t)(i - dconf->ss_slices_n[slice].fboffset) << 24;
								uint64_t sid = (uint64_t)(0 - (sum - prime));
								uint64_t polyid = (uint64_t)idx << 40;

								dconf->ss_slices_n[slice].element[dconf->ss_slices_n[slice].size++] =
									(polyid | pid | sid);

								if (dconf->ss_slices_n[slice].size >= dconf->ss_slices_n[slice].alloc)
								{
									dconf->ss_slices_n[slice].alloc *= 2;
									dconf->ss_slices_n[slice].element = (uint64_t*)xrealloc(
										dconf->ss_slices_n[slice].element, dconf->ss_slices_n[slice].alloc *
										sizeof(uint64_t));
									//printf("ss bucket allocation now %d elements\n", dconf->ss_slices_p[slice].alloc);
								}

								matchp1++;
							}
						}
					}
					if ((b - 1) >= 0)
					{
						for (k = 0; k < bins2[b - 1].size; k++)
						{
							int sum = bins1[ii].root[j] + bins2[b - 1].root[k];
							int polysum = bins1[ii].polynum[j] + bins2[b - 1].polynum[k];
							if ((sum >= prime) && ((sum - prime) < interval)) {

								int idx = polymap[polysum];

								if (idx == num_bpoly) continue;

								uint64_t pid = (uint64_t)(i - dconf->ss_slices_p[slice].fboffset) << 24;
								uint64_t sid = (uint64_t)(sum - prime);
								uint64_t polyid = (uint64_t)idx << 40;

								dconf->ss_slices_p[slice].element[dconf->ss_slices_p[slice].size++] =
									(polyid | pid | sid);

								if (dconf->ss_slices_p[slice].size >= dconf->ss_slices_p[slice].alloc)
								{
									dconf->ss_slices_p[slice].alloc *= 2;
									dconf->ss_slices_p[slice].element = (uint64_t*)xrealloc(
										dconf->ss_slices_p[slice].element, dconf->ss_slices_p[slice].alloc *
										sizeof(uint64_t));
									//printf("ss bucket allocation now %d elements\n", dconf->ss_slices_p[slice].alloc);
								}

								matchm1++;
							}
							else if ((0 - (sum - prime)) < interval)
							{

								int idx = polymap[polysum];

								if (idx == num_bpoly) continue;

								uint64_t pid = (uint64_t)(i - dconf->ss_slices_n[slice].fboffset) << 24;
								uint64_t sid = (uint64_t)(0 - (sum - prime));
								uint64_t polyid = (uint64_t)idx << 40;

								dconf->ss_slices_n[slice].element[dconf->ss_slices_n[slice].size++] =
									(polyid | pid | sid);

								if (dconf->ss_slices_n[slice].size >= dconf->ss_slices_n[slice].alloc)
								{
									dconf->ss_slices_n[slice].alloc *= 2;
									dconf->ss_slices_n[slice].element = (uint64_t*)xrealloc(
										dconf->ss_slices_n[slice].element, dconf->ss_slices_n[slice].alloc *
										sizeof(uint64_t));
									//printf("ss bucket allocation now %d elements\n", dconf->ss_slices_p[slice].alloc);
								}

								matchm1++;

							}
						}
					}
#endif
				}
			}

			//printf("found %d sieve hits matching bin %d to bin %d\n", nummatch, ii, b);
			//printf("found %d sieve hits matching bin %d to bin %d\n", matchp1, ii, b + 1);
			//printf("found %d sieve hits matching bin %d to bin %d\n", matchm1, ii, b - 1);
		}

		gettimeofday(&stop, NULL);
		t_match_roots += ytools_difftime(&start, &stop);

#if defined( DEBUG_SS_SEARCH) || defined(DEBUG_SS_SEARCH2)
		printf("found %d matches out of %d original hits for prime %d (%u)\n",
			nummatch, orighits.size, i, prime);

		//exit(1);
#endif

	}

	// now sort the bucket elements by polynum
	gettimeofday(&start, NULL);

	for (i = 0; i < sconf->factor_base->num_ss_slices; i++)
	{
		printf("sorting %d elements in +side bucket %d\n", dconf->ss_slices_p[i].size, i);
		qsort(dconf->ss_slices_p[i].element, dconf->ss_slices_p[i].size, sizeof(uint64_t), &cmp_ss_element);
		printf("sorting %d elements in -side bucket %d\n", dconf->ss_slices_n[i].size, i);
		qsort(dconf->ss_slices_n[i].element, dconf->ss_slices_n[i].size, sizeof(uint64_t), &cmp_ss_element);
	}

	gettimeofday(&stop, NULL);
	t_sort_buckets += ytools_difftime(&start, &stop);


	printf("enumerating roots: %1.4f seconds\n", t_enum_roots);
	printf("sorting roots: %1.4f seconds\n", t_sort_roots);
	printf("matching roots: %1.4f seconds\n", t_match_roots);
	printf("sorting buckets: %1.4f seconds\n", t_sort_buckets);

	//exit(0);

	mpz_clear(polyb1);
	mpz_clear(polyb2);
	free(polymap);
	free(polyv);
	free(polysign);
	free(polynums);

#ifdef DEBUG_SS_SEARCH
	free(orighits.root);
	free(orighits.polynum);
#endif


#endif

	return;
}
#endif

#ifdef USE_POLY_BUCKET_SS
void ss_search_poly_buckets(static_conf_t* sconf, dynamic_conf_t* dconf)
{
	// the subset-sum search algorithm using a bucket sort to
	// organize the sieve-hits per poly.  
	// bucket sorting makes sieving very fast because it is just
	// a linear dump from a bucket into the sieve, and
	// SIMD can even help a little with that.
	// But sorting itself is very slow because there are *many* 
	// polynomials and they are widely spaced in memory, making 
	// the sort slow.
	// An idea is sort groups of polynomials so that we aren't 
	// accessing so many buckets during the sort.  Then sieving
	// can use similar bookkeeping to the sorted-lists approach
	// within each poly-group bucket.
	siqs_poly* poly = dconf->curr_poly;
	fb_list* fb = sconf->factor_base;
	uint32_t numblocks = sconf->num_blocks;
	uint32_t interval;

	numblocks = sconf->num_blocks;
	interval = numblocks << 15;

#ifdef USE_SS_SEARCH

	int i, j, ii;
	int ss_num_poly_terms = poly->s / 2;
	ss_set_t ss_set1a, ss_set1b, ss_set2a, ss_set2b;

	mpz_t polyb1, polyb2;
	mpz_init(polyb1);
	mpz_init(polyb2);

	ss_set1a.root = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
	ss_set1a.polynum = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
	ss_set1a.size = ss_num_poly_terms;

	ss_set1b.root = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
	ss_set1b.polynum = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
	ss_set1b.size = ss_num_poly_terms;

	//gmp_printf("%Zd\n", dconf->Bl[0]);
	mpz_set(polyb1, dconf->Bl[0]);

	for (ii = 1; ii < ss_num_poly_terms; ii++) {
		//gmp_printf("%Zd\n", dconf->Bl[ii]);
		mpz_add(polyb1, polyb1, dconf->Bl[ii]);
	}
	mpz_tdiv_q_2exp(polyb1, polyb1, 1);

	ss_num_poly_terms = (poly->s - 1) - ss_num_poly_terms;

	ss_set2a.root = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
	ss_set2a.polynum = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
	ss_set2a.size = ss_num_poly_terms;

	ss_set2b.root = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
	ss_set2b.polynum = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
	ss_set2b.size = ss_num_poly_terms;

	// first poly b: sum of first n positive Bl
	mpz_set(polyb2, dconf->Bl[ii]);

	for (ii++; ii < poly->s; ii++) {
		mpz_add(polyb2, polyb2, dconf->Bl[ii]);
	}
	mpz_tdiv_q_2exp(polyb2, polyb2, 1);

	int numbins = 2 * fb->list->prime[fb->B - 1] / (interval)+1;
	int bindepth = 64;
	int binsize;

	ss_set_t* bins1a;
	ss_set_t* bins1b;
	ss_set_t* bins2;

	// init bins
	bins1a = (ss_set_t*)xmalloc(numbins * sizeof(ss_set_t));
	bins1b = (ss_set_t*)xmalloc(numbins * sizeof(ss_set_t));
	bins2 = (ss_set_t*)xmalloc(numbins * sizeof(ss_set_t));

	for (ii = 0; ii < numbins; ii++)
	{
		bins1a[ii].root = (int*)xmalloc(bindepth * sizeof(int));
		bins1b[ii].root = (int*)xmalloc(bindepth * sizeof(int));
		bins2[ii].root = (int*)xmalloc(bindepth * sizeof(int));
		bins1a[ii].polynum = (int*)xmalloc(bindepth * sizeof(int));
		bins1b[ii].polynum = (int*)xmalloc(bindepth * sizeof(int));
		bins2[ii].polynum = (int*)xmalloc(bindepth * sizeof(int));
		bins1a[ii].alloc = bindepth;
		bins1b[ii].alloc = bindepth;
		bins2[ii].alloc = bindepth;
		bins1a[ii].size = 0;
		bins1b[ii].size = 0;
		bins2[ii].size = 0;
	}

	// create a map from poly number back to 
	// gray code enumeration order.
	int polynum = 0;
	int* polymap = dconf->polymap; // (int*)xmalloc((1 << (poly->s - 1)) * sizeof(int));
	int* polynums = (int*)xmalloc((1 << (poly->s - 1)) * sizeof(int));
	int* polyv = (int*)xmalloc((1 << (poly->s - 1)) * sizeof(int));
	int* polysign = (int*)xmalloc((1 << (poly->s - 1)) * sizeof(int));

	//polymap[0] = 1;
	polymap[1] = 0;
	polynums[0] = 0;

	for (ii = 1, polynum = 0; ii < (1 << (poly->s - 1)); ii++)
	{
		int v;
		int tmp;
		int sign;

		// next poly enumeration
		v = 1;
		j = ii;
		while ((j & 1) == 0)
			v++, j >>= 1;
		tmp = ii + (1 << v) - 1;
		tmp = (tmp >> v);
		polyv[ii] = v;
		if (tmp & 1)
			polysign[ii] = -1;
		else
			polysign[ii] = 1;

		// next polynum
		polynum ^= (1 << (v - 1));
		polynums[ii] = polynum;
		//polymap[polynum] = ii + 1;
		polymap[ii + 1] = polynum;

		//printf("%d <- %d (%d)\n", ii + 1, polynum, polynums[ii]);
	}
	int invalid_polynum = polynum;

	int num_bpoly = 1 << (dconf->curr_poly->s - 1);

	if (num_bpoly > dconf->ss_slices_n[0].numbuckets)
	{
		printf("not enough buckets allocated (%d) for current number of b-polys (%d)\n",
			dconf->ss_slices_n[0].numbuckets, num_bpoly);
		exit(1);
	}

	int need_bucket_alloc = !dconf->poly_buckets_allocated;
	for (i = 0; i < dconf->num_ss_slices; i++)
	{
		for (ii = 0; ii <= num_bpoly; ii++)
		{
			dconf->ss_slices_p[i].buckets[ii].size = 0;
			dconf->ss_slices_n[i].buckets[ii].size = 0;

			if (need_bucket_alloc)
			{
				dconf->ss_slices_p[i].buckets[ii].alloc = 65536;
				dconf->ss_slices_n[i].buckets[ii].alloc = 65536;
				dconf->ss_slices_p[i].buckets[ii].element = (uint32_t*)xmalloc(
					dconf->ss_slices_p[i].buckets[ii].alloc * sizeof(uint32_t));
				dconf->ss_slices_n[i].buckets[ii].element = (uint32_t*)xmalloc(
					dconf->ss_slices_n[i].buckets[ii].alloc * sizeof(uint32_t));
			}
		}

		dconf->ss_slices_p[i].curr_poly_idx = 0;
		dconf->ss_slices_n[i].curr_poly_idx = 0;
		dconf->ss_slices_p[i].curr_poly_num = 0;
		dconf->ss_slices_n[i].curr_poly_num = 0;
	}
	if (need_bucket_alloc)
	{
		dconf->poly_buckets_allocated = 1;
	}

	double t_enum_roots;
	double t_sort_roots;
	double t_match_roots;
	double t_sort_buckets;
	struct timeval start, stop;


	// todo
	// allocate buckets for polysum hits, one for each b-poly and per sieve side.
	// X stop bucket sieving at the ss_start_B index.
	// X set the ss_start_B index in static poly init.
	// need a bucket dumping routine somewhere in tdiv.  maybe tdiv_large.
	// testing!

	// optimizations:
	// can still use vector root updates when enumerating roots
	// for each side of the sum.  just need to be careful to
	// do root matching per prime in that case.
	// probably can do vector match search since we take an element
	// on one side then sum and compare to many elements on the other side.
	// could either do this many primes at a time or maybe many elements per
	// prime at a time.
	// how many "other side" bins we need to search for matches needs to
	// be more exact than hardcoded +/- 1 bins from the center bin.
	int* rootupdates = dconf->rootupdates;
	update_t update_data = dconf->update_data;
	uint32_t* modsqrt = sconf->modsqrt_array;

	int maxbin1 = 0;
	int maxbin2 = 0;
	int nummatch = 0;
	int matchp1 = 0;
	int matchm1 = 0;

	double avg_size1 = 0.0;
	double avg_size2 = 0.0;
	double avg_size3 = 0.0;

	int totalbins1 = 0;
	int totalbins2 = 0;
	int totalbins3 = 0;

	for (i = fb->ss_start_B; i < fb->B; i++)
	{
		int numB;
		int maxB = 1 << (dconf->curr_poly->s - 1);
		int r1 = update_data.firstroots1[i], r2 = update_data.firstroots2[i];
		uint32_t bound = sconf->factor_base->B;
		uint32_t med_B = sconf->factor_base->med_B;
		int polynum = 0;
		int slice = (int)((i - sconf->factor_base->ss_start_B) / 65536);
		int fboffset = dconf->ss_slices_p[slice].fboffset;
		int prime = fb->list->prime[i];
		int p = prime;
		int root1 = modsqrt[i];
		int root2 = prime - root1;
		uint8_t logp = fb->list->logprime[i];

		int amodp = (int)mpz_tdiv_ui(poly->mpz_poly_a, prime);
		int bmodp = (int)mpz_tdiv_ui(poly->mpz_poly_b, prime);
		int inv;

		//find a^-1 mod p = inv(a mod p) mod p
		if (sconf->knmod8 == 1)
		{
			inv = modinv_1(2 * amodp, prime);
		}
		else
		{
			inv = modinv_1(amodp, prime);
		}

		// create set1, an enumeration of (+/-t - b)(a)^-1 mod p
		// for b in the set [+/-B1 +/- B2 +/- ... Bs/2]
		int ii, v, sign, n = poly->s / 2;
		int tmp;

#ifdef SS_TIMING
		gettimeofday(&start, NULL);
#endif

		prime = fb->list->prime[i];
		root1 = modsqrt[i];
		root2 = prime - root1;

		// first poly b: sum of first n positive Bl
		bmodp = (int)mpz_tdiv_ui(polyb1, prime);

		// under these conditions polya is added to polyb.
		// here we account for that by adding amodp into bmodp.
		// only here in set1 enumeration.
		if (sconf->knmod8 == 1)
		{
			if (sconf->charcount == 1)
			{
				bmodp += amodp;
				if (bmodp >= prime)
				{
					bmodp -= prime;
				}
			}
		}

		//COMPUTE_FIRST_ROOTS;
		r1 = (int)root1 - bmodp;
		r2 = (int)root2 - bmodp;
		if (r1 < 0) r1 += prime;
		if (r2 < 0) r2 += prime;

		r1 = (int)((uint64_t)inv * (uint64_t)r1 % (uint64_t)prime);
		r2 = (int)((uint64_t)inv * (uint64_t)r2 % (uint64_t)prime);

		int* ptr;

		// enumerate set1 roots
		for (ii = 1, polynum = 0; ii < (1 << n); ii++) {
			// these roots go into the set
			ss_set1a.root[ii - 1] = r1;
			ss_set1b.root[ii - 1] = r2;
			ss_set1a.polynum[ii - 1] = polynum;
			ss_set1b.polynum[ii - 1] = polynum;

			// next polynum
			polynum = polynums[ii];
			sign = polysign[ii];
			v = polyv[ii];

			// next roots
			ptr = &rootupdates[(v - 1) * bound + i];
			if (sign > 0)
			{
				r1 = (int)r1 - *ptr;
				r2 = (int)r2 - *ptr;
				r1 = (r1 < 0) ? r1 + prime : r1;
				r2 = (r2 < 0) ? r2 + prime : r2;
			}
			else
			{
				r1 = (int)r1 + *ptr;
				r2 = (int)r2 + *ptr;
				r1 = (r1 >= prime) ? r1 - prime : r1;
				r2 = (r2 >= prime) ? r2 - prime : r2;
			}
		}
		// these roots go into the set
		ss_set1a.root[ii - 1] = r1;
		ss_set1a.polynum[ii - 1] = polynum;
		ss_set1b.root[ii - 1] = r2;
		ss_set1b.polynum[ii - 1] = polynum;

		int size1 = n;

		// create set2, an enumeration of (+/-t - b)(a)^-1 mod p
		// for b in the set [+/-Bs/2+1 +/- Bs/2+2 ... + Bs]
		n = (poly->s - 1) - n;

		bmodp = (int)mpz_tdiv_ui(polyb2, prime);

		// - (B1 + B2 + B3 + ... Bs-1)
		bmodp = prime - bmodp;

		r1 = (uint32_t)((uint64_t)inv * (uint64_t)bmodp % (uint64_t)prime);

		// enumerate set2 roots
		for (ii = 1, polynum = 0; ii < (1 << n); ii++) {
			// these roots go into the set
			ss_set2a.root[ii - 1] = r1;
			ss_set2a.polynum[ii - 1] = polynum;

			polynum = polynums[ii] << size1;
			sign = polysign[ii];
			v = polyv[ii];

			// next roots
			ptr = &rootupdates[(size1 + v - 1) * bound + i];
			if (sign > 0)
			{
				r1 = (int)r1 - *ptr;
				r1 = (r1 < 0) ? r1 + prime : r1;
			}
			else
			{
				r1 = (int)r1 + *ptr;
				r1 = (r1 >= prime) ? r1 - prime : r1;
			}
		}

		// these roots go into the set
		ss_set2a.root[ii - 1] = r1;
		ss_set2a.polynum[ii - 1] = polynum;

#ifdef SS_TIMING
		gettimeofday(&stop, NULL);
		t_enum_roots += ytools_difftime(&start, &stop);

		gettimeofday(&start, NULL);
#endif

		// now sort the sets into a moderate number of bins over the range 0:p
		numbins = p / (16 * interval) + 1;
		bindepth = 128;
		binsize = p / numbins + 1;

		// initialize bin sizes
		for (ii = 0; ii < numbins; ii++)
		{
			bins1a[ii].size = 0;
			bins1b[ii].size = 0;
			bins2[ii].size = 0;
		}

		// sort root1 into set 1a
		for (ii = 0; ii < (1 << ss_set1a.size); ii++)
		{
			int binnum = ss_set1a.root[ii] / binsize;
			if (binnum < numbins)
			{
				//printf("bin %d <-- set1a root %d (L polynum %d)\n", binnum, ss_set1a.root[ii], ss_set1a.polynum[ii]);
				bins1a[binnum].root[bins1a[binnum].size] = ss_set1a.root[ii];
				bins1a[binnum].polynum[bins1a[binnum].size] = ss_set1a.polynum[ii];
				bins1a[binnum].size++;
				if (bins1a[binnum].size > bindepth)
				{
					printf("bin overflow\n");
					exit(1);
				}
			}
			else
			{
				printf("element %d of set 1, root %d, invalid bin %d [of %d]\n",
					ii, ss_set1a.root[ii], binnum, numbins);
			}
		}

		// sort root2 into set 1b
		for (ii = 0; ii < (1 << ss_set1b.size); ii++)
		{
			int binnum = ss_set1b.root[ii] / binsize;
			if (binnum < numbins)
			{
				// printf("bin %d <-- set1b root %d (polynum 0)\n", binnum, ss_set1a.root[ii]);
				bins1b[binnum].root[bins1b[binnum].size] = ss_set1b.root[ii];
				bins1b[binnum].polynum[bins1b[binnum].size] = ss_set1b.polynum[ii];
				bins1b[binnum].size++;
				if (bins1b[binnum].size > bindepth)
				{
					printf("bin overflow\n");
					exit(1);
				}
			}
			else
			{
				printf("element %d of set 1, root %d, invalid bin %d [of %d]\n",
					ii, ss_set1b.root[ii], binnum, numbins);
			}
		}

		// sort set 2
		for (ii = 0; ii < (1 << ss_set2a.size); ii++)
		{
			int binnum = ss_set2a.root[ii] / binsize;
			if (binnum < numbins)
			{
				//printf("bin %d <-- set2a root %d (R polynum %d)\n", binnum, ss_set2a.root[ii], ss_set2a.polynum[ii]);
				bins2[binnum].root[bins2[binnum].size] = ss_set2a.root[ii];
				bins2[binnum].polynum[bins2[binnum].size] = ss_set2a.polynum[ii];
				bins2[binnum].size++;
				if (bins2[binnum].size > bindepth)
				{
					printf("bin overflow\n");
					exit(1);
				}
			}
			else
			{
				printf("element %d of set 2, root %d, invalid bin %d [of %d]\n",
					ii, ss_set2a.root[ii], binnum, numbins);
			}
		}


#ifdef SS_TIMING
		gettimeofday(&stop, NULL);
		t_sort_roots += ytools_difftime(&start, &stop);

		gettimeofday(&start, NULL);
#endif

		// commence matching
		uint32_t pid = (uint32_t)(i - fboffset) << 16;

		for (ii = 0; ii < numbins; ii++)
		{
			avg_size1 += bins1a[ii].size;
			avg_size2 += bins1b[ii].size;
			avg_size3 += bins2[ii].size;
		}
		totalbins1 += numbins;

		__m512i vp = _mm512_set1_epi32(prime);
		__m512i vi = _mm512_set1_epi32(interval);
		__m512i vz = _mm512_setzero_epi32();

		for (ii = 0; ii < numbins; ii++)
		{
			int x = binsize * ii + binsize / 2;
			int px = p - x;
			int b = px / binsize;

			int k;
			for (k = 0; k < bins2[b].size; k++)
			{
				__m512i vb2root = _mm512_set1_epi32(bins2[b].root[k]);
				__m512i vb2poly = _mm512_set1_epi32(bins2[b].polynum[k]);
				int bin2root = bins2[b].root[k];

#if 0
				//if (j = 0) //
				//for (j = 0; j < bins1a[ii].size; j += 16)
				j = 0;
				if (bins1a[ii].size > 4)
				{
					__mmask16 loadmask;
					
					if ((bins1a[ii].size - j) >= 16)
						loadmask = 0xffff;
					else
						loadmask = (1 << (bins1a[ii].size - j)) - 1;

					__m512i vb1root = _mm512_mask_loadu_epi32(vz, loadmask,
						&bins1a[ii].root[j]);

					__m512i vsum = _mm512_add_epi32(vb1root, vb2root);
					__mmask16 mpos = loadmask & _mm512_cmpge_epi32_mask(vsum, vp);
					__mmask16 mneg = loadmask & (~mpos);

					__m512i vdiffp = _mm512_mask_sub_epi32(vp, mpos, vsum, vp);
					__m512i vdiffn = _mm512_mask_sub_epi32(vp, mneg, vp, vsum);

					__mmask16 mhitp = _mm512_cmplt_epi32_mask(vdiffp, vi);
					__mmask16 mhitn = _mm512_cmplt_epi32_mask(vdiffn, vi);

					__m512i vb1poly = _mm512_mask_loadu_epi32(vz, loadmask,
						&bins1a[ii].polynum[j]);

					if (mhitp > 0)
					{
						uint32_t sum[16], poly[16];

						_mm512_mask_storeu_epi32(sum, mhitp, vdiffp);
						_mm512_mask_storeu_epi32(poly, mhitp, _mm512_add_epi32(vb1poly, vb2poly));

						while (mhitp > 0)
						{
							int pos = _trail_zcnt(mhitp);
							int idx = poly[pos];

							//printf("p-hit at bins1a[%d] with poly %d: (%d | %d)\n",
							//	pos, idx, pid, sum[pos]);

							uint32_t bsz = dconf->ss_slices_p[slice].buckets[idx].size;
							dconf->ss_slices_p[slice].buckets[idx].element[bsz] = (pid | sum[pos]);
							dconf->ss_slices_p[slice].buckets[idx].size++;

							if ((dconf->ss_slices_p[slice].buckets[idx].size) >=
								dconf->ss_slices_p[slice].buckets[idx].alloc)
							{
								dconf->ss_slices_p[slice].buckets[idx].alloc *= 2;
								dconf->ss_slices_p[slice].buckets[idx].element = (uint32_t*)xrealloc(
									dconf->ss_slices_p[slice].buckets[idx].element,
									dconf->ss_slices_p[slice].buckets[idx].alloc *
									sizeof(uint32_t));
							}

							mhitp = _reset_lsb(mhitp);
						}
					}
					
					if (mhitn > 0)
					{
						uint32_t sum[16], poly[16];

						_mm512_mask_storeu_epi32(sum, mhitn, vdiffn);
						_mm512_mask_storeu_epi32(poly, mhitn, _mm512_add_epi32(vb1poly, vb2poly));

						while (mhitn > 0)
						{
							int pos = _trail_zcnt(mhitn);
							int idx = poly[pos];

							//printf("n-hit at bins1a[%d] with poly %d: (%d | %d)\n",
							//	pos, idx, pid, sum[pos]);

							uint32_t bsz = dconf->ss_slices_n[slice].buckets[idx].size;

							dconf->ss_slices_n[slice].buckets[idx].element[bsz] = (pid | sum[pos]);
							dconf->ss_slices_n[slice].buckets[idx].size++;

							if (dconf->ss_slices_n[slice].buckets[idx].size >=
								dconf->ss_slices_n[slice].buckets[idx].alloc)
							{
								dconf->ss_slices_n[slice].buckets[idx].alloc *= 2;
								dconf->ss_slices_n[slice].buckets[idx].element = (uint32_t*)xrealloc(
									dconf->ss_slices_n[slice].buckets[idx].element,
									dconf->ss_slices_n[slice].buckets[idx].alloc *
									sizeof(uint32_t));
							}

							mhitn = _reset_lsb(mhitn);
						}
					}

					j += 16;
				}

				for (; j < bins1a[ii].size; j++)
				{
					int sum1 = bin2root + bins1a[ii].root[j];
					int sign1 = (sum1 >= prime);
					sum1 = sign1 ? sum1 - prime : prime - sum1;

					if (sum1 < interval)
					{
						int polysum = bins1a[ii].polynum[j] + bins2[b].polynum[k];
						int idx = polysum;

						if (sign1)
						{
							uint32_t bsz = dconf->ss_slices_p[slice].buckets[idx].size;

							dconf->ss_slices_p[slice].buckets[idx].element[bsz] = (pid | sum1);
							dconf->ss_slices_p[slice].buckets[idx].size++;

							if (dconf->ss_slices_p[slice].buckets[idx].size >=
								dconf->ss_slices_p[slice].buckets[idx].alloc)
							{
								dconf->ss_slices_p[slice].buckets[idx].alloc *= 2;
								dconf->ss_slices_p[slice].buckets[idx].element = (uint32_t*)xrealloc(
									dconf->ss_slices_p[slice].buckets[idx].element,
									dconf->ss_slices_p[slice].buckets[idx].alloc *
									sizeof(uint32_t));
							}
						}
						else
						{
							uint32_t bsz = dconf->ss_slices_n[slice].buckets[idx].size;

							dconf->ss_slices_n[slice].buckets[idx].element[bsz] = (pid | sum1);
							dconf->ss_slices_n[slice].buckets[idx].size++;

							if (dconf->ss_slices_n[slice].buckets[idx].size >=
								dconf->ss_slices_n[slice].buckets[idx].alloc)
							{
								dconf->ss_slices_n[slice].buckets[idx].alloc *= 2;
								dconf->ss_slices_n[slice].buckets[idx].element = (uint32_t*)xrealloc(
									dconf->ss_slices_n[slice].buckets[idx].element,
									dconf->ss_slices_n[slice].buckets[idx].alloc *
									sizeof(uint32_t));
							}
						}

						nummatch++;
					}

				}

				//for (j = 0; j < bins1b[ii].size; j += 16)
				j = 0;
				if (bins1b[ii].size > 4) 
				{
					__mmask16 loadmask;

					if ((bins1b[ii].size - j) >= 16)
						loadmask = 0xffff;
					else
						loadmask = (1 << (bins1b[ii].size - j)) - 1;

					__m512i vb1root = _mm512_mask_loadu_epi32(vz, loadmask,
						&bins1b[ii].root[j]);

					__m512i vsum = _mm512_add_epi32(vb1root, vb2root);
					__mmask16 vpos = loadmask & _mm512_cmpge_epi32_mask(vsum, vp);
					__mmask16 vneg = loadmask & (~vpos);

					__m512i vdiffp = _mm512_mask_sub_epi32(vp, vpos, vsum, vp);
					__m512i vdiffn = _mm512_mask_sub_epi32(vp, vneg, vp, vsum);

					__mmask16 vhitp = _mm512_cmplt_epi32_mask(vdiffp, vi);
					__mmask16 vhitn = _mm512_cmplt_epi32_mask(vdiffn, vi);

					__m512i vb1poly = _mm512_mask_loadu_epi32(vz, loadmask,
						&bins1b[ii].polynum[j]);

					if (vhitp > 0)
					{
						uint32_t sum[16], poly[16];

						_mm512_mask_storeu_epi32(sum, vhitp, vdiffp);
						_mm512_mask_storeu_epi32(poly, vhitp, _mm512_add_epi32(vb1poly, vb2poly));

						while (vhitp > 0)
						{
							int pos = _trail_zcnt(vhitp);
							int idx = poly[pos];

							uint32_t bsz = dconf->ss_slices_p[slice].buckets[idx].size;
							dconf->ss_slices_p[slice].buckets[idx].element[bsz] = (pid | sum[pos]);
							dconf->ss_slices_p[slice].buckets[idx].size++;

							if ((dconf->ss_slices_p[slice].buckets[idx].size) >=
								dconf->ss_slices_p[slice].buckets[idx].alloc)
							{
								dconf->ss_slices_p[slice].buckets[idx].alloc *= 2;
								dconf->ss_slices_p[slice].buckets[idx].element = (uint32_t*)xrealloc(
									dconf->ss_slices_p[slice].buckets[idx].element,
									dconf->ss_slices_p[slice].buckets[idx].alloc *
									sizeof(uint32_t));
							}

							vhitp = _reset_lsb(vhitp);
						}
					}

					if (vhitn > 0)
					{
						uint32_t sum[16], poly[16];

						_mm512_mask_storeu_epi32(sum, vhitn, vdiffn);
						_mm512_mask_storeu_epi32(poly, vhitn, _mm512_add_epi32(vb1poly, vb2poly));

						while (vhitn > 0)
						{
							int pos = _trail_zcnt(vhitn);
							int idx = poly[pos];

							uint32_t bsz = dconf->ss_slices_n[slice].buckets[idx].size;

							dconf->ss_slices_n[slice].buckets[idx].element[bsz] = (pid | sum[pos]);
							dconf->ss_slices_n[slice].buckets[idx].size++;

							if (dconf->ss_slices_n[slice].buckets[idx].size >=
								dconf->ss_slices_n[slice].buckets[idx].alloc)
							{
								dconf->ss_slices_n[slice].buckets[idx].alloc *= 2;
								dconf->ss_slices_n[slice].buckets[idx].element = (uint32_t*)xrealloc(
									dconf->ss_slices_n[slice].buckets[idx].element,
									dconf->ss_slices_n[slice].buckets[idx].alloc *
									sizeof(uint32_t));
							}

							vhitn = _reset_lsb(vhitn);
						}
					}

					j += 16;
				}

				for (; j < bins1b[ii].size; j++)
				{
					int sum1 = bin2root + bins1b[ii].root[j];
					int sign1 = (sum1 >= prime);
					sum1 = sign1 ? sum1 - prime : prime - sum1;

					if (sum1 < interval)
					{
						int polysum = bins1b[ii].polynum[j] + bins2[b].polynum[k];
						int idx = polysum;

						if (sign1)
						{
							uint32_t bsz = dconf->ss_slices_p[slice].buckets[idx].size;

							dconf->ss_slices_p[slice].buckets[idx].element[bsz] = (pid | sum1);
							dconf->ss_slices_p[slice].buckets[idx].size++;

							if (dconf->ss_slices_p[slice].buckets[idx].size >=
								dconf->ss_slices_p[slice].buckets[idx].alloc)
							{
								dconf->ss_slices_p[slice].buckets[idx].alloc *= 2;
								dconf->ss_slices_p[slice].buckets[idx].element = (uint32_t*)xrealloc(
									dconf->ss_slices_p[slice].buckets[idx].element,
									dconf->ss_slices_p[slice].buckets[idx].alloc *
									sizeof(uint32_t));
							}
						}
						else
						{
							uint32_t bsz = dconf->ss_slices_n[slice].buckets[idx].size;

							dconf->ss_slices_n[slice].buckets[idx].element[bsz] = (pid | sum1);
							dconf->ss_slices_n[slice].buckets[idx].size++;

							if (dconf->ss_slices_n[slice].buckets[idx].size >=
								dconf->ss_slices_n[slice].buckets[idx].alloc)
							{
								dconf->ss_slices_n[slice].buckets[idx].alloc *= 2;
								dconf->ss_slices_n[slice].buckets[idx].element = (uint32_t*)xrealloc(
									dconf->ss_slices_n[slice].buckets[idx].element,
									dconf->ss_slices_n[slice].buckets[idx].alloc *
									sizeof(uint32_t));
							}
						}

						nummatch++;
					}

				}

#else

				for (j = 0; j < MIN(bins1a[ii].size, bins1b[ii].size); j++)
				{
					int sum1 = bin2root + bins1a[ii].root[j];
					int sum2 = bin2root + bins1b[ii].root[j];
					int sign1 = (sum1 >= prime);
					int sign2 = (sum2 >= prime);
					sum1 = sign1 ? sum1 - prime : prime - sum1;
					sum2 = sign2 ? sum2 - prime : prime - sum2;

					if (sum1 < interval)
					{
						int polysum = bins1a[ii].polynum[j] + bins2[b].polynum[k];
						int idx = polysum;

						if (sign1)
						{
							uint32_t bsz = dconf->ss_slices_p[slice].buckets[idx].size;

							dconf->ss_slices_p[slice].buckets[idx].element[bsz] = (pid | sum1);
							dconf->ss_slices_p[slice].buckets[idx].size++;

							if (dconf->ss_slices_p[slice].buckets[idx].size >=
								dconf->ss_slices_p[slice].buckets[idx].alloc)
							{
								dconf->ss_slices_p[slice].buckets[idx].alloc *= 2;
								dconf->ss_slices_p[slice].buckets[idx].element = (uint32_t*)xrealloc(
									dconf->ss_slices_p[slice].buckets[idx].element,
									dconf->ss_slices_p[slice].buckets[idx].alloc *
									sizeof(uint32_t));
							}
						}
						else
						{
							uint32_t bsz = dconf->ss_slices_n[slice].buckets[idx].size;

							dconf->ss_slices_n[slice].buckets[idx].element[bsz] = (pid | sum1);
							dconf->ss_slices_n[slice].buckets[idx].size++;

							if (dconf->ss_slices_n[slice].buckets[idx].size >=
								dconf->ss_slices_n[slice].buckets[idx].alloc)
							{
								dconf->ss_slices_n[slice].buckets[idx].alloc *= 2;
								dconf->ss_slices_n[slice].buckets[idx].element = (uint32_t*)xrealloc(
									dconf->ss_slices_n[slice].buckets[idx].element,
									dconf->ss_slices_n[slice].buckets[idx].alloc *
									sizeof(uint32_t));
							}
						}

						nummatch++;
					}

					if (sum2 < interval)
					{
						int polysum = bins1b[ii].polynum[j] + bins2[b].polynum[k];
						int idx = polysum;

						if (sign2)
						{
							uint32_t bsz = dconf->ss_slices_p[slice].buckets[idx].size;

							dconf->ss_slices_p[slice].buckets[idx].element[bsz] = (pid | sum2);
							dconf->ss_slices_p[slice].buckets[idx].size++;

							if (dconf->ss_slices_p[slice].buckets[idx].size >=
								dconf->ss_slices_p[slice].buckets[idx].alloc)
							{
								dconf->ss_slices_p[slice].buckets[idx].alloc *= 2;
								dconf->ss_slices_p[slice].buckets[idx].element = (uint32_t*)xrealloc(
									dconf->ss_slices_p[slice].buckets[idx].element,
									dconf->ss_slices_p[slice].buckets[idx].alloc *
									sizeof(uint32_t));
							}
						}
						else
						{
							uint32_t bsz = dconf->ss_slices_n[slice].buckets[idx].size;

							dconf->ss_slices_n[slice].buckets[idx].element[bsz] = (pid | sum2);
							dconf->ss_slices_n[slice].buckets[idx].size++;

							if (dconf->ss_slices_n[slice].buckets[idx].size >=
								dconf->ss_slices_n[slice].buckets[idx].alloc)
							{
								dconf->ss_slices_n[slice].buckets[idx].alloc *= 2;
								dconf->ss_slices_n[slice].buckets[idx].element = (uint32_t*)xrealloc(
									dconf->ss_slices_n[slice].buckets[idx].element,
									dconf->ss_slices_n[slice].buckets[idx].alloc *
									sizeof(uint32_t));
							}
						}

						nummatch++;
					}
				}

				if (bins1a[ii].size > bins1b[ii].size)
				{
					for (; j < bins1a[ii].size; j++)
					{
						int sum1 = bin2root + bins1a[ii].root[j];
						int sign1 = (sum1 >= prime);
						sum1 = sign1 ? sum1 - prime : prime - sum1;

						if (sum1 < interval)
						{
							int polysum = bins1a[ii].polynum[j] + bins2[b].polynum[k];
							int idx = polysum;

							if (sign1)
							{
								uint32_t bsz = dconf->ss_slices_p[slice].buckets[idx].size;

								dconf->ss_slices_p[slice].buckets[idx].element[bsz] = (pid | sum1);
								dconf->ss_slices_p[slice].buckets[idx].size++;

								if (dconf->ss_slices_p[slice].buckets[idx].size >=
									dconf->ss_slices_p[slice].buckets[idx].alloc)
								{
									dconf->ss_slices_p[slice].buckets[idx].alloc *= 2;
									dconf->ss_slices_p[slice].buckets[idx].element = (uint32_t*)xrealloc(
										dconf->ss_slices_p[slice].buckets[idx].element,
										dconf->ss_slices_p[slice].buckets[idx].alloc *
										sizeof(uint32_t));
								}
							}
							else
							{
								uint32_t bsz = dconf->ss_slices_n[slice].buckets[idx].size;

								dconf->ss_slices_n[slice].buckets[idx].element[bsz] = (pid | sum1);
								dconf->ss_slices_n[slice].buckets[idx].size++;

								if (dconf->ss_slices_n[slice].buckets[idx].size >=
									dconf->ss_slices_n[slice].buckets[idx].alloc)
								{
									dconf->ss_slices_n[slice].buckets[idx].alloc *= 2;
									dconf->ss_slices_n[slice].buckets[idx].element = (uint32_t*)xrealloc(
										dconf->ss_slices_n[slice].buckets[idx].element,
										dconf->ss_slices_n[slice].buckets[idx].alloc *
										sizeof(uint32_t));
								}
							}

							nummatch++;
						}

					}
				}
				else if (bins1a[ii].size < bins1b[ii].size)
				{
					for (; j < bins1b[ii].size; j++)
					{
						int sum1 = bin2root + bins1b[ii].root[j];
						int sign1 = (sum1 >= prime);
						sum1 = sign1 ? sum1 - prime : prime - sum1;

						if (sum1 < interval)
						{
							int polysum = bins1b[ii].polynum[j] + bins2[b].polynum[k];
							int idx = polysum;

							if (sign1)
							{
								uint32_t bsz = dconf->ss_slices_p[slice].buckets[idx].size;

								dconf->ss_slices_p[slice].buckets[idx].element[bsz] = (pid | sum1);
								dconf->ss_slices_p[slice].buckets[idx].size++;

								if (dconf->ss_slices_p[slice].buckets[idx].size >=
									dconf->ss_slices_p[slice].buckets[idx].alloc)
								{
									dconf->ss_slices_p[slice].buckets[idx].alloc *= 2;
									dconf->ss_slices_p[slice].buckets[idx].element = (uint32_t*)xrealloc(
										dconf->ss_slices_p[slice].buckets[idx].element,
										dconf->ss_slices_p[slice].buckets[idx].alloc *
										sizeof(uint32_t));
								}
							}
							else
							{
								uint32_t bsz = dconf->ss_slices_n[slice].buckets[idx].size;

								dconf->ss_slices_n[slice].buckets[idx].element[bsz] = (pid | sum1);
								dconf->ss_slices_n[slice].buckets[idx].size++;

								if (dconf->ss_slices_n[slice].buckets[idx].size >=
									dconf->ss_slices_n[slice].buckets[idx].alloc)
								{
									dconf->ss_slices_n[slice].buckets[idx].alloc *= 2;
									dconf->ss_slices_n[slice].buckets[idx].element = (uint32_t*)xrealloc(
										dconf->ss_slices_n[slice].buckets[idx].element,
										dconf->ss_slices_n[slice].buckets[idx].alloc *
										sizeof(uint32_t));
								}
							}

							nummatch++;
						}

					}
				}
#endif

			}

#if 0
			if ((b + 1) < numbins)
			{
				for (k = 0; k < bins2[b + 1].size; k++)
				{
					int bin2root = bins2[b + 1].root[k];

					for (j = 0; j < MIN(bins1a[ii].size, bins1b[ii].size); j++)
					{
						int sum1 = bin2root + bins1a[ii].root[j];
						int sum2 = bin2root + bins1b[ii].root[j];
						int sign1 = (sum1 >= prime);
						int sign2 = (sum2 >= prime);
						sum1 = sign1 ? sum1 - prime : prime - sum1;
						sum2 = sign2 ? sum2 - prime : prime - sum2;

						if (sum1 < interval)
						{
							int polysum = bins1a[ii].polynum[j] + bins2[b + 1].polynum[k];
							int idx = polysum;

							if (sign1)
							{
								uint32_t bsz = dconf->ss_slices_p[slice].buckets[idx].size;

								dconf->ss_slices_p[slice].buckets[idx].element[bsz] = (pid | sum1);
								dconf->ss_slices_p[slice].buckets[idx].size++;

								if (dconf->ss_slices_p[slice].buckets[idx].size >=
									dconf->ss_slices_p[slice].buckets[idx].alloc)
								{
									dconf->ss_slices_p[slice].buckets[idx].alloc *= 2;
									dconf->ss_slices_p[slice].buckets[idx].element = (uint32_t*)xrealloc(
										dconf->ss_slices_p[slice].buckets[idx].element,
										dconf->ss_slices_p[slice].buckets[idx].alloc *
										sizeof(uint32_t));
								}
							}
							else
							{
								uint32_t bsz = dconf->ss_slices_n[slice].buckets[idx].size;

								dconf->ss_slices_n[slice].buckets[idx].element[bsz] = (pid | sum1);
								dconf->ss_slices_n[slice].buckets[idx].size++;

								if (dconf->ss_slices_n[slice].buckets[idx].size >=
									dconf->ss_slices_n[slice].buckets[idx].alloc)
								{
									dconf->ss_slices_n[slice].buckets[idx].alloc *= 2;
									dconf->ss_slices_n[slice].buckets[idx].element = (uint32_t*)xrealloc(
										dconf->ss_slices_n[slice].buckets[idx].element,
										dconf->ss_slices_n[slice].buckets[idx].alloc *
										sizeof(uint32_t));
								}
							}

							matchp1++;
						}

						if (sum2 < interval)
						{
							int polysum = bins1b[ii].polynum[j] + bins2[b + 1].polynum[k];
							int idx = polysum;

							if (sign2)
							{
								uint32_t bsz = dconf->ss_slices_p[slice].buckets[idx].size;

								dconf->ss_slices_p[slice].buckets[idx].element[bsz] = (pid | sum2);
								dconf->ss_slices_p[slice].buckets[idx].size++;

								if (dconf->ss_slices_p[slice].buckets[idx].size >=
									dconf->ss_slices_p[slice].buckets[idx].alloc)
								{
									dconf->ss_slices_p[slice].buckets[idx].alloc *= 2;
									dconf->ss_slices_p[slice].buckets[idx].element = (uint32_t*)xrealloc(
										dconf->ss_slices_p[slice].buckets[idx].element,
										dconf->ss_slices_p[slice].buckets[idx].alloc *
										sizeof(uint32_t));
								}
							}
							else
							{
								uint32_t bsz = dconf->ss_slices_n[slice].buckets[idx].size;

								dconf->ss_slices_n[slice].buckets[idx].element[bsz] = (pid | sum2);
								dconf->ss_slices_n[slice].buckets[idx].size++;

								if (dconf->ss_slices_n[slice].buckets[idx].size >=
									dconf->ss_slices_n[slice].buckets[idx].alloc)
								{
									dconf->ss_slices_n[slice].buckets[idx].alloc *= 2;
									dconf->ss_slices_n[slice].buckets[idx].element = (uint32_t*)xrealloc(
										dconf->ss_slices_n[slice].buckets[idx].element,
										dconf->ss_slices_n[slice].buckets[idx].alloc *
										sizeof(uint32_t));
								}
							}

							matchp1++;
						}
					}

					if (bins1a[ii].size > bins1b[ii].size)
					{
						for (; j < bins1a[ii].size; j++)
						{
							int sum1 = bin2root + bins1a[ii].root[j];
							int sign1 = (sum1 >= prime);
							sum1 = sign1 ? sum1 - prime : prime - sum1;

							if (sum1 < interval)
							{
								int polysum = bins1a[ii].polynum[j] + bins2[b + 1].polynum[k];
								int idx = polysum;

								if (sign1)
								{
									uint32_t bsz = dconf->ss_slices_p[slice].buckets[idx].size;

									dconf->ss_slices_p[slice].buckets[idx].element[bsz] = (pid | sum1);
									dconf->ss_slices_p[slice].buckets[idx].size++;

									if (dconf->ss_slices_p[slice].buckets[idx].size >=
										dconf->ss_slices_p[slice].buckets[idx].alloc)
									{
										dconf->ss_slices_p[slice].buckets[idx].alloc *= 2;
										dconf->ss_slices_p[slice].buckets[idx].element = (uint32_t*)xrealloc(
											dconf->ss_slices_p[slice].buckets[idx].element,
											dconf->ss_slices_p[slice].buckets[idx].alloc *
											sizeof(uint32_t));
									}
								}
								else
								{
									uint32_t bsz = dconf->ss_slices_n[slice].buckets[idx].size;

									dconf->ss_slices_n[slice].buckets[idx].element[bsz] = (pid | sum1);
									dconf->ss_slices_n[slice].buckets[idx].size++;

									if (dconf->ss_slices_n[slice].buckets[idx].size >=
										dconf->ss_slices_n[slice].buckets[idx].alloc)
									{
										dconf->ss_slices_n[slice].buckets[idx].alloc *= 2;
										dconf->ss_slices_n[slice].buckets[idx].element = (uint32_t*)xrealloc(
											dconf->ss_slices_n[slice].buckets[idx].element,
											dconf->ss_slices_n[slice].buckets[idx].alloc *
											sizeof(uint32_t));
									}
								}

								matchp1++;
							}

						}
					}
					else if (bins1a[ii].size < bins1b[ii].size)
					{
						for (; j < bins1b[ii].size; j++)
						{
							int sum1 = bin2root + bins1b[ii].root[j];
							int sign1 = (sum1 >= prime);
							sum1 = sign1 ? sum1 - prime : prime - sum1;

							if (sum1 < interval)
							{
								int polysum = bins1b[ii].polynum[j] + bins2[b + 1].polynum[k];
								int idx = polysum;

								if (sign1)
								{
									uint32_t bsz = dconf->ss_slices_p[slice].buckets[idx].size;

									dconf->ss_slices_p[slice].buckets[idx].element[bsz] = (pid | sum1);
									dconf->ss_slices_p[slice].buckets[idx].size++;

									if (dconf->ss_slices_p[slice].buckets[idx].size >=
										dconf->ss_slices_p[slice].buckets[idx].alloc)
									{
										dconf->ss_slices_p[slice].buckets[idx].alloc *= 2;
										dconf->ss_slices_p[slice].buckets[idx].element = (uint32_t*)xrealloc(
											dconf->ss_slices_p[slice].buckets[idx].element,
											dconf->ss_slices_p[slice].buckets[idx].alloc *
											sizeof(uint32_t));
									}
								}
								else
								{
									uint32_t bsz = dconf->ss_slices_n[slice].buckets[idx].size;

									dconf->ss_slices_n[slice].buckets[idx].element[bsz] = (pid | sum1);
									dconf->ss_slices_n[slice].buckets[idx].size++;

									if (dconf->ss_slices_n[slice].buckets[idx].size >=
										dconf->ss_slices_n[slice].buckets[idx].alloc)
									{
										dconf->ss_slices_n[slice].buckets[idx].alloc *= 2;
										dconf->ss_slices_n[slice].buckets[idx].element = (uint32_t*)xrealloc(
											dconf->ss_slices_n[slice].buckets[idx].element,
											dconf->ss_slices_n[slice].buckets[idx].alloc *
											sizeof(uint32_t));
									}
								}

								matchp1++;
							}

						}
					}
				}
			}

			if ((b - 1) >= 0)
			{
				for (k = 0; k < bins2[b - 1].size; k++)
				{
					int bin2root = bins2[b - 1].root[k];

					for (j = 0; j < MIN(bins1a[ii].size, bins1b[ii].size); j++)
					{
						int sum1 = bin2root + bins1a[ii].root[j];
						int sum2 = bin2root + bins1b[ii].root[j];
						int sign1 = (sum1 >= prime);
						int sign2 = (sum2 >= prime);
						sum1 = sign1 ? sum1 - prime : prime - sum1;
						sum2 = sign2 ? sum2 - prime : prime - sum2;

						if (sum1 < interval)
						{
							int polysum = bins1a[ii].polynum[j] + bins2[b - 1].polynum[k];
							int idx = polysum;

							if (sign1)
							{
								uint32_t bsz = dconf->ss_slices_p[slice].buckets[idx].size;

								dconf->ss_slices_p[slice].buckets[idx].element[bsz] = (pid | sum1);
								dconf->ss_slices_p[slice].buckets[idx].size++;

								if (dconf->ss_slices_p[slice].buckets[idx].size >=
									dconf->ss_slices_p[slice].buckets[idx].alloc)
								{
									dconf->ss_slices_p[slice].buckets[idx].alloc *= 2;
									dconf->ss_slices_p[slice].buckets[idx].element = (uint32_t*)xrealloc(
										dconf->ss_slices_p[slice].buckets[idx].element,
										dconf->ss_slices_p[slice].buckets[idx].alloc *
										sizeof(uint32_t));
								}
							}
							else
							{
								uint32_t bsz = dconf->ss_slices_n[slice].buckets[idx].size;

								dconf->ss_slices_n[slice].buckets[idx].element[bsz] = (pid | sum1);
								dconf->ss_slices_n[slice].buckets[idx].size++;

								if (dconf->ss_slices_n[slice].buckets[idx].size >=
									dconf->ss_slices_n[slice].buckets[idx].alloc)
								{
									dconf->ss_slices_n[slice].buckets[idx].alloc *= 2;
									dconf->ss_slices_n[slice].buckets[idx].element = (uint32_t*)xrealloc(
										dconf->ss_slices_n[slice].buckets[idx].element,
										dconf->ss_slices_n[slice].buckets[idx].alloc *
										sizeof(uint32_t));
								}
							}

							matchm1++;
						}

						if (sum2 < interval)
						{
							int polysum = bins1b[ii].polynum[j] + bins2[b - 1].polynum[k];
							int idx = polysum;

							if (sign2)
							{
								uint32_t bsz = dconf->ss_slices_p[slice].buckets[idx].size;

								dconf->ss_slices_p[slice].buckets[idx].element[bsz] = (pid | sum2);
								dconf->ss_slices_p[slice].buckets[idx].size++;

								if (dconf->ss_slices_p[slice].buckets[idx].size >=
									dconf->ss_slices_p[slice].buckets[idx].alloc)
								{
									dconf->ss_slices_p[slice].buckets[idx].alloc *= 2;
									dconf->ss_slices_p[slice].buckets[idx].element = (uint32_t*)xrealloc(
										dconf->ss_slices_p[slice].buckets[idx].element,
										dconf->ss_slices_p[slice].buckets[idx].alloc *
										sizeof(uint32_t));
								}
							}
							else
							{
								uint32_t bsz = dconf->ss_slices_n[slice].buckets[idx].size;

								dconf->ss_slices_n[slice].buckets[idx].element[bsz] = (pid | sum2);
								dconf->ss_slices_n[slice].buckets[idx].size++;

								if (dconf->ss_slices_n[slice].buckets[idx].size >=
									dconf->ss_slices_n[slice].buckets[idx].alloc)
								{
									dconf->ss_slices_n[slice].buckets[idx].alloc *= 2;
									dconf->ss_slices_n[slice].buckets[idx].element = (uint32_t*)xrealloc(
										dconf->ss_slices_n[slice].buckets[idx].element,
										dconf->ss_slices_n[slice].buckets[idx].alloc *
										sizeof(uint32_t));
								}
							}

							matchm1++;
						}
					}

					if (bins1a[ii].size > bins1b[ii].size)
					{
						for (; j < bins1a[ii].size; j++)
						{
							int sum1 = bin2root + bins1a[ii].root[j];
							int sign1 = (sum1 >= prime);
							sum1 = sign1 ? sum1 - prime : prime - sum1;

							if (sum1 < interval)
							{
								int polysum = bins1a[ii].polynum[j] + bins2[b - 1].polynum[k];
								int idx = polysum;

								if (sign1)
								{
									uint32_t bsz = dconf->ss_slices_p[slice].buckets[idx].size;

									dconf->ss_slices_p[slice].buckets[idx].element[bsz] = (pid | sum1);
									dconf->ss_slices_p[slice].buckets[idx].size++;

									if (dconf->ss_slices_p[slice].buckets[idx].size >=
										dconf->ss_slices_p[slice].buckets[idx].alloc)
									{
										dconf->ss_slices_p[slice].buckets[idx].alloc *= 2;
										dconf->ss_slices_p[slice].buckets[idx].element = (uint32_t*)xrealloc(
											dconf->ss_slices_p[slice].buckets[idx].element,
											dconf->ss_slices_p[slice].buckets[idx].alloc *
											sizeof(uint32_t));
									}
								}
								else
								{
									uint32_t bsz = dconf->ss_slices_n[slice].buckets[idx].size;

									dconf->ss_slices_n[slice].buckets[idx].element[bsz] = (pid | sum1);
									dconf->ss_slices_n[slice].buckets[idx].size++;

									if (dconf->ss_slices_n[slice].buckets[idx].size >=
										dconf->ss_slices_n[slice].buckets[idx].alloc)
									{
										dconf->ss_slices_n[slice].buckets[idx].alloc *= 2;
										dconf->ss_slices_n[slice].buckets[idx].element = (uint32_t*)xrealloc(
											dconf->ss_slices_n[slice].buckets[idx].element,
											dconf->ss_slices_n[slice].buckets[idx].alloc *
											sizeof(uint32_t));
									}
								}

								matchm1++;
							}

						}
					}
					else if (bins1a[ii].size < bins1b[ii].size)
					{
						for (; j < bins1b[ii].size; j++)
						{
							int sum1 = bin2root + bins1b[ii].root[j];
							int sign1 = (sum1 >= prime);
							sum1 = sign1 ? sum1 - prime : prime - sum1;

							if (sum1 < interval)
							{
								int polysum = bins1b[ii].polynum[j] + bins2[b - 1].polynum[k];
								int idx = polysum;

								if (sign1)
								{
									uint32_t bsz = dconf->ss_slices_p[slice].buckets[idx].size;

									dconf->ss_slices_p[slice].buckets[idx].element[bsz] = (pid | sum1);
									dconf->ss_slices_p[slice].buckets[idx].size++;

									if (dconf->ss_slices_p[slice].buckets[idx].size >=
										dconf->ss_slices_p[slice].buckets[idx].alloc)
									{
										dconf->ss_slices_p[slice].buckets[idx].alloc *= 2;
										dconf->ss_slices_p[slice].buckets[idx].element = (uint32_t*)xrealloc(
											dconf->ss_slices_p[slice].buckets[idx].element,
											dconf->ss_slices_p[slice].buckets[idx].alloc *
											sizeof(uint32_t));
									}
								}
								else
								{
									uint32_t bsz = dconf->ss_slices_n[slice].buckets[idx].size;

									dconf->ss_slices_n[slice].buckets[idx].element[bsz] = (pid | sum1);
									dconf->ss_slices_n[slice].buckets[idx].size++;

									if (dconf->ss_slices_n[slice].buckets[idx].size >=
										dconf->ss_slices_n[slice].buckets[idx].alloc)
									{
										dconf->ss_slices_n[slice].buckets[idx].alloc *= 2;
										dconf->ss_slices_n[slice].buckets[idx].element = (uint32_t*)xrealloc(
											dconf->ss_slices_n[slice].buckets[idx].element,
											dconf->ss_slices_n[slice].buckets[idx].alloc *
											sizeof(uint32_t));
									}
								}

								matchm1++;
							}

						}
					}
				}
			}
#endif

		}


#ifdef SS_TIMING
		gettimeofday(&stop, NULL);
		t_match_roots += ytools_difftime(&start, &stop);
#endif

	}


#ifdef SS_TIMING
	printf("found %d sieve hits matching bins\n", nummatch);
	printf("found %d sieve hits matching bins + 1\n", matchp1);
	printf("found %d sieve hits matching bins - 1\n", matchm1);
	printf("avg_size bins1 = %1.4f\n", avg_size1 / totalbins1);
	printf("avg_size bins2 = %1.4f\n", avg_size2 / totalbins1);
	printf("avg_size bins3 = %1.4f\n", avg_size3 / totalbins1);
	printf("enumerating roots: %1.4f seconds\n", t_enum_roots);
	printf("sorting roots: %1.4f seconds\n", t_sort_roots);
	printf("matching roots: %1.4f seconds\n", t_match_roots);
#endif

	//exit(0);

	mpz_clear(polyb1);
	mpz_clear(polyb2);
	free(polyv);
	free(polysign);
	free(polynums);

	for (ii = 0; ii < numbins; ii++)
	{
		free(bins1a[ii].root);
		free(bins1b[ii].root);
		free(bins2[ii].root);
		free(bins1a[ii].polynum);
		free(bins1b[ii].polynum);
		free(bins2[ii].polynum);
	}

	free(bins1a);
	free(bins1b);
	free(bins2);

	free(ss_set1a.root);
	free(ss_set1a.polynum);
	free(ss_set1b.root);
	free(ss_set1b.polynum);
	free(ss_set2a.root);
	free(ss_set2a.polynum);
	free(ss_set2b.root);
	free(ss_set2b.polynum);


#endif

	return;
}

#define P_PER_LOOP 16

void ss_search_poly_buckets_x16(static_conf_t* sconf, dynamic_conf_t* dconf)
{
	// the subset-sum search algorithm using a bucket sort to
	// organize the sieve-hits per poly.  
	// bucket sorting makes sieving very fast because it is just
	// a linear dump from a bucket into the sieve, and
	// SIMD can even help a little with that.
	// But sorting itself is very slow because there are *many* 
	// polynomials and they are widely spaced in memory, making 
	// the sort slow.
	// An idea is sort groups of polynomials so that we aren't 
	// accessing so many buckets during the sort.  Then sieving
	// can use similar bookkeeping to the sorted-lists approach
	// within each poly-group bucket.
	siqs_poly* poly = dconf->curr_poly;
	fb_list* fb = sconf->factor_base;
	uint32_t numblocks = sconf->num_blocks;
	uint32_t interval;

	numblocks = sconf->num_blocks;
	interval = numblocks << 15;

#ifdef USE_SS_SEARCH

	int i, j, ii;
	int ss_num_poly_terms = poly->s / 2;
	ss_set_t ss_set1a[16], ss_set1b[16], ss_set2a[16];

	mpz_t polyb1, polyb2;
	mpz_init(polyb1);
	mpz_init(polyb2);

#ifdef DEBUG_SS_SEARCH
	printf("allocating %u bytes for set1a/b of size %d\n",
		(1 << ss_num_poly_terms) * sizeof(int) * 2 * 2, ss_num_poly_terms);
#endif

	for (i = 0; i < 16; i++)
	{
		ss_set1a[i].root = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
		ss_set1a[i].polynum = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
		ss_set1a[i].size = ss_num_poly_terms;

		ss_set1b[i].root = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
		ss_set1b[i].polynum = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
		ss_set1b[i].size = ss_num_poly_terms;
	}

	//gmp_printf("%Zd\n", dconf->Bl[0]);
	mpz_set(polyb1, dconf->Bl[0]);

	for (ii = 1; ii < ss_num_poly_terms; ii++) {
		//gmp_printf("%Zd\n", dconf->Bl[ii]);
		mpz_add(polyb1, polyb1, dconf->Bl[ii]);
	}
	mpz_tdiv_q_2exp(polyb1, polyb1, 1);

	ss_num_poly_terms = (poly->s - 1) - ss_num_poly_terms;

#ifdef DEBUG_SS_SEARCH
	printf("allocating %u bytes for set2a/b of size %d\n",
		(1 << ss_num_poly_terms) * sizeof(int) * 2 * 2, ss_num_poly_terms);
#endif

	for (i = 0; i < 16; i++)
	{
		ss_set2a[i].root = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
		ss_set2a[i].polynum = (int*)xmalloc((1 << ss_num_poly_terms) * sizeof(int));
		ss_set2a[i].size = ss_num_poly_terms;
	}

	// first poly b: sum of first n positive Bl
	//gmp_printf("%Zd\n", dconf->Bl[size1]);
	mpz_set(polyb2, dconf->Bl[ii]);

	for (ii++; ii < poly->s; ii++) {
		//gmp_printf("%Zd\n", dconf->Bl[size1 + ii]);
		mpz_add(polyb2, polyb2, dconf->Bl[ii]);
	}
	mpz_tdiv_q_2exp(polyb2, polyb2, 1);

	int numbins = 2 * fb->list->prime[fb->B - 1] / (interval)+1;
	int bindepth = 128;
	int binsize; // = p / numbins + 1;

	ss_set_t** bins1;
	ss_set_t** bins2;

	// init bins
	bins1 = (ss_set_t**)xmalloc(16 * sizeof(ss_set_t*));
	bins2 = (ss_set_t**)xmalloc(16 * sizeof(ss_set_t*));

	for (i = 0; i < 16; i++)
	{
		bins1[i] = (ss_set_t*)xmalloc(numbins * sizeof(ss_set_t));
		bins2[i] = (ss_set_t*)xmalloc(numbins * sizeof(ss_set_t));

		for (ii = 0; ii < numbins; ii++)
		{
			bins1[i][ii].root = (int*)xmalloc(bindepth * sizeof(int));
			bins2[i][ii].root = (int*)xmalloc(bindepth * sizeof(int));
			bins1[i][ii].polynum = (int*)xmalloc(bindepth * sizeof(int));
			bins2[i][ii].polynum = (int*)xmalloc(bindepth * sizeof(int));
			bins1[i][ii].alloc = bindepth;
			bins2[i][ii].alloc = bindepth;
			bins1[i][ii].size = 0;
			bins2[i][ii].size = 0;
		}
	}

	//printf("allocated %d bins of depth %d\n", numbins, bindepth);
	//exit(0);

#ifdef DEBUG_SS_SEARCH
	mpz_add(dconf->gmptmp1, polyb2, polyb1);

	if (sconf->knmod8 == 1)
	{
		if (sconf->charcount == 1)
		{
			printf("including polya in polyb sum\n");
			mpz_add(dconf->gmptmp1, dconf->gmptmp1, poly->mpz_poly_a);
		}
	}

	if (mpz_cmp(dconf->gmptmp1, poly->mpz_poly_b) == 0)
	{
		prime = fb->list->prime[fb->B - 1];
		gmp_printf("polyb = %Zd + %Zd = %Zd\npolyb matches! (%u, %u, %u, %u)\n",
			polyb1, polyb2, dconf->gmptmp1,
			prime,
			mpz_tdiv_ui(polyb1, prime),
			mpz_tdiv_ui(polyb2, prime),
			mpz_tdiv_ui(poly->mpz_poly_a, prime));
	}
	else
	{
		gmp_printf("polyb = %Zd\npolyb doesn't match :(\n", dconf->gmptmp1);
	}
#endif

#if defined( DEBUG_SS_SEARCH) || defined(DEBUG_SS_SEARCH2)
	ss_set_t orighits;
	orighits.root = (int*)xmalloc(512 * sizeof(int));
	orighits.polynum = (int*)xmalloc(512 * sizeof(int));
#endif

	// create a map from poly number back to 
	// gray code enumeration order.
	int polynum = 0;
	int* polymap = dconf->polymap; // (int*)xmalloc((1 << (poly->s - 1)) * sizeof(int));
	int* polynums = (int*)xmalloc((1 << (poly->s - 1)) * sizeof(int));
	int* polyv = (int*)xmalloc((1 << (poly->s - 1)) * sizeof(int));
	int* polysign = (int*)xmalloc((1 << (poly->s - 1)) * sizeof(int));

	//polymap[0] = 1;
	polymap[1] = 0;
	polynums[0] = 0;

	//printf("%d <- %d (%d)\n", 1, 0, polynums[0]);

	for (ii = 1, polynum = 0; ii < (1 << (poly->s - 1)); ii++)
	{
		int v;
		int tmp;
		int sign;

		// next poly enumeration
		v = 1;
		j = ii;
		while ((j & 1) == 0)
			v++, j >>= 1;
		tmp = ii + (1 << v) - 1;
		tmp = (tmp >> v);
		polyv[ii] = v;
		if (tmp & 1)
			polysign[ii] = -1;
		else
			polysign[ii] = 1;

		// next polynum
		polynum ^= (1 << (v - 1));
		polynums[ii] = polynum;
		//polymap[polynum] = ii + 1;
		polymap[ii + 1] = polynum;

		//printf("%d <- %d (%d)\n", ii + 1, polynum, polynums[ii]);
	}

	//exit(0);
	int num_bpoly = 1 << (dconf->curr_poly->s - 1);

	if (num_bpoly > dconf->ss_slices_n[0].numbuckets)
	{
		printf("not enough buckets allocated (%d) for current number of b-polys (%d)\n",
			dconf->ss_slices_n[0].numbuckets, num_bpoly);
		exit(1);
	}

	int need_bucket_alloc = !dconf->poly_buckets_allocated;
	for (i = 0; i < dconf->num_ss_slices; i++)
	{
		for (ii = 0; ii <= num_bpoly; ii++)
		{
			dconf->ss_slices_p[i].buckets[ii].size = 0;
			dconf->ss_slices_n[i].buckets[ii].size = 0;

			if (need_bucket_alloc)
			{
				dconf->ss_slices_p[i].buckets[ii].alloc = 65536;
				dconf->ss_slices_n[i].buckets[ii].alloc = 65536;
				dconf->ss_slices_p[i].buckets[ii].element = (uint32_t*)xmalloc(
					dconf->ss_slices_p[i].buckets[ii].alloc * sizeof(uint32_t));
				dconf->ss_slices_n[i].buckets[ii].element = (uint32_t*)xmalloc(
					dconf->ss_slices_n[i].buckets[ii].alloc * sizeof(uint32_t));
			}
		}

		dconf->ss_slices_p[i].curr_poly_idx = 0;
		dconf->ss_slices_n[i].curr_poly_idx = 0;
		dconf->ss_slices_p[i].curr_poly_num = 0;
		dconf->ss_slices_n[i].curr_poly_num = 0;
	}
	if (need_bucket_alloc)
	{
		dconf->poly_buckets_allocated = 1;
	}

	double t_enum_roots;
	double t_sort_roots;
	double t_match_roots;
	double t_sort_buckets;
	struct timeval start, stop;


	// todo
	// allocate buckets for polysum hits, one for each b-poly and per sieve side.
	// X stop bucket sieving at the ss_start_B index.
	// X set the ss_start_B index in static poly init.
	// need a bucket dumping routine somewhere in tdiv.  maybe tdiv_large.
	// testing!

	// optimizations:
	// can still use vector root updates when enumerating roots
	// for each side of the sum.  just need to be careful to
	// do root matching per prime in that case.
	// probably can do vector match search since we take an element
	// on one side then sum and compare to many elements on the other side.
	// could either do this many primes at a time or maybe many elements per
	// prime at a time.
	// how many "other side" bins we need to search for matches needs to
	// be more exact than hardcoded +/- 1 bins from the center bin.
	int* rootupdates = dconf->rootupdates;
	update_t update_data = dconf->update_data;
	uint32_t* modsqrt = sconf->modsqrt_array;

	uint32_t loopbound;

	if (((fb->B - fb->ss_start_B) % P_PER_LOOP) != 0)
	{
		loopbound = fb->B - P_PER_LOOP;
	}
	else
	{
		loopbound = fb->B;
	}

	//printf("commencing subset-sum search from prime index %d to %d (%d)\n",
//		fb->ss_start_B, loopbound, fb->B);

	for (i = fb->ss_start_B; i < loopbound; i += P_PER_LOOP)
	{
		int ii, p, v, sign, n = poly->s / 2;
		int tmp;

		uint32_t bound = sconf->factor_base->B;
		int polynum = 0;
		int slice = (int)((i - sconf->factor_base->ss_start_B) / 65536);
		int fboffset = dconf->ss_slices_p[slice].fboffset;
		uint8_t logp = fb->list->logprime[i];

		int amodp[P_PER_LOOP];
		int bmodp[P_PER_LOOP];
		int inv[P_PER_LOOP];
		int prime[P_PER_LOOP];
		int r1[P_PER_LOOP];
		int r2[P_PER_LOOP];
		int root1[P_PER_LOOP];
		int root2[P_PER_LOOP];

		for (p = 0; p < P_PER_LOOP; p++)
		{
			prime[p] = fb->list->prime[i + p];
			amodp[p] = (int)mpz_tdiv_ui(poly->mpz_poly_a, prime[p]);

			//find a^-1 mod p = inv(a mod p) mod p
			if (sconf->knmod8 == 1)
			{
				inv[p] = modinv_1(2 * amodp[p], prime[p]);
			}
			else
			{
				inv[p] = modinv_1(amodp[p], prime[p]);
			}
		}

		// create set1, an enumeration of (+/-t - b)(a)^-1 mod p
		// for b in the set [+/-B1 +/- B2 +/- ... Bs/2]


#ifdef SS_TIMING
		gettimeofday(&start, NULL);
#endif

		for (p = 0; p < P_PER_LOOP; p++)
		{
			root1[p] = modsqrt[i + p];
			root2[p] = prime[p] - root1[p];
		}


#ifdef DEBUG_SS_SEARCH
			//printf("p,t1,t2 = %d,%d,%d\n", prime, root1, root2);
		printf("building set1 B\n");
#endif

		// first poly b: sum of first n positive Bl
		for (p = 0; p < P_PER_LOOP; p++)
		{
			bmodp[p] = (int)mpz_tdiv_ui(polyb1, prime[p]);

			// under these conditions polya is added to polyb.
			// here we account for that by adding amodp into bmodp.
			// only here in set1 enumeration.
			if (sconf->knmod8 == 1)
			{
				if (sconf->charcount == 1)
				{
					bmodp[p] += amodp[p];
					if (bmodp[p] >= prime[p])
					{
						bmodp[p] -= prime[p];
					}
				}
			}

			//COMPUTE_FIRST_ROOTS;
			r1[p] = (int)root1[p] - bmodp[p];
			r2[p] = (int)root2[p] - bmodp[p];
			if (r1[p] < 0) r1[p] += prime[p];
			if (r2[p] < 0) r2[p] += prime[p];

			r1[p] = (int)((uint64_t)inv[p] * (uint64_t)r1[p] % (uint64_t)prime[p]);
			r2[p] = (int)((uint64_t)inv[p] * (uint64_t)r2[p] % (uint64_t)prime[p]);
		}

		int* ptr;

#ifdef DEBUG_SS_SEARCH
		printf("enumerating %d elements of set 1 with %d polyb terms\n", (1 << n), n);
#endif

		// enumerate set1 roots
		for (ii = 1, polynum = 0; ii < (1 << n); ii++) {

			for (p = 0; p < P_PER_LOOP; p++)
			{
				// these roots go into the set
				ss_set1a[p].root[ii - 1] = r1[p];
				ss_set1b[p].root[ii - 1] = r2[p];
				ss_set1a[p].polynum[ii - 1] = polynum;
				ss_set1b[p].polynum[ii - 1] = polynum;
			}

			// next polynum
			polynum = polynums[ii];
			sign = polysign[ii];
			v = polyv[ii];

			for (p = 0; p < P_PER_LOOP; p++)
			{
				// next roots
				ptr = &rootupdates[(v - 1) * bound + i + p];
				if (sign > 0)
				{
					r1[p] = (int)r1[p] - *ptr;
					r2[p] = (int)r2[p] - *ptr;
					r1[p] = (r1[p] < 0) ? r1[p] + prime[p] : r1[p];
					r2[p] = (r2[p] < 0) ? r2[p] + prime[p] : r2[p];
				}
				else
				{
					r1[p] = (int)r1[p] + *ptr;
					r2[p] = (int)r2[p] + *ptr;
					r1[p] = (r1[p] >= prime[p]) ? r1[p] - prime[p] : r1[p];
					r2[p] = (r2[p] >= prime[p]) ? r2[p] - prime[p] : r2[p];
				}
			}
		}
		for (p = 0; p < P_PER_LOOP; p++)
		{
			// these roots go into the set
			ss_set1a[p].root[ii - 1] = r1[p];
			ss_set1b[p].root[ii - 1] = r2[p];
			ss_set1a[p].polynum[ii - 1] = polynum;
			ss_set1b[p].polynum[ii - 1] = polynum;
		}

		int size1 = n;

		// create set2, an enumeration of (+/-t - b)(a)^-1 mod p
		// for b in the set [+/-Bs/2+1 +/- Bs/2+2 ... + Bs]
		n = (poly->s - 1) - n;

#ifdef DEBUG_SS_SEARCH
		printf("building set2 B\n");
#endif

		for (p = 0; p < P_PER_LOOP; p++)
		{
			// - (B1 + B2 + B3 + ... Bs-1)
			bmodp[p] = (int)mpz_tdiv_ui(polyb2, prime[p]);
			bmodp[p] = prime[p] - bmodp[p];
			r1[p] = (uint32_t)((uint64_t)inv[p] * (uint64_t)bmodp[p] % (uint64_t)prime[p]);
		}

		// enumerate set2 roots
		for (ii = 1, polynum = 0; ii < (1 << n); ii++) {

			for (p = 0; p < P_PER_LOOP; p++)
			{
				// these roots go into the set
				ss_set2a[p].root[ii - 1] = r1[p];
				ss_set2a[p].polynum[ii - 1] = polynum;
			}

			polynum = polynums[ii] << size1;
			sign = polysign[ii];
			v = polyv[ii];

			for (p = 0; p < P_PER_LOOP; p++)
			{
				// next roots
				ptr = &rootupdates[(size1 + v - 1) * bound + i + p];
				if (sign > 0)
				{
					r1[p] = (int)r1[p] - *ptr;
					r1[p] = (r1[p] < 0) ? r1[p] + prime[p] : r1[p];
				}
				else
				{
					r1[p] = (int)r1[p] + *ptr;
					r1[p] = (r1[p] >= prime[p]) ? r1[p] - prime[p] : r1[p];
				}
			}
		}
		for (p = 0; p < P_PER_LOOP; p++)
		{
			// these roots go into the set
			ss_set2a[p].root[ii - 1] = r1[p];
			ss_set2a[p].polynum[ii - 1] = polynum;
		}

#ifdef SS_TIMING
		gettimeofday(&stop, NULL);
		t_enum_roots += ytools_difftime(&start, &stop);

		gettimeofday(&start, NULL);
#endif

#ifdef DEBUG_SS_SEARCH
		printf("sorting rootset1 of size %d into %d bins of size %d\n", (1 << ss_set1a.size), numbins, binsize);
#endif

		for (p = 0; p < P_PER_LOOP; p++)
		{
			for (ii = 0; ii < numbins; ii++)
			{
				bins1[p][ii].size = 0;
				bins2[p][ii].size = 0;
			}
		}

		// sort root 1 into set 1
		numbins = prime[0] / (1 * interval) + 1;
		bindepth = 128;
		binsize = prime[0] / numbins + 1;
		for (p = 0; p < P_PER_LOOP; p++)
		{
			// now sort the sets into a moderate number of bins over the range 0:p
			for (ii = 0; ii < (1 << ss_set1a[p].size); ii++)
			{
				int binnum = ss_set1a[p].root[ii] / binsize;
				if (binnum < numbins)
				{
					bins1[p][binnum].root[bins1[p][binnum].size] = ss_set1a[p].root[ii];
					bins1[p][binnum].polynum[bins1[p][binnum].size] = ss_set1a[p].polynum[ii];
					bins1[p][binnum].size++;
				}
			}
		}

#ifdef DEBUG_SS_SEARCH
		printf("sorting rootset2 of size %d into %d bins of size %d\n", (1 << ss_set2a.size), numbins, binsize);
#endif

		// sort set 2
		for (j = 0; j < P_PER_LOOP; j++)
		{
			// now sort the sets into a moderate number of bins over the range 0:p
			for (ii = 0; ii < (1 << ss_set2a[j].size); ii++)
			{
				int binnum = ss_set2a[j].root[ii] / binsize;
				if (binnum < numbins)
				{
					//printf("bin %d <-- set2a root %d (R polynum %d)\n", binnum, ss_set2a.root[ii], ss_set2a.polynum[ii]);
					bins2[j][binnum].root[bins2[j][binnum].size] = ss_set2a[j].root[ii];
					bins2[j][binnum].polynum[bins2[j][binnum].size] = ss_set2a[j].polynum[ii];
					bins2[j][binnum].size++;
				}
			}
		}

		// now match all numbers in bin x with the other set's bin that
		// contains p - x
#ifdef DEBUG_SS_SEARCH
		printf("matching root1\n");
#endif

#ifdef SS_TIMING
		gettimeofday(&stop, NULL);
		t_sort_roots += ytools_difftime(&start, &stop);

		gettimeofday(&start, NULL);
#endif

		//uint16_t loadmasks[17] = { 0, 0x1, 0x3, 0x7, 0xf,
		//	0x1f, 0x3f, 0x7f, 0xff,
		//	0x1ff, 0x3ff, 0x7ff, 0xfff,
		//	0x1fff, 0x3fff, 0x7fff, 0xffff };
		uint32_t offset_pids[16];
		for (p = 0; p < P_PER_LOOP; p++)
		{
			offset_pids[p] = (i + p - fboffset) << 16;
		}

		__m512i vprime = _mm512_loadu_epi32(prime);
		__m512i vinterval = _mm512_set1_epi32(interval);
		__m512i vpid = _mm512_loadu_epi32(offset_pids);

		// match set1 holding root1 to set2
		int nummatch = 0;
		int matchp1 = 0;
		int matchm1 = 0;

		//for (p = 0; p < P_PER_LOOP; p++)
		{
			uint32_t pid = (uint32_t)(i + p - fboffset) << 16;

			for (ii = 0; ii < numbins; ii++)
			{
				int x = binsize * ii + binsize / 2;
				int px = prime[0] - x;
				int b = px / binsize; // px / binsize;

				nummatch = 0;
				matchp1 = 0;
				matchm1 = 0;

				//printf("prime %d: processing %d bin roots/polys from bin1 %d\n",
				//	prime, bins1[ii].size, ii);

				// there seem to be too few elements in each bin,
				// on average, for SIMD to do much good.  Instead lets
				// try to process many primes at once.  We'll need to 
				// sort into 16 sets of bins and gather load them, but
				// the scatter store into poly-buckets stays the same.

				// all primes in this loop use the same number and size of bins, so
				// ii and b (the bin numbers) will be the same for each prime.
				// however the bins won't necessarily have the same size across
				// all primes because that depends on how each prime gets sorted
				// which depends on its roots.

				// so we create masks for:
				// whether bins1[ii] has any more elements.
				// whether bins2[b] has any more elements.
				int b1sizes[P_PER_LOOP];
				int b2sizes[P_PER_LOOP];
				int b1index[P_PER_LOOP];
				int b2index[P_PER_LOOP];
				uint32_t memoffsets1r[P_PER_LOOP];
				uint32_t memoffsets2r[P_PER_LOOP];
				uint32_t memoffsets1p[P_PER_LOOP];
				uint32_t memoffsets2p[P_PER_LOOP];

				for (p = 0; p < P_PER_LOOP; p++)
				{
					b1sizes[p] = bins1[p][ii].size;
					b2sizes[p] = bins2[p][b].size;

					memoffsets1r[p] = (uint32_t)(&bins1[p][ii].root[0] - &bins1[0][ii].root[0]);
					memoffsets2r[p] = (uint32_t)(&bins2[p][b].root[0] - &bins2[0][b].root[0]);
					memoffsets1p[p] = (uint32_t)(&bins1[p][ii].polynum[0] - &bins1[0][ii].polynum[0]);
					memoffsets2p[p] = (uint32_t)(&bins2[p][b].polynum[0] - &bins2[0][b].polynum[0]);
				}

				__m512i vb1index = _mm512_setzero_epi32();
				__m512i vb2index = _mm512_setzero_epi32();
				__m512i vb1size = _mm512_loadu_epi32(b1sizes);
				__m512i vb2size = _mm512_loadu_epi32(b2sizes);
				__m512i vones = _mm512_set1_epi32(1);

				__mmask16 loadb1 = _mm512_cmplt_epi32_mask(vb1index, vb1size);
				__mmask16 loadb2 = _mm512_cmplt_epi32_mask(vb2index, vb2size);

				//__m512i vr1index = _mm512_add_epi32(_mm512_loadu_epi32(memoffsets1r), vb1index);
				//__m512i vr2index = _mm512_add_epi32(_mm512_loadu_epi32(memoffsets2r), vb2index);
				__m512i vp1index = _mm512_add_epi32(_mm512_loadu_epi32(memoffsets1p), vb1index);
				__m512i vp2index = _mm512_add_epi32(_mm512_loadu_epi32(memoffsets2p), vb2index);

				while (loadb1 > 0)
				{
					__m512i vbin1root = _mm512_mask_i32gather_epi32(vones,
						loadb1, _mm512_add_epi32(_mm512_loadu_epi32(memoffsets1r), vb1index), 
						&bins1[0][ii].root[0], 4);

					while (loadb2 > 0)
					{
						// in here, at least one prime still has an entry in
						// bins1[ii] and bins2[b].  The same prime needs to have
						// an entry in both in order for there to be a hit.
						// If anything passes this last check, do the gather load.
						uint32_t m = loadb1 & loadb2;

						printf("loadmask1: %04x\n", loadb1);
						printf("loadmask2: %04x\n", loadb2);
						printf("combined : %04x: ", m);
						for (p = 0; p < P_PER_LOOP; p++)
						{
							if ((1 << p) & (loadb1 & loadb2))
							{
								printf("%d ", p);
							}
						}
						printf("\n");

						if (m > 0)
						{
							uint32_t tmp_e[16];
							uint32_t tmp_poly[16];

							__m512i vbin2root = _mm512_mask_i32gather_epi32(vones,
								loadb2, _mm512_add_epi32(_mm512_loadu_epi32(memoffsets2r), vb2index), 
								&bins2[0][b].root[0], 4);

							// defer loading the polys, we only need them if any of
							// these roots sum to be less than the interval.
							vbin1root = _mm512_mask_add_epi32(vbin1root, m, vbin1root, vbin2root);

							__mmask16 mpos = _mm512_mask_cmpge_epi32_mask(m, vbin1root, vprime);
							__mmask16 mneg = m & (~mpos);
							__m512i vdiffp = _mm512_mask_sub_epi32(vbin1root, mpos, vbin1root, vprime);
							__m512i vdiffn = _mm512_mask_sub_epi32(vbin1root, mneg, vprime, vbin1root);

							// compare to interval...
							mpos = _mm512_mask_cmplt_epi32_mask(mpos, vdiffp, vinterval);
							mneg = _mm512_mask_cmplt_epi32_mask(mneg, vdiffn, vinterval);
							
							// at least one of the primes produced hit
							if (mpos | mneg)
							{
								// load and sum the polys
								__m512i vbin1poly = _mm512_mask_i32gather_epi32(vones,
									m, vp1index, &bins1[0][ii].polynum[0], 4);
								__m512i vbin2poly = _mm512_mask_i32gather_epi32(vones,
									m, vp2index, &bins2[0][b].polynum[0], 4);

								vbin1poly = _mm512_mask_add_epi32(vbin1poly, m, vbin1poly, vbin2poly);
								_mm512_storeu_epi32(tmp_poly, vbin1poly);
							}

							if (mpos)
							{
								// form the sieve element
								vdiffp = _mm512_or_epi32(vdiffp, vpid);
								_mm512_storeu_epi32(tmp_e, vdiffp);

								while (mpos > 0)
								{
									int idx = _trail_zcnt(mpos);
									int pidx = tmp_poly[idx];
									printf("positive-side gather sum[%d] = %d, pidx = %d\n",
										idx, tmp_e[idx], tmp_poly[idx]);
									mpos = _reset_lsb(mpos);
								}
							}

							// at least one of the primes produced a negative-side hit
							if (mneg)
							{
								// form the sieve element
								vdiffn = _mm512_or_epi32(vdiffn, vpid);
								_mm512_storeu_epi32(tmp_e, vdiffn);

								while (mneg > 0)
								{
									int idx = _trail_zcnt(mneg);
									int pidx = tmp_poly[idx];
									printf("negative-side gather sum[%d] = %d, pidx = %d\n",
										idx, tmp_e[idx], tmp_poly[idx]);
									mneg = _reset_lsb(mneg);
								}
							}
						}

						vb2index = _mm512_add_epi32(vb2index, vones);
						loadb2 = _mm512_cmplt_epi32_mask(vb2index, vb2size);
					}
					vb1index = _mm512_add_epi32(vb1index, vones);
					loadb1 = _mm512_cmplt_epi32_mask(vb1index, vb1size);
				}

				printf("reference:\n");
				for (p = 0; p < P_PER_LOOP; p++)
				{
					for (j = 0; j < bins1[p][ii].size; j++)
					{
						int k;
						for (k = 0; k < bins2[p][b].size; k++)
						{
							//if ((bins1[p][ii].size > 0) && (bins2[p][b].size > 0))
							{
								int sum = bins1[p][ii].root[j] + bins2[p][b].root[k];
								int polysum = bins1[p][ii].polynum[j] + bins2[p][b].polynum[k];
								int idx = polysum;

								if ((sum >= prime[p]) && ((sum - prime[p]) < interval))
								{
									uint32_t sid = (uint32_t)(sum - prime[p]);
									printf("positive-side sum[%d] = %u, pidx = %d\n", p,
										(sid) | ((i + p - fboffset) << 16), polysum);
								}
								else if ((sum < prime[p]) && ((prime[p] - sum) < interval))
								{
									uint32_t sid = (uint32_t)(prime[p] - sum);
									printf("negative-side sum[%d] = %u, pidx = %d\n", p,
										(sid) | ((i + p - fboffset) << 16), polysum);
								}
							}
						}
					}
				}

				exit(1);

				for (j = 0; j < bins1[p][ii].size; j++)
				{
					int k;

#if 0
					printf("loading %d bin roots/polys from bin2 %d of size %d\n",
						MIN(bins2[b].size, 16), b, bins2[b].size);

					__mmask16 vloadmask = loadmasks[MIN(bins2[b].size, 16)];

					__m512i bin1r = _mm512_set1_epi32(bins1[ii].root[j]);
					__m512i bin2r = _mm512_mask_loadu_epi32(_mm512_setzero_epi32(),
						vloadmask, bins2[b].root);
					__m512i bin1p = _mm512_set1_epi32(bins1[ii].polynum[j]);
					__m512i bin2p = _mm512_mask_loadu_epi32(_mm512_setzero_epi32(),
						vloadmask, bins2[b].polynum);

					__m512i vsum = _mm512_add_epi32(bin1r, bin2r);
					__m512i vpolysum = _mm512_add_epi32(bin1p, bin2p);

					__mmask16 vcmp = _mm512_cmpge_epi32_mask(vsum, vprime);
					__m512i vdiffp = _mm512_mask_sub_epi32(vsum, vcmp, vsum, vprime);
					__m512i vdiffn = _mm512_mask_sub_epi32(vsum, ~vcmp, vprime, vsum);

					vdiffp = _mm512_or_epi32(vdiffp, vpid);
					vdiffn = _mm512_or_epi32(vdiffn, vpid);

					__mmask16 vhitp = vloadmask & _mm512_mask_cmplt_epi32_mask(vcmp, vdiffp, vinterval);
					__mmask16 vhitn = vloadmask & _mm512_mask_cmplt_epi32_mask(~vcmp, vdiffn, vinterval);

					uint32_t elements[16], vpidx[16];
					_mm512_storeu_epi32(elements, vdiffp);
					_mm512_storeu_epi32(vpidx, vpolysum);

					printf("sorting %d p-side hits to poly buckets\n", _mm_popcnt_u32(vhitp));

					while (vhitp > 0)
					{
						int idx = _trail_zcnt(vhitp);
						int pidx = vpidx[idx];

						printf("element %u of poly idx %d\n", elements[idx], pidx);

						uint32_t bsz = dconf->ss_slices_p[slice].buckets[pidx].size;

						dconf->ss_slices_p[slice].buckets[pidx].element[bsz] = elements[idx];
						dconf->ss_slices_p[slice].buckets[pidx].size++;

						if (dconf->ss_slices_p[slice].buckets[pidx].size >=
							dconf->ss_slices_p[slice].buckets[pidx].alloc)
						{
							dconf->ss_slices_p[slice].buckets[pidx].alloc *= 2;
							dconf->ss_slices_p[slice].buckets[pidx].element = (uint32_t*)xrealloc(
								dconf->ss_slices_p[slice].buckets[pidx].element,
								dconf->ss_slices_p[slice].buckets[pidx].alloc *
								sizeof(uint32_t));
						}

						vhitp = _reset_lsb(vhitp);
					}

					_mm512_storeu_epi32(elements, vdiffn);

					printf("sorting %d n-side hits to poly buckets\n", _mm_popcnt_u32(vhitn));

					while (vhitn > 0)
					{
						int idx = _trail_zcnt(vhitn);
						int pidx = vpidx[idx];

						printf("element %u of poly idx %d\n", elements[idx], pidx);

						uint32_t bsz = dconf->ss_slices_p[slice].buckets[pidx].size;

						dconf->ss_slices_n[slice].buckets[pidx].element[bsz] = elements[idx];
						dconf->ss_slices_n[slice].buckets[pidx].size++;

						if (dconf->ss_slices_n[slice].buckets[pidx].size >=
							dconf->ss_slices_n[slice].buckets[pidx].alloc)
						{
							dconf->ss_slices_n[slice].buckets[pidx].alloc *= 2;
							dconf->ss_slices_n[slice].buckets[pidx].element = (uint32_t*)xrealloc(
								dconf->ss_slices_n[slice].buckets[pidx].element,
								dconf->ss_slices_n[slice].buckets[pidx].alloc *
								sizeof(uint32_t));
				}

						vhitn = _reset_lsb(vhitn);
			}
#endif

					for (k = 0; k < bins2[p][b].size; k++)
					{
						int sum = bins1[p][ii].root[j] + bins2[p][b].root[k];
						int polysum = bins1[p][ii].polynum[j] + bins2[p][b].polynum[k];
						int idx = polysum;

						if ((sum >= prime[p]) && ((sum - prime[p]) < interval))
						{
							uint32_t sid = (uint32_t)(sum - prime[p]);
							uint32_t bsz = dconf->ss_slices_p[slice].buckets[idx].size;

							dconf->ss_slices_p[slice].buckets[idx].element[bsz] = (pid | sid);
							dconf->ss_slices_p[slice].buckets[idx].size++;

							if (dconf->ss_slices_p[slice].buckets[idx].size >=
								dconf->ss_slices_p[slice].buckets[idx].alloc)
							{
								dconf->ss_slices_p[slice].buckets[idx].alloc *= 2;
								dconf->ss_slices_p[slice].buckets[idx].element = (uint32_t*)xrealloc(
									dconf->ss_slices_p[slice].buckets[idx].element,
									dconf->ss_slices_p[slice].buckets[idx].alloc *
									sizeof(uint32_t));
							}

							nummatch++;

						}
						else if ((sum < prime[p]) && ((prime[p] - sum) < interval))
						{
							uint32_t sid = (uint32_t)(prime[p] - sum);
							uint32_t bsz = dconf->ss_slices_n[slice].buckets[idx].size;

							dconf->ss_slices_n[slice].buckets[idx].element[bsz] = (pid | sid);
							dconf->ss_slices_n[slice].buckets[idx].size++;

							if (dconf->ss_slices_n[slice].buckets[idx].size >=
								dconf->ss_slices_n[slice].buckets[idx].alloc)
							{
								dconf->ss_slices_n[slice].buckets[idx].alloc *= 2;
								dconf->ss_slices_n[slice].buckets[idx].element = (uint32_t*)xrealloc(
									dconf->ss_slices_n[slice].buckets[idx].element,
									dconf->ss_slices_n[slice].buckets[idx].alloc *
									sizeof(uint32_t));
							}

							nummatch++;
						}
					}


#if 1
					if ((b + 1) < numbins)
					{
						for (k = 0; k < bins2[p][b + 1].size; k++)
						{
							int sum = bins1[p][ii].root[j] + bins2[p][b + 1].root[k];
							int polysum = bins1[p][ii].polynum[j] + bins2[p][b + 1].polynum[k];
							int idx = polysum;

							if ((sum >= prime[p]) && ((sum - prime[p]) < interval))
							{
								uint32_t sid = (uint32_t)(sum - prime[p]);
								uint32_t bsz = dconf->ss_slices_p[slice].buckets[idx].size;

								dconf->ss_slices_p[slice].buckets[idx].element[bsz] = (pid | sid);
								dconf->ss_slices_p[slice].buckets[idx].size++;

								if (dconf->ss_slices_p[slice].buckets[idx].size >=
									dconf->ss_slices_p[slice].buckets[idx].alloc)
								{
									dconf->ss_slices_p[slice].buckets[idx].alloc *= 2;
									dconf->ss_slices_p[slice].buckets[idx].element = (uint32_t*)xrealloc(
										dconf->ss_slices_p[slice].buckets[idx].element,
										dconf->ss_slices_p[slice].buckets[idx].alloc *
										sizeof(uint32_t));
								}

								matchp1++;
							}
							else if ((sum < prime[p]) && ((prime[p] - sum) < interval))
							{
								uint32_t sid = (uint32_t)(prime[p] - sum);
								uint32_t bsz = dconf->ss_slices_n[slice].buckets[idx].size;

								dconf->ss_slices_n[slice].buckets[idx].element[bsz] = (pid | sid);
								dconf->ss_slices_n[slice].buckets[idx].size++;

								if (dconf->ss_slices_n[slice].buckets[idx].size >=
									dconf->ss_slices_n[slice].buckets[idx].alloc)
								{
									dconf->ss_slices_n[slice].buckets[idx].alloc *= 2;
									dconf->ss_slices_n[slice].buckets[idx].element = (uint32_t*)xrealloc(
										dconf->ss_slices_n[slice].buckets[idx].element,
										dconf->ss_slices_n[slice].buckets[idx].alloc *
										sizeof(uint32_t));
								}

								matchp1++;
							}
						}
					}
					if ((b - 1) >= 0)
					{
						for (k = 0; k < bins2[p][b - 1].size; k++)
						{
							int sum = bins1[p][ii].root[j] + bins2[p][b - 1].root[k];
							int polysum = bins1[p][ii].polynum[j] + bins2[p][b - 1].polynum[k];
							int idx = polysum;

							if ((sum >= prime[p]) && ((sum - prime[p]) < interval)) {
								uint32_t sid = (uint32_t)(sum - prime[p]);
								uint32_t bsz = dconf->ss_slices_p[slice].buckets[idx].size;

								dconf->ss_slices_p[slice].buckets[idx].element[bsz] = (pid | sid);
								dconf->ss_slices_p[slice].buckets[idx].size++;

								if (dconf->ss_slices_p[slice].buckets[idx].size >=
									dconf->ss_slices_p[slice].buckets[idx].alloc)
								{
									dconf->ss_slices_p[slice].buckets[idx].alloc *= 2;
									dconf->ss_slices_p[slice].buckets[idx].element = (uint32_t*)xrealloc(
										dconf->ss_slices_p[slice].buckets[idx].element,
										dconf->ss_slices_p[slice].buckets[idx].alloc *
										sizeof(uint32_t));
								}

								matchm1++;
							}
							else if ((sum < prime[p]) && ((prime[p] - sum) < interval))
							{
								uint32_t sid = (uint32_t)(prime[p] - sum);
								uint32_t bsz = dconf->ss_slices_n[slice].buckets[idx].size;

								dconf->ss_slices_n[slice].buckets[idx].element[bsz] = (pid | sid);
								dconf->ss_slices_n[slice].buckets[idx].size++;

								if (dconf->ss_slices_n[slice].buckets[idx].size >=
									dconf->ss_slices_n[slice].buckets[idx].alloc)
								{
									dconf->ss_slices_n[slice].buckets[idx].alloc *= 2;
									dconf->ss_slices_n[slice].buckets[idx].element = (uint32_t*)xrealloc(
										dconf->ss_slices_n[slice].buckets[idx].element,
										dconf->ss_slices_n[slice].buckets[idx].alloc *
										sizeof(uint32_t));
								}

								matchm1++;
							}
						}
					}
#endif

				}

				//printf("found %d sieve hits matching bin %d to bin %d\n", nummatch, ii, b);
				//printf("found %d sieve hits matching bin %d to bin %d\n", matchp1, ii, b + 1);
				//printf("found %d sieve hits matching bin %d to bin %d\n", matchm1, ii, b - 1);

			}
		}

		// reset set1 so we can put root2 into it
		for (p = 0; p < P_PER_LOOP; p++)
		{
			for (ii = 0; ii < numbins; ii++)
			{
				bins1[p][ii].size = 0;
			}
		}

#ifdef DEBUG_SS_SEARCH
		printf("sorting rootset1b of size %d into %d bins of size %d\n", (1 << ss_set1b.size), numbins, binsize);
#endif

		// sort root2 into set1
		for (j = 0; j < P_PER_LOOP; j++)
		{
			for (ii = 0; ii < (1 << ss_set1b[j].size); ii++)
			{
				int binnum = ss_set1b[j].root[ii] / binsize;
				if (binnum < numbins)
				{
					// printf("bin %d <-- set1b root %d (polynum 0)\n", binnum, ss_set1a.root[ii]);
					bins1[j][binnum].root[bins1[j][binnum].size] = ss_set1b[j].root[ii];
					bins1[j][binnum].polynum[bins1[j][binnum].size] = ss_set1b[j].polynum[ii];
					bins1[j][binnum].size++;
				}
			}
		}


#ifdef DEBUG_SS_SEARCH
		printf("sorting rootset2b of size %d into %d bins of size %d\n", (1 << ss_set2b.size), numbins, binsize);
#endif

		// now match all numbers in bin x with the other set's bin that
		// contains p - x
#ifdef DEBUG_SS_SEARCH
		printf("matching root2\n");
#endif

		// match set1 holding root2 to set2
		for (p = 0; p < P_PER_LOOP; p++)
		{
			uint32_t pid = (uint32_t)(i + p - fboffset) << 16;

			for (ii = 0; ii < numbins; ii++)
			{
				int x = binsize * ii + binsize / 2;
				int px = prime[0] - x;
				int b = px / binsize;

				nummatch = 0;
				matchp1 = 0;
				matchm1 = 0;

				{
					for (j = 0; j < bins1[p][ii].size; j++)
					{
						int k;
						for (k = 0; k < bins2[p][b].size; k++)
						{
							int sum = bins1[p][ii].root[j] + bins2[p][b].root[k];
							int polysum = bins1[p][ii].polynum[j] + bins2[p][b].polynum[k];
							int idx = polysum;

							if ((sum >= prime[p]) && ((sum - prime[p]) < interval))
							{
								uint32_t sid = (uint32_t)(sum - prime[p]);
								uint32_t bsz = dconf->ss_slices_p[slice].buckets[idx].size;

								dconf->ss_slices_p[slice].buckets[idx].element[bsz] = (pid | sid);
								dconf->ss_slices_p[slice].buckets[idx].size++;

								if (dconf->ss_slices_p[slice].buckets[idx].size >=
									dconf->ss_slices_p[slice].buckets[idx].alloc)
								{
									dconf->ss_slices_p[slice].buckets[idx].alloc *= 2;
									dconf->ss_slices_p[slice].buckets[idx].element = (uint32_t*)xrealloc(
										dconf->ss_slices_p[slice].buckets[idx].element,
										dconf->ss_slices_p[slice].buckets[idx].alloc *
										sizeof(uint32_t));
								}

								nummatch++;

							}
							else if ((sum < prime[p]) && ((prime[p] - sum) < interval))
							{
								uint32_t sid = (uint32_t)(prime[p] - sum);
								uint32_t bsz = dconf->ss_slices_n[slice].buckets[idx].size;

								dconf->ss_slices_n[slice].buckets[idx].element[bsz] = (pid | sid);
								dconf->ss_slices_n[slice].buckets[idx].size++;

								if (dconf->ss_slices_n[slice].buckets[idx].size >=
									dconf->ss_slices_n[slice].buckets[idx].alloc)
								{
									dconf->ss_slices_n[slice].buckets[idx].alloc *= 2;
									dconf->ss_slices_n[slice].buckets[idx].element = (uint32_t*)xrealloc(
										dconf->ss_slices_n[slice].buckets[idx].element,
										dconf->ss_slices_n[slice].buckets[idx].alloc *
										sizeof(uint32_t));
								}

								nummatch++;
							}
						}

#if 1
						if ((b + 1) < numbins)
						{
							for (k = 0; k < bins2[p][b + 1].size; k++)
							{
								int sum = bins1[p][ii].root[j] + bins2[p][b + 1].root[k];
								int polysum = bins1[p][ii].polynum[j] + bins2[p][b + 1].polynum[k];
								int idx = polysum;

								if ((sum >= prime[p]) && ((sum - prime[p]) < interval)) {
									uint32_t sid = (uint32_t)(sum - prime[p]);
									uint32_t bsz = dconf->ss_slices_p[slice].buckets[idx].size;

									dconf->ss_slices_p[slice].buckets[idx].element[bsz] = (pid | sid);
									dconf->ss_slices_p[slice].buckets[idx].size++;

									if (dconf->ss_slices_p[slice].buckets[idx].size >=
										dconf->ss_slices_p[slice].buckets[idx].alloc)
									{
										dconf->ss_slices_p[slice].buckets[idx].alloc *= 2;
										dconf->ss_slices_p[slice].buckets[idx].element = (uint32_t*)xrealloc(
											dconf->ss_slices_p[slice].buckets[idx].element,
											dconf->ss_slices_p[slice].buckets[idx].alloc *
											sizeof(uint32_t));
									}

									matchp1++;

								}
								else if ((sum < prime[p]) && ((prime[p] - sum) < interval))
								{
									uint32_t sid = (uint32_t)(prime[p] - sum);
									uint32_t bsz = dconf->ss_slices_n[slice].buckets[idx].size;

									dconf->ss_slices_n[slice].buckets[idx].element[bsz] = (pid | sid);
									dconf->ss_slices_n[slice].buckets[idx].size++;

									if (dconf->ss_slices_n[slice].buckets[idx].size >=
										dconf->ss_slices_n[slice].buckets[idx].alloc)
									{
										dconf->ss_slices_n[slice].buckets[idx].alloc *= 2;
										dconf->ss_slices_n[slice].buckets[idx].element = (uint32_t*)xrealloc(
											dconf->ss_slices_n[slice].buckets[idx].element,
											dconf->ss_slices_n[slice].buckets[idx].alloc *
											sizeof(uint32_t));
									}

									matchp1++;
								}
							}
						}
						if ((b - 1) >= 0)
						{
							for (k = 0; k < bins2[p][b - 1].size; k++)
							{
								int sum = bins1[p][ii].root[j] + bins2[p][b - 1].root[k];
								int polysum = bins1[p][ii].polynum[j] + bins2[p][b - 1].polynum[k];
								int idx = polysum;

								if ((sum >= prime[p]) && ((sum - prime[p]) < interval)) {
									uint32_t sid = (uint32_t)(sum - prime[p]);
									uint32_t bsz = dconf->ss_slices_p[slice].buckets[idx].size;

									dconf->ss_slices_p[slice].buckets[idx].element[bsz] = (pid | sid);
									dconf->ss_slices_p[slice].buckets[idx].size++;

									if (dconf->ss_slices_p[slice].buckets[idx].size >=
										dconf->ss_slices_p[slice].buckets[idx].alloc)
									{
										dconf->ss_slices_p[slice].buckets[idx].alloc *= 2;
										dconf->ss_slices_p[slice].buckets[idx].element = (uint32_t*)xrealloc(
											dconf->ss_slices_p[slice].buckets[idx].element,
											dconf->ss_slices_p[slice].buckets[idx].alloc *
											sizeof(uint32_t));
									}

									matchm1++;
								}
								else if ((sum < prime[p]) && ((prime[p] - sum) < interval))
								{
									uint32_t sid = (uint32_t)(prime[p] - sum);
									uint32_t bsz = dconf->ss_slices_n[slice].buckets[idx].size;

									dconf->ss_slices_n[slice].buckets[idx].element[bsz] = (pid | sid);
									dconf->ss_slices_n[slice].buckets[idx].size++;

									if (dconf->ss_slices_n[slice].buckets[idx].size >=
										dconf->ss_slices_n[slice].buckets[idx].alloc)
									{
										dconf->ss_slices_n[slice].buckets[idx].alloc *= 2;
										dconf->ss_slices_n[slice].buckets[idx].element = (uint32_t*)xrealloc(
											dconf->ss_slices_n[slice].buckets[idx].element,
											dconf->ss_slices_n[slice].buckets[idx].alloc *
											sizeof(uint32_t));
									}

									matchm1++;

								}
							}
						}
#endif

					}
				}

			//printf("found %d sieve hits matching bin %d to bin %d\n", nummatch, ii, b);
			//printf("found %d sieve hits matching bin %d to bin %d\n", matchp1, ii, b + 1);
			//printf("found %d sieve hits matching bin %d to bin %d\n", matchm1, ii, b - 1);
			}
		}

#ifdef SS_TIMING
		gettimeofday(&stop, NULL);
		t_match_roots += ytools_difftime(&start, &stop);
#endif

#if defined( DEBUG_SS_SEARCH) || defined(DEBUG_SS_SEARCH2)
		printf("found %d matches out of %d original hits for prime %d (%u)\n",
			nummatch, orighits.size, i, prime);

		//exit(1);
#endif

	}

	//printf("\nmax bin sizes for root1 = %d,%d\n", maxbin1, maxbin2);

#ifdef SS_TIMING
	printf("enumerating roots: %1.4f seconds\n", t_enum_roots);
	printf("sorting roots: %1.4f seconds\n", t_sort_roots);
	printf("matching roots: %1.4f seconds\n", t_match_roots);
	printf("sorting buckets: %1.4f seconds\n", t_sort_buckets);
#endif

	//exit(0);

	mpz_clear(polyb1);
	mpz_clear(polyb2);
	//free(polymap);
	free(polyv);
	free(polysign);
	free(polynums);

	for (j = 0; j < 16; j++)
	{
		for (ii = 0; ii < numbins; ii++)
		{
			free(bins1[j][ii].root);
			free(bins2[j][ii].root);
			free(bins1[j][ii].polynum);
			free(bins2[j][ii].polynum);
		}
	}

	free(bins1);
	free(bins2);

	for (j = 0; j < 16; j++)
	{
		free(ss_set1a[j].root);
		free(ss_set1a[j].polynum);
		free(ss_set1b[j].root);
		free(ss_set1b[j].polynum);
		free(ss_set2a[j].root);
		free(ss_set2a[j].polynum);
	}

#ifdef DEBUG_SS_SEARCH
	free(orighits.root);
	free(orighits.polynum);
#endif


#endif

	return;
}
#endif

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
	uint32_t start_prime = 2;
	int *rootupdates = dconf->rootupdates;
	update_t update_data = dconf->update_data;
	sieve_fb_compressed *fb_p = dconf->comp_sieve_p;
	sieve_fb_compressed *fb_n = dconf->comp_sieve_n;
	lp_bucket *lp_bucket_p = dconf->buckets;
	uint32_t *modsqrt = sconf->modsqrt_array;

	//locals
	uint32_t i, interval;
	uint8_t logp;
    int root1, root2, nroot1, nroot2, prime, amodp, bmodp, inv, x, bnum, j, numblocks;
	int s = poly->s;
	int bound_index = 0, k;
	uint32_t bound_val = fb->med_B;
	uint32_t *bptr, *sliceptr_p, *sliceptr_n;
    uint64_t* bptr64, * sliceptr64_p, * sliceptr64_n;
    uint8_t* slicelogp_ptr = NULL;
    uint32_t* slicebound_ptr = NULL;
	uint32_t *numptr_p, *numptr_n;
	int check_bound = BUCKET_ALLOC/2 - 1, room;
    FILE *out;
    uint32_t shift = 24;
    int u = 0;

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

        slicelogp_ptr = lp_bucket_p->logp;
        slicebound_ptr = lp_bucket_p->fb_bounds;
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
		uint64_t q64, tmp, t2;

		prime = fb->tinylist->prime[i];
		root1 = modsqrt[i]; 
		root2 = prime - root1; 

		amodp = (int)mpz_tdiv_ui(poly->mpz_poly_a,prime);
		bmodp = (int)mpz_tdiv_ui(poly->mpz_poly_b,prime);

		//find a^-1 mod p = inv(a mod p) mod p
        if (sconf->knmod8 == 1)
        {
            inv = modinv_1(2 * amodp, prime);
        }
        else
        {
            inv = modinv_1(amodp, prime);
        }

		COMPUTE_FIRST_ROOTS
	
		// reuse integer inverse of prime that we've calculated for use
		// in trial division stage
		// inv * root1 % prime
		t2 = (uint64_t)inv * (uint64_t)root1;
		tmp = t2 + (uint64_t)fb->tinylist->correction[i];
		q64 = tmp * (uint64_t)fb->tinylist->small_inv[i];
		tmp = q64 >> 32; 
		root1 = t2 - tmp * prime;

		// inv * root2 % prime
		t2 = (uint64_t)inv * (uint64_t)root2;
		tmp = t2 + (uint64_t)fb->tinylist->correction[i];
		q64 = tmp * (uint64_t)fb->tinylist->small_inv[i];
		tmp = q64 >> 32; 
		root2 = t2 - tmp * prime;
	
		//we don't sieve these primes, so ordering doesn't matter
		update_data.firstroots1[i] = root1;
		update_data.firstroots2[i] = root2;

		fb_p->root1[i] = (uint16_t)root1;
		fb_p->root2[i] = (uint16_t)root2;
		fb_n->root1[i] = (uint16_t)(prime - root2);
		fb_n->root2[i] = (uint16_t)(prime - root1);
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
			t2 = (uint64_t)inv * (uint64_t)x;
			tmp = t2 + (uint64_t)fb->tinylist->correction[i];
			q64 = tmp * (uint64_t)fb->tinylist->small_inv[i];
			tmp = q64 >> 32; 
			x = t2 - tmp * prime;

			rootupdates[(j)*fb->B+i] = x;
		}
	}

	for (i=sconf->sieve_small_fb_start;i<fb->fb_15bit_B;i++)
	{
		uint64_t tmp, t2;

		prime = fb->list->prime[i];
		root1 = modsqrt[i]; 
		root2 = prime - root1; 

		amodp = (int)mpz_tdiv_ui(poly->mpz_poly_a,prime);
		bmodp = (int)mpz_tdiv_ui(poly->mpz_poly_b,prime);

		//find a^-1 mod p = inv(a mod p) mod p
        if (sconf->knmod8 == 1)
        {
            inv = modinv_1(2 * amodp, prime);
        }
        else
        {
            inv = modinv_1(amodp, prime);
        }

		COMPUTE_FIRST_ROOTS

        // Barrett reduction for inv * root
        t2 = (uint64_t)inv * (uint64_t)root1;
        tmp = (t2 * fb->list->binv[i]) >> 32;
        t2 -= tmp * prime;
        root1 = (t2 >= prime) ? t2 - prime : t2;

        t2 = (uint64_t)inv * (uint64_t)root2;
        tmp = (t2 * fb->list->binv[i]) >> 32;
        t2 -= tmp * prime;
        root2 = (t2 >= prime) ? t2 - prime : t2;

		if (root2 < root1)
		{
			update_data.sm_firstroots1[i] = (uint16_t)root2;
			update_data.sm_firstroots2[i] = (uint16_t)root1;

			fb_p->root1[i] = (uint16_t)root2;
			fb_p->root2[i] = (uint16_t)root1;
			fb_n->root1[i] = (uint16_t)(prime - root1);
			fb_n->root2[i] = (uint16_t)(prime - root2);
		}
		else
		{
			update_data.sm_firstroots1[i] = (uint16_t)root1;
			update_data.sm_firstroots2[i] = (uint16_t)root2;

			fb_p->root1[i] = (uint16_t)root1;
			fb_p->root2[i] = (uint16_t)root2;
			fb_n->root1[i] = (uint16_t)(prime - root2);
			fb_n->root2[i] = (uint16_t)(prime - root1);
		}

		//for this factor base prime, compute the rootupdate value for all s
		//Bl values.  amodp holds a^-1 mod p
		//the rootupdate value is given by 2*Bj*amodp
		//Bl[j] now holds 2*Bl
		for (j=0;j<s;j++)
		{
			x = (int)mpz_tdiv_ui(dconf->Bl[j],prime);

			// x * inv % prime
            t2 = (uint64_t)inv * (uint64_t)x;
            tmp = (t2 * fb->list->binv[i]) >> 32;
            t2 -= tmp * prime;
            x = (t2 >= prime) ? t2 - prime : t2;

			rootupdates[(j)*fb->B+i] = x;
			dconf->sm_rootupdates[(j)*fb->med_B+i] = (uint16_t)x;
		}
	}

	for (i=fb->fb_15bit_B;i<fb->med_B;i++)
	{
		uint64_t tmp, t2;

		prime = fb->list->prime[i];
		root1 = modsqrt[i]; 
		root2 = prime - root1; 

		amodp = (int)mpz_tdiv_ui(poly->mpz_poly_a,prime);
		bmodp = (int)mpz_tdiv_ui(poly->mpz_poly_b,prime);

		//find a^-1 mod p = inv(a mod p) mod p
        if (sconf->knmod8 == 1)
        {
            inv = modinv_1(2 * amodp, prime);
        }
        else
        {
            inv = modinv_1(amodp, prime);
        }

		COMPUTE_FIRST_ROOTS

        // Barrett reduction for inv * root
        t2 = (uint64_t)inv * (uint64_t)root1;
        tmp = (t2 * fb->list->binv[i]) >> 32;
        t2 -= tmp * prime;
        root1 = (t2 >= prime) ? t2 - prime : t2;

        t2 = (uint64_t)inv * (uint64_t)root2;
        tmp = (t2 * fb->list->binv[i]) >> 32;
        t2 -= tmp * prime;
        root2 = (t2 >= prime) ? t2 - prime : t2;

		if (root2 < root1)
		{
			update_data.sm_firstroots1[i] = (uint16_t)root2;
			update_data.sm_firstroots2[i] = (uint16_t)root1;

			fb_p->root1[i] = (uint16_t)root2;
			fb_p->root2[i] = (uint16_t)root1;
			fb_n->root1[i] = (uint16_t)(prime - root1);
			fb_n->root2[i] = (uint16_t)(prime - root2);
		}
		else
		{
			update_data.sm_firstroots1[i] = (uint16_t)root1;
			update_data.sm_firstroots2[i] = (uint16_t)root2;

			fb_p->root1[i] = (uint16_t)root1;
			fb_p->root2[i] = (uint16_t)root2;
			fb_n->root1[i] = (uint16_t)(prime - root2);
			fb_n->root2[i] = (uint16_t)(prime - root1);
		}

		//for this factor base prime, compute the rootupdate value for all s
		//Bl values.  amodp holds a^-1 mod p
		//the rootupdate value is given by 2*Bj*amodp
		//Bl[j] now holds 2*Bl
		for (j=0;j<s;j++)
		{
			x = (int)mpz_tdiv_ui(dconf->Bl[j],prime);

			// x * inv % prime
            t2 = (uint64_t)inv * (uint64_t)x;
            tmp = (t2 * fb->list->binv[i]) >> 32;
            t2 -= tmp * prime;
            x = (t2 >= prime) ? t2 - prime : t2;

			rootupdates[(j)*fb->B+i] = x;
			dconf->sm_rootupdates[(j)*fb->med_B+i] = (uint16_t)x;
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
                //uint32_t ii;
                //uint32_t bb;
                //for (bb = 0; bb < numblocks; bb++)
                //{
                //    printf("%u lp p-roots in slice %d, block %d:\n", numptr_p[bb], bound_index, bb);
                //    for (ii = 0; ii < numptr_p[bb]; ii++)
                //    {
                //        //bptr = sliceptr_p + ((uint64_t)bb << BLOCKBITS) + (uint64_t)ii;
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
                //        //bptr = sliceptr_n + ((uint64_t)bb << BLOCKBITS) + (uint64_t)ii;
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
        if (sconf->knmod8 == 1)
        {
            inv = modinv_1(2 * amodp, prime);
        }
        else
        {
            inv = modinv_1(amodp, prime);
        }

        COMPUTE_FIRST_ROOTS

        root1 = (uint32_t)((uint64_t)inv * (uint64_t)root1 % (uint64_t)prime);
        root2 = (uint32_t)((uint64_t)inv * (uint64_t)root2 % (uint64_t)prime);
		
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
			x = (int)((int64_t)x * (int64_t)inv % (int64_t)prime);
            rootupdates[(j)*fb->B + i] = x;
		}
	}
    
	logp = fb->list->logprime[fb->large_B-1];

    for (i = fb->large_B; i < fb->B; i++)
    {
        CHECK_NEW_SLICE(i);

        prime = fb->list->prime[i];
        root1 = modsqrt[i];
        root2 = prime - root1;
		uint8_t logp = fb->list->logprime[i];

        amodp = (int)mpz_tdiv_ui(poly->mpz_poly_a, prime);
        bmodp = (int)mpz_tdiv_ui(poly->mpz_poly_b, prime);

        //find a^-1 mod p = inv(a mod p) mod p
        if (sconf->knmod8 == 1)
        {
            inv = modinv_1(2 * amodp, prime);
        }
        else
        {
            inv = modinv_1(amodp, prime);
        }

        COMPUTE_FIRST_ROOTS;

        root1 = (uint32_t)((uint64_t)inv * (uint64_t)root1 % (uint64_t)prime);
        root2 = (uint32_t)((uint64_t)inv * (uint64_t)root2 % (uint64_t)prime);

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
        for (j = 0; j < s; j++)
        {
            x = (int)mpz_tdiv_ui(dconf->Bl[j], prime);
            x = (int)((int64_t)x * (int64_t)inv % (int64_t)prime);
            rootupdates[(j)*fb->B + i] = x;
        }
    }


#ifdef USE_SS_SEARCH
	//ss_search_linked_lists(sconf, dconf);
	//ss_search_sorted_lists(sconf, dconf);
	ss_search_poly_buckets(sconf, dconf);
	//ss_search_poly_buckets_x16(sconf, dconf);
#endif

    if (lp_bucket_p->list != NULL)
    {
#ifdef USE_BATCHPOLY_X2
        int pnum;
        pnum = (dconf->numB % dconf->poly_batchsize) - 1;
        if (pnum < 0)
            pnum += dconf->poly_batchsize;

        lp_bucket_p->num_slices_batch[pnum] = bound_index + 1;

#ifdef DEBUGPRINT_BATCHPOLY
        printf("num slices after first poly (%d,%d) bucket sieve: %u\n", 
			dconf->numB, pnum, bound_index + 1);
#endif

#else
        lp_bucket_p->num_slices = bound_index + 1;
#endif
    }

	return;
}

