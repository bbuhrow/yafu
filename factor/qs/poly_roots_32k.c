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

#define SS_TIMING

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
	int slicesz = sconf->factor_base->slice_size;

	numblocks = sconf->num_blocks;
	interval = numblocks << 15;

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
	int bindepth = 256;
	int binsize;

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

	// create a map from poly number back to 
	// gray code enumeration order.
	int polynum = 0;
	int* polymap = dconf->polymap;
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
		for (ii = 0; ii <= (2 * num_bpoly); ii++)
		{
			dconf->ss_slices_p[i].size[ii] = 0;
			dconf->ss_slices_n[i].size[ii] = 0;
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
		int slice = (int)((i - sconf->factor_base->ss_start_B) / slicesz);
		int fboffset = dconf->ss_slices_p[slice].fboffset;
		int prime = fb->list->prime[i];
		int p = prime;
		int root1 = modsqrt[i];
		int root2 = prime - root1;
		uint8_t logp = fb->list->logprime[i];

		int amodp = (int)mpz_tdiv_ui(poly->mpz_poly_a, prime);
		int bmodp = (int)mpz_tdiv_ui(poly->mpz_poly_b, prime);
		int inv;

		if (slice >= sconf->factor_base->num_ss_slices)
		{
			printf("error slice number %d too large, %d allocated, at fb index %d (prime %u)\n",
				slice, sconf->factor_base->num_ss_slices, i, prime);
			exit(0);
		}


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
		numbins = p / (2 * interval) + 1;
		binsize = p / numbins + 1;

		// initialize bin sizes
		for (ii = 0; ii < numbins; ii++)
		{
			bins1[ii].size = 0;
			bins2[ii].size = 0;
		}

		// sort root1 into set 1
		for (ii = 0; ii < (1 << ss_set1a.size); ii++)
		{
			int binnum = ss_set1a.root[ii] / binsize;
			if (binnum < numbins)
			{
				//printf("bin %d <-- set1a root %d (L polynum %d)\n", binnum, ss_set1a.root[ii], ss_set1a.polynum[ii]);
				bins1[binnum].root[bins1[binnum].size] = ss_set1a.root[ii];
				bins1[binnum].polynum[bins1[binnum].size] = ss_set1a.polynum[ii];
				bins1[binnum].size++;
				if (bins1[binnum].size >= bindepth)
				{
					printf("\nbin overflow\n\n");
					exit(1);
				}
			}
			else
			{
				printf("element %d of set 1, root %d, invalid bin %d [of %d]\n",
					ii, ss_set1a.root[ii], binnum, numbins);
			}
		}

		// sort root2 into set 1
		for (ii = 0; ii < (1 << ss_set1b.size); ii++)
		{
			int binnum = ss_set1b.root[ii] / binsize;
			if (binnum < numbins)
			{
				// printf("bin %d <-- set1b root %d (polynum 0)\n", binnum, ss_set1a.root[ii]);
				bins1[binnum].root[bins1[binnum].size] = ss_set1b.root[ii];
				bins1[binnum].polynum[bins1[binnum].size] = ss_set1b.polynum[ii];
				bins1[binnum].size++;
				if (bins1[binnum].size >= bindepth)
				{
					printf("\nbin overflow\n\n");
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
		uint32_t pid = (uint32_t)(i - fboffset) << 18;

		if ((i - fboffset) >= slicesz)
		{
			printf("pid %d too large for slice %d with fboffset %d\n",
				pid, slice, fboffset);
			exit(1);
		}

		for (ii = 0; ii < numbins; ii++)
		{
			avg_size1 += bins1[ii].size;
			avg_size3 += bins2[ii].size;
		}
		totalbins1 += numbins;

		__m512i vp = _mm512_set1_epi32(prime);
		__m512i vi = _mm512_set1_epi32(interval);
		__m512i vz = _mm512_setzero_epi32();
		__m512i vpid = _mm512_set1_epi32(pid);

		uint32_t bucketalloc = dconf->ss_slices_p[0].alloc;
		uint32_t* pslice_ptr = dconf->ss_slices_p[slice].elements;
		uint32_t* nslice_ptr = dconf->ss_slices_n[slice].elements;
		uint32_t* psize_ptr = dconf->ss_slices_p[slice].size;
		uint32_t* nsize_ptr = dconf->ss_slices_n[slice].size;

		// ideas to try:
		// * put p-side and n-side hits into the same elements array so that 
		//   all hits can be gather/scattered from one base address.
		// * sort the bins by poly prior to matching?  bins are relatively sparse
		//   so sorting them is hopefully fast.  Bin2 holds the upper bits
		//   of the poly and bin1 hold the lower bits, so when we match an element
		//   from bin2 to many from bin1, we create matches for many poly buckets
		//   that are near each other.  Buckets are still fairly big though, so
		//   even consecutive polys are not super close.  Still, may be helpful.
		// * process a couple b-bins at a time, the independent operations could 
		//   hide some latency.
		// * don't separate p and n side buckets, put all hits, p or n, into the
		//   same poly bucket.  sort p vs. n during sieving.  half the buckets
		//   to write to here trading off hopefully very slight hit to sieve time.

		printf("prime = %u, binsize = %d, interval = %d, num bpoly = %d\n", prime, binsize, interval, num_bpoly);
		for (ii = 0; ii < numbins; ii++)
		{
			int x = binsize * ii + binsize / 2;
			int px = p - x;
			int b = px / binsize;

			if (((b+1) < numbins) && ((b-1) >= 0))
				printf("pairing bin %d with bin %d (and %d,%d) of %d total bins.  sizes: %d,%d,%d,%d\n",
					ii, b, b + 1, b - 1, numbins, bins2[b].size, bins2[b+1].size, bins2[b-1].size, bins1[ii].size);
			else if ((b + 1) >= numbins)
				printf("pairing bin %d with bin %d (and %d) of %d total bins.  sizes: %d,%d,%d\n",
					ii, b, b - 1, numbins, bins2[b].size, bins2[b - 1].size, bins1[ii].size);
			else if ((b - 1) <= 0)
				printf("pairing bin %d with bin %d (and %d) of %d total bins.  sizes: %d,%d,%d\n",
					ii, b, b + 1, numbins, bins2[b].size, bins2[b + 1].size, bins1[ii].size);

			int k;
			for (k = 0; k < bins2[b].size; k++)
			{
				__m512i vb2root = _mm512_set1_epi32(bins2[b].root[k]);
				__m512i vb2poly = _mm512_set1_epi32(bins2[b].polynum[k]);
				int bin2root = bins2[b].root[k];

#if 1
				//if (j = 0) //
				//for (j = 0; j < bins1[ii].size; j += 16)
				j = 0;
				if (bins1[ii].size > 4)
				{
					__mmask16 loadmask;
					
					if ((bins1[ii].size - j) >= 16)
						loadmask = 0xffff;
					else
						loadmask = (1 << (bins1[ii].size - j)) - 1;

					__m512i vb1root = _mm512_mask_loadu_epi32(vz, loadmask,
						&bins1[ii].root[j]);

					__m512i vsum = _mm512_add_epi32(vb1root, vb2root);
					__mmask16 mpos = loadmask & _mm512_cmpge_epi32_mask(vsum, vp);
					__mmask16 mneg = loadmask & (~mpos);

					__m512i vdiffp = _mm512_mask_sub_epi32(vp, mpos, vsum, vp);
					__m512i vdiffn = _mm512_mask_sub_epi32(vp, mneg, vp, vsum);

					mpos = _mm512_cmplt_epi32_mask(vdiffp, vi);
					mneg = _mm512_cmplt_epi32_mask(vdiffn, vi);

					vdiffp = _mm512_or_epi32(vdiffp, vpid);
					vdiffn = _mm512_or_epi32(vdiffn, vpid);

					__m512i vb1poly = _mm512_mask_loadu_epi32(vz, loadmask,
						&bins1[ii].polynum[j]);

					if (mpos > 0)
					{
						uint32_t sum[16], poly[16];

						_mm512_mask_storeu_epi32(sum, mpos, vdiffp);
						_mm512_mask_storeu_epi32(poly, mpos, _mm512_add_epi32(vb1poly, vb2poly));

						while (mpos > 0)
						{
							int pos = _trail_zcnt(mpos);
							int idx = poly[pos];

							uint32_t bsz = psize_ptr[idx];
							pslice_ptr[idx * bucketalloc + bsz] = sum[pos];
							psize_ptr[idx]++;

							mpos = _reset_lsb(mpos);
							nummatch++;
						}
					}
					
					if (mneg > 0)
					{
						uint32_t sum[16], poly[16];

						_mm512_mask_storeu_epi32(sum, mneg, vdiffn);
						_mm512_mask_storeu_epi32(poly, mneg, _mm512_add_epi32(vb1poly, vb2poly));

						while (mneg > 0)
						{
							int pos = _trail_zcnt(mneg);
							int idx = poly[pos];

							uint32_t bsz = nsize_ptr[idx];
							nslice_ptr[idx * bucketalloc + bsz] = sum[pos];
							nsize_ptr[idx]++;

							mneg = _reset_lsb(mneg);
							nummatch++;
						}
					}

					j += 16;
				}

				for (; j < bins1[ii].size; j++)
				{
					int sum1 = bin2root + bins1[ii].root[j];
					int sign1 = (sum1 >= prime);
					sum1 = sign1 ? sum1 - prime : prime - sum1;
				
					if (sum1 < interval)
					{
						int polysum = bins1[ii].polynum[j] + bins2[b].polynum[k];
						int idx = polysum;
				
						if (sign1)
						{
							uint32_t bsz = psize_ptr[idx];
							pslice_ptr[idx * bucketalloc + bsz] = (pid | sum1);
							psize_ptr[idx]++;
						}
						else
						{
							uint32_t bsz = nsize_ptr[idx];
							nslice_ptr[idx * bucketalloc + bsz] = (pid | sum1);
							nsize_ptr[idx]++;
						}
				
						nummatch++;
					}
				
				}

#else

				for (j = 0; j < MIN(bins1[ii].size, bins1b[ii].size); j++)
				{
					uint32_t sum1 = bin2root + bins1[ii].root[j];
					uint32_t sum2 = bin2root + bins1b[ii].root[j];
					int sign1 = (sum1 >= prime);
					int sign2 = (sum2 >= prime);
					sum1 = sign1 ? sum1 - prime : prime - sum1;
					sum2 = sign2 ? sum2 - prime : prime - sum2;

					if (sum1 < interval)
					{
						int polysum = bins1[ii].polynum[j] + bins2[b].polynum[k];
						int idx = polysum;

						if (sign1)
						{
							uint32_t bsz = dconf->ss_slices_p[slice].buckets[idx].size;

							//printf("adding element (%d | %d) to poly-bucket %d of size %d (alloc %d) in slice %d\n",
							//	pid, sum1, idx, bsz, dconf->ss_slices_p[slice].buckets[idx].alloc, slice);
							//fflush(stdout);

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

				if (bins1[ii].size > bins1b[ii].size)
				{
					for (; j < bins1[ii].size; j++)
					{
						uint32_t sum1 = bin2root + bins1[ii].root[j];
						int sign1 = (sum1 >= prime);
						sum1 = sign1 ? sum1 - prime : prime - sum1;

						if (sum1 < interval)
						{
							int polysum = bins1[ii].polynum[j] + bins2[b].polynum[k];
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
				else if (bins1[ii].size < bins1b[ii].size)
				{
					for (; j < bins1b[ii].size; j++)
					{
						uint32_t sum1 = bin2root + bins1b[ii].root[j];
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

#if 1

			if ((b + 1) < numbins)
			{
				for (k = 0; k < bins2[b + 1].size; k++)
				{
					__m512i vb2root = _mm512_set1_epi32(bins2[b + 1].root[k]);
					__m512i vb2poly = _mm512_set1_epi32(bins2[b + 1].polynum[k]);
					int bin2root = bins2[b + 1].root[k];

					//if (j = 0) //
					//for (j = 0; j < bins1[ii].size; j += 16)
					j = 0;
					if (bins1[ii].size > 4)
					{
						__mmask16 loadmask;

						if ((bins1[ii].size - j) >= 16)
							loadmask = 0xffff;
						else
							loadmask = (1 << (bins1[ii].size - j)) - 1;

						__m512i vb1root = _mm512_mask_loadu_epi32(vz, loadmask,
							&bins1[ii].root[j]);

						__m512i vsum = _mm512_add_epi32(vb1root, vb2root);
						__mmask16 mpos = loadmask & _mm512_cmpge_epi32_mask(vsum, vp);
						__mmask16 mneg = loadmask & (~mpos);

						__m512i vdiffp = _mm512_mask_sub_epi32(vp, mpos, vsum, vp);
						__m512i vdiffn = _mm512_mask_sub_epi32(vp, mneg, vp, vsum);

						__mmask16 mhitp = _mm512_cmplt_epi32_mask(vdiffp, vi);
						__mmask16 mhitn = _mm512_cmplt_epi32_mask(vdiffn, vi);

						vdiffp = _mm512_or_epi32(vdiffp, vpid);
						vdiffn = _mm512_or_epi32(vdiffn, vpid);

						__m512i vb1poly = _mm512_mask_loadu_epi32(vz, loadmask,
							&bins1[ii].polynum[j]);

						if (mhitp > 0)
						{
							uint32_t sum[16], poly[16];

							_mm512_mask_storeu_epi32(sum, mhitp, vdiffp);
							_mm512_mask_storeu_epi32(poly, mhitp, _mm512_add_epi32(vb1poly, vb2poly));

							while (mhitp > 0)
							{
								int pos = _trail_zcnt(mhitp);
								int idx = poly[pos];

								uint32_t bsz = psize_ptr[idx];
								pslice_ptr[idx * bucketalloc + bsz] = sum[pos];
								psize_ptr[idx]++;

								mhitp = _reset_lsb(mhitp);
								matchp1++;
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

								uint32_t bsz = nsize_ptr[idx];
								nslice_ptr[idx * bucketalloc + bsz] = sum[pos];
								nsize_ptr[idx]++;

								mhitn = _reset_lsb(mhitn);
								matchp1++;
							}
						}

						j += 16;
					}

					for (; j < bins1[ii].size; j++)
					{
						int sum1 = bin2root + bins1[ii].root[j];
						int sign1 = (sum1 >= prime);
						sum1 = sign1 ? sum1 - prime : prime - sum1;

						if (sum1 < interval)
						{
							int polysum = bins1[ii].polynum[j] + bins2[b + 1].polynum[k];
							int idx = polysum;

							if (sign1)
							{
								uint32_t bsz = psize_ptr[idx];
								pslice_ptr[idx * bucketalloc + bsz] = (pid | sum1);
								psize_ptr[idx]++;
							}
							else
							{
								uint32_t bsz = nsize_ptr[idx];
								nslice_ptr[idx * bucketalloc + bsz] = (pid | sum1);
								nsize_ptr[idx]++;
							}

							matchp1++;
						}

					}
				}
			}

			if ((b - 1) >= 0)
			{
				for (k = 0; k < bins2[b - 1].size; k++)
				{
					__m512i vb2root = _mm512_set1_epi32(bins2[b - 1].root[k]);
					__m512i vb2poly = _mm512_set1_epi32(bins2[b - 1].polynum[k]);
					int bin2root = bins2[b - 1].root[k];

					//if (j = 0) //
					//for (j = 0; j < bins1[ii].size; j += 16)
					j = 0;
					if (bins1[ii].size > 4)
					{
						__mmask16 loadmask;

						if ((bins1[ii].size - j) >= 16)
							loadmask = 0xffff;
						else
							loadmask = (1 << (bins1[ii].size - j)) - 1;

						__m512i vb1root = _mm512_mask_loadu_epi32(vz, loadmask,
							&bins1[ii].root[j]);

						__m512i vsum = _mm512_add_epi32(vb1root, vb2root);
						__mmask16 mpos = loadmask & _mm512_cmpge_epi32_mask(vsum, vp);
						__mmask16 mneg = loadmask & (~mpos);

						__m512i vdiffp = _mm512_mask_sub_epi32(vp, mpos, vsum, vp);
						__m512i vdiffn = _mm512_mask_sub_epi32(vp, mneg, vp, vsum);

						__mmask16 mhitp = _mm512_cmplt_epi32_mask(vdiffp, vi);
						__mmask16 mhitn = _mm512_cmplt_epi32_mask(vdiffn, vi);

						vdiffp = _mm512_or_epi32(vdiffp, vpid);
						vdiffn = _mm512_or_epi32(vdiffn, vpid);

						__m512i vb1poly = _mm512_mask_loadu_epi32(vz, loadmask,
							&bins1[ii].polynum[j]);

						if (mhitp > 0)
						{
							uint32_t sum[16], poly[16];

							_mm512_mask_storeu_epi32(sum, mhitp, vdiffp);
							_mm512_mask_storeu_epi32(poly, mhitp, _mm512_add_epi32(vb1poly, vb2poly));

							while (mhitp > 0)
							{
								int pos = _trail_zcnt(mhitp);
								int idx = poly[pos];

								uint32_t bsz = psize_ptr[idx];
								pslice_ptr[idx * bucketalloc + bsz] = sum[pos];
								psize_ptr[idx]++;

								mhitp = _reset_lsb(mhitp);
								matchm1++;
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

								uint32_t bsz = nsize_ptr[idx];
								nslice_ptr[idx * bucketalloc + bsz] = sum[pos];
								nsize_ptr[idx]++;

								mhitn = _reset_lsb(mhitn);
								matchm1++;
							}
						}

						j += 16;
					}

					for (; j < bins1[ii].size; j++)
					{
						int sum1 = bin2root + bins1[ii].root[j];
						int sign1 = (sum1 >= prime);
						sum1 = sign1 ? sum1 - prime : prime - sum1;

						if (sum1 < interval)
						{
							int polysum = bins1[ii].polynum[j] + bins2[b - 1].polynum[k];
							int idx = polysum;

							if (sign1)
							{
								uint32_t bsz = psize_ptr[idx];
								pslice_ptr[idx * bucketalloc + bsz] = (pid | sum1);
								psize_ptr[idx]++;
							}
							else
							{
								uint32_t bsz = nsize_ptr[idx];
								nslice_ptr[idx * bucketalloc + bsz] = (pid | sum1);
								nsize_ptr[idx]++;
							}

							matchm1++;
						}

					}
				}
			}
#endif

		}

		exit(1);


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
	printf("avg_size bins2 = %1.4f\n", avg_size3 / totalbins1);
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
		free(bins1[ii].root);
		free(bins2[ii].root);
		free(bins1[ii].polynum);
		free(bins2[ii].polynum);
	}

	free(bins1);
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
	if ((sconf->factor_base->ss_start_B > 0) &&
		(sconf->factor_base->ss_start_B < sconf->factor_base->B))
	{
		//ss_search_linked_lists(sconf, dconf);
		//ss_search_sorted_lists(sconf, dconf);
		ss_search_poly_buckets(sconf, dconf);
		//ss_search_poly_buckets_x16(sconf, dconf);
	}
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

