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

#ifdef USE_SS_SEARCH
#define SS_TIMING
#endif

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


#ifdef USE_SS_SEARCH
#define SORT_SET1
//#define USE_QSORT
#define ENABLE_BUCKET_PROTECTION

void merge(ss_set_t* ss_set, ss_set_t* ap, ss_set_t* am, int sz)
{
	int i = 0, j = 0, k = 0;

	while ((i < sz) && (j < sz)) {
		if (ap->root[i] < am->root[j]) {
			ss_set->root[k] = ap->root[i];
			ss_set->polynum[k] = ap->polynum[i];
			k++;
			i++;
		}
		else if (ap->root[i] > am->root[j]) {
			ss_set->root[k] = am->root[j];
			ss_set->polynum[k] = am->polynum[j];
			k++;
			j++;
		}
		else {
			ss_set->root[k] = ap->root[i];
			ss_set->polynum[k] = ap->polynum[i];
			k++;
			i++;
			ss_set->root[k] = am->root[j];
			ss_set->polynum[k] = am->polynum[j];
			k++;
			j++;
		}
	}

	while (i < sz)
	{
		ss_set->root[k] = ap->root[i];
		ss_set->polynum[k] = ap->polynum[i];
		k++;
		i++;
	}

	while (j < sz)
	{
		ss_set->root[k] = am->root[j];
		ss_set->polynum[k] = am->polynum[j];
		k++;
		j++;
	}

	return;
}

void shift(ss_set_t* ss_set, ss_set_t* ap, ss_set_t* am, ss_set_t* at, int m,
	int offset, int rupdate, int prime, int do_sort)
{
	// Kleinjung's shift algorithm with the following differences:
	// 1) start with the first root element with the choice of
	// all positive b's (the 0-vector in binary-encoded form).
	// 2) with that first choice, each subsequent shift just needs
	// to append a 1-bit to the end of each element of the previous 
	// sorted list, which corresponds to a modular addition.
	//
	// as soon as the modulus is needed, all subsequent additions
	// will also need it and the results of these will all be smaller
	// than their corresponding source elements.  Thus, sorting works 
	// the same way, with an array swap at the point of the first 
	// modulus operation followed by a merge.
	int sz = (1 << m);
	int km = -1;
	int k;
	int s = (1 << (m + offset));

#ifdef USE_AVX512F
	if (sz >= 16)
	{
		__m512i vu = _mm512_set1_epi32(rupdate);
		__m512i vs = _mm512_set1_epi32(s);
		__m512i vprime = _mm512_set1_epi32(prime);

		for (k = 0; k < sz; k += 16)
		{
			__m512i vr = _mm512_load_epi32(ss_set->root + k);
			__m512i vp = _mm512_load_epi32(ss_set->polynum + k);

			_mm512_store_epi32(ap->root + k, vr);
			_mm512_store_epi32(ap->polynum + k, vp);

			vr = _mm512_add_epi32(vr, vu);
			__mmask16 mask = _mm512_cmpge_epi32_mask(vr, vprime);
			vr = _mm512_mask_sub_epi32(vr, mask, vr, vprime);

			if ((mask > 0) && (km == -1))
			{
				int pos = _trail_zcnt(mask);
				km = k + pos;
			}

			_mm512_store_epi32(am->root + k, vr);
			_mm512_store_epi32(am->polynum + k, _mm512_add_epi32(vp, vs));
		}
	}
	else
#endif
	{
		for (k = 0; k < sz; k++)
		{
			int poly = ss_set->polynum[k];
			int root = ss_set->root[k];
			int rm = root + rupdate;
			ap->root[k] = ss_set->root[k];
			ap->polynum[k] = ss_set->polynum[k];

			if (rm >= prime)
			{
				if (km == -1)
					km = k;
				rm = rm - prime;
			}

			am->root[k] = rm;
			am->polynum[k] = poly + s;
		}
	}

	if (do_sort)
	{

		if (km > 0)
		{
			for (k = 0; k < km; k++)
			{
				at->root[k] = am->root[k];
				at->polynum[k] = am->polynum[k];
			}

			for (k = km; k < sz; k++)
			{
				am->root[k - km] = am->root[k];
				am->polynum[k - km] = am->polynum[k];
			}

			for (k = 0; k < km; k++)
			{
				am->root[sz - km + k] = at->root[k];
				am->polynum[sz - km + k] = at->polynum[k];
			}
		}

		merge(ss_set, ap, am, sz);
	}
	else
	{
		for (k = 0; k < sz; k++)
		{
			ss_set->root[2 * k + 0] = ap->root[k];
			ss_set->root[2 * k + 1] = am->root[k];
			ss_set->polynum[2 * k + 0] = ap->polynum[k];
			ss_set->polynum[2 * k + 1] = am->polynum[k];
		}
	}

	return;
}

void shift_unsorted(ss_set_t* ss_set, int m, int offset, int rupdate, int prime, int index)
{
	int sz = (1 << m);
	int km = -1;
	int k;
	int s = (1 << (m + offset));

#ifdef USE_AVX512F
	if (sz >= 16)
	{
		__m512i vu = _mm512_set1_epi32(rupdate);
		__m512i vs = _mm512_set1_epi32(s);
		__m512i vprime = _mm512_set1_epi32(prime);

		for (k = 0; k < sz; k += 16)
		{
			__m512i vr = _mm512_load_epi32(ss_set->root + k + index);
			__m512i vp = _mm512_load_epi32(ss_set->polynum + k + index);

			vr = _mm512_add_epi32(vr, vu);
			__mmask16 mask = _mm512_cmpge_epi32_mask(vr, vprime);
			vr = _mm512_mask_sub_epi32(vr, mask, vr, vprime);

			_mm512_store_epi32(ss_set->root + sz + k + index, vr);
			_mm512_store_epi32(ss_set->polynum + sz + k + index, _mm512_add_epi32(vp, vs));
		}
	}
	else
#endif
	{
		for (k = 0; k < sz; k++)
		{
			int poly = ss_set->polynum[k + index];
			int root = ss_set->root[k + index];
			int rm = root + rupdate;

			ss_set->root[sz + k + index] = rm;
			ss_set->polynum[sz + k + index] = poly + s;
		}
	}

	return;
}

#ifdef USE_QSORT
int binary_find_ge(ss_sortable_set_t* ss_set, int set_sz, int element)
{
	int lo_idx = 0;
	int hi_idx = set_sz - 1;

	while ((hi_idx - lo_idx) > 1)
	{
		int mid_idx = lo_idx + (hi_idx - lo_idx) / 2;

		if (ss_set[mid_idx].root > element)
		{
			hi_idx = mid_idx;
		}
		else if (ss_set[mid_idx].root < element)
		{
			lo_idx = mid_idx;
		}
		else
		{
			return mid_idx;
		}
	}

	if (element > ss_set[hi_idx].root)
		return set_sz;
	else
		return hi_idx;
}
#else
int binary_find_ge(int* root, int set_sz, int element)
{
	int lo_idx = 0;
	int hi_idx = set_sz - 1;

	while ((hi_idx - lo_idx) > 1)
	{
		int mid_idx = lo_idx + (hi_idx - lo_idx) / 2;

		if ((uint32_t)root[mid_idx] > (uint32_t)element)
		{
			hi_idx = mid_idx;
		}
		else if ((uint32_t)root[mid_idx] < (uint32_t)element)
		{
			lo_idx = mid_idx;
		}
		else
		{
			return mid_idx;
		}
	}

	if ((uint32_t)element > (uint32_t)root[hi_idx])
		return set_sz;
	else
		return hi_idx;
}
#endif
#endif



#if defined( USE_POLY_BUCKET_SS )

#define TRY_GATHER_SCATTER

void ss_search_clear(static_conf_t* sconf, dynamic_conf_t* dconf)
{

	return;
}

void ss_search_setup(static_conf_t* sconf, dynamic_conf_t* dconf)
{
	siqs_poly* poly = dconf->curr_poly;
	fb_list* fb = sconf->factor_base;

	int i, j, ii;
	int ss_num_poly_terms = poly->s / 2;

	mpz_set(dconf->polyb1, dconf->Bl[0]);
	for (ii = 1; ii < ss_num_poly_terms; ii++) {
		mpz_add(dconf->polyb1, dconf->polyb1, dconf->Bl[ii]);
	}
	mpz_tdiv_q_2exp(dconf->polyb1, dconf->polyb1, 1);

	ss_num_poly_terms = (poly->s - 1) - ss_num_poly_terms;

	// first poly b: sum of first n positive Bl
	mpz_set(dconf->polyb2, dconf->Bl[ii]);
	for (ii++; ii < poly->s; ii++) {
		mpz_add(dconf->polyb2, dconf->polyb2, dconf->Bl[ii]);
	}
	mpz_tdiv_q_2exp(dconf->polyb2, dconf->polyb2, 1);

	// create a map from gray code enumeration order
	// to polynomial binary encoding.
	int polynum = 0;

	dconf->polymap[1] = 0;
	dconf->polynums[0] = 0;

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
		dconf->polyv[ii] = v;
		if (tmp & 1)
			dconf->polysign[ii] = -1;
		else
			dconf->polysign[ii] = 1;

		// next polynum
		polynum ^= (1 << (v - 1));
		dconf->polynums[ii] = polynum;
		dconf->polymap[ii + 1] = polynum;
	}

	int num_bpoly = 1 << (dconf->curr_poly->s - 1);

	if (num_bpoly > dconf->ss_slices_p[0].numbuckets)
	{
		printf("not enough buckets allocated (%d) for current number of b-polys (%d)\n",
			dconf->ss_slices_p[0].numbuckets, num_bpoly);
		exit(1);
	}

	if ((dconf->ss_slices_p[0].numbuckets > (4 * num_bpoly)) ||
		(dconf->ss_slices_p[0].numbuckets < num_bpoly))
	{
		printf("reallocating number of poly-buckets to %d\n", (2 * num_bpoly));
	}

	for (i = 0; i < dconf->num_ss_slices; i++)
	{
		if (dconf->ss_slices_p[i].numbuckets > (2 * num_bpoly))
		{
			int a = dconf->ss_slices_p[i].alloc;
			int nb = 2 * num_bpoly;
			dconf->ss_slices_p[i].numbuckets = nb;
			dconf->ss_slices_p[i].elements = 
				(uint32_t*)xrealloc(dconf->ss_slices_p[i].elements, a * nb * sizeof(uint32_t));
		}

		for (ii = 0; ii <= (2 * num_bpoly); ii++)
		{
			dconf->ss_slices_p[i].size[ii] = 0;
#ifndef USE_POLY_BUCKET_PN_COMBINED_VARIATION
			dconf->ss_slices_n[i].size[ii] = 0;
#endif
		}
	}

	return;
}

// for this version of subset-sum, we need a different type
// for the enumerated roots so it is easier to sort them
typedef struct
{
	uint32_t root;
	uint32_t polynum;
} ss_sortable_set_t;

int qsort_ss_set(const void* x, const void* y)
{
	ss_sortable_set_t* xx = (ss_sortable_set_t*)x;
	ss_sortable_set_t* yy = (ss_sortable_set_t*)y;

	if (xx->root > yy->root)
		return 1;
	else if (xx->root == yy->root)
		return 0;
	else
		return -1;
}

#define QS_COMPARE(a,b) ((a)-(b))  

void QuickSort(int* list, int* list2, int beg, int end)
{
	int piv, tmp;
	int tmp2, piv2;

	int  l, r, p;

	while (beg < end)    // This while loop will avoid the second recursive call  
	{
		l = beg; p = (beg + end) / 2; r = end;

		piv = list[p];
		piv2 = list2[p];

		while (1)
		{
			while ((l <= r) && (QS_COMPARE(list[l], piv) <= 0)) l++;
			while ((l <= r) && (QS_COMPARE(list[r], piv) > 0)) r--;

			if (l > r) break;

			tmp = list[l];  list[l] = list[r];   list[r] = tmp;
			tmp2 = list2[l]; list2[l] = list2[r]; list2[r] = tmp2;

			if (p == r) p = l;

			l++; r--;
		}

		list[p] = list[r]; list[r] = piv;
		list2[p] = list2[r]; list2[r] = piv2;
		r--;

		// Recursion on the shorter side & loop (with new indexes) on the longer  
		if ((r - beg) < (end - l))
		{
			QuickSort(list, list2, beg, r);
			beg = l;
		}
		else
		{
			QuickSort(list, list2, l, end);
			end = r;
		}
	}
}

// this version performs the search slower, but doesn't
// need AVX512 and efficiently finds all of the sieve events.
// the larger number of sieve hits results in an overall
// larger relation discovery rate (rels/sec).
void ss_search_poly_buckets_sorted(static_conf_t* sconf, dynamic_conf_t* dconf)
{
	// the subset-sum search algorithm using a bucket sort to
	// organize the sieve-hits per poly.  
	siqs_poly* poly = dconf->curr_poly;
	fb_list* fb = sconf->factor_base;
	uint32_t numblocks = sconf->num_blocks;
	int interval;
	int slicesz = sconf->factor_base->slice_size;
	int pid_offset = SS_SIGN_BIT + 1;

	numblocks = sconf->num_blocks;
	interval = (numblocks << 15);
	int interval2 = interval * 2;

	int i, j, k, ii;

#ifdef USE_QSORT
	ss_sortable_set_t *ss_set1;
	ss_sortable_set_t *ss_set2;
#else
	ss_set_t ss_set1;
	ss_set_t ss_set2;
	ss_set_t ss_set_t1;
	ss_set_t ss_set_t2;
	ss_set_t ss_set_t3;
	ss_set_t ss_tmp1;
	ss_set_t ss_tmp2;
#endif

	double t_enum_roots;
	double t_sort_roots;
	double t_match_roots;
	struct timeval start, stop;

	int* rootupdates = dconf->rootupdates;
	update_t update_data = dconf->update_data;
	uint32_t* modsqrt = sconf->modsqrt_array;

	uint64_t nummatch = 0;
	int nump = 0;
	
#ifdef USE_QSORT
	ss_set1 = (ss_sortable_set_t*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(ss_sortable_set_t));
	ss_set2 = (ss_sortable_set_t*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(ss_sortable_set_t));
#else
	ss_set1.root = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_set2.root = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_tmp1.root = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_tmp2.root = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_set1.polynum = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_set2.polynum = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_tmp1.polynum = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_tmp2.polynum = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_set_t1.root = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_set_t2.root = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_set_t3.root = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_set_t1.polynum = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_set_t2.polynum = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_set_t3.polynum = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
#endif

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
		int prime;
		int root1;
		int root2;
		uint8_t logp = fb->list->logprime[i];

		if (slice >= sconf->factor_base->num_ss_slices)
		{
			printf("error slice number %d too large, %d allocated, at fb index %d (prime %u)\n",
				slice, sconf->factor_base->num_ss_slices, i, prime);
			exit(0);
		}

		nump++;

		// create set1, an enumeration of (+/-t - b)(a)^-1 mod p
		// for b in the set [+/-B1 +/- B2 +/- ... Bs/2]
		int ii, v, sign, n = poly->s / 2 - 1;
		int tmp;

#ifdef SS_TIMING
		gettimeofday(&start, NULL);
#endif

		prime = fb->list->prime[i];
		r1 = dconf->firstroot1a[i];
		r2 = dconf->firstroot1b[i];

		int* ptr;

		// enumerate set1 roots
#if defined(SORT_SET1)

		ss_tmp1.root[0] = r1;
		ss_tmp1.polynum[0] = 0;

		ss_tmp1.root[1] = ss_tmp1.root[0] + rootupdates[(0)*bound + i];
		ss_tmp1.polynum[1] = ss_tmp1.polynum[0] + (1 << 0);

		if (ss_tmp1.root[1] >= prime)
		{
			ss_tmp1.root[1] = ss_tmp1.root[1] - prime;
			int rtmp = ss_tmp1.root[0];
			int ptmp = ss_tmp1.polynum[0];
			ss_tmp1.root[0] = ss_tmp1.root[1];
			ss_tmp1.polynum[0] = ss_tmp1.polynum[1];
			ss_tmp1.root[1] = rtmp;
			ss_tmp1.polynum[1] = ptmp;
		}

		for (ii = 1; ii < n; ii++)
		{
			shift(&ss_tmp1, &ss_set_t1, &ss_set_t2, &ss_set_t3, ii, 0,
				rootupdates[(ii)*bound + i], prime, 1);
		}

		ss_tmp2.root[0] = r2;
		ss_tmp2.polynum[0] = 0;

		ss_tmp2.root[1] = ss_tmp2.root[0] + rootupdates[(0) * bound + i];
		ss_tmp2.polynum[1] = ss_tmp2.polynum[0] + (1 << 0);

		if (ss_tmp2.root[1] >= prime)
		{
			ss_tmp2.root[1] = ss_tmp2.root[1] - prime;
			int rtmp = ss_tmp2.root[0];
			int ptmp = ss_tmp2.polynum[0];
			ss_tmp2.root[0] = ss_tmp2.root[1];
			ss_tmp2.polynum[0] = ss_tmp2.polynum[1];
			ss_tmp2.root[1] = rtmp;
			ss_tmp2.polynum[1] = ptmp;
		}

		for (ii = 1; ii < n; ii++)
		{
			shift(&ss_tmp2, &ss_set_t1, &ss_set_t2, &ss_set_t3, ii, 0,
				rootupdates[(ii)*bound + i], prime, 1);
		}

		merge(&ss_set1, &ss_tmp1, &ss_tmp2, 1 << n);

#else
		for (ii = 1, k = 0, polynum = 0; ii < (1 << n); ii++, k += 2) {
			// these roots go into the set
#ifdef USE_QSORT
			ss_set1[k + 0].root = r1;
			ss_set1[k + 1].root = r2;
			ss_set1[k + 0].polynum = polynum;
			ss_set1[k + 1].polynum = polynum;
#else
			ss_set1.root[k + 0] = r1;
			ss_set1.root[k + 1] = r2;
			ss_set1.polynum[k + 0] = polynum;
			ss_set1.polynum[k + 1] = polynum;
#endif

			// next polynum
			polynum = dconf->polynums[ii];
			sign = dconf->polysign[ii];
			v = dconf->polyv[ii];

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
#ifdef USE_QSORT
		ss_set1[k + 0].root = r1;
		ss_set1[k + 1].root = r2;
		ss_set1[k + 0].polynum = polynum;
		ss_set1[k + 1].polynum = polynum;
#else
		ss_set1.root[k + 0] = r1;
		ss_set1.root[k + 1] = r2;
		ss_set1.polynum[k + 0] = polynum;
		ss_set1.polynum[k + 1] = polynum;

#ifdef SORT_SET1
		// sort the 1st set.  This isn't necessary to do matching, its
		// an experiment to see if we can reduce the number of searches needed.
		// Search for the match to the largest root in a small range (say, of 16)
		// and then iterate until the smallest root in the range no longer matches
		// anything.  With SIMD we can check many matches at a time.  So, there
		// is an extra cost to the sorting with the tradeoff of 16x fewer
		// searches...
		QuickSort(ss_set1.root, ss_set1.polynum, 0, (2 << n) - 1);
#endif

#endif
#endif

		int size1 = n;

		// create set2, an enumeration of (+/-t - b)(a)^-1 mod p
		// for b in the set [+/-Bs/2+1 +/- Bs/2+2 ... + Bs]
		n = (poly->s - 1) - n;

		int size2 = n;
		r1 = dconf->firstroot2[i];

		// always build set2 sorted, using shift & merge
		int n0 = size1;

		ss_set2.root[0] = r1;
		ss_set2.polynum[0] = 0;

		ss_set2.root[1] = ss_set2.root[0] + rootupdates[(n0)*bound + i];
		ss_set2.polynum[1] = ss_set2.polynum[0] + (1 << n0);

		if (ss_set2.root[1] >= prime)
		{
			ss_set2.root[1] = ss_set2.root[1] - prime;
			int rtmp = ss_set2.root[0];
			int ptmp = ss_set2.polynum[0];
			ss_set2.root[0] = ss_set2.root[1];
			ss_set2.polynum[0] = ss_set2.polynum[1];
			ss_set2.root[1] = rtmp;
			ss_set2.polynum[1] = ptmp;
		}

		for (ii = n0 + 1; ii < poly->s - 1; ii++)
		{
			shift(&ss_set2, &ss_set_t1, &ss_set_t2, &ss_set_t3, ii - n0, n0,
				rootupdates[(ii)*bound + i], prime, 1);
		}

#ifdef SS_TIMING
		gettimeofday(&stop, NULL);
		t_enum_roots += ytools_difftime(&start, &stop);

		start.tv_sec = stop.tv_sec;
		start.tv_usec = stop.tv_usec;
#endif


		// admissible sieving events occur when r1 + r2 < interval
		// or p <= r1 + r2 < p + interval
		// which is equivalent to 
		// -r0 <= r1 < interval - r0 and
		// p - r0 <= r1 < p + interval - r0
		// from here there are two cases:
		// 1) where r0 >= interval, then the first inequality is impossible
		// and we only need to check the second
		// 2) r0 < interval and both inequalities are possible
		// for case 1) do two binary searches on the sorted set2 to find the
		// first root r1' such that r1' >= p - r0, and find the last r1'' such
		// that r1'' < p + interval - r0.  All r1 between these points are
		// sieving events.
		// the 2nd search could probably be replaced by a comparison during
		// a loop starting at r1', especially if that loop is vectorized.
		// (indeed, this is faster.)
		// for case 2) the first substep just iterates through set2 starting
		// from the beginning and stops once r1'' >= interval - r0.
		// The second substep finds by binary search the element r1' >= p - r0
		// and then from there iterates to the end of the set2 array. case 2b
		// is logically equivalent to case 1, and is combined below.
		
		uint32_t pid = (uint32_t)(i - fboffset) << pid_offset;

		if ((i - fboffset) >= slicesz)
		{
			printf("pid %u too large for slice %d with fboffset %d\n",
				pid, slice, fboffset);
			exit(1);
		}

		uint32_t bucketalloc = dconf->ss_slices_p[0].alloc;
		uint32_t* pslice_ptr = dconf->ss_slices_p[slice].elements;
		uint32_t* psize_ptr = dconf->ss_slices_p[slice].size;
		int num_set2 = 1 << size2;

#ifndef SORT_SET1
		// commence matching
		int max_id = num_set2;
		for (ii = 0; ii < (2 << size1); ii++)
		{
#ifdef USE_QSORT
			root1 = ss_set1[ii].root;
			int poly1 = ss_set1[ii].polynum;
#else
			root1 = ss_set1.root[ii];
			int poly1 = ss_set1.polynum[ii];
#endif

			// case 2a only occurs if the initial root is really small
			if (root1 < interval)
			{
				int startid = 0;
#ifdef USE_QSORT
				root2 = ss_set2[startid].root;
#else
				root2 = ss_set2.root[startid];
#endif

				// case 2a: beginning of set2 array matching
				while (((root1 + root2) < interval) && (startid < num_set2))
				{
#ifdef USE_QSORT
					int polysum = ss_set2[startid].polynum + poly1;
#else
					int idx = ss_set2.polynum[startid] + poly1;
#endif

#ifdef FULL_INTERVAL_NOTATION
					int sum = root1 + root2 + interval;
					int sign = sum >= interval;

					if (sign == 0) sum = interval - sum;
#else
					int sum = root1 + root2;
#endif

					uint32_t bsz = psize_ptr[idx];
					pslice_ptr[idx * bucketalloc + bsz] = (pid | sum);
					psize_ptr[idx]++;

#ifdef ENABLE_BUCKET_PROTECTION
					if (psize_ptr[idx] >= dconf->ss_slices_p[slice].alloc)
					{
						printf("error bucket overflow: size of poly bucket %d = %d\n",
							idx, psize_ptr[idx]);
						exit(1);
					}
#endif

					startid++;
#ifdef USE_QSORT
					root2 = ss_set2[startid].root;
#else
					root2 = ss_set2.root[startid];
#endif

#ifdef SS_TIMING
					nummatch++;
#endif
				}
			}

			// this case always occurs, either as case 1
			// or as case 2b
			{
				int startid;

#ifdef USE_QSORT
				int startid = binary_find_ge(ss_set2, num_set2, prime - root1 - interval);
				root2 = ss_set2[startid].root;
#else

#ifdef SORT_SET1

				// if set1 is sorted then we can gradually decrease the max index
				// of the search.
				startid = binary_find_ge(ss_set2.root, max_id, prime - root1 - interval);
				max_id = startid;
#else
				startid = binary_find_ge(ss_set2.root, num_set2, prime - root1 - interval);
#endif

				root2 = ss_set2.root[startid];

#endif

#ifndef USE_POLY_BUCKET_PN_COMBINED_VARIATION
				root1 = root1 - prime + interval;
#else
				root1 = root1 - prime;
#endif

				// case 1: middle of set2 array matching
#ifndef USE_POLY_BUCKET_PN_COMBINED_VARIATION
				while (((startid) < num_set2) &&
					((root1 + root2) < interval2))
#else
				while (((startid) < num_set2) &&
					((root1 + root2) < interval))
#endif
				{
#ifdef USE_QSORT
					int idx1 = ss_set2[startid].polynum + poly1;
#else
					int idx1 = ss_set2.polynum[startid] + poly1;
#endif
					uint32_t bsz = psize_ptr[idx1];

#ifndef USE_POLY_BUCKET_PN_COMBINED_VARIATION
#ifdef USE_QSORT
					int sum1 = root1 + ss_set2[startid].root;
#else
					int sum1 = root1 + ss_set2.root[startid];
#endif

					int sign1 = sum1 >= interval;

					// the smoothness checking code all assumes -side hits
					// are arranged with respect to the 0 point, symmetric 
					// to +side hits.  so if this is a -side hit we need
					// to flip its orientation relative to the 0 point so
					// that smoothness checking works.  Note that the
					// exact 0 point is ignored for now because it raises
					// errors in trial division... need to look into how to
					// properly handle hits at sum = interval.
					if (sign1 == 0) sum1 = interval - sum1;

#else
#ifdef USE_QSORT
					int sum1 = root1 + ss_set2[startid].root;
#else
					int sum1 = root1 + root2;
#endif
					int sign = sum1 >= 0;
					if (sign == 0)
					{
						sum1 = 0 - sum1;
						sum1 |= SS_MAX_ROOT;
					}
#endif

					pslice_ptr[idx1 * bucketalloc + bsz] = (pid | sum1);
					psize_ptr[idx1]++;

#ifdef ENABLE_BUCKET_PROTECTION
					if (psize_ptr[idx1] >= dconf->ss_slices_p[slice].alloc)
					{
						printf("error bucket overflow: size of poly bucket %d = %d\n",
							idx1, psize_ptr[idx1]);
						exit(1);
					}
#endif

					startid++;
#ifdef USE_QSORT
					root2 = ss_set2[startid].root;
#else
					root2 = ss_set2.root[startid];
#endif

#ifdef SS_TIMING
					nummatch++;
#endif
				}

				//printf("\n");
			}

			
		}

		//exit(1);

#else

		
		__m512i vpid = _mm512_set1_epi32(pid);
		__m512i vprime = _mm512_set1_epi32(prime);
		__m512i vinterval = _mm512_set1_epi32(interval);
		__m512i vz = _mm512_setzero_epi32();
		__m512i vsignbit = _mm512_set1_epi32(SS_MAX_ROOT);
		// with SIMD we can match much faster.  It doesn't find some hits compared
		// to non-SIMD, so there must be a bug somewhere.  But it doesn't really
		// matter because even with faster matching, subset-sum still isn't
		// better than normal.
		int max_id = num_set2;
		int bucket_alloc_bits = dconf->poly_buckets_allocated;
		for (ii = 0; ii < (2 << size1); ii += 16)
		{
			// this case always occurs, either as case 1
			// or as case 2b
			int startid;

			// the idea behind this SIMD is to find the match to the
			// largest root1 in this group, thus the smallest root2.
			// then as we increment root2, smaller root1s in the group
			// will start to hit.  Stop once the smallest root1 is no
			// longer possible to match.
			startid = binary_find_ge(ss_set2.root, max_id,
				prime - ss_set1.root[ii + 15] - interval);
			max_id = startid;

			__m512i vr0 = _mm512_load_epi32(ss_set1.root + ii);
			__m512i vr1 = _mm512_sub_epi32(vr0, vprime);
			__m512i vr2 = _mm512_set1_epi32(ss_set2.root[startid]);
			__m512i vp1 = _mm512_load_epi32(ss_set1.polynum + ii);
			__m512i vp2 = _mm512_set1_epi32(ss_set2.polynum[startid]);
			__mmask16 hits;

			// case1 and case2b
			while (((startid < num_set2)) && 
				((ss_set1.root[ii + 0] - prime + ss_set2.root[startid]) < interval))

			{
				vp2 = _mm512_add_epi32(vp1, vp2);
				vr2 = _mm512_add_epi32(vr1, vr2);
				hits = _mm512_cmplt_epi32_mask(_mm512_abs_epi32(vr2), vinterval);

				int numhits = _mm_popcnt_u32(hits);

				if (numhits > 0)
				{
					__mmask16 signs = _mm512_cmplt_epi32_mask(vr2, vz);
					vr2 = _mm512_mask_sub_epi32(vr2, signs, vz, vr2);
					vr2 = _mm512_mask_or_epi32(vr2, signs, vr2, vsignbit);
					vr2 = _mm512_or_epi32(vr2, vpid);

					__m512i vbsz = _mm512_mask_i32gather_epi32(vbsz, hits, vp2, psize_ptr, 4);
					__m512i vindex = _mm512_slli_epi32(vp2, bucket_alloc_bits);
					vindex = _mm512_add_epi32(vindex, vbsz);

					_mm512_mask_i32scatter_epi32(pslice_ptr, hits, vindex, vr2, 4);
					vbsz = _mm512_add_epi32(vbsz, _mm512_set1_epi32(1));
					_mm512_mask_i32scatter_epi32(psize_ptr, hits, vp2, vbsz, 4);

					nummatch += numhits;
				}

				startid++;
				vr2 = _mm512_set1_epi32(ss_set2.root[startid]);
				vp2 = _mm512_set1_epi32(ss_set2.polynum[startid]);

			}
		}

		// case 2a only occurs if the initial root is really small.
		// with sorted root1s, we can stop looking as soon as the
		// root comparison fails.
		for (ii = 0; ii < (2 << size1); ii++)
		{
			root1 = ss_set1.root[ii];
			if (root1 < interval)
			{
				int startid = 0;
				root2 = ss_set2.root[startid];

				// case 2a: beginning of set2 array matching
				while ((root1 + root2) < interval) // ((startid < num_set2) && 
				{
					int polysum = ss_set2.polynum[startid] + ss_set1.polynum[ii];
					int idx = polysum;
					int sum = root1 + root2;

					uint32_t bsz = psize_ptr[idx];
					pslice_ptr[idx * bucketalloc + bsz] = (pid | sum);
					psize_ptr[idx]++;

#ifdef ENABLE_BUCKET_PROTECTION
					if (psize_ptr[idx] >= dconf->ss_slices_p[slice].alloc)
					{
						printf("error bucket overflow: size of poly bucket %d = %d\n",
							idx, psize_ptr[idx]);
						exit(1);
					}
#endif

					startid++;
					root2 = ss_set2.root[startid];

#ifdef SS_TIMING
					nummatch++;
#endif
				}
			}
			else
			{
				break;
			}
		}

		//exit(1);

#endif


#ifdef SS_TIMING
		gettimeofday(&stop, NULL);
		t_match_roots += ytools_difftime(&start, &stop);
#endif
	}

#ifdef SS_TIMING
	printf("ran subset-sum on %d primes\n", nump);
	printf("found %lu sieve hits\n", nummatch);
	printf("enumerating roots: %1.4f seconds\n", t_enum_roots);
	printf("matching roots: %1.4f seconds\n", t_match_roots);
#endif

#ifdef USE_QSORT
	free(ss_set1);
	free(ss_set2);
#else
	free(ss_set1.root);
	free(ss_set2.root);
	free(ss_tmp1.root);
	free(ss_tmp2.root);
	free(ss_set1.polynum);
	free(ss_set2.polynum);
	free(ss_tmp1.polynum);
	free(ss_tmp2.polynum);
	free(ss_set_t1.root);
	free(ss_set_t2.root);
	free(ss_set_t3.root);
	free(ss_set_t1.polynum);
	free(ss_set_t2.polynum);
	free(ss_set_t3.polynum);
#endif
	

	return;
}

#endif

#if defined( USE_DIRECT_SIEVE_SS )
void ss_search_clear(static_conf_t* sconf, dynamic_conf_t* dconf)
{

	return;
}

void ss_search_setup(static_conf_t* sconf, dynamic_conf_t* dconf)
{
	siqs_poly* poly = dconf->curr_poly;
	fb_list* fb = sconf->factor_base;

	int i, j, ii;
	int ss_num_poly_terms = poly->s / 2;

	mpz_set(dconf->polyb1, dconf->Bl[0]);
	for (ii = 1; ii < ss_num_poly_terms; ii++) {
		mpz_add(dconf->polyb1, dconf->polyb1, dconf->Bl[ii]);
	}
	mpz_tdiv_q_2exp(dconf->polyb1, dconf->polyb1, 1);

	ss_num_poly_terms = (poly->s - 1) - ss_num_poly_terms;

	// first poly b: sum of first n positive Bl
	mpz_set(dconf->polyb2, dconf->Bl[ii]);
	for (ii++; ii < poly->s; ii++) {
		mpz_add(dconf->polyb2, dconf->polyb2, dconf->Bl[ii]);
	}
	mpz_tdiv_q_2exp(dconf->polyb2, dconf->polyb2, 1);

	// create a map from gray code enumeration order
	// to polynomial binary encoding.
	int polynum = 0;

	//polymap[0] = 1;
	dconf->polymap[1] = 0;
	dconf->polynums[0] = 0;

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
		dconf->polyv[ii] = v;
		if (tmp & 1)
			dconf->polysign[ii] = -1;
		else
			dconf->polysign[ii] = 1;

		// next polynum
		polynum ^= (1 << (v - 1));
		dconf->polynums[ii] = polynum;
		dconf->polymap[ii + 1] = polynum;

		//printf("%d <- %d (%d)\n", ii + 1, polynum, polynums[ii]);
	}

	return;
}

void ss_search_poly_buckets_sorted(static_conf_t* sconf, dynamic_conf_t* dconf)
{
	// the subset-sum search algorithm using a bucket sort to
	// organize the sieve-hits per poly.  
	siqs_poly* poly = dconf->curr_poly;
	fb_list* fb = sconf->factor_base;
	uint32_t numblocks = sconf->num_blocks;
	int interval;
	int slicesz = sconf->factor_base->slice_size;

	interval = dconf->ss_sieve_sz;
	int interval2 = interval * 2;
	int resieve = dconf->num_ss_slices;
	int invalid_report_count = 0;
	int i, j, k, ii;

	ss_set_t ss_set1;
	ss_set_t ss_set2;
	ss_set_t ss_set_t1;
	ss_set_t ss_set_t2;
	ss_set_t ss_set_t3;
	ss_set_t ss_tmp1;
	ss_set_t ss_tmp2;

	double t_enum_roots;
	double t_sort_roots;
	double t_match_roots;
	struct timeval start, stop;

	int* rootupdates = dconf->rootupdates;
	update_t update_data = dconf->update_data;
	uint32_t* modsqrt = sconf->modsqrt_array;

	int nummatch = 0;
	int nump = 0;

	ss_set1.root = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_set2.root = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_tmp1.root = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_tmp2.root = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_set1.polynum = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_set2.polynum = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_tmp1.polynum = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_tmp2.polynum = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_set_t1.root = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_set_t2.root = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_set_t3.root = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_set_t1.polynum = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_set_t2.polynum = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));
	ss_set_t3.polynum = (int*)xmalloc((2 << (MAX_A_FACTORS / 2)) * sizeof(int));

	for (i = fb->ss_start_B; i < fb->B; i++)
	{
		int numB;
		int maxB = 1 << (dconf->curr_poly->s - 1);
		int r1 = update_data.firstroots1[i], r2 = update_data.firstroots2[i];
		uint32_t bound = sconf->factor_base->B;
		uint32_t med_B = sconf->factor_base->med_B;
		int polynum = 0;
		int slice = (int)((i - sconf->factor_base->ss_start_B) / slicesz);
		int prime;
		int root1;
		int root2;
		uint8_t logp = fb->list->logprime[i];

		if (slice >= sconf->factor_base->num_ss_slices)
		{
			printf("error slice number %d too large, %d allocated, at fb index %d (prime %u)\n",
				slice, sconf->factor_base->num_ss_slices, i, prime);
			exit(0);
		}

		nump++;

		// create set1, an enumeration of (+/-t - b)(a)^-1 mod p
		// for b in the set [+/-B1 +/- B2 +/- ... Bs/2]
		int ii, v, sign, n = poly->s / 2 - 1;
		int tmp;

#ifdef SS_TIMING
		gettimeofday(&start, NULL);
#endif

		prime = fb->list->prime[i];
		r1 = dconf->firstroot1a[i];
		r2 = dconf->firstroot1b[i];

		int* ptr;

		// enumerate set1 roots
#if defined(SORT_SET1)

		ss_tmp1.root[0] = r1;
		ss_tmp1.polynum[0] = 0;

		ss_tmp1.root[1] = ss_tmp1.root[0] + rootupdates[(0) * bound + i];
		ss_tmp1.polynum[1] = ss_tmp1.polynum[0] + (1 << 0);

		if (ss_tmp1.root[1] >= prime)
		{
			ss_tmp1.root[1] = ss_tmp1.root[1] - prime;
			int rtmp = ss_tmp1.root[0];
			int ptmp = ss_tmp1.polynum[0];
			ss_tmp1.root[0] = ss_tmp1.root[1];
			ss_tmp1.polynum[0] = ss_tmp1.polynum[1];
			ss_tmp1.root[1] = rtmp;
			ss_tmp1.polynum[1] = ptmp;
		}

		for (ii = 1; ii < n; ii++)
		{
			shift(&ss_tmp1, &ss_set_t1, &ss_set_t2, &ss_set_t3, ii, 0,
				rootupdates[(ii)*bound + i], prime, 1);
		}

		ss_tmp2.root[0] = r2;
		ss_tmp2.polynum[0] = 0;

		ss_tmp2.root[1] = ss_tmp2.root[0] + rootupdates[(0) * bound + i];
		ss_tmp2.polynum[1] = ss_tmp2.polynum[0] + (1 << 0);

		if (ss_tmp2.root[1] >= prime)
		{
			ss_tmp2.root[1] = ss_tmp2.root[1] - prime;
			int rtmp = ss_tmp2.root[0];
			int ptmp = ss_tmp2.polynum[0];
			ss_tmp2.root[0] = ss_tmp2.root[1];
			ss_tmp2.polynum[0] = ss_tmp2.polynum[1];
			ss_tmp2.root[1] = rtmp;
			ss_tmp2.polynum[1] = ptmp;
		}

		for (ii = 1; ii < n; ii++)
		{
			shift(&ss_tmp2, &ss_set_t1, &ss_set_t2, &ss_set_t3, ii, 0,
				rootupdates[(ii)*bound + i], prime, 1);
		}

		merge(&ss_set1, &ss_tmp1, &ss_tmp2, 1 << n);

#else
		for (ii = 1, k = 0, polynum = 0; ii < (1 << n); ii++, k += 2) {
			// these roots go into the set

			ss_set1.root[k + 0] = r1;
			ss_set1.root[k + 1] = r2;
			ss_set1.polynum[k + 0] = polynum;
			ss_set1.polynum[k + 1] = polynum;

			// next polynum
			polynum = dconf->polynums[ii];
			sign = dconf->polysign[ii];
			v = dconf->polyv[ii];

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

		ss_set1.root[k + 0] = r1;
		ss_set1.root[k + 1] = r2;
		ss_set1.polynum[k + 0] = polynum;
		ss_set1.polynum[k + 1] = polynum;
#endif

		int size1 = n;

		// create set2, an enumeration of (+/-t - b)(a)^-1 mod p
		// for b in the set [+/-Bs/2+1 +/- Bs/2+2 ... + Bs]
		n = (poly->s - 1) - n;

		int size2 = n;
		r1 = dconf->firstroot2[i];

		// enumerate a sorted set2 using Kleinjung's shift & merge method.
		int n0 = size1;

		ss_set2.root[0] = r1;
		ss_set2.polynum[0] = 0;

		ss_set2.root[1] = ss_set2.root[0] + rootupdates[(n0)*bound + i];
		ss_set2.polynum[1] = ss_set2.polynum[0] + (1 << n0);

		if (ss_set2.root[1] >= prime)
		{
			ss_set2.root[1] = ss_set2.root[1] - prime;
			int rtmp = ss_set2.root[0];
			int ptmp = ss_set2.polynum[0];
			ss_set2.root[0] = ss_set2.root[1];
			ss_set2.polynum[0] = ss_set2.polynum[1];
			ss_set2.root[1] = rtmp;
			ss_set2.polynum[1] = ptmp;
		}

		for (ii = n0 + 1; ii < poly->s - 1; ii++)
		{
			shift(&ss_set2, &ss_set_t1, &ss_set_t2, &ss_set_t3, ii - n0, n0,
				rootupdates[(ii)*bound + i], prime, 1);
		}

#ifdef SS_TIMING
		gettimeofday(&stop, NULL);
		t_enum_roots += ytools_difftime(&start, &stop);

		start.tv_sec = stop.tv_sec;
		start.tv_usec = stop.tv_usec;
#endif



		// admissible sieving events occur when r1 + r2 < interval
		// or p <= r1 + r2 < p + interval
		// which is equivalent to 
		// -r0 <= r1 < interval - r0 and
		// p - r0 <= r1 < p + interval - r0
		// from here there are two cases:
		// 1) where r0 >= interval, then the first inequality is impossible
		// and we only need to check the second
		// 2) r0 < interval and both inequalities are possible
		// for case 1) do two binary searches on the sorted set2 to find the
		// first root r1' such that r1' >= p - r0, and find the last r1'' such
		// that r1'' < p + interval - r0.  All r1 between these points are
		// sieving events.
		// the 2nd search could probably be replaced by a comparison during
		// a loop starting at r1', especially if that loop is vectorized.
		// (indeed, this is faster.)
		// for case 2) the first substep just iterates through set2 starting
		// from the beginning and stops once r1'' >= interval - r0.
		// The second substep finds by binary search the element r1' >= p - r0
		// and then from there iterates to the end of the set2 array.

		int num_set2 = 1 << size2;

#if 1 //ndef SORT_SET1
		// commence matching
		int max_id = num_set2;
		int sieve_sz = dconf->ss_sieve_sz;

		for (ii = 0; ii < (2 << size1); ii++)
		{
			root1 = ss_set1.root[ii];
			int poly1 = ss_set1.polynum[ii];

			// case 2a only occurs if the initial root is really small
			if (root1 < interval)
			{
				int startid = 0;
				root2 = ss_set2.root[startid];

				// case 2a: beginning of set2 array matching
				while (((root1 + root2) < interval) && (startid < num_set2))
				{
					int idx = ss_set2.polynum[startid] + poly1;

#ifdef FULL_INTERVAL_NOTATION
					int sum = root1 + root2 + interval;
					int sign = sum >= interval;

					if (sign == 0) sum = interval - sum;
#else
					int sum = root1 + root2;
#endif

					if (resieve)
					{
						int report_num;
						uint8_t sieve_val;

						sieve_val = dconf->ss_sieve_p[idx * sieve_sz + sum];

						if (sieve_val < 0xff)
						{

							report_num = dconf->report_ht_p[hash64(idx * sieve_sz + sum) %
								dconf->report_ht_size];

							if (report_num == 0)
							{
								invalid_report_count++;
								startid++;
								root2 = ss_set2.root[startid];
								continue;
							}

							int smooth_num = dconf->smooth_num[report_num];
							while (mpz_tdiv_ui(dconf->Qvals[report_num], prime) == 0)
							{
								dconf->fb_offsets[report_num][smooth_num++] = i;
								mpz_tdiv_q_ui(dconf->Qvals[report_num],
									dconf->Qvals[report_num], prime);
							}
							dconf->smooth_num[report_num] = smooth_num;
						}
					}
					else
					{
						dconf->ss_sieve_p[idx * sieve_sz + sum] -= logp;
					}

					startid++;
					root2 = ss_set2.root[startid];

#ifdef SS_TIMING
					nummatch++;
#endif
				}
			}

			// this case always occurs, either as case 1
			// or as case 2b
			{
				int startid;
#ifdef SORT_SET1

				// if set1 is sorted then we can gradually decrease the max index
				// of the search.
				startid = binary_find_ge(ss_set2.root, max_id, prime - root1 - interval);
				max_id = startid;
#else
				startid = binary_find_ge(ss_set2.root, num_set2, prime - root1 - interval);
#endif

				root2 = ss_set2.root[startid];

#ifndef USE_POLY_BUCKET_PN_COMBINED_VARIATION
				root1 = root1 - prime + interval;
#else
				root1 = root1 - prime;
#endif

				// case 1: middle of set2 array matching
#ifndef USE_POLY_BUCKET_PN_COMBINED_VARIATION
				while (((startid) < num_set2) &&
					((root1 + root2) < interval2))
#else
				while (((startid) < num_set2) &&
					((root1 + root2) < interval))
#endif
				{
					int idx = ss_set2.polynum[startid] + poly1;

#ifndef USE_POLY_BUCKET_PN_COMBINED_VARIATION
#ifdef USE_QSORT
					int sum1 = root1 + ss_set2[startid].root;
#else
					int sum1 = root1 + ss_set2.root[startid];
#endif

					int sign1 = sum1 >= interval;

					// the smoothness checking code all assumes -side hits
					// are arranged with respect to the 0 point, symmetric 
					// to +side hits.  so if this is a -side hit we need
					// to flip its orientation relative to the 0 point so
					// that smoothness checking works.  Note that the
					// exact 0 point is ignored for now because it raises
					// errors in trial division... need to look into how to
					// properly handle hits at sum = interval.
					if (sign1 == 0) sum1 = interval - sum1;

#else
					int sum1 = root1 + root2;
					int sign = sum1 >= 0;
					if (sign == 0)
					{
						sum1 = 0 - sum1;
					}
#endif

					if (resieve)
					{
						if (resieve)
						{
							int report_num;
							uint8_t sieve_val;

							if (sign)
							{
								sieve_val = dconf->ss_sieve_p[idx * sieve_sz + sum1];
							}
							else
							{
								sieve_val = dconf->ss_sieve_n[idx * sieve_sz + sum1];
							}

							if (sieve_val < 0xff)
							{
								if (sign)
								{
									report_num = dconf->report_ht_p[hash64(idx * sieve_sz + sum1) %
										dconf->report_ht_size];
								}
								else
								{
									report_num = dconf->report_ht_n[hash64(idx * sieve_sz + sum1) %
										dconf->report_ht_size];
								}

								if (report_num == 0)
								{
									invalid_report_count++;
									startid++;
									root2 = ss_set2.root[startid];
									continue;
								}

								int smooth_num = dconf->smooth_num[report_num];
								while (mpz_tdiv_ui(dconf->Qvals[report_num], prime) == 0)
								{
									dconf->fb_offsets[report_num][smooth_num++] = i;
									mpz_tdiv_q_ui(dconf->Qvals[report_num],
										dconf->Qvals[report_num], prime);
								}

								dconf->smooth_num[report_num] = smooth_num;
							}
						}
					}
					else
					{
						if (sign)
						{
							dconf->ss_sieve_p[idx * sieve_sz + sum1] -= logp;
						}
						else
						{
							dconf->ss_sieve_n[idx * sieve_sz + sum1] -= logp;
						}
					}

					startid++;
					root2 = ss_set2.root[startid];

#ifdef SS_TIMING
					nummatch++;
#endif
				}

			}
		}

		//exit(1);

#else


		__m512i vpid = _mm512_set1_epi32(pid);
		__m512i vprime = _mm512_set1_epi32(prime);
		__m512i vinterval = _mm512_set1_epi32(interval);
		__m512i vz = _mm512_setzero_epi32();
		__m512i vmaxroot = _mm512_set1_epi32(SS_MAX_ROOT);
		// with SIMD we can match much faster.  It doesn't find some hits compared
		// to non-SIMD, so there must be a bug somewhere.  But it doesn't really
		// matter because even with faster matching subset-sum still isn't
		// better than normal.
		int max_id = num_set2;
		for (ii = 0; ii < (2 << size1); ii += 16)
		{
			// this case always occurs, either as case 1
			// or as case 2b

			int startid;

			startid = binary_find_ge(ss_set2.root, max_id,
				prime - ss_set1.root[ii + 15] - interval);
			max_id = startid;

			__m512i vr0 = _mm512_load_epi32(ss_set1.root + ii);
			__m512i vr1 = _mm512_sub_epi32(vr0, vprime);
			__m512i vr2 = _mm512_set1_epi32(ss_set2.root[startid]);
			__m512i vp1 = _mm512_load_epi32(ss_set1.polynum + ii);
			__m512i vp2 = _mm512_set1_epi32(ss_set2.polynum[startid]);
			__mmask16 hits;

			// case1 and case2b
			while (((startid < num_set2)) &&
				((ss_set1.root[ii + 0] - prime + ss_set2.root[startid]) < interval))

			{
				vp2 = _mm512_add_epi32(vp1, vp2);
				vr2 = _mm512_add_epi32(vr1, vr2);
				hits = _mm512_cmplt_epi32_mask(_mm512_abs_epi32(vr2), vinterval);

				int numhits = _mm_popcnt_u32(hits);

				if (numhits > 0)
				{
					__mmask16 signs = _mm512_cmplt_epi32_mask(vr2, vz);
					vr2 = _mm512_mask_sub_epi32(vr2, signs, vz, vr2);
					vr2 = _mm512_mask_or_epi32(vr2, signs, vr2, vmaxroot);
					vr2 = _mm512_or_epi32(vr2, vpid);

					__m512i vbsz = _mm512_mask_i32gather_epi32(vbsz, hits, vp2, psize_ptr, 4);
					__m512i vindex = _mm512_slli_epi32(vp2, 11);
					vindex = _mm512_add_epi32(vindex, vbsz);

					_mm512_mask_i32scatter_epi32(pslice_ptr, hits, vindex, vr2, 4);
					vbsz = _mm512_add_epi32(vbsz, _mm512_set1_epi32(1));
					_mm512_mask_i32scatter_epi32(psize_ptr, hits, vp2, vbsz, 4);

					nummatch += numhits;
				}

				startid++;
				vr2 = _mm512_set1_epi32(ss_set2.root[startid]);
				vp2 = _mm512_set1_epi32(ss_set2.polynum[startid]);

			}
		}

		// case 2a only occurs if the initial root is really small.
		// with sorted root1s, we can stop looking as soon as the
		// root comparison fails.
		for (ii = 0; ii < (2 << size1); ii++)
		{
			root1 = ss_set1.root[ii];
			if (root1 < interval)
			{
				int startid = 0;
				root2 = ss_set2.root[startid];

				// case 2a: beginning of set2 array matching
				while ((root1 + root2) < interval) // ((startid < num_set2) && 
				{
					int polysum = ss_set2.polynum[startid] + ss_set1.polynum[ii];
					int idx = polysum;
					int sum = root1 + root2;

					uint32_t bsz = psize_ptr[idx];
					pslice_ptr[idx * bucketalloc + bsz] = (pid | sum);
					psize_ptr[idx]++;

#ifdef ENABLE_BUCKET_PROTECTION
					if (psize_ptr[idx] >= dconf->ss_slices_p[slice].alloc)
					{
						printf("error bucket overflow: size of poly bucket %d = %d\n",
							idx, psize_ptr[idx]);
						exit(1);
					}
#endif

					startid++;
					root2 = ss_set2.root[startid];

#ifdef SS_TIMING
					nummatch++;
#endif
				}
			}
			else
			{
				break;
			}
		}

		//exit(1);

#endif


#ifdef SS_TIMING
		gettimeofday(&stop, NULL);
		t_match_roots += ytools_difftime(&start, &stop);
#endif
	}

#ifdef SS_TIMING
	printf("ran subset-sum on %d primes\n", nump);
	printf("found %d sieve hits\n", nummatch);
	printf("enumerating roots: %1.4f seconds\n", t_enum_roots);
	printf("matching roots: %1.4f seconds\n", t_match_roots);
#endif


	free(ss_set1.root);
	free(ss_set2.root);
	free(ss_tmp1.root);
	free(ss_tmp2.root);
	free(ss_set1.polynum);
	free(ss_set2.polynum);
	free(ss_tmp1.polynum);
	free(ss_tmp2.polynum);
	free(ss_set_t1.root);
	free(ss_set_t2.root);
	free(ss_set_t3.root);
	free(ss_set_t1.polynum);
	free(ss_set_t2.polynum);
	free(ss_set_t3.polynum);

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
	int using_ss_search = 0;

	if ((sconf->factor_base->ss_start_B > 0) &&
		(sconf->factor_base->ss_start_B < sconf->factor_base->B))
	{
		using_ss_search = 1;
	}

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

#ifdef USE_SS_SEARCH

		if (using_ss_search && (i >= sconf->factor_base->ss_start_B))
		{
			// first poly b: sum of first n positive Bl
			bmodp = (int)mpz_tdiv_ui(dconf->polyb1, prime);

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
			root1 = (int)(modsqrt[i]) - bmodp;
			root2 = (int)(prime - modsqrt[i]) - bmodp;
			if (root1 < 0) root1 += prime;
			if (root2 < 0) root2 += prime;

			dconf->firstroot1a[i] = (int)((uint64_t)inv * (uint64_t)root1 % (uint64_t)prime);
			dconf->firstroot1b[i] = (int)((uint64_t)inv * (uint64_t)root2 % (uint64_t)prime);

			bmodp = (int)mpz_tdiv_ui(dconf->polyb2, prime);

			// - (B1 + B2 + B3 + ... Bs-1)
			bmodp = prime - bmodp;

			dconf->firstroot2[i] = (uint32_t)((uint64_t)inv * (uint64_t)bmodp % (uint64_t)prime);
		}
#endif
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
        CHECK_NEW_SLICE(i);

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

#ifdef USE_SS_SEARCH

		if (using_ss_search && (i >= sconf->factor_base->ss_start_B))
		{
			// first poly b: sum of first n positive Bl
			bmodp = (int)mpz_tdiv_ui(dconf->polyb1, prime);

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
			root1 = (int)(modsqrt[i]) - bmodp;
			root2 = (int)(prime - modsqrt[i]) - bmodp;
			if (root1 < 0) root1 += prime;
			if (root2 < 0) root2 += prime;

			dconf->firstroot1a[i] = (int)((uint64_t)inv * (uint64_t)root1 % (uint64_t)prime);
			dconf->firstroot1b[i] = (int)((uint64_t)inv * (uint64_t)root2 % (uint64_t)prime);

			bmodp = (int)mpz_tdiv_ui(dconf->polyb2, prime);

			// - (B1 + B2 + B3 + ... Bs-1)
			bmodp = prime - bmodp;

			dconf->firstroot2[i] = (uint32_t)((uint64_t)inv * (uint64_t)bmodp % (uint64_t)prime);
		}
#endif

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

#ifdef USE_SS_SEARCH

		if (using_ss_search && (i >= sconf->factor_base->ss_start_B))
		{
			// first poly b: sum of first n positive Bl
			bmodp = (int)mpz_tdiv_ui(dconf->polyb1, prime);

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
			root1 = (int)(modsqrt[i]) - bmodp;
			root2 = (int)(prime - modsqrt[i]) - bmodp;
			if (root1 < 0) root1 += prime;
			if (root2 < 0) root2 += prime;

			dconf->firstroot1a[i] = (int)((uint64_t)inv * (uint64_t)root1 % (uint64_t)prime);
			dconf->firstroot1b[i] = (int)((uint64_t)inv * (uint64_t)root2 % (uint64_t)prime);

			bmodp = (int)mpz_tdiv_ui(dconf->polyb2, prime);

			// - (B1 + B2 + B3 + ... Bs-1)
			bmodp = prime - bmodp;

			dconf->firstroot2[i] = (uint32_t)((uint64_t)inv * (uint64_t)bmodp % (uint64_t)prime);
		}
#endif

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

