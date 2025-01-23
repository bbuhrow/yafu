/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: sieve_core.c 963 2014-05-03 02:04:54Z jasonp_sf $
--------------------------------------------------------------------*/

#include <ms_common.h>
#include "mpqs.h"

#if BLOCK_KB == 32
#define SIEVE_BLOCK_SIZE 32768
#define LOG2_SIEVE_BLOCK_SIZE 15
#elif BLOCK_KB == 64
#define SIEVE_BLOCK_SIZE 65536
#define LOG2_SIEVE_BLOCK_SIZE 16
#else
#error "unsupported sieve block size"
#endif

/*--------------------------------------------------------------------*/
static void add_to_hashtable(bucket_t *entry,
			     uint32 sieve_offset, 
			     uint32 mask, 
			     uint32 prime_index, 
			     uint32 logprime) {

	/* add a 'sieve update' to the hashtable bin pointed
	   to by 'entry'. Hashing is based on the top few
	   bits of sieve_offset, and the range of sieve values
	   represented by one hash bin is given by 'mask'.

	   Note that after the first polynomial it's unlikely
	   that any hash bin will need to grow in size. */

	bucket_entry_t *new_entry;
	uint32 i = entry->num_used++;

	if (i == entry->num_alloc) {
		entry->num_alloc = 2 * i;
		entry->list = (bucket_entry_t *)xrealloc(entry->list,
						2 * i * sizeof(bucket_entry_t));
	}
	new_entry = entry->list + i;
	new_entry->logprime = logprime;
	new_entry->prime_index = prime_index;
	new_entry->sieve_offset = sieve_offset & mask;
}
	
/*--------------------------------------------------------------------*/
static void fill_sieve_block(sieve_conf_t *conf,
			     bucket_t *hash_bucket) {

	/* Routine for conventional sieving
	   
	   Each successive call to fill_sieve_block() reuses the root
	   values that were computed in the previous call */

	uint32 i;
	uint8 *sieve_array = conf->sieve_array;
	uint32 sieve_small_fb_start = conf->sieve_small_fb_start;
	uint32 sieve_large_fb_start = conf->sieve_large_fb_start;
	packed_fb_t *packed_fb = conf->packed_fb;
	bucket_entry_t *list;
	uint32 num_large;

	/* First update the sieve block with values corresponding
	   to the small primes in the factor base. Assuming the
	   entire block will fit in cache, each factor base prime
	   will update many locations within the block, and this
	   phase will run very fast */

	TIME1(sieve_small_time)

	for (i = sieve_small_fb_start; i < sieve_large_fb_start; i++) {
		packed_fb_t *pfbptr = packed_fb + i;
		uint32 prime = pfbptr->prime;
		uint32 root1 = pfbptr->next_loc1;
		uint32 root2 = pfbptr->next_loc2;
		uint8 logprime = pfbptr->logprime;

		/* We assume that the starting positions for each 
		   root (i.e. 'next_loc1' and 'next_loc2') are 
		   stored in ascending order; thus, if next_loc2 
		   is within the sieve block then next_loc1 
		   automatically is too. This allows both root 
		   updates to be collapsed into a single loop, 
		   effectively unrolling by two. 
		   
		   Note that no exception is made for sieve roots
		   that are invalid; they skip the loop below
		   since INVALID_ROOT is a large number */

		while (root2 < SIEVE_BLOCK_SIZE) {
			sieve_array[root1] -= logprime;
			sieve_array[root2] -= logprime;
			root1 += prime;
			root2 += prime;
		}

		/* once next_loc2 exceeds the sieve block size,
		   next_loc1 may still have one update left. In
		   that case, perform the update and switch 
		   next_loc1 and next_loc2, since next_loc1 is now
		   the larger root */

		if (root1 < SIEVE_BLOCK_SIZE) {
			uint32 tmp = root2;
			sieve_array[root1] -= logprime;
			root2 = root1 + prime;
			root1 = tmp;
		}
		pfbptr->next_loc1 = root1 - SIEVE_BLOCK_SIZE;
		pfbptr->next_loc2 = root2 - SIEVE_BLOCK_SIZE;
	}

	TIME2(sieve_small_time)

	/* Now update the sieve block with the rest of the
	   factor base. All of the offsets to update have
	   previously been collected into a hash table.
	   While the updates are to random offsets in the
	   sieve block, the whole block was cached by the
	   previous loop, and memory access to the hashtable
	   entry is predictable and can be prefetched */

	TIME1(sieve_large_time)
	list = hash_bucket->list;
	num_large = hash_bucket->num_used;

	for (i = 0; i < (num_large & (uint32)(~7)); i += 8) {
#ifdef MANUAL_PREFETCH
		PREFETCH(list + i + 16);
#endif
		sieve_array[list[i  ].sieve_offset] -= list[i  ].logprime;
		sieve_array[list[i+1].sieve_offset] -= list[i+1].logprime;
		sieve_array[list[i+2].sieve_offset] -= list[i+2].logprime;
		sieve_array[list[i+3].sieve_offset] -= list[i+3].logprime;
		sieve_array[list[i+4].sieve_offset] -= list[i+4].logprime;
		sieve_array[list[i+5].sieve_offset] -= list[i+5].logprime;
		sieve_array[list[i+6].sieve_offset] -= list[i+6].logprime;
		sieve_array[list[i+7].sieve_offset] -= list[i+7].logprime;
	}
	for (; i < num_large; i++) {
		sieve_array[list[i].sieve_offset] -= list[i].logprime;
	}
	TIME2(sieve_large_time)
}

/*--------------------------------------------------------------------*/
#define PACKED_SIEVE_MASK ((uint64)0x80808080 << 32 | 0x80808080)

#if (defined(GCC_ASM32X) || defined(GCC_ASM64X)) && defined(NDEBUG)
	#if defined(HAS_SSE2)
		#define SCAN_SSE2
	#elif defined(HAS_AMD_MMX) || defined(HAS_SSE)
		#define SCAN_MMX
	#endif
#endif

static uint32 scan_sieve_block(sieve_conf_t *conf,
				mp_t *a, signed_mp_t *b, signed_mp_t *c,
				int32 block_start,
				uint32 cutoff1,
				uint32 poly_index,
				bucket_t *hashtable) {
				
	uint32 i, j;
	uint8 *sieve_array = conf->sieve_array;
	uint64 *packed_sieve = (uint64 *)conf->sieve_array;
	uint32 relations_found = 0;

	TIME1(tf_plus_scan_time)
	for (i = 0; i < SIEVE_BLOCK_SIZE / 8; i += 8) {

		/* test 64 sieve values at a time for large
		   logarithms. The sieve code actually
		   decrements by the logarithms of primes, so
		   if a value is smooth its top bit will be
		   set. We can easily test for this condition
		  
		   Note that for really big jobs the cutoff
		   itself will exceed 128, so that the test
		   below will have many false alarms */

#if defined(SCAN_SSE2)
		uint32 compare_result = 0;
		asm volatile (
			"movdqa (%1), %%xmm0          \n\t"
			"orpd 16(%1), %%xmm0          \n\t"
			"orpd 32(%1), %%xmm0          \n\t"
			"orpd 48(%1), %%xmm0          \n\t"
			"pmovmskb %%xmm0, %0          \n\t"
			: "=r"(compare_result)
			: "r"(packed_sieve + i), "0"(compare_result)
			: "%xmm0");

		if (compare_result == 0)
			continue;

#elif defined(SCAN_MMX)
		uint32 compare_result = 0;
		asm volatile (
			"movq (%1), %%mm0         \n\t"
			"por 8(%1), %%mm0        \n\t"
			"por 16(%1), %%mm0        \n\t"
			"por 24(%1), %%mm0        \n\t"
			"por 32(%1), %%mm0        \n\t"
			"por 40(%1), %%mm0        \n\t"
			"por 48(%1), %%mm0        \n\t"
			"por 56(%1), %%mm0        \n\t"
			"pmovmskb %%mm0, %0       \n\t"
			: "=r"(compare_result)
			: "r"(packed_sieve + i), "0"(compare_result)
			: "%mm0");

		if (compare_result == 0)
			continue;

		/* make it safe to perform floating point
		   again (required, since SQUFOF uses it) */

		asm volatile("emms");

#else
		if (((packed_sieve[i] | packed_sieve[i+1] |
		      packed_sieve[i+2] | packed_sieve[i+3] |
		      packed_sieve[i+4] | packed_sieve[i+5] |
		      packed_sieve[i+6] | packed_sieve[i+7]) & 
		      PACKED_SIEVE_MASK) == (uint64)(0))
			continue;
#endif
	
		/* one or more of the 64 values is probably
		   smooth. Test one at a time and attempt
		   to factor them */

		for (j = 0; j < 64; j++) {
			uint32 bits = sieve_array[8 * i + j];
			if (bits > cutoff1) {
				TIME1(tf_total_time)
				relations_found += 
					check_sieve_val(conf, 
						block_start + (int32)(8*i+j), 
						cutoff1 + 257 - bits,
						a, b, c, poly_index,
						hashtable);
				TIME2(tf_total_time)
			}
		}
	}
	TIME2(tf_plus_scan_time)

#if defined(SCAN_MMX)
	asm volatile("emms");
#endif

	return relations_found;
}

/*--------------------------------------------------------------------*/

uint32 ROUTINE_NAME(sieve_conf_t *conf, 
		  uint32 poly_start, uint32 num_poly) {

	/* The core of the sieving code. This routine performs
	   the sieving and relation collection for 'num_poly'
	   polynomials at a time.

	   Because this is the self-initializing quadratic sieve,
	   groups of polynomials are related. There are 'total_poly'
	   polynomials that share the same 'a' value; we number them
	   from 0 to total_poly-1. poly_start refers to the index
	   of the initial polynomial of the group of num_poly. 

	   Ordinarily, you initialize one polynomial, sieve with it,
	   initialize the next polynomial, sieve with it, etc. However,
	   we want the sieve interval to be small, and when it's very
	   small but the factor base is large we will spend a lot
	   of time initializing polynomials and only a little time 
	   sieving. Self-initialization improves that a lot, but when
	   the factor base exceeds the cache size, even for previous
	   versions of this program, polynomial initialization would
	   take 50% of the total sieving stage runtime when the sieve
	   interval is small. We could make the interval larger, but
	   that means smooth relations occur less often; sieve time
	   does not scale well but poly initialization time does.

	   This routine takes a novel approach: for small factor base
	   primes it sieves one polynomial at a time. However, for
	   all the other factor base primes (90+% of the factor base,
	   in a big factorization) it does not sieve at all, but fills
	   a hashtable with the locations of which prime updates which
	   sieve offset. There is a separate hashtable for each polynomial,
	   and a separate hashtable for positive and negative sieve
	   values. Because these hashtables do not interact with each other,
	   we can do all the sieving for a few factor base primes in one
	   polynomial, then switch to the next polynomial and do the 
	   sieving for the same primes. For these few primes all of the 
	   sieving for one poly happens at once, and so the roots of these 
	   primes can then be updated to become roots of the next polynomial.

	   The upshot is that we deal with a block of primes at a time,
	   and for that block do all the sieving for num_poly polynomials;
	   the roots for each prime get updated whenever sieving for a
	   polynomial has finished. This in effect performs a tiling of
	   the sieve interval, the factor base, and the auxiliary information
	   needed to switch polynomials. When the tile size of these three
	   quantities is chosen small enough, everything fits in cache
	   for the duration of the sieving.

	   Note that putting all the updates into a hashtable uses up
	   much more space than sieving with it; however, when filling up
	   a hashtable linearly only one cache line per hash bin is active
	   at any given time, and if the number of bins is small enough
	   the processor can fit all the updates into store buffers as
	   they happen.  Contrast that with filling up a conventional 
	   sieve interval, where all the cache lines in the interval
	   are active.

	   All of the hashtables are filled in first; then we proceed one
	   polynomial at a time. For each polynomial we sieve with the
	   small factor base primes, add in the sieve values from the 
	   appropriate set of hash bins (the bin size matches the sieve
	   block size for this implementation), then do trial factoring
	   on presumed smooth sieve values. There's one more perk to doing
	   things this way: since a hash bin can be made to contain the
	   actual primes as well as the log values of those primes, 
	   trial factoring can be sped up an order of magnitude by reusing
	   the hash bins to detect primes that divide sieve values.

	   In short: black magic, awful mess, really fast. */

	uint32 i, j, k, m, n;
	uint32 relations_found = 0;
	uint32 num_sieve_blocks = conf->num_sieve_blocks;
	uint32 sieve_size;
	uint32 fb_size = conf->fb_size;
	fb_t *factor_base = conf->factor_base;
	uint32 num_factors = conf->num_poly_factors;
	uint32 total_poly = 1 << (num_factors - 1);
	bucket_t *buckets;
	uint32 poly_index;
	uint32 *poly_b_array;

	/* all hash bins start off empty */

	for (i = 0; i < num_poly * num_sieve_blocks; i++) {
		conf->buckets[i].num_used = 0;
	}

	poly_b_array = conf->poly_b_array;
	sieve_size = num_sieve_blocks * SIEVE_BLOCK_SIZE;
	i = conf->sieve_large_fb_start;

	/* for each block of factor base primes */

	while (i < fb_size) {
		uint32 fb_block = MIN(conf->fb_block, fb_size - i);
		fb_t *fb_start = factor_base + i;

		buckets = conf->buckets;
		poly_index = poly_start;

		/* for polynomial j */

		for (j = 0; j < num_poly; j++) {

			uint32 next_action;
			uint32 *poly_b_start;

			TIME1(bucket_time)
			for (k = 0; k < fb_block; k++) {
				fb_t *fbptr = fb_start + k;
				uint32 prime = fbptr->prime;
				uint32 root1, root2;
				uint8 logprime = fbptr->logprime;

#ifdef MANUAL_PREFETCH
				if (k % 4 == 0)
					PREFETCH(fbptr + 4);
#endif

				/* grab the two roots of each prime in
				   the block, and do all the sieving. For each
				   sieve value, use the upper bits of the 
				   value to determine the hash bin of 
				   polynomial j to update */

				root1 = fbptr->root1;
				while (root1 < sieve_size) {
					add_to_hashtable(buckets + 
					       	(root1>>LOG2_SIEVE_BLOCK_SIZE),
						root1, SIEVE_BLOCK_SIZE - 1, 
						i + k, logprime);
					root1 += prime;
				}
	
				root2 = fbptr->root2;
				while (root2 < sieve_size) {
					add_to_hashtable(buckets + 
					       	(root2>>LOG2_SIEVE_BLOCK_SIZE),
						root2, SIEVE_BLOCK_SIZE - 1, 
						i + k, logprime);
					root2 += prime;
				}
			}
			TIME2(bucket_time)

			/* this block has finished sieving for polynomial
			   j; now select the new roots for polynomial j+1.
			   Note that there is no polynomial with index
			   'total_poly', so no update is needed after the
			   last polynomial has updated its hashtable */

			if (++poly_index == total_poly)
				continue;

			k = poly_index;
			next_action = conf->next_poly_action[k];
			poly_b_start = poly_b_array;
			n = next_action & 0x7f;
			
			/* Update the two roots for each prime in the block.
			   Do not sort them in ascending order; the loop 
			   above does not require it, and the unpredictable 
			   branch required to do so takes a large fraction 
			   of the total poly initialization time! */

			TIME1(next_poly_large_time)
			if (next_action & 0x80) {
				for (k = 0; k < fb_block; 
					k++, poly_b_start += num_factors) {
		
					fb_t *fbptr = fb_start + k;
					uint32 prime = fbptr->prime;
					uint32 root1 = fbptr->root1;
					uint32 root2 = fbptr->root2;
		
					m = poly_b_start[n];
					root1 = mp_modadd_1(root1, m, prime);
					root2 = mp_modadd_1(root2, m, prime);
					fbptr->root1 = root1;
					fbptr->root2 = root2;
				}
			}
			else {
				for (k = 0; k < fb_block; 
					k++, poly_b_start += num_factors) {
		
					fb_t *fbptr = fb_start + k;
					uint32 prime = fbptr->prime;
					uint32 root1 = fbptr->root1;
					uint32 root2 = fbptr->root2;
		
					m = poly_b_start[n];
					root1 = mp_modsub_1(root1, m, prime);
					root2 = mp_modsub_1(root2, m, prime);
					fbptr->root1 = root1;
					fbptr->root2 = root2;
				}
			}
			TIME2(next_poly_large_time)

			/* polynomial j is finished; point to the
			   hash bins for polynomial j+1 */

			buckets += num_sieve_blocks;
		}

		/* current block has updated all hashtables;
		   go on to the next block */

		i += fb_block;
		poly_b_array += fb_block * num_factors;
	}

	/* All of the hashtables have been filled; now proceed
	   one polynomial at a time and sieve with the small
	   factor base primes */

	buckets = conf->buckets;
	poly_index = poly_start;

	/* for each polynomial */

	for (i = 0; i < num_poly; i++) {

		mp_t *a;
		signed_mp_t *b;
		signed_mp_t c, tmp, tmp_n;
		uint32 cutoff1;
		uint32 block;
		packed_fb_t *packed_fb = conf->packed_fb;
		uint32 next_action;

		/* make temporary copies of all of the roots; sieving
		   cannot take place all at once like the large FB primes
		   got to do, so we must update the temporary roots as
		   we go from sieve block to sieve block */

		TIME1(plus_init_time)
		for (j = MIN_FB_OFFSET + 1; 
				j < conf->sieve_large_fb_start; j++) {
			fb_t *fbptr = factor_base + j;
			packed_fb_t *pfbptr = packed_fb + j;
			pfbptr->prime = fbptr->prime;
			pfbptr->logprime = fbptr->logprime;
			pfbptr->next_loc1 = fbptr->root1;
			pfbptr->next_loc2 = fbptr->root2;
		}

		/* Form 'c', equal to (b * b - n) / a (exact
		   division). Having this quantity available
		   makes trial factoring a little faster.

		   Also double the 'b' value, since nobody
		   else uses it and the trial factoring code
		   prefers 2*b */

		mp_copy(conf->n, &tmp_n.num);
		tmp_n.sign = POSITIVE;

		a = &conf->curr_a;
		b = conf->curr_b + poly_start + i;
		signed_mp_mul(b, b, &tmp);
		signed_mp_add(b, b, b);
		signed_mp_sub(&tmp, &tmp_n, &tmp);
		mp_div(&tmp.num, a, &c.num);
		c.sign = tmp.sign;

		/* choose the first sieving cutoff */

		cutoff1 = mp_bits(&c.num);
		if (cutoff1 >= conf->cutoff1)
			cutoff1 -= conf->cutoff1;
		else
			cutoff1 = 0;

		TIME2(plus_init_time)

		/* for each sieve block */

		for (block = 0; block < num_sieve_blocks; block++) {
	
			int32 block_start = -(sieve_size / 2) + 
					(block * SIEVE_BLOCK_SIZE);
	
			/* initialize the sieve array */

			memset(conf->sieve_array, (int8)(cutoff1 - 1), 
				(size_t)SIEVE_BLOCK_SIZE);
	
			/* do the sieving and add in the values from the
			   hash bin corresponding to this sieve block */

			fill_sieve_block(conf, buckets + block);
	
			
			/* trial factor sieve values whose logarithm is
			   sufficiently large */

			relations_found += scan_sieve_block(conf, a, b, &c,
							block_start, cutoff1,
							poly_index, 
							buckets + block);
		}

		/* update the roots for each of the small factor
		   base primes to correspond to the next polynomial */

		if (++poly_index == total_poly)
			break;

		k = poly_index;
		next_action = conf->next_poly_action[k];
		poly_b_array = conf->poly_b_small[next_action & 0x7f];
		k = conf->sieve_large_fb_start;

		TIME1(next_poly_small_time)
		if (next_action & 0x80) {
			for (j = MIN_FB_OFFSET + 1; j < k; j++) {
	
				fb_t *fbptr = factor_base + j;
				uint32 prime = fbptr->prime;
				uint32 root1 = fbptr->root1;
				uint32 root2 = fbptr->root2;
	
				if (root1 == INVALID_ROOT)
					continue;
	
				m = poly_b_array[j];
				root1 = mp_modadd_1(root1, m, prime);
				root2 = mp_modadd_1(root2, m, prime);
				if (root2 > root1) {
					fbptr->root1 = root1;
					fbptr->root2 = root2;
				}
				else {
					fbptr->root1 = root2;
					fbptr->root2 = root1;
				}
			}
		}
		else {
			for (j = MIN_FB_OFFSET + 1; j < k; j++) {
	
				fb_t *fbptr = factor_base + j;
				uint32 prime = fbptr->prime;
				uint32 root1 = fbptr->root1;
				uint32 root2 = fbptr->root2;
	
				if (root1 == INVALID_ROOT)
					continue;
	
				m = poly_b_array[j];
				root1 = mp_modsub_1(root1, m, prime);
				root2 = mp_modsub_1(root2, m, prime);
				if (root2 > root1) {
					fbptr->root1 = root1;
					fbptr->root2 = root2;
				}
				else {
					fbptr->root1 = root2;
					fbptr->root2 = root1;
				}
			}
		}
		TIME2(next_poly_small_time)
		buckets += num_sieve_blocks;
	}

	return relations_found;
}
