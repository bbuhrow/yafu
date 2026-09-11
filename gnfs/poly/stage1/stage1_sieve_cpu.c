/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: stage1_sieve_cpu.c 1084 2026-05-10 03:05:04Z jasonp_sf $
--------------------------------------------------------------------*/

#include <stage1.h>
#include <cpu_intrinsics.h>

/* CPU collision search; this code looks for self-collisions
   among arithmetic progressions, by finding k1 and k2 such that
   for two arithmetic progressions r1+k*p1^2 and r2+k*p2^2 we
   have

      r1 + k1*p1^2 = r2 + k2*p2^2

   such that
      - p1 and p2 are coprime and < 2^P_BITS
      - the value where they coincide is of size smaller
        than a fixed bound

   This code uses an algorithm from Robert Gerbicz to find collisions 
   across all the p1 and p2 in the set simultaneously. The size of the 
   range is comparatively very large, and the algorithm as described 
   above is only efficient for pretty small problems. We cannot practically 
   check all the possibilities because for e.g. 512-bit GNFS the range 
   on k1 and k2 is typically around 10^6 and the set has ~10^6 progressions. 
   We can reduce the memory use by using a hashtable with blocks of 
   size > p_min^2 and filling each block individually, but that only 
   reduces the memory use and not the actual work required. Additionally, 
   p_min^2 very quickly becomes too large for a perfect hash function, 
   so the hashtable would be very nearly empty in most places.
   
   To scale up the idea, we further use a 'special-q' formulation 
   where all the inputs to the sorting are constrained to fall 
   on a third arithmetic progression r3 + k*q^2 for some k. We choose 
   a given q and for each of its roots run the complete collision
   search; forcing all the p1 and p2 to lie on the q progression 
   limits the size of numbers in the colision search to max(p)^2. This is 
   analogous to lattice sieving across the interval.
   
   This allows us to choose q so that the collision search problem is of 
   reasonable size but the collisions found are still over the original, 
   impractically large range. 

   A side effect of this is that the complete range must be < 
   2^(64 + 2 * P_BITS) in size, though this appears to be sufficient 
   for very large problems, e.g. 1024-bit GNFS */

/* #define HAVE_PROF */
#ifdef HAVE_PROF
#define SHOW_PROF __attribute__((noinline))
#else
#define SHOW_PROF /* nothing */
#endif

/*---- packed arrays of arithmetic progressions p --------------------------*/

/* structure for a single progression */

typedef struct {
	uint32 p;
	uint32 num_roots;
	uint64 mont_t;
	uint64 mont_w;
	uint64 mont_r;
	uint64 pp;
	uint64 ss_mod_pp;
	uint64 roots[MAX_ROOTS];
} p_packed_t;

/* the number of uint64 sized elements
   preceding hash_list_t roots in the above structure */

#define P_PACKED_HEADER_WORDS 6

/* structure for a collection of the above structures. The
   storage is in compressed format, where unused roots in the
   above are squeezed out. We only ever expect to iterate
   linearly through the structure so this saves a great deal 
   of memory. */

typedef struct {
	uint32 num_p;
	uint32 num_roots;
	uint32 p_size_alloc;
	uint64 sieve_size;
	uint32 curr_hashiter;
	p_packed_t *curr;
	p_packed_t *packed_array;
} p_packed_var_t;

/* list manipulation functions */

static void 
p_packed_init(p_packed_var_t *s)
{
	memset(s, 0, sizeof(p_packed_var_t));

	s->p_size_alloc = 100;
	s->packed_array = s->curr = (p_packed_t *)xmalloc(s->p_size_alloc *
						sizeof(p_packed_t));
}

static void 
p_packed_free(p_packed_var_t *s)
{
	free(s->packed_array);
}

static void
p_packed_reset(p_packed_var_t *s)
{
	s->num_p = s->num_roots = 0;
	s->curr = s->packed_array;
}

static p_packed_t * 
p_packed_next(p_packed_t *curr)
{
	return (p_packed_t *)((uint64 *)curr + 
			P_PACKED_HEADER_WORDS + curr->num_roots);
}

static SHOW_PROF void 
store_p_packed(uint64 p, uint32 num_roots, mpz_t *roots, void *extra)
{
	/* used to insert a new arithmetic progression and all its
	   roots into the list */

	uint32 i;
	uint64 pp, t;
	p_packed_var_t *s = (p_packed_var_t *)extra;
	p_packed_t *curr;

	if ((void *)s->curr >=
	    (void *)(s->packed_array + s->p_size_alloc - 1)) {
		uint32 *oldptr = (uint32 *)s->packed_array;

		/* we have to be careful here because reallocating
		   the array will cause memory to change out from
		   under s->curr, and we cannot index into the new
		   array because its entries are stored in
		   compressed format */
	   
		s->p_size_alloc *= 2;
		s->packed_array = (p_packed_t *)xrealloc(
						s->packed_array,
						s->p_size_alloc *
						sizeof(p_packed_t));
		s->curr = (p_packed_t *)((uint32 *)s->packed_array +
					((uint32 *)s->curr - oldptr));
	}

	pp = (uint64)p * p;

	curr = s->curr;
	curr->p = (uint32)p;
	curr->pp = pp;
	curr->mont_w = montmul64_w(pp);

	/* form 2^128 mod pp and 2^192 mod pp */
	t = -pp % pp;
	curr->mont_r = mp_modmul_2(t, t, pp);
	curr->mont_t = mp_modmul_2(curr->mont_r, t, pp);

	/* the sieve is centered at 0 so its first offset is
	   at -sieve_size / 2 */
	curr->num_roots = num_roots;
	curr->ss_mod_pp = pp - (s->sieve_size / 2) % pp;
	for (i = 0; i < num_roots; i++)
		curr->roots[i] = gmp2uint64(roots[i]);

	s->num_p++;
	s->num_roots += num_roots;
	s->curr = p_packed_next(s->curr);
}

/*---- packed arrays of special-q ------------------------------------------*/

typedef struct {
	uint64 q;
	uint32 num_roots;
	uint32 pad;
	uint128 roots[MAX_ROOTS];
} q_packed_t;

/* the number of uint64 sized elements
   preceding the roots in the above structure */

#define Q_PACKED_HEADER_WORDS 2

/* structure for a collection of the above structures. The
   storage is in compressed format, where unused roots in the
   above are squeezed out. We only ever expect to iterate
   linearly through the structure so this saves a great deal 
   of memory. */

typedef struct {
	uint32 num_q;
	uint32 num_roots;
	uint32 q_size_alloc;
	q_packed_t *curr;
	q_packed_t *packed_array;
} q_packed_var_t;

/* list manipulation functions */

static void 
q_packed_init(q_packed_var_t *s)
{
	memset(s, 0, sizeof(q_packed_var_t));

	s->q_size_alloc = 100;
	s->packed_array = s->curr = (q_packed_t *)xmalloc(s->q_size_alloc *
						sizeof(q_packed_t));
}

static void 
q_packed_free(q_packed_var_t *s)
{
	free(s->packed_array);
}

static void
q_packed_reset(q_packed_var_t *s)
{
	s->num_q = s->num_roots = 0;
	s->curr = s->packed_array;
}

static q_packed_t * 
q_packed_next(q_packed_t *curr)
{
	return (q_packed_t *)((uint64 *)curr + 
			Q_PACKED_HEADER_WORDS + 2 * curr->num_roots);
}

static SHOW_PROF void 
store_q_packed(uint64 q, uint32 num_roots, mpz_t *roots, void *extra)
{
	/* used to insert a new arithmetic progression and all its
	   roots into the list */

	uint32 i;
	q_packed_var_t *s = (q_packed_var_t *)extra;
	q_packed_t *curr;

	if ((void *)s->curr >=
	    (void *)(s->packed_array + s->q_size_alloc - 1)) {
		uint32 *oldptr = (uint32 *)s->packed_array;

		/* we have to be careful here because reallocating
		   the array will cause memory to change out from
		   under s->curr, and we cannot index into the new
		   array because its entries are stored in
		   compressed format */
	   
		s->q_size_alloc *= 2;
		s->packed_array = (q_packed_t *)xrealloc(
						s->packed_array,
						s->q_size_alloc *
						sizeof(q_packed_t));
		s->curr = (q_packed_t *)((uint32 *)s->packed_array +
					((uint32 *)s->curr - oldptr));
	}

	curr = s->curr;
	curr->q = q;
	curr->pad = 0;
	curr->num_roots = num_roots;
	memset(curr->roots, 0, 
		num_roots * sizeof(curr->roots[0]));
	for (i = 0; i < num_roots; i++)
		mpz_export(curr->roots + i, NULL, -1, 
				sizeof(uint32), 0, 0, roots[i]);

	s->num_q++;
	s->num_roots += num_roots;
	s->curr = q_packed_next(s->curr);
}

/*------------------------------------------------------------------------*/
#define MAX_RADIX_BITS 8

#define BUCKET_BITS 8
#define BUCKET_SIZE (1 << BUCKET_BITS)
#define BUCKET_MASK (BUCKET_SIZE - 1)

typedef struct {
	uint64 key;
	uint32 data;
} match_t;

typedef struct {
	/* initial state */
	uint32 num_sort;
	uint32 num_sort_alloc;
	uint64 *sort_key;
	uint32 *sort_data;

	/* modular inverses */
	uint32 num_inv_alloc;
	uint64 *invtable;

	/* aprog factories */
	p_packed_var_t p_array;
	q_packed_var_t q_array;

	uint64 curr_q;
	uint128 curr_q_root;

	/* state for bucket sorting */
	uint32 num_bucket_alloc;
	uint32 * bucket_data;

	uint32 * bucket_len;
	uint32 * bucket_offset;
	uint32 * bucket_prefix;
	uint32 * bucket_next;

	/* state for pruning false positives */
	uint32 num_dup_alloc;
	uint32 * dup_data;

	uint32 * dup0;
	uint32 * dup1;
	uint32 * dup2;
	uint32 * bin_prefix;

	/* state for putative duplicate keys */
	uint32 num_survivors;
	uint32 num_survivors_alloc;
	uint64 * survivors;
	uint64 * sort_survivors;

	uint32 num_match_alloc;
	match_t * match;
} cpu_thread_data_t;

void
cpu_thread_data_init(void *data, int threadid)
{
	stage1_sieve_data_t *sd = (stage1_sieve_data_t *)data;
	cpu_thread_data_t *t = (cpu_thread_data_t *)xcalloc(1,
					sizeof(cpu_thread_data_t));

	t->num_survivors_alloc = 1000;
	t->survivors = (uint64 *)xmalloc(t->num_survivors_alloc *
					sizeof(uint64));
	t->sort_survivors = (uint64 *)xmalloc(t->num_survivors_alloc *
					sizeof(uint64));
	t->num_match_alloc = 1000;
	t->match = (match_t *)xmalloc(t->num_match_alloc *
					sizeof(match_t));

	/* set up root generation arrays */

	p_packed_init(&t->p_array);
	q_packed_init(&t->q_array);
	sd->threads[threadid].hw_thread_data = t;
}

void
cpu_thread_data_free(void *data, int threadid)
{
	stage1_sieve_data_t *sd = (stage1_sieve_data_t *)data;
	stage1_sieve_thread_data_t *t0 = sd->threads + threadid;
	cpu_thread_data_t *t = (cpu_thread_data_t *)t0->hw_thread_data;

	p_packed_free(&t->p_array);
	q_packed_free(&t->q_array);
	free(t->sort_key);
	free(t->sort_data);
	free(t->invtable);
	if (t->bucket_data)
		aligned_free(t->bucket_data);
	if (t->dup_data)
		aligned_free(t->dup_data);
	free(t->survivors);
	free(t->sort_survivors);
	free(t->match);
	free(t);
}

/*------------------------------------------------------------------------*/
static void 
grow_sort(cpu_thread_data_t * td, uint32 new_size)
{
	td->num_sort_alloc = new_size;
	td->sort_key = (uint64 *)xrealloc(td->sort_key, 
				new_size * sizeof(uint64));
	td->sort_data = (uint32 *)xrealloc(td->sort_data, 
				new_size * sizeof(uint32));
}

/*------------------------------------------------------------------------*/
static SHOW_PROF uint32 
bucket_sort(cpu_thread_data_t * td, uint32 radix_bits)
{
	/* perform a single radix sort pass of size
	   1<<radix_bits; this breaks each radix bin up
	   into buckets of size BUCKET_SIZE and 
	   allocates buckets on the fly, so keys are sorted
	   in a single pass */

	uint32 i;
	uint32 radix_size = 1 << radix_bits;
	uint32 radix_mask = radix_size - 1;
	uint32 num_sort = td->num_sort;
	uint32 num_buckets = radix_size + 1 + num_sort / BUCKET_SIZE;
	uint32 last_bucket;
	uint32 max_bin_size;
	uint64 * key = td->sort_key;
	uint32 * bucket_len;
	uint32 * bucket_offset;
	uint32 * bucket_prefix;
	uint32 * bucket_next;
	uint32 bucket_words;

	/* resize if necessary; allocate a single array and
	   concatenate all the bucket sort state to avoid
	   associative conflicts in L1 */

	bucket_words = 2 * radix_size + num_buckets * (BUCKET_SIZE + 1);
	if (bucket_words > td->num_bucket_alloc) {
		uint32 t = bucket_words * 3 / 2;
		td->num_bucket_alloc = t;
		aligned_free(td->bucket_data);
		td->bucket_data = (uint32 *)aligned_malloc(
					t * sizeof(uint32), 64);
	}
	bucket_len = td->bucket_len = td->bucket_data;
	bucket_offset = td->bucket_offset = bucket_len + radix_size;
	bucket_next = td->bucket_next = bucket_offset + radix_size;
	bucket_prefix = td->bucket_prefix = bucket_next + num_buckets;

	/* initialize; every radix bin gets one bucket. Also remember
	   the index of the last bucket allocated. Note that we only
	   store the bottom 32+radix_bits bits of keys */

	last_bucket = radix_size - 1;
	memset(bucket_len, 0, radix_size * sizeof(uint32));
	for (i = 0; i < radix_size; i++)
		bucket_offset[i] = i << BUCKET_BITS;

	/* move each key */

	for (i = 0; i < num_sort; i++) {
		uint64 k = key[i];
		uint32 bin = (uint32)(k & radix_mask);

		bucket_prefix[bucket_offset[bin]++] = (uint32)(k >> radix_bits);

		if ((bucket_offset[bin] & BUCKET_MASK) == 0) {

			/* bucket is full; allocate the next
			   one and extend the chain of buckets
			   for this bin. Also undo the increment
			   so the offset in the next bucket 
			   starts at zero */

			uint32 curr_bucket = (bucket_offset[bin] - 1) >> 
							BUCKET_BITS;
			bucket_next[curr_bucket] = ++last_bucket;
			bucket_offset[bin] = last_bucket << BUCKET_BITS;
			bucket_len[bin] += BUCKET_SIZE;
		}
	}

	/* account for partially filled bins and return the
	   largest bin size */

	max_bin_size = 0;
	for (i = 0; i < radix_size; i++) {
		bucket_len[i] += bucket_offset[i] & BUCKET_MASK;
		max_bin_size = MAX(max_bin_size, bucket_len[i]);
	}
	bucket_next[last_bucket] = 0;
	return max_bin_size;
}

/*------------------------------------------------------------------------*/
static const double d_extra_bits = 4.5;

static SHOW_PROF uint32
prune_dup_arrays(cpu_thread_data_t *td, uint32 match_bits, 
		uint32 dup_array_bits, uint32 num_prefix,
		uint32 estimated_dups)
{
	const uint32 max_iter = 20;
	uint32 i, j;
	uint32 *T = td->dup0;
	uint32 *T2 = td->dup1;
	uint32 *T3 = td->dup2;
	uint32 *bin_prefix = td->bin_prefix;
	uint32 num_survivors[max_iter];
	uint32 dup_shift = 0;
	uint32 dup_mask = (1 << dup_array_bits) - 1;

	/* Using hash table find possible repeated keys, if f(x) 
	   is not repeated then x is also not repeated.
           For a hash function use f(x)=floor(x/2^k) mod (2^l) 
	   for different k,l values.
	   Using a table of size 2^l we can get the possible 
	   repeated keys, in another table store the possible 
	   repeated hash values. In another pass extract the 
	   possible repeated keys.
	   In each ping-pong round we can reduce the size of 
	   the array (possible number of repeated keys). */

	for(i = 0; i < max_iter; i++) {
        
		uint32 *U, *U2, *U3;
		uint32 dup_shift2;
		uint32 dup_mask2;
		uint32 dup_array_size2;
		uint32 curr_num_prefix = 0;
		uint32 ilog2 = MAX(5,(int)(log(estimated_dups + 1) / 
				M_LN2 + d_extra_bits));

		ilog2 = MIN(ilog2, match_bits);
        
		if (i % 2 == 0) {
			/* hash data is in T array */
			U = T; U2 = T2; U3 = T3;
		}
		else {
			/* hash data is in T3 array */
			U = T3; U2 = T2; U3 = T;
		}
        
		dup_shift2 = dup_shift + 6;
		if (dup_shift2 + ilog2 > match_bits)
			dup_shift2 = 0;
        
		dup_mask2 = (1 << ilog2) - 1;
		dup_array_size2 = 1 << (ilog2 - 5);
		memset(U2, 0, dup_array_size2 * sizeof(uint32));
		memset(U3, 0, dup_array_size2 * sizeof(uint32));
        
		for (j = 0; j < num_prefix; j++) {
			/* look for set bit in previous dup array,
			   fill current dup array */
			uint32 t = bin_prefix[j];
			uint32 bit = (t >> dup_shift) & dup_mask;
			if (U[bit >> 5] & (1 << (bit & 31))) {
				uint32 bit2 = (t >> dup_shift2) & dup_mask2;
				if (U2[bit2 >> 5] & (1 << (bit2 & 31))) 
					U3[bit2 >> 5] |= 1 << (bit2 & 31); 
				else
					U2[bit2 >> 5] |= 1 << (bit2 & 31);
				bin_prefix[curr_num_prefix++] = t;
			}
		}
        
		num_prefix = num_survivors[i] = curr_num_prefix;

		/* fast paths */
		if (num_prefix == 0)
			return 0;
		else if (num_prefix == 2)
			return (bin_prefix[0] == bin_prefix[1]) ? 2 : 0;
		else if (i >= 3 && num_survivors[i-3] == num_survivors[i])
			break;

		dup_mask = dup_mask2;
		dup_shift = dup_shift2;
	}
	return num_prefix;
}

/*------------------------------------------------------------------------*/
static void
find_dup_one_bin(cpu_thread_data_t *td, 
		uint32 match_bits, uint32 dup_array_bits,
	        uint32 radix_bits, uint32 which_bin)
{
	uint32 dup_array_mask = (1 << dup_array_bits) - 1;
	uint32 dup_array_size = 1 << (dup_array_bits - 5);
	uint32 bucket = which_bin; /* first bucket at known position */
	uint32 bucket_off = bucket << BUCKET_BITS;
	uint32 *bucket_prefix = td->bucket_prefix;
	uint32 *bucket_next = td->bucket_next;
	uint32 max_bucket_off = td->bucket_offset[which_bin];
	uint32 num_dups = 0;

	uint32 *bin_prefix = td->bin_prefix;
	uint32 i = 0;
	uint32 num_survivors;
	uint32 *T = td->dup0;
	uint32 *T2 = td->dup1;

	/* initialize duplicate arrays */

	memset(T, 0, dup_array_size * sizeof(uint32));
	memset(T2, 0, dup_array_size * sizeof(uint32));
    
	/* proceed one bucket at a time, concatenating
	   all the buckets into bin_prefix and filling
	   the duplicate array T */

	while (1) {
		uint32 end = MIN(max_bucket_off, 
				bucket_off + BUCKET_SIZE);

		/* use two separate bit arrays to remember 
		   key groups with matching bottom bits */

		while (bucket_off < end) {
			uint32 curr = bucket_prefix[bucket_off++];
			uint32 bit = curr & dup_array_mask;
			bin_prefix[i++] = curr;
			if (T2[bit >> 5] & (1 << (bit & 31))) {
				T[bit >> 5] |= 1 << (bit & 31);
				num_dups++;
			}
			else {
				T2[bit >> 5] |= 1 << (bit & 31);
			}
		}
		/* find the next bucket for which_bin */
		if (bucket_off == max_bucket_off)
			break;
		bucket = bucket_next[bucket];
		if (bucket == 0)
			break;
		bucket_off = bucket << BUCKET_BITS;
	}

	/* initial bit array is written, iteratively prune the
	   list of false positives */

	if (num_dups == 0)
		return;

	num_survivors = prune_dup_arrays(td, match_bits, 
				dup_array_bits, i, 2 * num_dups);

	/* save any survivors */

	if (num_survivors > 0) {
		uint32 c = td->num_survivors;
		if (c + num_survivors >= td->num_survivors_alloc) {
			uint32 t = 2 * (c + num_survivors);
			td->num_survivors_alloc = t;
			td->survivors = (uint64 *)xrealloc(td->survivors,
							t * sizeof(uint64));
			td->sort_survivors = (uint64 *)xrealloc(
							td->sort_survivors,
							t * sizeof(uint64));
		}
		for (i = 0; i < num_survivors; i++) {
			td->survivors[c + i] = (uint64)bin_prefix[i] << 
						radix_bits | which_bin;
		}
		td->num_survivors += num_survivors;
	}
}

/*------------------------------------------------------------------------*/
static SHOW_PROF void 
find_dup_all(cpu_thread_data_t *td, uint32 match_bits,
		uint32 dup_array_bits)
{
	uint32 i;
	uint32 num_sort = td->num_sort;
	uint64 *key = td->sort_key;
	uint32 dup_array_mask = (1 << dup_array_bits) - 1;
	uint32 dup_array_size = 1 << (dup_array_bits - 5);
	uint32 num_dups = 0;

	uint32 *bin_prefix = td->bin_prefix;
	uint32 num_survivors;
	uint32 *T = td->dup0;
	uint32 *T2 = td->dup1;

	/* initialize duplicate arrays */

	memset(T, 0, dup_array_size * sizeof(uint32));
	memset(T2, 0, dup_array_size * sizeof(uint32));
    
	/* save prefix of each key directly */

	for (i = 0; i < num_sort; i++) {
		uint32 curr = (uint32)key[i];
		uint32 bit = curr & dup_array_mask;
		bin_prefix[i] = curr;
		if (T2[bit >> 5] & (1 << (bit & 31))) {
			T[bit >> 5] |= 1 << (bit & 31);
			num_dups++;
		}
		else {
			T2[bit >> 5] |= 1 << (bit & 31);
		}
	}

	/* initial bit array is written, iteratively prune the
	   list of false positives */

	if (num_dups == 0)
		return;

	num_survivors = prune_dup_arrays(td, match_bits, 
				dup_array_bits, num_sort,
				2 * num_dups);

	/* save any survivors */

	if (num_survivors > td->num_survivors_alloc) {
		uint32 t = 2 * num_survivors;
		td->num_survivors_alloc = t;
		td->survivors = (uint64 *)xrealloc(td->survivors,
						t * sizeof(uint64));
		td->sort_survivors = (uint64 *)xrealloc(
						td->sort_survivors,
						t * sizeof(uint64));
	}
	for (i = 0; i < num_survivors; i++)
		td->survivors[i] = bin_prefix[i];
	td->num_survivors = num_survivors;
}

/*------------------------------------------------------------------------*/
static int 
compare64(const void* a, const void* b) 
{
     uint64 a64 = *((uint64*)a);
     uint64 b64 = *((uint64*)b);
     
     if (a64 == b64)
	     return 0;
     else if (a64 < b64) 
	     return -1;
     return 1;
}

static int 
compare_match(const void* a, const void* b)
{
     match_t *a0 = (match_t *)a;
     match_t *b0 = (match_t *)b;
     
     if (a0->key == b0->key)
	     return 0;
     else if (a0->key < b0->key) 
	     return -1;
     return 1;
}

static SHOW_PROF void
finish_search(task_data_t *task, uint32 threadid, uint32 match_bits)
{
	stage1_sieve_data_t *d = task->d;
	stage1_sieve_thread_data_t *t = d->threads + threadid;
	cpu_thread_data_t *td = (cpu_thread_data_t *)t->hw_thread_data;

	uint32 i, j;
	uint32 num_sort = td->num_sort;
	uint64 * sort_key = td->sort_key;
	uint32 * sort_data = td->sort_data;
	uint32 num_survivors = td->num_survivors;
	uint64 * survivors = td->survivors;
	uint64 * sort_survivors = td->sort_survivors;
	uint64 match_mask = ((uint64)1 << match_bits) - 1;
	uint32 num_match = 0;
	p_packed_var_t *p_array = &td->p_array;
	uint64 sieve_size = p_array->sieve_size;

	/* presence array state */

	uint32 p_bits;
	uint32 p_size;
	uint32 p_mask;
	uint32 p_alloc;
	uint32 * p_count;
	uint32 * p_bit_array;

	/* find repeated keys and store without multiplicity */

	qsort(survivors, num_survivors, sizeof(uint64), compare64);
	for (i = 1, j = 0; i < num_survivors; i++) {
		if (survivors[i] == survivors[i-1] &&
		    (i == num_survivors - 1 || 
				survivors[i] != survivors[i+1]))
			survivors[j++] = survivors[i];
	}
	if (j == 0)
		return;
	num_survivors = j;

	/* list out entries for keys whose bottom match_bits
	   bits match one of the entries in survivors[] 
	   
	   With a small number of survivors we can just check
	   every key against every survivor. However, when a
	   large array is radix sorted we can get hundreds of
	   survivors in all bins combined, even if each bin
	   would only have one survivor. For these cases we 
	   need to radix sort the list of survivors to avoid
	   iterating through all of them */

	p_bits = MAX(5, (int)(log(num_survivors + 1) / 
				M_LN2 + d_extra_bits));
	p_size = 1 << (p_bits - 5);
	p_mask = (1 << p_bits) - 1;
	p_alloc = p_size + (1 << p_bits) + 1;
	if (p_alloc > td->num_dup_alloc) {
		uint32 t = 2 * p_alloc;
		td->num_dup_alloc = t;
		td->dup_data = (uint32 *)xrealloc(td->dup_data, 
					t * sizeof(uint32));
	}
	p_bit_array = td->dup_data;
	p_count = td->dup_data + p_size;
	memset(p_bit_array, 0, p_alloc * sizeof(uint32));

	/* increment the count in the *next* bucket for survivor i */

	for (i = 0; i < num_survivors; i++) {
		uint32 bit = (uint32)survivors[i] & p_mask;
		p_count[bit + 1]++;
		p_bit_array[bit >> 5] |= 1 << (bit & 31);
	}

	/* prefix sum */

	for (i = 1; i < ((uint32)1 << p_bits) + 1; i++)
		p_count[i] += p_count[i-1];

	/* scatter the survivor list */

	for (i = 0; i < num_survivors; i++) {
		uint32 bit = (uint32)survivors[i] & p_mask;
		sort_survivors[p_count[bit]++] = survivors[i];
	}

	/* iterate through keys */

	for (i = 0; i < num_sort; i++) {
		uint64 key = sort_key[i];
		uint32 bit = (uint32)key & p_mask;
		uint32 s_start, s_end;

		if (!(p_bit_array[bit >> 5] & (1 << (bit & 31))))
			continue;

		s_start = (bit == 0) ? 0 : p_count[bit - 1];
		s_end = p_count[bit];

		for (j = s_start; j < s_end; j++) {
			if ((key & match_mask) == sort_survivors[j]) {
				if (num_match >= td->num_match_alloc) {
					uint32 t = 2 * num_match;
					td->num_match_alloc = t;
					td->match = (match_t *)xrealloc(
							td->match, 
							t * sizeof(match_t));
				}
				td->match[num_match].key = key;
				td->match[num_match].data = sort_data[i];
				num_match++;
				break;
			}
		}
	}

	/* perform final collision search using entirety of keys */

	qsort(td->match, num_match, sizeof(match_t), compare_match);
	for (i = 0; i < num_match - 1; i++) {
		for (j = i + 1; j < num_match && 
				td->match[i].key == td->match[j].key; j++) {

			uint32 p0 = td->match[i].data;
			uint32 p1 = td->match[j].data;
			if (mp_gcd_1(p0, p1) != 1)
				continue;

			handle_collision(task, (uint64)p0 * p1,
				   td->curr_q, td->curr_q_root,
				   (int64)td->match[i].key - sieve_size / 2);
		}
	}
}

/*------------------------------------------------------------------------*/
static SHOW_PROF void
collision_search(task_data_t *task, uint32 threadid, uint32 key_bits)
{
	/* return the number of pairs of keys that are
	   identical in the input data. We assume the array
	   size is small to medium (less than maybe 1M 
	   keys for the largest problems) but key_bits is 
	   potentially very large (at least 32 bits, up 
	   to maybe 56-64 bits) */

	stage1_sieve_data_t *d = task->d;
	stage1_sieve_thread_data_t *t = d->threads + threadid;
	cpu_thread_data_t *td = (cpu_thread_data_t *)t->hw_thread_data;

	uint32 i; 
	uint32 max_bin_size;
	uint32 dup_array_bits;
	uint32 dup_array_size;
	uint32 match_bits;
	uint32 num_sort = td->num_sort;
	uint32 radix_bits = 0;
	uint32 dup_words;

	/* small problems use the input keys directly;
	   otherwise perform a single bucket sort pass 
	   to get future tables to approximately L1 
	   cache size. In either case, for large keys
	   we limit the size of internal data structures
	   by matching only the bottom (32+radix_bits)
	   bits of keys, and rely on the number of keys
	   to be small enough that false positive matches
	   are unlikely */

	if (num_sort < 55000) {
		max_bin_size = num_sort;
	}
	else {
		const uint32 target_dupsize = 5000;
		radix_bits = (uint32)(log((double)num_sort / 
					target_dupsize) / M_LN2 + 0.5);
		radix_bits = MIN(radix_bits, MAX_RADIX_BITS);
		max_bin_size = bucket_sort(td, radix_bits);
	}

	/* resize duplicate detection arrays if necessary,
	   and concatenate the arrays to avoid associative
	   conflicts in L1. 
	   
	   The duplicate detection arrays have 1-bit type 
	   so the actual memory allocated is much smaller */

	dup_array_bits = MAX(5,(int)((double)log(max_bin_size + 1) / 
				M_LN2 + d_extra_bits));
	match_bits = MIN(32, key_bits - radix_bits);
	dup_array_bits = MIN(dup_array_bits, match_bits);
	dup_array_size = 1 << (dup_array_bits - 5);
	dup_words = 3 * dup_array_size + max_bin_size;
	if (dup_words > td->num_dup_alloc) {
		uint32 t = dup_words * 3 / 2;
		td->num_dup_alloc = t;
		aligned_free(td->dup_data);
		td->dup_data = (uint32 *)aligned_malloc(
					t * sizeof(uint32), 64);
	}
	td->dup0 = td->dup_data;
	td->dup1 = td->dup_data + dup_array_size;
	td->dup2 = td->dup_data + 2 * dup_array_size;
	td->bin_prefix = td->dup_data + 3 * dup_array_size;
    
	/* for small problems, detect duplicates directly; 
	   otherwise radix sorted data proceeds one bin at a time */

	if (radix_bits) {
		td->num_survivors = 0;
		for (i = 0; i < ((uint32)1 << radix_bits); i++) {
			if (td->bucket_len[i] == 0)
				continue;

			find_dup_one_bin(td, match_bits,
					dup_array_bits, radix_bits, i);
		}
	}
	else {
		find_dup_all(td, match_bits, dup_array_bits);
	}
    
	if (td->num_survivors > 0)
		finish_search(task, threadid, radix_bits + match_bits);
}

/*------------------------------------------------------------------------*/
static SHOW_PROF void
special_q_trivial(task_data_t *task, uint32 threadid, uint32 key_bits)
{
	uint32 i, j;
	stage1_sieve_data_t *d = task->d;
	stage1_sieve_thread_data_t *t = d->threads + threadid;
	cpu_thread_data_t *td = (cpu_thread_data_t *)t->hw_thread_data;
	p_packed_var_t *p_array = &td->p_array;
	uint32 num_entries = p_array->num_p;
	uint64 sieve_size = p_array->sieve_size;
	p_packed_t *tmp = p_array->packed_array;
	uint32 num_sort = 0;
	uint32 num_sort_alloc = td->num_sort_alloc;

	for (i = 0; i < num_entries; i++) {
		uint32 p = tmp->p;
		uint64 pp = tmp->pp;
		uint32 num_roots = tmp->num_roots;

		for (j = 0; j < num_roots; j++) {
			uint64 proot = tmp->roots[j];

			proot = mp_modsub_2(proot, tmp->ss_mod_pp, pp);

			if (num_sort + 100 >= num_sort_alloc) {
				num_sort_alloc = 3 * num_sort_alloc / 2;
				grow_sort(td, num_sort_alloc);
			}
			while (proot < sieve_size) {
				td->sort_key[num_sort] = proot;
				td->sort_data[num_sort] = p;
				num_sort++;
				proot += pp;
			}
		}

		tmp = p_packed_next(tmp);
	}
	td->num_sort = num_sort;
	collision_search(task, threadid, key_bits);
}

/*------------------------------------------------------------------------*/
static SHOW_PROF void
special_q_nontrivial(task_data_t *task, uint32 threadid,
		     uint64 *inv_array, uint32 key_bits)
{
	/* queue up the roots for a single special-q

	   inv_array stores the modular inverse of q^2 mod p^2 for
	   each progression p in hash_array, in Montgomery
	   form */

	uint32 i, j;
	stage1_sieve_data_t *d = task->d;
	stage1_sieve_thread_data_t *t = d->threads + threadid;
	cpu_thread_data_t *td = (cpu_thread_data_t *)t->hw_thread_data;
	p_packed_var_t *p_array = &td->p_array;
	uint32 num_entries = p_array->num_p;
	uint64 sieve_size = p_array->sieve_size;
	p_packed_t *tmp = p_array->packed_array;
	uint32 num_sort = 0;
	uint32 num_sort_alloc = td->num_sort_alloc;

	/* for each progression p */

	for (i = 0; i < num_entries; i++) {
		/* skip if p and q have a factor in common */

		if (inv_array[i]) {
			uint32 p = tmp->p;
			uint64 pp = tmp->pp;
			uint32 num_roots = tmp->num_roots;
			uint64 pp_t = tmp->mont_t;
			uint64 pp_w = tmp->mont_w;
			uint64 qinv = inv_array[i];
			uint64 sq_root;

			/* for each root R mod p^2, we need the first
			   value of R + k * p^2 that also falls on
			   special_q_root + m * special_q^2. This is 
			   a standard arithmetic problem, finding the
			   intersection of two arithmetic progressions */

			sq_root = mod128_64(td->curr_q_root, pp_t, pp, pp_w);
			for (j = 0; j < num_roots; j++) {
				uint64 proot = tmp->roots[j];

				proot = mp_modsub_2(proot, sq_root, pp);
				proot = montmul64(proot, qinv, pp, pp_w);
				proot = mp_modsub_2(proot, tmp->ss_mod_pp, pp);

				if (num_sort + 100 >= num_sort_alloc) {
					num_sort_alloc = 3 * num_sort_alloc / 2;
					grow_sort(td, num_sort_alloc);
				}
				while (proot < sieve_size) {
					td->sort_key[num_sort] = proot;
					td->sort_data[num_sort] = p;
					num_sort++;
					proot += pp;
				}
			}
		}
		tmp = p_packed_next(tmp);
	}
	td->num_sort = num_sort;
	collision_search(task, threadid, key_bits);
}

/*------------------------------------------------------------------------*/
static SHOW_PROF void
batch_invert(p_packed_var_t * p_array, q_packed_var_t * q_array,
		uint64 * invtable)
{
	/* Assume invtable is a num_q x num_p array; for each p,
	   use Montgomery's batch modular inversion algorithm
	   to compute the inverse of q^2 mod p^2 for each q.
	   For one more Montgomery multiply, leave all the outputs 
	   in Montgomery form even though they didn't start that way */

	uint32 i, j;
	uint32 num_p = p_array->num_p;
	uint32 num_q = q_array->num_q;
	uint64 qsave[SPECIALQ_BATCH_SIZE];
	uint64 invlist[SPECIALQ_BATCH_SIZE];
	p_packed_t * p_ptr = p_array->packed_array;

	for (i = 0; i < num_p; i++, invtable++, p_ptr = p_packed_next(p_ptr)) {

		uint32 p = p_ptr->p;
		uint64 pp = p_ptr->pp;
		uint64 pp_r = p_ptr->mont_r;
		uint64 pp_w = p_ptr->mont_w;
		uint64 invprod = 0;
		q_packed_t * q_ptr = q_array->packed_array;
		uint64 * inv_ptr = invtable + (num_q - 1) * num_p;
		uint32 first = 0;

		for (j = 0; j < num_q; j++, q_ptr = q_packed_next(q_ptr)) {

			uint64 q = q_ptr->q;

			/* current p and current q must be coprime */

			uint32 qmodp = (q >> 32) ? mp_mod64(q, p) :
						(uint32)q % p;
			if (mp_gcd_1(qmodp, p) != 1) {
				qsave[j] = 0;
			}
			else {
				/* save q and accumulate without converting
				   to Montgomery form. This means invprod
				   is scaled by an increasing negative power 
				   of 2 */

				qsave[j] = q;
				if (invprod == 0) {
					first = j;
					invprod = q;
				}
				else {
					invprod = montmul64(invprod, q, pp, pp_w);
				}
			}
			invlist[j] = invprod;
		}

		j = num_q - 1;
		if (invprod) {
			/* invert the accumulated product and convert
			   (once) to Montgomery form. This turns the
			   negative power of 2 scale factor positive */

			invprod = mp_modinv_2(invprod, pp);
			invprod = montmul64(invprod, pp_r, pp, pp_w);

			/* peel off results in reverse order */

			for (; j > first; j--, inv_ptr -= num_p) {
				uint64 inv_q;

				if (qsave[j] == 0) {
					*inv_ptr = 0;
					continue;
				}

				/* multiplying by the previous invprod
				   will cancel out all the previous q's 
				   and almost all the powers of 2, leaving 
				   the inverse of the j_th q in Montgomery 
				   form. Square that to compute the j_th 
				   result, then remove the inverse of the 
				   current q from invprod */

				inv_q = montmul64(invprod, invlist[j-1], pp, pp_w);
				*inv_ptr = montmul64(inv_q, inv_q, pp, pp_w);
				invprod = montmul64(invprod, qsave[j], pp, pp_w);
			}

			/* initial nonzero result */
			*inv_ptr = montmul64(invprod, invprod, pp, pp_w);
			j--; inv_ptr -= num_p;
		}
		for (; (int32)j >= 0; j--, inv_ptr -= num_p)
			*inv_ptr = 0;
	}
}

/*------------------------------------------------------------------------*/
void
stage1_specialq_cpu(task_data_t *task, uint32 threadid,
		uint64 special_q_min, uint64 special_q_max, 
		uint32 p_min, uint32 p_max)
{
	/* core search code */

	uint32 i, j;
	uint32 num_q_roots = 0;
	uint32 num_p;
	uint32 num_p_roots;
	uint32 sort_key_bits;
	uint64 sieve_size;
	double cpu_start_time = get_cpu_time();

	msieve_obj * obj = task->obj;
	stage1_sieve_data_t *d = task->d;
	stage1_sieve_thread_data_t *t = d->threads + threadid;
	cpu_thread_data_t *td = (cpu_thread_data_t *)t->hw_thread_data;
	poly_coeff_t *c = task->c;
	p_packed_var_t *p_array = &td->p_array;
	q_packed_var_t *q_array = &td->q_array;

	/* Each of the p values needs at least one entry 
	   in the list to search, so the 'sieve size' implicit 
	   in the dataset is the square of the largest p, 
	   and the sieve starts at -sieve_size / 2 */

	sieve_size = (uint64)p_max * p_max & (uint64)(~1);
	td->p_array.sieve_size = sieve_size;

	/* find the size of keys for collision search */

	sort_key_bits = 64;
	while (!(sieve_size & ((uint64)1 << (sort_key_bits - 1))))
		sort_key_bits--;

	/* build all the p arithmetic progressions */

	p_packed_reset(p_array);
	sieve_fb_reset(t->sieve_p_fb, p_min, p_max, 1, MAX_ROOTS);
	while (sieve_fb_next(t->sieve_p_fb, c, store_p_packed,
			p_array) != P_SEARCH_DONE) {
		;
	}

	num_p = p_array->num_p;
	num_p_roots = p_array->num_roots;
#if 0
	printf("pprogs: %u entries, %u x %u-bit roots for sieve_size %lu\n", 
			num_p, num_p_roots, sort_key_bits, sieve_size);
#endif

	/* set up cache for the inverse of each special-q 
	   modulo each p in p_array */
	   
	if (td->num_inv_alloc < num_p * SPECIALQ_BATCH_SIZE) {
		td->num_inv_alloc = num_p * SPECIALQ_BATCH_SIZE * 3 / 2;
		td->invtable = (uint64 *)xrealloc(td->invtable,
					td->num_inv_alloc * sizeof(uint64));
	}

	/* guess how many roots will participate in a single
	   instance of the collision search. The smallest p
	   will contribute up to 4 roots and the largest will
	   contribute only one */

	i = 3 * num_p_roots;
	td->num_sort = 0;
	if (td->num_sort_alloc < i)
		grow_sort(td, i * 5 / 4);

	/* handle trivial lattice */

	if (special_q_min == 1) {
		td->curr_q = 1;
		memset(&td->curr_q_root, 0, sizeof(uint128));

		special_q_trivial(task, threadid, sort_key_bits);
		num_q_roots++;
		if (special_q_max == 1)
			goto finished;
	}

	sieve_fb_reset(t->sieve_q_fb, special_q_min, 
			special_q_max, 1, MAX_ROOTS);
	while (1) {
		q_packed_t * qptr;
		uint64 * inv_array;

		/* allocate the next batch of special-q and all the
		   roots they use */

		q_packed_reset(q_array);
		while (sieve_fb_next(t->sieve_q_fb, c, store_q_packed,
					q_array) != P_SEARCH_DONE) {
			if (q_array->num_q == SPECIALQ_BATCH_SIZE)
				break;
		}
#if 0
		printf("qprogs: %u entries, %u roots\n", 
				q_array->num_q, q_array->num_roots);
#endif
		if (q_array->num_q == 0)
			break;

		/* invert all the special-q at once modulo each p in 
		   hash_array.
		 
		   Note that we need one inverse for each (p,special_q)
		   pair, *not* one for each special-q root. For composite 
		   special-q that have many roots, this speeds up 
		   the modular inverse phase by much more than 
		   SPECIALQ_BATCH_SIZE */

		batch_invert(p_array, q_array, td->invtable);

		/* for each root of each special-q */

		qptr = q_array->packed_array;
		inv_array = td->invtable;
		for (i = 0; i < q_array->num_q; 
			i++, inv_array += num_p, qptr = q_packed_next(qptr)) {

			td->curr_q = qptr->q;
			for (j = 0; j < qptr->num_roots; 
						j++, num_q_roots++) {

				td->curr_q_root = qptr->roots[j];

				special_q_nontrivial(task, threadid, 
						inv_array, sort_key_bits);

				/* check for interrupt, check less often 
				   for timeout */

				if ((obj->flags & MSIEVE_FLAG_STOP_SIEVING) ||
				    (task->coeff_deadline && 
				     (num_q_roots+1) % 1000 == 0 && 
				      (get_cpu_time() - cpu_start_time) > 
						task->coeff_deadline)) {
					goto finished;
				}
			}
		}
	}

finished:
	//logprintf(obj, "searched %u special-q\n", num_q_roots);
	t->cumulative_elapsed += get_cpu_time() - cpu_start_time;
}
