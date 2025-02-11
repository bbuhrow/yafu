/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: stage1_sieve_cpu.c 1025 2018-08-19 02:20:28Z jasonp_sf $
--------------------------------------------------------------------*/

#include <stage1.h>
#include <cpu_intrinsics.h>

/* CPU collision search; this code looks for self-collisions
   among arithmetic progressions, by finding k1 and k2 such that
   for two arithmetic progressions r1+k*p1^2 and r2+k*p2^2 we
   have

      r1 + k1*p1^2 = r2 + k2*p2^2

   such that
      - p1 and p2 are coprime and < 2^32
      - the value where they coincide is of size smaller
        than a fixed bound

   This code uses a hashtable to find collisions across all the
   p1 and p2 in the set simultaneously. The size of the range is
   comparatively very large, and the algorithm as described above
   is only efficient for pretty small problems. We cannot
   practically fill up a sigle hashtable with all the possibilities
   because for e.g. 512-bit GNFS the range on k1 and k2 is typically 
   around 10^6 and the set has ~10^6 progressions. We can reduce the
   memory use by breaking the hashtable into blocks of size > p_min
   and filling each block individually, but that only reduces the 
   memory use and not the actual work required.
   
   To scale up the idea, we further use a 'special-q' formulation 
   where all the inputs to the hashtable are constrained to fall 
   on a third arithmetic progression r3 + k*q^2 for some k. We choose 
   a given q and for each of its roots run the complete hashtable search.
   This is analogous to lattice sieving across the interval.
   
   This allows us to choose q so that the hashtable problem is of 
   reasonable size but the collisions found are still over the 
   original, impractically large range. 

   A side effect of this is that the complete range must be < 2^96
   in size, though this appears to be sufficient for very large
   problems, e.g. 200-digit GNFS */

#define NEW_P_PACKED
#define P_BITS 27
#define HASHITER_BITS (32 - P_BITS)

/* one element of the hashtable. The table is keyed by the
   offset ito the sieve interval, i.e. by the first field below. */

typedef struct {
	int64 offset;
	uint32 p : P_BITS;
	uint32 iter : HASHITER_BITS;
} hash_entry_t;

/* every progression keeps a running count of the current sieve
   offset for each of its roots */

typedef struct {
	uint64 start_offset;
	int64 offset;
} hash_list_t;

/* structure for a single progression */

typedef struct {
	uint32 p;
	uint32 num_roots;
	uint32 pad;
	uint32 mont_w;
	uint64 mont_r;
	uint64 pp;
	uint64 ss_mod_pp;
#ifdef NEW_P_PACKED
	// refer to the start index of our roots within
	// a globally managed list of roots.
	uint32_t roots_idx;
#else
	// each progression struct keeps its own roots
	// which are packed into a larger structure of many progressions.
	hash_list_t roots[MAX_ROOTS];
#endif
	
} p_packed_t;

/* P_PACKED_HEADER_WORDS is the number of uint64 sized elements
   preceding hash_list_t roots in the above structure */

#define P_PACKED_HEADER_WORDS 5

/*------------------------------------------------------------------------*/

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
#ifndef NEW_P_PACKED
	p_packed_t *curr;
#endif
	p_packed_t *packed_array;

#ifdef NEW_P_PACKED
	uint32_t alloc_roots;
	hash_list_t* roots;
#endif
} p_packed_var_t;

/* list manipulation functions */
//static
void 
p_packed_init(p_packed_var_t *s)
{
	memset(s, 0, sizeof(p_packed_var_t));

	s->p_size_alloc = 100;
#ifdef NEW_P_PACKED
	s->packed_array = (p_packed_t *)xmalloc(s->p_size_alloc *
						sizeof(p_packed_t));
	// one globally managed list of roots
	s->alloc_roots = 2048;
	s->roots = (hash_list_t*)xmalloc(s->alloc_roots * sizeof(hash_list_t));
#else
	s->packed_array = s->curr = (p_packed_t *)xmalloc(s->p_size_alloc *
						sizeof(p_packed_t));
#endif
}

//static
void 
p_packed_free(p_packed_var_t *s)
{
#ifdef NEW_P_PACKED
	free(s->roots);
#endif
	free(s->packed_array);
}

//static
void
p_packed_reset(p_packed_var_t *s)
{
	s->num_p = s->num_roots = 0;
#ifndef NEW_P_PACKED
	s->curr = s->packed_array;
#endif
}

//static
p_packed_t * 
p_packed_next(p_packed_t *curr)
{
	return (p_packed_t *)((uint64 *)curr + 
			P_PACKED_HEADER_WORDS + 2 * curr->num_roots);
}

// gprof said that store_p_packed accounted for 80+% of the runtime of
// the entire polyselect, which is outrageous.  While I don't entirely
// believe that figure, I think this is a better and simpler way to
// do the packing.  Retests showed that this does improve things, but
// of course not by 80%.
#ifdef NEW_P_PACKED
void
store_p_packed(uint32 p, uint32 num_roots, uint64* roots, void* extra)
{
	/* used to insert a new arithmetic progression and all its
	   roots into the list */

	uint32 i;
	p_packed_var_t* s = (p_packed_var_t*)extra;
	p_packed_t* curr;

	if ((s->num_p + 1) > s->p_size_alloc)
	{
		// double the size of the packed array
		s->p_size_alloc *= 2;
		s->packed_array = (p_packed_t*)xrealloc(
			s->packed_array,
			s->p_size_alloc *
			sizeof(p_packed_t));
	}

	if ((s->num_roots + num_roots) >= s->alloc_roots)
	{
		// double the size of the roots array
		s->alloc_roots *= 2;
		s->roots = (hash_list_t*)xrealloc(
			s->roots,
			s->alloc_roots *
			sizeof(hash_list_t));
	}

	curr = &s->packed_array[s->num_p];
	curr->p = p;
	curr->pad = 0;
	curr->num_roots = num_roots;
	curr->roots_idx = s->num_roots;		// the start of our roots within the global list
	for (i = 0; i < num_roots; i++)
		s->roots[s->num_roots + i].start_offset = roots[i];
	curr->pp = (uint64)p * p;

	/* set up Montgomery arithmetic */

	curr->mont_w = montmul32_w((uint32)curr->pp);
	curr->mont_r = montmul64_r(curr->pp);

	curr->ss_mod_pp = s->sieve_size % curr->pp;

	s->num_p++;
	s->num_roots += num_roots;
}
#else
//static
void 
store_p_packed(uint32 p, uint32 num_roots, uint64 *roots, void *extra)
{
	/* used to insert a new arithmetic progression and all its
	   roots into the list */

	uint32 i;
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

	curr = s->curr;
	curr->p = p;
	curr->pad = 0;
	curr->num_roots = num_roots;
	for (i = 0; i < num_roots; i++)
		curr->roots[i].start_offset = roots[i];
	curr->pp = (uint64)p * p;

	/* set up Montgomery arithmetic */

	curr->mont_w = montmul32_w((uint32)curr->pp);
	curr->mont_r = montmul64_r(curr->pp);

	curr->ss_mod_pp = s->sieve_size % curr->pp;

	s->num_p++;
	s->num_roots += num_roots;
	s->curr = p_packed_next(s->curr);
}
#endif


#define SIEVE_MAX 9223372036854775807ULL

#define MAX_SPECIAL_Q ((uint32)(-1))
#define MAX_OTHER (((uint32)1 << P_BITS) - 1)

/*------------------------------------------------------------------------*/
//static
uint32
handle_special_q(msieve_obj *obj, poly_search_t *poly, poly_coeff_t *c,
		hash_entry_t *hashtable, uint32 hashtable_size_log2,
		p_packed_var_t *hash_array, uint32 special_q,
		uint64 special_q_root, uint64 block_size, uint64 *inv_array, uint32 *num_sizeopt)
{
	/* perform the hashtable search for a single special-q

	   inv_array stores the modular inverse of q^2 mod p^2 for
	   each progression p in hash_array, in Montgomery
	   form */

	uint32 i, j;
	uint32 quit = 0;
	p_packed_t *tmp;
	uint32 num_entries = hash_array->num_p;
	uint64 sieve_size = hash_array->sieve_size;
	int64 sieve_start = -sieve_size;
	uint32 num_blocks = 0;
	uint32 hashtable_count;
	uint32 hashtable_size = (uint32)1 << hashtable_size_log2;
	uint32 hashmask = hashtable_size - 1;

	if (sieve_start > 0) {
		printf("error: sieve size must fit in signed int "
			"in handle_special_q()\n");
		exit(1);
	}

	tmp = hash_array->packed_array;

	if (special_q == 1) {

#ifdef NEW_P_PACKED
		for (i = 0; i < num_entries; i++) {
			uint64 pp = tmp[i].pp;
			uint32 num_roots = tmp[i].num_roots;
			uint32 roffset = tmp[i].roots_idx;

			for (j = 0; j < num_roots; j++) {
				uint64 proot = hash_array->roots[roffset + j].start_offset;

				proot = mp_modsub_2(proot,
					pp - tmp[i].ss_mod_pp, pp);
				hash_array->roots[roffset + j].offset = proot - sieve_size;
			}
		}
#else
		for (i = 0; i < num_entries; i++) {
			uint64 pp = tmp->pp;
			uint32 num_roots = tmp->num_roots;

			for (j = 0; j < num_roots; j++) {
				uint64 proot = tmp->roots[j].start_offset;

				proot = mp_modsub_2(proot,
						pp - tmp->ss_mod_pp, pp);
				tmp->roots[j].offset = proot - sieve_size;
			}

			tmp = p_packed_next(tmp);
		}
#endif
	}
	else {
		/* for each progression p */

#ifdef NEW_P_PACKED
		for (i = 0; i < num_entries; i++) {
			uint64 pp = tmp[i].pp;
			uint32 num_roots = tmp[i].num_roots;
			uint32 pp_w = tmp[i].mont_w;
			uint64 qinv = inv_array[i];
			uint32 roffset = tmp[i].roots_idx;

			/* check if calling code determined that p and
			   q have a factor in common. Skip p if so;
			   actually we do this by pushing the starting
			   offset all the way to the end of the sieve
			   interval */

			if (qinv == 0) {
				for (j = 0; j < num_roots; j++)
					hash_array->roots[roffset + j].offset = SIEVE_MAX;
				continue;
			}

			/* for each root R mod p^2, we need the first
			   value of R + k * p^2 that also falls on
			   special_q_root + m * special_q^2. This is
			   a standard arithmetic problem, finding the
			   intersection of two arithmetic progressions */

			for (j = 0; j < num_roots; j++) {
				uint64 proot = hash_array->roots[roffset + j].start_offset;

				proot = mp_modsub_2(proot,
					special_q_root % pp, pp);
				proot = montmul64(proot, qinv, pp, pp_w);
				proot = mp_modsub_2(proot,
					pp - tmp[i].ss_mod_pp, pp);

				hash_array->roots[roffset + j].offset = proot - sieve_size;
			}
		}
#else
		for (i = 0; i < num_entries; i++) {
			uint64 pp = tmp->pp;
			uint32 num_roots = tmp->num_roots;
			uint32 pp_w = tmp->mont_w;
			uint64 qinv = inv_array[i];

			/* check if calling code determined that p and
			   q have a factor in common. Skip p if so;
			   actually we do this by pushing the starting
			   offset all the way to the end of the sieve
			   interval */

			if (qinv == 0) {
				for (j = 0; j < num_roots; j++)
					tmp->roots[j].offset = SIEVE_MAX;
				tmp = p_packed_next(tmp);
				continue;
			}

			/* for each root R mod p^2, we need the first
			   value of R + k * p^2 that also falls on
			   special_q_root + m * special_q^2. This is 
			   a standard arithmetic problem, finding the
			   intersection of two arithmetic progressions */

			for (j = 0; j < num_roots; j++) {
				uint64 proot = tmp->roots[j].start_offset;

				proot = mp_modsub_2(proot,
						special_q_root % pp, pp);
				proot = montmul64(proot, qinv, pp, pp_w);
				proot = mp_modsub_2(proot,
						pp - tmp->ss_mod_pp, pp);
				tmp->roots[j].offset = proot - sieve_size;
			}

			tmp = p_packed_next(tmp);
		}
#endif
	}


	/* sieve is ready to go; we proceed with the hashtable
	   one block at a time. The block size is chosen to be
	   the minimum value of p^2 in hash_array, so that for each
	   block a given root for a given p contributes at most
	   one entry to the hashtable. The hashtable is refilled
	   from scratch when the next block runs */

#ifdef NEW_P_PACKED
	while (sieve_start < (int64)sieve_size) {
		uint32 iter = hash_array->curr_hashiter;
		int64 sieve_end = sieve_start + MIN(block_size,
			sieve_size - sieve_start);

		tmp = hash_array->packed_array;
		hashtable_count = 0;

		/* only reset the hashtable memory when the iteration
			count loops back to zero. In other cases, we can
			ignore stale data in the hashtable by checking its
			iteration count, and save a memset here */

		if (iter == 0) {
			memset(hashtable, 0, sizeof(hash_entry_t) *
				hashtable_size);
		}

		for (i = 0; i < num_entries; i++) {

			uint32 num_roots = tmp[i].num_roots;
			uint32 roffset = tmp[i].roots_idx;

			for (j = 0; j < num_roots; j++) {
				int64 offset = hash_array->roots[roffset + j].offset;
				uint32 key;

				if (offset >= sieve_end)
					continue;

				key = offset & hashmask;

				while (hashtable[key].iter == iter &&
					hashtable[key].p != 0 &&
					hashtable[key].offset != offset) {

					key = (key + 1) & hashmask;
				}

				if (hashtable[key].iter == iter &&
					hashtable[key].p != 0) {

					if (mp_gcd_1(hashtable[key].p, tmp[i].p) == 1) {

						/* collision found! The sieve
							offset is in 'special q	coordinates' */

						uint64 p;

						p = tmp[i].p;
						p = p * hashtable[key].p;

						if (handle_collision(c, p, special_q,
							special_q_root, offset) != 0) {

							(*num_sizeopt)++;
							poly->callback(c->high_coeff,
								c->p, c->m,
								poly->callback_data);
						}
					}
				}
				else {
					/* no hit; insert this offset
						for root j into the table */

					hashtable[key].iter = iter;
					hashtable[key].p = tmp[i].p;
					hashtable[key].offset = offset;

					if (++hashtable_count == hashtable_size)
						goto next_hashtable;
				}

				/* compute the next offset
					for root j */

				offset += tmp[i].pp;

				if (offset < sieve_end)
					/* handle overflow */
					offset = SIEVE_MAX;

				hash_array->roots[roffset + j].offset = offset;
			}
		}

	next_hashtable:
		sieve_start = sieve_end;
		num_blocks++;

		hash_array->curr_hashiter++;
		hash_array->curr_hashiter %= ((uint32)1 << HASHITER_BITS);

		/* check for interrupt after every block */

		if (obj->flags & MSIEVE_FLAG_STOP_SIEVING) {
			quit = 1;
			break;
		}
	}
#else
	while (sieve_start < (int64)sieve_size) {
		uint32 iter = hash_array->curr_hashiter;
		int64 sieve_end = sieve_start + MIN(block_size,
						sieve_size - sieve_start);

		tmp = hash_array->packed_array;
		hashtable_count = 0;

		/* only reset the hashtable memory when the iteration
		   count loops back to zero. In other cases, we can
		   ignore stale data in the hashtable by checking its
		   iteration count, and save a memset here */

		if (iter == 0) {
			memset(hashtable, 0, sizeof(hash_entry_t) *
							hashtable_size);
		}

		for (i = 0; i < num_entries; i++) {

			uint32 num_roots = tmp->num_roots;

			for (j = 0; j < num_roots; j++) {
				int64 offset = tmp->roots[j].offset;
				uint32 key;

				if (offset >= sieve_end)
					continue;

				key = offset & hashmask;

				while (hashtable[key].iter == iter &&
				      hashtable[key].p != 0 &&
				      hashtable[key].offset != offset) {

					key = (key + 1) & hashmask;
				}

				if (hashtable[key].iter == iter &&
				   hashtable[key].p != 0) {

					if (mp_gcd_1(hashtable[key].p,
							tmp->p) == 1) {

						/* collision found! The sieve
						   offset is in 'special q
						   coordinates' */

						uint64 p;

						p = tmp->p;
						p = p * hashtable[key].p;

						if (handle_collision(c, p, special_q,
							special_q_root, offset) != 0) {

							(*num_sizeopt)++;
							poly->callback(c->high_coeff,
								c->p, c->m,
								poly->callback_data);
						}
					}
				}
				else {
					/* no hit; insert this offset
					   for root j into the table */

					hashtable[key].iter = iter;
					hashtable[key].p = tmp->p;
					hashtable[key].offset = offset;

					if (++hashtable_count == hashtable_size)
						goto next_hashtable;
				}

				/* compute the next offset 
				   for root j */

				offset += tmp->pp;

				if (offset < sieve_end)
					/* handle overflow */
					offset = SIEVE_MAX;

				tmp->roots[j].offset = offset;
			}

			tmp = p_packed_next(tmp);
		}

next_hashtable:
		sieve_start = sieve_end;
		num_blocks++;

		hash_array->curr_hashiter++;
		hash_array->curr_hashiter %= ((uint32)1 << HASHITER_BITS);

		/* check for interrupt after every block */

		if (obj->flags & MSIEVE_FLAG_STOP_SIEVING) {
			quit = 1;
			break;
		}
	}
#endif

//	printf("%u\n", num_blocks); 
	return quit;
}

/*------------------------------------------------------------------------*/
#define SPECIALQ_BATCH_SIZE 10

//static
void
batch_invert(uint32 *qlist, uint32 num_q, uint64 *invlist,
		uint32 p, uint64 pp_r, uint32 pp_w)
{
	/* use Montgomery's batch modular inversion algorithm
	   to compute the inverse of many special q modulo a
	   single p. For one extra Montgomery multiplication,
	   leave all the outputs in Montgomery form even though
	   they didn't start that way */

	uint32 i;
	uint64 qq[SPECIALQ_BATCH_SIZE];
	uint64 invprod;
	uint64 pp = (uint64)p * p;

	invlist[0] = invprod = (uint64)qlist[0] * qlist[0];
	for (i = 1; i < num_q; i++) {
		qq[i] = (uint64)qlist[i] * qlist[i];
		invlist[i] = invprod = montmul64(invprod, qq[i], pp, pp_w);
	}

	invprod = mp_modinv_2(invprod % pp, pp);
	invprod = montmul64(invprod, pp_r, pp, pp_w);
	for (i = num_q - 1; i; i--) {
		invlist[i] = montmul64(invprod, invlist[i-1], pp, pp_w);
		invprod = montmul64(invprod, qq[i], pp, pp_w);
	}
	invlist[i] = invprod;
}

/*------------------------------------------------------------------------*/
//static
uint32
sieve_specialq_64(msieve_obj *obj, poly_search_t *poly, 
		poly_coeff_t *c, int64 sieve_size,
		void *sieve_special_q, void *sieve_p, 
		uint32 special_q_min, uint32 special_q_max, 
		uint32 p_min, uint32 p_max, double deadline, double *elapsed)
{
	/* core search code */

	uint32 i, j;
	uint32 quit = 0;
	p_packed_t *tmp;
	p_packed_var_t specialq_array;
	p_packed_var_t hash_array;
	hash_entry_t *hashtable;
	uint32 num_p;
	uint32 num_roots;
	uint32 hashtable_size_log2;
	uint64 block_size;
	uint64 *invtable = NULL;
	double calc_hashtable_size;
	double cpu_start_time = get_cpu_time();
	struct timeval tstop;	// stop time of this job
	struct timeval tstart;	// start time of this job
	mpz_t qprod;
	uint32 num_sizeopt = 0;
	uint32 last_num_sizeopt = 0;

	//msieve_gettimeofday(&tstart, NULL);
	p_packed_init(&specialq_array);
	p_packed_init(&hash_array);
	mpz_init(qprod);
	hash_array.sieve_size = sieve_size;

	/* build all the arithmetic progressions */

	sieve_fb_reset(sieve_p, p_min, p_max, 1, MAX_ROOTS);
	while (sieve_fb_next(sieve_p, c, store_p_packed,
			&hash_array) != P_SEARCH_DONE) {
		;
	}

	num_p = hash_array.num_p;
	num_roots = hash_array.num_roots;

#if 1
	printf("aprogs: %u entries, %u roots\n", num_p, num_roots);
#endif

	block_size = (uint64)p_min * p_min;
	block_size = MIN(block_size, (uint64)2 * sieve_size);

	/* estimate how many entries each hashtable will contain,
	   in order to allocate just enough memory to hold all the
	   entries */

	calc_hashtable_size = 0;
#ifdef NEW_P_PACKED
	for (i = 0; i < num_p; i++) {

		calc_hashtable_size += (double)block_size / hash_array.packed_array[i].pp *
			hash_array.packed_array[i].num_roots;
	}
#else
	for (i = 0, tmp = hash_array.packed_array; i < num_p; i++) {

		calc_hashtable_size += (double)block_size / tmp->pp *
							tmp->num_roots;

		tmp = p_packed_next(tmp);
	}
#endif

	hashtable_size_log2 = ceil(log(calc_hashtable_size * 5. / 3.) / M_LN2);
	hashtable = (hash_entry_t *)xmalloc(sizeof(hash_entry_t) *
				((uint32)1 << hashtable_size_log2));

	/* handle trivial lattice */
	if (special_q_min == 1) {
		quit = handle_special_q(obj, poly, c, hashtable,
				hashtable_size_log2, &hash_array,
				1, 0, block_size, NULL, &num_sizeopt);
		if (quit || special_q_max == 1)
			goto finished;
	}

	/* invtable caches the inverse of each special-q modulo 
	   each p in hash_array */
	   
	invtable = (uint64 *)xmalloc(num_p * SPECIALQ_BATCH_SIZE * 
					sizeof(uint64));

	sieve_fb_reset(sieve_special_q, special_q_min, special_q_max, 
			1, MAX_ROOTS);

	uint32 total_q = 0;
	uint32 total_r = 0;

	while (1) {
		p_packed_t *qptr = specialq_array.packed_array;
		uint32 num_q;
		uint64 *invtmp;
		uint32 batch_q[SPECIALQ_BATCH_SIZE];
		uint64 batch_q_inv[SPECIALQ_BATCH_SIZE];

		/* allocate the next batch of special-q and all the
		   roots they use */

		p_packed_reset(&specialq_array);
		while (sieve_fb_next(sieve_special_q, c, store_p_packed,
					&specialq_array) != P_SEARCH_DONE) {
			if (specialq_array.num_p == SPECIALQ_BATCH_SIZE)
				break;
		}

		num_q = specialq_array.num_p;
		if (num_q == 0)
			break;

		total_q += num_q;
		total_r += specialq_array.num_roots;

#if 0
		printf("special q: %u entries, %u roots\n", num_q, 
					specialq_array.num_roots);
#endif

		/* multiply the special-q together; this, and the
		   need to reduce the size of invtable, is the main
		   reason the special-q batch size is not bigger */

		mpz_set_ui(qprod, 1);
#ifdef NEW_P_PACKED
		for (i = 0; i < num_q; i++) {
			batch_q[i] = qptr[i].p;
			mpz_mul_ui(qprod, qprod, qptr[i].p);
		}
#else
		for (i = 0, tmp = qptr; i < num_q; i++) {
			batch_q[i] = tmp->p;
			mpz_mul_ui(qprod, qprod, tmp->p);
			tmp = p_packed_next(tmp);
		}
#endif

		/* invert all the special-q at once modulo each p in 
		   hash_array. Because we use batch inversion, if a 
		   given p has factors in common with one of the special-q 
		   then the code below will produce garbage for all of
		   them, so skip all the entries for that one p 
		 
		   Note that we need one inverse for each (p,special_q)
		   pair, *not* one for each root. For degree 4 and 6,
		   when each special_q has many roots, this speeds up the
		   modular inverse phase by much more than
		   SPECIALQ_BATCH_SIZE */
#ifdef NEW_P_PACKED
		for (i = 0, tmp = hash_array.packed_array; i < num_p; i++) {

			if (mp_gcd_1(mpz_tdiv_ui(qprod, tmp[i].p), tmp[i].p) == 1)
				batch_invert(batch_q, num_q, 
						batch_q_inv, tmp[i].p,
						tmp[i].mont_r, tmp[i].mont_w);
			else
				memset(batch_q_inv, 0, sizeof(batch_q_inv));

			invtmp = invtable + i;
			for (j = 0; j < num_q; j++) {
				*invtmp = batch_q_inv[j];
				invtmp += num_p;
			}
			//tmp = p_packed_next(tmp);
		}
#else
		for (i = 0, tmp = hash_array.packed_array; i < num_p; i++) {

			if (mp_gcd_1(mpz_tdiv_ui(qprod, tmp->p), tmp->p) == 1)
				batch_invert(batch_q, num_q,
					batch_q_inv, tmp->p,
					tmp->mont_r, tmp->mont_w);
			else
				memset(batch_q_inv, 0, sizeof(batch_q_inv));

			invtmp = invtable + i;
			for (j = 0; j < num_q; j++) {
				*invtmp = batch_q_inv[j];
				invtmp += num_p;
			}
			tmp = p_packed_next(tmp);
		}
#endif

		/* process (each root of) each special-q in turn */

		invtmp = invtable;
		for (i = 0; i < num_q; i++) {

#ifdef NEW_P_PACKED
			for (j = 0; j < qptr[i].num_roots; j++) {
				quit = handle_special_q(obj, poly, c,
					hashtable, hashtable_size_log2,
					&hash_array, qptr[i].p,
					specialq_array.roots[qptr[i].roots_idx + j].start_offset,
					block_size, invtmp, &num_sizeopt);
				if (quit)
					goto finished;
			}

			//msieve_gettimeofday(&tstop, NULL);
			if (get_cpu_time() - cpu_start_time > deadline) {
				//if (msieve_difftime(&tstart, &tstop) > deadline) {
				quit = 1;
				goto finished;
			}

			invtmp += num_p;
#else
			for (j = 0; j < qptr->num_roots; j++) {
				quit = handle_special_q(obj, poly, c,
						hashtable, hashtable_size_log2,
						&hash_array, qptr->p, 
						qptr->roots[j].start_offset,
						block_size, invtmp, &num_sizeopt);
				if (quit)
					goto finished;
			}

			//msieve_gettimeofday(&tstop, NULL);
			if (get_cpu_time() - cpu_start_time > deadline) {
			//if (msieve_difftime(&tstart, &tstop) > deadline) {
				quit = 1;
				goto finished;
			}

			qptr = p_packed_next(qptr);
			invtmp += num_p;
#endif
		}

		if (((num_sizeopt & ((1 << 10) - 1)) == 0) && (num_sizeopt != last_num_sizeopt))
		{
			last_num_sizeopt = num_sizeopt;
			poly_sizeopt_t* data = (poly_sizeopt_t*)poly->callback_data;

			//msieve_gettimeofday(&tstop, NULL);

			double best_e = data->best_saved_combined_e;
			double time_elapsed = get_cpu_time() - cpu_start_time; //msieve_difftime(&tstart, &tstop);

			if (obj->flags & MSIEVE_FLAG_NFS_POLYSIZE)
			{
				// if running sizeopt as well, display those statistics too.
				gmp_printf("coeff: %Zd, entries: %u, sizeopt: %u, "
					"rootopt: %u, saved: %d, best_e: %1.4e, time: %1.2f sec\r",
					c->high_coeff, total_q, num_sizeopt, data->num_rootopt, data->num_saved,
					best_e, time_elapsed);
			}
			else
			{
				// just display the sieve stats
				gmp_printf("coeff: %Zd, entries: %u, saved for sizeopt: %u, time: %1.2f sec\r",
					c->high_coeff, total_q, num_sizeopt, time_elapsed);
			}
		}
	}

finished:

	printf("\n");
#if 0
	printf("hashtable: %u entries, %5.2lf MB\n", 
			(uint32)1 << hashtable_size_log2,
			(double)sizeof(hash_entry_t) *
				((uint32)1 << hashtable_size_log2) / 1048576);
#endif
	free(invtable);
	free(hashtable);
	p_packed_free(&specialq_array);
	p_packed_free(&hash_array);
	mpz_clear(qprod);
	*elapsed = get_cpu_time() - cpu_start_time;
	return quit;
}

/*------------------------------------------------------------------------*/
double
sieve_lattice_cpu(msieve_obj *obj, poly_search_t *poly, 
		poly_coeff_t *c, double deadline)
{
	uint32 i;
	uint32 degree = poly->degree;
	uint32 num_pieces;
	uint32 p_min, p_max;
	uint32 special_q_min, special_q_max;
	uint32 special_q_min2, special_q_max2;
	uint32 special_q_fb_max;
	uint32 num_passes;
	double p_size_max = c->p_size_max;
	double sieve_bound = c->coeff_max / c->m0;
	double elapsed_total = 0;
	void *sieve_p = sieve_fb_alloc(); 
	void *sieve_special_q = sieve_fb_alloc();

	/* size the problem; we choose p_min so that we will get
	   to use a small number of each progression's offsets in
	   the hashtable search, while respecting the bound on
	   a_{d-2}. Choosing larger p implies that we can use more
	   of its offsets while choosing smaller p means we must
	   use fewer (or sometimes none) of its offsets. Larger p
	   also implies that each pair of p are less likely to
	   yield a collision, which makes the search more
	   difficult */

	p_min = MIN(MAX_OTHER / P_SCALE, sqrt(2.5 / sieve_bound));
	p_min = MIN(p_min, pow((double)SIEVE_MAX / sieve_bound, 1. / 4.) - 1);
	p_min = MIN(p_min, sqrt(p_size_max) / P_SCALE);

	p_max = p_min * P_SCALE;

	special_q_min = 1;
	special_q_max = MIN(MAX_SPECIAL_Q, p_size_max / p_max / p_max);

	/* set up the special q factory; special-q may have 
	   arbitrary factors, but many small factors are 
	   preferred since that will allow for many more roots
	   per special q, so we choose the factors to be as 
	   small as possible */

	special_q_fb_max = MIN(100000, special_q_max);
	sieve_fb_init(sieve_special_q, c,
			5, special_q_fb_max,
			1, degree,
			0);

	/* because special-q can have any factors, we require that
	   the progressions we generate use p that have somewhat
	   large factors. This minimizes the chance that a given
	   special-q has factors in common with many progressions
	   in the set */

	sieve_fb_init(sieve_p, c, 
			35, 5000,
			1, degree,
		       	0);

	/* large search problems can be randomized so that
	   multiple runs over the same range of leading
	   a_d will likely generate different results */

	num_pieces = MIN(450, (double)special_q_max * p_max
				/ log(special_q_max) / log(p_max)
				/ 3e9);

	if (degree == 5)
		num_passes = 2;
	else
		num_passes = 1;

	for (i = 0; i < num_passes; i++) {
		uint32 quit;
		double elapsed = 0;
		int64 sieve_size = MIN(sieve_bound * pow((double)p_min, 4.),
					SIEVE_MAX);

		if (sieve_size == SIEVE_MAX && i > 0)
			break;

		if (num_pieces > 451) { /* randomize the special_q range */
			uint32 piece_length = (special_q_max - special_q_min)
					/ num_pieces;
			uint32 piece = get_rand(&obj->seed1, &obj->seed2)
					% num_pieces;

			printf("randomizing rational coefficient: "
				"using piece #%u of %u\n",
				piece + 1, num_pieces);

			special_q_min2 = special_q_min + piece * piece_length;
			special_q_max2 = special_q_min2 + piece_length;
		}
		else {
			special_q_min2 = special_q_min;
			special_q_max2 = special_q_max;
		}

#if 1
		gmp_printf("coeff %Zd specialq %u - %u other %u - %u\n",
				c->high_coeff,
				special_q_min2, special_q_max2,
				p_min, p_max);
#endif

		quit = sieve_specialq_64(obj, poly, c, sieve_size,
				sieve_special_q, sieve_p,
				special_q_min2, special_q_max2,
				p_min, p_max, deadline, &elapsed);

		elapsed_total += elapsed;
		deadline -= elapsed;

		if (quit || special_q_max < 250 ||
		    p_max >= MAX_OTHER / P_SCALE)
			break;

		p_min = p_max;
		p_max *= P_SCALE;
		special_q_max /= (P_SCALE * P_SCALE);
	}

	sieve_fb_free(sieve_special_q);
	sieve_fb_free(sieve_p);
	return elapsed_total;
}
