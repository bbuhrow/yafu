/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: hashtable.c 713 2012-03-10 12:56:24Z jasonp_sf $
--------------------------------------------------------------------*/

#include <ms_common.h>

/*--------------------------------------------------------------------*/
void hashtable_init(hashtable_t *h, uint32 blob_words, uint32 hash_words) {

	uint32 init_match_size = 5000;
	uint32 log2_hashtable_size = 14;

	h->blob_words = blob_words;
	h->hash_words = hash_words;
	if (hash_words == 0)
		h->hash_words = blob_words;

	h->log2_hashtable_size = log2_hashtable_size;
	h->hashtable = (uint32 *)xcalloc((size_t)1 << log2_hashtable_size,
					sizeof(uint32));
	h->num_used = 0;
	h->congestion_target = 0.8 * (1 << log2_hashtable_size);

	h->match_array_size = 1;
	h->match_array_alloc = init_match_size;
	h->match_array = (uint32 *)xmalloc(init_match_size * 
					(blob_words + 1) * sizeof(uint32));
}

/*--------------------------------------------------------------------*/
void hashtable_free(hashtable_t *h) {
	free(h->hashtable);
	free(h->match_array);
	h->hashtable = NULL;
	h->match_array = NULL;
}

/*--------------------------------------------------------------------*/
size_t hashtable_sizeof(hashtable_t *h) {
	return (sizeof(uint32) << h->log2_hashtable_size) +
		((h->blob_words + 1) * sizeof(uint32) * 
		 	h->match_array_alloc);
}

/*--------------------------------------------------------------------*/
void hashtable_close(hashtable_t *h) {
	free(h->hashtable);
	h->hashtable = NULL;

	h->match_array = (uint32 *)xrealloc(h->match_array,
					h->match_array_size *
					(h->blob_words + 1) *
					sizeof(uint32));
	h->match_array_alloc = h->match_array_size;
}

/*--------------------------------------------------------------------*/
void *hashtable_find(hashtable_t *h, void *blob, 
		     uint32 *ordinal_id, uint32 *present) {

	uint32 i;
	uint32 offset, hashval;
	uint32 *key = (uint32 *)blob;
	uint32 *entry = NULL;
	uint32 *hashtable = h->hashtable;
	uint32 *match_array = h->match_array;
	uint32 log2_hashtable_size = h->log2_hashtable_size;
	uint32 blob_words = h->blob_words;
	uint32 hash_words = h->hash_words;

	/* pre-emptively increase the size allocated
	   for hashtable matches */
	
	if (h->match_array_size + 1 >= h->match_array_alloc) {
		h->match_array_alloc *= 2;
		match_array = h->match_array = (uint32 *)xrealloc(
						h->match_array,
						sizeof(uint32) *
						h->match_array_alloc *
						(blob_words + 1));
	}

	/* pre-emptively grow the hashtable size if it's
	   getting congested */

	if (h->num_used + 1 >= h->congestion_target &&
			log2_hashtable_size < 28) {

		/* reset the hashtable array */

		log2_hashtable_size++;
		free(hashtable);
		hashtable = h->hashtable = (uint32 *)xcalloc((size_t)1 << 
						log2_hashtable_size,
						sizeof(uint32));
		h->log2_hashtable_size = log2_hashtable_size;
		h->congestion_target = 0.8 * (1 << log2_hashtable_size);
		h->num_used = 0;

		/* re-link all the previously found entries */

		key = h->match_array + blob_words + 1;

		for (i = 1; i < h->match_array_size; i++) {
			hashval = hash_function(key, hash_words);
			hashval = hashval >> (32 - log2_hashtable_size);

			key[blob_words] = hashtable[hashval];
			if (key[blob_words] == 0)
				h->num_used++;
			hashtable[hashval] = i;
			key += blob_words + 1;
		}
	}

	/* compute hash value; only the first one or two words
	   figure into it */

	key = (uint32 *)blob;
	hashval = hash_function(key, hash_words);
	hashval = hashval >> (32 - log2_hashtable_size);

	/* attempt to find it in the table; all hash_words
	   words must match for lookup to succeed */

	offset = hashtable[hashval];
	while (offset != 0) {
		entry = match_array + (size_t)offset * (blob_words + 1);
		for (i = 0; i < hash_words; i++) {
			if (entry[i] != key[i])
				break;
		}
		if (i == hash_words)
			break;
		offset = entry[blob_words];
	}

	if (ordinal_id) {
		if (offset == 0)
			*ordinal_id = h->match_array_size - 1;
		else
			*ordinal_id = offset - 1;
	}

	if (offset == 0) {

		/* not found; add it */

		entry = match_array + 
			(size_t)h->match_array_size * (blob_words + 1);

		for (i = 0; i < hash_words; i++) {
			entry[i] = key[i];
		}
		entry[blob_words] = hashtable[hashval]; /* add to linked list */
		if (entry[blob_words] == 0)
			h->num_used++;
		hashtable[hashval] = h->match_array_size++;
	}
	if (present)
		*present = offset;

	return entry;
}
