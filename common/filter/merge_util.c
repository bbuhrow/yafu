/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: merge_util.c 23 2009-07-20 02:59:07Z jasonp_sf $
--------------------------------------------------------------------*/

#include "filter_priv.h"
#include "merge_util.h"

/*--------------------------------------------------------------------*/
void ideal_list_init(ideal_list_t *ideal_list, 
		uint32 num_ideals, uint32 is_active) {

	uint32 i;

	ideal_list->num_ideals = num_ideals;
	ideal_list->list = (ideal_set_t *)xcalloc((size_t)num_ideals,
						sizeof(ideal_set_t));

	for (i = 0; i < num_ideals; i++) {
		ideal_set_t *entry = ideal_list->list + i;
		entry->active = is_active;
		entry->min_relset_size = 0xffff;
		entry->next = entry;
		entry->prev = entry;
	}
}

/*--------------------------------------------------------------------*/
void ideal_list_free(ideal_list_t *ideal_list) {

	uint32 i;

	for (i = 0; i < ideal_list->num_ideals; i++)
		free(ideal_list->list[i].relsets);
	free(ideal_list->list);
}

/*--------------------------------------------------------------------*/
size_t get_merge_memuse(relation_set_t *relsets, uint32 num_relsets,
			ideal_list_t *ideal_list) {

	uint32 i;
	size_t s = num_relsets * sizeof(relation_set_t) +
		   ideal_list->num_ideals * sizeof(ideal_set_t);

	for (i = 0; i < num_relsets; i++) {
		relation_set_t *r = relsets + i;
		s += sizeof(uint32) * (r->num_large_ideals +
					r->num_relations);
	}
	for (i = 0; i < ideal_list->num_ideals; i++) {
		s += sizeof(uint32) * 
			ideal_list->list[i].num_relsets_alloc;
	}
	return s;
}

/*--------------------------------------------------------------------*/
void heap_init(heap_t *heap) {

	uint32 i;

	heap->num_bins = 256;
	heap->num_ideals = 0;
	heap->next_bin = heap->num_bins;
	heap->worst_bin = (uint32)(-1);
	heap->hashtable = (ideal_set_t *)xcalloc((size_t)heap->num_bins,
						sizeof(ideal_set_t));

	for (i = 0; i < heap->num_bins; i++) {
		ideal_set_t *entry = heap->hashtable + i;
		entry->next = entry;
		entry->prev = entry;
	}
}

/*--------------------------------------------------------------------*/
void heap_free(heap_t *heap) {

	free(heap->hashtable);
}

/*--------------------------------------------------------------------*/
static uint32 heap_compute_key(ideal_set_t *ideal) {

	/* the heap is used to pick the next ideal to be merged,
	   and the key is the maximum amount of fill-in that can
	   occur when the merge takes place. The key computation
	   implements the Markowitz criterion: merging the least 
	   dense matrix row (of weight c) with r-1 other rows, in 
	   order to eliminate one ideal they all have in common, 
	   will cause at most (r-1)*(c-1) extra nonzero entries to 
	   appear in the matrix. This is a quantitative way of 
	   choosing the ideal to merge that will perturb the current
	   matrix the least, and is a general case of the minimum 
	   degree algorithm.

	   Using a heap to track the Markowitz value of every ideal
	   means that the next ideal to eliminate is chosen in
	   constant time, with no approximations. */

	if (ideal->num_relsets == 0)
		return (uint32)(-1);

	return ((uint32)ideal->num_relsets - 1) *
	       ((uint32)ideal->min_relset_size - 1);
}

/*--------------------------------------------------------------------*/
void heap_add_ideal(heap_t *heap, 
			ideal_list_t *ideal_list, 
			uint32 ideal) {

	ideal_set_t *head;
	ideal_set_t *new_node = ideal_list->list + ideal;
	uint32 key = heap_compute_key(new_node);

	if (key == (uint32)(-1)) {
		printf("error: attempted to heapify an empty ideal\n");
		exit(-1);
	}

	/* possibly increase the number of hash buckets in the
	   heap to hold the new key value */

	if (key >= heap->num_bins) {
		uint32 i;
		ideal_set_t *ideal_array = ideal_list->list;
		uint32 new_size = MAX(key + 100, 2 * heap->num_bins);
		ideal_set_t *new_hashtable = (ideal_set_t *)xmalloc(
						new_size * sizeof(ideal_set_t));

		/* transfer chains of ideals with the same
		   key value over to the new hashtable */
		for (i = 0; i < heap->num_bins; i++) {
			ideal_set_t *old_entry = heap->hashtable + i;
			ideal_set_t *new_entry = new_hashtable + i;
			new_entry->next = new_entry;
			new_entry->prev = new_entry;
			if (old_entry->next != old_entry) {
				uint32 j;

				j = old_entry->next - ideal_array;
				new_entry->next = old_entry->next;
				ideal_array[j].prev = new_entry;

				j = old_entry->prev - ideal_array;
				new_entry->prev = old_entry->prev;
				ideal_array[j].next = new_entry;
			}
		}
		
		/* initialize the rest of the entries */

		for (; i < new_size; i++) {
			ideal_set_t *entry = new_hashtable + i;
			entry->next = entry;
			entry->prev = entry;
		}
		free(heap->hashtable);
		heap->hashtable = new_hashtable;
		heap->num_bins = new_size;
	}

	/* add the new ideal */

	head = heap->hashtable + key;
	new_node->prev = head;
	new_node->next = head->next;
	head->next->prev = new_node;
	head->next = new_node;

	/* adjust the best and worst heap bucket pointers */

	heap->num_ideals++;
	heap->next_bin = MIN(heap->next_bin, key);
	if (heap->worst_bin == (uint32)(-1))
		heap->worst_bin = key;
	else
		heap->worst_bin = MAX(heap->worst_bin, key);
}
	

/*--------------------------------------------------------------------*/
void heap_remove_ideal(heap_t *heap, 
			ideal_list_t *ideal_list, 
			uint32 ideal) {

	ideal_set_t *node = ideal_list->list + ideal;
	uint32 key = heap_compute_key(node);
	ideal_set_t *head;

	/* do nothing if the ideal is not connected */

	if (node->next == node)
		return;

	/* remove the ideal from the chain of ideals with
	   this key value */

	heap->num_ideals--;
	node->next->prev = node->prev;
	node->prev->next = node->next;
	node->prev = node;
	node->next = node;

	/* adjust the best and worst heap bucket pointers */

	head = heap->hashtable + heap->next_bin;
	if (head->next == head) {
		uint32 i = key;
		while (i < heap->num_bins && head->next == head) {
			head++;
			i++;
		}
		heap->next_bin = i;
	}

	head = heap->hashtable + heap->worst_bin;
	if (head->next == head) {
		uint32 i = key;
		while ((int32)i >= 0 && head->next == head) {
			head--;
			i--;
		}
		heap->worst_bin = i;
	}
}

/*--------------------------------------------------------------------*/
uint32 heap_add_relset(heap_t *active_heap, 
			heap_t *inactive_heap,
			ideal_list_t *ideal_list, 
			relation_set_t *r,
			uint32 relset_num,
			uint32 min_ideal_weight) {
	
	uint32 i;
	uint32 weight = r->num_small_ideals + r->num_large_ideals;

	r->num_active_ideals = 0;

	/* add the relation set by adding each of its 
	   large ideals to the heap */

	for (i = 0; i < r->num_large_ideals; i++) {
		uint32 ideal = r->data[r->num_relations + i];
		ideal_set_t *ideal_set = ideal_list->list + ideal;
		heap_t *heap = ideal_set->active ? active_heap : inactive_heap;
		uint16 prev_weight = ideal_set->min_relset_size;

		/* pull ideal i off the heap temporarily */

		heap_remove_ideal(heap, ideal_list, ideal);

		/* add the relation set to the list of relation
		   sets containing ideal i */

		if (ideal_set->num_relsets == ideal_set->num_relsets_alloc) {
			ideal_set->num_relsets_alloc += 3;
			ideal_set->relsets = (uint32 *)xrealloc(
						ideal_set->relsets,
						ideal_set->num_relsets_alloc *
						sizeof(uint32));
		}

		ideal_set->relsets[ideal_set->num_relsets++] = relset_num;
		ideal_set->min_relset_size = MIN(prev_weight, weight);

		/* if there are enough relation sets still 
		   containing ideal i, put it back into the heap */

		if (ideal_set->num_relsets > min_ideal_weight)
			heap_add_ideal(heap, ideal_list, ideal);

		if (heap == active_heap)
			r->num_active_ideals++;
	}

	/* return the number of ideals that still 
	   need merging before r can become a cycle */

	return r->num_active_ideals;
}

/*--------------------------------------------------------------------*/
void heap_remove_relset(heap_t *active_heap, 
			heap_t *inactive_heap,
			ideal_list_t *ideal_list,
			relation_set_t *r,
			relation_set_t *relset_array,
			uint32 min_ideal_weight) {
	
	uint32 i, j;
	uint32 weight = r->num_small_ideals + r->num_large_ideals;
	uint32 relset_num = r - relset_array;

	/* remove the relation set by removing each of its 
	   large ideals from the heap */

	for (i = 0; i < r->num_large_ideals; i++) {
		uint32 ideal = r->data[r->num_relations + i];
		ideal_set_t *ideal_set = ideal_list->list + ideal;
		heap_t *heap = ideal_set->active ? active_heap : inactive_heap;
		uint16 prev_weight = ideal_set->min_relset_size;

		/* pull ideal i off the heap temporarily */

		heap_remove_ideal(heap, ideal_list, ideal);

		/* remove the relation set from the list of 
		   relation sets containing ideal i */

		for (j = 0; j < ideal_set->num_relsets; j++) {
			if (ideal_set->relsets[j] == relset_num)
				break;
		}
		ideal_set->relsets[j] = ideal_set->relsets[
					--ideal_set->num_relsets];

		if (ideal_set->num_relsets == 0) {

			/* no more relations with this ideal */

			ideal_set->min_relset_size = 0xffff;
			ideal_set->num_relsets_alloc = 0;
			free(ideal_set->relsets);
			ideal_set->relsets = NULL;
		}
		else {
			/* if r was the lightest relation set, we need 
			   to update the minimum relation set weight for
			   ideal i */

			if (weight == prev_weight) {
				uint32 curr_weight = 0xffff;
				for (j = 0; j < ideal_set->num_relsets; j++) {
					relation_set_t *r2 = relset_array + 
							ideal_set->relsets[j];
					curr_weight = MIN(curr_weight,
							r2->num_small_ideals +
							r2->num_large_ideals);
				}
				ideal_set->min_relset_size = curr_weight;
			}

			/* if there are enough relation sets still
			   containing ideal i, put it back into the heap */

			if (ideal_set->num_relsets > min_ideal_weight)
				heap_add_ideal(heap, ideal_list, ideal);
		}
	}
}

/*--------------------------------------------------------------------*/
uint32 heap_remove_best(heap_t *heap, ideal_list_t *ideal_list) {

	/* remove the heap element with the smallest
	   Markowitz value */

	uint32 ideal;
	uint32 key = heap->next_bin;

	if (key == heap->num_bins)
		return (uint32)(-1);
	
	ideal = heap->hashtable[key].next - ideal_list->list;
	heap_remove_ideal(heap, ideal_list, ideal);
	return ideal;
}
	
/*--------------------------------------------------------------------*/
uint32 heap_remove_worst(heap_t *heap, ideal_list_t *ideal_list) {

	/* remove the heap element with the largest
	   Markowitz value */

	uint32 ideal;
	uint32 key = heap->worst_bin;

	if ((int32)key < 0)
		return (uint32)(-1);
	
	ideal = heap->hashtable[key].next - ideal_list->list;
	heap_remove_ideal(heap, ideal_list, ideal);
	return ideal;
}
	
/*--------------------------------------------------------------------*/
void load_next_relset_group(merge_aux_t *aux,
			heap_t *active_heap, heap_t *inactive_heap, 
			ideal_list_t *ideal_list,
			relation_set_t *relset_array,
			uint32 ideal,
			uint32 min_ideal_weight) {

	uint32 i;
	ideal_set_t *ideal_set = ideal_list->list + ideal;

	/* make aux ready to receive a batch of relation sets */

	aux->num_relsets = ideal_set->num_relsets;
	if (aux->num_relsets >= MERGE_MAX_OBJECTS) {
		printf("error: number of relsets too large\n");
		exit(-1);
	}

	/* copy each relation set containing 'ideal', and
	   remember the index into the full array of
	   relation sets where these occur. We'll need this
	   information in order to put everything back when
	   processing of aux has finished */

	for (i = 0; i < ideal_set->num_relsets; i++) {
		uint32 relset_num = ideal_set->relsets[i];
		aux->tmp_relset_idx[i] = relset_num;
		aux->tmp_relsets[i] = relset_array[relset_num];
	}

	/* now that all the relation sets are safely copied,
	   go back and remove the originals from the heap,
	   wiping out the old relation set structures in
	   the process */

	for (i = 0; i < aux->num_relsets; i++) {
		uint32 relset_num = aux->tmp_relset_idx[i];
		relation_set_t *r = relset_array + relset_num;
		heap_remove_relset(active_heap, inactive_heap,
				ideal_list, r, relset_array,
				min_ideal_weight);
		memset(r, 0, sizeof(relation_set_t));
	}
}

/*--------------------------------------------------------------------*/
void bury_inactive_ideal(relation_set_t *relset_array,
			ideal_list_t *ideal_list, uint32 ideal) {

	uint32 i, j;
	ideal_set_t *ideal_set = ideal_list->list + ideal;

	for (i = 0; i < ideal_set->num_relsets; i++) {
		relation_set_t *r = relset_array + ideal_set->relsets[i];
		uint32 num_ideals = r->num_large_ideals;

		r->num_small_ideals++;
		if (--(r->num_large_ideals) > 0) {

			uint32 *r_array = r->data + r->num_relations;

			/* squeeze the ideal out of r. Because the total 
			   weight of the relation set is unchanged, we 
			   don't have to re-heapify any of the other 
			   ideals in r */

			for (j = 0; j < num_ideals; j++) {
				if (r_array[j] == ideal)
					break;
			}
			for (j++; j < num_ideals; j++) {
				r_array[j-1] = r_array[j];
			}
		}
	}

	/* reset the ideal structure */

	free(ideal_set->relsets);
	memset(ideal_set, 0, sizeof(ideal_set_t));
}

/*--------------------------------------------------------------------*/
void merge_two_relsets(relation_set_t *r1, relation_set_t *r2, 
			relation_set_t *r_out, merge_aux_t *aux) {

	uint32 i;
	uint32 max_relations = r1->num_relations + r2->num_relations;
	uint32 max_ideals = r1->num_large_ideals + r2->num_large_ideals - 2;

	memset(r_out, 0, sizeof(relation_set_t));

	r_out->num_small_ideals = r1->num_small_ideals + 
				  r2->num_small_ideals;

	/* combine the list of relations in r1 and r2
	   using a merge operation */

	if (max_relations >= MERGE_MAX_OBJECTS) {
		printf("error: relation list too large\n");
		exit(-1);
	}

	i = merge_relations(aux->tmp_relations,
				r1->data, r1->num_relations,
				r2->data, r2->num_relations);
	if (i == 0) {
		memset(r_out, 0, sizeof(relation_set_t));
		return;
	}
	r_out->num_relations = i;

	/* merge the ideal lists of r1 and r2 */

	if (max_ideals >= MERGE_MAX_OBJECTS) {
		printf("error: list of merged ideals too large\n");
		exit(-1);
	}

	i = merge_relations(aux->tmp_ideals,
				r1->data + r1->num_relations, 
				r1->num_large_ideals,
				r2->data + r2->num_relations, 
				r2->num_large_ideals);

	r_out->num_large_ideals = i;

	/* save the merged lists */

	r_out->data = (uint32 *)xmalloc(sizeof(uint32) *
					(r_out->num_relations + 
					 r_out->num_large_ideals));
	memcpy(r_out->data, 
	       aux->tmp_relations, 
	       r_out->num_relations * sizeof(uint32));
	memcpy(r_out->data + r_out->num_relations, 
	       aux->tmp_ideals, 
	       r_out->num_large_ideals * sizeof(uint32));
}

/*--------------------------------------------------------------------*/
uint32 estimate_new_weight(relation_set_t *r1,
			relation_set_t *r2) {

	/* guess the total number of ideals derived from
	   merging r1 and r2. This is much faster than
	   actually merging r1 and r2, but must occur much
	   more often */

	uint32 i, j, k;
	uint32 num_small;
	uint32 num_ideals1 = r1->num_large_ideals;
	uint32 num_ideals2 = r2->num_large_ideals;
	uint32 *ilist1 = r1->data + r1->num_relations;
	uint32 *ilist2 = r2->data + r2->num_relations;

	/* count the number of large ideals that would
	   be left if r1 and r2 were merged */

	i = j = k = 0;
	while (i < num_ideals1 && j < num_ideals2) {
		uint32 ideal1 = ilist1[i];
		uint32 ideal2 = ilist2[j];
		if (ideal1 < ideal2) {
			i++; k++;
		}
		else if (ideal1 > ideal2) {
			j++; k++;
		}
		else {
			i++; j++;
		}
	}

	num_small = r1->num_small_ideals + r2->num_small_ideals;
	return k + num_small + (num_ideals1 - i) + (num_ideals2 - j);
}

