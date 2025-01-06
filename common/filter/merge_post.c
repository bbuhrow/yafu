/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: merge_post.c 23 2009-07-20 02:59:07Z jasonp_sf $
--------------------------------------------------------------------*/

#include "filter_priv.h"
#include "merge_util.h"

/* This code implements a postprocessing phase for the collection
   of cycles that is the output of the merge phase of NFS filtering.
   It uses a highly simplified form of an algorithm by Denny and Muller
   that everyone else seems to have completely forgotten.

   Suppose you have two cycles, one containing relations (a,b,c) and 
   the other containing relations (a,b,d,e). These two cycles are 
   entirely equivalent to two other cycles, (a,b,c) and (c,d,e), except
   that the latter is one relation smaller than the former. The algorithm
   basically takes advantage of the fact that occaisionally two cycles
   have many relations in common, and if the common relations in the cycles
   outnumber the unique relations in one of the cycles, then the 'cycle
   basis' can be changed to skip the common relations some of the time.

   For NFS filtering that is not smart (i.e. that doesn't use Cavallar-
   like algorithms), the cycle basis reduction process can prune a huge
   number of relations from the collection of cycles. Here, however,
   most cycles are quite small, so the postprocessing amounts to pruning
   a small number of relations (2-4%) from the collection of cycles.

   Many of the tricks that Denny and Muller come up with to make the
   pruning more effective are completely unneeded here; we can instead
   use a variant of Cavallar's minimum spanning tree algorithm to rigorously
   find the basis with smallest relation weight */

/*--------------------------------------------------------------------*/
static int compare_uint32(const void *x, const void *y) {
	uint32 *xx = (uint32 *)x;
	uint32 *yy = (uint32 *)y;
	if (*xx > *yy)
		return 1;
	if (*xx < *yy)
		return -1;
	return 0;
}

/*--------------------------------------------------------------------*/
static void copy_relset(relation_set_t *src, relation_set_t *dst) {

	uint32 size = src->num_relations + src->num_large_ideals;

	*dst = *src;
	dst->data = (uint32 *)xmalloc(size * sizeof(uint32));
	memcpy(dst->data, src->data, size * sizeof(uint32));
}

/*--------------------------------------------------------------------*/
static uint32 merge_via_spanning_tree(merge_aux_t *aux) {
	
	/* given a collection of N cycles with at least one 
	   relation in common, use a minimum spanning tree 
	   process to produce N different cycles that contain
	   fewer total relations */

	uint32 i, j, k;
	uint32 edges[SPANNING_TREE_MAX_RELSETS]
		    [SPANNING_TREE_MAX_RELSETS];
	uint32 vertex[SPANNING_TREE_MAX_RELSETS];
	uint32 v_done;
	uint32 num_relsets = aux->num_relsets;
	relation_set_t *relsets = aux->tmp_relsets;
	relation_set_t *tmp_relsets = aux->tmp_relsets2;

	uint32 min_index = 0;
	uint32 min_weight = relsets[0].num_large_ideals;
	uint32 total_weight = relsets[0].num_large_ideals;
	uint32 new_total_weight;

	/* count the starting number of relations and find the
	   lightest cycle. This cycle forms the 'anchor' for the
	   minimum spanning tree. The choice of anchor cycle is
	   arbitrary, but because it will appear unmodified in 
	   the final list of cycles we want it to have few relations */

	for (i = 1; i < num_relsets; i++) {
		uint32 num_ideals = relsets[i].num_large_ideals;
		if (num_ideals < min_weight) {
			min_weight = num_ideals;
			min_index = i;
		}
		total_weight += num_ideals;
	}

	/* calculate the number of relations derived from
	   merging each cycle with all the others */

	for (i = k = 0; i < num_relsets - 1; i++) {
		for (j = i + 1; j < num_relsets; j++, k++) {
			edges[i][j] = edges[j][i] = estimate_new_weight(
						relsets + i, relsets + j);
		}
	}

	memcpy(tmp_relsets, relsets, 
			num_relsets * sizeof(relation_set_t));

	/* permute the list of cycles so that the anchor
	   cycle is first, and copy the anchor to the output */

	vertex[0] = min_index;
	for (i = 0, j = 1; i < num_relsets; i++) {
		if (i != min_index)
			vertex[j++] = i;
	}

	copy_relset(tmp_relsets + min_index, relsets + 0);
	new_total_weight = relsets[0].num_large_ideals;

	/* proceed with Prim's algorithm. Note that here there is
	   a twist: the objective is not to connect cycles via 
	   merges so that an ideal is eliminated, but to minimize 
	   the total weight. This means that we have an alternative to
	   merging two cycles: adding a cycle to the tree but not
	   combining it with another cycle. Hence we technically do
	   not form a spanning tree but a 'spanning archipelago'; in
	   fact, most often this process will not change the input
	   list of cycles at all */

	for (v_done = 1; v_done < num_relsets; v_done++) {
		uint32 index1 = 0;
	        uint32 index2 = 0;
		uint32 weight = (uint32)(-1);

		/* find the pair of cycles, one in the
		   tree and one not yet added, that has the 
		   minimum merge weight */

		for (i = 0; i < v_done; i++) {
			for (j = v_done; j < num_relsets; j++) {
				uint32 curr_weight = edges[vertex[i]]
							  [vertex[j]];
				if (curr_weight < weight) {
					index1 = i;
					index2 = j;
					weight = curr_weight;
				}
			}
		}

		/* merge the cycles only if the merged version has
		   fewer relations than the original. Either way,
		   add the new cycle to the tree */

		i = vertex[index1];
		j = vertex[index2];
		if (edges[i][j] < tmp_relsets[j].num_large_ideals) {
			merge_two_relsets(tmp_relsets + i, tmp_relsets + j,
					  relsets + v_done, aux);
		}
		else {
			copy_relset(tmp_relsets + j, relsets + v_done);
		}
		new_total_weight += relsets[v_done].num_large_ideals;

		i = vertex[index2];
		vertex[index2] = vertex[v_done];
		vertex[v_done] = i;
	}

	/* remove the original list of cycles */

	for (i = 0; i < num_relsets; i++)
		free(tmp_relsets[i].data);

	/* return the number of relations removed */

	return total_weight - new_total_weight;
}

/*--------------------------------------------------------------------*/
static uint32 merge_via_pivot(merge_aux_t *aux) {
	
	/* simplified form of the pevious routine */

	uint32 i;
	relation_set_t *relsets = aux->tmp_relsets;
	uint32 num_relsets = aux->num_relsets;
	uint32 weight;
	uint32 pivot_idx;
	uint32 total_weight = 0;
	uint32 new_total_weight = 0;

	/* find the lightest cycle */

	weight = relsets[0].num_large_ideals;
	total_weight = weight;
	pivot_idx = 0;
	for (i = 1; i < num_relsets; i++) {
		uint32 new_weight = relsets[i].num_large_ideals;
		total_weight += new_weight;
		if (new_weight < weight) {
			pivot_idx = i;
			weight = new_weight;
		}
	}
		
	/* merge it with any other cycle for which the total
	   merged weight would be reduced */

	for (i = 0; i < num_relsets; i++) {
		if (i != pivot_idx && 
		    estimate_new_weight(relsets + pivot_idx, relsets + i) < 
					relsets[i].num_large_ideals) {

			relation_set_t tmp = relsets[i];
			merge_two_relsets(relsets + pivot_idx, 
					&tmp, relsets + i, aux);

			free(tmp.data);
		}
		new_total_weight += relsets[i].num_large_ideals;
	}

	/* return number of relations pruned */

	return total_weight - new_total_weight;
}

/*--------------------------------------------------------------------*/
static uint32 do_merges_core(merge_aux_t *aux) {
	
	/* the same number of cycles are returned, along with
	   the number of relations that the shuffling of cycle
	   basis has saved */

	if (aux->num_relsets == 1) {
		return 0;
	}
	else if (aux->num_relsets >= 3 && 
		 aux->num_relsets <= SPANNING_TREE_MAX_RELSETS) { 
		return merge_via_spanning_tree(aux);
	}
	else {
		return merge_via_pivot(aux);
	}
}

/*--------------------------------------------------------------------*/
static void store_next_relset_group(merge_aux_t *aux,
			heap_t *active_heap, heap_t *inactive_heap, 
			ideal_list_t *ideal_list,
			relation_set_t *relset_array) {
	uint32 i;

	/* store all the cycles just processed. Relations get
	   added back into the heap if the relation occurs in
	   more than one cycle (i.e. is a candidate for later
	   optimization */

	for (i = 0; i < aux->num_relsets; i++) {
		uint32 relset_num = aux->tmp_relset_idx[i];
		relation_set_t *r = relset_array + relset_num;

		*r = aux->tmp_relsets[i];
		heap_add_relset(active_heap, inactive_heap,
				ideal_list, r, relset_num, 1);
	}
}

/*--------------------------------------------------------------------*/
#define LOG2_REL_HASHTABLE_SIZE 23
#define NUM_CYCLE_BINS 9

void filter_postproc_relsets(msieve_obj *obj, merge_t *merge) {

	/* external interface to the postprocessing phase */

	uint32 i, j;
	uint32 cycle_bins[NUM_CYCLE_BINS + 2] = {0};
	uint32 max_ideal_list_size;
	uint32 num_relations;
	relation_set_t *relset_array = merge->relset_array;
	uint32 num_relsets = merge->num_relsets;
	uint32 num_ideals;
	hashtable_t h;
	ideal_list_t ideal_list;
	heap_t active_heap;
	heap_t inactive_heap;
	merge_aux_t *aux;

	logprintf(obj, "commencing cycle optimization\n");

	/* we will call the relation sets 'cycles' here, but they
	   are treated the same as in the merge phase. Instead of a
	   list of ideals, each cycle has a list of relations that
	   occur in that cycle. We assign one 'ideal' to each relation
	   that appears in the collection of cycles, and minimize that.
	   Note that we assume that the merge phase has already 
	   cancelled out any relations that appear more than once in
	   a cycle */

	hashtable_init(&h, 1, 0);
	max_ideal_list_size = 0;
	for (i = num_relations = 0; i < num_relsets; i++) {
		relation_set_t *r = relset_array + i;
		uint32 curr_num_relations = r->num_relations;

		max_ideal_list_size = MAX(max_ideal_list_size, 
					  curr_num_relations);
		num_relations += curr_num_relations;

		/* initialize the cycle. There are no small ideals,
		   and cycles that have only one relation will never
		   be optimized to have zero relations, so skip those */

		r->num_small_ideals = 0;
		r->num_large_ideals = 0;
		r->num_active_ideals = 0;
		if (curr_num_relations == 1)
			continue;

		r->num_large_ideals = curr_num_relations;
		r->num_active_ideals = curr_num_relations;
		r->data = (uint32 *)xrealloc(r->data, 2 * curr_num_relations * 
						sizeof(uint32));

		/* every unique relation is assigned a number */

		for (j = 0; j < curr_num_relations; j++) {
			hashtable_find(&h, r->data + j, 
					r->data + curr_num_relations + j, 
					NULL);
		}

		/* relation lists occur in sorted order */

		qsort(r->data + curr_num_relations, 
			(size_t)curr_num_relations,
			sizeof(uint32), compare_uint32);
	}
	num_ideals = merge->num_ideals = hashtable_get_num(&h);
	hashtable_free(&h);
	logprintf(obj, "start with %u relations\n", num_relations);

	/* initialize the heaps for tracking cycles; all cycles
	   start out active. Then add each cycle to the heaps;
	   each relation that occurs in more than one cycle is
	   heapified */

	aux = (merge_aux_t *)xmalloc(sizeof(merge_aux_t));
	heap_init(&active_heap);
	heap_init(&inactive_heap);
	ideal_list_init(&ideal_list, num_ideals, 1);

	for (i = 0; i < num_relsets; i++) {
		relation_set_t *r = relset_array + i;

		if (r->num_large_ideals > 1) {
			heap_add_relset(&active_heap, &inactive_heap,
					&ideal_list, r, i, 1);
		}
	}

	/* for each relation, starting with the relations that
	   occur in the fewst cycles,
	   - remove all the cycles containing that relation
	   - optimize the collection of cycles
	   - put the cycles back into the heap
	   - make the relation inactive
	*/

	num_relations = 0;
	while (active_heap.num_ideals > 0) {

		uint32 ideal = heap_remove_best(&active_heap, 
						&ideal_list);

		load_next_relset_group(aux, &active_heap, 
				&inactive_heap, &ideal_list, 
				relset_array, ideal, 1);

		num_relations += do_merges_core(aux);

		ideal_list.list[ideal].active = 0;

		store_next_relset_group(aux, &active_heap, 
					&inactive_heap, &ideal_list, 
					relset_array);
	}

	logprintf(obj, "pruned %u relations\n", num_relations);
	logprintf(obj, "memory use: %.1f MB\n", (double)
			get_merge_memuse(relset_array, num_relsets,
						&ideal_list) / 1048576);

	/* print statistics on the final collection of cycles */

	for (i = j = 0; i < num_relsets; i++) {
		relation_set_t *r = relset_array + i;
		j = MAX(j, r->num_relations);

		if (r->num_relations > NUM_CYCLE_BINS)
			cycle_bins[NUM_CYCLE_BINS+1]++;
		else
			cycle_bins[r->num_relations]++;

		r->num_small_ideals = 0;
		r->num_large_ideals = 0;
	}
	logprintf(obj, "distribution of cycle lengths:\n");
	for (i = 1; i < NUM_CYCLE_BINS + 1; i++)
		logprintf(obj, "%u relations: %u\n", i, cycle_bins[i]);
	logprintf(obj, "%u+ relations: %u\n", i, cycle_bins[i]);
	logprintf(obj, "heaviest cycle: %u relations\n", j);

	heap_free(&active_heap);
	heap_free(&inactive_heap);
	ideal_list_free(&ideal_list);
	free(aux);
}
