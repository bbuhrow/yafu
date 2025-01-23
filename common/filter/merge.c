/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: merge.c 1020 2018-05-18 13:33:21Z jasonp_sf $
--------------------------------------------------------------------*/

#include "filter_priv.h"
#include "merge_util.h"

/*--------------------------------------------------------------------*/
/* structure representing the weights of all the cycles
   found in the merge phase. This includes an interface
   for retrieving the weight of the lightest few cycles */

typedef struct {
	uint32 num_weights_alloc;
	uint32 *weights;
} matrix_weight_t;

static void matrix_weight_init(matrix_weight_t *m) {
	m->num_weights_alloc = 100;
	m->weights = (uint32 *)xcalloc((size_t)m->num_weights_alloc,
					sizeof(uint32));
}

static void matrix_weight_free(matrix_weight_t *m) {
	free(m->weights);
}

static void matrix_weight_add(matrix_weight_t *m, relation_set_t *r) {

	uint32 w = r->num_small_ideals +
		   r->num_large_ideals;
	if (w >= m->num_weights_alloc) {
		uint32 i;
		m->weights = (uint32 *)xrealloc(m->weights, 
					(size_t)(w + 100) *
					sizeof(uint32));
		for (i = m->num_weights_alloc; i < w + 100; i++)
			m->weights[i] = 0;
		m->num_weights_alloc = w + 100;
	}
	m->weights[w]++;
}

static void matrix_weight_sub(matrix_weight_t *m, relation_set_t *r) {

	uint32 w = r->num_small_ideals +
		   r->num_large_ideals;
	m->weights[w]--;
}

static uint64 get_matrix_weight(matrix_weight_t *m, uint32 target_cycles) {

	uint32 i;
	uint32 total_cycles = 0;
	uint64 total_weight = 0;

	for (i = 0; total_cycles < target_cycles &&
				i < m->num_weights_alloc; i++) {

		uint32 curr_cycles = MIN(m->weights[i], 
					target_cycles - total_cycles);
		total_weight += i * curr_cycles;
		total_cycles += curr_cycles;
	}

	return total_weight;
}

/*--------------------------------------------------------------------*/
static int compare_relsets(const void *x, const void *y) {
	relation_set_t *xx = (relation_set_t *)x;
	relation_set_t *yy = (relation_set_t *)y;
	return (int)xx->num_small_ideals - (int)yy->num_small_ideals;
}

/*--------------------------------------------------------------------*/
static void merge_via_spanning_tree(merge_aux_t *aux) {
	
	/* merge the relation sets in aux using the minimum
	   spanning tree algorithm as described in Cavallar's
	   paper. 
	   
	   Minimum spanning trees are an active research area.
	   We use Prim's algorithm to find the spanning tree;
	   it's very easy to implement, and these spanning tree
	   problems are quite small so asymptotic behavior
	   really doesn't matter all that much */

	uint32 i, j;
	uint32 edges[SPANNING_TREE_MAX_RELSETS]
		    [SPANNING_TREE_MAX_RELSETS];
	uint32 vertex[SPANNING_TREE_MAX_RELSETS];
	uint32 v_done;
	uint32 num_relsets = aux->num_relsets;
	relation_set_t *relsets = aux->tmp_relsets;
	relation_set_t *tmp_relsets = aux->tmp_relsets2;

	/* we model the collection of relation sets as a fully
	   connected graph, with one vertex for each relation set
	   and an edge for each possible merge between two relation
	   sets. The objective is to find the (num_relsets-1) edges
	   with minimum weight that connect all num_relsets vertices */

	for (i = 0; i < num_relsets; i++)
		vertex[i] = i;

	/* estimate the cost of merging every relation set
	   with every other */

	for (i = 0; i < num_relsets - 1; i++) {
		for (j = i + 1; j < num_relsets; j++) {
			edges[i][j] = edges[j][i] = estimate_new_weight(
						relsets + i, relsets + j);
		}
	}
	memcpy(tmp_relsets, relsets, 
			num_relsets * sizeof(relation_set_t));

	/* start with relation set 0 and add every relation set
	   to the spanning tree in turn */

	for (v_done = 1; v_done < num_relsets; v_done++) {
		uint32 index1 = 0;
	        uint32 index2 = 0;
		uint32 weight = (uint32)(-1);

		/* find the pair of relation sets, one in the
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

		/* do the merge, add the new relation set to the
		   list of those that are in the tree */

		merge_two_relsets(tmp_relsets + vertex[index1],
				  tmp_relsets + vertex[index2],
				  relsets + v_done - 1, aux);

		i = vertex[index2];
		vertex[index2] = vertex[v_done];
		vertex[v_done] = i;
	}

	/* remove the original collection of relation sets */

	for (i = 0; i < num_relsets; i++)
		free(tmp_relsets[i].data);
}

/*--------------------------------------------------------------------*/
static void merge_via_pivot(merge_aux_t *aux) {
	
	/* a more conventional merge operation for relation
	   sets: find the lightest one and merge it into
	   all the other relation sets */

	uint32 i;
	relation_set_t *relsets = aux->tmp_relsets;
	uint32 num_relsets = aux->num_relsets;
	uint32 weight;
	relation_set_t pivot;
	uint32 pivot_idx;

	/* find the cheapest merge */

	weight = relsets[0].num_small_ideals +
		 relsets[0].num_large_ideals;
	pivot = relsets[0];
	pivot_idx = 0;
	for (i = 1; i < num_relsets; i++) {
		uint32 new_weight = relsets[i].num_small_ideals +
					relsets[i].num_large_ideals;
		if (new_weight < weight) {
			pivot = relsets[i];
			pivot_idx = i;
			weight = new_weight;
		}
	}

	/* squeeze the pivot relation set out of the list */

	relsets[pivot_idx] = relsets[i-1];
		
	/* merge with the other relation sets and free the pivot */

	for (i = 0; i < num_relsets - 1; i++) {
		relation_set_t tmp = relsets[i];
		merge_two_relsets(&pivot, &tmp, relsets + i, aux);
		free(tmp.data);
	}

	free(pivot.data);
}

/*--------------------------------------------------------------------*/
static void do_merges_core(merge_aux_t *aux) {
	
	/* main interface for merging an ideal so that it
	   does not appear in a group of relation sets */

	if (aux->num_relsets == 1) {
		/* relation set contains a singleton ideal; delete it */
		free(aux->tmp_relsets[0].data);
		memset(aux->tmp_relsets + 0, 0, sizeof(relation_set_t));
		return;
	}
	else if (aux->num_relsets >= 3 && 
		 aux->num_relsets <= SPANNING_TREE_MAX_RELSETS) { 
		merge_via_spanning_tree(aux);
	}
	else {
		/* use for 2-way merges, and merges that are too large */
		merge_via_pivot(aux);
	}
}

/*--------------------------------------------------------------------*/
static void toggle_ideal_state(ideal_list_t *ideal_list, uint32 ideal, 
				relation_set_t *relset_array,
				uint32 *num_cycles,
				matrix_weight_t *mat_weight) {
	
	/* make an inactive ideal active (i.e. allow it to be 
	   merged), or make an active ideal inert (do not allow
	   it to be merged) */

	uint32 i;
	ideal_set_t *ideal_set = ideal_list->list + ideal;
	uint32 active = ideal_set->active;

	/* for each relation set containing the ideal */

	for (i = 0; i < ideal_set->num_relsets; i++) {
		relation_set_t *r = relset_array + ideal_set->relsets[i];

		/* relation sets that contained this ideal, but
		   had no other ideals active, now are not cycles
		   anymore (i.e. ideal needs to be merged in
		   order to make r a cycle again). Conversely, if 
		   the last active ideal in r is this one, then
		   r is now a cycle */

		if (active) {
			if (--(r->num_active_ideals) == 0) {
				(*num_cycles)++;
				matrix_weight_add(mat_weight, r);
			}
		}
		else {
			if (r->num_active_ideals++ == 0) {
				(*num_cycles)--;
				matrix_weight_sub(mat_weight, r);
			}
		}
	}

	/* toggle the ideal state */

	ideal_set->active ^= 1;
}

/*--------------------------------------------------------------------*/
static uint32 store_next_relset_group(merge_aux_t *aux,
			heap_t *active_heap, heap_t *inactive_heap, 
			ideal_list_t *ideal_list,
			relation_set_t *relset_array,
			matrix_weight_t *mat_weight) {

	/* add a collection of relation sets back into the main
	   collection of ideals processed by the merge phase */

	uint32 i;
	uint32 num_cycles = 0;

	/* for each relation set */

	for (i = 0; i < aux->num_relsets - 1; i++) {
		uint32 relset_num = aux->tmp_relset_idx[i];
		relation_set_t *r = relset_array + relset_num;

		/* re-insert into the main list of relation sets */

		*r = aux->tmp_relsets[i];

		/* if the merging wiped out r (i.e. it had a singleton 
		   ideal), or there are too many relations in r, give up */

		if (r->num_relations == 0) {
			continue;
		}
		else if (r->num_relations > MAX_RELSET_SIZE) {
			free(r->data);
			memset(r, 0, sizeof(relation_set_t));
			continue;
		}

		/* add r into the merge heaps; if r does not have
		   any active ideals, then it is a cycle */

		if (heap_add_relset(active_heap, inactive_heap,
					ideal_list, r, relset_num, 0) == 0) {
			num_cycles++;
			matrix_weight_add(mat_weight, r);
		}
	}

	return num_cycles;
}

/*--------------------------------------------------------------------*/
#define NUM_CYCLE_BINS 9

int32 filter_merge_full(msieve_obj *obj, merge_t *merge, uint32 min_cycles) {

	/* The core of the NFS merge phase; this routine starts
	   with a collection of relation sets, each containing entries
	   drawn from a pool of large ideals, and produces another
	   collection of relation sets that can function as an
	   NFS matrix. More formally, given R relation sets containing
	   I large ideals, we can find at least R-I cycles, and probably
	   much more than that. The exact number will be near R-I+C, where
	   C is the number of connected components in the graph formed
	   by R; we don't know what C is, but it should be a small fraction
	   of R-I.

	   min_cycles is the minimum number of cycles that we absolutely
	   must have in order to make a usable matrix. It is the combined
	   number of ideals in the NFS rational and algebraic factor
	   bases that were never tracked in the rest of the filtering phase,
	   along with a little extra to allow discarding very heavy cycles.
	   Producing a matrix with exactly min_cycles cycles means we will
	   have to merge many ideals from I, which leads to a very dense
	   matrix. On the other hand, R-I can be much larger than min_cycles,
	   and we can instead attempt to produce a number of cycles close
	   to R-I. This reduces the number of merges we would have to do,
	   and leaves us with a very large but very light matrix. Between
	   these two extremes there is a range of 'compromise' matrices.
	   A compromise matrix will ignore some number S of the large ideals
	   in I, and consist of a moderate number (min_cycles+S) of cycles
	   that are moderately dense. The choice of S and resulting cycles
	   depends on how many excess relation sets we have to work with, 
	   and how dense we can tolerate the final matrix.

	   Unlike other filtering implementations, this code attempts to
	   adaptively and automatically determine a good matrix derived from
	   the considerations above, with no user tuning or user parameter
	   choices at all. The critical advantage is that we know at the
	   outset that the complete collection of relation sets will be
	   made to fit in memory, and this allows a global view of the
	   effects of merging choices as the merging runs. It also allows
	   choosing the next ideal to merge in a way that globally minimizes
	   the amount of fill-in added to the current collection of cycles.
	   It thus makes completely unnecessary the complex system of 
	   merge levels, shrinkage passes and relation set cutoff sizes that
	   is typical of filtering implementations derived from Cavallar's
	   paper.

	   We divide I into two sets, the set of active ideals (which we
	   will try to eliminate via merging) and inactive ideals (which
	   we do not try to merge but otherwise process identically). Each
	   large ideal is placed in one of two discrete priority heaps, the
	   active or inactive heap, and each heap is keyed to (an upper 
	   bound on) the amount of fill-in that merging that ideal would 
	   incur. We start with min_cycles ideals in the active heap, and 
	   keep merging ideals until the active heap is empty. A relation 
	   set is considered to be a cycle if it has no active large ideals. 
	   Ideals can move freely between the active and inactive heap,
	   dynamically changing the number of cycles that exist at any given 
	   time. In particular, the active heap will always contain the
	   ideals that will lead to the lightest merges. 

	   When the active heap is empty, if the current number of cycles
	   is at least (min_cycles+(ideals in the inactive heap)), merging
	   can finish. Otherwise, we transfer some ideals from the inactive
	   heap to the active heap, shrinking the number of cycles needed
	   but possibly spoiling previously created cycles, which then have
	   to be merged again. Ideals are also made active if the number
	   of cycles is sufficient but the average cycle weight is too low,
	   since in that case the matrix has room to grow more dense in
	   return for becoming smaller. Finally, whenever the active heap
	   is empty we also 'bury' some of the heaviest ideals in the inactive
	   heap; buried ideals are removed from the ideal lists of all the
	   relation sets containing those ideals, and become ineligible 
	   for merging (since they are not explicitly tracked anymore). 
	   The burying operation leads to a slightly denser matrix, but
	   saves significant amounts of memory as merging progresses.
	   Merging stops when all of these processes converge to a moderate 
	   size matrix with specified density and a sufficient number of 
	   cycles.
	   
	   The effect of this method is that when there are many excess
	   relations, the merge phase produces a usable but too-large matrix
	   as quickly as possible, then gradually reduces its size until
	   the number of nonzero matrix elements exceeds a threshold
	   (or no more merges are possible). */

	uint32 i;
	relation_set_t *relset_array;
	uint32 num_relsets;
	uint32 num_ideals;
	uint32 num_cycles;
	uint32 target_cycles;
	uint32 unmerged_ideals;
	matrix_weight_t mat_weight;
	double avg_cycle_weight;
	heap_t active_heap;
	heap_t inactive_heap;
	ideal_list_t ideal_list;
	merge_aux_t *aux;
	uint64 total_cycle_weight = 0;
	uint32 cycle_bins[NUM_CYCLE_BINS + 2] = {0};
	uint32 max_cycles;
	int32 status = 0;

	logprintf(obj, "commencing full merge\n");

	relset_array = merge->relset_array;
	num_relsets = merge->num_relsets;
	num_ideals = merge->num_ideals;

	/* initialize; all ideals start off inactive */

	aux = (merge_aux_t *)xmalloc(sizeof(merge_aux_t));
	heap_init(&active_heap);
	heap_init(&inactive_heap);
	ideal_list_init(&ideal_list, num_ideals, 0);
	matrix_weight_init(&mat_weight);

	/* add each relation set to the heaps, and count the
	   total relation set weight */

	for (i = num_cycles = 0; i < num_relsets; i++) {
		if (heap_add_relset(&active_heap, &inactive_heap,
				&ideal_list, relset_array + i, i, 0) == 0) {
			num_cycles++;
			matrix_weight_add(&mat_weight, relset_array + i);
		}
	}

	/* transfer the lightest min_cycles ideals to the
	   active heap, modifying the number of cycles and
	   the total cycle weight in the process */

	for (i = 0; i < min_cycles; i++) {
		uint32 ideal = heap_remove_best(&inactive_heap, &ideal_list);
		if (ideal == (uint32)(-1))
			break;

		toggle_ideal_state(&ideal_list, ideal, 
				relset_array, &num_cycles,
				&mat_weight);
		heap_add_ideal(&active_heap, &ideal_list, ideal);
	}
	unmerged_ideals = min_cycles + inactive_heap.num_ideals;
	target_cycles = unmerged_ideals + merge->num_extra_relations;
	avg_cycle_weight = (double)get_matrix_weight(&mat_weight, 
					num_cycles) / num_cycles;

	/* keep forming cycles as long as possible */

	while (1) {

		uint32 ideal;

		if (active_heap.num_ideals == 0) {

			/* we only recalculate the average cycle weight
			   when there are enough cycles to make a valid
			   matrix. Note that we do *not* care about the 
			   average weight of all the cycles we have, but 
			   the average weight of cycles that would 
			   actually appear in the matrix */
	
			if (num_cycles >= target_cycles) {
				avg_cycle_weight = (double)get_matrix_weight(
							&mat_weight, 
							target_cycles) / 
							target_cycles;
				if (avg_cycle_weight >= merge->target_density)
					break;
			}

			/* give up trying to merge a few of the heaviest
			   unmerged ideals. This makes the final matrix
			   somewhat more dense than it has to be, but
			   saves a great deal of memory as the merge
			   phase progresses */

			for (i = 0; i < 500; i++) {

				ideal = heap_remove_worst(&inactive_heap, 
							&ideal_list);
				if (ideal == (uint32)(-1))
					break;

				bury_inactive_ideal(relset_array, 
							&ideal_list, ideal);
				min_cycles++;
			}

			/* we cannot form more cycles, so transfer a 
			   few of the lightest unmerged ideals from 
			   the inactive to the active heap, then 
			   recalculate the target matrix size */

			for (i = 0; i < 2000; i++) {

				ideal = heap_remove_best(&inactive_heap, 
							&ideal_list);
				if (ideal == (uint32)(-1))
					break;

				toggle_ideal_state(&ideal_list, ideal, 
						relset_array, &num_cycles,
						&mat_weight);
				heap_add_ideal(&active_heap, 
						&ideal_list, ideal);
			}

			unmerged_ideals = min_cycles + inactive_heap.num_ideals;
			target_cycles = unmerged_ideals + 
					merge->num_extra_relations;
		}

		/* choose the next ideal to merge */

		ideal = heap_remove_best(&active_heap, &ideal_list);
		if (ideal == (uint32)(-1))
			break;

		/* remove all the relation sets that contain that
		   ideal from both heaps, merge them, and add them
		   back to the heaps, updating the number of cycles
		   formed and the total weight of all cycles */

		load_next_relset_group(aux, &active_heap, &inactive_heap,
					&ideal_list, relset_array, ideal, 0);

		do_merges_core(aux);

		num_cycles += store_next_relset_group(aux, 
					&active_heap, &inactive_heap,
					&ideal_list, relset_array, 
					&mat_weight);

		/* swap ideals between the active and inactive
		   heaps, until all the lightest ideals are in the
		   active heap */

		while (active_heap.num_ideals > 0 &&
		       inactive_heap.num_ideals > 0 &&
		       active_heap.worst_bin > inactive_heap.next_bin) {

			ideal = heap_remove_best(&inactive_heap, &ideal_list);
			toggle_ideal_state(&ideal_list, ideal, 
						relset_array, &num_cycles,
						&mat_weight);
			heap_add_ideal(&active_heap, &ideal_list, ideal);

			ideal = heap_remove_worst(&active_heap, &ideal_list);
			toggle_ideal_state(&ideal_list, ideal, 
						relset_array, &num_cycles,
						&mat_weight);
			heap_add_ideal(&inactive_heap, &ideal_list, ideal);
		}

		/* recalculate the targets to shoot for */

		unmerged_ideals = min_cycles + inactive_heap.num_ideals;
		target_cycles = unmerged_ideals + merge->num_extra_relations;
	}

	/* merging has finished; delete relation sets that
	   are not cycles, then sort in order of increasing
	   cycle weight */

	logprintf(obj, "memory use: %.1f MB\n", (double)
			get_merge_memuse(relset_array, num_relsets,
						&ideal_list) / 1048576);

	for (i = num_cycles = 0; i < num_relsets; i++) {
		relation_set_t *r = relset_array + i;

		if (r->data && r->num_active_ideals == 0) {
			relation_set_t *new_r = relset_array + num_cycles++;
			*new_r = *r;
			new_r->num_small_ideals += new_r->num_large_ideals;
			new_r->num_large_ideals = 0;
			new_r->data = (uint32 *)xrealloc(new_r->data,
							new_r->num_relations *
							sizeof(uint32));
		}
		else {
			free(r->data);
			r->data = NULL;
		}
	}
	logprintf(obj, "found %u cycles, need %u\n", 
				num_cycles, target_cycles);

	if (num_cycles < 0.8 * target_cycles) {
		logprintf(obj, "too few cycles, matrix probably "
				"cannot build\n");
		status = 1;
		goto clean_up;
	}

	qsort(relset_array, (size_t)num_cycles, 
			sizeof(relation_set_t), compare_relsets);

	/* keep only the minimum possible number of 
	   the lightest cycles */

	for (i = target_cycles; i < num_cycles; i++) {
		relation_set_t *r = relset_array + i;
		free(r->data);
		r->data = NULL;
	}
	num_cycles = MIN(target_cycles, num_cycles);
	merge->num_relsets = num_cycles;
	merge->num_ideals = min_cycles + inactive_heap.num_ideals;

	merge->relset_array = relset_array = 
			(relation_set_t *)xrealloc(relset_array, num_cycles * 
						sizeof(relation_set_t));

	/* print statistics on the final collection of cycles */

	for (i = max_cycles = 0; i < num_cycles; i++) {
		relation_set_t *r = relset_array + i;
		total_cycle_weight += r->num_small_ideals;
		max_cycles = MAX(max_cycles, r->num_relations);

		if (r->num_relations > NUM_CYCLE_BINS)
			cycle_bins[NUM_CYCLE_BINS+1]++;
		else
			cycle_bins[r->num_relations]++;
	}

	merge->avg_cycle_weight = (double)total_cycle_weight / num_cycles;
	merge->max_relations = max_cycles;
	logprintf(obj, "weight of %u cycles is about %" PRIu64 
			" (%.2f/cycle)\n",
			num_cycles, total_cycle_weight, 
			merge->avg_cycle_weight);
	logprintf(obj, "distribution of cycle lengths:\n");
	for (i = 1; i < NUM_CYCLE_BINS + 1; i++)
		logprintf(obj, "%u relations: %u\n", i, cycle_bins[i]);
	logprintf(obj, "%u+ relations: %u\n", i, cycle_bins[i]);
	logprintf(obj, "heaviest cycle: %u relations\n", max_cycles);

	/* clean up */
clean_up:
	matrix_weight_free(&mat_weight);
	heap_free(&active_heap);
	heap_free(&inactive_heap);
	ideal_list_free(&ideal_list);
	free(aux);
	return status;
}
