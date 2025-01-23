/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: merge_pre.c 23 2009-07-20 02:59:07Z jasonp_sf $
--------------------------------------------------------------------*/

#include "filter_priv.h"

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
#define NUM_IDEAL_BINS 6

void filter_merge_init(msieve_obj *obj, filter_t *filter) {

	/* start the merge process. Right now this only prints
	   histograms of the number of ideals per relation, and
	   sorts the ideals of each relation into ascending order */

	uint32 i;
	relation_ideal_t *curr_relation = filter->relation_array;
	uint32 num_relations = filter->num_relations;
	uint32 ideal_count_bins[NUM_IDEAL_BINS+2] = {0};

	for (i = 0; i < num_relations; i++) {
		uint32 num_ideals = curr_relation->ideal_count;
		if (num_ideals > NUM_IDEAL_BINS)
			ideal_count_bins[NUM_IDEAL_BINS+1]++;
		else
			ideal_count_bins[num_ideals]++;

		if (num_ideals > 1) {
			qsort(curr_relation->ideal_list, (size_t)num_ideals,
					sizeof(uint32), compare_uint32);
		}
		curr_relation = next_relation_ptr(curr_relation);
	}

	for (i = 0; i < NUM_IDEAL_BINS+1; i++) {
		logprintf(obj, "relations with %u large ideals: %u\n", 
					i, ideal_count_bins[i]);
	}
	logprintf(obj, "relations with %u+ large ideals: %u\n", 
				i, ideal_count_bins[i]);
}

/*--------------------------------------------------------------------*/
#define MAX_2WAY_RELATIONS 200
#define MAX_2WAY_IDEALS 1000

void filter_merge_2way(msieve_obj *obj, filter_t *filter,
			merge_t *merge) {

	/* performs 2-way merges and initializes the structures
	   needed by the rest of the merge code. 2-way merges are
	   handled separately because there are a lot of them to do
	   (typically 30% of the ideals participate in 2-way merges) 
	   and they can be performed very efficiently.

	   We can perform a 2-way merge by locating all the 2-way 
	   cliques in the input dataset, then collapsing each clique 
	   into a relation set. A relation belongs to at most one 
	   clique, and no optimization is possible in the combining 
	   process */

	uint32 i, j;
	ideal_map_t *ideal_map;
	relation_ideal_t *relation_array;
	relation_ideal_t *curr_relation;
	uint32 num_relations;
	uint32 num_ideals;
	uint32 num_deleted;

	ideal_relation_t *reverse_array;
	uint32 num_reverse;
	uint32 num_reverse_alloc;

	relation_set_t *relset_array;
	uint32 num_relset;
	uint32 num_relset_alloc;

	logprintf(obj, "commencing 2-way merge\n");

	/* set up the hashtable for ideal counts */

	relation_array = filter->relation_array;
	num_relations = filter->num_relations;
	num_ideals = filter->num_ideals;
	ideal_map = (ideal_map_t *)xcalloc((size_t)num_ideals, 
					sizeof(ideal_map_t));

	/* set up structure for linked lists of clique relations */

	num_reverse = 1;
	num_reverse_alloc = 10000;
	reverse_array = (ideal_relation_t *)xmalloc(num_reverse_alloc *
					sizeof(ideal_relation_t));

	/* count the number of times each ideal occurs in relations */

	curr_relation = relation_array;
	for (i = 0; i < num_relations; i++) {
		curr_relation->connected = 0;
		for (j = 0; j < curr_relation->ideal_count; j++) {
			uint32 ideal = curr_relation->ideal_list[j];
			ideal_map[ideal].payload++;
		}
		curr_relation = next_relation_ptr(curr_relation);
	}

	/* mark all the ideals with weight 2 as belonging
	   to a clique, and set the head of their linked
	   list of relations to empty */

	for (i = 0; i < num_ideals; i++) {
		if (ideal_map[i].payload == 2) {
			ideal_map[i].payload = 0;
			ideal_map[i].clique = 1;
		}
	}

	/* for each relation */

	curr_relation = relation_array;
	for (i = 0; i < num_relations; i++) {

		/* for each ideal in the relation */

		for (j = 0; j < curr_relation->ideal_count; j++) {
			uint32 ideal = curr_relation->ideal_list[j];

			if (!ideal_map[ideal].clique)
				continue;

			/* relation belongs in a clique because of this
			   ideal; add it to the ideal's linked list */

			if (num_reverse == num_reverse_alloc) {
				num_reverse_alloc *= 2;
				reverse_array = (ideal_relation_t *)xrealloc(
						reverse_array,
						num_reverse_alloc *
						sizeof(ideal_relation_t));
			}
			reverse_array[num_reverse].relation_array_word =
					(uint32 *)curr_relation -
					(uint32 *)relation_array;
			reverse_array[num_reverse].next = 
						ideal_map[ideal].payload;
			ideal_map[ideal].payload = num_reverse++;
		}

		curr_relation = next_relation_ptr(curr_relation);
	}

	num_relset = 0;
	num_relset_alloc = 10000;
	relset_array = (relation_set_t *)xmalloc(num_relset_alloc *
					sizeof(relation_set_t));
	num_deleted = 0;

	/* find all the cliques and convert to relation sets.
	   We perform breadth first search by iterating through
	   all of the relations.
	  
	   Note that the clique removal step iterated through
	   ideals looking for cliques. It could do that because
	   the only objective was to find cliques. Here the objective
	   is different: find all relations, performing extra work
	   on cliques. If we discovered relations through the ideals
	   they contained, we'd have to do a lot more traversing 
	   of linked lists */

	curr_relation = relation_array;
	num_deleted = 0;
	for (i = 0; i < num_relations; 
		i++, curr_relation = next_relation_ptr(curr_relation)) {

		uint32 tmp_relations[MAX_2WAY_RELATIONS];
		uint32 tmp_ideals[MAX_2WAY_IDEALS];
		uint32 accum_ideals[MAX_2WAY_IDEALS];
		uint32 num_tmp_relation = 0;
		uint32 num_tmp_ideal = 0;
		uint32 num_small_ideal;
		relation_set_t *curr_relset;

		if (curr_relation->connected)
			continue;

		/* relation hasn't been seen before. Start counting
		   the ideals in it, and list the large ideals that 
		   the relation contains. The clique is complete when
		   all the ideals in the list have been processed */

		tmp_relations[num_tmp_relation++] = curr_relation->rel_index;
		num_small_ideal = curr_relation->gf2_factors;
		num_tmp_ideal = curr_relation->ideal_count;
		for (j = 0; j < num_tmp_ideal; j++)
			tmp_ideals[j] = curr_relation->ideal_list[j];
		curr_relation->connected = 1;

		/* for each ideal in the combined ideal list */

		for (j = 0; j < num_tmp_ideal; j++) {
			uint32 ideal = tmp_ideals[j];
			uint32 offset;

			/* check if the ideal is not part of a clique */

			if (!ideal_map[ideal].clique)
				continue;

			/* we've found a clique, and have to find all the 
			   relations in it */

			offset = ideal_map[ideal].payload;
			while (offset) {

				ideal_relation_t *rev = reverse_array + offset;
				relation_ideal_t *r = (relation_ideal_t *)
						((uint32 *)relation_array +
						rev->relation_array_word);
				if (r->connected) {
					offset = rev->next;
					continue;
				}

				/* relation seen for the first time;
				   add its count of ideals to the totals for
				   the current clique, and merge its ideal
				   list with that of the current clique. The 
				   relation is now considered processed */

				/* check for overflow */

				if (num_tmp_ideal + r->ideal_count >=
							MAX_2WAY_IDEALS) {
					printf("error: clique merge requires "
						"too many ideals\n");
					exit(-1);
				}

				if (num_tmp_relation == MAX_2WAY_RELATIONS) {
					printf("error: clique merge requires "
						"too many relations\n");
					exit(-1);
				}

				/* perform the merge */

				tmp_relations[num_tmp_relation++] = 
							r->rel_index;
				num_small_ideal += r->gf2_factors;
				num_tmp_ideal = merge_relations(accum_ideals,
						tmp_ideals, num_tmp_ideal,
						r->ideal_list, r->ideal_count);
				memcpy(tmp_ideals, accum_ideals, 
						num_tmp_ideal * sizeof(uint32));
				r->connected = 1;
				offset = rev->next;
			}

			/* the ideals were merged into the existing list,
			   so that new ideals to examine can be anywhere in
			   tmp_ideals. Thus we have to start the loop for
			   examining ideals all over again. We have to
			   restart the loop once for every relation in the
			   clique, but this isn't a big deal because most
			   cliques have only two relations */

			j = (uint32)(-1);
		}

		/* clique is complete; throw it away if it 
		   contains too many relations */
		
		if (num_tmp_relation >= MAX_RELSET_SIZE) {
			num_deleted++;
			continue;
		}
		
		/* allocate a relation set for the clique */

		if (num_relset == num_relset_alloc) {
			num_relset_alloc *= 2;
			relset_array = (relation_set_t *)xrealloc(relset_array,
							num_relset_alloc *
							sizeof(relation_set_t));
		}
		curr_relset = relset_array + num_relset;
		curr_relset->num_relations = num_tmp_relation;
		curr_relset->num_small_ideals = num_small_ideal;
		curr_relset->num_large_ideals = num_tmp_ideal;

		/* sort the relation numbers in ascending order
		   (the ideal list is already ordered), then store */

		curr_relset->data = (uint32 *) xmalloc(sizeof(uint32) *
					(num_tmp_relation + num_tmp_ideal));
		if (num_tmp_relation > 1) {
			qsort(tmp_relations, (size_t)num_tmp_relation,
					sizeof(uint32), compare_uint32);
		}
		memcpy(curr_relset->data, tmp_relations,
				num_tmp_relation * sizeof(uint32));
		memcpy(curr_relset->data + num_tmp_relation, tmp_ideals,
				num_tmp_ideal * sizeof(uint32));
		num_relset++;
	}

	/* free unneeded objects */

	free(reverse_array);
	free(filter->relation_array);
	filter->relation_array = NULL;
	relset_array = (relation_set_t *)xrealloc(relset_array,
				num_relset * sizeof(relation_set_t));

	/* renumber the ideals to skip ideals that have been
	   completely merged */

	memset(ideal_map, 0, num_ideals * sizeof(ideal_map_t));
	for (i = 0; i < num_relset; i++) {
		relation_set_t *r = relset_array + i;
		uint32 *ideal_list = r->data + r->num_relations;
		for (j = 0; j < r->num_large_ideals; j++) {
			ideal_map[ideal_list[j]].payload++;
		}
	}
	for (i = j = 0; i < num_ideals; i++) {
		if (ideal_map[i].payload) {
			ideal_map[i].payload = j++;
		}
	}
	num_ideals = j;
	for (i = 0; i < num_relset; i++) {
		relation_set_t *r = relset_array + i;
		uint32 *ideal_list = r->data + r->num_relations;
		for (j = 0; j < r->num_large_ideals; j++) {
			ideal_list[j] = ideal_map[ideal_list[j]].payload;
		}
	}

	logprintf(obj, "reduce to %u relation sets and %u "
			"unique ideals\n", num_relset, num_ideals);
	if (num_deleted) {
		logprintf(obj, "ignored %u oversize "
				"relation sets\n", num_deleted);
	}
	merge->relset_array = relset_array;
	merge->num_relsets = num_relset;
	merge->num_ideals = num_ideals;
	free(ideal_map);
}
