/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: clique.c 1038 2021-02-07 20:21:24Z jasonp_sf $
--------------------------------------------------------------------*/

#include "filter_priv.h"

	/* Perform the clique removal phase of NFS filtering. This
	   operation does some really cool stuff, and it is not
	   described very thoroughly in references, so I'll do that
	   here.

	   Suppose that at the end of NFS sieving you have X relations.
	   Each relation has zero or more large ideals, and there is
	   a total of Y unique large ideals in the whole dataset.
	   Then this collection of relations can form X-Y cycles, 
	   allowing you to build a matrix for an NFS factor base 
	   of size X-Y. Suppose that the factor base size is Z but
	   X-Y is much larger than that. This means you can get rid
	   of relations until X-Y has been reduced to Z. The opti-
	   mization problem is then: which relations must you throw
	   away so that when X-Y = Z the resulting matrix is as sparse
	   as possible?

	   When you delete a relation, X is decremented; when that
	   relation contained the last instance of some ideal, Y is
	   decremented. If X and Y are both decremented by about the 
	   same amount, X-Y is about the same as before and you can
	   keep removing relations. We assume that when X-Y = Z, the
	   matrix will be more sparse when X is smaller. So the 
	   optimization problem now becomes: which relations must you
	   delete so that the largest number of ideals are also 
	   removed from the dataset? This is a much more specific
	   version of the original optimization problem.

	   Suppose that a particular large ideal occurs in exactly
	   two relations. Then if one of those relations is deleted
	   the other relation becomes a singleton, and will have to be
	   deleted as well. More generally suppose you have a collection
	   of relations, and the collection contains a bunch of ideals
	   that only occur in pairs of those relations. If each relation
	   in the group has at least one of those ideals that only occur 
	   in pairs, then deleting one of the relations means the 
	   whole group can be deleted as well. In graph theory terms, 
	   our collection of relations is a connected subgraph of 
	   the complete dataset. If it was *fully* connected, i.e. every
	   relation had an ideal in common with all other relations in the
	   component, it would have a shorter name: a clique. Even though
	   it is a misnomer to call the component a clique in all cases,
	   with NFS filtering it is traditional to do so. We therefore
	   call the process of finding and deleting connected components
	   a 'clique removal' phase.

	   In this specific case we only want cliques that contain ideals 
	   of 'weight 2', and that problem can be solved efficiently. We 
	   just need a hashtable for the ideals of weight 2, combined with 
	   a breadth-first traversal of the complete dataset. All the 
	   relations that form a clique are in a sense tied to each other,
	   so if relations have to be deleted then it makes sense to 
	   concentrate on deleting the 'heaviest' cliques.

	   A clique is 'heavy' if 1) it has many relations, and 2)
	   those relations contain many ideals that have low weight.
	   These conditions together mean that deleting the entire 
	   clique means reducing both X and Y by nearly the same 
	   (large) amount. 

	   A recent development on the scoring of cliques is due to
	   Cyril Bouvier (http://hal.inria.fr/hal-00734654): if an ideal
	   occurs n times in the dataset and m times in a clique, then
	   the weight of the clique is the sum of m/n for each of the
	   ideals of weight > 2. This actually produces slightly
	   worse results at the end of the clique removal compared to
	   Cavallar's weight function, but produces slightly smaller
	   matrices (2-3% smaller) at the end of the merge phase.
	   
	   The clique processing here locates all of the cliques in a 
	   list of NFS relations, saves a bunch of the heaviest ones, 
	   sorts those so the heaviest cliques are first, and then 
	   deletes cliques until either the value of X-Y is sufficiently 
	   reduced or the entire bunch is gone. The point at which clique 
	   deletion stops is only approximate; removing a clique can create 
	   singletons or other cliques elsewhere in the graph. So we can't 
	   delete too many cliques at once, because we could use up our 
	   budget for reducing X but may end up ignoring heavy cliques 
	   that are created later. However, building the list of cliques 
	   takes time, and we don't want to spend too much time squeezing 
	   out every last heavy clique. The size of the list of cliques 
	   to generate determines how aggressive the clique removal is: 
	   the smaller the size, the more aggressive the removal. Calling
	   code is expected to perform clique removal multiple times in
	   a row, and removing only a few cliques at a time means the number
	   of iterations needed to reduce X must increase.

	   Clique removal is cool because it's equivalent to doing a
	   little Gauss elimination on the complete dataset and deleting
	   the heaviest eliminated super-relations. Only this Gauss
	   elimination is very fast and memory efficient. Clique removal
	   is also cool because it greatly rewards having more relations
	   than you strictly need. The more you're willing to prune, 
	   the better the final matrix will turn out to be. 
	   
	   The above explanation applies to cliques containing ideals
	   that appear in exactly two relations. When there is a really
	   huge amount of initial excess, it is possible for the cliques
	   to run out before the excess does :) To deal with this situation,
	   this code is generalized to find cliques with ideals
	   of weight > 2. This is not ordinarily a good idea, but when 
	   there are no weight-2 ideals and there is still excess to burn,
	   it is the best algorithm available for continuing to prune
	   relations */

/* representation of one clique */

typedef struct {
	uint16 num_relations;       /* number of relations in clique */
	uint16 num_ideals;          /* number of weight-2 ideals in clique */
	float score;                /* measure of how 'heavy' the clique is
				       (higher implies heavier) */
	uint64 *relation_array_word;  /* list of word offsets within the
					 relation array of relations that
					 occur in this clique */
} clique_t;

/*--------------------------------------------------------------------*/

/* boilerplate code for managing a heap of cliques */

#define HEAP_SWAP(a,b) { tmp = a; a = b; b = tmp; }
#define HEAP_PARENT(i)  (((i)-1) >> 1)
#define HEAP_LEFT(i)    (2 * (i) + 1)
#define HEAP_RIGHT(i)   (2 * (i) + 2)

static void heapify(clique_t *h, uint32 index, uint32 size) {

	uint32 c;
	clique_t tmp;
	for (c = HEAP_LEFT(index); c < (size-1); 
			index = c, c = HEAP_LEFT(index)) {

		if( h[c].score > h[c+1].score )
			c++;

		if( h[index].score > h[c].score ) {
			HEAP_SWAP(h[index], h[c]);
		}
		else
			return;
	}
	if (c == (size-1) && h[index].score > h[c].score) {
		HEAP_SWAP(h[index], h[c]);
	}
}

static void make_heap(clique_t *h, uint32 size) {

	uint32 i;
	for (i = HEAP_PARENT(size); i; i--)
		heapify(h, i-1, size);
}

/*--------------------------------------------------------------------*/
static int compare_uint64(const void *x, const void *y) {
	uint64 *xx = (uint64 *)x;
	uint64 *yy = (uint64 *)y;
	if (*xx > *yy)
		return 1;
	if (*xx < *yy)
		return -1;
	return 0;
}

/*--------------------------------------------------------------------*/
static int compare_score_descending(const void *x, const void *y) {
	clique_t *xx = (clique_t *)x;
	clique_t *yy = (clique_t *)y;

	if (xx->score > yy->score)
		return -1;
	if (xx->score < yy->score)
		return 1;
	return 0;
}

/*--------------------------------------------------------------------*/
static void check_relations_array(filter_t *filter, uint32 location) {

	relation_ideal_t *curr_relation;
	uint32 num_relations;
	uint32 num_ideals, ideal;
	uint32 i, j;
	uint32 num_bad = 0;
	uint32 badrel = 0;

	printf("checking relations array at location %u\n", location);
	fflush(stdout);

	curr_relation = filter->relation_array;
	num_relations = filter->num_relations;
	num_ideals = filter->num_ideals;

	for (i = 0; i < num_relations; i++) {
	    	for (j = 0; j < curr_relation->ideal_count; j++) {
			ideal = curr_relation->ideal_list[j];
			if (ideal >= num_ideals) badrel = 1; 
	    	}
	    	if (badrel == 1) {
			printf("Loc %u: bad relation %u of %u, num_ideals is %u\n",
				location, i, num_relations, num_ideals);
			printf("rel_index: %u, ideal_count: %hhu, gf2_factors: %hhu, connected: %hhu\nIdeals: ",
				curr_relation->rel_index, curr_relation->ideal_count, curr_relation->gf2_factors, curr_relation->connected);
			for (j = 0; j < curr_relation->ideal_count; j++)
				printf("%u, ", curr_relation->ideal_list[j]);
			printf("\n");
			badrel = 0;
			num_bad++;
		}
	    	curr_relation = next_relation_ptr(curr_relation);
	}
		
	if (num_bad > 0) printf("Found %u bad relations in array.\n", num_bad);
	fflush(stdout);

	return;
}

/*--------------------------------------------------------------------*/
static void delete_relations(filter_t *filter,
			uint64 *delete_array, uint32 num_delete) {

	uint32 i, j;
	uint32 num_relations = filter->num_relations;
	relation_ideal_t *relation_array = filter->relation_array;
	relation_ideal_t *curr_relation;
	relation_ideal_t *old_relation;

	/* put relation pointers in sorted order */

	qsort(delete_array, (size_t)num_delete, 
			sizeof(uint64), compare_uint64);

	/* walk through the relation list and delete the
	   relations that appear in the delete list */

	curr_relation = relation_array;
	old_relation = relation_array;
	for (i = j = 0; i < num_relations; i++) {
		uint64 array_word = (uint64)((uint32 *)curr_relation -
						(uint32 *)relation_array);
		relation_ideal_t *next_relation = 
				next_relation_ptr(curr_relation);

		if (array_word == delete_array[j]) {

			/* relation is skipped; don't increment
			   the 'delete pointer' if there are no
			   more relations to delete */

			if (j < num_delete - 1)
				j++;
		}
		else {
			/* relation has survived */

			uint8 curr_num_ideals = curr_relation->ideal_count;
			uint32 k;
			old_relation->rel_index = curr_relation->rel_index;
			old_relation->gf2_factors = curr_relation->gf2_factors;
			old_relation->ideal_count = curr_num_ideals;
			for (k = 0; k < curr_num_ideals; k++) {
				old_relation->ideal_list[k] =
						curr_relation->ideal_list[k];
			}
			old_relation = next_relation_ptr(old_relation);
		}
		curr_relation = next_relation;
	}

	/* trim the relation array */

	filter->relation_array = (relation_ideal_t *)xrealloc(relation_array,
				(size_t)(old_relation + 1 - relation_array) *
				sizeof(relation_ideal_t));
	filter->num_relations = num_relations - num_delete;
	check_relations_array(filter, 1);
}

/*--------------------------------------------------------------------*/
static uint32 purge_cliques_core(msieve_obj *obj, 
				filter_t *filter,
				uint32 clique_heap_size,
				uint32 max_clique_relations,
				uint32 num_excess_relations) {

	uint32 i, j;
	ideal_map_t *ideal_map;
	relation_ideal_t *relation_array;
	relation_ideal_t *curr_relation;
	uint32 num_relations;
	uint32 num_ideals;
	uint32 num_ideals_delete;
	clique_t *clique_heap;
	uint32 num_clique;

	uint64 *delete_array;
	uint32 num_delete;
	uint32 num_delete_alloc;

	ideal_relation_t *reverse_array;
	uint32 num_reverse;
	uint32 num_reverse_alloc;

	uint64 *clique_relations;
	uint32 num_clique_relations;
	uint32 num_clique_relations_alloc;

	uint32 *clique_ideals;
	uint32 num_clique_ideals;
	uint32 num_clique_ideals_alloc;

	relation_array = filter->relation_array;
	num_relations = filter->num_relations;
	num_ideals = filter->num_ideals;

	/* set up the hashtable for ideal counts */

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

	/* mark all the ideals with small enough weight as 
	   belonging to a clique, and set the head of their 
	   linked list of relations to empty */

	for (i = 0; i < num_ideals; i++) {
		if (ideal_map[i].payload <= max_clique_relations) {
			ideal_map[i].payload = 0;
			ideal_map[i].clique = 1;
		}
	}

	/* for each relation */

	curr_relation = relation_array;
	for (i = 0; i < num_relations; i++) {

		uint64 relation_array_word = 
				((uint32 *)curr_relation -
				 (uint32 *)relation_array);

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
						relation_array_word;
			reverse_array[num_reverse].next = 
						ideal_map[ideal].payload;
			ideal_map[ideal].payload = num_reverse++;
		}

		curr_relation = next_relation_ptr(curr_relation);
	}

	num_clique = 0;
	clique_heap = (clique_t *)xmalloc(clique_heap_size * sizeof(clique_t));

	num_clique_relations_alloc = 500;
	clique_relations = (uint64 *)xmalloc(num_clique_relations_alloc *
					sizeof(uint64));

	num_clique_ideals_alloc = 500;
	clique_ideals = (uint32 *)xmalloc(num_clique_ideals_alloc *
					sizeof(uint32));

	/* find all the cliques and save the heaviest ones.
	   We perform breadth first search by iterating through
	   all of the ideals */

	for (i = 0; i < num_ideals; i++) {
		clique_t *next_clique;
		uint32 curr_clique_ideal;
		float clique_score = 0.0;

		/* check if the ideal is not part of a clique,
		   or is part of a clique but the clique has
		   already been visited */

		if (!ideal_map[i].clique || ideal_map[i].connected)
			continue;

		/* we've found a clique, and have to measure its
		   weight and find all the relations in it. 

		   This is a fairly standard breadth-first search.
		   We find all the unvisited relations that contain
		   this ideal, list them, and add any unvisited
		   clique ideals they contain to a queue. The clique
		   is complete when the queue of clique ideals is empty */

		num_clique_relations = 0;
		curr_clique_ideal = 0;
		clique_ideals[0] = i;
		ideal_map[i].connected = 1;
		num_clique_ideals = 1;

		/* for each clique ideal in the queue */

		while (curr_clique_ideal < num_clique_ideals) {
			uint32 idx = clique_ideals[curr_clique_ideal];
			uint32 offset = ideal_map[idx].payload;

			/* find the relations containing the ideal */

			while (offset != 0) {
				ideal_relation_t *rev = reverse_array + offset;
				relation_ideal_t *r = (relation_ideal_t *)(
						(uint32 *)relation_array +
						rev->relation_array_word);

				if (r->connected) {
					offset = rev->next;
					continue;
				}

				/* relation seen for the first time;
				   look for clique ideals it contains */

				for (j = 0; j < r->ideal_count; j++) {
					uint32 new_ideal = r->ideal_list[j];
					ideal_map_t *map = ideal_map+new_ideal;

					/* add the contribution of this ideal
					   to the score of the clique */

					if (!map->clique) {
						clique_score += 1.0 / 
							map->payload;
						continue;
					}

					if (map->connected)
						continue;

					/* ideal is a clique ideal and has
					   not been visited before; add it 
					   to the queue */

					if (num_clique_ideals ==
						num_clique_ideals_alloc) {
						num_clique_ideals_alloc *= 2;
						clique_ideals = (uint32 *)
							xrealloc(clique_ideals,
							num_clique_ideals_alloc
							* sizeof(uint32));
					}
					clique_ideals[num_clique_ideals++] =
							new_ideal;
					map->connected = 1;
				}

				/* save the relation and mark as visited */

				if (num_clique_relations ==
						num_clique_relations_alloc) {
					num_clique_relations_alloc *= 2;
					clique_relations = (uint64 *)xrealloc(
						clique_relations,
						num_clique_relations_alloc *
						sizeof(uint64));
				}
				clique_relations[num_clique_relations++] = 
						rev->relation_array_word;
				r->connected = 1;
				offset = rev->next;
			}

			/* current clique ideal is finished */

			curr_clique_ideal++;
		}

		/* clique is enumerated; throw it away if it is
		   too large to fit into a packed structure */

		if (num_clique_relations > 65535 ||
		     num_clique_ideals > 65535)
			continue;

		if (num_clique < clique_heap_size) {
			/* heap not full; append this clique */
			next_clique = clique_heap + num_clique;
		}
		else if (clique_score <= clique_heap[0].score) {
			/* all cliques in the heap are heavier
			   than this one; just skip it */
			continue;
		}
		else {
			/* this clique replaces the lowest-
			   scoring clique in the heap */
			next_clique = clique_heap;
			free(next_clique->relation_array_word);
		}
			
		/* save this clique */

		next_clique->num_relations = (uint16)num_clique_relations;
		next_clique->num_ideals = (uint16)num_clique_ideals;
		next_clique->score = clique_score;
		next_clique->relation_array_word = 
				(uint64 *)xmalloc(num_clique_relations *
						sizeof(uint64));
		memcpy(next_clique->relation_array_word,
			clique_relations, 
			num_clique_relations * sizeof(uint64));

		/* manage the heap */

		if (num_clique < clique_heap_size) {
			if (num_clique == clique_heap_size - 1)
				make_heap(clique_heap, clique_heap_size);
			num_clique++;
		}
		else {
			heapify(clique_heap, 0, clique_heap_size);
		}
	}

	free(reverse_array);
	free(ideal_map);
	free(clique_relations);
	free(clique_ideals);

	/* put the heaviest cliques first */

	qsort(clique_heap, (size_t)num_clique, sizeof(clique_t), 
				compare_score_descending);

	/* now figure out how many cliques to delete, and list all
	   the relations in those cliques */

	num_delete_alloc = 5000;
	num_delete = 0;
	num_ideals_delete = 0;
	delete_array = (uint64 *)xmalloc(num_delete_alloc * sizeof(uint64));

	for (i = 0; i < num_clique; i++) {
		clique_t *curr_clique = clique_heap + i;

		/* we stop deleting cliques either when all of
		   them are gone, or when the number of deleted
		   relations threatens to exceed the number of 
		   deleted ideals by too much.

		   Note that we do not count the number of deleted
		   ideals exactly; we know the clique ideals will
		   be eliminated, but have no way of estimating the
		   number of singletons that deleting the clique will
		   create */

		if (num_delete + curr_clique->num_relations >= 
			num_excess_relations + curr_clique->num_ideals)
			break;

		num_excess_relations += curr_clique->num_ideals;
		num_ideals_delete += curr_clique->num_ideals;

		/* append the list of relations in this clique
		   to the delete list. Since a relation can only
		   appear in one clique, the list will have no
		   redundancy */

		for (j = 0; j < curr_clique->num_relations; j++) {
			if (num_delete == num_delete_alloc) {
				num_delete_alloc *= 2;
				delete_array = (uint64 *)xrealloc(delete_array,
							num_delete_alloc *
							sizeof(uint64));
			}
			delete_array[num_delete++] = 
				curr_clique->relation_array_word[j];
		}
	}

	/* don't bother deleteing cliques if there are too few of
	   them; the rest of the filtering would hardly be affected */

	if (i > 8000) {
		logprintf(obj, "removing %u relations and %u ideals "
				"in %u cliques\n", 
				num_delete, num_ideals_delete, i);
		delete_relations(filter, delete_array, num_delete);
	}
	else {
		num_delete = 0;
	}

	for (i = 0; i < num_clique; i++)
		free(clique_heap[i].relation_array_word);
	free(clique_heap);
	free(delete_array);
	return num_delete;
}

/*--------------------------------------------------------------------*/
void filter_purge_cliques(msieve_obj *obj, filter_t *filter) {

	uint32 clique_heap_size;
	uint32 max_clique_relations = 2;

	/* iteratively delete cliques until there are
	   just enough excess relations for the merge 
	   phase to run. We choose the number of cliques
	   to delete so that there are always several
	   passes to perform, allowing us to discover new
	   cliques to prune */

	clique_heap_size = ((filter->num_relations - 
			filter->num_ideals) - filter->target_excess) / 2;
	if (clique_heap_size > 5000000)
		clique_heap_size = 2000000;
	else if (clique_heap_size > 2000000)
		clique_heap_size = 1000000;
	else
		clique_heap_size = MIN(clique_heap_size, 400000);

	/* iteratively delete relations containing ideals
	   appearing in exactly max_clique_relations relations */

	while (filter->num_relations - filter->num_ideals >
			filter->target_excess + 100) {
		if (purge_cliques_core(obj, filter, clique_heap_size,
				max_clique_relations,
				(filter->num_relations - 
				 filter->num_ideals) - 
				filter->target_excess) == 0) {
			break;
		}

		/* the above got rid of relations; now get rid
		   of the ideals from those relations and remove
		   any additional singletons */

		filter_purge_singletons_core(obj, filter);
	}
}
