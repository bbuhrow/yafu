/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: merge_util.h 23 2009-07-20 02:59:07Z jasonp_sf $
--------------------------------------------------------------------*/

#ifndef _COMMON_FILTER_MERGE_UTIL_H_
#define _COMMON_FILTER_MERGE_UTIL_H_

#include "filter.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MERGE_MAX_OBJECTS 500

/* structure for merging relations that all have an ideal
   in common */

typedef struct {
	uint32 num_relsets;            /* number of relation sets to merge */
	relation_set_t tmp_relsets[MERGE_MAX_OBJECTS];   /* array to hold 
							    relation sets */
	relation_set_t tmp_relsets2[MERGE_MAX_OBJECTS];  /* scratch array for 
							    tmp_relsets */
	uint32 tmp_relset_idx[MERGE_MAX_OBJECTS]; /* offsets into main 
						     relation set array 
						     where these relsets 
						     occurred (they're put 
						     back into these 
						     locations) */
	uint32 tmp_relations[MERGE_MAX_OBJECTS]; /* scratch buffer for 
						    merging relation numbers 
						    from 2 relsets */
	uint32 tmp_ideals[MERGE_MAX_OBJECTS]; /* scratch array for merging 
						 ideals from two relsets */
} merge_aux_t;

/* structure for tracking ideals during merging. Each ideal
   has one entry in each of the two arrays in this struct */

typedef struct {
	uint32 num_ideals;
	ideal_set_t *list; /* for circular linked lists of ideals 
			      connected together by a priority queue */
} ideal_list_t;

/* set up an ideal list. is_active gives the starting
   state of each ideal in the list (0 or 1) */

void ideal_list_init(ideal_list_t *list, uint32 num_ideals, uint32 is_active);

/* clean up an ideal list */

void ideal_list_free(ideal_list_t *list);

/* a discrete double-ended heap for ideals */

typedef struct {
	ideal_set_t *hashtable;  /* hashtable for locating linked lists of
			       ideals; the hash function is the key value
			       and there are num_bins possible entries */
	uint32 num_bins;    /* number of hash bins */
	uint32 next_bin;    /* index into hashtable of the best entry in
			       the priority queue (=num_bins if empty) */
	uint32 worst_bin;   /* index into hashtable of the worst entry in
			       the priority queue (= -1 if empty) */
	uint32 num_ideals;  /* number of ideals in the priority queue */
} heap_t;

void heap_init(heap_t *heap);
void heap_free(heap_t *heap);

/* add an ideal to the heap */

void heap_add_ideal(heap_t *heap, ideal_list_t *list, uint32 ideal);

/* remove an ideal from the heap */

void heap_remove_ideal(heap_t *heap, ideal_list_t *list, uint32 ideal);

/* add all the ideals in a relation set r, that appear in more 
   than min_ideal_weight relation sets, to active_heap or
   inactive_heap, depending on the ideal */

uint32 heap_add_relset(heap_t *active_heap, 
			heap_t *inactive_heap,
			ideal_list_t *list, 
			relation_set_t *r,
			uint32 relset_num,
			uint32 min_ideal_weight);
	
/* remove the ideals in a relation set r from active_heap and
   inactive_heap. If an ideal appears in min_ideal_weight or
   less relations after r is removed, then remove it from
   active_heap or inactive_heap entirely */

void heap_remove_relset(heap_t *active_heap, 
			heap_t *inactive_heap,
			ideal_list_t *list,
			relation_set_t *r,
			relation_set_t *relset_array,
			uint32 min_ideal_weight);

/* remove-min and remove-max operations for the heap */

uint32 heap_remove_best(heap_t *heap, ideal_list_t *list);
uint32 heap_remove_worst(heap_t *heap, ideal_list_t *list);

/* remove all the relation sets containing 'ideal' from 
   active_heap and inactive_heap, storing them in aux. The
   removal process also handles all the other ideals in
   the relation sets */

void load_next_relset_group(merge_aux_t *aux,
			heap_t *active_heap, heap_t *inactive_heap, 
			ideal_list_t *ideal_list,
			relation_set_t *relset_array,
			uint32 ideal, uint32 min_ideal_weight);

/* find all the relation sets that contain a given ideal,
   and remove the ideal from those relation sets. Each
   relation set has its count of small ideals incremented;
   basically this forces the ideal to never be merged */

void bury_inactive_ideal(relation_set_t *relset_array,
			ideal_list_t *ideal_list, uint32 ideal);

size_t get_merge_memuse(relation_set_t *relsets, uint32 num_relsets,
			ideal_list_t *ideal_list);

/* the largest number of relation sets to be merged
   via the spanning tree algorithm */

#define SPANNING_TREE_MAX_RELSETS 20

/* functions for merging relation sets */

void merge_two_relsets(relation_set_t *r1, relation_set_t *r2, 
			relation_set_t *r_out, merge_aux_t *aux);

uint32 estimate_new_weight(relation_set_t *r1, relation_set_t *r2);

#ifdef __cplusplus
}
#endif

#endif /* _COMMON_MERGE_UTIL_H_ */
