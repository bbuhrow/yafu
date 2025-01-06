/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: filter.h 1020 2018-05-18 13:33:21Z jasonp_sf $
--------------------------------------------------------------------*/

/* implementation of relation filtering, using an intermediate
   representation for relations that is amenable to either QS or NFS.
   These routines work on an input filter_t structure and produce
   an output merge_t structure containing the completed relation sets */

#ifndef _COMMON_FILTER_FILTER_H_
#define _COMMON_FILTER_FILTER_H_

#include <ms_common.h>

#ifdef __cplusplus
extern "C" {
#endif

/* the singleton removal phase uses a packed representation
   of each relation. The following maps relations to large ideals;
   the ideals themselves don't matter, only the unique number
   that each is assigned */

typedef struct {
	uint32 rel_index;      /* savefile line number where relation occurs */
	uint8 ideal_count;     /* number of large ideals */
	uint8 gf2_factors;     /* number of small factors that are ignored */
	uint8 connected;       /* scratch space used in clique formation */
	uint32 ideal_list[TEMP_FACTOR_LIST_SIZE];  /* the large ideals */
} relation_ideal_t;

/* relations can have between 0 and TEMP_FACTOR_LIST_SIZE 
   large ideals, and the average number is very small (2 or 3). 
   When processing large numbers of relations we do not 
   assume a full relation_ideal_t struct is allocated for each. 
   Instead all such structures are packed together as
   tightly as possible, and we always iterate sequentially
   through the array. The following abstracts away the process
   of figuring out the number of bytes used in a current
   relation_ideal_t in the list, then pointing beyond it */

static INLINE relation_ideal_t *next_relation_ptr(relation_ideal_t *r) {
	return (relation_ideal_t *)((uint8 *)r +
			sizeof(relation_ideal_t) - 
			sizeof(r->ideal_list[0]) * 
			(size_t)(TEMP_FACTOR_LIST_SIZE - r->ideal_count));
}

/* main structure controlling relation filtering */

typedef struct {
	relation_ideal_t *relation_array;  /* relations after singleton phase */
	uint32 num_relations;       /* current number of relations */
	uint32 num_ideals;          /* current number of unique large ideals */
	uint32 filtmin_r;           /* min. value a rational ideal needs 
				       to be tracked during filtering */
	uint32 filtmin_a;           /* min. value an algebraic ideal needs 
				       to be tracked during filtering */
	uint32 target_excess;      /* how many more relations than ideals
					are required for filtering to proceed */
	uint32 max_ideal_weight;   /* the largest number of relations that
				        contain the same ideal */
	uint64 lp_file_size;       /* number of bytes in LP file */
} filter_t;

/* representation of a 'relation set', i.e. a group
   of relations that will all participate in the same
   column when the final matrix is formed */

typedef struct {
	uint16 num_relations;     /* number of relations in this relation set */
	uint16 num_small_ideals;  /* number of ideals that are not tracked */
	uint16 num_large_ideals;  /* number of ideals that are tracked */
	uint16 num_active_ideals; /* number of ideals eligible for merging
				     (0 means relset is a cycle) */
	uint32 *data;             /* the first num_relations elements are
				     the relation numbers that participate in
				     this relation set, while the last 
				     num_large_ideals elements give the ideals
				     that are tracked for merging. Both lists
				     are assumed sorted in ascending order */
} relation_set_t;

/* structure controlling the merge process */

typedef struct {
	uint32 num_relsets;           /* current number of relation sets */
	uint32 num_ideals;            /* number of unique ideals that must
					 be merged */
	uint32 num_extra_relations;   /* number of excess relsets required */
	double target_density;        /* target number of nonzeros per cycle */
	double avg_cycle_weight;      /* the avg number of nonzeros per cycle */
	uint32 max_relations;         /* largest relations in a relation set */
	relation_set_t *relset_array; /* current list of relation sets */
} merge_t;

/* the large-prime-only versions of relations are assumed to 
   start off in a disk file. The following reads the relations
   into filter->relation_array, skipping over ideals that occur
   more and max_ideal_weight times. max_ideal_weight = 0 means
   read the entire file into memory */

void filter_read_lp_file(msieve_obj *obj, filter_t *filter,
				uint32 max_ideal_weight);

/* perform approximate singleton removal on relations packed
   into a disk file, then prune the singletons from the file.
   The relations are not read into memory */

void filter_purge_lp_singletons(msieve_obj *obj, filter_t *filter,
				uint64 ram_size);

#if 0
/* perform clique removal on relations packed into a disk file, 
   then prune the singletons from the file. The relations are 
   not read into memory. This code anticipates that later passes
   will need excess relations of their own, and so is a bit
   conservative about how much excess is left when clique removal
   completes */

void filter_purge_lp_cliques(msieve_obj *obj, filter_t *filter);
#endif

/* remove singletons from a memory-resident collection of relations */

void filter_purge_singletons_core(msieve_obj *obj, filter_t *filter);

/* create the matrix columns corresponding to a collection 
   of relations. min_cycles is the absolute minimum number of
   relation sets that are required */

int32 filter_make_relsets(msieve_obj *obj, filter_t *filter, 
			merge_t *merge, uint32 min_cycles);

/* perform post-processing optimizations on the collection of cycles
   found by the merge phase */

void filter_postproc_relsets(msieve_obj *obj, merge_t *merge);

void filter_free_relsets(merge_t *merge);

void filter_dump_relsets(msieve_obj *obj, merge_t *merge);

#ifdef __cplusplus
}
#endif

#endif /* _COMMON_FILTER_FILTER_H_ */
