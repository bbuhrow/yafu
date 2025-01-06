/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: filter_priv.h 1038 2021-02-07 20:21:24Z jasonp_sf $
--------------------------------------------------------------------*/

/* implementation of relation filtering, using an intermediate
   representation for relations that is amenable to either QS or NFS.
   These routines work on an input filter_t structure and produce
   an output merge_t structure containing the completed relation sets */

#ifndef _COMMON_FILTER_FILTER_PRIV_H_
#define _COMMON_FILTER_FILTER_PRIV_H_

#include "filter.h"

#ifdef __cplusplus
extern "C" {
#endif

/* structure for the mapping between large ideals 
   and relations (used during clique removal) */

typedef struct {
	uint64 relation_array_word;  /* 64-bit word offset into relation
					array where the relation starts */
	uint32 next;		     /* next relation containing this ideal */
} ideal_relation_t;

/* structure used to map between a large ideal and a 
   linked list of relations that use that ideal */

typedef struct {
	uint64 payload : 32;	/* offset in list of ideal_relation_t
				   structures where the linked list of
				   ideal_relation_t's for this ideal starts */
	uint64 clique : 1;      /* nonzero if this ideal can participate in
				   a clique */
	uint64 connected : 1;   /* nonzero if this ideal has already been
				   added to a clique under construction */
} ideal_map_t;

/* a relation_set_t simulates matrix rows; the following simulates
   the matrix columns, mapping ideals to the relation sets
   containing those ideals */

typedef struct ideal_set_t {
	uint16 num_relsets;     /* the number of relation sets
				   containing this ideal */
	uint16 num_relsets_alloc; /* maximum number of relset numbers the
				     'relsets' array can hold */
	uint16 active;          /* 1 if ideal is active, 0 if inactive */
	uint16 min_relset_size; /* the number of ideals in the
				   lightest relation set that
				   contains this ideal */
	uint32 *relsets;        /* list of members in an array of 
				   relation sets that contain this
				   ideal (no ordering assumed) */
	struct ideal_set_t *next;
	struct ideal_set_t *prev;  /* used to build circular linked lists */
} ideal_set_t;

/* relation sets with more than this many relations are deleted */

#define MAX_RELSET_SIZE 28

/* perform clique removal on the current set of relations */

void filter_purge_cliques(msieve_obj *obj, filter_t *filter);

/* initialize the merge process */

void filter_merge_init(msieve_obj *obj, filter_t *filter); 

/* perform all 2-way merges, converting the results into
   relation-sets that the main merge routine operates on */

void filter_merge_2way(msieve_obj *obj, filter_t *filter, merge_t *merge);

/* do the rest of the merging. min_cycles is the minimum number
   of cycles that the input collection of relation-sets must
   produce, corresponding to the smallest matrix that can be
   built (the actual matrix is expected to be much larger than 
   this). */

int32 filter_merge_full(msieve_obj *obj, merge_t *merge, uint32 min_cycles);

#ifdef __cplusplus
}
#endif

#endif /* _COMMON_FILTER_FILTER_PRIV_H_ */
