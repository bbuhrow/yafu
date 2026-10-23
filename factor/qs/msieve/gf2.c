/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

Some parts of the code (and also this header), included in this 
distribution have been reused from other sources. In particular I 
have benefitted greatly from the work of Jason Papadopoulos's msieve @ 
www.boo.net/~jasonp, Scott Contini's mpqs implementation, and Tom St. 
Denis Tom's Fast Math library.  Many thanks to their kind donation of 
code to the public domain.
       				   --bbuhrow@gmail.com 11/24/09
----------------------------------------------------------------------*/

#include "qs_impl.h"
#include "ytools.h"

/* below, the routines for building a qs matrix and solving it with
the block lanczos solver in common/lanczos (shared with NFS) are
implemented */

static void build_qs_matrix(uint32_t ncols, qs_la_col_t *cols,
		    	siqs_r *relation_list, uint32_t num_relations);

/* entry points into common/lanczos.  These are declared in
ms_include/ms_common.h in terms of la_col_t, which can't be included
alongside yafu's msieve headers.  la_col_t/la_cycle_t have exactly the
same layout as qs_la_col_t/qs_la_cycle_t (uint32 data pointer, uint32
weight, then {uint32 num_relations, uint32 *list}), and both sides
share the msieve_obj definition, so we declare them here directly. */

uint64_t * block_lanczos(msieve_obj *obj,
			uint32_t nrows, uint32_t max_nrows, uint32_t start_row,
			uint32_t num_dense_rows,
			uint32_t ncols, uint32_t max_ncols, uint32_t start_col,
			qs_la_col_t *cols, uint32_t *deps_found);

uint64_t count_matrix_nonzero(msieve_obj *obj,
			uint32_t nrows, uint32_t num_dense_rows,
			uint32_t ncols, qs_la_col_t *cols);

uint64_t reduce_matrix(msieve_obj *obj, uint32_t *nrows,
		uint32_t num_dense_rows, uint32_t *ncols,
		qs_la_col_t *cols, uint32_t num_excess);

/* number of times to re-run lanczos (from a new random start)
   before giving up and asking for more relations */
#define QS_LANCZOS_ATTEMPTS 3

/*------------------------------------------------------------------*/
static msieve_obj * qs_la_obj_new(fact_obj_t *obj) {

	/* make a temporary msieve_obj carrying the settings that
	   common/lanczos reads: threads, cache size, random state,
	   logging, and a name for checkpoint files (only written
	   for matrices with more than 1M rows) */

	msieve_obj *mobj;
	uint32_t flags = 0;
	uint32_t threads = (uint32_t)((obj->LATHREADS > 0) ?
		obj->LATHREADS : obj->THREADS);
	char *logname = NULL;

	if (obj->VFLAG > 0)
		flags |= MSIEVE_FLAG_LOG_TO_STDOUT;

	if ((obj->LOGFLAG) && (obj->flogname[0] != '\0')) {
		flags |= MSIEVE_FLAG_USE_LOGFILE;
		logname = obj->flogname;
		if (obj->logfile != NULL)
			fflush(obj->logfile);
	}

	mobj = msieve_obj_new("", flags, obj->qs_obj.siqs_savefile, logname,
		NULL, obj->seed1, obj->seed2, (uint32_t)0, 9,
		(uint32_t)obj->L1CACHE, (uint32_t)obj->L2CACHE,
		threads, (uint32_t)0, NULL);

	return mobj;
}

/*------------------------------------------------------------------*/
int qs_solve_linear_system(fact_obj_t *obj, uint32_t fb_size,
		    uint64_t **bitfield, siqs_r *relation_list,
		    qs_la_col_t *cycle_list, uint32_t *num_cycles) {

	/* Generate linear dependencies among the relations
	   in full_relations and partial_relations */

	qs_la_col_t *cols;
	uint32_t nrows, ncols;
	uint64_t *dependencies = NULL;
	uint32_t num_deps = 0;
	msieve_obj *mobj;
	int attempt;
	int status;

	ncols = *num_cycles;
	nrows = fb_size;
	cols = cycle_list;

	/* convert the list of relations from the sieving
	   stage into a matrix. */

	uint32_t max_rel_index = 0;
	int i;
	for (i = 0; i < *num_cycles; i++)
	{
		int j;
		for (j = 0; j < cycle_list[i].cycle.num_relations; j++)
		{
			if (cycle_list[i].cycle.list[j] > max_rel_index)
				max_rel_index = cycle_list[i].cycle.list[j];
		}
	}

	mobj = qs_la_obj_new(obj);

	build_qs_matrix(ncols, cols, relation_list, max_rel_index);
	count_matrix_nonzero(mobj, fb_size, 0, ncols, cols);

	/* reduce the matrix dimensions to ignore almost empty rows */

	reduce_matrix(mobj, &nrows, 0, &ncols, cols, NUM_EXTRA_QS_RELATIONS);

	if ((ncols == 0) || (ncols <= nrows)) {
		printf("matrix is corrupt (%u x %u); skipping linear algebra\n",
			nrows, ncols);
		*num_cycles = 0;
		msieve_obj_free(mobj);
		return -2;
	}

	/* solve the linear system.  block_lanczos frees the column
	   data (but not the cycle lists) when it finishes, so a retry
	   rebuilds the columns that survived reduce_matrix.  The random
	   state in mobj advances with each attempt. */

	for (attempt = 0; attempt < QS_LANCZOS_ATTEMPTS; attempt++) {

		if (attempt > 0) {
			uint32_t new_nrows = fb_size;

			printf("lanczos failed, retrying (attempt %d of %d)\n",
				attempt + 1, QS_LANCZOS_ATTEMPTS);
			build_qs_matrix(ncols, cols, relation_list, max_rel_index);
			reduce_matrix(mobj, &new_nrows, 0, &ncols, cols,
				NUM_EXTRA_QS_RELATIONS);
			nrows = new_nrows;
		}

		dependencies = block_lanczos(mobj, nrows, nrows, 0, 0,
			ncols, ncols, 0, cols, &num_deps);

		if ((dependencies != NULL) && (num_deps > 0))
			break;

		free(dependencies);
		dependencies = NULL;
		num_deps = 0;
	}

	/* keep yafu's random state moving forward */
	obj->seed1 = mobj->seed1;
	obj->seed2 = mobj->seed2;
	msieve_obj_free(mobj);

	if (num_deps == 0) {
		status = -1;
	}
	else {
		*bitfield = dependencies;
		*num_cycles = ncols;
		status = 0;
	}

    return status;
}

/*------------------------------------------------------------------*/
uint32_t qs_merge_relations(uint32_t *merge_array,
		  uint32_t *src1, uint32_t n1,
		  uint32_t *src2, uint32_t n2) {

	/* Given two sorted lists of integers, merge
	   the lists into a single sorted list with
	   duplicate entries removed. If a particular
	   entry occurs an even number of times in the
	   two input lists, don't add it to the final list
	   at all. Returns the number of elements in the
	   resulting list */

	uint32_t i1, i2, val1, val2, count1, count2;
	uint32_t num_merge;

	i1 = i2 = 0;
	num_merge = 0;

	while ((i1 < n1) && (i2 < n2)) {
		val1 = src1[i1];
		val2 = src2[i2];

		if (val1 < val2) {
			count1 = 0;
			do {
				i1++; count1++;
			} while (i1 < n1 && src1[i1] == val1);

			if (count1 & 1)
				merge_array[num_merge++] = val1;
		}
		else if (val1 > val2) {
			count2 = 0;
			do {
				i2++; count2++;
			} while (i2 < n2 && src2[i2] == val2);

			if (count2 & 1)
				merge_array[num_merge++] = val2;
		}
		else {
			count1 = count2 = 0;
			do {
				i1++; count1++;
			} while (i1 < n1 && src1[i1] == val1);
			do {
				i2++; count2++;
			} while (i2 < n2 && src2[i2] == val2);

			if ( (count1 + count2) & 1 )
				merge_array[num_merge++] = val1;
		}
	}

	if (i2 == n2) {
		src2 = src1;
		i2 = i1;
		n2 = n1;
	}

	while (i2 < n2) {
		count2 = 0; val2 = src2[i2];

		do {
			i2++; count2++;
		} while ((i2 < n2) && (src2[i2] == val2));

		if (count2 & 1)
			merge_array[num_merge++] = val2;
	}

	return num_merge;
}

/*------------------------------------------------------------------*/
//#define QS_MAX_COL_WEIGHT 10000

void build_qs_matrix(uint32_t ncols, qs_la_col_t *cols, 
			   siqs_r *relation_list, uint32_t num_relations) {

	/* Convert lists of relations from the sieving stage
	   into a sparse matrix. The matrix is stored by
	   columns, pointed to by 'cols'. The total number 
	   of nonzero entries in the matrix is returned */

	uint32_t i, j;
	qs_la_col_t *col;
    uint32_t* buf;
    uint32_t* accum;
    int QS_MAX_COL_WEIGHT = 10000;

	/* Cycles are assumed to be sorted in order of increasing
	   number of relations, so that any cycles that are
	   not used would have created the heaviest matrix columns
	   anyway */
	
	// printf("building matrix with %u columns\n", ncols);

    buf = (uint32_t*)xmalloc(QS_MAX_COL_WEIGHT * sizeof(uint32_t));
    accum = (uint32_t*)xmalloc(QS_MAX_COL_WEIGHT * sizeof(uint32_t));
	for (i = 0; i < ncols; i++) {
		
		uint32_t weight;

		/* merge each succeeding relation into the accumulated
		   matrix column */

		col = cols + i;

		for (j = weight = 0; j < col->cycle.num_relations; j++) {
			siqs_r *r = &relation_list[col->cycle.list[j]];
			
			//if (r->num_factors > 40)
			//	printf("warning merging large relation with %d factors "
			//		"with index %d\n", r->num_factors, col->cycle.list[j]);

			if ((weight + r->num_factors) > QS_MAX_COL_WEIGHT)
			{
                QS_MAX_COL_WEIGHT *= 2;
				printf("weight = %u: allocated new max weight vectors of size %u\n", 
					(weight + r->num_factors), QS_MAX_COL_WEIGHT); fflush(stdout);
                buf = (uint32_t*)xrealloc(buf, QS_MAX_COL_WEIGHT * sizeof(uint32_t));
                accum = (uint32_t*)xrealloc(accum, QS_MAX_COL_WEIGHT * sizeof(uint32_t));
			}

			weight = qs_merge_relations(accum, buf, weight,
						r->fb_offsets, r->num_factors);
			memcpy(buf, accum, weight * sizeof(uint32_t));
		}

		//printf("column %d has weight %u after merging %u relations\n",
		//	i, weight, col->cycle.num_relations);

		col->weight = weight;
		col->data = (uint32_t *)xmalloc(weight * sizeof(uint32_t));
		memcpy(col->data, buf, weight * sizeof(uint32_t));
	}

    free(buf);
    free(accum);
    return;
}

