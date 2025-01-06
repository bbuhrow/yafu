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

/* below, the routines for building and solving a qs matrix using
jasonp's block lanczos routines are implemented */

static void build_qs_matrix(uint32_t ncols, qs_la_col_t *cols, 
		    	siqs_r *relation_list, uint32_t num_relations);


/*------------------------------------------------------------------*/
int qs_solve_linear_system(fact_obj_t *obj, uint32_t fb_size, 
		    uint64_t **bitfield, siqs_r *relation_list, 
		    qs_la_col_t *cycle_list, uint32_t *num_cycles) {

	/* Generate linear dependencies among the relations
	   in full_relations and partial_relations */

	qs_la_col_t *cols;
	uint32_t nrows, ncols;
	uint64_t *dependencies;
	uint32_t num_deps;

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

	build_qs_matrix(ncols, cols, relation_list, max_rel_index);
	count_qs_matrix_nonzero(obj, fb_size, 0, ncols, cols);

	/* reduce the matrix dimensions to ignore almost empty rows */

	reduce_qs_matrix(obj, &nrows, 0, &ncols, cols, NUM_EXTRA_QS_RELATIONS);

	if (ncols == 0) {
		printf("matrix is corrupt; skipping linear algebra\n");
        cols = NULL;
		*num_cycles = 0;
		return -2;
	}

	/* solve the linear system */

	dependencies = qs_block_lanczos(obj, nrows, 0, ncols, cols, &num_deps);

    if (num_deps == (uint32_t)-1)
    {
        return -2;
    }

	if (num_deps == 0) {
		free(dependencies);
		return -1;
	}

	*bitfield = dependencies;
	*num_cycles = ncols;

    return 0;
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

