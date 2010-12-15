/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	
       				   --jasonp@boo.net 6/27/07

Modified:	Ben Buhrow
Date:		11/24/09
Purpose:	Port into Yafu-1.14.
--------------------------------------------------------------------*/

#include "qs.h"

static void build_matrix(uint32 ncols, la_col_t *cols, 
		    	siqs_r *relation_list);


/*------------------------------------------------------------------*/
void solve_linear_system(fact_obj_t *obj, uint32 fb_size, 
		    uint64 **bitfield, siqs_r *relation_list, 
		    la_col_t *cycle_list, uint32 *num_cycles) {

	/* Generate linear dependencies among the relations
	   in full_relations and partial_relations */

	la_col_t *cols;
	uint32 nrows, ncols;
	uint64 *dependencies;
	uint32 num_deps;

	ncols = *num_cycles;
	nrows = fb_size;
	cols = cycle_list;

	/* convert the list of relations from the sieving 
	   stage into a matrix. */

	build_matrix(ncols, cols, relation_list);
	count_matrix_nonzero(obj, fb_size, 0, ncols, cols);

	/* reduce the matrix dimensions to ignore almost empty rows */

	reduce_matrix(obj, &nrows, 0, &ncols, cols, NUM_EXTRA_RELATIONS);

	if (ncols == 0) {
		printf("matrix is corrupt; skipping linear algebra\n");
		free(cols);
		*num_cycles = 0;
		return;
	}

	/* solve the linear system */

	dependencies = block_lanczos(obj, nrows, 0, ncols, cols, &num_deps);

	if (num_deps == 0) {
		free(dependencies);
		return;
	}

	*bitfield = dependencies;
	*num_cycles = ncols;
}

/*------------------------------------------------------------------*/
uint32 merge_relations(uint32 *merge_array,
		  uint32 *src1, uint32 n1,
		  uint32 *src2, uint32 n2) {

	/* Given two sorted lists of integers, merge
	   the lists into a single sorted list with
	   duplicate entries removed. If a particular
	   entry occurs an even number of times in the
	   two input lists, don't add it to the final list
	   at all. Returns the number of elements in the
	   resulting list */

	uint32 i1, i2, val1, val2, count1, count2;
	uint32 num_merge;

	i1 = i2 = 0;
	num_merge = 0;

	while (i1 < n1 && i2 < n2) {
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
		} while (i2 < n2 && src2[i2] == val2);

		if (count2 & 1)
			merge_array[num_merge++] = val2;
	}

	return num_merge;
}

/*------------------------------------------------------------------*/
#define MAX_COL_WEIGHT 1000

static void build_matrix(uint32 ncols, la_col_t *cols, 
			   siqs_r *relation_list) {

	/* Convert lists of relations from the sieving stage
	   into a sparse matrix. The matrix is stored by
	   columns, pointed to by 'cols'. The total number 
	   of nonzero entries in the matrix is returned */

	uint32 i, j;
	la_col_t *col;

	/* Cycles are assumed to be sorted in order of increasing
	   number of relations, so that any cycles that are
	   not used would have created the heaviest matrix columns
	   anyway */
	
	for (i = 0; i < ncols; i++) {
		uint32 buf[MAX_COL_WEIGHT];
		uint32 accum[MAX_COL_WEIGHT];
		uint32 weight;

		/* merge each succeeding relation into the accumulated
		   matrix column */

		col = cols + i;

		for (j = weight = 0; j < col->cycle.num_relations; j++) {
			siqs_r *r = &relation_list[col->cycle.list[j]];
			weight = merge_relations(accum, buf, weight,
						r->fb_offsets, r->num_factors);
			memcpy(buf, accum, weight * sizeof(uint32));
		}

		col->weight = weight;
		col->data = (uint32 *)malloc(weight * sizeof(uint32));
		memcpy(col->data, buf, weight * sizeof(uint32));
	}
}
