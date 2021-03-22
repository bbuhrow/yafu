/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	
       				   --jasonp@boo.net 9/24/08

Modified:	Ben Buhrow
Date:		11/24/09
Purpose:	Port into Yafu-1.14.
--------------------------------------------------------------------*/

#include "lanczos.h"
#include "ytools.h"
#include "qs_impl.h"
#include "qs.h"

typedef struct {
	uint32_t count;
	uint32_t index;
} qs_row_count_t;

static int yafu_compare_row_count(const void *x, const void *y) {
	qs_row_count_t *xx = (qs_row_count_t *)x;
	qs_row_count_t *yy = (qs_row_count_t *)y;
	return yy->count - xx->count;
}

static int yafu_compare_row_index(const void *x, const void *y) {
	qs_row_count_t *xx = (qs_row_count_t *)x;
	qs_row_count_t *yy = (qs_row_count_t *)y;
	return xx->index - yy->index;
}

static int yafu_compare_uint32_t(const void *x, const void *y) {
	uint32_t *xx = (uint32_t *)x;
	uint32_t *yy = (uint32_t *)y;
	if (*xx > *yy)
		return 1;
	if (*xx < *yy)
		return -1;
	return 0;
}

static int yafu_compare_weight(const void *x, const void *y) {
	qs_la_col_t *xx = (qs_la_col_t *)x;
	qs_la_col_t *yy = (qs_la_col_t *)y;
	return xx->weight - yy->weight;
}

/*------------------------------------------------------------------*/
void count_qs_matrix_nonzero(fact_obj_t *obj,
			uint32_t nrows, uint32_t num_dense_rows,
			uint32_t ncols, qs_la_col_t *cols) {

	uint32_t i, j;
	uint32_t total_weight;
	uint32_t sparse_weight;
	size_t mem_use;
	
	mem_use = ncols * (sizeof(qs_la_col_t) +
		sizeof(uint32_t) * ((num_dense_rows + 31) / 32));

	for (i = total_weight = sparse_weight = 0; i < ncols; i++) {
		uint32_t w = cols[i].weight;
		total_weight += w;
		sparse_weight += w;
		mem_use += w * sizeof(uint32_t);
	}

	if (num_dense_rows > 0) {
		for (i = 0; i < ncols; i++) {
			uint32_t *dense_rows = cols[i].data + cols[i].weight;
			for (j = 0; j < num_dense_rows; j++) {
				if (dense_rows[j / 32] & (1 << (j % 32)))
					total_weight++;
			}
		}
	}

	if (obj->VFLAG > 0)
		printf("matrix is %u x %u (%.1f MB) with "
			"weight %u (%5.2f/col)\n", 
				nrows, ncols, 
				(double)mem_use / 1048576,
				total_weight, 
				(double)total_weight / ncols);
	if (obj->VFLAG > 0)
		printf("sparse part has weight %u (%5.2f/col)\n", 
				sparse_weight, 
				(double)sparse_weight / ncols);

	if (obj->logfile != NULL)
		logprint(obj->logfile, "matrix is %u x %u (%.1f MB) with "
			"weight %u (%5.2f/col)\n", 
				nrows, ncols, 
				(double)mem_use / 1048576,
				total_weight, 
				(double)total_weight / ncols);

	if (obj->logfile != NULL)
		logprint(obj->logfile, "sparse part has weight %u (%5.2f/col)\n", 
				sparse_weight, 
				(double)sparse_weight / ncols);
}

/*------------------------------------------------------------------*/
#define QS_MAX_COL_WEIGHT 4000

static void yafu_combine_cliques(uint32_t num_dense_rows, 
			uint32_t *ncols_out, qs_la_col_t *cols, 
			qs_row_count_t *counts) {

	uint32_t i, j;
	uint32_t ncols = *ncols_out;
	uint32_t dense_row_words = (num_dense_rows + 31) / 32;

	uint32_t num_merged;
	uint32_t merge_array[QS_MAX_COL_WEIGHT];

	/* for each row, mark the last column encountered 
	   that contains a nonzero entry in that row */

	for (i = 0; i < ncols; i++) {
		qs_la_col_t *c = cols + i;
		for (j = 0; j < c->weight; j++) {
			counts[c->data[j]].index = i;
		}
	}

	/* for each column */

	for (i = 0; i < ncols; i++) {
		qs_la_col_t *c0;
		qs_la_col_t *c1 = cols + i;
		uint32_t clique_base = (uint32_t)(-1);

		if (c1->data == NULL)
			continue;

		/* if the column hits to a row of weight 2, and the
		   other column containing this row is distinct and
		   not previously merged */

		for (j = 0; j < c1->weight; j++) {
			qs_row_count_t *curr_clique = counts + c1->data[j];
			if (curr_clique->count == 2) {
				clique_base = curr_clique->index;
				break;
			}
		}

		if (clique_base == (uint32_t)(-1) || clique_base == i)
			continue;

		c0 = cols + clique_base;
		if (c0->data == NULL || 
		    c0->weight + c1->weight >= QS_MAX_COL_WEIGHT)
			continue;

		/* remove c0 and c1 from the row counts */

		for (j = 0; j < c0->weight; j++)
			counts[c0->data[j]].count--;
		for (j = 0; j < c1->weight; j++)
			counts[c1->data[j]].count--;

		/* merge column c1 into column c0. First merge the
		   nonzero entries (do not assume they are sorted) */
		qsort(c0->data, (size_t)c0->weight, 
					sizeof(uint32_t), yafu_compare_uint32_t);
		qsort(c1->data, (size_t)c1->weight, 
					sizeof(uint32_t), yafu_compare_uint32_t);
		num_merged = qs_merge_relations(merge_array, 
						c0->data, c0->weight,
						c1->data, c1->weight);
		for (j = 0; j < dense_row_words; j++) {
			merge_array[num_merged + j] = c0->data[c0->weight+j] ^
						      c1->data[c1->weight+j];
		}
		free(c0->data);
		c0->data = (uint32_t *)xmalloc((num_merged + 
					dense_row_words) * sizeof(uint32_t));
		memcpy(c0->data, merge_array, (num_merged + 
					dense_row_words) * sizeof(uint32_t));
		c0->weight = num_merged;

		/* then combine the two lists of relation numbers */

		c0->cycle.list = (uint32_t *)xrealloc(c0->cycle.list, 
					(c0->cycle.num_relations +
					 c1->cycle.num_relations) *
					sizeof(uint32_t));
		memcpy(c0->cycle.list + c0->cycle.num_relations,
			c1->cycle.list, c1->cycle.num_relations * 
					sizeof(uint32_t));
		c0->cycle.num_relations += c1->cycle.num_relations;

		/* add c0 back into the row counts */

		for (j = 0; j < c0->weight; j++) {
			qs_row_count_t *curr_row = counts + c0->data[j];
			curr_row->count++;
			curr_row->index = clique_base;
		}

		/* kill off c1 */

		free(c1->data);
		c1->data = NULL;
		free(c1->cycle.list);
		c1->cycle.list = NULL;
	}

	/* squeeze out the merged columns from the list */

	for (i = j = 0; i < ncols; i++) {
		if (cols[i].data != NULL)
			cols[j++] = cols[i];
	}
	*ncols_out = j;
}

/*------------------------------------------------------------------*/
void reduce_qs_matrix(fact_obj_t *obj, uint32_t *nrows, 
		uint32_t num_dense_rows, uint32_t *ncols, 
		qs_la_col_t *cols, uint32_t num_excess) {

	/* Perform light filtering on the nrows x ncols
	   matrix specified by cols[]. The processing here is
	   limited to collapsing cliques, deleting columns that 
	   contain a singleton row, deleting empty rows, and then 
	   deleting the heaviest columns until the matrix has a 
	   few more columns than rows. Because deleting a column 
	   reduces the counts in several different rows, the process 
	   must iterate to convergence.
	   
	   Note that deleting singleton rows is not intended to 
	   make the Lanczos iteration run any faster (though it will); 
	   it's just that if we don't go to this trouble and the matrix
	   has many zero rows, then Lanczos iteration could fail 
	   to find any nontrivial dependencies. I've also seen cases
	   where cliques *must* be merged in order to find nontrivial
	   dependencies; this seems to happen for matrices that are large
	   and very sparse */

	uint32_t r, c, i, j, k;
	uint32_t passes;
	qs_row_count_t *counts;
	uint32_t reduced_rows;
	uint32_t reduced_cols;

	/* sort the columns in order of increasing weight */
	qsort(cols, (size_t)(*ncols), sizeof(qs_la_col_t), yafu_compare_weight);

	/* count the number of nonzero entries in each row */

	reduced_rows = *nrows;
	reduced_cols = *ncols;
	passes = 0;

	counts = (qs_row_count_t *)xcalloc((size_t)reduced_rows, 
					sizeof(qs_row_count_t));
	for (i = 0; i < reduced_cols; i++) {
		for (j = 0; j < cols[i].weight; j++)
			counts[cols[i].data[j]].count++;
	}

	do {
		r = reduced_rows;

		/* remove any columns that contain the only entry
		   in one or more rows, then update the row counts
		   to reflect the missing column. Iterate until
		   no more columns can be deleted */

		do {
			c = reduced_cols;

			/* delete columns that contain a singleton row */

			for (i = j = 0; i < reduced_cols; i++) {
				qs_la_col_t *col = cols + i;
				for (k = 0; k < col->weight; k++) {
					if (counts[col->data[k]].count < 2)
						break;
				}
	
				if (k < col->weight) {
					for (k = 0; k < col->weight; k++) {
						counts[col->data[k]].count--;
					}
					free(col->data);
					free(col->cycle.list);
				}
				else {
					cols[j++] = cols[i];
				}
			}
			reduced_cols = j;

			/* if the matrix is big enough, collapse most 
			   of the cliques that it contains */

			if (reduced_cols >= QS_MIN_NCOLS_TO_PACK) {
				yafu_combine_cliques(num_dense_rows, 
						&reduced_cols, 
						cols, counts);
			}
		} while (c != reduced_cols);
	
		/* count the number of rows that contain a
		   nonzero entry. Ignore the row indices associated
		   with the dense rows */

		for (i = reduced_rows = num_dense_rows; i < *nrows; i++) {
			if (counts[i].count)
				reduced_rows++;
		}

		/* Because deleting a column reduces the weight
		   of many rows, the number of nonzero rows may
		   be much less than the number of columns. Delete
		   more columns until the matrix has the correct
		   aspect ratio. Columns at the end of cols[] are
		   the heaviest, so delete those (and update the 
		   row counts again) */

		if (reduced_cols > reduced_rows + num_excess) {
			for (i = reduced_rows + num_excess;
					i < reduced_cols; i++) {

				qs_la_col_t *col = cols + i;
				for (j = 0; j < col->weight; j++) {
					counts[col->data[j]].count--;
				}
				free(col->data);
				free(col->cycle.list);
			}
			reduced_cols = reduced_rows + num_excess;
		}

		/* if any columns were deleted in the previous step,
		   then the matrix is less dense and more columns
		   can be deleted; iterate until no further deletions
		   are possible */

		passes++;

	} while (r != reduced_rows);

	/* if the linear system was underdetermined going
	   into this routine, the pruning above will likely
	   have destroyed the matrix. Linear algebra clearly
	   cannot proceed in this case */

	if (reduced_cols == 0) {
		free(counts);
		*nrows = reduced_rows;
		*ncols = reduced_cols;
		return;
	}

	if (obj->VFLAG > 0)
		printf("filtering completed in %u passes\n", passes);

	if (obj->logfile != NULL)
		logprint(obj->logfile, "filtering completed in %u passes\n", passes);
	count_qs_matrix_nonzero(obj, reduced_rows, num_dense_rows,
				reduced_cols, cols);

	/* permute the row indices to remove rows with zero
	   weight, put the heaviest row indices together,
	   and put each column in sorted order. The first
	   num_dense_rows rows are not affected */
	
	for (i = num_dense_rows; i < *nrows; i++)
		counts[i].index = i;
	qsort(counts + num_dense_rows, (size_t)(*nrows - num_dense_rows), 
			sizeof(qs_row_count_t), yafu_compare_row_count);
	for (i = num_dense_rows; i < *nrows; i++)
		counts[i].count = i;
	qsort(counts + num_dense_rows, (size_t)(*nrows - num_dense_rows), 
			sizeof(qs_row_count_t), yafu_compare_row_index);

	for (i = 0; i < reduced_cols; i++) {
		qs_la_col_t *col = cols + i;
		for (j = 0; j < col->weight; j++) {
			col->data[j] = counts[col->data[j]].count;
		}
		qsort(col->data, (size_t)col->weight, 
				sizeof(uint32_t), yafu_compare_uint32_t);
	}

	/* make heavy columns alternate with light columns; this
	   smooths out the distribution of nonzero entries of the
	   matrix. If we wanted to solve the matrix in parallel,
	   the best approach to achieve load balancing across
	   multiple CPUs is to use some sort of graph partitioning
	   scheme, but these consume huge amounts of memory */

	for (i = 1, j = reduced_cols - 2; i < j; i += 2, j -= 2) {
		qs_la_col_t tmp = cols[i];
		cols[i] = cols[j];
		cols[j] = tmp;
	}

	/* record the final matrix size */

	free(counts);
	*nrows = reduced_rows;
	*ncols = reduced_cols;
}
