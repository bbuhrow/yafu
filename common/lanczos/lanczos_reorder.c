/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: lanczos_reorder.c 348 2010-07-25 02:51:09Z jasonp_sf $
--------------------------------------------------------------------*/

#include "lanczos.h"

typedef struct {
	int16 cost;
	uint32 num_edges : 12;
	uint32 visited : 1;
	uint32 boundary : 1;
	uint32 partition : 1;
	uint32 is_col : 1;
	uint32 orig_index;
	uint32 edge_offset;

	uint32 next;
	uint32 prev;
} vertex_t;

typedef struct {
	uint32 num_rows;
	uint32 num_cols;
	uint32 num_rows_orig;
	uint32 num_cols_orig;

	vertex_t *vertex;
	uint32 *edges;

	uint32 stack_alloc;
	uint32 stack_num;
	uint32 *stack;
} graph_t;

typedef struct {
	int32 num_bins;
	int32 next_bin;
	int32 bias;

	vertex_t *hashtable;

	uint32 partition_max[2];
	uint32 partition_size[2][2];
} heap_t;

#define MIN_ROW_IDX 50000

#define INVALID_INDEX ((uint32)(-1))

/*--------------------------------------------------------------------*/
#define MAX_ROW_ENTRIES 1000
#define MAX_EDGES 40

static void graph_init(msieve_obj *obj, graph_t *graph, heap_t *heap) {

	uint32 i, j, k;
	uint32 nrows;
	uint32 dense_row_words;
	uint32 ncols;
	char buf[256];
	FILE *matrix_fp;

	uint32 num_edges;
	uint32 num_edges_alloc;
	uint32 *edges;
	vertex_t *vertex;
	uint32 num_vertex;

	sprintf(buf, "%s.mat", obj->savefile.name);
	matrix_fp = fopen(buf, "rb");
	if (matrix_fp == NULL) {
		logprintf(obj, "error: can't open matrix file\n");
		exit(-1);
	}

	fread(&nrows, sizeof(uint32), (size_t)1, matrix_fp);
	fread(&dense_row_words, sizeof(uint32), (size_t)1, matrix_fp);
	fread(&ncols, sizeof(uint32), (size_t)1, matrix_fp);
	dense_row_words = (dense_row_words + 31) / 32;

	vertex = (vertex_t *)xcalloc((size_t)ncols + nrows - MIN_ROW_IDX, 
					sizeof(vertex_t));

	num_vertex = 0;
	num_edges = 0;
	num_edges_alloc = 100000;
	edges = (uint32 *)xmalloc(num_edges_alloc * sizeof(uint32));

	for (i = 0; i < ncols; i++) {

		uint32 num_row_idx;
		uint32 row_entries[MAX_ROW_ENTRIES];
		vertex_t *v;

		fread(&num_row_idx, sizeof(uint32), 
				(size_t)1, matrix_fp);
		fread(row_entries, sizeof(uint32), 
				(size_t)(num_row_idx + dense_row_words), 
				matrix_fp);

		if (num_edges + num_row_idx >= num_edges_alloc) {
			num_edges_alloc *= 1.4;
			edges = (uint32 *)xrealloc(edges, num_edges_alloc *
						sizeof(uint32));
		}

		for (j = num_row_idx, k = 0; j; j--) {
			uint32 row = row_entries[j-1];

			if (row < MIN_ROW_IDX)
				break;

			v = vertex + ncols + row - MIN_ROW_IDX;
			if (v->num_edges < MAX_EDGES) {
				v->num_edges++;
				edges[num_edges + k] = row - MIN_ROW_IDX;
				k++;
			}
		}
		if (k == 0 || (num_row_idx - k) > MAX_EDGES)
			continue;

		v = vertex + num_vertex++;
		v->num_edges = k;
		v->is_col = 1;
		v->orig_index = i;
		v->edge_offset = num_edges;
		num_edges += k;
	}
	fclose(matrix_fp);

	for (i = num_edges = 0; i < num_vertex; i++) {
		vertex_t *v = vertex + i;
		uint32 curr_num_edges = v->num_edges;
		uint32 *row_list = edges + v->edge_offset;

		for (j = k = 0; j < curr_num_edges; j++) {
			uint32 row = row_list[j];
			vertex_t *row_v = vertex + ncols + row;

			if (row_v->num_edges < MAX_EDGES) {
				edges[num_edges + k] = row;
				k++;
			}
		}
		v->num_edges = k;
		v->edge_offset = num_edges;
		num_edges += k;
	}

	for (i = 0; i < nrows - MIN_ROW_IDX; i++) {
		vertex[ncols + i].num_edges = 0;
	}

	for (i = j = 0; i < num_vertex; i++) {
		vertex_t *v = vertex + i;
		if (v->num_edges > 0) {
			uint32 curr_num_edges = v->num_edges;
			uint32 *row_list = edges + v->edge_offset;

			for (k = 0; k < curr_num_edges; k++) {
				vertex_t *row_v = vertex + ncols + row_list[k];
				row_v->num_edges++;
			}
			vertex[j++] = *v;
		}
	}
	num_vertex = j;

	for (i = 0, j = num_edges, k = num_vertex; 
				i < nrows - MIN_ROW_IDX; i++) {
		vertex_t *v = vertex + ncols + i;
		if (v->num_edges > 0 && v->num_edges < MAX_EDGES) {
			v->orig_index = i + MIN_ROW_IDX;
			v->is_col = 0;
			v->edge_offset = j;
			v->next = k++;
			j += v->num_edges;
		}
	}

	for (i = 0; i < num_vertex; i++) {
		vertex_t *v = vertex + i;
		uint32 curr_num_edges = v->num_edges;
		uint32 *row_list = edges + v->edge_offset;

		for (j = 0; j < curr_num_edges; j++) {
			vertex_t *row_v = vertex + ncols + row_list[j];
			row_list[j] = row_v->next;
		}
	}

	for (i = j = 0; i < nrows - MIN_ROW_IDX; i++) {
		vertex_t *v = vertex + ncols + i;
		if (v->num_edges > 0 && v->num_edges < MAX_EDGES) {
			vertex[num_vertex + j] = *v;
			j++;
		}
	}

	graph->num_cols_orig = ncols;
	graph->num_rows_orig = nrows;
	logprintf(obj, "matrix is %u x %u\n", nrows, ncols);
	graph->num_cols = ncols = num_vertex;
	graph->num_rows = nrows = j;
	logprintf(obj, "sparse core is %u x %u\n", nrows, ncols);
	logprintf(obj, "graph has %u edges\n", num_edges);
	logprintf(obj, "memory use: %.2lf MB\n", 
			((nrows + ncols) * sizeof(vertex_t) +
			 2 * num_edges * sizeof(uint32)) / 1048576.0);

	edges = (uint32 *)xrealloc(edges, 2 * num_edges * sizeof(uint32));

	for (i = 0; i < ncols; i++) {
		vertex_t *v = vertex + i;
		uint32 curr_num_edges = v->num_edges;
		uint32 *row_list = edges + v->edge_offset;

		for (j = 0; j < curr_num_edges; j++) {
			vertex_t *row_v = vertex + row_list[j];
			uint32 *row_edges = edges + row_v->edge_offset;
			row_edges[row_v->prev++] = i;
		}
	}

	for (i = j = 0; i < nrows + ncols; i++) {
		vertex_t *v = vertex + i;
		j = MAX(j, v->num_edges);
	}
	logprintf(obj, "max edge weight is %u\n", j);
	graph->edges = edges;
	graph->vertex = (vertex_t *)xrealloc(vertex, 
					(nrows + ncols + 2 * j) *
					sizeof(vertex_t));
	heap->num_bins = 2 * j;
	heap->bias = j;
	heap->hashtable = graph->vertex + (nrows + ncols);

	graph->stack_alloc = 100;
	graph->stack = (uint32 *)xmalloc(graph->stack_alloc * sizeof(uint32));
}

/*--------------------------------------------------------------------*/
static void graph_free(graph_t *graph, uint32 **rowperm_out,
			uint32 **colperm_out) {

	uint32 i, j;
	uint32 nrows = graph->num_rows;
	uint32 ncols = graph->num_cols;
	uint32 nrows_orig = graph->num_rows_orig;
	uint32 ncols_orig = graph->num_cols_orig;
	uint32 rowgap = nrows_orig - nrows;
	uint32 colgap = ncols_orig - ncols;
	uint32 *rowperm;
	uint32 *colperm;

	free(graph->edges);
	rowperm = (uint32 *)xmalloc(nrows_orig * sizeof(uint32));
	colperm = (uint32 *)xmalloc(ncols_orig * sizeof(uint32));

	for (i = 0; i < nrows_orig; i++)
		rowperm[i] = INVALID_INDEX;

	for (i = 0; i < nrows; i++) {
		vertex_t *v = graph->vertex + ncols + i;
		rowperm[v->orig_index] = i + rowgap;
	}

	for (i = j = 0; i < nrows_orig; i++) {
		if (rowperm[i] == INVALID_INDEX)
			rowperm[i] = j++;
	}


	for (i = 0; i < ncols_orig; i++)
		colperm[i] = INVALID_INDEX;

	for (i = 0; i < ncols; i++) {
		vertex_t *v = graph->vertex + i;
		colperm[v->orig_index] = i + colgap;
	}

	for (i = j = 0; i < ncols_orig; i++) {
		if (colperm[i] == INVALID_INDEX)
			colperm[i] = j++;
	}

	free(graph->vertex);
	free(graph->stack);
	*rowperm_out = rowperm;
	*colperm_out = colperm;
}

/*--------------------------------------------------------------------*/
static void heap_insert_vertex(graph_t *graph, 
			heap_t *heap, uint32 v_offset) {

	vertex_t *base = graph->vertex;
	vertex_t *v = base + v_offset;
	int32 key = heap->bias + v->cost;
	vertex_t *head = heap->hashtable + key;
	uint32 head_offset = head - base;

	v->prev = head_offset;
	v->next = head->next;
	base[head->next].prev = v_offset;
	head->next = v_offset;

	heap->next_bin = MIN(heap->next_bin, key);
}

/*--------------------------------------------------------------------*/
static void heap_delete_vertex(graph_t *graph,
			heap_t *heap, uint32 v_offset) {

	vertex_t *base = graph->vertex;
	vertex_t *v = base + v_offset;
	int32 key = heap->bias + v->cost;

	base[v->next].prev = v->prev;
	base[v->prev].next = v->next;
	v->next = v_offset;
	v->prev = v_offset;

	if (key == heap->next_bin) {
		vertex_t *head = heap->hashtable + key;
		uint32 head_offset = head - base;
		while (key < heap->num_bins && head->next == head_offset) {
			head_offset++;
			head++;
			key++;
		}
		heap->next_bin = key;
	}
}

/*--------------------------------------------------------------------*/
static uint32 heap_delete_min(graph_t *graph, heap_t *heap) {

	vertex_t *base = graph->vertex;
	int32 key = heap->next_bin;

	if (key == heap->num_bins)
		return INVALID_INDEX;

	while (key < heap->num_bins) {
		vertex_t *head = heap->hashtable + key;
		uint32 head_offset = head - base;
		uint32 tmp_offset = head->next;

		while (tmp_offset != head_offset) {
			vertex_t *v = base + tmp_offset;
			uint32 partition = v->partition;

			if (heap->partition_size[v->is_col][partition ^ 1] <
					heap->partition_max[v->is_col]) {
				heap_delete_vertex(graph, heap, tmp_offset);
				return tmp_offset;
			}

			tmp_offset = v->next;
		}
		key++;
	}

	return INVALID_INDEX;
}

/*--------------------------------------------------------------------*/
static uint32 heap_init(msieve_obj *obj, 
			graph_t *graph, heap_t *heap, 
			uint32 row_start, uint32 row_end,
			uint32 col_start, uint32 col_end) {

	uint32 i, j, k;
	uint32 edge_cut = 0;
	
	heap->next_bin = heap->num_bins;
	for (i = 0; i < (uint32)heap->num_bins; i++) {
		vertex_t *head = heap->hashtable + i;
		head->next = head - graph->vertex;
		head->prev = head - graph->vertex;
	}

	heap->partition_max[0] = 1 + 1.01 * ((row_end - row_start + 1) / 2);
	heap->partition_max[1] = 1 + 1.01 * ((col_end - col_start + 1) / 2);
	heap->partition_size[0][0] = 0;
	heap->partition_size[0][1] = 0;
	heap->partition_size[1][0] = 0;
	heap->partition_size[1][1] = 0;

	for (i = 0; i < 2; i++) {
		uint32 vertex_lower = (i == 0) ? row_start : col_start;
		uint32 vertex_upper = (i == 0) ? row_end : col_end;

		j = get_rand(&obj->seed1, &obj->seed2);
		for (k = vertex_lower; k <= vertex_upper; k++) {
			vertex_t *v = graph->vertex + k;
			uint32 partition = j & 1;
	
			v->cost = 0;
			v->boundary = 0;
			v->visited = 1;
			v->partition = partition;
			heap->partition_size[v->is_col][partition]++;

			j >>= 1;
			if (k % 32 == 0)
				j = get_rand(&obj->seed1, &obj->seed2);
		}
	}

	for (i = 0; i < 2; i++) {
		uint32 vertex_lower1 = (i == 0) ? row_start : col_start;
		uint32 vertex_upper1 = (i == 0) ? row_end : col_end;
		uint32 vertex_lower2 = (i == 0) ? col_start : row_start;
		uint32 vertex_upper2 = (i == 0) ? col_end : row_end;

		for (j = vertex_lower1; j <= vertex_upper1; j++) {
			vertex_t *v = graph->vertex + j;
			uint32 num_edges = v->num_edges;
			uint32 *edges = graph->edges + v->edge_offset;
	
			for (k = 0; k < num_edges; k++) {
				uint32 v_offset2 = edges[k];
				vertex_t *v2 = graph->vertex + v_offset2;

				if (v_offset2 < vertex_lower2 ||
				    v_offset2 > vertex_upper2)
					continue;

				if (v->partition == v2->partition) {
					v->cost++;
				}
				else {
					v->cost--;
					edge_cut++;
				}
			}
			v->next = j;
			v->prev = j;
		}
	}

	return edge_cut;
}

/*--------------------------------------------------------------------*/
static uint32 repartition_vertex(graph_t *graph, heap_t *heap,
			uint32 v_offset, uint32 edge_cut, 
			uint32 row_start, uint32 row_end,
			uint32 col_start, uint32 col_end) {

	uint32 i;
	vertex_t *v = graph->vertex + v_offset;
	uint32 num_edges = v->num_edges;
	uint32 *edges = graph->edges + v->edge_offset;
	uint32 partition = v->partition;

	v->cost = -(v->cost);
	v->visited = 1;
	v->partition = partition ^ 1;

	heap->partition_size[v->is_col][partition]--;
	heap->partition_size[v->is_col][partition ^ 1]++;

	for (i = 0; i < num_edges; i++) {
		uint32 v2_offset = edges[i];
		vertex_t *v2 = graph->vertex + v2_offset;

		if (v2->is_col) {
			if (v2_offset < col_start || v2_offset > col_end)
				continue;
		}
		else {
			if (v2_offset < row_start || v2_offset > row_end)
				continue;
		}

		if (!v2->visited)
			heap_delete_vertex(graph, heap, v2_offset);

		if (v2->partition == partition) {
			v2->cost -= 2;
			edge_cut++;
		}
		else {
			v2->cost += 2;
			edge_cut--;
		}

		if (!v2->visited)
			heap_insert_vertex(graph, heap, v2_offset);
	}

	return edge_cut;
}

/*--------------------------------------------------------------------*/
#define MAX_HISTORY 100

static uint32 do_partition_core(graph_t *graph, 
				heap_t *heap, uint32 edge_cut,
				uint32 row_start, uint32 row_end,
				uint32 col_start, uint32 col_end) {

	uint32 i, j;
	uint32 min_edge_cut;
	uint32 num_history;
	uint32 history[MAX_HISTORY];

	for (i = 0; i < 2; i++) {
		uint32 vertex_lower = (i == 0) ? row_start : col_start;
		uint32 vertex_upper = (i == 0) ? row_end : col_end;

		for (j = vertex_lower; j <= vertex_upper; j++) {
			vertex_t *v = graph->vertex + j;
			if (v->visited) {
				v->visited = 0;
				heap_insert_vertex(graph, heap, j);
			}
		}
	}

	min_edge_cut = edge_cut;
	num_history = 0;

	while (1) {
		uint32 v_offset = heap_delete_min(graph, heap);
		if (v_offset == INVALID_INDEX)
			break;

		edge_cut = repartition_vertex(graph, heap, v_offset, edge_cut,
						row_start, row_end,
						col_start, col_end);

		if (edge_cut <= min_edge_cut) {
			num_history = 0;
			min_edge_cut = edge_cut;
		}
		else {
			history[num_history] = v_offset;
			if (++num_history == MAX_HISTORY)
				break;
		}
	}

	for (i = num_history - 1; (int32)i >= 0; i--) {
		edge_cut = repartition_vertex(graph, heap, 
					history[i], edge_cut,
					row_start, row_end,
					col_start, col_end);
	}

	return min_edge_cut;
}

/*--------------------------------------------------------------------*/
static void permute_graph_core(graph_t *graph,
			uint32 v_lower1, uint32 v_upper1,
			uint32 v_lower2, uint32 v_upper2,
			uint32 *middle1, uint32 *middle2) {

	int32 i, j, k;
	vertex_t *vertex_base = graph->vertex;
	uint32 *edge_base = graph->edges;
	int32 v_lower, v_upper;
	vertex_t tmp;

	for (i = v_lower1; i <= v_upper1; i++) {
		vertex_t *v = vertex_base + i;
		int32 num_edges = v->num_edges;
		uint32 *edges = edge_base + v->edge_offset;

		v->prev = i;
		for (j = k = 0; j < num_edges; j++) {
			int32 v2_offset = edges[j];
			if (v2_offset >= v_lower2 && 
			    v2_offset <= v_upper2) {
				k++;
			}
		}
		if (v->cost != k)
			v->boundary = 1;
	}

	v_lower = v_lower1;
	v_upper = v_upper1;
	i = v_lower - 1;
	j = v_upper + 1;
	while (1) {
		while (++i <= v_upper) {
			if (vertex_base[i].boundary == 0)
				break;
		}
		while (--j >= v_lower) {
			if (vertex_base[j].boundary == 1)
				break;
		}
		if (j < i)
			break;

		tmp = vertex_base[i];
		vertex_base[i] = vertex_base[j];
		vertex_base[j] = tmp;
	}

	*middle1 = v_lower = MIN(v_upper, i);
	i = v_lower - 1;
	j = v_upper + 1;
	while (1) {
		while (++i <= v_upper) {
			if (vertex_base[i].partition == 1)
				break;
		}
		while (--j >= v_lower) {
			if (vertex_base[j].partition == 0)
				break;
		}
		if (j < i)
			break;

		tmp = vertex_base[i];
		vertex_base[i] = vertex_base[j];
		vertex_base[j] = tmp;
	}
	*middle2 = MIN(v_upper, i);

	for (i = v_lower1; i <= v_upper1; i++) {
		vertex_t *v0 = vertex_base + i;
		vertex_t *v1 = vertex_base + v0->prev;
		v1->next = i;
	}

	for (i = v_lower2; i <= v_upper2; i++) {
		vertex_t *v = vertex_base + i;
		int32 num_edges = v->num_edges;
		uint32 *edges = edge_base + v->edge_offset;

		for (j = 0; j < num_edges; j++) {
			int32 v1_offset = edges[j];
			if (v1_offset >= v_lower1 && 
			    v1_offset <= v_upper1) {
				vertex_t *v0 = vertex_base + v1_offset;
				edges[j] = v0->next;
			}
		}
	}
}

/*--------------------------------------------------------------------*/
static void permute_graph(graph_t *graph,
			uint32 row_start, uint32 row_end,
			uint32 col_start, uint32 col_end,
			uint32 *row_start1, uint32 *row_start2,
			uint32 *col_start1, uint32 *col_start2) {

	permute_graph_core(graph, 
			row_start, row_end,
			col_start, col_end,
			row_start1, row_start2);

	permute_graph_core(graph, 
			col_start, col_end,
			row_start, row_end,
			col_start1, col_start2);
}

/*--------------------------------------------------------------------*/
#define MIN_BLOCKSIZE 4000

static void do_partition(msieve_obj *obj,
			graph_t *graph, heap_t *heap,
			uint32 row_start, uint32 row_end,
			uint32 col_start, uint32 col_end) {

	uint32 min_edge_cut, edge_cut;
	uint32 row_start1, row_start2;
	uint32 col_start1, col_start2;
#ifdef VERBOSE
	printf("optimize [%u,%u] x [%u,%u]\n",
			row_start, row_end, col_start, col_end);
#endif
	edge_cut = heap_init(obj, graph, heap,
		       		row_start, row_end,
				col_start, col_end);
	do {
		min_edge_cut = edge_cut;
#ifdef VERBOSE
		printf("edge cut = %u\n", edge_cut);
#endif
		edge_cut = do_partition_core(graph, heap, edge_cut,
						row_start, row_end,
						col_start, col_end);
	} while (min_edge_cut - edge_cut > 1000);

#ifdef VERBOSE
	printf("edge cut = %u\n", edge_cut);
#endif

	permute_graph(graph,
			row_start, row_end,
			col_start, col_end,
			&row_start1, &row_start2,
			&col_start1, &col_start2);

	if (col_end - col_start2 > MIN_BLOCKSIZE) {
		do_partition(obj, graph, heap,
				row_start2, row_end,
				col_start2, col_end);
	}
	if (col_start2 - col_start1 > MIN_BLOCKSIZE) {
		do_partition(obj, graph, heap,
				row_start1, row_start2 - 1,
				col_start1, col_start2 - 1);
	}
	if (col_start1 - col_start > MIN_BLOCKSIZE) {
		do_partition(obj, graph, heap,
				row_start, row_start1 - 1,
				col_start, col_start1 - 1);
	}
}

/*--------------------------------------------------------------------*/
void reorder_matrix(msieve_obj *obj, 
		    uint32 **rowperm, 
		    uint32 **colperm) {

	graph_t graph;
	heap_t heap;

	logprintf(obj, "permuting matrix for faster multiplies\n");

	graph_init(obj, &graph, &heap);

	do_partition(obj, &graph, &heap,
			graph.num_cols, 
			graph.num_cols + graph.num_rows - 1,
			0, 
			graph.num_cols - 1);

	graph_free(&graph, rowperm, colperm);
}
