/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: filter.c 1020 2018-05-18 13:33:21Z jasonp_sf $
--------------------------------------------------------------------*/

#include "filter_priv.h"

/*--------------------------------------------------------------------*/
void filter_free_relsets(merge_t *merge) {

	uint32 i;
	relation_set_t *relset_array = merge->relset_array;
	uint32 num_relsets = merge->num_relsets;

	for (i = 0; i < num_relsets; i++) {
		relation_set_t *r = relset_array + i;
		free(r->data);
	}
	free(merge->relset_array);
	merge->relset_array = NULL;
}

/*--------------------------------------------------------------------*/
void filter_dump_relsets(msieve_obj *obj, merge_t *merge) {

	uint32 i;
	relation_set_t *relset_array = merge->relset_array;
	uint32 num_relsets = merge->num_relsets;
	char buf[256];
	FILE *cycle_fp;

	sprintf(buf, "%s.cyc", obj->savefile.name);
	cycle_fp = fopen(buf, "wb");
	if (cycle_fp == NULL) {
		logprintf(obj, "error: can't open cycle file\n");
		exit(-1);
	}

	fwrite(&num_relsets, sizeof(uint32), (size_t)1, cycle_fp);

	for (i = 0; i < num_relsets; i++) {
		relation_set_t *r = relset_array + i;
		uint32 num = r->num_relations;

		fwrite(&num, sizeof(uint32), (size_t)1, cycle_fp);
		fwrite(r->data, sizeof(uint32),
				(size_t)num, cycle_fp);
	}
	fclose(cycle_fp);
}

/*--------------------------------------------------------------------*/
int32 filter_make_relsets(msieve_obj *obj, filter_t *filter,
				merge_t *merge, uint32 min_cycles) {

	filter_purge_cliques(obj, filter);
	filter_merge_init(obj, filter);
	filter_merge_2way(obj, filter, merge);
	return filter_merge_full(obj, merge, min_cycles);
}
