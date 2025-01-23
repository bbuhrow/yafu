/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: singleton.c 1049 2023-09-14 19:57:48Z jasonp_sf $
--------------------------------------------------------------------*/

#include "filter_priv.h"

/*--------------------------------------------------------------------*/
static void filter_read_lp_file_1pass(msieve_obj *obj, 
				filter_t *filter,
				uint32 max_ideal_weight) {

	uint32 i, j, k;
	FILE *fp;
	char buf[256];
	uint32 num_relations = filter->num_relations;
	uint32 num_ideals = filter->num_ideals;
	uint32 *counts;
	relation_ideal_t *r, *r_old, *r_next;

	logprintf(obj, "reading all ideals from disk\n");

	sprintf(buf, "%s.lp", obj->savefile.name);
	fp = fopen(buf, "rb");
	if (fp == NULL) {
		logprintf(obj, "error: can't open LP file\n");
		exit(-1);
	}

	filter->relation_array = (relation_ideal_t *)xmalloc(
					(size_t)filter->lp_file_size);
	fread(filter->relation_array, (size_t)1, 
			(size_t)filter->lp_file_size, fp);

	fclose(fp);
	logprintf(obj, "memory use: %.1f MB\n", 
			(double)filter->lp_file_size / 1048576);

	/* build a frequency table for the ideals */

	counts = (uint32 *)xcalloc((size_t)num_ideals, sizeof(uint32));
	r = filter->relation_array;
	for (i = 0; i < num_relations; i++) {

		for (j = 0; j < r->ideal_count; j++)
			counts[r->ideal_list[j]]++;

		r = next_relation_ptr(r);
	}

	/* renumber the ideals to ignore the ones that occur
	   too often */

	for (i = j = 0; i < num_ideals; i++) {
		if (counts[i] == 0 || counts[i] > max_ideal_weight)
			counts[i] = (uint32)(-1);
		else
			counts[i] = j++;
	}
	filter->target_excess += i - j;
	filter->num_ideals = j;

	if (i != j) {
		logprintf(obj, "keeping %u ideals with weight <= %u, "
				"target excess is %u\n",
				j, max_ideal_weight, 
				filter->target_excess);

		r = r_old = filter->relation_array;
		for (i = 0; i < num_relations; i++) {

			uint32 ideal_count = r->ideal_count;

			r_next = next_relation_ptr(r);
			r_old->rel_index = r->rel_index;
			r_old->gf2_factors = r->gf2_factors;

			for (j = k = 0; j < ideal_count; j++) {

				uint32 ideal = counts[r->ideal_list[j]];

				if (ideal != (uint32)(-1))
					r_old->ideal_list[k++] = ideal;
			}
			r_old->ideal_count = k;
			r_old->gf2_factors += ideal_count - k;
			r_old = next_relation_ptr(r_old);
			r = r_next;
		}

		/* trim the allocated memory */

		filter->relation_array = (relation_ideal_t *)xrealloc(
					filter->relation_array, 
					(size_t)((uint32 *)r_old - 
					(uint32 *)filter->relation_array) *
					sizeof(uint32));
	}

	free(counts);
}

/*--------------------------------------------------------------------*/
void filter_read_lp_file(msieve_obj *obj, filter_t *filter,
				uint32 max_ideal_weight) {
	uint32 i, j, k;
	FILE *fp;
	char buf[256];
	size_t header_words;
	relation_ideal_t tmp;
	relation_ideal_t *relation_array;
	uint64 curr_word;
	size_t num_relation_alloc;
	uint32 *counts;
	uint32 num_relations = filter->num_relations;
	uint32 num_ideals = filter->num_ideals;
	size_t mem_use;

	/* read in the relations from the lp file, as well as the large
	   ideals they contain. Do not save ideals that occur
	   more than max_ideal_weight times in the dataset */

	if (max_ideal_weight == 0) {
		filter_read_lp_file_1pass(obj, filter, 200);
		filter_purge_singletons_core(obj, filter);
		return;
	}

	logprintf(obj, "reading large ideals from disk\n");

	sprintf(buf, "%s.lp", obj->savefile.name);
	fp = fopen(buf, "rb");
	if (fp == NULL) {
		logprintf(obj, "error: singleton2 can't open LP file\n");
		exit(-1);
	}

	header_words = (sizeof(relation_ideal_t) - 
			sizeof(tmp.ideal_list)) / sizeof(uint32);
	counts = (uint32 *)xcalloc((size_t)num_ideals, sizeof(uint32));

	/* first build a frequency table for the large ideals */

	for (i = 0; i < num_relations; i++) {

		fread(&tmp, sizeof(uint32), header_words, fp);

		for (j = 0; j < tmp.ideal_count; j++) {
			uint32 curr_ideal;

			fread(&curr_ideal, sizeof(uint32), (size_t)1, fp);
			counts[curr_ideal]++;
		}
	}

	/* renumber the ideals to ignore the ones that occur
	   too often */

	for (i = j = 0; i < num_ideals; i++) {
		if (counts[i] <= max_ideal_weight)
			counts[i] = j++;
		else
			counts[i] = (uint32)(-1);
	}
	filter->target_excess += i - j;
	filter->num_ideals = j;
	logprintf(obj, "keeping %u ideals with weight <= %u, "
			"target excess is %u\n",
			j, max_ideal_weight, 
			filter->target_excess);

	/* reread the relation list, saving the sparse ideals */

	rewind(fp);
	num_relation_alloc = 10000;
	curr_word = 0;
	relation_array = (relation_ideal_t *)xmalloc(
					num_relation_alloc *
					sizeof(relation_ideal_t));
	for (i = 0; i < num_relations; i++) {

		relation_ideal_t *r;

		/* make sure the relation array has room for the
		   new relation. Be careful increasing the array
		   size, since this is probably the largest array
		   in the NFS code */

		if (curr_word * sizeof(uint32) >=
				(num_relation_alloc-1) * 
				sizeof(relation_ideal_t)) {

			num_relation_alloc = 1.4 * num_relation_alloc;
			relation_array = (relation_ideal_t *)xrealloc(
					relation_array, 
					num_relation_alloc *
					sizeof(relation_ideal_t));
		}

		r = (relation_ideal_t *)(
			(uint32 *)relation_array + curr_word);
		fread(r, sizeof(uint32), header_words, fp);

		for (j = k = 0; j < r->ideal_count; j++) {

			uint32 curr_ideal;

			fread(&curr_ideal, sizeof(uint32), (size_t)1, fp);
			curr_ideal = counts[curr_ideal];
			if (curr_ideal != (uint32)(-1))
				r->ideal_list[k++] = curr_ideal;
		}
		r->gf2_factors += j - k;
		r->ideal_count = k;
		curr_word += header_words + k;
	}

	/* finish up: trim the allocated relation array */

	filter->relation_array = (relation_ideal_t *)xrealloc(
						relation_array, 
						curr_word * 
						sizeof(uint32));
	free(counts);
	fclose(fp);
	mem_use = num_ideals * sizeof(uint32) +
			num_relation_alloc * 
			sizeof(relation_ideal_t);
	logprintf(obj, "memory use: %.1f MB\n", 
			(double)mem_use / 1048576);

	/* get rid of all the singletons now */

	filter_purge_singletons_core(obj, filter);
}

/*--------------------------------------------------------------------*/
void filter_purge_lp_singletons(msieve_obj *obj, 
				filter_t *filter,
				uint64 ram_size) {

	uint32 i, j, k, m;
	FILE *in_fp;
	FILE *out_fp;
	char buf[256];
	char buf2[256];
	size_t header_words;
	relation_ideal_t tmp;
	uint32 *relation_num;
	uint32 *counts;
	uint32 num_singletons;
	uint32 num_relations = filter->num_relations;
	uint32 num_ideals = filter->num_ideals;
	uint32 start_relations = filter->num_relations;
	uint32 num_passes = 0;
	uint64 new_file_size;

	logprintf(obj, "removing singletons from LP file\n");
	logprintf(obj, "start with %u relations and %u ideals\n",
			num_relations, num_ideals);

	sprintf(buf, "%s.lp", obj->savefile.name);
	in_fp = fopen(buf, "rb");
	if (in_fp == NULL) {
		logprintf(obj, "error: can't open LP file\n");
		exit(-1);
	}
	sprintf(buf2, "%s.lp0", obj->savefile.name);
	out_fp = fopen(buf2, "wb");
	if (out_fp == NULL) {
		logprintf(obj, "error: can't open LP output file\n");
		exit(-1);
	}

	header_words = (sizeof(relation_ideal_t) - 
			sizeof(tmp.ideal_list)) / sizeof(uint32);
	relation_num = (uint32 *)xmalloc(num_relations * sizeof(uint32));
	counts = (uint32 *)xcalloc((size_t)num_ideals, sizeof(uint32));

	/* first build a frequency table for the large ideals */

	for (i = 0; i < num_relations; i++) {

		uint32 *ideal_list = tmp.ideal_list;

		relation_num[i] = i;
		fread(&tmp, sizeof(uint32), header_words, in_fp);
		fread(ideal_list, sizeof(uint32), 
			(size_t)tmp.ideal_count, in_fp);

		for (j = 0; j < tmp.ideal_count; j++)
			counts[ideal_list[j]]++;
	}
	rewind(in_fp);

	/* iteratively ignore relations that contain singleton ideals;
	   we want to limit the number of passes over the disk file,
	   so stop when the number of relations in the file will easily
	   fit in memory, or the number of singletons starts decaying */

	do {
		new_file_size = 0;
		rewind(in_fp);

		for (i = j = k = 0; i < start_relations; i++) {

			uint32 *ideal_list = tmp.ideal_list;

			fread(&tmp, sizeof(uint32), header_words, in_fp);
			fread(ideal_list, sizeof(uint32), 
					(size_t)tmp.ideal_count, in_fp);

			if (i == relation_num[k]) {
				for (m = 0; m < tmp.ideal_count; m++) {
					if (counts[ideal_list[m]] < 2)
						break;
				}

				if (m == tmp.ideal_count) {
					relation_num[j++] = relation_num[k];
					new_file_size += (header_words +
							tmp.ideal_count) *
							sizeof(uint32);
				}
				else {
					for (m = 0; m < tmp.ideal_count; m++)
						counts[ideal_list[m]]--;
				}

				if (++k == num_relations)
					break;
			}
		}
		num_relations = j;
		num_singletons = k - j;
		logprintf(obj, "pass %u: found %u singletons\n",
				++num_passes, num_singletons);
		rewind(in_fp);

	} while (num_relations > 2000000 && 
			num_singletons > 500000 &&
			new_file_size >= ram_size / 2);


	/* renumber the ideals to squeeze out the removed ones */

	for (i = j = 0; i < num_ideals; i++) {
		if (counts[i] != 0)
			counts[i] = j++;
	}
	num_ideals = j;

	/* reread the relation list, saving relations that survived
	   singleton removal and renumbering their ideals */

	for (i = j = 0; i < start_relations; i++) {

		uint32 *ideal_list = tmp.ideal_list;

		fread(&tmp, sizeof(uint32), header_words, in_fp);
		fread(ideal_list, sizeof(uint32), 
				(size_t)tmp.ideal_count, in_fp);

		if (i == relation_num[j]) {
			for (k = 0; k < tmp.ideal_count; k++)
				ideal_list[k] = counts[ideal_list[k]];

			fwrite(&tmp, sizeof(uint32),
				header_words + tmp.ideal_count, out_fp);

			if (++j == num_relations)
				break;
		}
	}

	logprintf(obj, "pruned dataset has %u relations and "
			"%u large ideals\n", num_relations, num_ideals);

	filter->num_relations = num_relations;
	filter->num_ideals = num_ideals;
	filter->relation_array = NULL;
	free(counts);
	free(relation_num);

	fclose(in_fp);
	fclose(out_fp);
	if (remove(buf) != 0) {
		logprintf(obj, "error: can't delete LP file\n");
		exit(-1);
	}
	if (rename(buf2, buf) != 0) {
		logprintf(obj, "error: can't rename LP output file\n");
		exit(-1);
	}
	filter->lp_file_size = get_file_size(buf);
}

/*--------------------------------------------------------------------*/
void filter_purge_singletons_core(msieve_obj *obj, 
				filter_t *filter) {

	/* main routine for performing in-memory singleton
	   removal. We iterate until there are no more singletons */

	uint32 i, j;
	uint32 *freqtable;
	relation_ideal_t *relation_array;
	relation_ideal_t *curr_relation;
	relation_ideal_t *old_relation;
	uint32 orig_num_ideals;
	uint32 num_passes;
	uint32 num_relations;
	uint32 num_ideals;
	uint32 new_num_relations;

	logprintf(obj, "commencing in-memory singleton removal\n");

	num_relations = filter->num_relations;
	orig_num_ideals = num_ideals = filter->num_ideals;
	relation_array = filter->relation_array;
	freqtable = (uint32 *)xcalloc((size_t)num_ideals, sizeof(uint32));

	/* count the number of times each ideal occurs. Note
	   that since we know the exact number of ideals, we
	   don't need a hashtable to store the counts, just an
	   ordinary random-access array (i.e. a perfect hashtable) */

	curr_relation = relation_array;
	for (i = 0; i < num_relations; i++) {
		for (j = 0; j < curr_relation->ideal_count; j++) {
			uint32 ideal = curr_relation->ideal_list[j];
			freqtable[ideal]++;
		}
		curr_relation = next_relation_ptr(curr_relation);
	}

	logprintf(obj, "begin with %u relations and %u unique ideals\n", 
					num_relations, num_ideals);

	/* while singletons were found */

	num_passes = 0;
	new_num_relations = num_relations;
	do {
		num_relations = new_num_relations;
		new_num_relations = 0;
		curr_relation = relation_array;
		old_relation = relation_array;

		for (i = 0; i < num_relations; i++) {
			uint32 curr_num_ideals = curr_relation->ideal_count;
			uint32 ideal;
			relation_ideal_t *next_relation;

			/* the ideal count in curr_relation may get
			   overwritten when writing old_relation, so
			   cache the count and point to the next
			   relation now */

			next_relation = next_relation_ptr(curr_relation);

			/* check the count of each ideal */

			for (j = 0; j < curr_num_ideals; j++) {
				ideal = curr_relation->ideal_list[j];
				if (freqtable[ideal] <= 1)
					break;
			}

			if (j < curr_num_ideals) {

				/* relation is a singleton; decrement the
				   count of each of its ideals and skip it */

				for (j = 0; j < curr_num_ideals; j++) {
					ideal = curr_relation->ideal_list[j];
					freqtable[ideal]--;
				}
			}
			else {
				/* relation survived this pass; append it to
				   the list of survivors */

				old_relation->rel_index = 
						curr_relation->rel_index;
				old_relation->gf2_factors = 
						curr_relation->gf2_factors;
				old_relation->ideal_count = curr_num_ideals;
				for (j = 0; j < curr_num_ideals; j++) {
					old_relation->ideal_list[j] =
						curr_relation->ideal_list[j];
				}
				new_num_relations++;
				old_relation = next_relation_ptr(old_relation);
			}

			curr_relation = next_relation;
		}

		num_passes++;
	} while (new_num_relations != num_relations);

	/* find the ideal that occurs in the most
	   relations, and renumber the ideals to ignore
	   any that have a count of zero */

	num_ideals = 0;
	for (i = j = 0; i < orig_num_ideals; i++) {
		if (freqtable[i]) {
			j = MAX(j, freqtable[i]);
			freqtable[i] = num_ideals++;
		}
	}

	logprintf(obj, "reduce to %u relations and %u ideals in %u passes\n", 
				num_relations, num_ideals, num_passes);
	logprintf(obj, "max relations containing the same ideal: %u\n", j);
	
	/* save the current state */

	filter->max_ideal_weight = j;
	filter->num_relations = num_relations;
	filter->num_ideals = num_ideals;
	filter->relation_array = relation_array = 
			(relation_ideal_t *)xrealloc(relation_array,
				(old_relation - relation_array + 1) *
				sizeof(relation_ideal_t));

	curr_relation = relation_array;
	for (i = 0; i < num_relations; i++) {
		for (j = 0; j < curr_relation->ideal_count; j++) {
			uint32 ideal = curr_relation->ideal_list[j];
			curr_relation->ideal_list[j] = freqtable[ideal];
		}
		curr_relation = next_relation_ptr(curr_relation);
	}
	free(freqtable);
}
