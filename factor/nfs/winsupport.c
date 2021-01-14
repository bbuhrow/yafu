/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

       				   --bbuhrow@gmail.com 12/6/2012
----------------------------------------------------------------------*/

#include "gnfs.h"
#include "yafu_string.h"
#include "arith.h"
#include "factor.h"

#ifdef USE_NFS

// windows builds need to supply the following functions because they 
// are not available within gnfs.lib or common.lib, but are available
// within libmsieve.a
#if defined(WIN32)

// windows machines also need these declarations for functions supplied
// in winsupport.c, for NFS factorizations using msieve.
void logprintf(msieve_obj *obj, char *fmt, ...);
void free_cycle_list(la_col_t *cycle_list, uint32 num_cycles);
uint32 merge_relations(uint32 *merge_array,
			uint32 *src1, uint32 n1,
			uint32 *src2, uint32 n2);
uint32 factor_list_add(msieve_obj *obj, factor_list_t *list, mp_t *new_factor);
void factor_list_add_core(msieve_obj *obj, factor_list_t *list, mp_t *new_factor);

#endif

msieve_obj * msieve_obj_new(char *input_integer, uint32 flags,
			    char *savefile_name, char *logfile_name,
			    char *nfs_fbfile_name,
			    uint32 seed1, uint32 seed2, uint32 max_relations,
			    enum cpu_type cpu,
			    uint32 cache_size1, uint32 cache_size2,
			    uint32 num_threads, uint32 which_gpu, 
			    const char *nfs_args) {

	msieve_obj *obj = (msieve_obj *)xcalloc((size_t)1, sizeof(msieve_obj));

	obj->input = input_integer;
	obj->flags = flags;
	obj->seed1 = seed1;
	obj->seed2 = seed2;
	obj->max_relations = max_relations;
	obj->cpu = cpu;
	obj->cache_size1 = cache_size1;
	obj->cache_size2 = cache_size2;
	obj->num_threads = num_threads;
	obj->which_gpu = which_gpu;
	obj->logfile_name = MSIEVE_DEFAULT_LOGFILE;
	obj->nfs_args = nfs_args;
	if (logfile_name)
		obj->logfile_name = logfile_name;
	obj->nfs_fbfile_name = MSIEVE_DEFAULT_NFS_FBFILE;
	if (nfs_fbfile_name)
		obj->nfs_fbfile_name = nfs_fbfile_name;
	obj->mp_sprintf_buf = (char *)xmalloc(32 * MAX_MP_WORDS + 1);
	savefile_init(&obj->savefile, savefile_name);
	
	return obj;
}

/*--------------------------------------------------------------------*/
msieve_obj * msieve_obj_free(msieve_obj *obj) {

	msieve_factor *curr_factor;

	curr_factor = obj->factors;
	while (curr_factor != NULL) {
		msieve_factor *next_factor = curr_factor->next;
		free(curr_factor->number);
		free(curr_factor);
		curr_factor = next_factor;
	}

	savefile_free(&obj->savefile);
	free(obj->mp_sprintf_buf);
	free(obj);
	return NULL;
}


void factor_list_init(factor_list_t *list) {

	memset(list, 0, sizeof(factor_list_t));
}

void logprintf(msieve_obj *obj, char *fmt, ...) {

	va_list ap;

	/* do *not* initialize 'ap' and use it twice; this
	   causes crashes on AMD64 */

	if (LOGFLAG && (obj->flags & MSIEVE_FLAG_USE_LOGFILE)) {
		time_t t = time(NULL);
		char buf[64];
#ifdef HAVE_MPI
		char namebuf[256];
		sprintf(namebuf, "%s.mpi%02u", 
				obj->logfile_name, obj->mpi_rank);
		FILE *logfile = fopen(namebuf, "a");
#else
		FILE *logfile = fopen(obj->logfile_name, "a");
#endif

		if (logfile == NULL) {
			printf("fopen error: %s\n", strerror(errno));
			fprintf(stderr, "cannot open logfile\n");
			exit(-1);
		}

		va_start(ap, fmt);
		buf[0] = 0;
		strcpy(buf, ctime(&t));
		*(strchr(buf, '\n')) = 0;
		fprintf(logfile, "%s  ", buf);
		vfprintf(logfile, fmt, ap);
		fclose(logfile);
		va_end(ap);
	}
	if (obj->flags & MSIEVE_FLAG_LOG_TO_STDOUT) {
		va_start(ap, fmt);
		vfprintf(stdout, fmt, ap);
		va_end(ap);
	}
}

void free_cycle_list(la_col_t *cycle_list, uint32 num_cycles) {

	uint32 i;

	for (i = 0; i < num_cycles; i++)
		free(cycle_list[i].cycle.list);
	free(cycle_list);
}

/*------------------------------------------------------------------*/
uint32 merge_relations(uint32 *merge_array,
		  uint32 *src1, uint32 n1,
		  uint32 *src2, uint32 n2) {

	/* Given two sorted lists of integers, merge
	   the lists into a single sorted list. We assume
	   each list contains no duplicate entries.
	   If a particular entry occurs in both lists,
	   don't add it to the final list at all. Returns 
	   the number of elements in the resulting list */

	uint32 i1, i2;
	uint32 num_merge;

	i1 = i2 = 0;
	num_merge = 0;

	while (i1 < n1 && i2 < n2) {
		uint32 val1 = src1[i1];
		uint32 val2 = src2[i2];

		if (val1 < val2) {
			merge_array[num_merge++] = val1;
			i1++;
		}
		else if (val1 > val2) {
			merge_array[num_merge++] = val2;
			i2++;
		}
		else {
			i1++;
			i2++;
		}
	}

	while (i1 < n1)
		merge_array[num_merge++] = src1[i1++];
	while (i2 < n2)
		merge_array[num_merge++] = src2[i2++];

	return num_merge;
}

#define mp_is_zero(a) ((a)->nwords == 0)
#define mp_is_one(a) ((a)->nwords == 1 && (a)->val[0] == 1)

uint32 factor_list_add(msieve_obj *obj, factor_list_t *list, 
				mp_t *new_factor) {

	if (!mp_is_zero(new_factor) && !mp_is_one(new_factor))
		factor_list_add_core(obj, list, new_factor);

	return factor_list_max_composite(list);
}

uint32 factor_list_max_composite(factor_list_t *list) {

	uint32 i, bits;

	/* Find the number of bits in the largest composite factor, 
	   and return that (to give calling code an estimate of 
	   how much work would be left if it stopped trying to 
	   find new factors now) */

	for (i = bits = 0; i < list->num_factors; i++) {
		final_factor_t *curr_factor = list->final_factors[i];

		if (curr_factor->type == MSIEVE_COMPOSITE) {
			uint32 curr_bits = mp_bits(&curr_factor->factor);
			bits = MAX(bits, curr_bits);
		}
	}

	return bits;
}

#define PRECOMPUTED_PRIME_BOUND 100000
static INLINE void mp_clear(mp_t *a) {
	memset(a, 0, sizeof(mp_t));
}

static INLINE void mp_copy(mp_t *a, mp_t *b) {
	*b = *a;
}

static INLINE int32 mp_cmp(const mp_t *a, const mp_t *b) {

	uint32 i;

	if (a->nwords > b->nwords)
		return 1;
	if (a->nwords < b->nwords)
		return -1;

	for (i = a->nwords; i ; i--) {
		if (a->val[i-1] > b->val[i-1])
			return 1;
		if (a->val[i-1] < b->val[i-1])
			return -1;
	}

	return 0;
}


static void factor_list_add_core(msieve_obj *obj, 
				factor_list_t *list, 
				mp_t *new_factor) {

	/* recursive routine to do the actual adding of factors.
	   Upon exit, the factors in 'list' will be mutually coprime */

	uint32 i;
	mp_t tmp1, tmp2, common, q, r;
	uint32 num_factors = list->num_factors;

	mp_clear(&tmp1); tmp1.nwords = tmp1.val[0] = 1;
	mp_clear(&tmp2); tmp2.nwords = tmp2.val[0] = 1;
	mp_clear(&common);

	/* compare new_factor to all the factors 
	   already in the list */

	for (i = 0; i < num_factors; i++) {
		final_factor_t *curr_factor = list->final_factors[i];

		/* skip new_factor if element i of the current
		   list would duplicate it */

		if (mp_cmp(new_factor, &curr_factor->factor) == 0)
			return;

		/* if new_factor has a factor C in common with element
		   i in the list, remove all instances of C, remove 
		   factor i, compress the list and postprocess */

		mp_gcd(new_factor, &curr_factor->factor, &common);
		if (!mp_is_one(&common)) {

			mp_copy(new_factor, &tmp1);
			while (1) {
				mp_divrem(&tmp1, &common, &q, &r);
				if (mp_is_zero(&q) || !mp_is_zero(&r))
					break;
				mp_copy(&q, &tmp1);
			}

			mp_copy(&curr_factor->factor, &tmp2);
			while (1) {
				mp_divrem(&tmp2, &common, &q, &r);
				if (mp_is_zero(&q) || !mp_is_zero(&r))
					break;
				mp_copy(&q, &tmp2);
			}

			free(list->final_factors[i]);
			list->final_factors[i] = 
				list->final_factors[--list->num_factors];
			break;
		}
	}

	if (i < num_factors) {

		/* there is overlap between new_factor and one
		   of the factors previously found. In the worst
		   case there are three new factors to deal with */

		if (!mp_is_one(&tmp1))
			factor_list_add_core(obj, list, &tmp1);
		if (!mp_is_one(&tmp2))
			factor_list_add_core(obj, list, &tmp2);
		if (!mp_is_one(&common))
			factor_list_add_core(obj, list, &common);
	}
	else {
		/* list doesn't need to be modified, except
		   to append new_factor. We go to some trouble to
		   avoid unnecessary primality tests when new_factor
		   is small */

		i = list->num_factors++;
		list->final_factors[i] = (final_factor_t *)xmalloc(
						sizeof(final_factor_t));
		if (new_factor->nwords <= 2 &&
		    ((uint64)new_factor->val[1] << 32 | 
					new_factor->val[0]) <
		    ((uint64)PRECOMPUTED_PRIME_BOUND * 
		      			PRECOMPUTED_PRIME_BOUND)) {
			list->final_factors[i]->type = MSIEVE_PRIME;
		}
		else {
			list->final_factors[i]->type = (mp_is_prime(
						new_factor, 
						&obj->seed1, &obj->seed2)) ?
						MSIEVE_PROBABLE_PRIME : 
						MSIEVE_COMPOSITE;
		}
		mp_copy(new_factor, &(list->final_factors[i]->factor));
	}
}

#endif
