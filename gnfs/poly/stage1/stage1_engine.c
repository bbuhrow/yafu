/*--------------------------------------------------------------------
Stage-1 collision-engine registry (see stage1_engine.h).
--------------------------------------------------------------------*/

#include <stage1_engine.h>
#include <string.h>
#include <stdlib.h>

/*--------------------------------------------------------------------
Engine entry points implemented elsewhere, declared here so the table
can reference them.

  #2 (cpu_gerbicz)   : always compiled -- cpu_thread_data_init/free and
                       stage1_specialq_cpu live in stage1_sieve_cpu.c and
                       are already declared in stage1.h.

  #1 (cpu_hashtable) : the trunk wrapper (build item 2). Behind
                       HAVE_CPU_HASHTABLE until that wrapper lands.

  #3/#4 (gpu_*)      : gpu_sieve_data_init/free, gpu_thread_data_init/free
                       and stage1_specialq_gpu are declared in stage1.h and
                       provided by the GPU port (build item 3) only in a
                       CUDA build. #3 and #4 share the worker; the worker
                       selects CUB vs the collision engine from
                       task->d->engine->id.
--------------------------------------------------------------------*/

#ifdef HAVE_CPU_HASHTABLE
void cpu_hashtable_thread_data_init(void *data, int threadid);
void cpu_hashtable_thread_data_free(void *data, int threadid);
void stage1_specialq_cpu_hashtable(task_data_t *task, uint32 threadid,
			uint64 special_q_min, uint64 special_q_max,
			uint32 p_min, uint32 p_max);
#endif

/* envelope caps.
   trunk (#1/#3/#4): p is a 27-bit packed field, special-q is uint32.
   full  (#2)      : MAX_P (2^31-1) and MAX_SPECIAL_Q (2^63-1) per stage1.h */
#define ENV_TRUNK_P  ((uint32)((1u << 27) - 1u))
#define ENV_TRUNK_Q  ((uint64)0xFFFFFFFFu)
#define ENV_FULL_P   ((uint32)MAX_P)
#define ENV_FULL_Q   ((uint64)MAX_SPECIAL_Q)

static const stage1_engine_vtable_t stage1_engines[STAGE1_NUM_ENGINES] = {

	/* #2 -- the reference; always present */
	[STAGE1_ENGINE_CPU_GERBICZ] = {
		"cpu_gerbicz", STAGE1_ENGINE_CPU_GERBICZ,
		{ ENV_FULL_Q, ENV_FULL_P, 0 },
		STAGE1_OVERFLOW_NONE, 0,
		NULL, NULL,
		cpu_thread_data_init, cpu_thread_data_free,
		stage1_specialq_cpu,
	},

#ifdef HAVE_CPU_HASHTABLE
	/* #1 -- trunk hashtable, wrapped (build item 2) */
	[STAGE1_ENGINE_CPU_HASHTABLE] = {
		"cpu_hashtable", STAGE1_ENGINE_CPU_HASHTABLE,
		{ ENV_TRUNK_Q, ENV_TRUNK_P, 0 },
		STAGE1_OVERFLOW_NONE, 0,
		NULL, NULL,
		cpu_hashtable_thread_data_init, cpu_hashtable_thread_data_free,
		stage1_specialq_cpu_hashtable,
	},
#endif

#ifdef HAVE_CUDA_POLY
	/* #3 -- trunk CUB radix-sort (build item 3) */
	[STAGE1_ENGINE_GPU_CUBSORT] = {
		"gpu_cubsort", STAGE1_ENGINE_GPU_CUBSORT,
		{ ENV_TRUNK_Q, ENV_TRUNK_P, 1 },
		STAGE1_OVERFLOW_NONE, 4,
		gpu_sieve_data_init, gpu_sieve_data_free,
		gpu_thread_data_init, gpu_thread_data_free,
		stage1_specialq_gpu,
	},

	/* #4 -- kyleaskine collision engine; shares the GPU worker, may
	   overflow a dense cell -> treat as a graceful skip (build item 3) */
	[STAGE1_ENGINE_GPU_GERBICZ] = {
		"gpu_gerbicz", STAGE1_ENGINE_GPU_GERBICZ,
		{ ENV_TRUNK_Q, ENV_TRUNK_P, 1 },
		STAGE1_OVERFLOW_SKIP, 4,
		gpu_sieve_data_init, gpu_sieve_data_free,
		gpu_thread_data_init, gpu_thread_data_free,
		stage1_specialq_gpu,
	},
#endif
};

static int
engine_available(const stage1_engine_vtable_t *v)
{
	/* zero-initialized (compiled-out) rows have specialq == NULL */
	return v != NULL && v->specialq != NULL;
}

const stage1_engine_vtable_t *
stage1_engine_lookup(stage1_engine_id id)
{
	const stage1_engine_vtable_t *v;

	if ((int)id < 0 || id >= STAGE1_NUM_ENGINES)
		return NULL;
	v = &stage1_engines[id];
	return engine_available(v) ? v : NULL;
}

const stage1_engine_vtable_t *
stage1_engine_by_name(const char *name)
{
	uint32 i;

	for (i = 0; i < (uint32)STAGE1_NUM_ENGINES; i++) {
		const stage1_engine_vtable_t *v = &stage1_engines[i];
		if (engine_available(v) && strcmp(v->name, name) == 0)
			return v;
	}
	return NULL;
}

static void
list_available(msieve_obj *obj)
{
	uint32 i;

	for (i = 0; i < (uint32)STAGE1_NUM_ENGINES; i++) {
		const stage1_engine_vtable_t *v = &stage1_engines[i];
		if (engine_available(v))
			logprintf(obj, "  available stage 1 engine: %s%s\n",
				v->name,
				v->envelope.is_gpu ? " (gpu)" : " (cpu)");
	}
}

const stage1_engine_vtable_t *
stage1_engine_select(msieve_obj *obj)
{
	const stage1_engine_vtable_t *v;

	if (obj->nfs_args != NULL) {
		const char *tmp = strstr(obj->nfs_args, "stage1_engine=");

		if (tmp != NULL) {
			char name[64];
			uint32 n = 0;

			tmp += strlen("stage1_engine=");
			while (n < sizeof(name) - 1 && *tmp &&
					*tmp != ' ' && *tmp != ',')
				name[n++] = *tmp++;
			name[n] = 0;

			v = stage1_engine_by_name(name);
			if (v == NULL) {
				logprintf(obj, "error: stage 1 engine '%s' is "
					"not available in this build\n", name);
				list_available(obj);
				exit(-1);
			}
			logprintf(obj, "using stage 1 engine: %s\n", v->name);
			return v;
		}
	}

	/* default: the CPU Gerbicz reference (always compiled in) */
	v = stage1_engine_lookup(STAGE1_ENGINE_CPU_GERBICZ);
	if (v == NULL) {  /* should be unreachable */
		logprintf(obj, "error: no stage 1 engine available\n");
		exit(-1);
	}
	logprintf(obj, "using stage 1 engine: %s (default)\n", v->name);
	return v;
}

int
stage1_engine_cell_fits(const stage1_engine_vtable_t *v,
			uint32 p_max, uint64 *special_q_max)
{
	if (p_max > v->envelope.max_p)
		return 0;                                   /* skip this a_d */
	if (*special_q_max > v->envelope.max_special_q)
		*special_q_max = v->envelope.max_special_q; /* clamp window   */
	return 1;
}
