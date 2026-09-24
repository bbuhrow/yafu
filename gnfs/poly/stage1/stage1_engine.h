/*--------------------------------------------------------------------
Runtime registry of stage-1 collision engines.

All engines share the polysize contract:
    <engine>_thread_data_init/free   (thread_control init/shutdown shape)
    <engine> per-run hw context       (optional; GPU only)
    stage1_specialq_<engine>(task, threadid, q_min, q_max, p_min, p_max)
    -> handle_collision(task, ...) into the shared stage-2 pool

HAVE_CUDA still gates whether the GPU engines are *compiled*; this registry
gates which compiled engine is *selected* at run time. A CPU-only build
(no CUDA toolkit) registers only the CPU engines and links no GPU symbols.
--------------------------------------------------------------------*/

#ifndef _STAGE1_ENGINE_H_
#define _STAGE1_ENGINE_H_

#include <stage1.h>   /* task_data_t, stage1_sieve_data_t, msieve_obj, uint* */

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
	STAGE1_ENGINE_CPU_HASHTABLE = 0, /* #1 trunk hashtable (wrapped)      */
	STAGE1_ENGINE_CPU_GERBICZ,       /* #2 polysize Gerbicz-sort (ref)    */
	STAGE1_ENGINE_GPU_CUBSORT,       /* #3 trunk CUB radix-sort           */
	STAGE1_ENGINE_GPU_GERBICZ,       /* #4 kyleaskine collision engine    */
	STAGE1_NUM_ENGINES
} stage1_engine_id;

/* the cells an engine can accept. The driver routes around anything
   outside this instead of letting the engine crash: a_d whose p_max
   exceeds max_p is skipped; a q window above max_special_q is clamped. */
typedef struct {
	uint64 max_special_q;   /* inclusive cap on the special-q window */
	uint32 max_p;           /* inclusive cap on p                    */
	uint32 is_gpu;          /* device residency                      */
} stage1_envelope_t;

typedef enum {
	STAGE1_OVERFLOW_NONE = 0,  /* engine cannot fail a within-envelope cell */
	STAGE1_OVERFLOW_SKIP       /* engine may fail a cell; treat as a skip    */
} stage1_overflow_policy;

struct stage1_engine_vtable {
	const char *name;              /* selection token, e.g. "cpu_gerbicz" */
	stage1_engine_id id;
	stage1_envelope_t envelope;
	stage1_overflow_policy overflow;
	uint32 max_threads;            /* host-thread cap (GPU=4); 0 = uncapped */

	/* per-run hardware context -> stored in d->hw_data; NULL if none */
	void* (*sieve_data_init)(msieve_obj* obj, uint32 num_threads, uint32 id);
	void   (*sieve_data_free)(void *hw_data);

	/* per-thread working set; data = stage1_sieve_data_t *, allocate into
	   d->threads[threadid].hw_thread_data */
	void (*thread_data_init)(void *data, int threadid);
	void (*thread_data_free)(void *data, int threadid);

	/* the worker; ranges are forced by the caller. Engines that share a
	   worker (GPU #3/#4) read task->d->engine->id to pick their method. */
	void (*specialq)(task_data_t *task, uint32 threadid,
			uint64 special_q_min, uint64 special_q_max,
			uint32 p_min, uint32 p_max);
};
typedef struct stage1_engine_vtable stage1_engine_vtable_t;

/* NULL if the id is out of range or the engine was not compiled in */
const stage1_engine_vtable_t * stage1_engine_lookup(stage1_engine_id id);
const stage1_engine_vtable_t * stage1_engine_by_name(const char *name);

/* pick the engine for this run from obj->nfs_args ("stage1_engine=<name>"),
   defaulting to the CPU Gerbicz reference. Exits with the list of compiled
   engines if an unknown/unavailable engine is requested. */
const stage1_engine_vtable_t * stage1_engine_select(msieve_obj *obj);

/* 1 if the cell fits the engine (clamping *special_q_max into the envelope
   in place); 0 if p_max exceeds the engine's cap and the a_d must be skipped */
int stage1_engine_cell_fits(const stage1_engine_vtable_t *v,
			uint32 p_max, uint64 *special_q_max);

#ifdef __cplusplus
}
#endif

#endif /* !_STAGE1_ENGINE_H_ */
