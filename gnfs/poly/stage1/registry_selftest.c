/* Self-contained check of the stage1 engine-registry MECHANISM: the
   designated-initializer table with #ifdef holes, availability detection,
   name lookup, default selection, unavailable-engine rejection, and the
   envelope clamp/skip. Types and engine functions are mocked so this
   compiles without the msieve tree; the table + logic are structurally
   identical to stage1_engine.c. Build 4 ways: [-DHAVE_CUDA] x
   [-DHAVE_CPU_HASHTABLE]. */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

typedef uint32_t uint32;
typedef uint64_t uint64;
#define MAX_P          (((uint32)1u << 31) - 1u)
#define MAX_SPECIAL_Q  ((uint64)0x7FFFFFFFFFFFFFFFULL)

typedef struct { char *nfs_args; } msieve_obj;
typedef struct { int dummy; } task_data_t;

/* ---- mock engine entry points ---- */
static void   m_cpu_ti(void *d,int t){(void)d;(void)t;}
static void   m_cpu_tf(void *d,int t){(void)d;(void)t;}
static void   m_cpu_sq(task_data_t*a,uint32 b,uint64 c,uint64 e,uint32 f,uint32 g)
              {(void)a;(void)b;(void)c;(void)e;(void)f;(void)g;}
#ifdef HAVE_CPU_HASHTABLE
static void   m_ht_ti(void *d,int t){(void)d;(void)t;}
static void   m_ht_tf(void *d,int t){(void)d;(void)t;}
static void   m_ht_sq(task_data_t*a,uint32 b,uint64 c,uint64 e,uint32 f,uint32 g)
              {(void)a;(void)b;(void)c;(void)e;(void)f;(void)g;}
#endif
#ifdef HAVE_CUDA
static void  *m_gpu_sdi(msieve_obj*o,uint32 n){(void)o;(void)n;return (void*)1;}
static void   m_gpu_sdf(void *h){(void)h;}
static void   m_gpu_ti(void *d,int t){(void)d;(void)t;}
static void   m_gpu_tf(void *d,int t){(void)d;(void)t;}
static void   m_gpu_sq(task_data_t*a,uint32 b,uint64 c,uint64 e,uint32 f,uint32 g)
              {(void)a;(void)b;(void)c;(void)e;(void)f;(void)g;}
#endif

/* ==== begin logic copied from stage1_engine.{h,c} ==== */
typedef enum { STAGE1_ENGINE_CPU_HASHTABLE=0, STAGE1_ENGINE_CPU_GERBICZ,
	STAGE1_ENGINE_GPU_CUBSORT, STAGE1_ENGINE_GPU_GERBICZ,
	STAGE1_NUM_ENGINES } stage1_engine_id;
typedef struct { uint64 max_special_q; uint32 max_p; uint32 is_gpu; } stage1_envelope_t;
typedef enum { STAGE1_OVERFLOW_NONE=0, STAGE1_OVERFLOW_SKIP } stage1_overflow_policy;
typedef struct {
	const char *name; stage1_engine_id id; stage1_envelope_t envelope;
	stage1_overflow_policy overflow; uint32 max_threads;
	void *(*sieve_data_init)(msieve_obj*,uint32); void (*sieve_data_free)(void*);
	void (*thread_data_init)(void*,int); void (*thread_data_free)(void*,int);
	void (*specialq)(task_data_t*,uint32,uint64,uint64,uint32,uint32);
} stage1_engine_vtable_t;

#define ENV_TRUNK_P ((uint32)((1u<<27)-1u))
#define ENV_TRUNK_Q ((uint64)0xFFFFFFFFu)
#define ENV_FULL_P  ((uint32)MAX_P)
#define ENV_FULL_Q  ((uint64)MAX_SPECIAL_Q)

static const stage1_engine_vtable_t stage1_engines[STAGE1_NUM_ENGINES] = {
	[STAGE1_ENGINE_CPU_GERBICZ] = { "cpu_gerbicz", STAGE1_ENGINE_CPU_GERBICZ,
		{ENV_FULL_Q,ENV_FULL_P,0}, STAGE1_OVERFLOW_NONE, 0,
		NULL, NULL, m_cpu_ti, m_cpu_tf, m_cpu_sq },
#ifdef HAVE_CPU_HASHTABLE
	[STAGE1_ENGINE_CPU_HASHTABLE] = { "cpu_hashtable", STAGE1_ENGINE_CPU_HASHTABLE,
		{ENV_TRUNK_Q,ENV_TRUNK_P,0}, STAGE1_OVERFLOW_NONE, 0,
		NULL, NULL, m_ht_ti, m_ht_tf, m_ht_sq },
#endif
#ifdef HAVE_CUDA
	[STAGE1_ENGINE_GPU_CUBSORT] = { "gpu_cubsort", STAGE1_ENGINE_GPU_CUBSORT,
		{ENV_TRUNK_Q,ENV_TRUNK_P,1}, STAGE1_OVERFLOW_NONE, 4,
		m_gpu_sdi, m_gpu_sdf, m_gpu_ti, m_gpu_tf, m_gpu_sq },
	[STAGE1_ENGINE_GPU_GERBICZ] = { "gpu_gerbicz", STAGE1_ENGINE_GPU_GERBICZ,
		{ENV_TRUNK_Q,ENV_TRUNK_P,1}, STAGE1_OVERFLOW_SKIP, 4,
		m_gpu_sdi, m_gpu_sdf, m_gpu_ti, m_gpu_tf, m_gpu_sq },
#endif
};
static int engine_available(const stage1_engine_vtable_t *v){
	return v!=NULL && v->specialq!=NULL; }
static const stage1_engine_vtable_t *by_name(const char*n){
	uint32 i; for(i=0;i<(uint32)STAGE1_NUM_ENGINES;i++){
		const stage1_engine_vtable_t*v=&stage1_engines[i];
		if(engine_available(v)&&strcmp(v->name,n)==0)return v;} return NULL; }
static int cell_fits(const stage1_engine_vtable_t*v,uint32 p_max,uint64*q){
	if(p_max>v->envelope.max_p)return 0;
	if(*q>v->envelope.max_special_q)*q=v->envelope.max_special_q;
	return 1;
}
/* ==== end copied logic ==== */

int main(void){
	uint32 i;
#if defined(HAVE_CUDA) && defined(HAVE_CPU_HASHTABLE)
	const char *cfg="HAVE_CUDA + HAVE_CPU_HASHTABLE (all four)";
#elif defined(HAVE_CUDA)
	const char *cfg="HAVE_CUDA (cpu_gerbicz + gpu #3/#4)";
#elif defined(HAVE_CPU_HASHTABLE)
	const char *cfg="HAVE_CPU_HASHTABLE (cpu #1/#2)";
#else
	const char *cfg="cpu-only (cpu_gerbicz reference only)";
#endif
	printf("=== build: %s ===\n", cfg);

	printf("registered:");
	for(i=0;i<(uint32)STAGE1_NUM_ENGINES;i++)
		if(engine_available(&stage1_engines[i]))
			printf(" %s", stage1_engines[i].name);
	printf("\n");

	const char *probe[]={"cpu_gerbicz","cpu_hashtable","gpu_cubsort","gpu_gerbicz"};
	for(i=0;i<4;i++){
		const stage1_engine_vtable_t*v=by_name(probe[i]);
		printf("  select %-14s -> %s\n", probe[i],
			v ? "OK" : "REJECT (driver would exit + list availables)");
	}

	/* envelope: full engine accepts a big cell unchanged; a trunk-capped
	   engine skips an over-p a_d and clamps an over-q window */
	{
		const stage1_engine_vtable_t*full=by_name("cpu_gerbicz");
		uint64 q=(uint64)1e18; uint32 pmax=1000000u;
		int ok=cell_fits(full,pmax,&q);
		printf("  cell_fits cpu_gerbicz p=%u q=1e18 -> fit=%d q'=%llu%s\n",
			pmax,ok,(unsigned long long)q,q==(uint64)1e18?" (unchanged)":"");

		/* synthesize a trunk-capped envelope to exercise skip+clamp
		   regardless of which engines this config compiled */
		stage1_engine_vtable_t capped=*full;
		capped.envelope.max_p=ENV_TRUNK_P; capped.envelope.max_special_q=ENV_TRUNK_Q;
		q=(uint64)5e9;
		int skip=cell_fits(&capped,200000000u,&q);   /* p over 2^27 */
		printf("  cell_fits trunk-cap  p=2e8 -> fit=%d (expect 0, skip)\n",skip);
		q=(uint64)5e9;
		int clamp=cell_fits(&capped,100000u,&q);      /* q over 2^32 */
		printf("  cell_fits trunk-cap  p=1e5 q=5e9 -> fit=%d q'=%llu (clamped to 2^32-1)\n",
			clamp,(unsigned long long)q);
	}
	return 0;
}
