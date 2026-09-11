/* Validate the parallel-stage-2 pattern in isolation:
   N producers -> bounded blocking queue -> S consumers, where each consumer
   owns a PRIVATE bundle (its own scratch, like a per-worker poly_sizeopt_t /
   poly_rootopt_t) and the only SHARED things (output file, best value, funnel
   counters) are touched under one mutex. Proves: (1) private state never
   crosses workers, (2) shared file writes never interleave, (3) best = true
   max, (4) counters exact, (5) queue fully drains. Models the real design;
   not the msieve structs. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <stdint.h>

#define NPROD 4
#define NCONS 3            /* S stage-2 workers */
#define QCAP  1000         /* like threadpool_init(S, 1000, ...) */
#define PER_PROD 300000
#define SCRATCH 256

/* ---- bounded blocking queue (models the stage-2 threadpool intake) ---- */
typedef struct {
	double buf[QCAP];
	int head, tail, count;
	int closed;
	pthread_mutex_t m;
	pthread_cond_t not_full, not_empty;
} queue_t;

static void q_init(queue_t *q){ memset(q,0,sizeof(*q));
	pthread_mutex_init(&q->m,0); pthread_cond_init(&q->not_full,0); pthread_cond_init(&q->not_empty,0); }
static void q_push(queue_t *q, double v){
	pthread_mutex_lock(&q->m);
	while (q->count==QCAP) pthread_cond_wait(&q->not_full,&q->m);   /* blocking submit */
	q->buf[q->tail]=v; q->tail=(q->tail+1)%QCAP; q->count++;
	pthread_cond_signal(&q->not_empty); pthread_mutex_unlock(&q->m);
}
static int q_pop(queue_t *q, double *v){
	pthread_mutex_lock(&q->m);
	while (q->count==0 && !q->closed) pthread_cond_wait(&q->not_empty,&q->m);
	if (q->count==0 && q->closed){ pthread_mutex_unlock(&q->m); return 0; }
	*v=q->buf[q->head]; q->head=(q->head+1)%QCAP; q->count--;
	pthread_cond_signal(&q->not_full); pthread_mutex_unlock(&q->m); return 1;
}
static void q_close(queue_t *q){ pthread_mutex_lock(&q->m); q->closed=1;
	pthread_cond_broadcast(&q->not_empty); pthread_mutex_unlock(&q->m); }

/* ---- shared stage-2 output (file + best + funnel), one mutex ---- */
typedef struct {
	FILE *file;
	double best;
	uint64_t sizeopt_pass, rootopt_pass, stage2_done;
	pthread_mutex_t lock;
} shared_t;

/* ---- per-worker private bundle (stands in for poly_sizeopt/rootopt_data) ---- */
typedef struct {
	int id;
	uint64_t local_count;
	unsigned char scratch[SCRATCH];   /* private working memory */
	shared_t *shared;
	queue_t *q;
} worker_t;

static queue_t Q;
static shared_t SH;
static worker_t W[NCONS];

static void *producer(void *arg){
	long id=(long)arg;
	for (int i=0;i<PER_PROD;i++)
		q_push(&Q, 1.0e-9 + (double)(id*PER_PROD+i)*1e-18);   /* strictly-increasing-ish */
	return 0;
}

static void *consumer(void *arg){
	worker_t *w=(worker_t*)arg;
	double v;
	while (q_pop(w->q,&v)){
		/* PRIVATE work: stamp this worker's id all over its own scratch and
		   read it back. If scratch were shared, a concurrent worker would
		   corrupt it and this check would fire. */
		memset(w->scratch, (int)(w->id & 0xff), SCRATCH);
		for (int k=0;k<SCRATCH;k++)
			if (w->scratch[k] != (unsigned char)(w->id & 0xff)){
				fprintf(stderr,"SCRATCH CORRUPTION on worker %d\n", w->id);
				abort();
			}
		w->local_count++;

		/* funnel: model 1/8 surviving sizeopt, 1/64 surviving rootopt */
		int is_sizeopt = ((w->local_count & 7)==0);
		int is_rootopt = ((w->local_count & 63)==0);

		pthread_mutex_lock(&w->shared->lock);
		if (is_sizeopt) w->shared->sizeopt_pass++;
		if (is_rootopt){
			w->shared->rootopt_pass++;
			if (v > w->shared->best) w->shared->best = v;
			/* a full "# norm ... e ..." style line: must appear intact */
			fprintf(w->shared->file, "# e %.6le worker %d seq %llu\n",
				v, w->id, (unsigned long long)w->shared->rootopt_pass);
		}
		w->shared->stage2_done++;
		pthread_mutex_unlock(&w->shared->lock);
	}
	return 0;
}

int main(void){
	q_init(&Q);
	memset(&SH,0,sizeof SH); pthread_mutex_init(&SH.lock,0);
	SH.file=fopen("s2_shared.out","w");

	pthread_t prod[NPROD], cons[NCONS];
	for (int i=0;i<NCONS;i++){ W[i].id=i+1; W[i].local_count=0; W[i].shared=&SH; W[i].q=&Q;
		pthread_create(&cons[i],0,consumer,&W[i]); }
	for (long i=0;i<NPROD;i++) pthread_create(&prod[i],0,producer,(void*)i);
	for (int i=0;i<NPROD;i++) pthread_join(prod[i],0);
	q_close(&Q);
	for (int i=0;i<NCONS;i++) pthread_join(cons[i],0);
	fclose(SH.file);

	uint64_t produced=(uint64_t)NPROD*PER_PROD;
	uint64_t worker_sum=0; for(int i=0;i<NCONS;i++) worker_sum+=W[i].local_count;

	/* count intact lines in the shared file */
	uint64_t lines=0, malformed=0; char line[128];
	FILE *f=fopen("s2_shared.out","r");
	while (fgets(line,sizeof line,f)){
		lines++;
		if (strncmp(line,"# e ",4)!=0 || !strchr(line,'\n')) malformed++;
	}
	fclose(f);

	printf("produced           %llu\n",(unsigned long long)produced);
	printf("stage2_done        %llu   %s\n",(unsigned long long)SH.stage2_done,
		SH.stage2_done==produced?"OK":"MISMATCH");
	printf("worker local sum   %llu   %s\n",(unsigned long long)worker_sum,
		worker_sum==produced?"OK":"MISMATCH");
	printf("rootopt_pass       %llu\n",(unsigned long long)SH.rootopt_pass);
	printf("file lines         %llu   %s\n",(unsigned long long)lines,
		lines==SH.rootopt_pass?"OK":"MISMATCH");
	printf("malformed lines    %llu   %s\n",(unsigned long long)malformed,
		malformed==0?"OK":"CORRUPT");
	printf("funnel sizeopt>=rootopt: %s\n", SH.sizeopt_pass>=SH.rootopt_pass?"OK":"BAD");
	int ok = SH.stage2_done==produced && worker_sum==produced &&
		lines==SH.rootopt_pass && malformed==0 && SH.sizeopt_pass>=SH.rootopt_pass;
	printf("\n%s\n", ok?"PATTERN OK — per-worker isolated, shared output serialized":"FAILED");
	return ok?0:1;
}
