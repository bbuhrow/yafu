#include "gnfs.h"
#include "yafu_string.h"
#include "arith.h"
#include "factor.h"
#include "qs.h"

//----------------------- LOCAL DATA TYPES -----------------------------------//

/* used to place a deadline on how long polynomial 
   selection will run. Note that the time budget is
   independent of CPU speed; faster CPUs will simply
   search more of the polynomial space */

typedef struct {
	uint32 bits;
	uint32 seconds;
} poly_deadline_t;

static const poly_deadline_t time_limits[] = {
	{264, 4 * 60},
	{304, 8 * 60},
	{320, 15 * 60},
	{348, 30 * 60},
	{365, 1 * 3600},
	{383, 2 * 3600},
	{399, 4 * 3600},
	{416, 8 * 3600},
	{433, 16 * 3600},
	{449, 32 * 3600},
	{466, 64 * 3600},
	{482, 100 * 3600},
	{498, 200 * 3600},
	{514, 300 * 3600},
};

#define NUM_TIME_LIMITS sizeof(time_limits)/sizeof(time_limits[0])

enum nfs_thread_command {
	NFS_COMMAND_INIT,
	NFS_COMMAND_WAIT,
	NFS_COMMAND_RUN,
	NFS_COMMAND_RUN_POLY,
	NFS_COMMAND_END
};

enum nfs_state_e
{
	NFS_STATE_INIT,
	NFS_STATE_POLY,
	NFS_STATE_SIEVE,
	NFS_STATE_FILTER,
	NFS_STATE_LINALG,
	NFS_STATE_SQRT,
	NFS_STATE_CLEANUP,
	NFS_STATE_FILTCHECK,
	NFS_STATE_STARTNEW,
	NFS_STATE_RESUMESIEVE,
	NFS_STATE_RESUMEPOLY,
	NFS_STATE_DONE
};

enum param_flag_e {
	PARAM_FLAG_NONE = 0,
	PARAM_FLAG_FBLIM = 0x1,
	PARAM_FLAG_LPB = 0x2,
	PARAM_FLAG_MFB = 0x4,
	PARAM_FLAG_LAMBDA = 0x8
};

typedef struct
{
	uint32 fblim;
	uint32 lpb;
	uint32 mfb;
	float lambda;
	uint32 siever; 
	uint32 qrange;
	char sievername[1024];
	uint32 startq;
	uint32 min_rels;
	uint32 current_rels;
	uint32 type;			// 0==GNFS, 1==SNFS
	uint32 size;			// used for SNFS difficulty
	uint32 poly_time;
	uint32 last_leading_coeff;
} ggnfs_job_t;

typedef struct {
	// stuff for parallel ggnfs sieving
	char outfilename[80];
	ggnfs_job_t job;

	// stuff for parallel msieve poly select
	char *polyfilename, *logfilename, *fbfilename;
	uint64 poly_lower;
	uint64 poly_upper;
	msieve_obj *obj;
	mp_t *mpN;
	factor_list_t *factor_list;
	struct timeval thread_start_time;
	fact_obj_t *fobj;

	int tindex;
	int is_poly_select;

	/* fields for thread pool synchronization */
	volatile enum nfs_thread_command command;
	volatile int *thread_queue, *threads_waiting;

#if defined(WIN32) || defined(_WIN64)
	HANDLE thread_id;
	HANDLE run_event;
	HANDLE finish_event;

	HANDLE *queue_event;
	HANDLE *queue_lock;
#else
	pthread_t thread_id;
	pthread_mutex_t run_lock;
	pthread_cond_t run_cond;

	pthread_mutex_t *queue_lock;
	pthread_cond_t *queue_cond;
#endif

} nfs_threaddata_t;

//----------------------- LOCAL FUNCTIONS -------------------------------------//
void *lasieve_launcher(void *ptr);
void *polyfind_launcher(void *ptr);
void find_best_msieve_poly(fact_obj_t *fobj, ggnfs_job_t *job, int write_jobfile);
void msieve_to_ggnfs(fact_obj_t *fobj, ggnfs_job_t *job);
void ggnfs_to_msieve(fact_obj_t *fobj, ggnfs_job_t *job);
void get_ggnfs_params(fact_obj_t *fobj, ggnfs_job_t *job);
int check_for_sievers(fact_obj_t *fobj);
//void parse_job_file(fact_obj_t *fobj, ggnfs_job_t *job);
void parse_job_file(fact_obj_t *fobj, ggnfs_job_t *job, uint32* missing_params);
void fill_job_file(fact_obj_t *fobj, ggnfs_job_t *job, uint32 missing_params);

enum nfs_state_e check_existing_files(fact_obj_t *fobj, uint32 *last_spq, ggnfs_job_t *job);
void extract_factors(factor_list_t *factor_list, fact_obj_t *fobj);
uint32 get_spq(char **lines, int last_line, fact_obj_t *fobj);
uint32 do_msieve_filtering(fact_obj_t *fobj, msieve_obj *obj, ggnfs_job_t *job);
void do_msieve_polyselect(fact_obj_t *fobj, msieve_obj *obj, ggnfs_job_t *job, mp_t *mpN, factor_list_t *factor_list);
void get_polysearch_params(fact_obj_t *fobj, uint64 *start, uint64 *range);
void init_poly_threaddata(nfs_threaddata_t *t, msieve_obj *obj, 
	mp_t *mpN, factor_list_t *factor_list, int tid, uint32 flags, uint64 start, uint64 stop);
void do_sieving(fact_obj_t *fobj, ggnfs_job_t *job);
void savefile_concat(char *filein, char *fileout, msieve_obj *mobj);
void win_file_concat(char *filein, char *fileout);
void nfs_stop_worker_thread(nfs_threaddata_t *t,
				uint32 is_master_thread);
void nfs_start_worker_thread(nfs_threaddata_t *t, 
				uint32 is_master_thread);
void nfsexit(int sig);
#if defined(WIN32) || defined(_WIN64)
DWORD WINAPI nfs_worker_thread_main(LPVOID thread_data);
#else
void *nfs_worker_thread_main(void *thread_data);
#endif

int NFS_ABORT;

