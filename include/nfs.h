#include "gnfs.h"
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
//  bits, seconds
    {248, 1 * 60},		// 74 digits
    {264, 2 * 60},		// 80 digits
	{304, 6 * 60},		// 92 digits
	{320, 15 * 60},		// 97 digits
	{348, 30 * 60},		// 105 digits
	{365, 1 * 3600},	// 110 digits
	{383, 2 * 3600},	// 116 digits
	{399, 4 * 3600},	// 120 digits
	{416, 8 * 3600},	// 126 digits
	{433, 16 * 3600},	// 131 digits
	{449, 32 * 3600},	// 135 digits
	{466, 64 * 3600},	// 140 digits
	{482, 100 * 3600},	// 146 digits
	{498, 200 * 3600},	// 150 digits
	{514, 300 * 3600},	// 155 digits
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

enum param_flag_e
{
	PARAM_FLAG_NONE = 0,

	PARAM_FLAG_RLIM = 0x1,
	PARAM_FLAG_ALIM = 0x2,
	PARAM_FLAG_FBLIM = 0x1 + 0x2,

	PARAM_FLAG_LPBR = 0x4,
	PARAM_FLAG_LPBA = 0x8,
	PARAM_FLAG_LPB = 0x4 + 0x8,

	PARAM_FLAG_MFBR = 0x10,
	PARAM_FLAG_MFBA = 0x20,
	PARAM_FLAG_MFB = 0x10 + 0x20,

	PARAM_FLAG_RLAMBDA = 0x40,
	PARAM_FLAG_ALAMBDA = 0x80,
	PARAM_FLAG_LAMBDA = 0x40 + 0x80,

	PARAM_FLAG_ALL = 0xFF
};

enum snfs_form_e
{
	SNFS_NONE,
	SNFS_BRENT,
	SNFS_H_CUNNINGHAM,
	SNFS_XYYXF
};

enum special_q_e
{
	NEITHER_SPQ,
	RATIONAL_SPQ,
	ALGEBRAIC_SPQ
};

// gnfs.h has an mpz_poly_t struct:
/* typedef struct {
	uint32 degree;
	mpz_t coeff[MAX_POLY_DEGREE + 1];

	// scratch quantities for evaluating the homogeneous form of poly 
	mpz_t tmp1, tmp2, tmp3;
} mpz_poly_t; */

typedef struct
{
	mpz_poly_t rat; // linear (usually)
	mpz_poly_t alg;
	double skew;
	double murphy; // murphy e score
	double size;
	double alpha;
	int rroots;
	mpz_t m; // common root mod n
	enum special_q_e side;
} mpz_polys_t;

#define NUM_SNFS_POLYS 3
#define MAX_SNFS_BITS 1024

typedef struct
{
	// input integer
	mpz_t n; // the cofactor, what goes in the job file
	mpz_t primitive;
	// algebraic representation of the snfs form:
	// n divides c1*b1^e1 + c2*b2^e2
	mpz_t base1;
	mpz_t base2;
	int exp1;
	int exp2;
	int coeff1;
	int coeff2;
	// type of form
	enum snfs_form_e form_type;


	mpz_polys_t* poly;
	mpz_t c[MAX_POLY_DEGREE + 1]; // scratch space -- converted to mpz_poly_t
				      // in check_poly()

	// other useful parameters
	double difficulty;
	double sdifficulty;
	double anorm;
	double rnorm;
	int rank;
	int valid;
	int siever;
} snfs_t;

typedef struct
{
	mpz_polys_t* poly; // the idea is that job->snfs->poly == job->poly
	uint32 rlim, alim;
	uint32 lpbr, lpba;
	uint32 mfbr, mfba;
	double rlambda, alambda;
	uint32 qrange; // how large a sieving block is
	char sievername[1024];
	uint32 startq;
	uint32 min_rels;
	uint32 current_rels;
	uint32 poly_time;
	uint32 last_leading_coeff;
	uint32 use_max_rels;

	snfs_t* snfs; // NULL if GNFS
} nfs_job_t;

typedef struct {
	// stuff for parallel ggnfs sieving
	char outfilename[80];
	nfs_job_t job;
	uint32 siever;

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
double find_best_msieve_poly(fact_obj_t *fobj, nfs_job_t *job, int write_jobfile);
void msieve_to_ggnfs(fact_obj_t *fobj, nfs_job_t *job);
void ggnfs_to_msieve(fact_obj_t *fobj, nfs_job_t *job);
void get_ggnfs_params(fact_obj_t *fobj, nfs_job_t *job);
int check_for_sievers(fact_obj_t *fobj, int revert_to_siqs);
void print_poly(mpz_polys_t* poly, FILE *out);
void print_job(nfs_job_t *job, FILE *out);
uint32 parse_job_file(fact_obj_t *fobj, nfs_job_t *job);
void fill_job_file(fact_obj_t *fobj, nfs_job_t *job, uint32 missing_params);

enum nfs_state_e check_existing_files(fact_obj_t *fobj, uint32 *last_spq, nfs_job_t *job);
void extract_factors(factor_list_t *factor_list, fact_obj_t *fobj);
uint32 get_spq(char **lines, int last_line, fact_obj_t *fobj);
uint32 do_msieve_filtering(fact_obj_t *fobj, msieve_obj *obj, nfs_job_t *job);
void do_msieve_polyselect(fact_obj_t *fobj, msieve_obj *obj, nfs_job_t *job, mp_t *mpN, factor_list_t *factor_list);
void get_polysearch_params(fact_obj_t *fobj, uint64 *start, uint64 *range);
void init_poly_threaddata(nfs_threaddata_t *t, msieve_obj *obj, 
	mp_t *mpN, factor_list_t *factor_list, int tid, uint32 flags, uint64 start, uint64 stop);
void do_sieving(fact_obj_t *fobj, nfs_job_t *job);
void trial_sieve(fact_obj_t* fobj); // external test sieve frontend
int test_sieve(fact_obj_t* fobj, void* args, int njobs, int are_files);
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

// snfs stuff
void find_brent_form(fact_obj_t *fobj, snfs_t *poly);
void find_hcunn_form(fact_obj_t *fobj, snfs_t *poly);
void find_xyyxf_form(fact_obj_t *fobj, snfs_t *poly);
snfs_t* gen_brent_poly(fact_obj_t *fobj, snfs_t *poly, int* npolys); // the workhorse
snfs_t* gen_xyyxf_poly(fact_obj_t *fobj, snfs_t *poly, int* npolys);
int snfs_choose_poly(fact_obj_t* fobj, nfs_job_t* job);
void check_poly(snfs_t *poly, int VFLAG);
void print_snfs(snfs_t *poly, FILE *out);
void snfs_copy_poly(snfs_t *src, snfs_t *dest);
void approx_norms(snfs_t *poly);
void snfs_scale_difficulty(snfs_t *polys, int npoly, int VFLAG);
int snfs_rank_polys(fact_obj_t *fobj, snfs_t *polys, int npoly);
int qcomp_snfs_sdifficulty(const void *x, const void *y);
int qcomp_snfs_murphy(const void *x, const void *y);
nfs_job_t *snfs_test_sieve(fact_obj_t *fobj, snfs_t *polys, int npoly, nfs_job_t *jobs);
void snfs_make_job_file(fact_obj_t *fobj, nfs_job_t *job);
void snfs_init(snfs_t* poly);
void snfs_clear(snfs_t* poly);
void skew_snfs_params(fact_obj_t *fobj, nfs_job_t *job);
void find_primitive_factor(snfs_t *poly, int VFLAG);
void nfs_set_min_rels(nfs_job_t *job);
void copy_job(nfs_job_t *src, nfs_job_t *dest);
void copy_mpz_polys_t(mpz_polys_t *src, mpz_polys_t *dest);
void analyze_one_poly_xface(snfs_t *poly);
int est_gnfs_size(nfs_job_t *job);
int est_gnfs_size_via_poly(snfs_t *job);
snfs_t * snfs_find_form(fact_obj_t *fobj);
int tdiv_mpz(mpz_t x, int *factors);


int NFS_ABORT;
int IGNORE_NFS_ABORT;

static INLINE void mpz_polys_init(mpz_polys_t * poly) {
	mpz_poly_init(&poly->rat);
	mpz_poly_init(&poly->alg);
	mpz_init(poly->m);
	poly->skew = 0;
	poly->murphy = 0.;
	poly->size = 0.;
	poly->rroots = 0;
	poly->side = ALGEBRAIC_SPQ; // snfs routines will override if necessary
}

static INLINE void mpz_polys_free(mpz_polys_t * poly) {
	mpz_poly_free(&poly->rat);
	mpz_poly_free(&poly->alg);
	mpz_clear(poly->m);
}
