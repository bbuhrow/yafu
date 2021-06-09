/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may
benefit from your work.

Some parts of the code (and also this header), included in this
distribution have been reused from other sources. In particular I
have benefitted greatly from the work of Jason Papadopoulos's msieve @
www.boo.net/~jasonp, Scott Contini's mpqs implementation, and Tom St.
Denis Tom's Fast Math library.  Many thanks to their kind donation of
code to the public domain.
                       --bbuhrow@gmail.com 11/24/09
----------------------------------------------------------------------*/

#ifndef _NFS_IMPL_H_
#define _NFS_IMPL_H_

#pragma once

#include "nfs.h"
#include "factor.h"
#include "gnfs.h"
#include "arith.h"
#include "qs.h"
#include <stdint.h>

/* used to place a deadline on how long polynomial
   selection will run. Note that the time budget is
   independent of CPU speed; faster CPUs will simply
   search more of the polynomial space */

typedef struct {
    uint32_t bits;
    uint32_t seconds;
} poly_deadline_t;

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

typedef struct
{
    mpz_polys_t* poly; // the idea is that job->snfs->poly == job->poly
    uint32_t rlim, alim;
    uint32_t lpbr, lpba;
    uint32_t mfbr, mfba;
    double rlambda, alambda;
    uint32_t qrange; // how large a sieving block is
    char sievername[1024];
    uint32_t startq;
    uint32_t min_rels;
    uint32_t current_rels;
    uint32_t poly_time;
    uint32_t last_leading_coeff;
    uint32_t use_max_rels;
    double test_score;      // used to remember score if test-sieved.

    snfs_t* snfs; // NULL if GNFS
} nfs_job_t;

typedef struct {
    // stuff for parallel ggnfs sieving
    char outfilename[80];
    nfs_job_t job;
    uint32_t siever;

    // stuff for parallel msieve poly select
    char* polyfilename, * logfilename, * fbfilename;
    uint64_t poly_lower;
    uint64_t poly_upper;
    msieve_obj* obj;
    mp_t* mpN;
    factor_list_t* factor_list;
    struct timeval thread_start_time;
    fact_obj_t* fobj;

    int tindex;
    int is_poly_select;

    /* fields for thread pool synchronization */
    volatile enum nfs_thread_command command;
    volatile int* thread_queue, * threads_waiting;

#if defined(WIN32) || defined(_WIN64)
    HANDLE thread_id;
    HANDLE run_event;
    HANDLE finish_event;

    HANDLE* queue_event;
    HANDLE* queue_lock;
#else
    pthread_t thread_id;
    pthread_mutex_t run_lock;
    pthread_cond_t run_cond;

    pthread_mutex_t* queue_lock;
    pthread_cond_t* queue_cond;
#endif

} nfs_threaddata_t; 


//----------------------- NFS FUNCTIONS -------------------------------------//
void* lasieve_launcher(void* ptr);
void* polyfind_launcher(void* ptr);
double find_best_msieve_poly(fact_obj_t* fobj, nfs_job_t* job, int write_jobfile);
void msieve_to_ggnfs(fact_obj_t* fobj, nfs_job_t* job);
void ggnfs_to_msieve(fact_obj_t* fobj, nfs_job_t* job);
void get_ggnfs_params(fact_obj_t* fobj, nfs_job_t* job);
int check_for_sievers(fact_obj_t* fobj, int revert_to_siqs);
void print_poly(mpz_polys_t* poly, FILE* out);
void print_job(nfs_job_t* job, FILE* out);
uint32_t parse_job_file(fact_obj_t* fobj, nfs_job_t* job);
void fill_job_file(fact_obj_t* fobj, nfs_job_t* job, uint32_t missing_params);

enum nfs_state_e check_existing_files(fact_obj_t* fobj, uint32_t* last_spq, nfs_job_t* job);
void extract_factors(factor_list_t* factor_list, fact_obj_t* fobj);
uint32_t get_spq(char** lines, int last_line, fact_obj_t* fobj);
uint32_t do_msieve_filtering(fact_obj_t* fobj, msieve_obj* obj, nfs_job_t* job);
void do_msieve_polyselect(fact_obj_t* fobj, msieve_obj* obj, nfs_job_t* job, mp_t* mpN, factor_list_t* factor_list);
void get_polysearch_params(fact_obj_t* fobj, uint64_t* start, uint64_t* range);
void init_poly_threaddata(nfs_threaddata_t* t, msieve_obj* obj,
    mp_t* mpN, factor_list_t* factor_list, int tid, uint32_t flags, uint64_t start, uint64_t stop);
void do_sieving(fact_obj_t* fobj, nfs_job_t* job);
void trial_sieve(fact_obj_t* fobj); // external test sieve frontend
int test_sieve(fact_obj_t* fobj, void* args, int njobs, int are_files);
void savefile_concat(char* filein, char* fileout, msieve_obj* mobj);
void win_file_concat(char* filein, char* fileout);
void nfs_stop_worker_thread(nfs_threaddata_t* t,
    uint32_t is_master_thread);
void nfs_start_worker_thread(nfs_threaddata_t* t,
    uint32_t is_master_thread);
void nfsexit(int sig);
#if defined(WIN32) || defined(_WIN64)
DWORD WINAPI nfs_worker_thread_main(LPVOID thread_data);
#else
void* nfs_worker_thread_main(void* thread_data);
#endif


//----------------------- SNFS FUNCTIONS -------------------------------------//
void find_brent_form(fact_obj_t* fobj, snfs_t* poly);
void find_hcunn_form(fact_obj_t* fobj, snfs_t* poly);
void find_xyyxf_form(fact_obj_t* fobj, snfs_t* poly);
void find_direct_form(fact_obj_t* fobj, snfs_t* poly);
snfs_t* gen_brent_poly(fact_obj_t* fobj, snfs_t* poly, int* npolys); // the workhorse
snfs_t* gen_xyyxf_poly(fact_obj_t* fobj, snfs_t* poly, int* npolys);
int snfs_choose_poly(fact_obj_t* fobj, nfs_job_t* job);
void check_poly(snfs_t* poly, int VFLAG);
void compute_difficulty_from_poly(snfs_t* poly, int VFLAG);
void print_snfs(snfs_t* poly, FILE* out);
void snfs_copy_poly(snfs_t* src, snfs_t* dest);
void approx_norms(snfs_t* poly);
void snfs_scale_difficulty(snfs_t* polys, int npoly, int VFLAG);
int snfs_rank_polys(fact_obj_t* fobj, snfs_t* polys, int npoly);
int qcomp_snfs_sdifficulty(const void* x, const void* y);
int qcomp_snfs_murphy(const void* x, const void* y);
nfs_job_t* snfs_test_sieve(fact_obj_t* fobj, snfs_t* polys, int npoly, nfs_job_t* jobs, int force_test);
void snfs_make_job_file(fact_obj_t* fobj, nfs_job_t* job);
void snfs_init(snfs_t* poly);
void snfs_clear(snfs_t* poly);
void skew_snfs_params(fact_obj_t* fobj, nfs_job_t* job);
void nfs_set_min_rels(nfs_job_t* job);
void copy_job(nfs_job_t* src, nfs_job_t* dest);
void copy_mpz_polys_t(mpz_polys_t* src, mpz_polys_t* dest);
void analyze_one_poly_xface(snfs_t* poly);
int est_gnfs_size(nfs_job_t* job);
int est_gnfs_size_via_poly(snfs_t* job);
int tdiv_mpz(mpz_t x, int* factors, uint64_t *primes, uint64_t num_p);

void mpz_polys_init(mpz_polys_t* poly);
void mpz_polys_free(mpz_polys_t* poly);


extern int NFS_ABORT;
extern int IGNORE_NFS_ABORT;
#define GGNFS_TABLE_ROWS 23
extern double ggnfs_table[GGNFS_TABLE_ROWS][8];



#endif


