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

#include "yafu.h"
#include "factor.h"
#include "ytools.h"
#include "arith.h"
#include <gmp_xface.h>
#include <ecm.h>

/* 
declarations for types and functions used in ECM and other
group theoretic factorization routines
*/

typedef struct {
	mpz_t gmp_n, gmp_factor;
    mpz_t t, d; // scratch bignums
	ecm_params params;
	uint32 sigma;
	int stagefound;
	fact_obj_t *fobj;
	int thread_num;
	int curves_run;
    int *total_curves_run;
	char tmp_output[80];
    int factor_found;
    int *ok_to_stop;
    int *curves_in_flight;

    // timing data for ETA
    struct timeval stop;
    struct timeval start;
    double *total_time;

} ecm_thread_data_t;

//local function declarations
void *ecm_do_one_curve(void *ptr);
int print_B1B2(fact_obj_t *fobj, FILE *fid);
void ecmexit(int sig);
void ecm_process_init(fact_obj_t *fobj);
void ecm_process_free(fact_obj_t *fobj);
int ecm_check_input(fact_obj_t *fobj);
int ecm_get_sigma(ecm_thread_data_t *thread_data);
int ecm_deal_with_factor(ecm_thread_data_t *thread_data);
void ecm_stop_worker_thread(ecm_thread_data_t *t, uint32 is_master_thread);
void ecm_start_worker_thread(ecm_thread_data_t *t, uint32 is_master_thread);
void ecm_thread_free(ecm_thread_data_t *tdata);
void ecm_thread_init(ecm_thread_data_t *tdata);

#if defined(WIN32) || defined(_WIN64)
DWORD WINAPI ecm_worker_thread_main(LPVOID thread_data);
#else
void *ecm_worker_thread_main(void *thread_data);
#endif

// "local" globals
int ECM_ABORT;
int PM1_ABORT;
int PP1_ABORT;

int TMP_THREADS;
uint64 TMP_STG2_MAX;
