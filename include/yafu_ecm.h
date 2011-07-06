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
#include "arith.h"

/* 
declarations for types and functions used in ECM and other
group theoretic factorization routines
*/

//brent's rho
void brent_loop(fact_obj_t *fobj);

//pollard's p-1
void pollard_loop(fact_obj_t *fobj);

//williams p+1
void williams_loop(int trials, fact_obj_t *fobj);

//ecm
int ecm_loop(z *n, int numcurves, fact_obj_t *fobj);

//state save/recover for p-1,p+1
void recover_stg1(z *n, z *m, int method);
void recover_stg2(z *n, z *cc, z *m, z *mm, int method);
void save_state(int stage, z *n, z *m, z *mm, z *cc, int i, int method);

enum ecm_thread_command {
	ECM_COMMAND_INIT,
	ECM_COMMAND_WAIT,
	ECM_COMMAND_RUN,
	ECM_COMMAND_RUN_TRANS,
	ECM_COMMAND_END
};


#include <gmp_xface.h>
#include <ecm.h>

typedef struct {
	mpz_t gmp_n, gmp_factor;
	z n, factor;
	ecm_params params;
	uint32 sigma;
	int stagefound;

	/* fields for thread pool synchronization */
	volatile enum ecm_thread_command command;

#if defined(WIN32) || defined(_WIN64)
	HANDLE thread_id;
	HANDLE run_event;
	HANDLE finish_event;
#else
	pthread_t thread_id;
	pthread_mutex_t run_lock;
	pthread_cond_t run_cond;
#endif

} ecm_thread_data_t;


//globals
int ECM_ABORT;
int PM1_ABORT;
int PP1_ABORT;

