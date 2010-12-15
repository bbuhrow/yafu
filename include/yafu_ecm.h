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

#if !defined(HAVE_GMP) || !defined(HAVE_GMP_ECM)
	typedef struct 
	{
		z X;
		z Z;
	} ecm_pt;

	typedef struct 
	{
		z sum1;
		z diff1;
		z sum2;
		z diff2;
		z tt1;
		z tt2;
		z tt3;
		z s;
		ecm_pt *A;
		ecm_pt *B;
		ecm_pt *C;
		ecm_pt *tmp1;
		ecm_pt *tmp2;
	} ecm_work;

	typedef struct {
		ecm_work *work;
		uint8 *marks;
		uint8 *nmarks;
		z n;
		z u;
		z v;
		z factor;
		z bb;
		z acc;
		z Paprod;
		z *Pbprod;		//vector of D big ints
		ecm_pt *P;
		ecm_pt *Pb;		//vector of D ecm_pts
		ecm_pt Pa;
		ecm_pt Pd;
		ecm_pt Pad;
		int D;
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

	void ecm_thread_free(ecm_thread_data_t *tdata);
	void ecm_thread_init(z *n, int D, ecm_thread_data_t *tdata);
	void ecm_work_init(ecm_work *work);
	void ecm_work_free(ecm_work *work);
	void ecm_pt_init(ecm_pt *pt);
	void ecm_pt_free(ecm_pt *pt);
	void add(ecm_pt *Pin, ecm_pt *Pout, ecm_work *work, z *n);
	void duplicate(z *insum, z *indiff, ecm_pt *P, ecm_work *work, z *n);
	void next_pt(ecm_pt *P, int c, ecm_work *work, z *n);
	void prac (ecm_pt *A, unsigned long k, z *b, ecm_work *work, z *n);
	void pracdup(ecm_pt *P, z *diff, z *sum, ecm_work *work, z *b, z *n);
	void pracadd(ecm_pt *P, ecm_pt *P0, ecm_work *work, z *n);
	void pt_diff(ecm_pt *P, z *diff, z *n);
	void pt_sum(ecm_pt *P, z *sum, z *n);
	void pt_copy(ecm_pt *src, ecm_pt *dest);
	void ecm_pt_clear(ecm_pt *pt);
	int check_factor(z *Z, z *n, z *f);
#else
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

#endif

//globals
int ECM_ABORT;
int PM1_ABORT;
int PP1_ABORT;

