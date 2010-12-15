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
       				   --bbuhrow@gmail.com 11/2/10
----------------------------------------------------------------------*/

#include "yafu.h"
#include "yafu_ecm.h"
#include "soe.h"
#include "factor.h"
#include "monty.h"
#include "util.h"
#include "calc.h"

//local function declarations
void *ecm_do_one_curve(void *ptr);
int print_B1B2(void);
void ecmexit(int sig);
void ecm_process_init(z *n);
void ecm_process_free();

void ecm_stop_worker_thread(ecm_thread_data_t *t,
				uint32 is_master_thread);
void ecm_start_worker_thread(ecm_thread_data_t *t, 
				uint32 is_master_thread);
#if defined(WIN32) || defined(_WIN64)
DWORD WINAPI ecm_worker_thread_main(LPVOID thread_data);
#else
void *ecm_worker_thread_main(void *thread_data);
#endif

void ecm_start_worker_thread(ecm_thread_data_t *t, 
				uint32 is_master_thread) {

	/* create a thread that will process a polynomial The last poly does 
	   not get its own thread (the current thread handles it) */

	if (is_master_thread) {
		//matrix_thread_init(t);
		return;
	}

	t->command = ECM_COMMAND_INIT;
#if defined(WIN32) || defined(_WIN64)
	t->run_event = CreateEvent(NULL, FALSE, TRUE, NULL);
	t->finish_event = CreateEvent(NULL, FALSE, FALSE, NULL);
	t->thread_id = CreateThread(NULL, 0, ecm_worker_thread_main, t, 0, NULL);

	WaitForSingleObject(t->finish_event, INFINITE); /* wait for ready */
#else
	pthread_mutex_init(&t->run_lock, NULL);
	pthread_cond_init(&t->run_cond, NULL);

	pthread_cond_signal(&t->run_cond);
	pthread_mutex_unlock(&t->run_lock);
	pthread_create(&t->thread_id, NULL, ecm_worker_thread_main, t);

	pthread_mutex_lock(&t->run_lock); /* wait for ready */
	while (t->command != ECM_COMMAND_WAIT)
		pthread_cond_wait(&t->run_cond, &t->run_lock);
#endif
}

void ecm_stop_worker_thread(ecm_thread_data_t *t,
				uint32 is_master_thread)
{
	if (is_master_thread) {
		//matrix_thread_free(t);
		return;
	}

	t->command = ECM_COMMAND_END;
#if defined(WIN32) || defined(_WIN64)
	SetEvent(t->run_event);
	WaitForSingleObject(t->thread_id, INFINITE);
	CloseHandle(t->thread_id);
	CloseHandle(t->run_event);
	CloseHandle(t->finish_event);
#else
	pthread_cond_signal(&t->run_cond);
	pthread_mutex_unlock(&t->run_lock);
	pthread_join(t->thread_id, NULL);
	pthread_cond_destroy(&t->run_cond);
	pthread_mutex_destroy(&t->run_lock);
#endif
}

#if defined(WIN32) || defined(_WIN64)
DWORD WINAPI ecm_worker_thread_main(LPVOID thread_data) {
#else
void *ecm_worker_thread_main(void *thread_data) {
#endif
	ecm_thread_data_t *t = (ecm_thread_data_t *)thread_data;

	while(1) {

		/* wait forever for work to do */
#if defined(WIN32) || defined(_WIN64)
		WaitForSingleObject(t->run_event, INFINITE);		
#else
		pthread_mutex_lock(&t->run_lock);
		while (t->command == ECM_COMMAND_WAIT) {
			pthread_cond_wait(&t->run_cond, &t->run_lock);
		}
#endif
		/* do work */

		if (t->command == ECM_COMMAND_RUN)
			ecm_do_one_curve(t);
		//else if (t->command == COMMAND_RUN_TRANS)
		//	mul_trans_packed_core(t);
		//else if (t->command == COMMAND_INIT)
		//	matrix_thread_init(t);
		else if (t->command == ECM_COMMAND_END)
			break;

		/* signal completion */

		t->command = ECM_COMMAND_WAIT;
#if defined(WIN32) || defined(_WIN64)
		SetEvent(t->finish_event);		
#else
		pthread_cond_signal(&t->run_cond);
		pthread_mutex_unlock(&t->run_lock);
#endif
	}

	//matrix_thread_free(t);

#if defined(WIN32) || defined(_WIN64)
	return 0;
#else
	return NULL;
#endif
}


// declarations/definitions when using YAFU's ECM
#if !defined(HAVE_GMP) || !defined(HAVE_GMP_ECM)

	void ecm_work_init(ecm_work *work)
	{
		zInit(&work->diff1);
		zInit(&work->diff2);
		zInit(&work->sum1);
		zInit(&work->sum2);
		zInit(&work->tt1);
		zInit(&work->tt2);
		zInit(&work->tt3);
		zInit(&work->s);
		work->A = (ecm_pt *)malloc(sizeof(ecm_pt));
		work->B = (ecm_pt *)malloc(sizeof(ecm_pt));
		work->C = (ecm_pt *)malloc(sizeof(ecm_pt));
		work->tmp1 = (ecm_pt *)malloc(sizeof(ecm_pt));
		work->tmp2 = (ecm_pt *)malloc(sizeof(ecm_pt));
		ecm_pt_init(work->A);
		ecm_pt_init(work->B);
		ecm_pt_init(work->C);
		ecm_pt_init(work->tmp1);
		ecm_pt_init(work->tmp2);

		return;
	}

	void ecm_work_free(ecm_work *work)
	{
		zFree(&work->diff1);
		zFree(&work->diff2);
		zFree(&work->sum1);
		zFree(&work->sum2);
		zFree(&work->tt1);
		zFree(&work->tt2);
		zFree(&work->tt3);
		zFree(&work->s);
		ecm_pt_free(work->A);
		ecm_pt_free(work->B);
		ecm_pt_free(work->C);
		ecm_pt_free(work->tmp1);
		ecm_pt_free(work->tmp2);
		free(work->A);
		free(work->B);
		free(work->C);
		free(work->tmp1);
		free(work->tmp2);

		return;
	}

	void ecm_process_init(z *n)
	{
		//initialize things which all threads will need.

		//initialize montgomery arithmatic for this modulus, and cache needed primes using
		//the sieve of erathostenes
		GetPRIMESRange(0,10001000);
		if (ECM_STG2_MAX > PRIMES[NUM_P-1])
			GetPRIMESRange(0,ECM_STG2_MAX+1000);
		monty_init(n);

		return;
	}

	void ecm_process_free()
	{
		monty_free();
		return;
	}

	void ecm_thread_init(z *n, int D, ecm_thread_data_t *tdata)
	{
		//initialize everything for all threads doing ecm
		int j;

		zInit(&tdata->u);
		zInit(&tdata->v);
		zInit(&tdata->bb);
		zInit(&tdata->acc);
		zInit(&tdata->Paprod);
		zInit(&tdata->n);
		zInit(&tdata->factor);

		zCopy(n,&tdata->n);
		tdata->D = D;

		tdata->work = (ecm_work *)malloc(sizeof(ecm_work));
		ecm_work_init(tdata->work);
		tdata->P = (ecm_pt *)malloc(sizeof(ecm_pt));
		ecm_pt_init(tdata->P);
		ecm_pt_init(&tdata->Pa);
		ecm_pt_init(&tdata->Pd);
		ecm_pt_init(&tdata->Pad);

		//build an array to hold values of f(b)
		tdata->Pb = (ecm_pt *)malloc(tdata->D * sizeof(ecm_pt));
		for (j=0;j<tdata->D;j++)
			ecm_pt_init(&tdata->Pb[j]);
		
		//for storage of products
		tdata->Pbprod = (z *)malloc(tdata->D * sizeof(z));
		for (j=0;j<tdata->D;j++)
			zInit(&tdata->Pbprod[j]);

		tdata->marks = (uint8 *)malloc(tdata->D * sizeof(uint8));
		tdata->nmarks = (uint8 *)malloc(tdata->D * sizeof(uint8));

		return;
	}

	void ecm_thread_free(ecm_thread_data_t *tdata)
	{
		int j;

		zFree(&tdata->u);
		zFree(&tdata->v);
		zFree(&tdata->bb);
		zFree(&tdata->acc);
		zFree(&tdata->Paprod);
		zFree(&tdata->n);
		zFree(&tdata->factor);

		ecm_work_free(tdata->work);
		free(tdata->work);

		ecm_pt_free(tdata->P);
		free(tdata->P);

		ecm_pt_free(&tdata->Pa);
		ecm_pt_free(&tdata->Pd);
		ecm_pt_free(&tdata->Pad);

		//build an array to hold values of f(b)
		for (j=0;j<tdata->D;j++)
			ecm_pt_free(&tdata->Pb[j]);
		free(tdata->Pb);
		
		//for storage of products
		for (j=0;j<tdata->D;j++)
			zFree(&tdata->Pbprod[j]);
		free(tdata->Pbprod);

		free(tdata->marks);
		free(tdata->nmarks);

		return;
	}

	void ecm_pt_init(ecm_pt *pt)
	{
		zInit(&pt->X);
		zInit(&pt->Z);
	}

	void ecm_pt_clear(ecm_pt *pt)
	{
		zClear(&pt->X);
		zClear(&pt->Z);
	}

	void ecm_pt_free(ecm_pt *pt)
	{
		zFree(&pt->X);
		zFree(&pt->Z);
	}

	int check_factor(z *Z, z *n, z *f)
	{
		zLEGCD(Z,n,f);
		if (zCompare(f,&zOne) > 0)
		{
			if (zCompare(f,n) == 0)
			{
				zCopy(&zZero,f);
				return 0;
			}
			return 1;
		}
		return 0;
	}

	void add(ecm_pt *Pin, ecm_pt *Pout, ecm_work *work, z *n)
	{
		//x+ = z- * [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
		//z+ = x- * [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
		//x- = original x
		//z- = original z

		monty_mul(&work->diff1,&work->sum2,&work->tt1,n);	//U
		monty_mul(&work->sum1,&work->diff2,&work->tt2,n);	//V
		
		monty_add(&work->tt1,&work->tt2,&Pout->X,n);		//U + V
		monty_sub(&work->tt1,&work->tt2,&Pout->Z,n);		//U - V
		monty_sqr(&Pout->X,&work->tt1,n);					//(U + V)^2
		monty_sqr(&Pout->Z,&work->tt2,n);					//(U - V)^2
		monty_mul(&work->tt1,&Pin->Z,&Pout->X,n);			//Z * (U + V)^2
		monty_mul(&work->tt2,&Pin->X,&Pout->Z,n);			//x * (U - V)^2

		return;
	}

	void duplicate(z *insum, z *indiff, ecm_pt *P, ecm_work *work, z *n)
	{
		monty_sqr(indiff,&work->tt1,n);				//(x1 - z1)^2
		monty_sqr(insum,&work->tt2,n);					//(x1 + z1)^2
		monty_mul(&work->tt1,&work->tt2,&P->X,n);

		monty_sub(&work->tt2,&work->tt1,&work->tt3,n);	//w1 = T
		monty_mul(&work->tt3,&work->s,&work->tt2,n);
		monty_add(&work->tt2,&work->tt1,&work->tt2,n);
		monty_mul(&work->tt2,&work->tt3,&P->Z,n);
		
		return;
	}

	void next_pt(ecm_pt *P, int c, ecm_work *work, z *n)
	{
		uint8 t,ci[32];
		int i;
		z *x1,*z1,*x2,*z2,*s1,*s2,*d1,*d2;

		x1 = &work->A->X;
		z1 = &work->A->Z;
		x2 = &work->B->X;
		z2 = &work->B->Z;
		s1 = &work->sum1;
		s2 = &work->sum2;
		d1 = &work->diff1;
		d2 = &work->diff2;

		//get c in binary form
		i=0;
		while (c>0)
		{
			i++;
			ci[i] = c % 2;
			c = c/2;
		}
		t=i;
		
		//goal is to compute x_c, z_c using montgomery's addition
		//and multiplication formula's for x and z
		//the procedure will be similar to p+1 lucas chain formula's
		//but this time we are simultaneously incrementing both x and z
		//rather than just Vn.  In each bit of the binary expansion of
		//c, we need just one addition and one duplication for both x and z.

		//initialize
		zCopy(&P->X,x1);
		zCopy(&P->Z,z1);
		monty_sub(&P->X,&P->Z,d1,n);
		monty_add(&P->X,&P->Z,s1,n);
		duplicate(s1,d1,work->B,work,n);

		//compute loop
		//for each bit of M to the right of the most significant bit
		for(i=t-1;i>=1;i--)
		{
			monty_sub(x1,z1,d1,n);
			monty_add(x1,z1,s1,n);
			monty_sub(x2,z2,d2,n);
			monty_add(x2,z2,s2,n);

			//if the bit is 1
			if (ci[i])
			{
				//add x1,z1, duplicate x2,z2
				add(P,work->A,work,n);
				duplicate(s2,d2,work->B,work,n);
			}
			else
			{
				//add x2,z2, duplicate x1,z1
				add(P,work->B,work,n);
				duplicate(s1,d1,work->A,work,n);
			}
		}
		zCopy(x1,&P->X);
		zCopy(z1,&P->Z);

		return;
	}

	void pt_diff(ecm_pt *P, z *diff, z *n)
	{
		//compute diff = P_x - P_z
		monty_sub(&P->X,&P->Z,diff, n);
		return;
	}
		
	void pt_sum(ecm_pt *P, z *sum, z *n)
	{
		//compute sum = P_x + P_z
		monty_add(&P->X,&P->Z,sum,n);
		return;
	}

	void pt_copy(ecm_pt *src, ecm_pt *dest)
	{
		//copy src to dest
		zCopy(&src->X,&dest->X);
		zCopy(&src->Z,&dest->Z);
		return;
	}

	void pracadd(ecm_pt *P, ecm_pt *P0, ecm_work *work, z *n)
	{
		//add two points (x1,z1) and (x2,z2) with respect to the "origin" P0, 
		//where (x1 + z1), (z1 - z1), (x2 - z2), and (x2 + z2) have been precomputed 
		//and stored in work.  output result in P. 
		
		monty_mul(&work->diff2,&work->sum1,&work->tt1,n);		// (x2-z2)*(x1+z1) 
		monty_mul(&work->sum2,&work->diff1,&work->tt2,n);		// (x2+z2)*(x1-z1) 
		monty_add(&work->tt1,&work->tt2,&work->tt3,n);		// 2*(x1*x2-z1*z2) 
		monty_sub(&work->tt1,&work->tt2,&work->tt2,n);		// 2*(x2*z1-x1*z2) 
		monty_sqr(&work->tt3,&work->tt3,n);					// 4*(x1*x2-z1*z2)^2 
		monty_sqr(&work->tt2,&work->tt2,n);					// 4*(x2*z1-x1*z2)^2 

		if (&P0->X == &P->X)
		{
			monty_mul(&work->tt3,&P0->Z,&P->Z,n);
			monty_mul(&P0->X,&work->tt2,&P->X,n);
			swap(&P->X,&P->Z);
		}
		else
		{ 
			monty_mul(&work->tt3,&P0->Z,&P->X,n);				// 4*z*(x1*x2-z1*z2)^2 
			monty_mul(&P0->X,&work->tt2,&P->Z,n);				// 4*x*(x2*z1-x1*z2)^2
		}	
		return;
	}

	void pracdup (ecm_pt *P, z *diff, z *sum, ecm_work *work, z *b, z *n)
	{	
		//duplicate a point (x1,z1), where (x1 + z1) and (z1 - z1) have been
		//precomputed and stored in work.
		//output result in (x2,z2)
		//uses b = (A+2)/4

		monty_sqr(sum,&work->tt1,n);					// (x1+z1)^2 
		monty_sqr(diff,&work->tt2,n);					// (x1-z1)^2 
		monty_mul(&work->tt1,&work->tt2,&P->X,n);		// tt1*tt2 = (x1^2 - z1^2)^2 
		monty_sub(&work->tt1,&work->tt2,&work->tt3,n);// tt1-tt2 = 4*x1*z1 
		monty_mul(&work->tt3,b,&work->tt1,n);			// tt3*b = ((A+2)/4*(4*x1*z1)) 
		monty_add(&work->tt1,&work->tt2,&work->tt1,n);// (x1-z1)^2+(A+2)/4*(4*x1*z1) 
		monty_mul(&work->tt3,&work->tt1,&P->Z,n);		// ((4*x1*z1)*((x1-z1)^2+(A+2)/4*(4*x1*z1))) 

		return;
	}

	void prac (ecm_pt *A, unsigned long k, z *b, ecm_work *work, z *n)
	{
		//compute kA from A
		//k > 2

		unsigned long d, e, r;
		ecm_pt *tmp_pt, *B, *C, *t1, *t2;

		B = work->B;
		C = work->C;
		t1 = work->tmp1;
		t2 = work->tmp2;

		//we could try different values of phi here...
		d = k;
		r = (unsigned long) ((double) d / 1.61803398875 + 0.5);

		// Condition 3, then a swap
		d = k - r;
		e = 2 * r - k;

		pt_copy(A,B);
		pt_copy(A,C);
		
		pt_diff(A,&work->diff1,n);
		pt_sum(A,&work->sum1,n);
		pracdup (A, &work->diff1, &work->sum1, work, b,n); 
		while (d != e)
		{
			if (d < e)
			{
				r = d;
				d = e;
				e = r;
				swap(&A->X,&B->X);
				swap(&A->Z,&B->Z);
			}

			/* Do the first line of (Montgomery's) Table 4 whose condition qualifies */
			if (4 * d <= 5 * e && ((d + e) % 3) == 0)
			{ 
				// condition 1
				r = (2 * d - e) / 3;
				e = (2 * e - d) / 3;
				d = r;
				pt_diff(A,&work->diff2,n);
				pt_sum(A,&work->sum2,n);
				pt_diff(B,&work->diff1,n);
				pt_sum(B,&work->sum1,n);
				pracadd (t1, C, work,n); 

				pt_diff(t1,&work->diff1,n);
				pt_sum(t1,&work->sum1,n);
				pracadd (t2, B, work,n); 

				pt_diff(B,&work->diff2,n);
				pt_sum(B,&work->sum2,n);
				pracadd (B,A, work,n); 

				pt_copy(t2,A);
			}
			else if (4 * d <= 5 * e && (d - e) % 6 == 0)
			{ 
				// condition 2 
				d = (d - e) / 2;
				
				pt_diff(A,&work->diff2,n);
				pt_sum(A,&work->sum2,n);
				pt_diff(B,&work->diff1,n);
				pt_sum(B,&work->sum1,n);
				pracadd (B,C, work,n); 
				pracdup (A, &work->diff2, &work->sum2, work, b,n); 
			}
			else if (d <= (4 * e))
			{ 
				// condition 3
				d -= e;

				pt_diff(B,&work->diff2,n);
				pt_sum(B,&work->sum2,n);
				pt_diff(A,&work->diff1,n);
				pt_sum(A,&work->sum1,n);
				pracadd (t1, C, work,n); 

				tmp_pt = B;
				B = t1;
				t1 = C;
				C = tmp_pt;		
			}
			else if ((d + e) % 2 == 0)
			{ 
				// condition 4 
				d = (d - e) / 2;
				
				pt_diff(B,&work->diff2,n);
				pt_sum(B,&work->sum2,n);
				pt_diff(A,&work->diff1,n);
				pt_sum(A,&work->sum1,n);
				pracadd (B,C, work,n); 
				pracdup (A, &work->diff1, &work->sum1, work, b,n); 
			}
			/* now d+e is odd */
			else if (d % 2 == 0)
			{ 
				// condition 5 
				d /= 2;

				pt_diff(C,&work->diff2,n);
				pt_sum(C,&work->sum2,n);
				pt_diff(A,&work->diff1,n);
				pt_sum(A,&work->sum1,n);
				pracadd (C,B, work,n); 
				pracdup (A, &work->diff1, &work->sum1, work, b,n); 
			}
			/* now d is odd, e is even */
			else if (d % 3 == 0)
			{ 
				// condition 6 
				d = d / 3 - e;

				pt_diff(A,&work->diff1,n);
				pt_sum(A,&work->sum1,n);
				pracdup (t1, &work->diff1, &work->sum1, work, b,n);
				
				pt_diff(B,&work->diff2,n);
				pt_sum(B,&work->sum2,n);
				pracadd (t2, C, work,n);
				
				pt_diff(t1,&work->diff2,n);
				pt_sum(t1,&work->sum2,n);
				pracadd (A,A,work,n); 

				pt_diff(t2,&work->diff1,n);
				pt_sum(t2,&work->sum1,n);
				pracadd (t1,C,work,n); 

				tmp_pt = B;
				B = t1;
				t1 = C;
				C = tmp_pt;
			}
			else if ((d + e) % 3 == 0)
			{ 
				// condition 7 
				d = (d - 2 * e) / 3;

				pt_diff(A,&work->diff2,n);
				pt_sum(A,&work->sum2,n);
				pt_diff(B,&work->diff1,n);
				pt_sum(B,&work->sum1,n);
				pracadd (t1,C,work,n);

				pt_diff(t1,&work->diff1,n);
				pt_sum(t1,&work->sum1,n);
				pracadd (B,B,work,n); 
				pracdup (t1, &work->diff2, &work->sum2, work, b,n);

				pt_diff(t1,&work->diff1,n);
				pt_sum(t1,&work->sum1,n);
				pracadd (A,A,work,n); 
			}
			else if ((d - e) % 3 == 0)
			{ 
				// condition 8 
				d = (d - e) / 3;

				pt_diff(A,&work->diff2,n);
				pt_sum(A,&work->sum2,n);
				pt_diff(B,&work->diff1,n);
				pt_sum(B,&work->sum1,n);
				pracadd (t1,C,work,n);

				pt_diff(C,&work->diff1,n);
				pt_sum(C,&work->sum1,n);
				pracadd (C,B,work,n); 
				
				pt_copy(t1,B);
				pracdup (t1, &work->diff2, &work->sum2, work, b,n);

				pt_diff(t1,&work->diff1,n);
				pt_sum(t1,&work->sum1,n);
				pracadd (A,A, work,n);
			}
			else 
			{ 
				// condition 9 - necessarily e is even here
				e /= 2;

				pt_diff(C,&work->diff2,n);
				pt_sum(C,&work->sum2,n);
				pt_diff(B,&work->diff1,n);
				pt_sum(B,&work->sum1,n);
				pracadd (C,A,work,n); 
				pracdup (B, &work->diff1, &work->sum1, work, b,n); 
			}
		}

		pt_diff(A,&work->diff2,n);
		pt_sum(A,&work->sum2,n);
		pt_diff(B,&work->diff1,n);
		pt_sum(B,&work->sum1,n);
		pracadd (A,C,work,n);

		return;
	}

	void *ecm_do_one_curve(void *ptr)
	//uint32 ecm_do_one_curve(z *f, uint32 sigma)
	{
		//attempt to factor n with the elliptic curve method
		//following brent and montgomery's papers, CP's book, 
		//and with inspiration from other implementations including GMP-ECM.

		//unpack the data structure and stuff inside it
		ecm_thread_data_t *thread_data = (ecm_thread_data_t *)ptr;
		ecm_work *work = thread_data->work;
		z *acc = &thread_data->acc;
		z *bb = &thread_data->bb;
		uint8 *marks = thread_data->marks;
		uint8 *nmarks = thread_data->nmarks;
		ecm_pt *P = thread_data->P;
		ecm_pt *Pa = &thread_data->Pa;
		ecm_pt *Pad = &thread_data->Pad;
		z *Paprod = &thread_data->Paprod;
		ecm_pt *Pb = thread_data->Pb;
		z *Pbprod = thread_data->Pbprod;
		ecm_pt *Pd = &thread_data->Pd;
		z *u = &thread_data->u;
		z *v = &thread_data->v;
		z *n = &thread_data->n;

		uint32 sigma = thread_data->sigma;
		int D = thread_data->D;
		z *t1,*t2,*t3,*t4,*t5;
		
		int i,j;
		uint32 q, r;

		int a,b;
		int paired=0;

		//to build the curve, point to some scratch space
		t1 = &work->tt1;
		t2 = &work->tt2;
		t3 = &work->tt3;
		t4 = &work->sum1;
		t5 = &work->sum2;	

		sp2z(sigma,u);
		to_monty(u,n);
		sp2z(4,t1);
		to_monty(t1,n);
		monty_mul(u,t1,v,n);		//v = 4*sigma

		monty_sqr(u,u,n);
		sp2z(5,t1);
		to_monty(t1,n);
		monty_sub(u,t1,u,n);		//u = sigma^2 - 5

		monty_sqr(u,t1,n);
		monty_mul(t1,u,&P->X,n);	//x = u^3

		monty_sqr(v,t1,n);
		monty_mul(t1,v,&P->Z,n);	//z = v^3

		//compute parameter A
		monty_sub(v,u,t1,n);	//(v-u)

		monty_sqr(t1,t2,n);	
		monty_mul(t2,t1,t4,n);	//(v-u)^3

		sp2z(3,t1);
		to_monty(t1,n);
		monty_mul(t1,u,t2,n);	//3u
		monty_add(t2,v,t3,n);	//3u + v

		monty_mul(t3,t4,t1,n);	//a = (v-u)^3 * (3u + v)
		
		sp2z(16,t2);
		to_monty(t2,n);
		monty_mul(&P->X,t2,t3,n);	//16*u^3
		monty_mul(t3,v,t4,n);	//16*u^3*v

		//u holds the denom, t1 holds the numer
		//accomplish the division by multiplying by the modular inverse
		//of the denom, which we find using the gcd on the non-monty
		//representation of the denom and n
		zREDC(t4,n);
		xGCD(t4,n,t2,t3,t5);	//inverse is in &t2

		//gcd should be 1
		if (zCompare(t5,&zOne) != 0)
		{
			printf("WARNING: gcd != 1, aborting this curve\n");
			printf("t1 = %s\n",z2decstr(t1,&gstr1));
			printf("t2 = %s\n",z2decstr(t2,&gstr1));
			printf("t3 = %s\n",z2decstr(t3,&gstr1));
			printf("t4 = %s\n",z2decstr(t4,&gstr1));
			printf("t5 = %s\n",z2decstr(t5,&gstr1));
			printf("n = %s\n",z2decstr(n,&gstr1));
			printf("u = %s\n",z2decstr(u,&gstr1));
			printf("v = %s\n",z2decstr(v,&gstr1));
			printf("sigma = %u\n",sigma);
			return 0;
		}

		//verify inverse
		zModMul(t4,t2,n,t3);
		if (zCompare(t3,&zOne) != 0)
		{
			printf("WARNING: inverse incorrect: %s, aborting this curve\n",z2decstr(t3,&gstr1));
			return 0;
		}

		//compute division via multiplication by inverse
		zREDC(t1,n);
		zMul(t1,t2,t3);
		zDiv(t3,n,t4,t1);
		to_monty(t1,n);
		zCopy(t1,&work->s);
		zCopy(t1,bb);			//b = (A+2)/4 = [(v-u)^3 * (3u + v)]/16*u^3*v

		i=0;
		while (PRIMES[i] < ECM_STG1_MAX)
		{
			q = PRIMES[i];

			//these two cases are easy to handle manually, since the binary pattern is 
			//regular.
			if (q == 2)
			{
				for (r = 2; r <= ECM_STG1_MAX; r *= 2)
				{
					pt_diff(P,&work->diff1,n);
					pt_sum(P,&work->sum1,n);
					pracdup (P, &work->diff1, &work->sum1, work, bb,n);
				}
			}
			else if (q == 3)
			{
				for (r = 3; r <= ECM_STG1_MAX; r *= 3)
				{
					pt_diff(P,&work->diff1,n);
					pt_sum(P,&work->sum1,n);
					pracdup (work->B, &work->diff1, &work->sum1, work, bb,n);

					pt_diff(work->B,&work->diff2,n);
					pt_sum(work->B,&work->sum2,n);
					pracadd (P,P,work,n);
				}
			}
			else
			{
				//do this prime
				for (r = q; r <= ECM_STG1_MAX; r *= q)
					prac (P, (unsigned long) q, bb, work,n);

				//check for abort signal
				if (ECM_ABORT)
					goto done;

			}
			i++;
		}

		j = check_factor(&P->Z,n,&thread_data->factor);
		if (j)
		{
			thread_data->stagefound = 1;
			goto done;
		}

		//goto done;

		//stage 2 init
		//Q = P = result of stage 1
		//compute 2dQ for 0 < d <= D
		
		zCopy(&P->Z,&Pb[1].Z);
		zCopy(&P->X,&Pb[1].X);
		next_pt(&Pb[1],1,work,n);
		a=1;
			
		for (j=2;j<D;j++)
		{
			if (spGCD(j,D) == 1)
			{
				//I should be able to get away with this shortcut...
				//but for some reason it doesn't initialize things correctly
				//for stage 2 so the slower method is still used.

				//zCopy(&Pb[a].Z,&Pb[j].Z);
				//zCopy(&Pb[a].X,&Pb[j].X);
				//next_pt(&Pb[j],j-a);

				zCopy(&P->Z,&Pb[j].Z);
				zCopy(&P->X,&Pb[j].X);
				next_pt(&Pb[j],j,work,n);

				//store Pb[j].X * Pb[j].Z as well
				monty_mul(&Pb[j].X,&Pb[j].Z,&Pbprod[j],n);
				a=j;
			}
		}

		//first a value
		a = PRIMES[i]/D + (PRIMES[i]%D != 0);
		a *= D;

		//compute Pa
		zCopy(&P->Z,&Pa->Z);
		zCopy(&P->X,&Pa->X);
		next_pt(Pa,a,work,n);	
		//and Paprod
		monty_mul(&Pa->X,&Pa->Z,Paprod,n);

		//initialize info needed for giant step
		zCopy(&P->Z,&Pd->Z);
		zCopy(&P->X,&Pd->X);
		next_pt(Pd,D,work,n);

		zCopy(&P->Z,&Pad->Z);
		zCopy(&P->X,&Pad->X);
		next_pt(Pad,a-D,work,n);
				
		//initialize accumulator
		zCopy(&montyconst.one,acc);

		memset(marks,0,D*sizeof(uint8));
		memset(nmarks,0,D*sizeof(uint8));
		paired=0;

		//begin stage 2
		while (PRIMES[i]<ECM_STG2_MAX)
		{
			b = a - PRIMES[i];

			if (!marks[b])
			{
				//not marked, so doesn't have a match on the other side of the previous a
				//accumulate it, and mark the next range of a
				//we accumulate XrZd - XdZr = (Xr - Xd)*(Zr+Zd) + XdZd - XrZr
				//in CP notation, Pa -> (Xr,Zr), Pb -> (Xd,Zd)
				
				if (spGCD(b,D) != 1)
				{
					printf("error spGCD(%d,%d) == 1, please report this bug",b,D);
					exit(-1);
				}

				monty_sub(&Pa->X,&Pb[b].X,&work->tt1,n);
				monty_add(&Pa->Z,&Pb[b].Z,&work->tt2,n);
				monty_mul(&work->tt1,&work->tt2,&work->tt3,n);
				monty_add(&work->tt3,&Pbprod[b],&work->tt1,n);
				monty_sub(&work->tt1,Paprod,&work->tt2,n);
				monty_mul(acc,&work->tt2,&work->tt1,n);
				zCopy(&work->tt1,acc);

				nmarks[D-b]=1;
			}
			else
			{
				//it's marked, so don't need to accumulate it - it has a pair that is already
				//accumulated.
				paired++;
			}


			i++;
			
			if (PRIMES[i] > (uint32)a)
			{
				//check for abort signal
				if (ECM_ABORT)
					goto done;

				//set marks = nextmarks, then clear nextmarks
				memcpy(marks,nmarks,D*sizeof(uint8));
				memset(nmarks,0,D*sizeof(uint8));
				
				//giant step - use the addition formula for ECM
				zCopy(&Pa->X,&work->A->X);			
				zCopy(&Pa->Z,&work->A->Z);		
				
				//Pa + Pd
				//x+ = z- * [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
				//z+ = x- * [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
				//x- = original x (Pad.x)
				//z- = original z (Pad.z)
				monty_add(&Pa->X,&Pa->Z,&work->sum1,n);
				monty_add(&Pd->X,&Pd->Z,&work->sum2,n);
				monty_sub(&Pa->X,&Pa->Z,&work->diff1,n);
				monty_sub(&Pd->X,&Pd->Z,&work->diff2,n);
				add(Pad,Pa,work,n);

				//and Paprod
				monty_mul(&Pa->X,&Pa->Z,Paprod,n);
				
				//Pad holds the previous Pa
				zCopy(&work->A->X,&Pad->X);	
				zCopy(&work->A->Z,&Pad->Z);	

				//next range
				a += D;
			}
		}

		j = check_factor(acc,n,&thread_data->factor);
		if (j)
			thread_data->stagefound = 2;

	done:
		
		return 0;
	}


#else

	int TMP_THREADS;
	uint64 TMP_STG2_MAX;

	void ecm_process_init(z *n)
	{
		//initialize things which all threads will need when using
		//GMP-ECM
		TMP_THREADS = THREADS;
		TMP_STG2_MAX = ECM_STG2_MAX;
		if (THREADS > 1)
		{
			if (VFLAG >= 2)
				printf("GMP-ECM does not support multiple threads... running single threaded\n");
			THREADS = 1;
		}

		return;
	}

	void ecm_thread_init(z *n, int D, ecm_thread_data_t *tdata)
	{
		//initialize everything for all threads using GMP-ECM
		zInit(&tdata->n);
		zInit(&tdata->factor);
		mpz_init(tdata->gmp_n);
		mpz_init(tdata->gmp_factor);
		ecm_init(tdata->params);
		gmp_randseed_ui(tdata->params->rng, get_rand(&g_rand.low, &g_rand.hi));
		zCopy(n,&tdata->n);
		//gmp_randseed_ui(tdata->params->rng, 
		//	get_rand(&obj->seed1, &obj->seed2));

		tdata->params->method = ECM_ECM;
		
		return;
	}

	void ecm_thread_free(ecm_thread_data_t *tdata)
	{
		ecm_clear(tdata->params);
		mpz_clear(tdata->gmp_n);
		mpz_clear(tdata->gmp_factor);
		zFree(&tdata->n);
		zFree(&tdata->factor);

		return;
	}

	void ecm_process_free()
	{
		THREADS = TMP_THREADS;
		ECM_STG2_MAX = TMP_STG2_MAX;
		return;
	}

	void *ecm_do_one_curve(void *ptr)
	//uint32 ecm_do_one_curve(z *f, uint32 sigma)
	{
		int status;
		size_t count;

		//unpack the data structure and stuff inside it
		ecm_thread_data_t *thread_data = (ecm_thread_data_t *)ptr;

		thread_data->params->B1done = 1.0 + floor (1 * 128.) / 134217728.;
		//thread_data->params->verbose = 2;		
		mpz_set_ui(thread_data->params->x, (unsigned long)0);
		mpz_set_ui(thread_data->params->sigma, (unsigned long)0);

#if defined(_WIN64) && BITS_PER_DIGIT == 32
		mpz_import(thread_data->gmp_n, (size_t)(abs(thread_data->n.size)), -1, sizeof(uint32), 
			0, (size_t)0, thread_data->n.val);
#else
		//wrapper for YAFU bigints and call to gmp-ecm
		mp2gmp(&thread_data->n, thread_data->gmp_n);
#endif

		if (ECM_STG2_ISDEFAULT == 0)
		{
			//not default, tell gmp-ecm to use the requested B2
			//printf("using requested B2 value\n");
			sp642z(ECM_STG2_MAX,&thread_data->factor);
			mp2gmp(&thread_data->factor,thread_data->params->B2);
			zClear(&thread_data->factor);
		}

		status = ecm_factor(thread_data->gmp_factor, thread_data->gmp_n,
				ECM_STG1_MAX, thread_data->params);

#if defined(_WIN64) && BITS_PER_DIGIT == 32
		zClear(&thread_data->n);
		mpz_export(thread_data->n.val, &count, -1, sizeof(uint32),
				0, (size_t)0, thread_data->gmp_n);
		thread_data->n.size = count;
#else
		//update n: not sure if gmp-ecm modifies it
		gmp2mp(thread_data->gmp_n,&thread_data->n);
#endif

#if defined(_WIN64) && BITS_PER_DIGIT == 32
		zClear(&thread_data->factor);
		mpz_export(thread_data->factor.val, &count, -1, sizeof(uint32),
				0, (size_t)0, thread_data->params->sigma);
		thread_data->factor.size = count;
#else
		//pull a couple things out of params
		gmp2mp(thread_data->params->sigma,&thread_data->factor);
		thread_data->sigma = (uint32)thread_data->factor.val[0];
#endif

		//printf ("used B2: ");
		//mpz_out_str (stdout, 10, thread_data->params->B2);
		//printf ("\n");

		//NOTE: this required a modification to the GMP-ECM source code in ecm.c
		//in order to get the automatically computed B2 value out of the
		//library
		//gmp2mp(thread_data->params->B2,&thread_data->factor);
		//ECM_STG2_MAX = z264(&thread_data->factor);

#if defined(_WIN64) && BITS_PER_DIGIT == 32
		zClear(&thread_data->factor);
		mpz_export(thread_data->factor.val, &count, -1, sizeof(uint32),
				0, (size_t)0, thread_data->gmp_factor);
		thread_data->factor.size = count;
#else
		//pull out any factor found
		gmp2mp(thread_data->gmp_factor,&thread_data->factor);
#endif

		//the return value is the stage the factor was found in, if no error
		thread_data->stagefound = status;

		return 0;
	}

#endif

// function definitions
void ecmexit(int sig)
{
	printf("\nAborting...\n");
	ECM_ABORT = 1;
	return;
}

int ecm_loop(z *n, int numcurves, fact_obj_t *fobj)
{
	//thread data holds all data needed during sieving
	ecm_thread_data_t *thread_data;		//an array of thread data objects
	z d,t,nn;
	int D;
	FILE *flog;
	int curves_run = 0;
	int i,j;
	//maybe make this an input option: whether or not to stop after
	//finding a factor in the middle of running a requested batch of curves
	int bail_on_factor = 1;
	int bail = 0;
	int charcount = 0, charcount2 = 0;

	//open the log file and annouce we are starting ECM
	flog = fopen(fobj->logname,"a");
	if (flog == NULL)
	{
		printf("could not open %s for writing\n",fobj->logname);
		fclose(flog);
		return 0;
	}

	if (ECM_STG2_MAX > 0xEE6B2800)
	{
		fprintf(stderr,"primes greater than 4e9 not supported in ECM, reducing.\n");
		ECM_STG2_MAX = 0xEE6B2800;
	}

	//check for trivial cases
	if (isZero(n))
	{
		n->type = COMPOSITE;	
		logprint(flog,"Trivial input == 0 in ECM\n");
		fclose(flog);
		return 0;
	}
	if (isOne(n))
	{
		n->type = COMPOSITE;
		logprint(flog,"Trivial input == 1 in ECM\n");
		fclose(flog);
		return 0;
	}
	if (zShortMod(n,3) == 0)
	{
		zShortDiv(n,3,n);
		add_to_factor_list(&zThree);
		logprint(flog,"Trivial factor of 3 found in ECM\n");
		fclose(flog);
		return 0;
	}
	if ((n->val[0] & 0x1) != 1)
	{
		zShortDiv(n,2,n);
		add_to_factor_list(&zTwo);
		logprint(flog,"Trivial factor of 2 found in ECM\n");
		fclose(flog);
		return 0;
	}
	if (isPrime(n))
	{
		//maybe have an input flag to optionally not perform
		//PRP testing (useful for really big inputs)
		n->type = PRP;
		add_to_factor_list(n);
		logprint(flog,"prp%d = %s\n",ndigits(n),z2decstr(n,&gstr1));
		zCopy(&zOne,n);
		fclose(flog);
		return 0;
	}

	//close the log file for until we have something further to report
	fclose(flog);

	//ok, having gotten this far we are now ready to run the requested
	//curves.  initialize the needed data structures, then split
	//the curves up over N threads.  The main thread will farm out different
	//sigmas to the worker threads.

	//init ecm process
	ecm_process_init(n);

	//find the D value used in stg2 ECM.  thread data initialization
	//depends on this value
	if (ECM_STG1_MAX <= 2310)
	{
		if (ECM_STG1_MAX <= 210)
		{
			printf("ECM_STG1_MAX too low\n");
			return 0;
		}
		D = 210;
	}
	else
		D = 2310;

	thread_data = (ecm_thread_data_t *)malloc(THREADS * sizeof(ecm_thread_data_t));
	for (i=0; i<THREADS; i++)
	{
		ecm_thread_init(n,D,&thread_data[i]);
	}

	//init local big ints
	zInit(&d);
	zInit(&t);
	zInit(&nn);

	//round numcurves up so that each thread has something to do every iteration.
	//this prevents a single threaded "cleanup" round after the multi-threaded rounds
	//to finish the leftover requested curves.
	numcurves += numcurves % THREADS;

	if (VFLAG >= 0)
	{
		for (i=0;i<charcount+charcount2;i++)
			printf("\b");

		charcount = printf("ecm: %d curves on C%d input, at "
			,curves_run,ndigits(n));
		charcount2 = print_B1B2();
		fflush(stdout);
	}

	/* activate the threads one at a time. The last is the
	   master thread (i.e. not a thread at all). */

	for (i = 0; i < THREADS - 1; i++)
		ecm_start_worker_thread(thread_data + i, 0);

	ecm_start_worker_thread(thread_data + i, 1);

	//split the requested curves up among the specified number of threads. 
	for (j=0; j < numcurves / THREADS; j++)
	{
		//do work on different sigmas
		for (i=0; i<THREADS; i++)
		{

			if (SIGMA != 0)
			{
				if (curves_run == 0 &&  numcurves > 1)
					printf("WARNING: work will be duplicated with sigma fixed and numcurves > 1\n");
				thread_data[i].sigma = SIGMA;
			}
			else if (get_uvar("sigma",&t))
				thread_data[i].sigma = spRand(6,MAX_DIGIT);
			else
			{
				if (curves_run == 0 &&  numcurves > 1)
					printf("WARNING: work will be duplicated with sigma fixed and numcurves > 1\n");
				thread_data[i].sigma = (uint32)t.val[0];
			}

			if (i == THREADS - 1) {
				ecm_do_one_curve(&thread_data[i]);
			}
			else {
				thread_data[i].command = ECM_COMMAND_RUN;
#if defined(WIN32) || defined(_WIN64)
				SetEvent(thread_data[i].run_event);
#else
				pthread_cond_signal(&thread_data[i].run_cond);
				pthread_mutex_unlock(&thread_data[i].run_lock);
#endif
			}
		}

		//wait for threads to finish
		for (i=0; i<THREADS; i++)
		{
			if (i < THREADS - 1) {
#if defined(WIN32) || defined(_WIN64)
				WaitForSingleObject(thread_data[i].finish_event, INFINITE);
#else
				pthread_mutex_lock(&thread_data[i].run_lock);
				while (thread_data[i].command != ECM_COMMAND_WAIT)
					pthread_cond_wait(&thread_data[i].run_cond, &thread_data[i].run_lock);
#endif
			}
		}

		for (i=0; i<THREADS; i++)
		{
			//look at the result of each curve and see if we're done
			if (zCompare(&thread_data[i].factor,&zOne) > 0 && 
				zCompare(&thread_data[i].factor,n) < 0)
			{
				//non-trivial factor found
				//since we could be doing many curves in parallel,
				//it's possible more than one curve found this factor at the
				//same time.  We divide out factors as we find them, so
				//just check if this factor still divides n
				zCopy(n,&nn);
				zDiv(&nn,&thread_data[i].factor,&t,&d);
				if (zCompare(&d,&zZero) == 0)
				{
					//yes, it does... proceed to record the factor
					
					flog = fopen(fobj->logname,"a");
					if (flog == NULL)
					{
						printf("could not open %s for writing\n",fobj->logname);
						return 0;
					}

					if (isPrime(&thread_data[i].factor))
					{
						thread_data[i].factor.type = PRP;
						add_to_factor_list(&thread_data[i].factor);
						logprint(flog,"prp%d = %s (found in stg%d of curve %d (thread %d) with sigma = %u)\n",
							ndigits(&thread_data[i].factor),
							z2decstr(&thread_data[i].factor,&gstr1),
							thread_data[i].stagefound,curves_run+1,i,thread_data[i].sigma);
					}
					else
					{
						thread_data[i].factor.type = COMPOSITE;
						add_to_factor_list(&thread_data[i].factor);
						logprint(flog,"c%d = %s (found in stg%d of curve %d (thread %d) with sigma = %u)\n",
							ndigits(&thread_data[i].factor),
							z2decstr(&thread_data[i].factor,&gstr1),
							thread_data[i].stagefound,curves_run+1,i,thread_data[i].sigma);
					}

					fclose(flog);
			
					//reduce input
					zDiv(n,&thread_data[i].factor,&t,&d);
					zCopy(&t,n);
					zCopy(n,&thread_data[i].n);

					//we found a factor and might want to stop,
					//but should check all threads output first
					if (bail_on_factor)
						bail = 1;
					else if (isPrime(n))
						bail = 1;
					else
					{
						//found a factor and the cofactor is composite.
						//the user has specified to keep going with ECM until the 
						//curve counts are finished thus:
						//we need to re-initialize with a different modulus.  this is
						//independant of the thread data initialization
						ecm_process_free();
						ecm_process_init(n);
					}
				}
			}

			curves_run++;
		}

		if (bail)
			goto done;

		if (VFLAG >= 0)
		{
			for (i=0;i<charcount+charcount2;i++)
				printf("\b");

			charcount = printf("ecm: %d curves on C%d input, at "
				,curves_run,ndigits(n));
			charcount2 = print_B1B2();
			fflush(stdout);
		}

	}

done:
	if (VFLAG >= 0)
		printf("\n");

	flog = fopen(fobj->logname,"a");
	if (flog == NULL)
	{
		printf("could not open %s for writing\n",fobj->logname);
		fclose(flog);
		return 0;
	}

	logprint(flog,"Finished %d curves using Lenstra ECM method on C%d input, B1 = %u, B2 = %lu\n",
		curves_run,ndigits(n),ECM_STG1_MAX,ECM_STG2_MAX);

	fclose(flog);

	//stop worker threads
	for (i=0; i<THREADS - 1; i++)
	{
		ecm_stop_worker_thread(thread_data + i, 0);
	}
	ecm_stop_worker_thread(thread_data + i, 1);

	for (i=0; i<THREADS; i++)
		ecm_thread_free(&thread_data[i]);
	free(thread_data);

	zFree(&d);
	zFree(&t);
	zFree(&nn);
	
	ecm_process_free();

	return curves_run;
}

int print_B1B2(void)
{
	int i;
	char suffix;
	char stg1str[20];
	char stg2str[20];

	if (ECM_STG1_MAX % 1000000000 == 0)
	{
		suffix = 'B';
		sprintf(stg1str,"%u%c",ECM_STG1_MAX / 1000000000, suffix);
	}
	else if (ECM_STG1_MAX % 1000000 == 0)
	{
		suffix = 'M';
		sprintf(stg1str,"%u%c",ECM_STG1_MAX / 1000000, suffix);
	}
	else if (ECM_STG1_MAX % 1000 == 0)
	{
		suffix = 'K';
		sprintf(stg1str,"%u%c",ECM_STG1_MAX / 1000, suffix);
	}
	else
	{
		sprintf(stg1str,"%u",ECM_STG1_MAX);
	}

	if (ECM_STG2_MAX % 1000000000 == 0)
	{
		suffix = 'B';
#if defined(__unix__) && (BITS_PER_DIGIT == 64)
		sprintf(stg2str,"%lu%c",ECM_STG2_MAX / 1000000000, suffix);
#elif defined(__unix__) && (BITS_PER_DIGIT == 32)
		sprintf(stg2str,"%llu%c",ECM_STG2_MAX / 1000000000, suffix);
#else
		sprintf(stg2str,"%I64u%c",ECM_STG2_MAX / 1000000000, suffix);
#endif
	}
	else if (ECM_STG2_MAX % 1000000 == 0)
	{
		suffix = 'M';
#if defined(__unix__) && (BITS_PER_DIGIT == 64)
		sprintf(stg2str,"%lu%c",ECM_STG2_MAX / 1000000, suffix);
#elif defined(__unix__) && (BITS_PER_DIGIT == 32)
		sprintf(stg2str,"%llu%c",ECM_STG2_MAX / 1000000, suffix);
#else
		sprintf(stg2str,"%I64u%c",ECM_STG2_MAX / 1000000, suffix);
#endif
	}
	else if (ECM_STG2_MAX % 1000 == 0)
	{
		suffix = 'K';
#if defined(__unix__) && (BITS_PER_DIGIT == 64)
		sprintf(stg2str,"%lu%c",ECM_STG2_MAX / 1000, suffix);
#elif defined(__unix__) && (BITS_PER_DIGIT == 32)
		sprintf(stg2str,"%llu%c",ECM_STG2_MAX / 1000, suffix);
#else
		sprintf(stg2str,"%I64u%c",ECM_STG2_MAX / 1000, suffix);
#endif
	}
	else
	{
#if defined(__unix__) && (BITS_PER_DIGIT == 64)
		sprintf(stg2str,"%lu",ECM_STG2_MAX);
#elif defined(__unix__) && (BITS_PER_DIGIT == 32)
		sprintf(stg2str,"%llu",ECM_STG2_MAX);
#else
		sprintf(stg2str,"%I64u",ECM_STG2_MAX);
#endif
	}

	i = printf("B1 = %s, B2 = %s",stg1str,stg2str);

	return i;
}



