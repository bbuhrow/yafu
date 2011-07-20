/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

       				   --bbuhrow@gmail.com 7/1/10
----------------------------------------------------------------------*/

#include "soe.h"


void start_soe_worker_thread(thread_soedata_t *t, uint32 is_master_thread) 
{

	/* create a thread that will process a sieve line. The last line does 
	   not get its own thread (the current thread handles it) */

	if (is_master_thread) {
		return;
	}

	t->command = SOE_COMMAND_INIT;
#if defined(WIN32) || defined(_WIN64)
	t->run_event = CreateEvent(NULL, FALSE, TRUE, NULL);
	t->finish_event = CreateEvent(NULL, FALSE, FALSE, NULL);
	t->thread_id = CreateThread(NULL, 0, soe_worker_thread_main, t, 0, NULL);

	WaitForSingleObject(t->finish_event, INFINITE); /* wait for ready */
#else
	pthread_mutex_init(&t->run_lock, NULL);
	pthread_cond_init(&t->run_cond, NULL);

	pthread_cond_signal(&t->run_cond);
	pthread_mutex_unlock(&t->run_lock);
	pthread_create(&t->thread_id, NULL, soe_worker_thread_main, t);

	pthread_mutex_lock(&t->run_lock); /* wait for ready */
	while (t->command != SOE_COMMAND_WAIT)
		pthread_cond_wait(&t->run_cond, &t->run_lock);
#endif
}

void stop_soe_worker_thread(thread_soedata_t *t, uint32 is_master_thread)
{
	if (is_master_thread) {
		return;
	}

	t->command = SOE_COMMAND_END;
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
DWORD WINAPI soe_worker_thread_main(LPVOID thread_data) {
#else
void *soe_worker_thread_main(void *thread_data) {
#endif
	thread_soedata_t *t = (thread_soedata_t *)thread_data;

	while(1) {

		/* wait forever for work to do */
#if defined(WIN32) || defined(_WIN64)
		WaitForSingleObject(t->run_event, INFINITE);
#else
		pthread_mutex_lock(&t->run_lock);
		while (t->command == SOE_COMMAND_WAIT) {
			pthread_cond_wait(&t->run_cond, &t->run_lock);
		}
#endif
		/* do work */

		if (t->command == SOE_COMMAND_SIEVE_AND_COUNT)
		{
			sieve_line(t);
#ifdef DO_SPECIAL_COUNT
			count_line_special(t);
#else
			count_line(t);
#endif
		}
		else if (t->command == SOE_COMMAND_SIEVE_AND_COMPUTE)
		{
			sieve_line(t);
			primes_from_lineflags(t);
		}
		else if (t->command == SOE_COMMAND_END)
			break;

		/* signal completion */

		t->command = SOE_COMMAND_WAIT;
#if defined(WIN32) || defined(_WIN64)
		SetEvent(t->finish_event);
#else
		pthread_cond_signal(&t->run_cond);
		pthread_mutex_unlock(&t->run_lock);
#endif
	}

#if defined(WIN32) || defined(_WIN64)
	return 0;
#else
	return NULL;
#endif
}
