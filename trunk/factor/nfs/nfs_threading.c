#include "nfs.h"

void nfs_start_worker_thread(nfs_threaddata_t *t, 
				uint32 is_master_thread) {

	/* create a thread that will process a polynomial. The last poly does 
	   not get its own thread (the current thread handles it) */

	if (is_master_thread == 1) {
		return;
	}

	t->command = NFS_COMMAND_INIT;
#if defined(WIN32) || defined(_WIN64)
		
	// specific to different structure of poly selection threading
	if (is_master_thread == 2)
	{
		t->run_event = CreateEvent(NULL, FALSE, FALSE, NULL);
		t->finish_event = CreateEvent(NULL, FALSE, FALSE, NULL);
		*t->queue_event = CreateEvent(NULL, FALSE, FALSE, NULL);
	}
	else
	{
		t->run_event = CreateEvent(NULL, FALSE, TRUE, NULL);
		t->finish_event = CreateEvent(NULL, FALSE, FALSE, NULL);
	}

	t->thread_id = CreateThread(NULL, 0, nfs_worker_thread_main, t, 0, NULL);

	WaitForSingleObject(t->finish_event, INFINITE); /* wait for ready */
#else
	pthread_mutex_init(&t->run_lock, NULL);
	pthread_cond_init(&t->run_cond, NULL);

	if (is_master_thread == 0)
	{
		pthread_cond_signal(&t->run_cond);
		pthread_mutex_unlock(&t->run_lock);
	}

	pthread_create(&t->thread_id, NULL, nfs_worker_thread_main, t);

	pthread_mutex_lock(&t->run_lock); /* wait for ready */
	while (t->command != NFS_COMMAND_WAIT)
		pthread_cond_wait(&t->run_cond, &t->run_lock);

	if (is_master_thread == 2)
		pthread_mutex_unlock(&t->run_lock);

#endif

}

void nfs_stop_worker_thread(nfs_threaddata_t *t,
				uint32 is_master_thread)
{
	if (is_master_thread == 1) {
		return;
	}

	t->command = NFS_COMMAND_END;
#if defined(WIN32) || defined(_WIN64)
	SetEvent(t->run_event);
	WaitForSingleObject(t->thread_id, INFINITE);
	CloseHandle(t->thread_id);
	CloseHandle(t->run_event);
	CloseHandle(t->finish_event);

	// specific to different structure of poly selection threading
	if (is_master_thread == 2)
		CloseHandle(*t->queue_event);
#else
	if (is_master_thread == 2)
		pthread_mutex_lock(&t->run_lock);

	pthread_cond_signal(&t->run_cond);
	pthread_mutex_unlock(&t->run_lock);
	pthread_join(t->thread_id, NULL);
	pthread_cond_destroy(&t->run_cond);
	pthread_mutex_destroy(&t->run_lock);
#endif

}

#if defined(WIN32) || defined(_WIN64)
DWORD WINAPI nfs_worker_thread_main(LPVOID thread_data) {
#else
void *nfs_worker_thread_main(void *thread_data) {
#endif
	nfs_threaddata_t *t = (nfs_threaddata_t *)thread_data;

	/*
    * Respond to the master thread that we're ready for work. If we had any thread-
    * specific initialization which needed to be done, it would go before this signal.
    */

	// specific to different structure of poly selection threading
	if (t->is_poly_select)
	{
#if defined(WIN32) || defined(_WIN64)
		t->command = NFS_COMMAND_WAIT;
		SetEvent(t->finish_event);
#else
		pthread_mutex_lock(&t->run_lock);
		t->command = NFS_COMMAND_WAIT;
		pthread_cond_signal(&t->run_cond);
		pthread_mutex_unlock(&t->run_lock);
#endif
	}

	while(1) {

		/* wait forever for work to do */
#if defined(WIN32) || defined(_WIN64)
		WaitForSingleObject(t->run_event, INFINITE);		
#else
		pthread_mutex_lock(&t->run_lock);
		while (t->command == NFS_COMMAND_WAIT) {
			pthread_cond_wait(&t->run_cond, &t->run_lock);
		}
#endif
		/* do work */

		if (t->command == NFS_COMMAND_RUN)
			lasieve_launcher(t);
		else if (t->command == NFS_COMMAND_RUN_POLY)
			polyfind_launcher(t);
		else if (t->command == NFS_COMMAND_END)
			break;

		/* signal completion */
		t->command = NFS_COMMAND_WAIT;

		if (t->is_poly_select)
		{
#if defined(WIN32) || defined(_WIN64)

			WaitForSingleObject( 
				*t->queue_lock,    // handle to mutex
				INFINITE);  // no time-out interval
 			
			t->thread_queue[(*(t->threads_waiting))++] = t->tindex;

			SetEvent(*t->queue_event);			

			ReleaseMutex(*t->queue_lock);
		
#else
			pthread_mutex_unlock(&t->run_lock);

			// lock the work queue and insert my thread ID into it
			// this tells the master that my results should be collected
			// and I should be dispatched another polynomial
			pthread_mutex_lock(t->queue_lock);
			t->thread_queue[(*(t->threads_waiting))++] = t->tindex;
			pthread_cond_signal(t->queue_cond);
			pthread_mutex_unlock(t->queue_lock);
#endif
		}
		else
		{

#if defined(WIN32) || defined(_WIN64)
			SetEvent(t->finish_event);		
#else
			pthread_cond_signal(&t->run_cond);
			pthread_mutex_unlock(&t->run_lock);
#endif

		}
	}

#if defined(WIN32) || defined(_WIN64)
	return 0;
#else
	return NULL;
#endif
}
