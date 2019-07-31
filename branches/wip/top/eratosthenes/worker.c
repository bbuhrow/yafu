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


#ifdef USE_SOE_THREADPOOL
void start_soe_worker_thread(thread_soedata_t *t) {

    //create a thread that will process a polynomial     
    t->command = SOE_COMMAND_INIT;
#if defined(WIN32) || defined(_WIN64)
    t->run_event = CreateEvent(NULL, FALSE, FALSE, NULL);
    t->finish_event = CreateEvent(NULL, FALSE, FALSE, NULL);
    *t->queue_event = CreateEvent(NULL, FALSE, FALSE, NULL);
    t->thread_id = CreateThread(NULL, 0, soe_worker_thread_main, t, 0, NULL);
    WaitForSingleObject(t->finish_event, INFINITE); /* wait for ready */
#else
    pthread_mutex_init(&t->run_lock, NULL);
    pthread_cond_init(&t->run_cond, NULL);

    pthread_create(&t->thread_id, NULL, soe_worker_thread_main, t);

    pthread_mutex_lock(&t->run_lock); /* wait for ready */
    while (t->command != SOE_COMMAND_WAIT)
        pthread_cond_wait(&t->run_cond, &t->run_lock);
    pthread_mutex_unlock(&t->run_lock);
#endif
}

void stop_soe_worker_thread(thread_soedata_t *t)
{

#if defined(WIN32) || defined(_WIN64)
    t->command = SOE_COMMAND_END;
    SetEvent(t->run_event);
    WaitForSingleObject(t->thread_id, INFINITE);
    CloseHandle(t->thread_id);
    CloseHandle(t->run_event);
    CloseHandle(t->finish_event);
    CloseHandle(*t->queue_event);
#else
    pthread_mutex_lock(&t->run_lock);
    t->command = SOE_COMMAND_END;
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

    /*
    * Respond to the master thread that we're ready for work. If we had any thread-
    * specific initialization which needed to be done, it would go before this signal.
    */
#if defined(WIN32) || defined(_WIN64)
    t->command = SOE_COMMAND_WAIT;
    SetEvent(t->finish_event);
#else
    pthread_mutex_lock(&t->run_lock);
    t->command = SOE_COMMAND_WAIT;
    pthread_cond_signal(&t->run_cond);
    pthread_mutex_unlock(&t->run_lock);
#endif

    while(1) {
        uint32 i;

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
            t->sdata.lines[t->current_line] = 
                (uint8 *)xmalloc_align(t->sdata.numlinebytes * sizeof(uint8));
            sieve_line_ptr(t);
            t->linecount = count_line(&t->sdata, t->current_line);
            free(t->sdata.lines[t->current_line]);
        }
        else if (t->command == SOE_COMMAND_SIEVE_AND_COMPUTE)
        {
            sieve_line_ptr(t);
        }
        else if (t->command == SOE_COMPUTE_ROOTS)
        {
            if (VFLAG > 2)
            {
                printf("starting root computation over %u to %u\n", t->startid, t->stopid);
            }

            if (t->sdata.sieve_range == 0)
            {
                for (i = t->startid; i < t->stopid; i++)
                {
                    uint32 inv;
                    uint32 prime = t->sdata.sieve_p[i];

                    // slightly optimized modinv when prime >> prodN
                    inv = modinv_1c(t->sdata.prodN, prime);
                    t->sdata.root[i] = prime - inv;

                    //t->sdata.lower_mod_prime[i - t->sdata.bucket_start_id] = 
                    t->sdata.lower_mod_prime[i] =
                        (t->sdata.lowlimit + 1) % prime;
                }
            }
            else
            {
                mpz_t tmpz;
                mpz_init(tmpz);

                mpz_add_ui(tmpz, *t->sdata.offset, t->sdata.lowlimit + 1);
                for (i = t->startid; i < t->stopid; i++)
                {
                    uint32 inv;
                    uint32 prime = t->sdata.sieve_p[i];

                    // slightly optimized modinv when prime >> prodN
                    inv = modinv_1c(t->sdata.prodN, prime);
                    t->sdata.root[i] = prime - inv;

                    //t->sdata.lower_mod_prime[i - t->sdata.bucket_start_id] = 
                    t->sdata.lower_mod_prime[i] =
                        mpz_tdiv_ui(tmpz, prime);
                }
            }

        }
        else if (t->command == SOE_COMPUTE_PRIMES)
        {
            t->linecount = 0;

            for (i = t->startid; i < t->stopid; i+=8)
            {
                t->linecount = compute_8_bytes(&t->sdata, t->linecount, t->ddata.primes, i);		
            }
        }
        else if (t->command == SOE_COMPUTE_PRPS)
        {
            t->linecount = 0;
            for (i = t->startid; i < t->stopid; i++)
            {
                mpz_add_ui(t->tmpz, t->offset, t->ddata.primes[i - t->startid]);
                if ((mpz_cmp(t->tmpz, t->lowlimit) >= 0) && (mpz_cmp(t->highlimit, t->tmpz) >= 0))
                {
                    //if (is_mpz_prp(t->tmpz))
                    if (mpz_extrastrongbpsw_prp(t->tmpz))
                        t->ddata.primes[t->linecount++] = t->ddata.primes[i - t->startid];
                }
            }
        }
        else if (t->command == SOE_COMMAND_END)
            break;

        /* signal completion */

        t->command = SOE_COMMAND_WAIT;
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

#if defined(WIN32) || defined(_WIN64)
    return 0;
#else
    return NULL;
#endif
}


#else

void start_soe_worker_thread(thread_soedata_t *t) 
{

	// create a thread that will process a sieve line.

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

void stop_soe_worker_thread(thread_soedata_t *t)
{

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
		uint32 i;

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
			t->sdata.lines[t->current_line] = 
				(uint8 *)malloc(t->sdata.numlinebytes * sizeof(uint8));
			sieve_line_ptr(t);
			t->linecount = count_line(&t->sdata, t->current_line);
			free(t->sdata.lines[t->current_line]);
		}
		else if (t->command == SOE_COMMAND_SIEVE_AND_COMPUTE)
		{
			sieve_line_ptr(t);
		}
		else if (t->command == SOE_COMPUTE_ROOTS)
		{
			if (VFLAG > 2)
				printf("starting root computation over %u to %u\n", t->startid, t->stopid);

			if (t->sdata.sieve_range == 0)
			{
				for (i = t->startid; i < t->stopid; i++)
				{
					uint32 inv;
					uint32 prime = t->sdata.sieve_p[i];

                    // slightly optimized modinv when prime >> prodN
                    inv = modinv_1c(t->sdata.prodN, prime);
					t->sdata.root[i] = prime - inv;

					//t->sdata.lower_mod_prime[i - t->sdata.bucket_start_id] = 
                    t->sdata.lower_mod_prime[i] =
						(t->sdata.lowlimit + 1) % prime;
				}
			}
			else
			{
				mpz_t tmpz;
				mpz_init(tmpz);

				mpz_add_ui(tmpz, *t->sdata.offset, t->sdata.lowlimit + 1);
				for (i = t->startid; i < t->stopid; i++)
				{
					uint32 inv;
					uint32 prime = t->sdata.sieve_p[i];

                    // slightly optimized modinv when prime >> prodN
                    inv = modinv_1c(t->sdata.prodN, prime);
					t->sdata.root[i] = prime - inv;
		
					//t->sdata.lower_mod_prime[i - t->sdata.bucket_start_id] = 
                    t->sdata.lower_mod_prime[i] =
						mpz_tdiv_ui(tmpz, prime);
				}
			}

		}
		else if (t->command == SOE_COMPUTE_PRIMES)
		{
			t->linecount = 0;

			for (i = t->startid; i < t->stopid; i+=8)
			{
				t->linecount = compute_8_bytes(&t->sdata, t->linecount, t->ddata.primes, i);		
			}
		}
		else if (t->command == SOE_COMPUTE_PRPS)
		{
			t->linecount = 0;
			for (i = t->startid; i < t->stopid; i++)
			{
				mpz_add_ui(t->tmpz, t->offset, t->ddata.primes[i - t->startid]);
				if ((mpz_cmp(t->tmpz, t->lowlimit) >= 0) && (mpz_cmp(t->highlimit, t->tmpz) >= 0))
				{
					//if (is_mpz_prp(t->tmpz))
                    if (mpz_extrastrongbpsw_prp(t->tmpz))
						t->ddata.primes[t->linecount++] = t->ddata.primes[i - t->startid];
				}
			}
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

#endif
