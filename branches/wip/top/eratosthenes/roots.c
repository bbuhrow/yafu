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
#include "arith.h"

void getRoots(soe_staticdata_t *sdata, thread_soedata_t *thread_data)
{
    int prime, prodN;
    uint64 startprime;
    uint64 lblk_b, ublk_b, blk_b_sqrt;
    uint64 i;
    int j, k;
    uint32 range, lastid;

    //timing
    double t;
    struct timeval tstart, tstop;
    TIME_DIFF *	difference;

#ifdef USE_SOE_THREADPOOL
    // thread work-queue controls
    int threads_working = 0;
    int *thread_queue;
    int *threads_waiting;
#if defined(WIN32) || defined(_WIN64)
    HANDLE queue_lock;
    HANDLE *queue_events = NULL;
#else
    pthread_mutex_t queue_lock;
    pthread_cond_t queue_cond;
#endif

#endif

    prodN = (int)sdata->prodN;
    startprime = sdata->startprime;

    gettimeofday(&tstart, NULL);

    lblk_b = sdata->lowlimit;
    ublk_b = sdata->blk_r + lblk_b - sdata->prodN;
    blk_b_sqrt = (uint32)(sqrt(ublk_b + sdata->prodN)) + 1;

    for (i = startprime; i < sdata->bucket_start_id; i++)
    {
        uint32 inv;
        prime = sdata->sieve_p[i];

        //sieving requires that we find the offset of each sieve prime in each block 
        //that we sieve.  We are more restricted in choice of offset because we
        //sieve residue classes.  A good way to find the offset is the extended 
        //euclidean algorithm, which reads ax + by = gcd(a,b),
        //where a = prime, b = prodN, and therefore gcd = 1.  
        //since a and b are coprime, y is the multiplicative inverse of prodN modulo prime.
        //This value is a constant, so compute it here in order to facilitate 
        //finding offsets later.

        if (sdata->sieve_p[i] > blk_b_sqrt)
        {
            lblk_b = ublk_b + prodN;
            ublk_b += sdata->blk_r;
            blk_b_sqrt = (uint64)(sqrt((int64)(ublk_b + prodN))) + 1;
        }

        //solve prodN ^ -1 % p 
        inv = modinv_1(prodN, prime);
        sdata->root[i] = prime - inv;

        // as for larger primes we pull this division out of the
        // get_offsets loop.  This requires that we update the
        // lower block boundary (above) in which sieve primes are 
        // first used as we go.
        sdata->lower_mod_prime[i] =
            (lblk_b + 1) % prime;
    }

    gettimeofday(&tstop, NULL);

    difference = my_difftime(&tstart, &tstop);
    t = ((double)difference->secs + (double)difference->usecs / 1000000);
    free(difference);

    if (VFLAG > 2)
    {
        printf("time to compute linear sieve roots = %1.2f\n", t);
    }

    gettimeofday(&tstart, NULL);

    range = (sdata->pboundi - sdata->bucket_start_id) / THREADS;
    lastid = sdata->bucket_start_id;

    // divvy up the primes
    for (j = 0; j < THREADS; j++)
    {
        thread_soedata_t *t = thread_data + j;

        t->sdata = *sdata;
        t->startid = lastid;
        t->stopid = t->startid + range;
        lastid = t->stopid;

        if (VFLAG > 2)
        {
            printf("thread %d starting root computation over %u to %u\n",
                j, t->startid, t->stopid); fflush(stdout);
        }
    }

    // the last one gets any leftover
    if (thread_data[THREADS - 1].stopid != sdata->pboundi)
    {
        thread_data[THREADS - 1].stopid = sdata->pboundi;
    }


#ifdef USE_SOE_THREADPOOL

    // allocate the queue of threads waiting for work
    thread_queue = (int *)malloc(THREADS * sizeof(int));
    threads_waiting = (int *)malloc(sizeof(int));

    if (THREADS > 1)
    {
#if defined(WIN32) || defined(_WIN64)
        queue_lock = CreateMutex(
            NULL,              // default security attributes
            FALSE,             // initially not owned
            NULL);             // unnamed mutex
        queue_events = (HANDLE *)malloc(THREADS * sizeof(HANDLE));
#else
        pthread_mutex_init(&queue_lock, NULL);
        pthread_cond_init(&queue_cond, NULL);
#endif
    }

    for (i = 0; i < THREADS; i++)
    {
        thread_data[i].tindex = i;
        thread_data[i].tstartup = 1;
        // assign all thread's a pointer to the waiting queue.  access to 
        // the array will be controlled by a mutex
        thread_data[i].thread_queue = thread_queue;
        thread_data[i].threads_waiting = threads_waiting;

        if (THREADS > 1)
        {
#if defined(WIN32) || defined(_WIN64)
            // assign a pointer to the mutex
            thread_data[i].queue_lock = &queue_lock;
            thread_data[i].queue_event = &queue_events[i];
#else
            thread_data[i].queue_lock = &queue_lock;
            thread_data[i].queue_cond = &queue_cond;
#endif
        }
    }

    if (THREADS > 1)
    {
        // Activate the worker threads one at a time. 
        // Initialize the work queue to say all threads are waiting for work
        for (i = 0; i < THREADS; i++)
        {
            start_soe_worker_thread(thread_data + i);
            thread_queue[i] = i;
        }
    }
    *threads_waiting = THREADS;

    if (THREADS > 1)
    {
#if defined(WIN32) || defined(_WIN64)
        // nothing
#else
        pthread_mutex_lock(&queue_lock);
#endif
    }

    k = 0;
    //printf("=== starting threaded sieve\n");
    while (1)
    {

        // Process threads until there are no more waiting for their results to be collected
        while (*threads_waiting > 0)
        {
            int tid;

            if (THREADS > 1)
            {
                // Pop a waiting thread off the queue (OK, it's stack not a queue)
#if defined(WIN32) || defined(_WIN64)
                WaitForSingleObject(
                    queue_lock,    // handle to mutex
                    INFINITE);  // no time-out interval
#endif

                tid = thread_queue[--(*threads_waiting)];

#if defined(WIN32) || defined(_WIN64)
                ReleaseMutex(queue_lock);
#endif
            }
            else
            {
                tid = 0;
            }

            // if not in startup...
            if (thread_data[tid].tstartup == 0)
            {

                // this thread is done, so decrement the count of working threads
                threads_working--;
            }
            else
            {
                thread_data[tid].tstartup = 0;
            }

            // don't really need a threadpool for this, since we just assign
            // one chunk of work to each thread, but we need to use the existing
            // infrastructure.  In the future, process smaller batches in a pool.
            if (k < THREADS)
            {
                thread_soedata_t *t = thread_data + tid;

                t->command = SOE_COMPUTE_ROOTS;

                if (THREADS > 1)
                {
                    // send the thread a signal to start processing the poly we just generated for it
#if defined(WIN32) || defined(_WIN64)
                    SetEvent(thread_data[tid].run_event);
#else
                    pthread_mutex_lock(&thread_data[tid].run_lock);
                    pthread_cond_signal(&thread_data[tid].run_cond);
                    pthread_mutex_unlock(&thread_data[tid].run_lock);
#endif
                }

                // this thread is now busy, so increment the count of working threads
                threads_working++;

                // and keep track of where we are overall.
                k++;
            }

            if (THREADS == 1)
                *threads_waiting = 0;

        } // while (*threads_waiting > 0)

        // if all threads are done, break out
        if (threads_working == 0)
            break;

        if (THREADS > 1)
        {
            // wait for a thread to finish and put itself in the waiting queue
#if defined(WIN32) || defined(_WIN64)
            j = WaitForMultipleObjects(
                THREADS,
                queue_events,
                FALSE,
                INFINITE);
#else
            pthread_cond_wait(&queue_cond, &queue_lock);
#endif
        }
        else
        {
            //do some work
            thread_soedata_t *t = thread_data + 0;

            // run in the current thread
            // bucket sieved primes need more data
            if (sdata->sieve_range == 0)
            {
                for (i = t->startid; i < t->stopid; i++)
                {
                    uint32 inv;
                    prime = t->sdata.sieve_p[i];

                    //sieving requires that we find the offset of each sieve prime in each block 
                    //that we sieve.  We are more restricted in choice of offset because we
                    //sieve residue classes.  A good way to find the offset is the extended 
                    //euclidean algorithm, which reads ax + by = gcd(a,b),
                    //where a = prime, b = prodN, and therefore gcd = 1.  
                    //since a and b are coprime, y is the multiplicative inverse of prodN modulo prime.
                    //This value is a constant, so compute it here in order to facilitate 
                    //finding offsets later.

                    //solve prodN ^ -1 % p 
                    // slightly optimized modinv when prime >> prodN
                    inv = modinv_1c(prodN, prime);
                    t->sdata.root[i] = prime - inv;

                    //we can also speed things up by computing and storing the residue
                    //mod p of the first sieve location in the first residue class.  This provides
                    //a speedup by pulling this constant (involving a division) out of a critical loop
                    //when finding offsets of bucket sieved primes.
                    //these are only used by bucket sieved primes.
                    // - t->sdata.bucket_start_id] = 
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
                    prime = t->sdata.sieve_p[i];

                    //sieving requires that we find the offset of each sieve prime in each block 
                    //that we sieve.  We are more restricted in choice of offset because we
                    //sieve residue classes.  A good way to find the offset is the extended 
                    //euclidean algorithm, which reads ax + by = gcd(a,b),
                    //where a = prime, b = prodN, and therefore gcd = 1.  
                    //since a and b are coprime, y is the multiplicative inverse of prodN modulo prime.
                    //This value is a constant, so compute it here in order to facilitate 
                    //finding offsets later.

                    //solve prodN ^ -1 % p 
                    // slightly optimized modinv when prime >> prodN
                    inv = modinv_1c(prodN, prime);
                    t->sdata.root[i] = prime - inv;

                    //we can also speed things up by computing and storing the residue
                    //mod p of the first sieve location in the first residue class.  This provides
                    //a speedup by pulling this constant (involving a division) out of a critical loop
                    //when finding offsets of bucket sieved primes.
                    //these are only used by bucket sieved primes.			
                    //t->sdata.lower_mod_prime[i - t->sdata.bucket_start_id] = 
                    t->sdata.lower_mod_prime[i] =
                        mpz_tdiv_ui(tmpz, prime);
                }
            }

            *threads_waiting = 1;
        }
    }

    //printf("=== computing finished\n");

    //stop worker threads
    for (i = 0; i < THREADS; i++)
    {
        if (THREADS > 1)
            stop_soe_worker_thread(thread_data + i);
    }

    //printf("=== threading stopped\n");

    free(thread_queue);
    free(threads_waiting);

#if defined(WIN32) || defined(_WIN64)
    if (THREADS > 1)
        free(queue_events);
#endif

#else

    // start the threads
    for (i = 0; i < THREADS - 1; i++)
    {
        start_soe_worker_thread(thread_data + i);
    }

    // now run with the threads
    for (j = 0; j < THREADS; j++)
    {
        thread_soedata_t *t = thread_data + j;

        if (j == (THREADS - 1)) 
        {	
            // run in the current thread
            // bucket sieved primes need more data
            if (sdata->sieve_range == 0)
            {				
                for (i = t->startid; i < t->stopid; i++)
                {
                    uint32 inv;
                    prime = t->sdata.sieve_p[i];

                    //sieving requires that we find the offset of each sieve prime in each block 
                    //that we sieve.  We are more restricted in choice of offset because we
                    //sieve residue classes.  A good way to find the offset is the extended 
                    //euclidean algorithm, which reads ax + by = gcd(a,b),
                    //where a = prime, b = prodN, and therefore gcd = 1.  
                    //since a and b are coprime, y is the multiplicative inverse of prodN modulo prime.
                    //This value is a constant, so compute it here in order to facilitate 
                    //finding offsets later.

                    //solve prodN ^ -1 % p 
                    // slightly optimized modinv when prime >> prodN
                    inv = modinv_1c(prodN, prime);
                    t->sdata.root[i] = prime - inv;

                    //we can also speed things up by computing and storing the residue
                    //mod p of the first sieve location in the first residue class.  This provides
                    //a speedup by pulling this constant (involving a division) out of a critical loop
                    //when finding offsets of bucket sieved primes.
                    //these are only used by bucket sieved primes.
                    // - t->sdata.bucket_start_id] = 
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
                    prime = t->sdata.sieve_p[i];

                    //sieving requires that we find the offset of each sieve prime in each block 
                    //that we sieve.  We are more restricted in choice of offset because we
                    //sieve residue classes.  A good way to find the offset is the extended 
                    //euclidean algorithm, which reads ax + by = gcd(a,b),
                    //where a = prime, b = prodN, and therefore gcd = 1.  
                    //since a and b are coprime, y is the multiplicative inverse of prodN modulo prime.
                    //This value is a constant, so compute it here in order to facilitate 
                    //finding offsets later.

                    //solve prodN ^ -1 % p 
                    // slightly optimized modinv when prime >> prodN
                    inv = modinv_1c(prodN, prime);
                    t->sdata.root[i] = prime - inv;

                    //we can also speed things up by computing and storing the residue
                    //mod p of the first sieve location in the first residue class.  This provides
                    //a speedup by pulling this constant (involving a division) out of a critical loop
                    //when finding offsets of bucket sieved primes.
                    //these are only used by bucket sieved primes.			
                    //t->sdata.lower_mod_prime[i - t->sdata.bucket_start_id] = 
                    t->sdata.lower_mod_prime[i] =
                        mpz_tdiv_ui(tmpz, prime);
                }
            }
        }
        else 
        {
            t->command = SOE_COMPUTE_ROOTS;

#if defined(WIN32) || defined(_WIN64)
            SetEvent(t->run_event);
#else
            pthread_cond_signal(&t->run_cond);
            pthread_mutex_unlock(&t->run_lock);
#endif
        }		
    }


    //wait for each thread to finish
    for (i = 0; i < THREADS; i++) 
    {
        thread_soedata_t *t = thread_data + i;

        if (i < (THREADS - 1)) 
        {
#if defined(WIN32) || defined(_WIN64)
            WaitForSingleObject(t->finish_event, INFINITE);
#else
            pthread_mutex_lock(&t->run_lock);
            while (t->command != SOE_COMMAND_WAIT)
                pthread_cond_wait(&t->run_cond, &t->run_lock);
#endif
        }
    }

    //stop the worker threads
    for (i=0; i<THREADS - 1; i++)
    {
        stop_soe_worker_thread(thread_data + i);
    }

#endif

	gettimeofday(&tstop, NULL);

	difference = my_difftime(&tstart, &tstop);
	t = ((double)difference->secs + (double)difference->usecs / 1000000);
	free(difference);

    if (VFLAG > 2)
    {
        printf("time to compute bucket sieve roots = %1.2f\n", t);
    }

#ifdef INPLACE_BUCKET
	gettimeofday(&tstart, NULL);

	// inplace primes have special requirements because they operate on
	// the normal number line, and not in residue space
	for (; i < sdata->pboundi;i++)
	{
		uint64 starthit;
		uint32 startclass;
		uint64 startbit;
		uint32 rclass, bnum, rclassid;
		uint32 index = i - sdata->inplace_start_id;
		int a;

		// copy the prime into the special data structure
		//sdata->inplace_data[index].prime = sdata->sieve_p[i];

		// pull some computations involving a division out of the inner loop.
		// we need to know what prime/prodN and prime%prodN are.
		sdata->inplace_data[index].p_div = 
			sdata->sieve_p[i] / prodN;
		rclass = sdata->sieve_p[i] % prodN;
		rclassid = resID_mod30[rclass];
		sdata->inplace_data[index].p_mod = rclass;

		// now compute the starting hit in our sieve interval...
		starthit = (sdata->lowlimit / sdata->sieve_p[i] + 1) * sdata->sieve_p[i];

		// ... that is in one of our residue classes
		startclass = starthit % prodN;

		// using a lookup table
		startclass = next_mod30[rclassid][startclass];

		starthit += ((uint64)sdata->sieve_p[i] * (uint64)(startclass >> 8));
		startclass = startclass & 0xff;

		// the starting accumulated error is equal to the starting class
		sdata->inplace_data[index].eacc = startclass;

		// now compute the starting bit and block location for this starting hit
		startbit = (starthit - sdata->lowlimit - (uint64)startclass) / (uint64)prodN;

		// sanity check
		if (((starthit - sdata->lowlimit - (uint64)startclass) % (uint64)prodN) != 0)
			printf("starting bit is invalid!\n");

		sdata->inplace_data[index].bitloc = startbit & FLAGSIZEm1;
		bnum = startbit >> FLAGBITS;

		// finally, add the prime to a linked list
		// if the next hit is within our interval
		if (bnum < sdata->blocks)
		{
			//then reassign this prime to its next hit
			if (sdata->inplace_ptrs[bnum][resID_mod30[startclass]] == -1)
			{
				// this is the first hit in this block and rclass, so set the pointer
				// to this prime, and set next_pid = 0 so that we know to stop here
				// when we sieve
				sdata->inplace_ptrs[bnum][resID_mod30[startclass]] = index;
				sdata->inplace_data[index].next_pid = 0;
			}
			else
			{
				// add this prime to a listed list within the inplace sieve array.
				// this is done by first setting the next id to the current prime
				// at the end of the list
				sdata->inplace_data[index].next_pid = sdata->inplace_ptrs[bnum][resID_mod30[startclass]];

				// and then setting the end of the list to this prime
				sdata->inplace_ptrs[bnum][resID_mod30[startclass]] = index;						
			}
		}					

	}

	gettimeofday(&tstop, NULL);

	difference = my_difftime(&tstart, &tstop);
	t = ((double)difference->secs + (double)difference->usecs / 1000000);
	free(difference);

	if (VFLAG > 2)
		printf("time to compute inplace sieve roots = %1.2f\n", t);

#endif

	return;
}
