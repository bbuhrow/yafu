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

#include "yafu_ecm.h"
#include "factor.h"
#include "yafu.h"
#include "calc.h"

#if defined(FORK_ECM)
int ecm_loop(fact_obj_t *fobj)
{
	//thread data holds all data needed during sieving
	ecm_thread_data_t *thread_data, *my_thread_data;
	z *n = &fobj->ecm_obj.n;
	z * return_factor, *my_return_factor;
    int *my_curves_run;
    int *curves_run, *factor_found;
	int *shmid_list, num_shmids=0;
    int master_thread;
	int total_curves_run;
	int numcurves = fobj->ecm_obj.num_curves;
	z d,t,nn;
	FILE *flog;
	int i,j;
	//maybe make this an input option: whether or not to stop after
	//finding a factor in the middle of running a requested batch of curves
	int bail_on_factor = 1;
	int bail = 0;
	int charcount = 0, charcount2 = 0;
	int input_digits = ndigits(n);

	if (ecm_check_input(fobj) == 0)
		return 0;

	//ok, having gotten this far we are now ready to run the requested
	//curves.  initialize the needed data structures, then split
	//the curves up over N threads.  The main thread will farm out different
	//sigmas to the worker threads.

	//initialize the flag to watch for interrupts, and set the
	//pointer to the function to call if we see a user interrupt
	ECM_ABORT = 0;
	signal(SIGINT,ecmexit);

	//init ecm process
	ecm_process_init(fobj);

	//init local big ints
	zInit(&d);
	zInit(&t);
	zInit(&nn);

	// Allocate room for a list of all shared memory segments -- one per thread plus 4 more
	// We need this to correctly destroy all the segments on exit
	shmid_list = (int *)malloc((THREADS + 4) * sizeof(int));

	// Allocate shared memory for a variable saying someone found a factor, and an
    // array for each thread to record the number of curves it has run.
    factor_found = (int *)malloc_shared(sizeof(int), &shmid_list[num_shmids++]);
    curves_run = (int *)malloc_shared(THREADS * sizeof(int), &shmid_list[num_shmids++]);

    // We make the thread_data structures shared so that the master thread can
    // log a few pieces of information (e.g. the sigma and stage where the factor
    // was found) correctly. Note that many things in this structure are pointers
    // which are later dynamically allocated, and those allocations will not be
    // shared. The master thread must be careful to access only the statically
    // declared fields in other threads' structures.
    thread_data = malloc_shared(THREADS * sizeof(ecm_thread_data_t), &shmid_list[num_shmids++]);

    // Allocate shared memory for each thread to return the factor it found.
    // I don't know how to support dynamic resizing of shared memory, so I
    // just allocate plenty of space and assume it will never overflow.
    return_factor = malloc_shared(THREADS * sizeof(z), &shmid_list[num_shmids++]);
    for (i=0; i<THREADS; i++) {
        return_factor[i].size = 1;
        return_factor[i].alloc = 4096;
        return_factor[i].type = UNKNOWN;
        return_factor[i].val = (fp_digit *)malloc_shared(4096 * sizeof(fp_digit), &shmid_list[num_shmids++]);
    }

    // Set things up for the master thread
    master_thread = 1;
    my_return_factor = &return_factor[THREADS-1];
    my_curves_run = &curves_run[THREADS-1];
    my_thread_data = &thread_data[THREADS-1];

	// Divide the curves among the threads, rounding up so that all threads get the
    // same number of curves. Also keep track of the total number of curves for logging.
    numcurves = (numcurves + THREADS - 1) / THREADS;
    fobj->ecm_obj.num_curves = numcurves * THREADS;

	/* activate the threads one at a time. The last is the
	   master thread (i.e. not a thread at all). */
    for (i=0; i < THREADS-1; i++)
    {
        // We need to ensure that each child has different RNG state, otherwise they'll
        // all get the same stream of random numbers. So we have the master thread generate
        // a number and then copy that into the internal RNG state of the child.
        uint64 child_lcgstate = spRand(2,MAX_DIGIT);

        int pid = fork();
        if (pid == 0) {
            // I am a child thread. Set up pointers to my shared buffers accordingly.
            master_thread = 0;
            my_return_factor = &return_factor[i];
            my_curves_run = &curves_run[i];
            my_thread_data = &thread_data[i];
			my_thread_data->thread_num = i;
			my_thread_data->curves_run = *my_curves_run;
            LCGSTATE = child_lcgstate;			
            break;
        }
		my_thread_data->thread_num = i;
		my_thread_data->curves_run = *my_curves_run;
    }

    // Now we have THREADS number of separate processes running
	// Initialize the internal pieces of thread_data
	my_thread_data->fobj = fobj;	
    ecm_thread_init(my_thread_data);

	// TBD: make this a subroutine
	if (master_thread && VFLAG >= 0)
	{
		for (i=0;i<charcount+charcount2;i++)
			printf("\b");

        for (i=0, total_curves_run=0; i<THREADS; i++)
            total_curves_run += curves_run[i];

		charcount = printf("ecm: %d/%d curves on C%d input, at ",
			total_curves_run, fobj->ecm_obj.num_curves, input_digits);
		charcount2 = print_B1B2(fobj, NULL);
		fflush(stdout);
	}

	// Now each thread runs its assigned number of curves
	for (j=0; j < numcurves; j++)
	{
		//watch for an abort
		if (ECM_ABORT)
		{
            // On an abort, the master thread will wait for all the children to
            // complete and then log any found factors. The children just exit
            // right away.
            if (master_thread) 
            {
                for (i=0; i<THREADS-1; i++)
                    wait(NULL);
                print_factors(fobj);
            }

			// We still have to mark shared memory segments for destruction
            destroy_shm_segments(shmid_list, num_shmids);

            exit(1);
		}

		ecm_get_sigma(my_thread_data);

        // Run a curve!
        ecm_do_one_curve(my_thread_data);
        (*my_curves_run)++;
		my_thread_data->curves_run = *my_curves_run;

        if (zCompare(&my_thread_data->factor,&zOne) > 0 && zCompare(&my_thread_data->factor,n) < 0) {

            // non-trivial factor found
            // copy it into the shared return buffer so the master thread can see it, 
            // and set the shared flag which tells everyone the game is over.
            zCopy(&my_thread_data->factor, my_return_factor);
            *factor_found = 1;

            //printf("thread %d found a factor\n", my_thread_data->thread_num);
		}

        if (*factor_found) {

            // If someone found a factor, the child threads just exit
            if (!master_thread)
                exit(0);

            // The master thread waits for all the children to complete
            for (i=0; i<THREADS-1; i++)
				wait(NULL);

            // Now that all the children are done, look for and log any factors they found
            for (i=0; i<THREADS; i++) {

                // Check thread i's return buffer. If it doesn't contain a factor, just move on
                if (!(zCompare(&return_factor[i],&zOne) > 0 && zCompare(&return_factor[i],n) < 0))
					continue;

				//since we could be doing many curves in parallel,
				//it's possible more than one curve found this factor at the
				//same time.  We divide out factors as we find them, so
				//just check if this factor still divides n
				zCopy(n,&nn);
				zDiv(&nn,&return_factor[i],&t,&d);
				if (zCompare(&d,&zZero) == 0)
				{
					zCopy(&t, n);

					//yes, it does... proceed to record the factor
					zCopy(&return_factor[i], &thread_data[i].factor);
					ecm_deal_with_factor(&thread_data[i]);			

					// factor will have been printed with a newline - reset the backspaces.
					if (VFLAG > 0)
						charcount = charcount2 = 0;

					//we found a factor and might want to stop,
					//but should check all threads output first
					if (bail_on_factor)
						bail = 1;
					else if (isPrime(n))
						bail = 1;
					else
					{
                        // Note: this path no longer works since the child threads have exited.

						//found a factor and the cofactor is composite.
						//the user has specified to keep going with ECM until the 
						//curve counts are finished thus:
						//we need to re-initialize with a different modulus.  this is
						//independant of the thread data initialization
						ecm_process_free(fobj);
						ecm_process_init(fobj);
					}
				}
			}

        }

        // Note that bail can only be set on the master thread, since if a factor was
        // found all other threads will already have exited.
        if (bail)
			goto done;

        // Only the master thread prints updates, otherwise they all step on each other.
        if (master_thread && VFLAG >= 0)
        {
			for (i=0;i<charcount+charcount2;i++)
			printf("\b");

            for (i=0, total_curves_run=0; i<THREADS; i++)
				total_curves_run += curves_run[i];

			charcount = printf("ecm: %d/%d curves on C%d input, at ",
				total_curves_run, fobj->ecm_obj.num_curves, input_digits);
			charcount2 = print_B1B2(fobj, NULL);
			fflush(stdout);
        }


	}

done:

    // Either a factor was found, or this thread has completed all of its requested curves.
    // At this point the child threads can just exit.
    if (!master_thread)
        exit(0);

    // The master thread waits for all children to complete
    for (i=0; i<THREADS-1; i++)
        wait(NULL);

    // print status message again so it shows any work done after the master thread finished its loop
    if (VFLAG >= 0)
    {
        for (i=0;i<charcount+charcount2;i++)
            printf("\b");

        for (i=0, total_curves_run=0; i<THREADS; i++)
            total_curves_run += curves_run[i];

        charcount = printf("ecm: %d/%d curves on C%d input, at ",
			total_curves_run, fobj->ecm_obj.num_curves, input_digits);
        charcount2 = print_B1B2(fobj, NULL);
        fflush(stdout);
    }

	if (VFLAG >= 0)
		printf("\n");

	flog = fopen(fobj->logname,"a");
	if (flog == NULL)
	{
		printf("could not open %s for writing\n",fobj->logname);
		fclose(flog);
		return 0;
	}

    for (i=0, total_curves_run=0; i<THREADS; i++)
		total_curves_run += curves_run[i];

	logprint(flog,"Finished %d curves using Lenstra ECM method on C%d input, ",
		total_curves_run,input_digits);
	print_B1B2(fobj, flog);
	fprintf(flog, "\n");

	fclose(flog);

    // Detach from all shared memory segments allocated above
    ecm_thread_free(my_thread_data);
    for (i=0; i<THREADS; i++) {
        shmdt(return_factor[i].val);
    }
    shmdt(return_factor);
    shmdt(factor_found);
    shmdt(curves_run);
    shmdt(thread_data);

	// And mark them for destruction
    destroy_shm_segments(shmid_list, num_shmids);
    free(shmid_list);

	zFree(&d);
	zFree(&t);
	zFree(&nn);
	signal(SIGINT,NULL);
	
	ecm_process_free(fobj);

	return total_curves_run;
}

#else


int ecm_loop(fact_obj_t *fobj)
{
	ecm_thread_data_t *thread_data;		//an array of thread data objects
	z *n = &fobj->ecm_obj.n;
	z d,t,nn;
	FILE *flog;
	int i,j;
	//maybe make this an input option: whether or not to stop after
	//finding a factor in the middle of running a requested batch of curves
	int total_curves_run;
	int bail_on_factor = 1;
	int bail = 0;
	int charcount = 0, charcount2 = 0;
	int input_digits = ndigits(n);

	if (ecm_check_input(fobj) == 0)
		return 0;

	//ok, having gotten this far we are now ready to run the requested
	//curves.  initialize the needed data structures, then split
	//the curves up over N threads.  The main thread will farm out different
	//sigmas to the worker threads.

	//initialize the flag to watch for interrupts, and set the
	//pointer to the function to call if we see a user interrupt
	ECM_ABORT = 0;
	signal(SIGINT,ecmexit);

	//init ecm process
	ecm_process_init(fobj);

	//init local big ints
	zInit(&d);
	zInit(&t);
	zInit(&nn);

	thread_data = (ecm_thread_data_t *)malloc(THREADS * sizeof(ecm_thread_data_t));
	for (i=0; i<THREADS; i++)
	{
		thread_data[i].fobj = fobj;
		thread_data[i].thread_num = i;
		ecm_thread_init(&thread_data[i]);
	}

	//round numcurves up so that each thread has something to do every iteration.
	//this prevents a single threaded "cleanup" round after the multi-threaded rounds
	//to finish the leftover requested curves.
	fobj->ecm_obj.num_curves += ((fobj->ecm_obj.num_curves % THREADS) ? 
		(THREADS - (fobj->ecm_obj.num_curves % THREADS)) : 0);

	if (VFLAG >= 0)
	{
		for (i=0;i<charcount+charcount2;i++)
			printf("\b");

		for (i=0, total_curves_run=0; i<THREADS; i++)
            total_curves_run += thread_data[i].curves_run;

		charcount = printf("ecm: %d/%d curves on C%d input, at ",
			total_curves_run, fobj->ecm_obj.num_curves, ndigits(n));
		charcount2 = print_B1B2(fobj, NULL);
		fflush(stdout);
	}

	/* activate the threads one at a time. The last is the
	   master thread (i.e. not a thread at all). */

	for (i = 0; i < THREADS - 1; i++)
		ecm_start_worker_thread(thread_data + i, 0);

	ecm_start_worker_thread(thread_data + i, 1);

	//split the requested curves up among the specified number of threads. 
	for (j=0; j < fobj->ecm_obj.num_curves / THREADS; j++)
	{
		//watch for an abort
		if (ECM_ABORT)
		{
			print_factors(fobj);
			exit(1);
		}

		//do work on different sigmas
		for (i=0; i<THREADS; i++)
		{
			ecm_get_sigma(&thread_data[i]);

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
				zDiv(&nn,&thread_data[i].factor,&t, &d);
				if (zCompare(&d,&zZero) == 0)
				{
					//yes, it does... proceed to record the factor
					zCopy(&t, n);
					ecm_deal_with_factor(&thread_data[i]);
					
					// factor will have been printed with a newline - reset the backspaces.
					if (VFLAG > 0)
						charcount = charcount2 = 0;

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
						ecm_process_free(fobj);
						ecm_process_init(fobj);
					}
				}
			}

			thread_data[i].curves_run++;
		}

		if (bail)
			goto done;

		if (VFLAG >= 0)
		{
			for (i=0;i<charcount+charcount2;i++)
				printf("\b");

			for (i=0, total_curves_run=0; i<THREADS; i++)
				total_curves_run += thread_data[i].curves_run;

			charcount = printf("ecm: %d/%d curves on C%d input, at ",
				total_curves_run, fobj->ecm_obj.num_curves, ndigits(n));

			charcount2 = print_B1B2(fobj, NULL);
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

	for (i=0, total_curves_run=0; i<THREADS; i++)
		total_curves_run += thread_data[i].curves_run;

	logprint(flog,"Finished %d curves using Lenstra ECM method on C%d input, ",
		total_curves_run,input_digits);

	print_B1B2(fobj, flog);
	fprintf(flog, "\n");

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
	signal(SIGINT,NULL);
	
	ecm_process_free(fobj);

	return total_curves_run;
}

#endif

int ecm_deal_with_factor(ecm_thread_data_t *thread_data)
{
	FILE *flog;
	fact_obj_t *fobj = thread_data->fobj;
	z *factor = &thread_data->factor;
	z *n = &thread_data->n;
	int curves_run = thread_data->curves_run;
	int thread_num = thread_data->thread_num;
	z tmp1, tmp2;

	zInit(&tmp1);
	zInit(&tmp2);

	flog = fopen(fobj->logname,"a");
	if (flog == NULL)
	{
		printf("could not open %s for writing\n",fobj->logname);
		return 0;
	}

	if (isPrime(factor))
	{
		factor->type = PRP;
		add_to_factor_list(fobj, factor);

		if (VFLAG > 0)
			printf("\necm: found prp%d factor = %s\n",ndigits(factor),z2decstr(factor,&gstr1));

		logprint(flog,"prp%d = %s (curve %d stg%d B1=%u sigma=%u thread=%d)\n",
			ndigits(factor),
			z2decstr(factor,&gstr1),
			curves_run+1, thread_data->stagefound,
			fobj->ecm_obj.B1, thread_data->sigma, thread_num);
	}
	else
	{
		factor->type = COMPOSITE;
		add_to_factor_list(fobj, factor);
		if (VFLAG > 0)
			printf("\necm: found c%d factor = %s\n",ndigits(factor),z2decstr(factor,&gstr1));

		logprint(flog,"c%d = %s (curve %d stg%d B1=%u sigma=%u thread=%d)\n",
			ndigits(factor),
			z2decstr(factor,&gstr1),
			curves_run+1, thread_data->stagefound,
			fobj->ecm_obj.B1, thread_data->sigma, thread_num);
	}

	fclose(flog);
			
	//reduce input
	zDiv(n,factor,&tmp1,&tmp2);
	zCopy(&tmp1,n);
	zCopy(n,&thread_data->n);

	zFree(&tmp1);
	zFree(&tmp2);
	return 1;
}

int ecm_get_sigma(ecm_thread_data_t *thread_data)
{
	fact_obj_t *fobj = thread_data->fobj;
	z tmp;

	zInit(&tmp);
	if (fobj->ecm_obj.sigma != 0)
	{
		if (thread_data->curves_run == 0 &&  fobj->ecm_obj.num_curves > 1)
			printf("WARNING: work will be duplicated with sigma fixed and numcurves > 1\n");
		thread_data->sigma = fobj->ecm_obj.sigma;
	}
	else if (get_uvar("sigma",&tmp))
	{
		thread_data->sigma = spRand(6,MAX_DIGIT);
	}
	else
	{
		if (thread_data->curves_run == 0 &&  fobj->ecm_obj.num_curves > 1)
			printf("WARNING: work will be duplicated with sigma fixed and numcurves > 1\n");
		thread_data->sigma = (uint32)tmp.val[0];
	}

	zFree(&tmp);
	return 0;
}

int ecm_check_input(fact_obj_t *fobj)
{
	FILE *flog;
	z *n = &fobj->ecm_obj.n;

	//open the log file and annouce we are starting ECM
	flog = fopen(fobj->logname,"a");
	if (flog == NULL)
	{
		printf("could not open %s for writing\n",fobj->logname);
		fclose(flog);
		return 0;
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
		add_to_factor_list(fobj, &zThree);
		logprint(flog,"Trivial factor of 3 found in ECM\n");
		fclose(flog);
		return 0;
	}
	if ((n->val[0] & 0x1) != 1)
	{
		zShortDiv(n,2,n);
		add_to_factor_list(fobj, &zTwo);
		logprint(flog,"Trivial factor of 2 found in ECM\n");
		fclose(flog);
		return 0;
	}
	if (isPrime(n))
	{
		//maybe have an input flag to optionally not perform
		//PRP testing (useful for really big inputs)
		n->type = PRP;
		add_to_factor_list(fobj, n);
		logprint(flog,"prp%d = %s\n",ndigits(n),z2decstr(n,&gstr1));
		zCopy(&zOne,n);
		fclose(flog);
		return 0;
	}

	//close the log file for until we have something further to report
	fclose(flog);

	return 1;
}

void ecm_start_worker_thread(ecm_thread_data_t *t, uint32 is_master_thread) {

	/* create a thread that will process a polynomial The last poly does 
	   not get its own thread (the current thread handles it) */

	if (is_master_thread) {
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

void ecm_stop_worker_thread(ecm_thread_data_t *t, uint32 is_master_thread)
{
	if (is_master_thread) {
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

#if defined(WIN32) || defined(_WIN64)
	return 0;
#else
	return NULL;
#endif
}


#if defined(FORK_ECM)
void *malloc_shared(size_t bytes, int *save_shmid)
{
    int shmid = shmget(IPC_PRIVATE, bytes, SHM_R | SHM_W);
    if (shmid == -1) {
        printf("Couldn't allocated shared memory segment in ECM\n");
        exit(1);
    }
	*save_shmid = shmid;
    return shmat(shmid, 0, 0);
}

void destroy_shm_segments(int *shmid_list, int num_shmids)
{
    int i;

    for (i=0; i<num_shmids; i++)
		shmctl(shmid_list[i], IPC_RMID, NULL);
}
#endif

void ecm_process_init(fact_obj_t *fobj)
{
	//initialize things which all threads will need when using
	//GMP-ECM
	TMP_THREADS = THREADS;
	TMP_STG2_MAX = fobj->ecm_obj.B2;

#if defined(FORK_ECM)
	// For small curves, or if we're only running one curve, don't
    // bother spawning multiple threads
    if (fobj->ecm_obj.B1 < 10000 || fobj->ecm_obj.num_curves == 1)
		THREADS = 1;
#else
		
	if (THREADS > 1)
	{
		if (VFLAG >= 2)
			printf("GMP-ECM does not support multiple threads... running single threaded\n");
		THREADS = 1;
	}
#endif

	return;
}

void ecm_thread_init(ecm_thread_data_t *tdata)
{
	//initialize everything for all threads using GMP-ECM
	zInit(&tdata->n);
	zInit(&tdata->factor);
	mpz_init(tdata->gmp_n);
	mpz_init(tdata->gmp_factor);
	ecm_init(tdata->params);
	gmp_randseed_ui(tdata->params->rng, get_rand(&g_rand.low, &g_rand.hi));
	zCopy(&tdata->fobj->ecm_obj.n, &tdata->n);
	tdata->params->method = ECM_ECM;
	tdata->curves_run = 0;
		
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

void ecm_process_free(fact_obj_t *fobj)
{
	THREADS = TMP_THREADS;
	fobj->ecm_obj.B2 = TMP_STG2_MAX;
	return;
}

void *ecm_do_one_curve(void *ptr)
{
	int status;
#if defined(_WIN64) && BITS_PER_DIGIT == 32
	size_t count;
#endif

	//unpack the data structure and stuff inside it
	ecm_thread_data_t *thread_data = (ecm_thread_data_t *)ptr;
	fact_obj_t *fobj = thread_data->fobj;

	thread_data->params->B1done = 1.0 + floor (1 * 128.) / 134217728.;
	if (VFLAG >= 3)
		thread_data->params->verbose = VFLAG - 2;		
	mpz_set_ui(thread_data->params->x, (unsigned long)0);
	mpz_set_ui(thread_data->params->sigma, thread_data->sigma);

#if defined(_WIN64) && BITS_PER_DIGIT == 32
	mpz_import(thread_data->gmp_n, (size_t)(abs(thread_data->n.size)), -1, sizeof(uint32), 
		0, (size_t)0, thread_data->n.val);
#else
	//wrapper for YAFU bigints and call to gmp-ecm
	mp2gmp(&thread_data->n, thread_data->gmp_n);
#endif

	if (fobj->ecm_obj.stg2_is_default == 0)
	{
		//not default, tell gmp-ecm to use the requested B2
		//printf("using requested B2 value\n");
		sp642z(fobj->ecm_obj.B2,&thread_data->factor);
		mp2gmp(&thread_data->factor,thread_data->params->B2);
		zClear(&thread_data->factor);
	}

	status = ecm_factor(thread_data->gmp_factor, thread_data->gmp_n,
			fobj->ecm_obj.B1, thread_data->params);

#if defined(_WIN64) && BITS_PER_DIGIT == 32
	zClear(&thread_data->n);
	mpz_export(thread_data->n.val, &count, -1, sizeof(uint32),
			0, (size_t)0, thread_data->gmp_n);
	thread_data->n.size = count;
#else
	//update n: not sure if gmp-ecm modifies it
	gmp2mp(thread_data->gmp_n,&thread_data->n);
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

// function definitions
void ecmexit(int sig)
{
	printf("\nAborting...\n");
	ECM_ABORT = 1;
	return;
}

int print_B1B2(fact_obj_t *fobj, FILE *fid)
{
	int i;
	char suffix;
	char stg1str[20];
	char stg2str[20];

	if (fobj->ecm_obj.B1 % 1000000000 == 0)
	{
		suffix = 'B';
		sprintf(stg1str,"%u%c",fobj->ecm_obj.B1 / 1000000000, suffix);
	}
	else if (fobj->ecm_obj.B1 % 1000000 == 0)
	{
		suffix = 'M';
		sprintf(stg1str,"%u%c",fobj->ecm_obj.B1 / 1000000, suffix);
	}
	else if (fobj->ecm_obj.B1 % 1000 == 0)
	{
		suffix = 'K';
		sprintf(stg1str,"%u%c",fobj->ecm_obj.B1 / 1000, suffix);
	}
	else
	{
		sprintf(stg1str,"%u",fobj->ecm_obj.B1);
	}

	if (fobj->ecm_obj.stg2_is_default == 0)
	{
		if (fobj->ecm_obj.B2 % 1000000000 == 0)
		{
			suffix = 'B';
			sprintf(stg2str,"%" PRIu64 "%c",fobj->ecm_obj.B2 / 1000000000, suffix);
		}
		else if (fobj->ecm_obj.B2 % 1000000 == 0)
		{
			suffix = 'M';
			sprintf(stg2str,"%" PRIu64 "%c",fobj->ecm_obj.B2 / 1000000, suffix);
		}
		else if (fobj->ecm_obj.B2 % 1000 == 0)
		{
			suffix = 'K';
			sprintf(stg2str,"%" PRIu64 "%c",fobj->ecm_obj.B2 / 1000, suffix);
		}
		else
		{
			sprintf(stg2str,"%" PRIu64 "",fobj->ecm_obj.B2);
		}
	}
	else
		sprintf(stg2str, "gmp-ecm default");

	if (fid == NULL)
		i = printf("B1 = %s, B2 = %s",stg1str,stg2str);
	else
		i = fprintf(fid,"B1 = %s, B2 = %s",stg1str,stg2str);

	return i;
}



