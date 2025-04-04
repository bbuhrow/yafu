/*----------------------------------------------------------------------
MIT License

Copyright (c) 2021 bbuhrow

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
----------------------------------------------------------------------*/

#include <stdint.h>
#include <stdio.h>
#include "threadpool.h"


#if (defined(_WIN32) || defined(_WIN64))  // && (!defined(__clang__))
DWORD WINAPI tpool_worker_main(LPVOID thread_data);
#else
#include <errno.h>
void *tpool_worker_main(void *thread_data);
#endif
void tpool_start(tpool_t *t);
void tpool_stop(tpool_t *t);

void tpool_start(tpool_t *t) 
{

    // create a thread and execute any user defined start function
    t->state = TPOOL_STATE_INIT;

    if (t->debug > 1)
    {
        printf("tpool: starting thread %d\n", t->tindex);
    }

    if (t->tpool_start_fcn != NULL)
    {
        (t->tpool_start_fcn)(t);
    }

#if (defined(_WIN32) || defined(_WIN64))  // && (!defined(__clang__))
    t->run_event = CreateEvent(NULL, FALSE, FALSE, NULL);
    t->finish_event = CreateEvent(NULL, FALSE, FALSE, NULL);
    *t->queue_event = CreateEvent(NULL, FALSE, FALSE, NULL);
    t->thread_id = CreateThread(NULL, 0, tpool_worker_main, t, 0, NULL);
    WaitForSingleObject(t->finish_event, INFINITE); /* wait for ready */
#else
    pthread_mutex_init(&t->run_lock, NULL);
    pthread_cond_init(&t->run_cond, NULL);

#ifdef USE_TPOOL_AFFINITY
    CPU_ZERO(t->cpus);
    CPU_SET(t->tindex, t->cpus);
    pthread_attr_setaffinity_np(t->attr, sizeof(cpu_set_t), t->cpus);

    pthread_create(&t->thread_id, t->attr, tpool_worker_main, t);
#else

    pthread_create(&t->thread_id, NULL, tpool_worker_main, t);
#endif

    pthread_mutex_lock(&t->run_lock); /* wait for ready */
    while (t->state != TPOOL_STATE_WAIT)
        pthread_cond_wait(&t->run_cond, &t->run_lock);
    pthread_mutex_unlock(&t->run_lock);
#endif

    if (t->debug > 1)
    {
        printf("tpool: thread %d started\n", t->tindex);
    }
}

void tpool_stop(tpool_t* t)
{
    t->state = TPOOL_STATE_END;

    if (t->debug > 1)
    {
        printf("tpool: stopping thread %d\n", t->tindex);
    }

    if (t->tpool_stop_fcn != NULL)
    {
        (t->tpool_stop_fcn)(t);
    }

#if (defined(_WIN32) || defined(_WIN64))  // && (!defined(__clang__))

    SetEvent(t->run_event);
    WaitForSingleObject(t->thread_id, INFINITE);
    CloseHandle(t->thread_id);
    CloseHandle(t->run_event);
    CloseHandle(t->finish_event);
    CloseHandle(*t->queue_event);
#else

    //pthread_mutex_lock(&t->run_lock);
    int count = 0;
    while (1)
    {
        int retval = pthread_mutex_trylock(&t->run_lock);

        if (retval == EDEADLK)
        {
            printf("deadlock in tpool_stop, thread %d attempted to re-lock mutex run_lock\n",
                t->tindex);
            printf("attempting to ignore problem...\n");
            return;
            //exit(10);
        }
        else if (retval == EBUSY)
        {
            count++;

            //printf("run_lock is busy in thread %d, (count = %d)\n",
            //    t->tindex, count);
#ifdef _MSC_VER
            // don't sleep and hope for the best?

#else
            usleep(1);
#endif

            if (count > 100)
            {
                printf("too many attempts to lock mutex run_lock in thread %d\n",
                    t->tindex);
                printf("attempting to ignore problem...\n");
                return;
            }
        }
        else
        {
            break;
        }
    }
    pthread_cond_signal(&t->run_cond);
    pthread_mutex_unlock(&t->run_lock);
    pthread_join(t->thread_id, NULL);
    pthread_cond_destroy(&t->run_cond);
    pthread_mutex_destroy(&t->run_lock);
#endif

    if (t->debug > 1)
    {
        printf("tpool: thread %d stopped\n", t->tindex);
    }

}


#if (defined(_WIN32) || defined(_WIN64))  // && (!defined(__clang__))
DWORD WINAPI tpool_worker_main(LPVOID thread_data) {
#else
void *tpool_worker_main(void *thread_data) {
#endif
    tpool_t *t = (tpool_t *)thread_data;

    /*
    * Respond to the master thread that we're ready for work. If we had any thread-
    * specific initialization which needed to be done, it would go before this signal.
    */
#if (defined(_WIN32) || defined(_WIN64))  // && (!defined(__clang__))
    t->state = TPOOL_STATE_WAIT;
    SetEvent(t->finish_event);
#else
    pthread_mutex_lock(&t->run_lock);
    t->state = TPOOL_STATE_WAIT;
    pthread_cond_signal(&t->run_cond);
    pthread_mutex_unlock(&t->run_lock);
    //omp_get_num_procs();
#endif

    while (1) {

        if (t->debug > 1)
        {
            printf("tpool: thread %d waiting for work\n", t->tindex);
        }

        /* wait forever for work to do */
#if (defined(_WIN32) || defined(_WIN64))  // && (!defined(__clang__))
        WaitForSingleObject(t->run_event, INFINITE);
#else
        pthread_mutex_lock(&t->run_lock);
        while (t->state == TPOOL_STATE_WAIT) {
            pthread_cond_wait(&t->run_cond, &t->run_lock);
        }
#endif

        /* do work */
        if (t->state == TPOOL_STATE_WORK)
        {
            if (t->debug > 1)
            {
                printf("tpool: thread %d executing work_fcn %d\n", t->tindex, t->work_fcn_id);
            }

            (t->tpool_work_fcn[t->work_fcn_id])(t);
        }
        else if (t->state == TPOOL_STATE_END)
        {
            if (t->debug > 1)
            {
                printf("tpool: thread %d state end\n", t->tindex);
            }

            break;
        }

        /* signal completion */
        t->state = TPOOL_STATE_WAIT;
#if (defined(_WIN32) || defined(_WIN64))  // && (!defined(__clang__))

        WaitForSingleObject(
            *t->queue_lock,    // handle to mutex
            INFINITE);  // no time-out interval

        t->thread_queue[(*(t->threads_waiting))++] = t->tindex;
        SetEvent(*t->queue_event);
        ReleaseMutex(*t->queue_lock);

#else
        pthread_mutex_unlock(&t->run_lock);

        // lock the work queue and insert my thread ID into it
        // this tells the master that my results should be synced
        // and I should be dispatched more work
        pthread_mutex_lock(t->queue_lock);
        t->thread_queue[(*(t->threads_waiting))++] = t->tindex;
        pthread_cond_signal(t->queue_cond);
        pthread_mutex_unlock(t->queue_lock);
#endif

    }

#if (defined(_WIN32) || defined(_WIN64))  // && (!defined(__clang__))
    return 0;
#else
    return NULL;
#endif
}


void tpool_go(tpool_t *thread_data)
{
    // thread work-queue controls
    int threads_working = 0;
    int *thread_queue;
    int *threads_waiting;
    int i;

#if (defined(_WIN32) || defined(_WIN64))  // && (!defined(__clang__))
    HANDLE queue_lock;
    HANDLE *queue_events = NULL;
#else
    pthread_mutex_t queue_lock;
    pthread_cond_t queue_cond;
    
#ifdef USE_TPOOL_AFFINITY
    pthread_attr_t attr;
    cpu_set_t cpus;
#endif
#endif

#ifdef USE_TPOOL_AFFINITY
    pthread_attr_init(&attr);
#endif

    // allocate the queue of threads waiting for work
    thread_queue = (int *)malloc(thread_data->num_threads * sizeof(int));
    threads_waiting = (int *)malloc(sizeof(int));

#if (defined(_WIN32) || defined(_WIN64))  // && (!defined(__clang__))
    queue_lock = CreateMutex(
        NULL,              // default security attributes
        FALSE,             // initially not owned
        NULL);             // unnamed mutex
    queue_events = (HANDLE *)malloc(thread_data->num_threads * sizeof(HANDLE));
#else
    pthread_mutex_init(&queue_lock, NULL);
    pthread_cond_init(&queue_cond, NULL);
#endif

    

    for (i = 0; i<thread_data->num_threads; i++)
    {
        if (thread_data[i].debug > 1)
        {
            printf("tpool: thread %d go init\n", thread_data[i].tindex);
        }

        thread_data[i].tindex = i;
        thread_data[i].tstartup = 1;
        // assign all threads a pointer to the waiting queue.  access to 
        // the array will be controlled by a mutex
        thread_data[i].thread_queue = thread_queue;
        thread_data[i].threads_waiting = threads_waiting;

#if (defined(_WIN32) || defined(_WIN64))  // && (!defined(__clang__))
        // assign a pointer to the mutex
        thread_data[i].queue_lock = &queue_lock;
        thread_data[i].queue_event = &queue_events[i];
#else
        thread_data[i].queue_lock = &queue_lock;
        thread_data[i].queue_cond = &queue_cond;

#ifdef USE_TPOOL_AFFINITY
        thread_data[i].attr = &attr;
        thread_data[i].cpus = &cpus;
#endif
#endif
    }

    // Activate the worker threads one at a time. 
    // Initialize the work queue to say all threads are waiting for work
    for (i = 0; i < thread_data->num_threads; i++)
    {
        tpool_start(thread_data + i);
        thread_queue[i] = i;        
    }

    *threads_waiting = thread_data->num_threads;

#if (defined(_WIN32) || defined(_WIN64))  // && (!defined(__clang__))
    // nothing
#else
    pthread_mutex_lock(&queue_lock);
#endif

    if (thread_data[0].debug > 1)
    {
        printf("tpool: main thread starting threadpool\n");
    }
    // we enter this loop with the main thread and num_threads other
    // threads running in parallel. This loop in the main thread is
    // responsible for calling the sync and dispatch functions as threads
    // in the pool finish, one at a time.  The other threads will write
    // to an array to indicate when they are done.
    while (1)
    {
        int tid = 0;

        // Process threads until there are no more waiting for their results to be collected
        while (*threads_waiting > 0)
        {            
            // Pop a waiting thread off the queue (OK, it's stack not a queue)
#if (defined(_WIN32) || defined(_WIN64))  // && (!defined(__clang__))
            WaitForSingleObject(
                queue_lock,    // handle to mutex
                INFINITE);  // no time-out interval
#endif

            tid = thread_queue[--(*threads_waiting)];

#if (defined(_WIN32) || defined(_WIN64))  // && (!defined(__clang__))
            ReleaseMutex(queue_lock);
#endif
            
            if (thread_data[tid].debug > 1)
            {
                printf("tpool: main thread sync %d\n", thread_data[tid].tindex);
            }
            
            if (thread_data[tid].tstartup == 0)
            {
                // if not in startup...
                if (thread_data[tid].tpool_sync_fcn != NULL)
                    (thread_data[tid].tpool_sync_fcn)(&thread_data[tid]);

                // this thread is done, so decrement the count of working threads
                threads_working--;
            }
            else
            {
                thread_data[tid].tstartup = 0;
            }

            // user dispatch function
            if (thread_data[tid].debug > 1)
            {
                printf("tpool: main thread dispatch to %d\n", thread_data[tid].tindex);
            }

            (thread_data[tid].tpool_dispatch_fcn)(&thread_data[tid]);

            // in the dispatch function the user will determine
            // what, if any, work function needs to be executed.
            if (thread_data[tid].work_fcn_id < thread_data[tid].num_work_fcn)
            {
                thread_data[tid].state = TPOOL_STATE_WORK;

                if (thread_data[tid].debug > 1)
                {
                    printf("tpool: main thread unlocking %d to execute work_fcn %d\n",
                        thread_data[tid].tindex, thread_data[tid].work_fcn_id);
                }

                // send the thread a signal to start processing the poly we just generated for it
#if (defined(_WIN32) || defined(_WIN64))  // && (!defined(__clang__))
                SetEvent(thread_data[tid].run_event);
#else
                pthread_mutex_lock(&thread_data[tid].run_lock);
                pthread_cond_signal(&thread_data[tid].run_cond);
                pthread_mutex_unlock(&thread_data[tid].run_lock);
#endif
                // this thread is now busy, so increment the count of working threads
                threads_working++;
            }
            else if (thread_data[tid].debug > 1)
            {
                printf("tpool: main thread no work for %d, now %d threads waiting and %d threads working\n",
                    thread_data[tid].tindex, *threads_waiting, threads_working);
            }

        }

        // if all threads are done, break out.
        // I believe this is the source of a possible race condition lock-up.
        // threads_working is modified in this main loop, so it is possible this
        // will be zero and we break out of this loop before the thread
        // actually stops (by seeing the TPOOL_STATE_END flag and exiting out of
        // its tpool_worker_main loop).
        if (threads_working == 0)
            break;

        if (thread_data[0].debug > 1)
        {
            printf("tpool: main thread entering wait state\n");
        }

        // wait for a thread to finish and put itself in the waiting queue
#if (defined(_WIN32) || defined(_WIN64))  // && (!defined(__clang__))
        WaitForMultipleObjects(
            thread_data[tid].num_threads,
            queue_events,
            FALSE,
            INFINITE);
#else
        pthread_cond_wait(&queue_cond, &queue_lock);
#endif

        if (thread_data[0].debug > 1)
        {
            printf("tpool: main thread leaving wait state\n");
        }
        
    }

    if (thread_data[0].debug > 1)
    {
        printf("tpool: main thread work finished\n");
    }

    //stop worker threads
    for (i = 0; i<thread_data[0].num_threads; i++)
    {
        tpool_stop(thread_data + i);
    }

    if (thread_data[0].debug > 1)
    {
        printf("tpool: main thread threadpool stopped\n");
    }

    free(thread_queue);
    free(threads_waiting);

#if (defined(_WIN32) || defined(_WIN64))  // && (!defined(__clang__))
    free(queue_events);
#endif


    return;
}

tpool_t * tpool_setup(int num_threads, void *start_fcn, void *stop_fcn, 
    void *sync_fcn, void *dispatch_fcn, void *udata)
{
    tpool_t *t = (tpool_t *)malloc(num_threads * sizeof(tpool_t));
    int i;

    //int numberOfProcessors = sysconf(_SC_NPROCESSORS_ONLN);
    //printf("Number of processors: %d\n", numberOfProcessors);

    for (i = 0; i < num_threads; i++)
    {
        t[i].num_threads = num_threads;
        t[i].tpool_dispatch_fcn = dispatch_fcn;
        t[i].tpool_start_fcn = start_fcn;
        t[i].tpool_stop_fcn = stop_fcn;
        t[i].tpool_sync_fcn = sync_fcn;
        t[i].num_work_fcn = 0;
        t[i].tindex = i;
        t[i].user_data = udata;
        t[i].debug = 0;
    }

    return t;
}

void tpool_add_work_fcn(tpool_t *tdata, void *work_fcn)
{
    int i;

    // to make it easier on the user, could add a "call by name" feature
    // where the user can pass in a name associated with this function.
    // then in the dispatch function they can set a name instead of work_fcn_id;
    for (i = 0; i < tdata->num_threads; i++)
    {
        tdata[i].num_work_fcn++;
        tdata[i].tpool_work_fcn[tdata[i].num_work_fcn - 1] = work_fcn;
    }

    return;
}

