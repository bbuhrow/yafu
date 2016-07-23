#include "threadpool.h"
#include <stdint.h>

#if defined(WIN32) || defined(_WIN64)
DWORD WINAPI tpool_worker_main(LPVOID thread_data);
#else
void *tpool_worker_main(void *thread_data);
#endif
void tpool_start(tpool_t *t);
void tpool_stop(tpool_t *t);

void tpool_start(tpool_t *t) 
{

    // create a thread and execute any user defined start function
    t->state = TPOOL_STATE_INIT;

    if (t->tpool_start_fcn != NULL)
    {
        (t->tpool_start_fcn)(t);
    }

#if defined(WIN32) || defined(_WIN64)
    t->run_event = CreateEvent(NULL, FALSE, FALSE, NULL);
    t->finish_event = CreateEvent(NULL, FALSE, FALSE, NULL);
    *t->queue_event = CreateEvent(NULL, FALSE, FALSE, NULL);
    t->thread_id = CreateThread(NULL, 0, tpool_worker_main, t, 0, NULL);
    WaitForSingleObject(t->finish_event, INFINITE); /* wait for ready */
#else
    pthread_mutex_init(&t->run_lock, NULL);
    pthread_cond_init(&t->run_cond, NULL);

    pthread_create(&t->thread_id, NULL, tpool_worker_main, t);

    pthread_mutex_lock(&t->run_lock); /* wait for ready */
    while (t->state != TPOOL_STATE_WAIT)
        pthread_cond_wait(&t->run_cond, &t->run_lock);
    pthread_mutex_unlock(&t->run_lock);
#endif
}

void tpool_stop(tpool_t *t)
{
    t->state = TPOOL_STATE_END;

    if (t->tpool_stop_fcn != NULL)
    {
        (t->tpool_stop_fcn)(t);
    }

#if defined(WIN32) || defined(_WIN64)
    
    SetEvent(t->run_event);
    WaitForSingleObject(t->thread_id, INFINITE);
    CloseHandle(t->thread_id);
    CloseHandle(t->run_event);
    CloseHandle(t->finish_event);
    CloseHandle(*t->queue_event);
#else
    pthread_mutex_lock(&t->run_lock);
    pthread_cond_signal(&t->run_cond);
    pthread_mutex_unlock(&t->run_lock);
    pthread_join(t->thread_id, NULL);
    pthread_cond_destroy(&t->run_cond);
    pthread_mutex_destroy(&t->run_lock);
#endif
}


#if defined(WIN32) || defined(_WIN64)
DWORD WINAPI tpool_worker_main(LPVOID thread_data) {
#else
void *tpool_worker_main(void *thread_data) {
#endif
    tpool_t *t = (tpool_t *)thread_data;

    /*
    * Respond to the master thread that we're ready for work. If we had any thread-
    * specific initialization which needed to be done, it would go before this signal.
    */
#if defined(WIN32) || defined(_WIN64)
    t->state = TPOOL_STATE_WAIT;
    SetEvent(t->finish_event);
#else
    pthread_mutex_lock(&t->run_lock);
    t->state = TPOOL_STATE_WAIT;
    pthread_cond_signal(&t->run_cond);
    pthread_mutex_unlock(&t->run_lock);
#endif

    while (1) {

        /* wait forever for work to do */
#if defined(WIN32) || defined(_WIN64)
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
            (t->tpool_work_fcn[t->work_fcn_id])(t);
        }
        else if (t->state == TPOOL_STATE_END)
        {
            break;
        }

        /* signal completion */
        t->state = TPOOL_STATE_WAIT;
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


void tpool_go(tpool_t *thread_data)
{
    // thread work-queue controls
    int threads_working = 0;
    int *thread_queue;
    int *threads_waiting;
    int i;

#if defined(WIN32) || defined(_WIN64)
    HANDLE queue_lock;
    HANDLE *queue_events = NULL;
#else
    pthread_mutex_t queue_lock;
    pthread_cond_t queue_cond;
#endif

    // allocate the queue of threads waiting for work
    thread_queue = (int *)malloc(thread_data->num_threads * sizeof(int));
    threads_waiting = (int *)malloc(sizeof(int));

#if defined(WIN32) || defined(_WIN64)
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
        thread_data[i].tindex = i;
        thread_data[i].tstartup = 1;
        // assign all thread's a pointer to the waiting queue.  access to 
        // the array will be controlled by a mutex
        thread_data[i].thread_queue = thread_queue;
        thread_data[i].threads_waiting = threads_waiting;

#if defined(WIN32) || defined(_WIN64)
        // assign a pointer to the mutex
        thread_data[i].queue_lock = &queue_lock;
        thread_data[i].queue_event = &queue_events[i];
#else
        thread_data[i].queue_lock = &queue_lock;
        thread_data[i].queue_cond = &queue_cond;
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

#if defined(WIN32) || defined(_WIN64)
    // nothing
#else
    pthread_mutex_lock(&queue_lock);
#endif

    printf("=== starting threadpool\n");
    while (1)
    {
        int tid;

        // Process threads until there are no more waiting for their results to be collected
        while (*threads_waiting > 0)
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

            // if not in startup...
            if (thread_data[tid].tstartup == 0)
            {
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
            (thread_data[tid].tpool_dispatch_fcn)(&thread_data[tid]);

            // in the dispatch function the user will determine
            // what, if any, work function needs to be executed.
            if (thread_data[tid].work_fcn_id < thread_data[tid].num_work_fcn)
            {
                thread_data[tid].state = TPOOL_STATE_WORK;
                // send the thread a signal to start processing the poly we just generated for it
#if defined(WIN32) || defined(_WIN64)
                SetEvent(thread_data[tid].run_event);
#else
                pthread_mutex_lock(&thread_data[tid].run_lock);
                pthread_cond_signal(&thread_data[tid].run_cond);
                pthread_mutex_unlock(&thread_data[tid].run_lock);
#endif

                // this thread is now busy, so increment the count of working threads
                threads_working++;
            }

        } // while (*threads_waiting > 0)

        // if all threads are done, break out
        if (threads_working == 0)
            break;

        // wait for a thread to finish and put itself in the waiting queue
#if defined(WIN32) || defined(_WIN64)
        j = WaitForMultipleObjects(
            thread_data[tid].num_threads,
            queue_events,
            FALSE,
            INFINITE);
#else
        pthread_cond_wait(&queue_cond, &queue_lock);
#endif
        
    }

    printf("=== work finished\n");

    //stop worker threads
    for (i = 0; i<thread_data[0].num_threads; i++)
    {
        tpool_stop(thread_data + i);
    }

    printf("=== threadpool stopped\n");

    free(thread_queue);
    free(threads_waiting);

#if defined(WIN32) || defined(_WIN64)
    free(queue_events);
#endif


    return;
}

tpool_t * tpool_setup(int num_threads, void *start_fcn, void *stop_fcn, 
    void *sync_fcn, void *dispatch_fcn, void *udata)
{
    tpool_t *t = (tpool_t *)malloc(num_threads * sizeof(tpool_t));
    int i;

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
