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

#ifndef _THREADPOOL_H
#define _THREADPOOL_H


#include <stdlib.h>
#include <malloc.h>

#ifdef __cplusplus
extern "C" {
#endif

//#define USE_TPOOL_AFFINITY

#ifdef USE_TPOOL_AFFINITY
#define __USE_GNU
#include <sched.h>
#include <unistd.h>
#endif

#if (defined(_WIN32) || defined(_WIN64))  // && (!defined(__clang__))
    // MINGW defines these too: mingw will use windows 
    // threading instead of pthreads.
#include <windows.h>
#include <process.h>

#else /* !WIN32 */

#include <pthread.h>

#endif /* WIN32 */



enum tpool_state {
    TPOOL_STATE_INIT,
    TPOOL_STATE_WAIT,
    TPOOL_STATE_WORK,
    TPOOL_STATE_END
};

typedef struct
{
    int tindex;
    int tstartup;
    int num_threads;
    int done;
    int debug;

    volatile enum tpool_state state;

    int num_work_fcn;   // number of work functions
    int work_fcn_id;   // index of work function to execute

    // pointers to user-defined functions to execute during
    // the various states of the threadpool.
    void(*tpool_start_fcn)(void *);
    void(*tpool_stop_fcn)(void *);
    void(*tpool_sync_fcn)(void *);
    void(*tpool_dispatch_fcn)(void *);
    void(*tpool_work_fcn[16])(void *);

    // a user-defined structure to carry data throughout the threadpool.
    // user will have to cast this to their appropriate type before use.
    void *user_data;

    /* fields for thread pool synchronization */
    volatile int *thread_queue;
    volatile int *threads_waiting;

#if (defined(_WIN32) || defined(_WIN64))  // && (!defined(__clang__))
    HANDLE thread_id;
    HANDLE run_event;

    HANDLE finish_event;
    HANDLE *queue_event;
    HANDLE *queue_lock;

#else
    pthread_t thread_id;
    pthread_mutex_t run_lock;
    pthread_cond_t run_cond;

    pthread_mutex_t *queue_lock;
    pthread_cond_t *queue_cond;

#ifdef USE_TPOOL_AFFINITY
    pthread_attr_t *attr;
    cpu_set_t *cpus;
#endif

#endif

} tpool_t;


void tpool_go(tpool_t *thread_data);
void tpool_add_work_fcn(tpool_t *tdata, void *work_fcn);
tpool_t * tpool_setup(int num_threads, void *start_fcn, void *stop_fcn,
    void *sync_fcn, void *dispatch_fcn, void *udata);


#ifdef __cplusplus
}
#endif


#endif // #ifndef _THREADPOOL_H
