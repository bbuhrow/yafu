/**
 * based on threadpool.c (BSD licensed)
 * modified for use within Msieve by Jason Papadopoulos
 *
 * $Id: thread.c 892 2013-06-20 00:33:24Z jasonp_sf $
 *
 *  http://sourceforge.net/projects/cthreadpool
 *
 *  Created on: Dec 11, 2010
 *      Author: Tomer Heber (heber.tomer@gmail.com).
 */

#include <thread.h>

#include <pthread.h>

#define THREAD_POOL_DEBUG

#ifdef THREAD_POOL_DEBUG
#define REPORT_ERROR(...) fprintf (stderr,"line %d - ",__LINE__); fprintf (stderr, __VA_ARGS__); fprintf (stderr,"\n")
#else
#define REPORT_ERROR(...)
#endif /* THREAD_POOL_DEBUG */

struct threadpool_queue
{
	unsigned int head;
	unsigned int tail;
	unsigned int num_tasks;
	unsigned int max_tasks;
	void **tasks;
};

struct thread_init 
{
	int thread_num;
	struct threadpool *pool;
	thread_control_t control;
};

struct threadpool
{
	struct threadpool_queue tasks_queue;
	struct threadpool_queue free_tasks_queue;

	task_control_t *tasks;

	struct thread_init *thr_init;
	pthread_t *thr_arr;

	unsigned short num_of_threads;
	volatile unsigned short stop_flag;

	pthread_mutex_t free_tasks_mutex;
	pthread_cond_t free_tasks_cond;
	pthread_cond_t tasks_done_cond;

	pthread_mutex_t mutex;
	pthread_cond_t new_tasks_cond;
};

static void threadpool_queue_init(struct threadpool_queue *queue,
				int max_tasks)
{
	queue->head = 0;
	queue->tail = 0;
	queue->num_tasks = 0;
	queue->max_tasks = max_tasks;
	queue->tasks = (void **)xcalloc(max_tasks, sizeof(void *));
}

static void threadpool_queue_free(struct threadpool_queue *queue)
{
	free(queue->tasks);
}

static int threadpool_queue_enqueue(struct threadpool_queue *queue, void *data)
{
	if (queue->num_tasks == queue->max_tasks) {
		REPORT_ERROR("The queue is full, unable to add data to it.");
		return -1;
	}

	if (queue->tasks[queue->tail] != NULL) {
		REPORT_ERROR("A problem was detected in the queue (expected NULL, but found a different value).");
		return -1;
	}

	queue->tasks[queue->tail] = data;

	queue->num_tasks++;
	queue->tail++;

	if (queue->tail == queue->max_tasks) {
		queue->tail = 0;
	}

	return 0;
}

/**
 * This function removes and returns the head data element in the queue.
 *
 * @param queue The queue structure.
 * @return On success a data element is returned, on failure NULL is returned.
 */
static void *threadpool_queue_dequeue(struct threadpool_queue *queue)
{
	void *data;

	if (queue->num_tasks == 0) {
			REPORT_ERROR("Tried to dequeue from an empty queue.");
			return NULL;
	}

	data = queue->tasks[queue->head];

	queue->tasks[queue->head] = NULL;
	queue->num_tasks--;

	if (queue->num_tasks == 0) {
		queue->head = 0;
		queue->tail = 0;
	}
	else {
		queue->head++;
		if (queue->head == queue->max_tasks) {
			queue->head = 0;
		}
	}

	return data;
}

/**
 * This function checks if a given queue is empty.
 *
 * @param queue The queue structure.
 * @return 1 if the queue is empty, else 0.
 */
static int threadpool_queue_is_empty(struct threadpool_queue *queue)
{
	if (queue->num_tasks == 0) {
		return 1;
	}

	return 0;
}

/**
 * This function checks if a given queue is full.
 *
 * @param queue The queue structure.
 * @return 1 if the queue is full, else 0.
 */
static int threadpool_queue_is_full(struct threadpool_queue *queue)
{
	if (queue->num_tasks == queue->max_tasks) {
		return 1;
	}

	return 0;
}

/**
 * This function queries for the size of the given queue argument.
 *
 * @param queue The queue structure.
 * @return The size of the queue.
 */
static int threadpool_queue_getsize(struct threadpool_queue *queue)
{
	return queue->num_tasks;
}

static void threadpool_task_clear(task_control_t *task)
{
	memset(task, 0, sizeof(task_control_t));
}

/**
 * This function obtains a queued task from the pool and returns it.
 * If no such task is available the operation blocks.
 *
 * @param pool The thread pool structure.
 * @return A task or NULL on error (or if thread pool should shut down).
 */
static task_control_t * threadpool_task_get_task(struct threadpool *pool)
{
	task_control_t * task;

	if (pool->stop_flag) {
		/* The pool should shut down return NULL. */
		return NULL;
	}

	/* Obtain a task */
	if (pthread_mutex_lock(&(pool->mutex))) {
		perror("pthread_mutex_lock: ");
		return NULL;
	}

	while (threadpool_queue_is_empty(&(pool->tasks_queue)) && !pool->stop_flag) {
		/* Block until a new task arrives. */
		if (pthread_cond_wait(&(pool->new_tasks_cond),&(pool->mutex))) {
			perror("pthread_cond_wait: ");
			if (pthread_mutex_unlock(&(pool->mutex))) {
				perror("pthread_mutex_unlock: ");
			}

			return NULL;
		}
	}

	if (pool->stop_flag) {
		/* The pool should shut down return NULL. */
		if (pthread_mutex_unlock(&(pool->mutex))) {
			perror("pthread_mutex_unlock: ");
		}
		return NULL;
	}

	if ((task = (task_control_t *)threadpool_queue_dequeue(&(pool->tasks_queue))) == NULL) {
		/* Since task is NULL returning task will return NULL as required. */
		REPORT_ERROR("Failed to obtain a task from the jobs queue.");
	}

	if (pthread_mutex_unlock(&(pool->mutex))) {
		perror("pthread_mutex_unlock: ");
		return NULL;
	}

	return task;
}

/**
 * This is the routine the worker threads do during their life.
 *
 * @param data Contains a pointer to the startup data
 * @return NULL.
 */
#if defined(__GNUC__) && (__GNUC__ > 4 || __GNUC__ == 4 && __GNUC_MINOR__>1)

/* gcc on win32 needs to force 16-byte stack alignment on 
   thread entry, as this exceeds what windows may provide; see

   http://sourceware.org/ml/pthreads-win32/2008/msg00053.html
*/
__attribute__((force_align_arg_pointer))
#endif
static void *worker_thr_routine(void *data)
{
	struct thread_init *init = (struct thread_init *)data;
	int my_id = init->thread_num;
	struct threadpool *pool = init->pool;
	thread_control_t *t = &init->control;
	task_control_t *task;

	/* initialize thread-local state */

	if (t->init != NULL) {
		t->init(t->data, my_id);
	}

	while (1) {
		task = threadpool_task_get_task(pool);
		if (task == NULL) {
			if (pool->stop_flag) {
				/* Worker thr needs to exit (thread pool was shutdown). */
				break;
			}
			else {
				/* An error has occurred. */
				REPORT_ERROR("Warning an error has occurred when trying to obtain a worker task.");
				REPORT_ERROR("The worker thread has exited.");
				break;
			}
		}

		/* Execute task */

		if (task->init != NULL)
			task->init(task->data, my_id);

		if (task->run != NULL)
			task->run(task->data, my_id);

		if (task->shutdown != NULL)
			task->shutdown(task->data, my_id);

		/* Release the task by returning it to the free_task_queue. */
		threadpool_task_clear(task);
		if (pthread_mutex_lock(&(pool->free_tasks_mutex))) {
			perror("pthread_mutex_lock: ");
			REPORT_ERROR("The worker thread has exited.");
			break;
		}

		if (threadpool_queue_enqueue(&(pool->free_tasks_queue),task)) {
			REPORT_ERROR("Failed to enqueue a task to free tasks queue.");
			if (pthread_mutex_unlock(&(pool->free_tasks_mutex))) {
				perror("pthread_mutex_unlock: ");
			}

			REPORT_ERROR("The worker thread has exited.");
			break;
		}

		if (threadpool_queue_getsize(&(pool->free_tasks_queue)) == 1) {
			/* Notify all waiting threads that new tasks can added. */
			if (pthread_cond_broadcast(&(pool->free_tasks_cond))) {
				perror("pthread_cond_broadcast: ");
				if (pthread_mutex_unlock(&(pool->free_tasks_mutex))) {
					perror("pthread_mutex_unlock: ");
				}

				break;
			}
		}

		if (threadpool_queue_is_full(&(pool->free_tasks_queue)) == 1) {
			/* Notify any waiting threads that threadpool is not busy */

			if (pthread_cond_broadcast(&(pool->tasks_done_cond))) {
				perror("pthread_cond_broadcast: ");
				if (pthread_mutex_unlock(&(pool->free_tasks_mutex))) {
					perror("pthread_mutex_unlock: ");
				}

				break;
			}
		}

		if (pthread_mutex_unlock(&(pool->free_tasks_mutex))) {
			perror("pthread_mutex_unlock: ");
			REPORT_ERROR("The worker thread has exited.");
			break;
		}
	}

	/* tear down thread-local state */

	if (t->shutdown != NULL) {
		t->shutdown(t->data, my_id);
	}

	return NULL;
}

/**
 * This function does the following steps:
 * 1. It raises a flag that notifies the worker threads to stop working.
 * 2. It waits until all worker threads are done with their execution.
 * 3. It frees all the allocated memory of the threadpool struct.
 *
 * @param ptr The pool to stop its worker threads.

 * @return 0.
 */
void threadpool_free(struct threadpool *pool)
{
	int i;

	pool->stop_flag = 1;

	/* Wakeup all worker threads (broadcast operation). */
	if (pthread_mutex_lock(&(pool->mutex))) {
		perror("pthread_mutex_lock: ");
		REPORT_ERROR("Warning: Memory was not released.");
		REPORT_ERROR("Warning: Some of the worker threads may have failed to exit.");
		return;
	}

	if (pthread_cond_broadcast(&(pool->new_tasks_cond))) {
		perror("pthread_cond_broadcast: ");
		REPORT_ERROR("Warning: Memory was not released.");
		REPORT_ERROR("Warning: Some of the worker threads may have failed to exit.");
		return;
	}

	if (pthread_mutex_unlock(&(pool->mutex))) {
		perror("pthread_mutex_unlock: ");
		REPORT_ERROR("Warning: Memory was not released.");
		REPORT_ERROR("Warning: Some of the worker threads may have failed to exit.");
		return;
	}

	/* Wait until all worker threads are done. */
	for (i = 0; i < pool->num_of_threads; i++) {
		if (pthread_join(pool->thr_arr[i],NULL)) {
			perror("pthread_join: ");
		}
	}

	/* shut down any tasks that are still waiting */
	while (threadpool_queue_getsize(&(pool->tasks_queue))) {

		task_control_t *task = (task_control_t *)
				threadpool_queue_dequeue(&(pool->tasks_queue));

		if (task != NULL && task->shutdown != NULL) {
			task->shutdown(task->data, 0);
		}
	}

	/* Free all allocated memory. */
	threadpool_queue_free(&(pool->tasks_queue));
	threadpool_queue_free(&(pool->free_tasks_queue));
	free(pool->tasks);
	free(pool->thr_arr);
	free(pool->thr_init);
	free(pool);
}

struct threadpool* threadpool_init(int num_of_threads,
				int queue_size,
				thread_control_t *t)
{
	int i;
	struct threadpool *pool = (struct threadpool *)xcalloc(1,
					sizeof(struct threadpool));

	/* Init the mutex and cond vars. */
	if (pthread_mutex_init(&(pool->free_tasks_mutex),NULL)) {
		perror("pthread_mutex_init: ");
		free(pool);
		return NULL;
	}
	if (pthread_mutex_init(&(pool->mutex),NULL)) {
		perror("pthread_mutex_init: ");
		free(pool);
		return NULL;
	}
	if (pthread_cond_init(&(pool->free_tasks_cond),NULL)) {
		perror("pthread_mutex_init: ");
		free(pool);
		return NULL;
	}
	if (pthread_cond_init(&(pool->tasks_done_cond),NULL)) {
		perror("pthread_mutex_init: ");
		free(pool);
		return NULL;
	}
	if (pthread_cond_init(&(pool->new_tasks_cond),NULL)) {
		perror("pthread_mutex_init: ");
		free(pool);
		return NULL;
	}

	/* Init the queues. */
	threadpool_queue_init(&(pool->tasks_queue), queue_size);
	threadpool_queue_init(&(pool->free_tasks_queue), queue_size);
	pool->tasks = (task_control_t *)xmalloc(queue_size *
					sizeof(task_control_t));

	/* Add all the free tasks to the free tasks queue. */
	for (i = 0; i < queue_size; i++) {
		threadpool_task_clear((pool->tasks) + i);
		if (threadpool_queue_enqueue(&(pool->free_tasks_queue),(pool->tasks) + i)) {
			REPORT_ERROR("Failed to a task to the free tasks queue during initialization.");
			return NULL;
		}
	}

	/* Create the thr_arr. */
	if ((pool->thr_arr = malloc(sizeof(pthread_t) * num_of_threads)) == NULL) {
		perror("malloc: ");
		free(pool);
		return NULL;
	}

	if ((pool->thr_init = malloc(sizeof(struct thread_init) * num_of_threads)) == NULL) {
		perror("malloc: ");
		free(pool);
		return NULL;
	}

	/* Start the worker threads. */
	for (i = 0; i < num_of_threads; i++) {

		pool->thr_init[i].thread_num = i;
		pool->thr_init[i].pool = pool;
		pool->thr_init[i].control = *t;

		pool->num_of_threads = i + 1;
		if (pthread_create(&(pool->thr_arr[i]),NULL,
				worker_thr_routine,
				&(pool->thr_init[i]))) {
			perror("pthread_create:");

			threadpool_free(pool);

			return NULL;
		}
	}

	return pool;
}

int threadpool_add_task(struct threadpool *pool, task_control_t *new_task, int blocking)
{
	task_control_t *task;

	if (pool == NULL) {
		REPORT_ERROR("The threadpool received as argument is NULL.");
		return -1;
	}

	if (pthread_mutex_lock(&(pool->free_tasks_mutex))) {
		perror("pthread_mutex_lock: ");
		return -1;
	}

	/* Check if the free task queue is empty. */
	while (threadpool_queue_is_empty(&(pool->free_tasks_queue))) {
		if (!blocking) {
			/* Return immediately if the command is non blocking. */
			if (pthread_mutex_unlock(&(pool->free_tasks_mutex))) {
				perror("pthread_mutex_unlock: ");
				return -1;
			}

			return -2;
		}

		/* blocking is set to 1, wait until free_tasks queue has a task to obtain. */
		if (pthread_cond_wait(&(pool->free_tasks_cond),&(pool->free_tasks_mutex))) {
			perror("pthread_cond_wait: ");
			if (pthread_mutex_unlock(&(pool->free_tasks_mutex))) {
				perror("pthread_mutex_unlock: ");
			}

			return -1;
		}
	}

	/* Obtain an empty task. */
	if ((task = (task_control_t *)threadpool_queue_dequeue(&(pool->free_tasks_queue))) == NULL) {
		REPORT_ERROR("Failed to obtain an empty task from the free tasks queue.");
		if (pthread_mutex_unlock(&(pool->free_tasks_mutex))) {
			perror("pthread_mutex_unlock: ");
		}

		return -1;
	}

	if (pthread_mutex_unlock(&(pool->free_tasks_mutex))) {
		perror("pthread_mutex_unlock: ");
		return -1;
	}

	*task = *new_task;

	/* Add the task, to the tasks queue. */
	if (pthread_mutex_lock(&(pool->mutex))) {
		perror("pthread_mutex_lock: ");
		return -1;
	}

	if (threadpool_queue_enqueue(&(pool->tasks_queue),task)) {
		REPORT_ERROR("Failed to add a new task to the tasks queue.");
		if (pthread_mutex_unlock(&(pool->mutex))) {
			perror("pthread_mutex_unlock: ");
		}
		return -1;
	}

	if (threadpool_queue_getsize(&(pool->tasks_queue)) == 1) {
		/* Notify all worker threads that there are new jobs. */
		if (pthread_cond_broadcast(&(pool->new_tasks_cond))) {
			perror("pthread_cond_broadcast: ");
			if (pthread_mutex_unlock(&(pool->mutex))) {
				perror("pthread_mutex_unlock: ");
			}

			return -1;
		}
	}

	if (pthread_mutex_unlock(&(pool->mutex))) {
		perror("pthread_mutex_unlock: ");
		return -1;
	}

	return 0;
}

int threadpool_drain(struct threadpool *pool, int blocking)
{
	if (pthread_mutex_lock(&(pool->free_tasks_mutex))) {
		perror("pthread_mutex_lock: ");
		return -1;
	}

	/* Check if the free task queue is full. */
	while (!threadpool_queue_is_full(&(pool->free_tasks_queue))) {
		if (!blocking) {
			/* Return immediately if the command is non blocking. */
			if (pthread_mutex_unlock(&(pool->free_tasks_mutex))) {
				perror("pthread_mutex_unlock: ");
				return -1;
			}

			return 1;
		}

		/* blocking is set to 1, wait until free_tasks queue is full */
		if (pthread_cond_wait(&(pool->tasks_done_cond),&(pool->free_tasks_mutex))) {
			perror("pthread_cond_wait: ");
			if (pthread_mutex_unlock(&(pool->free_tasks_mutex))) {
				perror("pthread_mutex_unlock: ");
			}

			return -1;
		}
	}

	if (pthread_mutex_unlock(&(pool->free_tasks_mutex))) {
		perror("pthread_mutex_unlock: ");
		return -1;
	}

	return 0;
}

