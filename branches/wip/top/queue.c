#include "queue.h"

Queue_t * newQueue(uint32_t sz)
{
	Queue_t *Q = (Queue_t *)malloc(sizeof(Queue_t));
	Q->Q = (uint32_t *)malloc(sz * sizeof(uint32_t));
	Q->sz = sz;
	Q->head = 0;
	Q->tail = 0;
	Q->len = 0;
	return Q;
}

void enqueue(Queue_t *Q, uint32_t e)
{
	Q->Q[Q->tail++] = e;
	Q->len++;

	if (Q->tail == Q->sz)
	{
		Q->tail = 0;
	}

	if (Q->len >= Q->sz)
	{
		printf("warning: Q overflowed\n");
		exit(1);
	}
	return;
}

uint32_t dequeue(Queue_t *Q)
{
	uint32_t e = -1;

	if (Q->len > 0)
	{
		e = Q->Q[Q->head];
		Q->head++;
		if (Q->head == Q->sz)
		{
			Q->head = 0;
		}
		Q->len--;
	}
	else
	{
		printf("warning: attempted to dequeue from an empty queue\n");
	}

	return e;
}

uint32_t peekqueue(Queue_t *Q)
{
	uint32_t e = -1;
	if (Q->len > 0)
	{
		e = Q->Q[Q->head];
	}
	return e;
}

void clearQueue(Queue_t *Q)
{
	free(Q->Q);
	Q->len = 0;
	Q->sz = 0;
	Q->head = 0;
	Q->tail = 0;
}
