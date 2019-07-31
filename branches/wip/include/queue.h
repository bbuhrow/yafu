
#ifndef QUEUE_H
#define QUEUE_H


#include <stdint.h>

typedef struct
{
	uint32_t *Q;
	uint32_t len;
	uint32_t sz;
	uint32_t head;
	uint32_t tail;
} Queue_t;

void clearQueue(Queue_t *Q);
uint32_t peekqueue(Queue_t *Q);
uint32_t dequeue(Queue_t *Q);
void enqueue(Queue_t *Q, uint32_t e);
Queue_t * newQueue(uint32_t sz);

#endif
