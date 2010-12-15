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
       				   --bbuhrow@gmail.com 11/24/09
----------------------------------------------------------------------*/

#include "yafu.h"
#include "yafu_stack.h"

/*
	implements a custom stack/queue data structure
	the type parameter specifies which it is - a stack or a queue
	the only difference is that the queue pops from node 0 and
	then adjusts all the pointers down 1, while the stack
	pops from the top and no adjustment is necessary.
	  
	the type of element in either structure is a str_t, which
	implements a variable length string
*/

int stack_init(int num, bstack_t *stack, int type)
{
	int i;
	stack->elements = (str_t **)malloc(num*sizeof(str_t*));		//array of elements
	//space for each element (a str_t)
	for (i=0;i<num;i++)
	{
		stack->elements[i] = (str_t *)malloc(sizeof(str_t));
		//init each element (char array);
		sInit(stack->elements[i]);
	}
	stack->size = num;				//number of allocated stack elements
	stack->num = 0;					//number of currently occupied stack elements
	stack->top = 0;
	stack->type = type;

	return 0;
}

int stack_free(bstack_t *stack)
{
	int i;

	for (i=0;i<stack->size;i++)
	{
		sFree(stack->elements[i]);	//first free any occupied stack elements
		free(stack->elements[i]);
	}
	free(stack->elements);			//then free the stack

	return 0;
}

void push(str_t *str, bstack_t *stack)
{
	//str_t *newstr;

	//add an element to the stack, growing the stack if necessary
	if (stack->num >= stack->size)
	{
		stack->size *=2;
		stack->elements = (str_t **)realloc(stack->elements,
			stack->size * sizeof(str_t*));
		if (stack->elements == NULL)
		{
			printf("error allocating stack space\n");
			return;
		}
	}

	//create a new string
	//newstr = (str_t *)malloc(sizeof(str_t));
	//sInit(newstr);
	//sCopy(newstr,str);
	sCopy(stack->elements[stack->num],str);
	stack->num++;

	//both stacks and queues push to the same side of the array
	//the top element and the number of elements are the same
	stack->top = stack->num - 1;
	//store the pointer to it in the stack
	//stack->elements[stack->top] = newstr;

	return;
}

int pop(str_t *str, bstack_t *stack)
{
	//take an element off the stack.  return 0 if there are no elements
	//pass in a pointer to a string.  if necessary, this routine will 
	//reallocate space for the string to accomodate its size.  If this happens
	//the pointer to the string's (likely) new location is automatically
	//updated and returned.
	int i;

	//copy out the string at the top of the stack
	//then free the stack's copy.
	if (stack->num != 0)
	{
		stack->num--;
		if (stack->type == QUEUE)
		{
			//for queues, the top element is always node 0
			sCopy(str,stack->elements[0]);
			sFree(stack->elements[0]);
			free(stack->elements[0]);
			stack->top--;
			//now we need to adjust all the pointers down 1
			for (i=1;i<stack->num;i++)
				stack->elements[i-1] = stack->elements[i];
		}
		else
		{
			sCopy(str,stack->elements[stack->top]);
			//sFree(stack->elements[stack->top]);
			//free(stack->elements[stack->top]);
			stack->top--;
		}
		return 1;
	}
	else
		return 0;
}

