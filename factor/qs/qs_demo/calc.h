/*
MIT License

Copyright (c) 2021 Ben Buhrow

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
 */

/*
    Implements an extendable arbitrary precision calculator.
    Supports function calls with optional arguments.
    Supports interactive use by keeping internal state (user variables).
    Depends on GMP for arbitrary arithmetic.
    Library Usage:
        calc_init();                                        // set up some internal storage
        process_expression(expression_string, metadata);    // process expression_string
        calc_finalize();                                    // free internal storage

*/

#ifndef YCALC_H
#define YCALC_H

#ifdef __cplusplus
extern "C" {
#endif


#include "gmp.h"
#include "factor.h"
#include "soe.h"

// maximum length of strings, or characters in a string to process
// as a chunk
#define GSTR_MAXSIZE 1024

// string stuff
typedef struct
{
    char* s;		// pointer to beginning of s
    int nchars;		// number of valid characters in s (including \0)
    int alloc;		// bytes allocated to s
} str_t;

// user variables
typedef struct
{
    char name[40];
    mpz_t data;
} uvar_t;

typedef struct
{
    uvar_t* vars;
    int num;
    int alloc;
} uvars_t;

// string variables, used in for/if commands, for instance.
typedef struct
{
    char name[40];
    char* data;
} strvar_t;

typedef struct
{
    strvar_t* vars;
    int num;
    int alloc;
} strvars_t;

// some basic functions for the string type
void sInit(str_t* s);
void sFree(str_t* s);
void sClear(str_t* s);
void sCopy(str_t* dest, str_t* src);
void sGrow(str_t* s, int size);
void toStr(char* src, str_t* dest);
void sAppend(const char* src, str_t* dest);

// the expression processor needs a simple stack for strings
typedef struct
{
    str_t** elements;	//an array of pointers to elements
    int num;			//number of elements
    int size;			//allocated number of elements in stack
    int top;			//the top element
    int type;			//is this a stack (0) or a queue (1)?
} bstack_t;

int stack_init(int num, bstack_t * stack, int type);
int stack_free(bstack_t* stack);
void push(str_t* str, bstack_t* stack);
int pop(str_t* str, bstack_t* stack);

// User defined functions that process other data can 
// use this structure for getting that data in and out 
// of the calculator.
typedef struct
{
    fact_obj_t* fobj;
    soe_staticdata_t* sdata;
} meta_t;

// arbitrary precision calculator interface
extern char* process_expression(char* input_exp, meta_t* metadata,
    int force_quiet, int no_convert_result);
extern void reset_preprocessor(void);
extern void calc_finalize();
extern int calc_init(uint64_t rand_seed);


#ifdef __cplusplus
}
#endif

#endif /* #ifndef YCALC_H */
