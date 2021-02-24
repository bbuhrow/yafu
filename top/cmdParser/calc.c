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


#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <ctype.h>  // isdigit, isspace
#include <math.h>   // log, sqrt
#include "calc.h"
#include "gmp.h"
#include "mpz_aprcl.h"
#include "factor.h"
#include "autofactor.h"
#include "qs.h"
#include "yafu_ecm.h"
#include "nfs.h"

// define this for debug or a verbose interface
#define CALC_VERBOSE 0

// the number of functions defined
#define NUM_FUNC 76

// symbols in calc
#define EOE 1
#define IMM 2
#define NUM 3
#define OP 4
#define RP 5
#define LP 6
#define CH 7
#define AMBIG 8
#define COMMA 9
#define SPACE 10
#define POSITIVE 0
#define NEGATIVE 1

// operator associativity
#define RIGHT 1
#define LEFT 0

// bases
#define DEC 10
#define HEX 16
#define OCT 8
#define BIN 2

// useful definitions
#define MIN(a,b) ((a) < (b)? (a) : (b))
#define MAX(a,b) ((a) > (b)? (a) : (b))

// local function declarations
void handle_singleop(char* arg1, int op);
int single_op(char s);
int dual_op(char s);
str_t* preprocess(str_t* str, int* num);
int get_el_type(char s);
int processIMM(int opcode, str_t* str);
int calc(str_t* str, meta_t* metadata);
int isInt(char s);
int getIMM(char s);
int op_precedence(char* s1, char* s2, int assoc);
int getAssoc(char* s);
int processOP(char* s, str_t* n1, str_t* n2);
int getOP(char s);
int isEOE(char s);
int getFunc(char* s, int* nargs);
int feval(int func, int nargs, meta_t* metadata);
int new_uvar(const char* name, mpz_t data);
int set_uvar(const char* name, mpz_t data);
int get_uvar(const char* name, mpz_t data);
void free_uvars();
int new_strvar(const char* name, char* data);
int set_strvar(const char* name, char* data);
int get_strvar(const char* name, char* data);
char* get_strvarname(const char* data);
int is_strvar(const char* name);
void free_strvars();
int invalid_dest(char* dest);
int invalid_num(char* num);
char** tokenize(char* in, int* token_types, int* num_tokens);
int get_el_type2(char s);
int is_new_token(int el_type, int el_type2);
int invalid_dest(char* dest);
void calc_with_assignment(str_t* in, meta_t* metadata, int force_quiet);

#ifndef _MSC_VER
#define strtok_s strtok_r
#endif

// local data
char opchar[9] = { '=', '<', '>', '+', '-', '*', '/', '%', '^' }; // , '='};
char imms[3] = {'!','#','-'};
const int numopchars = 9;
char choperands[5][GSTR_MAXSIZE];
mpz_t operands[5];
int for_cnt = 0;
int forp_cnt = 0;
int forf_cnt = 0;
int if_cnt = 0;

// calculator data
int IBASE;
int OBASE;
str_t gstr1, gstr2, gstr3;
gmp_randstate_t gmp_randstate;
uvars_t uvars;
strvars_t strvars;

static char function_names[NUM_FUNC][11] = {
    "fib", "luc", "expr", "gcd", "jacobi",
    "rand", "lg2", "log", "ln", "size",
    "issquare", "isprime", "sqrt", "modinv", "modexp",
    "nroot", "shift", "ispow", "randb", "+",
    "-", "*", "/", "!", "#",
    "eq", "<<", ">>", "%", "^",
    "redc", "bitxor", "bitand", "bitor", "onecomp",
    "lte", "gte", "<", ">", "popcnt",
    "nextprime", "print", "lcm", "exit", "abs",
    "extgcd", "fac2", "facm", "binom", "randp",
    "hamdist", "snfs", "rsa", "factor", "pm1",
    "pp1", "rho", "trial", "shanks", "siqs",
    "primes", "torture", "ecm", "llt", "siqsbench",
    "sigma", "totient", "smallmpqs", "testrange", "sieverange",
    "fermat", "nfs", "tune", "bpsw", "aprcl",
    "semiprimes"};

static int function_nargs[NUM_FUNC] = {
    1, 1, 1, 2, 2, 
    1, 1, 1, 1, 1, 
    1, 1, 1, 2, 3,
    2, 2, 1, 1, 2,
    2, 2, 2, 1, 1, 
    2, 2, 2, 2, 2, 
    3, 2, 2, 2, 1,
    2, 2, 2, 2, 1, 
    1, 1, 2, 0, 1,
    2, 1, 2, 2, 1,
    2, 2, 1, 1, 1, 
    3, 1, 2, 1, 1, 
    3, 2, 2, 1, 0, 
    1, 1, 1, 4, 4, 
    3, 1, 0, 1, 1,
    2};


// =====================================================================
// basic string functions
// =====================================================================

void sInit(str_t* s)
{
    s->s = (char*)malloc(GSTR_MAXSIZE * sizeof(char));
    if (s->s == NULL)
    {
        printf("couldn't allocate str_t in sInit\n");
        exit(-1);
    }
    s->s[0] = '\0';
    s->nchars = 1;
    s->alloc = GSTR_MAXSIZE;
    return;
}

void sFree(str_t* s)
{
    free(s->s);
    return;
}

void sClear(str_t* s)
{
    s->s[0] = '\0';
    s->nchars = 1;
    return;
}

void toStr(char* src, str_t* dest)
{
    if ((int)strlen(src) > dest->alloc)
    {
        sGrow(dest, strlen(src) + 10);
        dest->alloc = strlen(src) + 10;
    }
    memcpy(dest->s, src, strlen(src) * sizeof(char));
    dest->s[strlen(src)] = '\0';
    dest->nchars = strlen(src) + 1;

    return;
}

void sGrow(str_t* str, int size)
{
    //printf("growing str_t size...\n");
    str->s = (char*)realloc(str->s, size * sizeof(char));
    if (str->s == NULL)
    {
        printf("unable to reallocate string in sGrow\n");
        exit(-1);
    }
    str->alloc = size;

    return;
}

void sAppend(const char* src, str_t* dest)
{
    if (((int)strlen(src) + dest->nchars) >= dest->alloc)
    {
        sGrow(dest, strlen(src) + dest->nchars + 10);
        dest->alloc = strlen(src) + dest->nchars + 10;
    }

    memcpy(dest->s + dest->nchars - 1, src, strlen(src) * sizeof(char));
    dest->nchars += strlen(src);	//already has a null char accounted for
    dest->s[dest->nchars - 1] = '\0';

    return;
}

void sCopy(str_t* dest, str_t* src)
{
    if (dest->alloc < src->nchars + 2)
    {
        dest->s = (char*)realloc(dest->s, (src->nchars + 2) * sizeof(char));
        dest->alloc = src->nchars + 2;
    }
    memcpy(dest->s, src->s, src->nchars * sizeof(char));
    dest->nchars = src->nchars;
    return;
}

// =====================================================================
// basic stack functions
// =====================================================================

int stack_init(int num, bstack_t* stack, int type)
{
    int i;
    stack->elements = (str_t * *)malloc(num * sizeof(str_t*));		//array of elements
    //space for each element (a str_t)
    for (i = 0; i < num; i++)
    {
        stack->elements[i] = (str_t*)malloc(sizeof(str_t));
        //init each element (char array);
        sInit(stack->elements[i]);
    }
    stack->size = num;				//number of allocated stack elements
    stack->num = 0;					//number of currently occupied stack elements
    stack->top = 0;
    stack->type = type;

    return 0;
}

int stack_free(bstack_t* stack)
{
    int i;

    for (i = 0; i < stack->size; i++)
    {
        sFree(stack->elements[i]);	//first free any occupied stack elements
        free(stack->elements[i]);
    }
    free(stack->elements);			//then free the stack

    return 0;
}

void push(str_t* str, bstack_t* stack)
{
    // add an element to the stack, growing the stack if necessary
    if (stack->num >= stack->size)
    {
        stack->size *= 2;
        stack->elements = (str_t * *)realloc(stack->elements,
            stack->size * sizeof(str_t*));
        if (stack->elements == NULL)
        {
            printf("error allocating stack space\n");
            return;
        }
    }

    sCopy(stack->elements[stack->num], str);
    stack->num++;

    // both stacks and queues push to the same side of the array
    // the top element and the number of elements are the same
    stack->top = stack->num - 1;

    return;
}

int pop(str_t* str, bstack_t* stack)
{
    // take an element off the stack.  return 0 if there are no elements
    // pass in a pointer to a string.  if necessary, this routine will 
    // reallocate space for the string to accomodate its size.  If this happens
    // the pointer to the string's (likely) new location is automatically
    // updated and returned.
    int i;

    // copy out the string at the top of the stack
    // then free the stack's copy.
    if (stack->num != 0)
    {
        stack->num--;
        if (stack->type == 1)   // this is a queue
        {
            // for queues, the top element is always node 0
            sCopy(str, stack->elements[0]);
            sFree(stack->elements[0]);
            free(stack->elements[0]);
            stack->top--;
            // now we need to adjust all the pointers down 1
            for (i = 1; i < stack->num; i++)
            {
                stack->elements[i - 1] = stack->elements[i];
            }
        }
        else
        {
            sCopy(str, stack->elements[stack->top]);
            stack->top--;
        }
        return 1;
    }
    else
    {
        return 0;
    }
}


// =====================================================================
// calc
// =====================================================================

void reset_preprocessor(void) {
    for_cnt = 0;
    forp_cnt = 0;
    forf_cnt = 0;
    if_cnt = 0;
    return;
}

int calc_init(uint64_t rand_seed)
{
	int i;
	// user variables space
	uvars.vars = (uvar_t *)malloc(10 * sizeof(uvar_t));
	uvars.alloc = 10;
	for (i=0;i<uvars.alloc;i++)
		mpz_init(uvars.vars[i].data);
	strcpy(uvars.vars[0].name,"ans");
	uvars.num = 1;

    // string variable space
    strvars.vars = (strvar_t *)malloc(10 * sizeof(strvar_t));
    strvars.alloc = 10;
    for (i = 0; i<strvars.alloc; i++)
        strvars.vars[i].data = (char *)malloc(GSTR_MAXSIZE * sizeof(char));
    strvars.num = 0;

    // mpz operands to functions
    for (i = 0; i < 5; i++)
    {
        mpz_init(operands[i]);
    }

    //global strings, used mostly for logprint stuff
    sInit(&gstr1);
    sInit(&gstr2);
    sInit(&gstr3);

    //global i/o base
    IBASE = DEC;
    OBASE = DEC;

    gmp_randinit_default(gmp_randstate);
    gmp_randseed_ui(gmp_randstate, rand_seed);

    reset_preprocessor();
	return 1;
}

void calc_finalize()
{
    int i;
	free_uvars();
    free_strvars();
    for (i = 0; i < 5; i++)
    {
        mpz_clear(operands[i]);
    }

    sFree(&gstr1);
    sFree(&gstr2);
    sFree(&gstr3);

    return;
}

int get_el_type2(char s)
{
	// there are several types of characters in an expression.  
	// decide which type this is
	if (isdigit(s) || (s <= 90 && s >= 65))
		return NUM;
	else if (s == '(')
		return LP;
	else if (s == ')')
		return RP;
	else if (s == '-')
		return AMBIG;
	else if (getIMM(s) >= 0)
		return IMM;
	else if (getOP(s) >= 0)
		return OP;
	else if (isEOE(s))
		return EOE;
	else if (s == ',')
		return COMMA;
	else if ((s <= 122 && s >= 95) || s == 39)
		return CH;
	else if (isspace(s))
		return SPACE;
	else
		return -1;
}

int isEOE(char s)
{
	if (s == 0)
		return 1;
	else
		return 0;
}

int getIMM(char s)
{
	int i;
    for (i = 0; i < 3; i++)
    {
        if (s == imms[i])
            return i;
    }

	return -1;
}

int getOP(char s)
{
	//return >=0 if this char is a opchar
	int i;

    for (i = 0; i < numopchars; i++)
    {
        if (opchar[i] == s)
            return i;
    }

	return -1;
}

int getAssoc(char *s)
{
	if (strcmp(s,"^") == 0)
		return RIGHT;
	else
		return LEFT;
}

int getPrecedence(char *s)
{
    if (strcmp(s, "=") == 0)  return -1;
    if (strcmp(s, "<<") == 0) return 0;
    if (strcmp(s, ">>") == 0) return 0;
    if (strcmp(s, "+") == 0)  return 1;
    if (strcmp(s, "-") == 0)  return 1;
    if (strcmp(s, "*") == 0)  return 2;
    if (strcmp(s, "/") == 0)  return 2;
    if (strcmp(s, "%") == 0)  return 2;
    if (strcmp(s, "^") == 0)  return 3;
    if (strcmp(s, "\\") == 0) return 4;
    return 0;
}

int op_precedence(char *s1, char *s2, int assoc)
{
	//if associativity is RIGHT, then use strictly greater than
	//else use greater than or equal to
	int p1=0,p2=0;

    p1 = getPrecedence(s1);
    p2 = getPrecedence(s2);

	if (assoc == LEFT)
		return p1 >= p2;
	else 
		return p1 > p2;
}

int is_closed(char *line, char *stopptr)
{
    int i;
    int openp, closedp, openb, closedb;

    openp = openb = closedp = closedb = 0;
    for (i = 0; i < strlen(line) && (&line[i] != stopptr); i++)
    {
        if (line[i] == '(') openp++;
        if (line[i] == ')') closedp++;
        if (line[i] == '{') openb++;
        if (line[i] == '}') closedb++;
    }
    if ((openp == closedp) && (openb == closedb))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

int exp_is_closed(char *start, char *stop)
{
    int i;
    int openp, closedp, openb, closedb;

    openp = openb = closedp = closedb = 0;

    i = 0;
    while ((i < strlen(start)) && ((start + i) < stop))
    {
        if (start[i] == '(') openp++;
        if (start[i] == ')') closedp++;
        if (start[i] == '{') openb++;
        if (start[i] == '}') closedb++;
        i++;
    }
    if ((openp == closedp) && (openb == closedb))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

int find_offset_matching_brace(char *ptr, char type)
{
    // find the matching brace or paren.
    // assume input string pointer starts just
    // beyond the opening brace or paren.
    int openp = 0, closedp = 0;
    char openchar = '(';
    int i;

    if (type == ')')
        openchar = '(';
    if (type == '}')
        openchar = '{';

    for (i = 0; i < strlen(ptr); i++)
    {
        if (ptr[i] == openchar) openp++;
        if (ptr[i] == type) {
            closedp++;
            if (closedp > openp)
                break;
        }
    }
    return i;
}

str_t* preprocess(str_t* in, int* numlines)
{
    // preprocess the expression in 'in'.
    // the expression evaluator is good at evaluating single-line
    // expressions (not assignments) that can be converted into
    // a post-fix notation consisting of operands, numbers, or variables.
    // The expression evaluator is therefore not able to handle assigments,
    // multiple lines, or compound expressions.
    // This preprocessor converts such complex syntax into a notation that
    // the expression evaluator can deal with.
    // We return the number of lines found and an array of strings,
    // one line per string.  user must free 'out'.
    int i, j, k;
    char* ptr;
    char str[1024];
    int openp, closedp, openb, closedb;
    str_t* out;
    str_t* current;

    if (in->s[0] == '{')
    {
        // if the first character is a open brace '{', then we need to
        // reformat this block expression by removing one level of
        // braces and removing leading or trailing commas.  The driver
        // will have placed commas between separate lines within the braces.
        strcpy(in->s, in->s + 1);
        in->s[strlen(in->s) - 1] = '\0';
        i = 0;
        while (in->s[i] == ',')
            in->s[i++] = ' ';

        i = 1;
        while ((in->s[strlen(in->s) - i] == ',') ||
            (in->s[strlen(in->s) - i] == '}') ||
            (in->s[strlen(in->s) - i] == ')'))
        {
            if ((in->s[strlen(in->s) - i] == ',') ||
                (in->s[strlen(in->s) - i] == '}'))
                in->s[strlen(in->s) - i] = ' ';
            i++;
        }
    }

    j = 0;
    // remove white space from the input
    for (i = 0; i < in->nchars; i++)
    {
        if (isspace(in->s[i]))
            continue;

        in->s[j] = in->s[i];
        j++;
    }
    in->s[j] = '\0';

    // copy input to first output location
    *numlines = 1;
    out = (str_t*)malloc(sizeof(str_t));
    current = &out[0];
    sInit(current);
    sCopy(current, in);

    // algebraic simplification (this would be cool...)    

    // reformat 'for' and 'if' tokens as functions 
    // taking string arguments, where the arguments are the various text
    // components of the function.  The actual looping is handled by 
    // recursive calls to process_expression from within the function evaluator.
    if (((ptr = strstr(current->s, "for(")) != NULL) &&
        exp_is_closed(current->s, ptr))
    {
        // new for loop
        char pre[8], start[80], vname[20];
        sprintf(pre, "for%d", for_cnt++);

        // save the beginning part of the command, if any
        if (ptr != current->s)
        {
            int n = (int)(ptr - current->s);
            strncpy(start, current->s, MIN(n, 79));
            start[n] = '\0';
        }
        else
            strcpy(start, "");

        // tokenize the loop
        //ptr = strtok(&current->s[4], ";");
        ptr = strtok(ptr, ";");
        if (ptr == NULL)
        {
            printf("badly formatted for loop: for(init; test; iter; body)\n");
            exit(3);
        }
        sprintf(vname, "%s_init", pre);
        if (set_strvar(vname, ptr + 4))
            new_strvar(vname, ptr + 4);
        ptr = strtok(NULL, ";");
        if (ptr == NULL)
        {
            printf("badly formatted for loop: for(init; test; iter; body)\n");
            exit(3);
        }
        sprintf(vname, "%s_test", pre);
        if (set_strvar(vname, ptr))
            new_strvar(vname, ptr);
        ptr = strtok(NULL, ";");
        if (ptr == NULL)
        {
            printf("badly formatted for loop: for(init; test; iter; body)\n");
            exit(3);
        }
        sprintf(vname, "%s_iter", pre);
        if (set_strvar(vname, ptr))
            new_strvar(vname, ptr);

        // this can't just find any old ')', it has to find the matching one.
        ptr = ptr + strlen(ptr) + 1;
        openp = closedp = 0;
        for (i = 0; i < strlen(ptr); i++)
        {
            if (ptr[i] == '(') openp++;
            if (ptr[i] == ')') {
                closedp++;
                if (closedp > openp)
                    break;
            }
        }

        if (ptr[i] == '\0')
        {
            printf("badly formatted for loop: for(init; test; iter; body)\n");
            exit(3);
        }

        ptr[i] = '\0';
        sprintf(vname, "%s_body", pre);
        if (set_strvar(vname, ptr))
            new_strvar(vname, ptr);

        if (ptr[i + 1] != '\0')
        {
            sprintf(str, "%s", ptr + i + 1);
            sprintf(current->s, "%sfor(%s_init, %s_test, %s_iter, %s_body);%s",
                start, pre, pre, pre, pre, str);
        }
        else
        {
            sprintf(current->s, "%sfor(%s_init, %s_test, %s_iter, %s_body);",
                start, pre, pre, pre, pre);
        }
    }

    if (0 && ((ptr = strstr(current->s, "forprime(")) != NULL) &&
        exp_is_closed(current->s, ptr))
    {
        // new for loop
        char pre[8], start[80], vname[20];
        sprintf(pre, "forp%d", forp_cnt++);

        // save the beginning part of the command, if any
        if (ptr != current->s)
        {
            int n = (int)(ptr - current->s);
            strncpy(start, current->s, MIN(n, 79));
            start[n] = '\0';
        }
        else
            strcpy(start, "");

        // tokenize the loop
        ptr = strtok(ptr, ";");
        if (ptr == NULL)
        {
            printf("badly formatted forprime loop: forprime(var=start; stop; body)\n");
            exit(3);
        }
        sprintf(vname, "%s_start", pre);
        if (set_strvar(vname, ptr + 9))
            new_strvar(vname, ptr + 9);
        ptr = strtok(NULL, ";");
        if (ptr == NULL)
        {
            printf("badly formatted forprime loop: forprime(var=start; stop; body)\n");
            exit(3);
        }
        sprintf(vname, "%s_stop", pre);
        if (set_strvar(vname, ptr))
            new_strvar(vname, ptr);

        // this can't just find any old ')', it has to find the matching one.
        ptr = ptr + strlen(ptr) + 1;
        openp = closedp = 0;
        for (i = 0; i < strlen(ptr); i++)
        {
            if (ptr[i] == '(') openp++;
            if (ptr[i] == ')') {
                closedp++;
                if (closedp > openp)
                    break;
            }
        }

        if (ptr[i] == '\0')
        {
            printf("badly formatted for loop: forprime(var=start; stop; body)\n");
            exit(3);
        }

        ptr[i] = '\0';
        sprintf(vname, "%s_body", pre);
        if (set_strvar(vname, ptr))
            new_strvar(vname, ptr);

        if (ptr[i + 1] != '\0')
        {
            sprintf(str, "%s", ptr + i + 1);
            sprintf(current->s, "%sforprime(%s_start, %s_stop, %s_body);%s",
                start, pre, pre, pre, str);
        }
        else
            sprintf(current->s, "%sforprime(%s_start, %s_stop, %s_body);", start, pre, pre, pre);
    }

    if (0 && ((ptr = strstr(current->s, "forfactors(")) != NULL) &&
        exp_is_closed(current->s, ptr))
    {
        // new for loop
        char pre[8], start[80], vname[20];
        sprintf(pre, "forf%d", forf_cnt++);

        // save the beginning part of the command, if any
        if (ptr != current->s)
        {
            int n = (int)(ptr - current->s);
            strncpy(start, current->s, MIN(n, 79));
            start[n] = '\0';
        }
        else
            strcpy(start, "");

        // tokenize the loop
        ptr = strtok(ptr, ";");
        if (ptr == NULL)
        {
            printf("badly formatted forfactors loop: forfactors(init, body)\n");
            exit(3);
        }
        sprintf(vname, "%s_init", pre);
        if (set_strvar(vname, ptr + 11))
            new_strvar(vname, ptr + 11);

        // this can't just find any old ')', it has to find the matching one.
        ptr = ptr + strlen(ptr) + 1;
        openp = closedp = 0;
        for (i = 0; i < strlen(ptr); i++)
        {
            if (ptr[i] == '(') openp++;
            if (ptr[i] == ')') {
                closedp++;
                if (closedp > openp)
                    break;
            }
        }

        if (ptr[i] == '\0')
        {
            printf("badly formatted for loop: forfactors(init, body)\n");
            exit(3);
        }

        ptr[i] = '\0';
        sprintf(vname, "%s_body", pre);
        if (set_strvar(vname, ptr))
            new_strvar(vname, ptr);

        if (ptr[i + 1] != '\0')
        {
            sprintf(str, "%s", ptr + i + 1);
            sprintf(current->s, "%sforfactors(%s_init, %s_body);%s",
                start, pre, pre, str);
        }
        else
            sprintf(current->s, "%sforfactors(%s_init, %s_body);", start, pre, pre);
    }

    if (((ptr = strstr(current->s, "if(")) != NULL) &&
        exp_is_closed(current->s, ptr))
    {
        // new if statement
        char pre[8], start[80], vname[20];
        sprintf(pre, "if%d", if_cnt++);

        // save the beginning part of the command, if any
        if (ptr != current->s)
        {
            int n = (int)(ptr - current->s);
            strncpy(start, current->s, MIN(n, 79));
            start[n] = '\0';
        }
        else
            strcpy(start, "");

        // tokenize the branch
        char* eptr;
        ptr = strtok(&current->s[3], ";");
        if (ptr == NULL)
        {
            printf("badly formatted if statement: if(condition; true-body; [false-body])\n");
            exit(3);
        }
        sprintf(vname, "%s_cond", pre);
        if (set_strvar(vname, ptr))
            new_strvar(vname, ptr);
        eptr = strtok(NULL, ";");
        if (eptr == NULL)
        {
            printf("badly formatted if statement: if(condition; true-body; [false-body])\n");
            exit(3);
        }
        else if (eptr[strlen(eptr) + 1] == '\0')
        {
            // no else statement and no output suppression character
            strncpy(str, eptr, strlen(eptr) - 1);
            str[strlen(eptr) - 1] = '\0';
            sprintf(vname, "%s_body", pre);
            if (set_strvar(vname, str))
                new_strvar(vname, str);

            sprintf(current->s, "%sif(%s_cond, %s_body);", start, pre, pre);
        }
        else
        {
            // either an else statement or an output suppression character or both
            if (eptr[strlen(eptr) + 1] == ';')
            {
                // both
                sprintf(str, "%s;", eptr);
                sprintf(vname, "%s_body", pre);
                if (set_strvar(vname, str))
                    new_strvar(vname, str);

                ptr = strtok(NULL, "\0");
                strncpy(str, ptr, strlen(ptr) - 1);
                str[strlen(ptr) - 1] = '\0';
                sprintf(vname, "%s_elsebody", pre);
                if (set_strvar(vname, str))
                    new_strvar(vname, str);
                sprintf(current->s, "%sif(%s_cond, %s_body, %s_elsebody);",
                    start, pre, pre, pre);
            }
            else if (eptr[strlen(eptr) + 1] == ')')
            {
                // just the if, with an output suppression character
                sprintf(str, "%s;", eptr);
                sprintf(vname, "%s_body", pre);
                if (set_strvar(vname, str))
                    new_strvar(vname, str);

                sprintf(current->s, "%sif(%s_cond, %s_body);", start, pre, pre);
            }
            else
            {
                // an else with no output suppression character
                sprintf(str, "%s", eptr);
                sprintf(vname, "%s_body", pre);
                if (set_strvar(vname, str))
                    new_strvar(vname, str);

                ptr = strtok(NULL, "\0");
                strncpy(str, ptr, strlen(ptr) - 1);
                str[strlen(ptr) - 1] = '\0';
                sprintf(vname, "%s_elsebody", pre);
                if (set_strvar(vname, str))
                    new_strvar(vname, str);

                sprintf(current->s, "%sif(%s_cond, %s_body, %s_elsebody);",
                    start, pre, pre, pre);
            }
        }
    }

    // search for commas within 'closed' areas and separate them into
    // a sequence of individual expressions.  A 'closed' area is text
    // that is not inside any parenthesis or brace.
    // for example, "i=0,j=0" should be parsed as the two expressions
    // i=0
    // j=0
    // but "for(i=0,j=0;..." is ignored because the comma is inside
    // an open parenthesis (it will eventually show up here without
    // the surrounding 'for')
    // to do this we go through the input string one character at a 
    // time and count open/closed parens/braces, while looking for
    // commas.
    openp = openb = closedp = closedb = 0;
    k = strlen(current->s);
    for (i = 0, j = 0; i < k; i++, j++)
    {
        if (current->s[i] == '(') openp++;
        if (current->s[i] == ')') closedp++;
        if (current->s[i] == '{') openb++;
        if (current->s[i] == '}') closedb++;
        if (current->s[i] == ',') {
            if ((openp == closedp) && (openb == closedb))
            {
                (*numlines)++;
                out = (str_t*)realloc(out, *numlines * sizeof(str_t));
                out[*numlines - 2].s[i] = '\0';
                sInit(&out[*numlines - 1]);
                current = &out[*numlines - 1];
                toStr(&out[*numlines - 2].s[i + 1], &out[*numlines - 1]);
                k = current->nchars;
                i = 0;
            }
        }
    }
    //toStr(current->s, &out[*numlines - 1]);

    return out;
}

int is_new_token(int el_type, int el_type2)
{

	if (el_type == EOE || el_type == LP || el_type == RP)
		return 1;

	if (el_type != el_type2)
	{
		// types are different
		if (el_type == CH && el_type2 == NUM)
		{
			// but this could be a function or variable name
			// so not different
			return 0;
		}
		else if (el_type == NUM && el_type2 == CH)
		{
			// but this could be a function or variable name
			// so not different
			return 0;
		}
		else
			return 1;
	}
	return 0;
}

char** tokenize(char *in, int *token_types, int *num_tokens)
{
	// take a string as input
	// break it into tokens
	// create an array of strings for each token
	// return the pointer to the array and the number of elements
	// in the array.  this will all have to be freed later
	// by the caller

	// a token in this context is one of the following things:
	//   a number, possibly including a base prefix (0x, 0d, 0b, etc)
	//   a variable name
	//   a function name
	//   an operator string (includes parens, commas)

	// read the string one character at a time
	// for each character read, decide if we've found the start of a new token

	int inpos, i, el_type, el_type2, token_alloc, tmpsize = GSTR_MAXSIZE;
	int len = strlen(in);
	char ch;
	char *tmp;
	char **tokens;

	token_alloc = 100;		//100 tokens
	tokens = (char **)malloc(token_alloc * sizeof(char *));
	*num_tokens = 0;

	tmp = (char *)malloc(GSTR_MAXSIZE * sizeof(char));

	// get the first character and check the type
	inpos = 0;
	i=1;
	ch = in[inpos];
	tmp[i-1] = ch;
	el_type = get_el_type2(ch);

	// when an expression gets cast into postfix, it aquires a leading
	// space which we can skip here
	if (el_type == SPACE)
	{
		inpos = 1;
		i=1;
		ch = in[inpos];
		tmp[i-1] = ch;
		el_type = get_el_type2(ch);
	}

	// ambiguous types:
	// a "-" can be either a num (if a negative sign) or an operator
	// a number can be a number or a string (num or func/var name)
	// a letter can be a string or a number (hex or func/var name)
	// we can tell them apart from the surrounding context
	// 
	// negative signs never have a num type before (or anything that
	// can be evaluated as a num , i.e. ")"
	// 
	// if we are reading CH's don't stop interpreting them as CH's until
	// we find a non-CH or non-NUM (use a flag)
	// 
	// watch for magic combinations '0x' '0d', etc.  set a flag to 
	// interpret what follows as the appropriate kind of num, and
	// discard the '0x', etc.
	// once this is fixed, change the final stack evaluation to print hex
	// strings to save some conversion time.
	if (el_type == AMBIG)
	{
		if (get_el_type2(in[inpos+1]) == NUM)
			el_type = NUM;
		else
			el_type = OP;
	}
	while (inpos < len)
	{
		// get another character and check the type
		inpos++;
		// if el_type == EOE, then no reason to keep reading.  This bug didn't seem to cause
		// any crashes, but couldn't have been healthy...
		if (el_type == EOE)
			break;
		ch = in[inpos];
		el_type2 = get_el_type2(ch);
		if (el_type2 == AMBIG)
		{
			switch (get_el_type2(in[inpos-1]))
			{
			case OP:
				el_type2 = NUM;
				break;
			case LP:
				el_type2 = NUM;
				break;
			case RP:
				el_type2 = OP;
				break;
			case CH:
				el_type2 = OP;
				break;
			case IMM:
				el_type2 = OP;
				break;
			case NUM:
				el_type2 = OP;
				break;
			case COMMA:
				el_type2 = NUM;
				break;
			case SPACE:
				// when processing postfix strings, we need this
				el_type2 = OP;
				break;
			default:
				printf("misplaced - sign\n");
				for (i=0;i< *num_tokens; i++)
					free(tokens[i]);
				free(tokens);
				free(tmp);
				return NULL;
			}
		}

		if (is_new_token(el_type,el_type2) || el_type == EOE)
		{
			if (el_type == EOE)
				break;

			if (el_type == -1)
			{
				// unrecognized character.  clear all tokens and return;
				printf("unrecognized character in input: %d\n", el_type);
				for (i=0;i< *num_tokens; i++)
					free(tokens[i]);
				free(tokens);
				free(tmp);
				return NULL;
			}

			if (el_type != SPACE)
			{
				// create a new token
				tmp[i] = '\0';
				tokens[*num_tokens] = (char *)malloc((strlen(tmp) + 2) * sizeof(char));
				strcpy(tokens[*num_tokens],tmp);
				token_types[*num_tokens] = el_type;
				*num_tokens = *num_tokens + 1;

				if (*num_tokens >= token_alloc)
				{
					tokens = (char **)realloc(tokens, token_alloc * 2 * sizeof(char *));
					token_types = (int *)realloc(token_types, token_alloc * 2 * sizeof(int));
					token_alloc *= 2;
				}
			}
		
			// then cycle the types
			el_type = el_type2;
			i=1;
			strcpy(tmp,&ch);
		}
		else
		{
			if (i == (tmpsize - 1))
			{
				tmpsize += GSTR_MAXSIZE;
				tmp = (char *)realloc(tmp,tmpsize * sizeof(char));
			}
			tmp[i] = ch;
			i++;
		}
	}

	free(tmp);
	return tokens;
}

int isNumber(char *str)
{
	int i = 0;
	int base = 10;
	int first_char_is_zero = 0;

	for (i = 0; i < strlen(str); i++)
	{
		if ((i == 0) && (str[i] == '0'))
		{
			first_char_is_zero = 1;
			continue;
		}
		if ((i == 1) && first_char_is_zero)
		{
			if (str[i] == 'b') base = 2;
			else if (str[i] == 'o') base = 8;
			else if (str[i] == 'd' || (str[i] >= '0' && str[i] <= '9')) base = 10;
			else if (str[i] == 'x') base = 16;
			else return 0;
			continue;
		}
		if ((base == 2) && !(str[i] >= '0' && str[i] <= '1')) return 0;
		if ((base == 8) && !(str[i] >= '0' && str[i] <= '7')) return 0;
		if ((base == 10) && !(str[i] >= '0' && str[i] <= '9')) return 0;
		if ((base == 16) && !(isdigit(str[i]) || (str[i] >= 'a' && str[i] <= 'f') || (str[i] >= 'A' && str[i] <= 'F'))) return 0;
	}
	return 1;
}

int isOperator(char *str)
{
	if (str[0] == '+' || str[0] == '-' || str[0] == '*' || str[0] == '/' || str[0] == '%' || str[0] == '^') return 1;
    if (str[0] == '!' || str[0] == '#' || str[0] == '=') return 1;
	return 0;
}

/* check to see if a string contains an operator */
int hasOperator(char *str)
{
	int i = 0;

	for (i = 0; i < strlen(str); i++)
		if (isOperator(&str[i]))
			return 1;

	return 0;
}

/* check to see if str contains + or - (Addition or Subtraction)*/
int hasOperatorAS(char *str)
{
	int i = 0;

	for (i = 0; i < strlen(str); i++)
		if (str[i] == '+' || str[i] == '-')
			return 1;

	return 0;
}

void get_expression(char* in, str_t* out)
{
    char* tok;
    char delim[] = { ' ', '\0' };
    int invalid_string = 0;
    bstack_t stk;
    str_t* tmp, * tmp1, * tmp2;

    // if there is no input to process, we are done... 
    if (in == NULL) return;

    stack_init(20, &stk, 0);
    tmp = (str_t*)malloc(sizeof(str_t));
    sInit(tmp);
    tmp1 = (str_t*)malloc(sizeof(str_t));
    sInit(tmp1);
    tmp2 = (str_t*)malloc(sizeof(str_t));
    sInit(tmp2);

    tok = strtok(in, delim);
    while (tok != NULL)
    {
        // if token is a number, push it onto the stack... 
        if (isNumber(tok))
        {
            toStr(tok, tmp);
            push(tmp, &stk);
        }
        else if (isOperator(tok))
        {
            if (tok[0] == '!' || tok[0] == '#')
            {
                // factorial and primorial are unary operators
                if (pop(tmp1, &stk) == 0)
                {
                    invalid_string = 1;
                    break;
                }
                if (hasOperator(tmp1->s))
                {
                    sClear(tmp);
                    sAppend("(", tmp);
                    sAppend(tmp1->s, tmp);
                    sAppend(")", tmp);
                    sAppend(tok, tmp);
                }
                else
                {
                    sClear(tmp);
                    sAppend(tmp1->s, tmp);
                    sAppend(tok, tmp);
                }
                push(tmp, &stk);
            }
            else if (tok[0] == '+')
            {
                // pop off two elements a,b
                // push back on "(aOPb)" 
                if ((pop(tmp2, &stk) == 0) || (pop(tmp1, &stk) == 0))
                {
                    invalid_string = 1;
                    break; // bad input string, don't parse any more...
                }
                sClear(tmp);
                sAppend(tmp1->s, tmp);
                sAppend(tok, tmp);
                sAppend(tmp2->s, tmp);
                push(tmp, &stk);
            }
            else if (tok[0] == '-')
            {
                // pop off two elements a,b
                //   push back on "(aOPb)" 
                if ((pop(tmp2, &stk) == 0) || (pop(tmp1, &stk) == 0))
                {
                    invalid_string = 1;
                    break; // bad input string, don't parse any more... 
                }
                sClear(tmp);
                sAppend(tmp1->s, tmp);
                sAppend(tok, tmp);
                if (hasOperatorAS(tmp2->s))
                {
                    sAppend("(", tmp);
                    sAppend(tmp2->s, tmp);
                    sAppend(")", tmp);
                }
                else
                {
                    sAppend(tmp2->s, tmp);
                }
                push(tmp, &stk);
            }
            else
            {
                // pop off two elements a,b
                // push back on "(aOPb)"
                if ((pop(tmp2, &stk) == 0) || (pop(tmp1, &stk) == 0))
                {
                    invalid_string = 1;
                    break; // bad input string, don't parse any more... 
                }
                sClear(tmp);
                if (hasOperator(tmp1->s))
                {
                    sAppend("(", tmp);
                    sAppend(tmp1->s, tmp);
                    sAppend(")", tmp);
                }
                else
                    sAppend(tmp1->s, tmp);

                sAppend(tok, tmp);

                if (hasOperator(tmp2->s))
                {
                    sAppend("(", tmp);
                    sAppend(tmp2->s, tmp);
                    sAppend(")", tmp);
                }
                else
                    sAppend(tmp2->s, tmp);

                push(tmp, &stk);
            }
        }
        else
        {
            // non-number and non-operator encountered, we're done grabbing the input string
            break;
        }
        tok = strtok(NULL, delim);
    }

    if (!invalid_string) pop(out, &stk);

    stack_free(&stk);
    sFree(tmp);
    free(tmp);
    sFree(tmp1);
    free(tmp1);
    sFree(tmp2);
    free(tmp2);
}

char* process_expression(char* input_exp, meta_t* metadata,
    int force_quiet, int no_convert_result)
{
    // process the expression in input_exp, which can be multi-line.
    // force_quiet will suppress all output.
    // no_convert_result will skip any conversion of the result back
    // to a string (e.g., for scripts that don't need to display 
    // intermediate results).
    str_t str;
    str_t* out;
    int num;
    int i;
    mpz_t tmp;
    char* outstr = NULL;

    mpz_init(tmp);
    sInit(&str);
    toStr(input_exp, &str);

    // multi-line statement blocks:
    // have the preprocessor check for open '{' while parsing
    // statement bodies.  If we see an open bracket without
    // a closing one, create a new variable with the name of
    // the block.  E.g. if parsing a 2nd for loop, the name
    // would be 'for1_block'.  Then we return, and get some more
    // text.  A global will need to keep track of open/closed
    // brackets to know whether to keep appending to a block 
    // variable or to finialize it and increment to the next 
    // statement.
    // multi-statement lines:
    // separate with commas?
    out = preprocess(&str, &num);

    for (i = 0; i < num; i++)
    {
        calc_with_assignment(&out[i], metadata, force_quiet);
        sFree(&out[i]);
    }

    // return the last result.  Return a new string
    // since the input string may not be big enough
    // to hold the result.
    if (!no_convert_result)
    {
        get_uvar("ans", tmp);
        outstr = mpz_get_str(NULL, OBASE, tmp);
    }

    mpz_clear(tmp);
    sFree(&str);
    free(out);
    return outstr;
}

void calc_with_assignment(str_t* in, meta_t* metadata, int force_quiet)
{
    char* ptr;
    char varname[80];
    int offset = 0;
    int nooutput;
    str_t str;
    mpz_t tmp;

    mpz_init(tmp);
    sInit(&str);

    if ((ptr = strchr(in->s, '=')) != NULL)
    {
        offset = (ptr - in->s);
        strncpy(varname, in->s, offset);
        varname[offset++] = '\0';
    }
    else
    {
        strcpy(varname, "ans");
    }

    // look for a trailing semicolon
    if (in->s[strlen(in->s) - 1] == ';')
    {
        nooutput = 1;
        in->s[strlen(in->s) - 1] = '\0';
    }
    else
    {
        nooutput = 0;
    }

    toStr(in->s + offset, &str);
    if (!calc(&str, metadata))
    {
        if (strcmp(str.s, "") != 0)
        {           
            // always set the default variable to the new answer
            mpz_set_str(tmp, str.s, 0);
            sCopy(&str, in);
            set_uvar("ans", tmp);

            // and optionally any assigned variable as well.
            if (set_uvar(varname, tmp))
            {
                new_uvar(varname, tmp);
            }

            if ((nooutput == 0) && (force_quiet == 0))
            {
                if (OBASE == DEC)
                {
                    gmp_printf("\n%s = %Zd\n\n", varname, tmp);

                }
                else if (OBASE == HEX)
                {
                    gmp_printf("\n%s = %Zx\n\n", varname, tmp);
                }
            }
        }
    }

    mpz_clear(tmp);
    sFree(&str);
    return;
}

int calc(str_t *in, meta_t *metadata)
{

	/*
	Dijkstra's shunting algorithm 

	While there are tokens to be read: 
		* Read a token. 
		* If the token is a number, then add it to the output queue. 
		* If the token is a function token, then push it onto the stack. 
		* If the token is a function argument separator (e.g., a comma): 
			* Until the topmost element of the stack is a left parenthesis, pop the element onto the 
				output queue. If no left parentheses are encountered, either the separator was misplaced 
				or parentheses were mismatched. 
		* If the token is an operator, o1, then: 
			* while there is an operator, o2, at the top of the stack, and either 
				o1 is associative or left-associative and its precedence is less than 
				(lower precedence) or equal to that of o2, or o1 is right-associative and its precedence 
				is less than (lower precedence) that of o2,
				* pop o2 off the stack, onto the output queue; 
			* push o1 onto the operator stack. 
		* If the token is a left parenthesis, then push it onto the stack. 
		* If the token is a right parenthesis: 
			* Until the token at the top of the stack is a left parenthesis, pop operators off the stack 
				onto the output queue. 
			* Pop the left parenthesis from the stack, but not onto the output queue. 
			* If the token at the top of the stack is a function token, pop it and onto the output queue. 
			* If the stack runs out without finding a left parenthesis, then there are mismatched parentheses. 
	* When there are no more tokens to read: 
		* While there are still operator tokens in the stack: 
			* If the operator token on the top of the stack is a parenthesis, then there are mismatched 
				parenthesis. 
			* Pop the operator onto the output queue. 
	Exit. 

    In addition, we have immediate operators to deal with.  for post immediates (i.e. 4!, 9#) add
	to the output queue, just like a number

	a couple test cases
	trial(4+5^78*(81345-4),21345)

	*/

	//////////////////////////////////////////////
	// as this is primarily a factoring engine frontend, create a factorization 
	// object and pass it around through all the feval calls.  Those functions
	// which do factoring will know what to do with it.  If the expression turns 
	// out to have nothing to do with factorization, then the object will not be
	// used and can just be freed at the end of this routine.
	// ///////////////////////////////////////////

	int i,retval,na,func,j,k;
    bstack_t stk;		//general purpose stack
	str_t *tmp;			//temporary str_t
	str_t *post;		//post fix expression
	char **tokens;		//pointer to an array of strings holding tokens
	char *tok;
	char *strN;			/* use this to make a copy of the post-fix string */
	char delim[2];
	int *token_types;	//type of each token
	int num_tokens;		//number of tokens in the array.
	int varstate;
	mpz_t tmpz;
    char *tok_context;

	retval = 0;

	//initialize and find tokens
	token_types = (int *)malloc(100 * sizeof(int));
	tokens = tokenize(in->s, token_types, &num_tokens);
	if (tokens == NULL)
	{		
		free(token_types);
		return 1;
	}

    stack_init(20, &stk, 0);
	tmp = (str_t *)malloc(sizeof(str_t));
	sInit(tmp);	
	mpz_init(tmpz);
	post = (str_t *)malloc(sizeof(str_t));
	sInit(post);

	// run the shunting algorithm
	i=0;
	post->s[0] = '\0';
	while (i<num_tokens)
	{
		switch (token_types[i])
		{
		case 3:
			// NUM
			// to num output queue
			sAppend(" ",post);
			sAppend(tokens[i],post);
			break;
		case 6:
			// LP
			toStr(tokens[i],tmp);
			push(tmp,&stk);
			break;
		case 7:
			// string (function or variable name)
			varstate = get_uvar(tokens[i],tmpz);
			if (varstate == 0)
			{
				// found a variable with that name, copy it's value
				// to num output queue
				sAppend(" ",post);
				sAppend(tokens[i],post);
			}
			else if (varstate == 2)
			{
				// do nothing, special case
			}
			else if (getFunc(tokens[i],&na) >= 0) 
			{
				// valid function, push it onto stack
				toStr(tokens[i],tmp);
				push(tmp,&stk);
			}
			else
			{
				// a non-numeric string that is not a variable or a function.
                // the only thing that makes sense is for it to be a string
                // representing a new assignment (new variable).  Assume that's the
                // case.  if not, errors will be raised further down the line.
                sAppend(" ", post);
                sAppend(tokens[i], post);
			}
			break;
		case 9:
			//comma (function argument separator)
			while (1)
			{
				if (pop(tmp,&stk) == 0)
				{
					// stack empty and we are still looking for a LP
					printf("bad function separator position or mismatched parens\n");
					retval = 1;
					goto free;
				}

				if (strcmp(tmp->s,"(") == 0)
				{
					// found a left paren.  put it back and continue
					push(tmp,&stk);
					break;
				}
				else
				{
					// copy to output operator queue
					sAppend(" ",post);
					sAppend(tmp->s,post);
				}
			}
			break;
		case 5:
			// right paren
			while (1)
			{
				if (pop(tmp,&stk) == 0)
				{
					// stack empty and we are still looking for a LP
					printf("mismatched parens\n");
					retval = 1;
					goto free;
				}

				if (strcmp(tmp->s,"(") == 0)
				{
					// found a left paren.  ignore it.
					if (pop(tmp,&stk) != 0)
					{
						// is the top of stack a function?
						if ((getFunc(tmp->s,&na) >= 0) && (strlen(tmp->s) > 1))
						{
							// the extra check for strlen > 1 fixes
							// the case where the string is an operator, not
							// a function.  for multichar operators this won't work
							// instead should probably separate out the operators
							// from the functions in getFunc

							// yes, put it on the output queue as well
							sAppend(" ",post);
							sAppend(tmp->s,post);
						}
						else
						{
							// no, put it back
							push(tmp,&stk);
						}
					}
					break;
				}
				else
				{
					// copy to output queue
					sAppend(" ",post);
					sAppend(tmp->s,post);
				}
			}
			break;
		case 4:
			// operator
			while (pop(tmp,&stk))
			{
				if (strlen(tmp->s) == 1 && getOP(tmp->s[0]) > 0)
				{
					// its an operator
					// check the precedence
					if (op_precedence(tmp->s,tokens[i],getAssoc(tmp->s)))
					{
						// push to output op queue
						sAppend(" ",post);
						sAppend(tmp->s,post);
					}
					else
					{
						// put the tmp one back and bail
						push(tmp,&stk);
						break;
					}
				}
				else
				{
					// its not an operator, put it back and bail
					push(tmp,&stk);
					break;
				}
			}
			// push the current op onto the stack.
			toStr(tokens[i],tmp);
			push(tmp,&stk);

			break;
		case 2:
			// post unary operator
			// I think we can jush push these into the output operator queue
			toStr(tokens[i],tmp);
			sAppend(" ",post);
			sAppend(tmp->s,post);
			break;
		}
		i++;
	}

	// now pop all operations left on the stack to the output queue
	while (pop(tmp,&stk))
	{
		if (strcmp(tmp->s,"(") == 0 || strcmp(tmp->s,")") == 0)
		{
			printf("mismatched parens\n");
			retval = 1;
			goto free;
		}
		sAppend(" ",post);
		sAppend(tmp->s,post);
	}

	// free the input tokens
	for (i=0;i<num_tokens;i++)
		free(tokens[i]);
	free(tokens);

	// process the output postfix expression:
	// this can be done with a simple stack
	// all tokens are separated by spaces
	// all tokens consist of numbers or functions

	// now evaluate the RPN expression
    if (CALC_VERBOSE)
	    printf("processing postfix expression: %s\n",post->s);
	delim[0] = ' ';
	delim[1] = '\0';

    // need to use strtok_s: some functions need to
    // do their own tokenizing so we need to remember
    // this one.
    tok_context = NULL;
	tok = strtok_s(post->s,delim,&tok_context);
	if (tok == NULL)
	{
        sClear(in);
		goto free;
	}

	do
	{
        if (CALC_VERBOSE)
        {
            printf("stack contents: ");
            for (i = 0; i < stk.num; i++)
                printf("%s ", stk.elements[i]->s);
            printf("\n");
            printf("current token: %s\n\n", &tok[0]);
        }

		switch (get_el_type2(tok[0]))
		{
		case NUM:
			toStr(tok,tmp);
			push(tmp,&stk);
			break;
		case AMBIG:
			// could be a number or a function
			// if the next character is a number, this is too
			// else its an operator
			if (get_el_type2(tok[1]) == NUM)
			{
				toStr(tok,tmp);
				push(tmp,&stk);
				break;
			}

			// if not a num, proceed into the next switch (function handle)
		default:
			func = getFunc(tok,&na);
            if (CALC_VERBOSE)
            {
                printf("processing function %d: %s from token %s\n", func, function_names[func], tok);
            }

			if (func >= 0)
			{
				// pop those args and put them in a global array
				for (j=0;j<na;j++)
				{
					// right now we must get all of the operands.
					// somewhere in there we should make allowances
					// for getting a reduced number (i.e. for unary "-"
					// and for variable numbers of arguments
					int r;
					k = pop(tmp,&stk);

                    if (k == 0)
                    {
                        // didn't get the expected number of arguments
                        // for this function.  This may be ok, if the
                        // function accepts varable argument lists.
                        // feval will handle it.
                        break;
                    }

					// try to make a number out of it
					r = mpz_set_str(tmpz, tmp->s, 0);
                    if (r < 0)
                    {
                        if (CALC_VERBOSE)
                        {
                            printf("adding %s to choperands\n", tmp->s);
                        }
                        strcpy(choperands[na - j - 1], tmp->s);
                    }
                    else
                    {
                        // it is a number, put it in the operand pile
                        mpz_set(operands[na - j - 1], tmpz);
                    }
				}

				na = j;
				// call the function evaluator with the 
				// operator string and the number of args available
				na = feval(func,na,metadata);

				// put result back on stack
				for (j=0;j<na;j++)
				{
					int sz = mpz_sizeinbase(operands[j], 10) + 10;
					if (tmp->alloc < sz)
					{
						tmp->s = (char *)realloc(tmp->s, sz * sizeof(char));
						tmp->alloc = sz;
					}
					mpz_get_str(tmp->s, 10, operands[j]);
					tmp->nchars = strlen(tmp->s)+1;
					push(tmp,&stk);
				}
			}
            else if (get_uvar(tok, tmpz) == 0)
            {
                int sz;

                // the string token is not a function, check if it's a defined variable

                sz = mpz_sizeinbase(tmpz, 10) + 10;
                if (gstr1.alloc < sz)
                {
                    gstr1.s = (char *)realloc(gstr1.s, sz * sizeof(char));
                    gstr1.alloc = sz;
                }
                mpz_get_str(gstr1.s, 10, tmpz);
                gstr1.nchars = strlen(gstr1.s) + 1;
                sCopy(tmp, &gstr1);
                push(tmp, &stk);
            }
            else if (is_strvar(tok))
            {
                // is a string variable... push it onto the stack
                toStr(tok, tmp);
                push(tmp, &stk);
            }
            else
			{
                printf("unrecognized variable or function '%s'\n", tok);
                sClear(in);
                goto free;
			}			
		}

        tok = strtok_s((char *)0, delim, &tok_context);
	} while (tok != NULL);
	pop(in,&stk);

free:	
	free(token_types);
	stack_free(&stk);
	mpz_clear(tmpz);
	sFree(tmp);
	free(tmp);
	sFree(post);
	free(post);
	return retval;
}

int getFunc(char *s, int *nargs)
{
	// return the opcode associated with the function, and
	// the number of arguments it takes
	int i,j;

	for (i = 0; i < NUM_FUNC; i++)
	{
		j = strcmp(function_names[i],s);
		if (j == 0)
		{
			*nargs = function_nargs[i];
			return i;
		}
	}

	return -1;
}

int check_args(int funcnum, int nargs)
{
    if (nargs != function_nargs[funcnum])
    {
        printf("wrong number of arguments in %s, expected %d\n", 
            function_names[funcnum], function_nargs[funcnum]);
        fflush(stdout);
        return 1;
    }
    else
    {
        return 0;
    }
}

int feval(int funcnum, int nargs, meta_t *metadata)
{
	// evaluate the function 'func', with 'nargs' argument(s) located
	// in the mpz_t array 'operands'.
	// place return values in operands[0]
	mpz_t mp1, mp2, mp3, tmp1, tmp2;
	mpz_t gmpz;

	str_t str;
	uint32_t i=0;
	uint64_t n64;
	uint32_t j,k;
	double t;
	struct timeval tstart, tstop;
	uint64_t lower, upper, inc, count;
    fact_obj_t* fobj = metadata->fobj;
    soe_staticdata_t* sdata = metadata->sdata;

	mpz_init(mp1);
	mpz_init(mp2);
	mpz_init(mp3);
	mpz_init(tmp1);
	mpz_init(tmp2);
	sInit(&str);

	switch (funcnum)
	{
	case 0:
		// fib - one argument
        if (check_args(funcnum, nargs)) break;
		mpz_fib_ui(operands[0], mpz_get_ui(operands[0]));
		break;
	case 1:
		// luc - one argument
        if (check_args(funcnum, nargs)) break;
		mpz_lucnum_ui(operands[0], mpz_get_ui(operands[0]));
		break;

	case 2:
		// expr - one argument
		// this is used to evaluate numerical expressions from the command line,
		// now that the default behavior is to factor the input.
		// since whatever was in the function will have already been evaluated,
		// simply return the input.
        if (check_args(funcnum, nargs)) break;

		break;

	case 3:
		// gcd - two arguments
        if (check_args(funcnum, nargs)) break;			
		mpz_gcd(operands[0], operands[0], operands[1]);

		break;
	case 4:
		// jacobi - two arguments
        if (check_args(funcnum, nargs)) break;

		if (mpz_odd_p(operands[1]))
			mpz_set_si(operands[0], mpz_jacobi(operands[0], operands[1]));
		else
			printf("jacobi defined only for odd denominators!\n");

		break;

	case 5:
		// rand - one argument
        if (check_args(funcnum, nargs)) break;

		mpz_set_ui(operands[1], 10);
		mpz_pow_ui(operands[1], operands[1], mpz_get_ui(operands[0]));
		mpz_urandomm(operands[0], gmp_randstate, operands[1]);
		break;
	case 6:
		// lg2 - one argument
        if (check_args(funcnum, nargs)) break;
		mpz_set_ui(operands[0], mpz_sizeinbase(operands[0], 2));
		break;
	case 7:
		// log - one argument
        if (check_args(funcnum, nargs)) break;
		mpz_set_ui(operands[0], mpz_sizeinbase(operands[0], 10));
		break;
	case 8:
		// ln - one argument
        if (check_args(funcnum, nargs)) break;
		mpz_set_ui(operands[0], (uint32_t)((mpz_sizeinbase(operands[0], 2)-1) * log(2.0)));
		break;
	case 9:
		// size - one argument
        if (check_args(funcnum, nargs)) break;

		printf("%d digits, %d bits\n", (int)mpz_sizeinbase(operands[0], 10),
			(int)mpz_sizeinbase(operands[0], 2));

		mpz_set_ui(operands[0], mpz_sizeinbase(operands[0], 2));

		break;
	case 10:
		// issquare
        if (check_args(funcnum, nargs)) break;
		i = mpz_perfect_square_p(operands[0]);
		mpz_set_ui(operands[0], i);
		break;
	case 11:
		// isprime - one argument
        if (check_args(funcnum, nargs)) break;
        i = mpz_probab_prime_p(operands[0], 3);
		mpz_set_ui(operands[0], i);
		break;
	case 12:
		// sqrt - one argument
        if (check_args(funcnum, nargs)) break;
		mpz_root(operands[0], operands[0], 2);
		break;
	case 13:
		// modinv - two arguments
        if (check_args(funcnum, nargs)) break;
		mpz_invert(operands[0], operands[0], operands[1]);        
		break;
	case 14:
		// modexp - three arguments
        if (check_args(funcnum, nargs)) break;
		mpz_powm(operands[0], operands[0], operands[1], operands[2]);
		break;
	case 15:
		// nroot - two arguments
        if (check_args(funcnum, nargs)) break;
		mpz_root(operands[0], operands[0], mpz_get_ui(operands[1]));
		break;
	case 16:
		// shift - two arguments
        if (check_args(funcnum, nargs)) break;

		if (mpz_sgn(operands[1]) >= 0)
			mpz_mul_2exp(operands[0], operands[0], mpz_get_ui(operands[1]));
		else
			mpz_tdiv_q_2exp(operands[0], operands[0], -1*mpz_get_si(operands[1]));

		break;
	case 17:
		// ispow - one argument
        if (check_args(funcnum, nargs)) break;
		i = mpz_perfect_power_p(operands[0]);
		mpz_set_ui(operands[0], i);
		break;
	case 18:
		// randb - one argument
        if (check_args(funcnum, nargs)) break;
		mpz_urandomb(operands[0], gmp_randstate, mpz_get_ui(operands[0]));
		break;
	case 19:
		// add
        if (check_args(funcnum, nargs)) break;
		mpz_add(operands[0], operands[0], operands[1]);
		break;
	case 20:
		//subtract or negate
		if (nargs == 1)
		{
			mpz_neg(operands[0], operands[0]);
		}
		else if (nargs == 2)
		{
			mpz_sub(operands[0], operands[0], operands[1]);
		}
		else
		{
			printf("wrong number of arguments in sub/neg\n");
			break;
		}

		break;

	case 21:
		// mul
        if (check_args(funcnum, nargs)) break;
		mpz_mul(operands[0], operands[0], operands[1]);
		break;
	case 22:
		// div
        if (check_args(funcnum, nargs)) break;
		mpz_tdiv_q(operands[0], operands[0], operands[1]);
		break;
	case 23:
		// !
        if (check_args(funcnum, nargs)) break;
		mpz_fac_ui(operands[0], mpz_get_ui(operands[0]));
		break;
	case 24:
		// primorial
        if (check_args(funcnum, nargs)) break;
		mpz_primorial_ui(operands[0], mpz_get_ui(operands[0]));
		break;
	case 25:
		// eq
        if (check_args(funcnum, nargs)) break;
		mpz_set_ui(operands[0], mpz_cmp(operands[0], operands[1]) == 0);
		break;
	case 26:
		// <<
        if (check_args(funcnum, nargs)) break;
		mpz_mul_2exp(operands[0], operands[0], mpz_get_ui(operands[1]));
		break;
	case 27:
		// >>
        if (check_args(funcnum, nargs)) break;
		mpz_tdiv_q_2exp(operands[0], operands[0], mpz_get_ui(operands[1]));
		break;
	case 28:
		// mod
        if (check_args(funcnum, nargs)) break;
		mpz_mod(operands[0], operands[0], operands[1]);
		break;
	case 29:
		// exp
        if (check_args(funcnum, nargs)) break;
		mpz_pow_ui(operands[0], operands[0], mpz_get_ui(operands[1]));
		break;
	case 30:
        // REDC (redc)
        if (check_args(funcnum, nargs)) break;

        // arguments are:
        // 1: T - input to reduce
        // 2: n - reduction modulus
        // 3: r - bits in R
        mpz_set_ui(mp1, 1);
        j = mpz_get_ui(operands[2]);
        mpz_mul_2exp(mp1, mp1, j);
        mpz_invert(mp2, operands[1], mp1);
        mpz_sub(mp2, mp1, mp2);
        mpz_mul(mp3, operands[0], mp2);
        mpz_tdiv_r_2exp(mp3, mp3, j);
        mpz_mul(mp3, mp3, operands[1]);
        mpz_add(mp3, mp3, operands[0]);
        mpz_tdiv_q_2exp(operands[0], mp3, j);
        if (mpz_cmp(operands[0], operands[1]) >= 0)
        {
            mpz_sub(operands[0], operands[0], operands[1]);
        }
		break;

	case 31:
		// bitxor
        if (check_args(funcnum, nargs)) break;
        mpz_xor(operands[0], operands[0], operands[1]);
		break;
	case 32:
		// bitand
        if (check_args(funcnum, nargs)) break;
		mpz_and(operands[0], operands[0], operands[1]);
		break;
	case 33:
		// bitor
        if (check_args(funcnum, nargs)) break;
		mpz_ior(operands[0], operands[0], operands[1]);
		break;
	case 34:
		// onecomp
        if (check_args(funcnum, nargs)) break;
		mpz_com(operands[0], operands[0]);
		break;
	case 35:
		// lte
        if (check_args(funcnum, nargs)) break;
		mpz_set_ui(operands[0], mpz_cmp(operands[0], operands[1]) <= 0);
		break;
	case 36:
		// gte
        if (check_args(funcnum, nargs)) break;
		mpz_set_ui(operands[0], mpz_cmp(operands[0], operands[1]) >= 0);
		break;
	case 37:
		// lt
        if (check_args(funcnum, nargs)) break;
		mpz_set_ui(operands[0], mpz_cmp(operands[0], operands[1]) < 0);
		break;
	case 38:
		// gt
        if (check_args(funcnum, nargs)) break;
		mpz_set_ui(operands[0], mpz_cmp(operands[0], operands[1]) > 0);
		break;
    case 39:
        // popcnt
        if (check_args(funcnum, nargs)) break;
        j = mpz_popcount(operands[0]);
        mpz_set_ui(operands[0], j);
        break;
    case 40:
        // nextprime
        if (check_args(funcnum, nargs)) break;
        mpz_nextprime(operands[0], operands[0]);
        break;
    case 41:
        // print
        if (check_args(funcnum, nargs)) break;

        if (OBASE == DEC)
            gmp_printf("%Zd\n", operands[0]);
        else if (OBASE == HEX)
            gmp_printf("%Zx\n", operands[0]);

        fflush(stdout);
        break;

    case 42:
        // lcm
        if (check_args(funcnum, nargs)) break;
        mpz_lcm(operands[0], operands[0], operands[1]);
        break;
    case 43:
        // exit
        exit(0);
        break;
	case 44:
		// abs
		if (check_args(funcnum, nargs)) break;
		mpz_abs(operands[0], operands[0]);
		break;
    case 45:
        // extgcd
        if (check_args(funcnum, nargs)) break;
        mpz_gcdext(operands[0], mp1, mp2, operands[0], operands[1]);
        if (OBASE == DEC)
            gmp_printf("a*s + b*t = gcd\ns = %Zd\nt = %Zd\n", mp1, mp2);
        else
            gmp_printf("a*s + b*t = gcd\ns = %Zx\nt = %Zx\n", mp1, mp2);
        break;
    case 46:
        // fac2
        if (check_args(funcnum, nargs)) break;
        mpz_2fac_ui(operands[0], mpz_get_ui(operands[0]));
        break;
    case 47:
        // facm
        if (check_args(funcnum, nargs)) break;
        mpz_mfac_uiui(operands[0], mpz_get_ui(operands[0]), mpz_get_ui(operands[1]));
        break;
    case 48:
        // binom
        if (check_args(funcnum, nargs)) break;
        mpz_bin_ui(operands[0], operands[0], mpz_get_ui(operands[1]));
        break;
    case 49:
        // randp
        if (check_args(funcnum, nargs)) break;
        mpz_urandomb(operands[0], gmp_randstate, mpz_get_ui(operands[0]));
        mpz_nextprime(operands[0], operands[0]);
        break;
    case 50:
        // hamdist
        if (check_args(funcnum, nargs)) break;
        j = mpz_hamdist(operands[0], operands[1]);
        mpz_set_ui(operands[0], j);
        break;
    case 51:
        // snfs - two arguments
        if (check_args(funcnum, nargs)) break;

        // the first argument is the full input form and the second is 
        // the cofactor we use as the input.
        // fill the job's 'n' parameter now, and nfs will detect the form (or 
        // bail).
        fobj->nfs_obj.snfs = 1;
        mpz_set(fobj->N, operands[1]);
        mpz_set(fobj->nfs_obj.snfs_cofactor, operands[1]);

        if (mpz_sizeinbase(fobj->nfs_obj.snfs_cofactor, 10) < fobj->nfs_obj.min_digits)
        {
            printf("***** warning: possibly malformed co-factor (too small)!\n");
            printf("*****          co-factor expected to be input divided by known factors.\n");
            printf("*****          attempting snfs anyway.\n");
        }

        mpz_set(fobj->nfs_obj.gmp_n, operands[0]);
        nfs(fobj);
        mpz_set(operands[0], fobj->nfs_obj.gmp_n);
        print_factors(fobj->factors, fobj->N, fobj->VFLAG, fobj->NUM_WITNESSES);

        break;
    case 52:
        // rsa - one argument
        if (check_args(funcnum, nargs)) break;

        mpz_set_ui(mp1, 2048);
        mpz_set_ui(mp2, 4096);
        if (mpz_cmp(operands[0], mp2) > 0)
        {
            printf("bitlength too large");
            mpz_set_ui(operands[0], 1);
            break;
        }
        else if (mpz_cmp(operands[0], mp1) > 0)
        {
            printf("Paranoid, huh?  This might take a minute\n");
        }

        build_RSA(mpz_get_ui(operands[0]), operands[0], gmp_randstate);
        break;
    case 53:
        // factor - one argument
        if (check_args(funcnum, nargs)) break;

        mpz_set(fobj->N, operands[0]);
        factor(fobj);
        mpz_set(operands[0], fobj->N);
        print_factors(fobj->factors, fobj->N, fobj->VFLAG, fobj->NUM_WITNESSES);

        {
            char vname[20];
            for (i = 0; i < fobj->factors->num_factors; i++)
            {

                sprintf(vname, "_f%d", i);
                if (set_uvar(vname, fobj->factors->factors[i].factor))
                    new_uvar(vname, fobj->factors->factors[i].factor);
                sprintf(vname, "_fpow%d", i);
                mpz_set_ui(operands[4], fobj->factors->factors[i].count);
                if (set_uvar(vname, operands[4]))
                    new_uvar(vname, operands[4]);
            }
            sprintf(vname, "_fnum");
            mpz_set_ui(operands[4], fobj->factors->num_factors);
            if (set_uvar(vname, operands[4]))
                new_uvar(vname, operands[4]);
        }

        reset_factobj(fobj);

        break;
    case 54:
        // pm1 - one argument
        if (check_args(funcnum, nargs)) break;
        mpz_set(fobj->N, operands[0]);
        mpz_set(fobj->pm1_obj.gmp_n, operands[0]);
        pollard_loop(fobj);
        mpz_set(operands[0], fobj->pm1_obj.gmp_n);
        print_factors(fobj->factors, fobj->N, fobj->VFLAG, fobj->NUM_WITNESSES);
        break;
    case 55:
        // pp1 - two arguments, one optional
        mpz_set(fobj->N, operands[0]);
        if (nargs == 2)
        {
            mpz_set(fobj->pp1_obj.gmp_n, operands[0]);
            fobj->pp1_obj.numbases = mpz_get_ui(operands[1]);
            williams_loop(fobj);
            mpz_set(operands[0], fobj->pp1_obj.gmp_n);
        }
        else if (nargs == 1)
        {
            mpz_set(fobj->pp1_obj.gmp_n, operands[1]);
            fobj->pp1_obj.numbases = 1;
            williams_loop(fobj);
            mpz_set(operands[0], fobj->pp1_obj.gmp_n);
        }
        else
        {
            printf("wrong number of arguments in pp1\n");
            break;
        }
        print_factors(fobj->factors, fobj->N, fobj->VFLAG, fobj->NUM_WITNESSES);
        break;
    case 56:
        // rho - one argument
        if (check_args(funcnum, nargs)) break;

        mpz_set(fobj->N, operands[0]);
        mpz_set(fobj->rho_obj.gmp_n, operands[0]);
        brent_loop(fobj);
        mpz_set(operands[0], fobj->rho_obj.gmp_n);
        print_factors(fobj->factors, fobj->N, fobj->VFLAG, fobj->NUM_WITNESSES);
        break;
    case 57:
        // trial - two arguments
        mpz_set(fobj->N, operands[0]);
        if (nargs == 2)
        {
            mpz_set(fobj->div_obj.gmp_n, operands[0]);
            fobj->div_obj.print = 0;
            fobj->div_obj.limit = mpz_get_ui(operands[1]);
            zTrial(fobj);
            mpz_set(operands[0], fobj->div_obj.gmp_n);
        }
        else if (nargs == 1)
        {
            printf("using default trial division bound of 10000\n");
            mpz_set(fobj->div_obj.gmp_n, operands[1]);
            fobj->div_obj.print = 0;
            fobj->div_obj.limit = 10000;
            zTrial(fobj);
            mpz_set(operands[0], fobj->div_obj.gmp_n);
        }
        else
        {
            printf("wrong number of arguments in trial\n");
            break;
        }

        print_factors(fobj->factors, fobj->N, fobj->VFLAG, fobj->NUM_WITNESSES);
        break;
    case 58:
        // shanks - one argument
        if (check_args(funcnum, nargs)) break;
        mpz_init(gmpz);

        n64 = sp_shanks_loop(operands[0], fobj);
        print_factors(fobj->factors, fobj->N, fobj->VFLAG, fobj->NUM_WITNESSES);
        mpz_set_64(operands[0], n64);
        break;
    case 59:
        // siqs - one argument
        if (check_args(funcnum, nargs)) break;

        mpz_set(fobj->N, operands[0]);
        mpz_set(fobj->qs_obj.gmp_n, operands[0]);
        SIQS(fobj);
        mpz_set(operands[0], fobj->qs_obj.gmp_n);
        print_factors(fobj->factors, fobj->N, fobj->VFLAG, fobj->NUM_WITNESSES);
        break;

    case 60:
        // primes
        if (nargs == 2)
        {
            gettimeofday(&tstart, NULL);
            lower = mpz_get_64(operands[1]);
            upper = mpz_get_64(operands[2]);

            {
                soe_staticdata_t* sdata = soe_init(fobj->VFLAG, fobj->THREADS, 32768);
                soe_wrapper(sdata, lower, upper, 1, &n64, 0, 0);
                soe_finalize(sdata);
            }
            

            mpz_set_64(operands[0], n64);
            gettimeofday(&tstop, NULL);
            t = ytools_difftime(&tstart, &tstop);

            printf("elapsed time = %6.4f\n", t);
            break;
        }
        else if (nargs < 2)
        {
            printf("not enough arguments, please specify min and max of range\n");
            break;
        }
        else if (nargs > 3)
        {
            printf("wrong number of arguments in primes\n");
            break;
        }

        lower = mpz_get_64(operands[0]);
        upper = mpz_get_64(operands[1]);
        
        {
            uint64_t* PRIMES, NUM_P, P_MAX, P_MIN;
            soe_staticdata_t* sdata = soe_init(fobj->VFLAG, fobj->THREADS, 32768);
            PRIMES = soe_wrapper(sdata, lower, upper, mpz_get_ui(operands[2]), &NUM_P, 
                metadata->options->pfile, metadata->options->pscreen);

            if (PRIMES != NULL)
            {
                P_MIN = PRIMES[0];
                P_MAX = PRIMES[NUM_P - 1];
            }
            soe_finalize(sdata);
            mpz_set_64(operands[0], NUM_P);
        }

        break;
    case 61:
        // torture - two arguments
        if (check_args(funcnum, nargs)) break;

        i = mpz_get_ui(operands[0]);
        k = mpz_get_ui(operands[1]);
        for (j = 0; j < i; j++)
        {
            printf("***********RUN %d***********\n", j + 1);
            mpz_urandomb(operands[2], gmp_randstate, k);
            mpz_set(fobj->N, operands[2]);
            factor(fobj);
            mpz_set(operands[0], fobj->N);
            print_factors(fobj->factors, fobj->N, fobj->VFLAG, fobj->NUM_WITNESSES);
            clear_factor_list(fobj);
        }

        break;
    case 62:
        // ecm - two arguments
        mpz_set(fobj->N, operands[0]);
        if (nargs == 2)
        {
            k = mpz_get_ui(operands[1]);
            fobj->ecm_obj.num_curves = k;
            mpz_set(fobj->ecm_obj.gmp_n, operands[0]);
            ecm_loop(fobj);
            mpz_set(operands[0], fobj->ecm_obj.gmp_n);
        }
        else if (nargs == 1)
        {
            fobj->ecm_obj.num_curves = 1;
            mpz_set(fobj->ecm_obj.gmp_n, operands[1]);
            ecm_loop(fobj);
            mpz_set(operands[0], fobj->ecm_obj.gmp_n);
        }
        else
        {
            printf("wrong number of arguments in ecm\n");
            break;
        }

        print_factors(fobj->factors, fobj->N, fobj->VFLAG, fobj->NUM_WITNESSES);
        break;
    case 63:
        // lucas lehmer test
        if (check_args(funcnum, nargs)) break;

        if (llt(mpz_get_ui(operands[0]), fobj->VFLAG))
        {
            if (fobj->VFLAG > 1)
                gmp_printf("%Zd is prime!\n", operands[0]);
            mpz_set_ui(operands[0], 1);
        }
        else
        {
            if (fobj->VFLAG > 1)
                gmp_printf("%Zd is composite.\n", operands[0]);
            mpz_set_ui(operands[0], 0);
        }
        break;

    case 64:
        // siqsbench
        if (check_args(funcnum, nargs)) break;
        siqsbench(fobj);
        break;

    case 65:
        // sigma - sum of divisors function
        if (check_args(funcnum, nargs)) break;

        int oldvflag = fobj->VFLAG;
        if (mpz_sizeinbase(operands[0], 2) < 192)
        {
            fobj->VFLAG = -1;
        }

        mpz_set(mp2, operands[0]);
        mpz_set(fobj->N, operands[0]);
        k = mpz_get_ui(operands[1]);
        factor(fobj);

        mpz_set_ui(mp2, 1);
        for (i = 0; i < fobj->factors->num_factors; i++)
        {
            mpz_set_ui(mp1, 1);
            for (j = 1; j <= fobj->factors->factors[i].count; j++)
            {
                mpz_pow_ui(mp3, fobj->factors->factors[i].factor, j * k);
                mpz_add(mp1, mp1, mp3);
            }
            mpz_mul(mp2, mp2, mp1);
        }
        mpz_set(operands[0], mp2);
        fobj->VFLAG = oldvflag;

        break;
    case 66:
        // Euler's totient function
        if (check_args(funcnum, nargs)) break;

        oldvflag = fobj->VFLAG;
        if (mpz_sizeinbase(operands[0], 2) < 192)
        {
            fobj->VFLAG = -1;
        }

        mpz_set(mp2, operands[0]);
        mpz_set(fobj->N, operands[0]);
        factor(fobj);
        mpz_set(operands[0], mp2);
        mpz_tdiv_q(mp1, mp2, fobj->factors->factors[0].factor);
        mpz_sub(operands[0], operands[0], mp1);

        for (i = 1; i < fobj->factors->num_factors; i++)
        {
            mpz_tdiv_q(mp1, mp2, fobj->factors->factors[i].factor);
            mpz_sub(mp1, mp2, mp1);
            mpz_mul(operands[0], operands[0], mp1);
            mpz_tdiv_q(operands[0], operands[0], mp2);
        }
        fobj->VFLAG = oldvflag;

        break;

    case 67:
        // smallmpqs - 1 argument
        if (check_args(funcnum, nargs)) break;

        mpz_set(fobj->N, operands[0]);
        mpz_set(fobj->qs_obj.gmp_n, operands[0]);
        if (strlen(fobj->flogname) > 0)
        {
            fobj->logfile = fopen(fobj->flogname, "a");
            if (fobj->logfile == NULL)
                printf("fopen error: %s\n", strerror(errno));
        }
        smallmpqs(fobj);
        if (strlen(fobj->flogname) > 0)
        {
            if (fobj->logfile != NULL)
                fclose(fobj->logfile);
        }
        mpz_set(operands[0], fobj->qs_obj.gmp_n);
        print_factors(fobj->factors, fobj->N, fobj->VFLAG, fobj->NUM_WITNESSES);

        break;
    case 68:
        // testrange - 4 arguments (low, high, depth, witnesses)
        if (check_args(funcnum, nargs)) break;

        // move to soe library...
        {
            uint64_t num_found;
            uint64_t* primes;
            uint64_t range;
            mpz_t lowz, highz;

            primes = NULL;

            mpz_init(lowz);
            mpz_init(highz);
            mpz_set(lowz, operands[0]);
            mpz_set(highz, operands[1]);
            primes = sieve_to_depth(sdata, lowz, highz,
                0, mpz_get_ui(operands[3]), mpz_get_ui(operands[2]), 
                &num_found, fobj->options->pscreen, fobj->options->pfile);

            if (!NULL)
                free(primes);

            mpz_clear(lowz);
            mpz_clear(highz);
            mpz_set_ui(operands[0], num_found);
        }

        break;

    case 69:
        // sieverange - 4 arguments
        // range lo, range hi, sieve bound, count?
        if (check_args(funcnum, nargs)) break;

        // move to soe library...
        {
            uint64_t num_found;
            uint64_t* primes;
            uint64_t range;
            mpz_t lowz, highz;

            primes = NULL;

            mpz_init(lowz);
            mpz_init(highz);
            mpz_set(lowz, operands[0]);
            mpz_set(highz, operands[1]);
            primes = sieve_to_depth(sdata, lowz, highz,
                0, 0, mpz_get_ui(operands[2]),
                &num_found, fobj->options->pscreen, fobj->options->pfile);

            if (!NULL)
                free(primes);

            mpz_clear(lowz);
            mpz_clear(highz);
            mpz_set_ui(operands[0], num_found);
        }

        break;

    case 70:
        // fermat - three arguments
        if (nargs == 2)
        {
            mpz_set(fobj->N, operands[1]);
            mpz_set(fobj->div_obj.gmp_n, operands[1]);
            n64 = mpz_get_64(operands[2]);
            j = 1;
        }
        else if (nargs == 3)
        {
            mpz_set(fobj->N, operands[0]);
            mpz_set(fobj->div_obj.gmp_n, operands[0]);
            n64 = mpz_get_64(operands[1]);
            j = mpz_get_ui(operands[2]);
        }
        else
        {
            printf("wrong number of arguments in fermat\n");
            break;
        }

        zFermat(n64, j, fobj);
        mpz_set(operands[0], fobj->div_obj.gmp_n);
        print_factors(fobj->factors, fobj->N, fobj->VFLAG, fobj->NUM_WITNESSES);
        break;

    case 71:
        // nfs - one argument
        if (check_args(funcnum, nargs)) break;

        mpz_set(fobj->N, operands[0]);
        mpz_set(fobj->nfs_obj.gmp_n, operands[0]);
        nfs(fobj);
        mpz_set(operands[0], fobj->nfs_obj.gmp_n);
        print_factors(fobj->factors, fobj->N, fobj->VFLAG, fobj->NUM_WITNESSES);
        break;

    case 72:
        // tune, no arguments
        if (check_args(funcnum, nargs)) break;

        factor_tune(fobj);
        break;

    case 73:
        // bpsw, 1 argument
        if (check_args(funcnum, nargs)) break;

        i = mpz_bpsw_prp(operands[0]);
        // should we print out the input number again?
        if (i == PRP_COMPOSITE)
        {
            printf("Input is composite.  C%d\n", gmp_base10(operands[0]));
        }
        else if (i == PRP_PRP)
        {
            printf("Input is prp.  PRP%d\n", gmp_base10(operands[0]));
        }
        else if (i == PRP_PRIME)
        {
            printf("Input is prime.  P%d\n", gmp_base10(operands[0]));
        }

        break;

    case 74:
        // aprcl, one argument
        if (check_args(funcnum, nargs)) break;

        i = mpz_aprtcle(operands[0], APRTCLE_VERBOSE1);
        // should we print out the input number again?
        if (i == APRTCLE_COMPOSITE)
        {
            if (mpz_bpsw_prp(operands[0]) != PRP_COMPOSITE)
            {
                printf("\n");
                // if BPSW doesn't think this composite number is actually composite, then ask the user
                // to report this fact to the YAFU sub-forum at mersenneforum.org
                printf(" *** ATTENTION: BPSW issue found.  Please report the following number to:\n");
                printf(" *** ATTENTION: http://www.mersenneforum.org/forumdisplay.php?f=96\n");
                gmp_printf(" *** ATTENTION: n = %Zd\n", operands[0]);
            }
            printf("\nInput is composite.  C%d\n", gmp_base10(operands[0]));
        }
        else if (i == APRTCLE_PRP)
        {
            printf("\nInput is prp.  PRP%d\n", gmp_base10(operands[0]));
        }
        else if (i == APRTCLE_PRIME)
        {
            printf("\nInput is prime.  P%d\n", gmp_base10(operands[0]));
        }

        break;

    case 75:
        // semiprime list generation: two arguments (count, bits)
        if (check_args(funcnum, nargs)) break;

        //generate_semiprime_list(mpz_get_ui(operands[0]), mpz_get_ui(operands[1]), gmp_randstate);
        {
            // generate a list of 'num' semiprimes, each of size 'bits'
            // save to semiprimes.dat

            FILE* out;
            int i;
            mpz_t tmp1, tmp2, tmp3;
            uint32_t num = mpz_get_ui(operands[0]);
            uint32_t bits = mpz_get_ui(operands[1]);

            mpz_init(tmp1);
            mpz_init(tmp2);
            mpz_init(tmp3);

            out = fopen("semiprimes.dat", "w");
            if (out == NULL)
            {
                printf("couldn't open semiprimes.dat for writing\n");
                break;
            }

            for (i = 0; i < num; i++)
            {
                mpz_urandomb(tmp3, gmp_randstate, bits / 2);
                mpz_setbit(tmp3, bits / 2 - 1);
                mpz_nextprime(tmp1, tmp3);
                if (bits & 1)
                {
                    mpz_urandomb(tmp3, gmp_randstate, bits / 2 + 1);
                    mpz_setbit(tmp3, bits / 2);
                }
                else
                {
                    mpz_urandomb(tmp3, gmp_randstate, bits / 2);
                    mpz_setbit(tmp3, bits / 2 - 1);
                }
                mpz_nextprime(tmp2, tmp3);
                mpz_mul(tmp3, tmp2, tmp1);
                gmp_fprintf(out, "%Zd,%Zd,%Zd\n",
                    tmp3, tmp1, tmp2);
            }
            fclose(out);

            printf("generated %d semiprimes in file semiprimes.dat\n", num);
            mpz_clear(tmp1);
            mpz_clear(tmp2);
            mpz_clear(tmp3);
        }

        break;

	default:
		printf("unrecognized function code\n");
		mpz_set_ui(operands[0], 0);
		break;
	}

	sFree(&str);
	mpz_clear(mp1);
	mpz_clear(mp2);
	mpz_clear(mp3);
	mpz_clear(tmp1);
	mpz_clear(tmp2);
	return 1;
}

int new_uvar(const char *name, mpz_t data)
{
	int i;
	// create a new user variable with name 'name', and return
	// its location in the global uvars structure
	if (uvars.num == uvars.alloc)
	{
		//need more room for variables
		uvars.vars = (uvar_t *)realloc(uvars.vars, uvars.num * 2 * sizeof(uvar_t));
		uvars.alloc *= 2;
		for (i=uvars.num;i<uvars.alloc;i++)
			mpz_init(uvars.vars[i].data);
	}

	strcpy(uvars.vars[uvars.num].name,name);
	mpz_set(uvars.vars[uvars.num].data, data);
	uvars.num++;
	return uvars.num - 1;
}

int set_uvar(const char *name, mpz_t data)
{
	// look for 'name' in the global uvars structure
	// if found, copy in data and return 0
	// else return 1
	int i;

	i = mpz_get_ui(data);
	// first look if it is a global constant
	if (strcmp(name,"IBASE") == 0)
	{
		if (i != DEC && i != HEX && i != BIN && i != OCT)
		{
			printf("unknown base\n");
			return 1;
		}
		else
		{
			IBASE = i;
			return 0;
		}
	}
	else if (strcmp(name,"OBASE") == 0)
	{
		if (i != DEC && i != HEX && i != BIN && i != OCT)
		{
			printf("unknown base\n");
			return 0;
		}
		else
		{
			OBASE = i;
			return 1;
		}
	}

	for (i=0;i<uvars.num;i++)
	{
		if (strcmp(uvars.vars[i].name,name) == 0)
		{
			mpz_set(uvars.vars[i].data, data);
			return 0;
		}
	}
	return 1;
}

int get_uvar(const char *name, mpz_t data)
{
	// look for 'name' in the global uvars structure
	// if found, copy out data and return 0
	// else return 1 if not found
	int i;

	// first look if it is a global constant
	if (strcmp(name,"IBASE") == 0) {
		mpz_set_ui(data, IBASE); return 0;}
	else if (strcmp(name,"OBASE") == 0) {
		mpz_set_ui(data, OBASE); return 0;}
	

	for (i=0;i<uvars.num;i++)
	{
		if (strcmp(uvars.vars[i].name,name) == 0)
		{
			mpz_set(data, uvars.vars[i].data);
			return 0;
		}
	}

	if (strcmp(name,"vars") == 0) {
		printf("dumping variable name data:\n");
		printf("IBASE              %u\n",IBASE);
		printf("OBASE              %u\n",OBASE);		

		for (i=0;i<uvars.num;i++)
			printf("%s      %s\n",uvars.vars[i].name,mpz_get_str(gstr1.s, 10, uvars.vars[i].data));

		return 2;
	}

    if (strcmp(name, "strvars") == 0) {
        printf("dumping string variable name data:\n");

        for (i = 0; i<strvars.num; i++)
            printf("%s      %s\n", strvars.vars[i].name, strvars.vars[i].data);

        return 2;
    }

	return 1;
}

void free_uvars()
{
	int i;
	for (i=0;i<uvars.alloc;i++)
		mpz_clear(uvars.vars[i].data);
	free(uvars.vars);
}

int new_strvar(const char *name, char *data)
{
    int i;
    // create a new user variable with name 'name', and return
    // its location in the global uvars structure
    if (strvars.num == strvars.alloc)
    {
        // need more room for variables
        strvars.vars = (strvar_t *)realloc(strvars.vars, strvars.num * 2 * sizeof(strvar_t));
        strvars.alloc *= 2;
        for (i = strvars.num; i<strvars.alloc; i++)
            strvars.vars[i].data = (char *)malloc(GSTR_MAXSIZE * sizeof(char));
    }

    strcpy(strvars.vars[strvars.num].name, name);
    strcpy(strvars.vars[strvars.num].data, data);
    strvars.num++;
    return strvars.num - 1;
}

int set_strvar(const char *name, char *data)
{
    // look for 'name' in the global uvars structure
    // if found, copy in data and return 0
    // else return 1
    int i;

    for (i = 0; i<strvars.num; i++)
    {
        if (strcmp(strvars.vars[i].name, name) == 0)
        {
            strcpy(strvars.vars[i].data, data);
            return 0;
        }
    }
    return 1;
}

int is_strvar(const char *name)
{
    // look for 'name' in the global uvars structure
    // if found, return 1
    // else return 0 if not found
    int i;

    for (i = 0; i<strvars.num; i++)
    {
        if (strcmp(strvars.vars[i].name, name) == 0)
        {
            return 1;
        }
    }
    return 0;
}

char * get_strvarname(const char *data)
{
    // look for 'data' in the global uvars structure
    // if found, return the name of the variable
    int i;
    char *name = NULL;

    for (i = 0; i<strvars.num; i++)
    {
        if (strcmp(strvars.vars[i].data, data) == 0)
        {
            name = strvars.vars[i].name;
            break;
        }
    }
    return name;
}

int get_strvar(const char *name, char *data)
{
    // look for 'name' in the global uvars structure
    // if found, copy out data and return 0
    // else return 1 if not found
    int i;

    for (i = 0; i<strvars.num; i++)
    {
        if (strcmp(strvars.vars[i].name, name) == 0)
        {
            strcpy(data, strvars.vars[i].data);
            return 0;
        }
    }

    if (strcmp(name, "strvars") == 0) {
        printf("dumping string variable name data:\n");

        for (i = 0; i<strvars.num; i++)
            printf("%s      %s\n", strvars.vars[i].name, strvars.vars[i].data);

        return 2;
    }

    return 1;
}

void free_strvars()
{
    int i;
    for (i = 0; i < strvars.alloc; i++)
    {
        free(strvars.vars[i].data);
    }
    free(strvars.vars);
}

int invalid_dest(char* dest)
{
    //return 1 if invalid, 0 otherwise
    int i;

    if (getFunc(dest, &i) >= 0)
        return 1;	//is a function name

    //global vars are ok
    if (strcmp(dest, "POLLARD_STG1_MAX") == 0) {
        return 0;
    }
    else if (strcmp(dest, "POLLARD_STG2_MAX") == 0) {
        return 0;
    }
    else if (strcmp(dest, "WILL_STG1_MAX") == 0) {
        return 0;
    }
    else if (strcmp(dest, "WILL_STG2_MAX") == 0) {
        return 0;
    }
    else if (strcmp(dest, "ECM_STG1_MAX") == 0) {
        return 0;
    }
    else if (strcmp(dest, "ECM_STG2_MAX") == 0) {
        return 0;
    }
    else if (strcmp(dest, "BRENT_MAX_IT") == 0) {
        return 0;
    }
    else if (strcmp(dest, "IBASE") == 0) {
        return 0;
    }
    else if (strcmp(dest, "OBASE") == 0) {
        return 0;
    }
    else if (strcmp(dest, "QS_DUMP_CUTOFF") == 0) {
        return 0;
    }
    else if (strcmp(dest, "NUM_WITNESSES") == 0) {
        return 0;
    }
    else if (strcmp(dest, "LOGFLAG") == 0) {
        return 0;
    }
    else if (strcmp(dest, "VFLAG") == 0) {
        return 0;
    }
    else if (strcmp(dest, "PRIMES_TO_FILE") == 0) {
        return 0;
    }
    else if (strcmp(dest, "PRIMES_TO_SCREEN") == 0) {
        return 0;
    }

    //check starting char not lower case letter or _ or `
    if ((dest[0] < 95) || (dest[0] > 122) || (dest[0] == 96)) return 1;

    // check that dest string doesn't contain any invalid characters.
    // we allow non-leading characters to be a-z,A-Z,0-9,_
    for (i = 1; i < strlen(dest); i++)
    {
        if ((dest[i] < 48) || (dest[i] > 122) ||
            ((dest[i] > 90) && (dest[i] < 95)) ||
            ((dest[i] > 57) && (dest[i] < 65)) ||
            (dest[i] == 96))
        {
            return 1;
        }
    }

    return 0;
}