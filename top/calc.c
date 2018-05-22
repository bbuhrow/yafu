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

/*
Implements an arbitrary precision calculator.
*/

#include "yafu.h"
#include "soe.h"
#include "calc.h"
#include "yafu_stack.h"
#include "util.h"
#include "yafu_ecm.h"
#include "qs.h"
#include "factor.h"
#include "common.h"
#include "mpz_aprcl.h"

#define CALC_VERBOSE 0

char opchar[9] = { '=', '<', '>', '+', '-', '*', '/', '%', '^' }; // , '='};
char imms[3] = {'!','#','-'};
const int numopchars = 9;
char choperands[5][80];
mpz_t operands[5];
int for_cnt = 0;
int forp_cnt = 0;
int if_cnt = 0;

void calc_with_assignment(str_t *in, fact_obj_t *fobj, int force_quiet);

void reset_preprocessor(void) {
    for_cnt = 0;
    forp_cnt = 0;
    if_cnt = 0;
    return;
}

int calc_init()
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

    for (i = 0; i<5; i++)
        mpz_init(operands[i]);

    reset_preprocessor();
	return 1;
}

void calc_finalize()
{
    int i;
	free_uvars();
    free_strvars();
    for (i = 0; i<5; i++)
        mpz_clear(operands[i]);
}

int get_el_type2(char s)
{
	//there are several types of characters in an expression.  
	//decide which type this is
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
	for (i=0;i<3;i++)
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

	for (i=0;i<numopchars;i++)
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

int preprocess(str_t *in)
{
	//preprocess the expression in 'in'
	int i,j;
    char *ptr;
    char str[80];

	j=0;
	//remove white space
	for (i=0;i<in->nchars;i++)
	{
		if (isspace(in->s[i]))
			continue;
		
		in->s[j] = in->s[i];
		j++;
	}
	in->s[j] = '\0';

	// algebraic simplification (this would be cool...)

    // substitute loop/if body.
    // loop and if statements have bodies that may contain
    // expressions (assignments).  One way to handle these
    // correctly is to store the entire body in a 
    // temporary variable string for later processing.
    // this way, the existing expression evaluator doesn't
    // need to change.  Down in the feval function that will
    // evaluate the if/loop we can recover the body text and
    // recursively call process_expression.
    if ((ptr = strstr(in->s, "for(")) != NULL)
    {
        // new for loop
        char pre[8], vname[20];
        sprintf(pre, "for%d", for_cnt++);

        // tokenize the loop
        ptr = strtok(&in->s[4], ";");
        if (ptr == NULL)
        {
            printf("badly formatted for loop: for(init; test; iter; body)\n");
            exit(3);
        }
        sprintf(vname, "%s_init", pre);
        if (set_strvar(vname, ptr))
            new_strvar(vname, ptr);
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
        ptr = strtok(NULL, "\0");
        if (ptr == NULL)
        {
            printf("badly formatted for loop: for(init; test; iter; body)\n");
            exit(3);
        }
        strncpy(str, ptr, strlen(ptr) - 1);
        str[strlen(ptr) - 1] = '\0';
        sprintf(vname, "%s_body", pre);
        if (set_strvar(vname, str))
            new_strvar(vname, str);

        sprintf(in->s, "for(%s_init, %s_test, %s_iter, %s_body);",pre, pre, pre, pre);
    }

    if ((ptr = strstr(in->s, "forprime(")) != NULL)
    {
        // new for loop
        char pre[8], vname[20];
        sprintf(pre, "forp%d", forp_cnt++);

        // tokenize the loop
        ptr = strtok(&in->s[9], ";");
        if (ptr == NULL)
        {
            printf("badly formatted forprime loop: forprime(var=start; stop; body)\n");
            exit(3);
        }
        sprintf(vname, "%s_start", pre);
        if (set_strvar(vname, ptr))
            new_strvar(vname, ptr);
        ptr = strtok(NULL, ";");
        if (ptr == NULL)
        {
            printf("badly formatted forprime loop: forprime(var=start; stop; body)\n");
            exit(3);
        }
        sprintf(vname, "%s_stop", pre);
        if (set_strvar(vname, ptr))
            new_strvar(vname, ptr);
        ptr = strtok(NULL, "\0");
        if (ptr == NULL)
        {
            printf("badly formatted forprime loop: forprime(var=start; stop; body)\n");
            exit(3);
        }
        strncpy(str, ptr, strlen(ptr) - 1);
        str[strlen(ptr) - 1] = '\0';
        sprintf(vname, "%s_body", pre);
        if (set_strvar(vname, str))
            new_strvar(vname, str);

        sprintf(in->s, "forprime(%s_start, %s_stop, %s_body);", pre, pre, pre);
    }

    if ((ptr = strstr(in->s, "if(")) != NULL)
    {
        // new if statement
        char pre[8], vname[20];
        sprintf(pre, "if%d", if_cnt++);

        // tokenize the branch
        char *eptr;
        ptr = strtok(&in->s[3], ";");
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
        else if (eptr[strlen(eptr)+1] == '\0')
        {
            // no else statement and no output suppression character
            strncpy(str, eptr, strlen(eptr) - 1);
            str[strlen(eptr) - 1] = '\0';
            sprintf(vname, "%s_body", pre);
            if (set_strvar(vname, str))
                new_strvar(vname, str);

            sprintf(in->s, "if(%s_cond, %s_body);", pre, pre);
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
                sprintf(in->s, "if(%s_cond, %s_body, %s_elsebody);", pre, pre, pre);
            }
            else if (eptr[strlen(eptr) + 1] == ')')
            {
                // just the if, with an output suppression character
                sprintf(str, "%s;", eptr);
                sprintf(vname, "%s_body", pre);
                if (set_strvar(vname, str))
                    new_strvar(vname, str);

                sprintf(in->s, "if(%s_cond, %s_body);", pre, pre);
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
                sprintf(in->s, "if(%s_cond, %s_body, %s_elsebody);", pre, pre, pre);
            }
        }

    }

	return 0;
}

int is_new_token(int el_type, int el_type2)
{

	if (el_type == EOE || el_type == LP || el_type == RP)
		return 1;

	if (el_type != el_type2)
	{
		//types are different
		if (el_type == CH && el_type2 == NUM)
		{
			//but this could be a function or variable name
			//so not different
			return 0;
		}
		else if (el_type == NUM && el_type2 == CH)
		{
			//but this could be a function or variable name
			//so not different
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

	//  a token in this context is one of the following things:
	//    a number, possibly including a base prefix (0x, 0d, 0b, etc)
	//    a variable name
	//    a function name
	//    an operator string (includes parens, commas)

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
				//printf("growing tmpsize in tokenize...\n");
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

void get_expression(char *in, str_t *out)
{
	char *tok;
	char delim[] = {' ', '\0'};
	int invalid_string = 0;
	bstack_t stk;
	str_t *tmp,*tmp1,*tmp2;

	// if there is no input to process, we are done... 
	if (in == NULL) return;

	stack_init(20,&stk,STACK);
	tmp = (str_t *)malloc(sizeof(str_t));
	sInit(tmp);
	tmp1 = (str_t *)malloc(sizeof(str_t));
	sInit(tmp1);
	tmp2 = (str_t *)malloc(sizeof(str_t));
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

int process_expression(char *input_exp, fact_obj_t *fobj, int force_quiet)
{
	str_t str;
	char *ptr;

    sInit(&str);
    toStr(input_exp, &str);
	preprocess(&str);

	// new factorization
	//fobj->refactor_depth = 0;
    calc_with_assignment(&str, fobj, force_quiet);

	sFree(&str);
	return 0;
}

void calc_with_assignment(str_t *in, fact_obj_t *fobj, int force_quiet)
{
    char *ptr;
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
        strcpy(varname, "ans");

    // look for a trailing semicolon
    if (in->s[strlen(in->s) - 1] == ';')
    {
        nooutput = 1;
        in->s[strlen(in->s) - 1] = '\0';
    }
    else
        nooutput = 0;

    toStr(in->s + offset, &str);
    if (!calc(&str, fobj))
    {
        if (strcmp(str.s, "") != 0)
        {
            mpz_set_str(tmp, str.s, 0);
            sCopy(&str, in);

            // always set the default variable to the new answer
            set_uvar("ans", tmp);
            // and optionally any assigned variable as well.
            if (set_uvar(varname, tmp))
                new_uvar(varname, tmp);

            if ((nooutput == 0) && (force_quiet == 0))
            {
                if (OBASE == DEC)
                {
                    if (VFLAG > 0)
                        gmp_printf("\n%s = %Zd\n\n", varname, tmp);
                    else
                        gmp_printf("%Zd\n", tmp);
                }
                else if (OBASE == HEX)
                {
                    if (VFLAG > 0)
                        gmp_printf("\n%s = %Zx\n\n", varname, tmp);
                    else
                        gmp_printf("%Zx\n", tmp);
                }
            }
        }
    }

    mpz_clear(tmp);
    sFree(&str);
    return;
}

int calc(str_t *in, fact_obj_t *fobj)
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

	stack_init(20,&stk,STACK);
	tmp = (str_t *)malloc(sizeof(str_t));
	sInit(tmp);	
	mpz_init(tmpz);
	post = (str_t *)malloc(sizeof(str_t));
	sInit(post);

	//run the shunting algorithm
	i=0;
	post->s[0] = '\0';
	while (i<num_tokens)
	{
		switch (token_types[i])
		{
		case 3:
			//NUM
			//to num output queue
			sAppend(" ",post);
			sAppend(tokens[i],post);
			break;
		case 6:
			//LP
			toStr(tokens[i],tmp);
			push(tmp,&stk);
			break;
		case 7:
			//string (function or variable name)
			varstate = get_uvar(tokens[i],tmpz);
			if (varstate == 0)
			{
				//found a variable with that name, copy it's value
				//to num output queue
				sAppend(" ",post);
				sAppend(tokens[i],post);
			}
			else if (varstate == 2)
			{
				//do nothing, special case
			}
			else if (getFunc(tokens[i],&na) >= 0) 
			{
				//valid function, push it onto stack
				toStr(tokens[i],tmp);
				push(tmp,&stk);
			}
			else
			{
				// a non-numeric string that is not a variable or a function.
                // the only thing that makes sense is for it to be a string
                // representing a new assignment (new variable).  Assume that's the
                // case.  if not, errors will be raised further down the line.
				//printf("unrecognized token: %s\n",tokens[i]);
				//retval=1;
				//goto free;
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
					//stack empty and we are still looking for a LP
					printf("bad function separator position or mismatched parens\n");
					retval = 1;
					goto free;
				}

				if (strcmp(tmp->s,"(") == 0)
				{
					//found a left paren.  put it back and continue
					push(tmp,&stk);
					break;
				}
				else
				{
					//copy to output operator queue
					sAppend(" ",post);
					sAppend(tmp->s,post);
				}
			}
			break;
		case 5:
			//right paren
			while (1)
			{
				if (pop(tmp,&stk) == 0)
				{
					//stack empty and we are still looking for a LP
					printf("mismatched parens\n");
					retval = 1;
					goto free;
				}

				if (strcmp(tmp->s,"(") == 0)
				{
					//found a left paren.  ignore it.
					if (pop(tmp,&stk) != 0)
					{
						//is the top of stack a function?
						if ((getFunc(tmp->s,&na) >= 0) && (strlen(tmp->s) > 1))
						{
							//the extra check for strlen > 1 fixes
							//the case where the string is an operator, not
							//a function.  for multichar operators this won't work
							//instead should probably separate out the operators
							//from the functions in getFunc


							//yes, put it on the output queue as well
							sAppend(" ",post);
							sAppend(tmp->s,post);
						}
						else
						{
							//no, put it back
							push(tmp,&stk);
						}
					}
					break;
				}
				else
				{
					//copy to output queue
					sAppend(" ",post);
					sAppend(tmp->s,post);
				}
			}
			break;
		case 4:
			//operator
			while (pop(tmp,&stk))
			{
				if (strlen(tmp->s) == 1 && getOP(tmp->s[0]) > 0)
				{
					//its an operator
					//check the precedence
					if (op_precedence(tmp->s,tokens[i],getAssoc(tmp->s)))
					{
						//push to output op queue
						sAppend(" ",post);
						sAppend(tmp->s,post);
					}
					else
					{
						//put the tmp one back and bail
						push(tmp,&stk);
						break;
					}
				}
				else
				{
					//its not an operator, put it back and bail
					push(tmp,&stk);
					break;
				}
			}
			//push the current op onto the stack.
			toStr(tokens[i],tmp);
			push(tmp,&stk);

			break;
		case 2:
			//post unary operator
			//I think we can jush push these into the output operator queue
			toStr(tokens[i],tmp);
			sAppend(" ",post);
			sAppend(tmp->s,post);
			break;
		}
		i++;
	}

	//now pop all operations left on the stack to the output queue
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

	//free the input tokens
	for (i=0;i<num_tokens;i++)
		free(tokens[i]);
	free(tokens);

	// process the output postfix expression:
	// this can be done with a simple stack
	// all tokens are separated by spaces
	// all tokens consist of numbers or functions

	/* try to grab the input number as an expression, if found, store it in fobj->str_N */
	strN = strdup(post->s); /* make a copy of the postfix string */
	if (strN != NULL)
	{
		/* find input expression and store in str_N */
		get_expression(strN, &(fobj->str_N));
        if (CALC_VERBOSE)
		    printf("Found expression: %s\n", fobj->str_N.s);
		free(strN);
	}

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
		// printf("nothing to process\n");
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
			//printf("pushing %s onto stack\n",tok);
			toStr(tok,tmp);
			push(tmp,&stk);
			break;
		case AMBIG:
			// could be a number or a function
			// if the next character is a number, this is too
			// else its an operator
			//printf("peeking at tok + 1: %d\n",(int)tok[1]);
			if (get_el_type2(tok[1]) == NUM)
			{
				//printf("pushing %s onto stack\n",tok);
				toStr(tok,tmp);
				push(tmp,&stk);
				break;
			}

			// if not a num, proceed into the next switch (function handle)
		default:
			func = getFunc(tok,&na);
            if (CALC_VERBOSE)
			    printf("processing function %d\n",func);

			if (func >= 0)
			{
				//pop those args and put them in a global array
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
					//printf("looking at argument %s\n",tmp->s);
					r = mpz_set_str(tmpz, tmp->s, 0);
                    if (r < 0)
                    {
                        //printf("input is not a valid number in a discoverable base\n");
                        //printf("placing %s into character operand array\n", tmp->s);
                        strcpy(choperands[na - j - 1], tmp->s);
                    }
                    else
                    {
                        // it is a number, put it in the operand pile
                        //gmp_printf("found numerical argument %Zd\n", tmpz);
                        mpz_set(operands[na - j - 1], tmpz);
                    }
				}

				na = j;
				// call the function evaluator with the 
				// operator string and the number of args available
				na = feval(func,na,fobj);

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
			else
			{
				int sz;

                // the string token is not a function, check if it's a defined variable
                if (get_uvar(tok, tmpz) == 0)
                {
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
                else
                {
                    // not a defined variable... push the string onto the stack
                    // where we will later treat it as a character operand to a function.
                    //printf("pushing %s onto stack\n", tok);
                    toStr(tok, tmp);
                    push(tmp, &stk);
                }
			}
			break;
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
	//return the opcode associated with the function, and
	//the number of arguments it takes
	int i,j;

	char func[NUM_FUNC][11] = {"fib","luc","snfs","expr","rsa",
						"gcd","jacobi","factor","rand","lg2",
						"log","ln","pm1","pp1","rho",
						"trial","mpqs","nextprime","size","issquare",
						"isprime","squfof","sqrt","modinv","modexp",
						"nroot","shift","siqs","primes", "ispow", 
						"torture","randb", "ecm","+","-",
						"*","/","!","#","eq",
						"<<",">>","%","^","redc",
						"puzzle","sieve","algebraic","llt","siqsbench",
						"pullp","sftest","smallmpqs","testrange","siqstune",
						"ptable","sieverange","fermat","nfs","tune",
						"xor", "and", "or", "not", "frange",
						"bpsw","aprcl","lte", "gte", "<", 
						">","dummy","if","print", "for",
                        "forprime", "dummy", "dummy", "dummy", "dummy",
                        "dummy", "dummy", "dummy", "dummy", "dummy",
                        "dummy", "dummy", "dummy", "dummy", "dummy",
                        "dummy", "dummy", "dummy", "dummy", "dummy",
                        "dummy", "dummy", "dummy", "dummy", "dummy"};

	int args[NUM_FUNC] = {1,1,2,1,1,
					2,2,1,1,1,
					1,1,1,2,1,
					2,1,2,1,1,
					1,1,1,2,3,
					2,2,1,3,1, 
					2,1,2,2,2,
					2,2,1,1,2,
					2,2,2,2,2,
					2,2,2,1,0,
					0,5,1,4,1,
					0,4,3,1,0,
					2,2,2,1,2,
					1,1,2,2,2,
					2,2,3,1,4,
                    3,-1,-1,-1,-1,
                    -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1};

	for (i=0;i<NUM_FUNC;i++)
	{
		j = strcmp(func[i],s);
		if (j == 0)
		{
			*nargs = args[i];
			return i;
		}
	}

	return -1;
}

int feval(int func, int nargs, fact_obj_t *fobj)
{
	// evaluate the function 'func', with 'nargs' argument(s) located
	// in the mpz_t array 'operands'.
	// place return values in operands[0]
	mpz_t mp1, mp2, mp3, tmp1, tmp2;
	mpz_t gmpz;

	str_t str;
	uint32 i=0;
	uint64 n64;
	uint32 j,k;
	//uint32 *offsets;
	double t;
	struct timeval tstart, tstop;
	uint64 lower, upper, inc, count;

	mpz_init(mp1);
	mpz_init(mp2);
	mpz_init(mp3);
	mpz_init(tmp1);
	mpz_init(tmp2);
	sInit(&str);

	switch (func)
	{
	case 0:
		//fib - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in fib\n");
			break;
		}
		mpz_fib_ui(operands[0], mpz_get_ui(operands[0]));
		break;
	case 1:
		//luc - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in luc\n");
			break;
		}
		mpz_lucnum_ui(operands[0], mpz_get_ui(operands[0]));
		break;
		
	case 2:
		// snfs - two arguments
		if (nargs != 2)
		{
			printf("wrong number of arguments in snfs\n");
			break;
		}

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
		print_factors(fobj);

		break;
	case 3:
		// expr - one argument
		// this is used to evaluate numerical expressions from the command line,
		// now that the default behavior is to factor the input.
		// since whatever was in the function will have already been evaluated,
		// simply return the input.
		if (nargs != 1)
		{
			printf("wrong number of arguments in expr\n");
			break;
		}

		break;
	case 4:
		//rsa - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in gcd\n");
			break;
		}
		mpz_set_ui(mp1, 2048);
		mpz_set_ui(mp2, 4096);
		if (mpz_cmp(operands[0], mp2) > 0)
		{
			printf("bitlength too large");
			mpz_set_ui(operands[0], 1);
			break;
		}
		else if (mpz_cmp(operands[0],mp1) > 0)
		{
			printf("Paranoid, huh?  This might take a minute\n");
		}

		build_RSA(mpz_get_ui(operands[0]), operands[0]);
		break;
	case 5:
		//gcd - two arguments
		if (nargs != 2)
		{
			printf("wrong number of arguments in gcd\n");
			break;
		}
			
		mpz_gcd(operands[0], operands[0], operands[1]);

		break;
	case 6:
		//jacobi - two arguments
		if (nargs != 2)
		{
			printf("wrong number of arguments in jacobi\n");
			break;
		}

		if (mpz_odd_p(operands[1]))
			mpz_set_si(operands[0], mpz_jacobi(operands[0], operands[1]));
		else
			printf("jacobi defined only for odd b!\n");

		break;
	case 7:
		//factor - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in factor\n");
			break;
		}

		mpz_set(fobj->N, operands[0]);
		factor(fobj);
		mpz_set(operands[0], fobj->N);
		print_factors(fobj);
		break;
	case 8:
		//rand - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in rand\n");
			break;
		}
		mpz_set_ui(operands[1], 10);
		mpz_pow_ui(operands[1], operands[1], mpz_get_ui(operands[0]));
		mpz_urandomm(operands[0], gmp_randstate, operands[1]);
		break;
	case 9:
		//lg2 - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in lg2\n");
			break;
		}
		mpz_set_ui(operands[0], mpz_sizeinbase(operands[0], 2));
		break;
	case 10:
		//log - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in log\n");
			break;
		}
		mpz_set_ui(operands[0], mpz_sizeinbase(operands[0], 10));
		break;
	case 11:
		//ln - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in ln\n");
			break;
		}
		mpz_set_ui(operands[0], (uint32)((mpz_sizeinbase(operands[0], 2)-1) * log(2.0)));
		break;
	case 12:
		//pm1 - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in pm1\n");
			break;
		}
		mpz_set(fobj->N, operands[0]);
		mpz_set(fobj->pm1_obj.gmp_n, operands[0]);
		pollard_loop(fobj);
		mpz_set(operands[0], fobj->pm1_obj.gmp_n);
		print_factors(fobj);
		break;
	case 13:
		//pp1 - two arguments
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
		print_factors(fobj);
		break;
	case 14:
		//rho - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in rho\n");
			break;
		}
		mpz_set(fobj->N, operands[0]);
		mpz_set(fobj->rho_obj.gmp_n, operands[0]);
		brent_loop(fobj);
		mpz_set(operands[0], fobj->rho_obj.gmp_n);
		print_factors(fobj);
		break;
	case 15:
		//trial - two arguments
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
		
		print_factors(fobj);
		break;
	case 16:
		//mpqs - one argument
		printf("no longer supported\n");
		break;

	case 17:
		//next prime - two arguments
		if (nargs == 2)
		{
			zNextPrime(operands[0], operands[2], mpz_get_ui(operands[1]));
			mpz_set(operands[0], operands[2]);
		}
		else if (nargs == 1)
		{
			zNextPrime(operands[1], operands[2], 1);
			mpz_set(operands[0], operands[2]);
		}
		else
		{
			printf("wrong number of arguments in NextPrime\n");
			break;
		}
		
		break;
	case 18:
		//size - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in size\n");
			break;
		}

		printf("%d digits, %d bits\n", gmp_base10(operands[0]),
			(int)mpz_sizeinbase(operands[0], 2));

		mpz_set_ui(operands[0], mpz_sizeinbase(operands[0], 2));

		break;
	case 19:
		//issquare
		if (nargs != 1)
		{
			printf("wrong number of arguments in issquare\n");
			break;
		}

		i = mpz_perfect_square_p(operands[0]);

		if (i)
			printf("input is square\n");
		else
			printf("input is not square\n");

		mpz_set_ui(operands[0], i);

		break;
	case 20:
		//isprime - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in isprime\n");
			break;
		}

		i = is_mpz_prp(operands[0]);

        if (VFLAG > 2)
        {
            if (i)
                printf("probably prime\n");
            else
                printf("not prime\n");
        }

		mpz_set_ui(operands[0], i);

		break;
	case 21:
		//shanks - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in squfof\n");
			break;
		}
		mpz_init(gmpz);

		n64 = sp_shanks_loop(operands[0],fobj);
		print_factors(fobj);
		mpz_set_64(operands[0], n64);
		break;
	case 22:
		//sqrt - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in sqrt\n");
			break;
		}
		mpz_root(operands[0], operands[0], 2);
		break;
	case 23:
		//modinv - two arguments
		if (nargs != 2)
		{
			printf("wrong number of arguments in modinv\n");
			break;
		}
        printf("modinv_1 = %u\n", modinv_1(operands[0]->_mp_d[0], operands[1]->_mp_d[0]));
		mpz_invert(operands[0], operands[0], operands[1]);        
		break;
	case 24:
		//modexp - three arguments
		if (nargs != 3)
		{
			printf("wrong number of arguments in modexp\n");
			break;
		}
		mpz_powm(operands[0], operands[0], operands[1], operands[2]);
		break;
	case 25:
		//nroot - two arguments
		if (nargs != 2)
		{
			printf("wrong number of arguments in Nroot\n");
			break;
		}
		mpz_root(operands[0], operands[0], mpz_get_ui(operands[1]));
		break;
	case 26:
		//shift - two arguments
		if (nargs != 2)
		{
			printf("wrong number of arguments in shift\n");
			break;
		}

		if (mpz_sgn(operands[1]) >= 0)
			mpz_mul_2exp(operands[0], operands[0], mpz_get_ui(operands[1]));
		else
			mpz_tdiv_q_2exp(operands[0], operands[0], -1*mpz_get_si(operands[1]));

		break;
	case 27:
		//siqs - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in siqs\n");
			break;
		}
		mpz_set(fobj->N, operands[0]);
		mpz_set(fobj->qs_obj.gmp_n, operands[0]);
		SIQS(fobj);
		mpz_set(operands[0], fobj->qs_obj.gmp_n);
		print_factors(fobj);
		break;

	case 28:
		//primes
		if (nargs == 2)
		{
			gettimeofday(&tstart, NULL);
			lower = mpz_get_64(operands[1]);
			upper = mpz_get_64(operands[2]);

			soe_wrapper(spSOEprimes, szSOEp, lower, upper, 1, &n64);

			mpz_set_64(operands[0], n64);
			gettimeofday (&tstop, NULL);
			t = my_difftime (&tstart, &tstop);

			printf("elapsed time = %6.4f\n",t);
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
		free(PRIMES);
		PRIMES = soe_wrapper(spSOEprimes, szSOEp, lower, upper, mpz_get_ui(operands[2]), &NUM_P);
		if (PRIMES != NULL)
		{
			P_MIN = PRIMES[0];
			P_MAX = PRIMES[NUM_P-1];
		}
		mpz_set_64(operands[0], NUM_P);

		break;
	case 29:
		//ispow - one argument

		i = mpz_perfect_power_p(operands[0]);

		if (i)
			printf("input is a perfect power\n");
		else
			printf("input is not a perfect power\n");

		mpz_set_ui(operands[0], i);

		break;

	case 30:
		//torture - two arguments
		if (nargs != 2)
		{
			printf("wrong number of arguments in torture\n");
			break;
		}
		i = mpz_get_ui(operands[0]);
		k = mpz_get_ui(operands[1]);
		for (j=0; j<i; j++)
		{
			printf("***********RUN %d***********\n",j+1);
			mpz_urandomb(operands[2], gmp_randstate, k);
			mpz_set(fobj->N, operands[2]);
			factor(fobj);
			mpz_set(operands[0],fobj->N);
			print_factors(fobj);
			clear_factor_list(fobj);
		}

		break;
	case 31:
		//randb - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in randb\n");
			break;
		}
		mpz_urandomb(operands[0], gmp_randstate, mpz_get_ui(operands[0]));
		break;
	case 32:
		//ecm - two arguments
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
		
		print_factors(fobj);
		break;

	case 33:
		//add
		if (nargs != 2)
		{
			printf("wrong number of arguments in add\n");
			break;
		}
		mpz_add(operands[0], operands[0], operands[1]);
		break;

	case 34:
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

	case 35:
		//mul
		if (nargs != 2)
		{
			printf("wrong number of arguments in mul\n");
			break;
		}
		mpz_mul(operands[0], operands[0], operands[1]);
		break;

	case 36:
		//div
		if (nargs != 2)
		{
			printf("wrong number of arguments in div\n");
			break;
		}
		mpz_tdiv_q(operands[0], operands[0], operands[1]);
		break;

	case 37:
		//!
		if (nargs != 1)
		{
			printf("wrong number of arguments in factorial\n");
			break;
		}
		mpz_fac_ui(operands[0], mpz_get_ui(operands[0]));
		break;

	case 38:
		//primorial
		if (nargs != 1)
		{
			printf("wrong number of arguments in primorial\n");
			break;
		}
		else
		{
			mpz_primorial_ui(operands[0], mpz_get_ui(operands[0]));
		}
		break;

	case 39:
		// eq
		if (nargs != 2)
		{
			printf("wrong number of arguments in eq\n");
			break;
		}
		else
		{
			mpz_set_ui(operands[0], mpz_cmp(operands[0], operands[1]) == 0);
		}

		break;

	case 40:
		//<<
		if (nargs != 2)
		{
			printf("wrong number of arguments in <<\n");
			break;
		}
		mpz_mul_2exp(operands[0], operands[0], mpz_get_ui(operands[1]));
		break;

	case 41:
		//>>
		if (nargs != 2)
		{
			printf("wrong number of arguments in >>\n");
			break;
		}
		mpz_tdiv_q_2exp(operands[0], operands[0], mpz_get_ui(operands[1]));
		break;

	case 42:
		//mod
		if (nargs != 2)
		{
			printf("wrong number of arguments in mod\n");
			break;
		}
		mpz_mod(operands[0], operands[0], operands[1]);
		break;

	case 43:
		//exp
		if (nargs != 2)
		{
			printf("wrong number of arguments in exp\n");
			break;
		}
		mpz_pow_ui(operands[0], operands[0], mpz_get_ui(operands[1]));
		break;

	case 44:
        // REDC (redc)
        if (nargs != 2)
        {
            printf("wrong number of arguments in REDC\n");
            break;
        }

        {
            uint32 b, x, mhat;
            uint32 a[2 * 256 + 1];
            uint32 m[2 * 256 + 1];
            uint32 c[2 * 256 + 1], *_c, *tmpm, mu;
            int      oldused, xx, y, pa, j;
            char *strn;            

            // REDC
            // export the input as a 4-byte array, least significant word first,
            // native endianness, full words (0 nails bits)
            mpz_export(a, &oldused, -1, 4, 0, 0, operands[0]);
            mpz_export(m, &pa, -1, 4, 0, 0, operands[1]);

            printf("size of a = %d, size of m = %d\n", oldused, pa);

            b = m[0];

            // this first equation results from taking our starting point a^-1 = 1 mod 2
            // and lifting it twice to produce -a^3 + 4a^2 - 6a + 4 mod 16.
            // if you work out the above algebraicly using the bit representation of a mod 16,
            // you get this result.
            x = (((b + 2) & 4) << 1) + b; // here x*a==1 mod 2**4
            x *= 2 - b * x;               // here x*a==1 mod 2**8
            x *= 2 - b * x;               // here x*a==1 mod 2**16
            x *= 2 - b * x;               // here x*a==1 mod 2**32

            // mhat = -1/m mod b
            mhat = (uint32)(((uint64)1 << ((uint64)32)) - ((uint64)x));

            printf("mhat = %u\n", mhat);

            /* copy the input */
            for (xx = 0; xx < oldused; xx++)
            {
                c[xx] = a[xx];
            }
            for (; xx < 2 * pa + 1; xx++)
            {
                c[xx] = 0;
            }

            for (xx = 0; xx < pa; xx++)
            {
                uint32 cy = 0;

                // get Mu for this round
                mu = c[xx] * mhat;

                if (VFLAG > 0)
                    printf("round %d mu = %08x\n", xx, mu);

                // compute c += mu*n*b^x, where b=2^DIGITBITS
                _c = c + xx;          // b^x
                tmpm = m;      // n
                y = 0;

#ifdef __unix__
                for (; y < pa; y++)
                {
                    __asm__(\
                        "movl %5,%%eax \n\t"
                        "mull %4       \n\t"
                        "addl %1,%%eax \n\t"
                        "adcl $0,%%edx \n\t"
                        "addl %%eax,%0 \n\t"
                        "adcl $0,%%edx \n\t"
                        "movl %%edx,%1 \n\t"
                        :"=g"(_c[0]), "=r"(cy)
                        : "0"(_c[0]), "1"(cy), "r"(mu), "r"(*tmpm++)
                        : "%eax", "%edx", "memory", "%cc");

                    ++_c;
                }
#endif

                while (cy)
                {
                    do { uint32 t = _c[0] += cy; cy = (t < cy); } while (0);
                    ++_c;
                }

                if (VFLAG > 0)
                {
                    printf("round %d partial sum: \n", xx);
                    for (y = 2*pa - 1; y >= 0; y--)
                    {
                        printf("%08x", _c[y]);
                    }
                    printf("\n");
                }

            }

            /* now copy out */
            _c = c + pa;
            tmpm = a;
            for (xx = 0; xx < pa + 1; xx++)
            {
                *tmpm++ = *_c++;
            }

            for (; xx < oldused; xx++)
            {
                *tmpm++ = 0;
            }

            mpz_import(operands[0], pa + 1, -1, 4, 0, 0, a);

            if (mpz_cmp(operands[0], operands[1]) >= 0)
            {
                mpz_sub(operands[0], operands[0], operands[1]);
            }


        }

		break;

	case 45:
        // modmul

		break;

	case 46:
		//sieve
		//guess at the number of expected results
		break;	// not supported in official release

	case 47:
		//algebraic
		break;	// not supported in official release

	case 48:
		//lucas lehmer test
		if (nargs != 1)
		{
			printf("wrong number of arguments in llt\n");
			break;
		}

        if (llt(mpz_get_ui(operands[0])))
        {
            if (VFLAG > 1)
                gmp_printf("%Zd is prime!\n", operands[0]);
            mpz_set_ui(operands[0], 1);
        }
        else
        {
            if (VFLAG > 1)
                gmp_printf("%Zd is composite.\n", operands[0]);
            mpz_set_ui(operands[0], 0);
        }
		break;

	case 49:
		//siqsbench
		if (nargs != 0)
		{
			printf("wrong number of arguments in siqsbench\n");
			break;
		}
		siqsbench(fobj);
		break;
	case 50:
		

		break;
	case 51: 

		break;	// not supported in official release
		// small factorization test routine
		// 5 inputs

		//if (0)
		//{
		//	FILE *test;
		//	test = fopen("res.dat", "r");
		//	char str[100];
		//	fp_digit *n,r;

		//	n = (fp_digit *)malloc(500000 * sizeof(fp_digit));
		//	while (!feof(test))
		//	{
		//		if (fgets(str, 100, test) == NULL)
		//			break;

		//		sscanf(str, "%llu", &n[i++]);				
		//	}

		//	gettimeofday(&tstart, NULL);
		//	for (i=0; i<500000; i++)
		//		spModExp(2, n[i] - 1, n[i], &r);
		//	fclose(test);
		//	gettimeofday (&tstop, NULL);
		//	difference = my_difftime (&tstart, &tstop);
		//	t = ((double)difference->secs + (double)difference->usecs / 1000000);
		//	free(difference);
		//	printf("%1.4f sec to test 500000 inputs from a c85\n", t);

		//	free(n);
		//}

		{
			mpz_t *input;
			int test_squfof = 0;
			int test_lehman = 0;
			int test_tinysiqs = 0;
			int test_squfof_tf = 0;
			//int force_hard = 1;
			int do_tf = 0;
			int use_hartolf = 0;
			double lehTune = 2.5;
			int numin = 100000;
			uint64 *outputs;
			int bits = 115;
			double sum;
			int maxsz, minsz;
			fact_obj_t *fobj2;
			mpz_t gmptmp;
			int tmpv;
			int startbits = 90, stopbits = 120, stepbits = 5;
			//int startbits = 40, stopbits = 60, stepbits = 5;

			// get numin
			numin = mpz_get_ui(operands[0]);

			// get method
			if (mpz_get_ui(operands[1]) == 1)
			{ 
				printf("commencing %d tests of squfof method\n",numin);
				test_squfof = 1;
				test_lehman = 0;
				test_tinysiqs = 0;
				test_squfof_tf = 0;
			}
			else if (mpz_get_ui(operands[1]) == 2)
			{ 
				printf("commencing %d tests of lehman method\n",numin);
				test_squfof = 0;
				test_lehman = 1;
				test_tinysiqs = 0;
				test_squfof_tf = 0;
				//init_lehman();
			}
			else if (mpz_get_ui(operands[1]) == 3)
			{ 
				printf("commencing %d tests of squfof + tf method\n",numin);
				test_squfof = 0;
				test_lehman = 0;
				test_squfof_tf = 1;
				test_tinysiqs = 0;
				//init_lehman();
			}
			else if (mpz_get_ui(operands[1]) == 4)
			{ 
				printf("commencing %d tests of tinySIQS method\n",numin);
				test_squfof = 0;
				test_lehman = 0;
				test_squfof_tf = 0;
				test_tinysiqs = 1;
				//init_lehman();
			}
			else
			{ 
				printf("commencing %d tests of smallmpqs method\n",numin);
				test_squfof = 0;
				test_lehman = 0;
				test_squfof_tf = 0;
				test_tinysiqs = 0;
			}

			// get other parameters
			if (nargs >= 3)
				do_tf = mpz_get_ui(operands[2]);

			if (nargs >= 4)
				lehTune = mpz_get_ui(operands[3]);
			
			if (nargs == 5)
				use_hartolf = mpz_get_ui(operands[4]);

			printf("sweeping inputs with bit sizes from %d to %d in steps of %d\n",
				startbits, stopbits, stepbits);

			if (test_lehman)
				printf("lehman parameters:\n  do TF = %d\n  LehTune = %1.0f\n  Use Hart_OLF = %d\n",
				do_tf, lehTune, use_hartolf);

			fobj2 = (fact_obj_t *)malloc(sizeof(fact_obj_t));
			init_factobj(fobj2);
			input = (mpz_t *)malloc(numin * sizeof(mpz_t));
			for (i=0; i<numin; i++)
				mpz_init(input[i]);

			outputs = (uint64 *)malloc(numin * sizeof(uint64));

			mpz_init(gmptmp);

			for (bits=startbits; bits<=stopbits; bits+=stepbits)
			{				
				int nbits;
				int proper_factors;

				gettimeofday(&tstart, NULL);

				sum = 0;
				maxsz = 0;
				minsz = 999;
				for (i=0; i<numin; i++)
				{
					fp_digit b;

					b = bits-2;
					mpz_urandomb(mp1, gmp_randstate, b/2); //zRandb(&mp1, b/2);
					//if (force_hard) 
					//	mp1.val[(b/2)/BITS_PER_DIGIT] |= (fp_digit)(1 << ((b/2) % BITS_PER_DIGIT));
					mpz_urandomb(mp2, gmp_randstate, b/2); //zRandb(&mp2, b/2);
					//if (force_hard) 
					//	mp2.val[(b/2)/BITS_PER_DIGIT] |= (fp_digit)(1 << ((b/2) % BITS_PER_DIGIT));
					
					mpz_nextprime(mp1, mp1);
					mpz_nextprime(mp2, mp2);

					mpz_mul(input[i], mp1,mp2);
					nbits = mpz_sizeinbase(input[i],2);				
					sum += nbits;
					if (nbits < minsz)
						minsz = nbits;
					if (nbits > maxsz)
						maxsz = nbits;
				}

				gettimeofday (&tstop, NULL);
				t = my_difftime (&tstart, &tstop);

				if (VFLAG > 0)
				{
					printf("generation of %d %d bit inputs took = %6.4f\n",numin,bits,t);
					printf("average size was %6.4f bits (max = %d, min = %d)\n",
						sum/(double)numin, maxsz, minsz);
				}
		
				gettimeofday(&tstart, NULL);
				for (i=0; i<numin; i++)
				{
					if ((i % (numin/10) == 0) && (VFLAG > 0))
						printf("input %d\n",i); //, %d correct\n",i,correct);
				
					mpz_set(fobj2->qs_obj.gmp_n, input[i]);
					if (test_squfof)
					{
						mpz_set(gmptmp, fobj2->qs_obj.gmp_n);
						outputs[i] = sp_shanks_loop(gmptmp, fobj);
					}
					else if (test_lehman)
					{
						uint64 N = mpz_get_64(input[i]);
						outputs[i] = LehmanFactor(N, lehTune, do_tf, 0.1);
					}
					else if (test_squfof_tf)
					{
						int k = 1;
						mpz_set(gmptmp, fobj2->qs_obj.gmp_n);

						// tf
						outputs[i] = 0;
						while (PRIMES[k] < 50)
						{
							uint64 q = (uint64)PRIMES[k];
							uint64 r = mpz_tdiv_ui(gmptmp, q);
		
							if (r != 0)
								k++;
							else	
							{
								outputs[i] = q;
								break;
							}
						}

						if (outputs[i] == 0)
							outputs[i] = sp_shanks_loop(gmptmp, fobj);

						if (outputs[i] == 0)
						{
							while (PRIMES[k] < 65536)
							{
								uint64 q = (uint64)PRIMES[k];
								uint64 r = mpz_tdiv_ui(gmptmp, q);
		
								if (r != 0)
									k++;
								else	
								{
									outputs[i] = q;
									break;
								}
							}
						}

					}
					else if (test_tinysiqs)
					{
						fobj2->qs_obj.flags = 12345;
						tmpv = VFLAG;
						VFLAG = -1;
						//tinySIQS(fobj2->qs_obj.gmp_n, NULL, NULL);

						//if (mpz_cmp_ui(fobj2->qs_obj.gmp_n, 1) != 0)
						//{						
						//	gmp_printf("not fully factored!  residue = %Zd\n", fobj2->qs_obj.gmp_n);
						//	print_factors(fobj2);
						//}
						VFLAG = tmpv;
						//clear_factor_list(fobj2);
					}
					else
					{
						fobj2->qs_obj.flags = 12345;
						smallmpqs(fobj2); //SIQS(fobj2);

						if (mpz_cmp_ui(fobj2->qs_obj.gmp_n, 1) != 0)
						{						
							gmp_printf("not fully factored!  residue = %Zd\n", fobj2->qs_obj.gmp_n);
							print_factors(fobj2);
						}

						clear_factor_list(fobj2);
					}
				}
			
			
				gettimeofday (&tstop, NULL);
				t = my_difftime (&tstart, &tstop);

				if (test_lehman || test_squfof || test_squfof_tf)
				{
					proper_factors = 0;
					for (i=0; i<numin; i++)
					{
						if ((mpz_get_64(input[i]) % outputs[i] == 0) && 
							(outputs[i] > 1) &&
							(mpz_get_64(input[i]) != outputs[i]))
							proper_factors++;
					}
					if (VFLAG > 0)
						printf("factoring %d %d bit inputs took = %6.4f, found %d proper factors\n",
							numin,bits,t,proper_factors);
					else
						printf("%6.2f %6.4f %d\n", sum / (double)numin, t, proper_factors);
				}
				else
				{
					if (VFLAG > 0)
						printf("factoring %d %d bit inputs took = %6.4f\n",numin,bits,t);
					else
						printf("%6.2f %6.4f\n", sum / (double)numin, t);
				}
					
			}

			// free temps
			free_factobj(fobj2);
			free(fobj2);

			for (i=0; i<numin; i++)
				mpz_clear(input[i]);

			free(input);
			free(outputs);
			mpz_clear(gmptmp);
		}

		break;

	case 52:
		//smallmpqs - 1 argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in smallmpqs\n");
			break;
		}
		mpz_set(fobj->N, operands[0]);
		mpz_set(fobj->qs_obj.gmp_n, operands[0]);
		fobj->logfile = fopen(fobj->flogname,"a");
		if (fobj->logfile == NULL)
			printf("fopen error: %s\n", strerror(errno));
		smallmpqs(fobj);		
		//tinySIQS(fobj->qs_obj.gmp_n, NULL, NULL);
		if (fobj->logfile != NULL)
			fclose(fobj->logfile);
		mpz_set(operands[0], fobj->qs_obj.gmp_n);
		print_factors(fobj);

		break;
	case 53:
		//testrange - 4 arguments
		if (nargs != 4)
		{
			printf("wrong number of arguments in testrange\n");
			break;
		}

		{
			uint64 num_found;
			uint64 *primes;
			uint64 range;
			uint32 *sieve_p, num_sp;
			mpz_t lowz, highz;
			int val1, val2;
			
			range = mpz_get_64(operands[2]);
			val1 = PRIMES_TO_SCREEN;
			val2 = PRIMES_TO_FILE;
			PRIMES_TO_SCREEN = 0;
			PRIMES_TO_FILE = 0;
			primes = soe_wrapper(spSOEprimes, szSOEp, 0, range, 0, &num_found);
			PRIMES_TO_SCREEN = val1;
			PRIMES_TO_FILE = val2;
			sieve_p = (uint32 *)malloc(num_found * sizeof(uint32));
			for (i=0; i<num_found; i++)
				sieve_p[i] = (uint32)primes[i];
			num_sp = (uint32)num_found;
			free(primes);
			primes = NULL;

			mpz_init(lowz);
			mpz_init(highz);
			mpz_set(lowz, operands[0]);
			mpz_set(highz, operands[1]);
			primes = sieve_to_depth(sieve_p, num_sp, lowz, highz, 
				0, mpz_get_ui(operands[3]), &num_found);


			free(sieve_p);
			if (!NULL)
				free(primes);

			mpz_clear(lowz);
			mpz_clear(highz);
			mpz_set_ui(operands[0], num_found);
		}
		
		break;

	case 54:
		printf("siqstune not currently supported\n");
		break;

	case 55:
		//print a table of prime counts similar to http://www.trnicely.net/pi/pix_0000.htm
		//printf("ptable not currently supported\n");
		//break;
		//lower = 1000000000;
		//count = 50847534;
		lower = 10;
		count = 4;
		printf("%" PRIu64 ": %" PRIu64 "\n",lower,count);
		for (i = 1; i < 13; i++)
		{
			inc = (uint64)pow(10,i);
			//this increment may be too high.  break it into chunks no larger than 10e9.
			k = 10; 
			if (inc > 10000000000ULL)
			{
				k *= (uint32)(inc / 10000000000ULL);
				inc = 10000000000ULL;
			}

			for (j = 1; j < k; j++)
			{
				upper = lower + inc; 
				gettimeofday(&tstart, NULL);	
				soe_wrapper(spSOEprimes, szSOEp, lower, upper, 1, &n64);
				count += n64;
				gettimeofday (&tstop, NULL);
				t = my_difftime (&tstart, &tstop);
				lower = upper;
				printf("%" PRIu64 ": %" PRIu64 ", elapsed time = %6.4f\n",upper,count,t);
			}
		}

		break;

	case 56: 

		//sieverange - 4 arguments
        // range lo, range hi, sieve bound, count?
		if (nargs != 4)
		{
			printf("wrong number of arguments in sieverange\n");
			break;
		}

		{
			uint64 num_found;
			uint64 *primes;
			uint64 range;
			uint32 *sieve_p, num_sp;
			mpz_t lowz, highz;
			int val1, val2;
			
			range = mpz_get_64(operands[2]);
			val1 = PRIMES_TO_SCREEN;
			val2 = PRIMES_TO_FILE;
			PRIMES_TO_SCREEN = 0;
			PRIMES_TO_FILE = 0;
			primes = soe_wrapper(spSOEprimes, szSOEp, 0, range, 0, &num_found);
			PRIMES_TO_SCREEN = val1;
			PRIMES_TO_FILE = val2;
			sieve_p = (uint32 *)malloc(num_found * sizeof(uint32));
            for (i = 0; i < num_found; i++)
            {
                if ((num_found - i) < 10) printf("%u\n", primes[i]);
                sieve_p[i] = (uint32)primes[i];
            }
			num_sp = (uint32)num_found;
			free(primes);
			primes = NULL;

			mpz_init(lowz);
			mpz_init(highz);
			mpz_set(lowz, operands[0]);
			mpz_set(highz, operands[1]);
			primes = sieve_to_depth(sieve_p, num_sp, lowz, highz, 
				mpz_get_ui(operands[3]), 0, &num_found);

			free(sieve_p);
			if (!NULL)
				free(primes);

			mpz_clear(lowz);
			mpz_clear(highz);
			mpz_set_ui(operands[0], num_found);
		}
		
		break;

	case 57:
		//fermat - three arguments
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
		print_factors(fobj);
		break;

	case 58:
		//nfs - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in nfs\n");
			break;
		}
		mpz_set(fobj->N, operands[0]);
		mpz_set(fobj->nfs_obj.gmp_n, operands[0]);
		nfs(fobj);
		mpz_set(operands[0], fobj->nfs_obj.gmp_n);
		print_factors(fobj);
		break;

	case 59:
		//tune, no arguments
		factor_tune(fobj);
		break;

	case 60:
		// xor
		if (nargs != 2)
		{
			printf("wrong number of arguments in xor\n");
			break;
		}
		else
		{
			mpz_xor(operands[0], operands[0], operands[1]);
		}

		break;
	case 61:
		// and
		if (nargs != 2)
		{
			printf("wrong number of arguments in and\n");
			break;
		}
		else
		{
			mpz_and(operands[0], operands[0], operands[1]);
		}

		break;
	case 62:
		// or
		if (nargs != 2)
		{
			printf("wrong number of arguments in or\n");
			break;
		}
		else
		{
			mpz_ior(operands[0], operands[0], operands[1]);
		}

		break;
	case 63:
		// not
		if (nargs != 1)
		{
			printf("wrong number of arguments in not\n");
			break;
		}
		else
		{
			mpz_com(operands[0], operands[0]);
		}

		break;

	case 64:
		// factor a range of single precision numbers
		if (nargs != 2)
		{
			printf("wrong number of arguments in frange\n");
			break;
		}
		else
		{
			if (mpz_sizeinbase(operands[0],2) >= 64 || mpz_sizeinbase(operands[1],2) >= 64)
			{
				printf("inputs must be single precision\n");
				break;
			}
			
			if (mpz_cmp(operands[1], operands[0]) <= 0)
			{
				printf("second operand must be bigger than the first\n");
				break;
			}
			else
			{
				// split range into blocks of size N.
				// sieve each block with primes less than N/2 and store all factors found.
				// complete factorizations with factor() with trial division bound set to 0.
				// print stored factors together with those found with factor().
				uint64 b = mpz_get_64(operands[1]);
				uint64 a = mpz_get_64(operands[0]);
				uint64 diff = b-a;
                uint64 res;
				int j;
				int nblk, blk;
				int v;
				uint32 dlimit;
				int do_log;
				uint64 *n;
				uint64 **f;
				uint8 *nf;
				int np = 0;
                uint32 BLKSIZE = diff;
                mpz_t gmpn;
                ecm_params my_ecm_params;
                int num = 0, k = 0, numff = 0, numt = 0, numtd = 0, numtdff = 0, numsq = 0;
                
                ecm_init(my_ecm_params);
                mpz_init(gmpn);

				// remember original settings that will be changed
				v = VFLAG;
				dlimit = fobj->div_obj.limit;
				do_log = LOGFLAG;

				// change them
				VFLAG = -1;				
				LOGFLAG = 0;
				fobj->div_obj.limit = 0;

				// setup structures to hold inputs and factors of a block of numbers.
				// we could halve (and then some) the storage requirements by storing
				// unique prime factors and a count.  but for now this is good enough.
				n = (uint64 *)malloc(BLKSIZE * sizeof(uint64));
				f = (uint64 **)malloc(BLKSIZE * sizeof(uint64 *));
				nf = (uint8 *)malloc(BLKSIZE * sizeof(uint8));
				for (i=0; i<BLKSIZE; i++)
				{
					nf[i] = 0;
					f[i] = (uint64 *)malloc(64 * sizeof(uint64));
				}

				nblk = diff / BLKSIZE + (diff % BLKSIZE != 0);
				for (blk = 0; blk<nblk; blk++)
				{
					uint64 pid;
					uint64 p;          
                    uint64 state = 28953;

#define SPECIAL_PROBLEM
                    

					  
#ifdef SPECIAL_PROBLEM
                    // initialize the block
                    for (j = 0; j<BLKSIZE; j++)
                    {
                        nf[j] = 0;
                        state = 6364136223846793005ULL * state + 1442695040888963407ULL;
                        n[j] = state >> 1;
                    }
                    
                    // can't sieve, as the numbers are not in order.
                    // do trial division instead

                    for (j = 0; j < BLKSIZE; j++)
                    {
                        while ((n[j] & 1) == 0)
                        {
                            n[j] >>= 1;
                            nf[j]++;
                            numtd++;
                        }

                        k = 1;
                        while ((spSOEprimes[k] < 5000) && (n[j] > 1))
                        {
                            while ((n[j] % spSOEprimes[k]) == 0)
                            {
                                n[j] /= spSOEprimes[k];
                                nf[j]++;
                                numtd++;
                            }
                            k++;
                        }
                    }
                    
#else
                    // initialize the block
                    for (j = 0; j<BLKSIZE; j++)
                    {
                        nf[j] = 0;
                        n[j] = a + blk*BLKSIZE + j;
                    }

                    // sieve the block                  
					pid = 0;
					p = spSOEprimes[pid];
					while (p < 2*BLKSIZE)
					{
						uint32 j;
						uint32 offset;
						uint32 mod = (a + blk*BLKSIZE) % p;
						if (mod == 0)
							offset = 0;
						else
							offset = p - mod;

						for (j=offset; j<BLKSIZE; j += p)
						{
							do
							{
								n[j] /= p;
								f[j][nf[j]++] = p;
                                numtd++;
							} while ((n[j] % p) == 0);                            
						}
						pid++;
						p = spSOEprimes[pid];
					}
#endif

                    // parallel squfof on the remaining composites
#ifdef USE_VEC_SQUFOF
                    {
                        uint64 *m = (uint64 *)xmalloc_align(BLKSIZE * sizeof(uint64));
                        uint64 *fac = (uint64 *)xmalloc_align(BLKSIZE * sizeof(uint64));
                        int *idx = (int *)xmalloc_align(BLKSIZE * sizeof(int));
                        
                        for (k = 0, j = 0; j < BLKSIZE; j++)
                        {
                            if (n[j] > 1)
                            {
                                mpz_set_64(fobj->N, n[j]);
                                if (mpz_probab_prime_p(fobj->N, 1))
                                {
                                    // count primes... i.e., prime with no sieve factors
                                    if (nf[j] == 0)
                                        np++;
                                    nf[j]++;
                                    n[j] = 1;
                                    numtdff++;
                                }
                                else if (mpz_sizeinbase(fobj->N, 2) < 56)
                                {
                                    // queue up for batch squfof
                                    m[k] = n[j];
                                    idx[k++] = j;
                                }
                            }
                            else
                            {
                                numtdff++;
                            }
                        }

                        fprintf(stderr, "removed %d factors with sieving\n", numtd);

                        numsq = par_shanks_loop(m, fac, k);
                        fprintf(stderr, "factored %d out of %d with par_squfof\n", numsq, k);
                        for (j = 0; j < k; j++)
                        {
                            if (fac[j] > 1)
                            {
                                n[idx[j]] = 1;
                                nf[idx[j]] += 2;
                            }
                        }

                        align_free(m);
                        align_free(fac);
                        align_free(idx);
                    }       

#else
                    fprintf(stderr, "removed %d factors with sieving\n", numtd);
#endif

                    num = 0;
                    numff = 0;                    
                    numsq = 0;
                    state = 28953;

					// now complete the factorization and print all factors for each 
					// number in the block
					for (j=0; j<BLKSIZE; j++)
					{
						//if (a + blk*BLKSIZE + j > b)
						//	break;

                        mpz_set_64(fobj->N, n[j]);
#ifdef SPECIAL_PROBLEM
                        state = 6364136223846793005ULL * state + 1442695040888963407ULL;
                        printf("%" PRIu64 " ", state>>1);

#else                        
                        printf("%" PRIu64 " ", a + blk*BLKSIZE + j);
#endif
                        fflush(stdout);

#ifndef USE_VEC_SQUFOF
                       
                        // check if prime
                        if (n[j] == 1)
                        {
                            numtdff++;
                        }
                        else if (mpz_probab_prime_p(fobj->N, 1))
                        {
                            // count primes... i.e., prime with no sieve factors
                            if (nf[j] == 0)
                                np++;
                            nf[j]++;
                            n[j] = 1;
                            numtdff++;
                        }
                        else 
#endif                        
                        {
                            numt++;

                            // still unfactored?
                            if (n[j] > 1)
                            {                                
                                // do some rho
                                for (k = 0; k < 3; k++)
                                {
                                    res = spbrent(n[j], k+1, 8192);
                                    if (res > 0)
                                    {
                                        nf[j]++;
                                        n[j] /= res;
                                        
                                        mpz_set_ui(fobj->N, n[j]);
                                        if (mpz_probab_prime_p(fobj->N, 1))
                                        {
                                            nf[j]++;
                                            n[j] = 1;
                                            num++;
                                            break;
                                        }
                                    }                                    
                                }
                            }
                            
                            // still unfactored?
                            if (n[j] > 1)
                            {
                                // do some ecm
                                // printf("trying ecm... ");
                                for (k = 0; k < 20; k++)
                                {
                                    uint64 sigma = spRand(100, 1000000000);
                                    my_ecm_params->B1done = 1.0 + floor(1 * 128.) / 134217728.;
                                    mpz_set_ui(my_ecm_params->x, (unsigned long)0);
                                    mpz_set_ui(my_ecm_params->sigma, sigma);
                                    my_ecm_params->method = ECM_ECM;
                                    mpz_set_ui(fobj->ecm_obj.gmp_n, n[j]);
                                    ecm_factor(fobj->ecm_obj.gmp_f, fobj->ecm_obj.gmp_n, 2000, my_ecm_params);
                                    if ((mpz_cmp_ui(fobj->ecm_obj.gmp_f, 1) > 0) &&
                                        (mpz_cmp_ui(fobj->ecm_obj.gmp_f, n[j]) < 0))
                                    {
                                        nf[j]++;
                                        n[j] /= mpz_get_ui(fobj->ecm_obj.gmp_f);

                                        mpz_set_ui(fobj->ecm_obj.gmp_n, n[j]);

                                        if (mpz_probab_prime_p(fobj->ecm_obj.gmp_n, 1))
                                        {
                                            nf[j]++;
                                            n[j] = 1;
                                            numff++;
                                            break;
                                        }
                                    }
                                }
                            }

                            // still unfactored?
                            if (n[j] > 1)
                            {
                                // give up
                                printf(" giving up... ");
                                fflush(stdout);
                            }
                        }                        

                        printf("has %d factors\n", nf[j]);
                        fflush(stdout);


                    /*
                        printf("%" PRIu64 " = ", n[j]);
                        for (k=0; k<nf[j]; k++)
                        {
                            printf("%" PRIu64 "", f[j][k]);
                            if (k != nf[j] - 1)
                                printf(" * ");
                        }

                        if (n[j] == 1)
                        {
                            printf("\n");
                            continue;
                        }

						if (is_mpz_prp(fobj->N))
						{
							if (nf[j] == 0)
							{
								np++;
								printf("is prime\n");
							}
							else
								printf(" * %" PRIu64 "\n", n[j]);
						}
						else
						{
							// input not prime and contains no small factors.
							// try squfof first										
							uint64 res;
                            
							mpz_set_ui(gmpn, n[j]);
							res = sp_shanks_loop(gmpn, fobj);
							if (res > 1)
							{
								printf(" * %" PRIu64 "", res);
								n[j] /= res;
							}

							mpz_set_64(fobj->N, n[j]);
							if (is_mpz_prp(fobj->N))
							{
								printf(" * %" PRIu64 "\n", n[j]);
								continue;
							}

							// rarely, squfof+trial won't complete the job, so pull out 
							// the big guns.
							factor(fobj);
							if (fobj->num_factors > 0)
								printf(" * ");

							for (k=0; k<fobj->num_factors; k++)
							{
								gmp_printf("%Zd", fobj->fobj_factors[k].factor);
								if (fobj->fobj_factors[k].count > 1)
									printf("^%d", fobj->fobj_factors[k].count);
								if (k != fobj->num_factors-1)
									printf(" * ");
							}
							printf("\n");
							clear_factor_list(fobj);
							fobj->refactor_depth = 0;
						}

                        */
					}
				}
				//printf("found %d primes in range\n", np);
                fprintf(stderr, "factored %d out of %d with trial division\n", numtdff, diff);
#ifndef USE_VEC_SQUFOF
                fprintf(stderr, "factored %d out of %d with squfof\n", numsq, numt);
#endif
                fprintf(stderr, "factored %d out of %d with rho\n", num, numt);
                fprintf(stderr, "factored %d out of %d with ecm\n", numff, numt);

				// restore settings
				VFLAG = v;
				LOGFLAG = do_log;
				fobj->div_obj.limit = dlimit;

				// free memory
				free(nf);
				free(n);
				for (i=0; i<BLKSIZE; i++)
					free(f[i]);
				free(f);

                mpz_clear(gmpn);
			}

			mpz_set_ui(operands[0], 0);
		}

		break;

	case 65:
		/* bpsw */
		if (nargs != 1)
		{
			printf("wrong number of arguments in bpsw\n");
			break;
		}
		else
		{
			int ret = 0;
			ret = mpz_bpsw_prp(operands[0]);
			// should we print out the input number again?
			if (ret == PRP_COMPOSITE)
			{
				printf("Input is composite.  C%d\n", gmp_base10(operands[0]));
			}
			else if (ret == PRP_PRP)
			{
				printf("Input is prp.  PRP%d\n", gmp_base10(operands[0]));
			}
			else if (ret == PRP_PRIME)
			{
				printf("Input is prime.  P%d\n", gmp_base10(operands[0]));
			}
		}
		break;

	case 66:
		/* aprcl */
		if (nargs != 1)
		{
			printf("wrong number of arguments in aprcl\n");
			break;
		}
		else
		{
			int ret = 0;
			ret = mpz_aprtcle(operands[0], APRTCLE_VERBOSE1);
			// should we print out the input number again?
			if (ret == APRTCLE_COMPOSITE)
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
			else if (ret == APRTCLE_PRP)
			{
				printf("\nInput is prp.  PRP%d\n", gmp_base10(operands[0]));
			}
			else if (ret == APRTCLE_PRIME)
			{
				printf("\nInput is prime.  P%d\n", gmp_base10(operands[0]));
			}
		}
		break;

	case 67:
		// lte
		if (nargs != 2)
		{
			printf("wrong number of arguments in lte\n");
			break;
		}
		else
		{
			mpz_set_ui(operands[0], mpz_cmp(operands[0], operands[1]) <= 0);
		}
		break;

	case 68:
		// gte
		if (nargs != 2)
		{
			printf("wrong number of arguments in gte\n");
			break;
		}
		else
		{
			mpz_set_ui(operands[0], mpz_cmp(operands[0], operands[1]) >= 0);
		}
		break;

	case 69:
		// lt
		if (nargs != 2)
		{
			printf("wrong number of arguments in lt\n");
			break;
		}
		else
		{
			mpz_set_ui(operands[0], mpz_cmp(operands[0], operands[1]) < 0);
		}
		break;

	case 70:
		// gt
		if (nargs != 2)
		{
			printf("wrong number of arguments in gt\n");
			break;
		}
		else
		{
			mpz_set_ui(operands[0], mpz_cmp(operands[0], operands[1]) > 0);
		}

		break;
    case 71:
        // = (assignment)
        if (nargs != 2)
        {
            printf("wrong number of arguments in =\n");
            break;
        }
        else
        {
            // the first operands is a string and should therefore
            // be in the choperands array.
            if (invalid_dest(choperands[0]))
            {
                printf("invalid destination %s\n", choperands[0]);
                break;
            }

            if (set_uvar(choperands[0], operands[1]))
            {
                new_uvar(choperands[0], operands[1]);                
            }
            mpz_set(operands[0], operands[1]);
        }

        break;
    case 72:
        // if
        if (nargs == 2)
        {
            // just the if case
            char cond[80], ifbody[80];
            int nooutput;

            get_strvar(choperands[1], cond);
            get_strvar(choperands[2], ifbody);

            if (CALC_VERBOSE)
                fprintf(stderr, "Processing if with cond = %s, body = %s\n",
                cond, ifbody);

            // execute cond expression
            process_expression(cond, fobj, 1);
            get_uvar("ans", operands[4]);

            if (mpz_get_ui(operands[4]) != 0)
            {
                if (ifbody[strlen(ifbody) - 1] == ';')
                {
                    nooutput = 1;
                    ifbody[strlen(ifbody) - 1] = '\0';
                }
                else
                    nooutput = 0;

                process_expression(ifbody, fobj, nooutput);
            }
        }
        else if (nargs == 3)
        {
            char cond[80], ifbody[80], elsebody[80];
            int nooutput;

            get_strvar(choperands[0], cond);
            get_strvar(choperands[1], ifbody);
            get_strvar(choperands[2], elsebody);

            if (CALC_VERBOSE)
                fprintf(stderr, "Processing if with cond = %s, body = %s, elsebody = %s\n",
                cond, ifbody, elsebody);

            // execute cond expression
            process_expression(cond, fobj, 1);
            get_uvar("ans", operands[4]);

            if (mpz_get_ui(operands[4]) != 0)
            {
                if (ifbody[strlen(ifbody) - 1] == ';')
                {
                    nooutput = 1;
                    ifbody[strlen(ifbody) - 1] = '\0';
                }
                else
                    nooutput = 0;

                process_expression(ifbody, fobj, nooutput);
            }
            else
            {
                if (elsebody[strlen(elsebody) - 1] == ';')
                {
                    nooutput = 1;
                    elsebody[strlen(elsebody) - 1] = '\0';
                }
                else
                    nooutput = 0;

                process_expression(elsebody, fobj, nooutput);
            }
        }
        else
        {
            printf("wrong number of arguments in if\n");
            break;
        }

        break;
    case 73:
        // print
        // todo: honor the current OBASE setting
        gmp_printf("%Zd\n", operands[0]);

        break;

    case 74:
        // for
        if (nargs != 4)
        {
            printf("wrong number of arguments in for loop\n");
            break;
        }

        {
            char init[80], test[80], iter[80], body[80];
            int nooutput;

            get_strvar(choperands[0], init);
            get_strvar(choperands[1], test);
            get_strvar(choperands[2], iter);
            get_strvar(choperands[3], body);

            if (CALC_VERBOSE)
                fprintf(stderr, "Processing loop with init = %s, test = %s, iter = %s, body = %s\n",
                init, test, iter, body);

            // look for a trailing semicolon
            if (body[strlen(body) - 1] == ';')
            {
                nooutput = 1;
                body[strlen(body) - 1] = '\0';
            }
            else
                nooutput = 0;

            // execute init expression
            //toStr(init, &str);
            process_expression(init, fobj, 1);

            while (1)
            {
                // execute body expression.
                //toStr(body, &str);
                //calc_with_assignment(&str, fobj, nooutput);
                process_expression(body, fobj, nooutput);

                // execute iter expression.
                //toStr(iter, &str);
                //calc_with_assignment(&str, fobj, 1);
                process_expression(iter, fobj, 1);

                // execute test expression.
                //toStr(test, &str);
                //calc_with_assignment(&str, fobj, 1);
                process_expression(test, fobj, 1);
                get_uvar("ans", operands[4]);
                if (mpz_get_ui(operands[4]) == 0)
                    break;
            }
        }

        break;

    case 75:
        // forprime
        if (nargs != 3)
        {
            printf("wrong number of arguments in forprime loop\n");
            break;
        }

        {
            char start[80], stop[80], body[80], name[80], *ptr;
            uint64 ustart, ustop;
            uint64 *plist, nump;
            int nooutput;

            get_strvar(choperands[0], start);
            get_strvar(choperands[1], stop);
            get_strvar(choperands[2], body);

            if (CALC_VERBOSE)
                fprintf(stderr, "Processing forprime loop with start = %s, stop = %s, body = %s\n",
                start, stop, body);

            // the start expression should include a '=' and a valid variable name
            ptr = strtok(start, "=");
            if (ptr == NULL)
            {
                printf("badly formatted forprime loop: forprime(var=start; stop; body)\n");
                exit(3);
            }
            strcpy(name, ptr);
            ptr = strtok(NULL, "\0");
            if (ptr == NULL)
            {
                printf("badly formatted forprime loop: forprime(var=start; stop; body)\n");
                exit(3);
            }
            strcpy(start, ptr);
            mpz_set_str(operands[0], start, 0);
            mpz_set_str(operands[1], stop, 0);

            if (set_uvar(name, operands[0]))
                new_uvar(name, operands[0]);

            if ((mpz_cmp_ui(operands[0], 0xffffffffffffffffULL) > 0) ||
                (mpz_cmp_ui(operands[1], 0xffffffffffffffffULL) > 0))
            {
                printf("start and stop must be < 2^64 in forprime loop\n");
                exit(1);
            }

            ustart = mpz_get_ui(operands[0]);
            ustop = mpz_get_ui(operands[1]);

            plist = soe_wrapper(spSOEprimes, szSOEp, ustart, ustop, 0, &nump);

            // look for a trailing semicolon
            if (body[strlen(body) - 1] == ';')
            {
                nooutput = 1;
                body[strlen(body) - 1] = '\0';
            }
            else
                nooutput = 0;

            for (i = 0; i < nump; i++)
            {
                mpz_set_ui(operands[0], plist[i]);
                set_uvar(name, operands[0]);

                process_expression(body, fobj, nooutput);
            }
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
	//create a new user variable with name 'name', and return
	//its location in the global uvars structure
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
	//look for 'name' in the global uvars structure
	//if found, copy in data and return 0
	//else return 1
	int i;

	i = mpz_get_ui(data);
	//first look if it is a global constant
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
	else if (strcmp(name,"LOGFLAG") == 0) {
		LOGFLAG = i; return 0;}
	else if (strcmp(name,"VFLAG") == 0) {
		VFLAG = i; return 0;}
	else if (strcmp(name,"PRIMES_TO_FILE") == 0) {
		PRIMES_TO_FILE = i; return 0;}
	else if (strcmp(name,"PRIMES_TO_SCREEN") == 0) {
		PRIMES_TO_SCREEN = i; return 0;}
	else if (strcmp(name,"NUM_WITNESSES") == 0) {
		NUM_WITNESSES = i; return 0;}

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
	//look for 'name' in the global uvars structure
	//if found, copy out data and return 0
	//else return 1 if not found
	int i;

	//first look if it is a global constant
	if (strcmp(name,"IBASE") == 0) {
		mpz_set_ui(data, IBASE); return 0;}
	else if (strcmp(name,"OBASE") == 0) {
		mpz_set_ui(data, OBASE); return 0;}
	else if (strcmp(name,"NUM_WITNESSES") == 0) {
		mpz_set_ui(data, NUM_WITNESSES); return 0;}
	else if (strcmp(name,"LOGFLAG") == 0) {
		mpz_set_ui(data, LOGFLAG); return 0;}
	else if (strcmp(name,"VFLAG") == 0) {
		mpz_set_ui(data, VFLAG); return 0;}
	else if (strcmp(name,"PRIMES_TO_FILE") == 0) {
		mpz_set_ui(data, PRIMES_TO_FILE); return 0;}
	else if (strcmp(name,"PRIMES_TO_SCREEN") == 0) {
		mpz_set_ui(data, PRIMES_TO_SCREEN); return 0;}

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
		printf("NUM_WITNESSES      %u\n",NUM_WITNESSES);
		printf("LOGFLAG            %u\n",LOGFLAG);
		printf("VFLAG              %u\n",VFLAG);
		printf("PRIMES_TO_FILE     %u\n",PRIMES_TO_FILE);
		printf("PRIMES_TO_SCREEN   %u\n",PRIMES_TO_SCREEN);

		for (i=0;i<uvars.num;i++)
			printf("%s      %s\n",uvars.vars[i].name,mpz_get_str(gstr1.s, 10, uvars.vars[i].data));

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
    //create a new user variable with name 'name', and return
    //its location in the global uvars structure
    if (strvars.num == strvars.alloc)
    {
        //need more room for variables
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
    for (i = 0; i<strvars.alloc; i++)
        free(strvars.vars[i].data);
    free(strvars.vars);
}


