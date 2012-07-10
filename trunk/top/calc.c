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

char opchar[10] = {'=','<','>','+','-','*','/','%','^'};
char imms[3] = {'!','#','-'};
const int numopchars = 9;
z operands[5];

int calc_init()
{
	int i;
	//user variables space
	uvars.vars = (uvar_t *)malloc(10 * sizeof(uvar_t));
	uvars.alloc = 10;
	for (i=0;i<uvars.alloc;i++)
		zInit(&uvars.vars[i].data);
	strcpy(uvars.vars[0].name,"ans");
	uvars.num = 1;
	return 1;
}

void calc_finalize()
{
	free_uvars();
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

int op_precedence(char *s1, char *s2, int assoc)
{
	//if associativity is RIGHT, then use strictly greater than
	//else use greater than or equal to
	int p1=0,p2=0;

	if (strcmp(s1,"=") == 0) p1=-1;
	if (strcmp(s1,"<<") == 0) p1=0;
	if (strcmp(s1,">>") == 0) p1=0;
	if (strcmp(s1,"+") == 0) p1=1;
	if (strcmp(s1,"-") == 0) p1=1;
	if (strcmp(s1,"*") == 0) p1=2;
	if (strcmp(s1,"/") == 0) p1=2;
	if (strcmp(s1,"%") == 0) p1=2;
	if (strcmp(s1,"^") == 0) p1=3;
	if (strcmp(s1,"\\") == 0) p1=4;

	if (strcmp(s1,"=") == 0) p1=-1;
	if (strcmp(s2,"<<") == 0) p2=0;
	if (strcmp(s2,">>") == 0) p2=0;
	if (strcmp(s2,"+") == 0) p2=1;
	if (strcmp(s2,"-") == 0) p2=1;
	if (strcmp(s2,"*") == 0) p2=2;
	if (strcmp(s2,"/") == 0) p2=2;
	if (strcmp(s2,"%") == 0) p2=2;
	if (strcmp(s2,"^") == 0) p2=3;
	if (strcmp(s2,"\\") == 0) p2=4;

	if (assoc == LEFT)
		return p1 >= p2;
	else 
		return p1 > p2;
}

int preprocess(str_t *in)
{
	//preprocess the expression in 'in'
	int i,j;

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

	//look for a^b%c combinations, and replace with a function call
	//or (a^b)%c
	//should ignore stuff inside parens, so for instance
	//(3+5)^(8*9*7)%(134513%46) should reduce to
	//modexp(3+5,8*9*7,134513%46)

	//algebraic simplification (this would be cool...)

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
	//take a string as input
	//break it into tokens
	//create an array of strings for each token
	//return the pointer to the array and the number of elements
	//in the array.  this will all have to be freed later
	//by the caller

	//a token in this context is one of the following things:
	//	a number, possibly including a base prefix (0x, 0d, 0b, etc)
	//	a variable name
	//	a function name
	//	an operator string (includes parens, commas)

	//read the string one character at a time
	//for each character read, decide if we've found the start of a new token

	int inpos, i, el_type, el_type2, token_alloc, tmpsize = GSTR_MAXSIZE;
	char ch;
	char *tmp;
	char **tokens;

	token_alloc = 100;		//100 tokens
	tokens = (char **)malloc(token_alloc * sizeof(char *));
	*num_tokens = 0;

	tmp = (char *)malloc(GSTR_MAXSIZE * sizeof(char));

	//get the first character and check the type
	inpos = 0;
	i=1;
	ch = in[inpos];
	tmp[i-1] = ch;
	el_type = get_el_type2(ch);

	//when an expression gets cast into postfix, it aquires a leading
	//space which we can skip here
	if (el_type == SPACE)
	{
		inpos = 1;
		i=1;
		ch = in[inpos];
		tmp[i-1] = ch;
		el_type = get_el_type2(ch);
	}

	//ambiguous types:
	//a "-" can be either a num (if a negative sign) or an operator
	//a number can be a number or a string (num or func/var name)
	//a letter can be a string or a number (hex or func/var name)
	//we can tell them apart from the surrounding context
	//
	//negative signs never have a num type before (or anything that
	//can be evaluated as a num , i.e. ")"
	//
	//if we are reading CH's don't stop interpreting them as CH's until
	//we find a non-CH or non-NUM (use a flag)
	//
	//watch for magic combinations '0x' '0d', etc.  set a flag to 
	//interpret what follows as the appropriate kind of num, and
	//discard the '0x', etc.
	//once this is fixed, change the final stack evaluation to print hex
	//strings to save some conversion time.
	if (el_type == AMBIG)
	{
		if (get_el_type2(in[inpos+1]) == NUM)
			el_type = NUM;
		else
			el_type = OP;
	}
	while (1)
	{
		//get another character and check the type
		inpos++;
		//if el_type == EOE, then no reason to keep reading.  This bug didn't seem to cause
		//any crashes, but couldn't have been healthy...
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
				//when processing postfix strings, we need this
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
				//unrecognized character.  clear all tokens and return;
				printf("unrecognized character in input\n");
				for (i=0;i< *num_tokens; i++)
					free(tokens[i]);
				free(tokens);
				free(tmp);
				return NULL;
			}

			if (el_type != SPACE)
			{
				//create a new token
				tmp[i] = '\0';
				tokens[*num_tokens] = (char *)malloc((strlen(tmp) + 1) * sizeof(char));
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
		
			//then cycle the types
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
	if (str[0] == '!' || str[0] == '#') return 1;
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

	/* if there is no input to process, we are done... */
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
		/* if token is a number, push it onto the stack... */
		if (isNumber(tok))
		{
			toStr(tok, tmp);
			push(tmp, &stk);
		}
		else if (isOperator(tok))
		{
			if (tok[0] == '!' || tok[0] == '#')
			{
				/* factorial and primorial are unary operators */
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
				/* pop off two elements a,b
				   push back on "(aOPb)" */
				if ((pop(tmp2, &stk) == 0) || (pop(tmp1, &stk) == 0))
				{
					invalid_string = 1;
					break; /* bad input string, don't parse any more... */
				}
				sClear(tmp);
				sAppend(tmp1->s, tmp);
				sAppend(tok, tmp);
				sAppend(tmp2->s, tmp);
				push(tmp, &stk);
			}
			else if (tok[0] == '-')
			{
				/* pop off two elements a,b
				   push back on "(aOPb)" */
				if ((pop(tmp2, &stk) == 0) || (pop(tmp1, &stk) == 0))
				{
					invalid_string = 1;
					break; /* bad input string, don't parse any more... */
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
				/* pop off two elements a,b
				   push back on "(aOPb)" */
				if ((pop(tmp2, &stk) == 0) || (pop(tmp1, &stk) == 0))
				{
					invalid_string = 1;
					break; /* bad input string, don't parse any more... */
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
			/* non-number and non-operator encountered, we're done grabbing the input string */
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
	z tmpz;
	
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
	for (i=0;i<5;i++)
		zInit(&operands[i]);
	zInit(&tmpz);
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
			varstate = get_uvar(tokens[i],&tmpz);
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
				//not a variable and not valid num
				printf("unrecognized token: %s\n",tokens[i]);
				retval=1;
				goto free;
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

	//process the output postfix expression:
	//this can be done with a simple stack
	//all tokens are separated by spaces
	//all tokens consist of numbers or functions

	/* try to grab the input number as an expression, if found, store it in fobj->str_N */
	strN = strdup(post->s); /* make a copy of the postfix string */
	if (strN != NULL)
	{
		/* find input expression and store in str_N */
		get_expression(strN, &(fobj->str_N));
		//printf("Found expression: %s\n", fobj->str_N.s);
		free(strN);
	}

	//now evaluate the RPN expression
	//printf("processing postfix expression: %s\n",post->s);
	delim[0] = ' ';
	delim[1] = '\0';
	tok = strtok(post->s,delim);
	if (tok == NULL)
	{
		printf("nothing to process\n");
		goto free;
	}

	do
	{
		//printf("stack contents: ");
		//for (i=0; i<stk.num; i++)
		//	printf("%s ",stk.elements[i]->s);
		//printf("\n");
		//printf("current token: %s\n\n", &tok[0]);

		switch (get_el_type2(tok[0]))
		{
		case NUM:
			//printf("pushing %s onto stack\n",tok);
			toStr(tok,tmp);
			push(tmp,&stk);
			break;
		case AMBIG:
			//could be a number or a function
			//if the next character is a number, this is too
			//else its an operator
			//printf("peeking at tok + 1: %d\n",(int)tok[1]);
			if (get_el_type2(tok[1]) == NUM)
			{
				//printf("pushing %s onto stack\n",tok);
				toStr(tok,tmp);
				push(tmp,&stk);
				break;
			}

			//if not a num, proceed into the next switch (function handle)
		default:
			func = getFunc(tok,&na);
			//printf("processing function %d\n",func);

			if (func >= 0)
			{
				//pop those args and put them in a global array
				for (j=0;j<na;j++)
				{
					//right now we must get all of the operands.
					//somewhere in there we should make allowances
					//for getting a reduced number (i.e. for unary "-"
					//and for variable numbers of arguments
					k = pop(tmp,&stk);

					//try to make a number out of it
					//printf("looking at argument %s\n",tmp->s);
					str2hexz(tmp->s,&tmpz);
					//printf("found numerical argument %s\n",z2decstr(&tmpz,&gstr1));

					//tmpz with size zero means this thing isn't a number
					if (tmpz.size == 0 || k == 0)
					{
						//didn't get the expected number of arguments
						//for this function.  This may be ok, if the
						//function accepts varable argument lists.
						//feval will handle it.

						//put the non-num back
						if (k != 0)
							push(tmp,&stk);

						//and bail
						break;
					}
					else
					{
						//it is a number, put it in the operand pile
						str2hexz(tmp->s,&operands[na-j-1]);
					}
				}

				na = j;
				//call the function evaluator with the 
				//operator string and the number of args available
				na = feval(func,na,fobj);

				//put result back on stack
				for (j=0;j<na;j++)
				{
					z2hexstr(&operands[j],tmp);
					push(tmp,&stk);
				}
			}
			else
			{
				get_uvar(tok,&tmpz);
				z2decstr(&tmpz,&gstr1);
				sCopy(tmp,&gstr1);
				push(tmp,&stk);
			}
			break;
		}

		tok = strtok((char *)0,delim);
	} while (tok != NULL);
	pop(in,&stk);


free:
	for (i=0;i<5;i++)
		zFree(&operands[i]);
	free(token_types);
	stack_free(&stk);
	zFree(&tmpz);
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

	char func[NUM_FUNC][11] = {"fib","luc","dec2hex","hex2dec","rsa",
						"gcd","jacobi","factor","rand","lg2",
						"log","ln","pm1","pp1","rho",
						"trial","mpqs","nextprime","size","issquare",
						"isprime","squfof","sqrt","modinv","modexp",
						"nroot","shift","siqs","primes", "ispow", 
						"torture","randb", "ecm","+","-",
						"*","/","!","#","==",
						"<<",">>","%","^","test",
						"puzzle","sieve","algebraic","llt","siqsbench",
						"pullp","sftest","smallmpqs","testrange","siqstune",
						"ptable","sieverange","fermat","nfs","tune",
						"xor", "and", "or", "not"};

	int args[NUM_FUNC] = {1,1,1,1,1,
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
					0,4,2,1,0,
					2,2,2,1};

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
	//evaluate the function 'fname', with argument(s) 'in'
	z mp1, mp2, mp3, tmp1, tmp2;
	mpz_t gmpz;

	str_t str;
	uint32 i=0;
	uint64 n64;
	uint32 j,k;
	//clock_t start, stop;
	uint32 *offsets;
	double t;
	struct timeval tstart, tstop;
	TIME_DIFF *	difference;
	uint64 lower, upper, inc, count;
	//FILE *out;
	//uint64 powof2, powof2m1, summodp;


	zInit(&mp1);
	zInit(&mp2);
	zInit(&mp3);
	zInit(&tmp1);
	zInit(&tmp2);
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
		i = operands[0].val[0];
		fib(i,&operands[0]);
		break;
	case 1:
		//luc - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in luc\n");
			break;
		}
		i = operands[0].val[0];
		lucas(i,&operands[0]);
		break;
		
	case 2:
		//dec2hex - one argument
		//calc(in);
		//str2hexz(in->s,&mp2);
		//z2hexstr(&mp2,in);
		break;
	case 3:
		//hex2dec - one argument
		//calc(in);
		//str2hexz(in->s,&mp1);
		//zHex2Dec(&mp1,&mp2);
		//z2decstr(&mp2,in);
		break;
	case 4:
		//rsa - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in gcd\n");
			break;
		}
		sp2z(2048,&mp1);
		sp2z(4096,&mp2);
		if (zCompare(&operands[0],&mp2) > 0)
		{
			printf("bitlength too large");
			sp2z(1,&operands[0]);
			break;
		}
		else if (zCompare(&operands[0],&mp1) > 0)
		{
			printf("Paranoid, huh?  This might take a minute\n");
		}
		i = operands[0].val[0];
		build_RSA(i,&operands[0]);
		break;
	case 5:
		//gcd - two arguments
		if (nargs != 2)
		{
			printf("wrong number of arguments in gcd\n");
			break;
		}

		{
			mpz_t gmp1, gmp2;

			mpz_init(gmp1);
			mpz_init(gmp2);
			mp2gmp(&operands[0], gmp1);
			mp2gmp(&operands[1], gmp2);
			mpz_gcd(gmp1, gmp1, gmp2);
			gmp2mp(gmp1, &operands[0]);
			mpz_clear(gmp1);
			mpz_clear(gmp2);
		}

		break;
	case 6:
		//jacobi - two arguments
		if (nargs != 2)
		{
			printf("wrong number of arguments in jacobi\n");
			break;
		}
		mp1.val[0] = zJacobi(&operands[0],&operands[1]);
		if (mp1.val[0] == MAX_DIGIT)
		{
			mp1.size *= -1;
			mp1.val[0] = 1;
		}
		zCopy(&mp1,&operands[0]);
		break;
	case 7:
		//factor - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in factor\n");
			break;
		}
		zCopy(&operands[0],&fobj->N);
		factor(fobj);
		zCopy(&fobj->N,&operands[0]);
		print_factors(fobj);
		break;
	case 8:
		//rand - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in rand\n");
			break;
		}
		i = operands[0].val[0];
		zRand(&operands[0],i);
		break;
	case 9:
		//lg2 - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in lg2\n");
			break;
		}
		mp1.val[0] = zBits(&operands[0]);
		mp1.size = 1;
		zCopy(&mp1,&operands[0]);
		break;
	case 10:
		//log - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in log\n");
			break;
		}
		mp1.val[0] = ndigits(&operands[0]);
		mp1.size = 1;
		zCopy(&mp1,&operands[0]);
		break;
	case 11:
		//ln - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in ln\n");
			break;
		}
		mp1.val[0] = (uint32)((double)(zBits(&operands[0]) - 1) * log(2.0));
		mp1.size = 1;
		zCopy(&mp1,&operands[0]);
		break;
	case 12:
		//pm1 - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in pm1\n");
			break;
		}
		zCopy(&operands[0],&fobj->N);
		mp2gmp(&operands[0],fobj->pm1_obj.gmp_n);
		pollard_loop(fobj);
		gmp2mp(fobj->pm1_obj.gmp_n,&operands[0]);
		print_factors(fobj);
		break;
	case 13:
		//pp1 - two arguments
		zCopy(&operands[0],&fobj->N);
		if (nargs == 2)
		{
			mp2gmp(&operands[0],fobj->pp1_obj.gmp_n);
			fobj->pp1_obj.numbases = operands[1].val[0];
			williams_loop(fobj);
			gmp2mp(fobj->pp1_obj.gmp_n,&operands[0]);
		}
		else if (nargs == 1)
		{
			mp2gmp(&operands[1],fobj->pp1_obj.gmp_n);
			fobj->pp1_obj.numbases = 1;
			williams_loop(fobj);
			gmp2mp(fobj->pp1_obj.gmp_n,&operands[1]);
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
		zCopy(&operands[0],&fobj->N);
		mp2gmp(&operands[0],fobj->rho_obj.gmp_n);
		brent_loop(fobj);
		gmp2mp(fobj->rho_obj.gmp_n,&operands[0]);
		print_factors(fobj);
		break;
	case 15:
		//trial - two arguments
		zCopy(&operands[0],&fobj->N);
		if (nargs == 2)
		{
			mp2gmp(&operands[0],fobj->div_obj.gmp_n);
			fobj->div_obj.print = 0;
			fobj->div_obj.limit = operands[1].val[0];
			zTrial(fobj);
			gmp2mp(fobj->div_obj.gmp_n,&operands[0]);
		}
		else if (nargs == 1)
		{
			printf("using default trial division bound of 10000\n");
			mp2gmp(&operands[1],fobj->div_obj.gmp_n);
			fobj->div_obj.print = 0;
			fobj->div_obj.limit = 10000;
			zTrial(fobj);
			gmp2mp(fobj->div_obj.gmp_n,&operands[1]);
			//apparently this comes in as operand 1, but the calling function
			//expects the result in operand 0, so put it there.  This should be
			//cleaner
			zCopy(&operands[1],&operands[0]);
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
			zNextPrime(&operands[0],&mp3,operands[1].val[0]);
			zCopy(&mp3,&operands[0]);
		}
		else if (nargs == 1)
		{
			//assume larger
			zNextPrime(&operands[1],&mp3,1);
			zCopy(&mp3,&operands[0]);
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
		mpz_init(gmpz);
		mp2gmp(&operands[0], gmpz);

		printf("%d digits, %d bits\n",gmp_base10(gmpz),zBits(&operands[0]));
		zCopy(&mp3,&operands[0]);
		mpz_clear(gmpz);

		break;
	case 19:
		//issquare
		if (nargs != 1)
		{
			printf("wrong number of arguments in issquare\n");
			break;
		}
		
		mpz_init(gmpz);
		mp2gmp(&operands[0], gmpz);
		i = mpz_perfect_square_p(gmpz);

		if (i)
			printf("input is square\n");
		else
			printf("input is not square\n");

		sp2z(i,&operands[0]);
		mpz_clear(gmpz);

		break;
	case 20:
		//isprime - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in isprime\n");
			break;
		}

		i = isPrime(&operands[0]);

		if (i)
			printf("probably prime\n");
		else
			printf("not prime\n");

		sp2z(i,&operands[0]);

		break;
	case 21:
		//shanks - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in squfof\n");
			break;
		}
		mpz_init(gmpz);

		mp2gmp(&operands[0], gmpz);
		n64 = sp_shanks_loop(gmpz,fobj);
		print_factors(fobj);
		sp642z(n64,&operands[0]);
		break;
	case 22:
		//sqrt - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in sqrt\n");
			break;
		}
		zNroot(&operands[0],&mp1,2);
		zCopy(&mp1,&operands[0]);
		break;
	case 23:
		//modinv - two arguments
		if (nargs != 2)
		{
			printf("wrong number of arguments in modinv\n");
			break;
		}
		xGCD(&operands[0],&operands[1],&tmp1,&tmp2,&mp1);
		zCopy(&tmp1,&operands[0]);
		break;
	case 24:
		//modexp - three arguments
		if (nargs != 3)
		{
			printf("wrong number of arguments in modexp\n");
			break;
		}
		zModExp(&operands[0],&operands[1],&operands[2],&tmp1);
		zCopy(&tmp1,&operands[0]);
		break;
	case 25:
		//nroot - two arguments
		if (nargs != 2)
		{
			printf("wrong number of arguments in Nroot\n");
			break;
		}
		zNroot(&operands[0],&mp3,operands[1].val[0]);
		zCopy(&mp3,&operands[0]);
		break;
	case 26:
		//shift - two arguments
		if (nargs != 2)
		{
			printf("wrong number of arguments in shift\n");
			break;
		}
		j = operands[1].val[0];

		if (operands[1].size < 0)
			zShiftRight(&mp3,&operands[0],-1*j);
		else
			zShiftLeft(&mp3,&operands[0],j);
		zCopy(&mp3,&operands[0]);
		break;
	case 27:
		//siqs - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in siqs\n");
			break;
		}
		zCopy(&operands[0],&fobj->N);
		mp2gmp(&operands[0],fobj->qs_obj.gmp_n);
		SIQS(fobj);
		gmp2mp(fobj->qs_obj.gmp_n,&operands[0]);
		print_factors(fobj);
		break;

	case 28:
		//primes
		if (nargs == 2)
		{
			gettimeofday(&tstart, NULL);
			lower = z264(&operands[1]);
			upper = z264(&operands[2]);

			soe_wrapper(spSOEprimes, szSOEp, lower, upper, 1, &n64);

			sp642z(n64,&operands[0]);
			gettimeofday (&tstop, NULL);
			difference = my_difftime (&tstart, &tstop);

			t = ((double)difference->secs + (double)difference->usecs / 1000000);
			free(difference);
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

		lower = z264(&operands[0]);
		upper = z264(&operands[1]);
		free(PRIMES);
		PRIMES = soe_wrapper(spSOEprimes, szSOEp, lower, upper, operands[2].val[0], &NUM_P);
		if (PRIMES != NULL)
		{
			P_MIN = PRIMES[0];
			P_MAX = PRIMES[NUM_P-1];
		}
		sp642z(NUM_P,&operands[0]);

		break;
	case 29:
		//ispow - one argument
		mpz_init(gmpz);

		mp2gmp(&operands[0], gmpz);
		i = mpz_perfect_power_p(gmpz);

		if (i)
			printf("input is a perfect power\n");
		else
			printf("input is not a perfect power\n");

		sp2z(i,&operands[0]);

		mpz_clear(gmpz);
		break;

	case 30:
		//torture - two arguments
		if (nargs != 2)
		{
			printf("wrong number of arguments in torture\n");
			break;
		}
		i = operands[0].val[0];
		k = operands[1].val[0];
		for (j=0; j<i; j++)
		{
			printf("***********RUN %d***********\n",j+1);
			zRand(&mp2,k);
			zCopy(&mp2,&fobj->N);
			factor(fobj);
			zCopy(&fobj->N,&mp2);
			print_factors(fobj);
			clear_factor_list(fobj);
		}

		zClear(&mp2);
		zCopy(&mp2,&operands[0]);
		break;
	case 31:
		//randb - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in randb\n");
			break;
		}
		i = operands[0].val[0];
		zRandb(&operands[0],i);
		break;
	case 32:
		//ecm - two arguments
		zCopy(&operands[0],&fobj->N);
		if (nargs == 2)
		{
			k = operands[1].val[0];
			fobj->ecm_obj.num_curves = k;
			mp2gmp(&operands[0], fobj->ecm_obj.gmp_n);
			ecm_loop(fobj);
			gmp2mp(fobj->ecm_obj.gmp_n, &operands[0]);
		}
		else if (nargs == 1)
		{
			fobj->ecm_obj.num_curves = 1;
			mp2gmp(&operands[1], fobj->ecm_obj.gmp_n);
			ecm_loop(fobj);
			gmp2mp(fobj->ecm_obj.gmp_n, &operands[0]);
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
		zAdd(&operands[0],&operands[1],&mp1);
		zCopy(&mp1,&operands[0]);
		break;

	case 34:
		//subtract or negate
		if (nargs == 1)
		{
			//zNeg(&operands[0]);
			zSub(&zZero,&operands[0],&operands[0]);
		}
		else if (nargs == 2)
		{
			zSub(&operands[0],&operands[1],&mp1);
			zCopy(&mp1,&operands[0]);
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
		zMul(&operands[0],&operands[1],&mp1);
		zCopy(&mp1,&operands[0]);
		break;

	case 36:
		//div
		if (nargs != 2)
		{
			printf("wrong number of arguments in div\n");
			break;
		}
		zDiv(&operands[0],&operands[1],&mp1,&mp2);
		zCopy(&mp1,&operands[0]);
		break;

	case 37:
		//!
		if (nargs != 1)
		{
			printf("wrong number of arguments in factorial\n");
			break;
		}
		zFactorial(operands[0].val[0],&operands[0]);
		break;

	case 38:
		//primorial
		if (nargs != 1)
		{
			printf("wrong number of arguments in primorial\n");
			break;
		}
		zPrimorial(operands[0].val[0],&operands[0]);
		break;

	case 39:
		//==

		break;

	case 40:
		//<<
		if (nargs != 2)
		{
			printf("wrong number of arguments in <<\n");
			break;
		}
		zShiftLeft(&operands[0],&operands[0],operands[1].val[0]);
		break;

	case 41:
		//>>
		if (nargs != 2)
		{
			printf("wrong number of arguments in >>\n");
			break;
		}
		zShiftRight(&operands[0],&operands[0],operands[1].val[0]);
		break;

	case 42:
		//mod
		if (nargs != 2)
		{
			printf("wrong number of arguments in mod\n");
			break;
		}
		zDiv(&operands[0],&operands[1],&mp1,&mp2);
		zCopy(&mp2,&operands[0]);
		break;

	case 43:
		//exp
		if (nargs != 2)
		{
			printf("wrong number of arguments in exp\n");
			break;
		}
		zExp(operands[1].val[0],&operands[0],&mp1);
		zCopy(&mp1,&operands[0]);
		break;

	case 44:
		//test
		//race(operands[0].val[0]);
		//GenTest(operands[0].val[0],operands[1].val[0]);
		//zCopy(&mp1,&operands[0]);
		break;

	case 45:
		//puzzle
		//puzzle1(&operands[0],&operands[1]);
		//wordpuzzle();
		break;

	case 46:
		//sieve
		//guess at the number of expected results
		break;	// not supported in official release

		if (nargs != 2)
		{
			printf("wrong number of arguments in sieve\n");
			break;
		}

		zShortDiv(&operands[1],(fp_digit)((double)(zBits(&operands[1]) - 1) * log(2.0)),&tmp1);
		zShortDiv(&operands[0],(fp_digit)((double)(zBits(&operands[0]) - 1) * log(2.0)),&tmp2);
		zSub(&tmp1,&tmp2,&tmp1);
		printf("expect %u primes\n",(uint32)tmp1.val[0]);
		printf("allocating space for %u entries\n",
			(uint32)(tmp1.val[0] * zBits(&operands[1])/10));
		offsets = (uint32 *)malloc(
			(tmp1.val[0] * zBits(&operands[1])/10) * sizeof(uint32));
		k=0;
		zCopy(&operands[0],&tmp1);
		while(1)
		{
			k += isPrime(&tmp1);
			zShortAdd(&tmp1,1,&tmp1);
			if (zCompare(&tmp1,&operands[1]) > 0)
				break;
		}
		printf("%d actual primes in range\n",k);

		//i = sieve_to_bitdepth(&operands[0],&operands[1],1,offsets);
		k=0;
		for (j=0; j<i; j++)
		{
			zShortAdd(&operands[0],offsets[j],&tmp1);
			if (isPrime(&tmp1) == 0)
				k++;
		}
		printf("%d are not prime\n",k);
		break;

	case 47:
		//algebraic
		break;	// not supported in official release

		if (nargs != 2)
		{
			printf("wrong number of arguments in algebraic\n");
			break;
		}
		printf("base is %s\n",z2decstr(&operands[0],&gstr2));
		printf("exponent is %s\n",z2decstr(&operands[1],&gstr3));
		break;
	case 48:
		//lucas lehmer test
		if (nargs != 1)
		{
			printf("wrong number of arguments in llt\n");
			break;
		}
		k = operands[0].val[0];
		i = llt(k);
		if (i)
			printf("prime!\n");
		else
			printf("composite.\n");
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
		
		//test_dlp_composites();
		//pull_large_primes();
		break;
	case 51: 
		// small factorization test routine
		// 5 inputs

		{
			z *input;
			int test_squfof = 0;
			int test_lehman = 0;
			int test_tinysiqs = 0;
			int test_squfof_tf = 0;
			int force_hard = 1;
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
			int startbits = 90, stopbits = 110, stepbits = 5;
			//int startbits = 40, stopbits = 60, stepbits = 5;

			// get numin
			numin = operands[0].val[0];

			// get method
			if (operands[1].val[0] == 1)
			{ 
				printf("commencing %d tests of squfof method\n",numin);
				test_squfof = 1;
				test_lehman = 0;
				test_tinysiqs = 0;
				test_squfof_tf = 0;
			}
			else if (operands[1].val[0] == 2)
			{ 
				printf("commencing %d tests of lehman method\n",numin);
				test_squfof = 0;
				test_lehman = 1;
				test_tinysiqs = 0;
				test_squfof_tf = 0;
				//init_lehman();
			}
			else if (operands[1].val[0] == 3)
			{ 
				printf("commencing %d tests of squfof + tf method\n",numin);
				test_squfof = 0;
				test_lehman = 0;
				test_squfof_tf = 1;
				test_tinysiqs = 0;
				//init_lehman();
			}
			else if (operands[1].val[0] == 4)
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
				do_tf = operands[2].val[0];

			if (nargs >= 4)
				lehTune = operands[3].val[0];
			
			if (nargs == 5)
				use_hartolf = operands[4].val[0];

			printf("sweeping inputs with bit sizes from %d to %d in steps of %d\n",
				startbits, stopbits, stepbits);

			if (test_lehman)
				printf("lehman parameters:\n  do TF = %d\n  LehTune = %1.0f\n  Use Hart_OLF = %d\n",
				do_tf, lehTune, use_hartolf);

			fobj2 = (fact_obj_t *)malloc(sizeof(fact_obj_t));
			init_factobj(fobj2);
			input = (z *)malloc(numin * sizeof(z));
			for (i=0; i<numin; i++)
				zInit(&input[i]);

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
					zRandb(&mp1, b/2);
					if (force_hard) 
						mp1.val[(b/2)/BITS_PER_DIGIT] |= (fp_digit)(1 << ((b/2) % BITS_PER_DIGIT));
					zRandb(&mp2, b/2);
					if (force_hard) 
						mp2.val[(b/2)/BITS_PER_DIGIT] |= (fp_digit)(1 << ((b/2) % BITS_PER_DIGIT));
					
					zNextPrime(&mp1, &mp1, 1);
					zNextPrime(&mp2, &mp2, 1);

					zMul(&mp1,&mp2,&input[i]);
					nbits = zBits(&input[i]);				
					sum += nbits;
					if (nbits < minsz)
						minsz = nbits;
					if (nbits > maxsz)
						maxsz = nbits;
				}

				gettimeofday (&tstop, NULL);
				difference = my_difftime (&tstart, &tstop);
				t = ((double)difference->secs + (double)difference->usecs / 1000000);
				free(difference);
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
				
					mp2gmp(&input[i], fobj2->qs_obj.gmp_n);
					if (test_squfof)
					{
						mpz_set(gmptmp, fobj2->qs_obj.gmp_n);
						outputs[i] = sp_shanks_loop(gmptmp, fobj);
					}
					else if (test_lehman)
					{
						uint64 N = z264(&input[i]);
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
						//tinySIQS(fobj2); //SIQS(fobj2);

						if (mpz_cmp_ui(fobj2->qs_obj.gmp_n, 1) != 0)
						{						
							gmp_printf("not fully factored!  residue = %Zd\n", fobj2->qs_obj.gmp_n);
							print_factors(fobj2);
						}
						VFLAG = tmpv;
						clear_factor_list(fobj2);
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
				difference = my_difftime (&tstart, &tstop);
				t = ((double)difference->secs + (double)difference->usecs / 1000000);
				free(difference);

				if (test_lehman || test_squfof || test_squfof_tf)
				{
					proper_factors = 0;
					for (i=0; i<numin; i++)
					{
						if ((z264(&input[i]) % outputs[i] == 0) && 
							(outputs[i] > 1) &&
							(z264(&input[i]) != outputs[i]))
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
				zFree(&input[i]);

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
		zCopy(&operands[0],&fobj->N);
		mp2gmp(&operands[0],fobj->qs_obj.gmp_n);
		fobj->logfile = fopen(fobj->flogname,"a");
		if (fobj->logfile == NULL)
			printf("fopen error: %s\n", strerror(errno));
		smallmpqs(fobj);		
		//tinySIQS(fobj);
		if (fobj->logfile != NULL)
			fclose(fobj->logfile);
		gmp2mp(fobj->qs_obj.gmp_n,&operands[0]);
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
			
			range = z264(&operands[2]);
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
			mp2gmp(&operands[0], lowz);
			mp2gmp(&operands[1], highz);
			//printf("test range with %d witnesses\n", operands[3].val[0]);
			primes = sieve_to_depth(sieve_p, num_sp, lowz, highz, 
				0, operands[3].val[0], &num_found);

			free(sieve_p);
			if (!NULL)
				free(primes);

			mpz_clear(lowz);
			mpz_clear(highz);
			sp2z(num_found, &operands[0]);
		}
		
		break;

	case 54:
		printf("siqstune not currently supported\n");
		break;
		//siqstune(operands[0].val[0]);
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
				difference = my_difftime (&tstart, &tstop);
				t = ((double)difference->secs + (double)difference->usecs / 1000000);
				free(difference);
				lower = upper;
				printf("%" PRIu64 ": %" PRIu64 ", elapsed time = %6.4f\n",upper,count,t);
			}
		}

		break;

	case 56: 

		//sieverange - 4 arguments
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
			
			range = z264(&operands[2]);
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
			mp2gmp(&operands[0], lowz);
			mp2gmp(&operands[1], highz);
			primes = sieve_to_depth(sieve_p, num_sp, lowz, highz, operands[3].val[0], 0, &num_found);

			free(sieve_p);
			if (!NULL)
				free(primes);

			mpz_clear(lowz);
			mpz_clear(highz);
			sp2z(num_found, &operands[0]);
		}
		
		break;

	case 57:
		//fermat - three arguments
		if (nargs == 2)
		{
			zCopy(&operands[1],&fobj->N);
			mp2gmp(&operands[1],fobj->div_obj.gmp_n);
			n64 = z264(&operands[2]);
			j = 1;
		}
		else if (nargs == 3)
		{			
			zCopy(&operands[0],&fobj->N);
			mp2gmp(&operands[0],fobj->div_obj.gmp_n);
			n64 = z264(&operands[1]);
			j = operands[2].val[0];
		}
		else
		{
			printf("wrong number of arguments in fermat\n");
			break;
		}
		
		zFermat(n64, j, fobj);
		gmp2mp(fobj->div_obj.gmp_n,&operands[0]);
		print_factors(fobj);
		break;

	case 58:
		//nfs - one argument
		if (nargs != 1)
		{
			printf("wrong number of arguments in nfs\n");
			break;
		}
		zCopy(&operands[0],&fobj->N);
		mp2gmp(&operands[0],fobj->nfs_obj.gmp_n);
		nfs(fobj);
		gmp2mp(fobj->nfs_obj.gmp_n,&operands[0]);
		print_factors(fobj);
		break;

	case 59:
		//tune, no arguments
		factor_tune(fobj);
		break;

	case 60:
		// xor
		{
			mpz_t a,b;

			mpz_init(a);
			mpz_init(b);
			mp2gmp(&operands[0], a);
			mp2gmp(&operands[1], b);
			mpz_xor(a, a, b);
			gmp2mp(a, &operands[0]);
			mpz_clear(a);
			mpz_clear(b);
		}

		break;
	case 61:
		// and
		{
			mpz_t a,b;

			mpz_init(a);
			mpz_init(b);
			mp2gmp(&operands[0], a);
			mp2gmp(&operands[1], b);
			mpz_and(a, a, b);
			gmp2mp(a, &operands[0]);
			mpz_clear(a);
			mpz_clear(b);
		}

		break;
	case 62:
		// or
		{
			mpz_t a,b;

			mpz_init(a);
			mpz_init(b);
			mp2gmp(&operands[0], a);
			mp2gmp(&operands[1], b);
			mpz_ior(a, a, b);
			gmp2mp(a, &operands[0]);
			mpz_clear(a);
			mpz_clear(b);
		}

		break;
	case 63:
		// not
		{
			mpz_t a;

			mpz_init(a);
			mp2gmp(&operands[0], a);
			mpz_com(a, a);
			gmp2mp(a, &operands[0]);
			mpz_clear(a);
		}

		break;

	default:
		printf("unrecognized function code\n");
		zCopy(&zZero,&operands[0]);
		break;
	}

	sFree(&str);
	zFree(&mp1);
	zFree(&mp2);
	zFree(&mp3);
	zFree(&tmp1);
	zFree(&tmp2);
	return 1;
}

/*
int getArgs(str_t *args, str_t *in, int num)
{
	int numparen,argnum;
	char *startptr, *stopptr;
	int i;

	i=0;
	startptr = stopptr = in->s;
	for (argnum=0;argnum<num;argnum++)
	{
		//look for comma
		numparen = 0;
		i=0;
		while (1)
		{
			if (startptr[i] == ',' && numparen == 0)
				break;
			if (startptr[i] == ')')
				numparen--;
			if (startptr[i] == '(')
				numparen++;
			if (startptr[i] == '\0')
				break;
			i++;
		}
		stopptr = startptr + i;

		//get the argument
		args[argnum].nchars = stopptr - startptr;
		strncpy(args[argnum].s,startptr,args[argnum].nchars);
		args[argnum].s[args[argnum].nchars] = '\0';
		startptr = stopptr + 1;	//skip the comma

		//process it
		calc(&args[argnum]);
		i++;
	}

	return argnum + 1;
}
*/

int new_uvar(const char *name, z *data)
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
			zInit(&uvars.vars[i].data);
	}

	strcpy(uvars.vars[uvars.num].name,name);
	zCopy(data,&uvars.vars[uvars.num].data);
	uvars.num++;
	return uvars.num - 1;
}

int set_uvar(const char *name, z *data)
{
	//look for 'name' in the global uvars structure
	//if found, copy in data and return 0
	//else return 1
	int i;

	i = data->val[0];
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
			zCopy(data,&uvars.vars[i].data);
			return 0;
		}
	}
	return 1;
}

int get_uvar(const char *name, z *data)
{
	//look for 'name' in the global uvars structure
	//if found, copy out data and return 0
	//else return 1 if not found
	int i;

	//first look if it is a global constant
	if (strcmp(name,"IBASE") == 0) {
		sp2z(IBASE,data); return 0;}
	else if (strcmp(name,"OBASE") == 0) {
		sp2z(OBASE,data); return 0;}
	else if (strcmp(name,"NUM_WITNESSES") == 0) {
		sp2z(NUM_WITNESSES,data); return 0;}
	else if (strcmp(name,"LOGFLAG") == 0) {
		sp2z(LOGFLAG,data); return 0;}
	else if (strcmp(name,"VFLAG") == 0) {
		sp2z(VFLAG,data); return 0;}
	else if (strcmp(name,"PRIMES_TO_FILE") == 0) {
		sp2z(PRIMES_TO_FILE,data); return 0;}
	else if (strcmp(name,"PRIMES_TO_SCREEN") == 0) {
		sp2z(PRIMES_TO_SCREEN,data); return 0;}

	for (i=0;i<uvars.num;i++)
	{
		if (strcmp(uvars.vars[i].name,name) == 0)
		{
			zCopy(&uvars.vars[i].data,data);
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
			printf("%s      %s\n",uvars.vars[i].name,z2decstr(&uvars.vars[i].data,&gstr1));

		return 2;
	}

	return 1;
}

void free_uvars()
{
	int i;
	for (i=0;i<uvars.alloc;i++)
		zFree(&uvars.vars[i].data);
	free(uvars.vars);
}

