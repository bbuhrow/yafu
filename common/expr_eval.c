/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Brian Gladman. You may use it for any purpose, free of charge, without 
having to notify anyone. I disclaim any responsibility for any errors.

Modified for use within the msieve library by Jason Papadopoulos

$Id: expr_eval.c 937 2013-08-08 00:19:28Z jasonp_sf $
--------------------------------------------------------------------*/

#include <ms_common.h>

#define LEFT_BRKT	'('
#define RIGHT_BRKT	')'

#define RETURN_SUCCESS		0
#define OUT_OF_MEMORY		-1
#define STACK_UNDERFLOW		-2
#define STACK_OVERFLOW		-3
#define ILLEGAL_INPUT_CHARACTER	-4
#define UNKNOWN_OPERATOR	-5
#define MISSING_RIGHT_BRACKET	-6
#define MISSING_LEFT_BRACKET	-7
#define EXPRESSION_ERROR	-8
#define EVALUATION_ERROR	-9
#define INTEGER_OVERFLOW	-10
#define SIGN_ERROR		-11

#define STACK_SIZE 100

typedef struct {
	void *base[STACK_SIZE];
	int num_used;
} eval_stack_t;

enum char_type { is_space, is_digit, is_operator, is_token, is_invalid };

/* arithmetic operators */
static const char *operator_list = "+-%*/^()";   

/* precedence */
static const int precedence[9] = { 1, 1, 2, 3, 3, 4, 5, 5 };

/*---------------------------------------------------------------------*/
static int isoperator(int c) {

	if (strchr(operator_list, c))
		return 1;
	return 0;
}

/*---------------------------------------------------------------------*/
static int find_precedence(int op_symbol) {

	char *p = strchr(operator_list, op_symbol);
	if (p) {
		return precedence[p - operator_list];
	}
	return 0;
}

/*---------------------------------------------------------------------*/
static void stack_init(eval_stack_t *stack) {

	stack->num_used = 0;
}

/*---------------------------------------------------------------------*/
static void stack_free(eval_stack_t *stack) {

	int i;
	for (i = 0; i < stack->num_used; i++) {
		free(stack->base[i]);
	}
	stack->num_used = 0;
}

/*---------------------------------------------------------------------*/
static int stack_push(eval_stack_t *stack, 
		const void *symbol, int symbol_len) {

	char *curr_entry;

	if(stack->num_used >= STACK_SIZE)
		return STACK_OVERFLOW;
		
	curr_entry = (char *)xmalloc((size_t)(symbol_len + 1));
	if (curr_entry == NULL)
		return OUT_OF_MEMORY;
	memcpy(curr_entry, symbol, (size_t)symbol_len);
	curr_entry[symbol_len] = 0;
	stack->base[stack->num_used++] = curr_entry;
	return RETURN_SUCCESS;
}

/*---------------------------------------------------------------------*/
static int stack_pop(eval_stack_t *stack) {

	if (stack->num_used <= 0)
		return STACK_UNDERFLOW;

	free(stack->base[--stack->num_used]);
	return RETURN_SUCCESS;
}

/*---------------------------------------------------------------------*/
static void *eval_stack_top(eval_stack_t *stack, int pos) {

	if (stack->num_used)
		return stack->base[stack->num_used - pos - 1];
	return NULL;
}

/*---------------------------------------------------------------------*/
static void *stack_pos(eval_stack_t *stack, int pos)
{
	if (pos < stack->num_used)
		return stack->base[pos];
	return NULL;
}

/*---------------------------------------------------------------------*/
static int input_to_tokens(const char* input, eval_stack_t *tokens) {

	enum char_type curr_type = is_invalid; 
	enum char_type last_type;
	int len = (int)strlen(input);
	int i, pos, status;

	for (i = pos = 0; i < len; i++) {
		char c = input[i];
		last_type = curr_type;
		if (isdigit(c))
			curr_type = is_digit;
		else if (isoperator(c))
			curr_type = is_operator;
		else if (isspace(c))
			curr_type = is_space;
		else
			return ILLEGAL_INPUT_CHARACTER;

		if (last_type == is_operator || 
		    (i > pos && curr_type != last_type)) {
			status = stack_push(tokens, input + pos, i - pos);
			if (status < 0)
				return status;
			pos = i;
		}
		if (curr_type == is_space)
			pos++;
	}
	if (i > pos) {
		status = stack_push(tokens, input + pos, i - pos);
		if (status < 0)
			return status;
	}

	return tokens->num_used;
}

/*---------------------------------------------------------------------*/
static int has_higher_precedence(int first, int second) {

	int prec1 = find_precedence(first); 
	int prec2 = find_precedence(second);

	if (prec1 == 0 || prec2 == 0)
		return UNKNOWN_OPERATOR;
	return (prec1 > prec2);
}

/*---------------------------------------------------------------------*/
static int infix_to_postfix(const char *str, char *output) {

	int num_tokens, i;
	int status = 0;
	eval_stack_t tokens;
	eval_stack_t symbol_stack;

	stack_init(&tokens);
	stack_init(&symbol_stack);
	num_tokens = input_to_tokens(str, &tokens);
	if (num_tokens < 0)
		return num_tokens;
	else if (num_tokens == 0)
		return RETURN_SUCCESS;

	for (i = 0; i < num_tokens; i++) {

		char sym = *(char*)stack_pos(&tokens, i);
		if (!isoperator(sym)) {
			output += sprintf(output, "%s ", 
					(char *)stack_pos(&tokens, i));
			continue;
		}

		if (symbol_stack.num_used && sym == RIGHT_BRKT) {
			/* drop a level */
			while (symbol_stack.num_used && 
			       strcmp(eval_stack_top(&symbol_stack, 0), "(")) {

				output += sprintf(output, "%s ",
					      (char *)eval_stack_top(&symbol_stack, 0));
				status = stack_pop(&symbol_stack);
				if (status < 0)
					return status;
			}

			if (symbol_stack.num_used == 0) {
				return MISSING_LEFT_BRACKET;
			}
			else {
				status = stack_pop(&symbol_stack);
				if (status < 0)
					return status;
			}
		}
		else if (symbol_stack.num_used == 0 || 
			 sym == LEFT_BRKT || 
			 (status = has_higher_precedence(sym, 
					*(char *)eval_stack_top(&symbol_stack, 0)))) {
			if (status < 0)
				return status;

			status = stack_push(&symbol_stack, 
					stack_pos(&tokens, i),
					(int)strlen(stack_pos(&tokens, i)));
			if (status < 0)
				return status;
		}
		else {
			while(!(status = has_higher_precedence(sym, 
				*(char*)eval_stack_top(&symbol_stack, 0))) ||
			      strcmp(stack_pos(&tokens, i), 
				     eval_stack_top(&symbol_stack, 0)) == 0) {

				if (strcmp(eval_stack_top(&symbol_stack, 0), 
						"(" ) == 0)
					break;

				output += sprintf(output, "%s ", (char *)
					      eval_stack_top(&symbol_stack, 0));

				status = stack_pop(&symbol_stack);
				if (status < 0)
					return status;

				if(symbol_stack.num_used == 0)
					break;
			}
			if (status < 0)
				return status;

			status = stack_push(&symbol_stack, 
					stack_pos(&tokens, i),
					(int)strlen(stack_pos(&tokens, i)));
			if (status < 0)
				return status;
		}
	}

	while (symbol_stack.num_used) {
		if (strcmp(eval_stack_top(&symbol_stack, 0), "(") == 0)
			return MISSING_RIGHT_BRACKET;

		output += sprintf(output, "%s ", 
				(char *)eval_stack_top(&symbol_stack, 0));

		status = stack_pop(&symbol_stack);
		if (status < 0)
			return status;
	}

	stack_free(&symbol_stack);
	stack_free(&tokens);
	return RETURN_SUCCESS;
}

/*---------------------------------------------------------------------*/
static int do_binary_op(mp_t *a, mp_t *b, int op, mp_t *r)
{
	switch(op) {
	case '+':
		mp_add(a, b, r); 
		break;

	case '-':
		if (mp_cmp(a, b) < 0)
			return SIGN_ERROR;
		mp_sub(a, b, r); 
		break;

	case '*':
		if (mp_bits(a) + mp_bits(b) >= 32 * MAX_MP_WORDS - 1)
			return INTEGER_OVERFLOW;
		mp_mul(a, b, r); 
		break;

	case '/':
		mp_div(a, b, r); 
		break;

	case '%':
		mp_mod(a, b, r); 
		break;

	case '^':
		if (b->nwords > 1)
			return INTEGER_OVERFLOW;
		if (a->nwords && b->val[0] * mp_log(a) / M_LN2 >=
					32.0 * MAX_MP_WORDS - 2) {
			return INTEGER_OVERFLOW;
		}
		mp_pow(a, b, r); 
		break;
	}
	return RETURN_SUCCESS;
}

/*---------------------------------------------------------------------*/
static int mp_evaluate(char *str, mp_t *res) {	

	int status;
	mp_t *a, *b, r;
	int num_tokens, i;
	eval_stack_t tokens;
	eval_stack_t int_stack;

	stack_init(&tokens);
	stack_init(&int_stack);
	num_tokens = input_to_tokens(str, &tokens);
	if (num_tokens < 0)
		return num_tokens;
	else if (num_tokens == 0)
		return RETURN_SUCCESS;

	for(i = 0; i < num_tokens; i++) {

		char sym = *(char*)stack_pos(&tokens, i);

		if (!isoperator(sym)) {
			mp_str2mp(stack_pos(&tokens, i), &r, (uint32)10);
			status = stack_push(&int_stack, &r, (int)sizeof(mp_t));
			if (status < 0)
				return status;
		}
		else if (int_stack.num_used > 1) {
			a = (mp_t *)eval_stack_top(&int_stack, 0);
			b = (mp_t *)eval_stack_top(&int_stack, 1);
			if ((status = do_binary_op(b, a, sym, &r)) < 0)
				return status;
			if ((status = stack_pop(&int_stack)) < 0)
				return status;
			if ((status = stack_pop(&int_stack)) < 0)
				return status;
			if ((status = stack_push(&int_stack, &r, 
						(int)sizeof(mp_t))) < 0)
				return status;
		}
		else {
			return EXPRESSION_ERROR;
		}
	}
	mp_copy((mp_t *)eval_stack_top(&int_stack, 0), res);
	if(int_stack.num_used > 1)
		return EXPRESSION_ERROR;
	stack_free(&tokens);
	stack_free(&int_stack);
	return RETURN_SUCCESS;
}

/*---------------------------------------------------------------------*/
int32 evaluate_expression(char *expr, mp_t *res) {

	int status;
	char postfix[BIGNUM_BUF_SIZE];
	char *tmp;

	/* fast path: if expr is an ascii integer, convert
	   it immediately */
	
	mp_clear(res);
	tmp = expr;
	if (tmp[0] == '0' && tolower(tmp[1]) == 'x')
		tmp += 2;
	else if (tmp[0] == '0')
		tmp++;

	while (*tmp != 0) {
		if (!isxdigit(*tmp))
			break;
		tmp++;
	}
	if (*tmp == 0) {
		if (tmp - expr >= 311) {
			printf("input integers must be under 311 digits\n");
			return INTEGER_OVERFLOW;
		}
		mp_str2mp(expr, res, 0);
		return RETURN_SUCCESS;
	}

	status = infix_to_postfix(expr, postfix);
	if (status < 0)
		return status;

	status = mp_evaluate(postfix, res);
	if (status < 0)
		return status;
	return RETURN_SUCCESS;
}
