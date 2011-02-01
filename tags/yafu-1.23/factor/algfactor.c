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
#include "arith.h"
#include "calc.h"

/* 
not functional at the moment, but this is where string parsing
and algebraic factorization will eventually occur
*/

int algebraic(str_t *str)
{
	//scan the string for expressions that can be algebraicly reduced
	//right now this only handles expressions of the form
	//num^num ... etc.  for example, no parens or nested expressions in the
	//exponent or base

	//to do -> general algebraic polynomial factorization
/*
	char ch;
	char rexpo_str[80];
	char expo_str[80];
	int i, el_type, num_expo, numchars, k, j;
	const int MAX_EXPO = 10;
	z *exponents;

	exponents = (z *)malloc(MAX_EXPO * sizeof(z));
	for (i=0;i<MAX_EXPO;i++)
		zInit(&exponents[i]);

	//read until we find a number for the base
	for (i=0;i<strlen(str->s);i++)
	{
		if (get_el_type2(str->s[i]) != NUM)
			continue;
		break;
	}
	k=i;
	j=0;
	//how many digits in base
	for (;i<strlen(str->s);i++)
	{
		if (get_el_type2(str->s[i]) == NUM)
			j++;
		else
			break;
	}

	memcpy(expo_str,str->s+k,j*sizeof(char));
	expo_str[j] = '\0';
	str2hexz(expo_str,&exponents[0]);

	printf("base is %s\n",z2decstr(&exponents[0],&gstr1));

	//read until we find a ^
	for (;i<strlen(str->s);i++)
	{
		if (str->s[i] == '^')
			break;
	}

	j=0;
	i++;
	k=i;
	//how many digits in exponent
	for (;i<strlen(str->s);i++)
	{
		if (get_el_type2(str->s[i]) == NUM)
			j++;
		else
			break;
	}

	memcpy(expo_str,str->s+k,j*sizeof(char));
	expo_str[j] = '\0';
	str2hexz(expo_str,&exponents[1]);

	printf("exponent is %s\n",z2decstr(&exponents[1],&gstr1));


free:

	for (i=0;i<MAX_EXPO;i++)
		zFree(&exponents[i]);
	free(exponents);
*/
	return 1;
}
