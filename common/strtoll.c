/*--------------------------------------------------------------------
Copyright (c) 1988 The Regents of the University of California.
Copyright (c) 1994 Sun Microsystems, Inc.

Modified by Jason Papadopoulos for use within the Msieve library

$Id: strtoll.c 185 2010-01-31 19:03:51Z jaysonking $
--------------------------------------------------------------------*/

#if defined(_MSC_VER)

#include <util.h>

/*
 * The table below is used to convert from ASCII digits to a
 * numerical equivalent.  It maps from '0' through 'z' to integers
 * (100 for non-digit characters).
 */

static char cvt[] = {
0, 1, 2, 3, 4, 5, 6, 7, 8, 9,		/* '0' - '9' */
100, 100, 100, 100, 100, 100, 100,		/* punctuation */
10, 11, 12, 13, 14, 15, 16, 17, 18, 19,	/* 'A' - 'Z' */
20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
30, 31, 32, 33, 34, 35,
100, 100, 100, 100, 100, 100,		/* punctuation */
10, 11, 12, 13, 14, 15, 16, 17, 18, 19,	/* 'a' - 'z' */
20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
30, 31, 32, 33, 34, 35};

/*-----------------------------------------------------------------------*/
uint64
strtoull(const char *string, char **endptr, int base) {

	const char *p;
	uint64 result = 0;
	unsigned digit;
	uint64 shifted;
	int anydigits = 0, negative = 0;

	/* skip any leading blanks */

	p = string;
	while (isspace(*p))
		p++;

	/* check for a sign */

	if (*p == '-') {
		p++;
		negative = 1;
	} 
    else if (*p == '+') {
		p++;
	}

	/* if no base was provided, pick one from the 
	   leading characters of the string  */
	
	if (base == 0) {
		if (*p == '0') {
			p++;
			if (*p == 'x' || *p == 'X') {
				p++;
				base = 16;
			} 
			else {
				/* must set anyDigits here, otherwise "0" 
		     		   produces a "no digits" error */ 
				anydigits = 1;
				base = 8;
			}
		} 
		else {
			base = 10;
		}
	} 
	else if (base == 16) {
		/* skip a leading "0x" from hex numbers */
		if ((p[0] == '0') && (p[1] == 'x' || *p == 'X'))
			p += 2;
	}

	if (base == 8) {
		while (1) {
			digit = *p - '0';
			if (digit > 7)
				break;

			shifted = result << 3;
			if ((shifted >> 3) != result)
				goto overflow;

			result = shifted + digit;
			if (result < shifted)
				goto overflow;

			anydigits = 1;
			p++;
		}
	} 
	else if (base == 10) {
		while (1) {
			digit = *p - '0';
			if (digit > 9)
				break;

			shifted = 10 * result;
			if ((shifted / 10) != result)
				goto overflow;

			result = shifted + digit;
			if (result < shifted)
				goto overflow;

			anydigits = 1;
			p++;
		}
	} 
	else if (base == 16) {
		while (1) {
			digit = *p - '0';
			if (digit > ('z' - '0'))
				break;

			digit = cvt[digit];
			if (digit > 15)
				break;

			shifted = result << 4;
			if ((shifted >> 4) != result)
				goto overflow;

			result = shifted + digit;
			if (result < shifted)
				goto overflow;

			anydigits = 1;
			p++;
		}
	} 
	else if (base >= 2 && base <= 36) {
		while (1) {
			digit = *p - '0';
			if (digit > ('z' - '0'))
				break;

			digit = cvt[digit];
			if (digit >= (unsigned) base)
				break;

			shifted = result * base;
			if ((shifted / base) != result)
				goto overflow;

			result = shifted + digit;
			if ( result < shifted )
				goto overflow;

			anydigits = 1;
			p++;
		}
	}

	if (negative)
		result = -result;
	if (!anydigits)
		p = string;
	if (endptr != 0)
		*endptr = (char *) p;

	return result;

	/* on overflow generate the right output */

overflow:
	if (endptr != 0) {
		while (1) {
			digit = *p - '0';
			if (digit > ('z' - '0'))
				break;

			digit = cvt[digit];
			if (digit >= (unsigned) base)
				break;

			p++;
		}
		*endptr = (char *) p;
	}
	return (uint64)(-1);
}

/*-----------------------------------------------------------------------*/
int64
strtoll(const char *string, char **endptr, int base) {

	const char *p;
	int64 result = 0;
	uint64 uresult;

	/* skip any leading blanks */

	p = string;
	while (isspace(*p))
		p++;

	/* check for a sign */

	errno = 0;
	if (*p == '-') {
		p++;
		uresult = strtoull(p, endptr, base);
		if (uresult >= ((uint64)1 << 62))
			return (uint64)1 << 63;
		result = -uresult;
	} 
	else {
		if (*p == '+')
			p++;

		uresult = strtoull(p, endptr, base);
		if (uresult >= ((uint64)1 << 62))
			return ((uint64)1 << 62) - 1;
		result = uresult;
	}

	if ((result == 0) && (endptr != 0) && (*endptr == p))
		*endptr = (char *)string;
	return result;
}

#endif /* _MSC_VER */
