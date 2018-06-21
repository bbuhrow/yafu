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

void mp2gmp(z *src, mpz_t dest) {

	mpz_import(dest, (size_t)(abs(src->size)), -1, sizeof(fp_digit), 
			0, (size_t)0, src->val);

	if (src->size < 0)
		mpz_neg(dest, dest);
}

void gmp2mp(mpz_t src, z *dest) {

	size_t count;

	dest->size = 1;
	dest->val[0] = 0;
	mpz_export(dest->val, &count, -1, sizeof(fp_digit),
			0, (size_t)0, src);
	dest->size = count;
}

char * mpz_conv2str(char **in, int base, mpz_t n)
{
	// a wrapper for mpz_get_str that will grow the
	// input string if necessary to hold n, when 'in' is
	// not null (in which case mpz_get_str will create
	// a string of the appropriate size)

	if (*in == NULL)
	{
		return mpz_get_str(*in, base, n);
	}
	else
	{
		int nchars = mpz_sizeinbase(n, base) + 2;
		if (nchars >= GSTR_MAXSIZE)
			*in = (char *)realloc(*in, nchars * sizeof(char));
		return mpz_get_str(*in, base, n);	
	}
}

uint64 mpz_get_64(mpz_t src)
{

	uint64 out = mpz_getlimbn(src, 0);
#if GMP_LIMB_BITS == 32
	if (mpz_size(src) >= 2)
		out |= ((uint64)mpz_getlimbn(src, 1) << 32ULL);
#endif

	return out;

}

void mpz_set_64(mpz_t dest, uint64 src)
{
    //printf("%d %d %d\n", GMP_LIMB_BITS, sizeof(mp_limb_t), mp_bits_per_limb);
#if GMP_LIMB_BITS == 64
	dest->_mp_d[0] = src;
	dest->_mp_size = (src ? 1 : 0);
#else
	/* mpz_import is terribly slow */
	mpz_set_ui(dest, (uint32)(src >> 32));
	mpz_mul_2exp(dest, dest, 32);
	mpz_add_ui(dest, dest, (uint32)src);
#endif

}

void mpz_to_z32(mpz_t src, z32 *dest)
{
	int i;

#if GMP_LIMB_BITS == 32
	for (i=0; i < mpz_size(src); i++)
		dest->val[i] = mpz_getlimbn(src, i);

	dest->size = mpz_size(src);
#else
	int j = 0;
	for (i=0; i < mpz_size(src); i++)
	{
		uint64 tmp = mpz_getlimbn(src, i);
		dest->val[j] = (uint32)tmp;
		dest->val[j+1] = (uint32)(tmp >> 32);
		j += 2;
	}
	if (dest->val[j-1] == 0)
		dest->size = j-1;
	else
		dest->size = j;
#endif
	return;
}

void z32_to_mpz(z32 *src, mpz_t dest)
{
	int i;
#if GMP_LIMB_BITS == 32
	for (i=0; i < src->size; i++)
		dest->_mp_d[i] = src->val[i];
	dest->_mp_size = src->size;
#else
	int j=0;
	if (src->size & 1)
		src->val[src->size] = 0;

	for (i=0; i < src->size; i+=2)
		dest->_mp_d[j++] = (uint64)src->val[i] | ((uint64)src->val[i+1] << 32);
	dest->_mp_size = j;

#endif
	return;
}


void zCopy(z *src, z *dest)
{
	//physically copy the digits of u into the digits of v
	int i,su = abs(src->size);
	
	if (src == dest)
		return;

	if (dest->alloc < su)
		zGrow(dest,su);

	//memcpy(dest->val,src->val,su * sizeof(fp_digit));
    for (i = 0; i < su; i++)
    {
        dest->val[i] = src->val[i];
    }

	dest->size = src->size;
	dest->type = src->type;
	return;
}

void zCopy32(z32 *src, z32 *dest)
{
	//physically copy the digits of u into the digits of v
	int i,su = abs(src->size);
	
	if (src == dest)
		return;

	if (dest->alloc < su)
		zGrow32(dest,su);

	//memcpy(dest->val,src->val,su * sizeof(uint32));
    for (i = 0; i < su; i++)
    {
        dest->val[i] = src->val[i];
    }

	dest->size = src->size;
	dest->type = src->type;
	return;
}

void zInit(z *num)
{
	num->val = (fp_digit *)calloc(MAX_DIGITS, sizeof(fp_digit));
	if (num->val == NULL)
	{
		printf("couldn't allocate bignum in zInit\n");
		exit(1);
	}
	num->size = 1;
	num->alloc = MAX_DIGITS;
	num->type = UNKNOWN;
	return;
}

void zInit32(z32 *num)
{
	num->val = (uint32 *)calloc(MAX_DIGITS, sizeof(uint32));
	if (num->val == NULL)
	{
		printf("couldn't allocate bignum in zInit\n");
		exit(1);
	}
	num->size = 1;
	num->alloc = MAX_DIGITS;
	num->type = UNKNOWN;
	return;
}

void zGrow(z *num, int newsz)
{
	//printf("growing\n");
	num->val = (fp_digit *)realloc(num->val, (abs(newsz) + 2) * sizeof(fp_digit));
	if (num->val == NULL)
	{
		printf("realloc failure in grow, couldn't allocate %d bytes\n",
			(int)((abs(newsz)+2)*sizeof(fp_digit)));
		exit(1);
	}
	num->alloc = abs(newsz)+2;
	return;
}

void zGrow32(z32 *num, int newsz)
{
	//printf("growing\n");
	num->val = (uint32 *)realloc(num->val, (abs(newsz) + 2) * sizeof(uint32));
	if (num->val == NULL)
	{
		printf("realloc failure in grow, couldn't allocate %d bytes\n",
			(int)((abs(newsz)+2)*sizeof(uint32)));
		exit(1);
	}
	num->alloc = abs(newsz)+2;
	return;
}

void zFree(z *num)
{
	free(num->val);
	return;
}

void zFree32(z32 *num)
{
	free(num->val);
	return;
}

void zClear(z *num)
{
	if (num->val != NULL)
		memset(num->val, 0, num->alloc * sizeof(fp_digit));
	num->size = 1;
	num->type = UNKNOWN;
	return;
}

void zClear32(z32 *num)
{
	memset(num->val, 0, num->alloc * sizeof(uint32));
	num->size = 1;
	num->type = UNKNOWN;
	return;
}

char *z2decstr(z *n, str_t *s)
{
	//pass in a pointer to a string.  if necessary, this routine will 
	//reallocate space for the string to accomodate its size.  If this happens
	//the pointer to the string's (likely) new location is automatically
	//updated and returned.
	z a;
	int i,sza;
	char *tmp;

	//for really long inputs, a significant amount of time is spent here.
	//for instance, in computing 10000!, 0.047sec is spent on actually 
	//computing the factorial, while ~.5 sec is needed for the Hex2Dec conversion
	//and ~.8 sec is required to print it to a string.
	//maybe try to unroll the loop a bit?

	strcpy(s->s,"");
	s->nchars = 1;
	zInit(&a);

	//printf("starting hex 2 dec conversion\n");
	zHex2Dec(n,&a);

	sza = abs(a.size);

	if (s->alloc < DEC_DIGIT_PER_WORD*sza + 2)
	{
		s->s = (char *)realloc(s->s,(DEC_DIGIT_PER_WORD*sza + 10)*sizeof(char));
		s->alloc = (DEC_DIGIT_PER_WORD*sza + 10);
	}

	tmp = (char *)malloc(30);

	//print negative sign, if necessary
	if (n->size < 0)
	{
		 sprintf(s->s,"-");
		 s->nchars++;
	}

	//print first word
#if BITS_PER_DIGIT == 32
		sprintf(s->s,"%s%u",s->s,(uint32)a.val[sza - 1]);
		s->nchars += ndigits_1(a.val[sza-1]) - 1;

		//print the rest
		for (i=sza - 2; i>=0; i--)
		{
            int j;
			//sprintf(s->s,"%s%09u",s->s,a.val[i]);
			//s->nchars += 9;
			sprintf(tmp,"%09u",(uint32)a.val[i]);
			//memcpy(s->s + s->nchars, tmp, 9);
            for (j=0; j<9; j++)
                s->s[j] = tmp[j];
			s->nchars += 9;
		}
#else
		sprintf(s->s,"%s%" PRIu64,s->s,a.val[sza - 1]);
		s->nchars += ndigits_1(a.val[sza-1]) - 1;

		//print the rest
		for (i=sza - 2; i>=0; i--)
		{
            int j;
			//sprintf(s->s,"%s%09u",s->s,a.val[i]);
			//s->nchars += 9;
			sprintf(tmp,"%019" PRIu64,a.val[i]);
			//memcpy(s->s + s->nchars, tmp, 19);
            for (j = 0; j<19; j++)
                s->s[j] = tmp[j];
			s->nchars += 19;
		}
#endif

	s->s[s->nchars] = '\0';
	s->nchars++;

	zFree(&a);
	free(tmp);
	return s->s;
}

char *z2hexstr(z *n, str_t *s)
{
	//input bignum n, and an already allocated str_t, s
	//convert n into s, in hex

	int i,szn = abs(n->size);
	char *tmp;

	strcpy(s->s,"");
	s->nchars = 0;

	if (s->alloc < HEX_DIGIT_PER_WORD*szn + 3)
	{
		s->s = (char *)realloc(s->s,(HEX_DIGIT_PER_WORD*szn + 10)*sizeof(char));
		s->alloc = (HEX_DIGIT_PER_WORD*szn + 10);
	}

	tmp = (char *)malloc(30);

	if (n->size < 0)
	{
		 sprintf(s->s,"-");
		 s->nchars++;
	}
	sprintf(s->s,"%s0x",s->s);
	s->nchars += 2;

#if BITS_PER_DIGIT == 32
		//print first word
		sprintf(s->s,"%s%x",s->s,(uint32)n->val[szn - 1]);
		s->nchars = strlen(s->s);

		//print the rest
		for (i=szn - 2; i>=0; i--)
		{
            int j;
			//sprintf(s->s,"%s%09u",s->s,a.val[i]);
			//s->nchars += 9;
			sprintf(tmp,"%08x",(uint32)n->val[i]);
			//memcpy(s->s + s->nchars, tmp, 8);
            for (j = 0; j<8; j++)
                s->s[j] = tmp[j];
			s->nchars += 8;
		}
#else
		//print first word
		sprintf(s->s,"%s%" PRIx64,s->s,n->val[szn - 1]);
		s->nchars = strlen(s->s);

		//print the rest
		for (i=szn - 2; i>=0; i--)
		{
            int j;
			sprintf(tmp,"%016" PRIx64,n->val[i]);
			//memcpy(s->s + s->nchars, tmp, 16);
            for (j = 0; j<16; j++)
                s->s[j] = tmp[j];
			s->nchars += 16;
		}
#endif

	s->s[s->nchars] = '\0';
	s->nchars++;

	free(tmp);
	return s->s;
}

void str2hexz(char in[], z *u)
{
	mpz_t t;
	mpz_init(t);
	if (mpz_set_str(t, in, 0) < 0)
	{
		// not valid - couldn't determine base?  Try the ones we support:
		if (mpz_set_str(t, in, 2) < 0)
		{
			if (mpz_set_str(t, in, 8) < 0)
			{
				if (mpz_set_str(t, in, 10) < 0)
				{
					if (mpz_set_str(t, in, 16) < 0)
					{
						printf("str2hexz couldn't convert input string %s\n", in);
						exit(1);
					}
				}
			}
		}
	}
	gmp2mp(t, u);
	mpz_clear(t);
	return;
}

void mp_t2z(mp_t *src, z *dest)
{
#if BITS_PER_DIGIT == 64
	z32 tmp;
	int i, sz,j;

	zInit32(&tmp);	
	
	if (src->nwords & 0x1)
	{
		sz = src->nwords + 1;
		src->val[src->nwords]=0;
	}
	else
		sz = src->nwords;

	j = 0 ;
	for (i=0; i<sz; i+=2)
		dest->val[j++] = (uint64)src->val[i] | ((uint64)src->val[i+1] << 32);

	dest->size = j;
	dest->type = UNKNOWN;

	zFree32(&tmp);
#else
	int i;

	dest->size = src->nwords;
	for (i=0; i<src->nwords; i++)
		dest->val[i] = src->val[i];

#endif

	return;
}

void dbl2z(double n, z *a)
{
	char s[GSTR_MAXSIZE];

	sprintf(s,"%0.0f",n);
	str2hexz(s,a);

	return;
}

double z2dbl(z *a)
{
	double d;
	d = a->val[0];
	d = d + MAX_DIGIT*(double)a->val[1];
	if (a->size < 0)
		d *= -1;
	return d;
}

void sp2z(fp_digit sp, z *mp)
{
	mp->size = 1;
	//printf("mp->val[0] = %llu\n",mp->val[0]);
	mp->val[0] = sp;
	return;
}

void sp642z(uint64 sp, z *mp)
{
#if BITS_PER_DIGIT == 32
		mp->size = 1;
		mp->val[0] = (uint32)sp;
		mp->val[1] = (uint32)(sp >> 32);
		if (mp->val[1] != 0)
			mp->size++;

#else
		mp->size = 1;
		mp->val[0] = sp;
#endif

	return;
}

int isFive(z *n)
{
	return (abs(n->size) == 1 && n->val[0] == 5);
}

int isZero(z *n)
{
	return (abs(n->size) == 1 && n->val[0] == 0);
}

int isOne(z *n)
{
	return (abs(n->size) == 1 && n->val[0] == 1);
}

uint64 z264(z *n)
{
	//assumes n is only 2 or less digits long
	//sign information is lost.
	uint64 out = (uint64)n->val[0];

#if BITS_PER_DIGIT == 32
	if (n->size > 1)
		return (out | ((uint64)n->val[1] << 32));
	else
		return out;
#else
	return out;
#endif
}

int ndigits_1(fp_digit n)
{
	int i=0;
	while (n != 0)
	{
		n /= 10;
		i++;
	}
	if (i==0)
		i++;
	return i;
}

int ndigits(z *n)
{
	int i=0;
	z nn,tmp;

	//can get within one digit using zBits and logs, which would
	//be tons faster.  Any way to 'correct' the +/- 1 error?
	zInit(&nn);
	zInit(&tmp);
	zCopy(n,&tmp);
	while (tmp.size > 1)
	{
		zCopy(&tmp,&nn);
		zShortDiv(&nn,MAX_DEC_WORD,&tmp);
		i += DEC_DIGIT_PER_WORD;
	}
	i += ndigits_1(tmp.val[0]);
	zFree(&nn);
	zFree(&tmp);
	return i;
}

int zCompare(z *u, z *v)
{
	//return 1 if u > v, -1 if u < v, 0 if equal
	int i,j,su,sv;

	i = u->size < 0;
	j = v->size < 0;

	su = abs(u->size);
	sv = abs(v->size);
	//su = i ? -1*u->size : u->size;
	//sv = j ? -1*v->size : v->size;
	
	if (i > j) 
	{
		//v pos, u neg
		//make sure both are not zero
		if (u->val[0] == 0 && su == 1 && v->val[0] == 0 && sv == 1)
			return 0;
		else
			return -1;	
	}
	if (j > i) 
	{
		//u pos, v neg
		//make sure both are not zero
		if (u->val[0] == 0 && su == 1 && v->val[0] == 0 && sv == 1)
			return 0;
		else
			return 1;	
	}	

	//check obvious
	if (j)
	{	//both are negative
		if (su > sv) return -1;
		if (su < sv) return 1;
	}
	else
	{	//both are positive
		if (su > sv) return 1;
		if (su < sv) return -1;
	}

	//if the numbers are both negative, then we'll need to switch the return value
	for (i = su - 1; i>=0; --i)
	{
		if (u->val[i] > v->val[i]) 
			return (1 - 2*j);
		if (u->val[i] < v->val[i])
			return (-1 + 2*j);
	}

	//equal if got to here
	return 0;
}

int zCompare32(z32 *u, z32 *v)
{
	//return 1 if u > v, -1 if u < v, 0 if equal
	int i,j,su,sv;

	i = u->size < 0;
	j = v->size < 0;

	su = abs(u->size);
	sv = abs(v->size);
	
	if (i > j) 
	{
		//v pos, u neg
		//make sure both are not zero
		if (u->val[0] == 0 && su == 1 && v->val[0] == 0 && sv == 1)
			return 0;
		else
			return -1;	
	}
	if (j > i) 
	{
		//u pos, v neg
		//make sure both are not zero
		if (u->val[0] == 0 && su == 1 && v->val[0] == 0 && sv == 1)
			return 0;
		else
			return 1;	
	}	

	//check obvious
	if (j)
	{	//both are negative
		if (su > sv) return -1;
		if (su < sv) return 1;
	}
	else
	{	//both are positive
		if (su > sv) return 1;
		if (su < sv) return -1;
	}

	//if the numbers are both negative, then we'll need to switch the return value
	for (i = su - 1; i>=0; --i)
	{
		if (u->val[i] > v->val[i]) 
			return (1 - 2*j);
		if (u->val[i] < v->val[i])
			return (-1 + 2*j);
	}

	//equal if got to here
	return 0;
}

void zDec2Hex(z *u, z *v)
{
	//convert u[] in dec to v[] in hex by multiplying the ith digit by (1e9)*i
	//and adding to the previous digits

	z a,b,vv;
	int i,su = abs(u->size);

	zInit(&a);
	zInit(&b);
	zInit(&vv);

	if (v->alloc < su)
		zGrow(v,su);

	if (a.alloc < su)
	{
		zGrow(&a,su);
		zClear(&a);
	}

	if (b.alloc < su)
	{
		zGrow(&b,su);
		zClear(&b);
	}

	if (vv.alloc < su)
	{
		zGrow(&vv,su);
		zClear(&vv);
	}
	vv.size = su;

	//a holds the value of (1e9)*i
	a.size = 1;
	a.val[0] = 1;
	for (i=0;i<su;i++)
	{
		zShortMul(&a,u->val[i],&b);
		zAdd(&vv,&b,&vv);
		zShortMul(&a,MAX_DEC_WORD,&a);
	}

	//v may have unused high order limbs
	for (i=su-1;i>=0;i--)
	{
		if (vv.val[i] != 0)
			break;
	}
	vv.size = i+1;

	if (u->size < 0)
		vv.size *= -1;

	if (vv.size == 0)
		vv.size = 1;

	zCopy(&vv,v);

	zFree(&vv);
	zFree(&a);
	zFree(&b);
	return;
}

void zHex2Dec(z *u, z *v)
{
	//convert u[] in hex to v[] in decimal by repeatedly dividing
	//u by 1e9 = 0x3b9aca00
	//the remainder of the ith division is the ith decimal digit.
	//when the quotient = 0, stop

	z a,b;
	fp_digit r = 0;
	int su = abs(u->size);
	int approx_words = (int)((double)su * 1.5);	
	//because decimal takes more room than hex to store

	zInit(&a);
	zInit(&b);

	if (v->alloc < approx_words)
		zGrow(v,approx_words);
	zClear(v);

	if (a.alloc < approx_words)
	{
		zGrow(&a,approx_words);
		zClear(&a);
	}

	if (b.alloc < approx_words)
	{
		zGrow(&b,approx_words);
		zClear(&b);
	}
	
	zCopy(u,&a);
	v->size = 1;
	do
	{
		r = zShortDiv(&a,MAX_DEC_WORD,&b);
		v->val[v->size - 1] = r;
		v->size++;
		zCopy(&b,&a);
	} while (zCompare(&a,&zZero) != 0);
	v->size--;

	if (u->size < 0)
		v->size *= -1;

	zFree(&a);
	zFree(&b);
	return;
}

void swap(z *a, z *b)
{
	//do I actually have to physically copy here, or can I just swap pointers?
	z tmp;
	zInit(&tmp);
	zCopy(a,&tmp);
	zCopy(b,a);
	zCopy(&tmp,b);
	zFree(&tmp);
	return;
}

double rint(double x)
{
	 double i, r = modf(x, &i);
	 if (r < 0.0) {
		 r += 1.0; 
		 i -= 1.0;
	 }
	 return (r > 0.5 || (r == 0.5 && ((int)i & 1)) ? i + 1.0 : i);
}
