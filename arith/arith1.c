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
implements basic arithmetic operations
with a fundamental base of U32.  U64 types are used
for the operations that require a carry.
*/

#include "yafu.h"
#include "arith.h"
#include "common.h"

int zBits(z *n)
{
	if (n->size == 0)
		return 0;

	return BITS_PER_DIGIT*(abs(n->size)-1) + spBits(n->val[abs(n->size)-1]);
}

int gmp_base10(mpz_t x)
{
	mpz_t t;	//temp
	int g;		//guess: either correct or +1

	mpz_init(t);
	g = mpz_sizeinbase(x,10);
	mpz_set_ui(t,10);
	mpz_pow_ui(t, t, g - 1);
	g = g - (mpz_cmp(t, x) > 0 ? 1 : 0);
	mpz_clear(t);
	return g;
}
	
// borrowed from jasonp... 
double zlog(mpz_t x) {

#if GMP_LIMB_BITS == 32
	uint32 i = mpz_size(x);

	switch(i) {
	case 0:
		return 0;
	case 1:
		return log((double)((uint32)mpz_get_ui(x)));
	case 2:

		return log((double)(mpz_getlimbn(x,0)) + 
			MP_RADIX * mpz_getlimbn(x,1));

	default:

		return 32 * (i-3) * LN2 + 
			log((double)(mpz_getlimbn(x,i-3)) + MP_RADIX * (
		     	((double)mpz_getlimbn(x,i-2) + MP_RADIX * 
				mpz_getlimbn(x,i-1))));

	}

#else

	uint32 i = mpz_size(x);

	switch(i) {

	case 0:
		return 0;
	case 1:
		return log((double)((uint32)mpz_get_ui(x)));
	case 2:

		return log((double)(mpz_getlimbn(x,0)) + 
			18446744073709551616.0 * mpz_getlimbn(x,1));

	default:

		return 64 * (i-3) * LN2 + 
			log((double)(mpz_getlimbn(x,i-3)) + 18446744073709551616.0 * (
		     	((double)mpz_getlimbn(x,i-2) + 18446744073709551616.0 * 
				mpz_getlimbn(x,i-1))));

	}


#endif

}

fp_digit spBits(fp_digit n)
{	
	int i = 0;
	while (n != 0)
	{
		n >>= 1;
		i++;
	}
	
	return i;
}

int bits64(uint64 n)
{
	int i = 0;
	while (n != 0)
	{
		n >>= 1;
		i++;
	}
	return i;
}

void zAdd(z *u, z *v, z *w)
{
	//w = u + v
	//u and/or v can be negative
	//schoolbook addition.  see knuth TAOCP, vol. 2

	int i,su,sv,sw;
	fp_digit k;
	z *bigger, *smaller;

	if (u->size < 0)
	{
		if (v->size > 0)
		{
			//u is negative, v is not
			u->size *= -1;
			zSub(v,u,w);
			if (u != w)
				u->size *= -1;
			return;
		}
	}
	else if (v->size < 0)
	{
		//v is negative, u is not
		v->size *= -1;
		zSub(u,v,w);
		if (v != w)
			v->size *= -1;
		return;
	}

	//at this point, they are either both negative or positive
	//w is is same size as the greater of u or v, or at most
	//one limb bigger.
	//also, if u and v are different sizes, need to zero pad one or the other

	su = abs(u->size);
	sv = abs(v->size);

	if (su > sv)
	{
		sw = sv;
		bigger = u;
		smaller = v;
		//realloc w if necessary.  make it a fixed block of limbs bigger, to
		//prevent repeated reallocation if a bunch of adds are done in a row.
		if (w->alloc < (su + 1))
			zGrow(w,su + LIMB_BLKSZ);
	}
	else
	{
		sw = su;
		bigger = v;
		smaller = u;
		//realloc w if necessary.  make it a fixed block of limbs bigger, to
		//prevent repeated reallocation if a bunch of adds are done in a row.
		if (w->alloc < (sv + 1))
			zGrow(w,sv + LIMB_BLKSZ);
	}

#if defined(MSC_ASM32A) && !defined(ASM_ARITH_DEBUG)

	//add up to length of shorter input
	k=0;
	for (i=0;i<sw;++i)
		spAdd3(bigger->val[i],smaller->val[i],k,w->val+i,&k);

	/*
	k=0;
	for (i=0;i<sw;++i)
	{
		fp_digit s = smaller->val[i];
		fp_digit c = 0;
		fp_digit big = bigger->val[i];

		__asm
		{
			mov eax, big
			add eax, k
			adc c, 0
			add eax, s
			adc c, 0
		}

		w->val[i] = s;
		k=c;
	}
	*/

#elif defined(GCC_ASM32A) && !defined(ASM_ARITH_DEBUG)

	//add up to length of shorter input
	k=0;
	for (i=0;i<sw;++i)
	{
		fp_digit s = smaller->val[i];
		fp_digit c = 0;

		ASM_G ("movl %2, %%eax		\n\t"
			"addl %3, %%eax		\n\t"
			"adcl $0, %5		\n\t"
			"addl %%eax, %4		\n\t"
			"adcl $0, %5		\n\t"
			: "=r"(s), "=r"(c)
			: "r"(bigger->val[i]), "r"(k), "0"(s), "1"(c)
			: "eax", "memory", "cc");

		w->val[i] = s;
		k=c;
	}

#elif defined(GCC_ASM64X) && !defined(ASM_ARITH_DEBUG)

	//add up to length of shorter input
	k=0;
	for (i=0;i<sw;++i)
	{
		fp_digit s = smaller->val[i];
		fp_digit c = 0;

		ASM_G ("movq %2, %%rax		\n\t"
			"addq %3, %%rax		\n\t"
			"adcq $0, %5		\n\t"
			"addq %%rax, %4		\n\t"
			"adcq $0, %5		\n\t"
			: "=r"(s), "=r"(c)
			: "r"(bigger->val[i]), "r"(k), "0"(s), "1"(c)
			: "rax", "memory", "cc");

		w->val[i] = s;
		k=c;
	}
	
#else
	
	k=0;
	for (i=0;i<sw;++i)
		spAdd3(bigger->val[i],smaller->val[i],k,w->val+i,&k);

	
#endif

	if (abs(bigger->size) > sw)
	{
		//there are more words to process
		if (bigger->val[sw] != MAX_DIGIT)
		{
			//then just add any carry and copy the rest of the words
			i=sw;
			sw = abs(bigger->size);
			w->val[i] = bigger->val[i] + k;
			k=0;
			i++;
			for (; i<sw;i++)
				w->val[i] = bigger->val[i];
		}
		else
		{
			//there will be a carry we need to propagate, this will happen very rarely.
			i=sw;
			sw = abs(bigger->size);

#if defined(MSC_ASM32A)  && !defined(ASM_ARITH_DEBUG)

			for (;i<sw;++i)
				spAdd3(bigger->val[i],0,k,w->val+i,&k);


#elif defined(GCC_ASM32A) && !defined(ASM_ARITH_DEBUG)

			for (;i<sw;++i)
			{
				fp_digit c = 0;
				fp_digit s = 0;

				ASM_G ("movl %2, %%eax		\n\t"
					"addl %3, %%eax		\n\t"
					"movl %%eax, %1		\n\t"
					"adcl $0, %4		\n\t"
					: "=r"(c), "=r"(s)
					: "r"(bigger->val[i]), "r"(k), "0"(c)
					: "eax", "memory", "cc");

				w->val[i] = s;
				k=c;
			}

#elif defined(GCC_ASM64X) && !defined(ASM_ARITH_DEBUG)

			for (;i<sw;++i)
			{
				fp_digit c = 0;
				fp_digit s = 0;

				ASM_G ("movq %2, %%rax		\n\t"
					"addq %3, %%rax		\n\t"
					"movq %%rax, %1		\n\t"
					"adcq $0, %4		\n\t"
					: "=r"(c), "=r"(s)
					: "r"(bigger->val[i]), "r"(k), "0"(c)
					: "rax", "memory", "cc");

				w->val[i] = s;
				k=c;
			}
			
#else

			for (;i<sw;++i)
				spAdd3(bigger->val[i],0,k,w->val+i,&k);

#endif
		}
	}

	w->size = abs(bigger->size);

	//if there is still a carry, add a digit to w.
	if (k)
	{
		w->val[sw]=1;
		w->size++;
	}

	//if the inputs are negative, so is the output.
	if (u->size < 0)
		w->size *= -1;

	w->type = UNKNOWN;

	return;
}

void zShortAdd(z *u, fp_digit v, z *w)
{
	//w = u + v
	//schoolbook addition.  see knuth TAOCP, vol. 2

	int i,su;
	fp_digit k;

	if (u->size < 0)
	{
		//u is negative
		u->size *= -1;
		zShortSub(u,v,w);
		w->size *= -1;
		u->size *= -1;
		return;
	}

	su = abs(u->size);

	//realloc w if necessary.  make it a fixed block of limbs bigger, to
	//prevent repeated reallocation if a bunch of adds are done in a row.
	if (w->alloc < (su + 1))
		zGrow(w,su + LIMB_BLKSZ);

	zCopy(u,w);

	//add
	spAdd(u->val[0],v,w->val,&k);

	//add the carry
	spAdd(u->val[1],k,w->val+1,&k);

	if (k)
	{
		//only rarely will the carry propagate more than one place
		//special case this.
		for (i=2;i<su;++i)
			spAdd(u->val[i],k,w->val+i,&k);
	}

	//if there is still a carry, add a digit to w.
	if (k) 
	{
		w->val[su]=k;
		w->size++;
	}

	w->type = UNKNOWN;

	return;
}

int zSub(z *u, z *v, z *w)
{
	//w = u-v
	//u and/or v can be negative
	//schoolbook subtraction, see knuth TAOCP, vol. 2

	fp_digit k=0;
	int i,j,m,su,sv,sign=0;
	z *bigger, *smaller;

	if (u->size < 0)
	{
		if (v->size > 0)
		{
			//u is negative, v is not, so really an addition
			u->size *= -1;
			zAdd(u,v,w);
			if (u != w)
				u->size *= -1;
			w->size *= -1;
			return 0;
		}
		else
		{
			//both are negative, so we really have -u + v or v - u
			v->size *= -1;
			u->size *= -1;
			zSub(v,u,w);
			if (v != w)
				v->size *= -1;
			if (u != w)
				u->size *= -1;
			return 0;
		}
	}
	else if (v->size < 0)
	{
		if (u->size > 0)
		{
			//v is negative, u is not, so really an addition
			v->size *= -1;
			zAdd(u,v,w);
			if (v != w)
				v->size *= -1;
			return 0;
		}
	}

	su = u->size;
	sv = v->size;
	
	if (su > sv) 
	{	
		bigger = u;
		smaller = v;
		//realloc w if necessary.  make it a fixed block of limbs bigger, to
		//prevent repeated reallocation if a bunch of adds are done in a row.
		if (w->alloc < (su + 1))
			zGrow(w,su + LIMB_BLKSZ);
		goto beginsub;
	}
	if (su < sv)
	{
		bigger = v;
		smaller = u;
		//realloc w if necessary.  make it a fixed block of limbs bigger, to
		//prevent repeated reallocation if a bunch of adds are done in a row.
		if (w->alloc < (sv + 1))
			zGrow(w,sv + LIMB_BLKSZ);
		sign=1;
		goto beginsub;
	}

	//equal number of words to this point
	if (w->alloc < (su + 1))
		zGrow(w,su + LIMB_BLKSZ);

	for (i = su - 1; i>=0; --i)
	{
		if (u->val[i] > v->val[i]) 
		{	
			bigger = u;
			smaller = v;
			goto beginsub;
		}
		if (u->val[i] < v->val[i])
		{	
			bigger = v;
			smaller = u;
			sign=1;
			goto beginsub;
		}
	}

	//equal if got to here
	w->size = 1;
	w->val[0] = 0;
	return 1;

beginsub:
	//subtract up to the length of the smaller
	m = smaller->size;

#if defined(MSC_ASM32A) && !defined(ASM_ARITH_DEBUG)

	for (j=0;j<m;++j)
		spSub3(bigger->val[j],smaller->val[j],k,w->val+j,&k);

#elif defined(GCC_ASM32A) && !defined(ASM_ARITH_DEBUG)

	for (j=0; j<m; j++)
	{
		fp_digit s = smaller->val[j];
		fp_digit b = 0;

		ASM_G ("movl %2, %%eax		\n\t"
			"subl %4, %%eax		\n\t"
			"adcl $0, %5		\n\t"
			"subl %3, %%eax		\n\t"
			"adcl $0, %5		\n\t"
			"movl %%eax, %4		\n\t"
			: "=r"(s), "=r"(b)
			: "r"(bigger->val[j]), "r"(k), "0"(s), "1"(b)
			: "eax", "memory", "cc");

		w->val[j] = s;
		k = b;
	}

#elif defined(GCC_ASM64X) && !defined(ASM_ARITH_DEBUG)

	for (j=0; j<m; j++)
	{
		fp_digit s = smaller->val[j];
		fp_digit b = 0;

		ASM_G ("movq %2, %%rax		\n\t"
			"subq %4, %%rax		\n\t"
			"adcq $0, %5		\n\t"
			"subq %3, %%rax		\n\t"
			"adcq $0, %5		\n\t"
			"movq %%rax, %4		\n\t"
			: "=r"(s), "=r"(b)
			: "r"(bigger->val[j]), "r"(k), "0"(s), "1"(b)
			: "rax", "memory", "cc");

		w->val[j] = s;
		k = b;
	}

#else

	for (j=0;j<m;++j)
		spSub3(bigger->val[j],smaller->val[j],k,w->val+j,&k);

	
#endif
	
	//if there is a leftover word that is != 0, then subtract any
	//carry and simply copy any other leftover words

	//if there is a leftover word that is == 0, then subtract with
	//borrow for the rest of the leftover words.  this will happen rarely

	//leftover word?
	if (bigger->size > m)
	{
		//not equal to zero?
		if (bigger->val[m] != 0)
		{
			//subtract any carry and copy the rest
			w->val[m] = bigger->val[m] - k;
			j=m+1;
			m=bigger->size;
			for(; j<m; j++)
				w->val[j] = bigger->val[j];
		}
		else
		{
			//equal to zero, need to subtract with borrow for the rest
			//of the leftover words.  
			j=m;
			m=bigger->size;

#if defined(MSC_ASM32A)  && !defined(ASM_ARITH_DEBUG)

			for (;j<m;++j)
				spSub3(bigger->val[j],0,k,w->val+j,&k);

#elif defined(GCC_ASM32A) && !defined(ASM_ARITH_DEBUG)

			for (;j<m;++j)
			{
				fp_digit s = 0;
				fp_digit b = 0;

				ASM_G ("movl %2, %%eax		\n\t"
					"subl %3, %%eax		\n\t"
					"adcl $0, %5		\n\t"
					"movl %%eax, %0		\n\t"
					: "=r"(s), "=r"(b)
					: "r"(bigger->val[j]), "r"(k), "0"(s), "1"(b)
					: "eax", "memory", "cc");

				w->val[j] = s;
				k = b;
			}

#elif defined(GCC_ASM64X) && !defined(ASM_ARITH_DEBUG)

			for (;j<m;++j)
			{
				fp_digit s = 0;
				fp_digit b = 0;

				ASM_G ("movq %2, %%rax		\n\t"
					"subq %3, %%rax		\n\t"
					"adcq $0, %5		\n\t"
					"movq %%rax, %0		\n\t"
					: "=r"(s), "=r"(b)
					: "r"(bigger->val[j]), "r"(k), "0"(s), "1"(b)
					: "rax", "memory", "cc");

				w->val[j] = s;
				k = b;
			}

#else
			for (;j<m;++j)
				spSub3(bigger->val[j],0,k,w->val+j,&k);

			
#endif
		}
	}
	
	w->size = bigger->size;
	fp_clamp(w);
	if (sign)
		w->size *= -1;

	if (w->size == 0)
		w->size = 1;

	w->type = UNKNOWN;
	return 0;	
}

void zShortSub(z *u, fp_digit v, z *w)
{
	//w = u - v
	//schoolbook subtraction.  see knuth TAOCP, vol. 2

	int i,su;
	fp_digit k=0;

	su = abs(u->size);
	w->size = su;

	if (u->size < 0)
	{
		//u is negative, really an addition
		u->size *= -1;
		zShortAdd(u,v,w);
		u->size *= -1;
		w->size *= -1;
		return;
	}

	zCopy(u,w);

	//subtract
	spSub3(u->val[0],v,0,w->val,&k);

	//subtract the borrow
	spSub3(u->val[1],k,0,w->val+1,&k);
	
	if (k)
	{
		//propagate the borrow
		for (i=2;i<su;++i)
			spSub3(u->val[i],0,k,w->val+i,&k);
	}

	//check if we lost the high digit
	if ((w->val[su - 1] == 0) && (su != 1))
		su--;
	w->size = su;

	//check for u < v
	if (k)
	{
		//then u < v, and result is negative
		w->val[0] = ~w->val[0];
		w->val[0]++;
		w->size *= -1;
	}

	w->type = UNKNOWN;

	return;
}


fp_digit zShortDiv(z *u, fp_digit v, z *q)
{
	//q = u/v
	//return the remainder
	//schoolbook long division.  see knuth TAOCP, vol. 2

	int su = abs(u->size);
	int i;
	fp_digit rem = 0;

	q->size = su;

	//realloc q if necessary
	if (q->alloc < (su + 1))
		zGrow(q,su + LIMB_BLKSZ);

	i = su - 1;
	if (u->val[i] < v) 
	{
		rem = u->val[i];
		q->val[i--] = 0;
	}

#if defined(MSC_ASM32A) 

	while (i >= 0)
	{
		fp_digit quot1 = u->val[i];
		ASM_M  {
			mov eax, quot1
			mov edx, rem
			div v
			mov rem, edx
			mov quot1, eax
		}
		q->val[i] = quot1;
		i--;
	}


#elif defined(GCC_ASM32A)

	while (i >= 0)
	{
		fp_digit quot1;
		ASM_G ("divl %4"
			: "=a"(quot1),"=d"(rem)
			: "1"(rem), "0"(u->val[i]), "r"(v) );

		q->val[i] = quot1;
		i--;
	}


#elif defined(GCC_ASM64X)

	while (i >= 0)
	{
		fp_digit quot1;
		ASM_G ("divq %4"
			: "=a"(quot1),"=d"(rem)
			: "1"(rem), "0"(u->val[i]), "r"(v) );

		q->val[i] = quot1;
		i--;
	}


	
#else
	
	
	while (i >= 0)
	{
		fp_digit uu[2];
		uu[1] = (fp_digit)rem;
		uu[0] = u->val[i];

		spDivide(q->val + i, &rem, uu, v);
		i--;
	}

#endif

	//the quotient could be one limb smaller than the input
	if ((q->val[q->size - 1] == 0) && (q->size != 1))
		q->size--;

	if (u->size < 0)
		q->size *= -1;

	q->type = UNKNOWN;
	return rem;
}

uint32 zShortDiv32(z32 *u, uint32 v, z32 *q)
{
	//q = u/v
	//return the remainder
	//schoolbook long division.  see knuth TAOCP, vol. 2
	//only used in the trial division stage of QS

	int su = u->size;
	int i;
	uint32 rem = 0;

	q->size = su;

	i = su - 1;
	if (u->val[i] < v) 
	{
		rem = u->val[i];
		q->val[i--] = 0;
	}

#if defined(MSC_ASM32A) 

	while (i >= 0)
	{
		fp_digit quot1 = u->val[i];
		ASM_M  {
			mov eax, quot1
			mov edx, rem
			div v
			mov rem, edx
			mov quot1, eax
		}
		q->val[i] = quot1;
		i--;
	}


#elif defined(GCC_ASM32A)

	while (i >= 0)
	{
		uint32 quot1;
		ASM_G ("divl %4"
			: "=a"(quot1),"=d"(rem)
			: "1"(rem), "0"(u->val[i]), "r"(v) );

		q->val[i] = quot1;
		i--;
	}


#elif defined(GCC_ASM64X)

	while (i >= 0)
	{
		uint32 quot1;
		ASM_G ("divl %4"
			: "=a"(quot1),"=d"(rem)
			: "1"(rem), "0"(u->val[i]), "r"(v) );

		q->val[i] = quot1;
		i--;
	}


#else

	while (i >= 0)
	{
		fp_word acc = (fp_word)rem << 32 | (fp_word)u->val[i];
		q->val[i] = (uint32)(acc / v);
		rem = (uint32)(acc % v);
		i--;
	}
	
#endif

	//the quotient could be one limb smaller than the input
	if ((q->val[q->size - 1] == 0) && (q->size != 1))
		q->size--;

	q->type = UNKNOWN;

	return rem;
}

fp_digit zShortMod(z *u, fp_digit v)
{
	//q = u/v
	//return the remainder
	int i;

#if defined(MSC_ASM32A) 

	fp_digit rem = 0;

	i = abs(u->size) - 1;
	if (u->val[i] < v) 
		rem = u->val[i--];

	while (i >= 0)
	{
		fp_digit quot1 = u->val[i];
		ASM_M  {
			mov eax, quot1
			mov edx, rem
			div v
			mov rem, edx
		}
		i--;
	}

	return rem;

#elif defined(GCC_ASM32A)

	fp_digit rem = 0;

	i = abs(u->size) - 1;
	if (u->val[i] < v) 
		rem = u->val[i--];

	while (i >= 0)
	{
		fp_digit quot1;
		ASM_G ("divl %4"
			: "=a"(quot1),"=d"(rem)
			: "1"(rem), "0"(u->val[i]), "r"(v) );
		i--;
	}

	return rem;

#elif defined(GCC_ASM64X)

	fp_digit rem = 0;

	i = abs(u->size) - 1;
	if (u->val[i] < v) 
		rem = u->val[i--];

	while (i >= 0)
	{
		fp_digit quot1;
		ASM_G ("divq %4"
			: "=a"(quot1),"=d"(rem)
			: "1"(rem), "0"(u->val[i]), "r"(v) );
		i--;
	}

	return rem;

#elif defined(_WIN64)
	
	uint64 rem = 0;

	i = abs(u->size) - 1;
	if (u->val[i] < v) 
		rem = (uint64)u->val[i--];

	while (i >= 0)
	{
		uint64 d = (rem << 32) | (u->val[i]);
		rem = mod_64(d, v);
		i--;
	}

	return (uint32)rem;


#else
	
	fp_digit qq[2], q64, rem;
	
	if (u->size == 1)
		return u->val[0] % v;

	qq[1] = (fp_digit)u->val[u->size - 1];
	for (i=u->size - 2;i>=0;i--)
	{
		qq[0] = u->val[i];
		spDivide(&q64, &rem, qq, v);

		qq[1] = rem;
	}

	return (fp_digit)rem;

#endif

}

uint32 zShortMod32(z32 *u, uint32 v)
{
	//q = u/v
	//return the remainder
	//only used during trial division stage of QS
	int i;

#if defined(MSC_ASM32A) 

	uint32 rem = 0;

	//u will always be positive
	i = u->size - 1;
	if (u->val[i] < v) 
		rem = u->val[i--];

	while (i >= 0)
	{
		uint32 quot1 = u->val[i];
		ASM_M  {
			mov eax, quot1
			mov edx, rem
			div v
			mov rem, edx
		}
		i--;
	}

	return rem;

#elif defined(GCC_ASM32A)

	uint32 rem = 0;

	i = abs(u->size) - 1;
	if (u->val[i] < v) 
		rem = u->val[i--];

	while (i >= 0)
	{
		uint32 quot1;
		ASM_G ("divl %4"
			: "=a"(quot1),"=d"(rem)
			: "1"(rem), "0"(u->val[i]), "r"(v) );
		i--;
	}

	return rem;

#elif defined(GCC_ASM64X)

	uint32 rem = 0;

	i = abs(u->size) - 1;
	if (u->val[i] < v) 
		rem = u->val[i--];

	while (i >= 0)
	{
		uint32 quot1;
		ASM_G ("divl %4"
			: "=a"(quot1),"=d"(rem)
			: "1"(rem), "0"(u->val[i]), "r"(v) );
		i--;
	}

	return rem;


#else
	
	fp_word q64;

	if (u->size == 1)
		return u->val[0] % v;

	q64 = (fp_word)u->val[u->size - 1] << 32;
	for (i=u->size - 2;i>=0;i--)
	{
		q64 |= u->val[i];
		q64 = (q64 % v) << 32;
	}

	return (fp_digit)(q64 >> 32);
	
#endif

}

uint32 zShortEDiv32(z32 *u, uint32 v, uint32 inv)
{
	//u is overwritten with the quotient (exact) of u / v;
	uint64 q64;
	uint32 qhat;
	int i,su = u->size;

	//this may not be the fastest implementation of the exact division algorithm
	//but it works, and is over twice as fast as mpDiv_1
	for (i=0;i<su;i++)
	{
		qhat = inv*u->val[i];
		q64 = ((uint64)u->val[i+2]) << 32;
		q64 |= (uint64)u->val[i+1];
		q64 = q64 - (((uint64)qhat * (uint64)v) >> 32);
		u->val[i+2] = (uint32)(q64 >> 32);
		u->val[i+1] = (uint32)(q64 & 0xFFFFFFFF);
		u->val[i] = qhat;
	}
		
	if (u->val[su-1] == 0 && u->size > 1)
			u->size--;

	return 1;
}

void zDiv(z *u, z *v, z *q, z *r)
{
	/*
	q = u \ v
	r = u mod v
	u is overwritten

	schoolbook long division.  see knuth TAOCP, vol. 2
	*/
	fp_digit v1=0,v2=0,k,qhat,rhat,uj2,tt[2],pp[2];
	int i,j,m,su,sv;
	int s =0,cmp,sdd,sd;
	unsigned int shift;
	fp_digit bitmask;
	
	su = abs(u->size);
	sv = abs(v->size);
	m = su-sv;

	//Catch special cases 	
	//divide by zero
	if ((sv == 0) || ((sv == 1) && (v->val[0] == 0)))
		return;	

	if (q->alloc < (m + 1))
		zGrow(q,m + 1);

	if (r->alloc < su)
		zGrow(r,su);

	//Use short division instead
	if (sv == 1)
	{	
		r->val[0] = zShortDiv(u,v->val[0],q);
		r->size = 1;
		s = (v->size < 0);
		if (s)
		{
			q->size *= -1;
			r->size *= -1;
		}
		return;
	}

	//v > u, so just set q = 0 and r = u
	if (su < sv)
	{	
		zClear(q);
		q->size = 1;
		zCopy(u,r);
		s = (u->size < 0) ^ (v->size < 0);
		if (s)
		{
			r->size *= -1;
		}
		q->type = UNKNOWN;
		r->type = UNKNOWN;
		return;
	}

	//u and v are the same length
	if (su == sv)
	{	
		cmp = zCompare(u,v);
		//v > u, as above
		if (cmp < 0)
		{	
			zClear(q);
			q->size = 1;
			zCopy(u,r);
			s = (u->size < 0) ^ (v->size < 0);
			if (s)
			{
				q->size *= -1;
				r->size *= -1;
			}
			q->type = UNKNOWN;
			r->type = UNKNOWN;
			return;
		} 
		else if (cmp == 0)	//v == u, so set q = 1 and r = 0
		{	
			q->size = 1;
			q->val[0]=1;
			r->size = 1;
			r->val[0]=0;
			s = (u->size < 0) ^ (v->size < 0);
			if (s)
			{
				q->size *= -1;
			}
			q->type = UNKNOWN;
			r->type = UNKNOWN;
			return;
		}
	}

	//normalize v by left shifting until the high bit of v is set (v1 >= floor(2^31))
	bitmask = HIBITMASK;
	for (shift = 0; shift < BITS_PER_DIGIT; ++shift)
	{
		if (v->val[sv-1] & bitmask)
			break;
		bitmask >>= 1;
	}

	//normalize v by shifting left (x2) shift number of times
	//overflow should never occur to v during normalization
	zShiftLeft_x(v,v,shift);

	//left shift u the same amount - may get an overflow here
	zShiftLeft_x(u,u,shift);

	if (abs(u->size) == su)
	{	//no overflow - force extra digit
		if (u->size < 0)
			u->size--;
		else
			u->size++;
		u->val[su] = 0;
		su++;
	}
	else
		su++;

	//copy first two digits of v to local variables for quick access
	v1=v->val[sv-1];
	v2=v->val[sv-2];
	
	sdd=0;
	sd=0;
	//main loop
	for (j=0;j<=m;++j)
	{
		//calculate qhat
		tt[1] = u->val[su-j-1];		//first digit of normalized u
		tt[0] = u->val[su-j-2];		//second digit of normalized u
		if (tt[1] == v1)
			qhat = MAX_DIGIT;
		else
			spDivide(&qhat, &rhat, tt, v1);

		//quick check if qhat is too big based on our initial guess involving 
		//the first two digits of u and v.
		uj2 = u->val[su-j-3];

		while (1)
		{
			spMultiply(qhat,v1,&pp[0],&pp[1]);
			shortSubtract(tt,pp,tt);
			if (tt[1]) break;
			tt[1] = tt[0]; tt[0] = uj2; 
			spMultiply(qhat,v2,&pp[0],&pp[1]);
			i = shortCompare(pp,tt);  //p = v2*qhat, t = (uj*b+uj1-qhat*v1)*b + uj2

			if (i == 1)
				qhat--;
			else
				break;
		}
	
		//keep track of the significant digits
		if (qhat > 0)
		{
			sdd = sdd + 1 + sd;
			sd = 0;
		}
		else if (sdd != 0)
			sd++;

		//multiply and subtract, in situ
		k=0;
		for (i=0;i<sv;++i)
		{
			s=(int)(su-j-sv+i-1);
			spMultiply(v->val[i],qhat,&pp[0],&pp[1]);
			spAdd(pp[0],k,&tt[0],&tt[1]);
			u->val[s] = u->val[s] - tt[0]; 
			//check if this result is negative, remember the borrow for the next digit
			if (u->val[s] > (u->val[s] + tt[0]))
				k = pp[1] + tt[1] + 1;
			else
				k = pp[1] + tt[1];
		}
		
		//if the final carry is bigger than the most significant digit of u, then qhat
		//was too big, i.e. qhat[v1v2...vn] > [u0u1u2...un]
		if (k > u->val[su-j-1])
		{
			//correct by decrementing qhat and adding back [v1v2...vn] to [u0u1...un]
			qhat--;
			//first subtract the final carry, yielding a negative number for [u0u1...un]
			u->val[su-j-1] -= k;
			//then add back v
			k=0;
			for (i=0;i<sv;i++)
				spAdd3(u->val[su-j-sv+i-1],v->val[i],k,&u->val[su-j-sv+i-1],&k);
			u->val[su-j-1] += k;
		}
		else //else qhat was ok, subtract the final carry
			u->val[su-j-1] -= k;

		//set digit of q
		q->val[m-j] = qhat;
	}
	q->size = sdd+sd;
	zCopy(u,r);

	for (s=r->size - 1; s>=0; --s)
	{
		if ((r->val[s] == 0) && (r->size > 0))
			r->size--;
		else
			break;
	}

	//unnormalize.
	zShiftRight_x(v,v,shift);
	zShiftRight_x(r,r,shift);

	if (r->size == 0)
		r->size = 1;

	s = (u->size < 0) ^ (v->size < 0);
	if (s)
	{
		q->size *= -1;
		r->size *= -1;
	}

	q->type = UNKNOWN;
	r->type = UNKNOWN;
	return;
}

int shortCompare(fp_digit p[2], fp_digit t[2])
{
	//utility function used in zDiv
	int i;

	for (i=1;i>=0;--i)
	{
		if (p[i] > t[i]) return 1;
		if (p[i] < t[i]) return -1;
	}
	return 0;
}

int shortSubtract(fp_digit u[2], fp_digit v[2], fp_digit w[2])
{
	//utility function used in zDiv
	fp_digit j=0;

	w[0] = u[0] - v[0];
	if (w[0] > (MAX_DIGIT - v[0]))
	{
		j=1;
		w[0] = w[0] + MAX_DIGIT + 1;
	}
	w[1] = u[1] - v[1] - j;
	
	return 1;
}

void zMul(z *u, z *v, z *w)
{
	//w = u*v
	//schoolbook multiplication, see knuth TAOCP, vol. 2
	int su, sv;
	int signu, signv;
	int i,j;
	fp_digit k=0;
	fp_digit *wptr;
	z tmp;
	zInit(&tmp);

	su = abs(u->size);
	sv = abs(v->size);
	signu = u->size < 0;
	signv = v->size < 0;


	if (tmp.alloc < (su+sv))
		zGrow(&tmp,su+sv);

	zClear(&tmp);

	if (w->alloc < (su+sv))
		zGrow(w,su+sv);

	if (u->alloc < v->alloc)
		zGrow(u,v->alloc);

	if (v->alloc < u->alloc)
		zGrow(v,u->alloc);

	//use short multiplication instead
	if (su == 1) 
	{
		zShortMul(v,u->val[0],w);
		if (signu)
			w->size *= -1;

		zFree(&tmp);
		return;
	}

	//use short multiplication instead
	if (sv == 1)
	{
		zShortMul(u,v->val[0],w);
		if (signv)
			w->size *= -1;

		zFree(&tmp);
		return;
	}
	
	//for each digit of u
	for (i=0;i<su;++i)
	{
		//take an inner product and add in-situ with the previous inner products
		k=0;
		wptr = &tmp.val[i];
		for (j=0;j<sv;++j)
			spMulAdd(u->val[i],v->val[j],wptr[j],k,&wptr[j],&k);
		wptr[j] += k;
	}
	tmp.size = su+sv;
	fp_clamp(&tmp);
	zCopy(&tmp,w);
	zFree(&tmp);	
	
	//set sign
	if (signu ^ signv)
		w->size *= -1;

	w->type = UNKNOWN;

	return;
}

void zModMul(z *u, z *v, z *n, z *w)
{
	z t1,t2;
	
	zInit(&t1);
	zInit(&t2);

	zMul(u,v,&t1);
	zDiv(&t1,n,&t2,w);

	zFree(&t1);
	zFree(&t2);

	return;
}

void zShortMul(z *u, fp_digit v, z *w)
{
	//w = u * v
	//schoolbook multiplication, see knuth TAOCP, vol. 2
	fp_digit k=0;
	long i;
	long su;

	su = abs(u->size);

	//initialize w
	if (w->alloc < (su+1))
		zGrow(w,su + LIMB_BLKSZ);

	//inner product
	for (i=0;i<su;++i)
		spMulAdd(u->val[i],v,0,k,&w->val[i],&k);

	//if still have a carry, add a digit to w
	if (k)
	{
		w->val[su]=k;
		su++;
	}
	
	//check for significant digits.  only necessary if v or u = 0?
	for (i = su - 1;i>=0;--i)
	{
		if (w->val[i] != 0)
			break;
	}
	w->size = i+1;

	//only necessary if v or u = 0?
	if (w->size == 0)
		w->size = 1;

	if (u->size < 0)
		w->size *= -1;

	w->type = UNKNOWN;

	return;
}

void zSqr(z *x, z *y)
{
	//w = x * x
	fp_digit v,u[2],c[2],highbitv,k;
	fp_digit *wptr;
	int i,j,t;
	z *w, tmp;

	zMul(x,x,y);
	return;

	//from the handbook of applied cryptography:
	//t = size of x
	//for i from 0 to 2t-1 do
	//	w[i]=0;
	//for i from 0 to t-1 do:
	//	(uv)b = w[2i] + x[i] * x[i]
	//	w[2i]=v
	//	c=u
	//	for j from (i+1) to (t-1) do:
	//		(uv)b = w[i+j] + 2x[j]*x[i] + c
	//		w[i+j]=v
	//		c=u
	//	end
	//	w[i+t]=u
	//end

	//need to accomodate the 2*xj*xi operation, as it can take 3 digits to
	//represent

	if (x == y)
	{ 
		zInit(&tmp);
		w = &tmp;
	}
	else
		w = y;

	c[0] = 0; c[1] = 0; u[0] = 0; u[1] = 0;
	t = abs(x->size);

	//output is always positive
	w->size = 2*t;
	
	//initialize w
	if (w->alloc < (2*t))
	{
		w->val = (fp_digit *)realloc(w->val,(2*t) * sizeof(fp_digit));
		w->alloc = 2*t;
	}
	memset(w->val,0,(2*t)*sizeof(fp_digit));

	for (i=0;i<t;i++)
	{
		wptr = w->val + (i<<1);

		spMultiply(x->val[i],x->val[i],&v,&u[0]);
		spAdd(*wptr,v,wptr,&c[0]);
		u[0]+=c[0];
		c[0]=u[0];
		c[1]=0;

		wptr = w->val + i;
		for (j=(i+1);j<t;j++)
		{
			spMultiply(x->val[j],x->val[i],&v,&u[0]);
			//now the left shift of both v and u, keeping track of carries
			//highbitv = ((v & 0x80000000) > 0);
			highbitv = (v & HIBITMASK) >> (BITS_PER_DIGIT - 1);
			v <<=1;
			//u[1] = ((u[0] & 0x80000000) > 0);
			u[1] = (u[0] & HIBITMASK) >> (BITS_PER_DIGIT - 1);
			u[0] <<= 1;
			//u is at most 0xFFFFFFFE, so we can add the carry from the v shift
			//without overflow
			u[0] += highbitv;
			//now v and u hold 2*xj*xi.  
			spAdd3(wptr[j],c[0],v,wptr+j,&k);
			spAdd3(u[0],c[1],k,&u[0],&k);
			u[1] += k;
			c[0] = u[0];
			c[1] = u[1];
		}
		wptr[t] += u[0];
		wptr[t+1] += u[1];
	}

	for (i = w->size - 1;i>=0;--i)
	{
		if (w->val[i] != 0)
			break;
	}
	w->size = i+1;

	zCopy(w,y);

	if (x == y)
		zFree(&tmp);

	return;
}

void zNeg(z *u)
{
	//two's complement negation
	//set u = -u
	//-u = ~u + 1

	int i,s=1;
	int su;

	//keep track of significant digits with s and su
	su = abs(u->size);
	for (i= su-1;i>=0;--i)
	{
		u->val[i] = ~u->val[i];
		if (!(u->val[i]) && s)
			su--;
		else
			s=0;
	}
	u->size = su;

	if (u->val[0] < MAX_DIGIT)
		u->val[0]++;
	else
		zShortAdd(u,1,u);

	if (u->size == 0)
		u->size = 1;

	u->size *= -1;
	return;
}

void zShiftLeft_1(z *a, z *b)
{	
	/* Computes a = b << 1 */
	int sign;
	int y;
	int sb,j;
	fp_digit carry, nextcarry;

	if (b->size < 0)
	{
		sb = -1 * b->size;
		sign = 1;
	}
	else
	{
		sb = b->size;
		sign = 0;
	}
	a->size = sb;

	//for each digit, remember the highest x bits using the mask, then shift.
	//the highest x bits becomes the lowest x bits for the next digit
	y = BITS_PER_DIGIT - 1;
	carry = 0;
	for (j = 0; j < sb; ++j)
	{
		nextcarry = (b->val[j] & HIBITMASK) >> y;
		a->val[j] = b->val[j] << 1 | carry;
		carry = nextcarry;
	}
	
	if (carry)
	{
		a->val[sb] = carry;
		a->size++;
	}

	if (sign)
		a->size *= -1;

	a->type = COMPOSITE;
	return;
}

static uint64 LEFTSHIFTMASKS64[64] = {
	0x0000000000000000ULL, 0x8000000000000000ULL, 0xC000000000000000ULL, 0xE000000000000000ULL, 
	0xF000000000000000ULL, 0xF800000000000000ULL, 0xFC00000000000000ULL, 0xFE00000000000000ULL, 
	0xFF00000000000000ULL, 0xFF80000000000000ULL, 0xFFC0000000000000ULL, 0xFFE0000000000000ULL, 
	0xFFF0000000000000ULL, 0xFFF8000000000000ULL, 0xFFFC000000000000ULL, 0xFFFE000000000000ULL,
	0xFFFF000000000000ULL, 0xFFFF800000000000ULL, 0xFFFFC00000000000ULL, 0xFFFFE00000000000ULL, 
	0xFFFFF00000000000ULL, 0xFFFFF80000000000ULL, 0xFFFFFC0000000000ULL, 0xFFFFFE0000000000ULL, 
	0xFFFFFF0000000000ULL, 0xFFFFFF8000000000ULL, 0xFFFFFFC000000000ULL, 0xFFFFFFE000000000ULL, 
	0xFFFFFFF000000000ULL, 0xFFFFFFF800000000ULL, 0xFFFFFFFC00000000ULL, 0xFFFFFFFE00000000ULL,
	0xFFFFFFFF00000000ULL, 0xFFFFFFFF80000000ULL, 0xFFFFFFFFC0000000ULL, 0xFFFFFFFFE0000000ULL, 
	0xFFFFFFFFF0000000ULL, 0xFFFFFFFFF8000000ULL, 0xFFFFFFFFFC000000ULL, 0xFFFFFFFFFE000000ULL, 
	0xFFFFFFFFFF000000ULL, 0xFFFFFFFFFF800000ULL, 0xFFFFFFFFFFC00000ULL, 0xFFFFFFFFFFE00000ULL, 
	0xFFFFFFFFFFF00000ULL, 0xFFFFFFFFFFF80000ULL, 0xFFFFFFFFFFFC0000ULL, 0xFFFFFFFFFFFE0000ULL,
	0xFFFFFFFFFFFF0000ULL, 0xFFFFFFFFFFFF8000ULL, 0xFFFFFFFFFFFFC000ULL, 0xFFFFFFFFFFFFE000ULL, 
	0xFFFFFFFFFFFFF000ULL, 0xFFFFFFFFFFFFF800ULL, 0xFFFFFFFFFFFFFC00ULL, 0xFFFFFFFFFFFFFE00ULL, 
	0xFFFFFFFFFFFFFF00ULL, 0xFFFFFFFFFFFFFF80ULL, 0xFFFFFFFFFFFFFFC0ULL, 0xFFFFFFFFFFFFFFE0ULL, 
	0xFFFFFFFFFFFFFFF0ULL, 0xFFFFFFFFFFFFFFF8ULL, 0xFFFFFFFFFFFFFFFCULL, 0xFFFFFFFFFFFFFFFEULL};

static uint32 LEFTSHIFTMASKS32[32] = {
	0x00000000, 0x80000000, 0xC0000000, 0xE0000000,
	0xF0000000, 0xF8000000, 0xFC000000, 0xFE000000,
	0xFF000000, 0xFF800000, 0xFFC00000, 0xFFE00000,
	0xFFF00000, 0xFFF80000, 0xFFFC0000, 0xFFFE0000,
	0xFFFF0000, 0xFFFF8000, 0xFFFFC000, 0xFFFFE000,
	0xFFFFF000, 0xFFFFF800, 0xFFFFFC00, 0xFFFFFE00,
	0xFFFFFF00, 0xFFFFFF80, 0xFFFFFFC0, 0xFFFFFFE0,
	0xFFFFFFF0, 0xFFFFFFF8, 0xFFFFFFFC, 0xFFFFFFFE};

void zShiftLeft_x(z *a, z *b, int x)
{	
	/* Computes a = b << x, where x is less than a word */
	int sign;
	int y;
	int sb,j;
	fp_digit mask, carry, nextcarry;
	
#if BITS_PER_DIGIT == 64
	mask = LEFTSHIFTMASKS64[x];
#else
	mask = LEFTSHIFTMASKS32[x];
#endif

	if (b->size < 0)
	{
		sb = -1 * b->size;
		sign = 1;
	}
	else
	{
		sb = b->size;
		sign = 0;
	}
	a->size = sb;

	if (a->alloc < (sb + 1))
		zGrow(a, sb + LIMB_BLKSZ);

	//for each digit, remember the highest x bits using the mask, then shift.
	//the highest x bits becomes the lowest x bits for the next digit
	y = BITS_PER_DIGIT - x;
	carry = 0;
	for (j = 0; j < sb; ++j)
	{
		nextcarry = (b->val[j] & mask) >> y;
		a->val[j] = b->val[j] << x | carry;
		carry = nextcarry;
	}
	
	if (carry)
	{
		a->val[sb] = carry;
		a->size++;
	}

	if (sign)
		a->size *= -1;

	a->type = COMPOSITE;
	return;
}

void zShiftLeft(z *a, z *b, int x)
{	
	/* Computes a = b << x */
	int i,sign,wordshift;
	int y;
	int sb,j;
	fp_digit mask, carry, nextcarry;

	wordshift = x/BITS_PER_DIGIT;
	x = x%BITS_PER_DIGIT;

#if BITS_PER_DIGIT == 64
	mask = LEFTSHIFTMASKS64[x];
#else
	mask = LEFTSHIFTMASKS32[x];
#endif
	
	sb = abs(b->size);
	sign = (b->size < 0);
	a->size = sb;

	if (a->alloc < (sb + wordshift + 1))
		zGrow(a,sb + wordshift + LIMB_BLKSZ);

	//for each digit, remember the highest x bits using the mask, then shift.
	//the highest x bits becomes the lowest x bits for the next digit
	y = BITS_PER_DIGIT - x;
	carry = 0;
	for (j = 0; j < sb; ++j)
	{
		nextcarry = (b->val[j] & mask) >> y;
		a->val[j] = b->val[j] << x | carry;
		carry = nextcarry;
	}
	
	if (carry)
	{
		a->val[sb] = carry;
		a->size++;
	}

	if (wordshift)
	{
		//now shift by any full words
		for (i=a->size - 1;i>=0;i--)
			a->val[i+wordshift] = a->val[i];
		//zero out the ones that were shifted
		for (i=wordshift-1;i>=0;i--)
			a->val[i] = 0;
		a->size += wordshift;
	}	

	if (sign)
		a->size *= -1;

	a->type = COMPOSITE;
	return;
}

void zShiftLeft32(z32 *a, z32 *b, int x)
{	
	/* Computes a = b << x */
	int i,sign,wordshift;
	int y;
	int sb,j;
	fp_digit mask, carry, nextcarry;

	wordshift = x/32;
	x = x%32;

	mask = LEFTSHIFTMASKS32[x];
	
	sb = abs(b->size);
	sign = (b->size < 0);
	a->size = sb;

	if (a->alloc < (sb + wordshift + 1))
		zGrow32(a,sb + wordshift + LIMB_BLKSZ);

	//for each digit, remember the highest x bits using the mask, then shift.
	//the highest x bits becomes the lowest x bits for the next digit
	y = 32 - x;
	carry = 0;
	for (j = 0; j < sb; ++j)
	{
		nextcarry = (b->val[j] & mask) >> y;
		a->val[j] = b->val[j] << x | carry;
		carry = nextcarry;
	}
	
	if (carry)
	{
		a->val[sb] = carry;
		a->size++;
	}

	if (wordshift)
	{
		//now shift by any full words
		for (i=a->size - 1;i>=0;i--)
			a->val[i+wordshift] = a->val[i];
		//zero out the ones that were shifted
		for (i=wordshift-1;i>=0;i--)
			a->val[i] = 0;
		a->size += wordshift;
	}	

	if (sign)
		a->size *= -1;

	a->type = COMPOSITE;
	return;
}

static uint64 RIGHTSHIFTMASKS64[64] = {
	0x0000000000000000ULL, 0x0000000000000001ULL, 0x0000000000000003ULL, 0x0000000000000007ULL, 
	0x000000000000000FULL, 0x000000000000001FULL, 0x000000000000003FULL, 0x000000000000007FULL, 
	0x00000000000000FFULL, 0x00000000000001FFULL, 0x00000000000003FFULL, 0x00000000000007FFULL, 
	0x0000000000000FFFULL, 0x0000000000001FFFULL, 0x0000000000003FFFULL, 0x0000000000007FFFULL,
	0x000000000000FFFFULL, 0x000000000001FFFFULL, 0x000000000003FFFFULL, 0x000000000007FFFFULL, 
	0x00000000000FFFFFULL, 0x00000000001FFFFFULL, 0x00000000003FFFFFULL, 0x00000000007FFFFFULL, 
	0x0000000000FFFFFFULL, 0x0000000001FFFFFFULL, 0x0000000003FFFFFFULL, 0x0000000007FFFFFFULL, 
	0x000000000FFFFFFFULL, 0x000000001FFFFFFFULL, 0x000000003FFFFFFFULL, 0x000000007FFFFFFFULL,
	0x00000000FFFFFFFFULL, 0x00000001FFFFFFFFULL, 0x00000003FFFFFFFFULL, 0x00000007FFFFFFFFULL, 
	0x0000000FFFFFFFFFULL, 0x0000001FFFFFFFFFULL, 0x0000003FFFFFFFFFULL, 0x0000007FFFFFFFFFULL, 
	0x000000FFFFFFFFFFULL, 0x000001FFFFFFFFFFULL, 0x000003FFFFFFFFFFULL, 0x000007FFFFFFFFFFULL, 
	0x00000FFFFFFFFFFFULL, 0x00001FFFFFFFFFFFULL, 0x00003FFFFFFFFFFFULL, 0x00007FFFFFFFFFFFULL,
	0x0000FFFFFFFFFFFFULL, 0x0001FFFFFFFFFFFFULL, 0x0003FFFFFFFFFFFFULL, 0x0007FFFFFFFFFFFFULL, 
	0x000FFFFFFFFFFFFFULL, 0x001FFFFFFFFFFFFFULL, 0x003FFFFFFFFFFFFFULL, 0x007FFFFFFFFFFFFFULL, 
	0x00FFFFFFFFFFFFFFULL, 0x01FFFFFFFFFFFFFFULL, 0x03FFFFFFFFFFFFFFULL, 0x07FFFFFFFFFFFFFFULL, 
	0x0FFFFFFFFFFFFFFFULL, 0x1FFFFFFFFFFFFFFFULL, 0x3FFFFFFFFFFFFFFFULL, 0x7FFFFFFFFFFFFFFFULL};

static uint32 RIGHTSHIFTMASKS32[32] = {
	0x00000000, 0x00000001, 0x00000003, 0x00000007,
	0x0000000F, 0x0000001F, 0x0000003F, 0x0000007F,
	0x000000FF, 0x000001FF, 0x000003FF, 0x000007FF,
	0x00000FFF, 0x00001FFF, 0x00003FFF, 0x00007FFF,
	0x0000FFFF, 0x0001FFFF, 0x0003FFFF, 0x0007FFFF,
	0x000FFFFF, 0x001FFFFF, 0x003FFFFF, 0x007FFFFF,
	0x00FFFFFF, 0x01FFFFFF, 0x03FFFFFF, 0x07FFFFFF,
	0x0FFFFFFF, 0x1FFFFFFF, 0x3FFFFFFF, 0x7FFFFFFF};


void zShiftRight(z *a, z *b, int x)
{	/* Computes a = b >> x */
	int i, y, sign, wordshift;
	long sb;
	fp_digit mask, carry, nextcarry;
	//z tmp, *a;

	//test: likely will make this routine slower.  is it in a critical routine?
	//zInit(&tmp);
	//zCopy(b,&tmp);
	//a = &tmp;

	wordshift = x/BITS_PER_DIGIT;
	x = x%BITS_PER_DIGIT;

#if BITS_PER_DIGIT == 64
	mask = RIGHTSHIFTMASKS64[x];
#else
	mask = RIGHTSHIFTMASKS32[x];
#endif
	
	sb = b->size;
	sign = (b->size < 0);
	a->size = sb;

	if (a->alloc < (sb + 1))
		zGrow(a,sb + LIMB_BLKSZ);

	//for each digit, remember the lowest x bits using the mask, then shift.
	//the lowest x bits becomes the highest x bits for the next digit
	y = BITS_PER_DIGIT - x;
	carry = 0;
	for (i = sb - 1; i >= 0; --i)
	{
		nextcarry = (b->val[i] & mask) << y;
		a->val[i] = b->val[i] >> x | carry;
		carry = nextcarry;
	}

	if (wordshift)
	{
		//now shift by any full words
		for (i=0;i<a->size - 1;i++)
			a->val[i] = a->val[i+wordshift];
		//zero out the ones that were shifted
		a->size -= wordshift;
	}	

	for (i = a->size - 1; i >= 0; i--)
	{
		if (a->val[i] == 0)
			a->size--;
		else
			break;
	}

	if (a->size == 0)
		a->size = 1;

	if (sign)
		a->size *= -1;

	a->type = UNKNOWN;

	//zCopy(a,out);
	//zFree(&tmp);

	return;
}

void zShiftRight_x(z *a, z *b, int x)
{	/* Computes a = b >> x, where x is less than a word*/
	int i, y, sign;
	long sb;
	fp_digit mask, carry, nextcarry;

#if BITS_PER_DIGIT == 64
	mask = RIGHTSHIFTMASKS64[x];
#else
	mask = RIGHTSHIFTMASKS32[x];
#endif

	if (b->size < 0)
	{
		sb = -1 * b->size;
		sign = 1;
	}
	else
	{
		sb = b->size;
		sign = 0;
	}

	a->size = sb;

	//for each digit, remember the lowest x bits using the mask, then shift.
	//the lowest x bits becomes the highest x bits for the next digit
	y = BITS_PER_DIGIT - x;
	carry = 0;
	for (i = sb - 1; i >= 0; --i)
	{
		nextcarry = (b->val[i] & mask) << y;
		a->val[i] = b->val[i] >> x | carry;
		carry = nextcarry;
	}

	for (i = a->size - 1; i >= 0; i--)
	{
		if (a->val[i] == 0)
			a->size--;
		else
			break;
	}

	if (a->size == 0)
		a->size = 1;

	if (sign)
		a->size *= -1;

	a->type = UNKNOWN;

	return;
}


void zShiftRight_1(z *a, z *b)
{	/* Computes a = b >> 1 */
	int i, y, sign;
	long sb;
	fp_digit carry, nextcarry;

	if (b->size < 0)
	{
		sb = -1 * b->size;
		sign = 1;
	}
	else
	{
		sb = b->size;
		sign = 0;
	}

	//for each digit, remember the lowest x bits using the mask, then shift.
	//the lowest x bits becomes the highest x bits for the next digit
	y = BITS_PER_DIGIT - 1;
	carry = 0;
	for (i = sb - 1; i >= 0; --i)
	{
		nextcarry = (b->val[i] & 0x1) << y;
		a->val[i] = b->val[i] >> 1 | carry;
		carry = nextcarry;
	}

	for (i = a->size - 1; i >= 0; i--)
	{
		if (a->val[i] == 0)
			a->size--;
		else
			break;
	}

	if (a->size == 0)
		a->size = 1;

	if (sign)
		a->size *= -1;

	a->type = UNKNOWN;

	return;
}


void zShiftRight32(z32 *a, z32 *b, int x)
{	
	//Computes a = in >> 1 
	//only used in the trial division stage of QS

	int i, y, sign, wordshift;
	long sb;
	fp_digit mask, carry, nextcarry;
	//z tmp, *a;

	//test: likely will make this routine slower.  is it in a critical routine?
	//zInit(&tmp);
	//zCopy(b,&tmp);
	//a = &tmp;

	wordshift = x/32;
	x = x%32;

	mask = RIGHTSHIFTMASKS32[x];
	
	sb = b->size;
	sign = (b->size < 0);
	a->size = sb;

	if (a->alloc < (sb + 1))
		zGrow32(a,sb + LIMB_BLKSZ);

	//for each digit, remember the lowest x bits using the mask, then shift.
	//the lowest x bits becomes the highest x bits for the next digit
	y = 32 - x;
	carry = 0;
	for (i = sb - 1; i >= 0; --i)
	{
		nextcarry = (b->val[i] & mask) << y;
		a->val[i] = b->val[i] >> x | carry;
		carry = nextcarry;
	}

	if (wordshift)
	{
		//now shift by any full words
		for (i=0;i<a->size - 1;i++)
			a->val[i] = a->val[i+wordshift];
		//zero out the ones that were shifted
		a->size -= wordshift;
	}	

	for (i = a->size - 1; i >= 0; i--)
	{
		if (a->val[i] == 0)
			a->size--;
		else
			break;
	}

	if (a->size == 0)
		a->size = 1;

	if (sign)
		a->size *= -1;

	a->type = UNKNOWN;

	return;
}

void zShiftRight32_x(z32 *a, z32 *b, int x)
{	
	//Computes a = in >> x, where x is less than a word 
	//only used in the trial division stage of QS

	int i, y, sign;
	long sb;
	fp_digit mask, carry, nextcarry;

	mask = RIGHTSHIFTMASKS32[x];
	
	sb = b->size;
	sign = (b->size < 0);
	a->size = sb;

	//for each digit, remember the lowest x bits using the mask, then shift.
	//the lowest x bits becomes the highest x bits for the next digit
	y = 32 - x;
	carry = 0;
	for (i = sb - 1; i >= 0; --i)
	{
		nextcarry = (b->val[i] & mask) << y;
		a->val[i] = b->val[i] >> x | carry;
		carry = nextcarry;
	}

	for (i = a->size - 1; i >= 0; i--)
	{
		if (a->val[i] == 0)
			a->size--;
		else
			break;
	}

	if (a->size == 0)
		a->size = 1;

	if (sign)
		a->size *= -1;

	a->type = UNKNOWN;

	return;
}

#if defined(MSC_ASM32A) && !defined(ASM_ARITH_DEBUG)

void spAdd(fp_digit u, fp_digit v, fp_digit *sum, fp_digit *carry)
{
	fp_digit s,c;

	s = v;
	c = 0;
	ASM_M {
		mov eax, u
		add s, eax
		adc c, 0
	}

	*sum = s;
	*carry = c;

	return;
}

void spAdd3(fp_digit u, fp_digit v, fp_digit w, fp_digit *sum, fp_digit *carry)
{
	fp_digit s,c;

	s = v;
	c = 0;

	ASM_M {
		mov eax, u
		add eax, w
		adc c, 0
		add s, eax
		adc c, 0
	}

	*sum = s;
	*carry = c;

	return;
}

void spSub3(fp_digit u, fp_digit v, fp_digit w, fp_digit *sub, fp_digit *borrow)
{
	fp_digit s,b;

	s = v;
	b = 0;

	ASM_M {
		mov eax, u
		sub eax, s
		adc b, 0
		sub eax, w
		adc b, 0
		mov s, eax
	}

	*sub = s;
	*borrow = b;
 
	return;
}

void spSub(fp_digit u, fp_digit v, fp_digit *sub, fp_digit *borrow)
{
	fp_digit s,b;

	s = v;
	b = 0;

	ASM_M {
		mov eax, u
		sub eax, s
		adc b, 0
		mov s, eax
	}

	*sub = s;
	*borrow = b;
 
	return;
}

fp_digit spDivide(fp_digit *q, fp_digit *r, fp_digit u[2], fp_digit v)
{
	fp_digit tmpq, tmpr;
	tmpr = *r = u[1];
	tmpq = *q = u[0];

	ASM_M {
		mov eax, tmpq
		mov edx, tmpr
		div v
		mov tmpr, edx
		mov tmpq, eax
	}

	*q = tmpq;
	*r = tmpr;

	return 0;
}

void spMultiply(fp_digit u, fp_digit v, fp_digit *product, fp_digit *carry)
{
	ASM_M {
		mov eax, u
		mul v
		mov v, eax
		mov u, edx
	}

	*product = v;
	*carry = u;

	return;
}

void spMulAdd(fp_digit u, fp_digit v, fp_digit w, fp_digit t, fp_digit *lower, fp_digit *carry)
{
	//u*v + (w+t)  used a lot in multiplication
	//fp_word uu;
	//uu = (fp_word)u * (fp_word)v + (fp_word)w + (fp_word)t;
	//*lower = (fp_digit)uu;
	//*carry = (fp_digit)(uu >> BITS_PER_DIGIT);
	fp_digit tmp;

	spMultiply(u,v,lower,carry);
	spAdd3(*lower,w,t,lower,&tmp);
	spAdd(tmp,*carry,carry,&tmp);
	return;
}

void spMulMod(fp_digit u, fp_digit v, fp_digit m, fp_digit *w)
{
	fp_digit p[2];
	fp_digit q;

	spMultiply(u,v,&p[0],&p[1]);
	spDivide(&q,w,p,m);

	return;
}

#elif defined(GCC_ASM32A) && !defined(ASM_ARITH_DEBUG)

void spAdd(fp_digit u, fp_digit v, fp_digit *sum, fp_digit *carry)
{
	fp_word s,c;

	s = v;
	c = 0;

	ASM_G ("movl %2, %%eax		\n\t"
		"addl %%eax, %3		\n\t"
		"adcl $0, %4		\n\t"
		: "=r"(s), "=r"(c)
		: "r"(u), "0"(s), "1"(c)
		: "eax", "memory", "cc");

	*sum = s;
	*carry = c;

	return;
}

void spAdd3(fp_digit u, fp_digit v, fp_digit w, fp_digit *sum, fp_digit *carry)
{
	fp_word s,c;

	s = v;
	c = 0;

	ASM_G ("movl %2, %%eax		\n\t"
		"addl %3, %%eax		\n\t"
		"adcl $0, %5		\n\t"
		"addl %%eax, %4		\n\t"
		"adcl $0, %5		\n\t"
		: "=r"(s), "=r"(c)
		: "r"(u), "r"(w), "0"(s), "1"(c)
		: "eax", "memory", "cc");

	*sum = s;
	*carry = c;

	return;
}

void spSub3(fp_digit u, fp_digit v, fp_digit w, fp_digit *sub, fp_digit *borrow)
{
	fp_word s,b;

	s = v;
	b = 0;

	ASM_G ("movl %2, %%eax		\n\t"
		"subl %4, %%eax		\n\t"
		"adcl $0, %5		\n\t"
		"subl %3, %%eax		\n\t"
		"adcl $0, %5		\n\t"
		"movl %%eax, %4		\n\t"
		: "=r"(s), "=r"(b)
		: "r"(u), "r"(w), "0"(s), "1"(b)
		: "eax", "memory", "cc");

	*sub = s;
	*borrow = b;
 
	return;
}

void spSub(fp_digit u, fp_digit v, fp_digit *sub, fp_digit *borrow)
{
	fp_word s,b;

	s = v;
	b = 0;

	ASM_G ("movl %2, %%eax		\n\t"
		"subl %3, %%eax		\n\t"
		"adcl $0, %4		\n\t"
		"movl %%eax, %3		\n\t"
		: "=r"(s), "=r"(b)
		: "r"(u), "0"(s), "1"(b)
		: "eax", "memory", "cc");

	*sub = s;
	*borrow = b;
 
	return;
}


fp_digit spDivide(fp_digit *q, fp_digit *r, fp_digit u[2], fp_digit v)
{
	*r = u[1];
	*q = u[0];
	ASM_G ("divl %4"
		: "=a"(*q),"=d"(*r)
		: "1"(*r), "0"(*q), "r"(v) );

	return 0;
}

void spMultiply(fp_digit u, fp_digit v, fp_digit *product, fp_digit *carry)
{
	*product = v;
	*carry = u;

	ASM_G ("movl %2, %%eax	\n\t"
		"mull %3	\n\t"
		"movl %%eax, %0		\n\t"
		"movl %%edx, %1		\n\t"
		: "=r"(*product), "=r"(*carry)
		: "1"(*carry), "0"(*product)
		: "eax", "%rdx", "cc");

	return;
}

void spMulAdd(fp_digit u, fp_digit v, fp_digit w, fp_digit t, fp_digit *lower, fp_digit *carry)
{
	fp_digit k,p;
	spMultiply(u,v,&p,carry);
	spAdd3(p,w,t,lower,&k);
	*carry += k;
	return;
}

void spMulMod(fp_digit u, fp_digit v, fp_digit m, fp_digit *w)
{
	fp_digit p[2];
	fp_digit q;

	spMultiply(u,v,&p[0],&p[1]);
	spDivide(&q,w,p,m);

	return;
}

#elif defined(GCC_ASM64X) && !defined(ASM_ARITH_DEBUG)

void spAdd(fp_digit u, fp_digit v, fp_digit *sum, fp_digit *carry)
{
	//fp_word s,c;
	uint64 s,c;

	s = v;
	c = 0;

	ASM_G ("movq %2, %%rax		\n\t"
		"addq %%rax, %3		\n\t"
		"adcq $0, %4		\n\t"
		: "=r"(s), "=r"(c)
		: "r"(u), "0"(s), "1"(c)
		: "rax", "memory", "cc");

	*sum = s;
	*carry = c;

	return;
}

void spAdd3(fp_digit u, fp_digit v, fp_digit w, fp_digit *sum, fp_digit *carry)
{
	//fp_word s,c;
	uint64 s,c;

	s = v;
	c = 0;

	ASM_G ("movq %2, %%rax		\n\t"
		"addq %3, %%rax		\n\t"
		"adcq $0, %5		\n\t"
		"addq %%rax, %4		\n\t"
		"adcq $0, %5		\n\t"
		: "=r"(s), "=r"(c)
		: "r"(u), "r"(w), "0"(s), "1"(c)
		: "rax", "memory", "cc");

	*sum = s;
	*carry = c;

	return;
}

void spSub3(fp_digit u, fp_digit v, fp_digit w, fp_digit *sub, fp_digit *borrow)
{
	//fp_word s,b;
	uint64 s,b;

	s = v;
	b = 0;

	ASM_G ("movq %2, %%rax		\n\t"
		"subq %4, %%rax		\n\t"
		"adcq $0, %5		\n\t"
		"subq %3, %%rax		\n\t"
		"adcq $0, %5		\n\t"
		"movq %%rax, %4		\n\t"
		: "=r"(s), "=r"(b)
		: "r"(u), "r"(w), "0"(s), "1"(b)
		: "rax", "memory", "cc");

	*sub = s;
	*borrow = b;
 
	return;
}

void spSub(fp_digit u, fp_digit v, fp_digit *sub, fp_digit *borrow)
{
	//fp_word s,b;
	uint64 s,b;

	s = v;
	b = 0;

	ASM_G ("movq %2, %%rax		\n\t"
		"subq %3, %%rax		\n\t"
		"adcq $0, %4		\n\t"
		"movq %%rax, %3		\n\t"
		: "=r"(s), "=r"(b)
		: "r"(u), "0"(s), "1"(b)
		: "rax", "memory", "cc");

	*sub = s;
	*borrow = b;
 
	return;
}

fp_digit spDivide(fp_digit *q, fp_digit *r, fp_digit u[2], fp_digit v)
{
	*r = u[1];
	*q = u[0];
	ASM_G ("divq %4"
		: "=a"(*q),"=d"(*r)
		: "1"(*r), "0"(*q), "r"(v) );

	return 0;
}

void spMultiply(fp_digit u, fp_digit v, fp_digit *product, fp_digit *carry)
{
	*product = v;
	*carry = u;

	ASM_G ("movq %2, %%rax	\n\t"
		"mulq %3	\n\t"
		"movq %%rax, %0		\n\t"
		"movq %%rdx, %1		\n\t"
		: "=r"(*product), "=r"(*carry)
		: "1"(*carry), "0"(*product)
		: "rax", "rdx", "cc");

	return;
}

uint64 spPRP2(uint64 p)
{
	// do a base-2 prp test on the input, where p is greater than 2^32
	// i.e., compute 2^(p-1) % p.
	// since p is more than 32 bits we can do the accumulation division 
	// free for the first 5 iterations.  may not be much, but it's something.
	uint64 result;

	ASM_G (
		"xorq	%%rbx, %%rbx \n\t"
		"xorq	%%rdi, %%rdi \n\t"
		"addq	$1, %%rdi \n\t"		/* n = 1 */
		"0:	\n\t"					/* begin loop */					
		"test	$1, %%rcx \n\t"		/* exp & 0x1 */
		"je	2f		\n\t"			/* bit not set, skip accumulation into n */
		"movq	%%rax, %%rsi \n\t"	/* save acc */
		"mulq	%%rdi \n\t"			/* n * acc mod m */
		"movq	%%rax, %%rdi \n\t"	/* save n */
		"movq	%%rsi, %%rax \n\t"	/* restore acc */
		"2:			\n\t"			/* square acc stage */
		"shrq	$1, %%rcx \n\t"		/* base >>= 1 */
		"addq	$1, %%rbx \n\t"
		"mulq	%%rax \n\t"			/* acc = acc * acc*/
		"cmpq	$5, %%rbx \n\t"		/* 5 iterations? */
		"jb 0b \n\t"
		"3:	\n\t"					/* begin loop */					
		"test	$1, %%rcx \n\t"		/* exp & 0x1 */
		"je	4f		\n\t"			/* bit not set, skip accumulation into n */
		"movq	%%rax, %%rsi \n\t"	/* save acc */
		"mulq	%%rdi \n\t"			/* n * acc mod m */
		"divq	%3 \n\t"		
		"movq	%%rdx, %%rdi \n\t"	/* save n */
		"movq	%%rsi, %%rax \n\t"	/* restore acc */
		"4:			\n\t"			/* square acc stage */
		"shrq	$1, %%rcx \n\t"		/* base >>= 1 */
		"mulq	%%rax \n\t"			/* acc = acc * acc*/
		"divq	%3 \n\t"
		"cmpq	$0, %%rcx \n\t"		/* exp == 0? */
		"movq	%%rdx, %%rax \n\t"	/* mod m */
		"jne 3b \n\t"
		"movq	%%rdi, %0 \n\t"
		: "=r"(result)
		: "a"(2), "c"(p-1), "r"(p)
		: "rbx", "rdx", "rdi", "rsi", "cc");				/* return result */

	return result;

}

uint64 spModExp_asm(uint64 b, uint64 e, uint64 m)
{
	uint64 result;

	ASM_G (
		"xorq	%%rdi, %%rdi \n\t"
		"addq	$1, %%rdi \n\t"		/* n = 1 */
		"cmpq	$0, %%rcx \n\t"		/* exp == 0? */
		"je 1f \n\t"
		"0:	\n\t"					/* begin loop */					
		"test	$1, %%rcx \n\t"		/* exp & 0x1 */
		"je	2f		\n\t"			/* bit not set, skip accumulation into n */
		"movq	%%rax, %%rsi \n\t"	/* save acc */
		"mulq	%%rdi \n\t"			/* n * acc mod m */
		"divq	%3 \n\t"		
		"movq	%%rdx, %%rdi \n\t"	/* save n */
		"movq	%%rsi, %%rax \n\t"	/* restore acc */
		"2:			\n\t"			/* square acc stage */
		"shrq	$1, %%rcx \n\t"		/* base >>= 1 */
		"mulq	%%rax \n\t"			/* acc = acc * acc*/
		"divq	%3 \n\t"
		"cmpq	$0, %%rcx \n\t"		/* exp == 0? */
		"movq	%%rdx, %%rax \n\t"	/* mod m */
		"jne 0b \n\t"
		"1:			\n\t"			/* end loop */
		"movq	%%rdi, %0 \n\t"
		: "=r"(result)
		: "a"(b), "c"(e), "r"(m)
		: "rdx", "rdi", "rsi", "cc");				/* return result */

	return result;
}

void spMulAdd(fp_digit u, fp_digit v, fp_digit w, fp_digit t, fp_digit *lower, fp_digit *carry)
{
	fp_digit k,p;
	spMultiply(u,v,&p,carry);
	spAdd3(p,w,t,lower,&k);
	*carry += k;
	return;
}

void spMulMod(fp_digit u, fp_digit v, fp_digit m, fp_digit *w)
{
	fp_digit p[2];
	fp_digit q;

	spMultiply(u,v,&p[0],&p[1]);
	spDivide(&q,w,p,m);

	return;
}

#elif BITS_PER_DIGIT == 32

fp_digit spDivide(fp_digit *q, fp_digit *r, fp_digit u[2], fp_digit v)
{
	fp_word uu,qq;
	uu = (fp_word)u[1] << BITS_PER_DIGIT | u[0];
	qq = (uu/(fp_word)v);
	*r = (fp_digit)(uu%(fp_word)v);
	*q = (fp_digit)qq;
	return (fp_digit)(qq >> BITS_PER_DIGIT);
}

void spMultiply(fp_digit u, fp_digit v, fp_digit *product, fp_digit *carry)
{
	fp_word uu;
	uu = (fp_word)u * (fp_word)v;
	*product = (fp_digit)uu;
	*carry = (fp_digit)(uu >> BITS_PER_DIGIT);
	return;
}

void spMulAdd(fp_digit u, fp_digit v, fp_digit w, fp_digit t, fp_digit *lower, fp_digit *carry)
{
	//u*v + (w+t)  used a lot in multiplication
	//fp_word uu;
	//uu = (fp_word)u * (fp_word)v + (fp_word)w + (fp_word)t;
	//*lower = (fp_digit)uu;
	//*carry = (fp_digit)(uu >> BITS_PER_DIGIT);
	fp_digit tmp;

	spMultiply(u,v,lower,carry);
	spAdd3(*lower,w,t,lower,&tmp);
	spAdd(tmp,*carry,carry,&tmp);
	return;
}

void spAdd(fp_digit u, fp_digit v, fp_digit *sum, fp_digit *carry)
{
	fp_word w;
	w = (fp_word)u + (fp_word)v;
	*sum = (fp_digit)w;
	*carry = (fp_digit)(w >> BITS_PER_DIGIT);
	return;
}

void spAdd3(fp_digit u, fp_digit v, fp_digit w, fp_digit *sum, fp_digit *carry)
{
	fp_word uu;
	uu = (fp_word)u + (fp_word)v + (fp_word)w;
	*sum = (fp_digit)uu ;
	*carry = (fp_digit)(uu >> BITS_PER_DIGIT);
	return;
}

void spSub3(fp_digit u, fp_digit v, fp_digit w, fp_digit *sub, fp_digit *borrow)
{
	fp_word uu;
	uu = (fp_word)u - (fp_word)v - (fp_word)w;
	*sub = (fp_digit)uu ;
	*borrow = (fp_digit)((uu >> BITS_PER_DIGIT) && MAX_DIGIT);
	return;
}

void spSub(fp_digit u, fp_digit v, fp_digit *sub, fp_digit *borrow)
{
	fp_word uu;
	uu = (fp_word)u - (fp_word)v;
	*sub = (fp_digit)uu ;
	*borrow = (fp_digit)((uu >> BITS_PER_DIGIT) && MAX_DIGIT);
	return;
}

void spMulMod(fp_digit u, fp_digit v, fp_digit m, fp_digit *w)
{
	fp_digit p[2];
	fp_digit q;

	spMultiply(u,v,&p[0],&p[1]);
	spDivide(&q,w,p,m);

	return;
}

#elif BITS_PER_DIGIT == 64

// TBD: needs low level assembly routines similar to mod_64

fp_digit spDivide(fp_digit *q, fp_digit *r, fp_digit u[2], fp_digit v)
{

	
}

void spMultiply(fp_digit u, fp_digit v, fp_digit *product, fp_digit *carry)
{

	
	return;
}

void spMulAdd(fp_digit u, fp_digit v, fp_digit w, fp_digit t, fp_digit *lower, fp_digit *carry)
{

	return;
}

void spAdd(fp_digit u, fp_digit v, fp_digit *sum, fp_digit *carry)
{
	
	return;
}

void spAdd3(fp_digit u, fp_digit v, fp_digit w, fp_digit *sum, fp_digit *carry)
{

	return;
}

void spSub3(fp_digit u, fp_digit v, fp_digit w, fp_digit *sub, fp_digit *borrow)
{
	
	return;
}

void spSub(fp_digit u, fp_digit v, fp_digit *sub, fp_digit *borrow)
{
	
	return;
}

void spMulMod(fp_digit u, fp_digit v, fp_digit m, fp_digit *w)
{

	return;
}

#endif



