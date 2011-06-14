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
#include "monty.h"
#include "tfm.h"

/*
implements routines to perform computations with 
montgomery arithmatic
*/

static int TFM_MONTY;
void monty_mul_interleaved(z *a, z *b, z *c, z *n);

void zREDC(z *T, z *n)
{
	/* from handbook of applied cryptography, ch. 14
	INPUT: integers m = (mn-1 . . .m1m0)b with gcd(m; b) = 1, R = b^n,m' = -m^-1 mod
	b, and T = (t2n-1 . . . t1t0)b <mR.
	OUTPUT: TR^-1 mod m, the reduction of T mod m in montgomery representation...
	1. A=T . (Notation: A = (a2n-1 . . . a1a0)b.)
	2. For i from 0 to (n - 1) do the following:
	2.1 ui=ai*m' mod b.
	2.2 A=A + ui*m*b^i.
	3. A=A/b^n.
	4. If A > m then A=A-m.
	5. Return(A).
	*/
	int i,j,ix,su;
	fp_digit nhat = montyconst.nhat.val[0], ui,k;
	z mtmp3;

	if (TFM_MONTY == 1)
	{
		fp_montgomery_reduce(T,n,montyconst.nhat.val[0]);
		return;
	}
	
	//printf("shouldn't get to here\n");
	zInit(&mtmp3);

	if (mtmp3.alloc < n->size * 2)
		zGrow(&mtmp3,n->size * 2);

	//T needs to have allocated montyconst.n.size + T.size
	if (T->alloc < n->size + T->size)
		zGrow(T,n->size + T->size + 1);
	
	for (i=0;i<n->size;i++)
	{
		//the mod b happens automatically because only the 
		//lower 32 bits of the product is returned.
		ui = T->val[i] * nhat;						//ui = a1*nhat mod b 	
		//zShortMul(&montyconst.n,ui,&mtmp3);			//t1 = ui * n
		
		//short mul
		k=0;
		su = n->size;
		for (ix=0;ix<su;++ix)
			spMulAdd(n->val[ix],ui,0,k,&mtmp3.val[ix],&k);

		//if still have a carry, add a digit to w
		if (k)
		{
			mtmp3.val[su]=k;
			su++;
		}
		
		//check for significant digits.  only necessary if v or u = 0?
		for (ix = su - 1;ix>=0;--ix)
		{
			if (mtmp3.val[ix] != 0)
				break;
		}
		mtmp3.size = ix+1;
		

		for (j=mtmp3.size - 1;j>=0;j--)				//t1 *= b^i 
			mtmp3.val[j+i] = mtmp3.val[j];
		mtmp3.size += i;
		zAdd(T,&mtmp3,T);								//A += t1
	}

	for (j=0; j<T->size; j++)						//A /= b^n
		T->val[j] = T->val[j+n->size];
	T->size -= n->size;

	if (zCompare(T,n) > 0)				//if A > n, A = A-n
		zSub(T,n,T);

	if (T->size == 0)
		zCopy(n,T);

	zFree(&mtmp3);

	return;	
}

void to_monty(z *x, z *n)
{
	//given a number x in normal (hexadecimal) representation, 
	//find its montgomery representation

	//this uses some precomputed monty constants
	//xhat = (x * r) mod n
	z t1,t2;
	zInit(&t1);
	zInit(&t2);

	zMul(x,&montyconst.r,&t1);
	zDiv(&t1,n,&t2,x);

	zFree(&t1);
	zFree(&t2);

	return;
}

void monty_init(z *n)
{
	//for a input modulus n, initialize constants for 
	//montogomery representation
	//this assumes that n is relatively prime to 2, i.e. is odd.
	z g, b, q, r;

	//global montyconst structure
	zInit(&montyconst.nhat);
	zInit(&montyconst.r);
	zInit(&montyconst.rhat);
	zInit(&montyconst.one);

	
	if (abs(n->size) <= 16) 
	{
		fp_montgomery_setup(n,&montyconst.nhat.val[0]);
		fp_montgomery_calc_normalization(&montyconst.r,n);
		montyconst.one.val[0] = 1;
		montyconst.one.size = 1;
		to_monty(&montyconst.one,n);
		TFM_MONTY = 1;
		return;
	}
	else
		TFM_MONTY = 0;

	zInit(&g);
	zInit(&b);
	zInit(&q);
	zInit(&r);

	b.val[1]=1; b.size=2;

	//find r = b^t > N, where b = 2 ^32
	if (montyconst.r.alloc < n->size + 1)
		zGrow(&montyconst.r,n->size + 1);

	zClear(&montyconst.r);
	montyconst.r.size = n->size + 1;
	montyconst.r.val[montyconst.r.size - 1] = 1;

	//find nhat = -n^-1 mod b
	//nhat = -(n^-1 mod b) mod b = b - n^-1 mod b
	//since b is 2^32, this can be simplified, and made faster.
	xGCD(n,&b,&montyconst.nhat,&montyconst.rhat,&g);
	zSub(&b,&montyconst.nhat,&q);
	zCopy(&q,&montyconst.nhat);

	zCopy(&zOne,&montyconst.one);
	to_monty(&montyconst.one,n);

	zFree(&g);
	zFree(&b);
	zFree(&q);
	zFree(&r);
	return;
}

void monty_add(z *u, z *v, z *w, z *n)
{
	//add two numbers in the montgomery representation, returning their
	//sum in montgomery representation

	//work with the currently defined montyconst

	zAdd(u,v,w);
	if (zCompare(w,n) >= 0)
		zSub(w,n,w);

	return;
}

void monty_mul(z *u, z *v, z *w, z *n)
{
	//multiply two numbers in the montgomery representation, returning their
	//product in montgomery representation

	//work with the currently defined montyconst
	zMul(u,v,w);
	zREDC(w,n);
	//monty_mul_interleaved(u,v,w);

	return;
}


void monty_mul_interleaved(z *a, z *b, z *c, z *n)
{

	fp_digit nhat = montyconst.nhat.val[0], u;
	int i,j,t=n->size;
	int szb = abs(b->size);
	fp_digit k;
	z *t1,*t2;
	z s1,s2;

	zInit(&s1);
	zInit(&s2);
	t1 = &s1;
	t2 = &s2;
	zClear(t1);
	zClear(t2);

	for (i=0;i<t;i++)
	{
		u = (t1->val[0] + a->val[i] * b->val[0]) * nhat;	//truncation will provide mod b
		
		/****** short mul of b with ai, simultaneous with addition of A (in t1) ********/
		for (j=t1->size;j<szb;j++)
			t1->val[j] = 0;		//zero any unused words up to size of b, so we can add
		//mul and add up to size of b
		k=0;
		for (j=0;j<szb ;j++)
			spMulAdd(b->val[j],a->val[i],t1->val[j],k,t2->val + j,&k);
		//continue with add if A has more words
		for (;j<t1->size;j++)
			spAdd(t1->val[j],k,t2->val+j,&k);

		//adjust size
		if (t1->size > szb)
			t2->size = t1->size;
		else
			t2->size = szb;

		//account for carry
		if (k)
		{
			t2->val[t2->size]=k;
			t2->size++;
			j++;
		}
		/****** short mul of b with ai, simultaneous with addition of A (in t1) ********/


		/****** short mul of n with u, simultaneous with add. of prev step (in t2) 
		 and with right shift of one word                                       ********/

		for (;j<t;j++)
			t2->val[j] = 0;		//zero any unused words up to size of n, so we can add
		//mul and add up to size of n, store into one word previous
		k=0;
		//needs first mul to get k set right, answer gets shifted to oblivion
		spMulAdd(n->val[0],u,t2->val[0],k,t1->val,&k);
		for (j=1;j<t;j++)
			spMulAdd(n->val[j],u,t2->val[j],k,t1->val + j - 1,&k);
		//continue if t2 is bigger than n
		for (;j<t2->size;j++)
			spAdd(t2->val[j],k,t1->val+j-1,&k);

		//adjust size
		if (t2->size > t)
			t1->size = t2->size - 1;
		else
			t1->size = t - 1;

		//account for carry
		if (k)
		{
			t1->val[t1->size]=k;
			t1->size++;
		}
		/****** short mul of n with u, simultaneous with add. of prev step (in t2) 
		 and with right shift of one word                                       ********/

	}

	//almost done
	if (zCompare(t1,n) >= 0)
		zSub(t1,n,c);
	else
		zCopy(t1,c);

	zFree(&s1);
	zFree(&s2);
	return;
}

void monty_sqr(z *x, z *w, z *n)
{
	//square a number in the montgomery representation, returning their
	//product in montgomery representation

	//work with the currently defined montyconst
	zSqr(x,w);
	zREDC(w,n);

	return;
}

void monty_sub(z *u, z *v, z *w, z *n)
{
	//subtract two numbers in the montgomery representation, returning their
	//difference in montgomery representation

	//work with the currently defined montyconst
	zSub(u,v,w);
	if (w->size < 0)
	{
		w->size *= -1;
		zSub(n,w,w);
	}

	return;
}

void monty_free(void)
{
	zFree(&montyconst.nhat);
	zFree(&montyconst.r);
	zFree(&montyconst.rhat);
	zFree(&montyconst.one);

	return;
}

void zmModExp(z *a, z *b, z *u, z *nn)
{
	//computes a^b mod m = u using the right to left binary method
	//see, for instance, the handbook of applied cryptography
	//uses monty arith
	//a is already in monty rep, b doesn't need to be.
	z n,bb,aa,t;

	zInit(&aa);
	zInit(&bb);
	zInit(&n);
	zInit(&t);

	//overflow possibilities:
	//t ranges to 2x input 'a'
	//u needs at least as much space as modulus

	zCopy(&montyconst.one,&n);
	zCopy(a,&aa);

	zCopy(b,&bb);
	while (!isZero(&bb))
	{
		if (bb.val[0] & 0x1)
		{
			monty_mul(&n,&aa,&t,nn);
			zCopy(&t,&n);
		}
		zShiftRight(&bb,&bb,1);   //compute successive squares of a
		monty_sqr(&aa,&t,nn);
		zCopy(&t,&aa);
		if (aa.size < 0)
			aa.size *= -1;
	}
	zCopy(&n,u);

	zFree(&aa);
	zFree(&bb);
	zFree(&n);
	zFree(&t);
	return;
}

void zmModExpw(z *a, z *e, z *u, z *n, int k)
{
	//computes a^e mod m = u using the sliding window left to right binary method
	//see, for instance, the handbook of applied cryptography
	//uses monty arith
	//a is already in monty rep, b doesn't need to be.  k is the window size

	/*
	INPUT: g, e = (etet-1 . . . e1e0)2 with et = 1, and an integer k >= 1.
	OUTPUT: g^e.
	1. Precomputation.
	1.1 g1 = g, g2 = g^2.
	1.2 For i from 1 to (2^(k-1) - 1) do: g_{2i+1} =  g_{2i-1} * g2.
	2. A = 1, i = t.
	3. While i >= 0 do the following:
	3.1 If ei = 0 then do: A = A^2, i = i - 1.
	3.2 Otherwise (ei != 0), find the longest bitstring eiei-1 . . . el such that i-l+1 <= k
	and el = 1, and do the following:
	A = A^{2^{i-l+1}} * g_{eiei-1...el}2 , i = l - 1.
	4. Return(A).

	test -> 11749.  3 multiplications at i=7,4,0
	*/

	//need to allocate (2^(k-1) + 1) g's for precomputation.
	z *g, g2, ztmp;
	int numg, i, j, l, t, tmp1, tmp2;
	fp_digit utmp1;
	uint8 *bitarray;

	//overflow possibilities:
	//t ranges to 2x input 'a'
	//u needs at least as much space as modulus

	numg = (int)((1<<(k-1))+1);
	g = (z *)malloc(numg*sizeof(z));
	for (i=0;i<numg;i++)
		zInit(&g[i]);
	zInit(&g2);
	zInit(&ztmp);

	//precomputation
	zCopy(a,&g[0]);						//g[0] = a
	monty_sqr(a,&g2,n);					//g2 = a^2

	for (i=1;i<numg;i++)
		monty_mul(&g[i-1],&g2,&g[i],n);	//g[i] = g[i-1] * g2, where g[i] holds g^{2*i+1}

	zCopy(&montyconst.one,u);
	t = zBits(e);

	bitarray = (uint8 *)malloc(t * sizeof(uint8));
	//get e in one array
	for (i=0;i< e->size - 1;i++)
	{
		utmp1 = e->val[i];
		j=0;
		while (j<BITS_PER_DIGIT)
		{
			bitarray[BITS_PER_DIGIT*i+j] = (uint8)(utmp1 & 0x1);
			utmp1 >>= 1;
			j++;
		}
	}
	utmp1 = e->val[i];
	j=0;
	while (utmp1)
	{
		bitarray[BITS_PER_DIGIT*i+j] = (uint8)(utmp1 & 0x1);
		utmp1 >>= 1;
		j++;
	}

	i=t-1;
	while (i >= 0)
	{
		if (bitarray[i])
		{
			//find the longest bitstring ei,e1-1,...el such that i-l+1 <= k and el == 1
			l=i;
			if (i >= (k-1))
			{
				//protect against accessing bitarray past its boundaries
				for (j=k-1;j>0;j--)
				{
					if (bitarray[i-j])
					{
						//this is the longest possible string, exit
						l=i-j;
						break;
					}
				}
			}
			//now, bitarray[i] to bitarray[i-j] is the longest bitstring
			//figure out the g value to use corresponding to this bitstring
			tmp1 = 1;
			tmp2 = 0;
			for (j=l;j<=i;j++)
			{
				tmp2 += tmp1 * bitarray[j];
				tmp1 <<= 1;
			}
			tmp2 = (tmp2-1)/2;

			//do the operation A = A^{2^{i-l+1}} * g_{eiei-1...el}2
			for (j=0;j<(i-l+1);j++)
			{
				monty_sqr(u,&ztmp,n);
				zCopy(&ztmp,u);
			}
			monty_mul(u,&g[tmp2],&ztmp,n);
			zCopy(&ztmp,u);

			//decrement bit pointer
			i = l-1;
		}
		else
		{
			monty_sqr(u,&ztmp,n);
			zCopy(&ztmp,u);
			i--;
		}
	}

	for (i=0;i<numg;i++)
		zFree(&g[i]);
	free(g);
	zFree(&g2);
	zFree(&ztmp);
	free(bitarray);
	return;
}
