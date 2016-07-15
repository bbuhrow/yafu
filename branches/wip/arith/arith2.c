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
implements some slightly less basic arithmetic operations
with a fundamental base of U32.  U64 types are used
for the operations that require a carry.
*/

#include "yafu.h"
#include "soe.h"
#include "arith.h"
#include "gmp_xface.h"

void zShanksTonelli(z *a, fp_digit p, fp_digit *sq)
{
	//a is a quadratic residue mod p
	//p is an odd prime
	//find x where x^2 == a mod p
	//we assume p will always fit into an fp_digit, therefore x will as well.
	//see paper by Ezra Brown
	fp_digit x=0,b=0,g=0,n=0,s=0,r=0,e=0,b2m=0,tmp=0;
	int i;
	z t1;

	zInit(&t1);

	//factor p-1 = Q*2^S, where Q is odd and S >= 1.
	s = p-1;
	e=0;
	while (!(s & 1))
	{
		s >>= 1;
		e++;
	}

	//find a quadratic non-residue mod p.  keep it small to reduce the work of modexp
	n = 3;
	while (1)
	{
		if (jacobi_1(n,p) < 0)
			break;
		n++;
	}

	//approximate the root x = a^[(s+1)/2] mod p
	zModExp_1(a,(s+1)/2,p,&x);
	
	//guess at fudge factor b = a^s
	zModExp_1(a,s,p,&b);

	//initialize g = n^s
	sp2z(n,&t1);
	zModExp_1(&t1,s,p,&g);

	//initialize r = e
	r = e;

	while (1)
	{
		//find m such that b^(2^m) == 1 mod p with m between 0 and r-1
		b2m = b;
		for (i=0;i<(int)r;i++)
		{
			if (b2m==1) break;

			//successivly square b mod p
			spMulMod(b2m,b2m,p,&b2m);
		}

		if (i == 0)
		{
			*sq = x;
			goto free;
		}

		//replace x by x*g^(2^(r-m-1))
		sp2z(g,&t1);
		//zModExp_1(&t1,(fp_digit)pow(2,(int)r-i-1),p,&tmp);
		zModExp_1(&t1,1<<(r-i-1),p,&tmp);
		spMulMod(tmp,x,p,&x);

		//replace g by g^(2^(r-m)) and
		//replace b by b*g
		zModExp_1(&t1,1<<(r-i),p,&g);
		spMulMod(g,b,p,&b);
		
		r = i;
	}
		
free:
	//return the smallest solution always
	if (*sq > (p>>1)) 
		*sq = p - *sq;

	zFree(&t1);
	return;
}

void ShanksTonelli_1(fp_digit a, fp_digit p, fp_digit *sq)
{
	//a is a quadratic residue mod p
	//p is an odd prime
	//find x where x^2 == a mod p
	//we assume p will always fit into an fp_digit, therefore x will as well.
	//see paper by Ezra Brown
	fp_digit x=0,b=0,g=0,n=0,s=0,r=0,e=0,b2m=0,tmp=0;
	int i;

	//factor p-1 = Q*2^S, where Q is odd and S >= 1.
	s = p-1;
	e=0;
	while (!(s & 1))
	{
		s >>= 1;
		e++;
	}

	//find a quadratic non-residue mod p.  keep it small to reduce the work of modexp
	n = 3;
	while (1)
	{
		if (jacobi_1(n,p) < 0)
			break;
		n++;
	}

	//approximate the root x = a^[(s+1)/2] mod p
	spModExp(a,(s+1)/2,p,&x);
	
	//guess at fudge factor b = a^s
	spModExp(a,s,p,&b);

	//initialize g = n^s
	spModExp(n,s,p,&g);

	//initialize r = e
	r = e;

	while (1)
	{
		//find m such that b^(2^m) == 1 mod p with m between 0 and r-1
		b2m = b;
		for (i=0;i<(int)r;i++)
		{
			if (b2m==1) break;

			//successivly square b mod p
			spMulMod(b2m,b2m,p,&b2m);
		}

		if (i == 0)
		{
			*sq = x;
			goto free;
		}

		//replace x by x*g^(2^(r-m-1))
		spModExp(g,1<<(r-i-1),p,&tmp);
		spMulMod(tmp,x,p,&x);

		//replace g by g^(2^(r-m)) and
		//replace b by b*g
		spModExp(g,1<<(r-i),p,&g);
		spMulMod(g,b,p,&b);
		
		r = i;
	}
		
free:
	//return the smallest solution always
	if (*sq > (p>>1)) 
		*sq = p - *sq;

	return;
}

uint32 modinv_1(uint32 a, uint32 p) {

	/* thanks to the folks at www.mersenneforum.org */

	uint32 ps1, ps2, parity, dividend, divisor, rem, q, t;

	q = 1;
	rem = a;
	dividend = p;
	divisor = a;
	ps1 = 1;
	ps2 = 0;
	parity = 0;

	while (divisor > 1) {
		rem = dividend - divisor;
		t = rem - divisor;
		if (rem >= divisor) {
			q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) {
			q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) {
			q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) {
			q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) {
			q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) {
			q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) {
			q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) {
			q += ps1; rem = t;
		if (rem >= divisor) {
			q = dividend / divisor;
			rem = dividend % divisor;
			q *= ps1;
		} } } } } } } } }

		q += ps2;
		parity = ~parity;
		dividend = divisor;
		divisor = rem;
		ps2 = ps1;
		ps1 = q;
	}
	
	if (parity == 0)
		return ps1;
	else
		return p - ps1;
}

uint32 modinv_1b(uint32 a, uint32 p) {

	/* thanks to the folks at www.mersenneforum.org */

	/* modification: p is fixed at 2^32.  a is only valid if odd */
	
	uint64 dividend = (uint64)0x1 << 32;
	uint32 ps1, ps2, parity, divisor, rem, q, t;

	q = 1;
	rem = a;
	//dividend = p;
	divisor = a;
	ps1 = 1;
	ps2 = 0;
	parity = 0;

	while (divisor > 1) {
		rem = (uint32)(dividend - (uint64)divisor);
		t = rem - divisor;
		if (rem >= divisor) {
			q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) {
			q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) {
			q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) {
			q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) {
			q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) {
			q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) {
			q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) {
			q += ps1; rem = t;
		if (rem >= divisor) {
			q = (uint32)(dividend / (uint64)divisor);
			rem = (uint32)(dividend % (uint64)divisor);
			q *= ps1;
		} } } } } } } } }

		q += ps2;
		parity = ~parity;
		dividend = divisor;
		divisor = rem;
		ps2 = ps1;
		ps1 = q;
	}
	
	if (parity == 0)
		return ps1;
	else
		return 0xFFFFFFFF - ps1 + 1;
}



int zExp(uint32 e, z *u, z *w)
{
	//return u^e = w

	mpz_t result, base;

	mpz_init(result);
	mpz_init(base);

	mp2gmp(u, base);
	mpz_pow_ui(result, base, e);


#if ((GMP_LIMB_BITS == 64) && (BITS_PER_DIGIT == 32))
    zGrow(w, 2 * result->_mp_size);
#else
    zGrow(w, result->_mp_size);
#endif
	gmp2mp(result, w);

	mpz_clear(result);
	mpz_clear(base);
	w->type = COMPOSITE;
	return 1;
}

int zPrimorial(uint32 n, z *w)
{
	//return n# = p1 * p2 * p3 ... all the primes < n
	uint32 i;
	fp_digit q;
	int approx_words = (int)(10*n/DEC_DIGIT_PER_WORD);  
	//same approximation as in factorial.  this will be overkill...

	//allocate about 10x more digits than n
	if (w->alloc < approx_words)
		zGrow(w,approx_words);

	free(PRIMES);
	PRIMES = soe_wrapper(spSOEprimes, szSOEp, 0, n, 0, &NUM_P);
	P_MIN = PRIMES[0];
	P_MAX = PRIMES[NUM_P-1];
	w->size = 1;
	w->val[0] = 1;
	//naive (but simple) method
	for (i=0; PRIMES[i] <= n; i++)
	{
		q = (fp_digit)PRIMES[i];
		zShortMul(w,q,w);
	}

	w->type = COMPOSITE;

	return 1;
}

int zFactorial(uint32 n, z *w)
{
	//return n! = n*(n-1)*(n-2)*...*(1)
	mpz_t result;

	mpz_init(result);
	mpz_fac_ui(result, n);

#if ((GMP_LIMB_BITS == 64) && (BITS_PER_DIGIT == 32))
    zGrow(w, 2 * result->_mp_size);
#else
    zGrow(w, result->_mp_size);
#endif
	gmp2mp(result, w);

	mpz_clear(result);
	w->type = COMPOSITE;
	return 1;
}

void zModExp(z *a, z *b, z *m, z *u)
{
	//computes a^b mod m = u using the binary method
	//see, for instance, the handbook of applied cryptography
	z n,bb,aa,q,r,t;

	zInit(&aa);
	zInit(&bb);
	zInit(&n);
	zInit(&q);
	zInit(&r);
	zInit(&t);

	//overflow possibilities:
	//t ranges to 2x input 'a'
	//u needs at least as much space as m
	zCopy(&zOne,&n);
	zCopy(a,&aa);
	zCopy(b,&bb);
	while (!isZero(&bb))
	{
		if (bb.val[0] & 0x1)
		{
			zMul(&n,&aa,&t);  //n*a
			zDiv(&t,m,&q,&n);   //n*a mod m
		}
		zShiftRight_1(&bb,&bb);   //compute successive squares of a
		zSqr(&aa,&t);
		zDiv(&t,m,&q,&aa);
		if (aa.size < 0)
			aa.size *= -1;
	}
	zCopy(&n,u);
	u->type = UNKNOWN;

	zFree(&aa);
	zFree(&bb);
	zFree(&n);
	zFree(&q);
	zFree(&r);
	zFree(&t);
	return;
}

void zModExp_1(z *a, fp_digit b, fp_digit m, fp_digit *u)
{
	//computes a^b mod m = u using the binary method
	//see, for instance, the handbook of applied cryptography
	z aa,t;
	fp_digit n,bb,ut;

	zInit(&aa);
	zInit(&t);

	n=1;
	zCopy(a,&aa);
	bb = b;
	while (bb != 0)
	{
		if (bb & 0x1)
		{
			zShortMul(&aa,n,&t);  //n*a
			n = zShortMod(&t,m);   //n*a mod m
		}
		bb >>= 1;   
		//compute successive squares of a
		zSqr(&aa,&t);
		ut = zShortMod(&t,m);
		sp2z(ut,&aa);
	}
	*u = n;

	zFree(&aa);
	zFree(&t);
	return;
}

void spModExp(fp_digit a, fp_digit b, fp_digit m, fp_digit *u)
{
	//computes a^b mod m = u using the binary method
	//see, for instance, the handbook of applied cryptography
	fp_digit n,bb,aa,t,prod[2];

	n=1;
	aa = a;
	bb = b;
	while (bb != 0)
	{
		if (bb & 0x1)
		{
			spMultiply(aa,n,&prod[0],&prod[1]);		//n*a
			spDivide(&t,&n,prod,m);					//n*a mod m
		}
		bb >>= 1;   
		//compute successive squares of a
		spMultiply(aa,aa,&prod[0],&prod[1]);
		spDivide(&t,&aa,prod,m);
	}
	*u = n;

	return;
}


int zNroot(z *u, z *w, int n)
{
	mpz_t result, base;

	mpz_init(result);
	mpz_init(base);

	mp2gmp(u, base);
	if (n == 2)
		mpz_sqrt(result, base);
	else
		mpz_root(result, base, n);

#if ((GMP_LIMB_BITS == 64) && (BITS_PER_DIGIT == 32))
    zGrow(w, 2 * result->_mp_size);
#else
    zGrow(w, result->_mp_size);
#endif
	gmp2mp(result, w);

	mpz_clear(result);
	mpz_clear(base);
	w->type = COMPOSITE;

	return 1;
}

void lucas(uint32 n, z *L)
{
	//compute the nth lucas number
	//Ln = Ln-1 + Ln-2 for L0 = 2, L1 = 1
	z L1, L2;
	uint32 i;

	zInit(&L1);
	zInit(&L2);

	zCopy(&zOne,&L1);
	zCopy(&zTwo,&L2);

	if (n==1)
	{
		zCopy(&L1,L);
		goto free;
	}

	if (n==0)
	{
		zCopy(&L2,L);
		goto free;
	}

	for (i=2;i<=n;i++)
	{
		zAdd(&L1,&L2,L);
		zCopy(&L1,&L2);
		zCopy(L,&L1);
	}

free:
	zFree(&L1);
	zFree(&L2);
	L->type = UNKNOWN;
	return;
}

void fib(uint32 n, z *F)
{
	//compute the nth fibonacci number
	//Fn = Fn-1 + Fn-2 for F0 = 0, L1 = 1
	z F1, F2;
	uint32 i;

	zInit(&F1);
	zInit(&F2);

	zCopy(&zOne,&F1);
	zCopy(&zZero,&F2);

	if (n==1)
	{
		zCopy(&F1,F);
		goto free;
	}

	if (n==0)
	{
		zCopy(&F2,F);
		goto free;
	}

	for (i=2;i<=n;i++)
	{
		zAdd(&F1,&F2,F);
		zCopy(&F1,&F2);
		zCopy(F,&F1);
	}

free:
	zFree(&F1);
	zFree(&F2);

	F->type = UNKNOWN;
	return;
}

