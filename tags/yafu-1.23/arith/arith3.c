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

/*
implements some more advanced arithmatic and/or number
theoretic routines
*/
int d_pull_twos(double *n, int *j, double p);
int pull_twos(fp_digit *n, int *j, fp_digit p);
int z_pull_twos(z *n, int *j, z *p);

int d_pull_twos(double *n, int *j, double p)
{
	int c = 0;

	while (fmod(*n,2) == 0)
	{
		*n /= 2;
		c = 1 - c;
	}
	if ((c * fmod(p*p-1,16)) == 8)
		*j *= -1;
	return c;
}

int pull_twos(fp_digit *n, int *j, fp_digit p)
{
	int c = 0;

	while (!(*n & 1))
	{
		*n >>= 1;
		c = 1 - c;
	}
	if ((c * (p*p-1)%16) == 8)
		*j *= -1;
	return c;
}

int z_pull_twos(z *n, int *j, z *p)
{
	//n is overwritten
	int c = 0;
	z t1, t2;
	fp_digit r;

	zInit(&t1);
	zInit(&t2);

	while (!(n->val[0] & 1))
	{
		zShiftRight(n,n,1);
		c = 1 - c;
	}
	zExp(2,p,&t2);
	zSub(&t2,&zOne,&t1);
	r = zShortDiv(&t1,16,&t2);

	if ((c * r) == 8)
		*j *= -1;

	zFree(&t1);
	zFree(&t2);
	return c;
}


int isPrime(z *n)
{
	/*first check by trial division for the small primes.  if none divide n, 
	then all small primes checked in this way are relatively prime to n.
	then use fermat's observation that b^(n-1) mod n = 1 if n is prime and b and n 
	are relatively prime.  check for several small primes: b1, b2... bn.  if all are 1, 
	then n is probably a prime, or exceedingly rarely, is a pseudoprime for 
	bases b1, b2, ... bn

    see knuth TAOCP Vol 2.
	*/

	int i,j;
	fp_digit a;
	z nn_1, t, d, q;
	
	if (zCompare(n,&zOne) <= 0)
		return 0;

	//trust this info
	if (n->type == PRIME)
		return 1;

	if (n->type == PRP)
		return 1;

	if (n->type == COMPOSITE)
		return 0;

	//trial divide by the first 100 primes
	for(i=0;i<100 || i < (int)NUM_WITNESSES ;i++)
    {
		if (n->size == 1 && n->val[0] == spSOEprimes[i])
		{
			n->type = PRIME;
			return 1;
		}

		if (zShortMod(n,(fp_digit)spSOEprimes[i]) == 0)
		{
			n->type = COMPOSITE;
			return 0;
		}
    }

	zInit(&d);
	zInit(&t);
	zInit(&q);
	zInit(&nn_1);

	if (t.alloc < n->size)
	{
		zGrow(&t,n->size);
		zGrow(&nn_1,n->size);
		zGrow(&q,n->size);
	}
		
	//knuth's algorithm P: simplified version of the rabin-miller strong
	//pseudoprime test

	zSub(n,&zOne,&nn_1);
	zCopy(&nn_1,&t);
	a=0;

	//find t and a satisfying: n-1 = 2^a * t, t odd
	while (!(t.val[0] & 0x1))
	{
		zShiftRight(&t,&t,1);
		a++;
	}

	for(j=0;j<(int)NUM_WITNESSES;j++)
	{
		//use the first N primes as witnesses
		d.size = 1;
		d.val[0] = (fp_digit)spSOEprimes[j];

 		zModExp(&d,&t,n,&q);
		d.size = 1;
		d.val[0] = 2;
		i=0;
		while (1)
		{
			if (i > 0 && q.size == 1 && q.val[0] == 1)
			{
				zFree(&d);
				zFree(&t);
				zFree(&q);
				zFree(&nn_1);
				n->type = COMPOSITE;
				return 0;
			}
			if (i==0 && q.size == 1 && q.val[0] == 1)
				break;
			if (zCompare(&q,&nn_1) == 0)
				break;
			i++;
			if ((fp_digit)i >= a)
			{
				zFree(&d);
				zFree(&t);
				zFree(&q);
				zFree(&nn_1);
				n->type = COMPOSITE;
				return 0;
			}
			zModExp(&q,&d,n,&q);
		}
	}

	zFree(&d);
	zFree(&t);
	zFree(&q);
	zFree(&nn_1);
	n->type = PRP;
	return 1;
}

int isSquare(z *n)
{
	//thanks fenderbender @ mersenneforum.org
	unsigned long m;
	unsigned long largeMod;
	z w2,w3;
	int ans;

	// start with mod 128 rejection. 82% rejection rate
	// VERY fast, can read bits directly
	m=n->val[0] & 127; // n mod 128
	if ((m*0x8bc40d7d) & (m*0xa1e2f5d1) & 0x14020a) return 0; 

	//Other modulii share one BigInt modulus.
	largeMod=zShortMod(n,(63UL*25*11*17*19*23*31)); // SLOW, bigint modulus

	// residues mod 63. 75% rejection
	m=largeMod%63; // fast, all 32-bit math
	if ((m*0x3d491df7) & (m*0xc824a9f9) & 0x10f14008) return 0;

	// residues mod 25. 56% rejection
	m=largeMod%25; 
	if ((m*0x1929fc1b) & (m*0x4c9ea3b2) & 0x51001005) return 0;

	// residues mod 31. 48.4% rejection
	//  Bloom filter has a little different form to keep it perfect
	m=0xd10d829a*(largeMod%31); 
	if (m & (m+0x672a5354) & 0x21025115) return 0;

	// residues mod 23. 47.8% rejection
	m=largeMod%23; 
	if ((m*0x7bd28629) & (m*0xe7180889) & 0xf8300) return 0;

	// residues mod 19. 47.3% rejection
	m=largeMod%19; 
	if ((m*0x1b8bead3) & (m*0x4d75a124) & 0x4280082b) return 0;

	// residues mod 17. 47.1% rejection
	m=largeMod%17; 
	if ((m*0x6736f323) & (m*0x9b1d499) & 0xc0000300) return 0;

	// residues mod 11. 45.5% rejection
	m=largeMod%11; 
	if ((m*0xabf1a3a7) & (m*0x2612bf93) & 0x45854000) return 0;

	// Net nonsquare rejection rate: 99.92%

	// We COULD extend to another round, doing another BigInt modulus and
	// then followup rejections here, using
	// primes of  13 29 37 41 43 53.  That'd give 98% further rejection.
	// Empirical timing shows this second round would be useful for n>10^100 or so.

	// VERY expensive final definitive test

	zInit(&w2);
	zInit(&w3);
	zNroot(n,&w2,2);	//w2 = sqrt(w1)
	zSqr(&w2,&w3);		//w3 = w2^2
	ans = zCompare(n,&w3);
	zFree(&w2);
	zFree(&w3);
	return (ans == 0);
}

int rec_jacobi_1(uint32 n, uint32 p)
{
	//compute the jacobi symbol (n/p) recursively
	//p must be odd
	//based on routine in Bressoud's book

	//return an error condition if p is even
	if (!(p & 1))
		return -2;

	if (p==1) return 1;
	else if (n==0) return 0;
	else if (n & 1)
	{
		if (((n % 4) == 3) && ((p % 4) == 3))
			return -1*rec_jacobi_1(p % n,n);
		else
			return rec_jacobi_1(p % n,n);
	}
	else if (((p % 8) == 3) || ((p % 8) == 5))
		return -1*rec_jacobi_1(n>>1,p);
	else
		return rec_jacobi_1(n>>1,p);
}

int jacobi_1(fp_digit n, fp_digit p)
{
	//compute the jacobi symbol (n/p) for positive inputs
	//p must be odd
	//based on routine in Bressoud's book

	int j = 1;
	fp_digit t,nn = n;

	//return an error condition if p is even
	if (!(p & 1))
		return -2;

	nn = nn%p;	

	//if p divides n then (n/p) = 0
	if (nn==0)
		return 0;

	/*
	if (nn<0)
	{
		nn *= -1;
		if (p%4 == 3)
			j = -1;
	}
	*/

	pull_twos(&nn,&j,p);	
	while (nn>1)
	{
		if (((nn-1)*(p-1))%8 == 4)
			j = -1*j;
		t = nn;
		nn = p%nn;
		p = t;
		
		pull_twos(&nn,&j,p);
	}
	return j;
}

int d_jacobi(double n, double p)
{
	//compute the jacobi symbol (n/p)
	//p must be odd
	//based on routine in Bressoud's book

	int j = 1;
	double t,nn = n;

	//return an error condition if p is even
	if (fmod(p,2)==0)
		return -2;

	nn = fmod(nn,p);	

	//if p divides n then (n/p) = 0
	if (nn==0)
		return 0;

	if (nn<0)
	{
		nn *= -1;
		if (fmod(p,4) == 3)
			j = -1;
	}

	//printf("pulling two's in jacobi, outer loop\n");
	d_pull_twos(&nn,&j,p);	
	while (nn>1)
	{
		if (fmod((nn-1)*(p-1),8) == 4)
			j = -1*j;
		t = nn;
		nn = fmod(p,nn);
		p = t;
		//printf("pulling two's in jacobi, inner loop\n");
		d_pull_twos(&nn,&j,p);
	}

	return j;
}

int zJacobi(z *n, z *p)
{
	//compute the jacobi symbol (n/p) using the iterative method
	//p must be odd
	//based on routine in Bressoud's book
	
	int j = 1, sz;
	fp_digit rem;
	z nn, pp, t, t3, t4;

	//return an error condition if p is even
	if (!(p->val[0] & 1))
		return -2;

	zInit(&nn);
	zInit(&pp);
	zInit(&t);
	zInit(&t3);
	zInit(&t4);

	sz = abs(n->size);
	if (abs(p->size) > sz)
		sz = abs(p->size);

	if (t.alloc < sz)
	{
		zGrow(&t,sz);
		zGrow(&t3,sz);
		zGrow(&t4,sz);
	}

	zCopy(n,&nn);
	zCopy(p,&pp);

	zDiv(&nn,&pp,&t,&t3);
	//restore nn
	zCopy(n,&nn);

	//if p divides n then (n/p) = 0
	if (zCompare(&t3,&zZero) == 0)
	{
		zFree(&nn);
		zFree(&pp);
		zFree(&t);
		zFree(&t3);
		zFree(&t4);
		return 0;
	}

	z_pull_twos(&nn,&j,&pp);	
	while (zCompare(&nn,&zOne) > 0)	//nn > 1
	{
		zShortSub(&pp,1,&t);
		zShortSub(&nn,1,&t4);
		zMul(&t4,&t,&t3);
		rem = zShortMod(&t3,8);
		if (rem == 4)
			j = -1*j;
		zCopy(&nn,&t);
		zDiv(&pp,&nn,&t3,&t4);
		zCopy(&t4,&nn);
		zCopy(&t,&pp);
		z_pull_twos(&nn,&j,&pp);	
	}

	zFree(&nn);
	zFree(&pp);
	zFree(&t);
	zFree(&t3);
	zFree(&t4);
	return j;
}

int zGCD(z *u, z *v, z *w)
{
	//computes w = gcd(u,v) 
	//follows algorithm A.  p.320 Knuth Vol. 2
	//needs temp storage, so that u is not modified
	z uu, vv, q;

	int i=0, retcode,sz;

	zInit(&uu);
	zInit(&vv);
	zInit(&q);

	sz = abs(u->size);
	if (abs(v->size) > sz)
		sz = abs(v->size);

	if (w->alloc < sz)
		zGrow(w,sz);

	i = zCompare(u,v);
	if (i >= 0)
	{
		zCopy(u,&uu);
		zCopy(v,&vv);
	}
	else
	{
		zCopy(v,&uu);
		zCopy(u,&vv);
	}

	for (i=1;i<=10000;++i)		//safety net
	{
		//if v == 0, terminate with u as the answer
		if (zCompare(&vv, &zZero) == 0)
		{
			zCopy(&uu,w);
			retcode = i;
			goto free;
		}

		//set w = u mod v, u = v, v = w and loop
		//uu gets destroyed by the divide.  doesn't matter because we overwrite it 
		//during the step after anyway.
		//q is never used.  seems a shame to have a temp variable for it to just be ignored.
		zDiv(&uu,&vv,&q,w);

		//set uu = vv
		zCopy(&vv,&uu);

		//set vv = w
		zCopy(w,&vv);
	}
	if (i >= 10000) 
		retcode = 0;

free:
	zFree(&uu);
	zFree(&vv);
	zFree(&q);
	return retcode;
}

void xGCD_1(int a, int b, int *x, int *y, int *g)
{
	//compute the extended GCD of a, b, returning g = GCD(a,b) and x, y 
	//such that ax + by = GCD(a,b) if a,b are coprime
	int t1,t2,t3,u,v,r,R,q;

//	int i;
	/*

	Step 1: 
	if a < b then 
	Set u=0, v=1, and r=b 
	Set U=1, V=0, and R=a 
	else 
	Set u=1, v=0, and r=a 
	Set U=0, V=1, and R=b 

	Step 2: 
	if R = 0 then return r (for the gcd) and no inverses exist. 
	if R = 1 then return R (for the gcd), V (for the inverse a(mod b)) and U (for the inverse of b(mod a)). 
	
	Step 3: 
	Calculate q = int(r/R) 
	Calculate t1 = u - U*q 
	Calculate t2 = v - V*q 
	Calculate t3 = r - R*q 
	set u=U, v=V, r=R 
	set U=t1, V=t2, R=t3 
	goto Step 2. 
	*/

	*x = *y = 0;

	if (a < b)
	{
		u=0;
		v=1;
		r=b;
		*x=1;
		*y=0;
		R=a;
	}
	else
	{
		u=1;
		v=0;
		r=a;
		*x=0;
		*y=1;
		R=b;
	}

	while (1)
	{
		if (R == 0)
		{
			*g = r;
			*x = 0;
			*y = 0;
			break;
		}

		if (R == 1)
		{
			*g = R;
			break;
		}

		//Calculate q = int(r/R) 
		q = r / R;
		t3 = r % R;

		//Calculate t1 = u - U*q 
		t1 = u - q * *x;

		//Calculate t2 = v - V*q 
		t2 = v - *y * q;

		//set u=U, v=V, r=R 
		u = *x;
		v = *y;
		r = R;

		//set U=t1, V=t2, R=t3 
		*x = t1;
		*y = t2;
		R = t3;
	}

	return;
}

void xGCD(z *a, z *b, z *x, z *y, z *g)
{
	//compute the extended GCD of a, b, returning g = GCD(a,b) and x, y 
	//such that ax + by = GCD(a,b) if a,b are coprime
	z t1,t2,t3,u,v,r,R,q,tmp;

//	int i;
	/*

	Step 1: 
	if a < b then 
	Set u=0, v=1, and r=b 
	Set U=1, V=0, and R=a 
	else 
	Set u=1, v=0, and r=a 
	Set U=0, V=1, and R=b 

	Step 2: 
	if R = 0 then return r (for the gcd) and no inverses exist. 
	if R = 1 then return R (for the gcd), V (for the inverse a(mod b)) and U (for the inverse of b(mod a)). 
	
	Step 3: 
	Calculate q = int(r/R) 
	Calculate t1 = u - U*q 
	Calculate t2 = v - V*q 
	Calculate t3 = r - R*q 
	set u=U, v=V, r=R 
	set U=t1, V=t2, R=t3 
	goto Step 2. 
	*/


	zInit(&tmp);
	zInit(&t1);
	zInit(&t2);
	zInit(&t3);
	zInit(&q);
	zInit(&r);
	zInit(&R);
	zInit(&u);
	zInit(&v);

	//need to check for temp allocation

	zClear(x);
	zClear(y);


	if (zCompare(a,b) < 0)
	{
		u.val[0]=0;
		v.val[0]=1;
		zCopy(b,&r);
		x->val[0]=1;
		y->val[0]=0;
		zCopy(a,&R);
	}
	else
	{
		u.val[0]=1;
		v.val[0]=0;
		zCopy(a,&r);
		x->val[0]=0;
		y->val[0]=1;
		zCopy(b,&R);
	}

	while (1)
	{
		if (zCompare(&zZero,&R) == 0)
		{
			zCopy(&r,g);
			zCopy(&zZero,x);
			zCopy(&zZero,y);
			break;
		}

		if (zCompare(&zOne,&R) == 0)
		{
			zCopy(&R,g);
			break;
		}

		zCopy(&r,&tmp);
		zDiv(&tmp,&R,&q,&t3);		//q = int(r/R), t3 = r % R

		zMul(&q,x,&tmp);			//t1 = u - U*q
		zSub(&u,&tmp,&t1);

		zMul(&q,y,&tmp);			//t2 = v - V*q
		zSub(&v,&tmp,&t2);

		zCopy(x,&u);
		zCopy(y,&v);
		zCopy(&R,&r);

		zCopy(&t1,x);
		zCopy(&t2,y);
		zCopy(&t3,&R);
	}

	if (x->size < 0)
	{
		x->size *= -1;
		zSub(b,x,x);
	}

	if (y->size < 0)
	{
		y->size *= -1;
		zSub(a,y,y);
	}

	zFree(&tmp);
	zFree(&t1);
	zFree(&t2);
	zFree(&t3);
	zFree(&q);
	zFree(&r);
	zFree(&R);
	zFree(&u);
	zFree(&v);
	x->type = UNKNOWN;
	y->type = UNKNOWN;
	g->type = UNKNOWN;
	return;
}

int zBinGCD(z *u, z *v, z *w)
{
	//computes w = gcd(u,v) 
	//follows algorithm B.  p.321 Knuth Vol. 2

	z uu, vv, t;
	long i=0,j;
	int k,sz;

	sz = abs(u->size);
	if (abs(v->size) > sz)
		sz = abs(v->size);

	if (w->alloc < sz)
		zGrow(w,sz);

	w->type = UNKNOWN;

	i = zCompare(u,&zZero);
	j = zCompare(v,&zZero);
	if (i == 0) 
	{
		zCopy(v,w);
		return 1;
	}
	if (j == 0)
	{
		zCopy(u,w);
		return 1;
	}
		
	zInit(&uu);
	zInit(&vv);
	zInit(&t);

	if (t.alloc < sz)
		zGrow(&t,sz);

	zCopy(u,&uu);
	zCopy(v,&vv);

	//find power of 2 such that u and v are not both even
	k = 0;
	while(((uu.val[0] & 0x1) == 0) && ((vv.val[0] & 0x1) == 0))
	{
		zShiftRight(&uu,&uu,1);
		zShiftRight(&vv,&vv,1);
		k++;
	}

	j=0;
	do
	{
		if ((uu.val[0] & 0x1) == 0)
			zShiftRight(&uu,&uu,1);
		else if ((vv.val[0] & 0x1) == 0)
			zShiftRight(&vv,&vv,1);
		else
		{
			zSub(&uu,&vv,&t);
			t.size = abs(t.size);
			zShiftRight(&t,&t,1);
			if (zCompare(&uu,&vv) < 0)
				zCopy(&t,&vv);
			else
				zCopy(&t,&uu);
		}
		++j;
		if (j>= 10000)
			break;
	} while (zCompare(&uu,&zZero) > 0);

	zClear(w);
	w->val[0] = 1;
	zShiftLeft(w,w,k);
	zMul(w,&vv,&uu);
	zCopy(&uu,w);
	
	zFree(&uu);
	zFree(&vv);
	zFree(&t);
	return j;
}

int zLEGCD(z *u, z *v, z *w)
{
	//use the Lehman-Euclid algorithm to calculate GCD(u,v) = w
	//Algorithm L in Knuth, 4.5.2 p. 329
	//assumes u,v nonnegative

	fp_digit aa,bb,cc,dd;
	int i,j,k,it;
	fp_signdigit a,b,c,d,t;
	fp_digit up,vdp, q1, q2;
	fp_digit mask;
	z y,zz;	//t and w, in knuth
	z x;	//tmp variable
	z uu, vv;	//so u and v don't get destroyed
	z uh, vh;
	fp_digit udp[2],vp[2];

#if BITS_PER_DIGIT == 32
		mask = 0xff000000;
#else
		mask = 0xff00000000000000;
#endif

	w->type = UNKNOWN;

	i = zCompare(u,&zZero);
	j = zCompare(v,&zZero);
	if (i == 0) 
	{
		zCopy(v,w);
		return 1;
	}
	if (j == 0)
	{
		zCopy(u,w);
		return 1;
	}

	//temp variables should be twice as big as the input, to make room
	//for intermediate operations.  w should be as big as the biggest input.
	i = u->size;
	j = v->size;
	if (j > i)
		i = j;
	if (w->alloc < i)
		zGrow(w,i + LIMB_BLKSZ);

	zInit(&y);
	zInit(&zz);
	zInit(&x);
	zInit(&uu);
	zInit(&vv);
	zInit(&uh);
	zInit(&vh);

	if (y.alloc < i)
	{
		zGrow(&y,i + LIMB_BLKSZ);
		zGrow(&zz,i + LIMB_BLKSZ);
		zGrow(&x,i + LIMB_BLKSZ);
		zGrow(&uh,i + LIMB_BLKSZ);
		zGrow(&uu,i + LIMB_BLKSZ);
		zGrow(&vv,i + LIMB_BLKSZ);
		zGrow(&vh,i + LIMB_BLKSZ);
	}

	//put bigger number in uu, other in vv.
	i = zCompare(u,v);
	if (i >= 0)
	{
		zCopy(u,&uu);
		zCopy(v,&vv);
	}
	else
	{
		zCopy(v,&uu);
		zCopy(u,&vv);
	}

	j=0;
	while (vv.size > 1)
	{
		//Step L1
		for (it=vv.size;it<uu.size;it++)
			vv.val[it]=0;
		vv.size = uu.size;
		//get the most significant 32 bits of u and v, such that uhat >= vhat
		uh.val[1] = uu.val[uu.size - 1];
		uh.val[0] = uu.val[uu.size - 2];
		vh.val[1] = vv.val[vv.size - 1];
		vh.val[0] = vv.val[vv.size - 2];
		uh.size = vh.size = 2;

		//rightshift until uhat is a single word
		//0xff000000 is magic
		if ((uh.val[1] & mask) > 0)
		{
			uh.val[0] = uh.val[1];
			vh.val[0] = vh.val[1];
		}
		else
		{
			i=0;
			aa=uh.val[1];
			while ((aa & MAX_DIGIT) != 0)
			{
				aa >>= 1;
				i++;
			}
			zShiftRight(&uh,&uh,i);
			zShiftRight(&vh,&vh,i);
		}
		
		//make u',v',u'',v''
		up = uh.val[0];
		vdp = vh.val[0];
		if (up == MAX_DIGIT)
		{
			udp[0] = 0;
			udp[1] = 1;
		}
		else
		{
			udp[0] = up+1;
			udp[1] = 0;
		}

		if (vdp == MAX_DIGIT)
		{
			vp[0] = 0;
			vp[1] = 1;
		}
		else 
		{
			vp[0] = vdp+1;
			vp[1] = 0;
		}

		a=1; b=0; c=0; d=1; 
		
		k=0;
		while (1)
		{
			//Step L2:
			/*test quotient, protecting for overflow.  the conditions:
			0 <= uhat + a <= 2^32
			0 <= uhat + b < 2^32
			0 <= vhat + c < 2^32
			0 <= vhat + d <= 2^32
			will always hold.  hence only need to check for the case where
			uhat == MAX_DIGIT and a = 1 or vhat == MAX_DIGIT and d = 1
			*/

			//first check for /0
			if (((vp[0] == 0) && (vp[1] == 0)) || (vdp == 0))
				break;

			//u''/v''
			if (udp[1] == 1)
				spDivide(&q2,&t,udp,vdp);
			else
				q2 = udp[0]/vdp;

			//u'/v'
			if (vp[1] >= 1)
				q1=0;
			else
				q1 = up/vp[0];

			if (q1 != q2)
				break;

			//Step L3: Emulate Euclid
			t=a-q1*c;
			a=c;
			c=t;
			t=b-q1*d;
			b=d;
			d=t;
			t=up-q1*vp[0];
			up=vp[0]; 
			vp[0]=t;
			t=udp[0]-q2*vdp;
			udp[0]=vdp;
			vdp=t;
			k++;
			if (k>10000)
				goto free;
		}

		//Step L4: multiprecision step
		if (b==0)
		{
			for (i=vv.size-1;i>=0;i--)
			{
				if (vv.val[i] == 0)
					vv.size--;
				else
					break;
			}
			zDiv(&uu,&vv,&y,&x);	//y = u mod v
			zCopy(&vv,&uu);			//u = v
			zCopy(&x,&vv);			//v = y
		}
		else
		{
			//aa=abs(a);
			//bb=abs(b);
			//cc=abs(c);
			//dd=abs(d);
			if (a<0)
				aa = -a;
			else
				aa = a;
			if (b<0)
				bb = -b;
			else
				bb = b;
			if (c<0)
				cc = -c;
			else
				cc = c;
			if (d<0)
				dd = -d;
			else
				dd = d;
			zShortMul(&uu,aa,&y);	//y = A*u
			zShortMul(&vv,bb,&x);	//y = y + B*v
			if (a<0)
			{
				zSub(&x,&y,&x);
				zCopy(&x,&y);
			}
			else if (b<0)
				zSub(&y,&x,&y);
			else
				zAdd(&y,&x,&y);				
			
			zShortMul(&uu,cc,&zz);	//z = c*u
			zShortMul(&vv,dd,&x);	//z = z + d*v
			if (c<0)
			{
				zSub(&x,&zz,&x);
				zCopy(&x,&zz);
			}
			else // (d<0)
				zSub(&zz,&x,&zz);

			zCopy(&y,&uu);				//u = y;
			zCopy(&zz,&vv);				//v = z;
		}
		j++;
		if (j>10000)
			goto free;
	}

	//here, the size of v is 1, so finish up with regular GCD
	zBinGCD(&uu,&vv,w);

free:
	zFree(&y);
	zFree(&zz);
	zFree(&x);
	zFree(&uu);
	zFree(&vv);
	zFree(&uh);
	zFree(&vh);
	return 1;
}

fp_digit spGCD(fp_digit x, fp_digit y)
{
	fp_digit a,b,c;
	a=x; b=y;
	while (b != 0)
	{
		c=a%b;
		a=b;
		b=c;
	}
	return a;
}

fp_digit spBinGCD(fp_digit x, fp_digit y)
{
	int k;
	//fp_digit a,b,c;
	fp_signdigit t;
	k = 0;
	while (!((x & 0x1) || (y & 0x1)))
	{
		k++;
		x >>= 1;
		y >>= 1;
	}
	if (x & 0x1)
		t = -1 * y;
	else
		t = x;

	//fix this!
	return 0;

}

uint64 gcd64(uint64 x, uint64 y)
{
	uint64 a,b,c;
	a=x; b=y;
	while (b != 0)
	{
		c=a%b;
		a=b;
		b=c;
	}
	return a;
}

void dblGCD(double x, double y, double *w)
{
	double a,b,c;
	a=x; b=y;
	while (b != 0)
	{
		c=a-b*(floor(a/b));
		a=b;
		b=c;
	}
	*w=a;
	return;
}

int llt(uint32 exp)
{
	z tmp,tmp2,n;
	clock_t start, stop;
	double t;
	uint32 i,j,nchars;
	fp_digit d;

	zInit(&tmp);
	sp2z(exp,&tmp);
	if (!isPrime(&tmp))
	{
		zFree(&tmp);
		return 0;
	}

	start = clock();
	zInit(&n);
	sp2z(4,&tmp);
	zShiftLeft(&n,&zOne,exp);
	zShortSub(&n,1,&n);
	i = (uint32)n.val[0];
	//should vary the depth depending on the size of p
	for (i = 1; i< 100000; i++)
	{
		d = 2*i*(fp_digit)exp + 1;
		if (zShortMod(&n,d) == 0)
		{
			zFree(&n);
			zFree(&tmp);
			printf("2*%d*p+1 is a factor\n",i);
			return 0;
		}
	}
	printf("trial division to %u bits is complete\n",
		(uint32)spBits(2*100000*exp+1));

	/*
	t1 = POLLARD_STG1_MAX;
	t2 = POLLARD_STG2_MAX;
	POLLARD_STG1_MAX=1000;
	POLLARD_STG2_MAX=50000;
	i = (uint32)n.val[0];
	zInit(&tmp2);
	mpollard(&n,3,&tmp2);
	POLLARD_STG1_MAX=t1;
	POLLARD_STG2_MAX=t2;
	if (n.val[0] != i)
	{
		printf("pm1 found a factor %s\n",z2decstr(&tmp2,&gstr1));
		zFree(&n);
		zFree(&tmp);
		zFree(&tmp2);
		return 0;
	}
	*/
	
	zInit(&tmp2);
	nchars=0;
	//else do the ll test
	for (i=0;i<exp-2;i++)
	{
		zSqr(&tmp,&tmp);
		zShortSub(&tmp,2,&tmp);
		zDiv(&tmp,&n,&tmp2,&tmp);
		for (j=0;j<nchars;j++)
			printf("\b");
		nchars = printf("llt iteration %d",i);
		fflush(stdout);
	}
	printf("\n");
	
	stop = clock();
	t = (double)(stop - start)/(double)CLOCKS_PER_SEC;
	printf("elapsed time = %6.4f\n",t);

	zFree(&tmp2);
	zFree(&n);

	if (zCompare(&tmp,&zZero) == 0)
	{
		zFree(&tmp);
		return 1;
	}
	else
	{
		zFree(&tmp);
		return 0;
	}
}
