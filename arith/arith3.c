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
#include "gmp.h"
#include "gmp_xface.h"
#include "mpz_aprcl.h"

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

int is_mpz_prp(mpz_t n)
{
	return mpz_probab_prime_p(n,NUM_WITNESSES) && 
		mpz_strongbpsw_prp(n);
}

int isPrime(z *n)
{
	// translate into gmp's test
	mpz_t gmpz;
	int i;

	mpz_init(gmpz);
		
	mp2gmp(n, gmpz);
	i = is_mpz_prp(gmpz);

	mpz_clear(gmpz);
	return (i > 0);
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

// straight from wikipedia.
uint64 spBinGCD(uint64 u, uint64 v)
{
    // binary GCD for non-zero inputs.
    int shift;

    /* Let shift := lg K, where K is the greatest power of 2
    dividing both u and v. */
    for (shift = 0; ((u | v) & 1) == 0; ++shift) {
        u >>= 1;
        v >>= 1;
    }

    while ((u & 1) == 0)
        u >>= 1;

    /* From here on, u is always odd. */
    do {
        /* remove all factors of 2 in v -- they are not common */
        /*   note: v is not zero, so while will terminate */
        while ((v & 1) == 0)  /* Loop X */
            v >>= 1;

        /* Now u and v are both odd. Swap if necessary so u <= v,
        then set v = v - u (which is even). For bignums, the
        swapping is just pointer movement, and the subtraction
        can be done in-place. */
        if (u > v) {
            uint64 t = v; v = u; u = t;
        }  // Swap u and v.
        v = v - u;                       // Here v >= u.
    } while (v != 0);

    /* restore common factors of 2 */
    return u << shift;
}

// assume u is odd
uint64 spBinGCD_odd(uint64 u, uint64 v)
{
    /* From here on, u is always odd. */
    do {
        /* remove all factors of 2 in v -- they are not common */
        /*   note: v is not zero, so while will terminate */
        while ((v & 1) == 0)  /* Loop X */
            v >>= 1;

        /* Now u and v are both odd. Swap if necessary so u <= v,
        then set v = v - u (which is even). For bignums, the
        swapping is just pointer movement, and the subtraction
        can be done in-place. */
        if (u > v) {
            uint64 t = v; v = u; u = t;
        }  // Swap u and v.
        v = v - u;                       // Here v >= u.
    } while (v != 0);

    /* restore common factors of 2 */
    return u;
}

// much faster version: assuming x is odd
uint64 bingcd64(uint64 x, uint64 y)
{
    if (y) {
        y >>= _trail_zcnt64(y);
        while (x != y)
            if (x < y)
                y -= x, y >>= _trail_zcnt64(y);
            else
                x -= y, x >>= _trail_zcnt64(x);
    }
    return x;
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
	mpz_t tmp,tmp2,n;
	clock_t start, stop;
	double t;
	uint32 i,j,nchars;
	fp_digit d;

	mpz_init(tmp);
	mpz_set_ui(tmp, exp);
	if (!is_mpz_prp(tmp))
	{
		mpz_clear(tmp);
		printf("exponent is not prime\n");
		return 0;
	}

	start = clock();
	mpz_init(n);
    mpz_setbit(n, exp);
    mpz_sub_ui(n, n, 1);
	//should vary the depth depending on the size of p
	for (i = 1; i < MIN(sqrt(exp)-1, 1000000); i++)
	{
		d = 2*i*(fp_digit)exp + 1;
		if (mpz_tdiv_ui(n,d) == 0)
		{
			mpz_clear(n);
			mpz_clear(tmp);
            if (VFLAG > 1)
			    printf("2*%d*p+1 is a factor\n",i);
			return 0;
		}
	}
    if (VFLAG > 1)
	    printf("trial division to %u bits is complete\n",
		(uint32)spBits(2*1000000*exp+1));

    mpz_init(tmp2);
    mpz_set_ui(tmp, 4);
	nchars=0;
	//else do the ll test
	for (i=0;i<exp-2;i++)
	{
		mpz_mul(tmp, tmp, tmp); 
		mpz_sub_ui(tmp, tmp, 2);
        /* Adapted from http://rosettacode.org/wiki/Lucas-Lehmer_test#GMP */
        /* mpz_tdiv_r(tmp, tmp, n); but more efficiently done given mod 2^p-1 */
        if (mpz_sgn(tmp) < 0) mpz_add(tmp, tmp, n);
        /* while (n > mp) { n = (n >> p) + (n & mp) } if (n==mp) n=0 */
        /* but in this case we can have at most one loop plus a carry */
        mpz_tdiv_r_2exp(tmp2, tmp, exp);
        mpz_tdiv_q_2exp(tmp, tmp, exp);
        mpz_add(tmp, tmp, tmp2);
        while (mpz_cmp(tmp, n) >= 0) mpz_sub(tmp, tmp, n);
		
        if (VFLAG > 1)
        {
            if ((i & 511) == 0)
            {
                for (j = 0; j < nchars; j++)
                    printf("\b");
                nchars = printf("llt iteration %d", i);
                fflush(stdout);
            }
        }
	}
    if (VFLAG > 1)
	    printf("\n");
	
    if (VFLAG > 1)
    {
        stop = clock();
        t = (double)(stop - start) / (double)CLOCKS_PER_SEC;
        printf("elapsed time = %6.4f\n", t);
    }

	mpz_clear(n);
    mpz_clear(tmp2);

	if (mpz_cmp_ui(tmp, 0) == 0)
	{
		mpz_clear(tmp);
		return 1;
	}
	else
	{
		mpz_clear(tmp);
		return 0;
	}
}
