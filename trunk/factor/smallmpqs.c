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
       				   --bbuhrow@gmail.com 3/26/10
----------------------------------------------------------------------*/

#include "yafu.h"
#include "arith.h"
#include "factor.h"
#include "util.h"

/*
implements the multiple polynomial quadratic sieve.
inspired by Jason Papadopolous's msieve,
Scott Contini's siqs code, and tidbits from many many other papers
*/

/*
for use in splitting triple large prime composites in siqs.  these will all be
less than say 100 bits.  in the interest of speed, we can fix the factor base size
and the blocksize, and several other parameters.
*/

/************************* MPQS types and functions *****************/


#define MAX_FACTORS_SM 20
#define MAX_SMOOTH_PRIMES 50
#define NUM_EXTRA_RELS 32
uint32 SM_BLOCKSIZE;
uint32 SM_NUMBLOCKS;



typedef struct
{
	uint16 largeprime;		//large prime in the pd.
	uint16 offset;			//offset specifying Q (the quadratic polynomial)
	uint16 *fboffset;		//offsets of factor base primes dividing Q(offset).  max # of fb primes < 2^16 with this choice
	uint8 num_factors;		//number of factor base factors in the factorization of Q
	uint8 polynum;			//which poly this relation uses
	uint8 parity;			//the sign of the offset (x) 0 is positive, 1 is negative
} qs_r;

typedef struct
{
	uint32 num_r;
	uint32 act_r;
	uint32 allocated;
	qs_r **list;
} qs_rlist;

//holds all the info for a factor base element
typedef struct
{
	uint16 prime;
	uint16 proot1;
	uint16 proot2;
	uint16 nroot1;
	uint16 nroot2;
	uint16 c1;
	uint16 c2;
	uint8 logprime;
} fb_element_qs;

typedef struct
{
	z poly_a;
	z poly_b;
	z poly_c;
	uint32 poly_d;		//limits the upper end of usefullness, but ok for n < 128 bits or so
} qs_poly;

typedef struct
{
	uint32 B;
	uint32 small_B;
	uint32 med_B;
	fb_element_qs *list;
} fb_list_qs;

typedef struct
{
	uint16 prime;
	uint16 root1;
	uint16 root2;
	uint8 logprime;
} sieve_fb;

static void sm_sieve_block(uint8 *sieve, sieve_fb *fb, 
						uint8 s_init, fb_list_qs *fullfb);
static int sm_check_relations(uint16 sieve_interval, uint8 blocknum, uint8 sieve[],
							  z *n, 
								qs_poly *poly, uint8 closnuf,
								sieve_fb *fb, fb_list_qs *fullfb,
								qs_rlist *full, qs_rlist *partial,
								uint16 cutoff, 
								uint8 parity, 
								uint16 *num, int numpoly);
static void sm_trial_divide_Q(z *Q,sieve_fb *fb, qs_rlist *full, qs_rlist *partial,
						  uint8 *sieve, uint16 offset, uint16 j, uint8 sign, 
						  fb_list_qs *fullfb, uint16 cutoff,
						  int numpoly, uint8 parity,
						  uint8 closnuf);
static void sm_save_relation_qs(qs_rlist *list, uint16 offset, uint16 largeprime, 
							   uint8 num_factors, uint16 rnum, uint16 *fboffset, 
							   uint8 numpoly, uint8 parity);
void sm_make_fb_qs(fb_list_qs *fb, z *n);
int sm_qcomp_qs(const void *x, const void *y);
void sm_computeRoots(qs_poly *poly, fb_list_qs *fb, uint8 mul);
void sm_computeB(qs_poly *poly, z *n, z *work1, z *work2);
uint32 sm_nextD(uint32 poly_d, z *n, z *work);
int BlockGauss_smqs(qs_rlist *full, qs_rlist *partial, 
					z *apoly, z *bpoly,
					fb_list_qs *fb, z *n, int mul, 
					z *factors, uint32 *num_factors);
int checkpoly_qs(qs_poly *poly, z *n);

void sm_make_fb_qs(fb_list_qs *fb, z *n)
{
	//finds the factor base primes, and computes the solutions to the congruence x^2 = N mod p
	//for the QS, these are the starting positions of the sieve relative to the sqrt of N.
	//for the MPQS, additional work using the polynomial coefficents and these congruences 
	//needs to be done to compute the starting positions of the sieve.

	int i;
	uint32 b,j,r,k;
	uint16 prime, root1, root2;
	fp_digit f;
	uint8 logp;

	//the 0th element in the fb is always  2, so start searching with 3
	j=2; i=1;
	while (j<fb->B)
	{
		r = (uint32)zShortMod(n,(fp_digit)spSOEprimes[i]);
		if (r == 0)
		{
			//p divides n, which means it divides the multiplier.
			//we can still use it, but it only has one solution to x^2 == n mod p instead
			//of two.  just divide its logprime in half.
			//we also can't find the root using shanks-tonelli, but it will be very small
			//because the multiplier is very small, so just use brute force.
			prime = (uint16)spSOEprimes[i];
			b = (uint32)zShortMod(n,(fp_digit)prime);
			k=0;
			while (1)
			{
				if (((k*k) % prime) == b)
					break;
				k++;
			}
			root1 = k;
			root2 = prime - k;

			//compute logp
			logp = (uint8)(log((double)prime)/log(2.0) + .5)/2;

			//fill in factor base
			fb->list[j].prime = prime;

			if (root2 > root1)
			{
				fb->list[j].c1 = root1;
				fb->list[j].c2 = root2;
			}
			else
			{
				fb->list[j].c1 = root2;
				fb->list[j].c2 = root1;
			}
			fb->list[j].logprime = logp;
			
			j++;
			i++;
			continue;
		}

		b = jacobi_1((fp_digit)r,(fp_digit)spSOEprimes[i]);
		if (b==1)
		{
			//this prime works
			prime = (uint16)spSOEprimes[i];
			ShanksTonelli_1((fp_digit)r,(fp_digit)prime,&f);
			root1 = (uint32)f;
			root2 = prime - root1;

			//compute logp
			logp = (uint8)(log((double)prime)/log(2.0) + .5);

			//fill in factor base
			fb->list[j].prime = prime;
			if (root2 > root1)
			{
				fb->list[j].c1 = root1;
				fb->list[j].c2 = root2;
			}
			else
			{
				fb->list[j].c1 = root2;
				fb->list[j].c2 = root1;
			}
			fb->list[j].logprime = logp;
			
			j++;
		}
		i++;
	}

	return;
}

void smallmpqs(z *n, z *f1, z *f2, z *f3)
{
	qs_rlist *full, *partial;
	fb_list_qs *fb;
	sieve_fb *fb_sieve_p,*fb_sieve_n;
	qs_poly *poly;
	
	z *factors;
	z tmp, tmp2, tmp3, sqrt_n;
	z *apoly, *bpoly;

	uint32 numpoly, polyalloc;

	uint32 i,j;
	uint8 mul;
	uint16 pmax;							//largest prime in factor base
	uint16 cutoff;
	uint16 sieve_interval;
	uint16 num, max_f;
	int charcount, block;
	uint32 num_factors;
	uint8 *sieve;							//sieve values
	uint8 closnuf;
	uint32 nsieveloc;

	zClear(f1);
	zClear(f2);
	zClear(f3);

	if (n->val[0] == 1 && n->size == 1)
		return;

	if ((n->val[0] % 2) == 0)
		return;

	if (isSquare(n))
	{
		zNroot(n,f1,2);
		zCopy(f1,f2);
		return;
	}

	zInit(&tmp);
	zInit(&tmp2);
	zInit(&tmp3);

	//allocate the space for the factor base
	fb = (fb_list_qs *)malloc(sizeof(fb_list_qs));

	SM_BLOCKSIZE = 512;
	i = zBits(n);
	if (i < 60)
	{
		SM_BLOCKSIZE = 512;
		SM_NUMBLOCKS = 1;
		fb->B = 25;
	}
	else if (i < 80)
	{
		SM_NUMBLOCKS = 2;
		fb->B = 50;
	}
	else if (i < 100)
	{
		SM_NUMBLOCKS = 4;
		fb->B = 100;
	}
	else if (i < 110)
	{
		SM_NUMBLOCKS = 5;
		fb->B = 150;
	}
	else if (i < 120)
	{
		SM_NUMBLOCKS = 6;
		fb->B = 200;
	}
	else
	{
		SM_NUMBLOCKS = 8;
		fb->B = 300;
	}

	//allocate storage for relations based on the factor base size
	max_f = fb->B + NUM_EXTRA_RELS;	//num_f should never get this high, should find at least a few fulls from partials.
	full = (qs_rlist *)malloc((size_t)(sizeof(qs_rlist)));
	full->allocated = max_f;
	full->num_r = 0;
	full->act_r = 0;
	full->list = (qs_r **)malloc((size_t) (max_f * sizeof(qs_r *)));

	//we will typically also generate max_f/2 * 10 partials (empirically determined)
	partial = (qs_rlist *)malloc((size_t)(sizeof(qs_rlist)));
	partial->allocated = 10*fb->B;
	partial->num_r = 0;
	partial->act_r = 0;
	partial->list = (qs_r **)malloc((size_t) (10 * fb->B * sizeof(qs_r *)));

	//set the sieve interval.
	sieve_interval = SM_NUMBLOCKS * SM_BLOCKSIZE;

	//allocate the space for the factor base
	fb->list = (fb_element_qs *)malloc((size_t)(fb->B * sizeof(fb_element_qs)));
	fb_sieve_p = (sieve_fb *)malloc((size_t)(fb->B * sizeof(sieve_fb)));
	fb_sieve_n = (sieve_fb *)malloc((size_t)(fb->B * sizeof(sieve_fb)));
	
	//allocate the sieve
	sieve = (uint8 *)malloc((size_t) (SM_BLOCKSIZE * sizeof(uint8)));

	//allocate the current polynomial
	poly = (qs_poly *)malloc(sizeof(qs_poly));
	zInit(&poly->poly_a);
	zInit(&poly->poly_b);
	zInit(&poly->poly_c);

	//allocate the polynomial lists
	polyalloc = 32;
	apoly = (z *)malloc(polyalloc * sizeof(z));
	bpoly = (z *)malloc(polyalloc * sizeof(z));
	for (i=0;i<polyalloc;i++)
	{
		zInit(&apoly[i]);
		zInit(&bpoly[i]);
	}

	//find multiplier
	mul = (uint8)choose_multiplier(n,fb->B);
	zShortMul(n,mul,&tmp);
	zCopy(&tmp,n);

	//find new sqrt_n
	zInit(&sqrt_n);
	zNroot(n,&sqrt_n,2);

	//compute the first polynominal 'a' value.  we'll need it before creating the factor base in order
	//to find the first roots
	//'a' values should be as close as possible to sqrt(2n)/M, they should be a quadratic residue mod N (d/N) = 1,
	//and be a prime congruent to 3 mod 4.  this last requirement is so that b values can be computed without using the 
	//shanks-tonelli algorithm, and instead use faster methods.
	//since a = d^2, find a d value near to sqrt(sqrt(2n)/M)
	zShiftLeft(&tmp,n,1);
	zNroot(&tmp,&tmp2,2);
	zShortDiv(&tmp2,sieve_interval,&tmp2);

	zNroot(&tmp2,&tmp,2);
	poly->poly_d = (uint32)tmp.val[0];
	if (!(poly->poly_d & 1))
		poly->poly_d++;

	if (poly->poly_d < 5)
		poly->poly_d = 5;
 
	poly->poly_d = sm_nextD(poly->poly_d,n,&tmp);
	sm_computeB(poly,n,&tmp,&tmp2);

	checkpoly_qs(poly,n);

	fb->list[0].prime = 1;
	fb->list[1].prime = 2;

	//construct the factor base, and copy to the sieve factor base
	sm_make_fb_qs(fb,n);
	for (i=2;i<fb->B;i++)
	{
		fb_sieve_p[i].prime = fb->list[i].prime;
		fb_sieve_p[i].logprime = fb->list[i].logprime;
		fb_sieve_n[i].prime = fb->list[i].prime;
		fb_sieve_n[i].logprime = fb->list[i].logprime;
	}

	//find the root locations of the factor base primes for this poly
	sm_computeRoots(poly,fb,mul);

	pmax = fb->list[fb->B-1].prime;
	cutoff = pmax * 30;

	//compute the number of bits in M/2*sqrt(N/2), the approximate value
	//of residues in the sieve interval
	//sieve locations greater than this are worthy of trial dividing
	closnuf = (uint8)(double)((zBits(n) - 1)/2);
	closnuf += (uint8)(log((double)sieve_interval/2)/log(2.0));
	closnuf -= (uint8)(1.3 * log(cutoff) / log(2.0));
	closnuf -=4;	//empirical fudge factor
	
	numpoly = 0;
	num = 0;
	charcount=0;
	nsieveloc = 0;
	while (1)
	{
		//copy current poly into the poly lists
		if (numpoly < polyalloc)
		{
			zCopy(&poly->poly_a,&apoly[numpoly]);
			zCopy(&poly->poly_b,&bpoly[numpoly]);
		}
		else
		{
			polyalloc *= 2;
			apoly = (z *)realloc(apoly,polyalloc * sizeof(z));
			bpoly = (z *)realloc(bpoly,polyalloc * sizeof(z));
			for (i=numpoly; i<polyalloc; i++)
			{
				zInit(&apoly[i]);
				zInit(&bpoly[i]);
			}
			zCopy(&poly->poly_a,&apoly[numpoly]);
			zCopy(&poly->poly_b,&bpoly[numpoly]);
		}

		//copy root info into the fb_sieve
		for (i=2;i<fb->B;i++)
		{
			fb_sieve_n[i].root1 = fb->list[i].nroot1;
			fb_sieve_n[i].root2 = fb->list[i].nroot2;
		}
		for (i=2;i<fb->B;i++)
		{
			fb_sieve_p[i].root1 = fb->list[i].proot1;
			fb_sieve_p[i].root2 = fb->list[i].proot2;
		}
			
		for (block = 0; (uint32)block < SM_NUMBLOCKS; block++)
		{
			sm_sieve_block(sieve,fb_sieve_p,closnuf,fb);
			sm_check_relations(sieve_interval,block,sieve,n,poly,closnuf,
				fb_sieve_p,fb,full,partial,cutoff,0,&num,numpoly);
			nsieveloc += SM_BLOCKSIZE;

			j=0;
			if (partial->num_r > 0)
			{
				//check the partials for full relations
				qsort(partial->list,partial->num_r,sizeof(qs_r *),&sm_qcomp_qs);
				for (i=0;i<partial->num_r-1;i++)
				{
					if (partial->list[i]->largeprime == partial->list[i+1]->largeprime)
						j++;
				}
				partial->act_r = j;
			}
		
			if (j+(full->num_r) >= fb->B + NUM_EXTRA_RELS) 
			{
				//we've got enough total relations to stop
				goto done;
			}

			sm_sieve_block(sieve,fb_sieve_n,closnuf,fb);
			sm_check_relations(sieve_interval,block,sieve,n,poly,closnuf,
				fb_sieve_n,fb,full,partial,cutoff,1,&num,numpoly);
			nsieveloc += SM_BLOCKSIZE;

			j=0;
			if (partial->num_r > 0)
			{
				//check the partials for full relations
				qsort(partial->list,partial->num_r,sizeof(qs_r *),&sm_qcomp_qs);
				for (i=0;i<partial->num_r-1;i++)
				{
					if (partial->list[i]->largeprime == partial->list[i+1]->largeprime)
						j++;
				}
				partial->act_r = j;
			}

			if (j+(full->num_r) >= fb->B + NUM_EXTRA_RELS) 
			{
				//we've got enough total relations to stop
				goto done;
			}
		}

		//next polynomial
		poly->poly_d = sm_nextD(poly->poly_d,n,&tmp);
		sm_computeB(poly,n,&tmp,&tmp2);
		sm_computeRoots(poly,fb,mul);

		checkpoly_qs(poly,n);

		numpoly++;
	}

done:

	//can free sieving structures now
	free(sieve);
	free(fb_sieve_p);
	free(fb_sieve_n);
	zFree(&poly->poly_a);
	zFree(&poly->poly_b);
	zFree(&poly->poly_c);
	free(poly);
	
	num_factors=0;
	factors = (z *)malloc(MAX_FACTORS_SM * sizeof(z));
	for (i=0;i<MAX_FACTORS_SM;i++)
		zInit(&factors[i]);

	i = BlockGauss_smqs(full,partial,apoly,bpoly,fb,n,mul,factors,&num_factors);

	if (num_factors == 1)
	{
		zCopy(&factors[0],f1);
	}
	else if (num_factors == 2)
	{
		zCopy(&factors[0],f1);
		zCopy(&factors[1],f2);
	}
	else if (num_factors == 3)
	{
		zCopy(&factors[0],f1);
		zCopy(&factors[1],f2);
		zCopy(&factors[2],f3);
	}

	zFree(&tmp);
	zFree(&tmp2);
	zFree(&tmp3);
	zFree(&sqrt_n);

	for (i=0;i<MAX_FACTORS_SM;i++)
		zFree(&factors[i]);
	free(factors);

	for (i=0;i<full->num_r;i++)
	{
		free(full->list[i]->fboffset);
		free(full->list[i]);
	}
	free(full->list);
	free(full);

	for (i=0;i<partial->num_r;i++)
	{
		free(partial->list[i]->fboffset);
		free(partial->list[i]);
	}
	free(partial->list);
	free(partial);

	for (i=0;i<polyalloc;i++)
	{
		zFree(&apoly[i]);
		zFree(&bpoly[i]);
	}
	free(apoly);
	free(bpoly);

	free(fb->list);
	free(fb);
	return;
}

int checkpoly_qs(qs_poly *poly, z *n)
{
	//check that b^2 == N mod a
	//and that c == (b*b - n)/a
	z t1,t2,t3,t4;

	zInit(&t1);
	zInit(&t2);
	zInit(&t3);
	zInit(&t4);

	zCopy(n,&t1);
	zDiv(&t1,&poly->poly_a,&t2,&t3);

	zMul(&poly->poly_b,&poly->poly_b,&t2);
	zDiv(&t2,&poly->poly_a,&t1,&t4);

	if (zCompare(&t3,&t4) != 0)
	{
		printf("\nError in checkpoly: %s^2 !== %s mod %s\n",
			z2decstr(&poly->poly_b,&gstr1),
			z2decstr(n,&gstr2),
			z2decstr(&poly->poly_a,&gstr3));
		if (poly->poly_b.size < 0)
			printf("b is negative\n");
	}

	if (zJacobi(n,&poly->poly_a) != 1)
		printf("\nError in checkpoly: (a|N) != 1\n");

	zSqr(&poly->poly_b,&t1);
	zSub(&t1,n,&t2);
	zDiv(&t2,&poly->poly_a,&t1,&t3);

	if (zCompare(&t1,&poly->poly_c) != 0)
		printf("\nError in checkpoly: c != (b^2 - n)/a\n");
	
	zFree(&t1);
	zFree(&t2);
	zFree(&t3);
	zFree(&t4);
	return 0;
}

static int sm_check_relations(uint16 sieve_interval, uint8 blocknum, 
							  uint8 sieve[], z *n, qs_poly *poly, uint8 closnuf,
							  sieve_fb *fb, fb_list_qs *fullfb, qs_rlist *full, 
							  qs_rlist *partial, uint16 cutoff,  
							  uint8 parity, uint16 *num, int numpoly)
{
	z Q,t1,t2,t3;
	uint16 offset,j,k;
	uint8 neg;
	uint64 *sieveblock;
	uint64 mask = ( ((uint64)0x80808080 << 32) | 0x80808080);
	
	sieveblock = (uint64 *)sieve;
	zInit(&Q);
	zInit(&t1);
	zInit(&t2);
	zInit(&t3);
	//check for relations
	for (j=0;j<SM_BLOCKSIZE/8;j++)
	{
		if ((sieveblock[j] & mask) == (uint64)(0))
			continue;

		if (full->num_r > fullfb->B + NUM_EXTRA_RELS - 5)
			break;

		//else figure out which one's need to be checked
		for (k=0;k<8;k++)
		{
			if ((sieve[8*j + k] & 0x80) == 0)
				continue;

			(*num)++;

			//this one is close enough, compute 
			//Q(x) = (ax + b)^2 - N, where x is the sieve index
			//(a*x +/- 2b)*x + c;
			//offset = blocknum*BLOCKSIZE + 8*j+k;
			//printf("computing Q for block %d, offset %d\n",blocknum,8*j+k);
			offset = (blocknum<<16) + (j<<3) + k;
			zShiftLeft(&t2,&poly->poly_b,1);

			zShortMul(&poly->poly_a,offset,&t1);
			if (parity)
				zSub(&t1,&t2,&t3);
			else
				zAdd(&t1,&t2,&t3);

			zShortMul(&t3,offset,&t1);
			zAdd(&t1,&poly->poly_c,&Q);
			if (Q.size < 0)
				neg = 1;
			else
				neg = 0;
			Q.size = abs(Q.size);

			sm_trial_divide_Q(&Q,fb,full,partial,sieve,offset,8*j+k,neg,
				fullfb,cutoff,numpoly,parity,closnuf);
		}
	}
	zFree(&Q);
	zFree(&t1);
	zFree(&t2);
	zFree(&t3);
	return full->num_r;
}

static void sm_trial_divide_Q(z *Q, sieve_fb *fb, qs_rlist *full, qs_rlist *partial,
							  uint8 *sieve, uint16 offset, uint16 j, uint8 sign, 
							  fb_list_qs *fullfb, uint16 cutoff,
							  int numpoly, uint8 parity, uint8 closnuf)
{
	sieve_fb *fbptr;
	uint32 i,num_f,num_p;
	uint16 root1,root2,prime;
	uint32 r;
	int smooth_num;
	uint16 fboffset[MAX_SMOOTH_PRIMES];
	uint8 logp;
	num_f = full->num_r;
	num_p = partial->num_r;
	
	//we have two signs to worry about.  the sign of the offset tells us how to calculate ax + b, while
	//the sign of Q(x) tells us how to factor Q(x) (with or without a factor of -1)
	//the square root phase will need to know both.  fboffset holds the sign of Q(x).  the sign of the 
	//offset is stored standalone in the relation structure.
	if (sign)
		fboffset[0] = 1;
	else
		fboffset[0] = 0;

	smooth_num=0;

	//take care of powers of two
	while (!(Q->val[0] & 1))
	{
		zShiftRight(Q,Q,1);
		fboffset[++smooth_num] = 1;
	}

	i=2;
	while (i < fullfb->B)
	{
		fbptr = fb + i;
		root1 = fbptr->root1 + SM_BLOCKSIZE - j;		//after sieving a block, the root is updated for the start of the next block
		root2 = fbptr->root2 + SM_BLOCKSIZE - j;		//get it back on the current block's progression

		prime = fbptr->prime;
		logp = fbptr->logprime;

		if (Q->size == 1)
		{
			if (Q->val[0] < prime)
				break;
		}

		if (root2 >= prime)
		{
			//r2 is bigger than prime, it could be on the progression, check it.
			if (!(root2 % prime))
			{
				//it is, so it will divide Q(x).  do so as many times as we can.
				do
				{
					fboffset[++smooth_num] = (uint16)i;
					zShortDiv(Q,prime,Q);
					sieve[j] += logp;
					r = zShortMod(Q,prime);
				} while (r == 0);
				if (sieve[j] == closnuf)
					goto done;
			}
			else if ((root1 >= prime) && (!(root1 % prime)))
			{
				//r2 was a bust, but root1 met the criteria.  divide Q(x).		
				do
				{
					fboffset[++smooth_num] = (uint16)i;
					zShortDiv(Q,prime,Q);
					sieve[j] += logp;
					r = zShortMod(Q,prime);
				} while (r == 0);
				if (sieve[j] == closnuf)
					goto done;
			}
		}
		i++;
	}

done:

	//check if it completely factored by looking at the unfactored portion in tmp
	if ((Q->size == 1) && (Q->val[0] == 1))
		sm_save_relation_qs(full,offset,1,smooth_num+1,num_f,fboffset,numpoly,parity);
	else if ((Q->size == 1)  && (Q->val[0] < cutoff))
	{
		sm_save_relation_qs(partial,offset,Q->val[0],smooth_num+1,num_p,fboffset,
			numpoly,parity);
		if (partial->num_r == partial->allocated) 
		{
			partial->allocated *= 2;
			partial->list = (qs_r **)realloc(partial->list, 
					partial->allocated * sizeof(qs_r *));
		}
	}
	return;
}

static void sm_save_relation_qs(qs_rlist *list, uint16 offset, uint16 largeprime, 
								uint8 num_factors, uint16 rnum, uint16 *fboffset, 
								uint8 numpoly, uint8 parity)
{
	uint32 i;
	list->list[rnum] = (qs_r *)malloc(sizeof(qs_r));
	list->list[rnum]->fboffset = (uint16 *)malloc(num_factors*sizeof(uint16));
	for (i=0;i<num_factors;i++)
		list->list[rnum]->fboffset[i] = fboffset[i];
	
	list->list[rnum]->offset = offset;
	list->list[rnum]->largeprime = largeprime;
	list->list[rnum]->parity = parity;
	list->list[rnum]->num_factors = num_factors;
	list->list[rnum]->polynum = numpoly;
	list->num_r++;
	return;
}

static void sm_sieve_block(uint8 *sieve, sieve_fb *fb,
						   uint8 s_init, fb_list_qs *fullfb)
{
	uint32 prime, root1, root2, tmp;
	uint32 B=fullfb->B;
	uint32 j;
	uint8 *inner_sieve;
	uint8 logp;
	sieve_fb *fbptr;

	//initialize block
	memset(sieve,s_init,SM_BLOCKSIZE);

	inner_sieve = sieve;

	//for some small number of primes (those that fit in L1 cache/2
	for (j=2; j<B; j++)
	{
		fbptr = fb + j;
		prime = fbptr->prime;
		root1 = fbptr->root1;
		root2 = fbptr->root2;
		logp = fbptr->logprime;

		while (root2 < SM_BLOCKSIZE)
		{
			inner_sieve[root1] -= logp;
			inner_sieve[root2] -= logp;
			root1 += prime;
			root2 += prime;
		}

		//don't forget the last proot1[i], and compute the roots for the next block
		if (root1 < SM_BLOCKSIZE)
		{
			inner_sieve[root1] -= logp;
			root1 += prime;
			//root1 will be bigger on the next iteration, switch them now
			tmp = root2;
			root2 = root1;
			root1 = tmp;
		}

		fbptr->root1 = root1 - SM_BLOCKSIZE;
		fbptr->root2 = root2 - SM_BLOCKSIZE;
	}

	return;
}

uint32 sm_nextD(uint32 poly_d, z *n, z *work)
{
	//the current poly_d is passed in.
	//compare it to opt_d and find the next poly_d, return in poly_d.
	fp_digit p;

	//the current poly_d is prime, so check and increment to find the next one that meets all the criteria
	do 
	{
		//get the next prime number
		zNextPrime_1(poly_d,&p,work,1);	
		poly_d = (uint32)p;
		sp2z(p,work);
	} while ((zJacobi(n,work) != 1) || ((poly_d & 3) != 3));

	return poly_d;
}

void sm_computeB(qs_poly *poly, z *n, z *t2, z *t3)
{
	//using poly_d, compute poly_b and poly_a = poly_d^2
	fp_digit tmp, t0, h1, h2, t6;
	//poly_a = d^2.  we just found d.  also compute b using Hegel theorem and lifting.

	//n^(d-3)/4 mod d
	tmp = (poly->poly_d - 3) / 4;
	zModExp_1(n,tmp,poly->poly_d,&t0);	

	//h1 = n*t0 mod d
	zShortMul(n,t0,t3);		
	h1 = zShortMod(t3,poly->poly_d);

	//(n - h1^2)/d mod d
	sp2z(h1,t2);
	zSqr(t2,t2);		
	zSub(n,t2,t3);		
	zShortDiv(t3,poly->poly_d,t2);	
	tmp = zShortMod(t2,poly->poly_d);	

	//(2*h1)^-1 mod d = (2*h1)^(d-2) mod d
	sp2z(2*h1,t2);		
	zModExp_1(t2,poly->poly_d - 2,poly->poly_d,&t6);	

	//h2 = ((2*h1)^-1 * (n - h1^2)/d) mod d
	spMulMod(t6,tmp,poly->poly_d,&h2);

	//h1 + h2*D
	sp2z(h2,t2);
	zShortMul(t2,poly->poly_d,t2);
	zShortAdd(t2,h1,t2);

	//we're now done with d, so compute a = d^2
	sp2z(poly->poly_d,t3);
	zSqr(t3,&poly->poly_a);		

	//b = (h1 + h2*d) mod a
	zDiv(t2,&poly->poly_a,t3,&poly->poly_b);	

	//make sure b < a/2
	zShiftRight(t2,&poly->poly_a,1);
	if (zCompare(&poly->poly_b,t2) > 0)
		zSub(&poly->poly_a,&poly->poly_b,&poly->poly_b);

	//now that we have b, compute c = (b*b - n)/a
	zSqr(&poly->poly_b,t2);
	zSub(t2,n,t3);
	zDiv(t3,&poly->poly_a,&poly->poly_c,t2);

	return;
}

void sm_computeRoots(qs_poly *poly, fb_list_qs *fb, uint8 multiplier)
{
	//the roots are computed using a and b as follows:
	//(+/-t - b)(a)^-1 mod p
	//assume b > t
	//uint32 root1, root2, prime, amodp;
	int root1, root2, prime, amodp, bmodp, x;
	uint32 i;

	for (i=2;i<fb->B;i++)
	{
		//fast method of computing the inverse from lenstra...
		root1 = fb->list[i].c1;
		root2 = fb->list[i].c2;
		prime = fb->list[i].prime;

		//find a^-1 mod p = inv(a mod p) mod p
		amodp = zShortMod(&poly->poly_a,prime);
		if (amodp == 0 || multiplier % prime == 0)
		{
			fb->list[i].proot1 = (uint16)-1;
			fb->list[i].nroot1 = (uint16)-1;
			fb->list[i].proot2 = (uint16)-1;
			fb->list[i].nroot2 = (uint16)-1;
			continue;
		}
		amodp = modinv_1(amodp,prime);

		//find (t - b) mod p and (-t - b) mod p
		bmodp = zShortMod(&poly->poly_b,prime);
		x = (int)root1 - bmodp;
		if (x < 0) x += prime;
		root1 = x;
		x = (int)root2 - bmodp;
		if (x < 0) x += prime;
		root2 = x;	

		root1 = (uint32)((uint64)amodp * (uint64)root1 % (uint64)prime);
		root2 = (uint32)((uint64)amodp * (uint64)root2 % (uint64)prime);		

		if (root2 < root1)
		{
			fb->list[i].proot1 = root2;
			fb->list[i].proot2 = root1;
			fb->list[i].nroot1 = prime - root1;
			fb->list[i].nroot2 = prime - root2;
		}
		else
		{
			fb->list[i].proot1 = root1;
			fb->list[i].proot2 = root2;
			fb->list[i].nroot1 = prime - root2;
			fb->list[i].nroot2 = prime - root1;
		}
	}
	return;
}

int sm_qcomp_qs(const void *x, const void *y)
{
	qs_r **xx = (qs_r **)x;
	qs_r **yy = (qs_r **)y;
	
	if (xx[0]->largeprime > yy[0]->largeprime)
		return 1;
	else if (xx[0]->largeprime == yy[0]->largeprime)
		return 0;
	else
		return -1;
}


static uint64 bitValRead64(uint64 **m, int row, int col);

int BlockGauss_smqs(qs_rlist *full, qs_rlist *partial, z *apoly, z *bpoly,
			fb_list_qs *fb, z *n, int mul, 
			z *factors, uint32 *num_factor)
{
	int i,j,k,l,a,q,polynum;
	int *bl;
	uint8 **m;		//matrix of the powers of the prime decompositions of the relations over the factor base
	uint64 **m2_64;	//m mod 2, packed into 32 bit words
	uint64 **aug_64;	//matrix to store the permutations of the rows of m2, packed into 32 bit words
	uint32 largep, bool_val, B = fb->B;
	uint32 *partial_index;
	int num_f,num_p;
	int num_r,num_col,num_col_aug,set_continue;
	const int blocksz = 64;

	uint32 *pd;
	uint32 r;
	z zx, zy, tmp, tmp2, tmp3, tmp4, nn,poly_b,poly_d1,poly_d2,tmp_a;
	
	zInit(&zx);
	zInit(&zy);
	zInit(&tmp);
	zInit(&tmp2);
	zInit(&tmp3);
	zInit(&tmp4);
	zInit(&nn);
	zInit(&poly_b);
	zInit(&poly_d1);
	zInit(&poly_d2);
	zInit(&tmp_a);

	num_f = full->num_r;
	num_p = partial->act_r;

	num_r = full->num_r + partial->act_r;
	num_col = (uint32)((B/blocksz)+1);
	num_col_aug = (uint32)(num_r/blocksz+1);

	//allocate storage based on total number of relations.
	pd = (uint32 *)malloc(B * sizeof(uint32));
	partial_index = (uint32 *)malloc(num_p * sizeof(uint32));

	aug_64 = (uint64 **)malloc(num_r * sizeof(uint64 *));
	for (i=0; i<num_r; i++)
		aug_64[i] = (uint64 *)malloc(num_col_aug * sizeof(uint64));

	m2_64 = (uint64 **)malloc(num_r * sizeof(uint64 *));
	for (i=0; i<num_r; i++)
		m2_64[i] = (uint64 *)malloc(num_col * sizeof(uint64));

	m = (uint8 **)malloc(num_r * sizeof(uint8 *));
	for (i=0; i<num_r; i++)
		m[i] = (uint8 *)malloc(B * sizeof(uint8));

	bl = (int *)malloc(num_r * sizeof(int));

	//write fulls to m
	for (i=0;i<num_f;i++)
	{
		//Initialize
		for (j=0;j<(int)B;j++)
			m[i][j] = 0;

		//copy the pd's of the fboffsets to the correct location in m
		//offset 0 is special - indicates the parity of the offset
		m[i][0] = (uint8)full->list[i]->fboffset[0];
		j=1;
		while (j<full->list[i]->num_factors)
		{
			m[i][full->list[i]->fboffset[j]]++;
			j++;
		}
	}

	//write fulls from partials to m, probably also redundant to do it this way?
	largep = partial->list[0]->largeprime;
	j=num_f;
	for (i=1;i<(int)partial->num_r;i++)
	{
		if (partial->list[i]->largeprime == largep)
		{
			//this partial's largep is the same as the one before, add the pd's and copy to m
			for (k=0;k<(int)B;k++)
				m[j][k] = 0;

			//do the factor of -1
			m[j][0] = (uint8)partial->list[i-1]->fboffset[0];
			//then the rest
			k=1;
			while (k<partial->list[i-1]->num_factors) 
			{
				m[j][partial->list[i-1]->fboffset[k]]++;
				k++;
			}
			//factor of -1
			m[j][0] += partial->list[i]->fboffset[0];
			//the rest
			k=1;
			while (k<partial->list[i]->num_factors)
			{
				m[j][partial->list[i]->fboffset[k]]++;
				k++;
			}
			//remember the index of the partial that made this full relation, we'll need it later
			partial_index[j-num_f]=i;
			//increment the relation counter
			j++;
		}
		largep = partial->list[i]->largeprime;
	}

	//construct the bit matrix
	for (i=0;i<num_r;i++)
	{
		for (j=0;j<num_col;j++)
		{
			m2_64[i][j] = 0;
			for (k=0;k<blocksz;k++)
			{
				if ((blocksz*j+k) < (int)B)
					m2_64[i][j] |= ((uint64)((uint64)m[i][blocksz*j+k]%2) << k);
			}
		}
	}

	//construct augmented matrix
	for (i=0;i<num_r;i++)
	{
		for (j=0;j<num_col_aug;j++)
		{
			aug_64[i][j] = 0;
			for (k=0;k<blocksz;k++)
			{
				if ((blocksz*j+k)==i)
					aug_64[i][j] = ((uint64)(1) << (uint64)(k));
			}
		}
	}

	*num_factor=0;
	//initialize blacklist
	for (i=0;i<num_r;i++) bl[i] = 0;
	//search over all columns, right to left (more sparse on the right side)
	for (i=B-1;i>=0;i--)
	{
		//and all rows
		for (j=0;j<num_r;j++)
		{
			//if the j'th row, i'th bit is 1 and not blacklisted, continue
			bool_val = (bitValRead64(m2_64,j,i) != 0) && (bl[j] == 0);
			if (bool_val)
			{
				//add the j'th row mod 2 to all rows after it with a 1 in the ith column
				for (k=j+1;k<num_r;k++)
				{
					bool_val = (bitValRead64(m2_64,k,i) != 0) && (bl[k] == 0);
					if (bool_val)
					{
						//found one in the k'th row.  add to the j'th row starting at column i.
						//record the addition in the augmented matrix
						for (l=(uint32)(i/blocksz);l>=0;l--)
							m2_64[k][l] = m2_64[j][l] ^ m2_64[k][l];
						for (l=0;l<num_col_aug;l++)
							aug_64[k][l] = aug_64[k][l] ^ aug_64[j][l];
						
						//then check if the row is all zeros
						a=0;
						for (l=(uint32)(i/blocksz);l>=0;l--)
							a = a || m2_64[k][l];

						if (a==0)
						{
							//initialize solution vector
							for (l=0;l<(int)B;l++) pd[l] = 0;

							//found a potential solution. check it.
							for (l=0;l<num_r;l++)
							{
								bool_val = bitValRead64(aug_64,k,l) != 0;
								if (bool_val)
								{
									//then the l'th row of m was involved
									for (q=0;q<(int)B;q++)
										pd[q] += m[l][q];
								}
							}

							//compute x mod n
							zCopy(&zOne,&zy);
							zCopy(&zOne,&zx);
							for (l=0;l<num_r;l++)
							{
								bool_val = bitValRead64(aug_64,k,l) != 0;
								if (bool_val)
								{
									//then the l'th relation is involved in the product of relations
									if (l >= num_f)
									{
										//if l >= num_f, then this row refers to a relation generated from two partials.
										//we'll need to go back to the two partial locations to find the two offsets to
										//multiply together
										//luckily, we've remembered the index in the complete list of partials that 
										//created this full relation
										
										//our relation is of the form (ax + b)^2 == a(ax^2 + 2bx + c) mod n
										//(ax^2 + 2bx + c) is what we trial divided, and we remembered
										//a and b, so we can form the left hand side easily
										
										if (partial->list[partial_index[l-num_f]]->largeprime != partial->list[partial_index[l-num_f]-1]->largeprime)
											printf("ERROR, large primes not equal\n");

										//recreate poly_b from poly_a (store instead??)
										polynum = partial->list[partial_index[l-num_f]]->polynum;
										zCopy(&bpoly[polynum],&poly_b);
										
										//compute Q1(x)
										zShortMul(&apoly[polynum],partial->list[partial_index[l-num_f]]->offset,&tmp);
										zNroot(&apoly[polynum],&poly_d1,2);
										if (partial->list[partial_index[l-num_f]]->parity)
											zSub(&tmp,&poly_b,&tmp2);
										else
											zAdd(&tmp,&poly_b,&tmp2);
										zCopy(&tmp2,&tmp);

										//include 'a'
										zMul(&zy,&poly_d1,&tmp2);
										zDiv(&tmp2,n,&tmp3,&zy);

										//compute Q2(x)
										polynum = partial->list[partial_index[l-num_f]-1]->polynum;
										zCopy(&bpoly[polynum],&poly_b);
										
										zShortMul(&apoly[polynum],partial->list[partial_index[l-num_f]-1]->offset,&tmp2);
										zNroot(&apoly[polynum],&poly_d2,2);
										if (partial->list[partial_index[l-num_f]-1]->parity)
											zSub(&tmp2,&poly_b,&tmp3);
										else
											zAdd(&tmp2,&poly_b,&tmp3);
										zCopy(&tmp3,&tmp2);

										//compute Q(x1)*Q(x2)
										zMul(&tmp,&tmp2,&tmp4);	
										zDiv(&tmp4,n,&tmp2,&tmp3);	//mod n
										zMul(&zx,&tmp3,&tmp);	//accumulate with previous terms
										zDiv(&tmp,n,&tmp2,&zx);	//mod n
										//include the large prime in mp_y
										sp2z(partial->list[partial_index[l-num_f]]->largeprime,&tmp);
										zMul(&tmp,&zy,&tmp2);
										zDiv(&tmp2,n,&tmp3,&zy);

										//include 'a'
										zMul(&zy,&poly_d2,&tmp2);	
										zDiv(&tmp2,n,&tmp3,&zy);
									}
									else
									{
										polynum = full->list[l]->polynum;
										zCopy(&bpoly[polynum],&poly_b);
										
										//compute Q1(x)
										zShortMul(&apoly[polynum],full->list[l]->offset,&tmp);
										zNroot(&apoly[polynum],&poly_d1,2);
										if (full->list[l]->parity)
											zSub(&tmp,&poly_b,&nn);
										else
											zAdd(&tmp,&poly_b,&nn);

										zMul(&zx,&nn,&tmp);	//accumulate with previous terms
										zDiv(&tmp,n,&tmp2,&zx);	//mod n

										zMul(&zy,&poly_d1,&tmp2);	//sqrt(a) = d is part of mp_y
										zDiv(&tmp2,n,&tmp3,&zy);
									}
								}
							}

							//compute y mod n
							//ignore the factor of -1 in this operation
							for (l=1;l<(int)B;l++)
							{
								if (pd[l] > 0)
								{
									sp2z(fb->list[l].prime,&tmp);
									//pd tracks the exponents of the smooth factors.  we know they are all even
									//at this point.  we don't want to compute pd^2, so divide by 2.
									sp2z(pd[l]/2,&tmp2);
									zModExp(&tmp,&tmp2,n,&tmp3);
									zMul(&tmp3,&zy,&tmp4);
									zDiv(&tmp4,n,&tmp2,&zy);
								}
							}

							//split this off into a subroutine... also look for all non-trivial factors if one is composite
							//compute gcd(x-y,n)
							zSub(&zx,&zy,&tmp);
							zLEGCD(&tmp,n,&nn);

							if ((r = (uint32)zShortDiv(&nn,mul,&tmp)) == 0)
							{
								//mul divides this factor
								zCopy(&tmp,&nn);
							}

							if (!(nn.val[0] & 0x1))
							{
								//if it's not odd (factors of 2 creep in there, for some reason
								//remove the 2's
								while (!(nn.val[0] & 0x1))
									zShiftRight(&nn,&nn,1);
							}
							
							
							if ((zCompare(&nn,&zOne) > 0) && (zCompare(n,&nn) > 0))
							{
								//printf("non-trivial factor found = %s\n",z2decstr(&nn,&gstr1));
								zCopy(&nn,&tmp2);
								r = (uint32)zShortDiv(&tmp2,mul,&tmp3);
								if (r == 0)
									zCopy(&tmp3,&nn);

								if (isPrime(&nn))
								{
									//check that we havent' already found this one
									set_continue = 0;
									for (l=0;l<(int)*num_factor;l++)
									{
										if (zCompare(&nn,&factors[l]) == 0)
											set_continue = 1;
									}
									if (set_continue)
										continue;

									zCopy(&nn,&factors[*num_factor]);

									(*num_factor)++;
									if (*num_factor > MAX_FACTORS)
									{
										printf("max number of factors found in block gauss\n");
										goto free;
									}
									//check if we're done by accumulating all factors and comparing to n
									zCopy(&factors[0],&nn);
									for (l=1;l<(int)*num_factor;l++)
									{
										zCopy(&factors[l],&tmp);
										zMul(&tmp,&nn,&tmp2);
										zCopy(&tmp2,&nn);
									}
									if (zBits(&nn) + 10 >= zBits(n))
									{
										//+ 10 accounts for the multiplier in n
										//found all factors, done
										goto free;
									}
								}

								//check the other factor
								zCopy(n,&tmp);
								zDiv(&tmp,&nn,&tmp2,&tmp3);

								//remove the multiplier if necessary
								zCopy(&tmp2,&tmp);
								r = zShortDiv(&tmp,mul,&tmp3);
								if (r == 0)
									zCopy(&tmp3,&tmp);
								else
									zCopy(&tmp2,&tmp);

								if (isPrime(&tmp))
								{
									//check that we havent' already found this one
									set_continue = 0;
									for (l=0;l<(int)*num_factor;l++)
									{
										if (zCompare(&tmp,&factors[l]) == 0)
											set_continue = 1;
									}
									if (set_continue)
										continue;

									zCopy(&tmp,&factors[*num_factor]);

									(*num_factor)++;
									if (*num_factor > MAX_FACTORS)
									{
										printf("max number of factors found in block gauss\n");
										goto free;
									}
									//check if we're done by accumulating all factors and comparing to n
									zCopy(&factors[0],&nn);
									for (l=1;l<(int)*num_factor;l++)
									{
										zCopy(&factors[l],&tmp);
										zMul(&tmp,&nn,&tmp2);
										zCopy(&tmp2,&nn);
									}
									if (zBits(&nn) + 10 >= zBits(n))
									{
										//+ 10 accounts for the multiplier in n
										//found all factors, done
										goto free;
									}
								}
							} //if non-trivial factor
						} //if a == 0
					} //if found in k'th row
				} //add jth row mod 2 to all appropriate rows after it
				//blacklist the j'th row
				bl[j] = 1;
			} //if not blacklisted
		} //for all rows
	} //for all columns

	//printf("matrix exhausted\n");
	r = (uint32)zShortDiv(n,mul,&tmp);
	for (i=0;(uint32)i<*num_factor;i++)
	{
		zCopy(&tmp,&nn);
		zDiv(&nn,&factors[i],&tmp,&tmp2);
	}

free:
	free(pd);
	free(partial_index);
	for (i=0; i<num_r; i++)
		free(aug_64[i]);
	free(aug_64);
	for (i=0; i<num_r; i++)
		free(m2_64[i]);
	free(m2_64);
	for (i=0; i<num_r; i++)
		free(m[i]);
	free(m);
	free(bl);
	zFree(&zx);
	zFree(&zy);
	zFree(&tmp);
	zFree(&tmp2);
	zFree(&tmp3);
	zFree(&tmp4);
	zFree(&nn);
	zFree(&poly_b);
	zFree(&poly_d1);
	zFree(&poly_d2);
	zFree(&tmp_a);
	return 0;
}

static uint64 masks64[64] = {0x1,0x2,0x4,0x8,
							0x10,0x20,0x40,0x80,
							0x100,0x200,0x400,0x800,
							0x1000,0x2000,0x4000,0x8000,
							0x10000,0x20000,0x40000,0x80000,
							0x100000,0x200000,0x400000,0x800000,
							0x1000000,0x2000000,0x4000000,0x8000000,
							0x10000000,0x20000000,0x40000000ULL,0x80000000ULL,
							0x100000000ULL,0x200000000ULL,0x400000000ULL,0x800000000ULL,
							0x1000000000ULL,0x2000000000ULL,0x4000000000ULL,0x8000000000ULL,
							0x10000000000ULL,0x20000000000ULL,0x40000000000ULL,0x80000000000ULL,
							0x100000000000ULL,0x200000000000ULL,0x400000000000ULL,0x800000000000ULL,
							0x1000000000000ULL,0x2000000000000ULL,0x4000000000000ULL,0x8000000000000ULL,
							0x10000000000000ULL,0x20000000000000ULL,0x40000000000000ULL,0x80000000000000ULL,
							0x100000000000000ULL,0x200000000000000ULL,0x400000000000000ULL,0x800000000000000ULL,
							0x1000000000000000ULL,0x2000000000000000ULL,0x4000000000000000ULL,0x8000000000000000ULL};

static uint64 bitValRead64(uint64 **m, int row, int col)
{
	//col is the column in 0 to B-1 representation
	//read the bit in the packed 64 bit representation of the appropriate row
	//don't bother to check the bounds of m w.r.t row and col, assume caller knows what it's doing
	//return 0 if bit not set, 1 << bit offset otherwize
	
	int offset, mcol;
	mcol = col/64;
	offset = (col%64);
	return (m[row][mcol] & masks64[offset]);
}
