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
#include "qs.h"
#include "factor.h"
#include "util.h"

/*
implements the basic quadratic sieve.
inspired by Jason Papadopolous's msieve,
Scott Contini's siqs code, and tidbits from many many other papers
*/
#define SMBLOCKSIZE 32768
#define SM_INNERBLOCKSIZE 16384



static void sieve_block_QS(uint8 *sieve, sieve_fb *fb, uint8 s_init, uint32 B);
static int check_relations_QS(int parity, uint8 *sieve, z *n,uint32 k, z *sq_rt, uint8 closnuf,
					sieve_fb *fb, uint32 B, mpqs_rlist *full, mpqs_rlist *partial, uint32 cutoff, uint32 *num);
static void Qdivide_QS(z *Q, sieve_fb *fb, mpqs_rlist *full, mpqs_rlist *partial, z *sqrt_n,
						  uint8 *sieve, uint32 offset, uint32 j, uint32 sign, uint32 B, uint32 cutoff);
static void save_relation_QS(mpqs_rlist *list, uint32 offset, uint32 largeprime, uint32 num_factors, 
						  uint32 rnum, uint16 *fboffset, z *sqrt_n);
int BlockGauss_smQS(mpqs_rlist *full, mpqs_rlist *partial, z *apoly, z *bpoly,
			fb_list_mpqs *fb, z *n, int mul, 
			FILE *sieve_log,z *factors,uint32 *num_factor);

void pQS(fact_obj_t *fobj)
{
	z *n = &fobj->N;
	mpqs_rlist *full, *partial;
	fb_list_mpqs *fb;
	sieve_fb *fb_sieve_p,*fb_sieve_n;

	z *factors;
	z tmp, tmp2, tmp3, sqrt_n;

	double t_time,ln_n,pQS_fudge=1.1;
	clock_t start;
	uint32 mul,i,j,k;
	uint32 pmax;				//largest prime in factor base
	uint32 cutoff, pQS_large_mult = 30, max_blocks=16000;
	uint32 check_total, check_inc, num, max_f;
	int digits_n;
	int parity;
	uint32 num_factors;
	uint8 *sieve;				//sieve values
	uint8 s_init;				//initial sieve value
	uint8 closnuf;
	FILE *sieve_log;
	uint64 squfof_factor;
	struct timeval myTVstart, myTVend;
	TIME_DIFF *	difference;
	uint32 qsflags = fobj->qs_obj.flags;

	if (qsflags != 12345)
	{
		start = clock();	
		sieve_log = fopen(flogname,"a");
		gettimeofday(&myTVstart, NULL);
	}
	else
		sieve_log = NULL;

	if (isPrime(n))
	{
		n->type = PRP;
		add_to_factor_list(n);
		logprint(sieve_log,"prp%d = %s\n",ndigits(n),z2decstr(n,&gstr1));
		zCopy(&zOne,n);
		fclose(sieve_log);
		return;
	}
	
	zInit(&tmp);
	zInit(&tmp2);
	zInit(&tmp3);

	if ((zBits(n) < 58) && (qsflags != 12345))
	{
		logprint(sieve_log,"Starting SQUFOF on %s\n",z2decstr(n,&gstr1));
		fflush(sieve_log);
		squfof_factor = sp_shanks_loop(n,fobj);
		sp642z(squfof_factor,&tmp);

		if (zCompare(&tmp,&zOne) <= 0)
		{
			//squfof didn't find a factor, give up and don't alter input n.
			zFree(&tmp);
			zFree(&tmp2);
			zFree(&tmp3);
			fclose(sieve_log);
			return;
		}

		//else, report the factor...
		if (isPrime(&tmp))
		{
			tmp.type = PRP;
			logprint(sieve_log,
				"prp%d = %s\n",ndigits(&tmp),z2decstr(&tmp,&gstr1));
		}
		else
		{
			tmp.type = COMPOSITE;
			logprint(sieve_log,
				"C%d = %s\n",ndigits(&tmp),z2decstr(&tmp,&gstr1));
		}

		add_to_factor_list(&tmp);
		
		//... and divide it out of n
		zDiv(n,&tmp,&tmp2,&tmp3);
		zCopy(&tmp2,n);

		//now check if the result is a proper factor
		if (zCompare(n,&zOne) <= 0)
		{
			//it's not, bail
			zFree(&tmp);
			zFree(&tmp2);
			zFree(&tmp3);
			fclose(sieve_log);
			return;
		}

		//else it is, divide it out and report it
		if (isPrime(n))
		{
			n->type = PRP;
			logprint(sieve_log,
				"prp%d = %s\n",ndigits(n),z2decstr(n,&gstr1));
		}
		else
		{
			n->type = COMPOSITE;
			logprint(sieve_log,
				"C%d = %s\n",ndigits(n),z2decstr(n,&gstr1));
		}

		add_to_factor_list(n);

		zCopy(&zOne,n);
		zFree(&tmp);
		zFree(&tmp2);
		zFree(&tmp3);
		fclose(sieve_log);
		return;
	}

	//default mpqs parameters
	sieve_params.fudge_factor = 1.3;
	sieve_params.large_mult = 30;
	sieve_params.num_blocks = 40;
	sieve_params.num_extra_relations = 16;
	sieve_params.small_limit = 100;

	zNroot(n,&tmp,2);
	zMul(&tmp,&tmp,&tmp2);
	if (zCompare(&tmp2,n) == 0)
	{
		if (isPrime(&tmp))
		{
			tmp.type = PRP;
			add_to_factor_list(&tmp);
			add_to_factor_list(&tmp);
		}
		else
		{
			tmp.type = COMPOSITE;
			add_to_factor_list(&tmp);
			add_to_factor_list(&tmp);
		}
		zCopy(&zOne,n);
		zFree(&tmp);
		zFree(&tmp2);
		zFree(&tmp3);
		return;
	}

	if (ndigits(n) > 50)
	{
		printf("n too big for QS, try MPQS or SIQS\n");
		zFree(&tmp);
		zFree(&tmp2);
		zFree(&tmp3);
		fclose(sieve_log);
		return;
	}

	if (VFLAG > 0 && (qsflags != 12345))
		printf("\nstarting pQS on c(%d): %s\n",ndigits(n),z2decstr(n,&gstr1));

	if (qsflags != 12345)
		logprint(sieve_log,"starting pQS on c(%d): %s\n",ndigits(n),z2decstr(n,&gstr1));

	//allocate the space for the factor base
	fb = (fb_list_mpqs *)malloc(sizeof(fb_list_mpqs));

	//calculate approximate factor base bound.  Should do some experiments to determine if this really is optimal.
	//I divide by two because only approximately half of the primes up to the bound will be quadratic residues
	//mod n, and thus be eligible for the factor base.
	i = zBits(n) - 1;
	ln_n = (double)i * log(2.0);
	//calculate the optimal largest prime in the factor base. 
	fb->B = (uint32)(exp(.5 * sqrt(ln_n * log(ln_n))));
	//swag at better limit
	//fb->B = (uint32)((double)fb->B / 2.8);
	//determine the number of primes we'll need such that the largest prime is close to optimal
	for (i=0;i<fb->B;i++)
	{
		if (spSOEprimes[i] > fb->B)
			break;
	}
	//because only about half the primes will be suitable
	fb->B = i/2;

	//fb->B = 64;
	

	//B should always be < 2048 in this routine.  if not, n is too big and mpqs should be used 
	if (fb->B>2048)
		printf("warning: n too big, should use mpqs\n");

	//compute the number of digits in n and how often to check our list of partial relations
	//and update the gui.
	digits_n = ndigits(n);
	check_inc = fb->B/20;
	check_total = check_inc;

	//allocate storage for relations based on the factor base size
	max_f = 4 * fb->B + sieve_params.num_extra_relations;	//num_f should never get this high, should find at least a few fulls from partials.
	full = (mpqs_rlist *)malloc((size_t)(sizeof(mpqs_rlist)));
	full->allocated = max_f;
	full->num_r = 0;
	full->act_r = 0;
	full->list = (mpqs_r **)malloc((size_t) (max_f * sizeof(mpqs_r *)));

	//we will typically also generate max_f/2 * 10 partials (empirically determined)
	partial = (mpqs_rlist *)malloc((size_t)(sizeof(mpqs_rlist)));
	partial->allocated = 10*fb->B;
	partial->num_r = 0;
	partial->act_r = 0;
	partial->list = (mpqs_r **)malloc((size_t) (10*fb->B* sizeof(mpqs_r *)));

	//allocate the space for the factor base
	fb->list = (fb_element_mpqs *)malloc((size_t)(fb->B * sizeof(fb_element_mpqs)));
	fb_sieve_p = (sieve_fb *)malloc((size_t)(fb->B * sizeof(sieve_fb)));
	fb_sieve_n = (sieve_fb *)malloc((size_t)(fb->B * sizeof(sieve_fb)));
	
	//allocate the sieve
	sieve = (uint8 *)malloc((size_t) (SMBLOCKSIZE * sizeof(uint8)));

	//find multiplier
	mul = (uint32)choose_multiplier(n,fb->B);
	zShortMul(n,mul,&tmp);
	zCopy(&tmp,n);

	//find new sqrt_n
	zInit(&sqrt_n);
	zNroot(n,&sqrt_n,2);

	fb->list[0].prime = 1;
	fb->list[1].prime = 2;

	//construct the factor base
	make_fb(fb,fb->B,n,&sqrt_n);

	pmax = fb->list[fb->B-1].prime;
	cutoff = pmax * pQS_large_mult;

	//no need for small prime variation - too hard to tune cutoff adjustment for these small numbers
	s_init=0;
	
	closnuf = (uint8)(double)((zBits(n) - 1)/2);
	closnuf += 16;
	closnuf -= (uint8)(pQS_fudge * log(cutoff) / log(2.0));

	//print some info to the screen and the log file
	if (VFLAG > 0 && (qsflags != 12345))
	{
		printf("s_init = %d, pmax = %d, cutoff = %u(%d * pmax)\n",s_init,pmax,cutoff,pQS_large_mult);
		printf("inner_block_size = %d, blocksize = %d\n",SM_INNERBLOCKSIZE,SMBLOCKSIZE);
		printf("n = %d digits, %d bits\n",digits_n,zBits(n));
		printf("B = %d, T = %1.2f, k = %d\n",fb->B,pQS_fudge,mul);
		printf("base closnuf = %d\n",closnuf);
	}

	if (qsflags != 12345)
	{
		logprint(sieve_log,"s_init = %d, pmax = %d, cutoff = %u(%d * pmax)\n",s_init,pmax,cutoff,pQS_large_mult);
		logprint(sieve_log,"n = %d digits, %d bits\n",digits_n,zBits(n));
		logprint(sieve_log,"B = %d, T = %1.2f, k = %d, blocksize = %d\n",fb->B,pQS_fudge,mul,SMBLOCKSIZE);
	}
	k = 0;
	num = 0;

	//copy root info into the fb_sieve
	for (i=2;i<fb->B;i++)
	{
		fb_sieve_n[i].prime = fb->list[i].prime;
		fb_sieve_n[i].logprime = fb->list[i].logprime;
		fb_sieve_n[i].root1 = fb->list[i].nroot1;
		fb_sieve_n[i].root2 = fb->list[i].nroot2;
	}
	for (i=2;i<fb->B;i++)
	{
		fb_sieve_p[i].prime = fb->list[i].prime;
		fb_sieve_p[i].logprime = fb->list[i].logprime;
		fb_sieve_p[i].root1 = fb->list[i].proot1;
		fb_sieve_p[i].root2 = fb->list[i].proot2;
	}

	//the sieve and relation timings take almost no time compared to pQS, so keep em in there.
	while (k < max_blocks)
	{
		//positive direction
		parity=1;
		sieve_block_QS(sieve,fb_sieve_p,s_init,fb->B);
		i = check_relations_QS(parity,sieve,n,k,&sqrt_n,closnuf,fb_sieve_p,fb->B,full,partial,cutoff,&num);

		//negative direction
		parity=-1;
		sieve_block_QS(sieve,fb_sieve_n,s_init,fb->B);
		i = check_relations_QS(parity,sieve,n,k,&sqrt_n,closnuf,fb_sieve_n,fb->B,full,partial,cutoff,&num);

		if (i >= check_total || full->num_r > fb->B - 2)
		{
			//check the partials for full relations
			qsort(partial->list,partial->num_r,sizeof(mpqs_r *),&qcomp_mpqs);
			j=0;
			for (i=0;i<partial->num_r-1;i++)
			{
				if (partial->list[i]->largeprime == partial->list[i+1]->largeprime)
					j++;
			}
			partial->act_r = j;

			if ((j + full->num_r) >= fb->B + sieve_params.num_extra_relations) 
				break;
			else
			{
				//need to keep sieving
				if (((fb->B + sieve_params.num_extra_relations) - (j + full->num_r)) < check_inc) 
					check_total += ((fb->B+ sieve_params.num_extra_relations) - (j + full->num_r))/2;
				else
					check_total += check_inc;
			}
		}

		//next block
		k++;
	}

	if (qsflags != 12345)
	{
		logprint(sieve_log,"%d relations found: %d full + %d from %d partial, using %d blocks\n",partial->act_r+full->num_r,full->num_r,partial->act_r,partial->num_r,k);
		logprint(sieve_log,"%d relations checked for completeness\n",num);
	}

	if (VFLAG > 0 && (qsflags != 12345))
		printf("%d(%d) relations found: %d full + %d from %d partial, using %d blocks\n",partial->act_r+full->num_r,fb->B,full->num_r,partial->act_r,partial->num_r,k);

	num_factors=0;
	factors = (z *)malloc(MAX_FACTORS * sizeof(z));
	for (i=0;i<MAX_FACTORS;i++)
		zInit(&factors[i]);

	i = BlockGauss_smQS(full,partial,&sqrt_n,NULL,fb,n,mul,sieve_log,factors,&num_factors);
	//stop = clock();
	//t_time = (double)(stop - start)/(double)CLOCKS_PER_SEC;

	if (qsflags != 12345)
	{
		gettimeofday(&myTVend, NULL);
		difference = my_difftime (&myTVstart, &myTVend);

		t_time = 
			((double)difference->secs + (double)difference->usecs / 1000000);
		free(difference);

		if (VFLAG > 0)
			printf("Total pQS time = %6.4f seconds.\n",t_time);
		logprint(sieve_log,"Total QS time = %6.4f seconds.\n",t_time);
	}

	sp2z(mul,&tmp);
	zDiv(n,&tmp,&tmp2,&tmp3);
	zCopy(&tmp2,n);
	j=0;
	fobj->qs_obj.num_factors = 0;
	for(i=0;i<num_factors;i++)
	{
		if (zShortMod(&factors[i],mul) == 0)
			zShortDiv(&factors[i],mul,&factors[i]);

		sp2z(mul,&tmp);
		if (zCompare(&factors[i],&tmp) == 0)
			continue;

		zDiv(n,&factors[i],&tmp,&tmp2);
		zCopy(&tmp,n);

		if (qsflags != 12345)
		{
			//non-trivial factor found
			if (isPrime(&factors[i]))
			{
				factors[i].type = PRP;
				add_to_factor_list(&factors[i]);
			}
			else
			{
				factors[i].type = COMPOSITE;
				add_to_factor_list(&factors[i]);
			}
		}
		else
		{
			zInit(&fobj->qs_obj.factors[j]);
			zCopy(&factors[i],&fobj->qs_obj.factors[j++]);
			fobj->qs_obj.num_factors++;
		}
			
	}

	for (i=0;i<MAX_FACTORS;i++)
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
	free(fb->list);
	free(fb);
	free(fb_sieve_p);
	free(fb_sieve_n);
	free(sieve);
	zFree(&tmp);
	zFree(&tmp2);
	zFree(&tmp3);
	zFree(&sqrt_n);
	if (qsflags != 12345)
		fclose(sieve_log);
	
	return;
}

static int check_relations_QS(int parity, uint8 *sieve, z *n, uint32 k, z *sqrt_n, uint8 closnuf,
						sieve_fb *fb, uint32 B, mpqs_rlist *full, mpqs_rlist *partial, 
						uint32 cutoff, uint32 *num)
{
	z Q, t2, t3;
	uint32 offset,j;
	uint8 blockbits;
	uint32 neg;
	
	//choosing the close enough value is very important.  to small, and too much time is spent trial factoring, while too
	//big, and many smooths will be lost.
	//scale our fudge factor based on the largest prime in the factor base as well as how far we are in the sieving.
	
	zInit(&Q);
	zInit(&t2);
	zInit(&t3);

	//compensate for distance away from origin
	blockbits=0;
	while (((uint32)(1) << blockbits) < k)
		blockbits++;
	closnuf += blockbits;

	//check for relations
	for (j=0;j<SMBLOCKSIZE;j++)
	{
		//look for f(r)'s that are close enough
		if (sieve[j] < closnuf)
			continue;

		//this check is due to cases where for very small N, too many
		//fulls are found, and overflow.  For the size of N at which
		//this matters, this check is negligible for speed.
		if (full->num_r >= full->allocated)
			break;

		(*num)++;
		//this one is close enough, compute its offset from sqrt_n
		//nn = (sq_rt + i)^2 - n
		offset = k*SMBLOCKSIZE + j;
		if (parity < 0)
		{
			//found during negative branch search...
			zShortSub(sqrt_n,offset,&t2);
			zMul(&t2,&t2,&t3);
			zSub(&t3,n,&Q);
			neg=1;
		}
		else
		{
			//found during positive branch search...
			zShortAdd(sqrt_n,offset,&t2);
			zMul(&t2,&t2,&t3);
			zSub(&t3,n,&Q);
			neg=0;
		}
		Q.size = abs(Q.size);
		Qdivide_QS(&Q,fb,full,partial,sqrt_n,sieve,offset,j,neg,B,cutoff);
	}

	zFree(&Q);
	zFree(&t2);
	zFree(&t3);
	return full->num_r;
}

static void Qdivide_QS(z *Q, sieve_fb *fb, mpqs_rlist *full, mpqs_rlist *partial, z *sqrt_n,
						  uint8 *sieve, uint32 offset, uint32 j, uint32 sign, uint32 B, uint32 cutoff)
{
	sieve_fb *fbptr;
	uint32 i,num_f,num_p;
	uint32 root1,root2,prime;
	uint32 r;
	int smooth_num;
	uint16 fboffset[MAX_SMOOTH_PRIMES];
	uint8 logp;

	num_f = full->num_r;
	num_p = partial->num_r;

	if (sign)
		fboffset[0] = 1;
	else
		fboffset[0] = 0;

	//take care of powers of two first
	smooth_num=0;
	while (!(Q->val[0] & 1))
	{
		zShiftRight(Q,Q,1);
		fboffset[++smooth_num] = 1;
	}

	i=2;
	//do the primes less than the blocksize.  primes bigger than the blocksize can be handled
	//even more efficiently.
	//a couple of observations from jasonp:
	//if a prime divides Q(x), then this index (j) and either
	//root1 or root2 are on the same arithmetic progression.  this we can
	//test with a single precision mod operation
	while ((!((Q->size == 1) && (Q->val[0] == 1))) && (i < B))
	{
		fbptr = fb + i;
		root1 = fbptr->root1 + SMBLOCKSIZE - j;	
		root2 = fbptr->root2 + SMBLOCKSIZE - j;
		prime = fbptr->prime;
		logp = fbptr->logprime;

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
					sieve[j] -= logp;
					if (sieve[j] == 0)
						goto done;
					r = zShortMod(Q,prime);
				} while (r == 0);
				
			}
			else if ((root1 >= prime) && (!(root1 % prime)))
			{
				//r2 was a bust, but root1 met the criteria.  divide Q(x).		
				do
				{
					fboffset[++smooth_num] = (uint16)i;
					zShortDiv(Q,prime,Q);
					sieve[j] -= logp;
					if (sieve[j] == 0)
						goto done;
					r = zShortMod(Q,prime);
				} while (r == 0);
			}
		}
		i++;
	}

done:
	//check if it completely factored by looking at the unfactored portion in tmp
	if ((Q->size == 1) && (Q->val[0] == 1 ))
		save_relation_QS(full,offset,1,smooth_num+1,num_f,fboffset,sqrt_n);
	else if ((Q->size == 1) && (Q->val[0] < cutoff))
		save_relation_QS(partial,offset,Q->val[0],smooth_num+1,num_p,fboffset,sqrt_n);

	return;
}

static void save_relation_QS(mpqs_rlist *list, uint32 offset, uint32 largeprime, uint32 num_factors, 
						  uint32 rnum, uint16 *fboffset, z *sqrt_n)
{
	uint32 i;

	list->list[rnum] = (mpqs_r *)malloc(sizeof(mpqs_r));
	list->list[rnum]->fboffset = (uint16 *)malloc(num_factors*sizeof(uint16));
	
	for (i=0;i<num_factors;i++)
		list->list[rnum]->fboffset[i] = fboffset[i];
	
	list->list[rnum]->offset = offset;
	list->list[rnum]->largeprime = largeprime;
	list->list[rnum]->num_factors = (uint8)(num_factors);
	list->num_r++;

	return;
}

static void sieve_block_QS(uint8 *sieve, sieve_fb *fb, uint8 s_init, uint32 B)
{
	uint32 prime, root1, root2;
	uint32 i,j;
	uint8 *inner_sieve;
	uint8 logp;
	sieve_fb *fbptr;

	//initialize block
	memset(sieve,0,SMBLOCKSIZE);

	//get pointers pointed in the right direction
	inner_sieve = sieve;

	//break block into sections to keep in L1 cache
	//blocksize should be a multiple of inner_blocksize
	for (i=0;i<SMBLOCKSIZE;i+=SM_INNERBLOCKSIZE)
	{
		//for some small number of primes (those that fit in L1 cache/2)
		for (j=2; j<B; j++)
		{
			fbptr = fb + j;
			prime = fbptr->prime;
			root1 = fbptr->root1;
			root2 = fbptr->root2;
			logp = fbptr->logprime;

			while (root2 < SM_INNERBLOCKSIZE)
			{
				inner_sieve[root1] += logp;
				inner_sieve[root2] += logp;
				root1 += prime;
				root2 += prime;
			}

			//don't forget the last proot1[i], and compute the roots for the next block
			if (root1 < SM_INNERBLOCKSIZE)
			{
				inner_sieve[root1] += logp;
				root1 += prime;
				//root1 will be bigger on the next iteration, switch them now
				fbptr->root1 = root2 - SM_INNERBLOCKSIZE;
				fbptr->root2 = root1 - SM_INNERBLOCKSIZE;
			}
			else
			{
				fbptr->root1 = root1 - SM_INNERBLOCKSIZE;
				fbptr->root2 = root2 - SM_INNERBLOCKSIZE;
			}
		}
		//move inner_sieve to the right to paint a new stripe of the small fb primes in sieve
		inner_sieve += SM_INNERBLOCKSIZE;
	}

	return;
}

void make_fb(fb_list_mpqs *fb, uint32 B, z *n, z *sqrt_n)
{
	int i;
	uint32 b,j,r,k;
	uint32 prime, root1, root2, sqrtn_modp;
	fp_digit f;
	uint8 logp;
	
	//the 0th element in the fb is always -1, and the 1st is always 2, so start searching with 3
	j=2; i=1;
	while (j<B)
	{
		r = zShortMod(n,(fp_digit)spSOEprimes[i]);
		if (r == 0)
		{
			//p divides n, which means it divides the multiplier.
			//we can still use it, but it only has one solution to x^2 == n mod p instead
			//of two.  just divide its logprime in half.
			//we also can't find the root using shanks-tonelli, but it will be very small
			//because the multiplier is very small, so just use brute force.
			prime = (uint32)spSOEprimes[i];
			b = zShortMod(n,(fp_digit)prime);
			k=0;
			while (1)
			{
				if (((k*k) % prime) == b)
					break;
				k++;
			}
			root1 = k;
			root2 = prime - k;

			//compute starting roots from sqrt_n in the + and - directions
			sqrtn_modp = zShortMod(sqrt_n,(fp_digit)prime);
			sqrtn_modp = prime - sqrtn_modp;

			root1 += sqrtn_modp;
			if (root1 >= prime)
				root1 = root1 - prime;

			root2 += sqrtn_modp;
			if (root2 >= prime)
				root2 = root2 - prime;

			//compute logp
			logp = (uint8)(log((double)prime)/log(2.0) + .5)/2;

			//fill in factor base
			fb->list[j].prime = prime;
			if (root2 > root1)
			{
				fb->list[j].proot1 = root1;
				fb->list[j].proot2 = root2;
				fb->list[j].nroot1 = prime - root2;
				fb->list[j].nroot2 = prime - root1;
			}
			else
			{
				fb->list[j].proot1 = root2;
				fb->list[j].proot2 = root1;
				fb->list[j].nroot1 = prime - root1;
				fb->list[j].nroot2 = prime - root2;
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
			prime = (uint32)spSOEprimes[i];
			ShanksTonelli_1((fp_digit)r,(fp_digit)prime,&f);
			root1 = (uint32)f;
			root2 = prime - root1;

			//compute starting roots from sqrt_n in the + and - directions
			sqrtn_modp = zShortMod(sqrt_n,(fp_digit)prime);
			sqrtn_modp = prime - sqrtn_modp;

			root1 += sqrtn_modp;
			if (root1 >= prime)
				root1 = root1 - prime;

			root2 += sqrtn_modp;
			if (root2 >= prime)
				root2 = root2 - prime;

			//compute logp
			logp = (uint8)(log((double)prime)/log(2.0) + .5);

			//fill in factor base
			fb->list[j].prime = prime;
			if (root2 > root1)
			{
				fb->list[j].proot1 = root1;
				fb->list[j].proot2 = root2;
				fb->list[j].nroot1 = prime - root2;
				fb->list[j].nroot2 = prime - root1;
			}
			else
			{
				fb->list[j].proot1 = root2;
				fb->list[j].proot2 = root1;
				fb->list[j].nroot1 = prime - root1;
				fb->list[j].nroot2 = prime - root2;
			}
			fb->list[j].logprime = logp;
			j++;
		}
		i++;
	}

	return;
}


static uint64 qs_bitValRead64(uint64 **m, int row, int col);

int BlockGauss_smQS(mpqs_rlist *full, mpqs_rlist *partial, z *apoly, z *bpoly,
			fb_list_mpqs *fb, z *n, int mul, 
			FILE *sieve_log,z *factors,uint32 *num_factor)
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
			bool_val = (qs_bitValRead64(m2_64,j,i) != 0) && (bl[j] == 0);
			if (bool_val)
			{
				//add the j'th row mod 2 to all rows after it with a 1 in the ith column
				for (k=j+1;k<num_r;k++)
				{
					bool_val = (qs_bitValRead64(m2_64,k,i) != 0) && (bl[k] == 0);
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
								bool_val = qs_bitValRead64(aug_64,k,l) != 0;
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
								bool_val = qs_bitValRead64(aug_64,k,l) != 0;
								if (bool_val)
								{									
									if (l >= num_f)
									{
										//if l >= num_f, then this row refers to a relation generated from two partials.
										//we'll need to go back to the two partial locations to find the two offsets to
										//multiply together
										//luckily, we've remembered the index in the complete list of partials that 
										//created this full relation
										
										//sqrt_n is stored in poly_a
										//polynum = partial->list[partial_index[l-num_f]]->polynum;
										//zCopy(&partial->list[partial_index[l-num_f]]->poly_a,&tmp_a);

										if (partial->list[partial_index[l-num_f]]->fboffset[0] == 0)
											zShortAdd(&apoly[0],partial->list[partial_index[l-num_f]]->offset,&nn);		//compute Q(x1)
										else
											zShortSub(&apoly[0],partial->list[partial_index[l-num_f]]->offset,&nn);

										//sqrt_n is stored in poly_a
										//zCopy(&partial->list[partial_index[l-num_f]-1]->poly_a,&tmp_a);
										//polynum = partial->list[partial_index[l-num_f]-1]->polynum;

										if (partial->list[partial_index[l-num_f]-1]->fboffset[0] == 0)
											zShortAdd(&apoly[0],partial->list[partial_index[l-num_f]-1]->offset,&tmp3);		//compute Q(x2)
										else
											zShortSub(&apoly[0],partial->list[partial_index[l-num_f]-1]->offset,&tmp3);

										zMul(&nn,&tmp3,&tmp4);	//compute Q(x1)*Q(x2)
										zMul(&zx,&tmp4,&tmp);	//accumulate with previous terms
										zDiv(&tmp,n,&tmp2,&zx);	//mod n
										
										//include the large prime in mp_y
										zShortMul(&zy,partial->list[partial_index[l-num_f]]->largeprime,&tmp2);
										zDiv(&tmp2,n,&tmp3,&zy);
									}
									else
									{
										//polynum = full->list[l]->polynum;
										if (full->list[l]->fboffset[0] == 0)
											zShortAdd(&apoly[0],full->list[l]->offset,&nn);
										else
											zShortSub(&apoly[0],full->list[l]->offset,&nn);

										zMul(&zx,&nn,&tmp);	//accumulate with previous terms
										zDiv(&tmp,n,&tmp2,&zx);	//mod n
									}																	
								}
							}

							//printf("attempting to factor\n");
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
									if (sieve_log != NULL)
										logprint(sieve_log,"prp%d = %s\n",ndigits(&nn),z2decstr(&nn,&gstr1));

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
									if (sieve_log != NULL)
										logprint(sieve_log,"prp%d = %s\n",ndigits(&tmp),z2decstr(&tmp,&gstr1));

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

	printf("matrix exhausted\n");
	r = (uint32)zShortDiv(n,mul,&tmp);
	for (i=0;(uint32)i<*num_factor;i++)
	{
		zCopy(&tmp,&nn);
		zDiv(&nn,&factors[i],&tmp,&tmp2);
	}
	if (sieve_log != NULL)
		logprint(sieve_log,"c%d = %s\n",ndigits(&tmp),z2decstr(&tmp,&gstr1));

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
	
	/*
	free_lmatrix64(aug_64,0,num_r,0,num_col_aug-1);
	free_lmatrix64(m2_64,0,num_r,0,num_col-1);
	free_cmatrix(m,0,num_r,0,B-1);
	*/
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

static uint64 qs_masks64[64] = {0x1,0x2,0x4,0x8,
							0x10,0x20,0x40,0x80,
							0x100,0x200,0x400,0x800,
							0x1000,0x2000,0x4000,0x8000,
							0x10000,0x20000,0x40000,0x80000,
							0x100000,0x200000,0x400000,0x800000,
							0x1000000,0x2000000,0x4000000,0x8000000,
							0x10000000,0x20000000,0x40000000,0x80000000,
							0x100000000,0x200000000,0x400000000,0x800000000,
							0x1000000000,0x2000000000,0x4000000000,0x8000000000,
							0x10000000000,0x20000000000,0x40000000000,0x80000000000,
							0x100000000000,0x200000000000,0x400000000000,0x800000000000,
							0x1000000000000,0x2000000000000,0x4000000000000,0x8000000000000,
							0x10000000000000,0x20000000000000,0x40000000000000,0x80000000000000,
							0x100000000000000,0x200000000000000,0x400000000000000,0x800000000000000,
							0x1000000000000000,0x2000000000000000,0x4000000000000000,0x8000000000000000};

static uint64 qs_bitValRead64(uint64 **m, int row, int col)
{
	//col is the column in 0 to B-1 representation
	//read the bit in the packed 64 bit representation of the appropriate row
	//don't bother to check the bounds of m w.r.t row and col, assume caller knows what it's doing
	//return 0 if bit not set, 1 << bit offset otherwize
	int offset, mcol;
	mcol = col/64;
	offset = (col%64);
	return (m[row][mcol] & qs_masks64[offset]);
}
