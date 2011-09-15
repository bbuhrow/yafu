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
#include "factor.h"
#include "soe.h"
#include "gmp_xface.h"

void zTrial(fact_obj_t *fobj)
{
	//trial divide n using primes below limit. optionally, print factors found
	uint32 r,k=0;
	uint32 limit = fobj->div_obj.limit;
	int print = fobj->div_obj.print;
	fp_digit q;
	z *n = &fobj->div_obj.n;
	z tmp;

	zInit(&tmp);
	mp2gmp(n, fobj->div_obj.gmp_n);

	if (P_MAX < limit)
		GetPRIMESRange(0,limit);

	while ((mpz_cmp_ui(fobj->div_obj.gmp_n, 1) > 0) && (PRIMES[k] < limit))
	{
		if (k >= NUM_P)
		{
			GetPRIMESRange(PRIMES[k-1]-1,PRIMES[k-1]+10000000);
			k=0;
		}

		q = (fp_digit)PRIMES[k];
		r = mpz_tdiv_ui(fobj->div_obj.gmp_n, q);
		
		if (r != 0)
			k++;
		else
		{
			mpz_tdiv_q_ui(fobj->div_obj.gmp_n, fobj->div_obj.gmp_n, q);
			sp2z(q,&tmp);
			tmp.type = PRIME;
			add_to_factor_list(fobj, &tmp);
			if (print && (VFLAG > 0))
#if BITS_PER_DIGIT == 64
				printf("div: found prime factor = %" PRIu64 "\n",q);
#else
				printf("div: found prime factor = %u\n",q);
#endif
		}
	}

	gmp2mp(fobj->div_obj.gmp_n, n);

	zFree(&tmp);
}

void zFermat(fp_digit limit, fact_obj_t *fobj)
{
	//	  Fermat's factorization method (wikipedia psuedo-code)
	z *n = &fobj->div_obj.n;
	mpz_t a, b2, tmp, tmp2, tmp3, maxa;
	int i;
	int numChars;
	fp_digit reportIt, reportInc;
	fp_digit count;

	mp2gmp(n, fobj->div_obj.gmp_n);

	if (mpz_even_p(fobj->div_obj.gmp_n))
	{
		zShiftRight(n,n,1);
		add_to_factor_list(fobj, &zTwo);
		return;
	}

	if (mpz_perfect_square_p(fobj->div_obj.gmp_n))
	{
		z f;
		zInit(&f);

		mpz_sqrt(fobj->div_obj.gmp_n, fobj->div_obj.gmp_n);
		gmp2mp(fobj->div_obj.gmp_n, &f);
		add_to_factor_list(fobj, &f);
		add_to_factor_list(fobj, &f);
		mpz_set_ui(fobj->div_obj.gmp_n, 1);
		zFree(&f);
		gmp2mp(fobj->div_obj.gmp_n, n);
		return;
	}

	mpz_init(a);
	mpz_init(b2);
	mpz_init(tmp);
	mpz_init(maxa);

	mpz_add_ui(maxa, fobj->div_obj.gmp_n, 1);
	mpz_tdiv_q_2exp(maxa, maxa, 1);

	mpz_sqrt(a, fobj->div_obj.gmp_n);
	mpz_add_ui(a, a, 1);
	mpz_mul(tmp, a, a);
	mpz_sub(b2, tmp, fobj->div_obj.gmp_n);

	count = 0;
	numChars = 0;
	reportIt = limit / 100;
	reportInc = reportIt;
	while((mpz_size(b2) > 0) && (!mpz_perfect_square_p(b2)))
	{
		// todo: special case this...
		mpz_add_ui(a, a, 1); 

		if (mpz_cmp(maxa, a) > 0)
			break;	//give up

		//b2 = a*a - N = b2 + 2*a - 1
		mpz_mul_2exp(tmp, a, 1);
		mpz_add(b2, b2, tmp);

		// todo: special case this...
		mpz_sub_ui(b2, b2, 1);

		count++;
		if (count > limit)
			break;

		//progress report
		if ((count > reportIt) && (VFLAG > 1))
		{
			for (i=0; i< numChars; i++)
				printf("\b");
			numChars = printf("%" PRIu64 "%%",(uint64)((double)count / (double)limit * 100));
			fflush(stdout);
			reportIt += reportInc;
		}
	}

	if ((mpz_size(b2) > 0) && mpz_perfect_square_p(b2))
	{
		//printf("found square at count = %d: a = %s, b2 = %s",count,
		//	z2decstr(&a,&gstr1),z2decstr(&b2,&gstr2));
		z tmpz;
		zInit(&tmpz);

		mpz_sqrt(tmp, b2); 
		mpz_add(tmp, a, tmp);
		gmp2mp(tmp, &tmpz);

		if (mpz_probab_prime_p(tmp, NUM_WITNESSES))
			tmpz.type = PRP;			
		else
			tmpz.type = PRP;
		add_to_factor_list(fobj, &tmpz);

		mpz_tdiv_q(fobj->div_obj.gmp_n, fobj->div_obj.gmp_n, tmp);
		mpz_sqrt(tmp, b2);
		mpz_sub(tmp, a, tmp);
		gmp2mp(tmp, &tmpz);

		if (mpz_probab_prime_p(tmp, NUM_WITNESSES))
			tmpz.type = PRP;			
		else
			tmpz.type = PRP;
		add_to_factor_list(fobj, &tmpz);

		mpz_tdiv_q(fobj->div_obj.gmp_n, fobj->div_obj.gmp_n, tmp);
		zFree(&tmpz);
	}

	mpz_clear(tmp);
	mpz_clear(a);
	mpz_clear(b2);
	mpz_clear(maxa);

	gmp2mp(fobj->div_obj.gmp_n, n);

	return;

}

