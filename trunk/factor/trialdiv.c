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
	//trial divide n using primes below limit. optionally, print factors found.
	//input expected in the gmp_n field of div_obj.
	uint32 r,k=0;
	uint32 limit = fobj->div_obj.limit;
	int print = fobj->div_obj.print;
	fp_digit q;
	mpz_t tmp;
	mpz_init(tmp);

	if (P_MAX < limit)
	{
		free(PRIMES);
		PRIMES = soe_wrapper(spSOEprimes, szSOEp, 0, limit, 0, &NUM_P);
		P_MIN = PRIMES[0];
		P_MAX = PRIMES[NUM_P-1];
	}

	while ((mpz_cmp_ui(fobj->div_obj.gmp_n, 1) > 0) && (PRIMES[k] < limit))
	{
		q = (fp_digit)PRIMES[k];
		r = mpz_tdiv_ui(fobj->div_obj.gmp_n, q);
		
		if (r != 0)
			k++;
		else
		{			
			mpz_tdiv_q_ui(fobj->div_obj.gmp_n, fobj->div_obj.gmp_n, q);
			mpz_set_64(tmp, q);

			add_to_factor_list(fobj, tmp);
			if (print && (VFLAG > 0))
#if BITS_PER_DIGIT == 64
				printf("div: found prime factor = %" PRIu64 "\n",q);
#else
				printf("div: found prime factor = %u\n",q);
#endif
		}
	}

	mpz_clear(tmp);
}

void zFermat(fp_digit limit, fact_obj_t *fobj)
{
	//	  Fermat's factorization method (wikipedia psuedo-code)
	//input expected in the gmp_n field of div_obj.
	mpz_t a, b2, tmp, maxa;
	int i;
	int numChars;
	fp_digit reportIt, reportInc;
	fp_digit count;
	FILE *flog;

	if (mpz_even_p(fobj->div_obj.gmp_n))
	{
		mpz_init(tmp);
		mpz_set_ui(tmp, 2);
		mpz_tdiv_q_2exp(fobj->div_obj.gmp_n, fobj->div_obj.gmp_n, 1);
		add_to_factor_list(fobj, tmp);
		mpz_clear(tmp);
		return;
	}

	if (mpz_perfect_square_p(fobj->div_obj.gmp_n))
	{
		//open the log file
		flog = fopen(fobj->flogname,"a");
		if (flog == NULL)
		{
			printf("could not open %s for writing\n",fobj->flogname);
			return;
		}

		mpz_sqrt(fobj->div_obj.gmp_n, fobj->div_obj.gmp_n);
		if (mpz_probab_prime_p(fobj->div_obj.gmp_n, NUM_WITNESSES))
		{			
			logprint(flog, "Fermat method found perfect square factorization:\n");
			logprint(flog,"prp%d = %s\n",
				mpz_sizeinbase(fobj->div_obj.gmp_n, 10),
				mpz_get_str(gstr1.s, 10, fobj->div_obj.gmp_n));
			logprint(flog,"prp%d = %s\n",
				mpz_sizeinbase(fobj->div_obj.gmp_n, 10),
				mpz_get_str(gstr1.s, 10, fobj->div_obj.gmp_n));
		}
		else
		{
			logprint(flog, "Fermat method found perfect square factorization:\n");
			logprint(flog,"c%d = %s\n",
				mpz_sizeinbase(fobj->div_obj.gmp_n, 10),
				mpz_get_str(gstr1.s, 10, fobj->div_obj.gmp_n));
			logprint(flog,"c%d = %s\n",
				mpz_sizeinbase(fobj->div_obj.gmp_n, 10),
				mpz_get_str(gstr1.s, 10, fobj->div_obj.gmp_n));
		}
		add_to_factor_list(fobj, fobj->div_obj.gmp_n);
		add_to_factor_list(fobj, fobj->div_obj.gmp_n);
		mpz_set_ui(fobj->div_obj.gmp_n, 1);		
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
		mpz_sqrt(tmp, b2); 
		mpz_add(tmp, a, tmp);

		flog = fopen(fobj->flogname,"a");
		if (flog == NULL)
		{
			printf("could not open %s for writing\n",fobj->flogname);
			return;
		}
		logprint(flog, "Fermat method found factors:\n");

		add_to_factor_list(fobj, tmp);
		if (mpz_probab_prime_p(tmp, NUM_WITNESSES))
		{			
			logprint(flog,"prp%d = %s\n",
				mpz_sizeinbase(tmp, 10),
				mpz_get_str(gstr1.s, 10, tmp));
		}
		else
		{
			logprint(flog,"c%d = %s\n",
				mpz_sizeinbase(tmp, 10),
				mpz_get_str(gstr1.s, 10, tmp));
		}

		mpz_tdiv_q(fobj->div_obj.gmp_n, fobj->div_obj.gmp_n, tmp);
		mpz_sqrt(tmp, b2);
		mpz_sub(tmp, a, tmp);

		add_to_factor_list(fobj, tmp);
		if (mpz_probab_prime_p(tmp, NUM_WITNESSES))
		{			
			logprint(flog,"prp%d = %s\n",
				mpz_sizeinbase(tmp, 10),
				mpz_get_str(gstr1.s, 10, tmp));
		}
		else
		{
			logprint(flog,"c%d = %s\n",
				mpz_sizeinbase(tmp, 10),
				mpz_get_str(gstr1.s, 10, tmp));
		}

		mpz_tdiv_q(fobj->div_obj.gmp_n, fobj->div_obj.gmp_n, tmp);
	}

	mpz_clear(tmp);
	mpz_clear(a);
	mpz_clear(b2);
	mpz_clear(maxa);

	return;

}

