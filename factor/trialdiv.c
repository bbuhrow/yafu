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

#include "factor.h"
#include "soe.h"
#include "gmp_xface.h"

static uint64_t* td_primes;
static uint64_t td_nump;
static uint64_t td_maxp;
static int initialized = 0;

void zTrial(fact_obj_t *fobj)
{
	//trial divide n using primes below limit. optionally, print factors found.
	//input expected in the gmp_n field of div_obj.
	uint32_t r,k=0;
	uint32_t limit = fobj->div_obj.limit;
	int print = fobj->div_obj.print;
	FILE *flog;
	uint64_t q;
	mpz_t tmp;
	mpz_init(tmp);

    if (fobj->LOGFLAG)
    {
        flog = fopen(fobj->flogname, "a");
        if (flog == NULL)
        {
            printf("fopen error: %s\n", strerror(errno));
            printf("could not open %s for writing\n", fobj->flogname);
            return;
        }
    }

    if (!initialized)
    {
        soe_staticdata_t* sdata = soe_init(0, 1, 32768);
        td_primes = soe_wrapper(sdata, 0, 1000000, 0, &td_nump, 0, 0);
        td_maxp = td_primes[td_nump - 1];
        soe_finalize(sdata);
        initialized = 1;
    }

	if (td_maxp < limit)
	{
        soe_staticdata_t* sdata = soe_init(0, 1, 32768);
		free(td_primes);
        td_primes = soe_wrapper(sdata, 0, limit, 0, &td_nump, 0, 0);
        td_maxp = td_primes[td_nump -1];
        soe_finalize(sdata);
	}

	while ((mpz_cmp_ui(fobj->div_obj.gmp_n, 1) > 0) && 
		(td_primes[k] < limit) &&
		(k < (uint32_t)td_nump))
	{
		q = (uint64_t)td_primes[k];
		r = mpz_tdiv_ui(fobj->div_obj.gmp_n, q);
		
		if (r != 0)
			k++;
		else
		{			
			mpz_tdiv_q_ui(fobj->div_obj.gmp_n, fobj->div_obj.gmp_n, q);
			mpz_set_64(tmp, q);

			add_to_factor_list(fobj->factors, tmp, fobj->VFLAG, fobj->NUM_WITNESSES);

#if BITS_PER_DIGIT == 64
			logprint(flog,"div: found prime factor = %" PRIu64 "\n",q);
#else
			logprint(flog,"div: found prime factor = %u\n",q);
#endif

			if (print && (fobj->VFLAG > 0))
#if BITS_PER_DIGIT == 64
				printf("div: found prime factor = %" PRIu64 "\n",q);
#else
				printf("div: found prime factor = %u\n",q);
#endif
		}
	}

    if (fobj->LOGFLAG)
    {
        fclose(flog);
    }

	mpz_clear(tmp);
}

void factor_perfect_power(fact_obj_t *fobj, mpz_t b)
{
	// check if (b^1/i)^i == b for i = 2 to bitlen(b)
	uint32_t bits = mpz_sizeinbase(b,2);
	uint32_t i;
	FILE *flog;
	mpz_t base, ans;

	mpz_init(base);
	mpz_init(ans);

    if (fobj->LOGFLAG)
    {
        //open the log file
        flog = fopen(fobj->flogname, "a");
        if (flog == NULL)
        {
            printf("fopen error: %s\n", strerror(errno));
            printf("could not open %s for writing\n", fobj->flogname);
            return;
        }
    }

	for (i=2; i<bits; i++)
	{
		mpz_root(base, b, i);
		mpz_pow_ui(ans, base, i);
		if (mpz_cmp(ans, b) == 0)
		{
			// found a base.  				
			if (is_mpz_prp(base, fobj->NUM_WITNESSES))
			{
				uint32_t j;
                char* s = mpz_get_str(NULL, 10, base);
				//gmp_printf("\nAdding prime base %Zd to factor list...\n", base);
				for (j=0; j<i; j++)
				{
					add_to_factor_list(fobj->factors, base, fobj->VFLAG, fobj->NUM_WITNESSES);
					mpz_tdiv_q(b, b, base);
					logprint(flog,"prp%d = %s\n",
						gmp_base10(base), s);
				}
                free(s);
			}
			else
			{
				// if composite, factor it and then multiply
				// all factors by i (the exponent).
				fact_obj_t *fobj_refactor;
				uint32_t j;

				gmp_printf("\nFactoring composite base %Zd...\n", base);

				// load the new fobj with this number
				fobj_refactor = (fact_obj_t *)malloc(sizeof(fact_obj_t));
				init_factobj(fobj_refactor, fobj->options);
				mpz_set(fobj_refactor->N, base);

				// recurse on factor
				//factor(fobj_refactor);
                // FIXME: need to find a better way to do this or 
                // a different place to put this function, since the
                // "factor()" function is not part of factor_common... it
                // is now in autofactor.  For now just do more trial division.
                zTrial(fobj_refactor);

				// add all factors found during the refactorization
				for (j=0; j< fobj_refactor->factors->num_factors; j++)
				{
					int k, c;
                    char* s = mpz_get_str(NULL, 10, fobj_refactor->factors->factors[j].factor);

					for (k=0; k < fobj_refactor->factors->factors[j].count; k++)
					{
						// add i copies of it, since this was a perfect power
						for (c = 0; c < i; c++)
						{
							add_to_factor_list(fobj->factors, 
                                fobj_refactor->factors->factors[j].factor,
                                fobj->VFLAG, fobj->NUM_WITNESSES);
							mpz_tdiv_q(b, b, fobj_refactor->factors->factors[j].factor);
							logprint(flog,"prp%d = %s\n",
								gmp_base10(fobj_refactor->factors->factors[j].factor), s);
						}
					}
                    free(s);
				}

				// free temps
				free_factobj(fobj_refactor);
				free(fobj_refactor);
			}
			break;
		}
	}

	mpz_clear(base);
	mpz_clear(ans);

    if (fobj->LOGFLAG)
    {
        fclose(flog);
    }

	return;
}

#define setbit(a,b) (((a)[(b) >> 3]) |= (nmasks[(b) & 7])) 
#define getbit(a,b) (((a)[(b) >> 3]) & (nmasks[(b) & 7])) 

void zFermat(uint64_t limit, uint32_t mult, fact_obj_t *fobj)
{
	// Fermat's factorization method with a sieve-based improvement
	// provided by 'neonsignal'
	mpz_t a, b2, tmp, multN, a2;
	int i;
    int sqchecks = 0;
	int numChars;
	uint64_t reportIt, reportInc;
	uint64_t count;
	uint64_t i64;
	FILE *flog = NULL;
	uint32_t M = 2 * 2 * 2 * 2 * 3 * 3 * 5 * 5 * 7 * 7; //176400u
	uint32_t M1 = 11 * 17 * 23 * 31; //133331u
	uint32_t M2 = 13 * 19 * 29 * 37; //265031u
	uint8_t *sqr, *sqr1, *sqr2, *mod, *mod1, *mod2;
	uint16_t *skip;
	uint32_t m, mmn, s, d;
	uint8_t masks[8] = {0xfe, 0xfd, 0xfb, 0xf7, 0xef, 0xdf, 0xbf, 0x7f};
	uint8_t nmasks[8];
	uint32_t iM = 0, iM1 = 0, iM2 = 0;

	if (mpz_even_p(fobj->div_obj.gmp_n))
	{
		mpz_init(tmp);
		mpz_set_ui(tmp, 2);
		mpz_tdiv_q_2exp(fobj->div_obj.gmp_n, fobj->div_obj.gmp_n, 1);
		add_to_factor_list(fobj->factors, tmp, fobj->VFLAG, fobj->NUM_WITNESSES);
		mpz_clear(tmp);
		return;
	}

	if (mpz_perfect_square_p(fobj->div_obj.gmp_n))
	{
		//open the log file
        if (fobj->LOGFLAG)
        {
            flog = fopen(fobj->flogname, "a");
            if (flog == NULL)
            {
                printf("fopen error: %s\n", strerror(errno));
                printf("could not open %s for writing\n", fobj->flogname);
                return;
            }
        }

		mpz_sqrt(fobj->div_obj.gmp_n, fobj->div_obj.gmp_n);

        char* s = mpz_get_str(NULL, 10, fobj->div_obj.gmp_n);
		if (is_mpz_prp(fobj->div_obj.gmp_n, fobj->NUM_WITNESSES))
		{			
			logprint(flog, "Fermat method found perfect square factorization:\n");
			logprint(flog,"prp%d = %s\n",
				gmp_base10(fobj->div_obj.gmp_n), s);
			logprint(flog,"prp%d = %s\n",
				gmp_base10(fobj->div_obj.gmp_n), s);
		}
		else
		{
			logprint(flog, "Fermat method found perfect square factorization:\n");
			logprint(flog,"c%d = %s\n",
				gmp_base10(fobj->div_obj.gmp_n), s);
			logprint(flog,"c%d = %s\n",
				gmp_base10(fobj->div_obj.gmp_n), s);
		}
        free(s);

		add_to_factor_list(fobj->factors, fobj->div_obj.gmp_n, fobj->VFLAG, fobj->NUM_WITNESSES);
		add_to_factor_list(fobj->factors, fobj->div_obj.gmp_n, fobj->VFLAG, fobj->NUM_WITNESSES);
		mpz_set_ui(fobj->div_obj.gmp_n, 1);
        
        if (fobj->LOGFLAG)
        {
            fclose(flog);
        }
		return;
	}

	mpz_init(a);
	mpz_init(b2);
	mpz_init(tmp);
	mpz_init(multN);
	mpz_init(a2);

	// apply the user supplied multiplier
	mpz_mul_ui(multN, fobj->div_obj.gmp_n, mult);
	
	// compute ceil(sqrt(multN))
	mpz_sqrt(a, multN);

	// form b^2
	mpz_mul(b2, a, a);
	mpz_sub(b2, b2, multN);

	// test successive 'a' values using a sieve-based approach.
	// the idea is that not all 'a' values allow a^2 or b^2 to be square.  
	// we pre-compute allowable 'a' values modulo various smooth numbers and 
	// build tables to allow us to quickly iterate over 'a' values that are 
	// more likely to produce squares.
	// init sieve structures
	sqr = (uint8_t*)calloc((M / 8 + 1) , sizeof(uint8_t));
	sqr1 = (uint8_t*)calloc((M1 / 8 + 1) , sizeof(uint8_t));
	sqr2 = (uint8_t*)calloc((M2 / 8 + 1) , sizeof(uint8_t));
	mod = (uint8_t*)calloc((M / 8 + 1) , sizeof(uint8_t));
	mod1 = (uint8_t*)calloc((M1 / 8 + 1) , sizeof(uint8_t));
	mod2 = (uint8_t*)calloc((M2 / 8 + 1) , sizeof(uint8_t));
	skip = (uint16_t*)malloc(M * sizeof(uint16_t));

	// test it.  This will be good enough if |u*p-v*q| < 2 * N^(1/4), where
	// mult = u*v
	count = 0;
	if (mpz_perfect_square_p(b2))
		goto found;

	for (i=0; i<8; i++)
		nmasks[i] = ~masks[i];

	// marks locations where squares can occur mod M, M1, M2
	for (i64 = 0; i64 < M; ++i64)
		setbit(sqr, (i64*i64)%M);

	for (i64 = 0; i64 < M1; ++i64)
		setbit(sqr1, (i64*i64)%M1);

	for (i64 = 0; i64 < M2; ++i64)
		setbit(sqr2, (i64*i64)%M2);

	// for the modular sequence of b*b = a*a - n values 
	// (where b2_2 = b2_1 * 2a + 1), mark locations where
	// b^2 can be a square
	m = mpz_mod_ui(tmp, a, M);
	mmn = mpz_mod_ui(tmp, b2, M);
	for (i = 0; i < M; ++i)
	{
		if (getbit(sqr, mmn)) setbit(mod, i);
		mmn = (mmn+m+m+1)%M;
		m = (m+1)%M;
	}

	// we only consider locations where the modular sequence mod M can
	// be square, so compute the distance to the next square location
	// at each possible value of i mod M.
	s = 0;
	d = 0;
	for (i = 0; !getbit(mod,i); ++i)
		++s;
	for (i = M; i > 0;)
	{
		--i;
		++s;
		skip[i] = s;
		if (s > d) d = s;
		if (getbit(mod,i)) s = 0;
	}
	//printf("maxSkip = %u\n", d);

	// for the modular sequence of b*b = a*a - n values 
	// (where b2_2 = b2_1 * 2a + 1), mark locations where the
	// modular sequence can be a square mod M1.  These will
	// generally differ from the sequence mod M.
	m = mpz_mod_ui(tmp, a, M1);
	mmn = mpz_mod_ui(tmp, b2, M1);
	for (i = 0; i < M1; ++i)
	{
		if (getbit(sqr1, mmn)) setbit(mod1, i);
		mmn = (mmn+m+m+1)%M1;
		m = (m+1)%M1;
	}

	// for the modular sequence of b*b = a*a - n values 
	// (where b2_2 = b2_1 * 2a + 1), mark locations where the
	// modular sequence can be a square mod M2.  These will
	// generally differ from the sequence mod M or M1.
	m = mpz_mod_ui(tmp, a, M2);
	mmn = mpz_mod_ui(tmp, b2, M2);
	for (i = 0; i < M2; ++i)
	{
		if (getbit(sqr2, mmn)) setbit(mod2, i);
		mmn = (mmn+m+m+1)%M2;
		m = (m+1)%M2;
	}

	// loop, checking for perfect squares
	mpz_mul_2exp(a2, a, 1);
	count = 0;
	numChars = 0;
	reportIt = limit / 100;
	reportInc = reportIt;
    
	do
	{
		d = 0;
		i64 = 0;
		do
		{
			// skip to the next possible square residue of b*b mod M
			s = skip[iM];

			// remember how far we skipped
			d += s;

			// update the other residue indices
			if ((iM1 += s) >= M1) iM1 -= M1;
			if ((iM2 += s) >= M2) iM2 -= M2;
			if ((iM += s) >= M) iM -= M;

			// some multpliers can lead to infinite loops.  bail out if so.
			if (++i64 > M) goto done;

			// continue if either of the other residues indicates non-square.
		} while (!getbit(mod1,iM1) || !getbit(mod2,iM2));

		// form b^2 by incrementing by many factors of 2*a+1
		mpz_add_ui(tmp, a2, d);
		mpz_mul_ui(tmp, tmp, d);
		mpz_add(b2, b2, tmp);

		// accumulate so that we can reset d 
		// (and thus keep it single precision)
		mpz_add_ui(a2, a2, d*2);
        
		count += d;
		if (count > limit)
			break;

		//progress report
		if ((count > reportIt) && (fobj->VFLAG > 1))
		{
			for (i=0; i< numChars; i++)
				printf("\b");
			numChars = printf("%" PRIu64 "%%",(uint64_t)((double)count / (double)limit * 100));
			fflush(stdout);
			reportIt += reportInc;
		}
        sqchecks++;
        
	} while (!mpz_perfect_square_p(b2));

    if (fobj->VFLAG > 1)
        printf("fmt: performed %d perfect square checks\n", sqchecks);

found:

	// 'count' is how far we had to scan 'a' to find a square b
	mpz_add_ui(a, a, count);

	if ((mpz_size(b2) > 0) && mpz_perfect_square_p(b2))
	{
		mpz_sqrt(tmp, b2); 		
		mpz_add(tmp, a, tmp);
		mpz_gcd(tmp, fobj->div_obj.gmp_n, tmp);

        if (fobj->LOGFLAG)
        {
            flog = fopen(fobj->flogname, "a");
            if (flog == NULL)
            {
                printf("fopen error: %s\n", strerror(errno));
                printf("could not open %s for writing\n", fobj->flogname);
                goto done;
            }
            logprint(flog, "Fermat method found factors:\n");
        }

		add_to_factor_list(fobj->factors, tmp, fobj->VFLAG, fobj->NUM_WITNESSES);
        char* s = mpz_get_str(NULL, 10, tmp);
		if (is_mpz_prp(tmp, fobj->NUM_WITNESSES))
		{			
			logprint(flog,"prp%d = %s\n",
				gmp_base10(tmp), s);
		}
		else
		{
			logprint(flog,"c%d = %s\n",
				gmp_base10(tmp), s);
		}
        free(s);

		mpz_tdiv_q(fobj->div_obj.gmp_n, fobj->div_obj.gmp_n, tmp);
		mpz_sqrt(tmp, b2);
		mpz_sub(tmp, a, tmp);
		mpz_gcd(tmp, fobj->div_obj.gmp_n, tmp);

		add_to_factor_list(fobj->factors, tmp, fobj->VFLAG, fobj->NUM_WITNESSES);
        s = mpz_get_str(NULL, 10, tmp);

		if (is_mpz_prp(tmp, fobj->NUM_WITNESSES))
		{			
			logprint(flog,"prp%d = %s\n",
				gmp_base10(tmp), s);
		}
		else
		{
			logprint(flog,"c%d = %s\n",
				gmp_base10(tmp), s);
		}
        free(s);

		mpz_tdiv_q(fobj->div_obj.gmp_n, fobj->div_obj.gmp_n, tmp);
	}

done:
	mpz_clear(tmp);
	mpz_clear(a);
	mpz_clear(b2);
	mpz_clear(multN);
	mpz_clear(a2);
	free(sqr);
	free(sqr1);
	free(sqr2);
	free(mod);
	free(mod1);
	free(mod2);
	free(skip);
    if (fobj->LOGFLAG && (flog != NULL))
    {
        fclose(flog);
    }
	return;

}

int sptestsqr(uint64_t n)
{
	uint64_t t;
	t = n & 31;
	if (t == 0 || t == 1 || t == 4 ||
		t == 9 || t == 16 || t == 17 || t == 25)
	{
		t = (uint64_t)sqrt((int64_t)n);
		if (n == t * t)
			return 1;
	}
	return 0;
}

uint64_t spfermat(uint64_t n, uint64_t limit)
{
	//	  Fermat's factorization method (wikipedia psuedo-code)
	uint64_t a, b2, maxa, count;

	maxa = n + 1;
	maxa >>= 1;

	a = (uint64_t)ceil(sqrt((int64_t)n));
	b2 = a * a - n;

	count = 0;
	while(!sptestsqr(b2))
	{
		// todo: special case this...
		a++;

		if (a > maxa)
			break;	//give up

		//b2 = a*a - N = b2 + 2*a - 1
		b2 = a*a - n;

		count++;
		if (count > limit)
			break;
	}

	if (sptestsqr(b2))
		return a + (uint64_t)sqrt(b2);
	else
		return 1;

}
