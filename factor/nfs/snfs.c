/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

       				   --bbuhrow@gmail.com 12/6/2012
----------------------------------------------------------------------*/

#include <stdio.h>
#include "nfs_impl.h"
#include "ytools.h"
#include "gmp_xface.h"
#include <math.h>

#define POSITIVE 1
#define NEGATIVE 0

int tdiv_int(int x, int *factors, uint64_t *primes, uint64_t num_p);

void snfs_init(snfs_t* poly)
{	
	int i;
	memset(poly, 0, sizeof(snfs_t));
	poly->form_type = SNFS_NONE;
	poly->siever = 0;
	poly->poly = (mpz_polys_t*)malloc(sizeof(mpz_polys_t));
	if( !poly->poly )
	{
		printf("couldn't malloc!\n");
		exit(-1);
	}
	mpz_polys_init(poly->poly);
	poly->poly->side = RATIONAL_SPQ;
	poly->poly->rat.degree = 1;
	mpz_init(poly->n);
	mpz_init(poly->primitive);
	for (i=0; i<MAX_POLY_DEGREE + 1; i++)
		mpz_init(poly->c[i]);
	mpz_init(poly->base1);
	mpz_init(poly->base2);
}

void snfs_clear(snfs_t* poly)
{
	int i;
	mpz_clear(poly->n);
	mpz_clear(poly->primitive);
	mpz_polys_free(poly->poly);
	free(poly->poly);
	for (i=0; i<MAX_POLY_DEGREE + 1; i++)
		mpz_clear(poly->c[i]);
	mpz_clear(poly->base1);
	mpz_clear(poly->base2);
}

void snfs_copy_poly(snfs_t *src, snfs_t *dest)
{
	int i;
	mpz_set(dest->n, src->n);
	mpz_set(dest->primitive, src->primitive);

	mpz_set(dest->base1, src->base1);
	mpz_set(dest->base2, src->base2);

	dest->exp1 = src->exp1;
	dest->exp2 = src->exp2;
	dest->coeff1 = src->coeff1;
	dest->coeff2 = src->coeff2;
	dest->form_type = src->form_type;
	dest->siever = src->siever;

	for (i=0; i<MAX_POLY_DEGREE + 1; i++)
		mpz_set(dest->c[i], src->c[i]);

	dest->poly->rat.degree = src->poly->rat.degree;
	for(i = 0; i <= src->poly->rat.degree; i++)
		mpz_set(dest->poly->rat.coeff[i], src->poly->rat.coeff[i]);
	dest->poly->alg.degree = src->poly->alg.degree;
	for(i = 0; i <= src->poly->alg.degree; i++)
		mpz_set(dest->poly->alg.coeff[i], src->poly->alg.coeff[i]);
	dest->poly->skew = src->poly->skew;
	dest->poly->murphy = src->poly->murphy;
	dest->poly->size = src->poly->size;
	dest->poly->rroots = src->poly->rroots;
	mpz_set(dest->poly->m, src->poly->m);
	dest->poly->side = src->poly->side;

	dest->difficulty = src->difficulty;
	dest->sdifficulty = src->sdifficulty;
	dest->anorm = src->anorm;
	dest->rnorm = src->rnorm;
}

void check_poly(snfs_t *poly, int VFLAG)
{
	// make sure the generated poly is correct.  For instance if
	// the coefficents overflow an int the poly will be invalid, which
	// is fine since we wouldn't want to use coefficients that large anyway.
	// check that each polynomial mod n at m is zero
	mpz_t t;
	int i;
	mpz_init(t);

    if (VFLAG > 0)
        printf("nfs: checking degree %d poly\n", poly->poly->alg.degree);

	poly->valid = 1;
	mpz_set_ui(t, 0);
	for (i = poly->poly->alg.degree; i >= 0; i--)
	{
		mpz_mul(t, t, poly->poly->m);
		mpz_add(t, t, poly->c[i]);
		mpz_mod(t, t, poly->n);
	}

	// set mpz_poly_t alg appropriately
	for (i = poly->poly->alg.degree; i >= 0; i--)
		mpz_set(poly->poly->alg.coeff[i], poly->c[i]);

	mpz_mod(t, t, poly->n);
	if (mpz_cmp_ui(t,0) != 0)
	{
		poly->valid = 0;
		if (VFLAG > 0) 
		{
			gmp_fprintf(stdout, "Error: M=%Zd is not a root of f(x) % N\n"
				"n = %Zd\n", poly->poly->m, poly->n);
			fprintf (stderr, "f(x) = ");
			for (i = poly->poly->alg.degree; i >= 0; i--)
				gmp_fprintf (stdout, "%c %Zd*x^%d ",
				(mpz_cmp_ui(poly->c[i],0) < 0) ? '-' : '+', poly->poly->alg.coeff[i], i);
			gmp_fprintf (stdout, "\n""Remainder is %Zd\n\n", t);
		}
	}
	
	mpz_mul(t, poly->poly->rat.coeff[1], poly->poly->m);
	mpz_add(t, t, poly->poly->rat.coeff[0]);
	mpz_mod(t, t, poly->n);
	if (mpz_cmp_ui(t,0) != 0)
	{
		poly->valid = 0;
		if (VFLAG > 0)
		gmp_fprintf (stdout, "n = %Zd\n" "Error: M=%Zd is not a root of g(x) % N\n" "Remainder is %Zd\n\n", 
			poly->n, poly->poly->m, t);
	}
	
	//if (!isnormal(poly->anorm) || !isnormal(poly->rnorm)) // can they be negative?
	//{
	//	poly->valid = 0;
	//	if (fobj->VFLAG > 0)
	//		fprintf(stderr, "Error: invalid norms\n");
	//}
	//
	//if (!isnormal(poly->sdifficulty) || poly->sdifficulty <= 0 || !isnormal(poly->difficulty) || poly->difficulty <= 0)
	//{
	//	poly->valid = 0;
	//	if (fobj->VFLAG > 0)
	//		fprintf(stderr, "Error: invalid difficulties\n");
	//}

	mpz_clear(t);

	return;
}

void print_snfs(snfs_t *poly, FILE *out)
{
	// print the poly to stdout
	char c, side[80];
	int d = poly->difficulty; // round double to int

	if (poly->coeff2 < 0)
		c = '-';
	else
		c = '+';

	if (poly->poly->side == RATIONAL_SPQ)
		sprintf(side, "rational");
	else
		sprintf(side, "algebraic");

	gmp_fprintf(out, "n: %Zd\n", poly->n);

	if (poly->form_type == SNFS_H_CUNNINGHAM)
	{
		gmp_fprintf(out, "# %Zd^%d%c%Zd^%d, difficulty: %1.2f, anorm: %1.2e, rnorm: %1.2e\n", 
			poly->base1, poly->exp1, c, poly->base2, poly->exp1, poly->difficulty,
			poly->anorm, poly->rnorm);
	}
	else if (poly->form_type == SNFS_XYYXF)
	{
		gmp_fprintf(out, "# %Zd^%d+%Zd^%d, difficulty: %1.2f, anorm: %1.2e, rnorm: %1.2e\n", 
			poly->base1, poly->exp1, poly->base2, poly->exp2, poly->difficulty,
			poly->anorm, poly->rnorm);
	}
    else if (poly->form_type == SNFS_DIRECT)
    {
        gmp_fprintf(out, "# m=%Zd^%d, difficulty: %1.2f, anorm: %1.2e, rnorm: %1.2e\n",
            poly->base1, poly->exp1, poly->difficulty,
            poly->anorm, poly->rnorm);
    }
	else if (poly->form_type == SNFS_BRENT)
	{
		if (poly->coeff1 == 1)
			gmp_fprintf(out, "# %Zd^%d%c%d, difficulty: %1.2f, anorm: %1.2e, rnorm: %1.2e\n", 
				poly->base1, poly->exp1, c, abs(poly->coeff2), poly->difficulty,
				poly->anorm, poly->rnorm);
		else
			gmp_fprintf(out, "# %d*%Zd^%d%c%d, difficulty: %1.2f, anorm: %1.2e, rnorm: %1.2e\n", 
				abs(poly->coeff1), poly->base1, poly->exp1, c, abs(poly->coeff2), poly->difficulty,
				poly->anorm, poly->rnorm);
	}
    else
    {
        gmp_fprintf(out, "# difficulty: %1.2f, anorm: %1.2e, rnorm: %1.2e\n",
            poly->difficulty, poly->anorm, poly->rnorm);
    }

	if (poly->sdifficulty > 0)
		fprintf(out, "# scaled difficulty: %1.2f, suggest sieving %s side\n", poly->sdifficulty, side);

	// print recommended siever version, if known
	if (poly->siever > 0)
		fprintf(out, "# siever: %u\n", poly->siever);

	// msieve "analyze_one_poly" output, if known
    if (poly->poly->size > 0)
    {
        fprintf(out, "# size = %1.3e, alpha = %1.3f, combined = %1.3e, rroots = %d\n",
            poly->poly->size, poly->poly->alpha, poly->poly->murphy, poly->poly->rroots);
    }

	fprintf(out, "type: snfs\nsize: %d\n", d);

	print_poly(poly->poly, out);
}

void approx_norms(snfs_t *poly)
{
	// anorm ~= b^d * f(a/b) where f = the algebraic poly
	// rnorm ~= b * g(a/b) where g is the linear poly
	// a,b depend on the siever used.  according to the rsa768 paper, 
	// a,b were respectively about 3e9*sqrt(skew), 3e9/sqrt(skew).
	// here we use 1e6 instead of 3e9... it might not matter so much as
	// long as it is consistent between a/rnorm
	//int i;
	double a, b;
    int i;
    int found;
    double scale;
    double d = poly->difficulty;
    // pick an average size if for some reason we don't find one.
    int siever = 13;
	mpz_t res, tmp; // should be floats not ints perhaps

	// be sure poly->poly->alg is set properly
	if( !poly->valid )
		return; 

    // thanks to charybdis for suggesting that norms scale, 
    // at least with siever area.  Ideally Q would scale too; here
    // we just go with 1e6.  We approximate the norms during poly
    // generation, so we haven't picked a siever yet.  Siever selection
    // is based on difficulty, which we *have* computed, so run this
    // to pick up a siever estimate (scaling difficulty may later change 
    // it (rare) but this is an estimate anyway).
    for (i = 0; i < GGNFS_TABLE_ROWS - 1; i++)
    {
        if (d > ggnfs_table[i][0] && d <= ggnfs_table[i + 1][0])
        {
            scale = (double)(ggnfs_table[i + 1][0] - d) /
                (double)(ggnfs_table[i + 1][0] - ggnfs_table[i][0]);

            // pick closest entry
            if ((d - ggnfs_table[i][0]) < (ggnfs_table[i + 1][0] - d))
                siever = ggnfs_table[i][9];
            else
                siever = ggnfs_table[i + 1][9];

            found = 1;
        }
    }
    
    // https://www.mersenneforum.org/showpost.php?p=571762&postcount=3
	a = sqrt((double)(1ULL << (2 * siever - 1)) * 1000000.0 * poly->poly->skew);
    b = sqrt((double)(1ULL << (2 * siever - 1)) * 1000000.0 / poly->poly->skew);
    //a = sqrt(poly->poly->skew) * 1000000.;
	//b = 1000000. / (sqrt(poly->poly->skew));

    //printf("skew = %lf, a = %lf, b = %lf\n", poly->poly->skew, a, b);

	mpz_init(tmp);
	mpz_init(res);

	// use msieve functions to compute the norm size on each side.
	// these are used to skew the polynomial parameters, if the norms
	// are disparate enough
    //printf("alg degree: %d\n", poly->poly->alg.degree);
    //for (i = 0; i < 9; i++)
    //    gmp_printf("%Zd\n", poly->poly->alg.coeff[i]);

    //printf("rat degree: %d\n", poly->poly->rat.degree);
    //for (i = 0; i < 9; i++)
    //    gmp_printf("%Zd\n", poly->poly->rat.coeff[i]);

	eval_poly(res, (int64_t)a, (uint32_t)b, &poly->poly->alg);
	poly->anorm = mpz_get_d(res);

	eval_poly(res, (int64_t)a, (uint32_t)b, &poly->poly->rat);
	poly->rnorm = mpz_get_d(res);

    //printf("anorm = %le, bnorm = %le\n", poly->anorm, poly->rnorm);

	// extract from msieve (in a convoluted way) the Murphy score of the polynomial
	// along with several other values.  This will be used to rank the polynomials generated.
	analyze_one_poly_xface(poly);

	return;
}

void skew_snfs_params(fact_obj_t *fobj, nfs_job_t *job)
{
	// examine the difference between scaled difficulty and difficulty for snfs jobs
	// and skew the r/a parameters accordingly.
	// the input job struct should have already been filled in by get_ggnfs_params()	
	double oom_skew;
	double percent_skew;		
	if (job->snfs == NULL)
		return;

	// sdiff gets incremented by 1 for every 6 orders of magnitude difference
	// between the r/a norms. 
	oom_skew = job->snfs->sdifficulty - job->snfs->difficulty;

	if (job->snfs->poly->side == RATIONAL_SPQ)
	{
		// sieving on rational side means that side's norm is the bigger one
		if (oom_skew >= 4)
		{
			// for really big skew, increment the large prime bound as well
			job->lpbr++;
			job->mfbr += 2;
			// 5% for every 6 orders of magnitude difference between the norms	
			percent_skew = oom_skew * 0.05;
			job->alim -= percent_skew*job->alim/5;
			job->rlim += percent_skew*job->rlim;
		}

		if (oom_skew >= 5)
		{
			// for really big skew, increment the large prime bound as well
			job->lpbr++;
			job->mfbr += 2;
		}

		if (oom_skew >= 6)
		{
			// for really really big skew, use 3 large primes
			job->mfbr = job->lpbr*2.9;
			job->rlambda = 3.6;
		}

		// this *should* just change min_rels based on the new lpbr/a value
        if (oom_skew >= 4)
        {
            //get_ggnfs_params(fobj, job);
            nfs_set_min_rels(job);
        }

	}
	else
	{
		// sieving on algebraic side means that side's norm is the bigger one
		if (oom_skew >= 4)
		{
			// for really big skew, increment the large prime bound as well
			job->lpba++;
			job->mfba += 2;
			// 5% for every 6 orders of magnitude difference between the norms	
			percent_skew = oom_skew * 0.05;
			job->alim += percent_skew*job->alim;
			job->rlim -= percent_skew*job->rlim/5;
		}

		if (oom_skew >= 5)
		{
			// for really big skew, increment the large prime bound as well
			job->lpba++;
			job->mfba += 2;
		}

		if (oom_skew >= 6)
		{
			// for really really big skew, use 3 large primes
			job->mfba = job->lpba*2.9;
			job->alambda = 3.6;
		}

		// this *should* just change min_rels based on the new lpbr/a value
        if (oom_skew >= 4)
        {
            //get_ggnfs_params(fobj, job);
            nfs_set_min_rels(job);
        }
	}

	return;
}

void find_brent_form(fact_obj_t *fobj, snfs_t *form)
{
	int i,j,maxa,maxb,startb;
	mpz_t p, a, b, r, n, q;
	uint32_t inc = 1<<30;

	// cunningham numbers take the form a^n +- 1 with with a=2, 3, 5, 6, 7, 10, 11, 12
	// brent numbers take the form a^n +/- 1, where 13<=a<=99 and the product is less than 10^255.
	// generallized cullen/woodall numbers take the form a*b^a +/- 1
	// many oddperfect proof terms take the form a^n-1 for very large a
	// others of interest include k*2^n +/- 1, repunits, mersenne plus 2, etc.
	// All of the above forms (and more) can be represented by the generic form a*b^n +/- c
	// This routine will quickly discover any such generic form less than MAX_SNFS_BITS bits in size.

	maxa = 100;
	maxb = 100;
    startb = 31;

	mpz_init(p);
	mpz_init(a);
	mpz_init(b);
	mpz_init(r);
	mpz_init(n);
	mpz_init(q);

	mpz_set(n, fobj->nfs_obj.gmp_n);

	for (i=2; i<maxa; i++)
	{
		// ignore prime power bases... because e.g. 9^n == 3^(2*n) and thus
		// will have already been tested.
		if ((i==4) || (i==9) || (i==16) || (i==25) || (i==36) || (i==49) ||
			(i==64) || (i==81) || (i==8) || (i==27) || (i==32))
			continue;

		mpz_set_ui(b, i);
		mpz_pow_ui(p, b, startb-1); // p = i^31

		// limit the exponent so that the number is less than MAX_SNFS_BITS bits
		maxb = MAX_SNFS_BITS / log((double)i) + 1;

		if (fobj->VFLAG > 1)
			printf("nfs: checking a*b^x +/- c for 32 <= x <= %d\n", maxb);

		for (j=startb; j<maxb; j++)
		{
			// find any n = a*i^j + c, where a is arbitrary and c < 2^32
			mpz_mul(p, p, b);		// p = i^j
			mpz_add_ui(r, n, inc);	// r = n + 2^30
			mpz_mod(r, r, p);		// r = (n + 2^30) % i^j

			// now, if r is a single limb, then the input has a small coefficient
			if (mpz_sizeinbase(r,2) <= 32)
			{
				// get the second coefficient first
				int sign, c1, c2;
			
				c2 = mpz_get_ui(r);
				if (c2 > inc)
				{
					c2 = c2 - inc;
					sign = POSITIVE;
				}
				else
				{
					c2 = inc - c2;
					sign = NEGATIVE;
				}

				// now get any leading coefficient
				if (sign == POSITIVE)
					mpz_sub_ui(r, n, c2);	// r = n - c2
				else
					mpz_add_ui(r, n, c2);	// r = n + c2

				mpz_tdiv_qr(a, r, r, p);
				if (mpz_cmp_ui(r,0) != 0)
					continue;		// didn't divide, something's wrong

				if (mpz_sizeinbase(a,2) >= 32)
					continue;		// leading coefficient too big

				c1 = mpz_get_ui(a);


				// if the base divides the leading coefficient then we've just detected a 
				// degenerate form
				if (c1 % i == 0)
					continue;
				
				if (c1 > 1)
				{
					if (sign == POSITIVE)
					{
						if (fobj->VFLAG >= 0) printf("nfs: input divides %d*%d^%d + %d\n", c1, i, j, c2);
						logprint_oc(fobj->flogname, "a", "nfs: input divides %d*%d^%d + %d\n", c1, i, j, c2);
					}
					else
					{
						if (fobj->VFLAG >= 0) printf("nfs: input divides %d*%d^%d - %d\n", c1, i, j, c2);
						logprint_oc(fobj->flogname, "a", "nfs: input divides %d*%d^%d - %d\n", c1, i, j, c2);
					}
				}
				else
				{
					if (sign == POSITIVE)
					{
						if (fobj->VFLAG >= 0) printf("nfs: input divides %d^%d + %d\n", i, j, c2);
						logprint_oc(fobj->flogname, "a", "nfs: input divides %d^%d + %d\n", i, j, c2);
					}
					else
					{
						if (fobj->VFLAG >= 0) printf("nfs: input divides %d^%d - %d\n", i, j, c2);
						logprint_oc(fobj->flogname, "a", "nfs: input divides %d^%d - %d\n", i, j, c2);
					}
				}

				//printf("c1 = %d, i = %d, j = %d, c2 = %d\n", c1, i, j, c2);

				if ((c1 == 1) && (c2 == 1) && (sign == NEGATIVE) && (j % 2) == 0)
				{
					// with no leading coefficients, a -1 constant term, and
					// even exponent, the exponent can be divided by 2.
					j /= 2;
				}

				form->form_type = SNFS_BRENT;
				form->coeff1 = c1;
				mpz_set_ui(form->base1, i);
				mpz_set_ui(form->base2, 1);
				form->exp1 = j;
				form->coeff2 = sign ? c2 : -c2;
				if (mpz_cmp_ui(fobj->nfs_obj.snfs_cofactor, 1) > 0)
				{
                    if (fobj->LOGFLAG)
                    {
                        FILE* f = fopen(fobj->flogname, "a");
                        if (f != NULL)
                        {
                            logprint(f, "nfs: using supplied cofactor: ");
                            gmp_fprintf(f, "%Zd\n", fobj->nfs_obj.snfs_cofactor);
                            fclose(f);
                        }
                    }
					mpz_set(form->n, fobj->nfs_obj.snfs_cofactor);
				}
				else
					mpz_set(form->n, n);
				goto done;
			}

			// now find any *divisors* of n = i^j + c, where c < 2^32
			mpz_mod(r, p, n);
			mpz_sub(q, n, r);

			// now, if r is a single limb, then the input has a small coefficient
			if ((mpz_sizeinbase(r,2) <= 32) || (mpz_sizeinbase(q,2) <= 32))
			{
				// get the second coefficient first
				int sign, c1, c2;
				
				// if q is smaller, use it
				if (mpz_cmp(q, r) < 0)
				{
					//printf("n-(i^j) mod n is small: %d\n", mpz_get_ui(q));
					mpz_set(r, q);
					sign = POSITIVE;
				}
				else
					sign = NEGATIVE;

				c2 = mpz_get_ui(r);

				// this method doesn't detect a leading coefficient...
				c1 = 1;

				// if the base divides the leading coefficient then we've just detected a 
				// degenerate form
				if (c1 % i == 0)
					continue;
				
				if (c1 > 1)
				{
					if (sign == POSITIVE)
					{
						if (fobj->VFLAG >= 0) printf("nfs: input divides %d*%d^%d + %d\n", c1, i, j, c2);
						logprint_oc(fobj->flogname, "a", "nfs: input divides %d*%d^%d + %d\n", c1, i, j, c2);
					}
					else
					{
						if (fobj->VFLAG >= 0) printf("nfs: input divides %d*%d^%d - %d\n", c1, i, j, c2);
						logprint_oc(fobj->flogname, "a", "nfs: input divides %d*%d^%d - %d\n", c1, i, j, c2);
					}
				}
				else
				{
					if (sign == POSITIVE)
					{
						if (fobj->VFLAG >= 0) printf("nfs: input divides %d^%d + %d\n", i, j, c2);
						logprint_oc(fobj->flogname, "a", "nfs: input divides %d^%d + %d\n", i, j, c2);
					}
					else
					{
						if (fobj->VFLAG >= 0) printf("nfs: input divides %d^%d - %d\n", i, j, c2);
						logprint_oc(fobj->flogname, "a", "nfs: input divides %d^%d - %d\n", i, j, c2);
					}
				}

				if ((c1 == 1) && (c2 == 1) && (sign == NEGATIVE) && (j % 2) == 0)
					j /= 2;

				form->form_type = SNFS_BRENT;
				form->coeff1 = c1;
				mpz_set_ui(form->base1, i);
				mpz_set_ui(form->base2, 1);
				form->exp1 = j;
				form->coeff2 = sign ? c2 : -c2;
				if (mpz_cmp_ui(fobj->nfs_obj.snfs_cofactor, 1) > 0)
				{
                    if (fobj->LOGFLAG)
                    {
                        FILE* f = fopen(fobj->flogname, "a");
                        if (f != NULL)
                        {
                            logprint(f, "nfs: using supplied cofactor: ");
                            gmp_fprintf(f, "%Zd\n", fobj->nfs_obj.snfs_cofactor);
                            fclose(f);
                        }
                    }
					mpz_set(form->n, fobj->nfs_obj.snfs_cofactor);
				}
				else
					mpz_set(form->n, n);
				goto done;
			}

		}
	}


	// this checks x^n +- p, p small, up to a larger exponent value
	maxb = 1000;
	for (i = maxb; i>2; i--)
	{
		// now that we've reduced the exponent considerably, check for 
		// large bases by looking at the remaining possible exponents.
        if (fobj->VFLAG > 1)
        {
            printf("nfs: checking x^%d +/- 1\n", i);
        }
		
		// check -1 case:
		mpz_add_ui(a, n, 1);
		mpz_root(b, a, i);
		mpz_pow_ui(p, b, i);
		if (mpz_cmp(p, a) == 0)
		{
			char s[2048];

			if (fobj->VFLAG >= 0) gmp_printf("nfs: input divides %Zd^%d - 1\n", b, i);
			logprint_oc(fobj->flogname, "a", "nfs: input divides %s^%d - 1\n", mpz_get_str(s, 10, b), i);
			form->form_type = SNFS_BRENT;
			form->coeff1 = 1;
			mpz_set(form->base1, b);
			form->exp1 = i;
			mpz_set_ui(form->base2, 1);
			form->coeff2 = -1;
			// if the exponent is divisible by 2 in this case, then we can algebraically factor
			// as b^(2n) - 1 = (b^n + 1)(b^n - 1)
			if ((i & 0x1) == 0)
				form->exp1 /= 2;
			if (mpz_cmp_ui(fobj->nfs_obj.snfs_cofactor, 1) > 0)
			{
                if (fobj->LOGFLAG)
                {
                    FILE* f = fopen(fobj->flogname, "a");
                    if (f != NULL)
                    {
                        logprint(f, "nfs: using supplied cofactor: ");
                        gmp_fprintf(f, "%Zd\n", fobj->nfs_obj.snfs_cofactor);
                        fclose(f);
                    }
                }
				mpz_set(form->n, fobj->nfs_obj.snfs_cofactor);
			}
			else
				mpz_set(form->n, n);
			goto done;
		}

		// check +1 case:
		mpz_sub_ui(a, n, 1);
		mpz_root(b, a, i);
		mpz_pow_ui(p, b, i);
		if (mpz_cmp(p, a) == 0)
		{
			char s[2048];

			if (fobj->VFLAG >= 0) gmp_printf("nfs: input divides %Zd^%d + 1\n", b, i);
			logprint_oc(fobj->flogname, "a", "nfs: input divides %s^%d + 1\n", mpz_get_str(s, 10, b), i);
			form->form_type = SNFS_BRENT;
			form->coeff1 = 1;
			mpz_set(form->base1, b);
			form->exp1 = i;
			mpz_set_ui(form->base2, 1);
			form->coeff2 = 1;
			if (mpz_cmp_ui(fobj->nfs_obj.snfs_cofactor, 1) > 0)
			{
                if (fobj->LOGFLAG)
                {
                    FILE* f = fopen(fobj->flogname, "a");
                    if (f != NULL)
                    {
                        logprint(f, "nfs: using supplied cofactor: ");
                        gmp_fprintf(f, "%Zd\n", fobj->nfs_obj.snfs_cofactor);
                        fclose(f);
                    }
                }
				mpz_set(form->n, fobj->nfs_obj.snfs_cofactor);
			}
			else
				mpz_set(form->n, n);
			goto done;
		}

        // check other "+" cases:
        mpz_root(b, n, i);      // find b^i
        mpz_pow_ui(p, b, i);
        mpz_sub(p, n, p);       // and see if n - b^i is "small"

        if (mpz_sizeinbase(p, 2) < 31)
        {
			char* s = NULL; // [2048] ;
			char* s2 = NULL; // [2048] ;
			if (fobj->VFLAG >= 0)
				gmp_printf("nfs: input divides %Zd^%d + %Zd\n", b, i, p);

            logprint_oc(fobj->flogname, "a", "nfs: input divides %s^%d + %s\n",
                mpz_get_str(s, 10, b), i, mpz_get_str(s2, 10, p));

            form->form_type = SNFS_BRENT;
            form->coeff1 = 1;
            mpz_set(form->base1, b);
            form->exp1 = i;
            mpz_set_ui(form->base2, 1);
            form->coeff2 = mpz_get_ui(p);
            form->exp2 = 1;
			if (mpz_cmp_ui(fobj->nfs_obj.snfs_cofactor, 1) > 0)
			{
				if (fobj->LOGFLAG)
				{
					FILE* f = fopen(fobj->flogname, "a");
					if (f != NULL)
					{
						logprint(f, "nfs: using supplied cofactor: ");
						gmp_fprintf(f, "%Zd\n", fobj->nfs_obj.snfs_cofactor);
						fclose(f);
					}
				}
				mpz_set(form->n, fobj->nfs_obj.snfs_cofactor);
			}
			else
				mpz_set(form->n, n);

			free(s);
			free(s2);
            goto done;
        }

        // check other "-" cases:
        mpz_add_ui(b, b, 1);
        mpz_pow_ui(p, b, i);
        mpz_sub(p, p, n);       // and see if (b+1)^i - n is "small"

        if (mpz_sizeinbase(p, 2) < 31)
        {
			char* s = NULL; // [2048] ;
			char* s2 = NULL; // [2048] ;
            gmp_printf("nfs: input divides %Zd^%d - %Zd\n", b, i, p);

            logprint_oc(fobj->flogname, "a", "nfs: input divides %s^%d - %s\n", 
                mpz_get_str(s, 10, b), i, mpz_get_str(s2, 10, p));
            
			form->form_type = SNFS_BRENT;
			form->coeff1 = 1;
			mpz_set(form->base1, b);
			form->exp1 = i;
			mpz_set_ui(form->base2, 1);
			form->coeff2 = -(int)mpz_get_ui(p);
			form->exp2 = 1;
			if (mpz_cmp_ui(fobj->nfs_obj.snfs_cofactor, 1) > 0)
			{
				if (fobj->LOGFLAG)
				{
					FILE* f = fopen(fobj->flogname, "a");
					if (f != NULL)
					{
						logprint(f, "nfs: using supplied cofactor: ");
						gmp_fprintf(f, "%Zd\n", fobj->nfs_obj.snfs_cofactor);
						fclose(f);
					}
				}
				mpz_set(form->n, fobj->nfs_obj.snfs_cofactor);
			}
			else
				mpz_set(form->n, n);

			free(s);
			free(s2);
            goto done;
        }

	}


done:

	mpz_clear(p);
	mpz_clear(a);
	mpz_clear(b);
	mpz_clear(r);
	mpz_clear(q);
	mpz_clear(n);

	return;
}

void find_hcunn_form(fact_obj_t *fobj, snfs_t *form)
{
	int i,j,k,maxa,kmax;
	mpz_t pa, pb, a, b, r, g, n;

	// homogeneous cunninghams take the form a^n +/- b^n, where a,b <= 12 and gcd(a,b) == 1.
	// this routine finds inputs of the hcunninghams form, and considers
	// inputs up to MAX_SNFS_BITS bits in size.
	// once we have the form, we can create a polynomial for it for snfs processing.

	maxa = 51;

	mpz_init(pa);
	mpz_init(pb);
	mpz_init(a);
	mpz_init(g);
	mpz_init(b);
	mpz_init(r);
	mpz_init(n);

	mpz_set(n, fobj->nfs_obj.gmp_n);

	for (i=3; i<maxa; i++)
	{
		for (j=2; j<i; j++)
		{
			if (spGCD(i,j) != 1)
				continue;
			
			mpz_set_ui(a, i);
			mpz_pow_ui(pa, a, 19);
			mpz_set_ui(b, j);
			mpz_pow_ui(pb, b, 19);			

			// limit the exponent so that the number is less than MAX_SNFS_BITS bits
			kmax = MAX_SNFS_BITS / log((double)i) + 1;
			if (fobj->VFLAG > 1)
				printf("nfs: checking %d^x +/- %d^x for 20 <= x <= %d\n", i, j, kmax);

			for (k=20; k<kmax; k++)
			{
				mpz_mul(pa, pa, a);
				mpz_mul(pb, pb, b);

				mpz_add(g, pa, pb);
				mpz_mod(r, g, n);
				if (mpz_cmp_ui(r, 0) == 0)
				{
					if (fobj->VFLAG >= 0) printf("nfs: input divides %d^%d + %d^%d\n", i, k, j, k);
					logprint_oc(fobj->flogname, "a", "nfs: input divides %d^%d + %d^%d\n", i, k, j, k);
					form->form_type = SNFS_H_CUNNINGHAM;
					mpz_set_ui(form->base1, i);
					mpz_set_ui(form->base2, j);
					form->exp1 = k;
					form->coeff1 = 1;
					form->coeff2 = 1;
					if (mpz_cmp_ui(fobj->nfs_obj.snfs_cofactor, 1) > 0)
					{
                        if (fobj->LOGFLAG)
                        {
                            FILE* f = fopen(fobj->flogname, "a");
                            if (f != NULL)
                            {
                                logprint(f, "nfs: using supplied cofactor: ");
                                gmp_fprintf(f, "%Zd\n", fobj->nfs_obj.snfs_cofactor);
                                fclose(f);
                            }
                        }
						mpz_set(form->n, fobj->nfs_obj.snfs_cofactor);
					}
					else
						mpz_set(form->n, n);
					goto done;
				}

				mpz_sub(g, pa, pb);
				mpz_mod(r, g, n);
				if (mpz_cmp_ui(r, 0) == 0)
				{
					if (fobj->VFLAG >= 0) printf("nfs: input divides %d^%d - %d^%d\n", i, k, j, k);
					logprint_oc(fobj->flogname, "a", "nfs: input divides %d^%d - %d^%d\n", i, k, j, k);
					form->form_type = SNFS_H_CUNNINGHAM;
					mpz_set_ui(form->base1, i);
					mpz_set_ui(form->base2, j);
					form->exp1 = k;
					form->coeff1 = 1;
					form->coeff2 = -1;
					if (mpz_cmp_ui(fobj->nfs_obj.snfs_cofactor, 1) > 0)
					{
                        if (fobj->LOGFLAG)
                        {
                            FILE* f = fopen(fobj->flogname, "a");
                            if (f != NULL)
                            {
                                logprint(f, "nfs: using supplied cofactor: ");
                                gmp_fprintf(f, "%Zd\n", fobj->nfs_obj.snfs_cofactor);
                                fclose(f);
                            }
                        }
						mpz_set(form->n, fobj->nfs_obj.snfs_cofactor);
					}
					else
						mpz_set(form->n, n);
					goto done;
				}

			}
		}
	}

done:
	mpz_clear(pa);
	mpz_clear(pb);
	mpz_clear(g);
	mpz_clear(a);
	mpz_clear(b);
	mpz_clear(r);
	mpz_clear(n);
	return;
}

void find_xyyxf_form(fact_obj_t *fobj, snfs_t *form)
{
	int x,y,maxx;
	mpz_t xy, yx, r, g, n;

	// xyyxf numbers take the form x^y + y^x, where 1<y<x<151
	// this routine finds inputs of the xyyxf form, and considers
	// inputs up to MAX_SNFS_BITS bits in size.
	// once we have the form, we can create a polynomial for it for snfs processing.

	maxx = 201;

	mpz_init(xy);
	mpz_init(yx);
	mpz_init(g);
	mpz_init(r);
	mpz_init(n);

	mpz_set(n, fobj->nfs_obj.gmp_n);

	for (x=3; x<maxx; x++)
	{
		mpz_set_ui(xy, x);
		if (fobj->VFLAG > 1)
			printf("nfs: checking %d^y + y^%d\n", x, x);

		for (y=2; y<x; y++)
		{
			// with x fixed, we have computed x^y0.  successively multiply in x to this term
			// to avoid powering.
			// y^x we compute by powering.  we could probably create a 151x151 table of all powers of 
			// x and y without powering and then combine terms to check all possible x^y+y^x, but 
			// that is overkill for this routine that takes no noticable time anyway.
			mpz_mul_ui(xy, xy, x);
			mpz_set_ui(yx, y);
			mpz_pow_ui(yx, yx, x);

			mpz_add(g, xy, yx);
			mpz_mod(r, g, n);
			if (mpz_cmp_ui(r, 0) == 0)
			{
				if (fobj->VFLAG >= 0) printf("nfs: input divides %d^%d + %d^%d\n", x, y, y, x);
				logprint_oc(fobj->flogname, "a", "nfs: input divides %d^%d + %d^%d\n", x, y, y, x);
				form->form_type = SNFS_XYYXF;
				mpz_set_ui(form->base1, x);
				mpz_set_ui(form->base2, y);
				form->exp1 = y;
				form->exp2 = x;
				form->coeff1 = 1;
				form->coeff2 = 1;
				if (mpz_cmp_ui(fobj->nfs_obj.snfs_cofactor, 1) > 0)
				{
                    if (fobj->LOGFLAG)
                    {
                        FILE* f = fopen(fobj->flogname, "a");
                        if (f != NULL)
                        {
                            logprint(f, "nfs: using supplied cofactor: ");
                            gmp_fprintf(f, "%Zd\n", fobj->nfs_obj.snfs_cofactor);
                            fclose(f);
                        }
                    }
					mpz_set(form->n, fobj->nfs_obj.snfs_cofactor);
				}
				else
					mpz_set(form->n, n);
				goto done;
			}
		}
	}

done:
	mpz_clear(xy);
	mpz_clear(yx);
	mpz_clear(g);
	mpz_clear(r);
    mpz_clear(n);
	return;
}

void find_direct_form(fact_obj_t* fobj, snfs_t* form)
{
    int b, p, found = 0, i, deg;
    int c[10];
    mpz_t t, r, g, n, m;

    // find the following forms:
    // c1*m^n + c2*m^(n-1) + c3*m^(n-2) ... + cn
    // where c1...cn are all small, m=b^p, and nmax = 6.
    // then c1...cn are the coefficients of the algebraic polynomial,
    // Y1=-1, Y0=m is the rational polynomial, and n is the degree.
    // we search b=2..12 and p such that m<10^75 or so.

    mpz_init(t);
    mpz_init(g);
    mpz_init(r);
    mpz_init(n);
    mpz_init(m);

    mpz_set(n, fobj->nfs_obj.gmp_n);
    //mpz_set(n, fobj->N);

    for (b = 19; b > 1; b--)
    {
        mpz_set_ui(m, b);

        int maxp = 75 * (log(10) / log(b));

        if (fobj->VFLAG > 1)
            printf("nfs: checking m=%d^p forms, maxp = %d\n", b, maxp);

        for (p = maxp; p > 1; p--)
        {
            mpz_set_ui(m, b);
            mpz_pow_ui(m, m, p);

            if (fobj->VFLAG > 1)
                gmp_printf("checking m = %Zd = %d^%d\n", m, b, p);

            // if the remainder of n by m is small, that is a coefficient.
            // remove, reduce, and keep going until we either get a large
            // coefficient or n > 6
            mpz_set(t, n);
            found = 1;
            for (i = 0, deg = 0; i < 7; i++)
            {
                mpz_mod(r, t, m);
                if (mpz_cmp_ui(r, 10000) < 1)
                {
                    deg = i;
                    c[i] = mpz_get_ui(r);
                    mpz_sub_ui(t, t, c[i]);
                    mpz_tdiv_q(t, t, m);
                }
                else
                {
                    mpz_sub(r, m, r);
                    if (mpz_cmp_ui(r, 10000) < 1)
                    {
                        deg = i;
                        c[i] = -(int)mpz_get_ui(r);
                        mpz_add_ui(t, t, -c[i]);
                        mpz_tdiv_q(t, t, m);
                    }
                    else
                    {
                        found = 0;
                        break;
                    }
                }
            }

            // our test mpz 't' must now be 0 for the input to be represented by any
            // small-coefficient form we may have found.
            if (mpz_cmp_ui(t, 0) > 0)
            {
                found = 0;
            }

            // also, insist that we have at least degree 4.
            if (deg < 4)
            {
                found = 0;
            }
            
            if (found)
            {
                char nstr[128];
                char sign;
                strcpy(nstr, "");

                form->form_type = SNFS_DIRECT;
                mpz_set_ui(form->base1, b);
                form->exp1 = p;
                for (i = 8; i >= 0; i--)
                {
                    mpz_set_si(form->c[i], c[i]);
                }

                if (fobj->VFLAG >= 0)
                {
                    printf("nfs: input divides ");
                    for (i = 6; i > 0; i--)
                    {                       
                        if (abs(c[i]) > 0)
                        {
                            if (c[i] > 0) sign = '+'; else sign = '-';
                            printf("%c%d*(%d^%d)^%d ", sign, abs(c[i]), b, p, i);
                            sprintf(nstr, "%s%c%d*(%d^%d)^%d ", nstr, sign, abs(c[i]), b, p, i);
                        }
                    }
                    if (c[i] > 0) sign = '+'; else sign = '-';
                    printf("%c%d\n", sign, abs(c[i]));
                    sprintf(nstr, "%s%c%d", nstr, sign, abs(c[i]));
                }
                logprint_oc(fobj->flogname, "a", "nfs: input divides %s\n", nstr);

                if (mpz_cmp_ui(fobj->nfs_obj.snfs_cofactor, 1) > 0)
                {
                    if (fobj->LOGFLAG)
                    {
                        FILE* f = fopen(fobj->flogname, "a");
                        if (f != NULL)
                        {
                            logprint(f, "nfs: using supplied cofactor: ");
                            gmp_fprintf(f, "%Zd\n", fobj->nfs_obj.snfs_cofactor);
                            fclose(f);
                        }
                    }
                    mpz_set(form->n, fobj->nfs_obj.snfs_cofactor);
                }
                else
                {
                    mpz_set(form->n, n);
                }
                goto done;
            }
        }
    }

done:
    mpz_clear(t);
    mpz_clear(m);
    mpz_clear(g);
    mpz_clear(r);
    mpz_clear(n);

    return;
}

// see: http://home.earthlink.net/~elevensmooth/MathFAQ.html#PrimDistinct
void find_primitive_factor(fact_obj_t *fobj, snfs_t *poly, uint64_t* primes, uint64_t num_p, int VFLAG)
{
	// factor the exponent.  The algebraic reductions yafu knows how to handle are
	// for cunningham and homogenous cunningham inputs where the exponent is in
	// the exp1 field.
	int e = poly->exp1;
	int f[32];
	int nf, i, j, k, m, mult;
	// ranks of factors - we support up to 3 distinct odd factors of e
	int franks[4][32];		// and beans
	// and counts of the factors in each rank
	int cranks[4];			// doesn't that make you happy?   <-- bonus points if you get this...
	int nr, mrank;
	mpz_t n, term, t;

	mpz_set_ui(poly->primitive, 1);

	nf = tdiv_int(e, f, primes, num_p);	
	
	for (i=0; i<4; i++)
		cranks[i] = 0;

	// now arrange the factors into ranks of combinations of unique, distinct, and odd factors.
	// rank 0 is always 1
	franks[0][0] = 1;
	cranks[0] = 1;

	if (VFLAG > 2) printf("gen: finding primitive factor of exponent %d\n", e);
	// rank 1 is a list of the distinct odd factors.
	j=0;
	if (VFLAG > 2) printf("gen: rank 1 terms: ");
	for (i=0; i<nf; i++)
	{
		if (f[i] & 0x1)
		{
			// odd
			if (j==0 || f[i] != franks[1][j-1])
			{
				// distinct
				franks[1][j++] = f[i];
				if (VFLAG > 2) printf("%d ", f[i]);
			}
		}
	}
	if (VFLAG > 2) printf("\n");
	cranks[1] = j;
	nr = j + 1;

	if (j > 3)
	{
		printf("gen: too many distinct odd factors in exponent!\n");
		exit(1);
	}

	// ranks 2...nf build on the first rank combinatorially.
	// knuth, of course, has a lot to say on enumerating combinations:
	// http://www.cs.utsa.edu/~wagner/knuth/fasc3a.pdf
	// in which algorithm T might be sufficient since e shouldn't have too many
	// factors.
	// but since e shouldn't have too many factors and I don't feel like implementing
	// algorithm T from that reference right now, I will proceed to hardcode a bunch of
	// simple loops.
	// here is the second rank, if necessary
	if (cranks[1] == 2)
	{
		franks[2][0] = franks[1][0] * franks[1][1];
		cranks[2] = 1;
		if (VFLAG > 2) printf("gen: rank 2 term: %d\n", franks[2][0]);
	}
	else if (cranks[1] == 3)
	{
		// combinations of 2 primes
		m=0;
		if (VFLAG > 2) printf("gen: rank 2 terms: ");
		for (j=0; j<cranks[1]-1; j++)
		{
			for (k=j+1; k<cranks[1]; k++)
			{
				franks[2][m++] = franks[1][j] * franks[1][k];
				if (VFLAG > 2) printf("%d ", franks[2][m-1]);
			}
		}
		cranks[2] = m;
		if (VFLAG > 2) printf("\n");

		// combinations of 3 primes
		franks[3][0] = franks[1][0] * franks[1][1] * franks[1][2];
		cranks[3] = 1;
		if (VFLAG > 2) printf("gen: rank 3 term: %d\n", franks[3][0]);
	}

	// for exponents with repeated or even factors, find the multiplier
	mult = e;
	for (i=0; i<cranks[1]; i++)
		mult /= franks[1][i];
	if (VFLAG > 2) printf("gen: base exponent multiplier: %d\n", mult);

	// form the primitive factor, following the rank system of
	// http://home.earthlink.net/~elevensmooth/MathFAQ.html#PrimDistinct
	mpz_init(n);
	mpz_set_ui(n, 1);
	mpz_init(term);
	mpz_init(t);
	if ((nr & 0x1) == 1) 
		mrank = 0;
	else
		mrank = 1;
	for (i=nr-1; i >= 0; i--)
	{
		char c;
		if ((i & 0x1) == mrank)
		{
			// multiply by every other rank - do this before the division
			for (j=0; j<cranks[i]; j++)
			{			
				if (poly->form_type == SNFS_H_CUNNINGHAM)
				{
					mpz_set(term, poly->base1);
					mpz_pow_ui(term, term, franks[i][j] * mult);
					mpz_set(t, poly->base2);
					mpz_pow_ui(t, t, franks[i][j] * mult);
					if (poly->coeff2 < 0) {
						mpz_sub(term, term, t); c = '-';
					}
					else {
						mpz_add(term, term, t); c = '+';
					}
					if (VFLAG > 2) gmp_printf("gen: multiplying by %Zd^%d %c %Zd^%d = %Zd\n", 
						poly->base1, franks[i][j] * mult, c, poly->base2, franks[i][j] * mult, term);

				}
				else
				{
					mpz_set(term, poly->base1);
					mpz_pow_ui(term, term, franks[i][j] * mult);
					if (poly->coeff2 < 0) {
						mpz_sub_ui(term, term, 1); c = '-';
					}
					else {
						mpz_add_ui(term, term, 1); c = '+';
					}
					if (VFLAG > 2) gmp_printf("gen: multiplying by %Zd^%d %c 1 = %Zd\n", 
						poly->base1, franks[i][j] * mult, c, term);
				}
				mpz_mul(n, n, term);
			}
		}
	}
	for (i=nr-1; i >= 0; i--)
	{
		char c;
		if ((i & 0x1) == (!mrank))
		{
			// divide by every other rank
			for (j=0; j<cranks[i]; j++)
			{		
				if (poly->form_type == SNFS_H_CUNNINGHAM)
				{
					mpz_set(term, poly->base1);
					mpz_pow_ui(term, term, franks[i][j] * mult);
					mpz_set(t, poly->base2);
					mpz_pow_ui(t, t, franks[i][j] * mult);
					if (poly->coeff2 < 0) {
						mpz_sub(term, term, t); c = '-';
					}
					else {
						mpz_add(term, term, t); c = '+';
					}
					if (VFLAG > 2) gmp_printf("gen: dividing by %Zd^%d %c %Zd^%d = %Zd\n", 
						poly->base1, franks[i][j] * mult, c, poly->base2, franks[i][j] * mult, term);
				}
				else
				{
					mpz_set(term, poly->base1);
					mpz_pow_ui(term, term, franks[i][j] * mult);
					if (poly->coeff2 < 0) {
						mpz_sub_ui(term, term, 1); c = '-';
					}
					else {
						mpz_add_ui(term, term, 1); c = '+';
					}
					if (VFLAG > 2) gmp_printf("gen: dividing by %Zd^%d %c 1 = %Zd\n", 
						poly->base1, franks[i][j] * mult, c, term);
				}
				mpz_mod(t, n, term);
				if (mpz_cmp_ui(t, 0) != 0) printf("gen: error, term doesn't divide n!\n");
				mpz_tdiv_q(n, n, term);

				if (fobj->autofact_obj.autofact_active)
				{
					// does this term divide the input we are trying to autofactor?
					mpz_gcd(t, fobj->nfs_obj.gmp_n, term);
					if (mpz_cmp_ui(t, 1) > 0)
					{
						//gmp_printf("adding factor %Zd of autofactor input %Zd to factor list (gcd of term %Zd)\n", 
						//	t, fobj->nfs_obj.gmp_n, term);
						add_to_factor_list(fobj->factors, t, fobj->VFLAG, fobj->NUM_WITNESSES);
						mpz_tdiv_q(fobj->nfs_obj.gmp_n, fobj->nfs_obj.gmp_n, t);
					}
				}
			}
		}
	}

	mpz_tdiv_r(poly->primitive, poly->n, n);
	if (mpz_cmp_ui(poly->primitive, 0) != 0)
	{
		// we found a primitive factor that doesn't divide our
		// input cofactor.  This primitive factor is therefore 
		// a factor of some larger power of the target polynomial
		// that we were supplied a smaller cofactor of.
		// see if any part of the discovered primitive factor can 
		// be used
		mpz_gcd(t, poly->n, n);
		if (mpz_cmp_ui(t, 1) > 0)
		{
			// GCD of the primitive factor we discovered and our input discovered
			// a divisor.  Is it ever possible that this divisor is not useful?
			// i.e., that if we divide it out, what we are left with is a no-longer-snfsable
			// number that is more difficult than the original SNFS-able input?
			// we can't know unless we *don't* divide it out and continue with SNFS
			// poly generation, and then later compare everything.  
			if (VFLAG > 0) gmp_printf("gen: found primitive cofactor %Zd \n", n);
			if (VFLAG > 0) gmp_printf("gen: found factor of input with GCD: %Zd\n", t);
			// keep the discovered factor
			mpz_set(poly->primitive, t);
			// and don't reduce the input polynomial.
			// this will cause gen_brent_poly (if that is our caller)
			// to deal with the factor
		}
		else
		{
			mpz_set(poly->primitive, t);
		}
	}
	else
	{
		mpz_set(poly->primitive, n);

		if (mpz_cmp(n, poly->n) < 0)
		{
			if (VFLAG > 0) gmp_printf("gen: found primitive cofactor < input number:\ngen: %Zd\n", n);
			mpz_set(poly->n, n);
		}
	}

	if (fobj->autofact_obj.autofact_active)
	{
		// does this term divide the input we are trying to autofactor?
		mpz_mod(t, fobj->nfs_obj.gmp_n, poly->primitive);
		if (mpz_cmp_ui(t, 0) == 0)
		{
			//gmp_printf("adding primitive factor %Zd of autofactor input %Zd to factor list\n",
			//	poly->primitive, fobj->nfs_obj.gmp_n);
			add_to_factor_list(fobj->factors, poly->primitive, fobj->VFLAG, fobj->NUM_WITNESSES);
			mpz_tdiv_q(fobj->nfs_obj.gmp_n, fobj->nfs_obj.gmp_n, poly->primitive);
		}
	}

	mpz_clear(n);
	mpz_clear(term);
	mpz_clear(t);
	return;
}

// thanks to Alex Kruppa for his phi program, on which a lot
// of this routine is based.
snfs_t* gen_brent_poly(fact_obj_t *fobj, snfs_t *poly, int* npolys)
{
	int i, me;
	int e = poly->exp1;
	mpz_t b, b2;
	mpz_t n, m, t;
	double d, skew, k;
	int f[100];
	int numf = 0;
	snfs_t *polys = NULL;
	int npoly = 0;
	int apoly;	
	int algebraic = 0, halved = 0;

	mpz_init(n);
	mpz_init(m);
	mpz_init(t);

	mpz_init(b);
	mpz_init(b2);

	mpz_set(b, poly->base1);
	mpz_set(b2, poly->base2);

	// cunningham numbers take the form a^n +- 1 with with a=2, 3, 5, 6, 7, 10, 11, 12
	// brent numbers take the form a^n +/- 1, where 13<=a<=99 and the product is less than 10^255.
	// generallized cullen/woodall numbers take the form a*b^a +/- 1
	// many oddperfect proof terms take the form a^n-1 for very large a
	// others of interest include k*2^n +/- 1, repunits, mersenne plus 2, etc.
	// All of the above forms (and more) can be represented by the generic form a*b^n +/- c
	// This routine generates polynomials suitable for input into the gnfs lattice sievers from
	// kleinjung/franke and considers simple reductions by algebraic factors where possible.
	// Further, homogeneous cunninghams take the form a^n +/- b^n, where a,b <= 12 and gcd(a,b) == 1.
	// These can be rearranged to resemble the generic form by dividing through by b^n, so that
	// we have (a/b)^n +/- 1.  This requires straightforward modifications to the polynomial 
	// generation

	// first test for algebraic factors - we look for several specific but common forms:
	// 21*k, 15*k, 13*k, 11*k, 7*k, 5*k, 3*k
	// if any of these are available we immediately use them, because dividing out an algebraic factor
	// will always be lower difficulty then playing with exponents only, even if the degree
	// is sub-optimal.  The possibility of simple algebraic reduction occurs only when the a,c 
	// coefficients are 1.  More complex algebraic reductions like Aurifeuillian 
	// factorizations are not attempted here.
	if ((poly->form_type != SNFS_DIRECT) && (poly->coeff1 == 1) && (abs(poly->coeff2) == 1))
	{
		find_primitive_factor(fobj, poly, fobj->primes, fobj->num_p, fobj->VFLAG);
	}


    if (mpz_cmp(poly->primitive, poly->n) == 0)
    {
		// it's possible for this primitive factor to be prime, in which case we are done.
		if (is_mpz_prp(poly->primitive, 1))
		{

			if (fobj->LOGFLAG)
			{
				FILE* f = fopen(fobj->flogname, "a");

				if (f != NULL)
				{
					logprint(f, "nfs: snfs primitive factor is a probable prime\nnfs: ");
					gmp_fprintf(f, "%Zd\n", poly->primitive);
					fclose(f);
				}
			}

			if (fobj->VFLAG > 0)
			{
				gmp_printf("nfs: snfs primitive factor %Zd is a probable prime\n", poly->primitive);
			}

			char c[4];

			add_to_factor_list(fobj->factors, poly->primitive, fobj->VFLAG, fobj->NUM_WITNESSES);
			strncpy(c, "prp", 4);

			if (fobj->LOGFLAG)
			{
				FILE *logfile = fopen(fobj->flogname, "a");
				if (logfile == NULL)
				{
					printf("fopen error: %s\n", strerror(errno));
					printf("could not open yafu logfile for appending\n");
				}
				else
				{
					char* s = mpz_get_str(NULL, 10, poly->primitive);
					logprint(logfile, "%s%d = %s\n", c,
						gmp_base10(poly->primitive), s);
					fclose(logfile);
					free(s);
				}
			}

			mpz_tdiv_q(fobj->nfs_obj.gmp_n, fobj->nfs_obj.gmp_n, poly->primitive);
			gmp_printf("nfs: residue is %Zd\n", fobj->nfs_obj.gmp_n);
			*npolys = 0;
			return NULL;

		}
		else
		{
			// it's not prime.  If it has a non-trivial residue with the input
			// number, add the residue to the factor list.
			mpz_tdiv_r(t, fobj->nfs_obj.gmp_n, poly->primitive);
			if ((mpz_cmp_ui(t, 0) == 0) && (mpz_cmp_ui(poly->primitive, 1) > 0) &&
				(mpz_cmp(fobj->nfs_obj.gmp_n, poly->primitive) < 0))
			{
				mpz_tdiv_q(t, fobj->nfs_obj.gmp_n, poly->primitive);
				//gmp_printf("fac: adding residue %Zd to the factor list\n", t);
				add_to_factor_list(fobj->factors, t, fobj->VFLAG, fobj->NUM_WITNESSES);
				mpz_set(fobj->nfs_obj.gmp_n, poly->primitive);
			}

			if (fobj->LOGFLAG)
			{
				FILE* f = fopen(fobj->flogname, "a");

				if (f != NULL)
				{
					logprint(f, "nfs: commencing snfs on c%d primitive factor: ",
						gmp_base10(poly->primitive));
					gmp_fprintf(f, "%Zd\n", poly->primitive);
					fclose(f);
				}
			}
		}
    }
	else if (mpz_cmp_ui(poly->primitive, 0) > 0)
	{
		// gmp_printf("found factor %Zd of input %Zd\n", poly->primitive, fobj->nfs_obj.gmp_n);
		// we found a factor of the input that isn't a primitive factor
		add_to_factor_list(fobj->factors, poly->primitive, fobj->VFLAG, fobj->NUM_WITNESSES);

		// reduce the input and keep going
		mpz_tdiv_q(fobj->nfs_obj.gmp_n, fobj->nfs_obj.gmp_n, poly->primitive);
		//mpz_set(poly->n, fobj->nfs_obj.gmp_n);
		if (fobj->VFLAG > 0)
			gmp_printf("nfs: continuing with residue %Zd\n", fobj->nfs_obj.gmp_n);

		*npolys = 0;
		return NULL;
	}
    else
    {
		if (fobj->LOGFLAG)
		{
			FILE* f = fopen(fobj->flogname, "a");

			if (f != NULL)
			{
				logprint(f, "nfs: commencing snfs on c%d: ",
					gmp_base10(poly->n));
				gmp_fprintf(f, "%Zd\n", poly->n);
				fclose(f);
			}
		}
    }


    if (poly->form_type == SNFS_DIRECT)
    {
        // the poly form directly gives a poly.
        int deg = 0;
        double d;

        mpz_set(m, poly->base1);
        mpz_pow_ui(m, m, poly->exp1);

        for (i = 6; i >= 0; i--)
        {
            if ((mpz_cmp_ui(poly->c[i],0) > 0) && (deg == 0))
            {
                deg = i;
            }
            mpz_set(poly->poly->alg.coeff[i], poly->c[i]);
        }
        d = mpz_get_d(m);
        d = log10(d) * (double)deg;

        // leading coefficient contributes to the difficulty
        d += log10(mpz_get_d(poly->c[deg]));

        // compute skew
        skew = pow(fabs(mpz_get_d(poly->c[0])) / mpz_get_d(poly->c[deg]), 1. / (double)deg);

        //printf("degree is %d\n", deg);
        //printf("skew is %lf\n", skew);
        //printf("difficulty is %lf\n", d);

        poly->difficulty = d;
        mpz_set(poly->poly->m, m);
        poly->poly->skew = skew;
        poly->poly->alg.degree = deg;

        mpz_set_si(poly->poly->rat.coeff[1], -1);
        mpz_set(poly->poly->rat.coeff[0], m);
 
        algebraic = 0;

        polys = (snfs_t*)malloc(sizeof(snfs_t));
        snfs_init(polys);
        npoly = 1;
        snfs_copy_poly(poly, polys);		// copy algebraic form

        if (fobj->VFLAG > 0)
        {
            printf("gen: ========================================================\n"
                "gen: considering the following polynomials:\n"
                "gen: ========================================================\n\n");
        }
        logprint_oc(fobj->nfs_obj.logfile, "a", "gen: considering the following polynomials:\n");

        check_poly(&polys[0], fobj->VFLAG);
        approx_norms(&polys[0]);

        if (polys[0].valid)
        {
            if (fobj->VFLAG > 0) print_snfs(&polys[0], stdout);
        }

    }
	else if ((poly->exp1 % 15 == 0) && (poly->coeff1 == 1) && (abs(poly->coeff2) == 1))
	{
		polys = (snfs_t *)malloc(sizeof(snfs_t));
		snfs_init(polys);
		npoly = 1;
		snfs_copy_poly(poly, polys);		// copy algebraic form

		// a^(15k) +/- 1 has an algebraic factor which is an 8th degree symmetric polynomial.
		// check for this case before a^(3k) or a^(5k) since the algbraic reduction is greater.
		// the 8th degree poly can be halved as in the 11*k and 13*k cases, resulting in 
		// the following quartic:
		polys->poly->alg.degree = 4;
		fobj->nfs_obj.pref_degree = 4;
		k = poly->exp1 / 15;		
		mpz_set_ui(polys->c[4], 1);
		mpz_set_si(polys->c[3], poly->coeff2);
		mpz_set_si(polys->c[2], -4);
		mpz_set_si(polys->c[1], -poly->coeff2 * 4);
		mpz_set_ui(polys->c[0], 1);
		mpz_set(m, poly->base1);
		polys->difficulty = log10(mpz_get_d(m)) * 8. * k;
		mpz_pow_ui(polys->poly->m, m, k);
		polys->poly->skew = 1.;

		algebraic = 1;
		halved = 1;
	}
	else if ((poly->exp1 % 21 == 0) && (poly->coeff1 == 1) && (abs(poly->coeff2) == 1))
	{
		polys = (snfs_t *)malloc(sizeof(snfs_t));
		snfs_init(polys);
		npoly = 1;
		snfs_copy_poly(poly, polys);		// copy algebraic form

		// a^(21k) +/- 1 has an algebraic factor which is a 12th degree symmetric polynomial.
		// check for this case before a^(3k) or a^(7k) since the algbraic reduction is greater.
		// the 12th degree poly can be halved as in the 11*k and 13*k cases, resulting in 
		// the following sextic:
		polys->poly->alg.degree = 6;
		fobj->nfs_obj.pref_degree = 6;
		k = poly->exp1 / 21;
		mpz_set_ui(polys->c[6], 1);
		mpz_set_si(polys->c[5], poly->coeff2);
		mpz_set_si(polys->c[4], -6);
		mpz_set_si(polys->c[3], -poly->coeff2 * 6);
		mpz_set_ui(polys->c[2], 8);
		mpz_set_si(polys->c[1], poly->coeff2 * 8);
		mpz_set_ui(polys->c[0], 1);
		mpz_set(m, poly->base1);
		polys->difficulty = log10(mpz_get_d(m)) * 12. * k;
		mpz_pow_ui(polys->poly->m, m, k);
		polys->poly->skew = 1.;
		
		algebraic = 1;
		halved = 1;
	}
    else if ((poly->exp1 % 3 == 0) && (poly->coeff1 == 1) && (abs(poly->coeff2) == 1))
    {
        // todo:
        // one case that looks to see if exponent is divisible by 3.
        // then handle conversion to a quartic or sextic depending on input difficulty.
        double d4;
        double d6;
        int degree;

        mpz_set(m, poly->base1);

        if (poly->exp1 % 6 == 0)
        {
            k = (poly->exp1 / 6);
            d4 = log10(mpz_get_d(m)) * 4. * k;
        }
        else
        {
            k = (poly->exp1 - 3) / 6;
            d4 = log10(mpz_get_d(m)) * (4. * k + 2);
        }

        if (poly->exp1 % 9 == 0)
        {
            k = (poly->exp1 / 9);
            d6 = log10(mpz_get_d(m)) * 6. * k;
        }
        else if (poly->exp1 % 9 == 3)
        {
            k = (poly->exp1 - 3) / 9;
            d6 = log10(mpz_get_d(m)) * (6. * k + 2);
        }
        else
        {
            k = (poly->exp1 + 3) / 9;
            d6 = log10(mpz_get_d(m)) * (6. * k);
        }

        poly->poly->rat.degree = 1;
        if (d4 < 160)
        {
            poly->poly->alg.degree = 4;
            fobj->nfs_obj.pref_degree = 4;
            degree = 4;
        }
        else
        {
            poly->poly->alg.degree = 6;
            fobj->nfs_obj.pref_degree = 6;
            degree = 6;
        }

        printf("nfs: degree 4 difficulty = %1.2f, degree 6 difficulty = %1.2f\n", d4, d6);
        printf("nfs: choosing degree %d\n", degree);


        if ((poly->exp1 % 6 == 0) && (poly->coeff1 == 1) && (abs(poly->coeff2) == 1) && (degree == 4))
        {
            polys = (snfs_t *)malloc(sizeof(snfs_t));
            snfs_init(polys);
            npoly = 1;
            snfs_copy_poly(poly, polys);		// copy algebraic form

            // a^(3k) +/- 1, k even, is divisible by (a^k +/- 1) giving a quadratic in a^k.
            // the quadratic can be converted into a quartic...
            // see: http://www.mersennewiki.org/index.php/SNFS_Polynomial_Selection
            // todo: look into making degree 6 based on difficulty
            k = poly->exp1 / 6;
            mpz_set_ui(polys->c[4], 1);
            mpz_set_si(polys->c[2], -poly->coeff2);
            mpz_set_ui(polys->c[0], 1);
            mpz_set(m, poly->base1);
            polys->difficulty = log10(mpz_get_d(m)) * 4. * k;
            mpz_pow_ui(polys->poly->m, m, k);
            polys->poly->skew = 1.;
        }
        else if ((poly->exp1 % 6 == 3) && (poly->coeff1 == 1) && (abs(poly->coeff2) == 1) && (degree == 4))
        {
            polys = (snfs_t *)malloc(sizeof(snfs_t));
            snfs_init(polys);
            npoly = 1;
            snfs_copy_poly(poly, polys);		// copy algebraic form

            // a^(3k) +/- 1, k odd, is divisible by (a^k +/- 1) giving a quadratic in a^k.
            // the quadratic can be converted into a quartic...
            // see: http://www.mersennewiki.org/index.php/SNFS_Polynomial_Selection
            // todo: look into making degree 6 based on difficulty
            k = (poly->exp1 - 3) / 6;
            mpz_mul(polys->c[4], poly->base1, poly->base1);
            mpz_set_si(polys->c[2], -poly->coeff2);
            mpz_mul(polys->c[2], polys->c[2], poly->base1);
            mpz_mul(polys->c[2], polys->c[2], poly->base2);
            mpz_mul(polys->c[0], poly->base2, poly->base2);

            mpz_set(m, poly->base1);
            polys->poly->skew = pow(mpz_get_d(poly->base1) / mpz_get_d(poly->base2), -0.5);
            polys->difficulty = log10(mpz_get_d(m)) * (4. * k + 2);
            mpz_pow_ui(polys->poly->m, m, k);
        }
        else if ((poly->exp1 % 9 == 0) && (poly->coeff1 == 1) && (abs(poly->coeff2) == 1) && (degree == 6))
        {            
            polys = (snfs_t *)malloc(sizeof(snfs_t));
            snfs_init(polys);
            npoly = 1;
            snfs_copy_poly(poly, polys);		// copy algebraic form

            // a^(3k) +/- 1, k odd, is divisible by (a^k +/- 1) giving a quadratic in a^k.
            // the quadratic can be converted into a sextic...
            // see: http://www.mersennewiki.org/index.php/SNFS_Polynomial_Selection            
            k = (poly->exp1 / 9);

            mpz_set_ui(polys->c[6], 1);
            mpz_set_si(polys->c[3], -poly->coeff2);
            mpz_set_ui(polys->c[0], 1);

            polys->poly->skew = 1.;
            polys->difficulty = log10(mpz_get_d(m)) * 6. * k;
            mpz_pow_ui(polys->poly->m, m, k);

        }
        else if ((poly->exp1 % 9 == 3) && (poly->coeff1 == 1) && (abs(poly->coeff2) == 1) && (degree == 6))
        {
            polys = (snfs_t *)malloc(sizeof(snfs_t));
            snfs_init(polys);
            npoly = 1;
            snfs_copy_poly(poly, polys);		// copy algebraic form

            // a^(3k) +/- 1, k odd, is divisible by (a^k +/- 1) giving a quadratic in a^k.
            // the quadratic can be converted into a sextic in q with k = 3q + 1.
            // (this brings in factors of base2 in hcunn forms).
            k = (poly->exp1 - 3) / 9;

            mpz_set(polys->c[6], poly->base1);
            mpz_mul(polys->c[6], polys->c[6], poly->base1);
            mpz_set_si(polys->c[3], -poly->coeff2);
            mpz_mul(polys->c[3], polys->c[3], poly->base1);
            mpz_mul(polys->c[3], polys->c[3], poly->base2);
            mpz_mul(polys->c[0], poly->base2, poly->base2);

            polys->poly->skew = pow(mpz_get_d(poly->base1) / mpz_get_d(poly->base2), -1. / 3.);
            polys->difficulty = log10(mpz_get_d(m)) * (6. * k + 2);
            mpz_pow_ui(polys->poly->m, m, k);
        }
        else if ((poly->exp1 % 9 == 6) && (poly->coeff1 == 1) && (abs(poly->coeff2) == 1) && (degree == 6))
        {
            polys = (snfs_t *)malloc(sizeof(snfs_t));
            snfs_init(polys);
            npoly = 1;
            snfs_copy_poly(poly, polys);		// copy algebraic form

            // a^(3k) +/- 1, k odd, is divisible by (a^k +/- 1) giving a quadratic in a^k.
            // the quadratic can be converted into a sextic in q with k = 3q - 1.
            // (this brings in factors of base2 in hcunn forms).
            k = (poly->exp1 + 3) / 9;

            mpz_mul(polys->c[6], poly->base2, poly->base2);
            mpz_set_si(polys->c[3], -poly->coeff2);
            mpz_mul(polys->c[3], polys->c[3], poly->base1);
            mpz_mul(polys->c[3], polys->c[3], poly->base2);
            mpz_mul(polys->c[0], poly->base1, poly->base1);

            //polys->poly->skew = pow(mpz_get_d(poly->base1) / mpz_get_d(poly->base2), 1./3.);
            polys->poly->skew = pow(mpz_get_d(poly->base1), 1. / 3.);
            polys->difficulty = log10(mpz_get_d(m)) * (6. * k);
            mpz_pow_ui(polys->poly->m, m, k);
        }

        algebraic = 1;

    }
    else if ((poly->exp1 % 35 == 0) && (poly->coeff1 == 1) && (abs(poly->coeff2) == 1))
    {
        // a^(35k) +/- 1 is divisible by (a^k +/- 1) giving 
        // either a quartic or sextic in a^k.  Choose which one based on difficulty.
        double d4;
        double d6;
        double ratio;
        int degree;

        mpz_set(m, poly->base1);
        k = poly->exp1 / 5;
        d4 = log10(mpz_get_d(m)) * 4. * k;

        k = poly->exp1 / 7;
        d6 = log10(mpz_get_d(m)) * 6. * k;

        // compare and choose one or the other based on difficulty.
        // if d4 is over ~160 then d6 may be the way to go unless d6
        // is much larger.  The ratio d6/d4 that is acceptable probably
        // grows with increasing d4, because degree 4 becomes more and
        // more undesirable with increasing size.  For now use this ad-hoc
        // table:
        ratio = d6 / d4;

        printf("d4 difficulty = %1.4f, d6 difficult = %1.4f, ratio = %1.4f\n", d4, d6, ratio);

        if (d4 > 220)
        {
            if (ratio < 1.2)
            {
                printf("d4 cutoff = 1.2: choosing degree 6\n");
                degree = 6;
            }
            else
            {
                printf("d4 cutoff = 1.2: choosing degree 4\n");
                degree = 4;
            }
        }
        else if (d4 > 200)
        {
            if (ratio < 1.15)
            {
                printf("d4 cutoff = 1.15: choosing degree 6\n");
                degree = 6;
            }
            else
            {
                printf("d4 cutoff = 1.15: choosing degree 4\n");
                degree = 4;
            }
        }
        else if (d4 > 180)
        {
            if (ratio < 1.1)
            {
                printf("d4 cutoff = 1.1: choosing degree 6\n");
                degree = 6;
            }
            else
            {
                printf("d4 cutoff = 1.1: choosing degree 4\n");
                degree = 4;
            }
        }
        else if (d4 > 160)
        {
            if (ratio < 1.05)
            {
                printf("d4 cutoff = 1.05: choosing degree 6\n");
                degree = 6;
            }
            else
            {
                printf("d4 cutoff = 1.05: choosing degree 4\n");
                degree = 4;
            }
        }
        else
        {
            degree = 4;
        }

        if (degree == 4)
        {
            polys = (snfs_t *)malloc(sizeof(snfs_t));
            snfs_init(polys);
            npoly = 1;
            snfs_copy_poly(poly, polys);		// copy algebraic form

            // a^(5k) +/- 1 is divisible by (a^k +/- 1) giving a quartic in a^k
            polys->poly->alg.degree = 4;
            fobj->nfs_obj.pref_degree = 4;
            k = poly->exp1 / 5;
            mpz_set_ui(polys->c[4], 1);
            mpz_set_si(polys->c[3], -poly->coeff2);
            mpz_set_ui(polys->c[2], 1);
            mpz_set_si(polys->c[1], -poly->coeff2);
            mpz_set_ui(polys->c[0], 1);
            mpz_set(m, poly->base1);
            polys->difficulty = log10(mpz_get_d(m)) * 4. * k;
            mpz_pow_ui(polys->poly->m, m, k);
            polys->poly->skew = 1.;

            algebraic = 1;
        }
        else
        {
            polys = (snfs_t *)malloc(sizeof(snfs_t));
            snfs_init(polys);
            npoly = 1;
            snfs_copy_poly(poly, polys);		// copy algebraic form

            // a^(7k) +/- 1 is divisible by (a^k +/- 1) giving a sextic in a^k
            polys->poly->alg.degree = 6;
            fobj->nfs_obj.pref_degree = 6;
            k = poly->exp1 / 7;
            mpz_set_ui(polys->c[6], 1);
            mpz_set_si(polys->c[5], -poly->coeff2);
            mpz_set_ui(polys->c[4], 1);
            mpz_set_si(polys->c[3], -poly->coeff2);
            mpz_set_ui(polys->c[2], 1);
            mpz_set_si(polys->c[1], -poly->coeff2);
            mpz_set_ui(polys->c[0], 1);
            mpz_set(m, poly->base1);
            polys->difficulty = log10(mpz_get_d(m)) * 6. * k;
            mpz_pow_ui(polys->poly->m, m, k);
            polys->poly->skew = 1.;

            algebraic = 1;
        }

    }
	else if ((poly->exp1 % 5 == 0) && (poly->coeff1 == 1) && (abs(poly->coeff2) == 1))
	{
		polys = (snfs_t *)malloc(sizeof(snfs_t));
		snfs_init(polys);
		npoly = 1;
		snfs_copy_poly(poly, polys);		// copy algebraic form

		// a^(5k) +/- 1 is divisible by (a^k +/- 1) giving a quartic in a^k
		polys->poly->alg.degree = 4;
		fobj->nfs_obj.pref_degree = 4;
		k = poly->exp1 / 5;
		mpz_set_ui(polys->c[4], 1);
		mpz_set_si(polys->c[3], -poly->coeff2);
		mpz_set_ui(polys->c[2], 1);
		mpz_set_si(polys->c[1], -poly->coeff2);
		mpz_set_ui(polys->c[0], 1);
		mpz_set(m, poly->base1);
		polys->difficulty = log10(mpz_get_d(m)) * 4. * k;
		mpz_pow_ui(polys->poly->m, m, k);
		polys->poly->skew = 1.;		

		algebraic = 1;
	}
	else if ((poly->exp1 % 7 == 0) && (poly->coeff1 == 1) && (abs(poly->coeff2) == 1) && 
        (mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10) > 90))
	{
		polys = (snfs_t *)malloc(sizeof(snfs_t));
		snfs_init(polys);
		npoly = 1;
		snfs_copy_poly(poly, polys);		// copy algebraic form

		// a^(7k) +/- 1 is divisible by (a^k +/- 1) giving a sextic in a^k
		polys->poly->alg.degree = 6;
		fobj->nfs_obj.pref_degree = 6;
		k = poly->exp1 / 7;
		mpz_set_ui(polys->c[6], 1);
		mpz_set_si(polys->c[5], -poly->coeff2);
		mpz_set_ui(polys->c[4], 1);
		mpz_set_si(polys->c[3], -poly->coeff2);
		mpz_set_ui(polys->c[2], 1);
		mpz_set_si(polys->c[1], -poly->coeff2);
		mpz_set_ui(polys->c[0], 1);
		mpz_set(m, poly->base1);
		polys->difficulty = log10(mpz_get_d(m)) * 6. * k;
		mpz_pow_ui(polys->poly->m, m, k);
		polys->poly->skew = 1.;		

		algebraic = 1;
	}
	else if ((poly->exp1 % 11 == 0) && (poly->coeff1 == 1) && (abs(poly->coeff2) == 1))
	{
		polys = (snfs_t *)malloc(sizeof(snfs_t));
		snfs_init(polys);
		npoly = 1;
		snfs_copy_poly(poly, polys);		// copy algebraic form

		// a^(11k) +/- 1 is divisible by (a^k +/- 1) giving a poly in a^k of degree 10.
		// but the poly is symmetric and can be halved to degree 5... 
		// see http://www.mersennewiki.org/index.php/SNFS_Polynomial_Selection
		polys->poly->alg.degree = 5;
		fobj->nfs_obj.pref_degree = 5;
		k = poly->exp1 / 11;
		mpz_set_ui(polys->c[5], 1);
		mpz_set_si(polys->c[4], -poly->coeff2*1);
		mpz_set_si(polys->c[3], -4);
		mpz_set_si(polys->c[2], poly->coeff2*3);
		mpz_set_ui(polys->c[1], 3);
		mpz_set_si(polys->c[0], -poly->coeff2);
		mpz_set(m, poly->base1);
		polys->difficulty = log10(mpz_get_d(m)) * 10. * k;
		mpz_pow_ui(polys->poly->m, m, k);
		polys->poly->skew = 1.; // is this right? c[0]/c[6] != 1	
				
		algebraic = 1;
		halved = 1;
	}
	else if ((poly->exp1 % 13 == 0) && (poly->coeff1 == 1) && (abs(poly->coeff2) == 1))
	{
		polys = (snfs_t *)malloc(sizeof(snfs_t));
		snfs_init(polys);
		npoly = 1;
		snfs_copy_poly(poly, polys);		// copy algebraic form

		// a^(13k) +/- 1 is divisible by (a^k +/- 1) giving a poly in a^k of degree 12.
		// but the poly is symmetric and can be halved to degree 6... 
		// see http://www.mersennewiki.org/index.php/SNFS_Polynomial_Selection
		polys->poly->alg.degree = 6;
		fobj->nfs_obj.pref_degree = 6;
		k = poly->exp1 / 13;
		mpz_set_ui(polys->c[6], 1);
		mpz_set_si(polys->c[5], -poly->coeff2);
		mpz_set_si(polys->c[4], -5);
		mpz_set_si(polys->c[3], poly->coeff2*4);
		mpz_set_ui(polys->c[2], 6);
		mpz_set_si(polys->c[1], -poly->coeff2*3);
		mpz_set_si(polys->c[0], -1);
		mpz_set(m, poly->base1);
		polys->difficulty = log10(mpz_get_d(m)) * 12. * k;
		mpz_pow_ui(polys->poly->m, m, k);
		polys->poly->skew = 1.;	
		
		algebraic = 1;
		halved = 1;
	}
	else 
	{
		mpz_t c0, cd, tmp;			

		// No algebraic factor - play with powers and composite bases
		int start_deg, stop_deg;

		mpz_init(c0);
		mpz_init(cd);
		mpz_init(tmp);

		mpz_set(m, b);
		if (!mpz_probab_prime_p(m,10))
			numf = tdiv_mpz(b, f, fobj->primes, fobj->num_p);

		// initialize candidate polynomials now that we know how many we'll need.
		// each factor of the base generates two, plus 2 for
		// the whole base raised and lowered, for each degree
        if (numf > 1)
        {
            apoly = (numf * 2 + 2) * 3;
        }
        else
        {
            apoly = 6;
        }

		polys = (snfs_t *)malloc(apoly * sizeof(snfs_t));
        for (i = 0; i < apoly; i++)
        {
            snfs_init(&polys[i]);
        }

		if (fobj->VFLAG > 0)
		{
			printf( "gen: ========================================================\n"
				"gen: considering the following polynomials:\n"
				"gen: ========================================================\n\n");
		}
		logprint_oc(fobj->nfs_obj.logfile, "a", "gen: considering the following polynomials:\n");
		
		npoly = 0;
		// be a little smarter about what degrees we consider... 
		me = (e-e % 5) / 5;
		//mpz_set_si(m, b);
		mpz_set(m, b);
		mpz_pow_ui(m, m, me);
		d = mpz_get_d(m);
		d = log10(d) * 5.;
		if (d < 120)
		{
			fobj->nfs_obj.pref_degree = 4;
			fobj->nfs_obj.alt_degree = 5;
			start_deg = 4; 
			stop_deg = 5;
		}
		else if (d < 170)
		{
			fobj->nfs_obj.pref_degree = 5;
			fobj->nfs_obj.alt_degree = 4;
			start_deg = 4; 
			stop_deg = 5;
		}
		else if (d < 220)
		{
			fobj->nfs_obj.pref_degree = 5;
			fobj->nfs_obj.alt_degree = 6;
			start_deg = 5; 
			stop_deg = 6;
		}
		else
		{
			fobj->nfs_obj.pref_degree = 6;
			fobj->nfs_obj.alt_degree = 5;
			start_deg = 5;
			stop_deg = 6;
		}

		fflush(stdout);

		for (i=start_deg; i<=stop_deg; i++)
		{			
			if (e % i == 0)
			{
				// the degree divides the exponent - resulting polynomial is straightforward
				me = e / i;
				mpz_set(m, b);		// signed?
				mpz_pow_ui(m, m, me);
				d = mpz_get_d(m);
				d = log10(d) * (double)i;
				snfs_copy_poly(poly, &polys[npoly]);		// copy algebraic form
				polys[npoly].difficulty = d;
				polys[npoly].poly->skew = 1.0;
				mpz_set_si(polys[npoly].c[i], poly->coeff1);
				mpz_set_si(polys[npoly].c[0], poly->coeff2);
				polys[npoly].poly->alg.degree = i;
				if (poly->form_type == SNFS_H_CUNNINGHAM)
				{
					mpz_set(polys[npoly].poly->rat.coeff[1], b2);
					mpz_pow_ui(polys[npoly].poly->rat.coeff[1], 
						   polys[npoly].poly->rat.coeff[1], me);
					mpz_set(polys[npoly].poly->rat.coeff[0], m);
					mpz_set(polys[npoly].poly->m, m);
					mpz_invert(n, polys[npoly].poly->rat.coeff[1], polys[npoly].n);
					mpz_mul(polys[npoly].poly->m, polys[npoly].poly->m, n);
					mpz_mod(polys[npoly].poly->m, polys[npoly].poly->m, polys[npoly].n);

					mpz_gcd(n, polys[npoly].poly->rat.coeff[0], polys[npoly].poly->rat.coeff[1]);
					mpz_tdiv_q(polys[npoly].poly->rat.coeff[0], polys[npoly].poly->rat.coeff[0], n);
					mpz_tdiv_q(polys[npoly].poly->rat.coeff[1], polys[npoly].poly->rat.coeff[1], n);

					mpz_neg(polys[npoly].poly->rat.coeff[1], polys[npoly].poly->rat.coeff[1]);					
				}
				else
				{
					mpz_set_si(polys[npoly].poly->rat.coeff[1], -1);
					mpz_set(polys[npoly].poly->rat.coeff[0], m);
					mpz_set(polys[npoly].poly->m, m);
				}

				mpz_gcd(tmp, polys[npoly].c[i], polys[npoly].c[0]);
				mpz_tdiv_q(polys[npoly].c[i], polys[npoly].c[i], tmp);
				mpz_tdiv_q(polys[npoly].c[0], polys[npoly].c[0], tmp);
				
				check_poly(&polys[npoly], fobj->VFLAG);
				approx_norms(&polys[npoly]);
				
				if (polys[npoly].valid)
				{				
					if (fobj->VFLAG > 0) print_snfs(&polys[npoly], stdout);
					npoly++;
				}
				else
				{	// being explicit
					snfs_clear(&polys[npoly]);
					snfs_init(&polys[npoly]);
				}
			}
			else
			{
				// degree does not divide the exponent, try increasing the exponent
				int inc = (i - e % i);
				me = (e+inc) / i;
				mpz_set(m, b);		// signed?
				mpz_pow_ui(m, m, me);
				d = mpz_get_d(m);
				d = log10(d) * (double)i;

				// cd = (int64)pow((double)b2, inc) * (int64)poly->coeff1;
				mpz_set_si(cd, poly->coeff1);
				mpz_pow_ui(tmp, b2, inc);
				mpz_mul(cd, cd, tmp);
				 
				// c0 = (int64)pow((double)b,inc) * (int64)poly->coeff2;
				mpz_set_si(c0, poly->coeff2);
				mpz_pow_ui(tmp, b, inc);
				mpz_mul(c0, c0, tmp);

                // leading coefficient contributes to the difficulty
                d += log10(mpz_get_d(cd));

				//skew = pow((double)abs((int)c0)/(double)cd, 1./(double)i);
                // leading coefficient contributes to the difficulty
                d += log10(mpz_get_d(cd));

				skew = pow(fabs(mpz_get_d(c0)) / mpz_get_d(cd), 1./(double)i);
				snfs_copy_poly(poly, &polys[npoly]);		// copy algebraic form
				polys[npoly].difficulty = d;
				polys[npoly].poly->skew = skew;
				mpz_set(polys[npoly].c[i], cd);
				mpz_set(polys[npoly].c[0], c0);
				polys[npoly].poly->alg.degree = i;
				if (poly->form_type == SNFS_H_CUNNINGHAM)
				{
					mpz_set(polys[npoly].poly->rat.coeff[1], b2);		// signed?
					mpz_pow_ui(polys[npoly].poly->rat.coeff[1], 
						   polys[npoly].poly->rat.coeff[1], me);
					mpz_set(polys[npoly].poly->rat.coeff[0], m);
					mpz_set(polys[npoly].poly->m, m);
					mpz_invert(n, polys[npoly].poly->rat.coeff[1], polys[npoly].n);
					mpz_mul(polys[npoly].poly->m, polys[npoly].poly->m, n);
					mpz_mod(polys[npoly].poly->m, polys[npoly].poly->m, polys[npoly].n);

					mpz_gcd(n, polys[npoly].poly->rat.coeff[0], polys[npoly].poly->rat.coeff[1]);
					mpz_tdiv_q(polys[npoly].poly->rat.coeff[0], polys[npoly].poly->rat.coeff[0], n);
					mpz_tdiv_q(polys[npoly].poly->rat.coeff[1], polys[npoly].poly->rat.coeff[1], n);

					mpz_neg(polys[npoly].poly->rat.coeff[1], polys[npoly].poly->rat.coeff[1]);					
				}
				else
				{
					mpz_set_si(polys[npoly].poly->rat.coeff[1], -1);
					mpz_set(polys[npoly].poly->rat.coeff[0], m);
					mpz_set(polys[npoly].poly->m, m);
				}
		
				mpz_gcd(tmp, polys[npoly].c[i], polys[npoly].c[0]);
				mpz_tdiv_q(polys[npoly].c[i], polys[npoly].c[i], tmp);
				mpz_tdiv_q(polys[npoly].c[0], polys[npoly].c[0], tmp);

				check_poly(&polys[npoly], fobj->VFLAG);
				approx_norms(&polys[npoly]);
				
				if (polys[npoly].valid)
				{					
					if (fobj->VFLAG > 0) print_snfs(&polys[npoly], stdout);
					npoly++;
				}
				else
				{	// being explicit
					snfs_clear(&polys[npoly]);
					snfs_init(&polys[npoly]);
				}

				// and decreasing the exponent
				inc = e % i;
				me = (e-inc) / i;
				mpz_set(m, b);		// signed?
				mpz_pow_ui(m, m, me);
				d = mpz_get_d(m);

                // thanks to jyb for finding that the leading coefficient was
                // double-counted in this case.
                // https://www.mersenneforum.org/showpost.php?p=571720&postcount=1
                d = log10(d) * (double)i; // +log10(pow(mpz_get_d(b), inc));

				mpz_set_si(cd, poly->coeff1);
				mpz_pow_ui(tmp, b, inc);
				mpz_mul(cd, cd, tmp);
				
				mpz_set_si(c0, poly->coeff2);
				mpz_pow_ui(tmp, b2, inc);
				mpz_mul(c0, c0, tmp);
				
				//skew = pow((double)abs((int)c0)/(double)cd, 1./(double)i);
				skew = pow(fabs(mpz_get_d(c0)) / mpz_get_d(cd), 1./(double)i);

                //gmp_printf("cd: %Zd\nc0: %Zd\n", cd, c0);
                //gmp_printf("b1: %Zd\nb2: %Zd\nm: %Zd\n", b, b2, m);
                //printf("exp: %d\nmul: %d\nd  : %lf\n", e, i, d);
                //printf("coeff1: %d\ncoeff2: %d\n", poly->coeff1, poly->coeff2);

				// leading coefficient contributes to the difficulty
				d += log10(mpz_get_d(cd));
                //printf("skew: %lf\n", skew);
                //printf("diff: %lf\n", d);

				snfs_copy_poly(poly, &polys[npoly]);		// copy algebraic form
				polys[npoly].difficulty = d;
				polys[npoly].poly->skew = skew;
				mpz_set(polys[npoly].c[i], cd);
				mpz_set(polys[npoly].c[0], c0);
				polys[npoly].poly->alg.degree = i;
				if (poly->form_type == SNFS_H_CUNNINGHAM)
				{
					mpz_set(polys[npoly].poly->rat.coeff[1], b2);		// signed?
					mpz_pow_ui(polys[npoly].poly->rat.coeff[1],
						   polys[npoly].poly->rat.coeff[1], me);
					mpz_set(polys[npoly].poly->rat.coeff[0], m);
					mpz_set(polys[npoly].poly->m, m);
					mpz_invert(n, polys[npoly].poly->rat.coeff[1], polys[npoly].n);
					mpz_mul(polys[npoly].poly->m, polys[npoly].poly->m, n);
					mpz_mod(polys[npoly].poly->m, polys[npoly].poly->m, polys[npoly].n);

					mpz_gcd(n, polys[npoly].poly->rat.coeff[0], polys[npoly].poly->rat.coeff[1]);
					mpz_tdiv_q(polys[npoly].poly->rat.coeff[0], polys[npoly].poly->rat.coeff[0], n);
					mpz_tdiv_q(polys[npoly].poly->rat.coeff[1], polys[npoly].poly->rat.coeff[1], n);

					mpz_neg(polys[npoly].poly->rat.coeff[1], polys[npoly].poly->rat.coeff[1]);
				}
				else
				{
					mpz_set_si(polys[npoly].poly->rat.coeff[1], -1);
					mpz_set(polys[npoly].poly->rat.coeff[0], m);
					mpz_set(polys[npoly].poly->m, m);
				}

				mpz_gcd(tmp, polys[npoly].c[i], polys[npoly].c[0]);
				mpz_tdiv_q(polys[npoly].c[i], polys[npoly].c[i], tmp);
				mpz_tdiv_q(polys[npoly].c[0], polys[npoly].c[0], tmp);
				
				check_poly(&polys[npoly], fobj->VFLAG);
				approx_norms(&polys[npoly]);
				
				if (polys[npoly].valid)
				{					
					if (fobj->VFLAG > 0) print_snfs(&polys[npoly], stdout);
					npoly++;
				}
				else
				{	// being explicit
					snfs_clear(&polys[npoly]);
					snfs_init(&polys[npoly]);
				}

				// and playing with composite bases
				if (numf > 1)
				{
					int j;

					// multiply by powers of each factor individually to get that 
					// factor to a multiple of degree
					for (j=0; j<numf; j++)
					{
						int k, i1, i2, bb;

						// unique factor
						if (j > 0)
							if (f[j] == f[j-1])
								continue;

						// move it up
						i1 = (i - e % i);
						
						//c0 = (int64)pow((double)f[j], i1) * (int64)poly->coeff2;
						mpz_set_si(c0, poly->coeff2);
						mpz_set_si(tmp, f[j]);
						mpz_pow_ui(tmp, tmp, i1);
						mpz_mul(c0, c0, tmp);

						// since we moved the current factor up, move the other factors down
						// otherwise this would be the same as moving the whole composite base up.
						// also move the second term up...
						
						//cd = (int64)pow((double)b2, i1) * (int64)poly->coeff1;
						mpz_set_si(cd, poly->coeff1);
						mpz_pow_ui(tmp, b2, i1);
						mpz_mul(cd, cd, tmp);


						bb = 1;
						i2 = e % i;
						for (k=0; k<numf; k++)
						{
							if (k == j) continue;
							//cd *= (int64)pow((double)f[k], i2);
							mpz_set_si(tmp, f[k]);
							mpz_pow_ui(tmp, tmp, i2);
							mpz_mul(cd, cd, tmp);

							bb *= f[k];
						}
						// m is now a mix of powers of factors of b.
						// here is the contribution of the factor we increased
						me = (e+i1) / i;
						mpz_set_si(m, f[j]);
						mpz_pow_ui(m, m, me);
						// here is the contribution of the factors we decreased
						me = (e-i2) / i;
						mpz_set_si(n, bb);
						mpz_pow_ui(n, n, me);
						// combine them and compute the base difficulty
						mpz_mul(m, m, n);
						d = mpz_get_d(m);
						d = log10(d) * (double)i;
						// the power we moved up appears in the constant term and the powers we moved 
						// down appear as a coefficient to the high order term and thus contribute
						// to the difficulty.
						d += log10(mpz_get_d(cd));
						skew = pow(fabs(mpz_get_d(c0))/mpz_get_d(cd), 1./(double)i);
						snfs_copy_poly(poly, &polys[npoly]);		// copy algebraic form
						polys[npoly].difficulty = d;
						polys[npoly].poly->skew = skew;
						mpz_set(polys[npoly].c[i], cd);
						mpz_set(polys[npoly].c[0], c0);
						polys[npoly].poly->alg.degree = i;
						if (poly->form_type == SNFS_H_CUNNINGHAM)
						{
							mpz_set(polys[npoly].poly->rat.coeff[1], b2);		// signed?
							mpz_pow_ui(polys[npoly].poly->rat.coeff[1],
								   polys[npoly].poly->rat.coeff[1], (e+i1) / i);
							mpz_set(polys[npoly].poly->rat.coeff[0], m);
							mpz_set(polys[npoly].poly->m, m);
							mpz_invert(n, polys[npoly].poly->rat.coeff[1], polys[npoly].n);
							mpz_mul(polys[npoly].poly->m, polys[npoly].poly->m, n);
							mpz_mod(polys[npoly].poly->m, polys[npoly].poly->m, polys[npoly].n);

							mpz_gcd(n, polys[npoly].poly->rat.coeff[0], polys[npoly].poly->rat.coeff[1]);
							mpz_tdiv_q(polys[npoly].poly->rat.coeff[0], polys[npoly].poly->rat.coeff[0], n);
							mpz_tdiv_q(polys[npoly].poly->rat.coeff[1], polys[npoly].poly->rat.coeff[1], n);

							mpz_neg(polys[npoly].poly->rat.coeff[1], polys[npoly].poly->rat.coeff[1]);							
						}
						else
						{
							mpz_set_si(polys[npoly].poly->rat.coeff[1], -1);
							mpz_set(polys[npoly].poly->rat.coeff[0], m);
							mpz_set(polys[npoly].poly->m, m);
						}

						mpz_gcd(tmp, polys[npoly].c[i], polys[npoly].c[0]);
						mpz_tdiv_q(polys[npoly].c[i], polys[npoly].c[i], tmp);
						mpz_tdiv_q(polys[npoly].c[0], polys[npoly].c[0], tmp);

						check_poly(&polys[npoly], fobj->VFLAG);
						approx_norms(&polys[npoly]);
						
						if (polys[npoly].valid)
						{							
							if (fobj->VFLAG > 0) print_snfs(&polys[npoly], stdout);
							npoly++;
						}
						else
						{	// being explicit
							snfs_clear(&polys[npoly]);
							snfs_init(&polys[npoly]);
						}

						// move it down
						i1 = e % i;
						
						//cd = (int64)pow((double)f[j], i1) * (int64)poly->coeff1;
						mpz_set_si(cd, poly->coeff1);
						mpz_set_si(tmp, f[j]);
						mpz_pow_ui(tmp, tmp, i1);
						mpz_mul(cd, cd, tmp);

						// since we moved the current factor down, move the other factors up
						// otherwise this would be the same as moving the whole composite base down.
						
						//c0 = (int64)pow((double)b2, i1) * (int64)poly->coeff2;
						mpz_set_si(c0, poly->coeff2);
						mpz_pow_ui(tmp, b2, i1);
						mpz_mul(c0, c0, tmp);

						bb = 1;
						i2 = (i - e % i);
						for (k=0; k<numf; k++)
						{
							if (k == j) continue;
							//c0 *= (int64)pow((double)f[k], i2);
							mpz_set_si(tmp, f[k]);
							mpz_pow_ui(tmp, tmp, i2);
							mpz_mul(c0, c0, tmp);

							bb *= f[k];
						}
						// m is now a mix of powers of factors of b.
						// here is the contribution of the factor we cecreased
						me = (e-i1) / i;
						mpz_set_si(m, f[j]);
						mpz_pow_ui(m, m, me);
						// here is the contribution of the factors we increased
						me = (e+i2) / i;
						mpz_set_si(n, bb);
						mpz_pow_ui(n, n, me);
						// combine them and compute the base difficulty
						mpz_mul(m, m, n);
						d = mpz_get_d(m);
						d = log10(d) * (double)i;
						// the power we moved up appears in the constant term and the powers we moved 
						// down appear as a coefficient to the high order term and thus contribute
						// to the difficulty.
						d += log10(mpz_get_d(cd));
						skew = pow(fabs(mpz_get_d(c0))/mpz_get_d(cd), 1./(double)i);
						snfs_copy_poly(poly, &polys[npoly]);		// copy algebraic form
						polys[npoly].difficulty = d;
						polys[npoly].poly->skew = skew;
						mpz_set(polys[npoly].c[i], cd);
						mpz_set(polys[npoly].c[0], c0);
						polys[npoly].poly->alg.degree = i;
						if (poly->form_type == SNFS_H_CUNNINGHAM)
						{
							mpz_set(polys[npoly].poly->rat.coeff[1], b2);		// signed?
							mpz_pow_ui(polys[npoly].poly->rat.coeff[1],
								   polys[npoly].poly->rat.coeff[1], (e+i1) / i);
							mpz_set(polys[npoly].poly->rat.coeff[0], m);
							mpz_set(polys[npoly].poly->m, m);
							mpz_invert(n, polys[npoly].poly->rat.coeff[1], polys[npoly].n);
							mpz_mul(polys[npoly].poly->m, polys[npoly].poly->m, n);
							mpz_mod(polys[npoly].poly->m, polys[npoly].poly->m, polys[npoly].n);

							mpz_gcd(n, polys[npoly].poly->rat.coeff[0], polys[npoly].poly->rat.coeff[1]);
							mpz_tdiv_q(polys[npoly].poly->rat.coeff[0], polys[npoly].poly->rat.coeff[0], n);
							mpz_tdiv_q(polys[npoly].poly->rat.coeff[1], polys[npoly].poly->rat.coeff[1], n);

							mpz_neg(polys[npoly].poly->rat.coeff[1], polys[npoly].poly->rat.coeff[1]);
						}
						else
						{
							mpz_set_si(polys[npoly].poly->rat.coeff[1], -1);
							mpz_set(polys[npoly].poly->rat.coeff[0], m);
							mpz_set(polys[npoly].poly->m, m);
						}

						mpz_gcd(tmp, polys[npoly].c[i], polys[npoly].c[0]);
						mpz_tdiv_q(polys[npoly].c[i], polys[npoly].c[i], tmp);
						mpz_tdiv_q(polys[npoly].c[0], polys[npoly].c[0], tmp);

						check_poly(&polys[npoly], fobj->VFLAG);
						approx_norms(&polys[npoly]);
						
						if (polys[npoly].valid)
						{							
							if (fobj->VFLAG > 0) print_snfs(&polys[npoly], stdout);
							npoly++;
						}
						else
						{	// being explicit
							snfs_clear(&polys[npoly]);
							snfs_init(&polys[npoly]);
						}
					} // loop over factors of base
				} // composite base?
			} // degree divides exponent?
		} // for each degree

		mpz_clear(c0);
		mpz_clear(cd);
		mpz_clear(tmp);
	} // check for algebraic factors

	if (algebraic)
	{
		if (halved)
		{
			if (poly->form_type == SNFS_H_CUNNINGHAM)
			{
				// multiplying through by b = a^k gives g(x) = bx - (b^2 + 1).  but
				// for homogeneous polys a = a1/a2 so we have g(x) = (a1/a2)^k * x - ((a1^2/a2^2)^k + 1) = 0
				// multiplying through by a2^(2k) gives:
				// g(x) = (a1*a2)^k*x - (a1^(2k) + a2^(2k)) with m = (a1/a2)^k + (a2/a1)^k
				mpz_set(polys->poly->rat.coeff[1], b2);		// signed?
				mpz_mul(polys->poly->rat.coeff[1], polys->poly->rat.coeff[1], b);
				mpz_pow_ui(polys->poly->rat.coeff[1], polys->poly->rat.coeff[1], k);
				mpz_mul(polys->poly->rat.coeff[0], polys->poly->m, polys->poly->m);		// a1^(2k)
				mpz_set(n, poly->base2);
				mpz_pow_ui(n, n, 2*k);
				mpz_add(polys->poly->rat.coeff[0], polys->poly->rat.coeff[0], n);
				mpz_sqrt(n, n);
				mpz_set(t, n);													// a2^k
				mpz_invert(n, n, poly->n);										// a2^-k
				mpz_mul(n, polys->poly->m, n);									// 
				mpz_mod(n, n, poly->n);											// (a1/a2)^k
				mpz_invert(polys->poly->m, polys->poly->m, poly->n);			// a1^-k
				mpz_mul(polys->poly->m, polys->poly->m, t);						// 
				mpz_mod(polys->poly->m, polys->poly->m, poly->n);				// (a2/a1)^k
				mpz_add(polys->poly->m, polys->poly->m, n);						//
				mpz_mod(polys->poly->m, polys->poly->m, poly->n);				// (a1/a2)^k + (a2/a1)^k
				mpz_neg(polys->poly->rat.coeff[1], polys->poly->rat.coeff[1]);
			}
			else
			{
				// for halved degree polynomials, Y1 becomes -x^k and Y0 becomes x^(2k) + 1
				// such that y1*x + y0 evaluated at M = x^k + x^-k is 0.
				mpz_neg(polys->poly->rat.coeff[1], polys->poly->m);
				mpz_mul(polys->poly->rat.coeff[0], polys->poly->m, polys->poly->m);
				mpz_add_ui(polys->poly->rat.coeff[0], polys->poly->rat.coeff[0], 1);
				mpz_invert(m, polys->poly->m, polys->n);
				mpz_add(polys->poly->m, polys->poly->m, m);
			}
		}
		else
		{
			if (poly->form_type == SNFS_H_CUNNINGHAM)
			{
				mpz_set(polys->poly->rat.coeff[1], b2);		// signed?
				mpz_pow_ui(polys->poly->rat.coeff[1], polys->poly->rat.coeff[1], k);
				mpz_set(polys->poly->rat.coeff[0], polys->poly->m);
				mpz_invert(n, polys->poly->rat.coeff[1], poly->n);
				mpz_mul(polys->poly->m, polys->poly->m, n);
				mpz_mod(polys->poly->m, polys->poly->m, poly->n);
				mpz_neg(polys->poly->rat.coeff[1], polys->poly->rat.coeff[1]);
			}
			else
			{
				// Y1 = -1, Y0 = m such that y1*m + y0 = 0
				mpz_set(polys->poly->rat.coeff[0], polys->poly->m);
				mpz_set_si(polys->poly->rat.coeff[1], -1);
			}			
		}
		check_poly(polys, fobj->VFLAG);
		approx_norms(polys);

        if (!polys->valid)
        {	// being explicit
            snfs_clear(polys);
            npoly = 0;
        }
	}

	mpz_clear(m);
	mpz_clear(n);
	mpz_clear(t);
	mpz_clear(b);
	mpz_clear(b2);
	*npolys = npoly;
	return polys;
}

snfs_t* gen_xyyxf_poly(fact_obj_t *fobj, snfs_t *poly, int* npolys)
{
	int deg, i, j, nump1, nump2, me, base, e, b;
	// xyyxf bases never exceed a small number, by definition
	int x = mpz_get_ui(poly->base1), y = mpz_get_ui(poly->base2);
	mpz_t n, m;
	double d, skew;
	int f1[100];
	int numf1 = 0;
	int f2[100];
	int numf2 = 0;
	snfs_t *polys, *final_polys;
	int npoly = 0;
	int apoly;	
    int alloc_base_poly;
	FILE *f;
	double avg_diff;

	mpz_init(n);
	mpz_init(m);

	// xyyxf numbers take the form x^y + y^x, where 1<y<x<151
	// as far as I know there are no simple algebraic reductions of these types of polynomials.
	// The strategy for generating polynomials for them is therefore the following:
	// for each of the two terms, x^y and y^x, raise and lower the exponents to multiples of 
	// each of degrees 4, 5, and 6 (as in brent poly generation).  keep track of the coefficients
	// these actions create for each of the two terms.  Each can be raised/lowered independently,
	// so there will be quite a few possible polynomials.  We also work with composite bases
	// similarly to brent poly generation.  
	// Once we have a suitable degree, the polynomial will look like this:
	// a*(x^yd)^d + b*(y^xd)^d, where xd and yd are x and y reduced by a factor of d, 
	// and a,b are coefficients arising from adjusting x,y to multiples of d.
	// the algebraic poly f(z) is then a*z^d + b, where z = x^yd/y^xd
	// the rational poly g(z) is then -(y^xd)*z + (x^yd)
	// To handle the case where, with composite bases, we
	// have common factors, reduce Y1 and Y0 by GCD(Y1, Y0) and reduce cx, c0 by gcd(cx, c0).

	mpz_set_ui(m, x);
	if (!mpz_probab_prime_p(m,10))
		numf1 = tdiv_int(x, f1, fobj->primes, fobj->num_p);

	mpz_set_ui(m, y);
	if (!mpz_probab_prime_p(m,10))
		numf2 = tdiv_int(y, f2, fobj->primes, fobj->num_p);

	// initialize candidate polynomials now that we know how many we'll need.
	// each factor of the base generates two, plus 2 for
	// the whole base raised and lowered, for each degree, for each base
	nump1 = (numf1 * 2 + 2) * 3;
	nump2 = (numf2 * 2 + 2) * 3;
	apoly = nump1 + nump2;

    // so we can free it later.  apoly gets reused.
    alloc_base_poly = apoly;

	printf("number of factors: %d, %d, total polys = %d\n", numf1, numf2, apoly);

	polys = (snfs_t *)malloc(apoly * sizeof(snfs_t));
	for (i=0; i<apoly; i++)
	{
		snfs_init(&polys[i]);
		polys[i].valid = 0;
	}
		
    if (fobj->LOGFLAG)
    {
        f = fopen(fobj->flogname, "a");
        if (f != NULL)
        {
            logprint(f, "nfs: commencing snfs on c%d: ", gmp_base10(poly->n));
            gmp_fprintf(f, "%Zd\n", poly->n);
            fclose(f);
        }
    }

	npoly = 0;
	// form all polys for x^y + 1 and 1 + y^x separately.  Then combine them later.
	for (base = 0; base < 2; base++)
	{
		int *f, b2, numf;
		// set the exponent/base.  on the first pass 
		// we process x^y and the next we process y^x
		e = base == 0 ? y : x;
		b = base == 0 ? x : y;
		f = base == 0 ? f1 : f2;
		numf = base == 0 ? numf1 : numf2;
		b2 = 1;

		for (deg=4; deg<7; deg++)
		{
			int64_t c0, cd;

			if (e % deg == 0)
			{
				// the degree divides the exponent - resulting polynomial is straightforward
				me = e / deg;
				mpz_set_si(m, b);
				mpz_pow_ui(m, m, me);
				// remember the base and exponent
				mpz_set_si(polys[npoly].base1, b);
				polys[npoly].exp1 = me;
				mpz_set_si(polys[npoly].base2, 1);
				polys[npoly].exp2 = 1;
				d = mpz_get_d(m);
				d = log10(d) * (double)deg;
				mpz_set(polys[npoly].n, poly->n);
				polys[npoly].difficulty = d;
				polys[npoly].poly->skew = 1.0;
				mpz_set_si(polys[npoly].c[deg], poly->coeff1);
				mpz_set_si(polys[npoly].c[0], poly->coeff2);
				polys[npoly].poly->alg.degree = deg;
				npoly++;
			}
			else
			{
				// degree does not divide the exponent, try increasing the exponent
				int inc = (deg - e % deg);
				me = (e+inc) / deg;
				mpz_set_si(m, b);
				mpz_pow_ui(m, m, me);
				// remember the base and exponent
				mpz_set_si(polys[npoly].base1, b);
				polys[npoly].exp1 = me;
				mpz_set_si(polys[npoly].base2, 1);
				polys[npoly].exp2 = 1;
				d = mpz_get_d(m);
				d = log10(d) * (double)deg;
				cd = (int64_t)pow((double)b2, inc) * poly->coeff1;
				c0 = (int64_t)pow((double)b,inc) * poly->coeff2;
				skew = pow((double)abs(c0)/(double)cd, 1./(double)deg);
				mpz_set(polys[npoly].n, poly->n);
				polys[npoly].difficulty = d;
				polys[npoly].poly->skew = skew;
				mpz_set_si(polys[npoly].c[deg], cd);
				mpz_set_si(polys[npoly].c[0], c0);
				polys[npoly].poly->alg.degree = deg;
				npoly++;

				// and decreasing the exponent
				inc = e % deg;
				me = (e-inc) / deg;
				mpz_set_si(m, b);
				mpz_pow_ui(m, m, me);
				// remember the base and exponent
				mpz_set_si(polys[npoly].base1, b);
				polys[npoly].exp1 = me;
				mpz_set_si(polys[npoly].base2, 1);
				polys[npoly].exp2 = 1;
				d = mpz_get_d(m);
				d = log10(d) * (double)deg + log10(pow((double)b,inc));
				cd = (int64_t)pow((double)b,inc) * poly->coeff1;
				c0 = (int64_t)pow((double)b2, inc) * poly->coeff2;
				skew = pow((double)abs(c0)/(double)cd, 1./(double)deg);
				// leading coefficient contributes to the difficulty
				//d += log10((double)cd);
				mpz_set(polys[npoly].n, poly->n);
				polys[npoly].difficulty = d;
				polys[npoly].poly->skew = skew;
				mpz_set_si(polys[npoly].c[deg], cd);
				mpz_set_si(polys[npoly].c[0], c0);
				polys[npoly].poly->alg.degree = deg;
				npoly++;

				// and playing with composite bases
				if (numf > 1)
				{
					int j;

					// multiply by powers of each factor individually to get that 
					// factor to a multiple of degree
					for (j=0; j<numf; j++)
					{
						int k, i1, i2, bb;

						// unique factors...
						if (j > 0)
							if (f[j] == f[j-1])
								continue;

						// move it up
						i1 = (deg - e % deg);
						c0 = pow((double)f[j], i1) * poly->coeff2;
						// since we moved the current factor up, move the other factors down
						// otherwise this would be the same as moving the whole composite base up.
						// also move the second term up...
						cd = pow((double)b2, i1) * poly->coeff1;
						bb = 1;
						i2 = e % deg;
						for (k=0; k<numf; k++)
						{
							if (k == j) continue;
							cd *= pow((double)f[k], i2);
							bb *= f[k];
						}
						// m is now a mix of powers of factors of b.
						// here is the contribution of the factor we increased
						me = (e+i1) / deg;
						mpz_set_si(m, f[j]);
						mpz_pow_ui(m, m, me);
						// remember the base and exponent
						mpz_set_si(polys[npoly].base1, f[j]);
						polys[npoly].exp1 = me;
						// here is the contribution of the factors we decreased
						me = (e-i2) / deg;
						mpz_set_si(n, bb);
						mpz_pow_ui(n, n, me);
						// remember the base and exponent
						mpz_set_si(polys[npoly].base2, bb);
						polys[npoly].exp2 = me;
						// combine them and compute the base difficulty
						mpz_mul(m, m, n);
						d = mpz_get_d(m);
						d = log10(d) * (double)deg;
						// the power we moved up appears in the constant term and the powers we moved 
						// down appear as a coefficient to the high order term and thus contribute
						// to the difficulty.
						d += log10((double)cd);
						skew = pow((double)abs(c0)/(double)cd, 1./(double)deg);
						mpz_set(polys[npoly].n, poly->n);
						polys[npoly].difficulty = d;
						polys[npoly].poly->skew = skew;
						mpz_set_si(polys[npoly].c[deg], cd);
						mpz_set_si(polys[npoly].c[0], c0);
						mpz_set(polys[npoly].poly->m, m);
						polys[npoly].poly->alg.degree = deg;
						npoly++;

						// move it down
						i1 = e % deg;
						cd = pow((double)f[j], i1) * poly->coeff1;
						// since we moved the current factor down, move the other factors up
						// otherwise this would be the same as moving the whole composite base down.
						c0 = pow((double)b2, i1) * poly->coeff2;
						bb = 1;
						i2 = (deg - e % deg);
						for (k=0; k<numf; k++)
						{
							if (k == j) continue;
							c0 *= pow((double)f[k], i2);
							bb *= f[k];
						}
						// m is now a mix of powers of factors of b.
						// here is the contribution of the factor we increased
						me = (e-i1) / deg;
						mpz_set_si(m, f[j]);
						mpz_pow_ui(m, m, me);
						// remember the base and exponent
						mpz_set_si(polys[npoly].base1, f[j]);
						polys[npoly].exp1 = me;
						// here is the contribution of the factors we decreased
						me = (e+i2) / deg;
						mpz_set_si(n, bb);
						mpz_pow_ui(n, n, me);
						// remember the base and exponent
						mpz_set_si(polys[npoly].base2, bb);
						polys[npoly].exp2 = me;
						// combine them and compute the base difficulty
						mpz_mul(m, m, n);
						d = mpz_get_d(m);
						d = log10(d) * (double)deg;
						// the power we moved up appears in the constant term and the powers we moved 
						// down appear as a coefficient to the high order term and thus contribute
						// to the difficulty.
						d += log10((double)cd);
						skew = pow((double)abs(c0)/(double)cd, 1./(double)deg);
						mpz_set(polys[npoly].n, poly->n);
						polys[npoly].difficulty = d;
						polys[npoly].poly->skew = skew;
						mpz_set_si(polys[npoly].c[deg], cd);
						mpz_set_si(polys[npoly].c[0], c0);
						mpz_set(polys[npoly].poly->m, m);
						polys[npoly].poly->alg.degree = deg;
						npoly++;
					} // loop over factors of base
				} // composite base?
			} // degree divides exponent?
		} // for each degree

		// record how many polys this base actually generated (depends on whether the exp
		// is divisible by any of the various degrees or not)
		// this matters because we need the actual number of polys for this base
		// to use as an offset into the other base's polynomials below...
		if (base == 0) 
			nump1 = npoly;
		else
			nump2 = npoly - nump1;
	} // for each base (x, y)

	apoly = nump1 * nump2;

	printf("actual polys = %d, %d, total actual polys = %d\n", nump1, nump2, apoly);

	final_polys = (snfs_t *)malloc(apoly * sizeof(snfs_t));
    for (i = 0; i < apoly; i++)
    {
        snfs_init(&final_polys[i]);
    }

	if (fobj->VFLAG > 0)
	{
		printf( "\ngen: ========================================================\n"
			"gen: considering the following polynomials:\n"
			"gen: ========================================================\n\n");
	}

	// now mix together the possible forms of x^y + 1 and 1 + y^x	
	npoly = 0;
	avg_diff = 0.;
	for (i=0; i<nump1; i++)
	{
		snfs_t *p1, *p2;
		int64_t c0, cd;

		p1 = &polys[i];
		for (j=0; j<nump2; j++)
		{
			p2 = &polys[nump1+j];

			// don't mix degrees.
			if (p1->poly->alg.degree != p2->poly->alg.degree)
				continue;

			deg = p1->poly->alg.degree;

			cd = mpz_get_si(p1->c[deg]) * mpz_get_si(p2->c[0]);
			c0 = mpz_get_si(p1->c[0]) * mpz_get_si(p2->c[deg]);

			// whichever of these is smaller can be the leading coefficient
			if (c0 > cd) 
			{
				mpz_set(final_polys[npoly].n, poly->n);
				mpz_set_si(final_polys[npoly].c[deg], cd / spGCD(c0, cd));
				mpz_set_si(final_polys[npoly].c[0], c0 / spGCD(c0, cd));
				final_polys[npoly].poly->skew = pow(
					fabs(mpz_get_d(final_polys[npoly].c[0]))/mpz_get_d(final_polys[npoly].c[deg]), 1./(double)deg);

				final_polys[npoly].poly->alg.degree = deg;
				final_polys[npoly].difficulty = log(pow(2.71828, p1->difficulty) + pow(2.71828, p2->difficulty));

				// the algebraic poly f(z) is then cd*z^d + c0, where z = x^yd/y^xd
				// the rational poly g(z) is then -(y^xd)*z + (x^yd)
				// to handle composite bases, the first loop records the part that was raised and 
				// the part that was lowered separately
				mpz_set(final_polys[npoly].poly->rat.coeff[1], p2->base1);
				mpz_pow_ui(final_polys[npoly].poly->rat.coeff[1], 
						final_polys[npoly].poly->rat.coeff[1], p2->exp1);
				mpz_set(n, p2->base2);
				mpz_pow_ui(n, n, p2->exp2);
				mpz_mul(final_polys[npoly].poly->rat.coeff[1], final_polys[npoly].poly->rat.coeff[1], n);

				mpz_set(final_polys[npoly].poly->rat.coeff[0], p1->base1);
				mpz_pow_ui(final_polys[npoly].poly->rat.coeff[0], 
						final_polys[npoly].poly->rat.coeff[0], p1->exp1);
				mpz_set(n, p1->base2);
				mpz_pow_ui(n, n, p1->exp2);
				mpz_mul(final_polys[npoly].poly->rat.coeff[0], final_polys[npoly].poly->rat.coeff[0], n);

				mpz_gcd(n, final_polys[npoly].poly->rat.coeff[0], final_polys[npoly].poly->rat.coeff[1]);
				final_polys[npoly].difficulty -= log10(mpz_get_d(n))*deg;
				mpz_tdiv_q(final_polys[npoly].poly->rat.coeff[0], final_polys[npoly].poly->rat.coeff[0], n);
				mpz_tdiv_q(final_polys[npoly].poly->rat.coeff[1], final_polys[npoly].poly->rat.coeff[1], n);

				mpz_set(final_polys[npoly].poly->m, final_polys[npoly].poly->rat.coeff[0]);
				mpz_invert(n, final_polys[npoly].poly->rat.coeff[1], final_polys[npoly].n);
				mpz_mul(final_polys[npoly].poly->m, final_polys[npoly].poly->m, n);
				mpz_mod(final_polys[npoly].poly->m, final_polys[npoly].poly->m, final_polys[npoly].n);
				mpz_neg(final_polys[npoly].poly->rat.coeff[1], final_polys[npoly].poly->rat.coeff[1]);
			}
			else
			{
				mpz_set(final_polys[npoly].n, poly->n);
				mpz_set_si(final_polys[npoly].c[deg], c0 / spGCD(c0, cd));
				mpz_set_si(final_polys[npoly].c[0], cd / spGCD(c0, cd));
				final_polys[npoly].poly->skew = pow(
					fabs(mpz_get_d(final_polys[npoly].c[0]))/mpz_get_d(final_polys[npoly].c[deg]), 1./(double)deg);
				final_polys[npoly].poly->alg.degree = deg;
				final_polys[npoly].difficulty = log(pow(2.71828, p1->difficulty) + pow(2.71828, p2->difficulty));

				// the algebraic poly f(z) is then cd*z^d + c0, where z = y^xd/x^yd
				// the rational poly g(z) is then -(x^yd)*z + (y^xd)
				// to handle composite bases, the first loop records the part that was raised and 
				// the part that was lowered separately
				mpz_set(final_polys[npoly].poly->rat.coeff[1], p1->base1);
				mpz_pow_ui(final_polys[npoly].poly->rat.coeff[1], 
						final_polys[npoly].poly->rat.coeff[1], p1->exp1);
				mpz_set(n, p1->base2);
				mpz_pow_ui(n, n, p1->exp2);
				mpz_mul(final_polys[npoly].poly->rat.coeff[1], final_polys[npoly].poly->rat.coeff[1], n);

				mpz_set(final_polys[npoly].poly->rat.coeff[0], p2->base1);
				mpz_pow_ui(final_polys[npoly].poly->rat.coeff[0], 
						final_polys[npoly].poly->rat.coeff[0], p2->exp1);
				mpz_set(n, p2->base2);
				mpz_pow_ui(n, n, p2->exp2);
				mpz_mul(final_polys[npoly].poly->rat.coeff[0], final_polys[npoly].poly->rat.coeff[0], n);

				mpz_gcd(n, final_polys[npoly].poly->rat.coeff[0], final_polys[npoly].poly->rat.coeff[1]);
				final_polys[npoly].difficulty -= log10(mpz_get_d(n))*deg;
				mpz_tdiv_q(final_polys[npoly].poly->rat.coeff[0], final_polys[npoly].poly->rat.coeff[0], n);
				mpz_tdiv_q(final_polys[npoly].poly->rat.coeff[1], final_polys[npoly].poly->rat.coeff[1], n);

				mpz_set(final_polys[npoly].poly->m, final_polys[npoly].poly->rat.coeff[0]);
				mpz_invert(n, final_polys[npoly].poly->rat.coeff[1], final_polys[npoly].n);
				mpz_mul(final_polys[npoly].poly->m, final_polys[npoly].poly->m, n);
				mpz_mod(final_polys[npoly].poly->m, final_polys[npoly].poly->m, final_polys[npoly].n);
				mpz_neg(final_polys[npoly].poly->rat.coeff[1], final_polys[npoly].poly->rat.coeff[1]);
			}

			final_polys[npoly].difficulty += log10(mpz_get_d(final_polys[npoly].c[deg]));

			// copy algebraic form
			mpz_set_si(final_polys[npoly].base1, x); final_polys[npoly].exp1 = y;	
			mpz_set_si(final_polys[npoly].base2, y); final_polys[npoly].exp2 = x;
			final_polys[npoly].form_type = SNFS_XYYXF;

			// check correctness and evaluate
			check_poly(&final_polys[npoly], fobj->VFLAG);
			approx_norms(&final_polys[npoly]);

			// be a little smarter about what degrees we consider
			if (final_polys[npoly].difficulty < 120)
			{
				if (deg == 6)
					final_polys[npoly].valid = 0;
				else
					avg_diff += final_polys[npoly].difficulty;
			}
			else if (final_polys[npoly].difficulty < 170)
			{
				if (deg == 6)
					final_polys[npoly].valid = 0;
				else
					avg_diff += final_polys[npoly].difficulty;
			}
			else if (final_polys[npoly].difficulty < 220)
			{
				if (deg == 4)
					final_polys[npoly].valid = 0;
				else
					avg_diff += final_polys[npoly].difficulty;
			}
			else
			{
				if (deg == 4)
					final_polys[npoly].valid = 0;
				else
					avg_diff += final_polys[npoly].difficulty;
			}

			if (final_polys[npoly].valid)
			{
				if (fobj->VFLAG > 0) print_snfs(&final_polys[npoly], stdout);
				npoly++;
			}
			else
			{	// being explicit
				snfs_clear(&final_polys[npoly]);
				snfs_init(&final_polys[npoly]);
			}
		}
	}

	// be a little smarter about what degrees we consider
	avg_diff /= (double)npoly;
	if (avg_diff < 120)
	{
		if (deg == 6)
			final_polys[npoly].valid = 0;
		fobj->nfs_obj.pref_degree = 4;
		fobj->nfs_obj.alt_degree = 5;
	}
	else if (avg_diff < 170)
	{
		if (deg == 6)
			final_polys[npoly].valid = 0;
		fobj->nfs_obj.pref_degree = 5;
		fobj->nfs_obj.alt_degree = 4;
	}
	else if (avg_diff < 220)
	{
		if (deg == 4)
			final_polys[npoly].valid = 0;
		fobj->nfs_obj.pref_degree = 5;
		fobj->nfs_obj.alt_degree = 6;
	}
	else
	{
		if (deg == 4)
			final_polys[npoly].valid = 0;
		fobj->nfs_obj.pref_degree = 6;
		fobj->nfs_obj.alt_degree = 5;
	}

	printf("generated %d polynomials: average difficulty = %1.2f, preferred degree = %d, alternate degree = %d\n",
		npoly, avg_diff, fobj->nfs_obj.pref_degree, fobj->nfs_obj.alt_degree);
    fflush(stdout);

    // clean up
    for (i = 0; i<alloc_base_poly; i++)
    {
        snfs_clear(&polys[i]);
    }
    free(polys);

    for (i = npoly; i < apoly; i++)
    {
        snfs_clear(&final_polys[i]);
    }

	*npolys = npoly;
	return final_polys;
}

// we've now measured the difficulty for poly's of all common degrees possibly formed
// in several different ways.  now we have a decision to make based largely on difficulty and 
// degree.  We want to pick low difficulty, but only if the degree allows the norms on 
// both sides to be approximately equal.  Sometimes multiple degrees satisfy this requirement
// approximately equally in which case only test-sieving can really resolve the difference.
// if the difficulty is below a threshold, just pick one, else, do some test sieving.

// not only can test sieving pick the best poly, but it can reveal optimizations
// to the job file, such as adjustments to the siever version or lpb values
// to get  rels/q into the desired range (somewhere around 2-4 rels/q).
// Therefore, pass in job structures to be modified rather than using throwaway
// objects during the test sieving.
nfs_job_t* snfs_test_sieve(fact_obj_t *fobj, snfs_t *polys, int npoly, nfs_job_t* jobs, int force_test)
{
	int i, dotest, minscore_id;	

	// only one poly - don't bother test sieving it :)
	if ((npoly < 2) && (!force_test))
		return &jobs[0];

	// see if any poly is big enough to justify test sieving
	for (i=0, dotest = 0; i<npoly; i++)
		if (polys[i].sdifficulty > fobj->nfs_obj.snfs_testsieve_threshold) 
			dotest = 1;

	// input too small - test sieving not justified
	if( !dotest )
		return &jobs[0];

	// we don't want canceled test sieves to impact the overall nfs job...
	IGNORE_NFS_ABORT = 1;
	minscore_id = test_sieve(fobj, jobs, npoly, 0);
	IGNORE_NFS_ABORT = 0;

	if( minscore_id < 0 )
	{
		printf("gen: warning: test sieving failed, reverting to top ranked poly");
		minscore_id = 0;
	}

	return &jobs[minscore_id];
}

void snfs_make_job_file(fact_obj_t *fobj, nfs_job_t *job)
{
	FILE *out;
	snfs_t *poly = job->snfs;

	out = fopen(fobj->nfs_obj.job_infile, "w");
	if (out == NULL)
	{
		printf("could not create %s for writing\n", fobj->nfs_obj.job_infile);
		exit(0);
	}

	// print the header stuff and the poly
	print_snfs(poly, out);

	fclose(out);

	return;
}

void snfs_scale_difficulty(snfs_t *polys, int npoly, int VFLAG)
{
	// poly degrees yielding very unbalanced norms are less desireable for sieving.
	// reflect this by scaling the difficulty of these polys.  Since special-q can
	// bring the norm down on one side or another (at the expense of the other side),
	// we add one to the difficulty for every X orders of magnitude imbalance beyond 
	// between the two sides (until we figure out something more accurate)
	int i;
    double magic = 5.0;

	for (i=0; i<npoly; i++)
	{
		double ratio, absa = fabs(polys[i].anorm), absr = fabs(polys[i].rnorm);

		// slight preference to sieve on the algebraic side...
		if ((log10(absr) - magic) > log10(absa))
		{
			ratio = absr / absa;
			polys[i].poly->side = RATIONAL_SPQ;
		}
		else
		{
			ratio = absa / absr;
			polys[i].poly->side = ALGEBRAIC_SPQ;			
		}

		if (VFLAG > 1)
			printf("gen: anorm: %1.2e, rnorm: %1.2e, ratio: %1.2e, log10(ratio) = %1.2f\n",
				polys[i].anorm, polys[i].rnorm, ratio, log10(ratio));
		ratio = log10(ratio) / magic;

		if (ratio < 0.)
			ratio = 0.;
		
		polys[i].sdifficulty = polys[i].difficulty + ratio;
	}
	
	return;
}

void analyze_one_poly_xface(snfs_t *poly)
{
	char line[1024], *ptr;
	FILE *in;
	msieve_obj *obj;
	char outfile[80];

	sprintf(outfile, "YAFU_get_poly_score.out");

	obj = msieve_obj_new("", MSIEVE_FLAG_USE_LOGFILE, NULL, outfile,
		NULL, 0, 0, (uint32_t)0, (enum cpu_type)9, 0, 0, 0, (uint32_t)0, "");

	remove(outfile);
	poly->poly->murphy = 1e-99;
	analyze_one_poly(obj, &poly->poly->rat, &poly->poly->alg, poly->poly->skew);
	in = fopen(outfile, "r");
	if (in != NULL)
	{
		while (!feof(in))
		{
			ptr = fgets(line, 1024, in);
			if (ptr == NULL)
				break;

			ptr = strstr(line, "size");
			if (ptr == NULL)
				continue;

			poly->poly->size = strtod(ptr + 5, NULL);

			ptr = strstr(ptr, "alpha");
			if (ptr == NULL)
				continue;

			poly->poly->alpha = strtod(ptr + 6, NULL);

			ptr = strstr(ptr, "combined");
			if (ptr == NULL)
				continue;

			poly->poly->murphy = strtod(ptr + 11, NULL);

			ptr = strstr(ptr, "rroots");
			if (ptr == NULL)
				continue;

			sscanf(ptr + 8, "%d", &poly->poly->rroots);
		}
		fclose(in);			
	}
	msieve_obj_free(obj);

	return;
}

int qcomp_snfs_murphy(const void *x, const void *y)
{
	snfs_t *xx = (snfs_t *)x;
	snfs_t *yy = (snfs_t *)y;

	// sort descending
	if ((*xx).poly->murphy > (*yy).poly->murphy)
		return -1;
	else if ((*xx).poly->murphy == (*yy).poly->murphy)
		return 0;
	else
		return 1; 
	
	// this sometimes doesn't work right... double/int conversion?
	//return ((snfs_t*)x)->sdifficulty - ((snfs_t*)y)->sdifficulty;
}

int qcomp_snfs_sdifficulty(const void *x, const void *y)
{
	snfs_t *xx = (snfs_t *)x;
	snfs_t *yy = (snfs_t *)y;

	if ((*xx).sdifficulty > (*yy).sdifficulty)
		return 1;
	else if ((*xx).sdifficulty == (*yy).sdifficulty)
		return 0;
	else
		return -1; 
	
	// this sometimes doesn't work right... double/int conversion?
	//return ((snfs_t*)x)->sdifficulty - ((snfs_t*)y)->sdifficulty;
}

int snfs_rank_polys(fact_obj_t *fobj, snfs_t *polys, int npoly)
{
	// rank by scaled difficulty
	int i, j, k;
	int pref_count = 0;
	int alt_count = 0;

	// first eliminate duplicate 
	for (i=0, k=0; i<npoly; i++)
	{
		for (j=i+1; j<npoly; j++)
		{
			if ((mpz_cmp(polys[i].poly->m, polys[j].poly->m) == 0) &&
				(polys[i].poly->alg.degree == polys[j].poly->alg.degree))
			{
				//polys[i].sdifficulty = 99999999.;
				polys[i].poly->murphy = 1e-99;
				k++;
				break;
			}
		}
	}
	if (fobj->VFLAG > 0 && k > 0)
		printf("nfs: rejected %d duplicate polys out of %d\n", k, npoly);

	// then sort
	//qsort(polys, npoly, sizeof(snfs_t), &qcomp_snfs_sdifficulty);
	qsort(polys, npoly, sizeof(snfs_t), &qcomp_snfs_murphy);

	// keep the first two polys of the preferred degree, and the top
	// one from a competing degree, that haven't already been rejected.
	for (i=0, j=0; i<npoly; i++)
	{
		if ((polys[i].poly->alg.degree == fobj->nfs_obj.pref_degree) &&
			polys[i].poly->murphy > 1e-99)
		{
			pref_count++;
			if (pref_count > 2)
				polys[i].poly->murphy = 1e-99;
		}
		else if (polys[i].poly->alg.degree == fobj->nfs_obj.alt_degree)
		{
			alt_count++;
			if (alt_count > 1)
				polys[i].poly->murphy = 1e-99;
		}
		else
			polys[i].poly->murphy = 1e-99;
	}

	// then sort again
	//qsort(polys, npoly, sizeof(snfs_t), &qcomp_snfs_sdifficulty);
	qsort(polys, npoly, sizeof(snfs_t), &qcomp_snfs_murphy);

	j = MIN(pref_count,2) + MIN(alt_count,1);
	return MIN(j,npoly);
}

int tdiv_int(int x, int *factors, uint64_t* primes, uint64_t num_p)
{
	int numf = 0;
	int xx = x;
	int i;

	i=0;
	while ((xx > 1) && (primes[i] < 1000))
	{
		int q = (int)primes[i];
		
		if (xx%q != 0)
			i++;
		else
		{			
			xx /= q;
			factors[numf++] = q;
		}
	}

	return numf;
}

int tdiv_mpz(mpz_t x, int *factors, uint64_t* primes, uint64_t num_p)
{
	int numf = 0;
	mpz_t xx;
	int i, r;

	mpz_init(xx);
	mpz_set(xx, x);

	i=0;
	while ((mpz_cmp_ui(xx,1) > 0) && (primes[i] < 1000))
	{
		int q = (int)primes[i];
		
		r = mpz_tdiv_ui(xx, q);
		
		if (r != 0)
			i++;
		else
		{			
			mpz_tdiv_q_ui(xx, xx, q);
			factors[numf++] = q;
		}
	}

	return numf;
}
