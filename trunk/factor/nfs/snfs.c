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

#include "nfs.h"
#include "util.h"
#include "gmp_xface.h"

void snfs_init(snfs_t *poly)
{	
	poly->alambda = 0.;
	poly->alim = 0;
	poly->base1 = 0;
	poly->base2 = 0;
	poly->c[0] = poly->c[1] = poly->c[2] = poly->c[3] = 0;
	poly->c[4] = poly->c[5] = poly->c[6] = poly->c[7] = poly->c[8] = 0;
	poly->coeff1 = 0;
	poly->coeff2 = 0;
	poly->degree = 0;
	poly->difficulty = 0.;
	poly->exp1 = 0;
	poly->exp2 = 0;
	poly->form_type = SNFS_NONE;
	poly->lpba = 0;
	poly->lpbr = 0;
	poly->mfba = 0;
	poly->mfbr = 0;
	poly->rlambda = 0;
	poly->rlim = 0;
	poly->skew = 0;
	poly->anorm = 0;
	poly->rnorm = 0;
	poly->sdifficulty = 0.;
	poly->side = RATIONAL_SPQ;
	mpz_init(poly->y0);
	mpz_init(poly->y1);
	mpz_init(poly->m);
	mpz_init(poly->n);

	return;
}

void snfs_clear(snfs_t *poly)
{
	mpz_clear(poly->y0);
	mpz_clear(poly->y1);
	mpz_clear(poly->m);
	mpz_clear(poly->n);
	poly->form_type = SNFS_NONE;

	return;
}

void snfs_copy_poly(snfs_t *src, snfs_t *dest)
{
	dest->alambda = src->alambda;
	dest->alim = src->alim;
	dest->base1 = src->base1;
	dest->base2 = src->base2;
	dest->c[0] = src->c[0];
	dest->c[1] = src->c[1];
	dest->c[2] = src->c[2];
	dest->c[3] = src->c[3];
	dest->c[4] = src->c[4];
	dest->c[5] = src->c[5];
	dest->c[6] = src->c[6];
	dest->c[7] = src->c[7];
	dest->c[8] = src->c[8];
	dest->coeff1 = src->coeff1;
	dest->coeff2 = src->coeff2;
	dest->degree = src->degree;
	dest->difficulty = src->difficulty;
	dest->exp1 = src->exp1;
	dest->exp2 = src->exp2;
	dest->form_type = src->form_type;
	dest->lpba = src->lpba;
	dest->lpbr = src->lpbr;
	dest->mfba = src->mfba;
	dest->mfbr = src->mfbr;
	dest->rlambda = src->rlambda;
	dest->rlim = src->rlim;
	dest->skew = src->skew;
	dest->anorm = src->anorm;
	dest->rnorm = src->rnorm;
	dest->side = src->side;
	dest->sdifficulty = src->sdifficulty;
	mpz_set(dest->y0, src->y0);
	mpz_set(dest->y1, src->y1);
	mpz_set(dest->n, src->n);
	mpz_set(dest->m, src->m);

	return;
}

void check_poly(snfs_t *poly)
{
	// make sure the generated poly is correct.  For instance if
	// the coefficents overflow and int the poly will be invalid, which
	// is fine since we wouldn't want to use coefficients that large anyway.
	// check that each polynomial mod n at m is zero
	mpz_t t;
	int i;

	mpz_init(t);

	poly->valid = 1;
	mpz_set_ui(t, 0);
	for (i = poly->degree; i >= 0; i--)
	{
		mpz_mul(t, t, poly->m);
		if (poly->c[i] < 0)
			mpz_sub_ui(t, t, abs(poly->c[i]));
		else
			mpz_add_ui(t, t, abs(poly->c[i]));
		mpz_mod(t, t, poly->n);
	}
	if (mpz_cmp_ui(t,0) != 0)
	{
		poly->valid = 0;
		//gmp_fprintf (stderr, "Error: M=%Zd is not a root of f(x) % N\n", poly->m);
		//gmp_fprintf (stderr, "n = %Zd\n", poly->n);
		//fprintf (stderr, "f(x) = ");
		//for (i = poly->degree; i >= 0; i--)
		//	gmp_fprintf (stderr, "%d*x^%d %c %d", poly->c < 0 ? '-' : '+', abs(poly->c[i]), i);
		//gmp_fprintf (stderr, "\n""Remainder is %Zd\n", t);
	}

	mpz_mul(t, poly->y1, poly->m);
	mpz_add(t, t, poly->y0);
	mpz_mod(t, t, poly->n);
	if (mpz_cmp_ui(t,0) != 0)
	{
		poly->valid = 0;
		//gmp_fprintf (stderr, "n = %Zd\n", poly->n);
		//gmp_fprintf (stderr, "Error: M=%Zd is not a root of g(x) % N\n", poly->m);
		//gmp_fprintf (stderr, "Remainder is %Zd\n", t);
	}

	return;
}

void print_poly(snfs_t *poly, FILE *out)
{
	// print the poly to stdout
	char c;
	int i;
	char side[80];

	if (poly->coeff2 < 0)
		c = '-';
	else
		c = '+';

	if (poly->side == RATIONAL_SPQ)
		sprintf(side, "rational");
	else
		sprintf(side, "algebraic");

	if (poly->form_type == SNFS_H_CUNNINGHAM)
	{
		fprintf(out, "# %d^%d%c%d^%d, difficulty: %1.2f, anorm: %1.2e, rnorm: %1.2e\n", 
			poly->base1, poly->exp1, c, poly->base2, poly->exp1, poly->difficulty,
			poly->anorm, poly->rnorm);
	}
	else
	{
		if (poly->coeff1 == 1)
			fprintf(out, "# %d^%d%c%d, difficulty: %1.2f, anorm: %1.2e, rnorm: %1.2e\n", 
				poly->base1, poly->exp1, c, abs(poly->coeff2), poly->difficulty,
				poly->anorm, poly->rnorm);
		else
			fprintf(out, "# %d*%d^%d%c%d, difficulty: %1.2f, anorm: %1.2e, rnorm: %1.2e\n", 
				abs(poly->coeff1), poly->base1, poly->exp1, c, abs(poly->coeff2), poly->difficulty,
				poly->anorm, poly->rnorm);
	}
	if (poly->sdifficulty > 0)
		fprintf(out, "# scaled difficulty: %1.2f, suggest sieving %s side\n", poly->sdifficulty, side);
	gmp_fprintf(out, "n: %Zd\n", poly->n);
	fprintf(out, "type: snfs\nsize: %d\n", (int)poly->sdifficulty);
	fprintf(out, "skew: %1.2f\n", poly->skew);
	for (i=8; i>=0; i--)
		if (poly->c[i] != 0) fprintf(out, "c%d: %d\n", i, poly->c[i]);
	gmp_fprintf(out, "Y1: %Zd\n", poly->y1);
	gmp_fprintf(out, "Y0: %Zd\n", poly->y0);
	gmp_fprintf(out, "m: %Zd\n\n", poly->m);

	return;
}

void approx_norms(snfs_t *poly)
{
	// anorm ~= b^d * f(a/b) where f = the algebraic poly
	// rnorm ~= b * g(a/b) where g is the linear poly
	// a,b depend on the siever used.  according to the rsa768 paper, 
	// a,b were respectively about 3e9*sqrt(skew), 3e9/sqrt(skew).
	// here we use 1e6 instead of 3e9... it might not matter so much as
	// long as it is consistent between a/rnorm
	int i;
	double a,b;

	a = sqrt(poly->skew) * 1000000.;
	b = 1000000. / (sqrt(poly->skew));

	poly->anorm = 0;
	for (i=8; i>=0; i--)
		poly->anorm += abs(poly->c[i])*pow(a/b, i);
	poly->anorm *= pow(b, poly->degree);

	poly->rnorm = (fabs(mpz_get_d(poly->y1))*a/b + mpz_get_d(poly->y0)) * b;

	return;
}

void find_brent_form(fact_obj_t *fobj, snfs_t *form)
{
	int i,j,maxa,maxb;
	mpz_t p, a, b, r, n;
	uint32 inc = 1<<30;

	// brent numbers take the form a^n +/- 1, where 13<=a<=99 and the product is less than 10^255.
	// this routine finds inputs of the brent form as well as cunningham (a <= 12), 
	// and oddperfect forms (a > 99, n > 1) and considers inputs up to 1000 bits in size (~ 10^302).
	// once we have the form, we can create a polynomial for it for snfs processing.

	maxa = 100;
	maxb = 100;

	mpz_init(p);
	mpz_init(a);
	mpz_init(b);
	mpz_init(r);
	mpz_init(n);

	mpz_set(n, fobj->nfs_obj.gmp_n);

	for (i=2; i<maxa; i++)
	{
		// ignore prime power bases... because e.g. 9^n == 3^(2*n) and thus
		// will have already been tested.
		if ((i==4) || (i==9) || (i==16) || (i==25) || (i==36) || (i==49) ||
			(i==64) || (i==81) || (i==8) || (i==27) || (i==32))
			continue;

		mpz_set_ui(b, i);
		mpz_pow_ui(p, b, 31);

		// limit the exponent so that the number is less than 1000 bits
		maxb = 1000 / log((double)i) + 1;

		if (VFLAG > 1)
			printf("nfs: checking %d^x +/- 1 for 20 <= x <= %d\n", i, maxb);

		for (j=32; j<maxb; j++)
		{
			mpz_mul(p, p, b);		// p = i^j
			mpz_add_ui(r, n, inc);	// r = n + 2^30
			mpz_mod(r, r, p);		// r = (n + 2^30) % i^j

			// now, if r is a single limb, then the input has a small coefficient
			if (mpz_sizeinbase(r,2) <= 32)
			{
				// get the second coefficient first
				int sign, c1, c2 = mpz_get_ui(r);
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
				
				if (VFLAG > 0) 
				{
					if (c1 > 1)
					{
						if (sign == POSITIVE)
							printf("nfs: input divides %d*%d^%d + %d\n", c1, i, j, c2);
						else
							printf("nfs: input divides %d*%d^%d - %d\n", c1, i, j, c2);
					}
					else
					{
						if (sign == POSITIVE)
							printf("nfs: input divides %d^%d + %d\n", i, j, c2);
						else
							printf("nfs: input divides %d^%d - %d\n", i, j, c2);
					}
				}

				form->form_type = SNFS_BRENT;
				form->coeff1 = c1;
				form->base1 = i;
				form->base2 = 1;
				form->exp1 = j;
				form->coeff2 = sign ? -c2 : c2;
				mpz_set(form->n, n);
				gen_brent_poly(fobj, form);
				goto done;
			}
		}
	}

	for (i = maxb; i>1; i--)
	{
		// now that we've reduced the exponent considerably, check for 
		// large bases by looking at the remaining possible exponents.
		if (VFLAG > 0)
			printf("nfs: checking x^%d +/- 1\n", i);

		// check -1 case:
		mpz_add_ui(a, n, 1);
		mpz_root(b, a, i);
		mpz_pow_ui(p, b, i);
		if ((mpz_cmp(p, a) == 0) && (mpz_sizeinbase(b, 2) < 32))
		{
			if (VFLAG > 0) printf("nfs: input divides %d^%d - 1\n", (int)mpz_get_ui(b), i);
			form->form_type = SNFS_BRENT;
			form->base1 = mpz_get_ui(b);
			form->exp1 = i;
			form->coeff1 = -1;
			mpz_set(form->n, n);
			gen_brent_poly(fobj, form);
			goto done;
		}

		// check -1 case:
		mpz_sub_ui(a, n, 1);
		mpz_root(b, a, i);
		mpz_pow_ui(p, b, i);
		if ((mpz_cmp(p, a) == 0) && (mpz_sizeinbase(b, 2) < 32))
		{
			if (VFLAG > 0) printf("nfs: input divides %d^%d + 1\n", (int)mpz_get_ui(b), i);
			form->form_type = SNFS_BRENT;
			form->base1 = mpz_get_ui(b);
			form->exp1 = i;
			form->coeff1 = 1;
			mpz_set(form->n, n);
			gen_brent_poly(fobj, form);
			goto done;
		}

	}


done:
	mpz_clear(p);
	mpz_clear(a);
	mpz_clear(b);
	mpz_clear(r);
	return;
}

void find_hcunn_form(fact_obj_t *fobj, snfs_t *form)
{
	int i,j,k,maxa,kmax;
	mpz_t pa, pb, a, b, r, g, n;

	// homogeneous cunninghams take the form a^n +/- b^n, where a,b <= 12 and gcd(a,b) == 1.
	// this routine finds inputs of the hcunninghams form, and considers
	// inputs up to 1000 bits in size (~ 10^302).
	// once we have the form, we can create a polynomial for it for snfs processing.

	maxa = 13;

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

			// limit the exponent so that the number is less than 1000 bits
			kmax = 1000 / log((double)i) + 1;
			if (VFLAG > 1)
				printf("nfs: checking %d^x +/- %d^x for 20 <= x <= %d\n", i, j, kmax);

			for (k=20; k<kmax; k++)
			{
				mpz_mul(pa, pa, a);
				mpz_mul(pb, pb, b);

				mpz_add(g, pa, pb);
				mpz_mod(r, g, n);
				if (mpz_cmp_ui(r, 0) == 0)
				{
					if (VFLAG > 0) printf("nfs: input divides %d^%d + %d^%d\n", i, k, j, k);
					form->form_type = SNFS_H_CUNNINGHAM;
					form->base1 = i;
					form->base2 = j;
					form->exp1 = k;
					form->coeff1 = 1;
					mpz_set(form->n, n);
					gen_brent_poly(fobj, form);
					goto done;
				}

				mpz_sub(g, pa, pb);
				mpz_mod(r, g, n);
				if (mpz_cmp_ui(r, 0) == 0)
				{
					if (VFLAG > 0) printf("nfs: input divides %d^%d - %d^%d\n", i, k, j, k);
					form->form_type = SNFS_H_CUNNINGHAM;
					form->base1 = i;
					form->base2 = j;
					form->exp1 = k;
					form->coeff1 = -1;
					mpz_set(form->n, n);
					gen_brent_poly(fobj, form);
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
	return;
}

void find_xyyxf_form(fact_obj_t *fobj, snfs_t *form)
{
	int i,j,k,maxa,kmax;
	mpz_t pa, pb, a, b, r, g, n;

	//form->form_type = SNFS_NONE;
	// xyyxf numbers take the form x^y + y^x, where 1<y<x<151
	// this routine finds inputs of the xyyxf form, and considers
	// inputs up to 1000 bits in size (~ 10^302).
	// once we have the form, we can create a polynomial for it for snfs processing.

	maxa = 13;

	// TODO: finish this routine!

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

			// limit the exponent so that the number is less than 1000 bits
			kmax = 1000 / log((double)i) + 1;
			if (VFLAG > 0)
				printf("nfs: checking %d^x +/- %d^x for 20 <= x <= %d\n", i, j, kmax);

			for (k=20; k<kmax; k++)
			{
				mpz_mul(pa, pa, a);
				mpz_mul(pb, pb, b);

				mpz_add(g, pa, pb);
				mpz_mod(r, g, n);
				if (mpz_cmp_ui(r, 0) == 0)
				{
					if (VFLAG > 0) printf("nfs: input divides %d^%d + %d^%d\n", i, k, j, k);
					form->form_type = SNFS_H_CUNNINGHAM;
					form->base1 = i;
					form->base2 = j;
					form->exp1 = k;
					form->coeff1 = 1;
					goto done;
				}

				mpz_sub(g, pa, pb);
				mpz_mod(r, g, n);
				if (mpz_cmp_ui(r, 0) == 0)
				{
					if (VFLAG > 0) printf("nfs: input divides %d^%d - %d^%d\n", i, k, j, k);
					form->form_type = SNFS_H_CUNNINGHAM;
					form->base1 = i;
					form->base2 = j;
					form->exp1 = k;
					form->coeff1 = -1;
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
	return;
}

void gen_brent_poly(fact_obj_t *fobj, snfs_t *poly)
{
	int i, me;
	int e = poly->exp1;
	int b = poly->base1;
	int b2 = poly->base2;
	mpz_t n, m;
	double d, skew, k;
	int f[100];
	int numf = 0;
	snfs_t *polys;
	int npoly = 0;
	int apoly;	

	mpz_init(n);
	mpz_init(m);

	// test for algebraic factors - we look for several specific but common forms:
	// 21*k, 15*k, 13*k, 11*k, 7*k, 5*k, 3*k
	// if any of these are available we immediately use them, because dividing out an algebraic factor
	// will always be lower difficulty then playing with exponents only, even if the degree
	// is sub-optimal.  
	if (poly->exp1 % 15 == 0 && (poly->coeff1 == 1))
	{
		polys = (snfs_t *)malloc(sizeof(snfs_t));
		snfs_init(polys);
		npoly = 1;
		snfs_copy_poly(poly, polys);		// copy algebraic form

		// a^(15k) +/- 1 has an algebraic factor which is an 8th degree symmetric polynomial.
		// check for this case before a^(3k) or a^(5k) since the algbraic reduction is greater.
		// the 8th degree poly can be halved as in the 11*k and 13*k cases, resulting in 
		// the following quartic:
		polys->degree = 4;
		k = poly->exp1 / 15;
		polys->c[4] = 1;
		polys->c[3] = poly->coeff2;
		polys->c[2] = -4;
		polys->c[1] = -poly->coeff2 * 4;
		polys->c[0] = 1;
		mpz_set_ui(m, poly->base1);
		polys->difficulty = log10(mpz_get_d(m)) * 8. * k;
		mpz_pow_ui(polys->m, m, k);
		polys->skew = 1.;

		if (poly->form_type == SNFS_H_CUNNINGHAM)
		{
			// multiplying through by b = a^k gives g(x) = bx - (b^2 + 1).  but
			// for homogeneous polys a = a1/a2 so we have g(x) = (a1/a2)^k * x - ((a1^2/a2^2)^k + 1) = 0
			// multiplying through by a2^2 gives:
			// g(x) = (a1*a2)^k*x - (a1^(2k) + a2^(2k)) with m = (a1/a2)^k
			mpz_set_si(polys->y1, b2*b);
			mpz_pow_ui(polys->y1, polys->y1, k);
			mpz_mul(polys->y0, polys->m, polys->m);
			mpz_set_ui(n, poly->base2);
			mpz_pow_ui(n, n, 2*k);
			mpz_add(polys->y0, polys->y0, n);
			mpz_sqrt(n, n);
			mpz_invert(n, n, poly->n);
			mpz_mul(polys->m, polys->m, n);
			mpz_mod(polys->m, polys->m, poly->n);
			mpz_mul_si(polys->y1, polys->y1, -1);
		}
		else
		{
			// for halved degree polynomials, Y1 becomes -x^k and Y0 becomes x^(2k) + 1
			// such that y1*x + y0 evaluated at M = x^k + x^-k is 0.
			mpz_neg(polys->y1, polys->m);
			mpz_mul(polys->y0, polys->m, polys->m);
			mpz_add_ui(polys->y0, polys->y0, 1);
			mpz_invert(m, polys->m, polys->n);
			mpz_add(polys->m, polys->m, m);
		}
		
		check_poly(polys);
		approx_norms(polys);
	}
	else if (poly->exp1 % 21 == 0 && (poly->coeff1 == 1))
	{
		polys = (snfs_t *)malloc(sizeof(snfs_t));
		snfs_init(polys);
		npoly = 1;
		snfs_copy_poly(poly, polys);		// copy algebraic form

		// a^(21k) +/- 1 has an algebraic factor which is a 12th degree symmetric polynomial.
		// check for this case before a^(3k) or a^(7k) since the algbraic reduction is greater.
		// the 12th degree poly can be halved as in the 11*k and 13*k cases, resulting in 
		// the following sextic:
		polys->degree = 6;
		k = poly->exp1 / 21;
		polys->c[6] = 1;
		polys->c[5] = poly->coeff2;
		polys->c[4] = -6;
		polys->c[3] = -poly->coeff2 * 6;
		polys->c[2] = 8;
		polys->c[1] = poly->coeff2 * 8;
		polys->c[0] = 1;
		mpz_set_ui(m, poly->base1);
		polys->difficulty = log10(mpz_get_d(m)) * 12. * k;
		mpz_pow_ui(polys->m, m, k);
		polys->skew = 1.;
		
		if (poly->form_type == SNFS_H_CUNNINGHAM)
		{
			// multiplying through by b = a^k gives g(x) = bx - (b^2 + 1).  but
			// for homogeneous polys a = a1/a2 so we have g(x) = (a1/a2)^k * x - ((a1^2/a2^2)^k + 1) = 0
			// multiplying through by a2^2 gives:
			// g(x) = (a1*a2)^k*x - (a1^(2k) + a2^(2k)) with m = (a1/a2)^k
			mpz_set_si(polys->y1, b2*b);
			mpz_pow_ui(polys->y1, polys->y1, k);
			mpz_mul(polys->y0, polys->m, polys->m);
			mpz_set_ui(n, poly->base2);
			mpz_pow_ui(n, n, 2*k);
			mpz_add(polys->y0, polys->y0, n);
			mpz_sqrt(n, n);
			mpz_invert(n, n, poly->n);
			mpz_mul(polys->m, polys->m, n);
			mpz_mod(polys->m, polys->m, poly->n);
			mpz_mul_si(polys->y1, polys->y1, -1);
		}
		else
		{
			// for halved degree polynomials, Y1 becomes -x^k and Y0 becomes x^(2k) + 1
			// such that y1*x + y0 evaluated at M = x^k + x^-k is 0.
			mpz_neg(polys->y1, polys->m);
			mpz_mul(polys->y0, polys->m, polys->m);
			mpz_add_ui(polys->y0, polys->y0, 1);
			mpz_invert(m, polys->m, polys->n);
			mpz_add(polys->m, polys->m, m);
		}

		check_poly(polys);
		approx_norms(polys);
	}
	else if (poly->exp1 % 6 == 0 && (poly->coeff1 == 1))
	{
		polys = (snfs_t *)malloc(sizeof(snfs_t));
		snfs_init(polys);
		npoly = 1;
		snfs_copy_poly(poly, polys);		// copy algebraic form

		// a^(3k) +/- 1, k even, is divisible by (a^k +/- 1) giving a quadratic in a^k.
		// the quadratic can be converted into a quartic...
		// see: http://www.mersennewiki.org/index.php/SNFS_Polynomial_Selection
		// todo: look into making degree 6 based on difficulty
		polys->degree = 4;
		k = poly->exp1 / 6;
		polys->c[4] = 1;
		polys->c[2] = -1;
		polys->c[0] = 1;
		mpz_set_ui(m, poly->base1);
		polys->difficulty = log10(mpz_get_d(m)) * 4. * k;
		mpz_pow_ui(polys->m, m, k);
		polys->skew = 1.;

		if (poly->form_type == SNFS_H_CUNNINGHAM)
		{
			mpz_set_si(polys->y1, b2);
			mpz_pow_ui(polys->y1, polys->y1, k);
			mpz_set(polys->y0, polys->m);
			mpz_invert(n, polys->y1, poly->n);
			mpz_mul(polys->m, polys->m, n);
			mpz_mod(polys->m, polys->m, poly->n);
			mpz_mul_si(polys->y1, polys->y1, -1);
		}
		else
		{
			// Y1 = -1, Y0 = m such that y1*m + y0 = 0
			mpz_set(polys->y0, polys->m);
			mpz_set_si(polys->y1, -1);
		}

		check_poly(polys);
		approx_norms(polys);
	}
	else if (poly->exp1 % 6 == 3 && (poly->coeff1 == 1))
	{
		polys = (snfs_t *)malloc(sizeof(snfs_t));
		snfs_init(polys);
		npoly = 1;
		snfs_copy_poly(poly, polys);		// copy algebraic form

		// a^(3k) +/- 1, k odd, is divisible by (a^k +/- 1) giving a quadratic in a^k.
		// the quadratic can be converted into a quartic...
		// see: http://www.mersennewiki.org/index.php/SNFS_Polynomial_Selection
		// todo: look into making degree 6 based on difficulty
		polys->degree = 4;
		k = (poly->exp1 - 3) / 6;
		polys->c[4] = poly->base1 * poly->base1;
		polys->c[2] = poly->base1 * -poly->coeff2;
		polys->c[0] = 1;
		mpz_set_ui(m, poly->base1);
		polys->skew = pow(mpz_get_d(m), -0.5);
		polys->difficulty = log10(mpz_get_d(m)) * 4. * k;
		mpz_pow_ui(polys->m, m, k);		

		if (poly->form_type == SNFS_H_CUNNINGHAM)
		{
			mpz_set_si(polys->y1, b2);
			mpz_pow_ui(polys->y1, polys->y1, k);
			mpz_set(polys->y0, polys->m);
			mpz_invert(n, polys->y1, poly->n);
			mpz_mul(polys->m, polys->m, n);
			mpz_mod(polys->m, polys->m, poly->n);
			mpz_mul_si(polys->y1, polys->y1, -1);
		}
		else
		{
			// Y1 = -1, Y0 = m such that y1*m + y0 = 0
			mpz_set(polys->y0, polys->m);
			mpz_set_si(polys->y1, -1);
		}

		check_poly(polys);
		approx_norms(polys);
	}
	else if (poly->exp1 % 5 == 0 && (poly->coeff1 == 1))
	{
		polys = (snfs_t *)malloc(sizeof(snfs_t));
		snfs_init(polys);
		npoly = 1;
		snfs_copy_poly(poly, polys);		// copy algebraic form

		// a^(5k) +/- 1 is divisible by (a^k +/- 1) giving a quartic in a^k
		polys->degree = 4;
		k = poly->exp1 / 5;
		polys->c[4] = 1;
		polys->c[3] = -poly->coeff2;
		polys->c[2] = 1;
		polys->c[1] = -poly->coeff2;
		polys->c[0] = 1;
		mpz_set_ui(m, poly->base1);
		polys->difficulty = log10(mpz_get_d(m)) * 4. * k;
		mpz_pow_ui(polys->m, m, k);
		polys->skew = 1.;		

		if (poly->form_type == SNFS_H_CUNNINGHAM)
		{
			mpz_set_si(polys->y1, b2);
			mpz_pow_ui(polys->y1, polys->y1, k);
			mpz_set(polys->y0, polys->m);
			mpz_invert(n, polys->y1, poly->n);
			mpz_mul(polys->m, polys->m, n);
			mpz_mod(polys->m, polys->m, poly->n);
			mpz_mul_si(polys->y1, polys->y1, -1);
		}
		else
		{
			// Y1 = -1, Y0 = m such that y1*m + y0 = 0
			mpz_set(polys->y0, polys->m);
			mpz_set_si(polys->y1, -1);
		}

		mpz_set(polys->n, poly->n);
		check_poly(polys);
		approx_norms(polys);
	}
	else if (poly->exp1 % 7 == 0 && (poly->coeff1 == 1))
	{
		polys = (snfs_t *)malloc(sizeof(snfs_t));
		snfs_init(polys);
		npoly = 1;
		snfs_copy_poly(poly, polys);		// copy algebraic form

		// a^(7k) +/- 1 is divisible by (a^k +/- 1) giving a sextic in a^k
		polys->degree = 6;
		k = poly->exp1 / 7;
		polys->c[6] = 1;
		polys->c[5] = -poly->coeff2;
		polys->c[4] = 1;
		polys->c[3] = -poly->coeff2;
		polys->c[2] = 1;
		polys->c[1] = -poly->coeff2;
		polys->c[0] = 1;
		mpz_set_ui(m, poly->base1);
		polys->difficulty = log10(mpz_get_d(m)) * 6. * k;
		mpz_pow_ui(polys->m, m, k);
		polys->skew = 1.;		

		if (poly->form_type == SNFS_H_CUNNINGHAM)
		{
			mpz_set_si(polys->y1, b2);
			mpz_pow_ui(polys->y1, polys->y1, k);
			mpz_set(polys->y0, polys->m);
			mpz_invert(n, polys->y1, poly->n);
			mpz_mul(polys->m, polys->m, n);
			mpz_mod(polys->m, polys->m, poly->n);
			mpz_mul_si(polys->y1, polys->y1, -1);
		}
		else
		{
			// Y1 = -1, Y0 = m such that y1*m + y0 = 0
			mpz_set(polys->y0, polys->m);
			mpz_set_si(polys->y1, -1);
		}

		check_poly(polys);
		approx_norms(polys);
	}
	else if (poly->exp1 % 11 == 0 && (poly->coeff1 == 1))
	{
		polys = (snfs_t *)malloc(sizeof(snfs_t));
		snfs_init(polys);
		npoly = 1;
		snfs_copy_poly(poly, polys);		// copy algebraic form

		// a^(11k) +/- 1 is divisible by (a^k +/- 1) giving a poly in a^k of degree 10.
		// but the poly is symmetric and can be halved to degree 5... 
		// see http://www.mersennewiki.org/index.php/SNFS_Polynomial_Selection
		polys->degree = 5;
		k = poly->exp1 / 11;
		polys->c[5] = 1;
		polys->c[4] = -poly->coeff2*1;
		polys->c[3] = -4;
		polys->c[2] = poly->coeff2*3;
		polys->c[1] = 3;
		polys->c[0] = -poly->coeff2;
		mpz_set_ui(m, poly->base1);
		polys->difficulty = log10(mpz_get_d(m)) * 10. * k;
		mpz_pow_ui(polys->m, m, k);
		polys->skew = 1.;	
				
		if (poly->form_type == SNFS_H_CUNNINGHAM)
		{
			// multiplying through by b = a^k gives g(x) = bx - (b^2 + 1).  but
			// for homogeneous polys a = a1/a2 so we have g(x) = (a1/a2)^k * x - ((a1^2/a2^2)^k + 1) = 0
			// multiplying through by a2^2 gives:
			// g(x) = (a1*a2)^k*x - (a1^(2k) + a2^(2k)) with m = (a1/a2)^k
			mpz_set_si(polys->y1, b2*b);
			mpz_pow_ui(polys->y1, polys->y1, k);
			mpz_mul(polys->y0, polys->m, polys->m);
			mpz_set_ui(n, poly->base2);
			mpz_pow_ui(n, n, 2*k);
			mpz_add(polys->y0, polys->y0, n);
			mpz_sqrt(n, n);
			mpz_invert(n, n, poly->n);
			mpz_mul(polys->m, polys->m, n);
			mpz_mod(polys->m, polys->m, poly->n);
			mpz_mul_si(polys->y1, polys->y1, -1);
		}
		else
		{
			// for halved degree polynomials, Y1 becomes -x^k and Y0 becomes x^(2k) + 1
			// such that y1*x + y0 evaluated at M = x^k + x^-k is 0.
			mpz_neg(polys->y1, polys->m);
			mpz_mul(polys->y0, polys->m, polys->m);
			mpz_add_ui(polys->y0, polys->y0, 1);
			mpz_invert(m, polys->m, polys->n);
			mpz_add(polys->m, polys->m, m);
		}

		check_poly(polys);
		approx_norms(polys);
	}
	else if (poly->exp1 % 13 == 0 && (poly->coeff1 == 1))
	{
		polys = (snfs_t *)malloc(sizeof(snfs_t));
		snfs_init(polys);
		npoly = 1;
		snfs_copy_poly(poly, polys);		// copy algebraic form

		// a^(13k) +/- 1 is divisible by (a^k +/- 1) giving a poly in a^k of degree 12.
		// but the poly is symmetric and can be halved to degree 6... 
		// see http://www.mersennewiki.org/index.php/SNFS_Polynomial_Selection
		polys->degree = 6;
		k = poly->exp1 / 13;
		polys->c[6] = 1;
		polys->c[5] = -poly->coeff2;
		polys->c[4] = -5;
		polys->c[3] = poly->coeff2*4;
		polys->c[2] = 6;
		polys->c[1] = -poly->coeff2*3;
		polys->c[0] = -1;
		mpz_set_ui(m, poly->base1);
		polys->difficulty = log10(mpz_get_d(m)) * 12. * k;
		mpz_pow_ui(polys->m, m, k);
		polys->skew = 1.;	
		
		if (poly->form_type == SNFS_H_CUNNINGHAM)
		{
			// multiplying through by b = a^k gives g(x) = bx - (b^2 + 1).  but
			// for homogeneous polys a = a1/a2 so we have g(x) = (a1/a2)^k * x - ((a1^2/a2^2)^k + 1) = 0
			// multiplying through by a2^2 gives:
			// g(x) = (a1*a2)^k*x - (a1^(2k) + a2^(2k)) with m = (a1/a2)^k
			mpz_set_si(polys->y1, b2*b);
			mpz_pow_ui(polys->y1, polys->y1, k);
			mpz_mul(polys->y0, polys->m, polys->m);
			mpz_set_ui(n, poly->base2);
			mpz_pow_ui(n, n, 2*k);
			mpz_add(polys->y0, polys->y0, n);
			mpz_sqrt(n, n);
			mpz_invert(n, n, poly->n);
			mpz_mul(polys->m, polys->m, n);
			mpz_mod(polys->m, polys->m, poly->n);
			mpz_mul_si(polys->y1, polys->y1, -1);
		}
		else
		{
			// for halved degree polynomials, Y1 becomes -x^k and Y0 becomes x^(2k) + 1
			// such that y1*x + y0 evaluated at M = x^k + x^-k is 0.
			mpz_neg(polys->y1, polys->m);
			mpz_mul(polys->y0, polys->m, polys->m);
			mpz_add_ui(polys->y0, polys->y0, 1);
			mpz_invert(m, polys->m, polys->n);
			mpz_add(polys->m, polys->m, m);
		}

		check_poly(polys);
		approx_norms(polys);
	}
	else 
	{

		// No algebraic factor - play with powers and composite bases
		mpz_set_ui(m, b);
		if (!mpz_probab_prime_p(m,10))
		{
			int bb = b;
			// factor the base
			i=0;
			while ((bb > 1) && (spSOEprimes[i] < 100))
			{
				int q = (int)spSOEprimes[i];
		
				if (bb%q != 0)
					i++;
				else
				{			
					bb /= q;
					f[numf++] = q;
				}
			}
		}

		// initialize candidate polynomials now that we know how many we'll need.
		// each factor of the base generates one, plus 2 for
		// the whole base raised and lowered, for each degree
		if (numf > 1) apoly = (numf + 2) * 3;
		else apoly = 6;

		polys = (snfs_t *)malloc(apoly * sizeof(snfs_t));
		for (i=0; i<apoly; i++)
			snfs_init(&polys[i]);

		if (VFLAG > 0)
		{
			printf("gen: ========================================================\n");
			printf("gen: considering the following polynomials:\n");
			printf("gen: ========================================================\n");
		}

		for (i=4; i<7; i++)
		{
			int c0, cd;

			if (e % i == 0)
			{
				// the degree divides the exponent - resulting polynomial is straightforward
				me = e / i;
				mpz_set_si(m, b);
				mpz_pow_ui(m, m, me);
				d = mpz_get_d(m);
				d = log10(d) * (double)i;
				skew = 1.0;
				cd = poly->coeff1; c0 = poly->coeff2;
				snfs_copy_poly(poly, &polys[npoly]);		// copy algebraic form
				polys[npoly].difficulty = d;
				polys[npoly].skew = skew;
				polys[npoly].c[i] = cd;
				polys[npoly].c[0] = c0;
				if (poly->form_type == SNFS_H_CUNNINGHAM)
				{
					mpz_set_si(polys[npoly].y1, b2);
					mpz_pow_ui(polys[npoly].y1, polys[npoly].y1, me);
					mpz_set(polys[npoly].y0, m);
					mpz_set(polys[npoly].m, m);
					mpz_invert(n, polys[npoly].y1, polys[npoly].n);
					mpz_mul(polys[npoly].m, polys[npoly].m, n);
					mpz_mod(polys[npoly].m, polys[npoly].m, polys[npoly].n);
					mpz_mul_si(polys[npoly].y1, polys[npoly].y1, -1);
				}
				else
				{
					mpz_set_si(polys[npoly].y1, -1);
					mpz_set(polys[npoly].y0, m);
					mpz_set(polys[npoly].m, m);
				}
				polys[npoly].degree = i;			
				check_poly(&polys[npoly]);
				approx_norms(&polys[npoly]);
				if (VFLAG > 0 && polys[npoly].valid) print_poly(&polys[npoly], stdout);
				if (polys[npoly].valid) npoly++;
			}
			else
			{
				// degree does not divide the exponent, try increasing the exponent
				int inc = (i - e % i);
				me = (e+inc) / i;
				mpz_set_si(m, b);
				mpz_pow_ui(m, m, me);
				d = mpz_get_d(m);
				d = log10(d) * (double)i;
				cd = (int)pow((double)b2, inc) * poly->coeff1;
				c0 = (int)pow((double)b,inc) * poly->coeff2;
				skew = pow((double)abs(c0)/(double)cd, 1./(double)i);
				snfs_copy_poly(poly, &polys[npoly]);		// copy algebraic form
				polys[npoly].difficulty = d;
				polys[npoly].skew = skew;
				polys[npoly].c[i] = cd;
				polys[npoly].c[0] = c0;
				if (poly->form_type == SNFS_H_CUNNINGHAM)
				{
					mpz_set_si(polys[npoly].y1, b2);
					mpz_pow_ui(polys[npoly].y1, polys[npoly].y1, me);
					mpz_set(polys[npoly].y0, m);
					mpz_set(polys[npoly].m, m);
					mpz_invert(n, polys[npoly].y1, polys[npoly].n);
					mpz_mul(polys[npoly].m, polys[npoly].m, n);
					mpz_mod(polys[npoly].m, polys[npoly].m, polys[npoly].n);
					mpz_mul_si(polys[npoly].y1, polys[npoly].y1, -1);
				}
				else
				{
					mpz_set_si(polys[npoly].y1, -1);
					mpz_set(polys[npoly].y0, m);
					mpz_set(polys[npoly].m, m);
				}
				polys[npoly].degree = i;			
				check_poly(&polys[npoly]);
				approx_norms(&polys[npoly]);
				if (VFLAG > 0 && polys[npoly].valid) print_poly(&polys[npoly], stdout);
				if (polys[npoly].valid) npoly++;

				// and decreasing the exponent
				inc = e % i;
				me = (e-inc) / i;
				mpz_set_si(m, b);
				mpz_pow_ui(m, m, me);
				d = mpz_get_d(m);
				d = log10(d) * (double)i + log10(pow((double)b,inc));
				cd = (int)pow((double)b,inc) * poly->coeff1;
				c0 = (int)pow((double)b2, inc) * poly->coeff2;
				skew = pow((double)abs(c0)/(double)cd, 1./(double)i);
				// leading coefficient contributes to the difficulty
				d += log10((double)cd);
				snfs_copy_poly(poly, &polys[npoly]);		// copy algebraic form
				polys[npoly].difficulty = d;
				polys[npoly].skew = skew;
				polys[npoly].c[i] = cd;
				polys[npoly].c[0] = c0;
				if (poly->form_type == SNFS_H_CUNNINGHAM)
				{
					mpz_set_si(polys[npoly].y1, b2);
					mpz_pow_ui(polys[npoly].y1, polys[npoly].y1, me);
					mpz_set(polys[npoly].y0, m);
					mpz_set(polys[npoly].m, m);
					mpz_invert(n, polys[npoly].y1, polys[npoly].n);
					mpz_mul(polys[npoly].m, polys[npoly].m, n);
					mpz_mod(polys[npoly].m, polys[npoly].m, polys[npoly].n);
					mpz_mul_si(polys[npoly].y1, polys[npoly].y1, -1);
				}
				else
				{
					mpz_set_si(polys[npoly].y1, -1);
					mpz_set(polys[npoly].y0, m);
					mpz_set(polys[npoly].m, m);
				}
				polys[npoly].degree = i;
				check_poly(&polys[npoly]);
				approx_norms(&polys[npoly]);
				if (VFLAG > 0 && polys[npoly].valid) print_poly(&polys[npoly], stdout);
				if (polys[npoly].valid) npoly++;

				// and playing with composite bases
				if (numf > 1)
				{
					int j;

					// multiply by powers of each factor individually to get that 
					// factor to a multiple of degree
					for (j=0; j<numf; j++)
					{
						int k, i1, i2, bb;
						// move it up
						i1 = (i - e % i);
						c0 = pow((double)f[j], i1) * poly->coeff2;
						// since we moved the current factor up, move the other factors down
						// otherwise this would be the same as moving the whole composite base up.
						// also move the second term up...
						cd = pow((double)b2, i1) * poly->coeff1;
						bb = 1;
						i2 = e % i;
						for (k=0; k<numf; k++)
						{
							if (k == j) continue;
							cd *= pow((double)f[k], i2);
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
						// leading coefficient contributes to the difficulty
						d += log10((double)cd);
						skew = pow((double)abs(c0)/(double)cd, 1./(double)i);
						snfs_copy_poly(poly, &polys[npoly]);		// copy algebraic form
						polys[npoly].difficulty = d;
						polys[npoly].skew = skew;
						polys[npoly].c[i] = cd;
						polys[npoly].c[0] = c0;
						if (poly->form_type == SNFS_H_CUNNINGHAM)
						{
							mpz_set_si(polys[npoly].y1, b2);
							mpz_pow_ui(polys[npoly].y1, polys[npoly].y1, (e+i1) / i);
							mpz_set(polys[npoly].y0, m);
							mpz_set(polys[npoly].m, m);
							mpz_invert(n, polys[npoly].y1, polys[npoly].n);
							mpz_mul(polys[npoly].m, polys[npoly].m, n);
							mpz_mod(polys[npoly].m, polys[npoly].m, polys[npoly].n);
							mpz_mul_si(polys[npoly].y1, polys[npoly].y1, -1);
						}
						else
						{
							mpz_set_si(polys[npoly].y1, -1);
							mpz_set(polys[npoly].y0, m);
							mpz_set(polys[npoly].m, m);
						}
						polys[npoly].degree = i;
						check_poly(&polys[npoly]);
						approx_norms(&polys[npoly]);
						if (VFLAG > 0 && polys[npoly].valid) print_poly(&polys[npoly], stdout);
						if (polys[npoly].valid) npoly++;
					}
				}			
			}
		}
	}

	// we've now measured the difficulty for poly's of all common degrees possibly formed
	// in several different ways.  now we have a decision to make based largely on difficulty and 
	// degree.  We want to pick low difficulty, but only if the degree allows the norms on 
	// both sides to be approximately equal.  Sometimes multiple degrees satisfy this requirement
	// approximately equally in which case only test-sieving can really resolve the difference.
	// if the difficulty is below a threshold, just pick one, else, do some test sieving.
	snfs_scale_difficulty(polys, npoly);
	snfs_rank_polys(polys, npoly);

	if (VFLAG > 0)
	{
		printf("gen: ========================================================\n");
		printf("gen: filtered polynomials:\n");
		printf("gen: ========================================================\n");

		for (i=0; i<npoly; i++)
			print_poly(&polys[i], stdout);
	}

	snfs_test_sieve(fobj, polys, npoly);

	if (VFLAG > 0)
	{
		printf("gen: ========================================================\n");
		printf("gen: selected polynomial:\n");
		printf("gen: ========================================================\n");

		print_poly(&polys[0], stdout);
	}

	// output the polynomial to the default nfs job file
	if (npoly > 0)
	{
		FILE *out;

		out = fopen(fobj->nfs_obj.job_infile, "w");
		print_poly(&polys[0], out);
		fclose(out);
	}

	mpz_clear(n);
	mpz_clear(m);

	return;
}

void snfs_test_sieve(fact_obj_t *fobj, snfs_t *polys, int npoly)
{
	int i, dotest;
	uint32 count;
	int minscore_id;
	double score[3];
	double min_score = 999999999.;
	struct timeval stop, stop2;	// stop time of this job
	struct timeval start, start2;	// start time of this job
	TIME_DIFF *	difference;
	double t_time;

	// only one poly - don't bother test sieving it :)
	if (npoly == 1)
		return;	

	// see if any poly within the top three polys for this input is 
	// big enough to justify test sieving
	for (i=0, dotest = 0; i<npoly && i<3; i++)
		if (polys[i].sdifficulty > fobj->nfs_obj.snfs_testsieve_threshold) dotest = 1;

	// input too small - test sieving not justified
	if (dotest == 0)
		return;

	// check to make sure we can find ggnfs sievers
	if (check_for_sievers(fobj, 0) == 1)
	{
		printf("gen: can't find ggnfs lattice sievers - aborting test sieving\n");
		return;
	}

	// else, for at most the top 3 polys above the cutoff, sieve a small range
	// and compute sec/rel.
	minscore_id = 0;
	gettimeofday(&start2, NULL);
	for (i=0; i<3 && i<npoly; i++)
	{		
		ggnfs_job_t job;
		uint32 fill_params;
		char syscmd[1024], tmpbuf[1024];
		FILE *in;
		char side;

		//make a job file for each input
		// initialize some job parameters
		job.current_rels = 0;
		job.lpb = 0;
		job.mfb = 0;
		job.lambda = 0;
		job.type = 1;		// 0==GNFS, 1==SNFS
		job.fblim = 0;
		job.size = polys[i].sdifficulty;
		get_ggnfs_params(fobj, &job);
		sprintf(fobj->nfs_obj.job_infile, "test-%d.poly", i);
		snfs_make_poly_file(fobj, &polys[i]);
		fill_params = PARAM_FLAG_FBLIM | PARAM_FLAG_LPB | PARAM_FLAG_MFB | PARAM_FLAG_LAMBDA;
		fill_job_file(fobj, &job, fill_params);

		//create the afb/rfb - we don't want the time it takes to do this to
		//pollute the sieve timings
		sprintf(syscmd,"%s -b %s -k -c 0 -F", job.sievername, fobj->nfs_obj.job_infile);

		printf("gen: commencing construction of afb\n");
		gettimeofday(&start, NULL);
		system(syscmd);
		gettimeofday(&stop, NULL);
		difference = my_difftime (&start, &stop);
		t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
		free(difference);			
		printf("gen: afb generation took %6.4f seconds.\n", t_time);
		sprintf(tmpbuf, "%s.afb.0", fobj->nfs_obj.job_infile);
		remove(tmpbuf);
		MySleep(.1);

		//start the test
		job.startq = job.fblim / 2;
		side = polys[i].side == RATIONAL_SPQ ? 'r' : 'a';
		sprintf(syscmd,"%s -%c %s -f %u -c %u -o %s.out",
			job.sievername, side, fobj->nfs_obj.job_infile, job.startq, 5000, fobj->nfs_obj.job_infile);
		printf("gen: commencing test sieving of polynomial %d on side '%c' over range %u-%u\n", i, 
			side, job.startq, job.startq + 5000);
		gettimeofday(&start, NULL);
		system(syscmd);
		gettimeofday(&stop, NULL);
		difference = my_difftime (&start, &stop);
		t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
		free(difference);			

		//count relations
		sprintf(tmpbuf, "%s.out", fobj->nfs_obj.job_infile);
		in = fopen(tmpbuf,"r");
		if (in != NULL)
		{
			count = 0;
			while (fgets(tmpbuf, GSTR_MAXSIZE, in) != NULL)
				count++;
			fclose(in);
		}
		else
			count = 1;	//no divide by zero

		score[i] = t_time / (double)count;
		if (score[i] < min_score)
		{
			minscore_id = i;
			min_score = score[i];
			printf("gen: new best score of %1.6f sec/rel, estimated total sieving time = %1.2f\n", 
				score[i], score[i] * job.min_rels);
		}
		else
			printf("gen: score was %1.6f sec/rel, estimated total sieving time = %1.2f\n", 
				score[i], score[i] * job.min_rels);
	}
	gettimeofday(&stop2, NULL);
	difference = my_difftime (&start2, &stop2);
	t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
	free(difference);			
	printf("gen: test sieving took %1.2f seconds\n", t_time);

	// copy the winner into the top slot of the poly array (it will be picked
	// by the calling routine)
	snfs_copy_poly(&polys[minscore_id], &polys[0]);

	return;
}

void snfs_make_poly_file(fact_obj_t *fobj, snfs_t *poly)
{
	FILE *out;

	out = fopen(fobj->nfs_obj.job_infile, "w");
	if (out == NULL)
	{
		printf("could not create %s for writing\n", fobj->nfs_obj.job_infile);
		exit(0);
	}

	print_poly(poly, out);
	fclose(out);

	return;
}

void snfs_scale_difficulty(snfs_t *polys, int npoly)
{
	// poly degrees yielding very unbalanced norms are less desireable for sieving.
	// reflect this by scaling the difficulty of these polys.  Since special-q can
	// bring the norm down on one side or another (at the expense of the other side),
	// we add one to the difficulty for every order of magnitude imbalance beyond 6
	// between the two sides (until we figure out something more accurate)
	int i;

	for (i=0; i<npoly; i++)
	{
		double ratio;
		
		if (polys[i].anorm > polys[i].rnorm)
		{
			ratio = polys[i].anorm / polys[i].rnorm;
			polys[i].side = ALGEBRAIC_SPQ;
		}
		else
		{
			ratio = polys[i].rnorm / polys[i].anorm;
			polys[i].side = RATIONAL_SPQ;
		}

		printf("anorm: %1.2e, rnorm: %1.2e, ratio: %1.2e, log10(ratio) = %1.2f\n",
			polys[i].anorm, polys[i].rnorm, ratio, log10(ratio));
		ratio = log10(ratio) - 6.;

		if (ratio < 0.)
			ratio = 0.;

		polys[i].sdifficulty = polys[i].difficulty + ratio;		
	}

	return;
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
}

void snfs_rank_polys(snfs_t *polys, int npoly)
{
	// rank by scaled difficulty
	int i;

	qsort(polys, npoly, sizeof(snfs_t), &qcomp_snfs_sdifficulty);

	for (i=0; i<npoly; i++)
		polys[i].rank = i;

	return;
}


