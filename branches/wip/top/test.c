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
#include "util.h"
#include "qs.h"
#include "factor.h"
#include "monty.h"
#include <ecm.h>

double spAcc(int m);
double sqr_acc(int m, int sz);
double subadd_acc(int m, int sz);
double shortmuldiv_acc(int m, int sz);
double muldiv_acc(int m, int sz);
double sqrt_acc(int m, int sz);
double gcd_acc(int m, int sz);
void richard_guy_problem_e7(void);
int SQUFOF_alpertron(int64 N, int64 queue[]);
int squfof_rds(int64 n, int *fac1, int *fac2);
void enqu(int q, int *iter);


int qqueue[100];
int qpoint;

void enqu(int q, int *iter)
{
    qqueue[qpoint] = q;
    if (++qpoint > 100) *iter = -1;
}


int squfof_rds(int64 n, int *fac1, int *fac2)
{   /* start of squfof: factor n as fac1*fac2  faster in FP?????*/
    int64 temp, temp1;
    register int iq, ll, l2, p, pnext, q, qlast, r, s, t, i;
    int jter, iter;

    qlast = 1;
    s = (int)sqrt(n);

    p = s;
    temp1 = s * s;
    temp = n - temp1;                 /* temp = n - floor(sqrt(n))^2   */
    if (temp == 0)
    {                                   /* Here n is a square            */
        *fac1 = s;
        *fac2 = s;
        return(1);
    }

    q = (int)temp;              /* q = excess of n over next smaller square */
    ll = 1 + 2 * (int)sqrt((double)(p + p));
    l2 = ll / 2;
    qpoint = 0;

    /*   In the loop below, we need to check if q is a square right before   */
    /*  the end of the loop.  Is there a faster way? The current way is      */
    /*   EXPENSIVE! (many branches and double prec sqrt)                     */

    for (jter = 0; jter < 800000; jter++)      /* I see no way to speed this   */
    {                                     /*  main loop                   */
        iq = (s + p) / q;
        pnext = iq*q - p;
        if (q <= ll)
        {
            if ((q & 1) == 0) enqu(q / 2, &jter);
            else if (q <= l2) enqu(q, &jter);
            if (jter < 0)
            {
                return 0;
            }
        }
        t = qlast + iq*(p - pnext);
        qlast = q;
        q = t;
        p = pnext;                          /* check for square; even iter   */
        if (jter & 1) continue;             /* jter is odd:omit square test  */
        r = (int)sqrt((double)q);                 /* r = floor(sqrt(q))      */
        if (q != r*r) continue;
        if (qpoint == 0) goto gotit;
        for (i = 0; i<qpoint - 1; i += 2)      /* treat queue as list for simplicity*/
        {
            if (r == qqueue[i]) goto contin;
            if (r == qqueue[i + 1]) goto contin;
        }
        if (r == qqueue[qpoint - 1]) continue;
        goto gotit;
    contin:;
    }   /* end of main loop */

gotit:;
    qlast = r;
    p = p + r*((s - p) / r);
    temp = (int64)p * (int64)p;
    temp = n - temp;
    temp1 = temp / qlast;
    q = (int)temp1;					/* q = (n - p*p)/qlast (div is exact)*/
    for (iter = 0; iter<40000; iter++)
    {                              /* begin second main loop            */
        iq = (s + p) / q;                /* unroll it, of course              */
        pnext = iq*q - p;
        if (p == pnext) goto gotfac;
        t = qlast + iq*(p - pnext);
        qlast = q;
        q = t;
        p = pnext;
        iq = (s + p) / q;
        pnext = iq*q - p;
        if (p == pnext) goto gotfac;
        t = qlast + iq*(p - pnext);
        qlast = q;
        q = t;
        p = pnext;
        iq = (s + p) / q;
        pnext = iq*q - p;
        if (p == pnext) goto gotfac;
        t = qlast + iq*(p - pnext);
        qlast = q;
        q = t;
        p = pnext;
        iq = (s + p) / q;
        pnext = iq*q - p;
        if (p == pnext) goto gotfac;
        t = qlast + iq*(p - pnext);
        qlast = q;
        q = t;
        p = pnext;
    }


    return(0);                               /* this shouldn't happen      */

gotfac:; if ((q & 1) == 0) q /= 2;      /* q was factor or 2*factor   */
    *fac1 = q;
    temp = n / q;
    *fac2 = (int)temp;
    return(1);
}


/*Implementation of algorithm explained in Gower and Wagstaff paper */
int SQUFOF_alpertron(int64 N, int64 queue[])
{
    double sqrtn;
    int B, Q, Q1, P, P1, L, S;
    int i, r, s, t, q, u;
    int queueHead, queueTail, queueIndex;

    /* Step 1: Initialize */
    if ((N & 3) == 1)
    {
        N <<= 1;
    }
    sqrtn = sqrt(N);
    S = (int)sqrtn;
    if ((long)(S + 1)*(long)(S + 1) <= N)
    {
        S++;
    }
    if ((long)S*(long)S > N)
    {
        S--;
    }
    if ((long)S*(long)S == N)
    {
        return S;
    }
    Q1 = 1;
    P = S;
    Q = (int)N - P*P;
    L = (int)(2 * sqrt(2 * sqrtn));
    B = L << 1;
    queueHead = 0;
    queueTail = 0;

    /* Step 2: Cycle forward to find a proper square form */
    for (i = 0; i <= B; i++)
    {
        q = (S + P) / Q;
        P1 = q*Q - P;
        if (Q <= L)
        {
            if ((Q & 1) == 0)
            {
                queue[queueHead++] = Q >> 1;
                queue[queueHead++] = P % (Q >> 1);
                if (queueHead == 100)
                {
                    queueHead = 0;
                }
            }
            else if (Q + Q <= L)
            {
                queue[queueHead++] = Q;
                queue[queueHead++] = P % Q;
                if (queueHead == 100)
                {
                    queueHead = 0;
                }
            }
        }

        t = Q1 + q*(P - P1);
        Q1 = Q;
        Q = t;
        P = P1;
        if ((i & 1) == 0 && ((Q & 7) < 2 || (Q & 7) == 4))
        {
            r = (int)sqrt(Q);
            if (r*r == Q)
            {
                queueIndex = queueTail;
                for (;;)
                {
                    if (queueIndex == queueHead)
                    {
                        /* Step 3: Compute inverse square root of the square form */
                        Q1 = r;
                        u = (S - P) % r;
                        u += (u >> 31) & r;
                        P = S - u;
                        Q = (int)((N - (long)P*(long)P) / Q1);
                        /* Step 4: Cycle in the reverse direction to find a factor of N */
                        for (;;)
                        {
                            q = (S + P) / Q;
                            P1 = q*Q - P;
                            if (P == P1)
                            {
                                break;
                            }
                            t = Q1 + q*(P - P1);
                            Q1 = Q;
                            Q = t;
                            P = P1;
                        }

                        /* Step 5: Get the factor of N */
                        if ((Q & 1) == 0)
                        {
                            return Q >> 1;
                        }
                        return Q;
                    }
                    s = queue[queueIndex++];
                    t = queue[queueIndex++];
                    if (queueIndex == 100)
                    {
                        queueIndex = 0;
                    }
                    if ((P - t) % s == 0)
                    {
                        break;
                    }
                }
                if (r > 1)
                {
                    queueTail = queueIndex;
                }
                if (r == 1)
                {
                    queueIndex = queueTail;
                    for (;;)
                    {
                        if (queueIndex == queueHead)
                        {
                            break;
                        }
                        if (queue[queueIndex] == 1)
                        {
                            return 0;
                        }
                        queueIndex += 2;
                        if (queueIndex == 100)
                        {
                            queueIndex = 0;
                        }
                    }
                }
            }
        }
    }
    return 0;
}


void test_dlp_composites_par()
{
    FILE *in;
    uint64 *comp, f64;
    uint32 *f1;
    uint64 *f;
    uint32 *f2, totBits, minBits, maxBits;
    double t_time;
    clock_t start, stop;
    int i, j, num, correct, num_successes;
    z tmp, tmp2, t1, t2, t3;
    mpz_t gmptmp;
    int64 queue[100];
    struct timeval gstart;
    struct timeval gstop;
    int nf;
    int num_files;
    char filenames[15][80];

    mpz_init(gmptmp);
    comp = (uint64 *)malloc(2000000 * sizeof(uint64));
    f1 = (uint32 *)malloc(2000000 * sizeof(uint32));
    f2 = (uint32 *)malloc(2000000 * sizeof(uint32));
    f = (uint64 *)malloc(2000000 * sizeof(uint64));
    zInit(&tmp);
    zInit(&tmp2);
    zInit(&t2);
    zInit(&t1);
    zInit(&t3);

    i = 0;
    strcpy(filenames[i++], "pseudoprimes_37bit.dat");
    strcpy(filenames[i++], "pseudoprimes_40bit.dat");
    strcpy(filenames[i++], "pseudoprimes_41bit.dat");
    strcpy(filenames[i++], "pseudoprimes_48bit.dat");
    strcpy(filenames[i++], "pseudoprimes_51bit.dat");
    strcpy(filenames[i++], "pseudoprimes_55bit.dat");
    strcpy(filenames[i++], "pseudoprimes_59bit.dat");
    strcpy(filenames[i++], "pseudoprimes_60bit.dat");
    strcpy(filenames[i++], "pseudoprimes_61bit.dat");
    strcpy(filenames[i++], "pseudoprimes_62bit.dat");
    strcpy(filenames[i++], "pseudoprimes_63bit.dat");
    strcpy(filenames[i++], "pseudoprimes_64bit.dat");
    num_files = 11;

    for (nf = 0; nf < num_files; nf++)
    {
        in = fopen(filenames[nf], "r");

        start = clock();
        i = 0;
        totBits = 0;
        minBits = 999;
        maxBits = 0;
        //read in everything
        while (!feof(in))
        {
            fscanf(in, "%" PRIu64 ",%u,%u", comp + i, f1 + i, f2 + i);
            sp642z(comp[i], &tmp);
            j = zBits(&tmp);
            totBits += j;
            if ((uint32)j > maxBits)
                maxBits = j;
            if ((uint32)j < minBits && j != 0)
                minBits = j;
            i++;
        }
        num = i;
        num = 100000;
        fclose(in);
        stop = clock();
        t_time = (double)(stop - start) / (double)CLOCKS_PER_SEC;
        printf("data read in %2.4f sec\n", t_time);
        printf("average bits of input numbers = %.2f\n", (double)totBits / (double)i);
        printf("minimum bits of input numbers = %d\n", minBits);
        printf("maximum bits of input numbers = %d\n", maxBits);

        gettimeofday(&gstart, NULL);

        num_successes = par_shanks_loop(comp, f, num);

        gettimeofday(&gstop, NULL);
        t_time = my_difftime(&gstart, &gstop);

        printf("par_squfof reported %d successes\n", num_successes);
        printf("average squfof time per input = %1.4f ms\n", 1000 * t_time / (double)num);

        correct = 0;
        for (i = 0; i < num; i++)
        {
            if (((uint32)f[i] == f1[i]) || ((uint32)f[i] == f2[i]))
            {
                correct++;
            }
            else
            {
                int j;
                for (j = 0; j < 3; j++)
                {
                    f[i] = spbrent(comp[i], j + 1, 8192);
                    if ((f[i] == f1[i]) || (f[i] == f2[i]))
                    {
                        correct++;
                        break;
                    }
                }
            }
        }
        
        gettimeofday(&gstop, NULL);
        t_time = my_difftime(&gstart, &gstop);

        printf("overall got %d of %d correct in %2.2f sec\n", correct, num, t_time);
        printf("squfof percent correct = %.2f\n", 100.0*(double)num_successes / (double)num);
        printf("average time per input = %1.4f ms\n", 1000 * t_time / (double)num);
        
    }

    free(f);
    free(f1);
    free(f2);
    free(comp);
    zFree(&tmp);
    zFree(&tmp2);
    zFree(&t2);
    zFree(&t1);
    zFree(&t3);
    mpz_clear(gmptmp);
    return;
}


void gcd128(uint64 * u, uint64 * v, uint64 * w)
{
	// binary GCD for non-zero inputs.
	int shift;
	uint64 aa[4], bb[4];
	z a, b;

	a.val = aa;
	b.val = bb;
	a.val[0] = u[0];
	b.val[0] = v[0];
	a.val[1] = u[1];
	b.val[1] = v[1];
	a.size = 2;
	b.size = 2;
	a.alloc = 4;
	b.alloc = 4;

	/* Let shift := lg K, where K is the greatest power of 2
	dividing both u and v. */
	for (shift = 0; ((a.val[0] | b.val[0]) & 1) == 0; ++shift) {
		zShiftRight_1(&a, &a);
		zShiftRight_1(&b, &b);
	}
	
	while ((a.val[0] & 1) == 0)
	{
		zShiftRight_1(&a, &a);
	}

	/* From here on, u is always odd. */
	do {
		/* remove all factors of 2 in v -- they are not common */
		/*   note: v is not zero, so while will terminate */
		while ((b.val[0] & 1) == 0)  /* Loop X */
			zShiftRight_1(&b, &b);

		/* Now u and v are both odd. Swap if necessary so u <= v,
		then set v = v - u (which is even). For bignums, the
		swapping is just pointer movement, and the subtraction
		can be done in-place. */
		if (zCompare(&a, &b) > 0) {
			//uint64 t = v; v = u; u = t;
			uint64 *tmp = a.val;
			a.val = b.val;
			b.val = tmp;
		}  // Swap u and v.
		//v = v - u;                       // Here v >= u.
		zSub(&b, &a, &b);
	} while (zCompare(&b, &zZero) > 0);

	/* restore common factors of 2 */
	zShiftLeft(&a, &a, shift);
	w[0] = a.val[0];
	w[1] = a.val[1];

	return;
}

int brent128(monty128_t *mdata, uint64 * n, uint64 * f, uint32 a, uint32 imax)
{
	/*
	run pollard's rho algorithm on n with Brent's modification,
	returning the first factor found in f, or else 0 for failure.
	use f(x) = x^2 + c
	see, for example, bressoud's book.
	*/

	uint64 x[2];
	uint64 y[2];
	uint64 c[2];
	uint64 q[2];
	uint64 g[2];
	uint64 ys[2];
	uint64 t1[2];
	uint64 s[3];
	uint32 i = 0, k, r, m;
	int it;

	// starting state of algorithm.  
	r = 1;
	m = 16;
	i = 0;
	it = 0;

	q[0] = 1;
	y[0] = 0;
	g[0] = 1;
	c[0] = a;

	q[1] = 0;
	y[1] = 0;
	g[1] = 0;
	c[1] = 0;

	monty128_init(mdata, n);
	to_monty128(mdata, c);

	do
	{
		x[1] = y[1]; x[0] = y[0];
		for (i = 0; i <= r; i++)
		{
			sqrmod128(y, y, mdata);
			addmod128(y, c, y, mdata->n);
		}

		k = 0;
		do
		{
			ys[1] = y[1]; ys[0] = y[0];
			for (i = 1; i <= MIN(m, r - k); i++)
			{
				sqrmod128(y, y, mdata);
				addmod128(y, c, y, mdata->n);
				submod128(y, x, t1, mdata->n);
				mulmod128(t1, q, q, mdata);
			}
			//mpz_gcd(g, q, mdata->n);
			gcd128(q, mdata->n, g);
			k += m;
			it++;

			if (it > imax)
			{
				f[1] = 0; f[0] = 0;
				goto free;
			}

		} while (k < r && (g[0] == 1));
		r *= 2;
	} while (g[0] == 1);

	if ((g[1] == n[1]) && (g[0] == n[0]))
	{
		// back track
		do
		{
			sqrmod128(ys, ys, mdata);
			addmod128(ys, c, ys, mdata->n);
			submod128(ys, x, t1, mdata->n);
			gcd128(t1, mdata->n, g);
			//mpz_gcd(g, t1, mdata->n);
		} while ((g[1] == 0) && (g[0] == 1));

		if ((g[1] == n[1]) && (g[0] == n[0]))
		{
			f[1] = 0; f[0] = 0;
			goto free;
		}
		else
		{
			f[1] = g[1]; f[0] = g[0];
			goto free;
		}
	}
	else
	{
		f[1] = g[1]; f[0] = g[0];
		goto free;
	}

free:

	return it;
}

void test_dlp_composites()
{
	FILE *in;
	uint64 *comp, f64;
	uint32 *f1;
	uint32 *f2, bits, totBits, minBits, maxBits;
	double t_time;
	clock_t start, stop;
	int i, j, k, num, correct;
	z tmp, tmp2, t1, t2, t3;
	mpz_t gmptmp;
	int64 queue[100];
	struct timeval gstart;
	struct timeval gstop;
	int nf;
	int num_files;
	char filenames[20][80];


	fact_obj_t *fobj2;

	uint64 testLehman[29] = {
		5640012124823LL,
		7336014366011LL,
		19699548984827LL,
		52199161732031LL,
		73891306919159LL,
		112454098638991LL,

		32427229648727LL,
		87008511088033LL,
		92295512906873LL,
		338719143795073LL,
		346425669865991LL,
		1058244082458461LL,
		1773019201473077LL,
		6150742154616377LL,

		44843649362329LL,
		67954151927287LL,
		134170056884573LL,
		198589283218993LL,
		737091621253457LL,
		1112268234497993LL,
		2986396307326613LL,

		26275638086419LL,
		62246008190941LL,
		209195243701823LL,
		290236682491211LL,
		485069046631849LL,
		1239671094365611LL,
		2815471543494793LL,
		5682546780292609LL
	};


	if (0)
	{
		mpz_t gmp_comp, gmp_f, gmp_f2;
		char buf[1024];
		ecm_params my_ecm_params;
		tiny_qs_params *cosiqs = init_tinyqs();
		monty128_t mdata;
		int curves;
		int B1;
		int natt = 0;
		correct = 0;

		ecm_init(my_ecm_params);
		num = 0;
		minBits = 9999;
		maxBits = 0;
		totBits = 0;

		in = fopen("tlp_attempts_c95.txt", "r");

		mpz_init(gmp_comp);
		mpz_init(gmp_f);
		mpz_init(gmp_f2);

		goto cosiqs_start;
		goto tinyecm_start;

		gettimeofday(&gstart, NULL);

		while (!feof(in))
		{
			fgets(buf, 1024, in);
			gmp_sscanf(buf, "%Zd", gmp_comp);
			num++;

			B1 = 300;
			curves = 16;

			if (mpz_probab_prime_p(gmp_comp, 1))
				continue;

			natt++;
			bits = mpz_sizeinbase(gmp_comp, 2);
			if (bits < minBits)
				minBits = bits;
			if (bits > maxBits)
				maxBits = bits;
			totBits += bits;

			for (k = 0; k < curves; k++)
			{
				uint64 sigma = spRand(100, 1000000000);
				my_ecm_params->B1done = 1.0 + floor(1 * 128.) / 134217728.;
				//mpz_set_ui(my_ecm_params->B2, 0);
				mpz_set_ui(my_ecm_params->x, (unsigned long)0);
				mpz_set_ui(my_ecm_params->sigma, sigma);
				my_ecm_params->method = ECM_ECM;
				ecm_factor(gmp_f, gmp_comp, B1, my_ecm_params);
				if ((mpz_cmp_ui(gmp_f, 1) > 0) &&
				    (mpz_cmp(gmp_f, gmp_comp) < 0))
				{
					correct++;
					break;
				}
			}
		}

		fclose(in);

		gettimeofday(&gstop, NULL);
		t_time = my_difftime(&gstart, &gstop);

		printf("ecm found %d factors in %d inputs in %1.4f seconds\n", correct, num, t_time);
		printf("attempts had an average of %1.2f bits, %d min, %d max\n", 
			(double)totBits / natt, minBits, maxBits);
		printf("average time per input = %1.4f ms\n", 1000 * t_time / (double)num);

cosiqs_start:

		in = fopen("tlp_attempts_c95.txt", "r");

		gettimeofday(&gstart, NULL);
		correct = 0;
		num = 0;

		while (!feof(in))
		{
			fgets(buf, 1024, in);
			gmp_sscanf(buf, "%Zd", gmp_comp);
			num++;

			if (mpz_probab_prime_p(gmp_comp, 1))
				continue;

			tinyqs(cosiqs, gmp_comp, gmp_f, gmp_f2);
			if ((mpz_cmp_ui(gmp_f, 1) > 0) &&
				(mpz_cmp(gmp_f, gmp_comp) < 0))
			{
				correct++;
			}
		}

		fclose(in);

		gettimeofday(&gstop, NULL);
		t_time = my_difftime(&gstart, &gstop);

		printf("cosiqs found %d factors in %d inputs in %1.4f seconds\n", correct, num, t_time);
		printf("average time per input = %1.4f ms\n", 1000 * t_time / (double)num);

tinyecm_start:
		in = fopen("tlp_attempts_c95.txt", "r");

		gettimeofday(&gstart, NULL);
		correct = 0;
		num = 0;

		while (!feof(in))
		{
			fgets(buf, 1024, in);
			gmp_sscanf(buf, "%Zd", gmp_comp);
			num++;

			if (mpz_probab_prime_p(gmp_comp, 1))
				continue;

			tinyecm(gmp_comp, gmp_f, 200, 10000, 8, 0);

			if ((mpz_cmp_ui(gmp_f, 1) > 0) &&
				(mpz_cmp(gmp_f, gmp_comp) < 0) &&
				(mpz_tdiv_ui(gmp_comp, mpz_get_ui(gmp_f)) == 0))
			{
				correct++;
			}
		}

		fclose(in);

		gettimeofday(&gstop, NULL);
		t_time = my_difftime(&gstart, &gstop);

		printf("tinyecm found %d factors in %d inputs in %1.4f seconds\n", correct, num, t_time);
		printf("average time per input = %1.4f ms\n", 1000 * t_time / (double)num);

		mpz_clear(gmp_comp);
		mpz_clear(gmp_f);
		mpz_clear(gmp_f2);
	}

	


	fobj2 = (fact_obj_t *)malloc(sizeof(fact_obj_t));
	init_factobj(fobj2);


	mpz_init(gmptmp);
	comp = (uint64 *)malloc(2000000 * sizeof(uint64));
	f1 = (uint32 *)malloc(2000000 * sizeof(uint32));
	f2 = (uint32 *)malloc(2000000 * sizeof(uint32));
	zInit(&tmp);
	zInit(&tmp2);
	zInit(&t2);
	zInit(&t1);
	zInit(&t3);

	goto brent_marker;
	goto tinyecm_marker;

	
	for (i = 0; i < 29; i++)
	{
		f1[i] = spbrent(testLehman[i], 1, 8192);
		f2[i] = testLehman[i] / f1[i];
	}

	correct = 0;
	k = 0;
	for (i = 0; i < 29; i++)
	{
		f64 = LehmanFactor(testLehman[i], 3.5, 0, 0.1);
		printf("input %lu returned factor %lu, actual factors %u and %u (%lu)\n", 
			testLehman[i], f64, f1[i], f2[i], (uint64)f1[i] * (uint64)f2[i]);
		if ((f64 == f1[i]) || (f64 == f2[i]))
		{
			correct++;
		}
	}

	printf("Lehman got %d of %d correct of testvector inputs\n", correct, 29);

	//exit(0);
brent_marker:

    i = 0;
	strcpy(filenames[i++], "pseudoprimes_32bit.dat");
	strcpy(filenames[i++], "pseudoprimes_34bit.dat");
	strcpy(filenames[i++], "pseudoprimes_36bit.dat");
	strcpy(filenames[i++], "pseudoprimes_38bit.dat");
    strcpy(filenames[i++], "pseudoprimes_40bit.dat");
	strcpy(filenames[i++], "pseudoprimes_42bit.dat");
	strcpy(filenames[i++], "pseudoprimes_44bit.dat");
	strcpy(filenames[i++], "pseudoprimes_46bit.dat");
    //strcpy(filenames[i++], "pseudoprimes_41bit.dat");
    strcpy(filenames[i++], "pseudoprimes_48bit.dat");
	strcpy(filenames[i++], "pseudoprimes_50bit.dat");
	strcpy(filenames[i++], "pseudoprimes_52bit.dat");
	strcpy(filenames[i++], "pseudoprimes_54bit.dat");
	strcpy(filenames[i++], "pseudoprimes_56bit.dat");
	strcpy(filenames[i++], "pseudoprimes_58bit.dat");
    //strcpy(filenames[i++], "pseudoprimes_51bit.dat");
    //strcpy(filenames[i++], "pseudoprimes_55bit.dat");
    //strcpy(filenames[i++], "pseudoprimes_59bit.dat");
    strcpy(filenames[i++], "pseudoprimes_60bit.dat");
    //strcpy(filenames[i++], "pseudoprimes_61bit.dat");
    strcpy(filenames[i++], "pseudoprimes_62bit.dat");
    //strcpy(filenames[i++], "pseudoprimes_63bit.dat");
    strcpy(filenames[i++], "pseudoprimes_64bit.dat");
    num_files = 17;

	// lehman test
	for (nf = 0; nf < 0; nf++)
	{
		in = fopen(filenames[nf], "r");

		start = clock();
		i = 0;
		totBits = 0;
		minBits = 999;
		maxBits = 0;
		//read in everything
		while (!feof(in))
		{
			fscanf(in, "%" PRIu64 ",%u,%u", comp + i, f1 + i, f2 + i);
			sp642z(comp[i], &tmp);
			j = zBits(&tmp);
			totBits += j;
			if ((uint32)j > maxBits)
				maxBits = j;
			if ((uint32)j < minBits && j != 0)
				minBits = j;
			i++;
		}
		num = i;
		num = 100000;
		fclose(in);
		stop = clock();
		t_time = (double)(stop - start) / (double)CLOCKS_PER_SEC;
		printf("data read in %2.4f sec\n", t_time);
		printf("average bits of input numbers = %.2f\n", (double)totBits / (double)i);
		printf("minimum bits of input numbers = %d\n", minBits);
		printf("maximum bits of input numbers = %d\n", maxBits);

		gettimeofday(&gstart, NULL);

		correct = 0;
		k = 0;
		for (i = 0; i < num; i++)
		{
			f64 = LehmanFactor(comp[i], 3.5, 0, 0.1);
			if ((f64 == f1[i]) || (f64 == f2[i]))
			{
				correct++;
			}
		}

		gettimeofday(&gstop, NULL);
		t_time = my_difftime(&gstart, &gstop);

		printf("Lehman got %d of %d correct in %2.2f sec\n", correct, num, t_time);
		printf("percent correct = %.2f\n", 100.0*(double)correct / (double)num);
		printf("average time per input = %1.4f ms\n", 1000 * t_time / (double)num);
	}



	// spbrent test
    for (nf = 0; nf < 13; nf++)
    {
        in = fopen(filenames[nf], "r");

        start = clock();
        i = 0;
        totBits = 0;
        minBits = 999;
        maxBits = 0;
        //read in everything
        while (!feof(in))
        {
            fscanf(in, "%" PRIu64 ",%u,%u", comp + i, f1 + i, f2 + i);
            sp642z(comp[i], &tmp);
            j = zBits(&tmp);
            totBits += j;
            if ((uint32)j > maxBits)
                maxBits = j;
            if ((uint32)j < minBits && j != 0)
                minBits = j;
            i++;
        }
        num = i;
        num = 100000;
        fclose(in);
        stop = clock();
        t_time = (double)(stop - start) / (double)CLOCKS_PER_SEC;
        printf("data read in %2.4f sec\n", t_time);
        printf("average bits of input numbers = %.2f\n", (double)totBits / (double)i);
        printf("minimum bits of input numbers = %d\n", minBits);
        printf("maximum bits of input numbers = %d\n", maxBits);

        gettimeofday(&gstart, NULL);

        correct = 0;
        k = 0;
        for (i = 0; i < num; i++)
        {
            int p;
            
            for (p = 0; p < 1; p++)
            {
                //f64 = spbrent64(comp[i], p + 1, 8192);
				f64 = spbrent64(comp[i], 8192);
                if ((f64 == f1[i]) || (f64 == f2[i]))
                {
                    correct++;
                    break;
                }
            }           

            if (p == 1)
            {
                mpz_set_64(gmptmp, comp[i]);
                f64 = sp_shanks_loop(gmptmp, NULL);
                if ((f64 == f1[i]) || (f64 == f2[i]))
                {
                    correct++;
                    k++;
                }
            }
        }

        gettimeofday(&gstop, NULL);
        t_time = my_difftime(&gstart, &gstop);

        printf("rho got %d of %d correct in %2.2f sec (%d by squfof)\n", correct, num, t_time, k);
        printf("percent correct = %.2f\n", 100.0*(double)correct / (double)num);
        printf("average time per input = %1.4f ms\n", 1000 * t_time / (double)num);
    }

	goto tinyecm_marker;

	// sequential squfof test
	for (nf = 0; nf < 0; nf++)
	{
		in = fopen(filenames[nf], "r");

		start = clock();
		i = 0;
		totBits = 0;
		minBits = 999;
		maxBits = 0;
		//read in everything
		while (!feof(in))
		{
			fscanf(in, "%" PRIu64 ",%u,%u", comp + i, f1 + i, f2 + i);
			sp642z(comp[i], &tmp);
			j = zBits(&tmp);
			totBits += j;
			if ((uint32)j > maxBits)
				maxBits = j;
			if ((uint32)j < minBits && j != 0)
				minBits = j;
			i++;
		}
		num = i;
		num = 100000;
		fclose(in);
		stop = clock();
		t_time = (double)(stop - start) / (double)CLOCKS_PER_SEC;
		printf("data read in %2.4f sec\n", t_time);
		printf("average bits of input numbers = %.2f\n", (double)totBits / (double)i);
		printf("minimum bits of input numbers = %d\n", minBits);
		printf("maximum bits of input numbers = %d\n", maxBits);

		gettimeofday(&gstart, NULL);

		correct = 0;
		k = 0;
		for (i = 0; i < num; i++)
		{
			mpz_set_64(gmptmp, comp[i]);
			f64 = sp_shanks_loop(gmptmp, NULL);

			if (((uint32)f64 == f1[i]) || ((uint32)f64 == f2[i]))
				correct++;
			else
			{
				int p;
				for (p = 0; p < 3; p++)
				{
					f64 = spbrent(comp[i], p + 1, 8192);
					if ((f64 == f1[i]) || (f64 == f2[i]))
					{
						correct++;
						k++;
						break;
					}
				}
			}
		}

		gettimeofday(&gstop, NULL);
		t_time = my_difftime(&gstart, &gstop);

		printf("squfof got %d of %d correct in %2.2f sec (%d by rho)\n", correct, num, t_time, k);
		printf("percent correct = %.2f\n", 100.0*(double)correct / (double)num);
		printf("average time per input = %1.4f ms\n", 1000 * t_time / (double)num);
	}

	// smallmpqs test
    for (nf = 5; nf < 0; nf++)
    {
        in = fopen(filenames[nf], "r");

        start = clock();
        i = 0;
        totBits = 0;
        minBits = 999;
        maxBits = 0;
        //read in everything
        while (!feof(in))
        {
            fscanf(in, "%" PRIu64 ",%u,%u", comp + i, f1 + i, f2 + i);
            sp642z(comp[i], &tmp);
            j = zBits(&tmp);
            totBits += j;
            if ((uint32)j > maxBits)
                maxBits = j;
            if ((uint32)j < minBits && j != 0)
                minBits = j;
            i++;
        }
        num = i;
        num = 10000;
        fclose(in);
        stop = clock();
        t_time = (double)(stop - start) / (double)CLOCKS_PER_SEC;
        printf("data read in %2.4f sec\n", t_time);
        printf("average bits of input numbers = %.2f\n", (double)totBits / (double)i);
        printf("minimum bits of input numbers = %d\n", minBits);
        printf("maximum bits of input numbers = %d\n", maxBits);

        {
            gettimeofday(&gstart, NULL);

            correct = 0;
            k = 0;
            for (i = 0; i < num; i++)
            {
                mpz_set_64(fobj2->qs_obj.gmp_n, comp[i]);
                fobj2->qs_obj.flags = 12345;
                smallmpqs(fobj2);

                if (mpz_cmp_ui(fobj2->qs_obj.gmp_n, 1) == 0)
                {
                    correct++;
                }

                clear_factor_list(fobj2);
            }

            gettimeofday(&gstop, NULL);
            t_time = my_difftime(&gstart, &gstop);

            printf("smallmpqs got %d of %d correct in %2.2f sec\n", correct, num, t_time);
            printf("percent correct = %.2f\n", 100.0*(double)correct / (double)num);
            printf("average time per input = %.2f ms\n", 1000 * t_time / (double)num);
        }
    }

    i = 0;
    strcpy(filenames[i++], "pseudoprimes_55bit.dat");
    strcpy(filenames[i++], "pseudoprimes_60bit.dat");
    strcpy(filenames[i++], "pseudoprimes_64bit.dat");
    strcpy(filenames[i++], "pseudoprimes_70bit.dat");
    strcpy(filenames[i++], "pseudoprimes_80bit.dat");
    strcpy(filenames[i++], "pseudoprimes_90bit.dat");
    strcpy(filenames[i++], "pseudoprimes_100bit.dat");
    //strcpy(filenames[i++], "pseudoprimes_110bit.dat");
    //strcpy(filenames[i++], "pseudoprimes_120bit.dat");
    //strcpy(filenames[i++], "pseudoprimes_125bit.dat");
    num_files = 7;
    num = 10000;

	// smallmpqs test
    for (nf = 0; nf < 0; nf++)
    {
        mpz_t gmp_comp;
        char buf[1024];
        in = fopen(filenames[nf], "r");

        mpz_init(gmp_comp);

        gettimeofday(&gstart, NULL);

        correct = 0;
        k = 0;
        for (i = 0; i < num; i++)
        {
            fgets(buf, 1024, in);
            gmp_sscanf(buf, "%Zd", gmp_comp);
            mpz_set(fobj2->qs_obj.gmp_n, gmp_comp);
            fobj2->qs_obj.flags = 12345;
            smallmpqs(fobj2);

            if (mpz_cmp_ui(fobj2->qs_obj.gmp_n, 1) == 0)
            {
                correct++;
            }

            clear_factor_list(fobj2);
        }

        fclose(in);

        gettimeofday(&gstop, NULL);
        t_time = my_difftime(&gstart, &gstop);

        printf("smallmpqs got %d of %d correct in %2.2f sec\n", correct, num, t_time);
        printf("percent correct = %.2f\n", 100.0*(double)correct / (double)num);
        printf("average time per input = %.2f ms\n", 1000 * t_time / (double)num);

        mpz_clear(gmp_comp);
    }
        

    i = 0;
    strcpy(filenames[i++], "pseudoprimes_48bit.dat");
    strcpy(filenames[i++], "pseudoprimes_51bit.dat");
    strcpy(filenames[i++], "pseudoprimes_55bit.dat");
    strcpy(filenames[i++], "pseudoprimes_60bit.dat");
    strcpy(filenames[i++], "pseudoprimes_63bit.dat");
    strcpy(filenames[i++], "pseudoprimes_64bit.dat");
    strcpy(filenames[i++], "pseudoprimes_70bit.dat");
    strcpy(filenames[i++], "pseudoprimes_80bit.dat");
    strcpy(filenames[i++], "pseudoprimes_90bit.dat");
    strcpy(filenames[i++], "pseudoprimes_100bit.dat");
    num_files = 8;
    num = 10000;

	// arbitrary precision rho test
    for (nf = 0; nf < 0; nf++)
    {
        mpz_t gmp_comp;
        char buf[1024];
        in = fopen(filenames[nf], "r");

        mpz_init(gmp_comp);

        gettimeofday(&gstart, NULL);

        correct = 0;
        k = 0;
        for (i = 0; i < num; i++)
        {
            int p;

            fgets(buf, 1024, in);
            gmp_sscanf(buf, "%Zd", gmp_comp);            
            mpz_set(fobj2->rho_obj.gmp_n, gmp_comp);
            fobj2->rho_obj.iterations = 8192;
            
            for (p = 0; p < 3; p++)
            {
                fobj2->rho_obj.curr_poly = p;
                mbrent(fobj2);

                if (mpz_cmp_ui(fobj2->rho_obj.gmp_f, 1) > 0)
                {
                    correct++;
                    if (p > 0)
                        k++;
                    break;
                }
            }

            clear_factor_list(fobj2);
        }

        fclose(in);

        gettimeofday(&gstop, NULL);
        t_time = my_difftime(&gstart, &gstop);

        printf("mbrent got %d of %d correct in %2.2f sec (%d with backup polynomial)\n", 
            correct, num, t_time, k);
        printf("percent correct = %.2f\n", 100.0*(double)correct / (double)num);
        printf("average time per input = %.2f ms\n", 1000 * t_time / (double)num);

        mpz_clear(gmp_comp);
    }

    i = 0;
    strcpy(filenames[i++], "pseudoprimes_55bit.dat");
    strcpy(filenames[i++], "pseudoprimes_60bit.dat");
    strcpy(filenames[i++], "pseudoprimes_63bit.dat");
    strcpy(filenames[i++], "pseudoprimes_64bit.dat");
    strcpy(filenames[i++], "pseudoprimes_70bit.dat");
    strcpy(filenames[i++], "pseudoprimes_80bit.dat");
    strcpy(filenames[i++], "pseudoprimes_90bit.dat");
    strcpy(filenames[i++], "pseudoprimes_100bit.dat");
    strcpy(filenames[i++], "pseudoprimes_110bit.dat");
    strcpy(filenames[i++], "pseudoprimes_120bit.dat");
    strcpy(filenames[i++], "pseudoprimes_125bit.dat");
    num_files = 11;
    num = 10000;

	// gmp-ecm test
    for (nf = 0; nf < 7; nf++)
    {
        mpz_t gmp_comp;
        char buf[1024];
        ecm_params my_ecm_params;
        int curves;
        int B1;
        
        ecm_init(my_ecm_params);

        switch (nf)
        {
        case 0:
			// 55 bit
            B1 = 100;
            curves = 20;
            break;
        case 1:
			// 60 bit
            B1 = 200;
            curves = 20;
            break;
        case 2:
        case 3:
			// 64 bit
            B1 = 400;
            curves = 50;
            break;
        case 4:
			// 70 bit
            B1 = 600;
            curves = 50;
            break;
        case 5:
			// 80 bit
            B1 = 800;
            curves = 50;
            break;
        case 6:
			// 90 bit
            B1 = 1000;
            curves = 50;
            break;
        case 7:
            B1 = 1200;
            curves = 50;
            break;
        case 8:
            B1 = 1400;
            curves = 50;
            break;
        default:
            B1 = 2000;
            curves = 50;
            break;
        }

        in = fopen(filenames[nf], "r");

        mpz_init(gmp_comp);

        gettimeofday(&gstart, NULL);

        correct = 0;
        k = 0;
        for (i = 0; i < num; i++)
        {
            int p;
			uint64 f1, f2;

            fgets(buf, 1024, in);
            gmp_sscanf(buf, "%Zd,%lu,%lu", gmp_comp, &f1, &f2);
            
            for (k = 0; k < curves; k++)
            {
                uint64 sigma = spRand(100, 1000000000);
                my_ecm_params->B1done = 1.0 + floor(1 * 128.) / 134217728.;
				//mpz_set_ui(my_ecm_params->B2, 0);
                mpz_set_ui(my_ecm_params->x, (unsigned long)0);
                mpz_set_ui(my_ecm_params->sigma, sigma);
                my_ecm_params->method = ECM_ECM;
                ecm_factor(fobj2->ecm_obj.gmp_f, gmp_comp, B1, my_ecm_params);
                //if ((mpz_cmp_ui(fobj2->ecm_obj.gmp_f, 1) > 0) &&
                //    (mpz_cmp(fobj2->ecm_obj.gmp_f, gmp_comp) < 0))
				if ((mpz_cmp_ui(fobj2->ecm_obj.gmp_f, f1) == 0) ||
					(mpz_cmp_ui(fobj2->ecm_obj.gmp_f, f2) == 0))
                {
                    correct++;
                    break;
                }
            }
        }

        fclose(in);

        gettimeofday(&gstop, NULL);
        t_time = my_difftime(&gstart, &gstop);

        printf("ecm got %d of %d correct in %2.2f sec\n",
            correct, num, t_time);
        printf("percent correct = %.2f\n", 100.0*(double)correct / (double)num);
        printf("average time per input = %.2f ms\n", 1000 * t_time / (double)num);

        mpz_clear(gmp_comp);
    }


tinyecm_marker:
	i = 0;
	strcpy(filenames[i++], "pseudoprimes_42bit.dat");		// 70
	strcpy(filenames[i++], "pseudoprimes_44bit.dat");		// 70
	strcpy(filenames[i++], "pseudoprimes_46bit.dat");		// 70
	strcpy(filenames[i++], "pseudoprimes_48bit.dat");		// 70
	strcpy(filenames[i++], "pseudoprimes_50bit.dat");		// 70
	strcpy(filenames[i++], "pseudoprimes_52bit.dat");		// 85
	//strcpy(filenames[i++], "pseudoprimes_55bit.dat");
	strcpy(filenames[i++], "pseudoprimes_56bit.dat");		// 125
	strcpy(filenames[i++], "pseudoprimes_58bit.dat");		// 125
	strcpy(filenames[i++], "pseudoprimes_60bit.dat");		// 165
	strcpy(filenames[i++], "pseudoprimes_62bit.dat");		// 165
	strcpy(filenames[i++], "pseudoprimes_64bit.dat");		// 205
	//strcpy(filenames[i++], "pseudoprimes_70bit.dat");
	strcpy(filenames[i++], "semiprimes_tlp_32x32x32.txt");
	strcpy(filenames[i++], "semiprimes_tlp_32x64.txt");
	strcpy(filenames[i++], "semiprimes_tlp_48x48.txt");
	strcpy(filenames[i++], "pseudoprimes_80bit.dat");
	strcpy(filenames[i++], "pseudoprimes_90bit.dat");
	strcpy(filenames[i++], "pseudoprimes_100bit.dat");
	strcpy(filenames[i++], "pseudoprimes_110bit.dat");
	strcpy(filenames[i++], "pseudoprimes_120bit.dat");
	strcpy(filenames[i++], "pseudoprimes_125bit.dat");
	num_files = 11;
	num = 1000;

	// tinyecm test
	for (nf = 11; nf < 14; nf++)
	{
		mpz_t gmp_comp, gmp_f;
		char buf[1024];
		int curves;
		int B1;
		uint64 known1, known2, known3;

		switch (nf)
		{
		case 0:	// probably needs D=30 and smaller B1
		case 1: // probably needs D=30 and smaller B1
		case 2:
		case 3:
		case 4:				// 42 - 50 bit
			B1 = 70;
			curves = 24;
			break;
		case 5:				// 52 bit
			B1 = 85;
			curves = 24;
			break;
		case 6:				// 54-55 bit
		case 7:				
			B1 = 125;
			curves = 24;
			break;
		case 8:				// 60 bit
		case 9:
			B1 = 165;
			curves = 32;
			break;
		case 10:				// 63-64 bit
			B1 = 205;
			curves = 40;
			break;
		//case 4:	// 70 bit
		//	B1 = 400;
		//	curves = 16;
		//	break;
		//case 5:	// 80 bit
		//	B1 = 500;
		//	curves = 16;
		//	break;
		//case 6:	// 90 bit
		//	B1 = 600;
		//	curves = 32;
		//	break;
		//case 7:
		//	B1 = 700;
		//	curves = 32;
		//	break;
		//case 8:
		//	B1 = 800;
		//	curves = 32;
		//	break;
		default:
			B1 = 205;
			curves = 24;
			break;
		}

		in = fopen(filenames[nf], "r");
		printf("testing file: %s\n", filenames[nf]);

		mpz_init(gmp_comp);
		mpz_init(gmp_f);

		gettimeofday(&gstart, NULL);

		correct = 0;
		k = 0;
		if (nf < 0)
		{
			num = 100000;
			for (i = 0; i < num; i++)
			{
				uint64 in64;
				uint64 outf;
				
				fgets(buf, 1024, in);
				sscanf(buf, "%lu, %lu, %lu", &in64, &known1, &known2);

				microecm(in64, &outf, B1, 25 * B1, curves, 0);
				if ((outf == known1) ||
					(outf == known2))
				{
					correct++;
				}
			}

			fclose(in);

			gettimeofday(&gstop, NULL);
			t_time = my_difftime(&gstart, &gstop);

			printf("microecm got %d of %d correct in %2.2f sec\n",
				correct, num, t_time);
			printf("percent correct = %.2f\n", 100.0*(double)correct / (double)num);
			printf("average time per input = %.2f ms\n", 1000 * t_time / (double)num);
		}
		else
		{
			num = 10000;
			if (nf == 11)
			{
				for (i = 0; i < num; i++)
				{
					int p;

					fgets(buf, 1024, in);
					gmp_sscanf(buf, "%Zd, %lu, %lu, %lu", gmp_comp, &known1, &known2, &known3);

					tinyecm(gmp_comp, gmp_f, B1, 25 * B1, curves, 0);
					if ((mpz_cmp_ui(gmp_f, known1) == 0) ||
						(mpz_cmp_ui(gmp_f, known2) == 0) ||
						(mpz_cmp_ui(gmp_f, known3) == 0))
					{
						correct++;
					}
				}

			}
			else if (nf > 11)
			{
				for (i = 0; i < num; i++)
				{
					int p;

					fgets(buf, 1024, in);
					gmp_sscanf(buf, "%Zd, %lu, %lu", gmp_comp, &known1, &known2);

					tinyecm(gmp_comp, gmp_f, B1, 25 * B1, curves, 0);
					if ((mpz_cmp_ui(gmp_f, known1) == 0) ||
						(mpz_cmp_ui(gmp_f, known2) == 0))
					{
						correct++;
					}
				}

			}
			else
			{
				for (i = 0; i < num; i++)
				{
					int p;

					fgets(buf, 1024, in);
					gmp_sscanf(buf, "%Zd, %lu, %lu", gmp_comp, &known1, &known2);

					tinyecm(gmp_comp, gmp_f, B1, 25 * B1, curves, 0);
					if ((mpz_cmp_ui(gmp_f, known1) == 0) ||
						(mpz_cmp_ui(gmp_f, known2) == 0))
					{
						correct++;
					}
				}
			}

			fclose(in);

			gettimeofday(&gstop, NULL);
			t_time = my_difftime(&gstart, &gstop);

			printf("tinyecm got %d of %d correct in %2.2f sec\n",
				correct, num, t_time);
			printf("percent correct = %.2f\n", 100.0*(double)correct / (double)num);
			printf("average time per input = %.2f ms\n", 1000 * t_time / (double)num);
		}

		mpz_clear(gmp_comp);
		mpz_clear(gmp_f);	
	}

	exit(0);

	free(f1);
	free(f2);
	free(comp);
	zFree(&tmp);
	zFree(&tmp2);
	zFree(&t2);
	zFree(&t1);
	zFree(&t3);
	mpz_clear(gmptmp);
	return;
}

void test_qsort(void)
{
	//test the speed of qsort in  sorting a few million lists
	//of a few thousand random integers

	uint32 *list;
	int i,j,k;
	const int listsz = 512;
	double t_time;
	clock_t start, stop;

	list = (uint32 *)malloc(listsz * sizeof(uint32));

	for (k=1000;k<1000000;k*=10)
	{
		start = clock();

		for (j=0; j<k; j++)
		{
			for (i=0; i<listsz; i++)
				list[i] = (uint32)spRand(1000000,50000000);
		}

		stop = clock();
		t_time = (double)(stop - start)/(double)CLOCKS_PER_SEC;
		printf("baseline elapsed time for %d lists = %2.2f\n",k,t_time);

		start = clock();

		for (j=0; j<k; j++)
		{
			for (i=0; i<listsz; i++)
				list[i] = (uint32)spRand(1000000,50000000);
			qsort(list,listsz,4,&qcomp_uint32);
		}

		stop = clock();
		t_time = (double)(stop - start)/(double)CLOCKS_PER_SEC;
		printf("elapsed time for %d sorts = %2.2f\n",k,t_time);
	}

	free(list);
	return;
}
	
void arith_timing(int num)
{

	int i,sz;
	double t_time;
	clock_t start, stop;
	z a,b,c,d,e;

	zInit(&a);
	zInit(&b);
	zInit(&c);
	zInit(&d);
	zInit(&e);
	
	//spAcc(num);
	
	//sqrt_acc(num, int sz);
	//gcd_acc(num, int sz);
	goto muldiv;

	for (sz=50; sz<=500; sz += 50)
	{
		printf("baseline: generating %d random numbers with %d digits: ",num,sz);
		start = clock();
		for (i=0;i<num;i++)
		{
			zRand(&a,sz);
			zRand(&b,sz);
		}
		stop = clock();
		t_time = (double)(stop - start)/(double)CLOCKS_PER_SEC;
		printf("elapsed time = %2.3f\n",t_time);
	}

	for (sz=50; sz<=500; sz += 50)
	{
		subadd_acc(num, sz);
		/*
		zRand(&a,sz);
		zRand(&b,sz);
		printf("adding %d random numbers with %d digits: ",num,sz);
		start = clock();
		for (i=0;i<num;i++)
		{
			//zRand(&a,sz);
			//zRand(&b,sz);
			zAdd(&a,&b,&c);
		}
		stop = clock();
		t_time = (double)(stop - start)/(double)CLOCKS_PER_SEC;
		printf("elapsed time = %2.3f\n",t_time);
		*/
	}

	for (sz=50; sz<=500; sz += 50)
	{
		shortmuldiv_acc(num, sz);

		/*
		zRand(&a,sz);
		zRand(&b,sz);
		printf("subtracting %d random numbers with %d digits: ",num,sz);
		start = clock();
		for (i=0;i<num;i++)
		{
			//zRand(&a,sz);
			//zRand(&b,sz);
			zSub(&a,&b,&c);
		}
		stop = clock();
		t_time = (double)(stop - start)/(double)CLOCKS_PER_SEC;
		printf("elapsed time = %2.3f\n",t_time);
		*/
	}

muldiv:
	//for (sz=50; sz<=500; sz += 50)
	for (sz=5; sz<=50; sz += 1)
	{
		muldiv_acc(num, sz);

		/*
		zRand(&a,sz);
		zRand(&b,sz);
		printf("Multiplying %d random numbers with %d digits: ",num,sz);
		start = clock();
		for (i=0;i<num;i++)
		{
			//zRand(&a,sz);
			//zRand(&b,sz);
			zMul(&a,&b,&c);
		}
		stop = clock();
		t_time = (double)(stop - start)/(double)CLOCKS_PER_SEC;
		printf("elapsed time = %2.3f\n",t_time);
		*/
	}

	goto sqrt_test;

//sqr_test:
	for (sz=50; sz<=500; sz += 50)
	{
		sqr_acc(num, sz);

		/*
		zRand(&a,sz);
		zRand(&b,sz);
		printf("Squaring %d random numbers with %d digits: ",num,sz);
		start = clock();
		for (i=0;i<num;i++)
		{
			//zRand(&a,sz);
			//zRand(&b,sz);
			zSqr(&a,&c);
		}
		stop = clock();
		t_time = (double)(stop - start)/(double)CLOCKS_PER_SEC;
		printf("elapsed time = %2.3f\n",t_time);
		*/
	}

sqrt_test:	
	for (sz=5; sz<=50; sz += 1)
	//for (sz=50; sz<=150; sz += 10)
	{
		//gcd_acc(num/100,sz);
		sqrt_acc(num,sz);

		/*
		zRand(&a,sz);
		zRand(&b,sz/2);
		printf("Dividing %d random numbers with %d digits: ",num,sz);
		start = clock();
		for (i=0;i<num;i++)
		{
			//zRand(&a,sz);
			//zRand(&b,sz/2);
			zCopy(&a,&e);
			zDiv(&a,&b,&c,&d);
			zCopy(&e,&a);
		}
		stop = clock();
		t_time = (double)(stop - start)/(double)CLOCKS_PER_SEC;
		printf("elapsed time = %2.3f\n",t_time);
		*/
	}
	

	zFree(&a);
	zFree(&b);
	zFree(&c);
	zFree(&d);
	zFree(&e);

	return;
}

double spAcc(int m)
{
	z a,b,c,d;
	unsigned int shift;
	fp_digit bitmask,k,v;
	fp_digit tmpd[2];
	int j;

	clock_t start, stop;
	double t=0.0;

	zInit(&a);
	zInit(&b);
	zInit(&c);
	zInit(&d);

	printf("\nAccuracy test of single precision stuff, %d iterations\n",m);
	start = clock();
	for (j=1; j<m; ++j)
	{
		a.size=2;
		b.size=1;
		a.val[0] = rand()*rand();
		a.val[1] = rand()*rand();
		b.val[0] = 0;
		while (!(b.val[0]))
			b.val[0] = rand()*rand();

		bitmask = HIBITMASK;
		for (shift = 0; shift < BITS_PER_DIGIT; shift++)
		{
			if (b.val[0] & bitmask)
				break;
			bitmask >>= 1;
		}

		b.val[0] <<= shift;
		zShiftLeft(&a,&a,shift);

		tmpd[1] = a.val[1]; tmpd[0] = a.val[0];
		c.val[2] = spDivide(&c.val[1],&c.val[0],tmpd,b.val[0]);

		spMultiply(c.val[1],b.val[0],&d.val[0],&d.val[1]);
		if (c.val[2])
			d.val[1] += b.val[0];
		spAdd(d.val[0],c.val[0],&k,&v);
		spAdd(d.val[1],v,&d.val[1],&v);
		if ((k != a.val[0]) || (d.val[1] != a.val[1]))
		{
			printf("error at %d\na: ",j);
			break;
		}
	}
	stop = clock();
	t = (double)(stop - start)/(double)CLOCKS_PER_SEC;
	printf("%d numbers verified.  Elapsed time = %6.4f seconds.\n", m,t);

	zFree(&a);
	zFree(&b);
	zFree(&c);
	zFree(&d);

	return t;
}

double sqr_acc(int m, int sz)
{
	int i,j;
	clock_t start, stop;
	double t=0.0;

	z a,b,c,d;
	zInit(&a);
	zInit(&b);
	zInit(&c);
	zInit(&d);

	printf("\nAccuracy test, square and verify by mul:\n");
	start = clock();
	for (j=1; j<m; ++j)
	{
		zRand(&a,sz);
		
		if (rand() > 16384)
			a.size *= -1;

		zSqr(&a,&b);
		//zNroot(&b,&c,2);
		zMul(&a,&a,&c);

		//if (a.size < 0)
		//	c.size *= -1;

		i = zCompare(&c,&b);

		if (!(i==0))
		{
			printf("failed at %d\n",j);
			printf("a = %s\nb = %s\nc = %s\n",
				z2decstr(&a,&gstr1),z2decstr(&b,&gstr2),z2decstr(&c,&gstr3));
			exit(-1);
		}
	}
	stop = clock();
	t = (double)(stop - start)/(double)CLOCKS_PER_SEC;
	printf("%d numbers verified.  Elapsed time = %6.4f seconds.\n", m,t);

	zFree(&a);
	zFree(&b);
	zFree(&c);
	zFree(&d);

	return t;
}

double subadd_acc(int m, int sz)
{

	int i,j;
	z a,b,c,d;

	clock_t start, stop;
	double t=0.0;

	zInit(&a);
	zInit(&b);
	zInit(&c);
	zInit(&d);

	printf("\nAccuracy test, subtract and verify by addition:\n");
	start = clock();
	for (j=1; j<m; ++j)
	{
		zRand(&a,sz);
		zRand(&b,sz);
		
		if (rand() > 16384)
			a.size *= -1;

		if (rand() > 16384)
			b.size *= -1;

		zAdd(&a,&b,&d);
		zSub(&d,&b,&c);

		i = zCompare(&c,&a);

		if (!(i==0))
		{
			printf("failed at %d\n",j);
			printf("a+b=d\nd-b=c\nassert a == c failed\na: ");
			break;
		}
	}
	stop = clock();
	t = (double)(stop - start)/(double)CLOCKS_PER_SEC;
	printf("%d numbers verified.  Elapsed time = %6.4f seconds.\n", m,t);

	zFree(&a);
	zFree(&b);
	zFree(&c);
	zFree(&d);
	return t;
}

double shortmuldiv_acc(int m, int sz)
{

	int i,j;
	z a,b,c,d,q;
	fp_digit v;

	clock_t start, stop;
	double t=0.0;

	zInit(&a);
	zInit(&b);
	zInit(&c);
	zInit(&d);
	zInit(&q);

	printf("\nAccuracy test, short divide and verify by short multiplication and addition:\n");
	start = clock();
	for (j=1;j<m;++j)
	{
		zRand(&a,sz);
		
		//make a non-zero digit 'v'
		v = spRand(1,MAX_DIGIT);

		b.val[0] = zShortDiv(&a,v,&q); 
		b.size=1;

		zShortMul(&q,v,&c);
		zAdd(&c,&b,&d);
		zShortAdd(&d,b.val[0],&c);
		zShortSub(&c,b.val[0],&d);
		i = zCompare(&a,&d);
		
		if (!(i==0)) 
		{
			printf("failed at %d\na: ",j);
			break;
		}
	}
	stop = clock();
	t = (double)(stop - start)/(double)CLOCKS_PER_SEC;
	printf("%d numbers verified.  Elapsed time = %6.4f seconds.\n", m,t);

	zFree(&a);
	zFree(&b);
	zFree(&c);
	zFree(&d);
	zFree(&q);
	return t;
}

double muldiv_acc(int m, int sz)
{
	int i,j;
	z a,b,c,d,q,rem,tmp;

	clock_t start, stop;
	double t=0.0;

	zInit(&a);
	zInit(&b);
	zInit(&c);
	zInit(&d);
	zInit(&q);
	zInit(&rem);
	zInit(&tmp);

	printf("\nAccuracy test, divide and verify by multiplication and addition:\n");
	start = clock();
	for (j=1;j<m;++j)
	{	
		zRand(&a,sz*2);

		do 
		{
			zRandb(&b,(int)((double)sz * 3.33));
		} while (zCompare(&b,&zZero) == 0);

		zCopy(&a,&d);	//because a will be overwritten
		zDiv(&a,&b,&q,&rem);
		zMul(&q,&b,&c);
		zAdd(&c,&rem,&a);

		i = zCompare(&a,&d);
		
		if (i != 0) 
		{
			printf("\nfailed at %d\n",j);
			printf("a = %s\nb = %s\nq = %s\n",
				z2decstr(&d,&gstr1),z2decstr(&b,&gstr2),z2decstr(&q,&gstr3));
			printf("r = %s\nqb = %s\nqb+r = %s\n",
				z2decstr(&rem,&gstr1),z2decstr(&c,&gstr2),z2decstr(&a,&gstr3));
			exit(-1);
			break;
		}
	}
	stop = clock();
	t = (double)(stop - start)/(double)CLOCKS_PER_SEC;
	printf("%d numbers verified.  Elapsed time = %6.4f seconds.\n", m,t);

	zFree(&a);
	zFree(&b);
	zFree(&c);
	zFree(&d);
	zFree(&q);
	zFree(&rem);
	zFree(&tmp);
	return t;
}

double sqrt_acc(int m, int sz)
{
	int i,j,k;
	z a,b,c,d,tmp;

	clock_t start, stop;
	double t=0.0;

	zInit(&a);
	zInit(&b);
	zInit(&c);
	zInit(&d);
	zInit(&tmp);

	printf("\nAccuracy test, Newton square root:\n");
	start = clock();
	for (j=1;j<m;++j)
	{
		zRand(&a,sz);
		k = zNroot(&a,&b,2);
		
		if (k >= 10000)
		{
			printf("error, max iterations at %d\na:       ",j);
			break;
		}
		
		zMul(&b,&b,&c);
		i = zCompare(&c,&a);
		if (!((i==0) || (i==MAX_DIGIT)))
		{
			printf("error, not <= 0 at %d, iterations: %d\na:        ",j,k);			
			printf("a = %s\nb = %s\nc = %s\n",
				z2decstr(&a,&gstr1),z2decstr(&b,&gstr2),z2decstr(&c,&gstr3));

			break;
		}
		zShortAdd(&b,1,&c);
		zMul(&c,&c,&d);
		i = zCompare(&d,&a);
		if (i<=0)
		{
			printf("error, not > 0 at %d, iterations: %d\na:       ",j,k);
			printf("a = %s\nb = %s\nd = %s\n",
				z2decstr(&a,&gstr1),z2decstr(&b,&gstr2),z2decstr(&d,&gstr3));
			break;
		}
	}
	stop = clock();

	t = (double)(stop - start)/(double)CLOCKS_PER_SEC;
	printf("%d numbers verified.  Elapsed time = %6.4f seconds.\n", j,t);

	zFree(&a);
	zFree(&b);
	zFree(&c);
	zFree(&d);
	zFree(&tmp);
	return t;
}

void richard_guy_problem_e7(void)
{
	uint64 n,j;
	double sum = 0.;
		
	free(PRIMES);
	PRIMES = GetPRIMESRange(spSOEprimes, szSOEp, NULL, 0, 
		100000000, &NUM_P);

	P_MIN = PRIMES[0]; 
	P_MAX = PRIMES[(uint32)NUM_P-1];	

	spSOEprimes = (uint32 *)realloc(spSOEprimes, 
		(size_t) (NUM_P * sizeof(uint32)));
	for (n=0; n<NUM_P; n++)
		spSOEprimes[n] = (uint32)PRIMES[n];
	szSOEp = NUM_P;

	j=0;
	for (n=1; n<1000000000000; )
	{
		if (j >= NUM_P)
		{
			j=0;

			free(PRIMES);
			PRIMES = GetPRIMESRange(spSOEprimes, szSOEp, NULL, P_MAX+1, 
				(uint64)P_MAX + 10000000000, &NUM_P);

			P_MIN = PRIMES[0]; 
			P_MAX = PRIMES[(uint32)NUM_P-1];	

			// print the odd iterations
			printf("%" PRIu64 ", %" PRIu64 ", %1.9f\n", n, PRIMES[j], sum);

			if ((n & 0x1) == 0)
				sum += (double)n++ / (double)PRIMES[j++];
			else
				sum -= (double)n++ / (double)PRIMES[j++];
					
			printf("%" PRIu64 ", %" PRIu64 ", %1.9f\n", n, PRIMES[j], sum);
	
		}

		if (j < (NUM_P - 8))
		{
			if (n & 0x1)
			{
				sum -= (double)n++ / (double)PRIMES[j++];
				sum += (double)n++ / (double)PRIMES[j++];
				sum -= (double)n++ / (double)PRIMES[j++];
				sum += (double)n++ / (double)PRIMES[j++];
				sum -= (double)n++ / (double)PRIMES[j++];
				sum += (double)n++ / (double)PRIMES[j++];
				sum -= (double)n++ / (double)PRIMES[j++];
				sum += (double)n++ / (double)PRIMES[j++];
			}
			else
			{
				sum += (double)n++ / (double)PRIMES[j++];
				sum -= (double)n++ / (double)PRIMES[j++];
				sum += (double)n++ / (double)PRIMES[j++];
				sum -= (double)n++ / (double)PRIMES[j++];
				sum += (double)n++ / (double)PRIMES[j++];
				sum -= (double)n++ / (double)PRIMES[j++];
				sum += (double)n++ / (double)PRIMES[j++];
				sum -= (double)n++ / (double)PRIMES[j++];
			}
		}
		else
		{
			if ((n & 0x1) == 0)
				sum += (double)n++ / (double)PRIMES[j++];
			else
				sum -= (double)n++ / (double)PRIMES[j++];
		}
	}
	printf("%" PRIu64 ", %1.9f\n", n-1, sum);

	exit(0);
}