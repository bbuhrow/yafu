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
#include "cofactorize.h"
#include <ecm.h>

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
	struct timeval gstart;
	struct timeval gstop;
	int nf;
	int num_files;
	char filenames[30][80];
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

    goto tinyqs_marker;
	//goto brent_marker;
	//goto tinyecm_marker;

	
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
    // msvc max of 62 bits
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

	//goto tinyecm_marker;

	// sequential squfof test
	for (nf = 0; nf < 16; nf++)
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

tinyqs_marker:

    i = 0;
    strcpy(filenames[i++], "pseudoprimes_56bit.dat");
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

    // tinyqs test
    {
        tiny_qs_params *params;
        mpz_t fact1, fact2, gmp_comp;

        mpz_init(fact1);
        mpz_init(fact2);
        mpz_init(gmp_comp);
        params = init_tinyqs();

        for (nf = 0; nf < 7; nf++)
        {
            uint64 known1, known2;
            char buf[1024];
            in = fopen(filenames[nf], "r");

            gettimeofday(&gstart, NULL);

            correct = 0;
            k = 0;
            for (i = 0; i < num; i++)
            {
                fgets(buf, 768, in);

#ifdef _MSC_VER
                gmp_sscanf(buf, "%Zd, %llu, %llu",
                    gmp_comp, &known1, &known2);
#else
                gmp_sscanf(buf, "%Zd,%u,%u", gmp_comp, &known1, &known2);
#endif

                
                mpz_set(fobj2->qs_obj.gmp_n, gmp_comp);
                
                k = tinyqs(params, gmp_comp, fact1, fact2);

                if (k > 0)
                {
                    correct++;
                }
            }

            fclose(in);

            gettimeofday(&gstop, NULL);
            t_time = my_difftime(&gstart, &gstop);

            printf("tinyqs got %d of %d correct in %2.2f sec\n", correct, num, t_time);
            printf("percent correct = %.2f\n", 100.0*(double)correct / (double)num);
            printf("average time per input = %.2f ms\n", 1000 * t_time / (double)num);
        }

        params = free_tinyqs(params);
        mpz_clear(gmp_comp);
        mpz_clear(fact1);
        mpz_clear(fact2);
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
    for (nf = 0; nf < 0; nf++)
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
    strcpy(filenames[i++], "pseudoprimes_36bit.dat");
    strcpy(filenames[i++], "pseudoprimes_38bit.dat");
    strcpy(filenames[i++], "pseudoprimes_40bit.dat");
	strcpy(filenames[i++], "pseudoprimes_42bit.dat");		// 70
	strcpy(filenames[i++], "pseudoprimes_44bit.dat");		// 70
	strcpy(filenames[i++], "pseudoprimes_46bit.dat");		// 70
	strcpy(filenames[i++], "pseudoprimes_48bit.dat");		// 70
	strcpy(filenames[i++], "pseudoprimes_50bit.dat");		// 70
	strcpy(filenames[i++], "pseudoprimes_52bit.dat");		// 85
	strcpy(filenames[i++], "pseudoprimes_54bit.dat");       // 85
	strcpy(filenames[i++], "pseudoprimes_56bit.dat");		// 125
	strcpy(filenames[i++], "pseudoprimes_58bit.dat");		// 125
	strcpy(filenames[i++], "pseudoprimes_60bit.dat");		// 165
	strcpy(filenames[i++], "pseudoprimes_62bit.dat");		// 165
	strcpy(filenames[i++], "pseudoprimes_64bit.dat");		// 205
	strcpy(filenames[i++], "semiprimes_tlp_32x32x32.txt");
	strcpy(filenames[i++], "semiprimes_tlp_32x64.txt");
	strcpy(filenames[i++], "semiprimes_tlp_48x48.txt");
    strcpy(filenames[i++], "pseudoprimes_70bit.dat");
	strcpy(filenames[i++], "pseudoprimes_80bit.dat");
	strcpy(filenames[i++], "pseudoprimes_90bit.dat");
	strcpy(filenames[i++], "pseudoprimes_100bit.dat");
	strcpy(filenames[i++], "pseudoprimes_110bit.dat");
	strcpy(filenames[i++], "pseudoprimes_120bit.dat");
	strcpy(filenames[i++], "pseudoprimes_125bit.dat");
	num_files = 24;
	num = 100000;

	// tinyecm test
	for (nf = 15; nf < 18; nf++)
	{
		mpz_t gmp_comp, gmp_f;
		char buf[1024];
		int curves;
		int B1;
		uint64 known1, known2, known3;

		switch (nf)
		{
		case 0:	
		case 1:             // <= 40-bit
		case 2:             
            B1 = 70;        // tinyecm
            //B1 = 47;        // uecm
            curves = 16;
            break;
		case 3:             // 42-44 bit
        case 4:  
            B1 = 70;        // tinyecm
            //B1 = 59;        // uecm
            curves = 16;
            break;
        case 5:             // 46-50 bit
        case 6:
		case 7:				
			B1 = 70;
			curves = 24;
			break;
		case 8:				// 52 bit
			B1 = 85;
			curves = 24;
			break;
		case 9:				// 54-55 bit
		case 10:				
			B1 = 125;
			curves = 24;
			break;
		case 11:				// 60 bit
		case 12:
			B1 = 165;
			curves = 32;
			break;
		case 13:				// 63-64 bit
			B1 = 205;
			curves = 40;
			break;
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

				sscanf(buf, "%" PRIu64 ", %" PRIu64 ", %" PRIu64 "", 
                    &in64, &known1, &known2);

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
			
			if ((nf == 15) || (nf == 16) || (nf == 17))
			{
                num = 10000;
				for (i = 0; i < num; i++)
				{
					int p;

					fgets(buf, 1024, in);
#ifdef _MSC_VER
                    gmp_sscanf(buf, "%Zd, %llu, %llu, %llu",
                        gmp_comp, &known1, &known2, &known3);
#else
                    gmp_sscanf(buf, "%Zd, %" PRIu64 ", %" PRIu64 ", %" PRIu64 "",
                        gmp_comp, &known1, &known2, &known3);
#endif
					

					tinyecm(gmp_comp, gmp_f, B1, 25 * B1, curves, 0);
					if ((mpz_cmp_ui(gmp_f, known1) == 0) ||
						(mpz_cmp_ui(gmp_f, known2) == 0) ||
						(mpz_cmp_ui(gmp_f, known3) == 0))
					{
						correct++;
					}
				}

			}
			else if (nf > 17)
			{
                num = 10000;
				for (i = 0; i < num; i++)
				{
					int p;

					fgets(buf, 1024, in);
#ifdef _MSC_VER
                    gmp_sscanf(buf, "%Zd, %llu, %llu",
                        gmp_comp, &known1, &known2);
#else
                    gmp_sscanf(buf, "%Zd, %" PRIu64 ", %" PRIu64 "",
                        gmp_comp, &known1, &known2);
#endif
					

					tinyecm(gmp_comp, gmp_f, B1, 25 * B1, curves, 0);
					if ((mpz_cmp_ui(gmp_f, known1) == 0) ||
						(mpz_cmp_ui(gmp_f, known2) == 0))
					{
						correct++;
					}
                    //else
                    //{
                    //    gmp_printf("found incorrect factor %Zd (known: %" PRIu64 ", %" PRIu64 "\n",
                    //        gmp_f, known1, known2);
                    //}
				}

			}
			else
			{
                num = 100000;
				for (i = 0; i < num; i++)
				{
					int p;

					fgets(buf, 1024, in);
#ifdef _MSC_VER
                    uint32 k1, k2;
                    gmp_sscanf(buf, "%Zd, %u, %u",
                        gmp_comp, &k1, &k2);
                    known1 = k1;
                    known2 = k2;
#else
                    gmp_sscanf(buf, "%Zd, %" PRIu64 ", %" PRIu64 "",
                        gmp_comp, &known1, &known2);
#endif

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

    return;

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

