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

#include <inttypes.h>
#include "yafu.h"
#include "microecm.h"
#include "microecm_edw.h"
#include "tinyecm.h"
#include "tinyecm_edw.h"
#include "arith.h"
#include "ytools.h"
#include "qs.h"
#include "factor.h"
#include "monty.h"
#include "cofactorize.h"
#include "smallmpqs.h"
//#include "..\factor\kmpqs\kmpqs.h"
#include "kmpqs.h"
#ifdef HAVE_GMP_ECM
#include <ecm.h>
#endif
#include <math.h>

#ifdef __MINGW32__
#include <sys/time.h>
#endif

int test_smallmpqs(mpz_t* inputs, int num_inputs)
{
	struct timeval gstart;
	struct timeval gstop;
	mpz_t n, r;
	mpz_init(n);
	mpz_init(r);
	gettimeofday(&gstart, NULL);

	int i;
	int correct = 0;
	for (i = 0; i < num_inputs; i++)
	{
		int numf;
		mpz_t* factors = smallmpqs(inputs[i], &numf);

		if (numf > 0)
		{
			mpz_set(n, inputs[i]);
			int j;
			for (j = 0; j < numf; j++)
			{
				mpz_tdiv_r(r, n, factors[j]);
				if (mpz_cmp_ui(r, 0) == 0)
				{
					mpz_tdiv_q(n, n, factors[j]);
				}
			}

			if (mpz_cmp_ui(n, 1) == 0)
			{
				correct++;
			}
		}

	}

	gettimeofday(&gstop, NULL);
	double t_time = ytools_difftime(&gstart, &gstop);

	mpz_clear(n);
	mpz_clear(r);
	printf("smallmpqs got %d of %d correct in %2.2f sec\n", correct, num_inputs, t_time);
	printf("percent correct = %.2f\n", 100.0 * (double)correct / (double)num_inputs);
	printf("average time per input = %1.4f ms\n", 1000 * t_time / (double)num_inputs);
	return correct;
}

int test_tinysiqs(mpz_t* inputs, int num_inputs)
{
	struct timeval gstart;
	struct timeval gstop;
	mpz_t n, r, f1, f2, f3;
	mpz_init(n);
	mpz_init(r);
	mpz_init(f1);
	mpz_init(f2);
	mpz_init(f3);

	tiny_qs_params* params;
	params = init_tinysiqs();

	gettimeofday(&gstart, NULL);

	int i;
	int correct = 0;
	for (i = 0; i < num_inputs; i++)
	{
		int numf = (int)tinysiqs(params, inputs[i], f1, f2, f3, 0);

		if (numf > 0)
		{
			mpz_set(n, inputs[i]);

			if (mpz_cmp_ui(f1, 1) > 0)
			{
				mpz_tdiv_r(r, n, f1);
				if (mpz_cmp_ui(r, 0) == 0)
				{
					mpz_tdiv_q(n, n, f1);
				}
			}
			if (mpz_cmp_ui(f2, 1) > 0)
			{
				mpz_tdiv_r(r, n, f2);
				if (mpz_cmp_ui(r, 0) == 0)
				{
					mpz_tdiv_q(n, n, f2);
				}
			}
			if (mpz_cmp_ui(f3, 1) > 0)
			{
				mpz_tdiv_r(r, n, f3);
				if (mpz_cmp_ui(r, 0) == 0)
				{
					mpz_tdiv_q(n, n, f3);
				}
			}
			if ((mpz_cmp_ui(n, 1) == 0) || (mpz_probab_prime_p(n, 1) > 0))
			{
				correct++;
			}
		}
	}

	gettimeofday(&gstop, NULL);
	double t_time = ytools_difftime(&gstart, &gstop);

	free_tinysiqs(params);
	mpz_clear(n);
	mpz_clear(r);
	mpz_clear(f1);
	mpz_clear(f2);
	mpz_clear(f3);
	printf("tinysiqs got %d of %d correct in %2.2f sec\n", correct, num_inputs, t_time);
	printf("percent correct = %.2f\n", 100.0 * (double)correct / (double)num_inputs);
	printf("average time per input = %1.4f ms\n", 1000 * t_time / (double)num_inputs);
	return correct;
}

int test_kmpqs(mpz_t* inputs, int num_inputs)
{
	struct timeval gstart;
	struct timeval gstop;
	mpz_t n, r, *factors;
	mpz_init(n);
	mpz_init(r);

	kmpqs_ctx_t* ctx = kmpqs_ctx_new();

	gettimeofday(&gstart, NULL);

	int i;
	int correct = 0;
	for (i = 0; i < num_inputs; i++)
	{
		int numf = kmpqs_factor_ctx(ctx, inputs[i], 60, &factors);

		if (numf > 0)
		{
			mpz_set(n, inputs[i]);
			int j;

			for (j = 0; j < numf; j++) {
				if (mpz_cmp_ui(factors[j], 1) > 0)
				{
					mpz_tdiv_r(r, n, factors[j]);
					if (mpz_cmp_ui(r, 0) == 0)
					{
						mpz_tdiv_q(n, n, factors[j]);
					}
				}
			}

			if ((mpz_cmp_ui(n, 1) == 0) || (mpz_probab_prime_p(n, 1) > 0))
			{
				correct++;
			}
		}
	}

	gettimeofday(&gstop, NULL);
	double t_time = ytools_difftime(&gstart, &gstop);

	kmpqs_ctx_free(ctx);
	mpz_clear(n);
	mpz_clear(r);
	printf("kmpqs got %d of %d correct in %2.2f sec\n", correct, num_inputs, t_time);
	printf("percent correct = %.2f\n", 100.0 * (double)correct / (double)num_inputs);
	printf("average time per input = %1.4f ms\n", 1000 * t_time / (double)num_inputs);
	return correct;
}

void test_dlp_composites(int testnum)
{
	FILE *in;
	uint64_t *comp, f64;
	uint32_t *f1, *f2;
	uint64_t* f64a, *f64b, *outf;
	uint32_t bits, totBits, minBits, maxBits;
	double t_time;
	clock_t start, stop;
	int i, j, k, num, correct;
	mpz_t gmptmp;
	struct timeval gstart;
	struct timeval gstop;
	int nf;
	int num_files;
	char filenames[30][80];
	fact_obj_t *fobj2;
	uint64_t lcg = 42;

	fobj2 = (fact_obj_t *)malloc(sizeof(fact_obj_t));
	init_factobj(fobj2);

	mpz_init(gmptmp);
	comp = (uint64_t*)malloc(2000000 * sizeof(uint64_t));
	outf = (uint64_t*)malloc(2000000 * sizeof(uint64_t));
	f1 = (uint32_t *)malloc(2000000 * sizeof(uint32_t));
	f2 = (uint32_t *)malloc(2000000 * sizeof(uint32_t));
	f64a = (uint64_t*)malloc(2000000 * sizeof(uint64_t));
	f64b = (uint64_t*)malloc(2000000 * sizeof(uint64_t));

	switch (testnum)
	{
	case 0: goto tinyqs_marker;
	case 1: goto kmpqs_marker;
	case 2: goto tinyecm_104_list_marker;
	case 3: goto tinyecm_128_marker;
	case 4: goto microecm_marker;
	case 5: goto uecm_52_list_marker;
	case 6: goto brent_marker;
	case 7: goto mpqs_marker;
	case 8: goto microecm_edw_marker;
	case 9: goto tinyecm_104_edw_list_marker;
	case 10: goto tinyecm_128_edw_marker;
	}

brent_marker:

    i = 0;
	strcpy(filenames[i++], "semiprimes_32bit.dat");
	strcpy(filenames[i++], "semiprimes_34bit.dat");
	strcpy(filenames[i++], "semiprimes_36bit.dat");
	strcpy(filenames[i++], "semiprimes_38bit.dat");
    strcpy(filenames[i++], "semiprimes_40bit.dat");
	strcpy(filenames[i++], "semiprimes_42bit.dat");
	strcpy(filenames[i++], "semiprimes_44bit.dat");
	strcpy(filenames[i++], "semiprimes_46bit.dat");
    strcpy(filenames[i++], "semiprimes_48bit.dat");
	strcpy(filenames[i++], "semiprimes_50bit.dat");
	strcpy(filenames[i++], "semiprimes_52bit.dat");
	strcpy(filenames[i++], "semiprimes_54bit.dat");
	strcpy(filenames[i++], "semiprimes_56bit.dat");
	strcpy(filenames[i++], "semiprimes_58bit.dat");
    strcpy(filenames[i++], "semiprimes_60bit.dat");
    strcpy(filenames[i++], "semiprimes_62bit.dat");
    strcpy(filenames[i++], "semiprimes_64bit.dat");
    num_files = i;

	// lehman test
	for (nf = 0; nf < 5; nf++)
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
			if (fscanf(in, "%" PRIu64 ",%u,%u", comp + i, f1 + i, f2 + i) != 3) break;
            mpz_set_u64(gmptmp, comp[i]);
            j = (int)mpz_sizeinbase(gmptmp, 2);
			totBits += (uint32_t)j;
			if ((uint32_t)j > maxBits)
				maxBits = (uint32_t)j;
			if ((uint32_t)j < minBits && j != 0)
				minBits = (uint32_t)j;
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
		int j;

		for (j = 0; j < 10; j++)
		{
			for (i = 0; i < num; i++)
			{
				f64 = LehmanFactor(comp[i], 3.5, 0, 0.1);
				if ((f64 == f1[i]) || (f64 == f2[i]))
				{
					correct++;
				}
			}
		}

		gettimeofday(&gstop, NULL);
		t_time = ytools_difftime(&gstart, &gstop);

		printf("Lehman got %d of %d correct in %2.2f sec\n", correct, num, t_time);
		printf("percent correct = %.2f\n", 100.0*(double)correct / (double)num);
		printf("average time per input = %1.4f ms\n", 1000 * t_time / (double)num);
	}

	// spbrent test
    // msvc max of 62 bits
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
            if (fscanf(in, "%" PRIu64 ",%u,%u", comp + i, f1 + i, f2 + i) != 3) break;
            mpz_set_u64(gmptmp, comp[i]);
            j = (int)mpz_sizeinbase(gmptmp, 2);
            totBits += (uint32_t)j;
            if ((uint32_t)j > maxBits)
                maxBits = (uint32_t)j;
            if ((uint32_t)j < minBits && j != 0)
                minBits = (uint32_t)j;
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
        t_time = ytools_difftime(&gstart, &gstop);

        printf("rho got %d of %d correct in %2.2f sec (%d by squfof)\n", correct, num, t_time, k);
        printf("percent correct = %.2f\n", 100.0*(double)correct / (double)num);
        printf("average time per input = %1.4f ms\n", 1000 * t_time / (double)num);
    }

	// sequential squfof test
	for (nf = 0; nf < 5; nf++)
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
			if (fscanf(in, "%" PRIu64 ",%u,%u", comp + i, f1 + i, f2 + i) != 3) break;
            mpz_set_u64(gmptmp, comp[i]);
            j = (int)mpz_sizeinbase(gmptmp, 2);
			totBits += (uint32_t)j;
			if ((uint32_t)j > maxBits)
				maxBits = (uint32_t)j;
			if ((uint32_t)j < minBits && j != 0)
				minBits = (uint32_t)j;
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

			if (((uint32_t)f64 == f1[i]) || ((uint32_t)f64 == f2[i]))
				correct++;
			else
			{
				int p;
				for (p = 0; p < 3; p++)
				{
					f64 = spbrent(comp[i], (uint64_t)(p + 1), 8192);
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
		t_time = ytools_difftime(&gstart, &gstop);

		printf("squfof got %d of %d correct in %2.2f sec (%d by rho)\n", correct, num, t_time, k);
		printf("percent correct = %.2f\n", 100.0*(double)correct / (double)num);
		printf("average time per input = %1.4f ms\n", 1000 * t_time / (double)num);
	}

	goto done;

tinyqs_marker:

    i = 0;
	strcpy(filenames[i++], "semiprimes_32bit.dat");
	strcpy(filenames[i++], "semiprimes_40bit.dat");
	strcpy(filenames[i++], "semiprimes_46bit.dat");
	strcpy(filenames[i++], "semiprimes_50bit.dat");
	strcpy(filenames[i++], "semiprimes_52bit.dat");
	strcpy(filenames[i++], "semiprimes_54bit.dat");
    strcpy(filenames[i++], "semiprimes_56bit.dat");
	strcpy(filenames[i++], "semiprimes_58bit.dat");
    strcpy(filenames[i++], "semiprimes_60bit.dat");
    strcpy(filenames[i++], "semiprimes_64bit.dat");
    //strcpy(filenames[i++], "semiprimes_70bit.dat");
    //strcpy(filenames[i++], "semiprimes_80bit.dat");
    //strcpy(filenames[i++], "semiprimes_90bit.dat");
    //strcpy(filenames[i++], "semiprimes_100bit.dat");
    //strcpy(filenames[i++], "semiprimes_110bit.dat");
	//strcpy(filenames[i++], "semiprimes_120bit.dat");
    //strcpy(filenames[i++], "pseudoprimes_125bit.dat");
    num_files = i;
    num = 100000;

    // tinyqs test
    {
        mpz_t gmp_comp;
        mpz_init(gmp_comp);

        for (nf = 0; nf < 0; nf++)
        {
            uint64_t known1, known2;
            char buf[1024];
            in = fopen(filenames[nf], "r");

			if (in == NULL)
			{
				printf("couldn't find %s to open\n", filenames[nf]);
				continue;
			}

			//num = 1000;
			totBits = 0;
			minBits = 99999;
			maxBits = 0;

			mpz_t* inputs = (mpz_t*)xmalloc((unsigned long)num * sizeof(mpz_t));

			for (i = 0; i < num; i++)
			{
				if (fgets(buf, 1024, in) == NULL) break;
#ifdef _MSC_VER
				gmp_sscanf(buf, "%Zd, %"PRIu64", %"PRIu64"",
					gmp_comp, &known1, &known2);
#else
				gmp_sscanf(buf, "%Zd, %" PRIu64 ", %" PRIu64 "",
					gmp_comp, &known1, &known2);
#endif

				bits = (uint32_t)mpz_sizeinbase(gmp_comp, 2);
				totBits += bits;
				if (bits < minBits) minBits = bits;
				if (bits > maxBits) maxBits = bits;

				mpz_init(inputs[i]);
				mpz_set(inputs[i], gmp_comp);
			}

			fclose(in);

			printf("average bits of input numbers: %1.2f\n", (double)totBits / (double)num);
			printf("minimum bits of input numbers: %u\n", minBits);
			printf("maximum bits of input numbers: %u\n", maxBits);
			printf("commencing mpqs test\n");

			test_tinysiqs(inputs, num);

			for (i = 0; i < num; i++)
			{
				mpz_clear(inputs[i]);
			}
        }

        mpz_clear(gmp_comp);
    }   

	goto done;

kmpqs_marker:

	i = 0;
	strcpy(filenames[i++], "semiprimes_32bit.dat");
	strcpy(filenames[i++], "semiprimes_40bit.dat");
	strcpy(filenames[i++], "semiprimes_46bit.dat");
	strcpy(filenames[i++], "semiprimes_50bit.dat");
	strcpy(filenames[i++], "semiprimes_52bit.dat");
	strcpy(filenames[i++], "semiprimes_54bit.dat");
	strcpy(filenames[i++], "semiprimes_56bit.dat");
	strcpy(filenames[i++], "semiprimes_58bit.dat");
	strcpy(filenames[i++], "semiprimes_60bit.dat");
	strcpy(filenames[i++], "semiprimes_64bit.dat");
	//strcpy(filenames[i++], "semiprimes_70bit.dat");
	//strcpy(filenames[i++], "semiprimes_80bit.dat");
	//strcpy(filenames[i++], "semiprimes_90bit.dat");
	//strcpy(filenames[i++], "semiprimes_100bit.dat");
	//strcpy(filenames[i++], "semiprimes_110bit.dat");
	//strcpy(filenames[i++], "semiprimes_120bit.dat");
	//strcpy(filenames[i++], "pseudoprimes_125bit.dat");
	num_files = i;
	num = 100000;

	// kmpqs test
	{
		mpz_t gmp_comp;
		mpz_init(gmp_comp);

		for (nf = 0; nf < num_files; nf++)
		{
			uint64_t known1, known2;
			char buf[1024];
			in = fopen(filenames[nf], "r");

			if (in == NULL)
			{
				printf("couldn't find %s to open\n", filenames[nf]);
				continue;
			}

			//num = 10000;
			totBits = 0;
			minBits = 99999;
			maxBits = 0;

			mpz_t* inputs = (mpz_t*)xmalloc((unsigned long)num * sizeof(mpz_t));

			for (i = 0; i < num; i++)
			{
				if (fgets(buf, 1024, in) == NULL) break;
#ifdef _MSC_VER
				gmp_sscanf(buf, "%Zd, %"PRIu64", %"PRIu64"",
					gmp_comp, &known1, &known2);
#else
				gmp_sscanf(buf, "%Zd, %" PRIu64 ", %" PRIu64 "",
					gmp_comp, &known1, &known2);
#endif

				bits = (uint32_t)mpz_sizeinbase(gmp_comp, 2);
				totBits += bits;
				if (bits < minBits) minBits = bits;
				if (bits > maxBits) maxBits = bits;

				mpz_init(inputs[i]);
				mpz_set(inputs[i], gmp_comp);
			}

			fclose(in);

			printf("average bits of input numbers: %1.2f\n", (double)totBits / (double)num);
			printf("minimum bits of input numbers: %u\n", minBits);
			printf("maximum bits of input numbers: %u\n", maxBits);
			printf("commencing mpqs test\n");

			test_kmpqs(inputs, num);

			for (i = 0; i < num; i++)
			{
				mpz_clear(inputs[i]);
			}
		}

		mpz_clear(gmp_comp);
	}

	goto done;

tinyecm_104_list_marker:
	i = 0;
	//strcpy(filenames[i++], "semiprimes_tlp_32x32x32.txt");
	//strcpy(filenames[i++], "semiprimes_tlp_32x64.txt");
	//strcpy(filenames[i++], "semiprimes_tlp_48x48.txt");
	strcpy(filenames[i++], "semiprimes_50bit.dat");
	strcpy(filenames[i++], "semiprimes_60bit.dat");
	strcpy(filenames[i++], "semiprimes_64bit.dat");
    strcpy(filenames[i++], "semiprimes_70bit.dat");
	strcpy(filenames[i++], "semiprimes_80bit.dat");
	strcpy(filenames[i++], "semiprimes_90bit.dat");
	strcpy(filenames[i++], "semiprimes_100bit.dat");
	num_files = i;

	// tinyecm_104_list test
	for (nf = 0; nf < num_files; nf++)
	{
		mpz_t gmp_comp, gmp_f;
		char buf[1024];
		UNUSED_VAR int curves;
		UNUSED_VAR int B1;
		UNUSED_VAR uint64_t known1, known2, known3;

		in = fopen(filenames[nf], "r");
		printf("testing file: %s\n", filenames[nf]);

		mpz_init(gmp_comp);
		mpz_init(gmp_f);

		gettimeofday(&gstart, NULL);

		num = 10000;
		totBits = 0;
		minBits = 99999;
		maxBits = 0;
		for (i = 0; i < num; i++)
		{
			if (fgets(buf, 1024, in) == NULL) break;

			gmp_sscanf(buf, "%Zd, %" PRIu64 ", %" PRIu64 "",
				gmp_comp, &known1, &known2);

			bits = (uint32_t)mpz_sizeinbase(gmp_comp, 2);
			totBits += bits;
			if (bits < minBits) minBits = bits;
			if (bits > maxBits) maxBits = bits;

			comp[2 * i + 0] = mpz_get_ui(gmp_comp) & 0xffffffff;
			mpz_tdiv_q_2exp(gmp_comp, gmp_comp, 32);
			comp[2 * i + 0] |= (uint64_t)(mpz_get_ui(gmp_comp) & 0xfffff) << 32;
			mpz_tdiv_q_2exp(gmp_comp, gmp_comp, 20);
			comp[2 * i + 1] = mpz_get_ui(gmp_comp) & 0xffffffff;
			mpz_tdiv_q_2exp(gmp_comp, gmp_comp, 32);
			comp[2 * i + 1] |= (uint64_t)(mpz_get_ui(gmp_comp) & 0xfffff) << 32;

			f64a[i] = known1;
			f64b[i] = known2;
		}

		printf("average bits of input numbers: %1.2f\n", (double)totBits / (double)num);
		printf("minimum bits of input numbers: %u\n", minBits);
		printf("maximum bits of input numbers: %u\n", maxBits);
		printf("commencing getfactor_tecm_x8_list test\n");

		getfactor_tecm_x8_list(comp, outf, (int)(maxBits / 2), (uint32_t)num, &lcg);

		correct = 0;
		for (j = 0; j < num; j++)
		{
			if ((outf[j * 2 + 0] == f64a[j]) ||
				(outf[j * 2 + 0] == f64b[j]))
			{
				correct++;
			}

		}

		gettimeofday(&gstop, NULL);
		t_time = ytools_difftime(&gstart, &gstop);

		printf("tinyecm got %d of %d correct in %2.2f sec\n", correct, num, t_time);
		printf("percent correct = %.2f\n", 100.0 * (double)correct / (double)num);
		printf("average time per input = %.2f ms\n", 1000 * t_time / (double)num);

		mpz_clear(gmp_comp);
		mpz_clear(gmp_f);	
	}

	goto done;

tinyecm_104_edw_list_marker:
	i = 0;
	//strcpy(filenames[i++], "semiprimes_tlp_32x32x32.txt");
	//strcpy(filenames[i++], "semiprimes_tlp_32x64.txt");
	//strcpy(filenames[i++], "semiprimes_tlp_48x48.txt");
	strcpy(filenames[i++], "semiprimes_50bit.dat");
	strcpy(filenames[i++], "semiprimes_60bit.dat");
	strcpy(filenames[i++], "semiprimes_64bit.dat");
	strcpy(filenames[i++], "semiprimes_70bit.dat");
	strcpy(filenames[i++], "semiprimes_80bit.dat");
	strcpy(filenames[i++], "semiprimes_90bit.dat");
	strcpy(filenames[i++], "semiprimes_100bit.dat");
	num_files = i;

	// tinyecm_104_list test
	for (nf = 0; nf < num_files; nf++)
	{
		mpz_t gmp_comp, gmp_f;
		char buf[1024];
		UNUSED_VAR int curves;
		UNUSED_VAR int B1;
		UNUSED_VAR uint64_t known1, known2, known3;

		in = fopen(filenames[nf], "r");
		printf("testing file: %s\n", filenames[nf]);

		mpz_init(gmp_comp);
		mpz_init(gmp_f);

		gettimeofday(&gstart, NULL);

		num = 10000;
		totBits = 0;
		minBits = 99999;
		maxBits = 0;
		for (i = 0; i < num; i++)
		{
			if (fgets(buf, 1024, in) == NULL) break;
			gmp_sscanf(buf, "%Zd, %"PRIu64", %"PRIu64"",
				gmp_comp, &known1, &known2);

			bits = (uint32_t)mpz_sizeinbase(gmp_comp, 2);
			totBits += bits;
			if (bits < minBits) minBits = bits;
			if (bits > maxBits) maxBits = bits;

			comp[2 * i + 0] = mpz_get_ui(gmp_comp) & 0xffffffff;
			mpz_tdiv_q_2exp(gmp_comp, gmp_comp, 32);
			comp[2 * i + 0] |= (uint64_t)(mpz_get_ui(gmp_comp) & 0xfffff) << 32;
			mpz_tdiv_q_2exp(gmp_comp, gmp_comp, 20);
			comp[2 * i + 1] = mpz_get_ui(gmp_comp) & 0xffffffff;
			mpz_tdiv_q_2exp(gmp_comp, gmp_comp, 32);
			comp[2 * i + 1] |= (uint64_t)(mpz_get_ui(gmp_comp) & 0xfffff) << 32;

			f64a[i] = known1;
			f64b[i] = known2;
		}

		printf("average bits of input numbers: %1.2f\n", (double)totBits / (double)num);
		printf("minimum bits of input numbers: %u\n", minBits);
		printf("maximum bits of input numbers: %u\n", maxBits);
		printf("commencing getfactor_tecm_x8_list test\n");

		getfactor_tecm_edw_x8_list(comp, outf, (int)(maxBits / 2), (uint32_t)num, &lcg);

		correct = 0;
		for (j = 0; j < num; j++)
		{
			if ((outf[j * 2 + 0] == f64a[j]) ||
				(outf[j * 2 + 0] == f64b[j]))
			{
				correct++;
			}

		}

		gettimeofday(&gstop, NULL);
		t_time = ytools_difftime(&gstart, &gstop);

		printf("tinyecm got %d of %d correct in %2.2f sec\n", correct, num, t_time);
		printf("percent correct = %.2f\n", 100.0 * (double)correct / (double)num);
		printf("average time per input = %.2f ms\n", 1000 * t_time / (double)num);

		mpz_clear(gmp_comp);
		mpz_clear(gmp_f);
	}

	goto done;

tinyecm_128_marker:
	i = 0;
	//strcpy(filenames[i++], "semiprimes_tlp_32x32x32.txt");
	//strcpy(filenames[i++], "semiprimes_tlp_32x64.txt");
	//strcpy(filenames[i++], "semiprimes_tlp_48x48.txt");
	strcpy(filenames[i++], "semiprimes_50bit.dat");
	strcpy(filenames[i++], "semiprimes_60bit.dat");
	strcpy(filenames[i++], "semiprimes_64bit.dat");
	strcpy(filenames[i++], "semiprimes_70bit.dat");
	strcpy(filenames[i++], "semiprimes_80bit.dat");
	strcpy(filenames[i++], "semiprimes_90bit.dat");
	strcpy(filenames[i++], "semiprimes_100bit.dat");
	num_files = i;

	// tinyecm_128 test
	for (nf = 0; nf < num_files; nf++)
	{
		mpz_t gmp_comp, gmp_f;
		char buf[1024];
		UNUSED_VAR int curves;
		UNUSED_VAR int B1;
		uint64_t known1, known2, known3;

		in = fopen(filenames[nf], "r");
		printf("testing file: %s\n", filenames[nf]);

		mpz_init(gmp_comp);
		mpz_init(gmp_f);

		gettimeofday(&gstart, NULL);

		num = 10000;
		totBits = 0;
		minBits = 99999;
		maxBits = 0;
		correct = 0;
		printf("commencing getfactor_tecm test\n");
		for (i = 0; i < num; i++)
		{
			if (fgets(buf, 1024, in) == NULL) break;
			gmp_sscanf(buf, "%Zd, %"PRIu64", %"PRIu64"",
				gmp_comp, &known1, &known2);

			bits = (uint32_t)mpz_sizeinbase(gmp_comp, 2);
			totBits += bits;
			if (bits < minBits) minBits = bits;
			if (bits > maxBits) maxBits = bits;

			getfactor_tecm(gmp_comp, gmp_f, 0, &lcg);

			known3 = mpz_get_ui(gmp_f) & 0xffffffff;
			mpz_tdiv_q_2exp(gmp_f, gmp_f, 32);
			known3 |= (uint64_t)mpz_get_ui(gmp_f) << 32;

			if ((known3 == known1) ||
				(known3 == known2))
			{
				correct++;
			}
		}

		printf("average bits of input numbers: %1.2f\n", (double)totBits / (double)num);
		printf("minimum bits of input numbers: %u\n", minBits);
		printf("maximum bits of input numbers: %u\n", maxBits);
		
		gettimeofday(&gstop, NULL);
		t_time = ytools_difftime(&gstart, &gstop);

		printf("tinyecm got %d of %d correct in %2.2f sec\n", correct, num, t_time);
		printf("percent correct = %.2f\n", 100.0 * (double)correct / (double)num);
		printf("average time per input = %.2f ms\n", 1000 * t_time / (double)num);

		mpz_clear(gmp_comp);
		mpz_clear(gmp_f);
	}

	goto done;

tinyecm_128_edw_marker:
	i = 0;
	//strcpy(filenames[i++], "semiprimes_tlp_32x32x32.txt");
	//strcpy(filenames[i++], "semiprimes_tlp_32x64.txt");
	//strcpy(filenames[i++], "semiprimes_tlp_48x48.txt");
	strcpy(filenames[i++], "semiprimes_50bit.dat");
	strcpy(filenames[i++], "semiprimes_60bit.dat");
	strcpy(filenames[i++], "semiprimes_64bit.dat");
	strcpy(filenames[i++], "semiprimes_70bit.dat");
	strcpy(filenames[i++], "semiprimes_80bit.dat");
	strcpy(filenames[i++], "semiprimes_90bit.dat");
	strcpy(filenames[i++], "semiprimes_100bit.dat");
	num_files = i;

	// tinyecm_128 test
	for (nf = 0; nf < num_files; nf++)
	{
		mpz_t gmp_comp, gmp_f;
		char buf[1024];
		UNUSED_VAR int curves;
		UNUSED_VAR int B1;
		uint64_t known1, known2, known3;

		in = fopen(filenames[nf], "r");
		printf("testing file: %s\n", filenames[nf]);

		mpz_init(gmp_comp);
		mpz_init(gmp_f);

		gettimeofday(&gstart, NULL);

		num = 10000;
		totBits = 0;
		minBits = 99999;
		maxBits = 0;
		correct = 0;
		printf("commencing getfactor_tecm test\n");
		for (i = 0; i < num; i++)
		{
			if (fgets(buf, 1024, in) == NULL) break;
			gmp_sscanf(buf, "%Zd, %"PRIu64", %"PRIu64"",
				gmp_comp, &known1, &known2);

			bits = (uint32_t)mpz_sizeinbase(gmp_comp, 2);
			totBits += bits;
			if (bits < minBits) minBits = bits;
			if (bits > maxBits) maxBits = bits;

			getfactor_tecm_edw(gmp_comp, gmp_f, 0, &lcg);

			known3 = mpz_get_ui(gmp_f) & 0xffffffff;
			mpz_tdiv_q_2exp(gmp_f, gmp_f, 32);
			known3 |= (uint64_t)mpz_get_ui(gmp_f) << 32;

			if ((known3 == known1) ||
				(known3 == known2))
			{
				correct++;
			}
		}

		printf("average bits of input numbers: %1.2f\n", (double)totBits / (double)num);
		printf("minimum bits of input numbers: %u\n", minBits);
		printf("maximum bits of input numbers: %u\n", maxBits);

		gettimeofday(&gstop, NULL);
		t_time = ytools_difftime(&gstart, &gstop);

		printf("tinyecm got %d of %d correct in %2.2f sec\n", correct, num, t_time);
		printf("percent correct = %.2f\n", 100.0 * (double)correct / (double)num);
		printf("average time per input = %.2f ms\n", 1000 * t_time / (double)num);

		mpz_clear(gmp_comp);
		mpz_clear(gmp_f);
	}

	goto done;

uecm_52_list_marker:
	i = 0;
	strcpy(filenames[i++], "semiprimes_32bit.dat");
	strcpy(filenames[i++], "semiprimes_34bit.dat");
	strcpy(filenames[i++], "semiprimes_36bit.dat");
	strcpy(filenames[i++], "semiprimes_38bit.dat");
	strcpy(filenames[i++], "semiprimes_40bit.dat");
	strcpy(filenames[i++], "semiprimes_42bit.dat");
	strcpy(filenames[i++], "semiprimes_44bit.dat");
	strcpy(filenames[i++], "semiprimes_46bit.dat");
	strcpy(filenames[i++], "semiprimes_48bit.dat");
	strcpy(filenames[i++], "semiprimes_50bit.dat");
	strcpy(filenames[i++], "semiprimes_52bit.dat");
	num_files = i;

	// uecm_52_list test
	for (nf = 0; nf < num_files; nf++)
	{
		mpz_t gmp_comp, gmp_f;
		char buf[1024];
		UNUSED_VAR int curves;
		UNUSED_VAR int B1;
		UNUSED_VAR uint64_t known1, known2, known3;

		in = fopen(filenames[nf], "r");
		printf("testing file: %s\n", filenames[nf]);

		mpz_init(gmp_comp);
		mpz_init(gmp_f);

		gettimeofday(&gstart, NULL);

		num = 100000;
		totBits = 0;
		minBits = 99999;
		maxBits = 0;
		for (i = 0; i < num; i++)
		{
			if (fgets(buf, 1024, in) == NULL) break;
#ifdef _MSC_VER
			gmp_sscanf(buf, "%Zd, %u, %u",
				gmp_comp, &known1, &known2);
#else
			gmp_sscanf(buf, "%Zd, %u, %u",
				gmp_comp, &known1, &known2);
#endif

			bits = (uint32_t)mpz_sizeinbase(gmp_comp, 2);
			totBits += bits;
			if (bits < minBits) minBits = bits;
			if (bits > maxBits) maxBits = bits;

			comp[i] = mpz_get_ui(gmp_comp) & 0xfffffffffffffull;

			f1[i] = (uint32_t)known1;
			f2[i] = (uint32_t)known2;
		}

		printf("average bits of input numbers: %1.2f\n", (double)totBits / (double)num);
		printf("minimum bits of input numbers: %u\n", minBits);
		printf("maximum bits of input numbers: %u\n", maxBits);
		printf("commencing getfactor_uecm_x8_list test\n");

		getfactor_uecm_x8_list(comp, outf, (uint32_t)num, &lcg);

		correct = 0;
		for (j = 0; j < num; j++)
		{
			if ((outf[j] == f1[j]) ||
				(outf[j] == f2[j]))
			{
				correct++;
			}

		}

		gettimeofday(&gstop, NULL);
		t_time = ytools_difftime(&gstart, &gstop);

		printf("uecm got %d of %d correct in %2.2f sec\n", correct, num, t_time);
		printf("percent correct = %.2f\n", 100.0 * (double)correct / (double)num);
		printf("average time per input = %.2f ms\n", 1000 * t_time / (double)num);

		mpz_clear(gmp_comp);
		mpz_clear(gmp_f);
	}

	goto done;

microecm_marker:
	i = 0;
	strcpy(filenames[i++], "semiprimes_32bit.dat");
	strcpy(filenames[i++], "semiprimes_34bit.dat");
	strcpy(filenames[i++], "semiprimes_36bit.dat");
	strcpy(filenames[i++], "semiprimes_38bit.dat");
	strcpy(filenames[i++], "semiprimes_40bit.dat");
	strcpy(filenames[i++], "semiprimes_42bit.dat");		// 70
	strcpy(filenames[i++], "semiprimes_44bit.dat");		// 70
	strcpy(filenames[i++], "semiprimes_46bit.dat");		// 70
	strcpy(filenames[i++], "semiprimes_48bit.dat");		// 70
	strcpy(filenames[i++], "semiprimes_50bit.dat");		// 70
	strcpy(filenames[i++], "semiprimes_52bit.dat");		// 85
	strcpy(filenames[i++], "semiprimes_54bit.dat");       // 85
	strcpy(filenames[i++], "semiprimes_56bit.dat");		// 125
	strcpy(filenames[i++], "semiprimes_58bit.dat");		// 125
	strcpy(filenames[i++], "semiprimes_60bit.dat");		// 165
	strcpy(filenames[i++], "semiprimes_62bit.dat");		// 165
	strcpy(filenames[i++], "semiprimes_64bit.dat");		// 205

	num_files = i;
	num = 100000;

	lcg = 42;

	// uecm test
	for (nf = 0; nf < num_files; nf++)
	{
		mpz_t gmp_comp, gmp_f;
		UNUSED_VAR char buf[1024];
		UNUSED_VAR int curves;
		UNUSED_VAR int B1;
		UNUSED_VAR uint64_t known1, known2, known3;

		in = fopen(filenames[nf], "r");
		printf("testing file: %s\n", filenames[nf]);

		mpz_init(gmp_comp);
		mpz_init(gmp_f);

		gettimeofday(&gstart, NULL);

		correct = 0;
		k = 0;

		totBits = 0;
		minBits = 999;
		maxBits = 0;
		i = 0;
		while (!feof(in))
		{
			if (fscanf(in, "%" PRIu64 ",%u,%u", comp + i, f1 + i, f2 + i) != 3) break;
			mpz_set_u64(gmptmp, comp[i]);

			j = (int)mpz_sizeinbase(gmptmp, 2);
			totBits += (uint32_t)j;
			if ((uint32_t)j > maxBits)
				maxBits = (uint32_t)j;
			if ((uint32_t)j < minBits && j != 0)
				minBits = (uint32_t)j;
			i++;
		}
		printf("average bits of input numbers = %.2f\n", (double)totBits / (double)i);
		printf("minimum bits of input numbers = %d\n", minBits);
		printf("maximum bits of input numbers = %d\n", maxBits);

		fclose(in);

		for (i = 0; i < num; i++)
		{
			uint64_t result = getfactor_uecm(comp[i], 0, &lcg);

			if ((result == f1[i]) || (result == f2[i]))
			{
				correct++;
			}
		}

		gettimeofday(&gstop, NULL);
		t_time = ytools_difftime(&gstart, &gstop);

		printf("microecm got %d of %d correct in %2.2f sec\n",
			correct, num, t_time);
		printf("percent correct = %.2f\n", 100.0 * (double)correct / (double)num);
		printf("average time per input = %.2f ms\n", 1000 * t_time / (double)num);

		mpz_clear(gmp_comp);
		mpz_clear(gmp_f);
	}
	goto done;

microecm_edw_marker:
	i = 0;
	strcpy(filenames[i++], "semiprimes_32bit.dat");
	strcpy(filenames[i++], "semiprimes_34bit.dat");
	strcpy(filenames[i++], "semiprimes_36bit.dat");
	strcpy(filenames[i++], "semiprimes_38bit.dat");
	strcpy(filenames[i++], "semiprimes_40bit.dat");
	strcpy(filenames[i++], "semiprimes_42bit.dat");		// 70
	strcpy(filenames[i++], "semiprimes_44bit.dat");		// 70
	strcpy(filenames[i++], "semiprimes_46bit.dat");		// 70
	strcpy(filenames[i++], "semiprimes_48bit.dat");		// 70
	strcpy(filenames[i++], "semiprimes_50bit.dat");		// 70
	strcpy(filenames[i++], "semiprimes_52bit.dat");		// 85
	strcpy(filenames[i++], "semiprimes_54bit.dat");       // 85
	strcpy(filenames[i++], "semiprimes_56bit.dat");		// 125
	strcpy(filenames[i++], "semiprimes_58bit.dat");		// 125
	strcpy(filenames[i++], "semiprimes_60bit.dat");		// 165
	strcpy(filenames[i++], "semiprimes_62bit.dat");		// 165
	strcpy(filenames[i++], "semiprimes_64bit.dat");		// 205

	num_files = i;
	num = 100000;

	lcg = 42;

	// uecm test
	for (nf = 0; nf < num_files; nf++)
	{
		mpz_t gmp_comp, gmp_f;
		UNUSED_VAR char buf[1024];
		UNUSED_VAR int curves;
		UNUSED_VAR int B1;
		UNUSED_VAR uint64_t known1, known2, known3;

		in = fopen(filenames[nf], "r");
		printf("testing file: %s\n", filenames[nf]);

		mpz_init(gmp_comp);
		mpz_init(gmp_f);

		gettimeofday(&gstart, NULL);

		correct = 0;
		k = 0;

		totBits = 0;
		minBits = 999;
		maxBits = 0;
		i = 0;
		while (!feof(in))
		{
			if (fscanf(in, "%" PRIu64 ",%u,%u", comp + i, f1 + i, f2 + i) != 3) break;
			mpz_set_u64(gmptmp, comp[i]);

			j = (int)mpz_sizeinbase(gmptmp, 2);
			totBits += (uint32_t)j;
			if ((uint32_t)j > maxBits)
				maxBits = (uint32_t)j;
			if ((uint32_t)j < minBits && j != 0)
				minBits = (uint32_t)j;
			i++;
		}
		printf("average bits of input numbers = %.2f\n", (double)totBits / (double)i);
		printf("minimum bits of input numbers = %d\n", minBits);
		printf("maximum bits of input numbers = %d\n", maxBits);

		fclose(in);

		for (i = 0; i < num; i++)
		{
			uint64_t result = getfactor_uecm_edw(comp[i], 0, &lcg);

			if ((result == f1[i]) || (result == f2[i]))
			{
				correct++;
			}
		}

		gettimeofday(&gstop, NULL);
		t_time = ytools_difftime(&gstart, &gstop);

		printf("microecm got %d of %d correct in %2.2f sec\n",
			correct, num, t_time);
		printf("percent correct = %.2f\n", 100.0 * (double)correct / (double)num);
		printf("average time per input = %.2f ms\n", 1000 * t_time / (double)num);

		mpz_clear(gmp_comp);
		mpz_clear(gmp_f);
	}
	goto done;

mpqs_marker:
	i = 0;
	strcpy(filenames[i++], "semiprimes_70bit.dat");
	strcpy(filenames[i++], "semiprimes_80bit.dat");
	strcpy(filenames[i++], "semiprimes_90bit.dat");
	strcpy(filenames[i++], "semiprimes_100bit.dat");
	num_files = i;

	for (nf = 0; nf < num_files; nf++)
	{
		mpz_t gmp_comp;
		char buf[1024];
		UNUSED_VAR int curves;
		UNUSED_VAR int B1;
		UNUSED_VAR uint64_t known1, known2, known3;

		in = fopen(filenames[nf], "r");
		printf("testing file: %s\n", filenames[nf]);

		mpz_init(gmp_comp);

		gettimeofday(&gstart, NULL);

		num = 10000;
		totBits = 0;
		minBits = 99999;
		maxBits = 0;

		mpz_t* inputs = (mpz_t*)xmalloc((unsigned long)num * sizeof(mpz_t));

		for (i = 0; i < num; i++)
		{
			if (fgets(buf, 1024, in) == NULL) break;
#ifdef _MSC_VER
			gmp_sscanf(buf, "%Zd, %"PRIu64", %"PRIu64"",
				gmp_comp, &known1, &known2);
#else
			gmp_sscanf(buf, "%Zd, %" PRIu64 ", %" PRIu64 "",
				gmp_comp, &known1, &known2);
#endif

			bits = (uint32_t)mpz_sizeinbase(gmp_comp, 2);
			totBits += bits;
			if (bits < minBits) minBits = bits;
			if (bits > maxBits) maxBits = bits;

			mpz_init(inputs[i]);
			mpz_set(inputs[i], gmp_comp);
		}

		printf("average bits of input numbers: %1.2f\n", (double)totBits / (double)num);
		printf("minimum bits of input numbers: %u\n", minBits);
		printf("maximum bits of input numbers: %u\n", maxBits);
		printf("commencing mpqs test\n");

		test_smallmpqs(inputs, num);

		mpz_clear(gmp_comp);
		for (i = 0; i < num; i++)
		{
			mpz_clear(inputs[i]);
		}
	}

	goto done;

done:
	
	return;

	free(outf);
	free(f1);
	free(f2);
	free(f64a);
	free(f64b);
	free(comp);
	mpz_clear(gmptmp);
	return;
}
