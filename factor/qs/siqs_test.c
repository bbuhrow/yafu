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

#include "qs.h"
#include "qs_impl.h"
#include "ytools.h"
#include "tinyecm.h"

int check_relation(mpz_t a, mpz_t b, siqs_r* r, fb_list* fb, mpz_t n, int VFLAG, int knmod8)
{
	int offset, parity, num_factors;
	int j, retval;
	mpz_t Q, RHS;
	// unsigned!
	uint32_t lp[4];

	mpz_init(Q);
	mpz_init(RHS);

	offset = r->sieve_offset;
	lp[0] = r->large_prime[0];
	lp[1] = r->large_prime[1];
	lp[2] = r->large_prime[2];
	lp[3] = r->large_prime[3];
	parity = r->parity;
	num_factors = r->num_factors;

	mpz_set_ui(RHS, lp[0]);
	mpz_mul_ui(RHS, RHS, lp[1]);
	mpz_mul_ui(RHS, RHS, lp[2]);
	mpz_mul_ui(RHS, RHS, lp[3]);
	for (j = 0; j < num_factors; j++)
	{
		if (r->fb_offsets[j] > fb->B)
		{
			printf("error fb_offset is greater than the factor base size (%u > %u)\n",
				r->fb_offsets[j], fb->B);
		}
		mpz_mul_ui(RHS, RHS, fb->list->prime[r->fb_offsets[j]]);
	}

	if (knmod8 == 1)
	{
		// Q(x) = (2ax + b)^2 - N, where x is the sieve index
		mpz_mul_ui(Q, a, offset);
		mpz_mul_2exp(Q, Q, 1);
	}
	else
	{
		//Q(x)/a = (ax + b)^2 - N, where x is the sieve index
		mpz_mul_ui(Q, a, offset);
	}

	if (parity)
		mpz_sub(Q, Q, b);
	else
		mpz_add(Q, Q, b);
	mpz_mul(Q, Q, Q);
	mpz_sub(Q, Q, n);

	retval = 0;
	if (mpz_sgn(Q) < 0)
		mpz_neg(Q, Q);

	if (mpz_cmp(Q, RHS) != 0)
	{
		if (VFLAG > 1)
		{
			printf("error Q != RHS\n");
		}

		if (VFLAG > 0)
		{
			printf("error Q != RHS\n");
			gmp_printf("Q = %Zd, RHS = %Zd\n", Q, RHS);
			gmp_printf("A = %Zd\nB = %Zd\n", a, b);
			printf("offset = %u, parity = %d\n", offset, parity);
			printf("fb_offsets:primes:\n");
			for (j = 0; j < num_factors; j++)
				printf("%u:%u ", r->fb_offsets[j], fb->list->prime[r->fb_offsets[j]]);
			printf("\n");
			printf("lp: %u, %u, %u, %u\n", lp[0], lp[1], lp[2], lp[3]);
			printf("lp: %u, %u, %u, %u\n", r->large_prime[0], r->large_prime[1], r->large_prime[2], r->large_prime[3]);
		}
		retval = 1;
	}

	mpz_clear(Q);
	mpz_clear(RHS);

	return retval;
}

int check_specialcase(FILE* sieve_log, fact_obj_t* fobj)
{
	// check for some special cases of input number
	// sieve_log is passed in already open, and should return open
	if (mpz_even_p(fobj->qs_obj.gmp_n))
	{
		printf("input must be odd\n");
		return 1;
	}

	if (is_mpz_prp(fobj->qs_obj.gmp_n, fobj->NUM_WITNESSES))
	{
		char* s = mpz_get_str(NULL, 10, fobj->qs_obj.gmp_n);
		add_to_factor_list(fobj->factors, fobj->qs_obj.gmp_n,
			fobj->VFLAG, fobj->NUM_WITNESSES, 0);
		if (sieve_log != NULL)
			logprint(sieve_log, "prp%d = %s\n", gmp_base10(fobj->qs_obj.gmp_n), s);
		mpz_set_ui(fobj->qs_obj.gmp_n, 1);
		free(s);
		return 1;
	}

	if (mpz_perfect_square_p(fobj->qs_obj.gmp_n))
	{
		mpz_sqrt(fobj->qs_obj.gmp_n, fobj->qs_obj.gmp_n);

		char* s = mpz_get_str(NULL, 10, fobj->qs_obj.gmp_n);
		add_to_factor_list(fobj->factors, fobj->qs_obj.gmp_n,
			fobj->VFLAG, fobj->NUM_WITNESSES, 0);
		if (sieve_log != NULL)
			logprint(sieve_log, "prp%d = %s\n", gmp_base10(fobj->qs_obj.gmp_n), s);
		add_to_factor_list(fobj->factors, fobj->qs_obj.gmp_n,
			fobj->VFLAG, fobj->NUM_WITNESSES, 0);
		if (sieve_log != NULL)
			logprint(sieve_log, "prp%d = %s\n", gmp_base10(fobj->qs_obj.gmp_n), s);
		free(s);
		mpz_set_ui(fobj->qs_obj.gmp_n, 1);
		return 1;
	}

	if (mpz_perfect_power_p(fobj->qs_obj.gmp_n))
	{
		if (fobj->VFLAG > 0)
			printf("input is a perfect power\n");

		factor_perfect_power(fobj, fobj->qs_obj.gmp_n);

		if (sieve_log != NULL)
		{
			uint32_t j;
			logprint(sieve_log, "input is a perfect power\n");

			for (j = 0; j < fobj->factors->num_factors; j++)
			{
				char* s = mpz_get_str(NULL, 10, fobj->factors->factors[j].factor);
				uint32_t k;
				for (k = 0; k < fobj->factors->factors[j].count; k++)
				{
					logprint(sieve_log, "prp%d = %s\n",
						gmp_base10(fobj->factors->factors[j].factor), s);
				}
				free(s);
			}
		}
		return 1;
	}

	if (mpz_sizeinbase(fobj->qs_obj.gmp_n, 2) < 115)
	{
		// run MPQS, as the main SIQS doesn't work for smaller inputs
		int i;
		mpz_t f1, f2, f3, f;

		mpz_init(f1);
		mpz_init(f2);
		mpz_init(f3);
		mpz_init(f);

		if (fobj->qs_obj.flags != 12345)
		{
			if (fobj->logfile != NULL)
			{
				char* s = mpz_get_str(NULL, 10, fobj->qs_obj.gmp_n);
				logprint(fobj->logfile,
					"starting tinyqs on C%d = %s\n", gmp_base10(fobj->qs_obj.gmp_n), s);
				free(s);
			}
		}

		{
			tiny_qs_params* params;
			int j;

			// todo: need to define alternate routines if this isn't defined...
			params = init_tinysiqs();
			i = tinysiqs(params, fobj->qs_obj.gmp_n, f1, f2, f3, 0);
			params = free_tinysiqs(params);

			for (j = 0; j < i; j++)
			{
				if (j == 0)
				{
					if (mpz_divisible_p(fobj->qs_obj.gmp_n, f1) == 0)
						continue;

					mpz_set(f, f1);
				}
				else if (j == 1)
				{
					if (mpz_divisible_p(fobj->qs_obj.gmp_n, f2) == 0)
						continue;

					mpz_set(f, f2);
				}
				else if (j == 2)
				{
					if (mpz_divisible_p(fobj->qs_obj.gmp_n, f3) == 0)
						continue;

					mpz_set(f, f3);
				}

				add_to_factor_list(fobj->factors, f,
					fobj->VFLAG, fobj->NUM_WITNESSES, 0);

				if (fobj->qs_obj.flags != 12345)
				{
					if (fobj->logfile != NULL)
					{
						char* s = mpz_get_str(NULL, 10, f);
						logprint(fobj->logfile,
							"prp%d = %s\n", gmp_base10(f1), s);
						free(s);
					}
				}

				mpz_tdiv_q(fobj->qs_obj.gmp_n, fobj->qs_obj.gmp_n, f);
			}

			if (mpz_cmp_ui(fobj->qs_obj.gmp_n, 1) > 0)
			{
				if (mpz_probab_prime_p(fobj->qs_obj.gmp_n, 1))
				{
					add_to_factor_list(fobj->factors, fobj->qs_obj.gmp_n,
						fobj->VFLAG, fobj->NUM_WITNESSES, 0);

					if (fobj->qs_obj.flags != 12345)
					{
						if (fobj->logfile != NULL)
						{
							char* s = mpz_get_str(NULL, 10, fobj->qs_obj.gmp_n);
							logprint(fobj->logfile,
								"prp%d = %s\n", gmp_base10(fobj->qs_obj.gmp_n), s);
							free(s);
						}
					}

					mpz_set_ui(fobj->qs_obj.gmp_n, 1);
				}
			}
		}

		if ((i == 0) || (mpz_cmp_ui(fobj->qs_obj.gmp_n, 1) > 0))
		{
			// didn't find anything or not completely factored (rare).  
			// try a different method.
			if (fobj->flags != 12345)
			{
				if (fobj->logfile != NULL)
				{
					char* s = mpz_get_str(NULL, 10, fobj->qs_obj.gmp_n);
					logprint(fobj->logfile,
						"tiny_qs failed to fully factor; "
						"starting tiny-ecm on C%d = %s\n", 
						gmp_base10(fobj->qs_obj.gmp_n), s);
					free(s);
				}
			}

			// do ECM until we find a factor.
			int found = 0;
			uint64_t lcg = 42;
			do
			{
				found = getfactor_tecm(fobj->qs_obj.gmp_n, f1, 0, &lcg);
				mpz_tdiv_r(f2, fobj->qs_obj.gmp_n, f1);

				if (found)
				{
					add_to_factor_list(fobj->factors, f1,
						fobj->VFLAG, fobj->NUM_WITNESSES, 0);

					if (fobj->qs_obj.flags != 12345)
					{
						if (fobj->logfile != NULL)
						{
							char* s = mpz_get_str(NULL, 10, f1);
							logprint(fobj->logfile,
								"prp%d = %s\n", gmp_base10(f1), s);
							free(s);
						}
					}

					mpz_tdiv_q(fobj->qs_obj.gmp_n, fobj->qs_obj.gmp_n, f1);

					if (mpz_cmp_ui(fobj->qs_obj.gmp_n, 1) > 0)
					{
						if (mpz_probab_prime_p(fobj->qs_obj.gmp_n, 1))
						{
							add_to_factor_list(fobj->factors, fobj->qs_obj.gmp_n,
								fobj->VFLAG, fobj->NUM_WITNESSES, 0);

							if (fobj->qs_obj.flags != 12345)
							{
								if (fobj->logfile != NULL)
								{
									char* s = mpz_get_str(NULL, 10, fobj->qs_obj.gmp_n);
									logprint(fobj->logfile,
										"prp%d = %s\n", gmp_base10(fobj->qs_obj.gmp_n), s);
									free(s);
								}
							}

							mpz_set_ui(fobj->qs_obj.gmp_n, 1);
						}
					}
				}
			} while (!found);


		}

		mpz_clear(f1);
		mpz_clear(f2);
		mpz_clear(f3);
		mpz_clear(f);

		return 1;
	}

	return 0;
}

int checkpoly_siqs(siqs_poly* poly, mpz_t n)
{
	//check that b^2 == N mod a
	//and that c == (b*b - n)/a
	mpz_t t1, t2, t3, t4;

	mpz_init(t1);
	mpz_init(t2);
	mpz_init(t3);
	mpz_init(t4);

	mpz_set(t1, n);
	mpz_tdiv_r(t3, t1, poly->mpz_poly_a);

	mpz_mul(t2, poly->mpz_poly_b, poly->mpz_poly_b);
	mpz_tdiv_r(t4, t2, poly->mpz_poly_a);

	if (mpz_cmp(t3, t4) != 0)
	{
		char* s1 = mpz_get_str(NULL, 10, poly->mpz_poly_b);
		char* s2 = mpz_get_str(NULL, 10, poly->mpz_poly_a);
		printf("\nError in checkpoly: %s^2 !== N mod %s\n", s1, s2);
		if (mpz_sgn(poly->mpz_poly_b) < 0)
		{
			printf("b is negative\n");
		}
		free(s1);
		free(s2);
	}

	if (mpz_kronecker(n, poly->mpz_poly_a) != 1)
		printf("\nError in checkpoly: (a|N) != 1\n");

	mpz_mul(t2, poly->mpz_poly_b, poly->mpz_poly_b);
	mpz_sub(t2, t2, n);
	mpz_tdiv_q(t4, t2, poly->mpz_poly_a);

	if (mpz_cmp(t4, poly->mpz_poly_c) != 0)
		printf("\nError in checkpoly: c != (b^2 - n)/a\n");

	mpz_clear(t1);
	mpz_clear(t2);
	mpz_clear(t3);
	mpz_clear(t4);
	return 0;
}

int checkBl(mpz_t n, uint32_t* qli, fb_list* fb, mpz_t* Bl, int s)
{
	//check that Bl^2 == N mod ql and Bl == 0 mod qj for 1<=j<=s, j != l
	int i, j, p, q;
	mpz_t t1, t2;

	mpz_init(t1);
	mpz_init(t2);

	for (i = 0; i < s; i++)
	{
		if (i == (s - 1))
			p = qli[i]; 
		else
			p = fb->list->prime[qli[i]];

		mpz_tdiv_q_2exp(t2, Bl[i], 1); //zShiftRight(&t2,&Bl[i],1);
		mpz_mul(t1, t2, t2); //zSqr(&t2,&t1);
		if (mpz_tdiv_ui(t1, p) % mpz_tdiv_ui(n, p) != 0)
		{
			char* s = mpz_get_str(NULL, 10, Bl[i]);
			printf("%s^2 mod %u != N mod %u\n", s, p, p);
			free(s);
		}

		for (j = 0; j < s; j++)
		{
			if (j == i)
				continue;

			

			if (j == (s - 1))
				q = qli[j];
			else
				q = fb->list->prime[qli[j]];

			mpz_tdiv_q_2exp(t2, Bl[i], 1); //zShiftRight(&t2,&Bl[i],1);
			if (mpz_tdiv_ui(t2, q) != 0)
			{
				char* s = mpz_get_str(NULL, 10, Bl[i]);
				printf("%s mod %u != 0\n", s, q);
				free(s);
			}
		}
	}


	mpz_clear(t1);
	mpz_clear(t2);
	return 0;
}

void siqsbench(fact_obj_t* fobj, info_t* comp_info, uint32_t digit_limit)
{
	// Run through the list of benchmark SIQS factorizations and emit
	// a Markdown results row to siqs_bench.md for use on the YAFU wiki.
	//
	// Table format (one row per run):
	//   Col 1   : System — CPU, OS, Compiler, SIMD, thread count (stacked with <br>)
	//   Col 2   : YAFU version
	//   Col 3-11: qs_time (total_time) for each benchmark digit size
	//
	// The file is opened in append mode so successive runs accumulate rows.
	// A header block is written only when the file is empty / newly created.

	// -------------------------------------------------------------------------//
	//  Benchmark inputs (digit counts: 44 50 55 60 65 70 74 80 85 90 95 100)   //
	// -------------------------------------------------------------------------//
	const char* list[] = {
		"405461849292216354219321922871108605045931309",                             /* 44 */
		"29660734457033883936073030405220515257819037444591",                        /* 50 */
		"1592739380299790264815370870751381488173344970710983383",                   /* 55 */
		"349594255864176572614071853194924838158088864370890996447417",              /* 60 */
		"34053408309992030649212497354061832056920539397279047809781589871",         /* 65 */
		"6470287906463336878241474855987746904297564226439499503918586590778209",    /* 70 */
		"281396163585532137380297959872159569353696836686080935550459706878100362721",           /* 74 */
		"33727095233132290409342213138708322681737322487170896778164145844669592994743377",      /* 80 */
		"1877138824359859508015524119652506869600959721781289179190693027302028679377371001561",  /* 85 */
		"427351849650748515507228344120452096326780093349980867041485502247153375067354165128307841",
		"48404068520546498995797968938385187958997290617596242601254422967869040251141325866025672337021",
		"1802716097522165018257858828415111497060066282677325501816640492782221110851604465066510547671104729"
	};
	const int digit_sizes[] = { 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100 };
	int num_bench = 12;

	// ------------------------------------------------------------------ //
	//  Storage for per-input timings captured before reset_factobj()      //
	// ------------------------------------------------------------------ //
	double qs_time[12];
	double total_time[12];

	// ------------------------------------------------------------------ //
	//  Build the System cell content                                       //
	//                                                                      //
	//  Each piece of info on its own line via <br>:                        //
	//    Line 1 (bold): CPU model        (comp_info->idstr)               //
	//    Line 2:        OS               (compile-time macro)             //
	//    Line 3:        Compiler+version (mirrors print_splash logic)     //
	//    Line 4:        SIMD feature set (highest enabled; BMI2 shown     //
	//                   only for AVX2-and-below builds)                  //
	//    Line 5:        Thread count     (fobj->THREADS)                  //
	// ------------------------------------------------------------------ //
	char compiler_str[128];
	char simd_str[128];
	char os_str[32];
	char system_cell[512];

#if defined(_WIN32) || defined(_WIN64)
	strcpy(os_str, "Windows");
#elif defined(__linux__)
	// Detect WSL and, if found, include the distro name from WSL_DISTRO_NAME.
	// Fall back to /proc/version string search if the env var isn't set.
	{
		int is_wsl = 0;
		const char* wsl_distro = getenv("WSL_DISTRO_NAME");

		if (wsl_distro != NULL && wsl_distro[0] != '\0')
		{
			// WSL_DISTRO_NAME is set — use it directly for the distro label
			snprintf(os_str, sizeof(os_str), "Linux (WSL/%s)", wsl_distro);
			is_wsl = 1;
		}

		if (!is_wsl)
		{
			// Fall back: inspect /proc/version for "Microsoft" or "WSL"
			FILE* pv = fopen("/proc/version", "r");
			if (pv != NULL)
			{
				char pvbuf[256];
				if (fgets(pvbuf, sizeof(pvbuf), pv) != NULL)
				{
					if (strstr(pvbuf, "Microsoft") != NULL || strstr(pvbuf, "WSL") != NULL)
					{
						strcpy(os_str, "Linux (WSL)");
						is_wsl = 1;
					}
				}
				fclose(pv);
			}
		}

		if (!is_wsl)
			strcpy(os_str, "Linux");
	}
#elif defined(__APPLE__)
	strcpy(os_str, "macOS");
#else
	strcpy(os_str, "Unknown OS");
#endif

	// Compiler identification — mirrors the logic in print_splash()
#if defined(_MSC_VER)
#  if defined(__INTEL_COMPILER)
	snprintf(compiler_str, sizeof(compiler_str), "MSVC %d + ICC %d", _MSC_VER, __INTEL_COMPILER);
#  elif defined(__INTEL_LLVM_COMPILER)
	snprintf(compiler_str, sizeof(compiler_str), "MSVC %d + ICX %s", _MSC_VER, __clang_version__);
#  elif defined(__clang_version__)
	snprintf(compiler_str, sizeof(compiler_str), "MSVC %d + Clang %s", _MSC_VER, __clang_version__);
#  else
	snprintf(compiler_str, sizeof(compiler_str), "MSVC %d", _MSC_VER);
#  endif
#elif defined(__INTEL_COMPILER)
	snprintf(compiler_str, sizeof(compiler_str), "ICC %d", __INTEL_COMPILER);
#elif defined(__INTEL_LLVM_COMPILER)
	snprintf(compiler_str, sizeof(compiler_str), "ICX %s", __clang_version__);
#elif defined(__clang_version__)
	snprintf(compiler_str, sizeof(compiler_str), "Clang %s", __clang_version__);
#elif defined(__GNUC__)
	snprintf(compiler_str, sizeof(compiler_str), "GCC %d.%d", __GNUC__, __GNUC_MINOR__);
#else
	snprintf(compiler_str, sizeof(compiler_str), "Unknown compiler");
#endif

	// SIMD feature set — highest enabled level first.
	// BMI2 is only appended for AVX2-and-below builds because AVX512 implies BMI2.
	// Guards use the same USE_* / IFMA macros as the rest of the build so
	// the string reflects what was actually compiled in.
	simd_str[0] = '\0';
	int simd_is_avx512 = 0;
#ifdef IFMA
	if (comp_info->AVX512IFMA) { strcat(simd_str, "AVX512IFMA"); simd_is_avx512 = 1; }
	else
#endif
#ifdef USE_AVX512BW
		if (comp_info->AVX512BW) { strcat(simd_str, "AVX512BW");   simd_is_avx512 = 1; }
		else
#endif
#ifdef USE_AVX512F
			if (comp_info->AVX512F) { strcat(simd_str, "AVX512F");    simd_is_avx512 = 1; }
			else
#endif
#ifdef USE_AVX2
				if (comp_info->AVX2) { strcat(simd_str, "AVX2"); }
				else
#endif
#ifdef USE_AVX
					if (comp_info->AVX) { strcat(simd_str, "AVX"); }
					else
#endif
#ifdef USE_SSE41
						if (comp_info->bSSE41Extensions) { strcat(simd_str, "SSE41"); }
						else
#endif
						{
							strcat(simd_str, "SSE2");
						}

#ifdef USE_BMI2
	// Only show BMI2 explicitly when it adds information (AVX2 and below)
	if (!simd_is_avx512 && comp_info->BMI2) { strcat(simd_str, " BMI2"); }
#endif

	if (simd_str[0] == '\0')
		strcpy(simd_str, "generic");

	// Compose System cell: bold CPU name, then one line per attribute
	snprintf(system_cell, sizeof(system_cell),
		"**%s** <br>%s<br>%s<br>%s<br>%d threads",
		comp_info->idstr,
		os_str,
		compiler_str,
		simd_str,
		fobj->THREADS);

	// ------------------------------------------------------------------ //
	//  Run the benchmarks                                                  //
	// ------------------------------------------------------------------ //
	fact_obj_t f;
	init_factobj(&f);
	copy_factobj(&f, fobj, 1);
	strcpy(f.flogname, "bench.log");

	{
		FILE* log = fopen(f.flogname, "a");
		if (log == NULL)
		{
			printf("fopen error: %s\n", strerror(errno));
			printf("couldn't open %s for writing\n", f.flogname);
			exit(1);
		}
		fprintf(log, "commencing siqsbench on %s\n", comp_info->idstr);
		fclose(log);
	}

	for (int i = 0; i < num_bench; i++)
	{
		if (digit_sizes[i] > digit_limit)
		{
			break;
		}
		mpz_set_str(f.qs_obj.gmp_n, list[i], 10);
		mpz_set(f.input_N, f.qs_obj.gmp_n);
		mpz_set(f.N, f.qs_obj.gmp_n);

		SIQS(&f);

		// Capture timings before reset_factobj() clears them
		qs_time[i] = f.qs_obj.qs_time;
		total_time[i] = f.qs_obj.total_time;

		reset_factobj(&f);
	}

	strcpy(f.flogname, "bench.log");

	// ------------------------------------------------------------------ //
	//  Emit Markdown                                                       //
	// ------------------------------------------------------------------ //
	//
	// Rendered table layout:
	//
	//   | System                  | Version | 44d | 50d | ... | 85d |
	//   |:------------------------|:-------:|----:|----:|...|----:|
	//   | **CPU**<br>OS<br>...    | 2.11    | ... |
	//
	// The header is written only when the file has zero length (new file).
	// ------------------------------------------------------------------ //
	{
		const char* mdname = "siqs_bench.md";
		int write_header = 0;

		// Check whether the file is empty / does not yet exist
		FILE* probe = fopen(mdname, "r");
		if (probe == NULL)
		{
			write_header = 1;   // file absent
		}
		else
		{
			fseek(probe, 0, SEEK_END);
			if (ftell(probe) == 0)
				write_header = 1;   // file present but empty
			fclose(probe);
		}

		FILE* md = fopen(mdname, "a");
		if (md == NULL)
		{
			printf("fopen error: %s\n", strerror(errno));
			printf("couldn't open %s for writing\n", mdname);
			free_factobj(&f);
			return;
		}

		if (write_header)
		{
			fprintf(md,
				"## SIQS Benchmark Results\n\n"
				"Times are wall-clock seconds: **qs\\_time** (sieve phase only),\n"
				"with **total\\_time** (including pre/post-processing) in parentheses.\n"
				"Digit counts refer to the number of decimal digits in the input.\n\n"
			);

			// Header row: System, Version, then one column per digit size
			fprintf(md, "| System | Version |");
			for (int i = 0; i < num_bench; i++)
				fprintf(md, " %dd |", digit_sizes[i]);
			fprintf(md, "\n");

			// Separator: System left-aligned, Version centered, timings right-aligned
			fprintf(md, "|:-------|:-------:|");
			for (int i = 0; i < num_bench; i++)
				fprintf(md, "------:|");
			fprintf(md, "\n");
		}

		// Data row
		fprintf(md, "| %s | %s |", system_cell, fobj->yafu_version);
		for (int i = 0; i < num_bench; i++)
		{
			if (digit_sizes[i] > digit_limit)
			{
				fprintf(md, " |");
			}
			else
			{
				fprintf(md, " %.2f<br>(%.2f) |", qs_time[i], total_time[i]);
			}
		}
			
		fprintf(md, "\n");

		fclose(md);
		printf("Benchmark results appended to %s\n", mdname);
	}

	free_factobj(&f);
	return;
}

void siqsbench_old(fact_obj_t* fobj)
{
	//run through the list of benchmark siqs factorizations

	char list[10][200];
	enum cpu_type cpu;
	FILE* log;
	int i;

	strcpy(list[0], "405461849292216354219321922871108605045931309");
	strcpy(list[1], "29660734457033883936073030405220515257819037444591");
	strcpy(list[2], "1592739380299790264815370870751381488173344970710983383");
	strcpy(list[3], "349594255864176572614071853194924838158088864370890996447417");
	strcpy(list[4], "34053408309992030649212497354061832056920539397279047809781589871");
	strcpy(list[5], "6470287906463336878241474855987746904297564226439499503918586590778209");
	strcpy(list[6], "281396163585532137380297959872159569353696836686080935550459706878100362721");
	strcpy(list[7], "33727095233132290409342213138708322681737322487170896778164145844669592994743377");
	strcpy(list[8], "1877138824359859508015524119652506869600959721781289179190693027302028679377371001561");

	int num_bench = 9;

	fact_obj_t f;
	init_factobj(&f);
	copy_factobj(&f, fobj, 1);

	strcpy(f.flogname, "bench.log");
	log = fopen(f.flogname, "a");
	if (log == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("couldn't open %s for writing\n", f.flogname);
		exit(1);
	}

	fprintf(log, "commencing siqsbench on %s\n", f.CPU_ID_STR);
	fclose(log);

	for (i = 0; i < num_bench; i++)
	{
		mpz_set_str(f.qs_obj.gmp_n, list[i], 10);
		mpz_set(f.input_N, f.qs_obj.gmp_n);
		mpz_set(f.N, f.qs_obj.gmp_n);
		SIQS(&f);
		reset_factobj(&f);
	}

	strcpy(f.flogname, "bench.log");

	free_factobj(&f);

	return;
}

void get_dummy_params(int bits, uint32_t* B, uint32_t* M, uint32_t* NB)
{
	int i;
	double scale;

	//parameter table
	//bits, fb primes, lp mulitplier, 64k blocks
	int param_table[22][4] = {
		{140,	600,	40,	1},
		{149,	875,	40,	1},
		{165,	1228,	50,	1},
		{181,	2247,	50,	1},
		{198,	3485,	60,	2},
		{215,	6357,	60,	2},
		{232,	12132,	70,	3},
		{248,	26379,	80,	4},
		{265,	42871,	90,	6},
		{281,	55137,	100,	8},
		{298,	65244,	120,	10},
		{310,	78247,	120,	12},
		{320,	90678,	140,	14},
		{330,	105000, 150,    18},
		{340,	125000, 150,    21},
		{350,	155000, 150,    25},
		{360,	195000, 150,    29},
		{370,	250000, 150,    34},
		{380,	310000, 150,    40},
		{390,	380000, 150,    47},
		{400,	460000, 150,    55},
		{410,	550000, 150,    64},
	};

	//linear interpolation according to bit size to determine
	//factor base bound.  use the closest parameter for lp multiplier
	//and number of blocks.

	if (bits <= param_table[0][0])
	{
		scale = (double)bits / (double)param_table[0][0];
		*B = (uint32_t)(scale * (double)(param_table[0][1]));
		*M = 40;
		*NB = 1;
	}
	else
	{
		for (i = 0; i < 22; i++)
		{
			if (bits > param_table[i][0] && bits <= param_table[i + 1][0])
			{
				scale = (double)(param_table[i + 1][0] - bits) /
					(double)(param_table[i + 1][0] - param_table[i][0]);
				*B = param_table[i + 1][1] -
					(uint32_t)(scale * (double)(param_table[i + 1][1] - param_table[i][1]));

				*M = (uint32_t)((param_table[i + 1][2] + param_table[i][2]) / 2.0 + 0.5);
				*NB = (uint32_t)((param_table[i + 1][3] + param_table[i][3]) / 2.0 + 0.5);
			}
		}
	}

	if (*B == 0)
	{
		//off the end of the table, extrapolate based on the slope of 
		//the last two

		scale = (double)(param_table[21][1] - param_table[20][1]) /
			(double)(param_table[21][0] - param_table[20][0]);
		*B = (uint32_t)(((double)bits - param_table[21][0]) *
			scale + param_table[21][1]);
		*M = param_table[21][2];	//reuse last one

		scale = (double)(param_table[21][3] - param_table[20][3]) /
			(double)(param_table[21][0] - param_table[20][0]);
		*NB = (uint32_t)(((double)bits - param_table[21][0]) *
			scale + param_table[21][3]);
	}

	return;
}

int check_poly_at_loc(int loc, int parity, int prime, int polynum, int index,
	int* polyv, int* polysign, dynamic_conf_t* dconf, static_conf_t* sconf)
{
	int result = 0;

	// staring from the first b-poly, compute the indicated b-poly index, verify the polynum,
	// and verify that the provided x-coordinate of the poly divides N.

	uint32_t Bnum = 1;
	mpz_ptr n = sconf->n;
	mpz_t polyb, polyb1, polyb2, polyc, Q;
	siqs_poly* poly = dconf->curr_poly;
	int i, s = poly->s;

	mpz_init(polyb);
	mpz_init(polyb1);
	mpz_init(polyb2);
	mpz_init(polyc);
	mpz_init(Q);

	//initialize b
	mpz_set_ui(polyb, 0);

	// build up first b
	for (i = 0; i < s; i++)
	{
		mpz_add(polyb, polyb, dconf->Bl[i]);
	}
	mpz_tdiv_q_2exp(polyb, polyb, 1);

	mpz_set(polyb1, dconf->Bl[0]);

	for (i = 1; i < s / 2; i++) {
		//gmp_printf("%Zd\n", dconf->Bl[ii]);
		mpz_add(polyb1, polyb1, dconf->Bl[i]);
	}
	mpz_tdiv_q_2exp(polyb1, polyb1, 1);

	mpz_set(polyb2, dconf->Bl[i++]);

	for (; i < s; i++) {
		//gmp_printf("%Zd\n", dconf->Bl[ii]);
		mpz_add(polyb2, polyb2, dconf->Bl[i]);
	}
	mpz_tdiv_q_2exp(polyb2, polyb2, 1);


	// need to do this whether or not we compute 'c'.
	if (sconf->knmod8 == 1)
	{
		if (mpz_tstbit(polyb, 0) == 0)
		{
			mpz_add(polyb, polyb, poly->mpz_poly_a);
		}
	}

	int p = 0;
	for (; Bnum < index; Bnum++)
	{
		//compute the next b
		if (poly->gray[Bnum] < 0)
			mpz_sub(polyb, polyb, dconf->Bl[poly->nu[Bnum] - 1]);
		else
			mpz_add(polyb, polyb, dconf->Bl[poly->nu[Bnum] - 1]);

		// need to do this whether or not we compute 'c'.
		if (sconf->knmod8 == 1)
		{
			if (mpz_tstbit(polyb, 0) == 0)
			{
				mpz_add(polyb, polyb, poly->mpz_poly_a);
			}
		}

		int v;
		int tmp;
		int sign;
		int j;

		// next poly enumeration
		v = 1;
		j = Bnum;

		while ((j & 1) == 0)
			v++, j >>= 1;
		tmp = Bnum + (1 << v) - 1;
		tmp = (tmp >> v);
		if (tmp & 1)
			sign = -1;
		else
			sign = 1;

		// next polynum
		p ^= (1 << (v - 1));
	}

	if (p != polynum)
	{
		printf("polynum %04x at index %d does not match input polynum %04x", p, index, polynum);
	}

	// compute C
	if (sconf->knmod8 == 1)
	{
		mpz_mul(polyc, polyb, polyb);
		mpz_sub(polyc, polyc, n);
		mpz_tdiv_q(polyc, polyc, poly->mpz_poly_a);

		if (mpz_tdiv_ui(polyc, 4) != 0)
		{
			printf("c not divisible by 4 in Q2(x) variation!\n");
		}
		mpz_tdiv_q_ui(polyc, polyc, 4);
	}
	else
	{
		//now that we have b, compute c = (b*b - n)/a
		mpz_mul(polyc, polyb, polyb);
		mpz_sub(polyc, polyc, n);
		mpz_tdiv_q(polyc, polyc, poly->mpz_poly_a);
	}

	// now we have 'b' and 'c' at the provided poly index; compute Q
	if (sconf->knmod8 == 1)
	{
		// this one is close enough, compute 
		// Q(x) = (2ax + b)^2 - N, where x is the sieve index
		// Q(x)/4a = (ax + b)x + c;	
		mpz_mul_ui(dconf->gmptmp1, dconf->curr_poly->mpz_poly_a, loc);

		if (parity)
			mpz_sub(Q, dconf->gmptmp1, polyb);
		else
			mpz_add(Q, dconf->gmptmp1, polyb);

		mpz_mul_ui(Q, Q, loc);
		mpz_add(Q, Q, polyc);

		if (mpz_sgn(Q) < 0)
		{
			mpz_neg(Q, Q);
		}
	}
	else
	{
		// this one is close enough, compute 
		// Q(x) = (ax + b)^2 - N, where x is the sieve index
		// Q(x)/a = (ax + 2b)x + c;	
		mpz_mul_2exp(dconf->gmptmp2, polyb, 1);
		mpz_mul_ui(dconf->gmptmp1, dconf->curr_poly->mpz_poly_a, loc);

		if (parity)
			mpz_sub(Q, dconf->gmptmp1, dconf->gmptmp2);
		else
			mpz_add(Q, dconf->gmptmp1, dconf->gmptmp2);

		mpz_mul_ui(Q, Q, loc);
		mpz_add(Q, Q, polyc);

		if (mpz_sgn(Q) < 0)
		{
			mpz_neg(Q, Q);
		}
	}

	// finally make sure prime divides this Q.
	if (mpz_tdiv_ui(Q, prime) == 0)
	{
		//gmp_printf("pass: prime %d divides Q = %Zd at x = %d\n", prime, Q, loc);
		result = 1;
	}
	else
	{
		for (i = sconf->factor_base->x2_large_B; i < sconf->factor_base->B; i++)
		{
			if (prime == sconf->factor_base->list->prime[i])
				break;
		}

		int pidx = i;

		int root1 = sconf->modsqrt_array[pidx];
		int root2 = prime - root1;
		int amodp = mpz_tdiv_ui(dconf->curr_poly->mpz_poly_a, prime);
		int bmodp = mpz_tdiv_ui(polyb, prime);
		int inv;

		//find a^-1 mod p = inv(a mod p) mod p
		if (sconf->knmod8 == 1)
		{
			inv = modinv_1(2 * amodp, prime);
		}
		else
		{
			inv = modinv_1(amodp, prime);
		}

		root1 = (int)root1 - bmodp;
		root2 = (int)root2 - bmodp;
		if (root1 < 0) root1 += prime;
		if (root2 < 0) root2 += prime;

		root1 = (uint32_t)((uint64_t)inv * (uint64_t)root1 % (uint64_t)prime);
		root2 = (uint32_t)((uint64_t)inv * (uint64_t)root2 % (uint64_t)prime);

		gmp_printf("fail: prime %d does not divide Q = %Zd at x = %d\n", prime, Q, loc);
		gmp_printf("knmod8 = %d\nA = %Zd\nB = %Zd\nC = %Zd\nN = %Zd\n",
			sconf->knmod8, dconf->curr_poly->mpz_poly_a, polyb, polyc, n);
		printf("t = %d,%d; prime = %d\namodp = %d\nbmodp = %d\nainvp = %d\n",
			sconf->modsqrt_array[pidx], prime - sconf->modsqrt_array[pidx],
			sconf->factor_base->list->prime[pidx], amodp, bmodp, inv);
		printf("first roots = %d,%d (%d,%d)\n", dconf->update_data.firstroots1[pidx],
			dconf->update_data.firstroots2[pidx], root1, root2);
		gmp_printf("b1 = %Zd\nb2 = %Zd\n", polyb1, polyb2);

		root1 = sconf->modsqrt_array[pidx];
		root2 = prime - root1;

		int bmodp1 = mpz_tdiv_ui(polyb1, prime);
		int bmodp2 = mpz_tdiv_ui(polyb2, prime);

		int r1 = (int)root1 - bmodp1;
		int r2 = (int)root2 - bmodp1;
		if (r1 < 0) r1 += prime;
		if (r2 < 0) r2 += prime;

		r1 = (int)((uint64_t)inv * (uint64_t)r1 % (uint64_t)prime);
		r2 = (int)((uint64_t)inv * (uint64_t)r2 % (uint64_t)prime);

		printf("left  side r1 = %d, r2 = %d\n", r1, r2);

		r1 = prime - bmodp2;
		r1 = (int)((uint64_t)inv * (uint64_t)r1 % (uint64_t)prime);

		printf("right side r1 = %d\n", r1);

		result = 0;
	}

	mpz_clear(Q);
	mpz_clear(polyb);
	mpz_clear(polyc);

	return result;
}
