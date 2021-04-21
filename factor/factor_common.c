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
       				   --bbuhrow@gmail.com 3/26/10
----------------------------------------------------------------------*/

#include "factor.h"
#include "mpz_aprcl.h"
#include "soe.h"

void init_factobj(fact_obj_t* fobj)
{
    uint32_t seed1, seed2;
    int i;

    // initialize general stuff in fobj
    get_random_seeds(&seed1, &seed2);
    fobj->seed1 = seed1;
    fobj->seed2 = seed2;
    fobj->lcg_state = ((uint64_t)seed2 << 32) | (uint64_t)seed1;
    fobj->flags = 0;
    fobj->num_threads = 1;
    strcpy(fobj->flogname, "factor.log");
    fobj->do_logging = 1;   // not used...
    fobj->LOGFLAG = 1;
    fobj->NUM_WITNESSES = 1;

    // get space for everything
    alloc_factobj(fobj);

    // initialize stuff for rho	
    fobj->rho_obj.iterations = 1000;
    fobj->rho_obj.curr_poly = 0;

    // initialize stuff for pm1	
    fobj->pm1_obj.B1 = 100000;
    fobj->pm1_obj.B2 = 10000000;
    fobj->pm1_obj.stg2_is_default = 1;
    fobj->pm1_obj.pm1_exponent = 0;
    fobj->pm1_obj.pm1_multiplier = 0;
    fobj->pm1_obj.pm1_tune_freq = 0;
    fobj->pm1_obj.vecnum = 0;

    // initialize stuff for pp1	
    fobj->pp1_obj.B1 = 20000;
    fobj->pp1_obj.B2 = 1000000;
    fobj->pp1_obj.stg2_is_default = 1;
    fobj->pp1_obj.pp1_exponent = 0;
    fobj->pp1_obj.pp1_multiplier = 0;
    fobj->pp1_obj.pp1_tune_freq = 0;
    fobj->pp1_obj.vecnum = 0;

    // initialize stuff for ecm	
    fobj->ecm_obj.B1 = 11000;
    fobj->ecm_obj.B2 = 1100000;
    fobj->ecm_obj.stg2_is_default = 1;
    fobj->ecm_obj.sigma = 0;
    fobj->ecm_obj.num_curves = 90;
    fobj->ecm_obj.curves_run = 0;
    fobj->ecm_obj.ecm_exponent = 0;
    fobj->ecm_obj.ecm_multiplier = 0;
    fobj->ecm_obj.ecm_tune_freq = 0;
    fobj->ecm_obj.bail_on_factor = 1;
    fobj->ecm_obj.save_b1 = 0;

    // unlike ggnfs, ecm does not *require* external binaries.  
    // an empty string indicates the use of the built-in GMP-ECM hooks, while
    // a non-empty string (filled in by the user) will indicate the use of
    // an external binary
    strcpy(fobj->ecm_obj.ecm_path, "");
    fobj->ecm_obj.use_external = 0;
#ifdef USE_AVX512F
    fobj->ecm_obj.prefer_gmpecm = 0;
    fobj->ecm_obj.ecm_ext_xover = 40000000;
#else
    fobj->ecm_obj.prefer_gmpecm = 1;
    fobj->ecm_obj.ecm_ext_xover = 48000;
#endif

    // initialize stuff for squfof
    fobj->squfof_obj.num_factors = 0;

    // initialize stuff for qs	
    fobj->qs_obj.gbl_override_B_flag = 0;
    fobj->qs_obj.gbl_override_B = 0;
    fobj->qs_obj.gbl_override_small_cutoff_flag = 0;
    fobj->qs_obj.gbl_override_small_cutoff = 0;
    fobj->qs_obj.gbl_override_blocks_flag = 0;
    fobj->qs_obj.gbl_override_blocks = 0;
    fobj->qs_obj.gbl_override_lpmult_flag = 0;
    fobj->qs_obj.gbl_override_lpmult = 0;
    fobj->qs_obj.gbl_override_rel_flag = 0;
    fobj->qs_obj.gbl_override_rel = 0;
    fobj->qs_obj.gbl_override_tf_flag = 0;
    fobj->qs_obj.gbl_override_tf = 0;
    fobj->qs_obj.gbl_override_time_flag = 0;
    fobj->qs_obj.gbl_override_time = 0;
    fobj->qs_obj.gbl_override_mfbd = 0.;
    fobj->qs_obj.gbl_override_mfbt = 0.;
    fobj->qs_obj.gbl_override_lpb = 0;
    fobj->qs_obj.gbl_override_bdiv_flag = 0;
    fobj->qs_obj.gbl_override_bdiv = 3;
    fobj->qs_obj.gbl_override_3lp_bat = 0;
    fobj->qs_obj.gbl_btarget = 1000000;
    fobj->qs_obj.flags = 0;
    fobj->qs_obj.gbl_force_DLP = 0;
    fobj->qs_obj.gbl_force_TLP = 0;
    fobj->qs_obj.qs_exponent = 0;
    fobj->qs_obj.qs_multiplier = 0;
    fobj->qs_obj.qs_tune_freq = 0;
    fobj->qs_obj.no_small_cutoff_opt = 0;
    strcpy(fobj->qs_obj.siqs_savefile, "siqs.dat");
    init_lehman();

    // initialize stuff for trial division	
    fobj->div_obj.print = 0;
    fobj->div_obj.limit = 10000;
    fobj->div_obj.fmtlimit = 1000000;

    //initialize stuff for nfs
    fobj->nfs_obj.snfs = 0;
    fobj->nfs_obj.gnfs = 0;
    fobj->nfs_obj.gnfs_exponent = 0;
    fobj->nfs_obj.gnfs_multiplier = 0;
    fobj->nfs_obj.gnfs_tune_freq = 0;
    fobj->nfs_obj.min_digits = 85;
    fobj->nfs_obj.filter_min_rels_nudge = 1.05;	// raise min_rels bounds by 5%
                                                    // on unsuccessful filtering
    fobj->nfs_obj.siever = 0;					//default, use automatic selection
    fobj->nfs_obj.startq = 0;					//default, not used
    fobj->nfs_obj.rangeq = 0;					//default, not used
    fobj->nfs_obj.polystart = 0;					//default, not used
    fobj->nfs_obj.polyrange = 0;					//default, not used
    strcpy(fobj->nfs_obj.outputfile, "nfs.dat");			//default
    strcpy(fobj->nfs_obj.logfile, "nfs.log");			//default
    strcpy(fobj->nfs_obj.fbfile, "nfs.fb");				//default
    fobj->nfs_obj.sq_side = 0;					//default = algebraic
    fobj->nfs_obj.timeout = 0;					//default, not used
    strcpy(fobj->nfs_obj.job_infile, "nfs.job");			//default
    fobj->nfs_obj.poly_option = 0;					//default = fast search
                                        //1 = wide
                                    //2 = deep
    fobj->nfs_obj.restart_flag = 0;					//default = not a restart
    fobj->nfs_obj.nfs_phases = NFS_DEFAULT_PHASES;
    fobj->nfs_obj.snfs_testsieve_threshold = 160;
    *fobj->nfs_obj.filearg = '\0';

    fobj->nfs_obj.polybatch = 250;						//default	
#if defined(_WIN64)
    strcpy(fobj->nfs_obj.ggnfs_dir, ".\\");
#elif defined(WIN32)
    strcpy(fobj->nfs_obj.ggnfs_dir, ".\\");
#else
    strcpy(fobj->nfs_obj.ggnfs_dir, "./");
#endif

    //initialize autofactor object
    //whether we want to output certain info to their own files...
    fobj->autofact_obj.want_output_primes = 0;
    fobj->autofact_obj.want_output_factors = 0;
    fobj->autofact_obj.want_output_unfactored = 0;
    fobj->autofact_obj.want_output_expressions = 1;
    fobj->autofact_obj.qs_gnfs_xover = 95;
    fobj->autofact_obj.qs_snfs_xover = 75;
    // use xover even when timing info is available
    fobj->autofact_obj.prefer_xover = 0;
    fobj->autofact_obj.want_only_1_factor = 0;
    fobj->autofact_obj.no_ecm = 0;
    fobj->autofact_obj.target_pretest_ratio = 4.0 / 13.0;
    fobj->autofact_obj.initial_work = 0.0;
    fobj->autofact_obj.has_snfs_form = -1;		// not checked yet

    //pretesting plan used by factor()
    fobj->autofact_obj.yafu_pretest_plan = PRETEST_NORMAL;
    strcpy(fobj->autofact_obj.plan_str, "normal");
    fobj->autofact_obj.only_pretest = 0;
    fobj->autofact_obj.autofact_active = 0;

    // if a number is <= aprcl_prove_cutoff, we will prove it prime or composite
    fobj->factors->aprcl_prove_cutoff = 500;
    // if a number is >= aprcl_display_cutoff, we will show the APRCL progress
    fobj->factors->aprcl_display_cutoff = 200;

    return;
}

void free_factobj(fact_obj_t *fobj)
{
    int i;

	// free stuff in rho
	free(fobj->rho_obj.polynomials);
	mpz_clear(fobj->rho_obj.gmp_n);
	mpz_clear(fobj->rho_obj.gmp_f);

	// free stuff in pm1
	mpz_clear(fobj->pm1_obj.gmp_n);
	mpz_clear(fobj->pm1_obj.gmp_f);
    for (i = 0; i < fobj->pm1_obj.vecnum; i++)
    {
        mpz_clear(fobj->pm1_obj.vecn[i]);
    }
    fobj->pm1_obj.vecnum = 0;

	// free stuff in pp1
	mpz_clear(fobj->pp1_obj.gmp_n);
	mpz_clear(fobj->pp1_obj.gmp_f);
    for (i = 0; i < fobj->pp1_obj.vecnum; i++)
    {
        mpz_clear(fobj->pp1_obj.vecn[i]);
    }
    fobj->pp1_obj.vecnum = 0;

	// then free other stuff in ecm
	mpz_clear(fobj->ecm_obj.gmp_n);
	mpz_clear(fobj->ecm_obj.gmp_f);
    free(fobj->ecm_obj.lcg_state);

	// then free other stuff in squfof
	mpz_clear(fobj->squfof_obj.gmp_n);
	mpz_clear(fobj->squfof_obj.gmp_f);

	// then free other stuff in qs
	mpz_clear(fobj->qs_obj.gmp_n);

	// then free other stuff in div
	mpz_clear(fobj->div_obj.gmp_n);
	mpz_clear(fobj->div_obj.gmp_f);

	// then free other stuff in nfs
	mpz_clear(fobj->nfs_obj.gmp_n);
	mpz_clear(fobj->nfs_obj.snfs_cofactor);

	//free general fobj stuff
	mpz_clear(fobj->N);

	clear_factor_list(fobj->factors);
	free(fobj->factors);
    free(fobj->primes);
    fobj->min_p = 0;
    fobj->max_p = 0;
    fobj->num_p = 0;

	return;
}

void alloc_factobj(fact_obj_t *fobj)
{
	int i;
	
	mpz_init(fobj->N);

	fobj->rho_obj.num_poly = 3;
	fobj->rho_obj.polynomials = (uint32_t *)xmalloc(fobj->rho_obj.num_poly * sizeof(uint32_t));
	fobj->rho_obj.polynomials[0] = 3;
	fobj->rho_obj.polynomials[1] = 2;
	fobj->rho_obj.polynomials[2] = 1;
	mpz_init(fobj->rho_obj.gmp_n);
	mpz_init(fobj->rho_obj.gmp_f);

	mpz_init(fobj->pm1_obj.gmp_n);
	mpz_init(fobj->pm1_obj.gmp_f);
	
	mpz_init(fobj->pp1_obj.gmp_n);
	mpz_init(fobj->pp1_obj.gmp_f);

	mpz_init(fobj->ecm_obj.gmp_n);
	mpz_init(fobj->ecm_obj.gmp_f);
    fobj->ecm_obj.lcg_state = (uint64_t*)xmalloc(fobj->num_threads * sizeof(uint64_t));

	mpz_init(fobj->squfof_obj.gmp_n);
	mpz_init(fobj->squfof_obj.gmp_f);

	mpz_init(fobj->qs_obj.gmp_n);

	mpz_init(fobj->div_obj.gmp_n);
	mpz_init(fobj->div_obj.gmp_f);

	mpz_init(fobj->nfs_obj.gmp_n);
	mpz_init(fobj->nfs_obj.snfs_cofactor);

    fobj->factors = (yfactor_list_t*)xmalloc(1 * sizeof(yfactor_list_t));
	fobj->factors->alloc_factors = 256;
	fobj->factors->factors = (yfactor_t *)xmalloc(256 * sizeof(yfactor_t));
	for (i=0; i < fobj->factors->alloc_factors; i++)
	{
		fobj->factors->factors[i].type = UNKNOWN;
		fobj->factors->factors[i].count = 0;
	}

	fobj->ecm_obj.num_factors = 0;	
	fobj->qs_obj.num_factors = 0;	
	fobj->div_obj.num_factors = 0;	
	fobj->nfs_obj.num_factors = 0;	
	fobj->factors->num_factors = 0;

    soe_staticdata_t* sdata;
    // put a list of primes in the fobj; many algorithms use it
    sdata = soe_init(0, 1, 32768);
    fobj->primes = soe_wrapper(sdata, 0, 10000000, 0, &fobj->num_p, 0, 0);
    fobj->min_p = 2;
    fobj->max_p = fobj->primes[fobj->num_p - 1];
    soe_finalize(sdata);

	return;
}

void reset_factobj(fact_obj_t *fobj)
{
	// keep all of the settings in fobj, but do an init/free cycle on all
	// allocated structures
	free_factobj(fobj);
	alloc_factobj(fobj);

	return;
}

int add_to_factor_list(yfactor_list_t *flist, mpz_t n, int VFLAG, int NUM_WITNESSES)
{
	// stick the number n into the provided factor list.
    // return the index into which the factor was added.
	int i;
    int fid;
	int found = 0, v = 0;

	// look to see if this factor is already in the list
    for (i = 0; i < flist->num_factors && !found; i++)
    {
        if (mpz_cmp(n, flist->factors[i].factor) == 0)
        {
            found = 1;
            flist->factors[i].count++;
            return i;
        }
    }

    // otherwise, allocate another factor to the list and add it
    fid = flist->num_factors;
    if (flist->num_factors >= flist->alloc_factors)
    {
        flist->alloc_factors *= 2;
        flist->factors = (yfactor_t*)realloc(flist->factors,
            flist->alloc_factors * sizeof(yfactor_t));
    }

	mpz_init(flist->factors[fid].factor);
	mpz_set(flist->factors[fid].factor, n);
    flist->factors[fid].count = 1;
	if (gmp_base10(n) <= flist->aprcl_prove_cutoff) /* prove primality of numbers <= aprcl_prove_cutoff digits */
	{
		int ret = 0;

		if (VFLAG > 0)
			v = (gmp_base10(n) < flist->aprcl_display_cutoff) ? APRTCLE_VERBOSE0 : APRTCLE_VERBOSE1;
		else
			v = APRTCLE_VERBOSE0;

		if (v == APRTCLE_VERBOSE1)
			printf("\n");
		ret = mpz_aprtcle(n, v);
		if (v == APRTCLE_VERBOSE1)
			printf("\n");

		if (ret == APRTCLE_PRIME)
            flist->factors[fid].type = PRIME;
		else
		{
			if (mpz_bpsw_prp(n) != PRP_COMPOSITE)
			{
				// if BPSW doesn't think this composite number is actually composite, then ask the user
				// to report this fact to the YAFU sub-forum at mersenneforum.org
				printf(" *** ATTENTION: BPSW issue found.  Please report the following number to:\n");
				printf(" *** ATTENTION: http://www.mersenneforum.org/forumdisplay.php?f=96\n");
				gmp_printf(" *** ATTENTION: n = %Zd\n", n);
			}
            flist->factors[fid].type = COMPOSITE;
		}
	}
	else if (is_mpz_prp(n, NUM_WITNESSES))
	{
		if (mpz_cmp_ui(n, 100000000) < 0)
            flist->factors[fid].type = PRIME;
		else
            flist->factors[fid].type = PRP;
	}
	else
        flist->factors[fid].type = COMPOSITE;

    flist->num_factors++;
	return fid;
}

void delete_from_factor_list(yfactor_list_t* flist, mpz_t n)
{
	// remove the number n from the global factor list
	int i;

	// find the factor
    for (i = 0; i < flist->num_factors; i++)
    {
        if (mpz_cmp(n, flist->factors[i].factor) == 0)
        {
            int j;
            // copy everything above this in the list back one position
            for (j = i; j < flist->num_factors - 1; j++)
            {
                mpz_set(flist->factors[j].factor, flist->factors[j + 1].factor);
                flist->factors[j].count = flist->factors[j + 1].count;
            }
            // remove the last one in the list
            flist->factors[j].count = 0;
            mpz_clear(flist->factors[j].factor);

            flist->num_factors--;
            break;
        }
    }

	return;
}

void clear_factor_list(yfactor_list_t * flist)
{
	int i;

	//clear this info
    for (i = 0; i < flist->num_factors; i++)
    {
        flist->factors[i].count = 0;
        mpz_clear(flist->factors[i].factor);
    }
    flist->num_factors = 0;

	return;
}

int resume_check_input_match(mpz_t file_n, mpz_t input_n, mpz_t common_fact, int VFLAG)
{
	// see if the value we
	// read from the file is really the same as the input number.
	// if the gcd of the two is a divisor of the file, then
	// this is a resume.  if the gcd is non-trival,
	// set the input to the value found in the file.

	mpz_t r;
	int retval;

	mpz_init(r);
	if (VFLAG > 1)
	{
		gmp_printf("input from file = %Zd\n", file_n);
		gmp_printf("input to yafu = %Zd\n", input_n);
	}

	//mpz_gcd(common_fact, file_n, input_n);	
	if (mpz_cmp_ui(file_n, 0) == 0)
		return 0;

	mpz_tdiv_qr(common_fact, r, input_n, file_n);	

	if (mpz_cmp_ui(r, 0) == 0)
	{
		retval = 1;
		if (VFLAG > 1)
		{
			gmp_printf("input matches with multiple of %Zd\n", common_fact);
		}
	}
	else 
	{
		retval = 0;
		mpz_set_ui(common_fact, 1);
	}
		
	mpz_clear(r);
	return retval;
}

void print_factors(yfactor_list_t* flist, mpz_t N, int VFLAG, int NUM_WITNESSES)
{
	uint32_t i;
	int j, v;
	mpz_t tmp, tmp2;

	//always print factors unless complete silence is requested
	if (VFLAG >= 0)
	{
		printf("\n\n***factors found***\n\n");
		mpz_init(tmp);
		mpz_set_ui(tmp, 1);

		for (i=0; i< flist->num_factors; i++)
		{

			//if (fobj->fobj_factors[i].type == COMPOSITE)
			//{
			//	for (j=0;j<fobj->fobj_factors[i].count;j++)
			//	{
			//		mpz_mul(tmp, tmp, fobj->fobj_factors[i].factor);
			//		gmp_printf("C%d = %Zd\n", gmp_base10(fobj->fobj_factors[i].factor),
			//			fobj->fobj_factors[i].factor);
			//	}
			//}
			//else if (fobj->fobj_factors[i].type == PRP)
			//{
			//	for (j=0;j<fobj->fobj_factors[i].count;j++)
			//	{
			//		mpz_mul(tmp, tmp, fobj->fobj_factors[i].factor);
			//		gmp_printf("PRP%d = %Zd\n", gmp_base10(fobj->fobj_factors[i].factor),
			//			fobj->fobj_factors[i].factor);
			//	}
			//}
			if (flist->factors[i].type == PRIME)
			{
				// don't redo APR-CL calculations already performed by add_to_factor_list
				for (j=0;j< flist->factors[i].count;j++)
				{
					mpz_mul(tmp, tmp, flist->factors[i].factor);
					gmp_printf("P%d = %Zd\n", gmp_base10(flist->factors[i].factor),
                        flist->factors[i].factor);
				}
			}
			else
			{
				//type not set, determine it now
				/* prove primality of numbers <= aprcl_prove_cutoff digits */
				if (gmp_base10(flist->factors[i].factor) <=
                    flist->aprcl_prove_cutoff)
				{
					int ret = 0;
					if (VFLAG > 0)
						v = (gmp_base10(flist->factors[i].factor) <
                            flist->aprcl_display_cutoff) ? APRTCLE_VERBOSE0 : APRTCLE_VERBOSE1;
					else
						v = APRTCLE_VERBOSE0;

					if (v == APRTCLE_VERBOSE1)
						printf("\n");
					ret = mpz_aprtcle(flist->factors[i].factor, v);
					if (v == APRTCLE_VERBOSE1)
						printf("\n");

					if (ret == APRTCLE_PRIME)
					{
						for (j=0;j< flist->factors[i].count;j++)
						{
							mpz_mul(tmp, tmp, flist->factors[i].factor);
							gmp_printf("P%d = %Zd\n", 
                                gmp_base10(flist->factors[i].factor),
                                flist->factors[i].factor);
						}
					}
					else
					{
						if (mpz_bpsw_prp(flist->factors[i].factor) != PRP_COMPOSITE)
						{
							// if BPSW doesn't think this composite number is actually composite, then ask the user
							// to report this fact to the YAFU sub-forum at mersenneforum.org
							printf(" *** ATTENTION: BPSW issue found.  Please report the following number to:\n");
							printf(" *** ATTENTION: http://www.mersenneforum.org/forumdisplay.php?f=96\n");
							gmp_printf(" *** ATTENTION: n = %Zd\n", flist->factors[i].factor);
						}
						for (j=0;j< flist->factors[i].count;j++)
						{
							mpz_mul(tmp, tmp, flist->factors[i].factor);
							gmp_printf("C%d = %Zd\n", gmp_base10(flist->factors[i].factor),
                                flist->factors[i].factor);
						}
					}
				}
				else if (is_mpz_prp(flist->factors[i].factor, NUM_WITNESSES))
				{
					for (j=0;j< flist->factors[i].count;j++)
					{
						mpz_mul(tmp, tmp, flist->factors[i].factor);
						if (mpz_cmp_ui(flist->factors[i].factor, 100000000) < 0)
							gmp_printf("P%d = %Zd\n", gmp_base10(flist->factors[i].factor),
                                flist->factors[i].factor);
						else
							gmp_printf("PRP%d = %Zd\n", gmp_base10(flist->factors[i].factor),
                                flist->factors[i].factor);
					}
				}
				else
				{
					for (j=0;j< flist->factors[i].count;j++)
					{
						mpz_mul(tmp, tmp, flist->factors[i].factor);
						gmp_printf("C%d = %Zd\n", gmp_base10(flist->factors[i].factor),
                            flist->factors[i].factor);
					}
				}
			}
		}

		mpz_init(tmp2);
		mpz_set(tmp2, N);
		if (mpz_cmp_ui(tmp2, 1) > 0)
		{
			// non-trivial N remaining... compute and display the known co-factor
			mpz_tdiv_q(tmp2, tmp2, tmp);
			if (mpz_cmp_ui(tmp2, 1) > 0)
			{
				/* prove primality of numbers <= aprcl_prove_cutoff digits */
				if (gmp_base10(tmp2) <= flist->aprcl_prove_cutoff)
				{
					int ret = 0;
					if (VFLAG > 0)
						v = (gmp_base10(tmp2) < flist->aprcl_display_cutoff) ? APRTCLE_VERBOSE0 : APRTCLE_VERBOSE1;
					else
						v = APRTCLE_VERBOSE0;

					if (v == APRTCLE_VERBOSE1)
						printf("\n");
					ret = mpz_aprtcle(tmp2, v);
					if (v == APRTCLE_VERBOSE1)
						printf("\n");

					if (ret == APRTCLE_PRIME)
						gmp_printf("\n***co-factor***\nP%d = %Zd\n",
							gmp_base10(tmp2), tmp2);
					else
					{
						if (mpz_bpsw_prp(tmp2) != PRP_COMPOSITE)
						{
							// if BPSW doesn't think this composite number is actually composite, then ask the user
							// to report this fact to the YAFU sub-forum at mersenneforum.org
							printf(" *** ATTENTION: BPSW issue found.  Please report the following number to:\n");
							printf(" *** ATTENTION: http://www.mersenneforum.org/forumdisplay.php?f=96\n");
							gmp_printf(" *** ATTENTION: n = %Zd\n", tmp2);
						}
						gmp_printf("\n***co-factor***\nC%d = %Zd\n",
							gmp_base10(tmp2), tmp2);
					}
				}
				else if (is_mpz_prp(tmp2, NUM_WITNESSES))
				{
					gmp_printf("\n***co-factor***\nPRP%d = %Zd\n",
						gmp_base10(tmp2), tmp2);
				}
				else
				{
					gmp_printf("\n***co-factor***\nC%d = %Zd\n",
						gmp_base10(tmp2), tmp2);
				}
			}			
		}
		mpz_clear(tmp);
		mpz_clear(tmp2);
	}

	return;
}


