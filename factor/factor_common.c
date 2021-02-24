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


void init_factobj(fact_obj_t *fobj, options_t *options)
{
    uint32_t seed1, seed2;
    int i;

    fobj->lcg_state = options->rand_seed;
    fobj->flags = 0;
    fobj->num_threads = options->threads;
    strcpy(fobj->flogname, options->factorlog);
    fobj->do_logging = 1;   // not used...

	// get space for everything
	alloc_factobj(fobj);

    // keep a reference to the options structure that dictates
    // much of how factoring behaves
    fobj->options = options;

	// initialize global stuff in fobj
    if (options->rand_seed == 0)
    {
        // random seeds
        get_random_seeds(&seed1, &seed2);
        fobj->seed1 = (uint32_t)options->rand_seed = seed1;
        fobj->seed2 = seed2;
        options->rand_seed += (uint64_t)seed2 << 32;
    }
    else
    {
        fobj->seed1 = (uint32_t)options->rand_seed & 0xffffffff;
        fobj->seed2 = (uint32_t)(options->rand_seed >> 32);
    }

	// initialize stuff for rho	
	fobj->rho_obj.iterations = options->rhomax;		
	fobj->rho_obj.curr_poly = 0;	

	// initialize stuff for pm1	
	fobj->pm1_obj.B1 = options->B1pm1;
	fobj->pm1_obj.B2 = options->B2pm1;
	fobj->pm1_obj.stg2_is_default = (options->B2pm1 == 0);		
	fobj->pm1_obj.pm1_exponent = 0;
	fobj->pm1_obj.pm1_multiplier = 0;
	fobj->pm1_obj.pm1_tune_freq = 0;

	// initialize stuff for pp1	
	fobj->pp1_obj.B1 = options->B1pp1;
	fobj->pp1_obj.B2 = options->B2pp1;
	fobj->pp1_obj.stg2_is_default = (options->B2pp1 == 0);
	fobj->pp1_obj.pp1_exponent = 0;
	fobj->pp1_obj.pp1_multiplier = 0;
	fobj->pp1_obj.pp1_tune_freq = 0;
    fobj->pp1_obj.lcg_state = hash64(lcg_rand_64(&fobj->lcg_state));


	// initialize stuff for ecm	
	fobj->ecm_obj.B1 = options->B1ecm;
	fobj->ecm_obj.B2 = options->B2ecm;
	fobj->ecm_obj.stg2_is_default = (options->B2ecm == 0);
	fobj->ecm_obj.sigma = options->sigma;
	fobj->ecm_obj.num_curves = 90;
	fobj->ecm_obj.curves_run = 0;
	fobj->ecm_obj.ecm_exponent = 0;
	fobj->ecm_obj.ecm_multiplier = 0;
	fobj->ecm_obj.ecm_tune_freq = 0;
    fobj->ecm_obj.bail_on_factor = 1;
    fobj->ecm_obj.save_b1 = options->saveB1;
    fobj->ecm_obj.rand_seed1 = fobj->seed1;
    fobj->ecm_obj.rand_seed2 = fobj->seed2;
    for (i = 0; i < (int)options->threads; i++)
    {
        fobj->ecm_obj.lcg_state[i] = 
            hash64(lcg_rand_64(&fobj->lcg_state));
    }

	// unlike ggnfs, ecm does not *require* external binaries.  
	// an empty string indicates the use of the built-in GMP-ECM hooks, while
	// a non-empty string (filled in by the user) will indicate the use of
	// an external binary
	strcpy(fobj->ecm_obj.ecm_path, options->ecm_path);
	fobj->ecm_obj.use_external = 0;
    fobj->ecm_obj.prefer_gmpecm = options->prefer_gmpecm;
    fobj->ecm_obj.ecm_ext_xover = options->ext_ecm_xover;

	// initialize stuff for squfof
	fobj->squfof_obj.num_factors = 0;	

	// initialize stuff for qs	
	fobj->qs_obj.gbl_override_B_flag = (options->siqsB > 0);
	fobj->qs_obj.gbl_override_B = options->siqsB;
    fobj->qs_obj.gbl_override_small_cutoff_flag = (options->siqsTFSm > 0);
    fobj->qs_obj.gbl_override_small_cutoff = options->siqsTFSm;
	fobj->qs_obj.gbl_override_blocks_flag = (options->siqsNB > 0);
	fobj->qs_obj.gbl_override_blocks = options->siqsNB;
	fobj->qs_obj.gbl_override_lpmult_flag = (options->siqsM > 0);
	fobj->qs_obj.gbl_override_lpmult = options->siqsM;
	fobj->qs_obj.gbl_override_rel_flag = (options->siqsR > 0);
	fobj->qs_obj.gbl_override_rel = options->siqsR;
	fobj->qs_obj.gbl_override_tf_flag = (options->siqsTF > 0);
	fobj->qs_obj.gbl_override_tf = options->siqsTF;
	fobj->qs_obj.gbl_override_time_flag = (options->siqsT > 0);
	fobj->qs_obj.gbl_override_time = options->siqsT;
    fobj->qs_obj.inmem_cutoff = options->inmem_cutoff;
    if ((options->siqsMFBD < 1.0) || (options->siqsMFBD > 2.0))
    {
        if (options->siqsMFBD < 1.0)
            options->siqsMFBD = 1.0;
        if (options->siqsMFBD > 2.0)
            options->siqsMFBD = 2.0;
        printf("DLP Tdiv exponent should be between 1.0 and 2.0; setting MFBD = %1.2f\n", 
            options->siqsMFBD);
    }
    if ((options->siqsMFBT < 2.0) || (options->siqsMFBT > 3.0))
    {
        if (options->siqsMFBT < 2.0)
            options->siqsMFBT = 2.0;
        if (options->siqsMFBT > 3.0)
            options->siqsMFBT = 3.0;
        printf("TLP Tdiv exponent should be between 2.0 and 3.0; setting MFBT = %1.2f\n",
            options->siqsMFBT);
    }
	fobj->qs_obj.gbl_override_mfbd = options->siqsMFBD;
	fobj->qs_obj.gbl_override_mfbt = options->siqsMFBT;
	fobj->qs_obj.gbl_override_lpb = options->siqsLPB;
    if (options->siqsBDiv < 1.0)
    {
        printf("Batch divisor cannot be less than 1.0; setting BDiv = 1.0\n");
        options->siqsBDiv = 1.0;
    }
    fobj->qs_obj.gbl_override_bdiv = (float)options->siqsBDiv;
    fobj->qs_obj.gbl_override_3lp_bat = options->siqsNobat;
    fobj->qs_obj.gbl_btarget = options->siqsBT;
	fobj->qs_obj.flags = 0;
	fobj->qs_obj.gbl_force_DLP = options->siqsForceDLP;
	fobj->qs_obj.gbl_force_TLP = options->siqsForceTLP;
	fobj->qs_obj.qs_exponent = 0;
	fobj->qs_obj.qs_multiplier = 0;
	fobj->qs_obj.qs_tune_freq = 0;
    fobj->qs_obj.no_small_cutoff_opt = options->no_opt;
	strcpy(fobj->qs_obj.siqs_savefile, options->qssave);

	init_lehman();

	// initialize stuff for trial division	
	fobj->div_obj.print = 0;	
	fobj->div_obj.limit = 10000;
	fobj->div_obj.fmtlimit = options->fermat_max;
	
	//initialize stuff for nfs
	fobj->nfs_obj.snfs = 0;
	fobj->nfs_obj.gnfs = options->force_gnfs;
	fobj->nfs_obj.gnfs_exponent = 0;
	fobj->nfs_obj.gnfs_multiplier = 0;
	fobj->nfs_obj.gnfs_tune_freq = 0;
	fobj->nfs_obj.min_digits = 85;
    // raise min_rels bounds by a percentage
    // on unsuccessful filtering
	fobj->nfs_obj.filter_min_rels_nudge = 1.0 + options->filt_bump / 100.0;	
    
	fobj->nfs_obj.siever = options->ggnfs_siever;					
	fobj->nfs_obj.startq = options->sieveQstart;					
	fobj->nfs_obj.rangeq = options->sieveQstop;					
    if (fobj->nfs_obj.rangeq < fobj->nfs_obj.startq)
    {
        printf("ignoring sieve stop < sieve start\n");
    }
    else
    {
        fobj->nfs_obj.rangeq = fobj->nfs_obj.rangeq - fobj->nfs_obj.startq;
    }
	fobj->nfs_obj.polystart = options->polystart;
	fobj->nfs_obj.polyrange = options->polystop;
    if (fobj->nfs_obj.rangeq < fobj->nfs_obj.startq)
    {
        printf("ignoring poly stop < poly start\n");
    }
    else
    {
        fobj->nfs_obj.polyrange = fobj->nfs_obj.polyrange - fobj->nfs_obj.polystart;
    }

    char tmp[MAXARGLEN];
    char* cptr;
    strcpy(fobj->nfs_obj.outputfile, options->nfs_outfile);
    strcpy(tmp, fobj->nfs_obj.outputfile);
    cptr = strchr(tmp, 46);
    if (cptr == NULL)
    {
        //no . in provided filename
        sprintf(fobj->nfs_obj.outputfile, "%s.dat", fobj->nfs_obj.outputfile);
        sprintf(fobj->nfs_obj.logfile, "%s.log", fobj->nfs_obj.outputfile);
        sprintf(fobj->nfs_obj.fbfile, "%s.fb", fobj->nfs_obj.outputfile);
    }
    else
    {
        cptr[0] = '\0';
        sprintf(fobj->nfs_obj.logfile, "%s.log", tmp);
        sprintf(fobj->nfs_obj.fbfile, "%s.fb", tmp);
    }

    if (options->rat_side)
	    fobj->nfs_obj.sq_side = -1;					
    else if (options->alg_side)
        fobj->nfs_obj.sq_side = 1;					
    else
        fobj->nfs_obj.sq_side = 0;					// default: choose what makes sense

	fobj->nfs_obj.timeout = options->nfs_timeout;
	strcpy(fobj->nfs_obj.job_infile, options->nfs_jobfile);
	
    // default = fast search
    fobj->nfs_obj.poly_option = 0;					
    if (strcmp(options->poly_method, "wide") == 0)
        fobj->nfs_obj.poly_option = 1;
    else if (strcmp(options->poly_method, "deep") == 0)
        fobj->nfs_obj.poly_option = 2;
    else if (strcmp(options->poly_method, "fast") == 0)
        fobj->nfs_obj.poly_option = 0;
    else if (strcmp(options->poly_method, "min") == 0)
        fobj->nfs_obj.poly_option = 3;
    else if (strcmp(options->poly_method, "avg") == 0)
        fobj->nfs_obj.poly_option = 4;
    else if (strcmp(options->poly_method, "good") == 0)
        fobj->nfs_obj.poly_option = 5;

	fobj->nfs_obj.restart_flag = options->nfs_resume;
	fobj->nfs_obj.nfs_phases = NFS_DEFAULT_PHASES;
	fobj->nfs_obj.snfs_testsieve_threshold = options->snfs_testsieve_threshold;
	*fobj->nfs_obj.filearg = '\0';

	fobj->nfs_obj.polybatch = options->poly_batch;
	strcpy(fobj->nfs_obj.ggnfs_dir, options->ggnfs_dir);

	// initialize autofactor object
	// whether we want to output certain info to their own files...
    fobj->autofact_obj.want_output_primes = 0;
    fobj->autofact_obj.want_output_factors = 0;
    fobj->autofact_obj.want_output_unfactored = 0;
    fobj->autofact_obj.want_output_expressions = options->want_output_expr;
    if (strlen(options->opfile) > 0)
    {
        fobj->autofact_obj.want_output_primes = 1;
    }
    if (strlen(options->offile) > 0)
    {
        fobj->autofact_obj.want_output_factors = 1;
    }
    if (strlen(options->oufile) > 0)
    {
        fobj->autofact_obj.want_output_unfactored = 1;
    }
	
	fobj->autofact_obj.qs_gnfs_xover = options->xover;
    fobj->autofact_obj.qs_snfs_xover = options->qs_snfs_xover;
	// use xover even when timing info is available
    fobj->autofact_obj.prefer_xover = (options->xover != 95);
    fobj->autofact_obj.want_only_1_factor = options->one_factor;
	fobj->autofact_obj.no_ecm = options->no_ecm;
	fobj->autofact_obj.target_pretest_ratio = options->pretest_ratio;
    fobj->autofact_obj.initial_work = options->work;
	fobj->autofact_obj.has_snfs_form = -1;		// not checked yet

	//pretesting plan used by factor()
	fobj->autofact_obj.yafu_pretest_plan = PRETEST_NORMAL;
	strcpy(fobj->autofact_obj.plan_str, options->fact_plan);
    
    if (strcmp(options->fact_plan, "none") == 0)
        fobj->autofact_obj.yafu_pretest_plan = PRETEST_NONE;
    else if (strcmp(options->fact_plan, "noecm") == 0)
        fobj->autofact_obj.yafu_pretest_plan = PRETEST_NOECM;
    else if (strcmp(options->fact_plan, "light") == 0)
        fobj->autofact_obj.yafu_pretest_plan = PRETEST_LIGHT;
    else if (strcmp(options->fact_plan, "normal") == 0)
        fobj->autofact_obj.yafu_pretest_plan = PRETEST_NORMAL;
    else if (strcmp(options->fact_plan, "deep") == 0)
        fobj->autofact_obj.yafu_pretest_plan = PRETEST_DEEP;
    else if (strcmp(options->fact_plan, "custom") == 0)
        fobj->autofact_obj.yafu_pretest_plan = PRETEST_CUSTOM;
	fobj->autofact_obj.only_pretest = options->pretest;
	fobj->autofact_obj.autofact_active = 0;

	// if a number is <= aprcl_prove_cutoff, we will prove it prime or composite
	fobj->factors->aprcl_prove_cutoff = options->aprcl_p;
	// if a number is >= aprcl_display_cutoff, we will show the APRCL progress
	fobj->factors->aprcl_display_cutoff = options->aprcl_d;

    fobj->VFLAG = options->verbosity;
    if (strlen(options->factorlog) == 0)
    {
        fobj->LOGFLAG = 0;
    }
    else
    {
        fobj->LOGFLAG = 1;
    }

	return;
}

void free_factobj(fact_obj_t *fobj)
{
	// free stuff in rho
	free(fobj->rho_obj.polynomials);
	mpz_clear(fobj->rho_obj.gmp_n);
	mpz_clear(fobj->rho_obj.gmp_f);

	// free stuff in pm1
	mpz_clear(fobj->pm1_obj.gmp_n);
	mpz_clear(fobj->pm1_obj.gmp_f);

	// free stuff in pp1
	mpz_clear(fobj->pp1_obj.gmp_n);
	mpz_clear(fobj->pp1_obj.gmp_f);

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

void add_to_factor_list(yfactor_list_t *flist, mpz_t n, int VFLAG, int NUM_WITNESSES)
{
	//stick the number n into the provided factor list
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
            return;
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
	return;
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


