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

    // initialize general stuff in fobj
    get_random_seeds(&seed1, &seed2);
    fobj->seed1 = seed1;
    fobj->seed2 = seed2;
    fobj->lcg_state = ((uint64_t)seed2 << 32) | (uint64_t)seed1;
    fobj->flags = 0;
    fobj->num_threads = 1;
    strcpy(fobj->flogname, "factor.log");
    strcpy(fobj->factor_json_name, "factor.json");
    fobj->do_logging = 1;   // not used...
    fobj->LOGFLAG = 1;
    fobj->NUM_WITNESSES = 1;
    fobj->OBASE = 10;
    fobj->refactor_depth = 0;
    fobj->OUTPUT_TERSE = 0;
    fobj->VFLAG = 0;

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
    fobj->pm1_obj.ttime = 0.0;

    // initialize stuff for pp1	
    fobj->pp1_obj.B1 = 20000;
    fobj->pp1_obj.B2 = 1000000;
    fobj->pp1_obj.stg2_is_default = 1;
    fobj->pp1_obj.pp1_exponent = 0;
    fobj->pp1_obj.pp1_multiplier = 0;
    fobj->pp1_obj.pp1_tune_freq = 0;
    fobj->pp1_obj.vecnum = 0;
    fobj->pp1_obj.ttime = 0.0;

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
    fobj->ecm_obj.use_cgbn = 0;
    fobj->ecm_obj.gpucurves = 0;
    fobj->ecm_obj.use_gpuecm = 0;
    fobj->ecm_obj.use_gpudev = 0;
    fobj->ecm_obj.ttime = 0.0;


    // unlike ggnfs, ecm does not *require* external binaries.  
    // an empty string indicates the use of the built-in GMP-ECM hooks, while
    // a non-empty string (filled in by the user) will indicate the use of
    // an external binary
    strcpy(fobj->ecm_obj.ecm_path, "");
    fobj->ecm_obj.use_external = 0;
#ifdef USE_AVX512F
    fobj->ecm_obj.prefer_gmpecm = 0;
    fobj->ecm_obj.ecm_ext_xover = 300000000;
#else
    fobj->ecm_obj.prefer_gmpecm = 1;
    fobj->ecm_obj.ecm_ext_xover = 48000;
#endif
    fobj->ecm_obj.prefer_gmpecm_stg2 = 0;
    fobj->ecm_obj.prefer_avxecm_stg2 = 0;

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
    fobj->qs_obj.gbl_override_ssidx_flag = 0;
    fobj->qs_obj.gbl_override_ssidx = 0;
    fobj->qs_obj.gbl_override_ssidx_flag = 0;
    fobj->qs_obj.gbl_override_ssidx = 12;
    fobj->qs_obj.gbl_btarget = 1000000;
    fobj->qs_obj.flags = 0;
    fobj->qs_obj.gbl_force_DLP = 0;
    fobj->qs_obj.gbl_force_TLP = 0;
    fobj->qs_obj.gbl_force_QLP = 0;
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
    fobj->nfs_obj.min_digits = 80;
    fobj->nfs_obj.filter_min_rels_nudge = 1.05;	// raise min_rels bounds by 5%
                                                    // on unsuccessful filtering
    fobj->nfs_obj.siever = 0;					//default, use automatic selection
    fobj->nfs_obj.startq = 0;					//default, not used
    fobj->nfs_obj.rangeq = 0;					//default, not used
    fobj->nfs_obj.polystart = 0;					//default, not used
    fobj->nfs_obj.polyrange = 0;					//default, not used
    fobj->nfs_obj.np1 = 0;
    fobj->nfs_obj.nps = 0;
    fobj->nfs_obj.npr = 0;
    strcpy(fobj->nfs_obj.params_file, "");          // default: use built-in table
    strcpy(fobj->nfs_obj.outputfile, "nfs.dat");			//default
    strcpy(fobj->nfs_obj.logfile, "nfs.log");			//default
    strcpy(fobj->nfs_obj.fbfile, "nfs.fb");				//default
    fobj->nfs_obj.sq_side = 0;					//default = algebraic
    fobj->nfs_obj.timeout = 1<<31;					//default, not used
    strcpy(fobj->nfs_obj.job_infile, "nfs.job");			//default
    fobj->nfs_obj.poly_option = 0;					//default = fast search
                                        //1 = wide
                                    //2 = deep
    fobj->nfs_obj.restart_flag = 0;					//default = not a restart
    fobj->nfs_obj.nfs_phases = NFS_DEFAULT_PHASES;
    fobj->nfs_obj.snfs_testsieve_threshold = 160;
    *fobj->nfs_obj.filearg = '\0';
    fobj->nfs_obj.minrels = 0;                      // 0 = default stopping point (heuristic)
    fobj->nfs_obj.poly_testsieve = 390;
    fobj->nfs_obj.poly_percent_max = 8;
    fobj->nfs_obj.td = 0;
    fobj->nfs_obj.batch_3lp = 0;
    mpz_set_ui(fobj->nfs_obj.snfs_cofactor, 0);

    fobj->nfs_obj.polybatch = 250;						//default	
#if defined(_WIN64)
    strcpy(fobj->nfs_obj.ggnfs_dir, ".\\");
    strcpy(fobj->nfs_obj.cado_dir, ".\\");
#elif defined(WIN32)
    strcpy(fobj->nfs_obj.ggnfs_dir, ".\\");
    strcpy(fobj->nfs_obj.cado_dir, ".\\");
#else
    strcpy(fobj->nfs_obj.ggnfs_dir, "./");
    strcpy(fobj->nfs_obj.cado_dir, "./");
#endif

    fobj->nfs_obj.poly_time = 0.0;
    fobj->nfs_obj.filter_time = 0.0;
    fobj->nfs_obj.sieve_time = 0.0;
    fobj->nfs_obj.la_time = 0.0;
    fobj->nfs_obj.sqrt_time = 0.0;
    fobj->nfs_obj.ttime = 0.0;

    fobj->nfs_obj.cadoMsieve = 0;
    strcpy(fobj->nfs_obj.convert_poly_path, "");
    fobj->nfs_obj.skip_snfs_check = 0;

    //initialize autofactor object
    //whether we want to output certain info to their own files...
    fobj->autofact_obj.want_output_primes = 0;
    fobj->autofact_obj.want_output_factors = 0;
    fobj->autofact_obj.want_output_unfactored = 0;
    fobj->autofact_obj.want_output_expressions = 1;
    fobj->autofact_obj.qs_gnfs_xover = 95;
    fobj->autofact_obj.qs_snfs_xover = 95;
    fobj->autofact_obj.prefer_xover = 0;
    fobj->autofact_obj.want_only_1_factor = 0;
    fobj->autofact_obj.no_ecm = 0;
    fobj->autofact_obj.target_pretest_ratio = 4.0 / 13.0;
    fobj->autofact_obj.initial_work = 0.0;
    fobj->autofact_obj.has_snfs_form = -1;		// not checked yet
    fobj->autofact_obj.stopbase = 10;
    fobj->autofact_obj.stopeq = -1;
    fobj->autofact_obj.stoplt = -1;
    fobj->autofact_obj.stople = -1;
    fobj->autofact_obj.stopgt = -1;
    fobj->autofact_obj.stopge = -1;
    fobj->autofact_obj.stopprime = 0;
    fobj->autofact_obj.check_stop_conditions = 0;
    fobj->autofact_obj.max_siqs = -1;
    fobj->autofact_obj.max_nfs = -1;

    //pretesting plan used by factor()
    fobj->autofact_obj.yafu_pretest_plan = PRETEST_NORMAL;
    strcpy(fobj->autofact_obj.plan_str, "normal");
    fobj->autofact_obj.only_pretest = 0;
    fobj->autofact_obj.autofact_active = 0;
    fobj->autofact_obj.json_pretty = 0;

    // if a number is <= aprcl_prove_cutoff, we will prove it prime or composite
    fobj->factors->aprcl_prove_cutoff = 500;
    // if a number is >= aprcl_display_cutoff, we will show the APRCL progress
    fobj->factors->aprcl_display_cutoff = 200;
    fobj->factors->total_factors = 0;

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
    mpz_clear(fobj->input_N);
    free(fobj->input_str);
    fobj->input_str_alloc = 0;

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
    mpz_init(fobj->input_N);
    mpz_set_ui(fobj->N, 0);
    mpz_set_ui(fobj->input_N, 0);
    fobj->input_str = (char*)xmalloc(1024 * sizeof(char));
    fobj->input_str_alloc = 1024;

	fobj->rho_obj.num_poly = 3;
	fobj->rho_obj.polynomials = (uint32_t *)xmalloc(fobj->rho_obj.num_poly * sizeof(uint32_t));
	fobj->rho_obj.polynomials[0] = 3;
	fobj->rho_obj.polynomials[1] = 2;
	fobj->rho_obj.polynomials[2] = 1;
	mpz_init(fobj->rho_obj.gmp_n);
	mpz_init(fobj->rho_obj.gmp_f);
    mpz_set_ui(fobj->rho_obj.gmp_n, 0);
    mpz_set_ui(fobj->rho_obj.gmp_f, 0);

	mpz_init(fobj->pm1_obj.gmp_n);
	mpz_init(fobj->pm1_obj.gmp_f);
    mpz_set_ui(fobj->pm1_obj.gmp_n, 0);
    mpz_set_ui(fobj->pm1_obj.gmp_f, 0);
	
	mpz_init(fobj->pp1_obj.gmp_n);
	mpz_init(fobj->pp1_obj.gmp_f);
    mpz_set_ui(fobj->pp1_obj.gmp_n, 0);
    mpz_set_ui(fobj->pp1_obj.gmp_f, 0);

	mpz_init(fobj->ecm_obj.gmp_n);
	mpz_init(fobj->ecm_obj.gmp_f);
    mpz_set_ui(fobj->ecm_obj.gmp_n, 0);
    mpz_set_ui(fobj->ecm_obj.gmp_f, 0);
    fobj->ecm_obj.lcg_state = (uint64_t*)xmalloc(fobj->num_threads * sizeof(uint64_t));

	mpz_init(fobj->squfof_obj.gmp_n);
	mpz_init(fobj->squfof_obj.gmp_f);
    mpz_set_ui(fobj->squfof_obj.gmp_n, 0);
    mpz_set_ui(fobj->squfof_obj.gmp_f, 0);

	mpz_init(fobj->qs_obj.gmp_n);
    mpz_set_ui(fobj->qs_obj.gmp_n, 0);

	mpz_init(fobj->div_obj.gmp_n);
	mpz_init(fobj->div_obj.gmp_f);
    mpz_set_ui(fobj->div_obj.gmp_n, 0);
    mpz_set_ui(fobj->div_obj.gmp_f, 0);

    mpz_init(fobj->nfs_obj.full_n);
	mpz_init(fobj->nfs_obj.gmp_n);
	mpz_init(fobj->nfs_obj.snfs_cofactor);
    mpz_set_ui(fobj->nfs_obj.full_n, 0);
    mpz_set_ui(fobj->nfs_obj.gmp_n, 0);
    mpz_set_ui(fobj->nfs_obj.snfs_cofactor, 0);

    fobj->factors = (yfactor_list_t*)xmalloc(1 * sizeof(yfactor_list_t));
	fobj->factors->aprcl_display_cutoff = 200;
	fobj->factors->aprcl_prove_cutoff = 500;
	fobj->factors->alloc_factors = 256;
	fobj->factors->factors = (yfactor_t *)xmalloc(256 * sizeof(yfactor_t));
	for (i=0; i < fobj->factors->alloc_factors; i++)
	{
		fobj->factors->factors[i].type = UNKNOWN;
		fobj->factors->factors[i].count = 0;
	}
    fobj->factors->factorization_str = NULL;
    fobj->factors->str_alloc = 0;

	fobj->ecm_obj.num_factors = 0;	
	fobj->qs_obj.num_factors = 0;	
	fobj->div_obj.num_factors = 0;	
	fobj->nfs_obj.num_factors = 0;	
	fobj->factors->num_factors = 0;
    fobj->factors->total_factors = 0;

    soe_staticdata_t* sdata;
    // put a list of primes in the fobj; many algorithms use it
    sdata = soe_init(0, 1, 32768);
    fobj->primes = soe_wrapper(sdata, 0, 10000000, 0, &fobj->num_p, 0, 0);
    fobj->min_p = 2;
    fobj->max_p = fobj->primes[fobj->num_p - 1];
    soe_finalize(sdata);

	return;
}

void copy_factobj(fact_obj_t* dest, fact_obj_t* src)
{
    int i;

    dest->seed1 = src->seed1;
    dest->seed2 = src->seed2;
    dest->lcg_state = src->lcg_state;
    dest->flags = src->flags;
    dest->num_threads = src->num_threads;
    strcpy(dest->flogname, src->flogname);
    dest->do_logging = src->do_logging;   // not used...
    dest->LOGFLAG = src->LOGFLAG;
    dest->NUM_WITNESSES = src->NUM_WITNESSES;
    dest->OBASE = src->OBASE;
    dest->VFLAG = src->VFLAG;
    dest->OUTPUT_TERSE = src->OUTPUT_TERSE;
    dest->THREADS = src->THREADS;
    dest->LATHREADS = src->LATHREADS;

    copy_factor_list(dest->factors, src->factors);

    // initialize stuff for rho	
    dest->rho_obj.iterations = src->rho_obj.iterations;
    dest->rho_obj.curr_poly = src->rho_obj.curr_poly;

    // initialize stuff for pm1	
    mpz_set(dest->pm1_obj.gmp_n, src->pm1_obj.gmp_n);
    dest->pm1_obj.B1 = src->pm1_obj.B1;
    dest->pm1_obj.B2 = src->pm1_obj.B2;
    dest->pm1_obj.stg2_is_default = src->pm1_obj.stg2_is_default;
    dest->pm1_obj.pm1_exponent = src->pm1_obj.pm1_exponent;
    dest->pm1_obj.pm1_multiplier = src->pm1_obj.pm1_multiplier;
    dest->pm1_obj.pm1_tune_freq = src->pm1_obj.pm1_tune_freq;
    dest->pm1_obj.vecnum = src->pm1_obj.vecnum;

    // initialize stuff for pp1	
    mpz_set(dest->pp1_obj.gmp_n, src->pp1_obj.gmp_n);
    dest->pp1_obj.B1 = src->pp1_obj.B1;
    dest->pp1_obj.B2 = src->pp1_obj.B2;
    dest->pp1_obj.stg2_is_default = src->pp1_obj.stg2_is_default;
    dest->pp1_obj.pp1_exponent = src->pp1_obj.pp1_exponent;
    dest->pp1_obj.pp1_multiplier = src->pp1_obj.pp1_multiplier;
    dest->pp1_obj.pp1_tune_freq = src->pp1_obj.pp1_tune_freq;
    dest->pp1_obj.vecnum = src->pp1_obj.vecnum;

    // initialize stuff for ecm	
    mpz_set(dest->ecm_obj.gmp_n, src->ecm_obj.gmp_n);
    dest->ecm_obj.B1 = src->ecm_obj.B1;
    dest->ecm_obj.B2 = src->ecm_obj.B2;
    dest->ecm_obj.stg2_is_default = src->ecm_obj.stg2_is_default;
    dest->ecm_obj.sigma = src->ecm_obj.sigma;
    dest->ecm_obj.num_curves = src->ecm_obj.num_curves;
    dest->ecm_obj.curves_run = src->ecm_obj.curves_run;
    dest->ecm_obj.ecm_exponent = src->ecm_obj.ecm_exponent;
    dest->ecm_obj.ecm_multiplier = src->ecm_obj.ecm_multiplier;
    dest->ecm_obj.ecm_tune_freq = src->ecm_obj.ecm_tune_freq;
    dest->ecm_obj.bail_on_factor = src->ecm_obj.bail_on_factor;
    dest->ecm_obj.save_b1 = src->ecm_obj.save_b1;
    dest->ecm_obj.gpucurves = src->ecm_obj.gpucurves;
    dest->ecm_obj.use_cgbn = src->ecm_obj.use_cgbn;
    dest->ecm_obj.use_gpudev = src->ecm_obj.use_gpudev;
    dest->ecm_obj.use_gpuecm = src->ecm_obj.use_gpuecm;

    // unlike ggnfs, ecm does not *require* external binaries.  
    // an empty string indicates the use of the built-in GMP-ECM hooks, while
    // a non-empty string (filled in by the user) will indicate the use of
    // an external binary
    strcpy(dest->ecm_obj.ecm_path, src->ecm_obj.ecm_path);
    dest->ecm_obj.use_external = src->ecm_obj.use_external;
#ifdef USE_AVX512F
    dest->ecm_obj.prefer_gmpecm = src->ecm_obj.prefer_gmpecm;
    dest->ecm_obj.ecm_ext_xover = src->ecm_obj.ecm_ext_xover;
#else
    dest->ecm_obj.prefer_gmpecm = src->ecm_obj.prefer_gmpecm;
    dest->ecm_obj.ecm_ext_xover = src->ecm_obj.ecm_ext_xover;
#endif

    dest->ecm_obj.lcg_state = (uint64_t*)xrealloc(dest->ecm_obj.lcg_state,
        src->num_threads * sizeof(uint64_t));
    for (i = 0; i < (int)src->num_threads; i++)
    {
        dest->ecm_obj.lcg_state[i] =
            hash64(lcg_rand_64(&dest->lcg_state));
    }


    // initialize stuff for squfof
    dest->squfof_obj.num_factors = src->squfof_obj.num_factors;

    // initialize stuff for qs	
    mpz_set(dest->qs_obj.gmp_n, src->qs_obj.gmp_n);
    dest->qs_obj.gbl_override_B_flag = src->qs_obj.gbl_override_B_flag;
    dest->qs_obj.gbl_override_B = src->qs_obj.gbl_override_B;
    dest->qs_obj.gbl_override_small_cutoff_flag = src->qs_obj.gbl_override_small_cutoff_flag;
    dest->qs_obj.gbl_override_small_cutoff = src->qs_obj.gbl_override_small_cutoff;
    dest->qs_obj.gbl_override_blocks_flag = src->qs_obj.gbl_override_blocks_flag;
    dest->qs_obj.gbl_override_blocks = src->qs_obj.gbl_override_blocks;
    dest->qs_obj.gbl_override_lpmult_flag = src->qs_obj.gbl_override_lpmult_flag;
    dest->qs_obj.gbl_override_lpmult = src->qs_obj.gbl_override_lpmult;
    dest->qs_obj.gbl_override_rel_flag = src->qs_obj.gbl_override_rel_flag;
    dest->qs_obj.gbl_override_rel = src->qs_obj.gbl_override_rel;
    dest->qs_obj.gbl_override_tf_flag = src->qs_obj.gbl_override_tf_flag;
    dest->qs_obj.gbl_override_tf = src->qs_obj.gbl_override_tf;
    dest->qs_obj.gbl_override_time_flag = src->qs_obj.gbl_override_time_flag;
    dest->qs_obj.gbl_override_time = src->qs_obj.gbl_override_time;
    dest->qs_obj.gbl_override_mfbd = src->qs_obj.gbl_override_mfbd;
    dest->qs_obj.gbl_override_mfbt = src->qs_obj.gbl_override_mfbt;
    dest->qs_obj.gbl_override_lpb = src->qs_obj.gbl_override_lpb;
    dest->qs_obj.gbl_override_bdiv_flag = src->qs_obj.gbl_override_bdiv_flag;
    dest->qs_obj.gbl_override_bdiv = src->qs_obj.gbl_override_bdiv;
    dest->qs_obj.gbl_override_3lp_bat = src->qs_obj.gbl_override_3lp_bat;
    dest->qs_obj.gbl_btarget = src->qs_obj.gbl_btarget;
    dest->qs_obj.gbl_override_ssidx_flag = src->qs_obj.gbl_override_ssidx_flag;
    dest->qs_obj.gbl_override_ssidx = src->qs_obj.gbl_override_ssidx;
    dest->qs_obj.gbl_override_ssalloc_flag = src->qs_obj.gbl_override_ssalloc_flag;
    dest->qs_obj.gbl_override_ssalloc = src->qs_obj.gbl_override_ssalloc;
    dest->qs_obj.flags = src->qs_obj.flags;
    dest->qs_obj.gbl_force_DLP = src->qs_obj.gbl_force_DLP;
    dest->qs_obj.gbl_force_TLP = src->qs_obj.gbl_force_TLP;
    dest->qs_obj.qs_exponent = src->qs_obj.qs_exponent;
    dest->qs_obj.qs_multiplier = src->qs_obj.qs_multiplier;
    dest->qs_obj.qs_tune_freq = src->qs_obj.qs_tune_freq;
    dest->qs_obj.no_small_cutoff_opt = src->qs_obj.no_small_cutoff_opt;
    strcpy(dest->qs_obj.siqs_savefile, src->qs_obj.siqs_savefile);
    init_lehman();

    // initialize stuff for trial division	
    mpz_set(dest->div_obj.gmp_n, src->div_obj.gmp_n);
    dest->div_obj.print = src->div_obj.print;
    dest->div_obj.limit = src->div_obj.limit;
    dest->div_obj.fmtlimit = src->div_obj.fmtlimit;

    //initialize stuff for nfs
    mpz_set(dest->nfs_obj.gmp_n, src->nfs_obj.gmp_n);
    mpz_set(dest->nfs_obj.full_n, src->nfs_obj.full_n);
    mpz_set(dest->nfs_obj.snfs_cofactor, src->nfs_obj.snfs_cofactor);
    dest->nfs_obj.snfs = src->nfs_obj.snfs;
    dest->nfs_obj.gnfs = src->nfs_obj.gnfs;
    dest->nfs_obj.gnfs_exponent = src->nfs_obj.gnfs_exponent;
    dest->nfs_obj.gnfs_multiplier = src->nfs_obj.gnfs_multiplier;
    dest->nfs_obj.gnfs_tune_freq = src->nfs_obj.gnfs_tune_freq;
    dest->nfs_obj.min_digits = src->nfs_obj.min_digits;
    dest->nfs_obj.filter_min_rels_nudge = src->nfs_obj.filter_min_rels_nudge;
    dest->nfs_obj.siever = src->nfs_obj.siever;
    dest->nfs_obj.startq = src->nfs_obj.startq;
    dest->nfs_obj.rangeq = src->nfs_obj.rangeq;
    dest->nfs_obj.polystart = src->nfs_obj.polystart;
    dest->nfs_obj.polyrange = src->nfs_obj.polyrange;
    dest->nfs_obj.np1 = src->nfs_obj.np1;
    dest->nfs_obj.nps = src->nfs_obj.nps;
    dest->nfs_obj.npr = src->nfs_obj.npr;
    strcpy(dest->nfs_obj.params_file, src->nfs_obj.params_file);
    strcpy(dest->nfs_obj.outputfile, src->nfs_obj.outputfile);
    strcpy(dest->nfs_obj.logfile, src->nfs_obj.logfile);
    strcpy(dest->nfs_obj.fbfile, src->nfs_obj.fbfile);
    dest->nfs_obj.sq_side = src->nfs_obj.sq_side;
    dest->nfs_obj.timeout = src->nfs_obj.timeout;
    strcpy(dest->nfs_obj.job_infile, src->nfs_obj.job_infile);
    dest->nfs_obj.poly_option = src->nfs_obj.poly_option;
    dest->nfs_obj.restart_flag = src->nfs_obj.restart_flag;
    dest->nfs_obj.nfs_phases = src->nfs_obj.nfs_phases;
    dest->nfs_obj.snfs_testsieve_threshold = src->nfs_obj.snfs_testsieve_threshold;
    strcpy(dest->nfs_obj.filearg, src->nfs_obj.filearg);
    dest->nfs_obj.skip_snfs_check = src->nfs_obj.skip_snfs_check;
    dest->nfs_obj.poly_percent_max = src->nfs_obj.poly_percent_max;
    dest->nfs_obj.poly_testsieve = src->nfs_obj.poly_testsieve;
    dest->nfs_obj.td = src->nfs_obj.td;
    dest->nfs_obj.batch_3lp = src->nfs_obj.batch_3lp;

    dest->nfs_obj.polybatch = src->nfs_obj.polybatch;
    strcpy(dest->nfs_obj.ggnfs_dir, src->nfs_obj.ggnfs_dir);

    dest->nfs_obj.cadoMsieve = src->nfs_obj.cadoMsieve;
    strcpy(dest->nfs_obj.cado_dir, src->nfs_obj.cado_dir);
    strcpy(dest->nfs_obj.convert_poly_path, src->nfs_obj.convert_poly_path);

    //initialize autofactor object
    //whether we want to output certain info to their own files...
    dest->autofact_obj.want_output_primes = src->autofact_obj.want_output_primes;
    dest->autofact_obj.want_output_factors = src->autofact_obj.want_output_factors;
    dest->autofact_obj.want_output_unfactored = src->autofact_obj.want_output_unfactored;
    dest->autofact_obj.want_output_expressions = src->autofact_obj.want_output_expressions;

    if (src->autofact_obj.want_output_primes && 
        (strlen(src->autofact_obj.op_str) > 0) &&
        (strlen(src->autofact_obj.op_str) < 256))
    {
        strcpy(dest->autofact_obj.op_str, src->autofact_obj.op_str);
    }

    if (src->autofact_obj.want_output_factors &&
        (strlen(src->autofact_obj.of_str) > 0) &&
        (strlen(src->autofact_obj.of_str) < 256))
    {
        strcpy(dest->autofact_obj.of_str, src->autofact_obj.of_str);
    }

    if (src->autofact_obj.want_output_unfactored &&
        (strlen(src->autofact_obj.ou_str) > 0) &&
        (strlen(src->autofact_obj.ou_str) < 256))
    {
        strcpy(dest->autofact_obj.ou_str, src->autofact_obj.ou_str);
    }

    dest->autofact_obj.qs_gnfs_xover = src->autofact_obj.qs_gnfs_xover;
    dest->autofact_obj.qs_snfs_xover = src->autofact_obj.qs_snfs_xover;
    // use xover even when timing info is available
    dest->autofact_obj.prefer_xover = src->autofact_obj.prefer_xover;
    dest->autofact_obj.want_only_1_factor = src->autofact_obj.want_only_1_factor;
    dest->autofact_obj.no_ecm = src->autofact_obj.no_ecm;
    dest->autofact_obj.target_pretest_ratio = src->autofact_obj.target_pretest_ratio;
    dest->autofact_obj.initial_work = src->autofact_obj.initial_work;
    dest->autofact_obj.has_snfs_form = src->autofact_obj.has_snfs_form;
    dest->autofact_obj.stopbase = src->autofact_obj.stopbase;
    dest->autofact_obj.stopge = src->autofact_obj.stopge;
    dest->autofact_obj.stopgt = src->autofact_obj.stopgt;
    dest->autofact_obj.stopeq = src->autofact_obj.stopeq;
    dest->autofact_obj.stoplt = src->autofact_obj.stoplt;
    dest->autofact_obj.stople = src->autofact_obj.stople;
    dest->autofact_obj.stopprime = src->autofact_obj.stopprime;
    dest->autofact_obj.check_stop_conditions = src->autofact_obj.check_stop_conditions;
    dest->autofact_obj.max_siqs = src->autofact_obj.max_siqs;
    dest->autofact_obj.max_nfs = src->autofact_obj.max_siqs;
    dest->autofact_obj.initial_work = src->autofact_obj.initial_work;
    dest->autofact_obj.ecm_total_work_performed = src->autofact_obj.ecm_total_work_performed;

    //pretesting plan used by factor()
    dest->autofact_obj.yafu_pretest_plan = src->autofact_obj.yafu_pretest_plan;
    strcpy(dest->autofact_obj.plan_str, src->autofact_obj.plan_str);
    dest->autofact_obj.only_pretest = src->autofact_obj.only_pretest;
    dest->autofact_obj.autofact_active = src->autofact_obj.autofact_active;

    // if a number is <= aprcl_prove_cutoff, we will prove it prime or composite
    dest->factors->aprcl_prove_cutoff = src->factors->aprcl_prove_cutoff;
    // if a number is >= aprcl_display_cutoff, we will show the APRCL progress
    dest->factors->aprcl_display_cutoff = src->factors->aprcl_display_cutoff;

    dest->MEAS_CPU_FREQUENCY = 42;  // not used anymore
    strcpy(dest->CPU_ID_STR, src->CPU_ID_STR);
    dest->HAS_AVX2 = src->HAS_AVX2;
    dest->HAS_AVX = src->HAS_AVX;
    dest->HAS_SSE41 = src->HAS_SSE41;
    dest->NUM_WITNESSES = src->NUM_WITNESSES;
    dest->cache_size1 = src->cache_size1;
    dest->cache_size2 = src->cache_size2;
    dest->LOGFLAG = src->LOGFLAG;
    dest->THREADS = src->THREADS;
    dest->HAS_BMI2 = src->HAS_BMI2;
    dest->HAS_AVX512F = src->HAS_AVX512F;
    dest->HAS_AVX512BW = src->HAS_AVX512BW;
    dest->VFLAG = src->VFLAG;
    dest->LATHREADS = src->LATHREADS;

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

int find_in_factor_list(yfactor_list_t* flist, mpz_t n)
{
    int i;

    // look to see if this factor is already in the list
    for (i = 0; i < flist->num_factors; i++)
    {
        if (mpz_cmp(n, flist->factors[i].factor) == 0)
        {
            return i;
        }
    }
    return -1;
}

int add_to_factor_list(yfactor_list_t *flist, mpz_t n, int VFLAG, int NUM_WITNESSES)
{
	// stick the number n into the provided factor list.
    // return the index into which the factor was added.
	int i;
    int fid;
	int found = 0, v = 0;

    flist->total_factors++;

    //gmp_printf("adding %Zd to factor list\n", n);

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
    flist->factors[fid].type = UNKNOWN;

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

        //printf("aprtcle returned %d\n", ret);

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
    else if (NUM_WITNESSES > 0)
    {

        if (is_mpz_prp(n, NUM_WITNESSES))
        {
            if (mpz_cmp_ui(n, 100000000) < 0)
                flist->factors[fid].type = PRIME;
            else
                flist->factors[fid].type = PRP;
        }
        else
        {
            flist->factors[fid].type = COMPOSITE;
        }
    }
    else
    {
        flist->factors[fid].type = UNKNOWN;
    }

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
            flist->total_factors -= flist->factors[i].count;

            int j;
            // copy everything above this in the list back one position
            for (j = i; j < flist->num_factors - 1; j++)
            {
                mpz_set(flist->factors[j].factor, flist->factors[j + 1].factor);
                flist->factors[j].count = flist->factors[j + 1].count;
                flist->factors[j].type = flist->factors[j + 1].type;
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

void copy_factor_list(yfactor_list_t* dest, yfactor_list_t* src)
{
    // don't "add_to_factor_list" because it will re-run apr-cl or is_mpz_prp on the factors.

    return;
}

char* append_factor(char* in, mpz_t factor)
{
    int numb = (strlen(in) + gmp_base10(factor) + 4);
    char* out = (char*)xmalloc(numb * sizeof(char));
    if (out == NULL)
    {
        printf("could not allocate %d bytes for output buffer\n", numb);
        exit(1);
    }
    gmp_sprintf(out, "%s%Zd*", in, factor);
    free(in);
    return out;
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
    free(flist->factors);
    flist->num_factors = 0;
    flist->total_factors = 0;

    if ((flist->str_alloc > 0) && (flist->factorization_str != NULL))
    {
        free(flist->factorization_str);
    }

	return;
}

void generate_factorization_str(yfactor_list_t* flist)
{
    int i;
    int j;
    mpz_t prod;
    char* tersebuf = NULL;
    int maxlen = 0;

    if (flist->num_factors == 0)
    {
        if (flist->factorization_str != NULL)
        {
            free(flist->factorization_str);
        }

        flist->factorization_str = (char*)xmalloc(2 * sizeof(char));
        strcpy(flist->factorization_str, "");
        return;
    }

    mpz_init(prod);

    get_prod_of_factors(flist, prod);

    maxlen = 3 * gmp_base10(prod) + 4;
    tersebuf = (char*)xmalloc(maxlen);
    strcpy(tersebuf, "");

    for (i = 0; i < flist->num_factors; i++)
    {
        for (j = 0; j < flist->factors[i].count; j++)
        {
            tersebuf = append_factor(tersebuf, flist->factors[i].factor);
        }
    }

    if (flist->factorization_str != NULL)
    {
        free(flist->factorization_str);
    }

    flist->factorization_str = (char*)xmalloc((gmp_base10(prod) + strlen(tersebuf) + 10) * sizeof(char));
    
    flist->str_alloc = gmp_sprintf(flist->factorization_str, "%Zd=%s", prod, tersebuf) + 1;
    free(tersebuf);

    mpz_clear(prod);

    return;
}

void get_prod_of_factors(yfactor_list_t* flist, mpz_t prod)
{
    int i;
    mpz_set_ui(prod, 1);
    for (i = 0; i < flist->num_factors; i++)
    {
        int j;
        for (j = 0; j < flist->factors[i].count; j++)
        {
            mpz_mul(prod, prod, flist->factors[i].factor);
        }
    }

    return;
}

int contains_any_composite_factor(yfactor_list_t* flist, int num_witnesses)
{
    int p = 0;
    int i;

    for (i = 0; i < flist->num_factors; i++)
    {
        if (flist->factors[i].type == UNKNOWN)
        {
            int status = is_mpz_prp(flist->factors[i].factor, num_witnesses);
            if (status)
            {
                flist->factors[i].type = PRP;
            }
            else
            {
                flist->factors[i].type = COMPOSITE;
                p = 1;
                break;
            }
        }
        else if (flist->factors[i].type == COMPOSITE)
        {
            p = 1;
            break;
        }
    }

    return p;
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

void print_factor(const char* prefix, mpz_t factor, int OBASE)
{
    // print the size in base-10 digits and the factor in the requested base.
    // could also print the size in base-X digits but that doesn't seem right...
    if (OBASE == 8)
    {
        gmp_printf("%s%do = 0o%Zo\n", prefix, mpz_sizeinbase(factor, 8), factor);
    }
    else if (OBASE == 10)
    {
        gmp_printf("%s%d = %Zd\n", prefix, gmp_base10(factor), factor);
    }
    else if (OBASE == 16)
    {
        gmp_printf("%s%dx = 0x%Zx\n", prefix, mpz_sizeinbase(factor, 16), factor);
    }
    else if (OBASE == 2)
    {
        int i;
        printf("%s%db = 0b", prefix, mpz_sizeinbase(factor, 2));
        for (i = mpz_sizeinbase(factor, 2) - 1; i >= 0; i--)
        {
            printf("%d", mpz_tstbit(factor, i));
        }
        printf("\n");
    }
    return;
}



void print_factors(fact_obj_t *fobj, yfactor_list_t* flist, mpz_t N, int VFLAG, int NUM_WITNESSES, int OBASE)
{
	uint32_t i;
	int j, v;
	mpz_t tmp, tmp2;
    char* tersebuf = NULL;
    int maxlen = 0;
    int print_terse = fobj->OUTPUT_TERSE;

    if (print_terse)
    {    
        maxlen = MAX(gmp_base10(fobj->input_N) + 4, 1024);
        tersebuf = (char*)xmalloc(maxlen * 4);
        strcpy(tersebuf, "");
    }

	// always print factors unless complete silence is requested
	if ((VFLAG >= 0) || print_terse)
	{
		if (VFLAG >= 0) printf("\n\n***factors found***\n");
		mpz_init(tmp);
		mpz_set_ui(tmp, 1);

		for (i=0; i< flist->num_factors; i++)
		{
			if (flist->factors[i].type == PRIME)
			{
				// don't redo APR-CL calculations already performed by add_to_factor_list
				for (j=0;j< flist->factors[i].count;j++)
				{
					mpz_mul(tmp, tmp, flist->factors[i].factor);
                    if (print_terse) {
                        tersebuf = append_factor(tersebuf, flist->factors[i].factor);
                    }
                    print_factor("P", flist->factors[i].factor, OBASE);
				}
			}
			else if (flist->factors[i].type == PRP)
            {
                // don't redo mpz_isprobab_prime calculations already performed by add_to_factor_list
                for (j = 0; j < flist->factors[i].count; j++)
                {
                    mpz_mul(tmp, tmp, flist->factors[i].factor);
                    if (print_terse) {
                        tersebuf = append_factor(tersebuf, flist->factors[i].factor);
                    }
                    print_factor("PRP", flist->factors[i].factor, OBASE);
                }
            }
            else if (flist->factors[i].type == COMPOSITE)
            {
                // don't redo mpz_isprobab_prime calculations already performed by add_to_factor_list
                for (j = 0; j < flist->factors[i].count; j++)
                {
                    mpz_mul(tmp, tmp, flist->factors[i].factor);
                    if (print_terse) {
                        tersebuf = append_factor(tersebuf, flist->factors[i].factor);
                    }
                    print_factor("C", flist->factors[i].factor, OBASE);
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

                    //gmp_printf("in print_factors aprtcle returned status %d on input %Zd\n", ret, 
                    //    flist->factors[i].factor);

					if (v == APRTCLE_VERBOSE1)
						printf("\n");

					if (ret == APRTCLE_PRIME)
					{
						for (j=0;j< flist->factors[i].count;j++)
						{
							mpz_mul(tmp, tmp, flist->factors[i].factor);
							//gmp_printf("P%d = %Zd\n", 
                            //    gmp_base10(flist->factors[i].factor),
                            //    flist->factors[i].factor);
                            if (print_terse) {
                                tersebuf = append_factor(tersebuf, flist->factors[i].factor);
                            }
                            print_factor("P", flist->factors[i].factor, OBASE);
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
							//gmp_printf("C%d = %Zd\n", gmp_base10(flist->factors[i].factor),
                            //    flist->factors[i].factor);
                            if (print_terse) {
                                tersebuf = append_factor(tersebuf, flist->factors[i].factor);
                            }
                            print_factor("C", flist->factors[i].factor, OBASE);
						}
					}
				}
                else if (NUM_WITNESSES > 0)
                {
                    if (is_mpz_prp(flist->factors[i].factor, NUM_WITNESSES))
                    {
                        for (j = 0; j < flist->factors[i].count; j++)
                        {
                            mpz_mul(tmp, tmp, flist->factors[i].factor);
                            if (mpz_cmp_ui(flist->factors[i].factor, 100000000) < 0)
                            {
                                if (print_terse) {
                                    tersebuf = append_factor(tersebuf, flist->factors[i].factor);
                                }
                                print_factor("P", flist->factors[i].factor, OBASE);
                            }
                            else
                            {
                                if (print_terse) {
                                    tersebuf = append_factor(tersebuf, flist->factors[i].factor);
                                }
                                print_factor("PRP", flist->factors[i].factor, OBASE);
                            }
                        }
                    }
                    else
                    {
                        for (j = 0; j < flist->factors[i].count; j++)
                        {
                            mpz_mul(tmp, tmp, flist->factors[i].factor);
                            if (print_terse) {
                                tersebuf = append_factor(tersebuf, flist->factors[i].factor);
                            }
                            print_factor("C", flist->factors[i].factor, OBASE);
                        }
                    }
                }
                else
                {
                    for (j = 0; j < flist->factors[i].count; j++)
                    {
                        mpz_mul(tmp, tmp, flist->factors[i].factor);
                        if (print_terse) {
                            tersebuf = append_factor(tersebuf, flist->factors[i].factor);
                        }
                        print_factor("U", flist->factors[i].factor, OBASE);
                    }
                }
			}
		}

        //printf("after processing factors, tersebuf = %s\n", tersebuf);

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
                    {
                        if (OBASE == 8)
                        {
                            gmp_printf("\n***co-factor***\nP%d = %Zo\n",
                                mpz_sizeinbase(tmp2, 8), tmp2);
                        }
                        else if (OBASE == 10)
                        {
                            gmp_printf("\n***co-factor***\nP%d = %Zd\n",
                                gmp_base10(tmp2), tmp2);
                        }
                        else if (OBASE == 16)
                        {
                            gmp_printf("\n***co-factor***\nP%d = %Zx\n",
                                mpz_sizeinbase(tmp2, 16), tmp2);
                        }
                    }
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
                        if (OBASE == 8)
                        {
                            gmp_printf("\n***co-factor***\nC%d = %Zo\n",
                                mpz_sizeinbase(tmp2, 8), tmp2);
                        }
                        else if (OBASE == 10)
                        {
                            gmp_printf("\n***co-factor***\nC%d = %Zd\n",
                                gmp_base10(tmp2), tmp2);
                        }
                        else if (OBASE == 16)
                        {
                            gmp_printf("\n***co-factor***\nC%d = %Zx\n",
                                mpz_sizeinbase(tmp2, 16), tmp2);
                        }
					}
				}
                else if (NUM_WITNESSES > 0)
                {
                    if (is_mpz_prp(tmp2, NUM_WITNESSES))
                    {
                        if (OBASE == 8)
                        {
                            gmp_printf("\n***co-factor***\nPRP%d = %Zo\n",
                                mpz_sizeinbase(tmp2, 8), tmp2);
                        }
                        else if (OBASE == 10)
                        {
                            gmp_printf("\n***co-factor***\nPRP%d = %Zd\n",
                                gmp_base10(tmp2), tmp2);
                        }
                        else if (OBASE == 16)
                        {
                            gmp_printf("\n***co-factor***\nPRP%d = %Zx\n",
                                mpz_sizeinbase(tmp2, 16), tmp2);
                        }
                    }
                    else
                    {
                        if (OBASE == 8)
                        {
                            gmp_printf("\n***co-factor***\nC%d = %Zo\n",
                                mpz_sizeinbase(tmp2, 8), tmp2);
                        }
                        else if (OBASE == 10)
                        {
                            gmp_printf("\n***co-factor***\nC%d = %Zd\n",
                                gmp_base10(tmp2), tmp2);
                        }
                        else if (OBASE == 16)
                        {
                            gmp_printf("\n***co-factor***\nC%d = %Zx\n",
                                mpz_sizeinbase(tmp2, 16), tmp2);
                        }
                    }
                }
                else
                {
                    if (OBASE == 8)
                    {
                        gmp_printf("\n***co-factor***\nU%d = %Zo\n",
                            mpz_sizeinbase(tmp2, 8), tmp2);
                    }
                    else if (OBASE == 10)
                    {
                        gmp_printf("\n***co-factor***\nU%d = %Zd\n",
                            gmp_base10(tmp2), tmp2);
                    }
                    else if (OBASE == 16)
                    {
                        gmp_printf("\n***co-factor***\nU%d = %Zx\n",
                            mpz_sizeinbase(tmp2, 16), tmp2);
                    }
                }
			}			
		}

        //gmp_printf("after processing cofactors, tmp2 = %Zd\n", tmp2);
        //printf("after processing cofactors, tersebuf = %s,maxlen = %d\n", tersebuf,maxlen);

        if (print_terse)
        {
            if (mpz_cmp_ui(tmp2, 1) > 0)
            {
                tersebuf = append_factor(tersebuf, tmp2);
            }

            tersebuf[strlen(tersebuf) - 1] = '\0';
            gmp_printf("\n***factorization:***\n%Zd=%s\n", fobj->input_N, tersebuf);
            free(tersebuf);
        }

		mpz_clear(tmp);
		mpz_clear(tmp2);
	}

	return;
}


