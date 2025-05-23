/*
Copyright (c) 2024, Ben Buhrow
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.
*/


#include "avx_ecm.h"
#include "omp.h"
#include "soe.h"
#include "gmp.h"
#include "calc.h"
#include "factor.h"
#include "yafu_ecm.h"
#include "cmdOptions.h"

int main(int argc, char **argv)
{
    mpz_t gmpn;
	uint64_t b1;
	uint32_t i;
	int threads;    
    str_t input;
    fact_obj_t* fobj;
    uint32_t numcurves;
    uint64_t B1, B2;
    int* numfactors;
    int verbose;
    int save_b1;
    uint32_t* curves_run;
    options_t* options;

    mpz_init(gmpn);

    // process the command line
    processOpts(argc, argv, options);

    sInit(&input);
    calc_init(options->rand_seed);
    toStr(argv[1], &input);
    calc(&input);
    calc_finalize();

    mpz_set_str(gmpn, input.s, 10);

    // a factorization object that gets passed around to any factorization routine
    // called out in the input expression.  if no factorization routine is specified,
    // this is not used.  initialize and pass in all of the options.
    fobj = (fact_obj_t*)malloc(sizeof(fact_obj_t));
    init_factobj(fobj);

    options_to_factobj(fobj, options);
    
    mpz_set(fobj->ecm_obj.gmp_n, gmpn);
    mpz_set(fobj->N, gmpn);

    ecm_loop(fobj);

    // deal with factors
    print_factors(fobj->factors, fobj->N, fobj->VFLAG, fobj->NUM_WITNESSES, fobj->OBASE);

    reset_factobj(fobj);
    free_factobj(fobj);
    free(fobj);
    for (i = 0; i < options->num_tune_info; i++)
    {
        free(options->tune_info[i]);
    }
    free(options->tune_info);
    free(options);
    mpz_clear(gmpn);
    sFree(&input);

	return 0;
}


void options_to_factobj(fact_obj_t* fobj, options_t* options)
{
    // set any parameters of fobj changed by user options
    uint32_t seed1, seed2;
    int i;

    fobj->num_threads = options->threads;
    strcpy(fobj->flogname, options->factorlog);
    strcpy(fobj->factor_json_name, options->jsonlog);

    // initialize global stuff in fobj
    if (options->rand_seed != 0)
    {
        fobj->seed1 = (uint32_t)options->rand_seed & 0xffffffff;
        fobj->seed2 = (uint32_t)(options->rand_seed >> 32);
        fobj->lcg_state = options->rand_seed;
    }

    fobj->NUM_WITNESSES = options->num_prp_witnesses;

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
    strncpy(fobj->pm1_obj.vpm1_work_file, options->vpm1_work_file, 256);
    strncpy(fobj->pm1_obj.resume_file, options->resume_file, 256);

    // initialize stuff for pp1	
    fobj->pp1_obj.B1 = options->B1pp1;
    fobj->pp1_obj.B2 = options->B2pp1;
    fobj->pp1_obj.stg2_is_default = (options->B2pp1 == 0);
    fobj->pp1_obj.pp1_exponent = 0;
    fobj->pp1_obj.pp1_multiplier = 0;
    fobj->pp1_obj.pp1_tune_freq = 0;
    fobj->pp1_obj.lcg_state = hash64(lcg_rand_64(&fobj->lcg_state));
    strncpy(fobj->pp1_obj.vpp1_work_file, options->vpp1_work_file, 256);
    strncpy(fobj->pp1_obj.resume_file, options->resume_file, 256);

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
    fobj->ecm_obj.lcg_state = (uint64_t*)xrealloc(fobj->ecm_obj.lcg_state,
        options->threads * sizeof(uint64_t));
    for (i = 0; i < (int)options->threads; i++)
    {
        fobj->ecm_obj.lcg_state[i] =
            hash64(lcg_rand_64(&fobj->lcg_state));
    }
    fobj->ecm_obj.use_cgbn = options->use_cgbn;
    fobj->ecm_obj.use_gpudev = options->use_gpudev;
    fobj->ecm_obj.use_gpuecm = options->use_gpuecm;
    fobj->ecm_obj.gpucurves = options->gpucurves;

    // unlike ggnfs, ecm does not *require* external binaries.  
    // an empty string indicates the use of the built-in GMP-ECM hooks, while
    // a non-empty string (filled in by the user) will indicate the use of
    // an external binary
    strcpy(fobj->ecm_obj.ecm_path, options->ecm_path);
    fobj->ecm_obj.use_external = 0;
    fobj->ecm_obj.prefer_gmpecm = options->prefer_gmpecm;
    fobj->ecm_obj.prefer_gmpecm_stg2 = options->prefer_gmpecm_stg2;
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
    fobj->qs_obj.gbl_override_ssidx_flag = (options->siqsSSidx > 0);
    fobj->qs_obj.gbl_override_ssidx = options->siqsSSidx;
    fobj->qs_obj.gbl_override_ssalloc_flag = (options->siqsSSalloc > 0);
    fobj->qs_obj.gbl_override_ssalloc = options->siqsSSalloc;
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
    fobj->nfs_obj.minrels = options->minrels;
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
    if (options->sieveQstart > 0)
    {
        if ((options->sieveQstart == 1) && (options->sieveQstop == 1))
        {
            fobj->nfs_obj.startq = fobj->nfs_obj.rangeq = 0;
        }
        fobj->nfs_obj.nfs_phases |= NFS_PHASE_SIEVE;
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
    if (options->polystart > 0)
    {
        if ((options->polystart == 1) && (options->polystop == 1))
        {
            fobj->nfs_obj.polystart = fobj->nfs_obj.polyrange = 0;
        }
        fobj->nfs_obj.nfs_phases |= NFS_PHASE_POLY;
    }

    if (options->nc)
    {
        fobj->nfs_obj.nfs_phases |= NFS_PHASE_FILTER;
        fobj->nfs_obj.nfs_phases |= NFS_PHASE_LA;
        fobj->nfs_obj.nfs_phases |= NFS_PHASE_SQRT;
    }
    if (options->ncr)
    {
        fobj->nfs_obj.nfs_phases |= NFS_PHASE_LA_RESUME;
    }
    if (options->nc1)
    {
        fobj->nfs_obj.nfs_phases |= NFS_PHASE_FILTER;
    }
    if (options->nc2)
    {
        fobj->nfs_obj.nfs_phases |= NFS_PHASE_LA;
    }
    if (options->nc3)
    {
        fobj->nfs_obj.nfs_phases |= NFS_PHASE_SQRT;
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
    fobj->nfs_obj.poly_option = 4;
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
    fobj->nfs_obj.snfs_testsieve_threshold = options->snfs_testsieve_threshold;
    *fobj->nfs_obj.filearg = '\0';

    fobj->nfs_obj.polybatch = options->poly_batch;
    strcpy(fobj->nfs_obj.ggnfs_dir, options->ggnfs_dir);

    fobj->nfs_obj.cadoMsieve = options->cadoMsieve;
    strcpy(fobj->nfs_obj.cado_dir, options->cado_dir);
    strcpy(fobj->nfs_obj.convert_poly_path, options->convert_poly_path);
    fobj->nfs_obj.skip_snfs_check = options->skip_snfscheck;

    // initialize autofactor object
    // whether we want to output certain info to their own files...
    fobj->autofact_obj.want_output_primes = 0;
    fobj->autofact_obj.want_output_factors = 0;
    fobj->autofact_obj.want_output_unfactored = 0;
    fobj->autofact_obj.want_output_expressions = options->want_output_expr;
    if (strlen(options->opfile) > 0)
    {
        fobj->autofact_obj.want_output_primes = 1;
        strcpy(fobj->autofact_obj.op_str, options->opfile);
    }
    if (strlen(options->offile) > 0)
    {
        fobj->autofact_obj.want_output_factors = 1;
        strcpy(fobj->autofact_obj.of_str, options->offile);
    }
    if (strlen(options->oufile) > 0)
    {
        fobj->autofact_obj.want_output_unfactored = 1;
        strcpy(fobj->autofact_obj.ou_str, options->oufile);
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
    fobj->autofact_obj.json_pretty = options->json_pretty;
    fobj->autofact_obj.stopbase = options->stopbase;
    fobj->autofact_obj.stopeq = options->stopeq;
    fobj->autofact_obj.stopge = options->stopge;
    fobj->autofact_obj.stopgt = options->stopgt;
    fobj->autofact_obj.stople = options->stople;
    fobj->autofact_obj.stoplt = options->stoplt;
    fobj->autofact_obj.stopk = options->stopk;
    fobj->autofact_obj.stop_strict = options->strict;
    fobj->autofact_obj.stopprime = options->stopprime;
    fobj->autofact_obj.check_stop_conditions = options->check_stop_conditions;


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
