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
       				   --bbuhrow@gmail.com 11/2/10
----------------------------------------------------------------------*/

#include "yafu_ecm.h"
#include "factor.h"
#include "calc.h"
#include "threadpool.h"
#include <ecm.h>
#include <signal.h>
#include <time.h>

//local function declarations
typedef struct {
    mpz_t gmp_n, gmp_factor;
    mpz_t t, d; // scratch bignums
    ecm_params params;
    uint32_t sigma;
    int stagefound;
    fact_obj_t* fobj;
    int thread_num;
    int curves_run;
    int* total_curves_run;
    char tmp_output[80];
    int factor_found;
    int* ok_to_stop;
    int* curves_in_flight;
    int resume_avx_ecm;

    // timing data for ETA
    struct timeval stop;
    struct timeval start;
    double* total_time;

} ecm_thread_data_t;

void* ecm_do_one_curve(void* ptr);
int print_B1B2(fact_obj_t* fobj, FILE* fid);
void ecmexit(int sig);
void ecm_process_init(fact_obj_t* fobj);
void ecm_process_free(fact_obj_t* fobj);
int ecm_check_input(fact_obj_t* fobj);
int ecm_get_sigma(ecm_thread_data_t* thread_data);
int ecm_deal_with_factor(ecm_thread_data_t* thread_data);
void ecm_stop_worker_thread(ecm_thread_data_t* t, uint32_t is_master_thread);
void ecm_start_worker_thread(ecm_thread_data_t* t, uint32_t is_master_thread);
void ecm_thread_free(ecm_thread_data_t* tdata);
void ecm_thread_init(ecm_thread_data_t* tdata);

int ECM_ABORT;
int TMP_THREADS;
uint64_t ECM_TMP_STG2_MAX;

#if defined(WIN32) || defined(_WIN64)
DWORD WINAPI ecm_worker_thread_main(LPVOID thread_data);
#else
void* ecm_worker_thread_main(void* thread_data);
#endif


void ecm_sync_fcn(void *ptr)
{
    tpool_t *tdata = (tpool_t *)ptr;
    ecm_thread_data_t *udata = (ecm_thread_data_t *)tdata->user_data;
    ecm_thread_data_t *thread_data = &udata[tdata->tindex];
    fact_obj_t *fobj = thread_data->fobj;

    // usually the threads just run one curve at a time, but 
    // nothing in the infrastructure says it can't run more.
    *thread_data->total_curves_run += thread_data->curves_run;
    thread_data->curves_run = 0;

    (*thread_data->curves_in_flight)--;

    // progress report
    if (fobj->VFLAG >= 0)
    {
        printf("ecm: %d/%d curves on C%d, ",
            *thread_data->total_curves_run, fobj->ecm_obj.num_curves,
            (int)gmp_base10(fobj->ecm_obj.gmp_n));

        print_B1B2(fobj, NULL);

        // stop counter for this batch of curves
        // for larger B1s
        if (fobj->ecm_obj.B1 > 48000)
        {
            double curve_time;
            double avg_curve_time;
            double est_time;
            double curves_left;

            curve_time = ytools_difftime(&thread_data->start, &thread_data->stop);
            *thread_data->total_time += curve_time;
            avg_curve_time = *thread_data->total_time / (double)(*thread_data->total_curves_run);
            curves_left = (int)fobj->ecm_obj.num_curves - *thread_data->total_curves_run;
            if (curves_left < 0) curves_left = 0.;
            est_time = (double)(curves_left / fobj->THREADS) * avg_curve_time;

            if (est_time > 3600)
                printf(", ETA: %1.2f hrs ", est_time / 3600);
            else if (est_time > 60)
                printf(", ETA: %1.1f min ", est_time / 60);
            else
                printf(", ETA: %1.0f sec ", est_time);
        }
        printf("\r");
        fflush(stdout);
    }

    if (thread_data->factor_found == 1)
    {
        // since we could be doing many curves in parallel,
        // it's possible more than one curve found this factor at the
        // same time.  We divide out factors as we find them, so
        // just check if this factor still divides n
        mpz_tdiv_qr(thread_data->t, thread_data->d,
            fobj->ecm_obj.gmp_n, thread_data->gmp_factor);
        if (mpz_cmp_ui(thread_data->d, 0) == 0)
        {
            // yes, it does... proceed to record the factor
            // and reduce the active number by the factor.
            mpz_set(fobj->ecm_obj.gmp_n, thread_data->t);
            ecm_deal_with_factor(thread_data);
        }
    }

    // watch for an abort
    if (ECM_ABORT)
    {
        print_factors(fobj->factors, fobj->N, fobj->VFLAG, fobj->NUM_WITNESSES);
        fobj->ecm_obj.exit_cond = ECM_EXIT_ABORT;
        //exit(1);
    }

    return;
}

void ecm_dispatch_fcn(void *ptr)
{
    tpool_t *tdata = (tpool_t *)ptr;
    ecm_thread_data_t *udata = (ecm_thread_data_t *)tdata->user_data;
    ecm_thread_data_t *thread_data = &udata[tdata->tindex];

    // handle found factors
    // we found a factor and might want to stop,
    // but should check all threads output first
    if (thread_data->factor_found == 1)
    {
        if (thread_data->fobj->ecm_obj.bail_on_factor)
        {
            *thread_data->ok_to_stop = 1;
        }
        else if (is_mpz_prp(thread_data->fobj->ecm_obj.gmp_n, 
            thread_data->fobj->NUM_WITNESSES))
        {
            // todo: should check primality in the work function so
            // this prp check doesn't bottleneck dispatching to
            // other threads.  report results to the thread structure
            // and handle here.
            // always bail if the cofactor is prime
            *thread_data->ok_to_stop = 1;
        }
        else
        {
            // set the thread local copy of n.
            // todo: we should not be setting or reading from fobj.gmp_n
            // in the work function - not thread safe.
            mpz_set(thread_data->gmp_n, thread_data->fobj->ecm_obj.gmp_n);
        }
    }

    if (*thread_data->ok_to_stop == 1)
    {
        // found a factor - ok to stop              
        tdata->work_fcn_id = tdata->num_work_fcn;
    }
    else if (thread_data->fobj->ecm_obj.exit_cond == ECM_EXIT_ABORT)
    {
        // user abort - ok to stop  
        tdata->work_fcn_id = tdata->num_work_fcn;
    }
    else if (*thread_data->total_curves_run < 
        (thread_data->fobj->ecm_obj.num_curves - *thread_data->curves_in_flight))
    {
        // more curves to run
        (*thread_data->curves_in_flight)++;
        tdata->work_fcn_id = 0;
    }
    else
    {
        // done running curves
        tdata->work_fcn_id = tdata->num_work_fcn;
    }

    return;
}

void ecm_start_fcn(tpool_t *ptr)
{
    tpool_t *tpool = (tpool_t *)ptr;
    ecm_thread_data_t *udata = (ecm_thread_data_t *)tpool->user_data;
    ecm_thread_data_t *tdata = &udata[tpool->tindex];
    
    //initialize everything for all threads using GMP-ECM
    mpz_init(tdata->gmp_n);
    mpz_init(tdata->gmp_factor);
    ecm_init(tdata->params);
    gmp_randseed_ui(tdata->params->rng,
        lcg_rand_32(&tdata->fobj->ecm_obj.lcg_state[tpool->tindex]));
    mpz_set(tdata->gmp_n, tdata->fobj->ecm_obj.gmp_n);
    tdata->params->method = ECM_ECM;
    tdata->curves_run = 0;
    mpz_init(tdata->d);
    mpz_init(tdata->t);

    return;
}

void ecm_stop_fcn(tpool_t *ptr)
{
    tpool_t *tpool = (tpool_t *)ptr;
    ecm_thread_data_t *udata = (ecm_thread_data_t *)tpool->user_data;
    ecm_thread_data_t *tdata = &udata[tpool->tindex];

    ecm_clear(tdata->params);
    mpz_clear(tdata->gmp_n);
    mpz_clear(tdata->gmp_factor);
    mpz_clear(tdata->d);
    mpz_clear(tdata->t);

    if (tdata->fobj->ecm_obj.use_external)
    {
        // remove temp file specific to this thread num
        sprintf(tdata->tmp_output, "_yafu_ecm_tmp%d.out", tdata->thread_num);
        remove(tdata->tmp_output);
    }

    return;
}

int ecm_loop(fact_obj_t *fobj)
{
	//expects the input in ecm_obj->gmp_n
	FILE *flog = NULL;
	int i,j;
    double total_time = 0;
	int num_batches;    

	//maybe make this an input option: whether or not to stop after
	//finding a factor in the middle of running a requested batch of curves
	int total_curves_run;
	int bail_on_factor = 0;
	int bail = 0;
    int curves_in_flight = 0;
	int input_digits = gmp_base10(fobj->ecm_obj.gmp_n);

    // threading structures
    tpool_t *tpool_data;
    ecm_thread_data_t *thread_data;		//an array of thread data objects

    // if nothing else happens to overwrite it, set the normal exit condition.
    fobj->ecm_obj.exit_cond = ECM_EXIT_NORMAL;

    // check if we need to continue ecm
	if (ecm_check_input(fobj) == 0)
		return 0;
	
    // AVX-ECM won't work without AVX512F and the MINGW build appears to not work either (compiler bug?)
#if !defined(USE_AVX512F) || defined(__MINGW32__)
    if (1)
#else
    if ((fobj->ecm_obj.prefer_gmpecm) || 
        (fobj->ecm_obj.B1 > fobj->ecm_obj.ecm_ext_xover) ||
        (!fobj->HAS_AVX512F))
#endif
    {
        // initialize the flag to watch for interrupts, and set the
        // pointer to the function to call if we see a user interrupt
        ECM_ABORT = 0;
        signal(SIGINT, ecmexit);

        //init ecm process
        ecm_process_init(fobj);

        thread_data = (ecm_thread_data_t*)malloc(fobj->THREADS * sizeof(ecm_thread_data_t));
        for (i = 0; i < fobj->THREADS; i++)
        {
            // several things in the thread structure need to be tied
            // to one master copy for syncronization purposes.
            // there should maybe be a separate field in the threadpool
            // structure for this...
            // e.g. user_data and user_shared_data structures
            thread_data[i].fobj = fobj;
            thread_data[i].thread_num = i;  // todo: remove and replace with tpool->tindex
            thread_data[i].total_curves_run = &total_curves_run;
            thread_data[i].total_time = &total_time;
            thread_data[i].factor_found = 0;
            thread_data[i].ok_to_stop = &bail;
            thread_data[i].curves_in_flight = &curves_in_flight;
            thread_data[i].resume_avx_ecm = 0;
        }

        tpool_data = tpool_setup(fobj->THREADS, &ecm_start_fcn, &ecm_stop_fcn, &ecm_sync_fcn,
            &ecm_dispatch_fcn, thread_data);

        total_curves_run = 0;
        tpool_add_work_fcn(tpool_data, &ecm_do_one_curve);
        tpool_go(tpool_data);

        free(tpool_data);

        if (fobj->VFLAG >= 0)
            printf("\n");

        if (fobj->LOGFLAG && (strcmp(fobj->flogname, "") != 0))
        {
            flog = fopen(fobj->flogname, "a");
            if (flog == NULL)
            {
                printf("fopen error: %s\n", strerror(errno));
                printf("could not open %s for appending\n", fobj->flogname);
                return 0;
            }
            else
            {
                logprint(flog, "Finished %d curves using GMP-ECM method on C%d input, ",
                    total_curves_run, input_digits);

                print_B1B2(fobj, flog);
                fprintf(flog, "\n");

                fclose(flog);
            }
        }

        // this is how we tell factor() to stop running curves at this level
        if (((bail_on_factor == 1) && (bail == 1)) ||
            (mpz_cmp_ui(fobj->ecm_obj.gmp_n, 1) == 0))
        {
            total_curves_run = fobj->ecm_obj.num_curves;
        }

        free(thread_data);
        signal(SIGINT, NULL);
        ecm_process_free(fobj);
    }
    else
    {
        mpz_t F;
        int numfactors = 0;
        uint64_t B2;
        uint32_t curves_run;
        int save_b1 = fobj->ecm_obj.save_b1;
        uint64_t ecm_ext_xover;

        if (fobj->ecm_obj.stg2_is_default)
        {
            B2 = 100 * (uint64_t)fobj->ecm_obj.B1;
        }
        else
        {
            B2 = fobj->ecm_obj.B2;
        }

        if (fobj->ecm_obj.prefer_gmpecm_stg2)
        {
            save_b1 = 2;
            B2 = fobj->ecm_obj.B1;
            ecm_ext_xover = fobj->ecm_obj.ecm_ext_xover;
            fobj->ecm_obj.ecm_ext_xover = 0;
        }

        vec_ecm_main(fobj, fobj->ecm_obj.num_curves, fobj->ecm_obj.B1,
            B2, fobj->THREADS, &numfactors, fobj->VFLAG, save_b1,
            &curves_run);

        if (fobj->ecm_obj.prefer_gmpecm_stg2)
        {
            // now resume each curve with gmp-ecm


            // initialize the flag to watch for interrupts, and set the
            // pointer to the function to call if we see a user interrupt
            ECM_ABORT = 0;
            signal(SIGINT, ecmexit);

            //init ecm process
            ecm_process_init(fobj);

            thread_data = (ecm_thread_data_t*)malloc(fobj->THREADS * sizeof(ecm_thread_data_t));
            for (i = 0; i < fobj->THREADS; i++)
            {
                // several things in the thread structure need to be tied
                // to one master copy for syncronization purposes.
                // there should maybe be a separate field in the threadpool
                // structure for this...
                // e.g. user_data and user_shared_data structures
                thread_data[i].fobj = fobj;
                thread_data[i].thread_num = i;  // todo: remove and replace with tpool->tindex
                thread_data[i].total_curves_run = &total_curves_run;
                thread_data[i].total_time = &total_time;
                thread_data[i].factor_found = 0;
                thread_data[i].ok_to_stop = &bail;
                thread_data[i].curves_in_flight = &curves_in_flight;
                thread_data[i].resume_avx_ecm = 1;
            }

            tpool_data = tpool_setup(fobj->THREADS, &ecm_start_fcn, &ecm_stop_fcn, &ecm_sync_fcn,
                &ecm_dispatch_fcn, thread_data);

            total_curves_run = 0;
            tpool_add_work_fcn(tpool_data, &ecm_do_one_curve);
            tpool_go(tpool_data);

            free(tpool_data);

            if (fobj->VFLAG >= 0)
                printf("\n");

            if (fobj->LOGFLAG && (strcmp(fobj->flogname, "") != 0))
            {
                flog = fopen(fobj->flogname, "a");
                if (flog == NULL)
                {
                    printf("fopen error: %s\n", strerror(errno));
                    printf("could not open %s for appending\n", fobj->flogname);
                    return 0;
                }
                else
                {
                    logprint(flog, "Finished %d curves using GMP-ECM method on C%d input, ",
                        total_curves_run, input_digits);

                    print_B1B2(fobj, flog);
                    fprintf(flog, "\n");

                    fclose(flog);
                }
            }

            free(thread_data);
            signal(SIGINT, NULL);
            ecm_process_free(fobj);

            fobj->ecm_obj.ecm_ext_xover = ecm_ext_xover;
        }

        // this is how we tell factor() to stop running curves at this level
        if (((bail_on_factor == 1) && (bail == 1)) ||
            (mpz_cmp_ui(fobj->ecm_obj.gmp_n, 1) == 0))
        {
            total_curves_run = MAX(fobj->ecm_obj.num_curves, curves_run);
        }
        else
        {
            total_curves_run = curves_run;
        }

    }

	return total_curves_run;
}

int ecm_deal_with_factor(ecm_thread_data_t *thread_data)
{
	fact_obj_t *fobj = thread_data->fobj;
	int curves_run = thread_data->curves_run;
	int thread_num = thread_data->thread_num;

	if (is_mpz_prp(thread_data->gmp_factor, fobj->NUM_WITNESSES))
	{
		add_to_factor_list(fobj->factors, thread_data->gmp_factor, 
            fobj->VFLAG, fobj->NUM_WITNESSES);

        if (fobj->VFLAG > 0)
        {
            gmp_printf("\necm: found prp%d factor = %Zd\n",
                gmp_base10(thread_data->gmp_factor),
                thread_data->gmp_factor);
        }

        char* s = mpz_get_str(NULL, 10, thread_data->gmp_factor);
		logprint_oc(fobj->flogname, "a", "prp%d = %s (curve %d stg%d B1=%u sigma=%u thread=%d)\n",
			gmp_base10(thread_data->gmp_factor), s,
			*thread_data->total_curves_run + 1, thread_data->stagefound,
			fobj->ecm_obj.B1, thread_data->sigma, thread_num);
        free(s);
	}
	else
	{
		add_to_factor_list(fobj->factors, thread_data->gmp_factor, 
            fobj->VFLAG, fobj->NUM_WITNESSES);
		
        if (fobj->VFLAG > 0)
        {
            gmp_printf("\necm: found c%d factor = %Zd\n",
                gmp_base10(thread_data->gmp_factor),
                thread_data->gmp_factor);
        }

        char* s = mpz_get_str(NULL, 10, thread_data->gmp_factor);
		logprint_oc(fobj->flogname, "a", "c%d = %s (curve %d stg%d B1=%u sigma=%u thread=%d)\n",
			gmp_base10(thread_data->gmp_factor), s,
            *thread_data->total_curves_run + 1, thread_data->stagefound,
			fobj->ecm_obj.B1, thread_data->sigma, thread_num);
	}

	return 1;
}

int ecm_get_sigma(ecm_thread_data_t *thread_data)
{
	fact_obj_t *fobj = thread_data->fobj;
    ecm_obj_t* ecmobj = &fobj->ecm_obj;
	mpz_t tmp;

	mpz_init(tmp);
	if (fobj->ecm_obj.sigma != 0)
	{
		//if (thread_data->curves_run == 0 &&  fobj->ecm_obj.num_curves > 1)
	//		printf("WARNING: work will be duplicated with sigma fixed and numcurves > 1\n");
		thread_data->sigma = fobj->ecm_obj.sigma;
	}
	else //if (get_uvar("sigma",tmp))
	{
		thread_data->sigma = lcg_rand_32_range(6, 0xffffffff, 
            &ecmobj->lcg_state[thread_data->thread_num]);
	}
	//else
	//{
	//	if (thread_data->curves_run == 0 &&  fobj->ecm_obj.num_curves > 1)
	//		printf("WARNING: work will be duplicated with sigma fixed and numcurves > 1\n");
	//	thread_data->sigma = (uint32_t)mpz_get_ui(tmp);
	//}

	mpz_clear(tmp);
	return 0;
}

int ecm_check_input(fact_obj_t *fobj)
{
	//check for trivial cases
	if (mpz_cmp_ui(fobj->ecm_obj.gmp_n, 0) == 0)
	{
		logprint_oc(fobj->flogname, "a","Trivial input == 0 in ECM\n");
		return 0;
	}

	if (mpz_cmp_ui(fobj->ecm_obj.gmp_n, 1) == 0)
	{
		logprint_oc(fobj->flogname, "a","Trivial input == 1 in ECM\n");
		return 0;
	}

	if (mpz_tdiv_ui(fobj->ecm_obj.gmp_n, 3) == 0)
	{
		mpz_t tmp;
		mpz_init(tmp);
		mpz_set_ui(tmp, 3);
		mpz_tdiv_q_ui(fobj->ecm_obj.gmp_n, fobj->ecm_obj.gmp_n, 3);
		add_to_factor_list(fobj->factors, tmp,
            fobj->VFLAG, fobj->NUM_WITNESSES);
		logprint_oc(fobj->flogname, "a","Trivial factor of 3 found in ECM\n");
		mpz_clear(tmp);
		return 0;
	}

	if (mpz_tdiv_ui(fobj->ecm_obj.gmp_n, 2) == 0)
	{
		mpz_t tmp;
		mpz_init(tmp);
		mpz_set_ui(tmp, 2);
		mpz_tdiv_q_ui(fobj->ecm_obj.gmp_n, fobj->ecm_obj.gmp_n, 2);
		add_to_factor_list(fobj->factors, tmp,
            fobj->VFLAG, fobj->NUM_WITNESSES);
		logprint_oc(fobj->flogname, "a","Trivial factor of 2 found in ECM\n");
		mpz_clear(tmp);
		return 0;
	}

	if (is_mpz_prp(fobj->ecm_obj.gmp_n, fobj->NUM_WITNESSES))
	{
		//maybe have an input flag to optionally not perform
		//PRP testing (useful for really big inputs)
		add_to_factor_list(fobj->factors, fobj->ecm_obj.gmp_n,
            fobj->VFLAG, fobj->NUM_WITNESSES);
        char* s = mpz_get_str(NULL, 10, fobj->ecm_obj.gmp_n);
		logprint_oc(fobj->flogname, "a","prp%d = %s\n", gmp_base10(fobj->ecm_obj.gmp_n), s);		
		mpz_set_ui(fobj->ecm_obj.gmp_n, 1);
        free(s);
		return 0;
	}

	return 1;
}

void ecm_process_init(fact_obj_t *fobj)
{
	// initialize things which all threads will need when using
	// GMP-ECM
	TMP_THREADS = fobj->THREADS;
	ECM_TMP_STG2_MAX = fobj->ecm_obj.B2;

    if (strcmp(fobj->ecm_obj.ecm_path, "") != 0)
    {
        // check that the file actually exists... print a warning if not and
        // resort to internal ecm.
        if (NULL == fopen(fobj->ecm_obj.ecm_path, "rb"))
        {
            printf("ecm: ECM executable does not exist at %s\n"
                "ecm: using internal single threaded ECM...\n",
                fobj->ecm_obj.ecm_path);
            fobj->ecm_obj.use_external = 0;
        }
        else
        {
            fobj->ecm_obj.use_external = 1;
        }
    }
	
	if (fobj->THREADS > 1)
	{
        float ver;

#ifdef ECM_VERSION
        sscanf(ECM_VERSION, "%f", &ver);
#else
        ver = 6;
#endif

		if (!fobj->ecm_obj.use_external)
		{                       
            if (ver > 7.0)
            {
                // ok to use multiple threads in ecm 7+

            }
            else
            {
                if (fobj->VFLAG >= 2)
                {
                    printf("GMP-ECM does not support multiple threads... running single threaded\n");
                }
                fobj->THREADS = 1;
            }
		}
		else
		{
			if ((fobj->ecm_obj.B1 < fobj->ecm_obj.ecm_ext_xover) || (fobj->ecm_obj.num_curves == 1))
			{
                if (ver > 7.0)
                {
                    // ok to use multiple threads in ecm 7+
                }
                else
                {
                    fobj->THREADS = 1;                    
                }
                fobj->ecm_obj.use_external = 0;
			}
		}
	}
	else
	{
        if (fobj->ecm_obj.B1 < fobj->ecm_obj.ecm_ext_xover)
        {
            fobj->ecm_obj.use_external = 0;
        }
	}

	return;
}

void ecm_process_free(fact_obj_t *fobj)
{
	fobj->THREADS = TMP_THREADS;
	fobj->ecm_obj.B2 = ECM_TMP_STG2_MAX;
	return;
}

void ecm_thread_init(ecm_thread_data_t *tdata)
{
    // initialize everything for all threads using GMP-ECM
    mpz_init(tdata->gmp_n);
    mpz_init(tdata->gmp_factor);
    ecm_init(tdata->params);
    gmp_randseed_ui(tdata->params->rng, 
        lcg_rand_32(&tdata->fobj->ecm_obj.lcg_state[tdata->thread_num]));
    mpz_set(tdata->gmp_n, tdata->fobj->ecm_obj.gmp_n);
    tdata->params->method = ECM_ECM;
    tdata->curves_run = 0;
    mpz_init(tdata->d);
    mpz_init(tdata->t);

    return;
}

void ecm_thread_free(ecm_thread_data_t *tdata)
{
    ecm_clear(tdata->params);
    mpz_clear(tdata->gmp_n);
    mpz_clear(tdata->gmp_factor);
    mpz_clear(tdata->d);
    mpz_clear(tdata->t);

    if (tdata->fobj->ecm_obj.use_external)
    {
        // remove temp file specific to this thread num
        sprintf(tdata->tmp_output, "_yafu_ecm_tmp%d.out", tdata->thread_num);
        remove(tdata->tmp_output);
    }

    return;
}

void *ecm_do_one_curve(void *ptr)
{
	// unpack the data structure and stuff inside it
    tpool_t *tpool = (tpool_t *)ptr;
    ecm_thread_data_t *udata = (ecm_thread_data_t *)tpool->user_data;
    ecm_thread_data_t *thread_data = &udata[tpool->tindex];
	fact_obj_t *fobj = thread_data->fobj;

    ecm_get_sigma(thread_data);

    // start a counter for this batch of curves
    // for larger B1s
    if (fobj->ecm_obj.B1 > 48000)
    {
        gettimeofday(&thread_data->start, NULL);
    }

	if (!fobj->ecm_obj.use_external)
	{
		int status;

		thread_data->params->B1done = 1.0 + floor (1 * 128.) / 134217728.;
		if (fobj->VFLAG >= 3)
			thread_data->params->verbose = fobj->VFLAG - 2;		
		mpz_set_ui(thread_data->params->x, (unsigned long)0);
		mpz_set_ui(thread_data->params->sigma, thread_data->sigma);

		if (fobj->ecm_obj.stg2_is_default == 0)
		{
			//not default, tell gmp-ecm to use the requested B2
			//printf("using requested B2 value\n");
			uint64_2gmp(fobj->ecm_obj.B2, thread_data->params->B2);
		}

		status = ecm_factor(thread_data->gmp_factor, thread_data->gmp_n,
				fobj->ecm_obj.B1, thread_data->params);		

		//printf ("used B2: ");
		//mpz_out_str (stdout, 10, thread_data->params->B2);
		//printf ("\n");

		//NOTE: this required a modification to the GMP-ECM source code in ecm.c
		//in order to get the automatically computed B2 value out of the
		//library
		//gmp2mp(thread_data->params->B2,&thread_data->factor);
		//ECM_STG2_MAX = z264(&thread_data->factor);

		//the return value is the stage the factor was found in, if no error
		thread_data->stagefound = status;
	}
	else
	{
		char *cmd;
		FILE *fid;
		char line[1024];
		char *ptr;
		char *tmpstr = NULL;
		int retcode;

		// let mpz figure out and allocate the string
		tmpstr = mpz_get_str(tmpstr, 10, thread_data->gmp_n);

		// allocate the appropriately sized command string
		cmd = (char *)malloc((strlen(tmpstr) + strlen(fobj->ecm_obj.ecm_path) + 256) 
			* sizeof(char));

		// external executable was specified
		sprintf(thread_data->tmp_output, "_yafu_ecm_tmp%d.out", thread_data->thread_num);

        // todo: add command line input of arbitrary argument string to append to this command
		// build system command
        if (thread_data->resume_avx_ecm)
        {
            fid = fopen("avx_ecm_resume.txt", "r");
            if (fid != NULL)
            {
                fclose(fid);
                sprintf(cmd, "echo %s | %s -resume `pwd`/avx_ecm_resume.txt %u | tee %s\n",
                    tmpstr, fobj->ecm_obj.ecm_path, fobj->ecm_obj.B1,
                    thread_data->tmp_output);
            }
            else
            {
                thread_data->curves_run++;
                free(tmpstr);
                free(cmd);
                return;
            }
        }
        else
        {
            sprintf(cmd, "echo %s | %s -sigma %u %u > %s\n",
                tmpstr, fobj->ecm_obj.ecm_path, thread_data->sigma, fobj->ecm_obj.B1,
                thread_data->tmp_output);
        }

		// run system command
        //printf("sys: %s\n", cmd);
		retcode = system(cmd);

		free(tmpstr);
		free(cmd);

		// this is what I observed ecm returning on ctrl-c.  hopefully it is portable.
		if (retcode == 33280)
			ECM_ABORT = 1;
		
		// parse output file
		fid = fopen(thread_data->tmp_output, "r");
		while ((fid != NULL) && (!feof(fid)))
		{
			char fact[1024];

			fgets(line, 1024, fid);
			if (line == NULL)
				break;

			ptr = strstr(line, "**********");
			if (ptr == NULL)
				continue;

			// found a factor.  search for the :
			ptr = strstr(line, ":");
			if (ptr == NULL)
				continue;

			// the character prior to this is the stage, and the rest of the line
			// after it is the factor
			sscanf(ptr-2,"%d",&thread_data->stagefound);

			strcpy(fact, ptr+1);

			mpz_set_str(thread_data->gmp_factor, fact, 10);
			//str2hexz(fact, &thread_data->factor);

			break;
		}
		fclose(fid);

        if (thread_data->resume_avx_ecm)
        {
            time_t rawtime;
            struct tm* info;
            time(&rawtime);
            info = localtime(&rawtime);

            char destfile[80];

            sprintf(destfile, "avx_ecm_resume_%d_%d_%d.txt", 
                info->tm_mon, info->tm_mday, info->tm_year + 1900);
            // archive the resume file
            rename("avx_ecm_resume.txt", destfile);
        }
	}

    if (fobj->ecm_obj.B1 > 48000)
    {
        gettimeofday(&thread_data->stop, NULL);
    }

    //look at the result of each curve and see if we're done
    if ((mpz_cmp_ui(thread_data->gmp_factor, 1) > 0)
        && (mpz_cmp(thread_data->gmp_factor, fobj->ecm_obj.gmp_n) < 0))
    {
        // non-trivial factor found.
        // flag it to be dealt with by the master process.
        thread_data->factor_found = 1;
    }

    thread_data->curves_run++;

	return 0;
}

// function definitions
void ecmexit(int sig)
{
	printf("\nAborting...\n");
	ECM_ABORT = 1;
	return;
}

int print_B1B2(fact_obj_t *fobj, FILE *fid)
{
	int i;
	char suffix;
	char stg1str[20];
	char stg2str[20];

	if (fobj->ecm_obj.B1 % 1000000000 == 0)
	{
		suffix = 'B';
		sprintf(stg1str,"%u%c",fobj->ecm_obj.B1 / 1000000000, suffix);
	}
	else if (fobj->ecm_obj.B1 % 1000000 == 0)
	{
		suffix = 'M';
		sprintf(stg1str,"%u%c",fobj->ecm_obj.B1 / 1000000, suffix);
	}
	else if (fobj->ecm_obj.B1 % 1000 == 0)
	{
		suffix = 'k';
		sprintf(stg1str,"%u%c",fobj->ecm_obj.B1 / 1000, suffix);
	}
	else
	{
		sprintf(stg1str,"%u",fobj->ecm_obj.B1);
	}

	if (fobj->ecm_obj.stg2_is_default == 0)
	{
		if (fobj->ecm_obj.B2 % 1000000000 == 0)
		{
			suffix = 'B';
			sprintf(stg2str,"%" PRIu64 "%c",fobj->ecm_obj.B2 / 1000000000, suffix);
		}
		else if (fobj->ecm_obj.B2 % 1000000 == 0)
		{
			suffix = 'M';
			sprintf(stg2str,"%" PRIu64 "%c",fobj->ecm_obj.B2 / 1000000, suffix);
		}
		else if (fobj->ecm_obj.B2 % 1000 == 0)
		{
			suffix = 'k';
			sprintf(stg2str,"%" PRIu64 "%c",fobj->ecm_obj.B2 / 1000, suffix);
		}
		else
		{
			sprintf(stg2str,"%" PRIu64 "",fobj->ecm_obj.B2);
		}
	}
	else
		sprintf(stg2str, "gmp-ecm default");

	if (fid == NULL)
		i = printf("B1=%s, B2=%s",stg1str,stg2str);
	else
		i = fprintf(fid,"B1=%s, B2=%s",stg1str,stg2str);

	return i;
}



