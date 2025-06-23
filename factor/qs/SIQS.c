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

#include "factor.h"
#include "qs.h"
#include "qs_impl.h"
#include "ytools.h"
#include "gmp_xface.h"
#include "threadpool.h"
#include "cofactorize.h"
#include <immintrin.h>

#ifdef USE_BATCH_FACTOR
#include "batch_factor.h"
#include "prime_sieve.h"
#endif

#include "soe.h"

#if defined( USE_SS_SEARCH ) && defined(USE_POLY_BUCKET_SS)
#define SS_TIMING
#endif

//#define VDEBUG
// opt debug will print out some info relevant to the process used to 
// optimize the tf_small_cutoff value.  This value is varied slightly during the
// first few polynomials and the value which maximizes relation discover rate
// is chosen.
// #define OPT_DEBUG

// define the fat-binary function pointers
void (*testRoots_ptr)(static_conf_t*, dynamic_conf_t*);
void (*nextRoots_ptr)(static_conf_t*, dynamic_conf_t*);
void (*firstRoots_ptr)(static_conf_t*, dynamic_conf_t*);
void (*resieve_med_ptr)(uint8_t, uint32_t, uint32_t,
    static_conf_t*, dynamic_conf_t*);
void (*tdiv_med_ptr)(uint8_t, uint32_t, uint32_t,
    static_conf_t*, dynamic_conf_t*);
void (*tdiv_LP_ptr)(uint32_t, uint8_t, uint32_t,
    static_conf_t* , dynamic_conf_t* );
int (*scan_ptr)(uint32_t, uint8_t, static_conf_t*, dynamic_conf_t*);
void (*lp_sieveblock_ptr)(uint8_t* , uint32_t , uint32_t ,
    lp_bucket* , int , dynamic_conf_t* );
void (*med_sieve_ptr)(uint8_t*, sieve_fb_compressed*, fb_list*, uint32_t, uint8_t);

typedef struct
{
    thread_sievedata_t *thread_data;

    // adaptive tf_small_cutoff variables needed across syncronizations
    int orig_value;
    int num_meas;
    double results[3];
    int averaged_polys;
    double rels_per_sec_avg;
#ifdef OPT_DEBUG
    FILE *optfile;
#endif

    // other stuff that we need a global copy of
    uint32_t num_needed;
    uint32_t num_found;
    int updatecode;
    uint32_t avg_rels_per_acoeff;

    // semaphores for running batch_run_override
    uint32_t ack_batch_override[256];
    uint32_t done_batch_override[256];

} siqs_userdata_t;


typedef struct
{
    uint32_t* hitloc;
    uint32_t* prime;
    uint32_t numhits;
    uint32_t hitalloc;
} poly_bucket_t;


uint64_t* siqs_primes;
uint64_t siqs_nump;
uint64_t siqs_minp;
uint64_t siqs_maxp;
static int siqs_primes_initialized = 0;
static int SIQS_ABORT;

void siqs_start(void *vptr)
{
    tpool_t *tdata = (tpool_t *)vptr;
    siqs_userdata_t *udata = tdata->user_data;
    static_conf_t *static_conf = udata->thread_data[0].sconf;
    fact_obj_t *fobj = static_conf->obj;
    int i;

#ifdef OPT_DEBUG
    udata->optfile = fopen("optfile.csv", "a");
    fprintf(udata->optfile, "Optimization Debug File\n");
    fprintf(udata->optfile, "Detected cpu %d, with L1 = %d bytes, L2 = %d bytes\n",
        ytools_get_cpu_type(), L1CACHE, L2CACHE);
    gmp_fprintf(udata->optfile, "Starting SIQS on c%d: %s\n\n",
        fobj->digits, mpz_conv2str(&gstr1.s, 10, fobj->gmp_n));
    fprintf(udata->optfile, "Meas #,Poly A #, Avg Rels/Poly/Sec, small_tf_cutoff\n");
	printf("Meas #,Poly A #, Avg Rels/Poly/Sec, small_tf_cutoff\n");
#endif


    udata->results[0] = 0.0;
    udata->results[1] = 0.0;
    udata->results[2] = 0.0;

    udata->averaged_polys = 0;
    udata->rels_per_sec_avg = 0.0;
    udata->updatecode = 0;

    // get best parameters, multiplier, and factor base for the job
    // initialize and fill out the static part of the job data structure
    siqs_static_init(udata->thread_data[0].sconf, 0);

    // allocate structures for use in sieving with threads
    for (i = 0; i < fobj->THREADS; i++)
    {
        siqs_dynamic_init(udata->thread_data[i].dconf, static_conf);
    }

    if (fobj->qs_obj.gbl_override_small_cutoff_flag)
    {
        printf("overriding small TF cutoff of %u to %d\n",
            static_conf->tf_small_cutoff, fobj->qs_obj.gbl_override_small_cutoff);
        fobj->qs_obj.no_small_cutoff_opt = 1;
        static_conf->tf_small_cutoff = fobj->qs_obj.gbl_override_small_cutoff;
    }

    // check if a savefile exists for this number, and if so load the data
    // into the master data structure
    siqs_check_restart(udata->thread_data[0].dconf, static_conf);
    print_siqs_splash(udata->thread_data[0].dconf, static_conf);

    // start the process
	udata->num_needed = static_conf->factor_base->B + static_conf->num_extra_relations;
    static_conf->num_needed = udata->num_needed;
    udata->num_found = static_conf->num_r;
    static_conf->total_poly_a = 0;
    udata->num_meas = 0;
    udata->orig_value = static_conf->tf_small_cutoff;

    return;
}

void siqs_sync(void *vptr)
{
    tpool_t *tdata = (tpool_t *)vptr;
    siqs_userdata_t *udata = tdata->user_data;
    thread_sievedata_t *t = udata->thread_data;
    static_conf_t *static_conf = udata->thread_data[0].sconf;
    fact_obj_t *fobj = static_conf->obj;
    
    int tid = tdata->tindex;
    
    // adaptive tf_small_cutoff variables
    int num_avg = 10;    
    int orig_value = udata->orig_value;
    int j;

	if (static_conf->bits > 300)
		num_avg = 5;
	else if (static_conf->bits > 320)
		num_avg = 2;
    
    // Check whether the thread has any results to collect and merge them if so.
    if (t[tid].dconf->buffered_rels)
    {
        udata->num_found = siqs_merge_data(t[tid].dconf, static_conf);        
        static_conf->num_found = udata->num_found;

        if (fobj->qs_obj.no_small_cutoff_opt == 0)
        {
            int	poly_start_num = 0;

            if (udata->num_meas < 3)
            {
                if (udata->averaged_polys >= num_avg)
                {
                    poly_start_num = static_conf->total_poly_a;
                    udata->results[udata->num_meas] = udata->rels_per_sec_avg /
                        (double)udata->averaged_polys;

#ifdef OPT_DEBUG
                    fprintf(udata->optfile, "%d,%d,%f,%d\n", udata->num_meas, 
                        static_conf->total_poly_a,
                        udata->results[udata->num_meas], 
                        static_conf->tf_small_cutoff);
					printf("%d,%d,%f,%d\n", udata->num_meas,
						static_conf->total_poly_a,
						udata->results[udata->num_meas],
						static_conf->tf_small_cutoff);
#endif

                    udata->rels_per_sec_avg = 0.0;
                    udata->averaged_polys = 0;

                    if (udata->num_meas == 0)
                    {
                        static_conf->tf_small_cutoff = orig_value + 5;
                    }
                    else if (udata->num_meas == 1)
                    {
                        static_conf->tf_small_cutoff = orig_value - 5;
                    }
                    else
                    {
                        //we've got our three measurements, make a decision.
                        //the experimental results need to be convincingly better
                        //in order to switch to a different value.  2% for numbers
                        //with one LP, 5% for DLP.
                        if (static_conf->use_dlp)
                            udata->results[0] *= 1.05;
                        else
                            udata->results[0] *= 1.02;

                        if (udata->results[0] > udata->results[1])
                        {
                            if (udata->results[0] > udata->results[2])
                                static_conf->tf_small_cutoff = orig_value;
                            else
                                static_conf->tf_small_cutoff = orig_value - 5;
                        }
                        else
                        {
                            if (udata->results[1] > udata->results[2])
                                static_conf->tf_small_cutoff = orig_value + 5;
                            else
                                static_conf->tf_small_cutoff = orig_value - 5;
                        }

                        //printf("\nfinal value for tf_small_cutoff = %u\n",
                        //    static_conf->tf_small_cutoff);

#ifdef OPT_DEBUG
                        fprintf(udata->optfile, "final value = %d\n", static_conf->tf_small_cutoff);
                        fprintf(udata->optfile, "\n\n");
						printf("final value = %d\n", static_conf->tf_small_cutoff);
						printf("\n\n");
                        fclose(udata->optfile);
#endif
                    }

                    udata->num_meas++;
                }
                else
                {
                    udata->rels_per_sec_avg += t[tid].dconf->rels_per_sec;
                    udata->averaged_polys++;
                }
            }
        }

        // free sieving structure
        for (j = 0; j < t[tid].dconf->buffered_rels; j++)
            t[tid].dconf->relation_buf[j].num_factors = 0;

        t[tid].dconf->num = 0;
        t[tid].dconf->tot_poly = 0;
        t[tid].dconf->buffered_rels = 0;
        t[tid].dconf->attempted_squfof = 0;
        t[tid].dconf->failed_squfof = 0;
        t[tid].dconf->dlp_outside_range = 0;
		t[tid].dconf->num_slp = 0;
		t[tid].dconf->num_full = 0;
		t[tid].dconf->dlp_prp = 0;
        t[tid].dconf->dlp_useful = 0;
		t[tid].dconf->attempted_cosiqs = 0;
		t[tid].dconf->failed_cosiqs = 0;
		t[tid].dconf->tlp_outside_range = 0;
		t[tid].dconf->tlp_prp = 0;
		t[tid].dconf->tlp_useful = 0;
        t[tid].dconf->total_blocks = 0;
        t[tid].dconf->total_reports = 0;
        t[tid].dconf->total_surviving_reports = 0;
        t[tid].dconf->lp_scan_failures = 0;
        t[tid].dconf->num_64bit_residue = 0;
        t[tid].dconf->num_scatter_opp = 0;
        t[tid].dconf->num_scatter = 0;
    }

    return;
}

void siqs_dispatch(void *vptr)
{
    tpool_t *tdata = (tpool_t *)vptr;
    siqs_userdata_t *udata = tdata->user_data;
    thread_sievedata_t *t = udata->thread_data;
    static_conf_t *static_conf = t[0].sconf;
    fact_obj_t *fobj = static_conf->obj;
    int tid = tdata->tindex;    
    double avg_rels_per_acoeff;
    double in_flight_rels;

    if (0) //static_conf->use_dlp == 2)
    {
        uint32_t total_rels = static_conf->num_full + static_conf->num_slp +
            static_conf->dlp_useful + static_conf->tlp_useful;

        if (total_rels > static_conf->check_total)
        {
            // this means we will be running a TLP filtering.  Force 
            // concurrently running threads to perform batch filtering so that
            // we run it with our full data set.
            int i;
            for (i = 0; i < fobj->THREADS; i++);
            {
                t[tid].dconf->batch_run_override = 1;
            }

            // wait a second 

            // but this won't work because even if we can make them run batch
            // factoring, the results only go into a buffer.  We need them
            // in the results file for them to do any good and that would
            // go against the spirit of the threadpool to force synchronous 
            // dumps to the file.
        }
    }

    // check whether to continue or not, and update the screen
    udata->updatecode = update_check(static_conf);

    if ((udata->num_found > 0) && ((int)static_conf->total_poly_a > 0))
    {
        avg_rels_per_acoeff = (double)udata->num_found / (double)(static_conf->total_poly_a + 1);

        if ((static_conf->is_restart == 1) || (static_conf->use_dlp >= 2))
        {

            // we'd have to separately track the relations found since the
            // restart for this to work right.  that impacts many things,
            // so for now just don't worry about overshoot due to multi-threading.
            in_flight_rels = 0;
        }
        else
        {

            in_flight_rels = (tdata->num_threads - 1) * 0.725 * avg_rels_per_acoeff;
            // test: exaggerate this to force more bad matrices
            //in_flight_rels = (tdata->num_threads - 1) * 1.25 * avg_rels_per_acoeff;

            if (static_conf->do_batch)
            {
                // get an estimate of batch conversion rate across all threads
                int i;
                double ratio = 0;
                uint32_t rels = 0;
                for (i = 0; i < fobj->THREADS; i++)
                {
                    ratio += t[i].dconf->rb.conversion_ratio;
                    rels += t[i].dconf->rb.num_relations;
                }

                // now an estimate of the relations in-flight inside the 
                // batch queues (with a fudge factor aiming low).
                in_flight_rels += ratio * 1.1 * (double)rels;
            }
        }
        //printf("found %u rels so far in %d 'A' polys, avg rels_per_poly = %f, rels in flight = %f\n",
        //    udata->num_found, static_conf->total_poly_a + 1, avg_rels_per_acoeff, in_flight_rels);
    }
    else
    {
        in_flight_rels = 0;
    }

    // if we have enough relations, or if there was a break signal, stop dispatching
    // any more threads
    if ((udata->updatecode == 0) && 
        ((udata->num_found + in_flight_rels) < udata->num_needed))
    {
        relation_batch_t *src_rb;
        // generate a new poly A value for the thread we pulled out of the queue
        // using its dconf.  this is done by the master thread because it also 
        // stores the coefficients in a master list
        new_poly_a(static_conf, t[tid].dconf);
        tdata->work_fcn_id = 0;

        if (static_conf->do_batch)
        {
            src_rb = &static_conf->rb[static_conf->batch_buffer_id - 1];

            // if this thread is signalling that it's done processing 
            // a batch of relations, then prepare for the next batch.
            if ((static_conf->use_dlp >= 2) &&
                (t[tid].dconf->batch_run_override == -1))
            {
                if (fobj->VFLAG > 0)
                {
                    printf("thread %d done processing batch %d\n",
                        tid, t[tid].dconf->batch_run_override);
                }
                static_conf->batch_run_override = 0;
                t[tid].dconf->batch_run_override = 0;
                static_conf->num_active_rb--;
            }

            // if we are batching relations and we have enough to process
            // and we are not already processing somewhere, then
            // instruct this thread to start processing batched relations.
            if ((static_conf->use_dlp >= 2) &&
                (src_rb->num_relations > src_rb->target_relations) &&
                (static_conf->num_active_rb < static_conf->num_alloc_rb) &&
                (fobj->THREADS > 1))
                //(static_conf->batch_run_override == 0)) 
            {
                int i;
                int count = 0;
                int first_free = -1;
                uint8_t* rb_active;

                rb_active = (uint8_t*)xmalloc(fobj->THREADS * sizeof(uint8_t));
                memset(rb_active, 0, fobj->THREADS * sizeof(uint8_t));

                // tell this thread to process the batched relations
                // and switch to another batch buffer for incoming batched
                // relations from other threads.
                static_conf->batch_run_override = 1;

                if (fobj->VFLAG > 1)
                {
                    printf("signalling thread %d to process batch %d\n",
                        tid, static_conf->batch_buffer_id);
                }

                if (mpz_cmp_ui(static_conf->rb[static_conf->batch_buffer_id].prime_product, 0) == 0)
                {
                    if (fobj->VFLAG > 1)
                    {
                        printf("copying prime_product to batch buffer %d, mem use +%u bytes\n",
                            static_conf->batch_buffer_id,
                            mpz_sizeinbase(static_conf->rb[0].prime_product, 2) / 8);
                    }
                    mpz_set(static_conf->rb[static_conf->batch_buffer_id].prime_product,
                        static_conf->rb[0].prime_product);
                }

                t[tid].dconf->batch_run_override = static_conf->batch_buffer_id;

                for (i = 0; i < fobj->THREADS; i++)
                {
                    if (t[i].dconf->batch_run_override > 0)
                    {
                        rb_active[t[i].dconf->batch_run_override] = 1;
                        count++;
                    }
                }
                static_conf->num_active_rb = count;
                static_conf->max_active_rb = MAX(count, static_conf->max_active_rb);

                for (i = 1; i < fobj->THREADS; i++)
                {
                    if (rb_active[i] == 0)
                    {
                        first_free = i;
                        break;
                    }
                }

                if (first_free < 0)
                {
                    printf("error no free buffers for batch processing!\n");
                    first_free = 0;
                }

                if (fobj->VFLAG > 0)
                {
                    printf("there are now %d active batch processing threads\n", count);
                }

                static_conf->batch_buffer_id = first_free;

                static_conf->rb[static_conf->batch_buffer_id].target_relations =
                    static_conf->obj->qs_obj.gbl_btarget;

                //static_conf->batch_buffer_id++;
                //if (static_conf->batch_buffer_id == static_conf->num_alloc_rb)
                //{
                //    static_conf->batch_buffer_id = 1;
                //}

                if (fobj->VFLAG > 0)
                {
                    printf("now gathering relations into batch buffer %d\n",
                        static_conf->batch_buffer_id);
                }

                free(rb_active);
            }
            else
            {
                t[tid].dconf->batch_run_override = 0;
            }
        }
        else
        {
            t[tid].dconf->batch_run_override = 0;
        }

    }
    else
    {
        static_conf->flag = 1;
		//printf("thread %d stopping work\n", tid);
        //printf("thread %d stopping work: found %d, %d in flight, need %d\n",
        //    tid, udata->num_found, (int)in_flight_rels, udata->num_needed);
        tdata->work_fcn_id = tdata->num_work_fcn;
        //printf("stopping work with %d poly-a's\n", static_conf->total_poly_a + 1);
    }


    return;
}

void siqsexit(int sig)
{
    if (SIQS_ABORT == 1)
    {
        exit(1);
    }
    printf("\nAborting... threads will finish their current poly\n");
    printf("Press Ctrl-C again to exit immediately and lose in-progress data\n");
    SIQS_ABORT = 1;
    return;
}


// rare issues:
// 1) when thread count is high, yafu can decide to stop because
// of in-flight relations.  If these prove to not be enough to build
// a matrix, wh enter a re-build loop which fails.
// 2) on windows, I've seen the non-trivial dependencies count
// be fairly low (10 or 11).  In some cases no factors are found with these.
// then we enter a loop where the savefile is continually re-scanned
// and the matrix re-fails, forever.
//
// I think both cases can be fixed by improving the ability to go back
// and find more relations if something goes wrong during post-processing
// and no factors are found.


void SIQS(fact_obj_t *fobj)
{
	// the input fobj->N and this 'n' are pointers to memory which holds
	// the input number to this function.  a copy in different memory
	// will also be made, and modified by a multiplier.  don't confuse
	// these two.
	// input expected in fobj->gmp_n	

	// thread data holds all data needed during sieving
	tpool_t *tpool_data;
	siqs_userdata_t udata;
	thread_sievedata_t *thread_data;		//an array of thread data objects

	// master control structure
	static_conf_t *static_conf;

	// stuff for lanczos
	qs_la_col_t *cycle_list;
	uint32_t num_cycles = 0;
	uint64_t *bitfield = NULL;
	siqs_r *relation_list;

	// log file
	FILE *sieve_log;

	// some locals	
	int i, j, k;
	clock_t start, stop;
	double t_time;
	struct timeval myTVend;
    int max_retries = 5;

	// checking savefile
	FILE *data;
	char tmpstr[GSTR_MAXSIZE];

	// logfile for this factorization
	// must ensure it is only written to by main thread
	if (fobj->LOGFLAG && (fobj->flags != 12345))
	{
		fobj->logfile = fopen(fobj->flogname, "a");
		sieve_log = fobj->logfile;
	}
	else
	{
		fobj->logfile = NULL;
		sieve_log = fobj->logfile;
	}

    if (!siqs_primes_initialized)
    {
        soe_staticdata_t* sdata = soe_init(0, 1, 32768);
        siqs_primes = soe_wrapper(sdata, 0, 100000000, 0, &siqs_nump, 0, 0);
        siqs_maxp = siqs_primes[siqs_nump - 1];
        soe_finalize(sdata);
        siqs_primes_initialized = 1;
    }

	fobj->bits = mpz_sizeinbase(fobj->qs_obj.gmp_n, 2);
	fobj->digits = gmp_base10(fobj->qs_obj.gmp_n);

	// print the "starting" line before we check for special cases.
	// other applications (e.g., aliqeit) check for the starting
	// line when parsing the factor.log file.  Should have always 
	// printed it anyway, for consistency.
	if (fobj->VFLAG >= 0)
	{
		gmp_printf("\nstarting SIQS on c%d: %Zd\n",
			fobj->digits, fobj->qs_obj.gmp_n);
	}

    // test/debug: force a particular random seed
    // fobj->lcg_state = 1946832636881765315ULL;
    // printf("warning: test lcg_state is enabled!\n");

    fobj->seed1 = (uint32_t)(fobj->lcg_state & 0xffffffff );
    fobj->seed2 = (uint32_t)(fobj->lcg_state >> 32);

	if (sieve_log != NULL)
	{
        char* s = mpz_get_str(NULL, 10, fobj->qs_obj.gmp_n);
		logprint(sieve_log, "starting SIQS on c%d: %s\n", fobj->digits, s);
        free(s);
		logprint(sieve_log, "random seed: %" PRIu64 "\n", fobj->lcg_state);
		fflush(sieve_log);
	}

	// check for special cases and bail if there is one
	if ((i = check_specialcase(fobj->logfile, fobj)) > 0)
	{
		if (i == 1)
		{
			if (sieve_log != NULL)
				fclose(sieve_log);
		}
		return;
	}

	// check to see if a siqs savefile exists for this input	
	data = fopen(fobj->qs_obj.siqs_savefile, "r");

	if (data != NULL)
	{
		char *substr;
		mpz_t tmpz;
		mpz_t g;

		//read in the number from the savefile
		mpz_init(tmpz);
		mpz_init(g);

		fgets(tmpstr, GSTR_MAXSIZE, data);
		substr = tmpstr + 2;
		mpz_set_str(tmpz, substr, 0);	//auto detect the base

		if (resume_check_input_match(tmpz, fobj->qs_obj.gmp_n, g, fobj->VFLAG))
		{
			// remove any common factor so the input exactly matches
			// the file.  
			// TODO: we should also factor the common factor 'g' 
			// somehow, since the user may be expecting the complete
			// factorization of the input instead of just the 
			// qs-resumed bit.
			if (mpz_cmp_ui(g, 1) > 0)
			{
				add_to_factor_list(fobj->factors, g, fobj->VFLAG, fobj->NUM_WITNESSES);
			}
			mpz_tdiv_q(fobj->qs_obj.gmp_n, fobj->qs_obj.gmp_n, g);
			//mpz_set(fobj->N, fobj->gmp_n);
		}
		mpz_clear(tmpz);
		mpz_clear(g);
		fclose(data);
	}

	// At this point, we are committed to doing qs on the input
	// we need to:
	// 1.) notify the screen and log
	// 2.) start watching for an abort signal
	// 3.) initialize data objects
	// 4.) get ready to find the factor base

	// fill in the factorization object	
	fobj->qs_obj.savefile.name = (char *)malloc(80 * sizeof(char));
	strncpy(fobj->savefile_name, fobj->qs_obj.siqs_savefile, 80);

    // The above leaks a little memory:
    // == 148041 == 80 bytes in 1 blocks are definitely lost in loss record 2 of 2
    // == 148041 == at 0x4C2B0F7: malloc(vg_replace_malloc.c:381)
    // == 148041 == by 0x42E5E7 : SIQS(SIQS.c:687)
    // == 148041 == by 0x420046 : feval(calc.c:2575)
    // == 148041 == by 0x41F175 : calc(calc.c:1947)
    // == 148041 == by 0x41CBE2 : calc_with_assignment(calc.c:1527)
    // == 148041 == by 0x41CBE2 : process_expression(calc.c:1473)
    // == 148041 == by 0x405FD2 : main(driver.c:401)


	// initialize the data objects both shared (static) and 
	// per-thread (dynamic)
	static_conf = (static_conf_t *)malloc(sizeof(static_conf_t));
	static_conf->obj = fobj;

	thread_data = (thread_sievedata_t *)malloc(fobj->THREADS * sizeof(thread_sievedata_t));
	for (i = 0; i < fobj->THREADS; i++)
	{
		thread_data[i].dconf = (dynamic_conf_t *)malloc(sizeof(dynamic_conf_t));
        thread_data[i].dconf->tid = i;
		thread_data[i].sconf = static_conf;
	}

	// initialize the flag to watch for interrupts, and set the
	// pointer to the function to call if we see a user interrupt
	SIQS_ABORT = 0;
	signal(SIGINT, siqsexit);

	// start a counter for the whole job
	gettimeofday(&static_conf->totaltime_start, NULL);

	udata.thread_data = thread_data;

	tpool_data = tpool_setup(fobj->THREADS, NULL, NULL,
		&siqs_sync, &siqs_dispatch, &udata);
	tpool_add_work_fcn(tpool_data, &process_poly);

	// this function is not run by tpool. It puts things into the user
	// data portion of the tpool structure and runs a few initialization 
	// routines in an attempt to clean up this toplevel function.
	siqs_start(tpool_data);

#if 0
    if (0)
    {
        // hack in a loop to examine the size characteristics for 
        // the minimum offset for many different polynomials.
        mpz_t Q;
        mpz_init(Q);
        dynamic_conf_t* dconf = thread_data[0].dconf;
        int bins[24];
        int more = 0;
        gmp_printf("sqrtN = %Zd\n", static_conf->sqrt_n);
        int threshold;

        for (i = 0; i < 24; i++)
            bins[i] = 0;


        int num_small = 0;
        dconf->tot_poly = 0;
        for (i = 0; i < 100000; i++)
        {
            static_conf->total_poly_a++;
            new_poly_a(static_conf, dconf);

            // update the gray code
            get_gray_code(dconf->curr_poly);

            // update roots, etc.
            dconf->maxB = 1 << (dconf->curr_poly->s - 1);
            dconf->numB = 1;
            computeBl(static_conf, dconf, 1);

            for (; dconf->numB < dconf->maxB; dconf->numB++, dconf->tot_poly++)
            {
                mpz_sub(dconf->gmptmp1,
                    static_conf->sqrt_n, dconf->curr_poly->mpz_poly_b);
                mpz_tdiv_q(dconf->gmptmp1, dconf->gmptmp1, dconf->curr_poly->mpz_poly_a);
                if (mpz_sgn(dconf->gmptmp1) < 0)
                    mpz_neg(dconf->gmptmp1, dconf->gmptmp1);

                uint64_t minoffset = mpz_get_ui(dconf->gmptmp1);
                int parity;

                for (parity = 0; parity <= 1; parity++);
                {
                    mpz_mul_2exp(dconf->gmptmp2, dconf->curr_poly->mpz_poly_b, 1);
                    mpz_mul_ui(dconf->gmptmp1, dconf->curr_poly->mpz_poly_a, minoffset);

                    if (parity)
                        mpz_sub(Q, dconf->gmptmp1, dconf->gmptmp2);
                    else
                        mpz_add(Q, dconf->gmptmp1, dconf->gmptmp2);

                    mpz_mul_ui(Q, Q, minoffset);
                    mpz_add(Q, Q, dconf->curr_poly->mpz_poly_c);

                    if (mpz_sgn(Q) < 0)
                    {
                        mpz_neg(Q, Q);
                    }

                    mpz_tdiv_q(dconf->gmptmp2, static_conf->sqrt_n, Q);
                    if (mpz_sizeinbase(dconf->gmptmp2, 2) > 23)
                    {
                        printf("rare size %d bits!\n", mpz_sizeinbase(dconf->gmptmp2, 2));
                        gmp_printf("sqrtN = %Zd\n", static_conf->sqrt_n);
                        gmp_printf("Q     = %Zd\n", Q);
                        more++;
                    }
                    else if (mpz_sizeinbase(dconf->gmptmp2, 2) >= 1)
                    {
                        bins[mpz_sizeinbase(dconf->gmptmp2, 2) - 1]++;
                    }
                }

                nextB(dconf, static_conf, 1);
            }
        }

        for (i = 0; i < 24; i++)
        {
            printf("found %d small Q's (%1.2f percent) with Q >= %d bits below sqrtN\n", bins[i],
                (double)bins[i] / dconf->tot_poly * 100, i + 1);
        }
        if (more > 0)
            printf("found %d small Q's (%1.2f percent) with Q >= %d bits below sqrtN\n", more,
                (double)bins[i] / dconf->tot_poly * 100, i + 1);

        mpz_clear(Q);
        exit(1);
    }
#endif


#ifdef GATHER_RESIDUE_STATS
    char fname[1024];
    
    for (i = 0; i < fobj->THREADS; i++)
    {
        sprintf(fname, "residues_B_%u_MFBT_%1.2f_TF_%u_LPB_%u_tid_%d.txt",
            static_conf->factor_base->B, static_conf->tlp_exp, static_conf->tf_closnuf,
            static_conf->large_prime_max, i);
        static_conf->residue_files[i] = fopen(fname, "a");
        if (static_conf->residue_files[i] == NULL)
        {
            printf("problem opening %s to append\n", fname);
            exit(1);
        }
    }

#endif

    uint32_t factors_found = 0;
    int num_retries = 0;
    do
    {


        if (fobj->THREADS == 1)
        {
            // it is noticably faster to remove the tpool overhead
            // if we just have one thread.  This is basically what
            // tpool_go() does without all of the threading overhead.
            // todo: maybe this *should* be what tpool_go() does when
            // num_threads == 1...
            while (1)
            {
                siqs_sync(tpool_data);
                siqs_dispatch(tpool_data);
                if (tpool_data->work_fcn_id == 0)
                {
                    process_poly(tpool_data);
                }
                else
                {
                    break;
                }
            }
        }
        else
        {
            tpool_go(tpool_data);
        }


        uint32_t missed_batch_rels = 0;
        for (i = 0; i < fobj->THREADS; i++)
        {
            missed_batch_rels += thread_data[i].dconf->rb.num_relations;
        }

        if ((static_conf->do_batch) && (fobj->VFLAG > 0) && (missed_batch_rels > 0))
        {
            printf("threw away %u batched relations\n", missed_batch_rels);
        }

        if (!static_conf->in_mem)
        {
            // finalize savefile
            savefile_flush(&static_conf->obj->qs_obj.savefile);
            savefile_close(&static_conf->obj->qs_obj.savefile);
        }

        update_final(static_conf);

        if (udata.updatecode == 2)
        {
            // a stop condition, other than finding enough relations to finish, was met.
            // e.g., specified elasped time or relation count or a CTRL-C.  In these
            // cases don't run filtering, just go to cleanup.
            goto done;
        }

        gettimeofday(&myTVend, NULL);
        t_time = ytools_difftime(&static_conf->totaltime_start, &myTVend);

        if (fobj->VFLAG > 0)
        {
            if (static_conf->use_dlp >= 2)
                printf("TLP filtering time = %6.4f seconds.\n", static_conf->t_time4);
            printf("QS elapsed time = %6.4f seconds.\n", t_time);
            printf("\n==== post processing stage (msieve-1.38) ====\n");
        }

        if (sieve_log != NULL)
        {
            if (static_conf->use_dlp >= 2)
                logprint(sieve_log, "TLP filtering time = %6.4f seconds.\n", static_conf->t_time4);
            logprint(sieve_log, "QS elapsed time = %6.4f seconds.\n", t_time);
        }

        fobj->qs_obj.qs_time = t_time;

        start = clock();

        // filtering will reload (if using a savefile) and possibly re-order
        // the a-polys.  If using a save file we clear our current list
        // so that filtering can do it's thing.  If in-mem, keep them; filtering
        // knows how to deal with that case too.
        if ((!static_conf->in_mem) && (static_conf->total_poly_a > 0))
        {
            for (i = 0; (uint32_t)i < static_conf->total_poly_a; i++)
            {
                mpz_clear(static_conf->poly_a_list[i]);
            }
            free(static_conf->poly_a_list);
        }

        // load and filter relations and polys.
        yafu_qs_filter_relations(thread_data[0].sconf);

        cycle_list = static_conf->cycle_list;
        num_cycles = static_conf->num_cycles;
        relation_list = static_conf->relation_list;

        // solve the system of equations
        j = qs_solve_linear_system(static_conf->obj, static_conf->factor_base->B,
            &bitfield, relation_list, cycle_list, &num_cycles);

        if ((j == -1) || (j == -2))
        {
            // not enough dependencies or corrupt matrix.
            if (num_retries >= max_retries)
            {
                printf("Problem with filtering and too many retries, giving up\n");
                logprint(sieve_log, "Problem with filtering and too many "
                    "retries, giving up\n");
            }
            else
            {
                uint32_t old_needed = udata.num_needed;
                udata.num_needed = udata.num_needed * 1.05;
                static_conf->num_extra_relations += (udata.num_needed - old_needed);
                printf("Problem with filtering (code %d), re-trying with new relation "
                    "target of %u\n", j, udata.num_needed);
                logprint(sieve_log, "Problem with filtering (code %d), re-trying with new relation "
                    "target of % u\n", j, udata.num_needed);

                // also update the RNG as that could fix it too
                lcg_rand_64(&fobj->lcg_state);
                fobj->seed1 = (uint32_t)(fobj->lcg_state & 0xffffffff);
                fobj->seed2 = (uint32_t)(fobj->lcg_state >> 32 );

                //initialize the bookkeeping for tracking partial relations
                free(static_conf->cycle_hashtable);
                free(static_conf->cycle_table);
                static_conf->components = 0;
                static_conf->vertices = 0;
                static_conf->num_cycles = 0;
                static_conf->num_relations = 0;

                static_conf->cycle_hashtable = (uint32_t*)xcalloc(
                    (size_t)(1 << QS_LOG2_CYCLE_HASH),
                    sizeof(uint32_t));
                static_conf->cycle_table_size = 1;
                static_conf->cycle_table_alloc = 10000;
                static_conf->cycle_table = (qs_cycle_t*)xmalloc(
                    static_conf->cycle_table_alloc * sizeof(qs_cycle_t));

                if (!static_conf->in_mem)
                {
                    restart_siqs(static_conf, thread_data[0].dconf);
                    savefile_open(&static_conf->obj->qs_obj.savefile, SAVEFILE_APPEND);
                }
                else
                {
                    printf("rebuilding graph with %u buffered in-memory relations\n",
                        static_conf->buffered_rels);
                    rebuild_graph(static_conf, static_conf->in_mem_relations,
                        static_conf->buffered_rels);
                }

                num_retries++;
                continue;
            }
        }

        stop = clock();
        static_conf->t_time1 = (double)(stop - start) / (double)CLOCKS_PER_SEC;

        start = clock();
        //sqrt stage
        if (bitfield != NULL && num_cycles > 0)
        {

            factors_found = yafu_find_factors(static_conf->obj,
                static_conf->n, static_conf->factor_base->list,
                static_conf->factor_base->B, cycle_list, num_cycles,
                relation_list, bitfield, static_conf->multiplier,
                static_conf->poly_a_list, static_conf->poly_list,
                &static_conf->factor_list);
        }

        stop = clock();
        static_conf->t_time2 = (double)(stop - start) / (double)CLOCKS_PER_SEC;

        gettimeofday(&myTVend, NULL);
        static_conf->t_time3 = ytools_difftime(&static_conf->totaltime_start, &myTVend);

        if (fobj->VFLAG > 0)
        {
            printf("Lanczos elapsed time = %6.4f seconds.\n", static_conf->t_time1);
            printf("Sqrt elapsed time = %6.4f seconds.\n", static_conf->t_time2);
        }

        if (factors_found == 0)
        {
            // we failed to find any factors... too few dependencies,
            // none of them worked, or something else went wrong.
            // There are a couple things to try: either take away some
            // relations or add some more.  Most often adding more
            // seems to help.  We can/should also keep track of this
            // failure and somehow revert back to ecm or maybe nfs if
            // QS keeps failing.
            if (num_retries >= max_retries)
            {
                printf("Failed to find factors and too many retries, giving up\n");
                logprint(sieve_log, "Failed to find factors and too many "
                    "retries, giving up\n");
            }
            else
            {
                uint32_t old_needed = udata.num_needed;
                udata.num_needed = udata.num_needed * 1.05;
                static_conf->num_extra_relations += (udata.num_needed - old_needed);
                printf("Failed to find factors, re-trying with new relation target of %u\n",
                    udata.num_needed);
                logprint(sieve_log, "Failed to find factors, re-trying with "
                    "new relation target of % u\n", udata.num_needed);

                // also update the RNG as that could fix it too
                lcg_rand_64(&fobj->lcg_state);
                fobj->seed1 = (uint32_t)(fobj->lcg_state & 0xffffffff);
                fobj->seed2 = (uint32_t)(fobj->lcg_state >> 32 );

                //initialize the bookkeeping for tracking partial relations
                free(static_conf->cycle_hashtable);
                free(static_conf->cycle_table);
                static_conf->components = 0;
                static_conf->vertices = 0;
                static_conf->num_cycles = 0;
                static_conf->num_relations = 0;

                static_conf->cycle_hashtable = (uint32_t*)xcalloc(
                    (size_t)(1 << QS_LOG2_CYCLE_HASH),
                    sizeof(uint32_t));
                static_conf->cycle_table_size = 1;
                static_conf->cycle_table_alloc = 10000;
                static_conf->cycle_table = (qs_cycle_t*)xmalloc(
                    static_conf->cycle_table_alloc * sizeof(qs_cycle_t));

                if (!static_conf->in_mem)
                {
                    restart_siqs(static_conf, thread_data[0].dconf);
                    savefile_open(&static_conf->obj->qs_obj.savefile, SAVEFILE_APPEND);
                }
                else
                {
                    printf("rebuilding graph with %u buffered in-memory relations\n",
                        static_conf->buffered_rels);
                    rebuild_graph(static_conf, static_conf->in_mem_relations, 
                        static_conf->buffered_rels);
                }
                num_retries++;
            }
        }
    } while ((factors_found == 0) && (num_retries < max_retries));

	if (fobj->VFLAG >= 0)
	{
		printf("SIQS elapsed time = %6.4f seconds.\n", static_conf->t_time3);
	}

	fobj->qs_obj.total_time = static_conf->t_time3;

	if (sieve_log != NULL)
	{
		logprint(sieve_log,"Lanczos elapsed time = %6.4f seconds.\n",static_conf->t_time1);
		logprint(sieve_log,"Sqrt elapsed time = %6.4f seconds.\n",static_conf->t_time2);
		logprint(sieve_log,"SIQS elapsed time = %6.4f seconds.\n",static_conf->t_time3);
		logprint(sieve_log,"\n");
		logprint(sieve_log,"\n");
	}

	static_conf->cycle_list = cycle_list;
	static_conf->num_cycles = num_cycles;

    // free stuff used during filtering
    free_filter_vars(static_conf);

    // and during sqrt
    if ((bitfield != NULL) && (num_cycles > 0))
    {
        free(bitfield);
    }

done:

#ifdef GATHER_RESIDUE_STATS

    for (i = 0; i < fobj->THREADS; i++)
    {
        fclose(static_conf->residue_files[i]);
    }

#endif

    free(tpool_data);

    // free dynamic thread data
    for (i = 0; i < fobj->THREADS; i++)
    {
#if defined(USE_SS_SEARCH)
        if ((static_conf->factor_base->ss_start_B > 0) &&
            (static_conf->factor_base->ss_start_B < static_conf->factor_base->B))
        {
            free(thread_data[i].dconf->ss_sieve_p);
            free(thread_data[i].dconf->ss_sieve_n);
        }
        else
        {
            free_sieve(thread_data[i].dconf);
            
        }
#else
        free_sieve(thread_data[i].dconf);
#endif
        
        free(thread_data[i].dconf->relation_buf);
    }

	//free everything else
	free_siqs(thread_data[0].sconf);

    if (sieve_log != NULL)
    {
        fclose(sieve_log);
    }

	for (i=0; i<fobj->THREADS; i++)
	{
		free(thread_data[i].dconf);
	}
	free(static_conf);
	free(thread_data);

	//reset signal handler to default (no handler).
	signal(SIGINT,NULL);

    if (num_retries > 0)
    {
        // during testing for successful retries, stop here
        //exit(0);
    }

	return;
}

#if defined( USE_SS_SEARCH ) && defined( USE_DIRECT_SIEVE_SS )

int check_Qval(static_conf_t* sconf, dynamic_conf_t* dconf,
    int polyid, int offset, int parity)
{

    if (sconf->knmod8 == 1)
    {
        // this one is close enough, compute 
        // Q(x) = (2ax + b)^2 - N, where x is the sieve index
        // Q(x)/4a = (ax + b)x + c;	
        mpz_mul_ui(dconf->gmptmp1, dconf->curr_poly->mpz_poly_a, offset);

        if (parity)
            mpz_sub(dconf->gmptmp1, dconf->gmptmp1, dconf->curr_poly->mpz_poly_b);
        else
            mpz_add(dconf->gmptmp1, dconf->gmptmp1, dconf->curr_poly->mpz_poly_b);

        mpz_mul_ui(dconf->gmptmp1, dconf->gmptmp1, offset);
        mpz_add(dconf->gmptmp1, dconf->gmptmp1, dconf->curr_poly->mpz_poly_c);

        if (mpz_sgn(dconf->gmptmp1) < 0)
        {
            mpz_neg(dconf->gmptmp1, dconf->gmptmp1);
        }
    }
    else
    {
        // this one is close enough, compute 
        // Q(x) = (ax + b)^2 - N, where x is the sieve index
        // Q(x)/a = (ax + 2b)x + c;	
        mpz_mul_2exp(dconf->gmptmp2, dconf->curr_poly->mpz_poly_b, 1);
        mpz_mul_ui(dconf->gmptmp1, dconf->curr_poly->mpz_poly_a, offset);

        if (parity)
            mpz_sub(dconf->gmptmp1, dconf->gmptmp1, dconf->gmptmp2);
        else
            mpz_add(dconf->gmptmp1, dconf->gmptmp1, dconf->gmptmp2);

        mpz_mul_ui(dconf->gmptmp1, dconf->gmptmp1, offset);
        mpz_add(dconf->gmptmp1, dconf->gmptmp1, dconf->curr_poly->mpz_poly_c);

        if (mpz_sgn(dconf->gmptmp1) < 0)
        {
            mpz_neg(dconf->gmptmp1, dconf->gmptmp1);
        }
    }

    uint32_t fboffsets[100];
    uint8_t sieveval = sconf->blockinit;
    int i;
    int j = 0;
    mpz_set(dconf->gmptmp2, dconf->gmptmp1);
    for (i = sconf->sieve_small_fb_start; i < sconf->factor_base->B; i++)
    {
        uint32_t prime = sconf->factor_base->list->prime[i];
        uint8_t logp = sconf->factor_base->list->logprime[i];

        while (mpz_tdiv_ui(dconf->gmptmp2, prime) == 0)
        {
            sieveval -= logp;
            fboffsets[j++] = i;
            mpz_tdiv_q_ui(dconf->gmptmp2, dconf->gmptmp2, prime);
        }
    }

    if ((sieveval & 0x80))
    {
        gmp_printf("Q(%d,%d), logp = %02x: %Zd = %u ", polyid, offset, 
            sieveval, dconf->gmptmp1,
            sconf->factor_base->list->prime[fboffsets[0]]);
        for (i = 1; i < j; i++)
        {
            printf("* %u ", sconf->factor_base->list->prime[fboffsets[i]]);
        }

        while (mpz_even_p(dconf->gmptmp2))
            mpz_tdiv_q_2exp(dconf->gmptmp2, dconf->gmptmp2, 1);

        for (i = 2; i < sconf->sieve_small_fb_start; i++)
        {
            uint32_t prime = sconf->factor_base->list->prime[i];
            if (mpz_tdiv_ui(dconf->gmptmp1, prime) == 0)
                mpz_tdiv_q_ui(dconf->gmptmp2, dconf->gmptmp2, prime);
        }
        gmp_printf("* %Zd\n", dconf->gmptmp2);
    }

    return ((sieveval & 0x80) > 0);
}

void init_Qval(static_conf_t* sconf, dynamic_conf_t* dconf,
    int polyid, int offset, int parity, int report_num)
{
    dconf->reports[report_num] = offset;
    dconf->num++;
    dconf->total_reports++;

    if (sconf->knmod8 == 1)
    {
        // this one is close enough, compute 
        // Q(x) = (2ax + b)^2 - N, where x is the sieve index
        // Q(x)/4a = (ax + b)x + c;	
        mpz_mul_ui(dconf->gmptmp1, dconf->curr_poly->mpz_poly_a, offset);

        if (parity)
            mpz_sub(dconf->Qvals[report_num], dconf->gmptmp1, dconf->curr_poly->mpz_poly_b);
        else
            mpz_add(dconf->Qvals[report_num], dconf->gmptmp1, dconf->curr_poly->mpz_poly_b);

        mpz_mul_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], offset);
        mpz_add(dconf->Qvals[report_num], dconf->Qvals[report_num], dconf->curr_poly->mpz_poly_c);

        if (mpz_sgn(dconf->Qvals[report_num]) < 0)
        {
            mpz_neg(dconf->Qvals[report_num], dconf->Qvals[report_num]);
        }

        // minimum Q(x) should occur when (2ax + b)^2 = N
        // 2ax + b = sqrt(N)
        // 2ax = sqrt(N) - b
        // x = (sqrt(N) - b) / 2a

        //mpz_sqrt(dconf->gmptmp1, sconf->n);
        //mpz_sub(dconf->gmptmp1, dconf->gmptmp1, dconf->curr_poly->mpz_poly_b);
        //mpz_tdiv_q(dconf->gmptmp1, dconf->gmptmp1, dconf->curr_poly->mpz_poly_a);
        //mpz_tdiv_q_2exp(dconf->gmptmp1, dconf->gmptmp1, 1);
        //
        //gmp_printf("Min x occurs at offset %Zd, this offset: %d\n", dconf->gmptmp1, offset);
        //gmp_printf("This Q: %Zd\n", dconf->Qvals[report_num]);
        //
        //mpz_mul(dconf->gmptmp1, dconf->curr_poly->mpz_poly_a, dconf->gmptmp1);
        //
        //if (parity)
        //    mpz_sub(dconf->gmptmp1, dconf->gmptmp1, dconf->curr_poly->mpz_poly_b);
        //else
        //    mpz_add(dconf->gmptmp1, dconf->gmptmp1, dconf->curr_poly->mpz_poly_b);
        //
        //mpz_mul_ui(dconf->gmptmp1, dconf->gmptmp1, offset);
        //mpz_add(dconf->gmptmp1, dconf->gmptmp1, dconf->curr_poly->mpz_poly_c);
        //
        //if (mpz_sgn(dconf->gmptmp1) < 0)
        //{
        //    mpz_neg(dconf->gmptmp1, dconf->gmptmp1);
        //}
        //
        //gmp_printf("Min Q : %Zd\n", dconf->gmptmp1);
        //
        //mpz_mul_ui(dconf->gmptmp1, dconf->curr_poly->mpz_poly_a, 0);
        //
        //if (parity)
        //    mpz_sub(dconf->gmptmp1, dconf->gmptmp1, dconf->curr_poly->mpz_poly_b);
        //else
        //    mpz_add(dconf->gmptmp1, dconf->gmptmp1, dconf->curr_poly->mpz_poly_b);
        //
        //mpz_mul_ui(dconf->gmptmp1, dconf->gmptmp1, offset);
        //mpz_add(dconf->gmptmp1, dconf->gmptmp1, dconf->curr_poly->mpz_poly_c);
        //
        //if (mpz_sgn(dconf->gmptmp1) < 0)
        //{
        //    mpz_neg(dconf->gmptmp1, dconf->gmptmp1);
        //}
        //
        //gmp_printf("Q(0)  : %Zd\n", dconf->gmptmp1);
        //exit(1);
    }
    else
    {
        // this one is close enough, compute 
        // Q(x) = (ax + b)^2 - N, where x is the sieve index
        // Q(x)/a = (ax + 2b)x + c;	
        mpz_mul_2exp(dconf->gmptmp2, dconf->curr_poly->mpz_poly_b, 1);
        mpz_mul_ui(dconf->gmptmp1, dconf->curr_poly->mpz_poly_a, offset);

        if (parity)
            mpz_sub(dconf->Qvals[report_num], dconf->gmptmp1, dconf->gmptmp2);
        else
            mpz_add(dconf->Qvals[report_num], dconf->gmptmp1, dconf->gmptmp2);

        mpz_mul_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], offset);
        mpz_add(dconf->Qvals[report_num], dconf->Qvals[report_num], dconf->curr_poly->mpz_poly_c);

        if (mpz_sgn(dconf->Qvals[report_num]) < 0)
        {
            mpz_neg(dconf->Qvals[report_num], dconf->Qvals[report_num]);
        }
    }
}

int td_small_p(static_conf_t* sconf, dynamic_conf_t* dconf,
    int polyid, int offset, int parity, int report_num, poly_bucket_t *pbucket)
{
    uint8_t bits = 0, logp;
    int smooth_num;
    int i;
    tiny_fb_element_siqs* fullfb_ptr, * fullfb = sconf->factor_base->tinylist;
    sieve_fb_compressed* fbc;
    uint32_t tmp1, tmp2, tmp3, tmp4;
    uint64_t q64;
    uint32_t tmp, prime, root1, root2;
    siqs_poly* poly = dconf->curr_poly;

    fullfb_ptr = fullfb;
    if (parity)
    {
        fbc = dconf->comp_sieve_n;
    }
    else
    {
        fbc = dconf->comp_sieve_p;
    }

    dconf->smooth_num[report_num] = 0;
    smooth_num = 0;

    //take care of powers of two
    while (mpz_even_p(dconf->Qvals[report_num]))
    {
        mpz_tdiv_q_2exp(dconf->Qvals[report_num], dconf->Qvals[report_num], 1);
        bits++;
    }

    i = 2;
    // explicitly trial divide by small primes which we have not
    // been sieving.  because we haven't been sieving, their progressions
    // have not been updated and thus we can't use faster methods.
    // fortunately, there shouldn't be many of these to test.

    int k = 0;
    int index_to_skip = poly->qlisort[k];

    // finish up the rest of the small primes
    while ((uint32_t)i < sconf->sieve_small_fb_start)
    {
        uint64_t q64;

        if (i == index_to_skip)
        {
            k++;
            if (k < poly->s)
                index_to_skip = poly->qlisort[k];
            else
                index_to_skip = -1;
            continue;
        }

        prime = fbc->prime[i];
        //root1 = fbc->root1[i];
        //root2 = fbc->root2[i];
        logp = fbc->logp[i];
        //
        //// this is just offset % prime (but divisionless!)
        //tmp = offset + fullfb_ptr->correction[i];
        //q64 = (uint64_t)tmp * (uint64_t)fullfb_ptr->small_inv[i];
        //tmp = q64 >> 32;
        //tmp = offset - tmp * prime;

        // if offset % prime == either root, it's on the progression.  also
        // need to check for the case if root1 or root2 == prime at the same
        // time as offset mod prime = 0.  for small primes, this happens fairly
        // often.  the simple offset % prime check will miss these cases.
        //if (tmp == root1 || tmp == root2)
        {
            while (mpz_tdiv_ui(dconf->Qvals[report_num], prime) == 0)
            {
                mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], prime);
                bits += logp;
            }
        }
        i++;
    }

    if (bits < sconf->tf_small_cutoff)
        dconf->valid_Qs[report_num] = 0;
    else
        dconf->valid_Qs[report_num] = 1;

    if (dconf->valid_Qs[report_num])
    {

        while ((uint32_t)i < sconf->factor_base->med_B)
        {
            uint64_t q64;

            if (i == index_to_skip)
            {
                k++;
                if (k < poly->s)
                    index_to_skip = poly->qlisort[k];
                else
                    index_to_skip = -1;
                continue;
            }

            prime = fbc->prime[i];

            while (mpz_tdiv_ui(dconf->Qvals[report_num], prime) == 0)
            {
                dconf->fb_offsets[report_num][smooth_num++] = i;
                mpz_tdiv_q_ui(dconf->Qvals[report_num],
                    dconf->Qvals[report_num], prime);
            }
            i++;
        }

#ifdef USE_AVX512F

        for (i = 0; i < pbucket->numhits; i++)
        {
            uint32_t pid = pbucket->prime[i];
            prime = sconf->factor_base->list->prime[pid];

            while (mpz_tdiv_ui(dconf->Qvals[report_num], prime) == 0)
            {
                dconf->fb_offsets[report_num][smooth_num++] = pid;
                mpz_tdiv_q_ui(dconf->Qvals[report_num],
                    dconf->Qvals[report_num], prime);
            }

        }
#endif
    }
   
    dconf->smooth_num[report_num] = smooth_num;
    return dconf->valid_Qs[report_num];
}

#endif

void* process_poly(void* vptr)
{
    // top-level sieving function which performs all work for a single
    // new a coefficient.  has pthread calling conventions, meant to be
    // used in a multi-threaded environment
    tpool_t* tdata = (tpool_t*)vptr;
    siqs_userdata_t* udata = tdata->user_data;
    thread_sievedata_t* thread_data = &udata->thread_data[tdata->tindex];
    static_conf_t* sconf = thread_data->sconf;
    dynamic_conf_t* dconf = thread_data->dconf;

    // unpack stuff from the job data structure
    sieve_fb_compressed* fb_sieve_p = dconf->comp_sieve_p;
    sieve_fb_compressed* fb_sieve_n = dconf->comp_sieve_n;
    siqs_poly* poly = dconf->curr_poly;
    uint8_t* sieve = dconf->sieve;
    fb_list* fb = sconf->factor_base;
    lp_bucket* buckets = dconf->buckets;
    uint32_t start_prime = sconf->sieve_small_fb_start;
    uint32_t num_blocks = sconf->num_blocks;
    uint8_t blockinit = sconf->blockinit;
    int using_ss_search = 0;

    // locals
    uint32_t i;

    // to get relations per second
    double t_time;
    struct timeval start, stop, st;

    double t_sieve_ss_buckets = 0.0;
    double t_update_buckets = 0.0;
    double t_ss_search_total = 0.0;
    struct timeval start1, stop1;

    sconf->t_time4 = 0.0;

    if ((sconf->factor_base->ss_start_B > 0) &&
        (sconf->factor_base->ss_start_B < sconf->factor_base->B))
    {
        using_ss_search = 1;
    }


    // this routine is handed a dconf structure which already has a
    // new poly a coefficient (and some supporting data).  continue from
    // there, first initializing the gray code...
    gettimeofday(&start, NULL);

    // used to print a little more status info for huge jobs.
    if (sconf->digits_n > 110)
        gettimeofday(&st, NULL);

    // update the gray code
    get_gray_code(dconf->curr_poly);

    // update roots, etc.
    dconf->maxB = 1 << (dconf->curr_poly->s - 1);
    dconf->numB = 1;
    computeBl(sconf, dconf, 1);

    


#if defined( USE_SS_SEARCH )
    if (using_ss_search)
    {
        ss_search_setup(sconf, dconf);
    }
#endif

    firstRoots_ptr(sconf, dconf);

#if defined( USE_SS_SEARCH ) && defined( USE_POLY_BUCKET_SS )

    if (using_ss_search)
    {
#ifdef SS_POLY_BUCKET_SMALL_GROUPS
        //ss_search_sort_set_1(sconf, dconf);
        //for (i = 1; i <= (1 << dconf->ss_set2.size); i++)
        //    ss_search_poly_buckets_2(sconf, dconf, 1);
#else
        gettimeofday(&start1, NULL);
        ss_search_poly_buckets_sorted(sconf, dconf);
        gettimeofday(&stop1, NULL);
        t_ss_search_total += ytools_difftime(&start1, &stop1);
#endif
    }

#endif

#if defined( USE_SS_SEARCH ) && defined( USE_DIRECT_SIEVE_SS )

    if (using_ss_search)
    {
        int p, s;
        int num_bpoly = 1 << (dconf->curr_poly->s - 1);
        int locs_to_resieve = 0;
        int sieve_sz = dconf->ss_sieve_sz;
        update_t* update_data = &dconf->update_data;
        int* rootupdates = dconf->rootupdates;
        uint32_t bound = sconf->factor_base->B;

        poly_bucket_t polybucket;
        polybucket.hitalloc = 1 << 24;
        polybucket.hitloc = (uint32_t*)xmalloc(polybucket.hitalloc * sizeof(uint32_t));
        polybucket.prime = (uint32_t*)xmalloc(polybucket.hitalloc * sizeof(uint32_t));
        polybucket.numhits = 0;

        dconf->ss_sieve_p = (uint8_t*)xrealloc(dconf->ss_sieve_p, num_bpoly *
            dconf->ss_sieve_sz * sizeof(uint8_t));
        dconf->ss_sieve_n = (uint8_t*)xrealloc(dconf->ss_sieve_n, num_bpoly *
            dconf->ss_sieve_sz * sizeof(uint8_t));

        memset(dconf->report_ht_p, 0, dconf->report_ht_size * sizeof(uint16_t));
        memset(dconf->report_ht_n, 0, dconf->report_ht_size * sizeof(uint16_t));
        memset(dconf->ss_sieve_p, blockinit, num_bpoly * sieve_sz * sizeof(uint8_t));
        memset(dconf->ss_sieve_n, blockinit, num_bpoly * sieve_sz * sizeof(uint8_t));

        // for primes up to this bound, do a normal sieve process
        // into the sieve region for all polys.
        // p is the gray-code enumeration.

#ifdef USE_AVX512F
        uint32_t inc[16] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
        __m512i vinterval = _mm512_set1_epi32(sieve_sz);
        __m512i vinc = _mm512_loadu_epi32(inc);
#endif

        // this does the bulk of the sieving.  All polynomials and
        // all primes > 15-bit
        dconf->num_ss_slices = 0;   // indicator for resieveing
        ss_search_poly_buckets_sorted(sconf, dconf);

        gettimeofday(&stop, NULL);
        t_time = ytools_difftime(&start, &stop);

#ifdef USE_AVX512F
        printf("subset-sum sieving of %u primes took: %1.4f seconds\n",
            fb->B - fb->x2_large_B, t_time);
#else
        printf("subset-sum sieving of %u primes took: %1.4f seconds\n",
            fb->B - fb->med_B, t_time);
#endif

        gettimeofday(&start, NULL);

        // last poly is invalid but subset-sum will still hit it.
        for (s = 0; s < dconf->ss_sieve_sz; s++)
        {
            dconf->ss_sieve_p[dconf->polymap[num_bpoly] * sieve_sz + s] = 0xff;
            dconf->ss_sieve_n[dconf->polymap[num_bpoly] * sieve_sz + s] = 0xff;
        }

        int report_num = 1;
        int valid_report_num = 0;
        int collisions = 0;
        dconf->numB = 1;
        for (p = 1; p < num_bpoly; p++)
        {
            int k = 0;
            int index_to_skip = poly->qlisort[k];
            int sign = dconf->polysign[p];
            int v = dconf->polyv[p];

            // map the gray-code enumeration to binary-polynum encoding,
            // because that's the way the subset-sum code wants it and
            // we need to agree with that here.
            int pidx = dconf->polymap[p];

            polybucket.numhits = 0;

            for (i = sconf->sieve_small_fb_start; i < fb->med_B; i++)
            {
                if (i == index_to_skip)
                {
                    k++;
                    if (k < poly->s)
                        index_to_skip = poly->qlisort[k];
                    else
                        index_to_skip = -1;
                    continue;
                }

                uint16_t logp = sconf->factor_base->list->logprime[i];
                uint16_t prime = fb_sieve_p->prime[i];
                uint16_t r1 = fb_sieve_p->root1[i];
                uint16_t r2 = fb_sieve_p->root2[i];

                while (r2 < sieve_sz)
                {
                    dconf->ss_sieve_p[pidx * sieve_sz + r2] -= logp;
                    dconf->ss_sieve_p[pidx * sieve_sz + r1] -= logp;
                    r1 += prime;
                    r2 += prime;
                }
                if (r1 < sieve_sz)
                {
                    dconf->ss_sieve_p[pidx * sieve_sz + r1] -= logp;
                }

                r1 = fb_sieve_n->root1[i];
                r2 = fb_sieve_n->root2[i];

                while (r2 < sieve_sz)
                {
                    dconf->ss_sieve_n[pidx * sieve_sz + r2] -= logp;
                    dconf->ss_sieve_n[pidx * sieve_sz + r1] -= logp;
                    r1 += prime;
                    r2 += prime;
                }
                if (r1 < sieve_sz)
                {
                    dconf->ss_sieve_n[pidx * sieve_sz + r1] -= logp;
                }

                r1 = (uint32_t)update_data->sm_firstroots1[i];
                r2 = (uint32_t)update_data->sm_firstroots2[i];

                int sign = dconf->polysign[p];
                int v = dconf->polyv[p];
                int* ptr = &rootupdates[(v - 1) * bound + i];

                if (sign > 0)
                {
                    r1 = r1 - *ptr;
                    r2 = r2 - *ptr;
                    r1 = (r1 >> 15) ? r1 + prime : r1;
                    r2 = (r2 >> 15) ? r2 + prime : r2;
                }
                else
                {
                    r1 = r1 + *ptr;
                    r2 = r2 + *ptr;
                    r1 = (r1 >= prime) ? r1 - prime : r1;
                    r2 = (r2 >= prime) ? r2 - prime : r2;
                }

                if (r2 < r1)
                {
                    update_data->sm_firstroots1[i] = (uint16_t)r2;
                    update_data->sm_firstroots2[i] = (uint16_t)r1;

                    fb_sieve_p->root1[i] = (uint16_t)r2;
                    fb_sieve_p->root2[i] = (uint16_t)r1;
                    fb_sieve_n->root1[i] = (uint16_t)(prime - r1);
                    fb_sieve_n->root2[i] = (uint16_t)(prime - r2);
                }
                else
                {
                    update_data->sm_firstroots1[i] = (uint16_t)r1;
                    update_data->sm_firstroots2[i] = (uint16_t)r2;

                    fb_sieve_p->root1[i] = (uint16_t)r1;
                    fb_sieve_p->root2[i] = (uint16_t)r2;
                    fb_sieve_n->root1[i] = (uint16_t)(prime - r2);
                    fb_sieve_n->root2[i] = (uint16_t)(prime - r1);
                }
            }

#ifdef USE_AVX512F
            // if AVX512 is available, use it to do some normal sieving of
            // this polynomial up to an intermediate bound, much larger than the
            // sieve interval but still much smaller than the factor-base bound.

            for (i = fb->med_B; i < fb->x2_large_B; i += 16)
            {
                int* ptr = &rootupdates[(v - 1) * bound + i];

                __m512i vprime;
                __m512i vroot1;
                __m512i vroot2;
                __m512i vpval;
                __m512i vnroot1;
                __m512i vnroot2;
                __mmask16 mask1;
                __mmask16 mask2;
                __m512i vindex = _mm512_set1_epi32(i);
                uint32_t logp = update_data->logp[i];

                vprime = _mm512_load_epi32((__m512i*)(&update_data->prime[i]));
                vroot1 = _mm512_load_epi32((__m512i*)(&update_data->firstroots1[i]));
                vroot2 = _mm512_load_epi32((__m512i*)(&update_data->firstroots2[i]));

                __mmask16 m1 = _mm512_cmplt_epi32_mask(vroot1, vinterval);
                __mmask16 m2 = _mm512_cmplt_epi32_mask(vroot2, vinterval);

                _mm512_mask_compressstoreu_epi32(polybucket.hitloc +
                    polybucket.numhits, m1, vroot1);
                _mm512_mask_compressstoreu_epi32(polybucket.prime +
                    polybucket.numhits, m1,
                    _mm512_add_epi32(vindex, vinc));

                polybucket.numhits += _mm_popcnt_u32(m1);

                _mm512_mask_compressstoreu_epi32(polybucket.hitloc +
                    polybucket.numhits, m2, vroot2);
                _mm512_mask_compressstoreu_epi32(polybucket.prime +
                    polybucket.numhits, m2,
                    _mm512_add_epi32(vindex, vinc));

                polybucket.numhits += _mm_popcnt_u32(m2);
                vnroot1 = _mm512_sub_epi32(vprime, vroot1);
                vnroot2 = _mm512_sub_epi32(vprime, vroot2);

                m1 = _mm512_cmplt_epi32_mask(vnroot1, vinterval);
                m2 = _mm512_cmplt_epi32_mask(vnroot2, vinterval);

                _mm512_mask_compressstoreu_epi32(polybucket.hitloc +
                    polybucket.numhits, m1, _mm512_add_epi32(vprime, vnroot1));
                _mm512_mask_compressstoreu_epi32(polybucket.prime +
                    polybucket.numhits, m1,
                    _mm512_add_epi32(vindex, vinc));

                polybucket.numhits += _mm_popcnt_u32(m1);

                _mm512_mask_compressstoreu_epi32(polybucket.hitloc +
                    polybucket.numhits, m2, _mm512_add_epi32(vprime, vnroot2));
                _mm512_mask_compressstoreu_epi32(polybucket.prime +
                    polybucket.numhits, m2,
                    _mm512_add_epi32(vindex, vinc));

                polybucket.numhits += _mm_popcnt_u32(m2);

                if (polybucket.numhits + 64 > polybucket.hitalloc)
                {
                    printf("polybucket almost full\n");
                }

                if (sign > 0)
                {
                    vpval = _mm512_load_epi32((__m512i*)ptr);
                    mask1 = _mm512_cmp_epu32_mask(vpval, vroot1, _MM_CMPINT_GT);
                    mask2 = _mm512_cmp_epu32_mask(vpval, vroot2, _MM_CMPINT_GT);
                    vroot1 = _mm512_sub_epi32(vroot1, vpval);
                    vroot2 = _mm512_sub_epi32(vroot2, vpval);
                    vroot1 = _mm512_mask_add_epi32(vroot1, mask1, vroot1, vprime);
                    vroot2 = _mm512_mask_add_epi32(vroot2, mask2, vroot2, vprime);
                    _mm512_store_epi32((__m512i*)(&update_data->firstroots1[i]), vroot1);
                    _mm512_store_epi32((__m512i*)(&update_data->firstroots2[i]), vroot2);
                }
                else
                {
                    vpval = _mm512_load_epi32((__m512i*)ptr);
                    vroot1 = _mm512_add_epi32(vroot1, vpval);
                    vroot2 = _mm512_add_epi32(vroot2, vpval);
                    mask1 = _mm512_cmp_epu32_mask(vroot1, vprime, _MM_CMPINT_GE);
                    mask2 = _mm512_cmp_epu32_mask(vroot2, vprime, _MM_CMPINT_GE);
                    vroot1 = _mm512_mask_sub_epi32(vroot1, mask1, vroot1, vprime);
                    vroot2 = _mm512_mask_sub_epi32(vroot2, mask2, vroot2, vprime);
                    _mm512_store_epi32((__m512i*)(&update_data->firstroots1[i]), vroot1);
                    _mm512_store_epi32((__m512i*)(&update_data->firstroots2[i]), vroot2);
                }
            }


            for (i = 0; i < polybucket.numhits; i++)
            {
                uint32_t prime = update_data->prime[polybucket.prime[i]];
                uint32_t logp = update_data->logp[polybucket.prime[i]];

                if (polybucket.hitloc[i] > sieve_sz)
                {
                    dconf->ss_sieve_n[pidx * sieve_sz +
                        (polybucket.hitloc[i] - prime)] -= logp;
                }
                else
                {
                    dconf->ss_sieve_p[pidx * sieve_sz +
                        polybucket.hitloc[i]] -= logp;
                }
            }

#endif

            // search for hits and assign report numbers to hits.
            // initialize Q's for each hit and trial divide up to 
            // the subset-sum starting bound.  Hits can be disqualified 
            // due to the small-prime variation.
            for (s = 0; s < dconf->ss_sieve_sz; s++)
            {
                if (dconf->ss_sieve_p[pidx * sieve_sz + s] & 0x80)
                {
                    // mark this location for re-sieving
                    int ht_idx = hash64(pidx * sieve_sz + s) %
                        dconf->report_ht_size;

                    if (dconf->report_ht_p[ht_idx] > 0)
                        collisions++;

                    dconf->report_ht_p[ht_idx] = report_num++;
                    int rnum = dconf->report_ht_p[ht_idx];

                    init_Qval(sconf, dconf, p - 1, s, 0, rnum);
                    int valid = td_small_p(sconf, dconf, p - 1, s, 0, rnum, &polybucket);

                    if (valid)
                    {
                        dconf->ss_sieve_p[pidx * sieve_sz + s] = 1;
                        valid_report_num++;
                    }
                    else
                    {
                        dconf->ss_sieve_p[pidx * sieve_sz + s] = 0xff;
                    }
                }
                else
                {
                    dconf->ss_sieve_p[pidx * sieve_sz + s] = 0xff;
                }

                if (dconf->ss_sieve_n[pidx * sieve_sz + s] & 0x80)
                {
                    // mark this location for re-sieving
                    int ht_idx = hash64(pidx * sieve_sz + s) %
                        dconf->report_ht_size;

                    dconf->ss_sieve_n[pidx * sieve_sz + s] = 1;
                    if (dconf->report_ht_n[ht_idx] > 0)
                        collisions++;

                    dconf->report_ht_n[ht_idx] = report_num++;
                    int rnum = dconf->report_ht_n[ht_idx];

                    init_Qval(sconf, dconf, p - 1, s, 1, rnum);
                    int valid = td_small_p(sconf, dconf, p - 1, s, 1, rnum, &polybucket);

                    if (valid)
                    {
                        valid_report_num++;
                    }
                    else
                    {
                        dconf->ss_sieve_n[pidx * sieve_sz + s] = 0xff;
                    }
                }
                else
                {
                    dconf->ss_sieve_n[pidx * sieve_sz + s] = 0xff;
                }
            }

            // next polynomial
            // use the stored Bl's and the gray code to find the next b
            nextB(dconf, sconf, 1);
            dconf->numB++;


        }

        gettimeofday(&stop, NULL);
        t_time = ytools_difftime(&start, &stop);

#ifdef USE_AVX512F
        printf("sieving and trial division of %d primes took: %1.4f seconds\n",
            fb->x2_large_B - sconf->sieve_small_fb_start, t_time);
#else
        printf("sieving and trial division of %d primes took: %1.4f seconds\n",
            fb->med_B - sconf->sieve_small_fb_start, t_time);
#endif

        gettimeofday(&start, NULL);

        // and finally process for any large primes and store discovered relations.
        printf("found %d locations to resieve, %d hashtable collisions\n",
            report_num, collisions);

        // now initialize Q and divide out the small primes at the marked locations.
        // here we need to actually increment through the b-polys.
        dconf->total_surviving_reports += valid_report_num;

        // reset the poly so we can scan over them again later.
        computeBl(sconf, dconf, 1);
        dconf->numB = 1;

        // do the resieve.  Same as before but now we divide out the primes
        // for hits at the marked locations.
        if (valid_report_num > 0)
        {
            //printf("commencing subset-sum resieve\n");
            dconf->num_ss_slices = 1;   // indicator for resieveing
            ss_search_poly_buckets_sorted(sconf, dconf);
        }
        
        gettimeofday(&stop, NULL);
        t_time = ytools_difftime(&start, &stop);
        printf("subset-sum re-sieving took: %1.4f seconds\n", t_time);
        gettimeofday(&start, NULL);

        int invalid_report_count = 0;
        for (p = 1; p < num_bpoly; p++)
        {
            // map the gray-code enumeration to binary-polynum encoding
            int pidx = dconf->polymap[p];

            for (s = 0; s < dconf->ss_sieve_sz; s++)
            {
                //int result;
                //result = check_Qval(sconf, dconf, pidx, s, 0);
                //if (result && (dconf->ss_sieve_p[pidx * sieve_sz + s] == 0xff))
                //{
                //    printf("missed potential relation at pidx=%d (p=%d), s = %d, p-side\n",
                //        pidx, p, s);
                //}
                //result = check_Qval(sconf, dconf, pidx, s, 1);
                //if (result && (dconf->ss_sieve_n[pidx * sieve_sz + s] == 0xff))
                //{
                //    printf("missed potential relation at pidx=%d (p=%d), s = %d, n-side\n",
                //        pidx, p, s);
                //}

                if (dconf->ss_sieve_p[pidx * sieve_sz + s] < 0xff)
                {
                    // initialize Qval for poly p at offset s
                    int ht_idx = hash64(pidx * sieve_sz + s) %
                        dconf->report_ht_size;
                    int rnum = dconf->report_ht_p[ht_idx];

                    if (rnum == 0)
                    {
                        invalid_report_count++;
                        continue;
                    }

                    trial_divide_Q_siqs(rnum, 0, p - 1, 0, sconf, dconf);
                }
                if (dconf->ss_sieve_n[pidx * sieve_sz + s] < 0xff)
                {
                    // initialize Qval for poly p at offset -s
                    int ht_idx = hash64(pidx * sieve_sz + s) %
                        dconf->report_ht_size;
                    int rnum = dconf->report_ht_n[ht_idx];

                    if (rnum == 0)
                    {
                        invalid_report_count++;
                        continue;
                    }

                    trial_divide_Q_siqs(rnum, 1, p - 1, 0, sconf, dconf);
                }
            }

            // next polynomial
            // use the stored Bl's and the gray code to find the next b
            nextB(dconf, sconf, 1);
            dconf->numB++;
            dconf->tot_poly++;
        }

        gettimeofday(&stop, NULL);
        t_time = ytools_difftime(&start, &stop);
        if (invalid_report_count > 0)
        {
            printf("%d invalid reports pulled during trial division\n", invalid_report_count);
        }
        printf("large prime trial division took: %1.4f seconds\n", t_time);

        dconf->total_blocks += num_bpoly * 2;

        free(polybucket.hitloc);
        free(polybucket.prime);

        goto done;
    }

#endif


    // loop over each possible b value, for the current a value
    for (; dconf->numB < dconf->maxB; dconf->numB++, dconf->tot_poly++)
    {
        uint32_t invalid_root_marker = 0xFFFFFFFF;
        uint32_t minblock;

        mpz_sub(dconf->gmptmp1, sconf->sqrt_n, dconf->curr_poly->mpz_poly_b);
        mpz_tdiv_q(dconf->gmptmp1, dconf->gmptmp1, dconf->curr_poly->mpz_poly_a);
        if (sconf->knmod8 == 1)
            mpz_tdiv_q_2exp(dconf->gmptmp1, dconf->gmptmp1, 1);

        if (mpz_sgn(dconf->gmptmp1) < 0)
            mpz_neg(dconf->gmptmp1, dconf->gmptmp1);

        mpz_set(dconf->gmptmp3, dconf->gmptmp1);
        mpz_tdiv_q_2exp(dconf->gmptmp1, dconf->gmptmp1, sconf->qs_blockbits);
        minblock = mpz_get_ui(dconf->gmptmp1);

#ifdef SS_POLY_BUCKET_SMALL_GROUPS

        if (using_ss_search &&
            ((dconf->numB % (1 << dconf->ss_set1.size)) == 0))
        {
        //    //printf("running search roots on set2 instance %d\n",
        //    //    dconf->numB / (1 << dconf->ss_set1.size));
            ss_search_poly_buckets_2(sconf, dconf, dconf->numB / (1 << dconf->ss_set1.size));
        }
#endif

        gettimeofday(&start1, NULL);

#if defined( USE_SS_SEARCH ) //&& defined( USE_LINKED_LIST_SS )

        if (using_ss_search)
        {
            for (i = 0; i < num_blocks; i++)
            {
                if (i == minblock)
                    memset(dconf->ss_sieve_p + i * 32768, blockinit - 4, 32768);
                else
                    memset(dconf->ss_sieve_p + i * 32768, blockinit, 32768);
            }

            //for ( ; i < 2 * num_blocks; i++)
            //{
            //    if ((i - num_blocks) == minblock)
            //        memset(dconf->ss_sieve_p + i * 32768, blockinit - 4, 32768);
            //    else
            //        memset(dconf->ss_sieve_p + i * 32768, blockinit, 32768);
            //}

            lp_sieve_ss(dconf->ss_sieve_p, 0, dconf);
        }

#endif

#if 0 //defined( USE_SS_SEARCH )
        if (using_ss_search)
        {
            // test if anything in the block sieve modifies values outside the block
            memcpy(dconf->ss_sieve_n, dconf->ss_sieve_p, (num_blocks + 1) * 32768);
        }
#endif

        

        for (i = 0; i < num_blocks; i++)
        {
#if defined( USE_SS_SEARCH )
            if (using_ss_search)
            {
                //dconf->sieve = dconf->ss_sieve_p + i * 32768 + num_blocks * 32768;
                dconf->sieve = dconf->ss_sieve_p + i * 32768;
                sieve = dconf->sieve;
            }
#endif

            // set the roots for the factors of a such that
            // they will not be sieved.  we haven't found roots for them
            set_aprime_roots(sconf, invalid_root_marker, poly->qlisort, poly->s, fb_sieve_p, 1);
            if (i == minblock)
            {
#if !defined( USE_SS_SEARCH )
                memset(sieve, blockinit - 4, 32768);
#else
                if (!using_ss_search)
                {
                    memset(sieve, blockinit - 4, 32768);
                }
#endif
                med_sieve_ptr(sieve, fb_sieve_p, fb, start_prime, blockinit - 4);
            }
            else
            {
#if !defined( USE_SS_SEARCH )
                memset(sieve, blockinit, 32768);
#else
                if (!using_ss_search)
                {
                    memset(sieve, blockinit, 32768);
                }
#endif
                med_sieve_ptr(sieve, fb_sieve_p, fb, start_prime, blockinit);
            }

#if 0 //defined( USE_SS_SEARCH )
            if (using_ss_search)
            {
                // test if anything in the block sieve modifies values outside the block
                if ((i + 1) < num_blocks)
                {
                    int j;
                    int a = 0;
                    for (j = 0; j < 32768; j++)
                    {
                        if (sieve[32768 + j] != dconf->ss_sieve_n[(i + 1) * 32768 + j])
                        {
                            printf("med_sieve modified next block location %d from %u to %u\n",
                                j, dconf->ss_sieve_n[(i + 1) * 32768 + j], sieve[32768 + j]);
                            a = 1;
                        }
                    }
                    if (a == 1)
                        exit(1);
                }
            }
#endif

            lp_sieveblock_ptr(sieve, i, num_blocks, buckets, 0, dconf);

            // set the roots for the factors of a to force the following routine
            // to explicitly trial divide since we haven't found roots for them
            set_aprime_roots(sconf, invalid_root_marker, poly->qlisort, poly->s, fb_sieve_p, 0);
            scan_ptr(i, 0, sconf, dconf);
        }

        

#if defined( USE_SS_SEARCH )


        if (using_ss_search)
        {
            for (i = 0; i < num_blocks; i++)
            {
                if (i == minblock)
                    memset(dconf->ss_sieve_p + i * 32768, blockinit - 4, 32768);
                else
                    memset(dconf->ss_sieve_p + i * 32768, blockinit, 32768);
            }
        
            lp_sieve_ss(dconf->ss_sieve_p, 1, dconf);
        }
#endif

        

        for (i = 0; i < num_blocks; i++)
        {
#if defined( USE_SS_SEARCH )
            if (using_ss_search)
            {
                dconf->sieve = dconf->ss_sieve_p + i * 32768;
                sieve = dconf->sieve;
            }
#endif

            // set the roots for the factors of a such that
            // they will not be sieved.  we haven't found roots for them
            set_aprime_roots(sconf, invalid_root_marker, poly->qlisort, poly->s, fb_sieve_n, 1);
            if (i == minblock)
            {
#if !defined( USE_SS_SEARCH )
                memset(sieve, blockinit - 4, 32768);
#else
                if (!using_ss_search)
                {
                    memset(sieve, blockinit - 4, 32768);
                }
#endif
                med_sieve_ptr(sieve, fb_sieve_n, fb, start_prime, blockinit - 4);
            }
            else
            {
#if !defined( USE_SS_SEARCH )
                memset(sieve, blockinit, 32768);
#else
                if (!using_ss_search)
                {
                    memset(sieve, blockinit, 32768);
                }
#endif
                med_sieve_ptr(sieve, fb_sieve_n, fb, start_prime, blockinit);
            }
            lp_sieveblock_ptr(sieve, i, num_blocks, buckets, 1, dconf);

            // set the roots for the factors of a to force the following routine
            // to explicitly trial divide since we haven't found roots for them
            set_aprime_roots(sconf, invalid_root_marker, poly->qlisort, poly->s, fb_sieve_n, 0);
            scan_ptr(i, 1, sconf, dconf);
        }


        gettimeofday(&stop1, NULL);
        t_sieve_ss_buckets += ytools_difftime(&start1, &stop1);
        gettimeofday(&start1, NULL);


        // print a little more status info for huge jobs.
        if (sconf->digits_n > 110)
        {
            gettimeofday(&stop, NULL);
            t_time = ytools_difftime(&st, &stop);

            if (t_time > 5)
            {
				if (tdata->tindex == 0)
				{
					uint32_t dlp = sconf->dlp_useful + dconf->dlp_useful;
					uint32_t tlp = sconf->tlp_useful + dconf->tlp_useful;
					uint32_t att = sconf->attempted_cosiqs + dconf->attempted_cosiqs;

					// print some status
					gettimeofday(&stop, NULL);
					t_time = ytools_difftime(&start, &stop);

					if (sconf->use_dlp >= 2)
					{
						//printf("B: %u of %u, %u ppr, %u tpr, %u tlp attempts"
						//	" (%1.2f rels/sec)\n",
						//	dconf->numB, dconf->maxB, dlp, tlp, att,
						//	(double)dconf->buffered_rels / t_time);


                        if (sconf->do_batch)
                        {
                            char suffix;
                            float batched = dconf->rb.num_relations;

                            if (batched > 10000000)
                            {
                                batched /= 1000000.;
                                suffix = 'M';
                            }
                            else
                            {
                                batched /= 1000.;
                                suffix = 'k';
                            }

                            t_time = ytools_difftime(&sconf->totaltime_start, &stop);

                            //printf("thread: %u full, %u slp, "
                            //    "%u dlp, %u tlp (%1.1f%c batched), B = %u of %u\n",
                            //    dconf->num_full, dconf->num_slp,
                            //    dconf->dlp_useful, dconf->tlp_useful, (float)batched, suffix,
                            //    dconf->numB, dconf->maxB);

                            if (sconf->num_lp == 3)
                            {
                                printf("%u full, (%u,%u,%u) lp, "
                                    "(%1.1f%c raw-GCD), B = %u of %u \r",
                                    sconf->num_full + dconf->num_full, sconf->num_slp + dconf->num_slp,
                                    sconf->dlp_useful + dconf->dlp_useful,
                                    sconf->tlp_useful + dconf->tlp_useful,
                                    (float)batched, suffix, dconf->numB, dconf->maxB); 
                                
                            }
                            else if(sconf->num_lp == 4)
                            {
                                printf("%u full, (%u,%u,%u,%u) lp, "
                                    "(%1.1f%c raw-GCD), B = %u of %u \r",
                                    sconf->num_full + dconf->num_full, sconf->num_slp + dconf->num_slp,
                                    sconf->dlp_useful + dconf->dlp_useful,
                                    sconf->tlp_useful + dconf->tlp_useful,
                                    sconf->qlp_useful + dconf->qlp_useful,
                                    (float)batched, suffix, dconf->numB, dconf->maxB);
                            }

                            //printf("%u full, %u slp, "
                            //    "%u dlp, %u tlp (%1.1f%c batched), B = %u of %u\n",
                            //    sconf->num_full + dconf->num_full, sconf->num_slp + dconf->num_slp,
                            //    sconf->dlp_useful + dconf->dlp_useful, 
                            //    sconf->tlp_useful + dconf->tlp_useful,
                            //    batched, suffix,
                            //    dconf->numB, dconf->maxB);

                            //printf("last: %u, now: B = %u of %u, %u full, %u slp, "
                            //    "%u dlp, %u tlp (%1.1f%c batched), (%1.2f r/sec)\n",
                            //    sconf->last_numfull + sconf->last_numpartial, dconf->numB, dconf->maxB,
                            //    dconf->num_full, dconf->num_slp,
                            //    dconf->dlp_useful, dconf->tlp_useful, batched, suffix,
                            //    (double)(dconf->num_full + dconf->num_slp + dconf->dlp_useful +
                            //        dconf->tlp_useful) / t_time);
                        }
                        else
                        {
                            printf("last: %u, now: B = %u of %u, %u full, %u slp, "
                                "%u dlp (%ukatt, %ukprp), "
                                "%u tlp (%ukatt, %ukprp), (%1.2f r/sec)\n",
                                sconf->last_numfull + sconf->last_numpartial, dconf->numB, dconf->maxB,
                                dconf->num_full, dconf->num_slp,
                                dconf->dlp_useful, dconf->attempted_squfof / 1000,
                                dconf->dlp_prp / 1000, dconf->tlp_useful, dconf->attempted_cosiqs / 1000,
                                dconf->tlp_prp / 1000,
                                (double)(dconf->num_full + dconf->num_slp + dconf->dlp_useful +
                                    dconf->tlp_useful) / t_time);
                        }
					}
					else
					{
						printf("Bpoly %u of %u: buffered %u rels, checked %u (%1.2f rels/sec)\n",
							dconf->numB, dconf->maxB, dconf->buffered_rels, dconf->num,
							(double)dconf->buffered_rels / t_time);
					}

					// reset the timer
					gettimeofday(&st, NULL);
				}
            }
        }

        // next polynomial
        // use the stored Bl's and the gray code to find the next b
        nextB(dconf, sconf, 1);

#ifdef USE_BATCHPOLY

#if defined(USE_AVX512F)

        if (sconf->obj->HAS_AVX512F)
        {
            // every iteration we update the small-med prime's roots
            nextRoots_32k_avx2_small(sconf, dconf);

            // every N iterations we do the bucket sieve
            if ((dconf->numB % dconf->poly_batchsize) == 1)
            {
                nextRoots_32k_knl_polybatch(sconf, dconf);
            }
        }
        else
        {
            // every iteration we update the small-med prime's roots
            nextRoots_32k_generic_small(sconf, dconf);

            // every N iterations we do the bucket sieve
            if ((dconf->numB % dconf->poly_batchsize) == 1)
            {
                nextRoots_32k_generic_polybatch(sconf, dconf);
            }
        }

#else
        // every iteration we update the small-med prime's roots
        nextRoots_32k_generic_small(sconf, dconf);

        // every N iterations we do the bucket sieve
        if ((dconf->numB % dconf->poly_batchsize) == 1)
        {
            nextRoots_32k_generic_polybatch(sconf, dconf);
        }
#endif

#elif defined(USE_XLBATCHPOLY)

        if (sconf->obj->HAS_AVX512F)
        {
            // every iteration we update the small-med prime's roots
            nextRoots_32k_avx2_small(sconf, dconf);
            // and the med-big prime's roots
            nextRoots_32k_knl_bucket(sconf, dconf);

            // every N iterations we do the extra big prime's bucket sieve.
            if ((dconf->numB % dconf->poly_batchsize) == 1)
            {
                nextBigRoots_32k_knl_polybatch(sconf, dconf);
            }
        }

#else

#if defined(USE_AVX512F)

        if (sconf->obj->HAS_AVX512F)
        {
            //nextRoots_32k_avx2_intrin(sconf, dconf);
            nextRoots_32k_avx2_small(sconf, dconf);
            nextRoots_32k_knl_bucket(sconf, dconf);
        }
        else
        {
            nextRoots_ptr(sconf, dconf);
        }

#else
        // and update the roots
        nextRoots_ptr(sconf, dconf);
#endif

#endif

        gettimeofday(&stop1, NULL);
        t_update_buckets += ytools_difftime(&start1, &stop1);

        if (sconf->obj->THREADS >= 32)
        {
            // other threads may have got us past the threshold while we
            // are still working on this 'A' poly... check if we can stop.
            if ((sconf->num_found > sconf->num_needed) || (sconf->flag == 1))
            {
                break;
            }
        }
    }

    //exit(0);

done:

	gettimeofday (&stop, NULL);
    t_time = ytools_difftime(&start, &stop);

	dconf->rels_per_sec = (double)dconf->buffered_rels / t_time;

#if defined( USE_SS_SEARCH ) && defined( USE_POLY_BUCKET_SS )
    if ((using_ss_search) || (sconf->obj->VFLAG > 1))
    {
        printf("ss search total    : %1.4f\n", t_ss_search_total);
        printf("sieve/tdiv  buckets: %1.4f\n", t_sieve_ss_buckets);
        printf("poly update buckets: %1.4f\n", t_update_buckets);
    }
    if (using_ss_search)
    {
        ss_search_clear(sconf, dconf);
    }
#endif

	return 0;
}

uint32_t siqs_merge_data(dynamic_conf_t *dconf, static_conf_t *sconf)
{
	// the sconf structure holds the master list of relations and cycles.
	// merge everything we found in the last round of sieving into
	// this table and save relations out to disk
	uint32_t i;
	siqs_r *rel;
	char buf[1024];
    uint32_t curr_poly_idx = dconf->curr_poly->index;

	// save the current A value.
	if (!sconf->in_mem)
	{
        curr_poly_idx = dconf->curr_poly->index;
		gmp_sprintf(buf,"A 0x%Zx\n", dconf->curr_poly->mpz_poly_a);
		savefile_write_line(&sconf->obj->qs_obj.savefile,buf);
	}

	// save the data and merge into master cycle structure
	for (i = 0; i < dconf->buffered_rels; i++)
	{
		rel = dconf->relation_buf + i;


        // check to see if this relation's a-index is the same
        // as the current A.  If so, nothing extra to do.  If not,
        // we need to print the A associated with this relation.
        // This means that A values could get printed multiple times;
        // this is necessary when using batch factoring.  I think
        // filtering code will be ok with that.
        if ((sconf->do_batch) && (!sconf->in_mem))
        {
            if (rel->apoly_idx != curr_poly_idx)
            {
                curr_poly_idx = rel->apoly_idx;
                gmp_sprintf(buf, "A 0x%Zx\n", sconf->poly_a_list[rel->apoly_idx]);
                savefile_write_line(&sconf->obj->qs_obj.savefile, buf);
            }
        }

		save_relation_siqs(rel->sieve_offset, rel->large_prime,
			rel->num_factors, rel->fb_offsets, rel->poly_idx, 
			rel->apoly_idx, rel->parity, sconf);
	}

	// update some progress indicators
    sconf->total_blocks += dconf->total_blocks;
    sconf->total_reports += dconf->total_reports;
    sconf->total_surviving_reports += dconf->total_surviving_reports;
	sconf->num += dconf->num;
	sconf->tot_poly += dconf->tot_poly;
	sconf->failed_squfof += dconf->failed_squfof;
	sconf->attempted_squfof += dconf->attempted_squfof;
	sconf->num_full += dconf->num_full;
	sconf->num_slp += dconf->num_slp;
	sconf->dlp_outside_range += dconf->dlp_outside_range;
	sconf->dlp_prp += dconf->dlp_prp;
	sconf->dlp_useful += dconf->dlp_useful;
	sconf->failed_cosiqs += dconf->failed_cosiqs;
	sconf->attempted_cosiqs += dconf->attempted_cosiqs;
	sconf->tlp_outside_range += dconf->tlp_outside_range;
	sconf->tlp_prp += dconf->tlp_prp;
	sconf->tlp_useful += dconf->tlp_useful;
    sconf->qlp_outside_range += dconf->qlp_outside_range;
    sconf->qlp_prp += dconf->qlp_prp;
    sconf->qlp_useful += dconf->qlp_useful;
    sconf->lp_scan_failures += dconf->lp_scan_failures;
    sconf->num_scatter_opp += dconf->num_scatter_opp;
    sconf->num_scatter += dconf->num_scatter;

	// compute total relations found so far
	if (sconf->use_dlp < 2)
	{
		sconf->num_r = sconf->num_relations +
			sconf->num_cycles +
			sconf->components - sconf->vertices;
	}
	else if (sconf->do_batch)
	{
        relation_batch_t *rb = &dconf->rb;
        relation_batch_t *dest_rb;
        mpz_t x, tmp;
        mpz_init(x);
        mpz_init(tmp);
        mpz_set_ui(tmp, 0);
        uint32_t memuse = 0;

        dest_rb = &sconf->rb[sconf->batch_buffer_id - 1];
        memuse += dest_rb->num_relations_alloc * sizeof(cofactor_t);
        memuse += dest_rb->num_factors_alloc * sizeof(uint32_t);

        // merge in the batched relations
        for (i = 0; i < rb->num_relations; i++)
        {
            cofactor_t *c = rb->relations + i;
            uint32_t *f = rb->factors + c->factor_list_word;
            uint32_t *lp1 = f + c->num_factors_r + c->num_factors_a;
            int j;

            mpz_set_ui(x, lp1[c->lp_r_num_words - 1]);
            for (j = c->lp_r_num_words - 2; j >= 0; j--)
            {
                mpz_mul_2exp(x, x, 32);
                mpz_add_ui(x, x, lp1[j]);
            }

            relation_batch_add(c->a, c->b, c->signed_offset, 
                f, c->num_factors_r,
                x, NULL, 0, tmp, NULL, dest_rb);
        }

        if (sconf->obj->VFLAG > 2)
        {
            printf("merged %d relations from thread %d into batch buffer %d\n",
                rb->num_relations, dconf->tid, sconf->batch_buffer_id);
            
            memuse = dest_rb->num_relations_alloc * sizeof(cofactor_t) + 
                dest_rb->num_factors_alloc * sizeof(uint32_t) - memuse;
            
            if (memuse > 0)
            {
                printf("memory usage increased by %u bytes\n", memuse);
            }
        }

        // these relations are now in the static structure.
        // start collecting fresh ones in the thread.
        rb->num_relations = 0;
        rb->num_success = 0;
        rb->num_factors = 0;

        mpz_clear(x);
        mpz_clear(tmp);

		// when using TLP, the total number of relations found is
		// only updated when we do a filtering step.
		sconf->num_r = sconf->last_numfull + sconf->last_numpartial;
	}
    else
    {
        sconf->num_r = sconf->last_numfull + sconf->last_numpartial;
    }

	return sconf->num_r;
}

int siqs_check_restart(dynamic_conf_t *dconf, static_conf_t *sconf)
{
	fact_obj_t *obj = sconf->obj;
	char buf[1024];
	int state = 0;
    sconf->is_restart = 0;

	// if we want to do an in-memory factorization, then 
	// ignore the current state of the savefile and don't
	// prepare it for use.
	if (sconf->in_mem)
		return 0;

	// we're now almost ready to start, but first
	// check if this number has had work done
	restart_siqs(sconf,dconf);
	
	if ((uint32_t)sconf->num_r >= (sconf->factor_base->B + sconf->num_extra_relations)) 
	{
		// we've got enough total relations to stop		
		savefile_open(&obj->qs_obj.savefile,SAVEFILE_APPEND);	
		dconf->buckets->list = NULL;
		// signal that we should proceed to post-processing
		state = 1;
	}
	else if (sconf->num_r > 0)
	{     
		// we've got some relations, but not enough to finish.
		// whether or not this is a big job, it needed to be resumed
		// once so treat it as if it will need to be again.  use the savefile.
		savefile_open(&obj->qs_obj.savefile,SAVEFILE_APPEND);

        // don't try to do optimization of the small cutoff... it is not 
        // designed to cope with starting with a bunch of relations and poly_a's
        obj->qs_obj.no_small_cutoff_opt = 1;
        sconf->is_restart = 1;
	}
	else
	{
		// no relations found, get ready for new factorization
		// we'll be writing to the savefile as we go, so get it ready
		savefile_open(&obj->qs_obj.savefile,SAVEFILE_WRITE);
		gmp_sprintf(buf,"N 0x%Zx\n", sconf->obj->qs_obj.gmp_n);
		savefile_write_line(&obj->qs_obj.savefile,buf);
		savefile_flush(&obj->qs_obj.savefile);
		savefile_close(&obj->qs_obj.savefile);
		// and get ready for collecting relations
		savefile_open(&obj->qs_obj.savefile,SAVEFILE_APPEND);
	}

	return state;
}

void print_siqs_splash(dynamic_conf_t *dconf, static_conf_t *sconf)
{
    //print some info to the screen and the log file
    char inst_set[16];

#if defined (USE_AVX512F)
    if (sconf->obj->HAS_AVX512BW) strcpy(inst_set, "AVX512BW");
    else if (sconf->obj->HAS_AVX512F) strcpy(inst_set, "AVX512F");
    else if (sconf->obj->HAS_AVX2) strcpy(inst_set, "AVX2");
    else if (sconf->obj->HAS_SSE41) strcpy(inst_set, "SSE41");
    else strcpy(inst_set, "SSE2");
#elif defined(USE_AVX2)
    if (sconf->obj->HAS_AVX2) strcpy(inst_set, "AVX2");
    else if (sconf->obj->HAS_SSE41) strcpy(inst_set, "SSE41");
    else strcpy(inst_set, "SSE2");
#elif defined(USE_SSE41)
    if (sconf->obj->HAS_SSE41) strcpy(inst_set, "SSE41");
    else strcpy(inst_set, "SSE2");
#else
    strcpy(inst_set, "SSE2");
#endif

    if (sconf->obj->VFLAG > 0)
    {
        printf("\n==== sieve params ====\n");
        printf("n = %d digits, %d bits\n", sconf->digits_n, sconf->bits);
        printf("factor base: %d primes (max prime = %u)\n", sconf->factor_base->B, sconf->pmax);

        if (sconf->obj->qs_obj.gbl_override_lpb > 0)
        {
            printf("single large prime cutoff: %u (2^%d)\n",
                sconf->large_prime_max, sconf->obj->qs_obj.gbl_override_lpb);
        }
        else if (sconf->large_mult > 1000)
        {
            int bits = 0;
            uint32_t lpb = sconf->large_mult;
            while (lpb > 1)
            {
                lpb >>= 1;
                bits++;
            }

            if (sconf->large_prime_max == (1 << bits))
            {
                printf("single large prime cutoff: %u (2^%d)\n",
                    sconf->large_prime_max, bits);
            }
            else
            {
                printf("single large prime cutoff: %u\n",
                    sconf->large_prime_max);
            }
        }
        else
        {
            printf("single large prime cutoff: %u (%d * pmax)\n",
                sconf->large_prime_max, sconf->large_mult);
        }

        if (sconf->use_dlp >= 1)
        {
            printf("double large prime range from %" PRIu64 " to %" PRIu64 "\n",
                sconf->max_fb2, sconf->large_prime_max2);
            printf("DLP MFB = %1.2f\n", sconf->dlp_exp);

            if (sconf->use_dlp >= 2)
            {
                printf("triple large prime range from %1.0f to %1.0f\n",
                    sconf->max_fb3, sconf->large_prime_max3);
                printf("TLP Batch Div = %1.2f, TLP MFB = %1.2f\n",
                    sconf->obj->qs_obj.gbl_override_bdiv, sconf->tlp_exp);
            }
            if (sconf->use_dlp >= 3)
            {
                printf("quad large prime range from %1.0f to %1.0f\n",
                    sconf->max_fb4, sconf->large_prime_max4);
                printf("QLP Batch Div = %1.2f, QLP MFB = %1.2f\n",
                    sconf->obj->qs_obj.gbl_override_bdiv, sconf->qlp_exp);
            }
        }
        if (dconf->buckets->list != NULL)
        {
            printf("allocating %d large prime slices of factor base\n",
                dconf->buckets->alloc_slices);
            printf("buckets hold %d elements\n", BUCKET_ALLOC);
#ifdef USE_BATCHPOLY
            printf("processing polynomials in batches of %d\n", dconf->poly_batchsize);
#endif
#ifdef USE_BATCHPOLY_X2
            printf("processing polynomials in batches of %d for %u XL primes\n",
                dconf->poly_batchsize, sconf->factor_base->B - sconf->factor_base->x2_large_B);
#endif
            printf("large prime hashtables have %d bytes\n",
                dconf->buckets->list_size * BUCKET_ALLOC * sizeof(uint32_t));
        }
        printf("using %s enabled 32k sieve core\n", inst_set);
        printf("sieve interval: %d blocks of size %d\n",
            sconf->num_blocks, sconf->qs_blocksize);
        printf("polynomial A has ~ %d factors\n",
            (int)mpz_sizeinbase(sconf->target_a, 2) / 11);

        if (sconf->knmod8 == 1)
        {
            printf("using multiplier of %u\n", sconf->multiplier);
            printf("using Q2(x) polynomials for kN mod 8 = %u\n", sconf->knmod8);
        }
        else
        {
            printf("using multiplier of %u (kn mod 8 == %u)\n", sconf->multiplier, sconf->knmod8);
        }

        if (sconf->in_mem)
        {
            printf("using in-memory storage, aborts will lose all progress\n");
        }

#ifdef USE_SS_SEARCH

        if ((sconf->factor_base->ss_start_B > 0) &&
            (sconf->factor_base->ss_start_B < sconf->factor_base->B))
        {
#if defined(USE_LINKED_LIST_SS)
            printf("using linked-list subsum-search at fb index > %d\n", sconf->factor_base->ss_start_B);
#elif defined(USE_SORTED_LIST_SS)
            printf("using sorted-list subsum-search at fb index > %d\n", sconf->factor_base->ss_start_B);
#elif defined(USE_POLY_BUCKET_SS)
            printf("using poly-bucket subsum-search at fb index > %d (prime %u), %d slices\n",
                sconf->factor_base->ss_start_B,
                sconf->factor_base->list->prime[sconf->factor_base->ss_start_B],
                sconf->factor_base->num_ss_slices);
#elif defined(USE_DIRECT_SIEVE_SS)
            printf("using direct-sieve subsum-search at fb index > %d (prime %u), %d slices\n",
                sconf->factor_base->ss_start_B,
                sconf->factor_base->list->prime[sconf->factor_base->ss_start_B],
                sconf->factor_base->num_ss_slices);
#endif
        }
#endif
        

        printf("using SPV correction of %d bits, starting at offset %d\n",
            sconf->tf_small_cutoff, sconf->sieve_small_fb_start);
		printf("trial factoring cutoff at %d bits\n",sconf->tf_closnuf);
	}

	if (sconf->obj->VFLAG >= 0)
	{
		if (sconf->obj->THREADS == 1)
		{
			printf("\n==== sieving in progress (1 thread): %7u relations needed ====\n",
				sconf->factor_base->B + sconf->num_extra_relations);
			printf(  "====           Press ctrl-c to abort and save state           ====\n");
		}
		else
		{
			printf("\n==== sieving in progress (%3d threads): %7u relations needed ====\n",
                sconf->obj->THREADS,sconf->factor_base->B + sconf->num_extra_relations);
			printf(  "====             Press ctrl-c to abort and save state            ====\n");
		}
	}

	if (sconf->obj->logfile != NULL)
	{
		logprint(sconf->obj->logfile,"==== sieve params ====\n");
		logprint(sconf->obj->logfile,"n = %d digits, %d bits\n",
			sconf->digits_n,sconf->bits);
		logprint(sconf->obj->logfile,"factor base: %d primes (max prime = %u)\n",
			sconf->factor_base->B,sconf->pmax);
		if (sconf->obj->qs_obj.gbl_override_lpb > 0)
		{
			logprint(sconf->obj->logfile, "single large prime cutoff: %u (2^%d)\n",
				sconf->large_prime_max, sconf->obj->qs_obj.gbl_override_lpb);
		}
        else if (sconf->large_mult > 1000)
        {
            int bits = 0;
            uint32_t lpb = sconf->large_mult;
            while (lpb > 1)
            {
                lpb >>= 1;
                bits++;
            }

            if (sconf->large_prime_max == (1 << bits))
            {
                logprint(sconf->obj->logfile, "single large prime cutoff: %u (2^%d)\n",
                    sconf->large_prime_max, bits);
            }
            else
            {
                logprint(sconf->obj->logfile, "single large prime cutoff: %u\n",
                    sconf->large_prime_max);
            }
        }
        else
        {
            logprint(sconf->obj->logfile, "single large prime cutoff: %u (%d * pmax)\n",
                sconf->large_prime_max, sconf->large_mult);
        }

		if (sconf->use_dlp >= 1)
		{
			logprint(sconf->obj->logfile, "double large prime range from %" PRIu64 " to %" PRIu64 "\n",
				sconf->max_fb2, sconf->large_prime_max2);
            logprint(sconf->obj->logfile, "DLP MFB = %1.2f\n", sconf->dlp_exp);

			if (sconf->use_dlp >= 2)
			{
				logprint(sconf->obj->logfile, "triple large prime range from %1.0f to %1.0f\n",
					sconf->max_fb3, sconf->large_prime_max3);
                logprint(sconf->obj->logfile, "TLP Batch Div = %1.2f, TLP MFB = %1.2f\n",
                    sconf->obj->qs_obj.gbl_override_bdiv, sconf->tlp_exp);
			}
            if (sconf->use_dlp >= 4)
            {
                logprint(sconf->obj->logfile, "quad large prime range from %1.0f to %1.0f\n",
                    sconf->max_fb4, sconf->large_prime_max4);
                logprint(sconf->obj->logfile, "QLP Batch Div = %1.2f, QLP MFB = %1.2f\n",
                    sconf->obj->qs_obj.gbl_override_bdiv, sconf->qlp_exp);
            }
		}
		if (dconf->buckets->list != NULL)
		{
			logprint(sconf->obj->logfile,"allocating %d large prime slices of factor base\n",
				dconf->buckets->alloc_slices);
			logprint(sconf->obj->logfile,"buckets hold %d elements\n",BUCKET_ALLOC);
#ifdef USE_BATCHPOLY
            logprint(sconf->obj->logfile, "processing polynomials in batches of %d\n", 
                dconf->poly_batchsize);
#endif
#ifdef USE_BATCHPOLY_X2
            logprint(sconf->obj->logfile, "processing polynomials in batches of %d for XL primes\n",
                dconf->poly_batchsize);
#endif
            logprint(sconf->obj->logfile, "large prime hashtables have %d bytes\n", 
                dconf->buckets->list_size * BUCKET_ALLOC * sizeof(uint32_t));
		}
        logprint(sconf->obj->logfile,"using %s enabled 32k sieve core\n", inst_set);
		logprint(sconf->obj->logfile,"sieve interval: %d blocks of size %d\n",
			sconf->num_blocks,sconf->qs_blocksize);
		logprint(sconf->obj->logfile,"polynomial A has ~ %d factors\n", 
			mpz_sizeinbase(sconf->target_a, 2) / 11);
		logprint(sconf->obj->logfile,"using multiplier of %u\n",sconf->multiplier);

        if (sconf->knmod8 == 1)
        {
            logprint(sconf->obj->logfile, "using multiplier of %u\n", sconf->multiplier);
            logprint(sconf->obj->logfile, "using Q2(x) polynomials for kN mod 8 = %u\n", sconf->knmod8);
        }
        else
        {
            logprint(sconf->obj->logfile, "using multiplier of %u (kn mod 8 == %u)\n", sconf->multiplier, sconf->knmod8);
        }

		logprint(sconf->obj->logfile,"using SPV correction of %d bits, starting at offset %d\n",
			sconf->tf_small_cutoff,sconf->sieve_small_fb_start);
		logprint(sconf->obj->logfile,"trial factoring cutoff at %d bits\n",
			sconf->tf_closnuf);
		if (sconf->obj->THREADS == 1)
			logprint(sconf->obj->logfile,"==== sieving started (1 thread) ====\n");
		else
			logprint(sconf->obj->logfile,"==== sieving started (%2d threads) ====\n",
                sconf->obj->THREADS);
	}

	return;
}

int siqs_dynamic_init(dynamic_conf_t *dconf, static_conf_t *sconf)
{
    //allocate the dynamic structure which hold scratch space for 
    //various things used during sieving.
    uint32_t i, memsize;

    if (sconf->obj->VFLAG > 2)
    {
        printf("memory usage during sieving:\n");
    }

    //workspace bigints
    mpz_init(dconf->gmptmp1);
    mpz_init(dconf->gmptmp2);
    mpz_init(dconf->gmptmp3);

    //this stuff changes with every new poly
    //allocate a polynomial structure which will hold the current
    //set of polynomial coefficients (a,b,c) and other info
    dconf->curr_poly = (siqs_poly *)malloc(sizeof(siqs_poly));
    mpz_init(dconf->curr_poly->mpz_poly_a);
    mpz_init(dconf->curr_poly->mpz_poly_b);
    mpz_init(dconf->curr_poly->mpz_poly_c);
    dconf->curr_poly->qlisort = (int *)malloc(MAX_A_FACTORS * sizeof(int));
    dconf->curr_poly->gray = (char *)malloc((1 << (MAX_A_FACTORS - 1)) * sizeof(char));
    dconf->curr_poly->nu = (char *)malloc((1 << (MAX_A_FACTORS - 1)) * sizeof(char));
    if (sconf->obj->VFLAG > 2)
    {
        memsize = MAX_A_FACTORS * sizeof(int);
        memsize += (1 << (MAX_A_FACTORS - 1)) * sizeof(char);
        memsize += (1 << (MAX_A_FACTORS - 1)) * sizeof(char);
        printf("\tcurr_poly structure: %d bytes\n", memsize);
    }

    //initialize temporary storage of relations
    dconf->relation_buf = (siqs_r *)malloc(32768 * sizeof(siqs_r));
    dconf->buffered_rel_alloc = 32768;
    dconf->buffered_rels = 0;

    if (sconf->obj->VFLAG > 2)
    {
        memsize = 32768 * sizeof(siqs_r);
        printf("\trelation buffer: %d bytes\n", memsize);
    }

    //allocate the sieving factor bases

    dconf->comp_sieve_p = (sieve_fb_compressed *)malloc(sizeof(sieve_fb_compressed));
    dconf->comp_sieve_n = (sieve_fb_compressed *)malloc(sizeof(sieve_fb_compressed));

    dconf->comp_sieve_p->prime = (uint16_t *)xmalloc_align(
        (size_t)(sconf->factor_base->med_B * sizeof(uint16_t)));
    dconf->comp_sieve_p->root1 = (uint16_t *)xmalloc_align(
        (size_t)(sconf->factor_base->med_B * sizeof(uint16_t)));
    dconf->comp_sieve_p->root2 = (uint16_t *)xmalloc_align(
        (size_t)(sconf->factor_base->med_B * sizeof(uint16_t)));
    dconf->comp_sieve_p->logp = (uint16_t *)xmalloc_align(
        (size_t)(sconf->factor_base->med_B * sizeof(uint16_t)));

    dconf->comp_sieve_n->prime = (uint16_t *)xmalloc_align(
        (size_t)(sconf->factor_base->med_B * sizeof(uint16_t)));
    dconf->comp_sieve_n->root1 = (uint16_t *)xmalloc_align(
        (size_t)(sconf->factor_base->med_B * sizeof(uint16_t)));
    dconf->comp_sieve_n->root2 = (uint16_t *)xmalloc_align(
        (size_t)(sconf->factor_base->med_B * sizeof(uint16_t)));
    dconf->comp_sieve_n->logp = (uint16_t *)xmalloc_align(
        (size_t)(sconf->factor_base->med_B * sizeof(uint16_t)));

    dconf->update_data.sm_firstroots1 = (uint16_t *)xmalloc_align(
        (size_t)(sconf->factor_base->med_B * sizeof(uint16_t)));
    dconf->update_data.sm_firstroots2 = (uint16_t *)xmalloc_align(
        (size_t)(sconf->factor_base->med_B * sizeof(uint16_t)));
    dconf->update_data.firstroots1 = (int *)xmalloc_align(
        (size_t)(sconf->factor_base->B * sizeof(int)));
    dconf->update_data.firstroots2 = (int *)xmalloc_align(
        (size_t)(sconf->factor_base->B * sizeof(int)));
    dconf->update_data.prime = (uint32_t *)xmalloc_align(
        (size_t)(sconf->factor_base->B * sizeof(uint32_t)));
    dconf->update_data.logp = (uint8_t *)xmalloc_align(
        (size_t)(sconf->factor_base->B * sizeof(uint8_t)));
    dconf->rootupdates = (int *)xmalloc_align(
        (size_t)(MAX_A_FACTORS * sconf->factor_base->B * sizeof(int)));
    dconf->sm_rootupdates = (uint16_t *)xmalloc_align(
        (size_t)(MAX_A_FACTORS * sconf->factor_base->med_B * sizeof(uint16_t)));


    if (sconf->obj->VFLAG > 2)
    {
        memsize = sconf->factor_base->med_B * sizeof(sieve_fb_compressed) * 2;
        printf("\tfactor bases: %d bytes\n", memsize);
    }

    if (sconf->obj->VFLAG > 2)
    {
        memsize = sconf->factor_base->B * sizeof(int) * 2;
        memsize += sconf->factor_base->B * sizeof(uint32_t);
        memsize += sconf->factor_base->B * sizeof(uint8_t);
        printf("\tupdate data: %d bytes\n", memsize);
    }


#ifdef USE_SS_SEARCH
    

    if ((sconf->factor_base->ss_start_B > 0) &&
        (sconf->factor_base->ss_start_B < sconf->factor_base->B))
    {
        dconf->using_ss_search = 1;

#if defined(USE_SS_SEARCH) && defined( USE_DIRECT_SIEVE_SS )

        // this should be scaled such that the total sieve 
        // area fits in L3 cache.
        dconf->ss_sieve_sz = 4096;

        // 65536 is max number of polys
        dconf->ss_sieve_p = (uint8_t*)xmalloc_align(2048 * 
            dconf->ss_sieve_sz * sizeof(uint8_t));
        dconf->ss_sieve_n = (uint8_t*)xmalloc_align(2048 *
            dconf->ss_sieve_sz * sizeof(uint8_t));

#else
        dconf->ss_sieve_p = (uint8_t*)xmalloc_align(4 * sconf->num_blocks *
            sconf->qs_blocksize * sizeof(uint8_t));
        dconf->ss_sieve_n = (uint8_t*)xmalloc_align(4 * sconf->num_blocks *
            sconf->qs_blocksize * sizeof(uint8_t));
#endif
    }
    else
    {
        dconf->using_ss_search = 0;

        dconf->ss_sieve_n = NULL;
        dconf->ss_sieve_p = NULL;

        dconf->sieve = (uint8_t*)xmalloc_align(
            (size_t)(sconf->qs_blocksize * 4 * sizeof(uint8_t)));
        dconf->sieve = &dconf->sieve[2 * sconf->qs_blocksize];
    }
#else
    // allocate the sieve with a guard region.  this lets
    // us be a little more sloppy with vector-based sieving
    // and resieving.
    dconf->sieve = (uint8_t*)xmalloc_align(
        (size_t)(sconf->qs_blocksize * 4 * sizeof(uint8_t)));
    dconf->sieve = &dconf->sieve[2 * sconf->qs_blocksize];
#endif

    if (sconf->obj->VFLAG > 2)
    {
        memsize = sconf->qs_blocksize * sizeof(uint8_t);
        printf("\tsieve: %d bytes\n", memsize);
    }

    //allocate the Bl array, space for MAX_Bl bigint numbers
    dconf->Bl = (mpz_t *)malloc(MAX_A_FACTORS * sizeof(mpz_t));
    for (i = 0; i < MAX_A_FACTORS; i++)
    {
        mpz_init(dconf->Bl[i]);
    }

    //copy the unchanging part to the sieving factor bases
    for (i = 2; i < sconf->factor_base->med_B; i++)
    {
        uint32_t p = sconf->factor_base->list->prime[i];
        uint32_t lp = sconf->factor_base->list->logprime[i];

        dconf->comp_sieve_p->logp[i] = (uint8_t)lp;
        dconf->comp_sieve_p->prime[i] = (uint16_t)p;
        dconf->comp_sieve_n->logp[i] = (uint8_t)lp;
        dconf->comp_sieve_n->prime[i] = (uint16_t)p;

        dconf->update_data.prime[i] = p;
        dconf->update_data.logp[i] = lp;
    }

    for (; i < sconf->factor_base->B; i++)
    {
        dconf->update_data.prime[i] = sconf->factor_base->list->prime[i];
        dconf->update_data.logp[i] = sconf->factor_base->list->logprime[i];
    }

    //check if we should use bucket sieving, and allocate structures if so
    if (sconf->factor_base->B > sconf->factor_base->med_B)
    {
        dconf->buckets = (lp_bucket *)malloc(sizeof(lp_bucket));

        //test to see how many slices we'll need.
        testRoots_ptr(sconf, dconf);

#if defined(USE_BATCHPOLY) || defined(USE_BATCHPOLY_X2)
        // this should be a function of the L2 size.
//        dconf->poly_batchsize = MAX(2, L2CACHE / 2 / fobj->THREADS /
//            (2 * sconf->num_blocks * dconf->buckets->alloc_slices *
//                BUCKET_ALLOC * sizeof(uint32_t)));

        dconf->poly_batchsize = 4;

#else
        dconf->poly_batchsize = 1;
#endif


#if defined(USE_SS_SEARCH) && defined(USE_POLY_BUCKET_SS)

        if ((sconf->factor_base->ss_start_B > 0) &&
            (sconf->factor_base->ss_start_B < sconf->factor_base->B))
        {
            int numslices = sconf->factor_base->num_ss_slices;
            int slicesz = sconf->factor_base->slice_size;

            dconf->firstroot1a = (int*)xmalloc(sconf->factor_base->B * sizeof(int));
            dconf->firstroot1b = (int*)xmalloc(sconf->factor_base->B * sizeof(int));
            dconf->firstroot2 = (int*)xmalloc(sconf->factor_base->B * sizeof(int));

            dconf->polymap = (int*)xmalloc((1 << MAX_A_FACTORS) * sizeof(int));
            dconf->polynums = (int*)xmalloc((1 << MAX_A_FACTORS) * sizeof(int));
            dconf->polyv = (int*)xmalloc((1 << MAX_A_FACTORS) * sizeof(int));
            dconf->polysign = (int*)xmalloc((1 << MAX_A_FACTORS) * sizeof(int));

            mpz_init(dconf->polyb1);
            mpz_init(dconf->polyb2);

            dconf->ss_slices_p = (ss_bucket_slice_t*)xmalloc(numslices * sizeof(ss_bucket_slice_t));
            dconf->ss_slices_n = (ss_bucket_slice_t*)xmalloc(numslices * sizeof(ss_bucket_slice_t));

            dconf->num_ss_slices = numslices;

            // the size of a poly bucket, in bits. This choice
            // will depend on the slice and sieve-region size choices
            // made at compile time.  A poly-bucket must be
            // large enough to store all sieve-hits for all primes treated
            // by subset-sum.  The probability of a sieve-hit for a prime
            // is sieve-region-size / p, per root.  So as more primes are 
            // treated by subset-sum, and with larger allowed sieve-region
            // size, the buckets will have to be correspondingly larger.
            // Making this as small as possible is best, because we need one
            // bucket per b-poly so it multiplies up quickly.
            if (sconf->obj->qs_obj.gbl_override_ssalloc_flag > 0)
            {
                if (sconf->obj->qs_obj.gbl_override_ssalloc >= 32)
                {
                    // invalid entry, use default
                    printf("Invalid subset-sum bucket alloc, should be a number of bits < 16\n");
                    printf("Using default value of 12 bits per bucket\n");
                    dconf->poly_buckets_allocated = 12;
                }
                else
                {
                    dconf->poly_buckets_allocated = sconf->obj->qs_obj.gbl_override_ssalloc;
                }
            }
            else
            {
                dconf->poly_buckets_allocated = 12;
            }
            

            for (i = 0; i < numslices; i++)
            {
                // larger slices means more primes per slice which means the buckets
                // need to be bigger per slice.  also of course, larger slices
                // means fewer maximum blocks allowed, because we need more bits to
                // store prime id's and so we have fewer bits for storing root locations.
                dconf->ss_slices_p[i].alloc = (1 << dconf->poly_buckets_allocated);
                dconf->ss_slices_p[i].numbuckets = 65536;

#ifndef USE_POLY_BUCKET_PN_COMBINED_VARIATION
                dconf->ss_slices_n[i].alloc = 2048;
                dconf->ss_slices_n[i].numbuckets = 131072;
#endif

                int a = dconf->ss_slices_p[i].alloc;
                int nb = dconf->ss_slices_p[i].numbuckets;

                dconf->ss_slices_p[i].elements = (uint32_t*)xmalloc(a * nb * sizeof(uint32_t));
                dconf->ss_slices_p[i].size = (uint32_t*)xmalloc(nb * sizeof(uint32_t));

#ifndef USE_POLY_BUCKET_PN_COMBINED_VARIATION
                dconf->ss_slices_n[i].elements = (uint32_t*)xmalloc(a * nb * sizeof(uint32_t));
                dconf->ss_slices_n[i].size = (uint32_t*)xmalloc(nb * sizeof(uint32_t));
#endif

                dconf->ss_slices_p[i].fboffset = i * slicesz + sconf->factor_base->ss_start_B;
                dconf->ss_slices_n[i].fboffset = i * slicesz + sconf->factor_base->ss_start_B;

                if ((i * slicesz + sconf->factor_base->ss_start_B) < sconf->factor_base->B)
                {
                    dconf->ss_slices_p[i].logp =
                        sconf->factor_base->list->logprime[i * slicesz + sconf->factor_base->ss_start_B];
                    dconf->ss_slices_n[i].logp =
                        sconf->factor_base->list->logprime[i * slicesz + sconf->factor_base->ss_start_B];
                }
                else
                {
                    dconf->ss_slices_p[i].logp =
                        sconf->factor_base->list->logprime[sconf->factor_base->B - 1];
                    dconf->ss_slices_n[i].logp =
                        sconf->factor_base->list->logprime[sconf->factor_base->B - 1];
                }
            }
        }

#elif defined(USE_SS_SEARCH) && defined( USE_DIRECT_SIEVE_SS )

        dconf->firstroot1a = (int*)xmalloc(sconf->factor_base->B * sizeof(int));
        dconf->firstroot1b = (int*)xmalloc(sconf->factor_base->B * sizeof(int));
        dconf->firstroot2 = (int*)xmalloc(sconf->factor_base->B * sizeof(int));
        dconf->report_ht_p = (uint16_t*)xcalloc((1<<24), sizeof(uint16_t));
        dconf->report_ht_n = (uint16_t*)xcalloc((1<<24), sizeof(uint16_t));
        dconf->report_ht_size = (1 << 24);

        dconf->polymap = (int*)xmalloc((1 << MAX_A_FACTORS) * sizeof(int));
        dconf->polynums = (int*)xmalloc((1 << MAX_A_FACTORS) * sizeof(int));
        dconf->polyv = (int*)xmalloc((1 << MAX_A_FACTORS) * sizeof(int));
        dconf->polysign = (int*)xmalloc((1 << MAX_A_FACTORS) * sizeof(int));

        mpz_init(dconf->polyb1);
        mpz_init(dconf->polyb2);
#endif

        //initialize the bucket lists and auxilary info.
        //printf("allocating space for buckets with %u blocks, %u slices, and %u polys\n",
        //    sconf->num_blocks, dconf->buckets->alloc_slices, dconf->poly_batchsize);

        // two entries per block (p and n sides), per slice, and per poly.
        dconf->buckets->num = (uint32_t *)xmalloc_align(dconf->poly_batchsize *
            2 * sconf->num_blocks * dconf->buckets->alloc_slices * sizeof(uint32_t));
        
        // one entry per slice to record the starting prime and logp indices per slice.
        // experimental variation: use one entry per slice and per poly, so that we
        // can sort xlarge primes in a batch of polys and other primes one poly at a time.
        dconf->buckets->fb_bounds = (uint32_t *)malloc(
            dconf->buckets->alloc_slices * dconf->poly_batchsize * sizeof(uint32_t));
        dconf->buckets->logp = (uint8_t *)calloc(
            dconf->buckets->alloc_slices * dconf->poly_batchsize, sizeof(uint8_t));

        // the total number of buckets (size of buckets->num)
        dconf->buckets->list_size = dconf->poly_batchsize * 2 *
            sconf->num_blocks * dconf->buckets->alloc_slices;

        dconf->buckets->num_slices_batch = (uint32_t*)xmalloc(dconf->poly_batchsize * sizeof(uint32_t));

        // now allocate space for the bucket entries.  Each entry is a 32-bit integer
        // split into 16 bits for prime id and 16 bits for block location.
        dconf->buckets->list = (uint32_t *)xmalloc_align(
            dconf->buckets->list_size *
            BUCKET_ALLOC * sizeof(uint32_t));
    }
    else
    {
        dconf->poly_batchsize = 1;
        dconf->buckets = (lp_bucket *)malloc(sizeof(lp_bucket));
        dconf->buckets->list = NULL;
        dconf->buckets->list_size = 0;
        dconf->buckets->alloc_slices = 0;
        dconf->buckets->num_slices = 0;
    }

    if (sconf->obj->VFLAG > 2)
    {
        memsize = dconf->buckets->list_size * sizeof(uint32_t);
        memsize += dconf->buckets->alloc_slices * sizeof(uint32_t);
        memsize += dconf->buckets->alloc_slices * sizeof(uint8_t);
        memsize += dconf->buckets->list_size * BUCKET_ALLOC * sizeof(uint32_t);
        printf("\tbucket data: %d bytes\n", memsize);
    }

#ifdef USE_XLBUCKET
    // allow room for every xl-prime to hit the interval once.
    i = sconf->factor_base->B - sconf->factor_base->x2_large_B;

    dconf->xl_nbucket.list = (uint32_t*)xmalloc_align(i * sizeof(uint32_t));
    dconf->xl_pbucket.list = (uint32_t*)xmalloc_align(i * sizeof(uint32_t));
    
    i = ((sconf->factor_base->B - sconf->factor_base->x2_large_B) / 2048) + 2;

    dconf->xl_nbucket.sliceid = (uint32_t*)xmalloc_align(i * sizeof(uint32_t));
    dconf->xl_pbucket.sliceid = (uint32_t*)xmalloc_align(i * sizeof(uint32_t));
    dconf->xl_nbucket.slicenum = (uint32_t*)xmalloc_align(i * sizeof(uint32_t));
    dconf->xl_pbucket.slicenum = (uint32_t*)xmalloc_align(i * sizeof(uint32_t));
    dconf->xl_nbucket.slicelogp = (uint8_t*)xmalloc_align(i * sizeof(uint32_t));
    dconf->xl_pbucket.slicelogp = (uint8_t*)xmalloc_align(i * sizeof(uint32_t));

    dconf->xl_pbucket.alloc_slices = i;
    dconf->xl_pbucket.alloc_slices = i;

    if (sconf->obj->VFLAG > 1)
    {
        printf("\txlbucket data: %d bytes\n", 2 * sizeof(uint32_t) *
            (sconf->factor_base->B - sconf->factor_base->x2_large_B));
    }
#endif

    //used in trial division to mask out the fb_index portion of bucket entries, so that
    //multiple block locations can be searched for in parallel using SSE2 instructions
    dconf->mask = (uint16_t *)xmalloc_align(16 * sizeof(uint16_t));
    dconf->mask2 = (uint32_t *)xmalloc_align(8 * sizeof(uint32_t));

    dconf->mask[1] = 0xFFFF;
    dconf->mask[3] = 0xFFFF;
    dconf->mask[5] = 0xFFFF;
    dconf->mask[7] = 0xFFFF;
    dconf->mask[9] = 0xFFFF;
    dconf->mask[11] = 0xFFFF;
    dconf->mask[13] = 0xFFFF;
    dconf->mask[15] = 0xFFFF;

    dconf->mask2[0] = 0x0000ffff;
    dconf->mask2[1] = 0x0000ffff;
    dconf->mask2[2] = 0x0000ffff;
    dconf->mask2[3] = 0x0000ffff;
    dconf->mask2[4] = 0x0000ffff;
    dconf->mask2[5] = 0x0000ffff;
    dconf->mask2[6] = 0x0000ffff;
    dconf->mask2[7] = 0x0000ffff;

    // used in SIMD optimized resiever
    dconf->corrections = (uint16_t *)xmalloc_align(32 * sizeof(uint16_t));

    // array of sieve locations scanned from the sieve block that we
    // will submit to trial division.  make it the size of a sieve block 
    // in the pathological case that every sieve location is a report
    dconf->reports = (uint32_t *)malloc(MAX_SIEVE_REPORTS * sizeof(uint32_t));
    dconf->num_reports = 0;
    dconf->Qvals = (mpz_t *)malloc(MAX_SIEVE_REPORTS * sizeof(mpz_t));
    for (i = 0; i < MAX_SIEVE_REPORTS; i++)
    {
        mpz_init(dconf->Qvals[i]); //, 2*sconf->bits);
        mpz_set_ui(dconf->Qvals[i], 0);
    }

    dconf->valid_Qs = (int *)malloc(MAX_SIEVE_REPORTS * sizeof(int));
    dconf->smooth_num = (int *)malloc(MAX_SIEVE_REPORTS * sizeof(int));
    dconf->failed_squfof = 0;
    dconf->attempted_squfof = 0;
    dconf->dlp_outside_range = 0;
    dconf->dlp_prp = 0;
    dconf->num_slp = 0;
    dconf->num_full = 0;
    dconf->dlp_useful = 0;
    dconf->failed_cosiqs = 0;
    dconf->attempted_cosiqs = 0;
    dconf->tlp_outside_range = 0;
    dconf->tlp_prp = 0;
    dconf->tlp_useful = 0;
    dconf->qlp_outside_range = 0;
    dconf->qlp_prp = 0;
    dconf->qlp_useful = 0;

    dconf->mdata = monty_alloc();

    // ?? can't remember what this was for but it's no longer needed.
    //dconf->fobj2 = (fact_obj_t *)malloc(sizeof(fact_obj_t));
    //init_factobj(dconf->fobj2, sconf->obj->options);
    // also not needed, except for maybe with mingw?
#if defined(__MINGW64__)
    dconf->cosiqs = init_tinyqs();
#endif

    dconf->batch_run_override = 0;

    dconf->do_batch = 0;
    if (sconf->do_batch)
    {       
        uint32_t pmax = (sconf->large_prime_max - sconf->pmax) / sconf->obj->qs_obj.gbl_override_bdiv;
        // sconf->large_prime_max / sconf->obj->qs_obj.gbl_override_bdiv;
        dconf->do_batch = 1;
        relation_batch_init(NULL, &dconf->rb, sconf->pmax, 
            pmax, sconf->large_prime_max, sconf->large_prime_max,
            NULL, 0);

        dconf->batch_run_override = 0;
    }

	// used in SIMD optimized versions of tdiv_med
#ifdef USE_8X_MOD_ASM
	dconf->bl_sizes = (uint16_t *)xmalloc_align(32 * sizeof(uint16_t));
	dconf->bl_locs = (uint16_t *)xmalloc_align(32 * sizeof(uint16_t));
#endif

    dconf->polyscratch = (uint32_t *)xmalloc_align(16 * sizeof(uint32_t));

	//initialize some counters
	dconf->tot_poly = 0;		//track total number of polys
	dconf->num = 0;				//sieve locations subjected to trial division
    dconf->lp_scan_failures = 0;
    dconf->total_reports = 0;
    dconf->total_surviving_reports = 0;
    dconf->total_blocks = 0;

    dconf->num_scatter_opp = 0;
    dconf->num_scatter = 0;

    dconf->lcg_state = hash64(lcg_rand_64(&sconf->obj->lcg_state));

	return 0;
}

int siqs_static_init(static_conf_t* sconf, int is_tiny)
{
    //find the best parameters, multiplier, and factor base
    //for the input.  This is an iterative job because we may find
    //a factor of the input while finding the factor base, in which case
    //we log that fact, and start over with new parameters, multiplier etc.
    //also allocate space for things that don't change during the sieving
    //process.  this scratch space is shared among all threads.

    fact_obj_t* obj = sconf->obj;
    uint32_t i, memsize;
    uint32_t closnuf;
    double sum, avg, sd;
    uint32_t dlp_cutoff;
    uint32_t tlp_cutoff;
    int VFLAG = sconf->obj->VFLAG;
    int THREADS = sconf->obj->THREADS;

    // this pretty much has to stay "8".  the reason is that many of the specialized routines
    // have picky requirements about how large or small the primes can be for them
    // to work correctly.  It seems we have trouble finding indices within the factor
    // base that simultaneously provide the correct size boundaries and are also
    // aligned to a multiple of 16.
    // When using AVX2, we'll maybe have to perform a round or two of SSE41 in a couple
    // places to cross over size boundaries within the factor base, if those boundaries
    // are not divisible by 16.
    int nump = 8;

    if (VFLAG > 2)
    {
        printf("static memory usage:\n");
    }

    // some things work different if the input is tiny
    sconf->is_tiny = is_tiny;

    //default parameters
    sconf->fudge_factor = 1.3;
    sconf->large_mult = 30;
    sconf->num_blocks = 40;
    sconf->num_extra_relations = 64;

    // can raise this if we also raise TFSm and lower TF bounds
    // and get about the same speed... 
    sconf->small_limit = 256;
    sconf->use_dlp = 0;

    // function pointer to the sieve array scanner
    scan_ptr = NULL;


    // relevant:
    // https://software.intel.com/en-us/node/523363
    // https://gcc.gnu.org/onlinedocs/gcc-4.9.2/gcc/X86-Built-in-Functions.html


    // fat binary assignments.  First assign the lowest level core functions.  
    // Sometimes these assume at least SSE2.

    // sieve core functions
    med_sieve_ptr = &med_sieveblock_32k;
    lp_sieveblock_ptr = &lp_sieveblock;

    // poly core function
    firstRoots_ptr = &firstRoots_32k;
    nextRoots_ptr = &nextRoots_32k;
    testRoots_ptr = &testfirstRoots_32k;

    // tdiv core functions
    tdiv_med_ptr = &tdiv_medprimes_32k;
    tdiv_LP_ptr = &tdiv_LP_sse2;
    resieve_med_ptr = &resieve_medprimes_32k;

    // now override the default assignments based on first, what was
    // available during compilation, and second, what is available 
    // at runtime.
#if defined(USE_AVX512F)
    if (obj->HAS_AVX512F)
    {
        tdiv_med_ptr = &tdiv_medprimes_32k_avx2;
        tdiv_LP_ptr = &tdiv_LP_avx512;
        if (sconf->obj->VFLAG > 1)
        {
            printf("assigning tdiv_LP_avx512 ptr\n");
            printf("assigning tdiv_medprimes_32k_avx2 ptr\n");
        }

#if defined(USE_AVX512BW)
        if (obj->HAS_AVX512BW)
        {
            resieve_med_ptr = &resieve_medprimes_32k_avx512bw;
            med_sieve_ptr = &med_sieveblock_32k_avx512bw;
            lp_sieveblock_ptr = &lp_sieveblock_avx512bw;
            if (sconf->obj->VFLAG > 1)
            {
                printf("assigning resieve_medprimes_32k_avx512bw ptr\n");
                printf("assigning lp_sieveblock_avx512bw ptr\n");
                printf("assigning med_sieveblock_32k_avx512bw ptr\n");
            }
        }
        else
        {
            resieve_med_ptr = &resieve_medprimes_32k_avx2;
            lp_sieveblock_ptr = &lp_sieveblock_avx512f;
            med_sieve_ptr = &med_sieveblock_32k_avx2;
            if (sconf->obj->VFLAG > 1)
            {
                printf("assigning resieve_medprimes_32k_avx2 ptr\n");
                printf("assigning lp_sieveblock_avx512f ptr\n");
                printf("assigning med_sieveblock_32k_avx2 ptr\n");
            }
        }

#else
        resieve_med_ptr = &resieve_medprimes_32k_avx2;
        lp_sieveblock_ptr = &lp_sieveblock_avx512f;
        med_sieve_ptr = &med_sieveblock_32k_avx2;
        if (sconf->obj->VFLAG > 1)
        {
            printf("assigning resieve_medprimes_32k_avx2 ptr\n");
            printf("assigning lp_sieveblock_avx512f ptr\n");
            printf("assigning med_sieveblock_32k_avx2 ptr\n");
        }
#endif


    }
    else if (obj->HAS_AVX2)
    {
        if (sconf->obj->VFLAG > 1)
        {
            printf("assigning tdiv_medprimes_32k_avx2 ptr\n");
            printf("assigning tdiv_LP_avx2 ptr\n");
            printf("assigning resieve_medprimes_32k_avx2 ptr\n");
            printf("assigning med_sieveblock_32k_avx2 ptr\n");
        }
        resieve_med_ptr = &resieve_medprimes_32k_avx2;
        tdiv_med_ptr = &tdiv_medprimes_32k_avx2;
        tdiv_LP_ptr = &tdiv_LP_avx2;
        med_sieve_ptr = &med_sieveblock_32k_avx2;

#if defined(_MSC_VER)
        // the avx2 code path for nextroots involves lots of inline 
        // ASM, so visual studio builds can't use it.
        if (obj->HAS_BMI2)
        {
            nextRoots_ptr = &nextRoots_32k_avx2_intrin;
            if (sconf->obj->VFLAG > 1)
            {
                printf("assigning nextRoots_32k_avx2_intrin ptr\n");
            }
        }
        else
        {
            nextRoots_ptr = &nextRoots_32k_sse41;
            if (sconf->obj->VFLAG > 1)
            {
                printf("assigning nextRoots_32k_sse41 ptr\n");
            }
        }
#else

        if (obj->HAS_BMI2)
        {
            if (VFLAG > 1)
            {
                printf("assigning nextRoots_32k_avx2_intrin ptr\n");
            }
            nextRoots_ptr = &nextRoots_32k_avx2_intrin;
        }
        else
        {
            if (VFLAG > 1)
            {
                printf("assigning nextRoots_32k_avx2 ptr\n");
            }
            nextRoots_ptr = &nextRoots_32k_avx2;
        }

#endif
    }
    else if (obj->HAS_SSE41)
    {
        nextRoots_ptr = &nextRoots_32k_sse41;
        med_sieve_ptr = &med_sieveblock_32k_sse41;
        if (sconf->obj->VFLAG > 1)
        {
            printf("assigning nextRoots_32k_sse41 ptr\n");
            printf("assigning med_sieveblock_32k_sse41 ptr\n");
        }

    }

#elif defined(USE_AVX2)
    if (obj->HAS_AVX2)
    {
        if (sconf->obj->VFLAG > 1)
        {
            printf("assigning tdiv_medprimes_32k_avx2 ptr\n");
            printf("assigning tdiv_LP_avx2 ptr\n");
            printf("assigning resieve_medprimes_32k_avx2 ptr\n");
            printf("assigning med_sieveblock_32k_avx2 ptr\n");
        }
        resieve_med_ptr = &resieve_medprimes_32k_avx2;
        tdiv_med_ptr = &tdiv_medprimes_32k_avx2;
        tdiv_LP_ptr = &tdiv_LP_avx2;
        med_sieve_ptr = &med_sieveblock_32k_avx2;

#if defined(_MSC_VER)
        // the avx2 code path for nextroots involves lots of inline 
        // ASM, so visual studio builds can't use it.

        if (obj->HAS_BMI2)
        {
            nextRoots_ptr = &nextRoots_32k_avx2_intrin;
            if (sconf->obj->VFLAG > 1)
            {
                printf("assigning nextRoots_32k_avx2_intrin ptr\n");
            }
        }
        else
        {
            nextRoots_ptr = &nextRoots_32k_sse41;
            if (sconf->obj->VFLAG > 1)
            {
                printf("assigning nextRoots_32k_sse41 ptr\n");
            }
        }
#else
        
        if (obj->HAS_BMI2)
        {
            if (VFLAG > 1)
            {
                printf("assigning nextRoots_32k_avx2_intrin ptr\n");
            }
            nextRoots_ptr = &nextRoots_32k_avx2_intrin;
        }
        else
        {
            if (VFLAG > 1)
            {
                printf("assigning nextRoots_32k_avx2 ptr\n");
            }
            nextRoots_ptr = &nextRoots_32k_avx2;
        }

        //if (VFLAG > 1)
        //{
        //    printf("assigning nextRoots_32k_sse41 ptr\n");
        //}
        //nextRoots_ptr = &nextRoots_32k_sse41;
#endif
    }
    else if (obj->HAS_SSE41)
    {
        nextRoots_ptr = &nextRoots_32k_sse41;
        med_sieve_ptr = &med_sieveblock_32k_sse41;
        if (sconf->obj->VFLAG > 1)
        {
            printf("assigning nextRoots_32k_sse41 ptr\n");
            printf("assigning med_sieveblock_32k_sse41 ptr\n");
        }

    }


#elif defined(USE_SSE41)

    nextRoots_ptr = &nextRoots_32k_sse41;
    med_sieve_ptr = &med_sieveblock_32k_sse41;
    if (sconf->obj->VFLAG > 1)
    {
        printf("assigning nextRoots_32k_sse41 ptr\n");
        printf("assigning med_sieveblock_32k_sse41 ptr\n");
    }

#endif
		
#if defined( __amd64__ ) && defined(USE_AVX512F)
    // amd eypc (zen4) was slower when using the avx512 variants of these
    lp_sieveblock_ptr = &lp_sieveblock;
    tdiv_LP_ptr = &tdiv_LP_avx2;
#endif

	sconf->qs_blocksize = 32768;
	sconf->qs_blockbits = 15;

	//allocate the space for the factor base structure
	sconf->factor_base = (fb_list *)malloc(sizeof(fb_list));

	//allocate space for a copy of input number in the job structure
	mpz_init(sconf->n); //zInit(&sconf->n);
	mpz_set(sconf->n, obj->qs_obj.gmp_n);

	//initialize some constants
	mpz_init(sconf->sqrt_n);
	mpz_init(sconf->target_a);

	//initialize the bookkeeping for tracking partial relations
	sconf->components = 0;
	sconf->vertices = 0;
	sconf->num_cycles = 0;
	sconf->num_relations = 0;
	//force this to happen for now, eventually should implement this flag
	if (1 || !(sconf->obj->flags & MSIEVE_FLAG_SKIP_QS_CYCLES)) {
		sconf->cycle_hashtable = (uint32_t *)xcalloc(
					(size_t)(1 << QS_LOG2_CYCLE_HASH),
					sizeof(uint32_t));
		sconf->cycle_table_size = 1;
		sconf->cycle_table_alloc = 10000;
		sconf->cycle_table = (qs_cycle_t *)xmalloc(
			sconf->cycle_table_alloc * sizeof(qs_cycle_t));
	}

	if (VFLAG > 2)
	{
		memsize = (1 << QS_LOG2_CYCLE_HASH) * sizeof(uint32_t);
		printf("\tinitial cycle hashtable: %d bytes\n",memsize);
		memsize = sconf->cycle_table_alloc * sizeof(qs_cycle_t);
		printf("\tinitial cycle table: %d bytes\n",memsize);
	}

	while (1)
	{
		//look up some parameters - tuned for 64k L1 cache
		get_params(sconf);

		//the number of primes in the factor base is rounded up to the 
		//next multiple of 16, so that we can use aligned moves to speed
		//computations of root updates
		sconf->factor_base->B += (16 - (sconf->factor_base->B % 16));

		//allocate the space for the factor base elements
		sconf->factor_base->list = (fb_element_siqs *)xmalloc_align(
			(size_t)(sizeof(fb_element_siqs)));
		sconf->factor_base->tinylist = (tiny_fb_element_siqs *)xmalloc_align(
			(size_t)(sizeof(tiny_fb_element_siqs)));

		sconf->modsqrt_array = (uint32_t *)xmalloc_align(
			sconf->factor_base->B * sizeof(uint32_t));
		sconf->factor_base->list->prime = (uint32_t *)xmalloc_align(
			(size_t)(sconf->factor_base->B * sizeof(uint32_t)));
        sconf->factor_base->list->binv = (uint64_t*)xmalloc_align(
            (size_t)(sconf->factor_base->B * sizeof(uint64_t)));


		sconf->factor_base->tinylist->prime = (uint32_t *)xmalloc_align(
			(size_t)(512 * sizeof(uint32_t)));
		sconf->factor_base->tinylist->small_inv = (uint32_t *)xmalloc_align(
			(size_t)(512 * sizeof(uint32_t)));
		sconf->factor_base->tinylist->correction = (uint32_t *)xmalloc_align(
			(size_t)(512 * sizeof(uint32_t)));
		sconf->factor_base->tinylist->logprime = (uint32_t *)xmalloc_align(
			(size_t)(512 * sizeof(uint32_t)));

#ifdef USE_8X_MOD_ASM
		sconf->factor_base->list->small_inv = (uint16_t *)xmalloc_align(
			(size_t)(sconf->factor_base->B * sizeof(uint16_t)));
		sconf->factor_base->list->correction = (uint16_t *)xmalloc_align(
			(size_t)(sconf->factor_base->B * sizeof(uint16_t)));
#else
		sconf->factor_base->list->small_inv = (uint32_t *)xmalloc_align(
			(size_t)(sconf->factor_base->B * sizeof(uint32_t)));
		sconf->factor_base->list->correction = (uint32_t *)xmalloc_align(
			(size_t)(sconf->factor_base->B * sizeof(uint32_t)));
#endif
		
		sconf->factor_base->list->logprime = (uint32_t *)xmalloc_align(
			(size_t)(sconf->factor_base->B * sizeof(uint32_t)));


		if (VFLAG > 2)
		{
			memsize = sconf->factor_base->B * sizeof(uint32_t) * 5;
			printf("\tfactor base: %d bytes\n",memsize);
		}

		//find multiplier
		sconf->multiplier = (uint32_t)choose_multiplier_siqs(sconf->factor_base->B, sconf->n);
		mpz_mul_ui(sconf->n, sconf->n, sconf->multiplier);
        sconf->knmod8 = mpz_tdiv_ui(sconf->n, 8);
        //printf("knmod8 = %u\n", sconf->knmod8);

		//sconf holds n*mul, so update its digit count and number of bits
		sconf->digits_n = gmp_base10(sconf->n);
		sconf->bits = mpz_sizeinbase(sconf->n, 2);

		//find sqrt_n
		mpz_sqrt(sconf->sqrt_n, sconf->n);

		//construct the factor base - divide out any primes which 
		//factor n other than mul
		//and go back to get_params if any are found
		if ((i = make_fb_siqs(sconf)) == 0)
			break;
		else
		{
			mpz_t tmpz;
			mpz_init(tmpz);

			//i is already divided out of n.  record the factor we found
			uint64_2gmp(i, tmpz);
			add_to_factor_list(sconf->obj->factors, tmpz, sconf->obj->VFLAG, 
                sconf->obj->NUM_WITNESSES);
			mpz_clear(tmpz);

			//and remove the multiplier we may have added, so that
			//we can try again to build a factor base.
			mpz_tdiv_q_ui(sconf->n, sconf->n, sconf->multiplier);
            align_free(sconf->modsqrt_array);
			align_free(sconf->factor_base->list->prime);
            align_free(sconf->factor_base->list->binv);
			align_free(sconf->factor_base->list->small_inv);
			align_free(sconf->factor_base->list->correction);
			align_free(sconf->factor_base->list->logprime);
			align_free(sconf->factor_base->list);	
		}
	}

#ifdef VDEBUG
	printf("found factor base of %u elements; max = %u\n",
		sconf->factor_base->B,sconf->factor_base->list->prime[sconf->factor_base->B-1]);
#endif

    if (0)
    {
        printf("multiplier: %u\nfb: ", sconf->multiplier);
        sum = 0.0;
        i = 2;
        while (sconf->factor_base->list->prime[i] < 1000)
        {
            printf("%u, ", sconf->factor_base->list->prime[i]);
            sum += log(sconf->factor_base->list->prime[i]);
            i++;
        }
        printf("\nsum of log(primes): %1.3lf\n", sum);
    }

	//the first two fb primes are always the same
	sconf->factor_base->list->prime[0] = 1;		//actually represents -1
	sconf->factor_base->list->prime[1] = 2;

    sconf->sieve_primes = (uint32_t *)xmalloc_align(
        (size_t)(sconf->factor_base->B * sizeof(uint32_t)));

    for (i = 2; i < sconf->factor_base->B; i++)
    {
        sconf->sieve_primes[i] = sconf->factor_base->list->prime[i];
    }

	//adjust for various architectures
	if (sconf->qs_blocksize == 32768)
		sconf->num_blocks *= 1;
	else if (sconf->qs_blocksize == 65536)
		sconf->num_blocks /= 2;
	else
	{
		printf("unknown block size!\n");
		exit(1);
	}

    sconf->in_mem_relations = (siqs_r*)malloc(32768 * sizeof(siqs_r));
    sconf->buffered_rel_alloc = 32768;
    sconf->buffered_rels = 0;

	if (sconf->num_blocks < 1)
		sconf->num_blocks = 1;

	//and adjust for really small jobs
	if (sconf->bits <= 140)
		sconf->num_blocks = 1;

	//set the sieve interval.  this many blocks on each side of 0
	sconf->sieve_interval = sconf->qs_blocksize * sconf->num_blocks;

	//compute sieving limits
	sconf->factor_base->small_B = MIN(
		sconf->factor_base->B,1024); //((INNER_BLOCKSIZE)/(sizeof(sieve_fb))));

	//test the contribution of the small primes to the sieve.  
	for (i = 2; i < sconf->factor_base->B; i++)
	{
		if (sconf->factor_base->list->prime[i] > sconf->small_limit)
			break;
	}
	sconf->sieve_small_fb_start = i;

    // ======================================================================
    // set thresholds at various bit levels.  Note: many places in the code
    // are extremely sensitive to these thresholds; the program can hang or
    // crash unexpectedly if set incorrectly.  Adjust carefully!
    // ======================================================================

	for (; i < sconf->factor_base->B; i++)
	{
		//find the point at which factor base primes exceeds 10 bits.  
		//wait until the index is a multiple of 8 so that we can enter
		//this region of primes aligned on a 16 byte boundary and thus be able to use
		//movdqa
		//don't let med_B grow larger than 1.5 * the blocksize
		if ((sconf->factor_base->list->prime[i] > 1024)  &&
			(i % nump == 0)) break;
	}
	sconf->factor_base->fb_10bit_B = i;

	for (; i < sconf->factor_base->B; i++)
	{
		//find the point at which factor base primes exceeds 11 bits.  
		//wait until the index is a multiple of 8 so that we can enter
		//this region of primes aligned on a 16 byte boundary and thus be able to use
		//movdqa
		//don't let med_B grow larger than 1.5 * the blocksize
		if ((sconf->factor_base->list->prime[i] > 2048)  &&
			(i % nump == 0)) break;
	}
	sconf->factor_base->fb_11bit_B = i;

	for (; i < sconf->factor_base->B; i++)
	{
		//find the point at which factor base primes exceeds 12 bits.  
		//wait until the index is a multiple of 8 so that we can enter
		//this region of primes aligned on a 16 byte boundary and thus be able to use
		//movdqa
        if ((sconf->factor_base->list->prime[i] > 4096) &&
            (i % nump == 0)) break;
	}
	sconf->factor_base->fb_12bit_B = i;

	for (; i < sconf->factor_base->B; i++)
	{
		//find the point at which factor base primes exceeds 13 bits.  
		//wait until the index is a multiple of 8 so that we can enter
		//this region of primes aligned on a 16 byte boundary and thus be able to use
		//movdqa
        if ((sconf->factor_base->list->prime[i] > 8192) &&
            (i % nump == 0)) break;
	}
	sconf->factor_base->fb_13bit_B = i;

	// the prime at which we want to start testing the asm code below is the point
	// at which we have at most 3 increments of the root.  this point starts when
	// prime > 32768 / 3 = 10922.
	for (; i < sconf->factor_base->B; i++)
	{
		//find the point at which factor base primes exceeds 13 bits.  
		//wait until the index is a multiple of 8 so that we can enter
		//this region of primes aligned on a 16 byte boundary and thus be able to use
		//movdqa
		//don't let med_B grow larger than 1.5 * the blocksize
		if ((sconf->factor_base->list->prime[i] > 10922)  &&
			(i % nump == 0)) break;
	}
	// put the upper bound just before prime exceeds blocksize/3, required by SSE2 sieving
	sconf->factor_base->fb_32k_div3 = i;

	for (; i < sconf->factor_base->B; i++)
	{
		//find the point at which factor base primes exceeds 14 bits.  
		//wait until the index is a multiple of 8 so that we can enter
		//this region of primes aligned on a 16 byte boundary and thus be able to use
		//movdqa
		if ((sconf->factor_base->list->prime[i] > 16384)  &&
            (i % nump == 0)) break;
	}	
	// put the upper bound just before prime exceeds 14 bits, required by SSE2 sieving
	sconf->factor_base->fb_14bit_B = i;

	for (; i < sconf->factor_base->B; i++)
	{
		//find the point at which factor base primes exceeds 15 bits.  
		//wait until the index is a multiple of 8 so that we can enter
		//this region of primes aligned on a 16 byte boundary and thus be able to use
		//movdqa
		if ((sconf->factor_base->list->prime[i] > 32768)  &&
            (i % 32 == 0)) break;
	}
	sconf->factor_base->fb_15bit_B = i;

	for (; i < sconf->factor_base->B; i++)
	{
		// find the point at which factor base primes exceed the blocksize.
        // this controls the point at which we switch from medsieve to 
        // lpsieve (and bucket sorting).
		// wait until the index is a multiple of 16 so that we can enter
		// this region of primes aligned on a 16 byte boundary and thus be able to use
		// movdqa
		// don't let med_B grow larger than ~1.5 * the blocksize
#if defined(USE_AVX512F) && !defined(_MSC_VER)
        if ((sconf->factor_base->list->prime[i] > (uint32_t)(1.5 * (double)sconf->qs_blocksize)) &&
            ((i % 16) == 0))
            break;
#else
        // bucket sieve not as efficient without compressstoreu, so 
        // do a little more work in medsieve before switching.
        if ((sconf->factor_base->list->prime[i] > (uint32_t)(1.75 * (double)sconf->qs_blocksize)) &&
            ((i % 16) == 0))
            break;
#endif
		

		//or 2^16, whichever is smaller
		if (sconf->factor_base->list->prime[i] > 65536)
		{
            printf("probably never get here\n");
			i -= i%16;
			break;
		}

		//or, of course, the entire factor base (loop condition)
	}
	sconf->factor_base->med_B = i;

	for (; i < sconf->factor_base->B; i++)
	{
		//find the point at which factor base primes exceed the size of the sieve 
		//interval.  wait until the index is a multiple of 16 so that we can enter
		//this region of primes aligned on a 16 byte boundary and thus be able to use
		//movdqa
		if ((sconf->factor_base->list->prime[i] > sconf->sieve_interval) &&
			((i % 16) == 0))
		{
			i -= 16;
			break;
		}
	}
	sconf->factor_base->large_B = i;

    // when p > sieve_interval we use a special bucket sieve tailored for primes that
    // will at most hit the interval once.  The AVX512 version of that method uses 
    // compressstoreu for each root on each block for a total of 2*numblocks calls per side.
    // (along with 2 popcnt_u32 calls as well for each block on each side.)
    // As primes continue to get larger and we have more blocks, at some point it
    // is better to switch back to the scatter loop.  The scatter loop uses 4 compressstoreu
    // calls followed by N writes, where N is the number of interval hits for both roots.
    // if we assume a compressstoreu is about the same as a normal write, then it makes
    // sense to switch when (N + 4) < (2 * numblocks).  To estimate N we find the area
    // under the triangle formed by the large sieve primes and their probability of
    // hitting the interval.
    for (; i < sconf->factor_base->B; i++)
    {
        if ((sconf->factor_base->list->prime[i] > 6 * sconf->sieve_interval) &&
            ((i % 16) == 0) && (((sconf->factor_base->B - i) % 16) == 0))
            break;
    }
    sconf->factor_base->x2_large_B = i; // sconf->factor_base->B;

#ifdef USE_SS_SEARCH
    if ((sconf->obj->qs_obj.gbl_override_ssidx_flag))
    {
        for (; i < sconf->factor_base->B; i++)
        {
            // when using subset-sum searching, use this as the switchover
            // point.  It occurs much higher up in the fb compared to 
            // standard poly enumeration.
            if ((sconf->obj->qs_obj.gbl_override_ssidx_flag) &&
                (sconf->obj->qs_obj.gbl_override_ssidx > 0))
            {
                if ((i > sconf->obj->qs_obj.gbl_override_ssidx) &&
                    ((i % 16) == 0) && (((sconf->factor_base->B - i) % 16) == 0))
                    break;
            }
            else if ((sconf->obj->qs_obj.gbl_override_ssidx_flag) &&
                (sconf->obj->qs_obj.gbl_override_ssidx < 0))
            {
                // -1 indicates default, current set to to top 1/3 of the factor base
                if ((i > sconf->factor_base->B * 2 / 3) &&
                    ((i % 16) == 0) && (((sconf->factor_base->B - i) % 16) == 0))
                    break;
            }
        }
        sconf->factor_base->ss_start_B = i;

#if defined(USE_DIRECT_SIEVE_SS)
        if ((sconf->obj->qs_obj.gbl_override_ssidx_flag) &&
            (sconf->obj->qs_obj.gbl_override_ssidx > 0))
        {
#ifdef USE_AVX512F
            sconf->factor_base->ss_start_B = sconf->factor_base->x2_large_B;
#else
            sconf->factor_base->ss_start_B = sconf->factor_base->med_B;
#endif
        }
        else
        {
            sconf->factor_base->ss_start_B = sconf->factor_base->B;
        }
#endif
    }
    else
    {
        sconf->factor_base->ss_start_B = sconf->factor_base->B;
    }

    // larger slices means more primes per slice which means the buckets
    // need to be bigger per slice.  also of course, larger slices
    // means fewer maximum blocks allowed, because we need more bits to
    // store prime id's and so we have fewer bits for storing root locations.
    sconf->factor_base->slice_size = SS_SLICE_SIZE;
    int numslices = (int)((sconf->factor_base->B - sconf->factor_base->ss_start_B) /
        sconf->factor_base->slice_size) + 1;
    sconf->factor_base->num_ss_slices = numslices;
#else
    sconf->factor_base->ss_start_B = sconf->factor_base->B;
#endif
    
    sconf->pmax = sconf->factor_base->list->prime[sconf->factor_base->B - 1];

    if (sconf->factor_base->large_B < sconf->factor_base->med_B)
        sconf->factor_base->large_B = sconf->factor_base->med_B;

    
	if (VFLAG > 1)
	{
		printf("fb bounds\n\tsmall: %u\n\tSPV: %u\n\t10bit: %u\n\t11bit: %u\n\t12bit: %u\n\t"
			"13bit: %u\n\t32k div 3: %u\n\t14bit: %u\n\t15bit: %u\n\tmed: %u\n\tlarge: %u\n"
            "\tlarge_x2: %u\n\tall: %u\n",
			sconf->factor_base->small_B,
			sconf->sieve_small_fb_start,
			sconf->factor_base->fb_10bit_B,
			sconf->factor_base->fb_11bit_B,
			sconf->factor_base->fb_12bit_B,
			sconf->factor_base->fb_13bit_B,
			sconf->factor_base->fb_32k_div3,
			sconf->factor_base->fb_14bit_B,
			sconf->factor_base->fb_15bit_B,
			sconf->factor_base->med_B,
			sconf->factor_base->large_B,
            sconf->factor_base->x2_large_B,
			sconf->factor_base->B);

		printf("start primes\n\tSPV: %u\n\t10bit: %u\n\t11bit: %u\n\t12bit: %u\n\t"
			"13bit: %u\n\t32k div 3: %u\n\t14bit: %u\n\t15bit: %u\n\tmed: %u\n\tlarge: %u\n"
            "\tlarge_x2: %u\n\tall: %u\n",
			sconf->factor_base->list->prime[sconf->sieve_small_fb_start - 1],
			sconf->factor_base->list->prime[sconf->factor_base->fb_10bit_B-1],
			sconf->factor_base->list->prime[sconf->factor_base->fb_11bit_B-1],
			sconf->factor_base->list->prime[sconf->factor_base->fb_12bit_B-1],
			sconf->factor_base->list->prime[sconf->factor_base->fb_13bit_B-1],
			sconf->factor_base->list->prime[sconf->factor_base->fb_32k_div3-1],
			sconf->factor_base->list->prime[sconf->factor_base->fb_14bit_B-1],
			sconf->factor_base->list->prime[sconf->factor_base->fb_15bit_B-1],
			sconf->factor_base->list->prime[sconf->factor_base->med_B-1],
			sconf->factor_base->list->prime[sconf->factor_base->large_B-1],
            sconf->factor_base->list->prime[sconf->factor_base->x2_large_B - 1],
            sconf->pmax);
	}

    //printf("large mult is %u, pmax is %u\n", sconf->large_mult, sconf->pmax);
    //printf("large prime max is %u, gbl_override_lpb is %u\n", 
    //    sconf->large_prime_max, sconf->obj->qs_obj.gbl_override_lpb);

	// lpmax from parameters and/or user specifications
    if ((sconf->large_mult > 1000) || (sconf->large_mult == 0))
    {
        // this is an actual LPB for a TLP job
        sconf->large_prime_max = sconf->large_mult;
        if (sconf->large_prime_max == 0)
            sconf->large_prime_max = 4294967295ULL;
    }
	else if ((4294967295ULL / sconf->large_mult) < sconf->pmax)
	{
		// job is so big that pmax * default large_mult won't fit in 32 bits
		// reduce large_mult accordingly
		sconf->large_mult = 4294967295ULL / sconf->pmax;
		sconf->large_prime_max = sconf->pmax * sconf->large_mult;
	}
	else
	{
		sconf->large_prime_max = sconf->pmax * sconf->large_mult;
	}

	if (sconf->obj->qs_obj.gbl_override_lpb > 0)
	{
        if (sconf->obj->qs_obj.gbl_override_lpb >= 32)
        {
            sconf->large_prime_max = 4294967295;
        }
        else
        {
            sconf->large_prime_max = (1 << sconf->obj->qs_obj.gbl_override_lpb);
        }
	}

	// based on the size of the input, determine how to proceed.
    // this should maybe be tuned based on machine type and/or other
    // factors as well, not just instruction set used during compile.
#if (defined(USE_AVX2) || defined(USE_SSE41))
    if (sconf->obj->HAS_SSE41 || sconf->obj->HAS_AVX2 || sconf->obj->HAS_AVX512F)
        dlp_cutoff = 70;
    else
        dlp_cutoff = 77;
#else
    dlp_cutoff = 77;
#endif

    // don't use tlp unless forced.  too much can go wrong still.
    tlp_cutoff = 221;

    // could maybe someday change this w.r.t input size... for now
    // just test it out...
    //sconf->poly_batch_size = 1;
    //printf("choosing scan algorithm based on HAS_AVX512F=%d, HAS_AVX2=%d, dlp_cutoff=%d\n",
    //    sconf->obj->HAS_AVX512F, sconf->obj->HAS_AVX2, dlp_cutoff);

	if (sconf->digits_n >= tlp_cutoff)
	{
		sconf->use_dlp = 2;
        sconf->num_lp = 3;
        // NOTE: for TF cutoffs above 127 bits, only the avx512 version
        // of scan_relations will work.
#ifdef USE_AVX512F
        if (sconf->obj->HAS_AVX512F) scan_ptr = &check_relations_siqs_16_avx512;
        else if (sconf->obj->HAS_AVX2) scan_ptr = &check_relations_siqs_16_avx2;
        else scan_ptr = &check_relations_siqs_16_sse2;
#elif defined (USE_AVX2)
        if (sconf->obj->HAS_AVX2) scan_ptr = &check_relations_siqs_16_avx2;
        else scan_ptr = &check_relations_siqs_16_sse2;
#else
        scan_ptr = &check_relations_siqs_16_sse2;
#endif
		sconf->scan_unrolling = 128;
	}
    else if (sconf->digits_n >= dlp_cutoff)
	{
		sconf->use_dlp = 1;
        sconf->num_lp = 2;
#ifdef USE_AVX512F
        if (sconf->obj->HAS_AVX512F) scan_ptr = &check_relations_siqs_16_avx512;
        else if (sconf->obj->HAS_AVX2) scan_ptr = &check_relations_siqs_16_avx2;
        else scan_ptr = &check_relations_siqs_16_sse2;
#elif defined (USE_AVX2)
        if (sconf->obj->HAS_AVX2) scan_ptr = &check_relations_siqs_16_avx2;
        else scan_ptr = &check_relations_siqs_16_sse2;
#else
        scan_ptr = &check_relations_siqs_16_sse2;
#endif
		sconf->scan_unrolling = 128;
	}
	else
	{
        sconf->num_lp = 1;
		if (sconf->digits_n < 30)
		{
#ifdef USE_AVX2
            if (sconf->obj->HAS_AVX2) scan_ptr = &check_relations_siqs_4_avx2;
            else scan_ptr = &check_relations_siqs_4_sse2;
#else
            scan_ptr = &check_relations_siqs_4_sse2;
#endif
			sconf->scan_unrolling = 32;
		}
		else if (sconf->digits_n < 60)
		{
#ifdef USE_AVX2
            if (sconf->obj->HAS_AVX2) scan_ptr = &check_relations_siqs_4_avx2;
            else scan_ptr = &check_relations_siqs_4_sse2;
#else
            scan_ptr = &check_relations_siqs_4_sse2;
#endif
			sconf->scan_unrolling = 32;
		}
		else
		{
#ifdef USE_AVX2
            if (sconf->obj->HAS_AVX2) scan_ptr = &check_relations_siqs_8_avx2;
            else scan_ptr = &check_relations_siqs_8_sse2;
#else
            scan_ptr = &check_relations_siqs_8_sse2;
#endif
			sconf->scan_unrolling = 64;
		}
		sconf->use_dlp = 0;
	}

    if (sconf->obj->qs_obj.gbl_force_QLP)
    {
        sconf->use_dlp = 3;
        sconf->num_lp = 4;
        // NOTE: for TF cutoffs above 127 bits, only the avx512 version
        // of scan_relations will work.
#ifdef USE_AVX512F
        if (sconf->obj->HAS_AVX512F) scan_ptr = &check_relations_siqs_16_avx512;
        else if (sconf->obj->HAS_AVX2) scan_ptr = &check_relations_siqs_16_avx2;
        else scan_ptr = &check_relations_siqs_16_sse2;
#elif defined(USE_AVX2)
        if (sconf->obj->HAS_AVX2) scan_ptr = &check_relations_siqs_16_avx2;
        else scan_ptr = &check_relations_siqs_16_sse2;
#else
        scan_ptr = &check_relations_siqs_16_sse2;
#endif
        sconf->scan_unrolling = 128;
    }

	if (sconf->obj->qs_obj.gbl_force_TLP)
	{
		sconf->use_dlp = 2;
        sconf->num_lp = 3;
        // NOTE: for TF cutoffs above 127 bits, only the avx512 version
        // of scan_relations will work.
#ifdef USE_AVX512F
        if (sconf->obj->HAS_AVX512F) scan_ptr = &check_relations_siqs_16_avx512;
        else if (sconf->obj->HAS_AVX2) scan_ptr = &check_relations_siqs_16_avx2;
        else scan_ptr = &check_relations_siqs_16_sse2;
#elif defined(USE_AVX2)
        if (sconf->obj->HAS_AVX2) scan_ptr = &check_relations_siqs_16_avx2;
        else scan_ptr = &check_relations_siqs_16_sse2;
#else
        scan_ptr = &check_relations_siqs_16_sse2;
#endif
		sconf->scan_unrolling = 128;
	}
	
    if (sconf->obj->qs_obj.gbl_force_DLP)
	{
		sconf->use_dlp = 1;
        sconf->num_lp = 2;
#ifdef USE_AVX512F
        if (sconf->obj->HAS_AVX512F) scan_ptr = &check_relations_siqs_16_avx512;
        else if (sconf->obj->HAS_AVX2) scan_ptr = &check_relations_siqs_16_avx2;
        else scan_ptr = &check_relations_siqs_16_sse2;
#elif defined(USE_AVX2)
        if (sconf->obj->HAS_AVX2) scan_ptr = &check_relations_siqs_16_avx2;
        else scan_ptr = &check_relations_siqs_16_sse2;
#else
        scan_ptr = &check_relations_siqs_16_sse2;
#endif
		sconf->scan_unrolling = 128;
	}

    sconf->do_batch = 1;
    if (sconf->obj->qs_obj.gbl_override_3lp_bat)
    {
        sconf->do_batch = 0;
    }
    sconf->batch_run_override = 0;

#if !defined( __MINGW64__)
    // not sure why, but batch factoring completely fails when using mingw64.
    if ((sconf->use_dlp >= 2) && (sconf->do_batch == 1))
    {
        uint32_t memuse = 0;
        struct timeval locstart, locstop;

        sconf->rb = (relation_batch_t *)xmalloc(obj->THREADS * sizeof(relation_batch_t));
        sconf->num_alloc_rb = obj->THREADS;
        sconf->num_active_rb = 0;
        memuse += obj->THREADS * sizeof(relation_batch_t);

        //relation_batch_init(stdout, &sconf->rb, sconf->pmax, sconf->large_prime_max,
        //    sconf->large_prime_max, sconf->large_prime_max, NULL, (print_relation_t)NULL, 1);

        // compute the product of primes only with one thread.
        gettimeofday(&locstart, NULL);

        uint32_t pmax = (sconf->large_prime_max - sconf->pmax) / sconf->obj->qs_obj.gbl_override_bdiv; 
        sconf->rb[0].num_uecm[0] = 0;
        sconf->rb[0].num_uecm[1] = 0;
        sconf->rb[0].num_uecm[2] = 0;
        sconf->rb[0].num_uecm[3] = 0;
        sconf->rb[0].num_tecm = 0;
        sconf->rb[0].num_tecm2 = 0;
        sconf->rb[0].num_qs = 0;
        sconf->rb[0].num_attempt = 0;
        sconf->rb[0].num_success = 0;
        for (i = 0; i < 8; i++)
        {
            sconf->rb[0].num_abort[i] = 0;
        }

        relation_batch_init(stdout, &sconf->rb[0], sconf->pmax, pmax,
            sconf->large_prime_max, sconf->large_prime_max, NULL, 1);
        memuse += sconf->rb[0].num_relations_alloc * sizeof(cofactor_t);
        memuse += sconf->rb[0].num_factors_alloc * sizeof(uint32_t);
        memuse += mpz_sizeinbase(sconf->rb[0].prime_product, 2) / 8;

        gettimeofday(&locstop, NULL);
        if (obj->VFLAG > 0)
        {
            printf("batch init took %1.4f sec\n", ytools_difftime(&locstart, &locstop));
        }

        sconf->rb[0].target_relations = sconf->obj->qs_obj.gbl_btarget;

        for (i = 1; i < sconf->num_alloc_rb; i++)
        {
            // allocate space for each thread to buffer relations.
            relation_batch_init(stdout, &sconf->rb[i], sconf->pmax, pmax,
                sconf->large_prime_max, sconf->large_prime_max, NULL, 0);
            memuse += sconf->rb[i].num_relations_alloc * sizeof(cofactor_t);
            memuse += sconf->rb[i].num_factors_alloc * sizeof(uint32_t);

            // copy the prime product from the first thread.
            // to save memory, we could maybe defer this to when a relation buffer is
            // first used?  We should only ever need a fraction of the threads
            // to act a merged relation buffers.
            //mpz_set(sconf->rb[i].prime_product, sconf->rb[0].prime_product);
            //memuse += mpz_sizeinbase(sconf->rb[0].prime_product, 2) / 8;
            mpz_set_ui(sconf->rb[i].prime_product, 0);

            sconf->rb[i].target_relations = sconf->obj->qs_obj.gbl_btarget;
        }

        sconf->num_active_rb = 0;
        sconf->max_active_rb = 0;
        sconf->batch_buffer_id = 1;
        sconf->batch_run_override = 0;
        sconf->do_batch = 1;

        sconf->do_periodic_tlp_filter = 1;
        
        if (obj->VFLAG > 0)
        {
            printf("memory use for relation batch structures: %u bytes\n", memuse);
        }
    }
#endif

	savefile_init(&obj->qs_obj.savefile, sconf->obj->qs_obj.siqs_savefile);

	// default values
    if (sconf->bits < 270)
    {
        sconf->dlp_exp = 1.75;
    }
    else if (sconf->bits < 290)
    {
        sconf->dlp_exp = 1.8;
    }
    else if (sconf->bits < 330)
    {
        sconf->dlp_exp = 1.9;
    }
    else
    {
        sconf->dlp_exp = 1.95;
    }
	sconf->tlp_exp = 2.7;
    sconf->qlp_exp = 3.4;

	// check for user overrides
	if (fabs(sconf->obj->qs_obj.gbl_override_mfbt) > 1e-9)     // - 2.8 
	{
		sconf->tlp_exp = sconf->obj->qs_obj.gbl_override_mfbt;
	}

	if (fabs(sconf->obj->qs_obj.gbl_override_mfbd) > 1e-9)      //  - 1.85
	{
		sconf->dlp_exp = sconf->obj->qs_obj.gbl_override_mfbd;
	}

    if (fabs(sconf->obj->qs_obj.gbl_override_mfbq) > 1e-9)      //  - 1.85
    {
        sconf->qlp_exp = sconf->obj->qs_obj.gbl_override_mfbq;
    }

	// if we're using dlp, compute the range of residues which will
	// undergo cofactorization (batch GCD, ECM)
	if (sconf->use_dlp >= 1)
	{
		sconf->max_fb2 = (uint64_t)sconf->pmax * (uint64_t)sconf->pmax;
		sconf->large_prime_max2 = (uint64_t)pow((double)sconf->large_prime_max, sconf->dlp_exp);
	}
	
	if (sconf->use_dlp >= 2)
	{
		sconf->max_fb3 = pow((double)sconf->pmax, 3.0);
		sconf->large_prime_max3 = pow((double)sconf->large_prime_max, sconf->tlp_exp);
	}

    if (sconf->use_dlp >= 3)
    {
        sconf->max_fb4 = pow((double)sconf->pmax, 4.0);
        sconf->large_prime_max4 = pow((double)sconf->large_prime_max, sconf->qlp_exp);
    }

	// 'a' values should be as close as possible to sqrt(2n)/M in order to make
	// values of g_{a,b}(x) as uniform as possible
	mpz_mul_2exp(sconf->target_a, sconf->n, 1);
	mpz_sqrt(sconf->target_a, sconf->target_a);
#if defined (USE_SS_SEARCH) && defined( USE_DIRECT_SIEVE_SS)
    mpz_tdiv_q_ui(sconf->target_a, sconf->target_a, 1024);
#else
    mpz_tdiv_q_ui(sconf->target_a, sconf->target_a, sconf->sieve_interval);
#endif
	
    if (sconf->knmod8 == 1)
    {
        mpz_tdiv_q_2exp(sconf->target_a, sconf->target_a, 1);
    }


	// compute the number of bits in M/2*sqrt(N/2), the approximate value
	// of residues in the sieve interval.  Then subtract some slack.
	// sieve locations greater than this are worthy of trial dividing
	closnuf = (uint8_t)(double)((sconf->bits - 1)/2);
	closnuf += (uint8_t)(log((double)sconf->sieve_interval/2)/log(2.0));
	closnuf -= (uint8_t)(sconf->fudge_factor * log(sconf->large_prime_max) / log(2.0));

	// contribution of all small primes we're skipping to a block's
	// worth of sieving... compute the average per sieve location
	sum = 0;
	for (i = 2; i < sconf->sieve_small_fb_start; i++)
	{
		uint32_t prime = sconf->factor_base->list->prime[i];

		sum += (double)sconf->factor_base->list->logprime[i] * 
			(double)sconf->qs_blocksize / (double)prime;
	}
	avg = 2*sum/sconf->qs_blocksize;

	// this was observed to be the typical standard dev. for one block of 
	// test sieving... wrap this magic number along with several others into
	// one empirically determined fudge factor...
	sd = sqrt(28);

	// this appears to work fairly well... paper mentioned doing it this
	// way... find out and reference here.
	sconf->tf_small_cutoff = (uint8_t)(avg + 2.5*sd);

	// it pays to trial divide a little more for big jobs...
	// "optimized" is used rather loosely here, but these corrections
	// were observed to make things faster.  

    if ((sconf->use_dlp == 1) || sconf->obj->qs_obj.gbl_force_DLP)
    {
        // empirically, these were observed to work fairly well.
        if(sconf->digits_n < 82)
            closnuf = sconf->digits_n + 7;
        else if (sconf->digits_n < 90)
            closnuf = sconf->digits_n + 5;
        else if (sconf->digits_n < 95)
            closnuf = sconf->digits_n + 3;
        else if (sconf->digits_n < 100)
            closnuf = sconf->digits_n + 1;
        else
            closnuf = sconf->digits_n;

#ifdef USE_AVX512F
        closnuf -= 2;
#endif
    }
    else if ((sconf->use_dlp >= 2) || 
        sconf->obj->qs_obj.gbl_force_TLP || 
        sconf->obj->qs_obj.gbl_force_QLP)
    {
        closnuf = sconf->digits_n - 8;
    }
    else
    {
        closnuf -= sconf->tf_small_cutoff;	//correction to the previous estimate
    }

    if (sconf->obj->qs_obj.gbl_override_tf_flag)
    {
        closnuf = sconf->obj->qs_obj.gbl_override_tf;
    }

    // need the highest bit to be clear in order to scan the sieve array efficiently
    // this means we may check more sieve reports than we otherwise would have, but
    // only for the very largest of jobs which might actually benefit from doing so anyway.
    if ((closnuf >= 128) && (sconf->use_dlp < 2))
        closnuf = 127;

#ifndef USE_AVX512BW
    // so far, detecting sieve hits when tf bound is > 127 only
    // works for AVX512BW enabled cpus.
    if (closnuf >= 128)
        closnuf = 127;
#endif

    sconf->blockinit = closnuf;
    sconf->tf_closnuf = closnuf;

    if (sconf->digits_n < sconf->obj->qs_obj.inmem_cutoff)
    {
        sconf->in_mem = 1;
    }
    else
    {
        sconf->in_mem = 0;
    }

	// needed during filtering
	mpz_init(sconf->curr_a);
	sconf->curr_poly = (siqs_poly *)malloc(sizeof(siqs_poly));
	mpz_init(sconf->curr_poly->mpz_poly_a);
	mpz_init(sconf->curr_poly->mpz_poly_b);
	mpz_init(sconf->curr_poly->mpz_poly_c);
	sconf->curr_poly->qlisort = (int *)malloc(MAX_A_FACTORS*sizeof(int));
	sconf->curr_poly->gray = (char *) malloc( (1 << (MAX_A_FACTORS-1) ) * sizeof(char));
	sconf->curr_poly->nu = (char *) malloc((1 << (MAX_A_FACTORS-1)) * sizeof(char));

	// compute how often to check our list of partial relations and update the gui.
    // for the TLP-variation, this determines when we run filtering... so it is
    // quite a bit larger in that case.  It is reduced as the factorization 
    // progresses (see update_check()).
    if (sconf->use_dlp >= 2)
    {
        // check frequently.  Sieving will take
        // so long that filtering costs are negligible, and good
        // data will emerge from the filtering runs that might hone in
        // cycle formation predictions.
        sconf->check_inc = 1000;        // run filtering after finding this many fulls
    }
    else
    {
        sconf->check_inc = sconf->factor_base->B / 10;
    }

	sconf->check_total = sconf->check_inc;
	sconf->update_time = 5;
    sconf->flag = 0;

	// get ready for the factorization and screen updating
	sconf->t_update=0;
	sconf->last_numfull = 0;
	sconf->last_numpartial = 0;
	sconf->last_numcycles = 0;
	sconf->last_cycrate = 0.0;
	sconf->last_fullrate = 0.0;
	sconf->num_cycles = 0;
	gettimeofday(&sconf->update_start, NULL);

	sconf->failed_squfof = 0;
	sconf->attempted_squfof = 0;
	sconf->failed_cosiqs = 0;
	sconf->attempted_cosiqs = 0;
	sconf->dlp_outside_range = 0;
	sconf->num_full = 0;
	sconf->num_slp = 0;
	sconf->dlp_prp = 0;
	sconf->dlp_useful = 0;
	sconf->tlp_outside_range = 0;
	sconf->tlp_prp = 0;
	sconf->tlp_useful = 0;
    sconf->qlp_outside_range = 0;
    sconf->qlp_prp = 0;
    sconf->qlp_useful = 0;
    sconf->num_found = 0;
	sconf->total_poly_a = 0;	//track number of A polys used
	sconf->num_r = 0;			//total relations found
	sconf->charcount = 0;		//characters on the screen

	sconf->t_time1 = 0;			//sieve time
	sconf->t_time2 = 0;			//relation scanning and trial division
	sconf->t_time3 = 0;			//polynomial calculations and large prime sieving
	sconf->t_time4 = 0;			//extra?
	
	sconf->tot_poly = 0;		//track total number of polys
	sconf->num = 0;				//sieve locations subjected to trial division
    sconf->total_reports = 0;
    sconf->total_surviving_reports = 0;
    sconf->total_blocks = 0;
    sconf->lp_scan_failures = 0;

	// no factors so far...
	sconf->factor_list.num_factors = 0;

    // diagnostics for how often AVX512 gather/scatter operations occur.
    sconf->num_scatter_opp = 0;
    sconf->num_scatter = 0;

    return 0;
}

int update_check(static_conf_t *sconf)
{
	// check to see if we should combine partial relations
	// found and check to see if we're done.  also update the screen.
	// this happens one of two ways, either we have found more than a 
	// certain amount of relations since the last time we checked, or if
	// a certain amount of time has elapsed.
	mpz_t tmp1;
	struct timeval update_stop;
	uint32_t num_full = sconf->num_relations;
	uint32_t check_total = sconf->check_total;
	uint32_t check_inc = sconf->check_inc;
    double update_time = sconf->update_time;
	double t_update, t_time;	
	//int i;
	fb_list *fb = sconf->factor_base;
	int retcode = 0;

    // if we are working on a large problem or have a lot of
    // threads, update as fast as we can.
    if ((sconf->digits_n >= 90) || (sconf->obj->THREADS >= 32))
        update_time = 0;

	mpz_init(tmp1);

	gettimeofday(&update_stop, NULL);
    t_update = ytools_difftime(&sconf->update_start, &update_stop);

	if ((num_full >= check_total) || (t_update > update_time))
	{
		// watch for an abort
		if (SIQS_ABORT)
		{
            // let the threads stop gracefully (merge rels).  
            return 2;
		}

        // if we are only collecting a specified number of relations
        if (sconf->obj->qs_obj.gbl_override_rel_flag)
        {
            if (((sconf->use_dlp < 2) &&
                ((num_full + sconf->num_cycles) > sconf->obj->qs_obj.gbl_override_rel)) ||
                ((sconf->use_dlp < 2) && ((sconf->num_full + sconf->num_slp +
                    sconf->dlp_useful + sconf->tlp_useful) > sconf->obj->qs_obj.gbl_override_rel)))
            {
                // this doesn't work when restarting... the relations loaded
                // count toward the total found so far.  Would need to separately count
                // relations loaded and relations found this run.
                printf("\nMax specified relations found\n");
                mpz_set_ui(tmp1, sconf->tot_poly);					//total number of polys
                mpz_mul_ui(tmp1, tmp1, sconf->num_blocks);	//number of blocks
                mpz_mul_2exp(tmp1, tmp1, 1);			//pos and neg sides
                mpz_mul_2exp(tmp1, tmp1, sconf->qs_blockbits);	//sieve locations per block	

                printf("\nsieve time = %6.4f, relation time = %6.4f, poly_time = %6.4f\n",
                    sconf->t_time1, sconf->t_time2, sconf->t_time3);
                gmp_printf("trial division touched %d sieve locations out of %Zd\n",
                    sconf->num, tmp1);
                fflush(stdout);
                fflush(stderr);

                sconf->obj->qs_obj.gbl_override_rel = num_full + sconf->num_cycles;

                return 2;
            }
        }

        // if we are only running a specified amount of time
        t_time = ytools_difftime(&sconf->totaltime_start, &update_stop);
		if (sconf->obj->qs_obj.gbl_override_time_flag &&
			(t_time > sconf->obj->qs_obj.gbl_override_time))
		{
			printf("\nMax specified time limit reached\n");

			mpz_set_ui(tmp1, sconf->tot_poly);					//total number of polys
			mpz_mul_ui(tmp1, tmp1, sconf->num_blocks);	//number of blocks
			mpz_mul_2exp(tmp1, tmp1, 1);			//pos and neg sides
			mpz_mul_2exp(tmp1, tmp1, sconf->qs_blockbits);	//sieve locations per block				

			printf("\nsieve time = %6.4f, relation time = %6.4f, poly_time = %6.4f\n",
				sconf->t_time1,sconf->t_time2,sconf->t_time3);
			gmp_printf("trial division touched %d sieve locations out of %Zd\n",
				sconf->num, tmp1);
			fflush(stdout);
			fflush(stderr);
			
			return 2;
		}


		if (sconf->use_dlp < 2)
		{
			sconf->num_r = sconf->num_relations +
				sconf->num_cycles +
				sconf->components - sconf->vertices;
		}
		else
		{
			// total number of relations found is only updated
			// when we do a filtering step (below)
			sconf->num_r = sconf->last_numfull + sconf->last_numpartial;
		}

		// update screen
		if (sconf->obj->VFLAG >= 0)
		{
            uint32_t total_rels = sconf->num_full + sconf->num_slp +
                sconf->dlp_useful + sconf->tlp_useful;
            uint32_t rels_left;

            if (check_total > sconf->num_full)
            {
                rels_left = check_total - sconf->num_full;
            }
            else
            {
                rels_left = 0;
            }

			if (sconf->use_dlp >= 2)
			{
                if (sconf->do_batch)
                {
                    char suffix;
                    float batched = sconf->attempted_cosiqs;

                    if (batched > 10000000)
                    {
                        batched /= 1000000.;
                        suffix = 'M';
                    }
                    else
                    {
                        batched /= 1000.;
                        suffix = 'k';
                    }

                    if (sconf->use_dlp == 3)
                    {
                        printf("%u full + %u partial from (%u,%u,%u,%u) lp, "
                            "(%1.1f%c raw-GCD), filt in %u (%1.2f r/sec)\r",
                            sconf->num_full, sconf->last_numpartial, sconf->num_slp,
                            sconf->dlp_useful, sconf->tlp_useful, sconf->qlp_useful, 
                            (float)batched, suffix, rels_left,
                            (double)(sconf->num_relations + sconf->dlp_useful +
                                sconf->tlp_useful) / t_time);
                    }
                    else
                    {
                        printf("%u full + %u partial from (%u,%u,%u) lp, "
                            "(%1.1f%c raw-GCD), filt in %u (%1.2f r/sec)\r",
                            sconf->num_full, sconf->last_numpartial, sconf->num_slp,
                            sconf->dlp_useful, sconf->tlp_useful, (float)batched, suffix,
                            rels_left,
                            (double)(sconf->num_relations + sconf->dlp_useful +
                                sconf->tlp_useful) / t_time);
                    }
                    
                }
                else
                {
                    printf("last: %u, now: %u full, %u slp, "
                        "%u dlp (%ukatt), "
                        "%u tlp (%ukatt), filt in %u (%1.2f r/sec)\r",
                        sconf->last_numfull + sconf->last_numpartial, sconf->num_full, sconf->num_slp,
                        sconf->dlp_useful, sconf->attempted_squfof / 1000,
                        sconf->tlp_useful, sconf->attempted_cosiqs / 1000,
                        rels_left,
                        (double)(sconf->num_relations + sconf->dlp_useful +
                            sconf->tlp_useful) / t_time);
                }
                fflush(stdout);
			}
			else
			{
                // https://www.mersenneforum.org/node/11822?p=1067569#post1067569
                // James Heinrich has determined that the following formula 
                // pretty closely predicts the number of full relations
                // needed to complete a factorization, as a function of input size,
                // as a percentage of the total number of relations needed.
                // from the info we are printing here and the total elasped time
                // so far, we can compute an ETA.
                float fractFullNeeded = 165.0 * pow(sconf->digits_n, -1.4);
                float fractFullHave = (float)sconf->num_relations / (float)sconf->factor_base->B;
                float percentDone = fractFullHave / fractFullNeeded;
                int eta = MAX(0, (int)(t_time / percentDone - t_time));
				printf("%d rels found: %d full + "
					"%d from %d partial, (%6.2f rels/sec, ETA %d sec)\r",
					sconf->num_r, sconf->num_relations,
					sconf->num_cycles +
					sconf->components - sconf->vertices,
					sconf->num_cycles,
					(double)(sconf->num_relations + sconf->num_cycles) / t_time, eta);
			}
            
			fflush(stdout);
		}

        sconf->update_start.tv_sec = update_stop.tv_sec;
        sconf->update_start.tv_usec = update_stop.tv_usec;
		sconf->t_update = 0;

		// TLP c135 leyland et. al testcase
		// 116915434112329921568236283928181979297762987646390347857868153872054154807376462439621333455331738807075404918922573575454310187518221

        // if we are using TLP, this is where we periodically filter to
        // test cycle formation.
		if (sconf->use_dlp >= 2)
		{
	        uint32_t total_rels = sconf->num_full + sconf->num_slp +
				sconf->dlp_useful + sconf->tlp_useful + sconf->qlp_useful;

            // use the number of full relations found as a gauge for
            // when to filter.  At the end of a TLP run with typical parameters we
            // have about 20% fulls and 80% cycles.
            if ((sconf->num_full > check_total) &&
                (sconf->num_found < (sconf->factor_base->B + sconf->num_extra_relations)))
			{
				fact_obj_t *obj = sconf->obj;
				siqs_r *relation_list;
				qs_la_col_t *cycle_list;
				uint32_t num_cycles;
				uint32_t *hashtable = sconf->cycle_hashtable;
				qs_cycle_t *table = sconf->cycle_table;
				uint32_t num_relations;
				uint32_t i = 0, passes;
				uint32_t curr_a_idx, curr_poly_idx, curr_rel;
				uint32_t curr_expected, curr_saved, curr_cycle;
				uint32_t all_relations;
				uint32_t total_poly_a = 0;
				uint32_t poly_saved;
				uint32_t *plist0;
				uint32_t *plist1;
				uint32_t *plist2;
                uint32_t* apolylist = NULL;
				uint32_t newrels;
				int j;
				char buf[LINE_BUF_SIZE];
				struct timeval filt_start, filt_stop;
                double tmp;

				gettimeofday(&filt_start, NULL);

                printf("\nreached deadline of %u (%u) full relations found\n",
                    sconf->num_full, check_total);
				t_time = ytools_difftime(&sconf->totaltime_start, &filt_start);
				printf("QS elasped time is now %1.2f sec\n", t_time);
				printf("reading relations\n");

				if (obj->logfile != NULL)
				{
                    printf("\nreached deadline of %u (%u) full relations found\n",
                        sconf->num_full, check_total);
					logprint(obj->logfile, "QS elasped time is now %1.2f sec\n", t_time);
				}

                if (sconf->in_mem)
                {
                    relation_list = (siqs_r*)xmalloc(sconf->buffered_rels * sizeof(siqs_r));
                    plist0 = (uint32_t*)xmalloc(sconf->buffered_rels * sizeof(uint32_t));
                    plist1 = (uint32_t*)xmalloc(sconf->buffered_rels * sizeof(uint32_t));
                    plist2 = (uint32_t*)xmalloc(sconf->buffered_rels * sizeof(uint32_t));
                    apolylist = (uint32_t*)xmalloc(sconf->buffered_rels * sizeof(uint32_t));
                    for (i = 0; i < sconf->buffered_rels; i++)
                    {
                        relation_list[i].poly_idx = i;
                        relation_list[i].apoly_idx = apolylist[i] = sconf->in_mem_relations[i].apoly_idx;
                        relation_list[i].large_prime[0] = plist0[i] = sconf->in_mem_relations[i].large_prime[0];
                        relation_list[i].large_prime[1] = plist1[i] = sconf->in_mem_relations[i].large_prime[1];
                        relation_list[i].large_prime[2] = plist2[i] = sconf->in_mem_relations[i].large_prime[2];
                    }
                    total_poly_a = sconf->total_poly_a;
                    all_relations = num_relations = sconf->buffered_rels;
                }
                else
                {
                    /* skip over the first line */
                    savefile_flush(&sconf->obj->qs_obj.savefile);
                    savefile_close(&sconf->obj->qs_obj.savefile);

                    savefile_open(&sconf->obj->qs_obj.savefile, SAVEFILE_READ);
                    savefile_read_line(buf, sizeof(buf), &sconf->obj->qs_obj.savefile);

                    // we don't know beforehand how many rels to expect, so start
                    // with some amount and allow it to increase as we read them
                    relation_list = (siqs_r*)xmalloc(10000 * sizeof(siqs_r));
                    plist0 = (uint32_t*)xmalloc(10000 * sizeof(uint32_t));
                    plist1 = (uint32_t*)xmalloc(10000 * sizeof(uint32_t));
                    plist2 = (uint32_t*)xmalloc(10000 * sizeof(uint32_t));
                    curr_rel = 10000;

                    while (!savefile_eof(&sconf->obj->qs_obj.savefile)) {
                        char* start;

                        switch (buf[0]) {
                        case 'A':
                            total_poly_a++;
                            break;

                        case 'R':
                            start = strchr(buf, 'L');
                            if (start != NULL) {
                                uint32_t primes[MAXLP];
                                int numlp = yafu_read_Nlp(start, primes);
                                if (i == curr_rel) {
                                    curr_rel = 3 * curr_rel / 2;
                                    relation_list = (siqs_r*)xrealloc(
                                        relation_list,
                                        curr_rel *
                                        sizeof(siqs_r));

                                    plist0 = (uint32_t*)xrealloc(plist0, curr_rel * sizeof(uint32_t));
                                    plist1 = (uint32_t*)xrealloc(plist1, curr_rel * sizeof(uint32_t));
                                    plist2 = (uint32_t*)xrealloc(plist2, curr_rel * sizeof(uint32_t));
                                }

                                //printf("found primes %u,%u,%u on line %u, current allocation: %u\n",
                                //	primes[0], primes[1], primes[2], i, curr_rel);

                                relation_list[i].poly_idx = i;
                                relation_list[i].large_prime[0] = primes[0];
                                relation_list[i].large_prime[1] = primes[1];
                                relation_list[i].large_prime[2] = primes[2];

                                plist0[i] = primes[0];
                                plist1[i] = primes[1];
                                plist2[i] = primes[2];
                                i++;
                            }
                            break;
                        }

                        savefile_read_line(buf, sizeof(buf), &sconf->obj->qs_obj.savefile);
                    }

                    savefile_flush(&sconf->obj->qs_obj.savefile);
                    savefile_close(&sconf->obj->qs_obj.savefile);
                    // get ready to collect more relations
                    savefile_open(&sconf->obj->qs_obj.savefile, SAVEFILE_APPEND);

                    all_relations = i;
                    num_relations = i;
                }

                gettimeofday(&filt_stop, NULL);
                tmp = ytools_difftime(&filt_start, &filt_stop);

				printf("read %u relations in %1.2f seconds\n", i, tmp);
				printf("building graph\n");

				if (obj->logfile != NULL)
				{
					logprint(obj->logfile, "building graph with %u relations\n", i);
				}

				newrels = num_relations - sconf->last_numcycles;
				sconf->last_numcycles = num_relations;
				rebuild_graph3(sconf, relation_list, num_relations);

				hashtable = sconf->cycle_hashtable;
				table = sconf->cycle_table;

				printf("commencing %dlp singleton removal\n", sconf->num_lp);
				//printf("cycle table size = %u, vertices = %u, components = %u\n",
				//	sconf->cycle_table_size, sconf->vertices, sconf->components);

				if (obj->logfile != NULL)
				{
					logprint(obj->logfile, "commencing %dlp singleton removal\n", sconf->num_lp);
					//logprint(obj->logfile, "cycle table size = %u, vertices = %u, components = %u\n",
					//	sconf->cycle_table_size, sconf->vertices, sconf->components);
				}

				num_relations = qs_purge_singletonsN(sconf->obj, relation_list, 
					num_relations, sconf->num_lp, table, hashtable);

				printf("%u relations survived singleton removal\n", num_relations);

				if (obj->logfile != NULL)
				{
					logprint(obj->logfile, "%u relations survived singleton removal\n", num_relations);
				}

				if (num_relations > 0)
				{
					double fullrate;
					double cycrate;
					//printf("commencing duplicate removal\n");
					//
					//num_relations = qs_purge_duplicate_relations3(sconf->obj,
					//	relation_list, num_relations);
					//
					//printf("%u relations survived duplicate removal\n", num_relations);

					printf("building reduced graph\n");

					if (obj->logfile != NULL)
					{
						logprint(obj->logfile, "building reduced graph\n");
					}

					rebuild_graph3(sconf, relation_list, num_relations);

					//printf("cycle table size = %u, edges = %u, vertices = %u, components = %u\n",
					//	sconf->cycle_table_size, sconf->num_cycles, 
					//	sconf->vertices, sconf->components);

					printf("commencing cycle find\n");

					if (obj->logfile != NULL)
					{
						//logprint(obj->logfile, "cycle table size = %u, edges = %u, vertices = %u, components = %u\n",
						//	sconf->cycle_table_size, sconf->num_cycles,
						//	sconf->vertices, sconf->components);
						logprint(obj->logfile, "commencing cycle find\n");
					}

					// before cycle find, num_cycles is the number of vertices 
					// (unique primes) and num_relations is all relations found.
					// after cycle find, sconf->num_r is the actual number of cycles
					// formed (full relations and cycles from partials),
					// sconf->num_relations is the number of fulls, and sconf->num_cycles 
					// is a count of all raw relations found.
					// the variable num_cycles is equal to num_r and the variable
					// num_relations is equal to sconf->num_cycles (all raw relations found).
					// this is very confusing and needs to get cleaned up.
					num_cycles = sconf->vertices;
					if (num_cycles > num_relations)
					{
						printf("not enough relations for cycle find\n");
						printf("restoring full %u relation set\n", all_relations);

						for (j = 0; j < all_relations; j++)
						{
							relation_list[j].large_prime[0] = plist0[j];
							relation_list[j].large_prime[1] = plist1[j];
							relation_list[j].large_prime[2] = plist2[j];
						}
						num_relations = all_relations;

						printf("rebuilding full graph\n");

						rebuild_graph3(sconf, relation_list, num_relations);
						sconf->num_r = 0;
						sconf->num_relations = 0;
						sconf->num_cycles = num_relations;
					}
					else
					{
						double next_cycrate;
                        uint32_t new_full, new_partial;

						cycle_list = find_cycles3(sconf->obj, sconf, relation_list,
							num_relations, &num_cycles, &passes);

                        new_full = sconf->num_relations - sconf->last_numfull;
                        new_partial = sconf->num_r - sconf->num_relations - sconf->last_numpartial;

						printf("added %u full and %u cycles from partials with %u "
							"new relations since last update\n",
                            new_full, new_partial, newrels);

						if (obj->logfile != NULL)
						{
							logprint(obj->logfile, "added %u full and %u cycles from partials with %u "
								"new relations since last update\n",
								sconf->num_relations - sconf->last_numfull,
								sconf->num_r - sconf->num_relations - sconf->last_numpartial, newrels);
						}

						fullrate = (double)(sconf->num_relations - sconf->last_numfull) /
							(double)(newrels);
						cycrate = (double)(sconf->num_r - sconf->num_relations - sconf->last_numpartial) /
							(double)(newrels);

						printf("cycle gathering rate is now ~%1.4f fulls/rel and ~%1.4f cycles/rel\n",
							fullrate, cycrate);

						if (obj->logfile != NULL)
						{
							logprint(obj->logfile, "cycle gathering rate is now ~%1.4f fulls/rel and ~%1.4f cycles/rel\n",
								fullrate, cycrate);
						}

						sconf->last_numfull = sconf->num_relations;
						sconf->last_numpartial = sconf->num_r - sconf->num_relations;

                        int have_cycle_growth_estimate = 0;
                        if (sconf->last_cycrate > 0)
                        {
                            uint32_t current_rels = sconf->num_r;
                            uint32_t needed_rels = (fb->B - sconf->num_r);
                            uint32_t batch_rels;
                            uint32_t raw_rels = 0;
                            double fr = fullrate;
                            double cr = cycrate;
                            double lr = sconf->last_cycrate;
                            double ri = (cr - lr);
                            uint32_t ci = check_inc;
                            uint32_t fulls_at_current_rate = 0;
                            uint32_t fulls_at_linear_rate = 0;

                            lr = cr;
                            cr = cr + ri;

                            fulls_at_current_rate = (uint32_t)(((double)needed_rels /
                                (fullrate + cycrate)) * fullrate);

                            printf("at the current cycle gathering rate we need %u "
                                "more raw relations (%u more fulls)\n",
                                (uint32_t)((double)needed_rels / (fullrate + cycrate)),
                                fulls_at_current_rate);

                            printf("cycle gathering rate increased by %1.4f "
                                "since last filter\n", ri);

                            while (current_rels < fb->B)
                            {
                                uint32_t new_fulls = ci;
                                uint32_t new_raws = ci / fr;
                                double new_cycle_rate = cr + ri;
                                uint32_t new_cycles = new_cycle_rate * new_raws;
                                current_rels += (new_cycles + new_fulls);
                                cr = new_cycle_rate;
                                raw_rels += new_raws;
                                ci = ci / 2;
                                if (ci < 16) ci = 16;
                            }

                            fulls_at_linear_rate = (uint32_t)(raw_rels * fr);

                            if (fulls_at_linear_rate > fulls_at_current_rate)
                            {
                                // this case occurs when the current increment (ci)
                                // is enough to get to the goal.  The current strategy:
                                // (fulls_at_current_rate + fulls_at_linear_rate) / 3
                                // actually seems to work pretty well, but maybe I've 
                                // been getting lucky in test cases.  Leave it for now.
                                // just don't print the linear growth case since it doesn't
                                // make much sense as-is.
                            }
                            else
                            {
                                printf("with linear growth of cycle rate we need %u "
                                    "more raw relations (%u more fulls)\n", raw_rels,
                                    fulls_at_linear_rate);
                            }

                            if ((fulls_at_current_rate < fb->B) &&
                                (fulls_at_linear_rate < fb->B))
                            {
                                have_cycle_growth_estimate = 1;
                                sconf->check_inc = (fulls_at_current_rate + 
                                    fulls_at_linear_rate) / 4;
                            }
                        }
                        
                        // use the number of full relations found as a gauge for
                        // when to filter.  At the end of a TLP run with typical parameters we
                        // have about 20% fulls and 80% cycles.  So we start filtering at 15% fulls
                        // and every 1% thereafter.
                        double percent_complete = 
                            ((double)sconf->last_numfull + (double)sconf->last_numpartial) / 
                            (double)sconf->factor_base->B;

                        // check frequently for very large inputs.  Sieving will take
                        // so long that filtering costs are negligible, and good
                        // data will emerge from the filtering runs that might hone in
                        // cycle formation predictions.
                        sconf->check_inc = 1000;

						// this is a better approximation than just assuming that 
						// cycle formation stays constant.
                        next_cycrate = cycrate + (cycrate - sconf->last_cycrate) / 1.75;

						sconf->last_fullrate = fullrate;
						sconf->last_cycrate = cycrate;

                        printf("predicting %u new fulls and %u new cycles in next batch of ~%u rels\n",
                            (uint32_t)((double)sconf->check_inc),
                            (uint32_t)((double)sconf->check_inc / fullrate * next_cycrate),
                            (uint32_t)((double)sconf->check_inc / fullrate));

						if (sconf->num_r >= (fb->B + sconf->num_extra_relations))
						{
							//we've got enough total relations to stop
							retcode = 1;
						}
						else
						{
							// as we get close enough that the next batch could possibly
							// be the last one, use the prediction to scale the batch
							// size better.
							//if (((uint32_t)((double)sconf->check_inc * fullrate) +
							//	(uint32_t)((double)sconf->check_inc * next_cycrate) +
							//	sconf->num_r) > (fb->B + sconf->num_extra_relations))
                            uint32_t predicted_rels = (uint32_t)((double)sconf->check_inc) +
                                (uint32_t)((double)sconf->check_inc / fullrate * cycrate);

                            if ((predicted_rels + sconf->num_r) > 
                                (fb->B + sconf->num_extra_relations))
							{
								uint32_t rels_needed = (fb->B + sconf->num_extra_relations) - sconf->num_r;
                                double batch_fraction = (double)rels_needed / (double)predicted_rels;
                                sconf->check_inc = (uint32_t)((double)sconf->check_inc * batch_fraction);
								printf("using cycle formation rate to set batch size to %u\n",
									sconf->check_inc);
							}

                            sconf->check_inc = MIN(sconf->check_inc, 1000);

							for (j = 0; j < num_cycles; j++)
							{
								free(cycle_list[j].cycle.list);
							}
							free(cycle_list);

							printf("found %u cycles in %u passes\n", num_cycles, passes);

							printf("restoring full %u relation set\n", all_relations);

							for (j = 0; j < all_relations; j++)
							{
                                if (apolylist != NULL)
                                {
                                    relation_list[j].apoly_idx = apolylist[j];
                                }
								relation_list[j].large_prime[0] = plist0[j];
								relation_list[j].large_prime[1] = plist1[j];
								relation_list[j].large_prime[2] = plist2[j];
							}
							num_relations = all_relations;

							printf("rebuilding full graph\n");

							rebuild_graph3(sconf, relation_list, num_relations);
						}
					}
				}
				else
				{
					printf("restoring full %u relation set\n", all_relations);

					for (j = 0; j < all_relations; j++)
					{
						relation_list[j].large_prime[0] = plist0[j];
						relation_list[j].large_prime[1] = plist1[j];
						relation_list[j].large_prime[2] = plist2[j];
					}
					num_relations = all_relations;

					printf("rebuilding full graph\n");

					rebuild_graph3(sconf, relation_list, num_relations);
					sconf->num_r = 0;
					sconf->num_relations = 0;
					sconf->num_cycles = num_relations;
				}

				free(plist0);
				free(plist1);
				free(plist2);
				free(relation_list);
                if (apolylist != NULL)
                {
                    free(apolylist);
                }

				sconf->check_total += sconf->check_inc;
                printf("new filtering deadline is %u full relations\n", sconf->check_total);

				gettimeofday(&filt_stop, NULL);
				sconf->t_time4 += ytools_difftime(&filt_start, &filt_stop);
			}
            else if (sconf->num_found >= (sconf->factor_base->B + sconf->num_extra_relations))
            {
                // we've got enough total relations to stop
                retcode = 1;
            }
            
		}
		else
		{

			if (sconf->num_r >= (fb->B + sconf->num_extra_relations))
			{
				//we've got enough total relations to stop
				mpz_clear(tmp1);
				return 1;
			}
			else
			{
				// need to keep sieving.  since the last time we checked, we've found
				// (full->num_r + partial->act_r) - (last_numfull + last_numpartial)
				// relations.  assume we'll find this many next time, and 
				// scale how much longer we need to sieve to hit the target.

				sconf->num_expected = sconf->num_r - sconf->last_numfull - sconf->last_numpartial;
				// if the number expected to be found next time puts us over the needed amount, scale the 
				// check total appropriately.  otherwise, just increment check_total by check_inc
				if ((sconf->num_r + sconf->num_expected) > (fb->B + sconf->num_extra_relations))
				{
					sconf->num_needed = (fb->B + sconf->num_extra_relations) - sconf->num_r;
					sconf->check_total +=
						(uint32_t)((double)check_inc * (double)sconf->num_needed / (double)sconf->num_expected);
					sconf->check_total += sconf->num_extra_relations;
					sconf->update_time *= (double)sconf->num_needed / (double)sconf->num_expected;

					// always go at least one more second.
					if (sconf->update_time < 1)
						sconf->update_time = 1;
				}
				else
					sconf->check_total += sconf->check_inc;
			}
			sconf->last_numfull = num_full;
			sconf->last_numcycles = sconf->num_cycles;
			sconf->last_numpartial = sconf->num_cycles + sconf->components - sconf->vertices;
		}
	}

	mpz_clear(tmp1);
	return retcode;
}

int update_final(static_conf_t *sconf)
{
	FILE *sieve_log = sconf->obj->logfile;
	mpz_t tmp1;
	struct timeval myTVend;
    double t_time;

	mpz_init(tmp1);
    gettimeofday(&myTVend, NULL);
    t_time = ytools_difftime(&sconf->totaltime_start, &myTVend);

	if (sconf->obj->VFLAG >= 0)
	{
        struct timeval update_stop;
        gettimeofday(&update_stop, NULL);
        t_time = ytools_difftime(&sconf->totaltime_start, &update_stop);

		if (sconf->use_dlp >= 2)
		{
			printf("last: %u, now: %u full, %u slp, "
				"%u dlp (%ukatt, %ukprp), "
				"%u tlp (%ukatt, %ukprp), (%1.0f r/sec)\r",
				sconf->last_numfull + sconf->last_numpartial, sconf->num_full, sconf->num_slp,
				sconf->dlp_useful, sconf->attempted_squfof / 1000,
				sconf->dlp_prp / 1000, sconf->tlp_useful, sconf->attempted_cosiqs / 1000,
				sconf->tlp_prp / 1000,
				(double)(sconf->num_relations + sconf->dlp_useful +
					sconf->tlp_useful) / t_time);
		}
		else
		{
			printf("%d rels found: %d full + "
				"%d from %d partial, (%6.2f rels/sec)\r",
				sconf->num_r, sconf->num_relations,
				sconf->num_cycles +
				sconf->components - sconf->vertices,
				sconf->num_cycles,
				(double)(sconf->num_relations + sconf->num_cycles) / t_time);
		}

		mpz_set_ui(tmp1, sconf->tot_poly);					//total number of polys
		mpz_mul_ui(tmp1, tmp1, sconf->num_blocks);	//number of blocks
		mpz_mul_2exp(tmp1, tmp1, 1);			//pos and neg sides
		mpz_mul_2exp(tmp1, tmp1, sconf->qs_blockbits);	//sieve locations per block	

		if (sconf->obj->VFLAG > 0)
		{
			printf("\n\nsieving required %u total polynomials (%u 'A' polynomials)\n",
				sconf->tot_poly, sconf->total_poly_a);
			gmp_printf("trial division touched %lu sieve locations out of %Zd\n",
				sconf->num, tmp1);

			if (sconf->use_dlp >= 2)
			{
				printf("tlp-ecm: %u failures, %u attempts, %u outside range, %u prp, %u useful\n",
					sconf->failed_cosiqs, sconf->attempted_cosiqs,
					sconf->tlp_outside_range, sconf->tlp_prp, sconf->tlp_useful);
			}

			if (sconf->use_dlp)
			{
				printf("dlp-ecm: %u failures, %u attempts, %u outside range, %u prp, %u useful\n",
					sconf->failed_squfof, sconf->attempted_squfof,
					sconf->dlp_outside_range, sconf->dlp_prp, sconf->dlp_useful);
			}

            printf("total reports = %lu, total surviving reports = %lu\ntotal blocks sieved = %lu, "
                "avg surviving reports per block = %1.2f\n", sconf->total_reports,
                sconf->total_surviving_reports, sconf->total_blocks,
                (float)sconf->total_surviving_reports / (float)sconf->total_blocks);
            printf("Elapsed time: %1.4f sec\n", t_time);

            //printf("scatter/gather processing: %lu of %lu opportunities\n",
            //    sconf->num_scatter, sconf->num_scatter_opp);

#ifdef USE_AVX2
            if (sconf->lp_scan_failures > 0)
            {
                printf("large prime scan failures = %u\n", sconf->lp_scan_failures);
            }
#endif

		}
		else
		{
			printf("\n\n");
		}

		if (sieve_log != NULL)
		{
            char* s = mpz_get_str(NULL, 10, tmp1);
			logprint(sieve_log, "trial division touched %u sieve locations out of %s\n",
				sconf->num, s);
            free(s);

			logprint(sieve_log, "total reports = %lu, total surviving reports = %lu\n", sconf->total_reports,
				sconf->total_surviving_reports);

			logprint(sieve_log, "total blocks sieved = %lu, avg surviving reports per block = %1.2f\n", 
				sconf->total_blocks, (float)sconf->total_surviving_reports / (float)sconf->total_blocks);
			
			if (sconf->use_dlp >= 2)
			{
				logprint(sieve_log, "tlp-ecm: %u failures, %u attempts, %u outside range, %u prp, %u useful\n",
					sconf->failed_cosiqs, sconf->attempted_cosiqs,
					sconf->tlp_outside_range, sconf->tlp_prp, sconf->tlp_useful);
			}

			if (sconf->use_dlp)
			{
				logprint(sieve_log, "dlp-ecm: %u failures, %u attempts, %u outside range, %u prp, %u useful\n",
					sconf->failed_squfof, sconf->attempted_squfof,
					sconf->dlp_outside_range, sconf->dlp_prp, sconf->dlp_useful);
			}
		}

		fflush(stdout);
		fflush(stderr);
	}

	if (sieve_log != NULL)
	{
		if (sconf->use_dlp >= 2)
		{
			uint32_t total_rels = sconf->num_full + sconf->num_slp +
				sconf->dlp_useful + sconf->tlp_useful;

			logprint(sieve_log, "%d relations found: %d full + "
				"%d from %d partial, using %d polys (%d A polys)\n",
				sconf->num_r, sconf->num_relations,
				sconf->num_r - sconf->num_relations,
				sconf->num_slp + sconf->dlp_useful + sconf->tlp_useful, 
				sconf->tot_poly, sconf->total_poly_a);
			logprint(sieve_log, "on average, sieving found %1.2f raw rels/poly and %1.2f raw rels/sec\n",
				(double)(total_rels) / (double)sconf->tot_poly,
				(double)(total_rels) / t_time);
		}
		else
		{
			logprint(sieve_log, "%d relations found: %d full + "
				"%d from %d partial, using %d polys (%d A polys)\n",
				sconf->num_r, sconf->num_relations,
				sconf->num_cycles +
				sconf->components - sconf->vertices,
				sconf->num_cycles, sconf->tot_poly, sconf->total_poly_a);
			logprint(sieve_log, "on average, sieving found %1.2f rels/poly and %1.2f rels/sec\n",
				(double)(sconf->num_relations + sconf->num_cycles) / (double)sconf->tot_poly,
				(double)(sconf->num_relations + sconf->num_cycles) / t_time);
		}

        char* s = mpz_get_str(NULL, 10, tmp1);
		logprint(sieve_log, "trial division touched %lu sieve locations out of %s\n",
			sconf->num, s);
        free(s);
		logprint(sieve_log,"==== post processing stage (msieve-1.38) ====\n");
	}

	sconf->obj->qs_obj.rels_per_sec =
        (double)(sconf->num_relations + sconf->num_cycles) / t_time;

	if (sieve_log != NULL)
	{
		fflush(sieve_log);
	}

	//sieve_log = fopen(flogname,"a");
	mpz_clear(tmp1);

	return 0;
}

int free_sieve(dynamic_conf_t* dconf)
{
    uint32_t i;

    //can free sieving structures now
#if defined(USE_SS_SEARCH)
    if (dconf->using_ss_search)
    {
        align_free(dconf->ss_sieve_n);
        align_free(dconf->ss_sieve_p);
    }
    else
    {
        dconf->sieve = dconf->sieve - 2 * 32768;
        align_free(dconf->sieve);
    }
    if (dconf->using_ss_search)
    {
        free(dconf->firstroot1a);
        free(dconf->firstroot1b);
        free(dconf->firstroot2);
        mpz_clear(dconf->polyb1);
        mpz_clear(dconf->polyb2);
        free(dconf->polyv);
        free(dconf->polysign);
        free(dconf->polynums);
        free(dconf->polymap);
    }
#ifdef USE_DIRECT_SIEVE_SS
    free(dconf->report_ht_n);
    free(dconf->report_ht_p);
#endif
#ifdef USE_POLY_BUCKET_SS
    if (dconf->using_ss_search)
    {
        for (i = 0; i < dconf->num_ss_slices; i++)
        {
            free(dconf->ss_slices_p[i].elements);
            free(dconf->ss_slices_p[i].size);
        }

        free(dconf->ss_slices_p);
        free(dconf->ss_slices_n);
    }
#endif
#else
    dconf->sieve = dconf->sieve - 2 * 32768;
    align_free(dconf->sieve);
#endif
	align_free(dconf->comp_sieve_p->prime);
	align_free(dconf->comp_sieve_p->root1);
	align_free(dconf->comp_sieve_p->root2);
	align_free(dconf->comp_sieve_p->logp);
	align_free(dconf->comp_sieve_n->prime);
	align_free(dconf->comp_sieve_n->root1);
	align_free(dconf->comp_sieve_n->root2);
	align_free(dconf->comp_sieve_n->logp);
	free(dconf->comp_sieve_p);
	free(dconf->comp_sieve_n);

	align_free(dconf->rootupdates);
	align_free(dconf->sm_rootupdates);

	align_free(dconf->update_data.sm_firstroots1);
	align_free(dconf->update_data.sm_firstroots2);
	align_free(dconf->update_data.firstroots1);
	align_free(dconf->update_data.firstroots2);
	align_free(dconf->update_data.prime);
	align_free(dconf->update_data.logp);

	monty_free(dconf->mdata);
    free(dconf->mdata);

	if (dconf->buckets->list != NULL)
	{
		align_free(dconf->buckets->list);
		free(dconf->buckets->fb_bounds);
		free(dconf->buckets->logp);
		align_free(dconf->buckets->num);
        free(dconf->buckets->num_slices_batch);
		free(dconf->buckets);
	}
    else
    {
        free(dconf->buckets);
    }

	//support data on the poly currently being sieved
	free(dconf->curr_poly->gray);
	free(dconf->curr_poly->nu);
	free(dconf->curr_poly->qlisort);
	mpz_clear(dconf->curr_poly->mpz_poly_a);
	mpz_clear(dconf->curr_poly->mpz_poly_b);
	mpz_clear(dconf->curr_poly->mpz_poly_c);
	free(dconf->curr_poly);

	for (i=0;i<MAX_A_FACTORS;i++)
		mpz_clear(dconf->Bl[i]);
	free(dconf->Bl);

	//workspace bigints
	mpz_clear(dconf->gmptmp1);
	mpz_clear(dconf->gmptmp2);
	mpz_clear(dconf->gmptmp3);
	
	align_free(dconf->mask);
    align_free(dconf->mask2);

	//free sieve scan report stuff
	free(dconf->reports);
	
	for (i=0; i<MAX_SIEVE_REPORTS; i++)
	{
		mpz_clear(dconf->Qvals[i]);
	}
	free(dconf->Qvals);	
	free(dconf->valid_Qs);
	free(dconf->smooth_num);

    if (dconf->do_batch)
    {
        relation_batch_free(&dconf->rb);
    }

#ifdef USE_8X_MOD_ASM
	align_free(dconf->bl_locs);
	align_free(dconf->bl_sizes);
#endif
	align_free(dconf->corrections);
    align_free(dconf->polyscratch);

	return 0;
}

void free_filter_vars(static_conf_t *sconf)
{
	int i;

	// cycle table stuff created at the beginning of the factorization
	free(sconf->cycle_hashtable);
	free(sconf->cycle_table);

	if (sconf->poly_list != NULL)
	{
		for (i=0;(uint32_t)i < sconf->poly_list_alloc;i++)
			mpz_clear(sconf->poly_list[i].b);
		free(sconf->poly_list);
	}

	if (sconf->relation_list != NULL)
	{
		// free post-processed relations
		free(sconf->relation_list);
	}

    if (sconf->cycle_list != NULL)
    {
        for (i = 0; i < sconf->num_cycles; i++)
        {
            if (&sconf->cycle_list[i] != NULL)
            {
                if (sconf->cycle_list[i].cycle.list != NULL)
                    free(sconf->cycle_list[i].cycle.list);
                if (sconf->cycle_list[i].data != NULL)
                    free(sconf->cycle_list[i].data);
            }
        }
        free(sconf->cycle_list);
    }

	return;
}

int free_siqs(static_conf_t *sconf)
{
	uint32_t i;

	//current poly info used during filtering
	free(sconf->curr_poly->gray);
	free(sconf->curr_poly->nu);
	free(sconf->curr_poly->qlisort);
	mpz_clear(sconf->curr_poly->mpz_poly_a);
	mpz_clear(sconf->curr_poly->mpz_poly_b);
	mpz_clear(sconf->curr_poly->mpz_poly_c);
	free(sconf->curr_poly);
	mpz_clear(sconf->curr_a);	
    align_free(sconf->modsqrt_array);
	align_free(sconf->factor_base->list->prime);
    align_free(sconf->factor_base->list->binv);
	align_free(sconf->factor_base->list->small_inv);
	align_free(sconf->factor_base->list->correction);
	align_free(sconf->factor_base->list->logprime);
	align_free(sconf->factor_base->tinylist->prime);
	align_free(sconf->factor_base->tinylist->small_inv);
	align_free(sconf->factor_base->tinylist->correction);
	align_free(sconf->factor_base->tinylist->logprime);
	align_free(sconf->factor_base->list);
	align_free(sconf->factor_base->tinylist);
	free(sconf->factor_base);
    align_free(sconf->sieve_primes);

    // list of a values used first to track all a coefficients
    // generated during sieving.
    for (i = 0; (uint32_t)i < sconf->total_poly_a; i++)
    {
        mpz_clear(sconf->poly_a_list[i]);
    }
    free(sconf->poly_a_list);

	//while freeing the list of factors, divide them out of the input
	for (i=0;i<sconf->factor_list.num_factors;i++)
	{
		mpz_t tmp;
		mpz_init(tmp);

		//convert the factor
		mp_t2gmp(&sconf->factor_list.final_factors[i]->factor,tmp);

		//divide it out
		mpz_tdiv_q(sconf->obj->qs_obj.gmp_n, sconf->obj->qs_obj.gmp_n, tmp);

		//log it
		add_to_factor_list(sconf->obj->factors, tmp, 
            sconf->obj->VFLAG, sconf->obj->NUM_WITNESSES);
		
		mpz_clear(tmp);
		free(sconf->factor_list.final_factors[i]);
	}

    free(sconf->in_mem_relations);

	mpz_clear(sconf->sqrt_n);
	mpz_clear(sconf->n);
	mpz_clear(sconf->target_a);

	//free(sconf->obj->savefile.name);
	savefile_free(&sconf->obj->qs_obj.savefile);
    

	return 0;
}

// used in multiplier selection
#define NUM_TEST_PRIMES 300
#define NUM_MULTIPLIERS (sizeof(mult_list)/sizeof(uint8_t))

static const uint8_t mult_list_orig[] =
{ 1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19,
 21, 22, 23, 26, 29, 30, 31, 33, 34, 35, 37, 38,
 39, 41, 42, 43, 46, 47, 51, 53, 55, 57, 58, 59,
 61, 62, 65, 66, 67, 69, 70, 71, 73 };

// more options.  Thanks Till for pointing out that my original list was too small!
static const uint8_t mult_list_many_even[] = {
1  ,   2  ,   3  ,   5  ,   6  ,   7  ,   9  ,  10  ,  11  ,  13  ,  14  ,
15 ,   17 ,   19 ,   21 ,   22 ,   23 ,   25 ,   26 ,   29 ,   30 ,   31 ,
33 ,   34 ,   35 ,   37 ,   38 ,   39 ,   41 ,   42 ,   43 ,   45 ,   46 ,
47 ,   49 ,   51 ,   53 ,   55 ,   57 ,   58 ,   59 ,   61 ,   62 ,   63 ,
65 ,   66 ,   67 ,   69 ,   70 ,   71 ,   73 ,   75 ,   77 ,   79 ,   83 ,
85 ,   87 ,   89 ,   91 ,   93 ,   95 ,   97 ,  101 ,  103 ,  105 ,  107 ,
109,   111,   113,   115,   119,   121,   123,   127,   129,   131,   133,
137,   139,   141,   143,   145,   147,   149,   151,   155,   157,   159,
161,   163,   165,   167,   173,   177,   179,   181,   183,   185,   187,
191,   193,   195,   197,   199,   201,   203,   205,   209,   211,   213,
215,   217,   219,   223,   227,   229,   231,   233,   235,   237,   239,
241,   249,   251,   253,   255 };

static const uint8_t mult_list[] = {
1  ,   2  ,   3  ,   5  ,   7  ,   9  ,   10 ,   11 ,   13 ,   14  ,
15 ,   17 ,   19 ,   21 ,   23 ,   25 ,   29 ,   31 ,
33 ,   35 ,   37 ,   39 ,   41 ,   43 ,   45 ,   
47 ,   49 ,   51 ,   53 ,   55 ,   57 ,   59 ,   61 ,   63 ,
65 ,   67 ,   69 ,   71 ,   73 ,   75 ,   77 ,   79 ,   83 ,
85 ,   87 ,   89 ,   91 ,   93 ,   95 ,   97 ,   101,   103,   105,   107,
109,   111,   113,   115,   119,   121,   123,   127,   129,   131,   133,
137,   139,   141,   143,   145,   147,   149,   151,   155,   157,   159,
161,   163,   165,   167,   173,   177,   179,   181,   183,   185,   187,
191,   193,   195,   197,   199,   201,   203,   205,   209,   211,   213,
215,   217,   219,   223,   227,   229,   231,   233,   235,   237,   239,
241,   249,   251,   253,   255 };

uint8_t choose_multiplier_siqs(uint32_t B, mpz_t n) 
{
	uint32_t i, j;
	uint32_t num_primes = MIN(2 * B, NUM_TEST_PRIMES);
	double best_score;
	uint8_t best_mult;
	double scores[NUM_MULTIPLIERS];
	uint32_t num_multipliers;
	double log2n = zlog(n);

	/* measure the contribution of 2 as a factor of sieve
	   values. The multiplier itself must also be taken into
	   account in the score. scores[i] is the correction that
	   is implicitly applied to the size of sieve values for
	   multiplier i; a negative score makes sieve values 
	   smaller, and so is better */

	for (i = 0; i < NUM_MULTIPLIERS; i++) {
		uint8_t curr_mult = mult_list[i];
		uint8_t knmod8 = (uint8_t)((curr_mult * mpz_get_ui(n)) % 8);
		double logmult = log((double)curr_mult);

		/* only consider multipliers k such than
		   k*n will not overflow an mp_t */

        // even k seem to be worthwhile, on average.  only occasionally 
        // are they a bit slower.
        //if ((curr_mult & 1) == 0)
        //{
        //    scores[i] = 1000;
        //    continue;
        //}

		if (log2n + logmult > (32 * MAX_DIGITS - 2) * LN2)
			break;

		scores[i] = 0.5 * logmult;
		switch (knmod8) {
		case 1:
            // when kn == 1 mod 8 then we use Q2(x) polys that
            // result in sieve locations being 1 bit smaller, so
            // we use 2.625 * log(2) as the contribution of 2.
            // the value 2.625 got 9 out of 10 borderline cases correct
            // in an ad-hoc test but could probably use more testing.
            scores[i] -= (2.625 * LN2); // (2 * LN2);
			break;
		case 5:
			scores[i] -= LN2;
			break;
		case 3:
		case 7:
			scores[i] -= 0.5 * LN2;
			break;
		/* even multipliers start with a handicap */
		}
	}
	num_multipliers = i;

	/* for the rest of the small factor base primes */

	for (i = 1; i < num_primes; i++) {
		uint32_t prime = (uint32_t)siqs_primes[i];
		double contrib = log((double)prime) / (prime - 1);
		uint32_t modp = (uint32_t)mpz_tdiv_ui(n, prime);

		for (j = 0; j < num_multipliers; j++) {
			uint8_t curr_mult = mult_list[j];
			//uint32_t knmodp = mp_modmul_1(modp, curr_mult, prime);
			uint32_t knmodp = (modp * curr_mult) % prime;

			/* if prime i is actually in the factor base
			   for k * n ... */

			if (knmodp == 0 || jacobi_1(knmodp, prime) == 1) {

				/* ...add its contribution. A prime p con-
				   tributes log(p) to 1 in p sieve values, plus
				   log(p) to 1 in p^2 sieve values, etc. The
				   average contribution of all multiples of p 
				   to a random sieve value is thus

				   log(p) * (1/p + 1/p^2 + 1/p^3 + ...)
				   = (log(p) / p) * 1 / (1 - (1/p)) 
				   = log(p) / (p-1)

				   This contribution occurs once for each
				   square root used for sieving. There are two
				   roots for each factor base prime, unless
				   the prime divides k*n. In that case there 
				   is only one root */

				if (knmodp == 0)
					scores[j] -= contrib;
				else
					scores[j] -= 2 * contrib;
			}
		}

	}

	/* use the multiplier that generates the best score */
    double best_score_below_73 = 1000.0;
    uint8_t best_mult_below_73 = 1;
	best_score = 1000.0;
	best_mult = 1;
	for (i = 0; i < num_multipliers; i++) {
		double score = scores[i];
        if ((score < best_score_below_73) && (mult_list[i] <= 73)) {
            best_score_below_73 = score;
            best_mult_below_73 = mult_list[i];
        }
		if (score < best_score) {
			best_score = score;
			best_mult = mult_list[i];
		}
	}
    if (best_mult_below_73 != best_mult)
    {
        //gmp_printf("n = %Zd\n", n);
        //printf("found better multiplier: %u, compare to original list's best: %u\n",
        //    best_mult, best_mult_below_73);
    }
	return best_mult;
}

void rebuild_graph(static_conf_t *sconf, siqs_r *relation_list, int num_relations)
{
    int i;
    for (i = 0; i < num_relations; i++) {
        siqs_r *r = relation_list + i;
        if (r->large_prime[0] != r->large_prime[1]) {
            yafu_add_to_cycles(sconf, sconf->obj->flags, r->large_prime[0], r->large_prime[1]);
            sconf->num_cycles++;
        }
        else {
            sconf->num_relations++;
        }
    }
    return;
}

void rebuild_graph3(static_conf_t *sconf, siqs_r *relation_list, int num_relations)
{
	int i;
	uint32_t *hashtable = sconf->cycle_hashtable;
	qs_cycle_t *table = sconf->cycle_table;

	memset(hashtable, 0, sizeof(uint32_t) * (1 << QS_LOG2_CYCLE_HASH));
	sconf->num_cycles = 0;
	sconf->vertices = 0;
	sconf->components = 0;
	sconf->cycle_table_size = 1;
	memset(table, 0, sconf->cycle_table_alloc * sizeof(qs_cycle_t));

	for (i = 0; i < num_relations; i++) {
		siqs_r *r = relation_list + i;
		uint32_t *p = r->large_prime;

		yafu_add_to_cyclesN(sconf, 0, p);
	}
	return;
}
