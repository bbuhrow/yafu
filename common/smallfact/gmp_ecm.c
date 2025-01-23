/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: gmp_ecm.c 915 2013-07-13 12:18:35Z brgladman $
--------------------------------------------------------------------*/

/* Conditionally compiled interface to the GMP-ECM library. This
   is used in order to plug the P-1, P+1 and ECM algorithms into the
   recursive framework used to completely factorize input numbers.
   GMP-ECM is a very advanced implementation of these algorithms,
   far more advanced than I'm able to do without an effort comparable 
   to reimplementing all of msieve */

#include <ms_common.h>

#ifndef HAVE_GMP_ECM

uint32 ecm_pp1_pm1(msieve_obj *obj, mp_t *n, mp_t *reduced_n, 
		   factor_list_t *factor_list) {

	logprintf(obj, "no P-1/P+1/ECM available, skipping\n");
	return 0;
}

#else

#include <ecm.h>

/* description of the work to do in order to find factors
   with a given number of digits. I should stress that this
   code is not a state-of-the-art ECM driver that can find 
   40- and 50-digit factors; we use ECM and friends here
   primarily to filter the inputs to the main QS and NFS code */

typedef struct {
   uint32 digits;
   double stage_1_bound;
   uint32 num_ecm_trials;
} work_t;

static const work_t work_table[] = {
   {15,    2000,   30},
   {20,   11000,   74},
   {25,   50000,  214},
   {30,  250000,  430},
   {35, 1000000,  904},
   {40, 3000000, 2350},
};

/* to save work on P+-1 runs, we remember how much
   stage 1 work was done in each instance and reuse
   the stage 1 output at higher digit levels */

typedef struct {
	double stage_1_done;
	mpz_t start_val;
} pm1_pp1_t;

#define NUM_TABLE_ENTRIES (sizeof(work_table) / sizeof(work_t))

#define NUM_PM1_TRIALS 1
#define NUM_PP1_TRIALS 1
#define NUM_NON_ECM (NUM_PM1_TRIALS + NUM_PP1_TRIALS)

/* the number of times we can tolerate any algorithm
   finding a factor of n that equals n. This typically
   occurs when n has many small factors, and it's better
   to just run QS instead of wasting time computing
   trivial factors here */

#define MAX_TRIVIAL 3

/*--------------------------------------------------------------------*/
static uint32 postprocess(msieve_obj *obj, mpz_t factor, 
			  mpz_t n, pm1_pp1_t *non_ecm_vals,
			  factor_list_t *factor_list) {

	/* handling for factors found. When complete, n
	   is changed to new input to be factored */

	uint32 i;
	mpz_t q, r;
	mp_t mp_factor;
	uint32 is_prime0, is_prime1;

	/* divide out all instances of factor from n */

	mpz_init(q);
	mpz_init(r);
	while (1) {
		mpz_tdiv_qr(q, r, n, factor);
		if (mpz_cmp_ui(q, (unsigned long)0) == 0 || 
		    mpz_cmp_ui(r, (unsigned long)0) != 0)
			break;

		mpz_set(n, q);
	}
	mpz_clear(q);
	mpz_clear(r);

	if (mpz_cmp_ui(n, (unsigned long)1) == 0) {
		mpz_set(n, factor);
	}
	else {
		/* in the worst case, n and factor are both composite.
     		   We would have to give up factoring one of them and
     		   continue on the other */

		is_prime0 = mpz_probab_prime_p(factor, 1);
		is_prime1 = mpz_probab_prime_p(n, 1);
		if ( (is_prime1 > 0 && is_prime0 == 0) ||
		     (is_prime1 == 0 && is_prime0 == 0 &&
				mpz_cmp(factor, n) > 0) ) {
			mpz_swap(factor, n);
		}

		/* switch the P+-1 stage 1 values to 
		   be modulo the new n */

		for (i = 0; i < NUM_NON_ECM; i++) {
			mpz_mod(non_ecm_vals[i].start_val,
				non_ecm_vals[i].start_val, n);
		}
	}

	gmp2mp(factor, &mp_factor);
	return factor_list_add(obj, factor_list, &mp_factor);
}

/*--------------------------------------------------------------------*/
static uint32 choose_max_digits(msieve_obj *obj, uint32 bits) {

	/* choose the amount of work to do. We want the
	   chosen digit level to be a small fraction of what
	   QS and NFS would need */

	uint32 max_digits = 15;

	if (bits == 0)
		return 0;

	if (obj->flags & MSIEVE_FLAG_DEEP_ECM) {
		if (bits > 220) {
			if (bits < 280)
				max_digits = 20;
			else if (bits < 320)
				max_digits = 25;
			else if (bits < 360)
				max_digits = 30;
			else if (bits < 400)
				max_digits = 35;
			else
				max_digits = 40;
		}
	}
	return max_digits;
}

/*--------------------------------------------------------------------*/
/* whenever a factor is found, recalculate the appropriate
   digit level for P+-1 and ECM, since smaller inputs make
   the QS code run faster too. Stop trying to find factors
   when it becomes faster to switch to QS */
 
#define HANDLE_FACTOR_FOUND(alg_name)				\
	if (status > 0) {					\
		logprintf(obj, alg_name " stage %d factor found\n",	\
				status);			\
		if (mpz_cmp(gmp_factor, gmp_n) == 0) {		\
			if (++num_trivial == MAX_TRIVIAL)	\
				goto clean_up;			\
		}						\
		else {						\
			factor_found = 1;			\
			bits = postprocess(obj, gmp_factor, 	\
					gmp_n, non_ecm_vals, 	\
					factor_list);		\
			max_digits = choose_max_digits(obj, bits); \
			if (curr_work->digits > max_digits)	\
				goto clean_up;			\
		}						\
	}							\
	if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)		\
		goto clean_up;

uint32 ecm_pp1_pm1(msieve_obj *obj, mp_t *n, mp_t *reduced_n, 
		   factor_list_t *factor_list) {

	uint32 i, j, max_digits;
	uint32 bits;
	mpz_t gmp_n, gmp_factor;
	pm1_pp1_t non_ecm_vals[NUM_NON_ECM];
	uint32 factor_found = 0;
	uint32 num_trivial = 0;
	ecm_params params;

	/* factorization will continue until any remaining
	   composite cofactors are smaller than the cutoff for
	   using other methods */

	/* initialize */

	mpz_init(gmp_n);
	mpz_init(gmp_factor);
	for (i = 0; i < NUM_NON_ECM; i++) {
		pm1_pp1_t *tmp = non_ecm_vals + i;
		tmp->stage_1_done = 1.001;
		mpz_init_set_ui(tmp->start_val, (unsigned long)0);
	}
	ecm_init(params);
	gmp_randseed_ui(params->rng, get_rand(&obj->seed1, &obj->seed2));
	mp2gmp(n, gmp_n);
	max_digits = choose_max_digits(obj, mp_bits(n));

	/* for each digit level */

	for (i = 0; i < NUM_TABLE_ENTRIES; i++) {
		const work_t *curr_work = work_table + i;
		uint32 log_rate, next_log;
		int status;

		if (curr_work->digits > max_digits)
			break;

		logprintf(obj, "searching for %u-digit factors\n",
				curr_work->digits);

		/* perform P-1. Stage 1 runs much faster for this
		   algorithm than it does for ECM, so crank up the
		   stage 1 bound
		
		   We save the stage 1 output to reduce the work at
		   the next digit level */

		for (j = 0; j < NUM_PM1_TRIALS; j++) {
			pm1_pp1_t *tmp = non_ecm_vals + j;

			params->method = ECM_PM1;
			params->B1done = tmp->stage_1_done;
			mpz_set(params->x, tmp->start_val);

			status = ecm_factor(gmp_factor, gmp_n,
					10 * curr_work->stage_1_bound, params);

			tmp->stage_1_done = params->B1done;
			mpz_set(tmp->start_val, params->x);
			HANDLE_FACTOR_FOUND("P-1");
		}

		/* perform P+1. Stage 1 runs somewhat faster for this
		   algorithm than it does for ECM, so crank up the
		   stage 1 bound
		
		   We save the stage 1 output to reduce the work at
		   the next digit level */

		for (j = 0; j < NUM_PP1_TRIALS; j++) {
			pm1_pp1_t *tmp = non_ecm_vals + NUM_PM1_TRIALS + j;

			params->method = ECM_PP1;
			params->B1done = tmp->stage_1_done;
			mpz_set(params->x, tmp->start_val);

			status = ecm_factor(gmp_factor, gmp_n,
					5 * curr_work->stage_1_bound, params);

			tmp->stage_1_done = params->B1done;
			mpz_set(tmp->start_val, params->x);
			HANDLE_FACTOR_FOUND("P+1");
		}

		/* run the specified number of ECM curves. Because
		   so many curves could be necessary, just start each
		   one from scratch */

		log_rate = 0;
		if (obj->flags & (MSIEVE_FLAG_USE_LOGFILE |
	    		   	  MSIEVE_FLAG_LOG_TO_STDOUT)) {
			if (curr_work->digits >= 35)
				log_rate = 1;
			else if (curr_work->digits >= 30)
				log_rate = 5;
			else if (curr_work->digits >= 25)
				log_rate = 20;
		}

		next_log = log_rate;
		for (j = 0; j < curr_work->num_ecm_trials; j++) {

			params->method = ECM_ECM;
			params->B1done = 1.0;
			mpz_set_ui(params->x, (unsigned long)0);
			mpz_set_ui(params->sigma, (unsigned long)0);

			status = ecm_factor(gmp_factor, gmp_n,
					curr_work->stage_1_bound, params);
			HANDLE_FACTOR_FOUND("ECM");

			if (log_rate && j == next_log) {
				next_log += log_rate;
				fprintf(stderr, "%u of %u curves\r",
						j, curr_work->num_ecm_trials);
				fflush(stderr);
			}
		}

		if (log_rate)
			fprintf(stderr, "\ncompleted %u ECM curves\n", j);
	}

clean_up:
	ecm_clear(params);
	gmp2mp(gmp_n, reduced_n);
	mpz_clear(gmp_n);
	mpz_clear(gmp_factor);
	for (i = 0; i < NUM_NON_ECM; i++)
		mpz_clear(non_ecm_vals[i].start_val);
	return factor_found;
}

#endif /* HAVE_GMP_ECM */
