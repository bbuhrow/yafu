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
#include "factor.h"
#include "common.h"
#include "util.h"
#include <gmp.h>

#define DEFINED 1
#define NUM_SQUFOF_MULT 38

#if defined (TARGET_KNC) || defined(USE_AVX512F)
#define NUM_LANES 16
#else
#define NUM_LANES 8
#endif

//#define PRINT_DEBUG

/*
implements shanks's squfof algorihm.  priceless help from the papers
of Stephen McMath, Danial Shanks, and Jason Gower
*/

typedef struct
{
    uint64 *N;
    uint64 *mN;
    uint32 *listref;
    uint32 *mult;
    uint32 *valid;
    uint32 *P;
    uint32 *bn;
    uint32 *Qn;
    uint32 *Q0;
    uint32 *b0;
    uint32 *it;
    uint32 *imax;
    uint32 *f;
    int *maxrounds;
    int *rounds;
    int *multnum;
    int *active;
} par_mult_t;

typedef struct
{
    uint32 mult;
    uint32 valid;
    uint32 P;
    uint32 bn;
    uint32 Qn;
    uint32 Q0;
    uint32 b0;
    uint32 it;
    uint32 imax;	
    uint32 maxrounds;
    uint32 rounds;
} mult_t;


// local functions
void par_shanks_mult_unit(par_mult_t *mult_save);
void par_shanks_mult_unit_asm(par_mult_t *mult_save);
void par_shanks_mult_unit_asm2(par_mult_t *mult_save);
void shanks_mult_unit(uint64 N, mult_t *mult_save, uint64 *f);
int init_multipliers(mult_t **savedata, par_mult_t batch_data, uint64 N, int lane, 
    int num_in, mpz_t gmptmp);
int init_next_multiplier(par_mult_t mult_save, int lane,
    int num_in, mpz_t gmptmp);
void copy_mult_save(par_mult_t batch_data, int dest_lane, int src_lane);
void save_multiplier_data(par_mult_t batch_data, mult_t **savedata, int lane);
void load_multiplier_data(par_mult_t batch_data, mult_t **savedata, int lane, int multnum);
int get_next_multiplier(par_mult_t batch_data, mult_t **savedata, int lane);

// larger list of square-free multipliers from Dana Jacobsen.  Together with fewer
// iterations per round and racing, this works faster on average.
const int multipliers[NUM_SQUFOF_MULT] = {
    3 * 5 * 7 * 11, 3 * 5 * 7, 3 * 5 * 7 * 11 * 13, 3 * 5 * 7 * 13, 3 * 5 * 7 * 11 * 17, 3 * 5 * 11,
    3 * 5 * 7 * 17, 3 * 5, 3 * 5 * 7 * 11 * 19, 3 * 5 * 11 * 13, 3 * 5 * 7 * 19, 3 * 5 * 7 * 13 * 17,
    3 * 5 * 13, 3 * 7 * 11, 3 * 7, 5 * 7 * 11, 3 * 7 * 13, 5 * 7,
    3 * 5 * 17, 5 * 7 * 13, 3 * 5 * 19, 3 * 11, 3 * 7 * 17, 3,
    3 * 11 * 13, 5 * 11, 3 * 7 * 19, 3 * 13, 5, 5 * 11 * 13,
    5 * 7 * 19, 5 * 13, 7 * 11, 7, 3 * 17, 7 * 13,
    11, 1 };


uint64 sp_shanks_loop(mpz_t N, fact_obj_t *fobj)
{
	// call shanks with multiple small multipliers
	int i, rounds,j;
	uint64 n64, nn64, f64, big1, big2;
	mult_t mult_save[NUM_SQUFOF_MULT];
	mpz_t gmptmp;
	
	big1 = 0xFFFFFFFFFFFFFFFFULL;
	big2 = 0x3FFFFFFFFFFFFFFFULL;

	if (mpz_sizeinbase(N,2) > 62)
	{
        if (VFLAG > 0)
		    printf("N too big (%d bits), exiting...\n", (int)mpz_sizeinbase(N,2));
		return 1;
	}	

	n64 = mpz_get_64(N);

	if (mpz_sizeinbase(N,2) <= 40)
		return LehmanFactor(n64, 3.5, 0, 0.1);

	//default return value
	f64 = 1;

	if (n64 <= 3)
		return n64;

	mpz_init(gmptmp);

	for (i=NUM_SQUFOF_MULT-1;i>=0;i--)
	{
		// can we multiply without overflowing 64 bits?
		if (big2/(uint64)multipliers[i] < n64)
		{
			//this multiplier makes the input bigger than 64 bits
			mult_save[i].mult = multipliers[i];
			mult_save[i].valid = 0;
			continue;
		}

		//form the multiplied input
		nn64 = n64 * (uint64)multipliers[i];

		mult_save[i].mult = multipliers[i];
		mult_save[i].valid = 1;

		//set imax = N^1/4
		mpz_set_64(gmptmp, nn64);
		mpz_sqrt(gmptmp, gmptmp);	
		mult_save[i].b0 = mpz_get_ui(gmptmp);
		mult_save[i].imax = (uint32)sqrt((double)mult_save[i].b0) / 16;

		//set up recurrence
		mult_save[i].Q0 = 1;
		mult_save[i].P = mult_save[i].b0;
		mult_save[i].Qn = (uint32)(nn64 - 
			(uint64)mult_save[i].b0 * (uint64)mult_save[i].b0);
			
		if (mult_save[i].Qn == 0)
		{
			//N is a perfect square
			f64 = (uint64)mult_save[i].b0;
			goto done;
		}
		mult_save[i].bn = (mult_save[i].b0 + mult_save[i].P)
			/ mult_save[i].Qn;
		mult_save[i].it = 0;

	}

	//now process the multipliers a little at a time.  this allows more
	//multipliers to be tried in order to hopefully find one that 
	//factors the input quickly
    if (mpz_sizeinbase(N, 2) < 50)
        rounds = 4;
    else if (mpz_sizeinbase(N, 2) < 55)
        rounds = 8;
    else if (mpz_sizeinbase(N, 2) < 58)
        rounds = 16;
    else if (mpz_sizeinbase(N, 2) < 61)
        rounds = 24;
    else
        rounds = 32;

	for (i = 0; i < rounds; i++)
	{
		for (j=0; j < NUM_SQUFOF_MULT; j++)
        //for (j = NUM_SQUFOF_MULT - 1; j >= 0; j--)
		{
			if (mult_save[j].valid == 0)
				continue;

			//form the input
			nn64 = n64 * multipliers[j];
			//try to factor
			shanks_mult_unit(nn64,&mult_save[j],&f64);

			//check the output for a non-trivial factor
			if (f64 == (uint64)-1)
			{
				//this is an error condition, stop processing this multiplier
				mult_save[j].valid = 0;
			}
			else if (f64 > 1)
			{
				if (f64 != multipliers[j])
				{
					//factor found.  check for and remove small multiplier if necessary.
					nn64 = gcd64(f64,multipliers[j]);
					f64 /= nn64;

					if (f64 != 1)
					{
						//found a non-trivial factor, return it;
						goto done;
					}
					else
					{
						//found trivial factor, stop processing this multiplier
						mult_save[j].valid = 0;
					}
				}
				else
				{
					//found trivial factor, stop processing this multiplier
					mult_save[j].valid = 0;
				}
			}
		}
	}
	//default return value
	f64 = 1;

	//if we've got to here, then the number is still unfactored.  returning
	//a value of 1 signifies this
done:

	mpz_clear(gmptmp);

	return f64;
}


int par_shanks_loop(uint64 *N, uint64 *f, int num_in)
{
    // this routine takes a list of input 64-bit integers and factors them in 
    // parallel using AVX2 enhanced racing-SQUFOF.  If a factor is found, it is 
    // placed into the corresponding location of the array 'f'.  Else the location 
    // in 'f' is set to 1. 
    int i, rounds, j, all_done, list_position, num_successes, num_processed, num_active;
    uint64 n64, nn64, f64;
    par_mult_t mult_batch;
    mpz_t gmptmp;

    mult_t **save_data;

    mpz_init(gmptmp);

    save_data = (mult_t **)xmalloc(NUM_LANES * sizeof(mult_t *));
    for (i = 0; i < NUM_LANES; i++)
    {
        save_data[i] = (mult_t *)xmalloc(NUM_SQUFOF_MULT * sizeof(mult_t));
    }
    mult_batch.active = (int *)xmalloc_align(NUM_LANES * sizeof(int));
    mult_batch.b0 = (uint32 *)xmalloc_align(NUM_LANES * sizeof(uint32));
    mult_batch.bn = (uint32 *)xmalloc_align(NUM_LANES * sizeof(uint32));
    mult_batch.f = (uint32 *)xmalloc_align(NUM_LANES * sizeof(uint32));
    mult_batch.imax = (uint32 *)xmalloc_align(NUM_LANES * sizeof(uint32));
    mult_batch.it = (uint32 *)xmalloc_align(NUM_LANES * sizeof(uint32));
    mult_batch.listref = (uint32 *)xmalloc_align(NUM_LANES * sizeof(uint32));
    mult_batch.mN = (uint64 *)xmalloc_align(NUM_LANES * sizeof(uint64));
    mult_batch.mult = (uint32 *)xmalloc_align(NUM_LANES * sizeof(uint32));
    mult_batch.multnum = (int *)xmalloc_align(NUM_LANES * sizeof(int));
    mult_batch.N = (uint64 *)xmalloc_align(NUM_LANES * sizeof(uint64));
    mult_batch.P = (uint32 *)xmalloc_align(NUM_LANES * sizeof(uint32));
    mult_batch.Q0 = (uint32 *)xmalloc_align(NUM_LANES * sizeof(uint32));
    mult_batch.Qn = (uint32 *)xmalloc_align(NUM_LANES * sizeof(uint32));
    mult_batch.valid = (uint32 *)xmalloc_align(NUM_LANES * sizeof(uint32));
    mult_batch.maxrounds = (int *)xmalloc_align(NUM_LANES * sizeof(int));
    mult_batch.rounds = (int *)xmalloc_align(NUM_LANES * sizeof(int));


    num_active = 0;
    num_successes = 0;
    num_processed = 0;
    all_done = 0;
    list_position = 0;
    for (j = 0; j < NUM_LANES; j++)
    {
        mult_batch.active[j] = 0;
    }

    while ((num_processed < num_in) || (num_active > 0))
    {

        // fill in any empty lane in the batch.  A lane can be empty if it
        // has been factored, or has had all multipliers exhaused with no factor
        for (j = 0; j < NUM_LANES; j++)
        {
            if (list_position == num_in)
            {
                // stop filling lanes if we've reached the end of the list.
                // parallel squfof will continue to run on the incomplete
                // list until all are finished.
                break;
            }

            if (mult_batch.active[j] == 0)
            {
                int result;

                //printf("assigning input %d = %lu to lane %d\n", list_position, N[list_position], j);
                // replace this lane.
                mult_batch.listref[j] = list_position;
                mult_batch.N[j] = N[list_position++];

                if (mult_batch.N[j] < 3)
                {
                    //printf("rejecting input N = %lu\n", mult_batch.N[j]);
                    f[mult_batch.listref[j]] = mult_batch.N[j];
                    j--;
                    continue;
                }

                // initialize all multipliers
                result = init_multipliers(save_data, mult_batch, 
                    mult_batch.N[j], j, num_in, gmptmp);

                if (result == 1)
                {
                    // a result code of 1 from init_multipliers means that we found
                    // a factor while trying out multipliers.  Can happen if the input
                    // is already a perfect square.
                    // record this (unlikely) success, and try to fill the lane again.
                    num_processed++;
                    num_successes++;
                    f[mult_batch.listref[j]] = (uint64)sqrt(N[j]);
                    j--;
                }
            }
        }

        if (num_active < NUM_LANES)
        {
            int active_lane = -1;

            // if we are at the end of the list, some lanes may be inactive.  Fill the
            // save_data for these lanes with copies of an active lane, to prevent
            // problems with divide-by-0.
            for (j = 0; j < NUM_LANES; j++)
            {
                if (mult_batch.active[j] == 1)
                {
                    active_lane = j;
                    break;
                }
            }

            if (active_lane < 0)
            {
                //printf("no active lanes found\n");
                active_lane = 0;
                break;
            }

            for (j = 0; j < NUM_LANES; j++)
            {
                if (mult_batch.active[j] == 0)
                {
                    copy_mult_save(mult_batch, j, active_lane);

                    // keep it inactive
                    mult_batch.active[j] = 0;
                }
            }

        }

        //for (j = 0; j < NUM_LANES; j++)
        //{
        //    printf("lane %d has input %lu with multiplier %d starting on iteration %u of %u\n", 
        //        j, mult_batch.mN[j], mult_batch.mult[j], mult_batch.it[j], mult_batch.imax[j]);
        //}


        // run parallel squfof
#if (USE_AVX2)
#if defined(__INTEL_COMPILER)
        par_shanks_mult_unit(&mult_batch);
#else
        par_shanks_mult_unit_asm(&mult_batch);
#endif
#else
        par_shanks_mult_unit(&mult_batch);
#endif

        // examine the batch and:
        // 1) flag completed numbers
        // 2a) change multpliers on any number that has more than imax iterations run.
        // 2b) or has found a trivial factor 
        // 3) check if we are all done
        for (j = 0; j < NUM_LANES; j++)
        {
            int result;

            if (mult_batch.active[j] == 1)
            {

                f64 = mult_batch.f[j];

                //check the output for a non-trivial factor
                if (f64 == (uint64)-1)
                {
                    //this is an error condition, stop processing this multiplier
                    save_data[j][mult_batch.multnum[j]].valid = 0;
                    get_next_multiplier(mult_batch, save_data, j);

#ifdef PRINT_DEBUG
                    printf("failed on input %lu \n", N[mult_batch.listref[j]]);
#endif
                }
                else if (f64 > 1)
                {
                    //found a non-trivial factor, return it;
                    f[mult_batch.listref[j]] = f64;
                    num_successes++;
                    num_processed++;

#ifdef PRINT_DEBUG
                    printf("found factor %u of input %lu on multiplier %u, round %u, iteration %u\n",
                        f64, N[mult_batch.listref[j]], mult_batch.multnum[j], mult_batch.rounds[j], mult_batch.it[j]);
#endif

                    // and flag this lane to be replaced by a new input
                    mult_batch.active[j] = 0;
                }
                else
                {

                    // if it's not an error or a factor, check to see if we should 
                    // continue or not.
                    if (mult_batch.it[j] >= mult_batch.imax[j])
                    {
                        int next;
                        int old_multnum;

                        // we've done enough iterations on this multiplier, switch
                        // to the next one.  First reset the number of performed iterations.
                        // The exact number of iterations doesn't really matter,
                        // but we have to maintain the parity of the iteration count.
                        while (mult_batch.it[j] >= mult_batch.imax[j])
                        {
                            if (mult_batch.imax[j] & 1)
                            {
                                mult_batch.it[j] -= (mult_batch.imax[j] - 1);
                            }
                            else
                            {
                                mult_batch.it[j] -= mult_batch.imax[j];
                            }
                        }
                        
                        // try to get the next multiplier
                        old_multnum = mult_batch.multnum[j];
                        save_multiplier_data(mult_batch, save_data, j);
                        next = get_next_multiplier(mult_batch, save_data, j);

                        if (next < 0)
                        {
                            // no more multipliers for this input, try another round
                            // on the same multiplier.
                            load_multiplier_data(mult_batch, save_data, j, old_multnum);
                            mult_batch.rounds[j]++;
                            save_data[j][old_multnum].rounds = mult_batch.rounds[j];

                            if (mult_batch.rounds[j] >= mult_batch.maxrounds[j])
                            {
                                // also done with rounds.  give up.
                                mult_batch.active[j] = 0;
                                num_processed++;

#ifdef PRINT_DEBUG
                                printf("giving up on input %lu \n", N[mult_batch.listref[j]]);
#endif
                            }
                            else
                            {
                                // start over with the first multiplier, to
                                // do another round of iterations on each.
                                mult_batch.multnum[j] = -1;

                                // try to get the next multiplier starting from the beginning                                
                                next = get_next_multiplier(mult_batch, save_data, j);
                                if (next < 0)
                                {
                                    // all multipliers invalid.  give up.
                                    mult_batch.active[j] = 0;
                                    num_processed++;

#ifdef PRINT_DEBUG
                                    printf("giving up on input %lu \n", N[mult_batch.listref[j]]);
#endif
                                }
                            }
                        }
                    }
                }
            }
        }

        // count the number of active lanes, to see if we should keep going.
        num_active = 0;
        for (j = 0; j < NUM_LANES; j++)
        {
            num_active += mult_batch.active[j];
        }

        //printf("now %d successes out of %d processed, %d currently active lanes\n", 
        //    num_successes, num_processed, num_active);

    }

    mpz_clear(gmptmp);

    align_free(mult_batch.active);
    align_free(mult_batch.b0);
    align_free(mult_batch.bn);
    align_free(mult_batch.f);
    align_free(mult_batch.imax);
    align_free(mult_batch.it);
    align_free(mult_batch.listref);
    align_free(mult_batch.mN);
    align_free(mult_batch.mult);
    align_free(mult_batch.multnum);
    align_free(mult_batch.N);
    align_free(mult_batch.P);
    align_free(mult_batch.Q0);
    align_free(mult_batch.Qn);
    align_free(mult_batch.valid);
    align_free(mult_batch.maxrounds);
    align_free(mult_batch.rounds);

    for (i = 0; i < NUM_LANES; i++)
    {
        free(save_data[i]);
    }
    free(save_data);

    return num_successes;
}

void copy_mult_save(par_mult_t batch_data, int dest_lane, int src_lane)
{
    // copy fields from src lane to dest lane.
    batch_data.active[dest_lane]    = batch_data.active[src_lane];
    batch_data.b0[dest_lane]        = batch_data.b0[src_lane];
    batch_data.bn[dest_lane]        = batch_data.bn[src_lane];
    batch_data.f[dest_lane]         = batch_data.f[src_lane];
    batch_data.imax[dest_lane]      = batch_data.imax[src_lane];
    batch_data.it[dest_lane]        = batch_data.it[src_lane];
    batch_data.listref[dest_lane]   = batch_data.listref[src_lane];
    batch_data.maxrounds[dest_lane] = batch_data.maxrounds[src_lane];
    batch_data.mN[dest_lane]        = batch_data.mN[src_lane];
    batch_data.mult[dest_lane]      = batch_data.mult[src_lane];
    batch_data.multnum[dest_lane]   = batch_data.multnum[src_lane];
    batch_data.N[dest_lane]         = batch_data.N[src_lane];
    batch_data.P[dest_lane]         = batch_data.P[src_lane];
    batch_data.Q0[dest_lane]        = batch_data.Q0[src_lane];
    batch_data.Qn[dest_lane]        = batch_data.Qn[src_lane];
    batch_data.rounds[dest_lane]    = batch_data.rounds[src_lane];
    batch_data.valid[dest_lane]     = batch_data.valid[src_lane];

    return;
}

void save_multiplier_data(par_mult_t batch_data, mult_t **savedata, int lane)
{
    // save current data from this lane
    int i = batch_data.multnum[lane];
    savedata[lane][i].b0 = batch_data.b0[lane];
    savedata[lane][i].bn = batch_data.bn[lane];
    savedata[lane][i].it = batch_data.it[lane];
    savedata[lane][i].P = batch_data.P[lane];
    savedata[lane][i].Q0 = batch_data.Q0[lane];
    savedata[lane][i].Qn = batch_data.Qn[lane];
    savedata[lane][i].valid = batch_data.valid[lane];
    savedata[lane][i].rounds = batch_data.rounds[lane];
    return;
}

void load_multiplier_data(par_mult_t batch_data, mult_t **savedata, int lane, int multnum)
{
    // save current data from this lane
    batch_data.b0[lane] = savedata[lane][multnum].b0;
    batch_data.bn[lane] = savedata[lane][multnum].bn;
    batch_data.it[lane] = savedata[lane][multnum].it;
    batch_data.P[lane] = savedata[lane][multnum].P;
    batch_data.Q0[lane] = savedata[lane][multnum].Q0;
    batch_data.Qn[lane] = savedata[lane][multnum].Qn;
    batch_data.valid[lane] = savedata[lane][multnum].valid;
    batch_data.rounds[lane] = savedata[lane][multnum].rounds;
    return;
}

int get_next_multiplier(par_mult_t batch_data, mult_t **savedata, int lane)
{
    
    int i;
    int found = -1;

    // restore infomation about the next valid multiplier for this lane
    batch_data.multnum[lane]++;
    for (i = batch_data.multnum[lane]; i < NUM_SQUFOF_MULT; i++)
    {
        if (savedata[lane][i].valid == 0)
        {
            continue;
        }

        batch_data.b0[lane] = savedata[lane][i].b0;
        batch_data.bn[lane] = savedata[lane][i].bn;
        batch_data.it[lane] = savedata[lane][i].it;
        batch_data.mN[lane] = batch_data.N[lane] * multipliers[i];
        batch_data.imax[lane] = savedata[lane][i].imax;
        batch_data.mult[lane] = multipliers[i];
        batch_data.multnum[lane] = i;
        batch_data.P[lane] = savedata[lane][i].P;
        batch_data.Q0[lane] = savedata[lane][i].Q0;
        batch_data.Qn[lane] = savedata[lane][i].Qn;
        batch_data.valid[lane] = 1;
        batch_data.maxrounds[lane] = savedata[lane][i].maxrounds;
        batch_data.rounds[lane] = savedata[lane][i].rounds;
        found = 1;
        break;
    }

    return found;
}

int init_multipliers(mult_t **savedata, par_mult_t batch_data, 
    uint64 N, int lane, int num_in, mpz_t gmptmp)
{
    int i;
    int rounds;
    int success = 0;
    uint64 nn64;
    uint64 big2 = 0x3FFFFFFFFFFFFFFFULL;


    if (batch_data.active[lane] == 0)
    {
        for (i = 0; i < NUM_SQUFOF_MULT; i++)
        {

            // can we multiply without overflowing 64 bits?
            if ((big2 / (uint64)multipliers[i]) < N)
            {
                //this multiplier makes the input bigger than 64 bits
                savedata[lane][i].mult = multipliers[i];
                savedata[lane][i].valid = 0;
                continue;
            }

            //form the multiplied input
            nn64 = N * (uint64)multipliers[i];

            savedata[lane][i].mult = multipliers[i];
            savedata[lane][i].valid = 1;

            //set imax = N^1/4
            mpz_set_64(gmptmp, nn64);

            if (mpz_sizeinbase(gmptmp, 2) < 50)
                rounds = 4;
            else if (mpz_sizeinbase(gmptmp, 2) < 55)
                rounds = 8;
            else if (mpz_sizeinbase(gmptmp, 2) < 58)
                rounds = 16;
            else if (mpz_sizeinbase(gmptmp, 2) < 61)
                rounds = 24;
            else
                rounds = 32;            

            savedata[lane][i].maxrounds = rounds;
            savedata[lane][i].rounds = 0;

            mpz_sqrt(gmptmp, gmptmp);
            savedata[lane][i].b0 = mpz_get_ui(gmptmp);
            savedata[lane][i].imax = (uint32)sqrt((double)savedata[lane][i].b0) / 16;

            //set up recurrence
            savedata[lane][i].Q0 = 1;
            savedata[lane][i].P = savedata[lane][i].b0;
            savedata[lane][i].Qn = (uint32)(nn64 -
                (uint64)savedata[lane][i].b0 * (uint64)savedata[lane][i].b0);

            if (savedata[lane][i].Qn == 0)
            {
                // N is a perfect square - this number is factored.
                //printf("perfect sqrt in init_multipliers with N = %lu, mN = %lu\n",
                //    N, nn64);
                success = 1;
                savedata[lane][i].valid = 0;
                continue;
            }

            savedata[lane][i].bn = (savedata[lane][i].b0 + savedata[lane][i].P) / savedata[lane][i].Qn;
            savedata[lane][i].it = 0;

            // copy the first valid multiplier to the batch data structure
            batch_data.active[lane] = 1;
            batch_data.b0[lane] = savedata[lane][i].b0;
            batch_data.bn[lane] = savedata[lane][i].bn;
            batch_data.it[lane] = savedata[lane][i].it;
            batch_data.mN[lane] = nn64;
            batch_data.imax[lane] = savedata[lane][i].imax;
            batch_data.f[lane] = 0;
            batch_data.mult[lane] = multipliers[i];
            batch_data.multnum[lane] = i;
            batch_data.P[lane] = savedata[lane][i].P;
            batch_data.N[lane] = N;
            batch_data.Q0[lane] = savedata[lane][i].Q0;
            batch_data.Qn[lane] = savedata[lane][i].Qn;
            batch_data.valid[lane] = savedata[lane][i].valid;
            batch_data.maxrounds[lane] = rounds;
            batch_data.rounds[lane] = 0;
        }
    }

    // whether we found a factor or not (overwhelming not... but it could happen).
    return success;
}


int init_next_multiplier(par_mult_t mult_save, int lane, 
    int num_in, mpz_t gmptmp)
{
    int i;
    int rounds;
    int success = 0;
    uint64 nn64;
    uint64 n64 = mult_save.N[lane];
    uint64 big2 = 0x3FFFFFFFFFFFFFFFULL;

    // initially deactivate this lane
    mult_save.active[lane] = 0;

    // then see if we can find another multiplier to use
    for (i = mult_save.multnum[lane]; i < NUM_SQUFOF_MULT; i++)
    {
        mult_save.multnum[lane] = i;

        // can we multiply without overflowing 64 bits?
        if ((big2 / (uint64)multipliers[i]) < n64)
        {
            //this multiplier makes the input bigger than 64 bits
            mult_save.mult[lane] = multipliers[i];
            mult_save.valid[lane] = 0;
            continue;
        }

        //form the multiplied input
        nn64 = n64 * (uint64)multipliers[i];
        mult_save.mN[lane] = nn64;

        mult_save.mult[lane] = multipliers[i];
        mult_save.valid[lane] = 1;

        //set imax = N^1/4
        mpz_set_64(gmptmp, nn64);

        if (mpz_sizeinbase(gmptmp, 2) < 50)
            rounds = 4;
        else if (mpz_sizeinbase(gmptmp, 2) < 55)
            rounds = 8;
        else if (mpz_sizeinbase(gmptmp, 2) < 58)
            rounds = 16;
        else if (mpz_sizeinbase(gmptmp, 2) < 61)
            rounds = 24;
        else
            rounds = 32;

        mpz_sqrt(gmptmp, gmptmp);
        mult_save.b0[lane] = mpz_get_ui(gmptmp);
        mult_save.imax[lane] = (uint32)sqrt((double)mult_save.b0[lane]) * rounds / 16;

        //set up recurrence
        mult_save.Q0[lane] = 1;
        mult_save.P[lane] = mult_save.b0[lane];
        mult_save.Qn[lane] = (uint32)(nn64 -
            (uint64)mult_save.b0[lane] * (uint64)mult_save.b0[lane]);

        if (mult_save.Qn[lane] == 0)
        {
            // N is a perfect square - this number is factored.
            mult_save.f[mult_save.listref[lane]] = (uint64)mult_save.b0[lane];
            success = 1;
            break;
        }
        mult_save.bn[lane] = (mult_save.b0[lane] + mult_save.P[lane]) / mult_save.Qn[lane];
        mult_save.it[lane] = 0;

        // found a valid multiplier, go on to the next lane.
        mult_save.active[lane] = 1;
        break;
    }

    // whether we found a factor or not (overwhelming not... but it could happen).
    return success;
}



void shanks_mult_unit(uint64 N, mult_t *mult_save, uint64 *f)
{
	//use shanks SQUFOF to factor N.  almost all computation can be done with longs
	//input N < 2^63
	//return 1 in f if no factor is found
	uint32 imax,i,Q0,b0,Qn,bn,P,bbn,Ro,S,So,t1,t2;
	int j=0;	

	//initialize output
	*f=0;

	//load previous save point
	P = mult_save->P;
	bn = mult_save->bn;
	Qn = mult_save->Qn;
	Q0 = mult_save->Q0;
	b0 = mult_save->b0;
	i = mult_save->it;
	imax = i + mult_save->imax;

    //printf("scheduled to run %d iterations on input %lu\n", imax, N);
	
	while (1)
	{
		j=0;

		//i must be even on entering the unrolled loop below
		if (i & 0x1)
		{
			t1 = P;		//hold Pn for this iteration
			P = bn*Qn - P;
			t2 = Qn;	//hold Qn for this iteration
			Qn = Q0 + bn*(t1-P);
			Q0 = t2;	//remember last Q
			bn = (b0 + P) / Qn;
			i++;
		}

		while (1)
		{
			//at the start of every iteration, we need to know:
			//	P from the previous iteration
			//	bn from the previous iteration
			//	Qn from the previous iteration
			//	Q0 from the previous iteration
			//	iteration count, i
			if (i >= imax)
			{
				//haven't progressed to the next stage yet.  let another
				//multiplier try for awhile.  save state so we
				//know where to start back up in the next round
				mult_save->P = P;
				mult_save->bn = bn;
				mult_save->Qn = Qn;
				mult_save->Q0 = Q0;
				mult_save->it = i;
				//signal to do nothing but continue looking for a factor
				*f = 0;	
				return;
			}

			t1 = P;		//hold Pn for this iteration
			P = bn*Qn - P;
			t2 = Qn;	//hold Qn for this iteration
			Qn = Q0 + bn*(t1-P);
			Q0 = t2;	//remember last Q

			// this happens fairly often, but special-casing the division
			// still makes things slower (at least on modern intel chips)
			//if (Qn > ((b0 + P) >> 1))
			//	bn = 1;
			//else
			bn = (b0 + P) / Qn;		
			i++;

			//even iteration
			//check for square Qn = S*S
			// try sse2: broadcast Qn, AND, parallel ==, bytebitmask
			t2 = Qn & 31;
			if (t2 == 0 || t2 == 1 || t2 == 4 ||
				t2 == 9 || t2 == 16 || t2 == 17 || t2 == 25)			
			{
				// extra squaritude tests are also slower on modern intel chips
				//t2 = Qn & 63;
				//if (t2 < 32 || t2 == 33 || t2 == 36 || 
				//	t2 == 41 ||t2 == 49 || t2 == 57)
				//{
					t1 = (uint32)sqrt(Qn);
					if (Qn == t1 * t1)
						break;
				//}
			}

			//odd iteration
			t1 = P;		//hold Pn for this iteration
			P = bn*Qn - P;
			t2 = Qn;	//hold Qn for this iteration
			Qn = Q0 + bn*(t1-P);
			Q0 = t2;	//remember last Q

			bn = (b0 + P) / Qn;
			i++;
					
		}

        //printf("checking symmetry point after %d iterations on input %lu\n", i, N);

		//reduce to G0
		S = (int)sqrt(Qn);
		Ro = P + S*((b0 - P)/S);
		t1 = Ro;
		So = (uint32)(((int64)N - (int64)t1*(int64)t1)/(int64)S);
		bbn = (b0+Ro)/So;

		//search for symmetry point
		while (1)
		{
			t1 = Ro;		//hold Ro for this iteration
			Ro = bbn*So - Ro;
			t2 = So;		//hold So for this iteration
			So = S + bbn*(t1-Ro);
			S = t2;			//remember last S
			bbn = (b0+Ro)/So;
			
			//check for symmetry point
			if (Ro == t1)
				break;

			t1 = Ro;		//hold Ro for this iteration
			Ro = bbn*So - Ro;
			t2 = So;		//hold So for this iteration
			So = S + bbn*(t1-Ro);
			S = t2;			//remember last S
			bbn = (b0+Ro)/So;
			
			//check for symmetry point
			if (Ro == t1)
				break;

			t1 = Ro;		//hold Ro for this iteration
			Ro = bbn*So - Ro;
			t2 = So;		//hold So for this iteration
			So = S + bbn*(t1-Ro);
			S = t2;			//remember last S
			bbn = (b0+Ro)/So;
			
			//check for symmetry point
			if (Ro == t1)
				break;

			t1 = Ro;		//hold Ro for this iteration
			Ro = bbn*So - Ro;
			t2 = So;		//hold So for this iteration
			So = S + bbn*(t1-Ro);
			S = t2;			//remember last S
			bbn = (b0+Ro)/So;
			
			//check for symmetry point
			if (Ro == t1)
				break;

		}

		*f = gcd64(Ro,N);
		//found a factor - don't know yet if it's trivial or not.
		//we don't need to remember any state info, as one way or the other
		//this multiplier will be invalidated		
		if (*f > 1)
			return;
	}
}

void par_shanks_mult_unit(par_mult_t *mult_save)
{
    //use shanks SQUFOF on 8 inputs simultaneously using AVX2 for the bulk of the calculations.
    //input N < 2^63
    //return 1 in f if no factor is found
    uint32 imax;
    
#if defined(__GNUC__)
    uint32 *iterations = mult_save->it;
    uint32 *P = mult_save->P;
    uint32 *Qn = mult_save->Qn;
    uint32 *Q0 = mult_save->Q0;
    uint32 *bn = mult_save->bn;
    uint32 *b0 = mult_save->b0;
    __attribute__((aligned(64))) uint32 bbn[NUM_LANES];
    __attribute__((aligned(64))) uint32 Ro[NUM_LANES];
    __attribute__((aligned(64))) uint32 S[NUM_LANES];
    __attribute__((aligned(64))) uint32 So[NUM_LANES];
    __attribute__((aligned(64))) uint32 t1[NUM_LANES];
    __attribute__((aligned(64))) uint32 t2[NUM_LANES];
    __attribute__((aligned(64))) uint32 success_vec[NUM_LANES];
#else
    uint32 *iterations = mult_save->it;
    uint32 *P = mult_save->P;
    uint32 *Qn = mult_save->Qn;
    uint32 *Q0 = mult_save->Q0;
    uint32 *bn = mult_save->bn;
    uint32 *b0 = mult_save->b0;
    __declspec(align(64)) uint32 bbn[NUM_LANES];
    __declspec(align(64)) uint32 Ro[NUM_LANES];
    __declspec(align(64)) uint32 S[NUM_LANES];
    __declspec(align(64)) uint32 So[NUM_LANES];
    __declspec(align(64)) uint32 t1[NUM_LANES];
    __declspec(align(64)) uint32 t2[NUM_LANES];
    __declspec(align(64)) uint32 success_vec[NUM_LANES];
#endif


    int j = 0;
    int i = 0;
    int k;
    int success;

    // find max iterations such that all active lanes will have met this multiplier's imax
    imax = 0xffffffff;
    for (k = 0; k < NUM_LANES; k++)
    {
        if (mult_save->active[k] == 1)
        {
            if ((mult_save->imax[k] - mult_save->it[k]) < imax)
            {
                imax = mult_save->imax[k] - mult_save->it[k];
            }
        }
    }

    imax = MAX(imax, 128);

#pragma ivdep
#pragma vector aligned
    for (k = 0; k < NUM_LANES; k++)
    {
        mult_save->f[k] = 0;
    }

    while (1)
    {
        j = 0;

        //i must be even on entering the unrolled loop below
#pragma ivdep
#pragma vector aligned
        for (k = 0; k < NUM_LANES; k++)
        {
            if (iterations[k] & 0x1)
            {
                t1[k] = P[k];
                P[k] = bn[k] * Qn[k] - P[k];
                t2[k] = Qn[k];
                Qn[k] = Q0[k] + bn[k] * (t1[k] - P[k]);
                Q0[k] = t2[k];
                bn[k] = (b0[k] + P[k]) / Qn[k];
                iterations[k]++;
            }
        }

#ifdef PRINT_DEBUG
        printf("starting values:\n");

        printf("P = [");
        for (k = 0; k < NUM_LANES; k++)
            printf("%u ", P[k]);
        printf("]\n");

        printf("Qn = [");
        for (k = 0; k < NUM_LANES; k++)
            printf("%u ", Qn[k]);
        printf("]\n");

        printf("Q0 = [");
        for (k = 0; k < NUM_LANES; k++)
            printf("%u ", Q0[k]);
        printf("]\n");

        printf("bn = [");
        for (k = 0; k < NUM_LANES; k++)
            printf("%u ", bn[k]);
        printf("]\n");

        printf("b0 = [");
        for (k = 0; k < NUM_LANES; k++)
            printf("%u ", b0[k]);
        printf("]\n");

        printf("iterations = [");
        for (k = 0; k < NUM_LANES; k++)
            printf("%u ", iterations[k]);
        printf("]\n\n");
#endif

        while (1)
        {
            //at the start of every iteration, we need to know:
            //	P from the previous iteration
            //	bn from the previous iteration
            //	Qn from the previous iteration
            //	Q0 from the previous iteration
            //	iteration count, i
            if (i >= imax)
            {
                //haven't progressed to the next stage yet.  let another
                //multiplier try for awhile.  save state so we
                //know where to start back up in the next round

#pragma ivdep
#pragma vector aligned
                for (k = 0; k < NUM_LANES; k++)
                {
                    //signal to do nothing but continue looking for a factor
                    mult_save->f[k] = 0;
                }

#ifdef PRINT_DEBUG
                printf("ending values:\n");

                printf("P = [");
                for (k = 0; k < NUM_LANES; k++)
                    printf("%u ", P[k]);
                printf("]\n");

                printf("Qn = [");
                for (k = 0; k < NUM_LANES; k++)
                    printf("%u ", Qn[k]);
                printf("]\n");

                printf("Q0 = [");
                for (k = 0; k < NUM_LANES; k++)
                    printf("%u ", Q0[k]);
                printf("]\n");

                printf("bn = [");
                for (k = 0; k < NUM_LANES; k++)
                    printf("%u ", bn[k]);
                printf("]\n");

                printf("b0 = [");
                for (k = 0; k < NUM_LANES; k++)
                    printf("%u ", b0[k]);
                printf("]\n");

                printf("iterations = [");
                for (k = 0; k < NUM_LANES; k++)
                    printf("%u ", iterations[k]);
                printf("]\n\n");

                printf("i >= imax\n");
#endif
                
                return;
            }

#pragma ivdep
#pragma vector aligned
            for (k = 0; k < NUM_LANES; k++)
            {
                t1[k] = P[k];
                P[k] = bn[k] * Qn[k] - P[k];
                t2[k] = Qn[k];
                Qn[k] = Q0[k] + bn[k] * (t1[k] - P[k]);
                Q0[k] = t2[k];
                bn[k] = (b0[k] + P[k]) / Qn[k];
                iterations[k]++;
            }
            i++;

#ifdef PRINT_DEBUG
            printf("P = [");
            for (k = 0; k < NUM_LANES; k++)
                printf("%u ", P[k]);
            printf("]\n");

            printf("Qn = [");
            for (k = 0; k < NUM_LANES; k++)
                printf("%u ", Qn[k]);
            printf("]\n");

            printf("Q0 = [");
            for (k = 0; k < NUM_LANES; k++)
                printf("%u ", Q0[k]);
            printf("]\n");

            printf("bn = [");
            for (k = 0; k < NUM_LANES; k++)
                printf("%u ", bn[k]);
            printf("]\n");

            printf("b0 = [");
            for (k = 0; k < NUM_LANES; k++)
                printf("%u ", b0[k]);
            printf("]\n");

            printf("iterations = [");
            for (k = 0; k < NUM_LANES; k++)
                printf("%u ", iterations[k]);
            printf("]\n\n");
#endif

            //even iteration
            //check for square Qn = S*S
            success = 0;

#pragma ivdep
#pragma vector aligned
#pragma unroll
            for (k = 0; k < NUM_LANES; k++)
            {
                success_vec[k] = 0;

                // much faster this way, can be autovectorized.
                t1[k] = (uint32)sqrt(Qn[k]);
                if (Qn[k] == (t1[k] * t1[k]))
                {
                    success_vec[k] = 1;
                    success = 1;
                }
            }

            if (success)
            {
#ifdef PRINT_DEBUG
                printf("found a square\n");
                printf("ending values:\n");

                printf("P = [");
                for (k = 0; k < NUM_LANES; k++)
                    printf("%u ", P[k]);
                printf("]\n");

                printf("Qn = [");
                for (k = 0; k < NUM_LANES; k++)
                    printf("%u ", Qn[k]);
                printf("]\n");

                printf("Q0 = [");
                for (k = 0; k < NUM_LANES; k++)
                    printf("%u ", Q0[k]);
                printf("]\n");

                printf("bn = [");
                for (k = 0; k < NUM_LANES; k++)
                    printf("%u ", bn[k]);
                printf("]\n");

                printf("b0 = [");
                for (k = 0; k < NUM_LANES; k++)
                    printf("%u ", b0[k]);
                printf("]\n");

                printf("iterations = [");
                for (k = 0; k < NUM_LANES; k++)
                    printf("%u ", iterations[k]);
                printf("]\n");

                printf("success_vec = [");
                for (k = 0; k < NUM_LANES; k++)
                    printf("%u ", success_vec[k]);
                printf("]\n\n");

                
#endif

                // found at least one square.  check to see if it produces a factorization.
                break;
            }

            //odd iteration
#pragma ivdep
#pragma vector aligned
            for (k = 0; k < NUM_LANES; k++)
            {
                t1[k] = P[k];
                P[k] = bn[k] * Qn[k] - P[k];
                t2[k] = Qn[k];
                Qn[k] = Q0[k] + bn[k] * (t1[k] - P[k]);
                Q0[k] = t2[k];                
                bn[k] = (b0[k] + P[k]) / Qn[k];
                iterations[k]++;
            }
            i++;

#ifdef PRINT_DEBUG
            printf("P = [");
            for (k = 0; k < NUM_LANES; k++)
                printf("%u ", P[k]);
            printf("]\n");

            printf("Qn = [");
            for (k = 0; k < NUM_LANES; k++)
                printf("%u ", Qn[k]);
            printf("]\n");

            printf("Q0 = [");
            for (k = 0; k < NUM_LANES; k++)
                printf("%u ", Q0[k]);
            printf("]\n");

            printf("bn = [");
            for (k = 0; k < NUM_LANES; k++)
                printf("%u ", bn[k]);
            printf("]\n");

            printf("b0 = [");
            for (k = 0; k < NUM_LANES; k++)
                printf("%u ", b0[k]);
            printf("]\n");

            printf("iterations = [");
            for (k = 0; k < NUM_LANES; k++)
                printf("%u ", iterations[k]);
            printf("]\n\n");
#endif

        }

        // we're about to check for a factorization.  
        // save progress in case we find one and exit.
#pragma ivdep
#pragma vector aligned
        for (k = 0; k < NUM_LANES; k++)
        {
            //set default signal: to do nothing but continue looking for a factor
            mult_save->f[k] = 0;
        }

        success = 0;
        for (k = 0; k < NUM_LANES; k++)
        {
            if ((success_vec[k] == 1) && (mult_save->active[k] == 1))
            {
                int failsafe = 10000;

                //printf("symmetry check on input %lu in lane %d after %d iterations\n", 
                //    mult_save->mN[k], k, i);

                // reduce to G0
                S[k] = (int)sqrt(Qn[k]);
                Ro[k] = P[k] + S[k] * ((b0[k] - P[k]) / S[k]);
                t1[k] = Ro[k];
                So[k] = (uint32)(((int64)mult_save->mN[k] - (int64)t1[k] * (int64)t1[k]) / (int64)S[k]);
                bbn[k] = (b0[k] + Ro[k]) / So[k];


                // search for symmetry point
                while (failsafe)
                {

                    t1[k] = Ro[k];
                    Ro[k] = bbn[k] * So[k] - Ro[k];
                    t2[k] = So[k];
                    So[k] = S[k] + bbn[k] * (t1[k] - Ro[k]);
                    S[k] = t2[k];
                    bbn[k] = (b0[k] + Ro[k]) / So[k];

                    if (Ro[k] == t1[k])
                    {
                        // eliminate trivial cases before running the expensive gcd
                        if (Ro[k] > mult_save->mult[k])
                        {
                            // these run at almost exactly the same speed.  use the normal
                            // gcd both because it is more general and because it will
                            // pick up any future advances in hardware integer division.
                            mult_save->f[k] = gcd64(Ro[k], mult_save->N[k]);
                            //mult_save->f[k] = spBinGCD_odd(mult_save->N[k], Ro[k]);

                            if (mult_save->f[k] > 1)
                            {
                                success = 1;
                            }
                        }
                        break;
                    }

                    failsafe--;
                }
            }
        }


        if (success)
        {
            // at least one of the lanes was successful.  drop out of this
            // routine so we can replace that lane.  else we'll keep
            // going until we hit imax (or another square).
            return;
        }
    }
}

#if defined(USE_AVX2)

void par_shanks_mult_unit_asm(par_mult_t *mult_save)
{
    //use shanks SQUFOF on 8 inputs simultaneously using AVX2 for the bulk of the calculations.
    //input N < 2^63
    //return 1 in f if no factor is found
    uint32 imax;

#if defined(__INTEL_COMPILER)
    uint32 *iterations = mult_save->it;
    uint32 *P = mult_save->P;
    uint32 *Qn = mult_save->Qn;
    uint32 *Q0 = mult_save->Q0;
    uint32 *bn = mult_save->bn;
    uint32 *b0 = mult_save->b0;
    __declspec(aligned(64)) uint32 bbn[NUM_LANES];
    __declspec(aligned(64)) uint32 Ro[NUM_LANES];
    __declspec(aligned(64)) uint32 S[NUM_LANES];
    __declspec(aligned(64)) uint32 So[NUM_LANES];
    __declspec(aligned(64)) uint32 t1[NUM_LANES];
    __declspec(aligned(64)) uint32 t2[NUM_LANES];
    __declspec(aligned(64)) uint32 success_vec[NUM_LANES];
    uint64 conversion[4];
    double fudge[4];
#else
    uint32 *iterations = mult_save->it;
    uint32 *P = mult_save->P;
    uint32 *Qn = mult_save->Qn;
    uint32 *Q0 = mult_save->Q0;
    uint32 *bn = mult_save->bn;
    uint32 *b0 = mult_save->b0;
    uint32 bbn[NUM_LANES];
    uint32 Ro[NUM_LANES];
    uint32 S[NUM_LANES];
    uint32 So[NUM_LANES];
    uint32 t1[NUM_LANES];
    uint32 t2[NUM_LANES];
    uint32 success_vec[NUM_LANES];
    uint64 conversion[4];
    double fudge[4];
#endif


    int j = 0;
    int i = 0;
    int k;
    int success;

    conversion[0] = 0x4330000000000000;
    conversion[1] = 0x4330000000000000;
    conversion[2] = 0x4330000000000000;
    conversion[3] = 0x4330000000000000;

    fudge[0] = 1/65536.0;
    fudge[1] = 1/65536.0;
    fudge[2] = 1/65536.0;
    fudge[3] = 1/65536.0;

    // find max iterations such that all active lanes will have met this multiplier's imax
    imax = 0xffffffff;
    for (k = 0; k < NUM_LANES; k++)
    {
        if (mult_save->active[k] == 1)
        {
            if ((mult_save->imax[k] - mult_save->it[k]) < imax)
            {
                imax = mult_save->imax[k] - mult_save->it[k];
            }
        }
    }

    //imax = MAX(imax, 128);
    //printf("scheduled to run %d iterations\n", imax);

    // load previous save point
#pragma novector
    for (k = 0; k < NUM_LANES; k++)
    {
        mult_save->f[k] = 0;
        success_vec[k] = 0;
    }

    while (1)
    {
        int kk;

        j = 0;

        //i must be even on entering the unrolled loop below
#pragma novector
        for (k = 0; k < NUM_LANES; k++)
        {
            if (iterations[k] & 0x1)
            {
                t1[k] = P[k];
                P[k] = bn[k] * Qn[k] - P[k];
                t2[k] = Qn[k];
                Qn[k] = Q0[k] + bn[k] * (t1[k] - P[k]);
                Q0[k] = t2[k];
                bn[k] = (b0[k] + P[k]) / Qn[k];
                iterations[k]++;
            }
        }

#ifdef PRINT_DEBUG
        printf("starting values:\n");

        printf("P = [");
        for (k = 0; k < NUM_LANES; k++)
            printf("%u ", P[k]);
        printf("]\n");

        printf("Qn = [");
        for (k = 0; k < NUM_LANES; k++)
            printf("%u ", Qn[k]);
        printf("]\n");

        printf("Q0 = [");
        for (k = 0; k < NUM_LANES; k++)
            printf("%u ", Q0[k]);
        printf("]\n");

        printf("bn = [");
        for (k = 0; k < NUM_LANES; k++)
            printf("%u ", bn[k]);
        printf("]\n");

        printf("b0 = [");
        for (k = 0; k < NUM_LANES; k++)
            printf("%u ", b0[k]);
        printf("]\n");

        printf("iterations = [");
        for (k = 0; k < NUM_LANES; k++)
            printf("%u ", iterations[k]);
        printf("]\n\n");
#endif

        //printf("entering loop at iteration %d (of %d)\n", i, imax);
        //fflush(stdout);
        //success = imax / 8;
        //for (kk = 0; kk < success; kk++)
        //{

        kk = i;

            __asm__
                (
                /* load quantities in xmm registers */
                "vmovdqa    (%1), %%ymm0 \n\t" /* P */
                "vmovdqa    (%2), %%ymm1 \n\t" /* Qn */
                "vmovdqa    (%3), %%ymm2 \n\t" /* Q0 */
                "vmovdqa    (%4), %%ymm3 \n\t" /* bn */
                "vmovdqa    (%5), %%ymm4 \n\t" /* b0 */
                "movl	    $1, %%r11d \n\t"		        /* fill utility register with 1's */
                "vmovd      %%r11d, %%xmm8 \n\t"            /* broadcast 1's to ymm8 */
                "vpshufd	    $0, %%xmm8, %%xmm8 \n\t"    /* broadcast 1's to ymm8 */
                "vinserti128    $1, %%xmm8, %%ymm8, %%ymm8 \n\t"
                "vmovd      %%r11d, %%xmm10 \n\t"               /* broadcast 1's to ymm10 */
                "vpshufd	    $0, %%xmm10, %%xmm10 \n\t"      /* broadcast 1's to ymm10 */
                "vpaddd     %%xmm10, %%xmm10, %%xmm10 \n\t"     /* make it an array of 2's */
                "vcvtdq2ps  %%xmm10, %%xmm10 \n\t"	            /* convert to single precision */
                "vcvtps2pd  %%xmm10, %%ymm5 \n\t"	            /* convert to double precision (2's array) */
                "vrcpps	%%xmm10, %%xmm10 \n\t"		            /* approx 1/2 */
                "vmulps	%%xmm10, %%xmm10, %%xmm10 \n\t"		    /* approx 1/4 */
                "vmulps	%%xmm10, %%xmm10, %%xmm10 \n\t"		    /* approx 1/16 */
                "vmulps	%%xmm10, %%xmm10, %%xmm10 \n\t"		    /* approx 1/256 */
                "vmulps	%%xmm10, %%xmm10, %%xmm10 \n\t"		    /* approx 1/65536 */
                "vcvtps2pd %%xmm10, %%ymm10 \n\t"	            /* convert to double precision */
                "vmovdqu    (%8), %%ymm11 \n\t"                 /* load conversion constant into ymm11 */
                "movl       %0, %%r10d \n\t"
                "movl       %7, %%r11d \n\t"
                "decl       %%r11d \n\t"

                /* top of while(1) loop */
                "0: \n\t"       

                /* if (i >= imax) do this then return */
                "cmpl       %%r10d, %%r11d \n\t"
                "jae 1f \n\t"
                "vmovdqa    %%ymm0, (%1)  \n\t" /* P */
                "vmovdqa    %%ymm1, (%2)  \n\t" /* Qn */
                "vmovdqa    %%ymm2, (%3)  \n\t" /* Q0 */
                "vmovdqa    %%ymm3, (%4)  \n\t" /* bn */
                "jmp 5f \n\t"                   /* jump out of loop and set return condition */

                "1: \n\t"
                /* even iteration */
                "vpmulld    %%ymm3, %%ymm1, %%ymm6 \n\t"    /* P = bn * Qn - P */
                "vpsubd     %%ymm0, %%ymm6, %%ymm6 \n\t"
                "vpsubd     %%ymm6, %%ymm0, %%ymm7  \n\t"   /* Qn = Q0 + bn * (t1 - P) */
                "vmovdqa    %%ymm6, %%ymm0 \n\t"            /* now ok to write new P */
                "vpmulld    %%ymm7, %%ymm3, %%ymm7  \n\t"
                "vpaddd     %%ymm2, %%ymm7, %%ymm7  \n\t"
                "vmovdqa    %%ymm1, %%ymm2 \n\t"            /* now ok to write new Q0 */
                "vmovdqa    %%ymm7, %%ymm1 \n\t"            /* now ok to write new Qn */
                "vpaddd     %%ymm4, %%ymm0, %%ymm6  \n\t"    /* bn = (b0 + P) / Qn */

#ifdef DEFINED //NOTDEFINED // DEFINED //
                /* first 4 divisions */
                "vpmovzxdq  %%xmm1, %%ymm7 \n\t"            /* (p5, L3) zero extend 32-bit to 64-bit */
                "vpmovzxdq  %%xmm6, %%ymm9 \n\t"            /* (p5, L3) zero extend 32-bit to 64-bit */
                "vpaddq     %%ymm7, %%ymm11, %%ymm7 \n\t"   /* (p1, L3) add magic constant as integer */
                "vextracti128   $1,%%ymm1,%%xmm12 \n\t"     /* (p5, L3) grab high half of inputs */
                "vpaddq     %%ymm9, %%ymm11, %%ymm9 \n\t"   /* (p1, L3) add magic constant as integer */
                "vextracti128   $1,%%ymm6,%%xmm13 \n\t"     /* (p5, L3) grab high half of inputs */
                "vsubpd     %%ymm11, %%ymm7, %%ymm7 \n\t"   /* (p1, L3) sub magic constant as double */
                "vpmovzxdq  %%xmm12, %%ymm12 \n\t"          /* (p5, L3) zero extend 32-bit to 64-bit */
                "vsubpd     %%ymm11, %%ymm9, %%ymm9 \n\t"   /* (p1, L3) sub magic constant as double */
                "vpmovzxdq  %%xmm13, %%ymm13 \n\t"          /* (p5, L3) zero extend 32-bit to 64-bit */
                "vpaddq     %%ymm12, %%ymm11, %%ymm12 \n\t" /* (p1, L3) add magic constant as integer */
                
                "vdivpd     %%ymm7, %%ymm9, %%ymm9 \n\t"
                "vpaddq     %%ymm13, %%ymm11, %%ymm13 \n\t" /* (p1, L3) add magic constant as integer */      
                "vsubpd     %%ymm11, %%ymm12, %%ymm12 \n\t" /* (p1, L3) sub magic constant as double */
                "vsubpd     %%ymm11, %%ymm13, %%ymm13 \n\t" /* (p1, L3) sub magic constant as double */
                "vcvttpd2dq	%%ymm9, %%xmm3 \n\t"	        /* (p1,p5, L6) truncate to uint32 */
                
                /* second 4 divisions */                
                "vdivpd     %%ymm12, %%ymm13, %%ymm13 \n\t"
                "vcvttpd2dq	%%ymm13, %%xmm13 \n\t"	        /* truncate to uint32 */
                "vinserti128    $1, %%xmm13, %%ymm3, %%ymm3 \n\t" /* insert into high half of bn */
#else

                // numerator is in ymm6
                // denominator is in ymm1
                // result needs to be in ymm3 at the end
                // used (unavailable) registers:
                // ymm0 (P), ymm1 (Qn), ymm2 (Q0), ymm4 (b0), ymm5 (double 2's), 
                // ymm10 (fudge), ymm11 (convert)
                // available registers:
                // ymm6 (after numerator is used)
                // ymm3, ymm7, ymm8, ymm9, ymm12, ymm13, ymm14, ymm15                

                "vextracti128   $1,%%ymm1,%%xmm13 \n\t"      /* grab high half of inputs */
                "vpmovzxdq  %%xmm1, %%ymm12 \n\t"            /* zero extend 32-bit to 64-bit */
                "vpmovzxdq  %%xmm13, %%ymm13 \n\t"            /* zero extend 32-bit to 64-bit */
                "vpaddq     %%ymm12, %%ymm11, %%ymm12 \n\t"   /* add magic constant as integer */
                "vpaddq     %%ymm13, %%ymm11, %%ymm13 \n\t"   /* add magic constant as integer */
                "vsubpd     %%ymm11, %%ymm12, %%ymm12 \n\t"   /* sub magic constant as double */
                "vsubpd     %%ymm11, %%ymm13, %%ymm13 \n\t"   /* sub magic constant as double */

                "vcvtpd2ps	    %%ymm12, %%xmm7 \n\t"	            /* convert double to float */    
                "vrcpps	        %%xmm7, %%xmm7 \n\t"	            /* approx 1/d */
                "vcvtpd2ps	    %%ymm13, %%xmm9 \n\t"	            /* convert double to float */
                "vrcpps	        %%xmm9, %%xmm9 \n\t"	            /* approx 1/d */ 
                "vcvtps2pd      %%xmm7, %%ymm7 \n\t"	            /* convert back to double precision */
                "vcvtps2pd      %%xmm9, %%ymm9 \n\t"	            /* convert back to double precision */
                

                "vmulpd	        %%ymm7, %%ymm12, %%ymm14 \n\t"	      /* d * 1/d */
                "vpmovzxdq      %%xmm6, %%ymm8 \n\t"                /* extend low half of numerator */
                "vmulpd	        %%ymm9, %%ymm13, %%ymm15 \n\t"	      /* d * 1/d */
                "vextracti128   $1,%%ymm6,%%xmm6 \n\t"              /* grab high half of numerator */
                "vmulpd	        %%ymm14, %%ymm7, %%ymm14 \n\t"	      /* d * 1/d * 1/d */
                "vpmovzxdq      %%xmm6, %%ymm6 \n\t"                /* extend high half of numerator */
                "vmulpd	        %%ymm15, %%ymm9, %%ymm15 \n\t"	      /* d * 1/d * 1/d */
                "vpaddq         %%ymm8, %%ymm11, %%ymm8 \n\t"       /* add magic constant as integer */
                "vfmsub132pd	%%ymm5, %%ymm14, %%ymm7 \n\t"	/* 1/d = 2 * 1/d - d * 1/d^2 */
                "vpaddq         %%ymm6, %%ymm11, %%ymm6 \n\t"       /* add magic constant as integer */
                "vfmsub132pd	%%ymm5, %%ymm15, %%ymm9 \n\t"	/* 1/d = 2 * 1/d - d * 1/d^2 */
                "vsubpd         %%ymm11, %%ymm8, %%ymm8 \n\t"   /* sub magic constant as double */
                "vmulpd	        %%ymm7, %%ymm12, %%ymm14 \n\t"	      /* d * 1/d */
                "vsubpd         %%ymm11, %%ymm6, %%ymm6 \n\t"   /* sub magic constant as double */
                "vmulpd	        %%ymm9, %%ymm13, %%ymm15 \n\t"	      /* d * 1/d */
                "vmulpd	        %%ymm14, %%ymm7, %%ymm14 \n\t"	      /* d * 1/d * 1/d */
                "vmulpd	        %%ymm15, %%ymm9, %%ymm15 \n\t"	      /* d * 1/d * 1/d */
                "vfmsub132pd	%%ymm5, %%ymm14, %%ymm7 \n\t"	/* 1/d = 2 * 1/d - d * 1/d^2 */
                "vfmsub132pd	%%ymm5, %%ymm15, %%ymm9 \n\t"	/* 1/d = 2 * 1/d - d * 1/d^2 */
                "vmovapd        %%ymm8, %%ymm14 \n\t"           /* copy numerator */
                "vmovapd        %%ymm6, %%ymm15 \n\t"           /* copy numerator */
                "vfmadd132pd	%%ymm7, %%ymm10, %%ymm14 \n\t"	/* n * 1/d + 1/65536 */
                "vfmadd132pd	%%ymm9, %%ymm10, %%ymm15 \n\t"	/* n * 1/d + 1/65536 */
                "vroundpd	    $3, %%ymm14, %%ymm14 \n\t"	/* truncate */
                "vroundpd	    $3, %%ymm15, %%ymm15 \n\t"	/* truncate */
                "vmulpd	        %%ymm14, %%ymm12, %%ymm12 \n\t"	      /* ans * denominators */
                "vmulpd	        %%ymm15, %%ymm13, %%ymm13 \n\t"	      /* ans * denominators */
                "vcmppd         $0xe, %%ymm8, %%ymm12, %%ymm7 \n\t"   /* if (ymm4[j] > numerator[j]) set ymm5[j] = 0xffffffffffffffff, for j=0:3 */
                "vcmppd         $0xe, %%ymm6, %%ymm13, %%ymm9 \n\t"   /* if (ymm4[j] > numerator[j]) set ymm5[j] = 0xffffffffffffffff, for j=0:3 */
                "vpand          %%ymm10, %%ymm7, %%ymm7 \n\t"            
                "vpand          %%ymm10, %%ymm9, %%ymm9 \n\t"
                "vsubpd         %%ymm7, %%ymm14, %%ymm14 \n\t"
                "vsubpd         %%ymm9, %%ymm15, %%ymm15 \n\t"
                "vcvttpd2dq	    %%ymm14, %%xmm3 \n\t"	        /* truncate to uint32 */
                "vcvttpd2dq	    %%ymm15, %%xmm7 \n\t"	        /* truncate to uint32 */
                "vinserti128   $1,%%xmm7, %%ymm3, %%ymm3 \n\t"      /* set high half of bn */
#endif

                "incl       %%r10d \n\t"

                /* square root check */
                "xorl       %%r8d, %%r8d \n\t"

                "vextracti128   $1,%%ymm1,%%xmm12 \n\t"      /* grab high half of inputs */
                "vpmovzxdq  %%xmm1, %%ymm6 \n\t"            /* zero extend 32-bit to 64-bit */
                "vmovdqa    %%xmm12, %%xmm13 \n\t"          /* copy high half */
                "vpaddq     %%ymm6, %%ymm11, %%ymm6 \n\t"   /* add magic constant as integer */
                "vpmovzxdq  %%xmm12, %%ymm12 \n\t"            /* zero extend 32-bit to 64-bit */
                "vsubpd     %%ymm11, %%ymm6, %%ymm6 \n\t"   /* sub magic constant as double */                

                "vsqrtpd    %%ymm6, %%ymm6	        	     \n\t"      /* take the square root */
                "vpaddq     %%ymm12, %%ymm11, %%ymm12 \n\t"   /* add magic constant as integer */
                "vsubpd     %%ymm11, %%ymm12, %%ymm12 \n\t"   /* sub magic constant as double */
                "vcvttpd2dq %%ymm6, %%xmm6 \n\t"
                "vsqrtpd    %%ymm12, %%ymm12	        	     \n\t"      /* take the square root */
                "vpmulld    %%xmm6, %%xmm6, %%xmm6	         \n\t"      /* multiply answer * answer */
                "vcvttpd2dq %%ymm12, %%xmm12 \n\t"
                "vpcmpeqd   %%xmm6, %%xmm1, %%xmm7           \n\t"      /* see if equal to starting value */
                "vpmulld    %%xmm12, %%xmm12, %%xmm12	         \n\t"      /* multiply answer * answer */
                "vmovdqa    %%xmm7, %%xmm9 \n\t"                        /* save the success mask */
                "vpcmpeqd   %%xmm12, %%xmm13, %%xmm13           \n\t"      /* see if equal to starting value */
                "vmovmskps  %%xmm7, %%r8d	                 \n\t"      /* set a bit if so */
                "vinserti128    $1, %%xmm13, %%ymm9, %%ymm9 \n\t"        /* insert into high half of success mask */
                "vmovmskps  %%xmm13, %%r9d	                 \n\t"      /* set a bit if so */
                "vpsrld     $31, %%ymm9, %%ymm6 \n\t"                   /* shift the mask to 1's in positions that pass */
                "orl        %%r9d, %%r8d \n\t"

                /* if (success) do this then exit loop */
                "testl      %%r8d, %%r8d \n\t"
                "jz 2f \n\t"
                "vmovdqa    %%ymm0, (%1)  \n\t" /* P */
                "vmovdqa    %%ymm1, (%2)  \n\t" /* Qn */
                "vmovdqa    %%ymm2, (%3)  \n\t" /* Q0 */
                "vmovdqa    %%ymm3, (%4)  \n\t" /* bn */
                "vmovdqa    %%ymm6, (%6)  \n\t" /* successes */
                "jmp 5f \n\t"                   /* jump out of loop and set break condition */

                "2: \n\t"
                /* odd iteration */ 
                "vpmulld    %%ymm3, %%ymm1, %%ymm6 \n\t"    /* P = bn * Qn - P */
                "vpsubd     %%ymm0, %%ymm6, %%ymm6 \n\t"
                "vpsubd     %%ymm6, %%ymm0, %%ymm7  \n\t"   /* Qn = Q0 + bn * (t1 - P) */
                "vmovdqa    %%ymm6, %%ymm0 \n\t"            /* now ok to write new P */
                "vpmulld    %%ymm7, %%ymm3, %%ymm7  \n\t"
                "vpaddd     %%ymm2, %%ymm7, %%ymm7  \n\t"
                "vmovdqa    %%ymm1, %%ymm2 \n\t"            /* now ok to write new Q0 */
                "vmovdqa    %%ymm7, %%ymm1 \n\t"            /* now ok to write new Qn */
                "vpaddd     %%ymm4, %%ymm0, %%ymm6  \n\t"    /* bn = (b0 + P) / Qn */

#ifdef DEFINED //NOTDEFINED // DEFINED //
                /* first 4 divisions */
                "vpmovzxdq  %%xmm1, %%ymm7 \n\t"            /* (p5, L3) zero extend 32-bit to 64-bit */
                "vpmovzxdq  %%xmm6, %%ymm9 \n\t"            /* (p5, L3) zero extend 32-bit to 64-bit */
                "vpaddq     %%ymm7, %%ymm11, %%ymm7 \n\t"   /* (p1, L3) add magic constant as integer */
                "vextracti128   $1,%%ymm1,%%xmm12 \n\t"     /* (p5, L3) grab high half of inputs */
                "vpaddq     %%ymm9, %%ymm11, %%ymm9 \n\t"   /* (p1, L3) add magic constant as integer */
                "vextracti128   $1,%%ymm6,%%xmm13 \n\t"     /* (p5, L3) grab high half of inputs */
                "vsubpd     %%ymm11, %%ymm7, %%ymm7 \n\t"   /* (p1, L3) sub magic constant as double */
                "vpmovzxdq  %%xmm12, %%ymm12 \n\t"          /* (p5, L3) zero extend 32-bit to 64-bit */
                "vsubpd     %%ymm11, %%ymm9, %%ymm9 \n\t"   /* (p1, L3) sub magic constant as double */
                "vpmovzxdq  %%xmm13, %%ymm13 \n\t"          /* (p5, L3) zero extend 32-bit to 64-bit */
                "vpaddq     %%ymm12, %%ymm11, %%ymm12 \n\t" /* (p1, L3) add magic constant as integer */

                "vdivpd     %%ymm7, %%ymm9, %%ymm9 \n\t"
                "vpaddq     %%ymm13, %%ymm11, %%ymm13 \n\t" /* (p1, L3) add magic constant as integer */   
                "vsubpd     %%ymm11, %%ymm12, %%ymm12 \n\t" /* (p1, L3) sub magic constant as double */
                "vsubpd     %%ymm11, %%ymm13, %%ymm13 \n\t" /* (p1, L3) sub magic constant as double */
                "vcvttpd2dq	%%ymm9, %%xmm3 \n\t"	        /* (p1,p5, L6) truncate to uint32 */                

                /* second 4 divisions */                
                "vdivpd     %%ymm12, %%ymm13, %%ymm13 \n\t"
                "vcvttpd2dq	%%ymm13, %%xmm13 \n\t"	        /* truncate to uint32 */
                "vinserti128    $1, %%xmm13, %%ymm3, %%ymm3 \n\t" /* insert into high half of bn */

#else

                // numerator is in ymm6
                // denominator is in ymm1
                // result needs to be in ymm3 at the end
                // used (unavailable) registers:
                // ymm0 (P), ymm1 (Qn), ymm2 (Q0), ymm4 (b0), ymm5 (double 2's), 
                // ymm10 (fudge), ymm11 (convert)
                // available registers:
                // ymm6 (after numerator is used)
                // ymm3, ymm7, ymm8, ymm9, ymm12, ymm13, ymm14, ymm15                

                "vextracti128   $1,%%ymm1,%%xmm13 \n\t"      /* grab high half of inputs */
                "vpmovzxdq  %%xmm1, %%ymm12 \n\t"            /* zero extend 32-bit to 64-bit */
                "vpmovzxdq  %%xmm13, %%ymm13 \n\t"            /* zero extend 32-bit to 64-bit */
                "vpaddq     %%ymm12, %%ymm11, %%ymm12 \n\t"   /* add magic constant as integer */
                "vpaddq     %%ymm13, %%ymm11, %%ymm13 \n\t"   /* add magic constant as integer */
                "vsubpd     %%ymm11, %%ymm12, %%ymm12 \n\t"   /* sub magic constant as double */
                "vsubpd     %%ymm11, %%ymm13, %%ymm13 \n\t"   /* sub magic constant as double */

                "vcvtpd2ps	    %%ymm12, %%xmm7 \n\t"	            /* convert double to float */
                "vpmovzxdq      %%xmm6, %%ymm8 \n\t"                /* extend low half of numerator */
                "vcvtpd2ps	    %%ymm13, %%xmm9 \n\t"	            /* convert double to float */
                "vextracti128   $1,%%ymm6,%%xmm6 \n\t"              /* grab high half of numerator */
                "vrcpps	        %%xmm7, %%xmm7 \n\t"	            /* approx 1/d */
                "vpmovzxdq      %%xmm6, %%ymm6 \n\t"                /* extend high half of numerator */
                "vpaddq         %%ymm8, %%ymm11, %%ymm8 \n\t"       /* add magic constant as integer */
                "vrcpps	        %%xmm9, %%xmm9 \n\t"	            /* approx 1/d */
                "vpaddq         %%ymm6, %%ymm11, %%ymm6 \n\t"       /* add magic constant as integer */
                "vcvtps2pd      %%xmm7, %%ymm7 \n\t"	            /* convert back to double precision */
                "vsubpd         %%ymm11, %%ymm8, %%ymm8 \n\t"   /* sub magic constant as double */
                "vcvtps2pd      %%xmm9, %%ymm9 \n\t"	            /* convert back to double precision */
                "vsubpd         %%ymm11, %%ymm6, %%ymm6 \n\t"   /* sub magic constant as double */

                "vmulpd	        %%ymm7, %%ymm12, %%ymm14 \n\t"	      /* d * 1/d */
                "vmulpd	        %%ymm9, %%ymm13, %%ymm15 \n\t"	      /* d * 1/d */
                "vmulpd	        %%ymm14, %%ymm7, %%ymm14 \n\t"	      /* d * 1/d * 1/d */
                "vmulpd	        %%ymm15, %%ymm9, %%ymm15 \n\t"	      /* d * 1/d * 1/d */
                "vfmsub132pd	%%ymm5, %%ymm14, %%ymm7 \n\t"	/* 1/d = 2 * 1/d - d * 1/d^2 */
                "vfmsub132pd	%%ymm5, %%ymm15, %%ymm9 \n\t"	/* 1/d = 2 * 1/d - d * 1/d^2 */
                "vmulpd	        %%ymm7, %%ymm12, %%ymm14 \n\t"	      /* d * 1/d */
                "vmulpd	        %%ymm9, %%ymm13, %%ymm15 \n\t"	      /* d * 1/d */
                "vmulpd	        %%ymm14, %%ymm7, %%ymm14 \n\t"	      /* d * 1/d * 1/d */
                "vmulpd	        %%ymm15, %%ymm9, %%ymm15 \n\t"	      /* d * 1/d * 1/d */
                "vfmsub132pd	%%ymm5, %%ymm14, %%ymm7 \n\t"	/* 1/d = 2 * 1/d - d * 1/d^2 */
                "vfmsub132pd	%%ymm5, %%ymm15, %%ymm9 \n\t"	/* 1/d = 2 * 1/d - d * 1/d^2 */
                "vmovapd        %%ymm8, %%ymm14 \n\t"           /* copy numerator */
                "vmovapd        %%ymm6, %%ymm15 \n\t"           /* copy numerator */
                "vfmadd132pd	%%ymm7, %%ymm10, %%ymm14 \n\t"	/* n * 1/d + 1/65536 */
                "vfmadd132pd	%%ymm9, %%ymm10, %%ymm15 \n\t"	/* n * 1/d + 1/65536 */
                "vroundpd	    $3, %%ymm14, %%ymm14 \n\t"	/* truncate */
                "vroundpd	    $3, %%ymm15, %%ymm15 \n\t"	/* truncate */
                "vmulpd	        %%ymm14, %%ymm12, %%ymm12 \n\t"	      /* ans * denominators */
                "vmulpd	        %%ymm15, %%ymm13, %%ymm13 \n\t"	      /* ans * denominators */
                "vcmppd         $0xe, %%ymm8, %%ymm12, %%ymm7 \n\t"   /* if (ymm4[j] > numerator[j]) set ymm5[j] = 0xffffffffffffffff, for j=0:3 */
                "vcmppd         $0xe, %%ymm6, %%ymm13, %%ymm9 \n\t"   /* if (ymm4[j] > numerator[j]) set ymm5[j] = 0xffffffffffffffff, for j=0:3 */
                "vpand          %%ymm10, %%ymm7, %%ymm7 \n\t"
                "vpand          %%ymm10, %%ymm9, %%ymm9 \n\t"
                "vsubpd         %%ymm7, %%ymm14, %%ymm14 \n\t"
                "vsubpd         %%ymm9, %%ymm15, %%ymm15 \n\t"
                "vcvttpd2dq	    %%ymm14, %%xmm3 \n\t"	        /* truncate to uint32 */
                "vcvttpd2dq	    %%ymm15, %%xmm7 \n\t"	        /* truncate to uint32 */
                "vinserti128   $1,%%xmm7, %%ymm3, %%ymm3 \n\t"      /* set high half of bn */
#endif

                "incl       %%r10d \n\t"

                "jmp    0b \n\t"

                /* exit marker */
                "5: \n\t"
                "movl   %%r10d, %0 \n\t"            /* export iteration count */

                : "+r"(i)
                : "r"(P), "r"(Qn), "r"(Q0), "r"(bn), "r"(b0),
                "r"(success_vec), "r"(imax), "r"(conversion)
                : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm9", "xmm10", "xmm11", 
                    "xmm12", "xmm13", "xmm14", "xmm15", "r8", "r9", "r10", "r11", "cc", "memory"
                );


                kk = i - kk;

#ifdef PRINT_DEBUG

                printf("values after i = %d\n", i);
                printf("P = [");
                for (k = 0; k < NUM_LANES; k++)
                    printf("%u ", P[k]);
                printf("]\n");

                printf("Qn = [");
                for (k = 0; k < NUM_LANES; k++)
                    printf("%u ", Qn[k]);
                printf("]\n");

                printf("Q0 = [");
                for (k = 0; k < NUM_LANES; k++)
                    printf("%u ", Q0[k]);
                printf("]\n");

                printf("bn = [");
                for (k = 0; k < NUM_LANES; k++)
                    printf("%u ", bn[k]);
                printf("]\n");

                printf("b0 = [");
                for (k = 0; k < NUM_LANES; k++)
                    printf("%u ", b0[k]);
                printf("]\n");

                printf("iterations = [");
                for (k = 0; k < NUM_LANES; k++)
                    printf("%u ", iterations[k]);
                printf("]\n\n");
#endif

        //}

        //printf("loop exited at iteration %d (of %d)\n", i, imax);
            
        for (k = 0; k < NUM_LANES; k++)
        {
            // set default signal: to do nothing but continue looking for a factor
            mult_save->f[k] = 0;
            iterations[k] += kk;
        }

        if (i >= imax)
        {

#ifdef PRINT_DEBUG
            printf("i > imax\n");
            printf("ending values:\n");

            printf("P = [");
            for (k = 0; k < NUM_LANES; k++)
                printf("%u ", P[k]);
            printf("]\n");

            printf("Qn = [");
            for (k = 0; k < NUM_LANES; k++)
                printf("%u ", Qn[k]);
            printf("]\n");

            printf("Q0 = [");
            for (k = 0; k < NUM_LANES; k++)
                printf("%u ", Q0[k]);
            printf("]\n");

            printf("bn = [");
            for (k = 0; k < NUM_LANES; k++)
                printf("%u ", bn[k]);
            printf("]\n");

            printf("b0 = [");
            for (k = 0; k < NUM_LANES; k++)
                printf("%u ", b0[k]);
            printf("]\n");

            printf("iterations = [");
            for (k = 0; k < NUM_LANES; k++)
                printf("%u ", iterations[k]);
            printf("]\n\n");

#endif

            return;
        }

#ifdef PRINT_DEBUG
        printf("found a square\n");
        printf("ending values:\n");

        printf("P = [");
        for (k = 0; k < NUM_LANES; k++)
            printf("%u ", P[k]);
        printf("]\n");

        printf("Qn = [");
        for (k = 0; k < NUM_LANES; k++)
            printf("%u ", Qn[k]);
        printf("]\n");

        printf("Q0 = [");
        for (k = 0; k < NUM_LANES; k++)
            printf("%u ", Q0[k]);
        printf("]\n");

        printf("bn = [");
        for (k = 0; k < NUM_LANES; k++)
            printf("%u ", bn[k]);
        printf("]\n");

        printf("b0 = [");
        for (k = 0; k < NUM_LANES; k++)
            printf("%u ", b0[k]);
        printf("]\n");

        printf("iterations = [");
        for (k = 0; k < NUM_LANES; k++)
            printf("%u ", iterations[k]);
        printf("]\n");

        printf("success_vec = [");
        for (k = 0; k < NUM_LANES; k++)
            printf("%u ", success_vec[k]);
        printf("]\n\n");

#endif

        success = 0;
        for (k = 0; k < NUM_LANES; k++)
        {
            if ((success_vec[k] == 1) && (mult_save->active[k] == 1))
            {
                int failsafe = 10000;

                //printf("symmetry check on input %lu in lane %d after %d iterations\n", 
                //    mult_save->mN[k], k, i);

                // reduce to G0
                S[k] = (int)sqrt(Qn[k]);
                Ro[k] = P[k] + S[k] * ((b0[k] - P[k]) / S[k]);
                t1[k] = Ro[k];
                So[k] = (uint32)(((int64)mult_save->mN[k] - (int64)t1[k] * (int64)t1[k]) / (int64)S[k]);
                bbn[k] = (b0[k] + Ro[k]) / So[k];

                if ((S[k] *  S[k]) != Qn[k])
                {
                    printf("error, reported square is not square: %u, %u\n", S[k], Qn[k]); 
                    fflush(stdout);
                    continue;                    
                }
                //printf("Searching for symmetry on Square %u, S = %u (%u)...", Qn[k], S[k], S[k] * S[k]);

                // search for symmetry point
                while (failsafe)
                {

                    t1[k] = Ro[k];
                    Ro[k] = bbn[k] * So[k] - Ro[k];
                    t2[k] = So[k];
                    So[k] = S[k] + bbn[k] * (t1[k] - Ro[k]);
                    S[k] = t2[k];
                    bbn[k] = (b0[k] + Ro[k]) / So[k];

                    // eliminate trivial cases before running the expensive gcd
                    if (Ro[k] == t1[k])
                    {
                        // eliminate trivial cases before running the expensive gcd
                        if (Ro[k] > mult_save->mult[k])
                        {
                            // these run at almost exactly the same speed.  use the normal
                            // gcd both because it is more general and because it will
                            // pick up any future advances in hardware integer division.
                            mult_save->f[k] = gcd64(Ro[k], mult_save->N[k]);
                            //mult_save->f[k] = spBinGCD(mult_save->N[k], Ro[k]);

                            if (mult_save->f[k] > 1)
                            {
                                success = 1;
                            }
                        }
                        break;
                    }

                    
                    failsafe--;
                }

                //if (success)
                //    printf("found result %u\n", mult_save->f[k]);
                //else
                //    printf("\n");
                
            }
        }


        if (success)
        {
            // at least one of the lanes was successful.  drop out of this
            // routine so we can replace that lane.  else we'll keep
            // going until we hit imax (or another square).
            return;
        }
    }
}

#endif

/****
Code by Warren D. Smith, Dec 2011, to implement Lehman integer-factorization
algorithm.  Rigorous O(N^(1/3)) step factoring
(or prime-proving) algorithm invented by RS Lehman 1974.
Compile:
 gcc LehmanClean.c -Wall -O6 -o lehcl
 gcc LehmanClean.c -Wall -DNDEBUG -O6 -o lehcl
*******/

static double sqr_tab[1024]; 

void make_sqr_tab()
{
	int i;
	for(i=0; i<1024; i++)
		sqr_tab[i] = sqrt((double)i);

	return;
}

static unsigned char issq1024[1024];
static unsigned char issq4199[4199];

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 1
#endif

void MakeIssq()
{
	int i;
	for(i=0; i<1024; i++){ issq1024[i] = 0; }
	for(i=0; i<1024; i++){ issq1024[(i*i)%1024] = 1; }
	for(i=0; i<4199; i++){ issq4199[i] = 0; }
	for(i=0; i<3465; i++){ issq4199[(i*i)%3465] |= 2; }
	for(i=0; i<4199; i++){ issq4199[(i*i)%4199] |= 1; }
	//printf("Square tables built.\n");
}

//the 6542 primes up to 65536=2^16, then sentinel 65535 at end
static uint16 prime[6543]; 

void MakePrimeTable(){
	uint32 i,j,k;
	prime[0]=2;
	prime[1]=3;
	prime[2]=5;
	k=3;
	for(i=7; i<65536; i+=2){
		for(j=0; prime[j]*(uint32)prime[j] <= i; j++){
			if(i%prime[j]==0) goto NONPRIME;
		}
		prime[k]=i;
		k++;
		NONPRIME: ;
	}

	prime[k] = 65535; //sentinel
	//printf("Prime table[0..%d] built: ", k);
	//for(i=0; i<20; i++){ printf("%d,", prime[i]); }
	//printf("%d,...,%d,(%d)\n", prime[20],prime[6541],prime[6542]);
}

/*** LehmanFactor  Testing:
Experimentally: Always returns N if N=prime<65536 and DoTrial=FALSE
and 0.1<=LehTune<=9.6.
  But: If LehTune<=0.8 then the   gcd(a+b, N)   can return N as its
mode of success if N=3.
  So the theorem that the special prime-detect return always will
fire for N prime can thus
  be violated if LehTune<0.8 and N=3 (still correctly returns N though), but
       valid if LehTune>=0.8  and N>2.
  Also valid if LehTune>=0.25 and N>31.
  Also valid if LehTune>=0.10 and N>83.

High end: if N < 2^43.000049 = 8796393022207  then apparently no
overflow issues.
Later note: sorry, the simple search loop I was using was too
simplistic a search approach;
The currently least known overflow-failure also occurs at   N =
5132209842943 = 2^42.2227.

Factor Validity: testing 100 million N with N=A*B with
 2611953 <= A <= B <= 2965920
LehmanFactor always succeeded in factoring.
However, this is deceptive, see above on weird rare failures.  It
appears that Lehman
can be sensitive to the distribution of N and there are "hard locales"
and "easy locales" which
are not obvious.  Also changing the tuning constants can alter the failure set.
***/

uint64 LehmanFactor(uint64 N, double Tune, int DoTrial, double CutFrac)
{
	uint32 b,p,k,r,B,U,Bred,inc,FirstCut,ip = 1;
	uint64 a,c,kN,kN4,B2;
	double Tune2, Tune3, x, sqrtn;
	mpz_t tmpz;

	if((N&1)==0) 
		return(2); //N is even

	if(Tune<0.1)
	{
		printf("Sorry0, Lehman only implemented for Tune>=0.1\n");
		return(1);
	}

	mpz_init(tmpz);
	mpz_set_64(tmpz, N);
	mpz_root(tmpz, tmpz, 3);
	B = Tune * (1+(double)mpz_get_ui(tmpz));

	FirstCut = CutFrac*B;

	//assures prime N will not activate "wrong" Lehman return
	if(FirstCut<84)
		FirstCut=84; 

	if(FirstCut>65535)
		FirstCut = 65535; 

	if(DoTrial)
	{
		for(ip=1; ; ip++)
		{                   
			p = prime[ip];
			if(p>=FirstCut) 
				break;
            if (N%p == 0)
            {
                mpz_clear(tmpz);
                return(p);
            }
		}
	}

	if(N>=8796393022207ull)
	{
		printf("Sorry1, Lehman only implemented for N<8796393022207\n");
        mpz_clear(tmpz);
		return(1);
	}

	Tune2 = Tune*Tune;
	Tune3 = Tune2*Tune;
	Bred = B / Tune3;

	B2 = B*B;
	kN = 0;

	//Lehman suggested (to get more average speed) trying highly-divisible k first. However,
	//my experiments on trying to to that have usually slowed things down versus this simple loop:
	sqrtn = sqrt((double)N);
	for(k=1; k<=Bred; k++)
	{
		if(k&1)
		{ 
			inc=4; 
			r=(k+N)%4; 
		}
		else
		{ 
			inc=2; 
			r=1; 
		} 

		kN += N;
		if(kN >= 1152921504606846976ull)
		{
			printf("Sorry2, overflow, N=%" PRIu64 " is too large\n", N);
            mpz_clear(tmpz);
			return(1);
		}

		//Actually , even if overflow occurs here, one could still use approximate
		//arithmetic to compute kN4's most-signif 64 bits only, then still exactly compute x...
		//With appropriate code alterations is should be possible to extent the range... but
		//I have not tried that idea for trying to "cheat" to gain more precision.
		kN4 = kN*4;
		if (k < 1024)
			x = sqrtn * sqr_tab[k];
		else
			x = sqrt((double)kN);

		a = x;
		if((uint64)a*(uint64)a==kN)
		{ 
			B2 = gcd64((uint64)a, N);
            mpz_clear(tmpz);
			return(B2);
		}

		x *= 2;
		a = x+0.9999999665; //very carefully chosen.
		//Let me repeat that: a = x+0.9999999665.  Really.
		b=a%inc;  
		b = a + (inc+r-b)%inc;   //b is a but adjusted upward to make b%inc=r.
		c = (uint64)b*(uint64)b - kN4;  //this is the precision bottleneck.
		//At this point, I used to do a test:
		//if( c+kN4 != (uint64)b*(uint64)b ) //overflow-caused failure: exit!
		//	printf("Sorry3, unrepairable overflow, N=%llu is too large\n", N);
		//  return(0);
		//However, I've now reconsidered.  I claim C language computes c mod 2^64 correctly.
		//If overflow happens, this is irrelevant because c is way smaller than 2^64, kN4, and b*b.
		//Hence c should be correct, despite the overflow. Hence I am removing this error-exit.

		U = x + B2/(2*x); 
		//old code was  U = sqrt((real)(B2+kN4+0.99));   and was 4% slower.

		//Below loop is: for(all integers a with 0<=a*a-kN4<=B*B and with a%inc==r)
		for(a=b;  a<=U;  c+=inc*(a+a+inc), a+=inc )
		{
			//again, even though this assert can fail due to overflow, that overflow should not matter:
			/** Programming trick:    naive code:     c = a*a-kN4;
			In the inner loop c is bounded between 0 and T^2*N^(2/3)
			and can be updated additively by adding inc*(anew+aold) to it
			when we update a to anew=a+inc. This saves a multiplication in
			the inner loop and/or allows us to reduce precision. **/
			if(issq1024[c&1023])
			{
				if(issq4199[c%3465]&2)
				{
					if(issq4199[c%4199]&1)
					{
						b = sqrt(c + 0.9);
						if(b*b==c)
						{ 
							//square found
							B2 = gcd64((uint64)(a+b), N);
							if(B2>=N)
								printf("theorem failure: B2=%" PRIu64 " N=%" PRIu64 "\n", B2,N); 
                            mpz_clear(tmpz);
							return(B2);
						}
					}
				}
			}
		}
	}

	//square-finding has failed so resume missing part of trial division: 
	if(DoTrial)
	{ 
		if(B>65535) 
			B = 65535;
		for( ; ; ip++)
		{ 
			p = prime[ip];
			if(p>=B) 
				break;
            if (N%p == 0)
            {
                mpz_clear(tmpz);
                return(p);
            }
		}
	}

    mpz_clear(tmpz);
	return(N); //N is prime
}

void init_lehman()
{
	MakeIssq();
	MakePrimeTable();
	make_sqr_tab();

	return;
}


#ifdef NOTDEF
void VecLehmanFactor(uint64 *N, double Tune, double CutFrac, uint64 *f, uint32 num)
{
    uint32 b, p, k, r, U, B, inc, FirstCut, ip = 1;
    uint32 Bred[NUM_LANES];
    uint64 a, c, kN, kN4;
    uint64 B2[NUM_LANES];
    double Tune2, Tune3, x, sqrtn;
    mpz_t tmpz;
    int vi;
    int num_success;
    int it;

    num_success = 0;
    mpz_init(tmpz);

    for (vi = 0; vi < num; vi += NUM_LANES)
    {

#pragma ivdep
#pragma vector aligned
        for (it = 0; it < NUM_LANES; it++)
        {
            mpz_set_64(tmpz, N[vi+it]);
            mpz_root(tmpz, tmpz, 3);
            B = Tune * (1 + (double)mpz_get_ui(tmpz));

            FirstCut = CutFrac*B;

            //assures prime N will not activate "wrong" Lehman return
            if (FirstCut < 84)
                FirstCut = 84;

            if (FirstCut > 65535)
                FirstCut = 65535;

            Tune2 = Tune*Tune;
            Tune3 = Tune2*Tune;
            Bred[it] = B / Tune3;

            B2[it] = B*B;
            kN = 0;

            //Lehman suggested (to get more average speed) trying highly-divisible k first. However,
            //my experiments on trying to to that have usually slowed things down versus this simple loop:
            sqrtn = sqrt((double)N);
            for (k = 1; k <= Bred; k++)
            {
                if (k & 1)
                {
                    inc = 4;
                    r = (k + N) % 4;
                }
                else
                {
                    inc = 2;
                    r = 1;
                }

                kN += N;
                kN4 = kN * 4;
                if (k < 1024)
                    x = sqrtn * sqr_tab[k];
                else
                    x = sqrt((double)kN);

                a = x;
                if ((uint64)a*(uint64)a == kN)
                {
                    B2 = gcd64((uint64)a, N);
                    return(B2);
                }

                x *= 2;
                a = x + 0.9999999665; //very carefully chosen.

                b = a%inc;
                b = a + (inc + r - b) % inc;   //b is a but adjusted upward to make b%inc=r.
                c = (uint64)b*(uint64)b - kN4;  //this is the precision bottleneck.

                U = x + B2 / (2 * x);

                //Below loop is: for(all integers a with 0<=a*a-kN4<=B*B and with a%inc==r)
                for (a = b; a <= U; c += inc*(a + a + inc), a += inc)
                {
                    b = sqrt(c + 0.9);
                    if (b*b == c)
                    {
                        //square found
                        B2 = gcd64((uint64)(a + b), N);
                        return(B2);
                    }
                }
            }
        }
    }

    return num_success;
}
#endif
