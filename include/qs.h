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

#ifndef _SIQS_H_
#define _SIQS_H_

#include <stdint.h>
#include "gmp.h"

typedef struct
{
    mpz_t gmp_n;
    qs_savefile_t savefile;		//savefile object
    char siqs_savefile[1024];

    double qs_exponent;
    double qs_multiplier;
    double qs_tune_freq;
    int no_small_cutoff_opt;	//1 is true - perform no optimization.  0 is false.

    int gbl_override_B_flag;
    uint32_t gbl_override_B;			//override the # of factor base primes
    int gbl_override_small_cutoff_flag;
    uint32_t gbl_override_small_cutoff;			//override the tf_small_cutoff value
    int gbl_override_tf_flag;
    uint32_t gbl_override_tf;			//extra reduction of the TF bound by X bits
    int gbl_override_time_flag;
    uint32_t gbl_override_time;		//stop after this many seconds
    int gbl_override_rel_flag;
    uint32_t gbl_override_rel;		//stop after collecting this many relations
    int gbl_override_blocks_flag;
    uint32_t gbl_override_blocks;		//override the # of blocks used
    int gbl_override_lpmult_flag;
    uint32_t gbl_override_lpmult;		//override the large prime multiplier
    int gbl_override_bdiv_flag;
    float gbl_override_bdiv;        // override the lpmax divider for batch GCD
    uint32_t gbl_btarget;             // the target number of batch relations
    int gbl_force_DLP;
    int gbl_force_TLP;
    uint32_t gbl_override_lpb;		// override the large prime bound (specified in bits)
    double gbl_override_mfbt;		// override the mfbt exponent
    double gbl_override_mfbd;		// override the mfbd exponent
    uint32_t gbl_override_3lp_bat;    // don't do 3lp batch factoring (default is do_batch)

    uint32_t inmem_cutoff;
    uint32_t num_factors;			//number of factors found in this method
    uint32_t flags;				//each bit corresponds to a location in the 
                                //flags enum
    double rels_per_sec;
    double qs_time;
    double total_time;

} qs_obj_t;



extern void SIQS(fact_obj_t* fobj);
extern void smallmpqs(fact_obj_t* fobj);
extern void siqsbench(fact_obj_t* fobj);

#endif /* _SIQS_H_ */


