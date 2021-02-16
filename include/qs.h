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

/* structure encapsulating the savefile used in a factorization */
typedef struct {

#if defined(WIN32) || defined(_WIN64)
    HANDLE file_handle;
    uint32 read_size;
    uint32 eof;
#else
    FILE* fp;
#endif
    char* name;
    char* buf;
    uint32 buf_off;
} qs_savefile_t;

typedef struct
{
    mpz_t gmp_n;
    char* savefile_name;
    char* flogname;
    qs_savefile_t savefile;		//savefile object
    FILE* logfile;
    char siqs_savefile[1024];

    double qs_exponent;
    double qs_multiplier;
    double qs_tune_freq;
    int no_small_cutoff_opt;	//1 is true - perform no optimization.  0 = optimize.

    // siqs parameters to override defaults
    uint32_t gbl_override_B;			//override the # of factor base primes
    uint32_t gbl_override_small_cutoff;			//override the tf_small_cutoff value
    uint32_t gbl_override_tf;			//extra reduction of the TF bound by X bits
    uint32_t gbl_override_time;		//stop after this many seconds
    uint32_t gbl_override_rel;		//stop after collecting this many relations
    uint32_t gbl_override_blocks;		//override the # of blocks used
    uint32_t gbl_override_lpmult;		//override the large prime multiplier
    float gbl_override_bdiv;        // override the lpmax divider for batch GCD
    uint32_t gbl_btarget;             // the target number of batch relations
    uint32_t gbl_override_lpb;		// override the large prime bound (specified in bits)
    double gbl_override_mfbt;		// override the mfbt exponent
    double gbl_override_mfbd;		// override the mfbd exponent
    uint32_t gbl_override_3lp_bat;    // don't do 3lp batch factoring (default is do_batch)
    uint32_t inmem_cutoff;          // below X digits, don't use savefile
    int gbl_force_DLP;
    int gbl_force_TLP;

    int gbl_override_B_flag;
    int gbl_override_small_cutoff_flag;
    int gbl_override_tf_flag;
    int gbl_override_time_flag;
    int gbl_override_rel_flag;
    int gbl_override_blocks_flag;
    int gbl_override_lpmult_flag;
    int gbl_override_bdiv_flag;
    
    uint32_t num_factors;			//number of factors found in this method
    uint32_t flags;				//each bit corresponds to a location in the 
                                //flags enum
    double rels_per_sec;
    double qs_time;
    double total_time;

    // configurable options
    int VFLAG;
    int THREADS;
    int LOGFLAG;
    int LATHREADS;
    int NUM_WITNESSES;

    // computer info
    double MEAS_CPU_FREQUENCY;
    char CPU_ID_STR[80];
    int HAS_SSE2;
    int HAS_SSE41;
    int HAS_AVX;
    int HAS_AVX2;
    int L1CACHE;
    int L2CACHE;
    int L3CACHE;
    uint32_t cache_size2;
    uint32_t num_threads;

    uint32_t seed1;
    uint32_t seed2;
    uint64_t lcg_state;

    uint32_t bits;
    uint32_t digits;

} qs_obj_t;



extern void SIQS(qs_obj_t* qsobj);
extern void smallmpqs(qs_obj_t* qsobj);
extern void siqsbench(qs_obj_t* qsobj);

#endif /* _SIQS_H_ */


