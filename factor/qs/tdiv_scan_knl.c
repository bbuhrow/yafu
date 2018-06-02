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
#include "qs.h"
#include "factor.h"
#include "util.h"
#include "common.h"

#ifdef USE_AVX512F
#include <immintrin.h>


//#define SIQSDEBUG 1

/*
We are given an array of bytes that has been sieved.  The basic trial 
division strategy is as follows:

1) Scan through the array and 'mark' locations that meet criteria 
indicating they may factor completely over the factor base.  

2) 'Filter' the marked locations by trial dividing by small primes
that we did not sieve.  These primes are all less than 256.  If after
removing small primes the location does not meet another set of criteria,
remove it from the 'marked' list (do not subject it to further trial
division).

3) Divide out primes from the factor base between 256 and 2^13 or 2^14, 
depending on the version (2^13 for 32k version, 2^14 for 64k).  

4) Resieve primes between 2^{13|14} and 2^16, max.  

5) Primes larger than 2^16 will have been bucket sieved.  Remove these
by scanning the buckets for sieve hits equal to the current block location.

6) If applicable/appropriate, factor a remaining composite with squfof

this file contains code implementing 1)


*/

int check_relations_siqs_4(uint32 blocknum, uint8 parity,
    static_conf_t *sconf, dynamic_conf_t *dconf)
{
    printf("check_relations_siqs_4 not implemented for KNL\n");
    exit(-1);
}

int check_relations_siqs_8(uint32 blocknum, uint8 parity,
    static_conf_t *sconf, dynamic_conf_t *dconf)
{
    printf("check_relations_siqs_8 not implemented for KNL\n");
    exit(-1);
}


int check_relations_siqs_16(uint32 blocknum, uint8 parity,
    static_conf_t *sconf, dynamic_conf_t *dconf)
{
    //unrolled x128; for large inputs
    uint32 i, j, it = sconf->qs_blocksize >> 3;
    uint32 thisloc;
    uint32 num_reports = 0;
    uint64 *sieveblock;
    uint32 *sieveblock32;

    __m512i vmask = _mm512_set1_epi32(0x80808080);

    dconf->num_reports = 0;
    sieveblock = (uint64 *)dconf->sieve;
    sieveblock32 = (uint32 *)dconf->sieve;

    for (j = 0; j < 4096; j += 8)	
    {		
        __mmask16 r_msk;
        int idx;
        int k;

        __m512i vsieve = _mm512_load_epi32((__m512i *)(&sieveblock[j]));
        r_msk = _mm512_test_epi32_mask(vsieve, vmask);

        thisloc = j * 8;
    
        while (r_msk > 0)
        {
            uint32 a_msk;
           
            idx = __builtin_ctzl(r_msk);
            
            // each lit bit identifies 4 possible bytes meeting our criteria
            // for a possible relation.
            a_msk = sieveblock32[(thisloc >> 2) + idx] & 0x80808080;

            do 
            {
                k = __builtin_ctzl(a_msk) >> 3;
                dconf->reports[dconf->num_reports++] = thisloc + k + idx*4;
                a_msk = _blsr_u32(a_msk);
            } while (a_msk > 0);

            r_msk = _blsr_u32(r_msk);
        }
      
    }

    if (dconf->num_reports >= MAX_SIEVE_REPORTS)
        dconf->num_reports = MAX_SIEVE_REPORTS-1;

    dconf->total_reports += dconf->num_reports;
    dconf->total_blocks++;

    //remove small primes, and test if its worth continuing for each report
    filter_SPV(parity, dconf->sieve, dconf->numB-1,blocknum,sconf,dconf);	
    tdiv_med_ptr(parity, dconf->numB - 1, blocknum, sconf, dconf);
    resieve_med_ptr(parity, dconf->numB - 1, blocknum, sconf, dconf);

    // factor all reports in this block
    for (j = 0; j<dconf->num_reports; j++)
    {
        if (dconf->valid_Qs[j])
        {
            dconf->total_surviving_reports++;
            tdiv_LP(j, parity, blocknum, sconf, dconf);
            trial_divide_Q_siqs(j, parity, dconf->numB - 1, blocknum, sconf, dconf);
        }
    }

 //   for (i = 0, j = 0; j < dconf->num_reports; j++)
 //   {
 //       if (dconf->valid_Qs[j])
 //       {
 //           survive_locs[i++] = j;
 //       }
 //   }
 //   resieve_medprimes_32k_knl(survive_locs, i, parity, dconf->numB - 1, blocknum, sconf, dconf);

	//// factor all reports in this block
	//for (j=0; j<dconf->num_reports; j++)
	//{
	//	if (dconf->valid_Qs[j])
	//	{
 //           //resieve_med_ptr(parity, dconf->numB - 1, blocknum, sconf, dconf);
 //           dconf->total_surviving_reports++;
	//		tdiv_LP(j, parity, blocknum, sconf, dconf);
	//		trial_divide_Q_siqs(j, parity, dconf->numB-1, blocknum,sconf,dconf);
	//	}
	//}

	return 0;
}


#endif

