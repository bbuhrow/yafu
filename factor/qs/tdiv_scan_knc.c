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
#include "ytools.h"
#include "common.h"

#if defined(TARGET_KNC)
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

#define SCAN_MASK 0x8080808080808080ULL

	//when we compress small primes into 16 bits of a 32 bit field, the
	//trick of fooling the sieve routine to not sieve those roots which
	//divide poly_a fails when the blocksize is 2^16, because we're doing this:
	//root1 = fbptr->roots & 0xFFFF;
	//root2 = fbptr->roots >> 16;
	//set_aprime_roots sets roots to all 1's, which then results in root1 and 
	//root2 being set to 65535 in the sieve routine.  this, of course, isn't right
	//so the sieve location 65535 is corrupted by many small prime hits when it
	//shouldn't be, and thus we might end up here more often then we should for
	//offset 65535.
	//
	//even if we do end up here when we shouldn't, often we'll fail to find many
	//small primes which actually divide this location, and we'll bail anyway.  this
	//is safe because we explicitly trial divide by these small primes.  
	//if we make it past the small prime test and go to check the progression
	//of a prime which divides poly_a then the roots we arrive at are false (65535 again)
	//but our computation of the progression will always be 65535 + prime - blocksize,
	//since we set the root to 65535 during the sieve step as well.  
	//65535 != 65535 + prime - blocksize, so we are safe here as well.
	//we may incur more trial division than necessary - is that better than always
	//throwing away block location 65535 - NO, empirically it is much better to just
	//always bail for location 65535 when the blocksize is 65536.

	// also need to bail on location 65534, because otherwise the 8x sse2 asm division
	// can fail.  this is because we add 1 to the block loc and then add the correction
	// factor on top of that, which can overflow if blockloc >= 65534.

	// I think that throwing away 3/1000th of 1 percent of the sieve hits in exchange
	// for the speedups associated with 8x sse2 asm division and compression of
	// small primes is worth it (on 64k builds only).

int check_relations_siqs_16(uint32_t blocknum, uint8_t parity,
    static_conf_t *sconf, dynamic_conf_t *dconf)
{
    //unrolled x128; for large inputs
    uint32_t i, j, it = sconf->qs_blocksize >> 3;
    uint32_t thisloc;
    uint32_t num_reports = 0;
    uint64_t *sieveblock;

    __m512i vmask = _mm512_set1_epi32(0x80808080);


    dconf->num_reports = 0;
    sieveblock = (uint64_t *)dconf->sieve;

    for (j = 0; j < 4096; j += 8)	
    {		
        __mmask16 r_msk;
        int idx;
        int k;

        __m512i vsieve = _mm512_load_epi32((__m512i *)(&sieveblock[j]));
        r_msk = _mm512_test_epi32_mask(vsieve, vmask);

        if (r_msk == 0)
            continue;

        thisloc = j * 8;
    
        // _mm_tzcnti_32 did not appear to work as advertised... at least,
        // the rels/sec figure of merit dropped noticably when I tried
        // it.  generally that means that something has gone wrong, when
        // otherwise valid relations are not being identified.
        while ((idx = _mm_tzcnt_32(r_msk)) < 16)
        {
            // each lit bit identifies 4 possible bytes meeting our criteria
            // for a possible relation.
            for (k = 0; k < 4; k++)
            {
                if (dconf->sieve[thisloc + k + idx*4] & 0x80)
                {
                    dconf->reports[dconf->num_reports++] = thisloc + k + idx*4;
                }
            }
            r_msk ^= (1 << idx);
        }
    }

	if (dconf->num_reports >= MAX_SIEVE_REPORTS)
		dconf->num_reports = MAX_SIEVE_REPORTS-1;

    dconf->total_reports += dconf->num_reports;
    dconf->total_blocks++;

	//remove small primes, and test if its worth continuing for each report
	filter_SPV(parity, dconf->sieve, dconf->numB-1,blocknum,sconf,dconf);
	tdiv_med_ptr(parity, dconf->numB-1,blocknum,sconf,dconf);
	resieve_med_ptr(parity, dconf->numB-1,blocknum,sconf,dconf);

	// factor all reports in this block
	for (j=0; j<dconf->num_reports; j++)
	{
		if (dconf->valid_Qs[j])
		{
            dconf->total_surviving_reports++;
			tdiv_LP(j, parity, blocknum, sconf, dconf);
			trial_divide_Q_siqs(j, parity, dconf->numB-1, blocknum,sconf,dconf);
		}
	}

	return 0;
}


#endif

