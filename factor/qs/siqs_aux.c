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

#include "qs_impl.h"
#include "soe.h"
#include "ytools.h"
#include "gmp_xface.h"

uint32_t make_fb_siqs(static_conf_t *sconf)
{
	// finds the factor base primes, and computes the solutions to the congruence x^2 = N mod p
	// for the QS, these are the starting positions of the sieve relative to the sqrt of N.
	// for the MPQS, additional work using the polynomial coefficents and these congruences 
	// needs to be done to compute the starting positions of the sieve.

	//locals
	int i;
	uint32_t b,j,r,k;
	uint32_t prime, root1, root2;
	uint8_t logp;
    uint32_t urange = sconf->factor_base->B * sconf->multiplier * 1.1; // 10000000;
	uint32_t lrange = 0;
	uint64_t f;
#ifdef USE_8X_MOD_ASM
	uint32_t shift = 24;
#endif
    soe_staticdata_t* sdata = soe_init(0, 1, 32768);


	//unpack stuff from static data structure
	fb_list *fb = sconf->factor_base;
	mpz_ptr n = sconf->n;
	uint32_t mul = sconf->multiplier;
	uint32_t *modsqrt = sconf->modsqrt_array;

    //printf("make_fb: summoning primes %u:%u with B=%u\n", lrange, urange, fb->B);
    if (siqs_maxp < urange)
    {
        free(siqs_primes);
        siqs_primes = soe_wrapper(sdata, 0, urange, 0, &siqs_nump, 0, 0);
        //printf("found %lu primes from %lu to %lu\n", NUM_P, P_MIN, P_MAX);

        //siqs_primes = soe_wrapper(siqs_primes, siqs_nump, lrange, urange, 0, &NUM_P);
        siqs_minp = siqs_primes[0];
        siqs_maxp = siqs_primes[siqs_nump - 1];
    }

	//the 0th and 1st elements in the fb are always 1 and 2, so start searching with 3
	j=2; i=1;
	while (j<fb->B)
	{
		if ((uint32_t)i >= siqs_nump)
		{
			lrange = urange + 1;
			urange = lrange + 10000000;
            if (siqs_maxp < urange)
            {
                free(siqs_primes);
                //printf("make_fb: summoning primes %u:%u with B=%u\n", lrange, urange, fb->B);
                siqs_primes = soe_wrapper(sdata, lrange, urange, 0, &siqs_nump, 0, 0);
                //siqs_primes = soe_wrapper(siqs_primes, siqs_nump, lrange, urange, 0, &NUM_P);
                siqs_minp = siqs_primes[0];
                siqs_maxp = siqs_primes[siqs_nump - 1];
                //printf("found %lu primes from %lu to %lu\n", NUM_P, P_MIN, P_MAX);
                i = 0;
            }
		}

		prime = (uint32_t)siqs_primes[i];
		r = mpz_tdiv_ui(n, prime);
		if (r == 0)
		{
			if (mul % prime != 0)
			{
				//prime doesn't divide the multiplier, so
				//this prime divides the input, divide it out and bail
				mpz_tdiv_q_ui(n, n, prime);
                soe_finalize(sdata);
				return prime;
			}

			//p divides n, which means it divides the multiplier.
			//we can still use it, but it only has one solution to x^2 == n mod p instead
			//of two.  just divide its logprime in half.
			//we also can't find the root using shanks-tonelli, but it will be very small
			//because the multiplier is very small, so just use brute force.
			b = mpz_tdiv_ui(n, prime);
			k=0;
			while (1)
			{
				if (((k*k) % prime) == b)
					break;
				k++;
			}
			root1 = k;
			root2 = prime - k;

			//compute logp
			logp = (uint8_t)(log((double)prime)/log(2.0) + .5)/2;

			//fill in factor base
			fb->list->prime[j] = prime;
			modsqrt[j] = root1;
			fb->list->logprime[j] = logp;

            fb->list->binv[j] = (1ULL << 32) / (uint64_t)prime;

			//store a couple things so we can replace single precision
			//mods with shifts and muls during trial division
			//the shift determines how big the dividends can be - it should be 
			//greater than the max possible dividend.  The shift should not be too
			//big though.  Max word size (inverse * dividend) >> shift should leave
			//more bits than the max divisor.
			if (prime < sconf->small_limit)
			{
				fb->tinylist->prime[j] = prime;
				fb->tinylist->logprime[j] = logp;

				fb->tinylist->small_inv[j] = (uint32_t)(((uint64_t)1 << 32) / (uint64_t)prime);
				if (floor(MP_RADIX / (double)prime + 0.5) ==
								(double)fb->tinylist->small_inv[j]) {
					fb->tinylist->correction[j] = 1;
				}
				else {
					fb->tinylist->correction[j] = 0;
					fb->tinylist->small_inv[j]++;
				}
			}
			else
			{
#ifdef USE_8X_MOD_ASM
				if ((shift == 24) && 
					(prime > 1024) && 
					(j % 8 == 0))
					shift = 26;

				if ((shift == 26) && 
					(prime > 4096) && 
					(j % 8 == 0))
					shift = 28;

				fb->list->small_inv[j] = (uint16_t)(((uint32_t)1 << shift) / prime);
				if (floor((double)(1 << shift) / (double)prime + 0.5) ==
								(double)fb->list->small_inv[j]) {
					fb->list->correction[j] = 1;
				}
				else {
					fb->list->correction[j] = 0;
					fb->list->small_inv[j]++;
				}
#else
				fb->list->small_inv[j] = (uint32_t)(((uint64_t)1 << FOGSHIFT) / (uint64_t)prime);
				if (floor((double)(1ULL << FOGSHIFT) / (double)prime + 0.5) ==
								(double)fb->list->small_inv[j]) {
					fb->list->correction[j] = 1;
				}
				else {
					fb->list->correction[j] = 0;
					fb->list->small_inv[j]++;
				}
#endif
			}

			j++;
			i++;
			continue;
		}

		b = jacobi_1((uint64_t)r,(uint64_t)prime);
		if (b==1)
		{
			//this prime works
			ShanksTonelli_1((uint64_t)r,(uint64_t)prime,&f);
			root1 = (uint32_t)f;
			root2 = prime - root1;

			//compute logp
			logp = (uint8_t)(log((double)prime)/log(2.0) + .5);

			//fill in factor base
			fb->list->prime[j] = prime;
			modsqrt[j] = root1;
			fb->list->logprime[j] = logp;

            fb->list->binv[j] = (1ULL << 32) / (uint64_t)prime;

			//store a couple things so we can replace single precision
			//mods with shifts and muls during trial division
			//this is very fragile... need better range checking.
			if (prime < sconf->small_limit)
			{
				fb->tinylist->prime[j] = prime;
				fb->tinylist->logprime[j] = logp;

				fb->tinylist->small_inv[j] = (uint32_t)(((uint64_t)1 << 32) / (uint64_t)prime);
				if (floor(MP_RADIX / (double)prime + 0.5) ==
								(double)fb->tinylist->small_inv[j]) {
					fb->tinylist->correction[j] = 1;
				}
				else {
					fb->tinylist->correction[j] = 0;
					fb->tinylist->small_inv[j]++;
				}
			}
			else
			{
#ifdef USE_8X_MOD_ASM

				if ((shift == 24) && 
					(prime > 1024) && 
					(j % 8 == 0))
					shift = 26;

				if ((shift == 26) && 
					(prime > 4096) && 
					(j % 8 == 0))
					shift = 28;

				fb->list->small_inv[j] = (uint16_t)(((uint32_t)1 << shift) / prime);
				if (floor((double)(1 << shift) / (double)prime + 0.5) ==
								(double)fb->list->small_inv[j]) {
					fb->list->correction[j] = 1;
				}
				else {
					fb->list->correction[j] = 0;
					fb->list->small_inv[j]++;
				}
				
#else
				fb->list->small_inv[j] = (uint32_t)(((uint64_t)1 << FOGSHIFT) / (uint64_t)prime);
				if (floor((double)(1ULL << FOGSHIFT) / (double)prime + 0.5) ==
								(double)fb->list->small_inv[j]) {
					fb->list->correction[j] = 1;
				}
				else {
					fb->list->correction[j] = 0;
					fb->list->small_inv[j]++;
				}
#endif
			}

			
			j++;
		}
		i++;
	}

    soe_finalize(sdata);
	return 0;
}

#define NUM_PARAM_ROWS 31
void get_params(static_conf_t *sconf)
{
	int bits,i;
	double scale;
	fb_list *fb = sconf->factor_base;

	//parameter table
	//bits, fb primes, lp mulitplier, 32k blocks, dlp_exp, tlp_exp, tlp_bdiv
	//adjustment in v1.27 - more primes and less blocks for numbers > ~80 digits
	//also different scaling for numbers bigger than 100 digits (constant increase
	//of 20% per line)
#if defined(SMALL_SIQS_INTERVALS)

    // does much better with much smaller sieve intervals.  maybe 
    // because smaller L2?
    int param_table[NUM_PARAM_ROWS][7] = {
        {50,	30,	30,	1, 1.8, 0, 0},
        { 60, 36, 40, 1  , 1.8, 0, 0},
        { 70, 50, 40, 1  , 1.8, 0, 0},
        { 80, 80, 40, 1  , 1.8, 0, 0},
        { 90, 120, 40, 1 , 1.8, 0, 0},
        { 100, 175, 50, 1, 1.8, 0, 0},
        { 110, 275, 50, 1, 1.8, 0, 0},
        { 120, 375, 50, 1, 1.8, 0, 0},

        { 140, 828, 50, 1     , 1.8, 0, 0},
        { 149, 1028, 60, 1    , 1.8, 0, 0},
        { 165, 1228, 60, 1    , 1.8, 0, 0},
        { 181, 2247, 70, 1    , 1.8, 0, 0},
        { 198, 3485, 70, 1    , 1.8, 0, 0},
        { 215, 6357, 80, 1    , 1.8, 0, 0},
        { 232, 12132, 80, 1   , 1.8, 0, 0},      // 70 digits
        { 248, 26379, 90, 1   , 1.8, 0, 0},
        { 265, 47158, 90, 1   , 1.8, 0, 0},      // 80 digits
        { 281, 60650, 100, 2  , 1.8, 0, 0},
        { 298, 71768, 120, 3  , 1.8, 0, 0},     // 90 digits
        { 310, 86071, 120, 3  , 1.8, 0, 0},
        { 320, 99745, 140, 4  , 1.8, 0, 0},
        { 330, 115500, 150, 4 , 1.8, 0, 0},     // 100 digits
        { 340, 138600, 150, 5 , 1.8, 0, 0},
        { 350, 166320, 150, 5 , 1.8, 0, 0},    // 105 digits
        { 360, 199584, 150, 6 , 1.8, 0, 0},
        { 370, 239500, 150, 6 , 1.8, 0, 0},    // 110 digits
        { 380, 287400, 175, 7 , 1.8, 0, 0},
        { 390, 344881, 175, 7 , 1.8, 0, 0},
        { 400, 413857, 175, 8 , 1.8, 0, 0},
        { 410, 496628, 175, 8 , 1.8, 0, 0},
        { 462, 750000, 175, 12, 1.8, 0, 0},        // 140 digits
};

#else

    // Haswell also likes smaller sieve intervals for the 
    // larger numbers... maybe this table needs to be 
    // processor dependent...
#if defined(USE_AVX512F)
    // another adjustment smaller (testing on skylake-x)
    double param_table[NUM_PARAM_ROWS][7] = {
        {50,	30,	    30,	    1,  1.0, 0, 0},
        {60,	36,	    40,	    1,  1.0, 0, 0},
        {70,	50,	    40,	    1,  1.0, 0, 0},
        {80,	80,	    40,	    1,  1.0, 0, 0},
        {90,	120,	40,	    1,  1.0, 0, 0},
        {100,	175,	50,	    2,  1.0, 0, 0},
        {110,	275,	50,	    2,  1.0, 0, 0},
        {120,	375,	50,	    2,  1.0, 0, 0},
        {140,	828,	50,	    2,  1.0, 0, 0},
        {149,	1028,	60,	    2,  1.0, 0, 0},
        {165,	1228,	60,	    2,  1.0, 0, 0},
        {181,	2247,	70,	    2,  1.0, 0, 0},
        {198,	3485,	70,	    2,  1.0, 0, 0},
        {215,	6357,	80,	    2,  1.0, 0, 0},         
        {232,	12132,	80,	    4,  1.75, 0, 0},         // 70 digits, DLP starts
        {248,	26379,	90,	    6,  1.75, 0, 0},         // 75 digits
        {265,	47158,	90,	    6,  1.75, 0, 0},         // 80 digits
        {281,	60650,	100,	8,  1.8, 0, 0},
        {298,	71768,	120,	8,  1.8, 0, 0},         // 90 digits
        {310,	86071 ,	120,	10, 1.85, 0, 0},
        {320,	99745 ,	140,	10, 1.85, 0, 0},        // 95 digits
        {330,	115500, 150,    12, 1.85, 0, 0},        // 100 digits 

        /* These are the DLP numbers (guesses) */
        {340,	139120, 150,    12, 1.85, 0, 2},
        {350,	166320, 150,    14, 1.85, 0, 2},        // 105 digits
        {360,   199584, 150,    14, 1.85, 0, 2},
        {370,	239500, 150,    16, 1.9, 0, 2},        // 110 digits
        {380,	287400, 175,    16, 1.9, 0, 2},
        {390,	344881, 175,    18, 1.9, 0, 2},
        {400,	413857, 175,    20, 1.9, 0, 2},        // 120 digits
        {410,	496628, 175,    22, 1.9, 0, 2},
        

        /* These are the TLP numbers (guesses)
        // row 22
        {340,	72000, 27,    6  , 1.8, 2.8, 2},          // here we start using TLP: 3rd column is LPB in bits
        {350,	80880,  27,    6 , 1.8, 2.8, 2},         // 105 digits
        {360,   139708, 28,    7 , 1.8, 2.8, 2},
        {370,	167650, 28,    8 , 1.8, 2.8, 2},        // 110 digits
        {380,	201180, 29,    8 , 1.8, 2.8, 2},
        {390,	241416, 29,    10, 1.8, 2.8, 2},
        {396,	289700, 30,    12, 1.8, 2.8, 2},        // 120 digits
        {429,	500000, 31,    14, 1.8, 2.8, 2},        // 130 digits
        {462,	750000, 32,    16, 1.8, 2.8, 2},        // 140 digits
        */

        
    }; 
#else
    int param_table[NUM_PARAM_ROWS][7] = {
        {50,	30,	    30,	    2 , 1.8, 0, 0},
        {60,	36,	    40,	    2 , 1.8, 0, 0},
        {70,	50,	    40,	    2 , 1.8, 0, 0},
        {80,	80,	    40,	    2 , 1.8, 0, 0},
        {90,	120,	40,	    2 , 1.8, 0, 0},
        {100,	175,	50,	    2 , 1.8, 0, 0},
        {110,	275,	50,	    2 , 1.8, 0, 0},
        {120,	375,	50,	    2 , 1.8, 0, 0},
        {140,	828,	50,	    2 , 1.8, 0, 0},
        {149,	1028,	60,	    2 , 1.8, 0, 0},
        {165,	1228,	60,	    2 , 1.8, 0, 0},
        {181,	2247,	70,	    2 , 1.8, 0, 0},
        {198,	3485,	70,	    4 , 1.8, 0, 0},
        {215,	6357,	80,	    4 , 1.8, 0, 0},
        {232,	12132,	80,	    6 , 1.8, 0, 0},      // 70 digits
        {248,	26379,	90,	    8 , 1.8, 0, 0},
        {265,	47158,	90,	    10, 1.8, 0, 0},     // 80 digits
        {281,	60650,	100,	12, 1.8, 0, 0},
        {298,	71768,	120,	12, 1.8, 0, 0},     // 90 digits
        {310,	86071 ,	120,	14, 1.8, 0, 0},
        {320,	99745 ,	140,	16, 1.8, 0, 0},     // 95 digits
        {330,	115500, 150,    16, 1.8, 0, 0},     // 100 digits 

        /* These are the DLP numbers (guesses) */
        {340,	139120, 150,    12, 1.85, 0, 2},
        {350,	166320, 150,    14, 1.85, 0, 2},        // 105 digits
        {360,   199584, 150,    14, 1.85, 0, 2},
        {370,	239500, 150,    16, 1.9, 0, 2},        // 110 digits
        {380,	287400, 175,    16, 1.9, 0, 2},
        {390,	344881, 175,    18, 1.9, 0, 2},
        {400,	413857, 175,    20, 1.9, 0, 2},        // 120 digits
        {410,	496628, 175,    22, 1.9, 0, 2},


        /* These are the TLP numbers (guesses)

        // row 22
        {340,	72000, 27,    6  , 1.9, 2.6, 3},          // here we start using TLP: 3rd column is LPB in bits
        {350,	80880,  27,    6 , 1.9, 2.6, 3},         // 105 digits
        {360,   139708, 28,    7 , 1.9, 2.7, 3},
        {370,	167650, 28,    8 , 1.9, 2.7, 3},        // 110 digits
        {380,	201180, 29,    8 , 1.9, 2.8, 2.5},
        {390,	241416, 29,    10, 1.9, 2.8, 2.5},
        {396,	289700, 30,    12, 1.95, 2.9, 2},        // 120 digits
        {429,	500000, 31,    14, 1.95, 2.9, 2},        // 130 digits
        {462,	750000, 32,    16, 1.95, 2.95, 2},        // 140 digits


        // These were the DLP numbers (guesses)

        {340,	138600, 150,    18},     // above 100 it is basically guesswork
        {350,	166320, 150,    18},     // 105 digits
        {360,   199584, 150,    20},
        {370,	239500, 150,    22},     // 110 digits
        {380,	287400, 175,    24},
        {390,	344881, 175,    26},
        {400,	413857, 175,    28},
        {410,	496628, 175,    30},
        */

        
    };
#endif
#endif

    int param_table_bkup[NUM_PARAM_ROWS][4] = {
        { 50, 30, 30, 1 },
        { 60, 36, 40, 1 },
        { 70, 50, 40, 1 },
        { 80, 80, 40, 1 },
        { 90, 120, 40, 1 },
        { 100, 175, 50, 1 },
        { 110, 275, 50, 1 },
        { 120, 375, 50, 1 },

        { 140, 828, 50, 1 },
        { 149, 1028, 50, 1 },
        { 165, 1228, 50, 1 },
        { 181, 2247, 50, 1 },
        { 198, 3485, 60, 2 },
        { 215, 6357, 60, 2 },
        { 232, 12132, 70, 3 },      // 70 digits
        { 248, 26379, 80, 4 },
        { 265, 47158, 90, 5 },      // 80 digits
        { 281, 60650, 100, 6 },
        { 298, 71768, 120, 7 },     // 90 digits
        { 310, 86071, 120, 8 },     
        { 320, 99745, 140, 9 },     // 95 digits (for benchmark c95)
        { 330, 115500, 150, 10 },   // 100 digits
        { 340, 138600, 150, 12 },
        { 350, 166320, 150, 14 },
        { 360, 199584, 150, 16 },   // 110 digits
        { 370, 239500, 150, 18 },
        { 380, 287400, 175, 22 },
        { 390, 344881, 175, 26 },
        { 400, 413857, 175, 30 },
        { 410, 496628, 175, 32 },
    };

	/*
	int param_table[22][4] = {
		{140,	600,	40,	    1},
		{149,	875,	40,	    1},
		{165,	1228,	50,	    1},
		{181,	2247,	50,	    1},
		{198,	3485,	60,	    2},
		{215,	6357,	60,	    2},	
		{232,	12132,	70,	    3},
		{248,	26379,	80,	    4},
		{265,	42871,	90,	    6},
		{281,	55137,	100,	8},
		{298,	65244,	120,	10},
		{310,	78247,	120,	12},
		{320,	90678,	140,	14},
		{330,	105000, 150,    18},
		{340,	125000, 150,    21},
		{350,	155000, 150,    25},
		{360,	195000, 150,    29},
		{370,	250000, 150,    34},
		{380,	310000, 150,    40},
		{390,	380000, 150,    47},
		{400,	460000, 150,    55},
		{410,	550000, 150,    64},
	};
	*/

	// linear interpolation according to bit size to determine
	// factor base bound.  use the closest parameter for lp multiplier
	// and number of blocks.

	bits = sconf->obj->bits;

	fb->B = 0;
	if (bits <= param_table[0][0])
	{
		scale = (double)bits / (double)param_table[0][0];
		fb->B = (uint32_t)(scale * (double)(param_table[0][1]));		
		sconf->large_mult = 40;
		sconf->num_blocks = 1;
	}
	else if (bits <= 1000)       // up to 100 digits, above that we use TLP
	{
        for (i = 0; i < NUM_PARAM_ROWS - 1; i++)
        {
            if (bits > param_table[i][0] && bits <= param_table[i + 1][0])
            {
                scale = (double)(param_table[i + 1][0] - bits) /
                    (double)(param_table[i + 1][0] - param_table[i][0]);
                fb->B = param_table[i + 1][1] -
                    (uint32_t)(scale * (double)(param_table[i + 1][1] - param_table[i][1]));

                sconf->large_mult = (uint32_t)((param_table[i + 1][2] + param_table[i][2]) / 2.0 + 0.5);
                sconf->num_blocks = (uint32_t)((param_table[i + 1][3] + param_table[i][3]) / 2.0 + 0.5);
            }
        }
	}
    else
    {
        // this can still be overridden with forceDLP.
        sconf->obj->qs_obj.gbl_force_TLP = 1;
        for (i = 21; i < NUM_PARAM_ROWS - 1; i++)
        {
            if ((bits > param_table[i][0]) && (bits <= param_table[i + 1][0]))
            {
                scale = (double)(param_table[i + 1][0] - bits) /
                    (double)(param_table[i + 1][0] - param_table[i][0]);
                fb->B = param_table[i + 1][1] -
                    (uint32_t)(scale * (double)(param_table[i + 1][1] - param_table[i][1]));

                if ((param_table[i + 1][0] - bits) < (bits - param_table[i][0]))
                {
                    sconf->large_mult = 1 << (int)param_table[i + 1][2];
                    sconf->num_blocks = param_table[i + 1][3];
                }
                else
                {
                    sconf->large_mult = 1 << (int)param_table[i][2];
                    sconf->num_blocks = param_table[i][3];
                }
            }
        }
    }

	if (fb->B == 0)
	{
		//off the end of the table, extrapolate based on the slope of 
		//the last two

		scale = (double)(param_table[NUM_PARAM_ROWS-1][1] - param_table[NUM_PARAM_ROWS-2][1]) /
			(double)(param_table[NUM_PARAM_ROWS-1][0] - param_table[NUM_PARAM_ROWS-2][0]);
		fb->B = (uint32_t)(((double)bits - param_table[NUM_PARAM_ROWS-1][0]) * 
			scale + param_table[NUM_PARAM_ROWS-1][1]);
		sconf->large_mult = param_table[NUM_PARAM_ROWS-1][2];	//reuse last one

		scale = (double)(param_table[NUM_PARAM_ROWS-1][3] - param_table[NUM_PARAM_ROWS-2][3]) /
			(double)(param_table[NUM_PARAM_ROWS-1][0] - param_table[NUM_PARAM_ROWS-2][0]);
		//sconf->num_blocks = param_table[NUM_PARAM_ROWS-1][3];	//reuse last one
		sconf->num_blocks = (uint32_t)(((double)bits - param_table[NUM_PARAM_ROWS-1][0]) * 
			scale + param_table[NUM_PARAM_ROWS-1][3]);

	}

    // make B divisible by 16
    while ((fb->B & 15) != 0)
    {
        fb->B++;
    }

	// minimum factor base - for use with really small inputs.
	// not efficient, but needed for decent poly selection
	//if (fb->B < 250)
    //		fb->B = 250;

	if (sconf->obj->qs_obj.gbl_override_B_flag)
		fb->B = sconf->obj->qs_obj.gbl_override_B;

	if (sconf->obj->qs_obj.gbl_override_blocks_flag)
		sconf->num_blocks = sconf->obj->qs_obj.gbl_override_blocks;

	if (sconf->obj->qs_obj.gbl_override_lpmult_flag)
		sconf->large_mult = sconf->obj->qs_obj.gbl_override_lpmult;

	return;
}

int qcomp_siqs(const void *x, const void *y)
{
	siqs_r **xx = (siqs_r **)x;
	siqs_r **yy = (siqs_r **)y;
	
	if (xx[0]->large_prime[0] > yy[0]->large_prime[0])
		return 1;
	else if (xx[0]->large_prime[0] == yy[0]->large_prime[0])
		return 0;
	else
		return -1;
}

void set_aprime_roots(static_conf_t *sconf, uint32_t val, int *qli, int s, 
	sieve_fb_compressed *fb, int action)
{
	int i;
	fb_list *fullfb = sconf->factor_base;

    /* invalid roots are currently marked by being set to 65535.  this works when we
    explicitly check the root against the blocksize before sieving it, but here where
    we've completely unrolled the loop it doesn't work.  what we could do is set roots
    and primes == 0 when roots are invalid instead of setting roots to 65535 for the
    range of primes that are going to be treated this way.  then only the first location
    in every block gets hosed and we can tell tdiv to always ignore that location.  the
    snippet of code below sets roots and primes = 0 for invalid roots */


    if (action == 1)
    {
        for (i = 0; i<s; i++)
        {
            if ((fullfb->list->prime[qli[i]] > 8192)) // && (fullfb->list->prime[qli[i]] < sconf->qs_blocksize))
            {
                fb->root1[qli[i]] = 0;
                fb->prime[qli[i]] = 0;
                fb->root2[qli[i]] = 0;
            }
            else
            {
                fb->root1[qli[i]] = 0xffff;
                fb->root2[qli[i]] = 0xffff;
            }
        }
    }
    else
    {
        for (i = 0; i < s; i++)
        {
            fb->root1[qli[i]] = 0xffff;
            fb->prime[qli[i]] = fullfb->list->prime[qli[i]];
            fb->root2[qli[i]] = 0xffff;
        }
    }

    /*
	for (i=0;i<s;i++)
	{		
		
		if ((fullfb->list->prime[qli[i]] > 8192)) // && (fullfb->list->prime[qli[i]] < sconf->qs_blocksize))
		{
			if (action == 1)
			{
				//printf("zeroing roots and primes at index %d, prime = %u\n", qli[i], fullfb->list->prime[qli[i]]);
				fb->root1[qli[i]] = 0;
				fb->prime[qli[i]] = 0;
				fb->root2[qli[i]] = 0;
			}
			else
			{
				//printf("restoring roots and primes at index %d, prime = %u\n", qli[i], fullfb->list->prime[qli[i]]);
				fb->root1[qli[i]] = 0xffff;
				fb->prime[qli[i]] = fullfb->list->prime[qli[i]];
				fb->root2[qli[i]] = 0xffff;
			}
		}
		else
		{
			fb->root1[qli[i]] = 0xffff;
			fb->root2[qli[i]] = 0xffff;
		}
	}
    */

	return;
}

void get_gray_code(siqs_poly *poly) 
{

	int i, v, j, n = poly->s;
	int tmp;

	for (i = 1; i < (1 << (n - 1)); i++) {
		v = 1;
		j = i;
		while ((j & 1) == 0)
			v++, j >>= 1;
		tmp = i + (1 << v) - 1;
		tmp = (tmp >> v);
		poly->nu[i] = v;
		if (tmp & 1)
			poly->gray[i] = -1;
		else
			poly->gray[i] = 1;
	}
	return;
}

uint32_t yafu_factor_list_add(fact_obj_t *obj, factor_list_t *list, 
				mpz_t new_factor) {

	uint32_t i, bitsleft;
	int isnew = 1;
	mpz_t tmpz;

	mpz_init(tmpz);

	//look to see if we've already included this one
	for (i=0; i<list->num_factors; i++)
	{
		mp_t2gmp(&list->final_factors[i]->factor, tmpz);
		isnew &= (mpz_cmp(tmpz,new_factor) != 0);
	}

	if (isnew)
	{
        if (obj->logfile != NULL)
        {
            char* s = mpz_get_str(NULL, 10, new_factor);
            logprint(obj->logfile,
                "prp%d = %s\n", gmp_base10(new_factor), s);
            free(s);
        }

		list->final_factors[list->num_factors] = (final_factor_t *)malloc(
			sizeof(final_factor_t));
		gmp2mp_t(new_factor, &list->final_factors[list->num_factors]->factor);
		list->num_factors++;
	}

	//now determine if we are done based on the bits of factors found compared to 
	//the bits in the original n
	bitsleft = obj->bits;
	for (i=0; i<list->num_factors; i++)
	{
		mp_t2gmp(&list->final_factors[i]->factor, tmpz);
		bitsleft -= mpz_sizeinbase(tmpz, 2);
	}

	mpz_clear(tmpz);
	return bitsleft;
}

