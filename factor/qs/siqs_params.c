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

// ---------------------------------------------------------------------------
// built-in SIQS parameter tables, keyed by large-prime kernel ISA, [large_B, B)
// strategy and kernel (see siqs_builtin_table_t in qs_impl.h).  New tables
// (per-ISA, per-kernel) go in siqs_param_tables[] below; siqs_find_param_table
// picks the most specific match, so a table keyed on more fields wins over a
// more general one.  The last entry matches everything and must stay last.
//
// columns: bits, B (fb primes), M (lp multiplier), NB (32k blocks),
// TF (tdiv bound), TFSm (small prime variation cutoff), MFBD, MFBT.
// 0 in TF/TFSm/MFBD/MFBT = use the formula/default in siqs_static_init.
//
// adjustment in v1.27 - more primes and less blocks for numbers > ~80 digits
// also different scaling for numbers bigger than 100 digits (constant increase
// of 20% per line).  Rows up to 148 bits are now handled by kmpqs.
//
// MFBD in the bucket tables is 1.85, the value these tables actually ran with:
// the yafu driver always passed its -siqsMFBD default of 1.85 as an override,
// so the dlp_exp column in the old tables (1.75-1.9) was never used.
// ---------------------------------------------------------------------------

// lp bucket algorithm, AVX512 builds (original tuning, adjusted smaller
// testing on skylake-x).  Haswell also likes smaller sieve intervals for the
// larger numbers.
static const siqs_builtin_row_t siqs_table_bucket_avx512[] = {
    {50,	30,	    30,	    1,  0, 0, 0,    0},
    {60,	36,	    40,	    1,  0, 0, 0,    0},
    {70,	50,	    40,	    1,  0, 0, 0,    0},
    {80,	80,	    40,	    1,  0, 0, 0,    0},
    {90,	120,	40,	    1,  0, 0, 0,    0},
    {100,	175,	50,	    2,  0, 0, 0,    0},
    {110,	275,	50,	    2,  0, 0, 0,    0},
    {120,	375,	50,	    2,  0, 0, 0,    0},
    {140,	828,	50,	    2,  0, 0, 0,    0},
    {148,	1028,	60,	    2,  0, 0, 0,    0},
    {165,	1228,	60,	    2,  0, 0, 0,    0},
    {181,	2247,	70,	    2,  0, 0, 0,    0},
    {198,	3485,	70,	    2,  0, 0, 0,    0},
    {215,	6357,	80,	    2,  0, 0, 1.85, 0},     // (MFBD: for interpolation; SLP here)
    {232,	12132,	80,	    4,  0, 0, 1.85, 0},     // 70 digits, DLP starts
    {248,	26379,	90,	    6,  0, 0, 1.85, 0},     // 75 digits
    {265,	47158,	90,	    6,  0, 0, 1.85, 0},     // 80 digits
    {281,	60650,	100,	8,  0, 0, 1.85, 0},
    {298,	71768,	120,	8,  0, 0, 1.85, 0},     // 90 digits
    {310,	86071 ,	120,	10, 0, 0, 1.85, 0},
    {320,	99745 ,	140,	10, 0, 0, 1.85, 0},     // 95 digits
    {330,	115500, 150,    12, 0, 0, 1.85, 0},     // 100 digits
    // above 100 digits: guesses
    {340,	139120, 150,    12, 0, 0, 1.85, 0},
    {350,	166320, 150,    14, 0, 0, 1.85, 0},     // 105 digits
    {360,   199584, 150,    14, 0, 0, 1.85, 0},
    {370,	239500, 150,    16, 0, 0, 1.85, 0},     // 110 digits
    {380,	287400, 175,    16, 0, 0, 1.85, 0},
    {390,	344881, 175,    18, 0, 0, 1.85, 0},
    {400,	413857, 175,    20, 0, 0, 1.85, 0},     // 120 digits
    {410,	496628, 175,    22, 0, 0, 1.85, 0},
};

// lp bucket algorithm, everything else (original tuning)
static const siqs_builtin_row_t siqs_table_bucket[] = {
    {50,	30,	    30,	    2 , 0, 0, 0,    0},
    {60,	36,	    40,	    2 , 0, 0, 0,    0},
    {70,	50,	    40,	    2 , 0, 0, 0,    0},
    {80,	80,	    40,	    2 , 0, 0, 0,    0},
    {90,	120,	40,	    2 , 0, 0, 0,    0},
    {100,	175,	50,	    2 , 0, 0, 0,    0},
    {110,	275,	50,	    2 , 0, 0, 0,    0},
    {120,	375,	50,	    2 , 0, 0, 0,    0},
    {140,	828,	50,	    2 , 0, 0, 0,    0},
    {148,	1028,	60,	    2 , 0, 0, 0,    0},
    {165,	1228,	60,	    2 , 0, 0, 0,    0},
    {181,	2247,	70,	    2 , 0, 0, 0,    0},
    {198,	3485,	70,	    4 , 0, 0, 0,    0},
    {215,	6357,	80,	    4 , 0, 0, 1.85, 0},     // (MFBD: for interpolation; SLP here)
    {232,	12132,	80,	    6 , 0, 0, 1.85, 0},     // 70 digits
    {248,	26379,	90,	    8 , 0, 0, 1.85, 0},
    {265,	47158,	90,	    10, 0, 0, 1.85, 0},     // 80 digits
    {281,	60650,	100,	12, 0, 0, 1.85, 0},
    {298,	71768,	120,	12, 0, 0, 1.85, 0},     // 90 digits
    {310,	86071 ,	120,	14, 0, 0, 1.85, 0},
    {320,	99745 ,	140,	16, 0, 0, 1.85, 0},     // 95 digits
    {330,	115500, 150,    16, 0, 0, 1.85, 0},     // 100 digits
    // above 100 digits: guesses
    {340,	139120, 150,    12, 0, 0, 1.85, 0},
    {350,	166320, 150,    14, 0, 0, 1.85, 0},     // 105 digits
    {360,   199584, 150,    14, 0, 0, 1.85, 0},
    {370,	239500, 150,    16, 0, 0, 1.85, 0},     // 110 digits
    {380,	287400, 175,    16, 0, 0, 1.85, 0},
    {390,	344881, 175,    18, 0, 0, 1.85, 0},
    {400,	413857, 175,    20, 0, 0, 1.85, 0},     // 120 digits
    {410,	496628, 175,    22, 0, 0, 1.85, 0},

    /* These are the TLP numbers (guesses).  M is the LPB in bits.
    {340,	72000,  27,    6 , 0, 0, 1.9,  2.6},
    {350,	80880,  27,    6 , 0, 0, 1.9,  2.6},    // 105 digits
    {360,   139708, 28,    7 , 0, 0, 1.9,  2.7},
    {370,	167650, 28,    8 , 0, 0, 1.9,  2.7},    // 110 digits
    {380,	201180, 29,    8 , 0, 0, 1.9,  2.8},
    {390,	241416, 29,    10, 0, 0, 1.9,  2.8},
    {396,	289700, 30,    12, 0, 0, 1.95, 2.9},    // 120 digits
    {429,	500000, 31,    14, 0, 0, 1.95, 2.9},    // 130 digits
    {462,	750000, 32,    16, 0, 0, 1.95, 2.95},   // 140 digits
    */
};

// lp direct sieve: new tuning (smaller intervals, smaller B).  Until we have
// better tuning, for DLP we generally reduce B, increase M, decrease TF
// threshold, decrease NB.  For SLP it's about the same.
static const siqs_builtin_row_t siqs_table_direct[] = {
    {50,	30,	    30,	    1,  0, 0, 0,    0},
    {60,	36,	    40,	    1,  0, 0, 0,    0},
    {70,	50,	    40,	    1,  0, 0, 0,    0},
    {80,	80,	    40,	    1,  0, 0, 0,    0},
    {90,	120,	40,	    1,  0, 0, 0,    0},
    {100,	175,	50,	    2,  0, 0, 0,    0},
    {110,	275,	50,	    2,  0, 0, 0,    0},
    {120,	375,	50,	    2,  0, 0, 0,    0},
    {140,	828,	50,	    2,  0, 0, 0,    0},
    {148,	1028,	60,	    2,  0, 0, 0,    0},
    {165,	1616,	66,	    2,  0, 0, 0,    0},     // 50 digits
    {181,	3008,	85,	    2,  0, 0, 0,    0},     // 55 digits
    {198,	4784,	105,    2,  0, 0, 0,    0},     // 60 digits
    {215,	6592,	250,	4,  0, 0, 1.72, 0},     // 65 digits (MFBD: for interpolation; SLP here)
    {232,	8064,	240,    2,  0, 0, 1.72, 0},     // 70 digits, DLP starts
    {248,	14928,	240,	3,  0, 0, 1.69, 0},     // 75 digits
    {265,	23376,	250,	3,  0, 0, 1.72, 0},     // 80 digits
    {281,	43008,	250,	4,  0, 0, 1.72, 0},     // 85 digits
    {298,	57000,	250,	6,  0, 0, 1.8,  0},     // 90 digits
    {310,	86071 ,	120,	10, 0, 0, 1.85, 0},
    {320,	99745 ,	140,	10, 0, 0, 1.85, 0},     // 95 digits
    {330,	115500, 150,    12, 0, 0, 1.85, 0},     // 100 digits
    // above 100 digits: guesses
    {340,	139120, 150,    12, 0, 0, 1.85, 0},
    {350,	166320, 150,    14, 0, 0, 1.85, 0},     // 105 digits
    {360,   199584, 150,    14, 0, 0, 1.85, 0},
    {370,	239500, 150,    16, 0, 0, 1.9,  0},     // 110 digits
    {380,	287400, 175,    16, 0, 0, 1.9,  0},
    {390,	344881, 175,    18, 0, 0, 1.9,  0},
    {400,	413857, 175,    20, 0, 0, 1.9,  0},     // 120 digits
    {410,	496628, 175,    22, 0, 0, 1.9,  0},
};

// lp direct sieve with AVX512 kernels (mid = sparse, kernel = cmem when
// measured).  Placeholder: c60-c80 are siqs_tune_v2.py pass-1 fitted optima
// (one input per size, i5-11600K, 1 thread, MFBD fixed at 1.72, 2026-10-25),
// about 0.72-0.90x the time of the defaults at the time (siqs_table_direct
// with the formula TF and -noopt TFSm).  B and NB optima were at the edges of
// the grid (B = 1.0x the default, NB at its smallest level at <= c70, TF at
// -8 at >= c70) so the true optima may lie beyond; M is nearly flat (+-0.3%).
// Other sizes are siqs_table_direct's.
static const siqs_builtin_row_t siqs_table_direct_avx512[] = {
    {50,	30,	    30,	    1,  0,  0,  0,    0},
    {60,	36,	    40,	    1,  0,  0,  0,    0},
    {70,	50,	    40,	    1,  0,  0,  0,    0},
    {80,	80,	    40,	    1,  0,  0,  0,    0},
    {90,	120,	40,	    1,  0,  0,  0,    0},
    {100,	175,	50,	    2,  0,  0,  0,    0},
    {110,	275,	50,	    2,  0,  0,  0,    0},
    {120,	375,	50,	    2,  0,  0,  0,    0},
    {140,	828,	50,	    2,  0,  0,  0,    0},
    {148,	1028,	60,	    2,  0,  0,  0,    0},
    {165,	1616,	66,	    2,  0,  0,  0,    0},     // 50 digits
    {181,	3008,	85,	    2,  0,  0,  0,    0},     // 55 digits
    {198,	4912,	120,    1,  64, 19, 0,    0},     // 60 digits, pass-1 fit 0.851
    {215,	6512,	300,	1,  69, 18, 1.70, 0},     // 65 digits, pass-1 fit 0.896 (MFBD: for interpolation; SLP here)
    {232,	9920,	200,    1,  69, 19, 1.70, 0},     // 70 digits, pass-1 fit 0.721, DLP starts
    //{248,	14944,	190,	2,  74, 22, 1.72, 0},     // 75 digits, pass-1 fit 0.810
    //best 0.964 mfbd 1.70 (B=1,M=225,TF=+0,TFSm=+2,NB=2)
    {248,	11208,	225,	2,  74, 22, 1.70, 0},     // 75 digits, pass-1 fit 0.810
    {265,	24624,	274,	2,  79, 20, 1.72, 0},     // 80 digits, pass-1 fit 0.813
    //{281,	43008,	250,	4,  82, 20, 1.72, 0},     // 85 digits (orig)
    //{281,	32256,	325,	3,  85, 20, 1.72, 0},     // 85 digits (testing)
    //best 0.945 (B=1,M=325,TF=+3,TFSm=+2,NB=3) (potential new with 1.72)
    //best 0.946 (B=1,M=325,TF=+3,TFSm=+2,NB=3) (potential new with 1.70)
    //best 0.946 (B=1,M=225,TF=+0,TFSm=+2,NB=3) (potential new with 1.74)
    {281,	32256,	225,	3,  85, 22, 1.74, 0},     // 85 digits (accepted after testing)
    {298,	43000,	250,	4,  89, 20, 1.74, 0},     // 90 digits
    {310,	56071,	250,	4,  94, 20, 1.74, 0},
    {320,	71745,	250,	5,  97, 20, 1.74, 0},     // 95 digits
    {330,	87500,  250,    6, 100, 20, 1.74, 0},     // 100 digits
    // above 100 digits: guesses
    {340,	139120, 150,    12, 0,  0,  1.85, 0},
    {350,	166320, 150,    14, 0,  0,  1.85, 0},     // 105 digits
    {360,   199584, 150,    14, 0,  0,  1.85, 0},
    {370,	239500, 150,    16, 0,  0,  1.9,  0},     // 110 digits
    {380,	287400, 175,    16, 0,  0,  1.9,  0},
    {390,	344881, 175,    18, 0,  0,  1.9,  0},
    {400,	413857, 175,    20, 0,  0,  1.9,  0},     // 120 digits
    {410,	496628, 175,    22, 0,  0,  1.9,  0},
};

#define SIQS_TABLE(name, isa, large, kernel, rows) \
    { name, isa, large, kernel, (int)(sizeof(rows) / sizeof(rows[0])), rows }

static const siqs_builtin_table_t siqs_param_tables[] = {
    SIQS_TABLE("direct-avx512", SIQS_ISA_AVX512, SIQS_LARGE_DIRECT, SIQS_STRAT_AUTO,
        siqs_table_direct_avx512),
    SIQS_TABLE("direct", SIQS_STRAT_AUTO, SIQS_LARGE_DIRECT, SIQS_STRAT_AUTO,
        siqs_table_direct),
    SIQS_TABLE("bucket-avx512", SIQS_ISA_AVX512, SIQS_LARGE_BUCKET, SIQS_STRAT_AUTO,
        siqs_table_bucket_avx512),
    // fallback: matches everything, keep last
    SIQS_TABLE("bucket", SIQS_STRAT_AUTO, SIQS_STRAT_AUTO, SIQS_STRAT_AUTO,
        siqs_table_bucket),
};

#define NUM_PARAM_TABLES (int)(sizeof(siqs_param_tables) / sizeof(siqs_param_tables[0]))

// the most specific table whose keys all match (each key either
// SIQS_STRAT_AUTO or equal); the first such in siqs_param_tables[] on ties
const siqs_builtin_table_t* siqs_find_param_table(uint32_t isa, int lp_direct,
    uint32_t kernel)
{
    uint32_t large = lp_direct ? SIQS_LARGE_DIRECT : SIQS_LARGE_BUCKET;
    const siqs_builtin_table_t* best = &siqs_param_tables[NUM_PARAM_TABLES - 1];
    int best_score = -1;
    int i;

    for (i = 0; i < NUM_PARAM_TABLES; i++)
    {
        const siqs_builtin_table_t* t = &siqs_param_tables[i];
        int score = 0;

        if (t->isa != SIQS_STRAT_AUTO)
        {
            if (t->isa != isa) continue;
            score++;
        }
        if (t->large != SIQS_STRAT_AUTO)
        {
            if (t->large != large) continue;
            score++;
        }
        if (t->kernel != SIQS_STRAT_AUTO)
        {
            if (t->kernel != kernel) continue;
            score++;
        }
        if (score > best_score)
        {
            best = t;
            best_score = score;
        }
    }
    return best;
}

// linear interpolation of an optional (0 = not tabulated) table value:
// 0 unless both bracketing rows have it
static double siqs_interp_opt(double lo, double hi, double scale)
{
    if ((lo <= 0) || (hi <= 0))
        return 0;
    return hi - scale * (hi - lo);
}

void get_params(static_conf_t *sconf)
{
	int bits,i;
	double scale;
	fb_list *fb = sconf->factor_base;
    const siqs_builtin_table_t* table =
        siqs_find_param_table((uint32_t)sconf->lp_isa, sconf->lp_direct,
            (uint32_t)sconf->lp_kernel);
    const siqs_builtin_row_t* t = table->rows;
    int nrows = table->num_rows;
    double tf = 0, tfsm = 0, mfbd = 0, mfbt = 0;

	// linear interpolation according to bit size to determine
	// factor base bound.  use the closest parameter for lp multiplier
	// and number of blocks.

	bits = (int)sconf->obj->bits;
    sconf->param_table_name = table->name;

	fb->B = 0;
	if (bits <= (int)t[0].bits)
	{
		scale = (double)bits / (double)t[0].bits;
		fb->B = (uint32_t)(scale * (double)(t[0].B));
		sconf->large_mult = 40;
		sconf->num_blocks = 1;
        tf = t[0].TF;
        tfsm = t[0].TFSm;
        mfbd = t[0].MFBD;
        mfbt = t[0].MFBT;
	}
	else
	{
        // TLP rows (large_mult = 2^M, M the LPB in bits) were here once, for
        // > 1000 bits; they're kept, commented out, in siqs_table_bucket
        if (bits > 1000)
        {
            // this can still be overridden with forceDLP.
            sconf->obj->qs_obj.gbl_force_TLP = 1;
        }

        for (i = 0; i < nrows - 1; i++)
        {
            if (bits > (int)t[i].bits && bits <= (int)t[i + 1].bits)
            {
                scale = ((double)t[i + 1].bits - bits) /
                    ((double)t[i + 1].bits - t[i].bits);
                // interpolate in double: B can decrease between rows, and a
                // negative double cast to uint32_t is undefined (with AVX512F
                // gcc saturates it to 0xffffffff, B wraps to 0 and SIQS hangs
                // in new_poly_a with a 16-prime factor base)
                fb->B = (uint32_t)((double)t[i + 1].B -
                    scale * ((double)t[i + 1].B - (double)t[i].B));

                sconf->large_mult = (uint32_t)((t[i + 1].M + t[i].M) / 2.0 + 0.5);
                sconf->num_blocks = (uint32_t)((t[i + 1].NB + t[i].NB) / 2.0 + 0.5);

                tf = siqs_interp_opt(t[i].TF, t[i + 1].TF, scale);
                tfsm = siqs_interp_opt(t[i].TFSm, t[i + 1].TFSm, scale);
                mfbd = siqs_interp_opt(t[i].MFBD, t[i + 1].MFBD, scale);
                mfbt = siqs_interp_opt(t[i].MFBT, t[i + 1].MFBT, scale);
            }
        }
	}

	if (fb->B == 0)
	{
		//off the end of the table, extrapolate based on the slope of
		//the last two
        const siqs_builtin_row_t* a = &t[nrows - 2];
        const siqs_builtin_row_t* b = &t[nrows - 1];

		scale = ((double)b->B - (double)a->B) / (double)(b->bits - a->bits);
		fb->B = (uint32_t)(((double)bits - b->bits) * scale + b->B);
		sconf->large_mult = b->M;	//reuse last one

		scale = ((double)b->NB - (double)a->NB) / (double)(b->bits - a->bits);
		sconf->num_blocks = (uint32_t)(((double)bits - b->bits) * scale + b->NB);

        // the TF bound grows with the input; the rest reuse the last row
        if ((a->TF > 0) && (b->TF > 0))
        {
            scale = ((double)b->TF - (double)a->TF) / (double)(b->bits - a->bits);
            tf = ((double)bits - b->bits) * scale + b->TF;
        }
        tfsm = b->TFSm;
        mfbd = b->MFBD;
        mfbt = b->MFBT;
	}

    // make B divisible by 16
    while ((fb->B & 15) != 0)
    {
        fb->B++;
    }

    sconf->param_tf = (uint32_t)(tf + 0.5);
    sconf->param_tfsm = (uint32_t)(tfsm + 0.5);
    sconf->param_mfbd = mfbd;
    sconf->param_mfbt = mfbt;

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
