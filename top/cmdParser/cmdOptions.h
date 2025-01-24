/*
MIT License

Copyright (c) 2021 Ben Buhrow

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

// ============================================================================
// a simple command line parser.
// supports defining required and options arguments to a program.
// supports defining command line switches with optional arguments.
// ============================================================================

#ifndef CMD_PARSE_H
#define CMD_PARSE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

    // the number of recognized command line options
#define NUMOPTIONS 115
// maximum length of command line option strings
#define MAXOPTIONLEN 20
// maximum length of help string for each option
#define MAXHELPLEN 128
// maximum length of an argument to an option
#define MAXARGLEN 256

typedef struct
{
    char usageHelp[MAXHELPLEN];

    // command line options (names to be preceeded by a -)
    char OptionArray[NUMOPTIONS][MAXOPTIONLEN];
    char OptionHelp[NUMOPTIONS][MAXHELPLEN];
    char LongOptionAliases[NUMOPTIONS][MAXOPTIONLEN];

    // indication of whether or not an option needs a corresponding argument
    // 0 = no argument
    // 1 = argument required
    // 2 = argument optional
    int needsArg[NUMOPTIONS];

    // ========================================================================
    // This defines how many required and optional arguments
    // there are.  These are things sent into the program that define its
    // behavior and don't have "-" or "-- in front of them.
    int numArguments;
    int numRequired;
    // ========================================================================


    // ========================================================================
    // These variables define the required and optional arguments
    char *inputExpr;
    // ========================================================================


    // ========================================================================
    // These variables define things that change program behavior and can be
    // set via command line switches.  I.e., these need a "-" or "--" in
    // front of them on the command line... add items here as necessary

    // general options
    uint64_t rand_seed;
    char sessionlog[MAXARGLEN];
    uint32_t threads;
    uint32_t verbosity;
    uint32_t vproc;
    int num_tune_info;
    char **tune_info;
    uint32_t yafu_idle;
    char expr[MAXARGLEN];
    char scriptfile[MAXARGLEN];
    uint32_t repeat;
    uint32_t no_clk_test;
    int json_pretty;
    int obase;

    // qs options
    uint32_t siqsB;
    uint32_t siqsTF;
    uint32_t siqsTFSm;
    uint32_t siqsNobat;
    uint32_t siqsR;
    uint32_t siqsT;
    uint32_t siqsNB;
    uint32_t siqsM;
    uint32_t siqsForceDLP;
    uint32_t siqsForceTLP;
    uint32_t siqsLPB;
    uint32_t siqsBT;
    uint32_t siqsSSidx;
    uint32_t siqsSSalloc;
    double siqsMFBD;
    double siqsMFBT;
    double siqsBDiv;
    uint32_t no_opt;
    char factorlog[MAXARGLEN];
    char batchfile[MAXARGLEN];
    char qssave[MAXARGLEN];
    uint32_t inmem_cutoff;

    // nfs options
    char ggnfs_dir[MAXARGLEN];
    char nfs_outfile[MAXARGLEN];
    char nfs_jobfile[MAXARGLEN];
    char poly_method[MAXARGLEN];
    char testsieve[MAXARGLEN];
    uint32_t nfs_resume;
    uint32_t poly_batch;
    uint32_t sieveQstart, sieveQstop;
    uint32_t polystart, polystop;
    double filt_bump;
    uint32_t ggnfs_siever;
    uint32_t lathreads;
    uint32_t nc;
    uint32_t nc1;
    uint32_t nc2;
    uint32_t nc3;
    uint32_t ncr;
    uint32_t alg_side;
    uint32_t rat_side;
    uint32_t nfs_timeout;
    uint32_t force_gnfs;
    uint32_t cadoMsieve;
    char cado_dir[MAXARGLEN];
    char convert_poly_path[MAXARGLEN];
    int skip_snfscheck;
    uint32_t minrels;

    // ecm/pp1/pm1/rho/tdiv options
    uint64_t B1pm1;
    uint64_t B1pp1;
    uint64_t B1ecm;
    uint32_t rhomax;
    uint64_t B2pm1;
    uint64_t B2pp1;
    uint64_t B2ecm;
    uint64_t sigma;
    uint32_t fermat_max;
    uint32_t saveB1;
    char ecm_path[MAXARGLEN];
    uint32_t prefer_gmpecm;
    uint32_t prefer_gmpecm_stg2;
    uint32_t prefer_avxecm_stg2;
    char vpp1_work_file[MAXARGLEN];
    char vpm1_work_file[MAXARGLEN];
    char resume_file[MAXARGLEN];
    int gpucurves;
    int use_cgbn;
    int use_gpudev;
    int use_gpuecm;

    // autofactor options
    uint32_t no_ecm;
    double pretest_ratio;
    double xover;
    uint32_t one_factor;
    char opfile[MAXARGLEN];
    char offile[MAXARGLEN];
    char oufile[MAXARGLEN];
    char fact_plan[MAXARGLEN];
    uint32_t pretest;
    uint32_t want_output_expr;
    uint64_t ext_ecm_xover;
    uint32_t snfs_testsieve_threshold;
    double work;
    double qs_snfs_xover;
    int stople;
    int stoplt;
    int stopgt;
    int stopge;
    int stopbase;
    int stopeq;
    int stopprime;
    int stopk;
    int strict;
    int check_stop_conditions;

    // prime finding options
    uint32_t pfile;
    uint32_t pscreen;
    uint32_t num_prp_witnesses;
    uint32_t aprcl_p;
    uint32_t aprcl_d;
    uint32_t soe_blocksize;


    // ========================================================================

} options_t;

extern options_t* initOpt(void);
extern void applyOpt(char* opt, char* arg, options_t* options);
extern int processOpts(int argc, char** argv, options_t* options);
extern int readINI(const char* filename, options_t* options);


#ifdef __cplusplus
}
#endif

#endif /* #ifndef CMD_PARSE_H */

