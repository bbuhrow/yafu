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
//
// gcc -O2 cmdOptions.h cmdOptions.c demo.c -o demo
//
// example:
// ./demo "hello" -v -v --str "world"
// ============================================================================

#include "cmdOptions.h"
#include "ytools.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

#ifdef __INTEL_LLVM_COMPILER
#include <ctype.h>
#endif

// this function prints the help information specified by usageHelp
// and OptionHelp.
void printUsage(options_t* options);

// ========================================================================
// In this section, define the options, any long form aliases for them,
// whether they have an argument or not, and help strings to display
// if the user tries -h

// define help text for usage of the program.  
char usageHelp[MAXHELPLEN] = "usage: [optional_command_string] [options]";

// command line options visible to driver.c
char OptionArray[NUMOPTIONS][MAXOPTIONLEN] = {
    "B1pm1", "B1pp1", "B1ecm", "rhomax", "B2pm1",
    "B2pp1", "B2ecm", "qssave", "siqsB", "siqsTF",
    "siqsR", "siqsT", "siqsNB", "siqsM", "logfile",
    "batchfile", "seed", "sigma", "session", "threads",
    "v", "silent", "pfile", "pscreen", "forceDLP",
    "fmtmax", "noopt", "vproc", "noecm", "ggnfs_dir",
    "tune_info", "pretest_ratio", "xover", "one", "op",
    "of", "ou", "plan", "pretest", "no_expr",
    "o", "a", "r", "ggnfsT", "job",
    "ns", "np", "nc", "psearch", "R",
    "pbatch", "ecm_path", "siever", "ncr", "lathreads",
    "nc2", "nc3", "p", "work", "nprp",
    "ext_ecm", "testsieve", "nt", "aprcl_p", "aprcl_d",
    "filt_bump", "nc1", "gnfs", "e", "repeat",
    "ecmtime", "no_clk_test", "siqsTFSm", "script", "degree",
    "snfs_xover", "soe_block", "forceTLP", "siqsLPB", "siqsMFBD",
    "siqsMFBT", "siqsBDiv", "siqsBT", "prefer_gmpecm", "saveB1",
    "siqsNobat", "inmem", "prefer_gmpecm_stg2", "vpp1_work_file", "vpm1_work_file",
    "resume", "jsonpretty", "cadoMsieve", "cado_dir", "convert_poly_path",
    "gpucurves", "cgbn", "use_gpuecm", "use_gpudev", "prefer_avxecm_stg2",
    "stoplt", "stople", "stopeq", "stopgt", "stopge", 
    "stopbase", "stopprime", "siqsSSidx", "siqsSSalloc", "skipSNFScheck",
    "obase"};

// help strings displayed with -h
// needs to be the same length as the above arrays, even if 
// some are left blank
char OptionHelp[NUMOPTIONS][MAXHELPLEN] = {
    "(Integer < 64-bit): B1 bound of P-1 algorithm",
    "(Integer < 64-bit): B1 bound of P+1 algorithm",
    "(Integer < 64-bit): B1 bound of ECM algorithm",
    "(Integer < 32-bit): Iteration limit of Rho algorithm",
    "(Integer < 64-bit): B2 bound of P-1 algorithm",
    "(Integer < 64-bit): B2 bound of P+1 algorithm",
    "(Integer < 64-bit): B2 bound of ECM algorithm",
    "(String)          : name of SIQS save file",
    "(Integer < 32-bit): Size of SIQS factor base (number of primes)",
    "(Integer < 32-bit): SIQS trial factoring bound", 
    "(Integer < 32-bit): SIQS Number of relations to gather", 
    "(Integer < 32-bit): SIQS timeout in seconds", 
    "(Integer < 32-bit): SIQS Number of sieve blocks",
    "(Integer < 32-bit): SIQS Large prime variation multiplier",
    "(String)          : Name of factoring logfile",
    "(String)          : Name of input batchfile",
    "(Integer < 64-bit): Seed for RNG", 
    "(Integer < 64-bit): ECM sigma",
    "(String)          : Name of session logfile",
    "(Integer < 32-bit): Number of threads",
    "                  : Increase verbosity",
    "                  : Set verbosity to -1 (no output)", 
    "                  : Output primes to file primes.dat",
    "                  : Output primes to screen",
    "                  : Force use of DLP variation in SIQS",
    "(Integer < 32-bit): Max iteration for Fermat factorization method",
    "                  : Do not optimize small prime variation cutoff in SIQS", 
    "                  : Verbose output of processor info",
    "                  : Do not use ECM in factor()",
    "(String)          : Directory containing ggnfs executables",
    "(String)          : Input tuning information for SIQS/NFS crossover",
    "(Floating point)  : ECM/NFS pretesting ratio: projected-difficulty * ratio = ECM depth", 
    "(Floating point)  : SIQS/NFS crossover point (decimal digits)",
    "                  : Stop after one factor found",
    "(String)          : Filename where primes found during factor() are output",
    "(String)          : Filename where factors found during factor() are output",
    "(String)          : Filename where unfactored residues remaining after factor() are output",
    "(String)          : ECM pretesting plan (light = 2/9, deep = 1/3, normal = 4/13, none, or custom (pretest_ratio)",
    "(Integer < 32-bit): Only ECM pretest (no NFS or SIQS) to specified digit level",
    "                  : For use with op, of, or ou, to output expressions",
    "(String)          : NFS output filename",
    "                  : Algebraic side GGNFS sieving", 
    "                  : Rational side GGNFS sieving",
    "(Integer < 32-bit): NFS sieving timeout in seconds",
    "(String)          : NFS input jobfile name",
    "(Start,Stop)      : Perform only sieving phase of NFS within specified 32-bit start,stop lattice Q's",
    "(Start,Stop)      : Perform only poly-search phase of NFS within specified 32-bit start,stop leading coefficients",
    "                  : Perform only post-processing phase of NFS, starting with filtering",
    "(String)          : NFS poly-search methology ('deep', 'wide', 'fast', 'min', 'avg', or 'good')",
    "                  : Resume NFS with an existing data/job file",
    "(Integer < 32-bit): Range of leading coefficients to distribute to threads during NFS poly-search",
    "(String)          : Path of GMP-ECM executable",
    "(Integer < 32-bit): Version of GGNFS to use (11,12,13,14,15,16)",
    "                  : Restart LA phase of NFS from an existing checkpoint",
    "(Integer < 32-bit): Number of threads to use in LA (if different from -threads)",
    "                  : Perform only LA phase of NFS",
    "                  : Perform only SQRT post-processing phase of NFS",
    "                  : Set yafu to idle priority",
    "(Floating point)  : Input work level (prior ECM work to specified number of digits)",
    "(Integer < 32-bit): Number of witnesses to use in Miller-Rabin PRP tests",
    "(Integer < 64-bit): B1 bound crossover point for switching to external ECM executable",
    "(Integer < 32-bit): Difficulty level above which SNFS testsieving is automatically performed to select best poly", 
    "(String)          : Comma-delimited list of job files to GNFS test sieve", 
    "(Integer < 32-bit): setting the threshold below which numbers are proved prime using APR-CL",
    "(Integer < 32-bit): setting the threshold above which numbers that are proved prime using APR-CL have additional verbosity",
    "(Floating point)  : If NFS filtering fails to produce a matrix, increment relations by specified percentage", 
    "                  : Perform only filtering phase of NFS ",
    "                  : Force GNFS instead of SNFS", 
    "(String)          : Supplies expression to execution", 
    "(Integer < 32-bit): Repeat input command line specified number of times",
    "(Integer < 32-bit): Not currently implemented", 
    "                  : Do not test clock speed when yafu initializes", 
    "(Integer < 32-bit): SIQS small prime variation trial factoring bound",
    "(String)          : Filename of script to execute", 
    "(Integer < 32-bit): Not currently implemented",
    "(Floating point)  : QS/SNFS crossover", 
    "(Integer < 32-bit): Sieve of Eratosthenes block size in bytes", 
    "                  : Force use of TLP variation in SIQS",
    "(Integer < 32-bit): SIQS large prime bound, in bits (2^bits)",
    "(Floating point)  : Exponent of SIQS DLP: attempt to split residues up to LPB^exponent",
    "(Floating point)  : Exponent of SIQS TLP: attempt to split residues up to LPB^exponent",
    "(Floating point)  : SIQS Batch Divisor: specifies max prime to include in batch GCD as divisor of max factorbase prime",
    "(Integer < 32-bit): SIQS Batch Target: how many residues to batch before doing batch GCD", 
    "                  : Uses external GMP-ECM instead of internal AVX-ECM", 
    "                  : Create savefile with B1 residues in AVX-ECM",
    "                  : Do not use SIQS Batch GCD",
    "(Integer < 32-bit): Digit level below which SIQS in-memory is used (no savefile)",
    "                  : Use GMP-ECM for stage 2",
    "(String)          : Filename for vecP+1 work input", 
    "(String)          : Filename for vecP-1 work input",
    "(String)          : Filename to resume parallel P+1 or P-1",
    "                  : Output JSON info pretty printed (one pair per line) ",
    "                  : cadoMsieve", 
    "                  : cado_dir", 
    "                  : convert_poly_path",
    "(Integer < 32-bit): Number of gpu curves to run in a batch",
    "                  : use cgbn with gpu-ecm (external ecm binary enabled with this option must exist)",
    "                  : use gpu-ecm (external ecm binary enabled with this option must exist)",
    "(Integer < 32-bit): gpu device number",
    "                  : use AVX-ECM for stage 2",
    "(Integer < 32-bit): stop if a factor is found of less than <n> digits",
    "(Integer < 32-bit): stop if a factor is found of less than or equal to <n> digits",
    "(Integer < 32-bit): stop if a factor is found of equal to <n> digits",
    "(Integer < 32-bit): stop if a factor is found of greater than <n> digits",
    "(Integer < 32-bit): stop if a factor is found of greater than or equal to <n> digits",
    "(Integer < 32-bit): Base to use for stopXY options(default 10, range: 2 <= b <= 62)",
    "                  : Use for stopXY options to add constraint that number is prime",
    "(Integer < 32-bit): Factor base index at which to start using subset sum algorithm (default 0: unused)",
    "(Integer < 32-bit): Size of a poly-bucket, in bits.  Smaller is better as long as it's big enough. 12 is default.",
    "                  : Set this to skip all checks for the existance of SNFS polynomials for the input",
    "(Integer==8,10,16): Output base in octal (8), decimal (default, 10), or hexadecimal (16)"
};

// indication of whether or not an option needs a corresponding argument.
// needs to be the same length as the above two arrays.
// 0 = no argument
// 1 = argument required
// 2 = argument optional
int needsArg[NUMOPTIONS] = {
    1,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,1,
    0,0,0,0,0,
    1,0,0,0,1,
    1,1,1,0,1,
    1,1,1,2,0,
    1,0,0,1,1,
    2,2,0,1,0,
    1,1,1,0,1,
    0,0,0,1,1,
    1,1,1,1,1,
    1,0,0,1,1,
    1,0,1,1,1,
    1,1,0,1,1,
    1,1,1,0,0,
    0,1,0,1,1,
    1,0,0,1,1,  // resume, json-pretty, new cado options
    1,0,0,1,0,   // gpucurves, cbgn, use gpu, gpu dev, prefer avxecm stg2
    1,1,1,1,1,  // "stoplt", "stople", "stopeq", "stopgt", "stopge", 
    1,0,1,1,0,  // "stopbase", "stopprime", "siqsSSidx", "siqsSSalloc", "skipSNFScheck"
    1
};

// command line option aliases, specified by '--'
// need the same number of strings here, even if
// some of them are blank (i.e., have no long form alias).
char LongOptionAliases[NUMOPTIONS][MAXOPTIONLEN] = {
    "", "", "", "", "", 
    "", "", "", "", "", 
    "", "", "", "", "", 
    "", "", "", "", "", 
    "", "", "", "", "", 
    "", "", "", "", "", 
    "", "", "", "", "", 
    "", "", "", "", "", 
    "", "", "", "", "", 
    "", "", "", "", "", 
    "", "", "", "", "", 
    "", "", "", "", "", 
    "", "", "", "", "", 
    "", "", "", "", "", 
    "", "", "", "", "", 
    "", "", "", "", "", 
    "", "", "", "", "",
    "", "", "", "", "",
    "", "", "", "", "",
    "", "", "", "", "",
    "", "", "", "", "",
    "", "", "", "", "",
    ""
};




// ========================================================================
int enforce_numeric(char* arg, char *opt)
{
    int i;
    for (i = 0; i < (int)strlen(arg); i++)
    {
        if (!isdigit((int)arg[i]))
        {
            printf("expected numeric input for option %s\n", opt);
            exit(1);
        }
    }
    return 0;
}


// ========================================================================
// modify this function to handle input arguments
// ========================================================================
void applyArg(char* arg, int argNum, options_t* options)
{
    if (argNum == 0)
    {
        options->inputExpr = (char*)realloc(options->inputExpr, 
            (strlen(arg) + 1) * sizeof(char));
        strcpy(options->inputExpr, arg);
    }

    return;
}

// ========================================================================
// modify this function to handle each option case, assigning values
// or arguments to the members of the options_t structure.
// ========================================================================
void applyOpt(char* opt, char* arg, options_t* options)
{
    int i;
    char** ptr = NULL;

    if (strcmp(opt, OptionArray[0]) == 0)
    {
        //"B1pm1"
        enforce_numeric(arg, opt);
        options->B1pm1 = strtoull(arg, ptr, 10);
    }
    else if (strcmp(opt, OptionArray[1]) == 0)
    {
        //"B1pp1"
        enforce_numeric(arg, opt);
        options->B1pp1 = strtoull(arg, ptr, 10);
    }
    else if (strcmp(opt, OptionArray[2]) == 0)
    {
        //"B1ecm"
        enforce_numeric(arg, opt);
        options->B1ecm = strtoull(arg, ptr, 10);
    }
    else if (strcmp(opt, OptionArray[3]) == 0)
    {
        //"rhomax"
        enforce_numeric(arg, opt);
        options->rhomax = strtoul(arg, ptr, 10);
    }
    else if (strcmp(opt, OptionArray[4]) == 0)
    {
        //"B2pm1"
        enforce_numeric(arg, opt);
        options->B2pm1 = strtoull(arg, ptr, 10);
    }
    else if (strcmp(opt, OptionArray[5]) == 0)
    {
        // "B2pp1"
        enforce_numeric(arg, opt);
        options->B2pp1 = strtoull(arg, ptr, 10);
    }
    else if (strcmp(opt, OptionArray[6]) == 0)
    {
        // "B2ecm"
        enforce_numeric(arg, opt);
        options->B2ecm = strtoull(arg, ptr, 10);
    }
    else if (strcmp(opt, OptionArray[7]) == 0)
    {
        // qssave
        if (strlen(arg) < MAXARGLEN)
            strcpy(options->qssave, arg);
        else
            printf("*** argument to savefile too long, ignoring ***\n");
    }
    else if (strcmp(opt, OptionArray[8]) == 0)
    {
        // siqsB
        enforce_numeric(arg, opt);
        options->siqsB = strtoul(arg, ptr, 10);
    }
    else if (strcmp(opt, OptionArray[9]) == 0)
    {
        // siqsTF
        enforce_numeric(arg, opt);
        options->siqsTF = strtoul(arg, ptr, 10);
    }
    else if (strcmp(opt, OptionArray[10]) == 0)
    {
        // siqsR
        enforce_numeric(arg, opt);
        options->siqsR = strtoul(arg, ptr, 10);
    }
    else if (strcmp(opt, OptionArray[11]) == 0)
    {
        // siqsT
        enforce_numeric(arg, opt);
        options->siqsT = strtoul(arg, ptr, 10);
    }
    else if (strcmp(opt, OptionArray[12]) == 0)
    {
        // siqsNB
        enforce_numeric(arg, opt);
        options->siqsNB = strtoul(arg, ptr, 10);
    }
    else if (strcmp(opt, OptionArray[13]) == 0)
    {
        // siqsM
        enforce_numeric(arg, opt);
        options->siqsM = strtoul(arg, ptr, 10);
    }
    else if (strcmp(opt, OptionArray[14]) == 0)
    {
        //argument is a string
        if ((strlen(arg) == 0) || (strcmp(arg, "NUL") == 0) ||
            (strcmp(arg, "NULL") == 0) || (strcmp(arg, "nul") == 0) ||
            (strcmp(arg, "null") == 0))
        {
            strcpy(options->factorlog, "");
        }
        else
        {
            if (strlen(arg) < MAXARGLEN)
                strcpy(options->factorlog, arg);
            else
                printf("*** argument to logfile too long, ignoring ***\n");
        }
    }
    else if (strcmp(opt, OptionArray[15]) == 0)
    {
        //argument is a string
        if (strlen(arg) < MAXARGLEN)
        {
            strcpy(options->batchfile, arg);
        }
        else
            printf("*** argument to batchfile too long, ignoring ***\n");
    }
    else if (strcmp(opt, OptionArray[16]) == 0)
    {
        //options->rand_seed = strtoull(arg, ptr, 10);
        printf("attempting to parse user rng seed from %s\n", arg);
        sscanf(arg, "%lu", &options->rand_seed);
        printf("read: % "PRIu64"\n", options->rand_seed);
    }
    else if (strcmp(opt, OptionArray[17]) == 0)
    {
        enforce_numeric(arg, opt);
        options->sigma = strtoull(arg, NULL, 10);
    }
    else if (strcmp(opt, OptionArray[18]) == 0)
    {
        //argument is a string
        if (strlen(arg) < MAXARGLEN)
        {
            strcpy(options->sessionlog, arg);
        }
        else
        {
            printf("*** argument to sessionname too long, ignoring ***\n");
        }
    }
    else if (strcmp(opt, OptionArray[19]) == 0)
    {
        enforce_numeric(arg, opt);
        options->threads = strtoul(arg, NULL, 10);
    }
    else if (strcmp(opt, OptionArray[20]) == 0)
    {
        options->verbosity++;
    }
    else if (strcmp(opt, OptionArray[21]) == 0)
    {
        options->verbosity = -1;
    }
    else if (strcmp(opt, OptionArray[22]) == 0)
    {
        options->pfile = 1;
    }
    else if (strcmp(opt, OptionArray[23]) == 0)
    {
        options->pscreen = 1;
    }
    else if (strcmp(opt, OptionArray[24]) == 0)
    {
        options->siqsForceDLP = 1;
    }
    else if (strcmp(opt, OptionArray[25]) == 0)
    {
        enforce_numeric(arg, opt);
        options->fermat_max = strtoul(arg, NULL, 10);
    }
    else if (strcmp(opt, OptionArray[26]) == 0)
    {
        options->no_opt = 1;
    }
    else if (strcmp(opt, OptionArray[27]) == 0)
    {
        options->vproc++;
    }
    else if (strcmp(opt, OptionArray[28]) == 0)
    {
        strcpy(options->fact_plan, "noecm");
    }
    else if (strcmp(opt, OptionArray[29]) == 0)
    {
        //argument is a string
        if (strlen(arg) < MAXARGLEN)
        {
            strcpy(options->ggnfs_dir, arg);
        }
        else
        {
            printf("*** argument to ggnfs_dir too long, ignoring ***\n");
        }
            
    }
    else if (strcmp(opt, OptionArray[30]) == 0)
    {
        if (strlen(arg) < MAXARGLEN)
        {
            options->num_tune_info++;
            options->tune_info = (char**)xrealloc(options->tune_info, options->num_tune_info * sizeof(char*));
            options->tune_info[options->num_tune_info - 1] = (char*)xmalloc(MAXARGLEN * sizeof(char));
            strcpy(options->tune_info[options->num_tune_info-1], arg);
        }
        else
        {
            printf("*** argument to tune_info too long, ignoring ***\n");
        }
    }
    else if (strcmp(opt, OptionArray[31]) == 0)
    {
        //argument "pretest_ratio"
        sscanf(arg, "%lf", &options->pretest_ratio);
    }
    else if (strcmp(opt, OptionArray[32]) == 0)
    {
        //argument "xover"
        sscanf(arg, "%lf", &options->xover);
        //fobj->nfs_obj.min_digits = fobj->autofact_obj.qs_gnfs_xover;
        //fobj->autofact_obj.prefer_xover = 1;
    }
    else if (strcmp(opt, OptionArray[33]) == 0)
    {
        //argument "one"
        options->one_factor = 1;
    }
    else if (strcmp(opt, OptionArray[34]) == 0)
    {
        //argument "op".  argument is a string
        if (strlen(arg) < MAXARGLEN)
        {
            strcpy(options->opfile, arg);
        }
        else
        {
            printf("*** argument to -op too long, ignoring ***\n");
        }
    }
    else if (strcmp(opt, OptionArray[35]) == 0)
    {
        //argument "of".  argument is a string
        if (strlen(arg) < MAXARGLEN)
        {
            strcpy(options->offile, arg);
        }
        else
        {
            printf("*** argument to -of too long, ignoring ***\n");
        }
    }
    else if (strcmp(opt, OptionArray[36]) == 0)
    {
        //argument "ou".  argument is a string
        if (strlen(arg) < MAXARGLEN)
        {
            strcpy(options->oufile, arg);
        }
        else
        {
            printf("*** argument to -ou too long, ignoring ***\n");
        }
    }
    else if (strcmp(opt, OptionArray[37]) == 0)
    {
        //argument "plan".  argument is a string
        if (strlen(arg) < MAXARGLEN)
        {
            // test for recognized options.  
            if (strcmp(arg, "none") == 0)
                strcpy(options->fact_plan, "none");
            else if (strcmp(arg, "noecm") == 0)
                strcpy(options->fact_plan, "noecm");
            else if (strcmp(arg, "light") == 0)
                strcpy(options->fact_plan, "light");
            else if (strcmp(arg, "deep") == 0)
                strcpy(options->fact_plan, "deep");
            else if (strcmp(arg, "normal") == 0)
                strcpy(options->fact_plan, "normal");
            else if (strcmp(arg, "custom") == 0)
                strcpy(options->fact_plan, "custom");
            else
            {
                printf("*** unknown plan option, ignoring ***\n");
                strcpy(options->fact_plan, "normal");
            }
        }
        else
        {
            printf("*** argument to -plan too long, ignoring ***\n");
        }
    }
    else if (strcmp(opt, OptionArray[38]) == 0)
    {
        //argument "pretest"
        if (arg == NULL)
        {
            // no argument, use the value "1" to signify doing 
            // pretesting to the bounds computed by factor()
            options->pretest = 1;
        }
        else
        {
            // an optional argument to pretest is interpreted as
            // a maximum t-level to pretest to
            options->pretest = strtoul(arg, NULL, 10);
        }
    }
    else if (strcmp(opt, OptionArray[39]) == 0)
    {
        //argument "no_expr"
        options->want_output_expr = 0;
    }
    else if (strcmp(opt, OptionArray[40]) == 0)
    {
        //argument "o".  Indicates output filename ggnfs sieving.
        if (strlen(arg) < MAXARGLEN)
        {
            strcpy(options->nfs_outfile, arg);
        }
        else
        {
            printf("*** argument to -o too long, ignoring ***\n");
        }
    }
    else if (strcmp(opt, OptionArray[41]) == 0)
    {
        //argument "a".  Indicates algebraic side special Q.
        options->alg_side = 1;
    }
    else if (strcmp(opt, OptionArray[42]) == 0)
    {
        //argument "r".  Indicates rational side special Q.
       options->rat_side = 1;
    }
    else if (strcmp(opt, OptionArray[43]) == 0)
    {
        //argument "ggnfsT".  Indicates timeout (in seconds) for NFS job.
        options->nfs_timeout = strtoul(arg, NULL, 10);
    }
    else if (strcmp(opt, OptionArray[44]) == 0)
    {
        //argument "job".  Indicates input .job file automated NFS.
        if (strlen(arg) < MAXARGLEN)
        {
            strcpy(options->nfs_jobfile, arg);
        }
        else
        {
            printf("*** argument to -job too long, ignoring ***\n");
        }
    }
    else if (strcmp(opt, OptionArray[45]) == 0)
    {
        // argument "ns".  do nfs sieving
        char** nextptr = &arg;

        if (arg != NULL)
        {
            // if an argument was supplied, parse the start and range of 
            //special Q in the format X,Y
            options->sieveQstart = strtoul(arg, nextptr, 10);

            if (*nextptr[0] != ',')
            {
                printf("format of sieving argument is START,STOP\n");
                exit(1);
            }
            options->sieveQstop = strtoul(*nextptr + 1, NULL, 10);
        }
        else
        {
            options->sieveQstart = 1;
            options->sieveQstop = 1;
        }

    }
    else if (strcmp(opt, OptionArray[46]) == 0)
    {
        char** nextptr = &arg;

        //argument "np".  do poly finding.

        if (arg != NULL)
        {
            // if an argument was supplied, parse the start and stop coefficient range in the
            // format X,Y
            options->polystart = strtoul(arg, nextptr, 10);

            if (*nextptr[0] != ',')
            {
                printf("format of poly select argument is START,STOP\n");
                exit(1);
            }
            options->polystop = strtoul(*nextptr + 1, NULL, 10);
        }
        else
        {
            options->polystart = 1;
            options->polystop = 1;
        }
    }
    else if (strcmp(opt, OptionArray[47]) == 0)
    {
        //argument "nc".  Do post processing, starting with filtering
        options->nc = 1;
    }
    else if (strcmp(opt, OptionArray[48]) == 0)
    {
        //argument "psearch".  modify poly search methodology
        if (strlen(arg) < MAXARGLEN)
        {
            strcpy(options->poly_method, arg);

            if (strcmp(arg, "wide") == 0)
                strcpy(options->poly_method, arg); // fobj->nfs_obj.poly_option = 1;
            else if (strcmp(arg, "deep") == 0)
                strcpy(options->poly_method, arg); //fobj->nfs_obj.poly_option = 2;
            else if (strcmp(arg, "fast") == 0)
                strcpy(options->poly_method, arg); //fobj->nfs_obj.poly_option = 0;
            else if (strcmp(arg, "min") == 0)
                strcpy(options->poly_method, arg); //fobj->nfs_obj.poly_option = 3;
            else if (strcmp(arg, "avg") == 0)
                strcpy(options->poly_method, arg); //fobj->nfs_obj.poly_option = 4;
            else if (strcmp(arg, "good") == 0)
                strcpy(options->poly_method, arg); //fobj->nfs_obj.poly_option = 5;
            else
            {
                printf("option -psearch recognizes arguments 'deep', 'wide', 'fast', 'min', 'avg', or 'good'.\n  see docfile.txt for details\n");
                exit(1);
            }

        }
        else
        {
            printf("*** argument to -psearch too long, ignoring ***\n");
        }

    }
    else if (strcmp(opt, OptionArray[49]) == 0)
    {
        //argument "R".  nfs restart flag
        options->nfs_resume = 1;
    }
    else if (strcmp(opt, OptionArray[50]) == 0)
    {
        //argument "pbatch".  Indicates size of blocks of leading coefficients to
        //distribute to each thread in threaded NFS poly selection.
        options->poly_batch = strtoul(arg, NULL, 10);
        if (options->poly_batch == 0)
            options->poly_batch = 250;
    }
    else if (strcmp(opt, OptionArray[51]) == 0)
    {
        // argument "ecm_path"
        //argument is a string
        if (strlen(arg) < MAXARGLEN)
            strcpy(options->ecm_path, arg);
        else
            printf("*** argument to ecm_path too long, ignoring ***\n");
    }
    else if (strcmp(opt, OptionArray[52]) == 0)
    {
        // argument "siever"
        enforce_numeric(arg, opt);
        options->ggnfs_siever = strtoul(arg, NULL, 10);
    }
    else if (strcmp(opt, OptionArray[53]) == 0)
    {
        //argument "ncr".  linear algebra restart flag
        options->ncr = 1;
    }
    else if (strcmp(opt, OptionArray[54]) == 0)
    {
        // argument "lathreads"
        enforce_numeric(arg, opt);
        options->lathreads = strtoul(arg, NULL, 10);
    }
    else if (strcmp(opt, OptionArray[55]) == 0)
    {
        //argument "nc2".  do linear algebra.
        options->nc2 = 1;
    }
    else if (strcmp(opt, OptionArray[56]) == 0)
    {
        //argument "nc3".  do nfs sqrt
        options->nc3 = 1;
    }
    else if (strcmp(opt, OptionArray[57]) == 0)
    {
        //argument "p".  set to idle priority.
        options->yafu_idle = 1;
    }
    else if (strcmp(opt, OptionArray[58]) == 0)
    {
        //argument "work"
        sscanf(arg, "%lf", &options->work);
    }
    else if (strcmp(opt, OptionArray[59]) == 0)
    {
        //argument "nprp"
        options->num_prp_witnesses = strtoul(arg, NULL, 10);
    }
    else if (strcmp(opt, OptionArray[60]) == 0)
    {
        // argument "ext_ecm"
        enforce_numeric(arg, opt);
        options->ext_ecm_xover = strtoull(arg, NULL, 10);
    }
    else if (strcmp(opt, OptionArray[61]) == 0)
    {
        //argument "testsieve"
        options->snfs_testsieve_threshold = strtoul(arg, NULL, 10);
    }
    else if (strcmp(opt, OptionArray[62]) == 0)
    {
        // argument "nt"
        if (arg == NULL)
        {
            printf("expected argument for option %s\n", opt);
            exit(1);
        }
        else if (strlen(arg) < MAXARGLEN)
        {
            printf("*** argument to nt too long, ignoring ***\n");
        }
        else
        {
            strcpy(options->testsieve, arg);
        }
    }
    else if (strcmp(opt, OptionArray[63]) == 0)
    {
        // argument "aprcl_p", setting the threshold below which numbers
        // are proved prime using APR-CL
        options->aprcl_p = strtoul(arg, NULL, 10);
        if (options->aprcl_p > 6021)
        {
            printf("APR-CL primality proving is possible only for numbers less"
                " than 6021 digits... setting limit to 6021 digits\n");
            options->aprcl_p = 6021;
        }
    }
    else if (strcmp(opt, OptionArray[64]) == 0)
    {
        // argument "aprcl_d", setting the threshold above which numbers
        // that are proved prime using APR-CL have additional verbosity enabled
        options->aprcl_d = strtoul(arg, NULL, 10);
    }
    else if (strcmp(opt, OptionArray[65]) == 0)
    {
        //argument "filt_bump"
        sscanf(arg, "%lf", &options->filt_bump);
    }
    else if (strcmp(opt, OptionArray[66]) == 0)
    {
        //argument "nc1".  do msieve filtering.
        options->nc1 = 1;
    }
    else if (strcmp(opt, OptionArray[67]) == 0)
    {
        //argument "gnfs"
        options->force_gnfs = 1;
    }
    else if (strcmp(opt, OptionArray[69]) == 0)			// NOTE: we skip -e because it is handled 
    {													// specially by process_flags
        //argument "repeat"
        options->repeat = strtoul(arg, NULL, 10);
    }
    else if (strcmp(opt, OptionArray[70]) == 0)
    {
        //argument "ecmtime"

    }
    else if (strcmp(opt, OptionArray[71]) == 0)
    {
        //argument "no_clk_test"
        options->no_clk_test = 1;
    }
    else if (strcmp(opt, OptionArray[72]) == 0)
    {
        //argument "siqsTFSm"
        options->siqsTFSm = atoi(arg);
    }
    else if (strcmp(opt, OptionArray[73]) == 0)
    {
        //argument "script"
        sscanf(arg, "%s", options->scriptfile);
    }
    else if (strcmp(opt, OptionArray[74]) == 0)
    {
        //argument "degree"

    }
    else if (strcmp(opt, OptionArray[75]) == 0)
    {
        //argument "snfs_xover"
        sscanf(arg, "%lf", &options->qs_snfs_xover);
        //fobj->autofact_obj.prefer_xover = 1;
    }
    else if (strcmp(opt, OptionArray[76]) == 0)
    {
        //argument "soe_block"
        options->soe_blocksize = strtoul(arg, NULL, 10);
    }
    else if (strcmp(opt, OptionArray[77]) == 0)
    {
        //argument "forceTLP"
        options->siqsForceTLP = 1;
    }
    else if (strcmp(opt, OptionArray[78]) == 0)
    {
        //argument "siqsLPB"
        // the maximum allowed large prime, in bits, in SIQS
        sscanf(arg, "%d", &options->siqsLPB);
    }
    else if (strcmp(opt, OptionArray[79]) == 0)
    {
        //argument "siqsMFBD"
        // the exponent of the large prime bound such that residues larger than
        // lpb^siqsMFBD are subjected to double large prime factorization attempts
        sscanf(arg, "%lf", &options->siqsMFBD);
    }
    else if (strcmp(opt, OptionArray[80]) == 0)
    {
        //argument "siqsMFBT"
        // the exponent of the large prime bound such that residues larger than
        // lpb^siqsMFBT are subjected to triple large prime factorization attempts
        sscanf(arg, "%lf", &options->siqsMFBT);
    }
    else if (strcmp(opt, OptionArray[81]) == 0)
    {
        //argument "siqsBDiv"
        // The divider of large_prime_max as the upper bound for
        // primes to multiply when using batch GCD in TLP factorizations.
        //fobj->qs_obj.gbl_override_bdiv_flag = 1;
        sscanf(arg, "%lf", &options->siqsBDiv);
    }
    else if (strcmp(opt, OptionArray[82]) == 0)
    {
        //argument "siqsBT" ("Batch Target")
        // How many relations to batch up before they are processed
        sscanf(arg, "%u", &options->siqsBT);
    }
    else if (strcmp(opt, OptionArray[83]) == 0)
    {
        // argument "prefer_gmpecm"
        options->prefer_gmpecm = 1;
    }
    else if (strcmp(opt, OptionArray[84]) == 0)
    {
        //argument "save_b1"
        options->saveB1 = 1;
    }
    else if (strcmp(opt, OptionArray[85]) == 0)
    {
        //argument "siqsNobat"
        // Whether or not to use batch factoring in 3LP
        options->siqsNobat = 1;
    }
    else if (strcmp(opt, OptionArray[86]) == 0)
    {
        // argument "inmem"
        // cutoff for processing in-memory
        sscanf(arg, "%u", &options->inmem_cutoff);
    }
    else if (strcmp(opt, OptionArray[87]) == 0)
    {
        // argument "prefer_gmpecm_stg2"
        options->prefer_gmpecm_stg2 = 1;
        options->prefer_avxecm_stg2 = 0;
    }
    else if (strcmp(opt, OptionArray[88]) == 0)
    {
        // argument "vpp1_work_file"
        if (strlen(arg) < MAXARGLEN)
            strcpy(options->vpp1_work_file, arg);
        else
            printf("*** argument to vpp1_work_file too long, ignoring ***\n");
    }
    else if (strcmp(opt, OptionArray[89]) == 0)
    {
        // argument "vpm1_work_file"
        if (strlen(arg) < MAXARGLEN)
            strcpy(options->vpm1_work_file, arg);
        else
            printf("*** argument to vpm1_work_file too long, ignoring ***\n");
    }
    else if (strcmp(opt, OptionArray[90]) == 0)
    {
        if (strlen(arg) < MAXARGLEN)
            strcpy(options->resume_file, arg);
        else
            printf("*** argument to resume too long, ignoring ***\n");
    }
    else if (strcmp(opt, OptionArray[91]) == 0)
    {
        options->json_pretty = 1;
    }
    else if (strcmp(opt, OptionArray[92]) == 0)
    {
        // argument "cadoMsieve"
        options->cadoMsieve = 1;
    }
    else if (strcmp(opt, OptionArray[93]) == 0)
    {
        // argument "cado_dir"
        if (strlen(arg) < MAXARGLEN)
            strcpy(options->cado_dir, arg);
        else
            printf("*** argument to cado_dir too long, ignoring ***\n");
    }
    else if (strcmp(opt, OptionArray[94]) == 0)
    {
        // argument "convert_poly_path"
        if (strlen(arg) < MAXARGLEN)
            strcpy(options->convert_poly_path, arg);
        else
            printf("*** argument to convert_poly_path too long, ignoring ***\n");
    }
    else if (strcmp(opt, OptionArray[95]) == 0)
    {
        // argument "gpucurves"
        options->gpucurves = atoi(arg);
    }
    else if (strcmp(opt, OptionArray[96]) == 0)
    {
        // argument "use_cgbn"
        options->use_cgbn = 1;
    }
    else if (strcmp(opt, OptionArray[97]) == 0)
    {
        // argument "use_gpuecm"
        options->use_gpuecm = 1;
    }
    else if (strcmp(opt, OptionArray[98]) == 0)
    {
        // argument "use_gpudev"
        options->use_gpudev = atoi(arg);
    }
    else if (strcmp(opt, OptionArray[99]) == 0)
    {
        // argument "prefer_avxecm_stg2"
        options->prefer_avxecm_stg2 = 1;
        options->prefer_gmpecm_stg2 = 0;
    }
    else if (strcmp(opt, OptionArray[100]) == 0)
    {
        // argument "stoplt"
        options->stoplt = atoi(arg);
        options->check_stop_conditions = 1;
    }
    else if (strcmp(opt, OptionArray[101]) == 0)
    {
        // argument "stople"
        options->stople = atoi(arg);
        options->check_stop_conditions = 1;
    }
    else if (strcmp(opt, OptionArray[102]) == 0)
    {
        // argument "stopeq"
        options->stopeq = atoi(arg);
        options->check_stop_conditions = 1;
    }
    else if (strcmp(opt, OptionArray[103]) == 0)
    {
        // argument "stopgt"
        options->stopgt = atoi(arg);
        options->check_stop_conditions = 1;
    }
    else if (strcmp(opt, OptionArray[104]) == 0)
    {
        // argument "stopge"
        options->stopge = atoi(arg);
        options->check_stop_conditions = 1;
    }
    else if (strcmp(opt, OptionArray[105]) == 0)
    {
        // argument "stopbase"
        options->stopbase = atoi(arg);
    }
    else if (strcmp(opt, OptionArray[106]) == 0)
    {
        // argument "stopprime"
        options->stopprime = 1;
    }
    else if (strcmp(opt, OptionArray[107]) == 0)
    {
        // argument "siqsSSidx"
        options->siqsSSidx = atoi(arg);
    }
    else if (strcmp(opt, OptionArray[108]) == 0)
    {
        // argument "siqsSSalloc"
        options->siqsSSalloc = atoi(arg);
    }
    else if (strcmp(opt, OptionArray[109]) == 0)
    {
        // argument "skipSNFScheck"
        options->skip_snfscheck = 1;
    }
    else if (strcmp(opt, OptionArray[110]) == 0)
    {
        // argument "OBASE"
        options->obase = atoi(arg);
        if ((options->obase != 8) && (options->obase != 10) && (options->obase != 16))
        {
            printf("*** argument obase must be either 8,10,or 16 ***\n");
            options->obase = 10;
        }
 
    }
    else
    {
        int i;

        printf("invalid option %s\n", opt);
        printf("supported options are: ");
        for (i = 0; i < NUMOPTIONS; i++)
        {
            if (i % 5 == 0) printf("\n");
            printf("%s ", options->OptionArray[i]);
        }
        exit(0);
    }

    return;
}

// ========================================================================
// modify this function to handle assign a default value to each 
// member of the options_t structure, if desired.
// ========================================================================
options_t* initOpt(void)
{
    options_t* options = (options_t*)malloc(sizeof(options_t));
    int i;

    for (i = 0; i < NUMOPTIONS; i++)
    {
        strcpy(options->OptionArray[i], OptionArray[i]);
        strcpy(options->OptionHelp[i], OptionHelp[i]);
        strcpy(options->LongOptionAliases[i], LongOptionAliases[i]);
        options->needsArg[i] = needsArg[i];
    }

    // ========================================================================
    // define how many required and optional arguments there are
    options->numArguments = 1;
    options->numRequired = 0;
    // ========================================================================


    // ========================================================================
    // default values assigned to optional arguments here
    options->inputExpr = (char*)malloc(MAXARGLEN * sizeof(char));
    strcpy(options->inputExpr, "");
    // ========================================================================


    // ========================================================================
    // default values assigned to switches here:
    
    // general options
    options->yafu_idle = 0;
    strcpy(options->factorlog, "factor.log");
    strcpy(options->batchfile, "");
    strcpy(options->sessionlog, "session.log");
    strcpy(options->scriptfile, "");
    options->num_tune_info = 0;
    options->tune_info = (char**)xmalloc(1 * sizeof(char*));
    options->rand_seed = 0;
    options->threads = 1;
    options->verbosity = 0;
    options->vproc = 0;
    strcpy(options->expr, "");
    options->repeat = 0;
    options->no_clk_test = 1;
    options->json_pretty = 0;
    options->obase = 10;

    // autofact options
    options->no_ecm = 0;
    options->pretest_ratio = 4.0 / 13.0;
    options->xover = 95;
    options->one_factor = 0;
    strcpy(options->opfile, "");
    strcpy(options->offile, "");
    strcpy(options->oufile, "");
    strcpy(options->fact_plan, "normal");
    options->pretest = 0;
    options->want_output_expr = 1;
    strcpy(options->tune_info, "");
    options->stopbase = 10;
    options->stopeq = -1;
    options->stople = -1;
    options->stoplt = -1;
    options->stopge = -1;
    options->stopgt = -1;
    options->stopprime = 0;
    options->check_stop_conditions = 0;
    
    // nfs options
    strcpy(options->nfs_outfile, "nfs.dat");
    options->alg_side = 0;
    options->rat_side = 0;
    options->nfs_timeout = 0;
    strcpy(options->nfs_jobfile, "nfs.job");
    options->sieveQstart = 0;
    options->sieveQstop = 0;
    options->polystart = 0;
    options->polystop = 0;
    options->nc = 0;
    options->nc1 = 0;
    options->nc2 = 0;
    options->nc3 = 0;
    options->ncr = 0;
    options->lathreads = 1;
    options->filt_bump = 5;
    options->force_gnfs = 0;
    options->qs_snfs_xover = 95;
    options->snfs_testsieve_threshold = 160;
    strcpy(options->testsieve, "");
    strcpy(options->poly_method, "avg");
    options->nfs_resume = 0;
    options->poly_batch = 250;
    options->ggnfs_siever = 0;
#if defined(_WIN64)
    strcpy(options->ggnfs_dir, ".\\");
#elif defined(WIN32)
    strcpy(options->ggnfs_dir, ".\\");
#else
    strcpy(options->ggnfs_dir, "./");
#endif
    options->skip_snfscheck = 0;
    
    // prime finding options
    options->soe_blocksize = 32768;
    options->aprcl_p = 500;
    options->aprcl_d = 200;
    options->num_prp_witnesses = 1;
    options->pfile = 0;
    options->pscreen = 0;

    // siqs options
    strcpy(options->qssave, "siqs.dat");
    options->siqsB = 0;
    options->siqsTF = 0;
    options->siqsR = 0;
    options->siqsT = 0;
    options->siqsNB = 0;
    options->siqsM = 0;
    options->siqsForceDLP = 0;
    options->siqsTFSm = 0;
    options->siqsNobat = 0;
    options->siqsSSidx = 0;
    options->siqsSSalloc = 12;
    options->siqsForceTLP = 0;
    options->siqsLPB = 0;
    options->siqsMFBD = 1.85;
    options->siqsMFBT = 2.9;
    options->siqsBDiv = 3.0;
    // larger batches of relations in TLP are more efficient
    // to run in the batch GCD:
    // 100k: 11314 rels/sec 
    // 500k: 19230 rels/sec
    // 1000k: 25784 rels/sec
    // but there are a couple tradeoffs.  One, the GCD uses
    // more memory and two, it both takes longer to gather
    // the batch and the processing takes longer.  These two
    // issues can significantly extend runtimes toward the
    // end of a factorization when not many more relations are
    // needed and when the cycle formation rate is very high.
    // so we start this out fairly high and gradually decrease it.
    options->siqsBT = 1000000;
    options->no_opt = 0;
    options->inmem_cutoff = 70;

    // ecm options
    strcpy(options->ecm_path, "");
    options->work = 0;
    options->sigma = 0;
    // if we can use AVX-ECM, do so, and save B1 checkpoints.
#ifdef USE_AVX512F
    options->prefer_gmpecm = 0;
    options->prefer_gmpecm_stg2 = 0;
    options->saveB1 = 0;
    options->ext_ecm_xover = 300000000;
#else
    options->prefer_gmpecm = 1;
    options->prefer_gmpecm_stg2 = 1;
    options->saveB1 = 0;
    options->ext_ecm_xover = 48000;
#endif
    options->prefer_avxecm_stg2 = 0;
    options->B1pm1 = 100000;
    options->B1pp1 = 20000;
    options->B1ecm = 11000;
    options->rhomax = 1000;
    options->B2pm1 = 0;
    options->B2pp1 = 0;
    options->B2ecm = 0;
    options->fermat_max = 1000000;
    strcpy(options->vpp1_work_file, "pp1_work.ini");
    strcpy(options->vpm1_work_file, "pm1_work.ini");
    strcpy(options->resume_file, "");
    
    // ========================================================================
    return options;
}


// ========================================================================
// this function should not need to be changed
// ========================================================================
int processOpts(int argc, char** argv, options_t* options)
{
    // this function takes the command line arguments and returns
    // an options_t, which will have any valid options specified
    // on the command line filled in.  Return the number of options
    // processed.
    int i, j, k = 0, valid, longswitch = 0;
    char optbuf[MAXOPTIONLEN];
    char argbuf[MAXARGLEN];
    int numOpt = 0;

    //argument loop
    i = 1;
    while (i < argc)
    {
        longswitch = 0;

        // read in the option
        if ((strlen(argv[i]) > 1) && (argv[i][0] == '-') && (argv[i][1] == '-'))
        {
            longswitch = 1;
        }
        else if (argv[i][0] == '-')
        {
            longswitch = 0;
        }
        else
        {
            if (numOpt >= (options->numArguments))
            {
                printf("Too many inputs\n\n");
                printUsage(options);
                exit(0);
            }
            else
            {
                // process argument
                applyArg(argv[i], numOpt, options);
                numOpt++;
                i++;
                continue;
            }
        }

        if (numOpt < options->numRequired)
        {
            if (options->numRequired == 1)
            {
                printf("Missing %d required argument \n\n", options->numRequired);
            }
            else
            {
                printf("Missing %d required arguments \n\n", options->numRequired);
            }
            printUsage(options);
            exit(0);
        }

        // check if the options is valid
        valid = 0;
        if (longswitch)
        {
            for (j = 0; j < NUMOPTIONS; j++)
            {
                if (strncmp(options->LongOptionAliases[j], &argv[i][2], MAXOPTIONLEN) == 0)
                {
                    valid = 1;
                    strncpy(optbuf, options->OptionArray[j], MAXOPTIONLEN);
                    break;
                }
            }
        }
        else
        {
            for (j = 0; j < NUMOPTIONS; j++)
            {
                if (strncmp(options->OptionArray[j], &argv[i][1], MAXOPTIONLEN) == 0)
                {
                    valid = 1;
                    strncpy(optbuf, &argv[i][1], MAXOPTIONLEN);
                    break;
                }
            }
        }

        if (valid == 0)
        {
            printUsage(options);
            exit(0);
        }

        //check to see if this option requires an argument
        if (options->needsArg[j] == 1)
        {
            i++;
            if ((i == argc) || argv[i][0] == '-')
            {
                printf("argument expected for %s%s\n\n", longswitch ? "--" : "-", optbuf);
                printUsage(options);
                exit(0);
            }
            strncpy(argbuf, argv[i], MAXARGLEN);

            //now apply -option argument
            applyOpt(optbuf, argbuf, options);
            k++;
        }
        else if (options->needsArg[j] == 2)
        {
            // check to see if an argument was supplied
            if (((i + 1) == argc) || argv[i + 1][0] == '-')
            {
                // no option supplied.  use default option
                applyOpt(optbuf, NULL, options);
                k++;
            }
            else
            {
                i++;
                // an option was supplied, pass it on
                strncpy(argbuf, argv[i], MAXARGLEN);

                //now apply -option argument
                applyOpt(optbuf, argbuf, options);
                k++;
            }
        }
        else
        {
            //apply -option
            //now apply -option argument
            applyOpt(optbuf, NULL, options);
            k++;
        }
        i++;
    }

    if (numOpt < options->numRequired)
    {
        if (options->numRequired == 1)
        {
            printf("Missing %d required argument \n\n", options->numRequired);
        }
        else
        {
            printf("Missing %d required arguments \n\n", options->numRequired);
        }
        printUsage(options);
        exit(0);
    }

    return k;
}

// ========================================================================
// this function should not need to be changed:
// edit OptionHelp and usageHelp strings.
// ========================================================================
void printUsage(options_t* options)
{
    int j;
    printf("Usage: %s\n", usageHelp);
    printf("supported options are: \n");
    for (j = 0; j < NUMOPTIONS; j++)
    {
        int n;
        if (strlen(options->LongOptionAliases[j]) > 0)
        {
            if (options->needsArg[j])
            {
                n = printf("%s <value>",
                    options->OptionArray[j]);
                for (; n < MAXOPTIONLEN + 8; n++)
                    printf(" ");
                printf(": %s (alias --%s)\n",
                    options->OptionHelp[j],
                    options->LongOptionAliases[j]);
            }
            else
            {
                n = printf("%s        ",
                    options->OptionArray[j]);
                for (; n < MAXOPTIONLEN + 8; n++)
                    printf(" ");
                printf(": %s (alias --%s)\n",
                    options->OptionHelp[j],
                    options->LongOptionAliases[j]);
            }
        }
        else
        {
            if (options->needsArg[j])
            {
                n = printf("%s <value>",
                    options->OptionArray[j]);
                for (; n < MAXOPTIONLEN + 8; n++)
                    printf(" ");
                printf(": %s\n", options->OptionHelp[j]);
            }
            else
            {
                n = printf("%s        ",
                    options->OptionArray[j]);
                for (; n < MAXOPTIONLEN + 8; n++)
                    printf(" ");
                printf(": %s\n", options->OptionHelp[j]);
            }
        }
    }
    return;
}

// ========================================================================
// this function should not need to be changed:
// parse options from a .ini (or other) file.
// ========================================================================
int readINI(const char* filename, options_t* options)
{
    FILE* doc;
    char* str;
    char* key;
    char* value;
    int len;

    doc = fopen(filename, "r");

    if (doc == NULL)
    {
        printf("warning: could not open %s, no options parsed\n", filename);
        return 0;
    }

    str = (char*)calloc(1024, sizeof(char));
    strcpy(str, "");
    while (fgets(str, 1024, doc) != NULL)
    {
        //if first character is a % sign, skip this line.
        if (str[0] == '%')
            continue;

        //if first character is a blank, skip this line.
        if (str[0] == ' ')
            continue;

        //if last character of line is newline, remove it
        len = strlen(str); 
        while (len > 0)
        {
            if (str[len - 1] == 10)
                str[len - 1] = '\0';
            else if (str[len - 1] == 13)
                str[len - 1] = '\0';
            else
                break;
            len = strlen(str);
        }

        //if line is now blank, skip it.
        if (strlen(str) == 0)
            continue;

        //read keyword by looking for an equal sign
        key = strtok(str, "=");

        if (key == NULL)
        {
            // no longer insist on having an argument
            key = str;
            value = NULL;
        }
        else
        {
            //read value
            value = strtok((char*)0, "=");
        }

        //apply the option... same routine command line options use
        applyOpt(key, value, options);
    }

    fclose(doc);
    free(str);

    return 1;
}
