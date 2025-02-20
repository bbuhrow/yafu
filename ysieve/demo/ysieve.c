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

/*

*/

#if defined(WIN32)

//#include <windows.h>
//#include <process.h>

#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include "calc.h"
#include "cmdOptions.h"
#include "soe.h"
#include "ytools.h"
#include "soe_impl.h"


int main(int argc, char** argv)
{
    options_t* options;
    uint64_t start;
    uint64_t stop;
    uint64_t *primes;
    int count;
    uint64_t num_found;
    int haveFile;
    // timing
    double t;
    struct timeval tstart, tstop;
    soe_staticdata_t* sdata;
    char* startStr;
    char* stopStr;

    options = initOpt();
    processOpts(argc, argv, options);


    calc_init();
    startStr = process_expression(options->startExpression, NULL, 1, 0);
    stopStr = process_expression(options->stopExpression, NULL, 1, 0);
    calc_finalize();

    if (strlen(options->outFile) > 0)
    {
        haveFile = 1;
    }
    else
    {
        haveFile = 0;
    }

    sdata = soe_init(options->verbosity, options->threads, options->blocksize);
    sdata->witnesses = options->num_witnesses;
    sdata->userclasses = options->numclasses;
    sdata->analysis = options->analysis_type;
    sdata->gapmin = options->gapmin;

    if (options->outScreen || haveFile || (options->num_witnesses > 0))
    {
        count = 0;
        // these options require calculating primes
        if (sdata->analysis < 1)
            sdata->analysis = 1;
    }
    else
    {
        count = 1;
    }

    gettimeofday(&tstart, NULL);
    
    if (options->sieve_primes_limit > 0)
    {
        mpz_t low, high;
        mpz_init(low);
        mpz_init(high);

        mpz_set_str(low, startStr, 10);
        mpz_set_str(high, stopStr, 10);

        primes = sieve_to_depth(sdata, low, high, count, options->num_witnesses,
            options->sieve_primes_limit, &num_found,
            haveFile, options->outScreen);

        mpz_sub(high, high, low);
        mpz_sub_ui(high, high, num_found);

        gmp_printf("Removed %Zd composites\n", high);
        if (sdata->analysis == 2)
        {
            if (sdata->witnesses > 0)
            {
                printf("Probable twins  : %" PRIu64 "\n", num_found);
            }
            else
            {
                printf("Candidate twins : %" PRIu64 "\n", num_found);
            }
            
        }
        else
        {
            if (sdata->witnesses > 0)
            {
                printf("Probable primes : %" PRIu64 "\n", num_found);
            }
            else
            {
                printf("Candidate primes: %" PRIu64 "\n", num_found);
            }
            
        }
        gettimeofday(&tstop, NULL);
        t = ytools_difftime(&tstart, &tstop);
        printf("Elapsed time    : %1.6f seconds\n", t);

        mpz_clear(low);
        mpz_clear(high);
    }
    else
    {
        sscanf(startStr, "%" PRIu64 "", &start);
        sscanf(stopStr, "%" PRIu64 "", &stop);

        printf("starting sieve on bounds %" PRIu64 " : %" PRIu64 "\n", start, stop);

        primes = soe_wrapper(sdata, start, stop, count, &num_found,
            haveFile, options->outScreen);

        if (sdata->analysis == 2)
        {
            printf("Num twins found : %" PRIu64 "\n", num_found);
        }
        else
        {
            printf("Num primes found: %" PRIu64 "\n", num_found);
        }
        gettimeofday(&tstop, NULL);
        t = ytools_difftime(&tstart, &tstop);
        printf("Elapsed time    : %1.6f seconds\n", t);
    }

    soe_finalize(sdata);
    free(startStr);
    free(stopStr);
    if (primes != NULL)
    {
        align_free(primes);
    }

    return 0;
}
