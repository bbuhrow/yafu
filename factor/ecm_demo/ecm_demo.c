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
#include "factor.h"
#include "yafu_ecm.h"

int main(int argc, char** argv)
{
    options_t* options;
    int count;
    int haveFile;
    // timing
    double t;
    struct timeval tstart, tstop;
    char* inputN;
    char* B1;
    char* B2;
    uint64_t b1;
    mpz_t b2;
    mpz_t N;
    fact_obj_t fobj;

    options = initOpt();
    processOpts(argc, argv, options);

    calc_init();
    inputN = process_expression(options->inputN, NULL, 1, 0);
    B1 = process_expression(options->B1, NULL, 1, 0);
    B2 = process_expression(options->B2, NULL, 1, 0);
    calc_finalize();

    mpz_set_str(N, inputN, 10);
    mpz_set_str(b2, B2, 10);
    b1 = strtoull(B1, NULL, 10);

    init_factobj(&fobj);
    new_factorization(&fobj, N);

    // customize for this method
    mpz_set(fobj.ecm_obj.gmp_n, N);
    fobj.ecm_obj.num_curves = options->numcurves;
    fobj.ecm_obj.B1 = b1;
    fobj.ecm_obj.prefer_avxecm_stg2 = options->prefer_gmpecm_stg2;
    fobj.ecm_obj.prefer_gmpecm = options->prefer_gmpecm;
    fobj.ecm_obj.use_gpuecm = options->gpu;
    fobj.ecm_obj.use_gpudev = options->gpudev;
    fobj.ecm_obj.save_b1 = options->save_b1;
    fobj.THREADS = options->threads;
    fobj.VFLAG = options->verbosity;
    fobj.ecm_obj.ecm_ext_xover = options->extecm_xover;
    strncpy(fobj.ecm_obj.ecm_path, options->extecm_path, 1024);
    
    if (options->sigma > 0)
    {
        fobj.ecm_obj.sigma = options->sigma;
    }

    if (mpz_cmp_ui(b2, 0) > 0)
    {
        fobj.ecm_obj.B2 = mpz_get_ui(b2);
    }
    else
    {
        fobj.ecm_obj.B2 = 100 * b1;
    }

    ecm_loop(&fobj);

    print_factors(&fobj);
   
    free(inputN);
    free(B1);
    free(B2);


    return 0;
}
