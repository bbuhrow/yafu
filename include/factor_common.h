#pragma once

#include <stdint.h>
#include "gmp.h"

// these are similar to things msieve defines.  Differences:
// the factor type contains more info about how the factor
// was found and the list is not fixed size.

typedef struct
{
    mpz_t factor;
    int count;
    int type;
    // new parameters to support returning factors
    // from avx-ecm.  In general I think it makes sense
    // to have this structure contain information not
    // only about the factor itself but how it was found.
    int method;             // factorization method used
    uint32_t curve_num;     // curve found
    uint64_t sigma;         // sigma value
    int tid;                // thread found
    int vid;                // vector position found
} yfactor_t;

typedef struct
{
    yfactor_t* factors;
    int num_factors;
    int alloc_factors;

    // the primitude of factors can be proven by aprcl: these
    // determine when and what is displayed while proving.
    int aprcl_prove_cutoff;
    int aprcl_display_cutoff;
} yfactor_list_t;



/*--------------DECLARATIONS FOR MANAGING FACTORS FOUND -----------------*/

void add_to_factor_list(yfactor_list_t *factors, mpz_t n, 
    int VFLAG, int NUM_WITNESSES);
void print_factors(yfactor_list_t* fobj);
void clear_factor_list(yfactor_list_t* fobj);
void delete_from_factor_list(yfactor_list_t* fobj, mpz_t n);

