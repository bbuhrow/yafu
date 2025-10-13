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

#ifndef _NFS_H_
#define _NFS_H_


#include "factor.h"
#include "gnfs-yafu.h"

enum special_q_e
{
    NEITHER_SPQ,
    RATIONAL_SPQ,
    ALGEBRAIC_SPQ
};

typedef struct
{
    mpz_poly_t rat; // linear (usually)
    mpz_poly_t alg;
    double skew;
    double murphy; // murphy e score
    double size;
    double alpha;
    int rroots;
    mpz_t m; // common root mod n
    enum special_q_e side;
} mpz_polys_t;

#define NUM_SNFS_POLYS 3
#define MAX_SNFS_BITS 1024

enum snfs_form_e
{
    SNFS_NONE,
    SNFS_BRENT,
    SNFS_H_CUNNINGHAM,
    SNFS_XYYXF,
    SNFS_DIRECT,
    SNFS_LUCAS
};

typedef struct
{
    // input integer
    mpz_t n; // the cofactor, what goes in the job file
    mpz_t primitive;
    // algebraic representation of the snfs form:
    // n divides c1*b1^e1 + c2*b2^e2
    mpz_t base1;
    mpz_t base2;
    int exp1;
    int exp2;
    int coeff1;
    int coeff2;
    // type of form
    enum snfs_form_e form_type;


    mpz_polys_t* poly;
    mpz_t c[MAX_POLY_DEGREE + 1]; // scratch space -- converted to mpz_poly_t
                      // in check_poly()

    // other useful parameters
    double difficulty;
    double sdifficulty;
    double anorm;
    double rnorm;
    int rank;
    int valid;
    int siever;
} snfs_t;

/* ============================ interface to nfs factoring ============================ */
extern void nfs(fact_obj_t* fobj);
extern snfs_t* snfs_find_form(fact_obj_t* fobj, mpz_t n);
extern void remove_algebraic_factors(fact_obj_t* fobj, 
    snfs_t* poly, uint64_t* primes, uint64_t num_p, int VFLAG);
extern void remove_algebraic_factors_lucas(fact_obj_t* fobj,
    snfs_t* poly, uint64_t* primes, uint64_t num_p, int VFLAG);

#endif
