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
#ifndef _YECM_H_
#define _YECM_H_

#include "factor.h"
#include "ytools.h"
#include "arith.h"
#include <gmp_xface.h>

/* ============================ interface to microecm ============================ */
extern void init_uecm(uint64_t lcg);
extern uint64_t do_uecm(uint64_t q);
extern void microecm(uint64_t n, uint64_t* f, uint32_t B1, uint32_t B2, uint32_t curves, int verbose);

/* ============================ interface to tinyecm ============================ */
extern void tinyecm(mpz_t n, mpz_t f, uint32_t B1, uint32_t B2, uint32_t curves,
    uint64_t* lcg_state, int verbose);

/* ============================ interface to gmpecm/pm1/pp1 ============================ */
extern void pollard_loop(fact_obj_t* fobj);
extern void williams_loop(fact_obj_t* fobj);
extern int ecm_loop(fact_obj_t* fobj);


/* ============================ interface to avxecm ============================ */
extern yfactor_t* vec_ecm_main(mpz_t N, uint32_t numcurves, uint64_t B1,
    uint64_t B2, int threads, int* numfactors, int verbose,
    int save_b1, uint32_t* curves_run);

#endif // #ifndef _YECM_H_
