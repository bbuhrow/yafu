#pragma once

#include "factor.h"
#include "gmp.h"

extern void test_dlp_composites();
extern void factor(fact_obj_t* fobj);

extern int factor_tiny(mpz_t in, mpz_t* out, 
	uint64_t* primes, uint64_t nump, uint64_t* prng);
