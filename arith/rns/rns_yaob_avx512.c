/*
 * RNS Montgomery modular multiplication, Yao regime, batch AVX-512 kernels
 * (IFMA, or FP-emulated products without it): rnsyb_*_avx512.  The kernel
 * body is rns_yaob_kern.h.
 */

#include "rns_yao.h"

#ifdef RNS_HAVE_AVX512
#define RNS_VW 512
#define RNSYB_NAME(f) f##_avx512
#define RNSYB_HALVES 1
#include "rns_yaob_kern.h"
#endif
