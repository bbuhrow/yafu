/*
 * RNS Montgomery modular multiplication, Yao regime, batch AVX2 + FMA
 * kernels: rnsyb_*_avx2, for builds without AVX-512.  Same 8-lane batch
 * layout as the AVX-512 kernel, processed as two 4-lane halves; FP-emulated
 * 52-bit products with MXCSR set to round toward zero for each call (build
 * with -frounding-math so the compiler keeps the FP ops between the MXCSR
 * writes).  The kernel body is rns_yaob_kern.h.
 */

#include "rns_yao.h"

#ifdef RNS_HAVE_AVX2
// the MXCSR switch: keep the FP ops in order (gcc/clang: -frounding-math)
#if defined(_MSC_VER) && !defined(__clang__)
#pragma fenv_access(on)
#elif defined(__clang__)
#pragma STDC FENV_ACCESS ON
#endif
#define RNS_VW 256
#define RNSYB_NAME(f) f##_avx2
#define RNSYB_HALVES 2
#include "rns_yaob_kern.h"
#endif
