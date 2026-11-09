/*
 * AVX-ECM stage 2, FFT continuation (avxecm_stage2fft.c).
 * Design: arith/avxecm_stage2_fft_design.md.
 */

#ifndef AVXECM_STAGE2FFT_H
#define AVXECM_STAGE2FFT_H

#include <stdint.h>
#include <gmp.h>

typedef struct avxecm_s2fft_s avxecm_s2fft_t;

#define AVXECM_BS_AUTO      64      // choose the Brent-Suyama polynomial
#define AVXECM_BS_DEFAULT   65      // option not given: AVXECM_BS_AUTO with
                                    // the automatic B2, else off

// Shared setup for one run: N = the arithmetic modulus (the backend's
// mdata->n), stage 2 range (B1, B2], per-thread memory limit in MB.
// Picks dF, D and the block count; the covered B2 (>= B2) is
// avxecm_s2fft_B2().  NULL if the FFT stage 2 is unavailable (no vector
// kernels, B2 <= B1, or B1 too small for the smallest D).
// auto_b2: ignore B2 and choose the most efficient one for the memory
// limit (design doc section 11.3: the least predicted time per factor at
// the factor size B1 is usually chosen for).  bs: Brent-Suyama polynomial
// (0 off, S > 1 x^S, S < -1 Dickson D_|S|(x, -1); design doc section 14;
// AVXECM_BS_AUTO: the least predicted time per factor among plain, D_3,
// D_6 and D_12).
avxecm_s2fft_t *avxecm_s2fft_setup(mpz_t N, uint64_t B1, uint64_t B2,
    int threads, double mem_mb, int auto_b2, int bs, int verbose);
uint64_t avxecm_s2fft_B2(avxecm_s2fft_t *s);
void avxecm_s2fft_free(avxecm_s2fft_t *s);

#endif
