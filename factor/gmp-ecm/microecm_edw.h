/*
microecm_edw: microecm (n < 2^64) with stage 1 on twisted Edwards curves.

Same interface and dispatch (B1 levels, curve counts, B2) as microecm.h's
getfactor_uecm.  Differences (see arith/micro_tiny_edwards_plan.md):
  - curves: Gelin-Kleinjung-Lenstra family (iv), a = -1, built mod n from
    [m]Q on E': y^2 = x^3 - x^2 - 9x + 9, Q = (5, 8), without inversions
    (consecutive m, random start per input from *pran);
  - stage 1: extended coordinates on precomputed signed-digit chains for
    the odd part of s, then the hand-off to the Montgomery image and the
    power of 2 as x-only doublings;
  - stage 2: microecm's x-only pairing stage 2, with projective A24.

This file is dual-licensed like microecm.h (FreeBSD or MPL 2.0).
*/

#ifndef MICROECM_EDW_H_INCLUDED
#define MICROECM_EDW_H_INCLUDED

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// Returns 1 if no factor was found, else a nontrivial factor of q64.
// is_arbitrary and *pran as for getfactor_uecm.
uint64_t getfactor_uecm_edw(uint64_t q64, int is_arbitrary, uint64_t *pran);

// ECM on a list of inputs of about the same size, 8 at a time with AVX-512
// (n < 2^52; otherwise, or without USE_AVX512F, one at a time with
// getfactor_uecm_edw).  f64[i] = a factor of q64[i], or 1.  B1 / curves by
// the size of q64[0], as getfactor_uecm_x8_list.  The x8 curve sequence per
// input is fixed (table curves, then GKL m = 3, 4, ...); *pran is used only
// by the scalar fallback.
void getfactor_uecm_edw_x8_list(uint64_t *q64, uint64_t *f64, uint32_t num_in, uint64_t *pran);

// test / research hooks
// A curve id >= 1 is GKL's [m]Q curve, m = id; UECM_EDW_TABLE | i is curve i
// of the small-constant table (edw_tiny_curves.h).
#define UECM_EDW_TABLE 0x80000000u
// One curve at a supported B1 (27, 47, 70, 85, 125, 165, 205) with
// microecm's B2 for that B1 (or B2 <= B1: stage 1 only).  Returns the factor
// found (1 if none); *stage = 0 (none), 1 or 2.
uint64_t uecm_edw_one_curve(uint64_t n, uint32_t id, uint32_t B1, uint32_t B2, int *stage);
// The curve, out of Montgomery form: X, Y, Z, T of the start point
// (extended coordinates) and d = dn / dd.
void uecm_edw_curve_params(uint64_t n, uint32_t id, uint64_t out[6]);
// The Montgomery (X : Z) after stage 1, out of Montgomery form.
void uecm_edw_stage1_point(uint64_t n, uint32_t id, uint32_t B1, uint64_t out[2]);
// How many table curves getfactor_uecm_edw runs before the GKL sequence
// (default: all; process-wide, for experiments only).
void uecm_edw_set_table(int ntab);
// Built with -DUEDW_STATS: copies and clears the histogram of curves run per
// input (hist[0..max]; the last bin counts "or more"); returns its size - 1,
// or -1 without UEDW_STATS.
int uecm_edw_hist(uint64_t *hist, int max);
// x8 list at a given B1 and curve count (B2 = 25 B1); returns the number of
// curves run (0 without USE_AVX512F).
uint64_t uecm_edw_x8_run(uint64_t *q64, uint64_t *f64, uint32_t num, uint32_t B1, uint32_t curves);

#ifdef __cplusplus
}
#endif

#endif
