/*
tinyecm_edw: tinyecm (n < 2^128; x8: n < 2^104) with stage 1 on twisted
Edwards curves.  Same interfaces, B1 levels, curve counts and stage-2
tables as tinyecm.h's getfactor_tecm / getfactor_tecm_x8 /
getfactor_tecm_x8_list; the curves, stage 1 and the curve build are those
of microecm_edw (see arith/micro_tiny_edwards_plan.md).

This file is dual-licensed like tinyecm.h (FreeBSD or MPL 2.0).
*/

#ifndef TINYECM_EDW_H_INCLUDED
#define TINYECM_EDW_H_INCLUDED

#include <stdint.h>
#include "gmp.h"

#ifdef __cplusplus
extern "C" {
#endif

// As getfactor_tecm: returns 1 and a factor in f, or 0.  target_bits = the
// expected factor size (0: unknown, run the small-factor pre-passes first).
int getfactor_tecm_edw(mpz_t n, mpz_t f, int target_bits, uint64_t *pran);

// As getfactor_tecm_x8: one n < 2^104, 8 curves at a time (AVX-512; else the
// scalar version).
int getfactor_tecm_edw_x8(mpz_t n, mpz_t f, int target_bits, uint64_t *pran);

// As getfactor_tecm_x8_list: q64[2i], q64[2i+1] = low / high 52 bits of
// input i (< 2^104); f64[2i], f64[2i+1] = low / high 64 bits of its factor
// (1 if none).
void getfactor_tecm_edw_x8_list(uint64_t *q64, uint64_t *f64, int target_bits,
    uint32_t num_in, uint64_t *pran);

// test hooks (curve ids as in microecm_edw.h: GKL m >= 1, or
// TECM_EDW_TABLE | i for table curve i)
#define TECM_EDW_TABLE 0x80000000u
// one curve, scalar; returns 1 and the factor in f if found, *stage = 0 / 1 / 2
int tecm_edw_one_curve(mpz_t n, uint32_t id, uint32_t B1, uint32_t B2, mpz_t f, int *stage);
// the Montgomery (X : Z) after stage 1, out of Montgomery form
void tecm_edw_stage1_point(mpz_t n, uint32_t id, uint32_t B1, mpz_t X, mpz_t Z);
// x8 list at a given B1 and curve count (B2 = 25 B1); q64 / f64 as in
// getfactor_tecm_edw_x8_list; returns the number of curves run (0 without
// USE_AVX512F)
uint64_t tecm_edw_x8_run(uint64_t *q64, uint64_t *f64, uint32_t num, uint32_t B1, uint32_t curves);
// x8 on one n, curves rounded up to a multiple of 8; returns 1 and the factor
// in f if found, and the number of curve rounds run in *rounds
int tecm_edw_x8_one(mpz_t n, mpz_t f, uint32_t B1, uint32_t curves, int *rounds);
// Parameter set (process-wide, for experiments): 1 = legacy levels, 2 (the
// default) = the re-tuned levels of edw_tiny_v2.h (same level choice).  All entry points and
// hooks follow it.
void tecm_edw_set_params(int v);

#ifdef __cplusplus
}
#endif

#endif
