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

#ifndef _SIQS_H_
#define _SIQS_H_

#include "factor.h"
#include <stdint.h>
#include "ytools.h"


extern void SIQS(fact_obj_t* qsobj);
extern void siqsbench(fact_obj_t* qsobj, info_t *comp_info, uint32_t digit_limit);

// SIQS parameter tables: the competing parameter sets exported by
// siqs_tune.py (siqs_candidates.txt, several rows per digit size) and the
// per-machine table siqstune() picks from them (siqs_params.txt, one row per
// digit size).  Both are whitespace-separated text whose columns are named
// by a "# digits bits rank ..." header line; "-" marks a missing value.
#define SIQSP_B         0x0001
#define SIQSP_M         0x0002
#define SIQSP_NB        0x0004
#define SIQSP_TF        0x0008
#define SIQSP_TFSM      0x0010
#define SIQSP_MFBD      0x0020
#define SIQSP_SPV       0x0040
#define SIQSP_MID       0x0080
#define SIQSP_LARGE     0x0100
#define SIQSP_KERNEL    0x0200
#define SIQSP_PARAMS    0x03ff

typedef struct
{
    int digits;
    int bits;
    int rank;               // 0 = yafu's built-in defaults
    double ratio;           // time relative to the defaults (where measured)
    uint32_t has;           // SIQSP_* bits of the parameters present
    uint32_t B, M, NB, TF, TFSm, SPV;
    double MFBD;
    uint32_t mid, large, kernel;    // SIQS_MID_*, SIQS_LARGE_*, SIQS_K_*
    // measurements, < 0 when missing
    double partials;        // partial relations collected at completion
    double partials_cv;
    double fulls;
    double needed;
    double t_sieve;
    double t_post;          // filtering + linear algebra + sqrt
    double est_time;        // siqstune(): estimated total time
    char n[256];            // the input the measurements were made on ("" = none)
} siqs_param_row_t;

extern int siqs_read_param_table(const char* fname, siqs_param_row_t** rows);
extern int siqs_write_param_table(const char* fname, siqs_param_row_t* rows,
    int num_rows, const char* comment);
extern int siqs_lookup_tuned(const char* fname, int bits, siqs_param_row_t* out);
extern void siqs_format_params(siqs_param_row_t* row, char* buf, int len);

#endif /* _SIQS_H_ */


