#ifndef SIQS_STRATEGY_H
#define SIQS_STRATEGY_H

// SIQS large-prime sieve strategy values, set with -siqsMid, -siqsLarge,
// -siqsKernel and -siqsSparseKernel (command line or yafu.ini) and
// resolved per factorization in siqs_static_init.  0 = automatic.
#define SIQS_STRAT_AUTO     0
#define SIQS_MID_BUCKET     1   // [med_B, large_B) bucket sieve
#define SIQS_MID_SPARSE     2   // [med_B, large_B) per-block sparse sieve
#define SIQS_LARGE_BUCKET   1   // [large_B, B) bucket sieve
#define SIQS_LARGE_DIRECT   2   // [large_B, B) direct interval sieve
#define SIQS_K_BITLOOP      1   // compare mask, trail_zcnt/reset_lsb loop
#define SIQS_K_SINK         2   // blend misses to a sink stripe, write all lanes
#define SIQS_K_CREG         3   // compress in register + unaligned store, batched
#define SIQS_K_CMEM         4   // compress-store to memory, batched

#endif
