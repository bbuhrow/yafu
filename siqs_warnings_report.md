# SIQS gcc warning cleanup — report (2026-11-06, branch L2sieve, uncommitted)

## Result
gcc 11.4 (WSL Ubuntu 22.04) with the Makefile's `-Wall -Wconversion -Wshadow -Wstrict-prototypes`:
**0 warnings** (baseline: ~1650 unique / ~2100 sites) in all six configs:

| config | flags |
|---|---|
| generic (SSE2) | all USE_* = 0 |
| sse41 | USE_SSE41=1 |
| avx2 | USE_AVX2=1 USE_BMI2=1 |
| avx512 | USE_AVX512=1 |
| ifma | USE_AVX512IFMA=1 |
| neon | `aarch64-linux-gnu-gcc` cross, IS_ARM=1 |

Files covered: every file in `YAFU_SIQS_SRCS` for each config, plus `factor/qs/qs_demo/{siqs_demo,cmdOptions,calc}.c` and `factor/qs/msieve/{gf2,sqrt}.c` (in the VS `siqs` project), and every header they pull in.

## How the fixes were made
* **Conversions (~90%)**: explicit casts to exactly the type gcc was converting to implicitly (assignment/return/argument), or to the signed operand of a binary op whose other operand already forces that type. Same values; in most files the generated code is byte-identical (see Verification).
* **Unused variables**: deleted only when the name appears nowhere else in the function under any `#ifdef` and in no macro body; otherwise marked `UNUSED_VAR` (the existing ytools.h macro, which MSVC also handles). No functions or static tables were deleted (`now_sec` → `UNUSED_FUNC`, `mult_list_orig`/`mult_list_many_even` → `UNUSED_VAR`).
* **Shadowing**: renamed the inner variable (`parity`→`c_parity` in tdiv.c, `s`→`str` in siqs_test.c, `prime`→`fbprime` in msieve/sqrt.c, outer `i`→`num_apoly` in filter.c `restart_siqs`) or dropped an inner redeclaration of a same-typed outer variable whose value is always re-assigned first (`int k` in the KNL/AVX2 slice loops, `int i` blocks in tdiv.c, `uint64_t q64` in tdiv_small.c).
* **Macro redefinitions**: `#undef` before the later definition (the later one already won), or made the duplicate definitions token-identical.

## Build system
* `Makefile`: GMP and GMP-ECM include dirs now passed with `-isystem` instead of `-I`, so warnings inside gmp.h/ecm.h inline code no longer show up in our builds.

## Shared headers touched (used outside SIQS too)
* `include/portable_carry.h`: the aarch64 carry/borrow helpers are now `extern __inline __attribute__((gnu_inline, always_inline))` (via `YAFU_PC_INLINE`, like gcc's own x86 intrinsic headers) instead of `static inline`. This fixes "static but used in inline function which is not static" from monty.h's `__inline mulredc*`. Only reached in non-MSVC branches.
* `include/monty.h`: `_umul128` on the portable path uses the same `YAFU_PC_INLINE`; AVX-512 non-IFMA 104-bit routines got `UNUSED_VAR` on unused locals plus a `#pragma GCC diagnostic` push/pop around the intentional shadowing of `dbias`/`vbias*` (the mul52 macros reference those names). Guarded by `__GNUC__`.
* `include/gmp_xface.h`: one cast.

## Flagged items — please review
1. **SIQS.c ~522** `for (i = 0; i < fobj->THREADS; i++);` had a stray `;`, so the block ran once. I removed the `;`, which doesn't change behavior because the body sets `t[tid].dconf->batch_run_override = 1` (the same thread every pass). **It probably should be `t[i]`.** I didn't change that.
2. **SIQS.c TLP status printf (~2936)**: a misplaced `)` left the `r/sec` argument out of the call, so `%1.2f` printed garbage (undefined behavior). Fixed. Only the printed output changes.
3. **printf formats**: `%u`/`%d` used with `uint64_t`/`size_t` arguments changed to `PRIu64` or given a matching cast (SIQS.c: "checked %…", "trial division touched %…", byte-count prints). Output only.
4. **String copies** (behavior only changes where the old code left a string unterminated):
   * SIQS.c: `strncpy(savefile_name, siqs_savefile, 80)` became a ≤80-char copy that is always NUL-terminated.
   * qs_demo/cmdOptions.c: `optbuf`/`argbuf` copies now use `N-1` and terminate explicitly.
   * qs_demo/calc.c: `strncpy(str, ptr, strlen(ptr)-1)` became `memcpy` (identical).
   * cmdOptions.c: `"% "PRIu64` typo fixed to `"%" PRIu64`.
5. **fgets** (SIQS.c resume check, filter.c `restart_siqs`): I kept the behavior exactly. The result is only tested so that `warn_unused_result` is satisfied, and the empty body says so. Separate pre-existing issue: the `while (!feof(data)) { fgets(...); ... }` loops in `restart_siqs` process the last line twice when the final `fgets` fails. Not changed.
6. **tdiv_large.c**: `DIVIDE_ONE_PRIME(i);` was invoking an object-like macro, which left a no-op `(i);` behind. Changed to `DIVIDE_ONE_PRIME;` (identical code).
7. **Pre-existing siqs_demo crash** (seen in both HEAD and new builds): `free(): invalid next size (fast)` when stdin is not a tty (batch mode) and the number on the command line is 40 digits (`1473679397014232648814484694383531642231`). Not investigated.
8. **tdiv_med_32k_avx2.c holds your uncommitted WIP.** My edits there are 16 `(uint16_t)block_loc` casts, 2 `(uint16_t)i` casts and `#undef TDIV_MED_CLEAN`. With and without them the file compiles to identical code. FYI: the WIP itself changes the relation yield compared with HEAD on a C49 (`9264451099144225512584183028022986703345379388409`, seed 42, IFMA: 7943 → 7920 relations read). It is deterministic and the factorization still succeeds.

## Verification
* **Disassembly diff, HEAD vs. new**, per file, configs generic/avx2/ifma/neon, built `-O2` without `-g`. Most objects are instruction-identical. In the ones that differ I checked the differences: the intentional fixes above, or register/scheduling changes where a cast sits on an operand instead of on the result. Those give the same values mod 2^32 except where the old code relied on signed overflow, which is undefined anyway. Two spots that only differed in scheduling (cofactorize `cutoff1`, new_poly `abs(a1-p)`) were rewritten so they compile identically again.
* **A/B factoring with siqs_demo** (seed 42, 1 thread, 7 inputs from C39 to C76, run under a pty):
  * generic build: output identical to HEAD (relation counts, cycles, matrix, factors).
  * IFMA build: identical except the C49 difference from item 8. Building HEAD plus only the current tdiv_med_32k_avx2.c reproduces it.
  * All factors were found.
* Full `make yafu` (IFMA) succeeds. The other users of `portable_carry.h` (microecm.c, micropm1.c, arith.c, factor_common.c) compile cleanly for aarch64.

## Not verified
* **MSVC**: no build available. I reviewed the new constructs only by reading: every `__attribute__`, `#pragma GCC` and `YAFU_PC_INLINE` sits in a non-MSVC branch, and `UNUSED_VAR`/`UNUSED_FUNC`/`PRIu64` are already defined for MSVC. The explicit casts should also clear many MSVC C4244/C4267 warnings.
* clang builds; running the NEON build on real ARM hardware.
