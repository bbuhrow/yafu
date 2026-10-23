# Porting lasieve5 mpqs / mpqs3 into yafu: plan

*Drafted 2026-10-23 on branch `L2sieve` (head `d52b9b7`). Status of each
phase is in section 0.*

## 0. Status

### Phase 1: pure-C import (done)

- `factor/kmpqs/`: `kmpqs1.c` (from `mpqs.c`), `kmpqs3.c` (from
  `mpqs3.c`), `kmpqs.c` (size dispatch), `kmpqs.h` (API),
  `kmpqs_common.h` (replaces lasieve's `if.h`, `gmp-aux.h`, `psp.c`),
  `kmpqs_bench.c` (driver).
- Generated mechanically: every `ASM_*`, `HAVE_XMM_MUL`, `A64_STYLE_TD`,
  `*_ZEIT` and `TOTAL_STAT` branch was resolved to its C fallback, the
  param tables were copied in from `asm/mpqs-config.h`, lasieve types were
  mapped to `<stdint.h>`, and all globals were made `static`. The
  `MPQS_STAT`/`MPQS3_STAT` counters are kept, for comparison against
  lasieve.
- lasieve's `Schlendrian()`/`complain()` used to `exit()`. They now longjmp
  back to the entry point, which returns `KMPQS_ERR_INTERNAL` (-100).
  Warnings print only if `kmpqs_verbose`.
- `psp()` (asm) is replaced by yafu's `MR_2sprp_64x1` / `MR_2sprp_128x1`,
  with `mpz_probab_prime_p` above 128 bits.
- API: `kmpqs1_factor`, `kmpqs3_factor`, `kmpqs_factor` (dispatch at 96
  bits, as in lasieve's strategy). These have the lasieve signature, but
  under new names, so they cannot collide with lasieve objects.
- Builds: Makefile `COMMON_SRCS`, both VS2022 `factor_common` projects,
  and `make kmpqs-bench`. Checked with WSL gcc 11, MinGW gcc (native
  Win64 run: same counters as Linux, no crash) and clang-cl
  (`-fsyntax-only`).
- Validation: `tools/kmpqs_compare.sh` builds lasieve as-is and with the
  asm switched off. Every counter (sieves, TD candidates/survivors, full,
  partial, combined relations) matches the no-asm lasieve exactly, up to
  100 bits.

**Bug found and fixed (lasieve's own, asm or not):** `mpqs3_convert` builds
the high word of a 96-bit value from `(op_dbl - (double)op_64) * 2^-64`. It
added 2^32 to a negative difference *before* rounding. Once `axb > 2^53`
(N above ~104 bits) the double rounding error can make that difference
slightly negative, and the high word came out `0xffffffff` instead of 0.
The relation was then corrupt, the matrix never yielded a split, and
mpqs3 returned -6 (matrix full) after seconds of sieving. In lasieve,
mpqs3 failed on 5/30 random 104-bit semiprimes, 28/30 at 108 bits and
30/30 at 110+. This is very likely the "lots of -6 errors with MFB=107"
noted at `MPQS3_GAUSS_MAX`. kmpqs3 now rounds first. 110- and 120-bit
inputs all succeed (about 2.3 ms and 4.6 ms, pure C). The same one-line
fix applies to `factor/lasieve5_64/mpqs3.c`, which the NFS siever uses.

Phase 1 timings (`tools/kmpqs_compare.sh`, WSL, i5-11600K; kmpqs built
with the Makefile's flags):

| bits | lasieve asm | kmpqs (C) | ratio |
|---|---|---|---|
| 50 | 0.063 s | 0.084 s | 1.33 |
| 64 | 0.137 s | 0.178 s | 1.30 |
| 80 | 0.356 s | 0.512 s | 1.44 |
| 90 | 0.673 s | 0.921 s | 1.37 |
| 96 (mpqs1) | 0.999 s | 1.435 s | 1.44 |
| 100 (mpqs3) | 1.351 s | 1.999 s | 1.48 |
| 110 / 120 (500 each) | fails | 1.127 / 2.321 s | |

(2000 inputs per size up to 100 bits.)

### Phase 2: context struct / thread safety (done)

- A script (`ctxify.py`, in the session notes, not committed) moved every
  mutable file-scope variable into `mpqs_ctx_t` / `mpqs3_ctx_t` (112 and
  109 variables; `mpqs_X` -> `qs->X`, `mpqs_2B` -> `qs->twoB`). It also
  added a `qs` parameter to every function that touches state, calls a
  function that does, or can fail. Pure helpers (`invert`, `sqrt`,
  `powmod_mm`, `inv_64` in kmpqs1, `convert`/`mod3`/`dbl64mod*` in kmpqs3)
  are unchanged. Only constant tables remain file-scope.
- The context also holds the prime and multiplier tables lasieve built
  once in `mpqs_init()`. Init runs once per context, not per call.
  Context size: about 0.3 MB (kmpqs1) and 1.3 MB (kmpqs3), plus
  allocations.
- API (`kmpqs.h`): `kmpqs_ctx_new/free` and `kmpqs[1|3]_factor_ctx`
  (explicit context; sub-contexts are created on first use).
  `kmpqs[1|3]_factor` use a thread-local context (`__thread` /
  `__declspec(thread)`); `kmpqs_thread_cleanup()` frees it.
  `KMPQS_ERR_INTERNAL` is now public.
- Other fixes: `C=(B^2-kN)/A` in kmpqs1 is now computed in `uint64_t`
  (lasieve relied on signed-overflow wraparound, which UBSan flags).
  kmpqs1 no longer allocates a `gauss_row[0]` that the C path overwrote
  (leaked in lasieve). kmpqs3's file-static `alog` is now local. Dead
  `powmod0`/`sqrt0` removed.
- Makefile: `KMPQS_WFLAGS` silences lasieve's `-Wconversion`/unused-variable
  noise for `factor/kmpqs` objects.
- Validation:
  - counters still match lasieve no-asm exactly (50..100 bits)
  - `kmpqs_bench -threads T` factors the inputs on T threads and compares
    every result with a single-threaded run: 0 differences on 3300 mixed
    50..130-bit inputs, with 4, 12 and 24 threads, on Linux and MinGW
  - ThreadSanitizer (4 threads) and ASan+UBSan: clean
  - speed unchanged within noise (±5%)

## 1. Goal

Move Kleinjung's small-number MPQS (`factor/lasieve5_64/mpqs.c` for
N <= ~96 bits, `mpqs3.c` for ~90..~150 bits) into yafu proper:

- Replace every m4/GAS `.asm` kernel with C intrinsics.
- Widen the SIMD kernels to AVX2 and AVX-512 (and later NEON, see the
  SIQS NEON plan) where it pays.
- Make it callable from yafu's threaded code (it is currently
  single-threaded global state).
- Use it wherever yafu needs to split ~50..~130-bit composites.

## 2. What exists today

### 2.1 lasieve5_64 sources

| file | role | lines |
|---|---|---|
| `mpqs.c` | MPQS, N <= 96 bits (`TINY` params, 14 size classes 44..96 bits) | ~2850 |
| `mpqs3.c` | MPQS, 3x32-bit words, N ~92..~150 bits (`TINY3`, 15 classes) | ~3150 |
| `asm/mpqs-config.h` | feature switches (`ASM_MPQS_*`, `HAVE_XMM_MUL`) + param tables | |
| `asm/mpqs3arith.c` | `mpqs3_mod3` (C) | |

The asm kernels (m4-preprocessed AT&T GAS, **SysV ABI only**):

| kernel | file | SIMD today | globals it touches |
|---|---|---|---|
| sieve (`asm_sieve`) | `mpqs_sieve.asm` | none (scalar, unrolled) | FB, FB_log, FB_start, sievearray, sievebegin, sievelen |
| sieve init | `mpqs_sieveinit.asm` | MMX/SSE2 | tiny prod / tiny array |
| eval (survivor scan) | `mpqs_eval.asm` | SSE2 | none (args only) |
| trial division | `mpqs_td.asm` | SSE2 16-bit (`pmullw/pmulhw`) | ~15 FB / inverse tables |
| next polynomial | `mpqs_nextpol.asm` | SSE2 16-bit | FB_disp, FB_mm_inv, FB_np_p(x), FB_start, nFB |
| gauss | `mpqs_gauss.asm` | MMX/SSE2 | gauss_* matrix state |
| mpqs3 sieve / sievea | `mpqs3_sieve.asm` | none (scalar) | mpqs3 FB tables + some mpqs_* |
| mpqs3 td | `mpqs3_td.asm` | SSE2 16-bit | ~16 mpqs3 tables |
| mpqs3 next pol | `mpqs3_nextpol.asm` | SSE2 16-bit | mpqs3 FB tables |
| mpqs3 gauss (`asm_re_strip`) | `mpqs3_gauss.asm` | SSE2 | none (args only) |

**Every asm block in both `.c` files has a C fallback** (`#else`
branch). The pure-C build works (section 3).

### 2.2 How yafu uses it now

- `factor/batch_factor.c` (8 call sites), `factor/autofactor.c`
  (lasieve cofactor strategy, 4 sites) and `factor/gpu_cofactorization*.c`
  call `mpqs_factor` / `mpqs3_factor` under `#ifdef HAVE_LASIEVE_MPQS`.
- The main yafu Makefile never defines `HAVE_LASIEVE_MPQS`.
  `gpu_cofactorization_cl.c` defines it for non-MinGW gcc, with the
  comment "mpqs builds in msys2 but crashes right away when run". That
  crash is the SysV-ABI asm being called under the Win64 ABI. The
  intrinsics port removes it.
- yafu's own small QS options: `tinyqs` (`factor/tinyqs/tinySIQS.c`,
  `init_tinyqs`/`tinyqs` in `cofactorize.h`), `smallmpqs`
  (`factor/qs/smallmpqs.c`), and msieve's `mpqs/`.

### 2.3 Prior benchmarks (`../tinyfactor/results-i5-11600k.txt`)

| bits | ggnfs mpqs (asm) | ggnfs mpqs3 (asm) | yafu tinyqs |
|---|---|---|---|
| 64 | 64 us | | 207 us |
| 80 | 174 us | | 506 us |
| 90 | 307 us | 336 us | 817 us |
| 100 | | 612 us | 1474 us |
| 120 | | 2539 us | |

So lasieve mpqs is about 2.5-3x faster than tinyqs over this whole
range. That is the reason to bring it in.

## 3. Baseline measured this session

WSL Ubuntu 22.04, gcc 11.4, i5-11600K, `-Ofast -march=native`. Inputs:
the first 2000 lines of `tinyfactor/semiprimes_<b>bit.dat`. These are
test builds (`-DMPQS_STAT -DMPQS_ZEIT`), so they include timer overhead
(`zeitA`/`zeitB` take ~3-5% of samples).

| test | asm build | pure-C build | C / asm |
|---|---|---|---|
| mpqs 50-bit | 0.124 s | 0.129 s | 1.04 |
| mpqs 64-bit | 0.237 s | 0.293 s | 1.24 |
| mpqs 80-bit | 0.674 s | 0.861 s | 1.28 |
| mpqs 90-bit | 1.253 s | 1.636 s | 1.31 |
| mpqs3 90-bit | 1.237 s | 1.867 s | 1.51 |
| mpqs3 100-bit | 2.507 s | 3.640 s | 1.45 |

Both builds factor all 2000 inputs at every size.

Profile (perf, mapped to symbols with `nm`, non-PIE relink):

- asm build, mpqs 80-bit: `mpqs_factor0` (inlined C) 21%, sieve
  `loop4` 13%, td `loop2` 11%, `mpqs_jacobi` 6%, `mpqs_invert` 4%,
  `mpqs_SI_init` 4%, sieve-init `innerloop16` 2%, gauss 2%, eval 2%.
- asm build, mpqs3 100-bit: sieve `loop4` 21%, `mpqs3_factor0` 16%, td
  `loop2` 13%, sieve-init 6%, jacobi 3%, gauss 3-4%.
- pure-C, mpqs3 100-bit: `mpqs3_factor0` 40% (td/decompose inlined),
  `mpqs3_sieve` 33%, `next_pol` 9%, `evaluate` 5%.

Takeaways:
1. Sieve and TD are the two kernels that matter. Next-pol, eval and
   sieve-init come after them. Gauss is small at these sizes.
2. The hot sieve kernel is *scalar* asm. It needs good C, not wide
   SIMD.
3. Setup costs (jacobi, invert, `SI_init`) are ~15% at 80 bits. These
   are scalar C and can be improved independently.

The harness is in the session scratchpad (not committed): copies of
`lasieve5_64` built as `mpqs_base` (asm) and `mpqs_noasm` (all
`ASM_MPQS*`, `HAVE_XMM_MUL`, `A64_STYLE_TD` off), plus `prof.sh`. Phase 1
should commit a reproducible equivalent.

## 4. Design decisions (proposed)

1. **Location.** New directory `factor/kmpqs/` (Kleinjung MPQS), with
   `kmpqs.c` (<= 96 bit) and `kmpqs3.c` (3-word). Per-ISA kernel files
   follow the SIQS pattern: `kmpqs_kernels_sse2.c`, `_avx2.c`,
   `_avx512.c`, later `_neon.c`. lasieve5_64 keeps its own copy for now
   and switches over once the port is proven.
2. **State.** Move the file-scope statics into a context struct
   (`kmpqs_ctx_t`) passed to every routine. Kernels get pointers and
   lengths as arguments, never globals. This is a precondition for
   thread safety, and the asm is what stops it today. Allocate once per
   thread and reuse across calls (the ~1 us per-call setup matters at
   50 bits).
3. **API.** Keep the lasieve signature as a thin wrapper so the existing
   `HAVE_LASIEVE_MPQS` call sites keep working:
   `int mpqs_factor(mpz_t N, size_t max_bits, mpz_t **factors)` plus a
   reentrant `kmpqs_factor(kmpqs_ctx_t*, mpz_t N, size_t max_bits,
   mpz_t *factors)`, and a dispatcher that picks mpqs / mpqs3 by size.
4. **ISA selection.** Compile-time, using the existing yafu build flags
   (`USE_AVX2`, `USE_AVX512F`/`BW`) like the rest of the cofactor code.
   SSE2 is the floor on x86-64, and there is a portable scalar fallback
   for ARM until NEON. Runtime dispatch can come later if the SIQS
   multi-ISA machinery is reused.
5. **Validation.** At each step, compare the factors found and the
   counts (relations, polynomials, TD candidates) against the asm build
   on fixed input sets. The algorithms are deterministic, so counters
   should match exactly until a deliberate parameter change.
6. **License.** lasieve5 is GPL (Franke/Kleinjung). yafu's own sources
   are public domain. Compiling these files into yafu makes the yafu
   *binary* GPL-encumbered, as linking lasieve's library already does
   for builds with `HAVE_LASIEVE_MPQS`. **This needs the owner's
   decision.** Keep the GPL headers on the ported files either way.

## 5. Phases

**Phase 1: pure-C import, no behavior change.**
- Copy `mpqs.c`/`mpqs3.c` and the needed helpers into `factor/kmpqs/`.
  Delete the asm branches. Keep the C fallbacks. Strip `MPQS_ZEIT`/stat
  instrumentation or put it behind one flag.
- Remove dependencies on lasieve's `if.h`, `siever-config.h`,
  `gmp-aux.h`, `zeit.h`. Use yafu types and `ytools`.
- Build it in the yafu Makefile and the VS2022 projects. Check that it
  runs on Win64 (MSVC/clang-cl and MinGW), Linux and WSL.
- Commit the benchmark/validation harness: a `kmpqs_test` driver in
  `test/` that reads semiprime files and prints time plus counters, and
  a script that runs it against a reference asm build.
- Gate: identical counters to the noasm baseline, and no Win64 crash.

**Phase 2: context struct / thread safety.**
- Mechanical move of the globals into `kmpqs_ctx_t`, with kernels
  taking explicit arguments.
- Gate: same counters, plus a multithreaded stress test (N threads
  factoring disjoint inputs, all results correct).

**Phase 3: kernels, in profile order.** Each kernel gets a scalar
reference and SSE2/AVX2/AVX-512 versions, plus a unit test against the
scalar version.
1. Sieve (`mpqs_sieve`, `mpqs3_sieve`/`sievea`): scalar, so the goal is
   to match or beat the asm in C (unrolling, two roots per prime,
   small-prime handling). It is worth trying AVX-512 scatter-free tricks
   for the smallest primes only.
2. TD (`mpqs_td`, `mpqs3_td`): 16-bit `mullo/mulhi` root tests. Going
   from SSE2 to AVX2 to AVX-512BW is 8 to 16 to 32 primes per
   instruction, with mask compares and compress.
3. Next-pol (`mpqs_nextpol`, `mpqs3_nextpol`): 16-bit root updates,
   same widening.
4. Eval (survivor scan): `vpcmpb`/`vpmovmskb`, or `kmask` on
   AVX-512BW.
5. Sieve init: memcpy of tiny-prime pattern. AVX2/512 broadcast stores.
6. Gauss / `re_strip`: 64-bit-word bit matrix, so AVX2/512 row XORs.
- Gate per kernel: counters unchanged, and the speedup recorded in this
  file.

**Phase 4: integration.**
- Turn on the `HAVE_LASIEVE_MPQS`-style paths by default against
  `factor/kmpqs` in `batch_factor.c`, `autofactor.c` and
  `gpu_cofactorization*.c`, and delete the MinGW workaround.
- Measure crossovers against `tinyqs`, `uecm`/`tecm` and SIQS for
  `factor()` of small inputs, and for SIQS DLP/TLP cofactorization.
  Wire it in where it wins.
- Optionally point lasieve5_64 at the shared sources and retire its
  copy and its mpqs `.asm` files.

**Phase 5: retune.** With wider TD and next-pol, the optimal
`mpqs_param`/`mpqs3_param` tables (FB size, sieve length, number of A
factors) will shift. Re-sweep them with a small tuner in the style of
`siqs_tune.py`.

## 6. Open questions

- License, see 4.6.
- Should mpqs3 grow a 4-word (~190-bit) variant, or does SIQS already
  win there? The data stops at 120 bits, so this needs measurement
  first.
- NEON: Phase 3 kernels should keep the same shape as the x86 ones so a
  NEON version slots in, as in the SIQS multi-ISA work.
