# Handoff: vector modular arithmetic (AVX-ECM template cores + RNS)

Status as of 2026-11-07, branch `wip_arith` (pushed, based on master
0a712aa). Read this first, then the documents it points to.

> **Two standing rules from the user:**
> 1. **Don't remove the old AVX-ECM cores from the builds until all of the
>    new work is confirmed working**, i.e. as the last step just before the
>    merge into master (§5).
> 2. **Never merge `wip_arith` into master automatically.** When the
>    merge is ready, stop and ask the user to review first. That means:
>    test results, benchmark results and the diff summary. Merge only on
>    their explicit go-ahead, and the same goes for pushing master.

## 1. Where things stand

### Code

| path | what |
|---|---|
| `arith/vecmont/vecarith52_tmpl.c` | Constant-S template IFMA Montgomery mul/sqr for AVX-ECM (8 numbers per vector, radix 2^52, R = 2^(52*S)). S = 2..32 fully unrolled (`vecmulmod52_s##S`), any larger S through the loop form `vecmulmod52_sN`. `vecmulmod52_tmpl_get(S, &mul, &sqr)`. Bit-identical to the old cores; 93-100% of the IFMA floor. IFMA only (`#if defined(USE_AVX512F) && defined(IFMA)`). |
| `arith/vecmont/vecmont.h` | Owns the digit-size macros (`DIGITBITS`, `base_t`, `VECLEN`), `vec_bignum_t`, `vec_monty_t`, the `vecmont_mul_fn` / `vecmont_sqr_fn` types and the `vecmulmod52_tmpl_get` prototype. `avx_ecm.h` includes it by relative path (`../../arith/vecmont/vecmont.h`), so no include-path changes were needed anywhere. |
| `factor/avx-ecm/avx_ecm.h` | `AVXECM_TMPL_KERNELS` (on by default for IFMA builds; `-DAVXECM_NO_TMPL_KERNELS` to A/B). `avxecm_redc_nwords(nbits)`: exact `ceil((nbits+3)/52)` words with the template, otherwise rounded up to a multiple of 4; always >= 3 spare bits (`AVXECM_SPARE_BITS`). Also the FP-emulation macros used when `IFMA` is undefined (§3.2). |
| `factor/avx-ecm/avx_ecm_main.c`, `avxppm1.c` | `assign_mul_ptrs` and the three P±1 sizing sites use the template cores and `avxecm_redc_nwords`. `insert_mpz_to_vec` now clears a lane's high words (fixed a real stage-2 miss bug). |
| `factor/avx-ecm/bench/` | `ae_bench`: a bit-exact AMM model vs GMP, worst cases, chains, powm, the add/sub range test (`-a`), a negative self-test (`-neg`), and benchmarks (`-b -p`). |
| `arith/rns/` | Standalone RNS library: small regime (`rns.h`, `rns_ref.c`, `rns_avx512.c`), Yao regime setup/reference (`rns_yao.h`, `rns_yao_ref.c`), shared vector helpers (`rns_yao_vec.h`), Yao one-number kernel (`rns_yao_avx512.c`), Yao batch kernel (`rns_yaob_avx512.c`), tests and `rns_bench.c`. IFMA only (`__AVX512F__ && __AVX512IFMA__`). |

### Documents

| path | contents |
|---|---|
| `arith/rns_design.md` | RNS design, the decisions log, all RNS benchmarks; §0d has the current standings. |
| `arith/avxecm_ifma_analysis.md` | Per-kernel analysis of the old AVX-ECM cores, the template design, end-to-end ECM A/B, the add/sub range finding, the stage-2 bug. |
| `arith/avx512_optimization_notes.md` | Machine model, measurement methods, what worked and what didn't. |

### Current speed

From `rns_bench`, Rocket Lake i5-11600K: powm with a full-size exponent,
µs per number, 8 numbers per vector. Full table in `rns_design.md` §0d.

| bits | AE-old | AE-tmpl | Yao1 | Yao8 | mpz_powm | fastest |
|---:|---:|---:|---:|---:|---:|---|
| 256 | 3.4 | **0.9** | 18.6 | 5.2 | 7.3 | AE-tmpl |
| 512 | 12.3 | **5.9** | 40.2 | 16.6 | 45.4 | AE-tmpl |
| 1024 | 57.5 | **43.8** | 97.8 | 55.0 | 264.9 | AE-tmpl |
| 1536 | 209.4 | 149.6 | 173.6 | **147.4** | 854.2 | tie |
| 2048 | 405.2 | 344.0 | 306.4 | **285.4** | 2005.0 | Yao8 |
| 4096 | 2853.7 | 2527.7 | 1500.8 | **1417.6** | 12637.4 | Yao8 |

The RNS small-regime kernel (`rns_v*`) is now dominated at every size.
Keep it as a reference, but don't integrate it.

End-to-end yafu ECM stage 1 with the template plus word-granular sizing:
1.10-1.38x faster than before. A 500-bit N runs about 1.9x faster than the
pre-template code (`avxecm_ifma_analysis.md` §8-9).

## 2. Build and test (WSL, gcc 11, `-O3 -march=native`)

GMP is at `../gmp-install/6.2.1` relative to the repo root (Linux build);
GMP-ECM at `../ecm-install/7.0.6-wsl`, both set in `config.mk`.

```
cd arith/rns && make all bench
./rns_test 50 && ./rns_yao_test 16 && ./rns_bench        # RNS tests + full comparison
make clean && make all CFLAGS="-O3 -g -Wall -DRNSY_COLUMN_FOLD=0" && ./rns_yao_test 16   # exact mode

cd factor/avx-ecm/bench && make && ./ae_bench -t          # AVX-ECM cores, bit-exact
./ae_bench -b -p                                          # benchmarks
```

The top-level `make` builds yafu; the IFMA path needs `USE_AVX512IFMA=1`
or `USE_NATIVE` on an IFMA machine. The template file is already in
`ECM_SRCS` and in both VS projects. **The MSVC build has never been
compiled** with the template.

## 3. Invariants and gotchas (learned the hard way)

* **>= 3 spare bits between N and R.** `vecaddmod52`/`vecsubmod52` act on
  AMM outputs (< R, not < N). With 1-2 spare bits, up to ~20% of lanes
  go wrong (1039-bit N in 1040 bits: ECM almost stopped finding factors).
  Always size through `avxecm_redc_nwords`.
* **Exact word counts are only valid with the template cores.** The
  fixed/BFIPS cores loop over 4-word blocks; `avxecm_redc_nwords` rounds
  up to a multiple of 4 when the template is off. Any code that assumes
  `NWORDS % 4 == 0` (NBLOCKS loops, add/sub, Mersenne paths) must be
  audited before it runs with template sizing (§9 of the analysis did
  this for ECM and P±1).
* **Clear unused high words** when inserting a shorter value into a
  reused vector (the `insert_mpz_to_vec` bug, which cost 2-5 of 40 sigma
  batches at 1249-1252 bits).
* **RNS representation is not radix-2^52.** RNS values are residues, in
  Montgomery form with respect to M_B, lazily reduced. Bounds: mulmod
  outputs < k*N (k = n + 2), mulmod inputs < 8*k*N (`RNS_HEADROOM`),
  `submod` adds 2*k*N. Callers must keep add/sub chains within those
  bounds. N must be coprime to M_B*M_C (the moduli are composite):
  `rns*_ctx_init` returns NULL otherwise, and for ECM that gcd is a found factor.
* **Odd n in Yao:** |C_j|_{b_i} carries a sign (-1)^(n-1) (`fx->cneg`).
  Any new kernel must handle both signs.
* **Benchmarks:** warm up 0.2 s before every timed loop (first-loop
  artifacts of 30-50%), and expect ±5-10% noise (WSL clock 4.2-4.6 GHz).
  Use `taskset` and min-of-N for small differences.
* **This shell collapses `\\` inside heredocs.** Write edit scripts to
  files instead of piping them through heredocs.

## 4. Future work

The order below is a suggestion: 4.4 and 4.5 are mechanical and make 4.1,
4.2, 4.3 and 4.6 easier. 4.7 is a benchmarking task that can run any time.

### 4.1 Replace the old AVX-ECM arithmetic with the template + RNS (user item 1)

**RNS part done (2026-11-07).** The template part (removing the old
fixed/BFIPS selection) is the pre-merge removal step (§5) and is not done.
* `factor/avx-ecm/avxecm_rns.c`: the RNS backend. `vec_monty_t.rns`
  (opaque) holds the shared `rnsy_fixed_t`, a per-thread `rnsy_ctx_t`, a
  single-number buffer and a batch temporary. The five `vec*mod_ptr`
  wrappers call the `rnsyb_*` batch kernels (lane = curve);
  `avxecm_rns_addsub` handles outputs aliasing inputs.
* Values cross between mpz and vectors only through backend-neutral
  helpers: `vec_lane_get` / `vec_lane_set` (lane as integer),
  `vec_set_raw_ui` (raw 1 for leaving Montgomery form), `vec_monty_mulR`
  (x*R mod N, R = 2^(52 NWORDS) or M_B), `vec_set_monty_one`. In
  `avxecm.c` they replace curve setup (`vec_build_one_curve`, the build
  work fcn), both stage-2 batch inversions, the gcd/save extracts and the
  raw-1 vectors (the dead `vececm_old` is untouched). In radix mode they
  do exactly what the old code did. `mdata->n` stays radix 2^52 (it is only
  read back as N).
* Selection (`vec_ecm_main`): non-Mersenne, `size_n >= avx_rns_bits > 0`,
  kernels built, and `avxecm_rns_fixed` succeeds (`rnsyb_choose_n` +
  `rnsy_fixed_init_n` + a ctx test; N sharing a factor with the composite
  moduli falls back to radix). Vectors are allocated with
  `max(NWORDS, n*s)` words so `vecCopy` / `vecClear` cover the layout.
  `assign_mul_ptrs` takes the RNS pointers first; `thread_init` attaches a
  ctx per thread; detach before `vec_monty_free` (owner last).
* Bounds audit: every add/sub feeds a mul or one more add/sub of mul
  outputs; worst is `CROSS_PRODUCT`, (tt3 + Pbprod) - Paprod < 4kN into a
  mul (limit 8kN); subtrahends are mul outputs (< kN < 2kN).
* Option `avxecm_rns_bits` (cmdOptions index 140, `ecm_obj.avx_rns_bits`;
  0 = never). Default `AVXECM_RNS_BITS_DEFAULT` in `include/factor.h`:
  1200 (IFMA) / 1500 (AVX-512F only). `yafu.ini` and the command line
  override it.
* `ecmtune()` (calc function 83, `ecm_rns_tune` in `factor/tune.c`; also
  run by `tune()`): one 8-curve batch, 1 thread, B1 = 20000, fixed sigma,
  radix vs forced RNS at 768..4096 bits (inputs: 256-bit primes times one
  of 257-512 bits, so no small factor), crossover where RNS/radix drops
  below 1 for good (linear interpolation), written to `yafu.ini`.
  Rocket Lake results (s per batch, RNS/radix):

| bits | IFMA radix | IFMA RNS | ratio | no-IFMA radix | no-IFMA RNS | ratio |
|---:|---:|---:|---:|---:|---:|---:|
| 768 | 0.116 | 0.150 | 1.30 | 0.193 | 0.267 | 1.39 |
| 1024 | 0.182 | 0.196 | 1.08 | 0.333 | 0.346 | 1.04 |
| 1280 | 0.273 | 0.257 | 0.94 | 0.423 | 0.462 | 1.09 |
| 1536 | 0.370 | 0.330 | 0.89 | 0.600 | 0.587 | 0.98 |
| 2048 | 0.619 | 0.468 | 0.76 | 0.976 | 0.817 | 0.84 |
| 3072 | 1.301 | 0.811 | 0.62 | 2.115 | 1.448 | 0.69 |
| 4096 | 2.241 | 1.201 | 0.54 | 3.479 | 2.301 | 0.66 |

  Crossover: 1169 bits with IFMA, 1489 bits without (lower than the powm
  benchmark's ~1500 / ~2000: in ECM the larger RNS vectors cost less than
  the powm comparison suggested).
* Verified: with a fixed starting sigma, radix and RNS report identical
  finds (lane, stage, sigma) for 64-curve runs at 520, 1036, 1600, 2048,
  3000, 4000 bits x 3 sigmas (IFMA) and 1036, 2048 (non-IFMA); p15 found
  with RNS at 520-4500 bits (RNS still works at 4500: 2 x 11 channels),
  with 4 threads, and in the non-IFMA build (no `vpmadd52` in it). MSVC
  v143 (avx512ifma) and ClangCL (avx512) compile `avxecm_rns.c`.
  `ae_bench -t`, `rns_test`, `rns_yao_test`, `rns_bench` still pass.
* Not done: P-1/P+1 (`avxppm1.c`) still radix only (unreachable from yafu,
  see §4.2); an RNS `vecmulmod_1` / PARAM1 path (PARAM1 is off);
  timings at large B1 (only the tune's B1 = 20000 batches were timed).

Original notes follow.

Goal: AVX-ECM uses only the template cores below the crossover and an RNS
backend above it. The old fixed/BFIPS cores are deleted or kept only for
A/B tests.

* Template part (straightforward): the template already overrides
  mul/sqr in `assign_mul_ptrs` for non-Mersenne inputs <= 16400 bits.
  - Remove the fixed416-1040/BFIPS selection.
  - Decide what happens above 16400 bits (karatsuba `veckmul_redc` today).
  - Remove the dead `printf("mul > R, reducing")` in
    `vecmulmod52_fixed1040_bfips` if the old cores stay.
  - Mersenne paths (`vecmulmod52_mersenne`, isMersenne != 0) are
    separate and untouched; decide whether they stay.
* RNS part (the real work): AVX-ECM's curve code works on `vec_bignum_t`
  through `vecmulmod_ptr`, `vecsqrmod_ptr`, `vecaddmod_ptr`, `vecsubmod_ptr`,
  `vecaddsubmod_ptr` and a few direct calls (about 180 call sites in
  `avxecm.c`, 140 in `avxppm1.c`). Options:
  1. An RNS "number" type with the same operation set, selected per run.
     The batch kernel's layout (8 numbers mod one N, lane = curve) matches
     AVX-ECM's 8 curves, so use `rnsyb_*`.
  2. Convert at the boundaries only: keep the curve state in RNS for all of
     stage 1 / stage 2, and convert to radix 2^52 / mpz for gcd and output.
  Needed in the RNS API, which doesn't exist yet:
  - combined add/sub (`vecaddsubmod`);
  - multiply by a small constant (the a24 / d constants in the ladder);
  - copies and lane select / masked copy if the curve code uses them;
  - the batch inversion in stage 2. RNS has no cheap comparison, so
    inversion must go through powm or a conversion out to mpz.
  - The headroom rules in §3: audit every add/sub chain between
    multiplies, e.g. (X+Z)(X-Z) and the ladder's differential add.
    Raising `RNS_HEADROOM` or `RNS_SUBMULT` costs a few channel bits.
* Old cores: remove them from the builds but keep the sources (§5); this
  requires splitting `vecarith52.c` first. **Do the removal last**, just
  before the merge, after everything else is confirmed working.
* Crossover: about 1500 bits on Rocket Lake (Yao8 ties AE-tmpl at 1536).
  Measure it inside ECM, not only powm: add/sub/copies cost linearly and
  are cheaper in radix 2^52. Make it a tunable.
* Done when: yafu ECM and P±1 select template or RNS by size; ECM finds
  factors at the same rate as before (use the Python reference / sigma
  batch comparison from `avxecm_ifma_analysis.md` §9 at several sizes,
  including just above and below the crossover); stage-1 and stage-2 timings
  are reported.

### 4.2 AVX-512F without IFMA (FP emulation) (user item 2)

**Done (2026-11-07): template and RNS kernels; see the end of this list.**
* `vecarith52_tmpl.c` builds for every `USE_AVX512F` build. Without
  `IFMA` (`AE_FP`) a product is `t = fma_rz(x, y, 2^104)` (hi, biased) and
  `fma_rz(x, y, 2^104 + 2^52 - t)` (lo, biased): 2 FMA + 1 sub + 2 integer
  adds. Raw bit patterns are summed and the biases removed once per column
  (`np * bias`, mod 2^64). A, B, N are converted to double once per call into
  `mdata->fpwork` (new `vec_monty_t` field, 3 * words vectors, allocated in
  `vec_monty_alloc`), m_k is stored as double in `s`. m_k itself comes out
  of the FP lo sequence as an exact double (no convert on the critical path).
  `AE_TMPL_NO_CARRY_TRICK` is IFMA-only.
* `avx_ecm.h` defines `AVXECM_TMPL_KERNELS` for all AVX-512F builds, so
  non-IFMA yafu now uses the template and exact word sizing. The kept
  routines (add/sub, `vecmul52_1`, `vec_simul_addsub52`) have no IFMA
  branches, so the exact-NWORDS audit of §9 carries over.
* `ae_bench`: `make NOIFMA=1` builds the non-IFMA bench (FMA-calibrated
  cycle; runs on IFMA machines too).
* Verified (Rocket Lake running the non-IFMA code): `ae_bench -t` bit-exact
  for all 75 template cases S = 2..80 (unrolled and `sN`), `-neg` fails as
  intended, `-a` fails only at 1-2 spare bits; IFMA `ae_bench -t` and the
  `rns_bench` AE-tmpl timings unchanged. MSVC v143 and ClangCL compile it
  for `avx512` (not run: `ae_bench` is gcc-only).
* Kernel speed vs the old non-IFMA cores, powm per number: 1.5x at 416
  bits, 1.6x at 1040, 1.4x at 1664 and 1872, 1.3x at 4160.
* yafu ECM, non-IFMA build (`make USE_AVX512=1 USE_AVX512IFMA=`), template
  vs `-DAVXECM_NO_TMPL_KERNELS`, B1=200000, 8 curves, `taskset -c 3`, min
  of 2, seconds:

| bits | words old/new | stage 1 old -> new | | stage 2 old -> new | |
|---:|---:|---:|---:|---:|---:|
| 400 | 8 / 8 | 0.474 -> 0.320 | 1.48x | 0.286 -> 0.218 | 1.31x |
| 520 | 12 / 11 | 0.948 -> 0.570 | 1.66x | 0.501 -> 0.350 | 1.43x |
| 700 | 16 / 14 | 1.604 -> 0.919 | 1.75x | 0.804 -> 0.520 | 1.55x |
| 1036 | 20 / 20 | 2.412 -> 1.881 | 1.28x | 1.088 -> 0.946 | 1.15x |
| 1600 | 32 / 31 | 5.164 -> 3.585 | 1.44x | 2.280 -> 1.659 | 1.37x |
| 2048 | 40 / 40 | 7.710 -> 5.839 | 1.32x | 3.302 -> 2.585 | 1.28x |

* ECM factor finding, non-IFMA template build: p15 found (64 curves,
  B1=3000) at 400, 520, 700, 1036, 1100, 1250, 1600, 2048, 3000 bits.
* P-1/P+1 (`avxppm1.c`) is **not reachable from yafu**: `pp1.c` runs
  GMP-ECM (`if (1)` before the `vecPP1` branch) and `vecPM1` is only called
  from the AVX-ECM demo (`factor/avx-ecm/demo`, VS project only). Its three
  sizing sites compile with the template but were not run.
* **RNS kernels: done too** (`rns_design.md` §0e). `RNS_HAVE_AVX512` (any
  AVX-512F target) now guards the kernels, `RNS_HAVE_IFMA` selects the
  native products; `rns_vec52.h` holds `LO52`/`HI52` (IFMA or exact FP
  emulation) and the fast FP accumulate helpers. `arith/rns/Makefile`
  takes `NOIFMA=1`. Verified non-IFMA: `rns_test`, `rns_yao_test` in both
  column-fold modes, `rns_bench`, on WSL and as MSVC builds on Windows;
  no `vpmadd52` in the non-IFMA yafu binary; IFMA results and timings
  unchanged. Non-IFMA yafu (template + FP RNS compiled in) builds clean
  and finds p15s at 700 / 1250 / 2048 bits.
* **Not tested on real non-IFMA hardware** (Skylake-X / Cascade Lake):
  all non-IFMA runs were the non-IFMA code on Rocket Lake. Port balance
  differs there (2 FMA units on server parts), so timings will differ.

Original notes follow.

Goal: everything still builds and works on AVX-512F machines without
IFMA (Skylake-X, Cascade Lake), using the FP-emulated 52-bit products.

* Existing mechanism: in `avx_ecm.h`, under `#ifndef IFMA` (from about line
  471), `VEC_MUL_ACCUM_LOHI_PD` and friends emulate `madd52lo/hi` with
  `_mm512_cvtepu64_pd` plus two `_mm512_fmadd_round_pd(..., _MM_FROUND_TO_ZERO)`
  and bias constants (`dbias`, `vbias1-3`). The old cores use these macros.
* Today: the template and all RNS kernels compile to nothing without
  IFMA, so a non-IFMA AVX-512 build uses the old cores with block sizing.
  **Not tested since these changes.** First step: build yafu with
  `USE_AVX512=1` and without `USE_AVX512IFMA`, run `ae_bench -t`
  (needs a non-IFMA build of the bench: drop `-DIFMA -DUSE_AVX512IFMA`)
  and an ECM run. Rocket Lake can run that code too; it just won't be
  representative for timing.
* Then port: give the template and the RNS kernels a lo/hi multiply
  abstraction (the RNS code already uses `LO52`/`HI52` macros in
  `rns_yao_vec.h` and `rns_avx512.c`; the template uses `AE_LO`/`AE_HI`)
  and implement it with the FP sequence.
  - The FP sequence yields lo and hi from one pair of FMAs, so a
    `MUL_LOHI` primitive is the natural unit, not separate lo and hi.
  - The cvt to double costs extra µops per input. Hoist conversions of
    reused operands (e.g. broadcast constants, row operands).
  - Port balance differs: FMA on p0 (+ p5 on 2-FMA server parts), so
    redo the bound from `avx512_optimization_notes.md` §1.
* Done when: non-IFMA builds pass `ae_bench -t`, `rns_test`,
  `rns_yao_test` (exact mode too), and ECM finds factors. Speed is
  compared to the old non-IFMA cores.

### 4.3 AVX2 versions of the template and RNS kernels (user item 3)

**Experiment done (2026-11-07): `arith/avx2_experiments.md`,
`arith/avx2mont/`.** Three 4-lane AVX2 Montgomery kernels with the
template's algorithm, all bit-exact vs a GMP AMM model at 100-5000 bits:
FP52 (52-bit digits, two FMAs per product, MXCSR round-toward-zero per
call), D28 (`vpmuludq`, 56-bit products summed whole), D32 (`vpmuludq`,
64-bit products with split high halves). Powm per number vs `mpz_powm`
(Rocket Lake running AVX2 code): **FP52 1.6-2.0x faster at 256-4096 bits**,
D28 ~1.0-1.3x (parity above 1 kbit), D32 0.83-1.03x. FP52 is the AVX2
format; 28/32-bit digits lose on product count. Full unrolling only pays
up to S = 12 on AVX2 (16 registers). Not measured on an AVX2-only CPU; RNS
on AVX2 and an AVX2 AVX-ECM are not started (next steps in the doc).

Original notes follow.

Goal: measure whether 4-lane AVX2 versions beat `mpz_powm`, and by how
much. Today AVX-ECM exists only for AVX-512 (`#ifdef USE_AVX512F` in
`avx_ecm.h` and `factor_common.c`).

* User direction: 32-bit digits, or 28-bit digits so lane accumulation has
  headroom. Experiment with both (§5), and with the FP FMA-based approach below.
* AVX2 has no 52-bit multiply and no 64x64 -> 128 vector multiply. Options:
  - The FP trick at 256 bits: `_mm256_fmadd_pd` with no embedded rounding
    control. Rounding mode must then come from MXCSR (set round-toward-zero
    around the kernel) or the error-free product needs a different formulation
    (Dekker / FMA-based two-product with round-to-nearest). Check the
    exactness argument carefully.
  - 32-bit digits with `_mm256_mul_epu32` (32x32 -> 64, 4 lanes). There's
    the `DIGITBITS == 32` path in `avx_ecm.h` (16 lanes on AVX-512) as a
    model, and `vecarith.c`.
* RNS fits AVX2 naturally with smaller channels (e.g. 26-31-bit moduli
  with `mul_epu32`, or 50-bit FP channels); the Yao constant-size
  analysis changes with the word size.
* Compare per-number powm against `mpz_powm` at 256-4096 bits. GMP's
  mulx/adx assembly is strong on AVX2-class cores, so the win isn't
  guaranteed. Measure before committing to a full port.
* The SIQS multi-ISA work in this repo (SSE4.1/AVX2/AVX-512/NEON kernel
  sets) is a pattern for dispatch and Makefile flags.

### 4.4 Move the template arithmetic into arith/ (user item 4)

**Done (2026-11-07).** `vecarith52_tmpl.c` moved to `arith/vecmont/`
(`git mv`, function names unchanged); the types moved to
`arith/vecmont/vecmont.h`, which `avx_ecm.h` includes; the template no
longer includes `avx_ecm.h` or ytools (own `AE_UNUSED` shim). Updated:
top-level Makefile `ECM_SRCS`, `factor/avx-ecm/bench/Makefile`,
`arith/rns/Makefile` (`VM_DIR`), and both VS ecm projects (+ filters;
MSVC still not compiled). Verified: yafu (USE_NATIVE, IFMA) builds,
`ae_bench -t`/`-a` pass, `-neg` fails as it should, `rns_test 50`,
`rns_yao_test 16` and `rns_bench` pass with unchanged timings, yafu ECM
finds a p15 at 400/1036/1600/2048 bits (64 curves, B1=3000) with the
template core selected. Original notes follow.

Goal: `arith/` holds the vector arithmetic (template + RNS); AVX-ECM is
a consumer.

* Move `vecarith52_tmpl.c` (and possibly the `vec_bignum_t` / `vec_monty_t`
  types it needs) to `arith/`, e.g. `arith/vecmont/` next to `arith/rns/`.
  Keep the function names, or add a thin header so `avx_ecm_main.c` and
  `avxppm1.c` keep compiling.
* The template depends on `avx_ecm.h` types and macros. Either move the
  minimal type definitions into an arith header that `avx_ecm.h`
  includes, or keep `avx_ecm.h` as the type owner and include it from
  arith (worse layering).
* Update the include paths in the Makefile (`-Ifactor/avx-ecm` today),
  `factor/avx-ecm/bench/Makefile`, `arith/rns/Makefile` (`AE_DIR`), and
  both VS projects (`build.vc22/ecm/ecm.vcxproj(.filters)` and
  `build.vc22/ecm-clang/ecm.vcxproj`).
* Done when: yafu, `ae_bench` and `rns_bench` all build from the new
  location and the tests pass. This is a no-behavior-change commit.

### 4.5 Integrate the template and RNS into yafu's Makefile (user item 5)

**Done (2026-11-07), except removing the old cores (last step, §5).**
* The five RNS sources are in `ECM_SRCS` and in both VS ecm projects
  (+ filters, headers listed). The references/setup always build; the
  kernels compile to empty objects without IFMA.
* One guard, `RNS_HAVE_AVX512IFMA` in `rns.h`: compiler IFMA macros **or**
  yafu's `USE_AVX512F` + `IFMA` (MSVC has no IFMA macro). The vector
  prototypes in `rns.h` / `rns_yao.h` and the vector parts of `rns_test` /
  `rns_yao_test` use it too. `<immintrin.h>` is included only inside it.
* Portability shims in `rns.h`: `RNS_INLINE`, `RNS_UNROLL(n)`,
  `RNS_LIKELY`. `rns_ref.c` uses `mul64` / `mulmod64` helpers: `__int128`
  where `__SIZEOF_INT128__` is defined, else `_umul128` / `_udiv128`.
  `rns_test.c` uses `timespec_get` on Windows. Allocation still uses the
  library's own `rns_aligned_alloc` (kept standalone, no yafu deps).
* yafu builds with `-Wconversion`; the RNS files are now warning-free
  under it (casts only, no behavior change).
* **Gotcha:** the local `config.mk` sets `USE_AVX512IFMA := 1`, so a
  plain `make USE_AVX2=1` is still an IFMA build. Pass
  `USE_AVX512IFMA=` (empty) to build lower ISA levels.
* Verified: WSL `make yafu` for generic / AVX2 / AVX-512F / IFMA / native
  (rc 0; only the two pre-existing `VEC_MAXDIGIT` / `VECLEN_I` macro
  warnings, now in `vecmont.h`). MSVC v143 compiles the ecm project for
  `avx512` and `avx512ifma` (IFMA kernels 330-450 KB objects, stubs
  otherwise); ClangCL `avx512` compiles. `rns_test` and `rns_yao_test`
  (both column-fold modes) built with `cl` against MPIR pass natively on
  Windows, vector kernels included. The ecm project's library step fails
  only on its hard-coded `libecm.lib` path (pre-existing, unrelated).
  ClangCL `avx512ifma` fails in the old `vecarith52.c` too: ClangCL never
  gets `-mavx512ifma` (`build_config.props` says IFMA is not in the
  Windows matrix). WSL tests: `rns_test`, `rns_yao_test` (both modes),
  `rns_bench`, `ae_bench -t` pass; `ae_bench -a` fails only at 1-2 spare
  bits. ECM p15 (64 curves, B1=3000) is found with the native IFMA build,
  the AVX-512F (old cores) build and the AVX2 (GMP-ECM) build; repeated
  runs on the native build found it in 46/48 runs (misses are the normal
  per-run chance at this B1).

Original notes follow.

* Template: already in `ECM_SRCS` (Makefile around line 950), now as
  `arith/vecmont/vecarith52_tmpl.c` (§4.4 done).
* RNS: currently only built by `arith/rns/Makefile`. Add `rns_ref.c`,
  `rns_avx512.c` (the small regime: keep it in, per §5), `rns_yao_ref.c`,
  `rns_yao_avx512.c`, `rns_yaob_avx512.c` to the yafu build.
* The old cores come out of `ECM_SRCS` and the VS projects only as the
  final pre-merge step (§5). Until then they stay in the build.
  - The kernels are guarded by `__AVX512F__ && __AVX512IFMA__` (compiler
    macros from `-march`), not by yafu's `USE_AVX512IFMA` / `IFMA`. Make
    those consistent: check both, or switch to yafu's defines.
  - The references and setup need GMP only. They can be built always and
    are needed for setup even when the vector kernels are compiled out.
  - The RNS code uses `posix_memalign` / `_aligned_malloc` itself. Consider
    yafu's `xmalloc_align` / `align_free` for consistency.
* Add the files to the VS projects too, and compile MSVC and clang-cl.
  RNS uses GCC-isms: `__attribute__((always_inline))`,
  `_Pragma("GCC unroll")`, `__builtin_expect`, `unsigned __int128` in
  `rns_ref.c`. The template already has MSVC macro shims; copy that
  approach. `__int128` needs `_umul128` / `__umulh` on MSVC.
* Done when: `make` (IFMA, AVX-512F only, AVX2, generic) builds clean on
  Linux/WSL and MSVC builds, with RNS compiled out where unsupported.

### 4.6 Keep AVX-ECM working with the new kernels (user item 6)

A test matrix to run after each of §4.1-4.5:

* `factor/avx-ecm/bench/ae_bench -t` (bit-exact cores, worst cases,
  1/2/26 spare bits) and `-a` (add/sub range), plus `-neg` (the negative
  self-test catches a corrupted kernel).
* `arith/rns`: `rns_test`, `rns_yao_test` in both `RNSY_COLUMN_FOLD`
  modes, and `rns_bench`, which checks every engine's powm against
  `mpz_powm` in every lane.
* yafu ECM end to end. `avxecm_ifma_analysis.md` §8-9 has the method:
  - p15 / p20 test inputs, B1 = 3000 / 11000, 64 curves;
  - count curves that find p, before vs after; the sigma-batch lane-for-lane
    comparison against the exact Python reference for stage 1 + stage 2;
  - sizes 400, 520, 700, 1000, 1038-1040, 1100, 1249-1252, 1500, 2048, 3000;
  - include sizes near every word boundary and near the template/RNS crossover.
* P-1 / P+1 (`avxppm1.c`, three sizing sites) at a few sizes.
* Mersenne / special inputs (2^n ± 1 and factors of them) still take the
  Mersenne paths.
* Thread counts > 1 (the RNS contexts have per-ctx scratch:
  **one ctx per thread**).
* Done when: all of the above pass on IFMA and non-IFMA AVX-512 builds,
  with factor-finding rates statistically unchanged.

### 4.7 Above 16400 bits (decided: benchmark three options)

AVX-ECM switches to karatsuba above `AVXECM_KARATSUBA_BITS` = 16400 today
(`veckmul_redc` / `vecksqr_redc` in `vecarith52_special.c`, built on
`kmuln` / `ksqrn` / `tmuln`). Benchmark these, at 16k-64k bits:
1. The existing karatsuba REDC (`veckmul_redc`, `vecksqr_redc`).
2. The Yao kernels as they are. `RNSY_MAXS` = 10 words per channel and
   the n <= ~11 constant limit cap the Yao setup near 4.1k bits today
   (`rnsy_fixed_init_n` returns NULL above that). Raise `RNSY_MAXS`, or
   accept larger s through a non-unrolled loop form like the template's `sN`.
3. Yao extended with karatsuba channel products. Channel products are
   5 per channel pair, so subquadratic channel multiplication pays once s
   is large; the extensions stay O(n^2 s).
   - The cost model in `rns_design.md` §2 becomes n * K(s) + 4 n^2 s, which
     pushes the optimum toward fewer, wider channels.
   - Squaring is especially cheap with karatsuba.
   - FFT/NTT channels are the eventual step beyond this.
* Also include the template `sN` loop form (no karatsuba) for reference,
  and `mpz_powm`. `rns_bench` can be extended: sizes above 4096, engines
  for karatsuba and Yao-large.

## 5. Decisions (user, 2026-11-07)

* **Above 16400 bits:** benchmark the karatsuba REDC, the Yao kernels, and
  Yao extended with karatsuba channel products (§4.7).
* **Old cores: don't delete them, but take them out of the builds** (Makefile
  and the MSVC / clang-cl projects). **Keep the RNS small-regime kernel and
  add it to the build files** for further experiments.
  - Complication: `vecarith52.c` mixes the old mul/sqr cores
    (`vecmulmod52_fixed*_bfips`, `vecsqrmod52_fixed*_bfips`, `vecmulmod52`,
    `vecsqrmod52`, the `_ifma` variants) with routines that must stay:
    `vecaddmod52`, `vecsubmod52`, `vecsignedaddmod52`,
    `vec_simul_addsub52`, `vecmul52_1` (multiply by a small constant),
    `vecsqrmod52_mul`, the Mersenne routines and `vecmodexp52`.
  - `vecarith52_special.c` holds the Mersenne and karatsuba routines,
    which stay in the build (§4.7).
  - So split `vecarith52.c`: the kept routines go to a new file that stays
    in `ECM_SRCS`, the old cores move to e.g. `vecarith52_oldcores.c`
    outside the builds, with the A/B harnesses (`ae_bench`, `rns_bench`)
    still compiling it. `assign_mul_ptrs` then needs a template-only path
    (remove the fixed/BFIPS selection, keep a clear error if
    `vecmulmod52_tmpl_get` fails).
  - **Timing (user decision):** keep the old cores in the builds until all
    new work is confirmed working. Their removal is the last change before
    the merge into master. The `vecarith52.c` split itself can happen
    earlier, as long as both halves stay in the build. Non-IFMA AVX-512 builds
    use the old cores today, because the template is IFMA-only, so the
    template's FP-emulation path (§4.2) must be done and confirmed before
    the removal. Otherwise keep `vecarith52_oldcores.c` in the build for
    `!IFMA` only.
  - `vecarith.c` (`DIGITBITS == 32`, 16 lanes) is a separate path; decide
    alongside the AVX2 work.
* **AVX2 digit size:** 32-bit, or possibly 28-bit digits so lane
  accumulation has headroom. This needs experiments (§4.3).  Also FP FMA.
  - With 32x32 -> 64 `vpmuludq` products, 28-bit digits give 56-bit
    products and 8 bits of in-lane headroom: about 256 accumulated terms
    before a carry pass.
  - 32-bit digits need a split or carry handling on every accumulation.
  - The same choice applies to the RNS channel size on AVX2 (28-31-bit
    moduli fit `vpmuludq` too; see the FCFS/Yao constant analysis in
    `rns_design.md`, which must be redone for small words).
* **Merge `wip_arith` to master** once AVX-ECM is confirmed working with
  the new kernels (§4.6) and shows a speedup. **Never merge
  automatically.** Pre-merge sequence:
  1. Every item is confirmed working: the §4.6 test matrix passes on IFMA
     and non-IFMA builds.
  2. Remove the old cores from the builds (Makefile, VS projects).
  3. Rebuild and rerun the full test matrix and benchmarks.
  4. **Stop and ask the user to review**: a test summary, the benchmark
     table vs master, `git log master..wip_arith` and `git diff --stat
     master..wip_arith`.
  5. Merge (and push master) only on the user's explicit go-ahead.
