# SIQS large-prime strategies across CPUs and ISAs: plan

*Drafted 2026-10-22. Planning only, no code yet. Companion to
`direct_sieve.md` and `Sink-Stripe Sieve Stores.md`.*

## 1. Goals

1. **Mix and match per prime range**, chosen at **run time**: buckets,
   sparse sieve or direct sieve, each with a selectable scatter kernel. One
   binary can then be A/B tested on any CPU without rebuilding.
2. **Port the direct sieve, the sparse sieve and the scatter kernels**
   (bit loop, sink stripe, compress store) to AVX2, SSE4.1, SSE2, generic C
   (any architecture), and later NEON.
3. **Fast, trustworthy comparison** on a new machine: one build script, one
   bench script, and a correctness gate that is bit-identical across
   strategies and ISAs.
4. **Defaults per CPU class** chosen from measurements, overridable from
   the command line and yafu.ini, and tunable by `siqs_tune.py`.

## 2. Where things stand (the L2sieve branch)

- **Compile time only.** `TRY_L2_SIEVE`, `TRY_SPARSE_SIEVE`,
  `TRY_DIRECT_SIEVE`, `DIRECT_MEDBUCKET`, `DIRECT_SCATTER_MODE`,
  `DIRECT_SPARSE_COMPRESS`, `DIRECT_PREFILTER_BITS` and
  `DIRECT_CMP_MAX_SURVIVORS` are all `#define`s. Trying a different
  strategy means rebuilding.
- **AVX-512 only.** The new code requires AVX-512F/BW:
  - `TRY_DIRECT_SIEVE` `#error`s without it.
  - `large_sieveblock_32k_avx512bw`, `L2_sieveblock_32k_avx512bw`,
    `direct_*` and `tdiv_large_32k_avx2` (AVX-512 despite its name) have no
    lower-ISA versions.
  - process_poly calls them without an ISA guard. So a USE_AVX2 or plain
    build of this branch probably fails to link or misbehaves (checked
    first in phase 0).
- **Deferred finalization is in only one scan.** It lives in
  `check_relations_siqs_16_avx512`. `_16_avx2` and `_16_sse2` still finalize
  each report as the block is scanned.
- **What already works at run time.** Kernel pointers (`med_sieve_ptr`,
  `lp_sieveblock_ptr`, `tdiv_*_ptr`, `nextRoots_ptr`, `scan_ptr`) are picked
  from `HAS_AVX512BW/F`, `HAS_AVX2` and `HAS_SSE41`. The stock bucket path
  exists for every ISA. That is the baseline everywhere.
- **Combinations that already exist in the code**, each behind ifdefs:

  | [med_B, large_B) | [large_B, B) | current switch |
  |---|---|---|
  | buckets | buckets | stock (no TRY_*) |
  | sparse | buckets | TRY_L2_SIEVE + TRY_SPARSE_SIEVE |
  | buckets | direct | TRY_DIRECT_SIEVE, DIRECT_MEDBUCKET 1 |
  | sparse | direct | TRY_DIRECT_SIEVE, DIRECT_MEDBUCKET 0 (default now) |

  All four combinations are already implemented. The main structural job
  is to make them selectable at run time.
- **Detection bug.** On the i5-11600K yafu reports L2 = 12582912 bytes,
  which is the L3 (the L2 is 512 KB). Any policy that uses L2 size needs
  this fixed first (`ytools/ytools.c`, `top/driver.c`).

## 3. How the problem decomposes

Primes split into ranges. Each range gets an **algorithm**, and some
algorithms take a **kernel** for their inner scatter or match loop.

| range | algorithms | kernels | notes |
|---|---|---|---|
| < med_B | med sieve (existing per-ISA versions) | n/a | out of scope; already per ISA |
| [med_B, large_B): 32 KB < p ≤ interval | **BUCKET** (medbucket + lp_sieveblock + tdiv_LP), **SPARSE** (per-block large_sieveblock + tdiv_large) | SPARSE: bitloop, sink, cstore-reg, cstore-mem | at most one hit per root per block |
| [large_B, B): p > interval | **BUCKET** (largebucket), **DIRECT** (whole-interval scatter + fused resieve/update/scatter) | DIRECT: bitloop, sink, cstore-reg, cstore-mem; resieve filter: gather, compare, scalar | at most one hit per root per side |

Kernel definitions (lane width W = 16, 8 or 4, or scalar):

- **bitloop**: compare → mask → `tzcnt`/`blsr` loop. One unpredictable
  exit per vector.
- **sink**: blend misses to a rotating stripe, then W unconditional
  read-modify-writes. No branches; wasted stores when the hit rate is low.
- **cstore-reg**: compress in a register (or emulate with a LUT shuffle),
  unaligned store, batch, then one sieve loop per batch.
- **cstore-mem**: `vpcompressd` straight to memory (AVX-512 only).
- **scalar-append** (new, for scalar and NEON): `buf[n] = x; n += (x < L);`
  is a branch-free scalar compress, so batching works with no SIMD at all.

## 4. Run-time strategy selection

### 4.1 Strategy record

Replace the compile-time switches with one struct in `static_conf_t`,
filled once per factorization:

```
typedef struct {
    uint8_t mid;        // MID_BUCKET | MID_SPARSE           [med_B, large_B)
    uint8_t large;      // LARGE_BUCKET | LARGE_DIRECT       [large_B, B)
    uint8_t scatter;    // K_BITLOOP | K_SINK | K_CSTORE_REG | K_CSTORE_MEM | K_SCALAR_APPEND
    uint8_t sparse_k;   // same kernel set, for the sparse sieve and tdiv_large
    uint8_t rs_filter;  // RS_GATHER | RS_COMPARE | RS_SCALAR   (direct resieve)
    uint16_t pf_bits;   // 0 / 256 / 512 / 1024 prefilter bits
    uint16_t cbuf_flush;
    uint16_t sink_bytes;
} siqs_strategy_t;
```

Command line: `-siqsMid bucket|sparse`, `-siqsLarge bucket|direct`,
`-siqsKernel bitloop|sink|creg|cmem|scalar`, plus `-siqsSparseKernel`,
`-siqsPF`, `-siqsCBuf`. The same keys work in yafu.ini. At `-v` yafu
prints the resolved strategy, as the EXP_LINE block does now, so logs
record what ran.

The existing `TRY_*` macros become **build-time capability** switches
("compile this code in"), not mode switches. Every ISA build includes all
four combinations.

### 4.2 What each choice affects (keep this table current)

| choice | fb bounds (SIQS.c) | root updates | per-block sieve | scan / tdiv | memory |
|---|---|---|---|---|---|
| mid = BUCKET | fb_blk8i = med_B | firstRoots middle loop + nextRoots medbucket | lp_sieveblock_ptr | tdiv_LP_ptr per report | bucket slices |
| mid = SPARSE | fb_blk8i = large_B | nextRoots L2 section | large_sieveblock | tdiv_large | sieve_p/n sized to fb_blk8i |
| large = BUCKET | large_B rounds down | largebucket | lp_sieveblock_ptr | tdiv_LP_ptr | bucket slices |
| large = DIRECT | large_B rounds up (first 16-aligned p > interval) | fused pass (no largebucket) | none per block (pass 1 / fused) | deferred: save_survivor → fused resieve → finalize | (2NB+4)·32 KB sieve, pool, bitmap |

Rules to enforce in an init-time validator:

- **DIRECT** requires a scan with the deferred-finalization hook (§5.6),
  and a-primes below large_B (below med_B when mid = BUCKET). Today that
  is a runtime `exit(1)`; turn it into a fallback to BUCKET with a warning.
- **SS search, batch polys and TLP** stay excluded with DIRECT until
  tested. The validator forces large = BUCKET for them.
- **process_poly has two shapes**: block-streaming (large = BUCKET) and
  whole-interval/deferred (large = DIRECT). Branch once per poly, or split
  it into two functions (`process_poly_stream`, `process_poly_direct`) and
  pick one with a pointer. Splitting is cleaner: the current ifdef nesting
  is the main source of bugs, and the two shapes differ in memset, block
  loop, the minblock trick, the save/restore of `blk_sieve[65535]`, and the
  end of each poly.

### 4.3 Dispatch granularity

Never branch on the kernel inside a hot loop. Instead:

- Write each kernel-parameterized pass once in a template include
  (`direct_pass_tmpl.h`) that uses `KERNEL_*` macros, and include it once
  per (ISA, kernel) pair to produce
  `direct_fused_pass_avx2_creg`, `..._avx512_cmem`, and so on.
- Select the function pointer at init from `(HAS_*, strategy.scatter)`.
  This is how yafu already handles the med sieve per ISA.
- If some (ISA, kernel) pairs are never worth having, leave them
  unimplemented; the validator maps them to the nearest available kernel
  and logs it.

## 5. Porting the kernels to lower ISAs

### 5.1 Design rule: bit-identical results across ISAs and kernels

**Phase 0 finding:** stock builds at different ISA levels do *not* do
identical work. On c77 at NB=4, the polynomial counts were 93182 (AVX2),
91142 (SSE4.1), 90887 (SSE2) and 90248 (generic), because their med-sieve
kernels approximate logp differently. Cross-ISA bit-identity is therefore
only possible for the large-prime kernels. The gate becomes:

- **Same binary:** the AVX-512 build also compiles in the AVX2 and scalar
  versions of every new large-prime kernel. A run-time option
  (`-siqsLPisa avx512|avx2|scalar`) selects them, so each kernel port is
  checked bit-for-bit against the AVX-512 kernel inside one binary, where
  everything else is identical.
- **Per build:** lower-ISA builds must factor, and their yield per poly
  must stay within a band of the same build's stock bucket strategy.

- **Keep the 16-prime logp group on every ISA.** AVX2 handles a group as
  2×8 lanes, SSE and NEON as 4×4. Keep large_B and fb_blk8i 16-aligned.
- **Keep the order of fb_offsets.** Matches are visited in ascending fb
  index, as the current compress code already does.
- **Then any (ISA, kernel, strategy-with-same-ranges)** gives the same
  sieve bytes, reports and relations under `-seed N -noopt`. That is the
  correctness gate for every port (§7.2).

### 5.2 Primitives and how each ISA provides them

| primitive | AVX-512 | AVX2 (W=8) | SSE4.1 (W=4) | SSE2 (W=4) | NEON (W=4) | scalar |
|---|---|---|---|---|---|---|
| u32 `a < b` → mask | `cmplt_epu32_mask` | `cmpgt_epi32` (all values < 2³¹) + `movemask_ps` | `cmpgt_epi32` + `movemask_ps` | same | `vcltq_u32`; mask via `vandq` {1,2,4,8} + `vaddvq_u32` | compare |
| blend | `mask_blend` | `blendv_ps` | `blendv_ps` | and/andnot/or | `vbslq_u32` | cmov |
| min/max u32 | native | `min/max_epu32` | `min/max_epu32` | xor-sign + `cmpgt` + blend | `vminq/vmaxq_u32` | cmov |
| compress | native (reg or mem) | 256-entry LUT → `cvtepu8_epi32` → `permutevar8x32` (avoid the pdep trick: microcoded on Zen 1/2) | 16-entry `pshufb` LUT (SSSE3) | none; use scalar-append | 16-entry `vqtbl1q_u8` LUT | scalar-append |
| popcount | `popcnt` | `popcnt` | `popcnt` | table or `popcnt` if present | `vcnt` / `__builtin_popcount` | builtin |
| gather (resieve bitmap) | `i32gather` | `i32gather` (slow on Zen 1–3) or scalar per surviving lane | scalar | scalar | scalar | scalar |
| prefilter table lookup | `permutex2var` (1024b) | `permutevar8x32` ×2 + blend (512b) | scalar 128-byte table | scalar | `vqtbl4q_u8`-based or scalar | scalar |
| byte `-= 4` (minblock) | `sub_epi8` | `sub_epi8` | `sub_epi8` | `sub_epi8` | `vsubq_u8` | autovectorized (`-O3`) or u64 SWAR |

Notes:
- **AVX2 signed compares** are safe: interval < 2³¹, and roots and primes
  are < 2³¹. This is what the sink doc's AVX2 note relies on.
- **Compress LUTs** are 2 KB (AVX2, 8-byte entries) and 64 bytes
  (SSE/NEON), so they stay in L1.
- **A thin `siqs_simd.h` vector layer** (`V_LOAD`, `V_CMPLT_MASK`,
  `V_BLEND`, `V_CSTORE`, …) for W = 16/8/4/1 would let the template
  includes in §4.3 be written once. Use it only for the new kernels; leave
  the existing hand-tuned med sieve and tdiv code as they are.

### 5.3 Direct scatter (pass 1 and the fused next-poly scatter)

- **bitloop and sink** port mechanically.
- **Sink for scalar and SSE2:** this is MPU's original form, one cmov per
  root. It is the natural default for the generic build and for SSE2.
- **cstore-reg on AVX2:** LUT compress, two ymm halves per 16-prime group.
- **cstore on SSE4.1 and NEON:** four 4-lane compresses per root vector,
  or just scalar-append. Measure both; at W = 4 the LUT shuffle may not
  beat the branch-free scalar append.
- **Packing** stays `logp<<24 | pos` (2·interval < 2²⁴ for NB < 256).

### 5.4 Fused resieve and root update (direct)

- **Root update** (sub/add mod p, unordered pair) is plain SIMD on every
  ISA.
- **Survivor filter.** Keep the order: in-range mask → prefilter → test.
  - Prefilter on AVX2: 512 bits (two `vpermd` + blend) or 256 bits (one
    `vpermd`).
  - On ISAs without gather (or where gather is slow), test the few lanes
    that survive the prefilter with scalar bitmap tests: get the mask with
    movemask and visit its bits. Lanes are rarely left, so this branch
    predicts well; unlike the sieve, it is not the hot case.
  - `RS_COMPARE` (broadcast and compare against each survivor position)
    ports directly and may beat gather-less bitmap tests at low survivor
    counts.

### 5.5 Sparse sieve and tdiv_large

- **large_sieveblock** with the "each root hits at most once" rule maps
  onto any W: compare root1/root2 < 32768 → kernel → masked add of p →
  subtract 32768 → min/max store. SSE2 has no min_epu32, so emulate it.
- **tdiv_large** compresses the matching fb indices. On scalar and SSE2
  use scalar-append of `i + lane`.

### 5.6 Scan hook for deferred finalization

- **Factor out the per-report tail** of the scans into one function
  (`siqs_report_tail(j, parity, blocknum, …)`). It runs tdiv_LP for
  mid = BUCKET, then either `direct_save_survivor` (large = DIRECT) or
  finalizes on the spot.
- **Call it from all three scans** (`_16_avx512`, `_16_avx2`, `_16_sse2`)
  and drop the AVX-512-only restriction.
- **Check the same block offset bugs on those scans** that bit the
  AVX-512 one: filter_SPV's blocknum offset, and the BW branch's sieve
  pointer (see L2_SIEVE_DEBUG_NOTES §4). The `dconf->sieve` re-pointing
  from process_poly already makes the scans layout-agnostic.

### 5.7 Generic and ARM builds

- **No intrinsics.** Use scalar kernels: sink for the direct scatter,
  scalar-append + batch for the direct scatter and the sparse sieve
  (measure both), scalar root update and scalar resieve filter.
- **Bit-identical to x86.** The generic build should match the x86
  builds under `-seed -noopt`, so it also serves as the reference
  implementation.
- **NEON later**: swap the scalar kernels for W = 4 NEON ones behind
  `USE_NEON`, using the same template includes. SVE2 would bring
  `svcompact` (native compress) and predicates, and can be added when
  hardware is available.

## 6. Choosing defaults

### 6.1 Policy

`siqs_strategy_t` is resolved in this order:

1. explicit command-line or yafu.ini values;
2. a **built-in defaults table** keyed by
   (ISA level, CPU vendor/family class, digits bucket, NB);
3. fallback: whatever the validator allows, closest to the table entry.

Seed values for the table, from what is known so far (all to be
re-measured):

| CPU class | mid | large | scatter | sparse kernel |
|---|---|---|---|---|
| Intel AVX-512 (Rocket/Ice/Sapphire Lake) | sparse | direct | cmem | cmem |
| AMD Zen 4/5 (AVX-512) | measure | direct | **creg** | **creg** |
| AVX2 (Intel and AMD) | measure (bucket likely) | measure | creg vs sink | creg vs bitloop |
| SSE4.1 / SSE2 | bucket | measure | sink | bitloop |
| generic / ARM | bucket | measure | sink vs scalar-append | scalar-append |

On Zen 4 the memory form of `vpcompressd` is reported to be microcoded
and much slower than the register form. That matches what you saw in the
bucket sort, so creg is the starting point for AMD. Confirm it in
phase 3/6.

### 6.2 Cache-aware rule (needs the L2 detection fix)

The direct scatter target is 2·NB·32 KB. On Rocket Lake (512 KB L2) the
direct sieve fell behind buckets as NB grew, before the medbucket and
compress work. A first heuristic to test: **DIRECT only when 2·NB·32 KB ≤
L2**, otherwise BUCKET for [large_B, B). Zen 4 (1 MB L2) and newer Intel
cores (2 MB L2) would then use DIRECT at larger NB. Measure this before
relying on it.

### 6.3 Tuning

`siqs_tune.py` already handles categorical parameters (SPV). Add `Mid`,
`Large`, `Kernel` and `SparseKernel` as categorical params. Its
best-per-size output then becomes a candidate row for the defaults table
on that machine. A `-siqsStrategyProbe` mode (sieve N polys under each
valid strategy and keep the fastest) could follow; it would only pay off
for long jobs.

## 7. Build, test and compare

### 7.1 Build matrix script (`tools/build_matrix.sh`)

Builds into `bin/<tag>/yafu`. Each tag is one Makefile ISA level:
`USE_AVX512=1`, `USE_AVX2=1`, `USE_SSE41=1`, plain x86-64 (SSE2), and a
generic build with no intrinsics (a new `GENERIC=1` flag, also forceable
on x86). Strategies and kernels are run-time choices, so this is at most
five builds per machine rather than one per experiment.

### 7.2 Correctness gate (run before any timing)

For each binary and each valid strategy/kernel, on two inputs with
`-seed 777 -noopt -siqsNB 4`:

- **Same ranges** (same mid and large algorithm) must give identical
  "trial division touched", total-blocks, relation and polynomial counts,
  and the correct factors, **across all ISAs and kernels**.
- **Different ranges** (bucket vs sparse, bucket vs direct) only need the
  correct factors, and yield within a stated band; the bucket logp
  convention differs (direct_sieve.md 6.4).
- Keep the reference numbers in `tools/bench_ref.txt`, so a new machine
  is checked against them automatically.

### 7.3 Bench script (`tools/bench_matrix.sh`, a generalized bench.sh)

- **Inputs:** C70, C77 (256-bit), C79 (261-bit), C85, C95, kept in the
  repo with their factors.
- **Axes:** binaries × strategies × kernels × NB {2, 4, 8, default} ×
  reps (default 2).
- **Run order:** one job at a time, varying the innermost axis fastest so
  that drift is spread across variants. Delete `siqs.dat` before each
  run, use a fixed seed, and refuse to start if other yafu processes are
  running.
- **Output:** CSV (`machine, cpu, isa, mid, large, kernel, sparse_k, NB,
  input, rep, qs_sec, polys, rels`), plus a small summary script that
  prints a mean table with % against a chosen baseline, like the tables in
  this session.
- **Optional `PERF=1`:** wraps runs in `perf stat -e
  cycles,instructions,branch-misses` on Linux. On Windows-native builds,
  skip it.

### 7.4 Machines

- Rocket Lake (have): AVX-512, plus AVX2/SSE/generic builds for
  cross-ISA correctness and "lower ISA on the same core" timings.
- An AMD Zen 4 or 5: AVX-512 cmem vs creg, and AVX2.
- An AMD Zen 2/3 and an older Intel AVX2 core (Skylake client).
- ARM: Apple M-series or Graviton, generic build now, NEON later.

## 8. Phases and exit criteria

0. **Baseline and fences.**
   - Record reference outputs and timings of the current AVX-512 default.
   - Try a USE_AVX2 build of this branch and list what breaks.
   - Guard the AVX-512-only calls so lower-ISA builds fall back to stock
     buckets.
   - Fix L2 detection.
   - *Exit:* every ISA build compiles and factors with the stock strategy.
1. **Runtime strategy on AVX-512.**
   - Add the strategy struct, flags and validator.
   - Split process_poly into stream/direct.
   - Template the kernels.
   - *Exit:* each of the four (mid, large) combinations × kernels is
     bit-identical to the matching compile-time build from this branch, at
     equal speed (±1%).
2. **Scan hook and generic kernels.**
   - Add `siqs_report_tail` to all scans.
   - Write scalar versions of every new pass (direct pass 1, fused pass,
     sparse sieve, tdiv_large) with the sink and scalar-append kernels.
   - *Exit:* the generic build runs all four combinations, bit-identical
     to the AVX-512 build under the same strategy; ARM build verified.
3. **AVX2.**
   - Kernels: bitloop, sink, creg (LUT); 512-bit prefilter; scalar or
     gather resieve test.
   - *Exit:* bit-identical; timing matrix on Rocket Lake (AVX2 build) and
     at least one AMD Zen.
4. **SSE4.1 / SSE2.** Port only where measured faster than the scalar
   kernels. *Exit:* bit-identical; timing matrix.
5. **NEON.** *Exit:* bit-identical on ARM; timing against generic.
6. **Defaults and tuning.**
   - Fill the defaults table from the matrices.
   - Add strategy params to siqs_tune.py.
   - Write the results into direct_sieve.md.
   - *Exit:* on each tested machine, the default strategy is within ~2% of
     the best measured.

Phases 2 and 3 can run in parallel once phase 1 lands. Phase 2 is the
priority, because the generic build becomes the reference for all the
other ports.

## 9. Risks and open questions

- **Relying on the compiler to optimize.** Templating by include keeps the
  inner loops free of branches. If gcc fails to unroll or spills in the
  template versions, compare with the current hand-written AVX-512 code in
  phase 1. That comparison is the ±1% exit criterion.
- **Yield differences between strategies** (the logp per slice vs per
  16-group) mix calibration into speed comparisons. Consider a common logp
  convention (or `DIRECT_LOGP_BIAS` equivalent for sparse) so that
  strategies differ only in speed, or always report yield per poly
  alongside time.
- **Lower lane widths change the trade-offs.** With W = 8 or 4 the
  bitloop exits are cheaper and better predicted (fewer bits per mask),
  so the sink and compress gains may shrink. That is why every port gets
  its own timing matrix rather than inheriting AVX-512 conclusions.
- **Gather-free resieve on AVX2 and below** may make the fused pass
  costlier. If so, the direct sieve may lose to buckets on those ISAs,
  and the table should say so; the plan does not assume DIRECT wins
  everywhere.
- **Threads.** All timings so far are `-threads 1`. Direct and sparse
  have smaller per-thread memory than buckets, which may matter more with
  many threads sharing L3. Add a `-threads N` column for one input in each
  matrix.
- **Out of scope for now:** SS search, batch polys and TLP with DIRECT,
  and 64k blocks.
