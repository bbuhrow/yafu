# AVX-512 IFMA arithmetic: optimization notes

Lessons from optimizing the RNS modmul kernels in `arith/rns/` (branch
`wip_arith`, Nov 2026), written as a starting point for the same kind of
work on AVX-ECM (`factor/avx-ecm/vecarith52*.c`). Design and full results
are in `arith/rns_design.md`. The numbers below are from an i5-11600K
(Rocket Lake), gcc 11.4 `-O3 -march=native` in WSL.

## 1. Machine model (client Ice/Rocket Lake)

Verify against uops.info before relying on it, especially on server parts
(ICX/SPR) and Zen 4, which differ.

| resource | what we relied on |
|---|---|
| `vpmadd52luq/huq` zmm | **port 0 only** on client cores: at most 1 IFMA per cycle. Latency ~4. |
| `vpaddq / vpandq / vpsrlq / vpsraq / vpsllvq` zmm | ports 0 and 5. When p0 is busy with IFMA, they all compete for p5. |
| shuffles (`vpermq`, `vpermi2q`, `vextracti64x*`, `vpsrldq`) | port 5 |
| loads | 2/cycle. A `{1to8}` embedded-broadcast memory operand (`_mm512_set1_epi64(mem)` folded into the IFMA) uses only the load port. |
| L1D | 48 KB. Working sets above that cost ~15% (see §4). |

Practical bound for an IFMA kernel:
`cycles >= max(IFMA, (IFMA + vector ALU + shuffles) / 2)`.
Once IFMA is roughly half of all vector µops, both ports are saturated
and *every* removed instruction helps. Below that, only removing IFMAs helps.

## 2. How we measured

* **Instruction mix:** dump the kernel and count mnemonics:
  ```
  objdump -d --no-show-raw-insn file.o | awk '/<sqr_5>:/,/^$/' | awk '{print $2}' | sort | uniq -c | sort -rn
  ```
  This is a *static* count. Rarely taken slow paths inflate it, so estimate
  the dynamic count by hand for the common path. Example, Yao s = 5 sqr:
  ~415 IFMA + ~300 ALU dynamic, ~560 cycles measured, against a 415-cycle
  IFMA floor.
* **Latency or throughput?** Run two independent dependency chains
  alternately (`f(a); f(b);` in a loop) and compare the per-op time with one
  chain. A big drop means latency-bound (interleaving or batching helps);
  no drop means throughput-bound (only fewer instructions help). Yao s = 2:
  -20% (latency-bound); s >= 3: ~0 (throughput-bound).
* **Rare paths:** count them with an `#ifdef` counter (we used
  `RNSY_COUNT_SLOW`). The fold slow path fired 0 times in 20k calls.
* **Timing hygiene:**
  - Warm up for ~0.2 s before each timed loop. The first loop on a new code
    path or size class ran 30-50% slow (cold code and tables, clock ramp),
    which produced a fake "squaring slower than multiplication" anomaly.
  - Expect +-5% run-to-run noise; repeat runs before acting on a 5% difference.
  - Use `clock_gettime(CLOCK_MONOTONIC)` and loops of >= 20k calls.
* **Cost models:** fit them to a measured sweep, don't derive them only from
  op counts. Our first model for picking the channel count missed by 13-18%
  at two sizes; a 6-parameter grid fit against one sweep got within 5% everywhere.

## 3. What worked (largest first)

1. **Fast path for rare carries (biggest single win).** A pseudo-Mersenne
   fold `X -> lo + hi*d` needs up to 3 rounds to guarantee `< 2^w`, each a
   full carry chain. After round 1 the remaining top is tiny, so rounds 2-3
   matter only if digit 1 overflows. One `_mm512_cmpge_epu64_mask` plus a
   `__builtin_expect` branch replaces two s-digit carry chains. Yao 2048-bit
   sqr went from 167 to 125 ns. General rule: when a branch-free worst-case
   sequence protects against an event with ~2^-40 probability, use a mask
   test and a predicted branch, and keep the exact slow path behind it.
2. **Batch layout: independent problems across lanes.** Putting 8 numbers in
   the 8 lanes (instead of 8 channels of one number) did three things:
   - It removed all horizontal reductions. Each `_mm512_reduce_add_epi64` is
     ~10 shuffle/extract µops on p5 plus latency.
   - Every constant became a broadcast operand.
   - The core got lots of independent work.
   Result: 2.5x at 512 bits and 1.8x at 1024; nothing at 2048, where the kernel
   was already throughput-bound. AVX-ECM already uses this layout (8 curves).
3. **Compile-time specialization plus full unrolling.** One `always_inline`
   template takes `const int S` (words per value). Per-S instances
   (`INSTANCES(1)..INSTANCES(10)`) go into a function-pointer table, with
   `_Pragma("GCC unroll 32")` on the digit loops. Accumulator arrays then
   live in registers. `-O3` vs `-O2` alone was **1.6x** on the
   small-regime kernel, because `-O2` didn't unroll the row-block loop. Check
   the flags first.
4. **Lazy column accumulation.** Accumulate the IFMA lo and hi halves of all
   products into 64-bit column accumulators. 52-bit values leave 12 spare
   bits, so up to 4096 terms fit. Normalize the carries once per output, not
   per product.
5. **Precompute per-modulus constant products in a better form.** For a
   product by a constant `f mod m`, store the rotations
   `F_i = f*2^(52i) mod m`. Then `a*f ≡ sum_i a_i*F_i` lands in s+1 columns
   with no high half to fold. That gained 5-9%. Precomputation that turns
   "multiply then reduce" into "multiply into already-reduced columns" is
   generally worth a look.
6. **Move work off the saturated port.** α' times a constant was replaced by
   a 16-entry table lookup (`_mm512_permutex2var_epi64` from two registers),
   moving 2s IFMA (p0) to s shuffles (p5).
7. **Squaring:** accumulate the off-diagonal products once, double the
   columns with one `vpaddq` each, then add the diagonal. In the RNS code
   this beat a separate off[] array (fewer live registers). AVX-ECM's
   sqr is already ~15-20% cheaper than its mul.
8. **Choose the decomposition per size.** When a size can be split several
   ways (here: channel count n vs words per channel s), sweep them.
   Capacity rounding (e.g. 3072 bits forced to s = 8 when s = 7.4 was
   needed) can cost more than the arithmetic differences.

## 4. What did not help (and why)

* **Reordering independent work in source** (issuing the C-base products
  alongside the B-base chain): no gain, since the out-of-order window
  already overlapped them. The extra live accumulators *cost* 4-10% at s = 10
  (register pressure / spills).
* **Register-block size** for matrix-vector loops (RB = 2, 4, 6, 8 rows):
  within noise once fully unrolled at -O3.
* **Trading ALU ops for IFMAs** (column-level fold: −4s ALU, +s IFMA per
  fold): only 1-4%, because the kernel was nearly IFMA-bound. Count per
  port before optimizing carry chains.
* **Interleaving two chains at s >= 3:** no gain (throughput-bound).
* **Large working sets:** the batch kernel at n = 8, s = 10 (~60 KB of
  operands, scratch and tables) was ~15% slower than a configuration that
  fit in the 48 KB L1.

## 5. Correctness methodology that made fast iteration safe

* **A scalar reference with the same reduction sequence** (GMP per value,
  the same fold rules), so the vector kernel can be checked **bit for bit**
  rather than only mod N. Most bugs show up immediately as a word mismatch
  at a specific step.
* **A `#define` switch** between the exact kernel and a faster inexact
  variant (`RNSY_COLUMN_FOLD`), plus a runtime query (`rnsy_vec_exact()`)
  so the tests pick the right checks automatically.
* **Inexact mode checks:** the same value mod N; every redundant component
  consistent with every other; digits normalized; value below its bound.
  A **negative self-test** (corrupt a word, check that the validator
  rejects it) proves the checks aren't vacuous.
* **Worst-case inputs** (max allowed value `h*k*N - 1`), in-place
  aliasing (`z == x`), long squaring chains, and powm against `mpz_powm`
  in every lane.
* Generalizing a parameter (n = 8 -> any n) exposed a hidden assumption: a
  sign `(-1)^(n-1)` that was constant at n = 8. Test the odd and edge cases
  of every parameter, not just the default.

## 6. Pointers for AVX-ECM

Baseline numbers measured this session (8 lanes, same N, per number):

| bits | AVX-ECM core | sqr ns | mul ns | powm (full exp) |
|---:|---|---:|---:|---:|
| 512 | fixed624 | 19-23 | 22-23 | 12.5 us |
| 1024 | fixed1040 | 46-49 | 46-47 | 57 us |
| 1536 | bfips (32 words) | 110 | 133 | 205 us |
| 2048 | bfips (40 words) | 163-167 | 200 | 404 us |
| 3072 | bfips (60 words) | 340 | 427 | 1250 us |
| 4096 | bfips (80 words) | 581-600 | 742-750 | 2840 us |

Observations and ideas to check:
* **Harness ready:** `arith/rns/rns_bench.c` already compiles
  `vec_common.c` and `vecarith52{,_common,_special}.c` standalone (flags
  `-DUSE_AVX512F -DUSE_AVX512BW -DSKYLAKEX -DUSE_AVX512IFMA -DIFMA`, plus an
  `xmalloc_align` shim). It mirrors `assign_mul_ptrs` and has a
  sliding-window powm that matches `vecslidingmodexp`. `vecmodexp` in
  avxppm1.c is unfinished (it frees `mdata->g`).
* **Fixed cores vs bfips:** the fixed416/624/832/1040 cores are much faster
  per word than the general BFIPS path. At 1536 bits (32 words) AVX-ECM is
  2.4x slower than at 1024 for 1.6x the words; quadratic scaling would
  predict 2.4x, so it is on trend, but more fixed-size instantiations (via
  a constant-S template as in §3.3) may help at 1248/1456/1664/2080 bits.
* **Block padding:** `vecmul52_n` / BFIPS work in 4-word blocks. For sizes
  that aren't a multiple of 4 words, padding wastes products. Check the
  padding and whether NWORDS is always rounded up to 4.
* **Count per port:** run the objdump count on `vecsqrmod52_fixed1040_bfips`
  and `vecmulmod52`. CIOS needs ~2n^2 lo/hi IFMA pairs per 8 numbers;
  compare the actual IFMA count and the add/and/shift overhead against
  that. If ALU >= IFMA, the carry handling is the lever (§3.1, §3.4); if
  IFMA dominates, look at Karatsuba (`veckmul_redc` exists; its crossover
  is quoted at ~16k bits on a Xeon 6254) or at reducing the reduction
  products.
* **Carry/normalization style:** look for per-word normalization inside the
  CIOS loop vs lazy column accumulation (§3.4), and for always-executed
  worst-case conditional subtractions that a mask test plus predicted
  branch could replace (§3.1).
* **Squaring:** AVX-ECM's sqr saves ~15-20% vs mul at 1.5k-4k bits.
  Check whether the theoretical ~40% for the multiply part is reached.
* **Working set:** at 80 words, one bignum is 80 x 64 B = 5 KB per value.
  Curve arithmetic keeps many temporaries plus `mdata->g` tables, so check
  L1/L2 behavior for ECM stage 1 and 2 at large sizes.
* **Higher level:** ECM op counts per bit (Edwards curves, alternate
  Montgomery parameterizations, stage-2 polynomial evaluation) are tracked
  in memory note `avx-ecm-edwards-idea`. Arithmetic speedups multiply with those.
* **Non-IFMA path:** the FP-emulated 52-bit path (`#ifndef IFMA`) uses FMA
  units. Port assignments differ (2 FMA ports on many server parts), so
  re-derive the bound in §1 for that path and target separately.
