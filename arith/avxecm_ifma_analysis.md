# AVX-ECM IFMA arithmetic: per-kernel analysis and the constant-S template

Branch `avxecm-ifma-opt` (off `wip_arith`), Nov 2026. Machine: i5-11600K
Rocket Lake (IFMA on port 0 only), WSL, gcc 11.4. Method and machine model
follow `arith/avx512_optimization_notes.md`.

## TL;DR

* The old fixed cores (`vecmulmod52_fixed416/624/832/1040_bfips` and
  their sqr versions) run at **54-75 % of the IFMA floor**. They are
  **latency-bound**, not port-bound: the in-block reduction chain
  serializes about 20-25 cycles per column. The general BFIPS
  `vecmulmod52`/`vecsqrmod52` reach 67-92 %.
* A new **constant-S template** (`factor/avx-ecm/vecarith52_tmpl.c`)
  uses column-wise Montgomery with lazy column accumulation. It is
  instantiated for S = 2..32 words plus a loop form for any larger S, and
  runs at **93-100 % of the IFMA floor** for S >= 8. Results are
  **bit-identical** to the old cores (the same AMM value), checked against
  GMP.
* powm per number, same word count, old -> new (-O2, yafu's flags):
  416 b **1.81x**, 624 b **1.49x**, 832 b **1.40x**, 1040 b **1.39x**,
  1248 b 1.33x, 1664 b 1.22x, 2080 b 1.17x, 3120 b 1.09x, 4160 b 1.09x.
* The template supports **any word count**, not only multiples of 4.
  Using it would remove the 4-word padding waste: a 512-bit N runs in
  10 words instead of 12, with powm 6.5 us vs 16.2 us today (**2.5x**).
  This needs driver changes (follow-up, see the last section).
* It is integrated into `assign_mul_ptrs` behind `AVXECM_TMPL_KERNELS`.
  That define is on by default for IFMA builds; build with
  `-DAVXECM_NO_TMPL_KERNELS` for the old cores.
* **Correctness finding (affects today's code):** when N is within 1-2
  bits of R = 2^(52*NWORDS), `vecaddmod52`/`vecsubmod52` on AMM outputs
  give wrong results for up to ~20 % of lanes. Today's sizing allows
  that, e.g. a 1038/1039-bit N in 1040 bits. With >= 3 spare bits there
  were no errors in 32k trials.
* Also found: `vecmulmod52_fixed1040_bfips` has a live
  `printf("mul > R, reducing")` in its final-subtraction path
  (vecarith52.c:1271). `vecsqrmod52_fixed1040_bfips` just calls the mul,
  and `vecsqrmod52_fixed416_bfips` has more IFMAs than the mul (272 vs 264).

## 1. Harness

`factor/avx-ecm/bench/` (`make GMP_PREFIX=...`; `./ae_bench -t` tests,
`-b -p` benchmarks with powm, `-a` add/sub range test, `-neg` negative
self-test). It compiles the AVX-ECM sources standalone, like
`arith/rns/rns_bench.c`.

* Correctness is **bit-exact** against a GMP model of the AMM product,
  `U = (ab + mN)/R`, minus N if `U >= R`. It checks all 8 lanes and:
  random inputs in [0, R); worst cases (R-1, N-1, 0, 1, sparse
  `rrandomb`); aliasing (c == a, c == b, sqr in place); 800-step
  dependent mul/sqr chains checked at every step; powm vs `mpz_powm` in
  every lane. Each size is run with three N: 1 spare bit, 2 spare bits,
  and 26 spare bits. The final subtraction fires on 0.5-4 % of random
  inputs with 1-2 spare bits, so that path is exercised.
* Negative self-test: a kernel wrapper that flips one bit every 50th call
  is caught at every size.
* Timing: min of 7 trials, pinned to one core with `taskset`. Cycles come
  from an IFMA throughput loop re-measured for every row, because the
  clock under WSL wanders between 4.2 and 4.6 GHz and single runs vary
  by about ±10 %. "dep" is a dependent chain (d = d*x); "ind" is
  independent consecutive calls.
* Floor = IFMA count at 1 IFMA/cycle: CIOS minimum mul 4S²+S, sqr 3S²+2S.
  The template with the carry trick needs S fewer (4S², 3S²+S); the
  tables keep the old floor.

## 2. Instruction mix (static, objdump) and bound

`max(IFMA, (IFMA + vector ALU)/2)` per §1 of the notes. "ALU" counts
vpadd/sub/and/or/xor/shift/ternlog/blend; "mov" counts loads, stores and
register moves.

| kernel | S | IFMA | ALU | mov | IFMA minimum | port bound | measured (cyc) | eff |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| vecmulmod52_fixed416 | 8 | 264 | 150 | 107 | 264 | 264 (IFMA) | 402-415 | 64-66 % |
| vecsqrmod52_fixed416 | 8 | 272 | 188 | 149 | 208 | 272 | 382-387 | 54 % |
| vecmulmod52_fixed624 | 12 | 588 | 226 | 192 | 588 | 588 | 831-871 | 67-71 % |
| vecsqrmod52_fixed624 | 12 | 456 | 268 | 142 | 456 | 456 | 718-719 | 63 % |
| vecmulmod52_fixed832 | 16 | 1040 | 302 | 289 | 1040 | 1040 | 1387-1659 | 63-75 % |
| vecsqrmod52_fixed832 | 16 | 864 | 376 | 401 | 800 | 864 | 1158-1203 | 67-69 % |
| vecmulmod52_fixed1040 | 20 | 1620 | 378 | 397 | 1620 | 1620 | 1732-1742 | 93-94 % |
| vecsqrmod52_fixed1040 | 20 | (calls mul) | | | 1240 | 1620 | 1789-2002 | 62-69 % |
| vecmulmod52 (BFIPS, loops) | 24-80 | dyn. 4S²+S | dyn. ~12/col | | 4S²+S | IFMA | | 74-92 % |
| vecsqrmod52 (BFIPS, loops) | 24-80 | dyn. 3S²+2S | | | 3S²+2S | IFMA | | 67-89 % |
| **vecmulmod52_s8 (tmpl)** | 8 | 264 | 132 | 359 | 264 | 264 | 265-283 | 93-100 % |
| **vecsqrmod52_s8** | 8 | 208 | 161 | 292 | 208 | 208 | 218-243 | 86-96 % |
| **vecmulmod52_s12** | 12 | 588 | 208 | 716 | 588 | 588 | 578-610 | 96-102 % |
| **vecsqrmod52_s12** | 12 | 456 | 261 | 567 | 456 | 456 | 469-492 | 93-97 % |
| **vecmulmod52_s16** | 16 | 1040 | 284 | 1169 | 1040 | 1040 | 1007-1064 | 98-103 % |
| **vecsqrmod52_s16** | 16 | 800 | 361 | 918 | 800 | 800 | 790-845 | 95-101 % |
| **vecmulmod52_s20** | 20 | 1620 | 360 | 1747 | 1620 | 1620 | 1585-1696 | 96-102 % |
| **vecsqrmod52_s20** | 20 | 1240 | 461 | 1356 | 1240 | 1240 | 1236-1297 | 96-100 % |
| **vecmulmod52_s24** | 24 | 2328 | 436 | 2480 | 2328 | 2328 | 2270-2416 | 96-103 % |
| **vecsqrmod52_s24** | 24 | 1776 | 561 | 1891 | 1776 | 1776 | 1821-1973 | 90-98 % |

(The template counts are without the carry trick. The trick removes S
IFMAs and adds 2S ALU.) The template's "mov" count is high mostly because
of `vmovdqa64 zero -> acc` copies that start each IFMA chain. They are
eliminated at rename and cost front-end slots only, at ~2.5 IPC.

Every kernel is **IFMA-bound by the port model**: ALU is always well
under IFMA, so `(IFMA+ALU)/2 < IFMA`. The fixed cores already have the
minimum IFMA count, except the sqr416/832/1040 cases above. Their 30-45 %
shortfall is therefore **latency**:

* dep ≈ ind for every old core, so the problem is not cross-call latency
  but a chain inside each call.
* In `CARRYPROP_ACCUM_n` (avx_ecm.h:~900-1000), each column of a 4-column
  block runs: add te_lo/te_hi into `acc_e0/acc_e1`; carry-propagate
  (srl, add, and); 0-3 IFMAs, *all into the same accumulator*; the m
  IFMA; the m*n0 IFMA; carry-propagate again; shift. That is a serial
  chain of ~20-28 cycles per column with ~2-3 IFMAs of independent work
  in it. S=12: 12 × ~24 ≈ 290 cycles exposed, which matches 871
  measured - 588 floor.
* fixed1040 mul (94 %) is the exception: at S=20 the out-of-order window
  finds enough block work (VEC_MUL4_ACCUM) to hide the chain.

## 3. The template (what worked)

`ae_montmul_t(C, A, B, N, rho, M, S, SQR, FULL, MREG)` is always-inline
with `const S`, instantiated by `AE_TMPL_INST(S)`.

Per column k = 0..2S-1:
```
lo = Σ lo(a_i b_{k-i}) + Σ_{i<k} lo(m_i n_{k-i})      (two accumulators each, ab and mn kept separate)
hi = Σ hi(a_i b_{k-i}) + Σ_{i<k} hi(m_i n_{k-i})
acc = lo + hi_{k-1} + carry
k <  S: m_k = lo52(acc * rho); hi += hi(m_k n0); carry = (acc>>52) + (acc mod 2^52 != 0)
k >= S: c_{k-S} = acc & mask;  carry = acc >> 52
final carry != 0 (rare, ~0.5-4 % of random inputs at 1-2 spare bits): masked c -= N, noinline slow path
```

What it does, against the lessons in the notes:

1. **Lazy column accumulation (§3.4).** One shift and 2-3 adds per column
   replace the old ~10 ALU ops per column. More importantly, the
   lo/hi/ab/mn **split accumulators** keep each column's dependency
   chains short and parallel. `m_{k-1}` (the newest value, on the
   critical path) is the last term of its chain and comes from a
   register, not a store-forward.
2. **Constant-S specialization (§3.3).** S = 2..32 instances; full unroll
   up to `AE_TMPL_FULL_MAX` = 24. Measured: the loop form is 5-10 %
   slower for sqr at S = 17-22 and equal at 24. Thanks to the unroll
   pragmas, **-O2 is as fast as -O3** (yafu builds at -O2).
3. **Squaring (§3.7):** off-diagonal products doubled with one add per
   lo/hi column, then the diagonal. 3S²+2S IFMAs, which is the minimum;
   sqr/mul ratio 0.73-0.80.
4. **Carry trick** (new; on by default, off with
   `AE_TMPL_NO_CARRY_TRICK`): acc + lo(m*n0) ≡ 0 mod 2^52, so its carry is
   `(acc>>52) + min(acc & mask, 1)`. This drops S IFMAs; measured 1-3 %.
5. **Rare final subtraction behind a mask test (§3.1).** The old cores
   always ran the masked subtraction loop: S × ~6 ops plus loads and
   stores.

## 4. Benchmarks (per 8-lane call; per-number ns = /8)

-O2 -march=icelake-client (yafu's flags), N with 1 spare bit:

| S | bits | old core | old sqr / mul ns | new sqr / mul ns | old powm us | new powm us | speedup |
|---:|---:|---|---:|---:|---:|---:|---:|
| 8 | 416 | fixed416 | 10.2 / 10.9 | 6.4 / 7.2 | 6.27 | 3.47 | 1.81x |
| 10 | 520 | (none: 624 used) | — | 9.4 / 11.1 | (16.15 at 624) | 6.50 | 2.48x vs today |
| 12 | 624 | fixed624 | 19.0 / 21.9 | 12.5 / 15.8 | 16.15 | 10.83 | 1.49x |
| 16 | 832 | fixed832 | 32.2 / 38.5 | 21.8 / 27.8 | 34.14 | 24.45 | 1.40x |
| 20 | 1040 | fixed1040 | 53.6 / 46.7 | 33.5 / 43.4 | 65.13 | 46.77 | 1.39x |
| 24 | 1248 | bfips | 66.8 / 78.7 | 49.5 / 63.8 | 107.65 | 81.23 | 1.33x |
| 32 | 1664 | bfips | 114.2 / 133.7 | 89.0 / 112.9 | 232.95 | 191.68 | 1.22x |
| 40 | 2080 | bfips | 165.8 / 204.1 | 141.0 / 176.9 | 431.84 | 368.36 | 1.17x |
| 60 | 3120 | bfips | 349.7 / 428.9 | 294.2 / 396.4 | 1320.7 | 1207.8 | 1.09x |
| 80 | 4160 | bfips | 593.8 / 761.4 | 520.9 / 721.7 | 3051.0 | 2807.7 | 1.09x |

Sizes that are not a multiple of 4 words (new only, -O3), sqr / mul ns per
number: S=5 2.9/3.1, 6 3.8/4.3, 7 6.0/6.2, 9 7.8/9.7, 11 12.6/14.8,
13 17.7/19.6, 14 16.8/21.4, 15 21.4/26.0, 17 24.6/31.7, 18 27.2/35.1,
22 41.4/52.6, 25 56.4/71.6, 27 63.3/83.7, 29 73.8/95.5, 30 84.7/103.1,
31 87.8/113.2, 33 98.7/125.3.

Comparison with the RNS Yao batch kernel (notes §6): at 2048 bits,
template mul is 177 ns vs the old bfips at 200, and at 4096 it is 722 vs
742-750. So the RNS advantage at 2-4k bits shrinks from 1.45x / 2.04x to
roughly 1.3x / 2.0x.

## 5. Candidates that did not help, or did not apply

* **Register-resident m (`AE_TMPL_MREG_MAX`)**: within noise at S <= 16
  (stores and reloads of m are off the critical path, since m_{k-1} comes
  from a register).
* **Full unroll beyond 24 words**: no gain at 28/32 and ~35 KB of code per
  instance. Left at 24.
* **Lookahead / software-pipelining of the ab sums**: not tried, because
  measured efficiency was already 93-100 %. The only kernels below 90 %
  are S <= 7, which are latency-bound (dep > ind by ~10 % at S=5). For
  those, the lever is interleaving two independent products (ECM always
  computes U and V in pairs), which would need a dual-product API.
* **4-word block padding**: the template has none of its own. The driver
  still rounds NWORDS to a multiple of 4 (follow-up).
* **Working set / L1**: at 80 words one operand is 5 KB, and A, B, N, M
  and C together are ~25 KB, inside the 48 KB L1. The loop form stays at
  94-98 % at S=80, so there was no L1 effect to fix. `vec_monty_alloc`
  allocates 1024 window tables (~10 MB at 80 words); that is unused
  memory, not a speed issue.
* **Karatsuba**: not attempted. At 2-4k bits the ab half is ~50 % of the
  IFMAs; one Karatsuba level saves ~12 % but needs 53-bit middle terms,
  which IFMA can't take. That means renormalizing to S/2+1 words. A
  possible follow-up for >= 2k bits.

## 6. Add/sub-mod range finding (`ae_bench -a`)

Values come from AMM mul/sqr chains started from residues < N, with sums
fed back in, as in curve arithmetic. Each run checks 32000 lanes:

| spare bits (52S - nbits) | addmod wrong | submod wrong |
|---:|---:|---:|
| 1 | 5120-7190 | 0-72 |
| 2 | 0-914 | 0-280 |
| >= 3 | 0 | 0 |

`vecaddmod52` subtracts N at most once and drops the carry if
a + b >= R + N. `vecsubmod52` adds N at most twice. Both depend on values
staying below ~2N, which AMM guarantees only with spare bits. The kernels
produce exactly the old values, so **current yafu is affected** whenever
nbits ≡ 206 or 207 (mod 208), e.g. 1038/1039-bit N in 1040 bits. The
likely effect is silently wasted curves rather than wrong factors (factors
are checked by gcd). Fix: size with at least 3 spare bits.

## 7. Integration status and follow-ups

* Done: `vecarith52_tmpl.c` in the Makefile and both VS projects (MSVC
  and clang-cl macros: `__forceinline`, `#pragma unroll`).
  `assign_mul_ptrs` overrides the mul/sqr pointers for non-Mersenne
  inputs with maxbits <= 16400, leaving add/sub selection unchanged.
  End-to-end yafu A/B: see §8.
* Follow-up 1 (correctness, small): size NWORDS with >= 3 spare bits in
  `avx_ecm_main.c` and the three sizing copies in `avxppm1.c`.
* Follow-up 2 (speed, larger): word-granular NWORDS
  (`ceil((nbits+3)/52)`) instead of multiples of 4. This is 1.2-2.5x on
  top of the template at most sizes. It needs an audit of the
  NBLOCKS/BLOCKWORDS assumptions in the stage-2, Mersenne and P-1 paths.
* Follow-up 3: the stray printf in fixed1040; the dead fixed cores could
  be removed once the template has been in use for a while.
* Not checked: the non-IFMA (FP-emulated) path, which is untouched; the
  template compiles to nothing without `IFMA`. MSVC build not compiled
  here.

## 8. End-to-end yafu (ECM, 1 thread, 8 curves)

Two WSL builds with yafu's own Makefile (-O2): default, and
`USER_CFLAGS=-DAVXECM_NO_TMPL_KERNELS`. `-v -v` confirms the switch
("Using constant-S template mul/sqr core (8 words)").

**Factor finding:** p15 · q at 400, 520, 600, 800, 1000, 1036, 1200, 1600
and 2048 bits, 64 curves at B1=3000. Both builds found the p15 at every
size.

**Timing:** B1=200000, B2=100·B1, one 8-curve batch, `taskset -c 3`,
mean of 2 runs, seconds:

| N bits | NWORDS | stage 1 old → new | speedup | stage 2 old → new | speedup |
|---:|---:|---:|---:|---:|---:|
| 400 | 8 | 0.302 → 0.219 | 1.38x | 0.254 → 0.202 | 1.26x |
| 520 | 12 | 0.569 → 0.430 | 1.32x | 0.374 → 0.321 | 1.16x |
| 600 | 12 | 0.553 → 0.424 | 1.31x | 0.384 → 0.328 | 1.17x |
| 800 | 16 | 0.914 → 0.730 | 1.25x | 0.576 → 0.502 | 1.15x |
| 1000 | 20 | 1.211 → 1.089 | 1.11x | 0.698 → 0.688* | ~1.0x |
| 1036 | 20 | 1.202 → 1.097 | 1.10x | 0.684 → 0.661 | 1.03x |
| 1200 | 24 | 1.799 → 1.521 | 1.18x | 0.968 → 0.866 | 1.12x |
| 1600 | 32 | 2.992 → 2.679 | 1.12x | 1.506 → 1.373 | 1.10x |
| 2048 | 40 | 4.465 → 4.042 | 1.10x | 2.141 → 1.962 | 1.09x |

(*) One new stage-2 run at 1000 bits was a 0.944 s outlier; the value
shown is the other run (0.688). End-to-end gains are below the kernel
gains, because add/sub, copies and ladder overhead are unchanged. At
1040 bits the old fixed1040 mul was already at 94 % (stage 2 is
mul-heavy), so there is little left there. A 520-bit N still runs in 12
words; word-granular sizing (follow-up 2) would give it 10.
