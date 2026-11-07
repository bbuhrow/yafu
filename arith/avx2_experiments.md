# AVX2 vector Montgomery: digit-format experiments (handoff 4.3)

Date: 2026-11-07. Code: `arith/avx2mont/` (standalone, not in yafu's build).

## Question

Can 4-lane AVX2 Montgomery multiplication beat `mpz_powm` per number, and
which digit format should an AVX2 port use? Candidates (user direction):
32-bit digits, 28-bit digits (lane accumulation headroom), and an FP FMA
approach.

## Variants

All three use the AVX-512 template's algorithm (column-wise product
scanning with lazy 64-bit column sums, fused Montgomery reduction, AMM
output with a rare final subtraction) on 4 lanes, with constant-S
specializations and a loop form. Each computes the same almost-Montgomery
product as AVX-ECM's cores, checked bit for bit against a GMP model.

| variant | digit | product | per product | notes |
|---|---:|---|---|---|
| FP52 | 52 bits | 2 x FP64 FMA: exact hi and lo | 2 FMA + 1 FP sub + 2 int adds | MXCSR set to round-toward-zero for the call; `x \| 0x433..` - 2^52 for uint64 -> double |
| D28 | 28 bits | `vpmuludq`, 56-bit product summed whole | 1 mul + 1 add | a*b and m*n sums in separate accumulators, so <= 255 terms each |
| D32 | 32 bits | `vpmuludq`, 64-bit product | 1 mul + 1 shift + 2 adds | full products summed mod 2^64 plus a sum of high halves; low sum recovered per column |

Uop estimate (relative cost ~ uops per product x number of products ~
1/digit^2): D28 2/28^2, D32 4/32^2 (1.5x D28), FP52 5/52^2 (0.72x D28).

FP52 details:
* AVX2 has no embedded rounding, so the kernel saves MXCSR, sets RZ, and
  restores it. Only the first FMA of each product depends on the mode;
  every other FP operation in the kernel is exact. The file is built with
  `-frounding-math` so gcc keeps the FP ops between the MXCSR writes.
* No `vcvtuqq2pd` and no `vpminuq` / `vpsraq` on AVX2: conversions use the
  2^52 bias trick, the low-column carry uses `vpcmpeqq`, the final
  subtraction (rare) is scalar.

## Correctness

`./a2m_bench -t`: every lane bit-exact against the GMP AMM model for all
three variants at 100 .. 5000 bits: random inputs in [0, R), R - 1, c == a,
500-step dependent mul/sqr chains, and sliding-window powm vs `mpz_powm`.
Passed with several seeds. The objects contain no zmm (AVX-512)
instructions.

## Speed

Rocket Lake i5-11600K running the AVX2 code (`-mavx2 -mfma`), `taskset`,
0.2 s warm-up, min of 3. Per number (4 numbers per call); powm with a
full-size exponent, sliding window.

| bits | FP52 S | FP52 powm us | D28 powm us | D32 powm us | mpz_powm us | FP52 vs GMP | D28 vs GMP |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 256 | 5 | 4.14 | 5.24 | 6.55 | 6.77 | 1.63x | 1.29x |
| 384 | 8 | 11.59 | 18.13 | 22.23 | 19.52 | 1.68x | 1.08x |
| 512 | 10 | 22.63 | 36.88 | 42.75 | 43.63 | 1.93x | 1.18x |
| 640 | 13 | 42.79 | 60.69 | 73.05 | 77.99 | 1.82x | 1.29x |
| 768 | 15 | 64.02 | 96.05 | 122.70 | 125.61 | 1.96x | 1.31x |
| 1024 | 20 | 137.31 | 267.78 | 301.29 | 257.95 | 1.88x | 0.96x |
| 1280 | 25 | 254.54 | 470.48 | 536.08 | 507.42 | 1.99x | 1.08x |
| 1536 | 30 | 422.63 | 772.35 | 867.85 | 809.68 | 1.92x | 1.05x |
| 1792 | 35 | 691.28 | 1199.73 | 1334.54 | 1303.63 | 1.89x | 1.09x |
| 2048 | 40 | 999.33 | 1738.20 | 1914.51 | 1892.63 | 1.89x | 1.09x |
| 2560 | 50 | 1847.54 | 3263.35 | 3580.10 | 3630.00 | 1.96x | 1.11x |
| 3072 | 60 | 3098.94 | 5421.24 | 6046.08 | 5519.29 | 1.78x | 1.02x |
| 3584 | 69 | 4663.70 | 9853.05 | 9578.00 | 9377.40 | 2.01x | 0.95x |
| 4096 | 79 | 6990.63 | 13320.73 | 14720.87 | 12749.04 | 1.82x | 0.96x |

For reference, same machine, powm per number (`rns_bench`): AVX-512 IFMA
template 43.8 us at 1024 bits, the AVX-512 FP-emulated template 84.1 us,
AVX2 FP52 137 us, mpz_powm 258 us.

Unrolling: with `A2M_FULL_MAX` 24 (the AVX-512 template's value), FP52 at
S = 20 was 25% slower than the loop form (16 ymm registers; a full unroll
at S = 20 is ~8000 instructions). 12 is best for D28 (S = 10), 8-12 for
FP52; the default is 12.

## Conclusions

1. **FP52 is the AVX2 format**: 1.6-2.0x faster than `mpz_powm` at every
   size from 256 to 4096 bits, and 1.6-1.9x faster than D28 at every size.
2. **D28** is about GMP parity (1.0-1.3x below 1 kbit, ~1.0-1.1x above):
   the 28-bit digits cost (52/28)^2 = 3.4x more products, which the cheaper
   product (2 uops instead of 5) doesn't make up.
3. **D32** is slower than GMP almost everywhere (0.83-1.03x): 4 uops per
   product and still 2.6x the products of FP52.
4. The 28-bit "lane headroom" argument holds (no carry passes needed up to
   S = 255, i.e. ~7100 bits), but headroom isn't the bottleneck; product
   count is.

## Caveats

* Measured on Rocket Lake executing AVX2 code, not on an AVX2-only CPU.
  On Zen 2/3 and Haswell..Skylake client, 256-bit FMA throughput is also 2
  per cycle, so the ranking should hold, but the ratios to GMP (whose
  `mulx`/`adx` code differs by core) need measuring there.
* MXCSR switching per call is cheap here (a few cycles against >= 4 us per
  powm per number) but it is a global-state dependency: callers must not
  run FP code that depends on the default rounding inside the kernel (none
  does), and MSVC needs `#pragma fenv_access(on)` in place of
  `-frounding-math`.

## Possible next steps (not done)

* An AVX2 AVX-ECM would need FP52 versions of the template and of every
  AVX-512-only routine AVX-ECM uses (add/sub, `vecmul52_1`, copies,
  Mersenne paths); AVX-ECM is AVX-512 only today (`USE_AVX512F` in
  `avx_ecm.h` and `factor_common.c`).
* RNS on AVX2: FP52 channels would reuse `rns_vec52.h`'s FP products at 4
  lanes; on AVX-512 the batch Yao kernel beats the template above ~1200
  bits in ECM, so a similar gain is plausible at large sizes.
